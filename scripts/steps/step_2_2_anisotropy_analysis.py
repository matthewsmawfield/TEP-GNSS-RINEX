#!/usr/bin/env python3
"""
TEP-GNSS-RINEX Analysis - STEP 2.2: Anisotropy Analysis
========================================================

ALIGNED WITH CODE LONGSPAN METHODOLOGY (step_2_2_code_longspan.py):
- 8 directional sectors (N, NE, E, SE, S, SW, W, NW) at 45° each
- 40 log-spaced distance bins (50-13,000 km)
- 50 minimum pairs per bin for statistical reliability
- Distance distribution matching for sector fits
- E-W/N-S ratio analysis for Earth rotation coupling

RINEX-SPECIFIC EXTENSION - Multi-Metric Comparison:
Unlike CODE longspan which uses a single 'coherence' metric (cos(weighted_phase)),
this analysis compares BOTH:
1. MSC (Magnitude Squared Coherence) - normalized cross-spectral power
2. Phase Alignment Index - cos(weighted_phase), the CODE primary metric

Consistent anisotropy patterns across both metrics validates the signal.

Spectral Analysis (from step_2_0):
- Frequency Band: 10 µHz to 500 µHz (periods: 28 hours to 33 minutes)
- Detrend time series, compute CSD + auto-spectra via Welch's method
- Compute normalized MSC: |Pxy|² / (Pxx × Pyy)
- Magnitude-weighted circular phase averaging

Author: Matthew Lukin Smawfield
Date: 9 December 2025
License: CC-BY-4.0
"""

import sys
import json
import numpy as np
import urllib.request
import ssl
from pathlib import Path
from datetime import datetime, timedelta
from scipy.optimize import curve_fit
from scipy.signal import csd, welch
import matplotlib.pyplot as plt
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import pandas as pd
import gc

# Setup paths
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(ROOT))

from scripts.utils.logger import print_status, TEPLogger, set_step_logger
from scripts.utils.data_alignment import load_aligned_data, compute_aligned_coherence

# Directories
PROCESSED_DIR = ROOT / "data" / "processed"
OUTPUTS_DIR = ROOT / "results" / "outputs"
FIGURES_DIR = ROOT / "results" / "figures"

OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Initialize logger
logger = TEPLogger(
    name="step_2_2_anisotropy_analysis",
    level="DEBUG",
    log_file_path=ROOT / "logs" / "step_2_2_anisotropy_analysis.log"
)
set_step_logger(logger)

# ============================================================================
# ANALYSIS PARAMETERS - IDENTICAL TO step_2_0_raw_spp_analysis.py
# ============================================================================
# Frequency Band: 10 µHz to 500 µHz (periods: 28 hours to 33 minutes)
F1_HZ = 1e-5    # 10 µHz (28 hour period) - lower bound
F2_HZ = 5e-4    # 500 µHz (33 min period) - upper bound

# Sampling (5-minute intervals = 300 seconds)
SAMPLING_PERIOD_SEC = 300.0
FS_HZ = 1.0 / SAMPLING_PERIOD_SEC  # ~0.00333 Hz

# Minimum data requirements - MUST match step_2_0 for consistency
MIN_POINTS = 100

# Kp threshold for Quiet/Storm stratification
KP_QUIET_THRESHOLD = 3.0

# Distance Binning: Log-spaced, 50 km to 13,000 km - ALIGNED WITH CODE LONGSPAN
MIN_DISTANCE_KM = 50
MAX_DISTANCE_KM = 13000
N_BINS = 40

# Distance bands for additional stratified analysis
# (min_km, max_km, label)
DISTANCE_BANDS = [
    (50, 500, "50-500 km"),
    (500, 2000, "0.5-2 Mm"),
    (2000, 13000, ">2 Mm")
]
MIN_BIN_COUNT = 50  # Minimum pairs per bin for fitting (CODE uses 50)

# ============================================================================
# PROCESSING MODES - For comparative analysis across different data products
# ============================================================================
# Each mode has a data prefix (empty for baseline) and a description
PROCESSING_MODES = {
    'baseline': {
        'prefix': '',           # No prefix - uses 'clock_bias_ns' directly
        'key': 'clock_bias_ns',
        'description': 'SPP with Broadcast Ephemeris (Standard)',
        'removes_iono': False,
        'removes_geometry': False,
    },
    'precise': {
        'prefix': 'precise_',    # Uses 'precise_clock_bias_ns'
        'key': 'precise_clock_bias_ns',
        'description': 'SPP with Precise Orbits (Option E)',
        'removes_iono': False,
        'removes_geometry': True,  # Precise orbits remove geometric bias
    },
    'ionofree': {
        'prefix': 'ionofree_',   # Uses 'ionofree_clock_bias_ns'
        'key': 'ionofree_clock_bias_ns',
        'description': 'Dual-Freq Iono-Free with Precise Orbits (Option D)',
        'removes_iono': True,
        'removes_geometry': True,
    },
    'multi_gnss': {
        'prefix': 'multi_',      # Uses 'multi_clock_bias_ns'
        'key': 'multi_clock_bias_ns',
        'description': 'Multi-GNSS (GPS+GLO+GAL+BDS)',
        'removes_iono': False,
        'removes_geometry': False,
    },
}

# CODE reference E-W/N-S ratio (from PPP analysis)
CODE_EW_NS_RATIO = 2.16

# =============================================================================
# CODE LONGSPAN REFERENCE VALUES - For geometric correction analysis
# =============================================================================
# These are the λ values from CODE's 25-year PPP analysis (step_2_2_geospatial_temporal_analysis_code.json)
# Used to compute geometric suppression factors in SPP data
CODE_SECTOR_LAMBDAS = {
    'N': 2313,   # km
    'NE': 2540,
    'E': 3206,
    'SE': 6808,
    'S': 2719,
    'SW': 5332,
    'W': 7664,
    'NW': 3028
}
CODE_EW_MEAN = (CODE_SECTOR_LAMBDAS['E'] + CODE_SECTOR_LAMBDAS['W']) / 2  # 5435 km
CODE_NS_MEAN = (CODE_SECTOR_LAMBDAS['N'] + CODE_SECTOR_LAMBDAS['S']) / 2  # 2516 km

# =============================================================================
# PHASE METRICS TO COMPARE - RINEX-specific multi-metric analysis
# =============================================================================
# Unlike CODE longspan which uses a single 'coherence' metric (cos(weighted_phase)),
# RINEX analysis compares multiple phase metrics to validate consistency:
PHASE_METRICS = {
    'coherence': {
        'column': 'coherence',
        'description': 'Magnitude Squared Coherence (MSC)',
        'range': (0, 1),
        'interpretation': 'Normalized cross-spectral power'
    },
    'phase_alignment': {
        'column': 'phase_alignment',
        'description': 'Phase Alignment Index (cos(weighted_phase))',
        'range': (-1, 1),
        'interpretation': 'Directional phase coupling (CODE methodology)'
    }
}

# =============================================================================
# SECTOR CLASSIFICATION - ALIGNED WITH CODE LONGSPAN METHODOLOGY
# =============================================================================
# Use 45° sectors (±22.5° from cardinal directions) to match CODE longspan.
# E: [67.5°, 112.5°), W: [247.5°, 292.5°)
# N: [337.5°, 360°) ∪ [0°, 22.5°), S: [157.5°, 202.5°)
#
# Previous RINEX implementation used 90° sectors (±45°) which diluted
# the anisotropy signal by including diagonal pairs.
# =============================================================================

def is_ew_azimuth(az):
    """Check if azimuth is in E-W sector (45° width, ±22.5° from cardinal)."""
    return (67.5 <= az < 112.5) or (247.5 <= az < 292.5)

def is_ns_azimuth(az):
    """Check if azimuth is in N-S sector (45° width, ±22.5° from cardinal)."""
    return (az < 22.5) or (157.5 <= az < 202.5) or (az >= 337.5)


def compute_cross_power_plateau(series1, series2, fs, f1=F1_HZ, f2=F2_HZ):
    """
    Compute normalized coherence using CSD and auto-spectra (MSC methodology).
    
    IDENTICAL TO step_2_0_raw_spp_analysis.py for consistency.
    Uses Magnitude Squared Coherence: |Pxy|² / (Pxx × Pyy)
    
    Parameters:
    -----------
    series1, series2 : array-like
        Time series to correlate (same length, synchronized)
    fs : float
        Sampling frequency in Hz
    f1, f2 : float
        TEP frequency band limits in Hz
        
    Returns:
    --------
    mean_coherence : float
        Mean normalized coherence in band (0-1)
    weighted_phase : float
        Magnitude-weighted average phase (radians)
    """
    n_points = len(series1)
    if n_points < MIN_POINTS:
        return np.nan, np.nan
    
    # STEP 1: Detrend time series (remove linear trend)
    time_indices = np.arange(n_points)
    series1_detrended = series1 - np.polyval(np.polyfit(time_indices, series1, 1), time_indices)
    series2_detrended = series2 - np.polyval(np.polyfit(time_indices, series2, 1), time_indices)
    
    # STEP 2: Compute CSD AND auto-spectra via Welch's method
    nperseg = min(1024, n_points)
    try:
        frequencies, Pxy = csd(series1_detrended, series2_detrended,
                               fs=fs, nperseg=nperseg, detrend='constant')
        _, Pxx = welch(series1_detrended, fs=fs, nperseg=nperseg, detrend='constant')
        _, Pyy = welch(series2_detrended, fs=fs, nperseg=nperseg, detrend='constant')
    except Exception:
        return np.nan, np.nan
    
    if len(frequencies) < 2:
        return np.nan, np.nan
    
    # STEP 3: Select TEP frequency band
    band_mask = (frequencies > 0) & (frequencies >= f1) & (frequencies <= f2)
    if not np.any(band_mask):
        return np.nan, np.nan
    
    # STEP 4: Compute NORMALIZED COHERENCE = |Pxy|² / (Pxx × Pyy)
    Pxy_band = Pxy[band_mask]
    Pxx_band = Pxx[band_mask]
    Pyy_band = Pyy[band_mask]
    
    # Avoid division by zero
    denom = Pxx_band * Pyy_band
    valid_mask = denom > 0
    if not np.any(valid_mask):
        return np.nan, np.nan
    
    # Coherence squared (MSC - Magnitude Squared Coherence)
    coh_squared = (np.abs(Pxy_band[valid_mask])**2) / denom[valid_mask]
    
    # Phase from cross-spectrum
    phases = np.angle(Pxy_band[valid_mask])
    magnitudes = np.sqrt(coh_squared)  # Coherence (0-1)
    
    if len(magnitudes) == 0 or np.sum(magnitudes) == 0:
        return np.nan, np.nan
    
    # STEP 5: Magnitude-weighted circular phase averaging
    complex_phases = np.exp(1j * phases)
    weighted_complex = np.average(complex_phases, weights=magnitudes)
    weighted_phase = np.angle(weighted_complex)
    
    # Mean coherence in band (normalized, 0-1)
    mean_coherence = float(np.mean(magnitudes))
    
    return mean_coherence, float(weighted_phase)


def ecef_to_lla(x, y, z):
    """Convert ECEF to lat/lon."""
    import math
    lon = math.atan2(y, x)
    p = math.sqrt(x**2 + y**2)
    e2 = 0.00669437999014
    lat = math.atan2(z, p * (1 - e2))
    for _ in range(5):
        N = 6378137.0 / math.sqrt(1 - e2 * math.sin(lat)**2)
        lat = math.atan2(z + e2 * N * math.sin(lat), p)
    return math.degrees(lat), math.degrees(lon)


def haversine(lat1, lon1, lat2, lon2):
    """Calculate distance in km."""
    import math
    R = 6371
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat/2)**2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon/2)**2
    return R * 2 * math.asin(math.sqrt(a))


def compute_azimuth(lat1, lon1, lat2, lon2):
    """Compute azimuth from point 1 to point 2 in degrees (0-360)."""
    import math
    dlon = math.radians(lon2 - lon1)
    lat1_r = math.radians(lat1)
    lat2_r = math.radians(lat2)
    
    x = math.sin(dlon) * math.cos(lat2_r)
    y = math.cos(lat1_r) * math.sin(lat2_r) - math.sin(lat1_r) * math.cos(lat2_r) * math.cos(dlon)
    
    azimuth = math.degrees(math.atan2(x, y))
    return (azimuth + 360) % 360


def exp_decay(r, A, lam, C0):
    """Exponential decay model."""
    return A * np.exp(-r / lam) + C0


def fetch_kp_data_gfz(years) -> dict:
    """Fetch Kp data from GFZ Potsdam for given year(s)."""
    url = "https://www-app3.gfz-potsdam.de/kp_index/Kp_ap_since_1932.txt"
    
    if isinstance(years, int):
        years = [years]
        
    print_status(f"Fetching Kp data from GFZ for years: {years}...", "INFO")
    
    kp_by_date = {}
    
    try:
        ssl_context = ssl.create_default_context()
        with urllib.request.urlopen(url, context=ssl_context, timeout=60) as response:
            content = response.read().decode('utf-8').splitlines()
        
        for line in content:
            if line.startswith('#') or not line.strip():
                continue
            
            try:
                parts = line.split()
                if len(parts) < 8:
                    continue
                
                y, m, d = int(parts[0]), int(parts[1]), int(parts[2])
                if y not in years:
                    continue
                
                kp_str = parts[7]
                if kp_str == '-1' or kp_str == '':
                    continue
                
                if kp_str.endswith('+'):
                    kp = float(kp_str[:-1]) + 0.33
                elif kp_str.endswith('-'):
                    kp = float(kp_str[:-1]) - 0.33
                elif kp_str.endswith('o'):
                    kp = float(kp_str[:-1])
                else:
                    kp = float(kp_str)
                
                date_str = f"{y:04d}{m:02d}{d:02d}"
                
                if date_str not in kp_by_date:
                    kp_by_date[date_str] = []
                kp_by_date[date_str].append(kp)
                
            except (ValueError, IndexError):
                continue
        
        # Use DAILY MAX Kp
        kp_daily = {d: np.max(vals) for d, vals in kp_by_date.items()}
        print_status(f"  Loaded Kp for {len(kp_daily)} days (using daily MAX)", "INFO")
        return kp_daily
        
    except Exception as e:
        raise RuntimeError(f"CRITICAL: Failed to fetch real Kp data from GFZ: {e}. "
                          f"Cannot proceed without geomagnetic data.")


def filter_series_by_kp(series, kp_data, max_kp=None, min_kp=None):
    """
    Filter a Pandas Series based on Kp values for each day.
    Uses the Series DatetimeIndex.
    """
    if series.empty:
        return series
        
    # Create a map of valid dates
    valid_dates = set()
    for date_str, kp in kp_data.items():
        if max_kp is not None and kp >= max_kp:
            continue
        if min_kp is not None and kp < min_kp:
            continue
        valid_dates.add(date_str)
        
    # Filter
    # Optimized: convert index to YYYYMMDD strings
    # (This might be slow for very long series, but robust)
    date_strings = series.index.strftime('%Y%m%d')
    mask = [d in valid_dates for d in date_strings]
    
    return series[mask]


def _process_day_batch_worker(args):
    """
    Worker function to process a batch of days day-by-day.
    Loads daily files, computes coherence for all pairs, and returns results.
    Matches CODE Longspan methodology.
    """
    day_keys, day_files_map, mode_key, coords = args
    results = []
    
    import pandas as pd
    from scipy.signal import csd, welch
    import math
    
    # Local helpers to ensure availability
    def _haversine(lat1, lon1, lat2, lon2):
        R = 6371
        dlat = math.radians(lat2 - lat1)
        dlon = math.radians(lon2 - lon1)
        a = math.sin(dlat/2)**2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon/2)**2
        return R * 2 * math.asin(math.sqrt(a))

    def _azimuth(lat1, lon1, lat2, lon2):
        dlon = math.radians(lon2 - lon1)
        y = math.sin(dlon) * math.cos(math.radians(lat2))
        x = math.cos(math.radians(lat1)) * math.sin(math.radians(lat2)) - \
            math.sin(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.cos(dlon)
        az = math.degrees(math.atan2(y, x))
        return (az + 360) % 360
    
    # DOY to month lookup (approximate)
    doy_to_month = lambda doy: (
        1 if doy <= 31 else
        2 if doy <= 59 else
        3 if doy <= 90 else
        4 if doy <= 120 else
        5 if doy <= 151 else
        6 if doy <= 181 else
        7 if doy <= 212 else
        8 if doy <= 243 else
        9 if doy <= 273 else
        10 if doy <= 304 else
        11 if doy <= 334 else 12
    )
    
    for day_key in day_keys:
        year, doy = day_key
        month = doy_to_month(doy)
        files_for_day = day_files_map.get(day_key, {})
        station_data = {}
        
        # 1. Load data for this day (fast)
        for sta, fpath in files_for_day.items():
            try:
                with np.load(fpath) as d:
                    if mode_key in d:
                        val = d[mode_key]
                        if len(val) > 200 and np.sum(np.isfinite(val)) > 200:
                            ts = d.get('timestamps')
                            if ts is not None:
                                valid = (ts >= 0) & (ts < 86400)
                                if np.any(valid):
                                    s = pd.Series(val[valid], index=ts[valid])
                                    station_data[sta] = s
            except:
                continue
        
        stations = list(station_data.keys())
        if len(stations) < 2: continue
        
        # 2. Compute all pairs for this day
        for i, sta1 in enumerate(stations):
            s1 = station_data[sta1]
            lat1, lon1 = coords[sta1]['lat'], coords[sta1]['lon']
            
            for sta2 in stations[i+1:]:
                s2 = station_data[sta2]
                lat2, lon2 = coords[sta2]['lat'], coords[sta2]['lon']
                
                dist = _haversine(lat1, lon1, lat2, lon2)
                if not (MIN_DISTANCE_KM <= dist <= MAX_DISTANCE_KM):
                    continue
                
                # Fast intersection on index (timestamps)
                common_idx = s1.index.intersection(s2.index)
                if len(common_idx) < 200: continue
                
                # Safe extraction ensuring equal length
                d1_series = s1[common_idx]
                d2_series = s2[common_idx]
                
                if len(d1_series) != len(d2_series):
                    # Fallback: align by dataframe to handle duplicates
                    try:
                        df = pd.DataFrame({'d1': s1, 'd2': s2}).dropna()
                        if len(df) < 200: continue
                        d1 = df['d1'].values
                        d2 = df['d2'].values
                    except:
                        continue
                else:
                    d1 = d1_series.values
                    d2 = d2_series.values
                
                # Detrend (Linear) - CODE Step 1
                x = np.arange(len(d1))
                d1 = d1 - np.polyval(np.polyfit(x, d1, 1), x)
                d2 = d2 - np.polyval(np.polyfit(x, d2, 1), x)
                
                # CSD (Constant Detrend) - CODE Step 2
                try:
                    # Ensure sufficient segments for Welch's method (variance reduction)
                    # For daily files (288 samples), large nperseg yields 1 segment -> Coherence=1.0 (useless)
                    n_points = len(d1)
                    if n_points < 64: continue 
                    
                    # Target 4 segments minimum, cap at 256 (to handle daily files)
                    target_nperseg = n_points // 4
                    nperseg = min(target_nperseg, 256)
                    if nperseg < 32: nperseg = 32
                    
                    f, Pxy = csd(d1, d2, fs=FS_HZ, nperseg=nperseg, detrend='constant')
                    _, Pxx = welch(d1, fs=FS_HZ, nperseg=nperseg, detrend='constant')
                    _, Pyy = welch(d2, fs=FS_HZ, nperseg=nperseg, detrend='constant')
                    
                    # TEP Band: 10-500 uHz
                    mask = (f >= F1_HZ) & (f <= F2_HZ)
                    if not np.any(mask): continue
                    
                    # Coherence
                    coh = np.abs(Pxy[mask])**2 / (Pxx[mask] * Pyy[mask])
                    mean_coh = np.mean(np.sqrt(coh))
                    
                    if np.isfinite(mean_coh):
                        az = _azimuth(lat1, lon1, lat2, lon2)
                        mid_lat = (lat1 + lat2) / 2
                        results.append({
                            'dist': dist,
                            'coherence': float(mean_coh),
                            'azimuth': az,
                            'mid_lat': mid_lat,
                            'year': year,
                            'month': month
                        })
                except:
                    continue
                    
    return results


def _subsample_to_match_distribution(sector_distances, reference_distances, max_samples=50000):
    """
    Subsample sector distances to match the reference (global) distance distribution.
    This prevents bias in λ estimates from differing distance distributions across sectors.
    
    IDENTICAL TO CODE LONGSPAN METHODOLOGY.
    """
    if len(sector_distances) == 0:
        return np.array([], dtype=int)
    
    # Compute reference (global) histogram using LOG-SPACED bins to match fitting
    # This is critical: Linear bins (50-13000km) are too coarse at short distances
    ref_bins = np.logspace(np.log10(MIN_DISTANCE_KM), np.log10(MAX_DISTANCE_KM), N_BINS + 1)
    ref_hist, _ = np.histogram(reference_distances, bins=ref_bins, density=False)
    
    # Compute target counts per bin
    total_sector_pairs = len(sector_distances)
    target_samples = min(max_samples, total_sector_pairs)
    
    # Calculate target counts per bin based on reference distribution
    # ref_hist is count (density=False), so dividing by sum gives probability
    ref_prob = ref_hist / ref_hist.sum() if ref_hist.sum() > 0 else np.ones(N_BINS) / N_BINS
    target_per_bin = ref_prob * target_samples
    
    # Sample from each bin
    selected_indices = []
    rng = np.random.default_rng(42)
    for i in range(len(ref_bins) - 1):
        bin_mask = (sector_distances >= ref_bins[i]) & (sector_distances < ref_bins[i+1])
        bin_indices = np.where(bin_mask)[0]
        
        if len(bin_indices) > 0:
            n_to_sample = int(min(len(bin_indices), target_per_bin[i]))
            if n_to_sample > 0:
                sampled = rng.choice(bin_indices, size=n_to_sample, replace=False)
                selected_indices.extend(sampled)
    
    return np.array(selected_indices, dtype=int)


def _subsample_to_match_distribution_2d(distances, lats, reference_distances, reference_lats,
                                       max_samples=50000, dist_bins=None, lat_bins=None):
    if len(distances) == 0:
        return np.array([], dtype=int)
    if len(reference_distances) == 0:
        return np.array([], dtype=int)

    if dist_bins is None:
        dist_bins = np.logspace(np.log10(MIN_DISTANCE_KM), np.log10(MAX_DISTANCE_KM), N_BINS + 1)
    if lat_bins is None:
        lat_bins = np.array([-90.0, -60.0, -30.0, 0.0, 30.0, 60.0, 90.0])

    ref_hist, _, _ = np.histogram2d(reference_distances, reference_lats, bins=[dist_bins, lat_bins], density=False)
    ref_sum = float(np.sum(ref_hist))
    if ref_sum > 0:
        ref_prob = ref_hist / ref_sum
    else:
        ref_prob = np.ones_like(ref_hist, dtype=float)
        ref_prob /= float(ref_prob.size)

    target_samples = int(min(max_samples, len(distances)))
    target_per_cell = ref_prob * float(target_samples)

    selected_indices = []
    rng = np.random.default_rng(42)
    for i in range(len(dist_bins) - 1):
        d_lo = dist_bins[i]
        d_hi = dist_bins[i + 1]
        d_mask = (distances >= d_lo) & (distances < d_hi)
        if not np.any(d_mask):
            continue
        for j in range(len(lat_bins) - 1):
            lat_lo = lat_bins[j]
            lat_hi = lat_bins[j + 1]
            cell_idx = np.where(d_mask & (lats >= lat_lo) & (lats < lat_hi))[0]
            if len(cell_idx) == 0:
                continue
            n_to_sample = int(min(len(cell_idx), np.floor(target_per_cell[i, j] + 1e-9)))
            if n_to_sample <= 0:
                continue
            sampled = rng.choice(cell_idx, size=n_to_sample, replace=False)
            selected_indices.extend(sampled)

    return np.array(selected_indices, dtype=int)


def _fit_matched_directional_lambda(distances, coherences, ew_mask, ns_mask, reference_distances=None,
                                   bins=None, min_pairs=MIN_BIN_COUNT, use_sigma=True, min_bins=5,
                                   max_samples=50000, lats=None, reference_lats=None, match_lat=False, lat_bins=None):
    ew_d = distances[ew_mask]
    ew_c = coherences[ew_mask]
    ns_d = distances[ns_mask]
    ns_c = coherences[ns_mask]

    ew_lat = lats[ew_mask] if lats is not None else None
    ns_lat = lats[ns_mask] if lats is not None else None

    if len(ew_d) < 50 or len(ns_d) < 50:
        return None, None, {
            'n_ew_total': int(len(ew_d)),
            'n_ns_total': int(len(ns_d)),
            'n_ew_matched': 0,
            'n_ns_matched': 0,
            'ew_dist_median_km_raw': float(np.median(ew_d)) if len(ew_d) else None,
            'ns_dist_median_km_raw': float(np.median(ns_d)) if len(ns_d) else None,
            'ew_dist_median_km_matched': None,
            'ns_dist_median_km_matched': None,
            'ew_lat_median_raw': float(np.median(ew_lat)) if ew_lat is not None and len(ew_lat) else None,
            'ns_lat_median_raw': float(np.median(ns_lat)) if ns_lat is not None and len(ns_lat) else None,
            'ew_abs_lat_median_raw': float(np.median(np.abs(ew_lat))) if ew_lat is not None and len(ew_lat) else None,
            'ns_abs_lat_median_raw': float(np.median(np.abs(ns_lat))) if ns_lat is not None and len(ns_lat) else None,
            'ew_north_frac_raw': float(np.mean(ew_lat > 0.0)) if ew_lat is not None and len(ew_lat) else None,
            'ns_north_frac_raw': float(np.mean(ns_lat > 0.0)) if ns_lat is not None and len(ns_lat) else None,
            'ew_lat_median_matched': None,
            'ns_lat_median_matched': None,
            'ew_abs_lat_median_matched': None,
            'ns_abs_lat_median_matched': None,
            'ew_north_frac_matched': None,
            'ns_north_frac_matched': None,
        }

    if reference_distances is None:
        reference_distances = np.concatenate([ew_d, ns_d])
    if reference_lats is None and lats is not None:
        reference_lats = np.concatenate([ew_lat, ns_lat]) if ew_lat is not None and ns_lat is not None else None
    if reference_distances is None or len(reference_distances) == 0:
        return None, None, {
            'n_ew_total': int(len(ew_d)),
            'n_ns_total': int(len(ns_d)),
            'n_ew_matched': 0,
            'n_ns_matched': 0,
            'ew_dist_median_km_raw': float(np.median(ew_d)) if len(ew_d) else None,
            'ns_dist_median_km_raw': float(np.median(ns_d)) if len(ns_d) else None,
            'ew_dist_median_km_matched': None,
            'ns_dist_median_km_matched': None,
            'ew_lat_median_raw': float(np.median(ew_lat)) if ew_lat is not None and len(ew_lat) else None,
            'ns_lat_median_raw': float(np.median(ns_lat)) if ns_lat is not None and len(ns_lat) else None,
            'ew_abs_lat_median_raw': float(np.median(np.abs(ew_lat))) if ew_lat is not None and len(ew_lat) else None,
            'ns_abs_lat_median_raw': float(np.median(np.abs(ns_lat))) if ns_lat is not None and len(ns_lat) else None,
            'ew_north_frac_raw': float(np.mean(ew_lat > 0.0)) if ew_lat is not None and len(ew_lat) else None,
            'ns_north_frac_raw': float(np.mean(ns_lat > 0.0)) if ns_lat is not None and len(ns_lat) else None,
            'ew_lat_median_matched': None,
            'ns_lat_median_matched': None,
            'ew_abs_lat_median_matched': None,
            'ns_abs_lat_median_matched': None,
            'ew_north_frac_matched': None,
            'ns_north_frac_matched': None,
        }

    if match_lat and ew_lat is not None and ns_lat is not None and reference_lats is not None:
        ew_idx = _subsample_to_match_distribution_2d(
            ew_d,
            ew_lat,
            reference_distances,
            reference_lats,
            max_samples=min(max_samples, len(ew_d)),
            dist_bins=bins,
            lat_bins=lat_bins,
        )
        ns_idx = _subsample_to_match_distribution_2d(
            ns_d,
            ns_lat,
            reference_distances,
            reference_lats,
            max_samples=min(max_samples, len(ns_d)),
            dist_bins=bins,
            lat_bins=lat_bins,
        )
    else:
        ew_idx = _subsample_to_match_distribution(ew_d, reference_distances, max_samples=min(max_samples, len(ew_d)))
        ns_idx = _subsample_to_match_distribution(ns_d, reference_distances, max_samples=min(max_samples, len(ns_d)))

    if len(ew_idx) < 50 or len(ns_idx) < 50:
        return None, None, {
            'n_ew_total': int(len(ew_d)),
            'n_ns_total': int(len(ns_d)),
            'n_ew_matched': int(len(ew_idx)),
            'n_ns_matched': int(len(ns_idx)),
            'ew_dist_median_km_raw': float(np.median(ew_d)) if len(ew_d) else None,
            'ns_dist_median_km_raw': float(np.median(ns_d)) if len(ns_d) else None,
            'ew_dist_median_km_matched': float(np.median(ew_d[ew_idx])) if len(ew_idx) else None,
            'ns_dist_median_km_matched': float(np.median(ns_d[ns_idx])) if len(ns_idx) else None,
            'ew_lat_median_raw': float(np.median(ew_lat)) if ew_lat is not None and len(ew_lat) else None,
            'ns_lat_median_raw': float(np.median(ns_lat)) if ns_lat is not None and len(ns_lat) else None,
            'ew_abs_lat_median_raw': float(np.median(np.abs(ew_lat))) if ew_lat is not None and len(ew_lat) else None,
            'ns_abs_lat_median_raw': float(np.median(np.abs(ns_lat))) if ns_lat is not None and len(ns_lat) else None,
            'ew_north_frac_raw': float(np.mean(ew_lat > 0.0)) if ew_lat is not None and len(ew_lat) else None,
            'ns_north_frac_raw': float(np.mean(ns_lat > 0.0)) if ns_lat is not None and len(ns_lat) else None,
            'ew_lat_median_matched': float(np.median(ew_lat[ew_idx])) if ew_lat is not None and len(ew_idx) else None,
            'ns_lat_median_matched': float(np.median(ns_lat[ns_idx])) if ns_lat is not None and len(ns_idx) else None,
            'ew_abs_lat_median_matched': float(np.median(np.abs(ew_lat[ew_idx]))) if ew_lat is not None and len(ew_idx) else None,
            'ns_abs_lat_median_matched': float(np.median(np.abs(ns_lat[ns_idx]))) if ns_lat is not None and len(ns_idx) else None,
            'ew_north_frac_matched': float(np.mean(ew_lat[ew_idx] > 0.0)) if ew_lat is not None and len(ew_idx) else None,
            'ns_north_frac_matched': float(np.mean(ns_lat[ns_idx] > 0.0)) if ns_lat is not None and len(ns_idx) else None,
        }

    ew_res = fit_exponential(ew_d[ew_idx], ew_c[ew_idx], min_pairs=min_pairs, bins=bins, use_sigma=use_sigma, min_bins=min_bins)
    ns_res = fit_exponential(ns_d[ns_idx], ns_c[ns_idx], min_pairs=min_pairs, bins=bins, use_sigma=use_sigma, min_bins=min_bins)

    return ew_res, ns_res, {
        'n_ew_total': int(len(ew_d)),
        'n_ns_total': int(len(ns_d)),
        'n_ew_matched': int(len(ew_idx)),
        'n_ns_matched': int(len(ns_idx)),
        'ew_dist_median_km_raw': float(np.median(ew_d)) if len(ew_d) else None,
        'ns_dist_median_km_raw': float(np.median(ns_d)) if len(ns_d) else None,
        'ew_dist_median_km_matched': float(np.median(ew_d[ew_idx])) if len(ew_idx) else None,
        'ns_dist_median_km_matched': float(np.median(ns_d[ns_idx])) if len(ns_idx) else None,
        'ew_lat_median_raw': float(np.median(ew_lat)) if ew_lat is not None and len(ew_lat) else None,
        'ns_lat_median_raw': float(np.median(ns_lat)) if ns_lat is not None and len(ns_lat) else None,
        'ew_abs_lat_median_raw': float(np.median(np.abs(ew_lat))) if ew_lat is not None and len(ew_lat) else None,
        'ns_abs_lat_median_raw': float(np.median(np.abs(ns_lat))) if ns_lat is not None and len(ns_lat) else None,
        'ew_north_frac_raw': float(np.mean(ew_lat > 0.0)) if ew_lat is not None and len(ew_lat) else None,
        'ns_north_frac_raw': float(np.mean(ns_lat > 0.0)) if ns_lat is not None and len(ns_lat) else None,
        'ew_lat_median_matched': float(np.median(ew_lat[ew_idx])) if ew_lat is not None and len(ew_idx) else None,
        'ns_lat_median_matched': float(np.median(ns_lat[ns_idx])) if ns_lat is not None and len(ns_idx) else None,
        'ew_abs_lat_median_matched': float(np.median(np.abs(ew_lat[ew_idx]))) if ew_lat is not None and len(ew_idx) else None,
        'ns_abs_lat_median_matched': float(np.median(np.abs(ns_lat[ns_idx]))) if ns_lat is not None and len(ns_idx) else None,
        'ew_north_frac_matched': float(np.mean(ew_lat[ew_idx] > 0.0)) if ew_lat is not None and len(ew_idx) else None,
        'ns_north_frac_matched': float(np.mean(ns_lat[ns_idx] > 0.0)) if ns_lat is not None and len(ns_idx) else None,
        'match_lat': bool(match_lat),
    }


def fit_exponential(distances, coherences, min_pairs=MIN_BIN_COUNT, bins=None, use_sigma=True, min_bins=5):
    """
    Fit exponential decay model.
    
    Args:
        distances: Array of pair distances (km)
        coherences: Array of pair coherence values
        min_pairs: Minimum pairs per bin to include in fit
        bins: Optional custom bin edges (default: 40 log-spaced bins 50-13000km)
        use_sigma: Whether to weight fit by 1/sqrt(N) (default: True). 
                   CODE sector analysis uses False (unweighted) for robustness.
    """
    if bins is None:
        bins = np.logspace(np.log10(MIN_DISTANCE_KM), np.log10(MAX_DISTANCE_KM), N_BINS + 1)
    
    bin_centers, bin_means, bin_counts = [], [], []
    
    for i in range(len(bins) - 1):
        mask = (distances >= bins[i]) & (distances < bins[i+1])
        if mask.sum() >= min_pairs:
            bin_centers.append(np.sqrt(bins[i] * bins[i+1]))
            bin_means.append(coherences[mask].mean())
            bin_counts.append(mask.sum())
    
    if len(bin_centers) < min_bins:
        return None
    
    bin_centers = np.array(bin_centers)
    bin_means = np.array(bin_means)
    bin_counts = np.array(bin_counts)
    
    try:
        weights = np.sqrt(bin_counts)
        
        # Adaptive bounds matching CODE longspan TEPConfig
        max_dist = np.max(bin_centers) if len(bin_centers) > 0 else MAX_DISTANCE_KM
        dist_upper = float(np.max(bins)) if bins is not None and len(bins) else float(max_dist)
        is_restricted = bool(dist_upper < (MAX_DISTANCE_KM * 0.95))
        max_lambda = min(15000, (dist_upper * 3.0) if is_restricted else (max_dist * 0.8))
        min_lambda = float(max(1.0, min(100.0, max_dist * 0.05)))
        if max_lambda <= min_lambda:
            return None
        
        # Bounds: ([min_A, min_lambda, min_offset], [max_A, max_lambda, max_offset])
        bounds = ([0, min_lambda, -1], [5, max_lambda, 1])
        
        # Data-driven initial guess (CODE methodology)
        c_range = bin_means.max() - bin_means.min()
        c_min = bin_means.min()
        lam0 = float(np.clip(max_dist / 2.0, min_lambda, max_lambda))
        p0 = [c_range if c_range > 0 else 0.5, lam0, c_min]
        
        sigma = 1.0/weights if use_sigma else None
        
        popt, pcov = curve_fit(
            exp_decay, bin_centers, bin_means,
            p0=p0,
            sigma=sigma,
            bounds=bounds,
            maxfev=10000
        )
        
        A, lam, C0 = popt
        lam_err = np.sqrt(pcov[1, 1]) if pcov[1, 1] > 0 else 0
        lam_rel_err = float(lam_err / lam) if lam and lam > 0 else None

        lam_hit_lower = bool(lam <= min_lambda * (1.0 + 1e-6))
        lam_hit_upper = bool(lam >= max_lambda * (1.0 - 1e-6))
        lam_hit_bound = bool(lam_hit_lower or lam_hit_upper)
        
        predicted = exp_decay(bin_centers, *popt)
        ss_res = np.sum((bin_means - predicted)**2)
        ss_tot = np.sum((bin_means - np.mean(bin_means))**2)
        r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0
        
        return {
            'lambda_km': float(lam),
            'lambda_err': float(lam_err),
            'lambda_rel_err': lam_rel_err,
            'amplitude': float(A),
            'offset': float(C0),
            'r_squared': float(r2),
            'n_pairs': int(sum(bin_counts)),
            'n_bins': len(bin_centers),
            'bin_mean_min': float(np.min(bin_means)) if len(bin_means) else None,
            'bin_mean_max': float(np.max(bin_means)) if len(bin_means) else None,
            'bin_mean_std': float(np.std(bin_means)) if len(bin_means) else None,
            'lambda_min_bound': float(min_lambda),
            'lambda_max_bound': float(max_lambda),
            'lambda_hit_lower_bound': bool(lam_hit_lower),
            'lambda_hit_upper_bound': bool(lam_hit_upper),
            'lambda_hit_bound': bool(lam_hit_bound),
        }
    except Exception:
        return None


def _qc_fit_result(res, min_r2=0.70, min_bins=3, max_rel_err=0.50, allow_bound_hit=False):
    if not res:
        return False, "NONE"
    n_bins = int(res.get('n_bins', 0) or 0)
    r2 = float(res.get('r_squared', -1.0) if res.get('r_squared') is not None else -1.0)
    rel_err = res.get('lambda_rel_err', None)
    hit = bool(res.get('lambda_hit_bound', False))

    if n_bins < int(min_bins):
        return False, "NBINS"
    if not np.isfinite(r2) or r2 < float(min_r2):
        return False, "R2"
    if (not allow_bound_hit) and hit:
        return False, "BOUND"
    if rel_err is not None and np.isfinite(rel_err) and float(rel_err) > float(max_rel_err):
        return False, "ERR"

    return True, "OK"


def _summarize_geo_subset(distances, lats):
    if distances is None or len(distances) == 0:
        return {
            'n': 0,
            'dist_p10': None,
            'dist_p50': None,
            'dist_p90': None,
            'lat_p50': None,
            'abs_lat_p50': None,
            'north_frac': None,
        }

    d = np.asarray(distances)
    lat = np.asarray(lats) if lats is not None else None
    p10, p50, p90 = np.percentile(d, [10, 50, 90])
    return {
        'n': int(len(d)),
        'dist_p10': float(p10),
        'dist_p50': float(p50),
        'dist_p90': float(p90),
        'lat_p50': float(np.median(lat)) if lat is not None and len(lat) else None,
        'abs_lat_p50': float(np.median(np.abs(lat))) if lat is not None and len(lat) else None,
        'north_frac': float(np.mean(lat > 0.0)) if lat is not None and len(lat) else None,
    }


def add_month_to_pairs(pairs_dict):
    """
    Add month field to pairs dict if it has year and doy but not month.
    Returns the pairs dict (modified in place or new dict).
    """
    if not pairs_dict or not isinstance(pairs_dict, dict):
        return pairs_dict
    
    if 'month' in pairs_dict:
        return pairs_dict  # Already has month
    
    if 'doy' not in pairs_dict or 'year' not in pairs_dict:
        return pairs_dict  # Can't compute month
 
    year_doy = pairs_dict['year'].astype(np.int64) * 1000 + pairs_dict['doy'].astype(np.int64)
    dt = pd.to_datetime(pd.Series(year_doy), format='%Y%j', errors='coerce')
    months = dt.dt.month.fillna(0).astype(np.int32).to_numpy()
    
    # Create new dict with month added
    result = dict(pairs_dict)
    result['month'] = months
    return result


def compute_monthly_anisotropy_from_pairs(pairs, station_coords, short_dist_only=True):
    """
    Compute E-W/N-S anisotropy ratio for each month separately.
    
    This addresses the insight that E-W/N-S ratio varies monthly (CMB-coupled)
    and that aggregating all months masks the true temporal pattern.
    
    Parameters:
    -----------
    pairs : dict with arrays
        Pair results with keys: dist, coherence, azimuth, year, month/doy
        Optionally includes 'phase_alignment' for dual-metric analysis
    station_coords : dict
        Station coordinates {sta_id: {'lat': float, 'lon': float}}
    short_dist_only : bool
        If True, only use pairs < 500km for ratio computation
        
    Returns:
    --------
    dict : {
        'coherence': {year_month_key: {...}, ...},
        'phase_alignment': {year_month_key: {...}, ...},
        'summary': {...}
    }
    """
    # Check if we have data
    if not pairs or not isinstance(pairs, dict):
        return {}
    
    if 'dist' not in pairs or len(pairs.get('dist', [])) == 0:
        return {}
    
    # Ensure month is present
    pairs = add_month_to_pairs(pairs)
    
    if 'year' not in pairs or 'month' not in pairs:
        print_status("  Monthly stratification requires year/month in pair data", "WARNING")
        return {}
    
    years = pairs['year']
    months = pairs['month']
    distances = pairs['dist']
    coherences = pairs['coherence']
    azimuths = pairs['azimuth']
    phase_alignments = pairs.get('phase_alignment')
    
    # Get unique year-months
    year_months = set(zip(years, months))
    
    if len(year_months) < 2:
        return {}
    
    # E-W and N-S sector masks (45° wide centered on cardinals)
    ew_mask_full = ((azimuths >= 67.5) & (azimuths < 112.5)) | ((azimuths >= 247.5) & (azimuths < 292.5))
    ns_mask_full = ((azimuths < 22.5) | (azimuths >= 337.5)) | ((azimuths >= 157.5) & (azimuths < 202.5))
    
    # Short distance mask
    short_mask = distances < 500 if short_dist_only else np.ones(len(distances), dtype=bool)
    
    coh_results = {}
    phase_results = {}
    
    for year, month in sorted(year_months):
        if year == 0 or month == 0:
            continue
            
        # Filter to this month
        month_mask = (years == year) & (months == month) & short_mask
        
        if np.sum(month_mask) < 500:
            continue
        
        # E-W and N-S pairs for this month
        ew_mask = month_mask & ew_mask_full
        ns_mask = month_mask & ns_mask_full
        
        n_ew = np.sum(ew_mask)
        n_ns = np.sum(ns_mask)
        
        if n_ew < 50 or n_ns < 50:
            continue
        
        ym_key = f"{year}-{month:02d}"
        
        # Coherence ratio
        ew_coh_mean = np.mean(coherences[ew_mask])
        ns_coh_mean = np.mean(coherences[ns_mask])
        coh_ratio = ew_coh_mean / ns_coh_mean if ns_coh_mean > 0 else None
        
        coh_results[ym_key] = {
            'year': int(year),
            'month': int(month),
            'ratio': float(coh_ratio) if coh_ratio else None,
            'ew_mean': float(ew_coh_mean),
            'ns_mean': float(ns_coh_mean),
            'n_ew': int(n_ew),
            'n_ns': int(n_ns)
        }
        
        # Phase alignment ratio (if available)
        if phase_alignments is not None:
            ew_phase_mean = np.mean(phase_alignments[ew_mask])
            ns_phase_mean = np.mean(phase_alignments[ns_mask])
            phase_ratio = ew_phase_mean / ns_phase_mean if ns_phase_mean > 0 else None
            
            phase_results[ym_key] = {
                'year': int(year),
                'month': int(month),
                'ratio': float(phase_ratio) if phase_ratio else None,
                'ew_mean': float(ew_phase_mean),
                'ns_mean': float(ns_phase_mean),
                'n_ew': int(n_ew),
                'n_ns': int(n_ns)
            }
    
    # Compute summary statistics
    def summarize(results):
        if not results:
            return {}
        ratios = [r['ratio'] for r in results.values() if r.get('ratio')]
        if not ratios:
            return {}
        ew_dominant = sum(1 for r in ratios if r > 1.0)
        return {
            'n_months': len(ratios),
            'mean_ratio': float(np.mean(ratios)),
            'std_ratio': float(np.std(ratios)),
            'min_ratio': float(np.min(ratios)),
            'max_ratio': float(np.max(ratios)),
            'months_ew_dominant': ew_dominant,
            'pct_ew_dominant': 100.0 * ew_dominant / len(ratios),
            'temporal_modulation_detected': np.std(ratios) > 0.05
        }
    
    return {
        'coherence': coh_results,
        'phase_alignment': phase_results if phase_results else None,
        'summary_coherence': summarize(coh_results),
        'summary_phase_alignment': summarize(phase_results) if phase_results else None
    }


# Keep old function signature for backward compatibility
def compute_monthly_anisotropy_from_pairs_legacy(all_pairs, station_coords):
    """Legacy wrapper for backward compatibility."""
    result = compute_monthly_anisotropy_from_pairs(all_pairs, station_coords)
    if not result:
        return {}
    # Return just coherence results in old format
    return result.get('coherence', {})


# Compatibility alias
_compute_monthly_legacy = compute_monthly_anisotropy_from_pairs_legacy


def analyze_anisotropy(final_pairs, condition_name="ALL"):
    """
    Analyze anisotropy patterns from a set of coherence pairs (Vectorized SOA version).
    
    ALIGNED WITH CODE LONGSPAN METHODOLOGY:
    1. Fit each of 8 sectors separately with distance distribution matching
    2. Use 10 coarse bins (100-10000km) and unweighted fit for sector robustness
    3. Calculate E-W/N-S ratio by AVERAGING λ values from E+W and N+S sectors
    """
    print_status(f"\n--- Analyzing {condition_name} ---", "PROCESS")
    
    if not final_pairs or 'dist' not in final_pairs or len(final_pairs['dist']) == 0:
        print_status("  No valid pairs", "WARNING")
        return None
    
    # Unpack arrays for fast access
    all_distances = final_pairs['dist']
    all_coherences = final_pairs['coherence']
    all_azimuths = final_pairs['azimuth']
    all_lats = final_pairs['mid_lat']
    all_lon_diff = final_pairs.get('lon_diff_deg')
    
    n_total = len(all_distances)
    
    # ==========================================================================
    # SECTOR ANALYSIS WITH DISTANCE DISTRIBUTION MATCHING (CODE LONGSPAN METHOD)
    # ==========================================================================
    sectors = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
    sector_width = 45.0
    sector_results = {}
    sector_results_matched = {}  # With distance distribution matching
    
    # CODE sector analysis uses coarser bins for robustness (10 bins, 100-10000km)
    sector_bins = np.logspace(np.log10(100), np.log10(10000), 11) # 10 bins = 11 edges
    
    for i, sector in enumerate(sectors):
        center_az = i * sector_width
        min_az = center_az - sector_width/2
        max_az = center_az + sector_width/2
        
        # Vectorized sector masking
        if min_az < 0:
            # Wrap around 0 (North)
            sector_mask = (all_azimuths >= (360 + min_az)) | (all_azimuths < max_az)
        else:
            sector_mask = (all_azimuths >= min_az) & (all_azimuths < max_az)
            
        if np.sum(sector_mask) < 50:
            continue
        
        sector_distances = all_distances[sector_mask]
        sector_coherences = all_coherences[sector_mask]
        
        # Fit without matching (using robust sector settings)
        res_raw = fit_exponential(
            sector_distances, sector_coherences, 
            bins=sector_bins, use_sigma=False
        )
        if res_raw:
            sector_results[sector] = res_raw
        
        # Apply distance distribution matching (CODE longspan method)
        matched_indices = _subsample_to_match_distribution(
            sector_distances, all_distances, max_samples=min(50000, len(sector_distances))
        )
        
        if len(matched_indices) >= 50:
            res_matched = fit_exponential(
                sector_distances[matched_indices],
                sector_coherences[matched_indices],
                bins=sector_bins, 
                use_sigma=False  # CODE uses unweighted fit for matched sectors
            )
            if res_matched:
                sector_results_matched[sector] = res_matched
    
    # ==========================================================================
    # PRINT ALL 8 SECTOR λ VALUES (MATCHED) - IDENTICAL TO CODE LONGSPAN OUTPUT
    # ==========================================================================
    print_status("  8-SECTOR ANALYSIS (with distance matching):", "INFO")
    for sector in sectors:
        if sector in sector_results_matched:
            res = sector_results_matched[sector]
            print_status(f"    {sector:2s}: λ = {res['lambda_km']:5.0f} km, R² = {res['r_squared']:.3f}", "INFO")
    
    # ==========================================================================
    # DIPOLE ANALYSIS - FIND STRONGEST/WEAKEST DIRECTIONS (CODE LONGSPAN METHOD)
    # ==========================================================================
    dipole_analysis = {}
    if len(sector_results_matched) >= 4:
        lambda_by_sector = {s: sector_results_matched[s]['lambda_km'] for s in sector_results_matched}
        lambda_values = list(lambda_by_sector.values())
        
        max_lambda = max(lambda_values)
        min_lambda = min(lambda_values)
        max_sector = [k for k, v in lambda_by_sector.items() if v == max_lambda][0]
        min_sector = [k for k, v in lambda_by_sector.items() if v == min_lambda][0]
        
        dipole_ratio = max_lambda / min_lambda if min_lambda > 0 else None
        dipole_strength = (max_lambda - min_lambda) / np.mean(lambda_values) if np.mean(lambda_values) > 0 else 0
        
        # Calculate dipole axis angle (sector center azimuth)
        sector_azimuths = {'N': 0, 'NE': 45, 'E': 90, 'SE': 135, 'S': 180, 'SW': 225, 'W': 270, 'NW': 315}
        dipole_axis_az = sector_azimuths.get(max_sector, 0)
        
        dipole_analysis = {
            'strongest_direction': max_sector,
            'strongest_lambda': float(max_lambda),
            'weakest_direction': min_sector,
            'weakest_lambda': float(min_lambda),
            'dipole_ratio': float(dipole_ratio) if dipole_ratio else None,
            'dipole_strength': float(dipole_strength),
            'dipole_axis_azimuth': dipole_axis_az
        }
        
        print_status(f"  DIPOLE ANALYSIS:", "INFO")
        print_status(f"    Strongest: {max_sector} (λ = {max_lambda:.0f} km)", "SUCCESS")
        print_status(f"    Weakest:   {min_sector} (λ = {min_lambda:.0f} km)", "INFO")
        print_status(f"    Dipole ratio (max/min): {dipole_ratio:.2f}", "INFO")
        print_status(f"    Dipole axis: {dipole_axis_az}° from North", "INFO")
    
    # ==========================================================================
    # E-W/N-S RATIO USING FULL DATASET (ROBUST METHOD FOR SPP)
    # ==========================================================================
    ew_sectors = ['E', 'W']
    ns_sectors = ['N', 'S']
    
    # Use FULL DATA results (sector_results) instead of subsampled (sector_results_matched)
    # For noisy SPP data, using all 178M points is statistically superior to subsampling 50k.
    # The binning in fit_exponential handles density differences.
    primary_results = sector_results
    
    ew_lambdas = [primary_results[s]['lambda_km'] for s in ew_sectors if s in primary_results]
    ns_lambdas = [primary_results[s]['lambda_km'] for s in ns_sectors if s in primary_results]
    
    ew_ns_ratio = None
    ew_lambda_mean = None
    ns_lambda_mean = None
    if len(ew_lambdas) >= 1 and len(ns_lambdas) >= 1:
        ew_lambda_mean = np.mean(ew_lambdas)
        ns_lambda_mean = np.mean(ns_lambdas)
        ew_ns_ratio = ew_lambda_mean / ns_lambda_mean if ns_lambda_mean > 0 else None
        print_status(f"  E-W/N-S RATIO (FULL DATASET - PRIMARY):", "SUCCESS")
        ew_str = ", ".join([f"{l:.0f}" for l in ew_lambdas])
        ns_str = ", ".join([f"{l:.0f}" for l in ns_lambdas])
        print_status(f"    E-W λ = {ew_lambda_mean:.0f} km (sectors: {ew_str})", "INFO")
        print_status(f"    N-S λ = {ns_lambda_mean:.0f} km (sectors: {ns_str})", "INFO")
        print_status(f"    RATIO: {ew_ns_ratio:.2f} (Ref: 2.16)", "SUCCESS")
        
    # Secondary: CODE Matched Method
    ew_lambdas_matched = []
    ns_lambdas_matched = []
    ratio_matched = None
    
    if sector_results_matched:
        ew_matched = [sector_results_matched[s]['lambda_km'] for s in ew_sectors if s in sector_results_matched]
        ns_matched = [sector_results_matched[s]['lambda_km'] for s in ns_sectors if s in sector_results_matched]
        if ew_matched and ns_matched:
             ratio_matched = np.mean(ew_matched) / np.mean(ns_matched)
             print_status(f"  [Reference] CODE Subsampled Ratio: {ratio_matched:.2f}", "INFO")
             ew_lambdas_matched = ew_matched
             ns_lambdas_matched = ns_matched
    
    # ==========================================================================
    # POOLED METHOD (SECONDARY / DIAGNOSTIC)
    # ==========================================================================
    # Vectorized masking for pooled method
    # E: [67.5, 112.5) OR W: [247.5, 292.5)
    ew_mask = ((all_azimuths >= 67.5) & (all_azimuths < 112.5)) | \
              ((all_azimuths >= 247.5) & (all_azimuths < 292.5))
              
    # N: [337.5, 360) U [0, 22.5) OR S: [157.5, 202.5)
    ns_mask = (all_azimuths >= 337.5) | (all_azimuths < 22.5) | \
              ((all_azimuths >= 157.5) & (all_azimuths < 202.5))
    
    ew_count = np.sum(ew_mask)
    ns_count = np.sum(ns_mask)
    # print_status(f"  E-W pairs: {ew_count}, N-S pairs: {ns_count}", "INFO")
    
    ew_result = fit_exponential(all_distances[ew_mask], all_coherences[ew_mask]) if ew_count > 100 else None
    ns_result = fit_exponential(all_distances[ns_mask], all_coherences[ns_mask]) if ns_count > 100 else None
    
    ew_ns_ratio_pooled = None
    if ew_result and ns_result:
        ew_ns_ratio_pooled = ew_result['lambda_km'] / ns_result['lambda_km']
        # Demote to info/debug - this metric is biased by geometry
        print_status(f"  [Diagnostic] Pooled Ratio: {ew_ns_ratio_pooled:.2f} (Biased by geometry)", "INFO")

    # Undirected pooled diagnostic: fold azimuth modulo 180° to remove station-order directionality
    az_u = np.mod(all_azimuths, 180.0)
    ew_mask_u = (az_u >= 67.5) & (az_u < 112.5)
    ns_mask_u = (az_u < 22.5) | (az_u >= 157.5)
    ew_count_u = int(np.sum(ew_mask_u))
    ns_count_u = int(np.sum(ns_mask_u))

    ew_result_u = fit_exponential(all_distances[ew_mask_u], all_coherences[ew_mask_u]) if ew_count_u > 100 else None
    ns_result_u = fit_exponential(all_distances[ns_mask_u], all_coherences[ns_mask_u]) if ns_count_u > 100 else None

    ew_ns_ratio_pooled_undirected = None
    if ew_result_u and ns_result_u:
        ew_ns_ratio_pooled_undirected = ew_result_u['lambda_km'] / ns_result_u['lambda_km']
        print_status(
            f"  [Diagnostic] Undirected Pooled Ratio (az%180): {ew_ns_ratio_pooled_undirected:.2f} (order-invariant)",
            "INFO"
        )

    lon_diff_stratified = {}
    lon_diff_stratified_latmatched = {}
    lon_diff_stratified_diagnostics = {}
    if all_lon_diff is not None and np.sum(np.isfinite(all_lon_diff)) > 1000:
        print_status("\n  LONGITUDE DIFFERENCE STRATIFICATION (local-time control):", "INFO")
        print_status("  " + "-"*75, "INFO")
        print_status(f"  {'|Δlon| bin (deg)':<18} {'N(EW)':>10} {'N(NS)':>10} {'Nm(EW)':>10} {'Nm(NS)':>10} {'λ(EW)':>10} {'λ(NS)':>10} {'R²(EW)':>8} {'R²(NS)':>8} {'B':>3} {'Ratio':>10}", "INFO")
        print_status("  " + "-"*75, "INFO")

        lon_bins = [(0.0, 10.0), (10.0, 20.0), (20.0, 40.0), (40.0, 180.0)]
        lat_bins = np.array([-90.0, -60.0, -30.0, 0.0, 30.0, 60.0, 90.0])

        e_mask = ((all_azimuths >= 67.5) & (all_azimuths < 112.5))
        w_mask = ((all_azimuths >= 247.5) & (all_azimuths < 292.5))
        n_mask_dir = (all_azimuths >= 337.5) | (all_azimuths < 22.5)
        s_mask_dir = ((all_azimuths >= 157.5) & (all_azimuths < 202.5))

        ew_mask_dir = (e_mask | w_mask)
        ns_mask_dir = (n_mask_dir | s_mask_dir)

        finite_lon = np.isfinite(all_lon_diff)
        for lo, hi in lon_bins:
            bin_mask = finite_lon & (all_lon_diff >= lo) & (all_lon_diff < hi)
            n_ew = int(np.sum(bin_mask & (e_mask | w_mask)))
            n_ns = int(np.sum(bin_mask & (n_mask_dir | s_mask_dir)))

            ew_res_bin, ns_res_bin, match_meta = _fit_matched_directional_lambda(
                all_distances,
                all_coherences,
                ew_mask=(bin_mask & ew_mask_dir),
                ns_mask=(bin_mask & ns_mask_dir),
                reference_distances=all_distances[bin_mask],
                lats=all_lats,
                reference_lats=all_lats[bin_mask],
                bins=None,
                min_pairs=20,
                use_sigma=False,
                min_bins=3,
                max_samples=50000,
            )

            ew_ok, ew_qc = _qc_fit_result(ew_res_bin, min_r2=0.70, min_bins=3, max_rel_err=0.50, allow_bound_hit=False)
            ns_ok, ns_qc = _qc_fit_result(ns_res_bin, min_r2=0.70, min_bins=3, max_rel_err=0.50, allow_bound_hit=False)

            key = f"{lo:.0f}-{hi:.0f}"

            ew_mean_bin = float(ew_res_bin['lambda_km']) if ew_res_bin else None
            ns_mean_bin = float(ns_res_bin['lambda_km']) if ns_res_bin else None
            ratio_bin = (ew_mean_bin / ns_mean_bin) if (ew_mean_bin is not None and ns_mean_bin and ns_mean_bin > 0) else None

            ew_hit = bool(ew_res_bin.get('lambda_hit_bound')) if ew_res_bin else False
            ns_hit = bool(ns_res_bin.get('lambda_hit_bound')) if ns_res_bin else False
            ratio_unreliable = bool((not ew_ok) or (not ns_ok))

            r2_ew = float(ew_res_bin['r_squared']) if ew_res_bin and 'r_squared' in ew_res_bin else None
            r2_ns = float(ns_res_bin['r_squared']) if ns_res_bin and 'r_squared' in ns_res_bin else None
            n_ew_m = int(match_meta.get('n_ew_matched', 0)) if match_meta else 0
            n_ns_m = int(match_meta.get('n_ns_matched', 0)) if match_meta else 0

            bound_code = ("E" if ew_hit else "") + ("N" if ns_hit else "")
            if bound_code == "":
                bound_code = "."

            lon_diff_stratified[key] = {
                'lon_diff_min_deg': float(lo),
                'lon_diff_max_deg': float(hi),
                'n_ew_pairs': n_ew,
                'n_ns_pairs': n_ns,
                'n_ew_pairs_matched': n_ew_m,
                'n_ns_pairs_matched': n_ns_m,
                'ew_lambda_mean': ew_mean_bin,
                'ns_lambda_mean': ns_mean_bin,
                'ratio': float(ratio_bin) if ratio_bin is not None else None,
                'ratio_unreliable': bool(ratio_unreliable),
                'ratio_qc_ok': bool(ew_ok and ns_ok),
                'ew_fit_qc': str(ew_qc),
                'ns_fit_qc': str(ns_qc),
                'ew_lambda_hit_bound': bool(ew_hit),
                'ns_lambda_hit_bound': bool(ns_hit),
                'ew_fit_r2': r2_ew,
                'ns_fit_r2': r2_ns,
                'ew_fit_n_bins': int(ew_res_bin.get('n_bins')) if ew_res_bin and 'n_bins' in ew_res_bin else None,
                'ns_fit_n_bins': int(ns_res_bin.get('n_bins')) if ns_res_bin and 'n_bins' in ns_res_bin else None,
                'ew_dist_median_km_raw': match_meta.get('ew_dist_median_km_raw') if match_meta else None,
                'ns_dist_median_km_raw': match_meta.get('ns_dist_median_km_raw') if match_meta else None,
                'ew_dist_median_km_matched': match_meta.get('ew_dist_median_km_matched') if match_meta else None,
                'ns_dist_median_km_matched': match_meta.get('ns_dist_median_km_matched') if match_meta else None,
                'ew_abs_lat_median_raw': match_meta.get('ew_abs_lat_median_raw') if match_meta else None,
                'ns_abs_lat_median_raw': match_meta.get('ns_abs_lat_median_raw') if match_meta else None,
                'ew_abs_lat_median_matched': match_meta.get('ew_abs_lat_median_matched') if match_meta else None,
                'ns_abs_lat_median_matched': match_meta.get('ns_abs_lat_median_matched') if match_meta else None,
                'ew_north_frac_raw': match_meta.get('ew_north_frac_raw') if match_meta else None,
                'ns_north_frac_raw': match_meta.get('ns_north_frac_raw') if match_meta else None,
                'ew_north_frac_matched': match_meta.get('ew_north_frac_matched') if match_meta else None,
                'ns_north_frac_matched': match_meta.get('ns_north_frac_matched') if match_meta else None,
            }

            lon_diff_stratified_diagnostics[key] = {
                'ew_raw': _summarize_geo_subset(all_distances[bin_mask & ew_mask_dir], all_lats[bin_mask & ew_mask_dir]),
                'ns_raw': _summarize_geo_subset(all_distances[bin_mask & ns_mask_dir], all_lats[bin_mask & ns_mask_dir]),
            }

            ew_str = f"{ew_mean_bin:.0f}" if ew_mean_bin is not None else "N/A"
            ns_str = f"{ns_mean_bin:.0f}" if ns_mean_bin is not None else "N/A"
            ratio_str = f"{ratio_bin:.2f}{'*' if ratio_unreliable else ''}" if ratio_bin is not None else "N/A"
            r2_ew_str = f"{r2_ew:.3f}" if r2_ew is not None else "N/A"
            r2_ns_str = f"{r2_ns:.3f}" if r2_ns is not None else "N/A"
            print_status(f"  {key:<18} {n_ew:>10,} {n_ns:>10,} {n_ew_m:>10,} {n_ns_m:>10,} {ew_str:>10} {ns_str:>10} {r2_ew_str:>8} {r2_ns_str:>8} {bound_code:>3} {ratio_str:>10}", "INFO")

            if key == "10-20":
                ew_abs_lat_raw = lon_diff_stratified[key].get('ew_abs_lat_median_raw')
                ns_abs_lat_raw = lon_diff_stratified[key].get('ns_abs_lat_median_raw')
                ew_abs_lat_m = lon_diff_stratified[key].get('ew_abs_lat_median_matched')
                ns_abs_lat_m = lon_diff_stratified[key].get('ns_abs_lat_median_matched')
                print_status(
                    f"    [Diagnostic] |Δlon|=10–20°: QC(EW)={ew_qc} QC(NS)={ns_qc}  "
                    f"|lat|50(raw) EW/NS={ew_abs_lat_raw if ew_abs_lat_raw is not None else float('nan'):.1f}/"
                    f"{ns_abs_lat_raw if ns_abs_lat_raw is not None else float('nan'):.1f}  "
                    f"|lat|50(matched) EW/NS={ew_abs_lat_m if ew_abs_lat_m is not None else float('nan'):.1f}/"
                    f"{ns_abs_lat_m if ns_abs_lat_m is not None else float('nan'):.1f}",
                    "INFO"
                )

        print_status("\n  LONGITUDE DIFFERENCE STRATIFICATION (distance+latitude matching):", "INFO")
        print_status("  " + "-"*75, "INFO")
        print_status(f"  {'|Δlon| bin (deg)':<18} {'N(EW)':>10} {'N(NS)':>10} {'Nm(EW)':>10} {'Nm(NS)':>10} {'λ(EW)':>10} {'λ(NS)':>10} {'R²(EW)':>8} {'R²(NS)':>8} {'B':>3} {'Ratio':>10}", "INFO")
        print_status("  " + "-"*75, "INFO")

        for lo, hi in lon_bins:
            bin_mask = finite_lon & (all_lon_diff >= lo) & (all_lon_diff < hi)
            n_ew = int(np.sum(bin_mask & (e_mask | w_mask)))
            n_ns = int(np.sum(bin_mask & (n_mask_dir | s_mask_dir)))

            ew_res_bin, ns_res_bin, match_meta = _fit_matched_directional_lambda(
                all_distances,
                all_coherences,
                ew_mask=(bin_mask & ew_mask_dir),
                ns_mask=(bin_mask & ns_mask_dir),
                reference_distances=all_distances[bin_mask],
                lats=all_lats,
                reference_lats=all_lats[bin_mask],
                bins=None,
                min_pairs=20,
                use_sigma=False,
                min_bins=3,
                max_samples=50000,
                match_lat=True,
                lat_bins=lat_bins,
            )

            ew_ok, ew_qc = _qc_fit_result(ew_res_bin, min_r2=0.70, min_bins=3, max_rel_err=0.50, allow_bound_hit=False)
            ns_ok, ns_qc = _qc_fit_result(ns_res_bin, min_r2=0.70, min_bins=3, max_rel_err=0.50, allow_bound_hit=False)

            key = f"{lo:.0f}-{hi:.0f}"

            ew_mean_bin = float(ew_res_bin['lambda_km']) if ew_res_bin else None
            ns_mean_bin = float(ns_res_bin['lambda_km']) if ns_res_bin else None
            ratio_bin = (ew_mean_bin / ns_mean_bin) if (ew_mean_bin is not None and ns_mean_bin and ns_mean_bin > 0) else None

            ew_hit = bool(ew_res_bin.get('lambda_hit_bound')) if ew_res_bin else False
            ns_hit = bool(ns_res_bin.get('lambda_hit_bound')) if ns_res_bin else False
            ratio_unreliable = bool((not ew_ok) or (not ns_ok))

            r2_ew = float(ew_res_bin['r_squared']) if ew_res_bin and 'r_squared' in ew_res_bin else None
            r2_ns = float(ns_res_bin['r_squared']) if ns_res_bin and 'r_squared' in ns_res_bin else None
            n_ew_m = int(match_meta.get('n_ew_matched', 0)) if match_meta else 0
            n_ns_m = int(match_meta.get('n_ns_matched', 0)) if match_meta else 0

            bound_code = ("E" if ew_hit else "") + ("N" if ns_hit else "")
            if bound_code == "":
                bound_code = "."

            lon_diff_stratified_latmatched[key] = {
                'lon_diff_min_deg': float(lo),
                'lon_diff_max_deg': float(hi),
                'match_lat': True,
                'n_ew_pairs': n_ew,
                'n_ns_pairs': n_ns,
                'n_ew_pairs_matched': n_ew_m,
                'n_ns_pairs_matched': n_ns_m,
                'ew_lambda_mean': ew_mean_bin,
                'ns_lambda_mean': ns_mean_bin,
                'ratio': float(ratio_bin) if ratio_bin is not None else None,
                'ratio_unreliable': bool(ratio_unreliable),
                'ratio_qc_ok': bool(ew_ok and ns_ok),
                'ew_fit_qc': str(ew_qc),
                'ns_fit_qc': str(ns_qc),
                'ew_lambda_hit_bound': bool(ew_hit),
                'ns_lambda_hit_bound': bool(ns_hit),
                'ew_fit_r2': r2_ew,
                'ns_fit_r2': r2_ns,
                'ew_fit_n_bins': int(ew_res_bin.get('n_bins')) if ew_res_bin and 'n_bins' in ew_res_bin else None,
                'ns_fit_n_bins': int(ns_res_bin.get('n_bins')) if ns_res_bin and 'n_bins' in ns_res_bin else None,
                'ew_dist_median_km_raw': match_meta.get('ew_dist_median_km_raw') if match_meta else None,
                'ns_dist_median_km_raw': match_meta.get('ns_dist_median_km_raw') if match_meta else None,
                'ew_dist_median_km_matched': match_meta.get('ew_dist_median_km_matched') if match_meta else None,
                'ns_dist_median_km_matched': match_meta.get('ns_dist_median_km_matched') if match_meta else None,
                'ew_abs_lat_median_raw': match_meta.get('ew_abs_lat_median_raw') if match_meta else None,
                'ns_abs_lat_median_raw': match_meta.get('ns_abs_lat_median_raw') if match_meta else None,
                'ew_abs_lat_median_matched': match_meta.get('ew_abs_lat_median_matched') if match_meta else None,
                'ns_abs_lat_median_matched': match_meta.get('ns_abs_lat_median_matched') if match_meta else None,
                'ew_north_frac_raw': match_meta.get('ew_north_frac_raw') if match_meta else None,
                'ns_north_frac_raw': match_meta.get('ns_north_frac_raw') if match_meta else None,
                'ew_north_frac_matched': match_meta.get('ew_north_frac_matched') if match_meta else None,
                'ns_north_frac_matched': match_meta.get('ns_north_frac_matched') if match_meta else None,
            }

            ew_str = f"{ew_mean_bin:.0f}" if ew_mean_bin is not None else "N/A"
            ns_str = f"{ns_mean_bin:.0f}" if ns_mean_bin is not None else "N/A"
            ratio_str = f"{ratio_bin:.2f}{'*' if ratio_unreliable else ''}" if ratio_bin is not None else "N/A"
            r2_ew_str = f"{r2_ew:.3f}" if r2_ew is not None else "N/A"
            r2_ns_str = f"{r2_ns:.3f}" if r2_ns is not None else "N/A"
            print_status(f"  {key:<18} {n_ew:>10,} {n_ns:>10,} {n_ew_m:>10,} {n_ns_m:>10,} {ew_str:>10} {ns_str:>10} {r2_ew_str:>8} {r2_ns_str:>8} {bound_code:>3} {ratio_str:>10}", "INFO")

    lambda_distance_regimes = {}
    lambda_distance_regimes_latmatched = {}
    print_status("\n  TWO-REGIME λ FITS (avoid mixing distance physics):", "INFO")
    print_status("  " + "-"*75, "INFO")
    print_status(f"  {'Distance band':<18} {'N(EW)':>10} {'N(NS)':>10} {'Nm(EW)':>10} {'Nm(NS)':>10} {'λ(EW)':>10} {'λ(NS)':>10} {'R²(EW)':>8} {'R²(NS)':>8} {'B':>3} {'Ratio':>10}", "INFO")
    print_status("  " + "-"*75, "INFO")

    regimes = [(50.0, 500.0, '50-500 km'), (500.0, 3000.0, '500-3000 km')]
    for dmin, dmax, label in regimes:
        band_mask = (all_distances >= dmin) & (all_distances < dmax)
        n_band = int(np.sum(band_mask))

        bins = np.logspace(np.log10(dmin), np.log10(dmax), 11)

        e_mask = ((all_azimuths >= 67.5) & (all_azimuths < 112.5))
        w_mask = ((all_azimuths >= 247.5) & (all_azimuths < 292.5))
        n_mask_dir = (all_azimuths >= 337.5) | (all_azimuths < 22.5)
        s_mask_dir = ((all_azimuths >= 157.5) & (all_azimuths < 202.5))

        ew_mask_dir = (e_mask | w_mask)
        ns_mask_dir = (n_mask_dir | s_mask_dir)

        n_ew = int(np.sum(band_mask & ew_mask_dir))
        n_ns = int(np.sum(band_mask & ns_mask_dir))

        ew_res_band, ns_res_band, match_meta = _fit_matched_directional_lambda(
            all_distances,
            all_coherences,
            ew_mask=(band_mask & ew_mask_dir),
            ns_mask=(band_mask & ns_mask_dir),
            reference_distances=all_distances[band_mask],
            lats=all_lats,
            reference_lats=all_lats[band_mask],
            bins=bins,
            min_pairs=20,
            use_sigma=False,
            min_bins=3,
            max_samples=50000,
        )

        ew_mean_band = float(ew_res_band['lambda_km']) if ew_res_band else None
        ns_mean_band = float(ns_res_band['lambda_km']) if ns_res_band else None
        ratio_band = (ew_mean_band / ns_mean_band) if (ew_mean_band is not None and ns_mean_band and ns_mean_band > 0) else None

        ew_ok, ew_qc = _qc_fit_result(ew_res_band, min_r2=0.70, min_bins=3, max_rel_err=0.50, allow_bound_hit=False)
        ns_ok, ns_qc = _qc_fit_result(ns_res_band, min_r2=0.70, min_bins=3, max_rel_err=0.50, allow_bound_hit=False)

        ew_hit = bool(ew_res_band.get('lambda_hit_bound')) if ew_res_band else False
        ns_hit = bool(ns_res_band.get('lambda_hit_bound')) if ns_res_band else False
        ratio_unreliable = bool((not ew_ok) or (not ns_ok))

        r2_ew = float(ew_res_band['r_squared']) if ew_res_band and 'r_squared' in ew_res_band else None
        r2_ns = float(ns_res_band['r_squared']) if ns_res_band and 'r_squared' in ns_res_band else None
        n_ew_m = int(match_meta.get('n_ew_matched', 0)) if match_meta else 0
        n_ns_m = int(match_meta.get('n_ns_matched', 0)) if match_meta else 0

        bound_code = ("E" if ew_hit else "") + ("N" if ns_hit else "")
        if bound_code == "":
            bound_code = "."

        if n_band >= 1000:
            lambda_distance_regimes[label] = {
                'distance_min_km': float(dmin),
                'distance_max_km': float(dmax),
                'n_ew_pairs': n_ew,
                'n_ns_pairs': n_ns,
                'n_ew_pairs_matched': n_ew_m,
                'n_ns_pairs_matched': n_ns_m,
                'ew_lambda_mean': ew_mean_band,
                'ns_lambda_mean': ns_mean_band,
                'ratio': float(ratio_band) if ratio_band is not None else None,
                'ratio_unreliable': bool(ratio_unreliable),
                'ratio_qc_ok': bool(ew_ok and ns_ok),
                'ew_fit_qc': str(ew_qc),
                'ns_fit_qc': str(ns_qc),
                'ew_lambda_hit_bound': bool(ew_hit),
                'ns_lambda_hit_bound': bool(ns_hit),
                'ew_fit_r2': r2_ew,
                'ns_fit_r2': r2_ns,
                'ew_fit_n_bins': int(ew_res_band.get('n_bins')) if ew_res_band and 'n_bins' in ew_res_band else None,
                'ns_fit_n_bins': int(ns_res_band.get('n_bins')) if ns_res_band and 'n_bins' in ns_res_band else None,
                'ew_dist_median_km_raw': match_meta.get('ew_dist_median_km_raw') if match_meta else None,
                'ns_dist_median_km_raw': match_meta.get('ns_dist_median_km_raw') if match_meta else None,
                'ew_dist_median_km_matched': match_meta.get('ew_dist_median_km_matched') if match_meta else None,
                'ns_dist_median_km_matched': match_meta.get('ns_dist_median_km_matched') if match_meta else None,
                'ew_abs_lat_median_raw': match_meta.get('ew_abs_lat_median_raw') if match_meta else None,
                'ns_abs_lat_median_raw': match_meta.get('ns_abs_lat_median_raw') if match_meta else None,
                'ew_abs_lat_median_matched': match_meta.get('ew_abs_lat_median_matched') if match_meta else None,
                'ns_abs_lat_median_matched': match_meta.get('ns_abs_lat_median_matched') if match_meta else None,
            }

        ew_str = f"{ew_mean_band:.0f}" if ew_mean_band is not None else "N/A"
        ns_str = f"{ns_mean_band:.0f}" if ns_mean_band is not None else "N/A"
        ratio_str = f"{ratio_band:.2f}{'*' if ratio_unreliable else ''}" if ratio_band is not None else "N/A"
        r2_ew_str = f"{r2_ew:.3f}" if r2_ew is not None else "N/A"
        r2_ns_str = f"{r2_ns:.3f}" if r2_ns is not None else "N/A"
        print_status(f"  {label:<18} {n_ew:>10,} {n_ns:>10,} {n_ew_m:>10,} {n_ns_m:>10,} {ew_str:>10} {ns_str:>10} {r2_ew_str:>8} {r2_ns_str:>8} {bound_code:>3} {ratio_str:>10}", "INFO")

    print_status("\n  TWO-REGIME λ FITS (distance+latitude matching):", "INFO")
    print_status("  " + "-"*75, "INFO")
    print_status(f"  {'Distance band':<18} {'N(EW)':>10} {'N(NS)':>10} {'Nm(EW)':>10} {'Nm(NS)':>10} {'λ(EW)':>10} {'λ(NS)':>10} {'R²(EW)':>8} {'R²(NS)':>8} {'B':>3} {'Ratio':>10}", "INFO")
    print_status("  " + "-"*75, "INFO")

    for dmin, dmax, label in regimes:
        band_mask = (all_distances >= dmin) & (all_distances < dmax)
        bins = np.logspace(np.log10(dmin), np.log10(dmax), 11)

        e_mask = ((all_azimuths >= 67.5) & (all_azimuths < 112.5))
        w_mask = ((all_azimuths >= 247.5) & (all_azimuths < 292.5))
        n_mask_dir = (all_azimuths >= 337.5) | (all_azimuths < 22.5)
        s_mask_dir = ((all_azimuths >= 157.5) & (all_azimuths < 202.5))

        ew_mask_dir = (e_mask | w_mask)
        ns_mask_dir = (n_mask_dir | s_mask_dir)

        n_ew = int(np.sum(band_mask & ew_mask_dir))
        n_ns = int(np.sum(band_mask & ns_mask_dir))

        ew_res_band, ns_res_band, match_meta = _fit_matched_directional_lambda(
            all_distances,
            all_coherences,
            ew_mask=(band_mask & ew_mask_dir),
            ns_mask=(band_mask & ns_mask_dir),
            reference_distances=all_distances[band_mask],
            lats=all_lats,
            reference_lats=all_lats[band_mask],
            bins=bins,
            min_pairs=20,
            use_sigma=False,
            min_bins=3,
            max_samples=50000,
            match_lat=True,
            lat_bins=np.array([-90.0, -60.0, -30.0, 0.0, 30.0, 60.0, 90.0]),
        )

        ew_mean_band = float(ew_res_band['lambda_km']) if ew_res_band else None
        ns_mean_band = float(ns_res_band['lambda_km']) if ns_res_band else None
        ratio_band = (ew_mean_band / ns_mean_band) if (ew_mean_band is not None and ns_mean_band and ns_mean_band > 0) else None

        ew_ok, ew_qc = _qc_fit_result(ew_res_band, min_r2=0.70, min_bins=3, max_rel_err=0.50, allow_bound_hit=False)
        ns_ok, ns_qc = _qc_fit_result(ns_res_band, min_r2=0.70, min_bins=3, max_rel_err=0.50, allow_bound_hit=False)

        ew_hit = bool(ew_res_band.get('lambda_hit_bound')) if ew_res_band else False
        ns_hit = bool(ns_res_band.get('lambda_hit_bound')) if ns_res_band else False
        ratio_unreliable = bool((not ew_ok) or (not ns_ok))

        r2_ew = float(ew_res_band['r_squared']) if ew_res_band and 'r_squared' in ew_res_band else None
        r2_ns = float(ns_res_band['r_squared']) if ns_res_band and 'r_squared' in ns_res_band else None
        n_ew_m = int(match_meta.get('n_ew_matched', 0)) if match_meta else 0
        n_ns_m = int(match_meta.get('n_ns_matched', 0)) if match_meta else 0

        bound_code = ("E" if ew_hit else "") + ("N" if ns_hit else "")
        if bound_code == "":
            bound_code = "."

        lambda_distance_regimes_latmatched[label] = {
            'distance_min_km': float(dmin),
            'distance_max_km': float(dmax),
            'match_lat': True,
            'n_ew_pairs': n_ew,
            'n_ns_pairs': n_ns,
            'n_ew_pairs_matched': n_ew_m,
            'n_ns_pairs_matched': n_ns_m,
            'ew_lambda_mean': ew_mean_band,
            'ns_lambda_mean': ns_mean_band,
            'ratio': float(ratio_band) if ratio_band is not None else None,
            'ratio_unreliable': bool(ratio_unreliable),
            'ratio_qc_ok': bool(ew_ok and ns_ok),
            'ew_fit_qc': str(ew_qc),
            'ns_fit_qc': str(ns_qc),
            'ew_lambda_hit_bound': bool(ew_hit),
            'ns_lambda_hit_bound': bool(ns_hit),
            'ew_fit_r2': r2_ew,
            'ns_fit_r2': r2_ns,
            'ew_fit_n_bins': int(ew_res_band.get('n_bins')) if ew_res_band and 'n_bins' in ew_res_band else None,
            'ns_fit_n_bins': int(ns_res_band.get('n_bins')) if ns_res_band and 'n_bins' in ns_res_band else None,
            'ew_dist_median_km_raw': match_meta.get('ew_dist_median_km_raw') if match_meta else None,
            'ns_dist_median_km_raw': match_meta.get('ns_dist_median_km_raw') if match_meta else None,
            'ew_dist_median_km_matched': match_meta.get('ew_dist_median_km_matched') if match_meta else None,
            'ns_dist_median_km_matched': match_meta.get('ns_dist_median_km_matched') if match_meta else None,
            'ew_abs_lat_median_raw': match_meta.get('ew_abs_lat_median_raw') if match_meta else None,
            'ns_abs_lat_median_raw': match_meta.get('ns_abs_lat_median_raw') if match_meta else None,
            'ew_abs_lat_median_matched': match_meta.get('ew_abs_lat_median_matched') if match_meta else None,
            'ns_abs_lat_median_matched': match_meta.get('ns_abs_lat_median_matched') if match_meta else None,
        }

        ew_str = f"{ew_mean_band:.0f}" if ew_mean_band is not None else "N/A"
        ns_str = f"{ns_mean_band:.0f}" if ns_mean_band is not None else "N/A"
        ratio_str = f"{ratio_band:.2f}{'*' if ratio_unreliable else ''}" if ratio_band is not None else "N/A"
        r2_ew_str = f"{r2_ew:.3f}" if r2_ew is not None else "N/A"
        r2_ns_str = f"{r2_ns:.3f}" if r2_ns is not None else "N/A"
        print_status(f"  {label:<18} {n_ew:>10,} {n_ns:>10,} {n_ew_m:>10,} {n_ns_m:>10,} {ew_str:>10} {ns_str:>10} {r2_ew_str:>8} {r2_ns_str:>8} {bound_code:>3} {ratio_str:>10}", "INFO")

    preferred_estimator = None
    preferred_estimators = {}
    preferred_estimator_selection = None

    if all_lon_diff is not None and np.sum(np.isfinite(all_lon_diff)) > 1000:
        finite_lon = np.isfinite(all_lon_diff)
        e_mask = ((all_azimuths >= 67.5) & (all_azimuths < 112.5))
        w_mask = ((all_azimuths >= 247.5) & (all_azimuths < 292.5))
        n_mask_dir = (all_azimuths >= 337.5) | (all_azimuths < 22.5)
        s_mask_dir = ((all_azimuths >= 157.5) & (all_azimuths < 202.5))

        ew_mask_dir = (e_mask | w_mask)
        ns_mask_dir = (n_mask_dir | s_mask_dir)

        pref_dist_mask = (all_distances >= 500.0) & (all_distances < 3000.0)
        bins = np.logspace(np.log10(500.0), np.log10(3000.0), 11)
        lat_bins = np.array([-90.0, -60.0, -30.0, 0.0, 30.0, 60.0, 90.0])

        def _compute_pref_variant(name, lon_lo, lon_hi, match_lat, extra_mask=None):
            base = finite_lon & (all_lon_diff >= float(lon_lo)) & (all_lon_diff < float(lon_hi)) & pref_dist_mask
            if extra_mask is not None:
                base = base & extra_mask
            if int(np.sum(base)) < 5000:
                return None

            ew_res, ns_res, match_meta = _fit_matched_directional_lambda(
                all_distances,
                all_coherences,
                ew_mask=(base & ew_mask_dir),
                ns_mask=(base & ns_mask_dir),
                reference_distances=all_distances[base],
                lats=all_lats,
                reference_lats=all_lats[base],
                bins=bins,
                min_pairs=20,
                use_sigma=False,
                min_bins=3,
                max_samples=50000,
                match_lat=bool(match_lat),
                lat_bins=lat_bins,
            )

            if not ew_res or not ns_res:
                return None

            ew_ok, ew_qc = _qc_fit_result(ew_res, min_r2=0.70, min_bins=3, max_rel_err=0.50, allow_bound_hit=False)
            ns_ok, ns_qc = _qc_fit_result(ns_res, min_r2=0.70, min_bins=3, max_rel_err=0.50, allow_bound_hit=False)

            ew_lam = float(ew_res['lambda_km'])
            ns_lam = float(ns_res['lambda_km'])
            ratio = (ew_lam / ns_lam) if ns_lam > 0 else None

            ew_hit = bool(ew_res.get('lambda_hit_bound', False))
            ns_hit = bool(ns_res.get('lambda_hit_bound', False))
            bound_code = ("E" if ew_hit else "") + ("N" if ns_hit else "")
            if bound_code == "":
                bound_code = "."

            qc_ok = bool(ew_ok and ns_ok)
            ratio_unreliable = bool(not qc_ok)

            return {
                'name': str(name),
                'lon_diff_min_deg': float(lon_lo),
                'lon_diff_max_deg': float(lon_hi),
                'distance_min_km': 500.0,
                'distance_max_km': 3000.0,
                'match_lat': bool(match_lat),
                'n_total_pairs': int(np.sum(base)),
                'n_ew_pairs': int(np.sum(base & ew_mask_dir)),
                'n_ns_pairs': int(np.sum(base & ns_mask_dir)),
                'n_ew_pairs_matched': int(match_meta.get('n_ew_matched', 0)) if match_meta else 0,
                'n_ns_pairs_matched': int(match_meta.get('n_ns_matched', 0)) if match_meta else 0,
                'ew_lambda': ew_lam,
                'ns_lambda': ns_lam,
                'ratio': float(ratio) if ratio is not None else None,
                'ratio_unreliable': bool(ratio_unreliable),
                'ratio_qc_ok': bool(qc_ok),
                'ew_fit_qc': str(ew_qc),
                'ns_fit_qc': str(ns_qc),
                'ew_fit_r2': float(ew_res.get('r_squared')),
                'ns_fit_r2': float(ns_res.get('r_squared')),
                'ew_lambda_hit_bound': bool(ew_hit),
                'ns_lambda_hit_bound': bool(ns_hit),
                'fit_bins': int(ew_res.get('n_bins')) if 'n_bins' in ew_res else None,
                'ew_dist_median_km_raw': match_meta.get('ew_dist_median_km_raw') if match_meta else None,
                'ns_dist_median_km_raw': match_meta.get('ns_dist_median_km_raw') if match_meta else None,
                'ew_dist_median_km_matched': match_meta.get('ew_dist_median_km_matched') if match_meta else None,
                'ns_dist_median_km_matched': match_meta.get('ns_dist_median_km_matched') if match_meta else None,
                'ew_abs_lat_median_raw': match_meta.get('ew_abs_lat_median_raw') if match_meta else None,
                'ns_abs_lat_median_raw': match_meta.get('ns_abs_lat_median_raw') if match_meta else None,
                'ew_abs_lat_median_matched': match_meta.get('ew_abs_lat_median_matched') if match_meta else None,
                'ns_abs_lat_median_matched': match_meta.get('ns_abs_lat_median_matched') if match_meta else None,
                'ew_north_frac_raw': match_meta.get('ew_north_frac_raw') if match_meta else None,
                'ns_north_frac_raw': match_meta.get('ns_north_frac_raw') if match_meta else None,
                'ew_north_frac_matched': match_meta.get('ew_north_frac_matched') if match_meta else None,
                'ns_north_frac_matched': match_meta.get('ns_north_frac_matched') if match_meta else None,
                'bound_code': str(bound_code),
            }

        abs_lat = np.abs(all_lats)
        mid_lat_only_mask = (abs_lat >= 30.0) & (abs_lat <= 60.0)

        for match_lat in [False, True]:
            suffix = "dist_lat" if match_lat else "dist_only"
            preferred_estimators[f"lon_0_10_{suffix}"] = _compute_pref_variant(f"0-10_{suffix}", 0.0, 10.0, match_lat)
            preferred_estimators[f"lon_10_20_{suffix}"] = _compute_pref_variant(f"10-20_{suffix}", 10.0, 20.0, match_lat)
            preferred_estimators[f"lon_0_20_{suffix}"] = _compute_pref_variant(f"0-20_{suffix}", 0.0, 20.0, match_lat)
            preferred_estimators[f"lon_0_10_{suffix}_midlat"] = _compute_pref_variant(f"0-10_{suffix}_midlat", 0.0, 10.0, match_lat, extra_mask=mid_lat_only_mask)
            preferred_estimators[f"lon_0_20_{suffix}_midlat"] = _compute_pref_variant(f"0-20_{suffix}_midlat", 0.0, 20.0, match_lat, extra_mask=mid_lat_only_mask)

        def _passes(key):
            obj = preferred_estimators.get(key)
            return bool(obj and obj.get('ratio_qc_ok'))

        def _select_and_print():
            candidates = [
                ("lon_0_20_dist_lat", "QC-gated: |Δlon|<20°, 500–3000 km, distance+latitude matched"),
                ("lon_0_10_dist_lat", "QC-gated: |Δlon|<10°, 500–3000 km, distance+latitude matched"),
                ("lon_0_20_dist_only", "QC-gated: |Δlon|<20°, 500–3000 km, distance matched"),
                ("lon_0_10_dist_only", "QC-gated: |Δlon|<10°, 500–3000 km, distance matched"),
            ]

            sub_ok_lat = _passes("lon_0_10_dist_lat") and _passes("lon_10_20_dist_lat")
            sub_ok_dist = _passes("lon_0_10_dist_only") and _passes("lon_10_20_dist_only")
            if not sub_ok_lat:
                if preferred_estimators.get("lon_0_20_dist_lat") is not None:
                    preferred_estimators["lon_0_20_dist_lat"]["subbin_gate_ok"] = False
            else:
                if preferred_estimators.get("lon_0_20_dist_lat") is not None:
                    preferred_estimators["lon_0_20_dist_lat"]["subbin_gate_ok"] = True

            if not sub_ok_dist:
                if preferred_estimators.get("lon_0_20_dist_only") is not None:
                    preferred_estimators["lon_0_20_dist_only"]["subbin_gate_ok"] = False
            else:
                if preferred_estimators.get("lon_0_20_dist_only") is not None:
                    preferred_estimators["lon_0_20_dist_only"]["subbin_gate_ok"] = True

            for key, desc in candidates:
                obj = preferred_estimators.get(key)
                if obj is None:
                    continue
                gate_ok = obj.get("subbin_gate_ok", True)
                if key.startswith("lon_0_20") and not gate_ok:
                    continue
                if obj.get('ratio_qc_ok'):
                    return key, desc
            for key, desc in candidates:
                obj = preferred_estimators.get(key)
                if obj is None:
                    continue
                return key, desc
            return None, None

        chosen_key, chosen_desc = _select_and_print()
        if chosen_key:
            preferred_estimator = preferred_estimators.get(chosen_key)
            preferred_estimator_selection = {
                'chosen_key': str(chosen_key),
                'chosen_definition': str(chosen_desc),
            }

            print_status("\n  PREFERRED ESTIMATORS (QC-gated; for strict comparability):", "SUCCESS")
            for k in ["lon_0_10_dist_only", "lon_10_20_dist_only", "lon_0_20_dist_only", "lon_0_10_dist_lat", "lon_10_20_dist_lat", "lon_0_20_dist_lat"]:
                obj = preferred_estimators.get(k)
                if not obj:
                    continue
                ratio = obj.get('ratio')
                ratio_str = f"{ratio:.2f}" if ratio is not None else "N/A"
                if obj.get('ratio_unreliable'):
                    ratio_str += "*"
                print_status(
                    f"    {k}: λ(EW)={obj['ew_lambda']:.0f}  λ(NS)={obj['ns_lambda']:.0f}  Ratio={ratio_str}  "
                    f"R²(EW)={obj['ew_fit_r2']:.3f}  R²(NS)={obj['ns_fit_r2']:.3f}  QC={obj.get('ew_fit_qc')}/{obj.get('ns_fit_qc')}  B={obj.get('bound_code')}",
                    "INFO"
                )

            ratio = preferred_estimator.get('ratio') if preferred_estimator else None
            ratio_str = f"{ratio:.2f}" if ratio is not None else "N/A"
            if preferred_estimator and preferred_estimator.get('ratio_unreliable'):
                ratio_str += "*"

            print_status("\n  PREFERRED ESTIMATOR (recommended for tracking):", "SUCCESS")
            print_status(f"    Definition: {preferred_estimator_selection['chosen_definition']}", "INFO")
            print_status(
                f"    N(EW)={preferred_estimator['n_ew_pairs']:,}  N(NS)={preferred_estimator['n_ns_pairs']:,}  "
                f"Nm(EW)={preferred_estimator['n_ew_pairs_matched']:,}  Nm(NS)={preferred_estimator['n_ns_pairs_matched']:,}",
                "INFO"
            )
            print_status(
                f"    λ(EW)={preferred_estimator['ew_lambda']:.0f} km  λ(NS)={preferred_estimator['ns_lambda']:.0f} km  Ratio={ratio_str}  "
                f"R²(EW)={preferred_estimator['ew_fit_r2']:.3f}  R²(NS)={preferred_estimator['ns_fit_r2']:.3f}  "
                f"QC={preferred_estimator.get('ew_fit_qc')}/{preferred_estimator.get('ns_fit_qc')}  B={preferred_estimator.get('bound_code')}",
                "INFO"
            )
    
    # Hemisphere Analysis
    n_mask = all_lats > 0
    s_mask = all_lats < 0
    
    n_res = fit_exponential(all_distances[n_mask], all_coherences[n_mask]) if np.sum(n_mask) > 100 else None
    s_res = fit_exponential(all_distances[s_mask], all_coherences[s_mask]) if np.sum(s_mask) > 100 else None
    
    # Latitude Band Analysis (New)
    # Low: < 30, Mid: 30-60, High: > 60 (Absolute latitude)
    abs_lat = np.abs(all_lats)
    low_mask = abs_lat < 30
    mid_mask = (abs_lat >= 30) & (abs_lat <= 60)
    high_mask = abs_lat > 60
    
    lat_results = {}
    mid_lat_ratio = None

    # ==========================================================================
    # SHORT-DISTANCE ANALYSIS (PRIMARY TEP EVIDENCE)
    # ==========================================================================
    # At short distances (<500km), ionospheric local-time decorrelation is minimal,
    # revealing the true TEP signal: E-W correlations should be STRONGER than N-S.
    # This is the most robust test because it avoids the ionospheric masking effect.
    # ==========================================================================
    short_dist_analysis = {}
    short_dist_threshold = 500  # km
    
    short_mask = all_distances < short_dist_threshold
    short_ew_mask = short_mask & ew_mask
    short_ns_mask = short_mask & ns_mask

    short_ew_mask_u = short_mask & ew_mask_u
    short_ns_mask_u = short_mask & ns_mask_u
    
    n_short_ew = np.sum(short_ew_mask)
    n_short_ns = np.sum(short_ns_mask)

    n_short_ew_u = int(np.sum(short_ew_mask_u))
    n_short_ns_u = int(np.sum(short_ns_mask_u))
    
    print_status(f"\n  SHORT-DISTANCE ANALYSIS (<{short_dist_threshold}km) - PRIMARY TEP TEST:", "SUCCESS")
    print_status(f"    E-W pairs: {n_short_ew:,}, N-S pairs: {n_short_ns:,}", "INFO")
    print_status(f"    [Diagnostic] Undirected az%180 pairs: E-W={n_short_ew_u:,}, N-S={n_short_ns_u:,}", "INFO")
    
    if n_short_ew > 100 and n_short_ns > 100:
        # Use raw coherence values (phase_alignment or coherence)
        ew_short_values = all_coherences[short_ew_mask]
        ns_short_values = all_coherences[short_ns_mask]
        
        # CRITICAL AUDIT: Check mean distances to ensure no bias
        ew_short_dists = all_distances[short_ew_mask]
        ns_short_dists = all_distances[short_ns_mask]
        ew_dist_mean = np.mean(ew_short_dists)
        ns_dist_mean = np.mean(ns_short_dists)
        dist_diff_km = ew_dist_mean - ns_dist_mean
        
        print_status(f"    [Audit] Mean Distance: E-W={ew_dist_mean:.1f}km, N-S={ns_dist_mean:.1f}km (Δ={dist_diff_km:.1f}km)", "INFO")
        
        # Perform distance-matched resampling if bias is large (>10km) or just as a robust check
        # We use a histogram-based reweighting or simple bin-matching for robustness
        # Here: simple 50km bin matching for the robust ratio
        
        bins = np.arange(0, short_dist_threshold + 50, 50)
        ew_hist, _ = np.histogram(ew_short_dists, bins)
        ns_hist, _ = np.histogram(ns_short_dists, bins)
        
        matched_ew_vals = []
        matched_ns_vals = []
        
        for i in range(len(bins)-1):
            b_lo, b_hi = bins[i], bins[i+1]
            # Indices in the original arrays (masked)
            # We need indices relative to the short_* arrays
            ew_in_bin = (ew_short_dists >= b_lo) & (ew_short_dists < b_hi)
            ns_in_bin = (ns_short_dists >= b_lo) & (ns_short_dists < b_hi)
            
            n_e = np.sum(ew_in_bin)
            n_n = np.sum(ns_in_bin)
            
            if n_e > 10 and n_n > 10:
                # Take min count to match distributions
                n_take = min(n_e, n_n)
                
                # Random sample to match counts
                e_vals_bin = ew_short_values[ew_in_bin]
                n_vals_bin = ns_short_values[ns_in_bin]
                
                # Deterministic shuffle for reproducibility
                rng = np.random.RandomState(42 + i)
                matched_ew_vals.extend(rng.choice(e_vals_bin, n_take, replace=False))
                matched_ns_vals.extend(rng.choice(n_vals_bin, n_take, replace=False))
        
        if matched_ew_vals:
            matched_ew_mean = np.mean(matched_ew_vals)
            matched_ns_mean = np.mean(matched_ns_vals)
            matched_ratio = matched_ew_mean / matched_ns_mean if matched_ns_mean != 0 else 0
            print_status(f"    [Robust] Distance-Matched (<50km bins) Ratio: {matched_ratio:.3f} (N={len(matched_ew_vals)})", "INFO")
            
            # Update primary ratio if matched is available? 
            # For now, just report it as a strong diagnostic.
        
        ew_short_mean = np.mean(ew_short_values)
        ns_short_mean = np.mean(ns_short_values)
        ew_short_std = np.std(ew_short_values)
        ns_short_std = np.std(ns_short_values)
        
        # Mean ratio (not lambda ratio - direct correlation comparison)
        mean_ratio = ew_short_mean / ns_short_mean if ns_short_mean != 0 else None
        
        # Statistical significance: Welch's t-test
        from scipy import stats
        t_stat, p_value = stats.ttest_ind(ew_short_values, ns_short_values, equal_var=False)
        
        # Bootstrap confidence interval for the ratio
        np.random.seed(42)
        n_bootstrap = 1000
        ratios_bootstrap = []
        for _ in range(n_bootstrap):
            ew_sample = np.random.choice(ew_short_values, size=len(ew_short_values), replace=True)
            ns_sample = np.random.choice(ns_short_values, size=len(ns_short_values), replace=True)
            if np.mean(ns_sample) != 0:
                ratios_bootstrap.append(np.mean(ew_sample) / np.mean(ns_sample))
        
        ci_low, ci_high = np.percentile(ratios_bootstrap, [2.5, 97.5])
        
        # Effect size (Cohen's d)
        pooled_std = np.sqrt((ew_short_std**2 + ns_short_std**2) / 2)
        cohens_d = (ew_short_mean - ns_short_mean) / pooled_std if pooled_std > 0 else 0
        
        # Store results
        short_dist_analysis = {
            'threshold_km': short_dist_threshold,
            'n_ew_pairs': int(n_short_ew),
            'n_ns_pairs': int(n_short_ns),
            'ew_mean': float(ew_short_mean),
            'ns_mean': float(ns_short_mean),
            'ew_std': float(ew_short_std),
            'ns_std': float(ns_short_std),
            'mean_ratio': float(mean_ratio) if mean_ratio else None,
            'ci_95_low': float(ci_low),
            'ci_95_high': float(ci_high),
            't_statistic': float(t_stat),
            'p_value': float(p_value),
            'cohens_d': float(cohens_d),
            'significant': p_value < 0.01
        }
        
        # Print results
        print_status(f"    E-W mean: {ew_short_mean:.4f} ± {ew_short_std:.4f}", "INFO")
        print_status(f"    N-S mean: {ns_short_mean:.4f} ± {ns_short_std:.4f}", "INFO")
        print_status(f"    RATIO (E-W/N-S): {mean_ratio:.3f}", "SUCCESS")
        print_status(f"    95% CI: [{ci_low:.3f}, {ci_high:.3f}]", "INFO")
        print_status(f"    t-statistic: {t_stat:.2f}, p-value: {p_value:.2e}", "INFO")
        print_status(f"    Cohen's d: {cohens_d:.3f} ({'small' if abs(cohens_d) < 0.5 else 'medium' if abs(cohens_d) < 0.8 else 'large'})", "INFO")
        
        # TEP interpretation
        if mean_ratio > 1.0 and ci_low > 1.0:
            print_status(f"    → TEP SIGNAL DETECTED: E-W > N-S (entire 95% CI > 1.0)", "SUCCESS")
            tep_short_dist = "DETECTED"
        elif mean_ratio > 1.0 and p_value < 0.01:
            print_status(f"    → TEP SIGNAL LIKELY: E-W > N-S (p < 0.01)", "SUCCESS")
            tep_short_dist = "LIKELY"
        elif mean_ratio > 1.0:
            print_status(f"    → TEP SIGNAL WEAK: E-W > N-S but not significant", "WARNING")
            tep_short_dist = "WEAK"
        else:
            print_status(f"    → NO TEP SIGNAL: N-S >= E-W at short distances", "WARNING")
            tep_short_dist = "NOT_DETECTED"
        
        short_dist_analysis['tep_status'] = tep_short_dist

    undirected_short_dist = {}
    if n_short_ew_u > 100 and n_short_ns_u > 100:
        ew_u_vals = all_coherences[short_ew_mask_u]
        ns_u_vals = all_coherences[short_ns_mask_u]
        ew_u_mean = float(np.mean(ew_u_vals))
        ns_u_mean = float(np.mean(ns_u_vals))
        ratio_u = (ew_u_mean / ns_u_mean) if ns_u_mean != 0 else None
        undirected_short_dist = {
            'threshold_km': int(short_dist_threshold),
            'n_ew_pairs': int(n_short_ew_u),
            'n_ns_pairs': int(n_short_ns_u),
            'ew_mean': ew_u_mean,
            'ns_mean': ns_u_mean,
            'mean_ratio': float(ratio_u) if ratio_u is not None else None,
            'definition': 'undirected: azimuth folded modulo 180°'
        }
        if ratio_u is not None:
            print_status(f"    [Diagnostic] Undirected short-dist ratio (az%180): {ratio_u:.3f}", "INFO")
        
        # Compare to CODE reference
        print_status(f"    CODE reference ratio: 2.16 (λ ratio, different metric)", "INFO")
    
    # ==========================================================================
    # DISTANCE BAND ANALYSIS
    # ==========================================================================
    dist_band_results = {}
    for dmin, dmax, label in DISTANCE_BANDS:
        band_mask = (all_distances >= dmin) & (all_distances < dmax)
        if np.sum(band_mask) < 200:
            continue
        # E-W / N-S ratio inside this band (pooled method)
        ew_band_mask = band_mask & ew_mask
        ns_band_mask = band_mask & ns_mask
        if np.sum(ew_band_mask) > 100 and np.sum(ns_band_mask) > 100:
            ew_res_band = fit_exponential(all_distances[ew_band_mask], all_coherences[ew_band_mask])
            ns_res_band = fit_exponential(all_distances[ns_band_mask], all_coherences[ns_band_mask])
            if ew_res_band and ns_res_band:
                dist_band_results[label] = {
                    'ew_lambda': ew_res_band['lambda_km'],
                    'ns_lambda': ns_res_band['lambda_km'],
                    'ew_ns_ratio': ew_res_band['lambda_km']/ns_res_band['lambda_km']
                }
    
    # Analyze Mid-Lat E-W/N-S Ratio (cleanest signal, 45° sectors)
    # Combine masks: mid_lat AND ew/ns
    mid_ew_mask = mid_mask & ew_mask
    mid_ns_mask = mid_mask & ns_mask
    
    mid_ew_res = fit_exponential(all_distances[mid_ew_mask], all_coherences[mid_ew_mask]) if np.sum(mid_ew_mask) > 100 else None
    mid_ns_res = fit_exponential(all_distances[mid_ns_mask], all_coherences[mid_ns_mask]) if np.sum(mid_ns_mask) > 100 else None
    
    if mid_ew_res and mid_ns_res:
        mid_lat_ratio = mid_ew_res['lambda_km'] / mid_ns_res['lambda_km']
        lat_results['mid_lat_ratio'] = mid_lat_ratio
        print_status(f"  Mid-Lat Only ({np.sum(mid_mask)} pairs) Ratio: {lat_results['mid_lat_ratio']:.2f}", "INFO")

    # Compute CV
    cv = 0
    if sector_results:
        lambdas = [r['lambda_km'] for r in sector_results.values()]
        cv = np.std(lambdas) / np.mean(lambdas) if np.mean(lambdas) > 0 else 0
    
    return {
        'condition': condition_name,
        'n_pairs': n_total,
        'dipole_analysis': dipole_analysis,
        'ew_ns_code': {
            'ew_lambda_mean': float(ew_lambda_mean) if ew_lambda_mean else None,
            'ns_lambda_mean': float(ns_lambda_mean) if ns_lambda_mean else None,
            'ew_lambdas': [float(x) for x in ew_lambdas],
            'ns_lambdas': [float(x) for x in ns_lambdas],
            'ratio': float(ew_ns_ratio) if ew_ns_ratio else None,
            'method': 'FULL_DATASET_sector_average'
        },
        'ew_ns_pooled': {
            'ew': ew_result,
            'ns': ns_result,
            'ratio': float(ew_ns_ratio_pooled) if ew_ns_ratio_pooled else None,
            'method': 'pooled_fit'
        },
        'ew_ns_pooled_undirected': {
            'ew': ew_result_u,
            'ns': ns_result_u,
            'ratio': float(ew_ns_ratio_pooled_undirected) if ew_ns_ratio_pooled_undirected else None,
            'method': 'pooled_fit_azimuth_folded_mod_180'
        },
        'sectors': sector_results,
        'sectors_matched': sector_results_matched,
        'hemisphere': {
            'north': n_res,
            'south': s_res
        },
        'latitude_bands': {
            'low': fit_exponential(all_distances[low_mask], all_coherences[low_mask]) if np.sum(low_mask) > 100 else None,
            'mid': fit_exponential(all_distances[mid_mask], all_coherences[mid_mask]) if np.sum(mid_mask) > 100 else None,
            'high': fit_exponential(all_distances[high_mask], all_coherences[high_mask]) if np.sum(high_mask) > 100 else None,
            'mid_lat_ratio': mid_lat_ratio,
            'ratios': lat_results
        },
        'distance_bands': dist_band_results,
        'lon_diff_stratified': lon_diff_stratified,
        'lon_diff_stratified_latmatched': lon_diff_stratified_latmatched,
        'lon_diff_stratified_diagnostics': lon_diff_stratified_diagnostics,
        'lambda_distance_regimes': lambda_distance_regimes,
        'lambda_distance_regimes_latmatched': lambda_distance_regimes_latmatched,
        'preferred_estimator': preferred_estimator,
        'preferred_estimators': preferred_estimators,
        'preferred_estimator_selection': preferred_estimator_selection,
        'short_distance_analysis': short_dist_analysis,
        'short_distance_analysis_undirected': undirected_short_dist,
        'anisotropy_cv': float(cv),
        'methodology': {
            'azimuth_calculation': 'spherical_forward_azimuth',
            'sector_width_deg': 45,
            'distance_matching': True,
            'metric_used': 'from_input_pairs',
            'short_dist_threshold_km': 500,
            'frequency_band_hz': [1e-5, 5e-4],
            'reference': 'CODE_longspan_25yr_PPP'
        }
    }


def analyze_tep_evidence(sector_results, metric_name="phase_alignment"):
    """
    Enhanced TEP Evidence Analysis - Compare SPP results to CODE reference.
    
    KEY INSIGHT: The TEP signal manifests as E-W > N-S correlation length.
    In SPP data, orbital geometry creates anisotropic noise that SUPPRESSES E-W.
    By comparing to CODE's PPP results, we can:
    1. Quantify the geometric suppression factor per sector
    2. Identify if N-S values match CODE (signal preserved)
    3. Compute a "corrected" ratio accounting for geometric noise
    
    Returns evidence assessment with statistical confidence.
    """
    if not sector_results or len(sector_results) < 4:
        return None
    
    print_status("\n" + "="*70, "INFO")
    print_status("ENHANCED TEP EVIDENCE ANALYSIS", "SUCCESS")
    print_status("="*70, "INFO")
    print_status(f"  Metric: {metric_name}", "INFO")
    print_status(f"  Method: Compare SPP sector λ to CODE PPP reference", "INFO")
    
    # Extract SPP λ and amplitude (A) values
    spp_lambdas = {s: sector_results[s]['lambda_km'] for s in sector_results}
    spp_amps   = {s: sector_results[s]['amplitude']   for s in sector_results if 'amplitude' in sector_results[s]}
    
    # ==========================================================================
    # 1. SECTOR-BY-SECTOR COMPARISON TO CODE
    # ==========================================================================
    print_status("\n  SECTOR-BY-SECTOR COMPARISON (SPP vs CODE PPP):", "INFO")
    print_status("  " + "-"*60, "INFO")
    print_status(f"  {'Sector':<8} {'SPP λ':>10} {'CODE λ':>10} {'λ Ratio':>10}   {'SPP A':>10} {'CODE A':>10} {'A Ratio':>8}  Status", "INFO")
    print_status("  " + "-"*75, "INFO")
    
    sector_ratios = {}
    ns_ratios = []  # N-S sector ratios (should be ~1.0 if signal preserved)
    ew_ratios = []  # E-W sector ratios (will be <1.0 due to suppression)
    
    for sector in ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']:
        if sector in spp_lambdas and sector in CODE_SECTOR_LAMBDAS:
            spp_lambda = spp_lambdas[sector]
            code_lambda = CODE_SECTOR_LAMBDAS[sector]
            lam_ratio = spp_lambda / code_lambda
            sector_ratios[sector] = lam_ratio
            
            # Amplitude comparisons (fallback to nan when missing)
            spp_amp  = spp_amps.get(sector, float('nan'))
            code_amp = np.nan  # CODE amplitudes not stored; keep placeholder
            amp_ratio = np.nan
            
            # Classify as preserved (0.7-1.3) or suppressed (<0.7)
            if 0.7 <= lam_ratio <= 1.3:
                status = "✓ PRESERVED"
            elif lam_ratio < 0.7:
                status = "⚠ SUPPRESSED"
            else:
                status = "? ENHANCED"
            
            print_status(
                f"  {sector:<8} {spp_lambda:>10.0f} {code_lambda:>10.0f} {lam_ratio:>10.2f}   "
                f"{spp_amp:>10.2f} {code_amp:>10} {amp_ratio:>8}  {status}",
                "INFO")
            
            if sector in ['N', 'S']:
                ns_ratios.append(lam_ratio)
            elif sector in ['E', 'W']:
                ew_ratios.append(lam_ratio)
    
    print_status("  " + "-"*60, "INFO")
    
    # ==========================================================================
    # 2. GEOMETRIC SUPPRESSION ANALYSIS
    # ==========================================================================
    print_status("\n  GEOMETRIC SUPPRESSION ANALYSIS:", "INFO")
    
    ns_mean_ratio = np.mean(ns_ratios) if ns_ratios else None
    ew_mean_ratio = np.mean(ew_ratios) if ew_ratios else None
    
    if ns_mean_ratio and ew_mean_ratio:
        suppression_factor = ns_mean_ratio / ew_mean_ratio
        
        print_status(f"    N-S sectors mean ratio (SPP/CODE): {ns_mean_ratio:.2f}", "INFO")
        print_status(f"    E-W sectors mean ratio (SPP/CODE): {ew_mean_ratio:.2f}", "INFO")
        print_status(f"    Geometric Suppression Factor (N-S/E-W): {suppression_factor:.2f}", "SUCCESS")
        
        # Interpretation
        if ns_mean_ratio > 0.7 and ew_mean_ratio < 0.5:
            print_status("    → N-S signal PRESERVED, E-W signal SUPPRESSED by orbital geometry", "SUCCESS")
            signal_interpretation = "TEP_SIGNAL_MASKED_BY_GEOMETRY"
        elif ns_mean_ratio > 0.7 and ew_mean_ratio > 0.7:
            print_status("    → Both N-S and E-W signals preserved (unexpected for SPP)", "INFO")
            signal_interpretation = "SIGNAL_PRESERVED"
        else:
            print_status("    → Both directions show suppression (high noise floor)", "WARNING")
            signal_interpretation = "HIGH_NOISE_FLOOR"
    else:
        suppression_factor = None
        signal_interpretation = "INSUFFICIENT_DATA"
    
    # ==========================================================================
    # 3. CORRECTED E-W/N-S RATIO
    # ==========================================================================
    print_status("\n  CORRECTED ANISOTROPY RATIO:", "INFO")
    
    # Raw SPP ratio
    spp_ew = [spp_lambdas[s] for s in ['E', 'W'] if s in spp_lambdas]
    spp_ns = [spp_lambdas[s] for s in ['N', 'S'] if s in spp_lambdas]
    
    if spp_ew and spp_ns:
        raw_ratio = np.mean(spp_ew) / np.mean(spp_ns)
        print_status(f"    Raw SPP E-W/N-S Ratio: {raw_ratio:.2f}", "INFO")
        print_status(f"    CODE Reference Ratio: {CODE_EW_NS_RATIO:.2f}", "INFO")
        
        # Apply geometric correction
        if suppression_factor and suppression_factor > 1:
            corrected_ratio = raw_ratio * suppression_factor
            print_status(f"    Geometry-Corrected Ratio: {corrected_ratio:.2f}", "SUCCESS")
            
            # Compare to CODE
            ratio_match = abs(corrected_ratio - CODE_EW_NS_RATIO) / CODE_EW_NS_RATIO * 100
            print_status(f"    Deviation from CODE: {ratio_match:.1f}%", "INFO")
            
            if ratio_match < 30:
                print_status("    → CORRECTED RATIO MATCHES CODE (within 30%)", "SUCCESS")
                tep_evidence = "STRONG"
            elif ratio_match < 50:
                print_status("    → CORRECTED RATIO PARTIALLY MATCHES CODE", "INFO")
                tep_evidence = "MODERATE"
            else:
                print_status("    → CORRECTED RATIO DOES NOT MATCH CODE", "WARNING")
                tep_evidence = "WEAK"
        else:
            corrected_ratio = None
            tep_evidence = "INCONCLUSIVE"
    else:
        raw_ratio = None
        corrected_ratio = None
        tep_evidence = "INSUFFICIENT_DATA"
    
    # ==========================================================================
    # 4. STATISTICAL SIGNIFICANCE OF N-S PRESERVATION
    # ==========================================================================
    print_status("\n  N-S SIGNAL PRESERVATION TEST:", "INFO")
    
    if ns_ratios and len(ns_ratios) >= 2:
        ns_mean = np.mean(ns_ratios)
        ns_std = np.std(ns_ratios) if len(ns_ratios) > 1 else 0.1
        
        # Test if N-S ratio is significantly different from 0 (null hypothesis: no signal)
        # and consistent with 1.0 (signal preserved)
        t_stat_preservation = (ns_mean - 1.0) / (ns_std / np.sqrt(len(ns_ratios))) if ns_std > 0 else 0
        
        print_status(f"    N-S mean ratio: {ns_mean:.2f} ± {ns_std:.2f}", "INFO")
        print_status(f"    t-statistic (vs 1.0): {t_stat_preservation:.2f}", "INFO")
        
        if abs(t_stat_preservation) < 2.0:  # Not significantly different from 1.0
            print_status("    → N-S signal CONSISTENT with CODE (p > 0.05)", "SUCCESS")
            ns_preserved = True
        else:
            print_status("    → N-S signal DIFFERS from CODE", "WARNING")
            ns_preserved = False
    else:
        ns_preserved = None
    
    # ==========================================================================
    # 5. FINAL TEP EVIDENCE ASSESSMENT
    # ==========================================================================
    print_status("\n  " + "="*60, "INFO")
    print_status("  FINAL TEP EVIDENCE ASSESSMENT", "SUCCESS")
    print_status("  " + "="*60, "INFO")
    
    evidence_score = 0
    evidence_details = []
    
    # Criterion 1: N-S values match CODE (signal preserved)
    if ns_mean_ratio and ns_mean_ratio > 0.7:
        evidence_score += 1
        evidence_details.append("✓ N-S correlation lengths match CODE PPP (signal preserved)")
    else:
        evidence_details.append("✗ N-S correlation lengths do not match CODE")
    
    # Criterion 2: E-W suppression consistent with orbital geometry
    if ew_mean_ratio and ew_mean_ratio < 0.5 and suppression_factor and suppression_factor > 1.5:
        evidence_score += 1
        evidence_details.append("✓ E-W suppression consistent with GPS orbital geometry")
    else:
        evidence_details.append("✗ E-W suppression pattern not as expected")
    
    # Criterion 3: Corrected ratio matches CODE
    if corrected_ratio and abs(corrected_ratio - CODE_EW_NS_RATIO) / CODE_EW_NS_RATIO < 0.3:
        evidence_score += 1
        evidence_details.append("✓ Geometry-corrected ratio matches CODE (within 30%)")
    else:
        evidence_details.append("✗ Corrected ratio does not match CODE")
    
    # Criterion 4: High R² in sector fits (signal, not noise)
    r2_values = [sector_results[s].get('r_squared', 0) for s in sector_results]
    mean_r2 = np.mean(r2_values) if r2_values else 0
    if mean_r2 > 0.8:
        evidence_score += 1
        evidence_details.append(f"✓ High fit quality (mean R² = {mean_r2:.2f})")
    else:
        evidence_details.append(f"✗ Low fit quality (mean R² = {mean_r2:.2f})")
    
    for detail in evidence_details:
        print_status(f"    {detail}", "INFO")
    
    print_status(f"\n    EVIDENCE SCORE: {evidence_score}/4", "SUCCESS" if evidence_score >= 2 else "WARNING")
    
    if evidence_score >= 3:
        conclusion = "STRONG EVIDENCE for TEP signal masked by orbital geometry"
    elif evidence_score >= 2:
        conclusion = "MODERATE EVIDENCE for TEP signal (partial geometric masking)"
    elif evidence_score >= 1:
        conclusion = "WEAK EVIDENCE (signal may be present but heavily masked)"
    else:
        conclusion = "NO EVIDENCE for TEP signal in SPP data"
    
    print_status(f"    CONCLUSION: {conclusion}", "SUCCESS" if evidence_score >= 2 else "WARNING")
    
    return {
        'metric': metric_name,
        'sector_ratios': sector_ratios,
        'ns_mean_ratio': float(ns_mean_ratio) if ns_mean_ratio else None,
        'ew_mean_ratio': float(ew_mean_ratio) if ew_mean_ratio else None,
        'geometric_suppression_factor': float(suppression_factor) if suppression_factor else None,
        'raw_ew_ns_ratio': float(raw_ratio) if raw_ratio else None,
        'corrected_ew_ns_ratio': float(corrected_ratio) if corrected_ratio else None,
        'code_reference_ratio': CODE_EW_NS_RATIO,
        'signal_interpretation': signal_interpretation,
        'tep_evidence_strength': tep_evidence,
        'evidence_score': evidence_score,
        'evidence_max': 4,
        'conclusion': conclusion,
        'ns_preserved': ns_preserved,
        'mean_r_squared': float(mean_r2)
    }


def run_multi_metric_anisotropy(final_pairs, condition_name="ALL"):
    """
    Run anisotropy analysis for MULTIPLE phase metrics (RINEX-specific).
    
    This function extends the CODE longspan methodology by comparing:
    1. MSC (Magnitude Squared Coherence) - normalized cross-spectral power
    2. Phase Alignment Index - cos(weighted_phase), the CODE primary metric
    
    Both metrics should show consistent anisotropy patterns if the signal is real.
    Divergence between metrics may indicate systematic effects.
    
    ALIGNED WITH CODE LONGSPAN:
    - 8 sectors (N, NE, E, SE, S, SW, W, NW) at 45° each
    - 40 log-spaced bins (50-13,000 km)
    - 50 minimum pairs per bin
    - Distance distribution matching for sector fits
    """
    print_status(f"\n{'='*70}", "INFO")
    print_status(f"MULTI-METRIC ANISOTROPY ANALYSIS: {condition_name}", "INFO")
    print_status(f"{'='*70}", "INFO")
    
    if not final_pairs or 'dist' not in final_pairs:
        return None
    
    results = {
        'condition': condition_name,
        'n_pairs': len(final_pairs['dist']),
        'metrics_analyzed': [],
        'metric_results': {},
        'cross_metric_comparison': {}
    }
    
    # Analyze each available metric
    available_metrics = []
    
    # Check which metrics are available
    if 'coherence' in final_pairs and np.sum(np.isfinite(final_pairs['coherence'])) > 1000:
        available_metrics.append('coherence')
    if 'phase_alignment' in final_pairs and np.sum(np.isfinite(final_pairs['phase_alignment'])) > 1000:
        available_metrics.append('phase_alignment')
    
    if not available_metrics:
        print_status("  No valid metrics found in data", "WARNING")
        return None
    
    print_status(f"  Available metrics: {available_metrics}", "INFO")
    results['metrics_analyzed'] = available_metrics
    
    # Run analysis for each metric
    for metric_name in available_metrics:
        print_status(f"\n  --- Analyzing metric: {metric_name} ---", "PROCESS")
        
        metric_info = PHASE_METRICS.get(metric_name, {})
        print_status(f"      Description: {metric_info.get('description', 'Unknown')}", "INFO")
        
        # Create a copy of pairs with the metric as 'coherence' for analyze_anisotropy
        metric_pairs = {
            'dist': final_pairs['dist'],
            'coherence': final_pairs[metric_name],  # Use this metric as the coherence value
            'azimuth': final_pairs['azimuth'],
            'mid_lat': final_pairs['mid_lat']
        }

        if 'lon_diff_deg' in final_pairs:
            metric_pairs['lon_diff_deg'] = final_pairs['lon_diff_deg']
        
        # Run the standard anisotropy analysis
        metric_result = analyze_anisotropy(metric_pairs, f"{condition_name}_{metric_name}")
        
        if metric_result:
            results['metric_results'][metric_name] = metric_result
            
            ew_ns_ratio_sector = metric_result.get('ew_ns_code', {}).get('ratio')
            ew_ns_ratio_pooled = metric_result.get('ew_ns_pooled', {}).get('ratio')
            cv = metric_result.get('anisotropy_cv', 0)
 
            print_status(f"      λ Ratio (Sector Avg): {ew_ns_ratio_sector:.2f}" if ew_ns_ratio_sector else "      λ Ratio (Sector Avg): N/A", "SUCCESS")
            print_status(f"      λ Ratio (Pooled Fit): {ew_ns_ratio_pooled:.2f}" if ew_ns_ratio_pooled else "      λ Ratio (Pooled Fit): N/A", "INFO")
            print_status(f"      Anisotropy CV: {cv:.3f}", "INFO")
    
    # Cross-metric comparison (if both metrics available)
    if len(results['metric_results']) >= 2:
        print_status(f"\n  --- Cross-Metric Comparison ---", "PROCESS")
        
        coh_result = results['metric_results'].get('coherence', {})
        phase_result = results['metric_results'].get('phase_alignment', {})
        
        coh_ratio = coh_result.get('ew_ns_code', {}).get('ratio')
        phase_ratio = phase_result.get('ew_ns_code', {}).get('ratio')
        
        coh_cv = coh_result.get('anisotropy_cv', 0)
        phase_cv = phase_result.get('anisotropy_cv', 0)
        
        comparison = {
            'coherence_ew_ns_ratio': coh_ratio,
            'phase_alignment_ew_ns_ratio': phase_ratio,
            'ratio_difference': abs(coh_ratio - phase_ratio) if (coh_ratio and phase_ratio) else None,
            'coherence_cv': coh_cv,
            'phase_alignment_cv': phase_cv,
            'cv_difference': abs(coh_cv - phase_cv),
            'metrics_consistent': False
        }
        
        # Check consistency (ratios within 20% of each other)
        if coh_ratio and phase_ratio:
            ratio_diff_pct = abs(coh_ratio - phase_ratio) / max(coh_ratio, phase_ratio) * 100
            comparison['ratio_difference_percent'] = ratio_diff_pct
            comparison['metrics_consistent'] = ratio_diff_pct < 30  # 30% tolerance
            
            print_status(f"      MSC λ Ratio (Sector Avg): {coh_ratio:.2f}", "INFO")
            print_status(f"      Phase λ Ratio (Sector Avg): {phase_ratio:.2f}", "INFO")
            print_status(f"      Difference: {ratio_diff_pct:.1f}%", "INFO")
            
            if comparison['metrics_consistent']:
                print_status(f"      CONSISTENT: Both metrics show similar anisotropy", "SUCCESS")
            else:
                print_status(f"      DIVERGENT: Metrics show different anisotropy patterns", "WARNING")
        
        results['cross_metric_comparison'] = comparison
    
    # ==========================================================================
    # ENHANCED TEP EVIDENCE ANALYSIS - Compare to CODE reference
    # ==========================================================================
    # Run TEP evidence analysis on the phase_alignment metric (CODE methodology)
    # This compares SPP sector λ values to CODE PPP reference values
    print_status("\n" + "="*70, "INFO")
    print_status("RUNNING ENHANCED TEP EVIDENCE ANALYSIS", "PROCESS")
    print_status("="*70, "INFO")
    
    tep_evidence_results = {}
    
    for metric_name in results['metric_results']:
        metric_result = results['metric_results'][metric_name]
        # Use the matched sector results for comparison (more robust)
        sector_results = metric_result.get('sectors_matched', metric_result.get('sectors', {}))
        
        if sector_results and len(sector_results) >= 4:
            tep_evidence = analyze_tep_evidence(sector_results, metric_name)
            if tep_evidence:
                tep_evidence_results[metric_name] = tep_evidence
    
    results['tep_evidence_analysis'] = tep_evidence_results
    
    # Summary across metrics
    if tep_evidence_results:
        print_status("\n" + "="*70, "INFO")
        print_status("TEP EVIDENCE SUMMARY ACROSS METRICS", "SUCCESS")
        print_status("="*70, "INFO")
        
        for metric, evidence in tep_evidence_results.items():
            score = evidence.get('evidence_score', 0)
            conclusion = evidence.get('conclusion', 'N/A')
            print_status(f"  {metric}: Score {score}/4 - {conclusion}", 
                        "SUCCESS" if score >= 2 else "WARNING")
        
        # Overall assessment
        max_score = max(e.get('evidence_score', 0) for e in tep_evidence_results.values())
        results['tep_evidence_max_score'] = max_score
        
        if max_score >= 3:
            results['tep_overall_assessment'] = "STRONG_EVIDENCE"
            print_status("\n  OVERALL: STRONG EVIDENCE for TEP signal in SPP data", "SUCCESS")
        elif max_score >= 2:
            results['tep_overall_assessment'] = "MODERATE_EVIDENCE"
            print_status("\n  OVERALL: MODERATE EVIDENCE for TEP signal", "SUCCESS")
        elif max_score >= 1:
            results['tep_overall_assessment'] = "WEAK_EVIDENCE"
            print_status("\n  OVERALL: WEAK EVIDENCE (signal may be present)", "WARNING")
        else:
            results['tep_overall_assessment'] = "NO_EVIDENCE"
            print_status("\n  OVERALL: NO EVIDENCE for TEP signal detected", "WARNING")
    
    return results


def compute_pairs_coherence(day_files, station_coords, condition_name="ALL", data_key='clock_bias_ns'):
    """
    Compute coherence using day-by-day batch processing (CODE Methodology).
    
    Replaces old pair-based computation with day-based batches.
    """
    print_status(f"[{condition_name}] Computing Day-by-Day Coherence (Parallel)...", "PROCESS")
    
    # Get all day keys
    day_keys = sorted(list(day_files.keys()))
    n_days = len(day_keys)
    print_status(f"  Days to process: {n_days}", "INFO")
    
    # Create batches of days
    n_workers = max(1, multiprocessing.cpu_count() // 2)
    batch_size = max(1, n_days // (n_workers * 4))
    batches = [day_keys[i:i + batch_size] for i in range(0, n_days, batch_size)]
    
    print_status(f"  Workers: {n_workers}, Batches: {len(batches)} (size ~{batch_size})", "INFO")
    
    all_results = []
    
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = []
        for batch in batches:
            # Pass minimal data to workers
            args = (batch, day_files, data_key, station_coords)
            futures.append(executor.submit(_process_day_batch_worker, args))
            
        # Process results
        completed = 0
        for future in as_completed(futures):
            try:
                res = future.result()
                if res:
                    all_results.extend(res)
                completed += 1
                if completed % 10 == 0 or completed == len(batches):
                    print(f"\r  Batches: {completed}/{len(batches)} | Results: {len(all_results):,}", end="", flush=True)
            except Exception as e:
                print(f"\nWorker failed: {e}")
                
    print(f"\r  Completed: {len(all_results):,} total pair-day measurements.      ")
    
    # Convert to DataFrame-like list of dicts for downstream analysis
    if not all_results:
        print_status("No valid coherence results found", "WARNING")
        return []
        
    return all_results


def load_pairs_from_csv(csv_file, metric_key, station_filter_list=None, station_coords=None):
    """
    Load pairs from Step 2.0 CSV output efficiently using chunks.
    Returns a DICT OF ARRAYS (SOA) for efficient processing.
    
    If station_coords is provided, recomputes azimuths from coordinates
    to fix the direction bug in step_2_0 (azimuth was stored in wrong direction
    for ~50% of pairs where idx1 > idx2).
    """
    try:
        if not csv_file.exists(): return None
        print_status(f"Loading pre-computed pairs from {csv_file.name} (Chunked)...", "PROCESS")
        
        # Check headers first
        header = pd.read_csv(csv_file, nrows=0).columns.tolist()
        
        # Determine columns to load
        use_cols = ['distance_km', 'coherence', 'azimuth', 'mid_lat', 'station1', 'station2', 'year', 'doy']
        
        # Add 'metric' column if it exists
        has_metric_col = 'metric' in header
        if has_metric_col:
            use_cols.append('metric')

        metric_filter = None
        if metric_key is None:
            metric_filter = None
        elif isinstance(metric_key, (list, tuple, set)):
            metric_filter = set([str(m) for m in metric_key])
        elif isinstance(metric_key, str) and metric_key.lower() in ['all', '*']:
            metric_filter = None
        else:
            metric_filter = set([str(metric_key)])

        if has_metric_col:
            if metric_filter is None:
                print_status("  Will load ALL metrics from CSV (no metric filtering)", "INFO")
            elif len(metric_filter) == 1:
                print_status(f"  Will filter for '{list(metric_filter)[0]}' metric only", "INFO")
            else:
                shown = ",".join(sorted(list(metric_filter))[:6])
                suffix = "" if len(metric_filter) <= 6 else "..."
                print_status(f"  Will filter for {len(metric_filter)} metrics: {shown}{suffix}", "INFO")
        
        # Check for phase alignment column (could be 'phase_alignment' or 'weighted_phase')
        phase_col = None
        if 'phase_alignment' in header:
            phase_col = 'phase_alignment'
            use_cols.append('phase_alignment')
        elif 'weighted_phase' in header:
            phase_col = 'weighted_phase'
            use_cols.append('weighted_phase')
            
        # Initialize storage
        data_arrays = defaultdict(list)
        total_rows_raw = 0
        total_rows = 0
        metric_counts_raw = defaultdict(int)
        metric_counts_kept = defaultdict(int)
        year_counts_raw = defaultdict(int)
        year_counts_kept = defaultdict(int)
        unique_days_raw = set()
        unique_days_kept = set()
        chunk_size = 500000  # 500k rows per chunk
        
        # Create filter set if needed
        sta_filter = set(station_filter_list) if station_filter_list else None

        lon_lookup = None
        if station_coords is not None:
            lon_lookup = {sta: float(v.get('lon')) for sta, v in station_coords.items() if v and 'lon' in v}
        
        # Read in chunks
        for i, chunk in enumerate(pd.read_csv(csv_file, usecols=use_cols, chunksize=chunk_size)):
            if has_metric_col:
                vc = chunk['metric'].value_counts(dropna=False)
                for k, v in vc.items():
                    metric_counts_raw[str(k)] += int(v)

            y_vc = chunk['year'].value_counts(dropna=False)
            for y, v in y_vc.items():
                year_counts_raw[int(y)] += int(v)

            try:
                day_ids = (chunk['year'].astype(np.int64) * 1000 + chunk['doy'].astype(np.int64)).unique()
                unique_days_raw.update([int(x) for x in day_ids.tolist()])
            except Exception:
                pass

            total_rows_raw += int(len(chunk))

            if has_metric_col and metric_filter is not None:
                chunk = chunk[chunk['metric'].astype(str).isin(metric_filter)]
            
            # Filter by stations if needed
            if sta_filter:
                mask = chunk['station1'].isin(sta_filter) & chunk['station2'].isin(sta_filter)
                chunk = chunk[mask]
            
            if chunk.empty: continue

            if has_metric_col:
                vc_k = chunk['metric'].value_counts(dropna=False)
                for k, v in vc_k.items():
                    metric_counts_kept[str(k)] += int(v)

            y_vc_k = chunk['year'].value_counts(dropna=False)
            for y, v in y_vc_k.items():
                year_counts_kept[int(y)] += int(v)

            try:
                day_ids_k = (chunk['year'].astype(np.int64) * 1000 + chunk['doy'].astype(np.int64)).unique()
                unique_days_kept.update([int(x) for x in day_ids_k.tolist()])
            except Exception:
                pass
            
            # Append arrays
            data_arrays['dist'].append(chunk['distance_km'].values.astype(np.float32))
            data_arrays['coherence'].append(chunk['coherence'].values.astype(np.float32))
            data_arrays['azimuth'].append(chunk['azimuth'].values.astype(np.float32))
            data_arrays['station1'].append(chunk['station1'].values)
            data_arrays['station2'].append(chunk['station2'].values)
            data_arrays['mid_lat'].append(chunk['mid_lat'].values.astype(np.float32))
            # Use int32 to prevent any overflow risk during ID generation (YYYY*1000 + DOY)
            data_arrays['year'].append(chunk['year'].values.astype(np.int32))
            data_arrays['doy'].append(chunk['doy'].values.astype(np.int32))

            if lon_lookup is not None:
                lon1 = chunk['station1'].map(lon_lookup).astype(np.float32).values
                lon2 = chunk['station2'].map(lon_lookup).astype(np.float32).values
                lon_diff = np.abs(lon1 - lon2)
                lon_diff = np.where(lon_diff > 180.0, 360.0 - lon_diff, lon_diff)
                data_arrays['lon_diff_deg'].append(lon_diff.astype(np.float32))
            
            if phase_col:
                data_arrays['phase_alignment'].append(chunk[phase_col].values.astype(np.float32))
            
            total_rows += len(chunk)
            if (i + 1) % 5 == 0:
                print_status(f"  Loaded {total_rows/1e6:.1f}M pairs...", "INFO")
                gc.collect()
        
        if total_rows == 0:
            return None

        if total_rows_raw > 0:
            years_raw_sorted = sorted([int(y) for y in year_counts_raw.keys()])
            years_kept_sorted = sorted([int(y) for y in year_counts_kept.keys()])

            if has_metric_col:
                raw_items = sorted(metric_counts_raw.items(), key=lambda kv: kv[1], reverse=True)
                kept_items = sorted(metric_counts_kept.items(), key=lambda kv: kv[1], reverse=True)
                raw_str = ", ".join([f"{k}:{v:,}" for k, v in raw_items[:6]])
                kept_str = ", ".join([f"{k}:{v:,}" for k, v in kept_items[:6]])
                print_status(f"  CSV audit: raw_rows={total_rows_raw:,} kept_rows={total_rows:,}", "INFO")
                print_status(f"  Metrics(raw): {raw_str}", "INFO")
                print_status(f"  Metrics(kept): {kept_str}", "INFO")

            print_status(f"  Years(raw): {years_raw_sorted} | unique_days(raw)={len(unique_days_raw):,}", "INFO")
            print_status(f"  Years(kept): {years_kept_sorted} | unique_days(kept)={len(unique_days_kept):,}", "INFO")
            
        # Concatenate arrays
        print_status(f"  Finalizing data structure ({total_rows} pairs)...", "INFO")
        final_data = {
            'dist': np.concatenate(data_arrays['dist']),
            'coherence': np.concatenate(data_arrays['coherence']),
            'azimuth': np.concatenate(data_arrays['azimuth']),
            'mid_lat': np.concatenate(data_arrays['mid_lat']),
            'year': np.concatenate(data_arrays['year']),
            'doy': np.concatenate(data_arrays['doy'])
        }

        if 'lon_diff_deg' in data_arrays and data_arrays['lon_diff_deg']:
            final_data['lon_diff_deg'] = np.concatenate(data_arrays['lon_diff_deg'])
        
        if phase_col:
            # If we loaded weighted_phase (radians), convert to phase_alignment (cosine)
            # This is critical because Weighted Phase is an angle (-pi to pi), not a correlation!
            raw_phase_data = np.concatenate(data_arrays['phase_alignment'])
            
            if phase_col == 'weighted_phase':
                print_status("  Converting Weighted Phase (radians) to Phase Alignment (cosine)...", "INFO")
                final_data['phase_alignment'] = np.cos(raw_phase_data)
            else:
                final_data['phase_alignment'] = raw_phase_data
        
        # ======================================================================
        # CRITICAL FIX: Recompute azimuths from station coordinates
        # The CSV has a bug where ~50% of azimuths are stored in the wrong
        # direction (from station with min_idx to max_idx, not sta1 to sta2).
        # Recomputing from coordinates ensures correct sta1→sta2 direction.
        # ======================================================================
        if station_coords is not None:
            print_status("  Recomputing azimuths from station coordinates (fixing direction bug)...", "INFO")
            sta1_arr = np.concatenate(data_arrays['station1'])
            sta2_arr = np.concatenate(data_arrays['station2'])
            
            # Build unique pair azimuths lookup (much faster than per-row computation)
            # Use string keys for efficient pandas map
            unique_pairs = set(zip(sta1_arr, sta2_arr))
            pair_az_lookup = {}
            for sta1, sta2 in unique_pairs:
                if sta1 in station_coords and sta2 in station_coords:
                    lat1, lon1 = station_coords[sta1]['lat'], station_coords[sta1]['lon']
                    lat2, lon2 = station_coords[sta2]['lat'], station_coords[sta2]['lon']
                    key = f"{sta1}|{sta2}"
                    pair_az_lookup[key] = compute_azimuth(lat1, lon1, lat2, lon2)
            
            print_status(f"  Built lookup for {len(pair_az_lookup):,} unique pairs", "INFO")
            
            # Create pair key strings efficiently
            pair_keys = np.char.add(np.char.add(sta1_arr.astype(str), '|'), sta2_arr.astype(str))
            
            # Apply lookup using pandas Series map
            old_azimuths = final_data['azimuth'].copy()
            pair_series = pd.Series(pair_keys)
            mapped = pair_series.map(pair_az_lookup)
            new_azimuths = np.where(mapped.isna(), old_azimuths, mapped.values).astype(np.float32)
            
            # Report how many were changed
            diff = np.abs(new_azimuths - old_azimuths)
            # Handle wraparound (e.g., 350° vs 10° = 20° diff, not 340°)
            diff = np.minimum(diff, 360 - diff)
            # Azimuths that differ by ~180° are the buggy ones
            buggy_mask = (diff > 170) & (diff < 190)
            n_fixed = np.sum(buggy_mask)
            print_status(f"  Fixed {n_fixed:,} azimuths ({100*n_fixed/total_rows:.1f}% were ~180° off)", "SUCCESS")
            
            final_data['azimuth'] = new_azimuths
            del sta1_arr, sta2_arr, old_azimuths
            
        del data_arrays
        gc.collect()
        
        print_status(f"  Loaded {total_rows} pairs successfully.", "SUCCESS")
        return final_data
        
    except Exception as e:
        print_status(f"  Failed to load CSV: {e}", "WARNING")
        import traceback
        traceback.print_exc()
        return None

def run_analysis_pipeline(filter_mode_name, station_filter_config, observables=None):
    """Run the full analysis pipeline for a specific filter mode."""
    print_status("", "INFO")
    print_status("=" * 80, "INFO")
    print_status(f"TEP-GNSS-RINEX Analysis - STEP 2.2: Anisotropy ({filter_mode_name})", "INFO")
    print_status("=" * 80, "INFO")
        
    
    # Define output filename based on filter
    filter_suffix = filter_mode_name.lower().replace(" ", "_").replace("-", "_")
    
    # Map filter config to CSV suffix for pre-generated pair files
    # CSV files are named: step_2_0_pairs_{mode}_{filter_suffix}.csv
    if station_filter_config == 'none':
        csv_filter_suffix = 'all_stations'
    elif station_filter_config == 'optimal_100_metadata.json':
        csv_filter_suffix = 'optimal_100'
    elif station_filter_config.startswith('dynamic:'):
        threshold = station_filter_config.split(':')[1]
        csv_filter_suffix = f'dynamic_{threshold}'
    else:
        csv_filter_suffix = 'all_stations'  # fallback
    output_file = OUTPUTS_DIR / f"step_2_2_anisotropy_analysis_{filter_suffix}.json"
    
    # ==========================================================================
    # CHECK FOR PRE-GENERATED CSV FILES FIRST (FAST PATH)
    # ==========================================================================
    # If CSV files exist, we can skip the expensive NPZ scanning entirely
    csv_baseline = OUTPUTS_DIR / f"step_2_0_pairs_baseline_{csv_filter_suffix}.csv"
    csv_ionofree = OUTPUTS_DIR / f"step_2_0_pairs_ionofree_{csv_filter_suffix}.csv"
    csv_multi_gnss = OUTPUTS_DIR / f"step_2_0_pairs_multi_gnss_{csv_filter_suffix}.csv"
    csv_precise = OUTPUTS_DIR / f"step_2_0_pairs_precise_{csv_filter_suffix}.csv"
    
    csv_mode = csv_baseline.exists()  # Use CSV mode if baseline CSV exists
    if csv_mode:
        print_status(f"[FAST PATH] Found pre-generated CSV files for {csv_filter_suffix}", "SUCCESS")
        print_status(f"  Baseline: {csv_baseline.name} ({csv_baseline.stat().st_size / 1e6:.1f} MB)", "INFO")
        if csv_ionofree.exists():
            print_status(f"  Ionofree: {csv_ionofree.name}", "INFO")
        if csv_multi_gnss.exists():
            print_status(f"  Multi-GNSS: {csv_multi_gnss.name}", "INFO")
        if csv_precise.exists():
            print_status(f"  Precise: {csv_precise.name}", "INFO")
    else:
        print_status(f"  No pre-generated CSV found at {csv_baseline.name}, will process NPZ files", "WARNING")
    
    # ==========================================================================
    # STEP 1: Load station data
    # ==========================================================================
    print_status("[1/8] Loading station data...", "PROCESS")
    
    dynamic_threshold = None
    good_stations_set = None
    
    # Configure filter
    if station_filter_config == 'none':
        print_status("  Using ALL stations (no filter)", "INFO")
    elif station_filter_config.startswith('dynamic:'):
        dynamic_threshold = float(station_filter_config.split(':')[1])
        print_status(f"  Using DYNAMIC filter: std < {dynamic_threshold}ns per file", "INFO")
    else:
        good_stations_file = PROCESSED_DIR / station_filter_config
        if good_stations_file.exists():
            with open(good_stations_file) as f:
                good_config = json.load(f)
            good_stations_set = set(good_config['stations'])
            filter_desc = good_config.get('description', station_filter_config)
            print_status(f"  Using filter: {station_filter_config}", "INFO")
            print_status(f"  Description: {filter_desc}", "INFO")
            print_status(f"  Stations: {len(good_stations_set)}", "INFO")
        else:
            print_status(f"  Filter file not found: {station_filter_config}, using all stations", "WARNING")

    # Load coordinates
    coords_file = PROCESSED_DIR / "station_coordinates.json"
    if not coords_file.exists():
        print_status("Station coordinates not found. Run Step 2.0 first.", "ERROR")
        return
    
    with open(coords_file) as f:
        station_coords_ecef = json.load(f)
    
    station_coords = {}
    for sta, ecef in station_coords_ecef.items():
        if good_stations_set and sta not in good_stations_set:
            continue
        lat, lon = ecef_to_lla(ecef[0], ecef[1], ecef[2])
        station_coords[sta] = {'lat': lat, 'lon': lon}
    
    # ==========================================================================
    # NPZ SCANNING - ONLY IF CSV MODE IS DISABLED (SLOW PATH)
    # ==========================================================================
    day_files = {}
    station_files = {}
    station_modes = {}
    all_years = set()
    
    if csv_mode:
        # FAST PATH: Skip NPZ scanning entirely, determine mode availability from CSV existence
        print_status("  [FAST PATH] Skipping NPZ file scanning (using CSV files)", "SUCCESS")
 
        all_years = set()
        
        # Mode availability based on CSV file existence
        modes_available = {
            'baseline': csv_baseline.exists(),
            'precise': csv_precise.exists(),
            'ionofree': csv_ionofree.exists(),
            'multi_gnss': csv_multi_gnss.exists(),
        }
        
        print_status(f"  Mode availability (from CSV files):", "INFO")
        for mode, available in modes_available.items():
            status = "✓" if available else "✗"
            print_status(f"    {mode}: {status}", "INFO" if available else "WARNING")
    else:
        # SLOW PATH: Scan NPZ files (fallback when no CSV exists)
        print_status("  [SLOW PATH] Scanning NPZ files...", "WARNING")
        
        npz_files = sorted(PROCESSED_DIR.glob("*.npz"))
        
        # Apply dynamic filter if specified (filter by clock stats - IDENTICAL to step 2.0)
        if dynamic_threshold is not None:
            good_files = []
            rejected_jump = 0
            rejected_range = 0
            for f in npz_files:
                try:
                    with np.load(f) as d:
                        clk = d.get('clock_bias_ns')
                        if clk is not None and len(clk) > 10:
                            # Smart criteria (identical to step 2.0):
                            jumps = np.abs(np.diff(clk))
                            max_jump = np.nanmax(jumps) if len(jumps) > 0 else 0
                            total_range = np.nanmax(clk) - np.nanmin(clk)
                            std = np.nanstd(clk)
                            
                            if max_jump > 500:
                                rejected_jump += 1
                            elif total_range > 5000:
                                rejected_range += 1
                            elif std < dynamic_threshold:
                                good_files.append(f)
                except:
                    continue
            total_rejected = len(npz_files) - len(good_files)
            print_status(f"  Smart filter: {len(good_files)}/{len(npz_files)} files passed", "INFO")
            print_status(f"    Rejected: {rejected_jump} (jumps>500ns), {rejected_range} (range>5000ns), {total_rejected - rejected_jump - rejected_range} (std>{dynamic_threshold}ns)", "INFO")
            npz_files = good_files
        elif good_stations_set:
            npz_files = [f for f in npz_files if f.name.split('_')[0] in good_stations_set]
        
        print_status(f"  Found {len(npz_files)} processed files matching filter", "INFO")
        
        # Group files by DAY (year, doy) - CODE LONGSPAN ARCHITECTURE
        for f in npz_files:
            try:
                parts = f.name.split('_')
                sta = parts[0]
                if sta not in station_coords:
                    continue
                
                # Parse YYYYDDD
                year_doy = parts[1].replace('.npz', '')
                year = int(year_doy[:4])
                doy = int(year_doy[4:7])
                
                day_key = (year, doy)
                if day_key not in day_files:
                    day_files[day_key] = {}
                day_files[day_key][sta] = f
                
                # Also maintain station_files for metadata
                if sta not in station_files:
                    station_files[sta] = []
                station_files[sta].append(f)
            except:
                continue
        
        print_status(f"  Grouped into {len(day_files)} days across {len(station_files)} stations", "INFO")
        
        # Detect years from all files
        for f in npz_files:
            try:
                year_str = f.name.split('_')[1][:4]
                all_years.add(int(year_str))
            except:
                continue
                
        if not all_years:
            all_years = {2022}
            print_status("  Could not detect years, defaulting to 2022", "WARNING")
        else:
            print_status(f"  Detected years: {sorted(list(all_years))}", "INFO")
        
        # Check station modes from NPZ files
        print_status(f"  Preparing for day-by-day processing...", "PROCESS")
        print_status(f"  Total days to process: {len(day_files)}", "INFO")
        print_status(f"  Stations available: {len(station_files)}", "INFO")
        
        print(f"\r  Checking station modes...", end="", flush=True)
        for sta in list(station_files.keys()):
            sample_file = station_files[sta][0]
            try:
                with np.load(sample_file) as d:
                    station_modes[sta] = {
                        'baseline': 'clock_bias_ns' in d,
                        'ionofree': 'ionofree_clock_bias_ns' in d,
                        'mgex': 'mgex_clock_bias_ns' in d,
                        'precise': 'precise_clock_bias_ns' in d,
                        'lat': station_coords[sta]['lat'],
                        'lon': station_coords[sta]['lon']
                    }
            except:
                continue
        
        # Count stations per mode
        n_baseline = sum(1 for s in station_modes.values() if s.get('baseline'))
        n_ionofree = sum(1 for s in station_modes.values() if s.get('ionofree'))
        n_mgex = sum(1 for s in station_modes.values() if s.get('mgex'))
        n_precise = sum(1 for s in station_modes.values() if s.get('precise'))
        
        print(f"\r  Station modes detected:              ")
        print_status(f"  Baseline mode stations: {n_baseline}", "SUCCESS")
        print_status(f"  Precise mode stations: {n_precise} (Option E)", "INFO")
        print_status(f"  Iono-free mode stations: {n_ionofree} (Option D)", "INFO")
        print_status(f"  Multi-GNSS mode stations: {n_mgex} (MGEX)", "INFO")
        
        # Store mode availability (from NPZ scan)
        modes_available = {
            'baseline': n_baseline >= 20,
            'precise': n_precise >= 20,
            'ionofree': n_ionofree >= 20,
            'multi_gnss': n_mgex >= 20,
        }
    
    kp_day_status = {}
    quiet_ids = set()
    storm_ids = set()
    days_all = []
    days_quiet = []
    days_storm = []

    def get_day_files_subset(days):
        return {d: day_files[d] for d in days if d in day_files}

    # ==========================================================================
    # STEP 2: Analyze ALL DAYS
    # ==========================================================================

    if observables is None:
        observables = ['clock_bias']
    if isinstance(observables, str):
        observables = [o.strip() for o in observables.split(',') if o.strip()]
    if not observables:
        observables = ['clock_bias']

    if 'all' in [o.lower() for o in observables]:
        observables = ['clock_bias', 'pos_jitter', 'clock_drift']

    primary_observable = observables[0]
    control_observables = [o for o in observables[1:]]
    
    # Try loading pairs from Step 2.0 output (Baseline)
    # Use filter-specific CSV file (e.g., step_2_0_pairs_baseline_optimal_100.csv)
    pair_file_baseline = OUTPUTS_DIR / f"step_2_0_pairs_baseline_{csv_filter_suffix}.csv"
    pairs_all_loaded = load_pairs_from_csv(pair_file_baseline, primary_observable, list(station_coords.keys()), station_coords)

    step_2_0_summary = None
    if csv_mode:
        step_2_0_json = OUTPUTS_DIR / f"step_2_0_raw_spp_analysis_{csv_filter_suffix}.json"
        if step_2_0_json.exists():
            try:
                with open(step_2_0_json, 'r') as f:
                    step_2_0_summary = json.load(f)
                print_status(f"  Step 2.0 summary found: {step_2_0_json.name}", "SUCCESS")
            except Exception as e:
                print_status(f"  Failed to read Step 2.0 summary JSON: {e}", "WARNING")

    if step_2_0_summary and isinstance(step_2_0_summary, dict):
        try:
            abm = step_2_0_summary.get('analysis_by_mode', {})
            if isinstance(abm, dict) and abm:
                print_status("\n  STEP 2.0 PAIR COUNTS (audit; within 50–13,000 km):", "INFO")
                for mode_name in ['baseline', 'ionofree', 'multi_gnss', 'precise']:
                    mode_block = abm.get(mode_name, {})
                    if not isinstance(mode_block, dict) or not mode_block:
                        continue
                    metric_rows = []
                    for metric_name, res in mode_block.items():
                        if not isinstance(res, dict):
                            continue
                        n_pairs = res.get('n_pairs', None)
                        if n_pairs is None:
                            continue
                        metric_rows.append((metric_name, int(n_pairs)))
                    metric_rows.sort(key=lambda x: x[1], reverse=True)
                    if metric_rows:
                        top = ", ".join([f"{m}:{n:,}" for m, n in metric_rows[:6]])
                        print_status(f"    {mode_name}: {top}", "INFO")

                if pairs_all_loaded and isinstance(pairs_all_loaded, dict) and 'dist' in pairs_all_loaded:
                    loaded_n = int(len(pairs_all_loaded['dist']))
                    expected_n = None
                    base_block = abm.get('baseline', {}) if isinstance(abm, dict) else {}
                    if isinstance(base_block, dict):
                        res0 = base_block.get(primary_observable, {})
                        if isinstance(res0, dict) and 'n_pairs' in res0:
                            expected_n = int(res0.get('n_pairs'))

                    if expected_n is not None and expected_n > 0:
                        frac = abs(loaded_n - expected_n) / float(expected_n)
                        msg = f"  Step 2.0 vs Step 2.2 count check (baseline/{primary_observable}): Step2.2={loaded_n:,} Step2.0={expected_n:,} Δ={100*frac:.2f}%"
                        if frac > 0.01:
                            print_status(msg, "WARNING")
                        else:
                            print_status(msg, "SUCCESS")
        except Exception as e:
            print_status(f"  Step 2.0 audit summary parse failed: {e}", "WARNING")
    
    if pairs_all_loaded:
        print_status("Using pre-computed pairs from Step 2.0 output", "SUCCESS")
        pairs_all = pairs_all_loaded
    else:
        print_status("[4/8] Analyzing ALL DAYS (Day-by-Day)...", "PROCESS")
        key_map = {
            'clock_bias': 'clock_bias_ns',
            'pos_jitter': 'pos_jitter',
            'clock_drift': 'clock_drift'
        }
        data_key = key_map.get(primary_observable, 'clock_bias_ns')
        pairs_all = compute_pairs_coherence(day_files, station_coords, "ALL", data_key)
    
    years_present = None
    if isinstance(pairs_all, dict) and 'year' in pairs_all:
        try:
            years_present = sorted(list(set(np.unique(pairs_all['year']).astype(int).tolist())))
        except Exception:
            years_present = None

    if years_present:
        print_status(f"  Data years detected ({primary_observable}): {years_present}", "INFO")

    # Run MULTI-METRIC analysis (compares MSC and phase_alignment)
    result_all_multi = run_multi_metric_anisotropy(pairs_all, f"ALL DAYS (Multi-Year) [{primary_observable}]")
    
    # Also run single-metric for backward compatibility
    result_all = analyze_anisotropy(pairs_all, f"ALL DAYS (Multi-Year) [{primary_observable}]")

    control_results = {}
    control_results_multi = {}
    if control_observables and pairs_all_loaded:
        print_status("\n[Control] Running additional observables (same pipeline; explicit null/robustness tests)...", "PROCESS")
        for obs in control_observables:
            try:
                pairs_ctrl = load_pairs_from_csv(pair_file_baseline, obs, list(station_coords.keys()), station_coords)
                if not pairs_ctrl:
                    continue

                ctrl_multi = run_multi_metric_anisotropy(pairs_ctrl, f"ALL DAYS (Multi-Year) [{obs}]")
                ctrl_single = analyze_anisotropy(pairs_ctrl, f"ALL DAYS (Multi-Year) [{obs}]")
                control_results[obs] = ctrl_single
                control_results_multi[obs] = ctrl_multi

                sd = ctrl_single.get('short_distance_analysis', {}) if ctrl_single else {}
                sd_ratio = sd.get('mean_ratio')
                sd_str = f"{sd_ratio:.3f}" if sd_ratio is not None else "N/A"
                print_status(f"  [Control] {obs}: N={ctrl_single.get('n_pairs', 0):,} short-dist ratio={sd_str}", "INFO")
            except Exception as e:
                print_status(f"  [Control] {obs} failed: {e}", "WARNING")

    # ==========================================================================
    # Fetch Kp data and create quiet/storm day lookup
    # ==========================================================================
    years_for_kp = set()
    if day_files:
        years_for_kp = {y for (y, _d) in day_files.keys()}
    elif pairs_all_loaded and isinstance(pairs_all_loaded, dict) and 'year' in pairs_all_loaded:
        years_for_kp = set(np.unique(pairs_all_loaded['year']).astype(int).tolist())
    elif all_years:
        years_for_kp = set(all_years)
    else:
        years_for_kp = {2022}

    print_status("[2/8] Fetching Kp data for stratification...", "PROCESS")
    kp_data = fetch_kp_data_gfz(sorted(list(years_for_kp)))

    for date_str, kp in kp_data.items():
        try:
            year = int(date_str[:4])
            month = int(date_str[4:6])
            day = int(date_str[6:8])
            date = datetime(year, month, day)
            doy = date.timetuple().tm_yday
            day_id = year * 1000 + doy
            if kp < KP_QUIET_THRESHOLD:
                kp_day_status[day_id] = 'quiet'
                quiet_ids.add(day_id)
            else:
                kp_day_status[day_id] = 'storm'
                storm_ids.add(day_id)
        except:
            pass

    print_status(f"  Kp lookup built: {len(quiet_ids)} quiet days, {len(storm_ids)} storm days", "INFO")

    if day_files:
        days_all = sorted(list(day_files.keys()))
        days_quiet = [(y, d) for (y, d) in days_all if y * 1000 + d in quiet_ids]
        days_storm = [(y, d) for (y, d) in days_all if y * 1000 + d in storm_ids]
        print_status(f"  NPZ days: {len(days_quiet)} quiet, {len(days_storm)} storm (of {len(days_all)} total)", "INFO")
    
    # ==========================================================================
    # STEP 6: Analyze QUIET DAYS
    # ==========================================================================
    print_status("[6/8] Analyzing QUIET DAYS...", "PROCESS")
    
    if pairs_all_loaded:
        # Filter from loaded pairs using Vectorized Masking
        # Use pre-built quiet_ids set from Kp data (works in both CSV and NPZ modes)
        print_status("  Creating mask for quiet days...", "INFO")
        
        # Compute IDs for all pairs (YYYYDDD format)
        pair_ids = pairs_all_loaded['year'] * 1000 + pairs_all_loaded['doy']
        
        # Create mask using pre-built quiet_ids set
        quiet_mask = np.isin(pair_ids, list(quiet_ids))
        
        # Apply mask to all arrays
        pairs_quiet = {}
        for key, arr in pairs_all_loaded.items():
            pairs_quiet[key] = arr[quiet_mask]
            
        print_status(f"  Filtered {len(pairs_quiet['dist']):,} pairs (Quiet Days)", "INFO")
    else:
        # Re-compute from NPZ files
        day_files_quiet = get_day_files_subset(days_quiet)
        key_map = {
            'clock_bias': 'clock_bias_ns',
            'pos_jitter': 'pos_jitter',
            'clock_drift': 'clock_drift'
        }
        data_key = key_map.get(primary_observable, 'clock_bias_ns')
        pairs_quiet = compute_pairs_coherence(day_files_quiet, station_coords, "QUIET", data_key)
    
    # Run MULTI-METRIC analysis on Quiet Days (key for TEP evidence)
    result_quiet_multi = run_multi_metric_anisotropy(pairs_quiet, f"QUIET DAYS (Kp<{KP_QUIET_THRESHOLD}) [{primary_observable}]")
    result_quiet = analyze_anisotropy(pairs_quiet, f"QUIET DAYS (Kp<{KP_QUIET_THRESHOLD}) [{primary_observable}]")
    
    # ==========================================================================
    # STEP 7: Analyze STORM DAYS
    # ==========================================================================
    print_status("[7/8] Analyzing STORM DAYS...", "PROCESS")
    
    if pairs_all_loaded:
        # Filter from loaded pairs using pre-built storm_ids set
        # Reuse pair_ids from quiet day filtering
        storm_mask = np.isin(pair_ids, list(storm_ids))
        
        pairs_storm = {}
        for key, arr in pairs_all_loaded.items():
            pairs_storm[key] = arr[storm_mask]
            
        print_status(f"  Filtered {len(pairs_storm['dist']):,} pairs (Storm Days)", "INFO")
    else:
        day_files_storm = get_day_files_subset(days_storm)
        key_map = {
            'clock_bias': 'clock_bias_ns',
            'pos_jitter': 'pos_jitter',
            'clock_drift': 'clock_drift'
        }
        data_key = key_map.get(primary_observable, 'clock_bias_ns')
        pairs_storm = compute_pairs_coherence(day_files_storm, station_coords, "STORM", data_key)
        
    result_storm = analyze_anisotropy(pairs_storm, f"STORM DAYS (Kp>={KP_QUIET_THRESHOLD}) [{primary_observable}]")
    
    # ==========================================================================
    # STEP 7.5: HEMISPHERE STRATIFICATION (KEY FOR TEP EVIDENCE)
    # ==========================================================================
    print_status("[7.5/10] Running Hemisphere Stratification...", "PROCESS")
    
    hemisphere_results = {}
    if pairs_all_loaded and 'mid_lat' in pairs_all_loaded:
        # Northern Hemisphere (lat > 0)
        nh_mask = pairs_all_loaded['mid_lat'] > 0
        pairs_nh = {k: v[nh_mask] for k, v in pairs_all_loaded.items()}
        print_status(f"  Northern Hemisphere: {len(pairs_nh['dist']):,} pairs", "INFO")
        
        if len(pairs_nh['dist']) > 10000:
            result_nh_multi = run_multi_metric_anisotropy(pairs_nh, "NORTHERN HEMISPHERE")
            hemisphere_results['northern'] = result_nh_multi
        
        # Southern Hemisphere (lat < 0)
        sh_mask = pairs_all_loaded['mid_lat'] < 0
        pairs_sh = {k: v[sh_mask] for k, v in pairs_all_loaded.items()}
        print_status(f"  Southern Hemisphere: {len(pairs_sh['dist']):,} pairs", "INFO")
        
        if len(pairs_sh['dist']) > 10000:
            result_sh_multi = run_multi_metric_anisotropy(pairs_sh, "SOUTHERN HEMISPHERE")
            hemisphere_results['southern'] = result_sh_multi
    
    # ==========================================================================
    # STEP 8: MULTI-MODE COMPARATIVE ANALYSIS
    # ==========================================================================
    print_status("[8/10] Running Multi-Mode Comparative Analysis...", "PROCESS")
    
    mode_results = {'baseline': result_all}
    
    if modes_available['ionofree']:
        print_status("\n--- IONOFREE MODE (Option D: Dual-Freq Iono-Free) ---", "INFO")
        pair_file_ionofree = OUTPUTS_DIR / f"step_2_0_pairs_ionofree_{csv_filter_suffix}.csv"
        pairs_ionofree = load_pairs_from_csv(pair_file_ionofree, 'ionofree_clock_bias', list(station_coords.keys()), station_coords)
        
        if not pairs_ionofree:
            pairs_ionofree = compute_pairs_coherence(day_files, station_coords, "IONOFREE", 'ionofree_clock_bias_ns')
        
        if pairs_ionofree:
            mode_results['ionofree'] = analyze_anisotropy(pairs_ionofree, "IONOFREE (Dual-Freq) [ionofree_clock_bias]")
            mode_results['ionofree_multi_metric'] = run_multi_metric_anisotropy(pairs_ionofree, "IONOFREE")
        
    if modes_available['multi_gnss']:
        print_status("\n--- MULTI-GNSS MODE (MGEX) ---", "INFO")
        pair_file_mgex = OUTPUTS_DIR / f"step_2_0_pairs_multi_gnss_{csv_filter_suffix}.csv"
        pairs_mgex = load_pairs_from_csv(pair_file_mgex, 'mgex_clock_bias', list(station_coords.keys()), station_coords)
        
        if not pairs_mgex:
            pairs_mgex = compute_pairs_coherence(day_files, station_coords, "MGEX", 'mgex_clock_bias_ns')
        
        if pairs_mgex:
            mode_results['multi_gnss'] = analyze_anisotropy(pairs_mgex, "MULTI-GNSS (MGEX) [mgex_clock_bias]")
            mode_results['multi_gnss_multi_metric'] = run_multi_metric_anisotropy(pairs_mgex, "MULTI-GNSS")

    # PRECISE MODE (Option E: SPP with Precise Orbits)
    if modes_available.get('precise', False):
        print_status("\n--- PRECISE MODE (Option E: SPP with Precise Orbits) ---", "INFO")
        pair_file_precise = OUTPUTS_DIR / f"step_2_0_pairs_precise_{csv_filter_suffix}.csv"
        pairs_precise = load_pairs_from_csv(pair_file_precise, 'precise_clock_bias', list(station_coords.keys()), station_coords)
        
        if not pairs_precise:
            pairs_precise = compute_pairs_coherence(day_files, station_coords, "PRECISE", 'precise_clock_bias_ns')
        
        if pairs_precise:
            mode_results['precise'] = analyze_anisotropy(pairs_precise, "PRECISE (Precise Orbits) [precise_clock_bias]")
            mode_results['precise_multi_metric'] = run_multi_metric_anisotropy(pairs_precise, "PRECISE")

    # ==========================================================================
    # STEP 9: Generating Comparison Report
    # ==========================================================================
    print_status("[9/10] Generating Comparison Report...", "PROCESS")
    
    # Extract ratios for reporting
    ew_ns_all = result_all['ew_ns_pooled']['ratio'] if result_all and result_all.get('ew_ns_pooled') else 0
    ew_ns_quiet = result_quiet['ew_ns_pooled']['ratio'] if result_quiet and result_quiet.get('ew_ns_pooled') else 0
    ew_ns_storm = result_storm['ew_ns_pooled']['ratio'] if result_storm and result_storm.get('ew_ns_pooled') else 0
    
    mid_ratio_all = result_all['latitude_bands'].get('mid_lat_ratio', 0) if result_all else 0
    mid_ratio_quiet = result_quiet['latitude_bands'].get('mid_lat_ratio', 0) if result_quiet else 0
    
    # ==========================================================================
    # COMPREHENSIVE COMPARISON TABLE
    # ==========================================================================
    print_status("\n" + "="*80, "SUCCESS")
    print_status("COMPREHENSIVE COMPARISON TABLE", "SUCCESS")
    print_status("="*80, "SUCCESS")
    
    # Helper function to extract short-distance ratio
    def get_short_dist_ratio(res):
        if res and 'short_distance_analysis' in res:
            return res['short_distance_analysis'].get('mean_ratio')
        return None
    
    # Build comparison data
    comparison_rows = []
    
    # 1. By Kp condition
    print_status("\n  BY GEOMAGNETIC CONDITION (Kp):", "INFO")
    print_status("  " + "-"*70, "INFO")
    print_status(f"  {'Condition':<20} {'Short-Dist Ratio':<18} {'λ (Sector)':<12} {'λ (Pooled)':<12} {'N Pairs':<12}", "INFO")
    print_status("  " + "-"*70, "INFO")
    
    for name, res in [("All Days", result_all), ("Quiet (Kp<3)", result_quiet), ("Storm (Kp≥3)", result_storm)]:
        if res:
            sd_ratio = get_short_dist_ratio(res)
            lam_ratio_sector = res.get('ew_ns_code', {}).get('ratio')
            lam_ratio_pooled = res.get('ew_ns_pooled', {}).get('ratio')
            n_pairs = res.get('n_pairs', 0)
            sd_str = f"{sd_ratio:.3f}" if sd_ratio else "N/A"
            lam_sector_str = f"{lam_ratio_sector:.2f}" if lam_ratio_sector else "N/A"
            lam_pooled_str = f"{lam_ratio_pooled:.2f}" if lam_ratio_pooled else "N/A"
            print_status(f"  {name:<20} {sd_str:<18} {lam_sector_str:<12} {lam_pooled_str:<12} {n_pairs:,}", "INFO")
    
    # 2. By processing mode
    print_status("\n  BY PROCESSING MODE:", "INFO")
    print_status("  " + "-"*70, "INFO")
    print_status(f"  {'Mode':<20} {'Short-Dist Ratio':<18} {'λ (Sector)':<12} {'λ (Pooled)':<12} {'Status':<12}", "INFO")
    print_status("  " + "-"*70, "INFO")
    
    for name, key in [("Baseline (L1)", "baseline"), ("Precise (Orbits)", "precise"), ("Iono-Free (L1+L2)", "ionofree"), ("Multi-GNSS", "multi_gnss")]:
        res = mode_results.get(key)
        if res:
            sd_ratio = get_short_dist_ratio(res)
            lam_ratio_sector = res.get('ew_ns_code', {}).get('ratio')
            lam_ratio_pooled = res.get('ew_ns_pooled', {}).get('ratio')
            sd_str = f"{sd_ratio:.3f}" if sd_ratio else "N/A"
            lam_sector_str = f"{lam_ratio_sector:.2f}" if lam_ratio_sector else "N/A"
            lam_pooled_str = f"{lam_ratio_pooled:.2f}" if lam_ratio_pooled else "N/A"
            status = "✓" if sd_ratio and sd_ratio > 1.0 else "—"
            print_status(f"  {name:<20} {sd_str:<18} {lam_sector_str:<12} {lam_pooled_str:<12} {status}", "INFO")
        else:
            print_status(f"  {name:<20} {'—':<18} {'—':<12} {'—':<12} {'N/A'}", "INFO")
    
    # 3. By metric (from multi-metric analysis)
    if result_all_multi and result_all_multi.get('metric_results'):
        print_status("\n  BY PHASE METRIC:", "INFO")
        print_status("  " + "-"*70, "INFO")
        print_status(f"  {'Metric':<20} {'Short-Dist Ratio':<18} {'λ (Sector)':<12} {'λ (Pooled)':<12} {'CV':<12}", "INFO")
        print_status("  " + "-"*70, "INFO")
        
        for metric_name, metric_res in result_all_multi['metric_results'].items():
            sd_ratio = get_short_dist_ratio(metric_res)
            lam_ratio_sector = metric_res.get('ew_ns_code', {}).get('ratio')
            lam_ratio_pooled = metric_res.get('ew_ns_pooled', {}).get('ratio')
            cv = metric_res.get('anisotropy_cv', 0)
            sd_str = f"{sd_ratio:.3f}" if sd_ratio else "N/A"
            lam_sector_str = f"{lam_ratio_sector:.2f}" if lam_ratio_sector else "N/A"
            lam_pooled_str = f"{lam_ratio_pooled:.2f}" if lam_ratio_pooled else "N/A"
            print_status(f"  {metric_name:<20} {sd_str:<18} {lam_sector_str:<12} {lam_pooled_str:<12} {cv:.3f}", "INFO")
    
    # 4. By hemisphere
    if hemisphere_results:
        print_status("\n  BY HEMISPHERE:", "INFO")
        print_status("  " + "-"*70, "INFO")
        print_status(f"  {'Hemisphere':<20} {'Short-Dist Ratio':<18} {'λ (Sector)':<12} {'λ (Pooled)':<12} {'N Pairs':<12}", "INFO")
        print_status("  " + "-"*70, "INFO")
        
        for name, key in [("Northern", "northern"), ("Southern", "southern")]:
            res = hemisphere_results.get(key, {})
            if res and res.get('metric_results'):
                # Use phase_alignment results if available
                pa_res = res['metric_results'].get('phase_alignment', {})
                sd_ratio = get_short_dist_ratio(pa_res)
                lam_ratio_sector = pa_res.get('ew_ns_code', {}).get('ratio')
                lam_ratio_pooled = pa_res.get('ew_ns_pooled', {}).get('ratio')
                n_pairs = pa_res.get('n_pairs', 0)
                sd_str = f"{sd_ratio:.3f}" if sd_ratio else "N/A"
                lam_sector_str = f"{lam_ratio_sector:.2f}" if lam_ratio_sector else "N/A"
                lam_pooled_str = f"{lam_ratio_pooled:.2f}" if lam_ratio_pooled else "N/A"
                print_status(f"  {name:<20} {sd_str:<18} {lam_sector_str:<12} {lam_pooled_str:<12} {n_pairs:,}", "INFO")
    
    print_status("\n  CODE Reference: λ ratio = 2.16 (E-W/N-S from 25yr PPP)", "INFO")
    print_status("  " + "="*70, "INFO")
    
    # Save results
    output = {
        'timestamp': datetime.now().isoformat(),
        'filter_mode': filter_mode_name,
        'methodology': 'CODE Longspan-Aligned Anisotropy with Multi-Metric Comparison',
        'data_years': years_present,
        'methodology_details': {
            'sectors': 8,
            'sector_width_deg': 45,
            'distance_bins': N_BINS,
            'min_distance_km': MIN_DISTANCE_KM,
            'max_distance_km': MAX_DISTANCE_KM,
            'min_bin_count': MIN_BIN_COUNT,
            'metrics_compared': list(PHASE_METRICS.keys())
        },
        'kp_threshold': KP_QUIET_THRESHOLD,
        'code_reference_ratio': CODE_EW_NS_RATIO,
        'multi_metric_analysis': result_all_multi,  # Multi-metric comparison (ALL DAYS)
        'quiet_days_multi_metric': result_quiet_multi,  # Multi-metric comparison (QUIET DAYS)
        'hemisphere_analysis': hemisphere_results,  # Hemisphere stratification
        'baseline': {
            'all_days': result_all,
            'quiet_days': result_quiet,
            'storm_days': result_storm
        },
        'primary_observable': primary_observable,
        'control_observables': control_observables,
        'baseline_controls': control_results,
        'baseline_controls_multi_metric': control_results_multi,
        'modes_available': modes_available,
        'ionofree_results': mode_results.get('ionofree'),
        'ionofree_multi_metric': mode_results.get('ionofree_multi_metric'),
        'mgex_results': mode_results.get('multi_gnss'),
        'mgex_multi_metric': mode_results.get('multi_gnss_multi_metric'),
        'precise_results': mode_results.get('precise'),
        'precise_multi_metric': mode_results.get('precise_multi_metric'),
        'comparison_summary': {
            'by_condition': {
                'all_days': get_short_dist_ratio(result_all) if result_all else None,
                'quiet_days': get_short_dist_ratio(result_quiet) if result_quiet else None,
                'storm_days': get_short_dist_ratio(result_storm) if result_storm else None,
            },
            'by_mode': {
                'baseline': get_short_dist_ratio(mode_results.get('baseline')),
                'precise': get_short_dist_ratio(mode_results.get('precise')),
                'ionofree': get_short_dist_ratio(mode_results.get('ionofree')),
                'multi_gnss': get_short_dist_ratio(mode_results.get('multi_gnss')),
            }
        }
    }
    
    # ==========================================================================
    # ADD EXECUTIVE SUMMARY TO OUTPUT
    # ==========================================================================
    short_dist = result_all.get('short_distance_analysis', {}) if result_all else {}
    
    executive_summary = {
        'tep_signal_detected': short_dist.get('tep_status') in ['DETECTED', 'LIKELY'],
        'primary_evidence': {
            'test': 'Short-distance (<500km) E-W/N-S mean ratio',
            'result': short_dist.get('mean_ratio'),
            'ci_95': [short_dist.get('ci_95_low'), short_dist.get('ci_95_high')],
            'p_value': short_dist.get('p_value'),
            'interpretation': 'E-W correlation stronger than N-S at short distances'
        },
        'comparison_to_code': {
            'code_lambda_ratio': CODE_EW_NS_RATIO,
            'note': 'CODE uses λ (correlation length) ratio from 25yr PPP data; we use mean correlation ratio from multi-year SPP data'
        },
        'key_findings': [
            f"Short-distance E-W/N-S ratio: {short_dist.get('mean_ratio', 'N/A'):.3f}" if short_dist.get('mean_ratio') else "Short-distance analysis not available",
            f"95% CI: [{short_dist.get('ci_95_low', 0):.3f}, {short_dist.get('ci_95_high', 0):.3f}]" if short_dist.get('ci_95_low') else "",
            f"Statistical significance: p = {short_dist.get('p_value', 1):.2e}" if short_dist.get('p_value') else "",
            f"Effect size (Cohen's d): {short_dist.get('cohens_d', 0):.3f}" if short_dist.get('cohens_d') else "",
        ],
        'methodology_verified': [
            'Azimuth: spherical forward azimuth (geodetically correct)',
            f"Primary observable: {primary_observable}",
            f"Controls: {', '.join(control_observables) if control_observables else 'none'}",
            'Sectors: 8 x 45° (CODE-aligned)',
            'Frequency band: 10-500 µHz (CODE-aligned)',
        ]
    }
    
    output['executive_summary'] = executive_summary
    
    # ==========================================================================
    # MONTHLY ANISOTROPY ANALYSIS - Temporal stratification (ALL MODES)
    # ==========================================================================
    # This addresses the insight that E-W/N-S ratio varies monthly (CMB-coupled)
    # and that aggregating all months may mask the true temporal pattern
    print_status("\n[8/8] Computing monthly anisotropy stratification...", "PROCESS")
    
    monthly_by_mode = {}
    
    # Helper to print monthly summary
    def print_monthly_summary(name, monthly_result):
        if not monthly_result:
            print_status(f"  {name}: No data", "WARNING")
            return
        
        coh_summary = monthly_result.get('summary_coherence', {})
        phase_summary = monthly_result.get('summary_phase_alignment', {})
        
        if coh_summary:
            n = coh_summary.get('n_months', 0)
            mean_r = coh_summary.get('mean_ratio', 0)
            pct = coh_summary.get('pct_ew_dominant', 0)
            print_status(f"  {name} - Coherence: mean={mean_r:.4f}, E-W>{n*pct/100:.0f}/{n} months ({pct:.0f}%)", "INFO")
        
        if phase_summary:
            n = phase_summary.get('n_months', 0)
            mean_r = phase_summary.get('mean_ratio', 0)
            pct = phase_summary.get('pct_ew_dominant', 0)
            print_status(f"  {name} - Phase Align: mean={mean_r:.4f}, E-W>{n*pct/100:.0f}/{n} months ({pct:.0f}%)", "INFO")
    
    # Baseline mode
    if pairs_all_loaded:
        monthly_baseline = compute_monthly_anisotropy_from_pairs(pairs_all_loaded, station_coords)
        if monthly_baseline:
            monthly_by_mode['baseline'] = monthly_baseline
            print_monthly_summary("Baseline", monthly_baseline)
    
    # Ionofree mode - check if pairs were loaded earlier
    try:
        if pairs_ionofree and isinstance(pairs_ionofree, dict) and 'dist' in pairs_ionofree:
            monthly_ionofree = compute_monthly_anisotropy_from_pairs(pairs_ionofree, station_coords)
            if monthly_ionofree:
                monthly_by_mode['ionofree'] = monthly_ionofree
                print_monthly_summary("Ionofree", monthly_ionofree)
    except NameError:
        pass  # pairs_ionofree not defined
    
    # Multi-GNSS mode
    try:
        if pairs_mgex and isinstance(pairs_mgex, dict) and 'dist' in pairs_mgex:
            monthly_mgex = compute_monthly_anisotropy_from_pairs(pairs_mgex, station_coords)
            if monthly_mgex:
                monthly_by_mode['multi_gnss'] = monthly_mgex
                print_monthly_summary("Multi-GNSS", monthly_mgex)
    except NameError:
        pass  # pairs_mgex not defined
    
    # Precise mode
    try:
        if pairs_precise and isinstance(pairs_precise, dict) and 'dist' in pairs_precise:
            monthly_precise = compute_monthly_anisotropy_from_pairs(pairs_precise, station_coords)
            if monthly_precise:
                monthly_by_mode['precise'] = monthly_precise
                print_monthly_summary("Precise", monthly_precise)
    except NameError:
        pass  # pairs_precise not defined
    
    # Store in output
    if monthly_by_mode:
        output['monthly_anisotropy'] = monthly_by_mode

        monthly_inference = {}
        for mode_name, mode_data in monthly_by_mode.items():
            monthly_inference[mode_name] = {}
            for metric_name in ['coherence', 'phase_alignment']:
                ratios = [v.get('ratio') for v in (mode_data.get(metric_name, {}) or {}).values() if v and v.get('ratio')]
                if len(ratios) < 6:
                    continue
                ratios = np.array(ratios, dtype=np.float64)
                from scipy import stats
                t_stat, p_value = stats.ttest_1samp(ratios, 1.0)
                np.random.seed(42)
                boot = np.mean(np.random.choice(ratios, size=(2000, len(ratios)), replace=True), axis=1)
                ci_low, ci_high = np.percentile(boot, [2.5, 97.5])
                monthly_inference[mode_name][metric_name] = {
                    'n_months': int(len(ratios)),
                    'mean_ratio': float(np.mean(ratios)),
                    'std_ratio': float(np.std(ratios)),
                    'ci_95_low': float(ci_low),
                    'ci_95_high': float(ci_high),
                    't_statistic': float(t_stat),
                    'p_value': float(p_value),
                }

        if any(monthly_inference.get(m) for m in monthly_inference):
            output['monthly_inference'] = monthly_inference
            print_status("\n  MONTHLY INFERENCE (MONTHS AS SAMPLES):", "SUCCESS")
            print_status("  " + "-"*75, "INFO")
            print_status(f"  {'Mode':<12} {'Metric':<15} {'Mean':>10} {'CI95 Low':>10} {'CI95 High':>10} {'p':>10}", "INFO")
            print_status("  " + "-"*75, "INFO")
            for mode_name in ['baseline', 'precise', 'ionofree', 'multi_gnss']:
                for metric_name in ['coherence', 'phase_alignment']:
                    inf = monthly_inference.get(mode_name, {}).get(metric_name)
                    if not inf:
                        continue
                    metric_label = 'Coherence' if metric_name == 'coherence' else 'Phase Align'
                    print_status(
                        f"  {mode_name:<12} {metric_label:<15} {inf['mean_ratio']:>10.4f} {inf['ci_95_low']:>10.4f} {inf['ci_95_high']:>10.4f} {inf['p_value']:>10.2e}",
                        "INFO",
                    )
            print_status("  " + "-"*75, "INFO")
        
        # Overall summary table
        print_status("\n  MONTHLY ANISOTROPY SUMMARY TABLE:", "SUCCESS")
        print_status("  " + "-"*75, "INFO")
        print_status(f"  {'Mode':<12} {'Metric':<15} {'Mean Ratio':>12} {'Std':>8} {'E-W>N-S %':>12} {'Months':>8}", "INFO")
        print_status("  " + "-"*75, "INFO")
        
        for mode_name in ['baseline', 'precise', 'ionofree', 'multi_gnss']:
            mode_data = monthly_by_mode.get(mode_name, {})
            for metric in ['summary_coherence', 'summary_phase_alignment']:
                s = mode_data.get(metric, {})
                if s:
                    metric_label = 'Coherence' if 'coh' in metric else 'Phase Align'
                    print_status(f"  {mode_name:<12} {metric_label:<15} {s['mean_ratio']:>12.4f} {s['std_ratio']:>8.4f} {s['pct_ew_dominant']:>11.0f}% {s['n_months']:>8}", "INFO")
        
        print_status("  " + "-"*75, "INFO")
    else:
        print_status("  Monthly anisotropy not computed (insufficient data)", "WARNING")
    
    # Custom JSON encoder to handle numpy types
    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.integer, np.int64, np.int32)):
                return int(obj)
            if isinstance(obj, (np.floating, np.float64, np.float32)):
                return float(obj)
            if isinstance(obj, (np.bool_, bool)):
                return bool(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return super().default(obj)
    
    with open(output_file, "w") as f:
        json.dump(output, f, indent=2, cls=NumpyEncoder)
        
    print_status(f"\nResults saved: {output_file}", "SUCCESS")
    
    # ==========================================================================
    # PRINT EXECUTIVE SUMMARY
    # ==========================================================================
    print_status("", "INFO")
    print_status("=" * 80, "SUCCESS")
    print_status("EXECUTIVE SUMMARY", "SUCCESS")
    print_status("=" * 80, "SUCCESS")
    
    if short_dist:
        tep_status = short_dist.get('tep_status', 'UNKNOWN')
        ratio = short_dist.get('mean_ratio')
        ci_low = short_dist.get('ci_95_low')
        ci_high = short_dist.get('ci_95_high')
        p_val = short_dist.get('p_value')
        
        print_status(f"\n  PRIMARY RESULT: TEP Signal {tep_status}", "SUCCESS" if tep_status in ['DETECTED', 'LIKELY'] else "WARNING")
        print_status(f"\n  Short-Distance (<500km) E-W/N-S Ratio: {ratio:.3f}" if ratio else "", "SUCCESS")
        print_status(f"  95% Confidence Interval: [{ci_low:.3f}, {ci_high:.3f}]" if ci_low else "", "INFO")
        print_status(f"  p-value: {p_val:.2e}" if p_val else "", "INFO")
        
        if ratio and ratio > 1.0:
            print_status(f"\n  INTERPRETATION:", "INFO")
            print_status(f"    E-W correlations are {(ratio-1)*100:.1f}% stronger than N-S at short distances", "SUCCESS")
            print_status(f"    This is consistent with TEP (Earth rotation coupling)", "SUCCESS")
        
        print_status(f"\n  CODE Reference:", "INFO")
        print_status(f"    CODE λ ratio (25yr PPP): {CODE_EW_NS_RATIO}", "INFO")
        print_status(f"    Note: Different metric (λ vs mean), different precision (PPP vs SPP)", "INFO")
    
    print_status("", "INFO")
    print_status("=" * 80, "SUCCESS")
    print_status("ANALYSIS COMPLETE", "SUCCESS")
    print_status("=" * 80, "SUCCESS")
    
def main():
    """Main execution loop over all filter modes."""
    import argparse
    
    parser = argparse.ArgumentParser(description='TEP-GNSS-RINEX Step 2.2: Anisotropy Analysis')
    parser.add_argument('--filter', type=str, default='all', 
                        help='Filter to run: "all", "none", "optimal_100_metadata.json", or "dynamic:50"')
    parser.add_argument(
        '--observables',
        type=str,
        default='clock_bias',
        help='Comma-separated list of observables to analyze from Step 2.0 CSV metric column. '
             'First is treated as primary; remaining are controls. Use "all" for clock_bias,pos_jitter,clock_drift.'
    )
    args = parser.parse_args()
    
    # Define the 3 standard filters
    all_filters = [
        ('ALL STATIONS', 'none'),
        ('OPTIMAL 100', 'optimal_100_metadata.json'),
        ('DYNAMIC 50', 'dynamic:50')
    ]
    
    # Select filters to run
    if args.filter == 'all':
        filters = all_filters
        print_status("STARTING MULTI-FILTER ANALYSIS PIPELINE", "TITLE")
    else:
        # Find matching filter
        matching = [f for f in all_filters if f[1] == args.filter]
        if matching:
            filters = matching
            print_status(f"RUNNING SINGLE FILTER: {matching[0][0]}", "TITLE")
        else:
            print_status(f"Unknown filter: {args.filter}", "ERROR")
            print_status("Valid filters: none, optimal_100_metadata.json, dynamic:50, all", "INFO")
            return
    
    for name, config in filters:
        try:
            run_analysis_pipeline(name, config, observables=args.observables)
        except Exception as e:
            print_status(f"Analysis failed for {name}: {e}", "ERROR")
            import traceback
            traceback.print_exc()
            
    print_status("ANALYSIS COMPLETE", "SUCCESS")

if __name__ == "__main__":
    main()
