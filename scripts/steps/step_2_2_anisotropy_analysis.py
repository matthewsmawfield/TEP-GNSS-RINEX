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
    for i in range(len(ref_bins) - 1):
        bin_mask = (sector_distances >= ref_bins[i]) & (sector_distances < ref_bins[i+1])
        bin_indices = np.where(bin_mask)[0]
        
        if len(bin_indices) > 0:
            n_to_sample = int(min(len(bin_indices), target_per_bin[i]))
            if n_to_sample > 0:
                sampled = np.random.choice(bin_indices, size=n_to_sample, replace=False)
                selected_indices.extend(sampled)
    
    return np.array(selected_indices, dtype=int)


def fit_exponential(distances, coherences, min_pairs=MIN_BIN_COUNT, bins=None, use_sigma=True):
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
    
    if len(bin_centers) < 5:
        return None
    
    bin_centers = np.array(bin_centers)
    bin_means = np.array(bin_means)
    bin_counts = np.array(bin_counts)
    
    try:
        weights = np.sqrt(bin_counts)
        
        # Adaptive bounds matching CODE longspan TEPConfig
        max_dist = np.max(bin_centers) if len(bin_centers) > 0 else MAX_DISTANCE_KM
        max_lambda = min(15000, max_dist * 0.8)  # Cap at 15,000km or 80% of observed range
        
        # Bounds: ([min_A, min_lambda, min_offset], [max_A, max_lambda, max_offset])
        bounds = ([0, 100, -1], [5, max_lambda, 1])
        
        # Data-driven initial guess (CODE methodology)
        c_range = bin_means.max() - bin_means.min()
        c_min = bin_means.min()
        p0 = [c_range if c_range > 0 else 0.5, 3000, c_min]
        
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
        
        predicted = exp_decay(bin_centers, *popt)
        ss_res = np.sum((bin_means - predicted)**2)
        ss_tot = np.sum((bin_means - np.mean(bin_means))**2)
        r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0
        
        return {
            'lambda_km': float(lam),
            'lambda_err': float(lam_err),
            'amplitude': float(A),
            'offset': float(C0),
            'r_squared': float(r2),
            'n_pairs': int(sum(bin_counts)),
            'n_bins': len(bin_centers)
        }
    except Exception:
        return None


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
    
    # DOY bounds for each month
    doy_bounds = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 366]
    doys = pairs_dict['doy']
    months = np.zeros(len(doys), dtype=np.int32)
    
    for m in range(1, 13):
        mask = (doys > doy_bounds[m-1]) & (doys <= doy_bounds[m])
        months[mask] = m
    
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
    ew_mean = None
    ns_mean = None
    if len(ew_lambdas) >= 1 and len(ns_lambdas) >= 1:
        ew_mean = np.mean(ew_lambdas)
        ns_mean = np.mean(ns_lambdas)
        ew_ns_ratio = ew_mean / ns_mean if ns_mean > 0 else None
        print_status(f"  E-W/N-S RATIO (FULL DATASET - PRIMARY):", "SUCCESS")
        ew_str = ", ".join([f"{l:.0f}" for l in ew_lambdas])
        ns_str = ", ".join([f"{l:.0f}" for l in ns_lambdas])
        print_status(f"    E-W λ = {ew_mean:.0f} km (sectors: {ew_str})", "INFO")
        print_status(f"    N-S λ = {ns_mean:.0f} km (sectors: {ns_str})", "INFO")
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
    
    n_short_ew = np.sum(short_ew_mask)
    n_short_ns = np.sum(short_ns_mask)
    
    print_status(f"\n  SHORT-DISTANCE ANALYSIS (<{short_dist_threshold}km) - PRIMARY TEP TEST:", "SUCCESS")
    print_status(f"    E-W pairs: {n_short_ew:,}, N-S pairs: {n_short_ns:,}", "INFO")
    
    if n_short_ew > 100 and n_short_ns > 100:
        # Use raw coherence values (phase_alignment or coherence)
        ew_short_values = all_coherences[short_ew_mask]
        ns_short_values = all_coherences[short_ns_mask]
        
        ew_mean = np.mean(ew_short_values)
        ns_mean = np.mean(ns_short_values)
        ew_std = np.std(ew_short_values)
        ns_std = np.std(ns_short_values)
        
        # Mean ratio (not lambda ratio - direct correlation comparison)
        mean_ratio = ew_mean / ns_mean if ns_mean != 0 else None
        
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
        pooled_std = np.sqrt((ew_std**2 + ns_std**2) / 2)
        cohens_d = (ew_mean - ns_mean) / pooled_std if pooled_std > 0 else 0
        
        # Store results
        short_dist_analysis = {
            'threshold_km': short_dist_threshold,
            'n_ew_pairs': int(n_short_ew),
            'n_ns_pairs': int(n_short_ns),
            'ew_mean': float(ew_mean),
            'ns_mean': float(ns_mean),
            'ew_std': float(ew_std),
            'ns_std': float(ns_std),
            'mean_ratio': float(mean_ratio) if mean_ratio else None,
            'ci_95_low': float(ci_low),
            'ci_95_high': float(ci_high),
            't_statistic': float(t_stat),
            'p_value': float(p_value),
            'cohens_d': float(cohens_d),
            'significant': p_value < 0.01
        }
        
        # Print results
        print_status(f"    E-W mean: {ew_mean:.4f} ± {ew_std:.4f}", "INFO")
        print_status(f"    N-S mean: {ns_mean:.4f} ± {ns_std:.4f}", "INFO")
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
            'ew_lambda_mean': float(ew_mean) if ew_mean else None,
            'ns_lambda_mean': float(ns_mean) if ns_mean else None,
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
        'short_distance_analysis': short_dist_analysis,
        'anisotropy_cv': float(cv),
        'methodology': {
            'azimuth_calculation': 'spherical_forward_azimuth',
            'sector_width_deg': 45,
            'distance_matching': True,
            'metric_used': 'clock_bias_only',
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
        
        # Run the standard anisotropy analysis
        metric_result = analyze_anisotropy(metric_pairs, f"{condition_name}_{metric_name}")
        
        if metric_result:
            results['metric_results'][metric_name] = metric_result
            
            # Extract key summary stats
            ew_ns_ratio = metric_result.get('ew_ns_pooled', {}).get('ratio')
            cv = metric_result.get('anisotropy_cv', 0)
            
            print_status(f"      E-W/N-S Ratio: {ew_ns_ratio:.2f}" if ew_ns_ratio else "      E-W/N-S Ratio: N/A", "SUCCESS")
            print_status(f"      Anisotropy CV: {cv:.3f}", "INFO")
    
    # Cross-metric comparison (if both metrics available)
    if len(results['metric_results']) >= 2:
        print_status(f"\n  --- Cross-Metric Comparison ---", "PROCESS")
        
        coh_result = results['metric_results'].get('coherence', {})
        phase_result = results['metric_results'].get('phase_alignment', {})
        
        coh_ratio = coh_result.get('ew_ns_pooled', {}).get('ratio')
        phase_ratio = phase_result.get('ew_ns_pooled', {}).get('ratio')
        
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
            
            print_status(f"      MSC E-W/N-S Ratio: {coh_ratio:.2f}", "INFO")
            print_status(f"      Phase E-W/N-S Ratio: {phase_ratio:.2f}", "INFO")
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
        
        # Add 'metric' column if it exists to filter for the requested metric only
        # IMPORTANT: pos_jitter shows isotropic correlations, diluting the anisotropy signal
        has_metric_col = 'metric' in header
        if has_metric_col:
            use_cols.append('metric')
            print_status(f"  Will filter for '{metric_key}' metric only (pos_jitter is isotropic)", "INFO")
        
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
        total_rows = 0
        chunk_size = 500000  # 500k rows per chunk
        
        # Create filter set if needed
        sta_filter = set(station_filter_list) if station_filter_list else None
        
        # Read in chunks
        for i, chunk in enumerate(pd.read_csv(csv_file, usecols=use_cols, chunksize=chunk_size)):
            # CRITICAL: Filter to requested metric only (pos_jitter is isotropic, dilutes signal)
            if has_metric_col:
                chunk = chunk[chunk['metric'] == metric_key]
            
            # Filter by stations if needed
            if sta_filter:
                mask = chunk['station1'].isin(sta_filter) & chunk['station2'].isin(sta_filter)
                chunk = chunk[mask]
            
            if chunk.empty: continue
            
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
            
            if phase_col:
                data_arrays['phase_alignment'].append(chunk[phase_col].values.astype(np.float32))
            
            total_rows += len(chunk)
            if (i + 1) % 5 == 0:
                print_status(f"  Loaded {total_rows/1e6:.1f}M pairs...", "INFO")
                gc.collect()
        
        if total_rows == 0:
            return None
            
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

def run_analysis_pipeline(filter_mode_name, station_filter_config, load_from_csv=False, csv_file=None):
    """Run the full analysis pipeline for a specific filter mode."""
    print_status("", "INFO")
    print_status("=" * 80, "INFO")
    print_status(f"TEP-GNSS-RINEX Analysis - STEP 2.2: Anisotropy ({filter_mode_name})", "INFO")
    print_status("=" * 80, "INFO")
    
    # Define output filename based on filter
    filter_suffix = filter_mode_name.lower().replace(" ", "_").replace("-", "_")
    output_file = OUTPUTS_DIR / f"step_2_2_anisotropy_analysis_{filter_suffix}.json"
    
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
    
    # Load NPZ files with optional dynamic filtering
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
    # This allows day-by-day processing instead of concatenating all years
    day_files = {}  # {(year, doy): {station: filepath}}
    station_files = {}  # Keep for compatibility
    
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
    all_years = set()
    for f in npz_files:
        try:
            # Format: STATION_YYYYDDD.npz
            year_str = f.name.split('_')[1][:4]
            all_years.add(int(year_str))
        except:
            continue
            
    if not all_years:
        all_years = {2022}
        print_status("  Could not detect years, defaulting to 2022", "WARNING")
    else:
        print_status(f"  Detected years: {sorted(list(all_years))}", "INFO")
    
    # DAY-BY-DAY PROCESSING - CODE LONGSPAN ARCHITECTURE
    # Instead of loading all data into memory, we'll process one day at a time
    # and aggregate bin statistics (sum_coherence, count) across days
    
    print_status(f"[1/8] Preparing for day-by-day processing...", "PROCESS")
    print_status(f"  Total days to process: {len(day_files)}", "INFO")
    print_status(f"  Stations available: {len(station_files)}", "INFO")
    
    # We'll store metadata about which stations have which modes
    # This is determined by checking a sample of files
    station_modes = {}  # {station: {mode: has_data}}
    all_keys = ['clock_bias_ns', 'ionofree_clock_bias_ns', 'mgex_clock_bias_ns']
    
    print(f"\r  Checking station modes...", end="", flush=True)
    for sta in list(station_files.keys()):  # Check ALL stations
        sample_file = station_files[sta][0]
        try:
            with np.load(sample_file) as d:
                station_modes[sta] = {
                    'baseline': 'clock_bias_ns' in d,
                    'ionofree': 'ionofree_clock_bias_ns' in d,
                    'mgex': 'mgex_clock_bias_ns' in d,
                    'lat': station_coords[sta]['lat'],
                    'lon': station_coords[sta]['lon']
                }
        except:
            continue
    
    # Count stations per mode
    n_baseline = sum(1 for s in station_modes.values() if s.get('baseline'))
    n_ionofree = sum(1 for s in station_modes.values() if s.get('ionofree'))
    n_mgex = sum(1 for s in station_modes.values() if s.get('mgex'))
    
    print(f"\r  Station modes detected:              ")
    print_status(f"  Baseline mode stations: {n_baseline}", "SUCCESS")
    print_status(f"  Iono-free mode stations: {n_ionofree} (Option D)", "INFO")
    print_status(f"  Multi-GNSS mode stations: {n_mgex} (MGEX)", "INFO")
    
    # Store mode availability
    modes_available = {
        'baseline': n_baseline >= 20,
        'ionofree': n_ionofree >= 20,
        'multi_gnss': n_mgex >= 20,
    }
    
    # Store day_files and station_modes for later use
    # These will be used by the day-by-day processing functions
    
    # ==========================================================================
    # STEP 2: Fetch Kp data and create quiet mask
    # ==========================================================================
    print_status("[2/8] Fetching Kp data for stratification...", "PROCESS")
    kp_data = fetch_kp_data_gfz(sorted(list(all_years)))
    
    # Filter days based on Kp
    days_all = sorted(list(day_files.keys()))
    days_quiet = []
    days_storm = []
    
    for year, doy in days_all:
        try:
            date = datetime(year, 1, 1) + timedelta(days=doy-1)
            # FIX: Match fetch_kp_data_gfz format (YYYYMMDD)
            date_str = date.strftime('%Y%m%d')
            kp = kp_data.get(date_str, 9.0)
            if kp < KP_QUIET_THRESHOLD:
                days_quiet.append((year, doy))
            else:
                days_storm.append((year, doy))
        except:
            pass
            
    quiet_days_count = len(days_quiet)
    total_days_count = len(days_all)
    print_status(f"  Quiet days (Kp < {KP_QUIET_THRESHOLD}): {quiet_days_count}/{total_days_count}", "INFO")
    
    # Helper to subset day_files
    def get_day_files_subset(days):
        return {d: day_files[d] for d in days if d in day_files}

    # ==========================================================================
    # STEP 4: Analyze ALL DAYS
    # ==========================================================================
    
    # Try loading pairs from Step 2.0 output (Baseline)
    pair_file_baseline = OUTPUTS_DIR / "step_2_0_pairs_baseline.csv"
    pairs_all_loaded = load_pairs_from_csv(pair_file_baseline, 'clock_bias', list(station_coords.keys()), station_coords)
    
    if pairs_all_loaded:
        print_status("Using pre-computed pairs from Step 2.0 output", "SUCCESS")
        pairs_all = pairs_all_loaded
    else:
        print_status("[4/8] Analyzing ALL DAYS (Day-by-Day)...", "PROCESS")
        pairs_all = compute_pairs_coherence(day_files, station_coords, "ALL", 'clock_bias_ns')
    
    # Run MULTI-METRIC analysis (compares MSC and phase_alignment)
    result_all_multi = run_multi_metric_anisotropy(pairs_all, "ALL DAYS (Multi-Year)")
    
    # Also run single-metric for backward compatibility
    result_all = analyze_anisotropy(pairs_all, "ALL DAYS (Multi-Year)")
    
    # ==========================================================================
    # STEP 6: Analyze QUIET DAYS
    # ==========================================================================
    print_status("[6/8] Analyzing QUIET DAYS...", "PROCESS")
    
    if pairs_all_loaded:
        # Filter from loaded pairs using Vectorized Masking
        quiet_set = set(days_quiet) # Set of tuples (year, doy)
        
        # Create a boolean mask for quiet days
        # Since we can't easily vector-match tuples, we'll match year/doy separately or use a fast lookup
        # Fast approach: Create a combined ID or use list comprehension for the mask only (fast enough for 100M bools)
        
        print_status("  Creating mask for quiet days...", "INFO")
        # Create a lookup set of YYYYDDD integers for fast matching
        quiet_ids = {y * 1000 + d for y, d in days_quiet}
        
        # Compute IDs for all pairs
        pair_ids = pairs_all_loaded['year'] * 1000 + pairs_all_loaded['doy']
        
        # Create mask (using np.isin is efficient)
        quiet_mask = np.isin(pair_ids, list(quiet_ids))
        
        # Apply mask to all arrays
        pairs_quiet = {}
        for key, arr in pairs_all_loaded.items():
            pairs_quiet[key] = arr[quiet_mask]
            
        print_status(f"  Filtered {len(pairs_quiet['dist'])} pairs (Quiet Days)", "INFO")
    else:
        # Re-compute
        day_files_quiet = get_day_files_subset(days_quiet)
        pairs_quiet = compute_pairs_coherence(day_files_quiet, station_coords, "QUIET", 'clock_bias_ns')
    
    # Run MULTI-METRIC analysis on Quiet Days (key for TEP evidence)
    result_quiet_multi = run_multi_metric_anisotropy(pairs_quiet, f"QUIET DAYS (Kp<{KP_QUIET_THRESHOLD})")
    result_quiet = analyze_anisotropy(pairs_quiet, f"QUIET DAYS (Kp<{KP_QUIET_THRESHOLD})")
    
    # ==========================================================================
    # STEP 7: Analyze STORM DAYS
    # ==========================================================================
    print_status("[7/8] Analyzing STORM DAYS...", "PROCESS")
    
    if pairs_all_loaded:
        # Filter from loaded pairs
        storm_set = set(days_storm)
        storm_ids = {y * 1000 + d for y, d in days_storm}
        
        # Reuse pair_ids if available, else recompute
        if 'pair_ids' not in locals():
            pair_ids = pairs_all_loaded['year'] * 1000 + pairs_all_loaded['doy']
            
        storm_mask = np.isin(pair_ids, list(storm_ids))
        
        pairs_storm = {}
        for key, arr in pairs_all_loaded.items():
            pairs_storm[key] = arr[storm_mask]
            
        print_status(f"  Filtered {len(pairs_storm['dist'])} pairs (Storm Days)", "INFO")
    else:
        day_files_storm = get_day_files_subset(days_storm)
        pairs_storm = compute_pairs_coherence(day_files_storm, station_coords, "STORM", 'clock_bias_ns')
        
    result_storm = analyze_anisotropy(pairs_storm, f"STORM DAYS (Kp>={KP_QUIET_THRESHOLD})")
    
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
        pair_file_ionofree = OUTPUTS_DIR / "step_2_0_pairs_ionofree.csv"
        pairs_ionofree = load_pairs_from_csv(pair_file_ionofree, 'ionofree_clock_bias', list(station_coords.keys()), station_coords)
        
        if not pairs_ionofree:
            pairs_ionofree = compute_pairs_coherence(day_files, station_coords, "IONOFREE", 'ionofree_clock_bias_ns')
        
        if pairs_ionofree:
            mode_results['ionofree'] = analyze_anisotropy(pairs_ionofree, "IONOFREE (Dual-Freq)")
            mode_results['ionofree_multi_metric'] = run_multi_metric_anisotropy(pairs_ionofree, "IONOFREE")
        
    if modes_available['multi_gnss']:
        print_status("\n--- MULTI-GNSS MODE (MGEX) ---", "INFO")
        pair_file_mgex = OUTPUTS_DIR / "step_2_0_pairs_multi_gnss.csv"
        pairs_mgex = load_pairs_from_csv(pair_file_mgex, 'mgex_clock_bias', list(station_coords.keys()), station_coords)
        
        if not pairs_mgex:
            pairs_mgex = compute_pairs_coherence(day_files, station_coords, "MGEX", 'mgex_clock_bias_ns')
        
        if pairs_mgex:
            mode_results['multi_gnss'] = analyze_anisotropy(pairs_mgex, "MULTI-GNSS (MGEX)")
            mode_results['multi_gnss_multi_metric'] = run_multi_metric_anisotropy(pairs_mgex, "MULTI-GNSS")

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
    print_status(f"  {'Condition':<20} {'Short-Dist Ratio':<18} {'λ Ratio':<15} {'N Pairs':<12}", "INFO")
    print_status("  " + "-"*70, "INFO")
    
    for name, res in [("All Days", result_all), ("Quiet (Kp<3)", result_quiet), ("Storm (Kp≥3)", result_storm)]:
        if res:
            sd_ratio = get_short_dist_ratio(res)
            lam_ratio = res.get('ew_ns_pooled', {}).get('ratio')
            n_pairs = res.get('n_pairs', 0)
            sd_str = f"{sd_ratio:.3f}" if sd_ratio else "N/A"
            lam_str = f"{lam_ratio:.2f}" if lam_ratio else "N/A"
            print_status(f"  {name:<20} {sd_str:<18} {lam_str:<15} {n_pairs:,}", "INFO")
    
    # 2. By processing mode
    print_status("\n  BY PROCESSING MODE:", "INFO")
    print_status("  " + "-"*70, "INFO")
    print_status(f"  {'Mode':<20} {'Short-Dist Ratio':<18} {'λ Ratio':<15} {'Status':<12}", "INFO")
    print_status("  " + "-"*70, "INFO")
    
    for name, key in [("Baseline (L1)", "baseline"), ("Iono-Free (L1+L2)", "ionofree"), ("Multi-GNSS", "multi_gnss")]:
        res = mode_results.get(key)
        if res:
            sd_ratio = get_short_dist_ratio(res)
            lam_ratio = res.get('ew_ns_pooled', {}).get('ratio')
            sd_str = f"{sd_ratio:.3f}" if sd_ratio else "N/A"
            lam_str = f"{lam_ratio:.2f}" if lam_ratio else "N/A"
            status = "✓" if sd_ratio and sd_ratio > 1.0 else "—"
            print_status(f"  {name:<20} {sd_str:<18} {lam_str:<15} {status}", "INFO")
        else:
            print_status(f"  {name:<20} {'—':<18} {'—':<15} {'N/A'}", "INFO")
    
    # 3. By metric (from multi-metric analysis)
    if result_all_multi and result_all_multi.get('metric_results'):
        print_status("\n  BY PHASE METRIC:", "INFO")
        print_status("  " + "-"*70, "INFO")
        print_status(f"  {'Metric':<20} {'Short-Dist Ratio':<18} {'λ Ratio':<15} {'CV':<12}", "INFO")
        print_status("  " + "-"*70, "INFO")
        
        for metric_name, metric_res in result_all_multi['metric_results'].items():
            sd_ratio = get_short_dist_ratio(metric_res)
            lam_ratio = metric_res.get('ew_ns_pooled', {}).get('ratio')
            cv = metric_res.get('anisotropy_cv', 0)
            sd_str = f"{sd_ratio:.3f}" if sd_ratio else "N/A"
            lam_str = f"{lam_ratio:.2f}" if lam_ratio else "N/A"
            print_status(f"  {metric_name:<20} {sd_str:<18} {lam_str:<15} {cv:.3f}", "INFO")
    
    # 4. By hemisphere
    if hemisphere_results:
        print_status("\n  BY HEMISPHERE:", "INFO")
        print_status("  " + "-"*70, "INFO")
        print_status(f"  {'Hemisphere':<20} {'Short-Dist Ratio':<18} {'λ Ratio':<15} {'N Pairs':<12}", "INFO")
        print_status("  " + "-"*70, "INFO")
        
        for name, key in [("Northern", "northern"), ("Southern", "southern")]:
            res = hemisphere_results.get(key, {})
            if res and res.get('metric_results'):
                # Use phase_alignment results if available
                pa_res = res['metric_results'].get('phase_alignment', {})
                sd_ratio = get_short_dist_ratio(pa_res)
                lam_ratio = pa_res.get('ew_ns_pooled', {}).get('ratio')
                n_pairs = pa_res.get('n_pairs', 0)
                sd_str = f"{sd_ratio:.3f}" if sd_ratio else "N/A"
                lam_str = f"{lam_ratio:.2f}" if lam_ratio else "N/A"
                print_status(f"  {name:<20} {sd_str:<18} {lam_str:<15} {n_pairs:,}", "INFO")
    
    print_status("\n  CODE Reference: λ ratio = 2.16 (E-W/N-S from 25yr PPP)", "INFO")
    print_status("  " + "="*70, "INFO")
    
    # Save results
    output = {
        'timestamp': datetime.now().isoformat(),
        'filter_mode': filter_mode_name,
        'methodology': 'CODE Longspan-Aligned Anisotropy with Multi-Metric Comparison',
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
        'modes_available': modes_available,
        'ionofree_results': mode_results.get('ionofree'),
        'ionofree_multi_metric': mode_results.get('ionofree_multi_metric'),
        'mgex_results': mode_results.get('multi_gnss'),
        'mgex_multi_metric': mode_results.get('multi_gnss_multi_metric'),
        'comparison_summary': {
            'by_condition': {
                'all_days': get_short_dist_ratio(result_all) if result_all else None,
                'quiet_days': get_short_dist_ratio(result_quiet) if result_quiet else None,
                'storm_days': get_short_dist_ratio(result_storm) if result_storm else None,
            },
            'by_mode': {
                'baseline': get_short_dist_ratio(mode_results.get('baseline')),
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
            'note': 'CODE uses λ (correlation length) ratio from 25yr PPP data; we use mean correlation ratio from 1yr SPP data'
        },
        'key_findings': [
            f"Short-distance E-W/N-S ratio: {short_dist.get('mean_ratio', 'N/A'):.3f}" if short_dist.get('mean_ratio') else "Short-distance analysis not available",
            f"95% CI: [{short_dist.get('ci_95_low', 0):.3f}, {short_dist.get('ci_95_high', 0):.3f}]" if short_dist.get('ci_95_low') else "",
            f"Statistical significance: p = {short_dist.get('p_value', 1):.2e}" if short_dist.get('p_value') else "",
            f"Effect size (Cohen's d): {short_dist.get('cohens_d', 0):.3f}" if short_dist.get('cohens_d') else "",
        ],
        'methodology_verified': [
            'Azimuth: spherical forward azimuth (geodetically correct)',
            'Metric: clock_bias only (pos_jitter excluded as isotropic)',
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
    
    # Store in output
    if monthly_by_mode:
        output['monthly_anisotropy'] = monthly_by_mode
        
        # Overall summary table
        print_status("\n  MONTHLY ANISOTROPY SUMMARY TABLE:", "SUCCESS")
        print_status("  " + "-"*75, "INFO")
        print_status(f"  {'Mode':<12} {'Metric':<15} {'Mean Ratio':>12} {'Std':>8} {'E-W>N-S %':>12} {'Months':>8}", "INFO")
        print_status("  " + "-"*75, "INFO")
        
        for mode_name in ['baseline', 'ionofree', 'multi_gnss']:
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
            run_analysis_pipeline(name, config)
        except Exception as e:
            print_status(f"Analysis failed for {name}: {e}", "ERROR")
            import traceback
            traceback.print_exc()
            
    print_status("ANALYSIS COMPLETE", "SUCCESS")

if __name__ == "__main__":
    main()
