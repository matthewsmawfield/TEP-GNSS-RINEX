#!/usr/bin/env python3
"""
TEP-GNSS-RINEX Analysis - STEP 2.0: Raw SPP Correlation Analysis
================================================================

Performs phase-coherent correlation analysis on Raw Single Point Positioning (SPP)
solutions derived from RINEX observation files using RTKLIB.

This is the "Fast Analysis" method that processes raw GNSS data without requiring
precise orbit/clock products from analysis centers. It serves as an independent
validation of TEP signatures detected in PPP-processed data.

Methodology:
    1. Process RINEX files with RTKLIB (rnx2rtkp) to get SPP solutions
    2. Extract position jitter and clock bias time series
    3. Apply band-pass filter (33 min - 12 hours) to isolate TEP-relevant frequencies
    4. Compute magnitude-weighted phase coherence for all station pairs
    5. Fit exponential decay model: C(r) = A*exp(-r/λ) + C₀
    6. Assess TEP consistency (λ in range 1000-3000 km, R² > 0.5)

Requirements: Step 1.0 complete (RINEX data downloaded)
Inputs:
    - data/rinex/*.crx, *.??o (RINEX observation files)
    - data/nav/*.??n (Broadcast navigation files)
Outputs:
    - results/outputs/step_2_0_raw_spp_analysis.json
    - results/figures/step_2_0_*.png

Author: Matthew Lukin Smawfield
Date: 9 December 2025
License: CC-BY-4.0
"""

import sys
from pathlib import Path

# Setup path for imports
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
import json
from datetime import datetime
from scipy.optimize import curve_fit
from scipy.signal import csd, welch, detrend
import matplotlib.pyplot as plt
import multiprocessing
from collections import defaultdict

from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.data_alignment import load_aligned_data, compute_aligned_coherence

# ============================================================================
# CONFIGURATION
# ============================================================================
DATA_DIR = PROJECT_ROOT / "data" / "rinex"
NAV_DIR = PROJECT_ROOT / "data" / "nav"
RESULTS_DIR = PROJECT_ROOT / "results" / "outputs"
FIGURES_DIR = PROJECT_ROOT / "results" / "figures"
LOGS_DIR = PROJECT_ROOT / "logs"
RTKLIB = PROJECT_ROOT / "RTKLIB" / "app" / "consapp" / "rnx2rtkp" / "gcc" / "rnx2rtkp"

# ============================================================================
# BATCH PROCESSING HELPERS
# ============================================================================
def get_day_files_subset(processed_dir, stations=None):
    """Group processed files by (Year, DOY)."""
    day_files = defaultdict(list)
    
    # Scan directory directly to avoid loading all files
    # Filename format: STAT_YYYYDOY.npz
    for f in processed_dir.glob("*.npz"):
        try:
            name_parts = f.stem.split('_')
            if len(name_parts) < 2: continue
            
            sta = name_parts[0]
            if stations and sta not in stations:
                continue
                
            # Parse YYYYDOY
            ts_str = name_parts[1]
            if len(ts_str) != 7: continue
            
            year = int(ts_str[:4])
            doy = int(ts_str[4:])
            
            day_files[(year, doy)].append(f)
        except:
            continue
            
    return day_files

# Ensure directories exist
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
LOGS_DIR.mkdir(parents=True, exist_ok=True)

# Initialize logger
logger = TEPLogger("step_2_0", log_file_path=LOGS_DIR / "step_2_0_raw_spp_analysis.log")
set_step_logger(logger)

# ============================================================================
# ANALYSIS PARAMETERS - IDENTICAL TO CODE-LONGSPAN (TEP-GNSS)
# ============================================================================
# Frequency Band: 10 µHz to 500 µHz (periods: 28 hours to 33 minutes)
# MUST MATCH CODE-LONGSPAN for valid cross-validation
F1_HZ = 1e-5   # 10 µHz (28 hour period) - lower bound
F2_HZ = 5e-4   # 500 µHz (33 min period) - upper bound
SAMPLING_PERIOD_SEC = 300.0   # 5-minute epoch interval (optimized processing)
FS_HZ = 1.0 / SAMPLING_PERIOD_SEC  # Sampling frequency (1/300 Hz)

# Distance Binning: Log-spaced, 50 km to 13,000 km - IDENTICAL TO CODE-LONGSPAN
MIN_DISTANCE_KM = 50
MAX_DISTANCE_KM = 13000
N_BINS = 40
MIN_BIN_COUNT = 10  # Minimum pairs per bin for fitting

# ============================================================================
# PROCESSING MODES - For comparative analysis across different data products
# ============================================================================
# SCIENTIFIC VALUE OF MULTI-MODE COMPARISON:
# If TEP is real: ALL modes should show similar λ (correlation length) and R²
# If TEP is artifact: Different error sources would produce different signatures
#
# Constellation comparison is particularly powerful because:
# - GPS: Rb/Cs clocks, MEO ~20,200 km
# - GLONASS: Cs clocks, MEO ~19,100 km, FDMA signals
# - Galileo: H-maser clocks (most precise!), MEO ~23,222 km
# - BeiDou: Rb clocks, MEO/GEO/IGSO mix

PROCESSING_MODES = {
    # === GPS BASELINE (reference) ===
    'baseline': {
        'clock_key': 'clock_bias_ns',
        'pos_key': 'pos_jitter',
        'drift_key': 'clock_drift',
        'description': 'GPS SPP with Broadcast Ephemeris',
    },
    
    # === ENHANCED GPS PROCESSING ===
    'precise': {
        'clock_key': 'precise_clock_bias_ns',
        'pos_key': 'precise_pos_jitter',
        'drift_key': 'precise_clock_drift',
        'description': 'GPS SPP with Precise Orbits',
    },
    'ionofree': {
        'clock_key': 'ionofree_clock_bias_ns',
        'pos_key': 'ionofree_pos_jitter',
        'drift_key': 'ionofree_clock_drift',
        'description': 'GPS Dual-Freq Iono-Free',
    },
    
    # === CONSTELLATION-SPECIFIC (cross-validation) ===
    'glonass': {
        'clock_key': 'glo_clock_bias_ns',
        'pos_key': 'glo_pos_jitter',
        'drift_key': 'glo_clock_drift',
        'description': 'GLONASS-only SPP (Cs clocks)',
    },
    'galileo': {
        'clock_key': 'gal_clock_bias_ns',
        'pos_key': 'gal_pos_jitter',
        'drift_key': 'gal_clock_drift',
        'description': 'Galileo-only SPP (H-maser clocks)',
    },
    'beidou': {
        'clock_key': 'bds_clock_bias_ns',
        'pos_key': 'bds_pos_jitter',
        'drift_key': 'bds_clock_drift',
        'description': 'BeiDou-only SPP (MEO/GEO/IGSO)',
    },
    
    # === COMBINED MULTI-GNSS ===
    'multi_gnss': {
        'clock_key': 'mgex_clock_bias_ns',
        'pos_key': 'mgex_pos_jitter',
        'drift_key': 'mgex_clock_drift',
        'description': 'All Constellations Combined',
    },
}

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def extract_coords_from_header(crx_file):
    """Fast coordinate extraction from RINEX header."""
    try:
        result = subprocess.run(['head', '-100', str(crx_file)], 
                               capture_output=True, text=True, timeout=5)
        for line in result.stdout.split('\n'):
            if 'APPROX POSITION XYZ' in line:
                parts = line.split()
                return float(parts[0]), float(parts[1]), float(parts[2])
    except Exception:
        pass
    return None


def ecef_to_lla(x, y, z):
    """Convert ECEF coordinates to Latitude/Longitude."""
    a = 6378137.0
    f = 1/298.257223563
    e2 = f * (2 - f)
    lon = np.arctan2(y, x)
    p = np.sqrt(x**2 + y**2)
    lat = np.arctan2(z, p * (1 - e2))
    for _ in range(5):
        N = a / np.sqrt(1 - e2 * np.sin(lat)**2)
        lat = np.arctan2(z + e2 * N * np.sin(lat), p)
    return np.degrees(lat), np.degrees(lon)


def haversine_km(lat1, lon1, lat2, lon2):
    """Calculate great-circle distance between two points."""
    R = 6371.0
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    return R * 2 * np.arcsin(np.sqrt(a))


def calculate_azimuth(lat1, lon1, lat2, lon2):
    """Compute azimuth from point 1 to point 2."""
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlon = lon2 - lon1
    y = np.sin(dlon) * np.cos(lat2)
    x = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    azimuth = np.degrees(np.arctan2(y, x))
    return (azimuth + 360) % 360


def apply_bandpass_filter(data, sampling_period_sec, low_period_sec, high_period_sec):
    """Apply Butterworth band-pass filter."""
    nyquist = 0.5 * (1 / sampling_period_sec)
    
    # Convert periods to frequencies
    low_freq = 1 / high_period_sec  # Longer period = Lower freq
    high_freq = 1 / low_period_sec  # Shorter period = Higher freq
    
    low = low_freq / nyquist
    high = high_freq / nyquist
    
    # 4th order Butterworth Bandpass
    b, a = signal.butter(4, [low, high], btype='band', analog=False)
    
    # Zero-phase filtering
    return signal.filtfilt(b, a, data)


def compute_cross_power_plateau(series1, series2, fs=FS_HZ, f1=F1_HZ, f2=F2_HZ):
    """
    Compute phase-coherent cross-power plateau using CSD (Welch's method).
    
    IDENTICAL TO CODE-LONGSPAN (TEP-GNSS) METHODOLOGY
    ================================================
    This is the exact same algorithm used in the 25-year CODE analysis.
    Using CSD ensures we measure the same physical quantity for valid comparison.
    
    Parameters:
    -----------
    series1, series2 : array-like
        Time series to correlate (same length, same sampling rate)
    fs : float
        Sampling frequency in Hz (default: 1/300 Hz for 5-minute epochs)
    f1, f2 : float
        Frequency band limits in Hz (default: 10µHz to 500µHz)
    
    Returns:
    --------
    coherence : float
        Phase-coherent correlation strength (comparable to CODE-Longspan)
    phase : float
        Representative phase difference in radians
    """
    from scipy.signal import csd, welch
    
    n_points = len(series1)
    if n_points < 100:
        return np.nan, np.nan
    
    # STEP 1: Detrend to remove systematic drifts (same as CODE-Longspan)
    time_indices = np.arange(n_points)
    series1_detrended = series1 - np.polyval(np.polyfit(time_indices, series1, 1), time_indices)
    series2_detrended = series2 - np.polyval(np.polyfit(time_indices, series2, 1), time_indices)
    
    # STEP 2: Compute cross-spectral density AND auto-spectra using Welch's method
    # nperseg=1024 matches CODE-Longspan for consistent frequency resolution
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
    
    # STEP 3: Select TEP frequency band (10µHz - 500µHz)
    band_mask = (frequencies > 0) & (frequencies >= f1) & (frequencies <= f2)
    if not np.any(band_mask):
        return np.nan, np.nan
    
    # STEP 4: Compute NORMALIZED COHERENCE = |Pxy|² / (Pxx × Pyy)
    # This gives values between 0 and 1
    Pxy_band = Pxy[band_mask]
    Pxx_band = Pxx[band_mask]
    Pyy_band = Pyy[band_mask]
    
    # Avoid division by zero
    denom = Pxx_band * Pyy_band
    valid_mask = denom > 0
    if not np.any(valid_mask):
        return np.nan, np.nan
    
    # Coherence squared (magnitude squared coherence)
    coh_squared = (np.abs(Pxy_band[valid_mask])**2) / denom[valid_mask]
    
    # Phase from cross-spectrum
    phases = np.angle(Pxy_band[valid_mask])
    magnitudes = np.sqrt(coh_squared)  # Coherence (0-1)
    
    if len(magnitudes) == 0 or np.sum(magnitudes) == 0:
        return np.nan, np.nan
    
    # STEP 5: Magnitude-weighted circular phase averaging (same as CODE-Longspan)
    complex_phases = np.exp(1j * phases)
    weighted_complex = np.average(complex_phases, weights=magnitudes)
    weighted_phase = np.angle(weighted_complex)
    
    # Mean coherence in band (normalized, 0-1)
    mean_coherence = np.mean(magnitudes)
    
    return float(mean_coherence), float(weighted_phase)


def _compute_coherence_phase(v1, v2, fs, f1, f2):
    """
    Compute coherence and phase for a single pair of time series.
    Returns (mean_coh, weighted_phase, phase_alignment) or (nan, nan, nan) on failure.
    
    OPTIMIZED: Single function call, minimal overhead.
    """
    n_points = len(v1)
    if n_points < 64:
        return np.nan, np.nan, np.nan
    
    # Detrend
    v1_d = detrend(v1, type='linear')
    v2_d = detrend(v2, type='linear')
    
    # nperseg must be < n_points to get multiple segments for coherence averaging
    # With nperseg = n_points, coherence is always 1.0 (single segment, no averaging)
    # Use n_points//2 to ensure at least 2 segments with 50% overlap (default)
    nperseg = min(256, n_points // 2)
    if nperseg < 32:
        return np.nan, np.nan, np.nan  # Too few points for reliable estimate
    
    try:
        # Compute spectra
        f, Pxy = csd(v1_d, v2_d, fs=fs, nperseg=nperseg, detrend='constant')
        _, Pxx = welch(v1_d, fs=fs, nperseg=nperseg, detrend='constant')
        _, Pyy = welch(v2_d, fs=fs, nperseg=nperseg, detrend='constant')
        
        # Band mask
        mask = (f >= f1) & (f <= f2)
        if not np.any(mask):
            return np.nan, np.nan, np.nan
        
        Pxy_band = Pxy[mask]
        Pxx_band = Pxx[mask]
        Pyy_band = Pyy[mask]
        
        # Normalized MSC
        denom = Pxx_band * Pyy_band
        valid_denom = denom > 0
        if not np.any(valid_denom):
            return np.nan, np.nan, np.nan
        
        coh_squared = np.abs(Pxy_band[valid_denom])**2 / denom[valid_denom]
        mean_coh = np.mean(np.sqrt(coh_squared))
        
        # Phase alignment (CODE longspan method)
        raw_magnitudes = np.abs(Pxy_band)
        phases = np.angle(Pxy_band)
        
        if np.sum(raw_magnitudes) > 0:
            complex_phases = np.exp(1j * phases)
            weighted_complex = np.average(complex_phases, weights=raw_magnitudes)
            weighted_phase = np.angle(weighted_complex)
            phase_alignment = np.cos(weighted_phase)
        else:
            weighted_phase = np.nan
            phase_alignment = np.nan
        
        return float(mean_coh), float(weighted_phase), float(phase_alignment)
    except:
        return np.nan, np.nan, np.nan


def _process_day_batch_worker(args):
    """
    Process a batch of days for multiple metrics.
    OPTIMIZED: Geometry pre-computed in main process and passed in.
    """
    day_batch, batch_files, metrics_config, pair_geometry = args
    # pair_geometry = (sta_to_idx, pair_info, sta_latlon) - pre-computed in main process
    sta_to_idx, pair_info, sta_latlon = pair_geometry
    results = defaultdict(list)
    
    processed_count = 0
    for year, doy in day_batch:
        daily_files = batch_files.get((year, doy), [])
        if not daily_files:
            continue
        
        # Load all station data for this day
        daily_data = {}
        for f in daily_files:
            try:
                sta = f.stem.split('_')[0]
                if sta not in sta_to_idx:
                    continue
                
                with np.load(f) as d:
                    ts = d.get('timestamps')
                    if ts is None or len(ts) < MIN_BIN_COUNT:
                        continue
                    
                    # Filter valid timestamps
                    valid = (ts >= 0) & (ts < 86400)
                    if np.sum(valid) < MIN_BIN_COUNT:
                        continue
                    
                    st_data = {'ts': ts[valid], 'idx': sta_to_idx[sta]}
                    has_data = False
                    
                    for metric_name, key in metrics_config:
                        val = d.get(key)
                        if val is not None and len(val) == len(ts):
                            st_data[metric_name] = val[valid]
                            has_data = True
                    
                    if has_data:
                        daily_data[sta] = st_data
            except:
                continue
        
        if len(daily_data) < 2:
            continue
        
        # Process valid pairs only
        stations = list(daily_data.keys())
        for i_s, sta1 in enumerate(stations):
            d1 = daily_data[sta1]
            idx1 = d1['idx']
            
            for sta2 in stations[i_s + 1:]:
                d2 = daily_data[sta2]
                idx2 = d2['idx']
                
                # Get pre-computed pair info (ensure i < j)
                pair_key = (idx1, idx2) if idx1 < idx2 else (idx2, idx1)
                if pair_key not in pair_info:
                    continue
                
                dist, _, mid_lat = pair_info[pair_key]  # Ignore stored azimuth
                
                # CRITICAL FIX: Compute azimuth directly from sta1 to sta2.
                # On a sphere, reversing direction is NOT simply +180° due to great circle curvature.
                lat1, lon1 = sta_latlon[sta1]
                lat2, lon2 = sta_latlon[sta2]
                az = calculate_azimuth(lat1, lon1, lat2, lon2)
                
                # Align timestamps once per pair
                common_ts, loc1, loc2 = np.intersect1d(d1['ts'], d2['ts'], return_indices=True)
                if len(common_ts) < MIN_BIN_COUNT:
                    continue
                
                # Process each metric
                for metric_name, _ in metrics_config:
                    if metric_name not in d1 or metric_name not in d2:
                        continue
                    
                    v1 = d1[metric_name][loc1]
                    v2 = d2[metric_name][loc2]
                    
                    mean_coh, weighted_phase, phase_alignment = _compute_coherence_phase(
                        v1, v2, FS_HZ, F1_HZ, F2_HZ
                    )
                    
                    if np.isfinite(mean_coh):
                        results[metric_name].append((
                            int(year), int(doy), dist, mean_coh, weighted_phase, 
                            phase_alignment, az, mid_lat, sta1, sta2
                        ))
    
    return results

def compute_pairs_coherence_stream(day_files, station_coords, mode_metrics, checkpoint_path=None):
    """Yields batch results for streaming processing with checkpoint support."""
    day_keys = sorted(list(day_files.keys()))
    n_days = len(day_keys)
    
    # PRE-COMPUTE GEOMETRY ONCE (not in each worker!)
    # This saves 144k distance calculations × 80 workers = 11.5M redundant calculations
    print_status("  Pre-computing pairwise geometry...", "INFO")
    sta_list = list(station_coords.keys())
    sta_to_idx = {s: i for i, s in enumerate(sta_list)}
    lats = np.array([station_coords[s]['lat'] for s in sta_list])
    lons = np.array([station_coords[s]['lon'] for s in sta_list])
    
    pair_info = {}  # (i, j) -> (dist, az, mid_lat)
    n_sta = len(sta_list)
    for i in range(n_sta):
        for j in range(i + 1, n_sta):
            dist = haversine_km(lats[i], lons[i], lats[j], lons[j])
            if MIN_DISTANCE_KM <= dist <= MAX_DISTANCE_KM:
                az = calculate_azimuth(lats[i], lons[i], lats[j], lons[j])
                mid_lat = (lats[i] + lats[j]) / 2
                pair_info[(i, j)] = (dist, az, mid_lat)
    
    # Include station lat/lon for correct azimuth computation in workers
    sta_latlon = {s: (station_coords[s]['lat'], station_coords[s]['lon']) for s in sta_list}
    pair_geometry = (sta_to_idx, pair_info, sta_latlon)
    print_status(f"  Pre-computed {len(pair_info)} valid station pairs", "SUCCESS")
    
    # Checkpoint handling
    completed_days = set()
    if checkpoint_path and checkpoint_path.exists():
        try:
            with open(checkpoint_path) as f:
                checkpoint_data = json.load(f)
                completed_days = set(tuple(d) for d in checkpoint_data.get('completed_days', []))
                print_status(f"  Resuming from checkpoint: {len(completed_days)}/{n_days} days already done", "SUCCESS")
        except Exception as e:
            print_status(f"  Checkpoint load failed: {e}, starting fresh", "WARNING")
    
    # Filter out completed days
    remaining_days = [d for d in day_keys if d not in completed_days]
    print_status(f"Processing {len(remaining_days)} remaining days (of {n_days} total)...", "INFO")
    
    if not remaining_days:
        print_status("  All days already processed!", "SUCCESS")
        return
    
    # OPTIMIZED FOR HIGH-CPU MACHINES (n1-highcpu-96: 96 vCPU, 86.4GB RAM)
    # Geometry is pre-computed in main process, so batch_size only affects I/O grouping
    # Smaller batches = faster feedback, same total work
    cpu_count = multiprocessing.cpu_count()
    if cpu_count >= 48:
        # High-CPU machine: maximize parallelism
        n_workers = 80  # Use most of 96 vCPUs
        batch_size = 1  # 1 day per batch for fastest feedback
    elif cpu_count >= 16:
        # Medium machine
        n_workers = min(24, cpu_count - 4)
        batch_size = 4
    else:
        # Small machine (original config)
        n_workers = min(12, max(1, cpu_count))
        batch_size = 4
    batches = [remaining_days[i:i + batch_size] for i in range(0, len(remaining_days), batch_size)]
    
    print_status(f"  Workers: {n_workers}, Batches: {len(batches)} (batch_size={batch_size})", "INFO")
    
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        # CRITICAL OPTIMIZATION: 
        # 1. Only pass files for each batch (not entire 400k file dict)
        # 2. Pass pre-computed geometry (not raw coords that need recomputation)
        futures = {}
        for i, batch_days in enumerate(batches):
            batch_files = {d: day_files[d] for d in batch_days if d in day_files}
            future = executor.submit(_process_day_batch_worker, (batch_days, batch_files, mode_metrics, pair_geometry))
            futures[future] = (i, batch_days)
        
        completed = len(completed_days)
        total_pairs = 0
        checkpoint_interval = 10  # Save checkpoint every 10 batches
        batches_since_checkpoint = 0
        
        for future in as_completed(futures):
            batch_idx, batch_days = futures[future]
            try:
                res = future.result()
                completed += len(batch_days)
                batch_pairs = sum(len(v) for v in res.values())
                total_pairs += batch_pairs
                
                # Track completed days
                for d in batch_days:
                    completed_days.add(d)
                
                batches_since_checkpoint += 1
                
                # Save checkpoint periodically
                if checkpoint_path and batches_since_checkpoint >= checkpoint_interval:
                    try:
                        checkpoint_data = {'completed_days': [list(d) for d in completed_days], 'total_pairs': total_pairs}
                        with open(checkpoint_path, 'w') as f:
                            json.dump(checkpoint_data, f)
                        batches_since_checkpoint = 0
                    except Exception as e:
                        print_status(f"  Checkpoint save failed: {e}", "WARNING")
                
                print_status(f"  Day {completed}/{n_days} complete: +{batch_pairs} pairs (total: {total_pairs})", "INFO")
                yield res
                
            except Exception as e:
                print_status(f"Batch failed: {e}", "ERROR")
                import traceback
                traceback.print_exc()
                yield {}
        
        # Final checkpoint
        if checkpoint_path:
            try:
                checkpoint_data = {'completed_days': [list(d) for d in completed_days], 'total_pairs': total_pairs, 'status': 'complete'}
                with open(checkpoint_path, 'w') as f:
                    json.dump(checkpoint_data, f)
                print_status(f"  Final checkpoint saved: {len(completed_days)} days", "SUCCESS")
            except Exception as e:
                print_status(f"  Final checkpoint save failed: {e}", "WARNING")

def exp_decay(r, A, lam, C0):
    return A * np.exp(-r/lam) + C0

def fit_exponential_model(distances, coherences):
    """Fit exponential decay model to pair data."""
    if len(distances) < 100: return None
    
    # Binning
    bin_edges = np.logspace(np.log10(MIN_DISTANCE_KM), np.log10(MAX_DISTANCE_KM), N_BINS + 1)
    bin_centers, bin_means, bin_counts = [], [], []
    
    for i in range(N_BINS):
        mask = (distances >= bin_edges[i]) & (distances < bin_edges[i+1])
        if np.sum(mask) >= MIN_BIN_COUNT:
            bin_centers.append((bin_edges[i] + bin_edges[i+1]) / 2)
            bin_means.append(np.mean(coherences[mask]))
            bin_counts.append(np.sum(mask))
            
    if len(bin_centers) < 5: return None
    
    x = np.array(bin_centers)
    y = np.array(bin_means)
    w = np.sqrt(np.array(bin_counts))
    
    try:
        popt, pcov = curve_fit(
            exp_decay, x, y,
            p0=[0.5, 2000, 0],
            sigma=1.0/w,
            bounds=([0, 100, -1], [2, 20000, 1]),
            maxfev=5000
        )
        
        # R2
        predicted = exp_decay(x, *popt)
        ss_res = np.sum((y - predicted)**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        return {
            'amplitude': popt[0],
            'correlation_length_km': popt[1],
            'offset': popt[2],
            'r_squared': r2,
            'success': True
        }
    except:
        return None

def process_station_rtklib(args):
    """Process one station with ALL its RINEX files using RTKLIB."""
    sta, obs_files, nav_files, station_pos = args
    
    out_file = RESULTS_DIR / f"rtklib_{sta}.pos"
    stat_file = RESULTS_DIR / f"rtklib_{sta}.pos.stat"
    
    # Build RTKLIB command with all observation and nav files
    cmd = [str(RTKLIB), '-p', '0', '-sys', 'G', '-y', '1', '-o', str(out_file)]
    cmd.extend([str(f) for f in obs_files])
    cmd.extend([str(f) for f in nav_files])
    
    try:
        result = subprocess.run(cmd, capture_output=True, timeout=300)
        
        # Parse output - extract position time series AND clock
        pos_data = []
        
        if out_file.exists():
            with open(out_file, 'r') as f:
                for line in f:
                    if line.startswith('%') or not line.strip():
                        continue
                    parts = line.split()
                    if len(parts) >= 5:
                        try:
                            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                            pos_data.append([x, y, z])
                        except Exception:
                            continue
        
        clock_data = []
        if stat_file.exists():
            with open(stat_file, 'r') as f:
                for line in f:
                    if line.startswith('$CLK'):
                        parts = line.split(',')
                        if len(parts) >= 6:
                            try:
                                clk_bias_m = float(parts[5])
                                clock_data.append(clk_bias_m)
                            except Exception:
                                continue
            stat_file.unlink()
            
        if out_file.exists():
            out_file.unlink()
        
        # Require minimum data points
        min_pts = 1000
        if len(pos_data) < min_pts or len(clock_data) < min_pts:
            return sta, None
            
        # Align lengths
        L = min(len(pos_data), len(clock_data))
        pos_data = np.array(pos_data[:L])
        
        # --- 1. Position Jitter Analysis (Proxy) ---
        mean_pos = pos_data.mean(axis=0)
        dr = np.sqrt(np.sum((pos_data - mean_pos)**2, axis=1))
        dr_filtered = apply_bandpass_filter(dr, SAMPLING_PERIOD_SEC, LOW_PERIOD_SEC, HIGH_PERIOD_SEC)
        
        # --- 2. True Clock Bias Analysis (Time) ---
        clock_ns = np.array(clock_data[:L]) * 3.33564095  # Convert meters to nanoseconds
        bias_filtered = apply_bandpass_filter(clock_ns, SAMPLING_PERIOD_SEC, LOW_PERIOD_SEC, HIGH_PERIOD_SEC)
        
        # --- 3. Clock Drift Analysis (Freq) ---
        clock_drift = np.diff(clock_ns) / SAMPLING_PERIOD_SEC
        clock_drift[np.abs(clock_drift) > 10] = 0  # Outliers to zero
        drift_filtered = apply_bandpass_filter(clock_drift, SAMPLING_PERIOD_SEC, LOW_PERIOD_SEC, HIGH_PERIOD_SEC)
        
        return sta, {'jitter': dr_filtered, 'bias': bias_filtered, 'drift': drift_filtered}
        
    except Exception:
        return sta, None




# ============================================================================
    
    return fit_results, metrics_data


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def run_analysis_pipeline(filter_mode_name, station_filter_config):
    """Run the full analysis pipeline for a specific filter mode."""
    print_status("", "INFO")
    print_status("=" * 80, "INFO")
    print_status(f"TEP-GNSS-RINEX Analysis - STEP 2.0: Multi-Metric SPP Correlation Analysis ({filter_mode_name})", "TITLE")
    print_status("=" * 80, "INFO")
    
    # Output filename suffix
    filter_suffix = filter_mode_name.lower().replace(" ", "_").replace("-", "_")
    
    # ========================================================================
    # STEP 1: Load Processed Data
    # ========================================================================
    print_status("[1/3] Loading processed .npz files...", "PROCESS")
    
    PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
    
    # Configure filtering
    dynamic_threshold = None
    good_stations_set = None
    
    if station_filter_config == 'none':
        print_status(f"  Using ALL stations (no filter)", "WARNING")
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

    npz_files = sorted(PROCESSED_DIR.glob("*.npz"))
    if good_stations_set:
        npz_files = [f for f in npz_files if f.name.split('_')[0] in good_stations_set]
    
    # Apply dynamic filtering if enabled
    if dynamic_threshold is not None:
        good_files = []
        rejected_jump = 0
        rejected_range = 0
        for f in npz_files:
            try:
                with np.load(f) as d:
                    clk = d.get('clock_bias_ns')
                    if clk is not None and len(clk) > 10:
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
        
    print_status(f"  Found {len(npz_files)} processed files matching filter", "INFO")
    
    # Load station coordinates
    coords_file = PROCESSED_DIR / "station_coordinates.json"
    if not coords_file.exists():
        print_status("Station coordinates file not found", "ERROR")
        return
    
    with open(coords_file) as f:
        station_coords_ecef = json.load(f)
    
    # Convert ECEF to lat/lon
    station_coords = {}
    for sta, ecef in station_coords_ecef.items():
        if len(ecef) == 3:
            lat, lon = ecef_to_lla(ecef[0], ecef[1], ecef[2])
            station_coords[sta] = {'lat': lat, 'lon': lon}
    
    print_status(f"  Loaded coordinates for {len(station_coords)} stations", "INFO")
    
    # ========================================================================
    # STEP 2: Batch Correlation Analysis - ALL MODES
    # ========================================================================
    print_status("[2/3] Computing CSD correlations using Day-by-Day Batch Processing...", "PROCESS")
    
    # Group files by Day (Year, DOY)
    valid_files_set = {f.name for f in npz_files}
    day_files_all = get_day_files_subset(PROCESSED_DIR)
    
    # Filter to only include files selected by smart filter
    day_files = defaultdict(list)
    for key, files in day_files_all.items():
        valid_subset = [f for f in files if f.name in valid_files_set]
        if valid_subset:
            day_files[key] = valid_subset
            
    print_status(f"  Grouped {len(npz_files)} files into {len(day_files)} day batches", "INFO")

    # Define Processing Modes - all 12 metrics (added precise)
    modes_config = {
        'baseline': [('clock_bias', 'clock_bias_ns'), ('pos_jitter', 'pos_jitter'), ('clock_drift', 'clock_drift')],
        'ionofree': [('ionofree_clock_bias', 'ionofree_clock_bias_ns'), ('ionofree_pos_jitter', 'ionofree_pos_jitter'), ('ionofree_clock_drift', 'ionofree_clock_drift')],
        'multi_gnss': [('mgex_clock_bias', 'mgex_clock_bias_ns'), ('mgex_pos_jitter', 'mgex_pos_jitter'), ('mgex_clock_drift', 'mgex_clock_drift')],
        'precise': [('precise_clock_bias', 'precise_clock_bias_ns'), ('precise_pos_jitter', 'precise_pos_jitter'), ('precise_clock_drift', 'precise_clock_drift')]
    }
    
    # OPTIMIZATION: Flatten metrics to run in a single batch pass (reads each file only once)
    all_metrics_config = []
    metric_to_mode = {}
    
    for mode, metrics in modes_config.items():
        for metric_name, key in metrics:
            all_metrics_config.append((metric_name, key))
            metric_to_mode[metric_name] = mode
            
    print_status(f"Running Single-Pass Optimized Batch Processing for {len(all_metrics_config)} metrics...", "PROCESS")
    
    # Check for existing checkpoint using unique filename for this filter mode
    checkpoint_file = RESULTS_DIR / f"step_2_0_checkpoint_{filter_suffix}.json"
    resuming = checkpoint_file.exists()
    
    # Initialize CSV outputs (append mode if resuming)
    csv_handles = {}
    for mode in modes_config.keys():
        # Use unique CSV filenames for this filter mode
        csv_path = RESULTS_DIR / f"step_2_0_pairs_{mode}_{filter_suffix}.csv"
        if resuming and csv_path.exists():
            f = open(csv_path, 'a')  # Append mode for resume
        else:
            f = open(csv_path, 'w')
            f.write("year,doy,metric,distance_km,coherence,phase_rad,phase_alignment,azimuth,mid_lat,station1,station2\n")
        csv_handles[mode] = f
        
    all_mode_results = defaultdict(dict)
    pair_counts = defaultdict(lambda: defaultdict(int))  # Track counts per metric
    
    try:
        # Stream processing with checkpoint support
        for batch_results in compute_pairs_coherence_stream(day_files, station_coords, all_metrics_config, checkpoint_path=checkpoint_file):
            for metric_name, pairs in batch_results.items():
                mode_name = metric_to_mode.get(metric_name)
                if not mode_name: continue
                
                # Batch write to CSV (much faster than per-line)
                f = csv_handles[mode_name]
                lines = []
                for p in pairs:
                    # p is tuple: (year, doy, dist, coh, phase, phase_align, az, mid_lat, sta1, sta2)
                    lines.append(f"{p[0]},{p[1]},{metric_name},{p[2]:.4f},{p[3]:.6f},{p[4]:.6f},{p[5]:.6f},{p[6]:.2f},{p[7]:.4f},{p[8]},{p[9]}\n")
                
                f.write(''.join(lines))  # Batch write
                pair_counts[mode_name][metric_name] += len(pairs)
                    
    except Exception as e:
        print_status(f"Processing failed: {e}", "ERROR")
    finally:
        for f in csv_handles.values():
            f.close()
    
    # Report pair counts
    for mode_name, metrics in pair_counts.items():
        for metric_name, count in metrics.items():
            print_status(f"    [{mode_name.upper()}] {metric_name}: {count:,} pairs written to CSV", "INFO")
            
    # Read ALL data from CSV for fitting (memory efficient chunked reading)
    # Now supports BOTH coherence metrics: MSC (normalized) and phase_alignment (cos(phase))
    print_status("Loading data from CSV for fitting (ALL pairs)...", "PROCESS")
    metrics_data = defaultdict(lambda: {'distances': [], 'coherences': [], 'phase_alignments': []})
    
    for mode_name in modes_config.keys():
        csv_path = RESULTS_DIR / f"step_2_0_pairs_{mode_name}_{filter_suffix}.csv"
        if not csv_path.exists():
            continue
            
        # Read CSV in chunks to avoid memory issues
        import csv as csv_module
        with open(csv_path, 'r') as f:
            reader = csv_module.DictReader(f)
            for row in reader:
                metric_name = row['metric']
                metrics_data[metric_name]['distances'].append(float(row['distance_km']))
                metrics_data[metric_name]['coherences'].append(float(row['coherence']))
                # Handle both old (no phase) and new (with phase) CSV formats
                if 'phase_alignment' in row and row['phase_alignment']:
                    try:
                        pa = float(row['phase_alignment'])
                        metrics_data[metric_name]['phase_alignments'].append(pa if np.isfinite(pa) else np.nan)
                    except:
                        metrics_data[metric_name]['phase_alignments'].append(np.nan)
                else:
                    metrics_data[metric_name]['phase_alignments'].append(np.nan)
    
    # Perform Fit on ALL data - BOTH coherence metrics
    for metric_name, data in metrics_data.items():
        dists = np.array(data['distances'])
        cohs = np.array(data['coherences'])
        phase_aligns = np.array(data['phase_alignments'])
        
        if len(dists) < 100:
            continue
        
        mode_name = metric_to_mode.get(metric_name, 'baseline')
        
        # Fit 1: Normalized MSC (coherence) - range [0, 1]
        res = fit_exponential_model(dists, cohs)
        if res:
            all_mode_results[mode_name][metric_name] = res
            res['n_pairs'] = len(dists)
            res['metric_type'] = 'normalized_msc'
            lam = res['correlation_length_km']
            r2 = res['r_squared']
            tep = "YES" if 500 < lam < 5000 and r2 > 0.5 else "NO"
            print_status(f"    [{mode_name.upper()}] {metric_name} (MSC): λ={lam:.0f}km, R²={r2:.3f}, N={len(dists):,}, TEP={tep}", 
                       "SUCCESS" if tep == "YES" else "INFO")
        
        # Fit 2: Phase Alignment cos(phase) - range [-1, 1] (like CODE longspan)
        valid_pa = np.isfinite(phase_aligns)
        if np.sum(valid_pa) > 100:
            res_pa = fit_exponential_model(dists[valid_pa], phase_aligns[valid_pa])
            if res_pa:
                pa_key = f"{metric_name}_phase_alignment"
                all_mode_results[mode_name][pa_key] = res_pa
                res_pa['n_pairs'] = int(np.sum(valid_pa))
                res_pa['metric_type'] = 'phase_alignment'
                lam = res_pa['correlation_length_km']
                r2 = res_pa['r_squared']
                tep = "YES" if 500 < lam < 5000 and r2 > 0.5 else "NO"
                print_status(f"    [{mode_name.upper()}] {metric_name} (Phase): λ={lam:.0f}km, R²={r2:.3f}, N={np.sum(valid_pa):,}, TEP={tep}", 
                           "SUCCESS" if tep == "YES" else "INFO")

    if not metrics_data:
        print_status("No valid metrics aggregated", "ERROR")
        return
        
    # Use baseline metrics for plotting
    fit_results = all_mode_results.get('baseline', {})
    
    # ========================================================================
    # STEP 3: Fit Exponential Decay Models (Publication Quality)
    # ========================================================================
    print_status("[3/3] Fitting exponential decay models...", "PROCESS")
    
    def exp_decay(r, A, lam, C0):
        return A * np.exp(-r/lam) + C0
    
    fit_results = {}
    
    # Metric display names
    metric_labels = {
        'pos_jitter': 'Position Jitter (Proxy)',
        'clock_bias': 'Clock Bias (Time)',
        'clock_drift': 'Clock Drift (Doppler)'
    }
    
    for metric_name, samples in metrics_data.items():
        if len(samples.get('distances', [])) < 100:
            print_status(f"  Skipping {metric_name}: too few pairs", "WARNING")
            continue
        
        # Use pre-collected arrays
        distances = np.array(samples['distances'])
        coherences = np.array(samples['coherences'])
        
        # LOG bins - captures exponential decay shape better at short distances
        bin_edges = np.logspace(np.log10(MIN_DISTANCE_KM), np.log10(MAX_DISTANCE_KM), N_BINS + 1)
        
        bin_centers = []
        bin_means = []
        bin_sems = []  # Standard Error of Mean
        bin_counts = []
        
        for i in range(N_BINS):
            mask = (bin_edges[i] <= distances) & (distances < bin_edges[i+1])
            bin_data = coherences[mask]
            
            if len(bin_data) >= MIN_BIN_COUNT:
                bin_centers.append((bin_edges[i] + bin_edges[i+1]) / 2)
                bin_means.append(np.mean(bin_data))
                bin_sems.append(np.std(bin_data) / np.sqrt(len(bin_data)))  # SEM
                bin_counts.append(len(bin_data))
        
        if len(bin_centers) < 5:
            print_status(f"  {metric_name}: too few bins", "WARNING")
            continue
        
        # Convert to numpy arrays
        bin_centers = np.array(bin_centers)
        bin_means = np.array(bin_means)
        bin_sems = np.array(bin_sems)
        bin_counts = np.array(bin_counts)
        
        try:
            # Weighted fit using inverse SEM
            weights = bin_counts  # More samples = more weight
            popt, pcov = curve_fit(
                exp_decay,
                bin_centers,
                bin_means,
                p0=[0.5, 2000, 0],
                sigma=1.0/np.sqrt(weights),
                bounds=([0, 100, -1], [2, 20000, 1]),
                maxfev=10000
            )
            
            A, lam, C0 = popt
            A_err, lam_err, C0_err = np.sqrt(np.diag(pcov))
            
            # Calculate R²
            predicted = exp_decay(bin_centers, A, lam, C0)
            ss_res = np.sum((bin_means - predicted)**2)
            ss_tot = np.sum((bin_means - np.mean(bin_means))**2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            
            fit_results[metric_name] = {
                'amplitude': float(A),
                'amplitude_err': float(A_err),
                'correlation_length_km': float(lam),
                'correlation_length_err_km': float(lam_err),
                'offset': float(C0),
                'offset_err': float(C0_err),
                'r_squared': float(r_squared),
                'success': True,
                'n_pairs': len(distances)
            }
            
            print_status(f"  {metric_name}:", "INFO")
            print_status(f"    λ: {lam:.0f} ± {lam_err:.0f} km, R²: {r_squared:.3f}", "INFO")
            
            # ================================================================
            # PUBLICATION-QUALITY PLOT (Linear scale, SEM errors)
            # ================================================================
            plt.style.use('seaborn-v0_8-whitegrid')
            fig, ax = plt.subplots(figsize=(10, 6))
            
            # 1. Individual pairs as faint scatter background
            if len(distances) > 10000:
                idx = np.random.choice(len(distances), 10000, replace=False)
                ax.scatter(distances[idx], coherences[idx], s=10, alpha=0.05,
                          c='#bdc3c7', edgecolors='none')
            else:
                ax.scatter(distances, coherences, s=10, alpha=0.05,
                          c='#bdc3c7', edgecolors='none')
            
            # 2. Fit line (smooth red) - log-spaced for log x-axis
            x_smooth = np.logspace(np.log10(MIN_DISTANCE_KM), np.log10(MAX_DISTANCE_KM), 500)
            y_smooth = exp_decay(x_smooth, A, lam, C0)
            ax.plot(x_smooth, y_smooth, color='#e74c3c', linewidth=2.5, zorder=5,
                   label=f'Fit: $\\lambda={lam:.0f}$ km ($R^2={r_squared:.2f}$)')
            
            # 3. Binned means with SEM error bars (only stable bins)
            plot_centers, plot_means, plot_errs = [], [], []
            for c, m, e, cnt in zip(bin_centers, bin_means, bin_sems, bin_counts):
                if cnt > 20:
                    plot_centers.append(c)
                    plot_means.append(m)
                    plot_errs.append(e)
            
            ax.errorbar(plot_centers, plot_means, yerr=plot_errs, fmt='o',
                       color='#2980b9', ecolor='#2980b9', elinewidth=1.5, capsize=2,
                       markersize=5, alpha=0.9, zorder=10, label='Binned Mean ± SEM')
            
            # Formatting
            ax.axhline(0, color='black', linewidth=0.8, linestyle='--', alpha=0.5)
            ax.set_xlabel('Distance (km)', fontsize=12, fontweight='bold')
            ax.set_ylabel(f'{metric_labels.get(metric_name, metric_name)} Coherence', 
                         fontsize=12, fontweight='bold')
            ax.set_title(f'TEP Analysis: {metric_labels.get(metric_name, metric_name)}',
                        fontsize=14, fontweight='bold', pad=15)
            
            # Stats box
            stats_text = (f"$\\lambda = {lam:.0f} \\pm {lam_err:.0f}$ km\n"
                         f"$A = {A:.3f}$\n"
                         f"$R^2 = {r_squared:.3f}$")
            props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='#bdc3c7')
            ax.text(0.97, 0.97, stats_text, transform=ax.transAxes, fontsize=11,
                   verticalalignment='top', horizontalalignment='right', bbox=props)
            
            ax.legend(loc='upper right', bbox_to_anchor=(0.97, 0.75), fontsize=10, frameon=True)
            ax.set_ylim(-0.1, 1.1)
            ax.set_xscale('log')
            ax.set_xlim(MIN_DISTANCE_KM, MAX_DISTANCE_KM)
            ax.grid(True, which='both', linestyle='-', alpha=0.3)
            
            # TEP detection badge
            tep_detected = 500 < lam < 5000 and r_squared > 0.5
            if tep_detected:
                ax.text(0.03, 0.97, 'TEP DETECTED', transform=ax.transAxes, fontsize=11,
                       fontweight='bold', verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.9))
            
            plt.tight_layout()
            fig_file = FIGURES_DIR / f'step_2_0_{metric_name}_{filter_suffix}.png'
            fig.savefig(fig_file, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print_status(f"    Figure: {fig_file}", "SUCCESS")
            
        except Exception as e:
            print_status(f"  {metric_name} fit failed: {e}", "WARNING")
            fit_results[metric_name] = {'success': False, 'error': str(e)}
    
    # ========================================================================
    # Save Results - ALL MODES
    # ========================================================================
    results = {
        'step': '2.0',
        'name': 'Multi-Metric SPP Correlation Analysis - All Processing Modes',
        'filter_mode': filter_mode_name,
        'completion_time': datetime.now().isoformat(),
        'status': 'completed',
        'n_stations': len(station_coords),
        'n_pairs': {metric: len(data.get('distances', [])) for metric, data in metrics_data.items()},
        'analysis_by_mode': all_mode_results,  # All 3 modes
        'analysis': fit_results,  # Baseline for backward compatibility
        'parameters': {
            'frequency_band_hz': [F1_HZ, F2_HZ],
            'frequency_band_period_sec': [1/F2_HZ, 1/F1_HZ],
            'distance_range_km': [MIN_DISTANCE_KM, MAX_DISTANCE_KM],
            'n_bins': N_BINS,
            'min_bin_count': MIN_BIN_COUNT
        },
        'modes_available': {mode: bool(all_mode_results.get(mode)) for mode in modes_config},
        'mode_station_counts': {mode: 'dynamic' for mode in modes_config}
    }
    
    output_file = RESULTS_DIR / f'step_2_0_raw_spp_analysis_{filter_suffix}.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print_status(f"Results saved: {output_file}", "SUCCESS")
    
    # ========================================================================
    # COMPARATIVE SUMMARY - ALL MODES
    # ========================================================================
    print_status("", "INFO")
    print_status("="*80, "INFO")
    print_status("COMPARATIVE ANALYSIS SUMMARY - ALL PROCESSING MODES", "INFO")
    print_status("="*80, "INFO")
    
    # Print table header
    print_status("", "INFO")
    print_status(f"  {'Mode':<12} {'Metric':<15} {'λ (km)':<12} {'R²':<8} {'TEP?':<6}", "INFO")
    print_status(f"  {'-'*12} {'-'*15} {'-'*12} {'-'*8} {'-'*6}", "INFO")
    
    tep_count = 0
    total_tests = 0
    
    for mode_name, mode_results in all_mode_results.items():
        for metric, result in mode_results.items():
            if result.get('success'):
                lam = result['correlation_length_km']
                r2 = result['r_squared']
                tep = "YES" if 500 < lam < 5000 and r2 > 0.5 else "NO"
                
                # Count TEP detections for clock metrics across all modes
                is_clock = 'clock_bias' in metric
                
                if is_clock and tep == "YES":
                    tep_count += 1
                if is_clock:
                    total_tests += 1
                
                status = "SUCCESS" if (is_clock and tep == "YES") else "INFO"
                print_status(f"  {mode_name:<12} {metric:<15} {lam:>8.0f}     {r2:>6.3f}   {tep:<6}", status)
    
    print_status("", "INFO")
    print_status("="*80, "INFO")
    
    # Final verdict
    if tep_count == total_tests and total_tests > 0:
        print_status(f"  TEP SIGNATURE CONFIRMED IN ALL {total_tests} MODES (clock_bias)", "SUCCESS")
    elif tep_count > 0:
        print_status(f"  TEP DETECTED in {tep_count}/{total_tests} modes (clock_bias)", "SUCCESS")
    else:
        print_status(f"  TEP NOT DETECTED in any mode", "WARNING")
    
    print_status("="*80, "INFO")
    print_status("Step 2.0 Analysis Complete", "TITLE")

def main():
    """Main execution loop over all filter modes."""
    print_status("STARTING MULTI-FILTER ANALYSIS PIPELINE", "TITLE")
    
    # Define the 3 standard filters
    filters = [
        ('ALL STATIONS', 'none'),
        ('OPTIMAL 100', 'optimal_100_metadata.json'),
        ('DYNAMIC 50', 'dynamic:50')
    ]
    
    for name, config in filters:
        try:
            run_analysis_pipeline(name, config)
        except Exception as e:
            print_status(f"Analysis failed for {name}: {e}", "ERROR")
            import traceback
            traceback.print_exc()
            
    print_status("MULTI-FILTER PIPELINE COMPLETE", "SUCCESS")

if __name__ == "__main__":
    main()
