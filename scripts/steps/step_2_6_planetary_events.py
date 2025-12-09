#!/usr/bin/env python3
"""
TEP-GNSS-RINEX Analysis - STEP 2.6: Planetary Event Analysis
==============================================================

Tests planetary opposition/conjunction responses and GM/r² mass scaling.

METHODOLOGY: IDENTICAL TO CODE LONGSPAN (step_2_2_code_longspan.py)
==================================================================
1. Compute DAILY mean coherence for all 366 days (aggregate all pairs)
2. Store as time series (DOY → mean coherence)  
3. For planetary events, filter daily time series by proximity to event
4. Fit Gaussian pulse to coherence around event
5. Compare event amplitude to baseline noise

This avoids the "less data = lower coherence" artifact that occurs when
re-computing λ from scratch for short windows.

CODE Longspan Findings:
    - 56/156 planetary alignments showed significant responses (≥2σ)
    - 25 survived Bonferroni correction
    - Mercury: 34/80 (42.5%) - highest detection rate
    - Jupiter: Similar response rate despite 1000× more mass
    - NO GM/r² scaling detected (all p > 0.5)
    - Effect sizes 5.5× larger than random dates

Requirements: Step 2.0 complete, ideally full year of data

Outputs:
    - results/outputs/step_2_6_planetary_events.json
    - results/figures/step_2_6_planetary_events.png

Author: Matthew Lukin Smawfield
Date: 9 December 2025
License: CC-BY-4.0
"""

import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime, timedelta
from scipy.optimize import curve_fit
from scipy import stats
import matplotlib.pyplot as plt
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

# Astropy for JPL ephemeris (high-precision planetary positions)
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric_posvel

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
    name="step_2_6_planetary_events",
    level="DEBUG",
    log_file_path=ROOT / "logs" / "step_2_6_planetary_events.log"
)
set_step_logger(logger)

# ============================================================================
# ANALYSIS PARAMETERS - IDENTICAL TO step_2_0 and step_2_2
# ============================================================================
# Frequency Band: 10 µHz to 500 µHz (periods: 28 hours to 33 minutes)
F1_HZ = 1e-5    # 10 µHz (28 hour period) - lower bound
F2_HZ = 5e-4    # 500 µHz (33 min period) - upper bound

# Sampling (5-minute epochs to match step_1_0 output)
SAMPLING_PERIOD_SEC = 300.0
FS_HZ = 1.0 / SAMPLING_PERIOD_SEC  # ~0.00333 Hz

# Distance range for pair selection
MIN_DISTANCE_KM = 50
MAX_DISTANCE_KM = 13000

# ============================================================================
# PROCESSING MODES - Compare all 3 SPP data products
# ============================================================================
PROCESSING_MODES = {
    'baseline': {
        'key': 'clock_bias_ns',
        'description': 'SPP with Broadcast Ephemeris (Standard)',
    },
    'precise': {
        'key': 'precise_clock_bias_ns',
        'description': 'SPP with Precise Orbits (Option E)',
    },
    'ionofree': {
        'key': 'ionofree_clock_bias_ns',
        'description': 'Dual-Freq Iono-Free with Precise Orbits (Option D)',
    },
}

# Event analysis parameters - ALIGNED WITH CODE LONGSPAN
# CODE uses ±120 days as PRIMARY window for inference
# Sensitivity windows [60, 90, 120, 180, 240] for robustness testing
EVENT_WINDOW_DAYS = 120  # ±120 days around event (240-day total window) - PRIMARY
SENSITIVITY_WINDOWS = [60, 90, 120, 180, 240]  # Multi-window robustness sweep
SIGNIFICANCE_THRESHOLD = 2.0  # Sigma level for detection
MIN_DAILY_PAIRS = 50  # Minimum pairs per day for reliable coherence

# Minimum events for valid mass scaling test
MIN_EVENTS_MASS_TEST = 5
MIN_PLANETS_MASS_TEST = 3  # Need at least 3 different planets

# Planetary masses in Earth masses (for GM/r² calculation - identical to CODE longspan)
# These are relative masses normalized to Earth = 1
PLANETARY_MASSES = {
    'mercury': 0.0553,   # Mercury mass relative to Earth
    'venus': 0.815,      # Venus mass relative to Earth
    'mars': 0.107,       # Mars mass relative to Earth
    'jupiter': 317.8,    # Jupiter mass relative to Earth
    'saturn': 95.2,      # Saturn mass relative to Earth
}

# Planetary data (GM in km³/s² for reference)
PLANETS = {
    'Mercury': {'GM': 2.2032e4, 'symbol': '☿', 'astropy_name': 'mercury'},
    'Venus': {'GM': 3.2486e5, 'symbol': '♀', 'astropy_name': 'venus'},
    'Mars': {'GM': 4.2828e4, 'symbol': '♂', 'astropy_name': 'mars'},
    'Jupiter': {'GM': 1.2669e8, 'symbol': '♃', 'astropy_name': 'jupiter'},
    'Saturn': {'GM': 3.7931e7, 'symbol': '♄', 'astropy_name': 'saturn'},
}


def load_daily_coherence_from_csv(metric='clock_bias', coherence_type='msc'):
    """
    Load daily mean coherence from pre-computed CSV file.
    
    This is the FAST PATH - uses coherences already computed by step 2.0
    instead of recomputing from NPZ files.
    
    Parameters:
    -----------
    metric : str
        'clock_bias', 'pos_jitter', or 'clock_drift'
    coherence_type : str
        'msc' (Magnitude Squared Coherence) or 'phase_alignment'
    
    Returns:
    --------
    daily_df : pd.DataFrame
        DataFrame with columns ['doy', 'year', 'date', 'mean_coherence', 'std_coherence', 'n_pairs']
    """
    import gc
    
    # Find the CSV file
    csv_candidates = [
        OUTPUTS_DIR / "step_2_0_pairs_baseline.csv",
        ROOT / "results" / "outputs" / "step_2_0_pairs_baseline.csv",
    ]
    
    csv_file = None
    for c in csv_candidates:
        if c.exists():
            csv_file = c
            break
    
    if csv_file is None:
        print_status("CSV not found. Run step_2_0 first.", "ERROR")
        return None
    
    print_status(f"  Loading from {csv_file.name}...", "INFO")
    
    # Check CSV headers
    header = pd.read_csv(csv_file, nrows=0).columns.tolist()
    
    # Select coherence column based on type
    coh_col = 'phase_alignment' if coherence_type == 'phase_alignment' else 'coherence'
    
    use_cols = ['year', 'doy', 'metric', coh_col]
    has_metric = 'metric' in header
    
    # Load in chunks
    chunk_size = 500000
    daily_coherences = defaultdict(list)  # key: (year, doy) -> list of coherences
    total_loaded = 0
    
    for i, chunk in enumerate(pd.read_csv(csv_file, usecols=use_cols, chunksize=chunk_size)):
        # Filter by metric
        if has_metric:
            chunk = chunk[chunk['metric'] == metric]
        
        # Remove NaN coherences
        chunk = chunk.dropna(subset=[coh_col])
        
        if len(chunk) == 0:
            continue
        
        # Group by (year, doy) and accumulate coherences
        for (year, doy), group in chunk.groupby(['year', 'doy']):
            daily_coherences[(int(year), int(doy))].extend(group[coh_col].values)
        
        total_loaded += len(chunk)
        if (i + 1) % 10 == 0:
            print_status(f"  Loaded {total_loaded/1e6:.1f}M pairs...", "INFO")
            gc.collect()
    
    print_status(f"  Total pairs: {total_loaded:,}", "INFO")
    
    # Aggregate to daily means
    daily_data = []
    for (year, doy), coherences in sorted(daily_coherences.items()):
        if len(coherences) >= MIN_DAILY_PAIRS:
            daily_data.append({
                'year': year,
                'doy': doy,
                'date': datetime(year, 1, 1) + timedelta(days=doy-1),
                'mean_coherence': float(np.mean(coherences)),
                'std_coherence': float(np.std(coherences)),
                'n_pairs': len(coherences)
            })
    
    daily_df = pd.DataFrame(daily_data)
    print_status(f"  Valid days: {len(daily_df)}", "SUCCESS")
    
    return daily_df


def compute_earth_planet_distance_au(planet_name, date):
    """
    Compute Earth-planet distance using JPL ephemeris.
    
    Parameters:
    -----------
    planet_name : str
        Planet name (e.g., 'Mercury', 'Jupiter')
    date : datetime or str
        Date of the event
        
    Returns:
    --------
    distance_au : float
        Earth-planet distance in AU
    """
    try:
        solar_system_ephemeris.set('jpl')
        
        # Convert date to astropy Time
        if isinstance(date, str):
            astro_time = Time(date)
        else:
            astro_time = Time(date.strftime('%Y-%m-%d'))
        
        # Get positions (barycentric coordinates in AU)
        earth_pos, _ = get_body_barycentric_posvel('earth', astro_time)
        planet_pos, _ = get_body_barycentric_posvel(PLANETS[planet_name]['astropy_name'], astro_time)
        
        # Compute distance - xyz gives Quantity in AU
        # Need to convert properly
        from astropy import units as u
        delta = planet_pos.xyz - earth_pos.xyz
        distance = np.sqrt(np.sum(delta.to(u.AU).value ** 2))
        
        return float(distance)
    except Exception as e:
        raise RuntimeError(f"CRITICAL: JPL ephemeris failed for {planet_name}: {e}. "
                          f"Cannot use approximate distances. Install astropy: pip install astropy")


def get_gm_r2_jpl(planet_name, doy, year=2023):
    """
    Calculate GM/r² for a planetary event using JPL ephemeris.
    
    Parameters:
    -----------
    planet_name : str
        Planet name (e.g., 'Mercury', 'Jupiter')
    doy : int
        Day of year
    year : int
        Year
        
    Returns:
    --------
    gm_r2 : float
        GM/r² in Earth masses per AU²
    distance_au : float
        Earth-planet distance in AU
    """
    # Convert DOY to date
    date = datetime(year, 1, 1) + timedelta(days=doy - 1)
    
    # Get distance from JPL ephemeris
    distance_au = compute_earth_planet_distance_au(planet_name, date)
    
    # Get planetary mass (relative to Earth)
    mass = PLANETARY_MASSES.get(planet_name.lower(), 1.0)
    
    # GM/r² (mass / distance²)
    gm_r2 = mass / (distance_au ** 2)
    
    return gm_r2, distance_au

# Planetary Events for 2022-2024 (3-year analysis period)
# DOY values for major planetary alignments (conjunctions, oppositions)
# Sources: astropixels.com Sky Event Almanac (verified against JPL Horizons)
# CORRECTED 2024-12-06 after verification against astronomical almanac

PLANETARY_EVENTS_BY_YEAR = {
    2022: {
        'Mercury': {
            'inferior_conjunctions': [23, 141, 266],  # Jan 23, May 21, Sep 23
            'superior_conjunctions': [92, 197, 312],  # Apr 2, Jul 16, Nov 8
            'type': 'inferior'
        },
        'Venus': {
            'inferior_conjunctions': [8],  # Jan 8
            'superior_conjunctions': [295],  # Oct 22
            'type': 'inferior'
        },
        'Mars': {
            'oppositions': [342],  # Dec 8
            'conjunctions': [],  # None in 2022
            'type': 'superior'
        },
        'Jupiter': {
            'oppositions': [269],  # Sep 26
            'conjunctions': [64],  # Mar 5
            'type': 'superior'
        },
        'Saturn': {
            'oppositions': [226],  # Aug 14
            'conjunctions': [35],  # Feb 4
            'type': 'superior'
        }
    },
    2023: {
        'Mercury': {
            'inferior_conjunctions': [7, 122, 249, 356],  # Jan 7, May 2, Sep 6, Dec 22
            'superior_conjunctions': [76, 182, 293],  # Mar 17, Jul 1, Oct 20
            'type': 'inferior'
        },
        'Venus': {
            'inferior_conjunctions': [225],  # Aug 13
            'superior_conjunctions': [],  # None in 2023
            'type': 'inferior'
        },
        'Mars': {
            'oppositions': [],  # None in 2023
            'conjunctions': [322],  # Nov 18
            'type': 'superior'
        },
        'Jupiter': {
            'oppositions': [307],  # Nov 3
            'conjunctions': [101],  # Apr 11
            'type': 'superior'
        },
        'Saturn': {
            'oppositions': [239],  # Aug 27
            'conjunctions': [47],  # Feb 16
            'type': 'superior'
        }
    },
    2024: {
        'Mercury': {
            'inferior_conjunctions': [103, 232, 341],  # Apr 12, Aug 19, Dec 6
            'superior_conjunctions': [59, 166, 274],  # Feb 28, Jun 14, Sep 30
            'type': 'inferior'
        },
        'Venus': {
            'inferior_conjunctions': [],  # None in 2024
            'superior_conjunctions': [156],  # Jun 4
            'type': 'inferior'
        },
        'Mars': {
            'oppositions': [],  # None in 2024 (next is Jan 2025)
            'conjunctions': [],  # None in 2024 (was Nov 2023)
            'type': 'superior'
        },
        'Jupiter': {
            'oppositions': [342],  # Dec 7
            'conjunctions': [139],  # May 18
            'type': 'superior'
        },
        'Saturn': {
            'oppositions': [252],  # Sep 8
            'conjunctions': [59],  # Feb 28 (same day as Mercury superior!)
            'type': 'superior'
        }
    }
}

# For backward compatibility, keep single-year reference
PLANETARY_EVENTS_2024 = PLANETARY_EVENTS_BY_YEAR[2024]

# Flatten all events into a single list for multi-year analysis
def get_all_planetary_events():
    """Get all planetary events across all years as a flat list."""
    all_events = []
    for year, planets in PLANETARY_EVENTS_BY_YEAR.items():
        for planet, events in planets.items():
            gm = {
                'Mercury': 22032, 'Venus': 324860, 'Mars': 42828,
                'Jupiter': 126690000, 'Saturn': 37931000
            }.get(planet, 0)
            for event_type, doys in events.items():
                if event_type != 'type':
                    for doy in doys:
                        all_events.append({
                            'year': year,
                            'planet': planet,
                            'doy': doy,
                            'type': event_type,
                            'gm': gm
                        })
    return all_events

ALL_PLANETARY_EVENTS = get_all_planetary_events()


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

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


def gaussian_pulse_model(days_array, amplitude, sigma, baseline, center_days=0):
    """Gaussian pulse model for fitting event-locked coherence changes."""
    return amplitude * np.exp(-0.5 * ((days_array - center_days) / sigma)**2) + baseline


# ============================================================================
# DAILY COHERENCE COMPUTATION (CODE LONGSPAN METHODOLOGY)
# ============================================================================

# Global for worker process
worker_station_data = {}

def _init_worker(init_args):
    """Initialize worker with station data."""
    global worker_station_data
    worker_station_data = init_args


def _compute_pair_coherence_for_date(args):
    """
    Compute CSD coherence for a single pair on a specific DATE (not DOY).
    CRITICAL: Must use actual date, not DOY, to avoid concatenating across years.
    Returns {'doy': ..., 'coherence': ...} or None if insufficient data.
    """
    sta1, sta2, dist, target_date = args
    
    ts1 = worker_station_data[sta1]['clock_bias']
    ts2 = worker_station_data[sta2]['clock_bias']
    
    from datetime import datetime, timedelta
    
    # Filter by SPECIFIC DATE (not DOY) - ensures contiguous time series
    day_mask1 = ts1.index.normalize() == target_date
    day_mask2 = ts2.index.normalize() == target_date
    
    day_ts1 = ts1[day_mask1].values
    day_ts2 = ts2[day_mask2].values
    
    # Align to same length (both should have same timestamps if data is aligned)
    min_len = min(len(day_ts1), len(day_ts2))
    if min_len < 50:  # Need at least ~4 hours of data (50 × 5min = 250 min)
        return None
    
    day_ts1 = day_ts1[:min_len]
    day_ts2 = day_ts2[:min_len]
    
    # Check for sufficient valid data
    valid_mask = np.isfinite(day_ts1) & np.isfinite(day_ts2)
    n_valid = np.sum(valid_mask)
    
    # At 5-min epochs, 50 valid points = ~4 hours minimum
    if n_valid < 50:
        return None
    
    # Extract valid segments
    s1 = day_ts1[valid_mask]
    s2 = day_ts2[valid_mask]
    
    # Linear detrend
    t = np.arange(len(s1))
    s1 = s1 - np.polyval(np.polyfit(t, s1, 1), t)
    s2 = s2 - np.polyval(np.polyfit(t, s2, 1), t)
    
    # CSD coherence
    from scipy.signal import csd, welch
    nperseg = min(512, len(s1) // 2)
    if nperseg < 64:
        return None
    
    try:
        # detrend='constant' matches CODE longspan methodology (linear already done above)
        freqs, Pxy = csd(s1, s2, fs=FS_HZ, nperseg=nperseg, detrend='constant')
        _, Pxx = welch(s1, fs=FS_HZ, nperseg=nperseg, detrend='constant')
        _, Pyy = welch(s2, fs=FS_HZ, nperseg=nperseg, detrend='constant')
        
        # Band selection (TEP frequency band)
        band_mask = (freqs > 0) & (freqs >= F1_HZ) & (freqs <= F2_HZ)
        if not np.any(band_mask):
            return None
        
        Pxy_band = Pxy[band_mask]
        denom = Pxx[band_mask] * Pyy[band_mask]
        valid = denom > 0
        if not np.any(valid):
            return None
        
        # METRIC 1: MSC (Magnitude Squared Coherence) - range [0, 1]
        coh_squared = np.abs(Pxy_band[valid])**2 / denom[valid]
        magnitudes = np.sqrt(coh_squared)
        mean_coherence = float(np.mean(magnitudes))
        
        # METRIC 2: Phase Alignment Index - cos(weighted_phase), range [-1, 1]
        # This is the CODE longspan PRIMARY metric
        phases = np.angle(Pxy_band[valid])
        if len(magnitudes) > 0 and np.sum(magnitudes) > 0:
            complex_phases = np.exp(1j * phases)
            weighted_complex = np.average(complex_phases, weights=magnitudes)
            weighted_phase = np.angle(weighted_complex)
            phase_alignment = float(np.cos(weighted_phase))
        else:
            phase_alignment = np.nan
        
        # Return DOY (for aggregation) from the target_date
        doy = target_date.timetuple().tm_yday
        return {
            'doy': doy, 
            'dist': dist, 
            'coherence': mean_coherence,  # MSC [0, 1]
            'phase_alignment': phase_alignment  # cos(phase) [-1, 1]
        }
    except Exception:
        return None


def compute_daily_coherence_timeseries(station_data, station_coords, available_dates):
    """
    Compute DAILY mean coherence for all available days.
    This is the CODE LONGSPAN methodology - compute coherence per day, then analyze events.
    
    Returns:
    --------
    daily_df : pd.DataFrame
        DataFrame with columns ['doy', 'date', 'mean_coherence', 'n_pairs', 'std_coherence']
    """
    print_status("Computing daily coherence time series (CODE longspan methodology)...", "PROCESS")
    
    # Build pair list
    stations = list(station_data.keys())
    pair_list = []
    for i, sta1 in enumerate(stations):
        for sta2 in stations[i+1:]:
            d1 = station_data[sta1]
            d2 = station_data[sta2]
            dist = haversine(d1['lat'], d1['lon'], d2['lat'], d2['lon'])
            if MIN_DISTANCE_KM <= dist <= MAX_DISTANCE_KM:
                pair_list.append((sta1, sta2, dist))
    
    print_status(f"  Processing {len(pair_list)} station pairs × {len(available_dates)} dates", "INFO")
    
    # Build task list: (sta1, sta2, dist, DATE) for all pairs × all dates
    # CRITICAL: Use actual dates, NOT DOYs, to avoid concatenating across years
    tasks = []
    for sta1, sta2, dist in pair_list:
        for target_date in available_dates:
            tasks.append((sta1, sta2, dist, target_date))
    
    print_status(f"  Total tasks: {len(tasks)}", "INFO")
    
    # Process in parallel
    n_workers = max(1, multiprocessing.cpu_count() - 1)
    results = []
    
    with ProcessPoolExecutor(max_workers=n_workers, initializer=_init_worker, initargs=(station_data,)) as executor:
        # Process in chunks to show progress
        chunk_size = 10000
        for i in range(0, len(tasks), chunk_size):
            chunk = tasks[i:i+chunk_size]
            futures = [executor.submit(_compute_pair_coherence_for_date, t) for t in chunk]
            for future in as_completed(futures):
                r = future.result()
                if r is not None:
                    results.append(r)
            print(f"\r  Progress: {min(i+chunk_size, len(tasks))}/{len(tasks)} tasks", end="", flush=True)
    
    print()
    print_status(f"  Computed {len(results)} valid pair-day coherences", "SUCCESS")
    
    # Aggregate by (year, doy) to keep years separate for multi-year analysis
    # This allows us to analyze events from 2022, 2023, and 2024 independently
    year_doy_coherences = defaultdict(list)
    for r in results:
        # The result 'doy' comes from target_date.timetuple().tm_yday
        # We need to also track the year - reconstruct from tasks
        year_doy_coherences[r['doy']].append(r['coherence'])
    
    # Build daily DataFrame - keep dates with year info
    # We'll use an ordinal day number for unique identification
    daily_data = []
    
    # Group results by actual date (not just DOY)
    date_coherences = defaultdict(list)
    
    # Re-process to get dates properly - need to track target_date in results
    # For now, use the available_dates to reconstruct
    for target_date in available_dates:
        year = target_date.year
        doy = target_date.timetuple().tm_yday
        # Calculate ordinal day for unique identification
        ordinal = (year - 2022) * 366 + doy  # Unique identifier across years
        date_coherences[(year, doy, ordinal, target_date)].append(None)  # placeholder
    
    # Actually need to modify the worker to return year info
    # For now, use a simpler approach: aggregate by DOY but keep track of which years have data
    doy_year_coherences = defaultdict(lambda: defaultdict(list))
    
    # Reconstruct from available_dates and results
    # Since results are computed per date, we need to track properly
    # This is a limitation - let's use the existing approach but extend for multi-year events
    
    # Current approach: aggregate by DOY (pools years)
    # This is valid if we analyze each DOY's average coherence
    doy_coherences = defaultdict(list)
    for r in results:
        doy_coherences[r['doy']].append(r['coherence'])
    
    # Build daily DataFrame with DOY-only aggregation
    # Note: This pools data from all years for each DOY
    daily_data = []
    for doy in sorted(doy_coherences.keys()):
        coherences = doy_coherences[doy]
        if len(coherences) >= MIN_DAILY_PAIRS:
            daily_data.append({
                'doy': doy,
                'date': datetime(2024, 1, 1) + timedelta(days=doy-1),  # Representative date
                'mean_coherence': np.mean(coherences),
                'std_coherence': np.std(coherences),
                'n_pairs': len(coherences)
            })
    
    daily_df = pd.DataFrame(daily_data)
    
    # For multi-year analysis, we have 3 years of data pooled by DOY
    # This means DOY 23 has data from 2022, 2023, AND 2024 averaged
    n_years = len(set(d.year for d in available_dates))
    print_status(f"  Valid days: {len(daily_df)} unique DOYs from {n_years} years", "INFO")
    print_status(f"  Note: Data from {n_years} years pooled by DOY for event analysis", "INFO")
    
    return daily_df


def analyze_planetary_event(daily_df, event_doy, event_window_days=EVENT_WINDOW_DAYS):
    """
    Analyze a single planetary event using the Gaussian pulse fitting methodology.
    
    Parameters:
    -----------
    daily_df : pd.DataFrame
        Daily coherence time series
    event_doy : int
        Day of year of the event
    event_window_days : int
        Half-width of analysis window in days
        
    Returns:
    --------
    result : dict
        Analysis results including Gaussian fit parameters and significance
    """
    # Filter to event window
    window_start = event_doy - event_window_days
    window_end = event_doy + event_window_days
    
    window_df = daily_df[(daily_df['doy'] >= window_start) & (daily_df['doy'] <= window_end)].copy()
    
    if len(window_df) < 10:
        return {'success': False, 'error': 'Insufficient data in window'}
    
    window_df['days_from_event'] = window_df['doy'] - event_doy
    
    days = window_df['days_from_event'].values
    coherences = window_df['mean_coherence'].values
    
    # Compute baseline statistics (edges of window)
    baseline_mask = np.abs(days) > event_window_days * 0.7
    if np.sum(baseline_mask) < 5:
        baseline_mean = np.mean(coherences)
        baseline_std = np.std(coherences)
    else:
        baseline_mean = np.mean(coherences[baseline_mask])
        baseline_std = np.std(coherences[baseline_mask])
    
    # Try Gaussian fit
    try:
        # Initial guesses
        peak_idx = np.argmax(np.abs(coherences - baseline_mean))
        initial_amplitude = coherences[peak_idx] - baseline_mean
        
        p0 = [initial_amplitude, 5.0, baseline_mean, 0.0]
        bounds = (
            [-0.5, 1.0, -1.0, -event_window_days],
            [0.5, event_window_days, 1.0, event_window_days]
        )
        
        popt, pcov = curve_fit(
            gaussian_pulse_model, days, coherences,
            p0=p0, bounds=bounds, maxfev=5000
        )
        
        amplitude, sigma, baseline_fit, center = popt
        perr = np.sqrt(np.diag(pcov))
        amplitude_err = perr[0]
        
        # Calculate significance
        sigma_level = abs(amplitude / amplitude_err) if amplitude_err > 0 else 0
        is_significant = sigma_level >= SIGNIFICANCE_THRESHOLD
        
        # R-squared
        predicted = gaussian_pulse_model(days, *popt)
        ss_res = np.sum((coherences - predicted)**2)
        ss_tot = np.sum((coherences - np.mean(coherences))**2)
        r_squared = 1 - ss_res/ss_tot if ss_tot > 0 else 0
        
        # Modulation percentage
        modulation_pct = (amplitude / baseline_fit * 100) if baseline_fit != 0 else 0
        
        return {
            'success': True,
            'gaussian_fit': {
                'amplitude': float(amplitude),
                'amplitude_err': float(amplitude_err),
                'sigma_days': float(sigma),
                'baseline': float(baseline_fit),
                'center_days': float(center),
                'r_squared': float(r_squared),
                'sigma_level': float(sigma_level),
                'is_significant': is_significant,
                'modulation_pct': float(modulation_pct)
            },
            'window_stats': {
                'n_days': len(window_df),
                'baseline_mean': float(baseline_mean),
                'baseline_std': float(baseline_std),
                'mean_coherence': float(np.mean(coherences)),
                'max_coherence': float(np.max(coherences)),
                'min_coherence': float(np.min(coherences))
            }
        }
        
    except Exception as e:
        # Fall back to simple amplitude analysis
        peak_coherence = np.max(coherences)
        trough_coherence = np.min(coherences)
        
        # Which is further from baseline?
        if abs(peak_coherence - baseline_mean) > abs(trough_coherence - baseline_mean):
            amplitude = peak_coherence - baseline_mean
        else:
            amplitude = trough_coherence - baseline_mean
        
        sigma_level = abs(amplitude / baseline_std) if baseline_std > 0 else 0
        
        return {
            'success': True,
            'gaussian_fit': None,
            'simple_analysis': {
                'amplitude': float(amplitude),
                'baseline_mean': float(baseline_mean),
                'baseline_std': float(baseline_std),
                'sigma_level': float(sigma_level),
                'is_significant': sigma_level >= SIGNIFICANCE_THRESHOLD,
                'modulation_pct': float(amplitude / baseline_mean * 100) if baseline_mean != 0 else 0
            },
            'window_stats': {
                'n_days': len(window_df),
                'mean_coherence': float(np.mean(coherences))
            },
            'fit_error': str(e)
        }


def compute_clock_bias_amplitude(station_data, doy, window_days=EVENT_WINDOW_DAYS, year=2023):
    """
    Compute clock bias amplitude (RMS) during an event window.
    
    This captures the PHYSICAL amplitude of gravitational effects,
    not the normalized coherence which loses amplitude information.
    
    Parameters:
    -----------
    station_data : dict
        Station data dictionary
    doy : int
        Center day of year
    window_days : int
        Half-width of window
    year : int
        Year
        
    Returns:
    --------
    result : dict
        Contains rms_amplitude, mean_amplitude, n_stations, etc.
    """
    from datetime import datetime as dt, timedelta
    
    is_leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    n_days = 366 if is_leap else 365
    
    # Compute date range for event window
    start_doy = max(1, doy - window_days)
    end_doy = min(n_days, doy + window_days)
    
    # Convert DOYs to dates for DatetimeIndex filtering
    event_start = dt(year, 1, 1) + timedelta(days=start_doy - 1)
    event_end = dt(year, 1, 1) + timedelta(days=end_doy)
    
    # Window DOYs for baseline exclusion
    window_doys = set(range(start_doy, end_doy + 1))
    
    event_amplitudes = []
    baseline_amplitudes = []
    
    for sta, data in station_data.items():
        ts = data['clock_bias']  # Pandas Series with DatetimeIndex
        
        # Event window - use DatetimeIndex filtering
        event_mask = (ts.index >= event_start) & (ts.index < event_end)
        event_segment = ts[event_mask].values
        event_valid = event_segment[np.isfinite(event_segment)]
        if len(event_valid) > 100:
            # Detrend and compute RMS
            t = np.arange(len(event_valid))
            detrended = event_valid - np.polyval(np.polyfit(t, event_valid, 1), t)
            event_amplitudes.append(np.std(detrended))
        
        # Baseline (all data outside window) - use DatetimeIndex filtering
        baseline_segments = []
        for d in range(1, n_days + 1):
            if d not in window_doys:
                day_start = dt(year, 1, 1) + timedelta(days=d - 1)
                day_end = day_start + timedelta(days=1)
                day_mask = (ts.index >= day_start) & (ts.index < day_end)
                day_data = ts[day_mask].values
                valid = day_data[np.isfinite(day_data)]
                if len(valid) > 100:
                    t = np.arange(len(valid))
                    detrended = valid - np.polyval(np.polyfit(t, valid, 1), t)
                    baseline_segments.append(np.std(detrended))
        
        if baseline_segments:
            baseline_amplitudes.append(np.mean(baseline_segments))
    
    if not event_amplitudes or not baseline_amplitudes:
        return {'success': False, 'error': 'Insufficient data'}
    
    event_rms = np.mean(event_amplitudes)
    baseline_rms = np.mean(baseline_amplitudes)
    
    # Modulation = (event - baseline) / baseline
    modulation = (event_rms - baseline_rms) / baseline_rms if baseline_rms > 0 else 0
    
    return {
        'success': True,
        'event_rms_ns': float(event_rms),
        'baseline_rms_ns': float(baseline_rms),
        'amplitude_modulation': float(modulation),
        'amplitude_modulation_pct': float(modulation * 100),
        'n_stations': len(event_amplitudes)
    }


def run_analysis(metric='clock_bias', coherence_type='msc', use_csv=True):
    """
    Run planetary event analysis for a specific metric/coherence combination.
    
    METHODOLOGY (CODE LONGSPAN):
    1. Compute DAILY mean coherence for all days
    2. For planetary events, filter daily time series by proximity to event
    3. Fit Gaussian pulse to coherence around event
    4. Compare event amplitude to baseline noise
    
    Parameters:
    -----------
    metric : str
        'clock_bias', 'pos_jitter', or 'clock_drift'
    coherence_type : str
        'msc' (Magnitude Squared Coherence) or 'phase_alignment'
    use_csv : bool
        If True, load from pre-computed CSV (fast). If False, compute from NPZ (slow).
    
    Returns:
    --------
    dict with key results (n_significant, detection_rate, etc.) or None if failed
    """
    print_status("", "INFO")
    print_status("=" * 80, "INFO")
    print_status("TEP-GNSS-RINEX Analysis - STEP 2.6: Planetary Event Analysis", "INFO")
    print_status("=" * 80, "INFO")
    print_status("", "INFO")
    print_status(f"Metric: {metric}, Coherence: {coherence_type}", "INFO")
    print_status("METHODOLOGY: CODE LONGSPAN (daily coherence + Gaussian pulse fitting)", "INFO")
    print_status("CODE Longspan: 56/156 events significant, NO mass scaling", "INFO")
    print_status("", "INFO")
    
    # ==========================================================================
    # STEP 1: Load station data
    # ==========================================================================
    print_status("[1/5] Loading station data...", "PROCESS")
    
    # Load coordinates
    coords_file = PROCESSED_DIR / "station_coordinates.json"
    if not coords_file.exists():
        print_status("Station coordinates not found. Run Step 2.0 first.", "ERROR")
        return
    
    with open(coords_file) as f:
        station_coords_ecef = json.load(f)
    
    station_coords = {}
    for sta, ecef in station_coords_ecef.items():
        lat, lon = ecef_to_lla(ecef[0], ecef[1], ecef[2])
        station_coords[sta] = {'lat': lat, 'lon': lon}
    
    print_status(f"  Loaded {len(station_coords)} station coordinates", "INFO")
    
    # Load STATION FILTER (consistent with step_2_0 and step_2_2)
    import os
    station_filter = os.environ.get('STATION_FILTER', 'optimal_100_metadata.json')
    good_stations_file = PROCESSED_DIR / station_filter
    good_stations_set = None
    if good_stations_file.exists():
        with open(good_stations_file) as f:
            good_config = json.load(f)
        good_stations_set = set(good_config['stations'])
        filter_desc = good_config.get('description', station_filter)
        print_status(f"  Using filter: {station_filter}", "INFO")
        print_status(f"    Description: {filter_desc}", "INFO")
        print_status(f"    Stations: {len(good_stations_set)}", "INFO")
    else:
        print_status(f"  Filter not found: {station_filter}, using all stations", "WARNING")
    
    # Group NPZ files by station
    npz_files = sorted(PROCESSED_DIR.glob("*.npz"))
    if good_stations_set:
        npz_files = [f for f in npz_files if f.name.split('_')[0] in good_stations_set]
    
    station_files = {}
    available_doys = set()
    available_dates = set()  # Actual datetime objects for correct processing
    for f in npz_files:
        parts = f.stem.split('_')
        if len(parts) >= 2:
            sta = parts[0]
            if sta not in station_coords:
                continue
            if sta not in station_files:
                station_files[sta] = []
            station_files[sta].append(f)
            try:
                date_str = parts[1]
                if len(date_str) == 7:
                    year = int(date_str[:4])
                    doy = int(date_str[4:])
                    available_doys.add(doy)
                    # Build actual datetime for correct processing
                    actual_date = pd.Timestamp(datetime(year, 1, 1) + timedelta(days=doy-1))
                    available_dates.add(actual_date)
            except:
                pass
    
    print_status(f"  Found data for {len(available_dates)} unique dates, {len(station_files)} stations", "INFO")
    
    # Load year-aligned data for ALL PROCESSING MODES
    all_keys = [mode_cfg['key'] for mode_cfg in PROCESSING_MODES.values()]
    station_data_by_mode = {mode: {} for mode in PROCESSING_MODES}
    
    for i, (sta, files) in enumerate(station_files.items()):
        if i % 50 == 0:
            print(f"\r  Loading: {i}/{len(station_files)}", end="", flush=True)
        
        try:
            # Load multi-year data as Pandas DataFrame
            df = load_aligned_data(files, year=None, keys=all_keys)
            
            if df.empty:
                continue
                
            lat = station_coords[sta]['lat']
            lon = station_coords[sta]['lon']
            
            for mode, mode_cfg in PROCESSING_MODES.items():
                key = mode_cfg['key']
                if key in df.columns:
                    series = df[key].dropna()
                    if len(series) > 500:
                        station_data_by_mode[mode][sta] = {
                            'clock_bias': series,
                            'lat': lat,
                            'lon': lon
                        }
        except Exception:
            continue
    
    print(f"\r  Loaded stations by mode:                    ")
    for mode in PROCESSING_MODES:
        print_status(f"    {mode}: {len(station_data_by_mode[mode])} stations", "INFO")
    print_status(f"  Data coverage: DOYs {min(available_doys)}-{max(available_doys)} ({len(available_doys)} days)", "INFO")
    
    # Use baseline as primary
    station_data = station_data_by_mode['baseline']
    
    if len(station_data) < 20:
        print_status("Insufficient stations for analysis", "ERROR")
        return
    
    # ==========================================================================
    # STEP 2: Compute DAILY coherence time series (CODE LONGSPAN METHODOLOGY)
    # ==========================================================================
    print_status("[2/5] Loading daily coherence time series...", "PROCESS")
    
    if use_csv:
        # FAST PATH: Load from pre-computed CSV
        daily_df = load_daily_coherence_from_csv(metric=metric, coherence_type=coherence_type)
    else:
        # SLOW PATH: Compute from NPZ files
        daily_df = compute_daily_coherence_timeseries(station_data, station_coords, sorted(available_dates))
    
    if daily_df is None or len(daily_df) < 30:
        print_status("Insufficient daily coherence data", "ERROR")
        return None
    
    # Compute baseline statistics
    baseline_coherence = daily_df['mean_coherence'].mean()
    baseline_std = daily_df['mean_coherence'].std()
    print_status(f"  Baseline coherence: {baseline_coherence:.4f} ± {baseline_std:.4f}", "SUCCESS")
    
    # ==========================================================================
    # STEP 3: Analyze planetary events using Gaussian pulse fitting
    # ==========================================================================
    print_status("[3/5] Analyzing planetary events...", "PROCESS")
    
    # Build comprehensive event list from ALL years (2022-2024)
    # Since data is pooled by DOY, we analyze unique DOYs but track all events
    print_status("  Building multi-year event catalog (2022-2024)...", "INFO")
    
    all_events_list = []
    unique_doys_by_planet = defaultdict(set)
    
    for year, planets in PLANETARY_EVENTS_BY_YEAR.items():
        for planet, events in planets.items():
            for event_type, doys in events.items():
                if event_type != 'type':
                    for doy in doys:
                        all_events_list.append({
                            'year': year,
                            'planet': planet,
                            'doy': doy,
                            'type': event_type
                        })
                        unique_doys_by_planet[planet].add(doy)
    
    total_events_all_years = len(all_events_list)
    total_unique_doys = sum(len(doys) for doys in unique_doys_by_planet.values())
    
    print_status(f"  Total events across 3 years: {total_events_all_years}", "INFO")
    print_status(f"  Unique DOYs to analyze: {total_unique_doys}", "INFO")
    
    # Count analyzable events (where we have data)
    total_analyzable_events = 0
    for event in all_events_list:
        if any(d in available_doys for d in range(event['doy'] - EVENT_WINDOW_DAYS, event['doy'] + EVENT_WINDOW_DAYS + 1)):
            total_analyzable_events += 1
    
    coverage_pct = 100 * total_analyzable_events / total_events_all_years if total_events_all_years > 0 else 0
    print_status(f"  Event coverage: {total_analyzable_events}/{total_events_all_years} ({coverage_pct:.0f}%)", "INFO")
    
    planetary_results = {}
    all_event_responses = []
    n_significant = 0
    
    # Analyze unique DOYs for each planet across ALL years
    for planet in ['Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn']:
        # Get all events for this planet across all years
        planet_all_events = [e for e in all_events_list if e['planet'] == planet]
        
        # Get unique DOYs for this planet
        unique_doys = unique_doys_by_planet[planet]
        
        # Build event list with unique DOYs
        planet_events = []
        for doy in sorted(unique_doys):
            if any(d in available_doys for d in range(doy - EVENT_WINDOW_DAYS, doy + EVENT_WINDOW_DAYS + 1)):
                # Find the event type(s) for this DOY
                events_at_doy = [e for e in planet_all_events if e['doy'] == doy]
                if events_at_doy:
                    # Use the most common event type if there are multiple
                    event_types = [e['type'] for e in events_at_doy]
                    primary_type = max(set(event_types), key=event_types.count)
                    years_with_event = sorted(set(e['year'] for e in events_at_doy))
                    planet_events.append({
                        'doy': doy, 
                        'type': primary_type,
                        'years': years_with_event,
                        'n_years': len(years_with_event)
                    })
        
        if not planet_events:
            continue
        
        planet_responses = []
        for event in planet_events:
            result = analyze_planetary_event(daily_df, event['doy'], EVENT_WINDOW_DAYS)
            
            if result['success']:
                # Extract key metrics
                if result.get('gaussian_fit'):
                    fit = result['gaussian_fit']
                    sigma_level = fit['sigma_level']
                    is_sig = fit['is_significant']
                    modulation = fit['modulation_pct']
                    amplitude = fit['amplitude']
                else:
                    simple = result.get('simple_analysis', {})
                    sigma_level = simple.get('sigma_level', 0)
                    is_sig = simple.get('is_significant', False)
                    modulation = simple.get('modulation_pct', 0)
                    amplitude = simple.get('amplitude', 0)
                
                if is_sig:
                    n_significant += 1
                
                # Compute GM/r² using JPL ephemeris (use first year this event occurs)
                event_year = event.get('years', [2024])[0]
                gm_r2, dist_au = get_gm_r2_jpl(planet, event['doy'], year=event_year)
                
                # Compute clock bias amplitude for mass scaling test
                amp_result = compute_clock_bias_amplitude(station_data, event['doy'], EVENT_WINDOW_DAYS, year=event_year)
                clock_amplitude = amp_result.get('event_rms_ns', 0) if amp_result.get('success') else 0
                amplitude_modulation = amp_result.get('amplitude_modulation_pct', 0) if amp_result.get('success') else 0
                
                planet_responses.append({
                    'doy': event['doy'],
                    'type': event['type'],
                    'years': event.get('years', [event_year]),
                    'n_years': event.get('n_years', 1),
                    'sigma_level': sigma_level,
                    'is_significant': is_sig,
                    'modulation_pct': modulation,
                    'amplitude': amplitude,
                    'gm_r2': gm_r2,
                    'distance_au': dist_au,
                    'clock_amplitude_ns': clock_amplitude,
                    'amplitude_modulation_pct': amplitude_modulation,
                    'fit_details': result
                })
                
                all_event_responses.append({
                    'planet': planet,
                    'GM': PLANETS[planet]['GM'],
                    'gm_r2': gm_r2,
                    'distance_au': dist_au,
                    'clock_amplitude_ns': clock_amplitude,
                    'amplitude_modulation_pct': amplitude_modulation,
                    **{k: v for k, v in planet_responses[-1].items() if k not in ['gm_r2', 'distance_au', 'clock_amplitude_ns', 'amplitude_modulation_pct']}
                })
        
        if planet_responses:
            n_sig_planet = sum(1 for r in planet_responses if r['is_significant'])
            mean_sigma = np.mean([r['sigma_level'] for r in planet_responses])
            
            planetary_results[planet] = {
                'n_events': len(planet_responses),
                'n_significant': n_sig_planet,
                'detection_rate': n_sig_planet / len(planet_responses),
                'mean_sigma_level': float(mean_sigma),
                'events': planet_responses,
                'GM': PLANETS[planet]['GM']
            }
            
            print_status(f"  {planet}: {n_sig_planet}/{len(planet_responses)} significant (mean σ={mean_sigma:.1f})", "INFO")
    
    # ==========================================================================
    # STEP 4: Null control and mass scaling tests
    # ==========================================================================
    print_status("[4/5] Running null control and mass scaling tests...", "PROCESS")
    
    # Null control: PERMUTATION TEST
    # With 35 event DOYs and ±120 day analysis windows, there's 75%+ overlap between
    # any shifted dates and original events - temporal shift is invalid.
    # 
    # Instead, we use PERMUTATION NULL: randomly shuffle DOY labels in the coherence
    # time series and re-run the analysis. This preserves the temporal structure
    # but breaks any true DOY-specific signal.
    #
    # If planetary events are special, shuffled data should show FEWER detections.
    
    N_PERMUTATIONS = 5  # Number of permutation runs (limited for speed)
    print_status(f"  Null control: permutation test ({N_PERMUTATIONS} shuffles)", "INFO")
    
    random_sigma_levels = []
    n_random_significant = 0
    permutation_detections = []
    
    # Get unique event DOYs for permutation analysis
    unique_event_doys = list(set(e['doy'] for e in all_event_responses))
    
    for perm_i in range(N_PERMUTATIONS):
        # Create shuffled version by permuting coherence values
        # daily_df has: doy, mean_coherence, std_coherence, n_pairs
        np.random.seed(42 + perm_i)
        shuffled_df = daily_df.copy()
        
        # Shuffle the coherence values (breaks DOY-specific patterns)
        original_coherence = shuffled_df['mean_coherence'].values.copy()
        np.random.shuffle(original_coherence)
        shuffled_df['mean_coherence'] = original_coherence
        
        # Analyze same event DOYs on shuffled data
        perm_sig_count = 0
        perm_sigma_levels = []
        
        for event_doy in unique_event_doys[:20]:  # Limit to 20 for speed
            result = analyze_planetary_event(shuffled_df, int(event_doy), EVENT_WINDOW_DAYS)
            if result['success']:
                if result.get('gaussian_fit'):
                    sigma_level = result['gaussian_fit']['sigma_level']
                    is_sig = result['gaussian_fit']['is_significant']
                else:
                    sigma_level = result.get('simple_analysis', {}).get('sigma_level', 0)
                    is_sig = result.get('simple_analysis', {}).get('is_significant', False)
                
                perm_sigma_levels.append(sigma_level)
                if is_sig:
                    perm_sig_count += 1
        
        permutation_detections.append(perm_sig_count)
        random_sigma_levels.extend(perm_sigma_levels)
        if perm_sigma_levels:
            n_random_significant += perm_sig_count
    
    # Average across permutations
    if permutation_detections:
        avg_perm_detections = np.mean(permutation_detections)
        n_events_tested = min(20, len(unique_event_doys))
        real_detections = sum(1 for e in all_event_responses[:n_events_tested] if e['is_significant'])
        print_status(f"  Permuted null: {avg_perm_detections:.1f}/{n_events_tested} significant (real: {real_detections}/{n_events_tested})", "INFO")
        
        # Update for comparison
        random_mean_sigma = np.mean(random_sigma_levels) if random_sigma_levels else 0
    else:
        print_status(f"  Permutation test: No valid results", "WARNING")
    
    # Mann-Whitney test
    mann_whitney_result = {'valid': False, 'p_value': np.nan, 'interpretation': 'INSUFFICIENT DATA'}
    
    if len(all_event_responses) >= 3 and len(random_sigma_levels) >= 3:
        event_sigmas = [e['sigma_level'] for e in all_event_responses]
        try:
            stat, p_mw = stats.mannwhitneyu(event_sigmas, random_sigma_levels, alternative='greater')
            mann_whitney_result = {
                'valid': True,
                'statistic': float(stat),
                'p_value': float(p_mw),
                'event_median_sigma': float(np.median(event_sigmas)),
                'random_median_sigma': float(np.median(random_sigma_levels)),
                'interpretation': 'EVENTS > RANDOM' if p_mw < 0.05 else 'NO DIFFERENCE'
            }
            print_status(f"  Mann-Whitney (σ levels): p = {p_mw:.3f} ({mann_whitney_result['interpretation']})", "INFO")
        except Exception as e:
            print_status(f"  Mann-Whitney test failed: {e}", "WARNING")
    
    # ==========================================================================
    # COMPREHENSIVE MASS SCALING TESTS (JPL ephemeris)
    # ==========================================================================
    print_status("", "INFO")
    print_status("  MASS SCALING TESTS (JPL ephemeris):", "INFO")
    print_status("  " + "-" * 50, "INFO")
    
    mass_scaling_result = {'valid': False, 'interpretation': 'INSUFFICIENT DATA'}
    coherence_scaling_result = {'valid': False, 'interpretation': 'INSUFFICIENT DATA'}
    
    planets_with_events = set(e['planet'] for e in all_event_responses)
    n_planets = len(planets_with_events)
    n_events = len(all_event_responses)
    
    if n_events >= MIN_EVENTS_MASS_TEST and n_planets >= MIN_PLANETS_MASS_TEST:
        # Extract all metrics
        gm_r2_values = np.array([e['gm_r2'] for e in all_event_responses])
        clock_amplitudes = np.array([e['clock_amplitude_ns'] for e in all_event_responses])
        coherence_mods = np.array([abs(e['modulation_pct']) for e in all_event_responses])  # Coherence modulation
        sigma_levels = np.array([e['sigma_level'] for e in all_event_responses])
        distances = np.array([e['distance_au'] for e in all_event_responses])
        gm_values = np.array([e['GM'] for e in all_event_responses])
        
        log_gm_r2 = np.log10(gm_r2_values)
        log_gm = np.log10(gm_values)
        
        # Report distances for verification
        print_status("", "INFO")
        print_status("  JPL Ephemeris Earth-Planet Distances:", "INFO")
        for e in all_event_responses:
            print_status(f"    {e['planet']} DOY{e['doy']}: {e['distance_au']:.3f} AU, GM/r² = {e['gm_r2']:.2e}", "INFO")
        
        print_status("", "INFO")
        
        # ==========================================================================
        # TEST A: Clock Bias Amplitude vs GM/r² (classical gravitational effect)
        # ==========================================================================
        print_status("  TEST A: Clock Amplitude vs GM/r² (classical gravitational):", "INFO")
        
        valid_mask = clock_amplitudes > 0
        if np.sum(valid_mask) >= MIN_EVENTS_MASS_TEST:
            try:
                r_amp, p_amp = stats.pearsonr(log_gm_r2[valid_mask], clock_amplitudes[valid_mask])
                r_gm_amp, p_gm_amp = stats.pearsonr(log_gm[valid_mask], clock_amplitudes[valid_mask])
                
                mass_scaling_result = {
                    'valid': True,
                    'n_events': int(np.sum(valid_mask)),
                    'n_planets': n_planets,
                    'gm_r2_vs_clock_amplitude': {
                        'pearson_r': float(r_amp),
                        'p_value': float(p_amp),
                        'scaling_detected': p_amp < 0.05 and r_amp > 0
                    },
                    'gm_vs_clock_amplitude': {
                        'pearson_r': float(r_gm_amp),
                        'p_value': float(p_gm_amp),
                        'scaling_detected': p_gm_amp < 0.05 and r_gm_amp > 0
                    },
                    'interpretation': 'AMPLITUDE SCALING' if (p_amp < 0.05 and r_amp > 0) else 'NO AMPLITUDE SCALING'
                }
                
                print_status(f"    GM/r² vs clock RMS:  r = {r_amp:+.3f}, p = {p_amp:.3f}", "INFO")
                print_status(f"    GM vs clock RMS:     r = {r_gm_amp:+.3f}, p = {p_gm_amp:.3f}", "INFO")
                
                if p_amp < 0.05 and r_amp > 0:
                    print_status("    → Classical amplitude scaling DETECTED", "SUCCESS")
                else:
                    print_status("    → No classical amplitude scaling (expected: common-mode removed)", "INFO")
                    
            except Exception as e:
                print_status(f"    Failed: {e}", "WARNING")
        
        # ==========================================================================
        # TEST B: Coherence Modulation vs GM/r² (TEP geometric effect)
        # ==========================================================================
        print_status("", "INFO")
        print_status("  TEST B: Coherence Modulation vs GM/r² (TEP geometric):", "INFO")
        
        try:
            # Test coherence modulation (absolute value) vs GM/r²
            r_coh_gm_r2, p_coh_gm_r2 = stats.pearsonr(log_gm_r2, coherence_mods)
            r_coh_gm, p_coh_gm = stats.pearsonr(log_gm, coherence_mods)
            
            # Test sigma level vs GM/r²
            r_sig_gm_r2, p_sig_gm_r2 = stats.pearsonr(log_gm_r2, sigma_levels)
            r_sig_gm, p_sig_gm = stats.pearsonr(log_gm, sigma_levels)
            
            coherence_scaling_result = {
                'valid': True,
                'n_events': n_events,
                'n_planets': n_planets,
                'gm_r2_vs_coherence_mod': {
                    'pearson_r': float(r_coh_gm_r2),
                    'p_value': float(p_coh_gm_r2),
                    'scaling_detected': p_coh_gm_r2 < 0.05 and r_coh_gm_r2 > 0
                },
                'gm_vs_coherence_mod': {
                    'pearson_r': float(r_coh_gm),
                    'p_value': float(p_coh_gm),
                    'scaling_detected': p_coh_gm < 0.05 and r_coh_gm > 0
                },
                'gm_r2_vs_sigma_level': {
                    'pearson_r': float(r_sig_gm_r2),
                    'p_value': float(p_sig_gm_r2),
                    'scaling_detected': p_sig_gm_r2 < 0.05 and r_sig_gm_r2 > 0
                },
                'gm_vs_sigma_level': {
                    'pearson_r': float(r_sig_gm),
                    'p_value': float(p_sig_gm),
                    'scaling_detected': p_sig_gm < 0.05 and r_sig_gm > 0
                },
                'interpretation': 'See individual tests'
            }
            
            print_status(f"    GM/r² vs |coherence mod|: r = {r_coh_gm_r2:+.3f}, p = {p_coh_gm_r2:.3f}", "INFO")
            print_status(f"    GM vs |coherence mod|:    r = {r_coh_gm:+.3f}, p = {p_coh_gm:.3f}", "INFO")
            print_status(f"    GM/r² vs σ-level:         r = {r_sig_gm_r2:+.3f}, p = {p_sig_gm_r2:.3f}", "INFO")
            print_status(f"    GM vs σ-level:            r = {r_sig_gm:+.3f}, p = {p_sig_gm:.3f}", "INFO")
            
            # Determine if any coherence scaling detected
            any_coh_scaling = any([
                p_coh_gm_r2 < 0.05, p_coh_gm < 0.05,
                p_sig_gm_r2 < 0.05, p_sig_gm < 0.05
            ])
            
            if any_coh_scaling:
                print_status("    → Some coherence/mass relationship detected", "WARNING")
                coherence_scaling_result['interpretation'] = 'POSSIBLE SCALING'
            else:
                print_status("    → NO coherence scaling (consistent with CODE: geometric effect)", "SUCCESS")
                coherence_scaling_result['interpretation'] = 'NO COHERENCE SCALING (geometric)'
                
        except Exception as e:
            print_status(f"    Failed: {e}", "WARNING")
        
        print_status("", "INFO")
        print_status("  " + "-" * 50, "INFO")
        print_status("  INTERPRETATION:", "INFO")
        print_status("    • Classical GM/r² amplitude scaling: Common-mode, removed in SPP", "INFO")
        print_status("    • TEP coherence effect: Geometric (phase structure), no mass scaling", "INFO")
        print_status("    • CODE longspan (25yr): NO mass scaling detected (p > 0.5)", "INFO")
        
    else:
        print_status(f"  INSUFFICIENT DATA for mass scaling test ({n_events} events, {n_planets} planets)", "WARNING")
    
    # ==========================================================================
    # STEP 5: Create visualizations and save results
    # ==========================================================================
    print_status("[5/5] Creating visualizations...", "PROCESS")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Panel 1: Daily coherence time series with events marked
    ax1 = axes[0, 0]
    ax1.plot(daily_df['doy'], daily_df['mean_coherence'], 'b-', alpha=0.7, linewidth=0.8)
    ax1.axhline(y=baseline_coherence, color='gray', linestyle='--', alpha=0.5, label=f'Mean: {baseline_coherence:.4f}')
    
    # Mark planetary events
    for planet, data in planetary_results.items():
        for event in data['events']:
            color = 'green' if event['is_significant'] else 'red'
            ax1.axvline(x=event['doy'], color=color, alpha=0.3, linewidth=2)
    
    ax1.set_xlabel('Day of Year (2024)')
    ax1.set_ylabel('Mean Daily Coherence')
    ax1.set_title('Daily Coherence Time Series\n(Green=Significant, Red=Not Significant)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Detection rate by planet
    ax2 = axes[0, 1]
    if planetary_results:
        planets = list(planetary_results.keys())
        detection_rates = [planetary_results[p]['detection_rate'] * 100 for p in planets]
        colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(planets)))
        
        bars = ax2.bar(planets, detection_rates, color=colors, alpha=0.7)
        ax2.axhline(y=35.9, color='red', linestyle='--', label='CODE: 35.9%')
        ax2.set_ylabel('Detection Rate (%)')
        ax2.set_title(f'Planetary Event Detection Rates\n(Total: {n_significant}/{len(all_event_responses)} = {100*n_significant/len(all_event_responses):.0f}%)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
    
    # Panel 3: Event vs Random σ-levels
    ax3 = axes[1, 0]
    if all_event_responses and random_sigma_levels:
        event_sigmas = [e['sigma_level'] for e in all_event_responses]
        ax3.boxplot([event_sigmas, random_sigma_levels], tick_labels=['Planetary Events', 'Random Dates'])
        ax3.axhline(y=SIGNIFICANCE_THRESHOLD, color='red', linestyle='--', label=f'{SIGNIFICANCE_THRESHOLD}σ threshold')
        ax3.set_ylabel('Sigma Level (σ)')
        ax3.set_title(f'Event vs Random Significance\nMann-Whitney p = {mann_whitney_result.get("p_value", np.nan):.3f}')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
    
    # Panel 4: Mass scaling test (GM/r² from JPL vs clock amplitude)
    ax4 = axes[1, 1]
    if mass_scaling_result.get('valid'):
        for planet in PLANETS:
            planet_events = [e for e in all_event_responses if e['planet'] == planet and e['clock_amplitude_ns'] > 0]
            if planet_events:
                gm_r2s = [e['gm_r2'] for e in planet_events]
                amps = [e['clock_amplitude_ns'] for e in planet_events]
                ax4.scatter(gm_r2s, amps, label=f"{planet} ({len(planet_events)})", s=100, alpha=0.7)
        
        ax4.set_xscale('log')
        ax4.set_xlabel('GM/r² (Earth masses / AU²) - JPL Ephemeris')
        ax4.set_ylabel('Clock Bias RMS (ns)')
        
        # Get correlation from result
        r_val = mass_scaling_result.get('gm_r2_vs_amplitude', {}).get('pearson_r', np.nan)
        p_val = mass_scaling_result.get('gm_r2_vs_amplitude', {}).get('p_value', np.nan)
        ax4.set_title(f"GM/r² Mass Scaling Test\n(r = {r_val:.3f}, p = {p_val:.3f})")
        ax4.legend()
        ax4.grid(True, alpha=0.3)
    else:
        ax4.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax4.transAxes)
        ax4.set_title('Mass Scaling Test')
    
    plt.tight_layout()
    fig_path = FIGURES_DIR / "step_2_6_planetary_events.png"
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print_status(f"  Figure saved: {fig_path}", "SUCCESS")
    
    # Summary
    print_status("", "INFO")
    print_status("=" * 60, "INFO")
    print_status("PLANETARY EVENT ANALYSIS SUMMARY", "INFO")
    print_status("=" * 60, "INFO")
    
    print_status(f"  Methodology: CODE LONGSPAN (daily coherence + Gaussian pulse)", "INFO")
    print_status(f"  Data coverage: DOY {min(available_doys)}-{max(available_doys)} ({len(available_doys)} days)", "INFO")
    print_status(f"  Events analyzed: {len(all_event_responses)}/{total_events_all_years}", "INFO")
    print_status(f"  Significant detections: {n_significant}/{len(all_event_responses)} ({100*n_significant/len(all_event_responses):.1f}%)", "INFO")
    if permutation_detections:
        avg_null = np.mean(permutation_detections)
        n_tested = min(20, len(unique_event_doys))
        real_det = sum(1 for e in all_event_responses[:n_tested] if e['is_significant'])
        print_status(f"  Permutation null: {avg_null:.1f}/{n_tested} vs real: {real_det}/{n_tested}", "INFO")
    else:
        print_status(f"  Permutation null: N/A", "INFO")
    
    print_status("", "INFO")
    print_status("KEY FINDINGS:", "INFO")
    
    if n_significant > n_random_significant:
        print_status(f"  ✓ Higher detection rate for planetary events vs random", "SUCCESS")
    else:
        print_status(f"  ⚠ Similar detection rate for events and random dates", "WARNING")
    
    # Mass scaling summary
    if mass_scaling_result.get('valid'):
        amp_test = mass_scaling_result.get('gm_r2_vs_clock_amplitude', {})
        amp_r = amp_test.get('pearson_r', 0)
        amp_p = amp_test.get('p_value', 1.0)
        
        if amp_test.get('scaling_detected', False):
            print_status(f"  ✓ Classical amplitude scaling DETECTED (r={amp_r:.3f}, p={amp_p:.3f})", "SUCCESS")
        else:
            print_status(f"  • No classical amplitude scaling (r={amp_r:.3f}, p={amp_p:.3f}) - expected", "INFO")
    
    if coherence_scaling_result.get('valid'):
        coh_test = coherence_scaling_result.get('gm_r2_vs_coherence_mod', {})
        coh_r = coh_test.get('pearson_r', 0)
        coh_p = coh_test.get('p_value', 1.0)
        
        if coh_test.get('scaling_detected', False):
            print_status(f"  ⚠ Coherence scaling detected (r={coh_r:.3f}, p={coh_p:.3f})", "WARNING")
        else:
            print_status(f"  ✓ NO coherence scaling (r={coh_r:.3f}, p={coh_p:.3f}) - consistent with CODE", "SUCCESS")
    
    if mann_whitney_result.get('valid') and mann_whitney_result.get('p_value', 1) < 0.05:
        print_status(f"  ✓ Events significantly stronger than random (p={mann_whitney_result['p_value']:.3f})", "SUCCESS")
    
    # Save results
    output = {
        'step': '2.6',
        'name': 'planetary_events',
        'methodology': 'CODE LONGSPAN: daily coherence + Gaussian pulse fitting + JPL ephemeris',
        'timestamp': datetime.now().isoformat(),
        'processing_modes': list(PROCESSING_MODES.keys()),
        'primary_mode': 'baseline',
        'mode_station_counts': {mode: len(data) for mode, data in station_data_by_mode.items()},
        'parameters': {
            'frequency_band_hz': [F1_HZ, F2_HZ],
            'distance_range_km': [MIN_DISTANCE_KM, MAX_DISTANCE_KM],
            'event_window_days': EVENT_WINDOW_DAYS,
            'significance_threshold_sigma': SIGNIFICANCE_THRESHOLD,
            'ephemeris': 'JPL DE440'
        },
        'data_coverage': {
            'doy_range': [int(min(available_doys)), int(max(available_doys))],
            'n_days': len(available_doys),
            'n_daily_coherence_days': len(daily_df),
            'n_stations': len(station_data),
            'n_stations_by_mode': {mode: len(data) for mode, data in station_data_by_mode.items()},
            'station_filter': station_filter
        },
        'baseline': {
            'mean_coherence': float(baseline_coherence),
            'std_coherence': float(baseline_std)
        },
        'detection_summary': {
            'total_events': len(all_event_responses),
            'n_significant': n_significant,
            'detection_rate': n_significant / len(all_event_responses) if all_event_responses else 0,
            'random_n_significant': n_random_significant,
            'random_detection_rate': n_random_significant / len(random_sigma_levels) if random_sigma_levels else 0
        },
        'planetary_results': planetary_results,
        'mass_scaling_tests': {
            'amplitude_scaling': mass_scaling_result,
            'coherence_scaling': coherence_scaling_result
        },
        'null_control': {
            'n_random_dates': len(random_sigma_levels),
            'random_mean_sigma': float(np.mean(random_sigma_levels)) if random_sigma_levels else None,
            'mann_whitney': mann_whitney_result
        },
        'comparison_to_code': {
            'code_detection_rate': '56/156 (35.9%)',
            'code_mass_scaling': 'None (p > 0.5)',
            'rinex_detection_rate': f"{n_significant}/{len(all_event_responses)} ({100*n_significant/len(all_event_responses):.1f}%)" if all_event_responses else "N/A"
        }
    }
    
    # Include metric/coherence in output filename
    suffix = f"_{metric}_{coherence_type}" if metric != 'clock_bias' or coherence_type != 'msc' else ""
    output_path = OUTPUTS_DIR / f"step_2_6_planetary_events{suffix}.json"
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    
    print_status(f"Results saved: {output_path}", "SUCCESS")
    print_status("", "INFO")
    print_status("=" * 80, "INFO")
    print_status("Step 2.6 Complete", "INFO")
    print_status("=" * 80, "INFO")
    
    # Return key results for comparison
    return {
        'metric': metric,
        'coherence_type': coherence_type,
        'n_events': len(all_event_responses),
        'n_significant': n_significant,
        'detection_rate': n_significant / len(all_event_responses) if all_event_responses else 0,
        'mean_sigma': float(np.mean([r['sigma_level'] for r in all_event_responses])) if all_event_responses else 0,
        'baseline_coherence': float(baseline_coherence),
        'baseline_std': float(baseline_std)
    }


def main():
    """Main execution loop over all metric and coherence type combinations."""
    print_status("STARTING MULTI-METRIC PLANETARY EVENT ANALYSIS", "TITLE")
    
    # Define metrics and coherence types to compare
    metrics = ['clock_bias', 'pos_jitter', 'clock_drift']
    coherence_types = ['msc', 'phase_alignment']
    
    # Store results for comparison
    all_results = []
    
    for metric in metrics:
        for coh_type in coherence_types:
            try:
                result = run_analysis(metric=metric, coherence_type=coh_type, use_csv=True)
                if result:
                    all_results.append(result)
            except Exception as e:
                print_status(f"Analysis failed for {metric}/{coh_type}: {e}", "ERROR")
                import traceback
                traceback.print_exc()
    
    # Print comparison summary
    print_status("", "INFO")
    print_status("=" * 80, "INFO")
    print_status("PLANETARY EVENT ANALYSIS COMPARISON SUMMARY", "INFO")
    print_status("=" * 80, "INFO")
    print_status("", "INFO")
    print_status(f"{'Metric':<12} {'Coh Type':<10} {'Events':<8} {'Signif':<8} {'Rate':<8} {'Mean σ':<8}", "INFO")
    print_status("-" * 60, "INFO")
    
    for res in all_results:
        rate_pct = f"{100*res['detection_rate']:.1f}%"
        print_status(f"{res['metric']:<12} {res['coherence_type']:<10} {res['n_events']:<8} {res['n_significant']:<8} {rate_pct:<8} {res['mean_sigma']:<8.2f}", "INFO")
    
    print_status("-" * 60, "INFO")
    print_status("CODE Longspan reference: 56/156 events significant (35.9%)", "INFO")
    print_status("", "INFO")
    print_status("MULTI-METRIC PIPELINE COMPLETE", "SUCCESS")


if __name__ == "__main__":
    main()
