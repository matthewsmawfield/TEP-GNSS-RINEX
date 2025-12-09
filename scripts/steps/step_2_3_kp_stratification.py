#!/usr/bin/env python3
"""
TEP-GNSS-RINEX - STEP 2.3a: Kp Geomagnetic Stratification (Memory Optimized)
============================================================================
Tests whether TEP signal differs between geomagnetically quiet (Kp<3) and
storm (Kp≥3) conditions.

Memory Optimized: Processes metrics sequentially.

Author: Matthew Lukin Smawfield
Date: 9 December 2025
License: CC-BY-4.0
"""

import sys
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import pandas as pd
import json
import os
from scipy.optimize import curve_fit
from datetime import datetime
import gc
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count

from scripts.utils.logger import TEPLogger, set_step_logger, print_status

# Initialize logger
logger = TEPLogger(
    name="step_2_3_kp_stratification",
    level="INFO",
    log_file_path=PROJECT_ROOT / "logs" / "step_2_3_kp_stratification.log"
)
set_step_logger(logger)

# Parallel processing
N_WORKERS = min(cpu_count(), 16)

# ============================================================================
# ANALYSIS PARAMETERS
# ============================================================================
MIN_DISTANCE_KM = 50
MAX_DISTANCE_KM = 13000
N_BINS = 40
MIN_BIN_COUNT = 50
KP_THRESHOLD = 3.0

RESULTS_DIR = PROJECT_ROOT / "results" / "outputs"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"

METRICS = ['clock_bias', 'pos_jitter', 'clock_drift']
COHERENCE_TYPES = ['msc', 'phase_alignment']

# Station filters - consistent with other Step 2.x scripts
STATION_FILTERS = [
    ('ALL_STATIONS', 'none'),
    ('OPTIMAL_100', 'optimal_100_metadata.json'),
    ('DYNAMIC_50', 'dynamic_50_metadata.json')
]

# Global station filter set
GOOD_STATIONS_SET = None

# ALL THREE PROCESSING MODES - consistent with other scripts
PROCESSING_MODES = {
    'baseline': {
        'csv_file': 'step_2_0_pairs_baseline.csv',
        'description': 'GPS L1 only'
    },
    'ionofree': {
        'csv_file': 'step_2_0_pairs_ionofree.csv',
        'description': 'Dual-frequency L1+L2 ionosphere-free'
    },
    'multi_gnss': {
        'csv_file': 'step_2_0_pairs_multi_gnss.csv',
        'description': 'GPS+GLO+GAL+BDS'
    }
}

def exp_decay(r, A, lam, C0):
    return A * np.exp(-r / lam) + C0

def fit_exponential(distances, coherences):
    if len(distances) < 1000: return None
    bin_edges = np.logspace(np.log10(MIN_DISTANCE_KM), np.log10(MAX_DISTANCE_KM), N_BINS + 1)
    bin_centers, bin_means, bin_counts = [], [], []
    
    distances = np.asarray(distances)
    coherences = np.asarray(coherences)
    
    for i in range(N_BINS):
        mask = (distances >= bin_edges[i]) & (distances < bin_edges[i+1])
        count = np.sum(mask)
        if count >= MIN_BIN_COUNT:
            bin_centers.append((bin_edges[i] + bin_edges[i+1]) / 2)
            bin_means.append(np.nanmean(coherences[mask]))
            bin_counts.append(count)
    if len(bin_centers) < 5: return None
    
    bin_centers, bin_means, bin_counts = np.array(bin_centers), np.array(bin_means), np.array(bin_counts)
    try:
        popt, pcov = curve_fit(exp_decay, bin_centers, bin_means, p0=[0.5, 2000, 0],
                               sigma=1.0/np.sqrt(bin_counts), bounds=([0, 100, -1], [2, 20000, 1]), maxfev=10000)
        # Weighted R² (consistent with weighted fit)
        weights = bin_counts.astype(float)
        predicted = exp_decay(bin_centers, *popt)
        weighted_mean = np.average(bin_means, weights=weights)
        ss_res = np.sum(weights * (bin_means - predicted)**2)
        ss_tot = np.sum(weights * (bin_means - weighted_mean)**2)
        r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0
        return {'lambda_km': float(popt[1]), 'r_squared': float(r2), 'n_pairs': int(sum(bin_counts))}
    except:
        return None

def load_kp_index():
    """
    Load REAL Kp index data from GFZ Potsdam (1932-present).
    Source: https://www-app3.gfz-potsdam.de/kp_index/Kp_ap_since_1932.txt
    """
    import urllib.request
    import ssl
    from datetime import datetime as dt
    
    url = "https://www-app3.gfz-potsdam.de/kp_index/Kp_ap_since_1932.txt"
    print_status(f"  Downloading real Kp data from GFZ Potsdam...")
    
    try:
        ssl_context = ssl.create_default_context()
        with urllib.request.urlopen(url, context=ssl_context, timeout=60) as response:
            content = response.read().decode('utf-8').splitlines()
    except Exception as e:
        raise RuntimeError(f"CRITICAL: Failed to download real Kp data from GFZ Potsdam: {e}. "
                         f"Cannot proceed with synthetic data. Check network connection.")
    
    # Parse the data - aggregate to daily mean Kp
    # Format: YYYY MM DD hh.h hh._m days days_m Kp ap D
    daily_kp_values = {}  # {(year, doy): [kp1, kp2, ...]}
    
    for line in content:
        if line.startswith('#') or not line.strip():
            continue
        try:
            parts = line.split()
            if len(parts) < 8:
                continue
            year = int(parts[0])
            month = int(parts[1])
            day = int(parts[2])
            kp_val = float(parts[7])
            
            # Only 2022-2024
            if year < 2022 or year > 2024:
                continue
            
            # Convert to DOY
            date_obj = dt(year, month, day)
            doy = date_obj.timetuple().tm_yday
            key = (year, doy)
            
            if key not in daily_kp_values:
                daily_kp_values[key] = []
            daily_kp_values[key].append(kp_val)
            
        except (ValueError, IndexError):
            continue
    
    # Compute daily mean Kp
    kp_data = {}
    for (year, doy), values in daily_kp_values.items():
        key = f"{year}_{doy:03d}"
        kp_data[key] = np.mean(values)
    
    print_status(f"  Loaded {len(kp_data)} days of REAL Kp data", level="SUCCESS")
    
    # Verify data quality
    kp_values = list(kp_data.values())
    print_status(f"  Kp range: {min(kp_values):.2f} - {max(kp_values):.2f}, mean: {np.mean(kp_values):.2f}")
    
    return kp_data

def process_metric(csv_file, target_metric, quiet_days, storm_days, station_set=None):
    print_status(f"\nProcessing metric: {target_metric}...")
    
    # Store only dist/coh for quiet/storm/all per coherence type
    # actually, all = quiet + storm, so we just store quiet and storm
    data = {c: {'quiet': {'dist': [], 'coh': []}, 'storm': {'dist': [], 'coh': []}} for c in COHERENCE_TYPES}
    
    chunk_size = 10_000_000  # Larger chunks for high-RAM systems
    for chunk in pd.read_csv(csv_file, chunksize=chunk_size):
        if 'metric' in chunk.columns:
            # Handle prefixed metric names (ionofree_*, multi_gnss_*)
            mc = chunk[chunk['metric'].str.endswith(target_metric)]
        else: continue
        if mc.empty: continue
        
        # Apply station filter if provided
        if station_set is not None:
            # Vectorized check: station1 IN set AND station2 IN set
            # This is slow with pandas .apply or list comp.
            # Faster: use .isin()
            s1_mask = mc['station1'].isin(station_set)
            s2_mask = mc['station2'].isin(station_set)
            mc = mc[s1_mask & s2_mask]
            if mc.empty: continue

        has_phase = 'phase_alignment' in mc.columns
        years = mc['year'].values
        doys = mc['doy'].values
        dist = mc['distance_km'].values
        msc = mc['coherence'].values
        phase = mc['phase_alignment'].values if has_phase else None
        
        # Vectorized mask
        # Combining year+doy to check efficiently?
        # A bit slow to loop. Let's make a set of (y,d).
        # Optimization: Map (y,d) to status 0=none, 1=quiet, 2=storm
        # ... assuming helper function done elsewhere or passed in.
        
        # Simple loop for now (vectorized is hard with set lookups)
        # Or build a map array if year/doy range is small
        # Actually, let's just loop over rows, it's Python but maybe fast enough for 5M chunk?
        # No, 5M loop is slow.
        # Efficient way:
        # Create a lookup array for [year][doy] -> status
        # status = lookup[years, doys]
        
        # Assuming quiet_days is a set of (year, doy)
        # We can make a pandas index or map
        
        # Fast way:
        # q_mask = [(y, d) in quiet_days for y, d in zip(years, doys)] is slow
        # Better: use lookup table
        pass # implemented in main loop logic below for speed
        
        # Re-implementing logic from original script but inside loop
        quiet_mask = np.array([(y, d) in quiet_days for y, d in zip(years, doys)])
        storm_mask = np.array([(y, d) in storm_days for y, d in zip(years, doys)])
        
        # MSC
        valid = ~np.isnan(msc)
        
        q_idx = valid & quiet_mask
        s_idx = valid & storm_mask
        
        if np.any(q_idx):
            data['msc']['quiet']['dist'].extend(dist[q_idx].tolist())
            data['msc']['quiet']['coh'].extend(msc[q_idx].tolist())
        if np.any(s_idx):
            data['msc']['storm']['dist'].extend(dist[s_idx].tolist())
            data['msc']['storm']['coh'].extend(msc[s_idx].tolist())
            
        # Phase
        if phase is not None:
            valid = ~np.isnan(phase)
            q_idx = valid & quiet_mask
            s_idx = valid & storm_mask
            if np.any(q_idx):
                data['phase_alignment']['quiet']['dist'].extend(dist[q_idx].tolist())
                data['phase_alignment']['quiet']['coh'].extend(phase[q_idx].tolist())
            if np.any(s_idx):
                data['phase_alignment']['storm']['dist'].extend(dist[s_idx].tolist())
                data['phase_alignment']['storm']['coh'].extend(phase[s_idx].tolist())

    results = {}
    for c in COHERENCE_TYPES:
        mode = f"{target_metric}/{c}"
        results[mode] = {}
        
        # Quiet
        qd = data[c]['quiet']
        q_res = fit_exponential(qd['dist'], qd['coh'])
        results[mode]['quiet'] = q_res
        
        # Storm
        sd = data[c]['storm']
        s_res = fit_exponential(sd['dist'], sd['coh'])
        results[mode]['storm'] = s_res
        
        # All (combine)
        all_dist = np.array(qd['dist'] + sd['dist'])
        all_coh = np.array(qd['coh'] + sd['coh'])
        all_res = fit_exponential(all_dist, all_coh)
        results[mode]['all'] = all_res
        
        print_status(f"  {mode}: Q={q_res['lambda_km'] if q_res else 'N/A'}, S={s_res['lambda_km'] if s_res else 'N/A'}")
        
    return results


def load_station_filter(filter_config):
    """Load station filter and return set of valid stations."""
    global GOOD_STATIONS_SET
    
    if filter_config == 'none':
        GOOD_STATIONS_SET = None
        return None, {'type': 'all'}
    
    # Load from metadata file (optimal_100_metadata.json or dynamic_50_metadata.json)
    filter_file = PROCESSED_DIR / filter_config
    if filter_file.exists():
        with open(filter_file) as f:
            meta = json.load(f)
        station_set = set(meta.get('stations', {}).keys())
        GOOD_STATIONS_SET = station_set
        return station_set, {'type': 'metadata', 'n_stations': len(station_set)}
    
    GOOD_STATIONS_SET = None
    return None, {'type': 'all'}


def main():
    import argparse
    parser = argparse.ArgumentParser(description='TEP-GNSS-RINEX Step 2.3: Kp Stratification')
    parser.add_argument('--filter', type=str, default='all',
                        help='Filter: "all", "none", "optimal_100_metadata.json", or "dynamic_50_metadata.json"')
    args = parser.parse_args()
    
    print_status('STEP 2.3a: Kp Geomagnetic Stratification (COMPREHENSIVE)', level="TITLE")
    print_status('Processing: 3 filters x 3 modes x 3 metrics x 2 coherence types')
    
    # Select filters
    if args.filter == 'all':
        filters = STATION_FILTERS
    else:
        matching = [f for f in STATION_FILTERS if f[1] == args.filter]
        filters = matching if matching else [('ALL_STATIONS', 'none')]
    
    # Load Kp data ONCE
    kp_data = load_kp_index()
    quiet_days = set()
    storm_days = set()
    for key, kp in kp_data.items():
        year, doy = map(int, key.split('_'))
        if kp < KP_THRESHOLD: quiet_days.add((year, doy))
        else: storm_days.add((year, doy))
        
    print_status(f"  Quiet days: {len(quiet_days)} ({100*len(quiet_days)/(len(quiet_days)+len(storm_days)):.1f}%)")
    print_status(f"  Storm days: {len(storm_days)} ({100*len(storm_days)/(len(quiet_days)+len(storm_days)):.1f}%)")
    
    all_filter_results = {}
    
    for filter_name, filter_config in filters:
        filter_key = filter_name.lower()
        print_status(f"\n{'#'*80}")
        print_status(f"# FILTER: {filter_name}")
        print_status(f"{'#'*80}")
        
        station_set, filter_meta = load_station_filter(filter_config)
        if station_set:
            print_status(f"  Using {len(station_set)} stations")
        
        all_mode_results = {}
        
        for mode_name, mode_config in PROCESSING_MODES.items():
            csv_file = RESULTS_DIR / mode_config['csv_file']
            
            print_status(f"\n  MODE: {mode_name.upper()} ({mode_config['description']})")
            
            if not csv_file.exists():
                print_status(f"    WARNING: {csv_file.name} not found, skipping", level="WARNING")
                continue
            
            mode_results = {}
            for metric in METRICS:
                res = process_metric(csv_file, metric, quiet_days, storm_days, station_set=station_set)
                mode_results.update(res)
                gc.collect()
            
            all_mode_results[mode_name] = mode_results
            
            # Print mode summary
            for key, val in mode_results.items():
                if val.get('quiet') and val.get('storm'):
                    q_lam = val['quiet']['lambda_km']
                    s_lam = val['storm']['lambda_km']
                    delta = 100 * (s_lam - q_lam) / q_lam
                    print_status(f"      {key}: Q={q_lam:.0f}km, S={s_lam:.0f}km, Δ={delta:+.1f}%")
        
        all_filter_results[filter_key] = {
            'filter_name': filter_name,
            'filter_metadata': filter_meta,
            'results_by_mode': all_mode_results
        }
        
        # Save per-filter JSON
        output = {
            'step': '2.3a', 
            'name': 'kp_stratification',
            'filter': filter_name,
            'timestamp': datetime.now().isoformat(),
            'kp_threshold': KP_THRESHOLD,
            'quiet_days': len(quiet_days),
            'storm_days': len(storm_days),
            'results_by_mode': all_mode_results
        }
        with open(RESULTS_DIR / f"step_2_3_kp_stratification_{filter_key}.json", 'w') as f:
            json.dump(output, f, indent=2, default=str)
        print_status(f"  Saved: step_2_3_kp_stratification_{filter_key}.json", level="SUCCESS")
    
    # Save summary
    summary = {
        'step': '2.3a', 
        'name': 'kp_stratification_comprehensive',
        'timestamp': datetime.now().isoformat(),
        'kp_threshold': KP_THRESHOLD,
        'filters_tested': [f[0] for f in filters],
        'modes_tested': list(PROCESSING_MODES.keys()),
        'results_by_filter': all_filter_results
    }
    with open(RESULTS_DIR / "step_2_3_kp_stratification_summary.json", 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    print_status(f"COMPLETE: {len(filters)} filters × {len(PROCESSING_MODES)} modes", level="SUCCESS")

if __name__ == '__main__':
    main()
