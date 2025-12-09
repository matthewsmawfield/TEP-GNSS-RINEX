#!/usr/bin/env python3
"""
TEP-GNSS-RINEX - STEP 2.4a: Seasonal Analysis (Memory Optimized)
================================================================
Tests whether TEP signal varies with season (proxy for diurnal/ionospheric).

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
    name="step_2_4_diurnal_analysis",
    level="INFO",
    log_file_path=PROJECT_ROOT / "logs" / "step_2_4_diurnal_analysis.log"
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

GOOD_STATIONS_SET = None

# ALL THREE PROCESSING MODES
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

SEASONS = {
    'winter': list(range(1, 80)) + list(range(356, 366)),
    'spring': list(range(80, 172)),
    'summer': list(range(172, 264)),
    'autumn': list(range(264, 356)),
}

def get_season_map():
    # Build array map 1-366 -> season index/string
    s_map = {}
    for s, doys in SEASONS.items():
        for d in doys: s_map[d] = s
    return s_map

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


def process_metric(csv_file, target_metric, season_map):
    print_status(f"\nProcessing metric: {target_metric}...")
    
    # Store dist/coh per season
    data = {c: {s: {'dist': [], 'coh': []} for s in SEASONS} for c in COHERENCE_TYPES}
    
    chunk_size = 10_000_000  # Larger chunks for high-RAM systems
    for chunk in pd.read_csv(csv_file, chunksize=chunk_size):
        # Apply station filter if set
        if GOOD_STATIONS_SET:
            chunk = chunk[chunk['station1'].isin(GOOD_STATIONS_SET) & chunk['station2'].isin(GOOD_STATIONS_SET)]
        if 'metric' in chunk.columns:
            # Handle prefixed metric names (ionofree_*, multi_gnss_*)
            mc = chunk[chunk['metric'].str.endswith(target_metric)]
        else: continue
        if mc.empty: continue
            
        has_phase = 'phase_alignment' in mc.columns
        doys = mc['doy'].values
        dist = mc['distance_km'].values
        msc = mc['coherence'].values
        phase = mc['phase_alignment'].values if has_phase else None
        
        # Map DOY to season efficiently
        # Vectorized map? Pandas map is slow.
        # Use simple list comprehension for now or numpy map
        # season_arr = np.array([season_map.get(d, 'winter') for d in doys])
        # This is slow. Faster:
        # Create a lookup array for 1-366.
        # But for now, loop is safer for 5M chunk than complexity.
        
        # Actually, let's optimize.
        # We can iterate by season mask?
        # SEASONS dict has lists of DOYs.
        # mask_winter = np.isin(doys, SEASONS['winter']) ...
        # This is fast!
        
        for s_name, s_doys in SEASONS.items():
            mask = np.isin(doys, s_doys)
            if not np.any(mask): continue
            
            # MSC
            valid = mask & ~np.isnan(msc)
            if np.any(valid):
                data['msc'][s_name]['dist'].extend(dist[valid].tolist())
                data['msc'][s_name]['coh'].extend(msc[valid].tolist())
                
            # Phase
            if phase is not None:
                valid = mask & ~np.isnan(phase)
                if np.any(valid):
                    data['phase_alignment'][s_name]['dist'].extend(dist[valid].tolist())
                    data['phase_alignment'][s_name]['coh'].extend(phase[valid].tolist())

    results = {}
    for c in COHERENCE_TYPES:
        mode = f"{target_metric}/{c}"
        results[mode] = {}
        for s in SEASONS:
            d = data[c][s]
            res = fit_exponential(d['dist'], d['coh'])
            results[mode][s] = res
            print_status(f"  {mode} {s}: {res['lambda_km'] if res else 'N/A'}")
            
    return results

def main():
    import argparse
    parser = argparse.ArgumentParser(description='TEP-GNSS-RINEX Step 2.4a: Seasonal Analysis')
    parser.add_argument('--filter', type=str, default='all',
                        help='Filter: "all", "none", "optimal_100_metadata.json", or "dynamic_50_metadata.json"')
    args = parser.parse_args()
    
    print_status('STEP 2.4a: Seasonal Analysis (COMPREHENSIVE)', level="TITLE")
    print_status('Processing: 3 filters x 3 modes x 3 metrics x 2 coherence types')
    
    if args.filter == 'all':
        filters = STATION_FILTERS
    else:
        matching = [f for f in STATION_FILTERS if f[1] == args.filter]
        filters = matching if matching else [('ALL_STATIONS', 'none')]
    
    s_map = get_season_map()
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
            print_status(f"\n  MODE: {mode_name.upper()}")
            
            if not csv_file.exists():
                print_status(f"    WARNING: {csv_file.name} not found", level="WARNING")
                continue
            
            mode_results = {}
            for metric in METRICS:
                res = process_metric(csv_file, metric, s_map)
                mode_results.update(res)
                gc.collect()
            
            all_mode_results[mode_name] = mode_results
        
        all_filter_results[filter_key] = {
            'filter_name': filter_name,
            'filter_metadata': filter_meta,
            'results_by_mode': all_mode_results
        }
        
        # Save per-filter JSON
        output = {
            'step': '2.4a', 'name': 'seasonal_analysis',
            'filter': filter_name,
            'timestamp': datetime.now().isoformat(),
            'results_by_mode': all_mode_results
        }
        with open(RESULTS_DIR / f"step_2_4_diurnal_analysis_{filter_key}.json", 'w') as f:
            json.dump(output, f, indent=2, default=str)
        print_status(f"  Saved: step_2_4_diurnal_analysis_{filter_key}.json", level="SUCCESS")
    
    # Save summary
    summary = {
        'step': '2.4a', 'name': 'seasonal_analysis_comprehensive',
        'timestamp': datetime.now().isoformat(),
        'filters_tested': [f[0] for f in filters],
        'modes_tested': list(PROCESSING_MODES.keys()),
        'results_by_filter': all_filter_results
    }
    with open(RESULTS_DIR / "step_2_4_diurnal_analysis_summary.json", 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    print_status(f"COMPLETE: {len(filters)} filters x {len(PROCESSING_MODES)} modes", level="SUCCESS")


if __name__ == '__main__':
    main()
