#!/usr/bin/env python3
"""
TEP-GNSS-RINEX - STEP 2.1: Regional Control Tests (Optimized)
=============================================================
Tests regional variations across 18 metric/mode combinations.
Memory optimized and consistent with Step 2.5 methodology.

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
import matplotlib.pyplot as plt
import gc
import os
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
from scripts.utils.tep_utils import bin_and_fit, ecef_to_lla, load_station_map
from scripts.utils.logger import TEPLogger, set_step_logger, print_status

# Initialize logger
logger = TEPLogger(
    name="step_2_1_control_tests",
    level="INFO",
    log_file_path=PROJECT_ROOT / "logs" / "step_2_1_control_tests.log"
)
set_step_logger(logger)

# Parallel processing settings
N_WORKERS = min(cpu_count(), 32)  # Cap at 32 to avoid memory issues

# Parameters
RESULTS_DIR = PROJECT_ROOT / "results" / "outputs"
FIGURES_DIR = PROJECT_ROOT / "results" / "figures"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"

METRICS = ['clock_bias', 'pos_jitter', 'clock_drift']
COHERENCE_TYPES = ['msc', 'phase_alignment']

# ALL THREE PROCESSING MODES - consistent across all Step 2.x scripts
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

# Default CSV file (will be updated per mode in main())
CSV_FILE = RESULTS_DIR / "step_2_0_pairs_baseline.csv"

# Station filters - consistent with other Step 2.x scripts
STATION_FILTERS = [
    ('ALL_STATIONS', 'none'),
    ('OPTIMAL_100', 'optimal_100_metadata.json'),
    ('DYNAMIC_50', 'dynamic_50_metadata.json')
]

REGIONS = ['global', 'europe', 'non_europe', 'northern', 'southern']

# Global variable for station filter
GOOD_STATIONS_SET = None

def get_station_regions(coords):
    """Map station to boolean flags."""
    # region_map: station -> (is_europe, is_northern)
    r_map = {}
    for sta, xyz in coords.items():
        lat, lon, _ = ecef_to_lla(*xyz)
        is_europe = (35 <= lat <= 72) and (-15 <= lon <= 45)
        is_northern = (lat > 0)
        r_map[sta] = (is_europe, is_northern)
    return r_map

def process_metric_worker(args):
    """Worker function for parallel processing - receives all args as tuple."""
    metric, csv_file, station_filter_set = args
    return process_metric(metric, csv_file, station_filter_set)

def process_metric(metric, csv_file=None, station_filter_set=None):
    """Process a single metric. Can be called directly or via worker."""
    print_status(f"[{os.getpid()}] Processing {metric}...")
    
    # Use passed csv_file or fall back to global
    if csv_file is None:
        csv_file = CSV_FILE
    
    data = {c: {r: {'dist': [], 'coh': []} for r in REGIONS} for c in COHERENCE_TYPES}
    
    chunk_size = 10_000_000  # Larger chunks for better I/O efficiency
    
    # Load Station Map
    coords = load_station_map(PROCESSED_DIR)
    if not coords:
        print_status("Error: Missing station coordinates.", level="ERROR")
        return None
    
    # Apply station filter if set
    if station_filter_set:
        coords = {k: v for k, v in coords.items() if k in station_filter_set}
        print_status(f"  [{os.getpid()}] Filtered to {len(coords)} stations")
    elif GOOD_STATIONS_SET:
        coords = {k: v for k, v in coords.items() if k in GOOD_STATIONS_SET}
        print_status(f"  [{os.getpid()}] Filtered to {len(coords)} stations")
    
    r_map = get_station_regions(coords)
    
    # Pre-compute station lookup for fast mapping?
    # Actually, dataframe map is fast enough.
    
    for chunk in pd.read_csv(csv_file, chunksize=chunk_size):
        if 'metric' in chunk.columns:
            # Handle prefixed metric names (ionofree_*, multi_gnss_*)
            mc = chunk[chunk['metric'].str.endswith(metric)]
        else: continue
        if mc.empty: continue
        
        # Extract columns
        s1 = mc['station1'].values
        s2 = mc['station2'].values
        dist = mc['distance_km'].values
        msc = mc['coherence'].values
        has_phase = 'phase_alignment' in mc.columns
        phase = mc['phase_alignment'].values if has_phase else None
        
        # Vectorized mapping using map/replace is slow on chunks?
        # Faster: list comprehension with dict lookup
        # (is_eur, is_north)
        flags1 = [r_map.get(s, (False, False)) for s in s1]
        flags2 = [r_map.get(s, (False, False)) for s in s2]
        
        # Unpack
        eur1 = np.array([f[0] for f in flags1])
        eur2 = np.array([f[0] for f in flags2])
        north1 = np.array([f[1] for f in flags1])
        north2 = np.array([f[1] for f in flags2])
        
        # Check if stations are in filtered set (r_map only contains filtered stations)
        in_filter1 = np.array([s in r_map for s in s1])
        in_filter2 = np.array([s in r_map for s in s2])
        in_filter = in_filter1 & in_filter2  # Both stations must be in the filtered set
        
        # Region Masks (BOTH stations must be in filter AND match region)
        # FIX: All masks now require both stations to be in the filter set
        mask_global = in_filter
        mask_europe = in_filter & eur1 & eur2
        mask_non_europe = in_filter & (~eur1) & (~eur2)
        mask_northern = in_filter & north1 & north2
        mask_southern = in_filter & (~north1) & (~north2)
        
        masks = {
            'global': mask_global,
            'europe': mask_europe,
            'non_europe': mask_non_europe,
            'northern': mask_northern,
            'southern': mask_southern
        }
        
        for coh_type in COHERENCE_TYPES:
            vals = msc if coh_type == 'msc' else phase
            if vals is None: continue
            
            valid = ~np.isnan(vals)
            
            for reg, mask in masks.items():
                final_mask = valid & mask
                if np.any(final_mask):
                    data[coh_type][reg]['dist'].extend(dist[final_mask].tolist())
                    data[coh_type][reg]['coh'].extend(vals[final_mask].tolist())
                    
    # Fit
    results = {}
    for coh_type in COHERENCE_TYPES:
        mode = f"{metric}/{coh_type}"
        results[coh_type] = {}
        print_status(f"  Fitting {mode}...")
        
        for reg in REGIONS:
            d = data[coh_type][reg]
            res = bin_and_fit(d['dist'], d['coh'])
            results[coh_type][reg] = res
            if res:
                print_status(f"    {reg:<12}: λ={res['lambda_km']:.0f}km, R²={res['r_squared']:.2f}")
            else:
                print_status(f"    {reg:<12}: Fit Failed", level="WARNING")
                
    gc.collect()
    return results

def plot_summary(final_results):
    # Plot Lambda comparisons
    fig, axes = plt.subplots(3, 2, figsize=(15, 18))
    fig.suptitle('TEP Regional Control Tests (Correlation Length λ)', fontsize=16)
    
    regions = REGIONS
    x = np.arange(len(regions))
    width = 0.35
    
    for i, metric in enumerate(METRICS):
        for j, coh in enumerate(COHERENCE_TYPES):
            ax = axes[i, j]
            if metric not in final_results: continue
            
            res = final_results[metric][coh]
            lambdas = [res[r]['lambda_km'] if res[r] else 0 for r in regions]
            errs = [res[r]['lambda_err'] if res[r] else 0 for r in regions]
            
            ax.bar(x, lambdas, yerr=errs, capsize=5, alpha=0.7, color='skyblue')
            ax.set_xticks(x)
            ax.set_xticklabels(regions, rotation=45)
            ax.set_title(f'{metric} / {coh}')
            ax.set_ylabel('Lambda (km)')
            ax.grid(axis='y', alpha=0.3)
            
            # Add value labels
            for k, v in enumerate(lambdas):
                ax.text(k, v + max(lambdas)*0.02, f"{v:.0f}", ha='center')
                
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(FIGURES_DIR / "step_2_1_control_tests.png")
    print_status(f"Saved figure to {FIGURES_DIR / 'step_2_1_control_tests.png'}", level="SUCCESS")

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
    
    parser = argparse.ArgumentParser(description='TEP-GNSS-RINEX Step 2.1: Control Tests')
    parser.add_argument('--filter', type=str, default='all',
                        help='Filter: "all", "none", "optimal_100_metadata.json", or "dynamic_50_metadata.json"')
    args = parser.parse_args()
    
    print_status("STEP 2.1: Control Tests (COMPREHENSIVE) - PARALLEL OPTIMIZED", level="TITLE")
    print_status(f"[PARALLEL] Using {N_WORKERS} workers (detected {cpu_count()} CPUs)")
    print_status("Processing: 3 filters × 3 modes × 3 metrics × 2 coherence types")
    
    # Select filters
    if args.filter == 'all':
        filters = STATION_FILTERS
    else:
        matching = [f for f in STATION_FILTERS if f[1] == args.filter]
        filters = matching if matching else [('ALL_STATIONS', 'none')]
    
    all_filter_results = {}
    
    for filter_name, filter_config in filters:
        filter_key = filter_name.lower()
        print_status(f"\n{'#'*80}")
        print_status(f"# FILTER: {filter_name}")
        print_status(f"{'#'*80}")
        
        station_set, filter_meta = load_station_filter(filter_config)
        if station_set:
            print_status(f"  Using {len(station_set)} stations from {filter_config}")
        
        all_mode_results = {}
        
        for mode_name, mode_config in PROCESSING_MODES.items():
            csv_file = RESULTS_DIR / mode_config['csv_file']
            
            print_status(f"\n{'='*60}")
            print_status(f"MODE: {mode_name.upper()} ({mode_config['description']})")
            print_status(f"{'='*60}")
            
            if not csv_file.exists():
                print_status(f"  WARNING: {csv_file.name} not found, skipping", level="WARNING")
                continue
            
            # Temporarily set global CSV_FILE for process_metric
            global CSV_FILE
            CSV_FILE = csv_file
            
            mode_results = {}
            
            # PARALLEL: Process all 3 metrics simultaneously
            # Pass csv_file and station_set explicitly for multiprocessing
            args_list = [(metric, csv_file, station_set) for metric in METRICS]
            
            with ProcessPoolExecutor(max_workers=min(3, N_WORKERS)) as executor:
                futures = {executor.submit(process_metric_worker, args): args[0] for args in args_list}
                for future in as_completed(futures):
                    metric = futures[future]
                    try:
                        mode_results[metric] = future.result()
                    except Exception as e:
                        print_status(f"  ERROR processing {metric}: {e}", level="ERROR")
                        mode_results[metric] = None
            
            all_mode_results[mode_name] = mode_results
        
        all_filter_results[filter_key] = {
            'filter_name': filter_name,
            'filter_config': filter_config,
            'filter_metadata': filter_meta,
            'results_by_mode': all_mode_results
        }
        
        # Save per-filter JSON
        output = {
            'step': '2.1',
            'name': 'control_tests',
            'filter': filter_name,
            'timestamp': datetime.now().isoformat(),
            'modes_processed': list(all_mode_results.keys()),
            'results_by_mode': all_mode_results,
            'results': all_mode_results.get('baseline', {})
        }
        output_file = RESULTS_DIR / f"step_2_1_control_tests_{filter_key}.json"
        with open(output_file, 'w') as f:
            json.dump(output, f, indent=2)
        print_status(f"Saved: {output_file.name}", level="SUCCESS")
    
    # Save comprehensive summary
    summary = {
        'step': '2.1',
        'name': 'control_tests_comprehensive',
        'timestamp': datetime.now().isoformat(),
        'filters_tested': [f[0] for f in filters],
        'modes_tested': list(PROCESSING_MODES.keys()),
        'results_by_filter': all_filter_results
    }
    with open(RESULTS_DIR / "step_2_1_control_tests_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Plot baseline results for first filter (for backward compat)
    first_filter = list(all_filter_results.keys())[0] if all_filter_results else None
    if first_filter and 'baseline' in all_filter_results[first_filter].get('results_by_mode', {}):
        plot_summary(all_filter_results[first_filter]['results_by_mode']['baseline'])
    
    print_status(f"COMPLETE: {len(filters)} filters × {len(PROCESSING_MODES)} modes", level="SUCCESS")

if __name__ == '__main__':
    main()
