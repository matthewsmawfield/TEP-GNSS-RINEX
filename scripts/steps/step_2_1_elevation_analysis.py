#!/usr/bin/env python3
"""
TEP-GNSS-RINEX - STEP 2.1b: Elevation Analysis (Optimized)
==========================================================
Tests systematic variation of correlation length (λ) with station altitude.
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
from datetime import datetime
from scripts.utils.tep_utils import bin_and_fit, ecef_to_lla, load_station_map
from scripts.utils.logger import TEPLogger, set_step_logger, print_status

# Initialize logger
logger = TEPLogger(
    name="step_2_1_elevation_analysis",
    level="INFO",
    log_file_path=PROJECT_ROOT / "logs" / "step_2_1_elevation_analysis.log"
)
set_step_logger(logger)

# Parameters
RESULTS_DIR = PROJECT_ROOT / "results" / "outputs"
FIGURES_DIR = PROJECT_ROOT / "results" / "figures"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"

METRICS = ['clock_bias', 'pos_jitter', 'clock_drift']
COHERENCE_TYPES = ['msc', 'phase_alignment']
QUINTILES = ['Q1', 'Q2', 'Q3', 'Q4', 'Q5']

# Station filters - consistent with other Step 2.x scripts
STATION_FILTERS = [
    ('ALL_STATIONS', 'none'),
    ('OPTIMAL_100', 'optimal_100_metadata.json'),
    ('DYNAMIC_50', 'dynamic_50_metadata.json')
]

# Global station filter set
GOOD_STATIONS_SET = None

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

def get_station_quintiles(coords):
    """Map station to quintile (Q1=Low, Q5=High)."""
    stas = []
    for s, xyz in coords.items():
        _, _, alt = ecef_to_lla(*xyz)
        stas.append((s, alt))
    
    stas.sort(key=lambda x: x[1])
    n = len(stas)
    q_size = n // 5
    
    q_map = {}
    q_stats = {}
    for i in range(5):
        q_name = f"Q{i+1}"
        start = i * q_size
        end = (i+1) * q_size if i < 4 else n
        subset = stas[start:end]
        alts = [x[1] for x in subset]
        
        q_stats[q_name] = {'min': min(alts), 'max': max(alts), 'count': len(subset)}
        for s, alt in subset:
            q_map[s] = q_name
            
    return q_map, q_stats

def process_metric(metric):
    print_status(f"\nProcessing {metric}...")
    
    data = {c: {q: {'dist': [], 'coh': []} for q in QUINTILES} for c in COHERENCE_TYPES}
    chunk_size = 10_000_000  # Larger chunks for high-RAM systems
    
    # Load Station Map
    coords = load_station_map(PROCESSED_DIR)
    if not coords: return None
    
    # Apply station filter if set
    if GOOD_STATIONS_SET:
        coords = {k: v for k, v in coords.items() if k in GOOD_STATIONS_SET}
        print_status(f"  Filtered to {len(coords)} stations")
    
    q_map, q_stats = get_station_quintiles(coords)
    
    for k, v in q_stats.items():
        print_status(f"  {k}: {v['min']:.0f}m - {v['max']:.0f}m ({v['count']} stations)")
    
    for chunk in pd.read_csv(CSV_FILE, chunksize=chunk_size):
        if 'metric' in chunk.columns:
            # Handle prefixed metric names (ionofree_*, multi_gnss_*)
            mc = chunk[chunk['metric'].str.endswith(metric)]
        else: continue
        if mc.empty: continue
        
        s1 = mc['station1'].values
        s2 = mc['station2'].values
        dist = mc['distance_km'].values
        msc = mc['coherence'].values
        has_phase = 'phase_alignment' in mc.columns
        phase = mc['phase_alignment'].values if has_phase else None
        
        # Map stations to quintiles
        # Optimization: q_map lookup
        q1 = [q_map.get(s, None) for s in s1]
        q2 = [q_map.get(s, None) for s in s2]
        
        # Filter: Both in same quintile
        # Create boolean masks for each Q
        for q in QUINTILES:
            # Vectorized check? List comprehension is reasonable for 5M items if simple string compare
            mask = np.array([qa == q and qb == q for qa, qb in zip(q1, q2)])
            
            if np.any(mask):
                for coh_type in COHERENCE_TYPES:
                    vals = msc if coh_type == 'msc' else phase
                    if vals is None: continue
                    
                    valid = ~np.isnan(vals) & mask
                    if np.any(valid):
                        data[coh_type][q]['dist'].extend(dist[valid].tolist())
                        data[coh_type][q]['coh'].extend(vals[valid].tolist())
                        
    # Fit
    results = {}
    for coh_type in COHERENCE_TYPES:
        mode = f"{metric}/{coh_type}"
        results[coh_type] = {}
        print_status(f"  Fitting {mode}...")
        
        for q in QUINTILES:
            d = data[coh_type][q]
            res = bin_and_fit(d['dist'], d['coh'])
            results[coh_type][q] = res
            if res:
                print_status(f"    {q:<4}: λ={res['lambda_km']:.0f}km, R²={res['r_squared']:.2f}")
            else:
                print_status(f"    {q:<4}: Fit Failed", level="WARNING")
                
    gc.collect()
    return results

def plot_summary(final_results):
    fig, axes = plt.subplots(3, 2, figsize=(15, 18))
    fig.suptitle('TEP Elevation Dependence (Correlation Length λ vs Altitude Quintile)', fontsize=16)
    
    x = np.arange(5) # Q1..Q5
    
    for i, metric in enumerate(METRICS):
        for j, coh in enumerate(COHERENCE_TYPES):
            ax = axes[i, j]
            if metric not in final_results: continue
            
            res = final_results[metric][coh]
            lambdas = [res[q]['lambda_km'] if res[q] else np.nan for q in QUINTILES]
            errs = [res[q]['lambda_err'] if res[q] else 0 for q in QUINTILES]
            
            # Plot line with error bars
            ax.errorbar(x, lambdas, yerr=errs, fmt='o-', capsize=5, linewidth=2, color='darkgreen')
            
            # Linear trend
            valid_idx = np.where(~np.isnan(lambdas))[0]
            if len(valid_idx) >= 3:
                z = np.polyfit(valid_idx, np.array(lambdas)[valid_idx], 1)
                p = np.poly1d(z)
                ax.plot(x, p(x), 'r--', alpha=0.5, label=f'Slope: {z[0]:.1f} km/quintile')
                ax.legend()
            
            ax.set_xticks(x)
            ax.set_xticklabels(QUINTILES)
            ax.set_title(f'{metric} / {coh}')
            ax.set_ylabel('Lambda (km)')
            ax.set_xlabel('Station Altitude Quintile (Low -> High)')
            ax.grid(True, alpha=0.3)
                
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(FIGURES_DIR / "step_2_1_elevation_analysis.png")
    print_status(f"Saved figure to {FIGURES_DIR / 'step_2_1_elevation_analysis.png'}", level="SUCCESS")

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
    parser = argparse.ArgumentParser(description='TEP-GNSS-RINEX Step 2.1b: Elevation Analysis')
    parser.add_argument('--filter', type=str, default='all',
                        help='Filter: "all", "none", "optimal_100_metadata.json", or "dynamic_50_metadata.json"')
    args = parser.parse_args()
    
    print_status("STEP 2.1b: Elevation Analysis (COMPREHENSIVE)", level="TITLE")
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
            for metric in METRICS:
                mode_results[metric] = process_metric(metric)
            
            all_mode_results[mode_name] = mode_results
        
        all_filter_results[filter_key] = {
            'filter_name': filter_name,
            'filter_config': filter_config,
            'filter_metadata': filter_meta,
            'results_by_mode': all_mode_results
        }
        
        # Save per-filter JSON
        output = {
            'step': '2.1b',
            'name': 'elevation_analysis',
            'filter': filter_name,
            'timestamp': datetime.now().isoformat(),
            'modes_processed': list(all_mode_results.keys()),
            'results_by_mode': all_mode_results,
            'results': all_mode_results.get('baseline', {})
        }
        output_file = RESULTS_DIR / f"step_2_1_elevation_analysis_{filter_key}.json"
        with open(output_file, 'w') as f:
            json.dump(output, f, indent=2)
        print_status(f"Saved: {output_file.name}", level="SUCCESS")
    
    # Save comprehensive summary
    summary = {
        'step': '2.1b',
        'name': 'elevation_analysis_comprehensive',
        'timestamp': datetime.now().isoformat(),
        'filters_tested': [f[0] for f in filters],
        'modes_tested': list(PROCESSING_MODES.keys()),
        'results_by_filter': all_filter_results
    }
    with open(RESULTS_DIR / "step_2_1_elevation_analysis_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Plot baseline results for first filter
    first_filter = list(all_filter_results.keys())[0] if all_filter_results else None
    if first_filter and 'baseline' in all_filter_results[first_filter].get('results_by_mode', {}):
        plot_summary(all_filter_results[first_filter]['results_by_mode']['baseline'])
    
    print_status(f"COMPLETE: {len(filters)} filters × {len(PROCESSING_MODES)} modes", level="SUCCESS")

if __name__ == '__main__':
    main()
