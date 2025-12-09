#!/usr/bin/env python3
"""
TEP-GNSS-RINEX - STEP 2.3b: Temporal Stability Analysis (Full 6-Mode)
=====================================================================
Tests whether TEP signal is stable across years (2022, 2023, 2024).

Consistent temporal structure would indicate a persistent physical phenomenon,
not a transient artifact.

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
from scipy import stats
from datetime import datetime
import matplotlib.pyplot as plt
import gc
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
from tqdm import tqdm

from scripts.utils.logger import TEPLogger, set_step_logger, print_status

# Initialize logger
logger = TEPLogger(
    name="step_2_3_temporal_analysis",
    level="INFO",
    log_file_path=PROJECT_ROOT / "logs" / "step_2_3_temporal_analysis.log"
)
set_step_logger(logger)

# Parallel processing
N_WORKERS = min(cpu_count(), 16)  # Cap workers for fitting tasks

# ============================================================================
# ANALYSIS PARAMETERS
# ============================================================================
MIN_DISTANCE_KM = 50
MAX_DISTANCE_KM = 13000
N_BINS = 40
MIN_BIN_COUNT = 50

RESULTS_DIR = PROJECT_ROOT / "results" / "outputs"
FIGURES_DIR = PROJECT_ROOT / "results" / "figures"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"

METRICS = ['clock_bias', 'pos_jitter', 'clock_drift']
COHERENCE_TYPES = ['msc', 'phase_alignment']
YEARS = [2022, 2023, 2024]

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

def exp_decay(r, A, lam, C0):
    return A * np.exp(-r / lam) + C0

def fit_exponential(distances, coherences):
    if len(distances) < 1000:
        return None
    bin_edges = np.logspace(np.log10(MIN_DISTANCE_KM), np.log10(MAX_DISTANCE_KM), N_BINS + 1)
    bin_centers, bin_means, bin_counts = [], [], []
    for i in range(N_BINS):
        mask = (distances >= bin_edges[i]) & (distances < bin_edges[i+1])
        if np.sum(mask) >= MIN_BIN_COUNT:
            bin_centers.append((bin_edges[i] + bin_edges[i+1]) / 2)
            bin_means.append(np.nanmean(coherences[mask]))
            bin_counts.append(np.sum(mask))
    if len(bin_centers) < 5:
        return None
    bin_centers, bin_means, bin_counts = np.array(bin_centers), np.array(bin_means), np.array(bin_counts)
    try:
        popt, pcov = curve_fit(exp_decay, bin_centers, bin_means, p0=[0.5, 2000, 0],
                               sigma=1.0/np.sqrt(bin_counts), bounds=([0, 100, -1], [2, 20000, 1]), maxfev=10000)
        predicted = exp_decay(bin_centers, *popt)
        ss_res, ss_tot = np.sum((bin_means - predicted)**2), np.sum((bin_means - np.mean(bin_means))**2)
        r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0
        return {'lambda_km': float(popt[1]), 'lambda_err': float(np.sqrt(pcov[1,1])),
                'amplitude': float(popt[0]), 'offset': float(popt[2]),
                'r_squared': float(r2), 'n_pairs': int(sum(bin_counts))}
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


def run_temporal_analysis(csv_file, filter_name):
    """Run temporal analysis for a single CSV file."""
    if not csv_file.exists():
        print_status(f"  WARNING: {csv_file.name} not found", level="WARNING")
        return None
    
    print_status(f"\n  Loading {csv_file.name}...")
    
    # Initialize: year + 'all'
    periods = YEARS + ['all']
    data = {m: {c: {p: {'dist': [], 'coh': []} for p in periods} for c in COHERENCE_TYPES} for m in METRICS}
    
    total = 0
    for chunk in pd.read_csv(csv_file, chunksize=10_000_000):  # Larger chunks for high-RAM systems
        # Apply station filter if set
        if GOOD_STATIONS_SET:
            chunk = chunk[chunk['station1'].isin(GOOD_STATIONS_SET) & chunk['station2'].isin(GOOD_STATIONS_SET)]
        has_phase = 'phase_alignment' in chunk.columns
        years = chunk['year'].values
        
        for metric in METRICS:
            if 'metric' in chunk.columns:
                # Handle prefixed metric names (ionofree_*, multi_gnss_*)
                metric_mask = chunk['metric'].str.endswith(metric)
                mc = chunk[metric_mask]
                yr = years[metric_mask]
            else:
                mc = chunk
                yr = years
            
            if mc.empty:
                continue
            
            dist = mc['distance_km'].values
            msc = mc['coherence'].values
            phase = mc['phase_alignment'].values if has_phase else None
            
            # By year
            for year in YEARS:
                year_mask = (yr == year)
                
                # MSC
                valid = year_mask & ~np.isnan(msc)
                data[metric]['msc'][year]['dist'].extend(dist[valid].tolist())
                data[metric]['msc'][year]['coh'].extend(msc[valid].tolist())
                
                # Phase
                if phase is not None:
                    valid = year_mask & ~np.isnan(phase)
                    data[metric]['phase_alignment'][year]['dist'].extend(dist[valid].tolist())
                    data[metric]['phase_alignment'][year]['coh'].extend(phase[valid].tolist())
            
            # All years
            valid = ~np.isnan(msc)
            data[metric]['msc']['all']['dist'].extend(dist[valid].tolist())
            data[metric]['msc']['all']['coh'].extend(msc[valid].tolist())
            if phase is not None:
                valid = ~np.isnan(phase)
                data[metric]['phase_alignment']['all']['dist'].extend(dist[valid].tolist())
                data[metric]['phase_alignment']['all']['coh'].extend(phase[valid].tolist())
        
        total += len(chunk)
        print(f"  {total/1e6:.1f}M rows...", end='\r')  # Progress indicator
        gc.collect()
    
    print_status(f"  Total: {total:,} rows")
    
    # Fit and report
    print_status("TEMPORAL STABILITY RESULTS", level="TITLE")
    
    all_results = {}
    for metric in METRICS:
        for coh_type in COHERENCE_TYPES:
            mode = f"{metric}/{coh_type}"
            print_status(f"\n--- {mode} ---")
            print_status(f"{'Period':<12} {'λ (km)':<12} {'R²':<10} {'Pairs':>12}")
            print_status("-"*50)
            
            all_results[mode] = {}
            lambdas = []
            for period in periods:
                d = data[metric][coh_type][period]
                result = fit_exponential(np.array(d['dist']), np.array(d['coh']))
                all_results[mode][str(period)] = result
                if result:
                    print_status(f"{period:<12} {result['lambda_km']:<12.0f} {result['r_squared']:<10.3f} {result['n_pairs']:>12,}")
                    if period != 'all':
                        lambdas.append(result['lambda_km'])
                else:
                    print_status(f"{period:<12} {'N/A':<12}")
            
            # Coefficient of variation
            if len(lambdas) >= 2:
                cv = np.std(lambdas) / np.mean(lambdas) * 100
                print_status(f"  Year-to-year CV: {cv:.1f}%")
                if cv < 20:
                    print_status(f"  → STABLE (CV < 20%)", level="SUCCESS")
                else:
                    print_status(f"  → VARIABLE (CV ≥ 20%)", level="WARNING")
    
    # Summary
    print_status("TEMPORAL STABILITY SUMMARY", level="TITLE")
    print_status(f"{'Mode':<30} {'2022 λ':<10} {'2023 λ':<10} {'2024 λ':<10} {'CV%':<8}")
    print_status("-"*70)
    
    for metric in METRICS:
        for coh_type in COHERENCE_TYPES:
            mode = f"{metric}/{coh_type}"
            lambdas = []
            row = [mode]
            for year in YEARS:
                r = all_results[mode].get(str(year))
                if r:
                    row.append(f"{r['lambda_km']:.0f}")
                    lambdas.append(r['lambda_km'])
                else:
                    row.append("N/A")
            cv = np.std(lambdas) / np.mean(lambdas) * 100 if len(lambdas) >= 2 else 0
            row.append(f"{cv:.1f}")
            print_status(f"{row[0]:<30} {row[1]:<10} {row[2]:<10} {row[3]:<10} {row[4]:<8}")
    
    return all_results


def main():
    import argparse
    parser = argparse.ArgumentParser(description='TEP-GNSS-RINEX Step 2.3b: Temporal Analysis')
    parser.add_argument('--filter', type=str, default='all',
                        help='Filter: "all", "none", "optimal_100_metadata.json", or "dynamic_50_metadata.json"')
    args = parser.parse_args()
    
    print_status('STEP 2.3b: Temporal Stability Analysis (COMPREHENSIVE)', level="TITLE")
    print_status('Processing: 3 filters x 3 modes x 3 metrics x 2 coherence types')
    
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
            print_status(f"  Using {len(station_set)} stations")
        
        all_mode_results = {}
        
        for mode_name, mode_config in PROCESSING_MODES.items():
            print_status(f"\n  MODE: {mode_name.upper()}")
            csv_file = RESULTS_DIR / mode_config['csv_file']
            result = run_temporal_analysis(csv_file, filter_name)
            if result:
                all_mode_results[mode_name] = result
        
        all_filter_results[filter_key] = {
            'filter_name': filter_name,
            'filter_metadata': filter_meta,
            'results_by_mode': all_mode_results
        }
        
        # Save per-filter JSON
        output = {
            'step': '2.3b', 'name': 'temporal_stability',
            'filter': filter_name,
            'timestamp': datetime.now().isoformat(),
            'years': YEARS,
            'results_by_mode': all_mode_results
        }
        with open(RESULTS_DIR / f"step_2_3_temporal_analysis_{filter_key}.json", 'w') as f:
            json.dump(output, f, indent=2, default=str)
        print_status(f"  Saved: step_2_3_temporal_analysis_{filter_key}.json", level="SUCCESS")
    
    # Save summary
    summary = {
        'step': '2.3b', 'name': 'temporal_analysis_comprehensive',
        'timestamp': datetime.now().isoformat(),
        'filters_tested': [f[0] for f in filters],
        'modes_tested': list(PROCESSING_MODES.keys()),
        'results_by_filter': all_filter_results
    }
    with open(RESULTS_DIR / "step_2_3_temporal_analysis_summary.json", 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    print_status(f"COMPLETE: {len(filters)} filters x {len(PROCESSING_MODES)} modes", level="SUCCESS")


if __name__ == '__main__':
    main()
