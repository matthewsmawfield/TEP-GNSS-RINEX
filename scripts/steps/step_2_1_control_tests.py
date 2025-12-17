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
from scripts.utils.tep_utils import (
    MIN_DISTANCE_KM,
    MAX_DISTANCE_KM,
    N_BINS,
    MIN_BIN_COUNT,
    fit_exponential_binned,
    ecef_to_lla,
    load_station_map,
)
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

# ALL FOUR PROCESSING MODES - consistent across all Step 2.x scripts
PROCESSING_MODES = {
    'precise': {
        'csv_file': 'step_2_0_pairs_precise.csv',
        'description': 'GPS SPP with Precise Orbits'
    },
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

def _resolve_pairs_csv(mode_name: str, filter_key: str) -> Path:
    """Resolve Step 2.0 pair CSV path for a given mode and station filter.

    Step 2.0 writes files with a filter suffix (e.g. *_all_stations.csv). Older
    runs may have unsuffixed filenames. This function supports both.
    """
    candidates = [
        RESULTS_DIR / f"step_2_0_pairs_{mode_name}_{filter_key}.csv",
        RESULTS_DIR / f"step_2_0_pairs_{mode_name}.csv",
    ]
    for p in candidates:
        if p.exists():
            return p
    return candidates[0]

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

def process_metric(metric, csv_file=None, station_filter_set=None):
    """Process a single metric. Can be called directly or via worker."""
    print_status(f"Processing {metric}...")
    
    # Use passed csv_file or fall back to global
    if csv_file is None:
        csv_file = CSV_FILE

    bin_edges = np.logspace(np.log10(MIN_DISTANCE_KM), np.log10(MAX_DISTANCE_KM), N_BINS + 1)
    bin_centers_all = np.sqrt(bin_edges[:-1] * bin_edges[1:])

    diagnostics = {
        coh_type: {reg: {} for reg in REGIONS}
        for coh_type in COHERENCE_TYPES
    }

    acc = {
        coh_type: {
            reg: {
                'sum': np.zeros(N_BINS, dtype=np.float64),
                'count': np.zeros(N_BINS, dtype=np.int64),
            }
            for reg in REGIONS
        }
        for coh_type in COHERENCE_TYPES
    }
    
    # Load Station Map
    coords = load_station_map(PROCESSED_DIR)
    if not coords:
        print_status("Error: Missing station coordinates.", level="ERROR")
        return None
    
    # Apply station filter if set
    if station_filter_set:
        coords = {k: v for k, v in coords.items() if k in station_filter_set}
        print_status(f"  Filtered to {len(coords)} stations")
    elif GOOD_STATIONS_SET:
        coords = {k: v for k, v in coords.items() if k in GOOD_STATIONS_SET}
        print_status(f"  Filtered to {len(coords)} stations")
    
    r_map = get_station_regions(coords)

    valid_stations = set(r_map.keys())
    is_europe = {k: v[0] for k, v in r_map.items()}
    is_northern = {k: v[1] for k, v in r_map.items()}

    header_cols = pd.read_csv(csv_file, nrows=0).columns.tolist()
    has_phase_col = 'phase_alignment' in header_cols

    usecols = ['metric', 'station1', 'station2', 'distance_km', 'coherence']
    if has_phase_col:
        usecols.append('phase_alignment')
    dtype = {
        'distance_km': 'float32',
        'coherence': 'float32',
    }
    if has_phase_col:
        dtype['phase_alignment'] = 'float32'

    for chunk in pd.read_csv(csv_file, chunksize=1_000_000, usecols=usecols, dtype=dtype):
        mc = chunk[chunk['metric'].str.endswith(metric)]
        if mc.empty:
            continue

        s1 = mc['station1']
        s2 = mc['station2']
        in_filter = s1.isin(valid_stations) & s2.isin(valid_stations)
        if not in_filter.any():
            continue

        dist = mc['distance_km'].to_numpy(dtype=np.float32, copy=False)
        msc = mc['coherence'].to_numpy(dtype=np.float32, copy=False)
        phase = mc['phase_alignment'].to_numpy(dtype=np.float32, copy=False) if has_phase_col else None

        eur1 = s1.map(is_europe).fillna(False).to_numpy(dtype=bool, copy=False)
        eur2 = s2.map(is_europe).fillna(False).to_numpy(dtype=bool, copy=False)
        north1 = s1.map(is_northern).fillna(False).to_numpy(dtype=bool, copy=False)
        north2 = s2.map(is_northern).fillna(False).to_numpy(dtype=bool, copy=False)

        in_filter_arr = in_filter.to_numpy(dtype=bool, copy=False)
        mask_global = in_filter_arr
        mask_europe = in_filter_arr & eur1 & eur2
        mask_non_europe = in_filter_arr & (~eur1) & (~eur2)
        mask_northern = in_filter_arr & north1 & north2
        mask_southern = in_filter_arr & (~north1) & (~north2)

        masks = {
            'global': mask_global,
            'europe': mask_europe,
            'non_europe': mask_non_europe,
            'northern': mask_northern,
            'southern': mask_southern,
        }

        bin_idx = np.searchsorted(bin_edges, dist, side='right') - 1
        in_range = (dist >= MIN_DISTANCE_KM) & (dist <= MAX_DISTANCE_KM)
        valid_bin = in_range & (bin_idx >= 0) & (bin_idx < N_BINS)
        if not np.any(valid_bin):
            continue

        for coh_type in COHERENCE_TYPES:
            vals = msc if coh_type == 'msc' else phase
            if vals is None:
                continue

            finite = np.isfinite(vals)
            for reg, reg_mask in masks.items():
                final_mask = valid_bin & reg_mask & finite
                if not np.any(final_mask):
                    continue

                idx = bin_idx[final_mask]
                v = vals[final_mask]
                acc[coh_type][reg]['sum'] += np.bincount(idx, weights=v, minlength=N_BINS)
                acc[coh_type][reg]['count'] += np.bincount(idx, minlength=N_BINS).astype(np.int64)
                    
    # Fit
    results = {}
    for coh_type in COHERENCE_TYPES:
        mode = f"{metric}/{coh_type}"
        results[coh_type] = {}
        print_status(f"  Fitting {mode}...")
        
        for reg in REGIONS:
            counts = acc[coh_type][reg]['count']
            valid = counts >= MIN_BIN_COUNT
            n_bins_used = int(np.sum(valid))
            n_pairs_total = int(np.sum(counts))
            n_pairs_used = int(np.sum(counts[valid]))

            diag = {
                'n_bins_used': n_bins_used,
                'n_pairs_total': n_pairs_total,
                'n_pairs_used': n_pairs_used,
                'has_fit': False,
            }

            if n_bins_used < 5:
                res = None
            else:
                sums = acc[coh_type][reg]['sum']
                bin_centers = bin_centers_all[valid]
                bin_means = (sums[valid] / counts[valid]).astype(np.float64)
                bin_counts = counts[valid].astype(np.int64)
                res = fit_exponential_binned(bin_centers, bin_means, bin_counts)

            results[coh_type][reg] = res
            if res:
                diag['has_fit'] = True
                diag['lambda_km'] = float(res.get('lambda_km'))
                diag['r_squared'] = float(res.get('r_squared'))
                diag['lambda_over_max_distance'] = float(res.get('lambda_km')) / float(MAX_DISTANCE_KM)
                diag['flat_decay_flag'] = bool(res.get('lambda_km') is not None and res.get('lambda_km') > MAX_DISTANCE_KM)
                print_status(f"    {reg:<12}: λ={res['lambda_km']:.0f}km, R²={res['r_squared']:.2f}")
            else:
                print_status(f"    {reg:<12}: Fit Failed", level="WARNING")

            diagnostics[coh_type][reg] = diag
                
    results['diagnostics'] = diagnostics
    gc.collect()
    return results

def _summarize_fit_block(metric_results, coh_type):
    reg_fits = metric_results.get(coh_type, {}) if metric_results else {}
    lambdas = []
    r2s = []
    for reg in REGIONS:
        res = reg_fits.get(reg)
        if not res:
            continue
        if res.get('lambda_km') is None or res.get('r_squared') is None:
            continue
        lambdas.append(float(res['lambda_km']))
        r2s.append(float(res['r_squared']))

    lambdas_arr = np.asarray(lambdas, dtype=np.float64)
    r2s_arr = np.asarray(r2s, dtype=np.float64)

    out = {
        'n_regions_fit': int(lambdas_arr.size),
        'global_lambda_km': None,
        'global_r_squared': None,
        'global_n_pairs_used': None,
        'lambda_median_km': None,
        'lambda_mean_km': None,
        'lambda_std_km': None,
        'lambda_cv': None,
        'r2_median': None,
        'r2_mean': None,
        'any_flat_decay_flag': False,
    }

    global_res = reg_fits.get('global') if reg_fits else None
    if global_res:
        out['global_lambda_km'] = float(global_res.get('lambda_km')) if global_res.get('lambda_km') is not None else None
        out['global_r_squared'] = float(global_res.get('r_squared')) if global_res.get('r_squared') is not None else None

    diag = metric_results.get('diagnostics', {}).get(coh_type, {}) if metric_results else {}
    if diag.get('global'):
        out['global_n_pairs_used'] = int(diag['global'].get('n_pairs_used')) if diag['global'].get('n_pairs_used') is not None else None

    if lambdas_arr.size:
        out['lambda_median_km'] = float(np.median(lambdas_arr))
        out['lambda_mean_km'] = float(np.mean(lambdas_arr))
        out['lambda_std_km'] = float(np.std(lambdas_arr, ddof=0))
        if out['lambda_mean_km'] and out['lambda_mean_km'] > 0:
            out['lambda_cv'] = float(out['lambda_std_km'] / out['lambda_mean_km'])

    if r2s_arr.size:
        out['r2_median'] = float(np.median(r2s_arr))
        out['r2_mean'] = float(np.mean(r2s_arr))

    out['any_flat_decay_flag'] = any(bool(diag.get(reg, {}).get('flat_decay_flag')) for reg in REGIONS)
    return out

def _summarize_mode(mode_results):
    summary = {
        'metrics': {},
        'bias_drift_consistency': {},
        'flags': {
            'any_flat_decay': False,
        },
    }

    for metric in METRICS:
        mr = mode_results.get(metric)
        if not mr:
            continue
        summary['metrics'][metric] = {}
        for coh_type in COHERENCE_TYPES:
            summary['metrics'][metric][coh_type] = _summarize_fit_block(mr, coh_type)
            if summary['metrics'][metric][coh_type].get('any_flat_decay_flag'):
                summary['flags']['any_flat_decay'] = True

    for coh_type in COHERENCE_TYPES:
        b = mode_results.get('clock_bias', {}).get(coh_type, {}).get('global') if mode_results.get('clock_bias') else None
        d = mode_results.get('clock_drift', {}).get(coh_type, {}).get('global') if mode_results.get('clock_drift') else None

        if b and d and b.get('lambda_km') is not None and d.get('lambda_km') is not None:
            lb = float(b['lambda_km'])
            ld = float(d['lambda_km'])
            denom = (lb + ld) / 2.0 if (lb + ld) else None
            rel = abs(lb - ld) / denom if denom else None
            summary['bias_drift_consistency'][coh_type] = {
                'lambda_bias_km': lb,
                'lambda_drift_km': ld,
                'delta_lambda_km': float(lb - ld),
                'relative_delta': float(rel) if rel is not None else None,
                'r2_bias': float(b.get('r_squared')) if b.get('r_squared') is not None else None,
                'r2_drift': float(d.get('r_squared')) if d.get('r_squared') is not None else None,
            }
        else:
            summary['bias_drift_consistency'][coh_type] = None

    return summary

def _print_mode_summary(mode_name, mode_config, mode_results, mode_summary):
    print_status(f"\n{'-'*60}")
    print_status(f"SUMMARY: {mode_name.upper()} ({mode_config['description']})")
    print_status(f"{'-'*60}")

    for metric in METRICS:
        if metric not in mode_summary.get('metrics', {}):
            continue
        for coh_type in COHERENCE_TYPES:
            s = mode_summary['metrics'][metric].get(coh_type, {})
            gl = s.get('global_lambda_km')
            gr2 = s.get('global_r_squared')
            gnp = s.get('global_n_pairs_used')
            med_l = s.get('lambda_median_km')
            cv = s.get('lambda_cv')
            med_r2 = s.get('r2_median')
            nreg = s.get('n_regions_fit')

            parts = []
            if gl is not None and gr2 is not None:
                parts.append(f"global λ={gl:.0f}km, R²={gr2:.2f}")
            if gnp is not None:
                parts.append(f"pairs_used={gnp:,}")
            if med_l is not None and med_r2 is not None:
                parts.append(f"median(λ)={med_l:.0f}km, median(R²)={med_r2:.2f}")
            if cv is not None:
                parts.append(f"CV(λ)={cv:.2f}")
            parts.append(f"regions_fit={nreg}")

            flag = "" if not s.get('any_flat_decay_flag') else " [FLAG: flat/very-long λ]"
            print_status(f"  {metric}/{coh_type}: " + "; ".join(parts) + flag)

    for coh_type, bd in (mode_summary.get('bias_drift_consistency') or {}).items():
        if not bd:
            continue
        rel = bd.get('relative_delta')
        r2b = bd.get('r2_bias')
        r2d = bd.get('r2_drift')
        tag = "" if rel is None or rel <= 0.20 else " [WARNING: bias–drift mismatch]"
        r2_part = "" if (r2b is None or r2d is None) else f", R²_bias={r2b:.2f}, R²_drift={r2d:.2f}"
        print_status(
            f"  bias–drift ({coh_type}/global): λ_bias={bd['lambda_bias_km']:.0f}km, λ_drift={bd['lambda_drift_km']:.0f}km, Δλ={bd['delta_lambda_km']:.0f}km, rel={rel:.2f}" + r2_part + tag
        )

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
    parser.add_argument('--parallel-metrics', action='store_true',
                        help='Process the three metrics in parallel per mode (higher I/O load).')
    args = parser.parse_args()
    
    print_status("STEP 2.1: Control Tests (COMPREHENSIVE) - STREAMING", level="TITLE")
    if args.parallel_metrics:
        print_status(f"[PARALLEL] Metrics parallelism enabled: up to {min(N_WORKERS, len(METRICS))} workers (detected {cpu_count()} CPUs)")
    else:
        print_status("[SEQUENTIAL] Metrics processed sequentially (recommended for low I/O contention)")
    print_status(f"Processing: 3 filters × {len(PROCESSING_MODES)} modes × {len(METRICS)} metrics × {len(COHERENCE_TYPES)} coherence types")
    
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
        all_mode_summaries = {}
        
        for mode_name, mode_config in PROCESSING_MODES.items():
            csv_file = _resolve_pairs_csv(mode_name, filter_key)
            
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

            if args.parallel_metrics:
                max_workers = min(N_WORKERS, len(METRICS))
                with ProcessPoolExecutor(max_workers=max_workers) as executor:
                    futures = {
                        executor.submit(process_metric, metric, csv_file, station_set): metric
                        for metric in METRICS
                    }
                    for fut in as_completed(futures):
                        metric = futures[fut]
                        try:
                            mode_results[metric] = fut.result()
                        except Exception as e:
                            print_status(f"  ERROR processing {metric}: {e}", level="ERROR")
                            mode_results[metric] = None
            else:
                for metric in METRICS:
                    try:
                        mode_results[metric] = process_metric(metric, csv_file, station_set)
                    except Exception as e:
                        print_status(f"  ERROR processing {metric}: {e}", level="ERROR")
                        mode_results[metric] = None

            mode_summary = _summarize_mode(mode_results)
            _print_mode_summary(mode_name, mode_config, mode_results, mode_summary)
            all_mode_summaries[mode_name] = mode_summary
            
            all_mode_results[mode_name] = mode_results
        
        all_filter_results[filter_key] = {
            'filter_name': filter_name,
            'filter_config': filter_config,
            'filter_metadata': filter_meta,
            'results_by_mode': all_mode_results,
            'summaries_by_mode': all_mode_summaries,
        }
        
        # Save per-filter JSON
        output = {
            'step': '2.1',
            'name': 'control_tests',
            'filter': filter_name,
            'timestamp': datetime.now().isoformat(),
            'modes_processed': list(all_mode_results.keys()),
            'results_by_mode': all_mode_results,
            'summaries_by_mode': all_mode_summaries,
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
