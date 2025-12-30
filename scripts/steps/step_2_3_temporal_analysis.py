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

# ALL FOUR PROCESSING MODES
PROCESSING_MODES = {
    'baseline': {
        'description': 'GPS L1 only'
    },
    'ionofree': {
        'description': 'Dual-frequency L1+L2 ionosphere-free'
    },
    'multi_gnss': {
        'description': 'GPS+GLO+GAL+BDS'
    },
    'precise': {
        'description': 'IGS Precise Orbits (L1+L2)'
    }

}


def _sanitize_tag(tag: str) -> str:
    tag = (tag or "").strip()
    if not tag:
        return ""
    safe = []
    for ch in tag:
        if ch.isalnum() or ch in ('-', '_'):
            safe.append(ch)
    return "_" + "".join(safe) if safe else ""

def exp_decay(r, A, lam, C0):
    return A * np.exp(-r / lam) + C0

def fit_exponential_binned(bin_centers, bin_means, bin_counts):
    if bin_centers.size < 5:
        return None
    try:
        lam_lower = 100.0
        lam_upper = 20000.0
        popt, pcov = curve_fit(
            exp_decay,
            bin_centers,
            bin_means,
            p0=[0.5, 2000, 0],
            sigma=1.0 / np.sqrt(bin_counts),
            bounds=([0, lam_lower, -1], [2, lam_upper, 1]),
            maxfev=10000,
        )

        predicted = exp_decay(bin_centers, *popt)

        weights = bin_counts.astype(float)
        weighted_mean = np.average(bin_means, weights=weights)
        ss_res = np.sum(weights * (bin_means - predicted) ** 2)
        ss_tot = np.sum(weights * (bin_means - weighted_mean) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

        lam = float(popt[1])
        lambda_hit_lower_bound = bool(lam <= lam_lower * 1.01)
        lambda_hit_upper_bound = bool(lam >= lam_upper * 0.99)
        lam_err = float(np.sqrt(pcov[1, 1])) if pcov.shape == (3, 3) and np.isfinite(pcov[1, 1]) else float('nan')
        lam_rel_err = float(lam_err / lam) if lam > 0 and np.isfinite(lam_err) else float('nan')

        max_bin_center_km = float(np.max(bin_centers))
        lambda_unconstrained = bool(
            lambda_hit_upper_bound
            or (np.isfinite(lam_rel_err) and lam_rel_err > 0.5)
            or (max_bin_center_km > 0 and lam > 1.5 * max_bin_center_km)
        )

        return {
            'lambda_km': lam,
            'lambda_err': lam_err,
            'lambda_rel_err': lam_rel_err,
            'lambda_hit_lower_bound': lambda_hit_lower_bound,
            'lambda_hit_upper_bound': lambda_hit_upper_bound,
            'lambda_unconstrained': lambda_unconstrained,
            'amplitude': float(popt[0]),
            'offset': float(popt[2]),
            'r_squared': float(r2),
            'n_pairs': int(np.sum(bin_counts)),
            'n_bins': int(bin_centers.size),
            'max_bin_center_km': max_bin_center_km,
        }
    except Exception:
        return None


def _init_binner():
    return {
        'sum': np.zeros(N_BINS, dtype=np.float64),
        'count': np.zeros(N_BINS, dtype=np.int64),
    }


def _accumulate_binner(binner, bin_idx, values):
    if values.size == 0:
        return
    binner['sum'] += np.bincount(bin_idx, weights=values, minlength=N_BINS)
    binner['count'] += np.bincount(bin_idx, minlength=N_BINS)


def _finalize_binner(binner, bin_edges):
    counts = binner['count']
    ok = counts >= MIN_BIN_COUNT
    if not np.any(ok):
        return None, None, None
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    means = binner['sum'][ok] / counts[ok]
    return bin_centers[ok], means, counts[ok]

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
    
    periods = YEARS + ['all']
    bin_edges = np.logspace(np.log10(MIN_DISTANCE_KM), np.log10(MAX_DISTANCE_KM), N_BINS + 1)
    binners = {
        m: {
            c: {p: _init_binner() for p in periods}
            for c in COHERENCE_TYPES
        }
        for m in METRICS
    }
    
    total = 0
    for chunk in pd.read_csv(csv_file, chunksize=10_000_000):
        if GOOD_STATIONS_SET:
            chunk = chunk[chunk['station1'].isin(GOOD_STATIONS_SET) & chunk['station2'].isin(GOOD_STATIONS_SET)]

        has_phase = 'phase_alignment' in chunk.columns

        for metric in METRICS:
            if 'metric' in chunk.columns:
                metric_mask = chunk['metric'].str.endswith(metric)
                mc = chunk[metric_mask]
            else:
                mc = chunk

            if mc.empty:
                continue

            yr = mc['year'].to_numpy(dtype=np.int64, copy=False)
            dist = mc['distance_km'].to_numpy(dtype=np.float64, copy=False)
            msc = mc['coherence'].to_numpy(dtype=np.float64, copy=False)
            phase = mc['phase_alignment'].to_numpy(dtype=np.float64, copy=False) if has_phase else None

            bin_idx = np.searchsorted(bin_edges, dist, side='right') - 1
            in_range = (bin_idx >= 0) & (bin_idx < N_BINS)

            valid_msc_all = in_range & ~np.isnan(msc)
            if np.any(valid_msc_all):
                _accumulate_binner(binners[metric]['msc']['all'], bin_idx[valid_msc_all], msc[valid_msc_all])

            if phase is not None:
                valid_phase_all = in_range & ~np.isnan(phase)
                if np.any(valid_phase_all):
                    _accumulate_binner(binners[metric]['phase_alignment']['all'], bin_idx[valid_phase_all], phase[valid_phase_all])

            for year in YEARS:
                year_mask = yr == year

                valid_msc = valid_msc_all & year_mask
                if np.any(valid_msc):
                    _accumulate_binner(binners[metric]['msc'][year], bin_idx[valid_msc], msc[valid_msc])

                if phase is not None:
                    valid_phase = (in_range & ~np.isnan(phase) & year_mask)
                    if np.any(valid_phase):
                        _accumulate_binner(binners[metric]['phase_alignment'][year], bin_idx[valid_phase], phase[valid_phase])

        total += len(chunk)
        print(f"  {total/1e6:.1f}M rows...", end='\r')
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
            unconstrained_years = []
            for period in periods:
                bc, bm, bn = _finalize_binner(binners[metric][coh_type][period], bin_edges)
                result = fit_exponential_binned(bc, bm, bn) if bc is not None else None
                all_results[mode][str(period)] = result
                if result:
                    print_status(f"{period:<12} {result['lambda_km']:<12.0f} {result['r_squared']:<10.3f} {result['n_pairs']:>12,}")
                    if period != 'all':
                        lambdas.append(result['lambda_km'])
                        if result.get('lambda_unconstrained'):
                            unconstrained_years.append(str(period))
                        else:
                            rel_err = result.get('lambda_rel_err', float('nan'))
                            if np.isfinite(rel_err) and rel_err > 0.5:
                                unconstrained_years.append(str(period))
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

                if unconstrained_years:
                    print_status(f"  → WARNING: Potentially unconstrained λ in years: {', '.join(unconstrained_years)} (large relative fit uncertainty)", level="WARNING")
    
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
    parser.add_argument('--tag', type=str, default='',
                        help='Optional tag appended to output filenames (for non-overwriting reruns).')
    args = parser.parse_args()
    
    print_status('STEP 2.3b: Temporal Stability Analysis (COMPREHENSIVE)', level="TITLE")
    run_tag = _sanitize_tag(args.tag)
    
    if args.filter == 'all':
        filters = STATION_FILTERS
    else:
        matching = [f for f in STATION_FILTERS if f[1] == args.filter or f[0] == args.filter]
        filters = matching if matching else [('ALL_STATIONS', 'none')]

    print_status(f"Processing: {len(filters)} filters x {len(PROCESSING_MODES)} modes x {len(METRICS)} metrics x {len(COHERENCE_TYPES)} coherence types")
    
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
            # Construct filename dynamically: step_2_0_pairs_{mode}_{filter}.csv
            csv_filename = f"step_2_0_pairs_{mode_name}_{filter_key}.csv"
            csv_file = RESULTS_DIR / csv_filename
            
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
        with open(RESULTS_DIR / f"step_2_3_temporal_analysis_{filter_key}{run_tag}.json", 'w') as f:
            json.dump(output, f, indent=2, default=str)
        print_status(f"  Saved: step_2_3_temporal_analysis_{filter_key}{run_tag}.json", level="SUCCESS")
    
    # Save summary
    summary = {
        'step': '2.3b', 'name': 'temporal_analysis_comprehensive',
        'timestamp': datetime.now().isoformat(),
        'filters_tested': [f[0] for f in filters],
        'modes_tested': list(PROCESSING_MODES.keys()),
        'results_by_filter': all_filter_results
    }
    with open(RESULTS_DIR / f"step_2_3_temporal_analysis_summary{run_tag}.json", 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    print_status(f"COMPLETE: {len(filters)} filters x {len(PROCESSING_MODES)} modes", level="SUCCESS")


if __name__ == '__main__':
    main()
