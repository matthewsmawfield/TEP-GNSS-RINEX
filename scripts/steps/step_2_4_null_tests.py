#!/usr/bin/env python3
"""
TEP-GNSS-RINEX - STEP 2.4b: Null Tests (Full 6-Mode) - Memory Optimized
=======================================================================
Tests that TEP signal does NOT correlate with known non-gravitational phenomena:
1. Solar Rotation (27-day cycle)
2. Lunar Month (29.5-day cycle)
3. Random Date Shuffling

Memory Optimized: Processes metrics sequentially to avoid OOM.

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
    name="step_2_4_null_tests",
    level="INFO",
    log_file_path=PROJECT_ROOT / "logs" / "step_2_4_null_tests.log"
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

# ALL FOUR PROCESSING MODES
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
    },
    'precise': {
        'csv_file': 'step_2_0_pairs_precise.csv',
        'description': 'IGS precise orbits/clocks (SP3)'
    }
}

# Cycle periods to test
SOLAR_ROTATION_DAYS = 27
LUNAR_MONTH_DAYS = 29.5

def exp_decay(r, A, lam, C0):
    return A * np.exp(-r / lam) + C0

def fit_exponential(distances, coherences):
    if len(distances) < 1000:
        return None
    bin_edges = np.logspace(np.log10(MIN_DISTANCE_KM), np.log10(MAX_DISTANCE_KM), N_BINS + 1)
    bin_centers, bin_means, bin_counts = [], [], []
    
    # Vectorized binning for speed
    distances = np.asarray(distances)
    coherences = np.asarray(coherences)
    
    for i in range(N_BINS):
        mask = (distances >= bin_edges[i]) & (distances < bin_edges[i+1])
        count = np.sum(mask)
        if count >= MIN_BIN_COUNT:
            bin_centers.append((bin_edges[i] + bin_edges[i+1]) / 2)
            bin_means.append(np.nanmean(coherences[mask]))
            bin_counts.append(count)
            
    if len(bin_centers) < 5:
        return None
        
    bin_centers = np.array(bin_centers)
    bin_means = np.array(bin_means)
    bin_counts = np.array(bin_counts)
    
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

def fit_exponential_binned(bin_centers, bin_means, bin_counts):
    if len(bin_centers) < 5:
        return None
    try:
        popt, pcov = curve_fit(
            exp_decay,
            bin_centers,
            bin_means,
            p0=[0.5, 2000, 0],
            sigma=1.0 / np.sqrt(bin_counts),
            bounds=([0, 100, -1], [2, 20000, 1]),
            maxfev=10000,
        )

        weights = bin_counts.astype(float)
        predicted = exp_decay(bin_centers, *popt)
        weighted_mean = np.average(bin_means, weights=weights)
        ss_res = np.sum(weights * (bin_means - predicted) ** 2)
        ss_tot = np.sum(weights * (bin_means - weighted_mean) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

        return {
            'lambda_km': float(popt[1]),
            'r_squared': float(r2),
            'n_pairs': int(np.sum(bin_counts)),
            'n_bins': int(len(bin_centers)),
        }
    except Exception:
        return None

def _init_bin_accumulators():
    return {
        'bin_counts': np.zeros(N_BINS, dtype=np.int64),
        'bin_sums': np.zeros(N_BINS, dtype=np.float64),
    }

def _accumulate_binned(distances, values, bin_edges, bin_counts, bin_sums):
    if distances.size == 0:
        return
    idx = np.searchsorted(bin_edges, distances, side='right') - 1
    valid = (idx >= 0) & (idx < N_BINS) & np.isfinite(values)
    if not np.any(valid):
        return
    idx = idx[valid]
    vals = values[valid]
    bin_counts += np.bincount(idx, minlength=N_BINS)
    bin_sums += np.bincount(idx, weights=vals, minlength=N_BINS)

def _accumulate_daily_means(years, doys, values, daily_sum, daily_count):
    valid = np.isfinite(values)
    if not np.any(valid):
        return
    keys = years[valid].astype(np.int64) * 1000 + doys[valid].astype(np.int64)
    vals = values[valid].astype(float)

    uniq, inv = np.unique(keys, return_inverse=True)
    sums = np.bincount(inv, weights=vals)
    counts = np.bincount(inv)
    for k, s, c in zip(uniq.tolist(), sums.tolist(), counts.tolist()):
        daily_sum[k] = daily_sum.get(k, 0.0) + float(s)
        daily_count[k] = daily_count.get(k, 0) + int(c)

def compute_cycle_phase_correlation(doys, coherences, period):
    phases = (doys % period) / period * 2 * np.pi
    sin_corr = np.corrcoef(coherences, np.sin(phases))[0, 1]
    cos_corr = np.corrcoef(coherences, np.cos(phases))[0, 1]
    return np.sqrt(sin_corr**2 + cos_corr**2)

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

def process_metric(csv_file, target_metric):
    """Process a single metric to save memory."""
    print_status(f"\nProcessing metric: {target_metric}...")

    bin_edges = np.logspace(np.log10(MIN_DISTANCE_KM), np.log10(MAX_DISTANCE_KM), N_BINS + 1)
    bin_centers_full = (bin_edges[:-1] + bin_edges[1:]) / 2
    rng = np.random.default_rng(42)

    accum_real = {c: _init_bin_accumulators() for c in COHERENCE_TYPES}
    accum_shuf = {c: _init_bin_accumulators() for c in COHERENCE_TYPES}
    daily_sum = {c: {} for c in COHERENCE_TYPES}
    daily_count = {c: {} for c in COHERENCE_TYPES}

    chunk_size = 5_000_000
    usecols_with_phase = ['metric', 'year', 'doy', 'distance_km', 'coherence', 'phase_alignment']
    usecols_no_phase = ['metric', 'year', 'doy', 'distance_km', 'coherence']
    try:
        reader = pd.read_csv(csv_file, chunksize=chunk_size, usecols=usecols_with_phase)
        phase_available = True
    except ValueError:
        reader = pd.read_csv(csv_file, chunksize=chunk_size, usecols=usecols_no_phase)
        phase_available = False

    total_rows = 0
    for chunk in reader:
        if 'metric' not in chunk.columns:
            continue

        mc = chunk[chunk['metric'].str.endswith(target_metric)]
        if mc.empty:
            continue

        dist = mc['distance_km'].to_numpy()
        years = mc['year'].to_numpy()
        doys = mc['doy'].to_numpy()

        msc = mc['coherence'].to_numpy()
        _accumulate_daily_means(years, doys, msc, daily_sum['msc'], daily_count['msc'])
        _accumulate_binned(
            dist,
            msc,
            bin_edges,
            accum_real['msc']['bin_counts'],
            accum_real['msc']['bin_sums'],
        )
        msc_valid = np.isfinite(msc)
        if np.any(msc_valid):
            msc_shuf = msc[msc_valid].astype(float).copy()
            rng.shuffle(msc_shuf)
            _accumulate_binned(
                dist[msc_valid],
                msc_shuf,
                bin_edges,
                accum_shuf['msc']['bin_counts'],
                accum_shuf['msc']['bin_sums'],
            )

        if phase_available and 'phase_alignment' in mc.columns:
            phase = mc['phase_alignment'].to_numpy()
            _accumulate_daily_means(
                years,
                doys,
                phase,
                daily_sum['phase_alignment'],
                daily_count['phase_alignment'],
            )
            _accumulate_binned(
                dist,
                phase,
                bin_edges,
                accum_real['phase_alignment']['bin_counts'],
                accum_real['phase_alignment']['bin_sums'],
            )
            phase_valid = np.isfinite(phase)
            if np.any(phase_valid):
                phase_shuf = phase[phase_valid].astype(float).copy()
                rng.shuffle(phase_shuf)
                _accumulate_binned(
                    dist[phase_valid],
                    phase_shuf,
                    bin_edges,
                    accum_shuf['phase_alignment']['bin_counts'],
                    accum_shuf['phase_alignment']['bin_sums'],
                )

        total_rows += len(mc)

    print_status(f"  Loaded {total_rows/1e6:.1f}M rows total.")
    
    results = {'solar': {}, 'lunar': {}, 'shuffle': {}}
    for coh_type in COHERENCE_TYPES:
        mode = f"{target_metric}/{coh_type}"
        print_status(f"  Analyzing {mode}...")

        # 1. Daily Means for Cycles
        d_keys = list(daily_sum[coh_type].keys())
        if len(d_keys) > 30:
            doys = np.array([int(k % 1000) for k in d_keys], dtype=float)
            vals = np.array(
                [daily_sum[coh_type][k] / max(daily_count[coh_type].get(k, 0), 1) for k in d_keys],
                dtype=float,
            )

            r_solar = compute_cycle_phase_correlation(doys, vals, SOLAR_ROTATION_DAYS)
            results['solar'][coh_type] = {
                'r': float(r_solar),
                'status': "NULL " if r_solar < 0.1 else "DETECTED ",
            }

            r_lunar = compute_cycle_phase_correlation(doys, vals, LUNAR_MONTH_DAYS)
            results['lunar'][coh_type] = {
                'r': float(r_lunar),
                'status': "NULL " if r_lunar < 0.1 else "DETECTED ",
            }

        # 2. Shuffle Test (streaming binned)
        real_counts = accum_real[coh_type]['bin_counts']
        real_sums = accum_real[coh_type]['bin_sums']
        shuf_counts = accum_shuf[coh_type]['bin_counts']
        shuf_sums = accum_shuf[coh_type]['bin_sums']

        real_valid = real_counts >= MIN_BIN_COUNT
        shuf_valid = shuf_counts >= MIN_BIN_COUNT

        real_res = None
        shuf_res = None

        if np.sum(real_valid) >= 5:
            real_centers = bin_centers_full[real_valid]
            real_means = (real_sums[real_valid] / real_counts[real_valid]).astype(float)
            real_res = fit_exponential_binned(real_centers, real_means, real_counts[real_valid].astype(float))

        if np.sum(shuf_valid) >= 5:
            shuf_centers = bin_centers_full[shuf_valid]
            shuf_means = (shuf_sums[shuf_valid] / shuf_counts[shuf_valid]).astype(float)
            shuf_res = fit_exponential_binned(shuf_centers, shuf_means, shuf_counts[shuf_valid].astype(float))

        if real_res and shuf_res:
            status = "PASS " if shuf_res['r_squared'] < 0.3 else "FAIL "
            results['shuffle'][coh_type] = {
                'real_r2': float(real_res['r_squared']),
                'shuffled_r2': float(shuf_res['r_squared']),
                'status': status,
            }

        gc.collect()
        
    return results

def run_null_tests(csv_file):
    """Run null tests for a single CSV file."""
    if not csv_file.exists():
        print_status(f"  WARNING: {csv_file.name} not found", level="WARNING")
        return None
    
    final_results = {'solar': {}, 'lunar': {}, 'shuffle': {}}
    
    for metric in METRICS:
        metric_results = process_metric(csv_file, metric)
        
        for c in COHERENCE_TYPES:
            mode = f"{metric}/{c}"
            if c in metric_results['solar']:
                final_results['solar'][mode] = metric_results['solar'][c]
                print_status(f"    Solar {mode}: r={metric_results['solar'][c]['r']:.3f}")
            if c in metric_results['lunar']:
                final_results['lunar'][mode] = metric_results['lunar'][c]
                print_status(f"    Lunar {mode}: r={metric_results['lunar'][c]['r']:.3f}")
            if c in metric_results['shuffle']:
                final_results['shuffle'][mode] = metric_results['shuffle'][c]
                print_status(f"    Shuffle {mode}: Real R²={metric_results['shuffle'][c]['real_r2']:.3f}, Shuffled R²={metric_results['shuffle'][c]['shuffled_r2']:.3f}")
        
        gc.collect()
    
    return final_results

def main():
    import argparse
    parser = argparse.ArgumentParser(description='TEP-GNSS-RINEX Step 2.4b: Null Tests')
    parser.add_argument('--filter', type=str, default='all',
                        help='Filter: "all", "none", "optimal_100_metadata.json", or "dynamic_50_metadata.json"')
    args = parser.parse_args()
    
    print_status('STEP 2.4b: Null Tests (COMPREHENSIVE)', level="TITLE")
    print_status('Processing: 3 filters x 4 modes x 3 metrics x 2 coherence types')
    
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
            csv_file = RESULTS_DIR / f"step_2_0_pairs_{mode_name}_{filter_key}.csv"
            print_status(f"\n  MODE: {mode_name.upper()}")
            
            result = run_null_tests(csv_file)
            if result:
                all_mode_results[mode_name] = result
        
        all_filter_results[filter_key] = {
            'filter_name': filter_name,
            'filter_metadata': filter_meta,
            'results_by_mode': all_mode_results
        }
        
        # Save per-filter JSON
        output = {
            'step': '2.4b', 'name': 'null_tests',
            'filter': filter_name,
            'timestamp': datetime.now().isoformat(),
            'results_by_mode': all_mode_results
        }
        with open(RESULTS_DIR / f"step_2_4_null_tests_{filter_key}.json", 'w') as f:
            json.dump(output, f, indent=2, default=str)
        print_status(f"  Saved: step_2_4_null_tests_{filter_key}.json", level="SUCCESS")
    
    # Save summary
    summary = {
        'step': '2.4b', 'name': 'null_tests_comprehensive',
        'timestamp': datetime.now().isoformat(),
        'filters_tested': [f[0] for f in filters],
        'modes_tested': list(PROCESSING_MODES.keys()),
        'results_by_filter': all_filter_results
    }
    with open(RESULTS_DIR / "step_2_4_null_tests_summary.json", 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    print_status(f"COMPLETE: {len(filters)} filters x {len(PROCESSING_MODES)} modes", level="SUCCESS")


if __name__ == '__main__':
    main()
