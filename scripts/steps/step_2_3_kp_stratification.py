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
DEFAULT_KP_THRESHOLDS = "3.0,4.0,5.0"

RESULTS_DIR = PROJECT_ROOT / "results" / "outputs"
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

# Global station filter set
GOOD_STATIONS_SET = None

# ALL FOUR PROCESSING MODES - consistent with other scripts
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


def _parse_kp_thresholds(raw: str) -> list[float]:
    raw = (raw or "").strip()
    if not raw:
        raw = DEFAULT_KP_THRESHOLDS
    parts = [p.strip() for p in raw.split(',') if p.strip()]
    thresholds = []
    for p in parts:
        thresholds.append(float(p))
    thresholds = sorted(set(thresholds))
    return thresholds


def _kp_tag(threshold: float) -> str:
    return f"_kp{threshold:.1f}".replace('.', 'p')

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

        weights = bin_counts.astype(float)
        predicted = exp_decay(bin_centers, *popt)
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

def process_metric(csv_file, target_metric, kp_values, kp_thresholds, station_set=None):
    print_status(f"\nProcessing metric: {target_metric}...")

    bin_edges = np.logspace(np.log10(MIN_DISTANCE_KM), np.log10(MAX_DISTANCE_KM), N_BINS + 1)

    binners_by_thr = {
        thr: {
            c: {
                'quiet': _init_binner(),
                'storm': _init_binner(),
            }
            for c in COHERENCE_TYPES
        }
        for thr in kp_thresholds
    }

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

        years = mc['year'].to_numpy(dtype=np.int64, copy=False)
        doys = mc['doy'].to_numpy(dtype=np.int64, copy=False)
        dist = mc['distance_km'].to_numpy(dtype=np.float64, copy=False)
        msc = mc['coherence'].to_numpy(dtype=np.float64, copy=False)
        phase = mc['phase_alignment'].to_numpy(dtype=np.float64, copy=False) if 'phase_alignment' in mc.columns else None

        year_idx = years - YEARS[0]
        valid_day = (year_idx >= 0) & (year_idx < len(YEARS)) & (doys >= 1) & (doys <= 366)
        kp = np.full(years.shape[0], np.nan, dtype=np.float64)
        if np.any(valid_day):
            kp[valid_day] = kp_values[year_idx[valid_day], doys[valid_day]]
        has_kp = np.isfinite(kp)

        bin_idx = np.searchsorted(bin_edges, dist, side='right') - 1
        in_range = (bin_idx >= 0) & (bin_idx < N_BINS)

        valid_msc = in_range & has_kp & ~np.isnan(msc)
        if np.any(valid_msc):
            for thr in kp_thresholds:
                q_idx = valid_msc & (kp < thr)
                s_idx = valid_msc & (kp >= thr)
                if np.any(q_idx):
                    _accumulate_binner(binners_by_thr[thr]['msc']['quiet'], bin_idx[q_idx], msc[q_idx])
                if np.any(s_idx):
                    _accumulate_binner(binners_by_thr[thr]['msc']['storm'], bin_idx[s_idx], msc[s_idx])

        if phase is not None:
            valid_phase = in_range & has_kp & ~np.isnan(phase)
            if np.any(valid_phase):
                for thr in kp_thresholds:
                    q_idx = valid_phase & (kp < thr)
                    s_idx = valid_phase & (kp >= thr)
                    if np.any(q_idx):
                        _accumulate_binner(binners_by_thr[thr]['phase_alignment']['quiet'], bin_idx[q_idx], phase[q_idx])
                    if np.any(s_idx):
                        _accumulate_binner(binners_by_thr[thr]['phase_alignment']['storm'], bin_idx[s_idx], phase[s_idx])

    results_by_thr = {thr: {} for thr in kp_thresholds}
    for thr in kp_thresholds:
        for c in COHERENCE_TYPES:
            mode = f"{target_metric}/{c}"
            results_by_thr[thr][mode] = {}

            qc, qm, qn = _finalize_binner(binners_by_thr[thr][c]['quiet'], bin_edges)
            q_res = fit_exponential_binned(qc, qm, qn) if qc is not None else None
            results_by_thr[thr][mode]['quiet'] = q_res

            sc, sm, sn = _finalize_binner(binners_by_thr[thr][c]['storm'], bin_edges)
            s_res = fit_exponential_binned(sc, sm, sn) if sc is not None else None
            results_by_thr[thr][mode]['storm'] = s_res

            all_sum = binners_by_thr[thr][c]['quiet']['sum'] + binners_by_thr[thr][c]['storm']['sum']
            all_count = binners_by_thr[thr][c]['quiet']['count'] + binners_by_thr[thr][c]['storm']['count']
            all_binner = {'sum': all_sum, 'count': all_count}
            ac, am, an = _finalize_binner(all_binner, bin_edges)
            all_res = fit_exponential_binned(ac, am, an) if ac is not None else None
            results_by_thr[thr][mode]['all'] = all_res

            print_status(f"  {mode} (Kp≥{thr:g}): Q={q_res['lambda_km'] if q_res else 'N/A'}, S={s_res['lambda_km'] if s_res else 'N/A'}")

    return results_by_thr


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
    parser.add_argument('--kp-thresholds', '--kp-threshold', dest='kp_thresholds', type=str, default=DEFAULT_KP_THRESHOLDS,
                        help='Comma-separated Kp thresholds for storm classification (storm if Kp ≥ threshold).')
    parser.add_argument('--tag', type=str, default='',
                        help='Optional tag appended to output filenames (for non-overwriting reruns).')
    args = parser.parse_args()
    
    print_status('STEP 2.3a: Kp Geomagnetic Stratification (COMPREHENSIVE)', level="TITLE")
    kp_thresholds = _parse_kp_thresholds(args.kp_thresholds)
    run_tag = _sanitize_tag(args.tag)
    
    # Select filters
    if args.filter == 'all':
        filters = STATION_FILTERS
    else:
        matching = [f for f in STATION_FILTERS if f[1] == args.filter or f[0] == args.filter]
        filters = matching if matching else [('ALL_STATIONS', 'none')]

    print_status(f"Processing: {len(filters)} filters x {len(PROCESSING_MODES)} modes x {len(METRICS)} metrics x {len(COHERENCE_TYPES)} coherence types x {len(kp_thresholds)} Kp thresholds")
    
    # Load Kp data ONCE
    kp_data = load_kp_index()

    kp_values = np.full((len(YEARS), 367), np.nan, dtype=np.float64)
    for key, kp in kp_data.items():
        year, doy = map(int, key.split('_'))
        if year not in YEARS or doy < 1 or doy > 366:
            continue
        y_idx = year - YEARS[0]
        kp_values[y_idx, doy] = float(kp)

    all_kp_vals = kp_values[np.isfinite(kp_values)]
    print_status(f"  Loaded {all_kp_vals.size} Kp-tagged days for analysis")
    for thr in kp_thresholds:
        q = int(np.sum(all_kp_vals < thr))
        s = int(np.sum(all_kp_vals >= thr))
        print_status(f"  Threshold Kp≥{thr:g}: Quiet={q} ({100*q/(q+s):.1f}%), Storm={s} ({100*s/(q+s):.1f}%)")
    
    all_filter_results_by_thr = {thr: {} for thr in kp_thresholds}
    
    for filter_name, filter_config in filters:
        filter_key = filter_name.lower()
        print_status(f"\n{'#'*80}")
        print_status(f"# FILTER: {filter_name}")
        print_status(f"{'#'*80}")
        
        station_set, filter_meta = load_station_filter(filter_config)
        if station_set:
            print_status(f"  Using {len(station_set)} stations")
        
        all_mode_results_by_thr = {thr: {} for thr in kp_thresholds}
        
        for mode_name, mode_config in PROCESSING_MODES.items():
            # Construct filename dynamically: step_2_0_pairs_{mode}_{filter}.csv
            csv_filename = f"step_2_0_pairs_{mode_name}_{filter_key}.csv"
            csv_file = RESULTS_DIR / csv_filename
            
            print_status(f"\n  MODE: {mode_name.upper()} ({mode_config['description']})")
            
            if not csv_file.exists():
                print_status(f"    WARNING: {csv_file.name} not found, skipping", level="WARNING")
                continue
            
            mode_results_by_thr = {thr: {} for thr in kp_thresholds}
            for metric in METRICS:
                res_by_thr = process_metric(csv_file, metric, kp_values, kp_thresholds, station_set=station_set)
                for thr in kp_thresholds:
                    mode_results_by_thr[thr].update(res_by_thr[thr])
                gc.collect()

            for thr in kp_thresholds:
                all_mode_results_by_thr[thr][mode_name] = mode_results_by_thr[thr]
            
            # Print mode summary
            for thr in kp_thresholds:
                print_status(f"\n    Threshold summary (Kp≥{thr:g}):")
                for key, val in mode_results_by_thr[thr].items():
                    if val.get('quiet') and val.get('storm'):
                        q_lam = val['quiet']['lambda_km']
                        s_lam = val['storm']['lambda_km']
                        delta = 100 * (s_lam - q_lam) / q_lam
                        print_status(f"      {key}: Q={q_lam:.0f}km, S={s_lam:.0f}km, Δ={delta:+.1f}%")

        for thr in kp_thresholds:
            all_filter_results_by_thr[thr][filter_key] = {
                'filter_name': filter_name,
                'filter_metadata': filter_meta,
                'results_by_mode': all_mode_results_by_thr[thr]
            }

            thr_tag = _kp_tag(thr)
            q = int(np.sum(all_kp_vals < thr))
            s = int(np.sum(all_kp_vals >= thr))
            output = {
                'step': '2.3a',
                'name': 'kp_stratification',
                'filter': filter_name,
                'timestamp': datetime.now().isoformat(),
                'kp_threshold': thr,
                'quiet_days': q,
                'storm_days': s,
                'results_by_mode': all_mode_results_by_thr[thr]
            }

            # Always write a threshold-tagged file
            with open(RESULTS_DIR / f"step_2_3_kp_stratification_{filter_key}{thr_tag}{run_tag}.json", 'w') as f:
                json.dump(output, f, indent=2, default=str)
            print_status(f"  Saved: step_2_3_kp_stratification_{filter_key}{thr_tag}{run_tag}.json", level="SUCCESS")

            # Also write legacy (unsuffixed) file for Kp≥3.0 for backward compatibility
            if abs(thr - KP_THRESHOLD) < 1e-9:
                with open(RESULTS_DIR / f"step_2_3_kp_stratification_{filter_key}{run_tag}.json", 'w') as f:
                    json.dump(output, f, indent=2, default=str)
                print_status(f"  Saved: step_2_3_kp_stratification_{filter_key}{run_tag}.json", level="SUCCESS")
    
    # Save summary (one per threshold)
    for thr in kp_thresholds:
        thr_tag = _kp_tag(thr)
        summary = {
            'step': '2.3a',
            'name': 'kp_stratification_comprehensive',
            'timestamp': datetime.now().isoformat(),
            'kp_threshold': thr,
            'filters_tested': [f[0] for f in filters],
            'modes_tested': list(PROCESSING_MODES.keys()),
            'results_by_filter': all_filter_results_by_thr[thr]
        }
        with open(RESULTS_DIR / f"step_2_3_kp_stratification_summary{thr_tag}{run_tag}.json", 'w') as f:
            json.dump(summary, f, indent=2, default=str)

        if abs(thr - KP_THRESHOLD) < 1e-9:
            with open(RESULTS_DIR / f"step_2_3_kp_stratification_summary{run_tag}.json", 'w') as f:
                json.dump(summary, f, indent=2, default=str)
    
    print_status(f"COMPLETE: {len(filters)} filters × {len(PROCESSING_MODES)} modes × {len(kp_thresholds)} thresholds", level="SUCCESS")

if __name__ == '__main__':
    main()
