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
import logging
import re
import gc
from scipy import stats
from datetime import datetime
from scripts.utils.tep_utils import bin_and_fit, ecef_to_lla, load_station_map
from scripts.utils.logger import TEPLogger, TEPFileFormatter, set_step_logger, print_status

LOG_FILE_PATH = PROJECT_ROOT / "logs" / "step_2_1_elevation_analysis.log"
logger = TEPLogger(
    name="step_2_1_elevation_analysis",
    level="INFO",
    log_file_path=LOG_FILE_PATH
)
set_step_logger(logger)

# Parameters
RESULTS_DIR = PROJECT_ROOT / "results" / "outputs"
FIGURES_DIR = PROJECT_ROOT / "results" / "figures"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"

METRICS = ['clock_bias', 'pos_jitter', 'clock_drift']
COHERENCE_TYPES = ['msc', 'phase_alignment']
QUINTILES = ['Q1', 'Q2', 'Q3', 'Q4', 'Q5']

LAT_BAND_DEG = None
REQUIRE_SAME_LAT_BAND = False
SHORT_RANGE_MAX_KM = 200.0
OUTPUT_TAG = None
OUTPUT_SUFFIX = ""
TAG_LOG_HANDLER = None

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

def get_station_quintiles(coords):
    """Map station to quintile (Q1=Low, Q5=High)."""

    return get_station_quintiles_lat_controlled(coords, lat_band_deg=None)

def get_station_quintiles_lat_controlled(coords, lat_band_deg=None):
    stas = []
    for s, xyz in coords.items():
        lat, _, alt = ecef_to_lla(*xyz)
        if lat_band_deg and lat_band_deg > 0:
            band = int(np.floor((lat + 90.0) / float(lat_band_deg)))
        else:
            band = 0
        stas.append((s, float(lat), float(alt), band))
    
    q_map = {}
    band_map = {}
    q_stats = {}

    by_band = {}
    for s, lat, alt, band in stas:
        by_band.setdefault(band, []).append((s, lat, alt))

    for band, band_stas in by_band.items():
        band_stas.sort(key=lambda x: x[2])
        n = len(band_stas)
        if n < 5:
            continue
        q_size = n // 5

        for i in range(5):
            q_name = f"Q{i+1}"
            start = i * q_size
            end = (i + 1) * q_size if i < 4 else n
            subset = band_stas[start:end]
            for s, lat, alt in subset:
                q_map[s] = q_name
                band_map[s] = band

    for q in QUINTILES:
        alts = []
        lats = []
        for s, lat, alt, _band in stas:
            if q_map.get(s) == q:
                alts.append(alt)
                lats.append(lat)
        if not alts:
            q_stats[q] = {'count': 0}
            continue

        alts_arr = np.asarray(alts, dtype=float)
        lats_arr = np.asarray(lats, dtype=float)
        q_stats[q] = {
            'min': float(np.min(alts_arr)),
            'max': float(np.max(alts_arr)),
            'mean': float(np.mean(alts_arr)),
            'median': float(np.median(alts_arr)),
            'count': int(len(alts_arr)),
            'lat_mean': float(np.mean(lats_arr)),
            'lat_median': float(np.median(lats_arr))
        }
            
    return q_map, q_stats, band_map

def _weighted_linregress(x, y, w):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    w = np.asarray(w, dtype=float)

    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    x = x[valid]
    y = y[valid]
    w = w[valid]

    n = int(len(x))
    if n < 3:
        return {
            'x': 'altitude_mean_m',
            'y': 'lambda_km',
            'n_quintiles_used': n,
            'trend_available': False
        }

    w_sum = float(np.sum(w))
    x_bar = float(np.sum(w * x) / w_sum)
    y_bar = float(np.sum(w * y) / w_sum)
    sxx = float(np.sum(w * (x - x_bar) ** 2))
    syy = float(np.sum(w * (y - y_bar) ** 2))
    sxy = float(np.sum(w * (x - x_bar) * (y - y_bar)))
    if sxx <= 0:
        return {
            'x': 'altitude_mean_m',
            'y': 'lambda_km',
            'n_quintiles_used': n,
            'trend_available': False
        }

    slope = float(sxy / sxx)
    intercept = float(y_bar - slope * x_bar)
    r_value = float(sxy / np.sqrt(sxx * syy)) if syy > 0 else 0.0

    y_hat = slope * x + intercept
    resid = y - y_hat
    sigma2 = float(np.sum(w * resid ** 2) / max(1, n - 2))
    stderr_slope = float(np.sqrt(sigma2 / sxx)) if sigma2 >= 0 else None
    if stderr_slope and np.isfinite(stderr_slope) and stderr_slope > 0:
        t_stat = slope / stderr_slope
        p_value = float(2.0 * stats.t.sf(np.abs(t_stat), df=n - 2))
    else:
        p_value = None

    return {
        'x': 'altitude_mean_m',
        'y': 'lambda_km',
        'n_quintiles_used': n,
        'slope_km_per_m': slope,
        'intercept_km': intercept,
        'r_value': r_value,
        'p_value': p_value,
        'stderr_slope': stderr_slope
    }

def _trend_from_quintiles(q_stats, q_results, y_key, y_label):
    altitudes = []
    yvals = []
    for q in QUINTILES:
        altitudes.append(q_stats.get(q, {}).get('mean', np.nan))
        yvals.append(q_results.get(q, {}).get(y_key, np.nan))
    altitudes = np.asarray(altitudes, dtype=float)
    yvals = np.asarray(yvals, dtype=float)

    valid = np.isfinite(altitudes) & np.isfinite(yvals)
    if int(np.sum(valid)) < 3:
        return {
            'x': 'altitude_mean_m',
            'y': y_label,
            'n_quintiles_used': int(np.sum(valid)),
            'trend_available': False
        }

    lr = stats.linregress(altitudes[valid], yvals[valid])
    return {
        'x': 'altitude_mean_m',
        'y': y_label,
        'n_quintiles_used': int(np.sum(valid)),
        'slope_per_m': float(lr.slope),
        'intercept': float(lr.intercept),
        'r_value': float(lr.rvalue),
        'p_value': float(lr.pvalue),
        'stderr_slope': float(lr.stderr) if lr.stderr is not None else None
    }

def process_metric(metric):
    print_status(f"\nProcessing {metric}...")
    
    data = {c: {q: {'dist': [], 'coh': []} for q in QUINTILES} for c in COHERENCE_TYPES}
    short_range = {c: {q: {'sum': 0.0, 'count': 0} for q in QUINTILES} for c in COHERENCE_TYPES}
    chunk_size = 10_000_000  # Larger chunks for high-RAM systems
    
    # Load Station Map
    coords = load_station_map(PROCESSED_DIR)
    if not coords: return None
    
    # Apply station filter if set
    if GOOD_STATIONS_SET:
        coords = {k: v for k, v in coords.items() if k in GOOD_STATIONS_SET}
        print_status(f"  Filtered to {len(coords)} stations")
    
    q_map, q_stats, band_map = get_station_quintiles_lat_controlled(coords, lat_band_deg=LAT_BAND_DEG)
    
    for k, v in q_stats.items():
        if v.get('count', 0) > 0:
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
        if REQUIRE_SAME_LAT_BAND:
            b1 = [band_map.get(s, None) for s in s1]
            b2 = [band_map.get(s, None) for s in s2]
        else:
            b1 = None
            b2 = None
        
        # Filter: Both in same quintile
        # Create boolean masks for each Q
        for q in QUINTILES:
            # Vectorized check? List comprehension is reasonable for 5M items if simple string compare
            if REQUIRE_SAME_LAT_BAND:
                mask = np.array([qa == q and qb == q and ba == bb for qa, qb, ba, bb in zip(q1, q2, b1, b2)])
            else:
                mask = np.array([qa == q and qb == q for qa, qb in zip(q1, q2)])
            
            if np.any(mask):
                for coh_type in COHERENCE_TYPES:
                    vals = msc if coh_type == 'msc' else phase
                    if vals is None: continue
                    
                    valid = ~np.isnan(vals) & mask
                    if np.any(valid):
                        data[coh_type][q]['dist'].extend(dist[valid].tolist())
                        data[coh_type][q]['coh'].extend(vals[valid].tolist())

                        sr_valid = valid & (dist <= SHORT_RANGE_MAX_KM)
                        if np.any(sr_valid):
                            short_range[coh_type][q]['sum'] += float(np.nansum(vals[sr_valid]))
                            short_range[coh_type][q]['count'] += int(np.sum(sr_valid))
                        
    # Fit
    results = {}
    for coh_type in COHERENCE_TYPES:
        mode = f"{metric}/{coh_type}"
        results[coh_type] = {}
        print_status(f"  Fitting {mode}...")
        
        for q in QUINTILES:
            d = data[coh_type][q]
            res = bin_and_fit(d['dist'], d['coh'])
            if res is not None:
                res = dict(res)
                res['n_samples'] = int(len(d['dist']))
            else:
                res = {
                    'lambda_km': np.nan,
                    'lambda_err': np.nan,
                    'r_squared': np.nan,
                    'n_pairs': 0,
                    'n_samples': int(len(d['dist'])),
                    'fit_failed': True
                }

            sr_count = short_range[coh_type][q]['count']
            sr_sum = short_range[coh_type][q]['sum']
            res['short_range_max_km'] = float(SHORT_RANGE_MAX_KM)
            res['short_range_n'] = int(sr_count)
            res['short_range_mean'] = float(sr_sum / sr_count) if sr_count > 0 else np.nan

            results[coh_type][q] = res
            if res:
                if 'lambda_km' in res and 'r_squared' in res:
                    print_status(f"    {q:<4}: λ={res['lambda_km']:.0f}km, R²={res['r_squared']:.2f}, N={res.get('n_pairs', 0):,}")
                else:
                    print_status(f"    {q:<4}: Fit Failed (N={res.get('n_samples', 0):,})", level="WARNING")
            else:
                print_status(f"    {q:<4}: Fit Failed", level="WARNING")

        altitudes = np.asarray([q_stats.get(q, {}).get('mean', np.nan) for q in QUINTILES], dtype=float)
        lambdas = np.asarray([results[coh_type][q].get('lambda_km', np.nan) for q in QUINTILES], dtype=float)
        valid = np.isfinite(altitudes) & np.isfinite(lambdas)
        if int(np.sum(valid)) >= 3:
            lr = stats.linregress(altitudes[valid], lambdas[valid])
            results[coh_type]['trend'] = {
                'x': 'altitude_mean_m',
                'y': 'lambda_km',
                'n_quintiles_used': int(np.sum(valid)),
                'slope_km_per_m': float(lr.slope),
                'intercept_km': float(lr.intercept),
                'r_value': float(lr.rvalue),
                'p_value': float(lr.pvalue),
                'stderr_slope': float(lr.stderr) if lr.stderr is not None else None
            }
        else:
            results[coh_type]['trend'] = {
                'x': 'altitude_mean_m',
                'y': 'lambda_km',
                'n_quintiles_used': int(np.sum(valid)),
                'trend_available': False
            }

        weights_n = np.asarray([results[coh_type][q].get('n_pairs', 0) for q in QUINTILES], dtype=float)
        results[coh_type]['trend_weighted_n_pairs'] = _weighted_linregress(altitudes, lambdas, weights_n)

        invvar = []
        for q in QUINTILES:
            le = results[coh_type][q].get('lambda_err', np.nan)
            if le and np.isfinite(le) and le > 0:
                invvar.append(1.0 / (le ** 2))
            else:
                invvar.append(0.0)
        weights_iv = np.asarray(invvar, dtype=float)
        results[coh_type]['trend_weighted_invvar'] = _weighted_linregress(altitudes, lambdas, weights_iv)

        results[coh_type]['trend_amplitude'] = _trend_from_quintiles(q_stats, results[coh_type], 'amplitude', 'amplitude')
        results[coh_type]['trend_offset'] = _trend_from_quintiles(q_stats, results[coh_type], 'offset', 'offset')
        results[coh_type]['trend_short_range_mean'] = _trend_from_quintiles(q_stats, results[coh_type], 'short_range_mean', 'short_range_mean')

        results[coh_type]['stratification'] = {
            'lat_band_deg': float(LAT_BAND_DEG) if LAT_BAND_DEG else None,
            'require_same_lat_band': bool(REQUIRE_SAME_LAT_BAND),
            'short_range_max_km': float(SHORT_RANGE_MAX_KM)
        }

        results[coh_type]['altitude_quintiles'] = q_stats
                
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
            lambdas = [res.get(q, {}).get('lambda_km', np.nan) for q in QUINTILES]
            errs = [res.get(q, {}).get('lambda_err', 0) for q in QUINTILES]
            
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
    fig_path = FIGURES_DIR / f"step_2_1_elevation_analysis{OUTPUT_SUFFIX}.png"
    plt.savefig(fig_path)
    print_status(f"Saved figure to {fig_path}", level="SUCCESS")

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

def _reset_master_log():
    try:
        LOG_FILE_PATH.parent.mkdir(parents=True, exist_ok=True)
        with open(LOG_FILE_PATH, 'w') as f:
            f.write("")
    except Exception:
        pass

def _configure_tag_log(tag: str):
    global TAG_LOG_HANDLER
    py_logger = logger.logger

    if TAG_LOG_HANDLER is not None:
        try:
            py_logger.removeHandler(TAG_LOG_HANDLER)
        except Exception:
            pass
        try:
            TAG_LOG_HANDLER.close()
        except Exception:
            pass
        TAG_LOG_HANDLER = None

    if not tag:
        return

    tag_path = LOG_FILE_PATH.parent / f"step_2_1_elevation_analysis_{tag}.log"
    tag_path.parent.mkdir(parents=True, exist_ok=True)

    fh = logging.FileHandler(tag_path, mode='w', encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(TEPFileFormatter())
    py_logger.addHandler(fh)
    TAG_LOG_HANDLER = fh

def _set_run_config(lat_band_deg: float, short_range_km: float, output_tag: str):
    global LAT_BAND_DEG
    global REQUIRE_SAME_LAT_BAND
    global SHORT_RANGE_MAX_KM
    global OUTPUT_TAG
    global OUTPUT_SUFFIX

    LAT_BAND_DEG = float(lat_band_deg) if lat_band_deg and lat_band_deg > 0 else None
    REQUIRE_SAME_LAT_BAND = bool(LAT_BAND_DEG is not None)
    SHORT_RANGE_MAX_KM = float(short_range_km)

    tag = (output_tag or "").strip()
    if tag:
        tag = re.sub(r"[^A-Za-z0-9_.-]+", "_", tag).strip("_")
    OUTPUT_TAG = tag if tag else None
    OUTPUT_SUFFIX = f"_{OUTPUT_TAG}" if OUTPUT_TAG else ""

def _run_analysis(filters):
    _configure_tag_log(OUTPUT_TAG)
    print_status("STEP 2.1b: Elevation Analysis (COMPREHENSIVE)", level="TITLE")
    print_status(
        f"Run tag: {OUTPUT_TAG if OUTPUT_TAG else '(none)'} | lat_band_deg={LAT_BAND_DEG if LAT_BAND_DEG else 0.0:g} | short_range_km={SHORT_RANGE_MAX_KM:g}"
    )
    print_status(
        f"Processing: {len(filters)} filters × {len(PROCESSING_MODES)} modes × {len(METRICS)} metrics × {len(COHERENCE_TYPES)} coherence types"
    )

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
            csv_file = _resolve_pairs_csv(mode_name, filter_key)
            print_status(f"\n{'='*60}")
            print_status(f"MODE: {mode_name.upper()} ({mode_config['description']})")
            print_status(f"{'='*60}")
            print_status("")

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

        output = {
            'step': '2.1b',
            'name': 'elevation_analysis',
            'filter': filter_name,
            'timestamp': datetime.now().isoformat(),
            'output_tag': OUTPUT_TAG,
            'modes_processed': list(all_mode_results.keys()),
            'results_by_mode': all_mode_results,
            'results': all_mode_results.get('baseline', {})
        }
        output_file = RESULTS_DIR / f"step_2_1_elevation_analysis_{filter_key}{OUTPUT_SUFFIX}.json"
        with open(output_file, 'w') as f:
            json.dump(output, f, indent=2)
        print_status(f"Saved: {output_file.name}", level="SUCCESS")

    summary = {
        'step': '2.1b',
        'name': 'elevation_analysis_comprehensive',
        'timestamp': datetime.now().isoformat(),
        'output_tag': OUTPUT_TAG,
        'filters_tested': [f[0] for f in filters],
        'modes_tested': list(PROCESSING_MODES.keys()),
        'results_by_filter': all_filter_results
    }
    with open(RESULTS_DIR / f"step_2_1_elevation_analysis_summary{OUTPUT_SUFFIX}.json", 'w') as f:
        json.dump(summary, f, indent=2)

    first_filter = list(all_filter_results.keys())[0] if all_filter_results else None
    if first_filter and 'baseline' in all_filter_results[first_filter].get('results_by_mode', {}):
        plot_summary(all_filter_results[first_filter]['results_by_mode']['baseline'])

    print_status(f"COMPLETE: {len(filters)} filters × {len(PROCESSING_MODES)} modes", level="SUCCESS")

    _configure_tag_log("")

def main():
    import argparse
    parser = argparse.ArgumentParser(description='TEP-GNSS-RINEX Step 2.1b: Elevation Analysis')
    parser.add_argument('--filter', type=str, default='all',
                        help='Filter: "all", "none", "optimal_100_metadata.json", or "dynamic_50_metadata.json"')
    parser.add_argument('--lat-band-deg', type=float, default=0.0,
                        help='If >0, compute altitude quintiles within latitude bands of this width (degrees) and only use pairs within the same band.')
    parser.add_argument('--short-range-km', type=float, default=200.0,
                        help='Maximum distance (km) for short-range coherence summary per quintile.')
    parser.add_argument('--output-tag', type=str, default='',
                        help='Optional suffix tag for output files (prevents overwriting previous runs).')
    parser.add_argument('--single', action='store_true',
                        help='Run a single configuration using the provided --lat-band-deg/--short-range-km settings. If omitted, runs the full recommended suite.')
    args = parser.parse_args()

    _reset_master_log()

    if args.filter == 'all':
        filters = STATION_FILTERS
    else:
        matching = [f for f in STATION_FILTERS if f[1] == args.filter]
        filters = matching if matching else [('ALL_STATIONS', 'none')]

    if args.single:
        _set_run_config(args.lat_band_deg, args.short_range_km, args.output_tag)
        _run_analysis(filters)
        return

    suite = [
        (0.0, 200.0, "global_sr200"),
        (10.0, 200.0, "lat10_sr200"),
        (20.0, 200.0, "lat20_sr200"),
        (10.0, 100.0, "lat10_sr100"),
        (20.0, 100.0, "lat20_sr100"),
    ]

    for lat_band_deg, short_range_km, tag in suite:
        _set_run_config(lat_band_deg, short_range_km, tag)
        _run_analysis(filters)

if __name__ == '__main__':
    main()
