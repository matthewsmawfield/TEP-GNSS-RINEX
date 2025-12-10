#!/usr/bin/env python3
"""
TEP-GNSS-RINEX - STEP 2.5: Orbital Velocity Coupling (COMPREHENSIVE)
=====================================================================
Validates the strongest TEP signature: Correlation between E-W/N-S anisotropy
and Earth's orbital velocity.

COMPREHENSIVE ANALYSIS:
- 3 Station Filters: ALL, OPTIMAL_100, DYNAMIC_50
- 3 Processing Modes: baseline, ionofree, multi_gnss
- 3 Metrics: clock_bias, pos_jitter, clock_drift
- 2 Coherence Types: msc, phase_alignment

Quality Control: Filters out lambda > 18,000 km (fitting failures).
Memory Optimized: Uses streaming binning to avoid OOM.
Outputs: Separate JSON per filter with monthly data for CMB analysis.

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
from scipy.optimize import curve_fit
from scipy import stats
from datetime import datetime
import math
import gc
import argparse

from scripts.utils.logger import print_status, TEPLogger, set_step_logger

# Initialize step logger for real-time output
LOGS_DIR = PROJECT_ROOT / "logs"
LOGS_DIR.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_2_5", log_file_path=LOGS_DIR / "step_2_5_orbital_coupling.log")
set_step_logger(logger)

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

# ============================================================================
# PROCESSING MODES - consistent with other Step 2.x scripts
# ============================================================================
PROCESSING_MODES = {
    'baseline': {
        'csv_file': 'step_2_0_pairs_baseline.csv',
        'description': 'SPP with Broadcast Ephemeris'
    },
    'ionofree': {
        'csv_file': 'step_2_0_pairs_ionofree.csv',
        'description': 'Dual-Freq Iono-Free with Precise Orbits'
    },
    'multi_gnss': {
        'csv_file': 'step_2_0_pairs_multi_gnss.csv',
        'description': 'Multi-GNSS (GPS+GLO+GAL+BDS)'
    }
}

# ============================================================================
# STATION FILTERS - consistent with Step 2.0/2.2
# ============================================================================
STATION_FILTERS = [
    ('ALL_STATIONS', 'none'),
    ('OPTIMAL_100', 'optimal_100_metadata.json'),
    ('DYNAMIC_50', 'dynamic_50_metadata.json')
]


def haversine_azimuth(lat1, lon1, lat2, lon2):
    """Compute distance and azimuth."""
    R = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat/2)**2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))
    dist = R * c
    
    y = math.sin(dlon) * math.cos(math.radians(lat2))
    x = math.cos(math.radians(lat1)) * math.sin(math.radians(lat2)) - math.sin(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.cos(dlon)
    az = math.degrees(math.atan2(y, x))
    return dist, (az + 360) % 360


def is_ew(az):
    """E-W direction: ±22.5° from East (90°) or West (270°)."""
    return (67.5 <= az < 112.5) or (247.5 <= az < 292.5)


def is_ns(az):
    """N-S direction: ±22.5° from North (0°) or South (180°)."""
    return (az < 22.5) or (157.5 <= az < 202.5) or (az >= 337.5)


def get_orbital_velocity(doy, year):
    """Get Earth's orbital velocity for a given day."""
    perihelion_doy = 3
    days_in_year = 366 if (year % 4 == 0) else 365
    theta = 2 * np.pi * (doy - perihelion_doy) / days_in_year
    v_mean = 29.78
    e = 0.0167
    return v_mean * (1 + e * np.cos(theta))


def ecef_to_lla(x, y, z):
    """Convert ECEF to lat/lon."""
    lon = math.atan2(y, x)
    p = math.sqrt(x**2 + y**2)
    lat = math.atan2(z, p * (1 - 0.00669438))
    v = 6378137 / math.sqrt(1 - 0.00669438 * math.sin(lat)**2)
    lat = math.atan2(z + 0.00669438 * v * math.sin(lat), p)
    return math.degrees(lat), math.degrees(lon)


def exp_decay(r, A, lam, C0):
    """Exponential decay function."""
    return A * np.exp(-r / lam) + C0


def fit_exponential_binned(bin_centers, bin_means, bin_counts):
    """Fit exponential decay to binned data."""
    try:
        if len(bin_centers) < 5:
            return None
        popt, pcov = curve_fit(
            exp_decay, bin_centers, bin_means, 
            p0=[0.5, 2000, 0],
            sigma=1.0/np.sqrt(bin_counts), 
            bounds=([0, 100, -1], [2, 20000, 1]), 
            maxfev=10000
        )
        
        # QC: Filter out fitting failures (hitting upper bound)
        if popt[1] > 18000:
            return None

        weights = bin_counts.astype(float)
        predicted = exp_decay(bin_centers, *popt)
        weighted_mean = np.average(bin_means, weights=weights)
        ss_res = np.sum(weights * (bin_means - predicted)**2)
        ss_tot = np.sum(weights * (bin_means - weighted_mean)**2)
        r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0
        return {
            'lambda_km': float(popt[1]), 
            'r_squared': float(r2), 
            'amplitude': float(popt[0]),
            'offset': float(popt[2]),
            'n_pairs': int(sum(bin_counts))
        }
    except:
        return None


def calculate_seasonal_stats(monthly_data):
    """
    Calculate seasonal anisotropy statistics to detect the 'Seasonal Breathing' signal.
    
    Equinox months: March (3), April (4), September (9), October (10)
    Solstice months: December (12), January (1), June (6), July (7)
    
    Returns:
        dict: Seasonal stats (means, counts, contrast ratio)
    """
    equinox_months = [3, 4, 9, 10]
    solstice_months = [12, 1, 6, 7]
    
    eq_ratios = [d['ratio'] for d in monthly_data if d['month'] in equinox_months]
    sol_ratios = [d['ratio'] for d in monthly_data if d['month'] in solstice_months]
    
    if not eq_ratios or not sol_ratios:
        return None
        
    eq_mean = np.mean(eq_ratios)
    sol_mean = np.mean(sol_ratios)
    
    return {
        'equinox_mean_ratio': eq_mean,
        'solstice_mean_ratio': sol_mean,
        'seasonal_contrast': eq_mean / sol_mean if sol_mean > 0 else 0,
        'n_equinox_months': len(eq_ratios),
        'n_solstice_months': len(sol_ratios),
        'modulation_detected': bool((eq_mean / sol_mean) > 1.2)  # Threshold for significant breathing
    }


def load_station_filter(filter_config, coords_lla):
    """Load station filter and return filtered coordinates + metadata."""
    if filter_config == 'none':
        return coords_lla, {'type': 'all', 'n_stations': len(coords_lla)}
    
    # Load from metadata file (including dynamic_50_metadata.json)
    filter_file = PROCESSED_DIR / filter_config
    if filter_file.exists():
        with open(filter_file) as f:
            meta = json.load(f)
        station_set = set(meta.get('stations', {}).keys())
        filtered = {k: v for k, v in coords_lla.items() if k in station_set}
        
        # Get hemisphere balance
        n_north = sum(1 for s in meta.get('stations', {}).values() 
                      if s.get('classification', {}).get('hemisphere') == 'N')
        n_south = sum(1 for s in meta.get('stations', {}).values() 
                      if s.get('classification', {}).get('hemisphere') == 'S')
        
        return filtered, {
            'type': 'metadata', 
            'file': filter_config,
            'n_stations': len(filtered),
            'n_north': n_north,
            'n_south': n_south,
            'hemisphere_ratio': f"{n_north}:{n_south}"
        }
    
    print_status(f"Filter file not found: {filter_config}, using all stations", "WARNING")
    return coords_lla, {'type': 'all', 'n_stations': len(coords_lla)}


def run_orbital_analysis(filter_name, filter_config, processing_mode, mode_config):
    """Run orbital coupling analysis for one filter + mode combination."""
    
    csv_file = RESULTS_DIR / mode_config['csv_file']
    if not csv_file.exists():
        print_status(f"CSV not found: {csv_file}", "WARNING")
        return None
    
    print_status(f"\n{'='*60}", "INFO")
    print_status(f"Filter: {filter_name} | Mode: {processing_mode}", "TITLE")
    print_status(f"{'='*60}", "INFO")
    
    # Load coordinates
    with open(PROCESSED_DIR / "station_coordinates.json") as f:
        coords_ecef = json.load(f)
    
    coords_lla_all = {k: ecef_to_lla(*v) for k, v in coords_ecef.items()}
    coords_lla, filter_meta = load_station_filter(filter_config, coords_lla_all)
    
    print_status(f"Stations in filter: {filter_meta['n_stations']}", "INFO")
    if 'hemisphere_ratio' in filter_meta:
        print_status(f"Hemisphere balance: {filter_meta['hemisphere_ratio']}", "INFO")
    
    bin_edges = np.logspace(np.log10(MIN_DISTANCE_KM), np.log10(MAX_DISTANCE_KM), N_BINS + 1)
    
    all_metric_results = {}
    
    for metric in METRICS:
        print_status(f"\nProcessing metric: {metric}...", "PROCESS")
        
        # Accumulators: [coh_type][year][month] -> direction bins
        acc = {c: {} for c in COHERENCE_TYPES}
        total_pairs = 0
        
        chunk_size = 2_000_000
        
        for chunk in pd.read_csv(csv_file, chunksize=chunk_size):
            if 'metric' in chunk.columns:
                # Handle prefixed metric names (e.g., ionofree_clock_bias, multi_gnss_clock_bias)
                # Match if metric column ends with the metric name
                mc = chunk[chunk['metric'].str.endswith(metric)]
            else:
                continue
            if mc.empty:
                continue
            
            # Filter by stations if using metadata filter
            if filter_meta['type'] == 'metadata':
                station_mask = mc['station1'].isin(coords_lla.keys()) & mc['station2'].isin(coords_lla.keys())
                mc = mc[station_mask]
                if mc.empty:
                    continue
            
            years = mc['year'].values
            doys = mc['doy'].values
            s1 = mc['station1'].values
            s2 = mc['station2'].values
            dists = mc['distance_km'].values
            
            total_pairs += len(mc)
            
            # Compute months (FIX: proper calendar month calculation)
            months = np.zeros(len(doys), dtype=int)
            doy_bounds = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 366]
            for m in range(1, 13):
                mask = (doys > doy_bounds[m-1]) & (doys <= doy_bounds[m])
                months[mask] = m
            
            # Compute azimuths for unique pairs
            pairs = list(zip(s1, s2))
            unique_pairs = set(pairs)
            
            chunk_az_map = {}
            for sta1, sta2 in unique_pairs:
                if sta1 in coords_lla and sta2 in coords_lla:
                    c1 = coords_lla[sta1]
                    c2 = coords_lla[sta2]
                    _, az = haversine_azimuth(c1[0], c1[1], c2[0], c2[1])
                    chunk_az_map[(sta1, sta2)] = az
            
            azimuths = np.array([chunk_az_map.get((a, b), np.nan) for a, b in zip(s1, s2)])
            
            # Classify directions
            valid_az = ~np.isnan(azimuths)
            ew_mask = np.zeros(len(azimuths), dtype=bool)
            ns_mask = np.zeros(len(azimuths), dtype=bool)
            
            az_v = azimuths[valid_az]
            ew_mask[valid_az] = ((az_v >= 67.5) & (az_v < 112.5)) | ((az_v >= 247.5) & (az_v < 292.5))
            ns_mask[valid_az] = ((az_v < 22.5) | (az_v >= 337.5)) | ((az_v >= 157.5) & (az_v < 202.5))
            
            # Bin data by year-month
            for coh_type in COHERENCE_TYPES:
                col_name = 'coherence' if coh_type == 'msc' else 'phase_alignment'
                if col_name not in mc.columns:
                    continue
                
                vals = mc[col_name].values
                bin_indices = np.digitize(dists, bin_edges) - 1
                valid = ~np.isnan(vals) & (bin_indices >= 0) & (bin_indices < N_BINS)
                
                df_subset = pd.DataFrame({
                    'y': years[valid], 'm': months[valid],
                    'ew': ew_mask[valid], 'ns': ns_mask[valid],
                    'b': bin_indices[valid], 'v': vals[valid]
                })
                
                # Group by Year, Month, Bin for E-W pairs
                for direction, dir_mask in [('ew', df_subset['ew']), ('ns', df_subset['ns'])]:
                    dir_df = df_subset[dir_mask]
                    if dir_df.empty:
                        continue
                    
                    g = dir_df.groupby(['y', 'm', 'b'])['v'].agg(['sum', 'count'])
                    for idx, row in g.iterrows():
                        y, m, b = idx
                        # FIX: Use (year, month) tuple as key to track 36+ unique months
                        ym_key = (int(y), int(m))
                        if ym_key not in acc[coh_type]:
                            acc[coh_type][ym_key] = {
                                'ew': np.zeros((N_BINS, 2)), 
                                'ns': np.zeros((N_BINS, 2)),
                                'year': int(y),
                                'month': int(m)
                            }
                        acc[coh_type][ym_key][direction][b, 0] += row['sum']
                        acc[coh_type][ym_key][direction][b, 1] += row['count']
            
            gc.collect()
        
        print_status(f"  Total pairs processed: {total_pairs:,}", "INFO")
        
        # Fit exponential decay for each month and compute correlation
        metric_results = {}
        
        for coh_type in COHERENCE_TYPES:
            monthly_data = []
            monthly_ratios = []
            monthly_velocities = []
            
            if not acc[coh_type]:
                continue
            
            # Sort by year-month
            sorted_ym = sorted(acc[coh_type].keys())
            
            for ym_key in sorted_ym:
                data = acc[coh_type][ym_key]
                y, m = data['year'], data['month']
                
                ew_counts = data['ew'][:, 1]
                ns_counts = data['ns'][:, 1]
                
                if np.sum(ew_counts) < 100 or np.sum(ns_counts) < 100:
                    continue
                
                ew_means = np.divide(data['ew'][:, 0], ew_counts, out=np.zeros(N_BINS), where=ew_counts>0)
                ns_means = np.divide(data['ns'][:, 0], ns_counts, out=np.zeros(N_BINS), where=ns_counts>0)
                
                valid_ew = ew_counts >= MIN_BIN_COUNT
                valid_ns = ns_counts >= MIN_BIN_COUNT
                
                b_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                
                fit_ew = fit_exponential_binned(b_centers[valid_ew], ew_means[valid_ew], ew_counts[valid_ew])
                fit_ns = fit_exponential_binned(b_centers[valid_ns], ns_means[valid_ns], ns_counts[valid_ns])
                
                if fit_ew and fit_ns and fit_ns['lambda_km'] > 0:
                    ratio = fit_ew['lambda_km'] / fit_ns['lambda_km']
                    mid_doy = 15 + (m-1)*30
                    vel = get_orbital_velocity(mid_doy, y)
                    
                    monthly_data.append({
                        'year': y,
                        'month': m,
                        'ratio': ratio,
                        'velocity': vel,
                        'ew_lambda': fit_ew['lambda_km'],
                        'ns_lambda': fit_ns['lambda_km'],
                        'ew_r2': fit_ew['r_squared'],
                        'ns_r2': fit_ns['r_squared'],
                        'n_ew_pairs': fit_ew['n_pairs'],
                        'n_ns_pairs': fit_ns['n_pairs']
                    })
                    monthly_ratios.append(ratio)
                    monthly_velocities.append(vel)
            
            # Compute correlation
            correlation = {}
            if len(monthly_ratios) >= 6:
                r, p = stats.pearsonr(monthly_velocities, monthly_ratios)
                
                # Monte Carlo significance test
                n_surrogates = 10000
                n_exceeded = 0
                for _ in range(n_surrogates):
                    shuffled = np.random.permutation(monthly_ratios)
                    r_null, _ = stats.pearsonr(monthly_velocities, shuffled)
                    if abs(r_null) >= abs(r):
                        n_exceeded += 1
                mc_p = (n_exceeded + 1) / (n_surrogates + 1)
                
                # Sigma
                sigma = abs(stats.norm.ppf(p/2)) if p > 0 else 0
                
                correlation = {
                    'pearson_r': r,
                    'p_value': p,
                    'monte_carlo_p': mc_p,
                    'sigma': sigma,
                    'n_months': len(monthly_ratios),
                    'n_surrogates': n_surrogates
                }
                
                print_status(f"  {coh_type}: r={r:.3f}, p={p:.2e}, σ={sigma:.1f} ({len(monthly_ratios)} months)", "SUCCESS")
            else:
                print_status(f"  {coh_type}: Insufficient data ({len(monthly_ratios)} months)", "WARNING")
            
            # Calculate seasonal stats
            seasonal_stats = calculate_seasonal_stats(monthly_data)
            if seasonal_stats and seasonal_stats['modulation_detected']:
                print_status(f"  SEASONAL SIGNAL: Equinox/Solstice Ratio = {seasonal_stats['seasonal_contrast']:.2f}", "SUCCESS")

            metric_results[coh_type] = {
                'monthly_anisotropy': {f"{d['year']}-{d['month']:02d}": d for d in monthly_data},
                'correlation': correlation,
                'seasonal_analysis': seasonal_stats,
                'n_months': len(monthly_data)
            }
        
        all_metric_results[metric] = metric_results
    
    return {
        'filter': filter_name,
        'filter_config': filter_config,
        'filter_metadata': filter_meta,
        'processing_mode': processing_mode,
        'mode_description': mode_config['description'],
        'results': all_metric_results,
        'timestamp': datetime.now().isoformat()
    }


def main():
    """Main entry point - comprehensive analysis."""
    parser = argparse.ArgumentParser(description='TEP-GNSS-RINEX Step 2.5: Orbital Coupling (Comprehensive)')
    parser.add_argument('--filter', type=str, default='all',
                        help='Filter: "all", "none", "optimal_100_metadata.json", or "dynamic_50_metadata.json"')
    parser.add_argument('--mode', type=str, default='all',
                        help='Processing mode: "all", "baseline", "ionofree", "multi_gnss"')
    args = parser.parse_args()
    
    print_status("", "INFO")
    print_status("=" * 80, "INFO")
    print_status("STEP 2.5: ORBITAL COUPLING - COMPREHENSIVE ANALYSIS", "TITLE")
    print_status("=" * 80, "INFO")
    print_status("Analyzing: 3 filters × 3 modes × 3 metrics × 2 coherence types", "INFO")
    print_status("=" * 80, "INFO")
    
    # Select filters
    if args.filter == 'all':
        filters = STATION_FILTERS
    else:
        matching = [f for f in STATION_FILTERS if f[1] == args.filter]
        if matching:
            filters = matching
        else:
            print_status(f"Unknown filter: {args.filter}", "ERROR")
            return
    
    # Select modes
    if args.mode == 'all':
        modes = PROCESSING_MODES
    elif args.mode in PROCESSING_MODES:
        modes = {args.mode: PROCESSING_MODES[args.mode]}
    else:
        print_status(f"Unknown mode: {args.mode}", "ERROR")
        return
    
    all_results = []
    summary = {}
    
    for filter_name, filter_config in filters:
        filter_key = filter_name.lower().replace(' ', '_')
        summary[filter_key] = {}
        
        for mode_name, mode_config in modes.items():
            result = run_orbital_analysis(filter_name, filter_config, mode_name, mode_config)
            
            if result:
                all_results.append(result)
                
                # Save individual file for each filter/mode combination
                output_file = RESULTS_DIR / f"step_2_5_orbital_coupling_{filter_key}_{mode_name}.json"
                with open(output_file, 'w') as f:
                    json.dump(result, f, indent=2)
                print_status(f"Saved: {output_file.name}", "SUCCESS")
                
                # Extract best result for summary
                best_r = 0
                best_combo = None
                for metric, metric_res in result['results'].items():
                    for coh_type, coh_res in metric_res.items():
                        corr = coh_res.get('correlation', {})
                        seas = coh_res.get('seasonal_analysis', {})
                        if abs(corr.get('pearson_r', 0)) > abs(best_r):
                            best_r = corr.get('pearson_r', 0)
                            best_combo = {
                                'metric': metric,
                                'coherence_type': coh_type,
                                'pearson_r': best_r,
                                'p_value': corr.get('p_value', 1),
                                'sigma': corr.get('sigma', 0),
                                'n_months': corr.get('n_months', 0),
                                'seasonal_contrast': seas.get('seasonal_contrast', 0) if seas else 0
                            }
                
                summary[filter_key][mode_name] = best_combo
    
    # Save comprehensive summary
    summary_output = {
        'step': '2.5',
        'name': 'orbital_coupling_comprehensive',
        'timestamp': datetime.now().isoformat(),
        'methodology': 'Monthly EW/NS lambda ratio vs orbital velocity correlation',
        'n_combinations': len(all_results),
        'filters_tested': [f[0] for f in filters],
        'modes_tested': list(modes.keys()),
        'summary': summary,
        'code_reference': {
            'code_r': -0.888,
            'code_sigma': 5.1,
            'code_n_months': 300
        }
    }
    
    # Find overall best result
    best_overall = None
    for filter_key, modes_dict in summary.items():
        for mode, result in modes_dict.items():
            if result and (not best_overall or abs(result.get('pearson_r', 0)) > abs(best_overall.get('pearson_r', 0))):
                best_overall = {**result, 'filter': filter_key, 'mode': mode}
    
    summary_output['best_result'] = best_overall
    
    summary_file = RESULTS_DIR / "step_2_5_orbital_coupling_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary_output, f, indent=2)
    
    # Print summary table
    print_status("\n" + "=" * 80, "INFO")
    print_status("SUMMARY: Best correlation per filter/mode", "TITLE")
    print_status("=" * 80, "INFO")
    print_status(f"{'Filter':<20} {'Mode':<12} {'Metric':<12} {'r':>8} {'σ':>6} {'SeasContrast':>12}", "INFO")
    print_status("-" * 80, "INFO")
    
    for filter_key, modes_dict in summary.items():
        for mode, result in modes_dict.items():
            if result:
                seas = result.get('seasonal_contrast', 0)
                seas_str = f"{seas:.2f}x" if seas > 0 else "N/A"
                print_status(
                    f"{filter_key:<20} {mode:<12} {result['metric']:<12} "
                    f"{result['pearson_r']:>8.3f} {result['sigma']:>6.1f} {seas_str:>12}",
                    "INFO"
                )
    
    print_status("-" * 80, "INFO")
    print_status(f"CODE Reference: r=-0.888, σ=5.1, 300 months", "INFO")
    
    if best_overall:
        print_status("", "INFO")
        print_status(f"BEST: {best_overall['filter']}/{best_overall['mode']} "
                    f"r={best_overall['pearson_r']:.3f} σ={best_overall['sigma']:.1f}", "SUCCESS")
    
    print_status(f"\nSaved summary: {summary_file.name}", "SUCCESS")
    print_status("=" * 80, "INFO")
    print_status("STEP 2.5 COMPLETE", "SUCCESS")
    print_status("=" * 80, "INFO")


if __name__ == '__main__':
    main()
