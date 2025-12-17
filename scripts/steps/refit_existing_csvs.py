#!/usr/bin/env python3
"""Re-run fitting on existing CSV files with fixed prebinned function."""

import sys
import json
import numpy as np
from pathlib import Path
from collections import defaultdict
from scipy.optimize import curve_fit

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.logger import print_status

RESULTS_DIR = Path(__file__).parent.parent.parent / "results" / "outputs"

def exp_decay(r, A, lam, C0):
    return A * np.exp(-r/lam) + C0

def fit_exponential_prebinned(bin_centers, bin_means, bin_counts=None):
    """Fit exponential decay to PRE-BINNED data (no re-binning)."""
    x = np.array(bin_centers)
    y = np.array(bin_means)
    
    if len(x) < 5:
        return None
    
    if bin_counts is not None:
        w = np.sqrt(np.array(bin_counts))
        w[w == 0] = 1
    else:
        w = np.ones_like(x)
    
    try:
        popt, pcov = curve_fit(
            exp_decay, x, y,
            p0=[0.5, 2000, 0.5],
            sigma=1.0/w,
            bounds=([0, 100, -1], [2, 20000, 1]),
            maxfev=5000
        )
        
        predicted = exp_decay(x, *popt)
        ss_res = np.sum((y - predicted)**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        return {
            'amplitude': popt[0],
            'correlation_length_km': popt[1],
            'offset': popt[2],
            'r_squared': r2,
            'success': True
        }
    except Exception as e:
        return None

def refit_csv(filter_suffix):
    """Re-run fitting on existing CSVs for a given filter."""
    import csv as csv_module
    
    modes_config = {
        'baseline': ['clock_bias', 'pos_jitter', 'clock_drift'],
        'ionofree': ['ionofree_clock_bias', 'ionofree_pos_jitter', 'ionofree_clock_drift'],
        'multi_gnss': ['mgex_clock_bias', 'mgex_pos_jitter', 'mgex_clock_drift'],
        'precise': ['precise_clock_bias', 'precise_pos_jitter', 'precise_clock_drift']
    }
    
    metric_to_mode = {}
    for mode, metrics in modes_config.items():
        for m in metrics:
            metric_to_mode[m] = mode
    
    print_status(f"\n{'='*80}", "INFO")
    print_status(f"Re-fitting CSVs for: {filter_suffix}", "INFO")
    print_status(f"{'='*80}\n", "INFO")
    
    # Stream and bin - USE LOG BINS LIKE CODE LONGSPAN
    N_BINS = 40  # Same as CODE longspan
    MIN_DIST = 50  # Same as CODE longspan  
    MAX_DIST = 13000  # Same as CODE longspan
    
    # Create log-spaced bin edges (same as CODE longspan)
    bin_edges = np.logspace(np.log10(MIN_DIST), np.log10(MAX_DIST), N_BINS + 1)
    
    metrics_bins = defaultdict(lambda: {
        'coh_sum': np.zeros(N_BINS),
        'coh_count': np.zeros(N_BINS, dtype=np.int64),
        'pa_sum': np.zeros(N_BINS),
        'pa_count': np.zeros(N_BINS, dtype=np.int64)
    })
    
    total_rows = 0
    skipped_out_of_range = 0
    
    for mode_name in modes_config.keys():
        csv_path = RESULTS_DIR / f"step_2_0_pairs_{mode_name}_{filter_suffix}.csv"
        if not csv_path.exists():
            print_status(f"  Skipping {csv_path.name} (not found)", "WARNING")
            continue
        
        print_status(f"  Streaming {csv_path.name}...", "INFO")
        with open(csv_path, 'r') as f:
            reader = csv_module.DictReader(f)
            for row in reader:
                total_rows += 1
                metric_name = row['metric']
                dist = float(row['distance_km'])
                coh = float(row['coherence'])
                
                # Skip pairs outside range (like CODE longspan)
                if dist < MIN_DIST or dist > MAX_DIST:
                    skipped_out_of_range += 1
                    continue
                
                # Find log bin index
                bin_idx = np.searchsorted(bin_edges, dist, side='right') - 1
                bin_idx = max(0, min(bin_idx, N_BINS - 1))
                
                if np.isfinite(coh):
                    metrics_bins[metric_name]['coh_sum'][bin_idx] += coh
                    metrics_bins[metric_name]['coh_count'][bin_idx] += 1
                
                if 'phase_alignment' in row and row['phase_alignment']:
                    try:
                        pa = float(row['phase_alignment'])
                        if np.isfinite(pa):
                            metrics_bins[metric_name]['pa_sum'][bin_idx] += pa
                            metrics_bins[metric_name]['pa_count'][bin_idx] += 1
                    except:
                        pass
                
                if total_rows % 10000000 == 0:
                    print_status(f"    Processed {total_rows/1e6:.0f}M rows...", "INFO")
    
    print_status(f"  Streamed {total_rows:,} rows ({skipped_out_of_range:,} out of range) into {N_BINS} log bins", "SUCCESS")
    
    # Convert bins to arrays - use LOG bin centers
    # Bin centers in log space: geometric mean of edges
    bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])
    
    metrics_data = {}
    for metric_name, bins in metrics_bins.items():
        
        valid_coh = bins['coh_count'] > 0
        distances = bin_centers[valid_coh]
        coherences = bins['coh_sum'][valid_coh] / bins['coh_count'][valid_coh]
        
        valid_pa = bins['pa_count'] > 0
        pa_distances = bin_centers[valid_pa]
        phase_alignments = bins['pa_sum'][valid_pa] / bins['pa_count'][valid_pa]
        
        metrics_data[metric_name] = {
            'distances': distances,
            'coherences': coherences,
            'coh_counts': bins['coh_count'][valid_coh],
            'n_pairs': int(bins['coh_count'].sum()),
            'pa_distances': pa_distances,
            'phase_alignments': phase_alignments,
            'pa_counts': bins['pa_count'][valid_pa],
            'n_pa_pairs': int(bins['pa_count'].sum())
        }
    
    # Fit and report
    all_results = defaultdict(dict)
    
    print_status("\nFitting Results:", "INFO")
    print_status("-" * 80, "INFO")
    
    for metric_name, data in metrics_data.items():
        dists = data['distances']
        cohs = data['coherences']
        coh_counts = data['coh_counts']
        n_pairs = data['n_pairs']
        
        if len(dists) < 10:
            continue
        
        mode_name = metric_to_mode.get(metric_name, 'baseline')
        
        # Fit MSC
        res = fit_exponential_prebinned(dists, cohs, coh_counts)
        if res:
            all_results[mode_name][metric_name] = res
            res['n_pairs'] = n_pairs
            res['n_bins'] = len(dists)
            res['metric_type'] = 'normalized_msc'
            lam = res['correlation_length_km']
            r2 = res['r_squared']
            tep = "YES" if 500 < lam < 5000 and r2 > 0.5 else "NO"
            print_status(f"  [{mode_name.upper()}] {metric_name} (MSC): λ={lam:.0f}km, R²={r2:.3f}, N={n_pairs:,}, TEP={tep}", 
                       "SUCCESS" if tep == "YES" else "INFO")
        
        # Fit Phase
        pa_dists = data['pa_distances']
        phase_aligns = data['phase_alignments']
        pa_counts = data['pa_counts']
        n_pa_pairs = data['n_pa_pairs']
        
        if len(pa_dists) >= 10:
            res_pa = fit_exponential_prebinned(pa_dists, phase_aligns, pa_counts)
            if res_pa:
                pa_key = f"{metric_name}_phase_alignment"
                all_results[mode_name][pa_key] = res_pa
                res_pa['n_pairs'] = n_pa_pairs
                res_pa['n_bins'] = len(pa_dists)
                res_pa['metric_type'] = 'phase_alignment'
                lam = res_pa['correlation_length_km']
                r2 = res_pa['r_squared']
                tep = "YES" if 500 < lam < 5000 and r2 > 0.5 else "NO"
                print_status(f"  [{mode_name.upper()}] {metric_name} (Phase): λ={lam:.0f}km, R²={r2:.3f}, N={n_pa_pairs:,}, TEP={tep}", 
                           "SUCCESS" if tep == "YES" else "INFO")
    
    # Save results
    output = {
        'filter_mode': filter_suffix,
        'total_rows': total_rows,
        'analysis_by_mode': dict(all_results)
    }
    
    output_file = RESULTS_DIR / f"step_2_0_refit_{filter_suffix}.json"
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2, default=float)
    
    print_status(f"\nResults saved: {output_file}", "SUCCESS")
    return all_results

if __name__ == "__main__":
    # Re-fit all available CSVs
    for suffix in ['all_stations', 'optimal_100']:
        refit_csv(suffix)
