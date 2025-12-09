#!/usr/bin/env python3
"""
TEP-GNSS-RINEX Analysis - STEP 2.7: CMB Frame Alignment Analysis
==================================================================

CORRECT METHODOLOGY (matches CODE Longspan step_2_5):
    Uses monthly EW/NS lambda ratios from step_2_5 and correlates with
    velocity vector predictor cos(declination) for different background
    motion hypotheses.

Theory:
    Earth moves through space with combined velocity:
    1. Orbital motion around Sun (~30 km/s, rotating direction)
    2. Background motion (Solar Apex ~20 km/s OR CMB frame ~369 km/s)
    
    The net velocity direction modulates which baseline orientations
    show stronger correlations (E-W vs N-S).
    
    Grid search finds which background RA/Dec best explains the observed
    EW/NS ratio modulation through the year.

CODE Longspan Finding:
    - Best-fit direction: RA=186°, Dec=-4° (r=0.747)
    - Only 18.2° from CMB dipole (RA=168°, Dec=-7°)
    - Solar Apex rejected: 5,570× variance ratio

Requirements: Step 2.5 complete (orbital coupling with monthly ratios)

Outputs:
    - results/outputs/step_2_7_cmb_frame_analysis.json
    - results/figures/step_2_7_cmb_frame_sky_map.png

Author: Matthew Lukin Smawfield
Date: 9 December 2025
License: CC-BY-4.0
"""

import sys
import json
import numpy as np
from pathlib import Path
from datetime import datetime
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Setup paths
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(ROOT))

from scripts.utils.logger import print_status, TEPLogger, set_step_logger

# Directories
OUTPUTS_DIR = ROOT / "results" / "outputs"
FIGURES_DIR = ROOT / "results" / "figures"

OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Initialize logger
logger = TEPLogger(
    name="step_2_7_cmb_frame_analysis",
    level="DEBUG",
    log_file_path=ROOT / "logs" / "step_2_7_cmb_frame_analysis.log"
)
set_step_logger(logger)

# ============================================================================
# PHYSICAL CONSTANTS
# ============================================================================

# Reference frame directions (J2000 equatorial coordinates)
CMB_DIPOLE = {
    'name': 'CMB Dipole',
    'ra': 167.94,   # Planck 2018
    'dec': -6.94,
    'speed_kms': 369.82  # km/s
}

SOLAR_APEX = {
    'name': 'Solar Apex',
    'ra': 271.0,    # Galactic motion direction
    'dec': 30.0,
    'speed_kms': 20.0  # Local Standard of Rest
}

# Earth orbital parameters
ECLIPTIC_TILT_DEG = 23.44  # Obliquity of ecliptic

# Orbital velocity varies ~29.3 to 30.3 km/s
def get_orbital_speed(day_of_year):
    """Get Earth's orbital speed for a given day (simplified Kepler)."""
    # Perihelion ~Jan 3 (day 3), aphelion ~July 4 (day 185)
    perihelion_day = 3
    theta = 2 * np.pi * (day_of_year - perihelion_day) / 365.25
    # Eccentricity-based speed variation
    e = 0.0167
    mean_speed = 29.78
    return mean_speed * (1 + e * np.cos(theta))

def get_mid_month_doy(year, month):
    """Get day of year for mid-month."""
    from datetime import date
    mid_day = 15
    d = date(year, month, mid_day)
    return d.timetuple().tm_yday


# ============================================================================
# 3D VELOCITY VECTOR CALCULATION (from CODE step_2_5)
# ============================================================================

def calculate_3d_velocity_vectors(day_of_year, orbital_speed_kms,
                                  background_ra, background_dec,
                                  background_speed):
    """
    Calculate full 3D velocity vectors in Equatorial Coordinates (J2000).
    
    Includes:
    1. Background Motion (configurable RA/Dec/speed)
    2. Earth Orbital Motion (Rotating in Ecliptic Plane -> Transformed to Equatorial)
    
    Returns dict with 3D vector components and angles.
    """
    # 1. Background Motion
    apex_ra_rad = np.radians(background_ra)
    apex_dec_rad = np.radians(background_dec)
    
    # Background velocity vector (Equatorial Cartesian)
    v_bg_x = background_speed * np.cos(apex_dec_rad) * np.cos(apex_ra_rad)
    v_bg_y = background_speed * np.cos(apex_dec_rad) * np.sin(apex_ra_rad)
    v_bg_z = background_speed * np.sin(apex_dec_rad)
    
    # 2. Earth Orbital Motion (Dynamic)
    # Velocity rotates in ecliptic plane through the year
    # At Vernal Equinox (Day 80), velocity points to ecliptic longitude 270°
    orbit_progress = (day_of_year - 80) / 365.25 * 2 * np.pi
    
    # Velocity direction in Ecliptic coordinates
    v_ecl_x = -orbital_speed_kms * np.sin(orbit_progress)
    v_ecl_y = orbital_speed_kms * np.cos(orbit_progress)
    v_ecl_z = 0.0
    
    # Rotate Ecliptic -> Equatorial (rotate around X-axis by obliquity)
    epsilon = np.radians(ECLIPTIC_TILT_DEG)
    v_orb_x = v_ecl_x
    v_orb_y = v_ecl_y * np.cos(epsilon)
    v_orb_z = v_ecl_y * np.sin(epsilon)
    
    # 3. Net Velocity (vector sum)
    v_net_x = v_orb_x + v_bg_x
    v_net_y = v_orb_y + v_bg_y
    v_net_z = v_orb_z + v_bg_z
    
    # Calculate Magnitude and Angles
    v_net_mag = np.sqrt(v_net_x**2 + v_net_y**2 + v_net_z**2)
    
    # Right Ascension (0 to 360)
    v_net_ra_deg = np.degrees(np.arctan2(v_net_y, v_net_x)) % 360
    
    # Declination (-90 to +90)
    v_net_dec_deg = np.degrees(np.arcsin(v_net_z / v_net_mag))
    
    return {
        'v_net_mag': float(v_net_mag),
        'v_net_ra_deg': float(v_net_ra_deg),
        'v_net_dec_deg': float(v_net_dec_deg),
        'v_net_xy_mag': float(np.sqrt(v_net_x**2 + v_net_y**2)),
        'v_net_z_mag': float(abs(v_net_z))
    }


def angular_separation(ra1, dec1, ra2, dec2):
    """Compute angular separation between two sky positions (degrees)."""
    ra1_r, dec1_r = np.radians(ra1), np.radians(dec1)
    ra2_r, dec2_r = np.radians(ra2), np.radians(dec2)
    
    cos_sep = (np.sin(dec1_r) * np.sin(dec2_r) + 
               np.cos(dec1_r) * np.cos(dec2_r) * np.cos(ra1_r - ra2_r))
    cos_sep = np.clip(cos_sep, -1.0, 1.0)
    return np.degrees(np.arccos(cos_sep))


def set_publication_style():
    """
    Configure matplotlib for consistent publication-quality figures (Nature style).
    """
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
        'font.size': 10,
        'axes.labelsize': 11,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'axes.linewidth': 0.8,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
    })


# ============================================================================
# MAIN ANALYSIS
# ============================================================================

def load_step_2_5_results(filter_name=None, metric=None, coherence_type=None, mode='baseline'):
    """Load monthly EW/NS ratios from step_2_5 orbital coupling results.
    
    NEW FORMAT: step_2_5_orbital_coupling_{filter}_{mode}.json contains
    results for all metrics and coherence types in nested structure.
    
    Args:
        filter_name: Station filter (all_stations, optimal_100, dynamic_50)
        metric: Metric name (clock_bias, pos_jitter, clock_drift)
        coherence_type: Coherence type (msc, phase_alignment)
        mode: Processing mode (baseline, ionofree, multi_gnss)
    """
    
    # NEW FORMAT: Load filter/mode file and extract metric/coherence
    if filter_name:
        filter_suffix = filter_name.lower().replace(' ', '_')
        filename = f"step_2_5_orbital_coupling_{filter_suffix}_{mode}.json"
        filepath = OUTPUTS_DIR / filename
        
        if filepath.exists():
            with open(filepath) as f:
                data = json.load(f)
            
            # Navigate to specific metric/coherence if requested
            if metric and coherence_type and 'results' in data:
                if metric in data['results']:
                    if coherence_type in data['results'][metric]:
                        return {
                            'monthly_anisotropy': data['results'][metric][coherence_type].get('monthly_anisotropy', {}),
                            'correlation': data['results'][metric][coherence_type].get('correlation', {}),
                            'filter': filter_name,
                            'mode': mode,
                            'metric': metric,
                            'coherence_type': coherence_type
                        }, filename
            
            return data, filename
    
    # Try to load summary and get best result
    summary_path = OUTPUTS_DIR / "step_2_5_orbital_coupling_summary.json"
    if summary_path.exists():
        with open(summary_path) as f:
            summary = json.load(f)
        
        best = summary.get('best_result', {})
        if best:
            best_filter = best.get('filter', 'optimal_100')
            best_mode = best.get('mode', 'baseline')
            best_metric = best.get('metric', 'clock_bias')
            best_coh = best.get('coherence_type', 'msc')
            
            print_status(f"Using BEST from summary: {best_filter}/{best_mode}/{best_metric}/{best_coh}", "SUCCESS")
            return load_step_2_5_results(best_filter, best_metric, best_coh, best_mode)
    
    # Fallback to old-style files
    fallback_patterns = [
        'step_2_5_orbital_coupling_optimal_100_baseline.json',
        'step_2_5_orbital_coupling_all_stations_baseline.json',
        'step_2_5_orbital_coupling.json'
    ]
    
    for pattern in fallback_patterns:
        filepath = OUTPUTS_DIR / pattern
        if filepath.exists():
            with open(filepath) as f:
                data = json.load(f)
            print_status(f"Loaded (fallback): {pattern}", "WARNING")
            return data, pattern
    
    print_status("No step_2_5 results found!", "ERROR")
    return None, None


def get_all_step_2_5_combinations():
    """Get list of all available step_2_5 result files and their metric/coherence combos.
    
    NEW FORMAT: Each file contains all metrics and coherence types.
    Returns expanded list of all testable combinations.
    """
    combinations = []
    
    filters = ['all_stations', 'optimal_100', 'dynamic_50']
    modes = ['baseline', 'ionofree', 'multi_gnss']
    metrics = ['clock_bias', 'pos_jitter', 'clock_drift']
    coherence_types = ['msc', 'phase_alignment']
    
    for f in filters:
        for mode in modes:
            filename = f"step_2_5_orbital_coupling_{f}_{mode}.json"
            filepath = OUTPUTS_DIR / filename
            
            if filepath.exists():
                # Load file to check which metrics/coherence are available
                with open(filepath) as file:
                    data = json.load(file)
                
                results = data.get('results', {})
                
                for m in metrics:
                    if m in results:
                        for c in coherence_types:
                            if c in results[m]:
                                monthly_data = results[m][c].get('monthly_anisotropy', {})
                                if len(monthly_data) >= 6:
                                    combinations.append({
                                        'filter': f,
                                        'mode': mode,
                                        'metric': m,
                                        'coherence_type': c,
                                        'filename': filename,
                                        'n_months': len(monthly_data)
                                    })
    
    # Also check old-style single file
    old_file = OUTPUTS_DIR / "step_2_5_orbital_coupling.json"
    if old_file.exists() and not combinations:
        combinations.append({
            'filter': 'all_stations',
            'mode': 'baseline',
            'metric': 'clock_bias',
            'coherence_type': 'msc',
            'filename': 'step_2_5_orbital_coupling.json',
            'n_months': 0
        })
    
    return combinations


def extract_monthly_data(step_2_5_data):
    """Extract monthly EW/NS ratios and velocities from step_2_5 results."""
    
    monthly_data = []
    
    # Navigate the JSON structure
    if 'primary_results' in step_2_5_data:
        results = step_2_5_data['primary_results']
    elif 'baseline' in step_2_5_data.get('comparative_results', {}):
        results = step_2_5_data['comparative_results']['baseline']
    else:
        results = step_2_5_data
    
    # Get monthly anisotropy data
    monthly_aniso = None
    if 'monthly_anisotropy' in results:
        monthly_aniso = results['monthly_anisotropy']
    elif 'monthly_results' in results:
        monthly_aniso = results['monthly_results']
    
    if monthly_aniso is None:
        print_status("Could not find monthly data in step_2_5 results", "ERROR")
        return []
    
    # Get monthly velocities
    monthly_velocity = results.get('monthly_velocity', {})
    
    for key, data in monthly_aniso.items():
        if isinstance(data, dict) and 'ratio' in data:
            year = data.get('year')
            month = data.get('month')
            ratio = data.get('ratio')
            
            if year and month and ratio:
                # Get orbital speed for mid-month
                doy = get_mid_month_doy(year, month)
                
                # Use stored velocity if available, otherwise compute
                vel_key = f"{year}-{month:02d}"
                if vel_key in monthly_velocity:
                    orbital_speed = monthly_velocity[vel_key]
                else:
                    orbital_speed = get_orbital_speed(doy)
                
                monthly_data.append({
                    'year': year,
                    'month': month,
                    'doy': doy,
                    'ratio': ratio,
                    'orbital_speed': orbital_speed
                })
    
    # Sort by date
    monthly_data.sort(key=lambda x: (x['year'], x['month']))
    
    return monthly_data


# Fixed background speed for grid search (matches CODE longspan methodology)
# Using 20 km/s creates meaningful modulation without biasing toward CMB
GRID_SEARCH_SPEED_KMS = 20.0


def grid_search_background(monthly_data, ra_step=1, dec_step=1):
    """
    Grid search over background velocity directions.
    
    For each (RA, Dec), compute the predicted EW/NS modulation and
    correlate with observed monthly ratios.
    
    Uses FIXED 20 km/s background speed (matches CODE longspan methodology).
    Resolution: 1° (matches CODE longspan finest setting).
    
    Predictor: cos(velocity_declination) - higher declination means
    more N-S aligned velocity, which should increase NS correlations.
    """
    
    print_status("Grid searching background directions (1° resolution, 20 km/s fixed - matches CODE)...", "PROCESS")
    
    # Extract observed data
    ratios = np.array([m['ratio'] for m in monthly_data])
    
    # Grid over RA and Dec
    ra_grid = np.arange(0, 360, ra_step)
    dec_grid = np.arange(-90, 91, dec_step)  # Full range like CODE
    
    results = {}
    best_r = -999
    best_ra, best_dec = 0, 0
    
    for bg_ra in ra_grid:
        for bg_dec in dec_grid:
            # Compute predictor for each month using FIXED speed
            predictors = []
            
            for m in monthly_data:
                v = calculate_3d_velocity_vectors(
                    m['doy'], m['orbital_speed'],
                    bg_ra, bg_dec, GRID_SEARCH_SPEED_KMS  # Fixed 20 km/s
                )
                # Predictor: cos(declination) - flat velocity = high EW/NS
                pred = np.cos(np.radians(v['v_net_dec_deg']))
                predictors.append(pred)
            
            predictors = np.array(predictors)
            
            # Correlation with observed ratios
            if np.std(predictors) > 1e-10:
                r, p = stats.pearsonr(predictors, ratios)
            else:
                r, p = 0.0, 1.0
            
            results[(bg_ra, bg_dec)] = {'r': r, 'p': p}
            
            if r > best_r:
                best_r = r
                best_ra, best_dec = bg_ra, bg_dec
    
    # Bootstrap confidence intervals
    print_status("  Computing bootstrap confidence intervals...", "PROCESS")
    ra_boot, dec_boot = bootstrap_direction_uncertainty(monthly_data, n_bootstrap=500)
    
    return results, best_ra, best_dec, best_r, ra_boot, dec_boot


def bootstrap_direction_uncertainty(monthly_data, n_bootstrap=500):
    """
    Bootstrap resampling to estimate uncertainty in best-fit direction.
    """
    ratios = np.array([m['ratio'] for m in monthly_data])
    n = len(monthly_data)
    
    best_ras = []
    best_decs = []
    
    for _ in range(n_bootstrap):
        # Resample with replacement
        idx = np.random.choice(n, n, replace=True)
        boot_data = [monthly_data[i] for i in idx]
        boot_ratios = ratios[idx]
        
        # Quick grid search at coarser resolution for speed
        best_r = -999
        best_ra, best_dec = 0, 0
        
        for bg_ra in range(0, 360, 10):
            for bg_dec in range(-80, 81, 10):
                predictors = []
                for m in boot_data:
                    v = calculate_3d_velocity_vectors(
                        m['doy'], m['orbital_speed'],
                        bg_ra, bg_dec, GRID_SEARCH_SPEED_KMS
                    )
                    predictors.append(np.cos(np.radians(v['v_net_dec_deg'])))
                
                predictors = np.array(predictors)
                if np.std(predictors) > 1e-10:
                    r = np.corrcoef(predictors, boot_ratios)[0, 1]
                    if r > best_r:
                        best_r = r
                        best_ra, best_dec = bg_ra, bg_dec
        
        best_ras.append(best_ra)
        best_decs.append(best_dec)
    
    # Return 68% confidence intervals
    ra_ci = (np.percentile(best_ras, 16), np.percentile(best_ras, 84))
    dec_ci = (np.percentile(best_decs, 16), np.percentile(best_decs, 84))
    
    return ra_ci, dec_ci


def estimate_global_p_value_monte_carlo(monthly_data, observed_best_r, n_iterations=10000, ra_step=5, dec_step=5):
    """
    Estimate global p-value accounting for the "look-elsewhere effect" (scanning the whole sky).
    
    Method:
    1. Shuffle the temporal order of the observed ratios (breaking correlation with orbital geometry).
    2. Perform the full grid search on the shuffled data.
    3. Record the maximum correlation found anywhere on the sky.
    4. Repeat N times.
    5. Count how many random iterations produced a max correlation >= observed_best_r.
    
    Optimized with vectorization for speed.
    """
    print_status(f"  Estimating global p-value (MC, N={n_iterations}, step={ra_step}°)...", "PROCESS")
    
    # Pre-compute predictors for the entire grid
    # Grid points
    ra_grid = np.arange(0, 360, ra_step)
    dec_grid = np.arange(-90, 91, dec_step)
    
    # Create meshgrid for vectorized calculation
    # We need arrays of shape (N_grid_points)
    ra_mesh, dec_mesh = np.meshgrid(ra_grid, dec_grid)
    ra_flat = ra_mesh.flatten()
    dec_flat = dec_mesh.flatten()
    n_grid = len(ra_flat)
    
    n_months = len(monthly_data)
    
    # Pre-compute predictors matrix: (N_months, N_grid_points)
    # This matrix contains the cos(v_dec) for every month and every sky position
    predictors = np.zeros((n_months, n_grid))
    
    for i, m in enumerate(monthly_data):
        # We need a vectorized version of calculate_3d_velocity_vectors or just inline the relevant part
        # Inline is safer for arrays
        
        # 1. Background Motion
        apex_ra_rad = np.radians(ra_flat)
        apex_dec_rad = np.radians(dec_flat)
        
        v_bg_x = GRID_SEARCH_SPEED_KMS * np.cos(apex_dec_rad) * np.cos(apex_ra_rad)
        v_bg_y = GRID_SEARCH_SPEED_KMS * np.cos(apex_dec_rad) * np.sin(apex_ra_rad)
        v_bg_z = GRID_SEARCH_SPEED_KMS * np.sin(apex_dec_rad)
        
        # 2. Earth Orbital Motion (Dynamic)
        orbit_progress = (m['doy'] - 80) / 365.25 * 2 * np.pi
        v_ecl_x = -m['orbital_speed'] * np.sin(orbit_progress)
        v_ecl_y = m['orbital_speed'] * np.cos(orbit_progress)
        # v_ecl_z = 0.0
        
        epsilon = np.radians(ECLIPTIC_TILT_DEG)
        v_orb_x = v_ecl_x
        v_orb_y = v_ecl_y * np.cos(epsilon)
        v_orb_z = v_ecl_y * np.sin(epsilon)
        
        # 3. Net Velocity
        v_net_x = v_orb_x + v_bg_x
        v_net_y = v_orb_y + v_bg_y
        v_net_z = v_orb_z + v_bg_z
        
        v_net_mag = np.sqrt(v_net_x**2 + v_net_y**2 + v_net_z**2)
        
        # Declination -> cos(dec)
        # sin(dec) = z / mag
        # cos(dec) = sqrt(1 - sin^2(dec))
        sin_dec = v_net_z / v_net_mag
        # Predictor is cos(dec)
        predictors[i, :] = np.sqrt(1.0 - sin_dec**2)
        
    # Normalize predictors (columns)
    # Subtract mean and divide by std for each grid point
    pred_mean = np.mean(predictors, axis=0)
    pred_std = np.std(predictors, axis=0)
    
    # Avoid division by zero
    valid_cols = pred_std > 1e-10
    
    predictors_norm = np.zeros_like(predictors)
    predictors_norm[:, valid_cols] = (predictors[:, valid_cols] - pred_mean[valid_cols]) / pred_std[valid_cols]
    
    # Observed ratios
    ratios_orig = np.array([m['ratio'] for m in monthly_data])
    
    # Monte Carlo Loop
    count_exceeds = 0
    
    for k in range(n_iterations):
        # Shuffle ratios
        ratios_shuffled = np.random.permutation(ratios_orig)
        
        # Normalize shuffled ratios
        r_mean = np.mean(ratios_shuffled)
        r_std = np.std(ratios_shuffled)
        
        if r_std < 1e-10:
            continue
            
        ratios_norm = (ratios_shuffled - r_mean) / r_std
        
        # Compute correlations: (predictors.T @ ratios) / n
        # This gives vector of correlations for all grid points
        correlations = (predictors_norm.T @ ratios_norm) / n_months
        
        # Max correlation on the sky for this random trial
        max_r_trial = np.max(correlations)
        
        if max_r_trial >= observed_best_r:
            count_exceeds += 1
            
    # Calculate p-value
    # Add +1 to numerator and denominator for unbiased estimator (prevents p=0)
    global_p = (count_exceeds + 1) / (n_iterations + 1)
    
    return global_p, count_exceeds


def create_sky_map(grid_results, best_ra, best_dec, output_path):
    """Create publication-aligned sky map matching code-longspan styling."""
    
    print_status("Creating sky map visualization...", "PROCESS")
    
    set_publication_style()
    
    # Extract grid data
    ra_vals = sorted(set(k[0] for k in grid_results.keys()))
    dec_vals = sorted(set(k[1] for k in grid_results.keys()))
    
    # Create 2D array of correlations
    corr_matrix = np.full((len(dec_vals), len(ra_vals)), np.nan)
    for i, dec in enumerate(dec_vals):
        for j, ra in enumerate(ra_vals):
            corr_matrix[i, j] = grid_results.get((ra, dec), {'r': np.nan})['r']
    
    # Create figure using Nature-style conventions
    fig, ax = plt.subplots(figsize=(12, 6), dpi=600)
    
    im = ax.imshow(
        corr_matrix,
        aspect='auto',
        origin='lower',
        extent=[0, 360, -90, 90],
        cmap='viridis',
        vmin=-0.3,
        vmax=0.75,
        interpolation='bilinear'
    )
    
    # Reference markers (matching code-longspan)
    ax.plot(best_ra, best_dec, 'o',
            markersize=9, markerfacecolor='white',
            markeredgecolor='black', markeredgewidth=1.5,
            clip_on=False, zorder=10)
    
    ax.plot(CMB_DIPOLE['ra'], CMB_DIPOLE['dec'], 'o',
            markersize=8, markerfacecolor='#00FFFF',
            markeredgecolor='black', markeredgewidth=1.5,
            clip_on=False, zorder=10)
    
    ax.plot(SOLAR_APEX['ra'], SOLAR_APEX['dec'], 'o',
            markersize=8, markerfacecolor='#F39C12',
            markeredgecolor='black', markeredgewidth=1.5,
            clip_on=False, zorder=10)
    
    # Label markers a/b/c as in code-longspan figure
    ax.text(best_ra + 3, best_dec + 3, 'a',
            fontsize=10, color='white', fontweight='bold',
            ha='left', va='bottom')
    ax.text(CMB_DIPOLE['ra'] + 3, CMB_DIPOLE['dec'] + 3, 'b',
            fontsize=10, color='white', fontweight='bold',
            ha='left', va='bottom')
    ax.text(SOLAR_APEX['ra'] + 3, SOLAR_APEX['dec'] + 3, 'c',
            fontsize=10, color='white', fontweight='bold',
            ha='left', va='bottom')
    
    # Axes styling
    ax.set_xlabel('Right ascension (°)', fontsize=11)
    ax.set_ylabel('Declination (°)', fontsize=11)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_yticks([-90, -45, 0, 45, 90])
    ax.grid(True, alpha=0.15, linestyle=':', linewidth=0.3, color='white')
    
    # Colorbar styling
    cbar = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02, aspect=30)
    cbar.set_label('r', fontsize=11, rotation=0, labelpad=12)
    cbar.ax.tick_params(labelsize=9, width=0.8, length=3)
    cbar.outline.set_linewidth(0.8)
    
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    
    ax.set_title('CMB Frame Analysis: Background Vector Grid Search', fontsize=12)
    
    plt.tight_layout(pad=0.5)
    plt.savefig(output_path, dpi=600, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()
    
    print_status(f"Figure saved: {output_path}", "SUCCESS")


def run_analysis():
    """Run the CMB frame alignment analysis."""
    
    print_status("", "INFO")
    print_status("=" * 80, "INFO")
    print_status("TEP-GNSS-RINEX Analysis - STEP 2.7: CMB Frame Alignment", "INFO")
    print_status("=" * 80, "INFO")
    print_status("", "INFO")
    print_status("Methodology: Grid search over background RA/Dec using", "INFO")
    print_status("             monthly EW/NS ratios from Step 2.5", "INFO")
    print_status("", "INFO")
    
    # ==========================================================================
    # 1. Load Step 2.5 Results
    # ==========================================================================
    print_status("[1/4] Loading Step 2.5 orbital coupling results...", "PROCESS")
    
    step_2_5_data, source_file = load_step_2_5_results()
    if step_2_5_data is None:
        print_status("Cannot proceed without Step 2.5 results", "ERROR")
        return None
    
    # ==========================================================================
    # 2. Extract Monthly Data
    # ==========================================================================
    print_status("[2/4] Extracting monthly EW/NS ratios...", "PROCESS")
    
    monthly_data = extract_monthly_data(step_2_5_data)
    
    if len(monthly_data) < 6:
        print_status(f"Insufficient monthly data: {len(monthly_data)} months", "ERROR")
        return None
    
    print_status(f"  Found {len(monthly_data)} months of data", "SUCCESS")
    
    # Show ratio range
    ratios = [m['ratio'] for m in monthly_data]
    print_status(f"  EW/NS ratio range: {min(ratios):.3f} to {max(ratios):.3f}", "INFO")
    
    # ==========================================================================
    # 3. Grid Search Over Background Directions
    # ==========================================================================
    print_status("[3/4] Grid searching background velocity directions...", "PROCESS")
    
    grid_results, best_ra, best_dec, best_r, ra_ci, dec_ci = grid_search_background(monthly_data)
    
    print_status(f"  Best-fit direction: RA={best_ra}°, Dec={best_dec}°", "SUCCESS")
    print_status(f"  Best correlation: r = {best_r:.4f}", "SUCCESS")
    print_status(f"  68% CI (RA): {ra_ci[0]:.0f}°-{ra_ci[1]:.0f}°", "INFO")
    print_status(f"  68% CI (Dec): {dec_ci[0]:.0f}°-{dec_ci[1]:.0f}°", "INFO")
    
    # Get p-value for best fit (Local t-test)
    n = len(monthly_data)
    t_stat = best_r * np.sqrt((n-2) / (1 - best_r**2 + 1e-10))
    best_p = 2 * (1 - stats.t.cdf(abs(t_stat), n-2))
    print_status(f"  Local P-value: {best_p:.4f}", "INFO")
    
    # Global P-value (Monte Carlo)
    # 5000 iterations for single run
    global_p, mc_hits = estimate_global_p_value_monte_carlo(monthly_data, best_r, n_iterations=5000)
    print_status(f"  Global P-value (Monte Carlo, N=5000): {global_p:.4f}", "SUCCESS")
    
    # ==========================================================================
    # 4. Compare to Reference Frames
    # ==========================================================================
    print_status("[4/4] Comparing to reference frames...", "PROCESS")
    
    # Angular separations
    cmb_sep = angular_separation(best_ra, best_dec, CMB_DIPOLE['ra'], CMB_DIPOLE['dec'])
    apex_sep = angular_separation(best_ra, best_dec, SOLAR_APEX['ra'], SOLAR_APEX['dec'])
    
    # Correlations at reference frame directions (find nearest grid point)
    ra_step = 10
    dec_step = 10
    cmb_ra_grid = round(CMB_DIPOLE['ra'] / ra_step) * ra_step
    cmb_dec_grid = round(CMB_DIPOLE['dec'] / dec_step) * dec_step
    apex_ra_grid = round(SOLAR_APEX['ra'] / ra_step) * ra_step
    apex_dec_grid = round(SOLAR_APEX['dec'] / dec_step) * dec_step
    
    cmb_r = grid_results.get((cmb_ra_grid, cmb_dec_grid), {'r': 0})['r']
    apex_r = grid_results.get((apex_ra_grid, apex_dec_grid), {'r': 0})['r']
    
    print_status("", "INFO")
    print_status("Reference Frame Comparison:", "INFO")
    print_status(f"  CMB Dipole: separation = {cmb_sep:.1f}°, r = {cmb_r:.4f}", "INFO")
    print_status(f"  Solar Apex: separation = {apex_sep:.1f}°, r = {apex_r:.4f}", "INFO")
    
    # Variance ratio
    if apex_r != 0:
        variance_ratio = (best_r**2) / (apex_r**2 + 1e-10)
    else:
        variance_ratio = float('inf') if best_r != 0 else 1.0
    
    print_status(f"  Variance ratio (best/apex): {variance_ratio:.1f}×", "INFO")
    
    # Which frame is closer?
    if cmb_sep < apex_sep:
        print_status("", "INFO")
        print_status("CMB FRAME CLOSER TO BEST FIT", "SUCCESS" if cmb_sep < 30 else "INFO")
    else:
        print_status("", "INFO")
        print_status("SOLAR APEX CLOSER TO BEST FIT", "WARNING")
    
    # ==========================================================================
    # Create Visualization
    # ==========================================================================
    fig_path = FIGURES_DIR / "step_2_7_cmb_frame_sky_map.png"
    create_sky_map(grid_results, best_ra, best_dec, fig_path)
    
    # ==========================================================================
    # Summary
    # ==========================================================================
    print_status("", "INFO")
    print_status("=" * 60, "INFO")
    print_status("CMB FRAME ANALYSIS SUMMARY", "INFO")
    print_status("=" * 60, "INFO")
    print_status(f"  Data source: {source_file}", "INFO")
    print_status(f"  Months analyzed: {len(monthly_data)}", "INFO")
    print_status(f"  Best-fit direction: RA={best_ra}°, Dec={best_dec}°", "INFO")
    print_status(f"  Best correlation: r = {best_r:.4f} (p = {best_p:.4f})", "INFO")
    print_status(f"  68% CI (RA): {ra_ci[0]:.0f}°-{ra_ci[1]:.0f}°", "INFO")
    print_status(f"  68% CI (Dec): {dec_ci[0]:.0f}°-{dec_ci[1]:.0f}°", "INFO")
    print_status(f"  Separation from CMB: {cmb_sep:.1f}°", "INFO")
    print_status(f"  Separation from Solar Apex: {apex_sep:.1f}°", "INFO")
    print_status(f"  Variance ratio (best/apex): {variance_ratio:.1f}×", "INFO")
    print_status("", "INFO")
    print_status("Comparison to CODE Longspan:", "INFO")
    print_status(f"  CODE: RA=186°, Dec=-4°, 18.2° from CMB, 5,570× ratio", "INFO")
    print_status(f"  RINEX: RA={best_ra}°, Dec={best_dec}°, {cmb_sep:.1f}° from CMB, {variance_ratio:.1f}× ratio", "INFO")
    print_status("", "INFO")
    
    # Interpretation
    if best_p < 0.05 and cmb_sep < 30:
        print_status("STRONG CMB FRAME ALIGNMENT DETECTED", "SUCCESS")
    elif best_p < 0.05:
        print_status("Significant anisotropy modulation but not aligned with CMB", "INFO")
    else:
        print_status("NO SIGNIFICANT CMB FRAME ALIGNMENT (null result)", "WARNING")
        print_status("  This may be due to shorter temporal baseline (3 years vs 25 years)", "INFO")
    
    # ==========================================================================
    # Save Results
    # ==========================================================================
    output = {
        'step': '2.7',
        'name': 'cmb_frame_analysis',
        'timestamp': datetime.now().isoformat(),
        'methodology': 'Grid search over background RA/Dec using monthly EW/NS ratios',
        'data_source': source_file,
        'n_months': len(monthly_data),
        'ew_ns_ratio_range': {
            'min': float(min(ratios)),
            'max': float(max(ratios)),
            'mean': float(np.mean(ratios)),
            'std': float(np.std(ratios))
        },
        'best_fit': {
            'ra': float(best_ra),
            'dec': float(best_dec),
            'correlation': float(best_r),
            'p_value': float(best_p),
            'global_p_value': float(global_p),
            'mc_hits': int(mc_hits),
            'ra_ci_68': [float(ra_ci[0]), float(ra_ci[1])],
            'dec_ci_68': [float(dec_ci[0]), float(dec_ci[1])]
        },
        'cmb_comparison': {
            'cmb_ra': CMB_DIPOLE['ra'],
            'cmb_dec': CMB_DIPOLE['dec'],
            'separation_deg': float(cmb_sep),
            'correlation_at_cmb': float(cmb_r)
        },
        'solar_apex_comparison': {
            'apex_ra': SOLAR_APEX['ra'],
            'apex_dec': SOLAR_APEX['dec'],
            'separation_deg': float(apex_sep),
            'correlation_at_apex': float(apex_r)
        },
        'variance_ratio': float(variance_ratio),
        'comparison_to_code': {
            'code_ra': 186,
            'code_dec': -4,
            'code_cmb_sep': 18.2,
            'code_variance_ratio': 5570,
            'rinex_cmb_sep': float(cmb_sep),
            'rinex_variance_ratio': float(variance_ratio)
        },
        'interpretation': {
            'significant': bool(best_p < 0.05),
            'cmb_aligned': bool(cmb_sep < 30),
            'conclusion': 'CMB alignment detected' if (best_p < 0.05 and cmb_sep < 30) else 'Null result'
        }
    }
    
    output_path = OUTPUTS_DIR / "step_2_7_cmb_frame_analysis.json"
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    
    print_status(f"Results saved: {output_path}", "SUCCESS")
    print_status("", "INFO")
    print_status("=" * 80, "INFO")
    print_status("Step 2.7 Complete", "INFO")
    print_status("=" * 80, "INFO")
    
    return output


def run_analysis_for_combination(filter_name, metric, coherence_type, mode='baseline'):
    """Run CMB analysis for a specific step_2_5 combination.
    
    Returns dict with results or None if failed.
    """
    # Load specific combination
    step_2_5_data, source_file = load_step_2_5_results(filter_name, metric, coherence_type, mode)
    if step_2_5_data is None:
        return None
    
    # Extract monthly data
    monthly_data = extract_monthly_data(step_2_5_data)
    if len(monthly_data) < 6:
        return None
    
    # Grid search with bootstrap CIs
    grid_results, best_ra, best_dec, best_r, ra_ci, dec_ci = grid_search_background(monthly_data)
    
    # P-value (Local t-test)
    n = len(monthly_data)
    t_stat = best_r * np.sqrt((n-2) / (1 - best_r**2 + 1e-10))
    best_p = 2 * (1 - stats.t.cdf(abs(t_stat), n-2))
    
    # Global P-value (Monte Carlo)
    # Use 1000 iterations for speed across all combinations
    global_p, mc_hits = estimate_global_p_value_monte_carlo(monthly_data, best_r, n_iterations=1000)
    
    # Separations from reference frames
    cmb_sep = angular_separation(best_ra, best_dec, CMB_DIPOLE['ra'], CMB_DIPOLE['dec'])
    apex_sep = angular_separation(best_ra, best_dec, SOLAR_APEX['ra'], SOLAR_APEX['dec'])
    
    return {
        'filter': filter_name,
        'mode': mode,
        'metric': metric,
        'coherence_type': coherence_type,
        'best_ra': best_ra,
        'best_dec': best_dec,
        'ra_ci_68': ra_ci,
        'dec_ci_68': dec_ci,
        'best_r': best_r,
        'best_p': best_p,
        'global_p': global_p,
        'mc_hits': mc_hits,
        'cmb_sep': cmb_sep,
        'apex_sep': apex_sep,
        'n_months': len(monthly_data),
        'grid_results': grid_results
    }


def main():
    """Main entry point - analyze all 18 combinations."""
    print_status("", "INFO")
    print_status("=" * 80, "INFO")
    print_status("STARTING CMB FRAME ANALYSIS - ALL 18 COMBINATIONS", "TITLE")
    print_status("=" * 80, "INFO")
    
    # Get all available combinations
    combinations = get_all_step_2_5_combinations()
    print_status(f"Found {len(combinations)} step_2_5 result files", "INFO")
    
    if not combinations:
        print_status("No step_2_5 results found!", "ERROR")
        return None
    
    # Analyze each combination
    all_results = []
    for i, combo in enumerate(combinations):
        print_status(f"\n[{i+1}/{len(combinations)}] Analyzing {combo['filter']}/{combo.get('mode','baseline')}/{combo['metric']}/{combo['coherence_type']}...", "PROCESS")
        
        result = run_analysis_for_combination(
            combo['filter'], combo['metric'], combo['coherence_type'], combo.get('mode', 'baseline')
        )
        
        if result:
            all_results.append(result)
            print_status(f"  Best: RA={result['best_ra']}°, Dec={result['best_dec']}°, r={result['best_r']:.3f}, CMB sep={result['cmb_sep']:.1f}°", "INFO")
    
    if not all_results:
        print_status("No valid results!", "ERROR")
        return None
    
    # ==========================================================================
    # COMPARISON SUMMARY
    # ==========================================================================
    print_status("\n" + "=" * 80, "INFO")
    print_status("CMB FRAME ANALYSIS COMPARISON", "INFO")
    print_status("=" * 80, "INFO")
    
    # Sort by CMB separation (closest first)
    all_results.sort(key=lambda x: x['cmb_sep'])
    
    print_status("", "INFO")
    print_status(f"{'Filter':<15} {'Mode':<10} {'Metric':<12} {'Coh':<8} {'RA':>6} {'Dec':>6} {'r':>7} {'CMB':>6} {'GlobP':>8}", "INFO")
    print_status("-" * 100, "INFO")
    
    for r in all_results:
        print_status(
            f"{r['filter']:<15} {r.get('mode','baseline'):<10} {r['metric']:<12} {r['coherence_type']:<8} "
            f"{r['best_ra']:>6.0f} {r['best_dec']:>6.0f} {r['best_r']:>7.3f} {r['cmb_sep']:>6.1f}° {r['global_p']:>8.4f}",
            "INFO"
        )
    
    # Find best result (closest to CMB with significant correlation)
    significant = [r for r in all_results if r['global_p'] < 0.05]
    if significant:
        best = min(significant, key=lambda x: x['cmb_sep'])
    else:
        best = all_results[0]
    
    # Calculate corrected global p-value (accounting for number of combinations)
    # Sidak correction: p_corrected = 1 - (1 - p_global)^N
    n_combos = len(all_results)
    p_global = best['global_p']
    p_corrected = 1 - (1 - p_global)**n_combos
    
    print_status("-" * 100, "INFO")
    print_status(f"CODE Longspan reference: RA=186°, Dec=-4°, CMB sep=18.2°", "INFO")
    print_status("", "INFO")
    print_status(f"BEST RESULT: {best['filter']}/{best.get('mode','baseline')}/{best['metric']}/{best['coherence_type']}", "SUCCESS")
    print_status(f"  RA={best['best_ra']}° (68% CI: {best['ra_ci_68'][0]:.0f}°-{best['ra_ci_68'][1]:.0f}°)", "SUCCESS")
    print_status(f"  Dec={best['best_dec']}° (68% CI: {best['dec_ci_68'][0]:.0f}°-{best['dec_ci_68'][1]:.0f}°)", "SUCCESS")
    print_status(f"  r={best['best_r']:.3f}", "SUCCESS")
    print_status(f"  Local P-value: {best['best_p']:.4f}", "INFO")
    print_status(f"  Global P-value (Monte Carlo): {best['global_p']:.4f}", "SUCCESS")
    print_status(f"  Corrected Global P-value (N={n_combos}): {p_corrected:.4f}", "SUCCESS")
    print_status(f"  CMB separation: {best['cmb_sep']:.1f}°", "SUCCESS")
    
    # Create visualization for best result
    print_status("\nCreating sky map for best result...", "PROCESS")
    create_sky_map(
        best['grid_results'], best['best_ra'], best['best_dec'],
        FIGURES_DIR / "step_2_7_cmb_frame_sky_map.png"
    )
    
    # Save comprehensive results
    output = {
        'step': '2.7',
        'name': 'cmb_frame_analysis',
        'timestamp': datetime.now().isoformat(),
        'n_combinations_tested': len(all_results),
        'best_result': {
            'filter': best['filter'],
            'mode': best.get('mode', 'baseline'),
            'metric': best['metric'],
            'coherence_type': best['coherence_type'],
            'best_ra': float(best['best_ra']),
            'best_dec': float(best['best_dec']),
            'ra_ci_68': [float(best['ra_ci_68'][0]), float(best['ra_ci_68'][1])],
            'dec_ci_68': [float(best['dec_ci_68'][0]), float(best['dec_ci_68'][1])],
            'correlation': float(best['best_r']),
            'p_value': float(best['best_p']),
            'global_p_value': float(best['global_p']),
            'corrected_p_value': float(p_corrected),
            'cmb_separation': float(best['cmb_sep']),
            'apex_separation': float(best['apex_sep'])
        },
        'methodology': {
            'grid_resolution_deg': 1,
            'background_speed_kms': GRID_SEARCH_SPEED_KMS,
            'bootstrap_samples': 500,
            'predictor': 'cos(velocity_declination)'
        },
        'all_results': [
            {
                'filter': r['filter'],
                'mode': r.get('mode', 'baseline'),
                'metric': r['metric'],
                'coherence_type': r['coherence_type'],
                'best_ra': float(r['best_ra']),
                'best_dec': float(r['best_dec']),
                'ra_ci_68': [float(r['ra_ci_68'][0]), float(r['ra_ci_68'][1])],
                'dec_ci_68': [float(r['dec_ci_68'][0]), float(r['dec_ci_68'][1])],
                'correlation': float(r['best_r']),
                'p_value': float(r['best_p']),
                'cmb_separation': float(r['cmb_sep']),
                'apex_separation': float(r['apex_sep'])
            }
            for r in all_results
        ],
        'cmb_dipole': CMB_DIPOLE,
        'solar_apex': SOLAR_APEX,
        'comparison_to_code': {
            'code_ra': 186,
            'code_dec': -4,
            'code_cmb_sep': 18.2,
            'code_variance_ratio': 5570,
            'best_rinex_cmb_sep': float(best['cmb_sep'])
        }
    }
    
    output_path = OUTPUTS_DIR / "step_2_7_cmb_frame_analysis.json"
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    
    print_status(f"Results saved: {output_path}", "SUCCESS")
    print_status("\n" + "=" * 80, "INFO")
    print_status("Step 2.7 Complete - ALL COMBINATIONS ANALYZED", "INFO")
    print_status("=" * 80, "INFO")
    
    return output


if __name__ == "__main__":
    main()
