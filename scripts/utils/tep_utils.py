#!/usr/bin/env python3
"""
TEP-GNSS-RINEX Core Mathematical Utilities
==========================================

Core mathematical functions for Temporal Equivalence Principle (TEP) analysis.
Provides exponential decay fitting, coordinate transformations, and binning
operations used throughout the analysis pipeline.

Key Functions:
    exp_decay: Exponential decay model C(r) = A*exp(-r/λ) + C₀
    fit_exponential_binned: Weighted least-squares fit to binned data
    bin_and_fit: Combined binning and fitting pipeline
    ecef_to_lla: ECEF to geodetic coordinate conversion
    load_station_map: Load station coordinate metadata

Constants:
    MIN_DISTANCE_KM: Minimum station pair distance (50 km)
    MAX_DISTANCE_KM: Maximum station pair distance (13,000 km)
    N_BINS: Number of logarithmic distance bins (40)
    MIN_BIN_COUNT: Minimum pairs per bin for fitting (50)

Author: Matthew Lukin Smawfield
Date: 9 December 2025
License: CC-BY-4.0
"""

import numpy as np
import json
import math
from scipy.optimize import curve_fit
from scipy import stats
import matplotlib.pyplot as plt
from pathlib import Path

# =============================================================================
# Constants
# =============================================================================
MIN_DISTANCE_KM = 50
MAX_DISTANCE_KM = 13000
N_BINS = 40
MIN_BIN_COUNT = 50

def exp_decay(r, A, lam, C0):
    """Exponential decay model: A * exp(-r/lambda) + C0"""
    return A * np.exp(-r / lam) + C0

def fit_exponential_binned(bin_centers, bin_means, bin_counts, p0=None):
    """
    Fit exponential decay to binned data.
    """
    if len(bin_centers) < 5:
        return None
    
    # Initial guess
    if p0 is None:
        p0 = [0.5, 2000, 0]
        
    try:
        # Weighted fit (sigma = 1/sqrt(N))
        # Points with more data have smaller sigma -> higher weight
        sigma = 1.0 / np.sqrt(bin_counts)
        
        # Bounds: A [0, 2], lambda [100, 20000], C0 [-1, 1]
        bounds = ([0, 100, -1], [2, 20000, 1])
        
        popt, pcov = curve_fit(
            exp_decay, bin_centers, bin_means, 
            p0=p0, sigma=sigma, bounds=bounds, maxfev=10000
        )
        
        A, lam, C0 = popt
        lambda_err = np.sqrt(pcov[1, 1])
        
        # Weighted R-squared (consistent with weighted fit)
        # Use bin_counts as weights (same as 1/sigma^2)
        weights = bin_counts.astype(float)
        predicted = exp_decay(bin_centers, *popt)
        weighted_mean = np.average(bin_means, weights=weights)
        ss_res = np.sum(weights * (bin_means - predicted)**2)
        ss_tot = np.sum(weights * (bin_means - weighted_mean)**2)
        r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0
        
        # Flag if parameters hit bounds
        boundary_hit = bool(A > 1.99 or A < 0.01 or 
                            C0 < -0.99 or C0 > 0.99 or 
                            lam > 19000 or lam < 101)
        
        return {
            'amplitude': float(A),
            'lambda_km': float(lam),
            'lambda_err': float(lambda_err),
            'offset': float(C0),
            'r_squared': float(r2),
            'n_pairs': int(np.sum(bin_counts)),
            'n_bins': len(bin_centers),
            'boundary_hit': boundary_hit
        }
    except Exception:
        return None

def bin_and_fit(distances, coherences, min_count=MIN_BIN_COUNT, n_bins=N_BINS, dist_range=(MIN_DISTANCE_KM, MAX_DISTANCE_KM)):
    """
    Bin raw data and fit exponential.
    """
    if len(distances) < 1000: return None
    
    distances = np.asarray(distances)
    coherences = np.asarray(coherences)
    
    # Logarithmic binning (better for exponential)
    bin_edges = np.logspace(np.log10(dist_range[0]), np.log10(dist_range[1]), n_bins + 1)
    
    bin_centers = []
    bin_means = []
    bin_counts = []
    
    # Vectorized binning is harder with variable edges, use digitize
    indices = np.digitize(distances, bin_edges) - 1
    
    for i in range(n_bins):
        mask = indices == i
        count = np.sum(mask)
        if count >= min_count:
            # Geometric mean center
            center = np.sqrt(bin_edges[i] * bin_edges[i+1])
            # Or arithmetic: (bin_edges[i] + bin_edges[i+1]) / 2
            
            bin_centers.append(center)
            bin_means.append(np.nanmean(coherences[mask]))
            bin_counts.append(count)
            
    return fit_exponential_binned(np.array(bin_centers), np.array(bin_means), np.array(bin_counts))

def ecef_to_lla(x, y, z):
    """Convert ECEF to Lat, Lon, Alt."""
    a = 6378137.0
    e2 = 0.00669437999014
    lon = math.atan2(y, x)
    p = math.sqrt(x**2 + y**2)
    lat = math.atan2(z, p * (1 - e2))
    # Iterative solution for lat
    for _ in range(5):
        N = a / math.sqrt(1 - e2 * math.sin(lat)**2)
        lat = math.atan2(z + e2 * N * math.sin(lat), p)
    
    N = a / math.sqrt(1 - e2 * math.sin(lat)**2)
    alt = p / math.cos(lat) - N
    return math.degrees(lat), math.degrees(lon), alt

def load_station_map(processed_dir):
    """Load station coordinates and return map."""
    path = Path(processed_dir) / "station_coordinates.json"
    if not path.exists(): return None
    with open(path) as f:
        data = json.load(f)
    return data
