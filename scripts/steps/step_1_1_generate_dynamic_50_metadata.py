#!/usr/bin/env python3
"""
TEP-GNSS-RINEX - STEP 1.1: Generate DYNAMIC_50 Station Selection
================================================================

Generates dynamic_50_metadata.json from NPZ files using quality filtering.
This creates an optimized subset of stations with the cleanest clock data,
ensuring high-quality input for correlation analysis.

Quality Criteria:
    - Clock bias standard deviation < 50 ns
    - No single-epoch jumps > 500 ns
    - Total range < 5000 ns

Output: data/processed/dynamic_50_metadata.json
    Station metadata with coordinates and quality metrics

Author: Matthew Lukin Smawfield
Date: 9 December 2025
License: CC-BY-4.0
"""

import sys
import json
import numpy as np
from pathlib import Path
from collections import defaultdict

# Paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status

# Initialize logger
logger = TEPLogger(
    name="step_1_1_generate_dynamic_50_metadata",
    level="INFO",
    log_file_path=PROJECT_ROOT / "logs" / "step_1_1_generate_dynamic_50_metadata.log"
)
set_step_logger(logger)
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
OUTPUT_FILE = PROCESSED_DIR / "dynamic_50_metadata.json"

# Quality thresholds (strict - same as step_2_0 dynamic50)
STD_THRESHOLD = 50.0  # nanoseconds (stricter = cleaner data)
JUMP_THRESHOLD = 500.0  # nanoseconds
RANGE_THRESHOLD = 5000.0  # nanoseconds


def analyze_npz_quality(npz_path):
    """Analyze NPZ file and return quality metrics."""
    try:
        data = np.load(npz_path, allow_pickle=True)
        if 'clock_bias_ns' not in data:
            return None
        
        clock = data['clock_bias_ns']
        if len(clock) < 100:  # Need minimum data
            return None
        
        # Calculate metrics
        std = np.nanstd(clock)
        total_range = np.nanmax(clock) - np.nanmin(clock)
        
        # Check for jumps
        diffs = np.abs(np.diff(clock))
        max_jump = np.nanmax(diffs) if len(diffs) > 0 else 0
        
        return {
            'std': std,
            'range': total_range,
            'max_jump': max_jump,
            'n_epochs': len(clock)
        }
    except Exception as e:
        return None


def get_station_hemisphere(lat):
    """Determine hemisphere from latitude."""
    return 'N' if lat >= 0 else 'S'


def main():
    print_status("GENERATING DYNAMIC_50_METADATA.JSON", level="TITLE")
    print_status(f"Quality criteria: std < {STD_THRESHOLD}ns, jumps < {JUMP_THRESHOLD}ns, range < {RANGE_THRESHOLD}ns")
    
    # Load station coordinates
    coords_file = PROCESSED_DIR / "station_coordinates.json"
    if not coords_file.exists():
        print_status(f"ERROR: {coords_file} not found", level="ERROR")
        return
    
    with open(coords_file) as f:
        coords_ecef = json.load(f)
    
    # Convert ECEF to LLA for hemisphere info
    def ecef_to_lla(x, y, z):
        a = 6378137.0
        f = 1/298.257223563
        e2 = 2*f - f**2
        lon = np.arctan2(y, x)
        p = np.sqrt(x**2 + y**2)
        lat = np.arctan2(z, p * (1 - e2))
        for _ in range(10):
            N = a / np.sqrt(1 - e2 * np.sin(lat)**2)
            lat = np.arctan2(z + e2 * N * np.sin(lat), p)
        h = p / np.cos(lat) - N
        return np.degrees(lat), np.degrees(lon), h
    
    station_lla = {}
    for sta, ecef in coords_ecef.items():
        lat, lon, alt = ecef_to_lla(ecef[0], ecef[1], ecef[2])
        station_lla[sta] = {'lat': lat, 'lon': lon, 'alt': alt}
    
    # Find all NPZ files
    npz_files = sorted(PROCESSED_DIR.glob("*.npz"))
    print_status(f"Found {len(npz_files)} NPZ files")
    
    # Track quality per station
    station_quality = defaultdict(list)
    
    rejected_jump = 0
    rejected_range = 0
    rejected_std = 0
    passed = 0
    
    for i, npz_path in enumerate(npz_files):
        if (i + 1) % 100 == 0:
            print(f"  Processing {i+1}/{len(npz_files)}...", end='\r')  # Progress indicator
        
        station = npz_path.stem.split('_')[0]
        quality = analyze_npz_quality(npz_path)
        
        if quality is None:
            continue
        
        # Apply quality filters
        if quality['max_jump'] > JUMP_THRESHOLD:
            rejected_jump += 1
            continue
        if quality['range'] > RANGE_THRESHOLD:
            rejected_range += 1
            continue
        if quality['std'] > STD_THRESHOLD:
            rejected_std += 1
            continue
        
        passed += 1
        station_quality[station].append(quality)
    
    print()  # Clear progress line
    print_status(f"Quality filtering results:")
    print_status(f"  Passed: {passed}")
    print_status(f"  Rejected (jumps > {JUMP_THRESHOLD}ns): {rejected_jump}")
    print_status(f"  Rejected (range > {RANGE_THRESHOLD}ns): {rejected_range}")
    print_status(f"  Rejected (std > {STD_THRESHOLD}ns): {rejected_std}")
    print_status(f"  Unique stations passing: {len(station_quality)}")
    
    # Build metadata
    stations_meta = {}
    for station, qualities in station_quality.items():
        if station not in station_lla:
            continue
        
        lla = station_lla[station]
        avg_std = np.mean([q['std'] for q in qualities])
        avg_range = np.mean([q['range'] for q in qualities])
        n_files = len(qualities)
        
        stations_meta[station] = {
            'lat': float(lla['lat']),
            'lon': float(lla['lon']),
            'alt': float(lla['alt']),
            'avg_std_ns': float(avg_std),
            'avg_range_ns': float(avg_range),
            'n_files_passed': n_files,
            'classification': {
                'hemisphere': get_station_hemisphere(lla['lat']),
                'quality': 'high' if avg_std < 50 else 'medium'
            }
        }
    
    # Count hemispheres
    n_north = sum(1 for s in stations_meta.values() if s['classification']['hemisphere'] == 'N')
    n_south = sum(1 for s in stations_meta.values() if s['classification']['hemisphere'] == 'S')
    
    print_status(f"Final station count: {len(stations_meta)}")
    print_status(f"  Northern hemisphere: {n_north}")
    print_status(f"  Southern hemisphere: {n_south}")
    
    # Save metadata
    output = {
        'name': 'DYNAMIC_50',
        'description': f'Stations with clock std < {STD_THRESHOLD}ns, no jumps > {JUMP_THRESHOLD}ns, range < {RANGE_THRESHOLD}ns',
        'quality_criteria': {
            'std_threshold_ns': STD_THRESHOLD,
            'jump_threshold_ns': JUMP_THRESHOLD,
            'range_threshold_ns': RANGE_THRESHOLD
        },
        'n_stations': len(stations_meta),
        'hemisphere_balance': f'{n_north}:{n_south}',
        'stations': stations_meta
    }
    
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(output, f, indent=2)
    
    print_status(f"Saved: {OUTPUT_FILE}", level="SUCCESS")


if __name__ == '__main__':
    main()
