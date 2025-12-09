#!/usr/bin/env python3
"""
TEP-GNSS-RINEX Geomagnetic Index Retrieval
==========================================

Downloads Kp geomagnetic index data from GFZ Helmholtz Centre Potsdam.
The Kp index is used for geomagnetic storm stratification in Step 2.3a,
allowing separation of quiet (Kp < 3) and disturbed (Kp >= 3) conditions.

Data Source: GFZ Potsdam Kp Index Service
    https://kp.gfz-potsdam.de/

Output: data/kp_index.json
    Daily mean Kp values keyed by YYYY_DOY format

Author: Matthew Lukin Smawfield
Date: 9 December 2025
License: CC-BY-4.0
"""

import requests
import json
import numpy as np
from datetime import datetime
from pathlib import Path

# GFZ Helmholtz Centre Potsdam - Kp Index
# Using their API/Download service
# URL format: https://kp.gfz-potsdam.de/app/json/\?download\=1\&start\=2022-01-01T00:00:00Z\&end\=2024-12-31T23:59:59Z

URL = "https://kp.gfz-potsdam.de/app/json/?download=1&start=2022-01-01T00:00:00Z&end=2024-12-31T23:59:59Z"
OUTPUT = Path("data/kp_index.json")

def main():
    print(f"Downloading Kp index from {URL}...")
    try:
        r = requests.get(URL, timeout=30)
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        raise RuntimeError(f"CRITICAL: Failed to download Kp data: {e}. "
                          f"Cannot proceed without real geomagnetic data.")

    # Parse
    timestamps = data.get('datetime', [])
    kp_values = data.get('Kp', [])
    
    if not timestamps:
        print("No data found")
        return

    # Aggregate by DOY
    daily_kp = {} # key: YYYY_DOY, val: list of Kp
    
    for ts_str, kp in zip(timestamps, kp_values):
        # ts_str: "2022-01-01T00:00:00Z"
        dt = datetime.strptime(ts_str, "%Y-%m-%dT%H:%M:%SZ")
        doy = dt.timetuple().tm_yday
        key = f"{dt.year}_{doy:03d}"
        
        if key not in daily_kp: daily_kp[key] = []
        daily_kp[key].append(kp)
        
    # Compute daily mean
    final_index = {}
    for key, vals in daily_kp.items():
        final_index[key] = round(float(np.mean(vals)), 2)
        
    print(f"Processed {len(final_index)} days.")
    
    with open(OUTPUT, 'w') as f:
        json.dump(final_index, f, indent=2)
    print(f"Saved to {OUTPUT}")

if __name__ == '__main__':
    main()
