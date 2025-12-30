#!/usr/bin/env python3
"""
TEP-GNSS-RINEX Data Alignment Utilities
=======================================

Time series alignment and spectral analysis functions for GNSS clock data.
Handles gap-aware loading, temporal alignment, and cross-spectral density
computation across multi-day datasets.

Key Functions:
    load_aligned_data: Load and align .npz files into pandas DataFrame
    compute_aligned_coherence: Compute CSD/coherence for aligned series

The spectral analysis methodology matches CODE longspan for cross-validation:
    1. Manual linear detrend followed by constant detrend in CSD
    2. Welch's method with 1024-sample segments (50% overlap)
    3. Magnitude-weighted circular phase averaging

Author: Matthew Lukin Smawfield
Date: 9 December 2025
License: CC-BY-4.0
"""

import numpy as np
import pandas as pd
from scipy.signal import csd, welch
from datetime import datetime, timedelta


def load_aligned_data(files, year=None, keys=['clock_bias_ns'], sampling_sec=None):
    """
    Load data from a list of .npz files into a Pandas DataFrame with DatetimeIndex.
    Supports multi-year loading and automatic gap handling.
    
    Parameters:
    -----------
    files : list of Path or str
        List of .npz files for a station.
    year : int or None
        Specific year to load. If None, loads all available years.
    keys : list of str
        List of keys to load from .npz files.
    sampling_sec : float, optional
        Sampling interval in seconds. If None, auto-detect from first file.
        
    Returns:
    --------
    full_series : pd.DataFrame or pd.Series
        Aligned data with DatetimeIndex.
    """
    # Auto-detect sampling from first file if not specified
    if sampling_sec is None:
        for f in files:
            try:
                with np.load(f) as d:
                    if 'timestamps' in d:
                        ts = d['timestamps']
                        if len(ts) > 1:
                            diffs = np.diff(ts[ts >= 0])
                            if len(diffs) > 0:
                                sampling_sec = float(np.median(diffs))
                                break
            except Exception:
                continue
        if sampling_sec is None:
            sampling_sec = 300.0
            
    frames = []
    
    for f in files:
        try:
            fname = str(f).split('/')[-1]
            parts = fname.split('_')
            if len(parts) < 2: continue
            
            doy_str = parts[1]
            if len(doy_str) < 7: continue
            
            file_year = int(doy_str[:4])
            file_doy = int(doy_str[4:7])
            
            if year is not None and file_year != year:
                continue
                
            with np.load(f) as d:
                if 'timestamps' not in d: continue
                
                ts = d['timestamps']
                valid_mask = (ts >= 0) & (ts < 86400)
                if not np.any(valid_mask): continue
                
                # Create datetime index for this day
                base_date = datetime(file_year, 1, 1) + timedelta(days=file_doy-1)
                
                # Convert seconds to timestamps
                # Round to nearest sampling interval to ensure alignment
                ts_valid = ts[valid_mask]
                ts_rounded = np.round(ts_valid / sampling_sec) * sampling_sec
                
                times = [base_date + timedelta(seconds=float(t)) for t in ts_rounded]
                
                data_dict = {}
                for k in keys:
                    if k in d:
                        val = d[k]
                        if len(val) == len(valid_mask):
                            data_dict[k] = val[valid_mask]
                            
                if data_dict:
                    df = pd.DataFrame(data_dict, index=times)
                    # Remove duplicates if any (shouldn't happen in standard processing)
                    df = df[~df.index.duplicated(keep='first')]
                    frames.append(df)
                    
        except Exception as e:
            continue
            
    if not frames:
        # Return empty structure with correct columns
        return pd.DataFrame(columns=keys) if len(keys) > 1 else pd.Series(dtype=float)
        
    # Concatenate all days and sort by time
    full_df = pd.concat(frames).sort_index()
    
    # Remove any duplicate timestamps that may exist across files
    full_df = full_df[~full_df.index.duplicated(keep='first')]
    
    if len(keys) == 1:
        return full_df[keys[0]]
    return full_df

def compute_aligned_coherence(s1, s2, fs, nperseg=1024, detrend='constant'):
    """
    Compute CSD and coherence for two aligned Pandas Series.
    Finds overlapping valid segments automatically.
    """
    # Ensure no duplicate indices before alignment
    if s1.index.has_duplicates:
        s1 = s1[~s1.index.duplicated(keep='first')]
    if s2.index.has_duplicates:
        s2 = s2[~s2.index.duplicated(keep='first')]
    
    # Align data on time index (intersection)
    # This automatically handles gaps: only times present in BOTH are kept
    common = pd.concat([s1, s2], axis=1, join='inner').dropna()
    
    if len(common) < nperseg:
        return None, None, None, None
        
    vals1 = common.iloc[:, 0].values
    vals2 = common.iloc[:, 1].values
    times = common.index
    
    # Find continuous segments based on sampling rate
    # A gap exists if time difference > 1.5 * sampling_period
    dt_seconds = np.diff(times.astype(np.int64)) // 10**9
    sampling_period = 1.0 / fs
    gap_indices = np.where(dt_seconds > 1.5 * sampling_period)[0]
    
    # Add start and end indices
    segment_starts = np.concatenate(([0], gap_indices + 1))
    segment_ends = np.concatenate((gap_indices + 1, [len(common)]))
    
    Pxy_sum = None
    Pxx_sum = None
    Pyy_sum = None
    total_weight = 0
    freqs = None
    
    for s, e in zip(segment_starts, segment_ends):
        # Skip short segments
        if e - s < nperseg:
            continue
            
        seg1 = vals1[s:e]
        seg2 = vals2[s:e]
        
        # Detrending - ALIGNED WITH CODE LONGSPAN METHODOLOGY
        # CODE does: (1) manual linear detrend, then (2) constant detrend in CSD
        if detrend == 'linear':
            # Step 1: Manual linear detrend (matches CODE line 1725-1726)
            x = np.arange(len(seg1))
            seg1 = seg1 - np.polyval(np.polyfit(x, seg1, 1), x)
            seg2 = seg2 - np.polyval(np.polyfit(x, seg2, 1), x)
            # Step 2: CSD will apply constant detrend (matches CODE line 1747)
            csd_detrend = 'constant'
        elif detrend == 'constant':
            seg1 = seg1 - np.mean(seg1)
            seg2 = seg2 - np.mean(seg2)
            csd_detrend = False
        else:
            csd_detrend = False
            
        try:
            f, Pxy = csd(seg1, seg2, fs=fs, nperseg=nperseg, detrend=csd_detrend)
            _, Pxx = welch(seg1, fs=fs, nperseg=nperseg, detrend=csd_detrend)
            _, Pyy = welch(seg2, fs=fs, nperseg=nperseg, detrend=csd_detrend)
            
            # Weight by number of windows
            n_windows = (len(seg1) - nperseg) // (nperseg // 2) + 1
            if n_windows < 1: n_windows = 1
            
            if Pxy_sum is None:
                Pxy_sum = Pxy * n_windows
                Pxx_sum = Pxx * n_windows
                Pyy_sum = Pyy * n_windows
                freqs = f
            else:
                Pxy_sum += Pxy * n_windows
                Pxx_sum += Pxx * n_windows
                Pyy_sum += Pyy * n_windows
            
            total_weight += n_windows
            
        except Exception:
            continue
            
    if total_weight == 0 or Pxy_sum is None:
        return None, None, None, None
        
    return freqs, Pxy_sum / total_weight, Pxx_sum / total_weight, Pyy_sum / total_weight
