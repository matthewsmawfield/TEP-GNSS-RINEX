#!/usr/bin/env python3
"""
TEP-GNSS-RINEX Analysis Configuration
=====================================

Centralized configuration management for the TEP-GNSS-RINEX analysis pipeline.
Supports environment variable overrides for all parameters, enabling flexible
deployment across local machines and cloud compute instances.

Configuration Hierarchy:
    1. Environment variables (highest priority)
    2. TEPConfig.DEFAULTS dictionary
    3. Function-level defaults (lowest priority)

Environment Variables:
    TEP_BINS: Number of distance bins for correlation analysis (default: 40)
    TEP_MAX_DISTANCE_KM: Maximum station pair distance in km (default: 13000)
    TEP_MIN_BIN_COUNT: Minimum pairs per bin for fitting (default: 50)
    TEP_BOOTSTRAP_ITER: Bootstrap iterations for uncertainty (default: 5000)
    TEP_NULL_ITERATIONS: Null hypothesis test iterations (default: 20)
    TEP_WORKERS: Number of parallel workers (default: CPU count)
    TEP_MIN_STATIONS: Minimum stations required (default: 0)

Example:
    >>> from scripts.utils.config import config
    >>> n_bins = config.get_int('TEP_BINS')
    >>> max_dist = config.get_float('TEP_MAX_DISTANCE_KM')

Author: Matthew Lukin Smawfield
Date: 9 December 2025
License: CC-BY-4.0
"""

import os
from typing import Optional, Union, List
from pathlib import Path
import multiprocessing as mp
import numpy as np


class TEPConfig:
    """Centralized configuration management for TEP GNSS Analysis"""
    
    # Default values for common configuration parameters
    DEFAULTS = {
        # Analysis parameters
        'TEP_BINS': 40,
        'TEP_MAX_DISTANCE_KM': 13000.0,
        'TEP_MIN_BIN_COUNT': 50,
        'TEP_BOOTSTRAP_ITER': 5000,
        'TEP_NULL_ITERATIONS': 20,
        
        # Processing parameters
        'TEP_WORKERS': 14,
        
        # Data parameters
        'TEP_MIN_STATIONS': 0,
    }

    @staticmethod
    def get_int(key: str, default: Optional[int] = None) -> int:
        if default is None:
            default = TEPConfig.DEFAULTS.get(key)
        value = os.getenv(key)
        if value is None:
            if default is None:
                raise ValueError(f"Required configuration {key} not found")
            return default
        try:
            return int(value)
        except ValueError:
            return default
    
    @staticmethod
    def get_float(key: str, default: Optional[float] = None) -> float:
        if default is None:
            default = TEPConfig.DEFAULTS.get(key)
        value = os.getenv(key)
        if value is None:
            if default is None:
                raise ValueError(f"Required configuration {key} not found")
            return default
        try:
            return float(value)
        except ValueError:
            return default
    
    @staticmethod
    def get_bool(key: str, default: Optional[bool] = None) -> bool:
        if default is None:
            default = TEPConfig.DEFAULTS.get(key, False)
        value = os.getenv(key)
        if value is None:
            return default
        return str(value).lower() in ('1', 'true', 'yes', 'on')
    
    @staticmethod
    def get_str(key: str, default: Optional[str] = None) -> str:
        if default is None:
            default = TEPConfig.DEFAULTS.get(key)
        value = os.getenv(key)
        if value is None:
            if default is None:
                raise ValueError(f"Required configuration {key} not found")
            return default
        return value

    @staticmethod
    def get_worker_count(env_var: str = 'TEP_WORKERS') -> int:
        default_workers = mp.cpu_count()
        try:
            user_workers = int(os.getenv(env_var, default_workers))
            max_workers = min(32, default_workers * 2)
            return max(1, min(user_workers, max_workers))
        except (ValueError, TypeError):
            return default_workers

# Convenience instance
config = TEPConfig()
