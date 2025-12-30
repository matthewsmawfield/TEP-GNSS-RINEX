"""
TEP-GNSS-RINEX Analysis Pipeline - Utility Modules
===================================================

This package provides shared utilities for the TEP-GNSS-RINEX analysis pipeline.

Modules:
    config: Centralized configuration management with environment variable support
    logger: Standardized logging with color-coded console output and file logging
    data_alignment: Time series alignment and spectral analysis functions
    tep_utils: Core mathematical functions for TEP correlation analysis
    fetch_kp: Geomagnetic Kp index data retrieval from GFZ Potsdam

Author: Matthew Lukin Smawfield
Date: 9 December 2025
License: CC-BY-4.0
"""

from .config import TEPConfig, config
from .logger import TEPLogger, print_status, set_step_logger
