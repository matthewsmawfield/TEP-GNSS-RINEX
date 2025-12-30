#!/usr/bin/env python3
"""
TEP-GNSS-RINEX Logging Utilities
================================

Standardized logging infrastructure for the TEP-GNSS-RINEX analysis pipeline.
Provides color-coded console output, file logging, and custom log levels for
consistent status reporting across all analysis steps.

Custom Log Levels:
    PROCESS (25): Blue - ongoing operations
    SUCCESS (26): Green - completed operations
    TEST (27): Magenta - test/validation results

Key Components:
    TEPLogger: Main logger class with console and file handlers
    TEPFormatter: Color-coded console formatter
    TEPFileFormatter: Clean formatter for log files (no ANSI codes)
    print_status: Global function for status messages

Usage:
    >>> from scripts.utils.logger import TEPLogger, print_status, set_step_logger
    >>> logger = TEPLogger("step_2_0", log_file_path=Path("logs/step_2_0.log"))
    >>> set_step_logger(logger)
    >>> print_status("Analysis complete", "SUCCESS")

Author: Matthew Lukin Smawfield
Date: 9 December 2025
License: CC-BY-4.0
"""

import logging
import sys
from pathlib import Path
import os
from typing import Optional

# Project root for relative path calculations
PACKAGE_ROOT = Path(__file__).resolve().parents[2]

class TEPFormatter(logging.Formatter):
    """Formatter with color support that matches step 3's print_status format."""

    # ANSI color codes matching step 3
    COLORS = {
        'SUCCESS': '\033[1;32m',  # Green bold
        'WARNING': '\033[1;33m',  # Yellow bold
        'ERROR': '\033[1;31m',    # Red bold
        'INFO': '\033[0;37m',     # White
        'DEBUG': '\033[0;90m',    # Dark gray
        'PROCESS': '\033[0;34m',  # Blue
        'TEST': '\033[1;35m',      # Magenta bold
        'TITLE': '\033[1;36m',     # Cyan bold for titles
        'CRITICAL': '\033[1;41m'  # White on red background for critical errors
    }
    RESET = '\033[0m'

    def __init__(self, fmt=None, datefmt=None, use_colors=True):
        # Use same timestamp format as step 3: %H:%M:%S
        super().__init__(fmt, datefmt='%H:%M:%S')
        self.use_colors = use_colors

    def format(self, record):
        message = record.getMessage()

        # Map log levels to display names and colors
        level_mapping = {
            25: ('PROCESS', self.COLORS['PROCESS']),   # Custom PROCESS level
            26: ('SUCCESS', self.COLORS['SUCCESS']),  # Custom SUCCESS level
            27: ('TEST', self.COLORS['TEST']),        # Custom TEST level
            logging.INFO: ('INFO', self.COLORS['INFO']),
            logging.WARNING: ('WARNING', self.COLORS['WARNING']),
            logging.ERROR: ('ERROR', self.COLORS['ERROR']),
            logging.DEBUG: ('DEBUG', self.COLORS['DEBUG']),
            logging.CRITICAL: ('CRITICAL', self.COLORS['CRITICAL']) # CRITICAL level
        }

        level_name, color = level_mapping.get(record.levelno, ('INFO', self.COLORS['INFO']))

        # Format exactly like step 3: [timestamp] [LEVEL] message
        timestamp = self.formatTime(record, self.datefmt)
        
        if self.use_colors:
            return f"{color}[{timestamp}] [{level_name}] {message}{self.RESET}"
        else:
            return f"[{timestamp}] [{level_name}] {message}"

class TEPFileFormatter(logging.Formatter):
    """Clean formatter for file output without ANSI color codes."""

    def __init__(self, fmt=None, datefmt=None):
        # Use same timestamp format as step 3: %H:%M:%S
        super().__init__(fmt, datefmt='%H:%M:%S')

    def format(self, record):
        message = record.getMessage()

        # Map log levels to display names (no colors for file)
        level_mapping = {
            25: 'PROCESS',   # Custom PROCESS level
            26: 'SUCCESS',  # Custom SUCCESS level
            27: 'TEST',     # Custom TEST level
            logging.INFO: 'INFO',
            logging.WARNING: 'WARNING',
            logging.ERROR: 'ERROR',
            logging.DEBUG: 'DEBUG',
            logging.CRITICAL: 'CRITICAL'
        }

        level_name = level_mapping.get(record.levelno, 'INFO')

        # Format exactly like step 3: [timestamp] [LEVEL] message (no colors)
        timestamp = self.formatTime(record, self.datefmt)
        return f"[{timestamp}] [{level_name}] {message}"

class TEPLogger:
    def __init__(self, name: str = "tep_gnss", level: str = "INFO", log_file_path: Optional[Path] = None, reset_log: bool = True):
        self.logger = logging.getLogger(name)
        self.logger.setLevel(self._get_log_level(level))
        
        # Clear any existing handlers to prevent duplication
        self.logger.handlers.clear()
        
        # Create a console handler with immediate flushing
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)  # Allow all custom levels
        # Ensure immediate output by setting stream to unbuffered
        ch.stream.reconfigure(line_buffering=True)

        # Create formatters
        console_formatter = TEPFormatter()  # With colors for console
        file_formatter = TEPFileFormatter()  # Clean format for file
        
        ch.setFormatter(console_formatter)

        self.logger.addHandler(ch)
        self.logger.propagate = False # Prevent messages from being passed to the root logger

        # Add file handler if log_file_path is provided, or use a default one
        if log_file_path is None:
            # Use a default log file if none is provided, e.g., a general_log.log in logs directory
            default_log_dir = PACKAGE_ROOT / "logs"
            default_log_dir.mkdir(parents=True, exist_ok=True)
            log_file_path = default_log_dir / "general_tep_gnss.log"
            reset_log = False  # Don't reset the general log file

        log_file_path.parent.mkdir(parents=True, exist_ok=True) # Ensure directory exists
        
        # Reset log file if requested (default for step-specific loggers)
        if reset_log:
            try:
                with open(log_file_path, 'w') as f:
                    f.write("")  # Clear the log file
            except Exception as e:
                print(f"Warning: Could not reset log file {log_file_path}: {e}")
        
        fh = logging.FileHandler(log_file_path, mode='a', encoding='utf-8')  # Use append mode after reset
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(file_formatter)  # Use clean formatter for file
        self.logger.addHandler(fh)

    def _get_log_level(self, level_name: str):
        return getattr(logging, level_name.upper(), logging.INFO)

    def info(self, message: str):
        self.logger.info(message)

    def warning(self, message: str):
        self.logger.warning(message)

    def error(self, message: str):
        self.logger.error(message)

    def debug(self, message: str):
        self.logger.debug(message)
        
    def process(self, message: str):
        # Use a custom logging level for PROCESS messages
        self.logger.log(25, message)  # 25 is between INFO(20) and WARNING(30)
        # Force flush to ensure real-time output
        for handler in self.logger.handlers:
            if hasattr(handler, 'stream') and hasattr(handler.stream, 'flush'):
                handler.stream.flush()

    def success(self, message: str):
        # Use a custom logging level for SUCCESS messages
        self.logger.log(26, message)  # 26 is between INFO(20) and WARNING(30)

    def test(self, message: str):
        # Use a custom logging level for TEST messages
        self.logger.log(27, message)  # 27 is between INFO(20) and WARNING(30)

    def debug_msg(self, message: str):
        self.logger.debug(message)

    def critical(self, message: str):
        # Use the standard CRITICAL logging level
        self.logger.critical(message)

# Global variable to track the current step logger (set by each step script)
# This is used by print_status to know which logger to use
_current_step_logger = None

def set_step_logger(logger: TEPLogger):
    """Set the current step logger that print_status will use."""
    global _current_step_logger
    _current_step_logger = logger

# Module-level logger instance for backward compatibility
_global_teplogger_instance = None

# Expose print_status and check_memory_usage as global functions
def print_status(message: str, level: str = "INFO") -> None:
    """Prints a status message to the console and logs to the current step logger.
    This function is exposed for consistent use across all scripts.
    """
    # Use the current step logger if set, otherwise just print to console
    if _current_step_logger is not None:
        if level == "SUCCESS":
            _current_step_logger.success(message)
        elif level == "ERROR":
            _current_step_logger.error(message)
        elif level == "WARNING":
            _current_step_logger.warning(message)
        elif level == "PROCESS":
            _current_step_logger.process(message)
        elif level == "DEBUG":
            _current_step_logger.debug(message)
        elif level == "TITLE":
            # For TITLE, log separator and message
            _current_step_logger.info("")
            _current_step_logger.info("="*80)
            _current_step_logger.info(message)
            _current_step_logger.info("="*80)
            _current_step_logger.info("")
        else:
            _current_step_logger.info(message)
    else:
        # No step logger set, just print to console with color
        import datetime
        timestamp = datetime.datetime.now().strftime("%H:%M:%S")
        formatter_instance = TEPFormatter()
        color = formatter_instance.COLORS.get(level, '')
        reset = formatter_instance.RESET
        
        if level == "TITLE":
            print(f"\n{color}{'='*80}{reset}")
            print(f"{color}{message}{reset}")
            print(f"{color}{'='*80}{reset}\n")
        else:
            print(f"{color}[{timestamp}] [{level}] {message}{reset}")
    
    # Force flush stdout to ensure real-time output
    sys.stdout.flush()

def reset_master_log() -> None:
    """Legacy function retained for backward compatibility. No-op."""
    pass

def check_memory_usage(context: str = "Unknown") -> None:
    """
    Logs the current memory usage of the process using the current step logger.
    This function is exposed for consistent use across all scripts.
    """
    import psutil
    import gc

    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    rss_mb = mem_info.rss / (1024 * 1024)
    vms_mb = mem_info.vms / (1024 * 1024)
    
    if _current_step_logger is not None:
        _current_step_logger.debug(f"Memory usage in {context}: RSS={rss_mb:.2f} MB, VMS={vms_mb:.2f} MB")
    else:
        print(f"[DEBUG] Memory usage in {context}: RSS={rss_mb:.2f} MB, VMS={vms_mb:.2f} MB")
    gc.collect()
