#!/usr/bin/env python3
"""
================================================================================
TEP-GNSS-RINEX Analysis Pipeline
STEP 1.0: OPTIMIZED Data Acquisition & Processing (ALL STATIONS)
================================================================================

STRATEGY:
  - Download ALL available RINEX stations (~300/day), not just 100
  - Filter to optimal stations in ANALYSIS steps (Step 2.x)
  - This provides FLEXIBILITY:
    * Station availability varies by year (2020 vs 2024)
    * Can dynamically select best stations per time period
    * Can try different selection criteria without re-downloading
    * More robust global coverage

OPTIMIZATIONS:
  - Batch downloads using CDDIS wildcards (1 request = 1 TAR with ~300 files)
  - ~100x fewer HTTP requests than individual downloads
  - Multiprocessing for CPU-bound RTKLIB processing
  - Parallel RTKLIB on all CPU cores

ESTIMATED OUTPUT:
  - ~300 stations × 1827 days = ~550,000 station-day files
  - ~50KB per file = ~27 GB total processed data

Author: Matthew Lukin Smawfield
Date: 9 December 2025
License: CC-BY-4.0
"""

import os
import sys
import gzip
import io
import subprocess
import json
import tarfile
import requests
import numpy as np
from datetime import datetime, timedelta
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from multiprocessing import cpu_count
from tqdm import tqdm
import tempfile
import time as _time  # Avoid conflict with time module
import shutil
import hatanaka
import unlzw3

# Force unbuffered output for real-time logging
os.environ['PYTHONUNBUFFERED'] = '1'

# Setup path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.logger import TEPLogger, set_step_logger

# ============================================================================
# CONFIGURATION
# ============================================================================

# Note: We download ALL stations, filtering happens in analysis steps
# This flag is kept for backwards compatibility but doesn't affect downloads
USE_OPTIMAL_100 = True  # Only affects logging, not downloads

PROJECT_ROOT = Path(__file__).parent.parent.parent
NAV_DIR = PROJECT_ROOT / "data" / "nav"
SP3_DIR = PROJECT_ROOT / "data" / "sp3"
RINEX_DIR = PROJECT_ROOT / "data" / "rinex"  # Temp storage for batch downloads
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
LOGS_DIR = PROJECT_ROOT / "logs"
OUTPUTS_DIR = PROJECT_ROOT / "results" / "outputs"

# Ensure directories exist
for d in [NAV_DIR, SP3_DIR, RINEX_DIR, PROCESSED_DIR, LOGS_DIR, OUTPUTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# RTKLIB path (configurable via env or CLI)
DEFAULT_RTKLIB = PROJECT_ROOT / "RTKLIB" / "app" / "consapp" / "rnx2rtkp" / "gcc" / "rnx2rtkp"

# Initialize logger
logger = TEPLogger("step_1_0", log_file_path=LOGS_DIR / "step_1_0.log")
set_step_logger(logger)

COORDS_FILE = PROCESSED_DIR / "station_coordinates.json"
CDDIS_BASE = "https://cddis.nasa.gov/archive/gnss/data/daily"

START_DATE = datetime(2022, 1, 1)  # 3 years: 2022-2024 (Galileo mature)
END_DATE = datetime(2024, 12, 31)

OPTIMAL_METADATA_FILE = PROCESSED_DIR / "optimal_100_metadata.json"


def get_auth():
    """Get CDDIS authentication credentials from .netrc."""
    import netrc
    try:
        auth = netrc.netrc().authenticators('urs.earthdata.nasa.gov')
        if auth:
            return (auth[0], auth[2])
    except Exception:
        pass
    user = os.getenv('CDDIS_USER')
    passwd = os.getenv('CDDIS_PASS')
    if user and passwd:
        return (user, passwd)
    return None


def load_coordinates():
    if COORDS_FILE.exists():
        with open(COORDS_FILE, 'r') as f:
            return json.load(f)
    return {}


def save_coordinates(coords):
    with open(COORDS_FILE, 'w') as f:
        json.dump(coords, f, indent=2)


def load_optimal_stations():
    """Load the curated 100 optimal stations."""
    if not OPTIMAL_METADATA_FILE.exists():
        logger.error(f"Optimal stations metadata not found: {OPTIMAL_METADATA_FILE}")
        return None
    
    with open(OPTIMAL_METADATA_FILE) as f:
        metadata = json.load(f)
    
    stations = []
    for sta4, data in metadata['stations'].items():
        igs_info = data.get('igs_info', {})
        code9 = igs_info.get('code_9char', f"{sta4}00XXX")
        stations.append(code9)
    
    logger.info(f"  Loaded {len(stations)} optimal stations")
    return stations


# ============================================================================
# BATCH DOWNLOAD FUNCTIONS (Using CDDIS wildcards)
# ============================================================================

def download_multi_gnss_nav(year, doy, auth):
    """
    Download multi-GNSS combined navigation file for a specific day.
    
    Tries (in order):
    1. BRDC00IGS combined nav (GPS+GLO+GAL+BDS+QZSS) - best option
    2. Individual constellation nav files as fallback
    
    Returns: Path to nav file or None
    """
    yy = str(year)[-2:]
    session = requests.Session()
    session.auth = auth
    
    # Output file
    nav_file = NAV_DIR / f"BRDC_{year}{doy:03d}.nav"
    if nav_file.exists():
        return nav_file
    
    # Try 1: IGS Combined Multi-GNSS Broadcast (RINEX3)
    # Format: BRDC00IGS_R_YYYYDOY0000_01D_MN.rnx.gz
    combined_urls = [
        f"https://cddis.nasa.gov/archive/gnss/data/daily/{year}/{doy:03d}/{yy}p/BRDC00IGS_R_{year}{doy:03d}0000_01D_MN.rnx.gz",
        f"https://cddis.nasa.gov/archive/gnss/data/daily/{year}/brdc/BRDC00IGS_R_{year}{doy:03d}0000_01D_MN.rnx.gz",
    ]
    
    for url in combined_urls:
        try:
            resp = session.get(url, timeout=60)
            if resp.status_code == 200:
                content = gzip.decompress(resp.content)
                with open(nav_file, 'wb') as f:
                    f.write(content)
                return nav_file
        except:
            continue
    
    # Try 2: GPS-only broadcast nav (fallback)
    gps_urls = [
        f"{CDDIS_BASE}/{year}/brdc/brdc{doy:03d}0.{yy}n.Z",
        f"{CDDIS_BASE}/{year}/brdc/brdc{doy:03d}0.{yy}n.gz",
        f"{CDDIS_BASE}/{year}/{doy:03d}/{yy}n/brdc{doy:03d}0.{yy}n.gz",
    ]
    
    for url in gps_urls:
        try:
            resp = session.get(url, timeout=30)
            if resp.status_code == 200:
                content = resp.content
                if url.endswith('.Z'):
                    content = unlzw3.unlzw(content)
                elif url.endswith('.gz'):
                    content = gzip.decompress(content)
                with open(nav_file, 'wb') as f:
                    f.write(content)
                return nav_file
        except:
            continue
    
    return None


def batch_download_nav_files(year, auth):
    """
    Download navigation files for a year.
    Tries multi-GNSS combined nav first, falls back to GPS-only.
    """
    yy = str(year)[-2:]
    
    # Check existing files
    existing_combined = list(NAV_DIR.glob(f"BRDC_{year}*.nav"))
    existing_gps = list(NAV_DIR.glob(f"brdc*.{yy}n"))
    total_existing = len(existing_combined) + len(existing_gps)
    n_days = 366 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 365
    
    if total_existing >= n_days:
        print(f"[NAV] {year}: {total_existing}/{n_days} files exist ✓ (skipping)", flush=True)
        return total_existing
    
    print(f"[NAV] {year}: Downloading multi-GNSS navigation files ({n_days} days)...", flush=True)
    
    # Download individually with multi-GNSS support
    count = 0
    multi_gnss = 0
    gps_only = 0
    failed = 0
    
    for doy in range(1, n_days + 1):
        nav_file = download_multi_gnss_nav(year, doy, auth)
        if nav_file:
            count += 1
            if 'BRDC_' in str(nav_file):
                multi_gnss += 1
            else:
                gps_only += 1
        else:
            failed += 1
        
        # Progress every 100 days
        if doy % 100 == 0:
            print(f"[NAV] {year}: {doy}/{n_days} days | {multi_gnss} multi-GNSS + {gps_only} GPS-only | {failed} failed", flush=True)
    
    print(f"[NAV] {year}: Complete! {count}/{n_days} files ({multi_gnss} multi-GNSS, {gps_only} GPS-only, {failed} failed)", flush=True)
    return count


def download_nav_individual(year, auth):
    """Fallback: download nav files one by one."""
    yy = str(year)[-2:]
    count = 0
    
    # Determine days in year
    if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0):
        n_days = 366
    else:
        n_days = 365
    
    session = requests.Session()
    session.auth = auth
    
    for doy in range(1, n_days + 1):
        nav_file = NAV_DIR / f"brdc{doy:03d}0.{yy}n"
        if nav_file.exists():
            count += 1
            continue
        
        urls = [
            f"{CDDIS_BASE}/{year}/brdc/brdc{doy:03d}0.{yy}n.Z",
            f"{CDDIS_BASE}/{year}/brdc/brdc{doy:03d}0.{yy}n.gz",
            f"{CDDIS_BASE}/{year}/{doy:03d}/{yy}n/brdc{doy:03d}0.{yy}n.gz",
        ]
        
        for url in urls:
            try:
                resp = session.get(url, timeout=30)
                if resp.status_code == 200:
                    content = resp.content
                    if url.endswith('.Z'):
                        content = unlzw3.unlzw(content)
                    elif url.endswith('.gz'):
                        content = gzip.decompress(content)
                    with open(nav_file, 'wb') as f:
                        f.write(content)
                    count += 1
                    break
            except:
                continue
    
    return count


def download_sp3_files(years, auth):
    """
    Download SP3 precise orbit files for all DOYs in the given years.
    Uses CODE final products (COD0MGXFIN) which are best for multi-GNSS.
    
    Returns: dict of (year, doy) -> sp3_path
    """
    import gzip
    
    session = requests.Session()
    session.auth = auth
    
    sp3_files = {}
    total_downloaded = 0
    total_skipped = 0
    
    for year in years:
        n_days = 366 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 365
        
        for doy in range(1, n_days + 1):
            sp3_file = SP3_DIR / f"sp3_{year}{doy:03d}.sp3"
            
            # Skip if already exists
            if sp3_file.exists():
                sp3_files[(year, doy)] = sp3_file
                total_skipped += 1
                continue
            
            # Convert year/doy to GPS week/day
            from datetime import datetime
            date = datetime(year, 1, 1) + timedelta(days=doy - 1)
            gps_epoch = datetime(1980, 1, 6)
            days_since_epoch = (date - gps_epoch).days
            gps_week = days_since_epoch // 7
            dow = days_since_epoch % 7
            
            # Try different SP3 products (CODE final is best for multi-GNSS)
            urls = [
                # CODE MGEX final
                f"https://cddis.nasa.gov/archive/gnss/products/{gps_week}/COD0MGXFIN_{year}{doy:03d}0000_01D_05M_ORB.SP3.gz",
                # CODE rapid
                f"https://cddis.nasa.gov/archive/gnss/products/{gps_week}/COD0OPSFIN_{year}{doy:03d}0000_01D_05M_ORB.SP3.gz",
                # IGS final
                f"https://cddis.nasa.gov/archive/gnss/products/{gps_week}/igs{year:02d}P{gps_week}{dow}.sp3.Z",
                # IGS rapid
                f"https://cddis.nasa.gov/archive/gnss/products/{gps_week}/igr{year:02d}P{gps_week}{dow}.sp3.Z",
            ]
            
            downloaded = False
            for url in urls:
                try:
                    resp = session.get(url, timeout=60)
                    if resp.status_code == 200:
                        content = resp.content
                        # Decompress
                        if url.endswith('.gz'):
                            content = gzip.decompress(content)
                        elif url.endswith('.Z'):
                            content = unlzw3.unlzw(content)
                        
                        with open(sp3_file, 'wb') as f:
                            f.write(content)
                        
                        sp3_files[(year, doy)] = sp3_file
                        total_downloaded += 1
                        downloaded = True
                        break
                except Exception as e:
                    continue
            
            if not downloaded:
                # SP3 not available for this day (normal for recent days)
                pass
        
        print(f"[SP3] {year}: {total_downloaded} downloaded, {total_skipped} existed", flush=True)
    
    return sp3_files


def get_available_files_for_doy(year, doy, session):
    """
    Get list of available RINEX files for a DOY using CDDIS directory listing.
    Uses *?list to get file listing without downloading.
    
    Returns: list of filenames available, or None on error
    """
    yy = str(year)[-2:]
    list_url = f"{CDDIS_BASE}/{year}/{doy:03d}/{yy}d/*_R_{year}{doy:03d}0000_01D_30S_MO.crx.gz?list"
    
    try:
        resp = session.get(list_url, timeout=30)
        if resp.status_code == 200:
            # Parse file listing - one file per line
            files = [f.strip() for f in resp.text.strip().split('\n') if f.strip()]
            return files
        else:
            return None
    except Exception as e:
        return None


def batch_download_rinex_for_doy(year, doy, auth):
    """
    Download ALL RINEX files for a specific DOY using Python requests.
    
    Uses CDDIS wildcard URL which returns a TAR file.
    Per CDDIS docs: requests.get(url) with .netrc authentication.
    
    Returns: list of rinex_paths (ALL stations)
    """
    start_time = _time.time()
    yy = str(year)[-2:]
    doy_str = f"{year}/{doy:03d}"
    doy_dir = RINEX_DIR / f"{year}" / f"{doy:03d}"
    doy_dir.mkdir(parents=True, exist_ok=True)
    
    # Use unique tar filename with PID and timestamp to avoid conflicts
    tar_file = doy_dir / f"rinex_{year}{doy:03d}_{os.getpid()}_{int(_time.time()*1000) % 100000}.tar"
    
    # Create session with auth (per CDDIS Python docs)
    session = requests.Session()
    session.auth = auth
    
    def log(status, msg=""):
        """Consistent log format"""
        timestamp = datetime.now().strftime("%H:%M:%S")
        if msg:
            print(f"[{timestamp}] [{doy_str}] {status}: {msg}", flush=True)
        else:
            print(f"[{timestamp}] [{doy_str}] {status}", flush=True)

    # -------------------------------------------------------------------------
    # CHECK IF ALREADY DOWNLOADED (check RINEX files in doy_dir)
    # We ALWAYS need RINEX files to process - check this first
    # -------------------------------------------------------------------------
    existing_rinex = list(doy_dir.glob("*.crx.gz")) + list(doy_dir.glob("*.rnx"))
    
    # -------------------------------------------------------------------------
    # CHECK IF ALREADY FULLY PROCESSED WITH ENHANCED MODES
    # PRIORITY: Check NPZ first - if complete, skip entirely (even if RINEX deleted)
    # -------------------------------------------------------------------------
    processed_dir = Path(__file__).parent.parent.parent / "data" / "processed"
    existing_npz = list(processed_dir.glob(f"*_{year}{doy:03d}.npz"))
    
    # PRIMARY CHECK: If we have 280+ NPZ with ALL modes including precise, we're DONE
    # No need to check RINEX - processing is complete
    if len(existing_npz) >= 280:
        try:
            sample_file = existing_npz[0]
            with np.load(sample_file) as data:
                keys = set(data.keys())
                # CRITICAL: Must have precise_clock_bias_ns to skip
                has_precise = 'precise_clock_bias_ns' in keys
                if has_precise:
                    log("SKIP", f"Found {len(existing_npz)} NPZ with precise mode (fully done)")
                    return []  # Fully processed - skip entirely
        except:
            pass
    
    # SECONDARY: Have RINEX but not enough/enhanced NPZ - return RINEX for processing
    if len(existing_rinex) >= 280:
        if len(existing_npz) >= 280:
            log("READY", f"Found {len(existing_rinex)} RINEX, {len(existing_npz)} NPZ need enhancement")
        else:
            log("READY", f"Found {len(existing_rinex)} RINEX, only {len(existing_npz)} NPZ - processing needed")
        return existing_rinex
    
    # -------------------------------------------------------------------------
    # STEP 1: Query CDDIS directory listing
    # -------------------------------------------------------------------------
    try:
        list_start = _time.time()
        expected_files = get_available_files_for_doy(year, doy, session)
        expected_count = len(expected_files) if expected_files else 0
        list_time = _time.time() - list_start
        
        if expected_count == 0:
            log("SKIP", f"No RINEX files on CDDIS")
            return []
        
        # Estimate size: ~2.5MB per file average
        est_mb = expected_count * 2.5
        log("START", f"Found {expected_count} files (~{est_mb:.0f}MB), downloading...")
        
    except Exception as e:
        log("ERROR", f"Failed to list directory: {e}")
        return []
    
    # -------------------------------------------------------------------------
    # STEP 2: Download TAR archive from CDDIS
    # -------------------------------------------------------------------------
    url = f"{CDDIS_BASE}/{year}/{doy:03d}/{yy}d/*_R_{year}{doy:03d}0000_01D_30S_MO.crx.gz"
    
    try:
        download_start = _time.time()
        resp = session.get(url, timeout=600, stream=True)
        
        if resp.status_code != 200:
            log("FAILED", f"HTTP {resp.status_code} from CDDIS")
            return []
        
        # Get expected size
        total_size = int(resp.headers.get('content-length', 0))
        total_mb = total_size / 1024 / 1024 if total_size else 0
        
        # Download with progress (less frequent - every 30s)
        downloaded = 0
        last_progress = _time.time()
        
        with open(tar_file, 'wb') as f:
            for chunk in resp.iter_content(chunk_size=1024*1024):
                f.write(chunk)
                downloaded += len(chunk)
                
                # Progress update every 30 seconds
                now = _time.time()
                if now - last_progress >= 30:
                    elapsed = now - download_start
                    speed = downloaded / elapsed / 1024 / 1024
                    log("DOWNLOADING", f"{downloaded/1024/1024:.0f}MB downloaded @ {speed:.1f}MB/s")
                    last_progress = now
        
        download_time = _time.time() - download_start
        
        # Verify download
        if not tar_file.exists():
            log("FAILED", "TAR file not created")
            return []
        
        tar_size = tar_file.stat().st_size
        tar_mb = tar_size / 1024 / 1024
        speed = tar_size / download_time / 1024 / 1024
        
        if tar_size < 1000:
            log("FAILED", f"TAR too small ({tar_size} bytes)")
            try:
                tar_file.unlink()
            except:
                pass
            return []
        
        log("DOWNLOADED", f"{tar_mb:.1f}MB in {download_time:.0f}s ({speed:.1f}MB/s)")
        
    except requests.exceptions.Timeout:
        log("FAILED", f"Download timeout after {_time.time() - start_time:.0f}s")
        try:
            tar_file.unlink()
        except:
            pass
        return []
    except Exception as e:
        log("FAILED", f"Download error: {type(e).__name__}: {e}")
        try:
            tar_file.unlink()
        except:
            pass
        return []
    
    # -------------------------------------------------------------------------
    # STEP 3: Extract RINEX files from TAR
    # -------------------------------------------------------------------------
    try:
        extract_start = _time.time()
        rinex_files = []
        
        # Verify TAR magic bytes
        with open(tar_file, 'rb') as f:
            f.seek(257)
            magic = f.read(5)
            if magic != b'ustar':
                log("FAILED", "Invalid TAR file (not ustar format)")
                tar_file.unlink()
                return []
        
        # Extract files
        with tarfile.open(tar_file, 'r') as tar:
            members = [m for m in tar.getmembers() if m.name.endswith('.crx.gz') and '_R_' in m.name]
            for member in members:
                fname = Path(member.name).name
                f = tar.extractfile(member)
                if f:
                    out_path = doy_dir / fname
                    with open(out_path, 'wb') as out:
                        out.write(f.read())
                    rinex_files.append(out_path)
        
        extract_time = _time.time() - extract_start
        
        # Clean up TAR
        try:
            tar_file.unlink()
        except:
            pass
        
        got_count = len(rinex_files)
        total_time = _time.time() - start_time
        
        # Final status - only warn if significant mismatch (>5% missing)
        if got_count == 0:
            log("FAILED", "No files extracted from TAR")
        elif expected_count > 0 and got_count < expected_count * 0.95:
            # Significant mismatch - more than 5% missing
            log("DONE", f"Extracted {got_count}/{expected_count} files ({total_time:.0f}s) ⚠️ {expected_count - got_count} MISSING")
        else:
            log("DONE", f"Extracted {got_count} files ({total_time:.0f}s) ✓")
        
        return rinex_files
        
    except Exception as e:
        log("FAILED", f"Extraction error: {type(e).__name__}: {e}")
        try:
            tar_file.unlink()
        except:
            pass
        return []


# ============================================================================
# PROCESSING MODES - Comprehensive multi-mode analysis
# ============================================================================
# Each mode provides different error characteristics for TEP validation:
# - baseline: Standard SPP - all errors present (ionosphere, orbit, clock)
# - precise: Precise orbits - removes orbit errors, isolates iono+clock
# - ionofree: Dual-freq - removes 99% ionosphere, isolates orbit+clock
# - Constellation-specific: Different clock types, orbital altitudes, systematics
#
# SCIENTIFIC VALUE OF CONSTELLATION COMPARISON:
# - GPS (G): MEO ~20,200 km, Rb/Cs clocks, oldest/most stable
# - GLONASS (R): MEO ~19,100 km, Cs clocks, FDMA (different signal structure)
# - Galileo (E): MEO ~23,222 km, H-maser clocks (most precise!)
# - BeiDou (C): MEO/GEO/IGSO mix, Rb clocks, newest constellation
#
# If TEP is real: ALL constellations should show similar λ and R²
# If TEP is artifact: Different systematics would produce different signatures

PROCESSING_MODES = {
    # === PRIMARY: GPS baseline (reference) ===
    # -ti 300 = 5-minute intervals (consistent with CODE longspan analysis)
    'baseline': {
        'flags': ['-p', '0', '-sys', 'G', '-ti', '300'],
        'needs_sp3': False,
        'prefix': '',  # Primary mode, no prefix
        'description': 'GPS SPP with broadcast ephemeris (5-min)',
    },
    
    # === ENHANCED PROCESSING: SPP with Iono-Free LC ===
    # SPP mode (-p 0) with ionosphere correction = Iono-Free LC
    # Per RTKLIB manual E.6: "If Ionosphere Correction set to Iono-Free LC,
    # the ionosphere-free LC pseudorange is used" - eliminates 1st order ionosphere!
    # Requires config file with pos1-ionoopt=iflc
    'ionofree': {
        'flags': ['-p', '0', '-sys', 'G', '-ti', '300'],
        'needs_sp3': False,
        'needs_config': True,  # Needs config file for iono-free LC
        'prefix': 'ionofree_',
        'description': 'GPS SPP dual-freq iono-free LC (5-min)',
    },
    
    # === GALILEO: H-maser clocks (most precise) ===
    'galileo': {
        'flags': ['-p', '0', '-sys', 'E', '-ti', '300'],
        'needs_sp3': False,
        'prefix': 'gal_',
        'description': 'Galileo-only SPP H-maser (5-min)',
    },
    
    # === COMBINED MULTI-GNSS ===
    'multi_gnss': {
        'flags': ['-p', '0', '-sys', 'G,R,E,C', '-ti', '300'],
        'needs_sp3': False,
        'prefix': 'mgex_',
        'description': 'All constellations combined (5-min)',
    },
    
    # === PRECISE: GPS SPP with Precise Orbits/Clocks ===
    # CRITICAL FOR TEP: Removes N-S orbital systematics
    'precise': {
        'flags': ['-p', '0', '-sys', 'G', '-ti', '300'],
        'needs_sp3': True,
        'needs_config': True,
        'prefix': 'precise_',
        'description': 'GPS SPP with IGS Precise Ephemeris',
    },
}


# ============================================================================
# RTKLIB PROCESSING (CPU-bound, uses multiprocessing)
# ============================================================================

def run_rtklib_mode(rtklib_path, rinex_file, nav_file, sp3_file, tmpdir, mode_name, mode_config):
    """
    Run RTKLIB for a single processing mode.
    Returns dict with clock_bias_ns, pos_jitter, clock_drift, or None if failed.
    """
    pos_file = tmpdir / f"output_{mode_name}.pos"
    stat_file = tmpdir / f"output_{mode_name}.pos.stat"
    
    # Build command
    cmd = [
        str(rtklib_path),
        '-e',       # ECEF output
        '-t',       # Time format
        '-y', '1',  # Solution status
        '-o', str(pos_file),
    ]
    
    # Create config file if needed
    if mode_config.get('needs_config'):
        config_file = tmpdir / f"config_{mode_name}.conf"
        
        if mode_name == 'ionofree':
            # Iono-Free LC Config
            config_content = """# RTKLIB config for Iono-Free LC SPP
pos1-posmode       =single
pos1-frequency     =l1+l2
pos1-soltype       =forward
pos1-elmask        =15
pos1-dynamics      =off
pos1-tidecorr      =off
pos1-ionoopt       =dual-freq
pos1-tropopt       =saas
pos1-sateph        =brdc
pos1-navsys        =1
"""
        elif mode_name == 'precise':
            # Precise Ephemeris Config (Dual-Freq Iono-Free + Precise Orbits)
            config_content = """# RTKLIB config for Precise Iono-Free SPP
pos1-posmode       =single
pos1-frequency     =l1+l2
pos1-soltype       =forward
pos1-elmask        =15
pos1-dynamics      =off
pos1-tidecorr      =off
pos1-ionoopt       =dual-freq
pos1-tropopt       =saas
pos1-sateph        =precise
pos1-navsys        =1
"""
        else:
            # Default/Fallback
            config_content = ""
            
        config_file.write_text(config_content)
        cmd.extend(['-k', str(config_file)])
    
    cmd.extend(mode_config['flags'])
    cmd.append(str(rinex_file))
    
    # Add SP3 if needed
    if mode_config.get('needs_sp3') and sp3_file and sp3_file.exists():
        cmd.append(str(sp3_file))
    elif mode_config.get('needs_sp3'):
        return None  # Skip mode if SP3 required but not available
    
    cmd.append(str(nav_file))
    
    try:
        subprocess.run(cmd, capture_output=True, timeout=120)
    except:
        return None
    
    if not pos_file.exists():
        return None
    
    # Parse position file
    timestamps, x_data, y_data, z_data = [], [], [], []
    nsat_data, q_data = [], []
    
    try:
        with open(pos_file, 'r') as f:
            for line in f:
                if line.startswith('%') or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) < 8:
                    continue
                try:
                    time_parts = parts[1].split(':')
                    sod = float(time_parts[0]) * 3600 + float(time_parts[1]) * 60 + float(time_parts[2])
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    q = int(parts[5]) if len(parts) > 5 else 5
                    ns = int(parts[6]) if len(parts) > 6 else 0
                    
                    if x == 0 and y == 0 and z == 0:
                        continue
                    
                    timestamps.append(sod)
                    x_data.append(x)
                    y_data.append(y)
                    z_data.append(z)
                    nsat_data.append(ns)
                    q_data.append(q)
                except:
                    continue
        
        if len(x_data) < 100:
            return None
        
        # Parse clock from stat file
        clock_data = []
        if stat_file.exists():
            with open(stat_file, 'r') as f:
                for line in f:
                    if line.startswith('$CLK'):
                        try:
                            parts = line.split(',')
                            if len(parts) >= 6:
                                clock_data.append(float(parts[5]))
                        except:
                            continue
        
        # Need actual clock data for meaningful results
        if not clock_data or len(clock_data) < 100:
            return None
            
        min_len = min(len(x_data), len(clock_data))
        if min_len < 100:
            return None
        
        # Reject if clock data is all zeros (failed processing)
        if all(c == 0 for c in clock_data[:min_len]):
            return None
        
        # Create arrays
        ts_arr = np.array(timestamps[:min_len], dtype=np.float32)
        x_arr = np.array(x_data[:min_len], dtype=np.float32)
        y_arr = np.array(y_data[:min_len], dtype=np.float32)
        z_arr = np.array(z_data[:min_len], dtype=np.float32)
        clk_arr = np.array(clock_data[:min_len], dtype=np.float32)
        nsat_arr = np.array(nsat_data[:min_len], dtype=np.uint8)
        q_arr = np.array(q_data[:min_len], dtype=np.uint8)
        
        # Compute metrics
        x_mean, y_mean, z_mean = np.mean(x_arr), np.mean(y_arr), np.mean(z_arr)
        dx = (x_arr - x_mean).astype(np.float32)
        dy = (y_arr - y_mean).astype(np.float32)
        dz = (z_arr - z_mean).astype(np.float32)
        dr = np.sqrt(dx**2 + dy**2 + dz**2).astype(np.float32)
        clk_ns = (clk_arr * 3.33564095).astype(np.float32)
        
        dt = np.diff(ts_arr)
        dt[dt == 0] = 30.0
        clk_drift = np.diff(clk_ns) / dt
        clk_drift = np.append(clk_drift, 0).astype(np.float32)
        
        return {
            'clock_bias_ns': clk_ns,
            'pos_jitter': dr,
            'clock_drift': clk_drift,
            'timestamps': ts_arr,
            'nsat': nsat_arr,
            'quality': q_arr,
            'dx': dx, 'dy': dy, 'dz': dz,
            'x_mean': x_mean, 'y_mean': y_mean, 'z_mean': z_mean,
            'n_epochs': len(ts_arr)
        }
        
    except:
        return None


def process_rinex_file(args):
    """
    Process a single RINEX file with RTKLIB - ALL processing modes.
    This is CPU-bound so we use multiprocessing.
    """
    rinex_gz_path, nav_file, sp3_file, rtklib_path, output_dir = args
    
    try:
        # Get station and date from filename
        fname = rinex_gz_path.name
        station = fname[:4].upper()
        
        # Parse date from filename: STATION_R_YYYYDOY0000_...
        parts = fname.split('_')
        if len(parts) >= 3:
            date_part = parts[2]  # YYYYDOY0000
            year = int(date_part[:4])
            doy = int(date_part[4:7])
        else:
            return None
        
        key = f"{station}_{year}{doy:03d}"
        output_file = output_dir / f"{key}.npz"
        
        # Skip if already processed WITH all 5 modes
        # Modes: baseline, ionofree, galileo, multi_gnss, precise
        if output_file.exists():
            try:
                with np.load(output_file) as data:
                    keys = set(data.keys())
                    has_baseline = 'clock_bias_ns' in keys
                    has_ionofree = 'ionofree_clock_bias_ns' in keys
                    has_galileo = 'gal_clock_bias_ns' in keys
                    has_multi = 'mgex_clock_bias_ns' in keys
                    has_precise = 'precise_clock_bias_ns' in keys
                    
                    if has_baseline and has_ionofree and has_galileo and has_multi and has_precise:
                        return {"key": key, "status": "exists"}
                    # File exists but missing modes - will reprocess
            except:
                pass
        
        # Verify file exists before processing
        if not rinex_gz_path.exists():
            return {"key": key, "status": "fail", "error": f"file_not_found: {rinex_gz_path.name}"}
        
        # Check file is not empty
        if rinex_gz_path.stat().st_size < 100:
            return {"key": key, "status": "fail", "error": "empty_file"}
        
        # Work in temp directory
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            
            # Decompress RINEX
            with gzip.open(rinex_gz_path, 'rb') as gz:
                content = gz.read()
            
            # Hatanaka decompress if needed
            if fname.endswith('.crx.gz'):
                content = hatanaka.decompress(content)
                rinex_file = tmpdir / fname.replace('.crx.gz', '.rnx')
            else:
                rinex_file = tmpdir / fname.replace('.gz', '')
            
            with open(rinex_file, 'wb') as f:
                f.write(content)
            
            # Extract coordinates
            coords = None
            with open(rinex_file, 'r', errors='ignore') as f:
                for line in f:
                    if 'APPROX POSITION XYZ' in line:
                        parts = line.split()
                        coords = [float(parts[0]), float(parts[1]), float(parts[2])]
                        break
                    if 'END OF HEADER' in line:
                        break
            
            if not coords:
                return {"key": key, "status": "no_coords"}
            
            # Run ALL processing modes
            output_data = {
                'station': station,
                'date': f"{year}-{doy:03d}",
                'doy': doy,
                'sampling_sec': 300,  # 5-minute intervals
                'approx_x': coords[0],
                'approx_y': coords[1],
                'approx_z': coords[2],
            }
            
            baseline_success = False
            modes_success = []
            
            for mode_name, mode_config in PROCESSING_MODES.items():
                result = run_rtklib_mode(
                    rtklib_path, rinex_file, nav_file, sp3_file, 
                    tmpdir, mode_name, mode_config
                )
                
                if result:
                    prefix = mode_config['prefix']
                    
                    # Add all keys with prefix
                    output_data[f'{prefix}clock_bias_ns'] = result['clock_bias_ns']
                    output_data[f'{prefix}pos_jitter'] = result['pos_jitter']
                    output_data[f'{prefix}clock_drift'] = result['clock_drift']
                    
                    # For baseline (no prefix), also add timestamps and other common data
                    if mode_name == 'baseline':
                        output_data['timestamps'] = result['timestamps']
                        output_data['nsat'] = result['nsat']
                        output_data['quality'] = result['quality']
                        output_data['dx'] = result['dx']
                        output_data['dy'] = result['dy']
                        output_data['dz'] = result['dz']
                        output_data['x_mean'] = result['x_mean']
                        output_data['y_mean'] = result['y_mean']
                        output_data['z_mean'] = result['z_mean']
                        output_data['n_epochs'] = result['n_epochs']
                        baseline_success = True
                    
                    modes_success.append(mode_name)
            
            # Must have at least baseline
            if not baseline_success:
                return {"key": key, "status": "rtklib_failed"}
            
            # Save with all modes
            np.savez_compressed(output_file, **output_data)
            
            return {
                "key": key,
                "status": "success",
                "station": station,
                "coords": coords,
                "n_epochs": output_data.get('n_epochs', 0),
                "modes": modes_success
            }
            
    except Exception as e:
        return {"key": str(rinex_gz_path), "status": f"error: {e}"}


# ============================================================================
# MAIN
# ============================================================================

def main():
    logger.info("")
    logger.info("=" * 80)
    logger.info("TEP-GNSS-RINEX - STEP 1.0: OPTIMIZED Batch Acquisition")
    logger.info("=" * 80)
    logger.info("")
    
    auth = get_auth()
    if not auth:
        logger.error("No CDDIS auth! Create ~/.netrc with urs.earthdata.nasa.gov credentials")
        return
    
    # We still load optimal stations for metadata, but download ALL stations
    _ = load_optimal_stations()  # Just for logging
    
    logger.info(f"  Auth user: {auth[0]}")
    logger.info(f"  Stations: ALL available (~300 per day)")
    logger.info(f"  Date range: {START_DATE.date()} to {END_DATE.date()}")
    
    years = list(range(START_DATE.year, END_DATE.year + 1))
    logger.info(f"  Years: {years}")
    
    if not DEFAULT_RTKLIB.exists():
        logger.error(f"RTKLIB not found: {DEFAULT_RTKLIB}")
        return
    
    all_coords = load_coordinates()
    
    # ========================================================================
    # PHASE 1: Batch download navigation files (one request per year)
    # ========================================================================
    logger.process("[1/3] Downloading navigation files (batch mode)...")
    
    nav_files = {}
    multi_gnss_count = 0
    gps_only_count = 0
    
    for year in years:
        count = batch_download_nav_files(year, auth)
        yy = str(year)[-2:]
        
        # Check for multi-GNSS nav files (BRDC_YEARDOY.nav)
        for nav_path in NAV_DIR.glob(f"BRDC_{year}*.nav"):
            # Parse DOY from filename: BRDC_2020001.nav -> doy=1
            try:
                doy = int(nav_path.stem.split('_')[1][4:7])
                nav_files[(year, doy)] = nav_path
                multi_gnss_count += 1
            except:
                continue
        
        # Also check for GPS-only nav files (brdc0010.20n) as fallback
        for nav_path in NAV_DIR.glob(f"brdc*.{yy}n"):
            # Parse DOY from filename: brdc0010.20n -> doy=1
            doy_str = nav_path.stem[4:7]
            try:
                doy = int(doy_str)
                if (year, doy) not in nav_files:  # Don't overwrite multi-GNSS
                    nav_files[(year, doy)] = nav_path
                    gps_only_count += 1
            except:
                continue
    
    logger.info(f"  Navigation files ready: {len(nav_files)} ({multi_gnss_count} multi-GNSS, {gps_only_count} GPS-only)")
    
    # ========================================================================
    # PHASE 2: Download SP3 precise orbit files (for ionofree mode)
    # ========================================================================
    logger.process("[2/4] Downloading SP3 precise orbit files...")
    
    sp3_files = download_sp3_files(years, auth)
    logger.info(f"  SP3 files ready: {len(sp3_files)}")
    
    # ========================================================================
    # PHASE 3: Batch download + process RINEX files (one request per DOY)
    # ========================================================================
    logger.process("[3/4] Downloading & processing RINEX files (batch mode)...")
    
    # Build list of ALL DOYs to process (don't skip - days may be incomplete)
    doys_to_process = []
    current = START_DATE
    while current <= END_DATE:
        year = current.year
        doy = current.timetuple().tm_yday
        if (year, doy) in nav_files:
            doys_to_process.append((year, doy))
        current += timedelta(days=1)
    
    # Count existing processed files and estimate DOYs completed
    existing_files = len(list(PROCESSED_DIR.glob("*.npz")))
    doys_completed = existing_files // 300  # ~300 stations per day
    doys_remaining = len(doys_to_process) - doys_completed
    expected_total = len(doys_to_process) * 300
    
    print(f"\n{'='*70}", flush=True)
    print(f"PROGRESS SUMMARY", flush=True)
    print(f"{'='*70}", flush=True)
    print(f"  Total DOYs:       {len(doys_to_process):,}", flush=True)
    print(f"  Already done:     ~{doys_completed:,} DOYs ({existing_files:,} files)", flush=True)
    print(f"  Remaining:        ~{doys_remaining:,} DOYs", flush=True)
    print(f"  Est. time:        ~{doys_remaining * 2 / 60:.1f} hours @ 2min/DOY", flush=True)
    print(f"{'='*70}\n", flush=True)
    
    # ========================================================================
    # PARALLEL PIPELINE: Download + Process simultaneously
    # ========================================================================
    # - 8 parallel download threads (network-bound, using requests)
    # - All CPUs for RTKLIB processing (CPU-bound)
    # - Producer-consumer pattern with queue
    
    from queue import Queue
    from threading import Thread
    import shutil
    
    # OPTIMIZED config for 7-mode processing (CPU-bound)
    # With 4 RTKLIB modes per file, processing is the bottleneck
    # MACHINE: e2-highcpu-32 (32 vCPUs)
    # Prioritize fast DOY completion over parallelism
    N_DOWNLOAD_WORKERS = 2   # 2 parallel downloads (processing is bottleneck anyway)
    N_DOY_PROCESSORS = 2     # 2 DOYs in parallel (faster per-DOY completion)
    WORKERS_PER_DOY = 16     # 16 RTKLIB workers per DOY (2×16=32 total)
    
    logger.info(f"  Download threads: {N_DOWNLOAD_WORKERS}")
    logger.info(f"  DOY processors: {N_DOY_PROCESSORS} (each with {WORKERS_PER_DOY} RTKLIB workers)")
    logger.info(f"  Total CPU utilization: {N_DOY_PROCESSORS * WORKERS_PER_DOY} workers")
    
    # Shared state
    total_success = 0
    total_failed = 0
    total_exists = 0
    download_queue = Queue()  # DOYs to download
    process_queue = Queue()   # Downloaded files ready to process
    
    # Fill download queue
    for year, doy in doys_to_process:
        download_queue.put((year, doy, nav_files.get((year, doy))))
    
    # Signal end of downloads
    for _ in range(N_DOWNLOAD_WORKERS):
        download_queue.put(None)
    
    def download_worker():
        """Download DOYs in parallel."""
        while True:
            item = download_queue.get()
            if item is None:
                process_queue.put(None)  # Signal this worker is done
                break
            
            year, doy, nav_file = item
            if nav_file is None:
                continue
            
            rinex_files = batch_download_rinex_for_doy(year, doy, auth)
            # Always put in queue (even if empty/skipped) so progress bar updates
            process_queue.put((year, doy, nav_file, rinex_files))
    
    # Start download workers
    download_threads = []
    for _ in range(N_DOWNLOAD_WORKERS):
        t = Thread(target=download_worker, daemon=True)
        t.start()
        download_threads.append(t)
    
    # Process files as they become available - 8 DOYs in parallel
    from threading import Lock
    
    workers_per_doy = WORKERS_PER_DOY
    
    workers_done = 0
    lock = Lock()
    pbar = tqdm(total=len(doys_to_process), desc="Processing", unit="day")
    
    def process_one_doy():
        """Process DOYs from queue until done."""
        nonlocal workers_done, total_success, total_failed, total_exists
        
        while True:
            item = process_queue.get()
            
            if item is None:
                with lock:
                    workers_done += 1
                break
            
            year, doy, nav_file, rinex_files = item
            
            # Log start of processing
            n_rinex = len(rinex_files) if rinex_files else 0
            msg = f"    [{year}/{doy:03d}] PROCESSING: {n_rinex} RINEX files..."
            print(msg, flush=True)
            logger.info(msg)
            
            if not rinex_files or n_rinex == 0:
                print(f"    [{year}/{doy:03d}] SKIPPED: No RINEX files to process", flush=True)
                with lock:
                    pbar.update(1)
                continue
            
            # Find SP3 file for this DOY (for precise/ionofree modes)
            sp3_file = SP3_DIR / f"sp3_{year}{doy:03d}.sp3"
            if not sp3_file.exists():
                sp3_file = None  # Will skip precise/ionofree modes
            
            # Build processing tasks
            tasks = []
            for rinex_path in rinex_files:
                tasks.append((
                    rinex_path,
                    nav_file,
                    sp3_file,  # SP3 file for precise/ionofree modes
                    DEFAULT_RTKLIB,
                    PROCESSED_DIR
                ))
            
            # Process in parallel
            doy_new, doy_exist, doy_fail = 0, 0, 0
            doy_errors = {}
            
            process_start = _time.time()
            
            with ProcessPoolExecutor(max_workers=workers_per_doy) as executor:
                futures = [executor.submit(process_rinex_file, t) for t in tasks]
                completed = 0
                total = len(futures)
                
                for future in as_completed(futures):
                    completed += 1
                    # Log progress every 50 files
                    if completed % 50 == 0:
                        print(f"    [{year}/{doy:03d}] Progress: {completed}/{total} files...", flush=True)
                        
                    try:
                        result = future.result()
                        if result:
                            status = result["status"]
                            if status == "success":
                                doy_new += 1
                                with lock:
                                    total_success += 1
                                    if result.get("coords"):
                                        all_coords[result["station"]] = result["coords"]
                            elif status == "exists":
                                doy_exist += 1
                                with lock:
                                    total_exists += 1
                            else:
                                doy_fail += 1
                                with lock:
                                    total_failed += 1
                                doy_errors[status] = doy_errors.get(status, 0) + 1
                    except Exception as e:
                        doy_fail += 1
                        with lock:
                            total_failed += 1
                        doy_errors["exception"] = doy_errors.get("exception", 0) + 1
            
            process_time = _time.time() - process_start
            
            # Log DOY summary - write directly to log file for reliability
            msg = f"    [{year}/{doy:03d}] PROCESSED: {doy_new} new, {doy_exist} exist, {doy_fail} fail ({process_time:.1f}s)"
            print(msg, flush=True)
            logger.info(msg)
            if doy_errors:
                error_str = ", ".join(f"{k}:{v}" for k, v in doy_errors.items())
                err_msg = f"             Failures: {error_str}"
                print(err_msg, flush=True)
                logger.warning(err_msg)
            
            # Clean up
            doy_dir = RINEX_DIR / f"{year}" / f"{doy:03d}"
            if doy_dir.exists():
                shutil.rmtree(doy_dir, ignore_errors=True)
            
            with lock:
                pbar.update(1)
                pbar.set_postfix({"new": total_success, "exist": total_exists, "fail": total_failed})
    
    # Start DOY processor threads
    processor_threads = []
    for _ in range(N_DOY_PROCESSORS):
        t = Thread(target=process_one_doy, daemon=True)
        t.start()
        processor_threads.append(t)
    
    # Wait for all download workers to signal done
    while workers_done < N_DOWNLOAD_WORKERS:
        _time.sleep(1)
    
    # Wait for processor threads to finish
    for t in processor_threads:
        t.join()
    
    pbar.close()
    
    # Wait for download threads
    for t in download_threads:
        t.join()
    
    save_coordinates(all_coords)
    
    # ========================================================================
    # PHASE 3: Summary
    # ========================================================================
    logger.process("[3/3] Generating summary...")
    
    npz_files = list(PROCESSED_DIR.glob("*.npz"))
    total_size = sum(f.stat().st_size for f in npz_files)
    unique_stations = set(f.stem.split('_')[0] for f in npz_files)
    
    logger.info("")
    logger.info("=" * 60)
    logger.info("Processing Summary:")
    logger.success(f"  New processed: {total_success}")
    logger.info(f"  Already existed: {total_exists}")
    if total_failed > 0:
        logger.warning(f"  Failed: {total_failed}")
    logger.info("")
    logger.info("Total in data/processed/:")
    logger.info(f"  Files: {len(npz_files)}")
    logger.info(f"  Size: {total_size / (1024*1024):.1f} MB")
    logger.info(f"  Stations: {len(unique_stations)}")
    logger.info(f"  Coordinates: {len(all_coords)}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
