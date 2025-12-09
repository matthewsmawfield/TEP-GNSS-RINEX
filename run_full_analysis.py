#!/usr/bin/env python3
"""
TEP-GNSS-RINEX Analysis Pipeline Orchestrator
=============================================

End-to-end orchestration of the complete TEP analysis pipeline.
Executes Step 1.0 (data acquisition) followed by all Step 2.x analyses
with configurable station filters.

Usage:
    python run_full_analysis.py --filters optimal_100 dynamic_50
    python run_full_analysis.py --skip-step1 --filters all

Arguments:
    --filters: Station filter(s) to use (all, optimal_100, dynamic_50)
    --skip-step1: Skip data acquisition (use existing NPZ files)
    --skip-step2: Skip analysis (run acquisition only)

Author: Matthew Lukin Smawfield
Date: 9 December 2025
License: CC-BY-4.0
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent
STEP_DIR = PROJECT_ROOT / "scripts" / "steps"

STEP_SCRIPTS = {
    "step_1_0": STEP_DIR / "step_1_0_data_acquisition.py",
    "step_1_1_dynamic50": STEP_DIR / "step_1_1_generate_dynamic_50_metadata.py",
    "step_2_0": STEP_DIR / "step_2_0_raw_spp_analysis.py",
    "step_2_1_control": STEP_DIR / "step_2_1_control_tests.py",
    "step_2_1_elevation": STEP_DIR / "step_2_1_elevation_analysis.py",
    "step_2_2": STEP_DIR / "step_2_2_anisotropy_analysis.py",
    "step_2_3_kp": STEP_DIR / "step_2_3_kp_stratification.py",
    "step_2_3_temporal": STEP_DIR / "step_2_3_temporal_analysis.py",
    "step_2_4_diurnal": STEP_DIR / "step_2_4_diurnal_analysis.py",
    "step_2_4_null": STEP_DIR / "step_2_4_null_tests.py",
    "step_2_5": STEP_DIR / "step_2_5_orbital_coupling.py",
    "step_2_6": STEP_DIR / "step_2_6_planetary_events.py",
    "step_2_7": STEP_DIR / "step_2_7_cmb_frame_analysis.py",
}

STEP2_SEQUENCE = [
    ("STEP 2.0 • Raw SPP Analysis", "step_2_0"),
    ("STEP 2.1 • Control Tests", "step_2_1_control"),
    ("STEP 2.1 • Elevation Analysis", "step_2_1_elevation"),
    ("STEP 2.2 • Anisotropy & Orbital Coupling", "step_2_2"),
    ("STEP 2.3 • Kp Stratification", "step_2_3_kp"),
    ("STEP 2.3 • Temporal Analysis", "step_2_3_temporal"),
    ("STEP 2.4 • Diurnal Analysis", "step_2_4_diurnal"),
    ("STEP 2.4 • Null Tests", "step_2_4_null"),
    ("STEP 2.5 • Orbital Coupling Deep Dive", "step_2_5"),
    ("STEP 2.6 • Planetary Events", "step_2_6"),
    ("STEP 2.7 • CMB Frame Analysis", "step_2_7"),
]

FILTER_ALIASES = {
    "all": "none",
    "none": "none",
    "all_stations": "none",
    "optimal": "optimal_100_metadata.json",
    "optimal_100": "optimal_100_metadata.json",
    "dynamic50": "dynamic_50_metadata.json",
    "dynamic_50": "dynamic_50_metadata.json",
}


def resolve_filter(token: str) -> Tuple[str, str]:
    """Return (display_name, env_value) for a filter token."""
    normalized = token.lower()
    env_value = FILTER_ALIASES.get(normalized, token)
    display = token
    return display, env_value


def run_step(label: str, cmd: List[str], env: Dict[str, str] | None = None) -> None:
    print("\n" + "=" * 80)
    print(f"{datetime.now():%Y-%m-%d %H:%M:%S} | {label}")
    print("=" * 80)
    subprocess.run(cmd, check=True, env=env)


def ensure_scripts_exist() -> None:
    missing = [name for name, path in STEP_SCRIPTS.items() if not path.exists()]
    if missing:
        pretty = ", ".join(missing)
        raise FileNotFoundError(f"Missing step scripts: {pretty}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Step 1.0 + Step 2.x pipeline locally (mirrors run_rinex_gcp.sh)."
    )
    parser.add_argument(
        "--filters",
        nargs="+",
        default=["optimal_100"],
        help="Station filters to run (aliases: all, optimal, dynamic50, or JSON filenames).",
    )
    parser.add_argument(
        "--skip-step1",
        action="store_true",
        help="Skip Step 1.0 data acquisition (assumes NPZ files already exist).",
    )
    parser.add_argument(
        "--skip-step2",
        action="store_true",
        help="Skip Step 2.x analyses (useful for running Step 1 only).",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python interpreter to use when invoking sub-scripts.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ensure_scripts_exist()

    if not args.skip_step1:
        run_step(
            "STEP 1.0 • Data Acquisition",
            [args.python, str(STEP_SCRIPTS["step_1_0"])],
        )
        run_step(
            "STEP 1.1 • Generate DYNAMIC_50 metadata",
            [args.python, str(STEP_SCRIPTS["step_1_1_dynamic50"])],
        )
    else:
        print("Skipping Step 1.0 (per --skip-step1)")

    if args.skip_step2:
        print("Skipping Step 2.x analyses (per --skip-step2)")
        return

    for filter_token in args.filters:
        display, env_value = resolve_filter(filter_token)
        env = os.environ.copy()
        env["STATION_FILTER"] = env_value

        if env_value not in {"none", ""}:
            filter_path = PROJECT_ROOT / "data" / "processed" / env_value
            if not filter_path.exists():
                print(
                    f"⚠️  Warning: {filter_path} not found. "
                    "Step scripts may fall back to ALL_STATIONS."
                )

        for label, key in STEP2_SEQUENCE:
            run_step(
                f"{label} ({display})",
                [args.python, str(STEP_SCRIPTS[key])],
                env=env,
            )

    print("\nPipeline complete ✅")


if __name__ == "__main__":
    main()
