#!/usr/bin/env python3
"""
TEP-GNSS Cross-Paper Synthesis — STEP 5.0
=========================================

Compares correlation-length results across the PPP pipeline (TEP-GNSS, Paper 1)
and the raw-RINEX SPP pipeline (TEP-GNSS-RINEX, Paper 3) to assess:

1. Whether the exponential decay signature survives the change from
   precise-orbit PPP residuals to raw broadcast-ephemeris SPP solutions.

2. Whether λ shifts between processing chains are directionally consistent
   with increased uncorrelated noise in raw data (noise dilution should
   shorten apparent λ and/or reduce R²).

3. Whether anisotropy ratios (E-W vs N-S) are preserved across chains.

Inputs
------
- ../TEP-GNSS/results/outputs/step_2_0_correlation_analysis_summary.json
- ../TEP-GNSS/results/outputs/step_2_2_geospatial_temporal_analysis_*.json
- results/outputs/step_2_0_raw_spp_analysis_*.json
- results/outputs/step_2_2_anisotropy_analysis_*.json

Outputs
-------
- results/outputs/step_5_0_cross_paper_synthesis.json
- results/figures/step_5_0_*.png (optional, if matplotlib available)

Author: Matthew Lukin Smawfield
Date: 4 June 2026
License: CC-BY-4.0
"""

import json
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any, Optional

import numpy as np

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
WORKSPACE_ROOT = PROJECT_ROOT.parents[0]
TEP_GNSS_ROOT = WORKSPACE_ROOT / "TEP-GNSS"

sys.path.insert(0, str(PROJECT_ROOT))

RESULTS_DIR = PROJECT_ROOT / "results" / "outputs"
FIGURES_DIR = PROJECT_ROOT / "results" / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_json(path: Path) -> Optional[Dict]:
    if not path.exists():
        return None
    with open(path, "r") as f:
        return json.load(f)


def load_jaipur_correlations() -> Dict[str, Dict]:
    """Load TEP-GNSS (Jaipur) step-2.0 correlation results per center."""
    path = TEP_GNSS_ROOT / "results" / "outputs" / "step_2_0_correlation_analysis_summary.json"
    data = load_json(path)
    if data is None:
        raise FileNotFoundError(f"Jaipur correlation summary not found: {path}")
    return data


def load_jaipur_anisotropy(center: str = "code") -> Optional[Dict]:
    """Load TEP-GNSS (Jaipur) step-2.2 geospatial temporal analysis."""
    path = TEP_GNSS_ROOT / "results" / "outputs" / f"step_2_2_geospatial_temporal_analysis_{center}.json"
    return load_json(path)


def load_kathmandu_spp(filter_mode: str) -> Optional[Dict]:
    """Load TEP-GNSS-RINEX (Kathmandu) step-2.0 raw SPP results."""
    path = RESULTS_DIR / f"step_2_0_raw_spp_analysis_{filter_mode}.json"
    return load_json(path)


def load_kathmandu_anisotropy(filter_mode: str) -> Optional[Dict]:
    """Load TEP-GNSS-RINEX (Kathmandu) step-2.2 anisotropy results."""
    path = RESULTS_DIR / f"step_2_2_anisotropy_analysis_{filter_mode}.json"
    return load_json(path)

# ---------------------------------------------------------------------------
# Extraction helpers
# ---------------------------------------------------------------------------

def extract_jaipur_lambda(data: Dict, use_exponential: bool = True) -> Dict[str, Any]:
    """
    Extract λ, R², amplitude, offset from a Jaipur center result.
    Jaipur results have both `best_fit` (AIC-best model) and `exponential_fit`.
    For cross-paper consistency we default to the *exponential* fit because
    Kathmandu always fits an exponential.
    """
    key = "exponential_fit" if use_exponential else "best_fit"
    fit = data.get(key, {})
    return {
        "model": fit.get("model", key),
        "lambda_km": fit.get("lambda_km"),
        "lambda_error": fit.get("lambda_error"),
        "amplitude": fit.get("amplitude"),
        "amplitude_error": fit.get("amplitude_error"),
        "offset": fit.get("offset"),
        "offset_error": fit.get("offset_error"),
        "r_squared": fit.get("r_squared"),
        "n_bins": fit.get("n_bins"),
        "n_pairs": data.get("data_summary", {}).get("total_pairs"),
        "distance_range": data.get("data_summary", {}).get("distance_range_km"),
    }


def extract_kathmandu_lambdas(data: Dict) -> List[Dict[str, Any]]:
    """
    Extract all (mode, metric, λ, R²) combinations from Kathmandu step-2.0 JSON.
    Returns a flat list for easier matching.
    """
    results = []
    analysis_by_mode = data.get("analysis_by_mode", {})
    for mode, metrics in analysis_by_mode.items():
        for metric_key, values in metrics.items():
            # Determine metric type from suffix or explicit field
            metric_type = values.get("metric_type", "unknown")
            # Distinguish phase_alignment vs normalized_msc
            is_phase = metric_type == "phase_alignment" or metric_key.endswith("_phase_alignment")
            results.append({
                "filter": data.get("filter_mode", "unknown"),
                "mode": mode,
                "metric_key": metric_key,
                "metric_type": "phase_alignment" if is_phase else "normalized_msc",
                "lambda_km": values.get("correlation_length_km"),
                "amplitude": values.get("amplitude"),
                "offset": values.get("offset"),
                "r_squared": values.get("r_squared"),
                "n_pairs": values.get("n_pairs"),
                "n_bins": values.get("n_bins"),
            })
    return results


def extract_kathmandu_anisotropy(data: Dict) -> Dict[str, Any]:
    """Extract E-W / N-S pooled ratios from Kathmandu anisotropy JSON."""
    out = {}
    mma = data.get("multi_metric_analysis", {})
    mr = mma.get("metric_results", {})
    for metric_name, metric_data in mr.items():
        pooled = metric_data.get("ew_ns_pooled", {})
        if pooled:
            ew = pooled.get("ew", {})
            ns = pooled.get("ns", {})
            ratio = pooled.get("ratio")
            out[metric_name] = {
                "ew_lambda_km": ew.get("lambda_km"),
                "ns_lambda_km": ns.get("lambda_km"),
                "ew_ns_ratio": ratio,
                "ew_r_squared": ew.get("r_squared"),
                "ns_r_squared": ns.get("r_squared"),
            }
    return out


def extract_jaipur_anisotropy(data: Dict) -> Dict[str, Any]:
    """Extract anisotropy from Jaipur geospatial temporal analysis."""
    out = {}
    # Jaipur step_2_2 stores results under keys like "w30", "w60", etc.
    # We look for the overall summary if present.
    for key, val in data.items():
        if not isinstance(val, dict):
            continue
        if "ew_ns_ratio" in val or "anisotropy" in key.lower():
            out[key] = val
    return out

# ---------------------------------------------------------------------------
# Comparison logic
# ---------------------------------------------------------------------------

def compute_comparison(jaipur: Dict, kathmandu: Dict) -> Dict[str, Any]:
    """Compute ratio, difference, and consistency flag for a matched pair."""
    j_lam = jaipur["lambda_km"]
    k_lam = kathmandu["lambda_km"]
    ratio = k_lam / j_lam if j_lam else None
    diff = k_lam - j_lam if j_lam is not None else None

    j_r2 = jaipur.get("r_squared", 0)
    k_r2 = kathmandu.get("r_squared", 0)

    # Consistent if both R² > 0.5 and λ within a factor of 3
    consistent = False
    if j_r2 and k_r2:
        consistent = (j_r2 > 0.5) and (k_r2 > 0.5) and (ratio is not None) and (0.33 < ratio < 3.0)

    return {
        "lambda_ratio_raw_to_ppp": ratio,
        "lambda_difference_km": diff,
        "pp_r_squared_ppp": j_r2,
        "pp_r_squared_raw": k_r2,
        "qualitatively_consistent": consistent,
    }


def build_primary_comparisons(jaipur_data: Dict, kathmandu_results: List[Dict]) -> List[Dict]:
    """
    Build the primary comparison table.
    Jaipur provides one λ per analysis center (phase_alignment_index).
    Kathmandu provides many λs (mode × metric × filter).
    We match Jaipur phase_alignment → Kathmandu phase_alignment metrics.
    """
    comparisons = []
    centers = ["code", "igs_combined", "esa_final"]

    for center in centers:
        if center not in jaipur_data:
            continue
        j = extract_jaipur_lambda(jaipur_data[center], use_exponential=True)

        for k in kathmandu_results:
            # Only compare phase_alignment metrics (same family as Jaipur)
            if k["metric_type"] != "phase_alignment":
                continue

            # The "precise" mode in Kathmandu uses precise ephemerides — closest to PPP
            proximity = {
                "precise": 1.0,
                "ionofree": 0.8,
                "multi_gnss": 0.6,
                "baseline": 0.4,
            }.get(k["mode"], 0.0)

            comp = compute_comparison(j, k)
            comparisons.append({
                "jaipur_center": center,
                "jaipur_metric": "phase_alignment_index (exponential)",
                "kathmandu_filter": k["filter"],
                "kathmandu_mode": k["mode"],
                "kathmandu_metric_key": k["metric_key"],
                "kathmandu_metric_type": k["metric_type"],
                "processing_chain_similarity": proximity,
                "jaipur_lambda_km": j["lambda_km"],
                "jaipur_lambda_error": j.get("lambda_error"),
                "kathmandu_lambda_km": k["lambda_km"],
                **comp,
            })

    # Sort by processing-chain similarity (most comparable first)
    comparisons.sort(key=lambda x: x["processing_chain_similarity"], reverse=True)
    return comparisons


def build_msc_comparisons(kathmandu_results: List[Dict]) -> List[Dict]:
    """
    Build a supplementary table for MSC (normalized_msc) metrics.
    Jaipur step-2.0 does not explicitly report MSC λs, but if they exist
    in other steps we could match them. For now we just catalogue Kathmandu
    MSC results as a baseline for future matching.
    """
    msc_records = [k for k in kathmandu_results if k["metric_type"] == "normalized_msc"]
    return [
        {
            "kathmandu_filter": k["filter"],
            "kathmandu_mode": k["mode"],
            "kathmandu_metric_key": k["metric_key"],
            "lambda_km": k["lambda_km"],
            "r_squared": k["r_squared"],
            "n_pairs": k["n_pairs"],
        }
        for k in msc_records
    ]

# ---------------------------------------------------------------------------
# Meta-analysis
# ---------------------------------------------------------------------------

def meta_analyse(comparisons: List[Dict]) -> Dict[str, Any]:
    """Aggregate statistics across all matched pairs."""
    ratios = [c["lambda_ratio_raw_to_ppp"] for c in comparisons if c["lambda_ratio_raw_to_ppp"] is not None]
    r2_ppp = [c["pp_r_squared_ppp"] for c in comparisons if c["pp_r_squared_ppp"] is not None]
    r2_raw = [c["pp_r_squared_raw"] for c in comparisons if c["pp_r_squared_raw"] is not None]

    if not ratios:
        return {"error": "No valid comparisons"}

    ratios_arr = np.array(ratios)
    r2_ppp_arr = np.array(r2_ppp)
    r2_raw_arr = np.array(r2_raw)

    # Directional hypothesis test: raw λ should be *shorter* than PPP λ
    # because raw SPP has higher uncorrelated noise, which dilutes correlations.
    # One-sided sign test: count ratios < 1.
    n_below_1 = int(np.sum(ratios_arr < 1.0))
    n_total = len(ratios_arr)
    # Binomial p-value for H0: P(ratio < 1) = 0.5
    try:
        from scipy.stats import binomtest
        sign_test_p = binomtest(n_below_1, n_total, p=0.5, alternative="two-sided").pvalue
    except ImportError:
        from scipy.stats import binom_test
        sign_test_p = binom_test(n_below_1, n_total, p=0.5, alternative="two-sided")

    return {
        "n_comparisons": len(comparisons),
        "n_valid_ratios": n_total,
        "lambda_ratio_mean": float(np.mean(ratios_arr)),
        "lambda_ratio_std": float(np.std(ratios_arr)),
        "lambda_ratio_median": float(np.median(ratios_arr)),
        "lambda_ratio_min": float(np.min(ratios_arr)),
        "lambda_ratio_max": float(np.max(ratios_arr)),
        "ratios_below_1_count": n_below_1,
        "ratios_below_1_fraction": float(n_below_1 / n_total) if n_total > 0 else None,
        "sign_test_pvalue_two_sided": float(sign_test_p),
        "interpretation_sign_test": (
            "Raw λ systematically shorter than PPP λ (p < 0.05)"
            if sign_test_p < 0.05 and n_below_1 > n_total / 2
            else "No significant systematic λ shift detected"
        ),
        "r2_ppp_mean": float(np.mean(r2_ppp_arr)),
        "r2_ppp_std": float(np.std(r2_ppp_arr)),
        "r2_raw_mean": float(np.mean(r2_raw_arr)),
        "r2_raw_std": float(np.std(r2_raw_arr)),
        "qualitatively_consistent_count": sum(1 for c in comparisons if c["qualitatively_consistent"]),
        "qualitatively_consistent_fraction": sum(1 for c in comparisons if c["qualitatively_consistent"]) / len(comparisons),
    }

# ---------------------------------------------------------------------------
# Anisotropy cross-check
# ---------------------------------------------------------------------------

def anisotropy_cross_check(kathmandu_aniso: Dict) -> Dict[str, Any]:
    """
    Evaluate whether Kathmandu anisotropy ratios are internally consistent
    and compare against the CODE benchmark ratio of 2.16.
    """
    code_ref = 2.16
    findings = {}
    for metric, data in kathmandu_aniso.items():
        ratio = data.get("ew_ns_ratio")
        if ratio is None:
            continue
        findings[metric] = {
            "ew_ns_ratio": ratio,
            "code_reference_ratio": code_ref,
            "ratio_vs_code": ratio / code_ref if code_ref else None,
            "ew_r_squared": data.get("ew_r_squared"),
            "ns_r_squared": data.get("ns_r_squared"),
        }
    return findings

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 60)
    print("STEP 5.0: Cross-Paper Synthesis (PPP vs Raw RINEX)")
    print("=" * 60)

    # ---- Load Jaipur (PPP) data ------------------------------------------------
    print("\n[1/4] Loading Jaipur (TEP-GNSS) PPP results ...")
    jaipur_corr = load_jaipur_correlations()
    print(f"      Centers available: {list(jaipur_corr.keys())}")

    jaipur_aniso = {}
    for center in ["code", "igs_combined", "esa_final"]:
        data = load_jaipur_anisotropy(center)
        if data:
            jaipur_aniso[center] = extract_jaipur_anisotropy(data)

    # ---- Load Kathmandu (Raw) data ---------------------------------------------
    print("\n[2/4] Loading Kathmandu (TEP-GNSS-RINEX) raw SPP results ...")
    filters = ["all_stations", "dynamic_50", "optimal_100"]
    kathmandu_results = []
    kathmandu_aniso = {}
    for filt in filters:
        spp = load_kathmandu_spp(filt)
        if spp:
            kathmandu_results.extend(extract_kathmandu_lambdas(spp))
            print(f"      {filt}: {len(extract_kathmandu_lambdas(spp))} metric entries")
        aniso = load_kathmandu_anisotropy(filt)
        if aniso:
            kathmandu_aniso[filt] = extract_kathmandu_anisotropy(aniso)

    # ---- Build comparisons -------------------------------------------------------
    print("\n[3/4] Building matched comparisons ...")
    primary_comparisons = build_primary_comparisons(jaipur_corr, kathmandu_results)
    msc_catalogue = build_msc_comparisons(kathmandu_results)
    meta = meta_analyse(primary_comparisons)
    aniso_findings = anisotropy_cross_check(
        kathmandu_aniso.get("all_stations", {})
    )

    # ---- Assemble output ---------------------------------------------------------
    output = {
        "step": "5.0",
        "name": "Cross-Paper Synthesis: PPP vs Raw RINEX",
        "timestamp": datetime.now().isoformat(),
        "jaipur_summary": {
            "paper": "1-TEP-GNSS-v0.25-Jaipur",
            "processing": "Precise Point Positioning (PPP) with CODE/ESA/IGS orbit-clock products",
            "primary_metric": "phase_alignment_index (exponential fit)",
            "centers": {
                c: {
                    "lambda_km": extract_jaipur_lambda(jaipur_corr[c], use_exponential=True)["lambda_km"],
                    "r_squared": extract_jaipur_lambda(jaipur_corr[c], use_exponential=True)["r_squared"],
                }
                for c in jaipur_corr.keys()
            },
        },
        "kathmandu_summary": {
            "paper": "3-TEP-GNSS-RINEX-v0.5-Kathmandu",
            "processing": "Single Point Positioning (SPP) with broadcast ephemerides",
            "filters": filters,
            "n_metric_entries_total": len(kathmandu_results),
        },
        "primary_comparisons": primary_comparisons,
        "msc_catalogue": msc_catalogue,
        "meta_analysis": meta,
        "anisotropy_findings": aniso_findings,
        "reviewer_relevance": {
            "forward_model_critique": (
                "If λ shifts systematically and predictably when moving from PPP to raw SPP, "
                "the space of possible forward models is constrained. "
                "A noise-dilution model predicts shorter λ in raw data; "
                "observing this directionally supports a physical (not software-artifact) origin."
            ),
            "multiplicity_critique": (
                "This synthesis treats all 72 Kathmandu metric combinations as consistency checks, "
                "not independent tests. The primary comparison is the phase-alignment index "
                "across the most similar processing chain (precise ephemerides)."
            ),
        },
    }

    out_path = RESULTS_DIR / "step_5_0_cross_paper_synthesis.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\n[4/4] Written: {out_path}")

    # ---- Console summary ---------------------------------------------------------
    print("\n" + "=" * 60)
    print("SYNTHESIS SUMMARY")
    print("=" * 60)
    print(f"Comparisons built: {len(primary_comparisons)}")
    print(f"Mean λ ratio (raw / PPP): {meta.get('lambda_ratio_mean', 'N/A'):.3f}")
    print(f"Median λ ratio: {meta.get('lambda_ratio_median', 'N/A'):.3f}")
    print(f"Ratios < 1 (raw shorter): {meta.get('ratios_below_1_count', 0)} / {meta.get('n_valid_ratios', 0)}")
    print(f"Sign-test p-value: {meta.get('sign_test_pvalue_two_sided', 'N/A'):.4f}")
    print(f"Qualitatively consistent: {meta.get('qualitatively_consistent_count', 0)} / {meta.get('n_comparisons', 0)}")
    print(f"Mean R² (PPP): {meta.get('r2_ppp_mean', 'N/A'):.3f}")
    print(f"Mean R² (raw): {meta.get('r2_raw_mean', 'N/A'):.3f}")
    print("\nTop 5 closest processing-chain matches:")
    for c in primary_comparisons[:5]:
        print(
            f"  {c['jaipur_center']:12s} → {c['kathmandu_filter']}/{c['kathmandu_mode']}/{c['kathmandu_metric_key']}: "
            f"λ_PPP={c['jaipur_lambda_km']:.0f} km, λ_raw={c['kathmandu_lambda_km']:.0f} km, "
            f"ratio={c['lambda_ratio_raw_to_ppp']:.3f}"
        )
    print("=" * 60)

    # Optional figure
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 6))
        labels = [f"{c['jaipur_center']}\n→ {c['kathmandu_mode'][:4]}.{c['kathmandu_metric_key'][-8:]}" for c in primary_comparisons[:15]]
        ratios = [c["lambda_ratio_raw_to_ppp"] for c in primary_comparisons[:15]]
        colors = ["#2ecc71" if c["qualitatively_consistent"] else "#e74c3c" for c in primary_comparisons[:15]]
        ax.barh(range(len(labels)), ratios, color=colors)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=7)
        ax.axvline(x=1.0, color="black", linestyle="--", linewidth=1, label="λ_raw = λ_PPP")
        ax.set_xlabel("λ_raw / λ_PPP")
        ax.set_title("Cross-Paper λ Comparison: Raw RINEX / PPP")
        ax.invert_yaxis()
        fig.tight_layout()
        fig_path = FIGURES_DIR / "step_5_0_lambda_ratio_comparison.png"
        fig.savefig(fig_path, dpi=150)
        print(f"Figure saved: {fig_path}")
    except Exception as e:
        print(f"Figure generation skipped: {e}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
