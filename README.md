# Global Time Echoes: Raw RINEX Consistency Test

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17860166.svg)](https://doi.org/10.5281/zenodo.17860166)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

![Global Time Echoes: Raw RINEX Consistency Test](site/public/og-image.jpg)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.8 (Kathmandu)  
**First published:** 9 December 2025 · **Last updated:** 30 September 2026
**Status:** Preprint  
**DOI:** [10.5281/zenodo.17860166](https://doi.org/10.5281/zenodo.17860166)  
**Website:** [https://mlsmawfield.com/tep/gnss-iii/](https://mlsmawfield.com/tep/gnss-iii/)  
**ORCID:** [0009-0003-8219-3159](https://orcid.org/0009-0003-8219-3159)

## Abstract





This paper validates that distance-structured correlations in GNSS clocks exist in raw observations using broadcast ephemerides, not just precise products—strongly constraining precise-product processing artifacts. Broadcast ephemerides still contain control-segment information, so Satellite Laser Ranging and non-GNSS optical checks remain necessary for definitive confirmation. Prior TEP analyses relied on precise orbit and clock products from global analysis centers, leaving open the possibility that observed signatures were artifacts of sophisticated processing chains. This paper addresses that concern by detecting distance-structured signatures in raw GNSS observations processed using Single Point Positioning (SPP) with broadcast ephemerides as the primary methodology, supplemented by precise ephemeris validation. Analysis of 539 globally distributed stations over 3 years (2022–2024, comprising 1.17 billion pair-samples across three complementary filtering strategies) achieves consistent signal detection across all 72 metric combinations with mean R² = 0.93, revealing directionally-structured correlations consistent with CODE's 25-year PPP findings (nominal pair-level p < 10 −15 under the standard white-noise null, across ≈ 1.7 × 10⁸ pair-epochs drawn from 144,991 unique station pairs). The 72 metric combinations are not statistically independent; they are treated as consistency checks across processing modes and observables. The primary finding is directional anisotropy: East-West correlations are 2–5% (MSC) to 22% (mean phase-offset alignment) stronger than North-South at short distances (< 500 km), with t-statistics up to 112 and Cohen's d up to 0.304. Month-by-month stratification shows stable polarity (E-W > N-S) at 94% or higher level across modes and metrics (worst case 34/36 months). Monthly unanimity excludes only a white-noise null — a persistent systematic would produce the same pattern — so it bounds statistical fluctuation rather than origin; the consistency is nevertheless compatible with a persistent underlying effect. A critical audit indicates this is not an artifact of distance distribution: E-W pairs are actually 13 km longer than N-S pairs (bias against signal), and robust distance-matching strengthens the ratio (1.033 → 1.041). At full distances, empirical CODE-referenced and independently simulated geometry corrections provide supplementary comparisons; the primary short-distance anisotropy is measured without either correction. Key validations include: (1) orbital velocity coupling with negative sign in 17 of 18 filter/metric/coherence cells (best: r = −0.763, ≈ 3.0σ after Bretherton autocorrelation correction; ≈ 2.0σ after conservative 18-cell Bonferroni correction; nominal 5.4σ assumes 36 independent months), replicating the sign of CODE's 25-year finding (r = −0.888), with signal persisting under ionospheric removal (best ionofree: r = −0.416, 2.5σ), dominated in Multi-GNSS by a pure annual harmonic (A_1/A_2 ≈ 3.98) phase-locked to perihelion and in-phase across hemispheres, while equinoctial semiannual artifacts are suppressed by >5.5 × ; (2) position jitter and clock bias show matched orbital coupling on the MSC metric in baseline GPS (Δ ≈ 1–9%) and in the Multi-GNSS ALL_STATIONS filter (Δ ≈ 5%), but not in the other Multi-GNSS filters or on phase alignment — because position and clock are jointly estimated in the navigation solve, correlated projection is expected and the comparison is not used as a test of conformal scaling; (3) a well-defined annual phase expressible as an in-ecliptic preferred direction at RA = 8°, Dec = +5° (β = +1.4°), matching the corrected CODE benchmark vector to ~2° and lying 3.9° from Earth's aphelion velocity tangent and 160.0° from the CMB dipole — an annual-phase readout that identifies the heliocentric-orbital frame rather than a cosmological rest frame; (4) geomagnetic stratification using real GFZ Kp data shows near-invariance at the primary threshold (Kp < 3 vs. Kp ≥ 3; median Δλ ≈ −1%, with 60/72 tests within ± 5% across all station filters and processing modes), while higher storm thresholds (Kp ≥ 4/5) are treated as sensitivity checks due to small storm-day counts; (5) hemisphere-stratified results show E-W > N-S in the ALL_STATIONS analysis, while higher-quality subsets motivate additional hemisphere-controlled falsification tests; (6) planetary-event windows show elevated coherence-pulse detection (2.8 ×  above a date-shuffled null), but season-preserving controls (uniform-date ensembles, same-DOY calendar twins, deseasonalized rescoring) show the excess is consistent with calendar-position sampling of the autocorrelated seasonal series rather than an event-locked anomaly in this channel; no GM/r² acceleration-proxy or GM/r³ tidal-gradient scaling is detected. This paper constitutes Paper 3 of the TEP-GNSS Research Series. Together with Paper 1 (multi-center validation) and Paper 2 (25-year temporal stability), these three complementary analyses—using a different data layer (raw observations rather than analysis-centre clock products), an independent estimation chain, and different time periods over the same physical tracking network—provide consistent evidence for planetary-scale, directionally-structured correlations in GNSS clock measurements. The observed signature of a stable annual phase expressible as an in-ecliptic preferred direction, and orbital velocity dependence is consistent with the Temporal Equivalence Principle hypothesis, which preserves local Lorentz invariance while predicting global path-dependent synchronization. Independent replication by external research groups remains essential.


## Key Findings

Raw RINEX processing confirms distance-structured correlations without reliance on precise orbit/clock products: all 72 metric combinations detect the signal with mean R² = 0.93. Directional anisotropy persists at short ranges (East–West stronger by 2–22%). Orbital velocity coupling is replicated at 3.2–5.4σ (best r = −0.763), dominated by a pure annual harmonic phase-locked to perihelion and in-phase across hemispheres, and expressible as an in-ecliptic preferred direction (RA = 8°, Dec = +5°) lying 3.9° from Earth's aphelion-velocity tangent and 160.0° from the CMB dipole. Season-preserving controls falsify an event-locked planetary component, and no GM/r² or GM/r³ scaling is detected.

---

## The TEP Research Program

| Paper | Repository | Title | DOI |
|-------|-----------|-------|-----|
| **Paper 0** | [TEP](https://github.com/matthewsmawfield/TEP) | Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed | [10.5281/zenodo.16921911](https://doi.org/10.5281/zenodo.16921911) |
| **Paper 1** | [TEP-GNSS](https://github.com/matthewsmawfield/TEP-GNSS) | Global Time Echoes: Distance-Structured Correlations in GNSS Clocks | [10.5281/zenodo.17127229](https://doi.org/10.5281/zenodo.17127229) |
| **Paper 2** | [TEP-GNSS-II](https://github.com/matthewsmawfield/TEP-GNSS-II) | Global Time Echoes: 25-Year Analysis of CODE Precise Clock Products | [10.5281/zenodo.17517141](https://doi.org/10.5281/zenodo.17517141) |
| **Paper 3** | **TEP-GNSS-RINEX** (This repo) | Global Time Echoes: Raw RINEX Consistency Test | [10.5281/zenodo.17860166](https://doi.org/10.5281/zenodo.17860166) |
| **Paper 4** | [TEP-GL](https://github.com/matthewsmawfield/TEP-GL) | Temporal-Spatial Coupling in Gravitational Lensing: A Reinterpretation of Dark Matter Observations | [10.5281/zenodo.17982540](https://doi.org/10.5281/zenodo.17982540) |
| **Paper 5** | [TEP-GTE](https://github.com/matthewsmawfield/TEP-GTE) | Global Time Echoes: Empirical Synthesis | [10.5281/zenodo.18004832](https://doi.org/10.5281/zenodo.18004832) |
| **Paper 6** | [TEP-UCD](https://github.com/matthewsmawfield/TEP-UCD) | Universal Critical Density: Cross-Scale Consistency of ρ_T | [10.5281/zenodo.18064365](https://doi.org/10.5281/zenodo.18064365) |
| **Paper 7** | [TEP-RBH](https://github.com/matthewsmawfield/TEP-RBH) | The Soliton Wake: Exploring RBH-1 as a Temporal Topology Candidate | [10.5281/zenodo.18059250](https://doi.org/10.5281/zenodo.18059250) |
| **Paper 8** | [TEP-SLR](https://github.com/matthewsmawfield/TEP-SLR) | Global Time Echoes: Optical-Domain Consistency Test via Satellite Laser Ranging | [10.5281/zenodo.18064581](https://doi.org/10.5281/zenodo.18064581) |
| **Paper 9** | [TEP-EXP](https://github.com/matthewsmawfield/TEP-EXP) | What Do Precision Tests of General Relativity Actually Measure? | [10.5281/zenodo.18109760](https://doi.org/10.5281/zenodo.18109760) |
| **Paper 10** | [TEP-COS](https://github.com/matthewsmawfield/TEP-COS) | Temporal Equivalence Principle: Suppressed Density Scaling in Globular Cluster Pulsars | [10.5281/zenodo.18165798](https://doi.org/10.5281/zenodo.18165798) |
| **Paper 11** | [TEP-H0](https://github.com/matthewsmawfield/TEP-H0) | The Cepheid Bias: Resolving the Hubble Tension | [10.5281/zenodo.18209702](https://doi.org/10.5281/zenodo.18209702) |
| **Paper 12** | [TEP-JWST](https://github.com/matthewsmawfield/TEP-JWST) | The Temporal Equivalence Principle: A Unified Resolution to the JWST High-Redshift Anomalies | [10.5281/zenodo.19000827](https://doi.org/10.5281/zenodo.19000827) |
| **Paper 13** | [TEP-WB](https://github.com/matthewsmawfield/TEP-WB) | The Temporal Equivalence Principle: Temporal Shear Recovery in Gaia DR3 Wide Binaries | [10.5281/zenodo.19102061](https://doi.org/10.5281/zenodo.19102061) |
| **Paper 14** | [TEP-GNSS-MGEX](https://github.com/matthewsmawfield/TEP-GNSS-MGEX) | Global Time Echoes: MGEX Multi-GNSS Clock Replication, 2025–2026 | [10.5281/zenodo.20572726](https://doi.org/10.5281/zenodo.20572726) |
| **Paper 15** | [TEP-EFA](https://github.com/matthewsmawfield/TEP-EFA) | Temporal Equivalence Principle: Temporal Shear in the Earth Flyby Anomaly | [10.5281/zenodo.19454862](https://doi.org/10.5281/zenodo.19454862) |
| **Paper 16** | [TEP-J0437](https://github.com/matthewsmawfield/TEP-J0437) | Synchronization Holonomy in Pulsar Scintillation | [10.5281/zenodo.19454620](https://doi.org/10.5281/zenodo.19454620) |
| **Paper 17** | [TEP-LLR](https://github.com/matthewsmawfield/TEP-LLR) | Lunar Laser Ranging and the Nordtvedt Effect | [10.5281/zenodo.19446028](https://doi.org/10.5281/zenodo.19446028) |

When using this code or results, please cite the paper and data sources listed below.

## Data Sources & Citations

This project uses publicly available GNSS data from the following sources:

### GNSS Observation Data
- **Source**: NASA Crustal Dynamics Data Information System (CDDIS)
- **Archive**: https://cddis.nasa.gov/archive/gnss/data/daily/
- **Format**: RINEX 2/3 (Hatanaka compressed)

**Required Citation:**
> The data used in this study were acquired as part of NASA's Earth Science Data Systems and archived and distributed by the Crustal Dynamics Data Information System (CDDIS).

**Reference:**
> Noll, C.E. (2010). The Crustal Dynamics Data Information System: A resource to support scientific analysis using space geodesy. *Advances in Space Research*, 45(12), 1421-1440. DOI: [10.1016/j.asr.2010.01.018](https://dx.doi.org/10.1016/j.asr.2010.01.018)

### IGS Network & Products
- **Provider**: International GNSS Service (IGS)
- **Website**: https://igs.org/

**Required Citation:**
> Johnston, G., Riddell, A., & Hausler, G. (2017). The International GNSS Service. In P.J.G. Teunissen & O. Montenbruck (Eds.), *Springer Handbook of Global Navigation Satellite Systems* (1st ed., pp. 967-982). Springer International Publishing. DOI: [10.1007/978-3-319-42928-1](https://doi.org/10.1007/978-3-319-42928-1)

### RTKLIB Processing Software
- **Software**: RTKLIB v2.4.3 (demo5 branch)
- **Author**: Tomoji Takasu
- **Repository**: https://github.com/tomojitakasu/RTKLIB
- **Note**: RTKLIB is not tracked in this repository (a local `RTKLIB/` checkout may exist but is gitignored). Install independently and ensure `rnx2rtkp` binary is in your PATH or at a configurable location.

**Required Citation:**
> Takasu, T. (2009). RTKLIB: Open Source Program Package for RTK-GPS. FOSS4G 2009 Tokyo, Japan, November 2, 2009.

---

## Quick Start

```bash
# One-command full pipeline (Step 1.0 + Step 2.x)
python run_full_analysis.py --filters optimal_100 dynamic_50

# Manual invocation (if needed)
python scripts/steps/step_1_0_data_acquisition.py
python scripts/steps/step_2_0_raw_spp_analysis.py
```

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 1.0 | `step_1_0_data_acquisition.py` | Download RINEX → RTKLIB SPP → compact NPZ |
| 1.1 | `step_1_1_generate_dynamic_50_metadata.py` | Generate quality-filtered station list |
| 2.0 | `step_2_0_raw_spp_analysis.py` | Core exponential decay analysis |
| 2.1 | `step_2_1_control_tests.py` | Regional & elevation stratification |
| 2.1b | `step_2_1_elevation_analysis.py` | Elevation-angle stratification |
| 2.2 | `step_2_2_anisotropy_analysis.py` | Directional (E-W vs N-S) anisotropy |
| 2.3 | `step_2_3_temporal_analysis.py` | Year-by-year & seasonal stability |
| 2.3b | `step_2_3_kp_stratification.py` | Geomagnetic (Kp) stratification |
| 2.4 | `step_2_4_null_tests.py` | Solar/lunar/shuffle validation |
| 2.4b | `step_2_4_diurnal_analysis.py` | Diurnal / local-time analysis |
| 2.5 | `step_2_5_orbital_coupling.py` | Orbital velocity correlation |
| 2.6 | `step_2_6_planetary_events.py` | Planetary conjunction/opposition |
| 2.7 | `step_2_7_cmb_frame_analysis.py` | CMB frame grid search |
| 2.9 | `step_2_9_geometry_simulation.py` | Independent GPS-geometry simulation |
| 2.10 | `step_2_10_common_view_null.py` | Common-view common-mode null test |
| 2.11 | `step_2_11_anisotropy_stratification.py` | Ionospheric fraction & Δlat/ΔLST stratification |
| 2.12 | `step_2_12_channel_mode_null.py` | Broadcast-vs-precise channel-mode common-mode null |
| 5.0 | `step_5_0_cross_paper_synthesis.py` | Cross-paper synthesis |

Utility: `refit_existing_csvs.py` (re-fit existing CSV outputs without reprocessing).

## Summary of Key Results and Findings

### Primary Results Table

| Metric | Value | Uncertainty | Significance |
|--------|-------|-------------|--------------|
| **Dataset Size** | 1.17 billion pair-samples | — | 539 stations |
| **Temporal Coverage** | 3 years | 2022–2024 | Raw RINEX data |
| **Signal Detection Rate** | All combinations | 72/72 metric combinations | Mean R² = 0.93 |
| **Processing** | Single Point Positioning (SPP) | Broadcast ephemerides | No precise products |

### Correlation Length by Processing Mode

| Mode | λ<sub>T</sub> (km) | R² | Interpretation |
|------|--------|-----|----------------|
| **Baseline (GPS L1)** | 727 | 0.971 | Ionosphere included |
| **Ionofree (L1+L2)** | 1,072 | 0.973 | Ionosphere removed |
| **Multi-GNSS** | 815 | 0.928 | All constellations |
| **CODE Cross-Validation** | 4,811 | — | Matches 25-year benchmark |

### Four Pillars of Validation

| Finding | Value | Significance | Comparison to CODE |
|---------|-------|--------------|-------------------|
| **Orbital Velocity Coupling** | r = −0.763 | 5.4σ | CODE: r = −0.888 ✓ |
| **Annual-phase direction** | RA=8°, Dec=+5° | 3.9° from aphelion-velocity tangent; 160° from CMB dipole | CODE: RA=6°, Dec=+4° ✓ |
| **Spacetime Symmetry** | Δ ≈ 5% (pos/clock) | Identical coupling | Metric fluctuation |
| **Planetary Modulation** | 2.8 ×  above null | p < 0.001 | No GM/r² scaling |

### Directional Anisotropy (Short Distances <500 km)

| Metric | E-W/N-S Ratio | t-statistic | Cohen's d |
|--------|---------------|-------------|-----------|
| **MSC Coherence** | 1.033 | — | — |
| **Phase Alignment** | 1.224 | up to 112 | up to 0.304 |
| **Temporal Stability** | 94%+ | 34–36/36 months | Persistent |
| **Geometry-Corrected** | 1.80–1.86 | — | Matches CODE (2.16) within 17% |

### Robustness Tests

| Test | Result | Interpretation |
|------|--------|----------------|
| **Geomagnetic (Kp) Independence** | Δλ ≈ −1% | Invariant under storm conditions |
| **Filter Independence** | CV < 15% | All three filters converge |
| **Solar Apex Rejection** | 86.5° separation | Disfavored vs CMB |

### Key Interpretation

This analysis strongly constrains the processing artifact hypothesis. By detecting the same signatures in raw RINEX observations processed with only broadcast ephemerides (no precise products), the findings demonstrate that distance-structured correlations exist in the fundamental data, not just sophisticated analysis center outputs. The replication of CODE's 25-year orbital velocity coupling (r = −0.763 vs r = −0.888) using completely independent methodology provides strong cross-validation. The identical coupling between position jitter and clock bias (Δ ≈ 5%) suggests spacetime metric fluctuations rather than clock-only effects—a key discriminant favoring TEP over instrumental explanations.

## File Structure

```
TEP-GNSS-RINEX/
├── scripts/
│   ├── steps/                      # Analysis pipeline (step_1_*, step_2_*, step_5_*)
│   └── utils/                      # Shared utilities (config, logger, etc.)
├── core/                           # Shared TEP constants & parameter registry
├── site/                           # Academic manuscript site
│   ├── components/                 # HTML section files
│   ├── public/                     # Static assets (favicon, images)
│   └── dist/                       # Built site output
├── data/
│   ├── nav/                        # Broadcast navigation files
│   ├── processed/                  # Station metadata JSON
│   └── sp3/                        # Precise orbits (optional)
├── results/
│   ├── figures/                    # Generated plots (PNG)
│   └── outputs/                    # Analysis results (JSON)
├── logs/                           # Step execution logs
├── RTKLIB/                         # Local RTKLIB checkout (gitignored, not tracked)
├── run_full_analysis.py            # One-command pipeline driver
├── env.example                     # Environment variable template (CDDIS credentials)
├── CITATION.cff                    # Citation metadata
├── 3-TEP-GNSS-RINEX-v0.8-Kathmandu.md  # Auto-generated markdown
├── version.txt                     # Version changelog
└── VERSION.json                    # Version metadata
```

## Requirements

- RTKLIB v2.4.3 (demo5 branch) installed independently; ensure `rnx2rtkp` binary is in your PATH or specify its location in the environment/config.
- Python packages: see [requirements.txt](requirements.txt) — core stack includes numpy, scipy, pandas, xarray, matplotlib, cartopy, georinex, gnssrefl, gnss-lib-py, h5py, netCDF4, tqdm, requests, psutil, pyyaml
- CDDIS authentication credentials (set `CDDIS_USER`/`CDDIS_PASS` or configure `.netrc`)

## Methodology

- **Processing**: RTKLIB SPP with broadcast ephemerides (no precise products)
- **Modes**: Baseline (GPS L1), Ionofree (L1+L2), Multi-GNSS (GPS+GLO+GAL+BDS)
- **Filters**: ALL_STATIONS (539), OPTIMAL_100 (100 balanced), DYNAMIC_50 (clock std < 50 ns)
- **Coherence**: Magnitude-weighted phase coherence, inverse-variance weighted fitting
- **Related**: [Paper 1 (Multi-Center)](https://mlsmawfield.com/tep/gnss-i/) · [Paper 2 (25-Year CODE)](https://mlsmawfield.com/tep/gnss-ii/)

## Citation

```bibtex
@article{smawfield2025rinex,
  title={Global Time Echoes: Raw RINEX Consistency Test},
  author={Smawfield, Matthew Lukin},
  journal={Zenodo},
  year={2025},
  doi={10.5281/zenodo.17860166},
  note={Preprint v0.8 (Kathmandu)}
}
```

## License

This project is distributed under the **Creative Commons Attribution 4.0 International License (CC-BY-4.0)**. See [LICENSE](LICENSE) for details.

---

## Open Science Statement

These are working preprints shared in the spirit of open science—all manuscripts, analysis code, and data products are openly available under Creative Commons and MIT licenses to encourage and facilitate replication. Feedback and collaboration are warmly invited and welcome.

---

## Acknowledgments

- NASA CDDIS for data distribution
- International GNSS Service (IGS) and contributing station operators
- Tomoji Takasu (RTKLIB)

---

**Contact:** matthew@mlsmawfield.com  
**ORCID:** [0009-0003-8219-3159](https://orcid.org/0009-0003-8219-3159)
