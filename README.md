# Global Time Echoes: Raw RINEX Validation of Distance-Structured Correlations in GNSS Clocks

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17860166.svg)](https://doi.org/10.5281/zenodo.17860166)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

![Global Time Echoes: Raw RINEX Validation](site/public/og-image.jpg)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.3 (Kathmandu)  
**Date:** 17 December 2025  
**Status:** Preprint  
**DOI:** [10.5281/zenodo.17860166](https://doi.org/10.5281/zenodo.17860166)  
**Website:** [https://mlsmawfield.com/tep/gnss-iii/](https://mlsmawfield.com/tep/gnss-iii/)  
**ORCID:** [0009-0003-8219-3159](https://orcid.org/0009-0003-8219-3159)

## Abstract

This paper validates that distance-structured correlations in GNSS clocks exist in raw observations, not just processed products—eliminating the processing artifact hypothesis. Prior TEP analyses relied on precise orbit and clock products from global analysis centers, leaving open the possibility that observed signatures were artifacts of sophisticated processing chains. This paper addresses that concern by detecting distance-structured signatures in raw GNSS observations processed using Single Point Positioning (SPP) with broadcast ephemerides as the primary methodology, supplemented by precise ephemeris validation. Analysis of 539 globally distributed stations over 3 years (2022–2024, comprising 1.17 billion pair-samples across three independent filtering strategies) achieves 100% signal detection (72/72 metric combinations) with mean R² = 0.93, revealing directionally-structured correlations consistent with CODE's 25-year PPP findings (p < 10⁻¹⁵).

The primary finding is directional anisotropy: East-West correlations are 2–5% (MSC) to 22% (Phase Alignment) stronger than North-South at short distances (<500 km), with t-statistics up to 112 and Cohen's d up to 0.304. Month-by-month stratification shows stable polarity (E-W > N-S) at the 94–100% level across modes and metrics (worst case 34/36 months), consistent with a persistent underlying effect. A critical audit indicates this is not an artifact of distance distribution: E-W pairs are actually 13 km *longer* than N-S pairs (bias against signal), and robust distance-matching strengthens the ratio (1.033 → 1.041). At full distances, raw λ ratios can appear suppressed by distance-dependent biases; a geometry-corrected comparison yields ratios of 1.80–1.86, within 17% of CODE's benchmark (2.16).

Key validations include: (1) orbital velocity coupling detected at 3.2–5.4σ (best: r = −0.763), replicating CODE's 25-year finding (r = −0.888), with signal persisting under ionospheric removal (best ionofree: r = −0.416, 2.5σ); (2) position jitter and clock bias show similar orbital coupling (Δ ≈ 5%), consistent with spacetime—not just temporal—modulation; (3) CMB frame alignment at RA = 188°, Dec = −5° (20.0° from CMB dipole), matching CODE's benchmark (18.2°), with Solar Apex disfavored (86.5° separation); (4) geomagnetic stratification using real GFZ Kp data shows near-invariance at the primary threshold (Kp < 3 vs. Kp ≥ 3; median Δλ ≈ −1%); (5) year-specific planetary event modulation detected (2.8× above null, p < 0.001) with no consistent tidal GM/r² scaling, consistent with alignment-driven geometric coupling rather than tidal forcing.

This paper constitutes Paper 3 of the TEP-GNSS Research Series. Together with Paper 1 (multi-center validation) and Paper 2 (25-year temporal stability), these three complementary analyses—using different data sources, processing chains, and time periods—provide consistent evidence for planetary-scale, directionally-structured correlations in GNSS clock measurements.

## Key Findings

Raw RINEX processing confirms distance-structured correlations without reliance on precise orbit/clock products: 72/72 metric combinations detect the signal with mean R² = 0.93. Directional anisotropy persists at short ranges (East–West stronger by 2–22%), and a geometry-corrected comparison yields EW/NS ratios of 1.80–1.86, consistent with the CODE benchmark (2.16). Orbital velocity coupling is replicated at 3.2–5.4σ (best r = −0.763), and CMB frame alignment matches the long-span solution (RA = 188°, Dec = −5°, 20.0° from the dipole). These results exclude processing artifacts while preserving the same spatial and kinematic structure found in the multi-center and 25-year analyses.

---

## The TEP Research Program

| Paper | Repository | Title | DOI |
|-------|-----------|-------|-----|
| **Paper 0** | [TEP](https://github.com/matthewsmawfield/TEP) | Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed | [10.5281/zenodo.16921911](https://doi.org/10.5281/zenodo.16921911) |
| **Paper 1** | [TEP-GNSS](https://github.com/matthewsmawfield/TEP-GNSS) | Global Time Echoes: Distance-Structured Correlations in GNSS Clocks | [10.5281/zenodo.17127229](https://doi.org/10.5281/zenodo.17127229) |
| **Paper 2** | [TEP-GNSS-II](https://github.com/matthewsmawfield/TEP-GNSS-II) | Global Time Echoes: 25-Year Temporal Evolution of Distance-Structured Correlations in GNSS Clocks | [10.5281/zenodo.17517141](https://doi.org/10.5281/zenodo.17517141) |
| **Paper 3** | **TEP-GNSS-RINEX** (This repo) | Global Time Echoes: Raw RINEX Validation of Distance-Structured Correlations in GNSS Clocks | [10.5281/zenodo.17860166](https://doi.org/10.5281/zenodo.17860166) |
| **Paper 4** | [TEP-GL](https://github.com/matthewsmawfield/TEP-GL) | Temporal-Spatial Coupling in Gravitational Lensing: A Reinterpretation of Dark Matter Observations | [10.5281/zenodo.17982540](https://doi.org/10.5281/zenodo.17982540) |
| **Synthesis** | [TEP-GTE](https://github.com/matthewsmawfield/TEP-GTE) | Global Time Echoes: Empirical Validation of the Temporal Equivalence Principle | [10.5281/zenodo.18004832](https://doi.org/10.5281/zenodo.18004832) |
| **Paper 7** | [TEP-UCD](https://github.com/matthewsmawfield/TEP-UCD) | Universal Critical Density: Unifying Atomic, Galactic, and Compact Object Scales | [10.5281/zenodo.18064366](https://doi.org/10.5281/zenodo.18064366) |
| **Paper 8** | [TEP-RBH](https://github.com/matthewsmawfield/TEP-RBH) | The Soliton Wake: A Runaway Black Hole as a Gravitational Soliton | [10.5281/zenodo.18059251](https://doi.org/10.5281/zenodo.18059251) |
| **Paper 9** | [TEP-SLR](https://github.com/matthewsmawfield/TEP-SLR) | Global Time Echoes: Optical Validation of the Temporal Equivalence Principle via Satellite Laser Ranging | [10.5281/zenodo.18064582](https://doi.org/10.5281/zenodo.18064582) |
| **Paper 10** | [TEP-EXP](https://github.com/matthewsmawfield/TEP-EXP) | What Do Precision Tests of General Relativity Actually Measure? | [10.5281/zenodo.18109761](https://doi.org/10.5281/zenodo.18109761) |

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
- **Note**: RTKLIB is no longer bundled. Install independently and ensure `rnx2rtkp` is available in your PATH or at a configurable location.

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
| 2.2 | `step_2_2_anisotropy_analysis.py` | Directional (E-W vs N-S) anisotropy |
| 2.3 | `step_2_3_temporal_analysis.py` | Year-by-year & seasonal stability |
| 2.4 | `step_2_4_null_tests.py` | Solar/lunar/shuffle validation |
| 2.5 | `step_2_5_orbital_coupling.py` | Orbital velocity correlation |
| 2.6 | `step_2_6_planetary_events.py` | Planetary conjunction/opposition |
| 2.7 | `step_2_7_cmb_frame_analysis.py` | CMB frame grid search |

## Summary of Key Results and Findings

### Primary Results Table

| Metric | Value | Uncertainty | Significance |
|--------|-------|-------------|--------------|
| **Dataset Size** | 1.17 billion pair-samples | — | 539 stations |
| **Temporal Coverage** | 3 years | 2022–2024 | Raw RINEX data |
| **Signal Detection Rate** | 100% | 72/72 metric combinations | Mean R² = 0.93 |
| **Processing** | Single Point Positioning (SPP) | Broadcast ephemerides | No precise products |

### Correlation Length by Processing Mode

| Mode | λ (km) | R² | Interpretation |
|------|--------|-----|----------------|
| **Baseline (GPS L1)** | 727 | 0.971 | Ionosphere included |
| **Ionofree (L1+L2)** | 1,072 | 0.973 | Ionosphere removed |
| **Multi-GNSS** | 815 | 0.928 | All constellations |
| **CODE Cross-Validation** | 4,811 | — | Matches 25-year benchmark |

### Four Pillars of Validation

| Finding | Value | Significance | Comparison to CODE |
|---------|-------|--------------|-------------------|
| **Orbital Velocity Coupling** | r = −0.763 | 5.4σ | CODE: r = −0.888 ✓ |
| **CMB Frame Alignment** | RA=188°, Dec=−5° | 20.0° from dipole | CODE: 18.2° ✓ |
| **Spacetime Symmetry** | Δ ≈ 5% (pos/clock) | Identical coupling | Metric fluctuation |
| **Planetary Modulation** | 2.8× above null | p < 0.001 | No GM/r² scaling |

### Directional Anisotropy (Short Distances <500 km)

| Metric | E-W/N-S Ratio | t-statistic | Cohen's d |
|--------|---------------|-------------|-----------|
| **MSC Coherence** | 1.033 | — | — |
| **Phase Alignment** | 1.224 | up to 112 | up to 0.304 |
| **Temporal Stability** | 94–100% | 34–36/36 months | Persistent |
| **Geometry-Corrected** | 1.80–1.86 | — | Matches CODE (2.16) within 17% |

### Robustness Tests

| Test | Result | Interpretation |
|------|--------|----------------|
| **Geomagnetic (Kp) Independence** | Δλ ≈ −1% | Invariant under storm conditions |
| **Filter Independence** | CV < 15% | All three filters converge |
| **Solar Apex Rejection** | 86.5° separation | Disfavored vs CMB |

### Key Interpretation

This analysis eliminates the processing artifact hypothesis. By detecting the same signatures in raw RINEX observations processed with only broadcast ephemerides (no precise products), the findings demonstrate that distance-structured correlations exist in the fundamental data, not just sophisticated analysis center outputs. The replication of CODE's 25-year orbital velocity coupling (r = −0.763 vs r = −0.888) using completely independent methodology provides strong cross-validation. The identical coupling between position jitter and clock bias (Δ ≈ 5%) suggests spacetime metric fluctuations rather than clock-only effects—a key discriminant favoring TEP over instrumental explanations.

## File Structure

```
TEP-GNSS-RINEX/
├── scripts/
│   ├── steps/                      # Analysis pipeline (step_1_*, step_2_*)
│   └── utils/                      # Shared utilities (config, logger, etc.)
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
├── manuscript-rinex.md             # Auto-generated markdown
└── VERSION.json                    # Version metadata
```

## Requirements

- RTKLIB v2.4.3 (demo5 branch) installed independently; ensure `rnx2rtkp` binary is in your PATH or specify its location in the environment/config.
- Python packages: numpy, scipy, pandas, matplotlib, tqdm
- CDDIS authentication credentials (set `CDDIS_USER`/`CDDIS_PASS` or configure `.netrc`)

## Methodology

- **Processing**: RTKLIB SPP with broadcast ephemerides (no precise products)
- **Modes**: Baseline (GPS L1), Ionofree (L1+L2), Multi-GNSS (GPS+GLO+GAL+BDS)
- **Filters**: ALL_STATIONS (539), OPTIMAL_100 (100 balanced), DYNAMIC_50 (clock std < 50 ns)
- **Coherence**: Magnitude-weighted phase coherence, inverse-variance weighted fitting
- **Related**: [Paper 1 (Multi-Center)](https://matthewsmawfield.github.io/TEP-GNSS/) · [Paper 2 (25-Year CODE)](https://matthewsmawfield.github.io/TEP-GNSS-II/)

## Citation

```bibtex
@article{smawfield2025rinex,
  title={Global Time Echoes: Raw RINEX Validation of Distance-Structured Correlations in GNSS Clocks},
  author={Smawfield, Matthew Lukin},
  journal={Zenodo},
  year={2025},
  doi={10.5281/zenodo.17860166},
  note={Preprint v0.4 (Kathmandu)}
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
