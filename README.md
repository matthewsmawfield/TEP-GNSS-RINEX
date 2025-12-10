# TEP-GNSS-RINEX: Raw RINEX Validation of Distance-Structured Correlations

Analysis of Temporal Equivalence Principle (TEP) signatures using raw RINEX data and RTKLIB Single Point Positioning (SPP). This repository contains the analysis code and site content for:

> **Global Time Echoes: Raw RINEX Validation of Distance-Structured Correlations in GNSS Clocks**  
> TEP-GNSS Paper 3 | v0.2 (Kathmandu) | DOI: [10.5281/zenodo.17860166](https://doi.org/10.5281/zenodo.17860166)

**Live Site**: [matthewsmawfield.github.io/TEP-GNSS-RINEX](https://matthewsmawfield.github.io/TEP-GNSS-RINEX/)

As the capstone of the TEP-GNSS trilogy, this work addresses the critical question of whether previously observed correlations are artifacts of precise post-processing (PPP). By processing raw pseudorange measurements with only broadcast ephemerides, this analysis demonstrates that TEP signatures exist in the fundamental unprocessed GNSS observables.

**The TEP-GNSS Trilogy:**
1. [**Paper 1: Multi-Center Validation**](https://matthewsmawfield.github.io/TEP-GNSS/) (CODE/ESA/IGS consistency)
2. [**Paper 2: 25-Year Temporal Evolution**](https://matthewsmawfield.github.io/TEP-GNSS/code-longspan/) (Decadal stability)
3. **Paper 3: Raw RINEX Validation** (This repository - Artifact exclusion)

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

## Key Results

**Dataset**: 440 stations × 3 years (2022–2024) = 172 million pairs

| Mode | λ (km) | R² | Interpretation |
|------|--------|-----|----------------|
| **Baseline (GPS L1)** | 727 | 0.971 | Ionosphere included |
| **Ionofree (L1+L2)** | 1,072 | 0.973 | Ionosphere removed |
| **Multi-GNSS** | 815 | 0.928 | All constellations |

**Key Findings (4 Pillars)**:
- **Orbital Velocity Coupling**: r = −0.752 (5.3σ), independently replicating CODE's 25-year finding (r = −0.888).
- **CMB Frame Alignment**: Best-fit at RA = 157°, 19.3° from CMB dipole (matches CODE's 18.2°). Solar Apex statistically excluded.
- **Spacetime Symmetry**: Position jitter and clock bias exhibit identical orbital coupling (Δ = 2–13%), consistent with spacetime metric fluctuation.
- **Planetary Modulation**: Coherence modulation detected around 37 planetary events (2.1× null rate) with no mass scaling, ruling out tidal forcing.

**Validation**:
- **Directional Anisotropy**: E-W > N-S at short distances (MSC 1.02–1.05, Phase 1.15–1.28, p < 10⁻³⁰⁰).
- **Seasonal Breathing**: Anisotropy oscillates seasonally, peaking at 1.35-1.51 during Equinoxes (Apr/Sept) and troughing at Solstices. Global average (~0.95) masks this dynamic signal.
- **Geometry-Corrected Ratio**: 1.79–1.86 (matches CODE's 2.16 within 17%).
- **Geomagnetic Independence**: Signal invariant under storm conditions (Δλ < 4%).
- **Filter Independence**: σ² = 0 across three station selection methods.

**CODE Cross-Validation**: Ionofree λ = 4,811 km (2024) matches CODE's 25-year benchmark (4,201 ± 1,967 km)

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
- **Filters**: ALL_STATIONS (440), OPTIMAL_100 (50N+50S balanced), DYNAMIC_50 (quality-filtered)
- **Coherence**: Magnitude-weighted phase coherence, inverse-variance weighted fitting
- **Related**: [Paper 1 (Multi-Center)](https://matthewsmawfield.github.io/TEP-GNSS/) · [Paper 2 (25-Year CODE)](https://matthewsmawfield.github.io/TEP-GNSS/code-longspan/)

## License

This project is open source. See [LICENSE](LICENSE) for details.

## Acknowledgments

- NASA CDDIS for data distribution
- International GNSS Service (IGS) and contributing station operators
- Tomoji Takasu (RTKLIB)
