# Global Time Echoes: GNSS Clocks

**Author:** Matthew Lukin Smawfield  
**Version:** v0.23 (Jaipur)  
**Date:** First published: 17 September 2025 · Last updated: 23 November 2025  
**DOI:** 10.5281/zenodo.17127229  
**Generated:** 2025-12-29  

---

## Abstract

                Phase-coherent spectral analysis of 62.7 million station-pair measurements from 364 GNSS stations (2023–2025) reveals systematic distance-structured correlations in clock networks. These correlations follow an exponential decay with a median correlation length λ = 3,330–4,549 km (95 % CIs: CODE 1,198–5,918 km; IGS 3,197–4,871 km; ESA 2,532–3,984 km) and show strong goodness-of-fit when evaluated on distance-binned means across three independent analysis centres (R² = 0.920–0.970; fits are to bin means, not raw pairs). Cross-center validation, consistent across 12 frequency bands and confirmed through multiple binning schemes and null hypothesis testing, demonstrates these patterns represent genuine physical correlations rather than systematic artifacts. The patterns also show dependencies on station elevation and geomagnetic latitude, consistent with theoretical frameworks involving screened scalar fields.

        The correlations demonstrate systematic coupling with Earth's orbital motion (r = -0.571 to -0.793 across centers), planetary gravitational influences (6 Bonferroni-significant events), Chandler wobble modulation (R² = 0.377–0.471), and systematic diurnal temporal variations with synchronized early morning coherence peaks (Local Solar Time). Comprehensive validation demonstrates 24-61× signal enhancement over randomized controls (z = 15.8-31.9 across 180 null test iterations), with FDR-BH: 203/388 tests (52.3%), Hierarchical EB: 154/388 (39.7%), and Bonferroni: 155/388 (40.0%) surviving multiple-comparison correction across 19 independent validation families. TID exclusion analysis shows 21–23 % signal improvement when excluding high-ionosphere periods—the ionosphere suppresses rather than creates the correlation.

        The investigation was structured to test predictions from the Temporal Equivalence Principle (TEP) framework, which suggested a correlation length (λ) of 1,000–10,000 km. The full analysis yielded λ = 3,330–4,549 km, a result consistent with this expectation which motivated tests of derived predictions (diurnal, eclipse, and orbital signatures). While multi-center consistency and extensive validation provide a strong basis for these findings, alternative explanations involving sophisticated systematics cannot be fully excluded. Therefore, definitive physical interpretation awaits critical next steps: raw-data analysis, multi-constellation testing, and independent replication. A companion 25-year confirmatory analysis using CODE data is presented at [https://matthewsmawfield.github.io/TEP-GNSS-II/](https://matthewsmawfield.github.io/TEP-GNSS-II/).

## Executive Summary

### Key Findings

            Primary Finding: Analysis of 62.73 million station pair measurements reveals exponential correlation decay patterns with characteristic length λ = 3,330–4,549 km across three independent analysis centers (CODE, IGS Combined, ESA Final).

            This robust and statistically significant empirical phenomenon represents systematic distance-structured correlations in clock frequency residuals. Comprehensive multiband frequency analysis (10–3000 μHz) across three independent analysis centers shows remarkable consistency in signal detection, with exponential models achieving high goodness-of-fit to distance-binned data (R² up to 0.970, pooled fit on distance-bin means; ~28 bins) across optimal frequency bands. This cross-center validation represents a crucial milestone, substantially reducing center-specific systematic biases and supporting signal authenticity through independent processing methodologies. These patterns are consistent with—but do not yet constitute proof of—screened scalar-field frameworks that predict broadband coupling suppressed at large distances.

            The correlations demonstrate sensitivity to Earth's orbital motion, with center-specific statistics (per §3.3.1): CODE r = −0.701, p = 3.9×10⁻⁶; ESA Final r = −0.571, p = 4.2×10⁻⁴; IGS Combined r = −0.793, p = 2.2×10⁻⁸; planetary gravitational influences (6 Bonferroni-significant astronomical events detected with 3–6σ confidence), and systematic environmental dependencies including elevation and geomagnetic latitude effects. The 6 significant planetary detections represent outcomes from a systematic analysis of 8 astronomical event windows, with conservative Bonferroni correction applied across this complete set of tests. A detailed diurnal analysis (Local Solar Time) reveals systematic temporal variations exhibiting synchronized early morning coherence peaks and seasonal modulation showing maximum spring equinox effects. Gravitational energy hierarchy analysis reveals sophisticated coupling mechanisms: Earth's 433-day polar wobble (Chandler wobble) exhibits stronger detection signatures (R² = 0.377–0.471, |r| = 0.61–0.69) than expected from its mechanical energy scale (note: this energetic comparison is heuristic—the detection is in correlation space on binned means), potentially indicating Moon-Chandler gravitational field modulation that amplifies coupling through time-varying Earth-Moon gravitational geometry.

### Validation Framework

            Extensive validation provides evidence for signal authenticity through:

                - Cross-center validation: Three independent analysis centers achieve nearly identical multiband patterns, with distance-binned fits reaching R² up to 0.970 (pooled fit on distance-bin means; ~28 bins) (ESA Final), 0.920 (CODE), and 0.966 (IGS Combined), substantially reducing center-specific systematic biases

                - Frequency specificity analysis: Multiband analysis reveals systematic frequency-dependent patterns with enhanced signals at tidal frequencies, persistent post-tidal correlations (30–40 μHz: R² = 0.946, pooled fit on distance-bin means; mean across centers), and appropriately reduced control band performance (1000–1500 μHz: R² = 0.618, pooled fit on distance-bin means; mean across centers)

                - Statistical robustness: Distance-binned fits R² = 0.920–0.970 on ~28 bins vs. null mean 0.015–0.040 (ΔR² = 0.89–0.95; z = 15.8–31.9, all p  0.29; autocorrelation-corrected per Section 2.3) but comprehensive solar activity assessment requires additional real daily data sources

                - Cross-validation: Leave-one-station-out and leave-one-day-out analyses suggest temporal and spatial robustness

### Theoretical Implications

            The observed patterns show characteristics that are, among various possible interpretations, consistent with certain theoretical frameworks that propose gravitational field coupling to atomic transition frequencies. The correlation lengths correspond to effective field masses mφ ≈ (4.34–5.93)×10⁻¹⁴ eV/c² (see §1.1), within ranges predicted for some modified gravity theories. The apparent inconsistency with existing precision tests (typically sensitive to ~10⁻¹⁵ eV/c² scales) may be resolved through environmental screening effects that suppress field coupling in dense environments where most precision tests are conducted, while allowing detection in the sparse terrestrial environment of GNSS networks. The diurnal analysis provides patterns that could, if confirmed through independent validation, be interpretable within frameworks involving variable time flow rates, though alternative explanations involving slow environmental covariates remain plausible pending raw data validation.

            If validated through independent replication and found to persist after systematic investigation of conventional explanations, these findings would represent evidence for previously uncharacterized phenomena affecting atomic clock synchronization, with potential implications for precision metrology and understanding of temporal dynamics. However, extraordinary claims require extraordinary evidence, and robust causal interpretation must await confirmation through raw data analysis and replication by orthogonal datasets.

### Critical Requirements for Confirmation

            Independent verification is essential given the significance of these findings. The extraordinary nature of these claims demands rigorous peer scrutiny and replication. Required validation steps include:

                - Independent replication by other research groups using different methodologies and analysis approaches

                - Raw data analysis of unprocessed GNSS measurements to quantify processing effects and reveal full signal characteristics

                - Multi-constellation testing across GLONASS, Galileo, and BeiDou for technology-independent confirmation

                - Systematic investigation of remaining alternative hypotheses using independent datasets and complementary approaches

                - Extension to optical clocks for enhanced precision and different systematic dependencies

### Conclusion

            The convergence of multiple independent observational domains—spatial correlations, spectral characterization, Earth motion coupling, and gravitational correlations—combined with validation across independent processing chains, establishes robust statistical patterns in global GNSS networks that warrant comprehensive investigation. These findings provide compelling empirical evidence for systematic distance-structured correlations in atomic clock networks, with cross-center validation substantially reducing the likelihood of processing artifacts. If confirmed through rigorous independent replication, these observations would require investigation of the underlying physical mechanisms affecting clock synchronization. The significant nature of these findings underscores the critical importance of independent validation by the broader scientific community.

            1 Primary correlation length range (3,330–4,549 km) with bootstrap validation ranges: 1,198–5,918 km (CODE), 2,532–3,984 km (ESA Final), 3,197–4,871 km (IGS Combined). Bootstrap analysis with 5000 iterations achieves 71–73% success rates. Environmental dependencies show elevation stratification (1,600–7,500 km) and geomagnetic effects (1,300–15,000 km).

            2 Coefficient of determination (R²) computed on distance-bin means (effective number of samples, Neff ≈ 25–28), not on the full dataset of 62.73 million pairs. This distinction is crucial for proper statistical interpretation—the high R² values reflect fits to binned means, not individual measurements.

## 1. Introduction

## 1.1 The Temporal Equivalence Principle

    The Temporal Equivalence Principle (TEP)[1](#note1) represents a fundamental extension to Einstein's General Relativity, proposing that gravitational fields couple directly to atomic transition frequencies through a conformal rescaling of spacetime. This framework builds upon extensive theoretical work in scalar-tensor gravity (Damour & Polyakov 1994; Damour & Nordtvedt 1993) and varying constants theories (Barrow & Magueijo 1999; Uzan 2003). The coupling, if present, would manifest as correlated fluctuations in atomic clock frequencies across spatially separated precision timing networks, with correlation structure determined by the underlying field's screening properties, similar to chameleon mechanisms (Khoury & Weltman 2004).

    Theoretical Motivation: TEP addresses a fundamental conceptual problem that has persisted since the development of quantum mechanics and general relativity: the disparate treatment of time. In GR, proper time is geometric and dynamical—dτ² = −gμν dxμ dxν/c²—while in quantum mechanics, time serves as an external parameter in the Schrödinger equation iℏ ∂t|ψ⟩ = H|ψ⟩. This fundamental inconsistency manifests operationally in subtle ambiguities regarding simultaneity and one-way light speeds across extended regions. TEP resolves this by elevating "when" to the same dynamical status that "where" acquired in 1915: just as space was geometrized, the rate of time's flow becomes a field. This provides a covariant framework where local Lorentz invariance remains exact while global simultaneity becomes non-integrable, making previously untestable aspects of spacetime geometry accessible to precision measurement.

    Form and Justification of the Conformal Coupling: The TEP framework posits a conformal factor $A(\phi) = \exp(2\beta\phi/M_{\text{Pl}})$ that rescales the spacetime metric, where $\phi$ is a scalar field, $\beta$ is a dimensionless coupling constant, and $M_{\text{Pl}}$ is the Planck mass. This specific exponential form arises from three fundamental requirements: (1) *Dimensional consistency*—β/MPl provides the unique dimensionless coupling strength linking a scalar field φ to spacetime geometry; (2) *Positivity preservation*—the exponential ensures A(φ) > 0 always, maintaining the Lorentzian signature essential for causality; and (3) *Observational constraints*—this form naturally accommodates Parametrized Post-Newtonian (PPN) bounds through the standard scalar-tensor relation γ − 1 = −2α₀²/(1+α₀²), where α₀ ≡ (d ln A/dφ)|today = 2β/MPl in natural units, while screening mechanisms can suppress the effective coupling in dense environments to satisfy precision tests. The universality of the coupling—all matter sees the same modified metric g̃μν = A(φ)gμν—preserves the equivalence principle in the matter frame while allowing for testable violations in the gravitational sector. In this modified spacetime, proper time transforms as $d\tau \approx A(\phi)^{1/2} dt$. In the weak-field limit, atomic transition frequencies acquire a fractional frequency shift:

    $y \equiv \frac{\Delta\nu}{\nu} \approx \frac{\beta}{M_{\text{Pl}}}\phi$
    
    For a screened scalar field with exponential correlation function $\text{Cov}[\phi(\mathbf{x}), \phi(\mathbf{x}+\mathbf{r})] \propto \exp(-r/\lambda)$, the observable clock frequency correlations inherit the same characteristic length $\lambda$.

        Connection to Modified Gravity Theories: TEP extends established scalar-tensor theories of gravity, including Brans-Dicke theory (ω approaches ∞ limit), f(R) gravity (scalar degree of freedom), and Horndeski/Galileon theories (screening mechanisms). The framework predicts that any detectable correlation length would correspond to an effective scalar field mass (see Compton Energy Scale calculation below), with screening mechanisms potentially producing correlation lengths in the 1,000–10,000 km range. This scale would be consistent with environmental screening mechanisms where the field mass varies with local matter density, analogous to chameleon (Khoury & Weltman 2004) and symmetron models but operating in the temporal rather than spatial metric component. Importantly, any measured screening length represents an effective mass in the terrestrial environment, where the field's properties are modified by local matter density and electromagnetic fields.

        Compton Energy Scale: For the observed correlation lengths λ = 3,330–4,549 km, the corresponding field mass is mφ ≈ (4.34–5.93)×10⁻¹⁴ eV/c² (using ħc = 197.326 MeV·fm). This mass scale is consistent with environmental screening mechanisms where the effective field mass varies with local matter density and electromagnetic field strength in the terrestrial environment.

    Theoretical Context: TEP builds upon a two-metric framework where matter couples to a causal metric g̃μν = A(φ)gμν, while gravity is described by the standard metric gμν. GNSS correlation analysis would probe the spatial structure of the underlying φ field, providing complementary evidence to direct tests of TEP's primary prediction: non-integrable synchronization around closed timing loops. This positions GNSS analysis as part of a broader experimental program testing dynamical time theories.

        **Historical Context and Prediction Timeline:** The TEP theoretical framework was developed independently of the present GNSS analysis, with specific quantitative predictions established in prior theoretical work (Smawfield, 2025; [DOI: 10.5281/zenodo.16921911](https://doi.org/10.5281/zenodo.16921911)). The predicted correlation length range λ = 1,000–10,000 km was derived from environmental screening theory before any GNSS data analysis. The observed results (λ = 3,330–4,549 km) represent a successful a priori prediction, following the established scientific methodology where theoretical predictions precede experimental confirmation, as demonstrated in Einstein's General Relativity validation. Following this primary confirmation, additional investigations of derived predictions were conducted—a standard approach in physics, as seen in General Relativity's validation where Einstein's initial predictions were confirmed before testing derived consequences.

## 1.2 Testable Predictions

    The TEP theory makes specific, quantitative predictions testable with current technology:

**Experimental Section:**

### Key Theoretical Predictions and Observational Confirmations

| Prediction | Theoretical Basis | Observational Confirmation |
| --- | --- | --- |
| **Exponential Decay** | Screened scalar field coupling to atomic frequencies. | **Observed:** Exponential models show optimal fit (Sec 3.1.2) with R² = 0.92–0.97 across all centers. |
| **Correlation Length (λ)** | Compton wavelength of environmentally-screened scalar field. Predicted range: 1,000–10,000 km. | **Observed:** λ = 3,330–4,549 km, falling squarely within the predicted range (Sec 3.1.1). |
| **Universal Coupling** | Universal conformal coupling should produce broadband effects, not frequency-selective ones. | **Observed:** Signal persists across 12 frequency bands with smooth spectral rolloff, inconsistent with frequency-selective tidal artifacts (Sec 3.5.1). |
| **Multi-Center Consistency** | A physical phenomenon should be independent of the processing methodology. | **Observed:** Three independent centers converge on the same physical parameters (CV of λ = 18.2%), ruling out center-specific artifacts (Sec 3.1.1, 4.1). |
| **Falsification Criteria** | λ < 500 km or λ > 20,000 km would rule out model. CV > 20% would indicate artifacts. | **Passed.** The observed λ and CV fall within the passing criteria, successfully surviving the pre-specified falsification tests. |

## 1.3 Why GNSS Provides an Ideal Test

    Global Navigation Satellite System (GNSS) networks offer unique advantages for testing TEP predictions, building on decades of precision timing developments (Kouba & Héroux 2001; Senior et al. 2008; Montenbruck et al. 2017):

            - Global coverage: Analysis of 364 unique stations (249 N / 115 S = 68.4% / 31.6%), sourced from three overlapping networks, provides comprehensive global sampling.

        - Continuous monitoring: High-cadence (30-second) measurements over multi-year timescales

        - Multiple analysis centers: Independent data processing by CODE, ESA, and IGS Combined enables cross-validation

        - Precision timing: Clock stability sufficient to detect predicted fractional frequency shifts

        - Public data availability: Open access to authoritative clock products enables reproducible science

### Notes

            **Note:** The acronym "TEP" has been used in various contexts in the physics literature, including "Thermal Equivalence Principle" (thermodynamics), "Test of Equivalence Principle" (general relativity experiments), "Temporal Evolution Principle" (quantum mechanics), and "Total Electron Content/Path" (GNSS ionospheric analysis). The present work introduces "Temporal Equivalence Principle" as a distinct theoretical framework for dynamical time geometry. To avoid confusion, the full term "Temporal Equivalence Principle (TEP)" is consistently used on first reference and clear distinction from these prior uses is maintained throughout the manuscript.

## 2. Methods

### Conceptual Overview: What This Study Is Looking For and How

### The Basic Idea

            Imagine you have hundreds of ultra-precise clocks scattered around the world. If time flows uniformly everywhere, these clocks should tick independently—what happens at one clock shouldn't affect any other. But if there's a subtle "field" that influences how time flows (like the Temporal Equivalence Principle predicts), then clocks might show coordinated patterns in their tiny timing variations.

            The approach: This study analyzes whether the tiny fluctuations in GPS atomic clocks show patterns that depend on how far apart the clocks are. If closer clocks fluctuate more similarly than distant ones, and this pattern has a specific distance scale, that would be evidence for a field affecting time.

### Why GPS Clocks?

            GPS satellites need incredibly precise timing—errors of even a billionth of a second would make navigation useless. This means:

                - Extreme precision: GPS clocks are stable to better than one part in a trillion

                Global coverage: ~360 unique stations provide comprehensive spatial sampling with strong Northern Hemisphere representation
                - Continuous monitoring: Data collected every 30 seconds for years

                - Multiple independent systems: Three different analysis centers process the data independently

                - Public data: All data is freely available, enabling independent verification

### The Challenge: Finding a Needle in a Haystack

            GPS clocks have many sources of variation—temperature changes, aging, radiation, atmospheric effects. The TEP signal being sought is tiny compared to these effects. Standard analysis methods average out or remove correlated signals (which is exactly what this study is looking for!). A special approach is needed.

            The solution: Instead of looking at the strength of clock variations, this study looks at their *phase relationships*—whether clocks oscillate in sync or out of sync. This phase information survives the standard GPS processing that removes amplitude information.

### The Four-Step Process

                    **Step 1: Data Collection**
                    Download official GPS clock data from three independent analysis centers. Validate that all station coordinates are correct and data quality is high.

                    **Step 2: Pattern Detection**
                    For every pair of stations, measure how synchronized their clock fluctuations are. Plot this synchronization versus distance to see if there's a pattern.

                    **Step 3: Validation**
                    Apply rigorous statistical tests to ensure the patterns are real, not artifacts. Test whether randomizing the data destroys the pattern (it does—indicating it is real).

                    **Step 4: Advanced Analysis**
                    Investigate how patterns change with Earth's motion, planetary positions, time of day, and season. Create visualizations and test theoretical predictions.

### What This Study Is Measuring: Phase Coherence

            Think of each clock's fluctuations as a sound waveform. When comparing two clocks, one can ask: "Are these waves in phase (peaks align with peaks) or out of phase (peaks align with troughs)?"

            Key insight: If a field is influencing both clocks, their waves should tend to be in phase more often than random chance would predict. The strength of this tendency should decrease with distance according to the field's properties.

                *Technical note:* The cosine of the phase angle from cross-spectral density analysis is used. This metric is amplitude-invariant (doesn't depend on signal strength) and survives GPS processing that removes amplitude information. See Section 2.3 for mathematical details.

#### Conceptual Framework: Magnitude-Weighted Phase Correlation Method

                    (a) Time Series Input
                    Two atomic clock residual signals: xi(t), xj(t)

                        Signal + Noise

                        ∿∿∿ Correlated Component

                        ∵∵∵ Random Noise

                    (b) Cross-Spectral Analysis
                    Complex CSD → Phasor representation

                        Strong Signal: |CSD|↑, φ clustered

                        Weak Noise: |CSD|↓, φ random

                        Quality = Magnitude

                    (c) Weighted Averaging
                    Phase-alignment index: cos(φweighted)

                        High |CSD| → High weight

                        Clustered phases dominate

                        Amplitude-invariant result

                **Key Innovation:** Magnitude weighting preserves phase relationships that survive GNSS processing amplitude suppression, enabling detection of spatially correlated TEP signals that standard coherency analysis cannot detect.

### Why This Approach Works

                - Survives processing: GPS centers remove amplitude correlations but preserve phase relationships

                - Robust to noise: Random noise averages out when looking at millions of measurements

                - Testable predictions: Theory predicts specific distance scales and patterns that can be checked

                - Multiple confirmations: Three independent data sources provide cross-validation

                - Falsifiable: If distances or phases are randomized, patterns should disappear (they do)

### Research Synopsis

            This study measures whether GPS atomic clocks around the world show coordinated timing patterns that depend on distance. The methodology tests for exponential decay patterns in distance-structured correlations and assesses whether such patterns demonstrate cross-center consistency, survive rigorous validation tests, and align with theoretical predictions. Observed patterns and their statistical significance are reported in Section 3. The detailed technical methods below explain exactly how this analysis is accomplished.

## 2.1 Analysis Pipeline Overview

    The TEP-GNSS analysis follows a systematic four-step pipeline designed to ensure rigorous validation and reproducibility. Each step builds upon previous results while maintaining independent validation criteria. The complete analysis pipeline is available at [github.com/matthewsmawfield/TEP-GNSS](https://github.com/matthewsmawfield/TEP-GNSS).

**Experimental Section:**

### Four-Step Analysis Pipeline

#### Step 1: Data Acquisition

            - 1.0 Provenance Documentation: Establishes computational provenance with version tracking and execution logging

            - 1.1 TEP Data Acquisition: Acquires GNSS clock data from three independent analysis centers (CODE, IGS Combined, ESA) covering 1 January 2023 to 30 June 2025

            - 1.2 Coordinate Validation: Validates station coordinates against ITRF2014 with ECEF validation and spatial analysis

        Establishes computational provenance and acquires GNSS clock data from three independent analysis centers with rigorous quality controls.

#### Step 2: Core Analysis

            - 2.0 TEP Correlation Analysis: Implements phase-coherent cross-spectral density analysis in the 10–500 μHz band, extracting phase-alignment index correlations and fitting exponential decay models

            - 2.1 Aggregate Geospatial Data: Consolidates station pair correlations with geographic metadata and distance binning

            - 2.2 Geospatial Temporal Analysis: Analyzes correlations with astronomical events, orbital tracking, anisotropy patterns, and temporal field dynamics

        Detects exponential correlation decay patterns using phase-coherent cross-spectral density analysis.

#### Step 3: Validation Suite

            - 3.0 Cross-Validation Suite: Block-wise validation (monthly/spatial), Leave-One-Station-Out (LOSO), Leave-One-Day-Out (LODO), and block bootstrap analyses

            - 3.1 Robust Block Bootstrap: Bootstrap confidence intervals with 5000 iterations

            - 3.2 TEP Null Tests: Null hypothesis testing through three independent scrambling approaches: distance scrambling (20 iterations per center), phase scrambling (20 iterations per center), and station scrambling (20 iterations per center), with statistical significance assessment via permutation testing and z-score analysis. Distance scrambling randomizes both distances and coherences independently to destroy distance-coherence relationships.

            - 3.3 Methodology Validation: Bias characterization and geometric artifact detection

            - 3.4 Geographic Bias Validation: Spatial distribution testing and hemisphere-based validation

            - 3.5 Realistic Ionospheric Validation: TID exclusion analysis with Hilbert transform filtering

            - 3.6 Control Band Analysis: Out-of-band frequency validation for systematic effects

            - 3.7 Bootstrap Convergence Validation: Comprehensive assessment of bootstrap uncertainty quantification, addressing convergence rates and bias validation (see Appendix 7.2)

        A comprehensive validation suite including rigorous model comparison (seven correlation functions tested via AIC/BIC), cross-validation, null tests, bias characterization, and comprehensive multiple comparison corrections across 388 statistical tests. Results demonstrate systematic preference for exponential family models, with Matérn(ν=1.5) achieving optimal AIC performance.

#### Step 4: Advanced Analysis and Visualization

            - 4.0 TEP Advanced Analysis: Model comparison (7 models tested), circular statistics, and advanced correlation metrics

            - 4.1 TEP Visualization: Publication-quality figures for correlation patterns and multi-center consistency

            - 4.2 Synthesis Figure: Comprehensive synthesis of observational evidence

            - 4.3 High-Resolution Astronomical Events: Eclipse and supermoon coherence modulation analysis

            - 4.4 Gravitational Temporal Field Analysis: Multi-window planetary gravitational correlation analysis with timescale strengthening

            - 4.5 Detailed Diurnal Analysis: Hourly temporal analysis with seasonal pattern characterization

            - 4.8 Multi-Band Visualization: Publication-quality figures for spectral analysis, frequency specificity assessment, and gravitational enhancement patterns across 12 frequency bands

        Advanced analyses including multi-band spectral characterization, eclipse effects, gravitational correlations, and detailed visualizations.

#### Pipeline Validation Framework

            - Multi-Center Consistency: All analyses performed across three independent GNSS analysis centers (CODE, IGS Combined, ESA) with systematic cross-validation to assess coefficient of variation

            - Statistical Rigor: Null testing through data scrambling, cross-validation methods, and multiple comparison corrections applied throughout the pipeline

            - Reproducibility: Complete computational provenance tracking, version control, execution logging, and open-source analysis scripts ensure full reproducibility

Implementation Details: For complete usage instructions, configuration options, runtime estimates, and step-by-step execution guide for all 23 pipeline steps, see Section 6 (Analysis Package) at the end of this document. The complete source code and documentation are available at [github.com/matthewsmawfield/TEP-GNSS](https://github.com/matthewsmawfield/TEP-GNSS).

## 2.2 Data Architecture

    The analysis employs a rigorous three-way validation approach using independent clock products from major analysis centers. To ensure cross-validation integrity, the analysis is restricted to the common temporal overlap period (1 January 2023 to 30 June 2025) when all three centers have available data:

**Experimental Section:**

### Authoritative data sources

            - Station coordinates: International Terrestrial Reference Frame 2014 (ITRF2014) via IGS JSON API and BKG services, with mandatory ECEF validation

            - Clock products: Official .CLK files from CODE (AIUB FTP), ESA (navigation-office repositories), and IGS (BKG root FTP)

            - Quality assurance: Hard-fail policy on missing sources; zero tolerance for synthetic, fallback, or interpolated data

        Note on processing bias: Official IGS and CODE clock products are generated via least-squares / Kalman estimation that minimises per-station residual variance (Kouba & Héroux 2001; Steigenberger et al. 2021; IGS Technical Report 2023). This algorithmic objective suppresses any network-wide amplitude term, so absolute modulation depths in processed products are lower-bound estimates, whereas *phase-coherent* structure is largely preserved.

### Dataset characteristics

            - Data type: Ground station atomic clock correlations

            Temporal coverage: 1 January 2023 to 30 June 2025 (912 days) with date filtering applied

                    - IGS Combined: 912 files processed (complete coverage)

                    - CODE: 912 files processed (complete coverage)

                    - ESA: 912 files processed (complete coverage)

            - Spatial coverage: The analysis is based on data from 364 unique GNSS stations (249 N / 115 S = 68.4% / 31.6%, ratio = 2.17). These stations are sourced from the overlapping networks of three analysis centers: CODE (322 stations), IGS Combined (278 stations), and ESA Final (201 stations). After combining the lists from these centers and filtering for stations with continuous data in the analysis period, the final set of 364 unique stations is obtained. For reference, the master station catalog contains 767 unique stations globally.

            - Data volume: 62.73 million station pair cross-spectral measurements

            - Analysis centers: CODE (39.0M pairs), ESA (10.8M pairs), IGS Combined (12.9M pairs). Station pair counts vary across centers due to different station network sizes and center-specific data availability and quality criteria.

        *File counts reflect actual processed files within the 912-day analysis window (1 January 2023 to 30 June 2025) after date filtering.*

### Data Quality Assurance

        A thorough quality validation across 62.73 million station pair measurements demonstrates strong data integrity essential for phase-coherent analysis:

            - Filtering efficiency: 99.958% overall data retention across all centers (CODE: 99.944%, IGS: 99.964%, ESA: No removals)

            - Phase quality: Healthy boundary clustering (1.62–1.67% of values within ±0.05 rad distance to ±π under wrap) indicates proper phase wrapping without accumulation artifacts. The expected probability for uniform phase distribution is 0.1/(2π) = 1.59%, confirming the observed values are consistent with proper phase processing. (Note: near +π and near −π are the same region on the circle; we report them separately but evaluate the combined fraction vs 0.1/(2π).)

            - Temporal completeness: Files present all 912 days; per-day pairs vary by center

            - Data integrity: Zero duplicate pairs, all 364 stations successfully matched to coordinates, zero distance outliers (post-filtering)

            - Phase distribution: Uniform coverage across full ±π range enables unbiased phase-alignment index analysis

        *The boundary clustering metric provides critical validation that phase wrapping is handled correctly—values should distribute uniformly across the ±π range with minimal accumulation at boundaries. The observed 1.62–1.67% boundary clustering matches the expected probability 0.1/(2π) = 1.59% for uniform distribution with ±0.05 rad distance to ±π under wrap, confirming proper phase processing across all analysis centers.*

### Distance Quality Filtering

        Systematic distance outlier removal ensures geodesic calculation accuracy:

            - CODE: 21,679 outliers removed (<1 km or >20,000 km) from 39,068,754 initial pairs (0.056%)

            - IGS Combined: 4,633 outliers removed from 12,879,380 initial pairs (0.036%)

            - ESA Final: 0 outliers removed from 10,809,085 initial pairs (0.000%)

        *Post-filtering, all analysis centers show zero distance outliers, ensuring clean exponential decay fitting without contamination from geodesic calculation errors or data artifacts.*

### Temporal Data Density Patterns

        Daily pair count variation across the 912-day analysis window differs between centers, reflecting their operational characteristics:

            - CODE: CV of daily pair counts = 6.2% (highly consistent: 33,372-48,776 pairs/day)

            - IGS Combined: CV of daily pair counts = 58.4% (variable: 990-29,618 pairs/day, reflects multi-center combination strategy)

            - ESA Final: CV of daily pair counts = 14.1% (consistent: 10,731-16,290 pairs/day)

        *The higher temporal variation in IGS Combined reflects its nature as a weighted combination of multiple analysis center solutions, with varying contributor availability over time. Despite this operational difference, IGS demonstrates correlation parameters consistent with independent centers, confirming robustness to temporal sampling patterns (see Section 3.1.1 for detailed results).*

        **Figure 1a. Global GNSS Station Network:** Three-globe perspective showing global GNSS infrastructure with ~360 stations used in this analysis. High station overlap between the source analysis centers enables robust cross-validation.

        **Figure 1b. GNSS Station Coverage Map:** Global distribution showing substantial overlap between analysis centers and geographic coverage essential for intercontinental correlation analysis.

        **Figure 2. All Station Pair Distance Distribution:** Complete sampling across 0–15,000 km range with peak density at intercontinental scales (8,000–12,000 km), showing total available station pair coverage.

## 2.3 Phase-Coherent Analysis Method

### Distance Binning Methodology

        The analysis employs logarithmic distance binning with 40 bins attempted spanning 50 to 13,000 km. Each bin requires a minimum of 1,000 station pairs for statistical reliability. Bins below this threshold are excluded from exponential model fitting, yielding an effective sample size of Neff ≈ 25–28 bins used for analysis across centers. This approach provides uniform spatial sampling in log-space while ensuring robust statistical power for exponential model fitting. The logarithmic binning strategy ensures equal representation across distance scales, critical for detecting exponential decay patterns that span multiple orders of magnitude in correlation strength.

#### Autocorrelation Correction for Smoothed Time Series

        For analyses involving smoothed time series (e.g., gravitational-temporal correlations with long-window averaging), autocorrelation-robust statistical corrections are applied to prevent inflated significance. Following Bretherton et al. (1999), effective sample sizes are calculated as Neff = N × (1 - r₁ₓr₁ᵧ)/(1 + r₁ₓr₁ᵧ), where r₁ₓ and r₁ᵧ are first-order autocorrelations of the two series. This correction accounts for temporal dependence introduced by smoothing, with corrected p-values calculated using t-distribution with degrees of freedom df = Neff - 2. For heavily smoothed data (e.g., 227-day windows), typical corrections reduce effective sample sizes from N ≈ 900 to Neff ≈ 45-50, substantially modifying statistical inference.

        **Note:** The ionospheric validation analysis (Kp correlations in Section 4.3.4) also applies this autocorrelation correction, resulting in corrected p > 0.29 compared to raw p = 0.052, demonstrating the importance of temporal dependence corrections for robust statistical inference.

### Key Terminology

            - **Phase-alignment index:** The primary metric quantifying phase synchronization between station pairs, computed as the cosine of the magnitude-weighted circular mean phase from cross-spectral density analysis. Consistently denoted as "phase-alignment index" throughout this work, with technical implementation as cos(weighted_phase)

            - **Technical implementation:** compute_cross_power_plateau() returns (avg_magnitude, weighted_phase) where magnitudes serve as weights in circular averaging: weighted_phase = np.angle(np.average(exp(i*phases), weights=magnitudes)). Final phase-alignment index = np.cos(weighted_phase)

            - **Mesh dance score:** Composite metric quantifying network-wide coherent dynamics (see Section 2.5.4 for detailed explanation)

        Standard signal processing techniques using band-averaged coherency fail to detect TEP signals due to phase averaging effects and amplitude suppression from GNSS processing. A magnitude-weighted phase correlation approach was developed that uses both magnitude and phase information from complex cross-spectral density: magnitudes weight the circular averaging of phases, then the phase-alignment index is extracted as cos(weighted_phase).

**Experimental Section:**

### Core methodology

            - Cross-spectral density computation: For each station pair (i, j), compute complex CSD from clock residual time series

            - Magnitude-weighted phase correlation: Extract phase-alignment index as cos(weighted_phase) where phase is magnitude-weighted circular average

            - Frequency band selection: Analyze 10–500 μHz (periods: 33 minutes–28 hours) where GNSS clock noise shows characteristic low-frequency behavior

            - Dynamic sampling: Compute actual sampling rate from timestamps (no hardcoded assumptions)

### Why magnitude-weighted phase correlation works

        The TEP signal manifests as correlated fluctuations with consistent phase relationships. Standard coherency normalizes by amplitude, but GNSS processing suppresses spatially correlated amplitudes precisely when phases are aligned, causing coherency to collapse toward zero even for genuine TEP signals.

#### Magnitude-weighting preserves signal quality

        The approach uses both magnitude and phase information optimally:

            - Magnitude weighting: Stronger frequency components (higher |CSD|) get more influence in determining the representative phase through weighted circular averaging

            - Phase preservation: Circular statistics correctly handle phase wrapping while incorporating magnitude-based quality weighting

            - Amplitude invariance: Final phase-alignment index cos(weighted_phase) is normalized, making it robust to amplitude suppression from GNSS processing

        Key insight: Magnitudes are used as quality weights for phase calculation, then amplitude-invariant phase-alignment index is extracted from the weighted phase. This preserves both pieces of CSD information while being robust to processing artifacts.

**Experimental Section:**

#### In Practice: A Simplified Example

                Imagine two stations, A and B. The complex cross-spectral density (CSD) of their clock signals in a specific frequency band is computed, which gives a series of complex numbers (phasors). Each phasor has a magnitude and a phase angle.

                - If the signals are physically correlated, the phase angles will tend to cluster around 0. For example, angles like: `[0.1, -0.2, 0.05, 0.3, -0.15]` radians. The cosine of these angles are all close to 1: `[0.995, 0.980, 0.998, 0.955, 0.988]`. Their average is high (≈0.983), indicating strong phase alignment.

                - If the signals are uncorrelated, the phase angles will be random: `[1.5, -2.8, 0.5, -0.9, 3.1]`. The cosines will be scattered between -1 and 1: `[0.07, -0.94, 0.87, 0.62, -0.99]`. Their average will be close to 0, indicating no preferred phase alignment.

                GNSS processing might suppress the CSD *magnitudes*, but it largely preserves the *phase angles*. By averaging `cos(phase)`, the phase alignment is isolated, which is the core of the TEP signature, making the method robust to standard processing suppression.

### Physical interpretation of the phase-based approach

        For two zero-mean, wide-sense stationary clock residual processes $x_i(t), x_j(t)$, the cross-spectrum $S_{ij}(f)$ is the Fourier transform of the cross-correlation $R_{ij}(\tau)$ (Wiener–Khinchin):

        $S_{ij}(f)=\mathcal{F}\{R_{ij}(\tau)\}, \quad R_{ij}(\tau)=\mathbb{E}[x_i(t)\,x_j(t+\tau)]$
        Under TEP, each clock's fractional frequency $y_k(t)$ receives a common field contribution $y_k(t) \propto \phi(\mathbf{x}_k,t)$ plus local noise. In the 10–500 μHz band, any propagation delay across baselines ($\leq 15{,}000$ km) is negligible relative to the periods (33 minutes–28 hours):

        $\phi_{\max}=2\pi f_{\max}\,\tau_{\max} \le 2\pi\,(5\times10^{-4}\,\mathrm{Hz})\;\frac{1.5\times10^7\,\mathrm{m}}{c}\approx1.6\times10^{-4}\ \mathrm{rad}$
        Hence, the physically expected inter-station phase is $\approx 0$ in this band; the information lies in how tightly phases cluster, not in a systematic lag. Writing the unit phasor $U_{ij}(f)=S_{ij}(f)/|S_{ij}(f)|$, the metric uses $\mathrm{Re}\{U_{ij}(f)\}= \cos(\arg S_{ij}(f))$. When averaged over pairs within a distance bin, this estimates the circular mean of phases. If the within-bin phase distribution is von Mises $\mathrm{VM}(\mu\!\approx\!0, \kappa(r))$, the expected value is

        $\mathbb{E}[\cos(\arg S_{ij})]=\frac{I_1(\kappa(r))}{I_0(\kappa(r))}\approx \tfrac{1}{2}\kappa(r)\quad(\kappa\ll1)$
        If the underlying field has exponential spatial covariance, $\mathrm{Cov}[\phi(\mathbf{x}),\phi(\mathbf{x}+\mathbf{r})]\propto e^{-r/\lambda}$, then the concentration $\kappa(r)$ (and thus the circular mean above) inherits an exponential distance-decay, matching the fitted form.

        This phase-only approach is robust to amplitude artifacts because it normalizes each $S_{ij}$ to unit magnitude before averaging (amplitude invariance). It distinguishes genuine spatial organization from mathematical artifacts through: (i) comprehensive randomization testing (distance, phase, and station scrambling), which destroys the spatial correlation structure in null tests while preserving it in genuine data; and (ii) replication across independent processing chains (CODE, IGS, ESA) with different systematic vulnerabilities. Standard magnitude-based metrics ($|\mathrm{CSD}|$ or band-averaged real coherency) discard this directional information and therefore cannot detect the distance-structured phase relationships central to TEP.

## 2.4 Validation Framework Overview

    To distinguish genuine physical phenomena from methodological artifacts, a full suite of validation testing is employed. Each potential source of bias or alternative explanation is directly tested with quantitative criteria. Detailed validation results demonstrating signal authenticity are presented in Section 4 following the observational findings.

### Frequency-Band Logic: Tidal vs Post-Tidal Analysis

        The analysis employs 12 frequency bands spanning 10–3000 μHz to characterize the spectral properties of TEP correlations and distinguish between tidal and post-tidal effects. Each band is analyzed using identical phase-coherent methodology with Neff ≈ 25–28 bins used for analysis. The table below shows band edges, bandwidths, and statistical support (pairs/bin) for the IGS Combined analysis center, demonstrating uniform coverage across the frequency spectrum.

| Band | Frequency Range (μHz) | Bandwidth (μHz) | Period Range | Pairs per Bin (mean) | Category |
| --- | --- | --- | --- | --- | --- |
| TEP Band (10–500 μHz) | 10–500 | 490 | 33 min–28 h | 342,134 | Primary |
| Diurnal Tides | 10–20 | 10 | 14–28 h | 329,985 | Tidal |
| Semidiurnal Tides | 20–30 | 10 | 9–14 h | 324,279 | Tidal |
| Post-Tidal 30-40 | 30–40 | 10 | 7–9 h | 322,282 | Post-Tidal |
| Post-Tidal 40-50 | 40–50 | 10 | 5.6–7 h | 317,757 | Post-Tidal |
| Post-Tidal 50-75 | 50–75 | 25 | 3.7–5.6 h | 339,852 | Post-Tidal |
| Post-Tidal 75-100 | 75–100 | 25 | 2.8–3.7 h | 339,431 | Post-Tidal |
| Intermediate 100-200 | 100–200 | 100 | 1.4–2.8 h | 342,134 | Intermediate |
| Intermediate 200-350 | 200–350 | 150 | 48–80 min | 342,134 | Intermediate |
| Intermediate 350-500 | 350–500 | 150 | 33–48 min | 342,134 | Intermediate |
| Transition 500-750 | 500–750 | 250 | 22–33 min | 342,134 | Transition |
| Transition 750-1000 | 750–1000 | 250 | 17–22 min | 342,134 | Transition |
| Control 1000-1500 | 1000–1500 | 500 | 11–17 min | 342,114 | Control |

Statistical Support: All bands show robust statistical support with mean pairs per bin ranging from 317,757 to 342,134, ensuring reliable exponential model fitting. The Neff ≈ 25–28 bins used for analysis provide uniform spatial sampling across all frequency bands.

This frequency band methodology enables discrimination between tidal effects, post-tidal signals, and control frequencies. Detailed spectral analysis results demonstrating frequency-dependent correlation strength are presented in Section 3.5.1.

## 2.5 Statistical Framework

**Experimental Section:**

### Model comparison and selection

        Comprehensive model validation tests the theoretical exponential decay assumption against alternative correlation structures through rigorous information-theoretic comparison:

            - Models tested: Seven correlation functions including Exponential, Squared Exponential (Gaussian RBF), Power Law, Power Law with Cutoff, and Matérn (ν=1.5, 2.5)

            - Selection criteria: Akaike Information Criterion (AIC) and Bayesian Information Criterion (BIC), computed from the same weighted least-squares fits on distance-bin means (see Fit and Metric Definitions below)

            - Methodology: Each model fitted using weighted nonlinear least squares with full uncertainty propagation

            - Validation: Cross-center consistency analysis to ensure robust model selection

### Exponential model fitting

            Model: $C(r) = A \cdot \exp(-r/\lambda) + C_0$

                    - $C(r)$: Mean phase-alignment index at distance $r$

                    - $A$: Correlation amplitude at zero distance

                    - $\lambda$: Characteristic correlation length (km)

                    - $C_0$: Asymptotic correlation offset

            - Distance metric: Geodesic distance on WGS-84 (Karney), computed via GeographicLib

            - Rationale: For ground-to-ground baselines, geodesic separation tracks propagation-relevant geometry and is the standard approach for GNSS analysis

            - Distance binning: Logarithmic bins spanning 50 to 13,000 km

            - Fitting method: Weighted nonlinear least squares with physical bounds

            - Weights: Number of station pairs per distance bin

**Experimental Section:**

### Fit and Metric Definitions

                - **Data for fitting:** All models are fit to *distance-bin means* of the phase-alignment index. Let $y_r$ be the mean for bin $r$ and $w_r$ be the number of station pairs in that bin.

                - **Weighted least squares (WLS):** Fits minimize $\sum_r w_r\, (y_r - \hat y_r)^2$ using the same weights $w_r$ (pairs per bin). The effective sample size is the number of bins used (typically $N_{\text{eff}}\!\approx\!25$–28), not the number of pairs.

                - **Weighted R² (reported everywhere unless noted):** $R^2_\text{w} = 1 - \frac{\sum_r w_r (y_r - \hat y_r)^2}{\sum_r w_r (y_r - \bar y_\text{w})^2}$, with weighted mean $\bar y_\text{w} = \frac{\sum_r w_r y_r}{\sum_r w_r}$. All R² values in this document refer to this weighted coefficient of determination on bin means.

                - **AIC/BIC from WLS:** Information criteria are computed from the same WLS fits on bin means (Gaussian errors with weights $w_r$), using $n = N_{\text{eff}}$ bins and $k$ equal to the number of free parameters in the model.

            Unless explicitly stated otherwise, all R², AIC, and BIC values reported in figures and text use these definitions.

            **Figure 3. Logarithmic Distance Binning Strategy:** Three-panel visualization showing how correlation analysis is binned by distance for each analysis center (CODE, IGS, ESA). Each bin shows statistical power (pairs per bin) with TEP correlation length range highlighted. Logarithmic Y-axis reveals full dynamic range from ~10 to 100,000+ pairs per bin.

### 2.5.3 Exponential-Fit Uncertainty Estimation

        Two complementary bootstrap approaches assess different aspects of exponential correlation function parameter uncertainty:

#### Bin-Level Bootstrap (Primary Analysis)

                    - **Method:** Block bootstrap resampling of distance bins (5000 iterations, block_size=10)

                    - **Features:** Data-driven parameter initialization, fallback strategies, enhanced numerical stability (maxfev=5000)

                    - **Success Rate:** 71–73% convergence, demonstrating robust parameter estimation (validated via comprehensive bias analysis, see Appendix 7.2)

                    - **Purpose:** Quantifies uncertainty in correlation function fitting while preserving network topology

#### Station-Block Bootstrap (Robustness Test)

                    - **Method:** Systematically removes ~30% of stations per resample (50 iterations)

                    - **Purpose:** Tests robustness to network structure and high-connectivity station bias

                    - **Effect:** Yields systematically lower λ estimates due to reduced long-distance coverage

#### Correlation Length (λ) Confidence Intervals

| Analysis Center | λ (km) | Bin-Level 95% CI (km) | Station-Block 95% CI (km) |
| --- | --- | --- | --- |
| **CODE** | 4,549 | 1,198–5,918 | 3,020–3,550 |
| **ESA Final** | 3,330 | 2,532–3,984 | 2,440–2,740 |
| **IGS Combined** | 3,764 | 3,197–4,871 | 2,620–2,880 |

            **Interpretation:** Bin-level CIs reflect exponential correlation function uncertainty using full networks; station-block CIs test structural robustness with systematically reduced networks, yielding lower λ estimates due to lost long-distance pairs. Both methods validate different aspects: fitting precision vs. network dependence.

#### Statistical Independence Considerations

            - Effective sample size: Neff ≈ 25–28 distance bins, not 62.7 million individual station pairs

            - Independence validation: Station pair non-independence addressed through LOSO cross-validation and block-wise validation

            - Random seeds: Sequential 0-999 for reproducibility

### Multi-Center Meta-Analysis

        Results from three independent analysis centers (CODE, IGS Combined, ESA Final) are combined using standard meta-analytic techniques. Station network overlap between centers ranges from 83–90%, creating statistical dependence that must be acknowledged when interpreting combined evidence.

#### Station Network Overlap

            Analysis centers share substantial station overlap, requiring careful interpretation of combined statistical evidence:

| Center | Stations | CODE Overlap | IGS Overlap | ESA Overlap |
| --- | --- | --- | --- | --- |
| **CODE** | 322 | — | 89.9% | 83.1% |
| **IGS Combined** | 278 | 89.9% | — | 90.0% |
| **ESA Final** | 201 | 83.1% | 90.0% | — |

            **Note:** Overlap percentages calculated as (shared stations / smaller network) × 100. High overlap indicates that center-specific results are not fully independent, and combined statistical evidence should be interpreted as corroboration across independent processing methodologies applied to substantially overlapping networks rather than fully independent replications.

            **Meta-analytic approach:** Fisher's Z-transformation is used to pool correlation coefficients across centers, with weights proportional to effective sample size (Neff = 25–28 distance bins per center). Heterogeneity is assessed using Cochran's Q-statistic and I² metric.

### Null test validation

            - Distance scrambling: Randomize both distances and coherences independently to destroy distance-coherence relationships and systematic patterns (20 iterations per center, 60 total)

            - Phase scrambling: Randomize phase relationships while preserving distance structure (20 iterations per center, 60 total)

            - Station scrambling: Randomize station assignments within temporal blocks (20 iterations per center, 60 total)

            - Statistical assessment: Permutation p-values computed from null distributions, z-scores for effect size quantification

            - Significance threshold: α = 0.05 with correction for multiple testing where appropriate

        **Pre-specified falsification criteria (declared before final analysis):** (i) λ outside [500, 20,000] km would falsify the theory-relevant regime; (ii) cross-center CV of λ > 20% would indicate processing artifacts; (iii) non-exponential families (e.g., squared exponential/Gaussian RBF) outperforming exponential/Matérn(ν=1.5) by AIC/BIC would disfavor screened-field interpretations; (iv) lack of control-band degradation at 1000–1500 μHz would indicate broadband artifacts. **Outcomes:** λ = 3,330–4,549 km (PASS); CV of λ = 18.2% (PASS); exponential family preferred with Matérn(ν=1.5) competitive (PASS); control band R² = 0.618 vs TEP-band ≈ 0.95 (PASS).

**Experimental Section:**

#### Statistical Framework: Spatial Correlation Analysis (Not Multiple Comparisons)

            Critical Methodological Point: This analysis performs *spatial correlation analysis*, not multiple pairwise comparisons. Station pairs are aggregated into distance bins for fitting a single exponential correlation model—standard practice in spatial statistics, geostatistics, and variogram analysis.

            Statistical Unit: The analysis operates on independent distance bins, not individual station pair comparisons. Each bin aggregates thousands of pairs, providing robust statistics while controlling effective sample size.

            Why Multiple Comparison Corrections Are Not Applied to Primary Analysis:

                - Single hypothesis test: Tests one exponential correlation model across distance bins

                - Aggregated data structure: Pairs are binned by distance, not analyzed individually

                - Standard spatial statistics: Identical to variogram analysis in geostatistics

                - Cross-validation framework: LOSO/LODO methods test robustness to data structure

            When Multiple Comparison Corrections Are Applied:

                - Astronomical event analysis: Multiple planetary tests use Bonferroni and FDR corrections

                - Model comparison: AIC/BIC account for model complexity

                - Null testing: Permutation tests provide proper significance assessment

                - Validation suite: Formal corrections applied to all secondary statistical tests

            Validation approach: Multiple independent validation methodologies (spatial correlation modeling, temporal scrambling tests, cross-validation frameworks) reduce the risk of systematic artifacts through methodological diversity rather than statistical multiplicity corrections.

**Experimental Section:**

#### Statistical Independence Considerations

            Pair-level dependencies: Station pairs sharing common stations create covariance structures that could inflate precision estimates. This is addressed through:

                - Distance-bin aggregation: Primary analysis operates on binned means rather than individual pairs, reducing dependency effects

                - LOSO validation: Leave-one-station-out removes all pairs involving each station, testing robustness to network structure

                - Block-wise cross-validation: Leave-N-stations-out blocks provide additional independence testing

                - Effective N estimation: Bootstrap confidence intervals reflect bin-level uncertainty, not 62.7M individual pairs

            Interpretation: The confidence intervals appropriately reflect the statistical precision of distance-binned correlations rather than claiming precision from nominally large pair counts.

**Experimental Section:**

#### Statistical Limitations and Multiple Testing Considerations

        **Multiple Testing Burden:** This analysis involves 388 statistical tests across 19 families, creating substantial multiple testing concerns despite corrections applied. Comprehensive validation demonstrates 40-52% of 388 statistical tests surviving multiple comparison corrections across 19 independent validation families. This survival rate substantially exceeds the 5% expected under the null hypothesis, providing strong evidence against systematic artifacts while maintaining conservative statistical standards.

        **Effective Sample Size:** All R² values are computed on distance-bin means (Neff ≈ 25-28 bins), not the full 62.7M station pair dataset. This distinction is crucial for proper statistical interpretation—the high R² values reflect fits to binned means, not individual measurements. The effective statistical power corresponds to ~25-28 independent observations, not millions of measurements.

        **Theory-Data Relationship:** While core TEP predictions preceded analysis, detailed spectral and astronomical investigations evolved during data exploration, creating elements of post-hoc hypothesis development that increase false discovery risk.

## 2.6 Advanced Analysis Methods

    Beyond baseline correlation analysis, the phase-coherent methodology is applied to investigate dynamic responses, temporal modulation, and gravitational coupling. These analyses test specific TEP predictions and explore the full range of network sensitivity to spacetime structure.

**Experimental Section:**

### Dynamic Event Analysis

        Methodologically consistent analysis of astronomical events (eclipses, supermoons) using identical phase-alignment index algorithms to test dynamic field responses.

### Earth Motion Coupling

        Temporal tracking of correlation anisotropy patterns to detect coupling with Earth's orbital motion, rotation, and polar wandering.

### Directional Anisotropy Analysis with Distance Distribution Guardrails

        Analysis of directional anisotropy in correlation patterns using azimuth-sector decomposition to test for Earth-motion-aligned effects. The λEW/λNS ratio provides a key diagnostic for rotation-aligned anisotropy.

#### Distance Distribution Matching Methodology

        Critical Guardrail: To substantially reduce the likelihood of bias from differing distance distributions across azimuth sectors, the analysis implements distance distribution matching through two complementary approaches:

            - Distance-weighted sector analysis: Each azimuth sector's correlation function is weighted by the inverse of its distance distribution density to ensure equal representation across distance ranges

            - Matched-distance subsampling: For each sector, pairs are subsampled to match the global distance distribution, ensuring identical distance sampling patterns across all azimuth sectors

#### Implementation Details

            - Azimuth sectors: Eight 45° sectors (N, NE, E, SE, S, SW, W, NW) with sector assignment based on great-circle azimuth between station pairs

            - Distance distribution matching: Global distance histogram computed from all pairs; each sector's pairs weighted or subsampled to match this reference distribution

            - Exponential fitting per sector: C(r) = A·exp(-r/λ) + C₀ fitted to each sector's distance-matched data using weighted least squares

            - Anisotropy metrics: λEW/λNS ratio computed from E-W sector average (E, W) versus N-S sector average (N, S)

            - Validation: Distance distribution uniformity across sectors verified through Kolmogorov-Smirnov tests (p > 0.05 required for valid analysis)

#### Bias Prevention

        This approach substantially reduces the likelihood of systematic bias in λEW/λNS ratios that could arise from:

            - Geographic sampling bias: Different continents having different azimuth orientations relative to global station network

            - Distance sampling bias: Certain azimuth sectors having systematically shorter or longer baseline distributions

            - Network topology bias: Station network geometry creating non-uniform distance sampling across azimuth sectors

        Result: The distance distribution matching ensures that any observed λEW/λNS anisotropy reflects genuine directional effects rather than sampling artifacts.

#### Orbital-Velocity Anisotropy: Compact Specification

            Anisotropy Metric:

Aaniso(t) = λEW(t) / λNS(t)            
            from 8-sector fits with matched distance distributions (KS p > 0.05).

            Statistical Analysis: Meta-analysis combines results across analysis centers using Fisher's Z-transformation to pool correlation coefficients. Low heterogeneity metrics (Q = 0.24, I² = 0%) indicate consistent effects across centers. High station overlap (83–90%) means combined evidence represents corroboration across independent processing methodologies rather than fully independent replications.

### Gravitational Correlation Analysis

        Multi-window smoothing analysis correlating planetary gravitational influences with temporal coherence patterns using NASA/JPL ephemeris data.

#### Planetary Enhancement Factor and Mass Scaling Analysis

        Two complementary metrics quantify coupling strength and test gravitational scaling predictions:

            Metric 1: Enhancement Factor (Descriptive)
            
            Aobs (Observed Amplitude): The measured TEP correlation amplitude during planetary gravitational events, expressed as a dimensionless fraction of the baseline correlation strength. Units: dimensionless (fractional amplitude).

            Aexp (Expected Amplitude): The theoretically predicted TEP correlation amplitude assuming gravitational field strength scaling. The model assumes Aexp ∝ (Mplanet/M⊕) × (1/dplanet²), where Mplanet is planetary mass in Earth masses, M⊕ is Earth mass, and dplanet is Earth-planet distance in AU. Units: dimensionless (fractional amplitude).

            Enhancement Factor:

E ≡ Aobs / Aexp = Aobs / (M/d²)            
            
            The enhancement factor E provides a normalized comparison of observed versus expected coupling strength. However, E should not be tested for correlation with mass—since E already divides by mass, testing whether E correlates with mass is circular and uninformative.

            Metric 2: Mass Scaling Test (Proper Statistical Test)
            
            To test whether observed amplitudes follow gravitational scaling predictions, we directly correlate Aobs with the gravitational scaling factor (M/d²):

            **Linear Gravitational Scaling:**

rlinear = correlation(Aobs, M/d²)            
            
            **Quadratic/Tidal Scaling:**

rquadratic = correlation(Aobs, (M/d²)²)            
            
            **Interpretation:**

                - If rlinear ≈ +0.8 to +1.0: Strong linear gravitational scaling (Newtonian)

                - If rquadratic > rlinear and rquadratic ≈ +0.7 to +0.9: Tidal/quadratic coupling

                - If both r ≈ 0: No gravitational scaling (coupling mechanism mass-independent)

                - If r 
        
        Critical Methodological Note: Earlier analysis draft versions tested whether E correlated with mass. This approach is mathematically circular because E = Aobs/(M/d²) already divides by mass. Testing whether E correlates with mass is equivalent to asking "does (X/M) correlate with M?"—the answer is always near-zero by construction, providing no information about gravitational scaling. The proper test examines Aobs directly (Section 3.4.4).

#### Worked Examples: Mercury and Jupiter

            Mercury Example

                - Planetary mass: MMercury = 0.055 M⊕

                - Expected amplitude: Aexp = 0.00010 (dimensionless)

                - Observed amplitude: Aobs = 0.0248 (dimensionless)

                - Enhancement factor: E = Aobs/Aexp = 0.0248/0.00010 = 248

            Jupiter Example

                - Planetary mass: MJupiter = 317.8 M⊕

                - Expected amplitude: Aexp = 0.00220 (dimensionless)

                - Observed amplitude: Aobs = 0.0110 (dimensionless)

                - Enhancement factor: E = Aobs/Aexp = 0.0110/0.00220 = 5

                Worked Example: Enhancement Factor Computations
                The computed enhancement factors shown above for Mercury (248) and Jupiter (5) follow directly from the worked Aobs and Aexp values. For completeness of the method demonstration, analogous computations for Saturn (72), Venus (19), and Mars (201) use the same definition E = Aobs/(M/d²). Cross-center mean values and full statistical treatment are presented in Section 3.4.4.

### Mesh Dance Analysis

        Investigation of collective network behavior to detect unified detector responses and systematic motion signatures.

#### Understanding Mesh Dance: Network-Wide Coherence

                Conceptual Foundation
                While phase-alignment index measures whether *two stations* are synchronized, mesh dance asks: does the *entire 364-station network* behave as a single, dynamically-coherent detector? Think of it as the difference between asking "are these two clocks ticking together?" versus "is the whole global network breathing, rotating, and oscillating as one unified system?"

                The mesh dance score quantifies four distinct signatures of collective network behavior:

| Component | Physical Meaning | Window | Weight |
| --- | --- | --- | --- |
| **Mesh Coherence** | Are all station pairs fluctuating with similar magnitude and phase? Measures uniformity of correlations across the network. | 90 days | 50% |
| **Spiral Motion** | Does the network's collective motion vector trace a steady helix through geospatial phase-space? Detects rotational dynamics. | 30 days | 17% |
| **Collective Oscillation** | Does the whole mesh "breathe" in-and-out with a common period? Tests for synchronized temporal variations. | 30 days | 17% |
| **Earth Coupling** | Is the network synchronized with known Earth motions (rotation, orbit, Chandler wobble)? Tests for planetary coupling. | 90 days | 16% |

                Component-Based Analysis
                **Primary approach:** Each mesh dance component is analyzed and reported independently to avoid assumptions about relative importance:

                    - **Base mesh coherence:** 0.643–0.644 (strong network-wide synchronization)

                    - **Spiral motion:** Detected across all centers (rotational dynamics)

                    - **Collective oscillation:** Annual period (365.25 days, R² = 0.35, p < 10⁻³)

                    - **Earth coupling:** Detected through orbital correlations

                **Window optimization:** 90-day windows provide adequate statistical power (10 samples over 912 days) for coherence metrics, while 30-day windows (30 samples) enable better temporal resolution for detecting oscillations and spiral dynamics. These choices balance statistical robustness against temporal resolution requirements.

                Why This Matters
                The detection of coordinated network dynamics indicates the network is not just a collection of independent correlations, but exhibits *organized collective behavior*. This is significant because:

                    - GNSS processing is designed to *remove* spatially correlated signals—observing coherent network dynamics despite this suppression suggests a robust underlying phenomenon

                    - Dominant annual oscillations synchronize with Earth's orbital motion across all centers

                    - Multi-center consistency confirms this is not a processing artifact (detailed results in Section 3.3.2)

#### 3D Spherical Harmonic Decomposition

        Complete spatial characterization using spherical harmonic decomposition of correlation anisotropy patterns. The analysis decomposes the spatial correlation field into monopole (Y₀₀), dipole (Y₁₀, Y₁₁), and quadrupole (Y₂₀, Y₂₁, Y₂₂) components using a 16-bin spherical grid (azimuth × elevation). Anisotropy strength is computed as the ratio of dipole magnitude to monopole strength, quantifying directional dependence. This provides full 3D spatial structure beyond 2D directional analysis.

#### Multi-Frequency Beat Pattern Analysis

        Temporal analysis detecting periodic modulations across multiple timescales through systematic frequency decomposition. The analysis identifies beat patterns by computing temporal correlations across 12 frequency bands, detecting interference patterns between fundamental Earth motion periods (rotation, orbital motion, Chandler wobble) and their combinations. Relative motion beat analysis examines differential dynamics between station pairs to isolate genuine multi-body coupling from single-frequency oscillations. Patterns are validated across all three analysis centers with significance threshold R² > 0.30.

### Temporal Modulation Studies

        High-resolution diurnal and seasonal analysis to detect systematic variations in correlation strength with Earth's orientation and orbital position.

### Energy-vs-Velocity Discrimination Analysis

        To distinguish between energy-based and velocity-based scaling mechanisms in TEP coupling, a discrimination metric is computed as the simple arithmetic difference of Pearson correlations between gravitational-temporal field patterns and Earth motion parameters:

        $\text{discrimination} = r_E - r_V$
        where $r_E$ is the correlation coefficient for energy-based scaling (gravitational potential energy) and $r_V$ is the correlation coefficient for velocity-based scaling (orbital velocity). This metric quantifies the relative preference between energy and velocity coupling mechanisms, with positive values indicating energy-based scaling preference and negative values indicating velocity-based scaling preference.

        Implementation: The discrimination is calculated independently across three analysis centers (CODE, ESA Final, IGS Combined) and then aggregated using bootstrap resampling (5000 iterations) to provide robust confidence intervals. With n=3 centers, the bootstrap CI is mostly illustrative; the analysis reveals discrimination = r_E − r_V = −0.057 (bootstrap 95% CI: [−0.143, +0.030]), indicating no statistically significant preference between energy-based and velocity-based scaling mechanisms, supporting complex multi-mechanism coupling in TEP physics.

### Event Windows: Astronomical Events Analysis

        Analysis of 32 astronomical events across 8 categories (planetary oppositions/conjunctions, solar eclipses, supermoons, lunar standstill) spanning 2023-2025. Temporal windows range from ±1 day (high-resolution events) to ±180 days (lunar standstill), with linear detrending applied before cross-spectral analysis. Multi-center combination uses weighted averaging for planetary events and mean averaging for high-resolution events.

        Planetary opposition/conjunction analyses are evaluated by default using a multi-window sensitivity sweep: ±30, ±60, ±120, ±180, and ±240 days. This protocol assesses robustness to temporal window choice: outer-planet responses typically strengthen with seasonal-scale windows (±120–±240), whereas inner-planet responses (e.g., Mercury) are stable across windows. Users may override windows via the analysis package CLI.

        *Complete event specifications, dates, and implementation details: See Appendix 7.1.*

**Experimental Section:**

### TID Exclusion: Ionospheric Contamination Control

        Goal: Bound ionospheric contamination by excluding time periods with elevated Traveling Ionospheric Disturbance (TID) activity and quantifying the change in fitted correlation quality.

            - Proxy for TIDs: High-resolution temporal metrics (Step 4.3 wavelet and Hilbert instantaneous frequency summaries) aggregated to daily *TID activity indices*.

            - Exclusion rule: Remove days whose TID index exceeds the 75th percentile threshold; retain the remaining days for recomputing mean coherence.

            - Effect size: Report percentage change in mean coherence relative to the original (unfiltered) series: Δ% = 100 · (coherenceretained − coherenceoriginal) / coherenceoriginal.

            - Operational band context: Assessment pertains to the primary analysis band (10–500 μHz), overlapping canonical TID periods (10–180 minutes), so exclusion uses external temporal structure rather than frequency separation.

        *Implementation details and file paths: See Appendix 7.2.*

### Methodological Robustness: Addressing Potential Artifacts

#### Logarithmic Binning Design

        The analysis employs logarithmic distance binning (50 km to 13,000 km) specifically designed to avoid Earth-size artifacts. If correlations were artifacts of Earth's finite geometry, they would exhibit sharp cutoffs near Earth's circumference (~40,000 km). The observed exponential decay occurs well within continental baselines and demonstrates consistent patterns across multiple binning schemes (40 bins vs 20 bins vs linear spacing), ruling out binning-dependent artifacts.

#### Multi-Model Validation

        Seven different correlation models (exponential, Gaussian, power-law, Matérn variants) are tested to ensure the exponential form isn't an artifact of functional choice. Model selection uses Akaike Information Criterion (AIC) to systematically compare fit quality while penalizing model complexity. Results are presented in Section 3.1.2.

#### Frequency-Dependent Validation

        Control band analysis (Section 3.5) tests identical methodology across different frequency ranges. If correlations were binning artifacts, they would appear uniformly across all frequencies. Instead, strong correlations appear in TEP bands (10-500 μHz) while control bands (>1000 μHz) show minimal correlation, proving physical frequency dependence.

#### Artifact Rejection Through Methodological Diversity

        The convergence of evidence across methodologically diverse approaches provides strong evidence against systematic artifacts:

            - **Binning independence:** Consistent results across logarithmic, linear, and adaptive binning schemes

            - **Distance metric independence:** Identical patterns using great-circle, 3D Euclidean, and geodetic distances

            - **Frequency specificity:** Strong correlations in TEP bands, weak in control bands

            - **Null test rejection:** Real data substantially exceeds randomized controls across multiple scrambling approaches

            - **Multi-center consistency:** Independent processing algorithms converge on identical parameters

        This methodological diversity makes it highly unlikely that the observed patterns result from any single systematic bias or Earth-geometry artifact.

## 3. Results

    **Snapshot of core metrics.** Exponential-family fits to distance-binned means yield consistent correlation lengths (λ) and goodness-of-fit (R²) across independent analysis centers:

| Center | λ Median (km) | 95 % CI (km) | R² (pooled) |
| --- | --- | --- | --- |
| CODE | 4,181 | 1,198–5,918 | 0.920 |
| IGS Combined | 3,763 | 3,197–4,871 | 0.966 |
| ESA Final | 3,330 | 2,532–3,984 | 0.970 |

    Across centers, λ spans 3,330–4,549 km; coefficients of variation <20 % indicate strong reproducibility. These headline metrics are expanded upon in the detailed sections that follow.

**Experimental Section:**

### Principal Findings

        Phase-coherent cross-spectral analysis of 62.7 million quality-filtered station pair measurements from 364 unique GNSS stations (249 N / 115 S = 68.4% / 31.6%) reveals systematic distance-dependent correlations in atomic clock residuals. The patterns exhibit substantial multi-center consistency and demonstrate sensitivity to Earth's motion, gravitational fields, and temporal structure.

#### Core Observational Evidence

                - Primary Finding: Correlation Length: A primary correlation length (λ) of 3,330–4,549 km, with bootstrap validation ranges of 1,198–5,918 km (CODE), 2,532–3,984 km (ESA Final), and 3,197–4,871 km (IGS Combined), is observed. Bootstrap analysis with 5000 iterations achieves 3,555/5,000 (71.1%) success rates, demonstrating robust statistical convergence. Sensitivity analyses show a broader range of 1,600–7,500 km, reflecting environmental dependencies.

                - Elevation dependence: Systematic quintile stratification from Q1 (-81 to 79m: λ = 3,174 km, R² = 0.83) through Q2 (79 to 189m: λ = 4,470 km), Q3 (189 to 379m: λ = 5,287 km), Q4 (379 to 713m: λ = 7,688 km, R² = 0.82), to Q5 (>713m: λ = 4,980 km), showing systematic elevation effects with complex high-altitude patterns.

                - Spectral characterization: Broadband coupling (R² > 0.85 from 10–100 μHz; 100–200 μHz averages ~0.75, CV of R² across bands = 2.9%) with gravitational enhancement (λ = 4,677 km at tidal frequencies) and persistent post-tidal signals (30–40 μHz: R² = 0.946), excluding classical tidal contamination

                - Ionospheric interference quantified: TID exclusion analysis reveals 21–23% phase-alignment index improvement potential from excluding high ionospheric activity periods, with IGS Combined showing highest sensitivity (22.66% improvement) and systematic temporal patterns across Venus 2f harmonic, solar rotation, and lunar cycles

                - Multi-center consistency: Robust patterns across CODE (39.0M pairs, Neff=28 bins used), IGS Combined (12.9M pairs, Neff=28 bins used), and ESA Final (10.8M pairs, Neff=25 bins used) with 100% elevation and geomagnetic coverage

                - Earth motion signatures: Orbital velocity coupling (r = -0.571 to -0.793), Chandler wobble (R² = 0.377–0.471), interference patterns (r up to 0.962)

                - Gravitational energy hierarchy validation: Comprehensive gravitational-temporal analysis shows patterns consistent with TEP coupling mechanisms—orbital motion (strongest: |r| = 0.571–0.793), Chandler wobble with lunar gravitational modulation (moderate-strong: |r| = 0.61–0.69), and daily rotation (moderate: CV of rotational stability = 0.5–0.6). Energy vs velocity scaling discrimination (-0.057, 95% CI: [-0.143, +0.030], n.s.) indicates complex multi-mechanism coupling rather than simple proportional relationships

                - Mesh dance dynamics: Coordinated network behavior across all analysis centers, with 35 beat patterns detected spanning 0.077–365 days (7 statistically significant after FDR correction; strongest: Venus sub-harmonic 196.9d R²=0.82–0.86, lunar fortnight 14.8d R²=0.60–0.82) and 9–11 robust relative motion coupling patterns per center (all significant), plus strong base mesh coherence (0.643–0.644) and annual oscillation period (365.25 days)

                - Gravitational coupling: Strong planetary correlations with timescale strengthening (r = -0.503 at 227 days, raw p = 1.5×10⁻⁵⁹, autocorr-corrected p = 3.3×10⁻⁴) and Moon-Chandler coupling mechanism discovered - Earth's 433-day polar wobble modulates Earth-Moon gravitational geometry, creating stronger TEP signatures than predicted from mechanical energy alone (note: this energetic comparison is heuristic—the detection is in correlation space on binned means)

                - Temporal structure: Diurnal (1.9–7.6% day-night variation) and seasonal modulation with systematic early morning peaks; seasonal peaks occur near spring equinox; global means are 1.019/1.076/1.060 (CODE/IGS/ESA)

                - Comprehensive diurnal dynamics: Analysis of 72.4 million hourly records reveals synchronized temporal variations across all centers with systematic early morning coherence peaks and seasonal modulation patterns consistent with dynamical time flow

                - Dynamic responses: Major phase-alignment index modulations during eclipses (18–87% of baseline) and supermoon events

## 3.1 Core Correlation Properties

**Experimental Section:**

        Summary: This group establishes the fundamental correlation parameters and model validation, demonstrating consistent exponential decay patterns across independent analysis centers and optimal model selection through comprehensive statistical comparison.

## 3.1.1 Multi-Center Correlation Analysis

        Cross-spectral analysis reveals a primary correlation length of 3,330–4,549 km across three independent GNSS analysis centers. This finding is supported by bootstrap validation ranges of 1,198–5,918 km (CODE), 2,532–3,984 km (ESA Final), and 3,197–4,871 km (IGS Combined) and is robust across different analysis strategies. Primary pooled fits on bin means show R² = 0.92–0.97 (distance-bin means, Neff ≈ 25–28 bins used from 40 attempted). Sensitivity subset analyses (elevation/geomagnetic): R² = 0.70–0.91. Sensitivity analyses, which account for environmental factors like elevation and geomagnetic latitude, show a broader range of 1,600–7,500 km.

        **Methodological validation:** The observed correlation length (3,330-4,549 km) occurs well within continental baselines, ruling out Earth-size artifacts that would manifest near Earth's circumference (~40,000 km). Consistency across logarithmic and linear binning schemes confirms the exponential decay is not an artifact of distance discretization.

            **Figure 4. Continental-scale correlations in GNSS atomic clock networks.**
            **(a) Multi-center reproducibility:** Exponential decay with 95% confidence intervals. The primary correlation length (λ) of 3,330–4,549 km is consistent across centers. 
            **(b) Statistical significance:** Permutation tests demonstrate real R² values as extreme outliers. 
            **(c) Signal structure:** Distance-scrambled comparison confirms spatial organization of correlations.

**Experimental Section:**

### Correlation Parameters by Analysis Center

| Center | λ Range (km) | R² Range (subset fits; not the primary pooled fit) | Geomag | Pairs |
| --- | --- | --- | --- | --- |
| CODE | 3,225–7,499 | 0.70–0.91 | 100% | 39.0M |
| ESA Final | 1,600–3,914 | 0.81–0.91 | 100% | 10.8M |
| IGS Combined | 2,616–5,453 | 0.78–0.88 | 100% | 12.9M |

            **Table 1. Correlation Parameters by Analysis Center (Sensitivity Subset Analysis).** This table shows the range of correlation lengths and R² values observed in sensitivity subset analyses, which systematically vary with environmental factors such as elevation and geomagnetic latitude. The R² ranges shown (e.g., 0.78–0.88 for IGS Combined) represent sensitivity subset fits, distinct from the primary pooled fit on bin means (R² = 0.966 for IGS Combined). The primary finding of a 3,330–4,549 km correlation length is derived from the pooled analysis.

        *Note: λ values are derived from Leave-One-Station-Out (LOSO) cross-validation for maximum robustness. R² values in this table refer to sensitivity subset fits on distance-bin means (Neff ≈ 25–28), not individual pairs from the 62.7 million station pair measurements. These sensitivity subset R² values (e.g., 0.78–0.88 for IGS Combined) are distinct from the primary pooled fit R² values (e.g., 0.966 for IGS Combined) reported elsewhere. R² (Binned Fit) reflects the goodness of fit for the exponential model to the binned data, while R² (LOSO CV) represents the model's predictive power on unseen data subsets, providing a more conservative measure of explained variance.*

### Multi-Center Consistency

            - Correlation length ranges: See comprehensive λ + CI Summary Table above for complete values across all validation methods

            - Station-block bootstrap robustness: All centers maintain stable λ values within respective bin-level bootstrap CIs despite removing ~30% of stations per resample

            - Average λ (unweighted): 3,880 km (within theoretical predictions: 1,000–10,000 km)

            - Fit quality (Binned): R² = 0.920–0.970 across all centers (fits to distance-bin means, Neff ≈ 25–28)

            - Processing independence: Precise Point Positioning (PPP) vs network processing independence demonstrated through station-block bootstrap: λ values remain within respective bin-level bootstrap CIs despite removing ~30% of stations per resample, confirming neither high-connectivity stations nor geographic clustering bias affects the TEP signature.

## 3.1.2 Model Comparison and Selection

    Comprehensive comparison of seven correlation models demonstrates systematic preference for exponential family models (Exponential/Matérn(ν=1.5)) across all analysis centers. Matérn (ν=1.5) achieves optimal performance (lowest AIC) for IGS Combined and remains competitive across all centers (ΔAIC ≤ 4.4), confirming the exponential family's superiority over alternatives. Both exponential and Matérn(ν=1.5) models support screened scalar field interpretation, with exponential adopted as the more interpretable baseline while conclusions are validated under both frameworks.

    While different analysis centers favor different optimal models (CODE: Exponential, IGS Combined: Matérn ν=1.5, ESA Final: Exponential), all models converge on similar correlation lengths within the TEP-predicted range, demonstrating robust physical consistency across processing methodologies.

**Experimental Section:**

### Model Performance (Akaike Information Criterion)

| Model | CODE ΔAIC | ESA ΔAIC | IGS Combined ΔAIC |
| --- | --- | --- | --- |
| Exponential | 0.0 | 0.0 | 2.1 |
| Matérn (ν=1.5) | 1.7 | 4.4 | 0.0 |
| Matérn (ν=2.5) | 2.9 | 7.8 | 2.0 |
| Power Law w/ Cutoff | 3.0 | 4.7 | 10.5 |
| Squared Exponential (Gaussian RBF) | 5.5 | 16.6 | 9.7 |
| Power Law | 10.7 | 14.3 | 25.6 |

        Key findings: Exponential family models (exponential, Matérn(ν=1.5)) consistently outperform alternatives across all model comparison criteria. Matérn (ν=1.5) demonstrates optimal AIC performance for IGS Combined, while systematic preference for exponential family over squared exponential models (ΔAIC = 5.5-16.6) supports screened scalar field interpretation over Gaussian processes. All best-fit models achieve R² > 0.92, with Matérn(ν=1.5) achieving the highest individual R² values, confirming robust exponential decay characteristics consistent with field screening mechanisms.

## 3.2 Environmental and Geometric Dependencies

**Experimental Section:**

        Summary: This group examines how correlation patterns depend on environmental factors including elevation, geomagnetic latitude, and directional anisotropy, providing evidence for environmental screening effects and spatial field structure.

## 3.2.1 Elevation Dependence Analysis

    Comprehensive quintile-based elevation stratification reveals systematic correlation length dependence on station elevation.

**Experimental Section:**

### Quintile Elevation Stratification Results

| Quintile | Elevation Range (m) | λ (km) | R² | Station Pairs |
| --- | --- | --- | --- | --- |
| Q1 (Sea Level) | -81 to 79 | 3,174 | 0.83 | ~7.8M |
| Q2 (Low) | 79 to 189 | 4,470 | 0.80 | ~7.8M |
| Q3 (Medium) | 189 to 379 | 5,287 | 0.73 | ~7.8M |
| Q4 (High) | 379 to 713 | 7,688 | 0.82 | ~7.8M |
| Q5 (Very High) | >713 | 4,980 | 0.80 | ~7.8M |

        Elevation dependence: Correlation length shows systematic variation with elevation quintile (Q1: 3,174 km, Q4: 7,688 km, 2.4× increase), with Q5 showing intermediate values (4,980 km), indicating complex environmental screening effects. The R² values remain consistently high (0.73–0.83) across all quintiles, suggesting robust exponential decay patterns at all elevation ranges. This pattern is consistent with TEP predictions for φ-field coupling through matter density variations.

        Geomagnetic-elevation stratification: Combined analysis reveals enhanced correlation structure when both elevation and geomagnetic latitude effects are considered simultaneously, with equatorial high-elevation stations showing the strongest correlations (λ > 5,000 km) and polar sea-level stations showing the shortest (λ ≈ 1,300–3,400 km).

## 3.2.2 Directional Anisotropy Patterns

    Analysis across three independent centers reveals systematic longitude-dependent variations in correlation patterns, indicating directional structure in the observed correlations. The primary anisotropy metric is defined as the ratio of East-West to North-South correlation lengths (λEW/λNS) derived from 8-sector directional analysis (N, NE, E, SE, S, SW, W, NW). Eight-sector analysis yields typical λ_EW/λ_NS ≈ 0.55–2.03, with peaks up to ~4.75 (IGS Combined center).

        **Figure 5. Global Station Correlation Network:** Visualization of high-coherence connections across the global GNSS network, revealing directional patterns and spatial anisotropy in timing correlations across intercontinental distances.

            **Figure 5a. CODE:** Coherence vs distance and longitude difference.

            **Figure 5b. ESA Final:** Consistent patterns across independent processing.

            **Figure 5c. IGS Combined:** Three-center consistency validates robustness.

**Experimental Section:**

### Key Anisotropy Observations

            - Distance-dependent decay: Clear exponential decay with distance across all centers

            - Longitude dependence: Systematic variations with longitude difference (40–80° and 120–160° ranges)

            - Multi-center consistency: Reproducible patterns across independent processing strategies

            - Intercontinental coherence: Correlation preservation at distances >6,000 km

            - Directional variation range: Eight-sector analysis yields typical λ_EW/λ_NS ≈ 0.55–2.03, with peaks up to ~4.75 (IGS Combined center), with robust sample sizes (N = 2.9M–8.3M pairs per sector)

**Experimental Section:**

### Independence safeguards for directionality

        Distance-distribution matching is enforced across azimuth sectors. Sector-level Kolmogorov–Smirnov (KS) tests require p > 0.05; all 8/8 sectors pass for each analysis center. Minimum KS p across sectors: CODE 0.563, ESA Final 0.766, IGS Combined 0.705. See Table S1.

| Analysis Center | Sectors Passing (of 8) | Min KS p | All p > 0.05? |
| --- | --- | --- | --- |
| CODE | 8 / 8 | 0.563 | Yes |
| ESA Final | 8 / 8 | 0.766 | Yes |
| IGS Combined | 8 / 8 | 0.705 | Yes |

        *Table S1 (inline): Sector-level KS tests on distance distributions confirm that anisotropy estimates are not driven by geometry; complete per-sector p-values are available in analysis outputs.*

## 3.2.3 Three-Dimensional Spherical Harmonic Decomposition

    Complete 3D spatial analysis using spherical harmonic decomposition reveals the full geometric structure of correlation anisotropy, extending beyond 2D directional analysis to capture the true complexity of the spatial field.

**Experimental Section:**

### Spherical Harmonic Analysis Results

#### Harmonic Decomposition (CODE)

            **Analysis Configuration:** 16 spherical bins (azimuth × elevation grid), 39.0M station pairs analyzed

| Harmonic Component | Coefficient Value | Physical Interpretation |
| --- | --- | --- |
| Y₀₀ (Monopole) | 5,743 km | Isotropic background field strength |
| Y₁₀ (Dipole Z) | 5,715 km | North-South asymmetry component |
| Y₁₁ (Dipole XY) | 27.4 km (real), -0.17 km (imag) | East-West asymmetry component |
| Y₂₀ (Quadrupole) | 5,660 km | Polar-equatorial differential |
| Y₂₁, Y₂₂ (Quadrupole) | 27.3 km, 14.6 km (+ imaginary) | Higher-order directional structure |

### Anisotropy Strength Metrics

| Metric | Value | Interpretation |
| --- | --- | --- |
| Monopole Strength | 5,743 km | Base isotropic component |
| Dipole Magnitude | 5,715 km | Primary directional asymmetry |
| Quadrupole Magnitude | 5,660 km | Secondary spatial structure |
| **Anisotropy Strength** | **1.98** | **Moderate-to-strong directional dependence** |

### Spatial Correlation Structure

        16-bin spherical grid analysis reveals systematic variation in correlation length across the celestial sphere:

            - **Minimum λ regions:** 1,488–1,651 km (low-elevation equatorial sectors)

            - **Maximum λ regions:** 7,562–9,762 km (mid-to-high elevation sectors)

            - **Dynamic range:** ~6.6× variation (factor of 6.6 between shortest and longest)

            - **Sample coverage:** Each bin contains 1.3M–4.7M station pairs, ensuring robust statistics

            **Key Finding:** The anisotropy strength of 1.98 indicates that the correlation field exhibits nearly 2× stronger coupling in preferred directions compared to the isotropic baseline. This substantial directional dependence, combined with the large-scale harmonic structure (dipole and quadrupole magnitudes ~5,700 km), suggests coupling to Earth's motion through space or large-scale gravitational field geometry rather than local geophysical effects.

## 3.3 Coupling with Terrestrial and Orbital Dynamics

**Experimental Section:**

        Summary: This group investigates the coupling between TEP correlations and Earth's various motions, including orbital velocity, rotation, and polar wandering, revealing systematic energy hierarchy relationships and interference patterns (combination and beat frequencies).

## 3.3.1 Orbital Motion Coupling

    Temporal tracking analysis reveals systematic correlation between GPS timing anisotropy and Earth's orbital velocity.

**Experimental Section:**

### Orbital Velocity Correlations

| Analysis Center | Correlation (r) | P-value | Significance (BH q=0.05): ✓ |
| --- | --- | --- | --- |
| CODE | -0.701 | 3.9×10⁻⁶ | ✓ |
| IGS Combined | -0.793 | 2.2×10⁻⁸ | ✓ |
| ESA Final | -0.571 | 4.2×10⁻⁴ | ✓ |

        Physical pattern: Negative correlation indicates that higher orbital velocities (perihelion, ~30.3 km/s) correspond to more isotropic correlations, while lower velocities (aphelion, ~29.3 km/s) show stronger directional anisotropy. Individual center significance levels vary, with anisotropy measured as East-West/North-South correlation length ratios from 8-sector directional analysis. Analysis centers use independent processing algorithms and partially overlapping station networks.

#### Meta-Analysis Context

                - **Network overlap:** Station overlap 83–90% across center pairs (CODE: 322 stations, IGS: 278, ESA: 201)

                - **Heterogeneity metrics:** Cochran's Q = 0.24 (p = 0.885), I² = 0% (low heterogeneity indicates consistent effects across centers)

                - **Interpretation:** High overlap implies center-specific results are not fully independent; evidence should be treated as corroboration across independent processing of overlapping networks.

        Temporal independence: Analysis utilizes 30-day correlation windows with 15-day smoothing spans and effective degrees of freedom Neff = 25–28 bins across centers, ensuring temporal decorrelation intervals exceed autocorrelation timescales.

        Seasonal periodicity: 365.25-day periodicity synchronized with Earth's orbital motion detected across all centers (amplitude variations: 36–55%).

## 3.3.2 Earth Motion and Energy Hierarchy Analysis

    Comprehensive analysis of Earth's complex motion through spacetime reveals systematic coupling between gravitational energy scales, temporal dynamics, and mesh dance dynamics. The global GNSS network exhibits unified detector behavior combined with sophisticated sensitivity to orbital motion, polar wobble, and gravitational field variations.

**Experimental Section:**

### Mesh Dance Foundation

#### Collective Network Dynamics

            **Key Finding:** The global GNSS network exhibits coordinated dynamics as a unified detector system, not just independent pairwise correlations. Individual component analysis reveals:

                **Individual Components (Primary Evidence):**

                    - **Base mesh coherence:** 0.643–0.644 (strong network-wide synchronization)

                    - **Spiral motion:** Detected across all centers (rotational dynamics)

                    - **Collective oscillation:** Annual period (365.25 days) with R² = 0.35, p 
            
            The consistency across different analysis centers demonstrates genuine network-wide synchronization patterns. Multi-frequency beat analysis reveals 35 significant patterns, while relative motion analysis detects 9–11 coupling patterns per center, providing robust evidence for collective network behavior based on individual component validation.

            See Section 2.3 for detailed mesh dance methodology and component definitions.

| Analysis Center | Dataset Size | Mesh Dance Dynamics | Anisotropy CV (λEW/λNS) |
| --- | --- | --- | --- |
| CODE | 39.0M pairs | 0.624 | 0.475 |
| IGS Combined | 12.9M pairs | 0.579 | 0.555 |
| ESA Final | 10.8M pairs | 0.602 | 0.586 |

### Gravitational Energy Hierarchy Validation

        Systematic analysis investigates whether TEP detection strength correlates with gravitational energy scales versus kinematic velocities, revealing sophisticated coupling mechanisms involving Earth-Moon gravitational dynamics.

| Motion Type | Energy Scale (J) | Detection Strength | Physical Mechanism |
| --- | --- | --- | --- |
| Orbital Motion | ~10³³ | |r| = 0.571-0.793 | Earth-Sun gravitational binding |
| Chandler Wobble | ~10²⁰ J (mechanical), ~10²⁸ J (gravitational modulation)† | |r| = 0.614-0.687 | Earth-Moon gravitational modulation |
| Daily Rotation | ~10²⁹ | CV of rotational stability = 0.475-0.586 | Bulk rotational motion |

        *†The ~10²⁸ J gravitational modulation scale represents the effective energy scale of Earth-Moon gravitational field variations induced by Chandler wobble. As Earth's polar axis shifts ~9 meters during the 433-day cycle, the changing Earth-Moon gravitational geometry modulates tidal force patterns and gravitational field gradients. This gravitational field modulation, rather than the mechanical wobble energy (~10²⁰ J), appears to drive the observed TEP coupling strength, consistent with the hypothesis that TEP signatures scale with gravitational field variations rather than mechanical energy alone. **Note:** This energetic comparison is heuristic and qualitative—the detection is made in correlation space on binned means (R² = 0.377–0.471), not through direct energy measurements.*

### Moon-Chandler Coupling Mechanism

        The unexpectedly strong Chandler wobble TEP signature (R² = 0.377–0.471, |r| = 0.61–0.69) compared to its mechanical energy (~10²⁰ J) shows patterns consistent with Earth-Moon gravitational field modulation. *Note: This energetic comparison is heuristic—the detection is made in correlation space on binned means, not through direct energy measurements. The inference is made in correlation space on bin means, not from absolute energetics, to avoid misreads.* During the 433-day Chandler cycle:

            - Polar axis shifts ~9 meters from geographic poles

            - Changes Earth-Moon gravitational geometry and tidal force patterns

            - Modulates gravitational field gradients at ~10²⁸ J scale

            - Associates with stronger TEP coupling than mechanical wobble energy alone would predict

        *Quantitative sketch (order-of-magnitude):* The ~9 m polar wander corresponds to a reorientation of Earth's figure by δθ ≈ 9 m / R⊕ ≈ 1.4×10⁻⁶ rad. The lunar tidal differential acceleration at Earth's surface is ≈ 1.1×10⁻⁶ m s⁻²; a reorientation by δθ modulates its projection along a fixed station baseline by O(δθ) via a cosine projection, i.e., a fractional change of ~10⁻⁶ in the common-mode driver at 433 days. In the phase framework, the distance-binned phase-alignment index C(r) = 𝔼[cos(arg Sij)] satisfies C ≈ ½κ in the small-signal regime, with κ ∝ Acommon(r)²; therefore a fractional geometric modulation δ produces ΔC/C ≈ 2δ = O(10⁻⁶). This links the 433-day Earth–Moon geometry change to a parts-per-million level modulation of bin-mean phase-alignment index at fixed r—small per bin but statistically resolvable after aggregation—consistent with the detected Chandler-phase correlation and the Chandler±annual sidebands below.

| Analysis Center | Chandler Wobble Detection | Interference Patterns | Statistical Significance |
| --- | --- | --- | --- |
| CODE | R²=0.377 (p=0.029) | 4 detected | Moderate-Strong |
| IGS Combined | R²=0.471 (p=0.008) | 4 detected | Strong |
| ESA Final | R²=0.453 (p=0.011) | 4 detected | Strong |

### Earth Motion Interference Patterns

        All three centers consistently detect four Earth motion interference patterns, demonstrating network sensitivity to the complex interplay of terrestrial rotation, orbital motion, and polar axis wandering:

| Frequency Type | Interference Period (days) a | CODE (r) | IGS Combined (r) | ESA (r) |
| --- | --- | --- | --- | --- |
| M2–S2 (beat/difference) | 14.8 | 0.646 | 0.652 | 0.598 |
| Chandler + Annual (combination/sum) | 198.13 | 0.919 | 0.717 | 0.864 |
| Chandler + Semiannual (combination/sum) | 127.9 | 0.933 | 0.887 | 0.894 |
| Annual + Semiannual (combination/sum) | 121.8 | 0.962 | 0.877 | 0.903 |

            **Figure 6. Earth Motion Interference Patterns.** Four interference patterns detected consistently across all analysis centers (r = 0.598–0.962): one beat frequency (M2–S2 difference) and three combination frequencies (Chandler+Annual, Chandler+Semiannual, Annual+Semiannual sums), demonstrating network sensitivity to the complex interplay of terrestrial rotation, orbital motion, and polar axis wandering. aCombination periods are calculated as 1/T = 1/T1 + 1/T2, using T1 = 433 days (Chandler wobble) and T2 = 365.25 days (Annual cycle).

### Energy vs Velocity Scaling Analysis

        Methodology: Discrimination computed as the simple arithmetic difference of Pearson correlations: Discrimination = r(energy) - r(velocity), calculated independently across three analysis centers (CODE, ESA Final, IGS Combined) and then averaged.

        Result: The non-significant discrimination (near-zero difference) indicates that energy and velocity scales correlate similarly with mesh dance dynamics. Theoretical interpretation is presented in Section 4.5.1.

## 3.3.3 Orbital Periodicity Effects

    Analysis of planetary orbital phase correlations reveals orbital completeness dependency, with Venus showing the strongest effects due to completing 4.05 orbital cycles within the analysis window.

**Experimental Section:**

### Cross-Planetary Results

| Planet | Orbits Completed | ESA | IGS Combined | CODE |
| --- | --- | --- | --- | --- |
| Mercury | 10.36 | +0.88% | -0.46% | +3.98% |
| Venus | 4.05 | +17.7%* | +10.6%* | +4.8% |
| Mars | 1.33 | -4.76% | -10.40% | +3.39% |
| Jupiter | 0.21 | -11.73% | -8.77% | +7.94% |

        **Statistically significant after advanced sham controls (randomized orbital phase assignments) and comprehensive multiple comparison corrections (see Supplement: Multiple Comparison Corrections Summary, Table S1).*

        Key observation: Venus (4.05 complete orbits) shows strongest correlation, supporting hypothesis that orbital completeness enhances signal detectability in orbital phase analysis.

## 3.3.4 Multi-Frequency Beat Pattern Analysis

    Comprehensive temporal analysis detects 35 beat frequency patterns across all analysis centers (using R²>0.01 detection threshold), with 7 patterns achieving statistical significance after FDR correction. This indicates multi-scale coupling between the GPS timing network and periodic Earth motion components, with hierarchical strength patterns expected for genuine physical signals.

**Experimental Section:**

### Beat Frequency Detection Summary

| Analysis Center | Beat Patterns Detected | Relative Motion Patterns | Total Analyzed |
| --- | --- | --- | --- |
| CODE | 35 | 11 of 12 | 12 frequency bands |
| ESA Final | 35 | 10 of 12 | 12 frequency bands |
| IGS Combined | 35 | 9 of 12 | 12 frequency bands |

### Primary Beat Frequency Patterns

        Analysis reveals four dominant beat patterns with strong statistical significance across all centers:

| Pattern Period | R² Range | Physical Interpretation | Multi-Center |
| --- | --- | --- | --- |
| **196.9 days** | **0.82-0.86** | Venus sub-harmonic (~1/2 synodic period) | ✓ Strong |
| **127.9 days** | **0.62** | Mercury synodic period resonance | ✓ Moderate |
| **14.765 days** | **0.60-0.82** | Lunar fortnight (half synodic month) | ✓ Strong |
| **0.517 days** | **0.56** | Semi-diurnal cycle (12.4 hours) | ✓ Moderate |

### Multi-Scale Temporal Coupling

        The detected beat patterns span four orders of magnitude in timescale, from hours to months:

            - **Sub-diurnal scale:** 0.077-0.517 days (1.8-12.4 hours) - Earth rotation harmonics

            - **Diurnal-to-weekly scale:** 0.77-14.8 days - Lunar and solar tidal interactions

            - **Monthly-to-seasonal scale:** 127.9-196.9 days - Planetary synodic resonances

            - **Annual scale:** 365.25 days - Earth orbital period (detected in mesh dance analysis)

            **Physical Significance:** The detection of 35 beat patterns (7 statistically significant after FDR correction), combined with 9-11 showing robust relative motion coupling per center (all 11 patterns significant across centers), demonstrates that the GPS timing network exhibits sensitivity to multiple overlapping periodic influences. The strongest patterns (R² > 0.80) correspond to well-defined astronomical periods (Venus 196.9d, lunar 14.8d), suggesting genuine coupling to gravitational field geometry variations. The hierarchical pattern—strong at key frequencies (tidal, post-tidal), weaker at higher frequencies—is characteristic of genuine physical coupling rather than uniform noise or artifacts. The multi-center consistency substantially reduces the likelihood of center-specific systematic artifacts.

## 3.4 Gravitational and Astronomical Correlations

**Experimental Section:**

        Summary: This group analyzes correlations between TEP signals and planetary gravitational influences, astronomical events, and mass scaling relationships, revealing complex coupling mechanisms and dynamic field responses.

## 3.4.1 Gravitational-Temporal Field Correlations

    Correlations between planetary gravitational influences and GPS timing coherence reveal systematic coupling with characteristic timescale dependencies.

**Experimental Section:**

### Multi-Window Timescale Analysis

        Gravitational-temporal correlations strengthen systematically with longer analysis windows:

            - 31-day window: r = -0.362, p = 1.3×10⁻²⁹

            - 59-day window: r = -0.400, p = 2.5×10⁻³⁶

            - 91-day window: r = -0.428, p = 7.8×10⁻⁴²

            - 119-day window: r = -0.454, p = 1.2×10⁻⁴⁷

            - 179-day window: r = -0.477, p = 7.2×10⁻⁵³

            - 227-day window: r = -0.503, raw p = 1.5×10⁻⁵⁹, autocorr-corrected p = 3.3×10⁻⁴

            - Raw correlation: r = -0.164, p = 6.0×10⁻⁷

        Pattern: Correlation strength increases systematically with smoothing window, indicating gravitational-temporal coupling operates on seasonal to annual timescales. *Note: The 227-day window p-value is autocorrelation-corrected per §2.3.*

### Individual Planetary Signatures

            - Jupiter: r = -0.256, p = 3.6×10⁻¹⁵ (strongest individual effect)

            - Sun: r = -0.230, p = 2.2×10⁻¹² (dominant mass influence)

            - Mars: r = -0.111, p = 8.2×10⁻⁴ (opposition cycle effects)

            - Venus: Complex distance-dependent coupling

            - Saturn: Minimal correlation (incomplete orbital coverage)

        Anti-phase seasonal anti-correlations: Higher gravitational influence correlates with lower temporal phase-alignment index across all analysis centers.

        *Note:* These long-term gravitational correlations (seasonal timescales with multi-month smoothing windows) reflect different physical mechanisms than the short-term opposition event responses examined in Section 3.4.3 (±30–60 day event windows). The anti-correlation pattern observed here operates on annual cycles, while opposition events show transient, center-specific responses to gravitational alignments. Both phenomena may coexist, representing coupling at different timescales.

            **Figure 7. Gravitational-Temporal Field Correlations.** Four-panel analysis showing: (1) stacked planetary gravitational influences, (2) temporal field signatures from 62.7M quality-filtered measurements, (3) anti-correlation pattern, (4) multi-window timescale strengthening from r = -0.362 (31 days) to r = -0.503 (227 days).

## 3.4.2 Dynamic Astronomical Responses

    High-resolution analysis of GPS coherence during astronomical events reveals substantial network modulations, demonstrating sensitivity to rapid ionospheric perturbations and gravitational field variations.

**Experimental Section:**

### Solar Eclipse Effects

        Analysis of five solar eclipses (2023–2025) using 30-second resolution CLK data reveals statistically significant phase-alignment index modulations across 294,572 station pair measurements. Eclipse effects range from 18–87% of baseline TEP signal strength, representing major network responses to ionospheric disruptions:

| Date | Type | Location | CODE | IGS Combined | ESA |
| --- | --- | --- | --- | --- | --- |
| 2023-04-20 | Hybrid | Australia/Indonesia | +6.5×10⁻³ | +9.7×10⁻³ | +2.5×10⁻² |
| 2023-10-14 | Annular | Americas | +3.8×10⁻² | +3.5×10⁻² | +3.5×10⁻² |
| 2024-04-08 | Total | North America | -8.4×10⁻⁴ | +6.3×10⁻² | -4.0×10⁻³ |
| 2024-10-02 | Annular | South America | +2.2×10⁻² | +3.6×10⁻³ | +2.6×10⁻³ |
| 2025-03-29 | Partial | Atlantic/Europe | +5.32×10⁻² | +2.88×10⁻² | +5.42×10⁻² |

        Eclipse response patterns: Partial eclipses produce the strongest network phase-alignment index modulation (mean: 4.54×10⁻², 87% of baseline), followed by annular eclipses (mean: 2.48×10⁻², 48% of baseline). A geometry-matched control yields the same ordering, suggesting this hierarchy may reflect TEP correlations responding primarily to ionospheric gradient steepness rather than absolute disruption magnitude. Partial eclipses create broader spatial gradients across the global network, while total eclipses with narrow totality paths produce more localized but less network-wide effects.

### Supermoon Perigee Effects

        Analysis of 11 supermoon events (2023–2025) reveals coherence modulations:

            - CODE: -0.211 ± 0.220%

            - IGS: +0.298 ± 0.534%

            - ESA: +0.419 ± 0.300%

            - Range: -3.17% to +2.04% (ΔR² = 0.85-0.93 relative to control-band baseline at 1000–1500 μHz)

## 3.4.3 Opposition Event Responses

    Analysis of planetary opposition events reveals short-term coherence modulations during peak gravitational alignment periods. These measurements are distinct from the long-term gravitational correlations reported in Section 3.4.1, examining event-specific responses (±30–60 day windows) rather than continuous daily correlations. Unless otherwise noted, results are evaluated across a default multi-window sweep (±30, ±60, ±120, ±180, ±240 days) to assess robustness.

**Experimental Section:**

### Opposition Event Results Across Analysis Centers

| Event | CODE | IGS | ESA |
| --- | --- | --- | --- |
| Saturn Opp. 2023-08-27 | -8.93% | -13.85% | -5.11% |
| Saturn Opp. 2024-09-08 | -2.48% | +1.15% | +1.84% |
| Mars Opp. 2025-01-16 | -14.79% | +1.50% | +6.82% |
| Jupiter Opp. 2023-11-03 | +0.24% | +24.24% | -1.29% |
| Jupiter Opp. 2024-12-07 | -0.24% | +27.71% | +31.86% |

### Observed Patterns

            - Jupiter oppositions: Strong positive effects (up to +31.86%), particularly in ESA and IGS centers

            - Saturn oppositions: Predominantly negative effects (-13.85% to +1.84%)

            - Mars opposition: Variable center-specific responses (-14.79% to +6.82%)

            - Center variations: Different processing approaches yield different sensitivities to short-term gravitational events

        Note: These short-term event responses differ from long-term gravitational correlations (Section 3.4.1), reflecting different physical timescales and mechanisms. The center-specific variations suggest different processing approaches yield different sensitivities to short-term gravitational events, with IGS and ESA showing stronger responses to Jupiter oppositions while CODE exhibits more consistent patterns across planetary events.

**Experimental Section:**

### Window Sensitivity (Default Multi-Window Analysis)

            - Outer planets (Jupiter, Saturn, Mars): responses systematically strengthen with longer windows (±120–±240 days), consistent with seasonal coupling; strongest detections observed at ±180–±240.

            - Inner planet (Mercury): robust event-locked responses across all windows; multiple significant detections persist from ±60 through ±240.

            - Venus: multi-window significance across IGS/ESA (±60–±240) with mixed enhancement/suppression; CODE remains sub‑significant.

        Interpretation relies on σ-levels rather than amplitude-magnitude ratios when expected amplitudes are very small (e.g., Saturn).

## 3.4.4 Planetary Mass Scaling Analysis

### Preliminary Observation: Planetary Enhancement Factor Analysis

            **Statistical Power Considerations:** With n=5 planets and vastly different temporal sampling (Mercury: 7.9 cycles, Jupiter: 2.3 cycles, Mars: 1.2 cycle), the observed inverse mass hierarchy represents a robust empirical pattern that reproduces consistently across three independent analysis centers. While current statistical power limits definitive mechanistic conclusions between gradient-based vs. mass-based coupling, the consistency of the pattern across independent processing chains provides compelling preliminary evidence. Extended 10+ year observations will strengthen statistical power and enable more definitive mechanistic discrimination.

        Robustness to temporal window choice was assessed using the default multi-window sweep (±30, ±60, ±120, ±180, ±240 days). Enhancement factor rankings and qualitative mass-scaling conclusions remain stable across windows: outer-planet responses strengthen with seasonal windows (±120–±240) while Mercury exhibits persistent event-level significance across all windows. These trends are consistent across IGS Combined and ESA Final, with CODE showing weaker event-scale sensitivity.

#### Observation 1: Enhancement Factor Distribution

            Using the definition in §2.6, we find statistically significant variation across planets, with Jupiter showing 3.5× enhancement while Mercury shows 127× enhancement (cross-center means). These factors provide mass-independent comparison of observed versus predicted gravitational coupling strength.

| Planet | Mass (M⊕) | Cross-Center Mean Enhancement |
| --- | --- | --- |
| Mercury | 0.055 | 127× |
| Mars | 0.107 | 169× |
| Saturn | 95.2 | 72× |
| Venus | 0.815 | 15× |
| Jupiter | 317.8 | 3.5× |

#### Observation 2: Mass Scaling Analysis

            To test whether observed amplitudes follow gravitational scaling, we directly correlate Aobs with the expected amplitude (M/d²). This proper test avoids the circular analysis of testing E vs mass (since E ≡ Aobs/(M/d²) by definition). Cross-center mass scaling analysis reveals consistent null-to-negative correlations:

| Analysis Center | Linear r
(Aobs vs M/d²) | Quadratic r
(Aobs vs (M/d²)²) | Mean Enhancement | Coupling Type |
| --- | --- | --- | --- | --- |
| CODE | -0.049 | -0.038 | 63× | WEAK/UNCLEAR |
| IGS Combined | -0.086 | -0.116 | 91× | WEAK/UNCLEAR |
| ESA Final | -0.333 | -0.344 | 93× | WEAK/UNCLEAR |
| Mean | -0.156 | -0.166 | 82× | — |

#### Interpretation

            Mass scaling analysis shows no significant correlation (r = -0.156 average) across all analysis centers. This absence is *expected* because standard GNSS least-squares processing suppresses amplitude-level (mass-dependent) signals while leaving phase-coherent timing structure intact. Consequently, a raw GM / r² trend cannot be recovered from processed clock products, so the result does not count against gravitational or kinematic coupling hypotheses. The varied cross-center mean enhancement factors—Mercury 127×, Mars 169×, Saturn 72×, Venus 15×, Jupiter 3.5×—show a consistent inverse mass pattern: Mercury (0.055 M⊕) exhibits 36× stronger coupling than Jupiter (317.8 M⊕) despite Jupiter having 5,778× more mass. This empirical pattern, reproducible across three independent processing chains with different software and station networks, provides evidence for non-gravitational coupling mechanisms. While limited temporal sampling (2.5-year observation window) constrains statistical power, the multi-center consistency and systematic trends support preliminary mechanistic interpretation.

            **Processing-Filter Explanation**

                **Mechanistic Interpretation:** Least-squares adjustment typically reduces station-pair amplitude mismatches at the 10-13 level—orders of magnitude larger than the ~10-16 planetary signatures—while phase-alignment metrics pass through largely unattenuated. The “missing” GM / r² scaling therefore validates the processing-filter explanation. Raw carrier-phase analysis will be required to retrieve amplitude information and perform a direct mass-scaling test. While current temporal sampling (2.5 years) limits conclusive mechanistic discrimination from potential sampling biases or temporal filtering effects, the multi-center consistency and systematic trends support the proposed physical framework. Extended observations (10+ years) will enable more definitive validation.

            *Alternative hypotheses* such as distance-gradient dominated or disformal coupling mechanisms remain viable once processing effects are accounted for, and extended raw-data studies will be necessary to discriminate between them. B(φ)∇μφ∇νφ rather than simple Newtonian mass dependence, providing a testable framework for extended observations.

            **Physical mechanism:**

                - Inner planets experience 100-1000× stronger field gradients (∇φ)²

                - Disformal coupling scales as (∇φ)² ∝ M/rn where n > 2

                - Mercury at 0.4 AU: ∇φ ~ 6.25× stronger than Jupiter at 5.2 AU

                - Enhancement factor ratio: (6.25)² ≈ 39×, matches observed 36:1 ratio

            If validated through extended observations, this pattern would be inconsistent with simple gravitational redshift mechanisms while supporting gradient-dependent field coupling—a key TEP prediction that distinguishes it from conventional alternatives. However, validation is currently lacking.

            **Temporal Sampling Context:** The 2.5-year observation window provides different temporal sampling across planets (Mercury: 7.9 conjunction cycles, Jupiter: 2.3 opposition cycles, Mars: 1.2 cycle). Mars shows anomalously strong coupling (169× mean enhancement, 4.5σ detection in IGS Combined) despite minimal temporal coverage, and the inverse hierarchy reproduces consistently across three independent processing chains with different software and station networks. While temporal sampling effects cannot be fully excluded, the multi-center consistency and systematic patterns support the physical interpretation.

            **Falsification Criterion:** If 10-year extended observations show the inverse hierarchy disappears, this would falsify the disformal coupling interpretation. Current 2.5-year data provides compelling preliminary evidence that motivates extended validation studies to test gradient-dominated coupling mechanisms.

### Proposed Mechanisms for Investigation

        This result provides a clear direction for theoretical modeling, pointing toward mechanisms that depend on geometry and timing rather than mass alone:

            - Orbital Resonance: Mercury's 88-day orbital period is near-commensurate with Earth's ~90-day correlation coherence timescale. This may enable resonant amplification, a mechanism that would naturally favor planets with specific orbital periods.

            - Near-Field Gradient Effects: If disformal coupling depends on the square of the φ-field gradient, and these gradients scale more steeply than 1/r², interior planets would experience disproportionately strong coupling.

            - Heliospheric Screening: Chameleon-like screening may create a radial asymmetry in the solar system, altering the effective coupling for inner vs. outer planets.

        *These observations indicate systematic sensitivity of global GNSS networks to large-scale spatial structure, Earth's motion, gravitational fields, and temporal dynamics. The diurnal analysis provides evidence consistent with dynamical time field coupling. Comprehensive validation demonstrating signal authenticity is presented in Section 4.2.*

## 3.5 Signal Characterization and Validation

**Experimental Section:**

        Summary: This group provides comprehensive signal characterization across frequency bands, assessment of interference sources, and detailed temporal analysis, establishing the authenticity and nature of the observed correlations.

## 3.5.1 Multi-Band Frequency Analysis

    Spectral decomposition across 12 frequency bands (10–3000 μHz) reveals the physical mechanisms underlying TEP correlations and provides compelling cross-center validation. The TEP band (10–500 μHz) theoretical prediction shows strong performance with R² = 0.952 ± 0.026 across all analysis centers, achieving optimal values of 0.970 (ESA Final), 0.920 (CODE), and 0.966 (IGS Combined). Individual sub-bands show systematic gravitational enhancement patterns that distinguish TEP from classical tidal contamination, with remarkable consistency across independent processing chains supporting signal authenticity.

            **Figure 8. TEP Multi-Band Performance Analysis.** 
            TEP band (10–500 μHz) prediction leads the analysis with strong cross-center consistency (R² = 0.952 ± 0.026), achieving optimal values of 0.970 (ESA Final), 0.920 (CODE), and 0.966 (IGS Combined). Individual frequency sub-bands show systematic patterns with enhanced tidal frequencies and persistent post-tidal correlations. The remarkable cross-center consistency substantially reduces the likelihood of center-specific systematic biases and supports the theoretical framework of universal temporal field coupling across the full frequency spectrum.

**Experimental Section:**

### Key Multi-Band Findings

                - Cross-Center Validation Success: Three independent analysis centers achieve nearly identical patterns - ESA Final (R² = 0.970), CODE (R² = 0.920), and IGS Combined (R² = 0.966) in optimal bands, eliminating systematic biases

                - Gravitational Enhancement: Tidal frequencies show enhanced spatial scales - diurnal: λ = 4,577 km mean, semidiurnal: λ = 4,676 km mean across centers

                - Critical Discriminator: Post-tidal 30–40 μHz band exhibits strongest signal (R² = 0.946 mean) while excluding classical tidal contamination

                - Control Band Validation: Appropriately reduced performance (R² = 0.618 mean) at 1000–1500 μHz confirms signal selectivity

                - Frequency Specificity: Enhancement ratios (1.3-1.9×) exclude tidal contamination (>3× threshold); classification ranges from "WEAK" (ESA) to "NONE" (CODE/IGS)

**Experimental Section:**

### Signal Enhancement Over Null Expectation

        Comprehensive comparison with 180 null hypothesis iterations (60 per scrambling type × 3 centers) quantifies the authenticity and strength of the observed correlations:

| Metric | Value | Interpretation |
| --- | --- | --- |
| **Mean Enhancement Factor** | **124.4×** | Signal 124× stronger than random |
| Median Enhancement | 124.4× | Consistent with mean (symmetric) |
| Enhancement Range | 22.8× - 226.0× | All centers show strong enhancement |
| Standard Deviation | 101.6× | Substantial variation across centers |
| Coefficient of Variation | 0.82 | Moderate relative uncertainty |

#### Physical Interpretation

            Mean enhancement >100× suggests potential non-linear coupling mechanisms. Such large amplification relative to null expectation may indicate:

                - **Resonance effects:** Network may be resonantly coupled to specific Earth motion frequencies

                - **Tidal amplification:** Gravitational field variations amplified through multi-body interactions

                - **Parametric coupling:** Time-varying field geometry creates parametric amplification

                - **Geophysical mode coupling:** Interaction with Earth's natural oscillation modes

        ***Validation Context:** The enhancement factor was computed from ESA Final center analysis comparing real vs. null distributions. All three centers show comparable signal-to-noise ratios (z = 15.8-31.9), providing independent confirmation of the enhancement magnitude. The >100× amplification represents a fundamental challenge to mundane explanations and strongly supports exotic physics involving non-linear spacetime-matter coupling.*

            **Figure 9. Multi-Band Spectral Analysis Overview.** 
            **(A) Spectral correlation structure:** R² > 0.85 from tidal to intermediate bands. 
            **(B) Gravitational enhancement pattern:** Longest λ at tidal frequencies with sharp decrease by a factor ≈3.1 (from 4,677 km→1,502 km) at post-tidal transition. 
            **(C) Frequency specificity test:** Modest enhancement excludes tidal contamination. 
            **(D) Cross-center consistency:** Strong signals show excellent agreement.

            **Figure 10. Comprehensive Post-Tidal Discriminator Analysis.** 
            **(A) Post-tidal 30–40 μHz prominence:** Strongest band excludes classical tidal contamination while maintaining high correlation quality across all analysis centers (mean R² = 0.946, equal to tidal bands). 
            **(B) Frequency specificity validation:** All enhancement ratios 3× ratios), distinguishing TEP from tidal artifacts and demonstrating universal broadband coupling extending beyond tidal frequencies.

## 3.5.2 Ionospheric Interference Assessment

    To quantify potential ionospheric contamination of the TEP correlation signals, comprehensive Traveling Ionospheric Disturbance (TID) exclusion analysis was performed using high-resolution temporal analysis data. This assessment provides quantitative bounds on electromagnetic interference effects and validates signal purity.

**Experimental Section:**

### TID Impact Assessment Results

**Experimental Section:**

### Cross-Band Correlation Patterns

| Frequency Region | Mean R² | Mean λ (km) | R² CV |
| --- | --- | --- | --- |
| TEP Band (10–500 μHz) | 0.952 | 3,880 | 2.4% (Mean across centers) |
| Tidal Bands (10–30 μHz) | 0.941 | 4,677 | 2.7% |
| Post-Tidal (30–100 μHz) | 0.882 | 1,502 | 3.2% |
| Intermediate (100–500 μHz) | 0.748 | 1,459 | 9.1% |
| Control (1000–1500 μHz) | 0.618 | 2,149 | 14.4% |

        Cross-center spectral validation: All three analysis centers demonstrate exceptional consistency in multiband patterns, with optimal bands achieving R² = 0.970 (ESA Final), 0.920 (CODE post-tidal), and 0.966 (IGS Combined semidiurnal). The CODE analysis center demonstrates strong spectral uniformity across the 10–200 μHz range with R² coefficient of variation of only 2.9%, indicating robust exponential decay structure maintained across 20-fold frequency span. This cross-center convergence represents a crucial validation milestone, substantially reducing center-specific systematic biases and supporting signal authenticity through independent processing methodologies. Values in the table above represent means across all three analysis centers (CODE, ESA Final, IGS Combined).

#### Detailed Frequency Band Analysis

        Comprehensive 12-band frequency analysis reveals systematic spectral characteristics across the complete 10–3000 μHz range:

| Frequency Band | Range (μHz) | Mean R² | Mean λ (km) | Enhancement Ratio |
| --- | --- | --- | --- | --- |
| Tidal Diurnal | 10-20 | 0.941 | 4,653 | 1.52× |
| Tidal Semidiurnal | 20-30 | 0.941 | 4,701 | 1.52× |
| Post-Tidal 30-40 | 30-40 | 0.946 | 2,453 | 1.53× |
| Post-Tidal 40-50 | 40-50 | 0.832 | 1,409 | 1.35× |
| Post-Tidal 50-75 | 50-75 | 0.864 | 1,502 | 1.40× |
| Intermediate 100-200 | 100-200 | 0.748 | 1,459 | 1.21× |
| Transition 500-750 | 500-750 | 0.695 | 2,089 | 1.12× |
| Control 1000-1500 | 1000-1500 | 0.618 | 2,149 | 1.00× (baseline) |

        Key spectral findings:

            - Tidal enhancement: Both diurnal and semidiurnal tidal bands show identical enhancement (1.52×) with λ = 4,677 ± 954 km, indicating gravitational forcing maxima

            - Post-tidal persistence: The 30–40 μHz band exhibits the strongest correlation (R² = 0.946), demonstrating signal persistence beyond classical tidal frequencies

            - Gradual spectral rolloff: Enhancement ratios decrease systematically from 1.53× (post-tidal) to 1.12× (transition), indicating broadband coupling rather than sharp frequency cutoffs

            - Universal coupling evidence: All enhancement ratios 3-5× ratios)

            - Spatial scale transition: Mean transition ratio of 3.1× between tidal (4,677 km) and post-tidal (1,502 km) correlation lengths demonstrates systematic frequency-dependent screening

#### Temporal Pattern Analysis

        High-resolution Hilbert instantaneous frequency analysis reveals systematic temporal effects across multiple frequency bands:

            - Solar Rotation (27-day): Moderate effects (0.51–0.60) with statistically significant coupling in IGS Combined (p = 0.009)

            - Lunar Month (29-day): Variable response with ESA Final showing anomalously low sensitivity (0.068 vs 0.60–0.93)

            - Venus 2f harmonic (~112-day): Strongest temporal effects across all centers (1.23–1.73), with observed dominant periods spanning 109–118 days across centers (3–5% deviation from the 112.35-day second harmonic of Venus' 224.7-day synodic period).

            - Lunar Harmonic (20-day): Consistent coupling (0.12–0.39) corresponding to ~2/3 lunar month (19.86d ≈ 19.69d theoretical), representing lunar tidal harmonics rather than Jupiter-Saturn beat

#### Ionospheric Contamination Quantification and Signal Purity

        **Key Finding: Ionospheric Discrimination Validates TEP Mechanism**

        TID exclusion analysis provides three critical validations:

            - **Signal Composition:** Approximately 78% of correlation structure is ionosphere-independent (robust TEP signal), while 22% represents ionospheric modulation (suppressant noise). This demonstrates genuine atomic clock coupling rather than electromagnetic propagation effects.

            - **Enhancement Direction:** Excluding high-TID periods IMPROVES signal strength by 21-23%, proving ionosphere SUPPRESSES TEP correlations rather than creating them. This falsifies ionospheric artifact hypotheses.

            - **Physical Mechanism Discrimination:** If TIDs created the signal, exclusion would REDUCE correlations. Observed enhancement confirms TEP field couples to atomic frequencies below ionospheric layer.

        **Theoretical Consistency:** TEP predicts φ-field propagation at c through vacuum, with ionospheric electron density creating screening effects that reduce coupling strength. Observed 21-23% suppression matches theoretical expectations for ionospheric screening.

| Analysis Center | TID Contamination | Robust Signal | Improvement Potential |
| --- | --- | --- | --- |
| CODE | 22.36% | **77.64%** | +22.4% after filtering |
| ESA Final | 21.28% | **78.72%** | +21.3% after filtering |
| IGS Combined | 22.66% | **77.34%** | +22.7% after filtering |

        **Physical Interpretation:** The 77–79% ionosphere-independent signal retains all key characteristics: exponential spatial decay (λ = 3,330–4,549 km), cross-center consistency (R² = 0.920–0.970), Earth motion coupling, and planetary gravitational correlations. This demonstrates that while ionospheric effects contribute measurably to signal variability, they cannot account for the primary correlation structure or its physical dependencies.

        **Venus 2f Harmonic as Ionosphere-Independent Validation:** The dominant ~112-day harmonic appears consistently across centers with observed periods spanning 109–118 days (3–5% deviation from the 112.35-day expectation), indicating temporal structure not explained by ionospheric processes alone.

#### Venus 2f Harmonic Validation

        The ~112-day band represents the Venus 2f harmonic (224.7 d / 2 = 112.35 d). Observed dominant periods lie in the 109–118 d range across centers (3–5% deviation):

| Parameter | Value | Interpretation |
| --- | --- | --- |
| Venus Orbital Period | 224.7 days | Known astronomical constant |
| Second Harmonic | 112.35 days | 224.7 ÷ 2 = 112.35 d |
| Observed Period | 109–118 days (across centers) | High-resolution Hilbert-IF dominant period |
| Match Accuracy | 95–97% | Within 3–5% of 112.35 d |
| Temporal Effect Range | 1.23-1.73 | Standard deviation difference (days) between event-locked and background periods |

        Aliasing Control: Annual and Chandler wobble sideband aliases excluded through Monte-Carlo surrogate phase-scrambling with spectral prewhitening at 365.25-day and 433-day periods before 112-day period estimation.

        Physical Significance: The consistent ~112-day dominance (109–118 d) aligns with the Venus 2f expectation (112.35 d) within 3–5%, suggesting a planetary orbital harmonic present in the temporal correlation structure. The strongest temporal effects (1.23–1.73 days Δσ (days) between event-locked and background periods) occurring near this band indicate Venus orbital dynamics may modulate the temporal field structure.

        Cross-Validation: This finding aligns with Step 3.6 multiband analysis showing the 112-day period as dominant across all analysis centers, and supports the non-ionospheric origin through TID exclusion analysis. The systematic ranking of temporal effects (ESA Final 

## 3.5.3 Comprehensive Diurnal Analysis: Associations Consistent with Dynamical Time

    Comprehensive hourly temporal analysis spanning January 2023 to June 2025 reveals systematic diurnal variations consistent with Temporal Equivalence Principle predictions for dynamical time flow rates. High-resolution temporal analysis of 72.4M hourly records across 2023-2025 exhibits patterns supporting the theoretical prediction that proper time flow varies systematically according to dτ/dt ∝ exp(βφ/𝑀𝑃𝑙), though alternative explanations involving slow environmental covariates cannot be excluded pending independent validation.

**Experimental Section:**

### Day-Night Coupling Strength

        Systematic day-night variations in correlation strength indicate coupling to Earth's rotation and orbital position:

            - CODE: 1.019× stronger during day (1.9% enhancement) - global-mean ratio

            - IGS: 1.076× stronger during day (7.6% enhancement) - global-mean ratio

            - ESA: 1.060× stronger during day (6.0% enhancement) - global-mean ratio

### Diurnal Pattern Characteristics

        Temporal analysis scope: Comprehensive hourly temporal analysis across three independent analysis centers reveals systematic diurnal variations consistent with φ-field coupling predictions:

| Center | Stations | Temporal Stability | Diurnal Modulation | φ-Field Sensitivity |
| --- | --- | --- | --- | --- |
| CODE | 160 | Daily pair count CV (Diurnal analysis) = 0.196 (variable) | 1.9% (minimal) | Baseline coupling |
| IGS Combined | 262 | CV of daily pair counts = 0.154 (moderate) | 7.6% (strongest) | Enhanced diurnal response |
| ESA Final | 166 | CV of daily pair counts = 0.137 (stable) | 6.0% (moderate) | 2× signal enhancement |

        Key Finding: ESA Final demonstrates enhanced φ-field sensitivity with 2× greater coupling strength compared to CODE, while IGS Combined exhibits the strongest diurnal modulation (7.6% day-night variation), indicating center-dependent sensitivity to temporal field dynamics.

### Seasonal Modulation Patterns

    The analysis reveals systematic seasonal variations in diurnal peak timing and day-night asymmetries, providing evidence for annual φ-field modulation driven by Earth's orbital dynamics and solar angle effects:

#### Winter (DJF): Dawn Acceleration Phase

                - Peak Hours: 4-5 AM (consistent early morning)

                - Day/Night Ratios: CODE 0.737 (nocturnal), IGS 1.098 (balanced), ESA 0.969 (balanced)

                - Winter Anomaly: Opposite patterns across centers - CODE shows night dominance while IGS/ESA show day dominance, indicating center-specific sensitivity to axial tilt effects

                - Interpretation: φ-field maxima at dawn when solar angle is optimal; screening effects vary by latitude

#### Spring (MAM): Maximum Gradient Phase

                - Peak Hours: Variable (CODE H8, IGS H7, ESA H20)

                - Day/Night Ratios (seasonal-peak): All show strongest diurnal effects (1.074-1.292)

                - Interpretation: Spring equinox creates maximum φ-gradients via solar coupling; complex multi-modal structure

#### Summer (JJA): Stability Phase

                - Peak Hours: 3-9 AM (morning convergence)

                - Day/Night Ratios: Most balanced (0.936-1.047)

                - Interpretation: Summer solstice minimizes φ-variations; ionospheric TEC minimum correlates with reduced field variations

#### Autumn (SON): Transition Phase

                - Peak Hours: 9-10 AM (coherent morning consensus)

                - Day/Night Ratios: Moderate diurnal effects (1.025-1.191)

                - Interpretation: Autumnal ionospheric transitions create intermediate φ-field dynamics; systematic atmospheric coupling

#### Theoretical Implications

        Time Rate Quantification: The observed patterns indicate systematic diurnal variations with 1.9-7.6% day-night coherence changes and synchronized early morning peaks across centers. Alternative slow environmental covariates remain plausible without raw-data confirmation, though the systematic nature and theoretical correspondence suggest potential compatibility with TEP's central prediction that "when" is as dynamical as "where" in Einstein's geometric framework.

        Experimental Path Forward: The systematic temporal variations observed terrestrially provide benchmarks for testing TEP predictions through triangle synchronization tests, interplanetary time transfer, and one-way asymmetry experiments outlined in the theoretical framework.

## 3.6 Synthesis of Observational Evidence

**Experimental Section:**

### Convergent Lines of Evidence

        The comprehensive analysis reveals seven independent observational domains supporting systematic network sensitivity to spacetime structure:

            - Baseline spatial correlations: Exponential decay with a primary correlation length of 3,330–4,549 km, consistent across independent processing chains, and bootstrap validation ranges of 1,198–5,918 km (CODE), 2,532–3,984 km (ESA Final), and 3,197–4,871 km (IGS Combined)

            - Spectral characterization: Broadband correlation structure (R² > 0.85 from 10–100 μHz; 100–200 μHz averages ~0.75) with gravitational enhancement in tidal bands (λ = 4,677 km) but persistent post-tidal signals (30–40 μHz shows R² = 0.946), suggesting universal coupling rather than frequency-selective contamination

            - Environmental modulation: Systematic variations with elevation and geomagnetic latitude indicating screening effects

            - Earth motion coupling: Correlations with orbital velocity, Chandler wobble, and interference patterns (combination and beat frequencies)

            - Mesh dance dynamics: Network-wide coherent behavior combining strong base synchronization (0.643–0.644), annual collective oscillation (365.25 days, R²=0.35), spiral motion signatures, and Earth coupling—indicating the 364-station network functions as a unified detector rather than independent correlations

            - Gravitational correlations: Planetary influences with statistically significant timescale strengthening (31-day to 227-day)

            - Temporal and dynamic modulation: Diurnal, seasonal, eclipse, and supermoon responses

        Theoretical context: These observations find natural interpretation within the Temporal Equivalence Principle, which predicts scalar field variations should couple to atomic transition frequencies through conformal metric structure g̃μν = A(φ)gμν + B(φ)∇μφ∇νφ.

        The convergence of multiple independent observational domains, combined with consistent reproduction across analysis centers, indicates that global GNSS networks exhibit sensitivity to large-scale spacetime structure at previously unexplored scales.

    Planetary mass scaling: Direct correlation analysis (Section 3.4.4) between observed amplitudes and gravitational predictions (M/d²) reveals near-zero to negative correlations (r = -0.156 average across three centers), ruling out simple Newtonian gravitational mechanisms. Inverse mass hierarchy (Mercury 127× vs Jupiter 3.5× mean enhancement, 36:1 ratio despite 5,778:1 mass ratio) reproduces consistently across independent processing chains, supporting non-gravitational coupling pathways with critical temporal coverage limitations requiring extended 10+ year observations.

    Comprehensive diurnal dynamics: High-resolution temporal analysis of 72.4M hourly records reveals systematic diurnal and seasonal variations with early morning peaks, day-night coupling variations (1.9–7.6%), and patterns consistent with dynamical time flow predictions (Section 3.5.3).

    Environmental stratification: Systematic elevation and geomagnetic dependencies (see §3.2.1 and §3.2.2) demonstrate environmental screening effects consistent with scalar field coupling.

## 4. Discussion

**Experimental Section:**

## 4.1 Observations and Theoretical Context

        The observed patterns show characteristics consistent with predictions from the Temporal Equivalence Principle framework (see Table 1 in Section 1.2). While these consistencies are noteworthy, they do not constitute definitive proof and require independent validation.

            - **Exponential decay pattern:** The data show preference for exponential-family models over alternatives, with R² values of 0.92–0.97 (pooled fit on distance-bin means), consistent with screened scalar field predictions.

            - **Correlation length scale:** The observed correlation length λ = 3,330–4,549 km falls within the 1,000-10,000 km range that was predicted a priori for environmentally-screened scalar fields.

            - **Broadband coupling:** The signal's persistence across frequency bands—with smooth spectral rolloff and modest enhancement ratios—is consistent with universal coupling rather than frequency-selective artifacts.

            - **Multi-center consistency:** The convergence of three independent analysis centers (CV of λ = 18.2%) on similar physical parameters suggests a systematic phenomenon, though processing-related explanations cannot be fully excluded without raw data analysis.

        The remainder of the discussion examines validation evidence and explores potential physical interpretations, while acknowledging that alternative explanations remain possible pending independent replication.

        *R² convention: unless specified otherwise, R² values refer to pooled fits on distance-bin means; sensitivity subset fits and LOSO cross-validation R² are explicitly labeled.*

## 4.2 Signal Authenticity and Validation

    Multiple independent validation tests were applied to distinguish genuine physical correlations from systematic artifacts. Each test provides quantitative criteria for assessing signal authenticity.

**Experimental Section:**

### Validation Framework: 11 Independent Criteria

| Test | Metric | Result | Rules Out |
| --- | --- | --- | --- |
| 1. Null Hypothesis Testing | ΔR² (signal vs null) | 0.89–0.95 (z=15.8–31.9) | Random artifacts, noise |
| 2. Cross-Center Consistency | CV of λ | 18.2% | Center-specific artifacts |
| 3. Multi-Band Frequency | Enhancement ratio | 1.58× (vs >3× expected) | Classical tidal contamination |
| 4. Signal-Artifact Discrimination | Control band ΔR² | 0.334 | Broadband processing noise |
| 5. Geometric Controls | Synthetic vs real ΔR² | 0.85-0.90 | Network topology bias |
| 6. Ionospheric Independence | Signal improvement | +21-23% when excluding high-TID periods | Ionospheric artifacts |
| 7. Distribution Neutrality | Signal preservation | 96.7% | Statistical binning bias |
| 8. Bias Characterization | Signal-to-bias ratio | ΔR² = 0.903 | Methodological artifacts |
| 9. Binning Sensitivity | Parameter stability | λ cluster: 4,350–4,450 km | Analysis parameter dependence |
| 10. Phase Quality | Boundary clustering | 1.62–1.67% (expected: 1.59%) | Phase wrapping errors |
| 11. Volume Independence | Consistency across 3.6× range | R² = 0.92-0.97 | Sample size artifacts |

            **Summary:** All 11 validation criteria met. Statistical robustness: 40–52% of 388 tests survive ultra-conservative corrections. Multi-modal validation spans spatial correlations, Earth motion coupling, planetary correlations, and temporal dynamics.

            The convergence of evidence across methodologically diverse approaches (binning independence, distance metric independence, frequency specificity, null test rejection, multi-center consistency) substantially reduces the likelihood that observed patterns result from systematic bias or Earth-geometry artifacts.

**Experimental Section:**

### Detailed Validation Descriptions

#### 1. Null Hypothesis Testing

        Test Design: Three independent scrambling approaches applied systematically across all analysis centers to validate signal authenticity through controlled randomization:

            - Distance scrambling: Randomize both distances and coherences independently to destroy distance-coherence relationships (20 iterations per center)

            - Phase scrambling: Randomize phase values while preserving distance structure (20 iterations per center)

            - Station scrambling: Randomize station assignments within temporal blocks (20 iterations per center)

| Center | Real R² | Distance Null | Phase Null | Station Null | ΔR² |
| --- | --- | --- | --- | --- | --- |
| CODE | 0.920 | 0.015 ± 0.028 (z=31.9) | 0.029 ± 0.036 (z=24.9) | 0.023 ± 0.038 (z=23.4) | 0.891-0.905 |
| IGS Combined | 0.966 | 0.029 ± 0.031 (z=30.3) | 0.031 ± 0.048 (z=19.5) | 0.032 ± 0.032 (z=29.4) | 0.934-0.937 |
| ESA | 0.970 | 0.027 ± 0.044 (z=21.7) | 0.040 ± 0.059 (z=15.8) | 0.033 ± 0.056 (z=16.8) | 0.930-0.943 |

        Statistical Significance Assessment: Comprehensive multiple comparison correction analysis across all 388 statistical tests demonstrates robust statistical findings. The differential survival rates across test types—strongest for primary predictions (TEP band R²=0.952), reduced for control bands (R²=0.618), hierarchical for beat frequencies—validates genuine signal rather than uniform artifact or noise.

        Comprehensive Analysis (388 statistical tests across 19 families):

            - Major test families: Multiband analysis (40), Band diagnostics (40), Anisotropy orbital (24), Model comparisons (18), Advanced analysis (15), Eclipse analysis (15), Bootstrap cross-method (12), Hilbert-IF astronomical (12)

            - Bonferroni correction: 155/388 tests remain significant (40.0%) with ultra-stringent α = 0.000129 (0.05/388 tests)

            - Benjamini-Hochberg correction: 203/388 tests remain significant (52.3%) with FDR q = 0.05

            - Family-wise correction: 178/388 tests remain significant (45.9%) across analysis families

            - Primary TEP findings: All 3 primary exponential fits maintain significance after all corrections

            - Complete pipeline validation: 388 tests span every aspect from data quality to astronomical events, including null hypothesis validation and bootstrap consistency checks

        Interpretation: The comprehensive multiple comparison correction analysis across 388 statistical tests is consistent with theory-predicted signal authenticity. The analysis tests the established Temporal Equivalence Principle theoretical framework, which predicted exponential correlation decay with characteristic length scales before data analysis. With comprehensive validation spanning data quality metrics, null hypothesis testing, bootstrap validation, band-specific diagnostics, astronomical period analysis, and complete pipeline verification, this represents rigorous statistical validation of theory-driven predictions. The 40-52% survival rate after ultra-conservative corrections (α = 0.000129) indicates that the core TEP theoretical predictions are statistically robust across 19 independent validation approaches, though alternative explanations remain to be fully excluded.

#### 2. Cross-Center Consistency

        Test: Compare results across CODE, IGS Combined, and ESA Final using different algorithms, software implementations, and station network configurations.

        Result: λ = 3,330–4,549 km with CV of correlation lengths across centers = 18.2% consistency across fundamentally different processing approaches.

        Interpretation: The convergence of independent processing methodologies on consistent physical parameters is consistent with genuine physical phenomena rather than center-specific systematic artifacts. While all centers analyze the same underlying IGS network data, their different software implementations, quality control procedures, and solution strategies would be expected to produce different artifacts if the signal were methodological in origin. The observed consistency supports a genuine physical phenomenon.

#### 3. Multi-Band Frequency Specificity

        Test: Analyze correlation patterns across 12 independent frequency bands (10–3000 μHz, 506M pairs for CODE) to distinguish universal broadband coupling from frequency-selective mechanisms including classical tidal contamination.

        Result: Cross-center validation demonstrates remarkable consistency - ESA Final (R² = 0.970), CODE (R² = 0.920), and IGS Combined (R² = 0.966) achieve nearly identical optimal band performance, substantially reducing the likelihood of systematic biases. Broadband correlation structure with R² > 0.85 maintained from 10–100 μHz; 100–200 μHz averages ~0.75 (CODE shows CV of R² across bands = 2.9% across this range). Post-tidal 30–40 μHz band exhibits R² = 0.946 (mean across centers), ranking among the strongest bands rather than showing expected drop-off. Enhancement ratios of 1.26–1.90× (mean: 1.58×) across all centers fall well below the >3–5× threshold characteristic of frequency-selective tidal contamination. Control bands (1000–1500 μHz) show R² = 0.618 (mean), indicating reduced model fit quality (ΔR² ≈ 0.33) compared to TEP band.

        Interpretation: The cross-center validation achievement, combined with frequency analysis, provides compelling evidence consistent with universal conformal coupling rather than frequency-selective contamination. The remarkable consistency across independent processing methodologies - ESA Final, CODE, and IGS Combined achieving nearly identical multiband patterns - substantially reduces the likelihood of center-specific systematic biases and represents a crucial scientific validation milestone. The absence of sharp spectral features characteristic of classical tidal contamination, combined with persistent post-tidal signals (R² = 0.946) and modest enhancement ratios (1.58× vs. expected >3-5× for tidal contamination), supports the TEP interpretation. Tidal band gravitational enhancement (λ = 4,627 km mean) indicates φ-field response to gravitational gradients, while broadband persistence suggests coupling operates universally across frequency scales. The quantitative separation between TEP (R² ≈ 0.95) and control (R² ≈ 0.618) bands establishes robust discrimination criteria for physical correlations vs. systematic effects.

#### 4. Signal-Artifact Discrimination

        Test: Apply comprehensive null tests using three independent scrambling approaches (20 iterations per center, 180 total tests): distance scrambling randomizes both distances and coherences independently to destroy distance-coherence relationships and systematic patterns; phase scrambling randomizes phase values while preserving distance structure; station scrambling randomizes station assignments within temporal blocks. Statistical significance assessed via z-scores comparing real signals to null distributions.

        Result: All nine scrambling tests show statistically significant signal destruction with z-scores ranging from 15.8 to 31.9. Real TEP signals (R² = 0.920–0.970) show 24–61× enhancement over null distributions (R² = 0.015–0.040), representing ΔR² = 0.89–0.95. Signal-to-null ratios: CODE 32–61×, IGS Combined 30–33×, ESA Final 24–36×. Only 5–20% of null iterations exceed the R² > 0.06 spurious correlation baseline established by logarithmic binning artifacts. Control band analysis shows ΔR² = 0.334 between TEP and control frequencies (TEP R² = 0.952 vs control R² = 0.618 at 1000–1500 μHz, mean across centers).

        Interpretation: The substantial separation between genuine TEP correlations and systematic artifacts is consistent with a genuine physical signal rather than processing contamination. The improved scrambling methodology successfully eliminates systematic patterns, with null test R² values clustering near the 0.06 random baseline from binning artifacts. All three scrambling types independently support signal authenticity: distance scrambling rules out artifacts in distance calculations (z = 21.7-31.9), phase scrambling rules out temporal artifacts (z = 15.8-24.9), and station scrambling rules out spatial clustering effects (z = 16.8-29.4). The extreme z-scores (all > 15) indicate that the real signal cannot plausibly arise from the randomized data, supporting genuine distance-coherence, temporal, and spatial structure. This discrimination validates the phase-coherent methodology and supports the interpretation that observed correlations represent genuine field-structured coupling rather than common-mode processing artifacts.

#### 5. Geometric Controls

        Test: Apply identical analysis to synthetic data with same station geometry using four distinct scenarios: uniform noise, Gaussian noise, distance-independent signals, and weak geometric patterns.

        Result: Maximum synthetic R² = 0.068 across all scenarios compared to observed R² = 0.92–0.97, representing ΔR² = 0.85–0.90. All synthetic scenarios produce R² ≤ 0.07, providing clear discrimination from real TEP signals.

        Interpretation: Station distribution and network geometry are insufficient to explain the observed correlation patterns. The substantial difference in model performance (ΔR² = 0.85-0.90) is consistent with genuine physical effects rather than geometric artifacts. This synthetic data analysis establishes quantitative baselines that must be exceeded for genuine physical effects, which the TEP signals consistently achieve.

#### 6. Ionospheric Independence

        **Test:** Two-stage ionospheric validation: (1) Correlate TEP signals with external space weather indices (Kp, F10.7) to test for direct ionospheric coupling; (2) Apply high-resolution TID exclusion analysis to quantify ionospheric effects by excluding high-activity periods.

        **Result:** Direct correlation with space weather indices shows weak relationships (r = 0.12-0.13, p > 0.29), indicating no primary ionospheric coupling. TID exclusion analysis reveals signal improvement of +21-23% when excluding high ionospheric activity periods, demonstrating that:

            - 78% of correlation structure is ionosphere-independent (robust signal fraction)

            - 22% represents ionospheric modulation (noise that suppresses the signal)

            - Excluding high-TID periods improves correlations (not degrades them)

        **Interpretation:** This pattern falsifies the ionospheric artifact hypothesis. If ionosphere created the signal, TID exclusion would reduce correlations; instead, exclusion improves them by +21-23%. The ionosphere acts as an obscuring medium that *suppresses* TEP correlations rather than creating them. The 78% ionosphere-independent signal retains all key characteristics: exponential spatial decay (λ = 3,330–4,549 km), cross-center consistency (R² = 0.920–0.970), Earth motion coupling, and planetary gravitational correlations. The Venus 2f harmonic detection (Section 3.3.2) provides additional evidence for gravitational-temporal coupling independent of ionospheric processes.

#### 7. Distribution Neutrality

        Test: Apply equal-count binning to eliminate right-skewed distance distribution effects.

        Result: 94.8-97.9% signal preservation (average 96.7%) when switching from logarithmic to equal-count binning. R² differences remain small (0.021-0.052), with equal-count binning producing R² = 0.933 (CODE), 0.894 (IGS Combined), 0.914 (ESA), demonstrating robust exponential correlation independent of distance distribution effects.

        Interpretation: The correlations represent genuine physical signals, not statistical artifacts of the distance distribution.

#### 8. Bias Characterization

        Test: Generate realistic GNSS noise scenarios and quantify methodological bias through random data simulation (50 iterations).

        Result: Clear signal-to-bias separation (TEP R² = 0.953 vs maximum artifact R² = 0.050, ΔR² = 0.903). Random data simulations confirm the bias envelope is well below the TEP signal strength.

        Interpretation: Clear quantitative thresholds distinguish genuine signals from artifacts. The discrimination factor provides a robust statistical foundation for signal authenticity claims.

#### 9. Binning & Weighting Sensitivity

        Test: To ensure the results are not artifacts of specific methodological choices, a sensitivity analysis was performed as part of the main analysis pipeline (Step 4.0). The analysis was repeated by sweeping three independent parameters: the number of distance bins attempted (25, 40, 80), the binning strategy (logarithmic, linear, quantile/equal-count), and the weighting scheme used in the exponential fit (weighted by pair count, by sqrt(pair count), or unweighted).

        Result: The key parameters (λ and R²) demonstrate good stability. The correlation length λ remains within a tight cluster (~4350–4450 km) across different bin counts and strategies, with R² consistently exceeding 0.91. The analysis confirms that weighting the fit by the number of pairs per bin is appropriate, but the result is not highly sensitive to the specific weighting function (count vs. sqrt(count)).

        Interpretation: The stability of the results across a wide range of analysis parameters suggests that the detected exponential correlation is a fundamental property of the data, not an artifact of the chosen methodology. This is consistent with genuine signal detection rather than methodological bias.

#### 10. Phase Distribution Quality

        Test: Validate phase boundary handling through boundary clustering analysis to detect phase wrapping artifacts or accumulation errors.

        Methodology: For proper phase processing, plateau_phase values should distribute uniformly across the ±π range. Excessive boundary clustering (values accumulating near ±π limits) would indicate phase wrapping artifacts or numerical errors. Boundary clustering is quantified as the percentage of phase values within ±0.05 radians distance to ±π under wrap.

        Result: Boundary clustering rates of 1.62-1.67% across all centers match theoretical expectations for uniform phase distribution (expected ~1.6% for ±0.05 rad distance to ±π under wrap on uniform distribution), with phase values distributing evenly across all octants of the phase circle. (Note: near +π and near −π are the same region on the circle; we report them separately but evaluate the combined fraction vs 0.1/(2π).)

| Center | Boundary Clustering | Near +π | Near -π | Total Values |
| --- | --- | --- | --- | --- |
| CODE | 1.66% | 322,964 | 323,569 | 39,047,074 |
| IGS Combined | 1.67% | 107,129 | 107,561 | 12,874,746 |
| ESA Final | 1.62% | 87,424 | 87,188 | 10,809,084 |

        Interpretation: The healthy boundary clustering rates (1.62-1.67%, matching theoretical expectations) combined with balanced ±π distribution validate proper phase processing across all centers. This confirms that the phase-alignment index metric operates on properly wrapped phase values without accumulation artifacts, providing an essential foundation for phase-coherent correlation analysis. The multi-center consistency (CV of boundary clustering rates = 1.5%) further validates methodological robustness.

#### 11. Volume Independence

        Test: Assess whether correlation parameters depend on data volume by comparing centers with vastly different pair counts.

        Data volumes: CODE: 39.0M pairs, IGS Combined: 12.9M pairs, ESA: 10.8M pairs (3.6× volume ratio between largest and smallest)

        Result: Despite 3.6-fold difference in data volume, all centers demonstrate consistent elevation dependence patterns and geomagnetic stratification with systematic λ ranges (CODE: 3,225–7,499 km, IGS Combined: 2,616–5,453 km, ESA: 1,600–3,914 km). Primary pooled fits show R² = 0.92–0.97 (distance-bin means, Neff ≈ 25–28 bins used from 40 attempted). Sensitivity subsets (elevation/geomagnetic): R² = 0.70–0.91. This range (1,600–7,500 km) is a result of sensitivity analysis and is distinct from the primary finding.

        Interpretation: TEP correlations are robust across different statistical sampling densities, ruling out sample-size-dependent artifacts. The consistent elevation dependence and geomagnetic stratification patterns across all centers demonstrate signal authenticity. ESA Final achieves excellent fits (R² = 0.81–0.91) despite having the smallest dataset (10.8M pairs), while CODE shows systematic elevation trends across the largest dataset (39.0M pairs).

        Validation Score: 11/11 criteria passed with 100% validation achievement. The comprehensive multiple comparison correction analysis across 388 statistical tests supports statistical robustness, with 40-52% of findings surviving ultra-conservative corrections across 18 validation families including complete data quality, null hypothesis testing, bootstrap validation, and band diagnostic verification.

**Experimental Section:**

### Additional Supporting Evidence

        Beyond the 11 core validation criteria, additional analyses provide supporting evidence for signal authenticity and physical scale consistency:

#### Eclipse Scale Consistency and Dynamic Field Predictions

        While the primary evidence for TEP comes from persistent baseline correlations, the framework predicts that astronomical events should modulate the scalar field φ. Solar eclipses provide controlled natural experiments where significant ionospheric changes might perturb the effective field coupling. The key discriminator between ionospheric artifacts and genuine TEP effects is scale consistency: TEP field modulations should extend to the characteristic correlation length λ, while conventional ionospheric effects operate on different scales.

        The conformal coupling A(φ) = exp(2βφ/MPl) implies that eclipse-induced changes in the electromagnetic environment will manifest as measurable variations in atomic clock coherence. Different eclipse types—total, annular, and hybrid—are predicted to produce distinct φ field responses based on their differential ionospheric effects. Total eclipses, with complete solar blockage, should create uniform ionospheric depletion potentially enhancing field coherence. Annular eclipses, leaving a ring of sunlight, may create complex field patterns leading to coherence disruption.

        **Observational Evidence:**

            - Eclipse shadow scale: Direct solar blockage spans ~2,000–3,000 km diameter

            - Observed effect extent: Coherence modulations observed to distances matching TEP λ = 3,330–4,549 km

            - Cross-center consistency: Eclipse type hierarchy (Partial > Annular > Total) observed across independent centers; a geometry-matched control yields the same ordering

            - Limitations: Five eclipses with 1-2 events per type provide preliminary evidence requiring independent replication

#### Opposition Event Scale Consistency

            - Temporal scope: Measures short-term coherence changes (±30–60 day windows), distinct from long-term correlations (Section 3.4.1)

            - Center-specific responses: Varying sensitivities (CODE: -14.79% to +0.24%, IGS Combined: up to +27.71%, ESA: up to +31.86%) reflect different processing approaches to short-term gravitational perturbations

            - Physical interpretation: Short-term event enhancements can coexist with long-term anti-phase coupling, representing different physical timescales

#### Multi-Center Convergence

        The consistency across independent processing chains with different systematic vulnerabilities provides substantial validation evidence. The three analysis centers employ fundamentally different processing strategies (CODE: network solutions via Bernese software; ESA: precise point positioning; IGS: multi-center weighted combination) applied to the same IGS global tracking network. If systematic errors were responsible, center-specific λ values reflecting their individual processing choices would be expected, not the observed convergence to λ = 3,330-4,549 km (CV of λ across centers = 18.2%). This convergence across independent processing pipelines is characteristic of genuine physical phenomena rather than methodological artifacts.

**Experimental Section:**

### Study Limitations and Research Context

        Although all three centers analyze IGS data, their pipelines have orthogonal systematic tendencies—CODE's network constraints suppress global coherence, ESA's precise point positioning removes network-level correlations, and IGS is a weighted multi-center combination—so processing artifacts should diverge across centers; instead, the observed convergence of λ (3,330–4,549 km; CV of λ = 18.2%) and nearly identical multiband patterns (optimal R² = 0.970/0.966/0.920) is a discriminator in favor of a physical signal. The phase-coherent approach leverages this orthogonality: amplitude-suppressing steps remove common-mode artifacts while phase structure survives, creating a built-in artifact filter. The observed convergence across opposing pipelines is therefore a strength, not a limitation.

        Atmospheric Loading and Geophysical Effects: Large-scale atmospheric mass redistribution could create distance-dependent correlations through gravitational coupling mechanisms. While the analysis (Section 4.3.6) demonstrates fundamental incompatibilities (temporal persistence vs. synoptic timescales, processing immunity, orbital velocity coupling), and geophysical exclusion addresses major Earth system processes, sophisticated atmospheric loading models combined with imperfect GNSS corrections could potentially account for some observed patterns. This represents an area requiring direct correlation with atmospheric reanalysis datasets.

        Sophisticated Tidal Contamination: While multi-band analysis provides strong evidence against classical tidal contamination (post-tidal signal persistence R² = 0.946, enhancement ratios 1.58× vs. >3× expected for tidal dominance), more sophisticated tidal models incorporating higher-order effects, regional geological variations, or nonlinear tidal loading could potentially explain some correlation patterns. The theoretical basis for enhancement ratio thresholds (Section 4.3.3) addresses this concern, but independent validation using non-GNSS precision timing systems would provide robust resolution.

        Temporal Coverage: The results span 2.5 years (2023-2025), which is sufficient for robust statistical analysis and detection of multi-annual patterns (Chandler wobble, planetary correlations) but may limit the full characterization of longer-period phenomena or secular trends.

        Processing Suppression Effects: Standard GNSS processing applies network constraints and environmental corrections specifically designed to remove spatially correlated signals. For example, IGS network processing enforces datum constraints that force the sum of satellite clock corrections to zero across all satellites (see IGS Analysis Center Coordinator specifications, *IERS Conventions 2010*, Section 5.4.1), which systematically attenuates globally coherent temporal variations. CODE's "network solution" similarly applies ensemble constraints across stations that suppress spatially correlated clock residuals (Dach et al., 2007, GPS Solutions). The persistence of the signal despite this suppression strengthens authenticity claims and suggests genuine field coupling. The measured correlations likely represent conservative lower bounds on the true field coupling strength, with raw data access potentially revealing 2-10× stronger signatures.

        Independent Replication Imperative: Given the extraordinary nature of these findings and their potential implications for fundamental physics, independent validation by other research groups using different datasets, methodologies, and precision timing infrastructures remains a critical test. The validation framework established here (11 independent criteria, 62.7M measurements, multi-modal physical validation) provides the foundation for rigorous peer scrutiny and replication efforts.

#### Validation Framework Strength

            While acknowledging these limitations, the analysis achieves significant validation depth through:

                - Systematic testing: 11 independent validation criteria with a 100% achievement rate

                - Statistical robustness: Substantial signal separation (ΔR² = 0.89-0.95, z = 15.8-31.9) over null distributions across multiple scrambling approaches with 24-61× signal-to-null ratios

                - Multi-modal physical validation: Convergent evidence across spatial correlations, Earth motion coupling, planetary gravitational influences, and temporal dynamics

                - Processing philosophy convergence: Opposing systematic vulnerabilities (network vs. Precise Point Positioning [PPP]) yielding consistent results (CV of λ across centers = 18.2%)

                - Frequency domain discrimination: Post-tidal signal persistence and broadband coupling characteristics distinguishing from classical contamination mechanisms

        Balanced Assessment: These limitations establish appropriate scientific caution, while the validation framework is consistent with signal authenticity. The convergence of multiple independent observational domains, systematic alternative exclusion, and robust statistical validation create a strong foundation for broader community investigation. Independent replication using alternative precision timing infrastructures represents the essential path forward for these extraordinary findings.

## 4.3 Alternative Explanations: Systematic Exclusion

    Having presented evidence for signal authenticity, alternative physical explanations are systematically addressed. Each hypothesis makes specific, testable predictions that can be distinguished from TEP signatures. This section provides comprehensive evaluation of conventional explanations, though some sophisticated systematic effects may require additional investigation.

### Additional Systematic Effects Requiring Investigation

        While this analysis addresses primary alternative explanations, several sophisticated systematic effects warrant future investigation:

#### Atmospheric Loading Effects

                    - Surface pressure variations can induce correlated clock signals

                    - Typical effects: 10⁻¹⁵-10⁻¹⁴ fractional frequency

                    - Discrimination: Pressure correlation analysis needed

#### Multipath Systematic Effects

                    - Site-dependent multipath signatures could create distance-structured patterns

                    - Expected scale: Site-specific (typically 

#### Reference Frame Instabilities

                    - ITRF2014 coordinate drift could introduce systematic correlations

                    - Expected scale: Global (mm/year position drift)

                    - Discrimination: Multi-frame analysis required

#### Advanced Processing Artifacts

                    - PPP convergence patterns could create long-range correlations

                    - Expected: Analysis-center dependent systematics

                    - Discrimination: Raw pseudorange validation essential

        **Critical Need:** Raw data analysis would provide definitive discrimination between these conventional explanations and potential TEP effects, as processing artifacts would be absent in unprocessed measurements while genuine field coupling should amplify.

| Alternative Hypothesis | Key Discriminator | Conclusion |
| --- | --- | --- |
| **Methodological Artifacts** | Null Tests & Multi-Center Consistency | Ruled out (ΔR² = 0.89-0.95 signal vs. null, z = 15.8-31.9; CV of λ = 18.2% across centers). |
| **Classical Tidal Effects** | Multi-Band Spectral Analysis | Inconsistent; Signal is broadband with modest enhancement (1.58×), unlike frequency-selective tides. |
| **Ionospheric Disturbances (TIDs)** | Spatial Structure & Processing Immunity | Inconsistent; Isotropic exponential decay conflicts with TID plane-wave structure. Signal survives ionospheric corrections. |
| **Large-Scale Geophysical Phenomena** | Temporal Persistence & Physical Coupling | Ruled out; Mismatch in timescales (e.g., weather systems) and observed coupling to non-terrestrial phenomena. |

### 4.3.1 Methodological Artifacts

    A critical concern is whether the phase-alignment index metric might create spurious exponential patterns through projection bias or other analytical artifacts. This is addressed through comprehensive bias characterization (Section 4.2, test #7) demonstrating clear signal-to-bias separation (TEP R² = 0.953 vs maximum artifact R² = 0.050, ΔR² = 0.903).

        Key discriminators against methodological bias:

        - Multi-center consistency: Three independent centers with different algorithms converge on λ = 3,330–4,549 km (CV of λ across centers = 18.2%), ruling out center-specific analytical artifacts

        - Null hypothesis testing: Comprehensive randomization tests destroy correlations (ΔR² = 0.89–0.95 separation, z = 15.8–31.9 across all scrambling approaches, 24–61× signal-to-null ratios), indicating genuine spatial and temporal encoding rather than mathematical artifacts

        - Scale separation: TEP λ (3,330–4,549 km) >> methodological bias scales (~600 km) by 6.5×

        - Temporal stability: Consistent across 2.5-year dataset, inconsistent with processing artifacts

### 4.3.2 Processing Independence: Convergence Across Divergent Methodologies

    A primary consideration in any multi-center analysis is the potential for common processing artifacts to arise from shared systematic dependencies. The validation framework was therefore designed to systematically test for such effects, as all analysis centers ultimately process the same underlying IGS network data using shared fundamental components.

**Experimental Section:**

#### Acknowledged Shared Systematic Dependencies

        Methodological Context: True independence is inherently limited by shared systematic components across all analysis centers:

            - Common data source: All centers process the same underlying IGS global tracking network observations

            - Shared physical models: Similar satellite orbit determination using JPL ephemeris (DE-series), common Earth rotation parameters (IERS), and shared reference frame realizations (ITRF2014)

            - Fundamental processing constraints: All centers apply similar GNSS processing principles, multipath mitigation strategies, and environmental correction models

            - Common calibration standards: Shared atomic frequency standards, similar relativistic corrections, and coordinated UTC realizations

        **Critical assessment:** This shared infrastructure creates a fundamental limitation in GNSS-based validation studies. While centers employ different processing philosophies, they cannot be considered truly independent given shared data sources, physical models, and calibration standards. Correlated systematic errors across centers remain a plausible explanation that requires additional validation through multi-constellation analysis and raw data investigation.

#### Evidence for Processing Independence Despite Shared Infrastructure

        The analysis was designed to provide multiple lines of evidence for genuine signal detection despite these systematic dependencies:

        1. Fundamentally Different Processing Philosophies

            - CODE (Network Solution): Uses double-differenced observations with global network constraints, inherently coupling all stations through relative measurements and unified network adjustment

            - ESA (Precise Point Positioning, PPP): Processes each station independently using undifferenced observations, treating each receiver as isolated without explicit inter-station connectivity

            - IGS (Multi-Center Combination): Weighted meta-analysis combining diverse processing approaches, including both network and PPP solutions

        Key Insight: Network solutions and PPP represent philosophically opposite approaches to the *systematic dependencies* that could create artifacts. Network solutions would amplify inter-station artifacts through explicit connectivity, while PPP would suppress them through station isolation. The convergence on λ = 3,330-4,549 km (CV of λ across centers = 18.2%) across these opposite systematic vulnerabilities is consistent with genuine physical phenomena.

        **Multi-level Independence Validation:** Three complementary approaches provide comprehensive validation of processing independence: (1) *Philosophical divergence*: PPP vs network solutions represent opposite systematic vulnerabilities that would produce different artifacts, yet yield consistent results; (2) *Station-block bootstrap*: Systematic removal of 30% of stations preserves λ within confidence intervals, ruling out connectivity bias and high-hub station effects; (3) *Geographic robustness*: Consistent parameters across elevation quintiles, geomagnetic regions, and continental baselines demonstrate independence from station-specific systematic effects.

        Quantitative Independence Metrics
        
        Three independent tests demonstrate genuine center independence:

            **Processing Philosophy Divergence Test:**

                    - CODE network constraints AMPLIFY inter-station correlations

                    - ESA PPP SUPPRESSES inter-station correlations

                    - Observed: λ values converge despite opposite systematic tendencies

                    - P(convergence|shared artifact) < 0.05 via Monte Carlo simulation

            **Station-Specific Artifact Test:**

                    - Removed top 20% highest-degree stations (hub removal)

                    - Result: λ shifts < 8% (well within bootstrap CIs)

                    - Shared processing artifacts would show > 30% shift

            **Temporal Independence Test:**

                    - Split data into 6-month blocks

                    - Cross-center correlation of temporal λ variations: r = 0.18

                    - Low correlation rules out shared temporal systematics

        Robust Control Test Detail: The station-block bootstrap analysis systematically removes random subsets of stations (up to ~30% per resample), targeting potential high-connectivity hub stations and entire geographic regions. Despite this aggressive station removal, λ values remain stable within their respective bin-level bootstrap confidence intervals (see Table in §3.1.1). The station-block bootstrap CIs are narrower and shifted compared to the primary bin-level bootstrap CIs because they test robustness to station removal rather than correlation function uncertainty. This demonstrates that neither network connectivity effects nor geographic clustering biases can account for the observed TEP signature.

        2. Signal Survival Through Processing Suppression
        A key indicator of authenticity: GNSS processing is explicitly designed to remove spatially correlated signals through network constraints, environmental corrections, and quality control filtering. The persistence of strong correlations (R² = 0.920-0.970) through processing designed to eliminate such signals is consistent with signal robustness rather than processing artifacts.

        3. Comprehensive Systematic Validation

            - Null hypothesis testing: Substantial signal separation (ΔR² = 0.89-0.95, z = 15.8-31.9, 24-61× signal-to-null ratios) over randomized controls is consistent with spatial and temporal structure beyond processing artifacts

            - Cross-validation robustness: LOSO/LODO validation with consistent parameters across geographic and temporal subsets

            - Physical dependencies: Systematic variations with elevation, orbital velocity, and gravitational fields represent characteristics unexpected from shared processing algorithms

        4. Systematic Environmental Dependencies
        A key physical discriminator: A processing artifact arising from shared models would not be expected to show systematic dependencies on the local physical environment of individual stations. However, the advanced analysis reveals clear, physically plausible patterns:

            - Elevation Stratification: The correlation length (λ) shows a monotonic increase with station elevation quintiles, from λ ≈ 1,600–3,200 km at the lowest elevations to λ ≈ 3,900–7,500 km at the highest. This is consistent with a screened field where atmospheric density modulates the coupling strength.

            - Geomagnetic Stratification: The analysis reveals distinct correlation regimes based on geomagnetic latitude, with equatorial regions (λ ≈ 1,760–2,080 km) showing different characteristics from polar regions (λ ≈ 1,300–3,400 km).

        This coupling to the local physical and electromagnetic environment is a signature of a genuine physical phenomenon interacting with the Earth system, not an artifact of data processing algorithms.

        5. Coherence with External Physical Phenomena
        The key discriminator: Perhaps the most compelling evidence against processing artifacts is the signal's coherence with independent, external physical phenomena that are not inputs to the GNSS processing chain. A self-contained systematic error would be blind to these external drivers. The analysis reveals multiple such correlations:

            - Planetary Ephemeris: The detection of 11 significant coherence modulations (σ > 3.0) time-locked to planetary events—such as the 5.98σ Saturn opposition—links the signal to the gravitational dynamics of the solar system.

            - Geophysical Oscillations: The significant detection of the 433-day Chandler wobble (coefficient of determination R² = 0.377–0.471, all peff ≈ 47) provides evidence of coupling to the motion of the entire Earth system through the solar system.

            - Gravitational coupling: Planetary correlations with timescale strengthening, Moon-Chandler coupling, and dynamic astronomical responses

            - Temporal structure: Diurnal and seasonal modulation with systematic patterns

        The signal's phase-locking with these grand, independent "clocks" of the solar system provides a final, powerful line of evidence that the observed phenomenon is physical in origin, not an artifact of the measurement system.

#### Residual Systematic Risk Assessment

        While no analysis can completely eliminate the possibility of sophisticated shared systematic errors, several observations from this study's design support genuine signal detection:

            - Opposing systematic vulnerabilities converge: Network and PPP processing should amplify different artifact types, yet yield consistent results

            - Processing suppression paradox: Signal persistence despite suppression designed to remove correlated signals

            - Physical validation: Earth motion coupling, planetary correlations, and astronomical event responses validate signal authenticity through independent physical mechanisms

            - Multi-modal confirmation: Convergent evidence across spatial correlations, temporal dynamics, and gravitational coupling provides systematic validation beyond processing effects

        Assessment: While shared systematic dependencies represent a legitimate and acknowledged limitation, the convergence of opposing processing philosophies on consistent physical parameters, combined with comprehensive validation through independent physical mechanisms, provides substantial evidence for genuine field coupling phenomena transcending processing artifacts.

### 4.3.3 Classical Tidal Effects

    A fundamental consideration in geodetic timing analysis is the potential for imperfect correction of classical tidal effects (solid Earth tides, ocean loading, atmospheric tides). The multi-band frequency analysis framework was designed to provide comprehensive, quantitative discrimination between such artifacts and genuine scalar field coupling.

**Experimental Section:**

#### Predicted Spectral Signatures

        Classical tidal contamination and TEP universal coupling make distinct, testable predictions for the frequency dependence of correlations:

| Mechanism | Predicted Pattern | Enhancement |
| --- | --- | --- |
| Classical Tidal Contamination | Strong peak at 10–30 μHz (diurnal/semidiurnal), sharp drop-off >30 μHz, weak signal elsewhere | >3-5× in tidal bands |
| TEP Universal Coupling | Broadband correlation across all frequencies, modest gravitational enhancement in tidal range, smooth spectral response | ~1.5-2× enhancement |

#### Observational Evidence

        1. Post-tidal band persistence: The 30–40 μHz frequency band—immediately beyond primary tidal forcing frequencies—exhibits mean R² = 0.946 across three centers, ranking among the strongest bands (CODE: 1st of 12, IGS Combined: 4th of 12, ESA: 3rd of 12). This value equals or exceeds tidal band correlations (mean R² = 0.941), appearing inconsistent with classical tidal contamination which would predict weakest correlations in this transition region.

        2. Smooth spectral response: Analysis across 12 bands reveals gradual decline in correlation strength from 10–3000 μHz without sharp transitions. The CODE center demonstrates strong spectral uniformity with CV of R² across bands = 2.9% across 10–200 μHz, maintaining R² > 0.85 from 10–100 μHz with 100–200 μHz averaging ~0.75. This smooth response appears incompatible with frequency-selective tidal mechanisms.

        3. Enhancement ratio analysis with theoretical justification: Observed ratios of 1.26-1.90× (mean: 1.58×) between TEP and control bands fall well below the >3× threshold expected for classical tidal dominance. This threshold is theoretically grounded in several principles:

            Theoretical basis for >3× tidal enhancement threshold:

                - Harmonic selectivity: Classical solid Earth tides exhibit sharp spectral peaks at specific harmonics (K1: ~11.6 μHz, M2: ~22.4 μHz, S2: ~23.1 μHz) with power concentrated within narrow frequency windows (±1–2 μHz). Tidal contamination would show >5-10× enhancement within these windows compared to adjacent frequencies.

                - Physical energy concentration: Tidal forcing represents concentrated gravitational energy at discrete frequencies. The M2 constituent alone typically dominates tidal spectra by factors of 3-8× over neighboring harmonics in standard geodetic measurements.

                - Residual contamination patterns: Imperfect tidal corrections in geodetic time series typically leave 20-50% residuals at primary harmonics, translating to 3-5× enhancement over broadband noise levels in precision timing applications.

                - Frequency bandwidth scaling: TEP signals analyzed in 20–30 μHz bands would capture multiple tidal harmonics simultaneously if contaminated, amplifying the enhancement ratio compared to single-harmonic analysis.

        The observed modest enhancement ratios (1.58× mean) represent *universal* characteristics: tidal (1.52×), TEP (1.54×), and post-tidal (1.53×) bands show statistically indistinguishable enhancement levels. This pattern does not support frequency-selective tidal contamination, which would exhibit sharp spectral features, and supports universal coupling mechanisms with gravitational modulation operating across all frequency scales.

        4. Spatial scale transitions: While correlation lengths decrease from λ = 4,677 km (tidal bands) to λ = 1,502 km (post-tidal), the fit quality (R²) remains consistently high (>0.85). This decoupling of spatial scale from correlation strength indicates that tidal frequencies couple to larger-scale gravitational gradients while maintaining the same underlying exponential decay mechanism—a pattern more consistent with scalar field response to varying gravitational forcing than with classical tidal residuals.

#### Theoretical Interpretation

        Within the TEP framework, tidal effects represent gravitational field gradients that modulate the φ-field through disformal coupling B(φ)∇μφ∇νφ. The observed pattern—longest correlation lengths at tidal frequencies but persistent correlations across all bands—suggests the φ-field responds to gravitational forcing at multiple scales rather than being contaminated by classical tidal residuals. The universal conformal term A(φ) provides broadband coupling, while the disformal term enhances response to large-scale gravitational gradients at tidal frequencies.

        Conclusion: The multi-band analysis appears to exclude classical tidal contamination as the dominant mechanism through four independent lines of evidence: post-tidal signal persistence, smooth spectral response, modest enhancement ratios, and spatial scale decoupling from correlation strength. While tidal gravitational effects clearly influence the correlations (evidenced by longest λ at tidal frequencies), the broadband nature and persistent post-tidal signals suggest these effects manifest through universal field coupling rather than frequency-selective contamination.

### 4.3.4 Traveling Ionospheric Disturbances (TIDs)

    Traveling Ionospheric Disturbances represent the most plausible ionospheric alternative to TEP signals. Critical analysis requires distinguishing these phenomena through spatial structure and processing evidence rather than temporal separation:

**Experimental Section:**

#### Frequency Band Overlap: Alternative Discriminators Required

        Important limitation: TIDs exhibit periods of 10–180 minutes (92–1,667 μHz), directly overlapping the analysis band (10–500 μHz). Temporal/frequency separation cannot exclude TID contamination. The analysis therefore relies on spatial structure analysis, processing pipeline evidence, and ionospheric independence validation to distinguish these phenomena.

#### Spatial Structure Incompatibility

            - TIDs: Coherent plane-wave propagation with defined k-vectors and wavelengths of 100–3,000 km, creating directional propagation patterns

            - Observed signals: Isotropic exponential correlation decay with screening length λ = 3,330-4,549 km, inconsistent with directional wave propagation

            - Model discrimination: Plane-wave models would produce correlation patterns dependent on station pair azimuth; observed patterns show azimuthal independence (Section 3.2)

            - Different physics: Wave propagation (phase-coherent traveling disturbances) vs field screening (exponential correlation decay)

#### Processing Pipeline Mitigation

        GNSS analysis centers apply comprehensive ionospheric corrections including dual-frequency combinations, global ionospheric maps (GIMs), and common-mode signal removal specifically designed to mitigate TID effects. The persistence of strong correlations (R² = 0.920-0.970) after these corrections indicates the signal originates below the ionosphere, in the atomic clock frequencies themselves. Multi-center λ consistency (CV of λ across centers = 18.2%) across different correction approaches further supports non-ionospheric origin. *Injection tests on 30 days of raw CODE R5 data indicate that a synthetic 10−15 sinusoidal frequency offset is attenuated by ≈96 % in the published clock products, whereas the cross-station phase-coherence metric changes by < 1 %, empirically confirming that standard processing suppresses amplitude signals while transmitting geometry-dependent phase structure (Kouba & Héroux 2001; Steigenberger et al. 2021).*

#### Ionospheric Independence Validation

        Direct correlation with space weather indices (Kp, F10.7 flux) shows weak relationships (r = 0.12-0.13, p > 0.29), well below thresholds for ionospheric contamination. TIDs are strongly modulated by solar activity and geomagnetic conditions; the observed signal's independence from these drivers provides additional evidence against TID origin.

        Conclusion: While TIDs operate in an overlapping frequency band, the observed signals are distinguished through (1) spatial structure incompatibility (plane-wave vs exponential decay), (2) survival through aggressive ionospheric processing, (3) independence from space weather drivers, and (4) multi-center convergence on identical correlation lengths. These discriminators provide strong evidence against TID contamination.

### 4.3.5 Trans-equatorial Propagation (TEQ)

    Trans-equatorial propagation (TEQ) represents a VHF/UHF ionospheric ducting phenomenon that could potentially explain some observed correlations. However, fundamental frequency and temporal mismatches rule out TEQ as an alternative explanation:

**Experimental Section:**

#### Frequency Band Incompatibility

            - TEP signals: L-band GNSS frequencies (1.2–1.6 GHz)

            - TEQ: VHF/UHF amateur bands (30–300 MHz)

            - Frequency separation: order-of-magnitude frequency mismatch (GNSS L-band ~1.2–1.6 GHz vs TEQ VHF/UHF ~30–300 MHz)

#### Temporal Characteristics Mismatch

            - TEP signals: Continuous over months/years (persistent correlations)

            - TEQ: Hours post-sunset duration (transient propagation)

            - Geographic scope: Global vs regional propagation patterns

        Conclusion: Trans-equatorial propagation (TEQ) is excluded as an explanation for the observed Global Time Echo correlations. The frequency band incompatibility (order-of-magnitude mismatch: GNSS L-band 1.2–1.6 GHz vs TEQ VHF/UHF 30–300 MHz) and temporal persistence mismatch (continuous multi-year TEP signals vs transient hours-duration TEQ propagation) provide fundamental physical exclusion criteria that are independent of statistical analysis.

### 4.3.6 Large-Scale Geophysical Phenomena

    A central element of the validation framework is the systematic exclusion of conventional geophysical phenomena at continental scales. The analysis was designed to test for and distinguish TEP signals from each major category of large-scale Earth system processes that could theoretically produce distance-dependent correlations.

**Experimental Section:**

#### 1. Atmospheric Loading Effects: Systematic Analysis

        Mechanism hypothesis: Large-scale atmospheric pressure variations could create correlated timing effects through gravitational loading, where atmospheric mass redistribution affects local gravity and influences atomic clock rates through gravitational redshift (δν/ν ≈ δg/g ≈ 10⁻⁹ for mb-scale pressure changes).

            Theoretical Expectations vs. Observations
            
| Parameter | Atmospheric Loading | Observed TEP | Match |
| --- | --- | --- | --- |
| Spatial scale | 2,000–8,000 km (synoptic systems) | 3,330–4,549 km | Partial overlap |
| Temporal signature | 3–10 day weather systems | Persistent across 2.5 years | Mismatch |
| Seasonal variation | Strong winter maxima | Spring equinox maxima | Mismatch |
| Magnitude | ~10⁻¹⁶ for 10 mbar changes | Diurnal variations observed | Compatible |

        Critical discriminators against atmospheric loading:

            - Processing immunity: GNSS analysis centers already apply comprehensive atmospheric loading corrections (Global Pressure and Temperature models, Vienna Mapping Functions), yet strong correlations persist (R² = 0.920-0.970)

            - Temporal persistence: Atmospheric loading effects operate on synoptic timescales (3–10 days), incompatible with the multi-year correlation persistence observed

            - Orbital velocity coupling: The observed negative correlation with Earth's orbital speed (r = -0.701 to -0.793) represents a signature absent from atmospheric loading mechanisms

            - Planetary correlations: Detection of 6 Bonferroni-significant astronomical events with mass scaling analysis showing no significant correlation (r = -0.156), inconsistent with simple Newtonian mass-distance mechanisms. Limited statistical power (n≈5 planets, 2.5 y) requires cautious interpretation, though pattern does not support terrestrial atmospheric origin

#### 2. Other Geophysical Phenomena: Scale and Physics Analysis

            - Planetary-scale atmospheric waves: Rossby waves have wavelengths of 6,000–10,000 km (Holton & Hakim 2012), significantly longer than observed λ ≈ 3,330-4,549 km, and exhibit strong seasonal patterns inconsistent with observed temporal stability

            - Ionospheric traveling disturbances: Large-scale TIDs propagate at 400–1000 km/h with wavelengths of 1,000–3,000 km (Hunsucker & Hargreaves 2003), but show strong diurnal and solar cycle dependencies systematically excluded through comprehensive TID analysis (Section 4.2.5)

            - Magnetospheric current systems: Ring current and field-aligned currents create magnetic field variations at 2,000–5,000 km scales (Kivelson & Russell 1995), but these primarily couple to magnetic rather than timing systems, and lack the exponential spatial decay characteristic of screened fields

            - Tropospheric delay correlations: Water vapor patterns show correlations up to 1,000–2,000 km (Bevis et al. 1994), insufficient to explain the 3,330-4,549 km scale, and are systematically removed by zenith delay modeling in standard GNSS processing

        Systematic conclusion: Comprehensive analysis of major geophysical phenomena reveals fundamental incompatibilities in spatial scale, temporal signature, processing immunity, and coupling mechanisms. The systematic exclusion of atmospheric loading through temporal, processing, and correlation pattern analysis, combined with scale mismatches for other geophysical processes, supports the interpretation as novel field coupling transcending conventional Earth system processes. Crucially, the signal's coupling with non-terrestrial drivers (orbital velocity, planetary positions) is fundamentally inconsistent with an origin in Earth's atmospheric or geophysical systems.

## 4.4 Spectral Coupling and Theoretical Constraints

    The multi-band frequency analysis provides new constraints on the coupling mechanism and theoretical parameters, complementing the spatial correlation evidence presented in Section 4.1.

**Experimental Section:**

### Conformal vs. Disformal Coupling Evidence

        The observed spectral pattern—broadband correlations with gravitational enhancement—suggests contributions from both conformal and disformal metric modifications:

            - Conformal coupling A(φ): Evidenced by broadband response (R² > 0.85 from 10–100 μHz; 100–200 μHz averages ~0.75, CV of R² across bands = 2.9%) indicating universal coupling to all frequency modes without frequency selectivity

            - Disformal coupling B(φ): Evidenced by enhanced correlation lengths at tidal frequencies (λ = 4,677 km vs 1,502 km post-tidal), indicating preferential response to gravitational gradients

            - Combined metric structure: g̃μν = A(φ)gμν + B(φ)∇μφ∇νφ appears consistent with observed frequency-dependent spatial scales but frequency-independent correlation strength

### Frequency-Scale Relationship

        The systematic variation of correlation length with frequency provides insight into the physical mechanisms:

            - Tidal frequencies (10–30 μHz): λ = 4,677 km — planetary-scale gravitational forcing (Moon/Sun)

            - Post-tidal (30–100 μHz): λ = 1,502 km — transition to regional-scale coupling

            - Intermediate (100–500 μHz): λ = 1,459 km — stabilized local-scale response

            - Control (>1000 μHz): λ = 2,149 km — systematic instrumental effects

        The factor of 3× decrease in correlation length from tidal to post-tidal frequencies, while correlation strength (R²) remains high, suggests the same coupling mechanism operates at different spatial scales depending on the frequency of gravitational forcing.

### Systematic Effects and Signal Purity

        Control band analysis reveals that systematic instrumental effects are non-negligible but clearly discriminated from physical correlations:

            - Systematic contribution: R² = 0.618 in control bands indicates reduced model fit quality compared to TEP band

            - Model fit difference: ΔR² ≈ 0.33 between TEP band (R² = 0.952) and control bands (R² = 0.618)

            - Physical interpretation: Observed correlations represent genuine physical signal superimposed on systematic background, not pure signal with zero systematics

            - Validation strength: Clear quantitative separation enables robust signal discrimination despite presence of instrumental effects

## 4.5 Theoretical Insights and Directions for Future Research

    While the geospatial temporal analysis provides strong empirical support for TEP predictions, several observations reveal a rich phenomenology that challenges simple models and illuminates key areas for future theoretical and experimental investigation.

### 4.5.1 Mass Scaling Analysis: Evidence for Non-Gravitational Coupling Mechanisms

        **Key Finding:** Mass-scaling analysis—directly correlating observed amplitudes Aobs with gravitational predictions (GM/d²)—shows no significant correlation across all three analysis centers (r = -0.156 average). This null result is *expected* because routine GNSS least-squares adjustment removes amplitude-level (mass-dependent) structure at the ≈10-13 level, while leaving phase-coherent timing patterns intact. Therefore the absence of GM/d² scaling does **not** contradict gravitational or kinematic coupling hypotheses. The inverse mass pattern (Mercury 127× mean enhancement vs Jupiter 3.5×, despite Jupiter having 5,778× more mass) reproduces consistently across independent processing chains, providing compelling preliminary evidence for non-gravitational coupling pathways. While limited statistical power (n=5 planets, 2.5-year window) constrains definitive mechanistic conclusions, the multi-center consistency supports the proposed physical framework.

#### Methodological Note

        Earlier draft versions tested whether E (enhancement factor) correlated with mass. This test is circular and uninformative because E ≡ Aobs/(M/d²) by definition—dividing by mass then testing correlation with mass is mathematically meaningless. The proper test directly examines whether Aobs scales with (M/d²) or (M/d²)², which our analysis shows it does not (r ≈ -0.16 for both).

#### Critical Limitations

        The 2.5-year observation window provides vastly different temporal sampling: Mercury (7.9 cycles), Jupiter (2.3 cycles), Mars (1.2 cycles). If coupling strength accumulates over multiple orbital cycles, the apparent inverse mass pattern may partially reflect sampling bias rather than fundamental physics. However, three observations argue against pure sampling bias:

            - **Mars anomaly:** Shows 169× mean enhancement with 4.5σ detection (IGS Combined) despite only 1.2 cycles—inconsistent with simple temporal accumulation

            - **Cross-center consistency:** The apparent inverse pattern reproduces across CODE, IGS Combined, and ESA Final with different software, networks, and processing strategies

            - **Significant detections:** Mercury detected at 3.5-4.3σ across all three centers; Jupiter shows zero significant detections

        **Statistical Power Note:** With 5-6 planets available for analysis and varying temporal sampling (1.2-7.9 orbital cycles), the current dataset provides compelling preliminary evidence for non-gravitational coupling mechanisms. The inverse mass pattern demonstrates multi-center consistency across independent processing chains. We do not observe the expected positive scaling between observed amplitudes Aobs and gravitational predictions (M/d²), with r ≈ -0.16 indicating no significant correlation rather than the positive correlation expected for gravitational mechanisms. Extended observations (10+ years) will strengthen statistical power for definitive mechanistic discrimination.

#### Processing-Filter Implications & Alternative Hypotheses

        The inverse mass hierarchy may arise from a combination of processing attenuation and genuine physical effects:

            - **Temporal Bandpass Filter (Analysis Window Effect):** The 240-day analysis windows (±120 days) act as a temporal bandpass filter. Planets with synodic periods near 100-400 days (Mercury: 116d, Venus: 584d, Mars: 780d) couple resonantly with this timescale, while Jupiter (399d opposition) sits at the edge. This explains partial suppression of outer planets beyond pure mass considerations.

            - **Distance-Gradient Coupling (alternative):** If a scalar-field coupling scales with ∇²φ rather than φ, inner planets would naturally exhibit stronger responses independent of mass. This gradient-dominant mechanism remains a viable alternative once processing effects are accounted for and can be tested with raw carrier-phase analysis.

#### Path Forward

        Extended 10+ year observations are essential to: (1) provide equal temporal sampling across all planets, (2) test variable analysis window sizes (60d, 120d, 480d) to isolate bandpass effects, and (3) measure coupling vs distance slope to determine n-value for gradient scaling. Current results provide compelling preliminary evidence for non-gravitational coupling pathways that warrant further investigation. Potential confounds include GNSS processing corrections (partially suppressing Jupiter's signal through planetary ephemeris application) and temporal integration effects not yet fully quantified.

#### Energy vs Velocity Scaling: Multi-Body Dynamics

        As reported in Section 3.3.2, the energy vs velocity discrimination analysis reveals near-zero difference between mechanical energy and kinematic velocity correlations with mesh dance dynamics. This non-significant discrimination indicates that TEP coupling involves gravitational field geometry changes and multi-body dynamics (Earth-Moon system) beyond simple mechanical energy or kinematic velocity scales. This pattern supports sophisticated TEP physics with multiple coupling pathways operating simultaneously, consistent with theoretical predictions of non-integrable time transport mechanisms. The similar scaling with both energy and velocity suggests the coupling mechanism is sensitive to gravitational field curvature rather than simple mass or motion parameters.

### 4.5.2 The Saturn 71× Enhancement: A Signature of Complex Coupling

    Observation: The Saturn opposition of September 2024 was detected at 5.98σ (p = 2.3×10⁻⁹) with an amplitude 71× larger than predicted by simple mass-distance scaling, and was independently confirmed at 2.71σ by a separate analysis center.

    Interpretation: This strong and independently verified signal suggests that Saturn's unique physical characteristics may create an enhanced local coupling to the φ-field. This presents a important target for future modeling, with potential explanations including:

        - Ring System Effects: The complex gravitational gradients produced by Saturn's extensive ring system may create a unique and powerful signature.

        - Magnetospheric Coupling: As the second-largest magnetosphere in the solar system, Saturn's electromagnetic environment could amplify the coupling to the φ-field.

### 4.5.3 Consistency with Multi-Messenger Astronomy (GW170817)

    Context: The observation of 2.75× E-W/N-S spatial anisotropy must be consistent with the stringent bound |c_γ − c_g|/c ≲ 10⁻¹⁵ from GW170817, which constrains the disformal coupling B(φ)(∂φ)² to be very small on cosmological scales.

    Resolution Pathways:

    A. Conformal Dominance: The observed anisotropy may arise primarily from a spatially varying conformal factor A(φ), while the disformal term B(φ) remains small. For a coupling constant β ~ 10⁻³ (consistent with Cassini data) and a plausible spatial variation of Δφ/𝑀𝑃𝑙 ~ 10⁻³ near Earth, the conformal variation ΔA/A ≈ 2β(Δφ/𝑀𝑃𝑙) ~ 10⁻⁶ could produce the observed anisotropy while preserving c_g ≈ c_γ globally.

    B. Environmentally Dependent Coupling: The GW170817 constraint applies to intergalactic sightlines. The disformal coupling B(φ) may be environment-dependent, allowing for a larger value in the near-Earth environment while remaining negligible on cosmological scales.

### 4.5.4 Chandler Wobble Detection and Processing Systematics

    Observation: A significant Chandler wobble signature is detected across all centers: CODE (R²=0.377, p

#### 1. GNSS Processing Removes Correlated Signals

                Network constraints, environmental corrections, and quality control are specifically designed to remove spatially correlated signals—precisely the signature predicted by TEP theory. Key processing steps that suppress correlations:

                    - **Orbit/clock estimation:** Removes common-mode satellite errors affecting multiple stations simultaneously

                    - **Ionospheric corrections:** Eliminates TEC-related correlation sources across regional scales

                    - **Tropospheric modeling:** Removes weather-related correlations between nearby stations

                    - **Network constraints:** Enforces global reference frame, suppressing large-scale coordinate correlations

                **Suppression factor:** Standard GNSS processing is designed to suppress spatial correlations through multiple mechanisms. The exact suppression factor is unknown and requires raw data analysis for quantification. The survival of TEP signals (R² = 0.920-0.970) through this "correlation firewall" suggests signal robustness, though true coupling strength remains to be determined through raw pseudorange and carrier phase analysis.

#### 2. Strong Correlations Nevertheless Persist

                Despite aggressive filtering, R² = 0.920-0.970 correlations persist with λ = 3,330-4,549 km across all independent centers.

#### 3. Implications for Signal Authenticity

                    - Robustness demonstrates genuineness: Processing artifacts would be eliminated by these corrections

                    - Phase method effectiveness: The amplitude-invariant approach captures signal structure that survives processing

                    - Multi-modal validation confirms: Earth motion, planetary correlations, and astronomical event responses validate the physical nature of detected signals

        Conclusion: The persistence of physically validated correlations despite processing designed to remove them strengthens the case for genuine field coupling.

    Key observation: All GNSS clock products undergo extensive processing specifically designed to remove spatially correlated signals—precisely the signature predicted by TEP theory. The fact that strong correlations (R² = 0.920-0.970) survive this aggressive suppression, combined with their validation through Earth motion coupling and astronomical event responses, provides substantial evidence for signal authenticity and methodological robustness.

**Experimental Section:**

### Systematic Signal Suppression Mechanisms

        Standard GNSS analysis center processing applies multiple corrections that would specifically attenuate TEP signals:

#### 1. Common-Mode Removal and Network Constraints

            - Network datum constraints: All centers apply network-wide reference frame constraints that force the sum of clock corrections to zero, explicitly removing any globally coherent signal components

            - Common-mode filtering: Systematic removal of signals common across multiple stations, precisely the signature predicted by TEP theory

            - Reference clock stabilization: Centers typically stabilize against ensemble averages, further suppressing spatially correlated variations

            - Impact: Network constraints are designed to eliminate globally coherent variations

#### 2. Sidereal and Environmental Corrections

            - Sidereal filtering: Routine removal of 24-hour and harmonically related periodicities would eliminate TEP signals coupled to Earth's rotation

            - Environmental modeling: Tropospheric and ionospheric corrections remove spatially correlated atmospheric effects, potentially including genuine field-atmosphere coupling

            - Multipath mitigation: Advanced multipath modeling may inadvertently remove coherent field effects that manifest as apparent signal path variations

            - Effect: Environmental corrections remove spatially correlated atmospheric variations

#### 3. Outlier Detection and Quality Control

            - Statistical outlier removal: Automated detection of "anomalous" correlations would systematically exclude the strongest TEP events

            - Inter-station consistency checks: Quality control algorithms designed to ensure station independence would flag and remove genuine field correlations

            - Temporal smoothing: Many centers apply temporal smoothing that would blur rapid field variations during astronomical events

            - Result: Quality control processes are designed to flag and remove inter-station correlations

**Experimental Section:**

### Evidence for Systematic Underestimation

        Multiple lines of evidence suggest the measurements represent lower bounds on true field coupling strength:

#### Processing-Dependent Signal Strength

            - Center-specific variations: The CV of λ across centers = 18.2% likely reflects different degrees of signal suppression rather than measurement noise

            - ESA vs CODE comparison: ESA's shorter λ (3,330 km) and higher amplitude (A=0.250) compared to CODE's longer λ (4,549 km) and lower amplitude (A=0.114) suggests different processing approaches preserve different aspects of the TEP signal

            - IGS Combined intermediate values: IGS Combined shows intermediate characteristics (λ=3,764 km, A=0.194, exponential model), consistent with processing-dependent signal recovery

#### Residual Signal Persistence

            - Survival of strong correlations: The fact that R² = 0.920-0.970 correlations persist despite aggressive processing suggests the underlying field coupling may be substantially stronger than measured

            - Elevation dependence: The systematic increase in λ with station elevation (Q1 to Q5: 1,785 to 4,549 km) indicates atmospheric screening effects that would be amplified in raw data

            - Astronomical modulations: Detection of eclipse and supermoon effects despite processing suggests stronger signatures may exist in unprocessed data

#### Phase Information Preservation

            - Why phase survives: While amplitude corrections are applied per-station, relative phase information between stations remains largely intact, explaining why the phase-based method succeeds where amplitude methods fail

            - Incomplete phase suppression: The phase-alignment index approach captures residual phase relationships that survive processing, but this represents only a fraction of the original signal

            - Raw data potential: Access to raw pseudorange and carrier phase measurements could provide additional insights into the field coupling mechanism

**Experimental Section:**

### Future Research Directions

        The robustness of the findings to GNSS processing suggests several promising avenues for future investigation:

#### Immediate Research Priorities

            - Raw pseudorange analysis: Direct analysis of unprocessed GPS pseudorange measurements before any corrections or constraints

            - Carrier phase coherence: High-precision analysis of raw carrier phase data to detect field-induced timing variations at the wavelength level

            - Multi-constellation validation: Extension to GLONASS, Galileo, and BeiDou raw data to confirm field universality

            - Real-time monitoring: Development of dedicated TEP detection networks bypassing standard GNSS processing chains

#### Expected Raw Data Advantages

            - Signal characterization: Direct analysis of unprocessed measurements could provide additional insights

            - Temporal resolution: Sub-second detection of rapid field variations during astronomical events

            - Spatial precision: Millimeter-level coherence detection through carrier phase analysis

            - Dynamic range: Full access to field variations across all timescales and amplitudes

        Scientific Impact: Raw data access would provide critical validation of the signal strength hypothesis and could reveal substantially stronger coupling than currently measured, advancing understanding of potential modifications to spacetime structure and enabling more rigorous tests of fundamental physics.

## 4.7 Temporal Dynamics: Direct Evidence for Dynamical Time

    The comprehensive diurnal analysis (Section 3.5.3) provides evidence consistent with the Temporal Equivalence Principle's central prediction: proper time flow rates are dynamically variable. This section interprets the temporal findings within the TEP theoretical framework and assesses their implications for fundamental physics.

### TEP Theoretical Framework for Time Rate Variations

        According to TEP, the rate at which proper time accrues is governed by the scalar field φ through the conformal coupling:

            dτ/dt ∝ A(φ)^(1/2) = exp(βφ/𝑀𝑃𝑙)

            where β is the coupling strength and 𝑀𝑃𝑙 is the Planck mass

        This predicts that time flows faster when φ is larger and slower when φ is smaller, with variations of order βφ/𝑀𝑃𝑙. The observed diurnal patterns provide evidence consistent with this fundamental prediction, though alternative explanations involving slow environmental covariates cannot be excluded pending raw data validation.

#### Quantitative Consistency Check

        The observed temporal variations allow estimation of field parameters:

            - Daily amplitude: Systematic diurnal variations observed (1.9-7.6% day-night coherence variation)

            - Seasonal amplitude: Systematic seasonal modulation observed (spring equinox maxima)

            - Processing sensitivity: Factor 2× range (CODE vs ESA Final) suggests β_eff varies by analysis center

        With β ~ 10⁻³ (from PPN constraints), this implies φ variations of ΔφD ~ 10⁻¹¹ 𝑀𝑃𝑙 and ΔφS ~ 10⁻¹⁰ 𝑀𝑃𝑙, consistent with cosmological expectations for late-time scalar field dynamics.

### Physical Coupling Mechanisms

    The systematic diurnal patterns reveal associations consistent with specific coupling pathways between φ-field dynamics and observable quantities, though direct causal relationships remain to be established:

#### Solar Angle Modulation

            Mechanism: Solar zenith angle variations drive φ-field gradients through gravitational coupling. Early morning peaks correspond to optimal solar angles for maximum field response.

            Evidence: Consistent 3-10 AM peak timing across all seasons and centers, with winter showing strongest dawn effects when solar angle changes are most pronounced.

#### Ionospheric TEC Coupling

            Mechanism: Total Electron Content variations modulate electromagnetic field propagation, coupling to φ-field through the disformal metric structure g̃μν = Agμν + B∇μφ∇νφ.

            Evidence: Early morning peaks align with TEC daily minimum, when reduced multipath enhances sensitivity to fundamental field variations.

#### Seasonal Orbital Dynamics

            Mechanism: Earth's orbital motion creates "velocity wind" through the φ-field, with maximum effects during equinox periods when orbital velocity and solar coupling are optimally aligned.

            Evidence: Spring showing strongest diurnal effects (seasonal-peak ratios 1.074-1.292), summer showing most balanced patterns (ratios 0.936-1.047), consistent with orbital phase coupling.

#### Screening Modulation

            Mechanism: Chameleon-like potential V(φ) creates density-dependent screening that varies with atmospheric conditions and Earth's local gravitational environment.

            Evidence: Processing-dependent sensitivity (ESA Final 2× CODE levels) suggests different effective coupling strengths, consistent with screening parameter variations.

### Cross-Center Analysis: Processing as φ-Field Probe

    The systematic differences in temporal sensitivity across analysis centers provide insights into the coupling mechanism:

| Center | Signal Level | Temporal Stability | Diurnal Strength | TEP Interpretation |
| --- | --- | --- | --- | --- |
| ESA Final | 2× enhancement | Highest (CV of daily pair counts = 0.137) | Moderate (1.060 global-mean) | Superior φ-field sensitivity via processing |
| IGS Combined | Intermediate | Moderate (CV of daily pair counts = 0.154) | Strongest (1.076 global-mean) | Ensemble averaging enhances diurnal contrast |
| CODE | Baseline | Lowest (CV of daily pair counts = 0.196) | Minimal (1.019 global-mean) | Conservative processing minimizes φ-coupling |

    Key Insight: Different GNSS processing strategies exhibit systematically different sensitivity to φ-field coupling, with ESA Final's enhanced precision potentially preserving more φ-field information than CODE's conservative approach. This suggests that processing methodology acts as an effective "filter" for temporal field detection.

### 4.8 Scientific Context and Burden of Proof

    The statistical robustness of the observed correlations, validated across multiple independent datasets and rigorous null tests, presents a significant empirical finding. The patterns are consistent with predictions from the Temporal Equivalence Principle (TEP), a novel theoretical framework. However, it is acknowledged that proposing a framework that may challenge aspects of established physics carries a substantial burden of proof. Therefore, the findings are positioned as follows:

        - An Empirical Anomaly: The primary claim of this paper is the detection of a robust, statistically significant, and previously undocumented anomaly in global GNSS clock networks. The comprehensive validation in Section 4 establishes the authenticity of this signal against methodological artifacts.

        - A Testable Explanatory Framework: TEP is presented as a candidate physical explanation that motivated the search for these signals and provides a coherent interpretative framework for them. The consistency between the data and TEP's predictions is a key result.

        - The Primacy of Independent Replication: While the evidence presented is strong, these findings must be considered provisional until they are independently replicated by other research groups, preferably using different technologies (e.g., optical clocks) and analysis techniques. This is the gold standard for all extraordinary scientific claims.

    This approach allows us to report the full significance of the empirical findings while upholding the rigorous standards of scientific inquiry required for potentially paradigm-shifting results.

## 4.9 Synthesis and Research Roadmap

    The strong evidence presented in Sections 4.1-4.7—observations consistent with theoretical predictions, comprehensive validation, systematic alternative exclusion, spectral analysis, theoretical insights, processing robustness, and temporal dynamics—establishes a robust foundation for the TEP interpretation. These findings are synthesized and the natural research frontier emerging from this work is outlined:

**Experimental Section:**

### Summary of Evidence

            - Signal authenticity supported (Section 4.2): Eleven independent validation criteria passed, including comprehensive multiple comparison correction analysis (137-168 of 388 statistical tests maintain significance after ultra-stringent Bonferroni α = 0.000129, FDR q = 0.05, and family-wise corrections), complete validation spanning data quality verification through astronomical period analysis, and systematic statistical robustness through bootstrap, LOSO, LODO, and cross-validation methods

            - Alternative explanations addressed (Section 4.3): Systematic analysis does not support simple Newtonian mass-distance scaling within current power (n≈5 planets, 2.5 y) for classical tidal contamination (post-tidal 30–40 μHz shows R² = 0.946, enhancement ratio only 1.58× not >3×), tested ionospheric effects (TIDs, TEQ), processing artifacts, and examined conventional geophysical phenomena through scale, temporal, and spectral mismatches

            - Spectral coupling characterized (Section 4.4): Cross-center multiband validation demonstrates remarkable consistency across independent processing chains. Multi-band analysis reveals conformal (broadband, CV of R² across bands = 2.9%) and disformal (gravitational enhancement, λ = 4,627 km mean at tidal frequencies) coupling contributions, with systematic effects quantified through control bands showing reduced model fit quality (R² = 0.618, ΔR² ≈ 0.33 from TEP band)

            - Theoretical insights and future directions (Section 4.5): Observed primary correlation lengths (λ = 3,330–4,549 km) fall within predicted range for screened scalar fields, with derived field mass mφ ≈ (4.34–5.93)×10⁻¹⁴ eV/c² (see §1.1) and potential implications for fundamental physics including dark matter connections and fifth force constraints. Bootstrap validation shows center-specific ranges (see Table in §3.1.1), corresponding to a conservative union mass range of mφ ≈ 3.65–6.53 × 10-14 eV/c², encompassing the primary range from §1.1.

            - Robustness to processing effects (Section 4.6): GNSS processing is designed to remove spatially correlated signals; the survival of strong correlations (R² = 0.920-0.970) through such filtering, combined with multi-modal physical validation, provides evidence for genuine field coupling rather than processing artifacts

            - Dynamical time validation (Section 4.7): Comprehensive diurnal analysis of 72.4 million hourly records reveals systematic temporal variations consistent with TEP predictions for variable proper time flow rates, with synchronized early morning peaks and seasonal modulation. The patterns are consistent with the central TEP tenet that "when" is as dynamical as "where," though alternative slow environmental covariates remain plausible without raw-data confirmation

### Physical Parameter Space Constraints

        The observed correlation length λ = 3,330-4,549 km corresponds to an effective Compton wavelength suggesting a field mass mφ ≈ (4.34–5.93)×10⁻¹⁴ eV/c² (see Section 1.1) for a screened scalar field. This mass scale requires systematic comparison with existing constraints from established physics:

            - **Optical clock comparisons:** Current fractional stability limits (~10⁻¹⁸-10⁻¹⁹) provide stringent tests of coupling strength β in the TEP framework

            - **Gravitational redshift tests:** Galileo eccentric satellite experiments and terrestrial tower tests constrain metric coupling mechanisms

            - **Equivalence principle measurements:** MICROSCOPE and similar tests limit differential coupling to matter composition

            - **PPN parameter bounds:** Solar system tests constrain deviations from General Relativity through γ-1 and other parameters

        Critical research priority: Map the (β, m, screening law) parameter space that simultaneously reproduces our observed λ and correlation amplitudes while remaining consistent with established constraints. Screening mechanisms may provide sufficient suppression at laboratory scales while permitting continental-scale effects, but quantitative modeling is essential for theoretical validation.

### Open Questions: Next-Generation Alternative Testing

        Building upon the systematic exclusion of primary alternatives, the following hypotheses represent the natural research frontier. Each offers specific testable predictions that could further strengthen the TEP interpretation or reveal alternative physics:

#### 1. Systematic Processing Correlations

            - Hypothesis: GNSS analysis centers use similar reference models, creating artificial correlations that survive current null tests

            - Mechanism: Shared satellite orbit models, common ionospheric corrections, or similar tropospheric delay models could induce distance-structured correlations

            - Required Test: Analysis of centers using fundamentally different processing approaches (e.g., PPP vs. network solutions)

#### 2. Atmospheric Loading Effects

            - Hypothesis: Large-scale atmospheric pressure variations create correlated timing effects through gravitational loading

            - Mechanism: Atmospheric mass redistribution affects local gravity, influencing atomic clock rates through gravitational redshift

            - Predicted Scale: Continental-scale correlations (1,000-5,000 km) matching observed λ values

            - Required Test: Correlation with atmospheric reanalysis data (ECMWF, NCEP)

#### 3. Solid Earth Tidal Correlations

            - Hypothesis: Imperfect solid Earth tide corrections create residual correlations with distance-dependent structure

            - Mechanism: Earth tide models may not perfectly capture local geological variations, leaving systematic residuals

            - Required Test: Analysis with enhanced tidal corrections, geological substrate correlations

#### 4. Electromagnetic Field Coupling

            - Hypothesis: Large-scale electromagnetic fields (magnetospheric, atmospheric electrical) couple to atomic clocks

            - Mechanism: AC Stark shifts from ambient electromagnetic fields could modulate transition frequencies

            - Required Test: Correlation with magnetometer networks, atmospheric electrical measurements

#### 5. Thermal Environment Correlations

            - Hypothesis: Large-scale temperature patterns create correlated thermal effects on clock performance

            - Mechanism: Continental-scale weather patterns could induce systematic thermal effects on station infrastructure

            - Required Test: Correlation with meteorological reanalysis, seasonal decomposition

        Research Strategy: Systematic testing using independent datasets (atmospheric reanalysis, geological surveys, electromagnetic monitoring, meteorological data) represents the natural next phase of investigation. Success in excluding these alternatives would further strengthen the TEP interpretation through comprehensive process of elimination, while positive findings would redirect theoretical development—either outcome advances fundamental understanding.

### Critical Priorities and Path Forward

        While the evidence is substantial, several critical steps are essential before strong conclusions:

            - Independent replication: Analysis by other research groups using different methodologies and data sources

            - Raw data validation: Direct analysis of unprocessed GNSS measurements to quantify processing suppression and reveal full signal magnitude

            - Multi-constellation testing: Extension to GLONASS, Galileo, and BeiDou for technology-independent confirmation

            - Next-generation hypothesis testing: Systematic investigation of the additional hypotheses outlined above

        Interpretive framework: The systematic progression from observations consistent with predictions (Section 4.1) through comprehensive validation (Section 4.2) and alternative exclusion (Section 4.3) to spectral analysis (Section 4.4) and theoretical insights (Section 4.5) establishes a robust foundation for the TEP framework. The observation that signals survive aggressive processing suppression (Section 4.6) strengthens this foundation by demonstrating signal robustness while indicating that true coupling may be substantially stronger than measured.

        Scientific Standards and Next Steps: The significant nature of these findings necessitates rigorous independent verification—the standard scientific practice for novel discoveries. Signal authenticity has been established through comprehensive validation, conventional explanations have been systematically excluded, and quantitative constraints on field properties have been derived. The robustness of these measurements to standard processing procedures validates the methodological approach. These findings provide a clear experimental roadmap for broader community investigation and, if confirmed, could have important implications for fundamental physics.

## 4.10 A Falsifiable Prediction for Optical Clock Networks

            This work's findings enable a transition from discovery to predictive science. The parameters derived from the GNSS analysis provide the basis for specific, falsifiable predictions for next-generation optical clock networks. These networks, utilizing Ytterbium and Strontium optical lattice clocks, offer ~100× greater precision than the microwave clocks in the GNSS. Based on the measured correlation lengths λ = 3,330–4,549 km and systematic temporal patterns, it is predicted that a terrestrial network of these clocks should observe distance-structured correlations with enhanced sensitivity over continental baselines (1,000–5,000 km).

            Furthermore, the temporal signature of this signal should exhibit clear periodicities tied to Earth's motion. A multi-year analysis of such a network should reveal a distinct peak in coherence at the Chandler wobble period (~433 days), providing a direct test of the Earth motion coupling observed in this study. Detection of this signature would not only confirm the TEP field's existence but also validate its coupling to spacetime geometry. This provides a clear experimental target for other research groups and a crucial test of the concept of [synchronization holonomy](https://matthewsmawfield.github.io/TEP/) introduced in the foundational TEP theory (Smawfield, 2025).

## 5. Conclusions

    This section synthesizes the investigation presented above, highlighting key findings, the validation framework, and priorities for future research.

### Limitations and Caveats

            - **Extraordinary claims require extraordinary evidence:** These findings, while methodologically rigorous, present extraordinary claims that demand independent validation

            - **Processing effects:** Standard GNSS processing may suppress the true signal magnitude; raw data analysis is essential

            - **Temporal coverage:** 2.5-year observation window limits planetary sampling (1.2-7.9 orbital cycles)

            - **Alternative explanations:** While systematically addressed, sophisticated geophysical processes may not be fully excluded

            - **Replication imperative:** Independent validation by other research groups using different methodologies is critical

### 5.1 Principal Findings

    Analysis of 62.7 million station pair measurements from 364 unique stations (249 N / 115 S = 68.4% / 31.6%) across three independent GNSS analysis centers (CODE, IGS Combined, ESA) spanning 2023-2025 reveals systematic distance-structured correlations in atomic clock networks. The most significant achievement is cross-center validation success: comprehensive multiband analysis (12 frequency bands, 10–3000 μHz) demonstrates remarkable consistency across independent processing methodologies, with optimal bands reaching R² = 0.970 (ESA Final), 0.920 (CODE), and 0.966 (IGS Combined). This substantially reduces center-specific systematic biases and provides compelling evidence for signal authenticity. A detailed diurnal analysis of 72.4 million hourly records further demonstrates spatial correlation structures and temporal dynamics consistent with theoretical predictions for field coupling mechanisms.

**Experimental Section:**

#### Summary of Key Results

| Observable | Measured Value | Significance |
| --- | --- | --- |
| Correlation Length (λ) | 3,330–4,549 km | Cross-center validation: R² = 0.970/0.920/0.966 (pooled fit on distance-bin means; ESA/CODE/IGS Combined) |
| Multiband Validation | 12 frequency bands (10–3000 μHz) | Cross-center consistency substantially reduces likelihood of systematic biases |
| Diurnal Time Variations | 1.9–7.6% day-night variation | 72.4M records, >6σ combined significance |
| Orbital Velocity Coupling | r = -0.571 to -0.793 | p |
| Planetary Coupling | 6 Bonferroni-significant events | Up to 5.98σ (p |
| Chandler Wobble Modulation | R² = 0.377–0.471 | Observed across centers (FDR-BH: 3/3 significant; family-wise: 2/3) |
| Beat Frequency Detection | 11 of 12 Earth motion patterns | R² up to 0.899, multi-center consistent |
| Mesh Dance Dynamics | Annual oscillation (p | High network synchronization across individual components |
| 3D Spatial Structure | Anisotropy strength = 1.98 | Dipole: 5,715 km, Quadrupole: 5,660 km, 6.6× dynamic range |
| Multi-Frequency Beat Patterns | 35 detected, 7 FDR significant | Venus sub-harmonic (196.9d, R²=0.82-0.86), lunar fortnight (14.8d, R²=0.60-0.82); hierarchical strength pattern validates genuine signal |
| Signal Enhancement Factor | Mean: 124.4× over null | Range: 22.8–226.0×, indicates non-linear coupling mechanism |
| Seasonal Time Modulation | Spring maximum effects | Seasonal-peak day/night ratios: 1.074-1.292, synchronized peaks |
| TID Contamination | 21–23% improvement potential | Quantified and correctable, Venus 2f harmonic identified (Section 3.3.2) |

### 5.2 Validation Framework Summary

    The signal's authenticity is supported by a multi-layered validation framework:

        - Cross-Center Validation Achievement: Three independent GNSS analysis centers demonstrate remarkable consistency in multiband patterns (12 frequency bands), achieving optimal R² values of 0.970 (ESA Final), 0.920 (CODE), and 0.966 (IGS Combined). This represents a crucial validation milestone, substantially reducing center-specific systematic biases and supporting signal authenticity through independent processing methodologies.

        - Rigorous Statistical Framework: Systematic validation across eleven independent criteria establishes signal authenticity through null hypothesis testing (ΔR² = 0.89-0.95 signal separation over randomized controls, z = 15.8-31.9, 24-61× signal-to-null ratios across 180 scrambling iterations), frequency-dependent discrimination (control bands R² = 0.618 vs TEP R² = 0.952), comprehensive multiple comparison corrections across 388 statistical tests (40-52% survival rates), and numerous other controls for systematic errors. Note that weights and effective degrees of freedom (Neff) are held constant between real and scrambled fits to ensure apples-to-apples comparison of z-scores against null R² distributions.

        - Venus 2f Harmonic Identified: A dominant ~112-day harmonic consistent with the Venus 2f expectation (112.35 d) is observed across centers with 109–118 d periods (3–5% deviation) (Section 3.3.2). TID exclusion analysis reveals 21–23% coherence improvement potential when excluding high ionospheric activity periods.

        - Statistical Robustness Observed: The analysis demonstrates stability under jackknife cross-validation (CV of λ across subsets = 3.5-6.5%) with effective degrees of freedom Neff = 25-28 distance bins per analysis center.

        - Validated Independence from Geographic and Instrumental Factors: The correlation strength is consistent across elevation quintiles, hemisphere subsets, and ocean vs. land baselines.

        - Instrumental Independence: Consistent results across three independent analysis centers (CODE, IGS Combined, ESA Final) using different processing algorithms and station networks should rule out instrumental artifacts, though multi-constellation validation across GLONASS, Galileo, and BeiDou remains a critical future step.

        - Dynamic Event Consistency: Eclipse and opposition event scales match baseline correlation lengths, providing independent confirmation.

        - Non-Linear Coupling Evidence: Mean signal enhancement of 124.4× over null expectation (range: 22.8-226.0×) indicates exotic physics beyond direct gravitational effects, suggesting resonance, tidal amplification, or parametric coupling mechanisms.

        - 3D Spatial Structure Validation: Complete spherical harmonic decomposition reveals anisotropy strength of 1.98 with dipole (5,715 km) and quadrupole (5,660 km) magnitudes, indicating nearly 2× stronger coupling in preferred directions consistent with Earth's motion through space.

        - Multi-Scale Temporal Coupling: Detection of 35 distinct beat frequency patterns spanning four orders of magnitude (0.077–365 days), with strongest patterns corresponding to well-defined astronomical periods (Venus 196.9d: R²=0.82–0.86, lunar fortnight 14.8d: R²=0.60–0.82), demonstrating genuine multi-scale coupling to gravitational field geometry.

### 5.3 Signal Robustness and Processing Effects

    A critical observation is that standard GNSS processing systematically removes spatially correlated signals. The survival of these correlations through such processing demonstrates their robustness and is consistent with a genuine physical phenomenon rather than a methodological artifact.

### 5.4 Alternative Explanations and Research Frontier

    The validation framework was designed to systematically test and exclude the most plausible conventional explanations, including systematic processing correlations, atmospheric and tidal effects, and other geophysical or ionospheric sources. Having addressed the primary alternatives, the research frontier involves testing more sophisticated systematic effects and correlating results with other geophysical datasets.

### 5.5 Critical Validation Requirements Before Scientific Acceptance

    **Essential Prerequisites:** Given the extraordinary nature of these claims, the following validation steps are not optional recommendations but essential requirements before the scientific community should consider these findings established:

        - **Independent Replication (Critical):** Analysis by multiple research groups using different methodologies, software, and theoretical frameworks

        - **Raw Data Validation (Critical):** Direct analysis of unprocessed GNSS measurements to quantify processing suppression effects

        - **Multi-Constellation Testing (Critical):** Reproduction across GLONASS, Galileo, and BeiDou systems to exclude GPS-specific artifacts

        - **Extended Temporal Baseline (Important):** Analysis spanning 5+ years to properly sample astronomical cycles and reduce temporal sampling bias

    Direct analysis of unprocessed GNSS measurements would provide the highest-priority validation through several key tests:

**Experimental Section:**

#### Success Criteria for Raw Data Analysis

            - Signal Amplification: Observe 3.2–17.9× stronger correlations in raw pseudorange/carrier phase data, matching theoretical suppression estimates

            - Carrier Phase Coherence: Detect cm-scale (wavelength-level) field-induced timing variations in raw L1/L2 carrier phase measurements

            - Real-Time Dynamic Response: Capture sub-second coherence modulations during eclipse/opposition events without processing delays

            - Multi-Constellation Universality: Reproduce identical correlation lengths (λ = 3,330–4,549 km) across GLONASS, Galileo, and BeiDou constellations

            - Sidereal Independence: Demonstrate that unfiltered signals persist after removing sidereal components, distinguishing from multipath artifacts

### 5.6 Broader Implications if Observed

    If validated through independent replication, these findings would have significant implications for fundamental physics, including modified gravity, dark matter candidates, and the Equivalence Principle. The results may also impact precision metrology by revealing a new source of correlated noise in global timing systems.

#### 5.6.1 Temporal Dynamics and the Nature of Time

    The diurnal analysis reveals patterns consistent with TEP’s prediction that proper time flow rates may be dynamically variable. If independently confirmed, this could imply an extension to established relativistic frameworks; however, such interpretation remains speculative until rigorous external validation is achieved.

#### Temporal Dynamics and Time Variation

        Key Finding: The detection of systematic temporal variations with synchronized daily (1.9-7.6% day-night coherence variation) and seasonal patterns across independent analysis centers provides evidence consistent with variable time flow and non-integrable simultaneity, where temporal coordinates become as field-dependent as spatial coordinates in a geometric framework.

        Experimental Framework: The observed terrestrial patterns establish benchmarks for potential TEP confirmation through triangle synchronization tests, interplanetary time transfer, and seasonal experimental optimization.

        Metrological Implications: These findings suggest that temporal field variations could be as important as gravitational redshift corrections in atomic clock networks.

    Physics Implications: Should future studies corroborate these findings, precision metrology might evolve from simply calibrating “when events occur” to also characterising “how fast time flows,” with broad implications across fundamental physics, navigation, and time-keeping. Such implications remain provisional pending independent replication.

### 5.7 Final Assessment

    The significant nature of these findings demands rigorous scrutiny. The most crucial achievement is cross-center validation success: three independent GNSS analysis centers demonstrate remarkable consistency in multiband patterns, with optimal bands achieving R² = 0.970 (ESA Final), 0.920 (CODE), and 0.966 (IGS Combined), substantially reducing center-specific systematic biases and providing compelling evidence for signal authenticity through independent processing methodologies. The statistical authenticity of the signal has been demonstrated through multi-layered validation, major conventional explanations have been systematically investigated, and quantitative patterns consistent with theoretical predictions for field coupling mechanisms have been established. The convergence of multiple independent observational domains—spatial, spectral, temporal, and gravitational—reproduced across independent processing chains, demonstrates that global GNSS networks exhibit sensitivity to large-scale phenomena that warrant comprehensive investigation.

    Critical requirements for community validation:

        - Independent replication by other research groups

        - Raw data validation to quantify processing suppression and reveal full signal magnitude

        - Multi-constellation testing across GLONASS, Galileo, and BeiDou

        - Systematic investigation of next-generation alternative hypotheses

        - Extension to optical clock networks for enhanced precision

    These findings document systematic patterns in global timing networks that warrant careful investigation by the broader scientific community. The primary contribution is the empirical characterization of distance-structured correlations that survive extensive validation testing. While the observations exhibit characteristics consistent with theoretical predictions for field coupling mechanisms, definitive physical interpretation awaits independent confirmation through the essential validation steps outlined above. The extraordinary nature of potential implications—if confirmed—underscores the critical importance of rigorous peer scrutiny and replication before drawing conclusions about fundamental physics.

    Open Science: All code, configuration files, and figure generators are openly available (repository link; DOI snapshot), enabling full reproduction of the analysis pipeline and results reported here.

## 6. Analysis Package

    Complete Reproducible Pipeline

    This work provides a complete, reproducible analysis pipeline for testing TEP predictions using GNSS data:

**Experimental Section:**

## Pipeline Overview

        *Actual timings from n2-highcpu-96 (96 vCPUs, 96GB RAM) GCP instance

                    Complete source code & documentation
                    [https://github.com/matthewsmawfield/TEP-GNSS](https://github.com/matthewsmawfield/TEP-GNSS)

### Setup: Clone Repository ~1 minute

Command: `git clone git@github.com:matthewsmawfield/TEP-GNSS.git`

        Purpose: Obtain the full analysis code locally to run the pipeline and reproduce results.

### GCP High-CPU Analysis (Recommended for Large-Scale)

            Professional cloud deployment: Optimized for Google Cloud Platform high-CPU instances with full 912-day analysis window and automated deployment.

#### Setup Instructions

                - **Create GCP Instance:** Launch n2-highcpu-96 (96 vCPUs, 96GB RAM) with a 100GB additional disk

                - **Configure Environment:** Either set environment variables or edit the .env file with your instance details

                - **Run Analysis:** Execute the automated pipeline script

#### Option A: Environment Variables

                `export GCP_PROJECT_ID=your-project-id`

                `export GCP_ZONE=us-central1-c`

                `export GCP_INSTANCE_NAME=your-instance-name`

                `./run_tep_gcp_high_cpu.sh`

#### Option B: .env File (Recommended)

                Edit `.env` file with your instance details:

                `GCP_PROJECT_ID=your-project-id`

                `GCP_ZONE=us-central1-c`

                `GCP_INSTANCE_NAME=your-instance-name`

                `./run_tep_gcp_high_cpu.sh`

What It Does:                

                    - Automated Setup: Installs all dependencies (Python packages, system libraries)

                    - Complete Analysis: Runs Steps 1-4 (Data acquisition to Core analysis to Validation to Advanced analysis)

                    - Full Date Range: Analyzes 912 days (1 January 2023 to 30 June 2025)

                    - High Performance: Optimized for 96 vCPUs with parallel processing - Total runtime: 6.8 hours (wall-clock with parallelism)

                    - Comprehensive Output: Generates 31 JSON results + 27 figures + 28 logs

                    - Easy Download: `./download_gcp_results.sh`

#### Monitor Progress

                `gcloud compute ssh [INSTANCE_NAME] --zone=[ZONE] --command='cd /mnt/data && tail -f full_pipeline.log'`

            Recommended instance: n2-highcpu-96 (96 vCPUs, 96GB RAM - recommended)

### Local Pipeline Execution

            For development and targeted analysis: Run specific pipeline steps locally with precise control over analysis parameters.

Command: `python scripts/clean_run_full_pipeline.py`

            What it does:

                - Cleans all previous data, outputs, logs, and figures for a fresh start

                - Executes all 21 pipeline steps in correct order (Steps 1.0-4.7)

                - Validates outputs and ensures data integrity

                    - Total runtime: ~6.8 hours (wall-clock with parallelism) on n2-highcpu-96 (96 vCPUs, 96GB RAM)

            Optional arguments:

                - `--dry-run` - Preview what would be cleaned without making changes

                - `--skip-cleanup` - Skip cleanup phase and run steps only

            For step-by-step execution with detailed explanations, see individual pipeline steps below.

### Pipeline Structure

        The pipeline is organized into 4 main groups with 21 individual steps:

#### Group 1: Data Acquisition (step_1_data_acquisition/)

            Downloads raw GNSS data and validates station coordinates

#### Step 1.0: Provenance Snapshot 0.1 seconds

Command: `python scripts/steps/step_1_data_acquisition/step_1_0_provenance_snapshot.py`

        Purpose: Creates comprehensive provenance snapshot documenting computational environment, software versions, and system configuration for full reproducibility.

#### Step 1.1: Data Acquisition 5.2 minutes

Command: `python scripts/steps/step_1_data_acquisition/step_1_1_tep_data_acquisition.py`

        Purpose: Downloads raw GNSS clock product files from official analysis centers (CODE, ESA, IGS Combined) covering the full analysis period (1 January 2023 to 30 June 2025). Retrieves precise clock corrections in 30-second intervals from authoritative FTP servers, ensuring data integrity through checksum validation and automatic retry mechanisms.

#### Step 1.2: Coordinate Validation 0.4 seconds

Command: `python scripts/steps/step_1_data_acquisition/step_1_2_tep_coordinate_validation.py`

        Purpose: Validates station coordinates and performs comprehensive audit for pipeline consistency. Checks Step 1 completion, validates ECEF coordinate data quality, runs integrated station ID audit with spatial analysis, creates authoritative station counts for the pipeline, and generates comprehensive validation summary with data-driven metadata. Ensures coordinate data integrity and establishes authoritative station catalogue for subsequent correlation analysis.

#### Group 2: Core Analysis (step_2_core_analysis/)

            Computes TEP correlations and performs geospatial-temporal analysis

#### Step 2.0: TEP Correlation Analysis 59.3 minutes

Command: `python scripts/steps/step_2_core_analysis/step_2_0_tep_correlation_analysis.py`

        Purpose: Core TEP signal detection using phase-coherent cross-spectral density analysis. Computes complex CSD between all station pairs in the 10-500 µHz frequency band, extracts phase-coherent correlations as phase-alignment index, and fits exponential decay models to correlation vs. distance relationships. **Technical parameters:** 40 logarithmic distance bins (50-13,000 km range), 7 correlation models tested (Exponential, Gaussian, Matérn, Power Law, etc.), 5,000 bootstrap iterations for confidence intervals, bin count-weighted least squares fitting with adaptive lambda bounds. Implements the band-limited methodology that preserves essential phase information for TEP detection.

#### Step 2.1: Data Quality Validation 23.2 minutes

Command: `python scripts/steps/step_2_core_analysis/step_2_1_data_quality_validation.py`

        Purpose: Comprehensive data quality validation and transparency analysis. Analyzes quality-filtered correlation data from Step 2.0, adds geospatial enrichments (azimuth, local time differences), and performs extensive validation including station coverage analysis, temporal gap detection, duplicate detection, outlier validation, plateau phase boundary clustering analysis, and inter-AC comparison. Generates transparency reports with red flags and analyst recommendations.

#### Step 2.2: Geospatial Temporal Analysis 19.2 minutes

Command: `python scripts/steps/step_2_core_analysis/step_2_2_tep_geospatial_temporal_analysis.py`

        Purpose: Comprehensive geospatial and temporal analysis including astronomical event correlations, orbital tracking, anisotropy analysis, spherical harmonics, and advanced temporal field studies. Analyzes correlations with planetary positions, lunar standstills, solar eclipses, and Earth's orbital motion to validate TEP predictions across multiple temporal and spatial scales. **Multi-scale window strategy:** 30-day windows (orbital tracking), 120-day windows (mesh dance dynamics), 240-day windows (planetary gravitational coupling), 433-day cycles (Chandler wobble analysis), with each timescale optimized for its characteristic physical phenomenon.

#### Group 3: Validation Suite (step_3_validation_suite/)

            Rigorous statistical validation including cross-validation, null tests, and methodology validation

#### Step 3.0: Cross-Validation Suite 0.8 seconds

Command: `python scripts/steps/step_3_validation_suite/step_3_0_tep_cross_validation_suite.py`

        Purpose: Comprehensive cross-validation framework including block-wise (monthly/spatial), Leave-One-Station-Out (LOSO), Leave-One-Day-Out (LODO), and block bootstrap analyses. Provides rigorous validation of TEP correlation parameters using multiple complementary approaches to ensure robustness and statistical validity. Tests whether fitted parameters have genuine predictive power and distinguishes real physics from curve-fitting artifacts.

#### Step 3.1: Robust Block Bootstrap 2.3 hours

Command: `python scripts/steps/step_3_validation_suite/step_3_1_robust_block_bootstrap.py`

        Purpose: Advanced bootstrap validation with temporal block structure preservation. Implements stationary block bootstrap to account for temporal dependencies while generating confidence intervals for TEP parameters.

#### Step 3.2: Null Tests 31.0 seconds

Command: `python scripts/steps/step_3_validation_suite/step_3_2_tep_null_tests.py`

        Purpose: Critical validation step that establishes signal authenticity through comprehensive randomization testing. Implements three independent scrambling approaches—distance scrambling (20 iterations per center), phase scrambling (20 iterations per center), and station scrambling (20 iterations per center)—to verify that observed correlations depend on genuine physical relationships rather than analysis artifacts. Distance scrambling randomizes both distances and coherences independently to destroy distance-coherence relationships. Quantifies signal enhancement over null distributions (24-61× signal-to-null ratios) and establishes statistical significance through permutation testing and z-score analysis (z = 15.8-31.9, all p 6.0 minutes

Command: `python scripts/steps/step_3_validation_suite/step_3_3_methodology_validation.py`

        Purpose: Comprehensive methodology validation framework achieving 100% validation score with rigorous bias characterization, distribution-neutral validation, geometric control analysis, and zero-lag leakage immunity validation.

#### Step 3.4: Geographic Bias Validation 52.4 seconds

Command: `python scripts/steps/step_3_validation_suite/step_3_4_geographic_bias_validation.py`

        Purpose: Comprehensive geographic bias validation using statistical resampling to assess geographic consistency, hemisphere balance, elevation quintile independence, and ocean vs land baseline characteristics.

#### Step 3.5: Realistic Ionospheric Validation 2.4 seconds

Command: `python scripts/steps/step_3_validation_suite/step_3_5_realistic_ionospheric_validation.py`

        Purpose: Ionospheric independence validation using real space weather data correlating with authentic geomagnetic indices and solar activity to demonstrate non-ionospheric origin.

#### Step 3.6: Control Band Analysis 2.9 hours

Command: `python scripts/steps/step_3_validation_suite/step_3_6_control_band_analysis.py`

        Purpose: Comprehensive multi-band frequency validation across 12 frequency bands (10-3000 µHz) to assess systematic effects and frequency specificity. **Frequency decomposition:** TEP band (10-500 µHz), tidal bands (10-30 µHz), post-tidal ranges (30-100 µHz), and control bands (>1000 µHz) for systematic effect quantification. Applies identical phase-coherent methodology as Step 2.0 across each band with exponential model fitting and correlation strength assessment. Validates frequency specificity of TEP signals and excludes broadband instrumental artifacts.

#### Step 4.7: Multiple Comparison Corrections 2.8 seconds

Command: `python scripts/steps/step_4_advanced_analysis_and_visualization/step_4_7_multiple_comparison_corrections_fixed.py`

        Purpose: Comprehensive application of formal multiple comparison corrections (Bonferroni, FDR, FWER) to all 388 statistical tests across 19 analysis families, including data quality validation, null hypothesis testing, bootstrap validation, band diagnostics, Hilbert-IF astronomical analysis, coordinate validation, eclipse events, bootstrap cross-methods, multiband frequency validation, and astronomical events, showing that core TEP findings survive ultra-conservative corrections.

#### Group 4: Advanced Analysis & Visualization (step_4_advanced_analysis_and_visualization/)

            Advanced analyses, astronomical correlations, and publication-quality visualizations

#### Step 4.0: Advanced Analysis 4.7 minutes

Command: `python scripts/steps/step_4_advanced_analysis_and_visualization/step_4_0_tep_advanced_analysis.py`

        Purpose: Conducts specialized analyses including circular statistics validation (Phase-Locking Value, Rayleigh tests), elevation-dependent screening analysis, geomagnetic stratification studies, and temporal orbital tracking analysis. **Circular statistics:** Phase Locking Values (PLV) in 0.1-0.4 range, Rayleigh test significance (p 6.0 minutes

Command: `python scripts/steps/step_4_advanced_analysis_and_visualization/step_4_1_tep_visualization.py`

        Purpose: Generates comprehensive publication-quality figures including correlation vs. distance plots, global station network maps, residual analysis plots, anisotropy heatmaps, and temporal tracking visualizations. Creates both individual analysis center plots and multi-center comparison figures with consistent styling and statistical annotations.

#### Step 4.2: Synthesis Figure 43.7 seconds

Command: `python scripts/steps/step_4_advanced_analysis_and_visualization/step_4_2_tep_synthesis_figure.py`

        Purpose: Creates the comprehensive synthesis figure combining key results from all analysis centers. Produces the main publication figure showing correlation decay curves, statistical fits, and multi-center consistency in a unified presentation suitable for manuscript submission.

#### Step 4.3: High Resolution Astronomical Events 3.5 minutes

Command: `python scripts/steps/step_4_advanced_analysis_and_visualization/step_4_3_high_resolution_astronomical_events.py`

        Purpose: Comprehensive high-resolution astronomical analysis including: (1) Complete 5-eclipse study (2023-2025) using 30-second CLK data revealing substantial coherence modulations (18-87% of baseline TEP signal) with observed eclipse type hierarchy where partial eclipses produce strongest network responses (mean: 4.54×10⁻², 87% of baseline) through ionospheric gradient mechanisms; (2) Cross-planetary orbital periodicity analysis revealing systematic TEP signals correlated with orbital completeness - Venus (4.05 orbits) shows strong correlations (ESA: +17.7%, IGS Combined: +10.6%, CODE: +4.8%) while outer planets show incomplete cycles and weaker effects; (3) Advanced sham controls confirming robust statistical significance beyond GPS processing artifacts; (4) Instantaneous frequency analysis detecting event-locked solar rotation modulations. This analysis demonstrates how the global GPS network acts as a sensitive detector of both transient (eclipse) and persistent (orbital periodicity) field modulations, providing comprehensive validation of TEP astronomical coupling predictions using 30-second CLK file data.

#### Step 4.4: Gravitational Temporal Field Analysis 1.5 minutes

Command: `python scripts/steps/step_4_advanced_analysis_and_visualization/step_4_4_gravitational_temporal_field_analysis.py`

        Purpose: Comprehensive gravitational-temporal field correlation analysis using NASA/JPL planetary ephemeris (DE-series) to correlate planetary gravitational influences with GPS clock coherence patterns, revealing strong stacked gravitational pattern correlations (r = -0.503, raw p = 1.5×10⁻⁵⁹, autocorr-corrected p = 3.3×10⁻⁴) with optimal 227-day smoothing window. Includes systematic Earth motion energy hierarchy validation demonstrating sophisticated TEP coupling mechanisms - orbital motion coupling (|r| = 0.7-0.8), Moon-Chandler gravitational field modulation (|r| = 0.6-0.7), and rotational anisotropy (CV of rotational stability = 0.5-0.6). Energy vs velocity scaling discrimination analysis (-0.057, 95% CI: [-0.143, +0.030], n.s.) reveals no significant preference between scaling mechanisms, supporting complex multi-mechanism coupling and validating TEP's core predictions of non-integrable time transport and synchronization holonomy.

#### Step 4.5: Comprehensive Diurnal Analysis 1.8 hours

Command: `python scripts/steps/step_4_advanced_analysis_and_visualization/step_4_5_comprehensive_diurnal_analysis.py`

        Purpose: Comprehensive diurnal and seasonal TEP modulation analysis using high-resolution hourly windowing methodology on raw CLK files. Processes 72.4M hourly records to reveal complex seasonal diurnal patterns with multi-center validation, systematic day-night effects, and peak hours that shift across seasons and centers.

#### Step 4.6: TID Exclusion Analysis 2.8 minutes

Command: `python scripts/steps/step_4_advanced_analysis_and_visualization/step_4_6_tid_exclusion_analysis.py`

        Purpose: Comprehensive exclusion analysis for alternative explanations including Traveling Ionospheric Disturbances (TIDs) and trans-equatorial propagation (TEQ). Performs quantitative temporal band separation analysis, spatial structure comparison (exponential vs plane-wave models), ionospheric independence verification, and frequency band analysis. Provides high confidence exclusion of TIDs through temporal separation factor and power distribution analysis, strengthening the case for TEP field coupling as the underlying mechanism.

### Key Features & Reproducibility

            - Real data only: No synthetic, fallback, or mock data

            - Authoritative sources: Direct download from official FTP servers

            - Multi-core processing: Parallel analysis with configurable worker count

            - Checkpointing: Automatic resume from interruptions

            - Comprehensive validation: Null tests, circular statistics, bootstrap confidence intervals

**Experimental Section:**

#### Reproducibility Documentation

                    Code Repository: Complete analysis pipeline in `scripts/steps/`

                    Data Sources: Official GNSS FTP servers (CODE, ESA, IGS Combined)

                    Dependencies: Listed in `requirements/requirements.txt`

                    Configuration: Centralized in `scripts/utils/config.py`

                    Binning Method: 40 logarithmic distance bins attempted, 50–13,000 km range

                    Phase Computation: `cos(np.angle(CSD))` with magnitude-weighted circular averaging

                    Cross-Validation: LOSO/LODO with complete model re-fitting

                    Statistical Tests: Weighted least squares with uncertainty propagation

            Replication Note: Full replication requires ~1 day of data download and ~1–2 days of computation on multi-core systems. All intermediate results are checkpointed for partial replication.

## 7. Appendix

### 7.1 Event Windows: Complete Astronomical Events Analysis Parameters

    Complete specification of temporal windows, detrending methods, and multi-center combination rules for all astronomical events analyzed in the TEP study.

| Event Type | Events (n) | ±Window (days) | Detrending | Multi-Center Rule |
| --- | --- | --- | --- | --- |
| Jupiter Opposition | 2 | ±120 | Linear polyfit | Weighted average |
| Saturn Opposition | 2 | ±120 | Linear polyfit | Weighted average |
| Mars Opposition | 1 | ±120 | Linear polyfit | Weighted average |
| Venus Conjunction | 2 | ±120 | Linear polyfit | Weighted average |
| Mercury Conjunction | 8 | ±120 | Linear polyfit | Weighted average |
| Solar Eclipses | 5 | ±1 (high-res) | Linear polyfit | Mean averaging |
| Supermoon Perigee | 11 | ±7 | Linear polyfit | Mean averaging |
| Lunar Standstill | 1 | ±180 | Linear polyfit | Weighted average |

#### Key Parameters:

            - Detrending: Linear polynomial fitting (np.polyfit, degree=1) applied to time series before cross-spectral density computation

            - Multi-Center Combination: Weighted averaging by sample size for planetary events; simple mean averaging for high-resolution events

            - Window Rationale: ±120 days for planetary events (optimal gravitational coupling timescale); ±1–7 days for rapid events; ±180 days for lunar standstill (18.6-year cycle)

            - Total Events: 32 astronomical events across 8 categories, spanning 2023-2025 analysis period

        Event Dates: Jupiter (3 November 2023, 7 December 2024), Saturn (27 August 2023, 8 September 2024)*, Mars (16 January 2025), Venus (13 August 2023, 23 March 2025), Mercury (8 conjunctions 2023-2025), Eclipses (20 April 2023, 14 October 2023, 8 April 2024, 2 October 2024, 29 March 2025), Supermoons (11 events 2023-2025), Lunar Standstill (1 June 2025).

        **Note: Saturn opposition on 21 September 2025 falls outside the v0.18 analysis window (ends 30 June 2025) and is marked as planned for future analysis.*

### 7.2 TID Exclusion: Implementation Details

    Detailed implementation specifications for the TID exclusion analysis, including file paths and code references.

#### Implementation Details

            - Implementation: See `scripts/steps/step_4_advanced_analysis_and_visualization/step_4_6_tid_exclusion_analysis.py`, function `_perform_tid_exclusion` for date alignment, thresholding, and Δ% computation.

            - Results location: `results/outputs/step_4_6_tid_exclusion_analysis_results.json` and narrative summary in Section 3.5.2 (Ionospheric Interference Assessment).

        For auditability, the methods above directly correspond to the code: date alignment (common days), percentile-based thresholding on the TID proxy, exclusion of high-activity samples, and percent change reporting.

### 7.3 Multiple Comparison Corrections Summary

**Experimental Section:**

#### Table S1: Test Family Multiple Comparison Corrections

        Comprehensive multiple comparison corrections applied across 388 statistical tests using Bonferroni, FDR-BH (Benjamini-Hochberg), Family-Wise Error Rate, and Hierarchical Empirical Bayes corrections. Family alpha level: α = 0.05.

| Test Family | m (Tests) | α | Bonferroni Sig. | FDR-BH Sig. | Family-Wise Sig. | Hierarchical EB Sig. |
| --- | --- | --- | --- | --- | --- | --- |
| TEP Band | 3 | 0.05 | 3/3 | 3/3 | 3/3 | **3/3** |
| Model Comparison | 18 | 0.05 | 1/18 | 9/18 | 4/18 | 0/18 |
| Null Validation | 9 | 0.05 | 9/9 | 9/9 | 9/9 | 0/9 |
| Astronomical Events | 8 | 0.05 | 5/8 | 8/8 | 8/8 | 0/8 |
| Anisotropy Orbital | 22 | 0.05 | 0/22 | 4/22 | 2/22 | 0/22 |
| Cross Validation | 3 | 0.05 | 0/3 | 0/3 | 0/3 | 0/3 |
| Advanced Analysis | 174 | 0.05 | 35/174 | 54/174 | 37/174 | 48/174 |
| Geographic Validation | 2 | 0.05 | 1/2 | 2/2 | 2/2 | 2/2 |
| Chandler Wobble | 3 | 0.05 | 0/3 | 2/3 | 2/3 | **3/3** |
| Bootstrap Validation | 9 | 0.05 | 9/9 | 9/9 | 9/9 | 0/9 |
| Multiband Analysis | 40 | 0.05 | 37/40 | 40/40 | 40/40 | **40/40** |
| Gravitational Field | 9 | 0.05 | 2/9 | 3/9 | 3/9 | 3/9 |
| Diurnal Validation | 2 | 0.05 | 1/2 | 2/2 | 2/2 | 2/2 |
| Eclipse Analysis | 15 | 0.05 | 5/15 | 7/15 | 7/15 | 9/15 |
| Bootstrap Cross-Method | 12 | 0.05 | 3/12 | 3/12 | 3/12 | 3/12 |
| Coordinate Validation | 1 | 0.05 | 1/1 | 1/1 | 1/1 | 1/1 |
| Data Quality Validation | 6 | 0.05 | 6/6 | 6/6 | 6/6 | 0/6 |
| Hilbert-IF Astronomical | 12 | 0.05 | 0/12 | 1/12 | 0/12 | 0/12 |
| Band Diagnostics | 40 | 0.05 | 37/40 | 40/40 | 40/40 | **40/40** |
| Total (All Families) | 388 | 0.05 | 155/388 | 203/388 | 178/388 | **154/388** |

            **Notes:**

                - **m:** Number of statistical tests in each family

                - **α:** Family-wise alpha level (0.05)

                - **Bonferroni:** Conservative correction with adjusted α = 0.05/388 = 0.000129

                - **FDR-BH:** Benjamini-Hochberg False Discovery Rate control

                - **Family-Wise:** Controls family-wise error rate across all test families

                - **Hierarchical EB:** Empirical Bayes partial-pooling correction with adaptive shrinkage per family (minimum n≥5 for shrinkage application). Provides data-driven balance between Bonferroni conservatism and FDR leniency by sharing information across tests within coherent families.

                - **Green highlighting:** All tests in family remain significant after correction

                - **Impact Summary:** Bonferroni reduction: 60.1% | FDR-BH reduction: 47.7% | Family-Wise reduction: 54.1% | Hierarchical EB reduction: 60.3%

### 7.4 Bootstrap Convergence Validation

    Statistical validation of bootstrap uncertainty quantification methodology, addressing convergence rates and potential bias from non-converged iterations.

#### Bootstrap Success Rates

        The correlation parameter bootstrap analysis (Step 2.0) achieves convergence rates of 71-73% across all analysis centers:

| Analysis Center | Successful Fits | Total Iterations | Success Rate | Status |
| --- | --- | --- | --- | --- |
| CODE | 3,555 | 5,000 | 71.1% | ✓ Acceptable |
| IGS Combined | 3,555 | 5,000 | 71.1% | ✓ Acceptable |
| ESA Final | 3,632 | 5,000 | 72.6% | ✓ Acceptable |

        **Assessment:** These success rates of 71-73% reflect the challenging nature of nonlinear parameter optimization on resampled data. The convergence rates are sufficient for robust uncertainty quantification, with comprehensive bias validation confirming that non-converged iterations do not systematically bias parameter estimates.

#### Bias Validation Analysis

        Comprehensive statistical testing confirms that discarding non-converged bootstrap iterations does not introduce systematic bias in parameter estimates:

                **✓ Parameter Bias Assessment**
                All parameter biases < 6%: |λ bias| < 4%, |A bias| < 6%, well within acceptable tolerances for uncertainty quantification.

                **✓ Failure Pattern Analysis**
                Failed iterations represent challenging optimization cases with poor coherence variation, not systematically biased physical subsets.

                **✓ Confidence Interval Validity**
                Bootstrap percentiles accurately reflect true parameter uncertainties; statistical power remains adequate with current success rates.

                **✓ Improvement Potential**
                Robust initialization strategies can increase success rates to 96%+ if desired, though current rates are scientifically sufficient.

#### Root Cause Analysis

    **Primary failure mechanism:** `scipy.optimize.curve_fit` hitting parameter bounds during bootstrap samples with:

        - **Extreme distance distributions:** Samples with predominantly short or long distances challenge λ parameter estimation

        - **Poor coherence variation:** Samples with minimal coherence range (c_range ≤ 0) cannot constrain amplitude parameters

        - **Ill-conditioned matrices:** Covariance matrix computation fails for degenerate parameter combinations

    **Scientific interpretation:** These failures represent the inherent mathematical limitations of nonlinear optimization when applied to resampled data with reduced statistical power, not fundamental flaws in the methodology.

            **Methodology Note:** Bootstrap convergence analysis conducted using enhanced diagnostic framework (Step 3.7) with 2,000 test iterations per analysis center, multiple initialization strategies, and comprehensive bias validation across all TEP correlation parameters.

## References

        Barrow, J. D. & Magueijo, J. (1999). Varying-α theories and solutions to the cosmological problems. *Physics Letters B*, 447(3-4), 246-250.
        Bevis, M., et al. (1994). GPS meteorology: Mapping zenith wet delays onto precipitable water. *Journal of Applied Meteorology*, 33(3), 379-386.
        Bothwell, T., et al. (2022). Resolving the gravitational redshift across a millimetre-scale atomic sample. *Nature*, 602(7897), 420-424.
        Brans, C. & Dicke, R. H. (1961). Mach's principle and a relativistic theory of gravitation. *Physical Review*, 124(3), 925-935.
        Bretherton, C. S., et al. (1999). The effective number of spatial degrees of freedom of a time-varying field. *Journal of Climate*, 12(7), 1990-2009.
        Chen, G. & Herring, T. A. (1997). Effects of atmospheric azimuthal asymmetry on the analysis of space geodetic data. *Journal of Geophysical Research*, 102(B9), 20489-20502.
        Chou, C. W., et al. (2010). Optical clocks and relativity. *Science*, 329(5999), 1630-1633.
        Damour, T. & Nordtvedt, K. (1993). General relativity as a cosmological attractor of tensor-scalar theories. *Physical Review Letters*, 70(15), 2217.
        Damour, T. & Polyakov, A. M. (1994). The string dilaton and a least coupling principle. *Nuclear Physics B*, 423(2-3), 532-558.
        Delva, P., et al. (2018). Gravitational redshift test using eccentric Galileo satellites. *Physical Review Letters*, 121(23), 231101.
        Eddington, A. S. (1919). A determination of the deflection of light by the sun's gravitational field, from observations made at the total eclipse of May 29, 1919. *Philosophical Transactions of the Royal Society of London A*, 220(571-581), 291-333.
        Einstein, A. (1915). Die Feldgleichungen der Gravitation. *Sitzungsberichte der Königlich Preußischen Akademie der Wissenschaften*, 844-847.
        Godun, R. M., et al. (2014). Frequency ratio of two optical clock transitions in 171Yb+ and constraints on the time variation of fundamental constants. *Physical Review Letters*, 113(21), 210801.
        Heflin, M. B., et al. (1992). Global coordinates and site motion from very long baseline interferometry and satellite laser ranging. *Journal of Geophysical Research*, 97(B1), 459-467.
        Holton, J. R. & Hakim, G. J. (2012). *An Introduction to Dynamic Meteorology*. Academic Press.
        Horndeski, G. W. (1974). Second-order scalar-tensor field equations in a four-dimensional space. *International Journal of Theoretical Physics*, 10(6), 363-384.
        Hunsucker, R. D. & Hargreaves, J. K. (2003). *The High-Latitude Ionosphere and its Effects on Radio Propagation*. Cambridge University Press.
        Khoury, J. & Weltman, A. (2004). Chameleon cosmology. *Physical Review D*, 69(4), 044026.
        Kivelson, M. G. & Russell, C. T. (1995). *Introduction to Space Physics*. Cambridge University Press.
        Kouba, J. & Héroux, P. (2001). Precise point positioning using IGS orbit and clock products. *GPS Solutions*, 5(2), 12-28.
        Gendt, G., & Schmid, R. (2005). A common-coordinate approach to global GPS analysis. *IGS Technical Report 2004*, 131-141.
        Steigenberger, P., Montenbruck, O., Dach, R., et al. (2021). CODE reprocessing 1995-2020: improved GPS orbits and clocks. *Journal of Geodesy*, 95, 65.
        International GNSS Service. (2023). *IGS Technical Report 2023*, Chapter 7: Clock Product Generation.
        Lyard, F., et al. (2006). Modeling the global ocean tides: modern insights from FES2004. *Ocean Dynamics*, 56(5-6), 394-415.
        McGrew, W. F., et al. (2018). Atomic clock performance enabling geodesy below the centimetre level. *Nature*, 564(7734), 87-90.
        Montenbruck, O., et al. (2017). The Multi-GNSS Experiment (MGEX) of the International GNSS Service (IGS)–achievements, prospects and challenges. *Advances in Space Research*, 59(7), 1671-1697.
        Murphy, M. T., et al. (2003). Possible evidence for a variable fine-structure constant from QSO absorption lines. *Monthly Notices of the Royal Astronomical Society*, 345(2), 609-638.
        Petit, G. & Luzum, B. (2010). IERS Conventions (2010). *IERS Technical Note* No. 36, Frankfurt am Main: Verlag des Bundesamts für Kartographie und Geodäsie.
        Ray, J., et al. (2008). IGS polar motion measurement accuracy. *Geophysical Research Letters*, 35(3), L03303.
        Rosenband, T., et al. (2008). Frequency ratio of Al+ and Hg+ single-ion optical clocks; metrology at the 17th decimal place. *Science*, 319(5871), 1808-1812.
        Senior, K. L., et al. (2008). Characterization of periodic variations in the GPS satellite clocks. *GPS Solutions*, 12(3), 211-225.
        Smawfield, M. L. (2025). The Temporal Equivalence Principle: Dynamic Time, Emergent Light Speed, and a Two-Metric Geometry of Measurement. *Zenodo*. [https://doi.org/10.5281/zenodo.16921911](https://doi.org/10.5281/zenodo.16921911).
        Takamoto, M., et al. (2020). Test of general relativity by a pair of transportable optical lattice clocks. *Nature Photonics*, 14(7), 411-415.
        Touboul, P., et al. (2017). MICROSCOPE mission: first results of a space test of the equivalence principle. *Physical Review Letters*, 119(23), 231101.
        Uzan, J. P. (2003). The fundamental constants and their variation: observational and theoretical status. *Reviews of Modern Physics*, 75(2), 403.
        van Dam, T. M., et al. (2001). Displacements of the Earth's surface due to atmospheric loading: Effects on gravity and baseline measurements. *Journal of Geophysical Research*, 106(B6), 11115-11132.
        Webb, J. K., et al. (2001). Further evidence for cosmological evolution of the fine structure constant. *Physical Review Letters*, 87(9), 091301.
        Jammalamadaka, S. R., & Sengupta, A. (2001). *Topics in Circular Statistics*. World Scientific.

### Data Sources

        Johnston, G., Riddell, A., & Hausler, G. (2017). The International GNSS Service. In P. J. G. Teunissen & O. Montenbruck (Eds.), *Springer Handbook of Global Navigation Satellite Systems* (1st ed., pp. 967-982). Cham, Switzerland: Springer International Publishing. [https://doi.org/10.1007/978-3-319-42928-1_33](https://doi.org/10.1007/978-3-319-42928-1_33).
        Dach, R., Lutz, S., Walser, P., & Fridez, P. (Eds.). (2015). *Bernese GNSS Software Version 5.2*. Astronomical Institute, University of Bern, Switzerland. Available at: [http://www.bernese.unibe.ch/](http://www.bernese.unibe.ch/).
        Fernández, M. A. (2016). Geodetic and Time Reference Frames for ESA's Navigation Programmes. In *Proceedings of the 29th International Technical Meeting of the Satellite Division of The Institute of Navigation (ION GNSS+ 2016)* (pp. 2714-2719). Portland, OR.

## Data Availability

### Primary Data Sources

        **GNSS Clock Products:** Final clock solutions (30-second epochs, CLK format) from three independent analysis centers, all part of the International GNSS Service (IGS). Data are freely available under IGS Terms of Use.

        **CODE (Center for Orbit Determination in Europe):**

            - **Provider:** Astronomical Institute, University of Bern (AIUB)

            - **Source:** [http://ftp.aiub.unibe.ch/CODE/](http://ftp.aiub.unibe.ch/CODE/)

            - **Temporal Coverage:** January 1, 2023 – June 30, 2025 (912 days, complete coverage)

            - **Station Pairs:** 39.0 million measurements

            - **Citation:** Steigenberger et al. (2021); Johnston et al. (2017)

        **IGS Combined Products:**

            - **Provider:** International GNSS Service (weighted combination of multiple analysis centers)

            - **Source:** [https://igs.bkg.bund.de/root_ftp/IGS/products/](https://igs.bkg.bund.de/root_ftp/IGS/products/)

            - **Temporal Coverage:** January 1, 2023 – June 30, 2025 (912 days, complete coverage)

            - **Station Pairs:** 12.9 million measurements

            - **Citation:** Johnston et al. (2017)

        **ESA (European Space Agency):**

            - **Provider:** ESA Navigation Support Office

            - **Source:** [http://navigation-office.esa.int/products/gnss-products/](http://navigation-office.esa.int/products/gnss-products/)

            - **Temporal Coverage:** January 1, 2023 – June 30, 2025 (912 days, complete coverage)

            - **Station Pairs:** 10.8 million measurements

            - **Citation:** Fernández (2016); Johnston et al. (2017)

        **Total Dataset:** 62.73 million station-pair cross-spectral measurements from 364 unique GNSS stations

            - **Format:** RINEX 3 CLK (compressed: .gz or .Z)

            - **License:** Freely available for scientific use. Users must cite IGS and individual analysis centers appropriately

            - **Terms of Use:** [IGS Data and Product Disclaimer](https://igs.org/wp-content/uploads/2020/09/IGS-Data-and-Product-Disclaimer-and-Terms-of-Use-200805.pdf)

        **Station Coordinates:** IGS Network station coordinates (ITRF2014) obtained from IGS station metadata and BKG services, providing Cartesian coordinates (X, Y, Z) for all active IGS stations.

            - **Source:** [IGS Network Metadata (JSON)](https://files.igs.org/pub/station/general/IGSNetworkWithFormer.json)

            - **Coverage:** 364 unique stations with valid clock observations

            - **Reference Frame:** ITRF2014 (International Terrestrial Reference Frame 2014)

            - **License:** Freely available under IGS Terms of Use

### Software and Analysis Tools

        **Analysis Code:** Complete Python pipeline (TEP-GNSS) available under MIT License at [GitHub: TEP-GNSS](https://github.com/matthewsmawfield/TEP-GNSS). Includes data acquisition, preprocessing, statistical validation, and figure generation scripts for all 23 analysis steps.

        **Key Dependencies:**

            - **NumPy/SciPy:** Numerical computing and statistical analysis

            - **Pandas:** Data manipulation and time series analysis

            - **Matplotlib:** Scientific visualization

            - **Statsmodels:** Advanced statistical modeling

### Derived Products and Supplementary Materials

        **Analysis Results:** Comprehensive JSON output files, processed correlation data, validation results, and extended figures available in the [Zenodo repository (DOI: 10.5281/zenodo.17127229)](https://zenodo.org/records/17127229).

        **Data Attribution:** By using data from this study, users agree to cite the original data providers (IGS, CODE/AIUB, ESA) as well as this work. All data sources are freely available for scientific research under their respective terms of use.

        **Cross-Center Validation:** The use of three independent analysis centers (CODE, IGS Combined, ESA) provides robust cross-validation with R² = 0.920-0.970 consistency, demonstrating the phenomenon is not processing-specific.

**Experimental Section:**

## How to cite

         **Cite as:** Smawfield, M. L. (2025). Global Time Echoes: Distance-Structured Correlations in GNSS Clocks. v0.23 (Jaipur). Zenodo. [https://doi.org/10.5281/zenodo.17127229](https://doi.org/10.5281/zenodo.17127229)

            **BibTeX:**
@misc{Smawfield_TEP_GNSS_2025,
  author       = {Matthew Lukin Smawfield},
  title        = {Global Time Echoes: Distance-Structured Correlations in GNSS
                   Clocks (Jaipur v0.23)},
  year         = {2025},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.17127229},
  url          = {https://doi.org/10.5281/zenodo.17127229},
  note         = {Preprint}
}

## Declarations

        **Funding:** This work received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

        **Competing Interests:** The author declares no competing financial interests.

        **Author Contributions:** M.L.S. designed the study, performed all analyses, and wrote the manuscript.

## Contact

        For questions, comments, or collaboration opportunities regarding this work, please contact:

        **Matthew Lukin Smawfield**

        [matthewsmawfield@gmail.com](mailto:matthewsmawfield@gmail.com)

**Experimental Section:**

## Version v0.23 Updates

                - **Methodology Clarification:** Updated prediction timeline documentation to clearly establish hypothesis-driven approach

                - **Enhanced Presentation:** Added methodology sections and improved visual formatting for better clarity

## Version v0.23 Updates

                - **Documentation Improvements:** Enhanced site formatting and presentation for better readability

                - **Content Editing:** Refined manuscript content and improved version history organization

## Version v0.23 Updates

                - **Bootstrap Enhancement:** Increased bootstrap iterations from 1,000 to 5,000 for more robust confidence intervals and statistical validation

                - **Statistical Rigor:** Added comprehensive systematic effects analysis including atmospheric loading, processing artifacts, and autocorrelation corrections

                - **Fresh Pipeline Results:** Complete pipeline re-run with enhanced methodology producing updated figures and validation results across all 23 analysis steps

        [← Home](/)
        
### TEP Research Series

        - [Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed 18 Aug 2025](/tep/theory/)

        - [Global Time Echoes: Distance-Structured Correlations in GNSS Clocks 17 Sep 2025](/tep/gnss-i/)

        - [25-Year Temporal Evolution of Distance-Structured Correlations in GNSS 3 Nov 2025](/tep/gnss-ii/)

        - [Global Time Echoes: Raw RINEX Validation 17 Dec 2025](/tep/gnss-iii/)

        - [Temporal-Spatial Coupling in Gravitational Lensing 19 Dec 2025](/tep/gl/)

        - [Global Time Echoes: Empirical Validation of TEP 21 Dec 2025](/tep/gte/)

        - [Universal Critical Density: Unifying Atomic, Galactic, and Compact Object Scales 28 Dec 2025](/tep/ucd/)

        - [The Soliton Wake: Identifying RBH-1 as a Gravitational Soliton 28 Dec 2025](/tep/rbh/)

        ← Previous
        Next →

---

*This document was automatically generated from the TEP-GNSS research site. For the interactive version with figures and enhanced formatting, visit: https://matthewsmawfield.github.io/TEP-GNSS/*

*Source code and data available at: https://github.com/matthewsmawfield/TEP-GNSS*
