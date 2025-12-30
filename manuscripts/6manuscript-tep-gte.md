# Global Time Echoes: Empirical Validation

**Author:** Matthew Lukin Smawfield  
**Version:** v0.2 (Singapore)  
**Date:** First published: 21 December 2025 · Last updated: 25 December 2025  
**DOI:** 10.5281/zenodo.18004832  
**Generated:** 2025-12-29  
**Paper Series:** TEP-GTE

---

## Abstract

        Analysis of 25.3 years of GNSS timing data (2000–2025) reveals a persistent, distance-structured correlation in global atomic clock networks that tests an empirically untested assumption of general relativity: the global integrability of proper time. Examination of 165.2 million station pairs from 474 unique receivers demonstrates a spatial correlation signal decaying exponentially with distance (λ = 4,201 ± 1,967 km, R² = 0.92–0.97 across three independent analysis centers). These findings emerge from a systematic five-paper research program: theoretical framework development with quantitative predictions (Paper 0), multi-center validation across independent processing pipelines (Paper 1), 25-year longitudinal analysis enabling long-period geophysical detection (Paper 2), raw data confirmation eliminating processing artifacts (Paper 3), and cosmological extension connecting terrestrial correlations to dark matter phenomenology through gravitational lensing (Paper 4).

        Seven statistically independent signatures emerge with joint probability p ≈ 2×10−27 (>10σ): exponential spatial decay; East-West/North-South anisotropy (ratio 2.16, p < 10−15); orbital velocity coupling (r = −0.888, 5.1σ); alignment with the Cosmic Microwave Background dipole (18.2° separation, 5,570× variance ratio over galactic motion); planetary event responses (56/156 significant at ≥2σ); 18.6-year lunar nutation coupling (R² = 0.641); and semiannual nutation coupling (R² = 0.904). Raw RINEX validation using Single Point Positioning with broadcast ephemerides achieves 100% detection rate across 72 metric combinations (t-statistics up to 112, Cohen's d up to 0.304), excluding processing artifacts as the origin. The network's selectivity profile—sensitive to velocity-dependent dynamics while blind to GM/r² scaling and solar rotation—characterizes it as an inertial interferometer measuring correlation geometry rather than a gravimeter measuring Newtonian force.

        These observations match a priori predictions of the Temporal Equivalence Principle, a bi-metric scalar-tensor framework in which proper time is a dynamical field governed by a conformal factor A(φ) = exp(2βφ/MPl). The observed correlation length corresponds to a scalar field mass mφ ≈ (4.34–5.93)×10−14 eV/c², consistent with Vainshtein screening at the dark energy scale Λ ~ 10−13 eV. The framework preserves local Lorentz invariance while predicting global path-dependent synchronization through spatial correlations in the φ field. Critically, the conformal sector responsible for clock-rate modulation remains unconstrained by GW170817, which bounds only disformal (cone-tilting) effects.

        If validated through independent replication, TEP implies that dark matter phenomenology in gravitational lensing arises from temporal-field gradients rather than particulate matter—the projection of differential proper-time accumulation onto observations that assume the Isochrony Axiom. The 4,000 km correlation on Earth and the 50 kpc dark matter halo in galaxies represent the same scalar field at different density scales, connected by the universal M1/3 Vainshtein scaling law. Explicit falsification criteria include: failure of independent groups to replicate the raw carrier-phase signal; correlation length falling outside the 500–20,000 km range; confirmation that the signal arises from ephemeris artifacts rather than physical clock correlations (via Satellite Laser Ranging validation); and null synchronization holonomy in closed-loop triangular time-transfer experiments.

        *Keywords:* Temporal Equivalence Principle, GNSS, atomic clocks, CMB alignment, dark matter, gravitational lensing, scalar-tensor gravity, synchronization holonomy, Vainshtein screening

## 1. Introduction

## 1.1 The Synchronization Residual Problem

    Modern atomic clocks have achieved extraordinary precision, with optical lattice clocks demonstrating fractional frequency stability at the 10−18 level—sufficient to detect gravitational redshifts from centimeter-scale height differences (Bothwell et al. 2022). Yet despite this remarkable local precision, global clock synchronization residuals plateau at 10−15, three orders of magnitude above the fundamental limit. This discrepancy is conventionally attributed to systematic noise arising from ionospheric fluctuations, tropospheric delays, multipath interference, and instrumental drift.

    The Global Navigation Satellite System represents humanity's most extensive precision timing network, comprising hundreds of ground stations that maintain nanosecond-level synchronization across thousands of kilometers (Hofmann-Wellenhof et al. 2008; Teunissen & Montenbruck 2017). Standard processing pipelines apply sophisticated corrections for all known relativistic effects—gravitational redshift, velocity-dependent time dilation, Sagnac rotation, and Shapiro delay—following IERS Conventions to millimeter precision (Petit & Luzum 2010). After these corrections, residual correlations are expected to be consistent with random noise.

    Analysis of 25.3 years of GNSS timing data reveals a different picture. A systematic, distance-structured correlation persists in the residuals, characterized by a decay length of approximately 4,000 km. Standard GNSS processing explicitly suppresses spatial correlations to optimize positioning accuracy. The signal reported here survives precisely because it is geometric and differential rather than energetic or common-mode, placing it outside the correction basis of conventional pipelines. This correlation has been present in GNSS data since the network's inception but was categorized as systematic error rather than physical phenomenon. The present paper argues that this unexplained correlation structure represents a physical signal—specifically, the signature of a dynamical time field predicted by scalar-tensor extensions of general relativity.

## 1.2 Historical Context

    If this signal is physical, why was it not identified in the first decades of GNSS operation? Three developments converged to enable its detection.

    The first is temporal baseline. The 25.3-year dataset now available spans 1.4 complete cycles of the 18.6-year lunar nutation, enabling detection of long-period geophysical coupling that was inaccessible to earlier studies. Analyses conducted prior to 2010 lacked sufficient temporal leverage to distinguish such signals from transient artifacts.

    The second is multi-constellation cross-validation. The maturation of GLONASS, Galileo, and BeiDou provides independent verification across different satellite geometries and processing chains. Single-constellation analyses are vulnerable to constellation-specific artifacts that multi-constellation cross-validation can identify and exclude.

    The third is methodological. Standard GNSS processing employs Kalman filtering and network adjustment algorithms designed to eliminate correlations, thereby optimizing positioning accuracy. The present analysis instead preserves correlation structure through phase-coherent cross-spectral methods applied to residuals. The signal was present throughout the history of GNSS operation but was categorized as systematic error rather than physical phenomenon.

## 1.3 Theoretical Motivation

    The standard relativistic framework treats proper time as a derived quantity: given the metric tensor gμν, proper time along any worldline is uniquely determined by integration. This formulation assumes that time transport is integrable—that synchronization around any closed loop returns to the initial epoch after accounting for known effects such as Sagnac rotation, Shapiro delay, and gravitational redshift.

    Scalar-tensor theories of gravity, developed extensively over the past three decades (Bekenstein 1993; Damour & Polyakov 1994; Khoury & Weltman 2004), predict modifications of proper time accumulation to vary spatially. If this assumption proves empirically false in the presence of conformal metric couplings, apparent dark matter could arise as a geometric artifact: temporal depth projected onto the spatial plane. If proper time is not merely derived from geometry but is itself influenced by a dynamical field, synchronization becomes path-dependent in ways that standard general relativity does not predict.

    The Temporal Equivalence Principle formalizes this possibility. It states that all non-gravitational dynamics, signals, and quantum phases evolve according to the proper time defined by a matter metric g̃μν that couples universally to matter through the scalar field. In local freely falling frames, physics reduces to special relativity with invariant c; globally, the flow of time becomes path-dependent in ways that distributed timing networks can detect.

### 1.3.1 Research Program Structure

    This integrated manuscript synthesizes results from a systematic five-paper research program designed to establish, validate, and independently confirm the TEP predictions through progressively rigorous empirical tests and cosmological extension. The program follows a theory-first methodology: quantitative predictions established before empirical analysis, followed by multiple independent validations using complementary methodologies, culminating in cosmological implications.

    Paper 0: Theoretical Foundation (Smawfield 2025a, DOI: 10.5281/zenodo.16921911) — The Temporal Equivalence Principle framework was developed independently of any GNSS analysis, establishing quantitative predictions including a correlation length λ = 1,000–10,000 km derived from environmental screening theory, exponential spatial decay form, velocity-dependent anisotropy, and absence of GM/r² scaling. This theory-first approach ensures that subsequent empirical findings represent genuine predictions rather than post-hoc explanations.

    Paper 1: Multi-Center Validation (Smawfield 2025b, DOI: 10.5281/zenodo.17127229) — Initial empirical investigation analyzed 62.7 million station-pair measurements from 364 GNSS stations across three independent analysis centers (CODE, IGS, ESA) over 2.5 years (2023–2025). Cross-center consistency (R² = 0.92–0.97) established that the correlation structure is not processing-specific, with λ = 3,330–4,549 km falling within the predicted range. This multi-center approach eliminated center-specific artifacts as a viable explanation.

    Paper 2: 25-Year Longitudinal Analysis (Smawfield 2025c, DOI: 10.5281/zenodo.17517141) — Extended temporal baseline to 25.3 years (2000–2025) using CODE products, analyzing 165.2 million station pairs from 474 unique receivers. This long-span analysis enabled detection of previously inaccessible signatures including 18.6-year lunar nutation coupling (R² = 0.641), CMB frame alignment (18.2° from dipole, 5,570× variance ratio over galactic motion), and confirmed temporal stability of the correlation structure with λ = 4,201 ± 1,967 km. The 25-year baseline provided access to long-period geophysical phenomena impossible to detect in shorter studies.

    Paper 3: Raw RINEX Validation (Smawfield 2025d, DOI: 10.5281/zenodo.17860166) — Definitive processing artifact elimination through analysis of raw GNSS observations using Single Point Positioning with broadcast ephemerides only. Analysis of 1.17 billion pair-samples across 539 stations (2022–2024) achieved 100% detection rate across 72 independent metric combinations (t-statistics up to 112, Cohen's d up to 0.304), confirming the signal exists in fundamental observables independent of sophisticated processing chains. This raw data validation eliminated the most plausible remaining artifact hypothesis.

    Paper 4: Cosmological Implications (Smawfield 2025e, DOI: 10.5281/zenodo.17982540) — Developed the connection between TEP and gravitational lensing, demonstrating how temporal-field gradients can produce apparent dark matter phenomenology. This paper introduced the Isochrony Axiom as a critical but untested assumption in standard lensing analysis, showed how its violation by TEP produces "Phantom Mass" through differential proper-time accumulation, and established the Earth-Galaxy scaling argument connecting the ~4,000 km GNSS correlation to ~50 kpc dark matter halos through Vainshtein screening (rV ∝ M1/3). This cosmological extension transforms TEP from a terrestrial metrology curiosity into a framework with profound implications for dark matter and cosmological tensions.

    The present manuscript integrates these five studies into a unified narrative, demonstrating that the TEP framework made specific quantitative predictions before empirical analysis, that three independent empirical investigations using progressively rigorous validation methods all confirm the predicted signatures, and that the framework's cosmological extension offers testable predictions for gravitational lensing and dark matter phenomenology. This theory-prediction-validation-cosmology sequence, with each empirical study addressing distinct artifact hypotheses (processing-specific, temporal transience, sophisticated processing) and the cosmological paper extending implications beyond terrestrial scales, follows the standard structure of scientific confirmation and theoretical development.

## 1.4 Prediction Timeline

    The TEP framework was developed independently of the GNSS analysis, with specific quantitative predictions established in prior theoretical work (Smawfield 2025a). The predicted correlation length range of λ = 1,000–10,000 km was derived from environmental screening theory before any GNSS data analysis was conducted.

    The chronology is significant: the TEP theory paper was published in August 2025 with quantitative predictions (DOI: 10.5281/zenodo.16921911); multi-center GNSS analysis was conducted from September through November 2025; the 25-year CODE longspan analysis was completed in November 2025; and raw RINEX validation confirmed the signal in fundamental observables in December 2025. The observed correlation length of λ = 4,201 ± 1,967 km falls within the predicted range, representing an a priori prediction confirmed by subsequent observation—the standard by which physics has validated theories from Einstein's perihelion precession to the Higgs boson mass.

## 1.5 The Temporal Equivalence Principle

    TEP extends Einstein's geometrization of space to the rate of time itself. The principle states that all non-gravitational dynamics, signals, and quantum phases evolve according to the proper time defined by a single causal metric g̃μν that couples universally to matter. The rate at which proper time accrues is governed by a scalar field φ through a conformal factor A(φ). In local freely falling frames, physics reduces to special relativity with invariant c; globally, the flow of time becomes path-dependent.

    This framework preserves all established physics locally—the speed of light remains exactly invariant in any freely falling laboratory—while predicting novel global effects that manifest as correlations in distributed timing networks. Crucially, the conformal sector of TEP, which governs clock-rate modulation, remains unconstrained by multi-messenger observations such as GW170817, which bound only disformal (cone-tilting) effects.

## 1.6 Scope and Non-Claims

    To avoid misinterpretation, the following claims are explicitly NOT made by this work:

    - Detection of a fifth force acting on test masses (the network responds to geometric configuration, not gravitational amplitude; GM/r² and GM/r³ scaling are both null)

    - Violation of local Lorentz invariance (the speed of light remains exactly invariant in any freely falling laboratory; TEP preserves local special relativity)

    - Breakdown of the equivalence principle in laboratory tests (universal coupling ensures the weak equivalence principle is preserved exactly; the Nordtvedt parameter η ~ 4 × 10⁻⁶ is well below LLR bounds)

    - Modification of light or gravitational wave propagation speeds (GW170817 constrains |cγ − cg|/c ≲ 10⁻¹⁵; TEP satisfies this through negligible disformal coupling B(φ) ≈ 0)

    - Violation of energy-momentum conservation (total stress-energy is conserved; apparent non-conservation in the Einstein frame is resolved by transformation to the Jordan frame where matter follows geodesics of g̃μν)

    TEP is a conformal scalar-tensor theory that modulates the rate of proper time accumulation while preserving all local physics. The GNSS signal reflects spatial correlations in clock rates arising from the scalar field's two-point correlation function, not violations of established physical principles.

## 1.7 Claim Hierarchy

#### Tier 1: Empirical Claim (Primary)

        GNSS clock residuals exhibit a persistent, distance-structured, anisotropic spatial correlation with exponential decay length λ = 4,201 ± 1,967 km. This correlation is robust across three independent processing centers (CODE, IGS, ESA; R² = 0.92–0.97), persists over 25.3 years spanning multiple satellite constellation changes, and appears in raw RINEX observations processed with Single Point Positioning (100% detection rate across 72 metric combinations, t-statistics up to 112). The correlation exhibits seven independent signatures: exponential spatial decay, East-West/North-South anisotropy (ratio 2.16), orbital velocity coupling (r = −0.888, 5.1σ), CMB frame alignment (18.2° from dipole), planetary event responses (56/156 significant), 18.6-year lunar nutation coupling (R² = 0.641), and semiannual nutation coupling (R² = 0.904).

#### Tier 2: Interpretive Claim (Secondary)

        The observed correlation structure is inconsistent with known atmospheric effects (ionosphere, troposphere), instrumental systematics (geomagnetic storms, altitude dependence), classical gravitational effects (GM/r² force scaling, GM/r³ tidal scaling), solar activity (r ≈ 0 with solar flux), and processing artifacts (persists in raw SPP data). The signal's selectivity profile—sensitive to velocity-dependent dynamics and CMB frame alignment while blind to mass scaling and solar rotation—distinguishes it from conventional systematic errors. However, the possibility of a presently unmodeled GNSS systematic tied to fundamental constellation geometry cannot be definitively excluded without independent non-GNSS validation (Section 3.6).

#### Tier 3: Theoretical Hypothesis (Tertiary)

        The observations are consistent with the Temporal Equivalence Principle (TEP), a bi-metric scalar-tensor framework in which proper time is modulated by a dynamical conformal field A(φ) = exp(2βφ/MPl). TEP made quantitative predictions before empirical analysis (λ = 1,000–10,000 km, exponential decay, velocity-dependent anisotropy, absence of mass scaling), which are confirmed by the data. The framework preserves local Lorentz invariance while predicting global path-dependent synchronization. This hypothesis is testable through independent experiments: SLR validation (Section 6.3.1.1), multi-constellation cross-validation, and triangle holonomy test (Section 6.3.3). The validity of Tier 1 (empirical observations) does not depend on the correctness of Tier 3 (TEP interpretation).

    This paper establishes Tier 1 through extensive validation, argues for Tier 2 through systematic exclusion of conventional explanations, and tests consistency with Tier 3 through a priori predictions. Independent replication is essential for all three tiers.

## 1.8 Paper Structure

    This paper presents evidence for TEP through systematic analysis of GNSS timing data. Section 2 characterizes the phenomenology of the observed signal. Section 3 validates the signal against artifact hypotheses through multiple independent tests. Section 4 presents the theoretical framework and its predictions. Section 5 develops the cosmological implications, including connections to gravitational lensing and dark matter phenomenology. Section 6 specifies explicit falsification criteria and outlines the experimental path to a definitive verdict.

    Falsification is central to this work. Concrete rejection criteria are specified, including failure of independent groups to replicate the signal in raw carrier-phase data, correlation length falling outside the 500–20,000 km range, confirmation that the signal arises from ephemeris artifacts rather than physical clock correlations (via Satellite Laser Ranging validation), and null synchronization holonomy in closed-loop triangular time-transfer experiments.

## 2. Phenomenology

    This section presents the empirical findings without theoretical interpretation. The signal characteristics are established through three complementary empirical studies, each designed to address distinct validation requirements:

    Paper 1 (Multi-Center Validation, 2023–2025): Cross-validation across three independent analysis centers (CODE, IGS, ESA) using 62.7 million station-pair measurements from 364 stations. This study established that the correlation structure is not processing-specific, achieving R² = 0.92–0.97 consistency across independent pipelines.

    Paper 2 (25-Year Longitudinal Analysis, 2000–2025): Extended temporal baseline using CODE products with 165.2 million station pairs from 474 unique receivers. This long-span analysis enabled detection of previously inaccessible signatures including 18.6-year lunar nutation coupling and CMB frame alignment, confirming temporal stability over 25.3 years.

    Paper 3 (Raw RINEX Validation, 2022–2024): Processing artifact elimination through analysis of raw GNSS observations using Single Point Positioning with broadcast ephemerides only. Analysis of 1.17 billion pair-samples across 539 stations achieved 100% detection rate across 72 independent metric combinations, confirming the signal exists in fundamental observables.

    The convergence of findings across these three independent methodologies—each addressing different potential artifact sources—provides robust empirical foundation for the reported signatures.

## 2.1 Spatial Correlation Structure

    Clock frequency residuals exhibit systematic spatial structure that persists after removal of all modeled relativistic effects. The phase coherence between station pairs decays exponentially with geodesic distance according to:

    $C(r) = A \cdot \exp(-r/\lambda) + C_0$

    where r denotes the inter-station distance, λ the characteristic correlation length, A the amplitude, and C0 the baseline offset. Fits to distance-binned means across approximately 28 bins achieve R² = 0.920–0.970 across all three analysis centers.

    The correlation metric employed is the magnitude-weighted phase alignment index, computed via cross-spectral density analysis in the 10–500 μHz band (corresponding to periods of 33 minutes to 28 hours). This approach measures whether clock fluctuations are in phase regardless of amplitude—information that survives GNSS processing because network adjustment removes common-mode offsets while preserving differential phase structure. Magnitude weighting ensures that frequency bins with stronger cross-spectral power contribute proportionally more to the phase average, so the metric reflects genuine correlated signals rather than noise.

    The primary result from the CODE 25-year analysis is a correlation length of λ = 4,201 ± 1,967 km. Cross-center validation yields consistent values: CODE reports λ = 3,330 km (95% CI: 1,198–5,918 km); IGS reports λ = 4,549 km (95% CI: 3,197–4,871 km); and ESA reports λ = 3,758 km (95% CI: 2,532–3,984 km). The coefficient of variation across centers is 18.2%, indicating consistency across independent processing pipelines and substantially reducing the likelihood of center-specific artifacts.

## 2.2 Seven Independent Signatures

    The 25-year longitudinal analysis identifies seven statistically independent signatures, summarized in Table 1. Their joint occurrence under the null hypothesis has probability p ≈ 2×10−27, exceeding 10σ significance.

| Signature | Observed Value | Significance |
| --- | --- | --- |
| Exponential Decay | λ = 4,201 ± 1,967 km | R² = 0.92–0.97 |
| Spatial Anisotropy | EW/NS = 2.16 | p < 10−15 |
| Orbital Velocity Coupling | r = −0.888 | 5.1σ (0/5M surrogates) |
| CMB Frame Alignment | 18.2° from CMB dipole | 5,570× variance ratio |
| Planetary Events | 56/156 significant (≥2σ) | 2.8× above null rate |
| 18.6-Year Nutation | R² = 0.641 | p < 10−8 |
| Semiannual Nutation | R² = 0.904 | p < 10−20 |

## 2.3 CMB Frame Alignment

    The most striking result concerns the alignment of the anisotropy vector with the Cosmic Microwave Background dipole. To guard against confirmation bias, the analysis was conducted blind to the locations of known cosmological vectors. A multi-resolution grid search tested 65,341 independent sky directions, identifying the best-fit direction (RA = 186°, Dec = −4°) solely by maximizing the anisotropy ratio R(θ,φ) in the GNSS data. Only after the vector was fixed was it compared to the CMB dipole (RA = 168°, Dec = −7°).

    The angular separation between the GNSS anisotropy vector and the CMB dipole is 18.2°. The raw p-value is p < 10−15, and after Bonferroni correction for the 65,341 trials (accounting for the look-elsewhere effect), the alignment remains significant at greater than 5.8σ.

    Bootstrap resampling analysis refines the significance estimate by accounting for spatial correlation in the grid search. The 65,341 tested directions are not statistically independent due to spatial smoothing in the anisotropy field. Resampling 10,000 bootstrap iterations with block sizes matched to the correlation length of the anisotropy ratio field (approximately 20° angular scale) yields an effective number of independent trials Neff ≈ 800–1,200, substantially smaller than the naive 65,341. With Neff ~ 1,000, the Bonferroni-corrected significance remains at 4.2σ (p ~ 10−5 after correction), confirming that the CMB alignment is robust to conservative multiple-comparison accounting. The 18.2° separation falls within the bootstrap 95% confidence interval (±22°), indicating that the observed alignment is consistent with the CMB dipole within statistical uncertainty.

    Within the TEP framework, this alignment admits a natural interpretation: the CMB frame corresponds to the cosmological rest frame where scalar field gradients ∇φ are minimized. Earth's motion at v ≈ 369 km/s through this frame creates a velocity-dependent anisotropy in the effective screening length λ, thereby modulating the clock correlation structure.

    Reference frame discrimination strongly favors the CMB interpretation over alternatives. The CMB frame yields r = 0.747 (R² = 55.7%), whereas the Solar Apex yields r = 0.007 (R² = 0.01%)—a variance ratio of 5,570×. Beyond statistical preference, the Solar Apex hypothesis is geometrically incompatible with the observations: the Solar Apex vector (+30° Dec) predicts North-South anisotropy, directly contradicting the observed East-West dominance. Only the CMB vector (−7° Dec) is geometrically consistent with the observed anisotropy plane. Hemisphere asymmetry corroborates this conclusion: Southern stations show r = −0.79 while Northern stations show r = +0.25, matching the geometry expected from the CMB velocity vector.

## 2.4 Planetary Event Responses

    Analysis of 156 planetary conjunction and opposition events reveals a pattern that constrains the underlying mechanism. The detection rate is 56/156 events showing significant response at the 2σ level or above, with 25 surviving Bonferroni correction (16.0%). Random dates show 5.5× smaller effect sizes (Mann-Whitney p < 10−17). Notably, however, neither GM/r² (force) nor GM/r³ (tidal gradient) scaling is detected (all p > 0.5).

    This absence of mass scaling is informative. Mercury produces comparable response rates (42.5%) to Jupiter (34.8%) despite approximately 106× smaller gravitational influence at Earth. Classical tidal forcing would produce mass-dependent signatures scaling with GM/r³ (tidal gradient), while direct gravitational acceleration would scale as GM/r² (force). The observed absence of both scalings (all p > 0.5) distinguishes the phenomenon from conventional gravitational effects—neither tidal deformation nor fifth-force acceleration. The pattern is consistent with a threshold-dependent or geometric (alignment-driven) mechanism rather than a continuous mass-scaling effect.

    The elevated detection rate (2.8× above null) combined with the absence of mass scaling constrains the mechanism to geometric rather than energetic coupling. Critically, the absence of GM/r² scaling persists in raw RINEX data processed with RTKLIB Single Point Positioning—a method that processes each station independently without network adjustment or common-mode filtering. This rules out the hypothesis that the null scaling is merely an artifact of CODE/IGS processing and supports the interpretation that the TEP coupling is intrinsically geometric: the network responds to planetary configuration geometry rather than gravitational field strength.

### 2.4.1 Proposed Mechanism: Velocity-Field Gradient Modulation

    The pattern of planetary event responses—elevated detection rate without mass scaling—suggests a specific physical mechanism consistent with TEP's velocity-dependent framework. The hypothesis is that planetary alignments modulate Earth's motion through the scalar field φ gradient, creating transient changes in the effective screening length that depend on geometric configuration rather than gravitational amplitude.

    **Physical Mechanism:** Earth moves through the φ field at velocity v⊕ ≈ 30 km/s (orbital) + 369 km/s (CMB frame). Planetary configurations create time-varying perturbations to the local field gradient ∇φ through their own screening halos. When planets align with Earth's velocity vector relative to the CMB frame, constructive interference enhances the velocity-dependent anisotropy; when perpendicular, destructive interference suppresses it. The effect depends on the alignment angle θ between the planetary position and Earth's CMB velocity vector:

    $\Delta\lambda(\theta) \propto \lambda_0 \left[1 + \epsilon \cos\theta\right]$

    where ε is a small coupling parameter (~10−3) and θ is measured from the CMB dipole direction. This predicts maximum modulation when planets are aligned with or opposite to the CMB velocity vector (θ = 0° or 180°), and null modulation when perpendicular (θ = 90°).

    **Why Mass Scaling is Absent:** The mechanism operates through geometric alignment rather than gravitational force. Each planet creates a local perturbation to ∇φ within its own screening radius rV ∝ M1/3. However, the Earth-planet distance (AU scale) far exceeds any planetary rV (Mercury: ~100 km; Jupiter: ~5,000 km), placing Earth in the unscreened regime where the field perturbation has decayed to a threshold-dependent residual. The detection depends on whether this residual exceeds a critical threshold for modulating Earth's local screening length, not on the perturbation's amplitude. This threshold behavior explains why Mercury (smaller M but closer) produces comparable detection rates to Jupiter (larger M but farther).

    **Testable Prediction:** Event detection rates should correlate with the alignment angle between planetary position and the CMB dipole direction. Specifically, events with |θ| < 30° (aligned) should show detection rates >60%, while events with 60° < |θ| < 120° (perpendicular) should show rates <30%. Analysis of the 156 planetary events stratified by alignment angle provides a falsifiable test of this geometric hypothesis.

## 2.5 Selectivity Profile: Inertial Interferometer

    The pattern of detections and non-detections characterizes the network's sensitivity and provides powerful constraints on viable physical mechanisms. Table 2 summarizes this selectivity profile.

| Observable | Result | Interpretation |
| --- | --- | --- |
| Orbital velocity | r = −0.888, 5.1σ | Sensitive to kinematic dynamics |
| Semiannual nutation | R² = 0.904 | Sensitive to inertial orientation |
| 18.6-year nutation | R² = 0.641 | Sensitive to long-period dynamics |
| CMB-frame motion | r = 0.747 | Sensitive to cosmic rest frame |
| Planetary alignments | 2.8× above null | Sensitive to geometric configuration |
| GM/r² and GM/r³ scaling | Both null (p > 0.5) | Insensitive to gravitational amplitude (force or tidal) |
| Solar rotation (27-day) | Null (r < 0.09) | Insensitive to solar surface phenomena |
| Lunar standstill | Null | Insensitive to static geometry |

    This selectivity profile is highly discriminating. The network responds to velocity-dependent and orientation-dependent phenomena while remaining blind to gravitational amplitude and solar surface effects. This pattern suggests that the GNSS network functions as an inertial interferometer—measuring velocity-dependent correlation geometry—rather than as a gravimeter measuring Newtonian force.

    The distinction is fundamental. A gravimeter would show GM/r² (force) scaling, while a differential tidal sensor would show GM/r³ (gradient) scaling; in both cases, Mercury's influence would be negligible compared to Jupiter's. An inertial interferometer responds to geometric configuration regardless of mass—precisely what is observed. The CMB frame alignment, with its 5,570× variance ratio over galactic motion, identifies the operative kinematic reference at cosmic scales: Earth's motion at approximately 369 km/s through the cosmological rest frame.

## 2.6 Statistical Power

    The raw RINEX validation analysis achieves extraordinary statistical power. Across 1.17 billion pair-samples, the directional anisotropy yields t-statistics up to 112 with Cohen's d up to 0.304. The signal is detected in 72 of 72 independent metric combinations (100% detection rate) with mean R² = 0.93.

    Hemisphere asymmetry provides independent corroboration. Southern stations show stronger orbital velocity coupling (r = −0.79) than Northern stations (r = +0.25), matching the geometry expected from Earth's motion through the CMB frame. Both hemispheres independently show East-West exceeding North-South correlation (Northern: 1.20×, Southern: 1.35×), ruling out hemisphere-specific artifacts.

    Monthly consistency is equally striking. East-West exceeds North-South in 94–100% of all 36 months analyzed, with short-distance coherence ratios showing coefficient of variation below 1%. The underlying signal is constant; the annual modulation in full-distance correlation lengths reflects atmospheric screening effects that vary seasonally.

## 3. Validation

    Extraordinary claims require extraordinary validation. This section systematically addresses the most plausible artifact hypotheses, demonstrating that the signal survives each challenge.

## 3.1 Processing Artifact Hypothesis

    The most immediate concern is that the signal arises from sophisticated processing chains—network adjustments, Kalman filtering, or integer ambiguity resolution. To address this possibility, raw RINEX observation files were processed using Single Point Positioning with broadcast ephemerides only, representing the simplest possible processing chain and one entirely independent of the network adjustments employed by CODE, IGS, and ESA.

    The raw RINEX validation dataset comprises 539 stations over three years (2022–2024), totaling 1.17 billion pair-samples. The signal is detected in all 72 metric combinations tested (100% detection rate), with mean R² = 0.93. Directional anisotropy matches the CODE findings, with East-West exceeding North-South by 2–22% at short distances (<500 km). A critical audit confirms this is robust to distance bias: E-W pairs are 13 km longer than N-S pairs (a bias *against* the signal), and distance-matching strengthens the ratio.

    The correlation length in raw data (λ ≈ 700–1,100 km) is shorter than in precise products (λ ≈ 4,000 km), consistent with ionospheric noise masking the long-range structure. When ionospheric effects are removed (Ionofree mode), the raw correlation length increases, supporting a "Ladder of Precision" where the signal scale converges toward the CODE benchmark as noise is mitigated. Orbital velocity coupling is detected at r = −0.763 (5.4σ) in multi-GNSS mode, and CMB alignment yields RA = 188°, Dec = −5° (20.0° from the CMB dipole), consistent with the 25-year analysis. The signal therefore exists in the fundamental observables and cannot be attributed to processing artifacts.

## 3.2 Ionospheric Hypothesis

    If the signal were ionospheric in origin—driven by solar activity modulating electron density along signal paths—it should correlate with solar cycle variations. The 25-year baseline spans Solar Cycles 23, 24, and 25, encompassing three solar maxima (2001, 2014, 2024) and two deep minima, providing substantial leverage for this test.

    The correlation between the signal and the Solar Flux Index is r ≈ 0, indicating no significant relationship. The correlation length λ remains stable across solar maxima within uncertainty bounds. Moreover, the signal persists in the dual-frequency ionosphere-free L1+L2 combination (λ = 1,069 km, R² = 0.969), and excluding high-ionosphere periods improves the signal by 21–23%. This last observation indicates that the ionosphere suppresses rather than creates the correlation. The signal is therefore not ionospheric in origin.

## 3.3 Geomagnetic Hypothesis

    If the signal were driven by geomagnetic storms affecting clock electronics or signal propagation, it should vary systematically with geomagnetic activity. Stratification by GFZ Kp index data across all 72 metric combinations yields median Δλ ≈ −1% at the primary threshold (Kp < 3 versus Kp ≥ 3), with 60 of 72 tests falling within ±5% across all station filters and processing modes. No systematic deviation appears at higher storm thresholds (Kp ≥ 4, Kp ≥ 5). The signal is therefore geomagnetically invariant.

## 3.4 Tropospheric Hypothesis

    If the signal were tropospheric in origin—driven by water vapor or pressure variations that correlate spatially—stations at high altitude (with thinner atmospheric columns) should show systematically different correlation lengths than sea-level stations.

    Analysis of 360 regressions spanning global and latitude-controlled altitude quintiles reveals that only 3.1% exhibit p < 0.05 slopes, consistent with chance expectation. The Q5/Q1 ratio (highest to lowest altitude quintile) has median 0.97 (IQR 0.76–1.27), indicating that low-altitude and high-altitude stations yield statistically indistinguishable λ values.

    The tropospheric column depth varies by approximately 20–30% between sea-level and high-altitude stations. If the correlation structure were tropospheric, one would expect systematic increase in λ at high altitude, latitude dependence correlated with tropopause height, and seasonal modulation correlated with monsoon patterns. None of these signatures are observed. The signal is altitude-invariant, consistent with a phenomenon operating at gravitational potential scales (~6,400 km) rather than atmospheric scales (~10 km).

    **Comprehensive Validation Summary:** The altitude independence test represents one component of a broader validation framework spanning 388 independent statistical tests across multiple analysis families. False Discovery Rate correction (FDR-BH: 52.3%, Hierarchical Empirical Bayes: 39.7%, Bonferroni: 40.0%) demonstrates robustness against multiple comparison artifacts. This extensive validation, combined with cross-center consistency (coefficient of variation 18.2% across CODE, IGS, ESA), establishes that the observed correlations survive rigorous scrutiny across independent processing pipelines, environmental conditions, and statistical frameworks.

## 3.5 Long-Period Stability

    If the signal were a transient artifact of the analysis window, it would not persist across the full 25.3-year baseline or couple coherently to long-period geophysical phenomena.

    The data reveal clear coupling to multiple geophysical cycles. The 18.6-year lunar nutation shows R² = 0.641 (p < 10−8), with 1.4 complete cycles observed. The semiannual nutation shows R² = 0.904 (p < 10−20), with over 50 complete cycles. The Chandler wobble (433 days) shows R² = 0.096 with phase stability of 0.72 across more than 21 cycles. Core signatures—correlation length, anisotropy, and orbital coupling—remain consistent across the full 25-year span. The signal is therefore temporally stable across multiple geophysical cycles.

## 3.6 Most Plausible Remaining Conventional Explanation

    The strongest remaining conventional hypothesis is that the observed correlations arise from a presently unmodeled, long-baseline GNSS systematic that survives both network adjustment and raw SPP processing. This section addresses this hypothesis with intellectual honesty: while extensive validation constrains conventional explanations, it cannot definitively exclude all possible systematics.

    **The Hypothesis:** A sophisticated, distance-dependent systematic error exists in GNSS observations that: (1) produces exponential decay with λ ≈ 4,000 km; (2) exhibits East-West/North-South anisotropy ratio of 2.16; (3) couples to Earth's orbital velocity with r = −0.888; (4) aligns with the CMB dipole to 18.2°; (5) responds to planetary configurations without mass scaling; (6) couples to 18.6-year lunar nutation and semiannual nutation; (7) survives independent processing by CODE, IGS, and ESA; (8) persists in raw RINEX data processed with broadcast ephemerides only; and (9) appears identically across 72 independent metric combinations with 100% detection rate.

    **Constraints from Validation:** Each validation step constrains this hypothesis without claiming impossibility:

    - **Multi-center consistency (R² = 0.92–0.97, CV = 18.2%):** The systematic must be present in three independent processing pipelines using different algorithms, software, and analysis centers. This excludes center-specific artifacts but does not exclude systematics inherent to GNSS observation geometry or satellite constellation design.

    - **Raw RINEX validation (100% detection, t-stats up to 112):** The systematic must survive Single Point Positioning with broadcast ephemerides, which processes each station independently without network adjustment. This excludes network-adjustment artifacts but does not exclude systematics in the broadcast ephemeris itself (the "ephemeris loophole" addressed in Section 6.3.1.1).

    - **25-year temporal stability:** The systematic must persist across satellite constellation changes, hardware upgrades, and processing methodology evolution spanning 1.4 complete cycles of the 18.6-year lunar nutation. This excludes transient artifacts but does not exclude systematics tied to fundamental GNSS geometry.

    - **CMB frame alignment (5,570× variance ratio):** The systematic must preferentially align with the cosmological rest frame rather than galactic motion, solar motion, or ecliptic plane. This is geometrically difficult to explain through satellite constellation geometry but not impossible if some unrecognized coupling exists between GNSS observation geometry and Earth's motion through the CMB frame.

    **Remaining Possibility:** The most plausible remaining conventional explanation is a systematic tied to the fundamental geometry of GNSS satellite constellations that couples to Earth's motion through space in a way not previously recognized. Such a systematic would need to produce all observed signatures through a single coherent mechanism—a requirement that becomes increasingly implausible as the number of independent signatures grows, but cannot be definitively excluded without independent validation through non-GNSS methods (SLR, triangle holonomy test).

    **Path to Resolution:** Definitive exclusion of this hypothesis requires: (1) SLR validation (Section 6.3.1.1) to close the ephemeris loophole; (2) multi-constellation cross-validation showing identical λ across GPS, GLONASS, Galileo, and BeiDou; and (3) triangle holonomy test (Section 6.3.3) providing direct measurement of synchronization non-closure. Until these tests are completed, the possibility of an unrecognized GNSS systematic, however implausible given the convergence of evidence, cannot be completely excluded.

## 3.7 Preliminary SLR Evidence

    Satellite Laser Ranging provides an independent validation pathway for TEP predictions. While dedicated SLR analysis remains a priority for future work (Section 6.3.1.1), existing anomalies in the SLR literature exhibit signatures consistent with the geometric screening effects predicted by TEP.

### 3.7.1 The Blue-Sky Effect

    SLR stations systematically measure Earth's crust as depressed during clear-sky (high atmospheric pressure) conditions compared to microwave-based techniques. This phenomenon, termed the "Blue-Sky effect," produces vertical displacements of 1.1–4.4 mm across the global SLR network, with the average effect of 1.1 mm corresponding to a scale discrepancy of approximately 0.2 ppb relative to Earth's radius. The maximum displacement of 4.4 mm occurs at inland stations in Eastern Europe and Central Asia.

    The standard explanation attributes this effect to atmospheric pressure loading—elastic deformation of Earth's crust under varying atmospheric mass. This explanation is well-established in the geodetic literature and successfully accounts for the observed magnitudes. The Blue-Sky effect is noted here not as evidence against the standard model, but as an observation whose phenomenology (systematic vertical offsets correlated with atmospheric conditions) is also consistent with TEP's prediction of path-dependent proper time accumulation. Distinguishing between these interpretations would require dedicated analysis beyond the scope of this preliminary survey.

### 3.7.2 Persistent Range Bias Residuals

    Despite decades of refinement, high-performance SLR stations exhibit persistent systematic residuals at the 10–20 mm level that resist conventional explanation. Historical analysis of the NASA MOBLAS network documented centimeter-level systematic biases through collocation analysis, with survey uncertainties and apparent system delays showing discrepancies up to 20 mm between calibration targets. These residuals persist even after corrections for tropospheric delay, relativistic effects, tidal loading, and hardware calibration.

    The International Laser Ranging Service Analysis Standing Committee continues to investigate station-satellite-specific range biases for LAGEOS satellites, noting that range biases correlate strongly with station height over short periods—a correlation that may influence SLR's ability to provide absolute scale. The magnitude of these unexplained residuals (10–20 mm) falls within the range predicted by TEP for proper time delays accumulated through Earth's gravitational screening field.

### 3.7.3 Elevation-Dependent Systematics

    SLR residuals exhibit elevation- and azimuth-dependent errors that standard tropospheric mapping functions fail to fully remove. At low elevation angles, these errors can reach 20 mm after mapping from zenith values. Recent work has demonstrated that estimating tropospheric biases absorbs these elevation-dependent errors but does not explain their physical origin.

    From the TEP perspective, elevation dependence maps directly to path length through Earth's screening field. A signal decaying with characteristic length λ ≈ 4,200 km would produce precisely the elevation-dependent curve that standard atmospheric models fail to capture—the residual structure represents not tropospheric error but geometric screening accumulated along the signal path.

### 3.7.4 Interpretation and Caveats

    These SLR anomalies are presented as preliminary supporting evidence rather than independent confirmation. The geodetic community has proposed conventional explanations for each phenomenon: atmospheric pressure loading for the Blue-Sky effect, hardware calibration issues for range biases, and tropospheric modeling limitations for elevation dependence. TEP does not claim these explanations are incorrect—rather, that the same observations admit an alternative interpretation consistent with geometric screening of proper time.

    Definitive SLR validation requires dedicated analysis of LAGEOS residuals stratified by elevation angle, with explicit comparison to TEP predictions. The publicly available data from stations such as Yarragadee (Australia) provide immediate opportunity for such analysis. If TEP is correct, plotting SLR residuals against elevation angle should reveal the characteristic exponential decay structure with λ consistent with GNSS-derived values.

    **Key References:** Sośnica et al. (2013), Journal of Geodesy 87(8):751-769; Luceri et al. (2019), Journal of Geodesy 93:2357-2379; Strugarek et al. (2022), EGU General Assembly, EGU22-9834.

## 3.8 Why the Signal Was Not Previously Identified

    The correlation structure documented here has been present in GNSS data since the network's inception but was systematically removed by processing algorithms optimized for positioning accuracy rather than fundamental physics characterization. Three factors explain why it was not identified earlier.

    Standard GNSS processing is designed to eliminate correlations, not to characterize them. The processing chain includes network adjustment (which removes common-mode signals), Kalman filtering (which smooths residuals toward zero-mean white noise), reference frame alignment (which forces consistency with ITRF, absorbing systematic offsets), and outlier rejection (which discards data points deviating from expected behavior).

    The correlation structure detected in this analysis is precisely what these algorithms are designed to remove. From the perspective of positioning accuracy, this is appropriate—the signal degrades positioning solutions. From the perspective of fundamental physics, however, it constitutes signal that reveals underlying structure.

    The present analysis examines the residuals that standard processing discards, employing phase-coherent cross-spectral methods that preserve correlation structure. The signal was present throughout the history of GNSS operation but was categorized as systematic error rather than physical phenomenon.

## 3.9 Summary of Excluded Hypotheses

    Table 3 summarizes the artifact hypotheses tested and their outcomes.

| Hypothesis | Test | Result |
| --- | --- | --- |
| Processing artifact | Raw RINEX/SPP | Signal persists (100% detection) |
| Ionospheric | 25-year solar cycle | No correlation (r ≈ 0) |
| Tropospheric | Altitude stratification | No altitude dependence (Q5/Q1 = 0.97) |
| Geomagnetic | Kp stratification | Near-invariant (Δλ ≈ −1%) |
| Tidal forcing | GM/r² and GM/r³ scaling | Neither force nor tidal scaling (p > 0.5) |
| Random noise | Shuffle test | 30× R² ratio (real vs. shuffled) |
| Center-specific bias | Multi-center validation | CV = 18.2% across CODE/IGS/ESA |

    Each row represents an independent test. The signal survives all challenges, constraining viable explanations to physical mechanisms operating at planetary scales with velocity-dependent, CMB-aligned geometry.

## 4. Theoretical Framework

    Having established that the signal is robust and not attributable to known artifacts, this section presents the theoretical framework that predicted these observations. The Temporal Equivalence Principle is a principled extension of general relativity within the well-established class of scalar-tensor theories, with specific quantitative predictions that preceded the GNSS analysis.

## 4.1 Two-Metric Geometry

    TEP posits that spacetime is endowed with two metrics related by a scalar field φ:

    $\tilde{g}_{\mu\nu} = A(\phi)\, g_{\mu\nu} + B(\phi)\, \nabla_\mu \phi\, \nabla_\nu \phi$

    The gravitational metric gμν determines geodesics of test masses and gravitational wave propagation. The matter metric g̃μν determines atomic clock rates, electromagnetic propagation, and all Standard Model physics. The conformal factor A(φ) = exp(2βφ/MPl) provides isotropic rescaling of proper time, while the disformal factor B(φ) introduces anisotropic light-cone tilts and is tightly constrained by GW170817.

    This structure belongs to the well-established class of scalar-tensor theories developed by Bekenstein (1993) and Damour & Polyakov (1994). The contribution of TEP lies in its operational interpretation: the conformal factor A(φ) governs the rate at which proper time accumulates, rendering time itself a dynamical field.

## 4.2 Physical Interpretation

    The framework is formulated as a scalar-tensor theory in the Einstein frame where gravity remains canonical. The scalar field φ couples to matter through the conformal factor A(φ), modulating the rate of proper time accumulation. The coupling strength is parameterized by β ~ 10−3, weak enough to evade existing constraints yet strong enough to produce observable GNSS correlations. The scalar field obeys a Klein-Gordon equation with source terms from matter density and a screening potential V(φ) that suppresses fifth forces in dense environments. (Full action and field equations are provided in Appendix A.)

## 4.3 Foundational Axioms

    The framework rests on four axioms. The first (two-metric structure) states that gravity is described by gμν while matter couples universally to g̃μν. The second (temporal equivalence) states that all non-gravitational processes evolve according to proper time dτ defined by g̃μν, and that in local freely falling frames, physics reduces to special relativity with invariant c. The third (causal safety) requires that cg = cγ within GW170817 bounds (|cg − cγ|/c ≲ 10−15), constraining B(φ) to be negligible at late times. The fourth (screening) posits a chameleon-like potential V(φ) that yields density-dependent effective mass, screening fifth forces in dense environments while permitting cosmological dynamics.

## 4.4 Compton Energy Scale

    The observed correlation length λ corresponds to a specific scalar field mass via the Compton relation λ = ℏ/(mφc). For λ = 3,330–4,549 km, the effective field mass is:

    $m_\phi \approx (4.34\text{--}5.93) \times 10^{-14}\, \text{eV}/c^2$

    This mass scale is consistent with environmental screening mechanisms where the effective field mass varies with local matter density. The apparent inconsistency with existing precision tests—which are typically sensitive to scales around 10−15 eV/c²—is resolved by screening: the field mass increases in dense environments where most precision tests are conducted, while remaining at the observed value in the sparse terrestrial environment sampled by the GNSS network.

### 4.4.1 Vainshtein Screening and the Terrestrial Correlation Length

    The Vainshtein mechanism provides a quantitative explanation for the observed GNSS correlation length. Derivative self-interactions in the scalar field Lagrangian (cubic Galileon terms) generate nonlinear screening that suppresses fifth forces near massive objects. The screening radius rV scales as:

    $r_V = \left(\frac{GM}{\Lambda^3}\right)^{1/3}$

    This M1/3 scaling is universal, applying from planetary to galactic masses without additional free parameters. For Earth with M = M⊕ ≈ 6 × 1024 kg and taking Λ ~ 10−13 eV (the dark energy scale, consistent with cosmological observations):

    $r_V \sim \left(\frac{G \times 6 \times 10^{24}\text{ kg}}{(10^{-13}\text{ eV})^3}\right)^{1/3} \approx 4{,}000\text{ km}$

    This prediction matches the observed GNSS correlation length λ = 4,201 ± 1,967 km with no adjustable parameters. The physical mechanism operates as follows:

    - **Inside rV:** Nonlinear derivative interactions dominate, screening the scalar field and suppressing fifth forces. The effective coupling to matter is reduced by factors of (r/rV)3, rendering TEP effects undetectable in laboratory tests and solar system observations.

    - **Outside rV:** The scalar field becomes weakly coupled and mediates long-range correlations. Clock networks sampling this region detect the φ-field structure through spatial covariance.

    - **Transition scale:** The correlation length λ ~ rV marks the boundary where screening becomes ineffective, allowing the scalar field's influence to manifest.

    **Falsifiable Predictions:** The M1/3 scaling law makes specific predictions for intermediate-mass systems:

    - **Jupiter** (M ~ 300 M⊕): Predicted λJupiter ~ (300)1/3 × 4,000 km ≈ 27,000 km. Testable with Juno mission clock data or future Jupiter atmospheric probe networks.

    - **Globular clusters** (M ~ 105 M☉): Predicted screening length ~ 1–10 pc, testable through pulsar timing array correlations within clusters.

    - **Galaxies** (M ~ 1012 M☉): Predicted rV ~ 50 kpc, matching the observed dark matter halo radius. This connection is developed in Section 5.5.

    The Vainshtein mechanism thus provides a unified explanation spanning 15 orders of magnitude in mass scale, from Earth to galaxy clusters, with a single fundamental parameter Λ set by cosmological observations.

## 4.5 Synchronization Holonomy and Spatial Correlations

    TEP predicts synchronization holonomy—the non-closure of proper time when transported around a closed loop:

    $\mathcal{H} \equiv \oint_C d\tau_{\text{prop}} \neq 0$

    In standard general relativity, after subtracting known effects (Sagnac, Shapiro, gravitational redshift), this holonomy vanishes. In TEP with disformal coupling B ≠ 0, it does not.

    An important subtlety arises in the conformal-only limit (B = 0): the time-transport one-form Ωμ = ∂μ(½ ln A) is exact, its curl vanishes, and holonomy is zero. Nonzero holonomy requires non-exact structure from disformal couplings or more general non-metricity. This raises a question: if GW170817 constrains B(φ) to be negligible, and conformal coupling alone produces no holonomy, what generates the observed GNSS correlations?

    The resolution lies in distinguishing two observables. Synchronization holonomy (H) represents the non-closure of proper time around closed loops; this requires B ≠ 0 and is the target of proposed triangle tests. Spatial correlation structure, by contrast, represents the statistical covariance of clock residuals as a function of distance; this arises from the spatial structure of the φ field itself—specifically, the two-point correlation function ⟨φ(x)φ(x+r)⟩ ∝ exp(−r/λ)—and persists even when B = 0.

    The GNSS signal reflects the latter mechanism: clocks at different locations sample different values of the conformal factor A(φ), and spatial correlations in φ induce correlations in clock rates. The exponential decay with distance reflects the Compton wavelength of the screened scalar field, while the velocity-dependent anisotropy arises because Earth's motion through the φ field creates a preferred direction in the correlation structure. The correlation structure is thus a statistical signature of the φ field's spatial configuration, whereas holonomy is a geometric signature of non-integrable time transport. Both are predictions of TEP, but they probe different aspects of the framework and have different observational requirements.

### 4.5.1 Resolving the Conformal Paradox

        **Figure 1:** Synchronization Holonomy: The Definitive Test of Non-Integrable Proper Time. Three atomic clocks (A, B, C) positioned at vertices of a triangle (1,000–3,000 km separation) exchange one-way time-transfer pulses around the loop. Panel A (Standard GR): After subtracting all known relativistic effects (Sagnac, Shapiro, gravitational redshift, Lense-Thirring), the accumulated time difference returns to zero—proper time is integrable. Panel B (TEP): A residual time gap ℋresid ~ 10−18–10−16 s remains, scaling as (B/A)|∇φ|² · Area. Target sensitivity: σ(ℋ)/Tloop ≤ 10−19.

    A critical objection must be addressed: if TEP employs only a conformal coupling A(φ), how can it produce synchronization holonomy? In standard differential geometry, a conformal transformation preserves angles and therefore cannot generate path-dependent effects in closed loops. The resolution lies in distinguishing synchronization holonomy from spatial correlation structure.

    **Observable 1: Synchronization Holonomy (ℋ)** — The non-closure of proper time around closed loops, defined as ℋ ≡ ∮C dτprop. This observable requires B(φ) ≠ 0 (disformal coupling) because only non-exact structure in the time-transport connection generates path-dependent accumulation. In the conformal-only limit (B = 0), the one-form Ωμ = ∂μ(½ ln A) is exact, its curl vanishes identically, and holonomy is zero by Stokes' theorem. The triangle holonomy test (Section 6.3.3) targets this observable and will definitively measure whether B(φ) is nonzero at the 10−19 fractional level.

    **Observable 2: Spatial Correlation Structure (C(r))** — The statistical covariance of clock frequency residuals as a function of inter-station distance. This observable requires only A(φ) ≠ 1 (conformal coupling) and persists even when B = 0. The mechanism is fundamentally different from holonomy: clocks at different spatial locations sample different values of the scalar field φ(x), and the two-point correlation function ⟨φ(x)φ(x+r)⟩ ∝ exp(−r/λ) of the underlying field induces correlations in clock rates through the conformal factor A(φ) = exp(2βφ/MPl). The exponential decay length λ = ℏ/(mφc) reflects the Compton wavelength of the screened scalar field.

    **Physical Analogy:** Consider a field of thermometers distributed across a landscape measuring a temperature field T(x) that varies spatially due to underground thermal sources. Even in the absence of heat flow (analogous to B = 0), nearby thermometers show correlated readings because they sample similar temperature values. The correlation length reflects the spatial structure of the temperature field itself, not thermal transport between locations. Similarly, GNSS clocks sample the φ-field "temperature," and their correlation structure reveals the field's spatial configuration.

    **Why GW170817 Does Not Constrain the GNSS Signal:** The multi-messenger constraint |cγ − cg|/c ≲ 10−15 bounds the disformal sector B(φ), which affects propagation speeds. The GNSS correlation structure arises from the conformal sector A(φ), which affects clock rates but preserves null cones. These are distinct physical effects constrained by different observables. GW170817 tells us that light and gravity travel at the same speed; it does not constrain whether clocks at different locations tick at correlated rates due to sampling a common scalar field.

    **Observational Distinction:** The GNSS signal is a "thermometer correlation"—a statistical signature of φ-field structure detected through spatial covariance analysis. Holonomy is a "transport non-closure"—a geometric signature requiring closed-loop, one-way time transfer. Both are TEP predictions, but they probe different sectors (conformal vs. disformal) and require different experimental geometries (distributed networks vs. closed triangles). The present work establishes the former; the latter awaits dedicated triangle tests.

## 4.6 Prediction Comparison

    TEP made specific quantitative predictions in August 2025 (Smawfield 2025a) before the GNSS analysis was conducted. Table 3 compares these forecasts to observations.

| Prediction | Forecast | Observed | Status |
| --- | --- | --- | --- |
| Correlation length | 1,000–10,000 km | 4,201 ± 1,967 km | Consistent |
| Spatial decay form | Exponential | ΔAIC = 12.8 vs. power-law | Consistent |
| Velocity-dependent anisotropy | Modulates with v/c | r = −0.888 (5.1σ) | Consistent |
| Geometric alignment | EW > NS in ecliptic frame | EW/NS = 2.16 | Consistent |
| Broadband coupling | Not frequency-selective | Detected across 10–500 μHz | Consistent |
| CMB frame alignment | Cosmological rest frame | 18.2° from CMB dipole | Consistent |
| GM/r² and GM/r³ scaling | Both absent (filtered) | No force or tidal scaling (p > 0.5) | Null confirmed |

    The null prediction regarding GM/r² scaling is particularly discriminating. Classical tidal forcing would produce mass-dependent signatures scaling with GM/r³, with Mercury contributing far less than Jupiter. The observed absence of such scaling (all p > 0.5 across multiple metrics) persists across both CODE network-adjusted products and raw RINEX data processed with RTKLIB Single Point Positioning. Since SPP processes each station independently without network adjustment or common-mode filtering, the persistence of the null in SPP data rules out the hypothesis that the absence of mass scaling is merely a processing artifact. The TEP coupling appears to be intrinsically geometric rather than energetic—responding to alignment configuration regardless of mass. Mercury (42.5% detection rate) produces comparable responses to Jupiter (34.8%) despite ~106× smaller gravitational influence, consistent with a threshold-dependent or alignment-driven mechanism rather than continuous force scaling.

## 4.7 Why Existing Tests Do Not Constrain TEP

    Existing precision tests probe different observables than TEP's primary signatures, as summarized in Table 4.

| Test | Observable | Reason TEP Is Unconstrained |
| --- | --- | --- |
| GW170817 | |cγ − cg|/c | Constrains disformal B(φ), not conformal A(φ) |
| Cassini (PPN-γ) | Two-way Shapiro delay | Reciprocity-even; does not bound loop holonomy |
| Resonator MM/KT | Closed-path anisotropy | Even-parity; blind to one-way non-reciprocity |
| GPS operation | Internal consistency | Uses explicit GR+Sagnac modeling; verifies model |
| Clock redshift tests | Pairwise comparisons | Confirms GR locally; closed loops required |

    The critical distinction is that existing tests probe two-way or pairwise measurements, whereas TEP's holonomy is a closed-loop, one-way observable to which these tests are structurally blind. Two-way measurements probe reciprocal paths where holonomy cancels by symmetry; pairwise clock comparisons measure local potential differences but cannot detect non-integrability, which requires closed loops. TEP effects are therefore structurally invisible to these geometries, explaining why decades of precision tests found nothing amiss while the signal was present in GNSS residuals.

#### Operational Measurement Protocol: What Existing Tests Probe

        The distinction between observables that constrain TEP and those that do not can be formalized through the synchronization holonomy (after GR subtraction), where σ is the time-transport one-form whose line integral equals the calibrated proper-time increment along a leg:

        $\mathcal{H}_{\text{resid}} = \oint_{\partial\Sigma} (\tilde{\sigma} - \sigma_{\text{GR}}) = \iint_{\Sigma} (\tilde{F} - F_{\text{GR}})$

        This represents the loop non-closure of time transport built from measured proper times (in time units; right-hand rule on Σ). It is gauge and synchronization-invariant and vanishes in standard relativity and in the conformal-only limit of TEP. GR subtraction includes Sagnac, Lense-Thirring gravito-magnetic effects, Shapiro delay, and gravitational redshift, with ITRF ephemerides and TT/TDB standards.

        How flagship constraints map to ℋresid:

        - GW170817 (GW-EM coincidence): |cγ − cg|/c ≲ 10−15 constrains global cone splits. In TEP, late-time conformal coupling preserves null cones, so electromagnetic and gravitational waves share causal structure; small disformal tilts today are allowed. This is a boundary condition, not a loop-holonomy test.

        - Cassini (PPN-γ): Two-way Doppler/Shapiro is reciprocity-even; it calibrates σGR to subtract but does not bound ℋresid.

        - Resonator MM/KT tests: Cavities bound closed-path, even-parity (two-way sums) anisotropy at 10−17–10−18, yet are blind to odd-parity (direction-reversing one-way differences) non-reciprocity and loop non-closure—the ingredients of ℋresid.

        - GPS operation: Network self-consistency uses explicit GR+Sagnac modeling and largely two-way/common-view calibration. This verifies internal consistency under the assumed GR model; it is not a direction-reversing one-way loop-closure null.

        - Clock redshift and pairwise A↔B tests: Exquisitely confirm GR locally; only closed loops (A→B→C→A with direction reversal) can reveal non-integrability captured by ℋresid.

        Why classical tests can be null while ℋresid ≠ 0: (1) Conformal null-cone invariance implies no large GW-EM kinematic delays (consistent with GW170817). (2) If ∂tφ = 0 over loop timescale and gradients are conservative, there is no first-order one-way anisotropy. Along a fixed path, forward and backward times cancel at order (B,∇φ); leading effects are second order or require time dependence and kinematics. Thus two-way and closed-path nulls can hold while a loop-holonomy test remains sensitive.

        Experimental falsifier (primary endpoints): Run a closed-loop, one-way time-transfer (and/or portable-clock) test and report: (1) Leg-wise antisymmetry: ΔtAB = tAB − tBA (and optionally ΞAB ≡ (tAB − tBA)/(tAB + tBA)). (2) Loop holonomy ℋresid after subtracting GR (Sagnac/gravito-magnetic/Shapiro). Use triangle/quadrilateral geometries with direction reversal; extend with interplanetary one-way links and multi-species clock networks.

### 4.7.1 Compatibility with Binary Pulsar Timing and Lunar Laser Ranging

    Two high-precision tests warrant explicit discussion: binary pulsar timing (PSR J0737-3039, testing GR to 0.05%) and Lunar Laser Ranging (testing the equivalence principle to 10−13). Both are compatible with TEP through screening and universal coupling.

    Binary Pulsar Timing: The Hulse-Taylor pulsar PSR B1913+16 and the double pulsar PSR J0737-3039 provide exquisite tests of general relativity through orbital decay measurements. The observed orbital period derivatives match GR predictions to better than 0.05%, constraining deviations in the strong-field regime. TEP compatibility requires that scalar field effects are screened in the vicinity of neutron stars.

    The Vainshtein screening radius for a neutron star (M ~ 1.4 M☉) is:

    $r_V^{\text{NS}} \sim \left(\frac{GM_{\text{NS}}}{\Lambda^3}\right)^{1/3} \sim 10^4 \text{ km}$

    This screening radius far exceeds the neutron star radius (RNS ~ 10 km) and encompasses the binary orbit (a ~ 106 km for PSR B1913+16). Within rV, the scalar field is strongly screened, suppressing fifth-force contributions by factors of (RNS/rV)3 ~ 10−9. The residual TEP effects on orbital dynamics are therefore below the 0.05% observational precision, rendering binary pulsar tests insensitive to the conformal coupling that generates GNSS correlations.

    Critically, pulsar timing measures orbital dynamics (gravitational sector, metric gμν), not clock correlations (matter sector, metric g̃μν). The two-metric structure allows the gravitational dynamics to remain GR-like while clock rates exhibit conformal modulation—these are distinct observables probing different sectors of the theory.

    Lunar Laser Ranging: LLR measures the Earth-Moon distance to millimeter precision, testing the equivalence principle through the Nordtvedt parameter η. The current bound η < 4.4 × 10−4 constrains deviations from the strong equivalence principle. TEP satisfies this constraint through universal coupling: the conformal factor A(φ) couples identically to all matter fields, preserving the equivalence principle in the matter frame.

    The key distinction is between weak and strong equivalence principles. The weak EP (universality of free fall for test masses) is preserved exactly in TEP because all matter couples to the same metric g̃μν. The strong EP (self-gravitating bodies fall identically to test masses) can show deviations in scalar-tensor theories through composition-dependent couplings. However, TEP's universal conformal coupling ensures that the Nordtvedt parameter:

    $\eta = 4\beta^2 - \gamma - 3$

    remains within observational bounds for the coupling strength β ~ 10−3 required to explain the GNSS correlation length. With γ ≈ 1 (Cassini constraint) and β ~ 10−3, the predicted Nordtvedt parameter is η ~ 4 × 10−6, well below the LLR bound.

    Furthermore, both Earth and Moon are within each other's screening radii (rVEarth ~ 4,000 km, Earth-Moon distance ~ 384,000 km), placing the system in a regime where screening effects partially suppress the scalar field influence. The combination of universal coupling and screening ensures LLR compatibility while preserving the GNSS signal in the terrestrial environment where screening is less effective.

## 4.8 Relation to Established Scalar-Tensor Theory

    The TEP framework is not a novel field theory but rather a phenomenological application of well-established scalar-tensor gravity to precision timing. TEP with B = 0 and A(φ) = exp(2βφ/MPl) reduces to Brans-Dicke gravity in the Jordan frame, with the coupling β related to the Brans-Dicke parameter ωBD. The screening mechanism employed (Vainshtein) arises naturally in the cubic Galileon sector of Horndeski gravity—the most general scalar-tensor theory with second-order equations of motion, ensuring freedom from Ostrogradsky instabilities. TEP's environmental screening is analogous to chameleon mechanisms but operates through derivative self-interactions rather than density-dependent mass. The two-metric structure is philosophically similar to Bekenstein's TeVeS, though TEP is more minimal, employing only a scalar field rather than scalar plus vector.

    The novelty of TEP lies not in the field theory but in the observational interpretation—specifically, the recognition that conformal coupling creates spatial correlations in clock rates that standard GNSS processing preserves while filtering amplitude information.

## 4.9 The Processing Filter Mechanism

    A critical insight explains both why the signal was missed for three decades and why it survives in the data. Standard GNSS processing is designed to optimize positioning accuracy by removing correlations, not to characterize them. The processing chain includes network adjustment (which removes common-mode signals), Kalman filtering (which smooths residuals toward zero-mean white noise), reference frame alignment (which forces consistency with ITRF, absorbing systematic offsets), and outlier rejection (which discards data points deviating from expected behavior).

    This architecture creates a natural filter: energetic (common-mode) signals are suppressed while geometric (differential) signals are preserved. Classical gravitational effects manifest as common-mode offsets—all clocks in a region experience similar tidal acceleration—and are therefore filtered. The TEP signal, by contrast, manifests as differential phase structure—the correlation between clock pairs depends on their geometric relationship—and is therefore preserved.

    This mechanism makes a specific prediction: GM/r² scaling should be absent in processed data because it is common-mode, while geometric signatures (anisotropy, CMB alignment) should survive because they are differential. The observations confirm both predictions. The processing filter is not a bug but a feature: it explains why the signal was categorized as systematic error rather than physical phenomenon, and why it can now be detected by examining the residuals that standard processing discards.

### 4.9.1 Mathematical Demonstration of Filter Selectivity

    The processing filter mechanism can be demonstrated mathematically by decomposing clock residuals into common-mode and differential components. Consider a network of N stations with clock residuals δti(t) after initial processing. Standard GNSS processing applies a datum constraint that forces residuals to have zero spatial mean, removing common-mode offsets. In Fourier space, this constraint eliminates the k = 0 (monopole) component while preserving all k ≠ 0 (dipole, quadrupole, and higher) components that have zero spatial mean.

    Classical gravitational effects manifest as either GM/r² (force) or GM/r³ (tidal gradient) scaling. Tidal effects are k ≈ 2 (quadrupole) in spherical harmonic decomposition. However, standard GNSS processing explicitly models and subtracts known tidal terms (solid Earth tides, ocean loading, pole tide) using IERS conventions with sub-millimeter precision. The datum constraint removes residual common-mode offsets, while explicit tidal modeling eliminates the structured quadrupole components.

    The TEP signal survives not merely because it is differential (k ≠ 0), but because it represents an unmodeled scalar correlation function with exponential spatial form C(r) ∝ exp(−r/λ). This functional form does not match the spherical harmonic basis (Ylm) used in standard tidal modeling. Tidal corrections are parameterized as sums of spherical harmonics with known frequencies (diurnal, semidiurnal, long-period); the TEP correlation is a distance-dependent phase coherence that is not decomposable into the standard tidal basis. The processing removes what it models (tidal harmonics) but preserves what it does not model (exponential distance correlations in phase space).

    The TEP correlation structure manifests as phase coherence between station pairs. The analysis uses cross-spectral methods that are invariant under common-mode amplitude rescaling because they depend only on relative phases, not absolute amplitudes. The datum constraint removes monopole amplitude information but preserves the phase structure encoded in higher-order spatial modes. The exponential decay C(r) ∝ exp(−r/λ) reflects the two-point correlation function of the underlying φ field, which has nonzero spatial gradients and therefore survives network adjustment. (Detailed mathematical formulation is provided in Appendix B.)

    **Prediction Verification:** This mathematical framework makes two testable predictions, both confirmed by observations:

    - **Null Prediction:** Classical gravitational effects (GM/r² force and GM/r³ tidal gradient) should be absent in processed data because they are predominantly k = 0 and k ≈ 0 (long-wavelength) modes that the datum constraint removes. Observed: Neither GM/r² nor GM/r³ scaling detected (all p > 0.5).

    - **Positive Prediction:** Geometric anisotropy (EW/NS ratio) should survive because it represents a k ≠ 0 quadrupole moment with zero spatial mean but nonzero phase structure. Observed: EW/NS = 2.16 detected at p < 10−15.

    The persistence of both predictions in raw RINEX data processed with Single Point Positioning (which applies no network adjustment) confirms that the filter mechanism is a property of the signal's spatial structure, not an artifact of CODE/IGS processing algorithms.

## 4.10 Systematic Error Budget for Triangle Test

    For the proposed triangle holonomy test to achieve definitive discrimination, systematic errors must be controlled below the target sensitivity of 10−19 fractional. Table 5 presents the error budget.

| Systematic | Mitigation | Residual |
| --- | --- | --- |
| Sagnac (Earth rotation) | IERS Earth orientation parameters | ≤ 10−21 |
| Shapiro (GR potentials) | Precise ephemerides + atmospheric models | ≤ 3 × 10−20 |
| Troposphere/ionosphere | Dual-frequency + water vapor radiometry | ≤ 5 × 10−20 |
| Fiber dispersion | Bidirectional calibration + thermal control | ≤ 5 × 10−20 |
| Link noise + clock instability | Optical clocks + day averaging | ≤ 3 × 10−19 |
| Total budget |  | ≤ 1 × 10−19 |

    For parameters saturating current bounds, the predicted signal is H/Tloop ~ 10−16–10−18. The error budget demonstrates that a definitive test is achievable with current technology.

## 5. Consistency with Cosmological Phenomenology

#### Epistemic Status

        The following section does not constitute an independent detection of cosmological effects. Instead, it addresses the question: if the terrestrial clock correlations reported in Sections 2–4 arise from a conformally coupled scalar field with Vainshtein screening, then what are the necessary consequences for astrophysical observations that already exist? The empirical GNSS result stands independently of the validity of this extension.

        This section demonstrates that the same scalar field producing λ ≈ 4,000 km correlations on Earth would generically contribute to gravitational lensing observables in a manner that is degenerate with particulate dark matter under the Isochrony Axiom—the assumption that all photons in a lensed image represent a synchronous snapshot. The framework makes specific, testable predictions that provide clean discriminators between particulate dark matter and conformal time-field interpretations at astrophysical scales.

    Within the TEP framework, the same conformal factor A(φ) that modulates clock correlations on Earth would affect light propagation across cosmological distances. The universal M1/3 Vainshtein scaling law connects terrestrial and galactic scales without additional free parameters: the screening radius rV = (GM/Λ³)1/3 scales from Earth (λ ≈ 4,000 km) to galaxies (≈50 kpc) across 15 orders of magnitude in mass. This quantitative continuity argument constrains the theory space and provides falsifiable predictions for gravitational lensing observations.

## 5.1 GW170817 Constraints

    The simultaneous detection of gravitational waves (GW170817) and gamma rays (GRB 170817A) from a binary neutron star merger at approximately 40 Mpc is widely cited as constraining modified gravity theories. The signals arrived within Δt = 1.74 ± 0.05 seconds after traveling some 130 million years, yielding:

    $\frac{|c_{\gamma} - c_g|}{c} \lesssim 10^{-15}$

    This constraint is often interpreted as ruling out theories in which light and gravity propagate differently. However, this interpretation conflates two distinct effects that must be carefully distinguished.

    The disformal sector B(φ) tilts light cones, causing cγ ≠ cg, and is indeed tightly constrained by GW170817. The conformal sector A(φ), by contrast, rescales proper time isotropically: photons and gravitational waves share the same null geodesic and experience identical time dilation. Because both messengers originated from the same coordinate and traversed the same path through the scalar field, any conformal time dilation is common-mode and cancels in the differential measurement. GW170817 therefore constrains cone tilting but remains blind to clock-rate modulation. The dark matter phenomenology in TEP arises entirely from the conformal sector, which remains unconstrained by multi-messenger observations.

    This can be seen explicitly in the conformal metric limit (B = 0):

    $ds^2 = A(\phi) g_{\mu\nu} dx^\mu dx^\nu = 0$

    Since A(φ) is an overall prefactor, the null condition ds² = 0 is identical for both the matter metric g̃μν (governing photons) and the gravitational metric gμν (governing gravitational waves). Their speeds are identical: cγ = cg = c/√A(φ). Both signals are delayed by the scalar field by exactly the same amount, and the differential arrival time is zero.

## 5.2 The Isochrony Axiom and Observational Degeneracy

    Standard gravitational lensing analysis assumes that all photons in a lensed image represent a synchronous snapshot of the source. This Isochrony Axiom—the assumption that differential time delays arise solely from geometric path length differences (Fermat potential), not from temporal-field gradients along the paths—is implicit in the lens equation and has never been empirically tested at the precision required to distinguish these scenarios.

    Within the TEP framework, this axiom would generically fail in the presence of conformal metric couplings. Photons traveling through regions with different scalar field values φ(x) would accumulate proper time at different rates. The observed arrival time difference would include both geometric delay and temporal-field contribution, creating an observational degeneracy with particulate dark matter under standard analysis assumptions:

    Consider two interpretations of the same Fermat potential surface. Under the standard interpretation, which assumes isochrony, an Einstein ring with apparent convergence κobs exceeding baryonic mass implies a dark matter halo with MDM = Mobs − Mbaryons. Under the TEP interpretation, which rejects isochrony, the observed image is a temporal composite: photons arriving simultaneously left the source at different epochs. For evolving sources, this temporal depth projects onto the image plane as apparent spatial distortion—apparent mass that does not correspond to matter.

    Both frameworks fit the data equally well. The difference is interpretive, depending on which axiom is taken as fundamental. The existence of dark matter as particulate matter is a conclusion contingent on the Isochrony Axiom.

## 5.3 The Phantom Mass Mechanism (Conditional)

    Within a two-metric framework where the Isochrony Axiom breaks down, a mechanism for apparent dark matter emerges naturally. If the conformal factor A(φ) varies spatially—forming halo-like configurations around mass concentrations—then photons taking different paths through the lens accumulate different proper times:

    $\Delta \tilde{\tau} = \frac{1}{c} \int (\sqrt{A(\phi)} - 1)\, dl$

    For an evolving source, this differential delay creates a temporal composite image: photons arriving simultaneously at the detector left the source at different epochs. Standard lens reconstruction, which assumes isochrony, interprets this temporal smearing as spatial distortion—apparent mass that does not exist.

    The amplification matrix in TEP decomposes as:

    $\mathcal{A}_{ij} = \underbrace{\mathcal{A}^{\text{GR}}_{ij}}_{\text{Baryonic}} + \underbrace{\mathcal{A}^{(\alpha,\text{static})}_{ij}}_{\text{Phantom Mass}} + \underbrace{\mathcal{A}^{(\alpha,\text{dyn})}_{ij}}_{\text{Shear Noise}}$

    Standard analyses attribute the sum of the first two terms to total mass. TEP identifies the second term as Phantom Mass—a geometric effect of the scalar field's refractive gradient rather than particulate matter.

    Order-of-magnitude estimates support this mechanism. For a coupling strength ε ≈ 10−6 and halo-scale integration (Lhalo ≈ 2 Mpc), the differential delay across an Einstein radius is:

    $\Delta \tilde{\tau} \sim \frac{\epsilon}{2} \frac{L_{\text{halo}}}{c} \approx 3 \text{ years}$

    This estimate is consistent with observed strong-lens time delays and resolves what might be called the "Refsdal Paradox": naive cosmological integration of conformal delays would predict millennia-long gaps between lensed images, which are not observed. The resolution requires careful distinction between static and dynamic delay components.

### 5.3.1 Resolving the Refsdal Paradox

    The apparent contradiction between TEP's predicted ~3 year delays and observed day-to-month scale time delays in strong lensing systems (e.g., SN Refsdal: Δt ~ 1 year between images) resolves through two complementary mechanisms that distinguish absolute from differential delays.

    **Mechanism 1: Static versus Dynamic Decomposition** — The amplification matrix decomposition (Equation above) separates temporal effects into distinct components with different observational signatures:

    - **Static Refraction (𝒜(α,static)):** The absolute proper-time delay Δτstatic ~ years accumulated across the entire halo depth affects all light rays equally, independent of which image they form. This component is source-independent and manifests as an overall convergence κ that standard lens reconstruction attributes to total mass. Because it is common to all images, it does not contribute to differential time delays between multiple images.

    - **Dynamic Shutter (𝒜(α,dyn)):** The gradient of the delay field ∇(Δτ) creates differential delays that depend on the ray's impact parameter and the source's evolutionary state. For transient sources with variability timescale tvar, the dynamic contribution scales as Δτdynamic ~ ∇(Δτ) × vsource × tvar, where vsource is the source's proper motion or internal velocity field.

    **Observational Consequence:** Time delays measured between multiple images of transient events (supernovae, AGN flares) probe only the differential dynamic component, not the absolute static delay. The ~3 year estimate applies to the total integrated delay affecting convergence measurements, while day-to-month delays reflect the gradient-driven differential component.

    **Mechanism 2: Screening Locality** — The conformal factor A(φ) varies significantly only within the Vainshtein screening radius rV ~ 100 kpc around the lens galaxy. Beyond this radius, A(φ) → 1 and delays cease accumulating. The integration path for the static delay is therefore limited to the halo depth Lhalo ~ 2 Mpc (the physical extent of the dark matter halo), not the full cosmological distance DL ~ Gpc to the source. This locality ensures that:

    $\Delta \tau_{\text{static}} \sim \frac{\epsilon}{2} \frac{L_{\text{halo}}}{c} \sim 3\text{ years} \quad \text{(not } \frac{D_L}{c} \sim 10^9\text{ years)}$

    **Quantitative Prediction for Time-Delay Cosmography:** TEP predicts systematic offsets in H0 measurements from time-delay cosmography of 5–10%. Standard analysis assumes that measured time delays Δtobs between images directly constrain the Fermat potential differences ΔΦ via Δt = (1+zL)ΔΦ/c². However, if the static refraction component contributes an unmodeled offset Δτstatic to the convergence κ, the inferred H0 is biased:

    $H_0^{\text{inferred}} = H_0^{\text{true}} \left(1 + \frac{\Delta\tau_{\text{static}}}{\Delta t_{\text{obs}}}\right)^{-1}$

    For Δτstatic ~ 3 years and typical strong-lens delays Δtobs ~ 50 days, this predicts H0inferred ~ 0.95 H0true, a 5% systematic offset. This is consistent with the observed H0 tension between local measurements (Cepheid-calibrated: 73 km/s/Mpc) and time-delay cosmography (H0LiCOW: 73.3 km/s/Mpc) versus CMB inference (Planck: 67.4 km/s/Mpc).

    **Falsification Test:** If TEP is correct, time-delay cosmography measurements should show systematic residuals correlated with lens halo mass (larger halos → larger Δτstatic → larger H0 bias). Analysis of the H0LiCOW and TDCOSMO samples stratified by lens velocity dispersion (a proxy for halo mass) should reveal this correlation at the 3σ level if the Extended Regime is operative.

## 5.4 The Variability-Mass Correlation

    A critical discriminator between TEP and standard dark matter emerges from the dynamic term in the amplification matrix. Because the stochastic shear contribution depends on source proper motion μs through the delay gradient, the inferred lens mass should correlate with source variability timescale. This correlation is forbidden in standard CDM (mass is source-independent) but required in TEP's Extended Regime.

    Specifically, lenses observed through rapidly variable sources (AGN, quasars) should show systematically different mass reconstructions than the same lenses observed through slowly evolving sources (elliptical galaxies). The predicted effect is 10–20% higher inferred mass for AGN-backlit lenses compared to galaxy-backlit lenses at the same redshift.

    This test can be performed with existing data from HST strong-lens archives and SDSS quasar lenses. If existing catalogs show no statistically significant correlation between inferred lens mass and source variability class at the 3σ level, the Extended Regime is disfavored. If such a correlation exists, it constitutes positive evidence for TEP independent of GNSS timing or FRB observations.

## 5.5 The Earth-Galaxy Scaling Argument

    The TEP framework suggests a connection between terrestrial metrology and galactic dynamics. If the Vainshtein screening mechanism is operative, the same scalar field that produces λ ≈ 4,000 km on Earth should produce the dark matter halo at galactic scales—representing the same physics at different density scales.

    In scalar-tensor theories with derivative self-interactions (the Galileon/Horndeski class), the scalar field is screened near massive objects. The screening radius rV scales as:

    $r_V \propto \left(\frac{M}{M_{\text{Pl}}}\right)^{1/3} \cdot \Lambda^{-1}$

    where M is the central mass and Λ is the strong-coupling scale. This M1/3 scaling is universal, applying from planetary to galactic masses without free parameters.

    If TEP is correct, the ~4,000 km correlation on Earth (λ = 4,201 ± 1,967 km) and the ~50 kpc dark matter halo in galaxies represent the same scalar field at different density scales. This constitutes a falsifiable prediction: dark matter halos would not be arbitrary structures but would follow a strict density-dependent scaling law determined by the field's coupling parameters.

    The Vainshtein screening radius scales as rV ∝ (GM/Λ³)1/3. For Earth mass M⊕ ≈ 6 × 1024 kg and typical galaxy mass Mgal ≈ 1012 M☉ ≈ 2 × 1042 kg, the mass ratio is approximately 3 × 1017, yielding a predicted screening radius ratio of (3 × 1017)1/3 ≈ 7 × 105.

    The GNSS correlation length λ reflects the environmental screening length in Earth's local gravitational environment, while rhalo reflects the transition radius where screening becomes ineffective. These are related but distinct quantities: the M1/3 scaling applies to the transition radius, whereas the correlation length involves additional factors including local matter density and field gradient structure. The appropriate comparison is between the ratio of screening lengths at similar relative positions, predicting that the ratio of correlation lengths should scale as M1/3 when measured at equivalent fractional radii—a testable prediction for intermediate-mass systems such as stellar clusters and dwarf galaxies.

| Scale | Observed λ | Interpretation |
| --- | --- | --- |
| Earth (GNSS) | ~4,000 km | Local screening length |
| Galaxy (rotation curves) | ~10–50 kpc | Dark matter halo radius |
| Ratio | ~104–105 | Consistent with rV ∝ M1/3 |

    Intermediate-scale tests provide a path to falsification. Jupiter (M ~ 300 ME) should exhibit a correlation length of approximately 25,000 km, while globular clusters should show intermediate screening lengths of approximately 1–10 pc.

## 5.6 Emergent MOND Scale

    A notable prediction of TEP with Vainshtein screening is the natural emergence of the MOND acceleration scale. The transition from screened (GR-like) to unscreened (modified) behavior occurs when the scalar-mediated acceleration equals the Newtonian acceleration:

    $a_0 \sim \frac{\Lambda^2}{M_{\text{Pl}}} \sim c H_0$

    For Λ ~ 10−13 eV (the dark energy scale):

    $a_0 \approx (3 \times 10^8 \text{ m/s}) \times (2.2 \times 10^{-18} \text{ s}^{-1}) \approx 6.6 \times 10^{-10} \text{ m/s}^2$

    This value lies within a factor of 5 of the empirical MOND value (a0 ≈ 1.2 × 10−10 m/s²). The cosmic coincidence that MOND requires as unexplained input emerges automatically in TEP from the scalar field's strong-coupling scale being set by the Hubble parameter.

## 5.7 Observational Predictions

    TEP makes quantitative predictions for existing and near-future observations, summarized in Table 7.

| Observable | TEP Prediction | Dataset |
| --- | --- | --- |
| Lensed quasar time delays | 5–15% systematic residuals | H0LiCOW, TDCOSMO |
| Galaxy rotation curves | rt ∝ M1/3 | SPARC, THINGS |
| Weak lensing S8 | Galaxy surveys 3–5% low vs. CMB | DES, KiDS, HSC vs. Planck |
| Strong lens mass-variability | AGN lenses 10–20% higher mass | SLACS, BELLS |
| FRB lensing residuals | Achromatic delays 0.1–10 ms | CHIME/FRB repeaters |

    These predictions are falsifiable within 2–5 years using archival data and ongoing surveys. A single null result does not exclude TEP, but a pattern of null results across all five observables would constitute strong evidence against the framework.

## 5.8 Connection to Cosmological Tensions

    TEP offers potential resolution to persistent tensions in precision cosmology.

    Regarding the H0 tension: time-delay cosmography assumes standard GR time transport. If TEP's temporal structure contributes to lens delays, the inferred H0 would be biased, potentially explaining the ~5 km/s/Mpc discrepancy between local and CMB-derived values.

    Regarding the S8 tension: galaxy weak lensing measures dynamic sources, whereas CMB lensing measures a static backlight. TEP predicts source-dependent kinematic noise that would suppress inferred S8 from galaxy surveys relative to CMB, matching the observed ~5% discrepancy.

    These are predictions of the framework rather than post-hoc explanations: if time is dynamical, measurements that assume isochrony will show systematic biases correlated with source properties.

## 5.9 Relation to ΛCDM

    The standard ΛCDM cosmological model fits observational data with remarkable precision. This framework attributes gravitational lensing anomalies to cold dark matter particles and explains cosmological structure formation through hierarchical assembly of particulate halos. ΛCDM is empirically successful and theoretically well-motivated.

    TEP does not claim that ΛCDM is incorrect. Rather, TEP proposes that the same observational phenomena can be reinterpreted under a different axiom set. The key distinction is operational, not observational at current precision:

    - **ΛCDM interpretation:** Gravitational lensing mass reconstructions directly measure matter distribution. Excess convergence beyond baryonic mass indicates particulate dark matter with density profile ρDM(r). The Isochrony Axiom (all photons in an image represent a synchronous snapshot) is implicitly assumed and never tested.

    - **TEP interpretation:** Gravitational lensing mass reconstructions measure the effective Fermat potential, which includes both matter distribution and temporal-field gradients. Excess convergence can arise from differential proper-time accumulation along different light paths. The Isochrony Axiom is empirically testable through achromatic timing residuals in lensed transients.

    Both frameworks fit existing lensing data because they make identical predictions for static sources observed under the assumption of isochrony. The degeneracy is broken by observations that probe temporal structure:

    - **Lensed fast radio bursts:** ΛCDM predicts achromatic residuals (after plasma correction) consistent with zero. TEP predicts achromatic residuals at the 0.1–10 ms level arising from differential temporal-field accumulation (Section 5.4).

    - **Variable source lensing:** ΛCDM predicts inferred lens mass is independent of source variability timescale. TEP predicts systematic correlation between inferred mass and source variability, with AGN showing higher apparent mass than elliptical galaxies for the same lens (Section 5.4).

    - **CMB vs. galaxy lensing S8:** ΛCDM predicts convergence of CMB-inferred and galaxy-inferred S8 as systematics are controlled. TEP predicts persistent tension arising from kinematic noise in galaxy surveys that is absent in CMB (static source).

    The frameworks compete rather than contradict. ΛCDM remains the default interpretation given its theoretical maturity and empirical success. TEP offers an alternative interpretation that becomes distinguishable through specific tests designed to probe temporal structure in lensing observations. Independent of the cosmological extension, the terrestrial GNSS correlation structure documented in Sections 2–4 requires explanation within any framework.

    **Falsification asymmetry:** Null results in the tests above would exclude TEP's cosmological extension while leaving ΛCDM intact. Positive results would require reinterpretation of dark matter phenomenology but would not exclude particulate dark matter in all contexts—only in regimes where temporal-field gradients dominate the lensing signal. This asymmetry reflects TEP's status as a challenger framework rather than an established paradigm.

## 6. Falsification and Experimental Program

    A theory that cannot be falsified is not science. TEP makes specific, testable predictions with explicit null-result thresholds. This section specifies the criteria that would exclude the framework and outlines the experimental path to a definitive verdict.

## 6.1 Explicit Falsification Criteria

    The TEP framework is falsified if any of the following observations are confirmed:

#### Tier 1: Terrestrial Falsification (GNSS Signal)

        - **Independent replication failure:** Independent research groups fail to reproduce exponential decay (λ ≈ 3,000–5,000 km) using either raw carrier-phase observables (RINEX data with Single Point Positioning) or independent clock solutions (CLK files from CODE, IGS, ESA, or JPL). This would exclude the entire GNSS evidence base.

        - **Correlation length out of bounds:** Independently replicated correlation length falls outside λ = 500–20,000 km range. This would exclude the TEP screening model, which predicts λ ~ 4,000 km for Earth based on Vainshtein screening at Λ ~ 10⁻¹³ eV.

        - **Ephemeris loophole confirmed:** Satellite Laser Ranging (SLR) residuals show no distance-structured correlation when analyzed using identical methods applied to GNSS data. This would indicate the GNSS signal arises from broadcast ephemeris artifacts rather than physical clock correlations (Section 6.3.1.1).

#### Tier 2: Synchronization Holonomy (Direct Test)

        - **Null holonomy at high precision:** Closed-loop one-way time-transfer around triangular geometry (1,000–3,000 km legs) yields null synchronization holonomy ℋ/Tloop < 10⁻¹⁹ after subtracting all known relativistic effects (Sagnac, Shapiro, gravitational redshift, Lense-Thirring). This would constrain disformal coupling B(φ) below cosmologically interesting levels and exclude the non-integrable proper time prediction (Section 6.3.3).

#### Tier 3: Cosmological Falsification (Dark Matter Connection)

        - **CMB-Galaxy S8 convergence:** Next-generation surveys (CMB-S4 + LSST) show CMB-inferred and galaxy lensing-inferred S8 agreeing to better than 1σ, with no systematic dependence on source redshift or variability. This would exclude TEP's prediction of kinematic noise in galaxy surveys arising from temporal-field gradients. The current S8 tension (4.5σ) between CMB and galaxy lensing measurements represents existing anomalous data consistent with TEP predictions.

    **Complete Framework Falsification:** The combination of Tier 1 (independent replication failure OR ephemeris loophole confirmed) AND Tier 2 (null holonomy at high precision) AND Tier 3 (S8 tension resolved without TEP mechanism) would exclude TEP as a viable framework rather than merely constraining one regime. However, failure of any single Tier 1 criterion is sufficient to exclude the GNSS evidence base, which forms the empirical foundation of the framework.

## 6.2 Consequences of Failed Replication

    If independent groups cannot replicate the raw RINEX signal, this would indicate one of several possibilities: a subtle processing artifact specific to the present pipeline that was not identified, a statistical fluctuation in the particular dataset selection, or an unidentified systematic in the cross-spectral analysis method.

    In such a scenario, TEP would lose its primary empirical support. The theoretical framework would remain mathematically consistent but empirically unmotivated, and the cosmological predictions would become untestable speculation. This outcome would represent the scientific method working as intended.

## 6.3 Three-Tier Experimental Strategy

    TEP can be confirmed or excluded through a coordinated experimental program organized by urgency and discriminating power.

### 6.3.1 Tier 1: Immediate Replication

    The most urgent priority is independent replication. If independent research groups cannot replicate the raw RINEX signal, the GNSS evidence base is excluded and TEP loses its primary empirical support.

    Required tests include independent clock analysis using either raw carrier-phase data (RINEX files with SPP processing) or processed clock solutions (CLK files from independent analysis centers), multi-constellation cross-check (verifying that the signal appears identically in GPS, GLONASS, Galileo, and BeiDou), and cross-center consistency (confirming CODE/IGS/ESA/JPL agreement with CV < 20%). Both RINEX and CLK approaches provide valid replication pathways; the key requirement is independence from the original analysis pipeline.

    These tests can be performed by any institution with GNSS processing capability, including IGS, NIST, PTB, SYRTE, and JPL. The publicly available data (both raw RINEX observations and processed CLK solutions) enable independent verification of the λ ≈ 4,000 km exponential decay, and the analysis can be completed in weeks using standard GNSS processing tools.

#### 6.3.1.1 Addressing the Ephemeris Loophole

    A legitimate concern exists regarding the Single Point Positioning validation: SPP relies on broadcast ephemerides generated by control segments using network adjustments. If the artifact originates in the broadcast orbit determination, SPP would reproduce rather than eliminate it. This represents a potential systematic loophole requiring independent validation.

    Satellite Laser Ranging (SLR) provides the definitive test. SLR is a purely geometric, optical measurement independent of the GNSS microwave clock and ephemeris loop. SLR stations measure satellite range using laser pulses with millimeter precision, providing orbit determination completely independent of GNSS processing chains. If the GNSS correlation structure arises from ephemeris artifacts, it should be absent in SLR residuals; if it reflects physical clock correlations, it should appear in both.

    Proposed Tier 1.5 Test: Correlate GNSS clock residuals with SLR range residuals for co-located or nearby stations. The International Laser Ranging Service (ILRS) operates approximately 40 SLR stations globally, with several co-located with IGS GNSS stations. Analysis would proceed as follows: (1) Extract SLR range residuals after precise orbit determination; (2) Compute cross-correlation between SLR residuals and GNSS clock residuals for station pairs as a function of distance; (3) Test whether the exponential decay structure (λ ≈ 4,000 km) appears in SLR data.

    Expected Outcome: If TEP is correct and the signal reflects physical clock-rate modulation, SLR residuals should show distance-structured correlations because the laser round-trip time depends on both the satellite orbit and the ground station clock. If the signal is an ephemeris artifact, SLR residuals (which use independent orbit determination) should show no correlation structure. This test definitively closes the ephemeris loophole and can be performed with existing ILRS data.

### 6.3.2 Tier 2: Direct Holonomy Test

    The triangle holonomy test provides the most direct discrimination between TEP and artifact hypotheses by measuring synchronization non-closure via closed-loop, multi-leg, one-way time-transfer. The proposed geometry involves two ground stations separated by 1,000–3,000 km plus a MEO satellite. The target sensitivity is σ(H)/Tloop ≤ 10−19 after averaging, with a signal forecast of H/Tloop ~ 10−16–10−18 for parameters saturating current bounds. This test requires a dedicated mission.

### 6.3.3 Tier 3: Cosmological Validation

    CMB-S4 combined with LSST weak lensing would test for source-dependent kinematic noise. TEP predicts that galaxy surveys (LSST) will infer S8 systematically 3–5% lower than CMB measurements (Planck, CMB-S4), with the discriminator being correlation with source redshift and variability. The current S8 tension (4.5σ) between CMB and galaxy lensing measurements represents existing anomalous data consistent with TEP predictions; resolution of this tension through conventional mechanisms would exclude TEP's kinematic noise explanation.

## 6.4 Path to Resolution

    Table 9 summarizes the experimental program and its decision structure.

| Tier | Test | Null-Result Threshold | What It Would Exclude |
| --- | --- | --- | --- |
| Tier 1 | Raw carrier-phase replication | Signal absent | Entire GNSS evidence base |
| Multi-constellation cross-check | Different λ across constellations | Universal coupling |
| SLR ephemeris validation | No correlation in SLR residuals | Physical clock origin |
| Tier 2 | Triangle holonomy test | Null < 10−19 | Holonomy mechanism |
| Tier 3 | CMB-S4 + LSST S8 | Agreement < 1σ | Kinematic noise |

    Tier 1 tests determine whether TEP survives as a viable hypothesis. If all three Tier 1 tests confirm the signal, Tier 2 tests discriminate among mechanisms. Tier 3 tests probe the cosmological regime. A definitive verdict—confirmation or exclusion—is achievable through this coordinated program.

## 7. Conclusions

    This paper presents empirical evidence for a persistent, distance-structured correlation in global atomic clock networks that probes an empirically untested assumption of general relativity: the global integrability of proper time.

## 7.1 Summary of Findings

    This integrated manuscript synthesizes results from a systematic five-paper research program that establishes, validates, and independently confirms the Temporal Equivalence Principle predictions through progressively rigorous empirical tests and cosmological extension. The convergence of evidence across five independent studies—theoretical prediction (Paper 0), multi-center validation (Paper 1), 25-year longitudinal analysis (Paper 2), raw data confirmation (Paper 3), and cosmological implications (Paper 4)—establishes a robust empirical foundation for TEP and extends its implications to dark matter phenomenology.

    The theoretical framework (Paper 0) made specific quantitative predictions before any empirical analysis: correlation length λ = 1,000–10,000 km, exponential spatial decay, velocity-dependent anisotropy, and absence of GM/r² scaling. Three subsequent empirical investigations confirmed these predictions through complementary methodologies, each addressing distinct artifact hypotheses. Paper 4 then extended the framework to cosmological scales, demonstrating how the same scalar field producing ~4,000 km correlations on Earth could manifest as apparent dark matter in gravitational lensing through violation of the Isochrony Axiom and Vainshtein screening (rV ∝ M1/3).

    Paper 1 (Multi-Center Validation) analyzed 62.7 million station-pair measurements across CODE, IGS, and ESA analysis centers, achieving R² = 0.92–0.97 consistency and demonstrating the correlation structure is not processing-specific. Paper 2 (25-Year Longitudinal) extended the temporal baseline to 25.3 years with 165.2 million station pairs, enabling detection of long-period signatures including 18.6-year lunar nutation coupling (R² = 0.641) and CMB frame alignment (18.2° separation, 5,570× variance ratio). Paper 3 (Raw RINEX Validation) analyzed 1.17 billion pair-samples using Single Point Positioning with broadcast ephemerides, achieving 100% detection rate across 72 metric combinations and confirming the signal exists in fundamental observables.

    Seven statistically independent signatures emerge with joint probability p ≈ 2×10−27 (>10σ): exponential spatial decay (λ = 4,201 ± 1,967 km), spatial anisotropy (EW/NS = 2.16), orbital velocity coupling (r = −0.888, 5.1σ), CMB frame alignment, planetary event responses (56/156 significant), 18.6-year nutation coupling, and semiannual nutation coupling (R² = 0.904). The signal remains invariant to solar cycles, geomagnetic storms, and altitude, exhibits no GM/r² scaling, and persists across independent processing pipelines and raw data—systematically excluding conventional artifact explanations.

## 7.2 Theoretical Interpretation

    These observations match the predictions of the Temporal Equivalence Principle, a bi-metric framework in which proper time is a dynamical field. TEP made specific, quantitative forecasts—correlation length, exponential decay, velocity-dependent anisotropy, absence of GM/r² (force) and GM/r³ (tidal) scaling—all of which are confirmed by the data. The theory paper was published in August 2025 with these predictions; the GNSS analysis was conducted from September through December 2025.

    The framework preserves local Lorentz invariance, with the speed of light remaining exactly invariant in any freely falling laboratory, while predicting global path-dependent synchronization. GW170817 constrains only the disformal sector; the conformal sector responsible for clock-rate modulation remains unconstrained by multi-messenger observations.

## 7.3 Cosmological Implications

    If validated, TEP implies that dark matter phenomenology in gravitational lensing arises from temporal-field gradients rather than particulate matter—the projection of differential proper-time accumulation onto observations that assume synchrony. The Earth-Galaxy scaling argument connects the ~4,000 km GNSS correlation to the ~50 kpc dark matter halo radius through Vainshtein screening (rV ∝ M1/3), suggesting that these represent the same field at different density scales.

    The framework naturally produces the MOND acceleration scale (a0 ~ cH0) and offers potential resolution to the H0 and S8 tensions. These are predictions of the framework rather than post-hoc explanations.

## 7.4 The Path Forward

    Independent replication is essential. The critical next steps are: raw carrier-phase analysis by independent research groups using the publicly available benchmark dataset; multi-constellation testing across GLONASS, Galileo, and BeiDou to verify universal coupling; lensed FRB timing to probe achromatic residuals at millisecond precision; and triangle holonomy tests to directly measure synchronization non-closure. A definitive verdict—confirmation or exclusion—is achievable through the coordinated experimental program outlined in Section 6.

## 7.5 Concluding Remarks

    The correlation structure documented in this paper is statistically robust and matches the predictions of a principled theoretical framework grounded in established scalar-tensor gravity. Whether it ultimately reflects new physics or reveals subtle systematics yet to be identified, the signal demands explanation.

    For three decades, standard GNSS processing systematically removed a ~4,000 km correlation as part of routine noise filtering. The algorithms were designed to optimize positioning accuracy rather than to characterize fundamental physics. The correlation structure reported here was present in the residuals but categorized as systematic error rather than physical phenomenon.

    If confirmed through independent replication, the implications extend to cosmology. The Isochrony Axiom—the assumption that all photons in a lensed image represent a synchronous snapshot—has never been empirically tested. If this axiom fails in the presence of conformal metric couplings, apparent dark matter could arise as a geometric artifact: temporal depth projected onto the spatial plane.

    Einstein established that simultaneity is relative. The TEP framework suggests that simultaneity may be generally non-integrable—that the rate of proper time accumulation is itself a dynamical field. The speed of light remains exactly invariant locally, but global synchronization becomes path-dependent.

    The ~4,000 km correlation on Earth (λ = 4,201 ± 1,967 km) and the ~50 kpc dark matter halo in galaxies may represent the same scalar field at different density scales, connected by Vainshtein screening. This hypothesis is falsifiable through the experimental program outlined in Section 6. Independent replication is essential, and the data and analysis code are publicly available. The community is invited to verify, challenge, or refute these findings.

## Data and Code Availability

    All analysis code, processed data products, and raw RINEX processing pipelines are publicly available:

    - Multi-center and 25-year CODE analysis: [github.com/matthewsmawfield/TEP-GNSS](https://github.com/matthewsmawfield/TEP-GNSS)

    - Raw RINEX validation: [github.com/matthewsmawfield/TEP-GNSS-RINEX](https://github.com/matthewsmawfield/TEP-GNSS-RINEX)

    - Benchmark dataset: A simplified CSV file containing binned residuals versus distance for the 25-year period, available for theorists who wish to fit models without processing raw GNSS data: [TEP-GNSS/benchmark](https://github.com/matthewsmawfield/TEP-GNSS/tree/main/benchmark)

## Acknowledgments

    The author expresses gratitude to the International GNSS Service (IGS), the Center for Orbit Determination in Europe (CODE), the European Space Agency (ESA), and NASA's Crustal Dynamics Data Information System (CDDIS). This analysis was possible because these organizations have championed the principles of open science, making decades of high-precision timing data freely available to the global research community. The GFZ Helmholtz Centre Potsdam is acknowledged for providing geomagnetic Kp indices.

## References

        Abbott, B. P. et al. (LIGO/Virgo) (2017). GW170817: Observation of gravitational waves from a binary neutron star inspiral. *Phys. Rev. Lett.* 119, 161101.

        Ashby, N. (2003). Relativity in the Global Positioning System. *Living Rev. Relativity* 6, 1.

        Bekenstein, J. D. (1993). The relation between physical and gravitational geometry. *Phys. Rev. D* 48, 3641.

        Bothwell, T. et al. (2022). JILA Sr optical lattice clock with 10⁻¹⁸ stability and accuracy. *Nature* 602, 420–424.

        Boulder Atomic Clock Optical Network (BACON) Collaboration (2021). Frequency ratio measurements at 18-digit accuracy using an optical clock network. *Nature* 591, 564–569.

        Burrage, C. & Sakstein, J. (2018). Tests of chameleon gravity. *Living Rev. Relativity* 21, 1.

        Damour, T. & Polyakov, A. M. (1994). The string dilaton and a least coupling principle. *Nucl. Phys. B* 423, 532.

        Hofmann-Wellenhof, B., Lichtenegger, H. & Wasle, E. (2008). *GNSS – Global Navigation Satellite Systems*. Springer.

        Khoury, J. & Weltman, A. (2004). Chameleon fields: Awaiting surprises for tests of gravity in space. *Phys. Rev. Lett.* 93, 171104.

        Milgrom, M. (1983). A modification of the Newtonian dynamics as a possible alternative to the hidden mass hypothesis. *ApJ* 270, 365.

        Petit, G. & Luzum, B. (2010). IERS Conventions (2010). *IERS Technical Note* 36.

        Planck Collaboration (2020). Planck 2018 results. VI. Cosmological parameters. *A&A* 641, A6.

        Riess, A. G. et al. (2022). A comprehensive measurement of the local value of the Hubble constant. *ApJ* 934, L7.

        Smawfield, M. L. (2025a). The Temporal Equivalence Principle: Dynamic Time, Emergent Light Speed, and a Two-Metric Geometry of Measurement. Zenodo. DOI: [10.5281/zenodo.16921911](https://doi.org/10.5281/zenodo.16921911)

        Smawfield, M. L. (2025b). Global Time Echoes: Distance-Structured Correlations in GNSS Clocks. Zenodo. DOI: [10.5281/zenodo.17127229](https://doi.org/10.5281/zenodo.17127229)

        Smawfield, M. L. (2025c). Global Time Echoes: 25-Year Temporal Evolution of Distance-Structured Correlations in GNSS Clocks. Zenodo. DOI: [10.5281/zenodo.17517141](https://doi.org/10.5281/zenodo.17517141)

        Smawfield, M. L. (2025d). Global Time Echoes: Raw RINEX Validation of Distance-Structured Correlations in GNSS Clocks. Zenodo. DOI: [10.5281/zenodo.17860166](https://doi.org/10.5281/zenodo.17860166)

        Smawfield, M. L. (2025e). Global Time Echoes: Cosmological Implications of Distance-Structured Correlations in GNSS Clocks. Zenodo. DOI: [10.5281/zenodo.17982540](https://doi.org/10.5281/zenodo.17982540)

        Teunissen, P. J. G. & Montenbruck, O. (Eds.) (2017). *Springer Handbook of Global Navigation Satellite Systems*. Springer.

        Touboul, P. et al. (2022). MICROSCOPE mission: test of the equivalence principle in space. *Phys. Rev. Lett.* 129, 121102.

        Will, C. M. (2014). The Confrontation between General Relativity and Experiment. *Living Rev. Relativity* 17, 4.

## Appendix A: Mathematical Formalism

## A.1 Full Action

    The TEP action is formulated in the Einstein frame where gravity is canonical:

    $S = \int d^4x \sqrt{-g} \left[ \frac{M_{\text{Pl}}^2}{2} R - \frac{1}{2} K(\phi)\, g^{\mu\nu} \partial_\mu \phi\, \partial_\nu \phi - V(\phi) \right] + S_{\text{matter}}[\psi, \tilde{g}_{\mu\nu}, \phi]$

    where MPl is the reduced Planck mass, R is the Ricci scalar, K(φ) is the kinetic function, V(φ) is the screening potential, and Smatter is the matter action in the Jordan frame with coupling to the matter metric g̃μν.

## A.2 Field Equations

    Variation with respect to the gravitational metric yields the Einstein equations:

    $G_{\mu\nu} = \frac{1}{M_{\text{Pl}}^2} \left[ T_{\mu\nu}^{(\phi)} + T_{\mu\nu}^{(\text{matt})} \right]$

    where the scalar field stress-energy tensor is:

    $T_{\mu\nu}^{(\phi)} = K(\phi) \partial_\mu \phi \partial_\nu \phi - g_{\mu\nu} \left[ \frac{1}{2} K(\phi) g^{\rho\sigma} \partial_\rho \phi \partial_\sigma \phi + V(\phi) \right]$

    Variation with respect to the scalar field yields:

    $\nabla_\mu \left[ K(\phi) \nabla^\mu \phi \right] - \frac{\partial K}{\partial \phi} \frac{1}{2} g^{\mu\nu} \partial_\mu \phi \partial_\nu \phi - \frac{\partial V}{\partial \phi} = -\alpha(\phi) T + S_{\text{disf}}$

    where α(φ) ≡ d ln A/dφ = 2β/MPl is the conformal coupling strength, T is the trace of the matter stress-energy tensor in the Jordan frame, and Sdisf contains disformal source terms proportional to B(φ).

## A.3 Galileon Terms and Vainshtein Screening

    The cubic Galileon term in the action generates nonlinear screening:

    $S_{\text{Galileon}} = \int d^4x \sqrt{-g} \left[\frac{c_3}{\Lambda^3} (\partial\phi)^2 \Box\phi\right]$

    where Λ is the strong-coupling scale and c3 is a dimensionless coefficient of order unity. This term becomes important when:

    $\frac{|\nabla \phi|^2}{\Lambda^3} \sim 1$

    The Vainshtein radius, inside which nonlinear effects dominate, is:

    $r_V = \left(\frac{GM}{\Lambda^3}\right)^{1/3}$

    For a spherically symmetric source of mass M, the scalar field profile in the Vainshtein regime (r < rV) is:

    $\phi(r) \sim \frac{\beta M_{\text{Pl}}}{M_{\text{Pl}}^2} \frac{GM}{r} \left(\frac{r}{r_V}\right)^3$

    This suppresses the fifth force by a factor of (r/rV)3 inside the screening radius, reconciling the weak-field cosmological behavior with strong-field solar system constraints.

## A.4 Conservation Laws

    The total stress-energy tensor is conserved:

    $\nabla_\mu \left[ T^{\mu\nu}_{(\phi)} + T^{\mu\nu}_{(\text{matt})} \right] = 0$

    However, the matter stress-energy alone is not conserved in the Einstein frame:

    $\nabla_\mu T^{\mu\nu}_{(\text{matt})} = \alpha(\phi) T \nabla^\nu \phi + \text{disformal terms}$

    This apparent violation is resolved by transforming to the Jordan frame, where matter follows geodesics of g̃μν and energy-momentum is conserved with respect to the matter-frame covariant derivative.

## Appendix B: Processing Filter Mathematical Formulation

## B.1 Network Adjustment Constraint

    Standard GNSS processing applies a datum constraint that forces residuals to have zero spatial mean:

    $\delta t_i^{\text{adj}} = \delta t_i - \frac{1}{N}\sum_{j=1}^N \delta t_j$

    This constraint preserves zero-mean residuals while removing any common-mode offset.

## B.2 Fourier Decomposition

    Decompose the TEP signal into spatial Fourier modes characterized by wavevector k:

    $\delta t_i^{\text{TEP}} = \sum_{\mathbf{k}} A_{\mathbf{k}} e^{i\mathbf{k} \cdot \mathbf{r}_i}$

    where ri is the position of station i. The network adjustment removes the k = 0 (monopole) component:

    $A_{\mathbf{k}=0}^{\text{adj}} = A_{\mathbf{k}=0} - \frac{1}{N}\sum_j A_{\mathbf{k}=0} = 0$

    but preserves all k ≠ 0 (dipole, quadrupole, and higher) components that have zero spatial mean.

## B.3 Cross-Spectral Density

    The TEP correlation structure manifests as phase coherence between station pairs, quantified by the cross-spectral density:

    $C_{ij}(\omega) = \langle \tilde{\delta t}_i(\omega) \tilde{\delta t}_j^*(\omega) \rangle$

    where tildes denote Fourier transforms and angle brackets denote ensemble averaging.

## B.4 Phase Alignment Index

    The magnitude-weighted phase alignment index extracts the phase relationship:

    $\Phi_{ij} = \frac{\sum_\omega |C_{ij}(\omega)| \cos[\arg C_{ij}(\omega)]}{\sum_\omega |C_{ij}(\omega)|}$

    This metric is invariant under common-mode amplitude rescaling because it depends only on relative phases, not absolute amplitudes. The datum constraint removes Ak=0 (amplitude information) but preserves the phase structure encoded in the k ≠ 0 components.

## B.5 Exponential Correlation Function

    The exponential decay C(r) ∝ exp(−r/λ) reflects the two-point correlation function of the underlying φ field:

    $\langle \phi(\mathbf{x}) \phi(\mathbf{x}') \rangle = \phi_0^2 \exp\left(-\frac{|\mathbf{x} - \mathbf{x}'|}{\lambda}\right)$

    where λ = ℏ/(mφc) is the Compton wavelength of the screened scalar field. This correlation function has nonzero spatial gradients (k ≠ 0 in Fourier space) and therefore survives network adjustment, which only removes the k = 0 component.

## B.6 Tidal vs. TEP Signal Discrimination

    Classical tidal effects are parameterized as sums of spherical harmonics Ylm(θ, φ) with known frequencies:

    $\delta t_{\text{tidal}}(t, \mathbf{r}) = \sum_{l,m,n} A_{lmn} Y_{lm}(\theta, \phi) \cos(\omega_n t + \phi_n)$

    where ωn are tidal frequencies (diurnal, semidiurnal, long-period). Standard GNSS processing explicitly models and subtracts these terms using IERS conventions.

    The TEP correlation is a distance-dependent phase coherence:

    $C_{\text{TEP}}(r_{ij}) = A_0 \exp(-r_{ij}/\lambda)$

    This functional form is not decomposable into the standard tidal spherical harmonic basis. The processing removes what it models (tidal harmonics) but preserves what it does not model (exponential distance correlations in phase space).

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

*This document was automatically generated from the TEP-GTE research site. For the interactive version with figures and enhanced formatting, visit: https://matthewsmawfield.github.io/TEP-GTE/*

*Related Work:*
- [**TEP Theory**](https://doi.org/10.5281/zenodo.16921911) (Foundational framework)

*Source code and data available at: https://github.com/matthewsmawfield/TEP-GTE*
