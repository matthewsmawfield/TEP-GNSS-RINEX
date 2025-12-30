# Global Time Echoes: Optical Validation of the Temporal Equivalence Principle via Satellite Laser Ranging

**Author:** Matthew Lukin Smawfield  
**Version:** v0.1 (Mombasa)  
**Date:** First published: 30 December 2025 · Last updated: 30 December 2025  
**DOI:** 10.5281/zenodo.18064582  
**Generated:** 2025-12-30  
**Paper Series:** TEP-SLR (Paper 9)

---

## Abstract

An independent, optical-domain test of the Temporal Equivalence Principle (TEP) is presented using 11 years (2015–2025) of Satellite Laser Ranging (SLR) data from passive ILRS geodetic satellites (LAGEOS-1/2 and Etalon-1/2). This analysis constrains "clock-artifact" explanations by employing two-way optical ranging to passive retroreflectors—a methodology orthogonal to the microwave measurements of active atomic clocks used in Global Navigation Satellite Systems (GNSS).

Under strict 5-minute contemporaneous binning, distance-binned mean pass-correlations fluctuate with high variance. However, widening the overlap window to 15 minutes (thereby increasing multi-station overlap) reveals statistically significant, distance-structured inter-station correlations (Fisher-combined $\chi^2=15.35$ with 4 d.o.f.; $p=0.0040$) under a family-wise circular-shift test.

This signal is driven primarily by LAGEOS-2 ($p=0.0005$), which exhibits a strong negative correlation ($r \approx -0.59$) in the 5,000–7,500 km distance bin, whereas LAGEOS-1 remains consistent with the null hypothesis ($p \approx 0.93$). Although observation counts and temporal overlap are comparable, this asymmetry likely reflects a combination of orbital geometry—LAGEOS-2's prograde $52.6^\circ$ orbit versus LAGEOS-1's retrograde $109.8^\circ$ orbit—and small-number statistics in the critical distance bin.

To validate this finding with more robust statistics, a daily-aggregation analysis ($N=190$ station pairs) was performed. This confirmed a subtler but statistically significant negative correlation at shorter ranges (500–1,000 km, $p=0.017$), suggesting a persistent global background structure independent of the high-amplitude LAGEOS-2 events.

The detection of matching low-frequency structure in a system devoid of active clocks and microwave propagation challenges receiver electronics, clock steering, and ionospheric modeling errors as complete explanations. While current network sparsity limits testing to the conformal sector, this work demonstrates SLR as an independent, technology-orthogonal line of evidence for TEP phenomenology.

## 1. Introduction

### 1.1 The Necessity of Independent Testing

The Temporal Equivalence Principle (TEP) posits that proper time is a dynamical field governed by a conformal factor $A(\phi) = \exp(2\beta\phi/M_{\text{Pl}})$, leading to path-dependent synchronization effects that manifest as spatial correlations in distributed timing networks. Previous analyses of the GNSS network (Smawfield 2025b, 2025c, 2025d; Papers 1-3) identified a persistent correlation structure with a length scale of $\lambda \approx 4,000$ km, providing evidence consistent with TEP's conformal-sector phenomenology. While robust across processing centers (CODE, IGS, ESA; $R^2 = 0.92-0.97$), temporally stable over 25 years, and present in raw RINEX observations, these findings relied on a single geodetic technique: one-way microwave transmission to active ground clocks.

This reliance leaves open a critical counter-hypothesis: that the signal arises from subtle couplings in receiver electronics (e.g., thermal responses), clock steering algorithms, ionospheric modeling errors specific to the microwave L-band, or pervasive "colored noise" (flicker noise) often treated as irreducible station systematics. To independently test this "clock artifact" hypothesis, the conformal sector predictions must be sought in a system with no active clocks and no microwave propagation. This study addresses this challenge by analyzing 11 years (2015–2025) of Satellite Laser Ranging (SLR) data from the LAGEOS-1 and LAGEOS-2 missions. Unlike GNSS, SLR relies on two-way optical pulses reflected off passive retroreflectors. The "clock" is a ground-based event timer, and the observable is the pure round-trip time of flight. If the TEP signal represents a fundamental metric propagation anomaly arising from the conformal sector, it is expected to manifest in these optical residuals—potentially offering a physical origin for the persistent "flicker noise" floor and scale drifts observed in geodetic time series.

### 1.2 TEP's Two Sectors: Conformal vs. Disformal

The TEP framework (Smawfield 2025a; Paper 0) distinguishes two physically distinct sectors with different observational signatures:

#### Conformal Sector (Clock-Rate Modulation)

**Coupling:** Universal conformal factor $A(\phi) = \exp(2\beta\phi/M_{\text{Pl}})$ modulates proper time rates: $\mathrm{d}\tau/\mathrm{d}t \propto A(\phi)^{1/2}$.

**Observable:** Spatial correlations in clock frequencies with characteristic length $\lambda \sim m_\phi^{-1}$ (screening length).

**Constraint Status:** *Unconstrained by GW170817* (multi-messenger bounds apply only to disformal sector).

**GNSS Evidence (Smawfield 2025b,c,d; Papers 1-3):** Supported by distance-structured correlations ($\lambda = 4,201 \pm 1,967$ km), orbital velocity coupling ($r = -0.888$, 5.1$\sigma$), and CMB frame alignment (18.2° from dipole).

#### Disformal Sector (Synchronization Holonomy)

**Coupling:** Disformal term $B(\phi) \nabla_\mu\phi \nabla_\nu\phi$ tilts photon light cones, creating path-dependent one-way time asymmetries.

**Observable:** Synchronization holonomy $H = \oint_C \Omega_\mu \mathrm{d}x^\mu \neq 0$ in closed-loop time transfer (after GR subtraction).

**Constraint Status:** *Tightly bounded by GW170817:* $|c_\gamma - c_g|/c \lesssim 10^{-15}$ requires $B(\phi)(\partial\phi)^2$ to be negligible today. This constraint applies specifically to the disformal sector because the disformal coupling tilts photon light cones independently of graviton propagation, creating a differential speed. In contrast, conformal coupling modulates proper time rates universally, affecting both photons and gravitons equally, leaving their relative speed unchanged and thus evading the GW170817 bound.

**GNSS Evidence:** *Untested* (no closed-loop holonomy experiments performed).

### 1.3 SLR Predictions: Conformal Sector Only

**Scope of This Analysis:** This paper tests TEP's conformal sector predictions only. The sparse ILRS network (47 stations) and two-way measurement geometry do not provide the closed-loop time-transfer topology and dense temporal coverage required to test disformal predictions (synchronization holonomy), including the proposed "Anti-Echo" sign-inversion signature.

The conformal sector makes three testable predictions for SLR:

- **Propagation Dependence:** Residual temporal structure should vary systematically with signal path length (stronger at lower elevations), distinguishing conformal propagation effects from purely local station systematics. *Prediction: order-unity to order-of-magnitude path dependence in a lag-1 diagnostic. Observed: a strong short-path vs long-path contrast in a gap-aware, binned lag-1 statistic, with estimator-sensitive uncertainty (§3.2).*

- **Spectral Structure:** Residuals should possess a non-white spectral signature concentrated in the low-frequency TEP band ($10-500 \mu\text{Hz}$), consistent with the transit time of Earth through the $\phi$ field. *Prediction: 2–3× enhancement relative to a full-spectrum mean PSD, with a larger contrast relative to a high-frequency broadband floor. Observed: 2.48× vs full-spectrum mean, and 14.00× vs broadband floor (§3.3).*

- **Frequency Independence:** The conformal coupling is achromatic; optical (SLR, $\sim 500$ THz) and microwave (GNSS, $\sim 1$ GHz) should show consistent correlation lengths when scaled by their respective screening environments (accounting for local matter density differences). *Prediction: Range-dependent coherence structure consistent with GNSS $\lambda \approx 4,000$ km. Observed: consistent scaling (§5.1).*

- **Spatial Coherence:** Residuals should exhibit distance-dependent spatial correlations consistent with the screening length $\lambda$. *Prediction: Negative correlations at intermediate ranges (Anti-Echo mechanism) or short-range coherence. Observed: Significant signals in daily aggregation ($p=0.017$) and high-cadence pass correlations ($p=0.004$) (§3.4).*

**What This Paper Does NOT Test:** The "Anti-Echo" sign inversion is a disformal prediction requiring experimental configurations described in Smawfield (2025a; Paper 0, §10.1): closed-loop optical time transfer, interplanetary one-way asymmetry measurements, and triangle synchronization holonomy tests. These experiments remain to be performed.

## 2. Methodology

### 2.1 The Target: LAGEOS

The Laser Geodynamics Satellites (LAGEOS-1 and LAGEOS-2) represent excellent test masses for this investigation. As dense, passive spheres covered in retroreflectors, they have the highest mass-to-area ratio of any satellite, minimizing non-gravitational perturbations (e.g., drag, radiation pressure). Crucially, they carry no active electronics and no clocks. The "time" measurement is performed entirely by the ground segment's event timer, measuring the round-trip flight time of a photon. This effectively isolates the propagation metric from onboard clock systematics.

### 2.2 Dataset & Processing

The complete International Laser Ranging Service (ILRS) dataset for passive geodetic satellites (LAGEOS-1/2 and Etalon-1/2) was analyzed over an 11-year period (Mar 2015 – Dec 2025), comprising 6,571,782 Normal Point observations from 47 global stations. Residuals were computed for 3,704,757 observations with available ephemerides and modeling inputs.

Residuals were computed relative to high-precision SP3 orbits (ASI/GFZ). For the 2025 reporting period, care was taken to utilize a consistent single-center orbit solution (ASI) to avoid systematic noise introduced by mixed-center product aggregation. The reduction strategy proceeded in two stages:

- **Geodetic Validation:** A broad 5-meter outlier rejection window retained $\approx 1.7$ million "valid" geodetic observations. The initial RMS ($\approx 2.77$ m) reflects the raw pre-fit state relative to the *a priori* orbit, preserving large-scale signal structures that are typically removed by aggressive orbital fitting.

- **Coherence Analysis Subset:** To isolate subtle timing correlations from gross interpolation and modeling errors, a strict 0.5-meter (50 cm) threshold was applied. This high-precision subset (192,561 residuals, $\approx 5.2\%$ of ephemeris-resolved residuals) forms the basis of the primary inter-station and propagation diagnostics reported in this work. This cut prioritizes epochs where orbit interpolation and environmental corrections remain within the sub-meter regime. Robustness is quantified by an explicit residual-threshold sweep (0.3, 0.5, 1.0 m) in the analysis outputs.

- **Corrections:** Standard Marini-Murray troposphere model, Shapiro delay, and Sagnac corrections were applied. No station-specific meteorological data was used, to avoid introducing local sensor systematics.

- **Parameter Estimation Strategy:** Crucially, this analysis avoids the standard practice of estimating frequent empirical accelerations or "geographically correlated" parameters, which are often used in precise orbit determination (POD) to whiten residuals. By fixing the orbit to the high-precision ASI solution and avoiding secondary empirical filtering, the "common mode" signals—typically discarded as noise—are preserved for analysis.

### 2.3 Inter-Station Metric: Contemporaneous Pass-Bin & Daily Aggregation

For sparse SLR networks, continuous, regularly sampled inter-station time series are generally unavailable. Therefore, the primary inter-station metric used here was based on *contemporaneous pass bins*. For each satellite and each time bin (5-minute and 15-minute windows), a pass-mean residual anomaly was computed at each station after subtracting that station’s global mean residual. Inter-station correlation was then computed for each station pair by correlating these pass-mean anomalies across bins.

Statistical significance was assessed using a family-wise circular-shift permutation test across distance bins (2000 permutations). Additionally, a daily-aggregation analysis was performed where residuals were averaged daily per station to maximize temporal overlap ($N=190$ station pairs), providing a robust check against short-term pass-geometry artifacts.

As a secondary check, an irregular-sampling phase-alignment statistic was computed on the same contemporaneous pass-bin series (without interpolation) and evaluated under an analogous family-wise circular-shift null test.

### 2.4 Spectral Diagnostic (TEP Band)

In addition to the pass-bin spatial test, the residuals were examined in the frequency domain to quantify the concentration of power in the predicted TEP band (10–500 $\mu$Hz). For SLR, this spectral diagnostic was treated as a station-level characterization of low-frequency structure, complementary to the contemporaneous pass-bin inter-station statistic.

## 3. Conformal Sector Assessment

The analysis of the 11-year LAGEOS dataset reveals four distinct signatures consistent with TEP's conformal sector phenomenology in the optical domain, challenging clock-artifact explanations.

### 3.1 Data Characterization: The Noise Floor

The raw SLR residuals exhibit a broad distribution. Within the strictly filtered subset ($N=192,561$, $|\Delta\rho|<0.5$ m), the RMS is 0.29 m (288 mm). The central question is whether there exists reproducible structure in this noise floor consistent with conformal propagation effects. In standard geodetic analysis, this floor is often characterized as "flicker noise" ($1/f$) or "colored noise" (Williams et al., 2004) and is typically attributed to unmodeled station systematics. TEP predicts this colored spectrum as a physical consequence of scalar field coupling.

### 3.2 Range-Dependent Coherence: The Propagation Signature

To distinguish propagation effects from local station systematics, a lag-1 residual statistic was analyzed as a function of signal path length (range to satellite). Here "lag-1" denotes the correlation between successive 5-minute binned residual means, computed on station-debiased residuals and made gap-aware by excluding time gaps exceeding 450 s. This diagnostic targets sub-hour temporal structure while avoiding trivial within-pass autocorrelation.

- **Shorter path ($\lesssim 6{,}500$ km):** mean lag-1 $\approx -0.042$ (95% CI: −0.120 to 0.028).

- **Longer path ($\gtrsim 8{,}000$ km):** mean lag-1 $\approx -0.275$ (95% CI: −0.534 to −0.018).

![Residual Coherence vs Path Length](results/figures/slr_residual_vs_elevation_full.png)

**Figure 3.1:** Lag-1 residual statistic as a function of signal path length (range to satellite), computed using 5-minute binned residual means with gap-aware pairing. Under a strict $|\Delta\rho|<0.5$ m filter, the long-minus-short contrast is $\Delta(\mathrm{low}-\mathrm{high})=-0.233$ (95% bootstrap CI: −0.466 to −0.010). The corresponding low/high ratio is 6.58 (95% bootstrap CI: −50.95 to 51.33), but is ill-conditioned because the short-path mean is near zero. Under a looser $|\Delta\rho|<1.0$ m filter, the contrast becomes $\Delta(\mathrm{low}-\mathrm{high})=-0.044$ (95% CI: −0.208 to 0.100) and the corresponding ratio becomes −1.31 (95% CI: −47.58 to 25.51). This threshold sensitivity motivates conservative interpretation of the path-length diagnostic as qualitative support rather than a precisely estimated amplitude.

The observed path dependence provides qualitative support for a propagation-linked component. While the magnitude of the contrast is sensitive to the outlier rejection threshold (see Figure 3.1), the sign structure remains consistent: longer paths accumulate greater decoherence. This sensitivity reflects estimator dependence and the near-zero short-path baseline; the path-length diagnostic is therefore treated conservatively as qualitative support for propagation effects, pending further robustness testing. The persistence of strong spectral concentration in the predicted TEP band (Section 3.3) provides the most significant quantitative evidence for a low-frequency, non-white process—a 14.00× enhancement that is threshold-stable and statistically robust.

### 3.3 Spectral Concentration: Primary Quantitative Signature

Frequency-domain analysis of the residuals reveals a significant concentration of power within the predicted TEP frequency band (10–500 $\mu$Hz). On 5-minute resampled station series, the station-averaged TEP-band mean PSD exceeds the full-spectrum mean PSD by $2.48\times$ (95% CI: 2.46–2.50; $N=46$ stations). Relative to a broadband floor defined by $f > 1$ mHz, the TEP-band mean exceeds broadband by $14.00\times$ (95% CI: 13.53–14.47; $N=46$ stations). This "spectral clumping" indicates that the signal is not white noise but a structured, low-frequency process consistent with the transit time of Earth through a scalar domain structure—matching the spectral characteristics observed in GNSS (Smawfield 2025b, 2025c, 2025d; Papers 1-3). The TEP band (10–500 $\mu$Hz) corresponds to periods of ~30 minutes to ~28 hours, consistent with Earth's motion through large-scale scalar field gradients.

### 3.4 Inter-Station Analysis: Expected Network Limitations

The ILRS network is sparse in time; most station pairs lack sufficient contiguous overlap to support the cross-spectral methods employed in GNSS analysis. Inter-station inference uses contemporaneous, same-satellite 5-minute bins: for each satellite and time bin, pass-mean residual anomalies are computed per station (after removing each station's global mean), and station-pair correlations are estimated across bins.

**Key Finding:**

#### Inter-Station Pass Correlations: Statistical Evidence with Geometric Complexity

Distance-binned mean pass-correlation estimates fluctuate with large variance under a strict 5-minute contemporaneous binning. However, when the contemporaneous overlap window is widened to 15-minute bins (increasing the number of multi-station overlap epochs), a family-wise circular-shift test yields statistically significant evidence for systematic distance-structured inter-station correlations (Fisher-combined $\chi^2=15.35$ with 4 d.o.f.; $p=0.0040$).

This signal is driven almost entirely by LAGEOS-2 ($p=0.0005$), which exhibits a strong negative correlation ($r \approx -0.59$) in the 5,000–7,500 km distance bin. In contrast, LAGEOS-1 remains consistent with the null hypothesis ($p \approx 0.93$). While observation counts and temporal overlap are nearly identical for both satellites, this asymmetry likely reflects a combination of orbital geometry—LAGEOS-2's prograde $52.6^\circ$ orbit versus LAGEOS-1's retrograde $109.8^\circ$ orbit—and small-number statistics in the critical distance bin ($N=3$ pairs for LAGEOS-2 vs $N=8$ for LAGEOS-1).

To examine this finding with improved statistics, a daily-aggregation analysis ($N=190$ station pairs) was performed. This identified a subtler but statistically significant negative correlation at shorter ranges (500–1,000 km, $p=0.017$), suggesting a persistent global background structure independent of the high-amplitude LAGEOS-2 events.

![SLR Pass-Correlation vs Distance](results/figures/slr_pass_correlation_decay.png)

**Figure 3.2:** Pass-based inter-station correlation of SLR residual anomalies as a function of baseline distance. While the primary 5-minute analysis shows limited structure, widening the window to 15 minutes reveals a statistically significant signal (Fisher combined $p=0.0040$), driven largely by strong anti-correlations in LAGEOS-2 residuals at 5,000–7,500 km baselines. A complementary daily-aggregation analysis ($N=190$ pairs) identifies significant negative correlations at 500–1,000 km ($p=0.017$).

The range-dependent and spectral signatures (Sections 3.2-3.3) offer a technology-orthogonal line of evidence consistent with the conformal-sector phenomenology reported in GNSS Papers 1-3. More stringent inter-station tests will require denser temporal overlap and experimental configurations beyond current ILRS cadence.

## 4. Future Experimental Directions

Based on the conformal-sector evidence presented, future experimental directions for testing the disformal sector are outlined below. Specifically, the "Anti-Echo" mechanism and the experimental requirements for future closed-loop time-transfer tests are described.

### 4.1 Conformal vs. Disformal: Two Distinct Predictions

The TEP bi-metric geometry $\tilde{g}_{\mu\nu} = A(\phi) g_{\mu\nu} + B(\phi) \nabla_\mu\phi \nabla_\nu\phi$ contains two physically distinct coupling mechanisms:

- **Conformal Coupling $A(\phi)$:** Modulates clock rates universally. Creates spatial correlations in timing residuals with length scale $\lambda \sim m_\phi^{-1}$. *This is the sector probed in GNSS and SLR analyses.*

- **Disformal Coupling $B(\phi)$:** Tilts photon light cones in directions transverse to $\nabla\phi$. Creates one-way time asymmetries and synchronization holonomy $H \neq 0$ in closed loops. *This requires closed-loop time transfer to test.*

The multi-messenger constraint from GW170817 requires $|c_\gamma - c_g|/c \lesssim 10^{-15}$, forcing $B(\phi)(\partial\phi)^2 \approx 0$ today. This bounds the disformal sector while leaving the conformal sector unconstrained. GNSS and SLR provide evidence consistent with conformal predictions; disformal predictions await dedicated holonomy experiments.

### 4.2 Estimator-Dependent Sign Structure: The "Anti-Echo" Mechanism

Within the disformal sector, TEP predicts *estimator dependence*: the same underlying common-mode propagation delay can map differently into post-fit residuals depending on whether the estimator absorbs global modes into clocks (kinematic GNSS positioning) or into orbital parameters (dynamic SLR orbit determination).

The mechanism operates as follows:

**Schematic:**

- GNSS (Kinematic PPP):  τTEP → δtreceiver → r > 0

- SLR (Dynamic OD):      τTEP → δaorbit   → r < 0 (Predicted)

- **GNSS (Kinematic PPP):** In Precise Point Positioning, the receiver coordinates and clock bias are solved epoch-by-epoch. The satellite orbit is fixed (from IGS products), but the receiver state is free. A common-mode TEP delay ($\bar{\tau}$) affecting a region is simply mapped into the receiver clock bias estimate. Since both stations measure this common delay, their residuals remain positively correlated ($r > 0$).

- **SLR (Dynamic Orbit Determination):** In SLR, a single orbital arc is fitted to observations from stations worldwide over several days. A persistent TEP delay acts phenomenologically like a scale error or an unmodelled drag force. The least-squares filter minimizes the global residual by adjusting the orbital parameters—typically the semi-major axis. This effectively absorbs the monopole term ($\bar{\tau}$) into the orbit solution, potentially leaving residuals that reflect deviations from the absorbed global average—manifesting as anti-correlation ($r < 0$) at regional scales.

### 4.3 Experimental Design Considerations

Monte Carlo simulations (see supplementary materials) demonstrate that estimator-dependent sign structure is plausible: when dynamic orbit fits absorb common-mode delays into orbital parameters, post-fit residual correlations can exhibit sign inversion at regional baselines. This provides the theoretical basis for future experimental designs.

### 4.4 Experimental Requirements for Testing the Anti-Echo

Detecting the "Anti-Echo" in real SLR data requires:

- **Dense Temporal Overlap:** Synchronous observations from multiple stations within correlation timescales (~hours), not available in current ILRS cadence.

- **Cross-Estimator Comparison:** Parallel processing with kinematic (epoch-by-epoch) and dynamic (arc-fit) estimators to isolate the sign-flip signature.

- **Alternative Orbit Solutions:** Sensitivity analysis across different Analysis Center products (ASI, GFZ, CSR) to test monopole absorption consistency.

The current SLR dataset does not meet these requirements. Future targeted SLR campaigns or dedicated optical time-transfer networks with denser coverage will be needed to test disformal predictions. The conformal-sector evidence presented in Section 3 represents the primary scientific contribution of this work.

## 5. Synthesis & Critical Analysis

### 5.1 The Ladder of Evidence: Conformal Sector Assessment

The TEP-SLR analysis provides an independent optical-domain test of the Temporal Equivalence Principle's conformal sector. By integrating these findings with the previous GNSS results (Smawfield 2025b, 2025c, 2025d; Papers 1-3), a "Ladder of Evidence" is proposed that supports the conformal interpretation while clearly distinguishing it from untested disformal predictions. Throughout, robust empirical observables are distinguished from interpretive inferences that remain contingent on additional systematics controls:

**Critical Analysis:**

#### 1. Universality (Frequency Independence)

The conformal signal affects optical frequencies ($\sim 10^{14}$ Hz, SLR) as well as microwave ($\sim 10^9$ Hz, GNSS). This vast frequency difference (factor of $\sim 10^5$) provides a strong argument against dispersive propagation effects, such as ionospheric delay or plasma dispersion, which scale with frequency ($1/f^2$). Furthermore, the spectral power concentration (14.00× enhancement in the TEP band relative to broadband, 95% CI: 13.53–14.47) is consistent with an achromatic, temporally structured anomaly. Taken together, the optical and microwave results support the universality of the conformal coupling $A(\phi)$ across the electromagnetic spectrum.

#### 2. Propagation Origin (Metric Modulation)

The residuals exhibit a path-length dependent temporal contrast in a gap-aware lag-1 diagnostic. In the primary $|\Delta\rho|<0.5$ m analysis, the long-minus-short contrast is $\Delta(\mathrm{low}-\mathrm{high})=-0.233$ (95% bootstrap CI: −0.466 to −0.010), while the short-path mean is near zero and therefore ratio-based summaries are ill-conditioned. This contrast is threshold-sensitive (e.g., $\Delta(\mathrm{low}-\mathrm{high})=-0.044$, 95% CI: −0.208 to 0.100 at $|\Delta\rho|<1.0$ m), and is therefore treated as qualitative support for a path-accumulated component rather than a precisely estimated amplitude. Crucially, standard tropospheric corrections (Marini-Murray) were applied prior to this analysis; the structure persists in the *post-correction* residuals, motivating deeper atmospheric-systematics testing rather than serving as a refutation of the baseline refraction model.

#### 3. Scale Consistency (The Screening Radius)

The GNSS correlation length ($\lambda = 4,201 \pm 1,967$ km) and the SLR path-length dependent contrast (over a $\sim 3{,}000$ km path difference between the high- and low-range selections) are both consistent with Vainshtein screening at the dark energy scale. The screening length is set by the scalar field mass: $\lambda = \hbar/(m_\phi c) \approx 4{,}000$ km for $m_\phi \approx 5 \times 10^{-14}\,\mathrm{eV}/c^2$. In natural-units order-of-magnitude terms, this mass scale is also in the range often associated with vacuum-energy scales (heuristically $\rho_\Lambda^{1/4} \sim m_\phi$), motivating the cosmological interpretation without requiring a direct identification. While the sparse ILRS network limits direct measurement of a continuous inter-station correlation length, the range-dependent diagnostic and TEP-band spectral concentration provide complementary constraints. The convergence of GNSS and SLR evidence at similar spatial scales supports the screening interpretation across two independent measurement systems.

#### 4. Inter-Station Correlation (Spatial Coherence)

While the sparse network limits synchronous overlap, advanced correlation techniques have successfully isolated a distance-structured signal. High-cadence (15-minute) pass-bin analysis reveals a statistically significant signal ($p=0.004$) driven by strong anti-correlations in LAGEOS-2 data at intermediate ranges (5,000–7,500 km). Furthermore, a robust daily-aggregation analysis ($N=190$ station pairs) identifies widespread negative correlations at shorter ranges (500–1,000 km, $p=0.017$). This spatial coherence—observed in an optical, passive system—mirrors the distance-structured correlations found in microwave GNSS networks, offering a technology-orthogonal test of the conformal sector phenomenology.

#### 5. Disformal Sector: Untested in Current Data

The "Anti-Echo" sign inversion is a *disformal prediction* ($B(\phi) \neq 0$) requiring closed-loop time transfer to test. Neither GNSS (Papers 1-3) nor SLR (this work) has performed such experiments with the topology and cadence required for synchronization-holonomy inference. The disformal sector remains a well-defined, falsifiable prediction for future experiments: closed-loop optical time transfer (Smawfield 2025a; Paper 0, §10.1), interplanetary one-way asymmetry measurements, and triangle synchronization holonomy tests. The simulation (Figure 4.1) demonstrates the mechanism's plausibility and provides quantitative forecasts for these future tests.

### 5.2 Limitations & Systematics

While the evidence supports a conformal-sector interpretation, several systematic effects warrant specific attention to assess the robustness of these findings.

#### Atmospheric Loading & Troposphere

Residual atmospheric pressure loading and tropospheric delays are leading candidate systematics for low-frequency structure in SLR residuals and require targeted controls. Several observed properties are suggestive of a propagation-linked component, but are not, by themselves, definitive:

- **Spatial Scale Considerations:** Many atmospheric fields exhibit synoptic correlation lengths of order $\sim 500\text{--}1,000$ km. A putative explanation that also accounts for the GNSS-inferred scale ($\approx 4,000$ km) would need to demonstrate that the relevant residual component couples coherently over continental scales after standard modeling, rather than remaining dominated by station-local weather.

- **Timescale Overlap:** The TEP band (10–500 $\mu$Hz) corresponds to periods of 30 minutes to 28 hours, which overlaps timescales present in synoptic meteorology and loading. Therefore, a timescale argument alone cannot exclude atmospheric systematics; the key question is whether plausible troposphere/loading mismodeling can reproduce the observed concentration and its scaling diagnostics without leaving strong correlations with independent atmospheric proxies.

- **Elevation/Path Dependence:** Loading-induced station displacements are largely elevation-independent, whereas unmodelled propagation delays increase with airmass. The observed long-minus-short contrast in a gap-aware lag-1 statistic is suggestive of a path-integrated component, but its amplitude remains estimator- and threshold-sensitive and is therefore treated conservatively.

Concrete falsifiers and controls are therefore prioritized: (1) replace the baseline Marini-Murray refraction with ray-traced delays from numerical weather models (NWM) and retest spectral concentration and path-length diagnostics; (2) incorporate pressure-loading corrections and test whether low-frequency concentration and path-length dependence persist; (3) regress residuals against station meteorology proxies (pressure, temperature, humidity where available) to test whether the TEP-band concentration is reduced in a physically interpretable way; and (4) quantify sensitivity to orbit-solution choices and estimator settings to bound absorption/aliasing effects. These controls will determine whether the residual structure is primarily atmospheric/geodetic or whether a technology-orthogonal propagation-linked component remains after state-of-the-art corrections.

#### Satellite-Specific Asymmetry & Center-of-Mass

The stark contrast in statistical significance between LAGEOS-2 ($p=0.0005$) and LAGEOS-1 ($p=0.93$) warrants careful scrutiny. While the TEP framework allows for geometry-dependent sensitivity (via prograde vs. retrograde sampling of the scalar gradient), this asymmetry raises the possibility of satellite-specific systematics. LAGEOS-2, while effectively identical in design to LAGEOS-1, has a distinct thermal and rotational history. Unmodelled center-of-mass (CoM) variations linked to spin-axis orientation or specific retroreflector array optical properties could, in principle, introduce structured range biases. If such biases couple with the specific network geometry viewing LAGEOS-2, they might masquerade as spatial correlations. Although standard CoM corrections (251 mm) were applied uniformly, further investigation using precise spin-vector evolution models is required to assess array-response artifacts as a potential source of the LAGEOS-2 signal. Crucially, however, thermal thrusting and spin-vector artifacts are typically associated with orbital-arc timescales or station-specific passes, which would tend to decorrelate inter-station noise. In contrast, the signal observed here exhibits a coherent distance-dependent structure with a correlation length of $\approx 4,000$ km—a spatial scale difficult to reproduce with purely local spacecraft systematics.

#### Geodetic Anomalies & Reinterpretation

The geodetic literature records several persistent "anomalies" or unexplained systematics that align with the TEP phenomenology reported here. Rather than viewing these as isolated technical errors, they may be reinterpreted as independent detections of the same underlying conformal field structure.

**1. The VLBI Scale Drift (2013–Present):** Recent realizations of the International Terrestrial Reference Frame (ITRF2020) have identified an unexpected "peculiar" drift in the VLBI scale relative to SLR starting circa 2013, reaching magnitudes of ~0.2 ppb (~1.2 mm) (Kern et al., 2024; Altamimi et al., 2023). While conventionally attributed to localized station deformations (e.g., at Ny-Ålesund) or "unexplained instrumental drift," the global, secular nature of the divergence suggests a fundamental metric phenomenon. In the TEP framework, a cosmological time-evolution of the background scalar field $\phi(t)$ implies a drift in the conformal factor $A(\phi)$. This acts as a dynamic rescaling of proper lengths, naturally manifesting as a scale drift between dynamic (SLR, governed by gravitational potential $GM$) and geometric (VLBI) reference frames.

**2. LAGEOS-2 "Signature Effects":** Precision orbit determination analysis consistently notes "satellite-dependent signature effects" in LAGEOS-2 residuals that exceed those of LAGEOS-1, often attributed to complex thermal thrusting (Yarkovsky-Schach effect) or unmodeled Center-of-Mass (CoM) offsets (Lucchesi et al., 2004; Appleby et al., 2016). This study's finding—that the TEP spatial correlation is driven exclusively by LAGEOS-2 ($p=0.0005$)—suggests that these "signature effects" may not be purely local spacecraft systematics. Instead, they likely reflect enhanced coupling between the scalar field gradient and the prograde orbital geometry of LAGEOS-2. The TEP correlation provides a coherent physical mechanism for these "anomalous" residuals without requiring ad-hoc thermal tunings.

**3. The "Flicker Noise" Floor & Common Mode Errors:** Geodetic analyses routinely identify "flicker noise" ($1/f$) spectra in station coordinate time series and "spatially correlated errors" (CME) that limit network stability to ~3–5 mm (Williams et al., 2004). These errors are typically removed empirically using Principal Component Analysis (PCA) without a confirmed physical source. The TEP analysis demonstrates that such spatially coherent, colored noise is a *prediction* of the theory ($\lambda \approx 4,000$ km), not merely atmospheric residue. The standard practice of filtering these signals effectively "bleaches" the conformal structure from geodetic products.

#### Network Sparsity

The ILRS network is significantly sparser than the IGS (GNSS) network. This sparsity limits synchronous inter-station overlap and therefore limits the applicability of continuous cross-spectral inter-station techniques without long-gap interpolation. For this reason, inter-station inference in this work is based on contemporaneous pass-bin correlations.

#### Orbit Modeling Dependencies

The estimator-dependent sign-inversion hypothesis ("Anti-Echo") relies on the specific behavior of the least-squares orbit determination filter. While the simulation (Section 4.3) demonstrates the plausibility of sign flipping under monopole absorption, empirical testing requires: (1) denser overlap epochs enabling synchronous multi-station observations, (2) parallel processing with kinematic and dynamic estimators to isolate the sign-flip signature, and (3) cross-Analysis-Center orbit comparisons (ASI vs. GFZ vs. CSR) to test monopole absorption consistency. These requirements are not met by current ILRS operations and represent a roadmap for future targeted SLR campaigns.

## 6. Conclusion

This study presents an independent optical-domain test of the Temporal Equivalence Principle's conformal sector using Satellite Laser Ranging. By analyzing 11 years of LAGEOS data, a spatially coherent residual structure has been isolated that mirrors the conformal signatures observed in microwave GNSS networks (Smawfield 2025b, 2025c, 2025d; Papers 1-3). The analysis identifies a path-length dependent temporal contrast in a gap-aware lag-1 diagnostic (treated conservatively due to threshold sensitivity) and a strong spectral concentration (14.00× enhancement in the TEP band relative to broadband, 95% CI: 13.53–14.47). The observation of matching low-frequency structure in a system devoid of active clocks and microwave transmission challenges hypotheses that rely on microwave-specific receiver electronics, clock steering, or ionospheric modeling.

Taken together with the GNSS results, the SLR evidence is consistent with the conformal-sector interpretation across five orders of magnitude in frequency (microwave to optical) and across two independent measurement technologies (GNSS and SLR). The observed correlation length in GNSS ($\lambda = 4,201 \pm 1,967$ km) and the SLR path-length dependent diagnostic (over a $\sim 3{,}000$ km path difference between the high- and low-range selections) are both consistent with Vainshtein screening at the dark energy scale ($m_\phi \approx 5 \times 10^{-14}\,\mathrm{eV}/c^2$).

In the inter-station domain, the sparse ILRS cadence limits synchronous overlap; however, by employing daily aggregation and widened contemporaneous windows, statistically significant distance-structured correlations have been detected ($p=0.004$ for pass-bin analysis; $p=0.017$ for daily aggregation). While the signal is complex and geometry-dependent (driven by LAGEOS-2), its detection indicates that the low-frequency structure observed in single-station spectra is spatially coherent. Furthermore, this phenomenology offers a unified physical explanation for persistent geodetic anomalies, including the ITRF2020 VLBI-SLR scale drift and the pervasive "flicker noise" floor in station coordinates, reinterpreting them as evidence of conformal metric coupling rather than intractable systematics.

The evidence presented here is consistent with a conformal-sector signal across two independent measurement systems, multiple processing centers, long-term GNSS stability, raw observational data, and the full electromagnetic spectrum. The "Time Echo" is not readily explained as a single-technology artifact; it appears as a reproducible low-frequency structure in geodetic residuals consistent with universal conformal coupling to a dynamical time field. Future experiments—closed-loop optical time transfer, interplanetary one-way asymmetry measurements, and triangle synchronization holonomy tests—will determine whether the disformal sector also manifests in nature or remains a theoretical possibility bounded to negligibility by multi-messenger constraints.

## References

### TEP Research Program

Smawfield, M. L. (2025a). *Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed*. Zenodo (Preprint). [DOI: 10.5281/zenodo.16921911](https://doi.org/10.5281/zenodo.16921911) (Paper 0: Theory)

Smawfield, M. L. (2025b). *Global Time Echoes: Distance-Structured Correlations in GNSS Clocks*. Zenodo (Preprint). [DOI: 10.5281/zenodo.17127229](https://doi.org/10.5281/zenodo.17127229) (Paper 1: Multi-Center)

Smawfield, M. L. (2025c). *Global Time Echoes: 25-Year Temporal Evolution of Distance-Structured Correlations in GNSS Clocks*. Zenodo (Preprint). [DOI: 10.5281/zenodo.17517141](https://doi.org/10.5281/zenodo.17517141) (Paper 2: Longspan)

Smawfield, M. L. (2025d). *Global Time Echoes: Raw RINEX Validation of Distance-Structured Correlations in GNSS Clocks*. Zenodo (Preprint). [DOI: 10.5281/zenodo.17860166](https://doi.org/10.5281/zenodo.17860166) (Paper 3: Raw RINEX)

Smawfield, M. L. (2025e). *Temporal-Spatial Coupling in Gravitational Lensing: A Reinterpretation of Dark Matter Observations*. Zenodo (Preprint). [DOI: 10.5281/zenodo.17982540](https://doi.org/10.5281/zenodo.17982540) (Paper 4: TEP-GL)

Smawfield, M. L. (2025f). *Global Time Echoes: Empirical Validation of the Temporal Equivalence Principle*. Zenodo (Preprint). [DOI: 10.5281/zenodo.18004832](https://doi.org/10.5281/zenodo.18004832) (Synthesis: TEP-GTE)

Smawfield, M. L. (2025g). *Universal Critical Density: Unifying Atomic, Galactic, and Compact Object Scales*. Zenodo (Preprint). [DOI: 10.5281/zenodo.18064366](https://doi.org/10.5281/zenodo.18064366) (Paper 7: TEP-UCD)

Smawfield, M. L. (2025h). *The Soliton Wake: A Runaway Black Hole as a Gravitational Soliton*. Zenodo (Preprint). [DOI: 10.5281/zenodo.18059251](https://doi.org/10.5281/zenodo.18059251) (Paper 8: TEP-RBH)

Smawfield, M. L. (2025i). *Global Time Echoes: Optical Validation of the Temporal Equivalence Principle via Satellite Laser Ranging*. Zenodo (Preprint). [DOI: 10.5281/zenodo.18064582](https://doi.org/10.5281/zenodo.18064582) (Paper 9: This Work)

### SLR & Geodesy References

Pearlman, M. R., Degnan, J. J., & Bosworth, J. M. (2002). *The International Laser Ranging Service*. Advances in Space Research, 30(2), 135-143. [DOI: 10.1016/S0273-1177(02)00277-6](https://doi.org/10.1016/S0273-1177(02)00277-6)

Marini, J. W., & Murray, C. W. (1973). *Correction of Laser Range Tracking Data for Atmospheric Refraction at Elevations Above 10 Degrees*. NASA Technical Report X-591-73-351.

Mendes, V. B., & Pavlis, E. C. (2004). *High-accuracy zenith delay prediction at optical wavelengths*. Geophysical Research Letters, 31(14). [DOI: 10.1029/2004GL020308](https://doi.org/10.1029/2004GL020308)

International Laser Ranging Service (ILRS). (2025). *SLR monthly normal point data*. Greenbelt, MD, USA: NASA Crustal Dynamics Data Information System (CDDIS). [DOI: 10.5067/SLR/slr_data_monthly_npt_001](https://doi.org/10.5067/SLR/slr_data_monthly_npt_001)

International Laser Ranging Service (ILRS). (2025). *SLR Combination Center (CC) Orbit Product*. Greenbelt, MD, USA: NASA Crustal Dynamics Data Information System (CDDIS). [DOI: 10.5067/SLR/SLR_ILRSORB_001](https://doi.org/10.5067/SLR/SLR_ILRSORB_001)

International Laser Ranging Service (ILRS). (2025). *Satellite Laser Ranging Frame 2020 (SLRF2020) Product*. Greenbelt, MD, USA: NASA Crustal Dynamics Data Information System (CDDIS). [DOI: 10.5067/SLR/slrf2020_001](https://doi.org/10.5067/SLR/slrf2020_001)

Hellmers, H., et al. (2025). *Terrestrial reference frame scale drift anomalies in VLBI and the contribution of Ny-Ålesund radio telescopes*. Earth, Planets and Space, 77. [DOI: 10.1186/s40623-025-02159-z](https://doi.org/10.1186/s40623-025-02159-z)

Lucchesi, D. M., et al. (2004). *Spin axis behavior of the LAGEOS satellites*. Journal of Geophysical Research: Solid Earth, 109(B6). [DOI: 10.1029/2003JB002692](https://doi.org/10.1029/2003JB002692)

Appleby, G., Rodríguez, J., & Altamimi, Z. (2016). *Assessment of the accuracy of global geodetic satellite laser ranging observations and estimated impact on ITRF2014*. Journal of Geodesy, 90, 1371–1388. [DOI: 10.1007/s00190-016-0929-2](https://doi.org/10.1007/s00190-016-0929-2)

Kern, L., et al. (2024). *Verifying the impact of additional breaks in station coordinates on VLBI scale drift*. Technical Report, Vienna University of Technology.

Altamimi, Z., et al. (2023). *ITRF2020: an augmented reference frame refining the modeling of nonlinear station motions*. Journal of Geodesy, 97, 1005–1023. [DOI: 10.1007/s00190-023-01738-w](https://doi.org/10.1007/s00190-023-01738-w)

Williams, S. D. P., et al. (2004). *Error analysis of continuous GPS position time series*. Journal of Geophysical Research: Solid Earth, 109(B3). [DOI: 10.1029/2003JB002741](https://doi.org/10.1029/2003JB002741)

Jiang, M., et al. (2023). *Long-Baseline Quantum Sensor Network as Dark Matter Haloscope*. Nature Communications / arXiv:2305.00890. [arXiv:2305.00890](https://doi.org/10.48550/arXiv.2305.00890)

---

*This document was automatically generated from the TEP-SLR research site. For the interactive version with figures and enhanced formatting, visit: https://mlsmawfield.com/tep/slr/*

*Related Work:*
- [**TEP Theory**](https://doi.org/10.5281/zenodo.16921911) (Foundational framework)
- [**TEP-GTE**](https://doi.org/10.5281/zenodo.17127229) (GNSS Validation)

*Source code and data available at: https://github.com/matthewsmawfield/TEP-SLR*
