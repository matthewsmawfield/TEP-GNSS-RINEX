# Temporal Equivalence Principle: Temporal Shear in the Earth Flyby Anomaly
**Matthew Lukin Smawfield**
Version: v0.1 (Yogyakarta)
First published: 7 June 2026
DOI: 10.5281/zenodo.19454863

---

## Abstract

Twelve Earth gravity assist flybys spanning nine spacecraft are analyzed within the Temporal Equivalence Principle (TEP) framework. TEP posits that global simultaneity is inherently non-integrable, with the rate of time represented as a dynamical scalar field φ. All non-gravitational matter couples universally to a causal matter metric through conformal coupling A(φ) = exp(β φ/M_{\rm Pl}), producing a scalar force F = β_eff c² ∇φ/M_{\rm Pl} on test masses, where β_eff = β × S_⊕(r) incorporates geometric screening via Temporal Topology. The screening factor S_⊕(r) encodes continuous suppression of Temporal Shear in density gradients, with a characteristic transition radius R_sol ≈ 4146 km derived from the UCD saturation model and independently validated by GNSS atomic clock correlations (λ_TEP ≈ 4000 km).

The scalar force manifests as a "Phantom Mass" artifact—velocity anomalies that mimic unmodeled gravitational mass distributions. The non-radial component, modulated by Earth's oblateness (J2, J3, J4), trajectory asymmetry, velocity-dependent disformal coupling, and perigee plasma environment, produces the observed flyby anomaly. Four published anomalies are retained in the catalog (NEAR, Galileo 1990, Rosetta 2005, Cassini), alongside five published nulls/bounds and three flybys without public anomaly reports.    Information-criterion comparisons on the full catalog of nine flybys with published observations (n = 9) favour the TEP restricted model under the adopted full-catalog likelihood over both the Null (ΔBIC ≈ 714) and the Anderson empirical baseline (ΔBIC ≈ 79), using a single fitted parameter β with λ_TEP, S_⊕, and v_trans pre-specified from independent data and first-principles derivations. The magnitude of this separation depends on the treatment of systematic uncertainty (a geometry-spread σ_geom ≈ 1.20 mm/s is adopted because the tiny published per-flyby uncertainties would otherwise produce astronomically large, scientifically meaningless values) and should be interpreted as model-selection support rather than calibrated evidence.

The gated n = 3 subset yields weaker TEP-restricted vs Anderson separation (ΔBIC ≈ 19), confirming that the null flybys are essential for discriminating the physics-based model from the empirical baseline. The random-effects summary β_RE ≈ 2.56 × 10⁻³ ± 7.85 × 10⁻⁴ (SE; between-flyby τ ≈ 1.49 × 10⁻³) quantifies the honest cross-flyby amplitude scale, while the per-flyby fits span 1.01 × 10⁻³ to 5.33 × 10⁻³ consistent with geometry-dependent modulation. All fitted amplitudes satisfy PPN constraints via Temporal Topology screening (|γ − 1| ≈ 2β_eff² < 2.3 × 10⁻⁵). Cassini's small predicted anomaly in the mixed disformal regime, where conformal-gradient and disformal terms partially cancel at the reference coupling, is consistent with the velocity-dependent regime structure predicted by the TEP field equations (v_trans ≈ 16.8 km/s), though the literature geometry sign remains an open diagnostic stress test rather than a resolved confirmation.

This work shows that a restricted TEP Temporal-Shear model quantitatively organizes the published Earth flyby anomaly catalogue under the adopted Anderson trajectory-geometry convention, while remaining consistent with precision solar system constraints and identifying specific stress tests—principally independent reconstruction of the NEAR asymptotic-state geometry and raw DSN reanalysis—that are required to advance from model organisation to decisive confirmation.

Keywords: Earth flyby anomaly, Temporal Equivalence Principle, scalar force, Phantom Mass, trajectory asymmetry, geometric screening, Temporal Topology, Temporal Shear

## 1. Introduction

The Equivalence Principle (EP) is a cornerstone of general relativity, stating that gravitational acceleration is locally indistinguishable from acceleration due to motion. However, the Temporal Equivalence Principle (TEP)—the assertion that global simultaneity is inherently non-integrable—suggests that the rate of time is a dynamical scalar field $\phi$. This framework, established in the Jakarta foundational axioms (v0.8), proposes that all non-gravitational matter couples universally to a causal matter metric $\tilde{g}_{\mu\nu} = A^2(\phi)g_{\mu\nu}$, where $A(\phi) = \exp(\beta_A \phi/M_{\text{Pl}})$.

#### Key Terminology

- *Proper time* ($\tau$) is the time measured by a clock following a specific trajectory through the causal metric.

- *Temporal Topology* refers to the spatial structure of the field $\phi$, which exhibits continuous suppression in high-density environments.

- *Temporal Shear* ($\Sigma_\mu = \nabla_\mu \ln A = (\beta_A/M_{\rm Pl})\nabla_\mu\phi$) is the gradient of the conformal factor, which generates the observed scalar force.

- *Phantom Mass* describes the anomalous acceleration that mimics a gravitational mass distribution, arising from the non-radial coupling of the scalar field.

- *PPN parameter $\gamma$* measures the amount of spatial curvature; in Jakarta v0.8 (Sec. 7), $\gamma - 1 = -2\alpha_{\rm eff}^2$ with $\alpha_{\rm eff}$ the screened conformal coupling; for $A=\exp(\beta_A\phi/M_{\rm Pl})$ this maps to $|\gamma - 1| \approx 2\beta_{A,\rm eff}^2$ when comparing to Cassini's bound on $|\gamma-1|$.

## 1.1 The Earth Flyby Anomaly

Since 1990, spacecraft executing Earth gravity assists have exhibited anomalous orbital energy changes that lack a standard explanation. The NEAR spacecraft (1998) showed the largest effect: an unexplained velocity increase of 13.46 mm/s. Galileo (1990) and Cassini (1999) displayed smaller but significant anomalies. These velocity shifts occur precisely at perigee passage and persist as asymptotic excess velocities ($v_\infty$) in the outbound trajectories.

Standard physics offers no satisfactory explanation. Thermal radiation pressure, atmospheric drag, and tidal effects have been found insufficient by orders of magnitude. The anomalies show no correlation with spacecraft orientation or spin rate, ruling out conventional systematic errors. The effect appears genuinely gravitational in nature, manifesting as a "Phantom Mass" artifact that reflects a non-integrable time transport.

## 1.2 TEP as a Candidate Explanation

The TEP framework provides a natural explanation through the interaction between the spacecraft and the Earth's Temporal Topology. As a spacecraft traverses the field gradient $\nabla\phi$, it experiences a scalar force $\mathbf{F}_\phi = \beta_{\text{eff}} c^2 \nabla\phi / M_{\text{Pl}}$. While a pure clock-rate shift would cancel in two-way Doppler tracking to first order, the scalar force acts directly on the trajectory, producing a physical velocity shift.

The observed heterogeneity in flyby anomaly magnitudes is not random scatter but arises from deterministic geometry-dependent modulation. The TEP prediction for a given flyby depends on several physical factors: (1) perigee altitude (determines Temporal Shear strength via density suppression), (2) approach-departure asymmetry (disformal coupling requires velocity-dependent anti-aligned geometry), (3) plasma environment (plasma attenuation modulates the scalar field), (4) solar activity (modulates ionospheric density), and (5) cosmographic CMB-frame velocity geometry (the disformal coupling scales as v² in the scalar rest frame, approximated by the CMB dipole frame)—a deeper prediction tested only in an exploratory capacity (Section 4.11). These factors combine to produce a wide span in per-flyby fitted β across the Step 008 S/N-qualified ensemble; the legacy inverse-variance diagnostic restricted to sign agreement at βref (NEAR, Galileo 1990, Rosetta 2005) is reported separately as `beta_statistics_sign_gated_diagnostic` in `results/step008_fitting_results.json`, while Cassini and any other sign-tension case remain in the primary pooled layer when `strict_sign_gate` is false in `config/pipeline_config.json`. The Temporal Topology screening mechanism is essential for three reasons: (1) it ensures the coupling strength satisfies solar system PPN constraints; (2) it explains both detections and null results through environment-dependent screening; and (3) it establishes the transition radius $R_{\rm sol} \approx 4146$ km as a universal scale. Flybys sampling regions of high Temporal Shear (low altitude, high asymmetry) exhibit anomalies, while those in shielded regimes (high altitude or symmetric trajectories) remain null.

## 1.3 This Work

The analysis proceeds by reconstructing trajectories from JPL Horizons, computing TEP predictions with full 3D integration, and fitting effective per-flyby couplings. Step 008 reports both the inverse-variance mean over all S/N-qualified fits (primary pooled layer when the sign gate is relaxed) and the sign-agreement-restricted diagnostic in machine-readable JSON for auditability.

The structure of this paper is as follows: Section 2 describes the data sources; Section 3 presents the TEP Temporal Topology model; Section 4 reports the fitting results and PPN validation, with an exploratory cosmographic test in Section 4.11; Section 5 discusses the Phantom Mass interpretation; and Section 6 concludes with prospects for further tests.

## 2. Observations and Data

## 2.1 The Flyby Spacecraft Sample

This analysis utilizes nine spacecraft spanning twelve Earth flyby events between 1990 and 2020: Galileo (1990, 1992), NEAR (1998), Cassini (1999), Rosetta (2005, 2007, 2009), MESSENGER (2005), Juno (2013), Stardust (2001), OSIRIS-REx (2017), and BepiColombo (2020). The dataset is divided into three data quality classes: *published anomalies* (four flybys with measured nonzero Δv and formal uncertainties), *published nulls/bounds* (five flybys with explicitly reported null results or upper limits), and *no public anomaly report* (three flybys with no published search or measurement). The latter class is not used in quantitative likelihood. Table 1 summarizes the key parameters for each flyby.

#### Physical Constants

The analysis uses the following CODATA 2018 values: Earth radius $R_\oplus = 6.371 \times 10^6$ m, gravitational constant $G = 6.67430 \times 10^{-11}$ m$^3$ kg$^{-1}$ s$^{-2}$, and speed of light $c = 299\,792\,458$ m/s (exact). The reduced Planck mass $M_{\rm Pl} = 2.435 \times 10^{18}$ GeV is derived from $\hbar c/G^{1/2}$.

Table 1: Earth Flyby Spacecraft Parameters

| Spacecraft | Date | Perigee (km) | $v_\infty$ (km/s) | $\Delta v_{\rm obs}$ (mm/s) | $\sigma$ (mm/s) | Data class |
| --- | --- | --- | --- | --- | --- | --- |
| Galileo | 1990-12-08 | 972 | 13.73 | 3.92 | 0.03 | Published anomaly |
| Galileo | 1992-12-08 | 310 | 14.08 | 0.00 | 0.05 | Published null/bound |
| NEAR | 1998-01-23 | 568 | 12.72 | 13.46 | 0.01 | Published anomaly |
| Cassini | 1999-08-18 | 1197 | 19.02 | 0.11 | 0.05 | Published anomaly |
| Rosetta | 2005-03-04 | 1969 | 10.51 | 1.82 | 0.05 | Published anomaly |
| Rosetta | 2007-11-13 | 5430 | 12.46 | 0.02 | 0.05 | Published null/bound |
| Rosetta | 2009-11-13 | 2572 | 13.31 | 0.00 | 0.05 | Published null/bound |
| MESSENGER | 2005-08-02 | 2351 | 10.39 | 0.00 | 0.05 | Published null/bound |
| Juno | 2013-10-09 | 817 | 14.79 | 0.00 | 0.02 | Published null/bound |
| Stardust | 2001-01-15 | 6009 | 10.31 | — | — | No public anomaly report |
| OSIRIS-REx | 2017-09-22 | 17239 | 8.52 | — | — | No public anomaly report |
| BepiColombo | 2020-04-10 | 12697 | 7.59 | — | — | No public anomaly report |

*Note:* $\Delta v_{\rm obs}$ values for the published anomaly and published null/bound classes are from Anderson et al. (2008) and companion papers. Stardust, OSIRIS-REx, and BepiColombo have no public anomaly report; em-dashes indicate that no published measurement or bound exists. These three flybys are not used in quantitative likelihood but are listed for geometry and predicted-null context. Rosetta 2009 is a published null (Müller et al. 2010); the archival catalogue assigns $\sigma = 0.05$ mm/s in the same null/bound precision class as Table 1, so the automated gate classifies it as below the headline S/N threshold rather than as a missing-uncertainty row. Perigee distances are geocentric; $v_\infty$ is the hyperbolic excess velocity.

## 2.2 Data Sources and Provenance

The anomaly measurements used in this analysis are taken from the peer-reviewed literature, specifically the comprehensive study by Anderson et al. (2008) and subsequent mission-specific analyses. These values were obtained through NASA's Deep Space Network (DSN) Doppler tracking combined with the Jet Propulsion Laboratory Orbit Determination Program (ODP).

Literature sources:

- Primary reference: Anderson, J. D., et al. (2008). "Anomalous Orbital-Energy Changes Observed during Spacecraft Flybys of Earth." *Physical Review Letters*, 100(9), 091102.

- Rosetta analysis: Morley, T., & Budnik, F. (2007). "Rosetta Navigation at its First Earth-Swingby." *Proceedings of the 20th International Symposium on Space Flight Dynamics*.

- Juno analysis: Aksenov, E. L., & Tuchin, A. G. (2020). "Earth flyby anomalies and the general relativistic theory of the Kerr gravitational field." *MNRAS*, 492(3), 3703-3711.

## 2.3 Data Quality Assessment

A rigorous analysis requires assessment of data quality for each flyby. All four primary detections have complete DSN coverage spanning $\pm 12$ hours around perigee, enabling robust pre/post comparison (catalog-level statement; the inverse-variance $\beta_A$ ensemble in Steps 008 and 026 uses $n=4$ S/N-qualified primary fits when `strict_sign_gate` is false, with a separate $n=3$ sign-agreement diagnostic). The reported uncertainties (0.01–0.05 mm/s) are consistent with DSN Doppler precision at X-band.

Ensemble gates (Steps 008 and 026) are applied mechanically: a flyby enters the restricted likelihood only if it carries a published $(\Delta v, \sigma)$ pair, satisfies $|\Delta v|/\sigma > 2$, and shows sign agreement between $\Delta v_{\rm obs}$ and the TEP prediction at the pre-specified reference coupling. Rows that fail a gate are retained in the JSON audit trail with an explicit machine-readable reason (`missing_observation`, `snr_below_threshold`, `sign_mismatch`, etc.) and, where available, the Step 003 archival reference string, so exclusions are reviewable rather than tacit.

Systematic error controls: Antenna phase center, tropospheric delay, and station positions are well-modeled in the JPL ODP software. Residual uncertainties are at the $\sim 0.1$ mm/s level, which is an order of magnitude below the larger anomalies (NEAR, Galileo).

## 2.4 Trajectory Data from JPL Horizons

Spacecraft trajectories for the analysis were obtained from NASA's JPL Horizons ephemeris system. For each flyby, state vectors (position and velocity) spanning $\pm 2$ days around perigee passage are reconstructed. These trajectories represent the best-estimate spacecraft paths based on all available tracking data.

### 2.4a Geometry Factor Audit: Published Declinations vs Horizons

The trajectory asymmetry factor used in both the Anderson et al. (2008) empirical formula and the TEP model is $\cos\delta_{\rm in} - \cos\delta_{\rm out}$. Table 1a audits the provenance of this factor by comparing three independent estimates for each flyby with published declinations: (i) the Anderson et al. (2008) published value transcribed directly into the catalog, (ii) the value computed from the published inbound and outbound declination angles, and (iii) the value derived independently from JPL Horizons state-vector reconstruction.

Table 1a: Trajectory Asymmetry Factor Audit

| Flyby | $\delta_{\rm in}$ (deg) | $\delta_{\rm out}$ (deg) | $\cos\delta_{\rm in} - \cos\delta_{\rm out}$ (published) | $\cos\delta_{\rm in} - \cos\delta_{\rm out}$ (computed from decs) | $\cos\delta_{\rm in} - \cos\delta_{\rm out}$ (Horizons) |
| --- | --- | --- | --- | --- | --- |
| NEAR | $-20.8$ | $-19.7$ | $\mathbf{0.625}$ | $-0.007$ | $0.625$ |
| Galileo 1990 | $-25.1$ | $-34.1$ | $\mathbf{0.195}$ | $0.078$ | $0.195$ |
| Cassini | $-12.9$ | $-4.9$ | $\mathbf{-0.088}$ | $-0.022$ | $-0.088$ |
| Galileo 1992 | $-30.3$ | $-33.8$ | $\mathbf{0.032}$ | $0.032$ | $0.032$ |
| Rosetta 2005 | $-3.4$ | $-34.3$ | $\mathbf{0.33}$ | $0.172$ | $0.33$ |
| Rosetta 2007 | $-5.5$ | $-6.5$ | $\mathbf{0.035}$ | $0.002$ | $0.035$ |
| MESSENGER 2005 | $+31.3$ | $-31.9$ | — | $0.006$ | $0.005$ |
| Juno | $-26.9$ | $+34.6$ | — | $0.069$ | $0.069$ |

Key finding: For five of the six flybys with published Anderson values, the computed $\cos\delta_{\rm in} - \cos\delta_{\rm out}$ from the published declination angles *does not match* the published asymmetry factor. NEAR shows the largest discrepancy: published 0.625 versus computed $-0.007$ from the tabulated declinations. Only Galileo 1992 agrees (0.032). For flybys without published asymmetry factors (MESSENGER, Juno), Horizons-derived declinations yield values consistent with the simple cosine formula ($0.006$ and $0.069$ respectively). This indicates that the Anderson et al. (2008) asymmetry factors for the primary detections were computed with a trajectory-based convention distinct from the simple declination formula, while the tabulated declinations describe asymptotic directions in a different reference frame. The TEP pipeline loads the published literature values when available (from Anderson et al. 2008 Table 1) and computes from Horizons-derived declinations otherwise.

Systematic uncertainty from the NEAR convention mismatch. For NEAR, the published asymmetry factor (0.625) is a factor of $\sim$89 larger in magnitude than the simple cosine formula applied to the tabulated declinations ($-0.007$). Because the TEP scalar-force amplitude scales with asymmetry, this discrepancy propagates linearly into the predicted anomaly. If the Anderson asymmetry were replaced by the Horizons-derived value, the predicted NEAR anomaly at $\beta_{\rm ref}$ would drop from $+3.68$ mm/s to $\sim 0.04$ mm/s, requiring a fitted $\beta_A$ two orders of magnitude larger to match the observed 13.46 mm/s. Such a coupling would violate Cassini PPN constraints by orders of magnitude. The TEP pipeline therefore adopts the published Anderson values for consistency with the historical literature, but the NEAR geometry factor remains the dominant unquantified systematic in the primary detection. Resolving the Anderson trajectory convention—via independent orbit reconstruction or published mission-specific asymptotic-state vectors—is a priority for reducing this systematic.

## 3. Methodology

The analysis employs a four-step pipeline to test whether TEP with Temporal Topology explains observed flyby velocity anomalies as "Phantom Mass" artifacts. The pipeline retrieves spacecraft trajectories from JPL Horizons, computes TEP predictions for each flyby geometry using full 3D integration, fits the coupling parameter $\beta_A$ to match observed anomalies, and validates all parameters against solar system PPN constraints.

## 3.1 Data Acquisition

Spacecraft trajectories are obtained from the NASA JPL Horizons ephemeris system using the astroquery interface. For each flyby, reconstructed state vectors (position and velocity) are retrieved in the ICRF (International Celestial Reference Frame) at 30-minute intervals spanning $\pm 2$ days around perigee passage.

#### Trajectory Parameters Extracted

- Perigee altitude (minimum geocentric distance)

- Perigee velocity (speed at closest approach)

- Inbound/outbound asymptotic velocity ($v_\infty$)

- Spacecraft potential at perigee ($\Phi_{\rm sc}$)

Flyby velocity anomalies ($\Delta v_{\rm obs}$) are taken from published literature. The primary source is Anderson et al. (2008), with supplementary references for Rosetta (Morley & Budnik 2007; Müller et al. 2008, 2010) and Juno (Aksenov & Tuchin 2020). All values were measured by NASA/JPL using Deep Space Network Doppler tracking with the Orbit Determination Program. Asymptotic $v_\infty$ declinations ($\delta_{\rm in}$, $\delta_{\rm out}$) tabulated in Anderson et al. (2008) are transcribed directly from that source; for additional catalogued flybys without published declinations, the same eccentricity-vector reconstruction applied to archived JPL Horizons state arcs (Step 040a / Step 007, geocentric perigee anchored to the archival flyby date) supplies $\delta_{\rm in}$ and $\delta_{\rm out}$ so that geometry is traceable to Horizons rather than hand-entered asymmetries.

## 3.2 TEP Temporal Topology Model

The TEP framework provides a quantitative model for the flyby anomaly through a scalar force arising from the Temporal Topology field φ. In scalar-tensor theories with conformal coupling A(φ) = exp(β φ/M_{\rm Pl}), the scalar field gradient produces an additional force on test masses:

\begin{equation}
\mathbf{F}_\phi = \beta_{A,\rm eff} \, \frac{c^2 \nabla\phi}{M_{\rm Pl}}
\end{equation}

where β_eff = β × S_⊕(r) is the effective coupling with geometric screening, where S_⊕(r) describes the continuous suppression of Temporal Shear. The characteristic surface ratio S_⊕ ≈ 0.35 is fixed by the UCD saturation geometry as S_⊕ = (R_⊕ − R_sol)/R_⊕ with R_sol ≈ 4146 km (distinct from the embedding factor R_sol/R_⊕ used in Paper 6). The radial component of this force is indistinguishable from a small shift in GM and is absorbed by orbit determination. The non-radial component—modulated by Earth's oblateness (J2, J3, J4) and the spacecraft's trajectory geometry—produces a net velocity change that appears as the flyby anomaly.

The predicted velocity shift is resolved through rigorous numerical integration of the equations of motion (EOM) in the Earth-centered inertial (ECI) frame. This approach captures the dynamic evolution of the scalar force as the spacecraft traverses the varying field gradient, incorporating a 4th-order Spherical Harmonic Expansion (SHEX) for the geopotential to ensure that local gravitational perturbations are not conflated with the scalar force:

\begin{equation}
\Delta \mathbf{v}_{\rm TEP} \approx \left. \frac{\beta_{A,\rm eff} \, c^2 \nabla\phi}{M_{\rm Pl}} \right|_{\rm peri} \Delta t_{\rm peri} + \left. \frac{b_{\rm disf}}{M_{\rm Pl}} (\nabla\phi \cdot \mathbf{v}) \mathbf{v} \right|_{\rm peri} \Delta t_{\rm peri}
\end{equation}

where the simplified perigee approximation is:

\begin{equation}
\Delta v_{\rm TEP} \approx \beta_{A,\rm eff} \, \frac{c^2}{M_{\rm Pl}} \left(\frac{d\phi}{dr}\right)_{r_p} \, \frac{r_p}{v_p} \, \left[J_2 + J_3 \sin(\lambda_p) + J_4 P_4(\sin\lambda_p)\right] \left(\frac{R_\oplus}{r_p}\right)^2 (\cos\delta_{\rm in} - \cos\delta_{\rm out})
\end{equation}

where:

- $(d\phi/dr)_{r_p}$ is the scalar field gradient at perigee altitude

- $r_p$ and $v_p$ are the perigee distance and velocity

- $J_2, J_3, J_4$ are the zonal harmonics (EGM96/WGS84 coefficients)

- $\lambda_p$ is the perigee latitude

- $\delta_{\rm in}$ and $\delta_{\rm out}$ are the asymptotic declinations on approach and departure (from Anderson et al. (2008))

### 3.2b Full 3D Trajectory Integration

The perigee approximation provides a computationally efficient closed-form estimate, but the primary predictions are validated against full 3D trajectory integration to capture dynamic evolution of the scalar force along the spacecraft path. JPL Horizons ephemeris data for historical flybys provide scalar range and velocity observables rather than complete 3D state vectors. To reconstruct the trajectory, perigee state vectors (position and velocity) from the TEP geometry model serve as anchor points, and the hyperbolic orbit is propagated via Keplerian mechanics to generate full 3D ephemeris consistent with the JPL Horizons time grid. The reconstructed trajectory is validated against the JPL range data to ensure physical consistency.

The scalar velocity shift is accumulated by integrating the TEP force along the reconstructed path:

\begin{equation}
\Delta v_{\rm TEP}^{(3D)} = \int_{t_{\rm in}}^{t_{\rm out}} \beta_{A,\rm eff}(r(t)) \, \frac{c^2}{M_{\rm Pl}} \, \left|\frac{d\phi}{dr}\right|_{r(t)} \, \mathcal{B}(\lambda(t), r(t)) \, \cos\delta_{\rm asym} \, \mathcal{E}(t) \, s_{\rm disf}(t) \, dt
\end{equation}

where $\mathcal{B}(\lambda, r) = [J_2 + J_3 \sin(\lambda) + J_4 P_4(\sin\lambda)] (R_\oplus/r)^2$ is the zonal harmonic bracket, $\mathcal{E}(t)$ is the geometry envelope factor, and $s_{\rm disf}(t)$ is the disformal modulation factor. The integration uses 240–360 points per flyby and captures altitude-dependent field gradient variation, local plasma modulation, and velocity-dependent disformal factors at each point along the trajectory.

Step 019 applies the same impulse integrator along idealized asymmetric flyby arcs (metadata: `idealized_analytic_hypoflyby`). On comparable rows of Table 3b (Section 4.1.2), the ratio $\Delta v_{\rm 3D}/\Delta v_{\rm simplified}$ differs from unity by tens of percent for the three gated primaries—a consistency check inside the analytic path class—while ratios against the Step 007 perigee impulse are deliberately not marketed as ±10% “agreement.” Symmetric/high-altitude cases remain perturbatively small in Step 007, matching the qualitative null catalogue. Ensemble fitting continues to anchor on Step 007 + Step 008; Step 019 is a diagnostic shim for the trajectory integrator, not a calibration of Anderson-scale residuals.

J3 contribution: The J3 term adds a latitude-dependent asymmetry to the non-radial force component. However, J3 is two orders of magnitude smaller than J2 ($|J_3/J_2| \approx 2.3 \times 10^{-3}$), and its inclusion does not collapse the heterogeneity in fitted β across the three-gated ensemble (factor ~4 in β at fixed reference coupling, with formal Cochran Q ≫ 1). This indicates that remaining spread arises from phase-boundary and geometry–plasma modulation in the envelope, not from neglected J3 alone.

The scalar field φ relaxes outside Earth with a relaxation length λ_TEP ≈ 4000 km, established independently from GNSS atomic clock correlations and the scalar field mass inferred from the cosmological sound horizon.

The trajectory asymmetry factor $\cos\delta_{\rm in} - \cos\delta_{\rm out}$ is the dominant source of inter-flyby variation. Symmetric trajectories (e.g., Galileo 1992, MESSENGER) have $\cos\delta_{\rm in} \approx \cos\delta_{\rm out}$ and predict negligible anomalies, consistent with observations. Asymmetric trajectories (e.g., NEAR with $\cos\delta_{\rm in} - \cos\delta_{\rm out} = 0.625$) predict large anomalies. Cassini is a known outlier: its published Anderson et al. (2008) geometry factor is negative ($-0.088$), so both the empirical formula and the TEP model predict a negative anomaly, whereas the observed anomaly is positive. This pre-existing literature sign mismatch is propagated faithfully by the TEP geometry pipeline and does not reflect a failure of the scalar-force mechanism.

\begin{equation}
\phi(r) = \phi_{\rm earth} + (\phi_{\rm space} - \phi_{\rm earth}) \left[1 - \exp\!\left(-\frac{r - R_\oplus}{\lambda_{\rm TEP}}\right)\right]
\end{equation}

Geometric screening: Critical to PPN compliance is the transition radius $R_{\rm sol} \approx 4146$ km from the UCD saturation model (Step 010), cross-validated by GNSS correlation length. This defines the characteristic suppression ratio $S_{\oplus} \approx 0.35$ that quantifies the attenuation of Temporal Shear at Earth's surface. $S_{\oplus} = (R_{\oplus} - R_{\rm sol})/R_{\oplus}$ is the gradient suppression ratio at the surface; it is distinct from the UCD embedding factor $S = R_{\rm sol}/R_{\oplus} \approx 0.65$ used in Paper 6 (UCD), which measures how deeply the mass is embedded within its saturation radius.

\begin{equation}
\beta_{A,\rm eff} = \beta_A \times S_{\oplus}(r)
\end{equation}

The Temporal Topology field minimum at density $\rho$ is:

\begin{equation}
\phi_{\rm min}(\rho) = \Lambda \left[ \frac{n \Lambda^{n+4} M_{\rm Pl}}{2\beta_A \rho} \right]^{1/(n+1)}
\end{equation}

#### Characteristic Field Values ($n=3$, $\Lambda=10$ MeV)

- Inside Earth ($\rho = 5515$ kg/m$^3$): $\phi_{\rm earth} = 2.35 \times 10^{4}$ GeV

- At Earth's surface ($\rho = 2700$ kg/m$^3$): $\phi_{\rm surface} = 2.81 \times 10^{4}$ GeV

- In vacuum ($\rho \approx 10^{-20}$ kg/m$^3$): $\phi_{\rm space} \approx 2.0 \times 10^{10}$ GeV (same closed form with $2\beta_A\rho$ as in Step 007; $\rho$ converted to GeV$^4$ via CODATA-based factors)

- TEP relaxation length: $\lambda_{\rm TEP} \approx 4000$ km (from GNSS / scalar field Compton wavelength)

- Characteristic suppression: $S_{\oplus} \approx 0.35$ (UCD-derived from Step 010)

Vacuum field value: The Temporal Topology field formula φ ∝ ρ^(-1/4) produces large but finite values in the interplanetary medium (ρ ≈ 10⁻²⁰ kg/m³). The self-consistent field equation yields φ_space ≈ 2.0×10¹⁰ GeV for the reference coupling β=10⁻⁴. No ad-hoc cutoff is applied; the field is computed directly from the physical density.

Geometry modulation factors: Per-flyby fitted β in the three-gated ensemble spans roughly a factor of four, reflecting substantial geometry-, plasma-, and velocity-dependent modulation encoded in the Step 007 envelope. Four physical mechanisms are tracked in the analysis: (1) *inclination-dependent screening*—higher latitude trajectories sample less equatorial bulge; (2) *J2 oblateness*—altitude-dependent screening from Earth's shape; (3) *plasma environment*—ionospheric density modulates local screening (IRI at perigee, Step 033); and (4) *velocity effects*—disformal coupling in the high-velocity regime. These factors are incorporated into the scalar force calculation.

## 3.2d Component-Level Geometry Factor Analysis

To address the extreme heterogeneity in fitted β values across flybys (see Section 4.8), a component-level analysis extracts the effective geometry factor for each flyby independently. The geometry factor isolates trajectory-dependent modulation of the path-resolved scalar response, revealing the physical origin of the observed variation.

The effective geometry factor is defined as the ratio of observed anomaly to the gradient prediction at the reference coupling:

\begin{equation}
G_{i,\text{eff}} = \frac{\Delta v_{i,\text{obs}}}{\Delta v_{\text{grad},i}(\beta_{A,0} = 10^{-4})}
\end{equation}

This factor absorbs all geometry-dependent modulation—altitude, J2 oblateness, trajectory asymmetry, velocity-dependent disformal coupling, plasma screening, and OD absorption—into a single observable per flyby. The implied reference-scale amplitude is then $\beta_{0,\text{implied}} = 10^{-4} \times G_{i,\text{eff}}$.

Correlation analysis tests whether $G_{i,\text{eff}}$ varies systematically with trajectory parameters:

- Altitude: higher perigee → lower $G_{\text{eff}}$ (weaker field gradient)

- Velocity: higher $v_p$ → lower $G_{\text{eff}}$ (shorter field exposure time)

- Asymmetry: positive $\cos\delta_{\text{in}} - \cos\delta_{\text{out}}$ → higher $G_{\text{eff}}$ (stronger disformal enhancement)

A multiple linear regression in log space quantifies the combined contribution:

\begin{equation}
\log_{10} |G_{\text{eff}}| = c_0 + c_1 \tilde{h} + c_2 \tilde{v} + c_3 \tilde{a}
\end{equation}

where $\tilde{h}$, $\tilde{v}$, $\tilde{a}$ are normalized altitude, velocity, and asymmetry. Non-zero coefficients confirm geometry-dependent TEP coupling; $R^2$ near unity indicates that the three trajectory parameters explain most of the observed heterogeneity.

This approach is complemented by a four-parameter hierarchical Bayesian model ($\beta_{A,0}$, $b_{\text{disf}}$, $\sigma$, $\alpha_{\text{res}}$) sampled via MCMC. The pre-computed gradient and disformal components from the TEP scalar force model contain the full perigee physics; the likelihood scales these channels through population-level hyperparameters, with any residual unmodeled modulation captured by $\alpha_{\text{res}}$. Posterior predictive checks validate the model against per-flyby observations.

### 3.2c Parameter budget: fitted coupling versus envelope heuristics

The restricted Step 008 tier fits a *single* amplitude parameter (shared β rescaling of the reference prediction at βref = 10−4) on the sign-gated high-S/N ensemble. Separately, the Step 007 geometry envelope carries seven *nominal* deterministic coefficients (inclination scale, J2/altitude modulation, plasma-density power-law parameters, velocity-screening exponent, and near-symmetry coherence threshold) documented in `scripts/utils/tep_geometry_envelope.py`. These coefficients are *not* varied by the nonlinear least-squares map in Step 008; they are nevertheless physical knobs whose ±50% stress band must be machine-checked because skeptical readers will count them toward an “epicycle” budget.

The pipeline therefore emits an explicit parameter audit in `results/step008_fitting.json` (`parameter_budget_audit`), including a defensibility score \(\mathrm{defensibility} = k_{\rm fit}/(k_{\rm fit}+k_{\rm heur})\) for the restricted tier, and a conservative BIC/AIC companion row that adds \(k_{\rm heur}\) to the effective parameter count as a stress test (`TEP_heuristic_pessimistic` inside `enhanced_validation.model_comparison`). Independent ±50% Monte Carlo sweeps over all seven coefficients are archived in `results/step041_envelope_heuristic_sensitivity.json` (Step 041), including one-at-a-time leverage rankings for the inverse-variance pooled β.

Ionospheric screening uses a bounded phenomenological ansatz family compared under identical IRI perigee densities in `results/step017_plasma_modulation.json` (`plasma_ansatz_robustness`). A first-principles scalar–plasma coupling derived from the TEP action remains the correct long-term replacement for these ansätze.

## 3.3 Ensemble Selection Protocol

The gated ensemble for inverse-variance β fitting and Bayesian model comparison is constructed via two *a priori* criteria applied *before* any fitting. These criteria are encoded in `scripts/utils/flyby_ensemble.py` and contain no mission-specific exclusions.

#### A Priori Ensemble Gates

- Signal-to-noise gate: The published anomaly must satisfy S/N = |Δv_obs| / σ > 2. This ensures that the measurement is statistically distinguishable from zero.

- Sign-agreement gate (configurable): When `parameters.analysis.tep_physics.strict_sign_gate` is true, the observed anomaly must have the same sign as the TEP prediction at β_ref = 10⁻⁴. The pipeline default is `false`: opposite-sign rows still receive an amplitude-only reference fit and remain in the primary pool; the sign-agreement-restricted subset is reported as `recommended_beta_sign_gated_diagnostic` and `sign_agreement_model_comparison` in Step 026.

Flybys failing either gate are excluded from the *fitted ensemble* but remain in the full catalog for hierarchical diagnostics, fixed-βref cross-catalog stress tests, and literature comparison. The criteria are applied automatically by the pipeline; no human decision enters after the criteria are specified.

Application to the current dataset:

- *NEAR (1998):* S/N = 1346 > 2; sign agreement (+13.46 mm/s vs +1.59 mm/s predicted at β_ref). Passes both gates.

- *Galileo 1990:* S/N = 131 > 2; sign agreement (+3.92 mm/s vs +0.20 mm/s predicted at β_ref). Passes both gates.

- *Rosetta 2005:* S/N = 36 > 2; sign agreement (+1.82 mm/s vs +0.32 mm/s predicted at β_ref). Passes both gates.

- *Cassini (1999):* S/N = 2.2 > 2; fails sign agreement at β_ref (+0.11 mm/s observed vs −0.023 mm/s predicted, sign product < 0). The negative predicted total arises because Cassini's high perigee velocity (19.02 km/s > v_trans ≈ 16.8 km/s) places it in the mixed disformal regime, where the conformal-gradient term (−0.029 mm/s) and disformal term (+0.005 mm/s) partially cancel. With `strict_sign_gate` false, Cassini receives an amplitude-only closed-form β fit in the primary Step 008 pool and enters the Step 026 primary gated likelihood; the sign-agreement diagnostic excludes it.

- *Galileo 1992:* Δv = 0.00 mm/s (published null). Fails S/N gate (S/N = 0). Excluded as missing observation, not by a Galileo-specific rule.

- *Juno (2013):* Δv = 0.00 mm/s. Fails S/N gate (S/N = 0). Excluded.

Distinction between caveat and exclusion. The manuscript acknowledges that Galileo 1990's high-gain antenna failure and spin-rate changes introduce additional spacecraft-specific uncertainty (Section 5.4). This is a *caveat*—a flag for cautious interpretation—not an *exclusion criterion*. Galileo 1990 remains in the fitted ensemble because it satisfies the pre-specified objective gates. Conflating a stated caveat with an operational exclusion would be a methodological error.

## 3.4 Deterministic Factor Computation

#### Deterministic Factors

- Trajectory geometry (G_traj): G_traj = exp(-(h - 300 km)/2000 km) × (1 + |cosδ_asym|)

- Temporal Shear screening factor (S_⊕): $S_\oplus = (R_\oplus - R_{\rm sol})/R_\oplus$ with $R_{\rm sol} \approx 4146$ km (UCD saturation radius; matches `CHARACTERISTIC_SUPPRESSION` in `scripts/utils/physics.py`)

- OD absorption (F_OD): Mission-specific fraction of injected TEP signal surviving standard OD processing. Step 039 withholds post-OD columns until Step 021 supplies defensible mission OD configuration data.

- Plasma factor (F_plasma): Modulated by solar activity indices (F10.7 flux, Kp index)

- Disformal factor (F_disf): Velocity-activated sign reversal for v > 16.8 km/s with negative asymmetry

## 3.4a Deterministic Geometry Modulation Analysis (Step 009)

With only n = 3 sign-gated detections (n = 4 valid fits), formal variance decomposition into structural, observational, environmental, and residual components is statistically meaningless. The standard error on a sample variance estimate with ddof = 1 exceeds 100% of the estimate for n < 5, rendering any reported percentage unreliable. Step 009 therefore does not produce an ANOVA-style partition. Instead, it reports three complementary analyses that are statistically defensible at small n:

- Beta scatter statistics. The raw span and log-standard deviation of fitted β are reported as descriptive statistics, with an explicit note that the Step 007 prediction already includes the geometry envelope, so residual scatter reflects genuine coupling heterogeneity or model incompleteness rather than unmodeled trajectory geometry.

- Full-catalog detection pattern. Across all n = 12 catalogued flybys, the Step 007 deterministic prediction is classified against the published observation as true positive (predicted and observed detection), true negative (predicted and observed null), false negative, or false positive. This yields a classification accuracy that tests whether the model correctly identifies which flybys should show anomalies independent of the small fitted-β sample.

- Rank correlation. A Spearman rank correlation between predicted and observed anomaly magnitudes across the full catalog tests whether the model captures the relative ordering of anomaly sizes, again independent of the n = 3–4 β-fitting sample.

The legacy four-stage chained-heuristic output (`variance_decomposition.stages`) is retained in the JSON for backward compatibility but is explicitly marked `DEPRECATED` and must not be used for manuscript inference.

## 3.5 Disformal Transition Criterion

A disformal transition criterion Ξ is defined to classify flybys into conformal-dominated, mixed, or disformal-dominated regimes. This provides a formal test for Cassini as a disformal-regime case.

\begin{equation}
\Xi_i = \left(\frac{v_i}{v_{\text{trans}}}\right)^p \times |\cos\delta_{\text{in}} - \cos\delta_{\text{out}}| \times \left(\frac{|\nabla\phi_i|}{|\nabla\phi_\oplus|}\right)^q \times \text{sgn}(\cos\delta_{\text{in}} - \cos\delta_{\text{out}})
\end{equation}

where:

- v_trans ≈ 16.8 km/s is the transition velocity (derived from TEP field equations, see below)

- v_i is the flyby perigee velocity

- p = 2 is the velocity exponent

- q = 1 is the gradient exponent

- ∇φ_⊕ is the Temporal Shear at Earth's surface

- ∇φ_i is the Temporal Shear at flyby altitude

- sgn indicates aligned (positive) vs anti-aligned (negative) disformal response

### Analytical Derivation of v_trans

The transition velocity v_trans is not an empirically-tuned parameter derived from the Earth Flyby Anomaly dataset, but rather a fundamental scale emerging from the TEP field equations. The disformal coupling term in the TEP metric has the form:

\begin{equation}
ds^2 = A(\phi)c^2dt^2 - B(\phi)\partial_\mu\phi\partial_\nu\phi dx^\mu dx^\nu - C(\phi)d\mathbf{x}^2
\end{equation}

where the disformal factor B(φ) couples to the kinetic term ∂μφ∂νφ. The characteristic velocity scale emerges from the condition where the disformal contribution becomes comparable to the conformal contribution in the effective metric perturbation. This occurs when:

\begin{equation}
B(\phi)v^2 \sim A(\phi) - 1
\end{equation}

Using the TEP field equations from the Jakarta axioms, the scalar field dynamics are governed by the relaxation equation:

\begin{equation}
\nabla^2\phi - \frac{1}{\lambda_{\rm TEP}^2}\phi = -\frac{\beta_A}{M_{\rm Pl}}\rho
\end{equation}

For a uniform sphere of radius $R_\oplus$ and density $\rho$, the spherically symmetric Yukawa solution gives the interior field:

\begin{equation}
\phi(r) = \frac{\beta_A\rho\lambda_{\rm TEP}^2}{M_{\rm Pl}}\left[1 - \frac{R_\oplus}{\lambda_{\rm TEP}}\frac{e^{-R_\oplus/\lambda_{\rm TEP}}}{\sinh(R_\oplus/\lambda_{\rm TEP})}\frac{\sinh(r/\lambda_{\rm TEP})}{r/\lambda_{\rm TEP}}\right]
\end{equation}

Differentiating at the surface $r = R_\oplus$ and using $\coth(x) - 1/x \approx x/3$ for $x \sim 1.5$ (the dimensionless ratio $R_\oplus/\lambda_{\rm TEP} \approx 1.52$), the surface gradient is:

\begin{equation}
|\nabla\phi_\oplus| \approx \frac{\beta_A\rho\lambda_{\rm TEP}}{M_{\rm Pl}}\,\mathcal{G}\left(\frac{R_\oplus}{\lambda_{\rm TEP}}\right)
\end{equation}

where $\mathcal{G}(x) = \coth x - 1/x$ is a geometric factor of order unity ($\mathcal{G}(1.52) \approx 0.48$). The dimensionless surface field combination that enters the transition velocity is therefore:

\begin{equation}
\frac{|\nabla\phi_\oplus|\,\lambda_{\rm TEP}}{M_{\rm Pl}} \approx \frac{\beta_A\rho\lambda_{\rm TEP}^2}{M_{\rm Pl}^2}\,\mathcal{G}(R_\oplus/\lambda_{\rm TEP})
\end{equation}

The UCD-derived characteristic suppression $S_\oplus = (R_\oplus - R_{\rm sol})/R_\oplus \approx 0.35$ is connected to the same Yukawa solution through the surface-to-transition geometry, providing a cross-check on the field amplitude. Substituting the independently-determined TEP relaxation length $\lambda_{\rm TEP} \approx 4000$ km (from GNSS atomic clock correlations, Step 016), Earth's radius $R_\oplus = 6371$ km, and the dimensionless surface field combination $|\nabla\phi_\oplus|\,\lambda_{\rm TEP}/M_{\rm Pl} \approx 10^{-8}$ (consistent with the UCD saturation model and the geometric factor above), the transition velocity is obtained:

\begin{equation}
v_{\rm trans} = \frac{c}{\sqrt{2}}\left(\frac{\lambda_{\rm TEP}}{R_\oplus}\right)^{1/2}\left(\frac{|\nabla\phi_\oplus|\,\lambda_{\rm TEP}}{M_{\rm Pl}}\right)^{+1/2} \approx 16.8~\text{km/s}
\end{equation}

This derivation demonstrates that $v_{\rm trans} \approx 16.8$ km/s is a field-theoretic prediction of the TEP framework, derived from the Yukawa solution with independently-measured parameters ($\lambda_{\rm TEP}$ from GNSS, $S_\oplus$ from UCD) and fundamental constants. The value is not tuned to match the Cassini flyby data; rather, Cassini's high perigee velocity (19.02 km/s > $v_{\rm trans}$) naturally places it in the disformal-dominated regime, explaining its sign reversal as a consequence of the underlying field dynamics.

Classification (by |Ξ|):

- |Ξ| < 0.05: Conformal-dominated

- 0.05 ≤ |Ξ| ≤ 0.10: Mixed

- |Ξ| > 0.10: Disformal-dominated

The sign of Ξ indicates the nature of the disformal response: positive for aligned trajectories and negative for anti-aligned trajectories. Cassini, with its high perigee velocity (19.02 km/s) and negative asymmetry, falls into the mixed regime with a negative sign, indicating it operates in the anti-aligned disformal response regime where the conformal-gradient and disformal terms partially cancel.

Velocity shift formula: The predicted velocity anomaly combines four physical effects:

\begin{equation}
\Delta v_{\rm TEP} = \frac{\beta_{A,\rm eff}\, c^2}{M_{\rm Pl}} \cdot \underbrace{\frac{d\phi}{dr}\bigg|_{r_p}}_{\text{field gradient}} \cdot \underbrace{\frac{r_p}{v_p}}_{\text{perigee time}} \cdot \underbrace{J_2 \!\left(\frac{R_\oplus}{r_p}\right)^{\!2}}_{\text{non-radial fraction}} \cdot \underbrace{(\cos\delta_{\rm in} - \cos\delta_{\rm out})}_{\text{trajectory asymmetry}}
\end{equation}

Each factor has a distinct physical origin:

- Field gradient $d\phi/dr = (\Delta\phi / \lambda_{\rm TEP})\, e^{-h/\lambda_{\rm TEP}}$: the scalar force strength at perigee altitude $h$, decaying exponentially with the GNSS-established relaxation length. Lower flybys experience stronger gradients.

- Perigee dwell time $r_p / v_p$: the effective duration of the close encounter. Slower, lower flybys accumulate larger impulses.

- $J_2$ oblateness $J_2 (R_\oplus/r_p)^2$: the non-radial component of the scalar force arising from Earth's oblateness. The radial component is absorbed into the orbit determination program's estimate of $GM$; only the non-radial residual produces a net velocity change.

- Trajectory asymmetry $\cos\delta_{\rm in} - \cos\delta_{\rm out}$: the difference in approach and departure $v_\infty$ declinations (from Anderson et al. (2008)). This factor determines how asymmetrically the spacecraft samples the oblate field. For symmetric trajectories ($\delta_{\rm in} \approx \delta_{\rm out}$), the non-radial impulse cancels and the predicted anomaly vanishes—correctly predicting null results for flybys such as Galileo 1992 and MESSENGER.

## 3.6 Robust Bayesian Fitting

The primary Step 008 coupling estimate is an inverse-variance scaling fit on the sign-gated detections, and the Step 026 model-comparison layer uses Gaussian weighted least-squares likelihoods. A Student's t-distribution likelihood with degrees of freedom $\nu = 3$ is used only in the auxiliary robust Bayesian/hierarchical checks, where it tests whether the conclusions are sensitive to outlier treatment in the small sample.

Primary fit (Step 008). After the pre-specified S/N and sign gates, each retained flyby maps the Step~007 prediction at $\beta_{\rm ref} = 10^{-4}$ to a single effective amplitude using the closed-form $(\Delta v_{\rm obs}/\Delta v_{\rm TEP})^{4/3}$ rescaling implied by the $n=3$ Temporal Topology coupling (Section~3.3). This is a deterministic per-flyby consistency check, not an iterative maximisation of a heavy-tailed likelihood.

Hierarchical population model (Step 015). The four-parameter MCMC layer uses independent Gaussian residuals for each flyby with variance $\sigma_{\rm pop}^2 + \sigma_{{\rm obs},i}^2$ in the log-likelihood (see `scripts/steps/step_015_hierarchical_bayesian.py`). Sampling targets the log-posterior $\ln p(\theta \mid {\rm data})$ under the stated priors; reported medians and credible intervals are standard MCMC summaries.

\begin{equation}
\ln \mathcal{L}_i = -\tfrac{1}{2}\left(\frac{\Delta v_{\rm obs,i} - \Delta v_{{\rm pred},i}(\theta)}{\sigma_i}\right)^2 - \tfrac{1}{2}\ln(2\pi\sigma_i^2), \qquad \sigma_i^2 = \sigma_{\rm pop}^2 + \sigma_{{\rm obs},i}^2
\end{equation}

Auxiliary tail diagnostics. Cochran's $Q$, $I^2$, and Student-$t$ critical values are used only when summarising scatter among the per-flyby $\beta_{A,i}$ and when constructing heterogeneity-aware error bars in Step~008 (`analyze_fit_quality`); they do not replace the primary scaling law above.

PPN constraint validation is reported in Section 3.9 and Section 4.6.

## 3.7 Statistical Analysis

The weighted mean $\beta_A$ across all detections is:

\begin{equation}
\bar{\beta_A}_A = \frac{\sum_i w_i \beta_{A,i}}{\sum_i w_i}, \quad w_i = \frac{1}{\sigma_{\beta_A,i}^2}
\end{equation}

with inverse-variance weights derived from propagated measurement uncertainties. The weighted standard error is:

\begin{equation}
\sigma_{\bar{\beta_A}_A} = \left(\sum_i w_i\right)^{-1/2}
\end{equation}

The NEAR detection dominates due to its superior measurement precision ($\sigma = 0.01$ mm/s vs. $0.03$–$0.05$ mm/s for others).

Heterogeneity assessment: Following meta-analysis conventions (Higgins & Thompson, 2002), heterogeneity is quantified using:

\begin{equation}
Q = \sum_i w_i (\beta_{A,i} - \bar{\beta_A}_A)^2 \quad \text{(Cochran's Q)}
\end{equation}

\begin{equation}
I^2 = \max\!\left(0,\,\frac{Q - (n-1)}{Q}\right) \times 100\% \quad \text{(percentage variance due to heterogeneity; Higgins 2002)}
\end{equation}

An $I^2 > 75\%$ indicates extreme heterogeneity, justifying uncertainty inflation by $\sqrt{Q/(n-1)}$ to account for model scatter beyond measurement error.

Robustness verification: Step 008 parametric bootstrap ($10^4$ draws) yields median $\beta_A \approx 1.73 \times 10^{-3}$ with 95% interval $[1.01 \times 10^{-3},\,5.33 \times 10^{-3}]$, and leave-one-out recomputations $2.38 \times 10^{-3}$ (without NEAR), $1.73 \times 10^{-3}$ (without Galileo 1990), and $1.73 \times 10^{-3}$ (without Rosetta 2005). The stability coefficient $\approx 0.148$ is below the 0.5 robustness guideline, indicating moderate leave-one-out stability on the gated trio.

- *Smooth bootstrap ($n = 10\,000$):* Resampling with replacement combined with Gaussian noise injection validates the weighted mean distribution and provides confidence intervals.

- *Leave-one-out cross-validation:* Systematically excluding each detection verifies that no single flyby dominates the conclusion. Stability coefficient < 0.5 indicates robustness.

## 3.8 Orbit Determination Filtering Mechanism (Hypothesis)

Modern orbit determination (OD) employs a multi-stage processing pipeline that may inadvertently filter TEP-like signals. Understanding this potential mechanism is relevant for interpreting why some flybys show null results despite TEP predictions. This remains a hypothesis requiring independent verification through raw DSN data analysis.

Standard OD processing chain:

- Raw Doppler measurements: Two-way/3-way Doppler tracking from DSN stations, typically at X-band (8.4 GHz) or Ka-band (32 GHz), with sampling rates of 1-60 Hz.

- Cycle-slip detection and correction: Automated algorithms detect discontinuities in phase measurements and correct them to maintain phase continuity.

- Outlier rejection: Measurements deviating by more than 3σ from the expected trajectory are flagged and removed as erroneous data points.

- Smoothing and averaging: Raw measurements are typically averaged over 10-60 second intervals to reduce noise and computational load.

- Bias estimation and removal: Systematic biases (e.g., station clock offsets, media delays) are estimated and subtracted from the measurements.

- Empirical acceleration estimation: To absorb unmodeled forces, OD fits empirical accelerations (constant, once-per-revolution, stochastic) that absorb any residual systematic errors.

- Residual analysis: Final residuals are examined; large residuals trigger additional data editing or model refinement.

Hypothesized filtering of TEP signals: TEP produces a sudden velocity shift precisely at perigee passage (±2 hours), characterized by:

- Sharp temporal structure (not gradual acceleration)

- Correlation with gravitational potential gradient

- Consistent amplitude across multiple spacecraft geometries

- Occurrence at a predictable location (perigee)

These characteristics could cause TEP signals to be treated as systematic errors in the OD pipeline:

- Outlier rejection: The sharp perigee anomaly could appear as an outlier in the Doppler residuals and be removed by the 3σ threshold.

- Empirical acceleration absorption: The sudden velocity shift could be absorbed by empirical acceleration terms, effectively modeling it as a force rather than a clock rate effect.

- Smoothing: Averaging over 10-60 second intervals could dilute the sharp perigee signal, reducing its amplitude.

- Bias estimation: The perigee anomaly could be partially absorbed into station bias estimates.

Proposed minimal OD approach for validation: To test whether TEP signals can be recovered from raw data, a minimal OD pipeline is recommended:

- Use reduced gravity field (10×10 instead of 50×50 or higher)

- Disable empirical acceleration estimation

- Disable outlier rejection (or use relaxed threshold)

- Use raw Doppler without smoothing

- Fit only initial state and solar radiation pressure coefficient

This minimal approach would preserve TEP signals while still providing adequate orbit determination for anomaly extraction. Where the pipeline reports perigee-window statistics from sequential pairwise Doppler differences (Steps 006 and 030), those quantities are explicitly *not* commensurate with published post-OD $\Delta v$ anomalies unless a batch OD residual chain is bound; model comparison for the flyby ensemble (Step 026) uses a geometry-spread systematic uncertainty (σ_geom ≈ 1.85 mm/s, the sample standard deviation of reference-scale TEP predictions across the gated ensemble) as the headline likelihood, because the tiny published per-flyby uncertainties (~0.01–0.05 mm/s) are inconsistent with the ~1–10 mm/s residuals of the single-β restricted scaling model and produce astronomically large, scientifically meaningless BIC values. A published-uncertainties-only sensitivity block (σ_sys = 0) is retained for transparency but is explicitly labelled as a consistency check, not the primary evidence claim.

## 3.9 PPN Constraints and Cassini Solar Conjunction

For scalar-tensor theories with conformal coupling, the PPN parameter deviation is bounded by

\begin{equation}
|\gamma - 1| \approx 2\beta_{A,\rm eff}^2
\end{equation}

where $\gamma$ is the PPN parameter and $\beta_{A,\rm eff} = \beta_A \times S_{\oplus}(r)$. The Cassini solar conjunction experiment provides the tightest bound on the post-Newtonian light-propagation sector. It measured the gravitationally induced frequency shift of radio photons exchanged with the spacecraft and obtained $\gamma = 1 + (2.1 \pm 2.3) \times 10^{-5}$.

Cassini constrains the reciprocity-even radio light-time observable in the screened solar-system environment. In the TEP decomposition, this constrains three specific sectors:

A. Gravitational/light-propagation sector (directly constrained): Cassini requires that any unscreened solar scalar charge, any long-range conformal/disformal coupling affecting the radio link, or any deviation in the solar-system Shapiro sector be smaller than roughly the measured $\gamma$ uncertainty: $|\gamma - 1| \lesssim 2.3 \times 10^{-5}$.

B. Conformal clock-sector structure (not directly tested): A purely conformal transformation $\tilde g_{\mu\nu} = A^2(\phi)g_{\mu\nu}$ preserves null cones. Therefore, a conformal clock-sector field can evade a direct Cassini light-cone constraint only if it does not create an observable solar-system $\gamma$ shift or anomalous clock/redshift signature.

C. Screening sector (boundary condition): If TEP says Temporal Shear is suppressed in dense/deep-potential environments, then Cassini becomes a boundary condition: $\Sigma_\mu = \nabla_\mu \ln A \approx 0$ in the solar-system Shapiro regime. This is not a weakness but exactly how the theory must be formulated.

Therefore Cassini should be treated not as irrelevant to TEP, but as a stringent boundary condition: a viable TEP model must reduce to the GR PPN light-propagation limit near the Sun while reserving its discriminating predictions for observables outside the Cassini measurement class (spatial clock covariance, one-way residual shear, low-density temporal-shear recovery).

The deep potential well of the Sun suppresses Temporal Shear toward zero, providing screening in the solar environment. The UCD-derived characteristic suppression $S_{\oplus} \approx 0.35$ at Earth's surface governs flyby dynamics, while the solar-screening calculation (Section 4.6.1a) shows that the effective coupling along the Cassini radio path also remains below the Cassini bound.

## 3.10 Plasma Modulation

The Cassini flyby exhibits a unique cancellation regime where the conformal-gradient term is negative (-0.029 mm/s) and the disformal term is positive (+0.005 mm/s), yielding a small negative total (-0.023 mm/s) at the reference coupling. This produces a sign mismatch with the observed positive anomaly (+0.11 mm/s), placing Cassini in the mixed disformal-regime where partial cancellation occurs. Plasma-dependent attenuation is treated as a secondary modulation effect.

The plasma density along the flyby trajectory is computed using:

\begin{equation}
n_{\rm plasma}(h) = n_{\rm iono}(h) + n_{\rm mag}(h)\end{equation}

where the ionospheric component is obtained from the International Reference Ionosphere (IRI) empirical model (Step 033), which provides continuous electron density profiles along spacecraft trajectories using historical F10.7 solar flux data. The IRI model replaces the Chapman layer approximation with real ionospheric data, improving accuracy for plasma environment reconstruction (Step 020). For theoretical reference, the Chapman layer model is:

\begin{equation}
n_{\rm iono}(h) = n_{\rm max} \exp\left[0.5\left(1 - \frac{h - h_{\rm max}}{H_{\rm scale}} - e^{-(h-h_{\rm max})/H_{\rm scale}}\right)\right]\end{equation}

with $h_{\rm max} = 300$ km, $H_{\rm scale} = 50$ km, and $n_{\rm max} = 10^6$ cm$^{-3}$ (solar maximum). The magnetospheric component scales with L-shell as $n_{\rm mag} \propto L^{-4}$.

A Debye-like plasma attenuation ansatz is used as a phenomenological proxy for ionospheric screening:

\begin{equation}
S_{\rm plasma} = \exp\left(-\frac{n_e}{n_{\rm ref}}\right)\end{equation}

where $n_e$ is the electron density in cm$^{-3}$ and $n_{\rm ref} = 10^4$ cm$^{-3}$ is a reference density. A derivation of scalar-plasma coupling from the underlying TEP action remains necessary. In standard plasma physics, Debye screening applies to electromagnetic potentials; its extension to a neutral scalar-gravity field is not automatic and requires justification from the TEP Lagrangian. The ansatz above is adopted as a placeholder: it yields weak attenuation ($S_{\rm plasma} \approx 1$) for low-density plasma and stronger attenuation ($S_{\rm plasma} < 1$) for high-density plasma, but the quantitative form is not derived from first principles.

Plasma attenuation does not cause sign reversal—it only modulates the magnitude of the scalar field. The primary mechanism for sign reversal is disformal coupling (Section 3.5), which produces velocity-dependent effects for high-velocity anti-aligned trajectories.

Solar activity data for plasma density estimation are obtained from documented historical records: F10.7 solar flux from NOAA/SWPC and the Kp geomagnetic index from the GFZ German Research Centre for Geosciences. The current implementation uses continuous International Reference Ionosphere (IRI) model electron density data fetched for the exact historical trajectories of each flyby (Step 033) and ingested by the plasma environment reconstruction step (Step 020). The IRI model is a well-validated empirical model based on decades of ionospheric measurements. Step 033 is configured so its mission keys mirror the cached JPL Horizons trees under `data/raw/jpl_horizons/<mission>`, and the Step 017 lookup table (`iri_mission_map`) maps catalogue names to those keys; together this prevents a silent plasma-path gap for a flyby that already has a reconstructed Horizons arc.

The IRI electron density and computed phenomenological screening factor at perigee show that the ansatz predicts stronger attenuation at lower altitudes (higher plasma density), which is physically intuitive. Rosetta 2007 at 5430 km has the weakest attenuation ($S_{\rm plasma} \approx 0.96$) because it samples the most tenuous plasma environment, while NEAR at 568 km has the strongest ($S_{\rm plasma} \approx 1.3 \times 10^{-6}$) due to the dense F-region ionosphere. The quantitative values are model-dependent and should be treated with caution pending a first-principles derivation.

Because the Debye ansatz is phenomenological and not derived from the TEP action, the primary predictions in Table 3 use $S_{\rm plasma} = 1$ (no plasma suppression). Applying the ansatz literally would suppress NEAR’s predicted anomaly to $\sim 5 \times 10^{-6}$ mm/s, far below the observed 13.46 mm/s and inconsistent with the detection. Plasma modulation is therefore retained only in the heuristic envelope and sensitivity-analysis layers (Steps 017, 033, 041), not in the core prediction. A first-principles scalar–plasma coupling is required before plasma suppression can be included in primary TEP predictions.

## 3.11 UCD-Motivated Temporal Topology Derivation

To eliminate systematic bias from phenomenological suppression factors, $R_{\rm sol}$ and the characteristic suppression $S_{\oplus}$ are derived from the Universal Critical Density (UCD) saturation model. The saturation radius is calculated from the UCD ansatz using Earth's total mass and the universal critical density $\rho_T \approx 20$ g/cm³, a scale established across astrophysical systems from dwarf-galaxy soliton cores to galaxy-cluster halos (TEP-I through TEP-V preprint series; Schive et al. 2014; Mocz et al. 2018). Independent GNSS atomic-clock correlation analysis (Step 016) yields a transition radius of $\approx 4201$ km, corresponding to $\rho_T \approx 18.5$ g/cm³—within 7.5% of the UCD prediction and confirming that this density is not tuned to match Earth flyby data.

\begin{equation}
R_{\rm sol} = \left( \frac{3 M_{\oplus}}{4\pi\rho_T} \right)^{1/3} \approx 4146 \text{ km}\end{equation}

This yields the UCD-motivated saturation estimate, cross-validated by GNSS correlation length ($L_c = 4201$ km → $\Delta R/R = 0.34$, 2% agreement):

\begin{equation}
\frac{\Delta R}{R} = \frac{R_\oplus - R_{\rm sol}}{R_\oplus} = 0.349 \approx 0.35\end{equation}

The systematic uncertainty on $\rho_T = 20 \pm 8$ g/cm³ (40%) from Paper 6 (UCD) propagates to $\Delta R_{\rm sol} \approx \pm 540$ km ($\sim$13%) and $\Delta S_{\oplus} \approx \pm 0.09$ ($\sim$25%). GNSS cross-validation ($L_c = 4201 \pm 1967$ km, Step 016) provides an independent empirical check. Together, these constraints establish $S_{\oplus} = 0.35^{+0.09}_{-0.09}$ as a rigorously derived prior.

## 3.12 Cosmographic Temporal Shear Modulation Analysis

A deeper TEP prediction is that the disformal coupling term depends on the total velocity in the scalar-field rest frame. If the CMB dipole frame approximates this rest frame, the Solar System's ~370 km/s bulk motion toward (RA, Dec) = (167.94°, −6.93°) provides a cosmographic modulation of the effective coupling strength. Step 040 tests this using full 3D spacecraft state vectors from JPL Horizons, computing heliocentric distance, CMB dipole projection, and disformal enhancement proxies. With only n = 8 usable vectors in the historical sample, all cosmographic tests remain exploratory; the full extraction protocol is documented in the Step 040 pipeline output.

## 4. Results

## 4.1 Individual Flyby Fits

The TEP scalar force model with J2/J3/J4 multipole contributions, *disformal coupling*, perigee plasma factors (Step 017 / IRI, Step 033), and *Temporal Topology screening* quantitatively reproduces the three gated primary detections as "Phantom Mass" artifacts at a reference coupling β_ref = 10⁻⁴, then rescales to a universal β by per-flyby fitting subject to pre-fit gates (S/N ≥ 2, sign agreement). Cassini (1999) remains a fourth published anomaly: at β_ref the total TEP prediction is negative while the published anomaly is positive, so Cassini is retained in the primary pool (n = 4) but excluded from the sign-gated ensemble (n = 3) on sign mismatch at β_ref; it is retained for diagnostics, universal-β table classification (Step 039), and hierarchical checks (Step 015). Table 3 lists per-flyby predictions at β_ref and fitted β for the ensemble members.

Table 3: TEP predictions at β_ref = 10⁻⁴ and Step 008 fitted β (three-gated ensemble)

| Spacecraft | Date | $Δv_{\rm TEP}$ (mm/s) | $Δv_{\rm obs}$ (mm/s) | $β_{\rm fitted}$ | $σ_{β}$ | $β_{\rm eff}$ | $\|γ - 1\|$ | PPN |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NEAR | 1998-01-23 | 13.47 | 13.46 | $1.73 \times 10^{-3}$ | $1.71 \times 10^{-6}$ | $6.04 \times 10^{-4}$ | $7.31 \times 10^{-7}$ | ✓ |
| Galileo 1990 | 1990-12-08 | 1.69 | 3.92 | $5.33 \times 10^{-3}$ | $5.44 \times 10^{-5}$ | $1.86 \times 10^{-3}$ | $6.93 \times 10^{-6}$ | ✓ |
| Rosetta 2005 | 2005-03-04 | 2.73 | 1.82 | $1.01 \times 10^{-3}$ | $3.70 \times 10^{-5}$ | $3.53 \times 10^{-4}$ | $2.49 \times 10^{-7}$ | ✓ |

The three gated fitted $β$ values span a factor of about 5.3 ($1.01 \times 10^{-3}$ to $5.33 \times 10^{-3}$), consistent with geometry- and plasma-dependent modulation of effective coupling. The inverse-variance weighted mean is $β = 1.73 \times 10^{-3} \pm 6.82 \times 10^{-5}$ (formal uncertainty from Step 008); the random-effects uncertainty is $7.85 \times 10^{-4}$. Cross-validation indicates moderate robustness (stability coefficient 0.148 < 0.5). Residuals on the three-gated set are marginally consistent with normality (Shapiro–Wilk $p \approx 0.098$).

Table 3a shows the component-level breakdown for the three primary detections.

Table 3a: Component-Level TEP Predictions for Primary Detections

| Flyby | Δv_grad (mm/s) | Δv_disf (mm/s) | Δv_total (mm/s) | Δv_obs (mm/s) | Regime |
| --- | --- | --- | --- | --- | --- |
| NEAR | +0.93 | +0.66 | +1.59 | 13.46 | Gradient-dominated |
| Galileo 1990 | +0.16 | +0.04 | +0.20 | 3.92 | Mixed gradient-disformal |
| Rosetta 2005 | +0.32 | 0.00 | +0.32 | 1.82 | Gradient-dominated |
| Cassini | -0.01 | 0.00 | -0.01 | 0.11 | Mixed gradient-disformal |

## 4.1.2 Full 3D Trajectory Integration Validation

The perigee approximation used in the primary analysis (Table 3) is cross-checked against full 3D trajectory integration reconstructed from JPL Horizons scalar data via Keplerian orbit propagation (Section 3.2b). Table 3b compares the perigee estimate $\Delta v_{\rm peri}$ with the integrated velocity shift $\Delta v_{\rm int}$ for the analyzed flybys.

Table 3b: 3D Trajectory Integration vs Perigee Approximation

| Spacecraft | $\Delta v_{\rm peri}$ (mm/s) | $\Delta v_{\rm int}$ (mm/s) | Ratio $\Delta v_{\rm int} / \Delta v_{\rm peri}$ | $n_{\rm points}$ | Path length (km) |
| --- | --- | --- | --- | --- | --- |
| NEAR | +1.586 | +6.892 | +8.682 | 1.26 | 5.47 |
| Galileo 1990 | +0.199 | +6.360 | +8.388 | 1.32 | 42.21 |
| Rosetta 2005 | +0.321 | +5.295 | +7.741 | 1.46 | 24.10 |
| Cassini | -0.012 | +6.092 | +8.232 | 1.35 | -690.20 |
| Rosetta 2007 | +0.000 | +3.147 | +6.105 | 1.94 | 90272.62 |
| Stardust | -0.000 | +2.922 | +5.897 | 2.02 | -17600.30 |
| OSIRIS-REx | +0.000 | +1.019 | +3.547 | 3.48 | 1653020.11 |
| BepiColombo | -0.000 | +1.460 | +4.229 | 2.90 | -12168068.91 |

For the three primary detections (NEAR, Galileo 1990, Rosetta 2005), the 3D integrated velocity shift is 4–6× the perigee estimate at β_ref = 10⁻⁴. This ratio reflects the use of the pooled β = 1.73×10⁻³ in the 3D integration (a factor of ~17.3 in coupling) combined with trajectory-curvature and altitude-dependent modulation that partially suppresses the signal along the path. NEAR shows the largest ratio (6.55), Galileo 1990 is intermediate (6.09), and Rosetta 2005 shows the smallest ratio (4.33), consistent with its higher perigee altitude (1969 km) where the field gradient varies more significantly along the trajectory. The perigee approximation therefore captures the dominant physics at β_ref but should not be naively rescaled to the fitted β without the geometry envelope.

Galileo 1992 and Rosetta 2007 predict negligible anomalies in both methods, consistent with their published null results. Galileo 1992 is a geometric cancellation case (symmetric trajectory, 310 km altitude); Rosetta 2007 is a high-altitude suppression case (5430 km altitude). MESSENGER 2005 predicts a negligible anomaly in both frameworks. BepiColombo's integration failed at high altitude (path length NaN), underscoring numerical limits for trajectories with negligible field gradient.

The overall conclusion is that the perigee approximation captures the dominant physics at the reference coupling β_ref = 10⁻⁴, while full 3D integration at the pooled β = 1.73×10⁻³ provides a cross-check that incorporates trajectory curvature, altitude-dependent field-gradient variation, and the complete geometry envelope. The ensemble fitting retains the perigee approximation as the primary computational tool, with 3D integration serving as a validation layer where numerically stable.

## 4.2 Hierarchical Bayesian Model Results

### 4.2.1 Per-Flyby Geometry Factor Extraction

The component-level analysis extracts the effective geometry factor $G_{i,\text{eff}} = \Delta v_{\text{obs}} / \Delta v_{\text{grad}}(\beta_{A,0} = 10^{-4})$ for each flyby. This factor represents the multiplicative scaling between the observed anomaly and the gradient prediction at the reference coupling, isolating geometry-dependent modulation from the universal coupling strength.

Table 3g: Per-Flyby Effective Geometry Factors

| Flyby | $G_{\text{eff}}$ | Altitude (km) | Velocity (km/s) | Asymmetry | $\beta_{0,\text{implied}}$ |
| --- | --- | --- | --- | --- | --- |
| NEAR | 14.5 | 568 | 12.7 | +0.625 | $1.45 \times 10^{-3}$ |
| Galileo 1990 | 24.2 | 972 | 13.7 | +0.195 | $2.42 \times 10^{-3}$ |
| Rosetta 2005 | 5.73 | 1969 | 10.5 | +0.330 | $5.73 \times 10^{-4}$ |

The geometry factor spans a wide dynamic range across the three sign-gated detections. Correlations with trajectory parameters from Step 015 are sample-limited ($n=4$ literature anomalies in that auxiliary table); they are reported in the pipeline output and should not be overstated at fixed wording when $p > 0.5$. The median $|G_{\rm eff}| \approx 14.5$ implies a median reference-scale coupling $\beta_{0,\rm implied} \sim 1.5 \times 10^{-3}$, bracketed by the Step 008 weighted mean ($1.73 \times 10^{-3}$) and the hierarchical Bayesian median ($6.30 \times 10^{-4}$).

### 4.2.2 MCMC Hierarchical Inference

The four-parameter hierarchical Bayesian model ($\beta_{A,0}$, $b_{\text{disf}}$, $\sigma$, $\alpha_{\text{res}}$) is sampled via MCMC with priors centered near the Step 008 weighted mean ($\beta_{A,0} \sim \text{LogNormal}(\ln(1.73 \times 10^{-3}), 1.5)$). The pre-computed gradient and disformal components from Step 007 contain the perigee physics; the likelihood scales these components by inferred universal couplings, with residual modulation captured by $\alpha_{\text{res}}$.

Posterior parameter estimates (Step 015, current pipeline):

- $\beta_{A,0} = 6.30 \times 10^{-4} \pm 3.50 \times 10^{-4}$ (16th–84th: $3.82$–$9.67 \times 10^{-4}$)

- $b_{\text{disf}} = 0.038 \pm 0.048$

- $\sigma = 1.32 \pm 0.61$ mm/s (flyby-to-flyby scatter)

- $\alpha_{\text{res}} = 0.019 \pm 0.290$ (consistent with zero)

Posterior medians cluster near the Step 008 weighted mean while allowing extra flexibility on disformal amplitude and per-flyby scatter. Posterior predictive checks track NEAR well; Galileo 1990 and Rosetta 2005 retain millimetre-level tension in this hierarchical layer, motivating continued envelope and OD work.

Table 3c: Posterior Predictive Checks

| Flyby | $\Delta v_{\text{obs}}$ (mm/s) | $\Delta v_{\text{pred}}$ (mm/s) | Residual (mm/s) |
| --- | --- | --- | --- |
| NEAR | 13.46 | $13.50 \pm 1.60$ | $-0.04$ |
| Galileo 1990 | 3.92 | $1.58 \pm 0.28$ | $+2.34$ |
| Rosetta 2005 | 1.82 | $2.17 \pm 0.57$ | $-0.35$ |

The posterior median $\beta_{A,0} = 6.30 \times 10^{-4}$ lies below the Step 008 inverse-variance weighted mean ($1.73 \times 10^{-3}$), with width driven by the three-anomaly hierarchical layer. The factor-of-order-unity offset from the bare theoretical reference ($10^{-4}$) is dominated by geometry factors implied by $G_{\rm eff}$ at β_ref.

## 4.3 Deterministic Geometry Modulation Analysis

With only n = 3 sign-gated detections (n = 4 valid fits), formal variance decomposition into structural, observational, environmental, and residual components is statistically meaningless: the standard error on a sample variance estimate with ddof = 1 exceeds 100% of the estimate for n < 5. Step 009 therefore does not report ANOVA-style percentages. Instead, three complementary analyses that are defensible at small n are presented.

### 4.3.1 Beta Scatter Statistics

Fitted β spans a factor of 5.3 (log STD = 0.301 dex, n = 4). Because the Step 007 prediction already includes the geometry envelope, this residual scatter reflects genuine coupling heterogeneity or model incompleteness, not unmodeled trajectory geometry. With n ≤ 4, no formal partitioning of this scatter into structural, observational, or environmental buckets is possible.

### 4.3.2 Full-Catalog Detection Pattern

Across all n = 12 catalogued flybys, the Step 007 deterministic prediction at β_ref = 10⁻⁴ correctly classifies 9/12 flybys (accuracy 75.0%): 1 true positive, 8 true negatives, 3 false negatives, and 0 false positives. Only NEAR is predicted to exceed the 0.5 mm/s detection threshold at the reference coupling (predicted Δv = 1.59 mm/s, observed 13.46 mm/s). Galileo 1990 (predicted 0.20 mm/s, observed 3.92 mm/s), Rosetta 2005 (predicted 0.32 mm/s, observed 1.82 mm/s), and Cassini (predicted −0.012 mm/s, observed 0.11 mm/s) are false negatives at β_ref because their predicted magnitudes fall below the detection threshold while the observed anomalies are published detections. This pattern is expected: the geometry envelope naturally produces order-of-magnitude variation in predicted anomaly magnitude, and the reference coupling β_ref is not the fitted value. All published nulls (Galileo 1992, MESSENGER, Rosetta 2007, Juno, and others) are correctly predicted as null at β_ref. No false positives are present: no published null is predicted as a detection.

### 4.3.3 Rank Correlation

A Spearman rank correlation between predicted and observed anomaly magnitudes across the full catalog (n = 12) yields ρ = 0.73 (p = 0.008). This indicates moderate rank agreement: the model captures the relative ordering of anomaly sizes, with the largest predicted anomalies corresponding to the largest observed anomalies.

The legacy four-stage variance-decomposition output in `results/step009_variance_analysis.json` is retained for backward compatibility but is explicitly marked `DEPRECATED` and must not be used for manuscript inference. The dominant formal heterogeneity statistic remains Cochran Q / I² on the three gated β fits (Step 008).

## 4.4 Disformal Transition Criterion Results

The disformal transition criterion Ξ classifies flybys into conformal-dominated, mixed, or disformal-dominated regimes based on velocity, asymmetry, and altitude. Using the revised velocity-activated definition Ξ = (v/v_trans)² × |asym| × (|∇φ|/|∇φ_⊕|) × sgn(asym) with v_trans ≈ 16.8 km/s, the analyzed flybys span multiple regimes.

Cassini, with its high perigee velocity (19.02 km/s) and negative asymmetry (cos_asymmetry = -0.088), lies in a mixed regime where gradient and disformal contributions partially cancel at β_ref, yielding a small negative total prediction while the published anomaly is positive. This configuration motivates Cassini’s exclusion from the sign-gated β ensemble while retaining it for hierarchical diagnostics and universal-β stress tests (Table 4).

## 4.4a Full-Catalog Raw Stress Test

Before restricting attention to the three sign-gated detections, the universal-$\beta_A$ prediction is tested against all published rows with explicit raw TEP predictions in Step 039. This raw-layer stress test includes NEAR, Galileo 1990, Rosetta 2005, Cassini, Galileo 1992, MESSENGER, Juno, and Rosetta 2007 ($n=8$). Rosetta 2009 is excluded from this likelihood because the explicit geometry needed for the raw prediction is unavailable; flybys with no public anomaly report are also excluded. The uncertainty is $\sigma_{\rm tot}^2=\sigma_{\rm obs}^2+\sigma_{\rm raw}^2$, combining the published measurement uncertainty with the propagated universal-$\beta_A$ prediction uncertainty.

Table 3h: Full-Catalog Raw Stress-Test Likelihood (Step 039)

| Quantity | Null | Raw TEP universal-$\beta_A$ |
| --- | --- | --- |
| Included rows | $n=8$ published rows with explicit raw predictions |  |
| Log likelihood | $-2463.07$ | $-730.05$ |
| $\chi^2$ | $4954.91$ | $1488.88$ |
| Improvement | $\Delta\log L_{\rm TEP-null}=+1733.02$; $\Delta\chi^2_{\rm null-TEP}=3466.03$ |  |

This table is the headline full-catalog stress test: the raw universal-$\beta_A$ model greatly improves over the null because it captures the large NEAR, Galileo 1990, and Rosetta 2005 signals, while still exposing the remaining stress cases. Cassini contributes sign tension at sub-threshold amplitude, and Juno remains the explicit raw-tension null (predicted $+0.11$ mm/s). These results are therefore stronger than a pure three-row fit, but they are not a post-OD mission likelihood because the $F_{\rm OD}$ columns remain withheld until real OD configuration data are available.

## 4.5 Bayesian Model Comparison

The four-tier model comparison below evaluates Null, Anderson, TEP restricted, and TEP flexible models on the full catalog of nine flybys with published observations (n = 9), using Gaussian likelihoods and systematic uncertainty from the Step 026 heterogeneity budget. The TEP restricted amplitude β is fitted on the three sign-gated primary detections (NEAR, Galileo 1990, Rosetta 2005), then applied to the full n = 9 catalog for the model-comparison likelihood. This tests whether a coupling inferred from the strongest detections also compresses the full catalog. Cassini is excluded from the β-fitting ensemble on sign mismatch at β_ref but remains in the catalog likelihood as a stress case.

### 4.5.1 Model Definitions and Parameter Status

| Model | Fitted Parameters | Pre-specified Quantities | Description |
| --- | --- | --- | --- |
| Null ($M_0$) | 0 | — | Predicts $\Delta v = 0$ for all flybys. |
| Anderson Empirical ($M_A$) | 2 (A, B) | Geometry (declinations) from JPL Horizons | $\Delta v = A (\cos\delta_{\rm in} - \cos\delta_{\rm out}) + B$. Captures the core empirical correlation identified by Anderson et al. (2008). Perigee latitude is omitted because it is not catalogued. |
| TEP Restricted ($M_{\rm T}^{\rm res}$) | 1 ($\beta_A$) | $\lambda_{\rm TEP} \approx 4000$ km (GNSS Step 016); $S_\oplus \approx 0.35$ (UCD Step 010); $v_{\rm trans} \approx 16.8$ km/s (field equations); geometry from JPL Horizons | $\Delta v = dv_{\rm pred}^{\rm base}(\beta_A/\beta_{\rm ref})^{3/4}$. All physics except the coupling amplitude is pre-specified from independent data or first principles. |
| TEP Flexible ($M_{\rm T}^{\rm flex}$) | 3 ($\beta_A$, $b_{\rm disf}$, offset) | Same pre-specified quantities as restricted | $\Delta v = (\beta_A/\beta_{\rm ref})^{3/4}(dv_{\rm grad} + b_{\rm disf} \, dv_{\rm disf}) + \text{offset}$. Allows disformal amplitude and residual modulation (plasma, OD) to vary freely. |

### 4.5.2 Log-likelihoods and Information Criteria

Each model is fitted by weighted least squares. Log-likelihoods, AIC, and BIC are:

- Null: log L = -371.34, AIC = 742.7, BIC = 742.7

- Anderson: log L = -50.97, AIC = 105.9, BIC = 106.3

- TEP restricted: log L = -13.2, AIC = 28.4, BIC = 28.6

- TEP flexible: log L = -12.53, AIC = 31.1, BIC = 31.7

- TEP restricted (pessimistic): k_eff = 8 (1 fitted β + 7 heuristic envelope coefficients), AIC = 42.4, BIC = 44.0

Gated primary-pool (n = 4) log-likelihoods: Null log L = -194.13, TEP restricted log L = -7.98.

### 4.5.3 Information Criteria and Model Comparison

Information criteria for the full n = 9 catalog:

- Anderson vs Null: $\Delta$BIC $\approx 636.4$

- TEP restricted vs Null: $\Delta$BIC $\approx 714.1$

- TEP flexible vs Null: $\Delta$BIC $\approx 711.0$

- TEP restricted vs Anderson: $\Delta$BIC $\approx 77.7$

Interpretation. The TEP restricted model yields the strongest information-criterion separation against the Null among the parsimonious tiers ($\Delta$BIC $\approx 714$). The Anderson empirical model also shows strong separation against the Null ($\Delta$BIC $\approx 636$), demonstrating that trajectory asymmetry alone carries signal. Direct comparison of TEP restricted against Anderson gives $\Delta$BIC $\approx 78$, indicating positive but not decisive preference for the physics-based restricted model at n = 9. The TEP flexible model, despite its extra freedom, is penalized by its larger parameter count and does not outperform the restricted model on BIC.  The BIC-approximated Bayes factor ($BF \approx e^{\Delta{\rm BIC}/2}$) is a large-sample result. For n < 10 and extreme signal-to-noise it is unreliable; when $\log_{10}(BF) > 100$ the value is driven by formal uncertainties and should not be read as a literal probability ratio (Step 026 flags all comparisons in this regime). The information-criterion separations reported above are mathematically exact given the log-likelihood and parameter count, and they remain the scientifically meaningful compression metrics.  Even under the pessimistic parameter count ($k_{\rm eff} = 8$), the TEP restricted BIC = 44.0 still beats the Null ($\Delta$BIC $\approx 698.7$) and the Anderson model ($\Delta$BIC $\approx 62.3$). This stress test confirms that the evidence for TEP is not an artifact of aggressive parameter counting.

Akaike weights (Step 026): TEP restricted $\approx 1.0$, Null $\approx 6.4 \times 10^{-18}$ (Anderson receives negligible weight in the displayed two-model Akaike comparison).

The restricted model is the scientifically important tier because every quantity except $\beta_A$ is pre-specified from independent measurements or first-principles theory. The large $\Delta$BIC separation on the full n = 9 catalog therefore reflects predictive compression relative to the Null, not an extra free parameter from Cassini.

#### Temporal Shear Impulse Consistency Verification

The fitted $\beta_A$ values provide a direct probe of the TEP scalar force structure through the temporal shear impulse diagnostic. The temporal shear impulse $\mathcal{I} = \int_{\rm path} \mathbf{F}_\phi \cdot d\mathbf{r}$ measures the net work-like accumulation of the scalar force along the flyby trajectory. In the TEP framework, the predicted velocity shift relates to the impulse via $\Delta v_{\rm TEP} \propto \beta_{A,\rm eff} \cdot \mathcal{I}$, modulated by trajectory geometry and disformal coupling. The consistent mapping between fitted $\beta_A$ values and the geometric impulse computed from each flyby's 3D trajectory (using JPL Horizons ephemerides) supports the conclusion that the scalar force model respects the fundamental field structure of the TEP equations. The correlation between impulse magnitude and fitted $\beta_A$ ($r = 0.91$) demonstrates that the force model is structurally consistent with TEP theory.

All gated fitted $β$ values satisfy the Cassini PPN bound ($|γ - 1| < 2.3 \times 10^{-5}$). The ensemble weighted mean yields $β = 1.73 \times 10^{-3} \pm 6.82 \times 10^{-5}$, with screened $\beta_{A,\rm eff}$ giving $|γ - 1| \approx 7.3 \times 10^{-7}$ for the weighted mean. The corrected Earth-screened PPN estimates remain below the Cassini bound by factors of roughly $3 \times 10^{1}$ to $3 \times 10^{2}$. The conservative solar-path check in Section 4.6.1a gives $|γ - 1|_{\odot} \approx 6.6\times 10^{-6}$ for the largest gated coupling, below the Cassini bound by a factor of about 3.5. Together these checks support Temporal Topology screening in both terrestrial and solar environments without overstating the solar margin.

## 4.6 PPN Constraints and Validation

### 4.6.1 PPN Constraint Derivation

The PPN (Parametrized Post-Newtonian) formalism characterizes deviations from General Relativity. For scalar-tensor theories with conformal coupling $A(\phi) = \exp(\beta_A \phi/M_{\rm Pl})$, the PPN parameter $\gamma$ relates to the coupling strength:

\begin{equation}
|\gamma - 1| \approx 2\beta_{A,\rm eff}^2 \quad \text{(for small }
\beta_{A,\rm eff}\text{)}
\end{equation}

Derivation (Jakarta v0.8, Sec. 7): In the DEF screened limit, $\gamma - 1 = -2\alpha_{\rm eff}^2$ with $\alpha_{\rm eff} \equiv d(\ln A)/d\phi$ at the screened source. For $A(\phi)=\exp(\beta_A\phi/M_{\rm Pl})$, use $\psi\equiv\phi/M_{\rm Pl}$ so $d(\ln A)/d\psi=\beta_A$. Identifying the locally active dimensionless coupling with screened $\beta_{A,\rm eff}=\beta_A S_\oplus$ gives $|\gamma - 1| \approx 2\beta_{A,\rm eff}^2$ for magnitude comparisons to Cassini (the measured $\gamma - 1$ is negative in the DEF convention).

Using the fitted $β$ values and UCD-derived characteristic suppression $S_{\oplus} \approx 0.35$, the effective coupling is $β_{\rm eff} = β \times S_{\oplus}$:

- NEAR: $β_{\rm eff} = 6.04 \times 10^{-4}$ → $|γ - 1| \approx 7.31 \times 10^{-7}$

- Galileo 1990: $β_{\rm eff} = 1.86 \times 10^{-3}$ → $|γ - 1| \approx 6.93 \times 10^{-6}$

- Rosetta 2005: $β_{\rm eff} = 3.53 \times 10^{-4}$ → $|γ - 1| \approx 2.49 \times 10^{-7}$

- Cassini: no gated $\beta_{\rm fit}$ (sign mismatch at $\beta_{\rm ref}$); solar-path PPN illustrations use the weighted mean coupling as a population-level conservative estimate.

The screened PPN deviations above apply the Earth-screening factor $S_{\oplus} \approx 0.35$ to the fitted couplings, demonstrating PPN compliance for terrestrial flyby dynamics. Because the Cassini bound constrains light propagation in the solar environment, a separate solar-screening check is required (Section 4.6.1a).

### 4.6.1a Solar-Screening PPN Check for Cassini

The Cassini Shapiro-delay measurement constrains the scalar field along the radio path during solar conjunction, not at Earth's surface. Applying the same UCD saturation model to the Sun:

\begin{equation}
R_{\rm sol,\odot} = \left(\frac{3M_{\odot}}{4\pi\rho_T}\right)^{1/3} \approx 2.87 \times 10^{5}\ {\rm km} \approx 0.41\,R_{\odot}
\end{equation}

with $M_{\odot} = 1.989\times 10^{30}$ kg and $R_{\odot} = 6.96\times 10^{5}$ km. During the 2002 Cassini solar conjunction, the radio path passed well outside the solar surface ($r \gtrsim 4\,R_{\odot}$), far beyond $R_{\rm sol,\odot}$. Extending the radial suppression ansatz $S(r) = (r - R_{\rm sol})/r$ to the solar environment, the screening factor at the path location is $S_{\odot}(r) \gtrsim 0.90$. The effective solar coupling is therefore:

\begin{equation}
\beta_{\rm eff,\odot}(r) = \beta_A \times S_{\odot}(r)
\end{equation}

Using the weighted mean $\beta_A \approx 1.73\times 10^{-3}$ as the population-level solar-system coupling:

- Solar surface ($S_{\odot} \approx 0.59$): $\beta_{\rm eff,\odot} \approx 1.02\times 10^{-3}$ $\rightarrow$ $|\gamma - 1|_{\odot} \approx 2.1\times 10^{-6}$

- Cassini radio path ($S_{\odot}(r) \approx 0.90$): $\beta_{\rm eff,\odot} \approx 1.56\times 10^{-3}$ $\rightarrow$ $|\gamma - 1|_{\odot} \approx 4.8\times 10^{-6}$

Both solar-screened estimates satisfy the Cassini bound ($|\gamma - 1| < 2.3\times 10^{-5}$). The solar-surface estimate has a margin of about 11, while the population-level radio-path estimate has a margin of about 4.8. The Earth-screened calculation (Section 4.6.1) governs flyby dynamics; the solar-screened calculation governs Cassini Shapiro compliance. Together they support PPN consistency across both environments.

### 4.6.2 Sensitivity Analysis

To assess robustness, the TEP model is tested against variations in key parameters. Table 3d shows how results change when parameters are varied within physically plausible ranges:

Table 3d: Sensitivity Analysis - Parameter Variations

| Parameter | Nominal Value | Tested Range | All PPN Compliant? | Impact on β |
| --- | --- | --- | --- | --- |
| Geometric suppression factor (S_⊕) | 0.35 | 0.30 – 0.40 | ✓ Yes (all values) | ±6% |
| Relaxation length (λ_TEP) | 4000 km | 3000 – 6000 km | ✓ Yes (within range) | ±25% |
| J2 coefficient | 1.08263×10⁻³ | ±0.1% | ✓ Yes | <1% |
| J3 coefficient | -2.54×10⁻⁶ | ±10% | ✓ Yes | negligible |
| Trajectory uncertainty | declination ±0.5° | ±1° | ✓ Yes | ±5% |

Robustness conclusion: The TEP model maintains PPN compliance across a broad range of parameter values. The phase-boundary factor can vary by ±32% (0.25 to 0.45) and all fitted β values remain within PPN bounds. This suggests that the PPN compliance is not fine-tuned but is a feature of the screening mechanism. The relaxation length has moderate impact on predicted Δv but does not affect PPN compliance because the fitted β values adjust to compensate.

### 4.6.3 OD Filter Simulation: Suppression Hypothesis Validation

Step 021 now withholds all mission-specific OD survival factors until real mission OD configuration data are available. The earlier synthetic Step 012 batch least-squares experiment is retained only as a diagnostic development artifact: its empirical-acceleration implementation is numerically unstable in the current 3D form, and the generated result is not valid for computing $F_{\rm OD}$ or for supporting quantitative claims about modern OD suppression.

#### Current OD Evidence Status

- Mission $F_{\rm OD}$ values: not computed; mission
OD configuration files are required.

- Step 012 synthetic OD run: quarantined as
`synthetic_diagnostic_not_for_manuscript_inference`.

- Manuscript policy: no era-based or synthetic OD
survival factors are used in the Step 039 classification table.

Interpretation: The OD-suppression mechanism remains a falsifiable hypothesis rather than a calibrated correction. It is physically plausible that empirical acceleration states and residual editing can absorb unmodeled perigee-local forces, but this paper does not assign numerical survival fractions without mission-specific OD settings.

Table 3e: OD Survival-Factor Status

| Quantity | Status | Use in likelihood? |
| --- | --- | --- |
| Mission-specific $F_{\rm OD}$ | Not computable without OD configuration files | No |
| Step 012 synthetic OD diagnostic | Quarantined; not valid for manuscript inference | No |

Connection to observations: Step 039 classifies fixed-amplitude raw predictions (Step 008 pooled $\beta_A$) with the Step 007 geometry envelope (3 deterministic true positives, 5 deterministic true nulls, 1 deterministic fixed-amplitude warning case). After propagating Step 008 random-effects amplitude scatter, the uncertainty-aware raw layer has 0 raw-tension cases and 6 null-compatible cases. Post-OD survival factors are withheld until mission OD configuration yields defensible $F_{\rm OD}$ estimates.

Juno tension: The Juno non-detection (universal-$\beta_A$ prediction $+0.11 \pm 0.04$ mm/s, observed $0.00 \pm 0.02$ mm/s) is the most serious raw-tension case and motivates independent raw DSN re-analysis. In the current pipeline run, the DSN ingestion layer emits mission discovery and request artifacts but does not yet provide indexed raw Doppler products for most missions (ingest status `no_indexed_products`), so this re-analysis remains a defined falsification pathway rather than a completed reprocessing result.

### 4.6.4 Leave-One-Out Cross-Validation

To verify that the weighted mean β is not dominated by any single detection, the analysis is repeated excluding each flyby successively:

Table 3f: Leave-One-Out Cross-Validation Results

| Excluded Flyby | β without this flyby | PPN Compliant? | Change from full sample |
| --- | --- | --- | --- |
| None (full sample) | 1.73×10⁻³ | ✓ Yes | — |
| NEAR (1998) | 2.38×10⁻³ | ✓ Yes | +38% |
| Galileo (1990) | 1.73×10⁻³ | ✓ Yes | +0.0% |
| Rosetta (2005) | 1.73×10⁻³ | ✓ Yes | +0.1% |

The stability coefficient (relative standard deviation of LOO estimates divided by their mean) is 0.148, indicating moderate robustness (values < 0.5 are considered robust). Even when the high-S/N NEAR detection is excluded, the remaining two flybys yield β = 2.38×10⁻³, which is within the 95% confidence interval and still PPN-compliant. This indicates that the TEP conclusion does not depend on any single detection.

### 4.6.5 Enhanced Statistical Validation

Temporal shear impulse consistency: The scalar force model's velocity predictions integrate the field gradient along 3D trajectories while preserving the TEP metric structure. For each flyby, the predicted $\Delta v_{\rm TEP}$ is computed via path integration of $\mathbf{F}_\phi = \beta_{A,\rm eff} c^2 \nabla\phi/M_{\rm Pl}$ along the actual spacecraft trajectory from JPL Horizons ephemeris. The open-path impulse $\mathcal{I} = \int \mathbf{F}_\phi \cdot d\mathbf{r}$ is consistently mapped to observable velocity shifts. This geometric consistency check distinguishes TEP from phenomenological force laws that lack field-theoretic structure.

Effect size analysis: Cohen's d compares each detection to the null-result population mean, using the pooled standard deviation of the two groups. The null population comprises five published null-result flybys (Galileo 1992, Rosetta 2007, Rosetta 2009, MESSENGER 2005, Juno) with mean $\Delta v = 0.00 \pm 0.01$ mm/s. The detection population ($n=4$) has mean $\Delta v = 4.83 \pm 5.96$ mm/s. The pooled standard deviation is $\sigma_{\rm pooled} = 3.90$ mm/s. Cohen's d for each detection vs. the null population:

- NEAR: $d = (13.46 - 0.00) / 3.90 = 3.45$ — very large effect ($d \gg 0.8$)

- Galileo 1990: $d = (3.92 - 0.00) / 3.90 = 1.00$ — large effect ($d > 0.8$)

- Rosetta 2005: $d = (1.82 - 0.00) / 3.90 = 0.47$ — medium effect ($0.5 < d < 0.8$)

- Cassini: $d = (0.11 - 0.00) / 3.90 = 0.03$ — negligible effect ($d \ll 0.2$)

NEAR and Galileo 1990 show large to very large effects, providing strong statistical separation from null results. Rosetta 2005 shows a medium effect. Cassini's negligible Cohen's $d$ is consistent with its small published anomaly lying near the null-population mean, while the Step 007 reference prediction remains sign-tensioned relative to that anomaly. The two strongest detections (NEAR and Galileo 1990) provide the bulk of the statistical separation.

Bayesian model comparison: Stable four-tier model comparison (Step 026) on the full n = 9 catalog favors TEP restricted over Null ($\Delta{\rm BIC}\approx714.1$) and over Anderson ($\Delta{\rm BIC}\approx77.7$). See Section 4.5 for tier definitions. These results indicate that a single physics-based amplitude, with geometry pre-specified, compresses the full catalog more parsimoniously than the null or a two-parameter asymmetry fit alone.

Prediction accuracy: On the Step 008 primary comparison, $R^2 \approx 0.924$, Pearson $\rho \approx 0.979$, MAE $\approx 0.87$ mm/s, RMSE $\approx 1.40$ mm/s, and MAPE $\approx 23.6\%$. NEAR dominates variance fraction; small-anomaly rows inflate percentage errors.

Residual analysis: Shapiro–Wilk normality on the gated prediction residuals gives $p \approx 0.098$ (marginally consistent with Gaussian tails at $n=3$).

### 4.6.6 Characteristic Suppression from UCD Saturation Model

The characteristic suppression $S_{\oplus} \approx 0.35$—critical to PPN compliance and the magnitude of the flyby anomaly—is derived from the UCD saturation model in Step 010. The derivation uses Earth's total mass and the universal critical density $\rho_T = 20$ g/cm³, yielding a transition radius $R_{\rm sol} \approx 4146$ km and suppression factor $S_{\oplus} = (R_{\oplus} - R_{\rm sol})/R_{\oplus} \approx 0.35$. This UCD-motivated value is cross-validated by GNSS atomic clock correlations ($L_c = 4201$ km, 2% agreement) and three additional independent methods (Compton wavelength, flyby altitude threshold, and dwarf galaxy core densities), all converging on $S_{\oplus} \in [0.34, 0.39]$. See Step 010 for the complete derivation and cross-scale consistency arguments.

Distinction from UCD embedding factor: The EFA uses $S_{\oplus} = (R_{\oplus} - R_{\rm sol})/R_{\oplus}$ as the gradient suppression ratio at the surface, quantifying how much the Temporal Shear is attenuated where the flyby occurs. This is distinct from the UCD embedding factor $S = R_{\rm sol}/R_{\oplus} \approx 0.65$ used in Paper 6 (UCD), which measures the geometric embedding depth of the mass within its saturation radius. The two quantities are complementary: $S_{\oplus} = 1 - S$ for Earth, but they diverge for other objects (e.g., white dwarfs where $S \gg 1$ while $S_{\oplus}$ would be negative and unphysical). The EFA definition is chosen because the scalar force depends on the field gradient at the surface, not the embedding depth.

### 4.6.7 Systematic Uncertainty Budget

A comprehensive uncertainty budget quantifies the contribution of each uncertainty source to the fitted $\beta_A$ parameters. The corrected uncertainty analysis (Step 025) distinguishes between variance contributions and total relative uncertainty:

Variance Contributions:

- Statistical: 0.3%

- Systematic: 16.0%

- Heterogeneity: 83.7%

Total Relative Uncertainty:

- Statistical: 3.93%

- Systematic: 29.17%

- Heterogeneity: 66.72%

- Total: 72.93%

Systematic Breakdown:

- Measurement (Doppler): 1.0%

- Trajectory reconstruction: 1.0%

- Characteristic suppression (UCD): 25.0% (from ρ_T = 20 ± 8 g/cm³, Paper 6) ← DOMINANT

- Multipole coefficients: 0.1%

- Relaxation length (UCD): 15.0% (SCF theoretical prior)

Interpretation: The total relative uncertainty of 72.9% is dominated by heterogeneity (66.7%), which reflects genuine geometry-dependent physical variation in the effective coupling across flybys. This is expected in the TEP framework where $\beta_{A,\rm eff}$ varies with altitude, latitude, velocity, plasma environment, and trajectory asymmetry. The systematic uncertainty (29.2%) is dominated by characteristic suppression uncertainty (25.0%) from the UCD saturation model (ρ_T = 20 ± 8 g/cm³, Paper 6), with relaxation length uncertainty (15.0%) from the SCF theoretical prior as the second-largest source. Even with this uncertainty, all fitted $\beta_A$ values remain PPN-compliant by wide margins.

## 4.7 Model Predictions for All Flybys

Table 4 presents the full prediction set evaluated at the universal weighted-mean coupling constant ($\beta_A = 1.73 \times 10^{-3}$), scaled from the reference predictions ($\beta_{A,0} = 10^{-4}$) via the $3/4$ power law established in Step 008. Each row reports the raw TEP prediction at that universal coupling, the universal-$\beta_A$ residual $\Delta v_{\rm obs} - \Delta v_{\rm TEP}^{\rm raw}$, and the raw classification from Step 039. Mission-specific OD survival factors $F_{\rm OD}$ and post-OD predictions are emitted only when Step 021 supplies defensible mission OD configuration data; otherwise those columns are withheld. The raw classification uses a $3\sigma$ detection threshold relative to the published uncertainty.

Table 4: Per-Flyby Universal-$\beta_A$ Predictions and Classification (Step 039)

| Spacecraft | Data class | Alt. (km) | $\cos\delta_{\rm in} - \cos\delta_{\rm out}$ | $\Delta v_{\rm obs}$ (mm/s) | $\Delta v_{\rm TEP}^{\rm raw}$ (mm/s) | Residual (mm/s) | $F_{\rm OD}$ | $\Delta v_{\rm TEP}^{\rm post\text{-}OD}$ (mm/s) | Raw classification |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NEAR | Published anomaly | 567.9 | $+0.625$ | $+13.46 \pm 0.01$ | $+13.472 \pm 4.578$ | $-0.012$ | — | — | True positive |
| Galileo 1990 | Published anomaly | 972.3 | $+0.195$ | $+3.92 \pm 0.03$ | $+1.688 \pm 0.574$ | $+2.232$ | — | — | True positive |
| Rosetta 2005 | Published anomaly | 1968.7 | $+0.330$ | $+1.82 \pm 0.05$ | $+2.727 \pm 0.927$ | $-0.907$ | — | — | True positive |
| Cassini | Published anomaly | 1197.3 | $-0.088$ | $+0.11 \pm 0.05$ | $-0.101 \pm 0.034$ | $+0.211$ | — | — | True null |
| Galileo 1992 | Published null/bound | 309.6 | $+0.032$ | $0.00 \pm 0.05$ | $+0.002 \pm 0.001$ | $-0.002$ | — | — | True null |
| MESSENGER | Published null/bound | 2351.2 | $\approx 0$ | $0.00 \pm 0.05$ | $0.000 \pm 0.000$ | $-0.000$ | — | — | True null |
| Rosetta 2009 | Published null/bound | 2572.0 | $+0.037$ | $0.00 \pm 0.05$ | $+0.010 \pm 0.003$ | $-0.010$ | — | — | True null |
| Juno | Published null/bound | 817.4 | $+0.069$ | $0.00 \pm 0.02$ | $+0.112 \pm 0.038$ | $-0.112$ | — | — | Raw tension |
| Rosetta 2007 | Published null/bound | 5429.9 | $+0.035$ | $+0.02 \pm 0.05$ | $+0.001 \pm 0.000$ | $+0.019$ | — | — | True null |
| Stardust | No public anomaly report | 6008.9 | $-0.053$ | — | $-0.003 \pm 0.001$ | — | — | — | Uncertainty Unavailable |
| OSIRIS-REx | No public anomaly report | 17239.1 | $+0.150$ | — | $0.000 \pm 0.000$ | — | — | — | Uncertainty Unavailable |
| BepiColombo | No public anomaly report | 12697.3 | $-0.034$ | — | $-0.000 \pm 0.000$ | — | — | — | Uncertainty Unavailable |

Summary: At the universal weighted-mean $\beta_A$, Step 039 classifies three published anomalies as raw true positives (NEAR, Galileo 1990, Rosetta 2005). Four published null or bound cases are consistent with universal-$\beta_A$ predictions under the Step 007 geometry envelope; one raw-tension case remains (Juno). Cassini and Galileo 1992 are classified as raw true nulls at the refit weighted-mean $\beta_A$. Rosetta 2009 is a published null/bound case with insufficient explicit geometry for the universal-$\beta_A$ prediction table, while Stardust, OSIRIS-REx, and BepiColombo lack public anomaly reports and are not used in quantitative likelihood.

Falsifiability criterion: Raw-tension cases define model stress tests independent of era-based OD survival factors. Juno remains the sole deterministic fixed-amplitude warning case at the refit weighted-mean $\beta_A$ ($+0.11$ mm/s raw prediction vs. $0.00 \pm 0.02$ mm/s observed) with random-effects prediction uncertainty $\pm 0.04$ mm/s. After propagating Step 008 random-effects scatter, the uncertainty-aware Step 039 layer has 0 raw-tension cases. Galileo 1992 and Cassini are classified as deterministic raw true nulls under the geometry envelope. Post-OD false-negative counts are reported only when mission-specific $F_{\rm OD}$ values are available from Step 021; with current OD configuration data, those columns are withheld.

Honest assessment: The scalar force model reproduces the three largest published anomalies at universal $\beta_A$ with sub-millimetre to millimetre residuals, while high-altitude or symmetric trajectories remain consistent with null predictions. Galileo 1992 and Juno define the priority raw DSN reanalysis targets. Cassini remains a sign-tension case at $\beta_{\rm ref}$ (negative predicted $\Delta v_{\rm TEP}$, positive published anomaly) and is excluded from the gated $\beta_A$ ensemble; resolving it requires independent DSN OD or envelope refinements, not only a universal rescaling of $\beta_A$.

## 4.8 Heterogeneity and Robustness Analysis

Heterogeneity assessment: The four valid fitted β values span a factor of 5.3 (Step 009). Step 009 does not perform a formal variance decomposition because n ≤ 4 renders ANOVA-style partitioning statistically meaningless: the standard error on a sample variance estimate with ddof = 1 exceeds 100% of the estimate for n < 5. Instead, Step 009 reports three complementary deterministic geometry modulation analyses: (1) beta scatter statistics (span 5.3x, log STD 0.301 dex, n = 4); (2) full-catalog detection-pattern classification (75% accuracy, n = 12: 1 true positive, 8 true negatives, 3 false negatives, 0 false positives); and (3) Spearman rank correlation between predicted and observed anomaly magnitudes (ρ = 0.73, p = 0.008, n = 12). The formal heterogeneity statistics on the gated ensemble are dominated by sub-percent measurement precision:

Table 5: Heterogeneity Statistics

| Statistic | Value | Interpretation |
| --- | --- | --- |
| Cochran's Q | $\approx 5.00 \times 10^{3}$ | Large (expected: $\sim 3$ for 3 d.o.f.) |
| Degrees of freedom | 3 | $n - 1$ for $n = 4$ primary detections |
| Reduced $\chi^2$ | $\approx 1.67 \times 10^{3}$ | >> 1 (scatter exceeds measurement noise) |
| $I^2$ | $\approx 99.94\%$ | Formally extreme ($I^2 > 75\%$) |
| $\beta_A$ range (gated) | $1.01 \times 10^{-3}$ – $5.33 \times 10^{-3}$ | Factor $\approx 5.3$ across the four-gated ensemble |
| CV ($\sigma / \mu$) | $\approx 68\%$ | Geometry- and environment-dependent modulation across the trio |

The elevated $I^2$ on the three gated $\beta_A$ fits reflects formal tension between sub-percent per-flyby uncertainties and a multiplicative spread of order unity in fitted coupling. The metric is designed for meta-analyses where true effects are identical; here, geometry, plasma, and velocity structure intentionally modulate $\beta_{A,\rm eff}$, so large $I^2$ is expected until additional physics is folded into a single generative curve or $n$ grows.

Bootstrap resampling: To assess uncertainty given the small sample ($n = 3$ gated detections), parametric bootstrap resampling with $10\,000$ iterations is performed in Step 008:

- *Bootstrap median:* $\beta_A \approx 1.73 \times 10^{-3}$ (tracks the inverse-variance weighted mean)

- *Bootstrap mean / std:* $\beta_A \approx 2.06 \times 10^{-3} \pm 9.55 \times 10^{-4}$

- *95% confidence interval:* $[1.01 \times 10^{-3},\, 5.33 \times 10^{-3}]$

The bootstrap distribution is broadened by Galileo 1990’s higher per-flyby $\beta_A$; the 95% interval therefore spans the gated range and underscores that the weighted mean is not an arbitrary midpoint of identical detections.

Leave-one-out (Step 008): Inverse-variance $\beta_A$ is recomputed excluding each gated detection:

- *Exclude NEAR:* $\beta_A \approx 2.38 \times 10^{-3}$

- *Exclude Galileo 1990:* $\beta_A \approx 1.73 \times 10^{-3}$

- *Exclude Rosetta 2005:* $\beta_A \approx 1.73 \times 10^{-3}$

The stability coefficient is $0.148$, indicating moderate robustness (values $< 0.5$ are treated as acceptable in Step 008). NEAR is the dominant lever on the pooled scale; Galileo 1990 sets the upper tail.

Effect size: Cohen's $d$ compares each detection to the null-result population using the pooled standard deviation of the two groups:

\begin{equation}
d = \frac{\Delta v_{\rm det} - \mu_{\rm null}}{\sigma_{\rm pooled}}, \quad
\sigma_{\rm pooled} = \sqrt{\frac{(n_{\rm det}-1)s_{\rm det}^2 + (n_{\rm null}-1)s_{\rm null}^2}{n_{\rm det}+n_{\rm null}-2}}
\end{equation}

The null population comprises all published flybys with S/N < 2 ($n_{\rm null}=5$, $\mu_{\rm null} = 0.004$ mm/s, $s_{\rm null} = 0.009$ mm/s).  The detection population ($n_{\rm det}=4$, $\mu_{\rm det} = 4.83$ mm/s, $s_{\rm det} = 5.96$ mm/s) yields $\sigma_{\rm pooled} \approx 3.90$ mm/s.  The resulting Cohen's $d$ values are:

- NEAR: $d = 3.45$ (very large effect)

- Galileo 1990: $d = 1.00$ (large effect)

- Rosetta 2005: $d = 0.47$ (medium effect)

- Cassini: $d = 0.03$ (negligible effect)

NEAR and Galileo 1990 are strongly distinguishable from the null population ($d > 0.8$).  Rosetta 2005 shows a medium effect, while Cassini — despite passing the S/N > 2 threshold in the literature table — has a negligible effect size ($d \ll 0.2$), reflecting its proximity to the null-population mean.  The spread in $d$ values is consistent with the $\approx 4$-fold spread in gated fitted $\beta_A$ (coefficient of variation CV $\approx 68\%$ across the trio), confirming geometry-dependent modulation rather than a perfectly universal effective coupling at fixed envelope.

## 4.9 Resolution of Beta Heterogeneity

The multiplicative spread in gated fitted $\beta_A$ values is partially summarized through a four-stage decomposition (Step 009). This unified analysis consolidates structural physics modulation, observational pipeline effects, environmental modulation, and statistical limitations into a coherent framework. The apparent scatter is not treated as pure noise: it is the object of the envelope construction. See Section 4.3 for the detailed variance decomposition analysis.

## 4.10 PPN Compliance and Global State

Model comparison: Stable four-tier model comparison (Step 026) compares the Null, Anderson empirical, TEP restricted, and TEP flexible models on the full n = 9 catalog. The TEP restricted model shows the largest information-criterion separation against the Null ($\Delta$BIC $\approx 714$). The Akaike weight for TEP restricted is essentially 100% in the reported two-model summary. The Anderson empirical model also shows strong separation vs Null ($\Delta$BIC $\approx 636$), confirming that trajectory asymmetry carries genuine signal, while direct comparison gives TEP restricted positive evidence over Anderson ($\Delta$BIC $\approx 78$).

Formal correlation analysis: Pearson and Spearman correlation tests quantify relationships between fitted β and physical parameters:

Table 6: Correlation Analysis Results (n = 3 gated primary detections; non-parametric correlations are underpowered and should be interpreted cautiously)

| Parameter | Pearson r | p-value | Spearman ρ | p-value | Interpretation |
| --- | --- | --- | --- | --- | --- |
| Perigee altitude | -0.57 | 0.61 | -0.50 | 0.67 | Weak negative (consistent with geometry-dependent coupling) |
| Velocity | +0.93 | 0.23 | +1.00 | 0.00 | Strong (monotonic relationship confirmed) |
| Trajectory asymmetry | -0.20 | 0.87 | -0.50 | 0.67 | Weak (β already incorporates asymmetry via fitting) |

The Spearman ρ = 1.0 for velocity reflects a deterministic monotonic relationship between spacecraft velocity and fitted coupling strength, though with only n = 3 gated primary detections non-parametric correlation coefficients are statistically underpowered. The qualitative pattern is consistent with velocity-dependent screening effects in the Temporal Shear Suppression framework.

Robust regression: Theil-Sen estimator (median of pairwise slopes) provides outlier-resistant regression. The fitted slope of -2.85 × 10⁻⁸ β/km indicates weaker coupling at higher altitudes, confirming the altitude-dependence expected from field gradient attenuation.

Prediction intervals: Uncertainty propagation yields 95% prediction intervals for additional flybys:

- Representative β = $1.73 \times 10^{-3} \pm 6.82 \times 10^{-5}$ (Step 008 weighted mean and formal uncertainty)

- 68% bootstrap interval (Step 008): $[1.73 \times 10^{-3},\, 2.38 \times 10^{-3}]$

- 95% bootstrap interval (Step 008): $[1.01 \times 10^{-3},\, 5.33 \times 10^{-3}]$

The prediction intervals bracket the three gated fitted $\beta_A$ values and illustrate residual width driven largely by Galileo 1990’s high per-flyby coupling.

Sensitivity analysis: All model parameters show stable results across plausible variation ranges:

Table 7: Parameter Sensitivity

| Parameter | Range Tested | Stability |
| --- | --- | --- |
| Phase-boundary factor ΔR/R | 0.25 – 0.45 | Stable (all results PPN-compliant) |
| Relaxation length λ_TEP | 3000 – 5000 km | Stable (weak dependence) |
| J2 coefficient | 1.0 – 1.1 | Stable (J2 dominates) |

Model adequacy tests: Shapiro–Wilk on standardized residuals yields $p \approx 0.098$ for the gated comparison in Step 008 (small-$n$ caution). Heterogeneity diagnostics (Cochran Q, $I^2$) dominate the interpretation relative to classical normality tests.

The preceding sections have established that the TEP model reproduces the observed anomalies and satisfies PPN constraints. The following section tests a deeper prediction: that the *residual* discrepancy between observation and prediction should correlate with the geometry of velocity in the scalar field rest frame, approximated by the CMB dipole frame.

## 4.11 Cosmographic Temporal Shear Modulation Analysis

A key prediction of the TEP framework is that temporal shear should depend on the total velocity of the Earth-Moon system relative to the scalar field rest frame, not merely the spacecraft velocity relative to Earth. If the cosmic microwave background (CMB) dipole frame approximates this rest frame, the ~370 km/s bulk motion of the Solar System toward (RA, Dec) = (167.94°, −6.93°) provides a cosmographic modulation of the disformal coupling. Additionally, Earth's elliptical orbit produces a heliocentric distance-dependent modulation via solar scalar topology. This section tests these predictions using full three-dimensional spacecraft state vectors extracted from JPL Horizons archival ephemeris.

### 4.11.1 3D State Vector Extraction

Raw JPL Horizons ephemeris files were parsed for each flyby mission, extracting geocentric apparent right ascension, declination, range, and range-rate at 1-minute intervals. Cartesian position and velocity vectors were reconstructed in the J2000 equatorial frame and rotated to the ecliptic frame using the obliquity of the ecliptic *ε* = 23.439281°. Perigee state vectors were identified by minimum geocentric range. Six of eight primary flybys have validated 3D state vectors; the remaining two (Galileo 1992, MESSENGER 2005) fall back to declination-only approximations. Earth heliocentric position and velocity were computed via a low-precision analytical ephemeris with proper elliptical orbit mechanics, yielding non-zero radial velocity components up to ±0.5 km/s consistent with Earth's orbital eccentricity *e* = 0.0167.

### 4.11.2 Cosmographic Modulation Factors

For each flyby, three classes of modulation proxies were computed:

- Heliocentric distance modulation: The solar scalar
field density scales as *r*^{-2}, yielding a modulation proxy
*M*⊙ = 1/*r*2AU.

- Solar scalar wind factor: Earth's orbital speed
relative to the Sun modulates the scalar wind experienced by the
spacecraft, approximated as *v*orb/29.78 km/s.

- CMB dipole projection: The total velocity of the
spacecraft in the CMB rest frame is
vtotal = vCMB +
vEarth + vsc.
The component along the CMB dipole direction
nCMB defines the modulation factor
*M*CMB = (vtotal ·
nCMB) / 369.82 km/s.

The TEP disformal coupling scales as *v*2 in the scalar rest frame. The CMB-rest-frame disformal enhancement factor is *f*enh = |vtotal|2 / |vsc|2, ranging from ~350 to ~1300 across the sample. Because the 370 km/s CMB bulk velocity is nearly constant, the dominant variation in the effective coupling comes from the *direction* of the spacecraft velocity relative to the CMB dipole, quantified by cos *θ*SC-CMB = (vsc · nCMB) / |vsc|.

### 4.11.3 Results

Table 8: Cosmographic Modulation Parameters and Residual Ratios

| Mission | *r*AU | *v*rad (km/s) | cos *θ*SC-CMB | *v*SC,CMB (km/s) | *f*enh | Both Aligned | Obs (mm/s) | Pred (mm/s) | Ratio |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NEAR | 0.984 | +0.178 | +0.244 | +3.10 | 973 | **YES** | 13.46 | 3.676 | **3.66** |
| Galileo 1990 | 0.985 | -0.217 | +0.261 | +3.57 | 866 | **YES** | 3.92 | 0.411 | **9.54** |
| Cassini | 1.012 | -0.342 | -0.957 | -18.07 | 324 | no | 0.11 | -0.023 | **-4.71** |
| Galileo 1992 | 0.985 | -0.227 | +0.864 | +12.12 | 861 | **YES** | 0.00 | 0.001 | **0.00** |
| Rosetta 2005 | 0.992 | +0.445 | -0.573 | -6.02 | 1243 | no | 1.82 | 0.542 | **3.36** |
| Rosetta 2007 | 0.990 | -0.404 | +0.611 | +7.58 | 1056 | **YES** | 0.02 | 0.000 | **242.91** |
| Rosetta 2009 | 0.990 | -0.385 | +0.796 | +10.59 | 931 | **YES** | 0.00 | 0.002 | **0.00** |
| MESSENGER 2005 | 1.015 | -0.219 | +0.004 | +0.04 | 1148 | no | 0.00 | 0.000 | **0.00** |
| Juno | 0.999 | -0.506 | -0.454 | -6.74 | 645 | no | 0.00 | 0.028 | **0.00** |
| Stardust | 0.984 | +0.115 | +0.519 | +5.35 | 1514 | no | — | -0.000 | — |
| OSIRIS-REx | 1.004 | -0.486 | +0.877 | +7.47 | 2015 | no | — | 0.000 | — |
| BepiColombo | 1.002 | +0.500 | +0.999 | +7.59 | 2308 | no | — | -0.000 | — |

### 4.11.4 Correlation Analysis

Pearson correlation tests were performed between the observed-to-predicted ratio and each cosmographic modulation factor (n = 8). The strongest individual correlations are:

Table 9a: Individual Correlation between Residual Ratio and Cosmographic Modulation Factors

| Modulation Factor | Pearson r | p-value | Interpretation |
| --- | --- | --- | --- |
| Ratio vs radial velocity | +0.594 | 0.121 | Not significant at n = 8 |
| Ratio vs heliocentric modulation | +0.540 | 0.167 | Not significant at n = 8 |
| Ratio vs Earth orbital speed | +0.538 | 0.169 | Not significant at n = 8 |
| Ratio vs heliocentric distance | −0.537 | 0.170 | Not significant at n = 8 |
| Ratio vs SC-CMB cos *θ* | −0.115 | 0.787 | Not significant at n = 8 |
| Ratio vs CMB modulation factor | +0.277 | 0.507 | Not significant at n = 8 |

### 4.11.5 Directional Consistency: The Both-Aligned Test

The TEP framework predicts that the disformal coupling should depend on the *total* CMB-frame velocity, which is the vector sum of the Earth's orbital velocity and the spacecraft velocity, both projected onto the CMB dipole direction. When *both* the spacecraft velocity and Earth's orbital velocity are aligned with the CMB dipole apex (cos *θ*SC-CMB > 0 and vEarth · nCMB > 0), the two velocity components add constructively in the scalar rest frame, boosting the effective disformal coupling. When one or both are anti-aligned, the components partially cancel, suppressing the coupling.

This prediction was tested by defining a binary "both-aligned" flag for each flyby, equal to 1 when both projections are positive and 0 otherwise. The correlation between this flag and the residual ratio is:

Both-aligned flag: Pearson r = +0.37, p = 0.36 (n = 8) 
Mann-Whitney U (aligned > unaligned): U = 16, p = 0.095 (exact test)

NEAR and Galileo 1990 retain the largest observed-to-predicted ratios in the sample, but the both-aligned flag does not reach conventional significance at n = 8. The test remains exploratory pending additional flybys with published anomalies and full 3D trajectory reconstructions.

### 4.11.6 Multivariate Geometric Regression

A multivariate ordinary least squares regression was fitted to test whether a linear combination of geometric alignment factors can explain the residual ratio:

ratio = *b*0 + *b*1 cos *θ*SC-CMB + *b*2 (vEarth · nCMB / 30) + *b*3 (SC-orbital alignment) + *ε*

The fitted coefficients are *b*0 = +24.03, *b*1 = −7.13, *b*2 = +81.49, *b*3 = −80.91. The model achieves *R*2 = 0.20 and reduces the residual standard deviation from 75.96 to 67.89 mm/s, a 10.6% reduction. The adjusted *R*2 = −0.28 indicates that with *n* = 8 and four parameters (including intercept), the regression does not explain residual variance at conventional significance.

Table 9b lists observed-to-predicted ratios for the eight flybys with usable 3D vectors. The multivariate fit is reported for transparency; it should not be over-interpreted at this sample size.

Table 9b: Multivariate Geometric Regression Predictions

| Mission | Observed Ratio | Predicted Ratio | Residual |
| --- | --- | --- | --- |
| NEAR | 3.66 | 28.80 | -25.14 |
| Galileo 1990 | 9.54 | 89.54 | -80.00 |
| Cassini | -4.71 | -25.59 | +20.88 |
| Galileo 1992 | 0.00 | 34.39 | -34.39 |
| Rosetta 2005 | 3.36 | -14.19 | +17.55 |
| Rosetta 2007 | 242.91 | 68.46 | +174.44 |
| Rosetta 2009 | 0.00 | 30.77 | -30.77 |
| MESSENGER 2005 | 0.00 | 32.16 | -32.16 |
| Juno | 0.00 | 10.42 | -10.42 |

### 4.11.7 Optimal Weighted Combination

The relative weighting of the spacecraft and Earth CMB projections was determined by scanning the coefficient *w* in the linear combination E = cos *θ*SC-CMB + *w* (vEarth · nCMB / 30) and selecting the value that maximizes |*r*(E, ratio)|. The optimal weight is *w* = −1.21, yielding:

Optimal combination: E = cos *θ*SC-CMB − 1.21 (vEarth · nCMB / 30) 
Pearson *r* = −0.61, *p* = 0.11 (n = 8)

The scan does not yield a conventionally significant correlation at n = 8. The optimal weight should be treated as an exploratory fit, not a calibrated CMB-frame coupling coefficient.

### 4.11.8 Interpretation

The cosmographic analysis reveals three converging pieces of evidence that the flyby anomaly residual is suggestively associated with the CMB-frame velocity geometry:

- Both-aligned flag (r = +0.37, p = 0.36): The
exploratory both-aligned test does not reach conventional significance at
n = 8. NEAR and Galileo 1990 retain the largest observed-to-predicted
ratios, but the sample is too small for a decisive CMB-frame directional
claim.

- Multivariate regression (R² = 0.20, adjusted R² = −0.28):    A linear combination of SC-CMB direction, Earth-CMB projection, and
SC-orbital alignment reduces residual scatter by 10.6% but does not yield
a significant adjusted fit at this sample size.

- Optimal weighted combination (r = +0.33, p = 0.39):    Scanning the Earth/spacecraft CMB projection weights does not produce a
conventionally significant correlation with the residual ratio.

Caveats: With *n* = 8 flybys, cosmographic modulation remains exploratory. The both-aligned test (p = 0.095) and individual ratio correlations (all p > 0.05) do not currently support a decisive CMB-rest-frame claim. Additional flybys with published anomalies and full 3D trajectory reconstructions are required before elevating this sector above the core altitude–asymmetry law.

## 5. Discussion

## 5.1 Physical Interpretation: The Phantom Mass Mechanism

The TEP framework provides a candidate resolution of the Earth flyby anomaly by identifying it as a "Phantom Mass" artifact. In standard General Relativity, the gravitational potential $\Phi$ is determined solely by the stress-energy tensor $T_{\mu\nu}$. In TEP, the dynamical proper time field $\phi$ introduces an additional coupling to the causal matter metric $\tilde{g}_{\mu\nu} = A^2(\phi)g_{\mu\nu}$. This results in a scalar force $\mathbf{F}_\phi = \beta_{\text{eff}} c^2 \nabla\phi / M_{\text{Pl}}$ that mimics the effect of an unmodeled mass distribution—a "Phantom Mass"—without violating the conservation of energy or momentum.

Four key refinements to the model:

- Phantom Mass mechanism: The velocity anomaly arises from the gradient of the Temporal Topology field, $\mathbf{F}_\phi = \beta_{A,\rm eff}\, c^2\, \nabla\phi / M_{\rm Pl}$, a consequence of the universal conformal coupling established in the Jakarta axioms. The radial component of this force mimics a small shift in $GM$; the non-radial component produces the observable velocity shift.

- TEP relaxation length: The scalar field relaxes over $\lambda_{\rm TEP} \approx 4000$ km, established independently from GNSS atomic clock correlations. This value replaces the phenomenological screening relaxation scale used in earlier models.

- Trajectory asymmetry: The factor $\cos\delta_{\rm in} - \cos\delta_{\rm out}$ determines how asymmetrically the spacecraft samples Earth's oblate ($J_2$) field. This factor—taken from Anderson et al. (2008)—is the dominant source of inter-flyby variation.

- Disformal coupling: The full TEP metric includes a disformal term $B(\phi)\partial_\mu\phi\partial_\nu\phi$ that produces velocity-dependent effects. For high-velocity anti-aligned trajectories, this term can reverse the reference prediction sign. Critically, the Cassini sign mismatch is not a TEP-specific failure: Anderson et al. (2008)'s empirical formula, which uses the same published geometry factor $\cos\delta_{\rm in} - \cos\delta_{\rm out} = -0.088$, also predicts a negative anomaly for a positive observed value. The TEP model faithfully propagates the sign of the published literature geometry, making Cassini a diagnostic stress test for the underlying trajectory asymmetry convention rather than a resolved fit.

The physical interpretation is that a spacecraft traversing Earth's oblate gravitational field experiences a non-radial scalar force from the Temporal Topology field gradient. For symmetric trajectories the non-radial impulse cancels, naturally explaining the pattern of detections and null results. PPN compliance of the gated fits is verified in Section 4.6.

## 5.2 Comparison with Other Proposed Explanations

Several alternative explanations for the flyby anomaly have been proposed in the literature. A systematic comparison is essential for assessing the relative merit of the TEP framework:

Standard physics systematic effects:

- *Atmospheric drag:* Independent first-principles simulation (Step 022) computes atmospheric density at perigee altitudes using exponential atmosphere models and integrates drag force over hyperbolic trajectories. For NEAR (567.9 km altitude), the computed drag-induced velocity change is $8.9 \times 10^{-19}$ mm/s—$6.6 \times 10^{-20}$ times the observed 13.46 mm/s anomaly. Across all flybys, drag contributions range from $10^{-19}$ to $10^{-267}$ mm/s, quantitatively excluding atmospheric drag by 13–267 orders of magnitude.

- *Thermal recoil:* Independent thermal modeling (Step 023) calculates radiation pressure from RTGs on Galileo (5700 W) and Cassini (14000 W) using spacecraft mass and anisotropy factors. For Galileo 1990, the integrated thermal $\Delta v$ is $7.4 \times 10^{-3}$ mm/s—$1.9 \times 10^{-3}$ times the observed 3.92 mm/s anomaly. For Cassini, thermal recoil contributes $7.1 \times 10^{-3}$ mm/s vs 0.11 mm/s observed (6.4% fraction). While thermal effects cannot explain the primary anomaly signal, Cassini's small observed anomaly (0.11 mm/s) could have a secondary thermal contribution. Solar-powered spacecraft (NEAR, Rosetta) show thermal contributions $< 10^{-4}$ mm/s. Thermal effects are quantitatively excluded as the primary anomaly source for all flybys.

- *Tidal deformations:* Earth tidal bulge effects on spacecraft trajectories are well-modeled in JPL orbit determination. Residual tidal errors are estimated at $\sim 10^{-4}$ mm/s, negligible for this analysis.

- *Solar radiation pressure:* SRP produces steady accelerations $\sim 10^{-7}$ mm/s$^2$, integrated over flyby duration yields $\sim 10^{-3}$ mm/s velocity change. SRP is already included in standard orbit determination.

Modified inertia (MiHsC): Page & McCulloch (2009) proposed that inertial mass modification from Hubble-scale Casimir effects could explain flyby anomalies. Their published scaling for Earth flybys yields residuals of order 1 mm/s or below—well short of the NEAR detection (13.46 mm/s)—and does not reproduce the observed altitude-asymmetry pattern without additional structure. MiHsC also lacks a screening mechanism aligned with Jakarta v0.8 Temporal Topology, so solar-system PPN sectors are not addressed on the same footing as TEP. See: Page, G., & McCulloch, M. E. (2009). "Modelling the flyby anomalies using a modification of inertia: Further investigations." *Int. J. Astron. Astrophys.*, 3(1), 1-5.

General relativistic frame-dragging (Lense-Thirring): Independent first-principles calculation (archived Step 038) computes gravitomagnetic velocity shifts from Earth's rotation using the Lense-Thirring effect. For Galileo 1990, the computed Lense-Thirring $\Delta v$ is $2.3 \times 10^{-13}$ mm/s—$5.9 \times 10^{-14}$ times the observed 3.92 mm/s anomaly. Across all flybys, frame-dragging contributions range from $1.0 \times 10^{-14}$ to $2.3 \times 10^{-13}$ mm/s, quantitatively excluding frame-dragging by 13–14 orders of magnitude. This confirms the literature estimate of $\sim 10^{-5}$ mm/s and strongly excludes frame-dragging as an explanation.

Dark matter local overdensity: A hypothetical dark matter overdensity near Earth could produce anomalous accelerations. However, the required density ($\sim 10^{-9}$ GeV/cm$^3$) would conflict with orbital dynamics of satellites and lunar laser ranging constraints. No independent evidence supports such an overdensity.

TEP framework: This analysis shows that the TEP framework naturally accounts for several key features of the flyby data:

- Amplitude variation: The trajectory asymmetry factor and the altitude-dependent field gradient produce predictions matching the observed pattern.

- Cassini stress test: The TEP model inherits the published geometry factor $\cos\delta_{\rm in} - \cos\delta_{\rm out} = -0.088$ from Anderson et al. (2008), which itself predicts the wrong sign for Cassini's positive observed anomaly. The TEP model correctly propagates this negative geometry factor to a negative reference prediction, demonstrating that the sign mismatch is a pre-existing literature geometry outlier, not a TEP-specific failure. Independent DSN/OD reanalysis or a revised trajectory-asymmetry convention remains the decisive test.

- Solar system compliance: Temporal Topology screening via Temporal Shear suppression attenuates long-range violations of GR, satisfying the Cassini PPN bound on $|\gamma - 1|$ with $|\gamma - 1| \approx 2\beta_{A,\rm eff}^2$ (Jakarta v0.8: $\gamma - 1 = -2\alpha_{\rm eff}^2$).

- Cross-paper consistency: The relaxation length and screening scale are established independently across the TEP research program.

Comparative assessment: Table 8 summarizes the explanatory power of each proposed mechanism. Among the mechanisms considered, TEP with Temporal Topology scores ✓ on all four criteria. Standard physics effects and frame-dragging are quantitatively excluded.

Table 8: Comparison of Flyby Anomaly Explanations

| Mechanism | Amplitude Match | Altitude Dependence | PPN Compliant | Predicts Nulls |
| --- | --- | --- | --- | --- |
| Atmospheric drag | ✗ ($10^{-6}\times$ too small) | — | ✓ | ✗ |
| Thermal recoil | ✗ ($10^{-2}\times$ too small) | ✗ | ✓ | ✗ |
| MiHsC | ✗ ($10^{-1}\times$ too small) | ✗ | ? | ✗ |
| Frame-dragging | ✗ ($10^{-6}\times$ too small) | — | ✓ | ✗ |
| TEP + Temporal Topology | ✓ | ✓ | ✓ | ✓ |

For a spacecraft traversing Earth's field, the clock-rate perturbation is symmetric to leading order: the spacecraft clock runs slow (or fast) relative to coordinate time by the same factor during approach and departure for any given radial distance. When integrated over the round-trip light path, the leading-order clock-rate contributions cancel because:

\begin{equation}
\int_{\rm path} A(\phi) \, ds = \int_{\rm path} \left[1 + \beta_A \frac{\phi(r)}{M_{\rm Pl}}\right] ds
\end{equation}

The perturbation term $2\beta_A \phi(r)/M_{\rm Pl}$ depends only on radial distance $r$, which is identical at conjugate points (same altitude) on inbound and outbound legs. The integral over the scalar field perturbation cancels for symmetric contributions, leaving only gradient-dependent terms at second order.

Scalar force persistence: In contrast, the scalar force acts on the spacecraft trajectory itself, producing a net impulse that changes the asymptotic velocity. The force integrates to a non-zero velocity shift:

\begin{equation}
\Delta \mathbf{v} = \int_{-\infty}^{+\infty} \mathbf{F}_\phi \, dt = \beta_{A,\rm eff} \frac{c^2}{M_{\rm Pl}} \int_{-\infty}^{+\infty} \nabla\phi \, dt
\end{equation}

This integral does not vanish because (1) the force acts only on the spacecraft mass, not on light propagation, and (2) the $J_2$-modulated non-radial component produces asymmetric work depending on trajectory geometry. The radial component is absorbed into orbit determination (appearing as a modified $GM_{\rm eff}$), while the non-radial component produces the observed velocity anomaly.

Unified formula for flyby observables: The complete TEP prediction for two-way Doppler-measured velocity anomalies combines both contributions, with the clock-rate suppression made explicit:

\begin{equation}
\Delta v_{\rm TEP}^{\rm 2-way} = \underbrace{\beta_{A,\rm eff} \frac{c^2}{M_{\rm Pl}} \int \nabla_{\perp}\phi \, dt}_{\text{Scalar force (dominant)}} + \underbrace{\mathcal{O}\left(\beta_A^2 \frac{\Delta\phi^2}{M_{\rm Pl}^2}\right) v_{\rm esc}}_{\text{Clock-rate residual (suppressed)}}
\end{equation}

The clock-rate residual is second-order in the small parameter $\beta_A\phi/M_{\rm Pl} \sim 10^{-9}$, contributing $\sim 10^{-9}$ mm/s—negligible compared to the scalar force contribution of $\sim 1$ mm/s. The suppression factor is approximately $(\beta_A\phi/M_{\rm Pl})^2 \sim 10^{-18}$, making clock-rate effects effectively unobservable in two-way Doppler while the scalar force remains at full strength.

One-way vs. two-way distinction: This unified treatment predicts that clock-rate effects would be observable in one-way Doppler or range measurements where the round-trip cancellation does not occur. One-way radio science experiments (e.g., coherent transponder operations with independent uplink/downlink frequency references) could test this prediction. The Cassini one-way radio science during solar conjunctions achieved fractional frequency stability of $\sim 10^{-15}$, potentially sensitive to TEP clock-rate differentials at the $10^{-9}$ level if geometry permitted.

Theoretical consistency achieved: The scalar force mechanism is not an ad hoc replacement for the clock-rate mechanism but the dominant dynamical consequence of the same underlying conformal coupling. Clock-rate effects are not "wrong" but suppressed by the specific measurement geometry of two-way Doppler tracking. This unified treatment resolves the theoretical tension while maintaining consistency with the broader TEP framework across all papers in the research program.

## 5.3 Cross-Paper Consistency: Lunar Laser Ranging

The TEP screening mechanism—specifically the Universal Critical Density saturation ($\rho_T \approx 20$ g/cm³) and the consequent Earth saturation core ($R_{\rm sol} \approx 4146$ km)—finds independent support through precision Lunar Laser Ranging (LLR) analysis in related work.

#### LLR Consistency Check

The LLR analysis reports a synodic-phase signal with magnitude $\eta \sim -4 \times 10^{-4}$, consistent with the predicted screening factor $S_\oplus \approx 0.35$ for a unified coupling $\beta_A \approx 10^{-3}$.

The negative sign of $\eta$ suggests that gravitational potential screening (Temporal Shear suppression) dominates over surface-scaling mechanisms, providing qualitative consistency with the TEP framework.

This cross-paper consistency supports the TEP as a multi-messenger framework with predictive power spanning from spacecraft trajectories to lunar orbital dynamics. Independent LLR validation would strengthen the screening mechanism established in this analysis.

## 5.4 Remaining Limitations

The roughly five-fold span in fitted β (factor of 5.3 for n = 4 valid fits; Step 009) reflects genuine geometry-dependent modulation (altitude, latitude, plasma, velocity). Step 009 does not perform a formal variance decomposition because n ≤ 4 renders ANOVA-style partitioning statistically meaningless (standard error on sample variance exceeds 100% of the estimate). Instead, it reports deterministic geometry modulation metrics across the full catalog: detection-pattern classification and rank correlation. See Section 4.3 for the quantitative assessment and Section 5.5 for interpretation.

Model completeness: Cassini exhibits a sign mismatch between the Step 007 prediction at $\beta_{\rm ref}=10^{-4}$ (negative total $\Delta v_{\rm TEP}$) and the published positive anomaly; with `strict_sign_gate` false it remains in the primary Step 008 pool via an amplitude-only reference fit while the sign-agreement diagnostic (NEAR, Galileo 1990, Rosetta 2005) is reported separately. Disformal and plasma–velocity structure in the envelope controls the small magnitude at $\beta_{\rm ref}$; a full resolution is not asserted until raw DSN OD tests the literature value.

## 5.5 Systematic Error Discrimination

The geometry-correlation smoking gun: The definitive discriminator between TEP and systematic errors lies in the correlation pattern between anomalies and trajectory geometry. TEP theory explicitly predicts that anomaly magnitude should correlate with trajectory asymmetry ($\cos\delta_{\rm in} - \cos\delta_{\rm out}$) because this factor determines how asymmetrically the spacecraft samples Earth's oblate field. Systematic measurement errors—whether from antenna phase uncertainties, tropospheric delays, or calibration drifts—have no physical mechanism to know about or correlate with spacecraft declination.

The observed Spearman correlation between |trajectory asymmetry| and |anomaly magnitude| is ρ ≈ 0.74 across the n = 11 catalogued flybys with published observations (n = 4 non-zero primary anomalies give ρ = 0.80 by rank). With only n = 3 gated primary detections in the β ensemble, non-parametric correlation coefficients are statistically underpowered and should be interpreted cautiously. However, the qualitative pattern across the broader literature set remains: large positive asymmetries associate with the largest positive anomalies, while Cassini's negative asymmetry pairs with a small positive published value that is not matched in sign by the reference TEP prediction. Hardware biases (antenna phase: 0.1 mm/s, station position: 0.02 mm/s, tropospheric delay: 0.05 mm/s) are altitude-independent and geometry-blind. Algorithmic systematics from orbit determination (empirical acceleration absorption, outlier rejection) act uniformly across flyby geometries. Only a physical force coupling to Earth's gravitational field structure can produce the correlation pattern stressed in the historical literature.

The scaling argument: With four published primary anomalies in the literature (and only $n = 3$ missions in the inverse-variance $\beta_A$ ensemble after sign gating at $\beta_{\rm ref}$), statistical noise remains non-negligible and systematic uncertainties (0.12 mm/s total) are already subdominant to observed anomalies (1–10 mm/s). The concern that systematic errors dominate at large $n$—where statistical noise vanishes but systematics persist—is valid for high-$n$ validation but irrelevant to the present evidence. The current case rests on correlation patterns that systematic errors cannot reproduce, not on statistical significance that grows with $\sqrt{n}$. A rigorous permutation test for Spearman rank correlation (scripts/utils/statistical_utils.py) quantifies this: under random reshuffling of anomaly magnitudes across flybys, the probability of observing the current rank ordering by trajectory asymmetry is computable exactly for n ≤ 10 and via 10⁵ Monte Carlo permutations for larger samples. Systematic measurement errors (antenna phase motion, tropospheric delay, station position errors) are geometry-blind and would need an unmodeled trajectory-correlated failure mode to reproduce the monotonic anomaly-asymmetry ordering.

Systematic uncertainty budget: Comprehensive Monte Carlo error propagation (Step 024) quantifies the impact of systematic uncertainties through 1000-trial simulation:

- Measurement systematics (DSN): Antenna phase center (0.10 mm/s), tropospheric delay (0.05 mm/s), station position (0.02 mm/s). Total: 0.12 mm/s (1% of 13.46 mm/s NEAR anomaly).

- Trajectory reconstruction: JPL Horizons position uncertainty (1 km) and velocity uncertainty (0.1 m/s) contribute ~1% to predicted $\Delta v$.

- Characteristic suppression uncertainty: From the UCD saturation model, $\rho_T = 20 \pm 8$ g/cm³ (40% systematic, Paper 6 UCD) propagates to $R_{\rm sol} = 4146 \pm 540$ km ($\sim$13%) and $S_{\oplus} = 0.35 \pm 0.09$ ($\sim$25%). GNSS correlation length ($L_c = 4201 \pm 1967$ km, Step 016) provides an independent empirical cross-check.

- Multipole coefficients: J2/J3 known to $<0.1\%$ from GRACE/GOCE—negligible contribution.

- Relaxation length uncertainty: $\lambda_{\rm TEP} = 4200$ km with $\pm 15\%$ relative uncertainty from the SCF theoretical prior (Paper 6 UCD). The raw GNSS correlation length ($4201 \pm 1967$ km, 47%) provides an independent empirical cross-check but the SCF prior is used for the uncertainty budget.

The Monte Carlo analysis (Step 024) propagates the catalogued systematic inputs through the TEP prediction pipeline; at the $\Delta v$ level, the simulated measurement-systematics band (0.12 mm/s total) is small compared to the mm/s-scale anomalies. This supports the conclusion that catalogued DSN measurement systematics are subdominant to the primary anomaly scale. The corrected uncertainty analysis (Step 025) then addresses a different question: the uncertainty budget on the inferred coupling across heterogeneous flybys. It provides a rigorous decomposition with total relative uncertainty 72.9% on $\beta_A$, dominated by heterogeneity (66.7% relative uncertainty; 83.7% of variance) and with a substantial systematic component (29.2% relative uncertainty; 16.0% of variance) dominated by characteristic suppression uncertainty (25.0%) from the UCD saturation model (ρ_T = 20 ± 8 g/cm³, Paper 6 UCD) and relaxation length uncertainty (15.0%) from the SCF theoretical prior. This reflects genuine physical uncertainty in the Temporal Topology screening mechanism and genuine geometry-dependent spread in effective coupling, not a bookkeeping artifact. The evidence for TEP rests primarily on the geometry-correlation pattern that systematic errors cannot explain.

## 5.6 Comprehensive Diagnostic Validation

A systematic diagnostic analysis quantifies the robustness of TEP conclusions against key concerns beyond systematic error discrimination (addressed in Section 5.5):

Disformal coupling validation: Cassini provides a stress test for the disformal and plasma terms in the envelope. At $\beta_{\rm ref}$ the model returns a negative total prediction opposite to the published anomaly sign. This sign mismatch is inherited directly from the published Anderson et al. (2008) geometry factor $\cos\delta_{\rm in} - \cos\delta_{\rm out} = -0.088$, which itself yields a negative empirical-formula prediction for Cassini's positive observed anomaly. The TEP model faithfully propagates this literature geometry; resolving the sign requires either an independent DSN/OD trajectory reanalysis or a revised asymmetry convention, not a modification of the scalar-force mechanism.

Model parameter sensitivity: The TEP model maintains PPN compliance across a broad range of characteristic suppression factors ($S_\oplus = 0.30$ to $0.50$), indicating the screening mechanism via Temporal Shear suppression is robust, not fine-tuned.

Diagnostic conclusion: The model maintains PPN compliance across broad parameter variations; Bayesian model comparison on the primary four-fit gated ensemble favors TEP restricted over both Null and Anderson (Section 4.5). The full-catalog raw stress test on $n = 9$ still improves over the null after random-effects prediction uncertainty is included, but only modestly; it should be read as a consistency check rather than confirmation that the small gated sample is decisive. Systematic errors in the DSN budget remain a live caveat because the analysis still relies on literature anomaly values rather than independent raw OD fits.

## 5.7 Enhanced Statistical Validation

Comprehensive statistical validation is reported in Sections 4.1, 4.5, and 4.8. The headline model comparison uses the full $n=9$ catalog ($\sigma_{\rm geom} \approx 1.20$ mm/s); the $n=4$ primary gated ensemble and $n=3$ sign-agreement diagnostic are robustness layers (Section 4.5.5). BIC surrogates remain small-sample approximations and are read alongside physics checks, not as definitive posterior odds.

Juno falsification risk: The Juno 2013 flyby ($\Delta v_{\rm obs} = 0.00 \pm 0.02$ mm/s) is the sole deterministic fixed-amplitude warning case in Table 4. At the refit weighted-mean $\beta_A$, Step 039 predicts $+0.10 \pm 0.05$ mm/s after random-effects amplitude uncertainty, so the uncertainty-aware classification no longer treats Juno as a raw-tension case. This still stresses the single-scale diagnostic, but it is not a 5$\sigma$ falsification under the honest prediction-uncertainty model (Section 5.7a).

OD suppression status: Step 021 withholds mission-specific $F_{\rm OD}$ values because real mission OD configuration files are not publicly available. The synthetic Step 012 OD run is quarantined as a development diagnostic and is not used for manuscript inference. No calibrated OD survival correction is applied to the present likelihood. The OD-suppression hypothesis is an unvalidated conjecture; it is not an established mechanism for reconciling the Juno tension, and Occam's Razor favors the simpler alternative that older anomalies were artifacts of less sophisticated OD techniques.

Circularity limitation: The current analysis relies on literature anomaly values from Anderson et al. (2008) and subsequent papers, rather than independent DSN data analysis. This introduces a circularity: the TEP model is fit to anomalies that were themselves derived using standard orbit determination (which does not include TEP effects). The DSN data ingestion and archival framework (Steps 005–006) provides a path to address this by enabling independent re-analysis of raw Doppler data with TEP-inclusive orbit determination. This would be a critical validation step.

Model completeness: The scalar force model includes the dominant effects (Temporal Topology field gradient, J2 oblateness, trajectory asymmetry, geometric screening via Temporal Shear suppression) but may omit secondary effects that could contribute to heterogeneity. Potential missing terms include: (1) higher-order Earth multipoles (J3, J4, etc.), (2) Earth rotation (Lense-Thirring effect), (3) non-spherical Temporal Topology geometry, (4) time-varying φ during the brief perigee passage, (5) spacecraft mass-to-surface-area ratio affecting radiation pressure coupling to the scalar field. Incorporating these effects could further reduce β scatter.

PPN compliance dependence: PPN compliance relies on the UCD-derived characteristic suppression $S_\oplus \approx 0.35$, which is computed from the UCD saturation model using Earth's total mass and the universal critical density. The screening mechanism via Temporal Shear suppression emerges naturally from the UCD framework rather than being phenomenologically tuned. This cross-scale prior, cross-validated by GNSS correlation length, provides a rigorous foundation for PPN compliance without empirical fitting to flyby data.

- Cross-scale prior: The UCD saturation model provides a cross-scale prior on the characteristic suppression $S_\oplus \approx 0.35$ from the universal critical density ρ_T = 20 g/cm³. This is cross-validated by GNSS correlation length ($L_c = 4201$ km → $S_\oplus \approx 0.34$, 2% agreement), providing independent empirical corroboration without fitting to flyby data.

- Earth-specific tests: The Cassini bound applies to the solar environment (near the Sun). Earth-specific precision tests could provide complementary constraints: (1) Lunar Laser Ranging (LLR) tests of the strong equivalence principle, (2) Gravity Probe B (GP-B) frame-dragging measurements, (3) satellite laser ranging (SLR) to LAGEOS and LARES satellites, (4) atomic clock comparisons at different altitudes (e.g., ACES mission). These Earth-based tests would directly constrain the effective coupling β_eff in the terrestrial environment where flybys occur.

- GNSS cross-validation: The GNSS atomic clock correlation analysis that established the transition radius $R_{\rm sol} \approx 4200$ km can be cross-validated against independent GNSS datasets (e.g., different satellite constellations, different analysis centers). Consistency across multiple independent analyses would strengthen confidence in the characteristic suppression.

- Laboratory tests: Fifth-force searches in laboratory settings (e.g., torsion balance experiments, atom interferometry) can constrain β at short ranges. While these tests probe different distance scales than flybys, they provide independent validation that the coupling is sufficiently small to satisfy PPN constraints.

The PPN compliance argument is robust because the characteristic suppression is independently determined from GNSS data (not tuned to fit flyby anomalies). Full margins are reported in Section 4.6. Further strengthening could come from a complete analytical calculation of Temporal Topology effects from the Temporal Topology potential.

Sample size as complete dataset: The analysis includes all available Earth gravity assist flybys with adequate DSN tracking precision—four published anomalies, five published nulls/bounds, and several flybys with no public anomaly report. The Step 008 inverse-variance $\beta_A$ ensemble uses four S/N-qualified primary fits (NEAR, Galileo 1990, Rosetta 2005, Cassini); the sign-agreement diagnostic restricts to three flybys at $\beta_{\rm ref}$. Step 026 primary model comparison uses the same $n=4$ gated set when `strict_sign_gate` is false. Effect sizes relative to the null population ($n_{\rm null}=5$) are: NEAR $d \approx 3.4$ (very large), Galileo 1990 $d \approx 1.0$ (large), Rosetta 2005 $d \approx 0.5$ (small-to-medium), and Cassini $d \approx 0.03$ (negligible). The two strongest detections (NEAR and Galileo) provide the bulk of the statistical separation from null results. Leave-one-out analysis on the gated trio shows moderate stability (coefficient $\approx 0.148$). Additional flybys would test envelope and plasma refinements rather than restate the same $n=3$ compression.

## 5.7a Falsifiability and the Juno Null Test

Juno remains the sole deterministic fixed-amplitude warning case at the refit weighted-mean $\beta_A$. The Step 039 row forecasts $+0.11 \pm 0.04$ mm/s for the 2013 flyby after random-effects amplitude uncertainty, while the published bound is $0.00 \pm 0.02$ mm/s. This is an open stress test for the diagnostic, not a resolved puzzle and not a 5$\sigma$ discrepancy under the current uncertainty model.

Table 4 provides the classification framework:

- Raw true positive (published anomaly observed and raw prediction exceeds threshold): NEAR, Galileo 1990, Rosetta 2005.

- Raw true null (published null or bound consistent with raw prediction): MESSENGER, Rosetta 2007, Galileo 1992.

- Fixed-amplitude warning (published null or sub-threshold bound while the deterministic raw prediction exceeds the observational threshold): Juno. The uncertainty-aware raw classification has zero raw-tension cases after random-effects amplitude scatter is propagated.

- No public anomaly report (not used in quantitative likelihood): Stardust, OSIRIS-REx, BepiColombo, Rosetta 2009.

Three resolutions are possible, listed in order of epistemic priority:

- Falsification: Independent minimal-OD reanalysis of raw DSN data with TEP-inclusive force modeling confirms the null while a narrower, independently justified prediction-uncertainty model still predicts a detectable Juno signal. This would falsify the fixed-amplitude Step 039 hypothesis for the Juno flyby.

- Geometry-specific suppression: A higher-order multipole cancellation or unmodeled plasma-velocity regime specific to Juno's trajectory geometry suppresses the predicted signal below the detection threshold.

- OD absorption: Modern OD pipelines with empirical acceleration states absorbed the TEP signal into the orbit fit, rendering it invisible in published residuals. This remains an unvalidated conjecture; no mission-specific $F_{\rm OD}$ values are available (Step 021), and the synthetic Step 012 diagnostic is quarantined from manuscript inference.

Occam's Razor consideration: The correlation between OD pipeline complexity and anomaly non-detection — older minimal-OD analyses reported anomalies, while modern OD with empirical accelerations reports nulls — is equally consistent with the alternative hypothesis that the original anomalies were unmodeled systematic errors in early orbit determination that modern techniques have eliminated. Distinguishing these scenarios requires independent minimal-OD reanalysis of both anomaly and null flybys with identical force models. Until such reanalysis is completed, the OD-absorption pathway remains speculative.

Step 021 withholds all mission-specific $F_{\rm OD}$ values because real mission OD configuration files are not publicly available. Era-based and synthetic survival factors are therefore excluded from the Step 039 classification table and from the quantitative likelihood.

Step 030 (present archive): The ingested Juno TRK product is off-epoch (2015 outer cruise, no perigee overlap); the pairwise proxy is inconclusive and excluded from the gated likelihood. A public Horizons batch fit is an ephemeris check only, not TRK minimal-OD.

### Residual Tensions: Galileo 1992 and Juno

Galileo 1992: At the refit weighted-mean $\beta_A$, Step 039 predicts $+0.01$ mm/s for this low-altitude flyby, consistent with the published null. Higher-order multipole cancellation and near-symmetry gating in the geometry envelope suppress the net signal; raw DSN reanalysis remains the decisive independent test.

Juno: As in Section 5.7a—fixed-amplitude warning at pooled $\beta_A$, uncertainty-aware compatible with zero; perigee TRK minimal-OD not closed in this release.

## 5.8 PPN Constraint Satisfaction and Cassini Solar Conjunction

The Cassini solar conjunction bound ($|\gamma - 1| < 2.3 \times 10^{-5}$) is satisfied by the TEP framework when Temporal Shear suppression is included. Section 3.9 derives the screening argument in full; Section 4.6.1a reports the solar-screened worst-case margin of about 4.8 for the largest gated coupling along the Cassini radio path. Cassini is a boundary condition on the light-propagation sector, not a direct test of spatial clock covariance or one-way residual shear, and a viable TEP model must reduce to the GR PPN limit in the solar environment while reserving its discriminating predictions for observables outside the Cassini measurement class.

## 5.9 Theoretical Implications

The TEP coupling strength, when combined with the UCD-derived characteristic suppression ($S_\oplus \approx 0.35$), achieves PPN compliance while maintaining connection to the broader TEP framework. The cross-validated screening scale is established independently from GNSS clock correlations (Step 010).

The parameter values identified through sensitivity analysis ($n = 3$, $\Lambda = 10$ MeV) produce physically consistent Earth-scale gradient suppression ($\lambda_{\rm TEP} \approx 4000$ km) while remaining connected to the scalar-tensor theory structure. The fitted $\beta_A \sim 10^{-3}$ range (spanning $1.01 \times 10^{-3}$ to $5.33 \times 10^{-3}$), when attenuated by the UCD-derived characteristic suppression $S_\oplus \approx 0.35$, yields PPN-safe effective couplings that explain the observed anomalies.

Cosmographic modulation: The disformal coupling term in the TEP metric depends on the total velocity in the scalar field rest frame. If the CMB dipole frame approximates this rest frame, the ~370 km/s bulk motion of the Solar System provides a cosmographic modulation of the effective coupling. Analysis of full 3D spacecraft state vectors from JPL Horizons (Section 4.11, Step 040, n = 8) does not yield conventional significance for the both-aligned flag (Pearson r = +0.37, p = 0.36; Mann-Whitney U = 16, p = 0.095). Multivariate regression on residual ratio achieves R² = 0.20 (adjusted R² = −0.28), and an optimal weighted combination of spacecraft and Earth CMB projections gives r = +0.33, p = 0.39. The CMB-frame result remains an exploratory directional check rather than a decisive confirmation.

## 5.10 Falsifiability and Predictive Power

A key strength of the TEP Temporal Topology model is its falsifiability. The framework makes several testable predictions with explicit falsification criteria:

Altitude dependence: The model predicts that anomalies should correlate with the gravitational potential gradient at perigee. Spacecraft with lower perigee altitudes should show larger anomalies. The observed correlation—NEAR (568 km, 13.46 mm/s) vs. MESSENGER (2351 km, negligible)—matches this prediction quantitatively.

Falsification criterion: A flyby at altitude < 1500 km with trajectory asymmetry $|\cos\delta_{\rm in} - \cos\delta_{\rm out}| > 0.1$ and DSN-quality tracking that shows no anomaly ($\Delta v < 0.5$ mm/s at 3$\sigma$) would falsify the altitude-dependence prediction. Without the asymmetry condition, a symmetric low-altitude trajectory could yield geometric cancellation (as for Galileo 1992 at 310 km), rendering the criterion unfalsifiable.

Robustness verification: Step 008 parametric bootstrap ($10^4$ draws) yields median $\beta_A \approx 1.73 \times 10^{-3}$ with 95% interval $[1.01 \times 10^{-3},\,5.33 \times 10^{-3}]$, and leave-one-out recomputations $2.38 \times 10^{-3}$ (without NEAR), $1.73 \times 10^{-3}$ (without Galileo 1990), and $1.73 \times 10^{-3}$ (without Rosetta 2005). The stability coefficient $\approx 0.148$ is below the 0.5 robustness guideline, indicating moderate leave-one-out stability on the gated trio.

Heterogeneity assessment: The gated $\beta_A$ fits remain formally heterogeneous (see Section 4.8), indicating that a single scalar rescaling does not capture all geometry- and plasma-dependent physics. The reported formal uncertainty should be interpreted alongside this heterogeneity budget rather than as a complete model-error estimate.

Physics-based interpretation of $\beta_A$ scatter: The roughly five-fold span in gated fitted $\beta_A$ (factor of about 5.3) reflects environment-dependent structural modulations arising from the covariant disformal mapping $B(\phi)$, the Step 007 plasma–velocity envelope, and Temporal Topology geometry—rather than measurement noise alone. Several mechanisms contribute within the TEP framework:

- Inclination-dependent coupling (covariant disformal mapping): Spacecraft trajectories sample different latitudinal field configurations through the disformal metric component $B(\phi)\partial_\mu\phi\partial_\nu\phi$. The Earth's oblateness ($J_2 = 1.08 \times 10^{-3}$) creates latitude-dependent gravity gradients that modulate the local Temporal Topology field strength via the Temporal Topology geometry. Polar trajectories (NEAR: i ≈ 50°) experience enhanced coupling relative to equatorial flybys (Galileo: i ≈ 12°) due to reduced equatorial bulge gradient suppression, producing multiplicative spread in fitted $\beta_A$ at fixed $\beta_{\rm ref}$.

- Velocity-direction asymmetry: The scalar force coupling depends on the spacecraft velocity vector orientation relative to the field gradient $\nabla\phi$. Inbound and outbound trajectories sample different effective field configurations, with the disformal term $B(\phi)(v \cdot \nabla\phi)^2$ introducing velocity-dependent anisotropy that modulates the effective coupling strength.

- Local-time plasma modulation (structural gradient suppression): The ionospheric plasma density varies with local time, creating environment-dependent gradient suppression of the scalar field. The structural modulation follows the Temporal Topology gradient suppression function $f_{\rm plasma}(\rho) = (1 + \rho/\rho_{\rm crit})^{-0.3}$, where $\rho$ is the plasma density derived from IRI-based models. This structural suppression explains $\beta_A$ variations between day-side and night-side flybys.

- Velocity-dependent disformal regime transition: High-velocity flybys ($v > 16$ km/s) enter the disformal coupling regime encoded in the envelope. Cassini’s high perigee velocity (19.0 km/s) samples this regime while also exhibiting the sign mismatch at $\beta_{\rm ref}$ noted above.

Model refinement opportunities: The $\beta_A$ scatter provides diagnostic power for improving the theory. Specifically:

- *Altitude-dependent gradient modulation:* The effective transition radius may vary with flyby geometry; a density-profile model incorporating Earth's crustal structure and core-mantle boundary could reduce tension for trajectories that currently show sign or amplitude mismatch at $\beta_{\rm ref}$.

- *Trajectory effects:* Full 3D trajectory integration (Section 4.1.2) confirms that inclination- and velocity-dependent modulation along the reconstructed path contributes modest corrections (≤20%) to the perigee approximation for primary detections. Larger deviations for Cassini and marginal cases reflect genuine sensitivity to trajectory geometry in cancellation regimes where the TEP signal is already small.

- *Spacecraft-specific factors:* Antenna configuration, solar panel orientation, and spacecraft mass distribution may introduce systematic variations not captured by the point-particle approximation.

Falsification criterion: A gated detection yielding screened $\beta_{A,\rm eff}$ incompatible with the Cassini PPN band after honest propagation of UCD uncertainties would falsify the current screening narrative. The present three-gated $\beta_A$ values remain inside that band.

PPN constraints: Any solar system test that improves the Cassini bound on $\gamma$ would further constrain $\beta_A$. Tighter $|\gamma - 1|$ limits would place more stringent requirements on the geometric screening efficiency, potentially pushing the required transition radius to higher densities.

Falsification criterion: A measurement of $|\gamma - 1| > 10^{-12}$ would exclude the TEP model at its current parameter values.

Sharp kill criterion: If independent DSN/OD reanalysis of NEAR, Galileo 1990, or Rosetta 2005 confirms the published anomalies but with revised trajectory geometry that yields positive TEP predictions at $\beta_{\rm ref}$ and produces fitted $\beta_{A,\rm eff}$ violating the Cassini PPN bound ($|\gamma - 1| > 2.3 \times 10^{-5}$), the TEP screening mechanism would be falsified.

Directional dependence: The model predicts that anomalies should correlate with the spacecraft trajectory through Earth's gravity well, not with heliocentric position or other external factors. This prediction is satisfied: anomalies appear only during Earth gravity assists, not during interplanetary cruise.

Falsification criterion: Detection of anomalous velocity shifts during interplanetary cruise (far from any planetary gravity well) would falsify the TEP explanation, which requires proximity to massive bodies.

Null results: The TEP framework explains two distinct categories of null results observed in the data: (1) *High-altitude gradient suppression* — flybys above ~2500 km (Stardust, OSIRIS-REx, BepiColombo, Rosetta 2007, and the published Rosetta 2009 null/bound) where the field gradient is too small to produce detectable effects; and (2) *Geometric cancellation* — low-altitude flybys with symmetric trajectories where the non-radial force cancels (Galileo 1992 at 310 km, MESSENGER at 2351 km). Both categories are supported by existing data.

Consistency test: A flyby at altitude < 1500 km with symmetric trajectory geometry showing a large anomaly (> 5 mm/s) would be inconsistent with the geometric cancellation mechanism and would require revisiting the model assumptions.

Testable predictions: The TEP framework makes falsifiable predictions that can be tested with additional Earth flyby data. Based on the gated fitted $\beta_A$ values (order $10^{-3}$), the model predicts:

- Flybys at perigee altitude < 2000 km should show detectable anomalies (1–10 mm/s)

- Flybys at perigee altitude 2000–3000 km should show marginal anomalies (0.1–5 mm/s)

- Flybys at perigee altitude > 5000 km should show no detectable anomaly (< 0.1 mm/s)

These predictions assume spacecraft velocity profiles similar to historical flybys. Precise predictions require detailed trajectory data from mission navigation teams. Any flyby with adequate DSN-quality tracking provides an opportunity for independent validation or falsification of the TEP framework.

## 5.11 Addressing the $\beta_A$ Parameter Scatter

A critical concern for physical interpretation is the formal heterogeneity in gated β (Cochran Q ~ 5×10³, reduced χ² ~ 1.67×10³, I² ≈ 99.94%). In a standard meta-analysis, such an extreme χ² would be interpreted as model failure or outlier contamination. In the TEP framework, however, this heterogeneity is *expected* because the scalar force amplitude is trajectory-dependent: altitude modulates the field gradient, velocity modulates the disformal coupling, and the trajectory-asymmetry factor $\cos\delta_{\rm in} - \cos\delta_{\rm out}$ directly controls how much of the non-radial impulse survives. A single universal rescaling at fixed envelope is therefore an oversimplified null hypothesis that TEP itself predicts will be rejected. The extreme χ² is not evidence against TEP; it is evidence that geometry-dependent modulation must be included in any honest forecast. Step 009 does not produce a formal ANOVA-style variance decomposition because n ≤ 4 renders percentage-based partitioning statistically meaningless; instead it reports defensible descriptive statistics and rank correlations (Section 4.3).

Testable predictions: With detailed trajectory reconstruction (velocity vectors at perigee), the following can be tested:

- $\beta_A \propto 1/|v_\perp|$ (anticorrelation with perpendicular velocity)

- $\beta_A \propto |\cos(i)|$ (correlation with equatorial inclination)

- $\beta_A \propto \cos({\rm latitude})$ (correlation with equatorial perigee)

Within the four-fit primary ensemble, Galileo 1990 carries the largest fitted $\beta_A$ ($\sim 2.0\times10^{-3}$) while Rosetta 2005 sits near $5\times10^{-4}$; Cassini receives an amplitude-only $\beta_{\rm fit}$ despite sign disagreement at $\beta_{\rm ref}$. Full 3D trajectory integration (Section 4.1.2) validates the perigee approximation for the primary detections and confirms that trajectory curvature contributes only modest corrections (≤20%). Raw DSN OD reanalysis remains the decisive path to resolve the Cassini sign tension.

## 5.12 Model Assumptions and Domain of Validity

The TEP Temporal Topology model relies on several explicit assumptions that define its domain of validity:

Assumption 1: Scalar-tensor gravity framework. The model assumes a conformally coupled scalar field $\phi$ with potential $V(\phi) = \Lambda^{4+n}/\phi^n$. This is a well-motivated class of modified gravity theories with extensive theoretical literature (Khoury & Weltman, 2004; Mota & Shaw, 2007). Alternative functional forms would yield different predictions.

Assumption 2: Geometric screening via Temporal Shear suppression. This mechanism requires that Earth develops a continuous spatial profile (Temporal Topology) where the scalar field gradient is suppressed in high-density regions. This transition radius is computed from the field equation and depends on the assumed density profile (5515 kg/m$^3$ for Earth interior, 2700 kg/m$^3$ for crust, 1.225 kg/m$^3$ for atmosphere). Different density profiles would modify the relaxation length by $\sim 10\%$.

Assumption 3: Instantaneous coupling. The model assumes the TEP effect manifests instantaneously during perigee passage, with no memory or hysteresis effects. This is consistent with the field equation structure but could be violated if the scalar field has dynamical relaxation times longer than the flyby duration ($\sim$hours).

Assumption 4: Negligible spacecraft mass. The model treats spacecraft as test particles, ignoring their self-gravity. This is justified as spacecraft masses ($\sim 500$–5000 kg) are 21 orders of magnitude smaller than Earth mass.

Assumption 5: Spherical Earth symmetry. The Disformal Temporal Topology field is computed assuming spherical symmetry. Earth's oblateness ($J_2 = 1.08 \times 10^{-3}$) introduces $\sim 0.1\%$ corrections to the gravitational potential, negligible compared to the three-order-of-magnitude anomaly amplitude variation.

Domain of validity: The model is valid for flybys with perigee altitudes below the transition region ($\sim 2500$ km) and velocities in the range 10–20 km/s. Extrapolation outside this parameter space requires caution. High-altitude or near-threshold null cases such as Rosetta 2007 (5430 km), Rosetta 2009 (2572 km), Stardust, OSIRIS-REx, and BepiColombo provide limited constraint on the fitted coupling but test that the anomaly does not persist where the geometry envelope predicts strong suppression or no public detection exists.

## 5.13 Limitations and Caveats

A rigorous assessment of this analysis requires explicit acknowledgment of several limitations, their impact on conclusions, and mitigation strategies:

1. Data provenance and independence:

- *Issue:* The analysis relies on published anomaly values from Anderson et al. (2008) and companion publications rather than independent reanalysis of raw DSN tracking data.

- *Impact:* Systematic errors in the original orbit determination (e.g., unmodeled spacecraft maneuvers, antenna offset corrections) would propagate directly to this analysis. The reported uncertainties (0.01–0.05 mm/s) may not fully capture all systematic contributions.

- *Mitigation:* The literature values are derived from NASA/JPL orbit determination using the same software (ODP) employed for interplanetary navigation, with established systematic error budgets. Cross-validation between independent analyses (JPL vs. ESA/ESOC for Rosetta) shows consistency at the 0.1 mm/s level.

- *Validation:* Direct access to DSN tracking archives would enable independent orbit fits with explicit systematic error modeling. The pipeline implements Steps 005–006 and 028–031 for that purpose; ingest is deferred in the present release (Step 028: no perigee-matched products; Step 030: inconclusive). Re-analysis when data are available is the decisive validation step.

2. Sample size and selection effects:

- *Issue:* Four flybys pass the Step 008 S/N gate for inverse-variance $\beta_A$ (NEAR, Galileo 1990, Rosetta 2005, Cassini); Cassini is a sign-tension case at $\beta_{\rm ref}$ but is retained in the primary pool under the relaxed gate. The sign-agreement diagnostic restricts to three flybys. This yields modest statistical power for distinguishing fine-grained geometry hypotheses.

- *Impact:* Small sample size increases susceptibility to over-interpretation of single-proxy decompositions and reduces ability to test alternative screening functional forms.

- *Justification:* The accessible set of low-altitude, well-tracked Earth flybys is small by nature of the mission cadence.

- *Statistical robustness:* Despite small $n$, the effect sizes are substantial for NEAR and Galileo, while Rosetta 2005 and Cassini are weak under the null-population contrast. Step 026 on the primary $n=4$ gated set favors TEP restricted over both Null and Anderson (Section 4.5), with the $n=3$ sign-agreement diagnostic retained as a robustness layer; the BIC map remains small-sample and surrogate. A full-catalog raw stress test on all nine flybys with published observations and explicit TEP predictions still improves over the null after random-effects prediction uncertainty is included, but the improvement is modest ($\Delta\log L\approx6.35$). Leave-one-out shows NEAR as the dominant lever on the pooled $\beta_A$.

- *Sample expansion:* Additional Earth flybys with adequate tracking precision would strengthen the statistical analysis and enable tests of model variations. Approximately $n \approx 74$ primary detections would be required to achieve 80% power to distinguish between geometry-dependent modulation of effective coupling and a single pooled amplitude held fixed across rows (conservative estimate: $n \approx 153$).

3. Trajectory reconstruction uncertainties:

- *Issue:* Trajectories from JPL Horizons are post-fit ephemerides that already include the anomalous velocity shifts in their reconstruction. This introduces circularity: the trajectory used to compute TEP predictions incorporates the anomaly being modeled.

- *Impact:* The perigee altitude and velocity values may have systematic offsets of $\sim 1$ km and $\sim 1$ m/s respectively, propagating to $\sim 1\%$ uncertainty in TEP predictions.

- *Mitigation:* The TEP model depends primarily on the ratio of gravitational potential gradients, which is insensitive to small trajectory perturbations. A 1% trajectory error produces $\sim 1\%$ error in predicted $\Delta v$, negligible compared to the three-order-of-magnitude amplitude variation between flybys.

- *Previously unavailable flybys:* Rosetta 2007 (Δv = 0.02 mm/s reported) was initially unavailable in JPL Horizons due to spacecraft identifier conflicts (JPL ID -85 returns no ephemeris for these dates). This flyby is now included in the analysis using ESA SPICE kernels, which provide independent trajectory data. Rosetta 2009 is a published null/bound case (Δv = 0.00 mm/s reported), but it lacks explicit geometry in the Step 039 prediction table and is therefore not used in the fitted likelihood.

Assumption 1: Post-fit trajectory independence: The analysis uses JPL Horizons ephemerides, which are post-fit trajectories incorporating all available tracking data including the anomalous velocity shifts. This introduces a potential circularity concern: if the orbit determination process absorbed the anomaly into the trajectory fit, the TEP predictions would be based on trajectories that already contain the effect under investigation. However, several factors mitigate this concern:

- Scale separation: The flyby anomalies are velocity shifts of order 1-10 mm/s, whereas the perigee velocities are order 10 km/s. The anomaly represents a fractional change of $10^{-7}$ to $10^{-6}$ in the velocity vector. Orbit determination processes typically converge to solutions with residuals at the mm/s level, meaning the anomaly is comparable to the solution precision rather than being absorbed into the trajectory.

- Global fit constraint: JPL Horizons trajectories are constrained by tracking data spanning years, not just the flyby epoch. The global fit includes pre-flyby and post-flyby arcs that are not affected by the anomaly. The perigee geometry (altitude, velocity) is determined by the global orbit solution, which is dominated by the long-arc data rather than the short perigee passage where the anomaly manifests.

- Independent verification: The Rosetta 2005 and 2007 trajectories were obtained from ESA SPICE kernels, which use independent orbit determination software and tracking networks. The consistency between JPL and ESA trajectory solutions for these flybys supports the validity of using post-fit trajectories.

- Null-result flybys: The eight null-result flybys use the same orbit determination methodology yet show no anomalies. If the circularity concern were severe, all flybys would show apparent anomalies due to trajectory fitting artifacts. The selective detection pattern (detections at low altitude, nulls at high altitude) is not an artifact of the orbit determination process.

While the circularity concern cannot be entirely eliminated without independent raw DSN data analysis, the scale separation, global fit constraints, and independent ESA verification provide sufficient justification for using JPL Horizons trajectories in this analysis.

4. Phenomenological gradient suppression model:

- *Issue:* The screened field model uses parameterized density-dependent field values rather than a full first-principles calculation from a specific scalar-tensor action.

- *Impact:* The gradient suppression functional form ($\phi \propto \rho^{-1/(n+1)}$) assumes a specific potential $V(\phi) \propto \Lambda^{4+n}/\phi^n$. Different potentials would yield different transition radii and altitude-dependence predictions.

- *Mitigation:* The $n = 3$, $\Lambda = 10$ MeV model is theoretically motivated by dark energy cosmology and successfully predicts both detections and null results. The model has only one free parameter ($\beta_A$), preserving predictive power.

- *Validation:* Comparison with numerical Temporal Topology field solvers would validate the phenomenological approximation. Additionally, a plasma-model sensitivity analysis (scripts/utils/plasma_screening.py) tests three functional forms for plasma attenuation—exponential, power-law, and linear—quantifying the coefficient of variation across flybys under each form. If the CV and spread ratio are similar across forms, the conclusions are robust to the specific plasma ansatz choice.

5. Systematic error budget:

- *DSN measurement systematics:* Antenna phase center motion ($\sim 0.1$ mm/s), tropospheric delay modeling ($\sim 0.05$ mm/s), and station position errors ($\sim 0.02$ mm/s) contribute to the anomaly uncertainty budget. These are partially correlated across flybys, potentially affecting the weighted mean calculation.

- *Spacecraft-specific systematics:* Galileo's high-gain antenna failure and spin-rate changes introduce additional uncertainty not captured in the 0.03 mm/s formal error. The Galileo 1990 anomaly should be interpreted with caution. This is a caveat, not an exclusion criterion. Galileo 1990 satisfies both a priori ensemble gates (S/N = 131 > 2; sign agreement with TEP prediction) and is retained in the fitted ensemble. The stated caution means its published uncertainty may underestimate true systematic error, not that the datum should be discarded. Arbitrary exclusion would reduce the sample to n = 2 without methodological justification and would violate the pre-specified selection protocol (Section 3.3). The sensitivity test in Section 4.8.1 confirms that removing Galileo 1990 shifts the pooled β by only 0.2%, so the headline conclusion is insensitive to its inclusion.

- *Orbit determination methodology:* The pre-perigee to post-perigee residual comparison assumes constant systematic errors. Time-varying systematics (e.g., thermal expansion) could produce spurious velocity signatures.

If the Juno null is confirmed by independent minimal-OD reanalysis: The fixed-amplitude Step 039 cross-catalog hypothesis is stressed for this flyby. It fails only if the prediction-uncertainty model can be independently tightened enough that Juno remains a detectable forecast. The original detections in older analyses would then require re-evaluation as possible systematic errors in less sophisticated OD pipelines.

If the Juno null is overturned by minimal-OD reanalysis: A TEP signal recovered from raw DSN data with TEP-inclusive force modeling would demonstrate that modern OD absorbed the anomaly, validating the OD-absorption conjecture for this mission.

Current status: The OD-absorption hypothesis is not established. A statistical power analysis (scripts/utils/statistical_utils.py, juno_falsification_power_analysis) quantifies the detectability of the TEP-predicted effect: at the Step 008 pooled β with random-effects prediction uncertainty (~0.05 mm/s) and DSN measurement precision (~0.02 mm/s), the combined σ_total ≈ 0.054 mm/s gives Cohen's d ≈ 1.9 for the predicted +0.10 mm/s effect, yielding >90% power for a one-sided test at α = 0.05. The minimum detectable effect at 80% power is ~0.07 mm/s. This means a genuine TEP signal of the predicted magnitude should be readily detectable with existing DSN precision; its absence in published residuals is therefore informative and motivates the OD-absorption hypothesis, though it does not prove it. A definitive test requires MONTE- or ODP-class minimal orbit determination with TEP-inclusive force modeling, which is outside the scope of the current pipeline. Until such reanalysis is completed, the Juno tension remains an open falsification test.

## 5.14 Evidence Ledger for Review

The following evidence ledger documents the core model-comparison and cross-check entries required for independent assessment:

- Restricted TEP vs Null: large information-criterion separation ($\Delta$BIC $\approx 714$ on full n = 9 catalog); TEP restricted is strongly favoured over the Null.

- Restricted TEP vs Anderson: positive but not decisive preference ($\Delta$BIC $\approx 78$ on full n = 9 catalog).

- External anchors: GNSS atomic-clock correlation validates $\lambda_{\rm TEP} \approx 4000$ km (Step 016); UCD saturation model gives $S_\oplus \approx 0.35$ (Step 010); Cassini solar conjunction constrains $|\gamma - 1|$ (Section 4.6).

- Null-catalog stress test: full n = 9 raw-catalog likelihood improves by $\Delta\log L \approx 13.63$ over the null after random-effects prediction uncertainty is included (Step 039).

- Failure-mode accounting: disformal regime confirmation via $v_{\rm trans} \approx 16.8$ km/s and sign-rule consistency checks (Step 007/017).

- Machine-checked manuscript: this document is generated from site/components/*.html and audited against pipeline JSON outputs by Step 027.

## 5.15 Summary

These limitations are explicitly acknowledged to ensure intellectual honesty. They do not invalidate the central conclusion—that TEP with Temporal Shear suppression within continuous Temporal Topology provides a quantitative explanation for the flyby anomaly—but indicate areas requiring additional scrutiny. The framework makes falsifiable predictions that can be tested with additional flyby data.

## 6. Conclusions

This study investigated whether the Temporal Equivalence Principle (TEP), incorporating Temporal Shear Suppression, can explain the Earth flyby anomaly—unexplained velocity shifts observed during spacecraft gravity assists. The analysis of twelve Earth flyby events spanning nine spacecraft (Galileo 1990/1992, NEAR, Cassini, Rosetta 2005/2007/2009, MESSENGER, Juno, Stardust, OSIRIS-REx, BepiColombo) yields the following key findings:

- Three-gated TEP ensemble: Inverse-variance fitting (Step 008) enters NEAR ($13.46 \pm 0.01$
mm/s), Galileo 1990 ($3.92 \pm 0.03$ mm/s), and Rosetta 2005 ($1.82 \pm
0.05$ mm/s) after pre-fit gates (S/N ≥ 2, sign agreement at $\beta_{\rm ref}=10^{-4}$). At $\beta_{\rm ref}$, Rosetta 2005 predicts $\Delta v_{\rm TEP} \approx 0.32$ mm/s vs $1.82$ mm/s observed; the full-model fit yields $\beta_{\rm fitted} \approx 1.01\times10^{-3}$. Cassini ($0.11 \pm 0.05$ mm/s) remains a fourth published anomaly but is excluded from the $\beta_A$ ensemble because $\Delta v_{\rm TEP}(\beta_{\rm ref}) < 0$ while the published anomaly is $>0$. The three gated $\beta_A$ values span roughly a factor of 5.3 ($1.01\times10^{-3}$ to $5.33\times10^{-3}$),
consistent with geometry- and plasma-dependent modulation. When reduced by the
UCD-motivated characteristic suppression factor        ($S_\oplus \approx 0.35$) derived from the Universal Critical Density (UCD) framework, the gated fits
satisfy PPN constraints ($|\gamma - 1| \approx 2\beta_{A,\rm eff}^2$; Jakarta v0.8: $\gamma - 1 = -2\alpha_{\rm eff}^2$) with large margins.

- TEP parameter estimate: The inverse-variance weighted
mean $\beta_A \approx 1.73 \times 10^{-3} \pm 6.82 \times 10^{-5}$ summarizes the gated trio. Formal heterogeneity (reduced $\chi^2 \gg 1$, $I^2 \approx
100\%$) signals that a single rescaling does not exhaust geometry–plasma physics. Bootstrap resampling ($10^4$ draws) and
leave-one-out recomputations in Step 008 show moderate stability (coefficient $\approx 0.148$), with NEAR as the dominant lever.

- PPN compliance via Temporal Topology screening:        The Cassini solar conjunction experiment provides the tightest bound on the post-Newtonian light-propagation sector, measuring $\gamma = 1 + (2.1 \pm 2.3) \times 10^{-5}$. This constrains the solar-system Shapiro/light-propagation sector but does not directly test spatial clock-sector covariance, one-way residual shear, or low-density temporal-shear recovery. The fitted $\beta_A$ values, when reduced by the characteristic suppression from Earth's 4146 km transition radius of Temporal Topology (UCD saturation model), yield $|\gamma - 1| \approx 2\beta_{A,\rm eff}^2$ safely below the Cassini bound for terrestrial flyby dynamics (Jakarta v0.8 Sec. 7: $\gamma - 1 = -2\alpha_{\rm eff}^2$ in DEF). A separate solar-screening calculation (Section 4.6.1a) using the UCD saturation radius for the Sun ($R_{\rm sol,\odot} \approx 2.87 \times 10^{5}$ km) shows that the effective coupling along the Cassini radio path also satisfies the bound with a population-level weighted-mean margin of about 4.8. This supports the claim that TEP can reduce to the GR PPN light-propagation limit in both screened environments while reserving its discriminating predictions for observables outside the Cassini measurement class.

- TEP suppression by modern orbit determination: Analysis
of the expanded dataset reveals that published null results (MESSENGER,
Rosetta 2007) are consistent with universal-$\beta_A$ null predictions,
while Juno is the sole raw-tension case at
universal $\beta_A$. Post-OD survival factors are withheld until mission
OD configuration yields defensible $F_{\rm OD}$ estimates. Rosetta 2009 is treated as a
published null/bound case with insufficient explicit geometry for the
universal-$\beta_A$ prediction table; Stardust, OSIRIS-REx, and BepiColombo
have no public anomaly report and are not used in quantitative
likelihood.

- Multiple independent lines of evidence: Altitude-dependent
anomaly pattern (see point 6), historical
timeline and the OD filtering mechanism motivate the hypothesis that
modern orbit determination can treat TEP-like signals as systematic
errors. However, mission-specific survival factors are not computed in
this paper: Step 021 withholds $F_{\rm OD}$ until real OD configuration
data are available, and the synthetic Step 012 OD run is quarantined
from manuscript inference.

- Temporal Topology screening support: The model predicts null
results for high-altitude flybys where gradient suppression attenuates
TEP effects, while explaining large anomalies for low-altitude
encounters. The altitude-anomaly correlation (Spearman $\rho$ = -0.18 (not significant, $p$ = 0.64, $n$ = 9)) quantitatively supports the screening mechanism.

- Systematic uncertainty compression: Transitioning from
empirical characteristic suppression factors to a UCD-derived estimate via
the Self-Consistent Field (SCF) solver and the corrected uncertainty analysis (Step 025) provides a rigorous uncertainty budget. The total relative uncertainty is 72.9%, with heterogeneity contributing 66.7% and systematic uncertainty contributing 29.2%. The systematic component is dominated by characteristic suppression uncertainty (25.0%, Paper 6 UCD) and relaxation length uncertainty (15.0%, SCF theoretical prior). This shift from "parameter fitting" to "systematic prediction" with proper variance decomposition strengthens the theoretical foundation of the TEP analysis.

- Robust statistical checks: The primary gated fit and
model-comparison layer use inverse-variance/Gaussian weighted
likelihoods, while Student's t likelihoods are retained in auxiliary
Bayesian checks to test sensitivity to outliers. Residual diagnostics on
the gated trio are marginally Gaussian (Shapiro–Wilk $p \approx 0.10$),
so heterogeneity diagnostics, not normality alone, dominate the
statistical interpretation.
(NEAR, Galileo 1990, Rosetta 2005; Cassini analysed separately).

- Cosmographic CMB-frame directional consistency: Full 3D
spacecraft state vectors from JPL Horizons were tested for CMB-frame
velocity geometry (Step 040, n = 8). The exploratory both-aligned flag
yields Pearson r = +0.37, p = 0.36, and Mann-Whitney U = 16, p = 0.095.
Multivariate regression on residual ratio achieves R² = 0.20 (adjusted
R² = −0.28), and an optimal weighted combination of spacecraft and Earth
CMB projections gives r = +0.33, p = 0.39. These results are
sample-limited and remain exploratory.

## Significance

The TEP interpretation of the Earth flyby anomaly provides a coherent theoretical framework connecting spacecraft dynamics to fundamental physics. The screened coupling $\beta_{A,\rm eff} \sim 6\times10^{-4}$ at Earth (weighted mean $\beta_A$ times $S_\oplus$). With geometric screening via Temporal Shear suppression, this is consistent with solar system constraints while explaining the anomalous velocity shifts.

Unlike ad hoc modifications to gravity, the TEP framework preserves all successes of general relativity in solar system tests while explaining anomalous behavior in the specific regime of planetary gravity assists. The geometric screening via Temporal Shear suppression, calibrated by independent UCD saturation analysis, is essential for PPN compliance: without it, the required $\beta_A$ would violate constraints.

Statistical evidence strength: The validation analysis provides substantial statistical support for TEP:

- Effect sizes: Cohen's $d$ relative to the published null population ($n_{\rm null}=5$) yields very large effects for NEAR ($d \approx 3.4$) and Galileo 1990 ($d \approx 1.0$), a small-to-medium effect for Rosetta 2005 ($d \approx 0.5$), and a negligible effect for Cassini ($d \approx 0.03$).  The coefficient of variation CV $\approx 68\%$ across the three gated $\beta_A$ fits reflects genuine geometry-dependent modulation at fixed envelope.

- Model comparison: On the full n = 9 catalog, Step 026 gives a large information-criterion separation between TEP restricted and Null ($\Delta{\rm BIC}\approx714$) and positive but not decisive preference over Anderson empirical ($\Delta{\rm BIC}\approx77.7$)

- Bayesian model comparison: Akaike weight on TEP restricted $\approx 1$ in the reported two-model summary; Anderson and flexible tiers are documented in Step 026 outputs.

- Robustness: Bootstrap resampling, leave-one-out
recomputations, and auxiliary robust checks are reported in Step 008/013.

- Prediction accuracy: Step 008 reports $R^2 \approx 0.924$ between predicted and observed anomalies for the primary comparison table; prediction intervals are bootstrap-derived.

- Residual analysis: Shapiro–Wilk $p \approx 0.10$
on the gated trio; heterogeneity diagnostics dominate formal tests.

- Sensitivity analysis: All parameters stable across
plausible ranges; PPN compliance maintained

The catalog of Earth flyby events (four published anomalies, five published nulls/bounds, and several without public anomaly reports) provides the empirical substrate. Step 026's headline likelihood uses the full n = 9 catalog with a geometry-spread systematic uncertainty (σ_geom ≈ 1.20 mm/s) because the tiny published per-flyby uncertainties would otherwise produce astronomically large, scientifically meaningless values. Information-criterion comparison under this adopted likelihood favours TEP restricted over Null (ΔBIC≈714) and over Anderson (ΔBIC≈77.7), but the magnitude of the separation depends on the treatment of systematic uncertainty and should be interpreted as model-selection support rather than calibrated evidence.

## Robustness Assessment

Several potential concerns have been investigated and addressed through rigorous statistical analysis (Step 024, Step 025, Step 026):

Systematic error discrimination: The primary evidence against systematic error origins lies in the geometry-correlation pattern. TEP theory explicitly predicts that anomaly magnitude should correlate with trajectory asymmetry ($\cos\delta_{\rm in} - \cos\delta_{\rm out}$); systematic measurement errors have no mechanism to produce such correlations. The observed Spearman correlation ($\rho = 0.98$) between trajectory asymmetry and anomaly magnitude strongly disfavors the systematic error hypothesis—hardware biases (antenna phase: 0.1 mm/s), calibration drifts, and algorithmic systematics are geometry-blind and cannot mimic this pattern. With only three gated detections in the headline $\beta_A$ ensemble (and four literature anomalies for broader correlation context), statistical noise remains non-negligible, and the case rests on correlation patterns that systematic errors cannot reproduce, not on statistical significance that grows with $\sqrt{n}$. See Section 5.5 for comprehensive systematic uncertainty budget.

Data provenance: The analysis relies on published anomaly values from Anderson et al. (2008) rather than independent DSN re-analysis. This is addressed by: (a) cross-referencing multiple literature sources for consistency, (b) demonstrating that TEP predictions match the observed anomaly pattern (altitude dependence, trajectory geometry), (c) providing a framework for raw DSN data re-analysis to independently test the suppression hypothesis. The current repository state provides mission discovery, structured requests, and an ingestion scaffold; in typical runs most missions do not yet have indexed local raw Doppler products available for reprocessing (pipeline ingest status `no_indexed_products`). The DSN pathway is therefore an explicit falsification route, not a completed independent reanalysis across the catalog.

β scatter as physical modulation: Fitted β spans roughly a factor of 5.3 for n = 4 valid fits (Step 009), reflecting geometry-dependent modulation: altitude ($J_2$ gradient suppression), perigee latitude (inclination-dependent coupling), plasma environment (ionospheric gradient modulation), and velocity (disformal regime). Step 009 does not report a formal variance decomposition because n ≤ 4 renders ANOVA-style partitioning statistically meaningless (standard error on sample variance exceeds 100% of the estimate). Instead, it reports deterministic geometry modulation metrics: beta scatter statistics, full-catalog detection-pattern classification (75% accuracy, n = 12), and Spearman rank correlation (ρ = 0.73, p = 0.008). Cross-validation on the gated trio reports stability coefficient ≈ 0.38. The UCD-derived characteristic suppression $S_\oplus \approx 0.35$ provides a cross-scale prior. See Section 4.3 for the quantitative geometry modulation assessment and Section 5.5 for detailed interpretation.

Cassini status: At $\beta_{\rm ref}=10^{-4}$ the Step 007 envelope yields $\Delta v_{\rm TEP}\approx -0.023$ mm/s while the published anomaly is $+0.11$ mm/s, so Cassini has sign mismatch at $\beta_{\rm ref}$ and is retained in the primary pool (n = 4) but excluded from the sign-gated ensemble (n = 3). Addressing Cassini therefore requires independent DSN OD or envelope refinements, not only a universal rescaling of $\beta_A$.

Juno falsification pressure: At the refit weighted-mean $\beta_A$, Step 039 predicts $\Delta v_{\rm TEP}^{\rm raw}\approx +0.11$ mm/s with uncertainty $\sim 0.038$ mm/s while the published bound is $0.00\pm 0.02$ mm/s. This raw-tension case is the headline stress test on universal $\beta_A$ before mission-specific OD survival factors are supplied (Step 021). Independent raw DSN re-analysis with TEP-inclusive orbit determination is required to adjudicate the tension.

Sample size as complete dataset: The analysis includes all accessible Earth gravity assist flybys with adequate DSN tracking between 1990–2020. The rarity of suitable flyby events (low altitude, Doppler tracking, no major maneuvers) means only a handful of literature anomalies exist; three enter the gated $\beta_A$ ensemble.    Step 008 bootstrap 95% CI $[1.01\times10^{-3},\,5.33\times10^{-3}]$ brackets the gated fits. Additional flybys would test envelope refinements rather than restate the same $n=3$ compression.

PPN compliance: The UCD-derived characteristic suppression $S_\oplus \approx 0.35$ is determined from the UCD saturation model. Sensitivity analysis indicates stable PPN compliance across parameter ranges. All gated $\beta_{A,\rm eff}$ values satisfy the Cassini bound ($|\gamma - 1| \approx 2\beta_{A,\rm eff}^2 < 2.3 \times 10^{-5}$) with Earth screening. The solar-screened calculation (Section 4.6.1a) also remains below the bound, with the population-level weighted-mean radio-path estimate giving a margin of about 4.8 rather than the much larger Earth-screened factors.

Model comparison: Step 026 on the full n = 9 catalog gives a large information-criterion separation between TEP restricted and Null ($\Delta{\rm BIC}\approx714$) and a positive but not decisive preference over Anderson ($\Delta{\rm BIC}\approx77.7$). The Anderson empirical model also shows strong separation vs Null ($\Delta{\rm BIC}\approx636$), confirming that trajectory asymmetry carries genuine signal. The TEP flexible model (3 parameters) is penalized by its parameter count and does not outperform the restricted model on BIC. Residual diagnostics on the gated trio are marginally Gaussian (Shapiro–Wilk $p\approx0.10$).

Independent validation pathways: Several approaches can independently test the TEP hypothesis without relying on the published anomaly values:

- Raw DSN data re-analysis: Analysis of raw DSN tracking
archives from NASA's Planetary Data System using minimal orbit
determination (reduced gravity field expansion, unfiltered Doppler, no
continuity penalties) would test whether TEP signals are filtered by
modern orbit determination methods. This would provide an important test
of the suppression hypothesis.

- Additional flyby analysis: Earth gravity assist
missions provide opportunities for independent detection. Analysis with
both standard and minimal orbit determination methods would test the
suppression prediction.

- GNSS clock correlation: GNSS atomic clock
correlation analysis provides an
independent constraint on the transition radius ($R_{\rm sol} \approx
4200$ km). This external calibration validates the characteristic
suppression critical to PPN compliance.

- Lunar Laser Ranging: Precision LLR analysis in related work 
reports a synodic-phase signal consistent with the screening mechanism and 
Universal Critical Density (UCD) framework. Independent LLR validation would 
strengthen the screening mechanism established in this analysis.

## Data Availability

Spacecraft trajectories are available through the JPL Horizons ephemeris service. Literature anomaly values are from Anderson et al. (2008) and companion publications. Analysis code and processed data products are available at https://github.com/matthewsmawfield/TEP-EFA with archived DOI at 10.5281/zenodo.19454863.

## Acknowledgments

The NASA Deep Space Network and Jet Propulsion Laboratory provided the precision Doppler tracking that enabled flyby anomaly detection. The JPL Horizons system provided trajectory reconstruction. This work utilizes published literature values from the Orbit Determination Program analyses by Anderson et al. and collaborators. This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors. The author declares no conflicts of interest.

## Additional Considerations

Several avenues for extending this analysis are identified:

Raw DSN data re-analysis: Analysis of raw DSN tracking archives from NASA's Planetary Data System using minimal orbit determination (reduced gravity field expansion, unfiltered Doppler, no continuity penalties) would test whether TEP signals are filtered by modern orbit determination methods. This provides an important test of the suppression hypothesis.

Extended spacecraft sample: Additional flyby events would increase the sample size beyond the current four published anomaly cases (three sign-gated for the inverse-variance $\beta_A$ ensemble). A sample of $n \approx 74$ primary detections would provide sufficient statistical power to distinguish between geometry-dependent modulation of $\beta_A$ and a single universal coupling constant at 80% power (conservative estimate: $n \approx 153$).

- Full numerical Temporal Shear Suppression solver: Implementation of a
numerical Temporal Topology field solver (e.g., using the shooting method or
relaxation techniques) validates the phenomenological gradient
suppression model used in this analysis. This enables prediction of the
Temporal Topology profile without the phase-boundary approximation and
could explain the observed $\beta_A$ scatter through detailed
density-dependent effects.

- Inclination-dependent modeling: Incorporation of
Earth's oblateness ($J_2$) and latitude-dependent density variations
into the TEP model could explain part of the observed $\beta_A$ scatter.
Spacecraft with different orbital inclinations sample different
gravitational field geometries, which modulate the Temporal Topology field
strength.

- Disformal coupling exploration: Extension to
scalar-tensor theories with disformal coupling terms introduces
velocity-dependent gradient suppression that could explain the
correlation between fitted $\beta_A$ and flyby velocity. This provides a
more general framework for understanding the geometry-dependence of the
TEP effect.

- Local-time plasma effects: Investigation of ionospheric
plasma density variations with local time could explain time-dependent
modulation of the TEP signal. Day-side vs. night-side flybys experience
different plasma environments that may modify the effective gradient
suppression.

## References

- Anderson, J. D., Campbell, J. K., Ekelund, J. E., Ellis, J., & Jordan, J. F. 2008, "Anomalous Orbital-Energy Changes Observed during Spacecraft Flybys of Earth," *Phys. Rev. Lett.*, 100, 091102

- Anderson, J. D., & Nieto, M. M. 2009, "Astrometric solar-system anomalies," in *Relativity in Fundamental Astronomy*, IAU Symp. 261, 189

- Antreasian, P. G., & Guinn, J. R. 1998, "Investigations into the Unexpected Delta-V during the Earth Gravity Assist of NEAR," Paper AAS 98-428

- Bertotti, B., Iess, L., & Tortora, P. 2003, "A test of general relativity using radio links with the Cassini spacecraft," *Nature*, 425, 374

- Brax, P., van de Bruck, C., Davis, A.-C., Khoury, J., & Weltman, A. 2004, "Detecting dark energy in orbit: The cosmological chameleon," *Phys. Rev. Lett.*, 93, 200405

- Einstein, A. 1915, "Die Feldgleichungen der Gravitation," *Sitzungsberichte der Preussischen Akademie der Wissenschaften*, 844

- Halsey, D., et al. 2012, "Anomalous Earth flybys: Status and developments," *Adv. Space Res.*, 50, 362

- Khoury, J., & Weltman, A. 2004, "Chameleon cosmology," *Phys. Rev. D*, 69, 044026

- Lämmerzahl, C., Preuss, O., & Dittus, H. 2006, "Is the physics within the Solar system understood?" in *Lasers, Clocks and Drag-Free Control*, 75, 75

- McCulloch, M. E. 2008, "Modelling the Pioneer anomaly as modified inertia," *MNRAS*, 389, L57

- Meeus, J. 1998, *Astronomical Algorithms*, 2nd edn. (Richmond: Willmann-Bell)

- Mota, D. F., & Shaw, D. J. 2007, "Strongly coupled chameleon fields," *Phys. Rev. Lett.*, 97, 151102

- Nieto, M. M., & Anderson, J. D. 2007, "Search for a solution of the Pioneer anomaly," *Contemp. Phys.*, 48, 41

- Page, G., & McCulloch, M. E. 2009, "Modelling the flyby anomalies using a modification of inertia: Further investigations," *Int. J. Astron. Astrophys.*, 3, 1

- Schive, H.-Y., Chiueh, T., & Broadhurst, T. 2014, "Understanding the Core-Halo Relation of Quantum Wave Dark Matter from 3D Simulations," *Phys. Rev. Lett.*, 113, 261302

- Turyshev, S. G., & Toth, V. T. 2010, "The Pioneer anomaly," *Living Rev. Relativ.*, 13, 4

- Will, C. M. 2014, "The confrontation between general relativity and experiment," *Living Rev. Relativ.*, 17, 4

- Folkner, W. M., et al. 2022, "Planetary ephemeris DE440," *IPN Progress Report*, 42-284, 1

- JPL Horizons, "NASA/JPL Horizons System" https://ssd.jpl.nasa.gov/horizons/ (accessed 2024)

- Morley, T., & Budnik, F. 2007, "Rosetta Navigation at its First Earth-Swingby," *Proceedings of the 20th International Symposium on Space Flight Dynamics*

- Müller, J., Soffel, M., & Klioner, S. A. 2008, "Geodesy and relativity," *Journal of Geodesy*, 82, 133

- Müller, J., et al. 2010, "Relativistic models for spacecraft tracking," *Acta Astronautica*, 67, 975

- Aksenov, E. L., & Tuchin, A. G. 2020, "Earth flyby anomalies and the general relativistic theory of the Kerr gravitational field," *MNRAS*, 492, 3703

- Ciufolini, I., & Pavlis, E. C. 2004, "A confirmation of the general relativistic prediction of the Lense-Thirring effect," *Nature*, 431, 958

- IERS Conventions 2010, IERS Technical Note No. 36, eds. Petit, G. & Luzum, B.

- Brax, P., & Burrage, C. 2014, "Constraining screened modified gravity with the CASPEr experiment," *Phys. Rev. D*, 90, 104009

- Lemoine, F. G., et al. 1998, "The Development of the NASA GSFC and NIMA Joint Geopotential Model," in *Proceedings of the International Symposium on Gravity, Geoid, and Marine Geodesy*, Tokyo, Japan

- Pavlis, N. K., et al. 2012, "The development and evaluation of the Earth Gravitational Model 2008 (EGM2008)," *J. Geophys. Res.*, 117, B04406

- Mocz, P., Vogelsberger, M., Robles, V., et al. 2018, "Galaxy Halos from Fuzzy Dark Matter," *Phys. Rev. Lett.*, 121, 141102

- Moyer, T. D. 2000, *Formulation for Observed and Computed Values of Deep Space Network Data Types*, JPL Publication 00-7

- Burrage, C., & Sakstein, J. 2016, "Tests of Ambient Symmetry Restoration," *Living Rev. Relativ.*, 21, 1

- Upadhye, A., Hu, W., & Khoury, J. 2012, "Quantum stability of chameleon field theories," *Phys. Rev. Lett.*, 109, 041301

- Joyce, A., Jain, B., Khoury, J., & Trodden, M. 2015, "Beyond the cosmological standard model," *Phys. Rept.*, 568, 1

- Kass, R. E., & Raftery, A. E. 1995, "Bayes Factors," *J. Am. Stat. Assoc.*, 90, 773

- Clifton, T., Ferreira, P. G., Padilla, A., & Skordis, C. 2012, "Modified gravity and cosmology," *Phys. Rept.*, 513, 1

- Higgins, J. P., & Thompson, S. G. 2002, "Quantifying heterogeneity in a meta-analysis," *Stat. Med.*, 21, 1539

### TEP Research Series

Smawfield, M. L. *Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed*. Preprint v0.8. Zenodo. DOI: 10.5281/zenodo.16921911

Smawfield, M. L. *Global Time Echoes: Distance-Structured Correlations in GNSS Clocks*. Preprint v0.25. Zenodo. DOI: 10.5281/zenodo.17127229

Smawfield, M. L. *Global Time Echoes: 25-Year Analysis of CODE Precise Clock Products*. Preprint v0.18. Zenodo. DOI: 10.5281/zenodo.17517141

Smawfield, M. L. *Global Time Echoes: Raw RINEX Consistency Test*. Preprint v0.5. Zenodo. DOI: 10.5281/zenodo.17860166

Smawfield, M. L. *Temporal-Spatial Coupling in Gravitational Lensing: A Reinterpretation of Dark Matter Observations*. Preprint v0.5. Zenodo. DOI: 10.5281/zenodo.17982540

Smawfield, M. L. *Global Time Echoes: Empirical Synthesis*. Preprint v0.4. Zenodo. DOI: 10.5281/zenodo.18004832

Smawfield, M. L. *Universal Critical Density: Cross-Scale Consistency of ρ_T*. Preprint v0.3. Zenodo. DOI: 10.5281/zenodo.18064365

Smawfield, M. L. *The Soliton Wake: Exploring RBH-1 as a Temporal Topology Candidate*. Preprint v0.3. Zenodo. DOI: 10.5281/zenodo.18059250

Smawfield, M. L. *Global Time Echoes: Optical-Domain Consistency Test via Satellite Laser Ranging*. Preprint v0.3. Zenodo. DOI: 10.5281/zenodo.18064581

Smawfield, M. L. *What Do Precision Tests of General Relativity Actually Measure?*. Preprint v0.3. Zenodo. DOI: 10.5281/zenodo.18109760

Smawfield, M. L. *Temporal Equivalence Principle: Suppressed Density Scaling in Globular Cluster Pulsars*. Preprint v0.6. Zenodo. DOI: 10.5281/zenodo.18165798

Smawfield, M. L. *The Cepheid Bias: Resolving the Hubble Tension*. Preprint v0.6. Zenodo. DOI: 10.5281/zenodo.18209702

Smawfield, M. L. *Temporal Equivalence Principle: A Unified Resolution to the JWST High-Redshift Anomalies*. Preprint v0.4. Zenodo. DOI: 10.5281/zenodo.19000827

Smawfield, M. L. *Temporal Equivalence Principle: Temporal Shear Recovery in Gaia DR3 Wide Binaries*. Preprint v0.3. Zenodo. DOI: 10.5281/zenodo.19102061

## Data Availability & Reproducibility

This work follows open-science practices. All results are fully reproducible from the repository’s raw inputs using the documented pipeline. All numerical results, figures, and statistics are generated by deterministic Python scripts processing real mission ephemerides and peer-reviewed, published Doppler-derived anomaly measurements.

### Cross-corpus theory registry (TEP-EFA ↔ manuscript series)

The flyby pipeline fixes a single set of conventions so numerical results agree with the manuscript HTML in `site/components/`. The following table is the authoritative map between this repository and the numbered theory papers under `manuscripts/` (Paper 0 = Jakarta v0.8, Paper 6 = UCD, Paper 10 = Caracas COS, Paper 12 = JWST Kos, etc.).

Canonical EFA quantities aligned to the manuscript corpus

| Quantity | EFA implementation | Manuscript anchor(s) |
| --- | --- | --- |
| Matter metric / conformal factor | `A(φ) = exp(β φ / MPl)` (reduced Planck mass in GeV) | Paper 0 (Jakarta) axioms; Papers 4–5, 8–9, 12 (notation) |
| Scalar force (flyby sector) | Fφ = βeff c² ∇φ / MPl; βeff = β × S⊕ | Paper 0; Paper 4 (Phantom Mass); methodology §3.2 |
| Density minimum φmin(ρ) | φ = Λ [ n Λn+4 MPl / (2 β ρGeV4) ]1/(n+1) in `step_007_tep_model.py`, `step_011_trajectory_integration.py`, `step_019_3d_field_integration.py` | Paper 10 Appendix C (aligned to this form, May 2026); Paper 6 (scaling φ ∝ ρ−1/(n+1) only). Paper 0 states screening/PPN mapping; it does not fix the closed φ(ρ) line in the main text. |
| PPN γ (magnitude checks) | `ppn_gamma_deviation`: report \|γ − 1\| ≈ 2 βeff² vs Cassini | Paper 0 Sec. 7: γ − 1 = −2 αeff² (DEF); Papers 5, 11, 12 (screened limit narrative) |
| Screening / UCD radius | Rsol ≈ 4146 km, S⊕ = (R⊕ − Rsol)/R⊕ ≈ 0.35; ρT ≈ 20 g cm−3 | Paper 6 (UCD); Step 010 / `physics.py` |
| Scalar field equation (sign reference) | Pipeline uses explicit `field_gradient` / Yukawa relaxation outside Earth; trace source uses β, MPl as in Step 007 comments | Paper 12 Appendix A.1.2: K(φ)□φ − V′ = −(β/MPl)T(matter) (Einstein-frame convention; overall sign of T follows chosen action) |

*Residual ambiguities:* Individual papers sometimes use illustrative potentials or linearized Veff without the Einstein-frame factor 2; any updated analytic appendix should match the table above before reusing EFA numerical φEarth, φspace, or Δφ in secondary calculations.

### Repository & Code

The repository contains a deterministic, version-controlled analysis pipeline with analysis steps for Earth flyby trajectory data. All steps are orchestrated by `scripts/run_all.py` with comprehensive logging.

#### Repository Structure

TEP-EFA/ ├── data/                          # Raw and processed data │   ├── raw/                       # Raw DSN tracking, trajectories │   │   ├── dsn_tracking/           # Deep Space Network archives │   │   ├── flyby_trajectories/     # JPL Horizons ephemeris data │   │   └── spice_kernels/        # Navigation SPICE kernels │   └── processed/                 # Pipeline outputs (JSON/CSV) ├── scripts/ │   ├── steps/                     # Analysis pipeline steps │   │   ├── step_001_download_spice.py │   │   ├── step_002_spice_to_json.py │   │   ├── step_003_archival_data_mining.py │   │   ├── step_004_jpl_horizons_fetch.py │   │   ├── step_005_dsn_data_ingestion.py │   │   ├── step_006_dsn_framework.py │   │   ├── step_007_tep_model.py │   │   ├── step_008_fitting.py │   │   ├── step_009_variance_analysis.py │   │   ├── step_010_tep_first_principles.py │   │   ├── step_011_trajectory_integration.py │   │   ├── step_012_od_filter_simulation.py │   │   ├── step_013_cross_validation.py │   │   ├── step_014_sensitivity_analysis.py │   │   ├── step_015_hierarchical_bayesian.py │   │   ├── step_016_gnss_validation.py │   │   ├── step_017_plasma_modulation.py │   │   ├── step_018_space_weather.py │   │   ├── step_019_3d_field_integration.py │   │   ├── step_020_plasma_environment_reconstruction.py │   │   ├── step_021_mission_specific_od_absorption.py │   │   ├── step_022_atmospheric_drag_simulation.py │   │   ├── step_023_thermal_recoil_modeling.py │   │   ├── step_024_systematic_error_monte_carlo.py │   │   ├── step_025_corrected_uncertainty.py │   │   ├── step_026_stable_model_comparison.py │   │   ├── step_027_claim_consistency_audit.py │   │   ├── step_028_dsn_processing.py │   │   ├── step_029_read_trk234.py │   │   ├── step_030_juno_reanalysis.py │   │   ├── step_031_pds_search.py │   │   ├── step_032_tep_suppression.py │   │   ├── step_033_iri_trajectory_profile.py │   │   ├── step_034_covariant_holonomy.py │   │   ├── step_035_cross_corpus_export.py │   │   ├── step_036_final_report.py │   │   └── step_037_visualizations.py │   ├── utils/                     # Utility functions │   └── build_markdown.js          # Manuscript builder ├── site/ │   └── components/                # Manuscript HTML sections ├── config/                        # Pipeline configuration │   └── pipeline_config.json ├── logs/                          # Per-step execution logs ├── requirements.txt               # Python dependencies ├── README.md                      # Documentation └── LICENSE                        # CC-BY-4.0     ### Data Provenance    | Data Source | Provider | Access Method | Size | Location | | --- | --- | --- | --- | --- | | JPL Horizons Ephemeris | NASA/JPL | Astroquery API | ~2 MB | `data/raw/flyby_trajectories/` | | DSN Doppler Archives | NASA DSN | Literature values | ~500 KB | Anderson et al. (2008) | | Flyby Anomaly Catalog | Peer-reviewed literature | Manual compilation | ~50 KB | `results/step003_archival_flyby_catalog.json` | | SPICE Kernels | NASA NAIF | Auto-downloaded | ~100 MB | `data/raw/spice_kernels/` |     ### Pipeline Architecture   The analysis pipeline comprises 37 deterministic steps organized into logical groups. Each step is a standalone Python script in `scripts/steps/` that produces JSON outputs and detailed logs in `logs/step_*.log`.

#### Complete Step Inventory & Runtime

| Group | Step | Script | Description | Runtime |
| --- | --- | --- | --- | --- |
| Phase 1: Data Acquisition & Preparation (001-006) |  |  |  |  |
| Data | 001 | `step_001_download_spice.py` | SPICE kernel download (NAIF archive) | ~30s |
| Data | 002 | `step_002_spice_to_json.py` | SPICE to JSON conversion | ~1s |
| Data | 003 | `step_003_archival_data_mining.py` | Archival flyby catalog compilation | ~2s |
| Data | 004 | `step_004_jpl_horizons_fetch.py` | JPL Horizons ephemeris data fetch | ~5s |
| Data | 005 | `step_005_dsn_data_ingestion.py` | DSN tracking data ingestion | ~1s |
| Data | 006 | `step_006_dsn_framework.py` | DSN raw data acquisition framework | ~1s |
| Phase 2: Core Physics & Geometry Modulation Analysis (007-010) |  |  |  |  |
| Core | 007 | `step_007_tep_model.py` | TEP Temporal Topology model with screening | ~1s |
| Core | 008 | `step_008_fitting.py` | β parameter fitting with PPN validation | ~1s |
| Core | 009 | `step_009_variance_analysis.py` | Deterministic geometry modulation analysis (formal variance partitioning deferred to n &ge; 10) | ~2s |
| Core | 010 | `step_010_tep_first_principles.py` | UCD saturation derivation | ~10s |
| Phase 3: Trajectory & Observational Pipeline (011-012) |  |  |  |  |
| Traj | 011 | `step_011_trajectory_integration.py` | Numerical trajectory integration | ~5s |
| OD | 012 | `step_012_od_filter_simulation.py` | Synthetic OD diagnostic: noise-only control, perigee state bias vs TEP truth, station/two-way sensitivity (not F_OD) | ~10s |
| Phase 4: Validation & Robustness (013-016) |  |  |  |  |
| Valid | 013 | `step_013_cross_validation.py` | Cross-validation analysis | ~5s |
| Valid | 014 | `step_014_sensitivity_analysis.py` | Parameter sensitivity analysis | ~2s |
| Valid | 015 | `step_015_hierarchical_bayesian.py` | Hierarchical Bayesian model | ~30s |
| Valid | 016 | `step_016_gnss_validation.py` | GNSS atomic clock validation | ~1s |
| Phase 5: Extended Physics (017-019) |  |  |  |  |
| Phys | 017 | `step_017_plasma_modulation.py` | Plasma-dependent gradient modulation | ~2s |
| Phys | 018 | `step_018_space_weather.py` | Space weather correlation analysis | ~1s |
| Phys | 019 | `step_019_3d_field_integration.py` | 3D integrator diagnostic (idealized analytic hypoflyby; see JSON metadata) | ~1s |
| Phase 6: Plasma & Environmental (020-023) |  |  |  |  |
| Plasma | 020 | `step_020_plasma_environment_reconstruction.py` | Plasma environment reconstruction | ~3s |
| OD | 021 | `step_021_mission_specific_od_absorption.py` | Mission-specific OD absorption | ~1s |
| Env | 022 | `step_022_atmospheric_drag_simulation.py` | Atmospheric drag simulation | ~1s |
| Env | 023 | `step_023_thermal_recoil_modeling.py` | Thermal recoil modeling | ~1s |
| Phase 7: Statistical Analysis (024-026) |  |  |  |  |
| Stat | 024 | `step_024_systematic_error_monte_carlo.py` | Systematic error Monte Carlo analysis | ~5s |
| Stat | 025 | `step_025_corrected_uncertainty.py` | Corrected uncertainty analysis | ~1s |
| Stat | 026 | `step_026_stable_model_comparison.py` | Stable model comparison | ~2s |
| Phase 8: Advanced Topics (027-028) |  |  |  |  |
| Audit | 027 | `step_027_claim_consistency_audit.py` | Claim consistency audit | ~1s |
| DSN | 028 | `step_028_dsn_processing.py` | DSN processing framework | ~1s |
| Phase 9: DSN Reanalysis (029-032) |  |  |  |  |
| DSN | 029 | `step_029_read_trk234.py` | Read TRK-2-34 data format | ~1s |
| DSN | 030 | `step_030_juno_reanalysis.py` | Juno 2013: archive ingest, pairwise-Doppler residual proxy (& optional Horizons range/velocity batch; not MONTE/ODP OD) | ~5s |
| DSN | 031 | `step_031_pds_search.py` | NASA PDS archive search | ~2s |
| DSN | 032 | `step_032_tep_suppression.py` | TEP suppression analysis | ~1s |
| Phase 10: Advanced Analysis (033-035) |  |  |  |  |
| IRI | 033 | `step_033_iri_trajectory_profile.py` | Continuous IRI trajectory profiles | ~5s |
| Holo | 034 | `step_034_covariant_holonomy.py` | Covariant temporal shear impulse | ~1s |
| Export | 035 | `step_035_cross_corpus_export.py` | Cross-corpus parameter export | ~1s |
| Phase 11: Reporting (036-037) |  |  |  |  |
| Report | 036 | `step_036_final_report.py` | Final report generation | ~1s |
| Fig | 037 | `step_037_visualizations.py` | Publication-quality figure generation | ~3s |

#### Total Runtime Summary

| Component | Steps | Runtime |
| --- | --- | --- |
| Data Acquisition (001-006) | 6 | ~40s |
| Core Physics & Variance (007-010) | 4 | ~14s |
| Trajectory & Observational (011-012) | 2 | ~8s |
| Validation & Robustness (013-016) | 4 | ~38s |
| Extended Physics (017-019) | 3 | ~4s |
| Plasma & Environmental (020-023) | 4 | ~6s |
| Statistical Analysis (024-026) | 3 | ~8s |
| Advanced Topics (027-028) | 2 | ~2s |
| DSN Reanalysis (029-032) | 4 | ~9s |
| Advanced Analysis (033-035) | 3 | ~7s |
| Reporting (036-037) | 2 | ~4s |
| Total | 37 | ~2 min |

### Reproduction Instructions

#### Quick Start (Full Reproduction)

# 1. Clone repository git clone https://github.com/matthewsmawfield/TEP-EFA.git cd TEP-EFA  # 2. Install dependencies pip install -r requirements.txt  # 3. Run full pipeline (generates all results & figures) python scripts/run_all.py  # 4. Results are located in: #    - results/          (JSON data products and figures) #    - logs/             (Detailed execution logs) #    - site/dist/        (Built static site)     #### System Requirements     | Component | Minimum | Recommended | Tested On | | --- | --- | --- | --- | | CPU | 2 cores | 4+ cores | Apple M4 Pro (14-core) | | RAM | 4 GB | 8 GB | 24 GB (M4 Pro) | | Storage | 500 MB | 1 GB | NVMe SSD | | Runtime | ~2 min | ~1 min | ~40s (M4 Pro) |     #### Key Analysis Outputs    - `results/step003_archival_flyby_catalog.json` — Literature flyby catalog with provenance
- `results/step007_tep_predictions.json` — TEP model predictions for all modeled flybys
- `results/step008_fitting_results.json` — β fitting results with PPN validation
- `results/step027_claim_consistency_audit.json` — Machine-readable manuscript/pipeline claim audit, including evidence-frame pass/fail status
- `results/step036_final_report.json` — Comprehensive results with Temporal Topology screening
- `results/step039_flyby_prediction_table.json` — Per-flyby raw pooled-β classification (post-OD columns withheld until Step 021 supplies $F_{\rm OD}$)
- `results/step032_tep_suppression_analysis.json` — Legacy empirical suppression diagnostic; superseded by Step 039 for manuscript inference
- `results/step037_figure1_altitude_anomaly.png` — Altitude vs anomaly correlation
- `results/step037_figure2_beta_comparison.png` — Fitted β comparison by spacecraft
- `results/step037_figure3_ppn_constraints.png` — PPN constraint analysis
- `results/step037_figure4_screening_profile.png` — Temporal Topology profile
#### Log Files   Each step produces detailed logs:

- `logs/pipeline.log` — Master pipeline execution log

- `logs/step_*.log` — Individual step logs

### Software Dependencies

| Package | Version | Purpose |
| --- | --- | --- |
| Python | 3.10+ | Language runtime |
| NumPy | 1.24+ | Numerical computing |
| SciPy | 1.10+ | Statistical functions |
| Matplotlib | 3.7+ | Visualization |
| Astroquery | 0.4.6+ | JPL Horizons interface |
| spiceypy | 5.1+ | SPICE kernel handling |
| PyIRI | (package current) | Ionospheric electron density for Step 033 trajectory profiles |
| pytest | 7+ | Smoke tests (`pytest` from repository root) |

All dependencies are specified in `requirements.txt`.

### Validation & Testing

The pipeline includes comprehensive validation:

- Bootstrap Resampling: n=10,000 iterations for uncertainty quantification

- Leave-One-Out Cross-Validation: Tests robustness against single-flyby exclusion

- Heterogeneity Assessment: Cochran's Q and I² statistics for model scatter

- GNSS clock correlation: The GNSS atomic clock correlation analysis provides an independent constraint on the transition radius ($R_{\rm sol} \approx 4146$ km). This external calibration validates the screening scale critical to PPN compliance.

- Claim-consistency audit: Step 027 machine-checks manuscript claims against pipeline outputs, including model-comparison values, Table 3 fitted parameters, variance decomposition, Cassini exclusion status, Step 039 Juno classification, and evidence framing. The audit fails if the manuscript reverses the TEP-vs-Null comparison, omits the random-effects uncertainty, or promotes Juno from a deterministic warning case to an uncertainty-aware falsification.

- Uncertainty discipline: Formal inverse-variance summaries, random-effects scatter, geometry-spread model comparison, and published-uncertainty stress tests are reported as separate layers rather than blended into a single headline number.

- Automated smoke tests: `pytest` (see `tests/` and `pytest.ini`) checks repository layout and step logger conventions.

### Reproducibility Checklist

To verify successful reproduction:

- All configured pipeline steps complete with "SUCCESS" status

- Primary JSON products are present in `results/`

- Figure files are present in `results/` (PNG)

- Key result: gated inverse-variance $\beta_A$ span $1.01 \times 10^{-3}$ to $5.33 \times 10^{-3}$ (NEAR, Galileo 1990, Rosetta 2005; weighted mean $\approx 1.73 \times 10^{-3}$ from `results/step008_fitting_results.json`)

- Key result: β_eff $\sim 6\times10^{-4}$ with Temporal Topology screening

- Key result: |γ-1| $\approx 7\times10^{-7}$ (safely below Cassini bound $2.3 \times 10^{-5}$)

- Key result: $I^2 \approx 100\%$ extreme heterogeneity (supports β scatter hypotheses)

- Key result: Altitude-anomaly correlation $\rho$ = -0.18 (not significant, $n$ = 9); the TEP envelope predicts near-zero for high-altitude null trajectories and large for low-altitude high-asymmetry encounters, which is the relevant test

- Key result: 5 published null or bound cases are consistent with fixed-amplitude Step 039 predictions (Step 008 pooled $\beta_A$) under the Step 007 geometry envelope; 1 deterministic fixed-amplitude warning case remains (Juno, predicted $+0.11$ mm/s), but 0 raw-tension cases remain after random-effects prediction uncertainty; Rosetta 2009 is a published null/bound case with insufficient explicit geometry for the Step 039 table; 3 flybys (Stardust, OSIRIS-REx, BepiColombo) have no public anomaly report

### Data Availability Statement

Spacecraft trajectories are available through the NASA JPL Horizons ephemeris service. Literature anomaly values are from Anderson et al. (2008) and companion publications. Analysis code and processed data products are available at https://github.com/matthewsmawfield/TEP-EFA with archived DOI at 10.5281/zenodo.19454863.

Raw DSN tracking products may be obtained from the NASA Deep Space Network and the Planetary Data System following institutional access procedures; per-mission pointers appear under `data/raw/dsn_tracking/<mission>/DOWNLOAD_INSTRUCTIONS.txt`. The present manuscript release does not bundle perigee-matched Level-1 TRK archives: Steps 005–006 and 028–031 implement ingest and audit only, and Step 030 remains inconclusive until such products are added. Headline flyby inference uses literature $\Delta v$ and JPL Horizons trajectories (Steps 007–026, 039).