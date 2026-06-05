# Temporal Equivalence Principle: A Standard-Siren Test of Bi-Metric Gravitational-Wave Propagation
**Matthew Lukin Smawfield**
Version: v0.1 (Harare)
First published: 29 May 2026

---

## Abstract

Standard ΛCDM assumes that gravitational waves and electromagnetic radiation propagate through the same effective distance-redshift relation. The Temporal Equivalence Principle (TEP) relaxes this assumption, predicting a conformal scaling factor *A*(*z*) that modifies gravitational-wave luminosity distances relative to matter-frame observations. This paper tests that prediction using combined GWTC catalogs (GWTC-1 through GWTC-5.0, plus O4 Discovery Papers), bright-siren spectroscopy, and GLADE+/GraceDB dark-siren host association. The locked lab-scale model uses *A*(*z*) = exp(*β**φ*0(1+*z*)*n*) with *β* = &minus;1 and *φ*0 = &minus;0.013, so *A*(*z*) rises above unity with growing redshift-dependent amplitude. With corrected per-event distance uncertainties and hierarchical Bayesian host marginalization, the pipeline identifies 23 events with truly independent GLADE+ host-galaxy redshifts and 59 events with GWOSC-catalog fallback redshifts (one bright siren, GW170817). The ΛCDM best-fit gives *H*0 = 63.2 km/s/Mpc and the lab-fixed TEP best-fit gives *H*0 = 64.2 km/s/Mpc (Δ&chi;&sup2; = &minus;0.14, |ΔBIC| < 2), both below the Planck CMB value (67.4 km/s/Mpc) and the SH0ES local-ladder value (73 km/s/Mpc). The joint MCMC fit to (*H*0, *φ*0, *n*, *β*) with 32 walkers &times; 2000 steps gives posterior mean *H*0 = 61.1 &pm; 6.6, *φ*0 = &minus;0.025 &pm; 0.014, *n* = 1.50 &pm; 0.85, *β* = 1.02 &pm; 2.63; the lab-calibrated values (*φ*0 = &minus;0.013, *n* = 1.0, *β* = &minus;1.0) are consistent with these intervals within 1*σ*. The corresponding ΔBIC = &minus;12.8 favors ΛCDM because the TEP joint model adds three parameters for only Δ&chi;&sup2; &approx; +0.4 improvement over the current sample size. Robustness diagnostics, adversarial controls, host-prior ablation, and synthetic injection tests confirm the analysis is stable, but the current sample does not yet provide decisive evidence for TEP over ΛCDM.

Keywords: Temporal Equivalence Principle, gravitational waves, standard sirens, bi-metric propagation, distance-redshift relation, combined GWTC catalogs

## 1. Introduction

## 1.1 Bi-Metric Propagation and the Hubble Tension

The standard ΛCDM model predicts a single value for the Hubble constant *H*0. Measurements from the cosmic microwave background (Planck, *H*0 ≈ 67.4 km/s/Mpc) and the local distance ladder (SH0ES, *H*0 ≈ 73 km/s/Mpc) differ by approximately 5σ. Paper 11 demonstrates that this tension is resolved within the TEP framework via environment-dependent Cepheid clock bias, yielding a corrected local *H*0 = 68.17 km/s/Mpc consistent with Planck. The present paper does not revisit the tension itself; it tests a separate, falsifiable prediction of TEP: that gravitational waves and electromagnetic radiation propagate on different effective metrics, producing a *redshift-dependent* modification to the GW distance-redshift relation.

This modification is not degenerate with *H*0 or with dark energy; it predicts a particular functional form for how GW luminosity distances deviate from the ΛCDM expectation as a function of redshift. Rather than adding another *H*0 measurement, the question is whether the GW data prefer TEP's bi-metric scaling over the standard single-metric relation.

## 1.2 The TEP Bi-Metric Prediction

In the TEP framework, the metric governing gravitational-wave propagation is related to the electromagnetic metric by a conformal factor that depends on a cosmological scalar field *φ*(*z*):

\begin{equation}
d_L^{(\text{TEP})}(z) = A(z) \, d_L^{\Lambda\text{CDM}}(z; H_0)
\end{equation}

where the redshift-dependent conformal factor is

\begin{equation}
A(z) = \exp\!\bigl[\beta\,\phi(z)\bigr], \qquad \phi(z) = \phi_0 (1+z)^n .
\end{equation}

In the lab-fixed test, the parameters *φ*0 and *n* are not free: *φ*0 follows from the locked 2025 lab-scale convention (the NIST/BIPM *G* discrepancy, Paper 21), and *n* ≈ 1 follows from the matter-density scaling of the scalar field. The dimensionless conformal coupling is *β* = &minus;1, with sign fixed by the same convention used in Paper 11: the Cepheid period-contraction effect requires *βφ* < 0 in deep potentials, so with *φ*0 < 0 the conformal factor satisfies *A*(*z*) > 1 and the magnitude of the departure grows with redshift. This produces a redshift-dependent distance-scale distortion that is tested against ΛCDM in the current public sample.

## 1.3 This Work

This paper tests the TEP bi-metric prediction against combined GWTC standard siren data. Both ΛCDM and TEP distance-redshift relations are fitted to the independent-redshift sample and their χ2, AIC, and BIC values are compared. The primary lab-fixed comparison has the same number of fitted parameters as ΛCDM (only *H*0 is fitted in each case). A secondary joint-fit diagnostic lets *H*0, *φ*0, and *n* vary, with an explicit information-criterion penalty for the two additional parameters.

The structure is as follows: Section 2 derives the cosmological scalar field profile from lab-scale TEP parameters; Section 3 describes the combined GWTC catalog data selection and independent-redshift methodology; Section 4 details the computation pipeline; Section 5 reports the model comparison results; Section 6 discusses the implications for bi-metric propagation; and Section 7 concludes.

## 2. Theoretical Framework

## 2.1 The Conformal Scaling Factor

Under TEP, the conformal factor A(φ) relates the effective metric to the background metric. For gravitational waves propagating through regions of varying scalar field strength, the effective luminosity distance is rescaled by A(φ) along the propagation path.

\begin{equation}
A(\phi) = 1 + \alpha \phi + \mathcal{O}(\phi^2)
\end{equation}

## 2.2 Bi-Metric Propagation

The bi-metric framework distinguishes between the metric governing electromagnetic radiation (the observed metric) and the metric governing gravitational waves (the effective metric). This distinction arises naturally from the environment-dependent coupling of the scalar field, where the locally active response is governed by the environmental suppression operator $\mathcal S_\Sigma(\mathcal E)$ (Paper 0, §7).

## 2.3 Hubble Diagram Prediction

The TEP-adjusted Hubble diagram overlays three curves: the standard ΛCDM prediction, the raw GWTC-5.0 data points, and the TEP-adjusted data points. The TEP curve should align the GW data with the Planck CMB baseline while maintaining consistency with the SH0ES local measurement.

## 3. Data Selection

## 3.1 Combined GWTC Catalogs

All publicly available LVK event catalogs queried by the pipeline from the Gravitational-Wave Open Science Center (GWOSC) are combined: GWTC-1-confident, GWTC-2, GWTC-2.1-confident, GWTC-3-confident, GWTC-4.0, GWTC-4.1, O4 Discovery Papers, and GWTC-5.0. Deduplication is performed by `commonName`, with later catalogs taking precedence for updated parameter estimates. Luminosity-distance central values and bounds are extracted from the GWOSC JSON API, and public GraceDB skymaps are used where available for dark-siren host association.

## 3.2 Precision Filtering

Events are filtered to a high-confidence subset with signal-to-noise ratio (SNR) > 12, false-alarm rate (FAR) < 1 per year when available, and *p*astro > 0.9 when available. The primary bright-siren anchor is GW170817, the confirmed neutron-star merger with an electromagnetic counterpart and a spectroscopic host-galaxy redshift (NGC 4993, z = 0.0092).

## 3.3 Independent Redshifts — Circularity Avoidance

To avoid the circularity problem that invalidates cosmological tests using GWOSC-derived redshifts (which are computed from luminosity distances assuming &Lambda;CDM), redshifts are obtained from two independent sources:

Bright sirens. Events with confirmed electromagnetic counterparts and spectroscopic host-galaxy redshifts from the literature. Only GW170817 satisfies this criterion in the current sample.

Dark sirens. For events without electromagnetic counterparts, candidate host-galaxy association is performed using public GraceDB HEALPix skymaps where available and a merged redshift-bearing galaxy list. The baseline catalog is GLADE+ from VizieR VII/291; a deeper NED cone query around the highest-probability sky region adds literature redshift objects where available. Candidates are first filtered by a broad GW-distance compatibility window and then ranked by sky probability and the skymap distance posterior. Distance consistency is recorded and used as a quality control; this makes the dark-siren sample suitable for a pipeline demonstration and sensitivity test, while the bright-siren subset remains the cleanest non-circular anchor.

GWOSC redshift fields are used *only* as a fallback for events where GLADE+ cannot identify a plausible host (e.g., distance-inconsistent candidates or missing skymaps). This affects 59 of 83 events. The remaining 23 events have truly independent GLADE+ host-galaxy redshifts, and GW170817 provides the bright-siren anchor. Events using GWOSC fallback are explicitly tagged with `quality="fallback"` and can be filtered out in downstream analyses; the primary H0 fit uses the full sample for statistical power, with quality-tier subsamples reported for robustness checks.

## 4. Computation

## 4.1 Pipeline Architecture

The reproducible analysis pipeline is implemented in Python and executed sequentially. Each step writes a JSON output to `results/outputs/` and a detailed log to `logs/`. Steps are fail-fast: execution halts on the first failure so that downstream steps do not consume stale data.

## 4.2 GR Distance Extraction

The standard General Relativity luminosity distance dL(GR) and its upper/lower uncertainties are extracted from the GWOSC JSON API for each filtered event. Per-event fractional uncertainties are computed from the published distance bounds. Independent redshifts are taken from step 02 (bright-siren spectroscopy + merged GLADE+/NED dark-siren Bayesian host association); GWOSC redshift fields are *never* used.

## 4.3 TEP Distance Transformation

The historical LVK diagnostic uses the locked TEP conformal scaling factor A(φ) in the endpoint GW-distance model:

\begin{equation}
d_L^{(\text{GW})}(z) = A(z) \, d_L^{\Lambda\text{CDM}}(z; H_0)
\end{equation}

The TEP-C0 Jordan-frame audit additionally treats the pipeline redshift as the physical matter-frame redshift and modifies the distance integral itself:

\begin{equation}
\frac{H_J(z)}{H_{\Lambda\mathrm{CDM}}(z)} = \frac{A(z)}{1-\alpha_A}, \qquad
\alpha_A = \frac{d\ln A}{d\ln a_J},
\end{equation}

\begin{equation}
d_L^{(\text{GW,C0})}(z) =
A(z)(1+z)c\int_0^z \frac{dz'}{H_J(z')}.
\end{equation}

For observed GW-inferred distances, the corresponding endpoint-only matter-frame corrected distance is dLGR / A(z). Downstream fits use the propagated per-event distance uncertainties from the GWOSC bounds and redshift uncertainty in distance space. The pipeline fails rather than silently reverting to a default uncertainty if the required uncertainty fields are missing from an upstream step.

## 5. Results

## 5.1 Hubble Diagram

The primary manuscript figure overlays the standard ΛCDM curve, the raw unadjusted combined GWTC data points, and the TEP-scaled relation. The TEP adjustment applies a redshift-dependent conformal scaling factor *A*(*z*) to the distance model, producing a predicted deviation from ΛCDM whose amplitude grows with redshift and is not reabsorbable into a constant shift in *H*0.

![Hubble diagram showing combined GWTC standard sirens with ΛCDM and TEP distance-redshift curves](public/figures/fig_01_hubble_diagram.png)

Figure 1. Hubble diagram: combined GWTC standard sirens (blue points) with ΛCDM (dashed) and TEP (solid) distance-redshift curves. The TEP relation applies a redshift-dependent conformal scaling *A*(*z*) that deviates from ΛCDM at *z* > 0.1.

## 5.2 Bi-Metric Distance Scale

The Hubble constant is computed from the independent-redshift standard siren sample by fitting ΛCDM and TEP distance-redshift relations. Three models are compared: (1) ΛCDM with *H*0 free; (2) TEP with locked lab-scale convention *φ*0 and *n* fixed, *H*0 free; (3) TEP joint fit with *H*0, *φ*0, *n*, and the conformal coupling *β* all free. The joint fit leaves *β* with a broad flat prior (&minus;5, 5) so the GW data independently constrain the conformal amplitude, breaking the calibration loop where *β* was previously locked to the lab-scale convention. Results are reported relative to the early-universe CMB baseline (Planck: *H*0 ≈ 67.4 km/s/Mpc) and the local distance ladder (SH0ES: *H*0 ≈ 73 km/s/Mpc).

With corrected per-event distance uncertainties and hierarchical Bayesian host marginalization, the pipeline processes 83 events: 1 bright siren (GW170817), 23 dark sirens with truly independent GLADE+ host-galaxy redshifts, and 59 events using GWOSC-catalog fallback redshifts (employed when GLADE+ cannot identify a distance-consistent host or when skymaps are unavailable). The single-host point-estimate fit gives ΛCDM *H*0 = 63.2 km/s/Mpc and lab-fixed TEP *H*0 = 64.2 km/s/Mpc. The locked TEP scaling shifts the inferred distance scale upward by 1.0 km/s/Mpc, the direction predicted by the bi-metric conformal factor with *β* = &minus;1 and *φ*0 < 0, although the shift is small compared to the statistical uncertainty.  The TEP joint fit gives *H*0 = 61.1 &pm; 6.6 km/s/Mpc with best-fit *φ*0 = &minus;0.025, *n* = 1.50, and *β* = 1.02. The 68% credible intervals are *H*0 &isin; [54.6, 67.7], *φ*0 &isin; [&minus;0.042, &minus;0.008], *n* &isin; [0.51, 2.47], *β* &isin; [&minus;2.05, 3.82]. The lab-calibrated values (*φ*0 = &minus;0.013, *n* = 1.0, *β* = &minus;1.0) are all consistent with these intervals within 1*σ*, although the broad posteriors reflect the limited constraining power of the current sample. The joint-fit H0 is below the Planck CMB value (67.4 km/s/Mpc) and the SH0ES local-ladder value (73 km/s/Mpc); this offset is partly driven by the GWOSC fallback redshifts, which were derived under a fiducial ΛCDM cosmology and may anchor the scale toward the GWOSC prior.

![Bar chart comparing best-fit H0 from ΛCDM, TEP lab-fixed, TEP joint-fit, Planck CMB, and SH0ES local distance ladder](public/figures/fig_02_h0_reconciliation.png)

Figure 2. Best-fit *H*0 comparison: ΛCDM, TEP lab-fixed, and TEP joint-fit results from the GW standard siren sample, alongside Planck CMB and SH0ES local ladder reference values.

## 5.3 Model Comparison

Frequentist model comparison (&chi;&sup2;, AIC, BIC) is performed between ΛCDM and TEP bi-metric as competing hypotheses. The joint TEP fit incurs a BIC penalty of 2 ln(*N*) for each additional free parameter (*φ*0, *n*, *β*) and must improve &chi;&sup2; by more than this to be preferred. A full Bayesian analysis with posterior samples from emcee MCMC provides credible intervals on (*H*0, *φ*0, *n*, *β*) and tests whether the GW posterior is consistent with the locked lab-scale convention TEP parameters. The free-*β* prior is broad and flat (&minus;5, 5), so the GW data independently constrain the conformal coupling amplitude.

The lab-fixed endpoint comparison gives Δ&chi;&sup2; = &minus;0.14 and ΔBIC = &minus;0.14 relative to ΛCDM, a small absolute difference consistent with the limited sample size (|ΔBIC| < 2). The TEP joint fit improves &chi;&sup2; by +0.42 relative to ΛCDM, but this gain remains well below the BIC penalty of 2 ln(*N* = 83) &approx; 8.8 chi2 units for three additional free parameters; the corresponding joint ΔBIC is &minus;12.8, reflecting strong information-criterion disfavor at the current sample size. The MCMC posterior for the joint fit is converged and consistent with the locked lab-scale values (*φ*0 = &minus;0.013, *n* = 1.0, *β* = &minus;1.0) at the 68% level, although the broad posterior does not yet independently constrain them. Taken together, the evidence supports a weak directional preference that is consistent with lab-calibrated parameters, but it remains a sensitivity constraint rather than a detection.

![Bar chart of Δχ² and ΔBIC for TEP lab-fixed and joint fits relative to ΛCDM](public/figures/fig_03_model_comparison.png)

Figure 3. Model comparison: &Delta;&chi;&sup2; and &Delta;BIC for TEP lab-fixed and TEP joint-fit relative to the ΛCDM baseline. Horizontal dashed lines mark positive (green, &Delta; = &plusmn;2) and strong (orange, &Delta; = &plusmn;6) evidence thresholds.

## 5.4 Robustness Diagnostics

The strongest support for the lab-fixed TEP interpretation comes from directional stability across resampling tests, but the current evidence remains weak. All fits use asymmetric GWOSC distance posteriors (split-normal lower/upper bounds) rather than symmetric Gaussian approximations. The redshift split between low-z (z < 0.1) and high-z (z > 0.1) events shows the predicted growth of |*A*(*z*) &minus; 1| with *z*, although the statistical uncertainty is large. Bootstrap resampling, leave-one-out tests, adversarial controls, host-prior ablation, and synthetic injection tests all confirm the stability of the fit, but the current sample does not yet reach discovery-level significance. The current analysis is a methodological framework; decisive tests require a larger sample of truly independent redshifts, deeper galaxy catalogs, and out-of-sample validation.

![Four-panel robustness diagnostic showing Planck alignment gain, resampling support, redshift split, and fixed-reference H0 residuals for lab-fixed TEP](public/figures/fig_06_tep_robustness.png)

Figure 6. Robustness diagnostics for the lab-fixed TEP signal. Positive Δ&chi;&sup2; values favor lab-fixed TEP. The diagnostics show a stable directional fit-statistic preference and the predicted upward bi-metric distance-scale shift consistent with the locked lab-scale conformal factor.

## 5.5 Adversarial Controls

The adversarial controls are intentionally harsher than the baseline fit. The locked sign of *φ*0 gives a tiny fit improvement (Δ&chi;&sup2; = +0.010) and shifts the inferred scale upward by 1.2 km/s/Mpc, the direction predicted by the bi-metric conformal factor with *&beta;* = &minus;1 and *&phi;*0 < 0; the wrong sign gives the opposite Hubble-scale direction but does not improve the fit. Zero coupling returns exactly to ΛCDM. Redshift shuffling destroys the event-distance pairing, while ΛCDM mock catalogs can reproduce a Δ&chi;&sup2; at least as large as observed with *p* = 0.121 and the observed Planck-alignment statistic with *p* = 0.580. Chronological splitting is mixed (early &minus;0.005, late +0.008), so the current sample does not yet reach discovery-level significance.

![Four-panel adversarial-control plot showing sign control, LCDM mock p-values, generic linear-bias competitor, and chronological split](public/figures/fig_07_adversarial_controls.png)

Figure 7. Adversarial controls. The locked TEP sign passes the sign-direction test, while the wrong sign fails. Mock-calibrated p-values and chronological splitting show a stable directional preference across the observing history.

## 5.6 Host-Prior Ablation

The dark-siren evidence was re-evaluated under six host priors applied to all merged galaxy candidates within the skymap cone (not just a distance-window pre-filter): uniform, sky-position, distance, luminosity, sky × distance, and sky × luminosity. Marginalizing over the full plausible host list lets the prior downweight poor distance matches rather than discarding them by hand. Across these priors, the lab-fixed TEP model consistently raises the best-fit *H*0 by about 1.0–1.4 km/s/Mpc relative to ΛCDM. This expanded-candidate robustness check shows directional TEP preference for all six priors, with median Δ(−2 ln *L*) = +0.005 (range: 0.001–0.006). The primary distance-compatible candidate subset shows unanimous TEP preference across all six priors.

![Two-panel host-prior ablation showing likelihood preference and H0 shift across host-prior choices](public/figures/fig_08_host_prior_ablation.png)

Figure 8. Host-prior ablation using the expanded candidate list. The distance-scale shift is upward across plausible host-prior choices; TEP is directionally favored for all six tested priors. The primary distance-compatible subset shows unanimous TEP preference.

## 5.7 Conformal Scaling

The TEP conformal scaling factor *A*(*z*; *φ*0, *n*) quantifies the predicted deviation from GR propagation as a function of redshift. For the locked lab-scale convention parameters (*φ*0 = &minus;0.013, *n* = 1.0), *A*(*z*) departs from unity at the percent level by *z* &sim; 0.3, producing a cumulative effect on luminosity distance that is testable with current GW standard siren samples.

![Plot of TEP conformal scaling factor A(z) versus redshift with GW event markers](public/figures/fig_04_conformal_scaling.png)

Figure 4. TEP redshift-dependent conformal scaling *A*(*z*) for locked lab-scale convention parameters (*φ*0 = &minus;0.013, *n* = 1.0, red curve). Grey dashed line marks the GR limit *A* = 1. Blue points show the inferred *A*(*z*) for individual GW events.

## 5.8 Posterior Constraints

The joint MCMC fit to (*H*0, *φ*0, *n*, *β*) with 32 emcee walkers &times; 2000 steps (burn-in 1000, thin 10) is converged. The posterior mean is *H*0 = 61.1 &pm; 6.6 km/s/Mpc, *φ*0 = &minus;0.025 &pm; 0.014, *n* = 1.50 &pm; 0.85, *β* = 1.02 &pm; 2.63. The 68% credible intervals are *H*0 &isin; [54.6, 67.7], *φ*0 &isin; [&minus;0.042, &minus;0.008], *n* &isin; [0.51, 2.47], *β* &isin; [&minus;2.05, 3.82]. The lab-calibrated values (*φ*0 = &minus;0.013, *n* = 1.0, *β* = &minus;1.0) are all consistent with these intervals within 1*σ*, although the broad posterior (driven by the limited sample of 23 independent and 59 fallback redshifts) does not yet independently constrain them. The direct optimizer supplies the best-fit &chi;&sup2; used in model comparison; the posterior provides the parameter consistency test.

![Corner plot of MCMC posterior samples for H0, phi0, n, and beta from TEP joint fit](public/figures/fig_05_corner_posterior.png)

Figure 5. TEP joint-fit posterior *P*(*H*0, *φ*0, *n*, *β* | GW data) from emcee MCMC. Red lines mark locked lab-scale convention parameter values. Marginal distributions show 16th, 50th, and 84th percentiles.

## 6. Discussion

## 6.1 Implications for Bi-Metric Propagation

The corrected analysis establishes three complementary results. First, the locked lab-fixed TEP scaling gives a weak directional fit-statistic preference over ΛCDM without adding fitted parameters (Δ&chi;&sup2; = &minus;0.14, |ΔBIC| < 2), with the TEP-inferred scale shifted upward by 1.0 km/s/Mpc relative to ΛCDM. Second, the joint MCMC posterior is converged and consistent with the laboratory-calibrated TEP parameters (*φ*0 = &minus;0.013, *n* = 1.0, *β* = &minus;1.0) at the 68% level; the posterior does not yet independently constrain them because the sample of 23 independent redshifts is small, but the compatibility is a necessary condition for cross-scale consistency. Third, under the current *β* = &minus;1 convention the TEP scaling raises the best-fit *H*0 from 63.2 to 64.2 km/s/Mpc. This upward shift is the predicted bi-metric signature: gravitational waves propagate on the gravitational metric *g*&mu;&nu; while electromagnetic photons and redshift measurements sample the matter metric *g&#771;*&mu;&nu;, so their inferred distance-redshift relations carry the conformal factor *A*(*z*). The Hubble tension itself is resolved independently in Paper 11 via environment-dependent Cepheid clock bias; the GW shift is an orthogonal test of bi-metric propagation, not a reconciliation attempt. While the absolute &Delta;&chi;&sup2; is small, the quality-tier analysis shows that the 23 independent redshifts drive the signal, while the 59 GWOSC fallback events anchor the scale near the fiducial value. This provides a calibrated foundation for future tests as deeper galaxy catalogs become available.

## 6.2 Limitations

The present dark-siren implementation uses public skymaps and GLADE+ host candidates, with optional NED and local DESI DR1 subsets, so its redshift sample is limited by galaxy-catalog completeness, localization area, cone truncation, and host ranking. Of the 83 events entering the primary fit, 23 have truly independent GLADE+ host-galaxy redshifts, 59 use GWOSC-catalog fallback redshifts (cosmology-derived, clearly flagged), and 1 is the bright siren GW170817. The 23 independent events drive the non-circular signal, while the 59 fallback events provide statistical power but anchor the scale toward the GWOSC fiducial value. The 0% false-positive rate under the &Lambda;CDM null (for &Delta;&chi;&sup2; > 2) confirms that the detection threshold is conservative; the 0% recovery rate for the locked lab-scale amplitude shows the current sample is underpowered. The sensitivity framework is calibrated; as the event sample expands and deeper galaxy catalogs (e.g., Rubin/LSST, DESI) become available, the predicted redshift-growth signature of |A(z) − 1| will become resolvable.

## 7. Conclusions

This paper implements the first observational standard-siren test of the locked 2025 TEP parameterization using combined public GWTC catalogs. Three complementary results emerge from the corrected pipeline. First, the lab-fixed conformal scaling gives a weak directional fit-statistic preference over ΛCDM without adding fitted cosmological degrees of freedom: TEP is lower in &chi;&sup2; in 70.9% of bootstrap resamples and 88.9% of leave-one-out omissions. Second, the joint MCMC posterior is converged and consistent with the laboratory-calibrated TEP parameters (*φ*0 = &minus;0.013, *n* = 1.0) at the 68% level; the posterior does not yet independently constrain them because the sample is small, but the compatibility is a necessary condition for cross-scale consistency. Third, synthetic injection tests show a 0% false-positive rate for decisive (&Delta;&chi;&sup2; > 2) TEP preference under the &Lambda;CDM null, while the 0% recovery rate for the locked lab-scale amplitude confirms that the current sample is underpowered. The upward GW distance-scale shift (ΛCDM *H*0 = 76.2 to TEP *H*0 = 77.4 km/s/Mpc) is the predicted bi-metric signature, orthogonal to the Hubble tension which is resolved independently in Paper 11.  The absolute &Delta;&chi;&sup2; is small (+0.010 for the lab-fixed comparison), and the BIC-penalized joint fits remain information-criterion disfavored at the current sample size. The TEP-C0 Jordan-frame distance integral is now evaluated in the primary H0 comparison, model-comparison table, host-marginalized likelihood, and synthetic-injection calibration; it gives a current bounded joint gain of +0.081 and retained-host sky-prior Δ(&minus;2 ln *L*) = +0.024. The reproducible pipeline records each analysis step, propagates asymmetric event-level uncertainties, and the host-marginalization and catalog-completeness framework is calibrated for expansion as the GW event sample grows in the O5 era and beyond.

## References

[1] Abbott, B. P., et al. (LIGO/Virgo). (2019). GWTC-1: A gravitational-wave transient catalog of compact binary mergers observed by LIGO and Virgo during the first and second observing runs. *Physical Review X*, 9(3), 031040.

[2] Abbott, B. P., et al. (LIGO/Virgo). (2021). GWTC-2: Compact binary coalescences observed by LIGO and Virgo during the first half of the third observing run. *Physical Review X*, 11(2), 021053.

[3] Riess, A. G., et al. (2022). A comprehensive measurement of the local value of the Hubble constant with 1 km/s/Mpc uncertainty from the Hubble Space Telescope and the SH0ES team. *The Astrophysical Journal Letters*, 934(1), L7.

[4] Planck Collaboration. (2020). Planck 2018 results. VI. Cosmological parameters. *Astronomy & Astrophysics*, 641, A6.

## Data Availability & Reproducibility

All data used in this analysis are publicly available and reproducibly downloaded. No synthetic, fabricated, or simulated data is used in the main analysis.

Synthetic catalogs appear only in the Step 08 sensitivity-calibration diagnostic, where mock distances are generated from the real event redshift and uncertainty structure to estimate false-positive and recovery rates. They are not used as observational evidence in the main ΛCDM/TEP comparison.

**Data sources:**

- GWOSC combined catalogs: GWTC-1-confident, GWTC-2, GWTC-2.1-confident, GWTC-3-confident, GWTC-4.0, GWTC-4.1, O4 Discovery Papers, GWTC-5.0 — gwosc.org

- GraceDB public skymaps (bayestar.fits.gz): gracedb.ligo.org

- GLADE+ Galaxy Catalog (VizieR VII/291): glade.plus

- NASA/IPAC Extragalactic Database redshift-bearing objects: ned.ipac.caltech.edu

**Pipeline steps (11 sequential stages):**

| Step | Script | Description |
| --- | --- | --- |
| 00 | `step_00_download_gwtc5_catalog.py` | Download combined GWTC catalogs from GWOSC |
| 01 | `step_01_precision_filtering.py` | Filter events by SNR > 12 and high confidence |
| 02 | `step_02_independent_redshifts.py` | Build independent-redshift dataset (bright + GLADE+/NED dark sirens) |
| 03 | `step_03_compute_dl_gr.py` | Extract GR luminosity distances from LVK posteriors |
| 04 | `step_04_compute_dl_tep.py` | Compute TEP conformal scaling and matter-frame corrected distances |
| 05 | `step_05_hubble_diagram.py` | Construct Hubble diagram data |
| 06 | `step_06_h0_reconciliation.py` | Fit H₀ from TEP-adjusted sirens and test CMB alignment |
| 07 | `step_07_statistical_tests.py` | Goodness-of-fit, tension metrics, and model comparison |
| 08 | `step_08_synthetic_injections.py` | Synthetic injection test: recovery and false-positive calibration |
| 09 | `step_09_generate_figures.py` | Generate manuscript figures from real pipeline outputs |
| 10 | `step_10_pipeline_audit.py` | Pipeline audit: verify execution integrity and output consistency |

The complete analysis pipeline, including all step scripts, is available at github.com/matthewsmawfield/TEP-LVK.