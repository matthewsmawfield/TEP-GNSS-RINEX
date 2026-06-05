# Temporal Equivalence Principle: A Sidereal Modulation Audit of LHC Luminosity Data
**Matthew Lukin Smawfield**
Version: v0.1 (Hanoi)
First published: 28 May 2026

---

## Abstract

This paper audits CERN LHC luminosity data for sidereal modulation predicted by the Temporal Equivalence Principle (TEP). The pipeline analyzes uncalibrated CMS Run 2 luminosity-section data, fill-level multi-IP LHC SuperTable summaries, and optional LHCb/IP8 time-resolved Massi or NXCALS exports when available.

The current run does not establish proof of a TEP-LHC signal or a primary candidate. CMS/IP5 time-resolved and fill-level tests are mostly null: the primary decay-rate regression gives p=0.802, and many wrong-period controls exceed the sidereal term. After excluding zero-duration SuperTable rows, LHCb/IP8 peak luminosity is not significant (p=0.516). The public LPC Massi-derived cross-check also fails for LHCb peak luminosity (p=0.104). The strongest remaining lead is no configured primary analysis. The only nominal lead is the secondary CMS quality-filtered L0 test (p=0.027), which is not replicated by the primary lambda_fit test.

The analysis therefore frames TEP-LHC as a falsifiable candidate pathway rather than a confirmed detection. The next decisive target is proof-grade time-resolved LHCb/IP8 luminosity, with year, Run-era, wrong-period, delivered-luminosity, and independent-source replication gates frozen before interpretation.

Keywords: Temporal Equivalence Principle, LHC, luminosity decay, sidereal modulation, LHCb, CMS, claim audit, open data

## 1. Introduction

## 1.1 The Persistence of Environmental Systematics in Accelerators

The optimization of luminosity in large storage rings is constrained by beam degradation, environmental variation, and accelerator-configuration effects. In both the Large Hadron Collider (LHC) and its predecessor, the Large Electron-Positron Collider (LEP), operators have observed time-dependent effects that alter transverse beam dynamics and phase-space geometry. Standard non-linear optics models—such as recent 4D fixed-line studies in the CERN Super Proton Synchrotron (SPS)—primarily treat these effects as local Hamiltonian, machine, or environmental phenomena. This paper asks whether a residual sidereal component remains after conventional explanations and controls are applied.

## 1.2 The Limitations of Classical and Quantum Explanations

Time-correlated variations in particle accelerators are usually interpreted through mechanical terrestrial tides, operational covariates, environmental controls, or quantum-level Lorentz Invariance Violation (LIV) searches. TEP-LHC is framed here as an additional candidate hypothesis that must survive those explanations rather than replace them by assumption.

### Mechanical Earth Tides

The established explanation for many diurnal shifts in a 27-kilometer ring is geological. As the gravitational pull of the Sun and Moon stretches the Earth's crust, the accelerator tunnel deforms by roughly 1 millimeter. Because protons travel near the speed of light within a rigid radio-frequency (RF) bucket, this minute path-length alteration can shift the radial orbit and beam energy. Terrestrial tides therefore become a required control: any proposed sidereal residual must be separated from solar, lunar, thermal, operational, and fill-structure effects before it can be treated as physically meaningful.

### Lorentz Invariance Violation (LIV) and the SME

Conversely, the search for a cosmic vector—a spatial anisotropy aligned with the galactic background rather than the Earth—has been almost exclusively monopolized by the Standard-Model Extension (SME) pioneered by Alan Kostelecký. Experiments like MINOS, IceCube, and T2K have actively searched for a 23.93-hour sidereal modulation. Yet, these searches consistently return null results, placing ever-tighter constraints on LIV coefficients.

Those null results do not directly test the macroscopic beam-degradation channel proposed here. TEP-LHC therefore does not use SME limits as evidence for a signal; it uses them as a reminder that any new sidereal claim must be narrow, falsifiable, and heavily controlled.

## 1.3 The TEP Hypothesis

The Temporal Equivalence Principle (TEP) resolves this blind spot by proposing that proper time is not a passive geometric outcome, but a dynamical scalar field (ϕ) characterized by a background cosmic gradient, or Temporal Shear (∇ϕ). Crucially, TEP dictates that this field couples to matter through a proximity-based screening mechanism: the response depends on the geometric overlap of topological charge cores within a coherence volume, not on the ambient density of the medium.

The TEP hypothesis proposes that dilute probes and proximity-saturated hadronic matter could couple differently to a macroscopic temporal field. Proton bunches circulating in the LHC are therefore a candidate environment for testing this idea. Under TEP, as proximity-saturated protons traverse the LHC's 27-kilometer loop, they would sample a directional proper-time gradient. The Earth's rotation would then alter the ring's angle relative to this gradient every 23 hours and 56 minutes, creating a possible sidereal contribution to temporal decoherence and emittance growth.

The practical prediction is that, if the effect exists at observable scale, the uncalibrated decay curve of LHC instantaneous luminosity (L(t)=L₀ e^(−λt)) should contain a repeatable 23.93-hour residual that survives wrong-period, source, era, and delivered-luminosity controls.

## 1.4 This Work

This paper presents a deterministic forensic audit of uncalibrated CMS luminosity data from LHC Run 2, with fill-level cross-checks across LHC interaction points. By transforming strict UTC timestamps to Local Sidereal Time (LST), the analysis tests for a sidereal component of beam degradation and compares it against wrong-period, solar, source, and replication controls.

The structure of this paper is as follows: Section 2 presents the TEP theoretical framework, including the proximity saturation limit and Synchronization Holonomy mechanism; Section 3 describes the open-source forensic audit methodology; Section 4 reports the data analysis and results; Section 5 discusses the implications for reinterpreting accelerator literature; and Section 6 concludes with prospects for future colliders.

## 2. Theoretical Framework

## 2.1 The Proximity-Saturated Scalar Coupling

In standard General Relativity, proper time is a passive parameter—the integrated path length of a worldline through a pre-existing metric tensor. The Temporal Equivalence Principle (TEP) reformulates this by elevating proper time to a dynamical scalar field, ϕ(x^μ), which acts as an active physical medium.

A central axiom of TEP is that the coupling between matter and the temporal field is not strictly proportional to mass, but is governed by a non-linear proximity screening mechanism. The interaction strength is mediated by a coupling parameter κ(ξ) that asymptotically saturates at a critical proximity scale, observationally proxied by the density threshold:

\begin{equation}
\rho_c \approx 20~\text{g/cm}^3
\end{equation}

For macroscopic dilute objects (ξ < ξ_c), the coupling is weakly active, recovering standard weak-field gravitational approximations. Geometric core overlap suppresses observable Temporal Shear in the screened regime, but accelerator probes access an effective channel response κ_LHC that depends on beam energy, interaction topology, and the microscopic structure of the topological charge, not on the naive ambient density of the target.

The internal proximity of a single proton — the geometric overlap of its constituent quark topological charge cores within the confinement radius r_conf ~ 1 fm — yields a single-particle proximity parameter ξ_p ≈ (r_c/r_conf)^3 ~ 10^7–10^8. This is enormously larger than the astrophysical saturation scale ξ_c ≈ 10^–7 (corresponding to ρ_c ≈ 20 g/cm³). The hadron therefore sits deep in the single-particle topological regime where the conformal factor A(φ) is driven toward unity and Temporal Shear is suppressed at the core, regardless of the ambient vacuum density in the LHC beam pipe. The LHC beam probes the channel-specific response of the topological charge under high-energy momentum transfer, not a macroscopic density-screened coupling. The phenomenology is treated through an effective response coefficient κ_LHC, distinct from the many-body proximity screening function S(ξ) used in astrophysical and cosmological domains.

## 2.2 Temporal Shear and the Cosmic Vector

If the temporal field ϕ is globally isotropic, local proper time accumulation remains symmetric. However, TEP posits the existence of a background cosmic gradient—a topological "slope" to the temporal field. The local gradient of this field is defined as the Temporal Shear, Σ:

\begin{equation}
\Sigma_\mu = \nabla_\mu \ln A(\phi) = \beta \, \nabla_\mu \phi
\end{equation}

Because Σ is a directional vector anchored to a cosmic background (such as the galactic center or the rest frame of the cosmic microwave background), its projection onto the Earth is not static. As the Earth rotates, the physical orientation of the LHC ring—fixed in the lithosphere at Geneva—shifts relative to Σ.

The projection of the Temporal Shear onto the transverse plane of the accelerator modulates continuously, parameterized by the Earth's angular velocity ω_⊕ and the Local Sidereal Time angle θ_LST(t):

\begin{equation}
\Sigma_{\text{transverse}}(t) = |\Sigma| \cos(\omega_⊕ t - Φ_0) \hat{u}
\end{equation}

where ω_⊕ = 2π/23.934 hours is the sidereal frequency, and Φ_0 is the cosmic phase offset relative to the Geneva meridian.

## 2.3 Synchronization Holonomy and Beam Decoherence

The observable consequence of matter traversing a dynamical temporal field is phase-space decoherence. Under TEP, when a particle completes a closed orbit C through a fluctuating temporal landscape, it accumulates a residual proper-time phase shift. In the matter-metric 3+1 decomposition, the synchronization one-form $\tilde{\sigma}$ receives a conformal contribution proportional to $d\ln A$ that is exact and therefore integrates to zero around any closed loop. The disformal sector $B(\phi)$ introduces a non-exact correction $\delta\tilde{\sigma}$ whose exterior derivative does not vanish. The Synchronization Holonomy is the closed-loop integral of this non-exact disformal transport structure:

\begin{equation}
\Delta \phi_{\text{TEP}} = \oint_C \delta\tilde{\sigma} = \iint_S d(\delta\tilde{\sigma}),
\end{equation}

where $d(\delta\tilde{\sigma}) \neq 0$ because the disformal term $B(\phi)\nabla_\mu\phi\nabla_\nu\phi$ deforms the matter-frame simultaneity connection in a direction-dependent, non-integrable way. The LHC ring samples this non-exact structure as the Earth rotates, producing a sidereally modulated holonomy.

In the LHC, a proton bunch is not a point particle; it possesses a finite transverse emittance. Because the disformal correction to the synchronization connection exerts a local, direction-dependent deformation across the width of the beam pipe, protons at different transverse positions within the bunch accumulate proper time at microscopically different rates.

This differential temporal accumulation, integrated over millions of laps, forces the quantum wave packets of the individual protons to drift out of phase. This temporal decoherence manifests macroscopically as physical emittance growth ($\Delta\epsilon$).

Because instantaneous luminosity $L(t)$ is inversely proportional to the cross-sectional area of the beam ($\epsilon_x \epsilon_y$), the classical exponential decay rate $\lambda_0$ becomes perturbed by the sidereally oscillating holonomy:

\begin{equation}
L(t) = L_0 \exp[-(\lambda_0 + \lambda_{\text{TEP}} \cos(\omega_\oplus t - \Phi_0))t]
\end{equation}

It is this specific, unmodeled $\lambda_{\text{TEP}}$ oscillation that the forensic audit attempts to test against conventional beam, detector, and environmental explanations.

## 3. Methodology: Open-Source Forensic Audit

## 3.1 CMS Open Data and the brilcalc Infrastructure

To test for macroscopic temporal shear without relying on proprietary or heavily redacted accelerator metadata, this analysis utilizes the public data infrastructure of the Compact Muon Solenoid (CMS) experiment. The CMS Beam Radiation, Instrumentation, and Luminosity (BRIL) suite provides high-resolution instantaneous luminosity measurements for LHC Run 2 collisions.

The data is queried at the granularity of a Luminosity Section (LS). An LS is the fundamental time unit of LHC data-taking, defined precisely as 2^18 LHC orbits, which equates to approximately 23.3 seconds. This high-frequency binning provides a continuous, highly granular macroscopic proxy for beam degradation over standard 10-to-15-hour "STABLE BEAMS" operational fills, making it the ideal dataset for isolating long-wave temporal oscillations.

## 3.2 Bypassing Algorithmic Scrubbing

The critical vulnerability in standard searches for environmental systematics is the reliance on "physics quality" data. In standard CMS data-handling protocols, raw luminosity measurements are subjected to offline calibration files known as normtags. Because the standard model assumes spatial isotropy and invariant proper time, calibration algorithms inherently treat time-correlated anomalies—such as a 23.93-hour oscillation—as systematic hardware biases (e.g., thermal fluctuations in detector electronics or mechanical stretching). The normtags are mathematically designed to flatten these anomalies, artificially smoothing the decay curve to fit classical emittance models.

To reduce dependence on offline physics calibration choices, the extraction pipeline queries the brilcalc API via the public web cache while explicitly omitting the --normtag argument. This returns uncalibrated High-Level Trigger (HLT) luminosity suitable for a residual audit. The result is not assumed to be cleaner or more physical than calibrated data; it is treated as a complementary source that must pass independent controls.

## 3.3 The Sidereal Transformation Pipeline

The extraction of the uncalibrated CMS data yields a schema containing the delivered luminosity per LS and a strict UTC timestamp, synchronized by the LHC's sub-nanosecond White Rabbit timing network. To test the TEP hypothesis, this temporal axis must be transformed from a solar-centric coordinate system to a cosmic vector.

A fully deterministic Python pipeline is deployed, utilizing the astropy.time and astropy.coordinates libraries. The script selects data flagged as STABLE BEAMS to reduce active machine-tuning intervals. The rigid UTC timestamps are then mathematically mapped to Local Sidereal Time (LST) based on the exact geographical coordinates of the LHC ring in Geneva:

#### LHC Coordinates

- Latitude: 46.23° N

- Longitude: 6.05° E

This transformation re-indexes the instantaneous luminosity from the 24.00-hour solar day to the 23.93-hour sidereal cycle, aligning the data with the rotational orientation of the Earth relative to the galactic background.

## 3.4 Isolation of the Holonomy Residuals

For each independent continuous fill, the classical physics baseline must be subtracted to reveal the temporal scalar interaction. The pipeline calculates the relative elapsed time (t) from the start of the fill and fits a standard classical exponential decay envelope to the raw delivered luminosity:

\begin{equation}
L_{\text{classical}}(t) = L_0 e^{-\lambda_0 t}
\end{equation}

Where L_0 is the initial luminosity and λ_0 is the classical decay constant driven by proton burn-off and geometric emittance growth. The theoretical curve is then subtracted from the uncalibrated HLT data. The resulting luminosity residuals are treated as candidate observables, not as direct measurements of Synchronization Holonomy. These residuals across Run 2 fills are then folded into sidereal bins and tested against control periods.

## 4. Results

## 4.1 Claim Ladder

The pipeline now reports an explicit claim status. The current status is **null_with_secondary_lead**: No configured primary analysis currently satisfies the candidate gates. The CMS quality-filtered L0 test is nominally significant (p=0.027), but it is secondary and is not replicated by the primary lambda_fit test.

The ladder is: null (no robust primary evidence), candidate (a primary test survives initial gates), strong candidate (candidate plus wrong-period and replication controls), and proof candidate (independent source/era/outcome replication).

## 4.2 CMS/IP5 Time-Resolved and Fill-Level Tests

CMS is the public, time-resolved starting point, but the present CMS results do not support a proof claim. The row-level Lomb-Scargle scan is dominated by a 20.83 h period after winsorization; the sidereal power is only 0.086 of the dominant power, while the solar term is 0.335.

| Group | Outcome | n | delta R2 | p | Wrong-period controls stronger than sidereal |
| --- | --- | --- | --- | --- | --- |
| All CMS units | lambda_fit | 202 | 0.00234 | 0.802 | 18.0h, 20.0h, 21.0h, 22.0h, 23.0h, 25.0h, 28.0h, 30.0h |
| All CMS units | L0 | 202 | 0.00120 | 0.618 | 18.0h, 21.0h, 22.0h, 26.0h, 28.0h |
| Quality-filtered CMS | lambda_fit | 125 | 0.00225 | 0.584 | 18.0h, 20.0h, 23.0h, 25.0h, 28.0h, 30.0h |
| Quality-filtered CMS | L0 | 125 | 0.01347 | 0.027 | none |

## 4.3 Multi-IP Fill-Level Audit

The fill-level SuperTable sidecar extends the audit to ALICE/IP2, ATLAS/IP1, CMS/IP5, and LHCb/IP8. It is not time-resolved proof data, but it is useful for cross-IP triage and for identifying where time-resolved follow-up should be concentrated.

| Experiment/IP | n peak | Peak delta R2 | Peak p | Peak q | Delivered p | Claim audit |
| --- | --- | --- | --- | --- | --- | --- |
| ALICE/IP2 | 1357 | 0.00102 | 0.210 | 0.669 | 0.149 | candidate_not_proof |
| ATLAS/IP1 | 1398 | 1.72e-4 | 0.564 | 0.723 | 0.287 | candidate_not_proof |
| CMS/IP5 | 1399 | 2.67e-4 | 0.335 | 0.669 | 0.632 | candidate_not_proof |
| LHCb/IP8 | 1386 | 3.50e-4 | 0.516 | 0.723 | 0.860 | candidate_not_proof |

## 4.4 LHCb/IP8 Time-Resolved Pathway

LHCb/IP8 remains the decisive target because it is the interaction point most directly implicated by the earlier fill-level triage and by the available LHCb-specific luminosity products. However, after requiring positive stable-beams duration in the SuperTable audit, LHCb/IP8 no longer survives as a fill-level peak-luminosity candidate. The decisive next dataset is therefore time-resolved LHCb luminosity from Massi archives or CERN NXCALS exports, not another pass over zero-duration fill summaries.

The time-resolved LHCb pathway is now pre-registered in the pipeline: primary luminosity must survive the frozen sidereal test, a dense 10-50 h alias scan, leave-one-year-out stress tests, year-phase consistency, independent Run 2 and Run 3 replication with consistent sign, and delivered-luminosity replication.

| Source | Input exists | Rows | Fills | Primary outcome | Status | p | Failed gates |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Massi local archive | no | 0 | 0 | n/a | missing_input | n/a | n/a |
| NXCALS LHC_STATS export | no | 0 | 0 | n/a | missing_input | n/a | n/a |

## 4.5 Public LPC Massi-Derived Cross-Check

During source discovery, the pipeline found public CERN LPC Run 2 luminosity tables derived from Massi files. These tables include fill-level LHCb integrated and peak luminosity for 2015-2018. They are not time-resolved proof data, but they provide an independent public fill-level cross-check.

The public LPC cross-check does not reproduce a peak-luminosity candidate: LHCb peak luminosity gives p=0.104 and q=0.176. The LPC integrated-luminosity column is cumulative, so the pipeline converts it to per-fill increments before testing delivered luminosity; that delivered-increment test also fails (p=0.176). A 28.0h control exceeds the sidereal peak term.

| Source | LHCb rows | Peak p | Peak q | Delivered p | Wrong periods stronger | Status |
| --- | --- | --- | --- | --- | --- | --- |
| CERN LPC public Massi-derived tables | 780 | 0.104 | 0.176 | 0.176 | 28.0h | not_candidate |

## 4.6 Injection Recovery

Blind injection tests show that the pipeline can recover synthetic sidereal signals on the real schedules. For all CMS units, the null residual baseline has p=0.810; injected signals at 0.5 to 1.5 residual-sigma are recovered with p-values below 0.012. This validates sensitivity without converting the observed CMS data into a detection.

## 5. Discussion

## 5.1 Interpretation

The present analysis should be interpreted as a disciplined null-to-weak-candidate search, not a validation of TEP. CMS/IP5 is mostly null under the configured robust tests. The earlier LHCb/IP8 fill-level candidate does not survive the corrected stable-beams filter, and the public LPC Massi-derived table does not reproduce LHCb peak significance.

## 5.2 Alternative Explanations and Controls

The main risk is schedule aliasing: operational choices, fill duration, year structure, and machine settings can create apparent periods near a sidereal day. The pipeline therefore includes solar terms, wrong-period controls, year-preserving permutations, leave-one-year-out stress tests, source consistency checks, and delivered-luminosity replication. The current CMS wrong-period controls are strong evidence against over-interpreting CMS residuals as a sidereal discovery.

## 5.3 What Would Move the Claim Up the Ladder?

- Time-resolved LHCb/IP8 luminosity must show a sidereal term after elapsed-fill, solar/lunar, machine-state, and fill fixed-effect controls.

- The same phase and sign must replicate independently in Run 2 and Run 3.

- Peak and delivered luminosity should agree, or the discrepancy must be physically explained.

- A dense wrong-period scan from 10 to 50 h must not find nearby aliases stronger than the sidereal term.

- IP geometry predictions should be frozen before comparing phase offsets across IP1/IP2/IP5/IP8.

## 5.4 Theory Rigor

The Synchronization Holonomy is sourced by the non-exact disformal transport structure of the matter metric. The closed-loop integral of the conformal gradient $\Sigma_\mu = \nabla_\mu\ln A$ vanishes identically; the residual holonomy is generated by the disformal correction $\delta\tilde{\sigma}$ arising from $B(\phi) \neq 0$. This anchors the LHC closed-orbit observable to the correct sector of the TEP metric and removes the mathematical ambiguity previously noted in this section.

## 5.5 Limitations

The current public CMS analysis is limited by detector/source systematics and by the absence of verified beamline geometry projections. The multi-IP SuperTable and public LPC tables are fill-level rather than luminosity-section-level, so they can triage candidates but cannot by themselves prove time-resolved sidereal beam degradation. The corrected SuperTable filter and public LPC table both weaken the current LHCb peak-luminosity case. LHCb Massi and NXCALS inputs remain the primary missing proof-grade pathway.

## 6. Conclusion

## 6.1 Summary of Findings

The current TEP-LHC pipeline does not support a proof claim or a primary fill-level candidate. It produces a transparent **null_with_secondary_lead** status: No configured primary analysis currently satisfies the candidate gates. The CMS quality-filtered L0 test is nominally significant (p=0.027), but it is secondary and is not replicated by the primary lambda_fit test.

- CMS/IP5 time-resolved and fill-level tests are mostly null under robust controls.

- The earlier LHCb/IP8 SuperTable lead was an artifact of including zero-duration fill rows; after filtering stable-beams duration, LHCb peak luminosity is not significant.

- A public CERN LPC Massi-derived Run 2 cross-check does not reproduce the LHCb peak-luminosity candidate.

- The LHCb time-resolved proof test is now implemented with frozen alias, replication, phase, and delivered-luminosity gates.

- Blind injections demonstrate that the pipeline can recover imposed sidereal structure on real schedules.

- The next decisive step is time-resolved LHCb/IP8 analysis with frozen replication gates.

## 6.2 Significance

The value of this version is methodological. It defines a reproducible ladder for separating null results, candidates, strong candidates, and proof candidates. That ladder prevents the manuscript from converting a suggestive fill-level feature into a detection claim before independent time-resolved replication exists.

## 6.3 Future Work

- Ingest CERN Massi or NXCALS LHCb/IP8 time-resolved luminosity exports.

- Run the implemented time-resolved LHCb gate suite after Massi or NXCALS exports are available.

- Freeze IP geometry and predicted phase offsets before cross-IP comparison.

- Add calibrated-versus-uncalibrated source comparisons when matched LHCb calibration products are available.

## 6.4 Final Remarks

TEP-LHC is currently best presented as a falsifiable proof pathway rather than a detection. The public LPC cross-check makes the present empirical case more conservative, and the pipeline has narrowed the decisive question to time-resolved LHCb/IP8 luminosity with rigorous replication and control gates.

## 7. References

[1] CMS Collaboration. CMS Open Data Portal and luminosity information for public collision datasets. CERN Open Data Portal.

[2] CMS BRIL Project. *brilcalc* luminosity tooling and BRIL Work Suite documentation.

[3] CERN Beams Department. LHC fill summaries and SuperTable operational data products.

[4] Scargle, J. D. (1982). Studies in astronomical time series analysis. II. Statistical aspects of spectral analysis of unevenly spaced data. *The Astrophysical Journal*, 263, 835-853.

[5] Zechmeister, M., & Kürster, M. (2009). The generalised Lomb-Scargle periodogram. *Astronomy & Astrophysics*, 496, 577-584.

[6] Kostelecký, V. A., & Russell, N. (2011). Data Tables for Lorentz and CPT Violation. *Reviews of Modern Physics*, 83, 11.

## 8. Data Availability & Reproducibility

## 8.1 Data Access

The CMS luminosity data used in this analysis is publicly available through the CMS Open Data portal. The analysis uses uncalibrated High-Level Trigger (HLT) luminosity measurements obtained via the brilcalc API without normtag calibration.

#### Data Sources

- CMS Open Data: https://cms.cern.ch/

- brilcalc API: https://cms-service-dqm.web.cern.ch/brilcalc/

- LHC Run 2 fills: 2015-2018

## 8.2 Code Availability

The analysis pipeline is implemented in Python and will be made publicly available in the project repository. The code includes:

- Data extraction scripts for CMS brilcalc API

- Sidereal time transformation utilities

- Frequency analysis tools (Lomb-Scargle periodogram)

- Visualization and plotting routines

#### Example: Sidereal Transformation

```
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord
import astropy.units as u

# LHC coordinates
lhc_location = EarthLocation(lat=46.23*u.deg, lon=6.05*u.deg)

# Convert UTC to LST
def utc_to_lst(utc_timestamp):
t = Time(utc_timestamp, format='unix', location=lhc_location)
lst = t.sidereal_time('apparent')
return lst.hour
```

## 8.3 Reproducibility Checklist

#### To reproduce this analysis:

- Install Python 3.8+ with required dependencies (astropy, numpy, scipy, pandas)

- Download CMS luminosity data via brilcalc API

- Run the data extraction script: python scripts/steps/step_001_download_lhc_data.py

- Run the sidereal transformation: python scripts/steps/step_003_sidereal_analysis.py

- Run the frequency analysis: python scripts/steps/step_003_sidereal_analysis.py

- Results will be saved in the results/ directory

## 8.4 Computational Requirements

The analysis can be performed on a standard laptop or desktop computer. Typical computational requirements:

- RAM: 4 GB minimum

- Disk space: 10 GB for raw data and results

- Processing time: 1-2 hours for full Run 2 analysis

## 8.5 Version Control

The analysis code and manuscript are maintained in a Git repository. The repository includes:

- Complete analysis pipeline with version history

- Manuscript source files

- Documentation and README files

- Issue tracking for bug reports and feature requests

## 8.6 Open Science Statement

This analysis follows open science principles. All data used is publicly available, all code will be released under an open-source license, and the manuscript will be made available as a preprint prior to peer review.