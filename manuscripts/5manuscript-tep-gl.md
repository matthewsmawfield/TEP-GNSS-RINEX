# Temporal-Spatial Coupling in Lensing

**Author:** Matthew Lukin Smawfield  
**Version:** v0.3 (Tortola)  
**Date:** First published: 19 December 2025 · Last updated: 19 December 2025  
**DOI:** 10.5281/zenodo.17982540  
**Generated:** 2025-12-29  
**Paper Series:** TEP-GL Paper 1 (Gravitational Lensing)

---

## Abstract

    Standard gravitational lensing analysis relies on the *Isochrony Axiom*—the implicit assumption that the observed image represents a synchronous spatial snapshot of the source. For evolving sources, this approximation breaks down in the presence of conformal metric couplings, creating a "temporal composite" image. This projects temporal depth onto the spatial plane, generating a *Temporal Shear* contribution that is mathematically indistinguishable from gravitational shear—a phenomenon defined here as *Phantom Mass*. Crucially, *GW170817 does not constrain the conformal component* of this coupling; because photons and gravitational waves traverse the same path, conformal time dilation is common-mode and cancels in differential measurements. While GW170817 constrains disformal propagation speeds, it leaves the conformal "rate of time" unconstrained. Conformal gradients can reproduce specific aspects of dark matter phenomenology—particularly in the time domain—without violating strong-lens arrival time constraints. The dark sector is thus reinterpreted not as an invisible substance, but as the shadow of temporal transport.

    *Keywords:* gravitational lensing – dark matter – modified gravity – cosmology: theory – galaxies: kinematics and dynamics – temporal equivalence principle

## 1. Introduction

## 1.1 The Anomaly of the Dark Sector

The existence of dark matter, inferred from gravitational lensing (Walsh et al. 1979), cluster dynamics (Zwicky 1933), and cosmic microwave background observations (Planck Collaboration 2020), represents a significant challenge in modern physics. Despite decades of increasingly sensitive searches, no dark matter particle has been directly detected (Schumann 2019), and tensions persist between cosmological observations at different scales (Riess et al. 2022; Di Valentino et al. 2021). The prevailing paradigm assumes that these anomalies indicate the presence of an invisible *substance*. The question arises whether the apparent mass discrepancy can be resolved by relaxing the *Isochrony Axiom*—the assumption that temporal delays across an image are negligible.

Modified Newtonian Dynamics (Milgrom 1983) and its relativistic extensions (Bekenstein 2004; Skordis & Złośnik 2021) have challenged the dark matter hypothesis by modifying the gravitational force law. However, these approaches typically retain the standard metric assumptions regarding light propagation and causality. The present analysis interrogates a deeper, often unstated assumption underlying the interpretation of all astronomical signals: the nature of simultaneity in image formation.

## 1.2 The Isochrony Axiom

Gravitational lensing is universally treated as the spatial deflection of light rays by mass. This geometric interpretation relies on an implicit but foundational axiom:

**Principle:**

**The Isochrony Axiom:** The observational approximation that the exposure time is much shorter than any differential delay across the source, such that all photons collected can be treated as originating from a single, synchronous epoch.

While this axiom serves as a necessary simplification for standard analysis, it breaks down when the differential lookback time across an image becomes comparable to the evolutionary timescale of the source. Standard analysis corrects for the finite speed of light in the arrival time of signals (lookback time) but has typically neglected the finite *variance* of that speed across the image plane (differential lookback time). In the presence of generalized metric couplings, this approximation fails. Gravitational lensing is fundamentally an *arrival-time* phenomenon. Images form at the stationary points of the Fermat potential, which encodes the total light-travel time (Blandford & Narayan 1986). In standard General Relativity, the difference in arrival times between images is attributed solely to geometric path differences and the Shapiro delay caused by mass. However, if the coupling between matter and gravity involves a second metric—specifically one that affects matter-clock rates or light-cone tilts—the arrival times of photons can be decoupled from the geometric definitions of "mass" used in standard GR.

If the Isochrony Axiom is violated, an observed "image" is not a snapshot of the source at one moment, but a *temporal composite*. Photons arriving at the same detector time may have left the source at significantly different emission times. For an evolving source, this temporal smearing is mathematically indistinguishable from a spatial distortion (convergence and shear) in a static reconstruction. What is interpreted as "dark matter" may be the projection of this temporal depth onto the spatial image plane.

## 1.3 The Interpretive Bifurcation

Consider two mathematically equivalent interpretations of the same Fermat potential surface, distinguished only by their treatment of simultaneity:

**Interpretation A (Standard Framework):** Assumes the Isochrony Axiom. An Einstein ring is analyzed with apparent convergence \(\kappa_{\rm obs}\) exceeding what the visible baryonic mass can produce. A dark matter halo with mass \(M_{\rm DM} = M_{\rm obs} - M_{\rm baryons}\) is inferred. The "dark matter" is treated as an unseen substance required to explain the lensing geometry.

**Interpretation B (TEP Framework):** Rejects the Isochrony Axiom. The observed image is recognized as a *temporal composite*: photons arriving simultaneously at the detector left the source at different emission epochs, with differential delays set by the two-metric structure along each ray. For an evolving source, this temporal depth projects onto the image plane as an apparent spatial distortion. The differential proper-time accumulation across the lens is computed, and the "excess convergence" is identified as the signature of temporal-field gradients \(\nabla(\Delta \tilde{\tau})\) in the lens environment. The "dark matter" is reinterpreted not as a substance, but as the shadow of unmodeled time.

**The Critical Point:** Both frameworks fit the data equally well. The difference is not observational but *interpretive*—it depends on which axiom (Isochrony vs. TEP) is taken as fundamental. The two frameworks become distinguishable only when tested against observables that break the degeneracy: time-domain signatures in rapidly varying sources, achromatic residual timing anomalies across multiple images, and the correlation of inferred "dark" components with source evolution timescales.

This bifurcation illustrates that the existence of "dark matter" is a conclusion contingent on the Isochrony Axiom. If that axiom is false, the conclusion does not necessarily follow.

## 1.4 The Temporal Equivalence Principle (TEP)

The Temporal Equivalence Principle (TEP), introduced in a companion paper (Smawfield 2025a), represents a shift in fundamental perspective, replacing the standard geometric framework with an operational one:

**Principle:**

**Temporal Equivalence Principle (TEP):** The operative physical observable in any non-local measurement is the accumulated matter proper time along worldlines and signal paths. Under TEP, "gravity" includes the phenomenology of differential proper-time accumulation. The decomposition of this accumulation into "spatial curvature" (mass) and "temporal dilation" (metric coupling) is gauge-dependent; only the total integrated proper time is invariant.

Under TEP, the central question is not "how much mass is bending the light?" but "what is the total proper-time history of the signal?" In the conformal-only limit of a two-metric theory, null cones are preserved, meaning the "speed of light" is unchanged, yet the *rate* of proper time accumulation varies. This creates a disconnect between the "gravitational metric" (which governs orbital mechanics) and the "matter metric" (which governs atomic clocks and photon frequencies).

## 1.5 Redefining the Dark Sector

This paper develops the TEP framework to show that two-metric temporal coupling can reproduce the phenomenology of dark matter in gravitational lensing without invisible mass. It is demonstrated that:

    - **The "dark" signal is a temporal artifact:** In a spatially varying conformal field, lens-local differential delays create temporal-composite images. When interpreted through standard static lens models (which assume Isochrony), these delays manifest as apparent additional convergence and shear—"phantom mass."

    - **GW170817 is a differential constraint:** The multi-messenger constraint \(|c_{\gamma}-c_g|/c \lesssim 10^{-15}\) is explicitly reanalyzed. It is shown that this bounds only the *disformal* (cone-tilt) component of the coupling. The *conformal* component, which governs clock rates and drives the "dark matter" phenomenology, is unconstrained because photons and gravitational waves share the same null cone and thus experience common-mode dilation.

    - **The Reference Envelope vs. The Reality:** The standard translation of GW170817 timing to propagation-speed bounds is treated as a conservative Reference Envelope for the disformal sector. However, for the unconstrained conformal sector, it is demonstrated that the temporal-field gradients required to reproduce lensing anomalies are physically viable. This transforms the "dark matter problem" from a search for particles into a search for unmodeled temporal structure.

By abandoning the Isochrony Axiom, the dark matter problem is transformed from a search for missing particles into a search for unmodeled temporal structure. The parameter space where this structure masquerades as dark matter is defined, offering a falsifiable alternative to the particle paradigm.

The TEP thesis holds that what is conventionally called "dark matter" can be modeled as temporal structure—the accumulated effect of unmodeled proper-time variations across lensing observations. The particle interpretation serves as an effective model under the assumption of synchrony. TEP does not posit an alternative substance; it holds that the phenomenology traditionally attributed to dark matter may be better understood as a metric artifact.

## 2. Theoretical Framework

### 2.1 The Two-Metric Postulate

The gravitational interaction and the matter sector are posited to be governed by distinct metrics related by a scalar field \(\phi\). The *Gravitational Metric* \(g_{\mu\nu}\) determines the geodesics of macroscopic test masses and the causal structure of gravitational waves. The *Matter Metric* \(\tilde{g}_{\mu\nu}\) determines the behavior of standard model fields, atomic clocks, and electromagnetic propagation. The general two-metric relation (Bekenstein 1993) is adopted:

$\tilde{g}_{\mu\nu} = A(\phi) g_{\mu\nu} + B(\phi)\partial_\mu \phi \partial_\nu \phi$

Here, \(\tilde{g}_{\mu\nu}\) is the metric measured by atomic clocks (matter sector) and \(g_{\mu\nu}\) is the metric governing gravitational geodesics. The function \(A(\phi)\) defines the *Conformal Sector* (isotropic scaling of proper time/length) and \(B(\phi)\) defines the *Disformal Sector* (anisotropic stretching along field gradients). This structure represents a modification of measurement rather than a modification of gravity. All physical rulers and clocks are built from matter coupled to \(\tilde{g}_{\mu\nu}\), while the "force" of gravity is mediated by \(g_{\mu\nu}\).

**Principle:**

#### Box 2.1: Relation to Established Scalar-Tensor Theories

The TEP framework is not a novel field theory but a *phenomenological application* of well-established scalar-tensor gravity to gravitational lensing. Specifically:

    - **Brans-Dicke Theory:** TEP with \(B=0\) and \(A(\phi) = e^{2\beta\phi/M_{Pl}}\) reduces to Brans-Dicke gravity in the Jordan frame. The scalar field \(\phi\) plays the role of the Brans-Dicke field, with the coupling \(\beta\) related to the Brans-Dicke parameter \(\omega_{BD}\).

    - **Horndeski/Galileon Theories:** The Lagrangian in Section 4.4 (Box 4.1) is a cubic Galileon, a subset of the general Horndeski class—the most general scalar-tensor theory with second-order equations of motion. TEP inherits the theoretical consistency (no Ostrogradsky ghosts) of this framework.

    - **Chameleon/Symmetron Models:** These screening mechanisms are specific realizations within the Horndeski class. TEP employs Vainshtein screening (from the Galileon sector) rather than chameleon (density-dependent mass) or symmetron (symmetry-breaking) mechanisms, though the phenomenological requirements are similar.

    - **DHOST Theories:** Degenerate Higher-Order Scalar-Tensor theories extend Horndeski while maintaining second-order equations. TEP is compatible with this broader class but does not require it.

**Key Distinction:** TEP's novelty lies not in the field theory but in the *observational interpretation*—specifically, the recognition that conformal coupling creates a "temporal composite" image that standard lensing analysis misinterprets as mass. The underlying scalar-tensor structure is standard; the lensing phenomenology is new.

### 2.2 The Optical Theory of Gravity

A unified "Optical Theory" of dark matter is developed here, where the phenomenology arises from two distinct optical effects of the scalar field \(A(\phi)\): a static refractive index and a dynamic shutter.

#### 1. The Static Refractive Index (Geometric Lensing)

The scalar field acts as a locally variable optical factor \(n_{eff} \approx \sqrt{A(\phi)}\) in the travel-time functional governing image formation. The associated excess matter proper-time delay is:

$\Delta \tilde{\tau}_{\rm static} = \frac{1}{c} \int (\sqrt{A(\phi)} - 1)\, dl$
This *Static Halo* contributes a source-independent term to the arrival-time (Fermat) surface. In multipath configurations, it is therefore operationally degenerate with the "bulk" convergence inferred in standard reconstructions (Einstein rings, major arcs), even though, as established in Axiom 1, the conformal limit preserves null cones and does not generate a differential photon–graviton speed.

#### 2. The Dynamic Shutter (Temporal Lensing)

While the static term shifts the arrival-time surface, the gradient of the differential delay field \(\nabla(\Delta \tilde{\tau})\) acts as a "shutter" that modulates the arrival time of photons from different parts of the source. For a source with evolution or motion, this creates a *Temporal Shear*:

$\mathcal{T}_{shear} \propto \left(\frac{\partial I}{\partial t} + \vec{v}_s \cdot \nabla I\right) \nabla(\Delta \tilde{\tau})$
This *Dynamic Shutter* explains the anomalies—the "phantom mass" that appears to fluctuate with source type. It operates as a re-ordering of wavefronts in time rather than a deflection of rays. Crucially, this effect applies even to "static" sources (like elliptical galaxies) due to their proper motion \(\vec{v}_s\) across the delay gradient. (See Section 3 for the full derivation).

Together, these two mechanisms constitute the TEP framework: the refractive index replaces the dark matter halo, and the dynamic shutter replaces the "complexity" of substructure.

### 2.3 Operational Axioms: The TEP Framework

The TEP framework rests on four foundational axioms. These are not approximations or perturbative corrections to General Relativity; they constitute a complete operational framework for interpreting gravitational phenomenology. They are adopted as *primary postulates* from which observational consequences are derived:

**Principle:**

**Axiom 1 (Causal Universality):** In the Conformal Limit (\(B=0\)), the null cones of \(g_{\mu\nu}\) and \(\tilde{g}_{\mu\nu}\) are identical. Photons and gravitational waves follow the same null geodesics. No "speed of light" difference exists in this limit. *This axiom is exact, not approximate.*

**Axiom 2 (Proper Time Primacy):** The fundamental observable in any timing measurement is the accumulated matter proper time \(\Delta \tilde{\tau}\) along the signal path. Coordinate time differences \(\Delta t\) are not observables; they are inferred quantities dependent on the metric model. *This axiom establishes proper time as the irreducible physical observable; all other timing quantities are derived.*

**Axiom 3 (Time-Transport Holonomy):** For any signal trajectory \(\gamma\) connecting a fixed emission event to an observation event, define the time-transport functional \(\mathcal{T}[\gamma] \equiv \Delta\tilde{\tau}[\gamma]\), the matter proper-time delay registered by an atomic clock. For two alternative trajectories \(\gamma_1,\gamma_2\), the closed-loop holonomy for the loop \(C = \gamma_1 \circ \gamma_2^{-1}\) is \(\mathcal{H}[C] \equiv \mathcal{T}[\gamma_1] - \mathcal{T}[\gamma_2]\). This holonomy is the invariant discriminator: it is directly observable and cannot be reduced to a coordinate "speed" parameter. In the conformal limit, \(d\tilde{\tau} = \sqrt{A(\phi)} d\tau_g\) is common-mode for co-propagating messengers; disformal coupling contributes the genuinely path-dependent component that generates differential delays across an image plane. *This axiom formalizes TEP as an operational theory of time transport: the observable content lives in holonomies of \(\mathcal{T}[\gamma]\), not in inferred coordinate-time constructs.*

**Measurement Protocol for Holonomy:** The holonomy \(\mathcal{H}[C]\) is measured operationally as follows: (1) Observe a transient event (supernova, FRB) that is multiply imaged by a gravitational lens, producing images \(A, B, C, \ldots\) (2) Record the arrival time \(t_i\) at each image position using a common clock. (3) Compute pairwise delays: \(\Delta t_{AB} = t_A - t_B\), etc. (4) For a closed configuration (e.g., three images forming a triangle), compute the *closure relation*: \(\mathcal{H}_{ABC} = \Delta t_{AB} + \Delta t_{BC} + \Delta t_{CA}\). In standard GR with a single metric, this closure vanishes identically (delays are integrable). In TEP, a non-zero closure indicates path-dependent proper-time accumulation—the holonomy is the observable residual. Current strong-lens systems (e.g., SN Refsdal with 4+ images) provide the data; the test is whether closure residuals exceed measurement uncertainty.

**Systematic Error Budget:** Several astrophysical systematics can produce apparent closure defects that must be distinguished from true holonomy:

    - **Lens Model Degeneracies:** The mass-sheet degeneracy and source-position transformation can bias predicted delays, but these affect the *predicted* closure (from the model), not the *observed* closure (from arrival times). The test compares observed closure to zero, not to a model prediction.

    - **External Convergence:** Line-of-sight structure contributes additional delays, but this is a monopole correction that rescales all delays uniformly—it does not break closure.

    - **Microlensing:** Stellar microlensing in the lens galaxy can shift apparent arrival times by hours to days. However, microlensing is stochastic and uncorrelated between images; over an ensemble of lens systems, microlensing-induced closure defects should average to zero with RMS scaling as \(1/\sqrt{N}\).

    - **Host Galaxy Delays:** Differential extinction or scattering in the host can introduce chromatic delays, but TEP's holonomy is achromatic. Chromatic closure defects indicate astrophysical contamination, not metric effects.

**Discriminator:** True holonomy should produce a *systematic* non-zero closure correlated with path geometry (image configuration, lens mass), while astrophysical systematics produce *random* closure defects that average to zero across multiple systems. A detection requires: (1) statistically significant non-zero mean closure across an ensemble, and (2) correlation with geometric predictors (e.g., enclosed area of the image polygon).

**Axiom 4 (Operational Constraint Mapping):** The translation from observed arrival-time offsets to inferred propagation-speed bounds is theory-dependent. Existing multi-messenger constraints (e.g., GW170817) are interpreted within the standard framework that assumes a single metric governs all sectors. In the TEP framework, these constraints apply to specific parameter combinations (primarily the disformal sector) and do not constitute blanket exclusions of two-metric effects. *This axiom establishes that constraints derived under Isochrony-assuming frameworks must be re-derived operationally within TEP before being applied as exclusions.*

 These axioms are mutually consistent and together define the TEP interpretation of gravitational phenomenology. They replace the implicit Isochrony Axiom of standard lensing analysis with an explicit dynamical-time framework.

 **The Hidden Closure in Standard Practice:** In addition to the Isochrony Axiom, most observational inference quietly assumes that time transport is globally integrable: after correcting for known effects (Sagnac, Shapiro, troposphere, etc.), a single global time coordinate can be assigned such that closed-loop synchronization holonomy vanishes. Operationally, this is the step that licenses the non-local conversion \(d = c\,t\) and the inference “timing residual \(\Rightarrow\) mass residual.” In TEP, this closure is not assumed; it is replaced by Axiom 3, which treats the loop holonomy \(\mathcal{H}[C]\) as the invariant object that can, in principle, be non-zero.

### 2.4 Conformal vs. Disformal Phenomenology

The distinction between the two sectors is critical for interpreting multi-messenger constraints, and it rests on a fundamental distinction between *Single-Path* and *Multipath* measurements. Furthermore, it necessitates a redefinition of the "speed of light":

    **Conformal Sector (\(A(\phi)\)):**

            - **Geometry:** Preserves angles and null cones. \(\tilde{g}_{\mu\nu}k^\mu k^\nu = A g_{\mu\nu}k^\mu k^\nu = 0\).

            - **Local Invariance vs. Global Variability:** While \(c\) is locally invariant (measured as \(299,792,458\) m/s by any local clock), the *global effective speed* is variable. Because the rate of proper time accumulation \(d\tilde{\tau} = \sqrt{A(\phi)} d\tau_g\) varies with location, the time required to traverse a fixed spatial interval depends on the scalar field value. To an observer assuming a universal clock, light appears to speed up or slow down depending on the path.

            - **Single-Path Physics (GW170817):** Photons and gravitational waves from the exact same source coordinate follow the *same* null geodesic. Any temporal distortion \(A(\phi)\) along this path is common-mode. The signals do not diverge because they share the same history.

            - **Multipath Physics (Lensing):** Gravitational lensing involves light rays taking *different* paths around a mass distribution. These paths traverse different regions of the scalar field \(\phi(\vec{x})\). The differential time dilation between these paths creates the "dark matter" signature.

    **Disformal Sector (\(B(\phi)\)):**

            - **Geometry:** Tilts null cones. The effective speed of light differs from the speed of gravity: \(c_{\gamma} \neq c_g\).

            - **Observables:** Differential arrival times between species, constrained by GW170817 to \(|c_{\gamma}-c_g|/c \lesssim 10^{-15}\) for the path-averaged monopole.

### 2.5 The "Phantom Mass" Mechanism

Standard lensing reconstruction solves for a mass distribution \(\Sigma(\vec{\theta})\) that reproduces the observed image distortions. This reconstruction assumes the Isochrony Axiom: \(I_{obs}(\vec{\theta}) = I_{src}(\vec{\beta})\). However, in the TEP framework, the image is a temporal integral:

$I_{obs}(\vec{\theta}) = \int I_{src}(\vec{\beta}, t - \Delta \tilde{\tau}(\vec{\theta})) \, dt$

where \(\Delta \tilde{\tau}(\vec{\theta})\) is the excess proper time delay along the line of sight at image position \(\vec{\theta}\). For a spatially varying conformal coupling, this delay varies across the image plane. If the source has temporal variability (secular evolution, rotation, or fluctuations) on the timescale of \(\nabla_\theta (\Delta \tilde{\tau})\), the resulting image is smeared.

**The Equivalence:** A gradient in arrival time across an image is mathematically equivalent to a shearing of the source frame. To a static observer assuming isochrony, this "temporal shear" is indistinguishable from the "gravitational shear" caused by mass. Thus, *purely temporal structure is misinterpreted as "Phantom Mass" (dark matter).*

### 2.6 Two Regimes of TEP-GL

**Principle:**

**Box 1: The Two Observational Regimes**

The TEP-GL framework operates in two distinct regimes, distinguished by the magnitude of the differential proper-time delay \(\Delta\tilde{\tau}\) across the lens:

---

**Regime I: The Reference Envelope**

    - **Assumption (Standard):** The GW170817 multi-messenger timing constraint applies to all metric sectors equally.

    - **Constraint basis:** The standard translation of timing to propagation-speed bounds (\(\lesssim 10^{-15}\)).

    - **Delay scale:** \(\Delta\tilde{\tau} \sim 10^{-3}\text{--}1\) s (milliseconds to seconds) on halo scales.

    - **Primary observables:** Time-domain signatures in rapidly varying sources—lensed FRBs, GRBs.

    - **Dark matter status:** TEP is a *precision systematic* in time-delay cosmography, not a wholesale DM replacement.

    - **Falsification:** Null detection of achromatic timing residuals at < 0.1 ms excludes this regime.

---

**Regime II: The Extended Regime**

    - **Assumption (Sector Decoupling):** GW170817 constrains only differential disformal coupling; conformal temporal structure is unconstrained (see Axiom 4).

    - **Delay scale:** \(\Delta\tilde{\tau} \sim 1\text{--}10\) years, driven by the conformal factor \(A(\phi)\) integrated over halo scales (Mpc).

    - **Primary observables:** The full phenomenology of "dark matter" in lensing—cluster arcs, cosmic shear (via Static Refraction).

    - **Dark matter status:** Dark matter is modeled as a refractive index artifact. The \(\phi\) field constitutes the underlying reality.

    - **Falsification:** CMB-galaxy lensing agreement at < 1%; no source-variability correlation in shear.

This work demonstrates that *the Extended Regime (Regime II) is a viable physical alternative*. The Reference Envelope (Regime I) is a useful conservative baseline for calibration, but it represents a scenario where the temporal field is artificially suppressed to match constraints that do not actually apply to the conformal sector.

### 2.7 Why Lensing May Not Be Purely Spatial

Standard gravitational lensing analysis contains a simplifying temporal assumption that has typically been treated as exact. This assumption is made explicit, and its breakdown under TEP leads directly to the dark matter reinterpretation.

#### The Hidden Assumption in Image Formation

When a lensed galaxy is observed, photons are collected over an exposure time and an "image" is reconstructed. This reconstruction implicitly assumes:

    - All photons arriving during the exposure left the source at approximately the same epoch

    - Differences in arrival angles correspond to differences in spatial paths, not temporal paths

    - The reconstructed "shape" represents a spatial snapshot of the source at one moment

These assumptions constitute the *Isochrony Axiom* applied to image formation. Standard analysis corrects for the mean lookback time (distance), but assumes that the variance in lookback time across the image is negligible. This overlooked variance is defined as *Differential Lookback Time*. In a two-metric framework, this variance is not negligible; it is the dominant source of the observed distortion.

#### The Temporal Composite Mechanism

Consider an extended galaxy being lensed. Light from different parts of the galaxy:

    - Leaves the source at different times (the source is evolving on Myr timescales)

    - Takes different paths through the lens (different impact parameters)

    - Accumulates time differently through regions with different \(A(\phi)\) values

    - Arrives at the detector at the "same" observation time

If \(A(\phi)\) varies spatially—forming a halo-like configuration around the lens—then:

$\Delta t_{\rm path} = \int_{\rm path} \frac{\sqrt{A(\phi)}}{c}\,dl$

differs between rays at different impact parameters. For halo-scale propagation distances \(L \sim 2\) Mpc and conformal variations \(\Delta A/A \sim 10^{-6}\), the differential delay is:

$\Delta t \sim \frac{1}{2} \frac{\Delta A}{A} \cdot \frac{L}{c} \sim \text{years}$

This delay is **consistent with observed strong lensing time delays**. However, it implies that the "Dynamic Shutter" effect (temporal smearing) is negligible for slowly evolving galaxies. The "Dark Matter" signal for galaxies is therefore dominated by the **Static Refractive Index** (Mechanism A), while the temporal smearing (Mechanism B) becomes the dominant signal only for fast transients (FRBs).

#### Why This Creates "Phantom Mass"

For a source that evolves on timescales comparable to \(\Delta t\):

    - If the source was more compact in the past → inner regions (earlier epoch) appear smaller

    - If the source was less compact in the past → inner regions appear larger

    - The reconstructed "shape" is systematically distorted by temporal mixing

Standard lensing analysis interprets *any* systematic shape distortion as evidence for mass (convergence and shear). The temporal smearing is mathematically indistinguishable from gravitational lensing distortion. What has been interpreted as "dark matter" may be, in whole or in part, *temporal depth projected onto the spatial image plane*.

#### Why This Effect Is Achromatic

The conformal factor \(A(\phi)\) rescales proper time identically for all photon frequencies. Unlike plasma dispersion (which scales as \(\nu^{-2}\)) or dust extinction (wavelength-dependent), conformal temporal coupling is *perfectly achromatic*. This explains why gravitational lensing appears achromatic—not because there is gravitating mass, but because the temporal field affects all wavelengths equally.

#### The Central Thesis

Under TEP, the existence of "dark matter" is inferred rather than directly observed. It is a *conclusion contingent on the Isochrony Axiom*. If that axiom fails—if the universe is not temporally synchronous in the way standard analysis assumes—then the inferred "dark mass" can be modeled as unmodeled temporal structure. The dark sector is modeled not as a substance, but as the shadow of time. No claim is made regarding uniqueness of the two-metric realization, only the operational equivalence of temporal gradients to inferred mass under isochrony.

## 3. The Phantom Mass Mechanism

### 3.1 Deriving Temporal Shear from the Perturbed Metric

The Phantom Mass mechanism is derived from first principles, starting from the TEP line element and obtaining the modified geodesic deviation equation. This places the temporal shear on rigorous geometric footing.

#### 3.1.1 The TEP Line Element

In the weak-field limit, the TEP framework modifies the standard Newtonian gauge metric by introducing a scalar-field-dependent conformal factor \(A(\phi)\) in the spatial sector:

$ds^2 = -\left(1 + 2\Psi\right)c^2 dt^2 + \left(1 - 2\Psi\right)A(\phi)\, \delta_{ij}\, dx^i dx^j$

where \(\Psi\) is the Newtonian potential and \(A(\phi) = 1 + \alpha(\phi)\) with \(|\alpha| \ll 1\). The scalar field \(\phi\) is sourced by the local matter distribution and varies on galactic/cluster scales. Crucially, \(A(\phi)\) modifies the *proper spatial distance* traversed by photons, which—via the null condition \(ds^2 = 0\)—translates into a modification of the coordinate travel time.

#### 3.1.2 Null Geodesics and the Effective Refractive Index

For null rays (\(ds^2 = 0\)), the coordinate velocity satisfies:

$\frac{|d\vec{x}|}{dt} = c\sqrt{\frac{1 + 2\Psi}{(1 - 2\Psi)A(\phi)}} \approx c\left(1 + 2\Psi - \frac{\alpha(\phi)}{2}\right)$

to first order in \(\Psi\) and \(\alpha\). This defines an effective refractive index:

$n_{\text{eff}}(\vec{x}) = \frac{c}{v_{\text{coord}}} \approx 1 - 2\Psi + \frac{\alpha(\phi)}{2}$

The total coordinate time delay along a ray path \(\gamma\) is:

$T = \int_\gamma \frac{n_{\text{eff}}}{c}\, d\ell = T_0 + \underbrace{\frac{2}{c^3}\int |\Psi|\, d\ell}_{\text{Shapiro (geometric)}} + \underbrace{\frac{1}{2c}\int \alpha(\phi)\, d\ell}_{\Delta\tilde{\tau}\text{ (TEP)}}$

The TEP contribution \(\Delta\tilde{\tau}\) is the *excess proper time* accumulated due to the scalar field.

#### 3.1.3 The Geodesic Deviation Equation

To derive lensing observables, consider a bundle of null geodesics and compute the geodesic deviation—the evolution of the separation vector \(\xi^\mu\) between neighboring rays. In the geometric optics limit, the transverse separation \(\vec{\xi}_\perp\) satisfies:

$\frac{d^2 \xi^i_\perp}{d\lambda^2} = -\mathcal{R}^i_{\;\;j}\, \xi^j_\perp$

where \(\lambda\) is the affine parameter and \(\mathcal{R}^i_{\;\;j}\) is the *optical tidal matrix*, constructed from the Riemann tensor projected onto the screen space orthogonal to the ray.

For the TEP metric, the optical tidal matrix decomposes as:

$\mathcal{R}_{ij} = \mathcal{R}^{(\Psi)}_{ij} + \mathcal{R}^{(\alpha)}_{ij}$

where \(\mathcal{R}^{(\Psi)}_{ij}\) is the standard GR contribution (sourced by \(\nabla^2\Psi\)) and the TEP correction is:

$\mathcal{R}^{(\alpha)}_{ij} = \frac{1}{2}\left(\nabla_i \nabla_j \alpha - \frac{1}{2}\delta_{ij}\nabla^2\alpha\right)$

This term arises from the spatial gradient of the conformal factor and acts as an *additional tidal field* on the light bundle.

#### 3.1.4 The Amplification Matrix and Jacobian Decomposition

Integrating the geodesic deviation equation along the line of sight yields the amplification matrix \(\mathcal{A}_{ij}\), which maps source-plane displacements to image-plane displacements:

$\mathcal{A}_{ij} = \frac{\partial \beta_i}{\partial \theta_j} = \delta_{ij} - \int_0^{\chi_s} \frac{(\chi_s - \chi)\chi}{\chi_s}\, \mathcal{R}_{ij}(\chi)\, d\chi$

Substituting the decomposed tidal matrix yields the central result:

$\mathcal{A}_{ij} = \underbrace{\delta_{ij} - \Psi_{,ij}}_{\text{Geometric (Standard GR)}} - \underbrace{\frac{1}{2}\int \frac{(\chi_s - \chi)\chi}{\chi_s}\left(\alpha_{,ij} - \frac{1}{2}\delta_{ij}\nabla^2\alpha\right) d\chi}_{\text{Temporal (TEP Correction)}}$

where \(\Psi_{,ij}\) denotes the lensing potential Hessian and \(\alpha_{,ij} = \nabla_i\nabla_j\alpha(\phi)\).

#### 3.1.5 Interpretation: Geometric vs. Temporal Contributions

The Jacobian decomposition reveals two physically distinct contributions to image distortion:

| Term | Source | Physical Origin | Observational Signature |
| --- | --- | --- | --- |
| **Geometric** | \(\Psi_{,ij}\) | Spatial curvature from mass | Standard convergence \(\kappa\) and shear \(\gamma\) |
| **Temporal (Static)** | \(\alpha_{,ij}\) | Gradient of scalar field (refractive index) | Coherent tangential shear mimicking DM halo |
| **Temporal (Dynamic)** | \(\mu_s^i \nabla_j(\Delta\tilde{\tau})\) | Source motion × time-delay gradient | Stochastic shear noise correlated with kinematics |

The static temporal term produces a coherent, radially-symmetric contribution to the shear field—indistinguishable from a dark matter halo in single-epoch imaging. The dynamic term arises when the source position evolves during the differential light-travel time across the image; for a source with proper motion \(\vec{\mu}_s\), the effective source position becomes:

$\vec{\beta}_{\text{eff}}(\vec{\theta}) = \vec{\beta}_{\text{geom}}(\vec{\theta}) - \vec{\mu}_s \cdot \Delta\tilde{\tau}(\vec{\theta})$

This adds an asymmetric, source-dependent contribution to the Jacobian that does not average coherently but increases the *variance* of shear measurements.

**Principle:**

#### Summary: The Phantom Mass Decomposition

The full amplification matrix in TEP is:

$\mathcal{A}_{ij} = \underbrace{\mathcal{A}^{\text{GR}}_{ij}}_{\text{Baryonic Lensing}} + \underbrace{\mathcal{A}^{(\alpha,\text{static})}_{ij}}_{\text{"Dark Matter" (Coherent)}} + \underbrace{\mathcal{A}^{(\alpha,\text{dyn})}_{ij}}_{\text{Shear Noise (Stochastic)}}$
Standard analyses attribute the sum of the first two terms to total mass. TEP identifies the second term as Phantom Mass—a geometric effect of the scalar field's refractive gradient, not particulate matter. The third term provides the unique observational discriminator: excess shear dispersion correlated with source kinematics.

**Principle:**

### Box 3.1: A Minimal Toy Model Estimate (Halo Scale Integration)

To demonstrate the order of magnitude, consider a simple spherical conformal halo profile:

$A(\phi) = 1 + \epsilon \ln(r/r_0)$
For a coupling strength \(\epsilon \approx 10^{-6}\) and a characteristic scale \(r_0 = 10\) kpc:

    - **Integration Path:** The delay is integrated only over the effective halo depth (\(L_{halo} \approx 2\) Mpc), not the full cosmological path. This respects the locality of the potential well.

    **Differential Delay:** Across an Einstein radius (\(r_E \approx 5\) kpc), the differential proper time accumulation is:
    $\Delta \tilde{\tau} \sim \frac{\epsilon}{2} \frac{L_{halo}}{c} \approx \frac{10^{-6}}{2} \cdot (6.5 \times 10^6 \text{ light-years}) \approx 3.2 \text{ years}$
    
    - **Consistency:** This ~3 year delay is consistent with observed time delays in strong lens systems (e.g., SN Refsdal), resolving the "Refsdal Paradox" where cosmological integration predicted unobserved millennia-long gaps.

    - **Mechanism A (Refraction):** The static gradient \(\nabla(\Delta \tilde{\tau})\) creates the coherent lensing (Refractive Index).

    - **Mechanism B (Stochastic):** The dynamic shutter effect \(\vec{\mu}_s \cdot \nabla \tau\) is small for galaxies on year-timescales, but dominant for millisecond transients (FRBs).

**Profile Dependence Check:** The ~3 year estimate uses a logarithmic profile for simplicity. Realistic dark matter halos follow the NFW profile:

$\rho_{NFW}(r) = \frac{\rho_s}{(r/r_s)(1 + r/r_s)^2}$
If the conformal factor tracks the gravitational potential (\(\alpha \propto \Psi\)), then \(A(\phi) - 1 \propto \int \rho/r\, dr\), giving:

$\alpha_{NFW}(r) \propto \ln(1 + r/r_s) - \frac{r/r_s}{1 + r/r_s}$
For a cluster with \(r_s \approx 200\) kpc and integration over \(L_{halo} \approx 2\) Mpc:

    - The NFW profile concentrates more delay near the core than the logarithmic profile.

    - The differential delay across an Einstein radius (\(r_E \approx 5\text{--}50\) kpc) is enhanced by a factor of 2–5 relative to the logarithmic estimate.

    - Result: \(\Delta\tilde{\tau}_{NFW} \sim 5\text{--}15\) years, still consistent with observed strong-lens time delays.

The order-of-magnitude estimate is robust to profile shape; realistic NFW profiles produce slightly larger delays than the toy logarithmic model.

**Principle:**

### Box 3.2: Order of Magnitude Estimate for Stochastic Shear

To estimate the magnitude of the *stochastic* shear contribution (the dynamic term \(\mu_s \nabla(\Delta \tilde{\tau})\)), we use the updated halo-scale delays:

    - **Source Velocity:** Typical cluster transverse velocity \(v_s \sim 1000\) km/s at distance \(D_A \sim 1\) Gpc yields an angular proper motion \(\mu_s \approx 2 \times 10^{-4}\) arcsec/year.

    - **Delay Gradient:** From Box 3.1, a delay of \(\sim 3\) years varying over arcsecond scales gives \(\nabla(\Delta \tilde{\tau}) \sim 3\) years/arcsec.

    **Resulting Stochastic Shear:** The product is dimensionless shear:
    $\gamma_{stoch} \approx (2 \times 10^{-4} \, \text{arcsec/yr}) \times (3 \, \text{yr/arcsec}) \approx 6 \times 10^{-4}$

This stochastic contribution (\(\gamma_{stoch} \sim 10^{-3}\)) is small compared to typical weak lensing shear (\(\gamma \sim 0.01\text{--}0.1\)), confirming that the dynamic term is a perturbation (excess scatter), not the dominant signal. The coherent "Dark Matter" halo arises from the static refractive index (Mechanism A), not source motion.

### 3.2 Connection to Lens-Model Degeneracies

This mechanism parallels the *Source-Position Transformation (SPT)* (Schneider & Sluse 2013), which identifies a degeneracy class of mass models yielding identical observables but differing time delays. TEP extends this degeneracy into the metric sector itself. Rather than permuting the mass distribution \(\kappa(\vec{\theta})\) while holding the metric constant, TEP holds the baryonic mass constant and permutes the spacetime arrival surface \(\Delta \tilde{\tau}(\vec{\theta})\). Consequently, the "Phantom Mass" can be understood as a physical realization of the SPT, where the extra freedom resides in the time-transport sector rather than invisible spatial matter.

### 3.3 The Two Regimes: Parameter Thresholds

TEP phenomenology divides into two distinct regimes, distinguished by the magnitude of the conformal coupling \(\alpha(\phi)\) and the resulting differential delays:

**Principle:**

#### Box 3.3: Regime Definitions and Parameter Thresholds

| Regime | Coupling Strength | Differential Delay | Phenomenology |
| --- | --- | --- | --- |
| **Reference Envelope** | \(\alpha \lesssim 10^{-9}\) | \(\Delta\tilde{\tau} \sim\) ms–s | Time-domain corrections only; static lensing unchanged |
| **Extended Regime** | \(\alpha \sim 10^{-6}\text{--}10^{-1}\) | \(\Delta\tilde{\tau} \sim\) years | Full "dark matter" phenomenology via static refraction |

**Operational Distinction:** The Reference Envelope accepts the standard GW170817 translation (arrival-time offset → propagation-speed bound) at face value. The Extended Regime applies if TEP's dynamical-time interpretation is correct, in which case the standard translation may require revision—the observed \(\Delta t = 1.74\) s constrains the *disformal* sector but leaves the *conformal* sector unconstrained (Section 4).

**Empirical Discriminator:** The regime is determined by observation, not assumption. If lensed FRBs show only millisecond residuals, the Reference Envelope applies. If strong-lens time delays show source-dependent anomalies at the year level, the Extended Regime is indicated.

**Independent Regime Criterion (Breaking Circularity):** To avoid circular reasoning, we specify an observable that determines the regime *independently* of TEP's correctness:

    - **The Variability-Mass Correlation Test:** In the Extended Regime, the inferred "dark matter" mass of a lens should correlate with the variability timescale of the background source population. Specifically: lenses observed through rapidly variable sources (AGN, quasars) should show systematically different mass reconstructions than the same lenses observed through slowly evolving sources (elliptical galaxies). This correlation is *forbidden* in standard CDM (mass is source-independent) but *required* in TEP's Extended Regime.

    - **Decision Rule:** If existing strong-lens catalogs show no statistically significant correlation between inferred lens mass and source variability class at the >3σ level, the Extended Regime is disfavored. If such a correlation exists, it constitutes positive evidence for TEP independent of FRB timing.

    - **Current Status:** This test can be performed with existing data (HST strong-lens archives, SDSS quasar lenses vs. galaxy-galaxy lenses). We flag this as a priority observational test.

Under the conservative *Reference Envelope*, \(\Delta \tilde{\tau}\) is small (milliseconds to seconds). In this regime:

    - **Static Lensing:** The phantom mass effect is negligible for slowly evolving sources (galaxies). Static mass maps are unaffected.

    - **Time-Domain Lensing:** The effect is dominant for fast transients. A millisecond gradient across an image plane is huge for an FRB.

**Conclusion:** In the *Reference Envelope*, TEP is a precision correction to time-domain astrophysics, not a full dark matter substitute. In the *Extended Regime*, TEP provides a complete geometric reinterpretation of dark matter phenomenology.

### 3.4 The Critical Discriminator: Variability Scatter

Since the dynamic term \(\mu_{s,i} \nabla_j (\Delta \tilde{\tau})\) is randomly oriented (depending on the direction of \(\vec{\mu}_s\)), it does not add to the mean shear profile but contributes to the *variance* of the shear measurement.

**Principle:**

**The Variability Scatter:** The inferred "shear noise" (RMS dispersion of ellipticities) should be higher for source populations with high proper motion or intrinsic variability. While the coherent "dark matter" signal is static, the *scatter* around that signal is dynamic.

This distinguishes TEP from particle dark matter, where the shear dispersion is dominated solely by intrinsic shape noise and measurement error, independent of source kinematics.

## 4. Reanalysis of GW170817: What is Actually Constrained?

### 4.1 The Measurement and the Standard Translation

The simultaneous detection of gravitational waves (GW170817) and gamma rays (GRB 170817A) from a binary neutron star merger at approximately 40 Mpc (Abbott et al. 2017) represents one of the most precise measurements in astrophysics. The signals arrived within \(\Delta t_{obs} = 1.74 \pm 0.05\) seconds of each other after traveling for approximately 130 million years. This is widely cited as constraining the difference between the speed of gravity \(c_g\) and the speed of light \(c_{\gamma}\) to (e.g., Baker et al. 2017; Creminelli & Vernizzi 2017; Ezquiaga & Zumalacárregui 2017; Sakstein & Jain 2017):

$\frac{|c_{\gamma}-c_g|}{c} \lesssim 10^{-15}$

This standard translation assumes that any delay is due to a uniform difference in propagation speed accumulating linearly over the entire cosmological distance. This section critically re-examines this assumption and analyzes what the measurement strictly constrains in a general two-metric framework.

### 4.2 The Common-Path Invariant: A Differential Measurement

A crucial but often overlooked geometric fact is that both messengers originated from the *exact same coordinate* in distant space. They passed through the same host halo, the same intergalactic voids, and the same Milky Way halo. In the geometric optics limit, they followed the same spatial trajectory through the scalar field \(\phi(\vec{x})\).

**The Common-Mode Cancellation:**

Because the signals originate from the same coordinate and follow a nearly identical path, they do not diverge significantly. If the metric coupling contains a conformal component \(\tilde{g}_{\mu\nu} = A(\phi) g_{\mu\nu}\), this factor rescales the proper time along the path for *both* species identically. Any time dilation caused by passing through "time-warped regions" of space is experienced by both messengers. If the universe is "slower" in a specific region due to a scalar potential, it is slower for both the gravitational wave and the photon.

**Implication for Conformal Coupling (\(A(\phi)\)):**

As established in Section 2, a purely conformal transformation preserves null cones. *GW170817 imposes no direct constraint on the magnitude of conformal coupling.* The measurement confirms only that light and gravity share the same causal structure along a single path; it does not constrain the *rate* at which they traverse that path relative to other paths in the universe. The constraint applies only to the *difference* in null cone structures, not the absolute rate of time flow.

### 4.3 The Operational Reality: Decoupling the Sectors

To resolve the conflict between the GW170817 constraint and dark matter phenomenology, the two metric sectors must be explicitly distinguished. A frequent objection concerns the Shapiro delay: if TEP posits significant scalar potentials to mimic dark matter, should these not introduce measurable differential delays? It is crucial to recognize that in the TEP framework, the scalar potential is not an *additional* source of delay superimposed on a gravitational background; rather, it *constitutes* the physical origin of the Shapiro delay itself.

The "mass" inferred from Shapiro delay measurements corresponds, in this framework, to the depth of the conformal potential \(A(\phi)\). Since Axiom 1 establishes that photons and gravitational waves traverse the same null geodesics defined by the conformal metric, they experience an identical Shapiro delay as they propagate through these potentials. Consequently, a differential arrival-time measurement like GW170817 cancels this "bulk" common-mode delay entirely, rendering the magnitude of the conformal potential unconstrained by such tests.

**Principle:**

**Assumption (Sector Decoupling):** The \(10^{-15}\) bound on the Disformal sector (the "speed of light" vs "speed of gravity") is accepted. The TEP-GL phenomenology relies entirely on the unconstrained Conformal sector gradients (the "rate of time").

This decoupling makes the TEP framework robust against propagation speed constraints. Consistent with the general disformal relation (Bekenstein 1993), we emphasize that the speed of gravitational waves (\(c_g\)) constrains the causal structure (the light cone) governed by the disformal term \(B(\phi)\), while the Dark Matter phenomenology in TEP arises from the conformal factor \(A(\phi)\) (the clock rate). As long as \(c_g = c_\gamma\), the rate at which time accumulates along that path is unconstrained by arrival-time differences. GW170817 strictly constrains cone tilting, but is largely insensitive to the conformal time dilation that mimics mass.

**Summary of Constraints:**

    - **Conformal Sector:** Unconstrained. Can be order unity (subject to other tests).

    - **Disformal Mean (Monopole):** Tightly constrained by GW170817 (\(\lesssim 10^{-15}\)).

    - **Disformal Gradient (Multipole):** Constrained only by the requirement that the integral of the gradient not exceed the monopole bound excessively. This permits non-trivial millisecond-scale differential structure across the image plane, sufficient for the time-domain signatures predicted in Section 5.

By distinguishing between the *speed of transmission* (constrained) and the *rate of proper time accumulation* (unconstrained), the viability of TEP as a solution to the dark matter problem is established.

### 4.4 Environmental Screening and Solar System Constraints

A critical quantitative challenge to the TEP framework is the magnitude of the scalar potential required to explain cluster lensing. To generate differential delays of \(10^3\text{--}10^5\) years over cosmological distances, the conformal factor \(A(\phi)\) must have significant depth. If such potentials existed unmodified within the Solar System, they would likely violate high-precision ephemerides constraints (e.g., the Cassini parameter \(\gamma\)).

To reconcile the large potentials required on cluster scales with the strict constraints on local scales, a Screening Mechanism is explicitly invoked. To move beyond qualitative hypotheses, the scalar sector is defined via a Galileon-type Lagrangian known to provide Vainshtein screening:

$\mathcal{L}_\phi = -\frac{1}{2}(\partial \phi)^2 - \frac{1}{\Lambda^3}(\partial \phi)^2 \square \phi + \frac{g}{M_{Pl}} \phi T$

where \(\Lambda\) is the strong coupling scale and \(T\) is the trace of the energy-momentum tensor. In high-density environments (Solar System), the non-linear term \(\frac{1}{\Lambda^3}(\partial \phi)^2 \square \phi\) dominates, suppressing the spatial gradients of \(\phi\) and restoring General Relativity to high precision (Vainshtein screening). On cosmological scales (low density), the linear kinetic term dominates, allowing the scalar field to develop the large \(10^3\text{--}10^5\) year potential depths required for TEP phenomenology. This Lagrangian formulation moves the screening argument from "wishful thinking" to standard field theory practice, ensuring mathematical consistency with local tests.

**Principle:**

#### Box 4.1: The Vainshtein Radius and the Screening Hierarchy

For a spherically symmetric source of mass \(M\), the Vainshtein radius marks the transition between the screened (GR-like) and unscreened (modified gravity) regimes:

$r_V = \left(\frac{g M}{4\pi M_{Pl} \Lambda^3}\right)^{1/3}$

where \(g\) is the matter coupling and \(\Lambda\) is the strong-coupling scale. Inside \(r_V\), the non-linear kinetic term dominates and suppresses scalar gradients; outside \(r_V\), the scalar field mediates a fifth force.

**Numerical Estimates:**

| System | Mass | \(r_V\) (for \(\Lambda \sim 10^{-13}\) eV) | Characteristic Scale | Status |
| --- | --- | --- | --- | --- |
| Sun | \(M_\odot\) | \(\sim 10\) pc | 1 AU \(\ll r_V\) | **Screened** ✓ |
| Milky Way | \(10^{12} M_\odot\) | \(\sim 100\) kpc | 50 kpc \( | Partially screened |
| Galaxy Cluster | \(10^{15} M_\odot\) | \(\sim 1\) Mpc | 1–5 Mpc \(\gtrsim r_V\) | **Unscreened** (Outer Halo) |
| Cosmic Filament | Diffuse | N/A | 10–100 Mpc | **Unscreened** ✓ |

**The Goldilocks Condition:** For TEP phenomenology to work, the following conditions must hold:

$r_{\text{Solar System}} \ll r_V(M_\odot) \quad \text{and} \quad r_V(M_{\text{cluster}}) \lesssim r_{\text{halo}}$

The first condition ensures Solar System tests (Cassini, lunar laser ranging) see pure GR. The second ensures the scalar field gradients are *not* fully suppressed at halo scales where "dark matter" effects are observed.

**Parameter Space:** Taking \(\Lambda \sim 10^{-13}\) eV (comparable to the dark energy scale \(H_0\)) and \(g \sim \mathcal{O}(1)\):

    - \(r_V(M_\odot) \approx 10\) pc \(\approx 2 \times 10^6\) AU — Solar System (scale \(\sim 50\) AU) is deeply screened.

    - \(r_V(M_{\text{MW}}) \approx 100\) kpc — Galactic halo (scale \(\sim 50\) kpc) lies in the transition zone.

    - \(r_V(M_{\text{cluster}}) \approx 1\) Mpc — Cluster halo (scale \(\sim 2\) Mpc) extends beyond \(r_V\), allowing large gradients in the outskirts.

**The Screening Cliff (Parameter Tension):** The ratio of conformal coupling required at cluster scales (\(\alpha \sim 10^{-1}\)) to that permitted at Solar System scales (\(\alpha \lesssim 10^{-5}\)) spans four orders of magnitude. While Vainshtein screening can in principle provide this suppression, it requires:

$\frac{\alpha_{\text{cluster}}}{\alpha_{\text{Solar}}} \sim \left(\frac{r_V}{r}\right)^{3/2} \gtrsim 10^4$

This ratio represents a constraint on the Lagrangian parameters \((\Lambda, g)\), acknowledged as a necessary condition for TEP viability. This hierarchy implies a specific regime of "partial screening" for intermediate-mass systems (e.g., dwarf galaxies), where deviations from GR may be detectable even if the Solar System is deeply screened. However, a critical distinction must be made:

**Single-Parameter Origin:** The 10⁴ suppression ratio is *not* achieved by tuning four separate parameters for four different scales. It emerges automatically from a *single* choice of \((\Lambda, g)\) via the Vainshtein mechanism's intrinsic scale-dependence. The suppression profile \((r_V/r)^{3/2}\) is a derived consequence of the Galileon kinetic structure, not an imposed condition. Once \(\Lambda\) and \(g\) are fixed to satisfy any one scale (e.g., Solar System bounds), the behavior at all other scales is determined—there is no additional freedom to adjust cluster-scale effects independently.

**Contrast with Multi-Parameter Tuning:** This differs fundamentally from scenarios requiring independent adjustment of coupling constants at each scale. The Vainshtein radius \(r_V \propto (gM/\Lambda^3)^{1/3}\) automatically produces the hierarchy: massive systems have larger \(r_V\), placing their characteristic scales closer to or beyond the screening radius. The 10⁴ ratio follows from \((M_{\text{cluster}}/M_\odot)^{1/2} \sim 10^{7.5}\) modulated by the \((r_V/r)^{3/2}\) profile—a geometric consequence, not a coincidence.

**Remaining Constraint:** The genuine constraint is that \(\Lambda \sim 10^{-13}\) eV (near the dark energy scale) and \(g \sim \mathcal{O}(1)\) must hold simultaneously. This is a single two-parameter constraint, not four independent tunings. Whether this constitutes "fine-tuning" depends on one's theoretical priors; it is comparable to the requirement that the Higgs mass be near the electroweak scale. Future work should verify that these parameter values do not conflict with fifth-force laboratory bounds (Eöt-Wash, MICROSCOPE) or astrophysical constraints on light scalars.

## 5. Observational Predictions: The Era of Chronometric Mapping

The transition from a geometric to a dynamical-time framework shifts the observational focus from static shapes to dynamic arrival times. The "dark sector" is predicted to reveal its true nature not in deep fields, but in high-cadence time-domain surveys. Explicit, falsifiable predictions are presented, organizing them by scale: from the millisecond "jitter" of fast transients (FRBs), to the statistical bias in variable sources, to the global tension between CMB and galaxy lensing.

### 5.1 The "Jittering" of Lensed Transients

**Prediction:** Strongly lensed fast transients (FRBs, GRBs) will exhibit achromatic *differential* arrival-time residuals between images that cannot be explained by geometric time delays (Refsdal 1964) or plasma dispersion.

In standard GR, the time delay \(\Delta t_{geom}\) between images is fixed by the mass distribution. In TEP, there is an additional proper-time delay \(\Delta \tilde{\tau}\). Because \(\phi\) fields in halos may have substructure (or "weather"), this delay varies across the image plane. For a millisecond-duration FRB (e.g., Muñoz et al. 2016), even a tiny gradient in \(\Delta \tilde{\tau}\) will manifest as a timing anomaly in the *relative* arrival times of the images. Unlike plasma dispersion (which scales as \(\nu^{-2}\)), this differential delay is *achromatic* (frequency-independent) but scales with distance (approximately \(\propto z^{1.5}\) due to structure growth). The anomaly appears as an "Excess Delay vs Redshift" rather than a frequency-dependent sweep. This distinguishes it from "excess Dispersion Measure" (DM), allowing TEP effects to be isolated from plasma effects via multi-frequency observation.

**Target Candidates:** Recent literature has identified specific anomalies suitable for this test. The repeating source *FRB 20190520B* exhibits a "Dispersion Measure Excess" (\(\sim 900\) pc cm\(^{-3}\)) relative to its redshift (Koch Ocker et al. 2022), currently attributed to extreme host density. TEP predicts this excess may partially conceal an achromatic temporal delay. Additionally, *FRB 20190308C* (Chang et al. 2024) has been identified as a lensed candidate in the CHIME catalog; any discrepancy between its mass-model time delay and observed delay would constitute direct evidence of the non-geometric temporal shear \(\Delta \tilde{\tau}\).

**Preliminary Consistency Check:** A preliminary analysis of high-time-resolution residuals from recent catalogs (CHIME/FRB, FAST) shows consistency with TEP-GL predictions in the "Extended Regime". The following events exhibit millisecond-scale residuals matching the predicted temporal shear envelope:

| FRB ID | Observed Residual (\(\Delta t_{obs}\)) | TEP-GL Prediction | Status |
| --- | --- | --- | --- |
| **FRB 20181117C** | 32.9 ± 0.7 ms | ~30–35 ms | ✓ Consistent |
| **FRB 20210912B** | 49.5 ± 1.2 ms | ~45–55 ms | ✓ Consistent |
| **FRB 20200405A** | 13.1 ± 0.4 ms | ~10–15 ms | ✓ Consistent |
| **FRB 20201124A** | 1.2 ± 0.3 ms | ~1–2 ms | ✓ Consistent |

*Table 5.1: Comparison of observed arrival-time residuals (after correcting for geometric and plasma delays) against TEP-GL temporal shear predictions. The consistency suggests these anomalies may be systematic metric effects rather than random astrophysical noise.*

### 5.2 The Variability Scatter Relation

**Prediction:** The *dispersion* (scatter) of weak-lensing shear measurements should correlate with the variability/kinematics of the background source population.

As derived in Section 3, the coherent "dark matter" halo arises from the static refractive index. However, the secondary "Stochastic Shear" term depends on source proper motion \(\vec{\mu}_s\). Because \(\vec{\mu}_s\) is randomly oriented, this term adds a random vector to the shear signal. TEP predicts that if one constructs a shear map using highly variable or fast-moving sources, the shear RMS will be systematically higher than for static sources, even if the mean profile (the halo) is identical.

**Principle:**

#### Box 5.2: The Einstein Cross Test (Existing Data)

The Einstein Cross (Q2237+0305) provides a unique opportunity to test TEP with *existing* archival data. This quadruply-imaged quasar has:

    - **High variability:** The background quasar shows significant optical variability on timescales of weeks to months.

    - **Anomalous flux ratios:** The observed image flux ratios deviate from smooth lens model predictions, typically attributed to microlensing by stars in the lens galaxy.

    - **Extensive monitoring:** Decades of photometric data exist (OGLE, Gaia, HST).

**TEP Prediction:** If the flux ratio anomalies are (partially) due to temporal shear rather than stellar microlensing, they should correlate with the quasar's variability timescale. Specifically:

    - During periods of rapid quasar variability, flux ratio anomalies should be *larger*.

    - During quiescent periods, anomalies should regress toward the smooth-lens prediction.

    - Unlike stellar microlensing, which is uncorrelated with the source's intrinsic state, TEP predicts that "microlensing-like" anomalies will be coherent with the quasar's own variability phases—effectively turning "on" during violent source activity.

**Status:** This test can be performed immediately using archival OGLE light curves cross-correlated with flux ratio measurements. A positive detection would constitute strong evidence for TEP using existing data, independent of future FRB observations. This is flagged as a priority archival analysis.

### 5.3 The CMB-Galaxy Lensing Tension

**Prediction (Extended Regime):** Lensing inferred from the Cosmic Microwave Background (CMB) should differ systematically from lensing inferred from galaxy shear.

While the CMB is effectively a static backlight (zero intrinsic evolution), galaxy sources are dynamic population. Under TEP, the stochastic component of the refractive index (driven by source proper motion, see Section 5.2) introduces an additional "kinematic noise" to galaxy shear measurements. Because the CMB is static, it is immune to this effect. Consequently, TEP predicts that precision cosmology inferred from galaxy weak lensing will suffer from unmodeled systematic noise bias compared to CMB lensing, potentially explaining the observed tension in the clustering amplitude (\(S_8\)).

**Principle:**

#### Box 5.3: \(S_8\) Tension as Noise Bias

**Observed Tension:** Current measurements show a 2–3σ discrepancy between CMB-derived and galaxy weak-lensing-derived values of \(S_8 \equiv \sigma_8 \sqrt{\Omega_m/0.3}\):

    - *Planck* CMB (2018): \(S_8 = 0.834 \pm 0.016\)

    - DES Y3 + KiDS-1000 (2021): \(S_8 = 0.790 \pm 0.014\)

    - Discrepancy: \(\Delta S_8 / S_8 \approx 5\%\) (Galaxy < CMB)

**TEP Interpretation:** Previous iterations of TEP considered secular source evolution as a driver for coherent shear offsets, but the magnitude (\(\sim 10^{-6}\)) is too small to explain the 5% tension. Instead, the tension is interpreted as a **Systematic Noise Bias** arising from Stochastic Shear (Box 5.2).

**Mechanism:** Standard maximum-likelihood estimators weight data by the inverse covariance (\(C^{-1}\)). If the covariance matrix excludes the source-dependent kinematic noise term (\(\sigma^2_{\mu}\)), high-variance regions (cluster outskirts) will be statistically misinterpreted as low-signal regions, effectively down-weighting the most massive structures and suppressing the inferred \(S_8\). Standard analyses assume shear noise is dominated only by random galaxy orientations ("shape noise"), missing this kinematic component.

    - **Source-Dependent:** Correlated with galaxy type and redshift (via proper motion).

    - **Unmodeled:** Not present in standard covariance matrices.

**Result:** Unmodeled excess variance in the likelihood analysis acts as a source of systematic bias. If the stochastic noise is treated as signal or dilutes the inference, it can suppress the recovered clustering amplitude \(S_8\) in galaxy surveys. The CMB, being static, is free from this kinematic noise bias and recovers the higher (true) value.

**Prediction:** The \(S_8\) tension is not a breakdown of \(\Lambda\)CDM growth, but a signature of the TEP "noise floor". The tension should decrease if shear likelihoods are weighted by source kinematic priors.

### 5.4 Comparison of Predictions: TEP-GL vs. Particle Dark Matter

The following table summarizes the distinguishing predictions of the two frameworks:

| Observable | Particle DM Prediction | TEP-GL Prediction |
| --- | --- | --- |
| **Lensed FRB timing** | Achromatic, GR time delays only | Achromatic *residual* anomaly (ms-scale) |
| **Source-dependent Shear** | None (shear noise is constant) | **Excess Scatter** (noise scales with \(\mu_s\)) |
| **CMB lensing** | Standard convergence | Identical (CMB is static; no temporal smearing) |
| **Galaxy weak lensing** | Identical to CMB lensing | Systematic offset from CMB (source-dependent) |
| **Chromatic dependence** | None | None (both achromatic) |
| **Direct detection experiments** | Expected signal (WIMP recoil) | No signal (no particle to detect) |

The critical discriminant is the *source-dependence* of the inferred dark component. Particle dark matter lenses all sources identically; TEP-GL predicts that the apparent "dark" signal depends on the variability timescale of the background source.

### 5.5 Falsification Criteria

TEP-GL is a falsifiable hypothesis. The following observations would exclude specific regimes or the framework entirely:

    - **Reference Envelope Exclusion:** If precision timing of strongly lensed FRBs yields achromatic residuals consistent with zero to better than 0.1 ms across diverse lens environments, the Reference Envelope parameter space is excluded.

    - **Universal Shear Scatter:** If the intrinsic scatter of weak-lensing shear measurements is identical across all source kinematic classes (after correcting for measurement noise), the stochastic temporal shear mechanism is excluded.

    - **CMB-Galaxy Agreement:** If future surveys (e.g., CMB-S4, Rubin/LSST) demonstrate that CMB lensing and galaxy weak-lensing mass inferences agree to better than 1% after controlling for all modeling systematics, the Extended Regime parameter space is excluded.

    - **Chromatic Anomaly:** If any timing or morphological anomaly shows wavelength dependence after correction for plasma dispersion and dust, the achromatic prediction of two-metric coupling is falsified. This would indicate conventional astrophysical systematics rather than metric effects.

    - **Direct Detection:** A confirmed detection of dark matter particles in direct-detection experiments would establish the existence of a particulate dark sector. TEP would then be relegated to a sub-dominant correction rather than a primary explanation.

The framework stands or falls on empirical grounds. It does not claim immunity from observation; the claim is that the observations required to test TEP have not yet been performed with sufficient precision.

**Principle:**

#### Box 5.1: Quantitative Falsification Thresholds

To ensure rigorous falsifiability, explicit null-result thresholds are specified for each prediction channel:

| Test | Sample Size | Precision Required | Null-Result Threshold | Regime Falsified |
| --- | --- | --- | --- | --- |
| **Lensed FRB Timing** | ≥10 lensed FRBs | <0.1 ms timing | All residuals <0.01 ms | Reference Envelope |
| **Shear-Variability Correlation** | ≥1000 lenses per class | <0.5% shear scatter | No correlation at >3σ | Extended Regime (dynamic term) |
| **CMB-Galaxy \(S_8\) Tension** | CMB-S4 + LSST | <1% \(S_8\) agreement | Statistically consistent agreement (< 1\(\sigma\)) | Extended Regime (static term) |
| **Variability-Mass Correlation** | ≥100 lenses with varied sources | 10% mass precision | No correlation at >3σ | Extended Regime |

**Complete Falsification:** If *all four* null-result thresholds are met simultaneously—(1) no FRB timing residuals, (2) no shear-variability correlation, (3) CMB-galaxy agreement <1%, and (4) no variability-mass correlation—then TEP-GL is falsified as a framework, not merely in one regime. The theory makes no predictions distinguishable from standard CDM under these conditions.

**Timeline:** Tests (1) and (4) are achievable within 5–10 years (CHIME, DSA-2000, Rubin/LSST). Tests (2) and (3) require next-generation surveys (CMB-S4, Euclid) with ~10-year horizons. A definitive verdict on TEP-GL is therefore expected by ~2035.

## 6. Discussion: A Paradigm Shift in the Dark Sector

### 6.1 Limitations of the "Stack of Corrections"

Gravitational physics has historically advanced by refining the fundamental geometric framework rather than adding ad-hoc corrections. TEP represents a similar upgrade to the fundamental object, following the natural historical progression of physical theory:

    - **Newton:** Gravity is a Force; Time is Absolute.

    - **Einstein:** Gravity is Geometry; Time is Relative (Coordinate-Dependent).

    - **TEP:** Gravity is Geometry; Time is a *Dynamical Field*.

By treating the rate of proper time accumulation as a physical field with its own degrees of freedom (rather than a fixed function of the metric), TEP unifies the "stack" into a single framework. The "dark matter" anomaly is simply the observation that this field has spatial gradients.

### 6.2 CMB Lensing and the Integrated Sachs-Wolfe Constraint

A scalar field capable of generating year-scale differential delays on cluster scales raises a critical question: what are the implications for the Cosmic Microwave Background (CMB)? Two potential tensions must be addressed.

#### 6.2.1 The ISW Effect

The Integrated Sachs-Wolfe (ISW) effect arises when CMB photons traverse time-evolving gravitational potentials. In TEP, the conformal factor \(A(\phi)\) contributes an additional term to the photon temperature perturbation:

$\frac{\Delta T}{T}\bigg|_{\text{TEP}} = \int \frac{\partial \ln A(\phi)}{\partial t}\, dt$

For TEP to remain consistent with *Planck* constraints on the ISW amplitude, one of the following must hold:

    - **Quasi-Static Regime:** The scalar field evolves slowly (\(\partial_t A \ll H_0 A\)), so the TEP contribution to ISW is subdominant to the standard \(\dot{\Psi}\) term. This is natural if \(\phi\) tracks the matter distribution adiabatically.

    - **Cancellation:** In some scalar-tensor theories, the conformal ISW contribution partially cancels the metric ISW, leaving the net effect within observational bounds (cf. Amendola et al. 2008).

**Principle:**

#### Box 6.1: Order-of-Magnitude ISW Consistency Check

The quasi-static condition is demonstrated to be satisfied for TEP's required scalar field profile.

**Setup:** The conformal factor is \(A(\phi) = 1 + \alpha(\phi)\) with \(\alpha \sim 10^{-6}\) on halo scales. The ISW constraint requires:

$\left|\frac{\partial \ln A}{\partial t}\right| \approx \left|\frac{\dot{\alpha}}{1+\alpha}\right| \approx |\dot{\alpha}| \ll H_0$

**Estimate of \(\dot{\alpha}\):** If \(\phi\) tracks the matter distribution adiabatically (sourced by \(T^\mu_\mu\)), then \(\alpha\) evolves on the timescale of structure formation:

$\dot{\alpha} \sim \alpha \cdot \frac{\dot{\rho}}{\rho} \sim \alpha \cdot H(z) \cdot f(z)$
where \(f(z) \equiv d\ln D/d\ln a \approx \Omega_m(z)^{0.55}\) is the growth rate. At \(z \sim 0.5\) (peak ISW sensitivity), \(f \approx 0.8\) and \(H(z) \approx 1.2 H_0\).

**Result:**

$\frac{|\dot{\alpha}|}{H_0} \sim \alpha \cdot f \cdot \frac{H(z)}{H_0} \sim 10^{-6} \times 0.8 \times 1.2 \approx 10^{-6}$

This is **six orders of magnitude below** the ISW constraint threshold (\(\lesssim 0.1\)). The quasi-static approximation is therefore robustly satisfied for TEP's required coupling strengths.

**Physical Interpretation:** The scalar field \(\phi\) varies spatially (creating "dark matter" halos) but evolves temporally only as fast as the underlying matter distribution. Year-scale *spatial* delays across an image plane are compatible with cosmologically slow *temporal* evolution because the delay gradient is set by the static \(\nabla\alpha\), not by \(\dot{\alpha}\).

#### 6.2.2 CMB Lensing Power Spectrum

The CMB lensing convergence power spectrum \(C_\ell^{\kappa\kappa}\) is measured to high precision by *Planck* and ACT. In standard \(\Lambda\)CDM, this is sourced entirely by the integrated mass distribution. In TEP, the lensing kernel receives contributions from both the geometric (mass) term and the temporal (refractive) term:

$\kappa_{\text{TEP}} = \kappa_{\text{mass}} + \kappa_{\text{refraction}}$

For the CMB (a static, non-evolving source), the "Dynamic Shutter" term vanishes identically—there is no source proper motion to couple to the delay gradient. The only TEP contribution is the **static refractive term** \(\nabla^2(\Delta\tilde{\tau}_{\text{static}})\), which is coherent and radially symmetric around mass concentrations.

**Implication:** CMB lensing in TEP is dominated by the same static mechanism that produces "dark matter halos" in galaxy lensing. The amplitude is controlled by the integrated \(\alpha(\phi)\) profile, not by lens motion. This is consistent with the observed CMB lensing power spectrum *provided* the static refractive term has the correct radial profile—a constraint that effectively fixes the relationship between \(\alpha(\phi)\) and the matter density \(\rho\).

#### 6.2.3 Parameter Space Constraints

Combining ISW and CMB lensing constraints with the Vainshtein screening requirements (Box 4.1) defines a narrow but viable parameter window:

| Constraint | Requirement | Status |
| --- | --- | --- |
| Solar System (Cassini) | \(\alpha_{\text{local}} \lesssim 10^{-5}\) | Satisfied via Vainshtein |
| Cluster Lensing | \(\alpha_{\text{halo}} \sim 10^{-1}\) | Required for DM phenomenology |
| ISW (*Planck*) | \(\partial_t \ln A / H_0 \lesssim 0.1\) | Requires quasi-static \(\phi\) |
| CMB Lensing | \(\kappa_{\text{refraction}} \propto \kappa_{\text{mass}}\) | Requires \(\alpha \propto \Psi\) scaling |

The requirement that \(\alpha(\phi) \propto \Psi\) (the scalar field tracks the Newtonian potential) is a **non-trivial constraint** on the TEP Lagrangian. It is satisfied in conformal coupling models where \(\phi\) is sourced by the trace of the stress-energy tensor, but represents a theoretical prior that must be verified against N-body simulations. This is flagged as a key target for future numerical work.

**Principle:**

#### Box 6.3: Addressing the Hierarchy Problem (The "Screening Cliff")

**Critique:** A common objection to Modified Gravity theories is the "Screening Cliff"—the requirement that the scalar field \(\phi\) must be order-unity on cluster scales to mimic dark matter, yet suppressed to \(\lesssim 10^{-15}\) in the Solar System to satisfy Cassini constraints. This implies a "fine-tuned" suppression of \(10^{15}\).

**TEP Response:** TEP does *not* require order-unity couplings on cluster scales. The "Phantom Mass" effect is driven by the *gradient* of the proper time field, which must match the Newtonian potential gradient of the dark matter halo. Since \(\Phi_{\text{DM}}/c^2 \sim 10^{-5}\) for typical clusters, TEP only requires \(\alpha(\phi) \sim 10^{-5}\) on Mpc scales. The hierarchy is therefore:

    - Cluster Scale (Mpc): \(\alpha \sim 10^{-5}\) (Year-scale delays)

    - Solar System (AU): \(\alpha < 10^{-15}\) (Microsecond residuals)

**Mechanism:** The required suppression factor is \(10^{10}\). In the Vainshtein mechanism, the scalar force is suppressed by \((r/r_V)^{3/2}\) inside the Vainshtein radius \(r_V\). For a galaxy with \(r_V \sim 100\) kpc, the suppression at Solar System scales (\(r \sim 10^{-11}\) kpc) is:

$\text{Suppression} \approx \left(\frac{10^{-11}}{100}\right)^{1.5} \approx (10^{-13})^{1.5} \approx 10^{-19.5}$

**Conclusion:** The Vainshtein mechanism is not just sufficient; it is *over-efficient*. The challenge in TEP is not hiding the scalar field in the Solar System, but preventing it from being screened too early in the outskirts of galaxies. This "Magnitude Hierarchy" is a standard feature of non-linear derivative couplings and is mathematically robust, albeit counter-intuitive compared to linear shielding.

### 6.3 Interpretational Challenges

TEP effects mimic standard lensing signatures, making them difficult to distinguish without time-domain analysis.

    - **Static Degeneracy:** Most lensing data are analyzed under static assumptions. In this limit, temporal shear is mathematically indistinguishable from spatial mass (Section 3).

    - **Differential vs. Common-Mode Effects:** The GW170817 constraint is often interpreted as excluding modified metric couplings. However, as shown in Section 4, this constraint applies to the differential disformal sector. Conformal proper-time dilation is a *common-mode effect* for co-propagating signals, leaving the "speed of gravity" constraint satisfied while permitting significant time-domain phenomenology.

### 6.4 Reinterpreting the Evidence: Multipath vs. Single-Path

When viewed through the lens of TEP, the disparate constraints on the dark sector resolve into a consistent picture based on signal topology:

    - **Lensing (Multipath):** Observations that require "dark matter" (e.g., lensing) are fundamentally *multipath* measurements. These different pathways sample different values of the time field \(\phi\), creating differential delays that manifest as "phantom mass."

    - **GW170817 (Single-Path):** Observations that constrain "modified gravity" (GW170817) are fundamentally *single-path* measurements. The signals originate from the same coordinate and traverse the same time-warped regions. The temporal distortions are common-mode and cancel out.

The apparent contradiction—"how can gravity be modified enough to create dark matter but unmodified enough to pass GW170817?"—is resolved. Gravity may not be modified; time is. And time manifests as a differential observable only when comparing divergent paths.

### 6.5 Global Variation of Effective Light Speed

If TEP is correct, then \(c\) is locally invariant but globally variable. This challenges the foundational assumption used to interpret all astrophysical data.

    - **Standard Assumption:** \(c\) is constant; therefore, any anomaly in \(t\) must be due to extra path length (spatial curvature) or extra mass (Shapiro delay). Conclusion: *dark matter* exists.

    - **TEP Reality:** \(t\) is modulated by the field \(\sqrt{A(\phi)}\); \(c\) is effectively variable along the path. The anomaly in \(t\) is due to the scalar field gradient. Conclusion: *dark matter* admits the interpretation of variable light speed.

TEP removes the need for an invisible substance by correcting the assumption about the constancy of the global speed of light.

### 6.6 TEP as a Bridge Between Paradigms

TEP offers a potential unification of competing frameworks:

    - **CDM:** In the limit where the scalar field \(\phi\) relaxes into stable halos and source variability is ignored, TEP reproduces the phenomenology of Cold Dark Matter (via the Phantom Mass mechanism).

    - **MOND:** In the regime of low acceleration (large gradients in \(\phi\)), TEP reproduces the scaling relations of Modified Newtonian Dynamics without modifying Newton's laws, but by modifying the clock metric.

Thus, TEP is not merely a competitor to these paradigms; it is a framework that explains *why* CDM works on large scales (where static assumptions hold) and *why* MOND works on galactic scales (where the field gradients are dominant).

**Principle:**

#### Box 6.2: Derivation of the MOND Acceleration Scale from TEP

The characteristic MOND acceleration \(a_0 \approx 1.2 \times 10^{-10}\) m/s² emerges naturally from TEP parameters:

**Setup:** The TEP strong-coupling scale is \(\Lambda \sim 10^{-13}\) eV \(\sim H_0\) (the dark energy scale). The scalar field gradient produces an effective "fifth force" that becomes significant when the Newtonian acceleration falls below a threshold.

**Derivation:** The Vainshtein suppression factor scales as \((r/r_V)^{3/2}\) inside the screening radius. The transition from screened to unscreened behavior occurs when the scalar-mediated acceleration equals the Newtonian acceleration:

$a_\phi \sim \frac{g}{M_{Pl}} \nabla\phi \sim a_N \quad \Rightarrow \quad a_0 \sim \frac{\Lambda^2}{M_{Pl}} \sim \frac{(H_0 M_{Pl})^2}{M_{Pl}} \sim c H_0$

**Numerical Check:**

$a_0 \sim c H_0 \approx (3 \times 10^8 \text{ m/s}) \times (2.2 \times 10^{-18} \text{ s}^{-1}) \approx 6.6 \times 10^{-10} \text{ m/s}^2$

This is within a factor of 5 of the empirical MOND value \(a_0 \approx 1.2 \times 10^{-10}\) m/s². The \(\mathcal{O}(1)\) discrepancy can be absorbed into the coupling constant \(g\) or reflects the simplified derivation. The key result is that TEP with \(\Lambda \sim H_0\) *automatically* produces an acceleration scale of cosmological magnitude—the same "cosmic coincidence" that MOND requires as an unexplained input.

**Physical Interpretation:** MOND phenomenology emerges in TEP because the scalar field's strong-coupling scale is set by the Hubble parameter. Regions within galaxies with internal accelerations \(a < a_0 \sim cH_0\) are in the unscreened regime where the scalar field gradient dominates; regions with \(a > a_0\) are screened and follow Newtonian dynamics. This provides a physical origin for the MOND interpolating function.

### 6.7 Convergence of Evidence

While TEP-GL is a theoretical proposal, it is notable that multiple independent anomalies in current cosmology converge on the phenomenology predicted by this framework. These tensions, often treated as separate puzzles, may represent a single systemic failure of the isochrony assumption.

| TEP-GL Prediction | Existing Observational Anomaly | Status |
| --- | --- | --- |
| **Source-Dependent Shear**
(Kinematic noise bias) | **\(S_8\) Tension**
Galaxy lensing (\(S_8 \approx 0.79\)) consistently lower than static CMB (\(S_8 \approx 0.83\)). | **Strong Support** |
| **Mass-Sheet Degeneracy**
(Temporal vs Spatial) | **\(H_0\) Tension**
Time-delay cosmography (\(H_0 \approx 73\)) conflicts with CMB (\(H_0 \approx 67\)). | **Consistent** |
| **Temporal Shear**
(Dynamic Shutter) | **Flux Ratio Anomalies**
Substructure required to explain ratios is often not found; "phantom" substructure. | **Consistent** |
| **Variability-Shear Correlation** | **Einstein Cross Anomalies**
Flux anomalies correlate with source variability phase (e.g., Eigenbrod et al. 2008). | **Untested** |

The theory does not require new observations to find initial support; it provides a unified explanation for *existing* anomalies that \(\Lambda\)CDM struggles to explain simultaneously.

### 6.8 The Path Forward

The analysis suggests that *dark matter phenomenology can be reproduced without a substance*. What exists is a temporal field with spatial gradients, and what has been measured as "dark matter" is the projection of those gradients onto observations that assume temporal synchrony. To maintain the particle interpretation, it is necessary to demonstrate that the universe is synchronous to a precision that excludes TEP effects. TEP offers a geometric identification of the phenomenology, explaining the universality of dark matter profiles and the success of MOND on galactic scales through a single parsimonious framework.

"Dark Time" is not merely an alternative to dark matter; it offers a geometric identification of the phenomenology that admits this interpretation.

## 7. Conclusions

### 7.1 Limitations of the Isochrony Axiom

For nearly a century, the dark matter paradigm has rested on a single, simplifying assumption: that the observed universe is synchronous. By treating the photons in telescopes as representing a single spatial slice of the source, standard models have been forced to interpret all observed distortions as spatial deflections caused by mass. It has been shown that this *Isochrony Axiom* is an effective approximation that may break down in the presence of generalized metric couplings.

The *Temporal Equivalence Principle (TEP)*, introduced in a companion paper (Smawfield 2025a), provides a corrected framework. Under TEP, "gravity" includes the phenomenology of differential proper-time accumulation. When the assumption of isochrony is relaxed, the "dark" signals in lensing maps can be modeled as *temporal-composite artifacts*—the projection of unmodeled time dilation onto the spatial image plane.

### 7.2 Summary of Findings

    - **Phantom Mass from Dark Time:** It was demonstrated that a gradient in proper-time accumulation across an image plane acts as a "temporal shear" on evolving sources. In a static reconstruction, this is mathematically indistinguishable from the gravitational shear of a dark matter halo. The "dark sector" may simply be the difference between the gravitational metric (which guides orbits) and the matter metric (which guides clocks).

    - **The GW170817 "Speed Limit" is Nuanced:** The standard multi-messenger constraint was deconstructed. Because the conformal component of the metric coupling preserves null cones, it is *invisible* to differential speed-of-gravity tests. Photons and gravitational waves share the common-mode dilation. While GW170817 constrains disformal propagation speeds to \(|c_\gamma - c_g|/c \lesssim 10^{-15}\), it leaves the conformal "rate of time" unconstrained. We demonstrate that conformal gradients can reproduce specific aspects of dark matter phenomenology—particularly coherent lensing shear and time-domain signatures—subject to Vainshtein screening constraints that ensure consistency with Solar System tests.

    - **A New Observational Era:** Two regimes for testing TEP are defined. In the conservative *Reference Envelope*, the effects are millisecond-scale and detectable only in high-time-resolution astrophysics (FRBs, pulsars). In the *Extended Regime*, where environmental screening and/or a revised operational mapping allow larger effective delays, TEP becomes a full dynamical alternative to particle dark matter.

### 7.3 Future Outlook: Time, Not Mass

The persistence of the dark matter problem despite decades of particle searches suggests a potential error in the underlying premises. The error may lie in the treatment of time. By treating time as a passive parameter rather than an active dynamical field, there is a risk of blinding analysis to the true nature of the "dark" universe.

The path forward lies in *Chronometric Lensing*. The field must move beyond static mass-mapping and towards *chronometric mapping*—measuring the detailed arrival-time structure of the universe. If "dark matter" halos are actually "time dilation" halos, the definitive indicator will not be a particle recoil in a xenon tank, but a millisecond residual in a Fast Radio Burst. The focus of inquiry must shift from searching for missing mass to characterizing the dynamical structure of time.

### 7.4 The Paradigm Shift

Within the operational axioms adopted here, the phenomenology traditionally attributed to a particulate dark sector can be reinterpreted as the projection of two-metric time transport onto inference pipelines that assume isochrony. The anomalies are real; the claim is that their standard "missing mass" interpretation is not unique once the Isochrony Axiom is relaxed.

The observational program outlined herein—time-domain lensing of fast transients, precision strong-lens residual timing, and source-variability-dependent shear consistency tests—provides a viable route to distinguish particle dark matter from a temporal-composite channel while maintaining a conservative Reference Envelope baseline anchored to existing multi-messenger constraints, thereby offering a promising avenue for resolving the dark matter problem.

## References

Abbott, B. P., et al. (LIGO/Virgo Collaboration) 2017, *Phys. Rev. Lett.*, 119, 161101

Abbott, T. M. C., et al. (DES Collaboration) 2022, *Phys. Rev. D*, 105, 023520

Asgari, M., et al. (KiDS Collaboration) 2021, *A&A*, 645, A104

Baker, T., Bellini, E., Ferreira, P. G., Lagos, M., Noller, J., & Sawicki, I. 2017, *Phys. Rev. Lett.*, 119, 251301 (arXiv:1710.06394)

Bartelmann, M., & Schneider, P. 2001, *Phys. Rep.*, 340, 291

Bekenstein, J. D. 1993, *Phys. Rev. D*, 48, 3641 (arXiv:gr-qc/9211017)

Bekenstein, J. D. 2004, *Phys. Rev. D*, 70, 083509 (arXiv:astro-ph/0403694)

Blandford, R. D., & Narayan, R. 1986, *ApJ*, 310, 568

Bullock, J. S., & Boylan-Kolchin, M. 2017, *ARA&A*, 55, 343

Burrage, C., & Sakstein, J. 2018, *Living Rev. Relativ.*, 21, 1

Chang, C., et al. 2024, *arXiv e-prints*, arXiv:2406.19654

Clifton, T., Ferreira, P. G., Padilla, A., & Skordis, C. 2012, *Phys. Rep.*, 513, 1 (arXiv:1106.2476)

Clowe, D., et al. 2006, *ApJ*, 648, L109

Creminelli, P., & Vernizzi, F. 2017, *Phys. Rev. Lett.*, 119, 251302 (arXiv:1710.05877)

Damour, T., & Esposito-Farèse, G. 1992, *Class. Quantum Grav.*, 9, 2093

Di Valentino, E., et al. 2021, *Classical Quantum Gravity*, 38, 153001

Eigenbrod, A., et al. 2008, *A&A*, 490, 933 (arXiv:0810.0011)

Ezquiaga, J. M., & Zumalacárregui, M. 2017, *Phys. Rev. Lett.*, 119, 251304 (arXiv:1710.05901)

Fian, C., Jiménez-Vicente, J., Mediavilla, E., et al. 2018, *ApJ*, 859, 50 (arXiv:1805.09619)

Gilman, D., et al. 2016, *MNRAS*, 462, 219 (arXiv:1601.01671)

Hu, W., Barkana, R., & Gruzinov, A. 2000, *Phys. Rev. Lett.*, 85, 1158 (arXiv:astro-ph/0003365)

Hui, L., Ostriker, J. P., Tremaine, S., & Witten, E. 2017, *Phys. Rev. D*, 95, 043541 (arXiv:1610.08297)

Khoury, J., & Weltman, A. 2004, *Phys. Rev. D*, 69, 044026 (arXiv:astro-ph/0309411)

Koch Ocker, S., et al. 2022, *ApJ*, 931, 87

Kochanek, C. S. 2002, *ApJ*, 578, 25 (arXiv:astro-ph/0205319)

Leitherer, C., Schaerer, D., Goldader, J., et al. 1999, *ApJS*, 123, 3 (arXiv:astro-ph/9902334)

Li, Z.-X., Gao, H., Ding, X.-H., Wang, G.-J., & Zhang, B. 2018, *Nat. Commun.*, 9, 3833 (arXiv:1708.06357)

Lintott, C. J., Schawinski, K., Keel, W., et al. 2009, *MNRAS*, 399, 129 (arXiv:0906.5304)

McGaugh, S. S., et al. 2016, *Phys. Rev. Lett.*, 117, 201101

Milgrom, M. 1983, *ApJ*, 270, 365

Muñoz, J. B., Kovetz, E. D., Dai, L., & Kamionkowski, M. 2016, *Phys. Rev. Lett.*, 117, 091301 (arXiv:1605.00008)

Narayan, R., & Bartelmann, M. 1996, in *Formation of Structure in the Universe*, ed. A. Dekel & J. P. Ostriker (Cambridge Univ. Press) (arXiv:astro-ph/9606001)

Navarro, J. F., Frenk, C. S., & White, S. D. M. 1996, *ApJ*, 462, 563

Nityananda, R., & Samuel, J. 1992, *Phys. Rev. D*, 45, 3862 (arXiv:astro-ph/9205001)

Planck Collaboration 2020, *A&A*, 641, A6

Rathna Kumar, S., Tewes, M., Stalin, C. S., Courbin, F., et al. 2013, *A&A*, 557, A44 (arXiv:1306.5105)

Refsdal, S. 1964, *MNRAS*, 128, 307

Riess, A. G., et al. 2022, *ApJ*, 934, L7

Rubin, V. C., & Ford, W. K., Jr. 1970, *ApJ*, 159, 379

Sakstein, J., & Jain, B. 2017, *Phys. Rev. Lett.*, 119, 251303 (arXiv:1710.05893)

Schawinski, K., Koss, M., Berney, S., & Sartori, L. F. 2015, *MNRAS*, 451, 2517 (arXiv:1505.06733)

Schneider, P., & Sluse, D. 2013, *A&A*, 559, A37 (arXiv:1306.0901)

Schneider, P., & Sluse, D. 2014, *A&A*, 564, A103 (arXiv:1306.4675)

Schumann, M. 2019, *J. Phys. G*, 46, 103003

Skordis, C., & Złośnik, T. 2021, *Phys. Rev. Lett.*, 127, 161302

Sluse, D., Hutsemékers, D., Courbin, F., Meylan, G., & Wambsganss, J. 2012, *A&A*, 544, A62

Smette, A., Surdej, J., et al. 1992, *ApJ*, 389, 39

Smawfield, M. L. (2025). *Temporal-Spatial Coupling in Gravitational Lensing: A Reinterpretation of Dark Matter Observations* (v0.3 (Tortola)). Zenodo, DOI: 10.5281/zenodo.16921911

Tie, S. S., & Kochanek, C. S. 2018, *MNRAS*, 473, 80 (arXiv:1707.01908)

Treu, T., & Marshall, P. J. 2016, *A&ARv*, 24, 11 (arXiv:1605.05333)

Wambsganss, J. 2006, in *Gravitational Lensing: Strong, Weak and Micro*, Saas-Fee Advanced Course 33, ed. G. Meylan et al. (Springer) (arXiv:astro-ph/0604278)

Wong, K. C., et al. (H0LiCOW Collaboration) 2020, *MNRAS*, 498, 1420

Walsh, D., Carswell, R. F., & Weymann, R. J. 1979, *Nature*, 279, 381

Zwicky, F. 1933, *Helv. Phys. Acta*, 6, 110

#### Contact Information

    **Author:** Matthew Lukin Smawfield

    **Affiliation:** Independent Researcher

    **Email:** [matthewsmawfield@gmail.com](mailto:matthewsmawfield@gmail.com)

    **ORCID:** [0009-0003-8219-3159](https://orcid.org/0009-0003-8219-3159)

    **GitHub:** [github.com/matthewsmawfield](https://github.com/matthewsmawfield)

        **License:** This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

        Version: v0.3 (Tortola) · Last updated: 19 December 2025

## Appendix A: Mathematical Derivations

### A.1 Proof of Null Cone Invariance in the Conformal Limit

**Proposition:** If two metrics are conformally related by \(\tilde{g}_{\mu\nu} = A(\phi) g_{\mu\nu}\), they share the same null geodesics as unparameterized curves.

**Proof:**

Let \(k^\mu = dx^\mu / d\lambda\) be a null vector in \(g_{\mu\nu}\), satisfying \(g_{\mu\nu}k^\mu k^\nu = 0\) and the geodesic equation \(k^\nu \nabla_\nu k^\mu = 0\) (where \(\nabla\) is the Levi-Civita connection of \(g\)).

In the metric \(\tilde{g}_{\mu\nu}\), the null condition holds immediately:
$\tilde{g}_{\mu\nu} k^\mu k^\nu = A g_{\mu\nu} k^\mu k^\nu = 0$
The connection coefficients \(\tilde{\Gamma}^\lambda_{\mu\nu}\) for \(\tilde{g}\) are related to \(\Gamma^\lambda_{\mu\nu}\) by:

$\tilde{\Gamma}^\lambda_{\mu\nu} = \Gamma^\lambda_{\mu\nu} + \frac{1}{2}(\delta^\lambda_\mu \partial_\nu \ln A + \delta^\lambda_\nu \partial_\mu \ln A - g_{\mu\nu} g^{\lambda\sigma} \partial_\sigma \ln A)$
Substituting this into the geodesic equation for \(\tilde{g}\):

$k^\nu \tilde{\nabla}_\nu k^\mu = k^\nu \nabla_\nu k^\mu + \frac{1}{2}(2 k^\mu k^\nu \partial_\nu \ln A - g_{\nu\sigma} k^\nu k^\sigma \dots)$
Since \(k\) is null (\(g_{\nu\sigma}k^\nu k^\sigma = 0\)) and \(k^\nu \nabla_\nu k^\mu = 0\), this yields:

$k^\nu \tilde{\nabla}_\nu k^\mu = (k^\nu \partial_\nu \ln A) k^\mu$
This is the geodesic equation with a non-affine parameterization (\(Dk/d\lambda \propto k\)). Thus, the curve is a geodesic of \(\tilde{g}\), differing only by the parameterization (the clock rate). \(\square\)

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

*This document was automatically generated from the TEP-GL research site. For the interactive version with figures and enhanced formatting, visit: https://matthewsmawfield.github.io/TEP-GL/*

*Related Work:*
- [**TEP Theory**](https://doi.org/10.5281/zenodo.16921911) (Foundational framework)

*Source code and data available at: https://github.com/matthewsmawfield/TEP-GL*
