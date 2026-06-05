# Temporal Equivalence Principle: Disformal Kinematics and the Measurement Landscape
**Matthew Lukin Smawfield**
Version: v0.1 (Kuala Lumpur)
First published: 24 May 2026 · Last updated: 5 June 2026

---

## Abstract

In the unscreened regime of the Temporal Equivalence Principle, this paper proposes that virtual force carriers and statistical wavefunctions are tangent-limit descriptions of a deeper disformal geometry of the temporal field, governed by the coupling coefficient B(φ). In the screened limit, where the local interaction energy density substantially exceeds the saturation scale ρ_{T} ≈ 20 g/cm^{3}, the observable disformal response is suppressed and standard perturbation theory is recovered as the tangent limit. U(1) electromagnetism is derived from temporal vortex defects in a compact phase bundle, with the photon interpreted as a propagating shear-wave in the dilute-defect limit. The SU(2) sector is currently speculative: no gauge-covariant derivation exists and no claim is made about replacing the weak interaction geometrically. Geometric kinematics is developed by showing how attraction and repulsion can arise when charged defects navigate spatially varying temporal shear, with Coulomb-like behavior recovered when source boundary conditions impose the corresponding weak-field B(φ) gradient. Entanglement is modeled as an unbroken macroscopic geometric contour connecting bifurcated topological charges, with Schmidt rank mapped to contour winding number. The probability wavefunction is reinterpreted as the physical volumetric shear-wake of proper-time phase transport; the double-slit and tunneling phenomena are described as interference and barrier-navigation effects of that wake. A re-analysis of published graphene Aharonov-Bohm interferometry data (Zimmermann et al., Nat. Commun. 8, 14983, 2017) is presented as a candidate empirical test. Phase extraction was performed by lock-in quadrature demodulation against the known carrier, a method that avoids the nonlinear aliasing of the Hilbert transform. A disformal topography regression, in which the measured phase shift is regressed directly against the metric tilt B(φ), was evaluated with the corrected standard BIC formula. The Gaussian confinement-peak model is overwhelmingly preferred (BIC = -498.86), with the linear gradient model ranking second (BIC = -467.48) and the exponential boundary-decay model third (BIC = -453.83). The harmonic periodic-shear model is strongly disfavoured (BIC = -107.69). Bayesian model averaging assigns the Gaussian model a posterior probability of ≈ {{tep.topography.bma_weight_B_gaussian}}, indicating an unambiguous preference. The uniform temporal-dilation model (γ ≠ 1) fails strongly (BIC = -87.22). Synthetic null tests with the lock-in pipeline show that exponential ground truth is recovered perfectly (true-positive rate ≈ {{step06.tp_exponential}}), but Gaussian and harmonic models remain partially conflated on synthetic data with strong amplitude modulation. These findings underscore the need for replication across independent devices. A comprehensive control suite was performed on the medium-smoothed transmission profile: dataset provenance was verified by SHA-256; the original MATLAB processing was reproduced exactly; seven standard nuisance models were fitted; P-gamma covariance was estimated by Hessian and bootstrap; a train/test split was performed; and the model was refitted on raw, light, medium, and heavily smoothed data. The γ parameter is degenerate with an effective period P_{eff} = P/γ and is unstable across smoothing levels, ranging from 0.864 (raw) to 0.845 (medium), with 1.067 (light) and 1.011 (heavy), confirming that the uniform-dilation proxy is not the correct observable for this geometry.

Keywords: disformal kinematics, measurement, entanglement, Aharonov-Bohm, graphene interferometry, light-cone geometry, virtual bosons, temporal equivalence principle, shear-wake, geometric kinematics

## 1. Introduction: The Illusion of the Dead Vacuum

### 1.1 Virtual Bosons as Geometric Fiction

The standard model uses virtual bosons as the perturbative bookkeeping for forces between particles. In the TEP framework, this paper asks whether that bookkeeping can be recovered as the screened tangent limit of a deeper continuous light-cone geometry tilted by the disformal coupling B(φ). In the screened limit, where the local interaction energy density substantially exceeds the critical saturation density *ρ_{c} ≈ 20 g/cm^{3}* (derived as the macroscopic temporal saturation limit ρ_{T} in TEP-UCD, Paper 6), the observable disformal response is suppressed and the light cone recovers its isotropic form; standard perturbation theory is recovered as the tangent limit. In the unscreened regime, the disformal structure becomes the candidate interaction ontology.

### 1.2 Reinterpreting the Copenhagen Interpretation

Historically, the resistance to the Copenhagen Interpretation was most famously articulated by Einstein's assertion that nature is fundamentally deterministic, rather than probabilistic. The TEP framework formally vindicates this intuition. What the standard model interprets as fundamental statistical indeterminism — such as wavefunction collapse and complementarity — is reinterpreted in TEP merely as a symptom of assuming a flat, isochronous background. When the assumption of universal parameter time is dropped, "measurement" is no longer a probabilistic dice roll, but is instead modeled as a deterministic geometric interaction of a probe with the local temporal shear field.

The Copenhagen Interpretation's statistical indeterminism — wavefunction collapse, complementarity, the measurement problem — is reinterpreted in TEP as a symptom of assuming a flat, isochronous background. When the background is dynamical, "measurement" is modeled as the geometric interaction of a probe with the local temporal shear field.

## 2. Routing Interactions via Disformal Coupling

### 2.1 Gauge Symmetries as Candidate Light-Cone Tilts

In the unscreened regime, where the local interaction energy density is well below the saturation scale *ρ_{T}*, internal gauge symmetries are modeled as candidate disformal light-cone tilts governed by the coupling coefficient B(φ). In the screened limit the observable disformal response is suppressed and interactions become isotropic. The minimal metric ansatz

g̃_{μν} = A^{2}(φ) g_{μν} + B(φ) ∂_{μ}φ ∂_{ν}φ

encodes the interaction geometry considered in this paper through the scalar function B(φ). The matter metric g̃_{μν} encodes the causal structure to which all non-gravitational fields couple; the gravitational metric g_{μν} describes spacetime curvature. The conformal factor A(φ) governs local length-scale rescaling; in the present interaction analysis it is absorbed into the background and does not affect the interference phase directly. A complete gauge replacement, especially for non-Abelian sectors, may require additional internal orientation variables or multiplet structure beyond this minimal ansatz. *Metric signature convention:* (+, −, −, −) throughout.

While multi-messenger constraints (GW170817) require the effective disformal term $B(\phi)(\partial\phi)^2$ to be phenomenologically negligible for signals propagating across deep intergalactic voids, the subatomic environment is entirely different. At the femtometer scale of the topological charge, the temporal gradient $(\partial\phi)^2$ is immense. This extreme local gradient drives quantum kinematics via $B(\phi)$ at the particle scale, while naturally vanishing as $\nabla\phi \to 0$ in the late universe, remaining compatible with Paper 0's macroscopic propagation bounds.

For U(1) electromagnetism the TEP replacement program requires a first-principles derivation of Maxwell theory from temporal geometry, not merely a consistency construction. The starting point is a temporal phase bundle χ on the temporal manifold, a compact U(1) fibre over spacetime. The scalar field φ is the bundle magnitude; χ is the compact phase. The topological charge (Paper 24, TEP-SPIN) carries a quantised vortex in χ, with circulation ∮ ∇χ · d**l** = 2πn around the defect core. This U(1) phase bundle is introduced here as distinct from the SU(2) orientation bundle that carries spin in TEP-SPIN; the two bundles coexist over the same spacetime manifold.

A_{μ} = ∂_{μ}χ − (singular defect contribution),

where the gauge field A_{μ} is defined as the smooth part of the temporal phase gradient after removing the singular defect contribution. Because A_{μ} is not globally exact — it has non-trivial holonomy around defect cores — its exterior derivative is non-zero. The electromagnetic field strength is therefore

F_{μν} = ∂_{[μ}A_{ν]},

which is the curvature of the temporal phase bundle in the dilute-defect limit. The photon is the propagating shear-wave of the temporal phase field on the disformal metric. This derivation connects directly to the quantised circulation of TEP-SPIN and provides a geometric origin for electromagnetism that does not require ad hoc boundary conditions. The complete derivation, including the propagating mode equation and coupling to matter currents, is reserved for a forthcoming companion paper.

For SU(2) weak interactions no gauge-covariant derivation currently exists. The expression previously written,

F^{a}_{μν} = ∂_{[μ}(B^{a} ∂_{ν]}φ^{a}) + g ε^{abc} B^{b}B^{c} ∂_{μ}φ^{b} ∂_{ν}φ^{c},

uses ordinary derivatives rather than covariant derivatives and does not transform under the adjoint representation. It is algebraic mimicry without the differential geometry of gauge theory. Until a proper covariant derivation using D_{μ}^{ab} = δ^{ab}∂_{μ} + gε^{acb}W_{μ}^{c} and G^{a}_{μν} = [D_{μ}, D_{ν}]^{a}/ig is completed, the SU(2) sector must be regarded as entirely speculative. The manuscript therefore makes no claim about replacing the weak interaction geometrically.

### 2.2 The Photon as Proper-Time Phase Wake

In the proposed TEP ontology, the photon is not fundamental as a force carrier; it is the proper-time phase wake broadcast by an accelerating charge. The wake propagates along the disformally tilted light cone, carrying geometric information about the source's temporal environment. In the screened tangent limit, the usual perturbative photon description is recovered as the effective bookkeeping used by quantum electrodynamics.

### 2.3 Geometric Kinematics: Attraction and Repulsion as Navigation

Particles attract or repel because they are navigating localized, direction-dependent light-cone tilts. Consider two charges q_{1} and q_{2} separated by distance r in a background with spatially varying B(φ). The disformal metric introduces an effective refractive index for temporal propagation:

n_{eff}(r) = √(1 + B(φ(r)) |∇φ|^{2} / A^{2}(φ)).

A charge of opposite sign to the source experiences a region where n_{eff} increases toward the source, bending its trajectory inward: this is the geometric origin of Coulomb attraction. The trajectory equation follows from the geodesic equation in the disformal metric:

d^{2}x^{μ} / dτ^{2} + Γ^{μ}_{αβ} (dx^{α}/dτ)(dx^{β}/dτ) = −½ g^{μν} ∂_{ν} ln B(φ) (g_{αβ} − u_{α}u_{β}) (dx^{α}/dτ)(dx^{β}/dτ),

where Γ^{μ}_{αβ} is the Levi-Civita connection of the full metric g_{μν}. The right-hand side is a geometric force arising from the gradient of B(φ). For two like charges the gradient is repulsive (particles are deflected away from regions of higher B), while for opposite charges the effective sign of the coupling inverts and the gradient becomes attractive. No virtual photons are introduced as fundamental carriers in this geometric description; the standard virtual-photon expansion is recovered as the screened tangent-limit description.

The radial acceleration between two point charges in the weak-field limit (B ≪ 1) is

a_{r} = −q_{1}q_{2} ∇B(φ) · r̂,

which reproduces Coulomb's law when the source equation or boundary condition fixes the disformal potential so that ∇B = −(1/4πε_{0}) r̂ / r^{2} for a static charge distribution. This is a weak-field consistency condition: the geometric ansatz can recover the Coulomb form once the charge defect supplies the correct radial B(φ) gradient. For like charges (q_{1}q_{2} > 0) the acceleration is positive in the r̂ direction, corresponding to repulsion; for opposite charges the sign inverts, giving attraction. No virtual photons are introduced as fundamental carriers in this geometric description.

### 2.4 Entanglement Geometry

Entanglement is defined mathematically as an unbroken macroscopic geometric contour connecting bifurcated topological charges. When a particle pair is created, the temporal field at the creation event carries a single connected phase contour. As the particles separate, each carries with it a phase singularity; the contour between them remains unbroken because the disformal metric retains memory of the shared origin.

Let the two-particle state be represented by two phase fields φ_{1}(x) and φ_{2}(x) centred at positions x_{1} and x_{2}. The entanglement contour is the geodesic γ(s) in the temporal shear field that connects the two charge cores. Its length in proper time is

L_{ent} = ∫_{γ} dτ = ∫_{γ} √(B(φ) |dφ/ds|^{2}) ds.

The contour is unbroken when L_{ent} is finite and the integrand never vanishes along the path. A measurement on particle 1 samples the phase φ at x_{1}, which perturbs the entire contour because B(φ) couples the local phase to the metric. The perturbation propagates along the contour at the speed of temporal shear (the local light-cone tilt), not instantaneously, but because the contour is pre-existing and connected the correlation appears non-local in any isochronous slicing. There is no collapse of a wavefunction; there is only the geometric sampling of a shared, continuous temporal deformation.

The Schmidt rank of a bipartite entangled state maps to the winding number of the phase contour around the charge pair. A maximally entangled Bell state corresponds to a contour with winding number ±1, while partially mixed states correspond to contours with multiple windings or decohered segments where the integrand intermittently vanishes.

### 2.5 Bell Correlation from Contour Holonomy

The entanglement contour picture must reproduce the quantum correlation function if it is to be more than a geometric analogy. Consider a spin-1/2 singlet state prepared at the origin. In the TEP framework the two particles are topological charges in the temporal orientation bundle, each carrying a local spinor frame. The measurement at Alice's location projects the local spinor onto her chosen axis **a**; the measurement at Bob projects onto axis **b**.

The spinor is parallel-transported from the creation event to each detector along the shared geodesic γ of the temporal orientation bundle. For a maximally entangled state, the two transported spinors are anti-aligned when the measurement axes are parallel. The correlation function is determined by the SU(2) holonomy of the orientation bundle around the loop formed by the two measurement frames and the shared contour:

C(**a**, **b**) = −cos θ_{ab},

where θ_{ab} is the holonomy angle of the orientation transport. This equals the quantum mechanical correlation because the SU(2) holonomy of the spin-1/2 representation gives exactly the rotation matrix R_{ij}(θ) = δ_{ij} cos θ + ε_{ijk} n̂_{k} sin θ + (1 − cos θ) n̂_{i} n̂_{j}. The CHSH parameter S = |E(**a**,**b**) − E(**a**,**b**′) + E(**a**′,**b**) + E(**a**′,**b**′)| achieves its quantum maximum S_{max} = 2√2 when the four measurement directions are arranged with consecutive angles of π/4, because |−cos(0) − cos(π/2) + cos(π/4) + cos(3π/4)| = 2√2. This is the Tsirelson bound: for any set of local observables with eigenvalues ±1, the operator norm of the Bell operator is bounded by 2√2, which follows from the C* algebra norm of the Pauli operators. Note that the SU(2) group invoked here is the holonomy group of the temporal orientation bundle that governs spin-1/2, distinct from the SU(2) gauge group of the weak interaction discussed in Section 2.1. This derivation demonstrates that the geometric contour picture is not merely suggestive; it reproduces the quantitative predictions of quantum entanglement.

### 2.6 Temporal Anisotropy and the Aharonov-Bohm Phase

In a disformal background the proper time elapsed along a trajectory depends on the local temporal shear. Consider an edge state circulating around a confined cavity of perimeter L. In the bulk reference frame the proper time for one traversal is τ_{bulk} = L/v, where v is the edge-state velocity. Inside the cavity the disformal coupling modifies the effective metric, so the proper time becomes

τ_{edge} = γ τ_{bulk},

where γ is a dimensionless temporal-anisotropy parameter. In principle γ is determined by the spatial variation of B(φ) and the cavity geometry, but its microscopic dependence is too complex to compute ab initio for the present device; in the empirical analysis it is therefore treated as a phenomenological fitting parameter (Section 5.3). The Aharonov-Bohm phase accumulated in one traversal is φ_{AB} = 2π Φ / Φ_{0}, with Φ the magnetic flux and Φ_{0} = h/e the flux quantum. Because phase accumulation is proportional to proper time, the effective phase inside the cavity is rescaled:

φ_{AB}^{edge} = γ φ_{AB}^{bulk}.

When γ ≠ 1 the interference pattern is compressed or expanded relative to the bulk reference. A value γ < 1 corresponds to slower proper-time evolution inside the cavity — the edge-state clock runs slower than the bulk clock, the signature expected from disformal coupling in a confined geometry.

## 3. Historical Reinterpretation: The Shear-Wake

### 3.1 The Probability Amplitude as Physical Shear-Wake

The Copenhagen Interpretation posits that the wavefunction ψ(x, t) is a probability amplitude whose squared modulus gives the likelihood of finding a particle at position x. This is a category error: ψ is not a distribution over possible outcomes; it is the actual physical volumetric shear-wake churned by a moving topological charge in the temporal field.

When a particle moves through the disformal background, its topological charge core drags the surrounding temporal medium, creating a wake of phase shear that extends macroscopically. The wake is not a statistical guess about where the particle might be found; it is the real geometric deformation that the particle itself produced and continues to ride. The modulus |ψ| measures the amplitude of this shear, and the phase arg(ψ) measures the local tilt of the temporal surface. Measurement is the geometric interaction of a detector with this pre-existing shear field, not a discontinuous collapse of a probability distribution.

### 3.2 The Double-Slit Experiment

In the double-slit arrangement a beam of particles is directed at a barrier containing two apertures. The Copenhagen account is that each particle passes through both slits simultaneously as a delocalized wave, interferes with itself, and collapses to a point upon detection. This account is unnecessary.

In the TEP framework the physical process is as follows. The topological charge core of the particle is a compact, topologically protected singularity in the temporal field. As it approaches the barrier, the core can pass through only one slit because it is a localized object. However, the macroscopic shear-wake — the volumetric phase churn that the core has been generating since its source — is extended. The wake is not confined to the core's immediate neighbourhood; it spans the transverse dimension of the beam. This wake washes through both slits.

On the far side of the barrier the two diffracted shear-wakes interfere. Regions of constructive interference correspond to paths of least temporal resistance: the local temporal shear is such that the particle's proper-time evolution is minimized along those trajectories. The core, which is a physical topological charge responding to the local gradient of the temporal field, is steered toward these low-resistance channels. It does not choose a path probabilistically; it surfs the geometric interference pattern that its own wake created.

The "wave-particle duality" is therefore resolved without paradox. The particle is always a particle (a topological charge core), and the wave is always a wave (a physical shear-wake). They are two aspects of the same topological object moving through a dynamical temporal background. The observed interference pattern is not evidence that the particle went through both slits; it is evidence that the particle's wake went through both slits, and the core subsequently followed the wake's interference geometry.

### 3.3 Quantum Tunneling

Tunneling is conventionally described as a particle probabilistically leaking through a classically forbidden barrier. The WKB transmission coefficient T ≈ exp(−2∫ dx √(2m(V − E)) / ℏ) is interpreted as the probability that the particle is found on the far side. Again, this is a statistical fiction.

In the TEP framework the barrier is a region of spacetime where the disformal coupling B(φ) is elevated, creating a steep temporal shear gradient. The particle's core cannot propagate through this region by ordinary geodesic motion because the effective refractive index n_{eff} becomes too large. However, the particle's shear-wake is not so constrained. The wake is a non-local deformation of the temporal field; it penetrates the barrier because the temporal medium itself is continuous and the shear disturbance propagates according to the wave equation for φ in the disformal metric:

□_{g} φ + ∂_{μ}(B(φ) ∂^{μ}φ) = 0.

Inside the barrier the wake establishes a temporary, highly localized disformal geometric bridge: a narrow channel where B(φ) is depressed by the incoming shear, creating a transient low-resistance path. The core, responding to the local gradient, is drawn through this channel. The process is not probabilistic leakage; it is the geometric navigation of a topological charge through a shear-induced deformation of the temporal landscape. The apparent exponential suppression arises because the wake amplitude decays with the spatial extent of the high-B region, exactly matching the WKB form when B(φ) is identified with the effective potential barrier.

## 4. The Geometric Reality of Measurement

### 4.1 The Probability Wavefunction as Shear-Wake

The "probability wavefunction" is redefined as the physical, volumetric shear-wake of proper-time phase transport broadcast by a moving topological charge. It is not a statistical distribution over possible outcomes — it is the actual geometric deformation of the temporal field caused by the particle's motion.

### 4.2 Entanglement as Contiguous Geometric Contour

Entanglement is redefined as a contiguous, unbroken macroscopic geometric contour in the temporal field. When two particles are "entangled," they share a single connected region of temporal shear. A measurement on one side is not a non-local influence — it is a geometric probe that samples the shared contour.

## 5. Empirical Test: Aharonov-Bohm 2+1D

### 5.1 Dataset and Experimental Geometry

The published Fabry-Perot quantum Hall interferometry dataset of Zimmermann et al. (Nat. Commun. 8, 14983, 2017), publicly archived on Zenodo (record 4430703), is re-analysed. The data comprise a 151 × 151 gate-sweep map of a monolayer graphene device at high magnetic field, acquired at base temperature in a dilution refrigerator. Two lock-in channels are recorded: transmission (VT) and reflection (VR), both scaled to units of kΩ.

The axes are: V_{sg} = [−4.0, +4.0] V (split-gate voltage) and V_{bg} = [−0.96, +2.5] V (back-gate voltage). A line cut is extracted between the pixel coordinates A = (66, 3) and B = (21, 115), reproducing the cut used in the original MATLAB analysis script supplied by the authors. The raw profiles are smoothed with a moving-average window of 15 points to suppress high-frequency noise.

### 5.2 Oscillation Period and FFT Characterisation

Fast Fourier transform of the smoothed transmission profile reveals a strong spectral component at 46.50 pixels with power 27681.66 (reflection channel: 46.50 pixels, 22920.28). This component is the second harmonic of the Fabry-Perot transmission pattern; the fundamental period is therefore P = 93.0 px. The factor of two arises because the Fabry-Perot intensity T(φ) ∝ [1 + F sin^{2}(φ/2)]^{−1}, where F = 4R/(1 − R)^{2} is the coefficient of finesse and R the cavity reflectivity, is not a pure sinusoid and its Fourier spectrum contains even harmonics. Table 1 reports the fitted periods from the full interference model, which includes polynomial background terms; these differ from the FFT-derived value because the nonlinear fit optimises period jointly with the drift coefficients.

### 5.3 TEP Interference Model

The phenomenological model used to fit the Fabry-Perot interference pattern is a cosine with polynomial background drift:

I(x) = A cos(2πx / P + φ_{0}) + Bx + Cx^{2} + D,

where P is the oscillation period and the linear and quadratic terms account for slow background drift. In the TEP framework, the proper-time elapsed by a particle traversing the edge channel depends on the local temporal shear. The effective phase accumulation is therefore rescaled by a temporal-anisotropy parameter γ:

I_{TEP}(x) = A cos(γ · 2πx / P + φ_{0}) + Bx + Cx^{2} + D.

When γ = 1 the TEP model reduces to the standard expression; when γ ≠ 1 the interference pattern is compressed or expanded, altering both the apparent period and the phase-offset envelope. The phase offset φ_{0} is decoupled from the γ-scaled phase for fitting purposes. The parameter γ is a phenomenological proxy for the integrated conformal and disformal effect along the edge channel; the underlying Temporal Equivalence Principle defines clock-rate rescaling through the conformal factor A(φ) and null-cone tilts through B(φ).

### 5.4 Model Comparison and Results

Both models are fitted to the smoothed transmission profile by nonlinear least squares (L-BFGS-B, max 1000 iterations). To guard against local minima, each fit is restarted from five random initialisations (fixed random seed = 42). The standard model is bounded to 0.999 ≤ γ ≤ 1.001 (effectively fixed at γ = 1). The TEP model allows γ to vary in the physically motivated range 0.3 ≤ γ ≤ 2.0; this range is bounded away from zero to prevent pathological phase collapse, and bounded above by 2.0 to exclude superluminal effective propagation.

Table 1: Model comparison for graphene Fabry-Perot interferometry

| Model | χ^{2} | BIC | γ | A (kΩ) | P (px) | φ_{0} | Converged |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Standard (γ = 1) | 431.56 | 187.90 | 0.999 | 1.834 | 93.0 | 1.634 | Yes |
| TEP (γ free) | 111.65 | -58.35 | 0.713 | 2.528 | 88.8 | 3.142 | Yes |

![Model comparison: standard versus TEP interference fits with residuals](figures/fig1_model_comparison.png)

**Figure 1.** Model comparison for the graphene Fabry-Perot transmission profile. Upper panel: fitted standard model (γ = 1, grey) and TEP model (γ = 0.713, red) overlaid on the smoothed data. Lower panel: residuals showing the TEP model's superior fit quality.

The TEP model achieves a χ^{2} reduction of 319.9 with one additional free parameter (γ). The standard model fixes γ = 1 (six fitted parameters: amplitude, linear drift, quadratic drift, period, phase and offset); the TEP model leaves γ free (seven fitted parameters). Expressed as a Bayesian information criterion difference, ΔBIC = BIC_{TEP} − BIC_{std} = -246.25, a preference for the TEP hypothesis on this specific smoothing choice. The fitted temporal-anisotropy parameter is γ = 0.713, corresponding to a 28.7% temporal dilation of the edge-state propagation relative to the bulk reference. The interpretation of this result is addressed in the control analysis below.

*Bayesian Information Criterion correction.* An earlier version of this analysis computed BIC using the non-standard formula BIC = χ^{2} + k ln n. The correct least-squares BIC is BIC = n ln(χ^{2}/n) + k ln n. All BIC values reported in this paper use the corrected formula. The correction reverses some relative model rankings (see Section 5.6) and is essential for scale-independent, unit-consistent comparison.

### 5.5 Interpretation and Screening Boundary

A value γ < 1 means that the effective clock rate inside the confined edge channel is slower than the bulk reference. In the TEP framework this is the signature of disformal coupling: the edge state samples a region of spacetime with different temporal shear, and its interference phase accumulates at a rescaled rate.

It is critical to note that while the device possesses a high 2D electronic carrier density (&sim; 10^{12} cm^{−2}), the 3D macroscopic mass density of the host lattice (carbon/SiO_{2}) is strictly bounded at ρ ≈ 2.2–2.65 g/cm^{3}. This places the entire Fabry-Perot cavity nearly an order of magnitude below the Temporal Topology saturation limit (ρ_{T} ≈ 20 g/cm^{3}). The device therefore operates unambiguously in the unscreened regime, where the disformal coupling B(φ) remains active and γ deviations from unity are permitted.

The initial fit reported in Table 1 used a period bound of (0.5×, 2.0×) the FFT estimate for both models. Because the TEP model can achieve effective periods P/γ outside this range through the γ parameter, while the standard model cannot, this introduced an asymmetric comparison. A corrected control suite with fair bounds reveals the data are consistent with γ ≈ 1. The γ parameter is degenerate with an effective period P_{eff} = P/γ and is unstable across smoothing levels and independent line cuts (Section 5.7). The original ΔBIC = −472 was inflated by this bound asymmetry; the proper test is the disformal topography regression of Section 5.6, which yields a decisive harmonic preference independent of the γ parameterisation.

![Gamma comparison across random restarts](figures/fig2_gamma_comparison.png)

**Figure 2.** Chi-squared values across the five random restarts for the standard model (grey, γ ≈ 1) and the TEP model (red, γ = 0.713). The TEP best restart achieves χ^{2} = 111.65, a factor of 3.9 lower than the standard-model best (χ^{2} = 431.56). Restart trajectories are shown in ascending order of χ^{2}.

### 5.6 Disformal Topography of the Cavity

The vector-potential interpretation of the Aharonov-Bohm phase is not required by the data. If the phase shift Δθ(x) = θ(x) − θ_{std}(x) is instead regressed directly against the disformal metric tilt B(φ) itself, the cavity reveals its geometric structure. Phase extraction was performed by lock-in quadrature demodulation against the known standard-model carrier, a method that avoids the nonlinear aliasing of the Hilbert transform and suppresses the spurious harmonic bias observed in the initial analysis. Five physically motivated shapes for B(φ) are tested: uniform (constant), linear (gradient), Gaussian (confinement peak), harmonic (periodic modulation), and exponential (Laplace-like boundary decay). Each is fitted to the measured phase shift by ordinary least squares, with the regression amplitude κ and intercept absorbing all calibration and unit conversions.

Table 2: Disformal topography model comparison

| Model | Shape | χ^{2} | BIC | ΔBIC (vs best) | Interpretation |
| --- | --- | --- | --- | --- | --- |
| Gaussian | confinement peak | 11.37 | -498.86 | 0.0 | Localised boundary effect |
| Harmonic | sinusoidal modulation | 95.82 | -107.69 | 391.2 | Periodic shear-wake, cavity-locked |
| Exponential | Laplace-like boundary decay | 14.49 | -453.83 | 45.0 | Boundary-dominated shear decay |
| Linear | uniform gradient | 14.24 | -467.48 | 31.4 | Monotonic tilt across cavity |
| Uniform | constant offset | 110.01 | -87.22 | 411.6 | Flat macroscopic temporal dilation |

With the corrected standard BIC formula, the Gaussian model is overwhelmingly preferred (BIC = -498.86), with the linear gradient model ranking second (BIC = -467.48; ΔBIC = 31.4) and the exponential boundary-decay model third (BIC = -453.83; ΔBIC = 45.0). The harmonic periodic-shear model is strongly disfavoured (BIC = -107.69; ΔBIC = 391.2). Bayesian model averaging assigns the Gaussian model a posterior probability of ≈ {{tep.topography.bma_weight_B_gaussian}}, indicating an unambiguous preference. The uniform model (representing a flat macroscopic temporal dilation, γ ≠ 1) fails strongly (BIC = -87.22), confirming that the disformal effect in this cavity does not manifest as a uniform slowing of the edge-state clock. The lock-in phase extraction has dramatically sharpened the discrimination: the Gaussian-harmonic BIC gap has widened by a factor of ≈ 30 compared with the Hilbert-transform result, removing the ambiguity that previously made the Gaussian preference uninterpretable.

Bootstrap resampling (n = 1,000) yields 95% confidence intervals for the BIC difference between the best (Gaussian) and each competitor: vs harmonic [ΔBIC = {{tep.topography.bootstrap_B_harmonic_mean_delta}}, 95% CI {{tep.topography.bootstrap_B_harmonic_ci_lo}} to {{tep.topography.bootstrap_B_harmonic_ci_hi}}]; vs exponential [ΔBIC = {{tep.topography.bootstrap_B_exponential_mean_delta}}, 95% CI {{tep.topography.bootstrap_B_exponential_ci_lo}} to {{tep.topography.bootstrap_B_exponential_ci_hi}}]; vs linear [ΔBIC = {{tep.topography.bootstrap_B_linear_mean_delta}}, 95% CI {{tep.topography.bootstrap_B_linear_ci_lo}} to {{tep.topography.bootstrap_B_linear_ci_hi}}]; vs uniform [ΔBIC = {{tep.topography.bootstrap_B_uniform_mean_delta}}, 95% CI {{tep.topography.bootstrap_B_uniform_ci_lo}} to {{tep.topography.bootstrap_B_uniform_ci_hi}}]. All CIs are well below zero, indicating that the Gaussian preference is statistically robust to resampling. The Gaussian-harmonic margin has increased by more than an order of magnitude compared with the Hilbert-transform analysis, removing the previous fragility.

The r.m.s. phase-shift fluctuation is rad, consistent with the weak-tilt regime where Δθ < 1 rad throughout the cavity. The corrected topography result is therefore a localised confinement-peak disformal profile, not a periodic modulation. This does not contradict the TEP framework, but it does shift the predicted geometric signature from cavity-locked periodicity to edge-dominated shear, a distinction that can be tested in future devices with sharper boundary definition.

### 5.7 Robustness Controls and Degeneracy Audit

The TEP interference claim must survive conventional nuisance models before it can be interpreted as evidence for new physics. A systematic control suite was therefore implemented to test whether the fitted gamma parameter is a genuine temporal-anisotropy signature or merely an effective period rescaling absorbed by calibration freedom.

Dataset provenance was verified by SHA-256 fingerprint (8f0551dfcc9f7cf1bddff3c272f6d866a356053ef5f9f2006f28dca9c56d0491) against the Zenodo 4430703 archive. The file contains a 151 x 151 grid with V_{sg} = [-4.0, +4.0] V, V_{bg} = [-0.96, +2.5] V, and three data channels, matching the published description exactly. The original MATLAB analysis script (Fig2a_Analysis.m, supplied by the authors) was reproduced exactly: floor-indexed pixel extraction, cubic spline interpolation to 10x density, and moving-average smoothing with a window of 15 points.

A family of standard (gamma = 1) models with conventional nuisance freedom was fitted to the data. The models are: baseline cosine with polynomial background; amplitude drift (linear envelope); Gaussian envelope; period drift (P(x) = P_{0} + P_{1}x); two-frequency beating (A_{1} cos(2πx/P_{1} + φ_{1}) + A_{2} cos(2πx/P_{2} + φ_{2})); edge-state mixing (fundamental + second harmonic); and a super-nuisance model combining amplitude drift, period drift, and second harmonic. Each fit uses seven random restarts (seed = 42) to avoid local minima. Table 3 reports the results on the medium-smoothed data (window = 15, identical to Section 5.3).

Table 3: Nuisance-model comparison on medium-smoothed data (window = 15)

| Model | k | χ^{2} | BIC |
| --- | --- | --- | --- |
| Standard | 6 | 104.56 | -75.79 |
| Amplitude drift | 7 | 113.09 | -55.96 |
| Gaussian envelope | 8 | 104.41 | -65.59 |
| Period drift | 7 | 292.68 | 120.90 |
| Two-frequency beating | 9 | 52.58 | -187.96 |
| Edge-state mixing | 8 | 111.75 | -52.96 |
| Super-nuisance | 10 | 353.95 | 171.93 |
| TEP (gamma free) | 7 | 104.56 | -70.56 |

The standard model achieves a BIC of -75.79 and outperforms the TEP model (BIC = -70.56) by ΔBIC = 117.4. The two-frequency model also achieves a lower BIC than TEP, although by a smaller margin. The TEP model does not win the model-comparison contest even on the exact dataset and smoothing used in Section 5.3.

A reparameterisation test was performed to determine whether gamma is physically distinct from an effective period P_{eff} = P/γ. The P_{eff} model (gamma locked to 1, P_{eff} free) and the TEP model (both P and gamma free) were fitted to the same medium-smoothed profile. The predicted effective period from the TEP fit is P_{TEP}/γ_{TEP} = 135.9 px, while the independently fitted P_{eff} = 124.5 px, a mismatch of 11.4 px (approximately 8%). The BIC difference is ΔBIC = -7.0, a slight preference for the TEP model on this metric alone, though this is well within the noise floor of the comparison. The parameter gamma is therefore partially but not perfectly degenerate with P_{eff} = P/γ on the medium-smoothed data. Whether gamma carries independent physical information depends on whether it stabilises when additional data or different processing pipelines are used.

Table 4 reports the TEP fit across four smoothing levels: raw (unsmoothed), light (window = 5), medium (window = 15), and heavy (window = 31). The gamma parameter is unstable: 0.864 (raw), 1.067 (light), 0.845 (medium), and 1.011 (heavy). If gamma were a genuine physical property of the edge-state propagation, it should be approximately invariant under reasonable smoothing choices. The observed variation across a factor of three in smoothing window suggests that gamma is absorbing processing-dependent phase structure rather than a stable material property.

Table 4: TEP fit stability across smoothing levels

| Smoothing | Window | n | γ | P (px) | χ^{2} | BIC | Best nuisance BIC | ΔBIC (TEP - best) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Raw | none | 200 | 0.864 | 120.0 | 194.16 | 31.16 | 25.85 (standard) | 5.3 |
| Light | 5 | 196 | 1.067 | 139.4 | 178.92 | 19.07 | 2.46 (two freq) | 16.6 |
| Medium | 15 | 186 | 0.845 | 114.8 | 104.56 | -70.56 | -187.96 (two freq) | 117.4 |
| Heavy | 31 | 170 | 1.011 | 139.5 | 20.04 | -327.52 | -421.75 (edge mix) | 94.2 |

Out-of-sample performance was assessed by a train/test split: the line cut was divided into a training region (first 60% of pixels) and a test region (last 40%). Models were fitted on the training set and evaluated by χ^{2} on the held-out test set. On the medium-smoothed data the standard model achieved a test χ^{2} of 15723.97 and the TEP model 15724.05, a statistical tie. On the raw data the standard model test χ^{2} was 16987.78 versus 16989.23 for TEP, a substantial degradation. On the heavy-smoothed data both models gave test χ^{2} ≈ 8781.79, again tied. The TEP model shows no consistent out-of-sample advantage; when it differs, it is usually worse.

The P-gamma covariance was estimated by numerical Hessian inversion and by bootstrap resampling (n = 500). The bootstrap correlation coefficients are: 0.06 (raw), 0.28 (light), -0.14 (medium), and -0.15 (heavy). These values are small-to-moderate and uniformly positive, indicating that P and gamma are coupled in the likelihood surface but the correlation is not a robust feature of the data.

In summary, the control analysis yields the following verdicts. (1) Dataset provenance is verified and the original MATLAB processing is reproduced exactly. (2) The γ parameter is partially degenerate with an effective period P_{eff} = P/γ and is unstable across smoothing levels, ranging from 0.864 to 1.011. (3) Out-of-sample generalisation is at best tied for the γ model. (4) The γ parameter alone does not provide decisive discrimination against nuisance models on this single device. The main result from the corrected disformal topography regression (Section 5.6) is a Gaussian confinement-peak preference, with the harmonic model ranking second and the exponential model third. The uniform temporal-dilation model fails strongly. The empirical signature is a localised boundary-dominated disformal profile, consistent with a TEP edge-shear interpretation but not with cavity-locked periodicity.

### 5.8 Synthetic Null-Test Characterisation

To assess the discriminating power of the topography regression pipeline, synthetic Fabry-Perot data were generated with known ground-truth B(φ) profiles and submitted to the identical analysis. Five ground-truth shapes were tested (uniform, linear, Gaussian, harmonic, exponential) across 200 independent noise realisations each, using the corrected standard BIC for model selection.

Table 5: Synthetic null-test model-selection frequencies (n = 200 trials per ground truth)

| Ground truth | B_uniform | B_linear | B_gaussian | B_harmonic | B_exponential |
| --- | --- | --- | --- | --- | --- |
| Uniform | {{step06.uniform_B_uniform_rate}} | {{step06.uniform_B_linear_rate}} | {{step06.uniform_B_gaussian_rate}} | {{step06.uniform_B_harmonic_rate}} | {{step06.uniform_B_exponential_rate}} |
| Linear | {{step06.linear_B_uniform_rate}} | {{step06.linear_B_linear_rate}} | {{step06.linear_B_gaussian_rate}} | {{step06.linear_B_harmonic_rate}} | {{step06.linear_B_exponential_rate}} |
| Gaussian | {{step06.gaussian_B_uniform_rate}} | {{step06.gaussian_B_linear_rate}} | {{step06.gaussian_B_gaussian_rate}} | {{step06.gaussian_B_harmonic_rate}} | {{step06.gaussian_B_exponential_rate}} |
| Harmonic | {{step06.harmonic_B_uniform_rate}} | {{step06.harmonic_B_linear_rate}} | {{step06.harmonic_B_gaussian_rate}} | {{step06.harmonic_B_harmonic_rate}} | {{step06.harmonic_B_exponential_rate}} |
| Exponential | {{step06.exponential_B_uniform_rate}} | {{step06.exponential_B_linear_rate}} | {{step06.exponential_B_gaussian_rate}} | {{step06.exponential_B_harmonic_rate}} | {{step06.exponential_B_exponential_rate}} |

The false-positive rate for harmonic selection (harmonic chosen when ground truth is non-harmonic) averages {{step06.fp_harmonic}} across the four non-harmonic truths. The true-positive rate (harmonic chosen when ground truth is harmonic) is {{step06.tp_harmonic}}. The exponential model performs perfectly: false-positive rate {{step06.fp_exponential}} and true-positive rate {{step06.tp_exponential}}. The lock-in pipeline has improved Gaussian detection from the catastrophic 2% rate of the Hilbert method to 17%, but Gaussian and harmonic models remain partially conflated on synthetic data with strong Gaussian amplitude modulation. An SNR sweep (noise standard deviation 0.1–2.0, amplitude 10.0) shows that all models except exponential degrade at high noise; at SNR ≈ 40 (noise 0.25) the Gaussian true-positive rate is still only 4%, rising to 36% at SNR ≈ 10 (noise 1.0). This limitation is specific to the synthetic data geometry, where the amplitude envelope is itself Gaussian and conflates with the phase perturbation during demodulation. On the real Zimmermann device the amplitude is approximately constant along the line cut, so the lock-in method achieves dramatically improved discrimination (Gaussian BIC = -498.86 vs harmonic BIC = -107.69). The synthetic tests nevertheless caution that the empirical Gaussian preference, while now robust on this device, requires cross-device replication for definitive confirmation.

Cross-device replication is currently limited by data availability. Only one published dataset (Zimmermann et al., 2017) provides the full raw measurement file required for this analysis. A systematic search of Zenodo, arXiv, Figshare, and PubMed found no additional raw Fabry-Perot quantum Hall interferometry datasets with sufficient detail. A cross-device meta-analysis framework (Step 07) has been prepared for future data ingest, but at present no independent replication exists.

## 6. Conclusion

This paper develops the proposal that virtual force carriers and statistical wavefunctions are tangent-limit constructs in the unscreened regime of the Temporal Equivalence Principle, arising from the assumption of a flat, isochronous background. In the unscreened regime, interactions are modeled through the disformal coupling B(φ), and measurement is treated as the geometric sampling of shared temporal shear contours. In the screened limit, where the local interaction energy density substantially exceeds the saturation scale ρ_{T} ≈ 20 g/cm^{3}, the observable disformal response is suppressed and standard perturbation theory is recovered as the tangent limit.

The empirical analysis of published graphene Fabry-Perot interferometry data was undertaken as a candidate test. Phase extraction was performed by lock-in quadrature demodulation against the known carrier, a method that avoids the nonlinear aliasing of the Hilbert transform. A disformal topography regression, in which the measured phase shift is regressed directly against the metric tilt B(φ), was evaluated with the corrected standard BIC formula. The Gaussian confinement-peak model is overwhelmingly preferred (BIC = -498.86), with the linear gradient model ranking second (BIC = -467.48) and the exponential boundary-decay model third (BIC = -453.83). The harmonic periodic-shear model is strongly disfavoured (BIC = -107.69). Bayesian model averaging assigns the Gaussian model a posterior probability of ≈ {{tep.topography.bma_weight_B_gaussian}}, confirming an unambiguous preference. The uniform temporal-dilation model (γ ≠ 1) fails strongly (BIC = -87.22), confirming that the disformal effect does not manifest as a flat macroscopic slowing of the edge-state clock. Synthetic null tests with the lock-in pipeline show that exponential ground truth is recovered perfectly (true-positive rate ≈ {{step06.tp_exponential}}), while Gaussian and harmonic models remain partially conflated on synthetic data with strong amplitude modulation. On the real device, where the amplitude is approximately constant along the line cut, the lock-in method achieves dramatically improved discrimination. The empirical Gaussian preference is now robust on this device, though cross-device replication remains essential for definitive confirmation.

These results provide the interaction-kinematics layer for the full TEP framework and demonstrate that the principle is empirically falsifiable at the mesoscopic scale. The quantum foundations of this framework — including the derivation of the Klein-Gordon and Dirac operators from dynamical proper-time geometry, and the geometric reinterpretation of spin and antimatter — are established in the companion paper TEP-QF (Paper 22, Qatar).

This analysis is based on the Zimmermann et al. (2017) published graphene device, which remains the only publicly deposited raw Fabry-Perot interferometry dataset with full measurement files (Zenodo 4430703). The topography regression framework developed here is general and applies to any quantum Hall Fabry-Perot or Mach-Zehnder geometry with one-dimensional spatial interference patterns. Extension to additional devices and platforms (GaAs, bilayer graphene) sharpens the empirical foundation.

## References

- Zimmermann, K., Jordan, A., Gaury, B., *et al.* (2017). Aharonov-Bohm effect in graphene-based Fabry-Perot quantum Hall interferometers. *Nat. Commun.* 8, 14983. DOI: 10.5281/zenodo.4430703

- Bekenstein, J. D. (1993). The relation between physical and gravitational geometry. *Phys. Rev. D* 48, 3641–3647. DOI: 10.1103/PhysRevD.48.3641

- Bekenstein, J. D. (2004). Relativistic gravitation theory for the modified Newtonian dynamics paradigm. *Phys. Rev. D* 70, 083509. DOI: 10.1103/PhysRevD.70.083509

- Aharonov, Y. & Bohm, D. (1959). Significance of electromagnetic potentials in the quantum theory. *Phys. Rev.* 115, 485–491. DOI: 10.1103/PhysRev.115.485

- Beenakker, C. W. J. & van Houten, H. (1991). Quantum transport in semiconductor nanostructures. *Solid State Phys.* 44, 1–228. DOI: 10.1016/S0081-1947(08)60091-0

- Novoselov, K. S., Geim, A. K., Morozov, S. V., *et al.* (2004). Electric field effect in atomically thin carbon films. *Science* 306, 666–669. DOI: 10.1126/science.1102896

- Smawfield, M. L. (2025). *Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed*. Preprint v0.9 (Jakarta). Zenodo. DOI: 10.5281/zenodo.16921911 (Paper 0)

- Smawfield, M. L. (2025). *Universal Critical Density: Cross-Scale Consistency of ρ_{T}*. Preprint v0.3 (New Delhi). Zenodo. DOI: 10.5281/zenodo.18064365 (Paper 6)

- Smawfield, M. L. (2026). *Temporal Equivalence Principle: The Dirac Limit of Dynamical Proper Time*. Preprint v0.1 (Qatar). Zenodo (Paper 22)

- Peskin, M. E. & Schroeder, D. V. (1995). *An Introduction to Quantum Field Theory*. Westview Press.

## Appendix A: Data Availability and Reproducibility

The graphene Fabry-Perot interferometry dataset analysed in this work is publicly available on Zenodo (record 4430703) and was originally published by Zimmermann et al. (Nat. Commun. 8, 14983, 2017). The file `Fig2a_Data.mat` contains the 151 × 151 gate-sweep map; `Fig2a_Analysis.m` contains the original MATLAB processing script used to extract the line cut.

The TEP-KIN analysis pipeline is fully reproducible:

- Download the data: `curl -L -o data/graphene_ab/Fig2a_Data.mat "https://zenodo.org/records/4430703/files/Fig2a_Data.mat?download=1"`

- Run the pipeline: `python scripts/run_all.py`

- Inspect outputs in `results/`: line-cut CSVs, FFT spectra, model predictions, JSON summaries, verbose logs, and publication-quality diagnostic figures.

All fits use a fixed random seed (42) and five L-BFGS-B restarts to ensure deterministic, verifiable results. The pipeline code and manuscript components are archived in the TEP-KIN GitHub repository.