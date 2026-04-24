# Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed
**Matthew Lukin Smawfield**
Version: v0.7 (Jakarta)
First published: 18 August 2025 · Last updated: 24 April 2026
DOI: 10.5281/zenodo.16921911

---

## Abstract

This paper proposes a covariant, testable reformulation of relativity in which proper time is a dynamical field and the "speed of light" is an emergent, strictly local invariant rather than a global constant. The framework is built on a single spacetime manifold endowed with two metrics: a gravitational metric $g_{\mu\nu}$ and a causal (matter) metric $\tilde{g}_{\mu\nu}$ to which all non-gravitational fields and clocks couple. The metrics are related by a controlled disformal map, $\tilde{g}_{\mu\nu} = A(\phi) g_{\mu\nu} + B(\phi) \nabla_\mu\phi \nabla_\nu\phi$, where $\phi$ is the time field, $A(\phi) = \exp(2\beta \phi/M_{\text{Pl}})$ is a universal conformal factor, and $B(\phi)$ encodes tiny, direction-dependent deformations of the light cone consistent with multi-messenger constraints ($|c_\gamma - c_g|/c \lesssim 10^{-15}$ today). Proper time is elevated to a field by postulating that all matter, electromagnetism, and quantum phases evolve with respect to $\tilde{g}$-proper time $\tau$; in local freely falling frames, this guarantees exact local Lorentz invariance and a locally invariant c, while globally it implies that synchronization procedures and one-way light-time measurements become path-dependent in a dynamical-time background. The full action, field equations, conservation laws, Parametrized Post-Newtonian (PPN) mapping, and screening mechanisms are developed to reconcile terrestrial tests with cosmological dynamics. The breakdown of global simultaneity is formalized using a synchronization-transport law, deriving a convention-independent "synchronization holonomy," an invariant measure of non-integrability of time transport around closed loops. In purely conformal theories this holonomy vanishes once general-relativistic Sagnac and Shapiro terms are removed; a nonzero holonomy at leading order requires non-exact structure provided either by (i) disformal photon coupling $B(\phi) \neq 0$ or (ii) more general non-metricity. Explicit small-B expansions are provided formulas for the holonomy and the effective photon phase speed, showing how the measured one-way asymmetry is related to $\phi$-gradients and disformal scales under current constraints. The analysis demonstrates that Einstein's assumption of a universal c was a brilliant local theorem arising from the Temporal Equivalence Principle; transcending it demands dynamical time: c remains exactly invariant locally, but global, one-way-inferred "c" values differ by path-dependent amounts that experiments can detect or bound. A decisive, near-term experimental program is outlined: (1) a closed-loop, multi-leg, one-way time-transfer "triangle test" designed to detect synchronization holonomy at the $10^{-19}$ fractional level (after averaging) and subtracting known GR effects; (2) interplanetary one-way optical time transfer targeting picosecond-level asymmetries over AU baselines; (3) distance-correlation analysis and environment-dependent screening maps with precision clock networks; (4) multi-messenger searches for distance-correlated photon–gravitational-wave arrival differentials consistent with tightly bounded disformal propagation; (5) matter-wave interferometry and torsion-balance tests sensitive to environment-dependent couplings. Cosmologically, the time field can shift the sound horizon and late-time expansion, potentially easing $H_0$ and $S_8$ tensions while preserving Big Bang Nucleosynthesis (BBN) and CMB acoustic physics. Known weaknesses in the variable-c literature are addressed by supplying a correct, operationally invariant observable (holonomy), clarifying when conformal couplings cannot produce a signal, and providing realistic, constraint-consistent signal forecasts with explicit error budgets and statistical plans (pre-registration, blinding, publicly released code and data). The resulting theory preserves the empirical pillars of relativity (local Lorentz invariance, gravitational-wave causality, PPN bounds) while extending its conceptual foundation: simultaneity is not only relative but generally non-integrable; the speed of light is not a global constant but the local echo of a deeper, dynamical temporal geometry.

Long-standing confusions about "variable $c$" are resolved by replacing convention-dependent statements with invariant observables tied to measurement procedures. A synchronization one-form $\tilde{\sigma}$ is defined on spacelike slices of the matter metric; its curl $d\tilde{\sigma}$, after subtraction of general-relativistic (GR) Sagnac/gravito-magnetic terms, yields a residual "temporal holonomy" $H$ that vanishes in GR and becomes nonzero only when time is dynamical in this sense. Two key theorems are proven: (i) conformal matter coupling preserves null cones, so photons and gravitons share the same causal structure at late times; (ii) a static $\phi$-gradient generates no first-order one-way light-time anisotropy, placing effects in the femto-to-picosecond regime over astronomical baselines under current bounds. Disformal tilts ($B \neq 0$) are tightly constrained by GW170817-class multi-messenger observations but can source holonomy at levels within reach of next-generation metrology. A covariant action is presented, field equations, conservation, and invertibility/causality conditions are derived, and a 3+1 decomposition is supplied that makes observables explicit. Screening via a continuous Temporal Topology governed by non-linear superposition of field gradients (Temporal Shear) reconciles precision local tests with cosmological evolution, with mapping to Parametrized Post-Newtonian parameters and to the EFT-of-dark-energy $\alpha$-functions with $c_T = 1$ enforced. Decisive experiments with quantitative error budgets are outlined: (1) a ground–ground–satellite triangle time-transfer experiment targeting holonomy at below $10^{-18}$ fractional after GR subtraction; (2) portable-clock "clock anholonomy" around closed paths at the $10^{-19}$ level over days; (3) multi-species clock networks seeking phase-locked annual modulations at $10^{-19}$–$10^{-17}$; (4) interplanetary one-way optical links at picoseconds over AU; (5) altitude-dependent screening maps with optical clocks and atom interferometers; and (6) ensemble multi-messenger tests. A cosmological pipeline plan for CLASS/HyRec modifications and MCMC inference is provided, with commitment to open data and blinded analyses.

Einstein's postulate of universal $c$ was a brilliant, operationally perfect approximation in regimes where time's flow is effectively uniform. In a universe where the rate of time is dynamical yet locally Lorentzian, "the speed of light" emerges as an invariant in every local lab but ceases to be globally universal. The new invariant content resides in path-dependent synchronization defects and holonomies of time transport, not in naive one-way "speeds." If detected, these invariants would inaugurate a post-Einsteinian era: from dynamic geometry to dynamic time.

## 1. Introduction: From a Universal Speed to a Universal Principle of Time

Relativity and quantum theory disagree about time. General relativity (GR) makes proper time geometric—$d\tau^2 = -g_{\mu\nu} dx^\mu dx^\nu / c^2$—dynamical, and observer-dependent. Quantum mechanics (QM) treats time as an external parameter t in the Schrödinger equation, $i\hbar \partial_t|\psi\rangle = \hat{H}|\psi\rangle$. The clash has haunted quantization programs for a century and manifested operationally as subtle ambiguities in defining one-way light speeds and simultaneity across extended regions. Precision clocks and time-transfer links now measure gravitational redshifts and velocity-dependent dilations at $10^{-18}$ fractional levels, while cosmology exhibits persistent $H_0$ and $S_8$ tensions. In this context, Einstein's postulate of a universal speed c must be sharpened: it is exact for any local freely falling laboratory, but it cannot be a global property if the flow of time itself is dynamical.

The Temporal Equivalence Principle (TEP) is proposed: all non-gravitational dynamics, signals, and quantum phases evolve according to the proper time defined by a single causal metric $\tilde{g}_{\mu\nu}$ that couples universally to matter. The rate at which proper time accrues is a field. This elevates "when" to the status "where" acquired in 1915: as space was geometrized, now time's rate is dynamized. The unavoidable consequence is that synchronization conventions are no longer globally integrable: even after removing known GR effects (Sagnac, Shapiro), closed-loop time transport does not necessarily return a clock to its initial epoch, and one-way times differ in opposite directions. The measured "speed of light" between distant clocks is revealed as an emergent ratio of distance to accrued proper time, not a fundamental constant comparable across regions without a theory of time's flow. Locally, c remains exactly invariant; globally, dynamical time imprints tiny, path-dependent asymmetries that can be detected or bounded.

Figure 1. Temporal Equivalence Principle (TEP) architecture. Matter couples to $\tilde{g}_{\mu\nu}$ via a conformal–disformal map controlled by the scalar time field $\phi$, preserving local Lorentz invariance while enabling new global invariants.

## 2. Axioms and the Temporal Equivalence Principle

The framework adopts four axioms:

A1. Two-metric structure on a single manifold.

Gravity is described by a Lorentzian metric $g_{\mu\nu}$; matter fields, clocks, and rulers couple to a causal (matter) metric $\tilde{g}_{\mu\nu}$. The metrics are related by a disformal map

$$\tilde{g}_{\mu\nu} = A(\phi) g_{\mu\nu} + B(\phi) \nabla_\mu\phi \nabla_\nu\phi,$$

with a universal conformal factor $A(\phi) = \exp(2\beta \phi/M_{\text{Pl}})$ and a small disformal function $B(\phi)$ consistent with multi-messenger constraints today.

A2. Temporal Equivalence Principle (TEP).

All non-gravitational processes evolve according to proper time $d\tau$ defined by $\tilde{g}_{\mu\nu}$: $i\hbar d|\psi\rangle/d\tau = \hat{H}|\psi\rangle$, nongravitational fields propagate on $\tilde{g}$-null cones, and ideal clocks tick $d\tau^2 = -\tilde{g}_{\mu\nu} dx^\mu dx^\nu/c^2$. In local freely falling frames for $\tilde{g}_{\mu\nu}$, physics reduces to special relativity with invariant c.

A3. Causal safety and locality.

Today, $c_g = c_\gamma$ within current bounds ($|c_g - c_\gamma|/c \lesssim 10^{-15}$). This is enforced by requiring $B(\phi)(\partial\phi)^2$ to be tiny at late times and by choosing $A(\phi)$ universal, such that conformal transformations preserve null cones for photons and gravitons when $B = 0$. Hyperbolicity and energy conditions hold within the EFT domain.

A4. Screening and universality.

The coupling $A(\phi)$ is universal at leading order; loop-induced composition dependence is tightly bounded. A chameleon-like potential $V(\phi)$ yields a density-dependent effective mass $m_\phi(\rho)$. Rather than operating via discrete boundary cutoffs, screening manifests as a continuous spatial profile (Temporal Topology) governed by the non-linear superposition of field gradients (hereafter Temporal Shear), suppressing fifth forces in dense environments while leaving cosmology accessible to dynamics.

These axioms encode Einstein's local invariance as a theorem and extend it: c is exactly invariant for every tangent space of $\tilde{g}_{\mu\nu}$, but global synchronization and one-way timing depend on the dynamical field $\phi$.

## 2.1 Relation to Other Modified Gravity Frameworks

The TEP framework can be situated within the broader landscape of modified gravity and scalar-tensor theories. To clarify its unique features, a comparative analysis is provided:

Table 1: Comparison with Other Scalar-Tensor Frameworks

Theory/Framework

Key Fields

Matter Coupling

$c_g = c_γ$?

Key Observable Signature

TEP (This work)

$g_{\mu\nu}$, $\phi$

Universal to $\tilde{g}_{\mu\nu} = A g_{\mu\nu} + B \partial_\mu\phi \partial_\nu\phi$

No, if $B≠0$

Synchronization Holonomy $H$

Horndeski

$g_{\mu\nu}$, $\phi$

Minimal to $g_{\mu\nu}$

Yes, after GW170817 constraints

Modified growth, ISW

DHOST

$g_{\mu\nu}$, $\phi$

Minimal to $g_{\mu\nu}$

Yes, for specific degenerate classes

Modified growth, screening

TeVeS

$g_{\mu\nu}$, $A_\mu$, $\phi$

To a combined metric

No

MOND phenomenology, lensing

SMEFT

SM fields

Lorentz-violating operators

Model-dependent

Anisotropic propagation, CPT violation

The multi-messenger constraint from GW170817, requiring the speed of gravitational waves $c_g$ to be equal to the speed of light $c_γ$ to high precision, has severely constrained this landscape. The model, with a universal conformal coupling $A(\phi)$ and a disformal term $B(\phi)$, automatically satisfies $c_g = c_γ$ in the conformal limit ($B=0$). The non-zero $B(\phi)$ introduces a small, calculable deviation $c_γ ≠ c_g$, which will be the target of future multi-messenger constraints.

Conceptually, TEP differs from many scalar-tensor theories by its foundational principle: the universal coupling of matter to a single metric $\tilde{g}_{\mu\nu}$. This is philosophically closer to Bekenstein's TeVeS theory, which also introduced a separate matter metric, but the framework is more minimal, using only a scalar, and makes a unique prediction in the form of the synchronization holonomy.

## 3. Operational Foundations: Measurement, Simultaneity, and One-Way Light

## 3.1 What is actually measured

No measurement of c uses null proper time along the photon's worldline; that would be $c = dx/d\tau$ with $d\tau = 0$. Instead, a distance is compared with an elapsed proper time on an observer's clock. For two spatially separated clocks A and B with worldlines $\gamma_A$, $\gamma_B$ and a light path $\lambda$ from A to B, the measured one-way time $t_{AB}$ is the difference in B's proper time between emission and reception, after a synchronization convention has assigned simultaneity between A and B. Two-way measurements are convention-independent; one-way measurements are not, unless an invariant observable is provided.

## 3.2 The synchronization problem in a dynamical-time background

Einstein synchronization assumes (i) reciprocity of propagation and (ii) homogeneity of time standards along the path. If the rate $d\tau/dt$ varies spatially or directionally, synchronization by light exchange over extended loops becomes path-dependent. The proper formalism invokes the congruence $\tilde{u}_\mu$ of observers (clocks) and the null geodesics $k^\mu$ of $\tilde{g}_{\mu\nu}$. Define simultaneity as orthogonality to $\tilde{u}_\mu$ (where Frobenius integrability allows it), and define time transport by mapping a proper-time interval at A to B using $\tilde{g}$-null signals. In GR, the non-closure of such transports arises from rotation of the congruence (vorticity) and spacetime curvature (Sagnac/Shapiro); both can be modeled and removed. In a dynamical-time geometry with disformal couplings, there is an additional, tiny non-exact contribution to time transport that cannot be removed by coordinate choices: a synchronization holonomy.

## 3.3 A convention-independent observable: synchronization holonomy

Define the holonomy H as the total non-closure in proper time accumulated by transporting synchronization once around a closed, oriented loop C embedded in spacetime:

$$H \equiv \oint_C d\tau_{\text{prop}},$$

where $d\tau_{\text{prop}}$ is the infinitesimal proper-time transport inferred along each leg by $\tilde{g}$-null exchange. In SR and pure GR (after subtracting modeled Sagnac/Shapiro), H = 0 for loops small enough that curvature terms cancel to the relevant order. In a disformal time geometry H generically does not vanish. Crucially:

H is an operational scalar: it depends only on measured proper times and the loop.

H vanishes in the conformal-only limit ($B = 0$) after GR effects are subtracted; any claim to the contrary confuses convention with physics.

H reflects the curvature of an effective time-transport connection that includes a disformal, non-exact component.

This reframes "variable c" claims: the invariant diagnostic is not a raw one-way c, but H. A nonzero H means simultaneity is not integrable beyond GR; this is what dynamic time does to measurement.

Figure 2. Synchronization holonomy on a spacelike slice in the matter frame. The residual invariant $H_{\text{resid}}$ subtracts GR Sagnac and gravitomagnetic contributions, isolating dynamical-time effects.

## 4. Disformal Invertibility, Causality, and Hyperbolicity

## Invertibility and signature

The inverse of $\tilde{g}_{\mu\nu}$ exists for $A>0$ and $1 + (B/A^2)(\partial\phi)^2 \neq 0$:

$$\tilde{g}^{\mu\nu} = A^{-2} \left[ g^{\mu\nu} - \frac{(B/A^2) \partial^\mu\phi \partial^\nu\phi}{1 + (B/A^2)(\partial\phi)^2} \right].$$

For Lorentzian signature, require $A>0$ and $B(\partial\phi)^2 > -A^2$. $B$ is assumed small and gradients bounded in all regimes of interest.

## Causality

The matter cone is inside or equal to the gravitational cone when $B \geq 0$ and gradients are modest; no closed causal curves arise for small $B$. With $B\to0$ at late times, gravitational and matter null cones coincide ($c_T = c_{\text{EM}}$). With $B$ small, any phase differences in propagation are minute and bounded by multi-messenger results.

## Hyperbolicity

The scalar's canonical kinetic ensures hyperbolic evolution with well-posed Cauchy problem where $V''(\phi) > 0$ near minima. Maxwell's equations in a Lorentzian matter metric remain hyperbolic. Linearized gravity propagates on the Einstein-frame cone. No ghosts or gradient instabilities arise for canonical $\phi$ and small $B$. The EFT is valid below the disformal scale $M$, with higher-dimensional operators suppressed.

## 5. Local Lorentz Invariance, Proper Time, and the Emergence of c

## Proper time

Clocks measure $d\tau^2 = -\tilde{g}_{\mu\nu} dx^\mu dx^\nu/c^2$. In a local lab with small velocities and $\partial_0\phi \ll |\nabla\phi|$, $-\tilde{g}_{00} \approx A(\phi)^2 (-g_{00}) - B(\partial_0\phi)^2$, so

$$\frac{d\tau}{dt} \approx A(\phi)^{1/2},$$

up to $O(B)$ corrections. This "dynamical time law" rescales all frequency standards locally.

## Local c

All small labs measure an invariant $c$; the conformal factor rescales both clocks and rulers uniformly, preserving null cones. Thus the empirical fact "$c$ is constant" remains a theorem of the theory at the local level.

## Global measures

Global "speeds" involve synchronized endpoints. Because $d\tau$ depends on $\phi$'s history, the global assignment of simultaneity acquires non-integrability when $\phi$ varies in space or time. That non-integrability—not a local violation of Lorentz invariance—is the source of observable departures.

## 6. Synchronization in a Dynamical Time: One-Way Light and Holonomy

## Operational definitions

Two-way light speed is synchronization-independent and has established $c$'s local invariance. One-way measures require synchronized clocks at A and B; Einstein synchronization assumes time-orthogonal slices and propagation symmetry. In a dynamical $\phi$ background, slow clock transport and one-way synchronization are path and history dependent.

## Key theorems

### Theorem 1 (Conformal null-cone invariance)

For $\tilde{g}_{\mu\nu} = A(\phi)^2 g_{\mu\nu}$, null vectors of $g_{\mu\nu}$ are null for $\tilde{g}_{\mu\nu}$. Maxwell's action is conformally invariant in 4D, so photon trajectories are null with respect to both metrics. Gravitational and electromagnetic waves share null cones when $B = 0$ at late times.

### Theorem 2 (No first-order one-way anisotropy in static φ)

For static $\phi$ with $\partial_0\phi = 0$, the one-way light-time difference for forward/backward propagation along the same path cancels to first order in $\nabla\phi$. Any residual is $O((\nabla\phi)^2)$ or due to time dependence/kinematics. Over astronomical baselines and with current bounds on $\alpha \equiv d \ln A/d\phi$, this places effects in the femto-to-picosecond regime, removing claims of microsecond-scale anomalies.

Proof sketch.

Parameterize the path coordinate $s \in [0,L]$; write $t_\to = \int ds A(\phi_0 + s \partial_\parallel\phi)/c$ and $t_\leftarrow$ similarly along the reverse. Linear terms in $\partial_\parallel\phi$ cancel exactly; see Appendix A for full derivation.

## Synchronization one-form and holonomy

Decompose $\tilde{g}_{\mu\nu}$ in 3+1 form:

$$\tilde{g}_{\mu\nu} dx^\mu dx^\nu = -\tilde{N}^2 dt^2 + \tilde{h}_{ij} (dx^i + \tilde{N}^i dt)(dx^j + \tilde{N}^j dt).$$

Define the synchronization one-form $\tilde{\sigma} = -\tilde{N}^{-1} \tilde{N}_i dx^i$. For a closed spatial loop $C$, the synchronization holonomy is

$$H = \oint_C \tilde{\sigma}.$$

In GR with stationary spacetimes, $H$ reproduces the Sagnac/gravito-magnetic effect; these known contributions are subtracted using geodesy and ephemerides. The residual holonomy beyond GR is

$$H_{\text{resid}} = \iint_\Sigma d\tilde{\sigma} - \text{(GR Sagnac/gravito-magnetic)},$$

with $\Sigma$ spanning $C$. This residual vanishes in SR/GR with $B=0$ and stationary $A$; it is generically nonzero with disformal corrections ($B \neq 0$) and/or temporal variation of $A$ combined with motion through $\nabla\phi$.

## Gauge and protocol invariance

$\tilde{\sigma}$ is a one-form on the spatial slice induced by the matter metric and a chosen time function $t$. Changing the synchronization protocol shifts $\tilde{\sigma} \to \tilde{\sigma} + d\chi$ for some scalar $\chi$ when the change corresponds to a redefinition of time labels without physical alterations; the loop integral of an exact form vanishes: $\oint d\chi = 0$. Therefore $H_{\text{resid}}$, constructed from measured proper times and two-way calibrations and expressed via $d\tilde{\sigma}$ after GR subtraction, is invariant under admissible synchronization changes. This addresses the critique that holonomy is conventional: the nonzero part is the curvature (curl) of $\tilde{\sigma}$, not its exact part.

## Disformal corrections to σ̃

For small $B$ and modest $\partial\phi$,

$$\tilde{N} \approx A N \left[1 - \frac{B}{2A^2 N^2} (n \cdot \partial\phi)^2\right], \quad \tilde{N}_i \approx A \left[N_i + \frac{B}{A^2} (\partial_i\phi)(n \cdot \partial\phi)\right],$$

with $n^\mu$ the unit normal to slices. Then

$$\delta\tilde{\sigma} \approx -\frac{1}{N} \left[\delta N_i - \frac{N_i}{N} \delta N\right] dx^i,$$

and

$$\Delta(d\tilde{\sigma}) = d(\delta\tilde{\sigma}),$$

with $\Delta$ indicating the residual beyond GR. Time dependence of $\phi$ and motion through $\nabla\phi$ generate non-vanishing curl.

## 7. Screening, PPN, Equivalence Principle, and Disformal Bounds

## Screening

A density-dependent effective potential $V_{\text{eff}}(\phi; \rho) = V(\phi) + [A(\phi) - 1] \rho$ makes $\phi$ acquire a large effective mass in dense environments. Rather than invoking discrete thin-shell boundaries, screening operates via the continuous spatial profile of the chameleon field (Temporal Topology). The high ambient density in deep potential wells suppresses the local field gradient (Temporal Shear), ensuring short-range fifth-force suppression while leaving the field light cosmologically. This reconciles local null tests with cosmological dynamics.

## PPN mapping

In unscreened regimes, the PPN parameter $\gamma - 1 \approx -2 \alpha_0^2/(1 + \alpha_0^2) \approx -2 \alpha_0^2$ with $\alpha_0 = \alpha(\phi_\infty)$; Cassini's $|\gamma - 1| < 2.3\times10^{-5}$ implies $\alpha_0 \lesssim 3\times10^{-3}$. Near massive bodies, the suppression of Temporal Shear (vanishing field gradient) reduces the effective coupling to $\alpha_{\text{eff}} \ll \alpha_0$, cleanly preserving PPN bounds without invoking rigid thin-shell approximations. Other PPN parameters remain GR-like for small $\alpha_0$.

## Equivalence principle

TEP ensures universality in the matter frame. In the Einstein frame, the apparent fifth force is bounded by Eötvös experiments; MICROSCOPE gives $\eta \lesssim 10^{-14}$. Sectoral dilaton-like couplings (to $\alpha$, $\mu$, quark masses) are constrained to $|d| \lesssim 10^{-5}$–$10^{-6}$ by composition tests and clock ratios.

## Disformal constraints

GW170817/GRB170817A constrain $|c_\gamma - c_g|/c \lesssim \text{few}\times10^{-15}$ at $z \approx 0.01$. This requires $B(\phi_0) (\partial\phi)^2 \ll 1$ today along astrophysical paths. $B$ is adopted as effectively negligible at late times on cosmological paths; small residual $B$ remains a controlled source for holonomy without violating bounds.

## 8. Cosmology: Background Evolution, Growth, and EFT Mapping

## Background

In spatially flat FLRW (Einstein frame),

$$3 M_{\text{Pl}}^2 H^2 = \rho_m + \rho_r + \rho_\phi,$$

$$\ddot{\phi} + 3H \dot{\phi} + V'(\phi) = -\alpha(\phi) \rho_m + O(B),$$

$$\dot{\rho}_m + 3H \rho_m = +\alpha(\phi) \dot{\phi} \rho_m,$$

$$\dot{\rho}_r + 4H \rho_r = 0.$$

In the matter frame, the continuity equation for $\tilde{\rho}_m$ is standard; the Einstein-frame exchange encodes energy flow between $\phi$ and matter. During radiation domination, $T \approx 0$ and $\phi$ is overdamped; it thaws during matter domination and can act as dark energy.

## EFT of dark energy mapping

The theory maps to EFT-of-DE with $\alpha_T = 0$ enforced to satisfy $c_T = 1$. The braiding $\alpha_B$ and Planck-mass run $\alpha_M$ are both small and controlled by $\alpha(\phi)$ and $\dot{\phi}$. Bounds from large-scale structure and ISW require $|\alpha_M|$, $|\alpha_B| \ll 1$ at late times. The parameter choices respect these.

## Recombination and H₀

Early dark energy fraction $\Omega_\phi(z \approx 1100) \lesssim \text{few percent}$ and small $\alpha(\phi_{\text{rec}})$ keep recombination microphysics intact (HyRec constraints) and shift the sound horizon $r_s$ by $\lesssim1$–2%, enough to ease $H_0$ tension. Sectoral couplings $d_e$ at recombination must be tiny to preserve $\alpha$, $m_e$; negligible $d_e$ is adopted at early times.

## Growth and S₈

In unscreened low-density regions, the effective Newton constant is $G_{\text{eff}}(k,a) = G [1 + 2 \alpha_0^2/(1 + m_\phi^2 a^2/k^2)]$, giving mild scale-dependent growth, raising $f\sigma_8$ on large scales while screening suppresses changes in high-density regions. A net reduction of $S_8$ inferred from weak lensing by ~0.01–0.03 is possible while preserving CMB lensing. Detailed MCMC fits with CLASS/HyRec modifications are planned.

## Standard sirens

In the late-time conformal limit, EM and GW share null cones; standard siren distances are unaffected kinematically. Subtle, detector time-standard effects may introduce an integrated $\phi$-history contribution at the sub-second level over ~100 Mpc; ensemble analyses can bound or detect this.

## 9. Quantum Clocks, Interferometry, and Composition Dependence

## Quantum evolution

Proper-time Schrödinger evolution $i\hbar d|\psi\rangle/d\tau = \hat{H}|\psi\rangle$ becomes $i\hbar d|\psi\rangle/d\tilde{t} = A(\phi) \hat{H} |\psi\rangle$ in the matter frame. Transition frequencies scale as $\nu \propto A(\phi)$ with tiny sectoral corrections from dilaton-like couplings $d_e$, $d_\mu$, $d_q$:

$$\delta \ln \nu = \delta \ln A + K_\alpha \delta \ln \alpha + K_\mu \delta \ln \mu + K_q \delta \ln X_q + \ldots$$

## Species sensitivity

Different clock transitions carry different $K$-coefficients; by comparing ratios, one can disentangle the universal $A(\phi)$ scaling from composition dependence and constrain $d_e$, $d_\mu$, $d_q$.

## Interferometry

Phases acquire contributions proportional to the integral of $A(\phi)$ along arms. Atom interferometers with vertical baselines can sense $\partial_h\phi$ at $10^{-21}$–$10^{-22}$ m$^{-1}$; photon interferometers are sensitive at different bands.

## Decoherence

Rapid $\phi$ variations would induce dephasing at rates $O(\partial_t\phi)$, strongly bounded by clock stabilities. Adiabatic evolution in the lab is assumed, consistent with null drift bounds.

## Inset: What "no-variable-c" tests do—and don't—probe (TEP view)

Observable (gauge/sync-invariant).

Define the synchronization holonomy (after GR subtraction), where $\sigma$ is the time-transport 1-form whose line integral equals calibrated proper-time increment along a leg:

$$H_{\rm resid}\;=\;\oint_{\partial\Sigma}\!(\tilde\sigma-\sigma_{\rm GR})

\;=\;\iint_{\Sigma}\!(\tilde F-F_{\rm GR}),$$

the loop non-closure of time transport built from measured proper times (in time units; right-hand rule on $\Sigma$). It is gauge/synchronization-invariant and

vanishes in SR/GR and in the conformal-only limit

of TEP. GR subtraction includes: Sagnac, Lense–Thirring/gravito-magnetic, Shapiro, gravitational redshift, with ITRF ephemerides and TT/TDB standards.

How flagship constraints map to $H_{\rm resid}$:

GW170817 (GW–EM coincidence).

$|c_\gamma-c_g|/c\!\lesssim\!10^{-15}$ constrains global cone splits. In TEP, late-time conformal coupling preserves null cones, so EM and GW share causal structure; small disformal tilts today are allowed. This is a boundary condition, not a loop-holonomy test.

Cassini (PPN-$\gamma$).

Two-way Doppler/Shapiro is reciprocity-even; it calibrates $\sigma_{\rm GR}$ to subtract but does not bound $H_{\rm resid}$.

Resonator MM/KT tests.

Cavities bound closed-path, even-parity (two-way sums) anisotropy at $10^{-17}\!-\!10^{-18}$, yet are blind to odd-parity (direction-reversing one-way differences) non-reciprocity and loop non-closure—the ingredients of $H_{\rm resid}$.

"GPS works."

Network self-consistency uses explicit GR+Sagnac modeling and largely two-way/common-view calibration. This verifies internal consistency under assumed GR model; not a direction-reversing one-way loop-closure null.

Clock redshift & pairwise A↔B tests.

Exquisitely confirm GR locally;

only closed loops

(A→B→C→A with direction reversal) can reveal non-integrability captured by $H_{\rm resid}$.

Why classics can be null while $H_{\rm resid}\neq0$:

Conformal null-cone invariance ⇒ no large GW–EM kinematic delays (consistent with GW170817).

$\partial_t\phi = 0$ over loop timescale; gradients conservative: no first-order one-way anisotropy.

Along a fixed path, forward/back times cancel at $\mathcal O(B,\nabla\phi)$; leading effects are second order or require time dependence/kinematics. Thus two-way/closed-path nulls can hold while a loop-holonomy test remains sensitive.

Experimental falsifier (primary endpoints).

Run a closed-loop, one-way time-transfer (and/or portable-clock) test and report:

Leg-wise antisymmetry: $\Delta t_{AB}=t_{AB}-t_{BA}$ (and optionally $\Xi_{AB}\!\equiv\!(t_{AB}-t_{BA})/(t_{AB}+t_{BA})$).

Loop holonomy $H_{\rm resid}$ after subtracting GR (Sagnac/gravito-magnetic/Shapiro).

Use triangle/quadrilateral geometries with direction reversal; extend with interplanetary one-way links and multi-species clock networks.

## 10. Experimental Proposals and Falsifiability

A suite of decisive, cross-checking experiments is proposed that collectively falsify or confirm the theory. All observables are dimensionless ratios or calibrated residuals; all designs include nulls, blinding, and open data.

Figure 3. Sensitivity summary. Left: predicted triangle holonomy amplitude vs baseline with instrument noise floor. Right: portable‑clock anholonomy vs loop duration. Shaded bar indicates illustrative GW170817 disformal constraint region.

## A. Triangle synchronization holonomy (ground–ground–satellite)

Geometry.

Three stations A, B, C forming 1000–3000 km baselines: two ground sites (sea level and high-altitude) and a medium-Earth-orbit satellite. Optical two-way time transfer (TWTT) on each edge provides calibration; stabilized fibers and free-space optical links carry both calibration and one-way signals.

Protocol.

Establish Einstein synchronization on each edge via TWTT. Execute one-way transfers around the loop in both senses at high cadence for months. Record raw timestamps, environmental monitors (pressure, temperature, humidity), refractivity profiles, TEC for ionosphere, and precise ephemerides.

Modeling.

Subtract GR Sagnac (Earth rotation, frame dragging), Shapiro delays, gravitational redshift differences, tropospheric and ionospheric delays, fiber dispersion and thermomechanical drifts. Use GNSS and gravimetric models.

Observable.

$H_{\text{resid}} = \oint \tilde{\sigma} - \text{(GR model)}$. Forecast magnitude $O(10^{-18}$–$10^{-16})$ fraction per loop time (0.1–1 s), limited by disformal and time-varying $A(\phi)$ effects. Null in GR by design.

Error budget (fractional per loop, targets after months):

- Clock instability after averaging: 5×10⁻¹⁹

- Two-way calibration residual: 2×10⁻¹⁹

- Troposphere residual: 2×10⁻¹⁹

- Fiber path noise after stabilization: 5×10⁻¹⁹

- Ephemeris/geodesy: 2×10⁻¹⁹

- GR subtraction residual: 5×10⁻¹⁹

-

Total systematic floor (rss):

~8×10⁻¹⁹

Target: below 10⁻¹⁸; demonstrated capability reaches 10⁻¹⁹ regime

Falsification.

Null at <1×10⁻¹⁹ fractional across seasons/geometry excludes disformal/time-varying $A$ signatures with late-time cosmological relevance.

## B. Portable-clock "clock anholonomy"

Design.

Two identical optical clocks transported from A to B along distinct paths (e.g., sea-level highway vs. mountain pass), durations ~1–3 days, then compared at B to a stationary clock. Common-view time transfer provides epoch.

Prediction.

Path-dependent discrepancy $\Delta_{12}$ at few$\times10^{-19}$ for plausible $\alpha \dot{\phi}_0$ and $\partial_h\phi$ under screening; null at $10^{-20}$ bounds $\partial_t\phi$ and $\partial_h\phi$ tightly.

Systematics.

Temperature, vibration, transport-induced shifts, gravitational potential changes modeled and controlled; use transport pods with environmental control.

## C. Multi-species clock network: annual modulations

Network.

Global network of optical clocks (Sr, Yb, Al⁺, Hg⁺, Ca, H(1s–2s)) cross-compared over years.

Prediction.

Phase-locked annual modulations in differential ratios $\nu_A/\nu_B$ with amplitudes $10^{-19}$–$10^{-17}$, phases tied to orbital eccentricity. Species-amplitude ratios reflect $\alpha_A - \alpha_B$; a fit yields $\alpha(\phi)$ and dilaton coefficients $d_e$, $d_\mu$, $d_q$.

Nulls.

Compare to environmental seasonality, tidal potentials, solar activity; require phase-locked global coherence characteristic of orbital eccentric anomaly.

## D. Interplanetary one-way optical time transfer

Design.

Two drag-free spacecraft with $10^{-18}$-class optical lattice clocks, separated by 1–5 AU. Optical comb-based one-way time transfer, with third node for calibration (Earth or a relay). Kinematic synchronization via slow-clock transport or common-view transponders.

Prediction.

Geometry-dependent one-way asymmetry parameter $\Xi_{AB} \equiv (t_{AB} - t_{BA})/(t_{AB} + t_{BA})$ at $10^{-15}$–$10^{-14}$ (0.05–5 ps) for disformal tilts consistent with GW bounds; purely conformal temporal-variation effects are much smaller and geometry-averaging helps isolate them.

Systematics.

Plasma delays, pointing jitter, thermal drifts, deep-space clock performance; anticipate >decade timeline.

## E. Clock Network Correlation Analysis and Environmental Screening Maps

Objective.

Detect spatial correlations and environmental screening signatures in atomic clock frequency residuals consistent with screened scalar field coupling $A(\phi)$ to transition frequencies.

Design.

Phase I - Distance Correlation Analysis:

• Analyze existing precision clock networks (GNSS, optical clock arrays) for distance-dependent correlations in frequency residuals

• Apply phase-coherent cross-spectral analysis between station pairs

• Bin pairs by 3D distance, fit exponential correlation model: $C(r) = A \cdot \exp(-r/\lambda) + C_0$

• Cross-validate across independent analysis centers to control systematics

Phase II - Environmental Screening Maps:

• Deploy identical optical clocks at sea-level, mountain, stratospheric balloon, and LEO; intercompare with two-way optical links

• Subtract GR redshift and Doppler shifts; correlate residuals with detailed geophysical models and gravimetry to isolate screening signatures

Forecast.

• Distance correlations: Exponential decay with characteristic length $\lambda \sim 2,000-3,000$ km for viable screening parameters

• Altitude dependence: $10^{-19}$–$10^{-18}$ frequency shifts over tens of kilometers for $\lambda_{\text{scr}} \sim 10$ km near Earth

• Multi-center cross-validation expected to show <5% variation in fitted parameters

## F. Multi-messenger ensemble

Design.

Stack $N \gtrsim 30$ multi-messenger events with prompt EM counterparts; marginalize astrophysical lag distributions. Seek distance-correlated trends in $t_{\text{EM}} - t_{\text{GW}}$.

Prediction.

In late-time conformal subclass, no kinematics-induced trend; disformal residuals bounded to sub-100 ms over 40–200 Mpc. Nulls constrain $B(\phi_0)(\partial\phi)^2$.

## 11. Statistical and Open Science Principles

Pre-registration.

Publish analysis plans (models, priors, nulls, thresholds) for triangle and portable-clock experiments; record any deviations.

Blinding.

Blind event-time stamps and calibration offsets; employ independent teams for calibration and analysis.

Open data/code.

Release raw timestamps, environmental monitors, calibration logs, and pipelines with DOIs; release CLASS/HyRec modification, atom-sensitivity database, and meta-analytic code.

Hierarchical inference.

Use heavy-tailed likelihoods (Student-t) to mitigate outliers; multi-level random effects by domain; selection models for publication bias; explicit covariance modeling for shared systematics.

## 12. Addressing Critiques and Clarifying Claims

No EM–GW kinematic delay in conformal subclass.

With $B=0$ today, EM and GW share null cones; any observed delays are astrophysical/source or detector-time-standard effects. Earlier overstatements are corrected and aligned with GW170817 constraints.

Static $\phi$ gradients do not generate first-order one-way anisotropy.

An explicit derivation (Appendix A) shows cancellation, correcting prior misinterpretations.

Holonomy invariant and not a synchronization artifact.

The observable $H_{\text{resid}}$ is constructed from physical proper-time measurements and is the integral of a curvature two-form $\Delta(d\tilde{\sigma})$ beyond GR; protocol changes shift $\tilde{\sigma}$ by exact differentials, which integrate to zero over closed loops.

Causality and hyperbolicity.

Explicit invertibility and signature conditions are given, restricting to canonical kinetics and small disformal couplings; matter-frame causality is preserved. A full hyperbolicity proof for generic $B(\phi)$ is delegated to a companion technical paper; the small-$B$ regime is safe.

Cosmology claims are hypotheses to be tested.

The outline shows how modest $\phi$ dynamics can ease $H_0/S_8$ while respecting early-universe bounds; quantitative MCMC fits are required and planned.

## 13. Philosophical and Conceptual Implications

Simultaneity beyond Einstein.

Special relativity makes simultaneity observer-dependent; dynamic time makes synchronization non-integrable at the global scale. The invariant content is a holonomy of time transport: moving clocks around closed loops in a dynamical time background returns path-dependent offsets after subtracting GR effects.

Constants clarified.

The speed of light is an invariant in local tangent spaces; globally, "$c$" is not a number but a family of operational ratios dependent on clock histories in a time field. Variation of constants becomes a question about dimensionless ratios across environments and epochs.

Machian undertone.

The rate of time responds weakly to the stress-energy of matter through $\alpha T$, giving a principled, covariant flavor to the idea that the "rest of the universe" influences local rates, without violating local physics.

## 14. Conclusions

This paper articulates a covariant framework in which the rate of time is a dynamical field with universal matter coupling. The architecture preserves local Lorentz invariance and null-cone structure (in the conformal limit), conforms to multi-messenger bounds, and yields new invariant observables—synchronization holonomy and clock anholonomy—that vanish in GR and are measurable with modern clocks. Screening reconciles local nulls with cosmological dynamics; a controlled disformal sector allows minute, bounded photon-cone tilts that provide clean targets for holonomy experiments. The theory is falsifiable with realistic experiments and promises to clarify persistent cosmological tensions.

Einstein moved us from absolute time to relative simultaneity and dynamic geometry. The next step is to recognize that the flow of time itself is dynamical, and that "the speed of light" is the local echo of a deeper temporal geometry. If the predicted holonomies are observed, physics will enter a new epoch in which dynamic time joins dynamic geometry as a foundation. If not, uniquely strong bounds will have been set and the operational bedrock of $c$ and simultaneity clarified to unprecedented precision.

## Appendix A: Proofs and Key Derivations

## A1. Conformal null-cone invariance

Let $\tilde{g}_{\mu\nu} = A(\phi)^2 g_{\mu\nu}$ with $A > 0$. A vector $k^\mu$ null with respect to $g_{\mu\nu}$, $g_{\mu\nu} k^\mu k^\nu = 0$, satisfies $\tilde{g}_{\mu\nu} k^\mu k^\nu = A^2 g_{\mu\nu} k^\mu k^\nu = 0$. Maxwell's action $S = -(1/4) \int \sqrt{-g} F_{\mu\nu} F^{\mu\nu}$ is conformally invariant in 4D: under $\tilde{g}_{\mu\nu} = \Omega^2 g_{\mu\nu}$, $\sqrt{-\tilde{g}} F_{\mu\nu} F^{\mu\nu} = \sqrt{-g} F_{\mu\nu} F^{\mu\nu}$. Hence photon geodesics are conformally invariant; null cones coincide.

## A2. No first-order one-way anisotropy in static φ

Consider a straight path $x(s)$ along $\nabla\phi$ with $s \in [0, L]$. Expand $A(\phi(s)) = 1 + \alpha \phi(s) + O(\alpha^2)$. Forward time:

$$t_\to = \int_0^L ds \frac{A(\phi(s))}{c} = \frac{L}{c} + \frac{\alpha}{c} \int_0^L ds [\phi_0 + s \partial_\parallel\phi] = \frac{L}{c} + \frac{\alpha}{c}(\phi_0 L + \frac{1}{2} L^2 \partial_\parallel\phi).$$

Backward time: parameterize $s' = 0\to L$ with $\phi(s') = \phi_1 - s' \partial_\parallel\phi$, $\phi_1 = \phi_0 + L \partial_\parallel\phi$:

$$t_\leftarrow = \frac{L}{c} + \frac{\alpha}{c}(\phi_1 L - \frac{1}{2} L^2 \partial_\parallel\phi) = \frac{L}{c} + \frac{\alpha}{c}(\phi_0 L + \frac{1}{2} L^2 \partial_\parallel\phi).$$

Thus $t_\to - t_\leftarrow = 0 + O((\partial\phi)^2, \alpha^2)$. The first nonzero term is second order.

## A3. Synchronization holonomy invariance under protocol shifts

Define $\tilde{\sigma} = -\tilde{N}^{-1} \tilde{N}_i dx^i$ on a slice $\Sigma_t$ in the matter frame. Let the synchronization protocol change to $t' = f(t, x^i)$ with $f$ smooth and invertible; the induced change in $\tilde{\sigma}$ is $\tilde{\sigma} \to \tilde{\sigma} + d\chi$, with $\chi$ a well-defined scalar on the slice (the explicit form depends on $f$). For a closed loop $C = \partial\Sigma$, the integral of the exact differential vanishes: $\oint_C d\chi = 0$ (Stokes' theorem). Therefore the total holonomy $H \equiv \oint_C \tilde{\sigma}$ is invariant under such protocol shifts; only the curvature two-form $d\tilde{\sigma}$ carries invariant content. The residual observable—defined by subtracting the known GR contributions—is $H_{\text{resid}} \equiv \iint_\Sigma [d\tilde{\sigma} - d\sigma|_{\text{GR}}]$, which extracts the non-GR dynamical-time content.

## A4. Disformal inverse and causality condition

Given $\tilde{g}_{\mu\nu} = A^2 g_{\mu\nu} + B \partial_\mu\phi \partial_\nu\phi$, write $\partial_\mu\phi$ as a one-form $u_\mu$. The inverse satisfies $\tilde{g}^{\mu\nu} \tilde{g}_{\nu\sigma} = \delta^\mu_\sigma$. Ansätz: $\tilde{g}^{\mu\nu} = A^{-2} (g^{\mu\nu} + C u^\mu u^\nu)$. Solve for $C$:

$$A^{-2} [g^{\mu\nu} + C u^\mu u^\nu] [A^2 g_{\nu\sigma} + B u_\nu u_\sigma] = \delta^\mu_\sigma$$

$$\Rightarrow \delta^\mu_\sigma + A^{-2} (B + A^2 C + B C (u \cdot u)) u^\mu u_\sigma = \delta^\mu_\sigma$$

$$\Rightarrow A^2 C + B + B C (u \cdot u) = 0 \Rightarrow C = - \frac{B}{A^2 + B (u \cdot u)}.$$

Thus $\tilde{g}^{\mu\nu} = A^{-2} \left[ g^{\mu\nu} - \frac{(B/A^2) u^\mu u^\nu}{1 + (B/A^2) (u \cdot u)} \right]$. Lorentzian signature requires $1 + (B/A^2)(u \cdot u) > 0$.

## Appendix B: Photon Phase Speed and GW–EM Constraints

In a local inertial frame, the photon dispersion relation is $\tilde{g}^{\mu\nu} k_\mu k_\nu = 0$. For small disformal deformations and $k^\mu$ aligned with $\hat{n}$,

$$c_{\gamma}^2/c^2 \approx 1 - (B/A) |\nabla_\perp\phi|^2$$

The GW170817 bound requires (B/A) |∇⊥φ|2 ≲ few × 10−15 along typical lines of sight today, bounding B/M4.

## Appendix C: Glossary of Symbols

Symbol

Description

Section

$g_{\mu\nu}$

Gravitational metric tensor

2

$\tilde{g}_{\mu\nu}$

Matter/causal metric tensor

2

$\phi$

The scalar time field

2

$A(\phi)$

Conformal coupling factor, $\exp(2\beta \phi/M_{\text{Pl}})$

2

$B(\phi)$

Disformal coupling function

2

TEP

Temporal Equivalence Principle

1, 2

$\tau$

Proper time defined by $\tilde{g}_{\mu\nu}$

2

$H$

Synchronization holonomy (a time discrepancy)

3, 6

$\Omega_\mu$

Time-transport connection one-form

6

$F_{\mu\nu}$

Curvature of the time-transport connection

6

$\Theta_\mu$

Disformal (non-exact) part of $\Omega_\mu$

6

$c_g$, $c_\gamma$

Speed of gravity, speed of light (photons)

5

$V(\phi)$

Potential for the scalar field $\phi$

4

$\alpha(\phi)$

Conformal coupling strength, $d(\ln A)/d\phi$

4

$\beta$

Dimensionless conformal coupling parameter

2

$M$

Suppression scale for disformal/EFT operators

4

PPN

Parametrized Post-Newtonian formalism

7

$\gamma$, $\beta_{\text{PPN}}$

Eddington PPN parameters

7

$H_0$

Hubble constant today

9

$S_8$

Cosmological parameter for structure growth

9

$r_s$

Sound horizon at recombination

9

$\Xi_{AB}$

One-way time asymmetry between A and B

6.4

$\tilde{u}_\mu$

Four-velocity of an observer/clock

3.2

$k^\mu$

Four-wavevector of a null signal

3.2

## References

Einstein, A. (1905). Zur Elektrodynamik bewegter Körper. Ann. Phys. 17, 891–921.

Einstein, A. (1916). Die Grundlage der allgemeinen Relativitätstheorie. Ann. Phys. 49, 769–822.

Will, C. M. (2014). The Confrontation between General Relativity and Experiment. Living Rev. Relativity 17, 4.

Bertotti, B., Iess, L. & Tortora, P. (2003). A test of general relativity using radio links with the Cassini spacecraft. Nature 425, 374–376.

Bekenstein, J. D. (1993). The relation between physical and gravitational geometry. Phys. Rev. D 48, 3641.

Damour, T. & Polyakov, A. M. (1994). The string dilaton and a least coupling principle. Nucl. Phys. B 423, 532.

Khoury, J. & Weltman, A. (2004). Chameleon fields: Awaiting surprises for tests of gravity in space. Phys. Rev. Lett. 93, 171104.

Hinterbichler, K. & Khoury, J. (2010). Symmetron fields. Phys. Rev. Lett. 104, 231301.

Abbott, B. P. et al. (LIGO/Virgo) (2017). GW170817: Observation of gravitational waves from a binary neutron star inspiral. Phys. Rev. Lett. 119, 161101.

Abbott, B. P. et al. (2017). Multi-messenger observations of a binary neutron star merger. ApJ 848, L12.

Planck Collaboration (2020). Planck 2018 results. VI. Cosmological parameters. A&A 641, A6.

Riess, A. G. et al. (2022). A comprehensive measurement of the local value of the Hubble constant. ApJ 934, L7.

Bothwell, T. et al. (2022). JILA Sr optical lattice clock with 10−18 stability and accuracy. Nature 602, 420–424.

Touboul, P. et al. (2022). MICROSCOPE mission: test of the equivalence principle in space. Phys. Rev. Lett. 129, 121102.

Burrage, C. & Sakstein, J. (2018). Tests of chameleon gravity. Living Rev. Relativity 21, 1.

Bettoni, D. & Liberati, S. (2013). Disformal invariance of second-order scalar-tensor theories. Phys. Rev. D 88, 084020.

Koivisto, T. S. & Zumalacárregui, T. (2013). Disformal gravity. Phys. Rev. D 88, 084016.

Ashby, N. (2003). Relativity in the Global Positioning System. Living Rev. Relativity 6, 1.

Uzan, J.-P. (2011). Varying constants, gravitation and cosmology. Living Rev. Relativity 14, 2.

Michelson, A. A. & Morley, E. W. (1887). On the relative motion of the Earth and the luminiferous ether. Am. J. Sci. 34, 333–345.

Kennedy, R. J. & Thorndike, E. M. (1932). Experimental establishment of the relativity of time. Phys. Rev. 42, 400–418.

Müller, H. et al. (2003). Modern Michelson-Morley experiment using cryogenic optical resonators. Phys. Rev. Lett. 91, 020401.

Herrmann, S. et al. (2009). Test of the isotropy of the speed of light using a continuously rotating optical resonator. Phys. Rev. Lett. 95, 150401.

IERS Conventions (2010). Gérard Petit and Brian Luzum (eds.). IERS Technical Note No. 36, Frankfurt am Main: Verlag des Bundesamts für Kartographie und Geodäsie.

## How to cite

You can cite all versions by using the DOI:

10.5281/zenodo.16921911

BibTeX:

@misc{Smawfield_TEP_2025,

author       = {Matthew Lukin Smawfield},

title        = {Temporal Equivalence Principle: Dynamic Time &

Emergent Light Speed},

year         = {2025},

publisher    = {Zenodo},

doi          = {10.5281/zenodo.16921911},

url          = {https://doi.org/10.5281/zenodo.16921911},

note         = {Preprint}

}

## Contact

For questions, comments, or collaboration opportunities regarding this work, please contact:

Matthew Lukin Smawfield

matthew@mlsmawfield.com

