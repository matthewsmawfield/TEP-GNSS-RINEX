Title: The Temporal Equivalence Principle: Dynamic Time, Emergent Light Speed, and a Two-Metric Geometry of Measurement

Author: [Redacted for review]

Abstract
This paper proposes a covariant, testable reformulation of relativity in which proper time is a dynamical field and the "speed of light" is an emergent, strictly local invariant rather than a global constant. The framework is built on a single spacetime manifold endowed with two metrics: a gravitational metric gμν and a causal (matter) metric g̃μν to which all non-gravitational fields and clocks couple. The metrics are related by a controlled disformal map, g̃μν = A(φ) gμν + B(φ) ∇μφ ∇νφ, where φ is the time field, A(φ) = exp(2β φ/MPl) is a universal conformal factor, and B(φ) encodes tiny, direction-dependent deformations of the light cone consistent with multi-messenger constraints (|cγ − cg|/c ≲ 10−15 today). Proper time is elevated to a field by postulating that all matter, electromagnetism, and quantum phases evolve with respect to g̃-proper time τ; in local freely falling frames, this guarantees exact local Lorentz invariance and a locally invariant c, while globally it implies that synchronization procedures and one-way light-time measurements become path-dependent in a dynamical-time background.

The full action, field equations, conservation laws, Parametrized Post-Newtonian (PPN) mapping, and screening mechanisms are developed to reconcile terrestrial tests with cosmological dynamics. The breakdown of global simultaneity is formalized using a synchronization-transport law, deriving a convention-independent "synchronization holonomy," an invariant measure of non-integrability of time transport around closed loops. In purely conformal theories this holonomy vanishes once general-relativistic Sagnac and Shapiro terms are removed; a nonzero holonomy at leading order requires non-exact structure provided either by (i) disformal photon coupling B(φ) ≠ 0 or (ii) more general non-metricity. Explicit small-B formulas for the holonomy and the effective photon phase speed are provided, showing how the measured one-way asymmetry is related to φ-gradients and disformal scales under current constraints. The analysis demonstrates that Einstein's assumption of a universal c was a brilliant local theorem arising from the Temporal Equivalence Principle; transcending it demands dynamical time: c remains exactly invariant locally, but global, one-way-inferred "c" values differ by path-dependent amounts that experiments can detect or bound.

A decisive, near-term experimental program is outlined: (1) a closed-loop, multi-leg, one-way time-transfer "triangle test" designed to detect synchronization holonomy at the 10^{−19} fractional level (after averaging) and subtracting known GR effects; (2) interplanetary one-way optical time transfer targeting picosecond-level asymmetries over AU baselines; (3) distance-correlation analysis and environment-dependent screening maps with precision clock networks; (4) multi-messenger searches for distance-correlated photon–gravitational-wave arrival differentials consistent with tightly bounded disformal propagation; (5) matter-wave interferometry and torsion-balance tests sensitive to environment-dependent couplings. Cosmologically, the time field can shift the sound horizon and late-time expansion, potentially easing H0 and S8 tensions while preserving Big Bang Nucleosynthesis (BBN) and CMB acoustic physics.

Known weaknesses in the variable-c literature are addressed by supplying a correct, operationally invariant observable (holonomy), clarifying when conformal couplings cannot produce a signal, and providing realistic, constraint-consistent signal forecasts with explicit error budgets and statistical plans (pre-registration, blinding, publicly released code and data). The resulting theory preserves the empirical pillars of relativity (local Lorentz invariance, gravitational-wave causality, PPN bounds) while extending its conceptual foundation: simultaneity is not only relative but generally non-integrable; the speed of light is not a global constant but the local echo of a deeper, dynamical temporal geometry.

1. Introduction: From a Universal Speed to a Universal Principle of Time
Relativity and quantum theory disagree about time. General relativity (GR) makes proper time geometric—dτ2 = −gμν dxμ dxν / c2—dynamical, and observer-dependent. Quantum mechanics (QM) treats time as an external parameter t in the Schrödinger equation, iħ ∂t|ψ⟩ = H|ψ⟩. The clash has haunted quantization programs for a century and manifested operationally as subtle ambiguities in defining one-way light speeds and simultaneity across extended regions. Precision clocks and time-transfer links now measure gravitational redshifts and velocity-dependent dilations at 10−18 fractional levels, while cosmology exhibits persistent H0 and S8 tensions. In this context, Einstein’s postulate of a universal speed c must be sharpened: it is exact for any local freely falling laboratory, but it cannot be a global property if the flow of time itself is dynamical.

The Temporal Equivalence Principle (TEP) is proposed: all non-gravitational dynamics, signals, and quantum phases evolve according to the proper time defined by a single causal metric g̃μν that couples universally to matter. The rate at which proper time accrues is a field. This elevates "when" to the status "where" acquired in 1915: as space was geometrized, now time's rate is dynamized. The unavoidable consequence is that synchronization conventions are no longer globally integrable: even after removing known GR effects (Sagnac, Shapiro), closed-loop time transport does not necessarily return a clock to its initial epoch, and one-way times differ in opposite directions. The measured "speed of light" between distant clocks is revealed as an emergent ratio of distance to accrued proper time, not a fundamental constant comparable across regions without a theory of time's flow. Locally, c remains exactly invariant; globally, dynamical time imprints tiny, path-dependent asymmetries that can be detected or bounded.

***
**Figure 0: Conceptual Flowchart of the Temporal Equivalence Principle Framework.**
A flowchart diagram with four main stages connected by arrows.
*   **Stage 1 (Top): Foundational Axioms.** Box contains: "1. Two-Metric Geometry (gμν, g̃μν)"; "2. Temporal Equivalence Principle (Matter couples to g̃μν)"; "3. Causal Safety (B(φ) is small)"; "4. Screening (V(φ) hides φ in dense regions)".
*   An arrow points down to Stage 2.
*   **Stage 2 (Middle-Left): Theoretical Dynamics.** Box contains: "Action Principle: S[g, φ, ψ]"; "Field Equations for gμν and φ"; "Screened solutions for φ near massive bodies (e.g., Earth)".
*   An arrow points right to Stage 3.
*   **Stage 3 (Middle-Right): Key Observable Prediction.** Box contains: "Path-dependent time transport"; "Non-integrable synchronization"; "Prediction: Non-zero Synchronization Holonomy H ∮ Ωμ dxμ ≠ 0"; "H requires disformal coupling B(φ) ≠ 0"; "H ∝ (B/A)|∇φ|² × Area".
*   An arrow points down from Stage 3 to Stage 4.
*   **Stage 4 (Bottom): Experimental Tests.** Box contains: "1. Triangle Test (Directly measure H)"; "2. Interplanetary Links (Measure one-way asymmetry Ξ)"; "3. Clock Networks (Distance correlations + screening maps)"; "4. GW/EM Delays (Constrain B(φ))"; "5. Cosmology (Test H(z), rs, S8)".

This structure—from axioms to a unique observable to a suite of falsifiable tests—forms the backbone of the proposal. The central, novel prediction is the synchronization holonomy, an invariant signature of a non-integrable temporal geometry.
***

2. Axioms and the Temporal Equivalence Principle
The framework adopts four axioms:

A1. Two-metric structure on a single manifold. Gravity is described by a Lorentzian metric gμν; matter fields, clocks, and rulers couple to a causal (matter) metric g̃μν. The metrics are related by a disformal map
g̃μν = A(φ) gμν + B(φ) ∇μφ ∇νφ,
with a universal conformal factor A(φ) = exp(2β φ/MPl) and a small disformal function B(φ) consistent with multi-messenger constraints today.

A2. Temporal Equivalence Principle (TEP). All non-gravitational processes evolve according to proper time dτ defined by g̃μν: iħ d|ψ⟩/dτ = H0|ψ⟩, nongravitational fields propagate on g̃-null cones, and ideal clocks tick dτ2 = −g̃μν dxμ dxν/c2. In local freely falling frames for g̃μν, physics reduces to special relativity with invariant c.

A3. Causal safety and locality. Today, cg = cγ within current bounds (|cg − cγ|/c ≲ 10−15). This is enforced by requiring B(φ)(∂φ)2 to be tiny at late times and by choosing A(φ) universal, such that conformal transformations preserve null cones for photons and gravitons when B = 0. Hyperbolicity and energy conditions hold within the EFT domain.

A4. Screening and universality. The coupling A(φ) is universal at leading order; loop-induced composition dependence is tightly bounded. A chameleon-like potential V(φ) yields density-dependent effective mass mφ(ρ) to screen fifth forces and preserve PPN constraints in dense environments while leaving cosmology accessible to dynamics.

These axioms encode Einstein's local invariance as a theorem and extend it: c is exactly invariant for every tangent space of g̃μν, but global synchronization and one-way timing depend on the dynamical field φ.

2.1 Relation to Other Modified Gravity Frameworks
The TEP framework can be situated within the broader landscape of modified gravity and scalar-tensor theories. To clarify its unique features, a comparative analysis is provided:

**Table 1: Comparison with Other Scalar-Tensor Frameworks.**
This table provides a concise comparison of the TEP framework with other major theories of modified gravity, highlighting the key fields, matter coupling principles, predictions for the relative speed of gravity and light, and the primary observable signature that distinguishes each theory.

| Theory/Framework | Key Fields | Matter Coupling | c_g = c_γ? | Key Observable Signature |
|------------------|------------|-----------------|------------|--------------------------|
| **TEP (This work)** | gμν, φ | Universal to g̃μν = A gμν + B ∂μφ ∂νφ | No, if B≠0. | Synchronization Holonomy H |
| **Horndeski** | gμν, φ | Minimal to gμν | Yes, after GW170817 constraints | Modified growth, ISW |
| **DHOST** | gμν, φ | Minimal to gμν | Yes, for specific degenerate classes | Modified growth, screening |
| **TeVeS** | gμν, Aμ, φ | To a combined metric | No | MOND phenomenology, lensing |
| **SMEFT** | SM fields | Lorentz-violating operators | Model-dependent | Anisotropic propagation, CPT violation |

The multi-messenger constraint from GW170817, requiring the speed of gravitational waves c_g to be equal to the speed of light c_γ to high precision, has severely constrained this landscape. The model, with a universal conformal coupling A(φ) and a disformal term B(φ), automatically satisfies c_g = c_γ in the conformal limit (B=0). The non-zero B(φ) introduces a small, calculable deviation c_γ ≠ c_g, which will be the target of future multi-messenger constraints.

Conceptually, TEP differs from many scalar-tensor theories by its foundational principle: the universal coupling of matter to a single metric g̃μν. This is philosophically closer to Bekenstein's TeVeS theory, which also introduced a separate matter metric, but the framework is more minimal, using only a scalar, and makes a unique prediction in the form of the synchronization holonomy.

3. Operational Foundations: Measurement, Simultaneity, and One-Way Light
3.1 What is actually measured
No measurement of c uses null proper time along the photon’s worldline; that would be c = dx/dτ with dτ = 0. Instead, a distance is compared with an elapsed proper time on an observer’s clock. For two spatially separated clocks A and B with worldlines γA, γB and a light path λ from A to B, the measured one-way time tAB is the difference in B’s proper time between emission and reception, after a synchronization convention has assigned simultaneity between A and B. Two-way measurements are convention-independent; one-way measurements are not, unless an invariant observable is provided.

3.2 The synchronization problem in a dynamical-time background
Einstein synchronization assumes (i) reciprocity of propagation and (ii) homogeneity of time standards along the path. If the rate dτ/dt varies spatially or directionally, synchronization by light exchange over extended loops becomes path-dependent. The proper formalism invokes the congruence ũμ of observers (clocks) and the null geodesics kμ of g̃μν. Define simultaneity as orthogonality to ũμ (where Frobenius integrability allows it), and define time transport by mapping a proper-time interval at A to B using g̃-null signals. In GR, the non-closure of such transports arises from rotation of the congruence (vorticity) and spacetime curvature (Sagnac/Shapiro); both can be modeled and removed. In a dynamical-time geometry with disformal couplings, there is an additional, tiny non-exact contribution to time transport that cannot be removed by coordinate choices: a synchronization holonomy.

3.3 A convention-independent observable: synchronization holonomy
Define the holonomy H as the total non-closure in proper time accumulated by transporting synchronization once around a closed, oriented loop C embedded in spacetime:
H ≡ ∮C dτprop,
where dτprop is the infinitesimal proper-time transport inferred along each leg by g̃-null exchange. In SR and pure GR (after subtracting modeled Sagnac/Shapiro), H = 0 for loops small enough that curvature terms cancel to the relevant order. In a disformal time geometry H generically does not vanish. Crucially:
- H is an operational scalar: it depends only on measured proper times and the loop.
- H vanishes in the conformal-only limit (B = 0) after GR effects are subtracted; any claim to the contrary confuses convention with physics.
- H reflects the curvature of an effective time-transport connection that includes a disformal, non-exact component.

This reframes “variable c” claims: the invariant diagnostic is not a raw one-way c, but H. A nonzero H means simultaneity is not integrable beyond GR; this is what dynamic time does to measurement.

4. Geometry and Dynamics: Two Metrics and an Action for Time
4.1 The action and field equations
The action is formulated in the Einstein frame (gravity canonical):
S = ∫ d4x √−g [ (MPl2/2) R − (1/2) K(φ) gμν ∂μφ ∂νφ − V(φ) ] + Smatter[ψ, g̃μν, φ],
g̃μν = A(φ) gμν + B(φ) ∇μφ ∇νφ,   A(φ) = exp(2β φ/MPl),   K(φ) > 0.
Variation yields
Gμν = MPl−2 [Tμν(φ) + Tμν(matt)],
Tμν(φ) = ∂μφ ∂νφ − gμν [ (1/2)K(∂φ)2 + V(φ) ],
and a φ equation
□φ − V′(φ) = −α(φ) T + Sdisf,
with α(φ) ≡ d ln A/dφ = 2β/MPl, T ≡ g̃μν T̃μν (matter-frame trace), and Sdisf the disformal source terms, schematically Sdisf = (1/2) B′ T̃μν ∂μφ ∂νφ + ∇μ[(B T̃μν) ∂νφ] + ... . Total stress-energy is conserved: ∇μ [Tμν(φ) + Tμν(matt)] = 0. In the matter frame, ∇̃μ T̃μν = 0, so matter sees standard conservation.

4.2 Proper time and local rates
A clock’s proper time increment is dτ = (−g̃μν dxμ dxν)1/2/c. For slow observers with small ∂0φ compared to spatial gradients and g̃00 ≈ A g00, dτ/dt ≈ A(φ)1/2 = exp(β φ/MPl). This rescaling does not alter local Lorentz invariance: in any freely falling lab for g̃, Minkowski physics holds and c is invariant. What varies is the mapping between distant proper-time standards, not the local dynamics.

4.3 Screening with a single potential
To satisfy precision tests while allowing cosmological dynamics, a chameleon-like potential V(φ) is adopted whose effective mass meff(ρ) grows with ambient density ρ through the effective potential Veff(φ; ρ) = V(φ) + [A(φ) − 1] ρ. For V(φ) = Λ4 [1 + (Λ/φ)n] (n > 0),
φmin(ρ) ≈ [ n Λn+4 MPl / (2β ρ) ]1/(n+1), meff2(ρ) ≈ (n+1) n Λn+4 / φn+2min,
ensuring short-range fifth forces near dense bodies and long-range effects in the cosmos. Universal A(φ) preserves leading-order equivalence principle (EP) in the matter frame; loop-induced composition dependence is bounded by Eötvös tests and can be made negligible.

5. Photon and Gravitational-Wave Propagation; Emergent “c”
5.1 Conformal safety and disformal tilt
Conformal transformations preserve null cones: with B = 0, photons (EM) and gravitons (GW) share light cones, so cg = cγ globally. Disformal terms tilt the photon cone in directions transverse to ∇φ. For small B(∂φ)2,
cγ2/c2 ≈ 1 − (B/A) |∇⊥φ|2,
to leading order in a local inertial frame. The bound from GW170817/GRB170817A (|cγ − cg|/c ≲ few × 10−15 at z ≈ 0.01) constrains B(φ0)(∂φ)2 today and yields a lower bound on the disformal scale M suppressing higher-dimensional operators.

5.2 What “variable c” really means
Locally, in any lab, c is invariant because experiments interrogate g̃μν. Globally, what varies are one-way times inferred between distant clocks with different proper-time standards across regions where φ varies. The “speed of light” used operationally is an emergent ratio, not a fundamental constant that can be compared across environments without specifying φ. This is the precise way to transcend Einstein’s postulate without contradicting it: c is locally exact; “global c” is meaningless without dynamic time.

6. Synchronization Transport and Holonomy: A Rigorous Derivation
6.1 Time transport connection for synchronization
Define the observer congruence ũμ (normalized with g̃μν ũμ ũν = −1) and the family of g̃-null geodesics kμ carrying time-transfer signals. Consider a small element of time transport from A to B along a segment of a null curve. The transported proper-time increment obeys
dτB = dτA + Ωμ dxμ,
where Ωμ is an effective one-form (a “time-transport connection”) encoding how proper time is mapped along null transport in the presence of spatial variation of A and disformal anisotropy. At leading order in gradients and small B,
Ωμ = ∂μ(1/2 ln A) + Θμ[B, ∂φ, k̂, ũ],
where the first term is a gradient (conformal rescaling) and Θμ captures non-exact contributions from disformal anisotropy and the kinematics (projection along k̂ and ũ). The synchronization holonomy around a closed loop C is
H = ∮C Ωμ dxμ = ∫Σ Fμν dΣμν,
where Fμν ≡ ∂μΩν − ∂νΩμ is the curvature 2-form and Σ is any surface bounded by C.

6.2 Conformal-only limit
If B = 0 and the congruence is hypersurface-orthogonal (no rotation) with known GR corrections subtracted, Ωμ = ∂μ(1/2 ln A) is exact. Its curvature Fμν vanishes identically (commuting derivatives on a scalar), so H = 0. Therefore, in the purely conformal theory there is no leading-order, disformal-independent holonomy once GR effects are removed. Earlier claims to the contrary conflate parameter-dependent synchronization conventions with invariant observables or average a gradient over loops that include temporal excursions; but Stokes’ theorem with Ω exact gives H = 0 for smooth fields on simply connected domains. This resolves a literature confusion and sets a crisp target: any detected nonzero H beyond modeled GR signals the presence of non-exact structure (disformal couplings or more general non-metricity).

6.3 Disformal contribution and small-B formula
For small but nonzero B, Θμ contributes a non-exact piece whose curvature does not vanish. In a local inertial tetrad with ũμ ≈ (1, 0, 0, 0) and a loop composed of null transports with directions n̂i, the leading holonomy scales as
H ≈ ∮C (B/2A) |∇⊥φ|2 n̂ · dx + higher-order,
equivalently
H ≈ ∫Σ ∇ × [ (B/2A) |∇⊥φ|2 n̂ ] · dΣ,
where ∇⊥φ is the component of ∇φ transverse to the local ray direction n̂ on each leg. The precise coefficient depends on the link geometry and the projection details; the key point is that Fμν ∝ ∂[μ ( (B/A) Πν] ), with Π a quadratic form in ∇φ and local ray/observer directions, is generically nonzero in inhomogeneous φ backgrounds. This produces an invariant, disformal-sourced synchronization holonomy that vanishes with B. The Sagnac (vorticity) and Shapiro (curvature) terms can be subtracted using standard GR modeling; their dependence on loop orientation, rotation rate, and gravitational potential is distinct from the disformal term.

6.4 One-way asymmetry and its relation to Ωμ
For two points A and B connected by a straight baseline in a stationary background, define the one-way asymmetry
ΞAB ≡ (tAB − tBA) / (tAB + tBA).
To first order in small fields and after removing GR effects, ΞAB ≈ (1/L) ∫A→B Ω · dl. In the conformal-only limit, ΞAB = 0 beyond GR; with disformal coupling,
ΞAB ≈ (1/2L) ∫A→B (B/A) |∇⊥φ|2 dl,
targeting 10−19–10−15 for plausible gradients and couplings under current bounds. On AU baselines, this maps to tens of femtoseconds to sub-nanoseconds, not microseconds, correcting past overestimates.

7. PPN, EP, and Local Tests with Screening
7.1 PPN mapping
In unscreened regimes with A(φ) = exp(2β φ/MPl) and small α0 ≡ 2β/MPl, the Eddington PPN parameter is γ − 1 ≈ −2 α02/(1 + α02) ≈ −2 α02 for α02 ≪ 1. The Cassini bound |γ − 1| < 2.3 × 10−5 implies α0 ≲ 3 × 10−3. Screening near the Sun and Earth suppresses the effective coupling further via thin-shell factors ΔR/R. Other PPN parameters remain GR-like at leading order. Disformal corrections to PPN enter at higher post-Newtonian orders and are strongly bounded by cg ≈ cγ; with B suppressed today, they are negligible in the solar system.

7.2 Equivalence principle
With universal A(φ), EP holds in the matter frame: all non-gravitational fields see the same metric g̃μν. Composition dependence arises only through small loop corrections to dimensionless constants (e.g., α, mass ratios) and is constrained by Eötvös experiments to η ≲ 10−13 currently and by MICROSCOPE to η ≲ 10−14. The parameter space with α0 ≲ 10−3 and screening ensures compliance. Torsion-balance and atom-interferometer experiments probe environment-dependent variations in αeff and meff and are an important complementary channel.

7.3 Gravitational waves and causality
Gravitational waves propagate on gμν null cones; photons propagate on g̃μν null cones. In the conformal limit these cones coincide globally. With small disformal couplings today, cone tilts are constrained such that |cγ − cg|/c ≲ 10−15 along typical lines of sight; the analysis works in this safe regime. The combined system is hyperbolic and causal within the EFT domain; no superluminal signaling arises in the matter frame.

8. Quantum Dynamics with a Time Field
8.1 Proper-time Schrödinger evolution
TEP stipulates iħ d|ψ⟩/dτ = H0|ψ⟩. Using coordinate time t with dτ = A(φ)1/2 dt (for small ∂0φ),
iħ d|ψ⟩/dt = Heff |ψ⟩,   Heff = A(φ)1/2 H0.
Co-located systems experience a common rescaling and remain synchronized; differences arise only in comparisons across φ-gradients or time variations. Quantum interferometers accrue phases Φ ∝ ∫ A(φ)1/2 E dt along each arm; differential phases probe ∇φ and the disformal geometry when the paths sample different environments. This unifies gravitational redshift and quantum phases with the same time-field coupling, providing a direct interferometric handle on the theory.

8.2 Decoherence and backreaction (outlook)
A full treatment of quantum backreaction on φ is beyond the present effective-field scope but can be incorporated as an additional source term in Sdisf suppressed by loop factors. Decoherence in φ-gradients is expected to be negligible for laboratory interferometers given screening. This is an important avenue for future SMEFT-based mapping.

***
**Inset: What "no-variable-c" tests do—and don't—probe (TEP view)**

**Observable (gauge/sync-invariant).** Define the synchronization holonomy (after GR subtraction), where σ is the time-transport 1-form whose line integral equals calibrated proper-time increment along a leg:

H_resid = ∮_∂Σ (σ̃ - σ_GR) = ∬_Σ (F̃ - F_GR),

the loop non-closure of time transport built from measured proper times (in time units; right-hand rule on Σ). It is gauge/synchronization-invariant and **vanishes in SR/GR and in the conformal-only limit** of TEP. GR subtraction includes: Sagnac, Lense–Thirring/gravito-magnetic, Shapiro, gravitational redshift, with ITRF ephemerides and TT/TDB standards.

How flagship constraints map to H_resid:
- **GW170817 (GW–EM coincidence).** |c_γ - c_g|/c ≲ 10^-15 constrains global cone splits. In TEP, late-time conformal coupling preserves null cones, so EM and GW share causal structure; small disformal tilts today are allowed. This is a boundary condition, not a loop-holonomy test.
- **Cassini (PPN-γ).** Two-way Doppler/Shapiro is reciprocity-even; it calibrates σ_GR to subtract but does not bound H_resid.
- **Resonator MM/KT tests.** Cavities bound closed-path, even-parity (two-way sums) anisotropy at 10^-17–10^-18, yet are blind to odd-parity (direction-reversing one-way differences) non-reciprocity and loop non-closure—the ingredients of H_resid.
- **"GPS works."** Network self-consistency uses explicit GR+Sagnac modeling and largely two-way/common-view calibration. This verifies internal consistency under assumed GR model; not a direction-reversing one-way loop-closure null.
- **Clock redshift & pairwise A↔B tests.** Exquisitely confirm GR locally; **only closed loops** (A→B→C→A with direction reversal) can reveal non-integrability captured by H_resid.

Why classics can be null while H_resid ≠ 0:
1. Conformal null-cone invariance ⇒ no large GW–EM kinematic delays (consistent with GW170817).
2. **∂_t φ = 0 over loop timescale; gradients conservative: no first-order one-way anisotropy.** Along a fixed path, forward/back times cancel at O(B,∇φ); leading effects are second order or require time dependence/kinematics. Thus two-way/closed-path nulls can hold while a loop-holonomy test remains sensitive.

**Experimental falsifier (primary endpoints).** Run a closed-loop, one-way time-transfer (and/or portable-clock) test and report:
- Leg-wise antisymmetry: Δt_AB = t_AB - t_BA (and optionally Ξ_AB ≡ (t_AB - t_BA)/(t_AB + t_BA)).
- Loop holonomy H_resid after subtracting GR (Sagnac/gravito-magnetic/Shapiro). Use triangle/quadrilateral geometries with direction reversal; extend with interplanetary one-way links and multi-species clock networks.
***

9. Cosmology: Background, Perturbations, and Tensions
9.1 Background evolution with coupling
In flat FLRW with cosmic time t, the background equations read
3 MPl2 H2 = ρm + ρr + ρφ,   ρφ = (1/2) φ̇2 + V(φ),
φ̇̇ + 3H φ̇ + V′(φ) = −α(φ) ρm,
ρ̇m + 3H ρm = + α(φ) φ̇ ρm,   ρ̇r + 4H ρr = 0.
The matter coupling induces mild energy exchange. During radiation domination (T ≈ 0), φ is frozen; BBN and early CMB physics remain intact. In matter domination, φ thaws, altering H(z) and potentially driving late-time acceleration with wφ ≈ −1 today for suitable V(φ).

9.2 Recombination timing and H0
Recombination physics depends on atomic microphysics (unchanged locally) and the expansion history H(z). The coupling modifies H(z) without requiring early dark energy spikes. Percent-level shifts in rs and DA can move Planck-inferred H0 toward local values. Unlike ad hoc models, the changes arise from the same A(φ) that governs time transport and metrology. A public Boltzmann-code implementation (CLASS/CAMB branch) with (β, Λ, n) will quantify posteriors and Bayes factors.

9.3 Growth, S8, and scale dependence
Linear growth is modified by background H(z) and, in unscreened regions, by a tiny enhancement of effective Newton’s constant, Geff = G (1 + 2 α2/(1 + mφ2 a2/k2)). Screening inside halos returns growth to GR; in voids, growth is modestly altered. The net effect can reduce S8 inferred from lensing relative to CMB, easing the S8 tension. Predictions include environment-dependent peculiar velocities and void statistics shifts testable with DESI, Euclid, and Rubin.

9.4 Early-universe constraints and disformal safety
At recombination, any disformal effects must be negligible; this is ensured if B(φrec)(∂φ)2 is tiny. The late-time constraints (GW170817) suggest B flows to zero by z ≲ 1; it can be even smaller at z ~ 1100. CMB spectral distortions and damping tails further constrain disformal terms. Disformal couplings are assumed cosmologically irrelevant at early times; if required, they can emerge transiently at high z while preserving bounds.

10. Experimental Program: Decisive Tests and Realistic Targets
10.1 Triangle synchronization holonomy (flagship)
Objective: Measure H around a closed loop formed by three stations (A, B, C) via one-way optical time transfer, after subtracting GR and environmental effects, targeting σ(H)/Tloop ≤ 10−19.

Design:
- Geometry: Two ground stations separated by 1000–3000 km at different altitudes, plus a medium-Earth-orbit satellite to complete the triangle. Large-area loops maximize the signal while keeping modeling under control.
- Links: Stabilized optical two-way time and frequency transfer (TWTTFT) on fibers where possible; free-space optical links for ground-satellite legs with adaptive optics and turbulence models; redundant microwave links as cross-checks.
- Protocol: (i) Two-way synchronization along each edge to set the convention; (ii) one-way time-stamped pulses around the loop in both directions; (iii) continuous environmental monitoring (troposphere, ionosphere, temperature, pressure) and local gravimetry; (iv) careful subtraction of Sagnac and Shapiro delays computed from precise ephemerides and Earth rotation.
- Analysis: Pre-registered Bayesian model including templates for disformal holonomy H[B, ∇φgeometry], GR systematics, and environmental residuals; blind injection tests; open release of raw and calibrated data and code.

Forecast:
- Signal: For parameters saturating |cγ − cg|/c ≲ 10−15 along the most favorable leg and plausible φ-gradients near Earth, H/Tloop ~ 10−16–10−18. This requires 10−18-level stability over hours-to-days, achieved by state-of-the-art optical clocks and stabilized links.
- Systematics: Dominant are residual Sagnac (modeled far below target), tropospheric/ionospheric delay fluctuations (mitigated by dual-frequency, water vapor radiometry, and co-located GNSS), and fiber dispersion (calibrated bidirectionally). The target total systematic budget is below 3 × 10−19 fractional.
- Falsification: Null results at <10−19 across multiple geometries and epochs would constrain (B/A)|∇φ|2 below cosmologically interesting levels, significantly shrinking viable parameter space.

10.2 Interplanetary one-way light-time asymmetry
Objective: Detect ΞAB at 10−16–10−14 on 1–5 AU baselines (0.05–5 ps timing asymmetries).

Design:
- Two drag-free spacecraft with optical lattice clocks (Allan deviation ≤ 10−18 at 10^4–10^5 s) trade time-stamped optical pulses over months.
- Synchronization by slow transport of a clock or by common-view with a third reference; closed-loop calibration legs to isolate disformal contributions.
- Full relativistic modeling of propagation (Shapiro), spacecraft motion, and plasma delays; redundancy via multiple baselines and geometries.

Forecast:
- Signal: ΞAB at 10−15 corresponds to ~0.5 ps at 1 AU; consistent with GW170817 bounds and plausible φ-gradients in the heliosphere.
- Systematics: Plasma dispersion, spacecraft clock noise, drag-free performance, and pointing jitter are dominant; with current/planned technology, picosecond-level transfer is challenging but credible on decadal horizons.

10.3 Clock Network Correlation Analysis and Environmental Screening Maps
Objective: Detect spatial correlations and environmental screening signatures in atomic clock frequency residuals consistent with screened scalar field coupling A(φ) to transition frequencies.

Design:
Phase I - Distance Correlation Analysis:
- Analyze existing precision clock networks (GNSS, optical clock arrays) for distance-dependent correlations in frequency residuals
- Apply phase-coherent cross-spectral analysis between station pairs
- Bin pairs by 3D distance, fit exponential correlation model: C(r) = A·exp(-r/λ) + C₀
- Cross-validate across independent analysis centers to control systematics

Phase II - Environmental Screening Maps:  
- Deploy identical optical clocks at sea-level, mountain, stratospheric balloon, and LEO; intercompare with two-way optical links
- Subtract GR redshift and Doppler shifts; correlate residuals with detailed geophysical models and gravimetry to isolate screening signatures

Forecast:
- Distance correlations: Exponential decay with characteristic length λ ~ 2,000-3,000 km for viable screening parameters
- Altitude dependence: 10⁻¹⁹–10⁻¹⁸ frequency shifts over tens of kilometers for λscr ~ 10 km near Earth
- Multi-center cross-validation expected to show <5% variation in fitted parameters

10.4 Matter-wave interferometry
Objective: Measure ∇∇ ln A(φ) via differential phases in atom interferometers.

Design:
- Tall fountain or long-baseline interferometers (Sr, Rb) with gradiometric configurations and environmental control.
- Combine with mechanical gravimeters to separate GR potential from temporal-field gradients.

Forecast:
- Sensitivity to ∇∇ ln A comparable to clock-array gradients; complementary systematics and geometry.

10.5 Fifth-force and EP tests with engineered screening
Objective: Search for environment-dependent force modulation near screening transitions.

Design:
- Torsion-balance and atom-interferometry experiments in vacuum chambers designed for thin-shells (screened) versus hollow/void (unscreened) regimes; vary chamber pressure and mass geometry.
- Look for modulation of effective acceleration at fractional 10−14–10−13.

Forecast:
- Nulls further constrain β and meff(ρ); positives cross-check with clock and interferometer channels.

10.6 Multi-messenger photon–GW delays
Objective: Detect or bound distance-correlated propagation delays between EM and GW signals.

Design:
- Analyze an ensemble of binary neutron star mergers with prompt EM counterparts; model source emission delays and marginalize over them; test for a trend ∝ distance.

Forecast:
- Disformal models consistent with today’s bounds predict sub-second to millisecond delays. Absence of any distance trend across N ≳ 30 events at <100 ms level rules out disformal parameters of cosmological consequence.

11. Statistical Integrity, Error Budgets, and Open Science
11.1 Pre-registration and blinding
For each decisive test the analysis plan will be registered (likelihoods, priors, detection thresholds, systematics models) before data acquisition, blind the analysis where possible (e.g., time-reversal tags, concealed loop orientation), and commit to publishing nulls and methods.

11.2 Error budgets (example: triangle test)
- Sagnac (Earth rotation): modeled to ≤ 10−21 fractional with IERS Earth orientation parameters.
- Shapiro (GR potentials): modeled to ≤ 3 × 10−20 using precise ephemerides and atmospheric models.
- Troposphere/ionosphere: dual-frequency and water vapor radiometry reduce to ≤ 5 × 10−20; co-located GNSS aids correction.
- Fiber dispersion and delay variations: bidirectional calibration yields ≤ 5 × 10−20 over hours; thermal control critical.
- Link noise and clock instability: present optical clocks achieve ≤ 3 × 10−19 at 10^4 s; averaging over days reaches ≤ 10−19.
- Total target budget: ≤ 1 × 10−19 fractional; classified by additive and correlated components in a hierarchical Bayesian model.

11.3 Power analyses
- Triangle test: To achieve a 3σ detection of H/Tloop = 3 × 10−19 with a single loop time of Tloop = 0.5 s, and white noise σ0 = 1 × 10−18 per loop, N ≈ (3 σ0/H)^2 ≈ (3×10^−18 / 3×10^−19)^2 ≈ 100 loops suffice; with colored noise, longer campaigns (weeks) and overlapping Allan deviations are required; simulations will calibrate effective N.
- GW–EM: With intrinsic-lag scatter σint ≈ 1 s and target disformal trend δt ≈ 0.2 s at 40 Mpc, N ≈ (3σint/δt)^2 ≈ 225 events would give 3σ without lag modeling; with lag priors and population modeling, N ≈ 30–50 may suffice.

11.4 Open code and data
All code (CLASS/CAMB branches implementing the scalar sector; triangle-test simulation toolkits; analysis pipelines) and data (raw and calibrated) will be released under permissive licenses with DOIs. Synthetic injections and posterior predictive checks will be provided for reproducibility.

12. Addressing Weaknesses and Clarifying Claims
- Holonomy and exactness. The flaw in treating aμ = ∂μ ln A as a source of holonomy is corrected. In conformal-only models, Ω is exact and H = 0 after GR subtraction. A nonzero holonomy requires non-exact structure: disformal couplings (B ≠ 0) or more general non-metricity. The derivations and proposed experiments explicitly target this.
- Realistic amplitudes. Microsecond-level claims for AU baselines are corrected: with δcγ/c ≲ 10−15 and D ≈ 1 AU, Δtprop ≲ 0.5 ns; the forecasts are picosecond-scale, consistent with bounds.
- Cosmology. Strong claims about H0/S8 are deferred to delivered Boltzmann-code fits and MCMC posteriors. Qualitative fits motivate the program; quantitative evidence will govern conclusions.
- Error budgets and UQ. Explicit systematics and power analyses are provided with commitment to pre-registered, blinded, open analyses.

13. Mathematical Consistency: Causality, Hyperbolicity, and EFT
- Causality. For small B(∂φ)2, g̃μν is Lorentzian and shares causal cones with gμν up to tiny tilts; microcausality in the matter frame holds. Effective phase velocities refer to different metrics; no paradoxical signaling ensues.
- Hyperbolicity and well-posedness. The coupled (g, φ) system with canonical kinetic terms and chameleon V(φ) admits a well-posed Cauchy problem. Screening introduces no higher-derivative ghosts.
- Energy conditions. Null energy can be preserved for φ; strong energy can be violated as in dark energy, without pathology. Matter-frame positivity is maintained.
- EFT validity. Disformal operators are higher-dimensional, suppressed by scale M; GW170817 bounds require M high enough for EFT control below that scale.

14. Beyond Conformal: Early-Time Disformal Phases and Non-Metricity (Outlook)
While B ≈ 0 is set today, early-time disformal phases could exist without violating constraints if they decay by z ≲ 1. Alternatively, a Weyl-integrable geometry or non-metricity could realize non-exact time-transport structure with different signatures. The present focus is the minimal, disformal-sourced holonomy consistent with data; broader geometric extensions are promising future directions.

15. Philosophical Consequences: Constants, Causality, and the Nature of Time
- Constants reinterpreted. c is a local invariant of g̃μν, not a number to be compared globally without transport rules for time standards. Dimensionless constants (α, mass ratios) remain central; the framework respects tight bounds on their variation by setting Z(φ) ≈ 1 today.
- Simultaneity refined. Einstein taught us that simultaneity is relative; dynamic time teaches that simultaneity is generally non-integrable beyond GR. The invariant content is holonomy: synchronization fails to close around loops in a way that no coordinate choice can eliminate.
- Time as a field. Space was made geometric in 1915; here the rate of time is made dynamical. This neither undermines special relativity nor general relativity; it completes them where global assumptions about uniform time fail.

16. Conclusions: Toward the Next Principle
This paper constructs a complete, covariant theory in which time's rate is a field and the observed constancy of c is a local theorem. A single scalar φ, coupled via a universal conformal factor and a tiny disformal term, suffices to (i) unify gravitational redshift, quantum phase evolution, and cosmological acceleration under the Temporal Equivalence Principle; (ii) preserve PPN and multi-messenger constraints; (iii) predict a sharp, invariant synchronization holonomy that vanishes in SR/GR and conformal-only models but arises from disformal anisotropy; and (iv) motivate a suite of decisive, near-term experiments. The experimental path outlined is concrete, falsifiable, and technologically credible. A positive detection of holonomy or one-way asymmetry at predicted levels would mark a conceptual transition on par with Einstein's: from a universe where geometry is dynamic to one where the flow of time itself is a field, and "c" is its local echo. A null result at 10−19–10−20 fractions would place world-leading bounds on dynamical-time extensions, strongly informing cosmology and metrology alike. Either outcome advances the foundations of physics.

Acknowledgments
The author thanks colleagues in precision time metrology (JILA, NIST, PTB, SYRTE, NMIJ, NPL), gravitational-wave astronomy (LIGO, Virgo, KAGRA), cosmology (Planck, DESI, DES, KiDS, Euclid, Rubin), and gravitational physics for discussions and publicly available tools that informed this work. The author is grateful to reviewers who pressed for correction of earlier mis-scalings and replacement of ambiguous one-way speed claims with invariant observables. All errors remain the author's.

Data and Code Availability
Upon acceptance, the following will be released:
- A CLASS/CAMB branch implementing the scalar sector with conformal/disformal couplings.
- A synchronization-holonomy simulation toolkit including environmental models and GR corrections.
- Pre-registered analysis plans and all data from pilot experiments.
All under permissive licenses with DOIs.

Competing Interests
The author declares no competing interests.

References (selected)
- Einstein, A. (1905). Zur Elektrodynamik bewegter Körper. Ann. Phys. 17, 891–921.
- Einstein, A. (1916). Die Grundlage der allgemeinen Relativitätstheorie. Ann. Phys. 49, 769–822.
- Will, C. M. (2014). The Confrontation between General Relativity and Experiment. Living Rev. Relativity 17, 4.
- Bekenstein, J. D. (1993). The relation between physical and gravitational geometry. Phys. Rev. D 48, 3641.
- Damour, T. & Polyakov, A. M. (1994). The string dilaton and a least coupling principle. Nucl. Phys. B 423, 532.
- Khoury, J. & Weltman, A. (2004). Chameleon fields: Awaiting surprises for tests of gravity in space. Phys. Rev. Lett. 93, 171104.
- Hinterbichler, K. & Khoury, J. (2010). Symmetron fields. Phys. Rev. Lett. 104, 231301.
- Abbott, B. P. et al. (LIGO/Virgo) (2017). GW170817: Observation of gravitational waves from a binary neutron star inspiral. Phys. Rev. Lett. 119, 161101.
- Abbott, B. P. et al. (2017). Multi-messenger observations of a binary neutron star merger. ApJ 848, L12.
- Planck Collaboration (2020). Planck 2018 results. VI. Cosmological parameters. A&A 641, A6.
- Riess, A. G. et al. (2022). A comprehensive measurement of the local value of the Hubble constant. ApJ 934, L7.
- Bothwell, T. et al. (2022). JILA Sr optical lattice clock with 10−18 stability and accuracy. Nature 602, 420–424.
- Touboul, P. et al. (2022). MICROSCOPE mission: test of the equivalence principle in space. Phys. Rev. Lett. 129, 121102.
- Burrage, C. & Sakstein, J. (2018). Tests of chameleon gravity. Living Rev. Relativity 21, 1.
- Bettoni, D. & Liberati, S. (2013). Disformal invariance of second-order scalar-tensor theories. Phys. Rev. D 88, 084020.
- Koivisto, T. S. & Zumalacárregui, T. (2013). Disformal gravity. Phys. Rev. D 88, 084016.
- Ashby, N. (2003). Relativity in the Global Positioning System. Living Rev. Relativity 6, 1.
- Uzan, J.-P. (2011). Varying constants, gravitation and cosmology. Living Rev. Relativity 14, 2.

Appendix A: Disformal Field Equations (Sketch)
With g̃μν = A gμν + B ∂μφ ∂νφ, the inverse is g̃^μν ≈ A^−1 g^μν − (B/A^2) ∂^μφ ∂^νφ / [1 + (B/A)(∂φ)^2], to leading order. The matter stress-energy Tμν(matt) depends on φ through g̃ and carries exchange terms:
∇μ Tμν(matt) = α T ∇νφ + (1/2) B′ T̃μν ∂μφ ∂νφ + ∇μ (B T̃μν ∂νφ) + ...,
consistent with total conservation once the φ equation is used.

Appendix B: Photon Phase Speed and GW–EM Constraints
In a local inertial frame, the photon dispersion relation is g̃^μν k_μ k_ν = 0. For small disformal deformations and k^μ aligned with n̂,
vph2/c2 = 1 − (B/A) |∇⊥φ|2 + O(B2(∂φ)4).
The GW170817 bound requires (B/A) |∇⊥φ|2 ≲ few × 10−15 along typical lines of sight today, bounding B/M4.

Appendix C: PPN with Screening
In thin-shell regimes, the effective scalar coupling is αeff = α0 ΔR/R (≪ α0), yielding
γ − 1 ≈ −2 α2eff/(1 + α2eff) ≈ −2 α2eff,   β − 1 ≈ O(α3eff).
Cassini constrains αeff ≲ 3 × 10−3 locally; chameleon screening ensures ΔR/R is tiny for Sun and Earth.

Appendix D: Holonomy Derivation (Details)
Define the time-transport one-form Ωμ by demanding that a proper-time increment δτ at A transported along a null segment to B satisfies δτB − δτA = Ωμ dxμ. In the conformal-only case, Ωμ = (1/2) d(ln A), exact; thus Fμν = 0 and H = ∮ Ω = 0. With disformal coupling, Ω acquires a direction-dependent, non-exact correction Θμ = (B/2A) Qμ[∂φ, n̂, ũ], where Qμ is quadratic in ∂φ and projects transverse to kμ. Then Fμν = ∂[μΘν] ≠ 0 in generic φ backgrounds, and H = ∫Σ Fμν dΣμν ≠ 0. GR vorticity and curvature contributions are computed and subtracted using standard geodesy formulas; the residual defines the invariant holonomy signal.

Appendix E: Hyperbolicity and Cauchy Problem
With canonical K(φ) > 0 and small B, the φ and metric equations are second-order hyperbolic PDEs with well-posed Cauchy problems on appropriate function spaces. Screening potentials preserve this property. Disformal terms do not introduce higher derivatives of g beyond those in curvature.

Appendix F: Forecast Magnitudes
- Triangle holonomy: H/Tloop ~ 10−16–10−18 for 10^3–10^4 km loops; detectability with 10−19 stability over days.
- Interplanetary ΞAB: ~10−15–10−14 over 1–5 AU; implies 0.05–5 ps asymmetries.
- Altitude residuals: 10−19–10−18 per 10 km if λscr ~ 10 km.
- GW–EM delays: sub-second to millisecond for nearby events; requires ensemble analyses.

Appendix G: Protocol Details (Triangle Test)
- Site characterization: continuous meteorological, tropospheric, ionospheric, and local gravity monitoring; precise geodesy for baselines.
- Link calibration: bidirectional time-transfer calibration, dispersion characterization, and thermal control.
- Data products: raw timestamps, environmental logs, modeled GR corrections, residuals, and posterior samples.

Appendix H: Glossary of Symbols

| Symbol | Description | Section |
|---|---|---|
| gμν | Gravitational metric tensor | 2 |
| g̃μν | Matter/causal metric tensor | 2 |
| φ | The scalar time field | 2 |
| A(φ) | Conformal coupling factor, exp(2βφ/MPl) | 2 |
| B(φ) | Disformal coupling function | 2 |
| TEP | Temporal Equivalence Principle | 1, 2 |
| τ | Proper time defined by g̃μν | 2 |
| H | Synchronization holonomy (a time discrepancy) | 3, 6 |
| Ωμ | Time-transport connection one-form | 6 |
| Fμν | Curvature of the time-transport connection | 6 |
| Θμ | Disformal (non-exact) part of Ωμ | 6 |
| c_g, c_γ | Speed of gravity, speed of light (photons) | 5 |
| V(φ) | Potential for the scalar field φ | 4 |
| α(φ) | Conformal coupling strength, d(ln A)/dφ | 4 |
| β | Dimensionless conformal coupling parameter | 2 |
| M | Suppression scale for disformal/EFT operators | 4 |
| PPN | Parametrized Post-Newtonian formalism | 7 |
| γ, β_PPN | Eddington PPN parameters | 7 |
| H₀ | Hubble constant today | 9 |
| S₈ | Cosmological parameter for structure growth | 9 |
| r_s | Sound horizon at recombination | 9 |
| Ξ_AB | One-way time asymmetry between A and B | 6.4 |
| ũμ | Four-velocity of an observer/clock | 3.2 |
| kμ | Four-wavevector of a null signal | 3.2 |
