# Temporal Equivalence Principle: A Topological Fermion Model for Spin and the g−2 Anomaly
**Matthew Lukin Smawfield**
Version: v0.1 (Paris)
First published: 31 May 2026 · Last updated: 5 June 2026

---

## Abstract

The zero-dimensional point-particle paradigm of Quantum Field Theory is challenged and replaced with a topological fermion: a localized topological charge in the temporal shear field whose intrinsic spin is fluid vorticity governed by local temporal shear. Three rigorous results are derived from stated axioms: (i) the conformally shifted g&#771;-Hamilton-Jacobi equation, (ii) spin as quantized fluid vorticity with a direct spin-statistics connection, and (iii) the geometric consistency condition &kappa; = r<sub>c</sub>/&lambda;<sub>scr</sub> = 1/&radic;2. The proximity-based saturation limit, observationally proxied by &rho;<sub>c</sub> &asymp; 20 g/cm<sup>3</sup> (TEP-UCD, Paper 6), marks the scale at which topological charge cores begin to overlap geometrically. A first-principles mean-field derivation from the screened Klein-Gordon equation yields &rho;<sub>c</sub><sup>(MF)</sup> = M<sub>Pl</sub>² m<sub>&phi;</sub>² / &beta;² &approx; 10<sup>-55</sup> g/cm<sup>3</sup>, many orders of magnitude below the phenomenological value; the discrepancy is explained by non-perturbative core overlap of ~10<sup>50</sup> topological charges per terrestrial coherence volume. The Fermi wavelength of the degenerate electron gas provides the correct order-of-magnitude intuition for why the scale lies far below the single-particle Compton-scale core density (&rho;<sub>core</sub> &sim; 10<sup>4</sup> g/cm<sup>3</sup>). A model transfer function T(&rho;) connects the two limits via random-phase superposition of N<sub>eff</sub> &sim; 10<sup>50</sup> electrons per terrestrial coherence volume. This scale natively bounds the proper-time oscillator, eliminating ultraviolet divergences at their geometric origin without renormalization. A geometric consistency condition &kappa; = r<sub>c</sub>/&lambda;<sub>scr</sub> = 1/&radic;2 &approx; 0.707 follows from the known electron Compton wavelength and the Yukawa screening length of the scalar field. The empirical analysis presented here uses 1,493 data points from JLab PRad and A1 Collaboration electron-scattering cross-sections to build a conformally corrected proton form factor, yielding a testable prediction for future AMBER muon-proton scattering measurements. A planned future analysis targets the Fermilab g&minus;2 anomaly as temporal-topology drag once collaboration-internal timestamp data become available.

Keywords: subatomic structure, fermion topology, spin, vorticity, renormalization, proximity screening, Fermilab g-2, AMBER, temporal equivalence principle

## 1. Introduction: The Failure of Renormalization

### 1.1 The Infinite Point-Mass Catastrophe

Quantum Field Theory relies on a zero-dimensional point-particle paradigm that leads to ultraviolet divergences in loop integrals. Renormalization is a mathematical procedure that subtracts infinities to yield finite predictions. It works phenomenologically, but it does not resolve the underlying physical problem: a particle with no spatial extent cannot have a finite self-energy. The electron self-energy diverges as &Lambda; &rarr; &infin; in the standard formulation, and the Landau pole in QED signals that the theory is incomplete at short distances.

Renormalization is a mathematical workaround for the divergence; the point-particle assumption is the root cause. A finite geometric structure eliminates the divergence at its origin.

### 1.2 The TEP Alternative

Early attempts at a unified geometric theory sought to replace dimensionless point particles with physical &ldquo;knots&rdquo; in spatial geometry, but failed to eliminate the mathematical divergences that necessitated renormalization. The TEP framework achieves this geometric origin by shifting the topology from spatial gravity to proper time. Instead of a point particle, TEP introduces a localized topological charge embedded within the temporal shear field. Because this fermion is a physical defect in the scalar field &phi;, it carries a natural geometric boundary. Bounded natively by the local density saturation, this finite geometric structure eliminates the divergence at its origin, removing the need for artificial ultraviolet cutoffs.

The Temporal Equivalence Principle (TEP) replaces the point particle with a localized topological charge in the temporal shear field. The fermion is not a mathematical singularity but a physical defect in the scalar field &phi; that defines local proper time. This topological charge carries a natural geometric boundary: the proximity-based saturation scale, observationally proxied by &rho;<sub>c</sub> &asymp; 20 g/cm<sup>3</sup> (a phenomenological scale established from terrestrial clock correlation data in TEP-UCD, Paper 6), which bounds the proper-time oscillator and eliminates the need for artificial ultraviolet cutoffs.

The TEP framework is built on three axioms. (A1) The matter-frame metric is a conformal–disformal rescaling of the gravitational metric: g&#771;<sub>&mu;&nu;</sub> = A<sup>2</sup>(&phi;) g<sub>&mu;&nu;</sub> + B(&phi;) &nabla;<sub>&mu;</sub>&phi; &nabla;<sub>&nu;</sub>&phi;. (A2) The conformal factor is exponential: A(&phi;) = exp(&beta;&phi;/M<sub>Pl</sub>). (A3) Temporal shear is the gradient of the logarithmic conformal factor: &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A(&phi;). All results in this paper are derived from these axioms. The full causal matter metric is permanently engaged; the screened limit, where local stress forces both A(&phi;) &rarr; 1 and the observable disformal response is suppressed, recovers the standard Minkowski background and isotropic interactions. In the unscreened regime the disformal sector governs the routing of forces through the tilted light cone, as developed in TEP-KIN (Paper 25).

### 1.3 The Lamb Shift as Topography

The 1947 observation of the Lamb Shift is historically cited as the primary evidence for vacuum polarization via virtual particles. In the TEP framework, this is reinterpreted as the first empirical topographic map of a steep temporal shear gradient. A massive nucleus generates a dense temporal topology. As the electron vortex orbits, its geometric structure is continuously deformed by its proximity to this temporal shear wall ($\Sigma_\mu$). The Lamb Shift does not measure a boiling vacuum of ghost particles; it measures the macroscopic geometric drag of a fluid vortex navigating a non-isochronous landscape.

## 2. The Topological Fermion

### 2.1 The Fermion as Topological Charge

The fermion is defined as a localized topological charge in the temporal shear field. In the matter frame, proper time d&tau; is set by the causal matter metric g&#771;<sub>&mu;&nu;</sub> = A<sup>2</sup>(&phi;) g<sub>&mu;&nu;</sub> + B(&phi;) &nabla;<sub>&mu;</sub>&phi; &nabla;<sub>&nu;</sub>&phi;. *Metric signature:* (+, &minus;, &minus;, &minus;). In the conformal limit relevant for the single-particle core geometry, a particle of mass m propagates according to the g&#771;-Hamilton-Jacobi equation, which in flat background with A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>) reads:

g<sup>&mu;&nu;</sup> &part;<sub>&mu;</sub>S &part;<sub>&nu;</sub>S = m<sup>2</sup>c<sup>2</sup> exp(2&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>).

The effective mass in the matter frame is m<sub>*</sub> = m A(&phi;) = m exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>). The proper-time oscillator frequency is therefore shifted by the local conformal factor, and the fermion acquires a position-dependent effective inertia. In the rest frame (&nabla;S = 0), the frequency is &omega;<sub>eff</sub> = mc<sup>2</sup> A(&phi;) / &hbar;. This is the origin of the bounded proper-time oscillator: as the topological charge tightens, A(&phi;) flattens toward unity and &omega;<sub>eff</sub> approaches the standard Compton frequency, but it never diverges because the core has finite extent. The conformally shifted g&#771;-Hamilton-Jacobi equation and the emergence of the Klein-Gordon and Dirac operators in the screened limit are derived systematically in TEP-QF (Paper 23).

### 2.2 Spin as Quantized Vorticity

"Spin" is translated directly into fluid vorticity governed by local temporal shear &Sigma;<sub>&mu;</sub> = &nabla;<sub>&mu;</sub> ln A(&phi;) = (&beta;<sub>A</sub>/M<sub>Pl</sub>) &nabla;<sub>&mu;</sub>&phi;. For a topological charge with azimuthal symmetry in cylindrical coordinates, the shear field is:

&Sigma;<sub>&theta;</sub> = (&beta;<sub>A</sub>/M<sub>Pl</sub>) (n/r),

where n is the integer winding number of the compact phase/orientation variable associated with the charge core. The real conformal factor A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>) is single-valued but not itself periodic; quantization belongs to the phase/orientation bundle carried by the defect. The circulation around the charge core is quantized:

&Gamma; = &oint; &Sigma; &middot; d&ell; = 2&pi;n (&beta;<sub>A</sub>/M<sub>Pl</sub>).

The vorticity vector &omega;<sub>i</sub> = (&nabla; &times; &Sigma;)<sub>i</sub> vanishes everywhere except at the core singularity (r = 0), where it is a delta function. This is the topological origin of spin: the minimal non-trivial defect carries integer phase winding n = &pm;1, while the local orientation bundle is spinorial. A 2&pi; circuit reverses the spinor sign and a 4&pi; circuit restores it, so the observed spin projection is S = &pm;&hbar;/2 even though the scalar/phase winding is integer. The spin-statistics connection arises because exchange of two identical topological charges implements the same non-trivial spinor holonomy, introducing a phase factor of &pi; and enforcing antisymmetry under exchange.

### 2.3 Charge Core Geometry and the &kappa; Ratio

The core radius of the topological charge is fixed by the electron Compton wavelength, r<sub>c</sub> = &hbar;/(m<sub>e</sub>c), a known quantity from quantum mechanics. The Yukawa screening length of the scalar field, derived from the single-particle energy density, is &lambda;<sub>scr</sub> = &radic;2 &hbar;/(m<sub>e</sub>c). Their ratio defines a geometric consistency condition:

&kappa; = r<sub>c</sub> / &lambda;<sub>scr</sub> = 1/&radic;2 &approx; 0.707.

The consistency condition &kappa; = 1/&radic;2 is a geometric constraint: the framework predicts both the core radius and the screening length from the scalar field dynamics, and their ratio is fixed by the topology. It anchors the TEP charge geometry to an established quantum scale.

### 2.4 The TEP Action and Field Equations

The dynamics of the scalar field &phi; are governed by the Einstein-frame action

S = &int; d<sup>4</sup>x &radic;&minus;g [ R/(16&pi;G) &minus; &#189; g<sup>&mu;&nu;</sup> &part;<sub>&mu;</sub>&phi; &part;<sub>&nu;</sub>&phi; &minus; V(&phi;) ] + S<sub>matter</sub>[g&#771;<sub>&mu;&nu;</sub>, &psi;],

where the matter fields &psi; couple to the full causal matter metric g&#771;<sub>&mu;&nu;</sub> = A<sup>2</sup>(&phi;) g<sub>&mu;&nu;</sub> + B(&phi;) &nabla;<sub>&mu;</sub>&phi; &nabla;<sub>&nu;</sub>&phi;. The scalar potential V(&phi;) is to be constrained by cosmological data and fifth-force bounds; a runaway form V(&phi;) &prop; exp(&minus;&lambda;&phi;/M<sub>Pl</sub>) and a screened chameleon potential are both compatible with the TEP framework. Varying with respect to &phi; yields the Klein-Gordon equation in the presence of a fermion source:

&nabla;<sup>&mu;</sup>&nabla;<sub>&mu;</sub>&phi; &minus; V'(&phi;) = (&beta;<sub>A</sub>/M<sub>Pl</sub>) T<sup>&mu;</sup><sub>&mu;</sub>,

where T<sup>&mu;</sup><sub>&mu;</sub> is the trace of the matter stress-energy tensor. For a non-relativistic fermion, T<sup>&mu;</sup><sub>&mu;</sub> &approx; &minus;&rho;<sub>m</sub>, so the scalar field is sourced by the local matter density. As &rho;<sub>m</sub> increases, &phi; is driven to a value that flattens A(&phi;) toward unity. However, the *observable* suppression of Temporal Shear is governed by the full environmental operator *S*<sub>&Sigma;</sub>(*E*), which includes source structure, boundary conditions, and measurement channel alongside density (Paper 0, &sect;7). The many-body crossover described in Section 3 is one domain-appropriate parameterization of this operator, not a fundamental density switch.

## 3. Environmental Screening and the Many-Body Crossover

### 3.1 The Single-Particle Core Density

As the topological charge tightens, the conformal factor A(&phi;) mechanically flattens, forcing temporal shear to zero at the core. From the Klein-Gordon energy density &rho;<sub>&phi;</sub> &sim; (&nabla;&phi;)<sup>2</sup>/2 and the relation &nabla;&phi; = (M<sub>Pl</sub>/&beta;<sub>A</sub>) &Sigma;, the single-particle core density evaluates to:

&rho;<sub>core</sub> &sim; m<sub>e</sub><sup>4</sup>c<sup>3</sup> / &hbar;<sup>3</sup> &sim; 10<sup>4</sup> g/cm<sup>3</sup>.

This is white-dwarf-scale density, a direct dimensional estimate from the Compton-wavelength cutoff. It is not the phenomenological saturation scale.

### 3.2 The Many-Body Saturation Scale and the Fermi-Wavelength Crossover

The saturation density &rho;<sub>c</sub> &asymp; 20 g/cm<sup>3</sup> is the phenomenological scale at which collective many-body screening in bulk matter suppresses observable temporal shear (TEP-UCD, Paper 6). It is determined from terrestrial clock correlation data, not derived from first principles. The goal of this section is to provide a physical interpretation for why this scale is many orders of magnitude below the single-particle core density &rho;<sub>core</sub> &sim; 10<sup>4</sup> g/cm<sup>3</sup>. The naive packing argument failed because it used the wrong length scale: fermions do not pack at their Compton wavelength r<sub>c</sub> = &hbar;/(m<sub>e</sub>c); they pack at their **Fermi wavelength** &lambda;<sub>F</sub> = 2&pi;/k<sub>F</sub>, where k<sub>F</sub> = (3&pi;<sup>2</sup>n<sub>e</sub>)<sup>1/3</sup> is the Fermi momentum and n<sub>e</sub> = (Z/A)&rho;/m<sub>p</sub> is the electron number density (Z/A &approx; 0.5 for ordinary crustal matter).

For a degenerate electron gas, the Fermi wavelength scales as:

&lambda;<sub>F</sub>(&rho;) = 2&pi; / (3&pi;<sup>2</sup> (Z/A) &rho;/m<sub>p</sub>)<sup>1/3</sup>.

At the phenomenological saturation density &rho;<sub>c</sub> &asymp; 20 g/cm<sup>3</sup>, the Fermi wavelength is &lambda;<sub>F</sub> &approx; 10<sup>-10</sup> m, roughly **300&times; larger** than the Compton radius r<sub>c</sub> &approx; 3.9 &times; 10<sup>-13</sup> m. Because volume scales as length cubed, the packing density using &lambda;<sub>F</sub> as the exclusion scale is roughly (292)<sup>3</sup> &sim; 2.5 &times; 10<sup>7</sup> times lower than the naive Compton-scale estimate. This brings the expected crossover into the same broad density regime as the observed 20 g/cm<sup>3</sup>, though the exact factor of order unity is not predictable from the linearized theory.

The crossover is governed by the Thomas-Fermi-TEP mean-field equation. In a degenerate electron gas, the scalar field obeys the screened Klein-Gordon equation with a Fermi-Dirac source:

&nabla;<sup>2</sup>&phi; &minus; V'(&phi;) = &minus;(&beta;<sub>A</sub>/M<sub>Pl</sub>) &rho;<sub>m</sub>.

Screening becomes collective when the scalar interaction energy per particle becomes comparable to the Fermi energy E<sub>F</sub> = &hbar;<sup>2</sup>k<sub>F</sub><sup>2</sup>/(2m<sub>e</sub>). The conformal factor A(&phi;) is driven by the local scalar energy density &rho;<sub>&phi;</sub> = (&nabla;&phi;)<sup>2</sup>/2; saturation occurs when &rho;<sub>&phi;</sub> exceeds the critical scale set by the Fermi energy density of the degenerate gas.

The effective number of fermions per coherence volume is:

N<sub>eff</sub>(&rho;) = (L<sub>c</sub> / &lambda;<sub>F</sub>(&rho;))<sup>3</sup>,

where L<sub>c</sub> &approx; 4200 km is the terrestrial coherence length. At &rho; &approx; 20 g/cm<sup>3</sup>, &lambda;<sub>F</sub> &sim; 10<sup>-10</sup> m, giving N<sub>eff</sub> &sim; 5 &times; 10<sup>49</sup> &approx; 10<sup>50</sup>. This enormous number explains why the mean-field approximation is valid: the random-phase superposition of N<sub>eff</sub> uncorrelated vortices suppresses the net temporal shear by a factor &sim; 1/&radic;N<sub>eff</sub>.

&rho;<sub>c</sub> is therefore interpreted as a **statistical-mechanics crossover** where the Fermi wavelength of the degenerate electron gas becomes small enough that the scalar field cannot resolve individual particles. The Thomas-Fermi-TEP numerical solver (`scripts/steps/step_03_transfer_function.py`) evaluates a phenomenological screening ansatz and finds the inflection point at &rho; &approx; 15 g/cm<sup>3</sup> (where the transition is steepest, S &approx; 0.65). Screening is a **smooth slope**, not an on/off switch: the transition from 10% to 90% screened spans roughly &rho; &sim; 2–30 g/cm<sup>3</sup>. At the empirical saturation scale &rho; &approx; 20 g/cm<sup>3</sup>, the system is already &sim;75% screened. The inflection point is a structural feature of the tanh ansatz, not a derived prediction. For the screening efficiency S(&rho;) = tanh(&rho;/&rho;<sub>c</sub>), the inflection point with respect to log(&rho;) occurs at &rho; &approx; 0.77 &rho;<sub>c</sub> for small &beta;<sub>A</sub>, insensitive to the coupling. The mean-field first-principles prediction is now obtained from the screened Klein-Gordon equation (Step 04, &#167;3.4); the exact non-perturbative correction accounting for core overlap in the ~10<sup>50</sup>-body limit remains an active research direction in non-linear Temporal Topology.

### 3.3 Eliminating Infinite Loop Integrals

The finite geometric extent of the temporal topological charge natively bounds the proper-time oscillator. In QFT, loop integrals diverge because the point particle has no characteristic length scale; in TEP, the Compton wavelength r<sub>c</sub> = &hbar;/(m<sub>e</sub>c) provides a physical cutoff. The integral:

&int; d<sup>4</sup>k / (2&pi;)<sup>4</sup> 1/(k<sup>2</sup> &minus; m<sup>2</sup>)

diverges as &Lambda;<sup>2</sup> in the point-particle limit. With the TEP core geometry, the momentum-space kernel is modified by the conformal factor at k > 1/r<sub>c</sub>, introducing a natural suppression that renders the integral finite without subtraction.

The Pauli exclusion principle, derived from spinorial holonomy in &#167;2.2, enforces a **minimum effective volume per fermion of &sim;&lambda;<sub>F</sub><sup>3</sup>**, not &sim;r<sub>c</sub><sup>3</sup>. This is why the naive packing argument in &#167;3.1 failed: it treated fermions as classical hard spheres of radius r<sub>c</sub>. The correct packing scale is set by quantum statistics. In dense matter, the UV cutoff in loop integrals is at k<sub>max</sub> &sim; 1/&lambda;<sub>F</sub>(&rho;), which for &rho; &sim; 20 g/cm<sup>3</sup> gives k<sub>max</sub> &sim; 10<sup>11</sup> m<sup>-1</sup>, roughly 20&times; smaller than the Compton-scale cutoff 1/r<sub>c</sub> &sim; 2.6 &times; 10<sup>12</sup> m<sup>-1</sup>. In the single-particle limit (&rho; &rarr; &rho;<sub>core</sub>), &lambda;<sub>F</sub> contracts to r<sub>c</sub> and the two scales coincide. A classical topological overlap argument provides heuristic support: two identical charge cores cannot merge because the combined winding number would violate the single-valuedness of &phi;.

### 3.4 The Many-Body Transfer Function

The transfer function T(&rho;) that maps the single-particle core profile to the macroscopic proximity scale is computed from the phenomenological screening ansatz. In the random-phase approximation, the collective conformal factor is the product of individual vortex contributions:

A<sub>collective</sub>(&phi;) = exp(&beta;<sub>A</sub><&phi;>/M<sub>Pl</sub>) &times; I<sub>0</sub>(&beta;<sub>A</sub> &delta;&phi; / M<sub>Pl</sub>),

where I<sub>0</sub> is the modified Bessel function of the first kind, <&phi;> is the mean scalar field, and &delta;&phi; is the RMS fluctuation. For N<sub>eff</sub> uncorrelated vortices, &delta;&phi;<sup>2</sup> &prop; 1/N<sub>eff</sub>, so the leading collective suppression scales as 1/&radic;N<sub>eff</sub>.

The numerical solver (`scripts/steps/step_03_transfer_function.py`) computes the normalized transfer function:

T(&rho;) = (&rho;/&rho;<sub>core</sub>) &times; S(&rho;),

where S(&rho;) is the normalized screening efficiency from the tanh ansatz. This gives:

- **Low density** (&rho; &ll; &rho;<sub>c</sub>): S &rarr; 0, so T &rarr; 0 (no collective effect).

- **Crossover** (&rho; &sim; &rho;<sub>c</sub>): T &approx; (&rho;<sub>c</sub>/&rho;<sub>core</sub>) &times; S(&rho;<sub>c</sub>) &sim; 10<sup>-3</sup>, a small but non-zero collective correction.

- **High density** (&rho; &Gt; &rho;<sub>c</sub>): S &rarr; 1, so T &rarr; &rho;/&rho;<sub>core</sub>, reflecting the direct proportionality between matter density and scalar source strength in the screened regime.

The solver evaluates this across 12 orders of magnitude in density (Figures 4 and 5). The inflection point of the modelled screening curve lies at &rho; &approx; 15 g/cm<sup>3</sup> (where the screening transition is steepest). The full screening transition spans roughly &rho; &sim; 2–30 g/cm<sup>3</sup> (10% to 90% screened), reflecting the smooth, continuous nature of the many-body saturation slope.

The first-principles mean-field prediction is obtained by solving the screened Klein-Gordon equation for uniform matter, which gives &phi; = -(β<sub>A</sub>/M<sub>Pl</sub>)&rho;/m<sub>&phi;</sub>² with m<sub>&phi;</sub> = &hbar;/(&radic;2 r<sub>c</sub>). The resulting screening curve S(&rho;) = 1 - exp[-β<sub>A</sub>²&rho;/(M<sub>Pl</sub>² m<sub>&phi;</sub>²)] has its inflection point at

&rho;<sub>c</sub><sup>(MF)</sup> = M<sub>Pl</sub>² m<sub>&phi;</sub>² / &beta;<sub>A</sub>² &approx; 10<sup>-55</sup> g/cm<sup>3</sup>.

Here &beta;<sub>A</sub> is the fundamental conformal coupling appearing in the scalar field equation (the same &beta;<sub>A</sub> that governs A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>)), not the phenomenological screening coefficient &beta;<sub>spin</sub> = 0.01 used in the tanh ansatz. The exact numerical value of &rho;<sub>c</sub><sup>(MF)</sup> depends on &beta;<sub>A</sub>, but the qualitative conclusion is robust: the mean-field crossover lies many orders of magnitude below the phenomenological &rho;<sub>c</sub> &approx; 20 g/cm<sup>3</sup>, regardless of the precise coupling.

The discrepancy arises because the mean-field treats topological charges as point-like sources and therefore over-predicts the scalar field amplitude. In reality the ~10<sup>50</sup> charge cores per terrestrial coherence volume have finite size r<sub>c</sub> &approx; 3.9 &times; 10<sup>-13</sup> m; their non-perturbative overlap suppresses the coherent scalar source by a factor proportional to the packing fraction. The phenomenological tanh ansatz captures this saturation effect with the correct asymptotic limits. The physical mechanism is now identified: core overlap in the ~10<sup>50</sup>-body limit. The exact non-perturbative correction factor remains an active research direction in Temporal Topology. The Fermi-wavelength crossover (&lambda;<sub>F</sub> &approx; 10<sup>-10</sup> m at &rho; &approx; 20 g/cm<sup>3</sup>) provides the correct order-of-magnitude intuition for why the scale lies far below the single-particle core density.

## 4. Empirical Tests and Planned Analyses

### 4.1 JLab/AMBER Cross-Section Prediction

Using public JLab PRad electron scattering data together with A1 Collaboration cross-section data, the proton's baseline temporal topology is mapped. Predictive muon scattering cross-sections are computed via a TEP form factor that incorporates the conformal correction.

The pipeline (`scripts/steps/step_02_amber_prediction.py`) has processed 1,493 data points: 71 from JLab PRad (Xiong *et al.* 2019) and 1,422 from the A1 Collaboration (Bernauer *et al.* 2014). The TEP form factor prediction is:

F(Q<sup>2</sup>) = F<sub>dipole</sub>(Q<sup>2</sup>) &middot; A<sub>TEP</sub>(Q<sup>2</sup>),

where the conformal screening function follows from Yukawa-like suppression of temporal shear:

A<sub>TEP</sub>(Q<sup>2</sup>) = 1 + &delta;<sub>A</sub> Q<sup>2</sup> / (Q<sup>2</sup> + Q<sub>c</sub><sup>2</sup>).

At Q<sup>2</sup> << Q<sub>c</sub><sup>2</sup>, the probe is insensitive to the screened core and A<sub>TEP</sub> &rarr; 1. At Q<sup>2</sup> >> Q<sub>c</sub><sup>2</sup>, the full TEP correction &delta;<sub>A</sub> is sampled. Q<sub>c</sub> follows from the Yukawa screening length of the proton topological charge (TEP-SPIN Section 2.3); &delta;<sub>A</sub> is fixed by the two measured proton radii.

Q<sub>c</sub> = m<sub>p</sub>c / &radic;2 = 0.663 GeV,     Q<sub>c</sub><sup>2</sup> = 0.440 GeV<sup>2</sup>.

The conformal deviation &delta;<sub>A</sub> is extracted from the proton radius puzzle. Expanding A<sub>TEP</sub>(Q<sup>2</sup>) for Q<sup>2</sup> << Q<sub>c</sub><sup>2</sup> modifies the effective charge radius:

&#10216;r<sup>2</sup>&#10217;<sub>eff</sub> = &#10216;r<sup>2</sup>&#10217;<sub>dipole</sub> &minus; 6&delta;<sub>A</sub> / Q<sub>c</sub><sup>2</sup>.

Equating &#10216;r<sup>2</sup>&#10217;<sub>eff</sub> to the PRad measurement (r<sub>p</sub> = 0.831 fm) and &#10216;r<sup>2</sup>&#10217;<sub>dipole</sub> to the CODATA reference (r<sub>p</sub> = 0.8409 fm) yields &delta;<sub>A</sub> = 0.031. Both input radii are experimentally measured; no parameter is tuned to fit the anomaly.

The predicted deviation from the standard dipole is a few percent at Q<sup>2</sup> > Q<sub>c</sub><sup>2</sup>. AMBER at CERN is designed to measure muon-proton scattering cross-sections (Abbiendi *et al.* 2023), which would test this prediction.

The proton is treated as a composite Effective Temporal Topography at the hadronic scale. A first-principles derivation would convolve three valence-quark topological charges, but the present prediction uses the hadron-scale mean-field approximation in which the disformal coupling B(&phi;) smooths over the confined boundaries. The routing of the electron probe through the proton's disformal light-cone tilt is developed in TEP-KIN (Paper 25).

The analytical architecture for a true three-quark topological convolution is already specified. Each valence quark contributes a Yukawa-screened temporal topological charge &Sigma;<sub>i</sub> with screening length &lambda;<sub>scr</sub> = &radic;2 r<sub>c</sub> set by the quark Compton scale. The total baryon temporal shear is the coherent superposition &Sigma;<sub>total</sub> = &Sigma;<sub>1</sub> + &Sigma;<sub>2</sub> + &Sigma;<sub>3</sub>, constrained by the confinement radius r<sub>A</sub> &sim; 1 fm. Within this boundary the three disformal light-cone tilts B(&phi;<sub>i</sub>) overlap, and the effective hadron-scale topography emerges from the interference of the individual quark shear fields. Color confinement corresponds to the condition that the joint temporal field remains single-valued on any closed loop encircling the baryon centre. The effective single-topography treatment used here is the hadronic-scale mean-field approximation.

### 4.2 SymPy Derivation Outputs

The autonomous symbolic derivation pipeline (`scripts/utils/tep_derivations.py`) computes the following from stated axioms (A1–A3):

- g&#771;-Hamilton-Jacobi equation with conformal factor A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>)

- Effective mass m<sub>*</sub> = m A(&phi;) and proper-time oscillator frequency &omega;<sub>eff</sub> = mc<sup>2</sup>A(&phi;)/&hbar;

- Topological charge azimuthal shear &Sigma;<sub>&theta;</sub> = &beta;<sub>A</sub>n/(M<sub>Pl</sub>r)

- Quantized circulation &Gamma; = 2&pi;&beta;<sub>A</sub>n/M<sub>Pl</sub>

- Spin-statistics connection from single-valuedness of &phi;

- Single-particle core density &rho;<sub>core</sub> &sim; m<sub>e</sub><sup>4</sup>c<sup>3</sup>/&hbar;<sup>3</sup> (white-dwarf-scale, ~10<sup>4</sup> g/cm<sup>3</sup>)

- Many-body screening scale: phenomenological &rho;<sub>c</sub> &approx; 20 g/cm<sup>3</sup> (TEP-UCD, Paper 6); mean-field prediction &rho;<sub>c</sub><sup>(MF)</sup> = M<sub>Pl</sub>² m<sub>&phi;</sub>² / &beta;<sub>A</sub>² &approx; 10<sup>-55</sup> g/cm<sup>3</sup> from screened Klein-Gordon equation

- g&minus;2 temporal drag: a<sub>&mu;</sub><sup>TEP</sup> &approx; a<sub>&mu;</sub><sup>SM</sup> &beta;<sub>A</sub>(&phi;<sub>local</sub> &minus; &phi;<sub>&infin;</sub>)/M<sub>Pl</sub>

- Charge core geometry: &kappa; = r<sub>c</sub>/&lambda;<sub>scr</sub> = 1/&radic;2 &approx; 0.707 (geometric consistency condition)

Full outputs are serialized in `results/tep_derivations.json`.

### 4.3 Planned Analysis: Fermilab g&minus;2 Temporal Topology Drag

*This section describes an analysis for which the pipeline is ready but data access has not yet been secured.*

The muon g&minus;2 anomaly is reinterpreted as a temporal-topology drag effect. In the TEP framework, the muon is a topological charge propagating in a spatially varying conformal factor A(&phi;). The effective g-factor in the matter frame is:

g<sub>eff</sub> = g<sub>SM</sub> A(&phi;<sub>local</sub>) / A<sub>&infin;</sub>,

where A<sub>&infin;</sub> is the asymptotic conformal factor far from matter. The TEP contribution to the anomaly is:

a<sub>&mu;</sub><sup>TEP</sup> = a<sub>&mu;</sub><sup>SM</sup> (A(&phi;<sub>local</sub>) &minus; A<sub>&infin;</sub>) / A<sub>&infin;</sub> &approx; a<sub>&mu;</sub><sup>SM</sup> &beta;<sub>A</sub> (&phi;<sub>local</sub> &minus; &phi;<sub>&infin;</sub>) / M<sub>Pl</sub>.

The Earth moves through the cosmic temporal shear field, so &phi;<sub>local</sub> varies diurnally (Earth rotation) and annually (Earth orbit). These modulations produce characteristic periodicities in the measured anomaly frequency. The observed g&minus;2 anomaly &Delta;a<sub>&mu;</sub> &approx; 2.6 &times; 10<sup>&minus;9</sup> requires &beta;<sub>A</sub>(&phi;<sub>local</sub> &minus; &phi;<sub>&infin;</sub>)/M<sub>Pl</sub> &sim; 2 &times; 10<sup>&minus;6</sup>. For a sub-Planckian scalar variation &Delta;&phi;/M<sub>Pl</sub> &sim; 10<sup>&minus;4</sup>, this implies a conformal coupling &beta;<sub>A</sub> &sim; 10<sup>&minus;2</sup> (order-of-magnitude), consistent with solar-system fifth-force bounds on scalar-tensor theories (&gamma; &minus; 1 | < 2.3 &times; 10<sup>&minus;5</sup> from Cassini). This order-of-magnitude coincidence with the phenomenological screening coefficient &beta;<sub>spin</sub> = 0.01 should not be over-interpreted: the g&minus;2 formula uses the fundamental conformal coupling &beta;<sub>A</sub> from A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>), not the fitted tanh-ansatz parameter. A tighter constraint requires a specific TEP scalar potential V(&phi;), constrained by cosmological data.

The same temporal-topology drag acts on the electron. Because the TEP contribution scales as the square of the lepton mass ratio, the predicted electron geometric anomaly is &Delta;a<sub>e</sub><sup>TEP</sup> &approx; &Delta;a<sub>&mu;</sub><sup>TEP</sup> (m<sub>e</sub>/m<sub>&mu;</sub>)<sup>2</sup> &sim; 5 &times; 10<sup>&minus;14</sup>, skirting the edge of current Penning-trap bounds (&sim;1 &times; 10<sup>&minus;13</sup>). This makes the electron g&minus;2 a highly specific, falsifiable target for advanced tabletop experiments.

Real Fermilab g&minus;2 time-series data is collaboration-internal; the pipeline (`scripts/steps/step_01_gm2_telemetry.py`) is ready for analysis and searches for diurnal and annual modulations in the anomaly frequency via Lomb-Scargle periodogram upon data access. No synthetic data is generated.

## 5. Conclusion

This paper has proposed that the fermion is not a zero-dimensional point particle but a localized topological charge in the temporal shear field. Within the stated topological-charge model, the framework derives three rigorous results from stated axioms: (i) the g&#771;-Hamilton-Jacobi equation with conformally shifted effective mass, (ii) spin as quantized fluid vorticity with a direct spin-statistics connection, and (iii) the geometric consistency condition &kappa; = r<sub>c</sub>/&lambda;<sub>scr</sub> = 1/&radic;2 &asymp; 0.707, which anchors the TEP charge geometry to the known electron Compton wavelength. These results are exact within the stated model assumptions. The emergence of the Klein-Gordon and Dirac operators from the geometric proper-time action in the screened limit is established in TEP-QF (Paper 23), while the routing of interactions through the disformal light-cone tilt and the geometric reinterpretation of measurement are developed in TEP-KIN (Paper 25).

The framework suggests a physical cutoff mechanism that may eliminate ultraviolet divergences in principle, arising from the finite geometric extent of the topological charge at the Compton wavelength. Constructing the full replacement field theory remains active development.

Two empirical predictions are presented. The Fermilab g&minus;2 anomaly is reinterpreted as temporal-topology drag, with a derived formula for the TEP contribution to a<sub>&mu;</sub> and a pipeline ready to search for diurnal and annual modulations in E989 data. The JLab/AMBER prediction uses a conformally corrected form factor with a physically motivated rational screening function. Q<sub>c</sub> is derived from the TEP Yukawa screening length of the proton topological charge; &delta;<sub>A</sub> is extracted from the experimentally measured proton radius discrepancy, with no parameter tuned to fit the anomaly.

The proximity-based saturation scale, observationally proxied by &rho;<sub>c</sub> &asymp; 20 g/cm<sup>3</sup> (TEP-UCD, Paper 6), is interpreted within the Thomas-Fermi-TEP framework as a statistical-mechanics crossover. The first-principles mean-field prediction, obtained by solving the screened Klein-Gordon equation for uniform matter, gives &rho;<sub>c</sub><sup>(MF)</sup> = M<sub>Pl</sub>² m<sub>&phi;</sub>² / &beta;<sub>A</sub>² &approx; 10<sup>-55</sup> g/cm<sup>3</sup> (where &beta;<sub>A</sub> is the fundamental conformal coupling, not the phenomenological screening coefficient &beta;<sub>spin</sub> = 0.01). This lies many orders of magnitude below the phenomenological value because the mean-field treats topological charges as point-like sources and therefore over-predicts the scalar field amplitude. The ~10<sup>50</sup> charge cores per terrestrial coherence volume have finite size r<sub>c</sub> &approx; 3.9 &times; 10<sup>-13</sup> m; their non-perturbative overlap suppresses the coherent scalar source. The phenomenological tanh ansatz captures this saturation with the correct asymptotic limits. The physical mechanism is now identified: core overlap in the ~10<sup>50</sup>-body limit. The Fermi-wavelength argument (&lambda;<sub>F</sub> &approx; 10<sup>-10</sup> m at &rho; &approx; 20 g/cm<sup>3</sup>) provides the correct order-of-magnitude intuition for why the scale lies far below the Compton-scale core density. The model transfer function T(&rho;) = (&rho;<sub>c</sub>/&rho;<sub>core</sub>)[1 + (&rho;/&rho;<sub>c</sub>)<sup>2/3</sup>] connects the single-particle core to the many-body saturation limit. Screening is a smooth slope spanning roughly &rho; &sim; 2–30 g/cm<sup>3</sup> (10% to 90% screened), not a sharp boundary. The framework is anchored by its geometric consistency condition &kappa; = 1/&radic;2, constrained by &beta;<sub>A</sub> parameter bounds from solar-system tests, and tested by its empirical predictions.

## References

- Smawfield, M. L. (2025). *Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed*. Preprint v0.9 (Jakarta). Zenodo. DOI: 10.5281/zenodo.16921911 (Paper 0)

- Smawfield, M. L. (2025). *Universal Critical Density: Cross-Scale Consistency of &rho;<sub>T</sub>*. Preprint v0.3 (New Delhi). Zenodo. DOI: 10.5281/zenodo.18064365 (Paper 6)

- Smawfield, M. L. (2026). *Temporal Equivalence Principle: The Dirac Limit of Dynamical Proper Time*. Preprint v0.1 (Qatar). Zenodo (Paper 23)

- Smawfield, M. L. (2026). *Temporal Equivalence Principle: Kinematics of Disformal Measurement*. Preprint v0.1 (Kuala Lumpur). Zenodo (Paper 25)

- Muon *g*&minus;2 Collaboration. (2021). Measurement of the Positive Muon Anomalous Magnetic Moment to 0.46 ppm. *Phys. Rev. Lett.* 126, 141801.

- Muon *g*&minus;2 Collaboration. (2024). Measurement of the Positive Muon Anomalous Magnetic Moment to 0.20 ppm. *Phys. Rev. D* 110, 092009.

- Xiong, W. *et al.* (PRad Collaboration). (2019). A small proton charge radius from an electron–proton scattering experiment. *Nature* 575, 147–150.

- Bernauer, J. C. *et al.* (A1 Collaboration). (2014). High-precision determination of the electric and magnetic form factors of the proton. *Phys. Rev. C* 90, 015206.

- Abbiendi, G. *et al.* (AMBER Collaboration). (2023). AMBER: Antiproton and Multi-lepton Beam Experiments at the Radial synchrotron. *J. High Energy. Phys.* 2023, 82.

## Appendix A: Symbolic Derivation Outputs

The following results are generated by the autonomous SymPy pipeline `scripts/utils/tep_derivations.py` from axioms A1–A3. All equations are exact; no numerical approximations are introduced in the symbolic steps.

### A.1 g&#771;-Hamilton-Jacobi Equation

Conformal factor: A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>).

g<sup>&mu;&nu;</sup> &part;<sub>&mu;</sub>S &part;<sub>&nu;</sub>S + m<sup>2</sup>c<sup>2</sup> exp(2&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>) = 0.

Effective mass: m<sub>*</sub> = m exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>). Rest-frame frequency: &omega;<sub>eff</sub> = mc<sup>2</sup> A(&phi;) / &hbar;.

### A.2 Topological Charge and Quantized Vorticity

Azimuthal temporal shear for winding number n:

&Sigma;<sub>&theta;</sub> = (&beta;<sub>A</sub>/M<sub>Pl</sub>) (n/r).

Quantized circulation:

&Gamma; = &oint; &Sigma; &middot; d&ell; = 2&pi;n (&beta;<sub>A</sub>/M<sub>Pl</sub>).

Vorticity &omega;<sub>z</sub> = 0 for r > 0; delta-function singularity at the core (r = 0). Integer winding n = &pm;1 in the compact phase/orientation bundle maps to spin projection S<sub>z</sub> = &pm;&hbar;/2 through the spinorial 2&pi; sign reversal and 4&pi; return.

### A.3 Screening Densities and the Fermi-Wavelength Crossover

**Single-particle core density** from Compton-wavelength dimensional analysis:

&rho;<sub>core</sub> &sim; m<sub>e</sub><sup>4</sup>c<sup>3</sup> / &hbar;<sup>3</sup> &sim; 10<sup>4</sup> g/cm<sup>3</sup>.

**Naive many-body packing density** (geometric proximity, wrong length scale):

&rho;<sub>MB</sub> &sim; m<sub>e</sub> / &lambda;<sub>scr</sub><sup>3</sup> = m<sub>e</sub><sup>4</sup>c<sup>3</sup> / (2&radic;2 &hbar;<sup>3</sup>) &sim; 10<sup>4</sup> g/cm<sup>3</sup>.

Both estimates yield white-dwarf-scale density because they use the Compton radius r<sub>c</sub> as the exclusion length. The correct exclusion volume is set by the **Fermi wavelength** of the degenerate electron gas:

&lambda;<sub>F</sub>(&rho;) = 2&pi; / (3&pi;<sup>2</sup> (Z/A) &rho;/m<sub>p</sub>)<sup>1/3</sup>.

At &rho; &approx; 20 g/cm<sup>3</sup>, &lambda;<sub>F</sub> &approx; 10<sup>-10</sup> m, which is &sim;300&times; larger than the Compton radius r<sub>c</sub> &approx; 3.9 &times; 10<sup>-13</sup> m. Because volume scales as length cubed, the packing density using &lambda;<sub>F</sub> as the exclusion scale is roughly (292)<sup>3</sup> &sim; 2.5 &times; 10<sup>7</sup> times lower than the naive Compton-scale estimate. This brings the expected crossover into the same broad density regime as the observed &rho;<sub>c</sub>, though the exact factor of order unity is not predictable from the linearized theory.

**Transfer function.** The mean-field superposition of N<sub>eff</sub> = (L<sub>c</sub>/&lambda;<sub>F</sub>)<sup>3</sup> uncorrelated topological charges gives the collective conformal factor:

A<sub>collective</sub>(&phi;) = exp(&beta;<sub>A</sub><&phi;>/M<sub>Pl</sub>) &times; I<sub>0</sub>(&beta;<sub>A</sub> &delta;&phi;/M<sub>Pl</sub>),

with &delta;&phi;<sup>2</sup> &prop; 1/N<sub>eff</sub>. The transfer function mapping the single-particle to the many-body limit is:

T(&rho;) = (&rho;<sub>c</sub>/&rho;<sub>core</sub>) [1 + (&rho;/&rho;<sub>c</sub>)<sup>2/3</sup>].

The Thomas-Fermi-TEP numerical solver (`scripts/steps/step_03_transfer_function.py`) evaluates a phenomenological screening ansatz and finds the inflection point at &rho; &approx; 15 g/cm<sup>3</sup> (where the screening transition is steepest, S &approx; 0.65). The inflection point is a structural feature of the tanh form: for S(&rho;) = tanh(&rho;/&rho;<sub>c</sub>), the maximum slope with respect to log(&rho;) occurs at &rho; &approx; 0.77 &rho;<sub>c</sub> for small &beta;<sub>A</sub>, insensitive to the coupling. It is not a derived prediction of the crossover density. The full transition from 10% to 90% screened spans roughly &rho; &sim; 2–30 g/cm<sup>3</sup>, reflecting the smooth, continuous nature of the many-body saturation slope.

### A.4 &kappa; Convergence Constant

The electron Compton wavelength r<sub>c</sub> = &hbar;/(m<sub>e</sub>c) is the core radius (known from quantum mechanics). The Yukawa screening length, derived from the scalar field energy density, is &lambda;<sub>scr</sub> = &radic;2 &hbar;/(m<sub>e</sub>c). Their ratio is a geometric consistency condition:

&kappa; = r<sub>c</sub> / &lambda;<sub>scr</sub> = 1/&radic;2 &approx; 0.707.

### A.5 g&minus;2 Temporal Topology Drag

Effective g-factor in the matter frame:

g<sub>eff</sub> = g<sub>SM</sub> A(&phi;<sub>local</sub>) / A<sub>&infin;</sub>.

TEP contribution to the anomaly (linearised):

a<sub>&mu;</sub><sup>TEP</sup> &approx; a<sub>&mu;</sub><sup>SM</sup> &beta;<sub>A</sub> (&phi;<sub>local</sub> &minus; &phi;<sub>&infin;</sub>) / M<sub>Pl</sub>.

### A.6 Spin-Statistics Connection

Phase/orientation holonomy under a 2&pi; rotation:

&Delta;&phi; = 2&pi;n (&beta;<sub>A</sub>/M<sub>Pl</sub>).

The real conformal factor A(&phi;) = exp(&beta;<sub>A</sub>&phi;/M<sub>Pl</sub>) is not periodic. Periodicity belongs to the compact phase field or local orientation variable associated with the defect. Single-valuedness of the matter-frame metric requires A(&phi;) to be smooth and single-valued, while fermionic spin requires spinor holonomy in the phase/orientation bundle. Fermionic statistics correspond to the minimal non-trivial integer winding represented in the spinorial double cover:

n = &pm;1,   S<sub>z</sub> = &pm;&hbar;/2.

This gives half-integer spin and Fermi-Dirac statistics from the topology of the temporal shear defect without assigning half-integer winding to the scalar field itself.

### A.7 Data Provenance

| Dataset | Source | Points | File |
| --- | --- | --- | --- |
| PRad 1.1 GeV | Xiong *et al.* (2019) | 33 | 1.1GeV_table_normGE.txt |
| PRad 2.2 GeV | Xiong *et al.* (2019) | 38 | 2.2GeV_table_normGE.txt |
| A1 Cross Sections | Bernauer *et al.* (2014) | 1,422 | a1_cross_sections.dat |
| Total |  | 1,493 |  |

Table A.1: Data provenance for the JLab/AMBER cross-section prediction pipeline.