## TEP Parameter Registry

Generated from parameter_registry.yaml. Classification records provenance; a benchmark is not an empirical calibration. Numerical values and descriptive claims retain the YAML's declared scope; this table is not independent validation.

| Symbol | Meaning | Value | Units | Class | First introduced | Used in |
|---|---|---|---|---|---|---|
| $\beta_A$ | Conformal coupling strength in A(phi) = exp(beta_A * phi / M_Pl) | -1.0 | dimensionless | fundamental | Paper 0 | 0, 1, 2, 3, 4, 6, 10, 11, 12, 14, 15, 17, 18, 19, 22, 23, 24, 25, 26, 28, 29, 31, 32 |
| $\rho_T$ | Temporal Topology saturation scale (observationally proxied by density) | 20.0 | g cm^{-3} | calibrated | Paper 6 | 6, 14, 17, 18, 23, 24 |
| $\lambda_T$ | Clock-correlation length (Temporal Topology coherence length) | 4200.0 | km | calibrated | Paper 1 | 1, 2, 3, 6, 14, 17 |
| $\alpha_{\log}$ | Lab-scale density-sector coupling | -0.00766 | dimensionless | exploratory | retired (withdrawn programme) | none |
| $\beta_{\text{geom}}$ | Mass-sector geometric coupling | 0.00015 | dimensionless | exploratory | retired (withdrawn programme) | none |
| $M_{\text{ref}}$ | Reference mass scale for geometric coupling | 1.0e18 | kg | exploratory | retired (withdrawn programme) | none |
| $\kappa_{\text{canonical}}$ | Canonical galaxy-scale observable response coefficient | 9.6e5 | mag | benchmark | Paper 11 | 11, 12 |
| $\kappa_{\text{nested}}$ | Nested clock-ratio bound on the raw conformal channel | -0.00075 | mag | derived | Paper 11 | 0, 11 |
| $\kappa_{\text{Cep}}$ | Cepheid observable response coefficient | Not specified | dimensionless | observable_response | Paper 11 | 11, 12 |
| $\kappa_{\text{MSP}}$ | Pulsar in-equation acceleration coupling (screened) | 0.05 | dimensionless | observable_response | Paper 10 | 10, 17 |
| $\tilde\kappa_{\text{MSP}}$ | Full empirical pulsar response coefficient (pipeline inversion) | 29000.0 ± 45000.0 | dimensionless | observable_response | Paper 10 | 10, 11, 6 |
| $\kappa_{\text{WB}}$ | Wide-binary observable response coefficient | Not specified | dimensionless | observable_response | Paper 0 | 0, 16 |
| $\kappa_{\text{LLR}}$ | LLR/flyby observable response coefficient | 1.59e-14 | dimensionless | derived | Paper 0 | 0, 8, 15, 17 |
| $\Gamma_X$ | Channel-specific geometric projector | Not specified | dimensionless | derived | Paper 0 | 0, 1, 2, 6, 10, 11, 12, 14, 15, 16, 17 |
| $\text{strong-field}$ | Strong-field hyperbolicity status | Not specified | dimensionless | derived | Paper 0 | 0 |
| $\text{quantum}$ | Quantum formalization status | Not specified | dimensionless | derived | Paper 0 | 0, 24, 25 |
| $\epsilon_T$ | Cosmology shear amplitude | Not specified | dimensionless | exploratory | Paper 26 | 18, 26 |
| $\phi_0$ | Present-day scalar field amplitude (phenomenological) | Not specified | dimensionless | exploratory | Paper 22 | 18, 22, 26 |
| $n$ | Redshift exponent in phenomenological conformal factor | Not specified | dimensionless | exploratory | Paper 22 | 18, 22 |
| $\beta_{\text{spin}}$ | Phenomenological screening coefficient in tanh ansatz | 0.01 | dimensionless | calibrated | Paper 24 | 24 |
| $\beta_{\text{Cassini}}$ | Upper bound on the unscreened DEF-normalized coupling \|α_DEF\| from the Cassini PPN test; the binding constraint on the screened Solar source charge is S_Σ(1 AU) ≲ 5.75×10⁻⁶ | 0.0034 | dimensionless | observable_response | Paper 0 | 0, 6, 18, 22, 24 |
| $B(\phi)$ | Disformal coupling function | Not specified | mass^{-4} for a dimension-one canonical scalar; dimensionless only in explicitly nondimensionalized benchmark coordinates | fundamental | Paper 0 | 0, 4, 9, 22, 28 |
| $H_0$ | Hubble constant (TEP-corrected local value) | 66.65 ± 1.58 | km s^{-1} Mpc^{-1} | calibrated | Paper 11 | 11, 26 |
| $S_8$ | Structure growth amplitude | Not specified | dimensionless | exploratory | Paper 26 | 26 |
| $V(\phi)$ | Canonical scalar potential branch | V(u) = V_matter(u) e^{-(u/u_s)^4} + V_0 e^{-(u_s/u)^4} | M_Pl^4 (dimensionless-u form) | fundamental | Paper 0 | 0 |
| $\Lambda_X$ | Kinetic-completion scale (gradient/shear screening threshold) | 1.9068 | meV | derived | Paper 0 | 0, 4, 6, 7, 12, 13, 17, 26, 28 |
| $\lambda_{\rm quartic}$ | Quartic self-coupling of the amplitude-sector potential | 7.526e-66 | dimensionless | calibrated | Paper 0 | 0, 6, 7, 10, 22, 28 |
| $\Lambda_V$ | Equivalent potential-normalization scale of the quartic branch | 127.5 | GeV | derived | Paper 6 | 0, 6 |
| $\rho_{\rm sat}$ | Self-quenching density (u_min ~ 1 crossing; temporal-well regime) | 6.1e25 | g cm^{-3} | derived | Paper 6 | 6, 28 |
| $u_{\rm min}$ | Density-set scalar-field equilibrium amplitude | Not specified | dimensionless | derived | Paper 0 | 0, 6, 7, 28 |
| $u_s$ | Master-potential plateau knee (self-quenching scale of the field amplitude) | 10.0 | dimensionless | benchmark | Paper 0 | 0, 6, 28 |

### Class definitions

- **fundamental** — Axiomatic or theoretically fixed parameter; changes would alter the framework itself.
- **calibrated** — Determined from experimental data with well-defined uncertainty; shared across analyses.
- **benchmark** — Declared numerical reference or prior; not an empirical calibration or independently derived value.
- **observable_response** — Domain-level empirical transfer coefficient; specific to a measurement domain.
- **nuisance** — Required for analysis but not a TEP parameter of direct physical interest.
- **exploratory** — Under investigation; value or existence not yet firmly established.
- **derived** — Structure or value computed from the action/pipeline; where the value additionally requires an environmental input or is a calibration identity, the description must say so.

### Parameter provenance and scope

#### beta_A

Locked lab-scale convention. Universal conformal coupling governing clock-rate rescaling.

#### rho_T

Phenomenological proximity scale at which collective many-body screening suppresses observable temporal shear, observationally proxied by the density rho_T ~ 20 g/cm^3. Corpus convention (decided 2026-09-28): rho_T = 20 g/cm^3 is the formal normalization adopted on the IGS precise-product branch (lambda_T ~ 4,200 km), not an independently measured constant. The fitted correlation-length family is estimator-dependent — precise-product bracket ~1.9-4.5x10^3 km (rho_T ~ 8-70 g/cm^3 effective band, epoch/window scatter dominant), raw channels bound below at ~0.7-1.1x10^3 km, MGEX held-out corroboration 1,862 +/- 112 km — so the numerical rho_T inherits the band of the adopted branch. The geometric identification lambda_T = R_T inverts the adopted 4,200 km to rho_T = 19.24 g/cm^3 — a 3.8% self-consistency residual that is a calibration identity, not an independent derivation (step_03). The atomic-scale (R_T(m_p) ~ a_0) and magnetar (P_crit ~ 6.8 s) anchors are estimator-independent cross-checks on the adopted normalization. Density is an observable proxy, not the causal parameter; the fundamental scale is geometric overlap governed by the packing fraction eta = (r_c/lambda_F)^3. Note: rho_c is reserved for per-cluster central density (Paper 10).

#### lambda_T

Canonical terrestrial correlation length for all forward analysis (terrestrial clock FEM, rho_T chain, etc.), adopted at 4,200 km on the IGS precise-product branch of the estimator family (25-year multi-center GNSS, Papers 1-2/6; long-span code 4,201 +/- 1,967 km). The fitted scale is product-dependent: precise-product bracket ~1.9-4.5x10^3 km with epoch/window scatter dominant; raw SPP channels ~0.7-1.1x10^3 km; raw ionofree 4,767 +/- 835 km. Paper 14 MGEX held-out verification (~1 yr combined-clock product, lambda = 1,862 +/- 112 km; annual-direction scan projection-degenerate and does not localize a CMB frame) corroborates a distance-structured signal on a different product but is not adopted for dimensional normalization. Cross-scale closure (steps 03-04): lambda_T is identified with the geometric saturation radius R_T = (3M/(4pi*rho_T))^{1/3}, not the Compton wavelength; the conditional Green-function covariance has its 1/e crossing near R_T (~4,003-4,012 km) without inserting the measured scale. Lab-scale crustal value is LAB_COHERENCE_LENGTH_M = 50 km.

#### alpha_log

Retired exploratory value from the withdrawn Naivasha lab-sector programme (no live paper). Negative by field-equation sign. Magnitude determined from requirement that TEP reproduce correct order of magnitude for laboratory metrology shifts. Former production use in Papers 25 (KIN 2DEG predictors) and 34 (PSR compact-object field) is resolved: both pipelines use the canonical local-equilibrium density response (core/scalar_field.py: equilibrium_u; Paper 34 uses the screened exterior charge u_scr(R)/u_unscr(R)). Legacy keyword arguments remain accepted for API compatibility but are ignored.

#### beta_geom

Retired exploratory parameter from the withdrawn Naivasha lab-sector programme (no live paper). Geometric coupling for mass-dependent scalar field contributions. Former production use in Paper 34 (PSR) is resolved: the exterior channel uses the canonical screened profile phi_profile_spherical with charge u_scr(R)/u_unscr(R). Legacy keyword arguments remain accepted for API compatibility but are ignored.

#### M_ref

Retired exploratory parameter from the withdrawn Naivasha lab-sector programme (no live paper). Threshold mass where phi_mass ~ beta_geom.

#### kappa_canonical

Expected bare-scale coupling for Cepheid P-L shifts before environmental screening. Declared theory benchmark used as prior in multi-anchor regression; the magnitude unit records the observational channel, not a physical dimension of the coupling. Step_05 derives the transfer structure (Gamma_gal inherited from the Cepheid projector with redshift transfer factor), but the numerical benchmark retains its declared provenance. Paper 12 uses the alias kappa_gal = kappa_canonical.

#### kappa_nested

Raw clock-ratio coefficient from the nested spectroscopic-to-Cepheid rate ratio on the uniform-profile benchmark (Paper 11, step_61): kappa_nested = (b/ln10) * Delta ln q / X = -7.5e-4 mag, delta_mu ~ 3e-10 mag at X = V_rot^2/c^2 for 200 km/s. The common galactic rate cancels in the ratio; alternative reference-clock conventions (volume-mean light depth, direct Cepheid-to-Cepheid baseline comparison) leave |kappa| of order unity at most. Establishes that the measured kappa_Cep is a response-sector coefficient, not a conformal clock-ratio artifact.

#### kappa_Cep

Domain-level transfer coefficient mapping TEP structure to Cepheid period-luminosity shifts. Fitted, not fundamental. Step_05 derives the transfer structure: kappa_Cep = |beta_A| S_Cep Gamma_Cep, where Gamma_Cep is the differential host-galaxy ambient field two-point function from the linearized scalar fluctuation equation — derived up to the environmental input (host-galaxy ambient field).

#### kappa_MSP

In-equation coupling in the pulsar spin-down equation; unsuppressed geometric amplitude ~10^6, screened to kappa_MSP^emp ~ 0.05. These numbers and the density-scaling result (observed slope Gamma=0.393+/-0.079 vs Newtonian CMC baseline Gamma_N=0.748+/-0.039, 4.1-sigma rejection; TEP prediction Gamma_TEP=2*Gamma_N-1~0.50 matches at 0.9-sigma) are Paper 10's (TEP-COS) own real empirical+theoretical result from 550 real MSPs (199 globular-cluster + 351 field) and a 20-cluster/18,813-MSP N-body CMC null baseline -- NOT computed by step_05. Paper 0's step_05 (derive_msp) instead runs its own separate, illustrative order-of-magnitude estimate (generic placeholder cluster, Gamma_MSP=1 direct proper-time projection) that does not reproduce or connect to Paper 10's slope test; it should not be read as independent confirmation of it. Corrected 2026-09-16 -- an earlier version of this description incorrectly attributed Paper 10's real result to step_05.

#### tilde_kappa_MSP

Full response coefficient from step_44_kappa_msp_prior.json, obtained by inverting the 0.63 dex raw excess with real cluster parameters: κ̃_MSP = (2.9 ± 4.5) × 10⁴. This is the observable response; the pre-screening geometric amplitude is ~10⁶ (see kappa_MSP). Suitable for cross-probe comparison with kappa_Cep. Effective screening factor S_MSP ~ 0.027 relative to unsuppressed ~10⁶.

#### kappa_WB

Wide-binary channel. Only the transition radius is derived (step_05): R_s ~ 2609 AU from the cross-scale acceleration condition (M_binary, M_gal, r_gal, beta_A only, no wide-binary input), matching observed 2646 +/- 182 AU within 1.4%. kappa_WB itself is NOT derived: an earlier version asserted a screening fraction S(R_s)=0.5 without derivation to produce a velocity excess alpha=sqrt(2)-1~0.414 (~13% from observed 0.366), but the radial ODE's own computed screening at R_s is S_sigma~0.994 (essentially unscreened), which the same formula would instead put at alpha~0.73 -- the assumed 0.5 was doing the work, not the theory. Paper 13 (TEP-WB) does not independently derive alpha_sat=0.366 either: it reports a mass-dependent self-screening fit alpha_sat(M)=alpha_0*exp(-M/M_screen) (alpha_0=1.74+/-0.26, M_screen=0.46+/-0.04 Msun) fit to the same wide-binary sample used to validate it. No validated Gamma_WB exists; kappa_WB is left open pending a real derivation.

**Derivation status:** partial; transition-scale structure only, amplitude not derived

#### kappa_LLR

LLR/flyby channel response coefficient. Derived (step_05): kappa_LLR = |beta_A| S_LLR Gamma_LLR, with Gamma_LLR = 1 (orbital range projector; the fitted Keplerian monopole is absorbed) and S_LLR = S_Sigma(g_lunar) = 1.59e-14 from the corpus screening operator S_Sigma(g) = [1+(g/g_t)^2]^-1, g_t = cH_0/(2|beta_A|) (master-action closure, R11/step_27; equivalently the pairwise form [1+(R_s/s)^4]^-1). Small but nonzero at lunar distance; far below the current LLR residual floor. Earlier versions used a Yukawa propagator exp(-r/R_T)(1+r/R_T) treating the density-transition radius R_T ~ 4146 km as a Compton range — an inconsistent massive-scalar ansatz superseded by the derived P(X) screening operator. Flyby: conformal channel screened to delta_v/v ~ 1.7e-30 at the surface field; any real flyby residual requires the disformal sector (step_08). Corrected 2026-09-25.

#### Gamma_X

Geometric/kinematic projector for channel X (steps 04-05): the dimensionless channel-geometry factor in kappa_X = |beta_A| S_X Gamma_X. Derived values: Gamma_GNSS = |beta_A| (two-point clock covariance is quadratic in the coupling; the fluctuation Green's-function spectrum and correlation length lambda_T ~ R_T enter the channel functional F_X, not the projector), Gamma_Cep = 5/ln10 ~ 2.17, Gamma_MSP = 1 (proper-time projection), Gamma_gal = Gamma_Cep (1+z)^(beta_A Delta u) ~ 2.17, Gamma_LLR = 1 with the S_Sigma(g) suppression carried by S_X; Gamma_WB partial (alpha_sat open). Earlier versions conflated the unit-source Green's-function normalization with the projector (Gamma_GNSS ~ 1e-84, Gamma_LLR ~ 5e-39 via a Yukawa-on-R_T ansatz); corrected 2026-09-25 to the dimensionless channel-geometry convention.

#### strong_field_status

Strong-field analysis (step_06): perturbative Einstein-conformal-scalar sector is strongly hyperbolic; GW170817 bounds the observed propagation-sector tensor-cone split (~few x 10^-15) but does not directly bound the matter-dependent scalar principal coefficients (KMZ 2012). TEP EOM remain second-order (Horndeski-class, no Ostrogradsky); scalar-Gauss-Bonnet is likewise second-order, though its characteristics can fail or degenerate on nontrivial backgrounds (Papallo-Reall 2017). Necessary conditions: Z_t > 0, Z_s > 0, metric nondegeneracy — not sufficient: A=1, B=1, q=0, rho=3, p=2 is regular but elliptic (Z_t=4, Z_s=-1, c_s^2=-1/4). Full coupled strong-field proof open.

#### quantum_status

Quantum formalization (step_07): 6 kinematics items established from proper-time postulate (Schrodinger equation, unitarity, species sensitivity, phase accumulation, decoherence requires trace, Dirac conformal rescaling). 6 dynamics items remain candidate (path integral, spin-1/2, antiparticles, Born rule, renormalizability, Standard Model) requiring additional structure in Papers 24-25.

#### epsilon_T

Amplitude of large-scale temporal shear in cosmological evolution. Under investigation; value not yet firmly established.

#### phi0

Phenomenological amplitude parameter in redshift-dependent ansatz. Not the fundamental scalar field value.

#### n_phi

Exponent in A(z) phenomenological template. Under investigation.

#### beta_spin

Coefficient in the phenomenological screening ansatz A(ρ) = exp(-β tanh(ρ/ρ_c)). NOT the fundamental conformal coupling β_A = -1.0. Determined by requiring the screening transition span the observed density range.

#### beta_cassini_max

Cassini PPN bound (Bertotti et al. 2003), |γ_PPN − 1| < 2.3×10⁻⁵, expressed as the bound on an effective unscreened DEF-normalized coupling |α_DEF| < 3.4×10⁻³ (α_DEF² = 2β_A²). For the frozen β_A = −1 the binding constraint is on the screened Solar source charge, S_Σ(1 AU) ≲ 5.8×10⁻⁶, via the corrected screened-source mapping γ_PPN−1 = −4β_A²S_Σ/(1+2β_A²S_Σ) (linear in S_Σ because the photon probe is unscreened; Burrage & Sakstein 2018, Eq. 3.31). Screening is supplied by the density-dependent minimum of the field equation (environmental operator S_Σ(E); the quartic completion satisfies the bound for μ_0 ≳ 10^10 while retaining the GNSS-scale response up to μ_0 ~ 10^11). The earlier phenomenological screening function f(ρ) = [1+(ρ/ρ_half)^n]⁻¹ is retired.

#### B

Disformal coupling function B(phi). GW170817 constrains the path-integral combination B(phi)(dphi)^2 along observed late-time astrophysical paths: |c_gamma - c_g|/c < few x 1e-15. This does NOT require B to vanish identically in all regimes. B(phi) = 0 in the conformal limit (holonomy vanishes). The nonzero form B(phi) = B0 * |phi|^2 / (1 + |phi|^2) * exp(-phi^4 / (2*sigma_B^4)) (quartic-Gaussian ultra-damped shear bump) is implemented in Paper 28 bh_common.py, with B0, n_B, sigma_B as shape parameters. B activates the disformal shear near strong-field regions while decaying super-polynomially in the deep core. On spacelike-gradient configurations the signature margin is wide (Paper 0 step 52); in the timelike-gradient drift sector the bare envelope is ambient-lapse-supercritical at absorber epochs for phenomenological normalizations (realized cap B0 <~ 0.03 at z ~ 2.5-3.6, phibar ~ ln(1+z), du/dt ~ H_T(z); Paper 29 gate10c), so ambient closure through the environmental gate is lapse-mandated, not merely phenomenological. The path-dependent GW constraint is satisfied because B(phi)(dphi)^2 is negligible along weak-field late-time astrophysical paths. Admissibility: the null-cone condition selects B >= 0 as the unconditional branch — the signature condition B(dphi)^2 > -A^2 then holds wherever the matter metric is non-degenerate, the disformal perfect-fluid principal coefficients obey Z_t, Z_s >= 1, and the matter cone lies inside or on the gravitational cone. B < 0 is not excluded a priori but carries a pointwise ledger on each realized configuration — |B|(dphi)^2 < A^2 (signature), Z_t, Z_s > 0 (hyperbolicity), no matter-cone excursions beyond the gravitational cone (causality) — and sign-indefinite B(phi) inherits the same condition along the realized profile. Configurations computed to date fail it on the negative branch (the step_13 volume-balance reconstruction drives Q = 1+(B/A^2)(dphi)^2 -> -1, a matter-metric degeneracy with superluminal longitudinal photons; Paper 0 App. E R5). Eligibility vetoes (corpus registry, issue 29-4): any candidate form must satisfy (i) Cassini and GW170817 path bounds along realized profiles, (ii) ambient-lapse sub-criticality under the self-consistent map 1+z = N(0)/N(phi_bar) with N = A*sqrt(1-u), u = B(dphi)^2/A^2 — the veto that excludes the descent-tracking form (Paper 29 gate10: u(phi_bar) = 1 - 8e-5 ambient at all three absorber epochs, presenting z_obs = 2.5-3.6 at 1+z ~ 390-595; required phi_bar = -3.35 to -3.46 violates the sign convention), and (iii) three-sightline amplitude transfer without per-sightline retuning. Candidate status: Paper-0/28 quartic-Gaussian envelope = live corpus candidate (fails BBN closure under the imported LCDM rate — B0 dispersion factor ~3-4 and Cassini-excluded normalizations — but not ruled out as a form); power-law B0*phi^n = fails transfer (n=2 Cassini-excluded unscreened, n=5 loses Lorentzian signature above z~3.4); descent-tracking B0*phi^4/(phi_c^4+phi^4)*e^{-lam phi} = EXCLUDED by veto (ii) notwithstanding inverse-root concordance to 1.0001; environment-modulated B_eff = B(phi)*G(X_local), closed in the ambient and saturating at G ~ 1e-3 inside wells (required gate slope dG/d ln(dphi) ~ (4-7)e-4 per e-fold, consistent across the three BBN sightlines within a factor 1.6; Paper 29 gate10b) = the only admissible form for the absorber amplitude channel, a D2(phi,X)-class generalisation.

#### H0

TEP-corrected local Hubble constant from per-host cz/d after removing ~1 km/s/Mpc host-potential clock distortion. Consistent with the TEP CMB (66.70 ± 0.58, Paper 26) at 0.03σ. The SH0ES ladder value 73.04 is the Hubble-flow temporal shear (cosmic-web path integral), not a correction to be closed.

#### S8

S_8 = σ_8 √(Ω_m/0.3) from TEP cosmological fit. Under investigation in Paper 26.

#### V_scalar

Rule 12/18 anchor: one universal scalar sector, fixed once. Non-canonical V(phi) realizations present in live pipelines are named benchmarks bounding methodology only, not the canonical potential: the V_rec(u) = V_0 - (rho_bar/2) e^{-u} descent branch (Paper 0 sec.8 reconstruction, BBN Gate 10), and the BBN step_04g massive-KG 1/2 m^2 phi^2 plus V=0 benchmarks. (issues.md v-phi-benchmark-proliferation, 29-4 one-action registry)

#### Lambda_X

Scale entering the non-perturbative kinetic completion P(X,phi) = X - V + X|X|/Lambda_X^4, i.e. P_{,X} = 1 + 2|X|/Lambda_X^4 on the spacelike branch. Fixed by Lambda_X^4 = M_Pl^2 H_0^2, Lambda_X = sqrt(M_Pl H_0) = 1.9068 meV with H_0 = 70 km/s/Mpc the measured drift (constants.py: LAMBDA_KINETIC_MEV). Anchored to H_0, not fitted. Sets the shear threshold g_t = cH_0/(2|beta_A|) = 3.4e-10 m/s^2 and the single-body screening operator S_Sigma(g) = [1+(g/g_t)^2]^-1. Distinct from the potential-sector scale: Lambda_X governs the gradient/shear sector, not the amplitude sector. (TEP-UCD step_14 scale audit.)

#### lambda_quartic

Dimensionless coupling in the matter-hosting weak-field branch V(u) = (lambda/4) M_Pl^4 u^4 of the master potential. The operative branch lambda_Cassini = 7.526e-66 = 10^5 lambda_ref is selected by the Cassini screening bound (the lambda_ref = 7.526e-71 point fails the corrected Cassini evaluation; Paper 6 App. C, Paper 7). The quartic exponent is itself selected, not free: the corpus amplitude law S_A ∝ rho^{1/3} = rho^{1/(m-1)} forces m = 4 (TEP-UCD step_14 amplitude_law_check). Compact-object bound: pulsar interiors remaining below u_s require lambda >~ 1.11e-84 (LAMBDA_QUARTIC_NS_BOUND); the Cassini branch passes by ~6.8e18, while the retired rho_T-normalized reading lambda = rho_T/M_Pl^4 = 2.45e-90 is excluded by compact objects, Compton resolution and Solar-System screening. Symbol collision note: several papers write it bare 'lambda'; where the same file also uses the GNSS covariance length the coupling is disambiguated (Paper 6 writes lambda_T for the length).

#### Lambda_V

Lambda_V = lambda_quartic^{1/4} M_Pl — the energy scale of the potential normalization expressed in units of mass: 127.5 GeV on the Cassini branch, 7.2 GeV on the reference branch (constants.py: LAMBDA_V_GEV / LAMBDA_V_REF_GEV). NOT an independent scale — derived from lambda_quartic. Its role is to prevent the two-scale conflation: the identification rho_T = Lambda_V^4 (equivalently Lambda ~ 96 keV) presupposes u_min(rho_T) ~ 1 and is excluded (step_14 excluded_normalization: 3e24 below the Cassini coupling, u_min(NS) ~ 21 > u_s, lambda_c(rho_T) ~ 3.4e7 km ≫ any terrestrial body).

#### rho_sat

Order-of-magnitude density at which the quartic equilibrium reaches u_min ~ 1: rho_sat ~ lambda_quartic M_Pl^4 (step_14: 6.14e25 g/cm^3; with the e^{-u} source factor the strict u_min = 1 crossing sits at e lambda M_Pl^4 ~ 1.7e26). This is the temporal-well interior scale (Paper 28), where exponential source suppression and the master-potential knee take over — NOT terrestrial saturation. Terrestrial screening saturation is the geometric/Compton-resolution crossover lambda_c(rho) ~ R_T(M), with lambda_c(rho_T) ~ 2.5e3 km ~ 0.6 R_T(Earth). The two scales form one density ladder on a single branch: Compton crossover at ~4 g/cm^3, perturbative quartic equilibria through pulsar densities (u_min ~ 1.5e-4), self-quenching only at rho ~ 1e26.

#### u_min

Interior minimum of the coupled field equation for the quartic branch: lambda M_Pl^4 u_min^3 = rho e^{-u_min}, i.e. u_min = (rho/(lambda M_Pl^4))^{1/3} for u_min << 1 (core/scalar_field.py: equilibrium_u). Values on the Cassini branch (step_14): 4.48e-9 at Earth's mean density, 6.88e-9 at rho_T, 1.48e-4 at neutron-star-core density. The amplitude response is exact quartic bookkeeping: S_A(rho) = min[1, u_min(rho)/u_min(rho_T)] = min[1, (rho/rho_T)^{1/3}]. Complementary to the kinetic sector: P_{,X} suppresses the gradient excursion while the potential fixes the interior amplitude (Earth interior u(0) ~ 2e-15 on the suppressed profile vs u_min ~ 4.5e-9).

#### u_s

Nominal knee of the master family V(u) = V_matter(u) e^{-(u/u_s)^4} + V_0 e^{-(u_s/u)^4}: the dimensionless field amplitude at which the matter-branch damping engages. Declared nominal value u_s = 10 (master-potential ansatz; closure probes scanned 10--50). Bounds: pulsar interiors must satisfy u_min(rho_NS) < u_s, which yields the quartic-normalization lower bound LAMBDA_QUARTIC_NS_BOUND = 1.11e-84.
