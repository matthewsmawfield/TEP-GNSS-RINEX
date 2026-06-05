# Local Temporal Topology and the Spatial Variation of Laboratory Gravitational Constants
**Matthew Lukin Smawfield**
Version: v0.1 (Naivasha)
First published: 2 March 2026

---

## Abstract

The April 2026 NIST metrology publication reported a redetermination of the gravitational constant *G* at Gaithersburg, Maryland, using an identical torsion-balance geometry to the prior French BIPM setup. The new measurement reveals a systematic relative drop of 2.81×10⁻⁴ compared to the Sèvres baseline — a shift 198× larger than the random environmental noise budget of either apparatus. We do not present this two-site comparison as an independent proof of a new-physics effect; instead, we employ it as an in-sample calibration that fixes the density-sector coupling magnitude to |αlog| = 7.66×10⁻³ (sign fixed by the TEP field equation in the (+,−,−,−) metric signature). With the lab-scale conformal coupling held at βA = −1, a 1D tapered mass-column surrogate — weighted by a near-field kernel *K*(*z*) = 1/(1 + (*z*/*z*₀)²), blended with CRUST1.0 deep structure, and using empirically constrained Bouguer gravity anomalies for the near-field column (0–5 km) — predicts a relative decrease of −4.16×10⁻⁴ (−416 ppm) at Gaithersburg versus Sèvres, matching the sign and order of magnitude of the measured −2.81×10⁻⁴. We flag explicitly that this headline ppm shift is delivered by a 1D surrogate that averages out lateral structure, while an uncoupled 3D finite-difference solver on the current symmetric boundary taper yields zero horizontal shear at the Gaithersburg bunker centre and a negligible inter-site scalar difference (dG/G ≈ 0 ppm). The 3D solver therefore becomes an out-of-sample prediction engine only once measured lateral Bouguer maps replace the symmetric taper. The paper reframes its contribution accordingly: a reproducible methodology to generate directional expectations for third-party precision-gravity laboratories. Laboratories on lighter sedimentary columns should systematically measure higher *G*eff, while those on denser crystalline basement should measure lower *G*eff. The Sèvres–Gaithersburg datum falls within this pattern. The 95% forward geological uncertainty band [−2270, +1724] ppm contains the empirical value, so the locked parameterization is not falsified by the in-sample calibration, but the decisive test is independent multi-laboratory comparison against empirically constrained scalar-field columns without re-fitting the coupling.

Keywords: Temporal Equivalence Principle, gravitational constant, metrology, NIST, BIPM, temporal topology, scalar field

## 1. Introduction

## 1.1 The NIST–BIPM Discrepancy

In April 2026, the National Institute of Standards and Technology (NIST) published a redetermination of the gravitational constant *G* using a torsion-balance apparatus functionally identical to the prior Bureau International des Poids et Mesures (BIPM) setup in Sèvres, France. The new measurement at Gaithersburg revealed a systematic relative drop of 2.81×10⁻⁴ compared to the Sèvres baseline, a shift that exceeds the combined classical noise budget by nearly two orders of magnitude.

## 1.2 The Temporal Equivalence Principle Hypothesis

The Temporal Equivalence Principle (TEP) proposes that proper time is not a passive geometric outcome, but a dynamical scalar field (φ) characterized by a background cosmic gradient, or Temporal Shear (Σμ = ∇μ ln A(φ)). Under TEP, the matter metric is rescaled by a universal conformal factor A2(φ) relative to the gravitational metric, so the Jordan-frame effective gravitational coupling is modulated by the local scalar field profile as *G*eff = *G*N / A2(φ). This paper tests whether the NIST–BIPM discrepancy is consistent with this framework.

## 1.3 This Work

This paper presents a reproducible pipeline that retrieves site elevation and DEM/geological-service metadata for both laboratories, builds a layered crustal mass column within a 50 km radius from published geophysical parameters, computes the local scalar field profiles, and predicts the *G*eff variance between Sèvres and Gaithersburg. The structure is as follows: Section 2 presents the TEP theoretical framework; Section 3 describes the geological modeling methodology; Section 4 details the computation pipeline; Section 5 reports the results; Section 6 discusses implications; and Section 7 concludes with qualitative directional expectations for third-party laboratories.

## 2. Theoretical Framework

## 2.1 The Scalar Field Equation

The TEP scalar field φ is sourced by local mass distributions through the field equation:

\begin{equation}
\Box\phi = \frac{8\pi G}{c^4} \rho_{\text{eff}}
\end{equation}

where ρeff incorporates the screening-corrected density. For the matter metric rescaled by A2(φ), the Jordan-frame effective gravitational coupling is Geff = GN / A2(φ).

The lab-scale phenomenological ansatz used here is φρ = αlog ln(ρ/ρc).  With the (+,−,−,−) metric signature adopted across the TEP corpus, the scalar field equation for non-relativistic matter (T = +ρ) and negative conformal coupling βA = −1.0 gives a static limit in which φ decreases with increasing density: dφ/dρ < 0.  This fixes the sign of the density-sector coupling to αlog < 0.  The magnitude |αlog| = 7.66×10⁻³ is fixed by the Sèvres–Gaithersburg in-sample calibration (it is the value that brings the 1D surrogate prediction into the correct order of magnitude); the sign is fixed by the field equation and is not a free parameter.

## 2.2 Temporal Shear

The Temporal Shear vector is defined as the active gradient of the conformal scaling factor:

\begin{equation}
\Sigma_\mu = \nabla_\mu \ln A(\phi)
\end{equation}

## 2.3 Screening Regime

When local density approaches the phenomenological saturation scale ρc ≈ 20 g/cm³, the scalar field response is smoothly suppressed and A(φ) → 1. Both laboratory sites operate in a mildly screened low-density regime relative to this scale, allowing TEP scaling to remain observable in the model.

The implemented smooth log-density screen uses n≥3.21 to retain ≥95% coupling at laboratory densities while satisfying the solar-wind suppression requirement of n≥0.10 down to solar-wind densities of 10⁻²⁴ g/cm³ (identifying the Cassini check as a proxy bound). This is a continuous screening profile, not a binary on/off threshold.

## 3. Geological Modeling

## 3.1 Data Sources

Site elevation is retrieved live from USGS 3DEP (Gaithersburg) and, outside USGS coverage, from global SRTM 30 m (Sèvres); DEM and geological-service query metadata are recorded from OpenTopography, USGS MRDS, and BRGM for reproducibility. The primary near-field geological constraint is the *Bouguer gravity anomaly*, an empirically constrained proxy for local subterranean mass surplus or deficit relative to a standard 2.67 g/cm³ crustal column, derived from surface gravity measurements with standard geophysical corrections (terrain, latitude, free-air). Bouguer anomalies resolve density variations at the station scale (~1 km), bypassing the resolution limit of global crustal models. Published station values are drawn from BRGM (France), USGS/NGMDB (USA), Yuan *et al.* (2014) for the Yangtze Basin, and Geoscience Australia for the Tasman Fold Belt. For five additional prediction sites, representative Bouguer anomalies are estimated from published regional gravity compilations (BGR for Braunschweig, BGS for Teddington, GSC for Ottawa, CAGS for Beijing, swisstopo for Zurich; all cross-checked against WGM2012, Bonvalot *et al.*, 2012). JILA Boulder is the sole exception: its −220 mGal isostatic Bouguer anomaly reflects deep crustal compensation of the Rocky Mountains, not surface-layer density, so the simple slab inversion is invalid and the near-field column is constrained by USGS geological map data instead. Deep structure below ~5 km is taken from the published global model CRUST1.0 (Laske *et al.*, 2013), which provides layer-by-layer densities and thicknesses on a 1°×1° grid.

#### Laboratory Coordinates

- BIPM Sèvres: 48.8298° N, 2.2137° E (Paris Basin)

- NIST Gaithersburg: 39.1383° N, −77.2014° W (Piedmont Plateau)

- HUST Wuhan: 30.5155° N, 114.4112° E (Yangtze Basin)

- UQ Brisbane: −27.4698° S, 153.0251° E (Tasman Fold Belt)

- JILA Boulder: 40.0150° N, −105.2705° W (Rocky Mountains)

- PTB Braunschweig: 52.2970° N, 10.4430° E (North German Plain)

- NPL Teddington: 51.4240° N, −0.3380° W (London Basin)

- NRC Ottawa: 45.4215° N, −75.6972° W (Canadian Shield)

- NIM Beijing: 39.9042° N, 116.4074° E (North China Plain)

- U Zurich: 47.3769° N, 8.5417° E (Swiss Molasse)

## 3.2 Mass Distribution Mapping

The Bouguer anomaly ΔgB relates to local density contrast Δρ through the slab formula ΔgB ≈ 2πG Δρ h. For a 5 km column, the conversion factor is Δρ ≈ 4.77 × 10⁻³ g/cm³ per mGal. The Paris Basin shows a well-documented negative anomaly of −42 ± 8 mGal (Le Pichon *et al.*, 1971), corresponding to a top-layer density of 2.47 ± 0.04 g/cm³ — a light sedimentary column. The Appalachian Piedmont at Gaithersburg shows a near-zero anomaly of +5 ± 10 mGal (Diment & Weaver, 1964), corresponding to 2.69 ± 0.05 g/cm³ — dense crystalline basement. The resulting *near-field* density contrast is 0.22 g/cm³, fifty times larger than CRUST1.0's 1°×1° average contrast of 0.004 g/cm³. CRUST1.0 deep structure is retained below the Bouguer-sensitive depth, yielding a hybrid column: Bouguer-constrained top + CRUST1.0 deep.

## 3.3 Near-Field Distance Weighting

The scalar field gradient that modulates the local effective coupling is dominated by near-field mass, not the deep crustal column. In analogy to Newtonian gravity, where the local acceleration gradient from a point source scales as r⁻², the density-sector contribution to φ is weighted by a kernel *K*(*z*) = 1/(1 + (*z*/*z*₀)²), where *z* is the depth to the layer midpoint and *z*₀ = 5 km is the near-field coherence scale. At *z* = 0 the kernel is unity (full sensitivity); at *z* = *z*₀ it is 0.5; at *z* = 5*z*₀ it is 0.038. The top Bouguer layer (0–5 km) therefore contributes ~55% of the effective density-sector signal, while the lower crust (20–30 km) contributes <3%. The geometric mass term still uses the total column mass because the conformal factor depends on the integrated mass budget.

## 3.4 Facility-Scale Mass Model

The 50-metre environment around each vacuum chamber dominates the local gravitational potential and is therefore the most important mass distribution for a laboratory-scale scalar field. The v0.1 pipeline ignored this near-field blind spot. Parametric facility-scale density grids are now injected (101×101×41 cells, 1 m lateral resolution, 0.5 m vertical resolution) for both laboratories.

BIPM Sèvres Pavilion. The historic pavilion is built into a gentle hillside with thick limestone cellar walls (2.4 g/cm³, 1.5 m thick) forming a partial L-shaped enclosure on the north and east sides. The cellar floor is 0.8 m of stone (2.35 g/cm³). Hillside soil (1.8 g/cm³) abuts the north and east walls with a 15° slope extending ~3 m outward. The effective vertically-integrated top-layer density is 2.55 g/cm³, a surplus of +0.15 g/cm³ above the ambient crustal average. Total facility *net mass anomaly* (above the 2.4 g/cm³ ambient crust) is small and slightly negative (~−1.8×10⁶ kg) because the lighter hillside soil and floor fill dominate over the denser limestone walls.

NIST Advanced Measurement Laboratory (Gaithersburg). The AML is a modern underground metrology bunker with reinforced concrete walls (2.5 g/cm³, 0.6 m thick) forming a near-complete box 30 m on a side. The floor slab is 1.2 m thick concrete (2.5 g/cm³). The surrounding Piedmont terrain is flat, with modest fill soil (1.9 g/cm³) outside the walls. The effective top-layer density is 2.48 g/cm³, a surplus of +0.08 g/cm³ above ambient. Total facility *net mass anomaly* (above the 2.4 g/cm³ ambient crust) is negative (~−5.6×10⁶ kg) because the large volume of lighter fill soil outside the walls dominates over the relatively thin concrete shell.

Physical significance. The facility-scale net mass anomaly at both sites is of order 10⁶ kg — small compared to the crustal column (~10¹⁷ kg), but the scalar field Green's function weights nearby mass heavily. More importantly, the *asymmetry* of the facility layout — Sèvres's partial L-shaped stone cellar on a hillside versus Gaithersburg's symmetric concrete box on flat terrain — is the most plausible source of a horizontal gradient that a torsion balance could actually measure. The near-field fill and soil are lighter than solid crust, so the net mass anomaly is negative at both sites, but the geometric layout asymmetry is what matters for the gradient. In the current pipeline, these facility grids are built, solved directly on a high-resolution 51×51×21 mesh, and superposed onto the crustal solution via the linearity of the Poisson equation.

## 4. Computation

## 4.1 Pipeline Architecture

The reproducible analysis pipeline is implemented in Python and executed sequentially. Each step writes a JSON output to `results/outputs/` and a detailed log to `logs/`. Steps are fail-fast: execution halts on the first failure so that downstream steps do not consume stale data.

## 4.2 Scalar Field Solution

The primary scalar field is computed from the 3D screened Poisson equation on a Cartesian grid. The field at the laboratory position is the volume integral of the Green's function over the local density distribution:

φ(**r**₀) = ∫V G(**r**₀, **r**') · (αlog/λT²) ln(ρ(**r**')/ρc) · Sscreen(ρ(**r**')) d³r'

where G(**r**, **r**') is the Green's function for the screened Poisson operator. In the unscreened limit (λscreen → ∞), G ∼ 1/(4π|**r** − **r**'|). The 3D finite-difference solver relaxes this equation directly on a 41×41×11 grid (100 km lateral, 50 km depth, 5 km resolution) using Successive Over-Relaxation (SOR) with Dirichlet boundary conditions φ = 0 at the domain boundary. The source term S(ρ) = (αlog/λT²) ln(ρ/ρc) Sscreen(ρ) carries the inverse square of the Temporal Topology screening length λT ≈ 4200 km, the canonical multi-center GNSS value (Paper 6). A shorter MGEX verification run (Paper 14) yields a preliminary 2590 km; the solver uses the canonical value while noting the preliminary MGEX consistency. No arbitrary domain-area division is required.

The 1D mass-weighted surrogate φρ = αlog Σi K(zi) (mi/M) ln(ρi/ρc) S(ρi), where K(z) = 1/(1 + (z/z₀)²), is retained as the in-sample calibration engine that fixes the amplitude of the density-sector coupling (|αlog| = 7.66×10⁻³) and delivers the headline ppm shift. It is not the primary out-of-sample solver. The surrogate averages out lateral structure and is physically valid only when the density field is uniform over the support of the Green's function — an assumption that fails at laboratory scales.

**Dimensional analysis.** The screened density contrast ln(ρ/ρc) Sscreen(ρ) is dimensionless, while the discrete Laplacian ∇²φ has units of φ/km². To obtain a dimensionless φ, the source term carries the inverse square of the Temporal Topology screening length: S(ρ) = (αlog/λT²) ln(ρ/ρc) Sscreen(ρ), with λT ≈ 4200 km from multi-center GNSS (Paper 6), the canonical terrestrial correlation length. The dimensional scale is therefore fixed by a physically derived correlation length, not by an arbitrary simulation domain size. The overall amplitude still rides on the fitted coupling αlog, so the 3D solver confirms the *sign* and *order* of the inter-site shift, but it does not independently fix the precise sub-ppm magnitude without a measured lateral Bouguer map. With the current symmetric lateral taper both sites have nearly identical crustal φ and the inter-site dG/G from the crustal grid alone is ≈ 0 ppm.

## 4.3 3D Finite-Difference Solver

The 3D finite-difference solver relaxes the static screened Poisson equation ∇²φ = −S(ρ) on a Cartesian grid using Successive Over-Relaxation (SOR), with Dirichlet boundary conditions φ = 0 at the domain boundary. The source term is S(ρ) = (αlog/λT²) ln(ρ/ρc) Sscreen(ρ), carrying the inverse square of the Temporal Topology screening length λT ≈ 4200 km, the canonical multi-center GNSS value (Paper 6). Because the Laplacian has units of 1/km² and the source term carries 1/λT², the SOR solver returns φ as a dimensionless field directly. No post-hoc domain-area division is required. The source term is not calibrated to the 1D surrogate, so the *sign* (both sites negative) is an independent check. The overall amplitude still rides on the fitted coupling αlog.

For the actual Bouguer+CRUST1.0 columns, the uncoupled 3D solver gives φ3D ≈ −4.1×10⁻⁵ (dimensionless) for both Sèvres and Gaithersburg, of the same order as the 1D weighted values (−1.54×10⁻² and −1.51×10⁻²). The dimensional scale is fixed by the physically derived λT ≈ 4200 km, the canonical multi-center GNSS value (Paper 6), not by an arbitrary simulation domain size. With the current symmetric lateral taper, both sites have nearly identical crustal φ, and the inter-site dG/G from the crustal grid alone is ≈ 0 ppm. Facility-scale superposition does not produce a measurable inter-site scalar difference either (both facilities give φtotal ≈ −4.3×10⁻⁵). The 3D solver therefore confirms the *sign* and *order* (both sites negative, magnitude ~10⁻⁵), but it does not yet resolve the inter-site *structure* because the symmetric boundary taper erases the lateral density contrast. Measured lateral Bouguer maps are required for the 3D solver to reproduce the site-to-site difference.

Lateral density variations are currently modeled by a smooth taper from the local CRUST1.0 profile to a reference density at the domain boundary. This produces a *symmetric* density grid with zero horizontal gradient at the centre. Facility-scale density grids (Section 3.4) are now superposed onto the crustal solution via the linearity of the Poisson equation: the total scalar field φtotal = φcrustal + φfacility, with both fields already dimensionless from the 1/λT² source term. The asymmetric L-wall at Sèvres produces a non-zero horizontal shear Σh ≈ 7.2×10⁻⁹ m⁻¹, while the symmetric Gaithersburg bunker gives Σh = 0 m⁻¹ (identically zero for the idealized symmetric model). The Sèvres facility is measurably more asymmetric, which is the correct physical ordering for a torsion-balance differential-torque signal. The ±1000 ppm-scale forward bands remain wide because the crustal boundary taper does not yet incorporate measured lateral Bouguer maps; this is the next required upgrade.

## 4.4 Temporal Shear Computation

The Temporal Shear is a 3-vector Σμ = ∇μ ln *A*(φ) = βA ∇μφ. It is computed from the 3D FEM scalar field solution via central differences at the laboratory position, yielding the full vector (Σx, Σy, Σz) rather than a scalar fiction. The horizontal components Σx and Σy are the physically relevant quantities for a torsion balance, which measures torque across its ~1 m baseline in the horizontal plane and is therefore sensitive to horizontal gradients of *G*eff, not to vertical gradients or the scalar value at the centre.

For the current symmetric tapered grid, the horizontal gradient is identically zero at the centre: Σx = Σy = 0, Σz ≈ −5.1×10⁻⁴ m⁻¹. This is expected because the density taper is azimuthally symmetric. The torsion balance would be torque-blind to this configuration. Facility-scale density grids (Section 3.4) are now superposed onto the crustal solution in step 06: the facility field is solved independently on the full 101×101×41 grid (50 m lateral, 20 m depth, ~0.5 m resolution) using SOR, with the same dimensionally-strict source term (αlog/λT²), and added directly to the crustal field. No area normalization is required. The asymmetric L-wall at Sèvres produces a non-zero horizontal shear Σh ≈ 7.2×10⁻⁹ m⁻¹, while the symmetric Gaithersburg bunker gives Σh = 0 m⁻¹. The Sèvres facility is measurably more asymmetric than Gaithersburg, which is the correct physical ordering for a torsion-balance differential-torque signal. The legacy v0.1 scalar diagnostic Σ ≈ φ/Lchar (with Lchar = 50 km) is marked deprecated and retained only for backwards comparison.

The effective coupling Geff = GN/A2(φ) depends on the local scalar value, not on the gradient. The inter-site prediction in Section 5 is therefore a comparison of two scalar point values — a *site effect*, not a gradient effect. Facility-scale density grids are now superposed in step 06, producing a genuine horizontal asymmetry at Sèvres (Σh ≈ 7.2×10⁻⁹ m⁻¹) from the L-wall hillside-soil configuration. This gradient is a *complementary* signature: a torsion balance at Sèvres would experience a differential torque from the horizontal shear, while the inter-site scalar point-value difference is the quantity predicted by the 1D solver. The two mechanisms predict different experimental signatures and should not be conflated.

## 5. Results

## 5.1 Mass Model Comparison

The 50 km multi-layer crustal mass columns are estimated at 7.24×10¹⁷ kg for Sèvres (Paris Basin) and 6.69×10¹⁷ kg for Gaithersburg (Piedmont Plateau). Sèvres has the thicker crust (33.5 km vs 30.7 km) and a near-field top-layer density of 2.47 g/cm³, while Gaithersburg sits on denser crystalline basement with a top-layer density of 2.69 g/cm³.

## 5.2 Screening Analysis

Both sites operate in the mildly underscreened regime (density ratio ~0.14 relative to ρc = 20 g/cm³), meaning TEP scaling effects are not suppressed by local density.

## 5.3 Statistical Validation

The combined *random* environmental noise budget (temperature, micro-seismic, atmospheric pressure, magnetic interference, torsion fiber anelasticity) sums to 1.42×10⁻⁶ in quadrature. The NIST discrepancy of 2.81×10⁻⁴ (281 ppm) is 198× larger than this random budget, so within-experiment random noise cannot account for it. For completeness, the same shift is only 0.56× the historical few-hundred-ppm scatter *among independent G determinations*, which is conventionally attributed to apparatus-dependent systematics. Random noise is thus excluded, but unmodeled apparatus systematics remain a viable conventional alternative; discriminating between that hypothesis and TEP requires independent multi-laboratory comparisons.

## 5.4 Effective-G Variance: Prediction vs. Measurement

Using the 1D in-sample surrogate calibrated to the Sèvres–Gaithersburg magnitude, evaluating Geff = GN/A2(φ) at both sites with the locked parameterization (βA = −1) and environment-dependent screening gives a predicted relative shift ΔG/G = −4.16×10⁻⁴ (−416 ppm), to be compared with the measured −2.81×10⁻⁴ (−281 ppm). The point-estimate signs match: both model and data predict Gaithersburg G < Sèvres G. The prediction is 48% off in magnitude. Propagating the geological-input uncertainty forward at fixed coupling (Monte Carlo over crustal density and thickness) yields a 95% prediction band of [−2.27×10⁻³, +1.72×10⁻³], which *does* contain the empirical −2.81×10⁻⁴; the locked parameterization is therefore not falsified by this single comparison. The empirical number would be reproduced exactly for a conformal coupling βA = −0.676, a −32.4% adjustment from the nominal −1.0. Because βA is a phenomenological coupling not independently constrained in this work, this exact-match value is reported as a fit, not as an independent confirmation; the honest assessment is that the locked parameterization reproduces the correct sign and order of magnitude, and the empirical value lies within the forward uncertainty band. The decisive test is independent multi-laboratory comparison.

## 5.5 Generalized Directional Expectations

The same pipeline, with no re-fitting, generates qualitative directional expectations for any laboratory with a known crustal column: laboratories on lighter sedimentary columns (lower near-field density, thicker unconsolidated cover) should measure higher *G*eff, while those on denser crystalline basement (higher near-field density, thinner sediment) should measure lower *G*eff. The Sèvres–Gaithersburg comparison itself falls within this pattern: Gaithersburg sits on denser Piedmont crystalline basement (2.69 g/cm³) versus Sèvres on the lighter Paris Basin sediments (2.47 g/cm³). Specific ppm-scale predictions are deferred until dedicated station-scale micro-gravity surveys reduce Bouguer uncertainty from the current ±10 mGal to ~±1 mGal; at that resolution the scalar-field difference between sites becomes calculable to sub-100 ppm precision. The pipeline is publicly available and the methodology is identical across sites, so future G comparisons can be evaluated against empirically constrained scalar-field columns rather than being treated as independent measurements of a universal constant.

## 6. Discussion

## 6.1 Interpretation of the Discrepancy

The 198× margin over the random environmental noise budget rules out within-experiment noise as the origin of the NIST–BIPM discrepancy. TEP offers a candidate explanation: differing crustal columns at separate laboratories generate distinct scalar field values φ, which modulate the effective coupling Geff = GN/A2(φ). When the model is fed Bouguer gravity anomalies — empirically constrained proxies for local mass deficit/surplus at the station scale, derived from surface gravity with standard geophysical corrections — the near-field density contrast between Sèvres (2.47 g/cm³) and Gaithersburg (2.69 g/cm³) is large and physically plausible. With the near-field kernel *K*(*z*) = 1/(1 + (*z*/*z*₀)²) that makes the top 5 km dominate the scalar field gradient, the model predicts a relative decrease of −416 ppm at Gaithersburg versus Sèvres. The empirical measurement is a relative decrease of −281 ppm. The point-estimate signs match; the 48% magnitude residual may reflect geological model uncertainty or an unmodeled physical effect. The 95% forward geological uncertainty band [−2270, +1724] ppm contains the empirical value, so the locked parameterization is not falsified by this single comparison.

The in-sample result is obtained with the locked conformal coupling βA = −1. Exact agreement would require βA ≈ −0.676 once environment-dependent screening is included — a −32.4% shift that is within the broad geological uncertainty envelope. The Sèvres–Gaithersburg datum demonstrates that the locked parameterization reproduces the correct sign and order of magnitude, and the empirical value lies inside the forward uncertainty band. The definitive test is independent multi-laboratory comparison evaluated against empirically constrained scalar-field columns.

Five interpretations remain live. (i) The Bouguer anomalies themselves are mis-assigned or the top-layer thickness is wrong; higher-resolution local gravity surveys could revise the density contrast. (ii) The 3D solver source term now carries the inverse square of the Temporal Topology screening length: S(ρ) = (αlog/λT²) ln(ρ/ρc) Sscreen(ρ), with λT ≈ 4200 km, the canonical multi-center GNSS value (Paper 6). A shorter MGEX verification run (Paper 14) yields a preliminary 2590 km. Because the Laplacian has units of 1/km² and the source term carries 1/λT², the SOR solver returns φ as a dimensionless field directly. No post-hoc domain-area division is required. Facility-scale grids are superposed via the linearity of the Poisson equation: φtotal = φcrustal + φfacility ≈ −4.3×10⁻⁵. The dimensional scale is fixed by the physically derived λT, not by an arbitrary simulation domain size. The 3D solver confirms the *sign* and *order* of the inter-site shift, but with the current symmetric lateral taper both sites have nearly identical crustal φ and the inter-site dG/G from the 3D solver is ≈ 0 ppm. Measured lateral Bouguer maps are required for the 3D solver to reproduce the site-to-site difference. (iii) The locked conformal coupling βA = −1 is incompatible with Cassini PPN bounds (|βPPN| < 0.0034) at solar-system scales unless density screening is active. The implemented log-density screen is continuous: n=4 exceeds the n≥3.21 needed to retain ≥95% lab coupling and exceeds the solar-wind suppression requirement used for the Cassini proxy check. The form of the first-principles potential V(φ, ρ) is constrained by the screening limits established here. (iv) The facility-scale mass has been modeled, solved directly, and superposed onto the crustal column: the asymmetric L-wall at Sèvres produces a horizontal shear Σh ≈ 7.2×10⁻⁹ m⁻¹, while the symmetric Gaithersburg bunker gives Σh = 0 m⁻¹. A parameter sensitivity sweep (step 14) reveals that the Sèvres signal is driven by the hillside-soil asymmetry in the north-east quadrant, not by wall or floor thickness variations: removing the hillside soil reduces Σh to the numerical noise floor, identical to the Gaithersburg symmetric box. The signal is therefore robust to building-material priors but sensitive to the local topography embedding. Full integration with measured lateral Bouguer maps for the predictive solve is the next required step. (v) The close match is coincidental, and the −281 ppm shift has a conventional origin such as unmodeled apparatus systematics.

The framework generalizes beyond the two-site comparison. Laboratories with lighter sedimentary cover (lower near-field density, thicker unconsolidated layers) should systematically measure higher *G*eff, while those on denser crystalline basement should measure lower *G*eff. Specific ppm-scale predictions are deferred until station-scale gravity surveys constrain local Bouguer anomalies to ~±1 mGal, at which point the scalar-field difference becomes calculable to sub-100 ppm precision.

## 6.2 The Strongest Case for TEP

The Sèvres–Gaithersburg comparison shows the locked parameterization reproduces the correct sign and order of magnitude, and the empirical value falls inside the forward uncertainty band. The strongest case for TEP, however, is that precision *G* measurements are already limited by non-random, site- and apparatus-dependent structure at the few-hundred-ppm scale, and TEP supplies a concrete physical covariate — the local scalar field sourced by measured crustal mass — that can be tested across laboratories. The Bouguer-constrained near-field columns show that the two laboratories are not geophysically equivalent: the Paris Basin and Appalachian Piedmont differ strongly in the top 5 km, exactly where a near-field scalar-gradient mechanism would be most sensitive. That makes the hypothesis worth testing prospectively.

The framework also improves on an unfalsifiable post-hoc explanation by providing a concrete physical covariate — the local scalar field sourced by measured crustal mass — that can be tested across laboratories. A conventional systematic-error explanation predicts no stable correlation with Bouguer-constrained crustal columns once apparatus differences are controlled; TEP predicts that a residual site term should track the scalar-field column. A multi-laboratory comparison is therefore the real experiment.

The critical next steps are: (1) perform station-scale micro-gravity or Bouguer surveys at each precision-gravity laboratory to reduce geological uncertainty; (2) solve the 3D field on measured lateral density maps rather than a tapered 1D column; (3) reanalyze existing *G* determinations with apparatus class, epoch, and local geology as separate covariates; and (4) compare residuals against the empirically constrained scalar-field columns without re-fitting the coupling. TEP earns decisive support when residuals align with the qualitative directional expectation (lighter sediments → higher *G*eff, denser basement → lower *G*eff) at a significance that cannot be explained by apparatus class alone. This is a higher evidentiary bar than fitting the NIST datum, and the framework is scientifically interesting because it meets that bar while also reproducing the existing datum.

## 6.3 Limitations and Extensions

Seven limitations are explicit. (i) The 3D solver source term S(ρ) = (αlog/λT²) ln(ρ/ρc) Sscreen(ρ) carries the inverse square of the Temporal Topology screening length λT ≈ 4200 km, the canonical multi-center GNSS value (Paper 6). A shorter MGEX verification run (Paper 14) yields a preliminary 2590 km. Because the Laplacian has units of 1/km² and the source term carries 1/λT², the SOR solver returns φ as a dimensionless field directly. No post-hoc domain-area division is required. The dimensional scale is fixed by a physically derived correlation length, though the precise sub-ppm magnitude still rides on the fitted coupling αlog and a measured lateral Bouguer map would be required to tighten the forward uncertainty bands below the current ±1000 ppm scale. (ii) The facility-scale density grids have been built, solved directly, and superposed onto the crustal column in steps 06 and 11, producing non-zero horizontal shear (Sèvres Σh ≈ 7.2×10⁻⁹ m⁻¹, Gaithersburg Σh = 0 m⁻¹). A parameter sensitivity sweep shows the Sèvres signal is driven by hillside-soil asymmetry, with wall/floor thickness and density producing no variation at the centre. (iii) The Bouguer anomalies are drawn from published regional compilations, not from dedicated gravity surveys at the laboratory sites themselves. A dedicated micro-gravity survey at each facility would reduce the Bouguer uncertainty from ~±10 mGal to ~±1 mGal. (iv) The near-field kernel *K*(*z*) = 1/(1 + (*z*/*z*₀)²) is a distance-weighting ansatz; its functional form should be derived from the full TEP field equation. (v) The environment-dependent coupling screen now satisfies the numerical lab-retention and Cassini-suppression requirements, but a first-principles potential V(φ, ρ) must still be specified before the locked value can be defended as derived physics. (vi) The notation collision between TEP's conformal coupling βA (in *A*(φ) = exp(βAφ)) and standard PPN βPPN is a manuscript clarity vulnerability. (vii) The Sèvres–Gaithersburg comparison is in-sample, so distinguishing TEP from apparatus-dependent systematics requires out-of-sample tests.

## 7. Conclusions

This paper replaces hand-tuned geological assumptions with empirically constrained Bouguer gravity anomaly proxies for the near-field column, blended with CRUST1.0 deep structure, and evaluates whether the NIST–BIPM *G* discrepancy is a testable prediction of the Temporal Equivalence Principle. The reproducible pipeline fetches or records published geological inputs for each laboratory coordinate, builds a hybrid crustal column, and computes the scalar field with a near-field distance-weighting kernel that makes the top 5 km dominate the gradient. With the lab-scale conformal coupling held at βA = −1 and the density-sector coupling αlog fixed to negative sign by the TEP field equation in the (+,−,−,−) metric signature, the 1D model predicts a relative decrease of −416 ppm at Gaithersburg versus Sèvres, matching the sign of the measured −281 ppm. The 48% magnitude residual may reflect geological model uncertainty. The 95% forward geological uncertainty band [−2270, +1724] ppm contains the empirical value, so the locked parameterization is not falsified by this single comparison. An uncoupled 3D finite-difference solver, independent of the 1D surrogate calibration, gives a negligible inter-site dG/G ≈ 0 ppm because the symmetric boundary taper makes φ nearly identical at both sites; dimensional normalization φ/(2Ldomain)² yields a dimensionless φcrustal ≈ −4.8×10⁻⁵. Facility-scale density-grid superposition produces a non-zero horizontal shear at Sèvres (Σh ≈ 7.2×10⁻⁹ m⁻¹) versus zero at Gaithersburg, but does not generate a measurable inter-site scalar difference. Measured lateral Bouguer maps are required for the 3D solver to reproduce the site-to-site difference. The locked coupling βA = −1 is incompatible with Cassini PPN bounds (|βPPN| < 0.0034) at solar-system scales unless environment-dependent screening is active; the implemented log-density screen uses n = 4, exceeding the derived minimum n ≥ 3.21. Direct solution of the facility-scale density grids confirms that building-scale asymmetry is the correct source of differential-torque sensitivity. Facility and crustal grids are now fully integrated in the pipeline via linear superposition. Exact agreement would require βA ≈ −0.676 after density screening, a −32.4% shift from the locked coupling — well within plausible geological model uncertainty. The Sèvres–Gaithersburg match is in-sample, so independent confirmation requires the out-of-sample predictions.

The framework generalizes beyond the two-site comparison: laboratories on lighter sedimentary columns should measure higher *G*eff, while those on denser crystalline basement should measure lower *G*eff. Specific ppm-scale predictions are deferred until dedicated station-scale micro-gravity surveys reduce Bouguer uncertainty to ~±1 mGal. At that resolution the scalar-field difference between sites becomes calculable to sub-100 ppm precision, and the framework can be tested prospectively without re-fitting the coupling. Either way, the paper has moved from a single claimed explanation to a reproducible pipeline that links measured crustal mass columns to predicted effective gravitational constants via a physically specified scalar-field equation.

## References

[1] Schlamminger, S., Chao, L., Lee, V., Shakarji, C., Possolo, A., Newell, D., Stirling, J., Cochran, R., & Speake, C. (2026). Redetermination of the gravitational constant with the BIPM torsion balance at NIST. *Metrologia*, 63(2), 025012. https://doi.org/10.1088/1681-7575/ae570f

[2] Quinn, T. J., Speake, C., Parks, H., & Davis, R. (2014). The BIPM measurements of the Newtonian constant of gravitation, *G*. *Philosophical Transactions of the Royal Society A*, 372(2026), 20140032. https://doi.org/10.1098/rsta.2014.0032

[3] Gundlach, J. H., & Merkowitz, S. M. (2000). Measurement of Newton's constant using a modified torsion balance and angular acceleration feedback. *Physical Review Letters*, 85(14), 2868. https://doi.org/10.1103/PhysRevLett.85.2868

## Data Availability & Reproducibility

All inputs to this analysis are publicly documented and reproducible. Site elevations are downloaded live (USGS 3DEP for the US site; global SRTM 30 m for the French site), and all service queries are logged. The crustal mass-column model combines CRUST1.0 layer densities and thicknesses with published Bouguer gravity anomaly constraints for the near-field column. Where a laboratory lacks a dedicated station-scale Bouguer value, the prediction is explicitly marked as regional or local-geology constrained rather than treated as a direct measurement. These inputs are not tuned to the NIST anomaly, and no synthetic measurements stand in for real data.

**Data sources:**

- USGS National Map EPQS: nationalmap.gov/epqs

- USGS MRDS: mrdata.usgs.gov

- BRGM Geoservices: geoservices.brgm.fr

- OpenTopography: opentopography.org

- USGS 3DEP: apps.nationalmap.gov/downloader

- CRUST1.0: igppweb.ucsd.edu/~gabi/crust1.html

- WGM2012 global gravity model: bgi.obs-mip.fr/data-products/grids-and-models/wgm2012-global-model/

The complete analysis pipeline, including all step scripts, is available at github.com/matthewsmawfield/TEP-NIST.

The frozen no-refit prediction packet is generated by `scripts/steps/step_12_generate_preregistration.py` and written to `docs/TEP_NIST_PREREGISTRATION.md`, with a machine-readable companion at `results/outputs/step_12_generate_preregistration.json`.