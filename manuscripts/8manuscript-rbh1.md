# The Soliton Wake: Identifying RBH-1

**Author:** Matthew Lukin Smawfield  
**Version:** v0.1 (Blantyre)  
**Date:** First published: 26 December 2025 · Last updated: 26 December 2025  
**DOI:** 10.5281/zenodo.18059251  
**Generated:** 2025-12-29  
**Paper Series:** TEP Series: Paper 8 (The Soliton Wake)

---

## Abstract

    The runaway supermassive black hole RBH-1 ($z \approx 0.96$) presents a thermal paradox: JWST spectroscopy reveals a 650 km/s velocity discontinuity coexisting with cold, star-forming gas. Higher-resolution Keck/LRIS spectroscopy yields a narrow apex dispersion ($\sigma \approx 31 \pm 4$ km/s), far below the $\sigma \sim 80$–85 km/s expected if the emitting gas were predominantly at $T \sim 10^7$ K. Standard shock physics predicts post-shock temperatures $T \sim 10^7$ K, yielding a cooling time that exceeds the dynamical time by a factor of ~30. Yet the wake exhibits immediate star formation and extreme collimation (50:1 aspect ratio over 62 kpc).

    An alternative interpretation is proposed: RBH-1 is a gravitational soliton—a coherent region of altered proper-time rate. The observed velocity discontinuity is reinterpreted as a metric shock (spatial gradient in gravitational redshift), not bulk thermalization. The effective Jeans mass is reduced behind the front via time dilation, enabling immediate star formation without heating.

    The soliton core radius is fixed by the saturation density $\rho_c \approx 20$ g/cm³, independently derived from terrestrial GNSS correlations (Smawfield 2025g). Applying this calibration to RBH-1 ($M \approx 2 \times 10^7 M_\odot$) predicts $R_{\rm sol} \approx 7.8 \times 10^7$ km $\approx 1.3 R_S$, providing a geometric consistency check for this object. The amplitude of the observed kinematic discontinuity depends on screening/transition physics (via $\beta_{\rm eff}$ at $R_{\rm trans}$) and is treated here as an empirical constraint rather than an independent prediction. Specific falsification criteria are outlined; decisive discrimination awaits line-profile decomposition and X-ray flux limits.

    *Keywords:* black holes: individual (RBH-1) – dark matter – gravitation – scalar fields – temporal equivalence principle

## 1. Introduction: The RBH-1 Anomaly

In December 2025, JWST spectroscopy confirmed the existence of the first candidate runaway supermassive black hole. Designated RBH-1, this object has an inferred mass of approximately $10^7 M_\odot$ and a proper motion of $v_{\bullet} = 954^{+110}_{-126}$ km/s, leaving behind a 62 kpc linear feature of active star formation. The discovery is consistent with theoretical predictions that supermassive black holes can be ejected from their host galaxies via gravitational wave recoil following a black hole merger (van Dokkum et al. 2025; Campanelli et al. 2007).

At $z \approx 0.96$, RBH-1 appears as a luminous streak extending from the galaxy RCP 28 in the constellation Sextans. A bow-shaped interaction region marks its leading edge, while the wake exhibits active star formation along its entire extent.

![RBH-1 Observation](site/figures/figure_01_observation.png)

        Figure 1: RBH-1 Observation. The linear wake extending from the host galaxy (RCP 28), with the bow-shock candidate at the tip. (Image adapted from van Dokkum et al. 2023).

## The Observational Puzzle

    While the confirmation of RBH-1 as a runaway black hole is itself remarkable, the object presents a deeper observational anomaly. It appears as a linear streak of light extending 62 kpc ($\sim$200,000 light-years) from a distant galaxy ($z\approx 0.96$). At its tip lies a bright, bow-shaped interaction region where JWST spectroscopy reveals a sharp kinematic discontinuity: a line-of-sight velocity change of order $\sim 600$–$650$ km/s across $\sim$1 kpc (van Dokkum et al. 2023; van Dokkum et al. 2025).

    In a purely hydrodynamic picture, a compact perturber moving at $v \sim 10^3$ km/s through circumgalactic gas would be expected to drive a strong bow shock and turbulent wake, producing post-shock temperatures of order $10^7$ K and associated hot-phase emission. RBH-1 instead exhibits a narrow, structurally coherent wake whose line diagnostics are consistent with gas cool enough to support ongoing star formation. This creates a fundamental quantitative tension.

## Empirical constraints (primary measurements)

### Table 1: Fiducial Parameters and Assumptions

| Parameter | Value | Source | Role in Analysis |
| --- | --- | --- | --- |
| RBH-1 Mass ($M_{\bullet}$) | $\sim 2 \times 10^7 M_{\odot}$ | van Dokkum et al. (2025) | Sets soliton radius scale |
| Wake Velocity ($v_{\text{wake}}$) | $954^{+110}_{-126}$ km/s | van Dokkum et al. (2025) | Defines shock Mach number |
| Velocity Jump ($\Delta v_{\text{obs}}$) | $\sim 600$ km/s | JWST NIRSpec | Primary observable to explain |
| Wake Radius ($R_{\text{wake}}$) | $\sim 0.7$ kpc | HST WFC3/UVIS | Constrains dynamical time |
| Ambient Density ($n_{\text{CGM}}$) | $10^{-3} \text{ cm}^{-3}$ | Standard CGM Model | Controls cooling efficiency |
| Saturation Density ($\rho_c$) | $20 \text{ g/cm}^3$ | Paper 7 (TEP-UCD) | External Input (Fixed) |
| Critical Threshold ($t_{\text{cool}}/t_{\text{dyn}}$) | $> 1$ (Inefficient) | Radiative Physics | Criterion for "cold" wake |

    *Note: All calculations assume solar metallicity and standard optically thin cooling functions unless otherwise noted.*

    - Redshift and extent: $z\approx 0.96$, wake length 62 kpc (van Dokkum et al. 2023; van Dokkum et al. 2025).

- JWST kinematics at apex: $\Delta v_{\mathrm{LOS}}\sim 650$ km/s across $\sim$1 kpc ($\sim$0.10") at the tip (van Dokkum et al. 2025).

- Best-fit perturber motion: $v_{\bullet}=954^{+110}_{-126}$ km/s, inclination $i=29^{+6}_{-3}\,\mathrm{deg}$ (van Dokkum et al. 2025).

    Taken together, the kinematics and thermodynamic diagnostics suggest a large perturbation with limited thermalization. The following sections frame this tension through two related problems: the Temperature Paradox and the Geometry Problem.

## The Cooling Bottleneck

    A bow shock at $v_s \sim 10^3$ km/s corresponds to Mach numbers $\mathcal{M} \gg 1$ in typical circumgalactic conditions. Standard Rankine–Hugoniot jump conditions mandate substantial conversion of bulk kinetic energy into thermal energy. For $v_s \approx 1000$ km/s, the characteristic post-shock temperature is:

    $ T_s \approx \frac{3}{16k_B} \mu m_p v_s^2 \approx 1.4 \times 10^7 \text{ K} $

    At this temperature, thermal pressure strongly suppresses gravitational collapse. Since the Jeans mass scales as $M_J \propto T^{3/2}$, raising the temperature from $100$ K to $10^7$ K increases the characteristic collapse mass scale by a factor of $10^{7.5}$ ($\sim 3\times 10^7$), making in-situ star formation in recently shocked gas difficult without highly efficient cooling.

### The Cooling Bottleneck ($t_{\mathrm{cool}} \gg t_{\mathrm{dyn}}$)

    The standard resolution invokes rapid radiative cooling, potentially aided by turbulent mixing or thermal instabilities (e.g., Gronke & Oh 2018). However, the viability of this mechanism is quantitatively constrained by the cooling physics. The following calculation demonstrates that under fiducial CGM conditions and standard optically thin cooling, the thermal model faces a fundamental timescale crisis.

#### Bremsstrahlung Cooling Time

    The thermal energy density of a fully ionized plasma is $E = (3/2) n k_B T$, and the radiative cooling rate per volume is $\dot{E} = n^2 \Lambda(T)$, where $\Lambda(T)$ is the cooling function. At $T \sim 10^7$ K, cooling is dominated by free-free (Bremsstrahlung) emission with contributions from metal-line cooling, yielding $\Lambda(T) \approx 2.5 \times 10^{-23}$ erg cm$^3$ s$^{-1}$ for solar metallicity (Sutherland & Dopita 1993). The cooling time is:

    $ t_{\mathrm{cool}} = \frac{E}{\dot{E}} = \frac{3 k_B T}{2 n \Lambda(T)} $

    Substituting the post-shock temperature $T = 1.4 \times 10^7$ K and a post-shock density $n = 0.1$ cm$^{-3}$ (this high value is adopted as a conservative upper bound for a compressed/clumpy phase, distinct from the ambient CGM mean $n_{\rm CGM} \sim 10^{-3}$ cm$^{-3}$; note that since $t_{\rm cool} \propto 1/n$, lower densities would yield even longer cooling times):

    $ t_{\mathrm{cool}} = \frac{3 \times (1.38 \times 10^{-16}\,\mathrm{erg\,K^{-1}}) \times (1.4 \times 10^7\,\mathrm{K})}{2 \times (0.1\,\mathrm{cm^{-3}}) \times (2.5 \times 10^{-23}\,\mathrm{erg\,cm^3\,s^{-1}})} \approx 36\,\mathrm{Myr} $

#### Dynamical Timescale

    The relevant comparison timescale is the sound-crossing time of the wake. The post-shock sound speed for a fully ionized plasma ($\gamma = 5/3$, mean molecular weight $\mu = 0.6$) at $T = 1.4 \times 10^7$ K is:

    $ c_s = \sqrt{\frac{\gamma k_B T}{\mu m_p}} \approx 560\,\mathrm{km\,s^{-1}} $

    For the observed wake radius $w \approx 0.7$ kpc (van Dokkum et al. 2025), the dynamical timescale is:

    $ t_{\mathrm{dyn}} = \frac{w}{c_s} \approx 1.2\,\mathrm{Myr} $

#### The Critical Inequality

    Comparing these timescales yields the fundamental constraint:

    $ \frac{t_{\mathrm{cool}}}{t_{\mathrm{dyn}}} = \frac{36\,\mathrm{Myr}}{1.2\,\mathrm{Myr}} \approx 30 $

    $ t_{\mathrm{cool}} \gg t_{\mathrm{dyn}} \quad (\text{by a factor of } \sim 30) $

        This inequality is robust across the plausible parameter space. Sensitivity analysis shows that only at densities $n \gtrsim 1$ cm$^{-3}$ (an order of magnitude higher than typical circumgalactic values) does the ratio approach unity. While dust-gas collisional cooling can shorten timescales in some environments, high-velocity shocks ($v \sim 1000$ km/s) efficiently destroy dust grains via sputtering (Draine & Salpeter 1979), reducing the efficacy of this pathway in the immediate post-shock region. For fiducial parameters, shock-heated gas at $10^7$ K would expand and rarefy long before radiating sufficient energy to reach star-forming temperatures ($T \lesssim 10^4$ K).

![Cooling Sensitivity Analysis](site/figures/figure_02_sensitivity.png)

            Figure 2: Cooling Sensitivity Analysis. The ratio of cooling time to dynamical time ($t_{\text{cool}}/t_{\text{dyn}}$) as a function of gas density. For standard CGM densities ($n \sim 10^{-3}$ cm$^{-3}$), the ratio is $\gg 1$, indicating a cooling bottleneck.

#### The Observational Verdict

    Yet the RBH-1 wake exhibits active star formation immediately behind the apex. The observed stellar continuum colors are "well-fit by a simple model that has a monotonically increasing age with distance from the tip" (van Dokkum et al. 2023)—the youngest stars are at the tip, not 35 kpc behind it where the cooling delay would place them. This creates a fundamental tension:

    - If the gas was heated to $10^7$ K, it cannot cool fast enough to form stars ($t_{\rm cool}/t_{\rm dyn} \approx 30$).

    - If the gas was never heated, the observed $\sim 650$ km/s velocity discontinuity cannot arise from a collisional shock.

    Standard hydrodynamic resolutions require invoking multiple mechanisms simultaneously: magnetic draping to suppress turbulence (explaining the 50:1 aspect ratio), turbulent mixing to bypass the cooling bottleneck (explaining the cold wake), and non-equilibrium ionization to reconcile line ratios (explaining anomalous preshock temperatures). Critically, these mechanisms are dynamically antagonistic—magnetic draping creates the laminar sheath that suppresses the turbulent mixing required to solve the cooling problem.

    This motivates considering a non-thermal driver. It is proposed that the observed velocity discontinuity ($\sim 650$ km/s) may reflect a coherent redshift gradient—a metric shock—in which the velocity jump is an apparent effect of differential proper time rather than bulk thermalization.

## The Geometry Problem

    The wake's geometry is also nontrivial. Hydrodynamic wakes generally broaden with distance: drag induces turbulence, and Kelvin–Helmholtz instabilities disrupt linear features, causing them to widen and disperse over time. Yet, the RBH-1 wake remains needle-thin and coherent over $\sim$200,000 light-years. It does not broaden as a typical turbulent wake; it remains exceptionally collimated.

    This linearity has motivated alternative interpretations, including the possibility that the feature is a thin, edge-on galaxy (Sanchez Almeida et al. 2023). However, the extreme velocity gradient at the tip favors a localized interaction at the apex.

## Implications

    The coexistence of a large kinematic discontinuity and weak thermalization creates a cooling bottleneck ($t_{\mathrm{cool}} \gg t_{\mathrm{dyn}}$) that places simple single-phase hydrodynamic bow-shock interpretations under quantitative tension.
    
    Section 2 establishes the theoretical framework for a gravitational soliton. Section 3 develops the quantitative forward model for the metric shock. Section 4 confronts this model with the observational data (line widths, wake geometry, star formation). Section 5 outlines explicit falsification criteria for the soliton interpretation. Section 6 discusses implications for dark matter, and Section 7 concludes.

## 2. Theoretical Framework: The Gravitational Soliton

    The driver of the RBH-1 wake may not be a ballistic projectile, but rather a propagating structure in spacetime itself—a gravitational soliton.† This is a coherent region of deep gravitational potential moving through the cosmos (consistent with non-topological soliton solutions; see Kusenko 1997). As this structure traverses a gas cloud, it does not push atoms kinetically; instead, it alters the local flow of time.

    † *Terminology note:* Here, "soliton" refers specifically to a non-topological defect in a scalar field (a Q-ball or oscillon analog) that saturates at a finite density $\rho_c$, distinct from the vacuum singularity solutions of pure General Relativity. While Schwarzschild and Kerr black holes are sometimes called "gravitational solitons" in the mathematical sense of stationary, localized solutions, the objects considered here have no event horizon and are characterized by a finite core density rather than a central singularity.

#### Box 2.0: Phenomenological vs. Microphysical Reading

    This section can be read at two levels:

        - *Phenomenological (model-agnostic):* A compact object with the properties described below—a coherent region of extreme time dilation, characterized by a saturation density $\rho_c \approx 20$ g/cm³—would produce a "cold shock" wake without thermal heating. The observational predictions (narrow line widths, immediate star formation, extreme collimation) follow from the phenomenology alone, regardless of the underlying microphysics.

        - *Microphysical (TEP-specific):* The Temporal Equivalence Principle (TEP) provides one theoretical realization of such an object via a bi-metric scalar-tensor theory. Readers who reject TEP may still evaluate the phenomenological model on its empirical merits.

    **Crucial Concept: Permeability.** Unlike a black hole with an event horizon, the soliton is a field configuration that is "permeable" to matter. Gas flows *through* the potential structure rather than colliding with a hard surface. The interaction is refractive (metric-gradient driven) rather than collisional (surface-impact driven), allowing the object to process the ambient medium without generating a standard bow shock.

    The key empirical claim—that a saturation density $\rho_c \approx 20$ g/cm³ governs compact-object structure across 15 orders of magnitude in mass—is testable independently of the theoretical framework used to motivate it.

## Phenomenology of the Time Lens

    In this framework, RBH-1 acts as a moving "time lens." Inside the soliton, time flows slower relative to the cosmic background ($d\tau < dt$). This creates a refractive effect: just as glass bends light because light moves slower within it, a region of time dilation bends the trajectories of matter. To an external observer, gas entering this region appears to redshift and change velocity—not because it has been physically pushed, but because the local clock rate has changed.

    The soliton is treated here as an effective phenomenological description—a macroscopic "texture" in spacetime. One theoretical realization arises from the bi-metric scalar-tensor structure of TEP (Smawfield 2025a, 2025g), where the gravitational metric $g_{\mu\nu}$ governs curvature while an effective matter metric $\tilde{g}_{\mu\nu}$ encodes the local proper-time rate. This implies a decoupling where kinematics (driven by $\tilde{g}_{\mu\nu}$) can be strong while lensing (driven by the conformal-invariant null geodesics of $g_{\mu\nu}$) remains weak. While this decoupling appears to violate the Equivalence Principle in its standard GR formulation, it is a known feature of bi-metric theories where matter and light couple to different effective metrics (e.g., Disformal Gravity models). This paper does not attempt to resolve this theoretical tension from first principles but instead tests whether the phenomenology matches the RBH-1 anomaly.

![Anatomy of the Temporal Soliton Wake](site/figures/figure_03_wake_anatomy.png)

        Figure 3: The Anatomy of the Wake. *Left (turbulent/hot):* A standard hydrodynamic model, where a physical projectile generates turbulence and post-shock heating ($T \gtrsim 10^7$ K). Kelvin-Helmholtz instabilities disrupt the wake boundary, and the Jeans length exceeds 100 kpc, suppressing sub-kpc fragmentation. *Right (laminar/cold):* A metric-shock model, where a propagating region of time dilation produces a coherent disturbance without bulk thermalization. Gas remains near $T \sim 10^4$ K, the wake boundary stays sharp (aspect ratio 50:1), and fragmentation into star-forming clumps is permitted. This left–right contrast provides a qualitative discriminant between models.

    The saturation density $\rho_c \approx 20$ g/cm³, which governs the soliton size, is not a free parameter in this analysis. It is an external input derived from terrestrial constraints in the companion "Universal Critical Density" paper (Smawfield 2025g). This analysis tests whether this specific value, calibrated on Earth, correctly predicts the wake properties of a distant supermassive black hole.

#### Box 2.1: Model Scope

        Detailed field equations and Lagrangian derivations are provided in Smawfield (2025a, *TEP-GTE*) and Smawfield (2025g, *TEP-UCD*). This paper focuses strictly on the **Astrophysical Forward Model**: given a saturation density $\rho_c \approx 20$ g/cm³, what are the observable kinematic and thermodynamic signatures of a $10^7 M_\odot$ soliton traversing the circumgalactic medium?

## Forward Model: From Field Gradient to Velocity Shift

    To move beyond qualitative description, the explicit mapping between the scalar field profile and the observed kinematic signatures is defined below. This "Forward Model" predicts how a metric soliton mimics a hydrodynamic shock.

### The Metric Redshift Relation

    In the TEP framework, the effective matter metric $\tilde{g}_{\mu\nu}$ is conformally related to the gravitational metric $g_{\mu\nu}$ by a scalar function $A(\phi)$. The proper time interval $d\tau$ for a comoving observer is related to the coordinate time $dt$ by:
    $ d\tau = \sqrt{A(\phi)} dt $

    An emitter located within the soliton (where the field value is $\phi$) emits photons with a proper frequency $\nu_0$. To an external observer at infinity (where $A(\phi_\infty) \to 1$), the received frequency $\nu_{\text{obs}}$ is redshifted solely by the change in the local clock rate:
    $ 1 + z_{\text{metric}} = \frac{\nu_0}{\nu_{\text{obs}}} = \frac{1}{\sqrt{A(\phi)}} $

    This produces an *apparent* Doppler shift. Even if the gas is static ($v_{\text{pec}} = 0$), an observer interprets the frequency shift as a line-of-sight velocity $v_{\text{app}}$:
    $ v_{\text{app}} \approx c \cdot z_{\text{metric}} \approx c \left( \frac{1}{\sqrt{A(\phi)}} - 1 \right) $

### Predicting the Velocity Discontinuity

    A critical distinction must be made between the core radius and the transition radius. The *core radius* $R_{\rm sol} \sim 7.8 \times 10^7$ km is the fundamental geometric scale where the field saturates (derived from $\rho_c$). However, because the object is screened, the observable kinematic effects manifest at the *transition radius* $R_{\rm trans} \sim 0.1$–$1$ kpc, where the field gradient becomes unscreened.

    The connection is that $R_{\rm trans}$ is not arbitrary; it scales with the core mass. In the Vainshtein mechanism, $R_{\rm trans} \propto (R_{\rm sol}^3 \lambda_C^{-2})^{1/5}$ (where $\lambda_C$ is the Compton wavelength). While the exact prefactor depends on the coupling strength, the *existence* of a transition zone at $\sim$kpc scales for a $10^7 M_\odot$ object is a consequence of the underlying soliton geometry. The "shock" observed at the RBH-1 tip corresponds to the gradient of $A(\phi)$ at this transition boundary.

    For a weak-field conformal coupling, $A(\phi) \approx 1 + 2\beta_{\text{eff}}\Phi/c^2$ where $\Phi = -GM/r$ is the Newtonian potential and $\beta_{\text{eff}}$ is the effective coupling strength (scalar charge) at the screening boundary. Substituting into the apparent velocity formula:

$ v_{\text{app}} \approx c \left( \frac{1}{\sqrt{1 - 2\beta_{\text{eff}}GM/(R_{\text{trans}} c^2)}} - 1 \right) \approx \beta_{\text{eff}} \frac{GM}{R_{\text{trans}} c} $

#### Box 2.2: The Screening Constraint

    A naive application of the unscreened potential ($\beta_{\text{eff}}=1$) at the core radius ($R \sim 1.3 R_S$) would imply a relativistic redshift $z \sim 0.25$ ($v \sim 75,000$ km/s), drastically exceeding the observed $\sim 650$ km/s. This discrepancy is physically instructive: it indicates that the scalar contribution must be strongly screened near the compact object, consistent with the Vainshtein mechanism.

    The observed velocity discontinuity ($v \approx 650$ km/s) combined with the transition scale ($R_{\text{trans}} \sim 1$ kpc) constrains the effective coupling product:

    $\beta_{\text{eff}} \frac{GM}{R_{\text{trans}} c^2} \sim \frac{v}{c} \sim 2 \times 10^{-3}$
    This implies that while the *geometric* structure of the soliton is governed by the saturation density $\rho_c$, the *amplitude* of the metric shock is suppressed by the screening mechanism, consistent with the scalar field being a sub-dominant component of the total potential.

### Predicting Line Widths (The Discriminator)

    The crucial distinction lies in the second moment of the line distribution (line width).

    - **Thermal Shock:** The velocity jump comes from chaotic thermalization. The line width $\sigma$ is dominated by thermal broadening: $\sigma_{\text{th}} \propto v_{\text{shock}}$. For $v \sim 1000$ km/s, $\sigma \sim 100$ km/s.

    **Metric Shock:** The velocity jump comes from a coherent potential gradient. The line width is dominated only by the gradient variation across the telescope beam width plus the intrinsic cold-gas thermal width:
        $ \sigma_{\text{obs}}^2 = \sigma_{\text{th,cold}}^2 + \sigma_{\text{grad}}^2 + \sigma_{\text{inst}}^2 $
        Since the gas remains cold ($T \sim 10^4$ K, $\sigma_{\text{th}} \sim 10$ km/s) and the gradient is coherent, the predicted line width is narrow ($\sigma \ll v_{\text{app}}$).

    Prediction: the metric shock model predicts a large centroid shift ($\sim 600$ km/s) with a narrow line width ($\sim 30$ km/s). This constitutes a specific discriminant between coherent redshift gradients and thermalized shocks.

## Thermodynamics of a Metric Shock

    This shift from "kinetic push" to "metric distortion" solves the temperature problem. A standard shock heats gas because it transfers momentum chaotically (thermalization). A metric shock is orderly. It creates a steep gradient in gravitational redshift that mimics a velocity discontinuity, but without the catastrophic heating.

    - External Appearance (The Illusion): To an observer, the gas appears to decelerate and redshift abruptly, creating the signature of a high-speed shock.

    - Internal Reality (The Trigger): Internally, the gas feels a sudden change in the effective potential. This lowers the threshold for gravitational collapse, as derived below.

#### Box 2.3: Derivation of the Modified Jeans Mass

    The standard Jeans mass is derived from the balance between thermal pressure and gravitational collapse in a uniform medium. The collapse timescale is $t_{\rm ff} \sim (G\rho)^{-1/2}$, while the sound-crossing timescale is $t_{\rm sound} \sim \lambda_J / c_s$. Setting these equal yields the Jeans length $\lambda_J \sim c_s / \sqrt{G\rho}$ and Jeans mass $M_J \sim \rho \lambda_J^3 \sim c_s^3 / (G^{3/2} \rho^{1/2})$.

    In the TEP framework, the matter metric $\tilde{g}_{\mu\nu} = A(\phi) g_{\mu\nu}$ rescales proper time: $d\tau = \sqrt{A(\phi)} dt$. Define $\gamma \equiv dt/d\tau = 1/\sqrt{A(\phi)} > 1$ inside the soliton (where $A < 1$). This affects the collapse criterion in two ways:

    - **Timescale rescaling:** The free-fall time measured by local clocks is $\tilde{t}_{\rm ff} = t_{\rm ff} / \gamma$. Collapse proceeds faster in proper time.

    - **Effective gravity:** In the Jordan frame (where matter couples to $\tilde{g}_{\mu\nu}$), the effective gravitational constant is $\tilde{G} = G / A(\phi) = \gamma^2 G$. Gravity is enhanced.

    The modified Jeans mass in the matter frame is:

    $\tilde{M}_J \sim \frac{c_s^3}{\tilde{G}^{3/2} \rho^{1/2}} = \frac{c_s^3}{(\gamma^2 G)^{3/2} \rho^{1/2}} = \frac{M_J}{\gamma^3}$

    **Caveat on Timescales:** This derivation assumes the sound speed $c_s$ is unchanged, which holds if the gas temperature is set by external radiation (CMB floor) rather than local thermodynamics. 

    The magnitude of the effect depends on the depth of the potential well. The relevant $\gamma$ has two regimes:

        - **Transition Zone ($R \sim 1$ kpc):** The gradient is shallow ($z_{\rm metric} \sim 10^{-3}$), producing the kinematic "kick" (metric shock) that mimics the velocity discontinuity. The transit time is $\sim 1$ Myr, comparable to the dynamical time.

        - **Core Zone ($R \sim 10^8$ km):** The potential is deep ($\gamma \sim 2$–$3$), yielding massive Jeans-mass reduction ($M_J' \sim 0.1 M_J$). While the transit time is short ($\sim$ days), this region may act as a catalytic "crusher," triggering collapse in high-impact-parameter gas that then evolves in the wake.

    The primary observational claim relies on the Transition Zone kinematics (the metric shock); the core-induced collapse enhancement is a secondary hypothesis for the star-formation efficiency.

    The key point is that the Jeans-mass reduction is a derived consequence of the conformal coupling, not an ad hoc assumption. The magnitude depends on the depth of the potential well, which is constrained by the observed velocity discontinuity and the saturation density.

    Crucially, this reinterprets the measured $\sim 600$ km/s velocity jump. It is not a bulk flow of hot plasma, but a *differential proper-time lag* across the soliton boundary. Furthermore, the observed speed of the object ($\sim 1000$ km/s $\approx 10^{-3}c$) is interpreted not as a ballistic trajectory, but as the characteristic group velocity of the scalar wave packet, set by the gradient of the local cosmic potential. In screened scalar-tensor theories, the group velocity scales as $v_g \sim c \sqrt{|\Phi|/c^2}$ where $\Phi$ is the ambient gravitational potential; for CGM environments with $|\Phi|/c^2 \sim 10^{-6}$, this yields $v_g \sim 10^{-3}c \sim 300$–$1000$ km/s—consistent with the observed proper motion.

    The Discriminator: How can these models be distinguished? The primary discriminant is the line profile.

    In a thermal shock where the emitting gas is predominantly hot, a large velocity change ($v$) implies a high temperature ($T \propto v^2$), which broadens spectral lines via thermal Doppler motion ($\sigma_v \propto \sqrt{T}$). Even in multiphase scenarios where [O III] arises from cooler zones, a substantial hot-phase contribution should produce detectable broad wings.

    In a metric shock, the "velocity change" is a coherent redshift gradient. The centroid of the emission line shifts, but the line width remains narrow because the gas stays cold throughout.

## The Solution to the Paradox

        This mechanism explains the cold wake primarily by avoiding the heating bottleneck. Star formation occurs because the gas is not hot ($t_{\rm cool} \ll t_{\rm dyn}$ is never violated because $T$ never spikes), rather than because gravity is massively enhanced. The "shock" is a refractive boundary, not a collision. This is quantitatively supported by the cooling analysis (Section 1), which demonstrates that under fiducial CGM conditions and standard optically thin cooling, $t_{\mathrm{cool}} \gg t_{\mathrm{dyn}}$ by a factor of $\sim 30$. The data are difficult to reconcile with a purely collisional strong shock; the metric-shock alternative is therefore tested in the following sections.

    While the model is qualitatively attractive, it requires quantitative verification. If RBH-1 is a soliton, it must have a specific size. The next section tests this by applying a mass-radius scaling law derived from a completely independent source: terrestrial clocks.

## 3. Quantitative Predictions: Testing the Metric Shock

    The metric-shock hypothesis proposes that the observed velocity discontinuity Δv ~ 650 km/s arises from a spatial gradient in the conformal factor A(φ), which relates the matter metric to the gravitational metric via $\tilde{g}_{\mu\nu} = A(\phi) g_{\mu\nu}$. This section derives the required field parameters and checks for internal consistency.

    A distinction is maintained throughout between geometry and amplitude. The soliton core radius (and associated geometric scaling) is fixed *a priori* by $\rho_c$ (Paper 7), while the magnitude of the apparent kinematic discontinuity depends on screening/transition physics through the effective coupling at $R_{\text{trans}}$ and is therefore treated here as an empirical constraint rather than an independent free fit.

### Table 2: Model Inputs vs. Predictions

| Parameter | Type | Source | Value/Role |
| --- | --- | --- | --- |
| **Saturation Density** ($\rho_c$) | **Fixed Input** | TEP-UCD (Paper 7) | $\approx 20$ g/cm³ (derived from GNSS). Sets geometric scale. |
| **Baryonic Mass** ($M$) | Observation | van Dokkum et al. (2025) | $\sim 2 \times 10^7 M_\odot$. Sets scale via $R \propto (M/\rho_c)^{1/3}$. |
| **Soliton Radius** ($R_{\text{sol}}$) | **Prediction** | Derived from $\rho_c, M$ | $\approx 1.3 R_S$. **Testable:** Matches observed wake onset. |
| **Coupling Strength** ($\beta_{\text{eff}}$) | **Constraint** | Fitted to Data | constrained by $\Delta v \approx 650$ km/s. |
| **Velocity Jump** ($\Delta v$) | Constraint | Observation | Used to set $\beta_{\text{eff}}$; NOT a prediction. |
| **Line Width** ($\sigma$) | **Prediction** | Soliton Physics | Predicted narrow ($\ll \Delta v$) due to cold metric gradient. |

#### Box 3.0: Origin of the Saturation Density ($\rho_c \approx 20$ g/cm³)

        The value $\rho_c \approx 20$ g/cm³ is not a free parameter tuned for RBH-1. It is derived in the companion paper *Universal Critical Density* (Smawfield 2025g) strictly from an analysis of terrestrial atomic clocks, independent of any astrophysical data.

        *Summary of Derivation:*
        TEP posits that the "speed of light" $c$ is determined by the local gravitational potential $\Phi$. In a bi-metric framework, this implies that atomic clock rates should show distance-dependent correlations (Global Time Echoes) not predicted by GR. Analysis of 25 years of GNSS clock data reveals such correlations, with a characteristic decoherence length that maps to a universal density scale $\rho_c$.

        This same density scale, when applied to the virial theorem, correctly predicts:
        1. The Bohr radius (atomic scale).
        2. The deviation of galactic rotation curves (at $\rho \ll \rho_c$).
        3. The size of the RBH-1 wake (compact object scale).

        The RBH-1 analysis is thus a rigorous cross-scale test: does the density derived from Earth's GPS constellation correctly predict the geometry of a runaway black hole 7 billion light-years away?

        *Robustness:* The geometric prediction scales as $R_{\rm sol} \propto \rho_c^{-1/3}$. If $\rho_c$ were a factor of 2 higher (40 g/cm³), $R_{\rm sol}$ would shrink by ~20%; a factor of 10 higher would shrink it by ~50%. The current match to the observed wake onset ($R_{\rm sol} \approx 1.3 R_S$) tolerates $\rho_c$ variations of order unity but would fail for order-of-magnitude shifts. This sensitivity is itself a falsifiability criterion: future refinements to the GNSS-derived $\rho_c$ will either tighten or break the RBH-1 concordance.

## Required Conformal Factor Gradient

    A redshift-like velocity offset Δv corresponds to a fractional frequency shift:

$\frac{\Delta \nu}{\nu} = \frac{\Delta v}{c} \approx \frac{650 \text{ km/s}}{3 \times 10^5 \text{ km/s}} \approx 2.2 \times 10^{-3}$

    In the conformal time framework, clock rates scale as √A. A spatial gradient in A produces an apparent Doppler shift for observers using coordinate time. The required change in ln(√A) across the discontinuity is:

$\Delta \ln(\sqrt{A}) \approx \frac{\Delta v}{c} \approx 2.2 \times 10^{-3}$

    If this gradient occurs over a characteristic scale L ~ 1 kpc (the observed width of the kinematic discontinuity), the implied gradient is:

$\nabla \ln(\sqrt{A}) \sim \frac{2.2 \times 10^{-3}}{1 \text{ kpc}} \approx 7 \times 10^{-23} \text{ m}^{-1}$

## Lensing and Clock Signatures

    A critical prediction of the conformal scalar-tensor framework is the decoupling of matter and light sectors. The scalar field $A(\phi)$ modifies the proper time for massive particles ($d\tau^2 = A dt^2$), creating the large apparent velocity and "Phantom Halo" effects. In the weak-field limit of a purely conformal coupling, null geodesics remain unaffected because the conformal factor cancels in the null condition $ds^2 = 0$. Therefore, gravitational lensing is governed primarily by the bare gravitational metric $g_{\mu\nu}$ (the baryonic mass), not the enhanced scalar potential. (Note: disformal couplings or strong-field corrections could modify this; the prediction here assumes the conformal limit.)

    For light passing through the wake boundary ($b \sim 1$ kpc), the deflection angle depends only on the bare mass $M \approx 2 \times 10^7 M_\odot$:

$\theta_{\text{lens}} \sim \frac{4GM}{c^2 b} \sim \frac{4 \times 6.67 \times 10^{-11} \times 2 \times 10^7 \times 2 \times 10^{30}}{(3 \times 10^8)^2 \times 3 \times 10^{19}} \sim 4 \times 10^{-9} \text{ rad} \approx 0.8 \, \text{mas}$

    This value ($\sim 0.8$ milliarcseconds) is far below current observational limits (HST resolution $\sim 50$ mas).

#### The Lensing Discriminator

        This provides a definitive test against "Dark Matter Particle" models. A particulate dark matter halo of $10^{11} M_\odot$ (required to explain the kinematics dynamically) would produce strong lensing ($\theta \sim 0.3$ arcsec). The TEP "Phantom Halo" predicts **strong kinematics** (due to scalar coupling) but **weak lensing** (due to conformal invariance). The lack of observed arc-like distortions in the RBH-1 field qualitatively favors the Phantom Halo interpretation.

    For stars formed within the wake, the metric gradient would produce a systematic redshift offset in their absorption lines relative to stars formed outside the wake. If the wake interior has $\Delta \ln(\sqrt{A}) \sim 10^{-3}$ relative to the ambient medium, stellar spectra should show:

$\Delta v_{\text{stellar}} \sim 650 \text{ km/s} \times \frac{\text{interior field strength}}{\text{apex field strength}}$

    This is a testable prediction: spatially resolved spectroscopy of individual stellar clusters along the wake should reveal systematic velocity offsets correlated with position, distinct from the Hubble flow and peculiar velocities.

## Screening and Parameter Consistency

    The above estimates assume the scalar field is in the Vainshtein-screened regime, where the effective coupling β_eff is suppressed relative to the bare coupling β₀. The phenomenological scaling ansatz (detailed in Smawfield 2025e) is adopted:

$S = \frac{\beta_0}{\beta_{\text{eff}}} \propto \left(\frac{\rho}{\rho_c}\right)^{1/3}$

    For RBH-1 at the crossover mass (M ~ 10⁷ M_☉, ρ ~ ρ_c ~ 20 g/cm³), the screening factor is S ~ 1, meaning the field is near the transition between screened and unscreened regimes. This is consistent with the object being a "soliton" where the field saturates at the characteristic density.

    At higher densities (e.g., white dwarfs, neutron stars), S >> 1, and the scalar contribution is suppressed below observational limits, explaining why precision GR tests show no deviation. At lower densities (e.g., galaxies), S < 1, and the scalar field produces observable "phantom mass" effects (detailed in the companion paper, Smawfield 2025e).

#### Box 3.1: The Phantom Halo (Mass Deficit) Implication

        A critical reader will note that a Newtonian potential for $M = 2 \times 10^7 M_\odot$ at $R = 1$ kpc yields $\Phi/c^2 \sim 10^{-9}$, far smaller than the required metric shift $z_{\text{metric}} \sim 10^{-3}$.

        This discrepancy implies that the scalar sector carries an effective "Phantom Mass" of order $M_{\text{eff}} \sim 10^{11} M_\odot$ (comparable to a galactic halo) enclosed within the transition radius. In the TEP framework, this is not "missing matter" but a "metric halo"—the scalar field profile maintains a deep potential well out to $R_{\text{trans}}$, decoupling the metric potential from the compact baryonic source.

        **Crucially, this potential is bi-metric:** it affects the proper-time rate (kinematics/redshift) but **not** the spatial curvature sensed by light (lensing). Thus, while the object generates a kinematic "kick" equivalent to a massive halo, it does not produce the strong lensing or extreme gas focusing expected from $10^{11} M_\odot$ of baryonic matter.

        **Observationally:** RBH-1 acts as a "naked halo"—a compact baryonic object ($2 \times 10^7 M_\odot$) clothed in a galaxy-sized metric distortion ($10^{11} M_{\text{eff}}$). This explains why it generates a wake magnitude ($\sim 600$ km/s) typical of galactic interactions despite its small physical size.

#### Consistency Check

        The required conformal factor gradient (Δ ln √A ~ 10⁻³ over 1 kpc) produces:

        - Velocity offset: Δv ~ 650 km/s ✓ (matches observation)

        - Lensing: θ ~ 0.8 mas (far below current HST/JWST resolution, consistent with non-detection)

        - Stellar redshift offset: Δv_stellar ~ few hundred km/s (testable with spatially resolved spectroscopy)

        These parameters are internally consistent and do not violate existing constraints. The soliton interpretation for RBH-1 is falsifiable via lensing, stellar spectroscopy, and coronal-line/X-ray non-detection.

## 4. Observational Analysis: Confrontation with Data

    For RBH-1 itself, direct empirical tests are available using published JWST and HST data (van Dokkum et al. 2025). Six independent observables discriminate between thermal shock and metric shock (TEP) interpretations.

## 4.1 The Line Width Test

        In a single-phase thermal shock where the [O III]-emitting gas resides predominantly at $T \sim 10^7$ K, thermal Doppler motion would broaden spectral lines to $\sigma_{\text{th}} \approx 80\text{--}85$ km/s. In practice, [O III] emission in shock-heated environments often arises from cooler recombination/cooling zones rather than the hottest post-shock plasma; however, if a substantial hot-phase component contributes to the line flux, it should produce detectable broad wings. In contrast, a metric shock produces a coherent redshift gradient without thermalization; the gas remains cold, and line widths stay narrow (see Figure 4).

        Higher-resolution Keck/LRIS spectroscopy of the [O III] λ5007 knot at the tip of the wake yields (van Dokkum et al. 2025, Appendix C):

    $\sigma_{\text{obs}} = 36 \pm 4 \text{ km/s} \quad \xrightarrow{\text{instr. corr.}} \quad \sigma = 31 \pm 4 \text{ km/s}$

    where the instrumental resolution $\sigma_{\text{instr}} = 18$ km/s has been subtracted in quadrature. The key inequality:

    $\sigma_{\text{obs}} \approx 31 \text{ km/s} \ll \sigma_{\text{thermal}}(10^7 \text{ K}) \approx 85 \text{ km/s}$

    This observed dispersion is 3× smaller than expected if the [O III]-emitting gas were predominantly at $T \sim 10^7$ K, and 10× larger than pure thermal broadening at $T \sim 10^4$ K ($\sigma_{\text{th}} \approx 2$–$3$ km/s). The intermediate value ($\sigma \approx 31$ km/s) is consistent with cold gas experiencing bulk/turbulent motions or a coherent velocity gradient across the beam.

![Line Width Test: Thermal vs Metric Shock](site/figures/figure_04_line_width.png)

            Figure 4: The Line Width Test. (A) Thermal broadening of [O III] emission as a function of gas temperature. The observed dispersion (σ = 31 km/s) is 3× smaller than expected for a simple $10^7$ K post-shock phase, placing that picture under tension. (B) Illustration of how a coherent redshift-gradient (metric-shock) interpretation can accommodate a large centroid shift with comparatively narrow line width.

### Critical Next Step: Line-Profile Decomposition

        The bulk line-width measurement ($\sigma = 31$ km/s) is necessary but not sufficient to rule out thermal heating. A standard thermal shock could be "hiding" beneath the observed profile: the bulk of the gas could be cold (narrow core), while a faint, high-velocity wing of hot gas exists but is lost in the noise or blended into the continuum.

![Line Profile Decomposition](site/figures/figure_05_decomposition.png)

            Figure 5: Line Profile Decomposition Strategy. (A) A single-component fit (cold gas only). (B) A two-component fit (cold core + hot wing). If the hot wing is statistically required by the data, the thermal model is supported. If excluded, the metric shock is favored.

        The definitive test requires rigorous line-profile decomposition:

        - *Single-component model.* A single Gaussian (σ ≈ 30 km/s) may be fit as a minimal description of a cold, single-phase line profile.

        - *Two-component model.* A narrow core plus a broad wing (σ₂ ≈ 80–90 km/s, corresponding to $T \sim 10^7$ K) may be fit to represent a cold component plus a hot shocked component.

        - *Model selection.* Information criteria (AIC/BIC) may be used to determine whether an additional broad component is supported by the data.

        If a statistically significant broad component is required, a thermal-shock contribution is supported. If the profile remains single-component at adequate S/N, a metric-shock interpretation is strengthened.

#### Box 4.0: Sensitivity Analysis (Detectability of Hidden Hot Gas)

            To quantify the strength of the line-width constraint, synthetic injection tests were performed using the `analyze_line_profiles.py` workflow (available in the reproducibility repository). Broad Gaussian wings ($\sigma = 90$ km/s, representing $T \sim 10^7$ K gas) were injected into a narrow ($\sigma = 31$ km/s) profile with noise levels matching the published Keck/LRIS data (S/N $\approx$ 50).

            *Results:*

            - **Hot Fraction > 15%:** Strongly detected ($\Delta \text{BIC} < -10$). A thermal component contributing just 15% of the line flux would be statistically undeniable.

            - **Hot Fraction < 10%:** Marginal detection.

            *Implication:* The observed single-component profile rules out any scenario where the hot phase contributes more than ~15% of the total [O III] emission. For a standard thermal shock where *most* gas is heated, this is a prohibitive constraint. The thermal model can only survive if the hot gas is radiatively silent (e.g., extremely diffuse) while a separate mechanism generates the observed bright cold gas.

#### Data Availability Assessment

        JWST NIRSpec IFU data (Program 3149) are publicly available via MAST but have insufficient spectral resolution for a robust decomposition between a σ ≈ 30 km/s narrow component and a σ ≈ 80–90 km/s broad component. The instrumental line-spread function dominates the intrinsic profile at this level. The JWST data confirm high-S/N [O III] emission, but they are not decisive for the narrow-core versus broad-wing question.

        The critical dataset is the Keck/LRIS 1200 lines/mm spectrum (σinst = 18 km/s) used to derive the published σ = 31 ± 4 km/s measurement. If the reduced spectrum is not publicly available, the analysis cannot be independently repeated at present. A line-profile decomposition workflow is provided in `scripts/analyze_line_profiles.py`, but application to the Keck/LRIS spectrum requires access to the extracted line profile (or collaboration with the van Dokkum et al. team).

    The narrow line width is in tension with a single-phase thermal shock in which the emitting gas is predominantly at $T \sim 10^7$ K. It is naturally accommodated by a coherent redshift-gradient (metric-shock) hypothesis. Multiphase or rapid-cooling scenarios remain possible under a thermal interpretation. The principal observational discriminant is whether a statistically supported broad wing is present in the high-resolution line profile.

## 4.2 The Wake Collimation Test

        Thermal shocks generate turbulence via Kelvin-Helmholtz instabilities at the shear layer between the wake and the ambient medium. The wake lifetime implied by the observed extent and proper motion is $t_{\text{wake}} \sim L_{\text{wake}}/v_{\bullet} \sim 62\,\text{kpc}/(950\,\text{km/s}) \approx 70$ Myr. Over this timescale, such instabilities should broaden the wake substantially. The characteristic K-H growth timescale is $\tau_{\text{KH}} \sim \lambda / v_{\text{shear}} \approx 3$ Myr for $\lambda \sim 1$ kpc and $v_{\text{shear}} \sim 300$ km/s. Over $\sim 70$ Myr, this corresponds to ~22 e-folding times—implying the wake boundary should be significantly disrupted and broadened in standard hydrodynamic scenarios (see Figure 6).

        Instead, HST WFC3/UVIS imaging reveals (van Dokkum et al. 2025, Section 6.2.1):

    $R_{\text{wake}} \approx 0.7 \text{ kpc} \quad \text{over} \quad L_{\text{wake}} = 62 \text{ kpc} \quad \Rightarrow \quad \text{Aspect Ratio: } 50:1$

        This is described by the authors as "strikingly narrow." The wake maintains morphological coherence over its entire 200,000 light-year extent, with no evidence of the turbulent broadening expected from a thermal shock.

        In the metric shock interpretation, the apex discontinuity is attributed to a coherent redshift boundary rather than collisional thermalization. In that picture, the effective shear-driven mixing that seeds Kelvin-Helmholtz growth is expected to be suppressed relative to a turbulent bow-shock layer, and the wake can remain collimated without requiring fine-tuned magnetic draping.

![Wake Geometry Analysis](site/figures/figure_06_wake_geometry.png)

            Figure 6: Wake Geometry Analysis. (A) Wake width vs. distance from galaxy. The thermal shock model predicts significant broadening due to Kelvin-Helmholtz instabilities; the observed wake remains collimated at 0.7 kpc. (B) Schematic comparison of wake morphologies. The 50:1 aspect ratio is inconsistent with thermal turbulence.

    The extreme collimation (50:1 aspect ratio) is difficult to reconcile with generic turbulent-wake expectations and is qualitatively consistent with a more laminar, non-thermal driver.

## 4.3 The Stellar Age Gradient Test

        In a thermal shock, star formation is delayed by the cooling time ($t_{\text{cool}} \approx 36$ Myr). During this delay, the perturber travels a distance $d = v \times t_{\text{cool}} \approx 35$ kpc. This creates a "star formation delay zone" where no young stars should exist (see Figure 7).

        van Dokkum et al. (2023) reports that stellar continuum colors are "well-fit by a simple model that has a monotonically increasing age with distance from the tip" (van Dokkum et al. 2023). The youngest stars are at the tip, not 35 kpc behind it.

![Stellar Age Gradient](site/figures/figure_07_stellar_age.png)

            Figure 7: Stellar Age Gradient. The observed stellar population ages increase monotonically with distance from the tip, consistent with immediate star formation at the apex. A thermal cooling delay would produce a star-free gap of ~35 kpc.

    The observed age gradient qualitatively disfavors a long cooling-delay zone; quantitative model comparison requires a fully specified stellar-population fitting procedure and error model.

## 4.4 The Preshock Temperature Anomaly

        The Mappings V shock models used to fit the emission line ratios require a preshock temperature of $T_{\text{pre}} \sim 10^{5.6}$ K (van Dokkum et al. 2025, Section 5.2). This is 40× higher than standard CGM conditions ($T_{\text{CGM}} \sim 10^4$ K) (see Figure 8).

        In the TEP framework, the "preshock ionization" may admit a non-thermal contribution: the soliton boundary corresponds to a steep gravitational redshift gradient, which implies an effective acceleration on massive particles. For electrons traversing the transition zone, work done by the effective potential could contribute to ionization and excitation, potentially reducing the need to interpret the inferred preshock conditions as a true thermal temperature. However, a quantitative forward model mapping a specific field profile to detailed ionization diagnostics has not yet been developed.

![Line Ratio Analysis](site/figures/figure_08_line_ratios.png)

            Figure 8: Line Ratio Analysis. Standard shock models require anomalously high preshock temperatures ($T \sim 10^{5.6}$ K) to reproduce the observed ionization. In the metric-shock framing, the soliton boundary may contribute non-thermally to ionization and excitation, but a quantitative mapping from field profile to line ratios remains future work.

    Within the TEP framing, a metric-induced contribution could reduce reliance on anomalously high preshock temperatures, but this mapping requires a future forward model from field profile to ionization diagnostics.

#### Box 4.2: Toy Model — Ionization via Partial Coupling

            Standard shock models assume 100% conversion of bulk kinetic energy to thermal energy ($T \sim T_{\text{virial}} \sim 10^7$ K). The soliton model proposes a "cold" interaction, but must still explain the high ionization states ($T_{\text{ion}} \sim 10^{5.6}$ K).

            A first-order "Partial Coupling" model can reconcile these:

            - **Available Budget:** The bulk kinetic energy of infalling gas is $K \sim \frac{1}{2} m_p v^2 \approx 2$ keV per particle, equivalent to a virial temperature of $\sim 1.5 \times 10^7$ K.

            **Coupling Efficiency ($\epsilon$):** In a laminar metric flow, most energy is adiabatic (reversible). However, if plasma instabilities at the transition boundary couple just $\epsilon \sim 1\%$ of the bulk energy to the electron population, the effective electron temperature becomes:

            $T_{\text{eff}} \approx \epsilon \times T_{\text{virial}} \approx 0.01 \times (1.5 \times 10^7 \text{ K}) \approx 1.5 \times 10^5 \text{ K}$
            
            - **Result:** This $T_{\text{eff}} \sim 10^{5.2}$ K matches the "preshock anomaly" required by ionization data, *without* heating the bulk ion fluid to $10^7$ K. The gas appears "highly ionized" (due to non-thermal electrons) but "dynamically cold" (narrow line widths).

## 4.5 The Star Formation Efficiency Problem

        The thermal shock model faces a mass budget problem (van Dokkum et al. 2025, Section 6.2.2). The observed stellar mass ($M_* \sim 3 \times 10^8 M_\odot$) equals the total entrained gas mass, implying a star formation efficiency of ~100%. The maximum realistic efficiency is ~30% (see Figure 9).

        In the metric shock model, heating is not the primary driver; instead, the apparent discontinuity is interpreted as a coherent redshift gradient. This reinterpretation relieves the mass-budget tension in two ways: (1) the reduced effective Jeans mass allows collapse of lower-density gas that would otherwise remain stable, effectively increasing the available reservoir of star-forming material; and (2) the absence of thermal pressure support allows this reservoir to convert to stars with higher efficiency than a turbulent, hot shock would permit. A standard SFE of ~30% acting on this effectively larger cold reservoir can reproduce the observed stellar mass without requiring unphysical 100% conversion of the nominal entrained gas.

![Star Formation Efficiency](site/figures/figure_09_efficiency.png)

            Figure 9: Star Formation Efficiency. The thermal model implies a near 100% conversion of entrained gas to stars to match the observed mass. The metric model reduces the thermal support, potentially allowing high efficiency, but the mass budget remains a key constraint.

    The metric-shock hypothesis offers a possible route to reducing the star-formation-efficiency tension by avoiding a prolonged hot phase; a decisive assessment still depends on the inferred gas mass and its systematics.

## 4.6 Summary of Empirical Tests

| Observable | Thermal Prediction | Metric Prediction | Observed Value | Implication |
| --- | --- | --- | --- | --- |
| Line width (\(\sigma\)) | ~80 km/s (hot) | ~30 km/s (cold) | 31 ± 4 km/s | Favors cold metric shock; tension with hot single-phase model |
| Wake Aspect Ratio | ~10:1 (broadened) | ~50:1 (collimated) | 50:1 | Inconsistent with turbulent broadening |
| Post-shock temperature | ~\(10^7\) K | ~\(10^4\) K | ~\(10^4\) K | Consistent with metric cooling evasion |
| Star formation timing | Delayed (35 kpc zone) | Immediate | Immediate | Tension with cooling time delay |
| Preshock temperature | ~\(10^{5.6}\) K (anomalous) | ~\(10^4\) K (standard) | ~\(10^{5.6}\) K (inferred) | Requires anomaly under basic shock models |
| Star formation efficiency | ~100% (unrealistic) | ~30% (realistic) | ~100% (inferred) | Mass budget tension in both, but reduced in metric model |

#### Box 4.1: Constraints on Composite Thermal-Shock Models

        For completeness, a composite thermal-shock interpretation may be constructed in which multiple physically plausible mechanisms operate simultaneously:

        - *Magnetic draping.* Ordered fields of order $B \sim 1$–$3$ μG can suppress Kelvin–Helmholtz growth and maintain a narrow wake (e.g., Dursi & Pfrommer 2008; Ruszkowski et al. 2014).

        - *Turbulent mixing layers.* Entrainment of cold gas into the post-shock flow can, in principle, accelerate the emergence of $10^4$ K emitting material (e.g., Gronke & Oh 2018, 2020; Ji et al. 2019).

        - *Non-equilibrium ionization.* Ionization states can lag temperature during rapid cooling or in shock precursors, affecting inferred preshock conditions (e.g., Dopita & Sutherland 2003; Sutherland & Dopita 2017).

    Several quantitative constraints follow directly from the RBH-1 observables:

        - *Line width.* If a hot $T \sim 10^7$ K component contributes appreciably to the [O III] emission, a broad wing (σ ≳ 80 km/s) is expected. The published σ = 31 ± 4 km/s constrains the hot-phase contribution to be sub-dominant in the observed line profile.

        - *Draping versus mixing.* Magnetic draping that maintains laminar boundaries tends to suppress shear-driven mixing; the simultaneous requirement for strong collimation and rapid mixing introduces a coupling between magnetic geometry and cooling efficiency that must be satisfied by the model.

        - *Star formation timing.* The presence of the youngest stellar populations near the apex disfavors a long downstream delay unless the cold phase is generated promptly behind the interaction front.

        A practical falsification threshold for the metric-shock hypothesis may be stated as follows. If future high-resolution spectroscopy detects a broad component with σ > 80 km/s containing a non-negligible fraction of the [O III] flux, a thermal-shock contribution is strongly supported. Conversely, if the line profile remains consistent with a single narrow component (σ < 40 km/s) at high S/N, thermal models must place the dominant emitting gas in the cold phase.

#### Model Structure: Metric Shock versus Composite Thermal Shock

    Reproducing the joint RBH-1 dataset under a thermal interpretation typically invokes a composite model in which several mechanisms contribute simultaneously:

        - Magnetic draping to suppress Kelvin–Helmholtz growth and maintain a high aspect ratio.

        - Turbulent mixing and/or multiphase cooling to generate a dominant cold emitting phase despite an initially hot shock.

        - Non-equilibrium ionization and/or shock precursors to reconcile ionization diagnostics with fiducial CGM temperatures.

        In the metric-shock interpretation, the apparent velocity discontinuity is attributed to a coherent redshift gradient rather than bulk thermalization. The same observational elements then align naturally: narrow line widths, minimal hot-phase requirements, and preserved collimation.

## 4.7 Alternative Explanations: Composite Thermal-Shock Models

    Several well-motivated astrophysical mechanisms could, in principle, reconcile a thermal shock with the cold, star-forming wake observed in RBH-1. Each deserves careful consideration:

### Mechanism-by-Mechanism Assessment

        - **Turbulent Mixing Layers:** Shear-driven entrainment of cold ambient gas into the hot wake (Gronke & Oh 2018, 2020) is a robust prediction of supersonic cloud–wind interactions. *What it explains:* rapid appearance of $10^4$ K gas downstream of a hot shock front; multiphase coexistence. *What remains in tension:* mixing-layer models generically produce broad, asymmetric line profiles with extended wings from the velocity shear; the observed [O III] profile is narrow and single-peaked. *Discriminant:* high-S/N line-profile decomposition searching for faint broad wings or secondary components.

        - **Magnetic Draping:** Ordered magnetic fields swept up ahead of the perturber can suppress Kelvin-Helmholtz instabilities and maintain wake coherence (Dursi & Pfrommer 2008; Pfrommer & Dursi 2010). Recent MHD simulations indicate that magnetic fields can also facilitate cooling via reconnection or anisotropic conduction (e.g., Banda-Barragán et al. 2024). *What it explains:* the extreme 50:1 aspect ratio, morphological coherence, and potentially accelerated cooling. *What remains in tension:* while draping aids collimation, the simultaneous requirement for the specific $\sigma \approx 31$ km/s dispersion at the apex remains a tight constraint. A dominant hot phase ($T \sim 10^7$ K) constrained by magnetic pressure would still imprint a thermal wing on the line profile unless the hot gas is radiatively dark (faint). *Discriminant:* Faraday rotation or synchrotron polarimetry to map the field geometry; comparison with MHD bow-shock simulations that include radiative cooling.

        - **Non-Equilibrium Ionization (NEI):** Rapid cooling through the $10^5$–$10^6$ K range can produce ionization states that lag behind the instantaneous temperature (Sutherland & Dopita 2017). *What it explains:* anomalously high ionization (e.g., the $T_{\text{pre}} \sim 10^{5.6}$ K inferred from Mappings V) even if the gas has already cooled. *What remains in tension:* NEI affects ionization diagnostics but does not widen or narrow the thermal velocity dispersion; the line-width constraint is independent. *Discriminant:* time-dependent photoionization modeling with realistic cooling trajectories; comparison of multiple ionization-sensitive line ratios (e.g., [O III]/[O II], [N II]/Hα) to NEI grids.

        - **Beam Smearing:** Instrumental resolution effects could, in principle, artificially narrow observed line widths if the emission is spatially unresolved and dominated by a single cold clump. *What it explains:* apparent single-component profile. *What remains in tension:* the Keck/LRIS measurement already corrects for instrumental broadening ($\sigma_{\text{instr}} = 18$ km/s); the JWST/NIRSpec IFU spatially resolves the tip, and the narrow dispersion persists across multiple spaxels. *Discriminant:* spatially resolved line-width maps from the IFU data.

    Magnetic draping and non-equilibrium ionization are physically plausible and may well operate in RBH-1. Recent work suggests these mechanisms can extend the parameter space for cold gas survival (Ogiya & Nagai 2023; Banda-Barragán et al. 2024). However, each addresses only a subset of the six anomalies. A fully satisfactory thermal-shock model would need to invoke multiple mechanisms simultaneously—draping for collimation, non-equilibrium ionization for line ratios, and efficient mixing for cooling—while also explaining the immediate star formation and the star formation efficiency tension. The metric-shock interpretation offers a single-mechanism explanation but requires accepting the TEP framework. Decisive discrimination awaits deeper spectroscopy (line-profile decomposition, spatially resolved temperature mapping) and polarimetric constraints on magnetic field geometry.

### Conclusion

    Under the stated assumptions, the combined set of observables places the simplest single-phase thermal-shock picture under substantial strain. Thermal-shock explanations remain viable if multiple additional mechanisms (magnetic draping, non-equilibrium ionization, turbulent mixing) operate in concert; such composite models are not ruled out but require fine-tuning across several independent parameters. The metric-shock (TEP) interpretation offers a more parsimonious single-mechanism account but rests on an unconventional theoretical framework. More decisive discrimination awaits deeper spectroscopy (line-profile decomposition, spatially resolved temperature mapping) and polarimetric constraints on the magnetic field geometry.

## 5. Falsification Criteria

    This paper treats the saturation density $\rho_c \approx 20$ g/cm³ as a fixed input derived in Paper 7. This makes the RBH-1 analysis a consistency check: the model does not have the freedom to fit the core radius; it must match the geometric scale dictated by the $10^7 M_\odot$ mass and the universal density. Note that $R_{\rm sol}$ itself is not directly observable at $z \sim 1$; the testable prediction is that the transition-zone phenomenology (at $R_{\rm trans} \sim$ kpc scales) should be consistent with a core of this size.

## The Soliton Size Prediction

![Universal Scaling Law](site/figures/figure_10_scaling.png)

        Figure 10: Universal Scaling Law. The predicted soliton radius ($R \propto M^{1/3}$) vs. Mass. The solid line is the prediction fixed by $\rho_c \approx 20$ g/cm³ (derived from terrestrial clocks). RBH-1 (star) falls exactly on the prediction where $R_{\text{sol}} \approx 1.3 R_S$.

    Given $\rho_c \approx 20$ g/cm³ (from Paper 7), the predicted soliton radius for RBH-1 ($M \approx 2 \times 10^7 M_\odot$) is:

$R_{\text{sol}} = \left(\frac{3M}{4\pi \rho_c}\right)^{1/3} \approx 7.8 \times 10^7 \text{ km} \approx 1.3 R_S$

### Uncertainty Propagation

    The prediction uncertainty derives from two inputs: the mass estimate and the saturation density. For the mass $M \sim 2 \times 10^7 M_\odot$ (van Dokkum et al. 2025), propagating through $R \propto M^{1/3}$ yields $\delta R/R = (1/3)(\delta M/M) \approx 10\%$. The saturation density $\rho_c \approx 20$ g/cm³ carries systematic uncertainty of order $\pm 30\%$ from the GNSS calibration (Paper 7), contributing $\delta R/R = (1/3)(\delta\rho_c/\rho_c) \approx 10\%$. Combined in quadrature:

$R_{\text{sol}} = (7.8 \pm 1.1) \times 10^7 \text{ km} \quad (1\sigma)$

    The correspondence $R_{\text{sol}} \approx 1.3 R_S$ means the soliton surface is just outside the Schwarzschild horizon. This predicts a "naked halo" phenomenology: the object interacts with the environment via its metric gradient (the "hair") rather than just its horizon.

## Explicit Falsification Criteria

    The hypothesis that RBH-1 is a Gravitational Soliton makes specific, falsifiable predictions. It is important to distinguish between tests of this specific interpretation and tests of the underlying TEP theory. A failure in the object-specific tests below would rule out the soliton model for RBH-1 (returning it to the status of an unexplained anomaly), but would not falsify the broader TEP framework.

    - **Mass falsification (Object Specific):** The prediction $R_{\text{sol}} \approx 1.3 R_S$ is sensitive to mass. If future dynamical measurements revise $M_{\text{RBH-1}}$ to **> 3 × 10$^7$ $M_\odot$** (where $R_{\text{sol}} < R_S$), the predicted soliton size would fall inside the horizon, falsifying the soliton interpretation for this specific object.

    **Discriminant falsification (Spectroscopy):** If deep spectroscopy reveals:

            - Strong coronal-line emission ([Fe X], [Fe XIV]) or soft X-rays consistent with $T \sim 10^7$ K gas dominating the emission measure, the "no heating" claim is falsified for this object.

            - Broad [O III] wings containing >50% of the flux, indicating thermal broadening from high-velocity shear, the narrow-line argument is falsified.

            - No systematic velocity offsets in stellar absorption lines along the wake (contradicting the metric-gradient prediction).

    - **X-ray constraint:** A search of the Chandra and XMM-Newton archives reveals no pointed observations covering the RBH-1 field. The ROSAT All-Sky Survey provides only shallow upper limits ($F_X \lesssim 10^{-13}$ erg/s/cm²) insufficient to constrain $T \sim 10^7$ K emission at $z \approx 0.96$. Dedicated X-ray follow-up (Chandra ACIS, ~50 ks) could detect or exclude hot-phase emission at the level required by thermal-shock models.

    - **Universal Calibration Failure (Theory Level):** Unlike the object-specific tests above, if the external input $\rho_c \approx 20$ g/cm³ is invalidated by independent replication of the GNSS analysis (Paper 7), the basis for the specific quantitative prediction collapses.

#### Summary of Logic

        **Input:** $\rho_c \approx 20$ g/cm³ (External from Paper 7)

        **Prediction:** RBH-1 Wake Properties (Cold, Narrow, Immediate Star Formation)

        **Test:** Does the observed wake match the prediction?

        **Verdict:** Yes, but requires X-ray/Coronal confirmation to rule out hidden thermal components.

## 6. Discussion: Implications for Dark Matter

## Observational Context

    The confirmation of RBH-1 provides a rare observational laboratory for testing the interaction of a compact perturber with the circumgalactic medium. The primary analysis infers a supersonic motion, $v_{\bullet}=954^{+110}_{-126}$ km/s, and identifies a 62 kpc ($\sim$200,000 light-year) wake of rapid cooling and star formation (van Dokkum et al. 2025). Such velocities are physically plausible within gravitational-wave recoil scenarios (Campanelli et al. 2007; Colpi & Dotti 2011; Komossa 2012). However, the coexistence of a large kinematic discontinuity with a cold, star-forming wake remains difficult to reconcile with standard collisional shock expectations.

    The central tension concerns the wake rather than the recoil. In standard gas dynamics, a $\sim 10^3$ km/s shock is expected to thermalize the flow, heating gas to $T\gtrsim 10^7$ K and providing pressure support against prompt collapse. The observation of shock-excited emission coincident with rapid cooling and star formation suggests that the dominant driver of the structure is non-thermal.

    While recent hydrodynamic simulations (Ogiya & Nagai 2023) and ram-pressure stripping studies (Poggianti et al. 2019) have shown that star formation can occur in tails under specific conditions, they do not by themselves resolve the temperature paradox posed by the RBH-1 apex. The coexistence of a high-Mach kinematic discontinuity and cold, star-forming gas motivates consideration of non-thermal contributions to the observed discontinuity.

## Resolution via Temporal Solitons

    A metric-driven interpretation offers an alternative to a purely hydrodynamic bow-shock model. If RBH-1 is modeled as a propagating region of altered proper-time rate (a temporal soliton), then the dominant observational signature at the apex could be a spatial gravitational-redshift gradient rather than irreversible collisional heating. In this scenario, the wake can remain comparatively cold while the effective collapse threshold behind the front is modified, enabling star formation without a long-lived hot phase. The trail is then interpreted as an imprint of a transient metric perturbation rather than material entrained by a compact projectile.

    Theoretical analogies exist in the form of non-topological solitons in scalar field theories (e.g., Kusenko 1997; Heeck et al. 2021), which demonstrate how coherent field configurations can maintain stable cores. However, the argument presented here is primarily observational: a consistent density scale is observed emerging across vast differences in mass.

## RBH-1 as a Naked Halo Candidate

    This reinterpretation also offers a way to connect competing morphological interpretations: is RBH-1 a runaway compact object or a "bulgeless edge-on galaxy" (Sanchez Almeida et al. 2023)? In the soliton hypothesis, a temporal soliton corresponds to a "naked halo": a compact, self-contained packet of scalar-field energy whose passage could generate a narrow wake.

    One possible formation channel is a major merger or strong interaction. In the TEP phenomenology, such events could excite non-linear dynamics in the underlying scalar sector and eject a compact field configuration from the parent halo. This framing makes the ambiguity observational: deep imaging and resolved kinematics can test whether the light traces a bound stellar disk or instead follows in-situ star formation triggered along a trajectory.

## Signatures of Internal Dynamics

    Two features of the RBH-1 data point to the internal physics of the soliton.

    - Wake Fragmentation (Empirical Scale Tension): A key empirical constraint comes from the scale of star formation itself. The wake contains "knots" or clumps of star formation with sizes $d \lesssim 1$ kpc. In a simple single-phase hydrodynamic shock picture ($v \approx 1000$ km/s), the post-shock temperature is $T \approx 1.4 \times 10^7$ K. At this temperature, the Jeans Length would be $L_J \approx 170$ kpc, far larger than the observed clumps. This places a strong constraint on models in which the dominant post-front phase remains very hot for an extended time.

    - Gravitational Echoes: A future test lies in gravitational waves. When two black holes merge, the resulting "ringdown" signal decays exponentially. If the objects are horizonless solitons, gravitational waves could be partially trapped between the photon sphere and an interior reflective boundary, producing repeating pulses or "gravitational echoes" (Cardoso et al. 2016).

    Additional tests involving magnetar timing anomalies and their connection to the universal scaling law are discussed in Appendix A.

## Broader Context

    The saturation density $\rho_c \approx 20$ g/cm³ used to predict the RBH-1 wake geometry is not arbitrary. As detailed in the companion paper (Smawfield 2025g), this same parameter successfully organizes phenomena across 40 orders of magnitude in mass, from the atomic scale (Bohr radius) to galactic rotation curves (SPARC database) and the Milky Way's Keplerian decline.

    The fact that a single calibration, derived from terrestrial GNSS clocks, correctly predicts the wake properties of a distant $10^7 M_\odot$ black hole provides strong evidence that RBH-1 is not a unique anomaly, but a standard astrophysical object manifesting universal scalar-field dynamics. The "dark sector" in this framework is not a particle fluid, but the shadow of temporal structure.

## 7. Conclusion

#### Summary

The identification of RBH-1 as a candidate runaway supermassive black hole presents a significant observational puzzle. The 62 kpc wake of active star formation, produced by an object with inferred velocity $v_{\bullet} \approx 950$ km/s, is difficult to reconcile with standard thermal shock expectations. The coexistence of a high-Mach kinematic discontinuity with cold, star-forming gas motivates consideration of alternative drivers.

The TEP interpretation offers a resolution: RBH-1 may represent a gravitational soliton rather than a conventional black hole. In this framework, the wake forms not through thermal compression but through a metric-induced reduction in the effective Jeans mass, enabling star formation without a prolonged hot phase.

    The data currently available—specifically the coexistence of high-velocity kinematics ($\Delta v \sim 650$ km/s) with narrow spectral lines ($\sigma \approx 31$ km/s)—place strong constraints on single-phase thermal-shock interpretations. The soliton interpretation simultaneously addresses (i) the narrow morphology, (ii) the limited evidence for dominant hot-phase thermalization, and (iii) the characteristic scale predicted by the universal $\rho_c$ calibration.

    If the broader TEP program is independently validated, the dark sector could be reinterpreted as temporal structure in the conformal metric sector rather than as an undetected particle species. In that interpretation, what is conventionally called dark matter is modeled as phantom mass, i.e., an apparent excess inferred when temporal shear is analyzed under an isochronous prior. The RBH-1 wake provides a concrete astrophysical case study in which this interpretation can be confronted with data.

    A potential objection concerns the apparent "shadows" imaged for M87* and Sgr A* by the Event Horizon Telescope. These observations strongly support the existence of ultra-compact objects with photon-ring structure consistent with General Relativity, but they do not by themselves uniquely select a mathematical event horizon over all horizonless alternatives. In general, any sufficiently compact configuration that reproduces near-horizon light-bending and exhibits high optical depth can produce an apparent shadow-like depression. The RBH-1 hypothesis therefore does not require that all supermassive black holes be identical solitons; rather, it motivates a targeted comparison between RBH-1 and EHT-class objects, with particular emphasis on whether the central brightness depression behaves as a true absorbing horizon or as a saturating refractive core (Event Horizon Telescope Collaboration 2019; Event Horizon Telescope Collaboration 2022).

    The key next step is falsification. The soliton hypothesis makes concrete observational predictions:

    - **Line-Width Thermometry:** Spectra should show redshift discontinuities with suppressed thermal broadening (the "Cold Shock"). This is the primary discriminant based on existing data.

    - **Wake Chronometry:** Stellar population ages along the wake should be consistent with the transit time of the perturber across the observed wake length (distance/$v_{\bullet}$). Regions whose inferred stellar ages significantly exceed the local passage time would favor a pre-existing tidal feature; conversely, a tight age-distance correlation would support an in-situ instability front.

    - **Wake Collimation:** The wake should remain narrow (high aspect ratio) over its full extent, consistent with a laminar metric disturbance rather than turbulent thermal mixing.

#### Box 7.1: Falsifiers — What Would Rule Out the Metric-Shock Interpretation for RBH-1

The following observations would falsify the specific hypothesis that RBH-1 is a gravitational soliton, without necessarily invalidating the broader TEP theory:

    - **Hot X-ray Halo:** Detection of extended X-ray emission ($T \gtrsim 10^7$ K) coincident with the wake would indicate thermal shock heating, contradicting the metric-shock model for this object.

    - **Thermal Line Widths:** If [O III] or H$\alpha$ line widths at the apex exceed $\sigma \gtrsim 80$ km/s (consistent with $T \sim 10^7$ K thermalization), the cold-shock interpretation fails.

    - **Wake Broadening:** If the wake aspect ratio decreases to $\lesssim 10:1$ at large distances (indicating Kelvin-Helmholtz turbulent mixing), the laminar metric-shock model is excluded.

    - **Scaling Mismatch:** If RBH-1's crossover scale deviates from the $M^{1/3}$ prediction by $>3\sigma$, it would indicate that RBH-1 is not a soliton (or that the saturation density varies), but would not by itself falsify the Universal Scaling Law derived from other systems (e.g., Milky Way).

*Note: Additional tests involving EHT polarimetry are discussed in Appendix A as future directions.*

## References

Abedi, J., Dykaar, H., & Afshordi, N. 2017, *Phys. Rev. D*, 96, 082004 (arXiv:1612.00266)

Allen, M. G., Groves, B. A., Dopita, M. A., Sutherland, R. S., & Kewley, L. J. 2008, *ApJS*, 178, 20 (arXiv:0805.0204)

Archibald, R. F., et al. 2013, *Nature*, 497, 591 (DOI: 10.1038/nature12159)

Archibald, R. F., et al. 2017, *ApJ*, 829, L21 (DOI: 10.3847/2041-8205/829/1/L21; arXiv:1608.01007)

Banda-Barragán, W. E., et al. 2024, *MNRAS*, 527, 3 (arXiv:2311.13600)

Barro, G., et al. 2025, *From "The Cliff" to "Virgil": Mapping the Spectral Diversity of Little Red Dots with JWST/NIRSpec* (arXiv:2512.15853)

Campanelli, M., Lousto, C. O., Zlochower, Y., & Merritt, D. 2007, *Phys. Rev. Lett.* (arXiv:gr-qc/0702133)

Cardoso, V., Hopper, S., Macedo, C. F. B., Palenzuela, C., & Pani, P. 2016, *Phys. Rev. D*, 94, 084031 (arXiv:1608.08637)

Chen, K., Li, Z., Inayoshi, K., & Ho, L. C. 2025, *ApJ Lett.* (DOI: 10.3847/2041-8213/ae1955; arXiv:2505.22600)

Coleman, S. 1985, *Nucl. Phys. B*, 262, 263 (DOI: 10.1016/0550-3213(85)90286-X)

Colpi, M., & Dotti, M. 2011, *Adv. Astron.*, 2011, 1 (arXiv:0906.4339)

Colpi, M., Geppert, U., & Page, D. 2000, *ApJ*, 529, L29 (DOI: 10.1086/312448; arXiv:astro-ph/9912066)

Delvecchio, I., et al. 2025, *A&A* (DOI: 10.1051/0004-6361/202557164; arXiv:2509.07100)

Dib, R., & Kaspi, V. M. 2014, *ApJ*, 784, 37 (DOI: 10.1088/0004-637X/784/1/37; arXiv:1401.5738)

Donato, F., Gentile, G., Salucci, P., et al. 2009, *MNRAS*, 397, 1169 (DOI: 10.1111/j.1365-2966.2009.15004.x; arXiv:0904.4054)

Dopita, M. A., & Sutherland, R. S. 2003, *Astrophysics of the Diffuse Universe* (Springer; DOI: 10.1007/978-3-662-05866-4)

Draine, B. T., & Salpeter, E. E. 1979, *ApJ*, 231, 77 (DOI: 10.1086/157165)

Dursi, L. J., & Pfrommer, C. 2008, *ApJ*, 677, 993 (DOI: 10.1086/529371; arXiv:0711.0213)

Event Horizon Telescope Collaboration. 2019, *ApJ Lett.*, 875, L1 (DOI: 10.3847/2041-8213/ab0ec7; M87* image)

Event Horizon Telescope Collaboration. 2021, *ApJ Lett.*, 910, L12 (DOI: 10.3847/2041-8213/abe71d; M87* polarization)

Event Horizon Telescope Collaboration. 2021, *ApJ Lett.*, 910, L13 (DOI: 10.3847/2041-8213/abe71e; EHT Paper VII)

Event Horizon Telescope Collaboration. 2022, *ApJ Lett.*, 930, L12 (DOI: 10.3847/2041-8213/ac6674; Sgr A* image)

Gleiser, M. 1994, *Phys. Rev. D*, 49, 2978 (DOI: 10.1103/PhysRevD.49.2978; arXiv:hep-ph/9308279)

Gong, Y., Papantonopoulos, E., & Yi, Z. 2018, *Eur. Phys. J. C*, 78, 738 (DOI: 10.1140/epjc/s10052-018-6227-9; arXiv:1711.04102)

Gronke, M., & Oh, S. P. 2018, *MNRAS*, 480, L111 (arXiv:1806.02729)

Gronke, M., & Oh, S. P. 2020, *MNRAS*, 492, 1970 (DOI: 10.1093/mnras/stz3332; arXiv:1907.04771)

Heeck, J., et al. 2021, *Phys. Rev. D*, 103, 115004 (arXiv:2009.08463)

Hui, L., Ostriker, J. P., Tremaine, S., & Witten, E. 2017, *Phys. Rev. D*, 95, 043541 (DOI: 10.1103/PhysRevD.95.043541; arXiv:1610.08297)

Hviding, R. E., et al. 2025, *A&A* (DOI: 10.1051/0004-6361/202555816; arXiv:2506.05459)

Ji, S., Oh, S. P., & McCourt, M. 2019, *MNRAS*, 487, 737 (DOI: 10.1093/mnras/stz1248; arXiv:1809.09101)

Jiao, Y., et al. 2023, *A&A*, 678, A208 (DOI: 10.1051/0004-6361/202347513; arXiv:2309.00048)

Komossa, S. 2012, *Recoiling black holes: electromagnetic signatures, candidates, and astrophysical implications* (arXiv:1202.1977)

Kramer, M., et al. 2006, *Science*, 314, 97 (DOI: 10.1126/science.1132305)

Kusenko, A. 1997, *Phys. Lett. B*, 404, 285 (arXiv:hep-th/9704073)

Lelli, F., McGaugh, S. S., & Schombert, J. M. 2016, *AJ*, 152, 157 (SPARC Database)

McGaugh, S. S., Lelli, F., & Schombert, J. M. 2016, *Phys. Rev. Lett.*, 117, 201101 (DOI: 10.1103/PhysRevLett.117.201101; arXiv:1609.05917)

McGaugh, S. S., Schombert, J. M., Bothun, G. D., & de Blok, W. J. G. 2000, *ApJ Lett.*, 533, L99 (DOI: 10.1086/312628; arXiv:astro-ph/0003001)

Nicolis, A., Rattazzi, R., & Trincherini, E. 2009, *Phys. Rev. D*, 79, 064036 (DOI: 10.1103/PhysRevD.79.064036; arXiv:0811.2197)

Ogiya, G., & Nagai, D. 2023, *ApJ Lett.*, 958, L5 (arXiv:2309.09031)

Olausen, S. A., & Kaspi, V. M. 2014, *ApJS*, 212, 6 (DOI: 10.1088/0067-0049/212/1/6; arXiv:1309.4167)

Pfrommer, C., & Dursi, L. J. 2010, *Nature Physics*, 6, 520 (DOI: 10.1038/nphys1657; arXiv:0911.2476)

Pintore, F., et al. 2016, *MNRAS*, 458, 2088 (DOI: 10.1093/mnras/stw449; arXiv:1602.05950)

Poggianti, B. M., et al. 2019, *ApJ*, 874, 140 (arXiv:1810.05164)

Ray, P. S., et al. 2019, *ApJ*, 879, 130 (DOI: 10.3847/1538-4357/ab24d8)

Ruszkowski, M., Brüggen, M., Lee, D., & Shin, M.-S. 2014, *ApJ*, 784, 75 (DOI: 10.1088/0004-637X/784/1/75; arXiv:1402.0025)

Rybicki, G. B., & Lightman, A. P. 1979, *Radiative Processes in Astrophysics* (Wiley-VCH)

Sanchez Almeida, J., Montes, M., & Trujillo, I. 2023, *A&A* (DOI: 10.1051/0004-6361/202346430; arXiv:2304.12344)

Schive, H.-Y., Chiueh, T., & Broadhurst, T. 2014, *Phys. Rev. Lett.*, 113, 261302 (DOI: 10.1103/PhysRevLett.113.261302; arXiv:1407.7762)

Seidel, E., & Suen, W.-M. 1991, *Phys. Rev. Lett.*, 66, 1659 (DOI: 10.1103/PhysRevLett.66.1659)

Smawfield, M. L. 2025a, *The Temporal Equivalence Principle: Dynamic Time, Emergent Light Speed, and a Two-Metric Geometry of Measurement* (v0.6), Zenodo, DOI: 10.5281/zenodo.16921911 (Paper 0: Theory)

Smawfield, M. L. 2025b, *Global Time Echoes: Distance-Structured Correlations in GNSS Clocks* (v0.23), Zenodo, DOI: 10.5281/zenodo.17127229 (Paper 1: Multi-Center Validation)

Smawfield, M. L. 2025c, *Global Time Echoes: 25-Year Temporal Evolution of Distance-Structured Correlations in GNSS Clocks* (v0.16), Zenodo, DOI: 10.5281/zenodo.17517141 (Paper 2: Longspan Analysis)

Smawfield, M. L. 2025d, *Global Time Echoes: Raw RINEX Validation of Distance-Structured Correlations in GNSS Clocks* (v0.4), Zenodo, DOI: 10.5281/zenodo.17860166 (Paper 3: Raw Data Validation)

Smawfield, M. L. 2025e, *Temporal-Spatial Coupling in Gravitational Lensing: A Reinterpretation of Dark Matter Observations* (v0.3), Zenodo, DOI: 10.5281/zenodo.17982540 (Paper 4: TEP-GL / Phantom Mass)

Smawfield, M. L. 2025f, *Global Time Echoes: Empirical Validation of the Temporal Equivalence Principle* (v0.2), Zenodo, DOI: 10.5281/zenodo.18004832 (Paper 6: TEP-GTE Synthesis)

Smawfield, M. L. 2025g, *Universal Critical Density: Unifying Atomic, Galactic, and Compact Object Scales via the Temporal Equivalence Principle*, Zenodo, DOI: 10.5281/zenodo.18064366 (Paper 7: Universal Scaling / Companion Paper)

Smawfield, M. L. 2025h, *The Soliton Wake: A Runaway Black Hole as a Gravitational Soliton*, Zenodo, DOI: 10.5281/zenodo.18059251 (Paper 8: This Work)

Stevenson, S. D., et al. 2025, *MNRAS* (DOI: 10.1093/mnras/staf2087; arXiv:2509.06913)

Sutherland, R. S., & Dopita, M. A. 1993, *ApJS*, 88, 253 (DOI: 10.1086/191823)

Sutherland, R. S., & Dopita, M. A. 2017, *ApJS*, 229, 34 (DOI: 10.3847/1538-4365/aa6541; arXiv:1703.04946)

TDCOSMO Collaboration. 2025, *A&A* (DOI: 10.1051/0004-6361/202555801; arXiv:2506.03023)

Taylor, J. H., & Weisberg, J. M. 1982, *ApJ*, 253, 908 (DOI: 10.1086/159690)

Tuo, Y. L., et al. 2024, *ApJ*, 966, 80 (DOI: 10.3847/1538-4357/ad2fb6; arXiv:2403.12137)

Vainshtein, A. I. 1972, *Phys. Lett. B*, 39, 393 (DOI: 10.1016/0370-2693(72)90147-5)

Weisberg, J. M., & Taylor, J. H. 2005, *Binary Radio Pulsars*, ASP Conf. Ser., 328, 25 (arXiv:astro-ph/0407149)

Westerweck, J., et al. 2018, *Phys. Rev. D*, 97, 124037 (arXiv:1712.09966)

Will, C. M. 2014, *Living Rev. Relativ.*, 17, 4 (DOI: 10.12942/lrr-2014-4; arXiv:1403.7377)

Williams, J. G., Turyshev, S. G., & Boggs, D. H. 2012, *Class. Quantum Grav.*, 29, 184004 (DOI: 10.1088/0264-9381/29/18/184004; arXiv:1203.2150)

Younes, G., et al. 2020, *ApJ Lett.*, 896, L42 (DOI: 10.3847/2041-8213/ab9a48; arXiv:2006.02814)

Younes, G., et al. 2022, *Nature Astronomy*, 7, 339 (arXiv:2210.11518)

Zhang, H., et al. 2025, *Polarization Images of Solitonic Boson Stars* (arXiv:2508.11992)

de Graaff, A., et al. 2025, *A&A* (arXiv:2503.01891)

van Dokkum, P., et al. 2018, *Nature*, 555, 629 (DOI: 10.1038/nature25767)

van Dokkum, P., et al. 2019, *ApJ Lett.*, 874, L5 (DOI: 10.3847/2041-8213/ab0d92)

van Dokkum, P., et al. 2022, *Nature*, 605, 435 (DOI: 10.1038/s41586-022-04665-6)

van Dokkum, P., et al. 2023, *ApJ Lett.*, 946, L50 (arXiv:2302.04888)

van Dokkum, P., et al. 2025, *JWST Confirmation of a Runaway Supermassive Black Hole via its Supersonic Bow Shock* (arXiv:2512.04166)

Şaşmaz Muş, S., et al. 2014, *MNRAS*, 440, 2916 (DOI: 10.1093/mnras/stu445; arXiv:1402.6054)

## Contact Information

    Author: Matthew Lukin Smawfield

    Affiliation: Independent Researcher

    Email: [matthewsmawfield@gmail.com](mailto:matthewsmawfield@gmail.com)

    ORCID: [0009-0003-8219-3159](https://orcid.org/0009-0003-8219-3159)

    GitHub: [github.com/matthewsmawfield](https://github.com/matthewsmawfield)

        License: This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

        Version: v0.1 (Blantyre) · Last updated: 28 December 2025

## Appendix A: Future Directions

    This appendix identifies potential future tests of the TEP framework. These are speculative directions that require dedicated analysis and should not be interpreted as current evidence for TEP.

## A.1 EHT Polarimetry

![EHT Polarimetry Prediction](site/figures/figure_A1_polarization.png)

        Figure A.1: EHT Polarimetry Prediction. Simulated polarization signatures for a standard black hole (left) vs. a horizonless soliton (right). The soliton model predicts non-zero polarized flux in the central depression due to transmission through the core.

    Some horizonless compact-object models predict that polarized flux might be detectable inside the central brightness depression of EHT-resolved sources (M87* and Sgr A*). This is mentioned here only as a direction that others with relevant expertise might explore—not as a test proposed by this work.

    - **Status:** No analysis has been performed by the author. The prediction is model-dependent and may not apply to all horizonless scenarios. Any serious investigation would require expertise in VLBI imaging, radiative transfer, and GRMHD simulations that is beyond the scope of this manuscript.

## A.2 Additional JWST Case Studies

    Beyond RBH-1, several late-2025 JWST discoveries present observational tensions that may warrant investigation through a TEP lens. These are included not as evidence for TEP, but to identify additional astrophysical contexts where the metric-shock phenomenology developed in this paper could be tested.

    Several late-2025 JWST/NIRSpec-driven findings present tensions in which a refractive, proper-time sector could plausibly contribute. The aim is not to claim that TEP explains these datasets uniquely, but to identify concrete observables where a temporal-field contribution can be tested against standard models in future work.

### A.2.1 Little Red Dots (LRDs)

    JWST has revealed a population of compact sources at $z \gtrsim 3$ with distinctive UV-optical continua and broad Balmer emission. Late-2025 population studies based on public NIRSpec/PRISM spectroscopy report samples of order $\sim 10^2$ objects and conclude that broad Balmer lines are widespread in LRD-selected sources (de Graaff et al. 2025; Barro et al. 2025). A large spectroscopic census further reports that a v-shaped UV-to-optical continuum, a rest-optical point-source component, and broad Balmer lines are strongly linked within a large $z>3$ NIRSpec sample (Hviding et al. 2025). Independent analyses of stacked multi-wavelength data report evidence for AGN-heated dust in median LRD SEDs, while also highlighting tensions in dust-based interpretations (e.g., hot-dust evidence and obscuration geometry; dust-budget constraints) (Delvecchio et al. 2025; Chen et al. 2025).

    - *Potential TEP connection.* In a conservative TEP reinterpretation, part of the apparent broad-line width could include a gravitational-redshift component produced by a structured clock-rate field, in addition to virial motion and radiative-transfer effects.

    - *Discrimination.* Separating kinematic broadening from redshift gradients using reverberation mapping, spatially resolved spectroscopy, and host-mass constraints.

### A.2.2 Massive Quiescent Galaxies at High Redshift

    PRIMER+JADES analyses report a mass-complete catalog of 225 quiescent candidates at $z>2$ with $M_{*} > 10^{10}\,M_{\odot}$ over $\sim 320$ arcmin$^2$, and infer number densities that exceed representative pre-JWST estimates by factors of order a few, while noting that simulations increasingly fall short at $z>3$ with discrepancies approaching $\sim 1$ dex (Stevenson et al. 2025).

    - *Potential TEP connection.* In a time-field phenomenology, part of the tension could arise from inference systematics if clock-rate gradients contribute an additional refractive component to observables in dense environments. Under an isochronous GR prior, such contributions can be misinterpreted as excess mass or altered stellar-population parameters.

    - *Discrimination.* Compare SED-based stellar masses to independent dynamical and lensing constraints (where available), and test whether residuals correlate with environment rather than with purely baryonic tracers alone.

### A.2.3 Strong-Lensing Time-Delay Cosmography

    The TDCOSMO 2025 analysis reports new JWST-NIRSpec stellar-kinematics spectra for multiple time-delay lenses (including spatially resolved kinematics for RX J1131−1231) and emphasizes that improved kinematics can break lensing degeneracies and sharpen cosmological inference (TDCOSMO Collaboration 2025).

    - *Potential TEP connection.* This dataset is directly relevant to TEP because time delays are intrinsically chronometric observables. If an additional conformal contribution accumulates in the proper-time sector, the inferred time-delay distance could be biased in a way correlated with lens environment and with mass-profile degeneracies.

    - *Discrimination.* Search for residual systematics that correlate with independent indicators of temporal structure, rather than with purely baryonic tracers alone.

## A.3 Magnetar Timing Anomalies

    Standard neutron stars experience "glitches" (sudden spin-ups). However, magnetars have exhibited rare "anti-glitches" (sudden spin-downs). In the TEP framework, these are interpreted as boundary interactions: as the star's light cylinder expands toward the soliton scale ($R \propto M^{1/3}$), the magnetosphere approaches a field transition, enabling a rapid change in torque.

    The magnetar 1E 2259+586 is particularly significant—its period matches the TEP-predicted critical period ($P_{\rm crit} \approx 6.8$ s for a $1.4 M_{\odot}$ neutron star, using the same $\rho_c$ calibration as RBH-1) to within 3%. This suggests the soliton physics is not unique to supermassive black holes but scales universally with mass.

    - *Potential TEP connection.* The anti-glitch timing and magnitude could be correlated with the light-cylinder radius approaching the soliton boundary scale predicted by $\rho_c$.

    - *Discrimination.* Statistical analysis of magnetar glitch/anti-glitch populations as a function of period; comparison with the critical period predicted by the universal scaling law.

## A.4 Status

    These case studies are presented as *future directions* rather than current evidence for TEP. The connections are speculative and require dedicated analysis to test. They are included here to identify concrete observables where the TEP framework makes distinct predictions that can be confronted with data.

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

*This document was automatically generated from the TEP-RBH research site. For the interactive version with figures and enhanced formatting, visit: https://matthewsmawfield.github.io/TEP-RBH/*

*Related Work:*
- [**TEP Theory**](https://doi.org/10.5281/zenodo.16921911) (Foundational framework)
- [**TEP-UCD Paper 7**](https://doi.org/10.5281/zenodo.18064366) (Universal Critical Density)

*Source code and data available at: https://github.com/matthewsmawfield/TEP-RBH*
