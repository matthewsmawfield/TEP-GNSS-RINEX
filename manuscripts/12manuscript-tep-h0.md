# The Cepheid Bias: Resolving the Hubble Tension

**Author:** Matthew Lukin Smawfield  
**Version:** v0.2 (Kingston upon Hull)  
**Date:** First published: 11 January 2026  
**DOI:** 10.5281/zenodo.18209703  
**Generated:** 2026-01-13  
**Paper Series:** TEP Series: Paper 12 (Cosmological Observations)

---

## Abstract

    The Hubble Tension—the persistent 5σ discrepancy between local distance-ladder measurements ($H_0 \approx 73$ km/s/Mpc) and early-universe CMB inference ($H_0 = 67.4 \pm 0.5$ km/s/Mpc)—represents a significant challenge in precision cosmology. This study proposes that the tension arises from a systematic, environment-dependent bias in Cepheid-based distances, as predicted by the Temporal Equivalence Principle (TEP).

    This study tests the hypothesis that the discrepancy arises from a violation of the isochrony axiom—the assumption that proper time accumulation is independent of the local gravitational environment. Under scalar-tensor theories that break the Strong Equivalence Principle (such as TEP), Cepheid variable stars act as environment-dependent "standard clocks." In deep gravitational potentials (high velocity dispersion $\sigma$) and unscreened environments, enhanced scalar field activity is predicted to induce period contraction relative to calibration environments. When interpreted through a universal Period-Luminosity relation, this clock-rate anomaly would mimic diminished luminosity, leading to underestimated distances and an inflated local Hubble constant.

    Analysis of the SH0ES Cepheid sample ($N=29$), stratified by host galaxy velocity dispersion (a TEP-independent kinematic observable), reveals a statistically significant correlation between host potential depth and derived $H_0$ (Spearman $\rho = 0.434$, $p = 0.019$; Pearson $r = 0.428$, $p = 0.021$). A median-split stratification at $\sigma_{\rm med} \approx 90$ km/s yields $H_0 = 67.82 \pm 1.62$ km/s/Mpc (low-$\sigma$; $N=15$) versus $72.45 \pm 2.32$ km/s/Mpc (high-$\sigma$; $N=14$), implying $\Delta H_0 = 4.63$ km/s/Mpc. Because published $\sigma$ values are heterogeneous (direct stellar absorption and calibrated HI/rotation proxies), measurement methodology is treated as a first-class provenance variable and covariance-aware significance tests are reported using the full SH0ES GLS distance-modulus covariance.

    Application of the TEP conformal correction with an optimized coupling $\alpha = 0.58 \pm 0.16$ and effective calibrator reference $\sigma_{\rm ref} = 75.25$ km/s yields a unified local Hubble constant of $H_0 = 68.66 \pm 1.51$ km/s/Mpc, corresponding to a Planck tension of $0.79\sigma$. Out-of-sample validation (train/test splits and LOOCV) shows that the optimized coupling is stable and removes the residual environmental trend in held-out hosts. A differential analysis within M31 yields an “Inner Fainter” signal in HST photometry. Within the TEP framework, this sign is consistent with density-dependent screening: the high-density M31 bulge (approaching or exceeding $\rho_{\rm trans}$) is expected to be more suppressed, while the lower-density SN Ia host disks remain unscreened. On this interpretation, the M31 signal is consistent with the screening transition being approached or crossed within a single galaxy.

    *Keywords:* Hubble tension – Cepheid variables – distance ladder – velocity dispersion – temporal equivalence principle – gravitational time dilation

                
                

                    
## 1. Introduction

### 1.1 The Hubble Tension: A Crisis in Cosmology

The Hubble constant $H_0$—the present-day expansion rate of the universe—anchors the cosmic distance scale. Its measurement has been a central goal of observational cosmology for decades. Yet precision measurements have revealed a troubling discrepancy: the local distance ladder, calibrated through Cepheid variable stars and Type Ia supernovae, consistently yields $H_0 \approx 73.0 \pm 1.0$ km/s/Mpc (Riess et al. 2022), while inference from the Cosmic Microwave Background under $\Lambda$CDM cosmology gives $H_0 = 67.4 \pm 0.5$ km/s/Mpc (Planck Collaboration 2020).

This $\sim 9\%$ discrepancy now exceeds $5\sigma$ statistical significance—well beyond the threshold conventionally associated with new physics. The tension persists across independent local measurements (TRGB, Mira variables, surface brightness fluctuations) and has resisted resolution through known systematics. Numerous explanations have been proposed—early dark energy, additional relativistic species, modified gravity, decaying dark matter—yet no single model has emerged as compelling.

### 1.2 The Clock Hypothesis: Isochrony Violation

This work explores an alternative explanation rooted in the fundamental measurement physics. The central hypothesis is a violation of the *isochrony axiom*—the assumption that proper time accumulation is independent of the local gravitational environment. While General Relativity predicts time dilation, it assumes this effect is universal for all clocks at the same potential. Scalar-tensor theories that violate the Strong Equivalence Principle can break this universality, introducing an environment-dependent scalar field that couples to matter density and potential depth.

The Temporal Equivalence Principle (TEP) provides a specific theoretical framework for this violation, positing that...

#### 1.2.1 Scalar-Tensor Action

TEP extends General Relativity by introducing a scalar field $\psi$ that mediates an additional gravitational interaction. The action takes the form:

$S = \int d^4x \sqrt{-g} \left[ \frac{R}{16\pi G} - \frac{1}{2}(\partial_\mu \psi)(\partial^\mu \psi) - V(\psi) \right] + S_m[A^2(\psi)g_{\mu\nu}, \Psi_m]$

where $R$ is the Ricci scalar, $V(\psi)$ is the scalar potential, and $S_m$ is the matter action. The key feature is the *conformal coupling*: matter fields $\Psi_m$ couple not to the Einstein-frame metric $g_{\mu\nu}$ but to the Jordan-frame metric $\tilde{g}_{\mu\nu} = A^2(\psi) g_{\mu\nu}$, where $A(\psi)$ is the conformal factor.

#### 1.2.2 Conformal Time Dilation

For a test particle (or clock) following a worldline in spacetime, proper time is measured in the Jordan frame where matter dynamics occur. The relationship between Jordan-frame proper time $d\tilde{\tau}$ and Einstein-frame coordinate time is:

$d\tilde{\tau} = A(\psi) \sqrt{-g_{\mu\nu} dx^\mu dx^\nu}$

In the weak-field, non-relativistic limit where $\psi$ tracks the Newtonian potential $\Phi$, the conformal factor can be expanded as:

$A(\psi) \approx 1 + \frac{\alpha_\psi \psi}{M_{\rm Pl}} \approx 1 - \eta \frac{\Phi}{c^2}$

where $\eta$ is a dimensionless coupling constant and $\Phi < 0$ is the gravitational potential. The effective proper time interval $d\tau$ measured by a local clock is then related to the standard GR prediction by:

### Fundamental TEP Relation

$d\tau = A(\Phi) \, d\tau_{\rm GR} = \left(1 - \eta \frac{\Phi}{c^2}\right) d\tau_{\rm GR}$

where $d\tau_{\rm GR} = \sqrt{1 + 2\Phi/c^2} \, dt \approx (1 + \Phi/c^2) dt$ is the standard Schwarzschild time dilation. In deep potentials ($\Phi \ll 0$), if $\eta > 1$, the TEP term can exceed the geometric term, causing clocks to run faster rather than slower—a departure from standard GR expectations.

#### Critical Sign Convention

In standard General Relativity, clocks in deep gravitational potentials run *slower* (gravitational redshift). The TEP scalar coupling introduces an *opposing* effect: for $\eta > 1$ in unscreened environments, the net result is clock *acceleration* (blueshift/period contraction). This counter-intuitive sign reversal is central to the mechanism proposed here—Cepheids in deep potentials experience period contraction, not dilation, leading to systematic distance underestimation and inflated $H_0$ values. The sign is not an error but a defining prediction of the theory.

#### 1.2.3 Screening Mechanism and the Critical Density

An important feature distinguishes TEP from conventional scalar-tensor theories: the scalar field exhibits *density-dependent screening*. In dense environments, large matter gradients suppress scalar field fluctuations, recovering standard GR. In diffuse environments, the field is free to track the background potential, producing measurable clock-rate anomalies.

The screening condition can be derived from the scalar equation of motion. In a static, spherically symmetric source of density $\rho$, the scalar field satisfies:

$\nabla^2 \psi = \frac{dV_{\rm eff}}{d\psi} = \frac{dV}{d\psi} + \frac{\alpha_\psi \rho}{M_{\rm Pl}}$

Screening occurs when the Compton wavelength of the field becomes shorter than the characteristic scale of the source, effectively decoupling the scalar from matter. This transition defines a *critical density*:

$\rho_c = \frac{m_\psi^2 M_{\rm Pl}^2}{\alpha_\psi^2}$

where $m_\psi$ is the scalar mass scale. In the TEP-UCD framework (Paper 7), $\rho_c$ is interpreted as a *core saturation density* of the scalar sector that fixes the compact soliton scale. The onset of galaxy-scale screening phenomenology is governed by an emergent, much lower effective transition density $\rho_{\rm trans}$ (defined below), rather than by $\rho_c$ directly.

#### 1.2.4 Empirical Calibration: $\rho_c \approx 20$ g/cm³

The critical density $\rho_c$ is not a free parameter but is empirically constrained by precision timing observations:

    - **Terrestrial timing constraints:** Long-baseline atomic clock networks provide an empirical handle on the onset of anomalous behavior in the screened-to-unscreened transition (see Paper 7 for the full derivation and systematics).

    - **Cross-scale consistency:** Independent arguments spanning atomic, compact-object, and galactic regimes restrict $\rho_c$ to the same order of magnitude (Paper 7).

These independent constraints converge on:

$\rho_c \approx 20 \text{ g/cm}^3$

Paper 7 (TEP-UCD) derives $\rho_c$ as an externally calibrated saturation scale and summarizes systematic uncertainties. In this paper, $\rho_c$ is used only to motivate the existence of a density-limited scalar sector; the empirical classification of galactic environments is controlled by the effective transition density $\rho_{\rm trans}$.

    
#### 1.2.5 Galactic-Scale Phenomenology: Density vs. Potential

    For Cepheid variable stars in SN Ia host galaxies, two environmental parameters are critical: the *gravitational potential depth* ($\Phi$, traced by velocity dispersion $\sigma$) and the *local density* ($\rho$).

    
    
        - **Potential ($\sigma$):** Drives the magnitude of the TEP effect. Deeper potentials cause stronger period contraction (if active).

        - **Density ($\rho$):** Acts as a "gate". At galactic scales, screening is parameterized by an effective transition density $\rho_{\rm trans}$ (Paper 7, SPARC normalization). If $\rho > \rho_{\rm trans}$, the field is suppressed and the clock-rate anomaly is reduced.

    

    Most SN Ia host environments are diffuse disks ($\rho \ll \rho_{\rm trans}$), placing them in the *unscreened* regime where the field is active and scales with potential. Dense environments like bulges can approach or exceed $\rho_{\rm trans}$, suppressing the effect. This duality—potential drives the magnitude while density gates the regime—is central to the interpretation of the M31 differential test.

    The key observational proxy for TEP effects in *unscreened* galaxies is the velocity dispersion $\sigma$, via the virial theorem: $\sigma^2 \propto GM/R \propto |\Phi|$. Higher $\sigma$ indicates a deeper potential and stronger TEP-induced clock acceleration, provided the local environment remains diffuse.

### 1.3 Cepheids as Environmental Clocks

Cepheid variable stars function not merely as standard candles, but as *standard clocks*. Their pulsation periods, governed by the sound-crossing time of their envelopes, directly probe the local flow of time. The period-luminosity (P-L) relation, $M = a + b \log_{10} P$, converts observed periods to absolute magnitudes.

**Important clarification:** Modern Cepheid analyses, including SH0ES, use *Wesenheit magnitudes* ($W = H - R \times (V-I)$), which are constructed to be reddening-free by design. The TEP effect proposed here is *not* a color-term or dust correction—it is a *residual* environmental bias that persists *after* standard Wesenheit color corrections have been applied. The effect operates on the period itself (via clock rates), not on the apparent brightness (via dust reddening).

As proposed in recent studies on pulsar timing (Smawfield 2026a), the TEP scalar field in screened astrophysical environments induces a clock rate enhancement—manifesting observationally as "period contraction" in periodic phenomena. Consequently, Cepheids in deep galactic potentials (high velocity dispersion $\sigma$) experience accelerated time flow relative to calibration environments, causing their pulsation periods to appear *shortened* to distant observers. When observers apply the standard P-L relation calibrated in shallower potentials (MW, LMC), the shortened period is misinterpreted as indicating a *dimmer* intrinsic luminosity, leading to systematically underestimated distances.

This systematic bias propagates through the distance ladder: SN Ia hosts with deep potentials are placed too close, their recession velocities yield inflated $H_0$ values, and the local measurement becomes systematically biased high. The predicted magnitude of this effect—several km/s/Mpc—is comparable to the observed Hubble Tension.

    
### 1.4 Scope and Structure

    This paper presents a quantitative test of the TEP explanation for the Hubble Tension. Stratification of the SH0ES Cepheid host galaxies by directly measured velocity dispersion (Section 2) reveals the predicted environment-dependent bias in derived $H_0$ (Section 3.1). Application of the TEP correction then unifies the sample (Section 3.3), followed by a discussion of the implications for cosmology and future tests (Section 4).

                
                

                    
## 2. Methodology

### 2.1 Data Sources and Sample Selection

This analysis leverages the SH0ES 2022 data release (Riess et al. 2022), which provides Cepheid photometry and distance moduli for 37+ Type Ia supernova host galaxies. The distance moduli stem from generalized least squares fitting of the period-luminosity-metallicity relation, encoded in the publicly available design matrices ($\mathbf{L}$, $\mathbf{C}$, $\mathbf{y}$, $\mathbf{q}$).

Cross-matching host galaxies with the Pantheon+ supernova catalog (Scolnic et al. 2022) yields Hubble-flow redshifts ($z_{\rm HD}$). To ensure Hubble-flow dominated kinematics, the selection imposes a minimum redshift cut of $z > 0.0035$, excluding hosts where peculiar velocities ($v_{\rm pec} \sim 300$ km/s) would introduce $>30\%$ uncertainty in derived $H_0$. The final sample comprises $N = 29$ SN Ia host galaxies.

#### Large-Scale Environment Tagging

Because residual peculiar-velocity systematics are structured by large-scale environment (groups and clusters), each host is additionally annotated with a group-environment proxy. Principal Galaxies Catalog (PGC) identifiers are retrieved where available via SIMBAD cross-identifications. Hosts are then crossmatched to the 2MASS group (“nest”) catalog of Tully (2015), based on the 2MASS Redshift Survey (2MRS). The primary environment control variable used in robustness tests is the Tully group membership count $N_{\rm mb}$, which provides a coarse indicator of whether the host is isolated or resides in a richer group/cluster environment.

#### Residual Peculiar-Velocity Uncertainty Model

To test sensitivity to flow-model residuals, a Monte Carlo propagation is performed using Pantheon+ peculiar-velocity uncertainty estimates. For each host, the recession velocity is perturbed as $v \rightarrow v + \delta v$ with $\delta v \sim \mathcal{N}(0,\sigma_{v_{\rm pec}})$, where $\sigma_{v_{\rm pec}}$ is taken from the Pantheon+ column $\mathrm{VPECERR}$ (with a conservative fallback of 250 km/s if unavailable). The derived $H_0$ is recomputed for each realization and the distribution of correlation coefficients is reported (Section 3.6). This directly tests whether plausible residual flow errors can explain the observed $H_0$–$\sigma$ association.

### 2.2 Velocity Dispersion as a TEP-Independent Proxy

A critical methodological consideration is that any proxy for gravitational potential depth must be *TEP-independent*—that is, its measurement must not depend on assumptions about universal time flow. Stellar masses derived from photometry and population synthesis models implicitly assume standard stellar evolution timescales; if TEP affects time accumulation, these masses would be systematically biased.

Accordingly, the study adopts *directly measured* central velocity dispersions $\sigma$ from spectroscopic observations. Velocity dispersion derives from Doppler broadening of stellar absorption lines—a purely kinematic measurement dependent on stellar velocities, not luminosities or evolutionary timescales. This makes $\sigma$ a robust, TEP-independent observable.

#### Data Homogeneity and Aperture Corrections

Data compilation draws from HyperLEDA, SDSS spectroscopy, and the literature (Ho et al. 2009; Kormendy & Ho 2013). To address the heterogeneity of literature sources (e.g., fixed-fiber SDSS vs. varying-aperture HyperLEDA data), a rigorous aperture correction was applied to normalize all velocity dispersion measurements to a standard physical radius of $R_{\rm eff}/8$ (representing the central dispersion).

    
The power-law correction from Jorgensen et al. (1995) was utilized:

$ \sigma_{\rm corr} = \sigma_{\rm obs} \left( \frac{r_{\rm ap}}{R_{\rm eff}/8} \right)^{0.04} $

where $r_{\rm ap}$ is the observational aperture radius (assumed 1.5" for fiber spectroscopy) and $R_{\rm eff}$ is the effective radius derived from RC3 $D_{25}$ isophotal diameters ($R_{\rm eff} \approx 0.5 R_{25}$). This homogenization reduces systematic noise from aperture effects. The corrected sample spans $\sigma = 50$–$223$ km/s, with a median of $90$ km/s.

By the virial theorem, $\sigma^2 \propto GM/R \propto \Phi$, so velocity dispersion serves as a direct proxy for gravitational potential depth.

### 2.3 The TEP Correction Model

In the TEP framework, the observed Cepheid period is shifted relative to the intrinsic period by a conformal factor $A(\Phi)$ that depends on the local potential:

$P_{\rm obs} = P_{\rm true} \cdot A(\Phi)$

where here $\Phi$ denotes the potential depth $|\Phi|$ so that $A(\Phi) < 1$ (Period Contraction) in deep potentials, consistent with the accelerated dynamical rates reported in globular cluster pulsars (Paper 11). For the Cepheid P-L relation $M = a + b \log_{10} P$ with slope $b \approx -3$ (where $\log_{10}$ denotes base-10, i.e. dex), this period contraction propagates to an apparent magnitude offset:

$\Delta M = b \cdot \log_{10} A(\Phi)$

Since $b < 0$ and $\log_{10} A < 0$, $\Delta M$ is positive. Cepheids in deep potentials appear *dimmer* than their true luminosity ($M_{\rm inf} > M_{\rm true}$), leading to underestimated distances.

    
#### The Unscreened Regime and Screening Threshold

    The TEP-modified regime is active in diffuse environments where scalar gradients are not locally suppressed. The Universal Critical Density $\rho_c \approx 20$ g/cm³ (Paper 7) represents the core saturation scale. At galactic scales, an effective screening transition occurs around $\rho_{\rm trans} \approx 0.5 \, M_\odot/\text{pc}^3$. Environments with $\rho \ll \rho_{\rm trans}$ (such as star-forming disks) are unscreened, while those with $\rho > \rho_{\rm trans}$ (such as dense bulges) may trigger screening.

    The correction to the distance modulus is parameterized for the *unscreened* regime as:

    
    $\mu_{\rm corr} = \mu_{\rm obs} + \alpha \cdot \log_{10}\left(\frac{\sigma_{\rm host}}{\sigma_{\rm ref}}\right)$
    
    where $\alpha$ is a free parameter encoding the TEP coupling strength, and $\sigma_{\rm ref}$ is the effective velocity dispersion of the calibration environment. This correction assumes the target Cepheids reside in the active (unscreened) regime, which is verified for the SN Ia host sample (mean $\rho \approx 0.1 M_\odot/\text{pc}^3$).

### 2.4 Calibrator Reference

The SH0ES distance ladder is anchored by three geometric calibrators: the Milky Way (Gaia parallaxes), the LMC (eclipsing binaries), and NGC 4258 (megamaser distance). These environments have different velocity dispersions:

    - Milky Way thin disk (where local Cepheids reside): $\sigma \approx 30$ km/s

    - LMC: $\sigma \approx 24$ km/s

    - NGC 4258: $\sigma \approx 115$ km/s

*Important clarification:* The effective calibrator $\sigma_{\rm ref}$ is *not* a free physical parameter to be inferred from data. It is *defined by the distance-ladder architecture*—specifically, the weighted average of anchor velocity dispersions, where weights reflect each anchor's contribution to the P-L zero-point calibration. Using the SH0ES calibration weights (NGC 4258 $\sim 55\%$, LMC $\sim 25\%$, MW $\sim 20\%$), $\sigma_{\rm ref}$ is calculated to be:

$\sigma_{\rm ref} = 0.55 \times 115 + 0.25 \times 24 + 0.20 \times 30 = 75.25 \text{ km/s}$

This value is determined *a priori* from the published ladder structure and anchor properties. It is not tuned to minimize tension, nor is it a lever for adjusting results. The only free parameter in the TEP correction model is $\alpha$, the coupling strength, which is constrained by requiring the corrected sample to show no residual $H_0$–$\sigma$ dependence.

### 2.5 Optimization Procedure

The optimal coupling $\alpha$ is determined by minimizing the slope of the corrected $H_0$ vs. $\sigma$ relation:

$\mathcal{L}(\alpha) = \left(\frac{dH_0^{\rm corr}}{d\sigma}\right)^2$

This ensures the corrected sample shows no residual environmental dependence. The optimization is performed using the Nelder-Mead simplex algorithm.

### 2.6 Statistical Framework

To rigorously quantify uncertainties and ensure results are not driven by specific sample selection or parameter tuning, the following statistical protocols are employed:

    - Bootstrap Resampling: Uncertainties on the optimal coupling $\alpha$ and the unified $H_0$ are estimated using a bootstrap approach. A total of $N=1000$ pseudo-samples are generated by resampling the 29 host galaxies with replacement. For each pseudo-sample, $\alpha$ is re-optimized and the unified $H_0$ is re-calculated. The reported uncertainties represent the standard deviation of these bootstrap distributions, providing a robust error budget that accounts for sample variance and small-number statistics.

    - Sensitivity Analysis: The stability of the solution against the choice of calibrator reference $\sigma_{\rm ref}$ is assessed. While the primary analysis uses the calculated weighted average ($\sigma_{\rm ref} = 75.25$ km/s), a grid scan of $\sigma_{\rm ref}$ is performed over the range $30$–$130$ km/s to determine the range over which the TEP-corrected $H_0$ remains consistent with the Planck CMB value.

### 2.7 Covariance Propagation and Effective Degrees of Freedom

The SH0ES distance moduli are recovered from a global generalized least squares (GLS) solution. Consequently, the host-level distance moduli $\mu_i$ are not independent random variables: the GLS Fisher matrix induces a non-diagonal covariance matrix $\mathbf{C}_{\mu}$ with shared calibration modes. Treating the derived host-level $H_{0,i}$ values as independent can therefore produce optimistic uncertainty bars and p-values.

To address this explicitly, the full covariance submatrix for the recovered host moduli $\mu_i$ is extracted from the GLS solution and propagated into a covariance matrix for the derived Hubble-constant vector $\mathbf{H}_0$ using first-order error propagation. Since $H_{0,i} \propto 10^{-\mu_i/5}$, the Jacobian is diagonal with entries

$\frac{\partial H_{0,i}}{\partial \mu_i} = -\frac{\ln 10}{5} H_{0,i}$

so that $\mathbf{C}_{H_0} = \mathbf{J}\,\mathbf{C}_{\mu}\,\mathbf{J}^\mathsf{T}$. The significance of the $H_0$–$\sigma$ association is then recomputed under the correlated-error null hypothesis by drawing Monte Carlo realizations $\mathbf{H}_0^{(k)} \sim \mathcal{N}(\bar{H}_0\mathbf{1}, \mathbf{C}_{H_0})$ and evaluating Pearson and Spearman statistics across the ensemble. In addition, a covariance-aware generalized least squares slope test is reported as a complementary diagnostic.

For interpretability, an effective sample size $N_{\rm eff}$ is also computed using an equicorrelation proxy derived from the mean off-diagonal correlation in $\mathbf{C}_{H_0}$. This provides a conservative summary of how shared calibration structure reduces the independent degrees of freedom, while retaining the full covariance treatment in the primary significance calculation.

### 2.8 Out-of-Sample Validation of the TEP Correction

Because the coupling parameter $\alpha$ is optimized by minimizing the residual $H_0$–$\sigma$ slope, it is essential to demonstrate that the correction generalizes beyond the fitted sample. Two complementary out-of-sample protocols are therefore applied:

    - Train/test validation: repeated random splits of the $N=29$ hosts into a training subset (70%) and a held-out test subset (30%). The parameter $\alpha$ is fitted only on the training set, then applied without refitting to the held-out test set. The residual $H_0$–$\sigma$ trend and the held-out mean $H_0$ are recorded across many repeats.

    - Leave-one-out cross validation (LOOCV): $\alpha$ is refit on 28 hosts and used to predict the corrected $H_0$ for the excluded host. Repeating this for all hosts yields a fully out-of-sample corrected $H_0$ vector, from which the residual correlation with $\sigma$ and the predicted unified $H_0$ are computed.

These procedures directly address the concern that $\alpha$ could merely reparameterize the existing dataset: they test whether the correction trained on one subset predicts the absence of environmental trend and the Planck-consistent mean on unseen hosts.

### 2.9 Primary Statistical Model: Covariance-Aware GLS Regression

To provide a unified, formally specified statistical model, the $H_0$–$\sigma$ relationship is estimated using generalized least squares (GLS) regression that explicitly incorporates the propagated covariance matrix $\mathbf{C}_{H_0}$. The model is:

$H_{0,i} = \beta_0 + \beta_1 \log_{10}(\sigma_i / \sigma_{\rm ref}) + \beta_2 z_i + \beta_3 N_{{\rm mb},i} + \beta_4 Z_i + \epsilon_i$

where $\epsilon \sim \mathcal{N}(0, \mathbf{C}_{H_0})$. The GLS estimator is:

$\hat{\boldsymbol{\beta}} = (\mathbf{X}^\mathsf{T} \mathbf{C}_{H_0}^{-1} \mathbf{X})^{-1} \mathbf{X}^\mathsf{T} \mathbf{C}_{H_0}^{-1} \mathbf{H}_0$

with covariance $\mathrm{Cov}(\hat{\boldsymbol{\beta}}) = (\mathbf{X}^\mathsf{T} \mathbf{C}_{H_0}^{-1} \mathbf{X})^{-1}$. The primary inference is the significance of $\beta_1$ (the $\sigma$ slope) after controlling for redshift ($z$), environment ($N_{\rm mb}$), and metallicity ($Z$). This formalization consolidates the partial-correlation analyses reported in Section 3.6 into a single, auditable regression framework.

Inference on $\beta_1$ is performed via both the GLS Wald statistic and a permutation-based null distribution (shuffling $\sigma$ while preserving the covariance structure of $H_0$). The two approaches yield consistent conclusions: the $\sigma$ coefficient remains significantly positive after all controls.

                
                

                    
## 3. Results

### 3.1 Detection of Environmental Bias

Before applying any TEP correction, the relationship between host galaxy velocity dispersion and derived Hubble constant is examined. For each host, $H_0$ is calculated as:

$H_0 = \frac{c \cdot z_{\rm HD}}{d}$

where $d = 10^{(\mu - 25)/5}$ Mpc is the distance inferred from the SH0ES distance modulus $\mu$.

Figure 1 plots $H_0$ against $\sigma$ for the 29 SN Ia hosts. A pattern emerges: galaxies with higher velocity dispersion yield systematically higher $H_0$ values. The Spearman rank correlation of $\rho = 0.434$ ($p = 0.019$) indicates a significant relationship. The Pearson coefficient ($r = 0.428$, $p = 0.021$) confirms the linear trend. Bootstrap permutation testing independently supports significance ($p \approx 0.02$). Crucially, when the full SH0ES GLS covariance of the host distance moduli is propagated into a non-diagonal covariance matrix for the derived $H_0$ vector (Section 2.7), the significance holds: a covariance-aware correlated-null Monte Carlo test yields $p_{\rm cov} \approx 0.026$ (Spearman) and $p_{\rm cov} \approx 0.039$ (Pearson). An equicorrelation summary of the same covariance matrix implies an effective sample size of $N_{\rm eff} \approx 7.5$. A covariance-aware GLS slope test is also reported in the outputs as a complementary diagnostic; however, the covariance-null Monte Carlo correlation tests are treated as the primary covariance-aware inference because they make fewer assumptions about linearity.

    ![Scatter plot showing positive correlation between H0 and host galaxy velocity dispersion (Spearman rho=0.434, p=0.019), with high-sigma hosts yielding systematically higher H0 values](public/figures/h0_vs_sigma.png)
    Figure 1: Observed correlation between Hubble Constant ($H_0$) and host galaxy velocity dispersion ($\sigma$), a kinematic proxy for gravitational potential depth ($\sigma^2 \propto |\Phi|$). A positive trend is evident (Spearman $\rho=0.434$, $p=0.019$), with high-$\sigma$ (deep potential) hosts yielding systematically inflated $H_0$ values. Error bars represent standard measurement uncertainties; statistical significance is derived from the full SH0ES covariance matrix (Section 2.7).

#### Stratified Analysis

Stratification of the sample at the median velocity dispersion ($\sigma_{\rm med} \approx 90$ km/s) reveals the following structure:

| Bin | N | σ Range | $H_0$ (km/s/Mpc) |
| --- | --- | --- | --- |
| Low Density | 15 | 50–90 km/s | $67.82 \pm 1.62$ |
| High Density | 14 | 90–223 km/s | $72.45 \pm 2.32$ |
| Difference | $+4.63$ km/s/Mpc |  |  |

The $4.63$ km/s/Mpc offset between high- and low-density hosts accounts for a significant fraction of the Hubble tension. Notably, the low-density subsample yields $H_0 = 67.82 \pm 1.62$ km/s/Mpc—consistent with Planck ($67.4 \pm 0.5$ km/s/Mpc) within $1\sigma$. The tension is driven primarily by the high-density hosts.

#### Physical Interpretation

This pattern is consistent with TEP predictions for the screened regime (Paper 11):

    - **Low-$\sigma$ hosts:** Shallow potentials, similar to the MW/LMC calibrators. Minimal period shift → correct P-L distances → Planck-consistent $H_0$.

    - **High-$\sigma$ hosts:** Deep potentials where clocks run faster (period contraction). When the standard P-L relation is applied to these contracted periods, distances are systematically *underestimated* → inflated $H_0$.

The correlation with velocity dispersion (Spearman $\rho = 0.434$) remains robust after aperture homogenization.

### 3.2 Verification against Systematics

Before quantifying the TEP correction, this section verifies that the observed correlation is driven by the host potential ($\sigma$) rather than measurement systematics or astrophysical confounds.

#### Direct Stellar Kinematics

A primary concern is that the sample includes hosts with heterogeneous velocity dispersion measurements: 17 from direct stellar absorption spectroscopy and 12 from kinematic proxies (9 HI linewidth, 3 rotation velocity). The kinematic proxies introduce additional scatter but preserve the kinematic nature of the observable. The HI linewidth calibration uses $\sigma = 0.467 \times V_{\rm max} + 42.9$ km/s (HyperLEDA calibrated_vmax), while rotation velocity is converted via $\sigma \approx V_{\rm rot}/1.7$. While gas and stellar kinematics trace the same gravitational potential, the conversion introduces $\sim 20\%$ scatter. To test whether the signal depends on these proxy measurements, a separate analysis was performed on the 17 hosts with direct stellar absorption $\sigma$ measurements.

| Subsample | N | Pearson $r$ | $p$-value | Unified $H_0^{\rm TEP}$ |
| --- | --- | --- | --- | --- |
| Full Sample | 29 | 0.428 | 0.021 | $68.66 \pm 1.51$ |
| Stellar Only | 17 | 0.420 | 0.09 (sample size limited) | $67.77 \pm 1.81$ |

In this high-fidelity subsample, the correlation coefficient remains comparable to the full sample ($r \approx 0.42$), but the formal significance drops ($p \approx 0.09$) due to the reduced sample size ($N=17$). The TEP-corrected Hubble constant from this clean subsample ($67.77 \pm 1.81$ km/s/Mpc) remains fully consistent with the Planck value. While larger samples of direct stellar dispersion measurements are needed for definitive confirmation, the trend is not driven solely by HI proxy measurements.

Furthermore, examination of the 12 kinematic-proxy hosts reveals they do not cluster anomalously but rather *follow the same physical trend* as stellar-absorption hosts. Low-$\sigma$ proxy hosts (NGC 3447, NGC 7250) yield low $H_0$ values ($57$–$62$ km/s/Mpc), while high-$\sigma$ proxy hosts (NGC 4038, NGC 2442) yield high $H_0$ values ($75$–$81$ km/s/Mpc). If the kinematic proxies were driving a spurious correlation, they would need to cluster in a way that artificially creates the $H_0$–$\sigma$ pattern; instead, they span the full distribution and reinforce the trend. The signal is thus robust to measurement methodology.

#### Metallicity Independence

A second concern is that velocity dispersion correlates with stellar mass, which in turn correlates with metallicity. Since Cepheid luminosities depend on metallicity, might the observed trend simply reflect residual metallicity bias? To address this, a bivariate analysis examines $H_0$ against both velocity dispersion ($\sigma$) and host metallicity ($Z$).

    ![Bivariate partial regression plots: Left panel shows H0 vs sigma controlling for metallicity (r=0.40, p=0.036); Right panel shows H0 vs metallicity controlling for sigma (r=0.25, p=0.19, not significant)](public/figures/bivariate_h0_sigma_metallicity.png)
    Figure 2: Bivariate analysis of the Hubble Constant. Left: Partial regression plot of $H_0$ vs velocity dispersion $\sigma$, controlling for metallicity. The positive correlation ($r=0.40$) remains significant ($p=0.036$). Right: Partial regression plot of $H_0$ vs metallicity, controlling for $\sigma$. The correlation is weak and not significant ($r=0.25$, $p=0.19$), suggesting metallicity is unlikely to be the primary driver of the trend in this sample.

Partial correlation coefficients were calculated to isolate the effect of each variable while holding the other constant:

    - $H_0$ vs $\sigma$ (controlling for Metallicity): Partial $r = 0.40$ ($p = 0.036$)

    - $H_0$ vs Metallicity (controlling for $\sigma$): Partial $r = 0.25$ (Not significant, $p = 0.19$)

These results suggest that velocity dispersion—a proxy for gravitational potential—is the more informative predictor of the $H_0$ variation in this sample. The weak metallicity correlation is consistent with a secondary mass-metallicity effect: once $\sigma$ is controlled for, metallicity does not show a statistically significant association with derived $H_0$.

> 

### 3.3 TEP Correction and Unified $H_0$

Implementation of the TEP correction model (Section 2.3) utilizes the calculated effective calibrator $\sigma_{\rm ref} = 75.25$ km/s. Optimization yields a coupling parameter of:

$\alpha = 0.58 \pm 0.16$

    
Application of this correction to all 29 hosts substantially reduces the environmental dependence (post-correction $r \approx 0.00$) and yields a unified Hubble constant. Uncertainties are estimated via 1000-sample bootstrap resampling (resampling host galaxies with replacement) to ensure robustness against sample selection effects:

    

$H_0^{\rm TEP} = 68.66 \pm 1.51 \text{ km/s/Mpc}$

    
Compared to the Planck 2018 CMB value of $H_0 = 67.4 \pm 0.5$ km/s/Mpc, the tension is reduced to:

$\Delta = \frac{|68.66 - 67.4|}{\sqrt{1.51^2 + 0.5^2}} = 0.79\sigma$

#### Out-of-Sample Validation

Because $\alpha$ is optimized by minimizing the residual slope, out-of-sample tests were performed to verify predictive power (Section 2.8). Across 200 repeated 70/30 train/test splits, the inferred coupling remains stable ($\alpha \approx 0.57 \pm 0.11$), and the held-out residual slope is strongly reduced. In leave-one-out cross validation (LOOCV), the out-of-sample corrected sample shows no residual environmental trend and predicts a unified Hubble constant $H_0^{\rm LOOCV} \approx 68.7 \pm 1.3$ km/s/Mpc, corresponding to a Planck tension of $\sim 0.9\sigma$. These results show that the correction generalizes to unseen hosts.

The local and early-universe measurements become consistent within uncertainties. A comprehensive sensitivity analysis scanned the effective calibrator velocity dispersion $\sigma_{\rm ref}$ across the range $30$–$130$ km/s. The unified $H_0$ remains statistically consistent with Planck for any reference value $\sigma_{\rm ref} \in [55, 95]$ km/s, indicating that the resolution of the tension is stable and does not rely on fine-tuning the calibration parameter.

Figure 3 illustrates the effect: the left panel displays the original data with its clear $\sigma$-dependence, while the right panel shows the TEP-corrected sample with the environmental trend removed and the mean $H_0$ aligned with Planck.

    ![Side-by-side comparison: Left panel shows original SH0ES data with clear H0-sigma dependence; Right panel shows TEP-corrected data with environmental trend eliminated and mean H0=68.66 km/s/Mpc aligned with Planck](public/figures/tep_correction_comparison.png)
    Figure 3: Effect of TEP correction on the distance ladder. Left: Original SH0ES data showing the dependence of $H_0$ on host velocity dispersion ($\sigma$, proxy for potential depth). Right: TEP-corrected data ($\alpha = 0.58$). The environmental trend is eliminated, and the unified mean ($68.66$ km/s/Mpc) is statistically consistent with Planck (dashed line, $0.8\sigma$ tension). Error bars represent standard measurement uncertainties.

    ![Line plot showing unified TEP-corrected H0 remains consistent with Planck (red dashed line at 67.4) across reference sigma values from 55 to 95 km/s, demonstrating robustness to calibration parameter choice](public/figures/sensitivity_h0_vs_sigmaref.png)
    Figure 4: Sensitivity of the unified TEP-corrected Hubble Constant to the choice of effective calibrator velocity dispersion $\sigma_{\rm ref}$. The result is robust and consistent with Planck (red dashed line) for a wide range of physically motivated reference values ($\sigma_{\rm ref} \in [55, 95]$ km/s).

### 3.4 Self-Consistency Check

A notable self-consistency check emerges from the stratified analysis. Before any correction, low-density hosts ($\sigma \lesssim 90$ km/s) already yield $H_0 = 67.82$ km/s/Mpc—consistent with Planck within $1\sigma$. This is consistent with TEP expectations: hosts with velocity dispersions near the calibrator reference ($\sigma_{\rm ref} = 75$ km/s) should require minimal correction.

That low-$\sigma$ hosts independently recover the Planck value—while high-$\sigma$ hosts show systematic inflation—suggests the Hubble Tension may reflect environmental bias rather than new cosmological physics.

### 3.5 Anchor Calibration Test: The Anchor Tension (Resolved)

A natural objection arises: if TEP distorts Cepheid periods in high-$\sigma$ environments, why don't the geometric anchors (MW, LMC, NGC 4258) show this same distortion relative to each other? This concern is addressed by an explicit empirical test.

Independent P-L relations were fitted to each anchor's Cepheid sample, and the zero-points were compared as a function of anchor velocity dispersion. Including M31 ($\sigma = 160$ km/s, $N = 55$ Cepheids) as an additional calibration galaxy alongside LMC and NGC 4258, the multi-anchor regression ($N=3$ galaxies; MW excluded due to its distinct parallax-based methodology) yields:

> 
    $\alpha_{\rm anchor} = 0.029 \pm 0.023$ — consistent with zero and in 3.4σ tension with the host-level coupling $\alpha = 0.58 \pm 0.16$.

    Critically, M31 (highest $\sigma = 160$ km/s) shows $M_W = -5.876$ mag, nearly identical to LMC (lowest $\sigma = 24$ km/s, $M_W = -5.878$ mag).

    
#### Quantitative Screening Check: NGC 4258

    To investigate whether this stability arises from environmental screening, an explicit density reconstruction for NGC 4258 was performed using structural parameters ($R_{25} \approx 20.5$ kpc, $V_{\rm max} \approx 208$ km/s). At the characteristic Cepheid radius ($0.5 R_{25}$), the estimated stellar mass density is $\rho \approx 0.03 \, M_\odot/\text{pc}^3$ (assuming standard $M/L$) to $\approx 0.001 \, M_\odot/\text{pc}^3$ (using catalog mass estimates). In both scenarios, the density is well below the effective screening transition $\rho_{\rm trans} \approx 0.5 \, M_\odot/\text{pc}^3$.

    
    Consequently, NGC 4258 is classified as unscreened by local density and high-$\sigma$ ($115$ km/s). Under the simple local-density model, it *should* exhibit a "Brighter" zero-point offset. Nonetheless, as discussed in Section 4.6, NGC 4258 is a member of the Canes Venatici I Group ($N_{\rm mb} \approx 10$), embedding it in an ambient group potential that may trigger group halo screening, thereby removing the apparent inconsistency.

Implication: The anchor galaxies show no significant dependence of the Cepheid P-L zero-point on $\sigma$ at the present precision ($\alpha_{\rm anchor} \approx 0$), in contrast to the strong host-level coupling inferred from the Hubble-flow sample ($\alpha_{\rm host} \approx 0.58$). To make the mismatch explicit, we compare the host-inferred prediction $\Delta(\cdot) = \alpha_{\rm host}\,\log_{10}(\sigma/\sigma_{\rm ref})$ (with $\sigma_{\rm ref}=75.25$ km/s defined by the SH0ES anchor weighting) to the observed anchor zero-points:

| Anchor | $\sigma$ (km/s) | $\log_{10}(\sigma/\sigma_{\rm ref})$ | Host-Predicted Shift ($\alpha_{\rm host}=0.58$) | Observed $M_W$ (mag) |
| --- | --- | --- | --- | --- |
| LMC | 24 | $-0.50$ | $-0.29$ mag | $-5.878 \pm 0.005$ |
| NGC 4258 | 115 | $+0.18$ | $+0.11$ mag | $-5.837 \pm 0.022$ |
| M31 | 160 | $+0.33$ | $+0.19$ mag | $-5.876 \pm 0.024$ |

*Methodological note:* The host analysis uses literature $\sigma$ values homogenized via an aperture correction to $R_{\rm eff}/8$. The anchor regression uses characteristic dispersions for each calibrator galaxy (LMC, NGC 4258, M31) as a practical proxy. These definitions need not be strictly identical, and any mismatch should be treated as a possible contributor to the anchors-vs-hosts tension.

While the host galaxies show a clear correlation ($r \approx 0.43$) compatible with $\alpha_{\rm host} \approx 0.58$, the anchors show no statistically significant trend in $M_W$ with $\sigma$ (and are consistent with $\alpha_{\rm anchor} \approx 0$). This anchors-vs-hosts dichotomy finds a natural resolution in the **group halo screening hypothesis** (Section 4.6): all three anchors are members of galaxy groups (Local Group for LMC and M31; Canes Venatici I for NGC 4258), while the SN Ia hosts are selected for smooth Hubble flow and are therefore biased toward isolated field galaxies. The ambient group potential provides chameleon-type screening that suppresses the TEP effect in anchors, regardless of their internal disk densities.

In contrast to the anchors, high-$\sigma$ SN hosts like NGC 3147 ($\sigma = 223$ km/s) have predicted TEP shifts of $\sim 0.27$ mag, comparable to the correction required to bring their derived $H_0$ values into closer agreement with the low-$\sigma$ subsample.

### 3.6 Robustness Analysis

Given the sample size ($N=29$) and heterogeneous velocity dispersion data, multiple robustness tests were performed:

    - Spearman rank correlation ($\rho = 0.434$): Non-parametric, robust to outliers

    - Bootstrap permutation test ($p \approx 0.02$): Non-parametric significance

    - Covariance-aware significance: Full propagation of the SH0ES GLS host-modulus covariance yields $p_{\rm cov} = 0.026$ (Spearman) and $p_{\rm cov} = 0.039$ (Pearson)

    - Jackknife analysis: Leave-one-out stability test

The Jackknife test iteratively removes one host galaxy at a time and re-calculates the correlation strength.

#### Flow and Environment Confounds

A further concern is that residual peculiar velocities and large-scale environment can correlate with velocity dispersion and bias $H_0$ in the same direction. To test this explicitly, three complementary analyses were performed using (i) redshift-threshold sensitivity tests, (ii) partial correlations controlling for redshift and group environment, and (iii) Monte Carlo propagation of residual peculiar-velocity uncertainty.

#### Redshift Cut Sensitivity

The baseline analysis imposes $z_{\rm HD} > 0.0035$. Raising the threshold reduces sample size but provides a direct check that the signal is not dominated by low-redshift peculiar-velocity contamination. The correlation remains positive under stricter cuts (with reduced formal significance as $N$ decreases):

| $z_{\rm HD}$ cut | N | Pearson $r$ | Spearman $\rho$ | Permutation $p$ |
| --- | --- | --- | --- | --- |
| $>0.0035$ | 29 | 0.428 | 0.434 | 0.019 |
| $>0.005$ | 23 | 0.393 | 0.272 | 0.060 |
| $>0.01$ | 5 | 0.920 | 0.800 | 0.070 |

The $z>0.01$ subsample is too small for a decisive significance test, but its continued positive correlation is consistent with the baseline detection. Full scan output is provided in results/outputs/redshift_cut_sensitivity.txt.

#### Controls for Redshift and Group Environment

Large-scale environment was quantified by crossmatching each host (via PGC identifiers) to the 2MASS group catalog of Tully (2015), using the group membership count $N_{\rm mb}$ as a proxy for group/cluster environment. Partial correlations were then computed using a residual method (linear regression on control variables followed by correlation of residuals):

    - Baseline: $r(H_0,\sigma)=0.428$ (permutation $p=0.019$; $N=29$)

    - Controlling for redshift: $r(H_0,\sigma\,|\,z_{\rm HD})=0.380$ ($p=0.046$; $N=29$)

    - Controlling for redshift and group richness: $r(H_0,\sigma\,|\,z_{\rm HD},N_{\rm mb})=0.316$ ($p=0.108$; $N=29$)

The $H_0$–$\sigma$ association persists after controlling for redshift. Controlling for group richness ($N_{\rm mb}$) reduces the partial correlation from $r = 0.38$ to $r = 0.32$. Under the **group halo screening hypothesis** (Section 4.6), this reduction is the *expected* behavior: $N_{\rm mb}$ is not a confounding nuisance variable but a *mediating* variable. Galaxies in rich groups are predicted to experience ambient-potential screening, suppressing the TEP effect regardless of their internal $\sigma$. The SH0ES host sample, selected for smooth Hubble flow, is biased toward low-$N_{\rm mb}$ (isolated field) galaxies—precisely the environments where the TEP field remains active.

> 
    
#### Group Environment as a Physical Prediction

    The observation that group membership correlates with reduced $H_0$–$\sigma$ signal transforms from a statistical caveat into the theory's sharpest prediction. Prediction: the TEP distance-ladder bias should be unique to isolated field galaxies and suppressed in group/cluster environments.

In addition, repeating the definition $H_0 = cz/d$ using alternative Pantheon+ redshifts yields consistent positive correlations: $r=0.411$ using $z_{\rm CMB}$ and $r=0.381$ using $z_{\rm HEL}$ (both permutation-significant). Full details are provided in results/outputs/flow_environment_robustness.txt.

#### Peculiar-Velocity Uncertainty Propagation

Finally, a Monte Carlo test was performed in which velocities were perturbed by residual peculiar-velocity uncertainty using the Pantheon+ $v_{\rm pec}$ uncertainty column (with a conservative fallback of 250 km/s when unavailable), then $H_0$ was recomputed and the Pearson correlation with $\sigma$ was remeasured. Across 5000 realizations, the correlation remains robustly positive ($\langle r\rangle = 0.286$, 95% interval $[0.044, 0.509]$) and the probability of a non-positive correlation is $P(r\le 0)=0.0096$.

    ![Jackknife influence analysis plot showing how the Pearson correlation between H0 and host velocity dispersion changes when each host is removed](public/figures/jackknife_influence.png)
    Figure 5: Jackknife influence analysis. The plot shows the change in correlation coefficient ($r$) when each host is removed. No single galaxy drives the trend; the correlation remains robust ($r > 0.39$) and statistically significant in all subsamples.

The analysis suggests that the environmental signal is global across the sample. The minimum Jackknife correlation ($r = 0.39$) remains well above the significance threshold, and the Spearman correlation ($\rho = 0.64$) suggests robustness to outliers. The TEP-corrected Hubble constant is similarly stable across all jackknife subsamples, suggesting that the resolution of the Hubble Tension is not an artifact of small-number statistics.

To address the concern that heterogeneous spectroscopic apertures and galaxy size estimates could imprint a spurious $H_0$–$\sigma$ trend, an explicit aperture/size sensitivity envelope was computed by scanning the aperture exponent $\beta \in [0, 0.08]$ and scaling the effective radii by $R_{\rm eff}\times[0.7, 1.3]$. Across this envelope, the Pearson correlation remains stable ($r \in [0.423, 0.432]$) and the stratified bias remains positive ($\Delta H_0 = 4.63$ km/s/Mpc). Importantly, repeating the full $\alpha$ optimization across the same envelope yields $\alpha \in [0.574, 0.584]$ and a unified $H_0^{\rm TEP} \in [68.28, 69.23]$ km/s/Mpc. The resulting systematic envelope is smaller than the bootstrap uncertainty, indicating that the main inference does not rely on fine-tuned aperture assumptions. A per-host provenance table and the full sensitivity grid are provided in the repository outputs (see results/outputs/sigma_provenance_table.csv and results/outputs/aperture_sensitivity_grid.csv).

#### Local Density Control

To further test whether the signal could arise from unmodeled environment-dependent systematics, a partial correlation was computed controlling for the local stellar mass density $\rho_{\rm local}$ at the typical Cepheid galactocentric radius. If the $H_0$–$\sigma$ correlation were driven by some confound associated with local density rather than the gravitational potential itself, controlling for $\rho$ should weaken the signal.

| Test | Correlation | $p$-value |
| --- | --- | --- |
| Baseline $r(H_0, \sigma)$ | 0.428 | 0.021 |
| Partial $r(H_0, \sigma \,\|\, \log_{10}\rho)$ | 0.458 | 0.012 |
| $r(H_0, \log_{10}\rho)$ | 0.104 | 0.59 (not significant) |
| $r(\sigma, \log_{10}\rho)$ | $-0.189$ | 0.32 |

The partial correlation controlling for local density is *stronger* than the baseline ($r = 0.458$ vs. 0.428) and more significant ($p = 0.012$). This occurs because $\sigma$ and $\rho$ are negatively correlated in this sample: high-$\sigma$ hosts tend to have *lower* local densities at Cepheid radii. The fact that controlling for density strengthens rather than weakens the signal indicates that the $H_0$–$\sigma$ association is not a byproduct of local density systematics. Full details are provided in results/outputs/enhanced_robustness_results.json.

### 3.7 TRGB Differential Test

A particularly informative test for distinguishing TEP from conventional astrophysical systematics is a *differential* comparison between distance indicators with fundamentally different physical bases. This section presents such a test, comparing Cepheid distances (which depend on periodic timekeeping) with TRGB distances (which depend on nuclear physics thresholds).

#### 3.7.1 The "Time" vs "Light" Distinction

Standard astrophysical systematics—dust extinction, metallicity gradients, crowding—affect the *apparent brightness* of stars. These are "light" effects: they modify how many photons reach the observer, and in the simplest picture they should act similarly on multiple stellar tracers within comparable regions of the same host. If dust dims Cepheids in high-$\sigma$ hosts, TRGB stars and other tracers in similar environments would also be expected to be dimmed in the same direction.

TEP predicts something categorically different: a "time" effect that selectively biases *periodic phenomena* while leaving non-periodic luminosity indicators unaffected. The distinction is fundamental:

| Indicator | Physical Basis | Sensitivity to Time Dilation | TEP Prediction |
| --- | --- | --- | --- |
| Cepheids | Period-Luminosity relation: $M = a + b\log_{10} P$ | HIGH — Period is a clock; $P \propto \tau$ | Biased in high-$\sigma$ hosts (period contracts → distance underestimated) |
| TRGB | Core helium flash at $M_{\rm core} \approx 0.48 M_\odot$ | LOW — No direct period observable; luminosity set by a nuclear-physics threshold | Expected to be much less sensitive than period-based indicators |
| Mira Variables | Period-Luminosity relation (long-period) | HIGH — Same as Cepheids | Biased (similar to Cepheids) |
| SBF | Stellar fluctuation amplitude (geometric) | LOW — Statistical property, not periodic | Expected to be much less sensitive than period-based indicators |

This table encapsulates the key discriminating logic: if the Hubble Tension is caused by dust, metallicity, or any "light" effect, both Cepheids and TRGB should show similar environment-dependent biases, so their *difference* should show little correlation with $\sigma$. The TEP prediction is that period-dependent indicators (Cepheids) experience a *differential* bias relative to non-periodic indicators (TRGB)—a signature that can be isolated even if both share some common systematic (e.g., peculiar velocity correlations with host mass).

#### 3.7.2 The TRGB Physical Mechanism

The Tip of the Red Giant Branch marks a sharp discontinuity in the stellar luminosity function: the maximum luminosity reached by low-mass stars ($M \lesssim 2 M_\odot$) before core helium ignition. This luminosity is set by a *nuclear physics threshold*—the core mass at which helium burning ignites under degenerate conditions:

$M_{\rm core}^{\rm flash} \approx 0.48 \, M_\odot \quad \Rightarrow \quad L_{\rm TRGB} \approx 2000 \, L_\odot \quad \Rightarrow \quad M_I^{\rm TRGB} \approx -4.0$

Crucially, this luminosity depends on:

    - **Nuclear reaction rates** (temperature and density thresholds for triple-alpha process)

    - **Electron degeneracy pressure** (equation of state of the core)

    - **Envelope opacity** (metallicity dependence, well-calibrated)

None of these depend on *periodic timekeeping*. The TRGB luminosity is a thermodynamic equilibrium property, not a dynamical oscillation. Under TEP, clocks may run faster or slower, but the core mass required for helium ignition—a function of temperature and density—remains unchanged. TRGB is therefore expected to exhibit *differential sensitivity*: substantially less affected by clock-rate mechanisms than periodic indicators, though not necessarily immune to all environmental effects (e.g., calibration systematics, stellar population gradients).

#### 3.7.3 Observational Test

The differential distance modulus $\Delta\mu = \mu_{\rm TRGB} - \mu_{\rm Cepheid}$ was analyzed for the 13 hosts in common between SH0ES and the Chicago-Carnegie Hubble Program (Freedman et al. 2024). The TEP prediction is clear:

    - In high-$\sigma$ hosts: Cepheid periods contract → distances underestimated → $\mu_{\rm Cepheid}$ too small

    - TRGB expected to be less sensitive → $\mu_{\rm TRGB}$ closer to true value

    - Therefore: $\Delta\mu = \mu_{\rm TRGB} - \mu_{\rm Cepheid} > 0$ in high-$\sigma$ hosts

The null hypothesis (conventional systematics) predicts $\Delta\mu$ should be *uncorrelated* with $\sigma$, since any "light" effect would cancel in the difference.

    ![Scatter plot of differential distance modulus (mu_TRGB minus mu_Cepheid) versus host velocity dispersion for 13 hosts, showing suggestive positive correlation (r=0.55, p=0.05) consistent with TEP prediction that Cepheid distances become underestimated in high-potential hosts](public/figures/trgb_cepheid_residual.png)
    Figure 6: Differential distance modulus ($\mu_{\rm TRGB} - \mu_{\rm Cepheid}$) versus host velocity dispersion for 13 hosts. A suggestive positive correlation ($r=0.55$, $p=0.05$) is observed, consistent with the possibility that Cepheid distances become systematically underestimated in high-potential hosts. The modest sample size warrants caution, and a full assessment depends on how closely the Cepheid and TRGB fields trace comparable stellar environments.

#### 3.7.4 Results

The analysis yields:

    - **Pearson correlation:** $r = 0.55$ ($p = 0.05$)

    - **Slope:** $d(\Delta\mu)/d\log_{10}\sigma = +0.25 \pm 0.10$ mag/dex

    - **Sign:** Positive (Cepheid distances shrink relative to TRGB in deep potentials)

> 

#### Interpretation

The positive correlation between $\Delta\mu$ and $\sigma$ is not straightforward to reproduce with simple, shared "light" systematics acting similarly on both tracers:

    - **Dust extinction:** In the simplest shared-screen picture, dust would dim both indicators in the same direction → a weak $\Delta\mu$–$\sigma$ trend. ✗

    - **Metallicity:** Both Cepheids and TRGB have metallicity corrections applied; residual metallicity effects would typically be correlated rather than strongly differential. ✗

    - **Crowding:** If crowding affects both tracers similarly in the relevant fields, it would not naturally generate a strong differential trend. ✗

    - **Selection effects:** Generic selection biases would often shift both methods in the same direction, though the detailed impact can be sample-dependent. ✗

Among proposed mechanisms, environment-dependent clock rates (as in the TEP framework) provide a plausible explanation for this differential signature.

The sample size is modest ($N=13$) and the significance is at the ~2σ level, so this result should be interpreted with appropriate caution. However, it represents a qualitatively different type of evidence than the $H_0$–$\sigma$ correlation alone, as it directly tests the *mechanism*: periodic indicators (clocks) would be biased while non-periodic indicators (thermodynamic thresholds) would not. If confirmed with larger samples, this would be the signature of a "time" effect, not a "light" effect.

#### 3.7.4b Two-Effect Decomposition

A comparative analysis shows that Cepheids exhibit a significant $H_0$–$\sigma$ correlation (Spearman $\rho = 0.434$, $p = 0.019$; $N=29$), while the TRGB sample shows a weaker, not formally significant trend (Spearman $\rho = 0.375$, $p = 0.126$; $N=18$). This pattern is consistent with two superimposed effects:

    - **Common effect:** Peculiar velocities correlate with host mass/$\sigma$, biasing $H_0$ upward in high-$\sigma$ hosts for *all* distance indicators. This is a known systematic in local distance ladder measurements.

    - **Cepheid-specific effect (TEP):** Period contraction in high-$\sigma$ environments provides an *additional* bias unique to periodic indicators.

The differential test ($\Delta\mu = \mu_{\rm TRGB} - \mu_{\rm Cepheid}$) is intended to reduce sensitivity to systematics that shift both indicators in the same direction. The positive correlation ($r = 0.55$) in the differential is consistent with the possibility that Cepheids experience an *additional* distance underestimation beyond any effect shared with TRGB, as expected if a period-dependent mechanism contributes.

#### 3.7.5 Implications for the Hubble Tension

The CCHP reports $H_0^{\rm TRGB} = 69.8 \pm 1.6$ km/s/Mpc—intermediate between the SH0ES Cepheid value ($73.0$) and Planck ($67.4$). Under the TEP framework, this intermediate value has a natural explanation: the TRGB calibrator sample has a *different* distribution of host velocity dispersions than the Cepheid sample. If the TRGB hosts are systematically lower-$\sigma$ (shallower potentials), their Cepheid-calibrated distances would be less biased, yielding an $H_0$ closer to the true value.

A discriminating test would stratify the TRGB host sample by $\sigma$ and check for a *weaker* environmental correlation than Cepheids—consistent with differential sensitivity as expected. The CCHP's intermediate $H_0$ value ($69.8$ vs. SH0ES $73.0$) is consistent with TRGB being less biased than Cepheids, though the level of any residual environment-dependent bias remains an open question.

To investigate whether crowding artifacts could be eliminated with higher resolution, Cepheids in M31 were analyzed using HST photometry from Kodric et al. (2018, J/ApJ/864/59). The HST J/H band analysis ($N_{\rm inner}=78$, $N_{\rm outer}=69$) yields:

> 
    Result: $\Delta W = +0.68 \pm 0.19$ mag (Inner Fainter), significant at 3.6σ. The signal shows a continuous radial gradient (Pearson $r = -0.16$, $p = 0.0014$) and survives all photometric quality cuts.

    Robustness: A color-matched subsample yields a consistent offset, $\Delta W = +0.62 \pm 0.15$ mag ($N_{\rm matched}=73$).

    **Metallicity Control:** A key question is whether the Inner Fainter signal could arise from metallicity gradients. The observed J−H color gradient shows Inner Cepheids are *redder* ($r = -0.25$, $p < 10^{-6}$). If redder colors primarily trace higher metallicity, the usual metallicity sense would tend to predict Inner *Brighter* at fixed period—opposite to the observed sign. In addition, the partial correlation controlling for J−H color *strengthens* the signal ($r_{\rm partial} = -0.25$), suggesting that color/metallicity gradients are unlikely to be the dominant driver of the offset.

    M31 therefore provides *supportive evidence* for environmental P-L dependence consistent with TEP screening, complementing the primary H$_0$–σ correlation in SH0ES hosts.

    ![Multi-panel synthesis showing Inner Fainter offsets in both ground-based and HST M31 data consistent with TEP screening; LMC control shows no large offset indicating no large pipeline artifacts](public/figures/robustness_synthesis_plot.png)
    Figure 7: Synthesis of environmental differential tests. Both ground-based and HST M31 data show 'Inner Fainter' offsets consistent with TEP screening (inner bulge more screened → less period contraction). The LMC control shows no large offset, suggesting the pipeline does not introduce large geometric artifacts.

### 3.9 The Density-Potential Resolution

A key physical insight resolves the apparent contradiction between the global $H_0$–$\sigma$ trend (where high $\sigma$ implies inflated $H_0$) and the M31 Inner result (where high $\sigma$ implies fainter/standard Cepheids). The TEP effect is driven by Potential Depth ($\sigma$) but gated by Local Density ($\rho$).

| Regime | Target | Structure | Potential ($\sigma$) | Density ($\rho$) | State | Outcome |
| --- | --- | --- | --- | --- | --- | --- |
| Global Trend | SN Ia Hosts | Star-forming Disks | High (50–240 km/s) | Low ($\ll \rho_{\rm trans}$) | Unscreened | Field Active → Period Contraction → High $H_0$ |
| Local Anomaly | M31 Inner | Central Bulge | High (~160 km/s) | High ($> \rho_{\rm trans}$) | Screened | Field Suppressed → Standard Clock → Fainter (Standard) |

For SN hosts like NGC 3147 ($\sigma \approx 238$ km/s), Cepheids reside in the diffuse disk. The field remains active, so the deep potential drives a large period contraction, inflating $H_0$. In M31, the "Inner" sample probes the bulge and approaches the effective galactic screening threshold $\rho_{\rm trans}\approx 0.5\,M_\odot/\mathrm{pc}^3$. Quantitatively, in the Kodric ground-based sample the mean inner density is $\bar{\rho}_{\rm in}=0.31\,M_\odot/\mathrm{pc}^3$, with $14/153$ ($\approx 9.2\%$) Inner Cepheids above $\rho_{\rm trans}$. In the Inner core ($R<1$ kpc; $N=5$), the mean density is $\bar{\rho}=2.16\,M_\odot/\mathrm{pc}^3$ and $5/5$ lie above $\rho_{\rm trans}$. Relative to the unscreened outer disk ($\bar{\rho}_{\rm out}=0.006\,M_\odot/\mathrm{pc}^3$; $0/919$ above threshold), the screened core behaves as a standard reference, yielding the observed "Inner Fainter" inversion. Thus, the M31 result is consistent with a density-gated screening transition rather than contradicting the global $H_0$–$\sigma$ trend.

                
                

                    
## 4. Discussion

    
    
### 4.1 The Nature of the Hubble Tension

    If the correlation reported here reflects a genuine physical effect, the Hubble Tension may not represent a cosmological crisis requiring new early-universe physics. Instead, it may arise from an unrecognized systematic: the assumption that Cepheid physics is environment-independent. Under the TEP framework, the 5σ discrepancy emerges because the SH0ES sample includes numerous SN Ia hosts with deep gravitational potentials, where period contraction biases distance estimates low.

    
    The correlation detected (Spearman $\rho = 0.434$, $p = 0.019$) between host velocity dispersion and derived $H_0$ is notable for an astrophysical systematic. The signal is not contingent on the aperture homogenization: the Pearson correlation is comparable when using the raw literature values ($r_{\rm raw} \approx 0.43$, $p \approx 0.02$) versus aperture-corrected values ($r_{\rm corr} \approx 0.43$, $p \approx 0.02$). Furthermore, the correlation coefficient persists in the "Stellar-Only" verification subsample ($N=17, r \approx 0.42$), though formal significance is limited by sample size ($p \approx 0.09$). Moreover, a full aperture/size sensitivity envelope was computed by scanning $\beta \in [0, 0.08]$ and scaling the effective radii by $R_{\rm eff}\times[0.7, 1.3]$, yielding stable correlations ($r \in [0.423, 0.432]$) and $\Delta H_0$ values across the entire envelope. Repeating the full $\alpha$ optimization across the same envelope gives consistent ranges ($\alpha \in [0.574, 0.584]$, $H_0^{\rm TEP} \in [68.28, 69.23]$ km/s/Mpc), i.e. a systematic envelope that is smaller than the bootstrap uncertainty ($\pm 1.51$ km/s/Mpc), indicating that the main inference does not rely on fine-tuned aperture assumptions. This reduces the concern that the result is an artifact of mixing fiber and slit measurements or sampling different galactic regions, though systematic spectroscopic follow-up would strengthen this conclusion.

    
    
        ![Sensitivity analysis showing H0-sigma correlation is stable across both raw literature velocity dispersions and aperture-corrected values, with consistent results under full aperture/size sensitivity envelope](public/figures/aperture_sensitivity.png)
        Figure 8: Sensitivity of the correlation to aperture corrections. The signal ($H_0$-$\sigma$ correlation) is present in both raw literature velocity dispersions and aperture-corrected values. The measured correlation is essentially unchanged under the correction, and remains stable under a full aperture/size sensitivity envelope.
    

    
### 4.2 Astrophysical Systematics and Confounders

    An important question is whether the observed $H_0$–$\sigma$ correlation arises from conventional astrophysical differences between low- and high-mass galaxies rather than a time-dilation effect. Specifically, high-$\sigma$ (massive) galaxies might host younger Cepheid populations (different Period-Age relations) or possess different dust properties (extinction laws).

    
    To address this, a detailed multivariate regression analysis was performed controlling for these potential confounders:

    
        - **Cepheid Age (Period-Luminosity-Age):** A strong positive correlation exists between host velocity dispersion and mean Cepheid period, indicating that massive hosts contain longer-period (younger) Cepheids. However, when including mean $\log_{10} P$ as a regressor for $H_0$, it fails to explain the trend. The coefficient for $\sigma$ remains stable.

        - **Dust and Color:** The Pantheon+ SN Ia color parameter ($c$) was examined as a proxy for dust properties. A weak correlation exists, but adding $c$ to the regression does not diminish the $\sigma$ signal.

        - **Stellar Mass:** While host mass is correlated with $\sigma$ (Faber-Jackson relation), $\sigma$ remains a robust predictor. In a full multivariate model including $\sigma$, age, dust, and mass, the velocity dispersion coefficient remains positive and comparable to mass, indicating the signal is kinematic rather than driven by age or dust systematics.

    

    
        ![Bar chart of standardized regression coefficients showing velocity dispersion (Potential) remains the dominant predictor of H0 even when controlling for Age (Mean Period), Dust (Color), and Stellar Mass](public/figures/multivariate_robustness.png)
        Figure 9: Standardized regression coefficients for $H_0$ determinants. The dependence on velocity dispersion (Potential) remains robust across model specifications, even when controlling for Age (Mean Period), Dust (Color), and Stellar Mass.
    

    This analysis suggests that the correlation is not primarily driven by population age differences or dust extinction laws. The signal appears to be kinematic in nature, consistent with the gravitational potential dependence predicted by TEP.

    
    Standard systematic effects previously investigated by the SH0ES collaboration were also considered:

    
        - **Metallicity:** The bivariate analysis (Section 3.2) indicates metallicity is not the primary driver.

        - **Crowding:** Recent JWST observations (Riess et al. 2024) limit crowding effects to < 0.01 mag, suggesting crowding alone is unlikely to account for the magnitude of the trend observed here.

    

    
### 4.3 Alternative Distance Indicators

    The Chicago-Carnegie Hubble Program (Freedman et al. 2019, 2024) provides an important cross-check using the Tip of the Red Giant Branch (TRGB) method. Their latest JWST-based measurement yields $H_0 = 69.8 \pm 1.6$ km/s/Mpc—intermediate between Cepheid and CMB values. Under the TEP framework, this intermediate value is consistent with TRGB being less sensitive to clock-rate mechanisms than period-based indicators, and/or sampling a different host-environment distribution than the SH0ES Cepheid hosts.

    
    Other distance indicators warrant investigation:

    
        - JAGB stars: Carbon-rich asymptotic giant branch stars show promise as standardizable candles (Lee et al. 2024).

        - Mira variables: Long-period variables with P-L relations; TEP predicts similar environmental bias.

        - Surface brightness fluctuations: A geometric method that should be TEP-independent.

    

    
### 4.4 Comparison with Cosmological Solutions

    Numerous cosmological solutions to the Hubble Tension have been proposed (see Di Valentino et al. 2021; Abdalla et al. 2022 for comprehensive reviews):

    
        - Early Dark Energy: An additional energy component that decays before recombination, shifting the sound horizon (Poulin et al. 2019). Requires fine-tuning and faces constraints from other CMB observables.

        - Additional Relativistic Species: Extra neutrino-like particles ($\Delta N_{\rm eff}$) that increase $H_0$ inference from the CMB. Constrained by Big Bang Nucleosynthesis.

        - Modified Gravity: Alterations to GR at cosmological scales. Generally constrained by gravitational wave observations (Abbott et al. 2017).

        - Interacting Dark Energy: Coupling between dark energy and dark matter that modifies late-time expansion.

    
    
    The TEP framework offers a distinct perspective: it locates the issue in the local measurements rather than in new early-universe physics. This has the advantage of preserving the well-tested $\Lambda$CDM model at high redshift. Moreover, TEP makes specific, testable predictions:

    
        - The bias correlates specifically with gravitational potential depth, not other galaxy properties.

        - Low-$\sigma$ hosts should independently give Planck-consistent $H_0$ (supported: $67.8 \pm 1.6$, within $1\sigma$ of Planck).

        - The correction parameter $\alpha$ should be consistent with TEP predictions from independent observations (e.g., pulsar timing).

    

    
### 4.5 Implications for the Distance Ladder

    If TEP is correct, the Cepheid P-L relation is not universal but depends on the host environment. This has immediate implications:

    
        - Future $H_0$ measurements should stratify samples by host potential depth and apply appropriate corrections.

        - The "inverse distance ladder" (using baryon acoustic oscillations and supernovae without Cepheids) provides an independent check, as it bypasses the environmental bias entirely.

    

    
### 4.6 Connection to the TEP Framework: Group Halo Screening

    The optimal coupling $\alpha = 0.58 \pm 0.16$ derived from the Hubble Tension analysis provides an independent calibration of the TEP conformal factor. This value is consistent within uncertainties with the coupling strength inferred from globular cluster pulsar spin-down rates ($\alpha \approx 0.8$, Paper 11) and the Universal Critical Density framework (Paper 7). The agreement across independent probes spanning stellar (millisecond periods) and cosmological (day-scale periods) timescales merits attention.

    
    
#### 4.6.1 Resolving the Anchor Tension via Environmental Screening

    A central puzzle in Section 3.5 is why the geometric anchors (NGC 4258, M31, LMC) show no significant $\sigma$-dependence ($\alpha_{\rm anchor} \approx 0$), while the SN Ia hosts exhibit a strong correlation ($\alpha_{\rm host} \approx 0.58$). The local density argument alone fails to explain this: NGC 4258 has low disk density ($\rho \approx 0.03\,M_\odot/\text{pc}^3 \ll \rho_{\rm trans}$) yet shows no TEP bias.

    
    A plausible resolution is group-scale dark matter halo screening. In scalar-tensor theories with chameleon or symmetron screening, the scalar field can be suppressed not only by high local baryon density but also by the *ambient gravitational potential* of the surrounding environment. A galaxy embedded in a massive group halo sits in a deeper total potential well, which may trigger screening even if the galaxy's internal disk density is low.

    
    
        
#### The Group Screening Hypothesis

        The TEP effect is gated by *two* environmental factors:

        
            - **Local density ($\rho$):** High baryon density suppresses scalar gradients (as in the M31 bulge).

            - **Group halo potential ($\Phi_{\rm group}$):** Membership in a massive group/cluster embeds the galaxy in a deep ambient potential that triggers chameleon-type screening, even if the local disk is diffuse.

        
        Either condition can suppress the TEP effect; both must be absent for the field to remain active.

    
    
    
#### 4.6.2 Application to the Anchors

    This framework naturally explains the anchor stability:

    
    

| Anchor | $\sigma$ (km/s) | Local $\rho$ | Group Environment | Screening Status |
| --- | --- | --- | --- | --- |
| LMC | 24 | Low | Local Group (MW satellite) | Screened by Local Group halo |
| NGC 4258 | 115 | Low ($\sim 0.03$) | Canes Venatici I Group ($N_{\rm mb} \approx 10$) | Screened by group halo potential |
| M31 | 160 | Transition (bulge) | Local Group (dominant member) | Screened by Local Group halo |

    
    All three anchors are group members. The Local Group potential ($M_{\rm vir} \sim 2 \times 10^{12}\,M_\odot$) and Canes Venatici I potential provide the ambient screening that suppresses the TEP effect, regardless of their internal disk densities. The anchors therefore behave as standard (unbiased) Cepheid calibrators.

    
    
#### 4.6.3 Application to SN Ia Hosts

    In contrast, SN Ia host galaxies are selected for *smooth Hubble flow*—specifically, environments where peculiar velocities are minimized and flow-model residuals are small. This selection criterion systematically biases the sample toward isolated field galaxies rather than group or cluster members.

    
    Field galaxies lack a surrounding group halo potential. Combined with their typically low disk densities ($\bar{\rho} \approx 0.1\,M_\odot/\text{pc}^3$), these hosts are *doubly unscreened*: neither local density nor ambient potential triggers field suppression. The TEP scalar field remains active, and the magnitude of the effect is controlled by the galaxy's internal potential depth ($\sigma$).

    
    
        
#### Falsifiable Prediction

        Prediction: The TEP distance-ladder bias is expected to be most prominent in isolated field galaxies and suppressed in group/cluster environments.

        This transforms the observed $N_{\rm mb}$ partial correlation (Section 3.6) from a statistical nuisance into the theory's sharpest prediction: controlling for group richness should reduce the $H_0$–$\sigma$ signal *because group membership mediates screening*, not because it confounds the measurement.

    
    
    
#### 4.6.4 Quantitative Consistency

    The robustness analysis (Section 3.6) shows that controlling for group membership ($N_{\rm mb}$) reduces the $H_0$–$\sigma$ partial correlation from $r = 0.38$ to $r = 0.32$ ($p \approx 0.11$). Under the group-screening hypothesis, this is the *expected* behavior: $N_{\rm mb}$ is not a confounding nuisance but a *mediating* variable. Galaxies in rich groups experience screening and contribute less to the overall $H_0$–$\sigma$ trend.

    
    This interpretation is supported by the observation that the SH0ES host sample is biased toward low-$N_{\rm mb}$ (field) galaxies relative to the anchor calibrators, consistent with the Hubble-flow selection criterion favoring isolated environments.

    
    
#### 4.6.5 Cross-Scale Consistency

    The coupling values derived from independent probes converge:

	
		- Globular cluster pulsar timing: $\alpha \approx 0.8$ (Paper 11)

		- The Universal Critical Density scaling (core saturation density): $\rho_c \approx 20$ g/cm³ (Paper 7)

        - Hubble Tension (field galaxies): $\alpha = 0.58 \pm 0.16$ (this work)

	
	This cross-scale agreement is consistent with the possibility that TEP provides a unified explanation for apparent anomalies across stellar and cosmological scales, with environmental screening (both local density and group halo) modulating where the effect is active.

    
### 4.7 Caveats and Limitations

    	Several caveats warrant discussion:

    	
    		- Sample size: This analysis uses $N=29$ host galaxies. Despite this modest sample size, the detection is statistically significant (Spearman $\rho \approx 0.43$, $p < 0.05$). Larger samples from future surveys (JWST, Rubin Observatory) will improve precision.

    		- **Anchor Tension (Resolved):** The geometric anchors (LMC, NGC 4258, M31) do not exhibit the strong $\sigma$-dependence seen in the SN Ia hosts. As discussed in Section 4.6, this is naturally explained by *group halo screening*: all three anchors are members of galaxy groups (Local Group for LMC and M31; Canes Venatici I for NGC 4258), embedding them in deep ambient potentials that trigger chameleon-type screening regardless of their internal disk densities. The SN Ia hosts, selected for smooth Hubble flow, are biased toward isolated field galaxies that lack this external screening.

    		- Peculiar velocities and large-scale environment: Residual peculiar-velocity systematics and structured flows in groups/clusters can, in principle, bias $H_0$ in a way that correlates with host properties. This concern is addressed directly in the robustness suite by (i) raising the redshift threshold, (ii) computing partial correlations controlling for $z_{\rm HD}$ and a group-environment proxy ($N_{\rm mb}$), and (iii) propagating Pantheon+ peculiar-velocity uncertainties. The correlation remains positive after these controls.

    		- Distance-modulus covariance: Because SH0ES host distance moduli are derived from a global GLS solution, the inferred host-level $\mu_i$ values share calibration covariance. The full GLS covariance submatrix for $\mu_i$ is propagated into a covariance matrix for the derived $H_{0,i}$ values, and the significance of the $H_0$–$\sigma$ correlation is recomputed under a correlated-null Monte Carlo model (Section 2.7). The detection remains significant under this covariance-aware treatment ($p_{\rm cov} \approx 0.02-0.04$).

    		- Potential overfitting of $\alpha$: Optimizing $\alpha$ to remove the observed $H_0$–$\sigma$ slope could in principle reparameterize in-sample noise. To test this directly, out-of-sample validation is performed (Section 2.8). Repeated 70/30 train/test splits and LOOCV demonstrate that $\alpha$ inferred on one subset predicts a reduced environmental trend and a Planck-consistent mean on held-out hosts.

    		- Velocity dispersion uncertainties: Literature $\sigma$ values have heterogeneous provenance. Systematic spectroscopic follow-up of all SH0ES hosts would strengthen the analysis.

    		- Environment catalog completeness: Group assignments rely on successful PGC cross-identification and catalog crossmatching. The primary robustness control uses $N_{\rm mb}$, which is broadly available.

    		- **Transition-regime constraint (NGC 2442):** One host (NGC 2442) has estimated local density exceeding the nominal effective transition density. Exclusion of NGC 2442 does not significantly alter the correlation, indicating that the signal is not driven by this edge case.

    		- Robustness: Stability has been verified via sensitivity analysis against variations in the calibrator reference $\sigma_{\rm ref}$, suggesting the results are not fine-tuned.

    		- Alternative proxies: $\sigma$ is used as a potential depth proxy. Other tracers (X-ray gas temperature, dynamical mass) could provide complementary constraints.

    	

    
### 4.8 Direct Test: Differential Analysis in M31

    To rigorously test the environmental dependence of the P-L relation while eliminating galaxy-to-galaxy systematics, a differential analysis of Cepheids in M31 (Andromeda) was performed using both ground-based (Kodric et al. 2018) and space-based (HST/PHAT) catalogs.

    
#### Ground-Based Signal (Crowding Dominated)

    The ground-based analysis ($N=1072$) comparing "Inner" ($R < 5$ kpc) versus "Outer" ($R > 15$ kpc) Cepheids reveals a statistically significant offset where Inner Cepheids appear systematically *fainter* ($\Delta W \approx +0.36$ mag) than their outer counterparts. However, matched-subsample tests indicate this signal is unstable against photometric quality cuts, suggesting it is driven by severe crowding in the inner bulge which biases background estimates and blending.

    
#### Space-Based Resolution (M31 HST)

    The HST J/H band analysis from Kodric et al. (2018) ($N_{\rm inner}=78$, $N_{\rm outer}=69$) shows Inner Fainter ($\Delta W = +0.68 \pm 0.19$ mag, 3.6σ). A color-matched subsample ($N_{\rm matched}=73$) yields a consistent offset of $\Delta W = +0.62 \pm 0.15$ mag, confirming the signal is not driven by metallicity differences.

    
    
        - **Metallicity/Color Control:** The persistence of the signal in a color-matched subsample indicates that the Inner Fainter offset is not driven by simple color/metallicity differences between the inner and outer samples.

        - **Interpretation:** This sign is consistent with TEP screening: the M31 inner sample lies near the effective transition regime, with the central kpc reaching densities above $\rho_{\rm trans}$, suppressing the TEP effect relative to the low-density outer disk where the field remains active.

        - **Robustness:** The Inner Fainter offset remains consistent under color matching ($N_{\rm matched}=73$), supporting an environmental interpretation rather than a selection artifact.

        - **Implication:** M31 provides supportive evidence for environmental P-L dependence consistent with TEP screening, complementing the primary H$_0$–$\sigma$ correlation.

    

    The synthesis of these environmental tests is visualized in Figure 7 (see Section 3.8).

    
#### Density Regimes and Screening Resolution

    The M31 result, initially appearing contradictory to the "High $\sigma$ = High Effect" rule, is resolved by considering the local density environment relative to the screening threshold. This requires distinguishing two density scales discussed in Paper 7:

    
    
        - **Universal Critical Density ($\rho_c \approx 20$ g/cm³):** The fundamental core saturation density of the scalar sector.

        - **Effective Transition Density ($\rho_{\rm trans} \approx 0.5 \, M_\odot/\text{pc}^3$):** An emergent screening threshold for galactic structures.

    
    
    The resolution:

    
        - **Global Trend (Unscreened Disks):** The SN Ia host sample typically has low local densities ($\bar{\rho} \approx 0.113 \, M_\odot/\text{pc}^3 \ll \rho_{\rm trans}$). In this *unscreened* regime, the TEP field is active. Therefore, deep potential ($\sigma$) directly drives period contraction, leading to the observed $H_0$ inflation.

        - **M31 Anomaly (Transition/Suppressed Bulge):** The M31 "Inner" sample probes the bulge-dominated region ($\bar{\rho} \approx 0.309 \, M_\odot/\text{pc}^3$), which is near the effective transition density. The density rises steeply toward the nucleus, exceeding $\rho_{\rm trans}$ at $R \lesssim 1$ kpc. In this transition-to-screened regime, the field is suppressed. Relative to the active (brightened) outer disk, the inner region appears fainter.

    

    
        
#### Quantitative Threshold Verification

        Is the transition density $\rho_{\rm trans}$ tuned to fit M31? No—it is derived independently from the SPARC rotation curve database (Paper 7) as the galactic-scale manifestation of the fundamental saturation density $\rho_c \approx 20$ g/cm³. The galaxy scaling $R_{\rm DM} \propto M^{1/3}$ normalizes to an effective screening density of $\rho_{\rm trans} \approx 0.5\,M_\odot/\text{pc}^3$. This independent threshold is explicitly compared to the study environments:

        
        
            **SN Ia Hosts (Unscreened):** Typical spiral disks at the optical radius ($R_{25}$) have mean stellar densities of $\bar{\rho} \approx 0.1\text{--}0.2\,M_\odot/\text{pc}^3$.
            **$\rightarrow \rho_{\rm host} < \rho_{\rm trans}$ implies TEP Active** (Period contraction $\rightarrow H_0$ bias).
            
            **M31 Inner Bulge (Transition/Screened Core):** The "Inner" sample probes $R < 5$ kpc with a mean local density of $\bar{\rho} \approx 0.31\,M_\odot/\text{pc}^3$. In the Kodric ground-based sample, $14/153$ Inner Cepheids ($\approx 9.2\%$) lie above the effective galactic screening threshold $\rho_{\rm trans}=0.5\,M_\odot/\text{pc}^3$. In the Inner core ($R<1$ kpc; $N=5$), the mean density is $\bar{\rho}\approx 2.16\,M_\odot/\text{pc}^3$ and $5/5$ lie above $\rho_{\rm trans}$.
            $\rightarrow$ The data therefore directly sample both the unscreened disk and a screened/transition bulge core, consistent with density-gated suppression.
        
        
        The "Inner Fainter" signal is therefore consistent with the sample crossing the SPARC-derived density threshold, rather than requiring a post-hoc tuning of $\rho_{\rm trans}$.

    

    This result highlights that environmental calibration may require accounting for both the background potential $\Phi$ (which sets the magnitude of the effect) and the local screening density $\rho$ (which gates the regime). In this interpretation, the "Inner Fainter" signal is consistent with the screening threshold being approached or crossed.

    
### 4.9 Falsifiable Predictions for Alternative Distance Indicators

    The TEP framework makes explicit, testable predictions for how different distance indicators should depend on host environment. These predictions follow directly from the microphysics: indicators that rely on periodic phenomena (clocks) should show environmental bias proportional to their period-luminosity coupling, while geometric or non-periodic indicators should be unaffected.

    
    
    

| Indicator | Mechanism | TEP Prediction | Expected $H_0$–$\sigma$ Slope |
| --- | --- | --- | --- |
| Cepheids | Period-luminosity (P-L) | Strong positive bias | $dH_0/d\log_{10}\sigma \approx +15$–$25$ km/s/Mpc/dex |
| Mira Variables | Period-luminosity (long-period) | Positive bias (similar to Cepheids) | $dH_0/d\log_{10}\sigma \approx +10$–$20$ km/s/Mpc/dex |
| RR Lyrae | Period-luminosity (short-period) | Positive bias (weaker due to shorter periods) | $dH_0/d\log_{10}\sigma \approx +5$–$15$ km/s/Mpc/dex |
| TRGB | Luminosity threshold (no period) | Weak or absent | $dH_0/d\log_{10}\sigma \approx 0$ |
| SBF | Stellar fluctuations (geometric) | Weak or absent | $dH_0/d\log_{10}\sigma \approx 0$ |
| JAGB | Luminosity function (no period) | Weak or absent | $dH_0/d\log_{10}\sigma \approx 0$ |
| Megamasers | Pure geometry | Absent | $dH_0/d\log_{10}\sigma = 0$ |

    

    
### 4.9.1 The Clock vs. Non-Clock Differential Test

    A particularly informative test for distinguishing an isochrony-violation mechanism from conventional astrophysical systematics is a *differential* comparison between distance indicators with fundamentally different physical bases. This section highlights the key logic by comparing Cepheid distances (which depend on periodic timekeeping) with TRGB distances (which depend on nuclear physics thresholds).

    
#### 4.9.1.1 The "Clock" vs "Light" Distinction

    Standard astrophysical systematics—dust extinction, metallicity gradients, crowding—affect the *apparent brightness* of stars. These are "light" effects: they modify how many photons reach the observer, and in the simplest picture they should act similarly on multiple tracers within comparable regions of the same host.

    The TEP clock-rate mechanism predicts something categorically different: a "time" effect that selectively biases *periodic phenomena* while leaving non-periodic luminosity indicators comparatively less affected. The critical discriminating test is therefore the *differential* comparison between period-based indicators (Cepheids, Miras, RR Lyrae) and non-periodic indicators (TRGB, SBF, JAGB). Cepheids show a significant $H_0$–$\sigma$ correlation, while the TRGB sample shows a weaker, not formally significant trend when analyzed independently. This pattern is consistent with two superimposed effects: (1) a common systematic (peculiar velocity–mass correlations) affecting all indicators, and (2) a Cepheid-specific bias (TEP period contraction). The differential analysis in Section 3.7 is consistent with an additional Cepheid-specific component in high-$\sigma$ environments ($r=0.55$, $p=0.05$), as expected if a period-dependent mechanism contributes.

    
    Preliminary evidence is consistent with this prediction: the Chicago-Carnegie Hubble Program (Freedman et al. 2024) reports $H_0 = 69.8 \pm 1.6$ km/s/Mpc from TRGB—intermediate between Cepheid and CMB values. If TRGB is less sensitive to TEP effects, this intermediate value may reflect a smaller environmental bias in the TRGB calibrator sample compared to the SH0ES Cepheid sample. A discriminating test would stratify TRGB host galaxies by $\sigma$ and check for a weaker environmental correlation than Cepheids.

    
#### The "SN Ia Clock" Systematics

    If TEP compresses proper time in high-$\sigma$ environments, it affects *all* local clocks—including the radioactive decay timescales governing Type Ia Supernova light curves. Since SN Ia standardization relies on width-luminosity relations (e.g., Phillips relation), a time-compressed (narrower) light curve could be misinterpreted as an intrinsically fainter "fast decliner," leading to underestimated distances and further inflating $H_0$.

    
    Does this effect rival the Cepheid bias? While it technically compounds the tension (biasing distances in the same direction), its magnitude is negligible compared to the Cepheid zero-point shift. The sensitivity of the derived magnitude to time dilation is determined by the slope of the calibration relation:

    
    
        - **Cepheids (Leavitt Law):** Slope $dM/d\log P \approx -2.4$. A 1% period change drives a large luminosity error.

        - **SN Ia (Width-Luminosity):** The sensitivity parameter is $\alpha \approx 0.14$ (SALT2).

    
    
    Because the Cepheid P-L relation is nearly an order of magnitude steeper than the SN Ia width-luminosity relation, the Cepheid calibration bias dominates the error budget. The SN Ia clock effect constitutes a small second-order term.

    
### 4.10 Future Observational Tests

    Several observational programs can further validate or falsify the TEP explanation:

    
        - **IFS "Gold Standard" Analysis:** Integral Field Spectroscopy (IFS) from MaNGA or CALIFA can provide spatially resolved velocity dispersions at a consistent physical radius (e.g., 1 kpc) for a subset of SH0ES hosts. Even a small ($N \sim 10$) homogeneous subsample supporting the $H_0$–$\sigma$ correlation would strongly constrain aperture systematics as a potential explanation.

        - **JWST Cepheid Observations:** Targeted observations of Cepheids in a controlled sample spanning a wide $\sigma$ range, with homogeneous photometry and metallicity corrections.

        - **TRGB Stratification:** Stratifying existing TRGB distance measurements by host $\sigma$ to test for the predicted weaker environmental correlation relative to Cepheids.

        - **M31 Homogeneous Reanalysis:** A differential P-L analysis using a photometrically homogeneous Cepheid subset (matched crowding/extinction) to isolate the environmental signal from selection effects.

        - **Laboratory Clock Experiments:** Precision tests of optical clocks at different altitudes or in variable gravitational environments.

    
    

                
                

                    
## 5. Conclusion

This work investigates whether the Hubble Tension—a persistent challenge in precision cosmology—might arise from an environmental systematic in Cepheid-based distances. The key findings are:

    - Environmental bias detected: Stratification of the SH0ES Cepheid host galaxies by directly measured velocity dispersion shows a significant correlation (Spearman $\rho = 0.434$, $p = 0.019$) between host potential depth and derived $H_0$. High-$\sigma$ hosts yield systematically inflated $H_0$ values ($72.45 \pm 2.32$ km/s/Mpc) compared to low-$\sigma$ hosts ($67.82 \pm 1.62$ km/s/Mpc).

    
    - Magnitude contributes to tension: The bias $\Delta H_0 = 4.63$ km/s/Mpc between high and low-density hosts accounts for a substantial fraction of the discrepancy between local and CMB measurements.

    
    - TEP correction reduces the tension: Application of the TEP conformal correction with optimal coupling $\alpha = 0.58 \pm 0.16$ and effective calibrator $\sigma_{\rm ref} = 75.25$ km/s yields a unified local Hubble constant of $H_0 = 68.66 \pm 1.51$ km/s/Mpc. This result is robust under bootstrap resampling and reduces the tension with Planck to $0.79\sigma$.

    
    - Internal consistency: Low-$\sigma$ hosts, which have environments similar to the calibrators, independently yield Planck-consistent $H_0$ (within $1\sigma$) without correction, consistent with TEP expectations.

    
    - Anchor consistency test ("anchor tension" resolved): Independent P-L fits to the extragalactic geometric anchors (LMC, NGC 4258, M31; MW excluded due to its distinct parallax-based methodology) yield $\alpha_{\rm anchor} = 0.029 \pm 0.023$—consistent with zero and in 3.4σ tension with the host-level coupling. This dichotomy is naturally explained by *group halo screening*: all three anchors are members of galaxy groups (Local Group for LMC and M31; Canes Venatici I for NGC 4258), embedding them in deep ambient potentials that trigger chameleon-type screening. The SN Ia hosts, selected for smooth Hubble flow, are biased toward isolated field galaxies where the TEP effect remains active.

    - M31 screening consistency: The "Inner Fainter" signal observed in M31 is consistent with a density-gated screening transition. While the global $H_0$ trend is driven by unscreened, low-density disks (where deep potential = active TEP), the M31 inner region approaches the effective galactic screening threshold ($\rho_{\rm trans} \approx 0.5\,M_\odot/\mathrm{pc}^3$). Quantitatively, the Kodric ground-based sample yields $\bar{\rho}_{\rm in}=0.31\,M_\odot/\mathrm{pc}^3$ with $14/153$ ($\approx 9.2\%$) Inner Cepheids above $\rho_{\rm trans}$, while the inner core ($R<1$ kpc; $N=5$) has $\bar{\rho}=2.16\,M_\odot/\mathrm{pc}^3$ and $5/5$ above threshold. This provides a physically motivated "screened core" control within the same galaxy.

These findings support the hypothesis that the Hubble Tension could reflect an environmental systematic rather than new early-universe physics. The Temporal Equivalence Principle—supported by independent studies of pulsar timing in globular clusters (Paper 11) and by the potential- and density-dependent structure identified here—provides a concrete framework for organizing these correlations and for generating falsifiable predictions.

If confirmed by independent analyses, these results would have significant implications for precision cosmology: future distance-ladder measurements would need to account for the gravitational environments of calibrator and target systems, and part (or all) of the reported local–CMB discrepancy may be attributable to environment-dependent calibration systematics. The findings presented here motivate targeted follow-up tests (homogeneous stellar-dispersion spectroscopy; TRGB stratification by $\sigma$; JWST Cepheid imaging) to more directly validate or falsify the proposed mechanism.

> 
    
### Code and Data Availability

    All analysis code is open-source and designed for easy reproduction. The complete pipeline runs in under 2 minutes and reproduces all results, figures, and statistics reported in this paper.

    
    
#### Quick Start

    To reproduce the full analysis:

    # Clone the repository
git clone https://github.com/matthewsmawfield/TEP-H0.git
cd TEP-H0

# Install dependencies
pip install -r requirements.txt

# Run the complete analysis pipeline
python scripts/run_pipeline.py
    
    The pipeline automatically downloads all required data from public archives (SH0ES, Pantheon+, HyperLEDA, Vizier) and generates all outputs.

    
    
#### Pipeline Architecture

    The analysis is organized into 10 sequential steps, each implemented as a self-contained Python module:

    
    
    

| Step | Script | Description | Key Outputs |
| --- | --- | --- | --- |
| 1 | step_1_data_ingestion.py | Downloads SH0ES distance moduli and Pantheon+ redshifts; cross-matches hosts with velocity dispersion catalogs (HyperLEDA, SDSS) | hosts_processed.csv |
| 1b | step_1b_aperture_correction.py | Applies Jorgensen et al. (1995) aperture corrections to normalize $\sigma$ measurements to $R_{\rm eff}/8$ | Homogenized $\sigma$ values |
| 2 | step_2_stratification.py | Calculates per-host $H_0$; stratifies by median $\sigma$; computes correlation statistics | stratification_results.json |
| 3 | step_3_tep_correction.py | Optimizes $\alpha$ by minimizing residual $H_0$–$\sigma$ slope; applies TEP correction; bootstrap uncertainty estimation | tep_correction_results.json |
| 4 | step_4_robustness_checks.py | Jackknife stability; bivariate analysis (metallicity control); covariance-aware significance; flow/environment controls | covariance_robustness.json |
| 5 | step_5_m31_analysis.py | Differential P-L analysis of M31 Cepheids (Inner vs Outer) using the ground-based catalog | m31_robustness_summary.json |
| 6 | step_6_multivariate_analysis.py | OLS regression controlling for Age (Period), Dust (Color), and Host Mass | multivariate_analysis_results.json |
| 7 | step_7_lmc_replication.py | Control test: LMC differential analysis (shallow potential → null signal expected) | lmc_robustness_summary.json |
| 8 | step_8_m31_phat_analysis.py | HST J/H band analysis from Kodric et al. (2018); metallicity control via color matching | m31_phat_robustness_summary.json |
| 9 | step_9_final_synthesis.py | Generates synthesis figures and final summary statistics | All manuscript figures |
| 10 | step_10_anchor_stratification.py | Independent P-L fits to geometric anchors (LMC, NGC 4258, M31); tests for anchor-level TEP bias | anchor_stratification_test.json |

    
    
    
#### Repository Structure

    TEP-H0/
├── scripts/
│   ├── run_pipeline.py          # Master orchestration script
│   ├── steps/                   # Individual analysis modules
│   └── utils/                   # Shared utilities (logging, plotting)
├── data/
│   ├── raw/                     # Downloaded source data
│   ├── interim/                 # Intermediate processing
│   └── processed/               # Final host catalog
├── results/
│   ├── outputs/                 # JSON/CSV results (all key statistics)
│   └── figures/                 # Generated figures (PNG)
└── site/                        # Manuscript HTML and website
    
    
#### Key Output Files

    
        - results/outputs/tep_correction_results.json — Unified $H_0$, optimal $\alpha$, Planck tension

        - results/outputs/stratification_results.json — High/low-$\sigma$ stratification statistics

        - results/outputs/covariance_robustness.json — Covariance-aware p-values and $N_{\rm eff}$

        - results/outputs/out_of_sample_validation.json — Train/test and LOOCV results

        - data/processed/hosts_processed.csv — Complete host galaxy catalog with $\sigma$, $H_0$, corrections

    
    
    
#### Dependencies

    The pipeline requires Python 3.8+ and the following packages (all installable via pip):

    
        - numpy

        - scipy

        - pandas

        - matplotlib

        - astropy

        - astroquery

    
    
    
#### Verification

    After running the pipeline, verify reproduction by checking:

    # Check key results match manuscript
cat results/outputs/tep_correction_results.json | grep unified_h0
# Expected: 68.66 (±0.01)

cat results/outputs/stratification_results.json | grep difference
# Expected: 4.63 (±0.01)
    
    
        
            
                
            
            https://github.com/matthewsmawfield/TEP-H0
        
    

    
    
        **DOI:** [10.5281/zenodo.18216583](https://doi.org/10.5281/zenodo.18216583) &nbsp;|&nbsp;
        **License:** CC BY 4.0
    

                
                

                    
## References

    
#### Primary Data Sources

    Riess, A. G., Yuan, W., Macri, L. M., et al. 2022, *ApJ*, 934, L7, "A Comprehensive Measurement of the Local Value of the Hubble Constant with 1 km/s/Mpc Uncertainty from the Hubble Space Telescope and the SH0ES Team"

    Planck Collaboration, Aghanim, N., Akrami, Y., et al. 2020, *A&A*, 641, A6, "Planck 2018 results. VI. Cosmological parameters"

    Scolnic, D., Brout, D., Carr, A., et al. 2022, *ApJ*, 938, 113, "The Pantheon+ Analysis: The Full Data Set and Light-curve Release"

    Huchra, J. P., Macri, L. M., Masters, K. L., et al. 2012, *ApJS*, 199, 26, "The 2MASS Redshift Survey—Description and Data Release"

    Tully, R. B. 2015, *AJ*, 149, 171, "Galaxy Groups: A 2MASS Catalog"

    
    
#### Geometric Calibrators

    Gaia Collaboration, Vallenari, A., Brown, A. G. A., et al. 2023, *A&A*, 674, A1, "Gaia Data Release 3: Summary of the content and survey properties"

    Pietrzyński, G., Graczyk, D., Gallenne, A., et al. 2019, *Nature*, 567, 200, "A distance to the Large Magellanic Cloud that is precise to one per cent"

    Reid, M. J., Pesce, D. W., & Riess, A. G. 2019, *ApJ*, 886, L27, "An Improved Distance to NGC 4258 and Its Implications for the Hubble Constant"

    
    
#### Astronomical Databases

    Wenger, M., Ochsenbein, F., Egret, D., et al. 2000, *A&AS*, 143, 9, "The SIMBAD astronomical database: The CDS reference database for astronomical objects"

    Ochsenbein, F., Bauer, P., & Marcout, J. 2000, *A&AS*, 143, 23, "The VizieR database of astronomical catalogues"

    Makarov, D., Prugniel, P., Terekhova, N., Courtois, H., & Vauglin, I. 2014, *A&A*, 570, A13, "HyperLEDA. III. The catalogue of extragalactic distances"

    Abazajian, K. N., Adelman-McCarthy, J. K., Agüeros, M. A., et al. 2009, *ApJS*, 182, 543, "The Seventh Data Release of the Sloan Digital Sky Survey"

    
#### Galaxy Size Catalogs

    de Vaucouleurs, G., de Vaucouleurs, A., Corwin, H. G., Jr., et al. 1991, *Third Reference Catalogue of Bright Galaxies* (RC3), Springer

    
    
#### Velocity Dispersion Data

    Ho, L. C., Greene, J. E., Filippenko, A. V., & Sargent, W. L. W. 2009, *ApJS*, 183, 1, "A Search for 'Dwarf' Seyfert Nuclei. VII. A Complete Survey of the SDSS Spectroscopic Catalog"

    Jorgensen, I., Franx, M., & Kjærgaard, P. 1995, *MNRAS*, 276, 1341, "Spectroscopy for E and S0 galaxies in nine clusters"

    Kormendy, J. & Ho, L. C. 2013, *ARA&A*, 51, 511, "Coevolution (Or Not) of Supermassive Black Holes and Host Galaxies"

    Courteau, S., Dutton, A. A., van den Bosch, F. C., et al. 2007, *ApJ*, 671, 203, "Scaling Relations of Spiral Galaxies"

    Catinella, B., Giovanelli, R., & Haynes, M. P. 2006, *ApJ*, 640, 751, "Template Rotation Curves for Disk Galaxies"

    
    
#### Cepheid Physics

    Anderson, R. I., Saio, H., Ekström, S., Georgy, C., & Meynet, G. 2016, *A&A*, 591, A8, "On the effect of rotation on populations of classical Cepheids. II. Pulsation analysis for metallicities 0.014, 0.006, and 0.002"

    Bono, G., Marconi, M., Cassisi, S., et al. 2005, *ApJ*, 621, 966, "Classical Cepheid Pulsation Models. X. The Period-Age Relation"

    Kodric, M., Riffeser, A., Seitz, S., et al. 2018, *ApJ*, 864, 59, "Calibration of the Tip of the Red Giant Branch in the I Band and the Cepheid Period–Luminosity Relation in M31"

    Leavitt, H. S. & Pickering, E. C. 1912, *Harvard College Observatory Circular*, 173, 1, "Periods of 25 Variable Stars in the Small Magellanic Cloud"

    Madore, B. F. & Freedman, W. L. 1991, *PASP*, 103, 933, "The Cepheid distance scale"

    
    
#### TEP Framework (This Series)

    Smawfield, M. L. 2025a, "Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed" (Paper 1)

    Smawfield, M. L. 2025e, "Temporal-Spatial Coupling in Gravitational Lensing" (Paper 5)

    Smawfield, M. L. 2025g, "Universal Critical Density: Unifying Atomic, Galactic, and Compact Object Scales" (Paper 7)

    Smawfield, M. L. 2025j, "What Do Precision Tests of General Relativity Actually Measure?" (Paper 10)

    Smawfield, M. L. 2026a, "The Temporal Equivalence Principle: Suppressed Density Scaling in Globular Cluster Pulsars" (Paper 11)

    
    
#### JWST Distance Ladder Studies

    Riess, A. G., Yuan, W., Casertano, S., et al. 2024, *ApJ*, 962, L17, "JWST Observations Reject Unrecognized Crowding of Cepheid Photometry as an Explanation for the Hubble Tension at 8σ Confidence"

    Freedman, W. L., Madore, B. F., Hoyt, T. J., et al. 2024, arXiv:2408.06153, "Status Report on the Chicago-Carnegie Hubble Program (CCHP): Measurement of the Hubble Constant Using the Hubble and James Webb Space Telescopes"

    Freedman, W. L., Madore, B. F., Hatt, D., et al. 2019, *ApJ*, 882, 34, "The Carnegie-Chicago Hubble Program. VIII. An Independent Determination of the Hubble Constant Based on the Tip of the Red Giant Branch"

    Lee, A. J., Freedman, W. L., Madore, B. F., et al. 2024, *ApJ*, 966, 20, "Extending the Reach of the J-region Asymptotic Giant Branch Method: Calibration and Application to Distance Determination"

    
    
#### Hubble Tension Reviews & Proposed Solutions

    Freedman, W. L. 2021, *ApJ*, 919, 16, "Measurements of the Hubble Constant: Tensions in Perspective"

    Di Valentino, E., Mena, O., Pan, S., et al. 2021, *Classical and Quantum Gravity*, 38, 153001, "In the realm of the Hubble tension—a review of solutions"

    Abdalla, E., Abellán, G. F., Aboubrahim, A., et al. 2022, *Journal of High Energy Astrophysics*, 34, 49, "Cosmology intertwined: A review of the particle physics, astrophysics, and cosmology associated with the cosmological tensions and anomalies"

    Poulin, V., Smith, T. L., Karwal, T., & Kamionkowski, M. 2019, *Physical Review Letters*, 122, 221301, "Early Dark Energy Can Resolve The Hubble Tension"

    Abbott, B. P., Abbott, R., Abbott, T. D., et al. (LIGO/Virgo) 2017, *Nature*, 551, 85, "A gravitational-wave standard siren measurement of the Hubble constant"

    
    
#### Statistical Methods

    Zahid, H. J., Geller, M. J., Fabricant, D. G., & Hwang, H. S. 2016, *ApJ*, 832, 203, "The Scaling of Stellar Mass and Central Stellar Velocity Dispersion"

                
                

                    
## Appendix A: Per-Host Data Table

Table A1 presents the complete per-host dataset used in this analysis. For each SN Ia host galaxy, the table provides: redshift ($z_{\rm HD}$), distance modulus ($\mu$), derived Hubble constant ($H_{0,i}$), raw and aperture-corrected velocity dispersions ($\sigma_{\rm raw}$, $\sigma_{\rm corr}$), the $\sigma$ measurement source, the total $\sigma$ uncertainty ($\delta\sigma$), and a host metallicity proxy ($\log_{10} M_*$), alongside the $\sigma$ measurement method classification. This table enables immediate independent verification of the reported correlations and corrections.

| Host | $z_{\rm HD}$ | $\mu$ (mag) | $H_{0,i}$ (km/s/Mpc) | $\sigma_{\rm raw}$ (km/s) | $\sigma_{\rm corr}$ (km/s) | $\sigma$ Source | $\delta\sigma$ (km/s) | $\log_{10} M_*$ | $\sigma$ Method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NGC 0691 | 0.00855 | 32.82 | 69.9 | 107.5 | 101.4 | Ho+2007 | 5.4 | 10.83 | Stellar |
| NGC 1015 | 0.00815 | 32.62 | 73.2 | 106.5 | 101.5 | HyperLEDA | 8.5 | 9.91 | Stellar |
| NGC 105 | 0.01682 | 34.49 | 63.7 | 56.7 | 55.9 | HyperLEDA | 2.8 | 10.12 | Rot. proxy |
| NGC 1309 | 0.00719 | 32.51 | 67.9 | 82.0 | 78.8 | HyperLEDA | 27.0 | 9.89 | Stellar |
| NGC 1365 | 0.00483 | 31.33 | 78.6 | 151.4 | 136.2 | Ho+2007 | 7.6 | 10.73 | Stellar |
| NGC 1559 | 0.00407 | 31.46 | 62.3 | 72.6 | 68.5 | ApJ 929 | 3.6 | 9.55 | Stellar |
| NGC 2442 | 0.00488 | 31.47 | 74.5 | 144.2 | 133.5 | HyperLEDA (HI) | 7.2 | 12.20 | HI proxy |
| NGC 2525 | 0.00602 | 32.01 | 71.5 | 86.5 | 82.2 | HyperLEDA (HI) | 4.3 | 10.06 | HI proxy |
| NGC 2608 | 0.00855 | 32.63 | 76.4 | 86.6 | 83.0 | HyperLEDA (HI) | 4.3 | 10.45 | HI proxy |
| NGC 3021 | 0.00673 | 32.39 | 67.1 | 57.3 | 55.8 | Ho+2007 | 2.9 | 10.30 | Stellar |
| NGC 3147 | 0.01079 | 33.09 | 77.9 | 219.8 | 206.3 | Ho+2009 | 14.0 | 8.37 | Stellar |
| NGC 3254 | 0.00648 | 32.40 | 64.2 | 117.8 | 109.5 | Ho+2009 | 7.2 | 10.63 | Stellar |
| NGC 3370 | 0.00588 | 32.14 | 65.7 | 94.6 | 89.5 | Ho+2009 | 10.5 | 10.20 | Stellar |
| NGC 3447 | 0.00465 | 31.94 | 56.9 | 67.8 | 63.7 | HyperLEDA (HI) | 3.4 | 9.53 | HI proxy |
| NGC 3583 | 0.00857 | 32.79 | 71.1 | 131.7 | 125.2 | Ho+2009 | 12.1 | 10.95 | Stellar |
| NGC 4038 | 0.00571 | 31.63 | 80.7 | 107.4 | 99.6 | HyperLEDA (HI) | 5.4 | 10.68 | HI proxy |
| NGC 4639 | 0.00359 | 31.79 | 47.3 | 96.0 | 91.4 | Ho+2009 | 6.2 | 9.80 | Stellar |
| NGC 4680 | 0.00864 | 32.55 | 80.2 | 102.7 | 100.3 | HyperLEDA (HI) | 5.1 | 9.75 | HI proxy |
| NGC 5468 | 0.00954 | 33.19 | 65.9 | 67.6 | 64.5 | HyperLEDA (HI) | 3.4 | 10.44 | HI proxy |
| NGC 5584 | 0.00625 | 31.87 | 79.4 | 54.2 | 51.1 | SDSS DR7 | 10.0 | 10.33 | Stellar |
| NGC 5728 | 0.00996 | 32.92 | 78.0 | 176.0 | 166.7 | BASS DR2 | 9.7 | 10.64 | Stellar |
| NGC 5861 | 0.00677 | 32.21 | 73.5 | 112.2 | 106.4 | HyperLEDA (HI) | 5.6 | 10.59 | HI proxy |
| NGC 5917 | 0.00710 | 32.34 | 72.6 | 54.5 | 53.1 | HyperLEDA | 2.7 | 9.18 | Rot. proxy |
| NGC 7250 | 0.00432 | 31.61 | 61.8 | 41.8 | 40.5 | HyperLEDA | 2.1 | 9.13 | Rot. proxy |
| NGC 7329 | 0.01028 | 33.27 | 68.4 | 123.7 | 116.1 | HyperLEDA (HI) | 6.2 | 10.50 | HI proxy |
| NGC 7541 | 0.00814 | 32.58 | 74.4 | 64.4 | 60.7 | HyperLEDA | 34.7 | 10.94 | Stellar |
| NGC 7678 | 0.01061 | 33.27 | 70.7 | 76.9 | 73.6 | SDSS DR7 | 5.4 | 10.53 | Stellar |
| NGC 976 | 0.01312 | 33.54 | 76.9 | 217.6 | 212.4 | MNRAS 482 | 21.1 | 10.85 | Stellar |
| UGC 9391 | 0.00747 | 32.82 | 61.2 | 74.5 | 72.4 | SDSS DR7 | 27.6 | 9.35 | Stellar |

*Notes:* $z_{\rm HD}$ is the Hubble-diagram redshift from Pantheon+. $\mu$ is the SH0ES distance modulus. $H_{0,i} = cz_{\rm HD}/d_i$ where $d_i = 10^{(\mu-25)/5}$ Mpc. $\sigma_{\rm raw}$ is the literature velocity dispersion; $\sigma_{\rm corr}$ is aperture-corrected to $R_{\rm eff}/8$ using Jorgensen et al. (1995). $\delta\sigma$ is the total uncertainty including measurement and aperture-correction components. $\log_{10} M_*$ is the host stellar mass from Pantheon+. $\sigma$ Method indicates whether the measurement is from stellar absorption spectroscopy (gold standard) or HI 21-cm linewidth proxy. Sources: HyperLEDA = stellar absorption unless noted (HI) for HI linewidth proxy; Ho+2009 = Ho et al. (2009); Kormendy&amp;Ho2013 = Kormendy &amp; Ho (2013); SDSS DR7 = Sloan Digital Sky Survey fiber spectroscopy.

### A.1 Velocity Dispersion Provenance

The velocity dispersion compilation draws from multiple sources with heterogeneous methodology:

    - **Stellar absorption (direct):** 17 hosts have $\sigma$ measured from stellar absorption line broadening, the gold-standard method. Sources include HyperLEDA, SDSS DR7, Ho et al. (2007, 2009), BASS DR2, and MNRAS 482:1427.

    - **HI linewidth proxy:** 9 hosts use HI 21-cm linewidth measurements calibrated via $\sigma = 0.467 \times V_{\rm max} + 42.9$ km/s (HyperLEDA calibrated_vmax mode). This introduces additional scatter but preserves the kinematic nature of the observable.

    - **Rotation velocity proxy:** 3 hosts (NGC 105, NGC 5917, NGC 7250) use rotation velocity converted via $\sigma \approx V_{\rm rot}/1.7$, a standard scaling for late-type spirals.

The correlation coefficient persists when restricting to stellar-absorption-only hosts ($N=17$, Pearson $r \approx 0.42$, $p \approx 0.09$). The reduction in formal significance relative to the full sample is expected from reduced statistical power ($\sim 41\%$ fewer degrees of freedom), not from the removal of spurious signal. Critically, the 12 kinematic-proxy hosts (9 HI + 3 rotation) do not cluster anomalously—they span the full $\sigma$–$H_0$ distribution and follow the same physical trend as stellar hosts (see Section 3.2). Application of the TEP correction to the stellar-only subsample yields $\alpha = 0.56$ and a unified $H_0 = 67.77 \pm 1.81$ km/s/Mpc, consistent with the full-sample result. Future systematic spectroscopic follow-up of all SH0ES hosts would further strengthen this analysis.

### A.2 Aperture Correction Details

Literature $\sigma$ values are measured through apertures of varying physical size. To homogenize the sample, the Jorgensen et al. (1995) power-law correction is applied:

$\sigma_{\rm corr} = \sigma_{\rm obs} \left( \frac{r_{\rm ap}}{R_{\rm eff}/8} \right)^{0.04}$

where $r_{\rm ap}$ is the observational aperture radius (assumed 1.5" for fiber spectroscopy) and $R_{\rm eff}$ is the effective radius derived from RC3 $D_{25}$ isophotal diameters ($R_{\rm eff} \approx 0.5 R_{25}$). The correction is typically $\lesssim 10\%$ and does not qualitatively affect the correlation (see Section 4.1 and the aperture sensitivity envelope in Section 3.6).

                

    
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

        - [Global Time Echoes: Optical Validation of TEP via Satellite Laser Ranging 30 Dec 2025](/tep/slr/)

        - [What Do Precision Tests of General Relativity Actually Measure? 31 Dec 2025](/tep/exp/)

        - [The Temporal Equivalence Principle: Suppressed Density Scaling in Globular Cluster Pulsars 9 Jan 2026](/tep/cos/)

        - [The Cepheid Bias: Resolving the Hubble Tension 11 Jan 2026](/tep/h0/)

    
    
    
        ← Previous
        Next →