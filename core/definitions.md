## Appendix: Canonical TEP Definitions

This appendix states the canonical definitions used throughout the TEP corpus. Every paper in the series adopts these conventions unless an explicit deviation is noted.

### A.1 Matter Metric

The TEP framework operates on a single manifold with two metrics. Gravity is described by a Lorentzian metric $g_{\mu\nu}$; matter fields, clocks, and rulers couple to a causal matter metric $\tilde{g}_{\mu\nu}$:

$$\tilde{g}_{\mu\nu} = A^2(\phi) g_{\mu\nu} + B(\phi) \nabla_\mu\phi \nabla_\nu\phi$$

where $\phi$ is the time field, $A(\phi) = \exp(\beta_A\phi/M_{\text{Pl}})$ is the universal conformal factor, and $B(\phi)$ is a small disformal function consistent with multi-messenger constraints.

### A.2 Sector Mapping

| Sector | Physical effect | Coupling | Key bound |
|---|---|---|---|
| Conformal $A(\phi)$ | Clock rates, proper-time rescaling | Universal ($\beta_A$) | Cassini PPN-$\gamma$, redshift tests |
| Disformal $B(\phi)$ | Null-cone tilts, direction-dependent propagation | Small, environment-dependent | GW170817: $\|c_\gamma - c_g\|/c \lesssim 10^{-15}$ |
| Temporal Topology | Spatial/covariance structure of $\ln A(\phi)$ | Field configuration + gradient suppression | Clock/covariance correlation length $\lambda_T$ |

The conformal sector governs clock rates; the disformal sector governs cone tilts.

### A.3 Temporal Shear

The Temporal Shear is the gradient of the conformal factor:

$$\Sigma_\mu = \nabla_\mu \ln A(\phi) = (\beta_A/M_{\text{Pl}}) \nabla_\mu\phi$$

It measures the local rate of change of the conformal clock-rate rescaling. It is not a separate interaction: matter follows geodesics of the matter metric, and in Einstein-frame variables the shear appears as the anomalous acceleration that fifth-force searches and PPN tests conventionally parameterize.

### A.4 Temporal Topology

The spatial profile and coherence structure of the scalar field $\phi$ (equivalently, of $\ln A$). It describes how the scalar field is organised in space and time, including its correlation lengths, domain boundaries, and topological defects. The coherence length $\lambda_T$ is the characteristic scale over which the field remains correlated; at smaller separations the field behaves coherently, while at larger scales random-phase superposition suppresses the net observable effect.

### A.5 Screening

Screening is the continuous suppression of the locally observable Temporal Shear/source-charge sector, expressed through the conformal factor $\ln A(\phi)$, its gradient $\Sigma_\mu$, and its covariance $C_A$. It is not a binary on/off switch but a smooth, environment-dependent suppression governed by the operator

$$\Sigma_\mu^{\text{obs}} = \mathcal S_\Sigma(\mathcal E) \, \Sigma_\mu,$$

where $\mathcal E = \{\rho, \Phi/c^2, \text{source structure}, \text{ambient environment}, \text{boundary conditions}, z, \text{measurement channel}\}$. The screening factor $\mathcal S_\Sigma$ is not a density threshold; it is a smooth suppression of the locally active shear.

Quantities such as $\rho_T$, $R_T(M)$, $S_\oplus(r)$, compactness $\Phi/c^2$, local stellar density, geometric coherence length, and channel-specific response coefficients are domain-specific projections of $\mathcal E$, not independent screening mechanisms and not interchangeable universal thresholds. Each is an observational transfer model that parameterizes the same underlying operator in a regime-appropriate form.

The saturation scale $\rho_T \approx 20$ g/cm$^3$ is the characteristic scale at which Temporal Topology effects saturate. It is **not** a local on/off condition of the form $\rho > \rho_T \Rightarrow$ GR and $\rho < \rho_T \Rightarrow$ active TEP. Rather, it is the scale at which the non-linear Temporal Topology response saturates. The subatomic core density ($\rho_{\text{core}} \sim 10^4$ g/cm$^3$, Paper 24), the macroscopic many-body suppression scale ($\rho_c \approx 20$ g/cm$^3$), and the galactic transition density ($\rho_{\text{half}} \approx 0.5 \, M_\odot/\text{pc}^3$, Paper 26) are different effective projections of the same non-linear Temporal Topology response. The first-principles transfer relation between them remains an open derivation.

Recovery of GR in local tests is controlled by suppression of the observable shear/source-charge sector. Source structure, environmental state, and boundary conditions suppress the locally active shear sector in screened regimes.

#### A.5.1 Screening Domain Table

The following canonical definitions are provided so that every paper in the TEP series references a single, unambiguous screening ontology. The underlying physics is one continuous suppression operator $\mathcal S_\Sigma(\mathcal E)$; the domain-specific forms listed below are observational transfer models, not separate fundamental laws.

| Domain | Screening variable | Observable response | Paper |
|---|---|---|---|
| GNSS / Earth | terrestrial topology, $L_c$, $\rho_T$ | clock covariance | 1, 2, 3, 4 |
| Cepheids / H$_0$ | $\sigma^2/c^2$ drives, $S(\rho)$ suppresses | period contraction, $H_0$ bias | 11 |
| Pulsars / clusters | cluster potential + suppressed density scaling | spin-down residual | 10 |
| Wide binaries | $R_s$, $\epsilon_{\text{env}}$, $\rho_T$ | velocity boost saturation | 13 |
| LLR | compactness + interior shielding | $\eta$ residual channel | 17 |
| JWST / high-$z$ | halo potential + density suppression | age/mass inflation | 12 |
| Cosmology / C0 | $S(\rho) \, \epsilon_T$, void vs mass-weighted average | redshift-distance transport | 26 |
| LHC / quantum | channel-specific response $\kappa_{\text{LHC}}$ | topological charge form factor | 20 |

Each entry uses a domain-appropriate parameterization of the same underlying operator $\mathcal S_\Sigma(\mathcal E)$. The LHC entry is explicitly marked as probing a channel-specific response under high-energy momentum transfer, not the macroscopic bulk-density screening that governs astrophysical and cosmological observables.

#### A.5.2 Screening Scales and Symbols

The following table lists the canonical screening-related symbols used across the TEP corpus, their meanings, and the papers in which they appear. This is provided to prevent symbol collisions and scale misattribution (e.g., $R_s$ for wide binaries vs $R_{\rm sol}$ for pulsar companion saturation).

| Symbol | Meaning | Defined in | Typical value / formula | Used in |
|---|---|---|---|---|
| $\rho_T$ | Temporal Topology saturation scale | Paper 6 | $\approx 20$ g/cm$^3$ | Papers 6, 10, 11, 13, 17 |
| $\rho_{\rm half}$ | Galactic stellar half-suppression density | Paper 11 | $\approx 0.5\,M_\odot/\text{pc}^3 \approx 3\times10^{-23}$ g/cm$^3$ | Paper 11 ($S(\rho_*)$ for Cepheids) |
| $\rho_c$ | Macroscopic many-body suppression scale | retired lab-sector estimate | $\approx 20$ g/cm$^3$ (equivalent to $\rho_T$) | none (superseded by $\rho_T$) |
| $R_T$ | Geometric saturation radius | Paper 6 | $\left(3M / 4\pi\rho_T\right)^{1/3}$ | Papers 6, 13 |
| $R_s$ | Effective screening radius (wide binaries) | Paper 13 | $\left(3M / 4\pi\epsilon_{\rm env}\rho_T\right)^{1/3}$ | Paper 13 |
| $R_{\rm sol}$ | Companion saturation radius (pulsars) | Paper 10 | $\sim$ companion orbital scale | Paper 10 |
| $\lambda_T$ | GNSS correlation length / relaxation scale | Papers 1–3 | $\sim 4{,}200$ km | Papers 1–3, 4, 14 |
| $\rho_{\rm amb}$ | Ambient halo density | Context-dependent | $\sim 10^{-18}$ g/cm$^3$ (GC halo) | Papers 10, 11, 13 |
| $\rho_*$ | Local stellar mass density | Context-dependent | $\sim 0.5\,M_\odot/\text{pc}^3$ (galactic disk) | Paper 11 |
| $\Lambda_X$ | Kinetic-completion gradient/shear scale, $P_{,X} = 1 + 2\vert X\vert/\Lambda_X^4$ | Paper 0 | $\sqrt{M_{\rm Pl}H_0} \approx 1.9$ meV | Papers 0, 4, 6, 7, 12, 13, 17, 26, 28 |
| $\lambda_{\rm quartic}$ | Quartic amplitude-sector coupling, $V = \lambda_{\rm quartic}M_{\rm Pl}^4 u^4/4$ | Paper 0 | $7.526\times10^{-66}$ (Cassini branch; reference $7.526\times10^{-71}$) | Papers 0, 6, 7, 10, 22, 28 |
| $\Lambda_V$ | Equivalent potential normalization $\lambda_{\rm quartic}^{1/4}M_{\rm Pl}$ (derived, not independent) | Paper 6 | $\approx 127.5$ GeV | Papers 0, 6 |
| $\rho_{\rm sat}$ | Self-quenching density ($u_{\min}\sim 1$ crossing; temporal-well regime) | Paper 6 | $\sim 6\times10^{25}$ g/cm$^3$ (order-of-magnitude) | Papers 6, 28 |
| $u_{\min}$ | Density-set scalar equilibrium, $\lambda_{\rm quartic}M_{\rm Pl}^4 u_{\min}^3 = \rho\,e^{-u_{\min}}$ | Paper 0 | $\left(\rho/\lambda_{\rm quartic}M_{\rm Pl}^4\right)^{1/3}$ | Papers 0, 6, 7, 28 |
| $u_s$ | Master-potential plateau knee | Paper 0 | $\sim 10$ (nominal; closure scans 10–50) | Papers 0, 6, 28 |

Key distinction: $\rho_{\rm half}$ is a **local stellar density** parameter for the continuous suppression factor $S(\rho_*)$ in galactic disks (Paper 11), not an ambient halo density. Globular clusters are active not because $\rho_{\rm amb} \ll \rho_{\rm half}$ but because their internal potential structure permits gradient coherence on scales larger than the cluster size. Conversely, $\rho_T \approx 20$ g/cm$^3$ is the asymptotic saturation scale of the conformal-factor sector; systems with $\rho \gg \rho_T$ (e.g., Solar System interiors) are in the saturated/GR-recovered regime.

Scale-sector distinction: the scalar sector carries three anchored scales that must not be conflated — $\Lambda_X$ (the kinetic gradient/shear scale, anchored to the measured drift $H_0$), $\lambda_{\rm quartic}\leftrightarrow\Lambda_V$ (the amplitude-sector potential normalization, anchored to the Cassini screening bound), and $\rho_T$ (the geometric calibration of $R_T(M)$, anchored to the GNSS correlation length). The energy-density reading $\rho_T^{1/4}\approx 96$ keV is not a field scale; the identification $\rho_T\equiv\Lambda_V^4$ is excluded by compact-object equilibria, Compton-resolution behaviour and Solar-System screening (Paper 6, Appendix C).

Symbol-overloading note: $\lambda$ carries several unrelated local meanings across the corpus — the quartic coupling $\lambda_{\rm quartic}$ above, the GNSS covariance length $\lambda_T$ (Papers 1–3, 6, 14), the $B(\phi)$ descent exponent in Paper 29 benchmark forms, and the affine parameter on geodesics (Paper 28). Manuscript-local definitions govern; where a single paper uses more than one meaning it disambiguates explicitly (Paper 6 writes the covariance length $\lambda_T$ alongside the quartic coupling).

### A.6 Observable Response Coefficient

An observable response coefficient, denoted $\kappa$ with an appropriate subscript (e.g. $\kappa_{\text{Cep}}$, $\kappa_{\text{MSP}}$), is a domain-level empirical transfer coefficient that maps a TEP field structure to a specific observable shift in a given experiment. It is not a bare microscopic coupling; it absorbs the combined effect of local screening, instrument bandwidth, calibration conventions, and systematic residuals. Different experiments constrain different response coefficients, and the same underlying $\beta_A$ can manifest through different $\kappa$ values in different domains.

### A.7 Synchronization Holonomy

The residual non-closure of synchronization transport around a closed loop, after subtraction of the full GR, kinematic, clock-scale, and reference-frame synchronization model:

$$H_{\rm resid}(C) \equiv \oint_C (\tilde{\sigma} - \sigma_{\rm GR}) = \iint_\Sigma (\tilde{F} - F_{\rm GR})$$

where $\sigma_{\rm GR}$ includes the standard Sagnac, gravitomagnetic/Lense--Thirring, Shapiro, gravitational-redshift, station-motion, clock-scale, and reference-frame contributions. This residual vanishes in SR/GR with $B=0$ and stationary $A$; it is generically nonzero with disformal corrections ($B \neq 0$) and/or temporal variation of $A$ combined with motion through $\nabla\phi$.

This quantity is invariant under admissible synchronization re-gaugings because the matter-frame connection $\tilde\sigma$ and the corresponding GR reference connection $\sigma_{\rm GR}$ shift by the same exact one-form, leaving the residual connection $\Delta\sigma=\tilde\sigma-\sigma_{\rm GR}$ unchanged. In the conformal-only limit ($B = 0$), the scalar contribution to closed-loop transport is generated by the exact one-form $d \ln A$, so $\oint_C d \ln A = 0$ and $H_{\rm resid}$ vanishes at leading order. A non-zero $H_{\rm resid}$ requires the disformal sector, non-metricity, or other explicitly non-exact transport structure.

Configurations in which $C$ does not bound a smooth surface — non-simply-connected transport, or a connection that is locally closed but not globally exact — are not described by this Stokes expression and are treated as a separate topological case.
