## TEP Parameter Registry

This table records the flow of parameters across the TEP corpus. Each parameter is classified as **fundamental**, **calibrated**, **observable response**, **nuisance**, or **exploratory**.

| Symbol | Meaning | Value | Class | First calibrated | Used in |
|---|---|---|---|---|---|
| $\beta_A$ | Conformal coupling strength in $A(\phi) = \exp(\beta_A\phi/M_{\text{Pl}})$ | $-1.0$ | fundamental | Paper 0 | 0, 1, 2, 3, 6, 10, 11, 14, 17, 18, 19, 21, 22, 23, 24, 25, 26 |
| $\rho_c$ | Temporal Topology saturation scale | $20.0$ g/cm$^3$ | calibrated | Paper 6 | 6, 14, 17, 18, 21, 23, 24 |
| $\lambda_T$ | Clock-correlation length (Temporal Topology coherence) | $4200$ km | calibrated | Papers 1--2 | 1, 2, 3, 6, 14, 17, 21 |
| $\alpha_{\log}$ | Lab-scale density-sector coupling | $-7.66 \times 10^{-3}$ | calibrated | Paper 21 | 21 |
| $\beta_{\text{geom}}$ | Mass-sector geometric coupling | $1.50 \times 10^{-4}$ | calibrated | Paper 21 | 21 |
| $M_{\text{ref}}$ | Reference mass scale for geometric coupling | $10^{18}$ kg | fundamental | Paper 21 | 21 |
| $\kappa_{\text{Cep}}$ | Cepheid observable response coefficient | -- | observable response | Paper 11 | 11, 12 |
| $\kappa_{\text{MSP}}$ | Pulsar observable response coefficient | -- | observable response | Paper 10 | 10, 17 |
| $\epsilon_T$ | Cosmology shear amplitude | -- | exploratory | Paper 26 | 18, 26 |
| $\phi_0$ | Present-day scalar field amplitude (phenomenological) | -- | exploratory | Paper 22 | 18, 22, 26 |
| $n$ | Redshift exponent in phenomenological conformal factor | -- | exploratory | Paper 22 | 18, 22 |
| $B(\phi)$ | Disformal coupling function | 0 (leading order) | fundamental | Paper 0 | 0, 4, 9, 22 |

**Class definitions:**

- **fundamental** — Axiomatic or theoretically fixed; changes would alter the framework itself.
- **calibrated** — Determined from experimental data with well-defined uncertainty; shared across analyses.
- **observable response** — Domain-level empirical transfer coefficient; specific to a measurement domain.
- **nuisance** — Required for analysis but not a TEP parameter of direct physical interest.
- **exploratory** — Under investigation; value or existence not yet firmly established.
