#!/usr/bin/env python3
"""
TEP Scalar Field
================

Version: TEP v0.15 (Jakarta)

Canonical scalar-sector implementation (Paper 0, master-action closure;
verified by scripts/steps/step_27, step_30, step_50 and tep_model.py).

All field values are dimensionless u = phi/M_Pl (reduced Planck
convention).  The scalar sector is:

  * Potential: the master family
        V(u) = V_matter(u) e^{-(u/u_s)^4} + V_0 e^{-(u_s/u)^4}
    which on the matter-hosting weak-field branch (u << u_s ~ 10) is the
    quartic V = lambda u^4 / 4 with matter-sourced equilibrium
        lambda u^3 e^u = rho / M_Pl^4
    (quartic approximation u_min ~ (rho/(lambda M_Pl^4))^{1/3}).
  * Kinetic: P(X, phi) = X - V + X|X|/Lambda^4, Lambda^4 = M_Pl^2 H0^2,
    giving the screening operators
        S_Sigma(g) = [1 + (g/g_t)^2]^-1 ,  g_t = c H0 / (2 |beta_A|)
        S_eff(s)   = [1 + (R_s/s)^4]^-1 ,  R_s = sqrt(G M / g_t)
    and the flux-conserving profile y(1 + y^2 (g/g_t)^2) = 1.
  * Amplitude sector: S_A = min[1, (rho_bar/rho_T)^{1/3}].
  * Disformal: B(u) = B0 u^2/(1+u^2) exp(-u^4/2), admissible branch
    B0 >= 0; canonical benchmark B0 = +1 (GW170817 bounds B0 <~ 78 on
    the galactic profile).

-----------------------------------------------------------------------------
LEGACY: the functions solve_scalar_field_* and scalar_field_logarithmic
below implement the retired Naivasha lab-scale logarithmic
ansatz
    phi = alpha_log * ln(rho/rho_T) + beta_geom * ln(M/M_ref).
Its parameters alpha_log/beta_geom are class 'exploratory' in
parameter_registry.yaml.  This ansatz is NOT the canonical field solution:
it is a saturation-law fit to laboratory metrology shifts whose density
scaling (phi ~ ln rho) differs from the canonical equilibrium
(u_min ~ rho^{1/3}) and whose sign convention measures the screened
excursion rather than the field value.  It is retained for backward
compatibility of historical pipelines only; new work must use the
canonical functions above.  Known consumer: TEP-C0 step_04_08 (flagged
for re-derivation, see issues.md `core-scalar-field-legacy-ansatz`).
-----------------------------------------------------------------------------
"""

import numpy as np
from scipy.optimize import brentq
from . import constants as tep_const
from .screening import screening_factor

RHO_T = tep_const.RHO_T


# ============================================================================
# CANONICAL SECTOR — derived operators and profile solvers
# ============================================================================

def screening_y(g, g_t=None):
    """Flux-conserving screening factor: real root of y(1+y^2(g/g_t)^2)=1.

    Cardano solution of the depressed cubic y^3 + y/q^2 - 1/q^2 = 0 with
    q = g/g_t.  At g >> g_t, y ~ q^{-2/3} (strong suppression); at
    g << g_t, y -> 1 (unscreened).  This is the exact profile-level
    operator; S_sigma below is its algebraic weak-gradient projection.
    """
    if g_t is None:
        g_t = tep_const.G_T_TRANSITION
    g = np.asarray(g, dtype=float)
    qq = (g / g_t) ** 2
    y = np.empty_like(g)
    small = qq < 1e-8
    y[small] = 1.0 - qq[small]
    q = qq[~small]
    disc = (0.5 / q) ** 2 + (1.0 / (3.0 * q)) ** 3
    s = np.sqrt(np.maximum(disc, 0.0))
    y[~small] = np.cbrt(0.5 / q + s) + np.cbrt(0.5 / q - s)
    return y


def S_sigma(g, g_t=None):
    """Algebraic screening operator S_Sigma(g) = [1+(g/g_t)^2]^-1.

    Derived from P_,X of the kinetic completion (step_27); controls the
    observable Temporal Shear (gradient projection).
    """
    if g_t is None:
        g_t = tep_const.G_T_TRANSITION
    g = np.asarray(g, dtype=float)
    return 1.0 / (1.0 + (g / g_t) ** 2)


def S_eff_pairwise(separation_m, mass_kg, g_t=None):
    """Pairwise screening S_eff(s) = [1+(R_s/s)^4]^-1, R_s = sqrt(GM/g_t)."""
    if g_t is None:
        g_t = tep_const.G_T_TRANSITION
    s = np.asarray(separation_m, dtype=float)
    R_s = np.sqrt(tep_const.G_NEWTON * mass_kg / g_t)
    return 1.0 / (1.0 + (R_s / s) ** 4)


def R_s_transition(mass_kg, g_t=None):
    """Shear transition radius R_s = sqrt(GM/g_t)."""
    if g_t is None:
        g_t = tep_const.G_T_TRANSITION
    return float(np.sqrt(tep_const.G_NEWTON * mass_kg / g_t))


def S_A_density(rho_bar_g_cm3, rho_T=None):
    """Clock-amplitude sector S_A = min[1, (rho_bar/rho_T)^{1/3}].

    Quartic-branch result: the matter-sourced equilibrium field and the
    response amplitude scale as rho^{1/3} below the saturation density.
    """
    if rho_T is None:
        rho_T = RHO_T
    rho = np.asarray(rho_bar_g_cm3, dtype=float)
    return np.minimum(1.0, (rho / rho_T) ** (1.0 / 3.0))


def equilibrium_u(rho_g_cm3, lam=None):
    """Equilibrium dimensionless field u_min on the quartic branch.

    Solves  lam u^3 e^u = rho/M_Pl^4  (exact conformal-source
    equilibrium; for u << 1 the e^u factor is ~1 and this reduces to
    u_min = (rho/(lam M_Pl^4))^{1/3}).
    """
    if lam is None:
        lam = tep_const.LAMBDA_QUARTIC_CASSINI
    rho = np.asarray(rho_g_cm3, dtype=float)
    rho_gev4 = rho * tep_const.G_CM3_TO_GEV4
    rhs = rho_gev4 / (lam * tep_const.M_PL_REDUCED_GEV**4)
    out = np.zeros_like(rho)
    pos = rhs > 0
    for i in np.ndindex(*rhs.shape):
        if not pos[i]:
            continue
        log_rhs = np.log(rhs[i])
        # solve 3 log u + u = log_rhs for u (monotonic for u > 0)
        x = brentq(lambda t: 3.0 * t + np.exp(t) - log_rhs, -1000.0, 100.0)
        out[i] = np.exp(x)
    return out


def m_eff2_quartic(u, rho_g_cm3, lam=None):
    """Effective scalar mass squared (GeV^2) on the quartic branch.

    m_eff^2 = V''(u) + matter term  =  3 lam M_Pl^2 u^2 + rho e^{-u}/M_Pl^2
    with rho converted to GeV^4.
    """
    if lam is None:
        lam = tep_const.LAMBDA_QUARTIC_CASSINI
    u = np.asarray(u, dtype=float)
    rho_gev4 = np.asarray(rho_g_cm3, dtype=float) * tep_const.G_CM3_TO_GEV4
    Mpl = tep_const.M_PL_REDUCED_GEV
    return 3.0 * lam * (Mpl * u) ** 2 + rho_gev4 / Mpl**2 * np.exp(-u)


def compton_wavelength_m(u, rho_g_cm3, lam=None):
    """Compton wavelength lambda_c = hbar c / m_eff in metres."""
    m2 = m_eff2_quartic(u, rho_g_cm3, lam)
    m2 = np.asarray(m2, dtype=float)
    out = np.full_like(m2, np.nan)
    pos = m2 > 0
    out[pos] = tep_const.HBAR_C_GEV_M / np.sqrt(m2[pos])
    return out


def R_T_geometric(mass_kg, rho_T=None):
    """Geometric saturation radius R_T = (3M/4 pi rho_T)^{1/3}."""
    if rho_T is None:
        rho_T = RHO_T
    return float((3.0 * mass_kg / (4.0 * np.pi * 1000.0 * rho_T)) ** (1.0 / 3.0))


def master_potential(u, lam=None, u_s=None, V0=None):
    """Master potential V(u) in M_Pl^4 units (dimensionless u = phi/M_Pl).

        V(u) = (lam/4) u^4 exp(-(u/u_s)^4) + V0 exp(-(u_s/u)^4)

    The quartic factor is the matter-hosting branch; the essential-
    singularity floor V0 e^{-(u_s/u)^4} activates only at u ~ u_s
    (temporal-horizon pileup) and underflows identically to zero at all
    screening densities.
    """
    if lam is None:
        lam = tep_const.LAMBDA_QUARTIC_CASSINI
    if u_s is None:
        u_s = tep_const.U_S_TRANSITION
    if V0 is None:
        V0 = tep_const.V0_PLATEAU
    u = np.asarray(u, dtype=float)
    quartic = (lam / 4.0) * u**4 * np.exp(-(u / u_s) ** 4)
    floor = np.zeros_like(u)
    mask = np.abs(u) > 1e-10
    floor[mask] = V0 * np.exp(-(u_s / u[mask]) ** 4)
    return quartic + floor


def B_disformal(u, B0=None):
    """Disformal envelope family B(u) = B0 u^2/(1+u^2) exp(-u^4/2).

    Admissible branch B0 >= 0 (B < 0 produces elliptic counterexamples in
    the perfect-fluid principal symbol; Paper 0 SS4).  Canonical benchmark
    B0 = +1; GW170817 path bound B0 <~ 78 on the galactic profile
    (step_50).
    """
    if B0 is None:
        B0 = tep_const.B0_DISFORMAL_CANONICAL
    u = np.asarray(u, dtype=float)
    return B0 * u**2 / (1.0 + u**2) * np.exp(-0.5 * u**4)


def phi_profile_spherical(r_m, mass_kg, radius_m, g_t=None):
    """Screened scalar profile u(r) for a uniform sphere, integrated inward
    from a flat infinity (u -> 0 as r -> infinity).

    du/dr = -y(g) g/c^2 with y the flux-conserving screening factor.
    Returns (u, u_unscreened, y) interpolated onto the input radii —
    the same machinery as step_50 Gate A.
    """
    if g_t is None:
        g_t = tep_const.G_T_TRANSITION
    r = np.asarray(r_m, dtype=float)
    order = np.argsort(r)
    rr = r[order]

    def g_newton(x):
        x = np.asarray(x, dtype=float)
        g = np.empty_like(x)
        inn = x < radius_m
        g[inn] = tep_const.G_NEWTON * mass_kg * x[inn] / radius_m**3
        g[~inn] = tep_const.G_NEWTON * mass_kg / x[~inn] ** 2
        return g

    r_out = max(rr[-1] * 20.0, radius_m * 50.0)
    grid = np.geomspace(max(rr[0], radius_m * 1e-4), r_out, 20000)
    g = g_newton(grid)
    y = screening_y(g, g_t)
    integ = y * g / tep_const.C_LIGHT**2
    integ_u = g / tep_const.C_LIGHT**2
    dr = np.diff(grid)
    u_grid = np.zeros_like(grid)
    uu_grid = np.zeros_like(grid)
    u_grid[:-1] = np.cumsum((0.5 * (integ[1:] + integ[:-1]) * dr)[::-1])[::-1]
    uu_grid[:-1] = np.cumsum((0.5 * (integ_u[1:] + integ_u[:-1]) * dr)[::-1])[::-1]
    u_r = np.interp(rr, grid, u_grid)
    uu_r = np.interp(rr, grid, uu_grid)
    y_r = screening_y(g_newton(rr), g_t)
    u_out = np.empty_like(r)
    uu_out = np.empty_like(r)
    y_out = np.empty_like(r)
    u_out[order] = u_r
    uu_out[order] = uu_r
    y_out[order] = y_r
    return u_out, uu_out, y_out


# ============================================================================
# LEGACY SECTOR — retired Naivasha lab-scale logarithmic ansatz.
# Do not use for new predictions; see module docstring.
# ============================================================================

ALPHA_LOG = tep_const.ALPHA_LOG
BETA_GEOM = tep_const.BETA_GEOM


def solve_scalar_field_cylinder(total_mass_kg, radius_m, height_m,
                                density_g_cm3, alpha_log=ALPHA_LOG,
                                beta_geom=BETA_GEOM):
    """
    Compute dimensionless scalar field phi using the TEP lab-scale
    logarithmic model.

    Parameters
    ----------
    total_mass_kg : float
        Total mass of the cylinder
    radius_m : float
        Cylinder radius
    height_m : float
        Cylinder height
    density_g_cm3 : float
        Average density
    alpha_log : float
        Locked density-sector coupling (dimensionless)
    beta_geom : float
        Locked geometric coupling (dimensionless)

    Returns
    -------
    phi : float
        Dimensionless scalar field value
    screen : float
        Screening suppression factor (0-1)
    """
    screen = screening_factor(density_g_cm3, RHO_T)

    # Density-sector contribution with screening
    phi_rho = alpha_log * np.log(density_g_cm3 / RHO_T) * screen

    # Mass-sector contribution (geometric)
    phi_mass = beta_geom * np.log(total_mass_kg / tep_const.M_REF)

    phi = phi_rho + phi_mass
    return float(phi), float(screen)


def solve_scalar_field_uniform_density(rho_kg_m3, radius_m, height_m,
                                       alpha_log=ALPHA_LOG,
                                       beta_geom=BETA_GEOM):
    """Cross-check using density-based mass estimate."""
    density_g_cm3 = rho_kg_m3 / 1000.0
    total_mass_kg = np.pi * radius_m**2 * height_m * rho_kg_m3
    return solve_scalar_field_cylinder(
        total_mass_kg, radius_m, height_m, density_g_cm3,
        alpha_log, beta_geom
    )


def solve_scalar_field_layered(layers, radius_m=50000.0, alpha_log=ALPHA_LOG, beta_geom=BETA_GEOM):
    """
    Compute dimensionless scalar field phi for a layered mass model.

    Computes total mass and volume-averaged density, then calls the
    cylinder solver.

    Parameters
    ----------
    layers : list of dict
        Each dict must have 'thickness_km' and 'density_g_cm3'
    radius_m : float
        Radius of the cylinder
    alpha_log, beta_geom : float
        Coupling constants

    Returns
    -------
    phi : float
    screen : float
    total_mass_kg : float
    """
    total_mass_kg = 0.0
    total_vol_m3 = 0.0

    for layer in layers:
        thickness_m = layer["thickness_km"] * 1000.0
        density_kg_m3 = layer["density_g_cm3"] * 1000.0
        if density_kg_m3 <= 0:
            raise ValueError("density_g_cm3 must be strictly positive")
        vol = np.pi * radius_m**2 * thickness_m
        mass = vol * density_kg_m3
        total_mass_kg += mass
        total_vol_m3 += vol

    avg_density_kg_m3 = total_mass_kg / total_vol_m3
    avg_density_g_cm3 = avg_density_kg_m3 / 1000.0
    height_m = total_vol_m3 / (np.pi * radius_m**2)

    phi, screen = solve_scalar_field_cylinder(
        total_mass_kg, radius_m, height_m, avg_density_g_cm3, alpha_log, beta_geom
    )
    return phi, screen, total_mass_kg


def solve_scalar_field_layered_weighted(layers, radius_m=50000.0, z0_km=5.0,
                                        alpha_log=ALPHA_LOG, beta_geom=BETA_GEOM):
    """
    Compute dimensionless scalar field phi for a layered mass model using
    a near-field depth-weighting kernel K(z) = 1 / (1 + (z/z0)^2).

    This solver accounts for the extreme sensitivity of the scalar field
    gradient to near-field density variations (e.g., local facility geology).
    The unweighted average smooths over these variations.

    Parameters
    ----------
    layers : list of dict
        Each dict must have 'thickness_km' and 'density_g_cm3'
    radius_m : float
        Radius of the cylinder
    z0_km : float
        Coherence scale of the near-field kernel (usually 5 km)
    alpha_log, beta_geom : float
        Coupling constants

    Returns
    -------
    phi : float
    screen : float
    total_mass_kg : float
    kernel_info : dict
        Includes z0_km and total_kernel_weight
    """
    total_mass_kg = 0.0
    weighted_density_sum = 0.0
    total_weight = 0.0
    current_z_km = 0.0

    for layer in layers:
        thickness_km = layer["thickness_km"]
        density_g_cm3 = layer["density_g_cm3"]

        if density_g_cm3 <= 0:
            raise ValueError("density_g_cm3 must be strictly positive")

        z_top = current_z_km
        z_bot = current_z_km + thickness_km

        # Integral of the weighting kernel across this layer
        w_i = z0_km * (np.arctan(z_bot / z0_km) - np.arctan(z_top / z0_km))

        weighted_density_sum += w_i * density_g_cm3
        total_weight += w_i

        thickness_m = thickness_km * 1000.0
        density_kg_m3 = density_g_cm3 * 1000.0
        vol = np.pi * radius_m**2 * thickness_m
        total_mass_kg += vol * density_kg_m3

        current_z_km = z_bot

    eff_density_g_cm3 = weighted_density_sum / total_weight

    # Calculate scalar field with the effective near-field weighted density
    screen = screening_factor(eff_density_g_cm3, RHO_T)
    phi_rho = alpha_log * np.log(eff_density_g_cm3 / RHO_T) * screen
    phi_mass = beta_geom * np.log(total_mass_kg / tep_const.M_REF)
    phi = phi_rho + phi_mass

    kernel_info = {
        "z0_km": z0_km,
        "total_kernel_weight": total_weight
    }

    return float(phi), float(screen), float(total_mass_kg), kernel_info


def scalar_field_logarithmic(density_g_cm3, total_mass_kg,
                             alpha_log=ALPHA_LOG, beta_geom=BETA_GEOM):
    """
    Scalar field from the TEP lab-scale logarithmic ansatz.

    phi = alpha_log * ln(rho / rho_T) + beta_geom * ln(M / M_ref)

    Parameters
    ----------
    density_g_cm3 : float or ndarray
        Local matter density in g/cm^3.
    total_mass_kg : float
        Total mass in kg.
    alpha_log : float
        Density-sector coupling.
    beta_geom : float
        Mass-sector geometric coupling.

    Returns
    -------
    float or ndarray
        Dimensionless scalar field phi.
    """
    if np.any(np.asarray(density_g_cm3) <= 0):
        raise ValueError("density_g_cm3 must be strictly positive")
    if total_mass_kg <= 0:
        raise ValueError("total_mass_kg must be strictly positive")

    phi_rho = alpha_log * np.log(np.asarray(density_g_cm3) / RHO_T)
    phi_geom = beta_geom * np.log(total_mass_kg / tep_const.M_REF)
    return phi_rho + phi_geom


def scalar_field_difference(density_1, density_2, mass_1, mass_2,
                            alpha_log=ALPHA_LOG, beta_geom=BETA_GEOM):
    """
    Scalar field difference between two sites.

    Delta_phi = phi(site_2) - phi(site_1)
              = alpha_log * ln(rho_2 / rho_1) + beta_geom * ln(M_2 / M_1)
    """
    if density_1 <= 0 or density_2 <= 0:
        raise ValueError("Densities must be strictly positive")
    if mass_1 <= 0 or mass_2 <= 0:
        raise ValueError("Masses must be strictly positive")

    delta_phi_rho = alpha_log * np.log(density_2 / density_1)
    delta_phi_geom = beta_geom * np.log(mass_2 / mass_1)
    return delta_phi_rho + delta_phi_geom


def compute_temporal_shear_from_mass_gradient(phi, mass_kg, density_g_cm3,
                                              radius_m, height_m,
                                              beta_A=tep_const.BETA_A):
    """
    Deprecated v0.1 scalar diagnostic: Σ ≈ β_A ∇φ with ∇φ ≈ φ / L_char.

    Retained for backwards comparison in NIST step 06. The 3D FEM gradient
    is the primary shear estimator.
    """
    if radius_m <= 0:
        raise ValueError("radius_m must be strictly positive")
    l_char = float(radius_m)
    grad_phi = float(phi) / l_char
    sigma = float(beta_A) * grad_phi
    return sigma, grad_phi, l_char
