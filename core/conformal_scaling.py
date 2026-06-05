#!/usr/bin/env python3
"""
TEP Conformal Scaling
=====================

Version: TEP v0.9 (Jakarta)

Computes the conformal factor A(phi), Temporal Shear Sigma_mu,
and effective gravitational coupling from the scalar field profile phi.

All TEP papers should import from this module to ensure consistent
definitions across the corpus.
"""

import numpy as np
from . import constants as tep_const

BETA_A = tep_const.BETA_A
BETA_CASSINI_MAX = tep_const.BETA_CASSINI_MAX
M_PLANCK = tep_const.M_PLANCK

DEFAULT_RHO_TRANSITION = 1.0
DEFAULT_SCREENING_STEEPNESS = 4.0


def _as_numeric_array(value, name):
    arr = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite numeric values")
    return arr


def _return_scalar_if_scalar(original, value):
    if np.isscalar(original):
        return float(value)
    return value


def conformal_factor(phi, beta_A=BETA_A):
    """
    Compute conformal factor A(phi) = exp(beta_A * phi).

    In the TEP framework, phi is dimensionless (measured in Planck-mass
    units: phi = phi_tilde / M_pl). The conformal factor relates the
    effective metric to the background metric.
    """
    phi_arr = _as_numeric_array(phi, "phi")
    if not np.isfinite(beta_A):
        raise ValueError("beta_A must be finite")
    factor = np.exp(beta_A * phi_arr)
    return _return_scalar_if_scalar(phi, factor)


def conformal_factor_small(phi, beta_A=BETA_A):
    """
    Linearized conformal factor for small phi:
    A(phi) ~ 1 + beta_A*phi
    """
    phi_arr = _as_numeric_array(phi, "phi")
    if not np.isfinite(beta_A):
        raise ValueError("beta_A must be finite")
    factor = 1.0 + beta_A * phi_arr
    return _return_scalar_if_scalar(phi, factor)


def temporal_shear_from_scalar_field(phi, grad_phi, beta_A=BETA_A):
    """
    Compute Temporal Shear vector Sigma_mu = nabla_mu ln A(phi).

    Parameters
    ----------
    phi : float
        Scalar field value (dimensionless)
    grad_phi : array-like
        Spatial gradient of phi (components in m^-1)

    Returns
    -------
    sigma : ndarray
        Temporal Shear vector (same shape as grad_phi)
    magnitude : float
        |Sigma| in m^-1
    """
    if not np.isfinite(beta_A):
        raise ValueError("beta_A must be finite")
    _as_numeric_array(phi, "phi")
    grad = _as_numeric_array(grad_phi, "grad_phi")
    sigma = beta_A * grad
    magnitude = float(np.linalg.norm(sigma))
    return sigma, magnitude


def effective_g(newton_g, phi, beta_A=BETA_A):
    r"""
    Effective gravitational constant under TEP (Jordan frame).

    For the matter metric \tilde{g}_{\mu\nu} = A^2(\phi) g_{\mu\nu},
    the Jordan-frame effective gravitational constant is:

        G_eff = G_N / A^2(\phi)

    This convention means that larger phi (stronger scalar field)
    -> larger A(phi) -> smaller G_eff.
    """
    if not np.isfinite(newton_g) or newton_g <= 0:
        raise ValueError("newton_g must be finite and strictly positive")
    a = conformal_factor(phi, beta_A)
    return newton_g / (a ** 2)


def g_eff_variance(phi_1, phi_2, beta_A=BETA_A, beta_2=None):
    r"""
    Predicted relative variance in G_eff between two sites.

    For the matter metric \tilde{g}_{\mu\nu} = A^2(\phi) g_{\mu\nu},
    the Jordan-frame effective gravitational constant is G_eff = G_N / A^2.
    Therefore:

    Delta_G / G = [G_eff(site_2) - G_eff(site_1)] / G_eff(site_1)
                = [A(phi_1)/A(phi_2)]^2 - 1

    For small phi differences with equal beta_A:
    Delta_G / G ~ -2*beta_A*(phi_2 - phi_1)

    Parameters
    ----------
    phi_1, phi_2 : float
        Scalar field at sites 1 and 2.
    beta_A : float
        Conformal coupling at site 1 (default BETA_A).
    beta_2 : float, optional
        Conformal coupling at site 2. If None, uses beta_A.
    """
    b2 = beta_2 if beta_2 is not None else beta_A
    _as_numeric_array(phi_1, "phi_1")
    _as_numeric_array(phi_2, "phi_2")
    if not np.isfinite(beta_A) or not np.isfinite(b2):
        raise ValueError("beta_A and beta_2 must be finite")
    a1 = conformal_factor(phi_1, beta_A)
    a2 = conformal_factor(phi_2, b2)
    return (a1 / a2) ** 2 - 1.0


def coupling_screening_factor(rho_local_g_cm3, rho_transition=DEFAULT_RHO_TRANSITION,
                              n=DEFAULT_SCREENING_STEEPNESS):
    """
    Dimensionless density-screening factor f(rho) for the coupling.

    The implemented coupling screen can be written:

        beta_eff(rho) = beta_A * f(rho)
        f(rho) = 1 / [1 + (rho_transition / rho)^n]

    This is the inverted power-law form of a density-threshold response.
    It exposes the exact assumptions needed for PPN compatibility:
    the transition density and steepness.
    """
    rho = _as_numeric_array(rho_local_g_cm3, "rho_local_g_cm3")
    if np.any(rho <= 0):
        raise ValueError("rho_local_g_cm3 must be strictly positive")
    if not np.isfinite(rho_transition) or rho_transition <= 0:
        raise ValueError("rho_transition must be finite and strictly positive")
    if not np.isfinite(n) or n <= 0:
        raise ValueError("n must be finite and strictly positive")

    factor = tep_const.universal_screening_function(rho, rho_transition, n=n, invert=True)
    return _return_scalar_if_scalar(rho_local_g_cm3, factor)


def beta_screened(rho_local_g_cm3, beta_A=BETA_A,
                  rho_transition=DEFAULT_RHO_TRANSITION,
                  n=DEFAULT_SCREENING_STEEPNESS):
    """
    Density-dependent conformal coupling (chameleon-like screening).

    The lab-scale locked coupling |beta_A| = 1.0 violates the Cassini PPN
    bound |beta_PPN| < 0.0034 (solar-system scale). A proximity-dependent
    coupling suppresses beta_A in dilute environments while preserving
    the full lab-scale strength inside dense matter:

        beta_eff(rho) = beta_A / [1 + (rho_transition / rho)^n]

    Parameters
    ----------
    rho_local_g_cm3 : float or ndarray
        Local matter density in g/cm^3.
    beta_A : float
        Bare conformal coupling (default BETA_A = -1.0).
    rho_transition : float
        Transition density in g/cm^3.
    n : float
        Steepness of the transition (default 4).

    Returns
    -------
    beta_eff : float or ndarray
        Screened coupling at the given density.
    """
    if not np.isfinite(beta_A):
        raise ValueError("beta_A must be finite")
    beta_eff = beta_A * coupling_screening_factor(
        rho_local_g_cm3,
        rho_transition=rho_transition,
        n=n,
    )
    return _return_scalar_if_scalar(rho_local_g_cm3, beta_eff)


def minimum_steepness_for_retention(rho_g_cm3, retention_fraction,
                                  rho_transition=DEFAULT_RHO_TRANSITION):
    """
    Minimum n such that f(rho) >= retention_fraction for rho > rho_transition.
    """
    if not np.isfinite(rho_g_cm3) or rho_g_cm3 <= rho_transition:
        raise ValueError("rho_g_cm3 must be finite and greater than rho_transition")
    if not np.isfinite(retention_fraction) or not (0 < retention_fraction < 1):
        raise ValueError("retention_fraction must be in (0, 1)")
    if not np.isfinite(rho_transition) or rho_transition <= 0:
        raise ValueError("rho_transition must be finite and strictly positive")

    odds = retention_fraction / (1.0 - retention_fraction)
    return float(np.log(odds) / np.log(rho_g_cm3 / rho_transition))


def minimum_steepness_for_suppression(rho_g_cm3, max_fraction,
                                      rho_transition=DEFAULT_RHO_TRANSITION):
    """
    Minimum n such that f(rho) <= max_fraction for rho < rho_transition.
    """
    if not np.isfinite(rho_g_cm3) or rho_g_cm3 <= 0 or rho_g_cm3 >= rho_transition:
        raise ValueError("rho_g_cm3 must be finite and between 0 and rho_transition")
    if not np.isfinite(max_fraction) or not (0 < max_fraction < 1):
        raise ValueError("max_fraction must be in (0, 1)")
    if not np.isfinite(rho_transition) or rho_transition <= 0:
        raise ValueError("rho_transition must be finite and strictly positive")

    max_odds = max_fraction / (1.0 - max_fraction)
    return float(np.log(max_odds) / np.log(rho_g_cm3 / rho_transition))


def screening_diagnostics(beta_A=BETA_A, ppn_beta_bound=BETA_CASSINI_MAX,
                          lab_density_g_cm3=2.5, solar_density_g_cm3=1e-24,
                          rho_transition=DEFAULT_RHO_TRANSITION,
                          n=DEFAULT_SCREENING_STEEPNESS,
                          lab_retention_target=0.95):
    """
    Summarize whether the density-screening choice satisfies lab and PPN needs.

    Returns a dict with diagnostics that can be used for validation and
    reporting across all TEP papers.
    """
    if not np.isfinite(beta_A) or beta_A == 0:
        raise ValueError("beta_A must be finite and non-zero")
    if not np.isfinite(ppn_beta_bound) or ppn_beta_bound <= 0:
        raise ValueError("ppn_beta_bound must be finite and strictly positive")

    required_solar_fraction = ppn_beta_bound / abs(beta_A)
    lab_factor = coupling_screening_factor(lab_density_g_cm3, rho_transition, n)
    solar_factor = coupling_screening_factor(solar_density_g_cm3, rho_transition, n)
    lab_min_n = minimum_steepness_for_retention(
        lab_density_g_cm3,
        lab_retention_target,
        rho_transition,
    )
    solar_min_n = minimum_steepness_for_suppression(
        solar_density_g_cm3,
        required_solar_fraction,
        rho_transition,
    )
    required_n = max(lab_min_n, solar_min_n)

    return {
        "model": "f(rho)=1/[1+(rho_transition/rho)^n]",
        "interpretation": "log-density logistic threshold for beta_eff=beta_A*f(rho)",
        "beta_A": float(beta_A),
        "rho_transition_g_cm3": float(rho_transition),
        "screening_steepness_n": float(n),
        "ppn_beta_bound": float(ppn_beta_bound),
        "required_solar_screening_fraction": float(required_solar_fraction),
        "lab_density_g_cm3": float(lab_density_g_cm3),
        "solar_density_g_cm3": float(solar_density_g_cm3),
        "lab_retention_target": float(lab_retention_target),
        "lab_screening_factor": float(lab_factor),
        "solar_screening_factor": float(solar_factor),
        "beta_eff_lab": float(beta_A * lab_factor),
        "beta_eff_solar": float(beta_A * solar_factor),
        "minimum_n_for_lab_retention": float(lab_min_n),
        "minimum_n_for_solar_suppression": float(solar_min_n),
        "minimum_n_required": float(required_n),
        "chosen_n_exceeds_requirement": bool(n >= required_n),
        "ppn_satisfied_at_solar_density": bool(abs(beta_A * solar_factor) <= ppn_beta_bound),
        "lab_retention_satisfied": bool(lab_factor >= lab_retention_target),
    }
