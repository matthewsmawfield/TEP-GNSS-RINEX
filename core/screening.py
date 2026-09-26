#!/usr/bin/env python3
"""
TEP Screening Module
====================

Version: TEP v0.14 (Jakarta)

Environment-dependent Temporal Shear suppression for the Temporal Equivalence Principle.

Screening is the continuous suppression of locally observable Temporal Shear,
Sigma_mu^obs = S_Sigma(E) Sigma_mu, where E includes density, potential, source
structure, boundary conditions, and measurement channel. The functions here
provide domain-appropriate parameterizations of the underlying operator; they
are transfer models, not separate fundamental laws.

All TEP papers should import from this module to ensure consistent screening
models across the corpus.
"""

import numpy as np
from . import constants as tep_const

RHO_C = tep_const.RHO_C


def universal_screening_function(rho, rho_scale, n=2.0, invert=False):
    """
    Universal TEP density screening function.

    Parameters
    ----------
    rho : float or ndarray
        Local matter density.
    rho_scale : float
        Transition density scale (same units as rho).
    n : float
        Steepness of the power-law transition. Default is 2.0.
    invert : bool
        If False (default): factor = 1 / [1 + (rho/rho_scale)^n].
        Used for source and cosmology screening (suppressed at high density).

        If True: factor = 1 / [1 + (rho_scale/rho)^n].
        Used for chameleon coupling screening (suppressed at low density).
    """
    rho = np.asarray(rho, dtype=float)
    if np.any(rho <= 0):
        raise ValueError("rho must be strictly positive")
    if rho_scale <= 0:
        raise ValueError("rho_scale must be strictly positive")
    if n <= 0:
        raise ValueError("n must be strictly positive")
    if invert:
        ratio = rho_scale / rho
    else:
        ratio = rho / rho_scale
    return 1.0 / (1.0 + ratio ** n)


def screening_factor(rho_local_g_cm3, rho_c=RHO_C):
    """
    Continuous Temporal Topology suppression factor for the scalar field source.

    When rho_local << rho_c: suppression -> 1 (full TEP effect)
    When rho_local -> rho_c: suppression -> 0.5 (transition)
    When rho_local >> rho_c: suppression -> 0 (saturated, A -> 1)

    WARNING: rho_c = 20.0 g/cm^3 is calibrated for lab/stellar-body densities.
    For galactic-scale densities (~1e-17 g/cm^3), rho/rho_c ~ 5e-19 and this
    function returns S ~ 1.0 for every galaxy, making it numerically useless as
    an environmental discriminant. For galaxy-scale Cepheid-host screening, use
    the galactic-onset density rho_half ~ 0.5 M_sun/pc^3 (see TEPCosmology).
    """
    return universal_screening_function(rho_local_g_cm3, rho_c, n=2.0, invert=False)


def coupling_screening_factor(rho_local_g_cm3, rho_transition=1.0, n=4.0):
    """
    Dimensionless density-screening factor for the coupling.

    f(rho) = 1 / [1 + (rho_transition / rho)^n]

    This is the inverted power-law form, used for chameleon-like
    coupling screening (suppressed at low density).
    """
    return universal_screening_function(rho_local_g_cm3, rho_transition, n=n, invert=True)


def physical_shear(grad_phi_solved, beta_A=tep_const.BETA_A):
    """
    Physical temporal shear is exactly the gradient of the conformal factor
    Sigma_mu = nabla_mu ln A(phi) = beta_A * nabla_mu phi
    for the exact, fully solved (and thus already screened) scalar field.

    Parameters
    ----------
    grad_phi_solved : float or ndarray
        The gradient of the fully solved (nonlinear) scalar field.
    beta_A : float
        Bare conformal coupling.

    Returns
    -------
    float or ndarray
        Physical Temporal Shear.
    """
    return beta_A * np.asarray(grad_phi_solved, dtype=float)


def screening_ratio(q_screened, q_reference):
    """
    Screening ratio S_Sigma describes how much the nonlinear solution
    differs from an unscreened reference solution.

    It should not automatically multiply the gradient of an already-screened field.

    Parameters
    ----------
    q_screened : float or ndarray
        Observable quantity (e.g., charge, force, shear) from the nonlinear solved field.
    q_reference : float or ndarray
        The same observable from an unscreened reference solution.

    Returns
    -------
    float or ndarray
        The ratio S_Sigma = Q_screened / Q_reference.
    """
    q_s = np.asarray(q_screened, dtype=float)
    q_r = np.asarray(q_reference, dtype=float)
    # Avoid division by zero warnings
    return np.divide(q_s, q_r, out=np.zeros_like(q_s), where=q_r!=0)


def matter_acceleration(grad_newton, physical_shear_val, c=tep_const.C_LIGHT):
    """
    Acceleration directly from the matter metric in the weak-field limit.
    In the appropriate limit, its scalar contribution is obtained from
    the physical nabla ln A, not from an additional arbitrary suppression multiplier.

    Parameters
    ----------
    grad_newton : float or ndarray
        The Newtonian gravitational acceleration (grad Phi_N).
    physical_shear_val : float or ndarray
        The physical temporal shear (nabla ln A).
    c : float
        Speed of light.

    Returns
    -------
    float or ndarray
        Total acceleration.
    """
    return - np.asarray(grad_newton, dtype=float) - (c**2) * np.asarray(physical_shear_val, dtype=float)
