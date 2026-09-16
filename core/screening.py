#!/usr/bin/env python3
"""
TEP Screening Module
====================

Version: TEP v0.10 (Jakarta)

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

RHO_T = tep_const.RHO_T


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


def screening_factor(rho_local_g_cm3, rho_t=RHO_T):
    """
    Continuous Temporal Topology suppression factor for the scalar field source.

    When rho_local << rho_t: suppression -> 1 (full TEP effect)
    When rho_local -> rho_t: suppression -> 0.5 (transition)
    When rho_local >> rho_t: suppression -> 0 (saturated, A -> 1)

    WARNING: rho_t = 20.0 g/cm^3 is calibrated for lab/stellar-body densities.
    For galactic-scale densities (~1e-17 g/cm^3), rho/rho_t ~ 5e-19 and this
    function returns S ~ 1.0 for every galaxy, making it numerically useless as
    an environmental discriminant. For galaxy-scale Cepheid-host screening, use
    the galactic-onset density rho_half ~ 0.5 M_sun/pc^3 (see TEPCosmology).
    """
    return universal_screening_function(rho_local_g_cm3, rho_t, n=2.0, invert=False)


def coupling_screening_factor(rho_local_g_cm3, rho_transition=1.0, n=4.0):
    """
    Dimensionless density-screening factor for the coupling.

    f(rho) = 1 / [1 + (rho_transition / rho)^n]

    This is the inverted power-law form, used for chameleon-like
    coupling screening (suppressed at low density).
    """
    return universal_screening_function(rho_local_g_cm3, rho_transition, n=n, invert=True)


def beta_screened(rho_local_g_cm3, beta_A=tep_const.BETA_A,
                  rho_transition=1.0, n=4.0):
    """
    Screened conformal coupling beta_eff(rho) = beta_A * f(rho).

    Parameters
    ----------
    rho_local_g_cm3 : float or ndarray
        Local matter density in g/cm^3.
    beta_A : float
        Bare conformal coupling.
    rho_transition : float
        Density at which coupling is half-screened.
    n : float
        Steepness of the transition.

    Returns
    -------
    float or ndarray
        Effective conformal coupling after density screening.
    """
    f = coupling_screening_factor(rho_local_g_cm3, rho_transition, n)
    return beta_A * f
