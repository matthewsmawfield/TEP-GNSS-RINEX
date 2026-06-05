#!/usr/bin/env python3
"""
TEP Screening Module
====================

Version: TEP v0.9 (Jakarta)

Density-dependent screening functions for the Temporal Equivalence Principle.

All TEP papers should import from this module to ensure consistent screening
models across the corpus.
"""

import numpy as np
from . import constants as tep_const

RHO_C = tep_const.RHO_C


def universal_screening_function(rho, rho_threshold, n=2.0, invert=False):
    """
    Universal TEP density screening function.

    Parameters
    ----------
    rho : float or ndarray
        Local matter density.
    rho_threshold : float
        Transition density threshold (same units as rho).
    n : float
        Steepness of the power-law transition. Default is 2.0.
    invert : bool
        If False (default): factor = 1 / [1 + (rho/rho_threshold)^n].
        Used for source and cosmology screening (suppressed at high density).

        If True: factor = 1 / [1 + (rho_threshold/rho)^n].
        Used for chameleon coupling screening (suppressed at low density).
    """
    rho = np.asarray(rho, dtype=float)
    if invert:
        ratio = rho_threshold / rho
    else:
        ratio = rho / rho_threshold
    return 1.0 / (1.0 + ratio ** n)


def screening_factor(rho_local_g_cm3, rho_c=RHO_C):
    """
    Continuous Temporal Topology suppression factor for the scalar field source.

    When rho_local << rho_c: suppression -> 1 (full TEP effect)
    When rho_local -> rho_c: suppression -> 0.5 (transition)
    When rho_local >> rho_c: suppression -> 0 (saturated, A -> 1)
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
