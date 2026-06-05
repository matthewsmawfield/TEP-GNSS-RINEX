#!/usr/bin/env python3
"""
TEP Core Constants
================

Canonical physical and phenomenological parameters for the Temporal Equivalence
Principle (TEP) framework.  All TEP papers should import from this module to ensure
consistency across the corpus.  Human-readable registry: parameter_registry.yaml
in this directory.  Do not duplicate these values in project scripts.

Version: TEP v0.9 (Jakarta)
"""

import numpy as np

VERSION = "0.9"
VERSION_CODENAME = "Jakarta"
VERSION_STRING = f"TEP v{VERSION} ({VERSION_CODENAME})"

# =============================================================================
# PHYSICAL CONSTANTS (CODATA 2018)
# =============================================================================
G_NEWTON = 6.67430e-11          # m^3 kg^-1 s^-2
C_LIGHT = 299792458.0            # m s^-1
M_PLANCK = 2.176434e-8           # kg (Planck mass m_P = sqrt(hbar*c/G))
M_SUN = 1.98847e30               # kg
MPC_TO_M = 3.08567758e22         # m

# =============================================================================
# TEP UNIVERSAL PARAMETERS
# =============================================================================

# Conformal coupling strength.
# phi is dimensionless (measured in reduced-Planck-mass units:
# phi = phi_tilde / M_pl), so the conformal factor is:
#     A(phi) = exp(beta_A * phi)
# with no further M_pl normalization in the code.
BETA_A = -1.0                 # Dimensionless conformal coupling (locked lab-scale convention)

# Phenomenological screening coefficient in the TEP-SPIN tanh ansatz.
# This is NOT the fundamental conformal coupling (BETA_A). It is a
# calibrated parameter of the density-dependent screening model.
BETA_SPIN = 0.01                 # Dimensionless; Paper 24

# Solar-system PPN bound on conformal coupling from Cassini time-delay test.
BETA_CASSINI_MAX = 0.0034        # Bertotti et al. 2003

# Phenomenological saturation proximity scale for Temporal Topology screening.
# When local proximity approaches rho_c (observationally proxied by density),
# the scalar field saturates and A(phi) -> 1, suppressing TEP effects.
RHO_C = 20.0                     # g cm^-3

# Coherence length for lab-scale scalar field
LAB_COHERENCE_LENGTH_M = 50000.0  # 50 km crustal column

# Reference mass scale for geometric coupling beta_geom
M_REF = 1.0e18                   # kg (threshold mass where phi_mass ~ beta_geom)

# Temporal Topology coherence length (long-duration GNSS analysis)
SCREENING_LENGTH_KM = 4200.0      # km; canonical multi-center baseline (Paper 6)

# MGEX held-out verification (Paper 14; TEP-GNSS-MGEX step_2_1_correlation_length.json)
LAMBDA_T_MGEX_KM = 1396.19
LAMBDA_T_MGEX_ERR_KM = 90.19
LAMBDA_T_MGEX_R2 = 0.486

# Multi-center GNSS exponential fits (Paper 1; TEP-GNSS step_2_0_correlation_analysis_summary.json)
GNSS_LAMBDA_T_LONGSPAN_CODE_KM = 4201
GNSS_LAMBDA_T_LONGSPAN_CODE_ERR_KM = 1967
GNSS_LAMBDA_T_EXPONENTIAL_BY_CENTER = {
    "CODE": {"lambda_km": 4549, "ci_low_km": 1198, "ci_high_km": 5918},
    "IGS": {"lambda_km": 3764, "ci_low_km": 3197, "ci_high_km": 4871},
    "ESA": {"lambda_km": 3330, "ci_low_km": 2532, "ci_high_km": 3984},
}

# Lab-scale coupling constants (TEP-NIST Paper 21)
# alpha_log sign is fixed by the TEP field equation in the (+,-,-,-) metric
# signature: nabla_mu[K(phi) nabla^mu phi] = -alpha(phi) T with alpha = beta_A/M_Pl < 0.
# For non-relativistic dust T = +rho, the static limit gives nabla^2 phi ~ +|alpha| rho,
# so phi decreases with increasing density: dphi/drho < 0.  Since the
# phenomenological ansatz is phi_rho = alpha_log * ln(rho/rho_c), this requires
# alpha_log < 0.  The magnitude |7.66e-3| was determined from the requirement
# that the TEP model reproduce the correct order of magnitude for laboratory
# metrology shifts.
ALPHA_LOG = -7.66e-3             # Density-sector coupling (negative by field-equation sign)
BETA_GEOM = 1.50e-4              # Mass-sector geometric coupling

# =============================================================================
# SCREENING MODEL
# =============================================================================

def universal_screening_function(rho, rho_threshold, n=2.0, invert=False):
    """
    Universal TEP density screening function.

    Replaces previously incompatible phenomenological forms across the TEP corpus.
    Parameters match the NIST standard (Paper 21) exactly.

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
