#!/usr/bin/env python3
"""
TEP Core Constants
================

Canonical physical and phenomenological parameters for the Temporal Equivalence
Principle (TEP) framework.  All TEP papers should import from this module to ensure
consistency across the corpus.  Human-readable registry: parameter_registry.yaml
in this directory.  Do not duplicate these values in project scripts.

Version: TEP v0.10 (Jakarta)
"""

import numpy as np

VERSION = "0.10"
VERSION_CODENAME = "Jakarta"
VERSION_STRING = f"TEP v{VERSION} ({VERSION_CODENAME})"

# =============================================================================
# PHYSICAL CONSTANTS (CODATA 2018)
# =============================================================================
G_NEWTON = 6.67430e-11          # m^3 kg^-1 s^-2
C_LIGHT = 299792458.0            # m s^-1
M_PLANCK = 2.176434e-8           # kg (Planck mass m_P = sqrt(hbar*c/G))
M_SUN = 1.98847e30               # kg
M_EARTH = 5.972e24               # kg
R_EARTH = 6.371e6                # m
R_SUN = 6.96e8                   # m
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
# calibrated parameter of the environment-dependent screening model.
BETA_SPIN = 0.01                 # Dimensionless; Paper 24

# Illustrative coupling used in numerical lattice solvers where |BETA_A| = 1
# would cause overflow in FFT-based solvers. This is a numerical convenience
# parameter (not a physical coupling). The lattice solver rescales results back
# to physical values using the mean-field ratio.
ILLUSTRATIVE_BETA_A = 0.01       # Dimensionless; used in step_09 and legacy scripts

# Solar-system PPN bound on conformal coupling from Cassini time-delay test.
BETA_CASSINI_MAX = 0.0034        # Bertotti et al. 2003

# Phenomenological saturation scale for Temporal Topology screening.
# When local proximity approaches rho_T (observationally proxied by density),
# the scalar field saturates and A(phi) -> 1, suppressing TEP effects.
# This is a candidate saturation scale, NOT a binary density threshold.
RHO_T = 20.0                     # g cm^-3

# Backward-compatible alias (deprecated; use RHO_T in new code)
RHO_C = RHO_T

# Coherence length for lab-scale scalar field
LAB_COHERENCE_LENGTH_M = 50000.0  # 50 km crustal column

# Reference mass scale for geometric coupling beta_geom
M_REF = 1.0e18                   # kg (threshold mass where phi_mass ~ beta_geom)

# Temporal Topology coherence length — canonical value for all forward analysis
# (25-year multi-center GNSS baseline, Papers 1–2/6). Do not substitute short-run
# verification estimates (e.g. Paper 14 MGEX ~1396 km on ~1 yr span).
SCREENING_LENGTH_KM = 4200.0

# MGEX held-out verification only (Paper 14; TEP-GNSS-MGEX step_2_1). Not used in
# NIST/UCD forward models or FEM dimensional normalization.
LAMBDA_T_MGEX_KM = 1396.19
LAMBDA_T_MGEX_ERR_KM = 90.19
LAMBDA_T_MGEX_R2 = 0.486

# Unit conversion: kg/m^3 <-> g/cm^3
# 1 kg/m^3 = 1000 g / 1,000,000 cm^3 = 10^-3 g/cm^3
# Therefore: g/cm^3 = kg/m^3 / 1000.0 and kg/m^3 = g/cm^3 * 1000.0
KG_M3_TO_G_CM3 = 1e-3   # multiply kg/m^3 by this to get g/cm^3
G_CM3_TO_KG_M3 = 1e3    # multiply g/cm^3 by this to get kg/m^3

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
# phenomenological ansatz is phi_rho = alpha_log * ln(rho/rho_T), this requires
# alpha_log < 0.  The magnitude |7.66e-3| was determined from the requirement
# that the TEP model reproduce the correct order of magnitude for laboratory
# metrology shifts.
ALPHA_LOG = -7.66e-3             # Density-sector coupling (negative by field-equation sign)
BETA_GEOM = 1.50e-4              # Mass-sector geometric coupling

# =============================================================================
# GALAXY-SCALE OBSERVABLE PARAMETERS
# =============================================================================

# Canonical galaxy-scale observable response coefficient (Paper 11).
KAPPA_GAL = 9.6e5                # mag
KAPPA_GAL_UNCERTAINTY = 4.0e5    # mag

# Stellar evolution index (M/L ~ t^n from stellar isochrones)
ALPHA_NUCLEAR = 0.7

# Reference halo mass for potential calculations
LOG_MH_REF = 12.0

# Dimensionless virial potential Phi/c^2 for 10^12 Msun halo at z=0
PHI_REF_0 = 1.6e-7

# Reference redshift for chronological enhancement
Z_REF = 5.5
