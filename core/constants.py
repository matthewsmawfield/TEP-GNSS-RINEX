#!/usr/bin/env python3
"""
TEP Core Constants
================

Canonical physical and phenomenological parameters for the Temporal Equivalence
Principle (TEP) framework.  All TEP papers should import from this module to ensure
consistency across the corpus.  Human-readable registry: parameter_registry.yaml
in this directory.  Do not duplicate these values in project scripts.

Version: TEP v0.15 (Jakarta)
"""

import numpy as np

VERSION = "0.15"
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
BETA_SPIN_PHEN = 0.01          # Phenomenological screening coefficient beta_spin (tanh ansatz), Paper 29 lattice demos only -- NOT the conformal coupling; opposite-sign to locked BETA_A = -1.0
ILLUSTRATIVE_BETA_A = BETA_SPIN_PHEN  # deprecated alias, kept for backcompat

# Solar-system PPN bound from the Cassini time-delay test (Bertotti et al.
# 2003, |gamma-1| < 2.3e-5). Under the canonical linear source-charge map
# gamma-1 = -4 beta_A^2 S_Sigma / (1 + 2 beta_A^2 S_Sigma), the bound applies
# to the screened solar source-charge fraction: |beta_A * S_Sigma| < 5.75e-6
# for beta_A = -1 (equivalently S_Sigma < 5.8e-6). It is NOT a bound on the
# bare coupling beta_A itself.
BETA_CASSINI_MAX = 5.75e-6       # bound on |beta_eff| = |beta_A * S_Sigma|

# Phenomenological saturation scale for Temporal Topology screening.
# rho_T is the geometric reference density of the R_T(M) law — the density
# argument of the amplitude response S_A = min[1, (rho/rho_T)^(1/3)]. At
# rho ~ rho_T the clock-amplitude response reaches its ceiling; the interior
# field stays at its local equilibrium u_min(rho) > 0 (A < 1, clocks slowed)
# while the shear gradient is pinned — saturation suppresses the response,
# it does not relax the field to vacuum. NOT a binary density threshold.
RHO_T = 20.0                     # g cm^-3

# Backward-compatible alias (deprecated; use RHO_T in new code)
RHO_C = RHO_T

# Coherence length for lab-scale scalar field
LAB_COHERENCE_LENGTH_M = 50000.0  # 50 km crustal column

# Reference mass scale for geometric coupling beta_geom
M_REF = 1.0e18                   # kg (threshold mass where phi_mass ~ beta_geom)

# Temporal Topology coherence length — canonical value for all forward analysis
# (25-year multi-center GNSS baseline, Papers 1–2/6). Do not substitute short-run
# verification estimates (e.g. Paper 14 MGEX ~1862 km on a ~15-month span).
SCREENING_LENGTH_KM = 4200.0

# MGEX held-out verification only (Paper 14; TEP-GNSS-MGEX step_2_1, corrected
# pooled 462-day estimate; the pre-correction ~1396 km figure for the original
# ~1 yr span is superseded). Short-baseline corroboration of the canonical
# anchor only — not used in NIST/UCD forward models or FEM dimensional
# normalization.
LAMBDA_T_MGEX_KM = 1861.63
LAMBDA_T_MGEX_ERR_KM = 111.75
LAMBDA_T_MGEX_R2 = 0.699

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

# Lab-scale coupling constants (retired lab-scale phenomenological ansatz; see
# the legacy notice in core/scalar_field.py). In the canonical (-,+,+,+)
# signature the dust trace is T = -rho and the static field equation is
# nabla^2 phi = -|alpha| rho: positive-density sources support a local field
# maximum (phi > 0 in wells, Rule 4). The ansatz's alpha_log < 0 instead
# encodes the screening direction: the observable excursion
# phi_rho = alpha_log * ln(rho/rho_T) is driven to zero as rho -> rho_T from
# below, so dense environments show a suppressed response. The magnitude
# |7.66e-3| was calibrated against laboratory metrology shift scales.
ALPHA_LOG = -7.66e-3             # Density-sector coupling (screening-direction sign; legacy ansatz)
BETA_GEOM = 1.50e-4              # Mass-sector geometric coupling (legacy ansatz)

# =============================================================================
# GALAXY-SCALE OBSERVABLE PARAMETERS
# =============================================================================

# Canonical galaxy-scale observable response coefficient (Paper 11).
# Declared channel-response benchmark (kappa_canonical in
# core/parameter_registry.yaml), applied corpus-wide without domain-specific
# refitting. It is NOT the raw clock ratio: the nested
# spectroscopic-to-Cepheid computation (TEP-H0 step_61) shows the clock
# channel contributes |kappa| of order unity at most (KAPPA_NESTED below),
# so the measured response is carried by the transfer/screening sector.
KAPPA_GAL = 9.6e5                # mag
KAPPA_GAL_UNCERTAINTY = 4.0e5    # mag

# Null-channel bound on the raw conformal clock ratio (TEP-H0 step_61). The
# spectroscopic and Cepheid clocks share each galactic potential, so the
# common rate cancels; the residual core-to-disk contrast difference gives
# kappa_nested = -7.5e-4 mag (delta_mu ~ 3e-10 mag). Alternative
# reference-clock conventions move the sign but leave |kappa| of order
# unity at most -- the raw clock channel cannot supply the measured
# kappa_Cep and the observed response is not a clock-ratio artifact.
KAPPA_NESTED = -7.5e-4           # mag

# Stellar evolution index (M/L ~ t^n from stellar isochrones)
ALPHA_NUCLEAR = 0.7

# Fixed response-transfer exponent (Paper 12). The magnitude-to-kernel
# transfer K_gal = kappa * ln(10) / (2.5 * n_ref) uses this fixed
# convention at every redshift; it is NOT the stellar-population
# exponent and must not be renormalised by population-dependent n_SPS.
RESPONSE_KERNEL_N_REF = 0.7

# Default stellar-population-synthesis M/L evolution index, applied only
# when the response is propagated into a stellar-mass correction.
SPS_ML_EXPONENT_DEFAULT = ALPHA_NUCLEAR

# Reference halo mass for potential calculations
LOG_MH_REF = 12.0

# Dimensionless virial potential Phi/c^2 for 10^12 Msun halo at z=0
PHI_REF_0 = 1.6e-7

# Reference redshift for chronological enhancement
Z_REF = 5.5

# =============================================================================
# CANONICAL SCALAR SECTOR (master-action closure, Paper 0 steps 03/27/50)
# =============================================================================
# All field values below are dimensionless u = phi / M_Pl (reduced Planck
# convention); convert to GeV via M_PL_REDUCED_GEV.

M_PL_REDUCED_GEV = 2.435e18      # reduced Planck mass in GeV
HBAR_C_GEV_M = 0.1973269804e-15  # hbar*c in GeV*m
KG_TO_GEV = 5.60958885e26        # 1 kg in GeV (natural units)
# 1 g/cm^3 = 1000 kg/m^3 -> GeV^4 :  kg -> GeV, 1/m^3 -> (hbar c)^3 GeV^3
G_CM3_TO_GEV4 = 1000.0 * KG_TO_GEV * HBAR_C_GEV_M**3

# Observed cosmic drift rate used to set the kinetic screening scale. This is
# the *measured* Hubble drift (the shear floor of the landscape), distinct from
# the TEP-corrected ladder output H0 = 66.65 reported in parameter_registry.
H0_DRIFT_KM_S_MPC = 70.0
# Transition shear/acceleration scale of the screening operator, derived from
# the kinetic completion P(X) = X - V + X|X|/Lambda^4 with Lambda^4 = M_Pl^2 H0^2
# (zero free parameters; Paper 0 Appendix E, step_27):
#     g_t = c H0 / (2 |beta_A|)
G_T_TRANSITION = C_LIGHT * (H0_DRIFT_KM_S_MPC * 1e3 / MPC_TO_M) / (2.0 * abs(BETA_A))
# Kinetic scale Lambda = sqrt(M_Pl H0) ~ 1.9 meV (the dark-energy scale).
# Conversion: H0 [s^-1] -> GeV via hbar = 6.582e-25 GeV s.
HBAR_GEV_S = 6.582119569e-25
LAMBDA_KINETIC_MEV = 1e12 * np.sqrt(
    M_PL_REDUCED_GEV * (H0_DRIFT_KM_S_MPC * 1e3 / MPC_TO_M) * HBAR_GEV_S
)

# Quartic self-coupling on the matter-hosting weak-field branch of the master
# potential V(u) = V_matter(u) e^{-(u/u_s)^4} + V_0 e^{-(u_s/u)^4}.
# Reference normalization (step_03, fixes lambda_c/R_T ~= 4.07 for Earth);
# the Cassini-compatible branch sits 1e5 higher (step_01 solar scan).
LAMBDA_QUARTIC_REF = 7.526e-71
LAMBDA_QUARTIC_CASSINI = 7.526e-66

# Master-potential family parameters (dimensionless-field convention).
# u_s is the plateau knee (manuscript nominal ~10; closure probes scanned
# 10--50); V_0 is the deep-well floor in M_Pl^4 units (fiducial).
U_S_TRANSITION = 10.0
V0_PLATEAU = 0.3

# Equivalent potential-normalization scale of the quartic branch,
# Lambda_V = lambda_quartic^(1/4) * M_Pl. Derived, not an independent
# parameter (Paper 6 App. C; TEP-UCD step_14 scale audit).
LAMBDA_V_GEV = LAMBDA_QUARTIC_CASSINI**0.25 * M_PL_REDUCED_GEV      # ~127.5 GeV
LAMBDA_V_REF_GEV = LAMBDA_QUARTIC_REF**0.25 * M_PL_REDUCED_GEV      # ~7.2 GeV

# Self-quenching density marker: the order-of-magnitude density at which the
# quartic equilibrium crosses u_min ~ 1, rho_sat ~ lambda * M_Pl^4
# (step_14: 6.1e25 g/cm^3; with the e^{-u} source factor the strict u_min = 1
# crossing sits at e * lambda * M_Pl^4 ~ 1.7e26). Temporal-well regime
# (Paper 28), NOT a terrestrial scale — terrestrial saturation is the
# Compton-resolution crossover lambda_c(rho) ~ R_T(M).
RHO_SAT_G_CM3 = (
    LAMBDA_QUARTIC_CASSINI * M_PL_REDUCED_GEV**4 / G_CM3_TO_GEV4
)

# Compact-object lower bound on the quartic normalization: neutron-star
# interiors (rho_NS ~ 2e14 g/cm^3) must remain below the master-potential
# knee, u_min(rho_NS) < u_s, i.e. lambda > rho_NS e^{-u_s} / (u_s^3 M_Pl^4).
# Step_14: 1.11e-84; the Cassini branch exceeds it by ~6.8e18 while the
# excluded rho_T-normalized branch (lambda ~ 2.45e-90) fails by ~5e5.
LAMBDA_QUARTIC_NS_BOUND = (
    (2e14 * G_CM3_TO_GEV4) * np.exp(-U_S_TRANSITION)
    / (U_S_TRANSITION**3 * M_PL_REDUCED_GEV**4)
)

# Disformal envelope family B(u) = B0 * u^2/(1+u^2) * exp(-u^4/2), admissible
# branch B0 >= 0. Canonical benchmark B0 = +1 (Paper 28 normalization);
# GW170817 path bound on the galactic profile gives B0 <~ 77.7 (step_50).
B0_DISFORMAL_CANONICAL = 1.0
B0_DISFORMAL_GW170817_MAX = 77.7
