#!/usr/bin/env python3
"""
TEP Scalar Field Solver
=======================

Version: TEP v0.9 (Jakarta)

Computes the dimensionless scalar field profile phi from local
mass distributions using the TEP lab-scale logarithmic model.

Phenomenological model (TEP Lab-Scale Logarithmic Form):
    phi = alpha_log * ln(rho_local / rho_c) + beta_geom * ln(M / M_ref)

where alpha_log is the lab-scale density coupling and
beta_geom is the geometric coupling. Both are dimensionless.

This formulation ensures phi remains O(10^-3) for laboratory
conditions, keeping A(phi) = exp(beta_A*phi) numerically stable
while reproducing the observed metrology shifts.

All TEP papers should import from this module to ensure consistent
scalar-field calculations across the corpus.
"""

import numpy as np
from . import constants as tep_const

RHO_C = tep_const.RHO_C
screening_factor = tep_const.screening_factor

# Lab-scale coupling constants (from TEP-NIST Paper 21)
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
    screen = screening_factor(density_g_cm3, RHO_C)

    # Density-sector contribution with screening
    phi_rho = alpha_log * np.log(density_g_cm3 / RHO_C) * screen

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
    screen = screening_factor(eff_density_g_cm3, RHO_C)
    phi_rho = alpha_log * np.log(eff_density_g_cm3 / RHO_C) * screen
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

    phi = alpha_log * ln(rho / rho_c) + beta_geom * ln(M / M_ref)

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

    phi_rho = alpha_log * np.log(np.asarray(density_g_cm3) / RHO_C)
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
