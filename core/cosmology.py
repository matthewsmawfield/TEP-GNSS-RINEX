#!/usr/bin/env python3
"""
TEP Cosmological Distance Models
=================================

Version: TEP v0.9 (Jakarta)

Provides the luminosity-distance--redshift relation for competing hypotheses,
used by standard-siren and cosmology analyses across the TEP corpus.

All TEP papers should import from this module to ensure consistent
cosmological calculations.

Conventions
-----------
- H0 in km/s/Mpc, distances in Mpc.
- Redshift z is the cosmological (CMB-frame, Hubble-flow) redshift,
  obtained independently of the GW distance.

Native TEP background (Paper 26 / TEP-C0, Paper 18 / TEP-HC)
------------------------------------------------------------
The hi_class ``tep_mode`` implementation uses:

    f_T(z) = ln(1+z) * exp(-(z/z_T)^n_T)
    A(z) = exp(epsilon_T * ln(1+z) * S(z)),  S(z) = exp(-(z/z_T)^n_T)
    M(z) = A(z) / (1 - alpha_A(z))
    H_TEP(z) = H_LCDM(z) * M(z)

Functions ``f_T_suppression``, ``f_T``, ``conformal_factor_native``,
``alpha_A_native``, and ``jordan_frame_M`` mirror ``background.c`` in
``external/patches/hiclass_tep_native.patch`` (TEP-HC).
"""

import numpy as np
from scipy.integrate import quad
from . import constants as tep_const

C_KM_S = 299792.458  # speed of light, km/s

_Z_CAP_FACTOR = 3.0
_A_MIN = 0.1


def _as_array(z):
    return np.atleast_1d(np.asarray(z, dtype=float))


def _scalar_result(z, arr):
    return arr if np.ndim(z) else float(arr[0])


def _validate_tep_params(z_T, n_T, epsilon_T=None):
    if not np.isfinite(z_T) or z_T <= 0.0:
        raise ValueError("z_T must be finite and strictly positive")
    if not np.isfinite(n_T) or n_T <= 0.0:
        raise ValueError("n_T must be finite and strictly positive")
    if epsilon_T is not None and (not np.isfinite(epsilon_T) or epsilon_T < 0.0):
        raise ValueError("epsilon_T must be finite and non-negative")


def _z_effective(z, z_T):
    z_arr = _as_array(z)
    cap = z_T * _Z_CAP_FACTOR
    return np.where(z_arr > cap, cap, z_arr)


def f_T_suppression(z, z_T, n_T):
    """
    Early-universe suppression factor S(z) = exp(-(z/z_T)^n_T).

    Matches hi_class ``tep_f_transition`` (suppression only; the full
    transition function is ``f_T`` = ln(1+z) * S(z)).
    """
    _validate_tep_params(z_T, n_T)
    z_arr = _as_array(z)
    out = np.zeros_like(z_arr, dtype=float)
    mask = z_arr > 0.0
    if np.any(mask):
        z_eff = _z_effective(z_arr[mask], z_T)
        out[mask] = np.exp(-np.power(z_eff / z_T, n_T))
    return _scalar_result(z, out)


def f_T(z, z_T, n_T):
    """
    Authoritative TEP transition function (Paper 26 / TEP-HC native tep_mode):

        f_T(z) = ln(1+z) * exp(-(z/z_T)^n_T)
    """
    _validate_tep_params(z_T, n_T)
    z_arr = _as_array(z)
    out = np.zeros_like(z_arr, dtype=float)
    mask = z_arr > 0.0
    if np.any(mask):
        S = f_T_suppression(z_arr[mask], z_T, n_T)
        out[mask] = np.log(1.0 + z_arr[mask]) * S
    return _scalar_result(z, out)


def conformal_factor_native(z, epsilon_T, z_T, n_T):
    """
    Covariant conformal factor A(z) = exp(epsilon_T * ln(1+z) * S(z)).

    Matches hi_class ``tep_gamma_factor``.
    """
    _validate_tep_params(z_T, n_T, epsilon_T)
    if epsilon_T == 0.0:
        z_arr = _as_array(z)
        return _scalar_result(z, np.ones_like(z_arr, dtype=float))
    z_arr = _as_array(z)
    out = np.ones_like(z_arr, dtype=float)
    mask = z_arr > 0.0
    if np.any(mask):
        S = f_T_suppression(z_arr[mask], z_T, n_T)
        out[mask] = np.exp(epsilon_T * np.log(1.0 + z_arr[mask]) * S)
        out[mask] = np.maximum(out[mask], _A_MIN)
    return _scalar_result(z, out)


def alpha_A_native(z, epsilon_T, z_T, n_T):
    """
    Jordan-frame coupling alpha_A = d ln A / d ln(1+z).

    Matches hi_class ``tep_M_factor`` intermediate.
    """
    _validate_tep_params(z_T, n_T, epsilon_T)
    z_arr = _as_array(z)
    out = np.zeros_like(z_arr, dtype=float)
    if epsilon_T == 0.0:
        return _scalar_result(z, out)
    mask = z_arr > 0.0
    if np.any(mask):
        z_pos = z_arr[mask]
        S = f_T_suppression(z_pos, z_T, n_T)
        dS = np.zeros_like(S)
        deriv_mask = (z_pos > 1e-10) & (z_pos <= z_T * _Z_CAP_FACTOR)
        if np.any(deriv_mask):
            z_d = z_pos[deriv_mask]
            S_d = S[deriv_mask]
            dS[deriv_mask] = (
                -S_d * n_T * np.power(z_d / z_T, n_T - 1.0) / z_T
            )
        out[mask] = -epsilon_T * (S + (1.0 + z_pos) * np.log(1.0 + z_pos) * dS)
    return _scalar_result(z, out)


def jordan_frame_M(z, epsilon_T, z_T, n_T):
    """
    Jordan-frame expansion modifier M(z) = A(z) / (1 - alpha_A(z)).

    Matches hi_class ``tep_M_factor``; H_TEP = M * H_LCDM.
    """
    _validate_tep_params(z_T, n_T, epsilon_T)
    if epsilon_T == 0.0:
        z_arr = _as_array(z)
        return _scalar_result(z, np.ones_like(z_arr, dtype=float))
    z_arr = _as_array(z)
    out = np.ones_like(z_arr, dtype=float)
    mask = z_arr > 0.0
    if np.any(mask):
        z_pos = z_arr[mask]
        A = conformal_factor_native(z_pos, epsilon_T, z_T, n_T)
        alpha = alpha_A_native(z_pos, epsilon_T, z_T, n_T)
        denom = 1.0 - alpha
        if np.any(denom <= 0.0):
            raise ValueError("Unphysical native TEP modifier: 1 - alpha_A must be positive")
        out[mask] = A / denom
    return _scalar_result(z, out)


def hubble_modifier_native(z, epsilon_T, z_T, n_T):
    """Alias for ``jordan_frame_M`` (H_TEP / H_LCDM)."""
    return jordan_frame_M(z, epsilon_T, z_T, n_T)


def E_of_z(z, Om=0.315):
    """Dimensionless Hubble parameter E(z) = H(z)/H0 for a flat universe."""
    if Om < 0.0 or Om > 1.0:
        raise ValueError("Om must be in the interval [0, 1] for a flat LCDM model")
    return np.sqrt(Om * (1.0 + z) ** 3 + (1.0 - Om))


def comoving_distance(z, H0, Om=0.315):
    """Line-of-sight comoving distance D_C(z) in Mpc (scalar z)."""
    if H0 <= 0.0:
        raise ValueError("H0 must be positive")
    if z < 0.0:
        raise ValueError("Redshift z must be non-negative")
    integral, _ = quad(lambda zp: 1.0 / E_of_z(zp, Om), 0.0, z)
    return (C_KM_S / H0) * integral


def luminosity_distance_lcdm(z, H0, Om=0.315):
    """
    LambdaCDM luminosity distance in Mpc for a flat universe.

    Accepts scalar or array-like z.
    """
    if H0 <= 0.0:
        raise ValueError("H0 must be positive")
    z_arr = np.atleast_1d(np.asarray(z, dtype=float))
    if np.any(z_arr < 0.0):
        raise ValueError("Redshift z must be non-negative")
    dc = np.array([comoving_distance(zi, H0, Om) for zi in z_arr])
    dl = (1.0 + z_arr) * dc
    return dl if np.ndim(z) else float(dl[0])


def phenomenological_redshift_factor(z, phi0, n=1.0, beta_A=tep_const.BETA_A):
    """
    Phenomenological redshift fit template.

    WARNING: This is a phenomenological redshift-dependent ansatz, not the fundamental
    scalar-field conformal factor A(phi) defined by the Jakarta axioms. It replaces
    phi/M_Pl with phi0 * (1+z)^n, and there is no scalar-field equation of motion that
    produces this exactly. Do not equate this A(z) with the true A(phi).

    Accepts scalar or array-like z.
    """
    z_arr = np.asarray(z, dtype=float)
    return np.exp(beta_A * phi0 * np.power(1.0 + z_arr, n))


def conformal_factor_derivative(z, phi0, n=1.0, beta_A=tep_const.BETA_A):
    """
    Analytical derivative dA/dz of the TEP conformal factor.

    A(z) = exp(beta_A * phi0 * (1+z)^n)
    dA/dz = A(z) * beta_A * phi0 * n * (1+z)^(n-1)

    For n == 0 the conformal factor is constant and the derivative is zero.
    Accepts scalar or array-like z.
    """
    z_arr = np.asarray(z, dtype=float)
    a_z = phenomenological_redshift_factor(z_arr, phi0, n, beta_A)
    if n == 0:
        return np.zeros_like(a_z)
    return a_z * beta_A * phi0 * n * np.power(1.0 + z_arr, n - 1.0)


def alpha_a_z(z, phi0, n=1.0, beta_A=tep_const.BETA_A):
    """
    Jordan-frame alpha_A = d ln A / d ln a_J.

    With A(z) = exp[beta_A * phi0 * (1+z)^n] and 1+z = a_J^-1,
    alpha_A = - beta_A * phi0 * n * (1+z)^n.
    """
    z_arr = np.asarray(z, dtype=float)
    return -beta_A * phi0 * n * np.power(1.0 + z_arr, n)


def hubble_modifier_tep_c0(z, phi0, n=1.0, beta_A=tep_const.BETA_A):
    """
    TEP-C0 Jordan-frame expansion modifier H_J / H_LCDM.

    This is the exact conformal-frame chain-rule factor for physical
    Jordan-frame redshift:

        H_J = A / (1 - alpha_A) * H_LCDM.
    """
    alpha = alpha_a_z(z, phi0, n, beta_A)
    denom = 1.0 - alpha
    if np.any(denom <= 0.0):
        raise ValueError("Unphysical TEP-C0 modifier: 1 - alpha_A must be positive")
    return phenomenological_redshift_factor(z, phi0, n, beta_A) / denom


def luminosity_distance_tep(z, H0, phi0, n=1.0, beta_A=tep_const.BETA_A, Om=0.315):
    """
    Endpoint-only TEP luminosity distance in Mpc: d_L^TEP = A(z) * d_L^LCDM.

    Accepts scalar or array-like z.
    """
    return phenomenological_redshift_factor(z, phi0, n, beta_A) * luminosity_distance_lcdm(z, H0, Om)


def luminosity_distance_tep_c0_jordan(
    z,
    H0,
    phi0,
    n=1.0,
    beta_A=tep_const.BETA_A,
    Om=0.315,
    include_gw_endpoint=True,
):
    """
    TEP-C0 Jordan-frame GW luminosity distance in Mpc.

    The matter-frame distance integral uses physical Jordan-frame redshift and
    H_J = [A / (1 - alpha_A)] H_LCDM. For standard-siren comparisons, the
    default also applies the GW endpoint factor A(z).
    """
    if H0 <= 0.0:
        raise ValueError("H0 must be positive")
    z_arr = np.atleast_1d(np.asarray(z, dtype=float))
    if np.any(z_arr < 0.0):
        raise ValueError("Redshift z must be non-negative")

    def integrand(zp):
        modifier = hubble_modifier_tep_c0(zp, phi0, n, beta_A)
        return 1.0 / (E_of_z(zp, Om) * modifier)

    integrals = np.array([
        quad(integrand, 0.0, zi)[0]
        for zi in z_arr
    ])
    dl = (1.0 + z_arr) * (C_KM_S / H0) * integrals
    if include_gw_endpoint:
        dl = phenomenological_redshift_factor(z_arr, phi0, n, beta_A) * dl
    return dl if np.ndim(z) else float(dl[0])


def redshift_at_distance(dl_target, H0, Om=0.315, z_max=3.0):
    """
    Invert the LambdaCDM d_L(z) relation to find z given d_L (Mpc).

    Used only for diagnostics / plotting limits, never to assign event
    redshifts in the likelihood.
    """
    from scipy.optimize import brentq
    if dl_target < 0.0:
        raise ValueError("Target luminosity distance must be non-negative")
    if dl_target == 0.0:
        return 0.0
    if H0 <= 0.0:
        raise ValueError("H0 must be positive")
    f = lambda z: luminosity_distance_lcdm(z, H0, Om) - dl_target
    if f(1e-4) > 0:
        return 1e-4
    while f(z_max) < 0.0 and z_max < 100.0:
        z_max *= 2.0
    if f(z_max) < 0.0:
        raise ValueError(f"Could not bracket redshift for d_L={dl_target} Mpc")
    return brentq(f, 1e-4, z_max)
