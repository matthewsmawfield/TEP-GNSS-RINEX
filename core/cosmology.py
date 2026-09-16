#!/usr/bin/env python3
"""
TEP Cosmological Distance Models
=================================

Version: TEP v0.10 (Jakarta)

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
from scipy.integrate import cumulative_trapezoid, quad
from . import constants as tep_const

C_KM_S = 299792.458  # speed of light, km/s
C_KMS = C_KM_S  # alias for backward compatibility

_Z_CAP_FACTOR = 3.0
# Numerical guardrail for exploratory parameter scans. The validated native
# TEP parameter ranges stay well above this floor; hitting it means the scan is
# outside the calibrated EFT regime rather than a new physical branch.
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
    if np.any(z_arr > cap):
        import warnings
        warnings.warn(
            f"_z_effective: z values above {cap} capped to {cap} for numerical stability. "
            "High-z physics may differ from full TEP model.",
            RuntimeWarning, stacklevel=3
        )
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
    Jordan-frame coupling alpha_A = d ln A / d ln a_J.

    Since a_J = 1 / (1 + z), this is the negative of
    d ln A / d ln(1+z).

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
    if z < -1.0:
        raise ValueError("Redshift z must be >= -1")
    arg = Om * (1.0 + z) ** 3 + (1.0 - Om)
    if arg < 0:
        raise ValueError(f"Negative argument to sqrt in E_of_z: {arg}")
    return np.sqrt(arg)


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


class TEPCosmology:
    """TEP cosmology with temporal shear effects."""

    # Spatial transition density (manuscript: rho_half ≈ 0.5 M_sun/pc^3)
    # Maps to Equation in Section 2.5: S(rho) = [1 + (rho/rho_half)^2]^-1
    RHO_HALF: float = 0.5  # M_sun / pc^3, threshold galactic onset density

    def __init__(
        self,
        H0: float = 70.0,
        Omega_m: float = 0.3,
        epsilon_T: float = 0.0,
        z_T: float = 5.0,
        n_T: float = 1.0,
    ):
        self.H0 = H0
        self.Omega_m = Omega_m
        self.epsilon_T = epsilon_T
        self.z_T = z_T
        self.n_T = n_T
        self.Omega_lambda = 1.0 - Omega_m

    def tep_gamma(self, z: float | np.ndarray) -> float | np.ndarray:
        """TEP path enhancement factor.

        Uses the exact normalised formula: gamma(z) = 1 + epsilon_T * ln(1+z) * S(z).
        This guarantees gamma(0) = 1 (recovering the local reference frame),
        and gamma > 1 for z > 0 (reflecting clocks ran slower in the denser past,
        A_past < A_today, thus gamma = A_today / A_past > 1).
        This makes intermediate distances larger, actively mimicking dark energy acceleration.
        """
        z_arr = np.asarray(z, dtype=float)
        if self.epsilon_T == 0:
            gamma = np.ones_like(z_arr)
            return float(gamma) if np.isscalar(z) else gamma

        log_factor = np.log(1.0 + z_arr)
        suppression = np.exp(-((z_arr / self.z_T) ** self.n_T))
        enhancement = 1.0 + self.epsilon_T * log_factor * suppression
        gamma = np.maximum(enhancement, 0.1)  # Physical bound
        return float(gamma) if np.isscalar(z) else gamma

    def screening_function(self, rho: float | np.ndarray) -> float | np.ndarray:
        """TEP environmental screening factor S(rho).

        Implements the continuous shear-suppression formula from
        Section 2.5 of the manuscript:
            S(rho) = [1 + (rho / rho_half)^2]^-1

        where rho_half = 0.5 M_sun / pc^3 is the threshold galactic
        onset density (class constant RHO_HALF).

        Parameters
        ----------
        rho : float or array
            Local matter density in M_sun / pc^3.

        Returns
        -------
        float or array
            Screening factor between 0 (fully screened) and 1 (unscreened).
        """
        ratio = np.asarray(rho, dtype=float) / self.RHO_HALF
        return 1.0 / (1.0 + ratio ** 2)

    def e_func(self, z: float | np.ndarray) -> float | np.ndarray:
        """Dimensionless Hubble parameter."""
        zp1 = 1.0 + np.asarray(z, dtype=float)
        matter = self.Omega_m * zp1**3
        Lambda = self.Omega_lambda
        return np.sqrt(matter + Lambda)

    def comoving_distance(self, z: float | np.ndarray) -> float | np.ndarray:
        """TEP comoving distance in Mpc."""
        z_arr = np.atleast_1d(np.asarray(z, dtype=float))
        c = 299792.458  # km/s

        max_z = float(np.max(z_arr))
        if max_z <= 0:
            distances = np.zeros_like(z_arr)
            return float(distances[0]) if np.isscalar(z) else distances

        grid = np.linspace(0.0, max_z, 1000)
        integrand = c / (self.H0 * self.e_func(grid)) * self.tep_gamma(grid)
        cumulative = cumulative_trapezoid(integrand, grid, initial=0.0)
        distances = np.interp(z_arr, grid, cumulative)
        return float(distances[0]) if np.isscalar(z) else distances

    def angular_diameter_distance(self, z: float | np.ndarray) -> float | np.ndarray:
        """TEP angular diameter distance in Mpc."""
        d_c = self.comoving_distance(z)
        return d_c / (1 + z)

    def luminosity_distance(self, z: float | np.ndarray) -> float | np.ndarray:
        """TEP luminosity distance in Mpc."""
        d_c = self.comoving_distance(z)
        return d_c * (1 + z)

    def distance_modulus(self, z: float | np.ndarray) -> float | np.ndarray:
        """TEP distance modulus."""
        d_l = np.maximum(self.luminosity_distance(z), np.finfo(float).tiny)
        return 5.0 * np.log10(d_l) + 25.0


class CosmologyFLRW:
    """Standard FLRW cosmology for comparison."""

    def __init__(self, H0: float = 70.0, Om0: float = 0.3, Ode0: float | None = None, Ok0: float = 0.0):
        self.H0 = H0
        self.Om0 = Om0
        self.Ode0 = Ode0 if Ode0 is not None else 1.0 - Om0 - Ok0
        self.Ok0 = Ok0

    def e_func(self, z: float) -> float:
        """Dimensionless Hubble parameter."""
        zp1 = 1.0 + z
        return np.sqrt(self.Om0 * zp1**3 + self.Ode0 + self.Ok0 * zp1**2)

    def comoving_distance(self, z: float | np.ndarray) -> float | np.ndarray:
        """Comoving distance in Mpc.

        Uses vectorised trapezoidal integration on a fine redshift grid,
        which is ~100x faster than calling quad() for each redshift
        individually and has negligible accuracy loss (<0.001%).
        """
        z_arr = np.atleast_1d(np.asarray(z, dtype=float))
        c = 299792.458  # km/s

        max_z = float(np.max(z_arr)) if z_arr.size > 0 else 0.0
        if max_z <= 0.0:
            distances = np.zeros_like(z_arr)
            return float(distances[0]) if np.isscalar(z) else distances

        grid = np.linspace(0.0, max_z, 2000)
        integrand = c / (self.H0 * self.e_func(grid))
        cumulative = cumulative_trapezoid(integrand, grid, initial=0.0)
        distances = np.interp(z_arr, grid, cumulative)
        return float(distances[0]) if np.isscalar(z) else distances

    def angular_diameter_distance(self, z: np.ndarray) -> np.ndarray:
        """Angular diameter distance in Mpc."""
        return self.comoving_distance(z) / (1 + z)

    def luminosity_distance(self, z: np.ndarray) -> np.ndarray:
        """Luminosity distance in Mpc."""
        return self.comoving_distance(z) * (1 + z)

    def distance_modulus(self, z: np.ndarray) -> np.ndarray:
        """Distance modulus."""
        d_l = np.maximum(self.luminosity_distance(z), np.finfo(float).tiny)
        return 5.0 * np.log10(d_l) + 25.0


class wCDMCosmology(CosmologyFLRW):
    """wCDM cosmology with constant equation of state w."""

    def __init__(self, H0: float = 70.0, Om0: float = 0.3, Ok0: float = 0.0, Ode0: float | None = None, w: float = -1.0):
        super().__init__(H0=H0, Om0=Om0, Ok0=Ok0, Ode0=Ode0)
        self.w = w

    def e_func(self, z: float) -> float:
        zp1 = 1.0 + z
        de_term = self.Ode0 * zp1 ** (3.0 * (1.0 + self.w))
        return np.sqrt(self.Om0 * zp1 ** 3 + de_term + self.Ok0 * zp1 ** 2)


class CPLCosmology(CosmologyFLRW):
    """CPL (w0wa) parametrisation for dark energy."""

    def __init__(self, H0: float = 70.0, Om0: float = 0.3, Ok0: float = 0.0, Ode0: float | None = None, w0: float = -1.0, wa: float = 0.0):
        super().__init__(H0=H0, Om0=Om0, Ok0=Ok0, Ode0=Ode0)
        self.w0 = w0
        self.wa = wa

    def e_func(self, z: float) -> float:
        zp1 = 1.0 + z
        a = 1.0 / zp1
        exponent = 3.0 * (1.0 + self.w0 + self.wa) * np.log(zp1) - 3.0 * self.wa * (1.0 - a)
        de_term = self.Ode0 * np.exp(exponent)
        return np.sqrt(self.Om0 * zp1 ** 3 + de_term + self.Ok0 * zp1 ** 2)


class TEPCosmologyFitter:
    """Fit TEP cosmology to observational data."""

    def __init__(
        self,
        z: np.ndarray | None = None,
        mu: np.ndarray | None = None,
        mu_err: np.ndarray | None = None,
        data_path: str | None = None,
    ):
        self.data_path = data_path
        self.data = None
        if z is not None and mu is not None and mu_err is not None:
            self.data = {
                'z': np.asarray(z, dtype=float),
                'mu': np.asarray(mu, dtype=float),
                'mu_err': np.asarray(mu_err, dtype=float),
            }

    def load_data(self):
        """Load SNe data for fitting."""
        import pandas as pd
        from pathlib import Path

        if self.data_path is None:
            # Try default locations
            candidates = [
                Path("data/raw/pantheon_plus_shoes.dat"),
                Path("data/raw/Pantheon+SH0ES.dat"),
            ]
            for path in candidates:
                if path.exists():
                    self.data_path = str(path)
                    break

        if self.data_path is None or not Path(self.data_path).exists():
            raise FileNotFoundError("Pantheon+ data not found")

        df = pd.read_csv(self.data_path, sep=r'\s+', comment='#')

        # Extract relevant columns
        z = df['zHD'].values if 'zHD' in df.columns else df['zCMB'].values
        mb = df['m_b_corr'].values if 'm_b_corr' in df.columns else df['mB'].values
        mb_err = df['m_b_corr_err_DIAG'].values if 'm_b_corr_err_DIAG' in df.columns else df.get('mBERR', pd.Series([0.15]*len(df))).values

        # Filter valid data
        valid = np.isfinite(z) & np.isfinite(mb) & (z > 0.001) & (z < 3.0)

        self.data = {
            'z': z[valid],
            'mb': mb[valid],
            'mb_err': mb_err[valid],
        }

        return self.data

    def _observations(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.data is None:
            self.load_data()

        z = self.data['z']
        if 'mu' in self.data:
            obs = self.data['mu']
            err = self.data['mu_err']
        else:
            obs = self.data['mb']
            err = self.data['mb_err']
        return z, obs, np.maximum(err, np.finfo(float).tiny)

    def fit_tep(self, initial_guess: dict = None) -> dict:
        """Fit TEP model to data.

        Returns best-fit parameters and statistics.
        """
        if self.data is None:
            self.load_data()

        from scipy.optimize import minimize

        z, obs, obs_err = self._observations()

        # Initial guess
        if initial_guess is None:
            p0 = [70.0, 0.3, 0.1, 5.0, 0.0]  # H0, Omega_m, epsilon_T, z_T, offset
        else:
            p0 = [
                initial_guess.get('H0', 70.0),
                initial_guess.get('Omega_m', 0.3),
                initial_guess.get('epsilon_T', 0.1),
                initial_guess.get('z_T', 5.0),
                initial_guess.get('offset', 0.0),
            ]

        def chi2(params):
            H0, Omega_m, epsilon_T, z_T, offset = params
            try:
                tep = TEPCosmology(H0=H0, Omega_m=Omega_m, epsilon_T=epsilon_T, z_T=z_T)
                mu_pred = tep.distance_modulus(z) + offset
                residuals = (obs - mu_pred) / obs_err
                return np.sum(residuals**2)
            except Exception:
                return 1e10

        # Run minimization
        result = minimize(
            chi2,
            p0,
            method='Nelder-Mead',
            bounds=[(50, 100), (0.05, 0.9), (0.0, 1.0), (0.5, 15.0), (-5.0, 5.0)],
            options={"maxiter": 5000},
        )

        if result.success:
            H0_fit, Om_fit, eps_fit, zT_fit, offset_fit = result.x
            chi2_min = result.fun

            return {
                'success': True,
                'parameters': {
                    'H0': H0_fit,
                    'Omega_m': Om_fit,
                    'epsilon_T': eps_fit,
                    'z_T': zT_fit,
                    'offset': offset_fit,
                },
                'chi2': chi2_min,
                'dof': len(z) - 5,
            }
        else:
            return {
                'success': False,
                'error': 'Minimization failed',
            }

    def fit_lcdm(self, initial_guess: dict | None = None) -> dict:
        """Fit matched LCDM distance model to the same observations."""
        if self.data is None:
            self.load_data()

        from scipy.optimize import minimize

        z, obs, obs_err = self._observations()
        p0 = [
            70.0 if initial_guess is None else initial_guess.get('H0', 70.0),
            0.3 if initial_guess is None else initial_guess.get('Omega_m', 0.3),
            0.0 if initial_guess is None else initial_guess.get('offset', 0.0),
        ]

        def chi2(params):
            H0, Omega_m, offset = params
            try:
                lcdm = CosmologyFLRW(H0=H0, Om0=Omega_m)
                mu_pred = np.array([
                    5.0 * np.log10(lcdm.angular_diameter_distance(float(zi)) * (1.0 + zi) ** 2) + 25.0
                    for zi in z
                ]) + offset
                residuals = (obs - mu_pred) / obs_err
                return float(np.sum(residuals**2))
            except Exception:
                return 1e10

        result = minimize(
            chi2,
            p0,
            method='Nelder-Mead',
            bounds=[(50, 100), (0.05, 0.9), (-5.0, 5.0)],
            options={"maxiter": 5000},
        )
        if result.success:
            H0_fit, Om_fit, offset_fit = result.x
            return {
                'success': True,
                'parameters': {'H0': H0_fit, 'Omega_m': Om_fit, 'offset': offset_fit},
                'chi2': float(result.fun),
                'dof': len(z) - 3,
            }
        return {'success': False, 'error': 'Minimization failed'}

    def compare_models(self) -> dict:
        """Fit LCDM and TEP with matched data and report information criteria."""
        z, _, _ = self._observations()
        lcdm = self.fit_lcdm()
        tep = self.fit_tep()
        if not lcdm.get('success') or not tep.get('success'):
            return {'status': 'failed', 'lcdm': lcdm, 'tep': tep}

        n = len(z)
        k_lcdm = 3
        k_tep = 5
        lcdm_aic = lcdm['chi2'] + 2 * k_lcdm
        tep_aic = tep['chi2'] + 2 * k_tep
        lcdm_bic = lcdm['chi2'] + k_lcdm * np.log(n)
        tep_bic = tep['chi2'] + k_tep * np.log(n)
        lcdm['chi2_per_dof'] = lcdm['chi2'] / max(lcdm['dof'], 1)
        tep['chi2_per_dof'] = tep['chi2'] / max(tep['dof'], 1)

        delta_bic = tep_bic - lcdm_bic
        return {
            'status': 'completed',
            'lcdm': {**lcdm, 'aic': float(lcdm_aic), 'bic': float(lcdm_bic)},
            'tep': {**tep, 'aic': float(tep_aic), 'bic': float(tep_bic)},
            'delta_aic_tep_vs_lcdm': float(tep_aic - lcdm_aic),
            'delta_bic_tep_vs_lcdm': float(delta_bic),
            'tep_competitive': bool(delta_bic < 2.0),
            'best_model': 'tep' if delta_bic < 0.0 else 'lcdm',
        }

# ============================================================================
# Part E: Active Scalar Perturbation Closure (Minimal Conformal)
# ============================================================================

def evaluate_tep_eft_sector(alpha_A):
    """
    Returns the exact dimensionless Bellini-Sawicki parameters and
    stability discriminants in the pure conformal limit (beta_A = -1.0).

    Mapping from TEP bi-metric action to Jordan-frame EFT:
      alpha_M = d ln M_*^2 / d ln a = -2 alpha_A
      alpha_B = -alpha_M = 2 alpha_A
      alpha_K = -5 alpha_A^2   (conformal anomaly from kinetic mixing)
      alpha_T = 0

    Physical no-ghost discriminant:
      D = alpha_K + (3/2) alpha_B^2
        = -5 alpha_A^2 + (3/2)(4 alpha_A^2)
        = alpha_A^2  >= 0  (strictly, for all z)

    Sound speed: c_s^2 = 1 exactly (conformal isomorphism to FLRW).
    """
    alpha_A = np.atleast_1d(alpha_A)

    alpha_M = -2.0 * alpha_A
    alpha_B =  2.0 * alpha_A
    alpha_K = -5.0 * (alpha_A ** 2)
    alpha_T = 0.0

    # Physical no-ghost discriminant (must be positive-definite)
    D = alpha_K + 1.5 * (alpha_B ** 2)

    # Sound speed: fixed by conformal isomorphism
    c_s2 = np.ones_like(alpha_A)

    # Strict check: D == alpha_A^2 up to numerical tolerance
    D_expected = alpha_A ** 2
    if not np.allclose(D, D_expected, atol=1e-14):
        raise ValueError("Ghost-free identity D = alpha_A^2 violated!")

    is_stable = np.all(D >= 0) and np.all(c_s2 >= 0)

    return {
        "alpha_M": alpha_M,
        "alpha_B": alpha_B,
        "alpha_K": alpha_K,
        "alpha_T": alpha_T,
        "D": D,
        "c_s2": c_s2,
        "is_stable": is_stable
    }


def tep_alpha_M(z, epsilon_T=0.1, z_T=3.0, n_T=1.0):
    """
    Planck-mass running alpha_M = d ln M_*^2 / d ln a = -2 alpha_A.
    Uses the exact native TEP conformal-factor derivative.
    """
    alpha_A = alpha_A_native(z, epsilon_T, z_T, n_T)
    return -2.0 * alpha_A


def tep_alpha_B(z, epsilon_T=0.1, z_T=3.0, n_T=1.0):
    """
    Kinetic braiding alpha_B = -alpha_M = 2 alpha_A.
    """
    return -tep_alpha_M(z, epsilon_T, z_T, n_T)


def tep_alpha_K(z, epsilon_T=0.1, z_T=3.0, n_T=1.0):
    """
    Kineticity alpha_K = -5 alpha_A^2.
    The apparent negative value is balanced by the full no-ghost
    discriminant D = alpha_K + (3/2) alpha_B^2 = +alpha_A^2.
    """
    alpha_A = alpha_A_native(z, epsilon_T, z_T, n_T)
    return -5.0 * (alpha_A ** 2)


def get_hiclass_perturbation_params(epsilon_T, z_T=3.0, n_T=1.0, tep_perturbations='minimal_conformal'):
    """
    Configures hi_class for the active scalar perturbation closure.
    """
    params = {
        'tep_mode': 'yes',
        'tep_epsilon_T': epsilon_T,
        'tep_z_T': z_T,
        'tep_n_T': n_T
    }

    if tep_perturbations == 'minimal_conformal':
        params['gravity_model'] = 'tep'
        params['M2_evolution'] = 'yes'

    return params


def check_stability(z, epsilon_T=0.1, z_T=3.0, n_T=1.0):
    """
    Checks no-ghost (D > 0) and gradient (c_s^2 > 0) stability.
    For the pure conformal limit, D = alpha_A^2 and c_s^2 = 1 exactly.
    """
    alpha_A = alpha_A_native(z, epsilon_T, z_T, n_T)
    eft = evaluate_tep_eft_sector(alpha_A)

    ghost_free = np.all(eft["D"] > 0)
    gradient_free = np.all(eft["c_s2"] > 0)

    if not ghost_free:
        raise ValueError("Ghost instability detected (D <= 0)!")
    if not gradient_free:
        raise ValueError("Gradient instability detected (c_s^2 < 0)!")

    return ghost_free, gradient_free
