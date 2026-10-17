"""Stationary metric null propagation, separate from endpoint clock rates.

Metric components use x^0=ct and a (-,+,+,+) signature. The returned
coordinate-time increment is seconds per kilometre of coordinate path.
Conformal cancellation is tested by constructing the scaled metric itself;
no conformal refractive-index approximation is introduced.
"""
import numpy as np
C_KM_S=299792.458

def null_seconds_per_km(metric,spatial_direction):
    metric=np.asarray(metric,dtype=float)
    direction=np.asarray(spatial_direction,dtype=float)
    if metric.shape!=(4,4) or direction.shape!=(3,):
        raise ValueError('A 4x4 metric and three-component direction are required')
    if not np.all(np.isfinite(metric)) or not np.all(np.isfinite(direction)):
        raise ValueError('Nonfinite metric or direction')
    scale=np.max(np.abs(metric))
    if scale==0:raise ValueError('Degenerate metric')
    metric=metric/scale
    if not np.allclose(metric,metric.T,atol=1e-14,rtol=1e-14):
        raise ValueError('Metric must be symmetric')
    norm=np.linalg.norm(direction)
    if norm==0:raise ValueError('Nonzero spatial direction required')
    direction=direction/norm
    a=metric[0,0];b=np.dot(metric[0,1:],direction)
    d=direction@metric[1:,1:]@direction
    if a>=0 or d<=0:raise ValueError('Stationary timelike coordinate and spacelike path required')
    discriminant=b*b-a*d
    if discriminant<=0:raise ValueError('No real null ray in this direction')
    root=np.sqrt(discriminant)
    # Positive future-directed root. Rationalize where b is negative to
    # avoid cancellation between b and sqrt(b^2-a*d).
    coordinate_slope=(b+root)/(-a) if b>=0 else d/(root-b)
    return float(coordinate_slope/C_KM_S)

def stationary_clock_rate(metric):
    """d(tau)/dt for a stationary endpoint; conformal scaling changes it."""
    metric=np.asarray(metric,dtype=float)
    if metric.shape!=(4,4) or not np.isfinite(metric[0,0]) or metric[0,0]>=0:
        raise ValueError('Stationary endpoint must be timelike')
    return float(np.sqrt(-metric[0,0]))


def unstretched_temporal_metric(metric, u, beta_A=-1.0):
    """Construct the unstretched-space temporal metric
        g_tilde_mn = g_mn - (A^2(u) - 1) U_m U_n
    in the cosmic temporal rest frame U_m = (-sqrt(-g_00), 0, 0, 0), where
    clock rate scales by A(u) = exp(beta_A u) while spatial rulers g_ij do
    not stretch. Preserves local Lorentz invariance dl_proper/dtau_proper = c
    at every point while retaining intermediate temporal depth u(x) along
    null rays in coordinate time.
    """
    metric = np.asarray(metric, dtype=float).copy()
    if metric.shape != (4, 4) or not np.all(np.isfinite(metric)) or metric[0, 0] >= 0:
        raise ValueError('Valid 4x4 stationary timelike metric required')
    a2 = float(np.exp(2.0 * float(beta_A) * float(u)))
    metric[0, 0] *= a2
    return metric


def unstretched_refractive_index(u, phi_N_over_c2=0.0, beta_A=-1.0):
    """Effective coordinate refractive index n(x) = c dt / dl along a null
    ray in the unstretched-space temporal metric:
        n(x) = exp(-beta_A u(x)) * sqrt((1 - 2 Phi_N/c^2) / (1 + 2 Phi_N/c^2))
             ~ 1 - 2 Phi_N(x)/c^2 - beta_A u(x).
    """
    u_arr = np.asarray(u, dtype=float)
    phi_arr = np.asarray(phi_N_over_c2, dtype=float)
    return np.exp(-beta_A * u_arr) * np.sqrt(
        (1.0 - 2.0 * phi_arr) / (1.0 + 2.0 * phi_arr)
    )

