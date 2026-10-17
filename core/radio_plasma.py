"""First-order multi-link dispersive Doppler calibration.

Aligned fractional-frequency residuals have a common nondispersive component
and separate uplink/downlink plasma components proportional to 1/f^2.
Frequency/path/clock convention certification and measurement covariance
must precede application to observations. This operator alone is not a fit.
"""
import numpy as np

def mixing_matrix(uplink_hz,downlink_hz):
    up=np.asarray(uplink_hz,dtype=float);down=np.asarray(downlink_hz,dtype=float)
    if up.ndim!=1 or up.shape!=down.shape or not len(up):
        raise ValueError('Matched one-dimensional frequency arrays required')
    if not np.all(np.isfinite(up)) or not np.all(np.isfinite(down)) or np.any(up<=0) or np.any(down<=0):
        raise ValueError('Positive finite calibrated carrier frequencies required')
    # Normalize each plasma coefficient to the lowest frequency on that
    # leg. This preserves 1/f^2 while avoiding ill-conditioned unit columns.
    return np.column_stack([np.ones(len(up)),(np.min(up)/up)**2,(np.min(down)/down)**2])

def nondispersive_weights(uplink_hz,downlink_hz,covariance=None):
    matrix=mixing_matrix(uplink_hz,downlink_hz);n=len(matrix)
    if n<3 or np.linalg.matrix_rank(matrix)<3:
        raise ValueError('Independent uplink/downlink plasma removal requires three independent link responses')
    noise=np.eye(n) if covariance is None else np.asarray(covariance,dtype=float)
    if noise.shape!=(n,n) or not np.all(np.isfinite(noise)) or not np.allclose(noise,noise.T,atol=1e-15,rtol=1e-12):
        raise ValueError('Finite symmetric link covariance required')
    try:np.linalg.cholesky(noise)
    except np.linalg.LinAlgError as exc:raise ValueError('Positive-definite covariance required') from exc
    whitened=np.linalg.solve(noise,matrix)
    normal=matrix.T@whitened
    if np.linalg.cond(normal)>1e12:raise ValueError('Carrier responses are too poorly conditioned')
    # Minimum-variance linear weights satisfying w A = (1,0,0).
    weights=np.linalg.solve(normal,[1.,0.,0.])@whitened.T
    if not np.allclose(weights@matrix,[1,0,0],atol=1e-10,rtol=0):
        raise ValueError('Plasma cancellation constraints failed')
    return weights,float(weights@noise@weights)

def nondispersive_residual(residuals,uplink_hz,downlink_hz,covariance=None):
    values=np.asarray(residuals,dtype=float)
    weights,variance=nondispersive_weights(uplink_hz,downlink_hz,covariance)
    if values.shape[-1:]!=(len(weights),) or not np.all(np.isfinite(values)):
        raise ValueError('Last residual dimension must contain all aligned links')
    return values@weights,variance
