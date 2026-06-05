#!/usr/bin/env python3
"""
TEP Evidence and Model-Comparison Utilities
==========================================

Version: TEP v0.9 (Jakarta)

Provides Bayesian evidence calculations and model-comparison statistics
used across the TEP corpus.

All TEP papers should import from this module to ensure consistent
statistical methodology.
"""

import numpy as np


def log_bayes_factor_from_bf(bf):
    """Convert a Bayes factor to log_10 Bayes factor."""
    if bf <= 0:
        return -np.inf
    return np.log10(bf)


def jeffreys_scale_interpretation(log10_bf):
    """
    Return Jeffreys scale verbal interpretation of log_10(BF).

    Parameters
    ----------
    log10_bf : float
        Log_10 of the Bayes factor.

    Returns
    -------
    str
        Verbal interpretation.
    """
    if log10_bf < 0:
        return "Negative (favours null / standard model)"
    elif log10_bf < 0.5:
        return "Not worth more than a bare mention"
    elif log10_bf < 1.0:
        return "Substantial"
    elif log10_bf < 1.5:
        return "Strong"
    elif log10_bf < 2.0:
        return "Very strong"
    else:
        return "Decisive"


def approximate_evidence_from_chi2(chi2, n_dof, n_params):
    """
    Approximate log-evidence from chi-squared fit.

    log Z ~ -0.5 * chi2 - 0.5 * n_params * log(N)

    where N is the effective number of data points.
    This is a Laplace approximation, not a full nested-sampling result.
    """
    if n_dof <= 0:
        raise ValueError("n_dof must be positive")
    if n_params < 0:
        raise ValueError("n_params must be non-negative")
    return -0.5 * chi2 - 0.5 * n_params * np.log(max(n_dof, 1))


def model_comparison_summary(evidence_null, evidence_alt, model_names=("Null", "TEP")):
    """
    Return a dict with Bayes factor, log10 BF, and Jeffreys interpretation.

    Parameters
    ----------
    evidence_null : float
        Log-evidence for the null / standard model.
    evidence_alt : float
        Log-evidence for the alternative (TEP) model.
    model_names : tuple
        (null_name, alt_name)

    Returns
    -------
    dict
        {'bf', 'log10_bf', 'strength', 'favoured_model'}
    """
    log_bf = evidence_alt - evidence_null
    bf = np.exp(log_bf)
    log10_bf = np.log10(bf)
    strength = jeffreys_scale_interpretation(log10_bf)
    favoured = model_names[1] if log_bf > 0 else model_names[0]

    return {
        "bf": float(bf),
        "log_bf": float(log_bf),
        "log10_bf": float(log10_bf),
        "strength": strength,
        "favoured_model": favoured,
    }


def angular_separation_deg(ra1, dec1, ra2, dec2):
    """Angular separation in degrees."""
    import numpy as np
    dra = np.radians(float(ra2) - float(ra1))
    ddec = np.radians(float(dec2) - float(dec1))
    a = (
        np.sin(ddec / 2.0) ** 2
        + np.cos(np.radians(float(dec1)))
        * np.cos(np.radians(float(dec2)))
        * np.sin(dra / 2.0) ** 2
    )
    return float(np.degrees(2.0 * np.arcsin(np.sqrt(min(1.0, a)))))


def p_value_from_sigma(sigma):
    """Two-tailed p-value from a standard-deviation threshold."""
    from scipy.stats import norm
    return 2.0 * (1.0 - norm.cdf(abs(sigma)))


def sigma_from_p_value(p):
    """Two-tailed sigma from a p-value."""
    from scipy.stats import norm
    if p <= 0 or p >= 1:
        return np.inf
    return float(norm.ppf(1.0 - p / 2.0))
