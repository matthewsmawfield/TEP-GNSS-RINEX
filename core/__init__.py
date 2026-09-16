"""TEP Core — canonical Python package for the Temporal Equivalence Principle.

Version: TEP v0.10 (Jakarta)

This package provides the shared physics layer used by all TEP papers:
  - constants: physical and phenomenological parameters
  - conformal_scaling: A(phi), Temporal Shear, effective coupling
  - cosmology: LCDM + TEP distance-redshift relations
  - scalar_field: lab-scale scalar-field solver
  - screening: environment-dependent screening functions
  - evidence: Bayesian evidence and model-comparison utilities
"""

from . import constants
from . import conformal_scaling
from . import cosmology
from . import scalar_field
from . import screening
from . import evidence

__all__ = [
    "constants",
    "conformal_scaling",
    "cosmology",
    "scalar_field",
    "screening",
    "evidence",
]

# Re-export commonly used symbols at package level for convenience
from .constants import (
    VERSION,
    VERSION_CODENAME,
    VERSION_STRING,
    G_NEWTON,
    C_LIGHT,
    M_PLANCK,
    M_SUN,
    MPC_TO_M,
    BETA_A,
    RHO_T,
    RHO_C,
    ILLUSTRATIVE_BETA_A,
    KG_M3_TO_G_CM3,
    G_CM3_TO_KG_M3,
    LAB_COHERENCE_LENGTH_M,
    M_REF,
    ALPHA_LOG,
    BETA_GEOM,
    SCREENING_LENGTH_KM,
    LAMBDA_T_MGEX_KM,
    LAMBDA_T_MGEX_ERR_KM,
    LAMBDA_T_MGEX_R2,
    GNSS_LAMBDA_T_LONGSPAN_CODE_KM,
    GNSS_LAMBDA_T_LONGSPAN_CODE_ERR_KM,
    GNSS_LAMBDA_T_EXPONENTIAL_BY_CENTER,
)
from .screening import (
    screening_factor,
    universal_screening_function,
    coupling_screening_factor,
    beta_screened,
)
from .conformal_scaling import (
    conformal_factor,
    conformal_factor_small,
    temporal_shear_from_scalar_field,
    effective_g,
    g_eff_variance,
    minimum_steepness_for_retention,
    minimum_steepness_for_suppression,
    screening_diagnostics,
)
from .scalar_field import (
    solve_scalar_field_cylinder,
    solve_scalar_field_uniform_density,
    solve_scalar_field_layered,
    solve_scalar_field_layered_weighted,
    scalar_field_logarithmic,
    scalar_field_difference,
    compute_temporal_shear_from_mass_gradient,
)
