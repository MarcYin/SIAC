"""
Aerosol retrieval solver module.

Implements multi-grid L-BFGS-B optimization for AOT and TCWV retrieval.

The solver minimizes a cost function with three components:
- Observation cost: (boa_model - boa_prior)^2 / sigma^2
- Prior cost: (param - param_prior)^2 / sigma^2
- Smoothness cost: gamma^2 * ||grad(param)||^2

The multi-grid approach progressively refines from coarse to fine resolution.
"""

from siac.solver.cost import (
    CostFunction,
    CostFunctionConfig,
    compute_laplacian_eigenvalues,
    apply_smoothness_filter,
    create_sparse_laplacian,
)
from siac.solver.multigrid import (
    MultiGridSolver,
    MultiGridConfig,
    SolverResult,
    solve_atmospheric_parameters,
)

__all__ = [
    # Cost function
    "CostFunction",
    "CostFunctionConfig",
    "compute_laplacian_eigenvalues",
    "apply_smoothness_filter",
    "create_sparse_laplacian",
    # Multi-grid solver
    "MultiGridSolver",
    "MultiGridConfig",
    "SolverResult",
    "solve_atmospheric_parameters",
]
