"""
Aerosol retrieval solver module.

Implements multi-grid L-BFGS-B optimization for AOT and TCWV retrieval.

The solver minimizes a cost function with three components:
- Observation cost: (boa_model - boa_prior)^2 / sigma^2
- Prior cost: (param - param_prior)^2 / sigma^2
- Smoothness cost: gamma^2 * ||grad(param)||^2

The multi-grid approach progressively refines from coarse to fine resolution.
"""

from siac.algorithms.solver.cost import (
    CostFunction as CostFunction,
)
from siac.algorithms.solver.cost import (
    CostFunctionConfig as CostFunctionConfig,
)
from siac.algorithms.solver.cost import (
    apply_smoothness_filter as apply_smoothness_filter,
)
from siac.algorithms.solver.cost import (
    compute_laplacian_eigenvalues as compute_laplacian_eigenvalues,
)
from siac.algorithms.solver.cost import (
    create_sparse_laplacian as create_sparse_laplacian,
)
from siac.algorithms.solver.multigrid import (
    MultiGridConfig as MultiGridConfig,
)
from siac.algorithms.solver.multigrid import (
    MultiGridSolver as MultiGridSolver,
)
from siac.algorithms.solver.multigrid import (
    SolverResult as SolverResult,
)
from siac.algorithms.solver.multigrid import (
    SolverStageConfig as SolverStageConfig,
)
from siac.algorithms.solver.multigrid import (
    StagedMultiGridSolver as StagedMultiGridSolver,
)
from siac.algorithms.solver.multigrid import (
    build_solver_valid_mask as build_solver_valid_mask,
)
from siac.algorithms.solver.multigrid import (
    solve_atmospheric_parameters as solve_atmospheric_parameters,
)
