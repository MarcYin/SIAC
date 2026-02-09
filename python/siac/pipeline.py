"""
Pipeline orchestration for SIAC atmospheric correction.

This module contains the ``run_pipeline()`` function that wires modules
M1–M6 together.  Each module is passed in as a plain callable — the
pipeline never instantiates providers itself.

Type aliases for each callable signature are defined here as well so that
type checkers can verify user-provided overrides.
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Protocol

import xarray as xr

from siac.core.aoi import AOI
from siac.core.types import (
    AtmosphericState,
    CorrectionResult,
    GeometryAngles,
    ObservationBundle,
    SensorBand,
    SensorConfig,
    SolvedAtmosphere,
    SolverInputBundle,
    SurfacePrior,
)
from siac.core.validation import (
    _validate_atmospheric_state,
    _validate_observation_bundle,
    _validate_surface_prior,
)

logger = logging.getLogger(__name__)

# ── Type aliases for module callables ──────────────────────────────────

# M1 preprocessor
PreprocessorFn = Callable[[Path, Any], ObservationBundle]  # (path, aoi|None) -> ObservationBundle

# M2 atmospheric prior provider
AtmoPriorFn = Callable[
    [tuple[float, float, float, float], str, datetime, float],
    AtmosphericState,
]

# M3 surface prior provider
SurfacePriorFn = Callable[
    [tuple[float, float, float, float], str, datetime, SensorConfig, GeometryAngles, float],
    SurfacePrior,
]

# M4 grid assembler
GridAssemblerFn = Callable[
    [ObservationBundle, AtmosphericState, SurfacePrior, Any, float, float],
    SolverInputBundle,
]

# M5 aerosol solver
SolverFn = Callable[[SolverInputBundle, Any], SolvedAtmosphere]

# M6 atmospheric corrector
CorrectorFn = Callable[[ObservationBundle, SolvedAtmosphere, Any], CorrectionResult]


# ── Pipeline orchestrator ──────────────────────────────────────────────

def run_pipeline(
    input_path: Path,
    aoi: AOI | None,
    config: Any,
    *,
    preprocessor: PreprocessorFn,
    atmo_provider: AtmoPriorFn,
    surface_prior_provider: SurfacePriorFn,
    grid_assembler: GridAssemblerFn,
    solver: SolverFn,
    corrector: CorrectorFn,
    rt_model: Any,
    max_workers: int = 4,
) -> CorrectionResult:
    """Orchestrate module execution with concurrent data sourcing.

    This is a plain function — no class, no state.  Each module callable
    is passed as an argument (either a bound method or a plain function).

    Phases:
        1. M1 (preprocess) — must complete first to get bounds/crs/time.
        2. M2 + M3 run concurrently (independent I/O).
        3. M4 → M5 → M6 run sequentially.
    """
    # Phase 1: Preprocess
    logger.info("M1: Preprocessing satellite data…")
    obs = preprocessor(input_path, aoi)
    _validate_observation_bundle(obs)

    bounds = obs.bounds
    crs = obs.crs
    obs_time = obs.metadata["observation_time"]
    resolution = getattr(config, "solver", None)
    if resolution is not None:
        resolution = getattr(resolution, "aerosol_resolution", 1000.0)
    else:
        resolution = 1000.0

    # Phase 2: Concurrent data sourcing
    logger.info("M2+M3: Fetching atmospheric & surface priors…")
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        f_m2 = executor.submit(atmo_provider, bounds, crs, obs_time, resolution)
        f_m3 = executor.submit(
            surface_prior_provider,
            bounds,
            crs,
            obs_time,
            obs.sensor_config,
            obs.geometry,
            resolution,
        )
        atmo = f_m2.result()
        surface = f_m3.result()

    _validate_atmospheric_state(atmo)
    _validate_surface_prior(surface)

    # Phase 3: Sequential processing
    logger.info("M4: Assembling solver grids…")
    solver_inputs = grid_assembler(obs, atmo, surface, rt_model)

    logger.info("M5: Solving for aerosol parameters…")
    solved = solver(solver_inputs, config)

    logger.info("M6: Applying atmospheric correction…")
    result = corrector(obs, solved, rt_model)

    logger.info("Pipeline complete.")
    return result
