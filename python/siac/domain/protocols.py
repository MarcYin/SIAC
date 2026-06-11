"""
Protocol definitions for SIAC pluggable modules.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Protocol, runtime_checkable

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence
    from datetime import datetime
    from pathlib import Path

    import xarray as xr

    from siac.algorithms.surface.brdf_monthly_composite import MonthlyCompositeCollection
    from siac.domain.sensors import SensorBand, SensorConfig
    from siac.runtime.models import (
        AtmosphericState,
        BRDFKernelWeights,
        CorrectionResult,
        GeometryAngles,
        ObservationBundle,
        RTCoefficients,
        SolvedAtmosphere,
        SolverInputBundle,
        SurfacePrior,
    )


# NOTE: runtime_checkable here matches attribute presence only,
# not property semantics — see REVIEW.md §2.3
@runtime_checkable
class SatellitePreprocessor(Protocol):
    """Protocol for satellite-specific data preprocessing."""

    @property
    def sensor_config(self) -> SensorConfig: ...

    def load_toa(self, input_path: Path | str) -> xr.Dataset: ...

    def extract_geometry(self, input_path: Path | str) -> GeometryAngles: ...

    def extract_cloud_mask(
        self,
        input_path: Path | str,
        toa: xr.Dataset | None = None,
    ) -> xr.DataArray: ...

    def get_metadata(self, input_path: Path | str) -> dict[str, Any]: ...


# NOTE: runtime_checkable here matches attribute presence only,
# not property semantics — see REVIEW.md §2.3
@runtime_checkable
class AtmosphericPriorProvider(Protocol):
    """Protocol for atmospheric state prior data providers."""

    def get_prior(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        resolution: float,
    ) -> AtmosphericState: ...

    @property
    def source_name(self) -> str: ...


# NOTE: runtime_checkable here matches attribute presence only,
# not property semantics — see REVIEW.md §2.3
@runtime_checkable
class BRDFProductProvider(Protocol):
    """Protocol for BRDF product data providers."""

    def get_brdf_parameters(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        target_resolution: float,
        bands: Sequence[int],
        temporal_window: int = 16,
    ) -> BRDFKernelWeights: ...

    @property
    def source_name(self) -> str: ...


@runtime_checkable
class SurfacePriorDeriver(Protocol):
    """Protocol for deriving surface reflectance priors."""

    def compute_surface_prior(
        self,
        brdf_weights: BRDFKernelWeights,
        geometry: GeometryAngles,
        psf_params: tuple[float, float] | None = None,
    ) -> SurfacePrior: ...


# NOTE: runtime_checkable here matches attribute presence only,
# not property semantics — see REVIEW.md §2.3
@runtime_checkable
class MonthlyCompositeProvider(Protocol):
    """Protocol for providers returning prepared monthly composite inputs."""

    @property
    def source_name(self) -> str: ...

    @property
    def source_bands(self) -> Sequence[SensorBand]: ...

    def get_monthly_composites(
        self,
        observation: ObservationBundle,
        resolution: float,
    ) -> MonthlyCompositeCollection: ...


# NOTE: runtime_checkable here matches attribute presence only,
# not property semantics — see REVIEW.md §2.3
@runtime_checkable
class RTModelBackend(Protocol):
    """Protocol for RT backends.

    This is the REQUIRED contract every backend implements. Backends may
    additionally provide the optional capabilities below, discovered via the
    ``rt_*`` accessor functions in this module (a Protocol cannot express
    optional methods, and scattering ``getattr`` probes across the solver and
    pipeline hid which hooks exist):

    - ``compute_coefficients_multi(geometry, atmo_state, bands, ...)`` —
      batch per-band evaluation that amortises shared scene work.
    - ``build_joint_grid_search_lut(geometry=..., atmo_state=..., aot_axis=...,
      tcwv_axis=..., bands=...)`` — one LUT spanning the solver's block-grid
      search range, served per candidate by interpolation (6S, libRadtran).
    - ``preload_joint_grid_search_lut(...)`` — background warm-up of the above.
    - ``preload_scene_subset(geometry, atmo_state, bands)`` — background
      materialisation of the remote-LUT scene subset.
    - ``set_observation_time(datetime | None)`` — observation-time hint for
      backends whose physics depend on it.
    """

    def compute_coefficients(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        band: SensorBand,
        compute_jacobian: bool = False,
    ) -> RTCoefficients: ...

    def supports_jacobian(self) -> bool: ...

    @property
    def backend_name(self) -> str: ...

    def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool: ...


def rt_optional_capability(rt_model: object, name: str) -> Callable[..., Any] | None:
    """Return the named optional RT-backend hook, or None when not provided.

    The single discovery point for the optional capabilities documented on
    :class:`RTModelBackend` — call sites use this instead of ad-hoc
    ``getattr`` so the set of recognised hooks stays auditable here.
    """
    fn = getattr(rt_model, name, None)
    return fn if callable(fn) else None


def rt_supports_jacobian(rt_model: object) -> bool:
    """Whether the backend can provide per-pixel RT Jacobians (False on error)."""
    fn = rt_optional_capability(rt_model, "supports_jacobian")
    if fn is None:
        return False
    try:
        return bool(fn())
    except Exception:
        return False


@runtime_checkable
class AerosolSolver(Protocol):
    """Protocol for aerosol solvers."""

    def solve(
        self,
        inputs: SolverInputBundle,
    ) -> SolvedAtmosphere: ...


@runtime_checkable
class OutputWriter(Protocol):
    """Protocol for writing output products."""

    def write(
        self,
        result: CorrectionResult,
        output_dir: str | Path,
    ) -> dict[str, Path]: ...
