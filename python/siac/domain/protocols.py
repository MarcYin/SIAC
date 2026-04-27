"""
Protocol definitions for SIAC pluggable modules.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol, runtime_checkable

if TYPE_CHECKING:
    from collections.abc import Sequence
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


@runtime_checkable
class SatellitePreprocessor(Protocol):
    """Protocol for satellite-specific data preprocessing."""

    @property
    def sensor_config(self) -> SensorConfig: ...

    def load_toa(self, input_path: str) -> xr.Dataset: ...

    def extract_geometry(self, input_path: str) -> GeometryAngles: ...

    def extract_cloud_mask(
        self,
        input_path: str,
        toa: xr.Dataset | None = None,
    ) -> xr.DataArray: ...

    def get_metadata(self, input_path: str) -> dict: ...


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


@runtime_checkable
class RTModelBackend(Protocol):
    """Protocol for RT backends."""

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
