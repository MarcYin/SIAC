"""
Protocol definitions for SIAC pluggable modules.

This module defines the interfaces (Protocols) that allow SIAC components
to be swapped at runtime. Each protocol specifies the required methods
and their signatures for a particular type of functionality.

Using Python's structural subtyping (typing.Protocol), implementations
don't need to explicitly inherit from these protocols - they just need
to implement the required methods with matching signatures.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol, runtime_checkable

if TYPE_CHECKING:
    from collections.abc import Sequence
    from datetime import datetime

    import xarray as xr

    from siac.core.types import (
        AtmosphericState,
        BRDFKernelWeights,
        GeometryAngles,
        RTCoefficients,
        SensorBand,
        SensorConfig,
        SurfacePrior,
    )


# =============================================================================
# Satellite Preprocessing Protocol
# =============================================================================


@runtime_checkable
class SatellitePreprocessor(Protocol):
    """
    Protocol for satellite-specific data preprocessing.

    Implementations handle the specifics of reading data from different
    satellite formats (e.g., Sentinel-2 SAFE, Landsat MTL) and extracting
    the information needed for atmospheric correction.

    Example implementations:
        - Sentinel2Preprocessor: Reads S2 SAFE format
        - Landsat8Preprocessor: Reads L8/L9 MTL format
        - GenericPreprocessor: Reads generic COG/NetCDF with metadata
    """

    @property
    def sensor_config(self) -> SensorConfig:
        """
        Return the sensor configuration for this preprocessor.

        Returns:
            SensorConfig with band specifications and calibration constants.
        """
        ...

    def load_toa(self, input_path: str) -> xr.Dataset:
        """
        Load top-of-atmosphere reflectance data.

        Args:
            input_path: Path to satellite data directory or file

        Returns:
            xr.Dataset with TOA reflectance for each band.
            Variables should be named by band (e.g., "B02", "B1").
            Dataset should include CRS and transform in attrs.
        """
        ...

    def extract_geometry(self, input_path: str) -> GeometryAngles:
        """
        Extract view and sun geometry from satellite data.

        Args:
            input_path: Path to satellite data directory or file

        Returns:
            GeometryAngles with sza, saa, vza, vaa in radians.
        """
        ...

    def extract_cloud_mask(
        self,
        input_path: str,
        toa: xr.Dataset | None = None,
    ) -> xr.DataArray:
        """
        Generate cloud mask from satellite data.

        Args:
            input_path: Path to satellite data directory or file
            toa: Optional preloaded TOA dataset

        Returns:
            xr.DataArray with boolean mask (True = cloudy/invalid).
        """
        ...

    def get_metadata(self, input_path: str) -> dict:
        """
        Extract observation metadata from satellite data.

        Args:
            input_path: Path to satellite data directory or file

        Returns:
            Dictionary with at minimum:
                - "observation_time": datetime
                - "tile_id": str (optional)
                - "processing_baseline": str (optional)
        """
        ...


# =============================================================================
# Atmospheric Prior Protocol
# =============================================================================


@runtime_checkable
class AtmosphericPriorProvider(Protocol):
    """
    Protocol for atmospheric state prior data providers.

    Implementations fetch atmospheric parameters (AOT, TCWV, TCO3) from
    various sources to use as priors in the aerosol retrieval.

    Example implementations:
        - CAMSProvider: ECMWF CAMS near real-time data
        - MERRA2Provider: NASA MERRA-2 reanalysis
        - ERA5Provider: ECMWF ERA5 reanalysis
        - UserProvider: User-supplied prior arrays
    """

    def get_prior(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        resolution: float,
    ) -> AtmosphericState:
        """
        Retrieve atmospheric priors for a given region and time.

        Args:
            bounds: Bounding box (xmin, ymin, xmax, ymax) in target CRS
            crs: Coordinate reference system (e.g., "EPSG:32632")
            obs_time: Observation datetime for temporal interpolation
            resolution: Target spatial resolution in meters

        Returns:
            AtmosphericState with AOT, TCWV, TCO3 and their uncertainties.
        """
        ...

    @property
    def source_name(self) -> str:
        """
        Name of the prior data source.

        Returns:
            String identifier (e.g., "CAMS", "MERRA-2", "ERA5")
        """
        ...


# =============================================================================
# BRDF Product Protocol
# =============================================================================


@runtime_checkable
class BRDFProductProvider(Protocol):
    """
    Protocol for BRDF product data providers.

    Implementations fetch BRDF kernel parameters from various sources
    to compute surface reflectance priors for atmospheric inversion.

    Example implementations:
        - MCD43Provider: MODIS Collection 6 MCD43A1
        - VNP43Provider: VIIRS VNP43A1
        - MCD19Provider: MODIS MCD19A3
        - GEEProvider: Google Earth Engine access
        - ZarrProvider: Pre-cached Zarr store
    """

    def get_brdf_parameters(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        target_resolution: float,
        bands: Sequence[int],
        temporal_window: int = 16,
    ) -> BRDFKernelWeights:
        """
        Retrieve BRDF kernel parameters for a given region and time.

        Args:
            bounds: Bounding box (xmin, ymin, xmax, ymax) in target CRS
            crs: Coordinate reference system
            obs_time: Observation datetime (center of temporal window)
            target_resolution: Target spatial resolution in meters
            bands: MODIS band numbers to retrieve (1-7)
            temporal_window: Days before/after obs_time to include

        Returns:
            BRDFKernelWeights with f0, f1, f2 and uncertainties.
            Shape should be (n_bands, y, x).
        """
        ...

    @property
    def source_name(self) -> str:
        """
        Name of the BRDF product source.

        Returns:
            String identifier (e.g., "MCD43", "VNP43", "MCD19")
        """
        ...


# =============================================================================
# Surface Prior Protocol
# =============================================================================


@runtime_checkable
class SurfacePriorDeriver(Protocol):
    """
    Protocol for deriving surface reflectance prior from BRDF parameters.

    Implementations compute the expected surface reflectance given BRDF
    kernel weights and observation geometry, optionally accounting for
    scale differences via PSF modeling.

    Example implementations:
        - KernelModelDeriver: Standard Ross-Thick Li-Sparse model
        - NeuralDeriver: Neural network-based BRDF interpolation
        - DirectDeriver: Direct BRDF sampling without PSF
    """

    def compute_surface_prior(
        self,
        brdf_weights: BRDFKernelWeights,
        geometry: GeometryAngles,
        psf_params: tuple[float, float] | None = None,
    ) -> SurfacePrior:
        """
        Compute surface reflectance prior from BRDF parameters.

        Args:
            brdf_weights: BRDF kernel coefficients (f0, f1, f2)
            geometry: Observation geometry (sza, vza, raa)
            psf_params: Optional (sigma_x, sigma_y) for PSF convolution
                        to match coarse BRDF resolution to fine image resolution

        Returns:
            SurfacePrior with BOA reflectance, uncertainty, and mask.
        """
        ...


# =============================================================================
# Radiative Transfer Model Protocol
# =============================================================================


@runtime_checkable
class RTModelBackend(Protocol):
    """
    Protocol for radiative transfer model backends.

    Implementations compute the atmospheric correction coefficients
    (xap, xbp, xcp) given viewing geometry and atmospheric state.

    Example implementations:
        - EmulatorBackend: Pre-trained neural network emulators (fast)
        - LUTBackend: Look-up table interpolation (medium)
        - Py6SBackend: Direct 6S simulation via Py6S (slow)
    """

    def compute_coefficients(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        band: SensorBand,
        compute_jacobian: bool = False,
    ) -> RTCoefficients:
        """
        Compute radiative transfer coefficients for atmospheric correction.

        The correction equation is:
            y = xap * toa - xbp
            boa = y / (1 + xcp * y)

        Args:
            geometry: Viewing geometry (sza, vza, raa)
            atmo_state: Atmospheric state (aot, tcwv, tco3, elevation)
            band: Sensor band specification
            compute_jacobian: Whether to compute d(coeff)/d(aot,tcwv)

        Returns:
            RTCoefficients with xap, xbp, xcp and optional Jacobians.
        """
        ...

    def supports_jacobian(self) -> bool:
        """
        Check if this backend supports analytical Jacobian computation.

        Emulator backends typically support analytical Jacobians via
        backpropagation. LUT and direct simulation backends may only
        support numerical differentiation.

        Returns:
            True if analytical Jacobians are available.
        """
        ...

    @property
    def backend_name(self) -> str:
        """
        Name of the RT backend.

        Returns:
            String identifier (e.g., "emulator", "lut", "py6s")
        """
        ...

    def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool:
        """
        Check if this backend has data/models for the specified sensor.

        Emulator backends may only support S2/L8, while LUT and direct
        simulation backends can support any sensor.

        Args:
            sensor_id: Sensor identifier (e.g., "MSI", "OLI")
            satellite_id: Satellite identifier (e.g., "S2A", "L8")

        Returns:
            True if sensor is supported.
        """
        ...


# =============================================================================
# Aerosol Solver Protocol
# =============================================================================


@runtime_checkable
class AerosolSolver(Protocol):
    """
    Protocol for aerosol retrieval solvers.

    Implementations solve for AOT and TCWV by minimizing the difference
    between observed and modeled TOA reflectance, subject to prior
    constraints and smoothness regularization.

    Example implementations:
        - MultiGridSolver: Multi-grid L-BFGS-B optimization
        - PixelwiseSolver: Independent pixel-by-pixel optimization
        - BayesianSolver: Full Bayesian inference
    """

    def solve(
        self,
        toa: xr.DataArray,
        surface_prior: SurfacePrior,
        geometry: GeometryAngles,
        atmo_prior: AtmosphericState,
        rt_model: RTModelBackend,
        cloud_mask: xr.DataArray,
        solver_config: dict,
    ) -> AtmosphericState:
        """
        Solve for atmospheric parameters (AOT, TCWV).

        The optimization minimizes:
            J = J_obs + J_prior + J_smooth

        where:
            J_obs = Σ (boa_model - boa_prior)² / σ²_boa
            J_prior = (aot - aot_prior)²/σ²_aot + (tcwv - tcwv_prior)²/σ²_tcwv
            J_smooth = γ_aot * ||∇aot||² + γ_tcwv * ||∇tcwv||²

        Args:
            toa: Top-of-atmosphere reflectance
            surface_prior: Surface reflectance prior from BRDF
            geometry: Viewing geometry
            atmo_prior: Prior atmospheric state from CAMS/MERRA-2
            rt_model: Radiative transfer model backend
            cloud_mask: Boolean mask (True = cloudy)
            solver_config: Configuration dictionary with:
                - aot_gamma: Smoothness weight for AOT
                - tcwv_gamma: Smoothness weight for TCWV
                - aot_bounds: (min, max) bounds for AOT
                - tcwv_bounds: (min, max) bounds for TCWV
                - max_iterations: Maximum optimizer iterations
                - aerosol_resolution: Grid resolution for solving (m)

        Returns:
            AtmosphericState with solved AOT, TCWV and uncertainties.
        """
        ...


# =============================================================================
# Output Writer Protocol
# =============================================================================


@runtime_checkable
class OutputWriter(Protocol):
    """
    Protocol for writing output products.

    Implementations handle the specifics of writing corrected imagery
    and auxiliary products in various formats.

    Example implementations:
        - GeoTiffWriter: Standard GeoTIFF output
        - COGWriter: Cloud-optimized GeoTIFF
        - ZarrWriter: Zarr store for analysis-ready data
        - NetCDFWriter: NetCDF-CF compliant output
    """

    def write_boa(
        self,
        boa: xr.Dataset,
        output_path: str,
        crs: str,
        transform: tuple,
    ) -> None:
        """
        Write BOA reflectance products.

        Args:
            boa: Dataset with BOA reflectance per band
            output_path: Output file or directory path
            crs: Coordinate reference system
            transform: Affine transform tuple
        """
        ...

    def write_auxiliary(
        self,
        aot: xr.DataArray,
        tcwv: xr.DataArray,
        output_path: str,
        crs: str,
        transform: tuple,
    ) -> None:
        """
        Write auxiliary products (AOT, TCWV maps).

        Args:
            aot: Aerosol optical thickness
            tcwv: Total column water vapor
            output_path: Output file or directory path
            crs: Coordinate reference system
            transform: Affine transform tuple
        """
        ...
