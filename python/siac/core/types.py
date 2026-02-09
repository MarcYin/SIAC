"""
Core data types for SIAC atmospheric correction.

This module defines the fundamental data structures used throughout SIAC,
including viewing geometry, atmospheric state, radiative transfer coefficients,
BRDF parameters, and sensor configurations.

All types are immutable (frozen dataclasses) to ensure data integrity and
enable caching/hashing where appropriate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Sequence

import numpy as np
import xarray as xr
from numpy.typing import NDArray


# =============================================================================
# Geometry Types
# =============================================================================


@dataclass(frozen=True)
class GeometryAngles:
    """
    View and sun geometry for atmospheric correction.

    All angles are stored in radians. The relative azimuth angle (raa)
    is computed as a property from the view and solar azimuths.

    Attributes:
        sza: Solar zenith angle in radians (0 = nadir, π/2 = horizon)
        saa: Solar azimuth angle in radians (0 = north, clockwise)
        vza: View zenith angle in radians (0 = nadir)
        vaa: View azimuth angle in radians (0 = north, clockwise)
    """

    sza: xr.DataArray
    saa: xr.DataArray
    vza: xr.DataArray
    vaa: xr.DataArray

    @property
    def raa(self) -> xr.DataArray:
        """Relative azimuth angle (view - solar) in radians."""
        return self.vaa - self.saa

    @property
    def cos_sza(self) -> xr.DataArray:
        """Cosine of solar zenith angle."""
        return np.cos(self.sza)

    @property
    def cos_vza(self) -> xr.DataArray:
        """Cosine of view zenith angle."""
        return np.cos(self.vza)

    @property
    def cos_raa(self) -> xr.DataArray:
        """Cosine of relative azimuth angle."""
        return np.cos(self.raa)

    def to_emulator_input(self) -> xr.Dataset:
        """
        Convert geometry to emulator input format.

        Returns:
            Dataset with cos(sza), cos(vza), cos(raa) variables.
        """
        return xr.Dataset(
            {
                "cos_sza": self.cos_sza,
                "cos_vza": self.cos_vza,
                "cos_raa": self.cos_raa,
            }
        )

    @classmethod
    def from_degrees(
        cls,
        sza_deg: xr.DataArray,
        saa_deg: xr.DataArray,
        vza_deg: xr.DataArray,
        vaa_deg: xr.DataArray,
    ) -> GeometryAngles:
        """Create GeometryAngles from angles in degrees."""
        deg_to_rad = np.pi / 180.0
        return cls(
            sza=sza_deg * deg_to_rad,
            saa=saa_deg * deg_to_rad,
            vza=vza_deg * deg_to_rad,
            vaa=vaa_deg * deg_to_rad,
        )


# =============================================================================
# Atmospheric State Types
# =============================================================================


@dataclass(frozen=True)
class AtmosphericState:
    """
    Atmospheric parameters for radiative transfer calculations.

    Attributes:
        aot: Aerosol Optical Thickness at 550nm (dimensionless)
        tcwv: Total Column Water Vapor (g/cm² or cm)
        tco3: Total Column Ozone (atm-cm or DU/1000)
        aot_unc: Uncertainty in AOT
        tcwv_unc: Uncertainty in TCWV
        tco3_unc: Uncertainty in TCO3
        elevation: Surface elevation in km
    """

    aot: xr.DataArray
    tcwv: xr.DataArray
    tco3: xr.DataArray
    aot_unc: xr.DataArray
    tcwv_unc: xr.DataArray
    tco3_unc: xr.DataArray
    elevation: xr.DataArray

    def to_emulator_input(self, geometry: GeometryAngles) -> xr.Dataset:
        """
        Combine atmospheric state with geometry for emulator input.

        The standard emulator input order is:
        [cos(sza), cos(vza), cos(raa), aot, tcwv, tco3, elevation]

        Args:
            geometry: Viewing geometry angles

        Returns:
            Dataset with all emulator input variables.
        """
        return xr.Dataset(
            {
                "cos_sza": geometry.cos_sza,
                "cos_vza": geometry.cos_vza,
                "cos_raa": geometry.cos_raa,
                "aot": self.aot,
                "tcwv": self.tcwv,
                "tco3": self.tco3,
                "elevation": self.elevation,
            }
        )

    def with_updated_aot_tcwv(
        self,
        aot: xr.DataArray,
        tcwv: xr.DataArray,
        aot_unc: xr.DataArray | None = None,
        tcwv_unc: xr.DataArray | None = None,
    ) -> AtmosphericState:
        """
        Create a new AtmosphericState with updated AOT and TCWV values.

        Used after aerosol retrieval to update the prior with solved values.
        """
        return AtmosphericState(
            aot=aot,
            tcwv=tcwv,
            tco3=self.tco3,
            aot_unc=aot_unc if aot_unc is not None else self.aot_unc,
            tcwv_unc=tcwv_unc if tcwv_unc is not None else self.tcwv_unc,
            tco3_unc=self.tco3_unc,
            elevation=self.elevation,
        )


# =============================================================================
# Radiative Transfer Types
# =============================================================================


@dataclass(frozen=True)
class RTCoefficients:
    """
    Radiative transfer coefficients from emulator or LUT.

    These coefficients are used in the atmospheric correction equation:
        y = xap * toa - xbp
        boa = y / (1 + xcp * y)

    Attributes:
        xap: Atmospheric path term (related to total transmittance)
        xbp: Background radiance / path reflectance term
        xcp: Adjacency / multiple scattering correction term
        d_xap: Jacobian of xap w.r.t. [aot, tcwv] (optional)
        d_xbp: Jacobian of xbp w.r.t. [aot, tcwv] (optional)
        d_xcp: Jacobian of xcp w.r.t. [aot, tcwv] (optional)
    """

    xap: xr.DataArray
    xbp: xr.DataArray
    xcp: xr.DataArray
    d_xap: xr.DataArray | None = None
    d_xbp: xr.DataArray | None = None
    d_xcp: xr.DataArray | None = None

    @property
    def has_jacobian(self) -> bool:
        """Check if Jacobians are available."""
        return self.d_xap is not None

    def apply_correction(self, toa: xr.DataArray) -> xr.DataArray:
        """
        Apply atmospheric correction to TOA reflectance.

        Args:
            toa: Top-of-atmosphere reflectance

        Returns:
            Bottom-of-atmosphere (surface) reflectance
        """
        y = self.xap * toa - self.xbp
        boa = y / (1.0 + self.xcp * y)
        return boa

    def compute_boa_jacobian(
        self, toa: xr.DataArray
    ) -> tuple[xr.DataArray, xr.DataArray]:
        """
        Compute Jacobian of BOA w.r.t. AOT and TCWV.

        Uses chain rule through the correction equation.

        Returns:
            Tuple of (d_boa/d_aot, d_boa/d_tcwv)
        """
        if not self.has_jacobian:
            raise ValueError("RTCoefficients does not have Jacobian information")

        y = self.xap * toa - self.xbp
        denom = 1.0 + self.xcp * y

        # dy/d_param = d_xap * toa - d_xbp
        # dboa/dy = 1/denom - xcp * y / denom^2 = (1 - xcp*y/denom) / denom
        #         = 1 / denom^2
        dboa_dy = 1.0 / (denom**2)

        # Chain rule for each parameter
        # d_xap, d_xbp, d_xcp have shape (..., 2) for [aot, tcwv]
        d_y_aot = self.d_xap.sel(param="aot") * toa - self.d_xbp.sel(param="aot")
        d_y_tcwv = self.d_xap.sel(param="tcwv") * toa - self.d_xbp.sel(param="tcwv")

        # Include xcp Jacobian contribution
        d_denom_aot = self.d_xcp.sel(param="aot") * y
        d_denom_tcwv = self.d_xcp.sel(param="tcwv") * y

        d_boa_aot = (d_y_aot * denom - y * d_denom_aot) / (denom**2)
        d_boa_tcwv = (d_y_tcwv * denom - y * d_denom_tcwv) / (denom**2)

        return d_boa_aot, d_boa_tcwv


# =============================================================================
# BRDF Types
# =============================================================================


@dataclass(frozen=True)
class BRDFKernelWeights:
    """
    Ross-Thick Li-Sparse BRDF kernel weights.

    The BRDF model is:
        ρ(θv, θs, φ) = f0 + f1 * K_vol + f2 * K_geo

    where K_vol is the Ross-Thick volumetric scattering kernel
    and K_geo is the Li-Sparse geometric scattering kernel.

    Attributes:
        f0: Isotropic kernel coefficient (wavelength dependent)
        f1: Volumetric scattering coefficient (Ross kernel)
        f2: Geometric scattering coefficient (Li kernel)
        f0_unc: Uncertainty in f0
        f1_unc: Uncertainty in f1
        f2_unc: Uncertainty in f2
    """

    f0: xr.DataArray  # Shape: (n_bands, y, x) or (y, x)
    f1: xr.DataArray
    f2: xr.DataArray
    f0_unc: xr.DataArray
    f1_unc: xr.DataArray
    f2_unc: xr.DataArray

    def compute_reflectance(
        self, k_vol: xr.DataArray, k_geo: xr.DataArray
    ) -> xr.DataArray:
        """
        Compute surface reflectance from BRDF kernels.

        Args:
            k_vol: Ross-Thick volumetric kernel values
            k_geo: Li-Sparse geometric kernel values

        Returns:
            Surface reflectance
        """
        return self.f0 + self.f1 * k_vol + self.f2 * k_geo

    def compute_reflectance_uncertainty(
        self, k_vol: xr.DataArray, k_geo: xr.DataArray
    ) -> xr.DataArray:
        """
        Compute uncertainty in surface reflectance.

        Assumes kernel uncertainties are negligible compared to
        coefficient uncertainties.
        """
        var = (
            self.f0_unc**2
            + (k_vol * self.f1_unc) ** 2
            + (k_geo * self.f2_unc) ** 2
        )
        return np.sqrt(var)


@dataclass(frozen=True)
class SurfacePrior:
    """
    Surface reflectance prior derived from BRDF model.

    Used as a constraint in the aerosol retrieval optimization.

    Attributes:
        boa: Prior surface reflectance (BOA)
        boa_unc: Uncertainty in surface reflectance
        kernels: BRDF kernel weights used to derive the prior
        mask: Valid pixel mask (True = valid)
    """

    boa: xr.DataArray
    boa_unc: xr.DataArray
    kernels: BRDFKernelWeights
    mask: xr.DataArray


# =============================================================================
# Sensor Configuration Types
# =============================================================================


@dataclass(frozen=True)
class SensorBand:
    """
    Specification for a single sensor band.

    Attributes:
        name: Band identifier (e.g., "B02", "B1")
        center_wavelength: Center wavelength in nanometers
        bandwidth: Full-width half-maximum in nanometers
        resolution: Native spatial resolution in meters
        band_index: Zero-based index in data arrays
    """

    name: str
    center_wavelength: float  # nm
    bandwidth: float  # nm
    resolution: float  # meters
    band_index: int

    @property
    def wavelength_um(self) -> float:
        """Center wavelength in micrometers."""
        return self.center_wavelength / 1000.0


@dataclass(frozen=True)
class SensorConfig:
    """
    Complete sensor configuration for a satellite instrument.

    Attributes:
        sensor_id: Sensor identifier (e.g., "MSI", "OLI")
        satellite_id: Satellite identifier (e.g., "S2A", "L8")
        bands: Tuple of band specifications
        default_ref_scale: Default reflectance scale factor
        default_ref_offset: Default reflectance offset
    """

    sensor_id: str
    satellite_id: str
    bands: tuple[SensorBand, ...]
    default_ref_scale: float = 1.0 / 10000.0
    default_ref_offset: float = 0.0

    def get_band(self, name: str) -> SensorBand:
        """Get band specification by name."""
        for band in self.bands:
            if band.name == name:
                return band
        raise KeyError(f"Band {name!r} not found in sensor {self.sensor_id}")

    def get_band_by_wavelength(
        self, wavelength_nm: float, tolerance_nm: float = 20.0
    ) -> SensorBand | None:
        """Find band closest to specified wavelength."""
        closest = None
        min_diff = float("inf")

        for band in self.bands:
            diff = abs(band.center_wavelength - wavelength_nm)
            if diff < min_diff and diff <= tolerance_nm:
                min_diff = diff
                closest = band

        return closest

    @property
    def band_names(self) -> tuple[str, ...]:
        """Tuple of all band names."""
        return tuple(b.name for b in self.bands)

    @property
    def band_wavelengths(self) -> tuple[float, ...]:
        """Tuple of all band center wavelengths."""
        return tuple(b.center_wavelength for b in self.bands)

    def select_bands_in_range(
        self, wl_min_nm: float, wl_max_nm: float
    ) -> list[SensorBand]:
        """Select bands whose center wavelength falls within [wl_min_nm, wl_max_nm]."""
        return [
            b for b in self.bands
            if wl_min_nm <= b.center_wavelength <= wl_max_nm
        ]

    def select_nearest_band(
        self, target_nm: float, tolerance_nm: float = 50.0
    ) -> SensorBand | None:
        """Find the band closest to target_nm within tolerance."""
        closest = None
        min_diff = float("inf")
        for band in self.bands:
            diff = abs(band.center_wavelength - target_nm)
            if diff < min_diff and diff <= tolerance_nm:
                min_diff = diff
                closest = band
        return closest

    @property
    def vis_bands(self) -> list[SensorBand]:
        """Visible bands (400-700 nm)."""
        return self.select_bands_in_range(400.0, 700.0)

    @property
    def nir_bands(self) -> list[SensorBand]:
        """Near-infrared bands (750-1000 nm)."""
        return self.select_bands_in_range(750.0, 1000.0)

    @property
    def swir_bands(self) -> list[SensorBand]:
        """Shortwave infrared bands (1000-2500 nm), excluding absorption windows."""
        return [
            b for b in self.bands
            if 1000.0 <= b.center_wavelength <= 2500.0
            and not (1350.0 <= b.center_wavelength <= 1420.0)  # cirrus
            and not (1800.0 <= b.center_wavelength <= 1950.0)  # WV absorption
        ]


# =============================================================================
# Predefined Sensor Configurations
# =============================================================================

SENTINEL2A_CONFIG = SensorConfig(
    sensor_id="MSI",
    satellite_id="S2A",
    bands=(
        SensorBand("B01", 443.0, 20.0, 60.0, 0),   # Coastal aerosol
        SensorBand("B02", 490.0, 65.0, 10.0, 1),   # Blue
        SensorBand("B03", 560.0, 35.0, 10.0, 2),   # Green
        SensorBand("B04", 665.0, 30.0, 10.0, 3),   # Red
        SensorBand("B05", 705.0, 15.0, 20.0, 4),   # Red edge 1
        SensorBand("B06", 740.0, 15.0, 20.0, 5),   # Red edge 2
        SensorBand("B07", 783.0, 20.0, 20.0, 6),   # Red edge 3
        SensorBand("B08", 842.0, 115.0, 10.0, 7),  # NIR
        SensorBand("B8A", 865.0, 20.0, 20.0, 8),   # NIR narrow
        SensorBand("B09", 945.0, 20.0, 60.0, 9),   # Water vapor
        SensorBand("B10", 1375.0, 30.0, 60.0, 10), # Cirrus
        SensorBand("B11", 1610.0, 90.0, 20.0, 11), # SWIR 1
        SensorBand("B12", 2190.0, 180.0, 20.0, 12), # SWIR 2
    ),
    default_ref_scale=1.0 / 10000.0,
    default_ref_offset=0.0,
)

SENTINEL2B_CONFIG = SensorConfig(
    sensor_id="MSI",
    satellite_id="S2B",
    bands=(
        SensorBand("B01", 442.0, 20.0, 60.0, 0),
        SensorBand("B02", 492.0, 65.0, 10.0, 1),
        SensorBand("B03", 559.0, 35.0, 10.0, 2),
        SensorBand("B04", 665.0, 30.0, 10.0, 3),
        SensorBand("B05", 704.0, 15.0, 20.0, 4),
        SensorBand("B06", 739.0, 15.0, 20.0, 5),
        SensorBand("B07", 780.0, 20.0, 20.0, 6),
        SensorBand("B08", 833.0, 115.0, 10.0, 7),
        SensorBand("B8A", 864.0, 20.0, 20.0, 8),
        SensorBand("B09", 943.0, 20.0, 60.0, 9),
        SensorBand("B10", 1377.0, 30.0, 60.0, 10),
        SensorBand("B11", 1610.0, 90.0, 20.0, 11),
        SensorBand("B12", 2186.0, 180.0, 20.0, 12),
    ),
    default_ref_scale=1.0 / 10000.0,
    default_ref_offset=0.0,
)

LANDSAT8_OLI_CONFIG = SensorConfig(
    sensor_id="OLI",
    satellite_id="L8",
    bands=(
        SensorBand("B1", 443.0, 16.0, 30.0, 0),    # Coastal aerosol
        SensorBand("B2", 482.0, 60.0, 30.0, 1),    # Blue
        SensorBand("B3", 561.5, 57.0, 30.0, 2),    # Green
        SensorBand("B4", 654.5, 37.0, 30.0, 3),    # Red
        SensorBand("B5", 865.0, 28.0, 30.0, 4),    # NIR
        SensorBand("B6", 1608.5, 85.0, 30.0, 5),   # SWIR 1
        SensorBand("B7", 2200.5, 187.0, 30.0, 6),  # SWIR 2
    ),
    default_ref_scale=2.75e-5,
    default_ref_offset=-0.2,
)

# Sensor registry for lookup
SENSOR_CONFIGS: dict[tuple[str, str], SensorConfig] = {
    ("MSI", "S2A"): SENTINEL2A_CONFIG,
    ("MSI", "S2B"): SENTINEL2B_CONFIG,
    ("OLI", "L8"): LANDSAT8_OLI_CONFIG,
}


# =============================================================================
# Pipeline Contract Types (Module Output Contracts)
# =============================================================================


@dataclass(frozen=True)
class ObservationBundle:
    """Complete observation data from satellite preprocessing.

    Output contract of M1 (Satellite Preprocessor).
    Contains everything the pipeline needs from the satellite input.
    """

    toa: xr.Dataset  # TOA reflectance, one var per band
    geometry: GeometryAngles  # SZA/SAA/VZA/VAA in radians
    cloud_mask: xr.DataArray  # bool, True = cloudy/invalid
    sensor_config: SensorConfig  # band definitions + scale/offset
    metadata: dict[str, Any]  # must include 'observation_time': datetime
    crs: str  # e.g. "EPSG:32632"
    bounds: tuple[float, float, float, float]  # (xmin, ymin, xmax, ymax)


@dataclass(frozen=True)
class SolverInputBundle:
    """All inputs to the aerosol solver, resampled to solver grids.

    Output contract of M4 (Grid Assembler). Everything in this bundle
    is spatially aligned and ready for the solver to consume directly.
    """

    # Observation (resampled to aux resolution)
    toa: xr.DataArray  # (bands, y, x) at aux resolution
    geometry: GeometryAngles  # at aux resolution
    cloud_mask: xr.DataArray  # (y, x) at aux resolution
    sensor_config: SensorConfig
    bands: list[SensorBand]  # solver bands (wavelength-selected)

    # Atmospheric prior (resampled to aux resolution)
    atmo_prior: AtmosphericState  # all fields at aux resolution

    # Surface prior (resampled to aux resolution)
    surface_prior: SurfacePrior  # at aux resolution

    # RT backend (not resampled — it's a model, not raster data)
    rt_model: Any  # RTModelBackend (Any to avoid circular import)

    # Grid metadata
    aux_resolution_m: float  # e.g. 500.0
    aerosol_resolution_m: float  # e.g. 1000.0


@dataclass(frozen=True)
class SolvedAtmosphere:
    """Solver output: retrieved atmospheric parameters + diagnostics.

    Output contract of M5 (Aerosol Solver). Contains the solved AOT/TCWV
    fields at the aerosol retrieval resolution, plus the full
    AtmosphericState (with solved values merged in) for use by the
    corrector.
    """

    atmo_state: AtmosphericState  # full state with solved AOT/TCWV
    aot: xr.DataArray  # solved AOT at aerosol resolution
    tcwv: xr.DataArray  # solved TCWV at aerosol resolution
    aot_unc: xr.DataArray  # posterior AOT uncertainty
    tcwv_unc: xr.DataArray  # posterior TCWV uncertainty

    # Diagnostics
    cost_final: float  # final cost function value
    n_iterations: int  # total optimizer iterations
    converged: bool  # did the solver converge?


@dataclass(frozen=True)
class CorrectionResult:
    """Final output of atmospheric correction.

    Output contract of M6 (Atmospheric Corrector).
    """

    boa: xr.Dataset  # BOA reflectance, one var per band, at native resolution
    boa_unc: xr.Dataset | None  # per-band uncertainty (optional)
    aot: xr.DataArray  # solved AOT map
    tcwv: xr.DataArray  # solved TCWV map
    cloud_mask: xr.DataArray  # final cloud mask
    metadata: dict[str, Any]  # processing metadata (timings, versions, etc.)


def get_sensor_config(sensor_id: str, satellite_id: str) -> SensorConfig:
    """
    Get sensor configuration by ID.

    Args:
        sensor_id: Sensor identifier (e.g., "MSI", "OLI")
        satellite_id: Satellite identifier (e.g., "S2A", "L8")

    Returns:
        SensorConfig for the specified sensor

    Raises:
        KeyError: If sensor is not found in registry
    """
    key = (sensor_id, satellite_id)
    if key not in SENSOR_CONFIGS:
        raise KeyError(
            f"Unknown sensor {sensor_id}/{satellite_id}. "
            f"Available: {list(SENSOR_CONFIGS.keys())}"
        )
    return SENSOR_CONFIGS[key]
