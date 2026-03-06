"""
Configuration system for SIAC using Pydantic.

This module provides a hierarchical configuration system that allows users
to specify processing options via YAML files, environment variables, or
programmatically. All configuration is validated at load time.

Example usage:
    >>> from siac.core.config import SIACConfig
    >>>
    >>> # Load from YAML file
    >>> config = SIACConfig.from_yaml("siac_config.yaml")
    >>>
    >>> # Or create programmatically
    >>> config = SIACConfig(
    ...     sensor="s2",
    ...     rt_model={"backend": "lut", "lut_path": "/path/to/lut.zarr"},
    ... )
    >>>
    >>> # Save to YAML
    >>> config.to_yaml("output_config.yaml")
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal
from urllib.parse import urlparse

import yaml
from pydantic import BaseModel, Field, field_validator, model_validator
from pydantic_settings import BaseSettings

DEFAULT_LUT_URL = (
    "https://gws-access.jasmin.ac.uk/public/nceo_isp/"
    "libradtran_continental_average_lut_1nm.zarr.zip"
)


# =============================================================================
# Sub-configuration Models
# =============================================================================


class AtmoPriorConfig(BaseModel):
    """Configuration for atmospheric prior data provider."""

    provider: Literal["cams", "merra2", "mcd19", "era5", "user"] = Field(
        default="cams",
        description=(
            "Atmospheric data source: 'cams' (ECMWF), 'merra2' (NASA), "
            "'mcd19' (Earthaccess AOD-focused), 'era5', or 'user'"
        ),
    )
    data_path: str | Path | None = Field(
        default=None,
        description="Path or HTTP URL to atmospheric data directory or files",
    )
    cache_dir: Path | None = Field(
        default=None,
        description="Directory for caching downloaded atmospheric data",
    )
    download_missing: bool = Field(
        default=True,
        description="Auto-download missing atmospheric prior files when the provider supports it.",
    )
    temporal_interpolation: Literal["nearest", "linear"] = Field(
        default="nearest",
        description="Temporal interpolation method for 3-hourly data",
    )

    # User-provided prior arrays (when provider="user")
    user_aot: Path | None = Field(
        default=None,
        description="Path to user-provided AOT array (GeoTIFF)",
    )
    user_tcwv: Path | None = Field(
        default=None,
        description="Path to user-provided TCWV array (GeoTIFF)",
    )
    user_tco3: Path | None = Field(
        default=None,
        description="Path to user-provided TCO3 array (GeoTIFF)",
    )

    @field_validator("data_path", mode="before")
    @classmethod
    def normalize_data_path(cls, value: Any) -> str | Path | None:
        """Preserve remote URLs and normalize local paths."""
        if value is None or isinstance(value, Path):
            return value

        text = str(value).strip()
        if not text:
            return None

        if urlparse(text).scheme.lower() in {"http", "https"}:
            return text
        return Path(text).expanduser()

    @model_validator(mode="after")
    def validate_user_provider(self) -> AtmoPriorConfig:
        """Ensure user-provided paths are set when provider='user'."""
        if self.provider == "user" and self.user_aot is None and self.user_tcwv is None:
            raise ValueError(
                "When provider='user', at least user_aot or user_tcwv must be provided"
            )
        return self


class BRDFConfig(BaseModel):
    """Configuration for BRDF product data provider."""

    provider: Literal["mcd43", "vnp43", "mcd19", "gee", "zarr", "user"] = Field(
        default="mcd43",
        description="BRDF product source: 'mcd43' (MODIS), 'vnp43' (VIIRS), 'mcd19', 'gee', 'zarr', or 'user'",
    )
    data_path: Path | None = Field(
        default=None,
        description="Path to BRDF data directory (HDF files or Zarr store)",
    )
    temporal_window: int = Field(
        default=16,
        ge=1,
        le=32,
        description="Days before/after observation to include in temporal composite",
    )
    use_cache: bool = Field(
        default=True,
        description="Whether to cache BRDF data for faster repeated access",
    )
    cache_dir: Path | None = Field(
        default=None,
        description="Directory for caching BRDF data",
    )

    # GEE-specific settings
    gee_project: str | None = Field(
        default=None,
        description="Google Earth Engine project ID (for gee provider)",
    )

    # Zarr-specific settings
    zarr_url: str | None = Field(
        default=None,
        description="URL to Zarr store (for zarr provider)",
    )


class SurfacePriorConfig(BaseModel):
    """Configuration for surface prior derivation from BRDF."""

    method: Literal["kernel_model", "neural", "direct"] = Field(
        default="kernel_model",
        description="Surface prior derivation method",
    )
    psf_sigma_x: float = Field(
        default=29.75,
        gt=0,
        description="PSF standard deviation in x direction (high-res pixels)",
    )
    psf_sigma_y: float = Field(
        default=39.0,
        gt=0,
        description="PSF standard deviation in y direction (high-res pixels)",
    )
    apply_psf: bool = Field(
        default=True,
        description="Whether to apply PSF convolution for scale matching",
    )


class RTModelConfig(BaseModel):
    """Configuration for radiative transfer model backend."""

    backend: Literal["emulator", "lut", "py6s"] = Field(
        default="emulator",
        description="RT model backend: 'emulator' (NN), 'lut' (Zarr), or 'py6s' (direct)",
    )
    emulator_dir: Path | None = Field(
        default=None,
        description="Directory containing emulator .npz files",
    )
    lut_path: Path | str | None = Field(
        default=DEFAULT_LUT_URL,
        description=(
            "Path or URL to LUT store. Supports local .zarr directory, local .zarr.zip, "
            "HTTP(S) URLs, and s3:// URLs."
        ),
    )
    lut_storage_options: dict[str, Any] = Field(
        default_factory=dict,
        description=(
            "Extra fsspec storage options for LUT loading. For S3, you can pass "
            "{region, endpoint_url, key, secret, anon}."
        ),
    )
    lut_interpolation: Literal["linear", "nearest", "cubic"] = Field(
        default="linear",
        description="Interpolation method for LUT lookup",
    )

    # Py6S settings
    py6s_aero_profile: str = Field(
        default="Continental",
        description="Py6S aerosol profile (Continental, Maritime, Urban, etc.)",
    )

    # Fallback behavior
    fallback_to_lut: bool = Field(
        default=True,
        description="Fall back to LUT if emulator unavailable for sensor",
    )
    fallback_to_py6s: bool = Field(
        default=True,
        description="Fall back to Py6S if LUT unavailable",
    )


class SolverConfig(BaseModel):
    """Configuration for aerosol retrieval solver."""

    # Smoothness regularization
    aot_gamma: float = Field(
        default=10.0,
        ge=0,
        description="Smoothness weight for AOT (higher = smoother)",
    )
    tcwv_gamma: float = Field(
        default=5.0,
        ge=0,
        description="Smoothness weight for TCWV (higher = smoother)",
    )

    # Spectral weighting
    alpha: float = Field(
        default=-1.6,
        description="Spectral weighting exponent (negative emphasizes shorter wavelengths)",
    )

    # Optimization settings
    max_iterations: int = Field(
        default=300,
        ge=1,
        description="Maximum L-BFGS-B iterations per grid level",
    )
    gtol: float = Field(
        default=1e-2,
        gt=0,
        description="Gradient tolerance for convergence",
    )
    ftol: float = Field(
        default=1e-7,
        gt=0,
        description="Function tolerance for convergence",
    )

    # Parameter bounds
    aot_bounds: tuple[float, float] = Field(
        default=(0.001, 2.5),
        description="(min, max) bounds for AOT",
    )
    tcwv_bounds: tuple[float, float] = Field(
        default=(0.0, 7.0),
        description="(min, max) bounds for TCWV (g/cm²)",
    )

    # Resolution
    aerosol_resolution: float = Field(
        default=1000.0,
        gt=0,
        description="Resolution for aerosol retrieval grid (meters)",
    )

    # Multi-grid settings
    use_multigrid: bool = Field(
        default=True,
        description="Use multi-grid optimization (coarse to fine)",
    )
    min_grid_size: int = Field(
        default=4,
        ge=2,
        description="Minimum grid dimension for multi-grid",
    )

    @field_validator("aot_bounds", "tcwv_bounds")
    @classmethod
    def validate_bounds(cls, v: tuple[float, float]) -> tuple[float, float]:
        """Ensure bounds are valid (min < max)."""
        if v[0] >= v[1]:
            raise ValueError(f"Bounds must have min < max, got {v}")
        return v


class CredentialConfig(BaseModel):
    """Credentials for remote data providers.

    All fields are optional.  When ``SIACConfig`` is loaded via
    ``pydantic_settings`` the env-var prefix is
    ``SIAC_CREDENTIALS__<FIELD>`` (double-underscore nested delimiter).
    """

    cdse_username: str | None = None
    cdse_password: str | None = None
    cds_api_key: str | None = None
    aws_access_key_id: str | None = None
    aws_secret_access_key: str | None = None
    earthdata_username: str | None = None
    earthdata_password: str | None = None
    gcs_credentials_file: Path | None = None


class S2DataAccessConfig(BaseModel):
    """Configuration for Sentinel-2 pre-M1 data access."""

    backend: Literal["cdse", "gcs", "local"] = Field(
        default="local",
        description=(
            "S2 data source backend: 'cdse' (Copernicus Data Space), "
            "'gcs' (Google public bucket), or 'local'."
        ),
    )
    cache_dir: Path | None = Field(
        default=None,
        description="Local cache directory for downloaded SAFE products",
    )
    max_cloud_cover: float = Field(
        default=80.0,
        ge=0.0,
        le=100.0,
        description="Maximum cloud cover filter applied during S2 search",
    )
    prefer_newest_baseline: bool = Field(
        default=True,
        description="Prefer newest processing baseline during product selection",
    )
    processing_level: Literal["L1C", "L2A"] = Field(
        default="L1C",
        description="Sentinel-2 processing level to target for search/fetch",
    )
    cdse_access_key: str | None = Field(
        default=None,
        description="CDSE access key/username override for S2 data access",
    )
    cdse_secret_key: str | None = Field(
        default=None,
        description="CDSE secret/password override for S2 data access",
    )


class CloudMaskConfig(BaseModel):
    """Configuration for cloud/cloud-shadow class generation."""

    mode: Literal["auto", "external_file", "user_callable", "none"] = Field(
        default="auto",
        description=(
            "Cloud-mask source mode: 'auto' (run default detector), "
            "'external_file', 'user_callable', or 'none'."
        ),
    )
    provider: Literal["omnicloudmask"] = Field(
        default="omnicloudmask",
        description="Default cloud detector provider for mode='auto'.",
    )
    external_mask_path: Path | None = Field(
        default=None,
        description="Path to existing cloud/shadow raster for mode='external_file'.",
    )
    class_mapping: dict[int, list[int]] | None = Field(
        default=None,
        description=(
            "Mapping from standardized class keys {0,1,2,3} to source class values. "
            "Multiple source classes can map to one target class."
        ),
    )
    unmapped_to_missing: bool = Field(
        default=True,
        description="When True, source classes not present in mapping are set to class 0.",
    )
    target_resolution_m: float = Field(
        default=10.0,
        gt=0.0,
        description="Target resolution (meters) for cloud-mask generation.",
    )
    resolution_policy: Literal["auto", "force"] = Field(
        default="auto",
        description=(
            "'auto': downsample finer-than-target data; keep coarser by default. "
            "'force': always resample to target_resolution_m."
        ),
    )
    allow_upsample_to_target: bool = Field(
        default=False,
        description=(
            "If True and resolution_policy='auto', coarser data can be upsampled to target resolution."
        ),
    )
    # Runtime-only callable hook for Python users; not serialized in YAML workflows.
    user_callable: Any | None = Field(
        default=None,
        description="Callable returning cloud classes for mode='user_callable'.",
    )


class ExecutionConfig(BaseModel):
    """Configuration for orchestration backend and concurrency controls."""

    backend: Literal["thread", "dask"] = Field(
        default="thread",
        description="Execution backend for pipeline orchestration.",
    )
    max_workers: int = Field(
        default=4,
        ge=1,
        description="Maximum parallel workers for orchestration tasks.",
    )
    retries: int = Field(
        default=2,
        ge=0,
        description="Retry count for transient M2/M3 provider failures.",
    )
    stage_timeout_s: float | None = Field(
        default=None,
        gt=0.0,
        description="Optional timeout per concurrent provider stage in seconds.",
    )
    dashboard: bool = Field(
        default=False,
        description="Enable dashboard when backend='dask'.",
    )
    dashboard_address: str | None = Field(
        default=None,
        description="Dashboard bind address for dask backend.",
    )
    performance_report_path: Path | None = Field(
        default=None,
        description="Optional dask performance report HTML output path.",
    )
    show_progress: bool = Field(
        default=False,
        description="Emit backend progress endpoints/log hints.",
    )


class OutputConfig(BaseModel):
    """Configuration for output products."""

    output_dir: Path | None = Field(
        default=None,
        description="Output directory (defaults to input directory)",
    )
    format: Literal["geotiff", "cog", "zarr", "netcdf"] = Field(
        default="cog",
        description="Output format: 'geotiff', 'cog' (cloud-optimized), 'zarr', or 'netcdf'",
    )
    compression: Literal["deflate", "lzw", "zstd", "none"] = Field(
        default="deflate",
        description="Compression algorithm for GeoTIFF/COG",
    )
    include_rgb: bool = Field(
        default=True,
        description="Generate RGB quicklook image",
    )
    include_uncertainty: bool = Field(
        default=True,
        description="Include per-pixel uncertainty bands",
    )
    include_auxiliary: bool = Field(
        default=True,
        description="Output AOT and TCWV maps",
    )

    # Data type settings
    boa_dtype: Literal["float32", "float64", "uint16"] = Field(
        default="float32",
        description="Data type for BOA reflectance output",
    )
    boa_scale: float = Field(
        default=10000.0,
        gt=0,
        description="Scale factor for uint16 output (BOA * scale)",
    )
    boa_nodata: float = Field(
        default=0.0,
        description="NoData value for BOA output",
    )


# =============================================================================
# Main Configuration Model
# =============================================================================


class SIACConfig(BaseSettings):
    """
    Complete SIAC configuration.

    This is the main configuration class that aggregates all sub-configurations.
    Can be loaded from YAML files, environment variables, or created programmatically.

    Environment variables are prefixed with SIAC_ and use double underscore for
    nested values (e.g., SIAC_RT_MODEL__BACKEND=lut).
    """

    # Sensor selection
    sensor: Literal["s2", "l8", "l9", "auto"] = Field(
        default="auto",
        description="Sensor type: 's2' (Sentinel-2), 'l8' (Landsat 8), 'l9' (Landsat 9), or 'auto' (detect)",
    )

    # Module configurations
    atmo_prior: AtmoPriorConfig = Field(
        default_factory=AtmoPriorConfig,
        description="Atmospheric prior provider configuration",
    )
    brdf: BRDFConfig = Field(
        default_factory=BRDFConfig,
        description="BRDF product provider configuration",
    )
    surface_prior: SurfacePriorConfig = Field(
        default_factory=SurfacePriorConfig,
        description="Surface prior derivation configuration",
    )
    rt_model: RTModelConfig = Field(
        default_factory=RTModelConfig,
        description="Radiative transfer model backend configuration",
    )
    solver: SolverConfig = Field(
        default_factory=SolverConfig,
        description="Aerosol retrieval solver configuration",
    )
    output: OutputConfig = Field(
        default_factory=OutputConfig,
        description="Output product configuration",
    )
    credentials: CredentialConfig = Field(
        default_factory=CredentialConfig,
        description="Credentials for remote data providers (CDSE, CDS, AWS, Earthdata, GCS)",
    )
    s2_data: S2DataAccessConfig = Field(
        default_factory=S2DataAccessConfig,
        description="Sentinel-2 data access configuration (query -> local SAFE resolution)",
    )
    cloud_mask: CloudMaskConfig = Field(
        default_factory=CloudMaskConfig,
        description="Cloud/cloud-shadow mask generation configuration",
    )
    execution: ExecutionConfig = Field(
        default_factory=ExecutionConfig,
        description="Execution backend and concurrency settings",
    )

    # Global settings
    aoi: Path | str | None = Field(
        default=None,
        description="Area of interest (GeoJSON path, WKT string, or bounding box)",
    )
    global_dem: Path | str | None = Field(
        default=None,
        description="Path or URL to global DEM VRT/COG",
    )
    global_water: Path | str | None = Field(
        default=None,
        description="Path or URL to global water mask",
    )

    # Processing settings
    log_level: Literal["DEBUG", "INFO", "WARNING", "ERROR"] = Field(
        default="INFO",
        description="Logging verbosity level",
    )
    n_jobs: int = Field(
        default=-1,
        description="Number of parallel jobs (-1 for all CPUs)",
    )
    chunk_size: int = Field(
        default=2048,
        gt=0,
        description="Chunk size for dask/xarray processing",
    )

    model_config = {
        "env_prefix": "SIAC_",
        "env_nested_delimiter": "__",
        "extra": "forbid",
    }

    @classmethod
    def from_yaml(cls, path: Path | str) -> SIACConfig:
        """
        Load configuration from a YAML file.

        Args:
            path: Path to YAML configuration file

        Returns:
            SIACConfig instance

        Example YAML:
            sensor: auto
            atmo_prior:
              provider: cams
            brdf:
              provider: mcd43
              temporal_window: 16
            rt_model:
              backend: emulator
        """
        path = Path(path)
        with path.open() as f:
            data = yaml.safe_load(f)

        return cls(**data)

    def to_yaml(self, path: Path | str) -> None:
        """
        Save configuration to a YAML file.

        Args:
            path: Output path for YAML file
        """
        path = Path(path)

        # Convert to dict, handling Path objects
        data = self.model_dump(mode="json")

        with path.open("w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)

    def to_dict(self) -> dict[str, Any]:
        """Convert configuration to dictionary."""
        return self.model_dump()

    def with_overrides(self, **kwargs: Any) -> SIACConfig:
        """
        Create a new config with specified overrides.

        Args:
            **kwargs: Configuration values to override

        Returns:
            New SIACConfig with overrides applied
        """
        data = self.model_dump()

        # Deep merge overrides
        for key, value in kwargs.items():
            if isinstance(value, dict) and key in data and isinstance(data[key], dict):
                data[key].update(value)
            else:
                data[key] = value

        return SIACConfig(**data)


# =============================================================================
# Default Configuration Templates
# =============================================================================


def get_default_config() -> SIACConfig:
    """Get default SIAC configuration."""
    return SIACConfig()


def get_jasmin_config() -> SIACConfig:
    """Get SIAC configuration optimized for JASMIN HPC."""
    return SIACConfig(
        atmo_prior=AtmoPriorConfig(
            provider="cams",
            data_path=Path("/work/scratch-pw3/marc/CAMS/"),
        ),
        brdf=BRDFConfig(
            provider="mcd43",
            data_path=Path("/gws/nopw/j04/nceo_ard/public/MCD43/"),
        ),
        global_dem="/vsicurl/https://raw.githubusercontent.com/MarcYin/Copernicus_GLO_30_DEM_VRT/refs/heads/main/copernicus_GLO_30_dem.vrt",
        n_jobs=8,
    )


def get_lut_config(lut_path: str | Path) -> SIACConfig:
    """
    Get SIAC configuration using LUT backend.

    Suitable for sensors without pre-trained emulators.

    Args:
        lut_path: Path or URL to LUT Zarr store
    """
    return SIACConfig(
        rt_model=RTModelConfig(
            backend="lut",
            lut_path=_coerce_lut_source(lut_path),
        ),
    )


def _coerce_lut_source(lut_path: str | Path) -> str | Path:
    """Preserve URL-like LUT sources as strings, local paths as Path."""
    if isinstance(lut_path, Path):
        return lut_path

    lut_str = str(lut_path)
    if lut_str.startswith(("http://", "https://", "s3://", "file://")):
        return lut_str

    return Path(lut_str)
