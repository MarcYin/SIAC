"""
Configuration system for SIAC using Pydantic.

This module provides a hierarchical configuration system that allows users
to specify processing options via YAML files, environment variables, or
programmatically. All configuration is validated at load time.

Example usage:
    >>> from siac.core.config import SIACConfig
    >>>
    >>> # Load from TOML or YAML
    >>> config = SIACConfig.from_file("siac_config.toml")
    >>>
    >>> # Or create programmatically
    >>> config = SIACConfig(
    ...     sensor="s2",
    ...     rt_model={"backend": "lut", "lut_path": "/path/to/lut.zarr"},
    ... )
    >>>
    >>> # Save to TOML
    >>> config.to_toml("config.toml")
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal
from urllib.parse import urlparse

import yaml
from pydantic import BaseModel, Field, PrivateAttr, field_validator, model_validator
from pydantic_settings import BaseSettings

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11
    import tomli as tomllib

DEFAULT_LUT_URL = (
    "https://gws-access.jasmin.ac.uk/public/nceo_isp/"
    "libradtran_continental_average_lut_1nm.zarr.zip"
)

_REMOTE_URI_SCHEMES = {"http", "https", "s3", "file"}
_SNAPSHOT_REDACTION_TOKEN = "<redacted>"
_TRACKED_ENV_PREFIXES = ("SIAC_",)
_TRACKED_ENV_KEYS = (
    "CDSAPI_KEY",
    "AWS_ACCESS_KEY_ID",
    "AWS_SECRET_ACCESS_KEY",
    "EARTHDATA_USERNAME",
    "EARTHDATA_PASSWORD",
    "GOOGLE_APPLICATION_CREDENTIALS",
    "RSRF_ROOT",
)
_CATEGORIZED_CONFIG_KEYS = frozenset(
    {
        "inputs",
        "paths",
        "runtime",
        "execution",
        "providers",
        "processing",
        "models",
        "auth",
        "output",
    }
)


def _coerce_pathlike(value: Any) -> Path | None:
    if value is None or isinstance(value, Path):
        return value

    text = str(value).strip()
    if not text:
        return None
    return Path(text).expanduser()


def _coerce_path_or_url(value: Any) -> str | Path | None:
    if value is None or isinstance(value, Path):
        return value

    text = str(value).strip()
    if not text:
        return None

    if urlparse(text).scheme.lower() in _REMOTE_URI_SCHEMES:
        return text
    return Path(text).expanduser()


def _normalize_snapshot_value(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, tuple):
        return [_normalize_snapshot_value(item) for item in value]
    if isinstance(value, list):
        return [_normalize_snapshot_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _normalize_snapshot_value(item) for key, item in value.items()}
    return value


def _path_state(value: Any) -> dict[str, Any]:
    if value is None:
        return {"value": None, "kind": "unset", "exists": None}

    normalized = _normalize_snapshot_value(value)
    if isinstance(normalized, str) and urlparse(normalized).scheme.lower() in _REMOTE_URI_SCHEMES:
        return {"value": normalized, "kind": "remote", "exists": None}

    path = Path(normalized).expanduser() if not isinstance(normalized, Path) else normalized.expanduser()
    resolved = path.resolve(strict=False)
    if resolved.exists():
        kind = "directory" if resolved.is_dir() else "file"
        exists = True
    else:
        kind = "missing"
        exists = False
    return {"value": str(resolved), "kind": kind, "exists": exists}


def _deep_merge(base: dict[str, Any], updates: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def _tracked_environment_keys() -> list[str]:
    keys = [key for key in os.environ if any(key.startswith(prefix) for prefix in _TRACKED_ENV_PREFIXES)]
    keys.extend(key for key in _TRACKED_ENV_KEYS if key in os.environ)
    return sorted(set(keys))


def _looks_like_categorized_config(data: dict[str, Any]) -> bool:
    return any(key in data for key in _CATEGORIZED_CONFIG_KEYS)


def _flatten_categorized_config(data: dict[str, Any]) -> dict[str, Any]:
    if not _looks_like_categorized_config(data):
        return data

    flattened: dict[str, Any] = {}

    inputs = data.get("inputs")
    if inputs is not None:
        if not isinstance(inputs, dict):
            raise TypeError("The 'inputs' config section must be a mapping/object.")
        for key in ("sensor", "aoi", "s2_data"):
            if key in inputs:
                flattened[key] = inputs[key]

    if "paths" in data:
        flattened["paths"] = data["paths"]

    runtime = data.get("runtime")
    if runtime is not None:
        if not isinstance(runtime, dict):
            raise TypeError("The 'runtime' config section must be a mapping/object.")
        runtime_settings = runtime.get("settings")
        if runtime_settings is not None:
            flattened["runtime"] = runtime_settings
        elif any(key in runtime for key in RuntimeSettingsConfig.model_fields):
            flattened["runtime"] = runtime

        runtime_execution = runtime.get("execution")
        if runtime_execution is not None:
            flattened["execution"] = runtime_execution

    if "execution" in data:
        flattened["execution"] = data["execution"]

    providers = data.get("providers")
    if providers is not None:
        if not isinstance(providers, dict):
            raise TypeError("The 'providers' config section must be a mapping/object.")
        for key in ("atmo_prior", "brdf"):
            if key in providers:
                flattened[key] = providers[key]

    processing = data.get("processing")
    if processing is not None:
        if not isinstance(processing, dict):
            raise TypeError("The 'processing' config section must be a mapping/object.")
        if "cloud_mask" in processing:
            flattened["cloud_mask"] = processing["cloud_mask"]

    models = data.get("models")
    if models is not None:
        if not isinstance(models, dict):
            raise TypeError("The 'models' config section must be a mapping/object.")
        for key in ("surface_prior", "rt_model", "solver"):
            if key in models:
                flattened[key] = models[key]

    auth = data.get("auth")
    if auth is not None:
        if not isinstance(auth, dict):
            raise TypeError("The 'auth' config section must be a mapping/object.")
        flattened["credentials"] = auth.get("credentials", auth)

    if "output" in data:
        flattened["output"] = data["output"]

    for key, value in data.items():
        if key not in _CATEGORIZED_CONFIG_KEYS:
            flattened[key] = value

    return flattened


def _toml_key(key: Any) -> str:
    text = str(key)
    if text and (text[0].isalpha() or text[0] in {"_", "-"}) and all(
        char.isalnum() or char in {"_", "-"} for char in text
    ):
        return text
    return json.dumps(text)


def _toml_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int) and not isinstance(value, bool):
        return str(value)
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, str):
        return json.dumps(value)
    if isinstance(value, list):
        return "[" + ", ".join(_toml_value(item) for item in value) + "]"
    raise TypeError(f"Unsupported TOML value type: {type(value)!r}")


def _toml_lines(data: dict[str, Any], prefix: tuple[str, ...] = ()) -> list[str]:
    lines: list[str] = []
    scalar_items: list[tuple[Any, Any]] = []
    table_items: list[tuple[str, dict[str, Any]]] = []

    for key, value in data.items():
        if value is None:
            continue
        if isinstance(value, dict):
            nested_lines = _toml_lines(value, prefix + (str(key),))
            if nested_lines:
                table_items.append((str(key), value))
        else:
            scalar_items.append((key, value))

    if prefix:
        lines.append("[" + ".".join(_toml_key(part) for part in prefix) + "]")

    for key, value in scalar_items:
        lines.append(f"{_toml_key(key)} = {_toml_value(value)}")

    if scalar_items and table_items:
        lines.append("")

    for index, (key, value) in enumerate(table_items):
        lines.extend(_toml_lines(value, prefix + (key,)))
        if index < len(table_items) - 1:
            lines.append("")

    return lines


def _dump_toml(data: dict[str, Any]) -> str:
    lines = _toml_lines(data)
    if not lines:
        return ""
    return "\n".join(lines) + "\n"


def _load_config_mapping(path: Path) -> dict[str, Any]:
    suffix = path.suffix.lower()
    if suffix == ".toml":
        with path.open("rb") as handle:
            loaded = tomllib.load(handle) or {}
    elif suffix in {".yaml", ".yml", ""}:
        with path.open(encoding="utf-8") as handle:
            loaded = yaml.safe_load(handle) or {}
    else:
        raise ValueError(f"Unsupported config format for '{path}'. Use .toml, .yaml, or .yml.")

    if not isinstance(loaded, dict):
        raise ValueError("Configuration file must contain a mapping/object at the top level.")
    return loaded


# =============================================================================
# Sub-configuration Models
# =============================================================================


class AtmoPriorConfig(BaseModel):
    """Configuration for atmospheric prior data provider."""

    provider: Literal["cams", "merra2", "mcd19", "vnp19", "era5", "user"] = Field(
        default="cams",
        description=(
            "Atmospheric data source: 'cams' (ECMWF), 'merra2' (NASA), "
            "'mcd19' (MODIS MAIAC), 'vnp19' (VIIRS MAIAC), 'era5', or 'user'"
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

        if urlparse(text).scheme.lower() in {"http", "https", "s3"}:
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
        description=(
            "BRDF product source: 'mcd43' (MODIS), 'vnp43' (VIIRS), "
            "'mcd19' (MAIAC kernels), 'gee', 'zarr', or 'user'"
        ),
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


class SurfaceSpectralMappingSettings(BaseModel):
    """Configuration for multisensor spectral mapping used by surface-prior flows."""

    enabled: bool = Field(
        default=True,
        description="Enable spectral mapping when BRDF/source and target sensor bands differ.",
    )
    siac_library_root: Path | None = Field(
        default=None,
        description="Path to the prepared SIAC hyperspectral library export root.",
    )
    rsrf_root: Path | None = Field(
        default=None,
        description="Path to the RSRF repository checkout/root.",
    )
    cache_dir: Path | None = Field(
        default=None,
        description="Cache directory for prepared spectral-mapping runtimes.",
    )
    k_neighbors: int = Field(
        default=5,
        ge=1,
        description="Number of hyperspectral neighbors used during spectral mapping.",
    )
    neighbor_estimator: str = Field(
        default="distance_weighted_mean",
        description="Neighbor combination strategy used by the spectral-library mapper.",
    )
    knn_backend: str = Field(
        default="numpy",
        description="KNN backend used by spectral-library.",
    )
    knn_eps: float = Field(
        default=0.0,
        ge=0.0,
        description="Approximate-neighbor epsilon passed to the spectral-library backend.",
    )
    min_valid_bands: int = Field(
        default=1,
        ge=1,
        description="Minimum count of valid source bands required for mapping.",
    )

    @field_validator("siac_library_root", "rsrf_root", "cache_dir", mode="before")
    @classmethod
    def normalize_local_paths(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)


class SurfacePriorConfig(BaseModel):
    """Configuration for surface prior derivation from BRDF."""

    method: Literal["kernel_model", "whittaker", "monthly_database", "neural", "direct"] = Field(
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
    whittaker_lambda: float = Field(
        default=10.0,
        gt=0,
        description="Temporal smoothness strength for Route-A Whittaker BRDF priors",
    )
    spectral_mapping: SurfaceSpectralMappingSettings = Field(
        default_factory=SurfaceSpectralMappingSettings,
        description="Spectral mapping configuration for cross-sensor BRDF/surface-prior transfers.",
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

    @field_validator("emulator_dir", mode="before")
    @classmethod
    def normalize_emulator_dir(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)

    @field_validator("lut_path", mode="before")
    @classmethod
    def normalize_lut_path(cls, value: Any) -> str | Path | None:
        return _coerce_path_or_url(value)


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

    All fields are optional explicit overrides. Provider-native environment
    variables are preferred for normal runtime configuration; these fields are
    primarily for programmatic or YAML-driven injection.

    When ``SIACConfig`` is loaded via ``pydantic_settings`` the nested env-var
    prefix is ``SIAC_CREDENTIALS__<FIELD>``.
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

    @field_validator("cache_dir", mode="before")
    @classmethod
    def normalize_cache_dir(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)


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


class RuntimeSettingsConfig(BaseModel):
    """General SIAC runtime settings that are not specific to one processing stage."""

    log_level: Literal["DEBUG", "INFO", "WARNING", "ERROR"] = Field(
        default="INFO",
        description="Logging verbosity level.",
    )
    n_jobs: int = Field(
        default=-1,
        description="Number of parallel jobs (-1 for all CPUs).",
    )
    chunk_size: int = Field(
        default=2048,
        gt=0,
        description="Chunk size for dask/xarray processing.",
    )


class PackagePathsConfig(BaseModel):
    """Centralized package asset/path configuration."""

    global_dem: str | Path | None = Field(
        default=None,
        description="Path or URL to the global DEM VRT/COG.",
    )
    global_water: str | Path | None = Field(
        default=None,
        description="Path or URL to the global water mask.",
    )
    emulator_dir: Path | None = Field(
        default=None,
        description="Default directory containing RT emulator weight files.",
    )
    lut_path: str | Path | None = Field(
        default=None,
        description="Default path or URL to the RT LUT store.",
    )
    spectral_library_root: Path | None = Field(
        default=None,
        description="Root path to the SIAC spectral library export.",
    )
    rsrf_root: Path | None = Field(
        default=None,
        description="Root path to the RSRF repository data checkout.",
    )
    spectral_mapping_cache_dir: Path | None = Field(
        default=None,
        description="Cache directory for prepared spectral-mapping runtimes.",
    )

    @field_validator("global_dem", "global_water", "lut_path", mode="before")
    @classmethod
    def normalize_path_or_url_fields(cls, value: Any) -> str | Path | None:
        return _coerce_path_or_url(value)

    @field_validator(
        "emulator_dir",
        "spectral_library_root",
        "rsrf_root",
        "spectral_mapping_cache_dir",
        mode="before",
    )
    @classmethod
    def normalize_local_path_fields(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)


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
    paths: PackagePathsConfig = Field(
        default_factory=PackagePathsConfig,
        description="Centralized package asset and path configuration.",
    )
    runtime: RuntimeSettingsConfig = Field(
        default_factory=RuntimeSettingsConfig,
        description="General runtime/logging/chunking configuration.",
    )

    # Global settings
    aoi: Path | str | None = Field(
        default=None,
        description=(
            "Optional default area of interest (GeoJSON path, WKT string, or bounding box). "
            "Prefer passing AOI per run when it varies by scene/tile."
        ),
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

    _config_file_path: Path | None = PrivateAttr(default=None)

    model_config = {
        "env_prefix": "SIAC_",
        "env_nested_delimiter": "__",
        "extra": "forbid",
    }

    @model_validator(mode="before")
    @classmethod
    def normalize_categorized_payload(cls, data: Any) -> Any:
        """Accept the canonical categorized YAML shape alongside flat kwargs."""
        if isinstance(data, dict):
            return _flatten_categorized_config(data)
        return data

    @model_validator(mode="after")
    def synchronize_categorized_settings(self) -> SIACConfig:
        """Keep categorized fields and legacy mirrors aligned.

        The categorized sections are treated as canonical when they are set away
        from their defaults; otherwise legacy mirrors populate them.
        """
        default_paths = PackagePathsConfig()
        if self.paths != default_paths:
            if self.paths.global_dem is not None:
                self.global_dem = self.paths.global_dem
            elif self.global_dem is not None:
                self.paths.global_dem = self.global_dem

            if self.paths.global_water is not None:
                self.global_water = self.paths.global_water
            elif self.global_water is not None:
                self.paths.global_water = self.global_water
        else:
            self.paths = PackagePathsConfig(
                global_dem=self.global_dem,
                global_water=self.global_water,
            )

        default_runtime = RuntimeSettingsConfig()
        if self.runtime != default_runtime:
            self.log_level = self.runtime.log_level
            self.n_jobs = self.runtime.n_jobs
            self.chunk_size = self.runtime.chunk_size
        else:
            self.runtime = RuntimeSettingsConfig(
                log_level=self.log_level,
                n_jobs=self.n_jobs,
                chunk_size=self.chunk_size,
            )
        return self

    def resolved_global_dem(self) -> str | Path | None:
        return self.paths.global_dem if self.paths.global_dem is not None else self.global_dem

    def resolved_global_water(self) -> str | Path | None:
        return self.paths.global_water if self.paths.global_water is not None else self.global_water

    def resolved_emulator_dir(self) -> Path | None:
        if self.rt_model.emulator_dir is not None:
            return self.rt_model.emulator_dir
        return self.paths.emulator_dir

    def resolved_lut_path(self) -> str | Path | None:
        if self.rt_model.lut_path not in (None, DEFAULT_LUT_URL):
            return self.rt_model.lut_path
        if self.paths.lut_path is not None:
            return self.paths.lut_path
        return self.rt_model.lut_path

    def resolved_surface_spectral_mapping(self) -> SurfaceSpectralMappingSettings:
        settings = self.surface_prior.spectral_mapping.model_copy(deep=True)
        if settings.siac_library_root is None:
            settings.siac_library_root = self.paths.spectral_library_root
        if settings.rsrf_root is None:
            settings.rsrf_root = self.paths.rsrf_root
        if settings.cache_dir is None:
            settings.cache_dir = self.paths.spectral_mapping_cache_dir
        return settings

    def _canonical_dump(self, *, mode: Literal["python", "json"] = "python") -> dict[str, Any]:
        return {
            "inputs": {
                "sensor": self.sensor,
                "aoi": _normalize_snapshot_value(self.aoi) if mode == "json" else self.aoi,
                "s2_data": self.s2_data.model_dump(mode=mode),
            },
            "paths": self.paths.model_dump(mode=mode),
            "runtime": {
                "settings": self.runtime.model_dump(mode=mode),
                "execution": self.execution.model_dump(mode=mode),
            },
            "providers": {
                "atmo_prior": self.atmo_prior.model_dump(mode=mode),
                "brdf": self.brdf.model_dump(mode=mode),
            },
            "processing": {
                "cloud_mask": self.cloud_mask.model_dump(mode=mode, exclude={"user_callable"}),
            },
            "models": {
                "surface_prior": self.surface_prior.model_dump(mode=mode),
                "rt_model": self.rt_model.model_dump(mode=mode),
                "solver": self.solver.model_dump(mode=mode),
            },
            "auth": self.credentials.model_dump(mode=mode),
            "output": self.output.model_dump(mode=mode),
        }

    def categorized_state(self, *, redact_secrets: bool = True) -> dict[str, Any]:
        config_state = self._canonical_dump(mode="json")
        auth_state = config_state["auth"]
        s2_state = config_state["inputs"]["s2_data"]
        if redact_secrets:
            for payload, secret_keys in (
                (
                    auth_state,
                    {
                        "cdse_username",
                        "cdse_password",
                        "cds_api_key",
                        "aws_access_key_id",
                        "aws_secret_access_key",
                        "earthdata_username",
                        "earthdata_password",
                        "gcs_credentials_file",
                    },
                ),
                (
                    s2_state,
                    {
                        "cdse_access_key",
                        "cdse_secret_key",
                    },
                ),
            ):
                for key in secret_keys:
                    if payload.get(key) is not None:
                        payload[key] = _SNAPSHOT_REDACTION_TOKEN

        package_paths_state = {
            **self.paths.model_dump(mode="json"),
            "global_dem_resolved": _normalize_snapshot_value(self.resolved_global_dem()),
            "global_water_resolved": _normalize_snapshot_value(self.resolved_global_water()),
            "emulator_dir_resolved": _normalize_snapshot_value(self.resolved_emulator_dir()),
            "lut_path_resolved": _normalize_snapshot_value(self.resolved_lut_path()),
        }
        spectral_mapping_state = self.resolved_surface_spectral_mapping().model_dump(mode="json")

        return {
            "load": {
                "config_file": None if self._config_file_path is None else str(self._config_file_path),
                "active_environment_keys": _tracked_environment_keys(),
            },
            "config": config_state,
            "resolved": {
                "paths": package_paths_state,
                "surface_spectral_mapping": spectral_mapping_state,
                "path_states": {
                    "paths.global_dem": _path_state(self.resolved_global_dem()),
                    "paths.global_water": _path_state(self.resolved_global_water()),
                    "paths.emulator_dir": _path_state(self.resolved_emulator_dir()),
                    "paths.lut_path": _path_state(self.resolved_lut_path()),
                    "paths.spectral_library_root": _path_state(spectral_mapping_state["siac_library_root"]),
                    "paths.rsrf_root": _path_state(spectral_mapping_state["rsrf_root"]),
                    "paths.spectral_mapping_cache_dir": _path_state(spectral_mapping_state["cache_dir"]),
                    "s2_data.cache_dir": _path_state(self.s2_data.cache_dir),
                    "atmo_prior.cache_dir": _path_state(self.atmo_prior.cache_dir),
                    "brdf.cache_dir": _path_state(self.brdf.cache_dir),
                    "output.output_dir": _path_state(self.output.output_dir),
                },
            },
        }

    def write_state_snapshot(self, path: Path | str, *, redact_secrets: bool = True) -> None:
        path = Path(path)
        snapshot = {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "state": self.categorized_state(redact_secrets=redact_secrets),
        }
        with path.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(snapshot, handle, default_flow_style=False, sort_keys=False)

    @classmethod
    def default_config_path(cls, format: Literal["toml", "yaml"] = "toml") -> Path:
        """Return the recommended user-level SIAC config path."""
        suffix = ".toml" if format == "toml" else ".yaml"
        return Path.home() / ".config" / "siac" / f"config{suffix}"

    @classmethod
    def write_default_config(
        cls,
        path: Path | str | None = None,
        *,
        format: Literal["toml", "yaml"] = "toml",
        overwrite: bool = False,
    ) -> Path:
        """Write a default categorized config file and return its path."""
        resolved = Path(path) if path is not None else cls.default_config_path(format=format)
        if not resolved.suffix:
            resolved = resolved.with_suffix(".toml" if format == "toml" else ".yaml")
        if resolved.exists() and not overwrite:
            raise FileExistsError(f"Config file already exists: {resolved}")
        resolved.parent.mkdir(parents=True, exist_ok=True)
        cls().to_file(resolved)
        return resolved

    @classmethod
    def load(cls, path: Path | str | None = None, **overrides: Any) -> SIACConfig:
        """Load config from categorized YAML/env sources and track the source path."""
        source_path: Path | None = None
        data: dict[str, Any] = {}

        env_config_path = os.getenv("SIAC_CONFIG_FILE")
        if path is None and env_config_path:
            path = env_config_path

        if path is not None:
            source_path = Path(path)
            loaded = _load_config_mapping(source_path)
            if isinstance(loaded.get("config"), dict):
                data = loaded["config"]
            elif isinstance(loaded.get("state"), dict) and isinstance(loaded["state"].get("config"), dict):
                data = loaded["state"]["config"]
            else:
                data = loaded

        if overrides:
            data = _deep_merge(data, overrides)

        config = cls(**data)
        config._config_file_path = source_path
        return config

    @classmethod
    def from_yaml(cls, path: Path | str) -> SIACConfig:
        """
        Load configuration from a YAML file.

        Args:
            path: Path to YAML configuration file

        Returns:
            SIACConfig instance

        Example YAML:
            inputs:
              sensor: auto
            providers:
              atmo_prior:
                provider: cams
              brdf:
                provider: mcd43
                temporal_window: 16
            models:
              rt_model:
                backend: emulator
        """
        return cls.load(path)

    @classmethod
    def from_toml(cls, path: Path | str) -> SIACConfig:
        """Load configuration from a TOML file."""
        return cls.load(path)

    @classmethod
    def from_file(cls, path: Path | str) -> SIACConfig:
        """Load configuration from a TOML, YAML, or YML file."""
        return cls.load(path)

    def to_yaml(self, path: Path | str) -> None:
        """
        Save configuration to a YAML file.

        Args:
            path: Output path for YAML file
        """
        path = Path(path)

        data = self._canonical_dump(mode="json")

        with path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(data, f, default_flow_style=False, sort_keys=False)

    def to_toml(self, path: Path | str) -> None:
        """Save configuration to a TOML file."""
        path = Path(path)
        data = self._canonical_dump(mode="json")
        with path.open("w", encoding="utf-8") as handle:
            handle.write(_dump_toml(data))

    def to_file(self, path: Path | str) -> None:
        """Save configuration using the file extension to select TOML or YAML."""
        path = Path(path)
        suffix = path.suffix.lower()
        if suffix == ".toml":
            self.to_toml(path)
            return
        if suffix in {".yaml", ".yml", ""}:
            self.to_yaml(path)
            return
        raise ValueError(f"Unsupported config format for '{path}'. Use .toml, .yaml, or .yml.")

    def to_dict(self) -> dict[str, Any]:
        """Convert configuration to dictionary."""
        return self._canonical_dump()

    def with_overrides(self, **kwargs: Any) -> SIACConfig:
        """
        Create a new config with specified overrides.

        Args:
            **kwargs: Configuration values to override

        Returns:
            New SIACConfig with overrides applied
        """
        data = self._canonical_dump()
        return SIACConfig(**_deep_merge(data, kwargs))


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
        paths=PackagePathsConfig(
            global_dem="/vsicurl/https://raw.githubusercontent.com/MarcYin/Copernicus_GLO_30_DEM_VRT/refs/heads/main/copernicus_GLO_30_dem.vrt",
        ),
        runtime=RuntimeSettingsConfig(n_jobs=8),
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
