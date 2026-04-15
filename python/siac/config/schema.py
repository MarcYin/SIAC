"""
Typed configuration schema for SIAC.

This module contains only data models. File I/O, environment overlays,
resolution, and snapshotting live in sibling modules.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal
from urllib.parse import urlparse

from pydantic import AliasChoices, BaseModel, Field, field_validator, model_validator

DEFAULT_LUT_URL = (
    "https://gws-access.jasmin.ac.uk/public/nceo_isp/"
    "libradtran_continental_average_lut_1nm.zarr.zip"
)

SensorName = Literal["s2", "l8", "l9", "auto"]
LogLevel = Literal["DEBUG", "INFO", "WARNING", "ERROR"]
AtmoProviderKind = Literal["cams", "merra2", "mcd19", "vnp19", "era5", "user"]
BRDFProviderKind = Literal["mcd43", "vnp43", "mcd19", "gee", "zarr", "user"]
MonthlyCompositeProviderKind = Literal["generated_brdf", "user_callable", "prepared_store"]
AtmosphericParameterName = Literal["aot", "tcwv", "tco3"]

_REMOTE_URI_SCHEMES = {"http", "https", "s3", "file", "gs"}


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
    if text.startswith("/vsi"):
        return text
    if urlparse(text).scheme.lower() in _REMOTE_URI_SCHEMES:
        return text
    return Path(text).expanduser()


class SIACBaseModel(BaseModel):
    model_config = {
        "extra": "forbid",
    }


class CachePathsConfig(SIACBaseModel):
    atmo: Path | None = None
    brdf: Path | None = None
    s2: Path | None = None

    @field_validator("atmo", "brdf", "s2", mode="before")
    @classmethod
    def normalize_paths(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)


class PathsConfig(SIACBaseModel):
    dem: str | Path | None = None
    water_mask: str | Path | None = None
    emulator_dir: Path | None = None
    lut_path: str | Path | None = DEFAULT_LUT_URL
    rsrf_root: Path | None = None
    cache_root: Path | None = None
    caches: CachePathsConfig = Field(default_factory=CachePathsConfig)

    @field_validator("dem", "water_mask", "lut_path", mode="before")
    @classmethod
    def normalize_path_or_url_fields(cls, value: Any) -> str | Path | None:
        return _coerce_path_or_url(value)

    @field_validator("emulator_dir", "rsrf_root", "cache_root", mode="before")
    @classmethod
    def normalize_local_paths(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)


class CDSEAuthConfig(SIACBaseModel):
    username: str | None = None
    password: str | None = None
    username_env: str = "SIAC_CDSE_USERNAME"
    password_env: str = "SIAC_CDSE_PASSWORD"


class CDSAuthConfig(SIACBaseModel):
    api_key: str | None = None
    api_key_env: str = "CDSAPI_KEY"


class AWSAuthConfig(SIACBaseModel):
    access_key_id: str | None = None
    secret_access_key: str | None = None
    access_key_id_env: str = "AWS_ACCESS_KEY_ID"
    secret_access_key_env: str = "AWS_SECRET_ACCESS_KEY"


class EarthdataAuthConfig(SIACBaseModel):
    username: str | None = None
    password: str | None = None
    username_env: str = "EARTHDATA_USERNAME"
    password_env: str = "EARTHDATA_PASSWORD"


class GCSAuthConfig(SIACBaseModel):
    credentials_file: Path | None = None
    credentials_file_env: str = "GOOGLE_APPLICATION_CREDENTIALS"

    @field_validator("credentials_file", mode="before")
    @classmethod
    def normalize_credentials_file(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)


class AuthConfig(SIACBaseModel):
    cdse: CDSEAuthConfig = Field(default_factory=CDSEAuthConfig)
    cds: CDSAuthConfig = Field(default_factory=CDSAuthConfig)
    aws: AWSAuthConfig = Field(default_factory=AWSAuthConfig)
    earthdata: EarthdataAuthConfig = Field(default_factory=EarthdataAuthConfig)
    gcs: GCSAuthConfig = Field(default_factory=GCSAuthConfig)


class AtmoProviderConfig(SIACBaseModel):
    kind: AtmoProviderKind = Field(
        default="cams",
        validation_alias=AliasChoices("kind", "provider"),
        serialization_alias="kind",
    )
    data_path: str | Path | None = None
    cache_dir: Path | None = None
    download_missing: bool = True
    temporal_interpolation: Literal["nearest", "linear"] = "nearest"
    user_aot: Path | None = None
    user_tcwv: Path | None = None
    user_tco3: Path | None = None

    @field_validator("data_path", mode="before")
    @classmethod
    def normalize_data_path(cls, value: Any) -> str | Path | None:
        return _coerce_path_or_url(value)

    @field_validator("cache_dir", "user_aot", "user_tcwv", "user_tco3", mode="before")
    @classmethod
    def normalize_local_paths(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)

    @model_validator(mode="after")
    def validate_user_provider(self) -> AtmoProviderConfig:
        if self.kind == "user" and self.user_aot is None and self.user_tcwv is None:
            raise ValueError(
                "When kind='user', at least one of user_aot or user_tcwv must be provided."
            )
        return self

    @property
    def provider(self) -> AtmoProviderKind:
        return self.kind


class BRDFProviderConfig(SIACBaseModel):
    kind: BRDFProviderKind = Field(
        default="mcd43",
        validation_alias=AliasChoices("kind", "provider"),
        serialization_alias="kind",
    )
    data_path: str | Path | None = None
    cache_dir: Path | None = None
    temporal_window: int = Field(default=16, ge=1, le=32)
    use_cache: bool = True
    gee_project: str | None = None
    zarr_url: str | None = None

    @field_validator("data_path", mode="before")
    @classmethod
    def normalize_data_path(cls, value: Any) -> str | Path | None:
        return _coerce_path_or_url(value)

    @field_validator("cache_dir", mode="before")
    @classmethod
    def normalize_cache_dir(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)

    @property
    def provider(self) -> BRDFProviderKind:
        return self.kind


class S2ProviderConfig(SIACBaseModel):
    backend: Literal["cdse", "gcs", "local"] = "local"
    cache_dir: Path | None = None
    max_cloud_cover: float = Field(default=80.0, ge=0.0, le=100.0)
    prefer_newest_baseline: bool = True
    processing_level: Literal["L1C", "L2A"] = "L1C"

    @field_validator("cache_dir", mode="before")
    @classmethod
    def normalize_cache_dir(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)


class MonthlyCompositeProviderConfig(SIACBaseModel):
    kind: MonthlyCompositeProviderKind = Field(
        default="generated_brdf",
        validation_alias=AliasChoices("kind", "provider"),
        serialization_alias="kind",
    )
    user_callable: Any | None = None
    store_path: Path | None = None
    strict_coverage: bool = True

    @field_validator("store_path", mode="before")
    @classmethod
    def normalize_store_path(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)

    @property
    def provider(self) -> MonthlyCompositeProviderKind:
        return self.kind


class ProvidersConfig(SIACBaseModel):
    atmo: AtmoProviderConfig = Field(default_factory=AtmoProviderConfig)
    brdf: BRDFProviderConfig = Field(default_factory=BRDFProviderConfig)
    s2: S2ProviderConfig = Field(default_factory=S2ProviderConfig)
    monthly_composites: MonthlyCompositeProviderConfig = Field(default_factory=MonthlyCompositeProviderConfig)


class SpectralMappingAlgorithmConfig(SIACBaseModel):
    enabled: bool = True
    k_neighbors: int = Field(default=5, ge=1)
    neighbor_estimator: str = "distance_weighted_mean"
    knn_backend: str = "scipy_ckdtree"
    knn_eps: float = Field(default=0.0, ge=0.0)
    min_valid_bands: int = Field(default=1, ge=1)


class MonthlyDatabaseQualityFilterConfig(SIACBaseModel):
    enabled: bool = True
    max_prediction_uncertainty: float = Field(default=0.05, gt=0.0)
    max_composite_quality: float = Field(default=0.05, gt=0.0)
    max_source_fit_rmse: float = Field(default=0.05, gt=0.0)
    max_knn_feature_distance: float = Field(default=0.05, gt=0.0)


class SurfacePriorAlgorithmConfig(SIACBaseModel):
    method: Literal["kernel_model", "whittaker", "monthly_database", "neural", "direct"] = (
        "kernel_model"
    )
    monthly_database_resolution_policy: Literal["provider_or_coarser", "aerosol"] = (
        "provider_or_coarser"
    )
    psf_sigma_x: float = Field(default=29.75, gt=0.0)
    psf_sigma_y: float = Field(default=39.0, gt=0.0)
    apply_psf: bool = True
    whittaker_lambda: float = Field(default=10.0, gt=0.0)
    spectral_mapping: SpectralMappingAlgorithmConfig = Field(
        default_factory=SpectralMappingAlgorithmConfig,
    )
    monthly_database_filter: MonthlyDatabaseQualityFilterConfig = Field(
        default_factory=MonthlyDatabaseQualityFilterConfig,
    )


class RTAlgorithmConfig(SIACBaseModel):
    backend: Literal["emulator", "lut", "py6s"] = "emulator"
    lut_interpolation: Literal["linear", "nearest", "cubic"] = "linear"
    lut_storage_options: dict[str, Any] = Field(default_factory=dict)
    py6s_aero_profile: str = "Continental"
    fallback_to_lut: bool = True
    fallback_to_py6s: bool = True


class SolverBoundsConfig(SIACBaseModel):
    aot: tuple[float, float] = (0.001, 2.5)
    tcwv: tuple[float, float] = (0.0, 7.0)

    @field_validator("aot", "tcwv")
    @classmethod
    def validate_bounds(cls, value: tuple[float, float]) -> tuple[float, float]:
        if value[0] >= value[1]:
            raise ValueError(f"Bounds must have min < max, got {value}")
        return value


class SolverStageConfig(SIACBaseModel):
    name: str = "default"
    solve: tuple[AtmosphericParameterName, ...] = ("aot", "tcwv")
    fixed: tuple[AtmosphericParameterName, ...] = ("tco3",)
    bands: tuple[str, ...] | None = None
    initial_state: Literal["prior", "previous"] = "previous"

    @field_validator("solve", "fixed", mode="before")
    @classmethod
    def normalize_parameters(cls, value: Any) -> Any:
        if value is None:
            return ()
        if isinstance(value, str):
            return (value,)
        return value

    @field_validator("bands", mode="before")
    @classmethod
    def normalize_bands(cls, value: Any) -> Any:
        if value is None or isinstance(value, tuple):
            return value
        if isinstance(value, str):
            return (value,)
        return tuple(value)

    @model_validator(mode="after")
    def validate_stage(self) -> SolverStageConfig:
        duplicated_solve = {name for name in self.solve if self.solve.count(name) > 1}
        duplicated_fixed = {name for name in self.fixed if self.fixed.count(name) > 1}
        duplicated = duplicated_solve | duplicated_fixed
        if duplicated:
            raise ValueError(
                f"Solver stage {self.name!r} duplicates parameter(s): "
                f"{', '.join(sorted(duplicated))}"
            )

        overlap = set(self.solve) & set(self.fixed)
        if overlap:
            raise ValueError(
                f"Solver stage {self.name!r} cannot both solve and fix "
                f"{', '.join(sorted(overlap))}"
            )
        return self


class SharpTransitionFilterConfig(SIACBaseModel):
    enabled: bool = False
    blur_kernel_pixels_native: int = Field(default=31, ge=3)
    residual_threshold_uint8: int = Field(default=12, ge=0, le=255)
    dilation_pixels: int = Field(default=1, ge=0)
    solver_cell_fraction_threshold: float = Field(default=0.03, ge=0.0, le=1.0)
    cloud_buffer_pixels: int = Field(default=2, ge=0)

    @field_validator("blur_kernel_pixels_native")
    @classmethod
    def normalize_blur_kernel(cls, value: int) -> int:
        return value if value % 2 == 1 else value + 1

    @model_validator(mode="before")
    @classmethod
    def normalize_legacy_fields(cls, data: Any) -> Any:
        if not isinstance(data, dict):
            return data

        normalized = dict(data)
        if "blur_kernel_pixels_native" not in normalized:
            for key in (
                "context_window_pixels_native",
                "window_pixels_native",
                "solver_local_window_cells",
            ):
                if key in normalized:
                    normalized["blur_kernel_pixels_native"] = normalized[key]
                    break

        legacy_keys = (
            "window_pixels_native",
            "context_window_pixels_native",
            "solver_local_window_cells",
            "residual_z_threshold",
            "gradient_z_threshold",
            "coherence_window_pixels_native",
            "road_std_z_threshold_native",
            "solver_road_std_z_threshold",
            "road_coherence_threshold_native",
            "solver_road_coherence_threshold",
            "road_std_floor_native",
            "solver_road_std_floor",
            "point_range_z_threshold_native",
            "solver_point_range_z_threshold",
            "point_outlier_fraction_max_native",
            "solver_point_outlier_fraction_max",
            "point_range_floor_native",
            "solver_point_range_floor",
            "outlier_sigma_native",
            "outlier_floor_native",
        )
        for key in legacy_keys:
            normalized.pop(key, None)
        return normalized


class SolverAlgorithmConfig(SIACBaseModel):
    aot_gamma: float = Field(default=10.0, ge=0.0)
    tcwv_gamma: float = Field(default=5.0, ge=0.0)
    alpha: float = -1.6
    smoothness_delta: float = Field(default=0.02, gt=0.0)
    max_iterations: int = Field(default=300, ge=1)
    gtol: float = Field(default=1e-2, gt=0.0)
    ftol: float = Field(default=1e-7, gt=0.0)
    aerosol_resolution: float = Field(default=120.0, gt=0.0)
    grid_search_aot_points: int = Field(default=11, ge=3)
    grid_search_tcwv_points: int = Field(default=11, ge=3)
    fixed_atmospheric_parameter: Literal["none", "aot", "tcwv"] = "none"
    stages: tuple[SolverStageConfig, ...] = Field(default_factory=tuple)
    quadratic_block_size: int = Field(default=1, ge=1)
    quadratic_block_min_valid_fraction: float = Field(default=0.5, ge=0.0, le=1.0)
    water_mask_buffer_pixels: int = Field(default=0, ge=0)
    use_multigrid: bool = True
    min_grid_size: int = Field(default=4, ge=2)
    bounds: SolverBoundsConfig = Field(default_factory=SolverBoundsConfig)
    sharp_transition_filter: SharpTransitionFilterConfig = Field(
        default_factory=SharpTransitionFilterConfig
    )

    @property
    def aot_bounds(self) -> tuple[float, float]:
        return self.bounds.aot

    @property
    def tcwv_bounds(self) -> tuple[float, float]:
        return self.bounds.tcwv


class CloudMaskAlgorithmConfig(SIACBaseModel):
    mode: Literal["auto", "external_file", "user_callable", "none"] = "auto"
    provider: Literal["omnicloudmask"] = "omnicloudmask"
    external_mask_path: Path | None = None
    class_mapping: dict[int, list[int]] | None = None
    unmapped_to_missing: bool = True
    target_resolution_m: float = Field(default=10.0, gt=0.0)
    resolution_policy: Literal["auto", "force"] = "auto"
    allow_upsample_to_target: bool = False
    user_callable: Any | None = None

    @field_validator("external_mask_path", mode="before")
    @classmethod
    def normalize_external_mask_path(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)


class AlgorithmsConfig(SIACBaseModel):
    surface_prior: SurfacePriorAlgorithmConfig = Field(default_factory=SurfacePriorAlgorithmConfig)
    rt: RTAlgorithmConfig = Field(default_factory=RTAlgorithmConfig)
    solver: SolverAlgorithmConfig = Field(default_factory=SolverAlgorithmConfig)
    cloud_mask: CloudMaskAlgorithmConfig = Field(default_factory=CloudMaskAlgorithmConfig)


class ExecutionRuntimeConfig(SIACBaseModel):
    backend: Literal["thread", "dask"] = "thread"
    max_workers: int = Field(default=4, ge=1)
    correction_max_workers: int | None = Field(default=None, ge=1)
    retries: int = Field(default=2, ge=0)
    stage_timeout_s: float | None = Field(default=None, gt=0.0)
    stage_timeouts: dict[str, float] = Field(default_factory=dict)
    dashboard: bool = False
    dashboard_address: str | None = None
    performance_report_path: Path | None = None
    show_progress: bool = False
    profiling_sample_interval_s: float = Field(default=5.0, gt=0.0)
    progress_heartbeat_s: float = Field(default=30.0, gt=0.0)
    max_scatter_points_per_band: int = Field(default=4096, ge=64)
    default_aux_resolution_m: float = Field(default=500.0, gt=0.0)

    @field_validator("performance_report_path", mode="before")
    @classmethod
    def normalize_report_path(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)


class RuntimeConfig(SIACBaseModel):
    log_level: LogLevel = "INFO"
    n_jobs: int = -1
    chunk_size: int = Field(default=2048, gt=0)
    execution: ExecutionRuntimeConfig = Field(default_factory=ExecutionRuntimeConfig)


class OutputDefaultsConfig(SIACBaseModel):
    output_dir: Path | None = None
    format: Literal["geotiff", "cog", "zarr", "netcdf"] = "cog"
    compression: Literal["deflate", "lzw", "zstd", "none"] = "deflate"
    include_rgb: bool = True
    include_uncertainty: bool = True
    include_auxiliary: bool = True
    skip_correction: bool = False
    boa_dtype: Literal["float32", "float64", "uint16"] = "float32"
    boa_scale: float = Field(default=10000.0, gt=0.0)
    boa_nodata: float = 0.0
    reopen_streamed_boa: bool = True

    @field_validator("output_dir", mode="before")
    @classmethod
    def normalize_output_dir(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)


class OutputConfig(SIACBaseModel):
    defaults: OutputDefaultsConfig = Field(default_factory=OutputDefaultsConfig)


class _ConfigShortcutsMixin:
    """Shared shortcut properties for both SystemConfig and ResolvedConfig.

    Both config classes have identical ``algorithms``, ``runtime``, ``output``,
    and ``paths`` sub-structures. This mixin deduplicates the accessor
    properties that navigate into those sub-structures.
    """

    algorithms: AlgorithmsConfig
    runtime: RuntimeConfig
    output: OutputConfig
    paths: PathsConfig | ResolvedPathsConfig

    @property
    def solver(self) -> SolverAlgorithmConfig:
        return self.algorithms.solver

    @property
    def cloud_mask(self) -> CloudMaskAlgorithmConfig:
        return self.algorithms.cloud_mask

    @property
    def execution(self) -> ExecutionRuntimeConfig:
        return self.runtime.execution

    @property
    def global_dem(self) -> str | Path | None:
        return self.paths.dem

    @property
    def global_water(self) -> str | Path | None:
        return self.paths.water_mask

    @property
    def log_level(self) -> LogLevel:
        return self.runtime.log_level

    @property
    def n_jobs(self) -> int:
        return self.runtime.n_jobs

    @property
    def chunk_size(self) -> int:
        return self.runtime.chunk_size


class SystemConfig(_ConfigShortcutsMixin, SIACBaseModel):
    paths: PathsConfig = Field(default_factory=PathsConfig)
    auth: AuthConfig = Field(default_factory=AuthConfig)
    providers: ProvidersConfig = Field(default_factory=ProvidersConfig)
    algorithms: AlgorithmsConfig = Field(default_factory=AlgorithmsConfig)
    runtime: RuntimeConfig = Field(default_factory=RuntimeConfig)
    output: OutputConfig = Field(default_factory=OutputConfig)

    @property
    def atmo_prior(self) -> AtmoProviderConfig:
        return self.providers.atmo

    @property
    def brdf(self) -> BRDFProviderConfig:
        return self.providers.brdf

    @property
    def surface_prior(self) -> SurfacePriorAlgorithmConfig:
        return self.algorithms.surface_prior

    @property
    def rt_model(self) -> RTAlgorithmConfig:
        return self.algorithms.rt

    @property
    def s2_data(self) -> S2ProviderConfig:
        return self.providers.s2

    @property
    def monthly_composites(self) -> MonthlyCompositeProviderConfig:
        return self.providers.monthly_composites


class RunRequest(SIACBaseModel):
    input_path: Path | None = None
    output_path: Path | None = None
    sensor: SensorName | None = None
    aoi: dict[str, Any] | Path | str | tuple[float, float, float, float] | list[float] | None = None
    s2_query: str | Path | None = None

    @field_validator("input_path", "output_path", mode="before")
    @classmethod
    def normalize_paths(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)

    @field_validator("s2_query", mode="before")
    @classmethod
    def normalize_s2_query(cls, value: Any) -> str | Path | None:
        return _coerce_path_or_url(value)


ResolvedCachePathsConfig = CachePathsConfig


class ResolvedPathsConfig(SIACBaseModel):
    dem: str | Path | None = None
    water_mask: str | Path | None = None
    emulator_dir: Path | None = None
    lut_path: str | Path | None = None
    rsrf_root: Path | None = None
    cache_root: Path | None = None
    caches: ResolvedCachePathsConfig = Field(default_factory=ResolvedCachePathsConfig)


class ResolvedCDSEAuthConfig(SIACBaseModel):
    username: str | None = None
    password: str | None = None


class ResolvedCDSAuthConfig(SIACBaseModel):
    api_key: str | None = None


class ResolvedAWSAuthConfig(SIACBaseModel):
    access_key_id: str | None = None
    secret_access_key: str | None = None


class ResolvedEarthdataAuthConfig(SIACBaseModel):
    username: str | None = None
    password: str | None = None


class ResolvedGCSAuthConfig(SIACBaseModel):
    credentials_file: Path | None = None


class ResolvedAuthConfig(SIACBaseModel):
    cdse: ResolvedCDSEAuthConfig = Field(default_factory=ResolvedCDSEAuthConfig)
    cds: ResolvedCDSAuthConfig = Field(default_factory=ResolvedCDSAuthConfig)
    aws: ResolvedAWSAuthConfig = Field(default_factory=ResolvedAWSAuthConfig)
    earthdata: ResolvedEarthdataAuthConfig = Field(default_factory=ResolvedEarthdataAuthConfig)
    gcs: ResolvedGCSAuthConfig = Field(default_factory=ResolvedGCSAuthConfig)


class ResolvedAtmoProviderConfig(AtmoProviderConfig):
    cache_dir: Path | None = None


class ResolvedBRDFProviderConfig(BRDFProviderConfig):
    cache_dir: Path | None = None


class ResolvedS2ProviderConfig(S2ProviderConfig):
    cache_dir: Path | None = None


ResolvedMonthlyCompositeProviderConfig = MonthlyCompositeProviderConfig


class ResolvedProvidersConfig(SIACBaseModel):
    atmo: ResolvedAtmoProviderConfig
    brdf: ResolvedBRDFProviderConfig
    s2: ResolvedS2ProviderConfig
    monthly_composites: ResolvedMonthlyCompositeProviderConfig


class ResolvedSpectralMappingConfig(SpectralMappingAlgorithmConfig):
    pass


class ResolvedSurfacePriorAlgorithmConfig(SurfacePriorAlgorithmConfig):
    spectral_mapping: ResolvedSpectralMappingConfig


class ResolvedRTAlgorithmConfig(RTAlgorithmConfig):
    emulator_dir: Path | None = None
    lut_path: str | Path | None = None


class ResolvedAlgorithmsConfig(SIACBaseModel):
    surface_prior: ResolvedSurfacePriorAlgorithmConfig
    rt: ResolvedRTAlgorithmConfig
    solver: SolverAlgorithmConfig
    cloud_mask: CloudMaskAlgorithmConfig


class ResolvedRunConfig(SIACBaseModel):
    input_path: Path | None = None
    output_path: Path | None = None
    sensor: SensorName = "auto"
    aoi: dict[str, Any] | Path | str | tuple[float, float, float, float] | list[float] | None = None
    s2_query: str | Path | None = None


class ResolvedConfig(_ConfigShortcutsMixin, SIACBaseModel):
    run: ResolvedRunConfig
    paths: ResolvedPathsConfig
    auth: ResolvedAuthConfig
    providers: ResolvedProvidersConfig
    algorithms: ResolvedAlgorithmsConfig
    runtime: RuntimeConfig
    output: OutputConfig

    @property
    def sensor(self) -> SensorName:
        return self.run.sensor

    @property
    def aoi(self) -> dict[str, Any] | Path | str | tuple[float, float, float, float] | list[float] | None:
        return self.run.aoi

    @property
    def atmo_prior(self) -> ResolvedAtmoProviderConfig:
        return self.providers.atmo

    @property
    def brdf(self) -> ResolvedBRDFProviderConfig:
        return self.providers.brdf

    @property
    def s2_data(self) -> ResolvedS2ProviderConfig:
        return self.providers.s2

    @property
    def monthly_composites(self) -> ResolvedMonthlyCompositeProviderConfig:
        return self.providers.monthly_composites

    @property
    def surface_prior(self) -> ResolvedSurfacePriorAlgorithmConfig:
        return self.algorithms.surface_prior

    @property
    def rt_model(self) -> ResolvedRTAlgorithmConfig:
        return self.algorithms.rt
