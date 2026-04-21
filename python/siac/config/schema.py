"""
Typed configuration schema for SIAC.

This module contains only data models. File I/O, environment overlays,
resolution, and snapshotting live in sibling modules.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal, cast
from urllib.parse import urlparse

from pydantic import AliasChoices, BaseModel, Field, field_validator, model_validator

from siac.sixs_outputs import SIXS_DEFAULT_OUTPUT_VARIABLES, SIXS_OUTPUT_VARIABLE_CHOICES

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
SixSMode = Literal["direct", "scene_lut", "auto"]
SixSParallelBackend = Literal["openmp", "worker_libraries"]
SixSBuildProfile = Literal["release", "parity"]
SixSAtmosphericColumnsMode = Literal["input_columns", "profile_default"]
SixSAtmosphericProfile = Literal[
    "no_gas",
    "tropical",
    "midlatitude_summer",
    "midlatitude_winter",
    "subarctic_summer",
    "subarctic_winter",
    "us_standard_62",
    "auto_latitude_date",
    "user_water_ozone",
    "user_profile",
]
SixSAerosolProfile = Literal[
    "none",
    "continental",
    "maritime",
    "urban",
    "desert",
    "biomass_burning",
    "stratospheric",
    "user_mixture",
    "multimodal_log_normal",
    "modified_gamma",
    "junge_power_law",
    "sun_photometer",
    "layered_profile",
    "user_model",
]
SixSBuiltinReflectance = Literal["green_vegetation", "clear_water", "sand", "lake_water"]
SixSSurfaceReflectanceKind = Literal["constant", "built_in", "spectrum"]
SixSSurfaceMode = Literal[
    "homogeneous_lambertian",
    "heterogeneous_lambertian",
    "homogeneous_brdf",
]
SixSBRDFModel = Literal[
    "user_defined",
    "hapke",
    "verstraete",
    "roujean",
    "walthall",
    "minnaert",
    "ocean",
    "iaquinta_pinty",
    "rahman",
    "kuusk",
    "modis",
    "ross_li_maignan",
]
SixSAtmosphericCorrectionMode = Literal[
    "none",
    "lambertian_reflectance",
    "lambertian_radiance",
    "brdf_reflectance",
    "brdf_radiance",
]

DEFAULT_SIXS_SOURCE_URL = "https://salsa.umd.edu/files/6S/6sV2.1.tar"

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


def _coerce_float_tuple(value: Any) -> tuple[float, ...] | None:
    if value is None:
        return None
    if isinstance(value, tuple):
        return tuple(float(item) for item in value)
    if isinstance(value, list):
        return tuple(float(item) for item in value)
    return (float(value),)


def _coerce_int_tuple(value: Any) -> tuple[int, ...] | None:
    if value is None:
        return None
    if isinstance(value, tuple):
        return tuple(int(item) for item in value)
    if isinstance(value, list):
        return tuple(int(item) for item in value)
    return (int(value),)


def _coerce_float_matrix(value: Any) -> tuple[tuple[float, ...], ...] | None:
    if value is None:
        return None
    rows = []
    for row in value:
        rows.append(tuple(float(item) for item in row))
    return tuple(rows)


def _coerce_float_series_with_broadcast(
    value: Any,
    *,
    target_length: int,
) -> tuple[float, ...] | None:
    series = _coerce_float_tuple(value)
    if series is None:
        return None
    if len(series) == 1:
        return tuple(series[0] for _ in range(target_length))
    return series


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


class SixSRadiosondeProfileConfig(SIACBaseModel):
    altitude_km: tuple[float, ...] = Field(validation_alias=AliasChoices("altitude_km", "altitude"))
    pressure_mb: tuple[float, ...] = Field(validation_alias=AliasChoices("pressure_mb", "pressure"))
    temperature_k: tuple[float, ...] = Field(
        validation_alias=AliasChoices("temperature_k", "temperature")
    )
    water_g_m3: tuple[float, ...] = Field(validation_alias=AliasChoices("water_g_m3", "water"))
    ozone_g_m3: tuple[float, ...] = Field(validation_alias=AliasChoices("ozone_g_m3", "ozone"))

    @field_validator(
        "altitude_km",
        "pressure_mb",
        "temperature_k",
        "water_g_m3",
        "ozone_g_m3",
        mode="before",
    )
    @classmethod
    def normalize_series(cls, value: Any) -> tuple[float, ...]:
        series = _coerce_float_tuple(value)
        if series is None:
            raise ValueError("Radiosonde profile fields cannot be null.")
        return series

    @model_validator(mode="after")
    def validate_lengths(self) -> SixSRadiosondeProfileConfig:
        lengths = {
            "altitude_km": len(self.altitude_km),
            "pressure_mb": len(self.pressure_mb),
            "temperature_k": len(self.temperature_k),
            "water_g_m3": len(self.water_g_m3),
            "ozone_g_m3": len(self.ozone_g_m3),
        }
        if len(set(lengths.values())) != 1:
            raise ValueError(
                "All radiosonde profile fields must contain the same number of samples."
            )
        count = next(iter(lengths.values()))
        if count != 34:
            raise ValueError("6S radiosonde profiles must contain exactly 34 levels.")
        return self


class SixSAerosolDistributionComponentConfig(SIACBaseModel):
    rmean: float = Field(gt=0.0)
    sigma: float = Field(gt=0.0)
    percentage_density: float = Field(gt=0.0)
    refr_real: tuple[float, ...]
    refr_imag: tuple[float, ...]

    @field_validator("refr_real", "refr_imag", mode="before")
    @classmethod
    def normalize_refractive_index(cls, value: Any) -> tuple[float, ...]:
        series = _coerce_float_tuple(value)
        if series is None:
            raise ValueError("Refractive-index arrays cannot be null.")
        return series

    @model_validator(mode="after")
    def validate_component(self) -> SixSAerosolDistributionComponentConfig:
        if len(self.refr_real) != 20 or len(self.refr_imag) != 20:
            raise ValueError("Aerosol distribution components require 20 refractive-index values.")
        return self


class SixSAerosolDistributionConfig(SIACBaseModel):
    rmin: float = Field(gt=0.0)
    rmax: float = Field(gt=0.0)
    components: tuple[SixSAerosolDistributionComponentConfig, ...]

    @field_validator("components", mode="before")
    @classmethod
    def normalize_components(cls, value: Any) -> tuple[SixSAerosolDistributionComponentConfig, ...]:
        if value is None:
            raise ValueError("Aerosol distribution components cannot be null.")
        return tuple(value)

    @model_validator(mode="after")
    def validate_distribution(self) -> SixSAerosolDistributionConfig:
        if self.rmax <= self.rmin:
            raise ValueError("Aerosol distribution rmax must be greater than rmin.")
        if not (1 <= len(self.components) <= 4):
            raise ValueError("6S aerosol distributions require between 1 and 4 components.")
        return self


class SixSSunPhotometerAerosolConfig(SIACBaseModel):
    radii_um: tuple[float, ...] = Field(validation_alias=AliasChoices("radii_um", "r"))
    dv_dlogr: tuple[float, ...] = Field(validation_alias=AliasChoices("dv_dlogr", "dvdlogr"))
    refr_real: tuple[float, ...]
    refr_imag: tuple[float, ...]

    @field_validator("radii_um", "dv_dlogr", mode="before")
    @classmethod
    def normalize_series(cls, value: Any) -> tuple[float, ...]:
        series = _coerce_float_tuple(value)
        if series is None:
            raise ValueError("Sun-photometer aerosol fields cannot be null.")
        return series

    @field_validator("refr_real", "refr_imag", mode="before")
    @classmethod
    def normalize_refractive_index(cls, value: Any) -> tuple[float, ...]:
        series = _coerce_float_series_with_broadcast(value, target_length=20)
        if series is None:
            raise ValueError("Sun-photometer aerosol fields cannot be null.")
        return series

    @model_validator(mode="after")
    def validate_lengths(self) -> SixSSunPhotometerAerosolConfig:
        if len(self.radii_um) != len(self.dv_dlogr):
            raise ValueError("Sun-photometer radii and dv_dlogr arrays must have the same length.")
        if not (1 <= len(self.radii_um) <= 50):
            raise ValueError("6S sun-photometer aerosol inputs support between 1 and 50 radius samples.")
        if len(self.refr_real) != 20 or len(self.refr_imag) != 20:
            raise ValueError("Sun-photometer aerosol refractive indices require 20 wavelengths.")
        return self


class SixSAerosolLayerConfig(SIACBaseModel):
    height_km: float = Field(gt=0.0, validation_alias=AliasChoices("height_km", "height"))
    optical_thickness: float = Field(ge=0.0)
    aerosol_type: Literal[
        "none",
        "continental",
        "maritime",
        "urban",
        "desert",
        "biomass_burning",
        "stratospheric",
    ] = "continental"


class SixSSpectralReflectanceConfig(SIACBaseModel):
    kind: SixSSurfaceReflectanceKind = "constant"
    constant: float | None = None
    built_in: SixSBuiltinReflectance | None = None
    values: tuple[float, ...] | None = None
    wavelengths_um: tuple[float, ...] | None = Field(
        default=None,
        validation_alias=AliasChoices("wavelengths_um", "wavelengths"),
    )

    @field_validator("values", "wavelengths_um", mode="before")
    @classmethod
    def normalize_series(cls, value: Any) -> tuple[float, ...] | None:
        return _coerce_float_tuple(value)

    @model_validator(mode="after")
    def validate_payload(self) -> SixSSpectralReflectanceConfig:
        if self.kind == "constant":
            if self.constant is None:
                raise ValueError("Constant surface reflectance requires `constant`.")
        elif self.kind == "built_in":
            if self.built_in is None:
                raise ValueError("Built-in surface reflectance requires `built_in`.")
        elif self.kind == "spectrum":
            if self.values is None:
                raise ValueError("Spectrum surface reflectance requires `values`.")
            if self.wavelengths_um is not None and len(self.wavelengths_um) != len(self.values):
                raise ValueError("Surface spectrum wavelengths and values must have the same length.")
        return self


class SixSBRDFConfig(SIACBaseModel):
    model: SixSBRDFModel
    parameters: dict[str, float] = Field(default_factory=dict)
    options: dict[str, int] = Field(default_factory=dict)
    table_solar_zenith: tuple[tuple[float, ...], ...] | None = None
    table_view_zenith: tuple[tuple[float, ...], ...] | None = None
    spherical_albedo: float | None = None
    directional_reflectance: float | None = None

    @field_validator("table_solar_zenith", "table_view_zenith", mode="before")
    @classmethod
    def normalize_tables(cls, value: Any) -> tuple[tuple[float, ...], ...] | None:
        return _coerce_float_matrix(value)

    @model_validator(mode="after")
    def validate_brdf(self) -> SixSBRDFConfig:
        if self.model == "user_defined":
            if self.table_solar_zenith is None or self.table_view_zenith is None:
                raise ValueError("User-defined BRDF requires both solar and view reflectance tables.")
            if len(self.table_solar_zenith) != 13 or len(self.table_view_zenith) != 13:
                raise ValueError("User-defined BRDF tables must have 13 azimuth rows.")
            if any(len(row) != 10 for row in self.table_solar_zenith + self.table_view_zenith):
                raise ValueError("User-defined BRDF tables must have 10 zenith samples per row.")
            if self.spherical_albedo is None or self.directional_reflectance is None:
                raise ValueError(
                    "User-defined BRDF requires `spherical_albedo` and `directional_reflectance`."
                )
        return self


class SixSSurfaceConfig(SIACBaseModel):
    mode: SixSSurfaceMode = "homogeneous_lambertian"
    target: SixSSpectralReflectanceConfig = Field(
        default_factory=lambda: SixSSpectralReflectanceConfig(kind="constant", constant=0.0)
    )
    environment: SixSSpectralReflectanceConfig | None = None
    radius_km: float = Field(default=1.0, gt=0.0, validation_alias=AliasChoices("radius_km", "radius"))
    brdf: SixSBRDFConfig | None = None

    @model_validator(mode="after")
    def validate_surface(self) -> SixSSurfaceConfig:
        if self.mode == "heterogeneous_lambertian" and self.environment is None:
            raise ValueError("Heterogeneous Lambertian surfaces require an `environment` reflectance.")
        if self.mode == "homogeneous_brdf" and self.brdf is None:
            raise ValueError("BRDF surfaces require a `brdf` configuration.")
        return self


class SixSAtmosphericCorrectionConfig(SIACBaseModel):
    mode: SixSAtmosphericCorrectionMode = "lambertian_reflectance"
    value: float | None = None

    @model_validator(mode="after")
    def validate_correction(self) -> SixSAtmosphericCorrectionConfig:
        if self.mode != "none" and self.value is None:
            raise ValueError("Atmospheric correction modes other than `none` require a `value`.")
        return self


class SixSAlgorithmConfig(SIACBaseModel):
    source_url: str = DEFAULT_SIXS_SOURCE_URL
    source_dir: Path | None = None
    build_dir: Path | None = None
    module_path: Path | None = None
    library_path: Path | None = None
    auto_build: bool = True
    compiler: str = "gfortran"
    build_profile: SixSBuildProfile = "release"
    mode: SixSMode = "auto"
    parallel_backend: SixSParallelBackend = "openmp"
    native_threads: int | None = Field(default=None, ge=1)
    worker_libraries: int | None = Field(default=None, ge=1)
    chunk_size: int = Field(default=4096, ge=1)
    scene_lut_min_pixels: int = Field(default=512, ge=1)
    scene_lut_max_nodes_per_axis: int = Field(default=4, ge=1)
    scene_lut_max_cases: int = Field(default=4096, ge=1)
    scene_lut_required_speedup: float = Field(default=1.5, gt=1.0)
    month: int = Field(default=1, ge=1, le=12)
    day: int = Field(default=1, ge=1, le=31)
    atmospheric_profile_latitude: float | None = Field(default=None, ge=-90.0, le=90.0)
    atmospheric_columns_mode: SixSAtmosphericColumnsMode = "input_columns"
    reference_reflectance: float = Field(default=0.1, gt=0.0, le=1.0)
    atmospheric_profile: SixSAtmosphericProfile = "user_water_ozone"
    radiosonde_profile: SixSRadiosondeProfileConfig | None = None
    aerosol_profile: SixSAerosolProfile = "continental"
    aerosol_mixture: tuple[float, float, float, float] | None = None
    aerosol_distribution: SixSAerosolDistributionConfig | None = None
    sun_photometer_aerosol: SixSSunPhotometerAerosolConfig | None = None
    aerosol_layer_profile: tuple[SixSAerosolLayerConfig, ...] | None = None
    aerosol_model_path: Path | None = Field(
        default=None,
        validation_alias=AliasChoices("aerosol_model_path", "mie_file_path"),
    )
    surface: SixSSurfaceConfig = Field(default_factory=SixSSurfaceConfig)
    atmospheric_correction: SixSAtmosphericCorrectionConfig = Field(
        default_factory=lambda: SixSAtmosphericCorrectionConfig(
            mode="lambertian_reflectance",
            value=0.1,
        )
    )
    output_variables: tuple[str, ...] = SIXS_DEFAULT_OUTPUT_VARIABLES

    @field_validator("source_dir", "build_dir", "module_path", "library_path", "aerosol_model_path", mode="before")
    @classmethod
    def normalize_local_paths(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)

    @field_validator("atmospheric_profile", mode="before")
    @classmethod
    def normalize_atmospheric_profile(cls, value: Any) -> Any:
        if value is None:
            return value
        key = str(value).strip().lower()
        aliases = {
            "from_latitude_and_date": "auto_latitude_date",
            "latitude_date": "auto_latitude_date",
            "user_water_and_ozone": "user_water_ozone",
            "radiosonde": "user_profile",
            "radiosonde_profile": "user_profile",
        }
        return aliases.get(key, value)

    @field_validator("atmospheric_columns_mode", mode="before")
    @classmethod
    def normalize_atmospheric_columns_mode(cls, value: Any) -> Any:
        if value is None:
            return value
        key = str(value).strip().lower()
        aliases = {
            "input": "input_columns",
            "inputs": "input_columns",
            "scene_inputs": "input_columns",
            "scene_input_columns": "input_columns",
            "profile": "profile_default",
            "profile_defaults": "profile_default",
            "profile_columns": "profile_default",
        }
        return aliases.get(key, value)

    @field_validator("aerosol_profile", mode="before")
    @classmethod
    def normalize_aerosol_profile(cls, value: Any) -> Any:
        if value is None:
            return value
        key = str(value).strip().lower()
        aliases = {
            "user": "user_mixture",
            "from_mie_file": "user_model",
            "mie_file": "user_model",
            "user_profile": "layered_profile",
            "multimodal_lognormal": "multimodal_log_normal",
            "multimodal_log_normal_distribution": "multimodal_log_normal",
            "modified_gamma_distribution": "modified_gamma",
            "junge_power_law_distribution": "junge_power_law",
            "sun_photometer_distribution": "sun_photometer",
        }
        return aliases.get(key, value)

    @field_validator("aerosol_mixture", mode="before")
    @classmethod
    def normalize_aerosol_mixture(cls, value: Any) -> tuple[float, float, float, float] | None:
        if value is None:
            return None
        if isinstance(value, dict):
            mapping = {str(key).strip().lower(): float(component) for key, component in value.items()}
            return (
                mapping.get("dust", 0.0),
                mapping.get("water", 0.0),
                mapping.get("oceanic", 0.0),
                mapping.get("soot", 0.0),
            )
        series = _coerce_float_tuple(value)
        if series is None:
            return None
        if len(series) != 4:
            raise ValueError("sixs.aerosol_mixture must contain four components: dust, water, oceanic, soot.")
        return (
            float(series[0]),
            float(series[1]),
            float(series[2]),
            float(series[3]),
        )

    @field_validator("output_variables", mode="before")
    @classmethod
    def normalize_output_variables(cls, value: Any) -> tuple[str, ...]:
        if value is None:
            return cast("tuple[str, ...]", SIXS_DEFAULT_OUTPUT_VARIABLES)
        items = [value] if isinstance(value, str) else list(value)
        normalized: list[str] = []
        for item in items:
            text = str(item).strip()
            if not text:
                continue
            if text not in SIXS_OUTPUT_VARIABLE_CHOICES:
                raise ValueError(
                    "Unknown 6S output variable "
                    f"{text!r}. Expected one of {', '.join(SIXS_OUTPUT_VARIABLE_CHOICES)}."
                )
            if text not in normalized:
                normalized.append(text)
        return tuple(normalized) or SIXS_DEFAULT_OUTPUT_VARIABLES

    @model_validator(mode="after")
    def validate_custom_aerosol(self) -> SixSAlgorithmConfig:
        if self.atmospheric_profile == "auto_latitude_date" and self.atmospheric_profile_latitude is None:
            raise ValueError(
                "sixs.atmospheric_profile_latitude must be provided when atmospheric_profile='auto_latitude_date'."
            )
        if self.atmospheric_profile == "user_profile" and self.radiosonde_profile is None:
            raise ValueError(
                "sixs.radiosonde_profile must be provided when atmospheric_profile='user_profile'."
            )
        if self.atmospheric_profile != "user_profile" and self.radiosonde_profile is not None:
            raise ValueError(
                "sixs.radiosonde_profile is only valid when atmospheric_profile='user_profile'."
            )
        if self.aerosol_profile == "user_mixture" and self.aerosol_mixture is None:
            raise ValueError(
                "sixs.aerosol_mixture must be provided when aerosol_profile='user_mixture'."
            )
        if self.aerosol_profile == "user_mixture" and self.aerosol_mixture is not None:
            mixture_total = sum(float(value) for value in self.aerosol_mixture)
            if abs(mixture_total - 1.0) > 1e-6:
                raise ValueError(
                    "sixs.aerosol_mixture must sum to 1.0 for user_mixture aerosol profiles."
                )
        if self.aerosol_profile != "user_mixture" and self.aerosol_mixture is not None:
            raise ValueError(
                "sixs.aerosol_mixture is only valid when aerosol_profile='user_mixture'."
            )
        if self.aerosol_profile in {"multimodal_log_normal", "modified_gamma", "junge_power_law"}:
            if self.aerosol_distribution is None:
                raise ValueError(
                    "sixs.aerosol_distribution must be provided for aerosol distribution profiles."
                )
            if (
                self.aerosol_profile in {"modified_gamma", "junge_power_law"}
                and len(self.aerosol_distribution.components) != 1
            ):
                raise ValueError(
                    "6S modified_gamma and junge_power_law profiles require exactly one component."
                )
        elif self.aerosol_distribution is not None:
            raise ValueError(
                "sixs.aerosol_distribution is only valid for multimodal_log_normal, modified_gamma, or junge_power_law."
            )
        if self.aerosol_profile == "sun_photometer" and self.sun_photometer_aerosol is None:
            raise ValueError(
                "sixs.sun_photometer_aerosol must be provided when aerosol_profile='sun_photometer'."
            )
        if self.aerosol_profile != "sun_photometer" and self.sun_photometer_aerosol is not None:
            raise ValueError(
                "sixs.sun_photometer_aerosol is only valid when aerosol_profile='sun_photometer'."
            )
        if self.aerosol_profile == "layered_profile" and not self.aerosol_layer_profile:
            raise ValueError(
                "sixs.aerosol_layer_profile must be provided when aerosol_profile='layered_profile'."
            )
        if self.aerosol_profile != "layered_profile" and self.aerosol_layer_profile is not None:
            raise ValueError(
                "sixs.aerosol_layer_profile is only valid when aerosol_profile='layered_profile'."
            )
        if self.aerosol_profile == "user_model" and self.aerosol_model_path is None:
            raise ValueError(
                "sixs.aerosol_model_path must be provided when aerosol_profile='user_model'."
            )
        if self.aerosol_profile != "user_model" and self.aerosol_model_path is not None:
            raise ValueError(
                "sixs.aerosol_model_path is only valid when aerosol_profile='user_model'."
            )
        if self.aerosol_profile == "layered_profile" and self.aerosol_layer_profile is not None:
            aerosol_types = {layer.aerosol_type for layer in self.aerosol_layer_profile}
            if len(aerosol_types) > 1:
                raise ValueError(
                    "6S layered aerosol profiles require the same aerosol_type for every layer."
                )
        if self.surface.mode != "homogeneous_brdf" and self.atmospheric_correction.mode.startswith("brdf_"):
            raise ValueError(
                "BRDF atmospheric correction modes require surface.mode='homogeneous_brdf'."
            )
        if (
            self.atmospheric_correction.mode == "lambertian_reflectance"
            and abs((self.atmospheric_correction.value or -1.0) - 0.1) < 1e-9
            and abs(self.reference_reflectance - 0.1) > 1e-9
        ):
            self.atmospheric_correction = SixSAtmosphericCorrectionConfig(
                mode="lambertian_reflectance",
                value=self.reference_reflectance,
            )
        if self.parallel_backend == "worker_libraries" and self.worker_libraries is not None:
            if self.worker_libraries < 1:
                raise ValueError("sixs.worker_libraries must be >= 1.")
            if self.chunk_size < 1:
                raise ValueError("sixs.chunk_size must be >= 1 for worker_libraries.")
        if self.scene_lut_max_cases < self.scene_lut_max_nodes_per_axis:
            raise ValueError(
                "sixs.scene_lut_max_cases must be >= sixs.scene_lut_max_nodes_per_axis."
            )
        return self


class RTAlgorithmConfig(SIACBaseModel):
    backend: Literal["emulator", "lut", "py6s", "sixs"] = "emulator"
    lut_interpolation: Literal["linear", "nearest", "cubic"] = "linear"
    lut_storage_options: dict[str, Any] = Field(default_factory=dict)
    py6s_aero_profile: str = "Continental"
    sixs: SixSAlgorithmConfig = Field(default_factory=SixSAlgorithmConfig)
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
