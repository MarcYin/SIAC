"""Algorithm configuration models."""

from __future__ import annotations

from pathlib import Path  # noqa: TC003
from typing import Any, cast

from pydantic import Field, field_validator, model_validator

from siac.config._base import SIACBaseModel
from siac.config.types import (
    DEFAULT_SIXS_SOURCE_URL,
    AtmosphericParameterName,
    CloudMaskMode,
    CloudMaskProvider,
    FixedAtmosphericParameter,
    LUTInterpolationMethod,
    MonthlyDatabaseResolutionPolicy,
    ResolutionPolicy,
    RTAerosolProfile,
    RTBackend,
    SixSAerosolLayerType,
    SixSAtmosphericColumnsMode,
    SixSAtmosphericCorrectionMode,
    SixSAtmosphericProfile,
    SixSBRDFModel,
    SixSBuildProfile,
    SixSBuiltinReflectance,
    SixSMode,
    SixSParallelBackend,
    SixSSurfaceMode,
    SixSSurfaceReflectanceKind,
    SolverStageInitialState,
    SurfacePriorMethod,
    _coerce_float_matrix,
    _coerce_float_series_with_broadcast,
    _coerce_float_tuple,
    _coerce_pathlike,
    _normalize_aerosol_mixture_payload,
)
from siac.sixs_outputs import SIXS_DEFAULT_OUTPUT_VARIABLES, SIXS_OUTPUT_VARIABLE_CHOICES


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
    method: SurfacePriorMethod = SurfacePriorMethod.KERNEL_MODEL
    monthly_database_resolution_policy: MonthlyDatabaseResolutionPolicy = (
        MonthlyDatabaseResolutionPolicy.PROVIDER_OR_COARSER
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
    altitude_km: tuple[float, ...]
    pressure_mb: tuple[float, ...]
    temperature_k: tuple[float, ...]
    water_g_m3: tuple[float, ...]
    ozone_g_m3: tuple[float, ...]

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
    radii_um: tuple[float, ...]
    dv_dlogr: tuple[float, ...]
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
            raise ValueError(
                "6S sun-photometer aerosol inputs support between 1 and 50 radius samples."
            )
        if len(self.refr_real) != 20 or len(self.refr_imag) != 20:
            raise ValueError("Sun-photometer aerosol refractive indices require 20 wavelengths.")
        return self


class SixSAerosolLayerConfig(SIACBaseModel):
    height_km: float = Field(gt=0.0)
    optical_thickness: float = Field(ge=0.0)
    aerosol_type: SixSAerosolLayerType = SixSAerosolLayerType.CONTINENTAL


class SixSSpectralReflectanceConfig(SIACBaseModel):
    kind: SixSSurfaceReflectanceKind = SixSSurfaceReflectanceKind.CONSTANT
    constant: float | None = None
    built_in: SixSBuiltinReflectance | None = None
    values: tuple[float, ...] | None = None
    wavelengths_um: tuple[float, ...] | None = None

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
                raise ValueError(
                    "Surface spectrum wavelengths and values must have the same length."
                )
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
                raise ValueError(
                    "User-defined BRDF requires both solar and view reflectance tables."
                )
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
    mode: SixSSurfaceMode = SixSSurfaceMode.HOMOGENEOUS_LAMBERTIAN
    target: SixSSpectralReflectanceConfig = Field(
        default_factory=lambda: SixSSpectralReflectanceConfig(
            kind=SixSSurfaceReflectanceKind.CONSTANT,
            constant=0.0,
        )
    )
    environment: SixSSpectralReflectanceConfig | None = None
    radius_km: float = Field(default=1.0, gt=0.0)
    brdf: SixSBRDFConfig | None = None

    @model_validator(mode="after")
    def validate_surface(self) -> SixSSurfaceConfig:
        if self.mode == "heterogeneous_lambertian" and self.environment is None:
            raise ValueError(
                "Heterogeneous Lambertian surfaces require an `environment` reflectance."
            )
        if self.mode == "homogeneous_brdf" and self.brdf is None:
            raise ValueError("BRDF surfaces require a `brdf` configuration.")
        return self


class SixSAtmosphericCorrectionConfig(SIACBaseModel):
    mode: SixSAtmosphericCorrectionMode = SixSAtmosphericCorrectionMode.LAMBERTIAN_REFLECTANCE
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
    auto_build: bool = True
    compiler: str = "gfortran"
    build_profile: SixSBuildProfile = SixSBuildProfile.RELEASE
    mode: SixSMode = SixSMode.AUTO
    parallel_backend: SixSParallelBackend = SixSParallelBackend.OPENMP
    native_threads: int | None = Field(default=None, ge=1)
    worker_libraries: int | None = Field(default=None, ge=1)
    chunk_size: int = Field(default=4096, ge=1)
    scene_lut_min_pixels: int = Field(default=512, ge=1)
    scene_lut_max_nodes_per_axis: int = Field(default=4, ge=1)
    scene_lut_max_cases: int = Field(default=4096, ge=1)
    scene_lut_required_speedup: float = Field(default=1.5, gt=1.0)
    #: Whether the native 6S extension was built with polarized radiative
    #: transfer enabled (``ipol=1`` in the Fortran source). The polarized
    #: branch computes Stokes Q and U via ``ospol_``/``kernelpol_``, which
    #: dominates the runtime at off-nadir geometry. SIAC's downstream
    #: pipeline uses only the scalar (I) component, so polarization
    #: contributes only via second-order corrections to the scalar
    #: result — typically <1% in BOA reflectance at typical S2 geometry.
    #:
    #: This is a **build-time** flag. Flipping it requires rebuilding the
    #: native extension (``pixi run -e rt6s build-6s-native`` after
    #: setting ``SIAC_SIXS_COMPUTE_POLARIZATION=1`` or ``0`` in the
    #: environment). Setting this field on the config alone is not
    #: enough — the build picks up the value via the env var at
    #: source-patching time.
    #:
    #: Default ``False`` (since wave 16) for ~3-5x faster runtime.
    compute_polarization: bool = False
    month: int = Field(default=1, ge=1, le=12)
    day: int = Field(default=1, ge=1, le=31)
    output_variables: tuple[str, ...] = SIXS_DEFAULT_OUTPUT_VARIABLES

    @field_validator("source_dir", "build_dir", "module_path", mode="before")
    @classmethod
    def normalize_local_paths(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)

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


class RTAtmosphereSetupConfig(SIACBaseModel):
    profile: SixSAtmosphericProfile | None = None
    columns_mode: SixSAtmosphericColumnsMode | None = None
    profile_latitude: float | None = Field(default=None, ge=-90.0, le=90.0)
    radiosonde_profile: SixSRadiosondeProfileConfig | None = None

    @model_validator(mode="after")
    def validate_profile_constraints(self) -> RTAtmosphereSetupConfig:
        if self.profile == "auto_latitude_date" and self.profile_latitude is None:
            raise ValueError(
                "rt.setup.atmosphere.profile_latitude is required when profile='auto_latitude_date'."
            )
        if self.profile == "user_profile" and self.radiosonde_profile is None:
            raise ValueError(
                "rt.setup.atmosphere.radiosonde_profile is required when profile='user_profile'."
            )
        if self.profile not in {None, "user_profile"} and self.radiosonde_profile is not None:
            raise ValueError(
                "rt.setup.atmosphere.radiosonde_profile is only valid when profile='user_profile'."
            )
        return self


class RTAerosolSetupConfig(SIACBaseModel):
    profile: RTAerosolProfile | None = None
    mixture: tuple[float, float, float, float] | None = None
    distribution: SixSAerosolDistributionConfig | None = None
    sun_photometer_aerosol: SixSSunPhotometerAerosolConfig | None = None
    layer_profile: tuple[SixSAerosolLayerConfig, ...] | None = None
    model_path: Path | None = None

    @field_validator("mixture", mode="before")
    @classmethod
    def normalize_mixture(cls, value: Any) -> tuple[float, float, float, float] | None:
        return _normalize_aerosol_mixture_payload(value)

    @field_validator("model_path", mode="before")
    @classmethod
    def normalize_model_path(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)

    @model_validator(mode="after")
    def validate_aerosol_constraints(self) -> RTAerosolSetupConfig:
        if self.profile == "user_mixture" and self.mixture is None:
            raise ValueError(
                "rt.setup.aerosol.mixture must be provided when profile='user_mixture'."
            )
        if self.profile == "user_mixture" and self.mixture is not None:
            mixture_total = sum(float(value) for value in self.mixture)
            if abs(mixture_total - 1.0) > 1e-6:
                raise ValueError(
                    "rt.setup.aerosol.mixture must sum to 1.0 for user_mixture aerosol profiles."
                )
        if self.profile != "user_mixture" and self.mixture is not None:
            raise ValueError("rt.setup.aerosol.mixture is only valid when profile='user_mixture'.")
        if self.profile in {"multimodal_log_normal", "modified_gamma", "junge_power_law"}:
            if self.distribution is None:
                raise ValueError(
                    "rt.setup.aerosol.distribution must be provided for aerosol distribution profiles."
                )
        elif self.distribution is not None:
            raise ValueError(
                "rt.setup.aerosol.distribution is only valid for multimodal_log_normal, modified_gamma, or junge_power_law."
            )
        if self.profile == "sun_photometer" and self.sun_photometer_aerosol is None:
            raise ValueError(
                "rt.setup.aerosol.sun_photometer_aerosol must be provided when profile='sun_photometer'."
            )
        if self.profile != "sun_photometer" and self.sun_photometer_aerosol is not None:
            raise ValueError(
                "rt.setup.aerosol.sun_photometer_aerosol is only valid when profile='sun_photometer'."
            )
        if self.profile == "layered_profile" and not self.layer_profile:
            raise ValueError(
                "rt.setup.aerosol.layer_profile must be provided when profile='layered_profile'."
            )
        if self.profile != "layered_profile" and self.layer_profile is not None:
            raise ValueError(
                "rt.setup.aerosol.layer_profile is only valid when profile='layered_profile'."
            )
        if self.profile == "user_model" and self.model_path is None:
            raise ValueError(
                "rt.setup.aerosol.model_path must be provided when profile='user_model'."
            )
        if self.profile != "user_model" and self.model_path is not None:
            raise ValueError("rt.setup.aerosol.model_path is only valid when profile='user_model'.")
        return self


class RTSurfaceSetupConfig(SIACBaseModel):
    mode: SixSSurfaceMode | None = None
    target: SixSSpectralReflectanceConfig | None = None
    environment: SixSSpectralReflectanceConfig | None = None
    radius_km: float | None = Field(default=None, gt=0.0)
    brdf: SixSBRDFConfig | None = None

    @model_validator(mode="after")
    def validate_surface(self) -> RTSurfaceSetupConfig:
        if self.mode == "heterogeneous_lambertian" and self.environment is None:
            raise ValueError(
                "Heterogeneous Lambertian surfaces require an `environment` reflectance."
            )
        if self.mode == "homogeneous_brdf" and self.brdf is None:
            raise ValueError("BRDF surfaces require a `brdf` configuration.")
        return self


class RTAtmosphericCorrectionSetupConfig(SIACBaseModel):
    mode: SixSAtmosphericCorrectionMode | None = None
    value: float | None = None

    @model_validator(mode="after")
    def validate_correction(self) -> RTAtmosphericCorrectionSetupConfig:
        if self.mode not in {None, "none"} and self.value is None:
            raise ValueError("Atmospheric correction modes other than `none` require a `value`.")
        return self


class RTSetupConfig(SIACBaseModel):
    atmosphere: RTAtmosphereSetupConfig | None = None
    aerosol: RTAerosolSetupConfig | None = None
    surface: RTSurfaceSetupConfig | None = None
    atmospheric_correction: RTAtmosphericCorrectionSetupConfig | None = None
    reference_reflectance: float | None = Field(default=None, gt=0.0, le=1.0)

    @model_validator(mode="after")
    def validate_setup_consistency(self) -> RTSetupConfig:
        if (
            self.surface is not None
            and self.atmospheric_correction is not None
            and self.surface.mode != "homogeneous_brdf"
            and self.atmospheric_correction.mode is not None
            and self.atmospheric_correction.mode.startswith("brdf_")
        ):
            raise ValueError(
                "BRDF atmospheric correction modes require rt.setup.surface.mode='homogeneous_brdf'."
            )
        return self

    def has_overrides(self) -> bool:
        return any(
            value is not None
            for value in (
                self.atmosphere,
                self.aerosol,
                self.surface,
                self.atmospheric_correction,
                self.reference_reflectance,
            )
        )


class RTAlgorithmConfig(SIACBaseModel):
    backend: RTBackend = RTBackend.EMULATOR
    lut_interpolation: LUTInterpolationMethod = LUTInterpolationMethod.LINEAR
    lut_storage_options: dict[str, Any] = Field(default_factory=dict)
    setup: RTSetupConfig = Field(default_factory=RTSetupConfig)
    sixs: SixSAlgorithmConfig = Field(default_factory=SixSAlgorithmConfig)
    fallback_to_lut: bool = True


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
    solve: tuple[AtmosphericParameterName, ...] = (
        AtmosphericParameterName.AOT,
        AtmosphericParameterName.TCWV,
    )
    fixed: tuple[AtmosphericParameterName, ...] = (AtmosphericParameterName.TCO3,)
    bands: tuple[str, ...] | None = None
    initial_state: SolverStageInitialState = SolverStageInitialState.PREVIOUS

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
                f"Solver stage {self.name!r} cannot both solve and fix {', '.join(sorted(overlap))}"
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
    fixed_atmospheric_parameter: FixedAtmosphericParameter = FixedAtmosphericParameter.NONE
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
    mode: CloudMaskMode = CloudMaskMode.AUTO
    provider: CloudMaskProvider = CloudMaskProvider.OMNICLOUDMASK
    external_mask_path: Path | None = None
    class_mapping: dict[int, list[int]] | None = None
    unmapped_to_missing: bool = True
    target_resolution_m: float = Field(default=10.0, gt=0.0)
    resolution_policy: ResolutionPolicy = ResolutionPolicy.AUTO
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
