"""Algorithm configuration models."""

from __future__ import annotations

from pathlib import Path  # noqa: TC003
from typing import Any, Literal, cast

from pydantic import Field, field_validator, model_validator

from siac.config._base import SIACBaseModel
from siac.config.types import (
    DEFAULT_LIBRADTRAN_OPTPROP_URL,
    DEFAULT_LIBRADTRAN_REPTRAN_URL,
    DEFAULT_LIBRADTRAN_SOURCE_URL,
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
    SolverMethod,
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
    #: kNN spectral-lookup key for the Route-B monthly database. ``"query"``
    #: (default) keys on ``[corrected NIR/SWIR | temporal median of NIR/SWIR]``.
    #: ``"all"`` additionally appends the per-pixel temporal median of the
    #: *visible* bands — a full-spectrum location fingerprint — so neighbours are
    #: matched on each pixel's multi-year climatology as well as its current
    #: SWIR signal. Default preserves the original behaviour (and goldens).
    monthly_database_median_key: Literal["query", "all"] = "query"
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
    #: Persistent on-disk cache of the native scene-LUT grid batch (the Fortran
    #: RT compute), keyed by the exact batch inputs + compiled-module identity. A
    #: re-run of the same scene reuses the grid coefficients instead of
    #: recomputing; the cheap per-pixel interpolation is never cached.
    run_cache_enabled: bool = True
    #: Cache directory; defaults to ``~/.cache/siac/rt6s/run_cache``.
    run_cache_dir: Path | None = None
    mode: SixSMode = SixSMode.AUTO
    #: Parallelism mode for 6S native calls.
    #:
    #: - ``"openmp"`` shares one library copy and uses OpenMP threads
    #:   within each batch — efficient for a single large batch but
    #:   leaves cores idle when several bands need to run.
    #: - ``"worker_libraries"`` loads N isolated library copies and
    #:   dispatches different batches (e.g. the per-band joint-LUT
    #:   batches) concurrently. Each worker uses
    #:   ``max(1, native_threads // worker_count)`` OpenMP threads so the
    #:   total stays inside the core budget. This is the default since
    #:   wave 18 because the joint-LUT band loop is the largest remaining
    #:   chunk of wall-clock that benefits from band-level parallelism.
    parallel_backend: SixSParallelBackend = SixSParallelBackend.WORKER_LIBRARIES
    native_threads: int | None = Field(default=None, ge=1)
    worker_libraries: int | None = Field(default=None, ge=1)
    chunk_size: int = Field(default=4096, ge=1)
    scene_lut_min_pixels: int = Field(default=512, ge=1)
    scene_lut_max_nodes_per_axis: int = Field(default=4, ge=1)
    scene_lut_max_cases: int = Field(default=4096, ge=1)
    scene_lut_required_speedup: float = Field(default=1.5, gt=1.0)
    #: Joint (aot × tcwv × geometry) LUT used by the solver's block-grid search.
    #: When enabled, the solver builds **one** large LUT spanning the entire
    #: grid-search range across (aot, tcwv) — plus the per-pixel geometric
    #: dimensions — before invoking the inner block-grid-search kernel. Each
    #: (aot_val, tcwv_val) candidate is then served by interpolation rather
    #: than a fresh 6S batch. This eliminates the ~N_aot × N_tcwv redundant
    #: 6S evaluations the previous per-candidate scene-LUT path incurred.
    #: ``joint_grid_search_lut_enabled`` toggles the optimization; the other
    #: two fields control the LUT size budget (analogous to the per-candidate
    #: scene-LUT fields). The default ``joint_grid_search_lut_max_cases`` is
    #: sized to fully preserve the per-candidate scene-LUT's geometric
    #: resolution (4 nodes per axis × 6 geometric axes = 4096 cases) on top
    #: of the typical 11×11 aot/tcwv grid: 121 × 4^6 ≈ 500 K. This keeps
    #: joint-LUT outputs numerically equivalent to the per-candidate
    #: scene-LUT path at the LUT nodes — only the geometric axes use linear
    #: interpolation, exactly as the per-candidate path does.
    joint_grid_search_lut_enabled: bool = True
    joint_grid_search_lut_max_nodes_per_axis: int = Field(default=4, ge=1)
    joint_grid_search_lut_max_cases: int = Field(default=524288, ge=1)
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

    @field_validator("source_dir", "build_dir", "module_path", "run_cache_dir", mode="before")
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
        if self.joint_grid_search_lut_max_cases < self.joint_grid_search_lut_max_nodes_per_axis:
            raise ValueError(
                "sixs.joint_grid_search_lut_max_cases must be >= "
                "sixs.joint_grid_search_lut_max_nodes_per_axis."
            )
        return self


class LibRadtranAlgorithmConfig(SIACBaseModel):
    """Configuration for the libRadtran (``uvspec``) RT backend.

    The backend compiles libRadtran from source (autotools) and drives the
    ``uvspec`` binary in batch over a scene-scoped grid, then interpolates
    per-pixel — mirroring 6S's ``scene_lut`` mode. It reproduces the libRadtran
    preset the public remote LUT encodes (US Standard atmosphere,
    ``continental_average`` aerosol, homogeneous Lambertian surface) via the
    two-albedo method, yielding the same ``TOA_rho1/2`` + ``Eg_rho1/2`` spectral
    terms the LUT stores.
    """

    #: libRadtran core source tarball plus the auxiliary-data archives. The
    #: reptran (1nm molecular absorption) and OPAC optical-property
    #: (``continental_average`` aerosol) archives are fetched separately and
    #: merged into the build's ``data/`` directory.
    source_url: str = DEFAULT_LIBRADTRAN_SOURCE_URL
    reptran_url: str = DEFAULT_LIBRADTRAN_REPTRAN_URL
    optprop_url: str = DEFAULT_LIBRADTRAN_OPTPROP_URL
    #: Pre-unpacked source tree (skips download).
    source_dir: Path | None = None
    #: Build/cache root. Defaults to ``~/.cache/siac/libradtran/<profile>``.
    build_dir: Path | None = None
    #: Explicit path to a prebuilt ``uvspec`` binary (skips the build).
    uvspec_path: Path | None = None
    #: Explicit path to the libRadtran ``data/`` directory (inferred from the
    #: build tree when unset).
    data_dir: Path | None = None
    auto_build: bool = True
    build_profile: SixSBuildProfile = SixSBuildProfile.RELEASE
    #: Persistent on-disk cache of ``uvspec`` run outputs (per-node spectra),
    #: keyed by the exact deck text. A cache hit skips the ``uvspec`` subprocess
    #: entirely, so re-running a scene - or reusing nodes across the retrieval /
    #: diagnostics / correction scene-LUT builds - costs no radiative-transfer
    #: time. The compiled engine is always reused when present (see
    #: ``ensure_libradtran``); this additionally caches the *runs*.
    run_cache_enabled: bool = True
    #: Directory for the run-output cache. Defaults to
    #: ``~/.cache/siac/libradtran/uvspec_run_cache``. Keys embed the engine data
    #: path (which carries the libRadtran version), so different builds never
    #: collide in a shared cache directory.
    run_cache_dir: Path | None = None
    #: Max concurrent ``uvspec`` subprocesses. ``None`` -> a small memory-bound
    #: default (never ``os.cpu_count()``); always further capped by memory.
    native_threads: int | None = Field(default=None, ge=1)
    #: Hard ceiling (GB) on TOTAL estimated concurrent ``uvspec`` memory. The
    #: worker pool is sized so ``heaviest_process_estimate x workers`` stays under
    #: ``min(0.7 x detected-available-RAM, memory_budget_gb)``, and cgroup/SLURM
    #: ``--mem`` limits are detected and respected. Default 30 GB is a
    #: conservative single-node budget; raise it on large-memory machines. A
    #: single process larger than the budget runs one-at-a-time with a warning.
    memory_budget_gb: float | None = Field(default=30.0, gt=0.0)
    #: Per-``uvspec`` subprocess timeout (seconds): a wedged run raises instead
    #: of blocking a pool worker forever. Generous by default.
    uvspec_timeout_s: float = Field(default=1800.0, gt=0.0)
    #: Scene-grid sizing (mirrors the 6S scene-LUT budget). Each grid node runs
    #: ``uvspec`` once per surface albedo.
    scene_lut_min_pixels: int = Field(default=512, ge=1)
    scene_lut_max_nodes_per_axis: int = Field(default=4, ge=1)
    scene_lut_max_cases: int = Field(default=512, ge=1)
    #: DISORT polar-stream count (even, >= 2).
    number_of_streams: int = Field(default=16, ge=2)
    #: Base molecular-absorption band model (``mol_abs_param``) for the spectral
    #: *windows*. Defaults to ``"reptran medium"``: it holds <0.5% vs ``fine``
    #: outside the deep water bands at ~12 GB/process, versus ~50-75 GB for
    #: ``fine`` over the full 340-2500 nm range (the cause of past >100 GB OOMs).
    #: The deep H2O absorption bands are upgraded to ``"reptran fine"``
    #: automatically (see ``adaptive_deep_water_fine``), so fine accuracy is kept
    #: where it actually matters. Set ``"reptran fine"`` to match the production
    #: LUT setting everywhere (far higher memory), or ``"reptran"`` (coarse,
    #: bundled, no download). reptran ``medium``/``fine`` tables are auto-fetched
    #: from the same ``reptran`` archive by the build harness.
    mol_abs_param: str = "reptran medium"
    #: Optional per-spectral-region band models. Each entry ``(lo_nm, hi_nm,
    #: model)`` overrides ``mol_abs_param`` over ``[lo_nm, hi_nm]``; ``mol_abs_param``
    #: covers everything else. ``uvspec`` is run once per resulting contiguous
    #: segment and the 1 nm pieces are stitched. Use this to spend ``"reptran
    #: fine"`` only in the deep H2O bands (see ``DEEP_WATER_H2O_BANDS_NM``) while a
    #: cheap ``"reptran"`` base covers the windows. Region boundaries should fall
    #: OUTSIDE every band's response support (a runtime warning fires otherwise),
    #: since the band models differ ~1-2% across a seam. ``None`` (default) uses a
    #: single global ``mol_abs_param`` for the whole window.
    mol_abs_regions: tuple[tuple[float, float, str], ...] | None = None
    #: When ``True`` (default) AND ``mol_abs_regions`` is unset, the deep H2O
    #: absorption bands (``DEEP_WATER_H2O_BANDS_NM``) overlapping the window are
    #: run at ``"reptran fine"`` over the ``mol_abs_param`` base - the adaptive
    #: scheme: cheap base in the windows, fine only where strong water-vapour
    #: absorption needs it, keeping each ``uvspec`` process small (~7-12 GB) so a
    #: scene build fits ``memory_budget_gb``. An explicit ``mol_abs_regions``
    #: takes precedence; set ``False`` for one global ``mol_abs_param``.
    adaptive_deep_water_fine: bool = True
    #: Output spectral window (nm) and resampling step for each node's spectrum.
    #: The reference LUT spans 340-2600 nm at 1 nm; 340 nm covers every S2 band's
    #: response support.
    wavelength_min_nm: float = Field(default=340.0, gt=0.0)
    wavelength_max_nm: float = Field(default=2500.0, gt=0.0)
    wavelength_step_nm: float = Field(default=1.0, gt=0.0)
    #: The two Lambertian surface albedos of the two-albedo method; persisted as
    #: dataset attrs so the coefficient derivation reads them back.
    rho1: float = Field(default=0.15, ge=0.0, le=1.0)
    rho2: float = Field(default=0.5, ge=0.0, le=1.0)
    #: Per-pixel interpolation method over the assembled scene LUT.
    interpolation_method: LUTInterpolationMethod = LUTInterpolationMethod.LINEAR

    @field_validator(
        "source_dir", "build_dir", "uvspec_path", "data_dir", "run_cache_dir", mode="before"
    )
    @classmethod
    def normalize_local_paths(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)

    @field_validator("source_url", "reptran_url", "optprop_url")
    @classmethod
    def require_https_urls(cls, value: str) -> str:
        # The build harness downloads + compiles these archives; require TLS so a
        # plaintext or non-http(s) override cannot inject build inputs.
        if not str(value).lower().startswith("https://"):
            raise ValueError("libradtran archive URLs must use https://.")
        return value

    @field_validator("mol_abs_regions", mode="before")
    @classmethod
    def normalize_mol_abs_regions(cls, value: Any) -> tuple[tuple[float, float, str], ...] | None:
        if value is None:
            return None
        regions: list[tuple[float, float, str]] = []
        for item in value:
            lo, hi, model = item
            regions.append((float(lo), float(hi), str(model)))
        regions.sort(key=lambda r: r[0])
        return tuple(regions)

    @model_validator(mode="after")
    def validate_libradtran(self) -> LibRadtranAlgorithmConfig:
        if self.number_of_streams % 2 != 0:
            raise ValueError("libradtran.number_of_streams must be even (DISORT requirement).")
        if self.wavelength_max_nm <= self.wavelength_min_nm:
            raise ValueError("libradtran.wavelength_max_nm must exceed wavelength_min_nm.")
        if self.scene_lut_max_cases < self.scene_lut_max_nodes_per_axis:
            raise ValueError(
                "libradtran.scene_lut_max_cases must be >= scene_lut_max_nodes_per_axis."
            )
        if not 0.0 <= self.rho1 < self.rho2 <= 1.0:
            raise ValueError("libradtran requires 0 <= rho1 < rho2 <= 1.")
        if self.mol_abs_regions:
            prev_hi: float | None = None
            for lo, hi, _model in self.mol_abs_regions:
                if hi <= lo:
                    raise ValueError("libradtran.mol_abs_regions entries require lo_nm < hi_nm.")
                if lo < self.wavelength_min_nm or hi > self.wavelength_max_nm:
                    raise ValueError(
                        "libradtran.mol_abs_regions must lie within "
                        "[wavelength_min_nm, wavelength_max_nm]."
                    )
                if prev_hi is not None and lo < prev_hi:
                    raise ValueError("libradtran.mol_abs_regions must not overlap.")
                prev_hi = hi
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
    libradtran: LibRadtranAlgorithmConfig = Field(default_factory=LibRadtranAlgorithmConfig)
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
    #: Solver method. ``multigrid`` (default) is the Bayesian multi-grid
    #: inversion; ``surface_driven`` sweeps a fixed AOT axis at the prior TCWV
    #: and picks the per-pixel arg-min of the surface-prior mismatch (Rust
    #: kernel). The surface-driven method ignores the multi-grid / L-BFGS-B
    #: fields below and uses ``grid_search_aot_points`` + ``aot_bounds`` to
    #: build its AOT axis and ``surface_driven_*`` for pooling/backstop.
    method: SolverMethod = SolverMethod.MULTIGRID
    #: (surface_driven) half-width of the spatial cost-pooling window, in metres.
    #: ~600 m → a ~1.2 km median window, the robustness scale validated in the
    #: surface-driven harness. Converted to solver pixels via aerosol_resolution.
    surface_driven_pool_radius_m: float = Field(default=600.0, ge=0.0)
    #: (surface_driven) minimum finite cost samples required in a pooling window
    #: for a node to be considered (else that node is skipped for the pixel).
    surface_driven_pool_min_count: int = Field(default=1, ge=1)
    #: (surface_driven) when True, widen/tighten the prior (CAMS) backstop by
    #: AOT regime (loose when clean and at the high-AOD tail, tight in the
    #: moderate band where the surface signal is shallow); else flat 50 %.
    surface_driven_backstop_calibrated: bool = True
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
    #: Torch device for OmniCloudMask inference. Defaults to ``"cpu"``
    #: because PyTorch's MPS (Apple GPU) and CUDA backends are NOT
    #: bit-deterministic across processes — Apple's Metal driver and
    #: NVIDIA's cuDNN can pick different shader kernels per process
    #: launch, producing ULP-level softmax drift at edge-of-cloud pixels.
    #: Wave 18g of REVIEW_FIXES.md traces how a 79-pixel cloud-mask
    #: difference cascades through M5's Voronoi-fill amplifier
    #: (``nearest_seed_fill``) into 100 % of AOT pixels drifting with a
    #: ~4 % multiplicative bias.
    #:
    #: Set to ``"auto"`` to let OmniCloudMask pick its own default
    #: (typically the fastest available device — non-deterministic on
    #: GPU). Set to ``"cuda"`` or ``"mps"`` for explicit GPU
    #: acceleration. Keep ``"cpu"`` for scientific reproducibility.
    inference_device: str = "cpu"

    @field_validator("external_mask_path", mode="before")
    @classmethod
    def normalize_external_mask_path(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)


class AlgorithmsConfig(SIACBaseModel):
    surface_prior: SurfacePriorAlgorithmConfig = Field(default_factory=SurfacePriorAlgorithmConfig)
    rt: RTAlgorithmConfig = Field(default_factory=RTAlgorithmConfig)
    solver: SolverAlgorithmConfig = Field(default_factory=SolverAlgorithmConfig)
    cloud_mask: CloudMaskAlgorithmConfig = Field(default_factory=CloudMaskAlgorithmConfig)
