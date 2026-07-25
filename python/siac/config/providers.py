"""Path, authentication, and provider configuration models."""

from __future__ import annotations

from pathlib import Path  # noqa: TC003
from typing import Any

from pydantic import Field, field_validator, model_validator

from siac.config._base import SIACBaseModel
from siac.config.types import (
    DEFAULT_LUT_URL,
    AODFusionOp,
    AtmoProviderKind,
    BRDFProviderKind,
    MonthlyCompositeProviderKind,
    S2Backend,
    S2ProcessingLevel,
    TemporalInterpolation,
    _coerce_path_or_url,
    _coerce_pathlike,
)


class CachePathsConfig(SIACBaseModel):
    atmo: Path | None = None
    brdf: Path | None = None
    s2: Path | None = None
    #: Content-addressed cache directory for cloud-mask model outputs.
    #: Skips the OmniCloudMask PyTorch inference (~20-25 s) on subsequent
    #: runs over the same TOA inputs. Set to ``None`` to disable caching.
    cloud: Path | None = None
    #: Content-addressed cache directory for reprojected priors and
    #: other geometric resamples (wave 18 opt 3). When set, the
    #: ``cached_reproject_match`` helper persists each reprojection
    #: output under a SHA-256 key over (target-grid signature,
    #: source-identity, resampling method, format version). On
    #: subsequent runs the reprojection is a single-digit-ms NetCDF
    #: load instead of a multi-second ``warp.reproject`` call. The MCD43
    #: monthly-composite reprojection — the largest single contributor
    #: to the wave-17 35 s of ``warp.py:reproject`` wall-clock — is the
    #: primary consumer. Set to ``None`` to disable.
    reproject: Path | None = None

    @field_validator("atmo", "brdf", "s2", "cloud", "reproject", mode="before")
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
    kind: AtmoProviderKind = AtmoProviderKind.CAMS
    data_path: str | Path | None = None
    cache_dir: Path | None = None
    download_missing: bool = True
    temporal_interpolation: TemporalInterpolation = TemporalInterpolation.NEAREST
    #: (MAIAC / MCD19 + VNP19 only) decode the ``AOD_QA`` bit field and keep only
    #: best-quality pixels (clear cloud mask + no adjacency + best AOD QA level)
    #: instead of the loose ``QA > 0`` filter. Matches the validated ``maiac_qa``
    #: best-quality recipe — the loose filter reads systematically higher at
    #: clean-coastal/polar sites (~0.63 vs ~0.54), which the surface-driven
    #: backstop then anchors to. ``False`` (default) preserves the legacy loose
    #: filter (and the multigrid atmo prior). No effect for CAMS / MERRA-2.
    maiac_best_quality_qa: bool = False
    #: Accept a fabricated constant AOD when the aerosol source has no usable
    #: data for a scene. Off by default: a surface-driven retrieval is
    #: prior-limited, so an invented prior silently becomes the answer instead
    #: of failing visibly, and reads far too low exactly where aerosol is thick.
    allow_default_prior: bool = False
    #: Additional AOD sources fused into the aerosol prior. ``kind`` remains the
    #: primary provider — it supplies water vapour, ozone, terrain and the grid —
    #: while these sources contribute aerosol optical depth only, combined under
    #: :attr:`fuse_aod_op`.
    #:
    #: The surface-driven retrieval defers to the aerosol prior wherever the
    #: visible bands constrain AOT weakly (most of the clean and moderate range),
    #: so the prior's centre is the dominant error term. Both routine sources
    #: under-estimate independently — MAIAC by ~0.05 (0.17 above AOD 0.5), CAMS by
    #: ~0.08 (0.27) — so ``max`` over the two removes most of the shared bias
    #: without fitting against reference data. Measured on the 152-matchup
    #: campaign: within-EE 79.2% → 84.6%, bias −0.021 → +0.001, and the
    #: prior-limited moderate band 69% → 84%.
    #:
    #: ``("cams",)`` with ``kind="mcd19"`` is the validated configuration.
    fuse_aod_with: tuple[AtmoProviderKind, ...] = ()
    #: How :attr:`fuse_aod_with` sources combine with the primary AOD. ``max`` is
    #: validated; ``mean`` preserves the shared under-estimate (it scored 27% in
    #: the high-AOD band against 54% for ``max``) and ``min`` compounds it.
    fuse_aod_op: AODFusionOp = "max"
    #: Read the aerosol optical depth from a prepared per-scene store instead of
    #: deriving it live. Live granule extraction over a small AOI can yield
    #: nothing (missing granule, cloud screening, QA rejection), in which case a
    #: default is substituted — and because the surface-driven retrieval is
    #: prior-limited wherever the visible bands constrain AOT weakly, that
    #: fabricated prior silently becomes the answer. Preparing the value offline
    #: removes the failure mode. The store holds one ``<scene_key>.json`` per
    #: scene with an ``aot`` (and optional ``aot_unc``); the remaining state
    #: still comes from :attr:`kind`.
    prepared_scalar_path: Path | None = None
    #: Entry to read from :attr:`prepared_scalar_path` for this run.
    prepared_scalar_scene_key: str | None = None
    #: Fail when the prepared entry is missing (default), rather than silently
    #: falling back to the live provider.
    prepared_scalar_required: bool = True
    #: Optional prepared total-column water vapour (cm precipitable water) for
    #: the scene, overriding the live provider's field the same way the
    #: prepared AOD does. The anchor bands (B8A/B11/B12) sit in water
    #: absorption, so the tau-dependent surface prior inherits any TCWV error.
    prepared_scalar_tcwv_cm: float | None = None

    @field_validator("fuse_aod_with", mode="before")
    @classmethod
    def normalize_fuse_sources(cls, value: Any) -> tuple[Any, ...]:
        if value is None:
            return ()
        if isinstance(value, str):
            return (value,)
        return tuple(value)

    @model_validator(mode="after")
    def validate_fusion_sources(self) -> AtmoProviderConfig:
        if self.kind in self.fuse_aod_with:
            raise ValueError(
                f"providers.atmo.fuse_aod_with must not repeat the primary provider "
                f"{self.kind.value!r}; list only the additional AOD sources."
            )
        if len(set(self.fuse_aod_with)) != len(self.fuse_aod_with):
            raise ValueError("providers.atmo.fuse_aod_with must not contain duplicates.")
        return self

    @field_validator("data_path", mode="before")
    @classmethod
    def normalize_data_path(cls, value: Any) -> str | Path | None:
        return _coerce_path_or_url(value)

    @field_validator("cache_dir", mode="before")
    @classmethod
    def normalize_local_paths(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)


class BRDFProviderConfig(SIACBaseModel):
    kind: BRDFProviderKind = BRDFProviderKind.MCD43
    data_path: str | Path | None = None
    cache_dir: Path | None = None
    temporal_window: int = Field(default=16, ge=1, le=32)
    use_cache: bool = True

    @field_validator("data_path", mode="before")
    @classmethod
    def normalize_data_path(cls, value: Any) -> str | Path | None:
        return _coerce_path_or_url(value)

    @field_validator("cache_dir", mode="before")
    @classmethod
    def normalize_cache_dir(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)


class S2ProviderConfig(SIACBaseModel):
    backend: S2Backend = S2Backend.LOCAL
    cache_dir: Path | None = None
    max_cloud_cover: float = Field(default=80.0, ge=0.0, le=100.0)
    prefer_newest_baseline: bool = True
    processing_level: S2ProcessingLevel = S2ProcessingLevel.L1C

    @field_validator("cache_dir", mode="before")
    @classmethod
    def normalize_cache_dir(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)


class MonthlyCompositeProviderConfig(SIACBaseModel):
    kind: MonthlyCompositeProviderKind = MonthlyCompositeProviderKind.GENERATED_BRDF
    user_callable: Any | None = None
    store_path: Path | None = None
    strict_coverage: bool = True

    # --- bestpixel provider (kind == "bestpixel") -----------------------
    #: STAC source the ``bestpixel`` package composites from. One of
    #: ``"auto" | "earth-search" | "pc" | "hls" | "hls-s30" | "mcd43a4"``.
    #: ``"hls"`` is mixed HLS L30+S30; ``"hls-s30"`` restricts the pool to
    #: Planetary Computer HLS S30.
    bestpixel_endpoint: str = "pc"
    #: Canonical bestpixel band names to request. ``None`` selects a
    #: default set spanning visible→SWIR (coastal, blue, green, red, nir,
    #: swir16, swir22) which covers SIAC's Route-B visible + query bands.
    bestpixel_bands: tuple[str, ...] | None = None
    #: Number of full calendar years before the observation year to build
    #: monthly composites for (e.g. 5 → the five years preceding the scene).
    bestpixel_lookback_years: int = 5
    #: Calendar months (1–12) to composite. ``None`` (default) composites only
    #: the scene's own acquisition month each year, so the surface prior is
    #: seasonally matched (e.g. a March scene -> March of each lookback year).
    #: Set explicitly to widen the seasonal window (e.g. [2, 3, 4]) or to
    #: composite a full year ([1, 2, ..., 12]).
    bestpixel_months: tuple[int, ...] | None = None
    #: When ``bestpixel_months`` is ``None``, widen the seasonal window to the
    #: scene month ± this many calendar months (wrapping 1–12): e.g. a May scene
    #: with ``1`` composites Apr/May/Jun of each lookback year. ``0`` (default)
    #: keeps the scene month only. The validated ``comp_ref`` harness used ±1
    #: (a 3-month seasonal window × the lookback years), tripling the realization
    #: count so the temporal median + RMSE are better constrained. Ignored when
    #: ``bestpixel_months`` is set explicitly.
    bestpixel_seasonal_window_months: int = 0
    #: Output CRS bestpixel composites in: ``"utm"`` or ``"native"``.
    bestpixel_output_crs: str = "utm"
    #: Fixed clearest-observations-per-chunk depth.
    bestpixel_top_k: int = 3
    #: Drop source scenes whose cloud cover exceeds this percent.
    bestpixel_max_cloud_cover: float = 90.0
    #: Composite build resolution in metres. ``None`` builds composites at
    #: the Route-B database resolution SIAC requests at runtime.
    bestpixel_resolution_m: float | None = None
    #: When set (finer than ``bestpixel_resolution_m``), fetch composites at
    #: this resolution and then spatially AREA-AVERAGE them down to
    #: ``bestpixel_resolution_m`` before use as the surface prior. Averaging
    #: N:1 reduces per-pixel noise ~sqrt(N) — trades resolution for a smoother,
    #: lower-outlier prior. ``None`` fetches directly at the build resolution
    #: (no averaging).
    bestpixel_fetch_resolution_m: float | None = None
    #: Optional on-disk cache directory passed straight to bestpixel.
    bestpixel_disk_cache: Path | None = None

    @field_validator("store_path", "bestpixel_disk_cache", mode="before")
    @classmethod
    def normalize_store_path(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)


class ProvidersConfig(SIACBaseModel):
    atmo: AtmoProviderConfig = Field(default_factory=AtmoProviderConfig)
    brdf: BRDFProviderConfig = Field(default_factory=BRDFProviderConfig)
    s2: S2ProviderConfig = Field(default_factory=S2ProviderConfig)
    monthly_composites: MonthlyCompositeProviderConfig = Field(
        default_factory=MonthlyCompositeProviderConfig
    )
