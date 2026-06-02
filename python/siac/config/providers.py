"""Path, authentication, and provider configuration models."""

from __future__ import annotations

from pathlib import Path  # noqa: TC003
from typing import Any

from pydantic import Field, field_validator

from siac.config._base import SIACBaseModel
from siac.config.types import (
    DEFAULT_LUT_URL,
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
    #: ``"auto" | "earth-search" | "pc" | "hls" | "mcd43a4"``.
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
