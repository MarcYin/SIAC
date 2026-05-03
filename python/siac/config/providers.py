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

    @field_validator("store_path", mode="before")
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
