"""Resolved runtime configuration models."""

from __future__ import annotations

from pathlib import Path  # noqa: TC003
from typing import Any

from pydantic import Field

from siac.config._base import SIACBaseModel
from siac.config.algorithms import (
    CloudMaskAlgorithmConfig,
    RTAlgorithmConfig,
    SolverAlgorithmConfig,
    SpectralMappingAlgorithmConfig,
    SurfacePriorAlgorithmConfig,
)
from siac.config.providers import (
    AtmoProviderConfig,
    BRDFProviderConfig,
    CachePathsConfig,
    MonthlyCompositeProviderConfig,
    S2ProviderConfig,
)
from siac.config.system import OutputConfig, RuntimeConfig  # noqa: TC001
from siac.config.types import SensorName  # noqa: TC001

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
    sensor: SensorName = SensorName.AUTO
    aoi: dict[str, Any] | Path | str | tuple[float, float, float, float] | list[float] | None = None
    s2_query: str | Path | None = None


class ResolvedConfig(SIACBaseModel):
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
    def aoi(
        self,
    ) -> dict[str, Any] | Path | str | tuple[float, float, float, float] | list[float] | None:
        return self.run.aoi
