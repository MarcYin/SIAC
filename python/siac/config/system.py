"""Top-level runtime, output, and system configuration models."""

from __future__ import annotations

from pathlib import Path  # noqa: TC003
from typing import Any

from pydantic import Field, field_validator

from siac.config._base import SIACBaseModel
from siac.config.algorithms import AlgorithmsConfig
from siac.config.providers import (
    AuthConfig,
    PathsConfig,
    ProvidersConfig,
)
from siac.config.types import (
    BOADType,
    ExecutionBackend,
    LogLevel,
    OutputCompression,
    OutputFormat,
    _coerce_pathlike,
)


class ExecutionRuntimeConfig(SIACBaseModel):
    backend: ExecutionBackend = ExecutionBackend.THREAD
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

    @field_validator("stage_timeouts")
    @classmethod
    def validate_stage_timeouts(cls, value: dict[str, float]) -> dict[str, float]:
        for stage_name, timeout_s in value.items():
            if not stage_name:
                raise ValueError("stage_timeouts keys must be non-empty stage names")
            if timeout_s <= 0.0:
                raise ValueError("stage_timeouts values must be > 0")
        return value


class RuntimeConfig(SIACBaseModel):
    log_level: LogLevel = LogLevel.INFO
    n_jobs: int = -1
    chunk_size: int = Field(default=2048, gt=0)
    execution: ExecutionRuntimeConfig = Field(default_factory=ExecutionRuntimeConfig)


class OutputDefaultsConfig(SIACBaseModel):
    output_dir: Path | None = None
    format: OutputFormat = OutputFormat.COG
    compression: OutputCompression = OutputCompression.DEFLATE
    include_rgb: bool = True
    include_uncertainty: bool = True
    include_auxiliary: bool = True
    #: Wave 18f: whether to render the decorative PNG previews
    #: (cloud-mask overlay, false-colour, AOT/TCWV colour maps, AOT
    #: scatter plots). These add ~12-16 s of wall-clock per scene on
    #: T33KWP because they run the per-pixel scaling + PIL render on
    #: the full-resolution BOA bands. Set to ``False`` for production
    #: runs that don't need the visual quicklooks — the RGB GeoTIFF
    #: quicklook (gated by ``include_rgb``) is still written.
    include_previews: bool = True
    #: Max edge length (pixels) for preview PNG outputs. Wave 18f
    #: downsamples the input BOA bands to this size before the
    #: per-pixel PIL render — the cloud-mask and false-colour previews
    #: were spending most of their wall-clock processing 10980×10980
    #: arrays for what eventually displays at ~1000-2000 px anyway. A
    #: 2048-px cap captures essentially all the visual fidelity at a
    #: fraction of the cost.
    preview_max_size_px: int = Field(default=2048, gt=0)
    skip_correction: bool = False
    boa_dtype: BOADType = BOADType.FLOAT32
    boa_scale: float = Field(default=10000.0, gt=0.0)
    boa_nodata: float = 0.0
    reopen_streamed_boa: bool = True

    @field_validator("output_dir", mode="before")
    @classmethod
    def normalize_output_dir(cls, value: Any) -> Path | None:
        return _coerce_pathlike(value)


class OutputConfig(SIACBaseModel):
    defaults: OutputDefaultsConfig = Field(default_factory=OutputDefaultsConfig)


class SystemConfig(SIACBaseModel):
    paths: PathsConfig = Field(default_factory=PathsConfig)
    auth: AuthConfig = Field(default_factory=AuthConfig)
    providers: ProvidersConfig = Field(default_factory=ProvidersConfig)
    algorithms: AlgorithmsConfig = Field(default_factory=AlgorithmsConfig)
    runtime: RuntimeConfig = Field(default_factory=RuntimeConfig)
    output: OutputConfig = Field(default_factory=OutputConfig)
