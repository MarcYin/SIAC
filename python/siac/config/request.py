"""Per-run request configuration model."""

from __future__ import annotations

from pathlib import Path  # noqa: TC003
from typing import Any

from pydantic import field_validator

from siac.config._base import SIACBaseModel
from siac.config.types import SensorName, _coerce_path_or_url, _coerce_pathlike


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
