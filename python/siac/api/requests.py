"""Typed public workflow requests."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from typing import TYPE_CHECKING, TypeAlias

if TYPE_CHECKING:
    from siac.adapters.auth import CredentialManager
    from siac.adapters.data.s2_data_source import S2Query
    from siac.config import SIACConfig
    from siac.domain.aoi import AOI


AOISpec: TypeAlias = "AOI | Path | str | tuple[float, float, float, float] | list[float] | None"
DateSpec: TypeAlias = date | datetime | str | None
PathLike: TypeAlias = str | Path
Sentinel2QuerySpec: TypeAlias = "S2Query | str | Path"


@dataclass(frozen=True)
class SceneProcessRequest:
    config: SIACConfig
    input_path: PathLike
    output_path: PathLike | None = None
    aoi: AOISpec = None
    auth: CredentialManager | None = None


@dataclass(frozen=True)
class Sentinel2ResolveRequest:
    config: SIACConfig
    query: Sentinel2QuerySpec
    auth: CredentialManager | None = None


@dataclass(frozen=True)
class Sentinel2SearchRequest:
    tile: str | None = None
    date: DateSpec = None
    date_value: DateSpec = None
    start_date: DateSpec = None
    end_date: DateSpec = None
    bbox: tuple[float, float, float, float] | None = None
    max_cloud_cover: float = 80.0
    backend: str = "cdse"
    config: SIACConfig | None = None
    auth: CredentialManager | None = None


@dataclass(frozen=True)
class Sentinel2ProcessRequest:
    config: SIACConfig
    query: Sentinel2QuerySpec
    output_path: PathLike | None = None
    aoi: AOISpec = None
    auth: CredentialManager | None = None
