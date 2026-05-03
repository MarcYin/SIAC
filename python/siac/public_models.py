"""Public SIAC data models shared across API layers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path


@dataclass(frozen=True)
class PreparedMonthlyCompositeBuildResult:
    store_path: Path
    period_count: int
    periods: tuple[str, ...]
    source_name: str | None
    source_band_names: tuple[str, ...]
    representation: str
