"""Month and sample-date helpers for SWIR Route-B surface priors."""

from __future__ import annotations

import calendar
from dataclasses import dataclass
from datetime import datetime
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import collections.abc

_HISTORY_YEARS = 5
_HISTORY_MONTH_OFFSETS = (-1, 0, 1)
_WEEKLY_STEP_DAYS = 7


@dataclass(frozen=True)
class _MonthSpec:
    year: int
    month: int
    center_time: datetime
    temporal_window: int
    sample_dates: tuple[datetime, ...]


def _iter_history_months(obs_time: datetime) -> list[tuple[int, int]]:
    months: list[tuple[int, int]] = []
    for year_offset in range(1, _HISTORY_YEARS + 1):
        base_year = obs_time.year - year_offset
        for month_offset in _HISTORY_MONTH_OFFSETS:
            month = obs_time.month + month_offset
            year = base_year
            while month < 1:
                month += 12
                year -= 1
            while month > 12:
                month -= 12
                year += 1
            months.append((year, month))
    return months


def _build_explicit_month_specs(
    year_months: collections.abc.Sequence[tuple[int, int]],
    *,
    template_time: datetime | None,
) -> list[_MonthSpec]:
    if not year_months:
        raise ValueError("year_months must not be empty")

    specs: list[_MonthSpec] = []
    seen: set[tuple[int, int]] = set()
    for year, month in sorted((int(year), int(month)) for year, month in year_months):
        if not 1 <= month <= 12:
            raise ValueError(f"month must be between 1 and 12, got {month}")
        key = (year, month)
        if key in seen:
            raise ValueError(f"Duplicate monthly composite selection: {year:04d}-{month:02d}")
        seen.add(key)
        center_template = template_time or datetime(year, month, 15, 12, 0, 0)
        specs.append(
            _MonthSpec(
                year=year,
                month=month,
                center_time=_month_center_datetime(year, month, center_template),
                temporal_window=max(1, calendar.monthrange(year, month)[1] // 2 + 1),
                sample_dates=_weekly_sample_dates(year, month),
            )
        )
    return specs


def _month_center_datetime(year: int, month: int, template_time: datetime) -> datetime:
    n_days = calendar.monthrange(year, month)[1]
    center_day = int(np.ceil(n_days / 2.0))
    return template_time.replace(year=year, month=month, day=center_day)


def _weekly_sample_dates(year: int, month: int) -> tuple[datetime, ...]:
    n_days = calendar.monthrange(year, month)[1]
    days = list(range(1, n_days + 1, _WEEKLY_STEP_DAYS))
    if days[-1] != n_days and (n_days - days[-1]) >= (_WEEKLY_STEP_DAYS // 2):
        days.append(n_days)
    return tuple(datetime(year, month, day) for day in days)


def _select_month_mask(time_values: np.ndarray, *, year: int, month: int) -> np.ndarray:
    month_strings = np.asarray(time_values, dtype="datetime64[D]").astype("datetime64[M]")
    target = np.datetime64(f"{year:04d}-{month:02d}", "M")
    return np.asarray(month_strings == target, dtype=np.bool_)
