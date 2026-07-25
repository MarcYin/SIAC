"""Aerosol-species candidates for native 6S surface-driven retrievals.

This ports the research-harness Aerosol_cci climatology path into package data.
For a scene centre and month it reads the monthly 1-degree climatology, finds
the nearest aerosol-type rows, and converts each row into a 6S
``multimodal_log_normal`` aerosol setup.
"""

from __future__ import annotations

import csv
import os
from functools import lru_cache
from importlib import resources
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import xarray as xr

from siac.config.algorithms import (
    RTAerosolSetupConfig,
    SixSAerosolDistributionComponentConfig,
    SixSAerosolDistributionConfig,
)

if TYPE_CHECKING:
    from collections.abc import Mapping

_AERO_DATA_PACKAGE = "siac.data.aero_species"
_DEFAULT_PERCENT_FRACTIONS = np.array([0.0, 0.0, 0.0, 100.0], dtype=np.float64)
_CCI_COMPONENTS: dict[str, dict[str, float]] = {
    "dust": {"rmean": 0.788, "sigma": 1.822, "n": 1.56, "k": 0.0018},
    "sea_salt": {"rmean": 0.788, "sigma": 1.822, "n": 1.40, "k": 1e-8},
    "fine_weak": {"rmean": 0.07, "sigma": 1.7, "n": 1.40, "k": 0.003},
    "fine_strong": {"rmean": 0.07, "sigma": 1.7, "n": 1.50, "k": 0.040},
}


def _resource(name: str) -> Any:
    return resources.files(_AERO_DATA_PACKAGE).joinpath(name)


@lru_cache(maxsize=1)
def _climatology_dataset() -> xr.Dataset:
    with resources.as_file(_resource("aerosolClimatology.nc")) as path:
        ds = xr.open_dataset(path)
        try:
            return ds.load()
        finally:
            ds.close()


@lru_cache(maxsize=1)
def _type_lut() -> np.ndarray:
    rows: list[tuple[float, float, float, float]] = []
    with _resource("aerosol_type_lut.csv").open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                (
                    float(row["dust"]),
                    float(row["sea_salt"]),
                    float(row["fine_mode_high_absorption"]),
                    float(row["fine_mode_low_absorption"]),
                )
            )
    if not rows:
        raise ValueError("Packaged aerosol_type_lut.csv contains no aerosol rows.")
    return cast("np.ndarray", np.asarray(rows, dtype=np.float64))


def _normalize_lon_for_dataset(lon: float, ds: xr.Dataset) -> float:
    lon_values = np.asarray(ds["lon"].values, dtype=np.float64)
    if lon_values.size == 0:
        return float(lon)
    lon_min = float(np.nanmin(lon_values))
    lon_max = float(np.nanmax(lon_values))
    if lon_max > 180.0 and lon < 0.0:
        return float(lon % 360.0)
    if lon_min >= -180.0 and lon_max <= 180.0 and lon > 180.0:
        return float(((lon + 180.0) % 360.0) - 180.0)
    return float(lon)


def _normalize_percent_fractions(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    if values.shape != (4,) or not np.all(np.isfinite(values)):
        return cast("np.ndarray", _DEFAULT_PERCENT_FRACTIONS.copy())
    clipped = np.clip(values, 0.0, None)
    total = float(np.sum(clipped))
    if total <= 0.0:
        return cast("np.ndarray", _DEFAULT_PERCENT_FRACTIONS.copy())
    return cast("np.ndarray", clipped * (100.0 / total))


def climatology_fraction_percentages(lon: float, lat: float, month: int) -> np.ndarray:
    """Return climatological percentages in CCI order.

    The returned order is ``[dust, sea_salt, fine_strong, fine_weak]`` and the
    values sum to 100.
    """

    month_int = int(month)
    if not 1 <= month_int <= 12:
        raise ValueError(f"month must be in 1..12, got {month!r}")
    ds = _climatology_dataset()
    lon_norm = _normalize_lon_for_dataset(float(lon), ds)
    pt = ds.sel(time=month_int).sel(lat=float(lat), lon=lon_norm, method="nearest")
    fine = float(pt["fine_of_total_fraction"])
    less_abs = float(pt["lessAbs_of_fine_fraction"])
    dust_coarse = float(pt["dust_of_coarse_fraction"])
    dust = dust_coarse * (1.0 - fine)
    sea_salt = (1.0 - dust_coarse) * (1.0 - fine)
    fine_strong = fine * (1.0 - less_abs)
    fine_weak = fine * less_abs
    return _normalize_percent_fractions(
        np.array([dust, sea_salt, fine_strong, fine_weak], dtype=np.float64) * 100.0
    )


def candidate_fraction_sets(
    lon: float,
    lat: float,
    month: int,
    *,
    n: int = 3,
) -> tuple[dict[str, float], ...]:
    """Return nearest CCI aerosol-type component fractions for a scene.

    Each item maps ``dust``, ``sea_salt``, ``fine_strong``, and ``fine_weak`` to
    fractions summing to 1.0. Matching follows the harness: nearest rows by L1
    distance in the four climatological percentage fields.
    """

    clim = climatology_fraction_percentages(lon, lat, month)
    lut = _type_lut()
    count = min(max(1, int(n)), int(lut.shape[0]))
    indices = np.argsort(np.abs(lut - clim[np.newaxis, :]).sum(axis=1))[:count]
    out: list[dict[str, float]] = []
    for idx in indices:
        dust, sea_salt, fine_strong, fine_weak = lut[int(idx)] / 100.0
        out.append(
            {
                "dust": float(dust),
                "sea_salt": float(sea_salt),
                "fine_strong": float(fine_strong),
                "fine_weak": float(fine_weak),
            }
        )
    return tuple(out)


def climatology_fraction_set(lon: float, lat: float, month: int) -> dict[str, float]:
    """Return the continuous monthly CCI component fractions for one scene."""

    percentages = climatology_fraction_percentages(lon, lat, month)
    return {
        name: float(percentages[index] / 100.0)
        for index, name in enumerate(("dust", "sea_salt", "fine_strong", "fine_weak"))
    }


def _component_config(
    name: str, percentage_density: float
) -> SixSAerosolDistributionComponentConfig:
    component = _CCI_COMPONENTS[name]
    return SixSAerosolDistributionComponentConfig(
        rmean=component["rmean"],
        sigma=component["sigma"],
        percentage_density=float(percentage_density),
        refr_real=tuple([component["n"]] * 20),
        refr_imag=tuple([component["k"]] * 20),
    )


def cci_distribution(
    fractions: Mapping[str, float],
    *,
    min_fraction: float = 0.01,
) -> SixSAerosolDistributionConfig:
    """Build a native-6S ``multimodal_log_normal`` distribution config."""

    kept = {
        name: float(value)
        for name, value in fractions.items()
        if name in _CCI_COMPONENTS and np.isfinite(float(value)) and float(value) > min_fraction
    }
    if not kept:
        kept = {"fine_weak": 1.0}
    total = sum(kept.values())
    components = tuple(
        _component_config(name, 100.0 * value / total) for name, value in kept.items()
    )
    return SixSAerosolDistributionConfig(rmin=0.005, rmax=20.0, components=components)


def candidate_cci_aerosol_setups(
    lon: float,
    lat: float,
    month: int,
    *,
    n: int = 3,
) -> tuple[RTAerosolSetupConfig, ...]:
    """Return native-6S aerosol setup configs for nearest CCI candidates."""

    return tuple(
        RTAerosolSetupConfig(
            profile="multimodal_log_normal",
            distribution=cci_distribution(fractions),
        )
        for fractions in candidate_fraction_sets(lon, lat, month, n=n)
    )


_DEFAULT_FINE_STRONG_CAP = 0.25


def climatology_cci_aerosol_setup(
    lon: float,
    lat: float,
    month: int,
) -> RTAerosolSetupConfig:
    """Return one native-6S setup using the continuous monthly CCI mixture.

    The soot-like ``fine_strong`` fraction is capped (excess moved to
    ``fine_weak``): large soot fractions in the CCI climatology over-absorb
    and inflate retrieved AOT, while capping at 0.25 leaves the sites that
    need absorbing aerosol unaffected. ``SIAC_SPECIES_FINE_STRONG_CAP``
    overrides the cap (``none`` disables it).
    """

    fractions = climatology_fraction_set(lon, lat, month)
    # Optional per-scene composition override from CAMS speciated AOD
    # (dust, sea_salt, fine_strong=black_carbon, fine_weak=organic+sulphate),
    # supplied as a normalized "d,s,fs,fw" string. Replaces the climatology
    # mixture with the actual day's mix; the fine-strong cap below still applies.
    cams_raw = os.environ.get("SIAC_CAMS_FRACTIONS", "")
    if cams_raw.strip():
        vals = [float(x) for x in cams_raw.split(",")]
        if len(vals) == 4 and sum(vals) > 0.0:
            total = sum(vals)
            fractions = dict(
                zip(("dust", "sea_salt", "fine_strong", "fine_weak"), [v / total for v in vals])
            )
    cap_raw = os.environ.get("SIAC_SPECIES_FINE_STRONG_CAP", "")
    if cap_raw.lower() in {"none", "off"}:
        cap = None
    else:
        cap = float(cap_raw) if cap_raw else _DEFAULT_FINE_STRONG_CAP
    if cap is not None:
        excess = fractions["fine_strong"] - cap
        if excess > 0.0:
            fractions["fine_strong"] = cap
            fractions["fine_weak"] += excess
    return RTAerosolSetupConfig(
        profile="multimodal_log_normal",
        distribution=cci_distribution(fractions),
    )
