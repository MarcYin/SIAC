"""Stage a tile-wide QA-filtered MAIAC prior on the native product grid.

This avoids resampling sparse best-quality pixels before aggregation.  The
operational provider can otherwise turn an AOI containing valid MCD19 pixels
into an all-missing target grid and silently replace it with the default AOD.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import tempfile
from collections import defaultdict
from datetime import date, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from pyproj import Transformer

from siac.adapters.atmo.mcd19_earthaccess import MCD19AODProvider
from siac.adapters.earthdata_common import parse_granule_date

DEFAULT_ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_HALF_DEG = 0.12

if TYPE_CHECKING:
    import xarray as xr


def _atomic_json(path: Path, record: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as stream:
        json.dump(record, stream, indent=2)
        stream.write("\n")
        temporary = Path(stream.name)
    temporary.replace(path)


def _projected_bounds(
    array: xr.DataArray,
    bounds_4326: tuple[float, float, float, float],
) -> tuple[float, float, float, float]:
    crs = array.rio.crs
    if crs is None:
        raise ValueError("MAIAC native tile has no CRS")
    west, south, east, north = bounds_4326
    transformer = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    projected = [transformer.transform(lon, lat) for lon in (west, east) for lat in (south, north)]
    x = [point[0] for point in projected]
    y = [point[1] for point in projected]
    return min(x), min(y), max(x), max(y)


def native_window_values(
    aot: xr.DataArray,
    aot_unc: xr.DataArray,
    bounds_4326: tuple[float, float, float, float],
) -> tuple[np.ndarray, np.ndarray]:
    """Return paired finite native pixels whose centres fall inside an AOI."""
    xmin, ymin, xmax, ymax = _projected_bounds(aot, bounds_4326)
    x = np.asarray(aot.coords["x"].values, dtype=np.float64)
    y = np.asarray(aot.coords["y"].values, dtype=np.float64)
    x_index = np.flatnonzero((x >= xmin) & (x <= xmax))
    y_index = np.flatnonzero((y >= ymin) & (y <= ymax))
    if not x_index.size or not y_index.size:
        return np.zeros(0, dtype=np.float64), np.zeros(0, dtype=np.float64)
    aot_values = np.asarray(aot.values[np.ix_(y_index, x_index)], dtype=np.float64)
    unc_values = np.asarray(aot_unc.values[np.ix_(y_index, x_index)], dtype=np.float64)
    valid = np.isfinite(aot_values) & np.isfinite(unc_values) & (aot_values >= 0.0)
    return aot_values[valid], unc_values[valid]


def select_day(
    values_by_day: dict[date, tuple[list[float], list[float]]],
    observation_day: date,
) -> tuple[date, np.ndarray, np.ndarray] | None:
    """Select the nearest day, preferring more QA pixels at equal distance."""
    candidates = []
    for day, (aot_values, _unc_values) in values_by_day.items():
        if not aot_values:
            continue
        candidates.append((abs((day - observation_day).days), -len(aot_values), day))
    if not candidates:
        return None
    _offset, _negative_count, day = min(candidates)
    aot_values, unc_values = values_by_day[day]
    return day, np.asarray(aot_values, dtype=np.float64), np.asarray(unc_values, dtype=np.float64)


def calibrated_uncertainty(aot: float, native_uncertainty: float, count: int) -> float:
    base = max(native_uncertainty, 0.04)
    high_aod = 0.18 * max(aot - 0.3, 0.0)
    sparse = 0.05 if count < 8 else 0.0
    return float(np.clip(base + high_aod + sparse, 0.04, 0.6))


def stage_one(
    matchup: dict[str, str],
    *,
    output_dir: Path,
    cache_dir: Path,
    half_deg: float = DEFAULT_HALF_DEG,
    aggregation: str = "median",
) -> dict[str, Any]:
    matchup_id = matchup["matchup_id"]
    lon = float(matchup["longitude"])
    lat = float(matchup["latitude"])
    observation_time = datetime.fromisoformat(matchup["sensing_time_utc"].replace("Z", "+00:00"))
    bounds = (lon - half_deg, lat - half_deg, lon + half_deg, lat + half_deg)
    provider = MCD19AODProvider(
        cache_dir=cache_dir,
        best_quality_qa=True,
        temporal_window_days=2,
        max_granules=40,
    )
    paths, short_name = provider._download_granules(bounds, "EPSG:4326", observation_time)
    values_by_day: dict[date, tuple[list[float], list[float]]] = defaultdict(lambda: ([], []))
    for path in paths:
        tile = provider._load_native_tile(path, obs_time=observation_time)
        aot_values, unc_values = native_window_values(tile["aot"], tile["aot_unc"], bounds)
        day = parse_granule_date(path).date()
        values_by_day[day][0].extend(aot_values.tolist())
        values_by_day[day][1].extend(unc_values.tolist())

    selected = select_day(values_by_day, observation_time.date())
    if selected is None:
        record = {
            "matchup_id": matchup_id,
            "site": matchup.get("site", ""),
            "status": "NO_VALID_AOD",
            "source": "MCD19A2",
            "short_name": short_name,
            "aot": None,
            "half_deg": half_deg,
            "aggregation": aggregation,
            "searched_granules": len(paths),
        }
    else:
        day, aot_values, unc_values = selected
        aot = float(np.median(aot_values) if aggregation == "median" else np.mean(aot_values))
        native_uncertainty = float(np.median(unc_values))
        if not math.isfinite(aot) or aot < 0.0:
            raise ValueError(f"Invalid aggregated MAIAC AOD for {matchup_id}: {aot!r}")
        record = {
            "matchup_id": matchup_id,
            "site": matchup.get("site", ""),
            "status": "OK",
            "source": "MCD19A2",
            "short_name": short_name,
            "aot": aot,
            "aggregation": aggregation,
            "half_deg": half_deg,
            "aot_mean": float(np.mean(aot_values)),
            "aot_median": float(np.median(aot_values)),
            "aot_std": float(np.std(aot_values)),
            "unc_nat": native_uncertainty,
            "n_good": int(aot_values.size),
            "day": day.isoformat(),
            "aot_unc": calibrated_uncertainty(aot, native_uncertainty, int(aot_values.size)),
            "searched_granules": len(paths),
        }
    _atomic_json(output_dir / f"{matchup_id}.json", record)
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matchup_id", nargs="+")
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--cache-dir", type=Path)
    parser.add_argument("--half-deg", type=float, default=DEFAULT_HALF_DEG)
    parser.add_argument("--aggregation", choices=("median", "mean"), default="median")
    args = parser.parse_args()
    output_dir = args.output_dir or args.root / "maiac_qa_native"
    cache_dir = args.cache_dir or args.root / "cache" / "mcd19_native_stage"
    with (args.root / "matchups" / "matchups.csv").open(encoding="utf-8", newline="") as stream:
        matchups = {row["matchup_id"]: row for row in csv.DictReader(stream)}
    for matchup_id in args.matchup_id:
        record = stage_one(
            matchups[matchup_id],
            output_dir=output_dir,
            cache_dir=cache_dir,
            half_deg=args.half_deg,
            aggregation=args.aggregation,
        )
        print(
            f"MAIAC_NATIVE {matchup_id} status={record['status']} "
            f"aot={record.get('aot')} n={record.get('n_good', 0)}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
