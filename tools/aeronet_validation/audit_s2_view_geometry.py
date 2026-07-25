"""Audit Sentinel-2 viewing geometry on the fixed low-cloud AERONET cohort."""

from __future__ import annotations

import argparse
import csv
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import numpy as np
from pyproj import Transformer

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
COHORT = ROOT / "campaign250_lowcloud20_mids.txt"
RESULTS = ROOT / "phaseD_results_lowcloud20_singleprior_b03_chi2_20260711"
SAFE_ROOT = ROOT / "phaseD_cache/s2"
SOLVE_BAND_IDS = (1, 2, 3)  # B02, B03, B04


def _local_name(element: ET.Element) -> str:
    return element.tag.rsplit("}", 1)[-1]


def _descendants(element: ET.Element, name: str) -> list[ET.Element]:
    return [child for child in element.iter() if _local_name(child) == name]


def _first_descendant(element: ET.Element, name: str) -> ET.Element:
    matches = _descendants(element, name)
    if not matches:
        raise ValueError(f"Missing {name} in Sentinel-2 tile metadata")
    return matches[0]


def _child(element: ET.Element, name: str) -> ET.Element:
    for candidate in element:
        if _local_name(candidate) == name:
            return candidate
    raise ValueError(f"Missing {name} below {_local_name(element)}")


def _float_child(element: ET.Element, name: str) -> float:
    child = _child(element, name)
    if child.text is None:
        raise ValueError(f"Empty {name} below {_local_name(element)}")
    return float(child.text)


def _angle_grid(element: ET.Element, angle_name: str) -> np.ndarray:
    angle = _child(element, angle_name)
    values_list = _child(angle, "Values_List")
    rows = []
    for row in values_list:
        if _local_name(row) == "VALUES" and row.text:
            rows.append([float(value) for value in row.text.split()])
    if not rows:
        raise ValueError(f"Empty {angle_name} grid")
    return np.asarray(rows, dtype=np.float64)


def _circular_mean_deg(values: np.ndarray, axis: int | tuple[int, ...] | None = None) -> Any:
    radians = np.deg2rad(values)
    finite = np.isfinite(radians)
    count = np.sum(finite, axis=axis)
    sine = np.sum(np.where(finite, np.sin(radians), 0.0), axis=axis)
    cosine = np.sum(np.where(finite, np.cos(radians), 0.0), axis=axis)
    angle = np.mod(np.rad2deg(np.arctan2(sine, cosine)), 360.0)
    return np.where(count > 0, angle, np.nan)


def _detector_composite(root: ET.Element, band_id: int) -> tuple[np.ndarray, np.ndarray]:
    zenith_grids = []
    azimuth_grids = []
    for element in _descendants(root, "Viewing_Incidence_Angles_Grids"):
        if int(element.attrib["bandId"]) != band_id:
            continue
        zenith_grids.append(_angle_grid(element, "Zenith"))
        azimuth_grids.append(_angle_grid(element, "Azimuth"))
    if not zenith_grids:
        raise ValueError(f"No detector geometry grids for bandId={band_id}")
    zenith_stack = np.stack(zenith_grids)
    azimuth_stack = np.stack(azimuth_grids)
    with np.errstate(invalid="ignore"):
        zenith_count = np.sum(np.isfinite(zenith_stack), axis=0)
        zenith = np.divide(
            np.nansum(zenith_stack, axis=0),
            zenith_count,
            out=np.full(zenith_count.shape, np.nan, dtype=np.float64),
            where=zenith_count > 0,
        )
        azimuth = _circular_mean_deg(azimuth_stack, axis=0)
    return zenith, np.asarray(azimuth, dtype=np.float64)


def _bilinear(grid: np.ndarray, row: float, col: float) -> float:
    row = float(np.clip(row, 0.0, grid.shape[0] - 1.0))
    col = float(np.clip(col, 0.0, grid.shape[1] - 1.0))
    r0 = int(math.floor(row))
    c0 = int(math.floor(col))
    r1 = min(r0 + 1, grid.shape[0] - 1)
    c1 = min(c0 + 1, grid.shape[1] - 1)
    samples = np.asarray([grid[r0, c0], grid[r0, c1], grid[r1, c0], grid[r1, c1]])
    weights = np.asarray(
        [
            (1.0 - (row - r0)) * (1.0 - (col - c0)),
            (1.0 - (row - r0)) * (col - c0),
            (row - r0) * (1.0 - (col - c0)),
            (row - r0) * (col - c0),
        ]
    )
    valid = np.isfinite(samples)
    if np.any(valid):
        return float(np.sum(samples[valid] * weights[valid]) / np.sum(weights[valid]))
    finite = np.argwhere(np.isfinite(grid))
    if finite.size == 0:
        return float("nan")
    nearest = finite[np.argmin((finite[:, 0] - row) ** 2 + (finite[:, 1] - col) ** 2)]
    return float(grid[tuple(nearest)])


def _bilinear_azimuth(grid: np.ndarray, row: float, col: float) -> float:
    sine = _bilinear(np.sin(np.deg2rad(grid)), row, col)
    cosine = _bilinear(np.cos(np.deg2rad(grid)), row, col)
    return float(np.mod(np.rad2deg(np.arctan2(sine, cosine)), 360.0))


def _angular_difference_deg(candidate: float, reference: float) -> float:
    return float((candidate - reference + 180.0) % 360.0 - 180.0)


def _folded_raa(vaa: float, saa: float) -> float:
    separation = abs(_angular_difference_deg(vaa, saa))
    return float(min(separation, 360.0 - separation))


def _mean_view_angles(root: ET.Element) -> dict[int, tuple[float, float]]:
    values: dict[int, tuple[float, float]] = {}
    for element in _descendants(root, "Mean_Viewing_Incidence_Angle"):
        band_id = int(element.attrib["bandId"])
        values[band_id] = (
            _float_child(element, "ZENITH_ANGLE"),
            _float_child(element, "AZIMUTH_ANGLE"),
        )
    return values


def _grid_position(root: ET.Element, lon: float, lat: float) -> tuple[float, float]:
    crs = _first_descendant(root, "HORIZONTAL_CS_CODE").text
    if crs is None:
        raise ValueError("Empty HORIZONTAL_CS_CODE")
    geoposition = next(
        element
        for element in _descendants(root, "Geoposition")
        if element.attrib.get("resolution") == "10"
    )
    ulx = _float_child(geoposition, "ULX")
    uly = _float_child(geoposition, "ULY")
    xdim = _float_child(geoposition, "XDIM")
    ydim = _float_child(geoposition, "YDIM")
    transformer = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    x, y = transformer.transform(lon, lat)
    angle_step = 5000.0
    col = (x - ulx) / angle_step
    row = (y - uly) / (-angle_step if ydim < 0.0 else angle_step)
    if xdim <= 0.0:
        raise ValueError(f"Unexpected Sentinel-2 XDIM={xdim}")
    return row, col


def _case_row(record: dict[str, Any]) -> dict[str, Any]:
    safe = SAFE_ROOT / f"{record['product_id']}.SAFE"
    xml_paths = list(safe.glob("GRANULE/*/MTD_TL.xml"))
    if len(xml_paths) != 1:
        raise ValueError(f"Expected one MTD_TL.xml for {safe}, found {len(xml_paths)}")
    root = ET.parse(xml_paths[0]).getroot()
    row, col = _grid_position(root, float(record["lon"]), float(record["lat"]))

    mean_angles = _mean_view_angles(root)
    current_vza = float(np.mean([value[0] for value in mean_angles.values()]))
    current_vaa = float(np.mean([value[1] for value in mean_angles.values()]))
    solve_mean_vza = float(np.mean([mean_angles[band_id][0] for band_id in SOLVE_BAND_IDS]))
    solve_mean_vaa = float(
        _circular_mean_deg(np.asarray([mean_angles[band_id][1] for band_id in SOLVE_BAND_IDS]))
    )

    sun = _first_descendant(root, "Sun_Angles_Grid")
    station_sza = _bilinear(_angle_grid(sun, "Zenith"), row, col)
    station_saa = _bilinear_azimuth(_angle_grid(sun, "Azimuth"), row, col)

    band_vza: dict[int, float] = {}
    band_vaa: dict[int, float] = {}
    for band_id in SOLVE_BAND_IDS:
        zenith, azimuth = _detector_composite(root, band_id)
        band_vza[band_id] = _bilinear(zenith, row, col)
        band_vaa[band_id] = _bilinear_azimuth(azimuth, row, col)
    grid_vza = float(np.mean(list(band_vza.values())))
    grid_vaa = float(_circular_mean_deg(np.asarray(list(band_vaa.values()))))
    current_raa = _folded_raa(current_vaa, station_saa)
    grid_raa = _folded_raa(grid_vaa, station_saa)

    truth = float(record["truth"])
    retrieved = float(record["retrieved"])
    return {
        "matchup_id": record["matchup_id"],
        "site": record["site"],
        "product_id": record["product_id"],
        "truth": truth,
        "retrieved": retrieved,
        "error": retrieved - truth,
        "within_ee": bool(record["within_ee"]),
        "station_sza_deg": station_sza,
        "station_saa_deg": station_saa,
        "current_all_band_mean_vza_deg": current_vza,
        "current_all_band_mean_vaa_deg": current_vaa,
        "solve_band_tile_mean_vza_deg": solve_mean_vza,
        "solve_band_tile_mean_vaa_deg": solve_mean_vaa,
        "solve_band_station_grid_vza_deg": grid_vza,
        "solve_band_station_grid_vaa_deg": grid_vaa,
        "current_raa_deg": current_raa,
        "station_grid_raa_deg": grid_raa,
        "station_minus_current_vza_deg": grid_vza - current_vza,
        "station_minus_current_vaa_deg": _angular_difference_deg(grid_vaa, current_vaa),
        "station_minus_current_raa_deg": grid_raa - current_raa,
        "B02_station_vza_deg": band_vza[1],
        "B03_station_vza_deg": band_vza[2],
        "B04_station_vza_deg": band_vza[3],
        "B02_station_vaa_deg": band_vaa[1],
        "B03_station_vaa_deg": band_vaa[2],
        "B04_station_vaa_deg": band_vaa[3],
    }


def _summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    def values(name: str) -> np.ndarray:
        return np.asarray([float(row[name]) for row in rows], dtype=np.float64)

    errors = values("error")
    vza_delta = values("station_minus_current_vza_deg")
    vaa_delta = values("station_minus_current_vaa_deg")
    raa_delta = values("station_minus_current_raa_deg")
    miss = np.asarray([not bool(row["within_ee"]) for row in rows])
    return {
        "cases": len(rows),
        "hits": int(np.sum(~miss)),
        "vza_delta_deg": _distribution(vza_delta),
        "absolute_vaa_delta_deg": _distribution(np.abs(vaa_delta)),
        "absolute_raa_delta_deg": _distribution(np.abs(raa_delta)),
        "vza_delta_error_correlation": _correlation(vza_delta, errors),
        "raa_delta_error_correlation": _correlation(raa_delta, errors),
        "median_absolute_raa_delta_hits_deg": float(np.median(np.abs(raa_delta[~miss]))),
        "median_absolute_raa_delta_misses_deg": float(np.median(np.abs(raa_delta[miss]))),
        "cases_abs_raa_delta_gt_5deg": int(np.sum(np.abs(raa_delta) > 5.0)),
        "cases_abs_vza_delta_gt_1deg": int(np.sum(np.abs(vza_delta) > 1.0)),
    }


def _distribution(values: np.ndarray) -> dict[str, float]:
    return {
        "min": float(np.min(values)),
        "p05": float(np.quantile(values, 0.05)),
        "median": float(np.median(values)),
        "p95": float(np.quantile(values, 0.95)),
        "max": float(np.max(values)),
        "mean": float(np.mean(values)),
    }


def _correlation(left: np.ndarray, right: np.ndarray) -> float | None:
    if np.std(left) == 0.0 or np.std(right) == 0.0:
        return None
    return float(np.corrcoef(left, right)[0, 1])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    matchup_ids = [line.strip() for line in COHORT.read_text().splitlines() if line.strip()]
    records = [
        json.loads((RESULTS / f"{matchup_id}.json").read_text()) for matchup_id in matchup_ids
    ]
    rows = [_case_row(record) for record in records]
    summary = _summary(rows)
    print(json.dumps(summary, indent=2, sort_keys=True))

    if args.output_dir is not None:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        with (args.output_dir / "cases.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    main()
