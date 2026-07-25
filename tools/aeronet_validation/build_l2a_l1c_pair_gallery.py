"""Render a small spatial gallery of exact same-day L2A/current-RT pairs.

The pair archive contains sampled clear-land pixels for fitting. This tool
re-fetches four selected acquisitions on the original 60 m grid so visual
inspection uses the same L2A product, L1C radiometry, saved MAIAC-conditioned
RT coefficients, and clear-land mask as the archived pairs.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
from tools.aeronet_validation.build_l2a_l1c_harmonization_pairs import (
    ELEVATION_CACHE,
    HALF_SIZE_DEGREES,
    MATCHUPS,
    ROOT,
    _load_pair_grids,
    _site_elevation_km,
)

PAIR_DIR = ROOT / "analysis/l2a_l1c_exact_pairs_mediumdev_20260713"
DEFAULT_OUTPUT = ROOT / "reports/l2a-harmonization-20260714/assets/pair-examples"
RGB_BANDS = (3, 2, 1)
VISIBLE_BANDS = (1, 2, 3)
EXAMPLES: tuple[tuple[str, str, str], ...] = (
    (
        "Medenine-IRA__T32SPC_20240613T100031",
        "S2B_32SPC_20200520_0_L1C",
        "arid North African plain",
    ),
    (
        "Windpoort__T33KWU_20241005T084729",
        "S2B_33KWU_20200906_0_L1C",
        "semi-arid inland vegetation",
    ),
    (
        "Beijing-CAMS__T50TMK_20240104T030119",
        "S2B_50TMK_20200115_1_L1C",
        "winter urban and peri-urban terrain",
    ),
    (
        "Mexico_City__T14QMG_20240331T165851",
        "S2A_14QMG_20200411_0_L1C",
        "highland urban basin",
    ),
)


def _scene_from_archive(matchup_id: str, scene_id: str) -> dict[str, Any]:
    archive_path = PAIR_DIR / f"{matchup_id}.npz"
    with np.load(archive_path, allow_pickle=False) as archive:
        scenes = json.loads(str(archive["scenes_json"]))
    for scene in scenes:
        if str(scene["scene_id"]) == scene_id:
            return scene
    raise KeyError(f"{scene_id} is not retained in {archive_path.name}")


def _pair_scene_input(scene: dict[str, Any]) -> dict[str, Any]:
    return {
        "l1c_id": scene["scene_id"],
        "day": scene["day"],
        "tile": scene["tile"],
        "maiac": scene["maiac_aot"],
        "wvp": scene["maiac_tcwv_cm"],
        "sza": scene["sza_deg"],
        "saa": scene["saa_deg"],
        "vza": scene["vza_deg"],
        "vaa": scene["vaa_deg"],
    }


def _scaled_rgb(
    l2a_surface: np.ndarray,
    siac_surface: np.ndarray,
    valid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    l2a_rgb = np.asarray(l2a_surface, dtype=np.float32)[list(RGB_BANDS)]
    siac_rgb = np.asarray(siac_surface, dtype=np.float32)[list(RGB_BANDS)]
    paired_values = np.concatenate(
        (l2a_rgb[:, valid], siac_rgb[:, valid]), axis=1
    )
    lower = np.nanpercentile(paired_values, 1.0, axis=1)
    upper = np.nanpercentile(paired_values, 99.0, axis=1)
    scale = np.maximum(upper - lower, 1.0e-5)

    def rescale(values: np.ndarray) -> np.ndarray:
        image = np.moveaxis(values, 0, -1)
        image = np.clip((image - lower) / scale, 0.0, 1.0)
        image = np.power(image, 1.0 / 2.2)
        image[~valid] = np.nan
        return image

    return rescale(l2a_rgb), rescale(siac_rgb)


def _visible_difference(
    l2a_surface: np.ndarray,
    siac_surface: np.ndarray,
    valid: np.ndarray,
) -> np.ndarray:
    difference = np.full(valid.shape, np.nan, dtype=np.float32)
    signed = (
        np.asarray(l2a_surface, dtype=np.float32)[list(VISIBLE_BANDS)]
        - np.asarray(siac_surface, dtype=np.float32)[list(VISIBLE_BANDS)]
    )
    difference[valid] = np.mean(signed[:, valid], axis=0)
    return difference


def _plot_row(
    axes: np.ndarray,
    *,
    l2a_rgb: np.ndarray,
    siac_rgb: np.ndarray,
    difference: np.ndarray,
    label: str,
    metadata: dict[str, Any],
) -> Any:
    axes[0].imshow(l2a_rgb, interpolation="nearest")
    axes[1].imshow(siac_rgb, interpolation="nearest")
    image = axes[2].imshow(
        difference,
        cmap="RdBu_r",
        vmin=-0.04,
        vmax=0.04,
        interpolation="nearest",
    )
    for axis in axes:
        axis.set_xticks([])
        axis.set_yticks([])
    axes[0].set_ylabel(
        f"{label}\n{metadata['day']} | MAIAC {metadata['maiac_aot']:.3f}",
        fontsize=9,
        rotation=0,
        ha="right",
        va="center",
    )
    return image


def _example_record(
    ee: Any,
    *,
    matchup_id: str,
    scene_id: str,
    setting: str,
    rows: dict[str, dict[str, str]],
    elevations: dict[str, float],
) -> tuple[dict[str, Any], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    saved_scene = _scene_from_archive(matchup_id, scene_id)
    row = rows[matchup_id]
    lon = float(row["longitude"])
    lat = float(row["latitude"])
    bbox = (
        lon - HALF_SIZE_DEGREES,
        lat - HALF_SIZE_DEGREES,
        lon + HALF_SIZE_DEGREES,
        lat + HALF_SIZE_DEGREES,
    )
    grids, metadata = _load_pair_grids(
        ee,
        bbox=bbox,
        sidecar_path=Path(saved_scene["sidecar"]),
        scene=_pair_scene_input(saved_scene),
        elevation_km=_site_elevation_km(row, elevations),
    )
    valid = np.asarray(grids["valid"], dtype=bool)
    l2a_surface = np.asarray(grids["l2a_surface"], dtype=np.float32)
    siac_surface = np.asarray(grids["siac_surface"], dtype=np.float32)
    l2a_rgb, siac_rgb = _scaled_rgb(l2a_surface, siac_surface, valid)
    difference = _visible_difference(l2a_surface, siac_surface, valid)
    b8a_difference = l2a_surface[4] - siac_surface[4]
    visible_values = difference[valid]
    b8a_values = b8a_difference[valid]
    return (
        {
            "matchup_id": matchup_id,
            "site": row["site"],
            "setting": setting,
            "scene_id": scene_id,
            "day": metadata["day"],
            "maiac_aot": metadata["maiac_aot"],
            "maiac_tcwv_cm": metadata["maiac_tcwv_cm"],
            "l2a_aot_median": float(np.nanmedian(np.asarray(grids["l2a_aot"])[valid])),
            "valid_pixels": int(np.count_nonzero(valid)),
            "valid_fraction": float(np.mean(valid)),
            "visible_delta_mean": float(np.nanmean(visible_values)),
            "visible_delta_median": float(np.nanmedian(visible_values)),
            "visible_delta_mae": float(np.nanmean(np.abs(visible_values))),
            "b8a_delta_median": float(np.nanmedian(b8a_values)),
            "grid": metadata["grid"],
            "l2a": metadata["l2a"],
            "l1c_asset": metadata["l1c_asset"],
            "cloud_score_asset": metadata["cloud_score_asset"],
            "sidecar": metadata["sidecar"],
        },
        (l2a_rgb, siac_rgb, difference),
    )


def _render_single(
    path: Path,
    *,
    record: dict[str, Any],
    arrays: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 3, figsize=(15.5, 5.6), facecolor="white")
    image = _plot_row(
        axes,
        l2a_rgb=arrays[0],
        siac_rgb=arrays[1],
        difference=arrays[2],
        label=f"{record['site']} ({record['setting']})",
        metadata=record,
    )
    for axis, title in zip(
        axes,
        (
            "Operational L2A surface reflectance",
            "Same-day L1C corrected with MAIAC/current-RT",
            "Mean B02-B04 delta: L2A - current-RT",
        ),
        strict=True,
    ):
        axis.set_title(title, fontsize=10)
    colorbar = figure.colorbar(image, ax=axes[2], fraction=0.046, pad=0.035)
    colorbar.set_label("Surface reflectance", fontsize=9)
    figure.suptitle(
        f"Exact same-day pair: {record['scene_id']}", x=0.03, ha="left", fontsize=13
    )
    figure.subplots_adjust(left=0.16, right=0.97, top=0.88, bottom=0.08, wspace=0.05)
    figure.savefig(path, dpi=155, bbox_inches="tight")
    plt.close(figure)


def _render_gallery(
    path: Path,
    examples: list[tuple[dict[str, Any], tuple[np.ndarray, np.ndarray, np.ndarray]]],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(len(examples), 3, figsize=(15.5, 18.6), facecolor="white")
    last_image = None
    for row_axes, (record, arrays) in zip(axes, examples, strict=True):
        last_image = _plot_row(
            row_axes,
            l2a_rgb=arrays[0],
            siac_rgb=arrays[1],
            difference=arrays[2],
            label=f"{record['site']} ({record['setting']})",
            metadata=record,
        )
    for axis, title in zip(
        axes[0],
        (
            "Operational L2A surface reflectance",
            "Same-day L1C corrected with MAIAC/current-RT",
            "Mean B02-B04 delta: L2A - current-RT",
        ),
        strict=True,
    ):
        axis.set_title(title, fontsize=11, pad=10)
    assert last_image is not None
    colorbar = figure.colorbar(last_image, ax=axes[:, 2], fraction=0.025, pad=0.025)
    colorbar.set_label("Surface reflectance", fontsize=9)
    figure.suptitle(
        "Selected exact same-day L2A and current-RT spatial pairs",
        x=0.04,
        y=0.99,
        ha="left",
        fontsize=15,
    )
    figure.text(
        0.04,
        0.968,
        "Each row uses its own paired RGB contrast stretch; blank pixels fail the same clear-land mask used for training.",
        ha="left",
        fontsize=9,
    )
    figure.subplots_adjust(left=0.20, right=0.93, top=0.94, bottom=0.035, hspace=0.16, wspace=0.045)
    figure.savefig(path, dpi=155, bbox_inches="tight")
    plt.close(figure)


def build(output: Path) -> dict[str, Any]:
    from bestpixel._gee import init_ee

    output.mkdir(parents=True, exist_ok=True)
    rows = {row["matchup_id"]: row for row in csv.DictReader(MATCHUPS.open())}
    elevations = json.loads(ELEVATION_CACHE.read_text(encoding="utf-8"))
    ee = init_ee()
    rendered: list[tuple[dict[str, Any], tuple[np.ndarray, np.ndarray, np.ndarray]]] = []
    for matchup_id, scene_id, setting in EXAMPLES:
        record, arrays = _example_record(
            ee,
            matchup_id=matchup_id,
            scene_id=scene_id,
            setting=setting,
            rows=rows,
            elevations=elevations,
        )
        filename = f"{matchup_id}__{scene_id}.png"
        _render_single(output / filename, record=record, arrays=arrays)
        record["image"] = filename
        rendered.append((record, arrays))
        print(
            f"PAIR_GALLERY {record['site']} {record['day']} "
            f"valid={record['valid_pixels']} aod={record['maiac_aot']:.3f}",
            flush=True,
        )
    gallery_path = output / "pair-gallery.png"
    _render_gallery(gallery_path, rendered)
    manifest = {
        "title": "Selected exact same-day L2A/current-RT spatial pairs",
        "selection": "Four representative settings were selected before rendering for visual diversity, not by AOD retrieval outcome.",
        "source_definition": "Operational COPERNICUS/S2_SR_HARMONIZED L2A paired on the same day and 60 m UTM grid with L1C corrected using the saved MAIAC-conditioned current libRadTran LUT coefficients.",
        "mask_definition": "Cloud Score+ cs > 0.60, SCL vegetation or bare soil, finite in all seven bands and L2A AOT/WVP, reflectance within (0.001, 0.8), followed by one-pixel binary erosion.",
        "rgb_definition": "B04/B03/B02, with a paired per-row 1st-99th percentile contrast stretch shared by L2A and current-RT.",
        "difference_definition": "Mean B02/B03/B04 surface reflectance: L2A minus current-RT. Red is positive and blue is negative; range is fixed at +/-0.04.",
        "gallery_image": "assets/pair-examples/pair-gallery.png",
        "examples": [record for record, _arrays in rendered],
    }
    (output / "pair-gallery-metadata.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(build(args.output), indent=2))


if __name__ == "__main__":
    main()
