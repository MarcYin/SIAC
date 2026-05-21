"""Refresh the captured goldens for a regression scene.

Usage::

    PYTHONPATH=python pixi run -e rt6s python \
        tests/regression/regenerate_goldens.py \
        --output-dir /Users/fengyin/siac_runs/output_t33kwp_sixs \
        --golden tests/regression/goldens/t33kwp_sixs_20260329.json \
        --scene-id S2B_MSIL1C_20260329T084559_N0512_R107_T33KWP_20260329T140503

The script reads the *AOT*, *TCWV*, *CLOUD*, and per-band *BOA_*\\* COG
outputs from a completed pipeline run plus the STAC sidecar JSON, and
writes a normalised goldens file. Run this **after** an intentional
numerical change has been validated by other means (e.g. visual review,
comparison against a reference scene). Do not use it to silence a real
regression — the regression test exists to catch those.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

# Standard set of STAC properties we lock in the goldens. Add to this list
# (and to the test) when a new scene-stable property becomes worth checking.
_STAC_KEYS_TO_CAPTURE: tuple[str, ...] = (
    "siac:aot_mean",
    "siac:tcwv_mean",
    "siac:satellite",
    "siac:tile_id",
    "siac:sensor",
    "siac:processing_baseline",
    "platform",
    "constellation",
    "view:sun_elevation",
    "view:sun_azimuth",
    "view:off_nadir",
    "proj:epsg",
    "datetime",
)

# Sentinel-2 band layout. Edit if a new scene captures a different sensor.
_S2_BANDS: tuple[str, ...] = (
    "B01",
    "B02",
    "B03",
    "B04",
    "B05",
    "B06",
    "B07",
    "B08",
    "B8A",
    "B09",
    "B10",
    "B11",
    "B12",
)


def _resolve_prefix(out_dir: Path) -> str:
    matches = sorted(out_dir.glob("*_AOT.tif"))
    if not matches:
        raise FileNotFoundError(f"No *_AOT.tif under {out_dir}")
    if len(matches) > 1:
        raise RuntimeError(f"Multiple *_AOT.tif under {out_dir}: {matches}")
    return matches[0].stem.removesuffix("_AOT")


def _stat_dict(path: Path) -> dict[str, Any]:
    """Read a single-band raster and return the same stats dict tests use."""
    # Imported here so the script doesn't pay rasterio startup unless invoked.
    import numpy as np
    import rasterio

    with rasterio.open(path) as ds:
        arr = ds.read(1, masked=True)
        finite = np.ma.compressed(arr)
        if finite.size == 0:
            return {"shape": list(arr.shape), "all_masked": True}
        return {
            "shape": list(arr.shape),
            "dtype": str(arr.dtype),
            "valid_fraction": float(finite.size / arr.size),
            "mean": float(np.mean(finite)),
            "std": float(np.std(finite)),
            "p01": float(np.percentile(finite, 1)),
            "p50": float(np.percentile(finite, 50)),
            "p99": float(np.percentile(finite, 99)),
            "min": float(np.min(finite)),
            "max": float(np.max(finite)),
        }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory containing the SIAC pipeline outputs (AOT/TCWV/BOA_*).",
    )
    parser.add_argument(
        "--golden",
        required=True,
        type=Path,
        help="Path to the goldens JSON to (over)write.",
    )
    parser.add_argument(
        "--scene-id",
        required=True,
        help="The Sentinel-2 product id (used only as metadata in the goldens).",
    )
    parser.add_argument(
        "--config-path",
        default=None,
        help="Optional: relative path to the config TOML used (kept as metadata).",
    )
    parser.add_argument(
        "--rt-backend",
        default="sixs",
        choices=("sixs", "lut", "emulator"),
        help="Which RT backend was used for the run.",
    )
    args = parser.parse_args(argv)

    out_dir = args.output_dir
    if not out_dir.is_dir():
        parser.error(f"--output-dir does not exist: {out_dir}")
    golden_path = args.golden

    prefix = _resolve_prefix(out_dir)
    products: dict[str, dict[str, Any]] = {}
    for label, file_suffix in (
        ("AOT", "AOT.tif"),
        ("TCWV", "TCWV.tif"),
        ("CLOUD", "CLOUD.tif"),
    ):
        path = out_dir / f"{prefix}_{file_suffix}"
        if path.exists():
            products[label] = _stat_dict(path)
    for band in _S2_BANDS:
        path = out_dir / f"{prefix}_BOA_{band}.tif"
        if path.exists():
            products[f"BOA_{band}"] = _stat_dict(path)

    stac_path = out_dir / f"{prefix}.json"
    if not stac_path.exists():
        parser.error(f"Missing STAC sidecar {stac_path}")
    stac = json.loads(stac_path.read_text())
    selected_props = {
        k: stac["properties"][k] for k in _STAC_KEYS_TO_CAPTURE if k in stac["properties"]
    }

    # Read existing meta (if any) to preserve non-numerical context fields.
    existing_meta: dict[str, Any] = {}
    if golden_path.exists():
        with golden_path.open() as fh:
            existing_meta = json.load(fh).get("_meta", {})
    new_meta = {
        **existing_meta,
        "scene_id": args.scene_id,
        "rt_backend": args.rt_backend,
        "captured_from": str(out_dir),
        "captured_at_utc": _utcnow_iso(),
    }
    if args.config_path is not None:
        new_meta["config_path_relative_to_repo"] = args.config_path

    doc = {
        "_meta": new_meta,
        "products": products,
        "stac": selected_props,
        "bbox": stac["bbox"],
    }

    golden_path.parent.mkdir(parents=True, exist_ok=True)
    with golden_path.open("w") as fh:
        json.dump(doc, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"Wrote {golden_path} ({golden_path.stat().st_size} bytes)")
    print(f"  products: {len(products)} ({sorted(products)})")
    print(f"  stac    : {len(selected_props)} keys")
    return 0


def _utcnow_iso() -> str:
    from datetime import datetime, timezone

    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


if __name__ == "__main__":
    sys.exit(main())
