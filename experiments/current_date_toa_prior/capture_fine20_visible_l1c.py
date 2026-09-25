#!/usr/bin/env python3
"""Add exact-product 20 m current-date L1C visible inputs to the mixed archive.

The original mixed archive already contains native 20 m red-edge/NIR/SWIR
planes and a 60 m wide-context tensor.  This companion capture reads only the
missing B01--B04 planes on the existing ``detail20`` grid, avoiding a second
download of data that are already frozen.  Reads remain windowed against the
anonymous Google Cloud Sentinel-2 public bucket; no SAFE is downloaded.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from experiments.current_date_toa_prior.capture_mixed_resolution_l1c import (
    CAPTURE_RECEIPT_SCHEMA,
    MIXED_SCHEMA,
    _all_band_finite_fraction,
    _radiometric_calibration,
    _read_multitile_bands,
)

FINE20_VISIBLE_BANDS = ("B01", "B02", "B03", "B04")
FINE20_SCHEMA = "siac_l1c_fine20_visible_v3"


def _scalar(value: np.ndarray) -> Any:
    return np.asarray(value).item()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_one(matchup_id: str, args: argparse.Namespace) -> dict[str, Any]:
    source = args.current_root / f"{matchup_id}.npz"
    target = args.output_root / f"{matchup_id}.npz"
    receipt_path = args.output_root / f"{matchup_id}.receipt.json"
    if target.is_file() and receipt_path.is_file() and not args.force:
        try:
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            domain = receipt.get("reflectance_domain") or {}
            with np.load(target, allow_pickle=False) as existing:
                schema = str(_scalar(existing["schema_version"]))
                existing_bands = existing["fine20_bands"] if "fine20_bands" in existing.files else ()
            if (
                schema == FINE20_SCHEMA
                and tuple(str(b) for b in existing_bands) == tuple(getattr(args, "bands", FINE20_VISIBLE_BANDS))
                and receipt.get("schema_version") == CAPTURE_RECEIPT_SCHEMA
                and domain.get("maximum") is None
                and domain.get("upper_clipping") is False
                and (receipt.get("radiometric_calibration") or {}).get("policy")
                == "exact_product_mtd_msil1c"
            ):
                return {"matchup_id": matchup_id, "status": "exists", "output": str(target)}
        except (OSError, ValueError, KeyError, EOFError, json.JSONDecodeError):
            pass
    if not source.is_file():
        raise FileNotFoundError(source)
    with np.load(source, allow_pickle=False) as payload:
        schema = str(_scalar(payload["schema_version"]))
        if schema != MIXED_SCHEMA:
            raise ValueError(f"{matchup_id}: unsupported source schema {schema!r}")
        if str(_scalar(payload["matchup_id"])) != matchup_id:
            raise ValueError(f"{matchup_id}: source matchup mismatch")
        product_id = str(_scalar(payload["product_id"]))
        shape = tuple(int(value) for value in np.asarray(payload["detail20_shape"]).ravel())
        transform = tuple(
            float(value) for value in np.asarray(payload["detail20_transform"]).ravel()[:6]
        )
        crs = str(_scalar(payload["crs"]))
    if len(shape) != 2 or len(transform) != 6:
        raise ValueError(f"{matchup_id}: malformed detail20 grid")

    bands = tuple(getattr(args, "bands", FINE20_VISIBLE_BANDS))
    planes, records, products = _read_multitile_bands(
        product_id,
        bands,
        crs=crs,
        transform=transform,
        height=shape[0],
        width=shape[1],
        minimum_finite_fraction=float(args.minimum_finite_fraction),
    )
    finite_fraction = _all_band_finite_fraction(planes)
    if finite_fraction < float(args.minimum_source_fraction):
        raise RuntimeError(
            f"{matchup_id}: fine20 visible coverage {finite_fraction:.4f} is below "
            f"{float(args.minimum_source_fraction):.4f}"
        )

    args.output_root.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.stem}.{os.getpid()}.partial.npz")
    try:
        with temporary.open("wb") as stream:
            np.savez_compressed(
                stream,
                schema_version=np.asarray(FINE20_SCHEMA),
                matchup_id=np.asarray(matchup_id),
                product_id=np.asarray(product_id),
                fine20_shape=np.asarray(shape, dtype=np.int64),
                fine20_transform=np.asarray(transform, dtype=np.float64),
                crs=np.asarray(crs),
                fine20_bands=np.asarray(bands),
                **{f"fine20_{band}": value for band, value in planes.items()},
            )
        temporary.replace(target)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    receipt = {
        "schema_version": CAPTURE_RECEIPT_SCHEMA,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "matchup_id": matchup_id,
        "product_id": product_id,
        "source_archive": {"path": str(source), "sha256": _sha256(source)},
        "source": "anonymous Google Cloud public Sentinel-2 L1C windowed JP2 reads",
        "reflectance_domain": {
            "minimum": 0.0,
            "maximum": None,
            "upper_clipping": False,
        },
        "radiometric_calibration": _radiometric_calibration(records),
        "source_products": products,
        "grid": {"shape": list(shape), "transform": list(transform), "crs": crs},
        "bands": records,
        "quality": {"fine20_visible_all_bands_finite_fraction": finite_fraction},
        "output": {"path": str(target), "bytes": target.stat().st_size, "sha256": _sha256(target)},
    }
    receipt_path.write_text(json.dumps(receipt, indent=2, allow_nan=False), encoding="utf-8")
    return {"matchup_id": matchup_id, "status": "ok", "output": str(target), **receipt["quality"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matchup_ids", nargs="+")
    parser.add_argument("--current-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--minimum-finite-fraction", type=float, default=0.98)
    parser.add_argument("--minimum-source-fraction", type=float, default=0.50)
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--bands",
        type=lambda text: tuple(name.strip() for name in text.split(",") if name.strip()),
        default=FINE20_VISIBLE_BANDS,
        help="bands to read on the detail20 grid (default the visible B01-B04; e.g. B08)",
    )
    args = parser.parse_args()
    lock_root = args.output_root / ".locks"
    lock_root.mkdir(parents=True, exist_ok=True)
    for matchup_id in args.matchup_ids:
        lock_path = lock_root / f"{matchup_id}.lock"
        with lock_path.open("a+", encoding="utf-8") as lock:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
            result = run_one(matchup_id, args)
        print(json.dumps(result, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
