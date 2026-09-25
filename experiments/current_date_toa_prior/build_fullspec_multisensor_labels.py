"""Same-date full-spectrum supervision for the S2 half of the joint Landsat/S2 corpus.

Extends ``build_hybrid_multisensor_labels`` from the visible triple to every land band.
For each S2 sample a dense query block is drawn from the acquisition's own
full-spectrum M5 teacher (``fullspec_training_labels_20260925/teacher_fullspec``):

* B01-B04: the ExtraTrees surface prediction (the label's visible construction);
* B05, B06, B07, B08, B8A, B11, B12: the 6S correction of the scene's own TOA at the
  M5-solved AOD (``extra_anchor_boa_at_solution`` / ``anchor_boa_at_solution``), which
  beats the tree readout of those bands against the AERONET measured-aerosol reference
  (B05/B06/B11/B12 better, B07/B8A tied).

The month-matched library labels on those bands are masked off, so the library block
only keeps its role in the native RT closure. Landsat samples are left untouched.
Uncertainty is the teacher's per-band sigma; eligibility is the same-date teacher's QA
weight. The S2B/S2C -> S2A constant offset is the measured visible one on B02-B04 and
zero elsewhere (not yet measured beyond the visible; recorded in the archive).

No AERONET measurement is read; the labels are model output, not truth.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import traceback
from pathlib import Path

import numpy as np
from experiments.current_date_toa_prior.build_hybrid_multisensor_labels import (
    BANDS,
    SCHEMA,
    _extend,
    digest,
    visible_correction,
)

LAND = ("B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12")
TREE_BANDS = ("B01", "B02", "B03", "B04")
SIXS_BANDS = ("B05", "B06", "B07", "B08", "B8A", "B11", "B12")
TARGET = [BANDS.index(b) for b in LAND]
FULLSPEC_SCHEMA = "siac_multisensor_fullspec_prepared_npz_v1"
CONTRACT = "siac_multisensor_same_date_fullspec_trees_visible_6s_extension_v1"


def dense_label(teacher_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(H, W, 11) values and sigma in LAND order, plus the label's own validity."""
    with np.load(teacher_path, allow_pickle=False) as z:
        names = [str(b) for b in z["surface_bands"]]
        surface, sigma = z["surface"], z["surface_uncertainty"]
        anchors = {str(b): z["anchor_boa_at_solution"][..., k] for k, b in enumerate(z["anchor_boa_bands"])}
        extra = {str(b): z["extra_anchor_boa_at_solution"][..., k] for k, b in enumerate(z["extra_anchor_bands"])}
        valid = np.asarray(z["valid"]).astype(bool)
    value = np.stack([surface[..., names.index(b)] if b in TREE_BANDS else {**anchors, **extra}[b] for b in LAND], -1)
    spread = np.stack([sigma[..., names.index(b)] for b in LAND], -1)
    return value.astype(np.float32), spread.astype(np.float32), valid


def build_one(record: dict, out_dir: Path, *, teacher_root: Path, queries: int, offsets: dict) -> dict:
    source = Path(record["prepared_path"])
    original = Path(record["inputs"]["original_prepared"]["path"])
    teacher_path = teacher_root / f"{record['matchup_id']}.npz"
    if not teacher_path.is_file():
        return {"status": "no_supervision", "reason": "no_fullspec_teacher"}
    with np.load(source, allow_pickle=False) as archive:
        data = {key: archive[key] for key in archive.files}
    with np.load(original, allow_pickle=False) as archive:
        weight_map = np.asarray(archive["weight"], dtype=np.float32)
    dense, dense_sigma, label_valid_map = dense_label(teacher_path)
    height, width = dense.shape[:2]
    if data["local_toa"].shape[:2] != (height, width) or weight_map.shape != (height, width):
        raise ValueError("The full-spectrum teacher and the joint sample use different grids")
    if data["target_band_names"].tolist() != list(BANDS):
        raise ValueError("Unexpected target band layout")
    library_queries = int(np.asarray(data["query_indices"]).size)

    eligible = (
        (weight_map > 0) & label_valid_map
        & np.isfinite(dense).all(-1) & np.isfinite(dense_sigma).all(-1) & (dense_sigma > 0).all(-1)
        & data["local_valid"][..., TARGET].all(-1)
    )
    flat = np.setdiff1d(np.flatnonzero(eligible.ravel()), np.asarray(data["query_indices"]).ravel())
    if flat.size < 256:
        return {"status": "no_supervision", "reason": "fewer_than_256_dense_pixels", "available": int(flat.size)}
    seed = int(hashlib.sha256((record["matchup_id"] + "|fullspec").encode()).hexdigest()[:16], 16)
    pick = np.sort(np.random.default_rng(seed).choice(flat, min(queries, flat.size), replace=False))

    visible_offset, visible_spread = visible_correction(offsets, str(record["sensor"]))
    offset = np.zeros(len(LAND), np.float32)
    spread = np.zeros(len(LAND), np.float32)
    for k, band in enumerate(("B02", "B03", "B04")):
        offset[LAND.index(band)], spread[LAND.index(band)] = visible_offset[k], visible_spread[k]
    values = dense.reshape(-1, len(LAND))[pick] + offset[None]
    sigma = np.sqrt(dense_sigma.reshape(-1, len(LAND))[pick] ** 2 + spread[None] ** 2)
    weight = np.clip(weight_map.reshape(-1)[pick], 0.0, 1.0)

    count, bands = pick.size, len(BANDS)
    teacher = np.full((count, bands), np.nan, np.float32)
    uncertainty = np.full((count, bands), np.nan, np.float32)
    quality = np.zeros((count, bands), np.float32)
    label_valid = np.zeros((count, bands), bool)
    good = np.isfinite(values).all(-1) & np.isfinite(sigma).all(-1) & (sigma > 0).all(-1) & (weight > 0)
    teacher[np.ix_(good, TARGET)] = values[good]
    uncertainty[np.ix_(good, TARGET)] = sigma[good]
    quality[np.ix_(good, TARGET)] = weight[good, None]
    label_valid[np.ix_(good, TARGET)] = True

    library_valid = np.asarray(data["target_label_valid"]).copy()
    library_dropped = int(library_valid[:, TARGET].sum())
    library_valid[:, TARGET] = False

    data["teacher"] = np.concatenate([data["teacher"], teacher], axis=0)
    data["uncertainty"] = np.concatenate([data["uncertainty"], uncertainty], axis=0)
    data["quality_weight"] = np.concatenate([data["quality_weight"], quality], axis=0)
    data["target_label_valid"] = np.concatenate([library_valid, label_valid], axis=0)
    data["native_quality_weight"] = np.concatenate(
        [data["native_quality_weight"], np.zeros((count, data["local_toa"].shape[-1]), np.float32)], axis=0)
    data["query_indices"] = np.concatenate(
        [np.asarray(data["query_indices"]).ravel(), pick.astype(data["query_indices"].dtype)])
    data["query_valid"] = np.concatenate([np.asarray(data["query_valid"]).ravel(), np.ones(count, bool)])
    handled = {"teacher", "uncertainty", "quality_weight", "target_label_valid",
               "native_quality_weight", "query_indices", "query_valid"}
    for key, value in list(data.items()):
        if key in handled or value.ndim == 0:
            continue
        if value.shape[0] == library_queries and key not in ("target_band_names", "input_band_names"):
            data[key] = _extend(value, count)
    data["schema_version"] = np.asarray(SCHEMA)
    data["fullspec_schema_version"] = np.asarray(FULLSPEC_SCHEMA)
    data["dense_label_bands"] = np.asarray(LAND)
    data["dense_label_source"] = np.asarray("same-date full-spectrum M5 teacher: trees B01-B04, 6S at M5 AOD B05-B12")
    data["dense_query_count"] = np.asarray(count, np.int32)
    data["library_query_count"] = np.asarray(library_queries, np.int32)
    data["platform_offset"] = offset
    data["platform_offset_spread"] = spread

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{record['matchup_id']}.npz"
    np.savez_compressed(path, **data)
    return {"status": "prepared", "prepared_path": str(path), "prepared_sha256": digest(path),
            "query_count": int(data["query_indices"].size), "dense_queries": count,
            "dense_labelled": int(good.sum()), "library_labels_dropped": library_dropped,
            "extension_labels": int(library_valid.sum()),
            "median_sigma": float(np.nanmedian(sigma)), "platform_offset": offset.tolist()}


def run(args: argparse.Namespace) -> int:
    release = json.loads(Path(args.release).read_text())
    offsets = json.loads(Path(args.platform_offsets).read_text())
    records = [r for r in release["records"] if "original_prepared" in r.get("inputs", {}) and r.get("prepared_path")]
    mine = records[args.shard :: args.shards]
    out_dir, status_dir = Path(args.out) / "prepared", Path(args.out) / "status"
    status_dir.mkdir(parents=True, exist_ok=True)
    failures, written = 0, []
    for record in mine:
        entry = {"matchup_id": record["matchup_id"], "sensor": record.get("sensor")}
        try:
            path = out_dir / f"{record['matchup_id']}.npz"
            if path.is_file() and not args.force:
                entry.update(status="cached", prepared_path=str(path), prepared_sha256=digest(path))
            else:
                entry.update(build_one(record, out_dir, teacher_root=Path(args.teacher_root),
                                       queries=args.queries, offsets=offsets))
        except Exception as exc:  # noqa: BLE001 - one scene must not sink the shard
            failures += 1
            traceback.print_exc()
            entry.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        written.append(entry)
        print(json.dumps(entry, sort_keys=True), flush=True)
    (status_dir / f"shard_{args.shard:04d}.json").write_text(json.dumps(written, indent=1))
    return 1 if failures else 0


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--release", required=True)
    value.add_argument("--teacher-root", required=True)
    value.add_argument("--out", required=True)
    value.add_argument("--platform-offsets", required=True)
    value.add_argument("--queries", type=int, default=4096)
    value.add_argument("--shards", type=int, default=1)
    value.add_argument("--shard", type=int, default=0)
    value.add_argument("--force", action="store_true")
    return value


def main() -> None:
    sys.exit(run(parser().parse_args()))


if __name__ == "__main__":
    main()
