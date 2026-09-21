"""Restore same-date visible supervision to the joint Landsat/S2 training samples.

The joint multisensor campaign supervises every band from a month-matched
library: for each scene ``fit_teacher`` refits the current TOA against
exact-month composites drawn from OTHER years, jointly solving AOD. That label
is 3.8x noisier per pixel than the same-date M5 teacher (sigma 0.0254 vs
0.0067) and covers only ~560 of 1024 fixed query pixels rather than the
teacher's 140k, so each scene carries roughly 60x less information than the
S2-only models received. Its own provenance declares
``simultaneous_surface_truth: False``.

The same-date teacher exists for every S2 scene in the corpus -- the release
records carry it as ``inputs.original_prepared`` -- and was deliberately not
copied ("The old RGB teacher is never copied"); ``fit_teacher`` reads its
``weight`` only as an eligibility screen. This rebuilds each S2 sample so that:

* B02/B03/B04 are supervised by the dense same-date M5 teacher, sampled at a
  large fresh query block, and
* the eight remaining land bands keep their month-matched library labels on the
  original query block, whose visible entries are masked off.

Landsat samples have no same-date teacher and are left untouched.

The joint model predicts canonical S2A reflectance while the teacher sits in
its own spacecraft's band space. That correction is small in the visible, but
it is not assumed away: for S2B/S2C it is a constant measured across a large
scene sample by ``estimate_platform_offsets``, applied with the scene-to-scene
spread folded into the label sigma. S2A scenes take the identity. See
:func:`visible_correction` for why the offset is global rather than per pixel.

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

BANDS = ("B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B09", "B11", "B12")
LAND_BANDS = tuple(b for b in BANDS if b != "B09")
VISIBLE = ("B02", "B03", "B04")
#: Target-band and input-band positions of the visible triple (both use BANDS).
VISIBLE_TARGET = [BANDS.index(b) for b in VISIBLE]
#: The same-date teacher stores exactly B02, B03, B04 on its last axis, in that
#: order. Verified against the library labels: pairing channel k with VISIBLE[k]
#: is ~3x closer than any permutation.
DENSE_BANDS = VISIBLE
#: The archive keeps the audited NPZ contract -- its arrays and dtypes are
#: unchanged, only the query block is longer -- so the released trainer reads it
#: without modification. The hybrid provenance rides in its own keys below.
SCHEMA = "siac_multisensor_prepared_npz_v1"
HYBRID_SCHEMA = "siac_multisensor_hybrid_prepared_npz_v1"
CONTRACT = "siac_multisensor_same_date_visible_plus_seasonal_extension_v1"


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def visible_correction(offsets: dict, sensor: str) -> tuple[np.ndarray, np.ndarray]:
    """Constant B02/B03/B04 offset from ``sensor`` to canonical S2A, and its spread.

    A per-pixel offset was tried first and rejected: bridging a single scene's
    own spectra puts the correction wherever the k-nearest-neighbour dictionary
    extrapolates, which on a dark scene produced -0.0099 in B02 against a
    population value of +0.0018 -- a 40% correction on a surface of reflectance
    0.024, and an artifact rather than physics. The offset is instead the median
    over a large scene sample (``estimate_platform_offsets``), applied as a
    constant with the scene-to-scene spread folded into the label sigma. It is
    small either way: every component is well under the teacher's own sigma.
    """
    if sensor == "S2A":
        return np.zeros(len(VISIBLE), np.float32), np.zeros(len(VISIBLE), np.float32)
    entry = offsets.get(sensor)
    if entry is None:
        raise ValueError(f"No measured {sensor} -> S2A visible offset is available")
    offset = np.asarray(entry["per_scene_median"], dtype=np.float32)
    spread = np.asarray(entry["per_scene_sd"], dtype=np.float32)
    if offset.shape != (len(VISIBLE),) or spread.shape != (len(VISIBLE),):
        raise ValueError("Platform offset must give one value and one spread per visible band")
    if not np.isfinite(offset).all() or not np.isfinite(spread).all() or (spread < 0).any():
        raise ValueError("Platform offset and spread must be finite and nonnegative")
    return offset, spread


def _extend(value: np.ndarray, count: int) -> np.ndarray:
    """Append ``count`` unsupervised rows to a per-query audit array."""
    tail = list(value.shape)
    tail[0] = count
    if value.dtype.kind == "f":
        block = np.full(tail, np.nan, dtype=value.dtype)
    elif value.dtype.kind == "b":
        block = np.zeros(tail, dtype=value.dtype)
    else:
        block = np.full(tail, -1, dtype=value.dtype)
    return np.concatenate([value, block], axis=0)


def build_one(record: dict, out_dir: Path, *, visible_queries: int, offsets: dict) -> dict:
    source = Path(record["prepared_path"])
    original = Path(record["inputs"]["original_prepared"]["path"])
    sensor = str(record["sensor"])
    with np.load(source, allow_pickle=False) as archive:
        data = {key: archive[key] for key in archive.files}
    with np.load(original, allow_pickle=False) as archive:
        dense = np.asarray(archive["teacher"], dtype=np.float32)
        dense_sigma = np.asarray(archive["uncertainty"], dtype=np.float32)
        dense_weight = np.asarray(archive["weight"], dtype=np.float32)
    height, width = dense.shape[:2]
    if data["local_toa"].shape[:2] != (height, width):
        raise ValueError("The same-date teacher and the joint sample use different grids")
    if data["target_band_names"].tolist() != list(BANDS):
        raise ValueError("Unexpected target band layout")
    # fit_teacher caps the query block at 1024 but takes fewer where the scene
    # offers fewer eligible pixels, so the existing width is read, never assumed.
    library_queries = int(np.asarray(data["query_indices"]).size)

    eligible = (
        (dense_weight > 0)
        & np.isfinite(dense).all(-1)
        & np.isfinite(dense_sigma).all(-1)
        & (dense_sigma > 0).all(-1)
        & data["local_valid"][..., VISIBLE_TARGET].all(-1)
    )
    flat = np.flatnonzero(eligible.ravel())
    flat = np.setdiff1d(flat, np.asarray(data["query_indices"]).ravel())
    if flat.size < 256:
        return {"status": "no_supervision", "reason": "fewer_than_256_dense_visible_pixels",
                "available": int(flat.size)}
    seed = int(hashlib.sha256((record["matchup_id"] + "|hybrid").encode()).hexdigest()[:16], 16)
    rng = np.random.default_rng(seed)
    pick = np.sort(rng.choice(flat, min(visible_queries, flat.size), replace=False))

    offset, spread = visible_correction(offsets, sensor)
    values = dense.reshape(height * width, len(DENSE_BANDS))[pick] + offset[None]
    sigma = np.sqrt(
        dense_sigma.reshape(height * width, len(DENSE_BANDS))[pick] ** 2 + spread[None] ** 2
    )
    weight = np.clip(dense_weight.reshape(height * width)[pick], 0.0, 1.0)

    count = pick.size
    bands = len(BANDS)
    teacher = np.full((count, bands), np.nan, dtype=np.float32)
    uncertainty = np.full((count, bands), np.nan, dtype=np.float32)
    quality = np.zeros((count, bands), dtype=np.float32)
    label_valid = np.zeros((count, bands), dtype=bool)
    good = np.isfinite(values).all(-1) & np.isfinite(sigma).all(-1) & (sigma > 0).all(-1) & (weight > 0)
    teacher[np.ix_(good, VISIBLE_TARGET)] = values[good]
    uncertainty[np.ix_(good, VISIBLE_TARGET)] = sigma[good]
    quality[np.ix_(good, VISIBLE_TARGET)] = weight[good, None]
    label_valid[np.ix_(good, VISIBLE_TARGET)] = True

    # The library's own visible labels are the noisy ones this replaces.
    library_valid = np.asarray(data["target_label_valid"]).copy()
    library_visible = int(library_valid[:, VISIBLE_TARGET].sum())
    library_valid[:, VISIBLE_TARGET] = False

    data["teacher"] = np.concatenate([data["teacher"], teacher], axis=0)
    data["uncertainty"] = np.concatenate([data["uncertainty"], uncertainty], axis=0)
    data["quality_weight"] = np.concatenate([data["quality_weight"], quality], axis=0)
    data["target_label_valid"] = np.concatenate([library_valid, label_valid], axis=0)
    # Native RT closure stays on the library block, whose scene-uniform state
    # screens it was built under; the dense visible queries never enter it.
    data["native_quality_weight"] = np.concatenate(
        [data["native_quality_weight"], np.zeros((count, data["local_toa"].shape[-1]), np.float32)],
        axis=0,
    )
    data["query_indices"] = np.concatenate(
        [np.asarray(data["query_indices"]).ravel(), pick.astype(data["query_indices"].dtype)]
    )
    data["query_valid"] = np.concatenate(
        [np.asarray(data["query_valid"]).ravel(), np.ones(count, dtype=bool)]
    )
    handled = {"teacher", "uncertainty", "quality_weight", "target_label_valid",
               "native_quality_weight", "query_indices", "query_valid"}
    for key, value in list(data.items()):
        if key in handled or value.ndim == 0:
            continue
        if value.shape[0] == library_queries and key not in ("target_band_names", "input_band_names"):
            data[key] = _extend(value, count)

    data["schema_version"] = np.asarray(SCHEMA)
    data["hybrid_schema_version"] = np.asarray(HYBRID_SCHEMA)
    data["visible_label_source"] = np.asarray("same-date M5 teacher (original_prepared)")
    data["visible_query_count"] = np.asarray(count, dtype=np.int32)
    data["library_query_count"] = np.asarray(library_queries, dtype=np.int32)
    data["visible_platform_offset"] = offset.astype(np.float32)
    data["visible_platform_offset_spread"] = spread.astype(np.float32)

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{record['matchup_id']}.npz"
    np.savez_compressed(path, **data)
    return {
        "status": "prepared",
        "prepared_path": str(path),
        "prepared_sha256": digest(path),
        "query_count": int(data["query_indices"].size),
        "visible_queries": count,
        "visible_labelled": int(good.sum()),
        "library_visible_labels_dropped": library_visible,
        "extension_labels": int(library_valid.sum()),
        "median_visible_sigma": float(np.nanmedian(sigma)),
        "platform_offset": offset.tolist(),
    }


def run(args: argparse.Namespace) -> int:
    release = json.loads(Path(args.release).read_text())
    offsets = json.loads(Path(args.platform_offsets).read_text())
    records = [r for r in release["records"]
               if "original_prepared" in r.get("inputs", {}) and r.get("prepared_path")]
    mine = records[args.shard :: args.shards]
    out_dir = Path(args.out) / "prepared"
    status_dir = Path(args.out) / "status"
    status_dir.mkdir(parents=True, exist_ok=True)
    failures = 0
    written = []
    for record in mine:
        entry = {"matchup_id": record["matchup_id"], "sensor": record.get("sensor")}
        try:
            path = out_dir / f"{record['matchup_id']}.npz"
            if path.is_file() and not args.force:
                entry.update(status="cached", prepared_path=str(path), prepared_sha256=digest(path))
            else:
                entry.update(build_one(record, out_dir, visible_queries=args.visible_queries,
                                       offsets=offsets))
        except Exception as exc:  # noqa: BLE001 - one scene must not sink the shard
            failures += 1
            traceback.print_exc()
            entry.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        written.append(entry)
        print(json.dumps(entry, sort_keys=True), flush=True)
    (status_dir / f"shard_{args.shard:04d}.json").write_text(json.dumps(written, indent=1))
    if failures:
        print(f"{failures} of {len(mine)} scenes failed", file=sys.stderr)
        return 1
    return 0


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--release", required=True)
    value.add_argument("--out", required=True)
    value.add_argument("--platform-offsets", required=True,
                       help="JSON from estimate_platform_offsets: sensor -> per_scene_median/per_scene_sd")
    value.add_argument("--visible-queries", type=int, default=8192)
    value.add_argument("--shards", type=int, default=1)
    value.add_argument("--shard", type=int, default=0)
    value.add_argument("--force", action="store_true")
    return value


def main() -> None:
    sys.exit(run(parser().parse_args()))


if __name__ == "__main__":
    main()
