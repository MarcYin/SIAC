#!/usr/bin/env python3
"""Add per-image sun and view geometry to a locally composed winner index.

``build_cloudscore_index_local`` writes ``day``, ``idx`` and ``ratio`` per
image, and declares the ``siac_l1c_cloudscore_winner_index_v1`` schema. Its
consumer -- ``live_l1c_library.read_mosaic_index`` -- also requires ``sza``,
``saa``, ``vza`` and ``vaa``, because the historical composites have to be
corrected through 6S at each contributing acquisition's own geometry. Declaring
a schema is not the same as satisfying it, and nothing checked.

The angles come from Earth Engine scene metadata rather than from per-image tile
XML: one listing per scene returns all of them, against roughly 111 XML reads
per scene the other way. View angles are B2's, matching the reference index.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Sequence

COLLECTION = "COPERNICUS/S2_HARMONIZED"
SERVER_URL = "https://earthengine-highvolume.googleapis.com"

#: Earth Engine property names, and the index field each becomes.
ANGLE_PROPERTIES = {
    "sza": "MEAN_SOLAR_ZENITH_ANGLE",
    "saa": "MEAN_SOLAR_AZIMUTH_ANGLE",
    "vza": "MEAN_INCIDENCE_ZENITH_ANGLE_B2",
    "vaa": "MEAN_INCIDENCE_AZIMUTH_ANGLE_B2",
}

REQUIRED_FIELDS = ("day", "idx", "sza", "saa", "vza", "vaa")


def needs_geometry(path: Path) -> bool:
    """True when an index is missing any angle the library reader requires."""

    archive = np.load(path, allow_pickle=False)
    if "image_table" not in archive.files:
        return True
    rows = [json.loads(str(row)) for row in archive["image_table"]]
    if not rows:
        return False
    return any(field not in rows[0] for field in REQUIRED_FIELDS)


def fetch_angles(image_ids: Sequence[str]) -> dict[str, dict[str, float]]:
    """``system:index -> angles`` for the given acquisitions, in one request."""

    import ee
    from edown.auth import initialize_earth_engine

    initialize_earth_engine(SERVER_URL)
    wanted = sorted({str(value).rsplit("/", 1)[-1] for value in image_ids if value})
    if not wanted:
        return {}
    collection = ee.ImageCollection(COLLECTION).filter(
        ee.Filter.inList("system:index", list(wanted))
    )
    properties = ["system:index", *ANGLE_PROPERTIES.values()]
    records = collection.reduceColumns(ee.Reducer.toList(len(properties)), properties).getInfo()
    resolved: dict[str, dict[str, float]] = {}
    for row in records.get("list", []) or []:
        index = str(row[0])
        values = {}
        for position, field in enumerate(ANGLE_PROPERTIES, start=1):
            value = row[position]
            if value is None:
                break
            values[field] = float(value)
        if len(values) == len(ANGLE_PROPERTIES):
            resolved[index] = values
    return resolved


def augment(path: Path, *, dry_run: bool = False) -> dict[str, Any]:
    """Rewrite one index with angles attached to every image row."""

    archive = np.load(path, allow_pickle=False)
    rows = [json.loads(str(row)) for row in archive["image_table"]]
    angles = fetch_angles([row["idx"] for row in rows])
    missing = [row["idx"] for row in rows if str(row["idx"]) not in angles]
    if missing:
        # A row without geometry cannot be corrected, and silently dropping it
        # would change the composite; fail loudly instead.
        raise RuntimeError(
            f"{path.name}: no Earth Engine geometry for {len(missing)} of "
            f"{len(rows)} images, first {missing[0]}"
        )
    updated = [{**row, **angles[str(row["idx"])]} for row in rows]
    if dry_run:
        return {"images": len(rows), "written": False}

    payload = {key: archive[key] for key in archive.files}
    payload["image_table"] = np.asarray([json.dumps(r, sort_keys=True) for r in updated])
    # The temp name must itself end in .npz: savez_compressed appends the
    # extension when it is absent, so a ".npz.tmp" target is written as
    # ".npz.tmp.npz" and the rename then finds nothing.
    temporary = path.with_name(f"{path.stem}.tmp.npz")
    np.savez_compressed(temporary, **payload)
    temporary.replace(path)
    return {"images": len(rows), "written": True}


def run(args: argparse.Namespace) -> dict[str, Any]:
    root = Path(args.index_root)
    paths = sorted(root.glob("*.npz"))
    if args.shard_count > 1:
        paths = [p for i, p in enumerate(paths) if i % args.shard_count == args.shard]
    if args.limit:
        paths = paths[: int(args.limit)]

    updated = skipped = 0
    failures: dict[str, str] = {}
    started = time.perf_counter()
    for position, path in enumerate(paths):
        try:
            if not args.force and not needs_geometry(path):
                skipped += 1
                continue
            augment(path, dry_run=args.dry_run)
            updated += 1
        except Exception as error:  # noqa: BLE001 - one bad index must not end the shard
            failures[path.stem] = f"{type(error).__name__}: {error}"
        if (position + 1) % 25 == 0:
            print(
                f"  {position + 1}/{len(paths)} indices, {updated} updated, "
                f"{skipped} already complete, {len(failures)} failed",
                flush=True,
            )

    summary = {
        "indices": len(paths),
        "updated": updated,
        "already_complete": skipped,
        "failed": len(failures),
        "failures": dict(list(failures.items())[:10]),
        "wall_seconds": round(time.perf_counter() - started, 1),
    }
    if args.summary:
        Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
        Path(args.summary).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--index-root", required=True)
    value.add_argument("--force", action="store_true")
    value.add_argument("--dry-run", action="store_true")
    value.add_argument("--limit", type=int, default=0)
    value.add_argument("--shard", type=int, default=0)
    value.add_argument("--shard-count", type=int, default=1)
    value.add_argument("--summary")
    return value


def main() -> None:
    print(json.dumps(run(parser().parse_args()), indent=2))


if __name__ == "__main__":
    main()
