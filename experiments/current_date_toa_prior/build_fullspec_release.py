"""Assemble the full-spectrum prepared release from the joint campaign's release.

S2 records point at the archives rebuilt by ``build_fullspec_multisensor_labels`` --
the acquisition's own M5 teacher on every land band (ExtraTrees B01-B04, 6S at the
M5-solved AOD B05-B12) -- and declare that provenance. Landsat records have no
same-date teacher and are carried through untouched, archive and digest included.
S2 scenes without a usable full-spectrum label keep their original archive.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

LAND = ("B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12")
CONTRACT = "siac_multisensor_same_date_fullspec_trees_visible_6s_extension_v1"


def fullspec_provenance(original: dict) -> dict:
    value = dict(original)
    value.update(
        month_matched=False,
        same_date=True,
        state_qa_applied=True,
        dense_label_bands=list(LAND),
        visible_label_source="same-date M5 teacher ExtraTrees surface (B01-B04)",
        extension_label_source="6S correction of the scene TOA at the M5-solved AOD (B05-B08, B8A, B11, B12)",
        platform_correction="measured constant to canonical S2A on B02-B04; zero elsewhere (unmeasured)",
        simultaneous_surface_truth=False,
        scientific_accuracy_validated=False,
        AERONET_used=False,
        supersedes="month-matched library labels on every land band",
    )
    return value


def run(args: argparse.Namespace) -> int:
    release = json.loads(Path(args.release).read_text())
    built: dict[str, dict] = {}
    for shard in sorted(Path(args.status).glob("shard_*.json")):
        for entry in json.loads(shard.read_text()):
            if entry.get("status") in ("prepared", "cached"):
                built[entry["matchup_id"]] = entry
    records = []
    rebuilt = carried = dropped = 0
    for record in release["records"]:
        if "original_prepared" not in record.get("inputs", {}):
            records.append(record)
            carried += 1
            continue
        entry = built.get(record["matchup_id"])
        if entry is None:
            # A handful of scenes have almost no usable same-date visible pixels
            # -- heavy cloud, mostly zero. Carrying them on their original
            # archive keeps the corpus identical to the campaign being compared
            # against, so the labels stay the only difference; their visible
            # supervision is simply the month-matched one it always was.
            records.append(record)
            dropped += 1
            continue
        value = dict(record)
        value["prepared_path"] = entry["prepared_path"]
        value["prepared_sha256"] = entry["prepared_sha256"]
        value["teacher_contract"] = CONTRACT
        value["teacher_provenance"] = fullspec_provenance(record.get("teacher_provenance", {}))
        value["query_count"] = entry.get("query_count", value.get("query_count"))
        value["dense_queries"] = entry.get("dense_queries")
        value["labelled_band_values"] = entry.get("extension_labels", 0) + len(LAND) * entry.get(
            "dense_labelled", 0
        )
        value["source_prepared_path"] = record["prepared_path"]
        records.append(value)
        rebuilt += 1
    if dropped:
        print(
            f"{dropped} S2 records kept their original archive (no usable full-spectrum label)",
            file=sys.stderr,
        )
    out = dict(release)
    out["records"] = records
    out["fullspec_label_build"] = {
        "rebuilt_s2_records": rebuilt,
        "carried_records": carried,
        "s2_records_kept_on_original_archive": dropped,
        "platform_offsets": json.loads(Path(args.platform_offsets).read_text()),
        "source_release": str(Path(args.release).absolute()),
    }
    for split in ("train", "development"):
        if not any(r["split"] == split for r in records):
            raise ValueError(f"No {split} records survived the rebuild")
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1))
    counts: dict[str, int] = {}
    for record in records:
        counts[record["split"]] = counts.get(record["split"], 0) + 1
    print(f"rebuilt {rebuilt} S2 records, carried {carried} Landsat records")
    print(f"splits: {counts}")
    print(f"wrote {path}")
    return 0


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--release", required=True)
    value.add_argument("--status", required=True)
    value.add_argument("--platform-offsets", required=True)
    value.add_argument("--out", required=True)
    return value


def main() -> None:
    sys.exit(run(parser().parse_args()))


if __name__ == "__main__":
    main()
