"""Assemble the hybrid-label prepared release from the joint campaign's release.

S2 records point at the archives rebuilt by ``build_hybrid_multisensor_labels``
-- same-date M5 teacher on B02/B03/B04, month-matched library labels on the
remaining land bands -- and declare both provenances explicitly. Landsat records
have no same-date teacher and are carried through untouched, archive and digest
included.

Everything the trainer audits is preserved: schema, the geographic audit
attestation, the global target definition, and the split/geography separation.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

VISIBLE = ("B02", "B03", "B04")
CONTRACT = "siac_multisensor_same_date_visible_plus_seasonal_extension_v1"


def hybrid_provenance(original: dict) -> dict:
    """Provenance for a record carrying two label families, naming each one."""
    value = dict(original)
    value.update(
        # True of the extension bands, whose labels are unchanged.
        month_matched=True,
        # True of the visible bands, which now carry the acquisition's own retrieval.
        same_date=True,
        state_qa_applied=True,
        visible_bands=list(VISIBLE),
        visible_label_source="same-date M5 teacher (original_prepared) at a dense query block",
        visible_label_platform_correction="measured constant to canonical S2A; spread added to sigma",
        extension_label_source="joint current native-SRF TOA fit over exact-month historical realizations",
        # The M5 teacher is a retrieval, not a measurement; neither family is truth.
        simultaneous_surface_truth=False,
        scientific_accuracy_validated=False,
        AERONET_used=False,
        supersedes="month-matched library labels on the visible bands",
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
            dropped += 1
            continue
        value = dict(record)
        value["prepared_path"] = entry["prepared_path"]
        value["prepared_sha256"] = entry["prepared_sha256"]
        value["teacher_contract"] = CONTRACT
        value["teacher_provenance"] = hybrid_provenance(record.get("teacher_provenance", {}))
        value["query_count"] = entry.get("query_count", value.get("query_count"))
        value["visible_queries"] = entry.get("visible_queries")
        value["labelled_band_values"] = entry.get("extension_labels", 0) + 3 * entry.get(
            "visible_labelled", 0
        )
        value["source_prepared_path"] = record["prepared_path"]
        records.append(value)
        rebuilt += 1
    if dropped:
        print(f"{dropped} S2 records had no hybrid archive and were dropped", file=sys.stderr)
    out = dict(release)
    out["records"] = records
    out["hybrid_label_build"] = {
        "rebuilt_s2_records": rebuilt,
        "carried_records": carried,
        "dropped_records": dropped,
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
