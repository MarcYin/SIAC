"""Build a fixed campaign cohort from retrieval cloud fractions."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

DEFAULT_ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_CLOUD_DIR = "phaseD_results_campaign250_masked_r2_l2awvp_6s_20260710"
DEFAULT_MIDS_FILE = "campaign250_mids.txt"
DEFAULT_OUTPUT_FILE = "campaign250_lowcloud20_mids.txt"


def _finite_float(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def build_cohort(
    campaign_ids: list[str],
    records: dict[str, dict[str, Any]],
    *,
    threshold: float,
    max_unclassified: int = 0,
) -> list[str]:
    """Return campaign IDs whose recorded cloud fraction is below threshold."""

    if not 0.0 <= threshold <= 1.0:
        raise ValueError("threshold must be between 0 and 1")
    if len(campaign_ids) != len(set(campaign_ids)):
        raise ValueError("campaign matchup IDs must be unique")

    missing = [matchup_id for matchup_id in campaign_ids if matchup_id not in records]
    if missing:
        raise ValueError(f"missing cloud records for {len(missing)} campaign matchups")

    invalid = [
        matchup_id
        for matchup_id in campaign_ids
        if _finite_float(records[matchup_id].get("cloud_frac")) is None
    ]
    if len(invalid) > max_unclassified:
        raise ValueError(f"invalid cloud_frac for {len(invalid)} campaign matchups")

    return [
        matchup_id
        for matchup_id in campaign_ids
        if (cloud_fraction := _finite_float(records[matchup_id].get("cloud_frac"))) is not None
        and cloud_fraction < threshold
    ]


def _load_records(directory: Path) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for path in sorted(directory.glob("*.json")):
        record = json.loads(path.read_text(encoding="utf-8"))
        matchup_id = str(record.get("matchup_id") or path.stem)
        if matchup_id in records:
            raise ValueError(f"duplicate matchup_id: {matchup_id}")
        records[matchup_id] = record
    return records


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--cloud-dir", default=DEFAULT_CLOUD_DIR)
    parser.add_argument("--campaign-mids", default=DEFAULT_MIDS_FILE)
    parser.add_argument("--output", default=DEFAULT_OUTPUT_FILE)
    parser.add_argument("--threshold", type=float, default=0.20)
    parser.add_argument("--expect-campaign", type=int, default=250)
    parser.add_argument("--expect-selected", type=int, default=152)
    parser.add_argument("--max-unclassified", type=int, default=1)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    root = args.root
    campaign_path = root / args.campaign_mids
    cloud_dir = root / args.cloud_dir
    output_path = root / args.output
    campaign_ids = [
        line.strip()
        for line in campaign_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(campaign_ids) != args.expect_campaign:
        raise SystemExit(f"expected {args.expect_campaign} campaign IDs, found {len(campaign_ids)}")

    records = _load_records(cloud_dir)
    cohort = build_cohort(
        campaign_ids,
        records,
        threshold=args.threshold,
        max_unclassified=args.max_unclassified,
    )
    if len(cohort) != args.expect_selected:
        raise SystemExit(f"expected {args.expect_selected} selected IDs, found {len(cohort)}")

    serialized = "\n".join(cohort) + "\n"
    output_path.write_text(serialized, encoding="utf-8")
    digest = hashlib.sha256(serialized.encode("utf-8")).hexdigest()
    unclassified = [
        matchup_id
        for matchup_id in campaign_ids
        if _finite_float(records[matchup_id].get("cloud_frac")) is None
    ]
    manifest = {
        "campaign_mids": str(campaign_path),
        "cloud_source_directory": str(cloud_dir),
        "cloud_field": "cloud_frac",
        "comparison": "strictly_less_than",
        "threshold_fraction": args.threshold,
        "campaign_count": len(campaign_ids),
        "selected_count": len(cohort),
        "unclassified_count": len(unclassified),
        "unclassified_matchup_ids": unclassified,
        "selected_mids_sha256": digest,
        "selected_mids_file": str(output_path),
    }
    manifest_path = output_path.with_suffix(".json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(
        f"selected {len(cohort)}/{len(campaign_ids)} matchups with cloud_frac "
        f"< {args.threshold:.3f}; sha256={digest}"
    )
    print(output_path)
    print(manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
