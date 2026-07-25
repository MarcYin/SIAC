"""Select deterministic non-campaign matchups for source-trust training.

The campaign-250 set intentionally oversamples high AOD. To fit source-selection
rules without using those 250 labels, this script samples matchups from the
larger AERONET pool while excluding the exact campaign matchup IDs. High-AOD
bins are kept dense because they are rare.
"""

from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
MATCHUPS = ROOT / "matchups" / "matchups.csv"
CAMPAIGN_MIDS = ROOT / "campaign250_mids.txt"
DEFAULT_OUT = ROOT / "external_train1000_mids.txt"

BIN_PLAN = (
    ("lt010", 0.0, 0.1, 143),
    ("010_020", 0.1, 0.2, 150),
    ("020_040", 0.2, 0.4, 200),
    ("040_060", 0.4, 0.6, 200),
    ("060_100", 0.6, 1.0, 191),
    ("100_150", 1.0, 1.5, 88),
    ("gt150", 1.5, float("inf"), 28),
)


def _read_campaign_mids(path: Path) -> set[str]:
    return {line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()}


def select_mids(seed: int, target: int) -> list[str]:
    campaign = _read_campaign_mids(CAMPAIGN_MIDS)
    rows = [row for row in csv.DictReader(MATCHUPS.open()) if row["matchup_id"] not in campaign]
    rng = random.Random(seed)
    selected: list[str] = []
    selected_set: set[str] = set()
    for _label, lo, hi, quota in BIN_PLAN:
        pool = [row for row in rows if lo <= float(row["aeronet_aod550_mean"]) < hi]
        rng.shuffle(pool)
        take = min(quota, len(pool), max(0, target - len(selected)))
        for row in pool[:take]:
            selected.append(row["matchup_id"])
            selected_set.add(row["matchup_id"])
    if len(selected) < target:
        remaining = [row for row in rows if row["matchup_id"] not in selected_set]
        rng.shuffle(remaining)
        for row in remaining[: target - len(selected)]:
            selected.append(row["matchup_id"])
    return selected[:target]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20260704)
    parser.add_argument("--target", type=int, default=1000)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    mids = select_mids(args.seed, args.target)
    args.out.write_text("\n".join(mids) + "\n", encoding="utf-8")
    print(f"Wrote {len(mids)} matchup IDs to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
