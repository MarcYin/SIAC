"""Build a deterministic bad-case plus AOD-stratified control subset."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
AOD_BINS = (
    (-math.inf, 0.1, "<0.1"),
    (0.1, 0.2, "0.1-0.2"),
    (0.2, 0.4, "0.2-0.4"),
    (0.4, 0.6, "0.4-0.6"),
    (0.6, 1.0, "0.6-1.0"),
    (1.0, 1.5, "1.0-1.5"),
    (1.5, math.inf, ">=1.5"),
)


def _ids(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def _rows(directory: Path) -> dict[str, dict[str, Any]]:
    rows = {}
    for path in directory.glob("*.json"):
        try:
            row = json.loads(path.read_text())
        except (OSError, ValueError):
            continue
        matchup_id = str(row.get("matchup_id") or "")
        if matchup_id:
            rows[matchup_id] = row
    return rows


def _truth(row: dict[str, Any]) -> float:
    return float(row.get("truth", row.get("aeronet_aod550_mean")))


def _hit(row: dict[str, Any]) -> bool:
    if row.get("within_ee") is not None:
        return bool(row["within_ee"])
    if row.get("flag") is not None:
        return str(row["flag"]).upper() == "OK"
    return abs(float(row["err"])) <= float(row["ee"])


def _aod_bin(value: float) -> str:
    return next(label for lower, upper, label in AOD_BINS if lower <= value < upper)


def _stable_order(matchup_id: str, seed: str) -> str:
    return hashlib.sha256(f"{seed}:{matchup_id}".encode()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--controls-per-bin", type=int, default=8)
    parser.add_argument("--seed", default="surface-e1-20260709")
    parser.add_argument("--output", type=Path, default=ROOT / "surface_e1_matched_mids.txt")
    parser.add_argument("--manifest", type=Path, default=ROOT / "surface_e1_matched_manifest.json")
    args = parser.parse_args()
    if args.controls_per_bin < 1:
        raise ValueError("controls-per-bin must be positive")

    campaign = _ids(ROOT / "campaign250_mids.txt")
    bad = set(_ids(ROOT / "campaign250_surface_bad89_mids.txt"))
    r2 = _rows(ROOT / "phaseD_results_campaign250_R2_full")
    controls_by_bin: dict[str, list[str]] = {label: [] for _, _, label in AOD_BINS}
    for matchup_id in campaign:
        row = r2.get(matchup_id)
        if row is None or matchup_id in bad or not _hit(row):
            continue
        controls_by_bin[_aod_bin(_truth(row))].append(matchup_id)

    controls = set()
    for candidates in controls_by_bin.values():
        candidates.sort(key=lambda matchup_id: _stable_order(matchup_id, args.seed))
        controls.update(candidates[: args.controls_per_bin])
    selected = [
        matchup_id for matchup_id in campaign if matchup_id in bad or matchup_id in controls
    ]
    manifest_rows = []
    for matchup_id in selected:
        row = r2.get(matchup_id, {})
        truth = _truth(row) if row else None
        manifest_rows.append(
            {
                "matchup_id": matchup_id,
                "role": "bad" if matchup_id in bad else "control_hit",
                "truth": truth,
                "aod_bin": _aod_bin(truth) if truth is not None else None,
                "existing_species_prior": (
                    ROOT / "l1c_seasonal_species_lut" / f"{matchup_id}.npz"
                ).exists(),
            }
        )
    args.output.write_text("\n".join(selected) + "\n", encoding="utf-8")
    args.manifest.write_text(
        json.dumps(
            {
                "seed": args.seed,
                "controls_per_bin": args.controls_per_bin,
                "count": len(selected),
                "bad_count": sum(row["role"] == "bad" for row in manifest_rows),
                "control_count": sum(row["role"] == "control_hit" for row in manifest_rows),
                "rows": manifest_rows,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(
        f"wrote {args.output}: {len(selected)} records "
        f"({len(bad & set(selected))} bad, {len(controls)} controls)"
    )


if __name__ == "__main__":
    main()
