"""Evaluate scene AOD choices derived from saved per-band surface minima."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")


def band_consensus_aod(
    row: dict[str, Any],
    *,
    spread_threshold: float,
    require_flat_delta: float | None = None,
) -> tuple[float, bool]:
    """Use the median winning-family band minimum on diagnosed disagreement."""

    candidates = row.get("species_candidates") or []
    if not candidates:
        return float(row["aod"]), False
    winning = min(candidates, key=lambda candidate: float(candidate["surface_cost"]))
    spread = float(winning["band_argmin_spread"])
    curve_delta = float(winning["curve_relative_second_delta"])
    values = np.asarray(list(winning["band_argmin_aod"].values()), dtype=float)
    finite = values[np.isfinite(values)]
    use = finite.size >= 2 and np.isfinite(spread) and spread >= spread_threshold
    if require_flat_delta is not None:
        use = use and np.isfinite(curve_delta) and curve_delta <= require_flat_delta
    return (float(np.median(finite)) if use else float(row["aod"])), bool(use)


def aggregation_consensus_aod(
    row: dict[str, Any],
    *,
    spread_max: float,
    deviation_min: float,
) -> tuple[float, bool]:
    """Use agreeing band minima when final spatial aggregation contradicts them."""

    candidates = row.get("species_candidates") or []
    if not candidates:
        return float(row["aod"]), False
    winning = min(candidates, key=lambda candidate: float(candidate["surface_cost"]))
    values = np.asarray(list(winning["band_argmin_aod"].values()), dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size < 2:
        return float(row["aod"]), False
    consensus = float(np.median(finite))
    spread = float(np.max(finite) - np.min(finite))
    current = float(row["aod"])
    use = spread <= spread_max and abs(current - consensus) >= deviation_min
    return (consensus if use else current), bool(use)


def _load_tag(tag: str) -> dict[str, dict[str, Any]]:
    records = {}
    for path in (ROOT / "seasonal_val").glob(f"*__seas_maiac_{tag}.json"):
        row = json.loads(path.read_text())
        records[str(row["matchup_id"])] = row
    return records


def _hit(aod: float, truth: float) -> bool:
    return abs(aod - truth) <= 0.05 + 0.15 * truth


def _summary(
    name: str,
    rows: list[dict[str, Any]],
    *,
    spread_threshold: float | None,
    require_flat_delta: float | None = None,
    aggregation_spread_max: float | None = None,
    aggregation_deviation_min: float | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    selected = []
    for entry in rows:
        if aggregation_spread_max is not None and aggregation_deviation_min is not None:
            aod, used = aggregation_consensus_aod(
                entry["record"],
                spread_max=aggregation_spread_max,
                deviation_min=aggregation_deviation_min,
            )
        elif spread_threshold is None:
            aod, used = float(entry["record"]["aod"]), False
        else:
            aod, used = band_consensus_aod(
                entry["record"],
                spread_threshold=spread_threshold,
                require_flat_delta=require_flat_delta,
            )
        truth = float(entry["record"]["truth"])
        selected.append(
            {
                "matchup_id": entry["matchup_id"],
                "role": entry["role"],
                "truth": truth,
                "aod": aod,
                "within_ee": _hit(aod, truth),
                "used_consensus": used,
            }
        )
    errors = np.asarray([row["aod"] - row["truth"] for row in selected], dtype=float)
    by_role = {}
    for role in sorted({row["role"] for row in selected}):
        subset = [row for row in selected if row["role"] == role]
        by_role[role] = {
            "hits": sum(row["within_ee"] for row in subset),
            "count": len(subset),
        }
    return (
        {
            "strategy": name,
            "hits": sum(row["within_ee"] for row in selected),
            "count": len(selected),
            "pct": 100.0 * sum(row["within_ee"] for row in selected) / len(selected),
            "rmse": float(np.sqrt(np.mean(np.square(errors)))),
            "bias": float(np.mean(errors)),
            "consensus_count": sum(row["used_consensus"] for row in selected),
            "by_role": by_role,
        },
        selected,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tag")
    parser.add_argument("--manifest", type=Path, default=ROOT / "surface_e1_matched_manifest.json")
    parser.add_argument("--output", type=Path, default=ROOT / "band_consensus_aod_screen.json")
    args = parser.parse_args()

    metadata = {row["matchup_id"]: row for row in json.loads(args.manifest.read_text())["rows"]}
    records = _load_tag(args.tag)
    rows = [
        {"matchup_id": matchup_id, "role": metadata[matchup_id]["role"], "record": record}
        for matchup_id, record in records.items()
        if matchup_id in metadata
    ]
    configs = [
        ("current", None, None, None, None),
        ("spread020_median", 0.20, None, None, None),
        ("spread030_median", 0.30, None, None, None),
        ("flat005_spread020_median", 0.20, 0.005, None, None),
        ("flat010_spread020_median", 0.20, 0.010, None, None),
        ("agree010_deviation010_median", None, None, 0.10, 0.10),
    ]
    summaries = []
    strategy_rows = {}
    for name, spread, flat, aggregation_spread, aggregation_deviation in configs:
        summary, selected = _summary(
            name,
            rows,
            spread_threshold=spread,
            require_flat_delta=flat,
            aggregation_spread_max=aggregation_spread,
            aggregation_deviation_min=aggregation_deviation,
        )
        summaries.append(summary)
        strategy_rows[name] = selected
    baseline = {row["matchup_id"]: row for row in strategy_rows["current"]}
    for summary in summaries[1:]:
        selected = strategy_rows[summary["strategy"]]
        summary["gains"] = [
            row["matchup_id"]
            for row in selected
            if not baseline[row["matchup_id"]]["within_ee"] and row["within_ee"]
        ]
        summary["losses"] = [
            row["matchup_id"]
            for row in selected
            if baseline[row["matchup_id"]]["within_ee"] and not row["within_ee"]
        ]
    result = {
        "tag": args.tag,
        "available": len(rows),
        "summaries": summaries,
        "rows": strategy_rows,
    }
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
