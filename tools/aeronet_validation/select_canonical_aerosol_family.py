"""Select among four saved canonical 6S aerosol-family retrievals."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from itertools import product
from pathlib import Path
from typing import Any

import numpy as np

from siac.algorithms.rt.aerosol_species import climatology_fraction_percentages

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
FAMILIES = ("continental", "biomass_burning", "desert", "maritime")


def family_priors_from_percentages(values: np.ndarray) -> dict[str, float]:
    """Map CCI component percentages to weak canonical-family probabilities."""

    dust, sea_salt, fine_strong, fine_weak = np.asarray(values, dtype=float) / 100.0
    raw = {
        "continental": fine_weak + 0.25 * dust + 0.10 * sea_salt,
        "biomass_burning": fine_strong,
        "desert": dust,
        "maritime": sea_salt,
    }
    floored = {family: max(0.02, float(value)) for family, value in raw.items()}
    total = sum(floored.values())
    return {family: value / total for family, value in floored.items()}


def selection_score(
    *,
    surface_cost: float,
    prior_probability: float,
    band_spread: float,
    curve_delta: float,
    prior_weight: float,
    spread_weight: float,
    flat_weight: float,
    edge_fraction: float = 0.0,
    edge_weight: float = 0.0,
    flat_target: float = 0.02,
) -> float:
    prior_term = -math.log(max(prior_probability, 1e-6))
    spread_term = band_spread if math.isfinite(band_spread) else 4.0
    flat_term = (
        max(0.0, flat_target - curve_delta) / flat_target if math.isfinite(curve_delta) else 1.0
    )
    return (
        surface_cost
        + prior_weight * prior_term
        + spread_weight * spread_term
        + flat_weight * flat_term
        + edge_weight * edge_fraction
    )


def band_edge_fraction(
    band_argmins: dict[str, float],
    *,
    lower: float = 0.01,
    upper: float = 4.0,
    tolerance: float = 1e-6,
) -> float:
    """Return the fraction of finite per-band minima clipped to the AOD axis."""

    values = [float(value) for value in band_argmins.values() if math.isfinite(float(value))]
    if not values:
        return 1.0
    return sum(value <= lower + tolerance or value >= upper - tolerance for value in values) / len(
        values
    )


def scale_surface_costs(costs: list[float], mode: str) -> list[float]:
    """Scale family costs without mixing evidence scales between scenes."""

    values = np.asarray(costs, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("surface costs must be finite")
    if mode == "raw":
        return values.tolist()
    if mode != "relative":
        raise ValueError(f"unsupported surface-cost mode {mode!r}")
    minimum = float(np.min(values))
    return ((values - minimum) / max(abs(minimum), 1e-6)).tolist()


def _load_tag(tag: str) -> dict[str, dict[str, Any]]:
    records = {}
    for path in (ROOT / "seasonal_val").glob(f"*__seas_maiac_{tag}.json"):
        row = json.loads(path.read_text())
        records[str(row["matchup_id"])] = row
    return records


def _diagnostic(row: dict[str, Any]) -> dict[str, Any]:
    candidates = row.get("species_candidates") or []
    if len(candidates) != 1:
        raise ValueError(f"expected one forced-family diagnostic for {row['matchup_id']}")
    return candidates[0]


def _fold(site: str, folds: int) -> int:
    return int.from_bytes(hashlib.sha256(site.encode()).digest()[:4], "big") % folds


def _within(aod: float, truth: float) -> bool:
    return abs(aod - truth) <= 0.05 + 0.15 * truth


def _summarize(rows: list[dict[str, Any]], folds: int) -> dict[str, Any]:
    errors = np.asarray([row["aod"] - row["truth"] for row in rows], dtype=float)
    by_role = {}
    for role in sorted({str(row["role"]) for row in rows}):
        selected = [row for row in rows if row["role"] == role]
        by_role[role] = {
            "hits": sum(row["within_ee"] for row in selected),
            "count": len(selected),
        }
    by_fold = []
    for fold in range(folds):
        selected = [row for row in rows if row["fold"] == fold]
        by_fold.append(
            {
                "fold": fold,
                "hits": sum(row["within_ee"] for row in selected),
                "count": len(selected),
            }
        )
    return {
        "hits": sum(row["within_ee"] for row in rows),
        "count": len(rows),
        "rmse": float(np.sqrt(np.mean(np.square(errors)))),
        "bias": float(np.mean(errors)),
        "by_role": by_role,
        "by_fold": by_fold,
        "family_counts": {
            family: sum(row["family"] == family for row in rows) for family in FAMILIES
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=ROOT / "surface_e1_matched_manifest.json")
    parser.add_argument("--continental-tag", required=True)
    parser.add_argument("--biomass-tag", required=True)
    parser.add_argument("--desert-tag", required=True)
    parser.add_argument("--maritime-tag", required=True)
    parser.add_argument("--output", type=Path, default=ROOT / "canonical_aerosol_selector_e4.json")
    parser.add_argument("--folds", type=int, default=5)
    args = parser.parse_args()

    tags = {
        "continental": args.continental_tag,
        "biomass_burning": args.biomass_tag,
        "desert": args.desert_tag,
        "maritime": args.maritime_tag,
    }
    records = {family: _load_tag(tag) for family, tag in tags.items()}
    manifest = json.loads(args.manifest.read_text())
    metadata = {row["matchup_id"]: row for row in manifest["rows"]}
    complete = [
        matchup_id
        for matchup_id in metadata
        if all(matchup_id in records[family] for family in FAMILIES)
    ]

    weight_grid = list(
        product(
            (0.0, 0.25, 0.5, 1.0, 2.0),
            (0.0, 0.5, 1.0, 2.0),
            (0.0, 0.25, 0.5),
            (0.0, 0.5, 1.0, 2.0),
        )
    )
    candidates = []
    candidate_rows: dict[str, list[dict[str, Any]]] = {}
    for cost_mode in ("raw", "relative"):
        for prior_weight, spread_weight, flat_weight, edge_weight in weight_grid:
            name = (
                f"{cost_mode}_prior{prior_weight:g}_spread{spread_weight:g}"
                f"_flat{flat_weight:g}_edge{edge_weight:g}"
            )
            selected_rows = []
            for matchup_id in complete:
                base = records["continental"][matchup_id]
                truth = float(base["truth"])
                site = str(base.get("site") or matchup_id.split("__")[0])
                matchup_parts = matchup_id.rsplit("_", 1)[-1]
                month = int(matchup_parts[4:6])
                if base.get("lon") is not None and base.get("lat") is not None:
                    percentages = climatology_fraction_percentages(
                        float(base["lon"]), float(base["lat"]), month
                    )
                    priors = family_priors_from_percentages(percentages)
                else:
                    priors = dict.fromkeys(FAMILIES, 0.25)
                raw_options = []
                for family in FAMILIES:
                    row = records[family][matchup_id]
                    diagnostic = _diagnostic(row)
                    edge_fraction = band_edge_fraction(diagnostic["band_argmin_aod"])
                    raw_options.append((family, row, diagnostic, edge_fraction))
                scaled_costs = scale_surface_costs(
                    [float(option[2]["surface_cost"]) for option in raw_options],
                    cost_mode,
                )
                options = []
                for scaled_cost, (family, row, diagnostic, edge_fraction) in zip(
                    scaled_costs, raw_options, strict=True
                ):
                    score = selection_score(
                        surface_cost=scaled_cost,
                        prior_probability=priors[family],
                        band_spread=float(diagnostic["band_argmin_spread"]),
                        curve_delta=float(diagnostic["curve_relative_second_delta"]),
                        prior_weight=prior_weight,
                        spread_weight=spread_weight,
                        flat_weight=flat_weight,
                        edge_fraction=edge_fraction,
                        edge_weight=edge_weight,
                    )
                    options.append((score, family, row, diagnostic, edge_fraction, scaled_cost))
                score, family, row, diagnostic, edge_fraction, scaled_cost = min(
                    options, key=lambda option: option[0]
                )
                aod = float(row["aod"])
                selected_rows.append(
                    {
                        "matchup_id": matchup_id,
                        "site": site,
                        "role": metadata[matchup_id]["role"],
                        "fold": _fold(site, args.folds),
                        "truth": truth,
                        "aod": aod,
                        "within_ee": _within(aod, truth),
                        "family": family,
                        "selection_score": score,
                        "surface_cost": diagnostic["surface_cost"],
                        "scaled_surface_cost": scaled_cost,
                        "band_spread": diagnostic["band_argmin_spread"],
                        "curve_delta": diagnostic["curve_relative_second_delta"],
                        "edge_fraction": edge_fraction,
                    }
                )
            summary = _summarize(selected_rows, args.folds)
            candidates.append(
                {
                    "candidate": name,
                    "cost_mode": cost_mode,
                    "prior_weight": prior_weight,
                    "spread_weight": spread_weight,
                    "flat_weight": flat_weight,
                    "edge_weight": edge_weight,
                    **summary,
                }
            )
            candidate_rows[name] = selected_rows

    candidates.sort(key=lambda row: (-row["hits"], row["rmse"], abs(row["bias"])))
    best = candidates[0]
    result = {
        "description": "Scene-level selection among four forced native-6S aerosol families",
        "tags": tags,
        "manifest_count": len(metadata),
        "complete_count": len(complete),
        "candidates": candidates,
        "best_candidate": best,
        "best_rows": candidate_rows[best["candidate"]],
    }
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(best, indent=2))
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
