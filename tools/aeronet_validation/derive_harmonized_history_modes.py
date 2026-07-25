"""Derive band-limited histories from aligned identity and all-band archives."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_HISTORY_ROOT = (
    ROOT
    / "analysis/l2a_l1c_harmonizer_clean015_mediumdev_20260714/daily_histories"
)
BAND_INDICES = {
    "visible": (0, 1, 2, 3),
    "solver": (1, 2, 3),
}


def _archive(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as source:
        return {name: np.asarray(source[name]) for name in source.files}


def derive(
    identity_path: Path,
    mapped_path: Path,
    output_path: Path,
    *,
    mode: str,
) -> None:
    if mode not in BAND_INDICES:
        raise ValueError(f"unknown mode {mode!r}")
    identity = _archive(identity_path)
    mapped = _archive(mapped_path)
    if identity["comp"].shape != mapped["comp"].shape:
        raise ValueError("identity and mapped comp shapes differ")
    if identity["comp"].ndim != 4 or identity["comp"].shape[1] != 7:
        raise ValueError("expected comp shape (realization, 7, y, x)")
    for name in ("epsg", "transform", "realizations", "scene_year", "month"):
        if name not in identity or name not in mapped:
            raise ValueError(f"missing aligned history field {name}")
        if not np.array_equal(identity[name], mapped[name]):
            raise ValueError(f"identity and mapped {name} differ")
    composite = np.asarray(identity["comp"], dtype=np.float32).copy()
    indices = BAND_INDICES[mode]
    composite[:, indices] = np.asarray(mapped["comp"], dtype=np.float32)[:, indices]
    provenance: dict[str, Any] = {
        "source_identity": str(identity_path),
        "source_all_band_mapping": str(mapped_path),
        "mode": mode,
        "mapped_band_indices": list(indices),
        "derivation": "band substitution after independent per-band temporal compositing",
        "uses_aeronet": False,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(".tmp.npz")
    np.savez_compressed(
        temporary,
        comp=composite,
        epsg=identity["epsg"],
        transform=identity["transform"],
        realizations=identity["realizations"],
        scene_year=identity["scene_year"],
        month=identity["month"],
        harmonization_json=np.asarray(json.dumps(provenance)),
    )
    temporary.replace(output_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matchup_id", nargs="+")
    parser.add_argument("--history-root", type=Path, default=DEFAULT_HISTORY_ROOT)
    parser.add_argument("--mapped-tag", default="full_a100_all_cap0p030")
    parser.add_argument("--mode", choices=tuple(BAND_INDICES), required=True)
    parser.add_argument("--output-tag")
    args = parser.parse_args()
    output_tag = args.output_tag or f"full_a100_{args.mode}_cap0p030"
    for matchup_id in args.matchup_id:
        derive(
            args.history_root / "identity_daily" / f"{matchup_id}.npz",
            args.history_root / args.mapped_tag / f"{matchup_id}.npz",
            args.history_root / output_tag / f"{matchup_id}.npz",
            mode=args.mode,
        )
        print(f"DERIVED {matchup_id} {args.mode}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
