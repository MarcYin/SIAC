#!/usr/bin/env python
"""Filter a built L1C seasonal prior to its globally cleanest windows."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tools.aeronet_validation.clean_day_sidecar import select_clean_windows


def filter_prior(
    source: Path,
    destination: Path,
    *,
    max_median_aod: float,
    min_windows: int,
) -> dict[str, object]:
    """Filter one prior while preserving its spatial and provenance fields."""

    with np.load(source, allow_pickle=False) as archive:
        payload = {name: archive[name] for name in archive.files}
    quality = json.loads(str(payload["clean_quality_json"].item()))
    composites = np.asarray(payload["comp"])
    realizations = np.asarray(payload["realizations"])
    if len(quality) != composites.shape[0] or realizations.shape[0] != composites.shape[0]:
        raise ValueError("clean quality, realizations, and composites are not aligned")
    keep, selection = select_clean_windows(
        quality,
        max_median_aod=max_median_aod,
        min_windows=min_windows,
    )
    payload["comp"] = composites[keep]
    payload["realizations"] = realizations[keep]
    payload["clean_quality_json"] = np.array(json.dumps([quality[index] for index in keep]))
    payload["clean_window_max_median_aod"] = np.array(max_median_aod)
    payload["clean_window_min"] = np.array(min_windows)
    payload["clean_window_selection_json"] = np.array(json.dumps(selection))
    destination.parent.mkdir(parents=True, exist_ok=True)
    np.savez(destination, **payload)
    return {
        "source": str(source),
        "destination": str(destination),
        "kept_realizations": [str(value) for value in realizations[keep]],
        "selection": selection,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--max-median-aod", type=float, default=0.15)
    parser.add_argument("--min-windows", type=int, default=3)
    args = parser.parse_args()
    result = filter_prior(
        args.source,
        args.destination,
        max_median_aod=args.max_median_aod,
        min_windows=args.min_windows,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
