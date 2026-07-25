"""Merge per-scene exact-WVP 6S coefficients into an AtmoSidecar."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

DEFAULT_BANDS = ["coastal", "blue", "green", "red", "nir08", "swir16", "swir22"]
LEGACY_BAND_NAMES = {
    "ultra_blue": "coastal",
    "blue": "blue",
    "green": "green",
    "red": "red",
    "nir": "nir08",
    "swir1": "swir16",
    "swir2": "swir22",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--meta", required=True)
    parser.add_argument("--coeffs", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    metadata = json.loads(Path(args.meta).read_text(encoding="utf-8"))["selected"]
    with np.load(args.coeffs, allow_pickle=False) as archive:
        coefficients = archive["coeffs"]
        tcwv_by_day = archive["tcwv_by_day"]
        raw_bands = archive["bands"].tolist() if "bands" in archive.files else DEFAULT_BANDS
    bands = [LEGACY_BAND_NAMES.get(str(name), str(name)) for name in raw_bands]
    if coefficients.shape != (len(metadata), 2, 1, 3, len(bands)):
        raise ValueError(
            f"unexpected exact-WVP coefficient shape {coefficients.shape!r} "
            f"for {len(metadata)} scenes"
        )
    if tcwv_by_day.shape != (len(metadata),):
        raise ValueError(f"unexpected tcwv_by_day shape {tcwv_by_day.shape!r}")

    scenes = {}
    for index, day in enumerate(metadata):
        values = coefficients[index, 0, 0]
        tcwv = float(tcwv_by_day[index])
        scenes[day["l1c_id"]] = {
            "maiac_aod": float(day["maiac"]),
            "wvp": tcwv,
            "tcwv_nodes": [tcwv],
            "xap": [[float(value) for value in values[0]]],
            "xbp": [[float(value) for value in values[1]]],
            "xcp": [[float(value) for value in values[2]]],
        }
    Path(args.out).write_text(json.dumps({"bands": bands, "scenes": scenes}), encoding="utf-8")
    print(f"wrote {args.out}: {len(scenes)} exact-WVP scenes", flush=True)


if __name__ == "__main__":
    main()
