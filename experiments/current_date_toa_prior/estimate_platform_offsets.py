"""Measure the S2B/S2C -> S2A visible offset used by the hybrid label build.

The joint multisensor model predicts canonical S2A reflectance. The same-date
M5 teacher that supervises its visible bands sits in its own spacecraft's band
space, so an S2B or S2C scene needs a small spectral correction before its
teacher can be used as an S2A target.

The correction is measured, not assumed, and it is measured once over a large
scene sample rather than per scene. Bridging a single scene's spectra puts the
answer wherever the k-nearest-neighbour dictionary happens to extrapolate: on a
dark scene that produced -0.0099 in B02 against a population value near
+0.002 -- a 40% correction on a surface of reflectance 0.024. Pooling across
scenes removes that sensitivity, and the scene-to-scene spread reported here is
what the builder folds into the label sigma.

Only month-matched library spectra are read. No AERONET measurement is used.
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np

BANDS = ("B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B09", "B11", "B12")
LAND_BANDS = tuple(b for b in BANDS if b != "B09")
VISIBLE = ("B02", "B03", "B04")
VISIBLE_LAND = [LAND_BANDS.index(b) for b in VISIBLE]
#: Only well-reconstructed spectra inform the offset; fit_teacher's own label
#: gate is 0.03, but an offset applied to thousands of scenes is held tighter.
MAX_BRIDGE_FIT = 0.01
MIN_ROWS = 200


def scene_offset(library_path: Path, sensor: str, rows: int, rng: np.random.Generator):
    """Median visible offset for one scene, or ``None`` where it cannot be measured."""
    from experiments.current_date_toa_prior.landsat_spectral_bridge import bridge_surface_rows

    with np.load(library_path, allow_pickle=False) as archive:
        names = [str(v) for v in archive["band_names"]]
        support = archive["comp"][:, [names.index(b) for b in LAND_BANDS]]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        spectra = np.nanmedian(support, axis=0).reshape(len(LAND_BANDS), -1).T.astype(np.float32)
    spectra = spectra[np.isfinite(spectra).all(-1)]
    if spectra.shape[0] < MIN_ROWS * 4:
        return None
    spectra = spectra[rng.choice(spectra.shape[0], min(rows, spectra.shape[0]), replace=False)]
    bridged = bridge_surface_rows(
        spectra,
        source_platform=sensor,
        source_band_names=LAND_BANDS,
        target_sensor="S2A",
        target_band_names=LAND_BANDS,
        source_uncertainty=np.full(spectra.shape, 0.01, dtype=np.float32),
    )
    canonical = np.asarray(bridged["reflectance"], dtype=np.float32)
    fit = np.asarray(bridged["source_fit_rmse"], dtype=np.float32)
    good = np.isfinite(canonical).all(-1) & np.isfinite(fit) & (fit <= MAX_BRIDGE_FIT)
    if good.sum() < MIN_ROWS:
        return None
    return np.median(canonical[good][:, VISIBLE_LAND] - spectra[good][:, VISIBLE_LAND], axis=0)


def run(args: argparse.Namespace) -> int:
    release = json.loads(Path(args.release).read_text())
    result: dict[str, dict] = {}
    for sensor in args.sensors:
        records = [
            r
            for r in release["records"]
            if r.get("sensor") == sensor
            and "expanded_library" in r.get("inputs", {})
            and r.get("split") != "holdout"
        ]
        if not records:
            continue
        rng = np.random.default_rng(args.seed)
        chosen = [records[i] for i in rng.choice(len(records), min(args.scenes, len(records)), replace=False)]
        per_scene = []
        for record in chosen:
            try:
                value = scene_offset(
                    Path(record["inputs"]["expanded_library"]["path"]), sensor, args.rows, rng
                )
            except Exception as exc:  # noqa: BLE001 - one library must not sink the estimate
                print(f"skip {record['matchup_id']}: {type(exc).__name__}: {exc}", file=sys.stderr)
                continue
            if value is not None:
                per_scene.append(value)
        if len(per_scene) < 20:
            raise ValueError(f"Too few usable {sensor} scenes to fix a global offset")
        values = np.asarray(per_scene, dtype=np.float32)
        result[sensor] = {
            "scenes": int(values.shape[0]),
            "bands": list(VISIBLE),
            "per_scene_median": np.median(values, axis=0).tolist(),
            "per_scene_sd": values.std(axis=0, ddof=1).tolist(),
            "per_scene_p5": np.percentile(values, 5, axis=0).tolist(),
            "per_scene_p95": np.percentile(values, 95, axis=0).tolist(),
            "max_bridge_fit": MAX_BRIDGE_FIT,
            "rows_per_scene": args.rows,
        }
        print(f"\n=== {sensor} -> S2A  ({values.shape[0]} scenes) ===")
        for i, band in enumerate(VISIBLE):
            print(
                f"  {band}: median {np.median(values[:, i]):+.5f}  sd {values[:, i].std(ddof=1):.5f}"
                f"   5-95% [{np.percentile(values[:, i], 5):+.5f}, {np.percentile(values[:, i], 95):+.5f}]"
            )
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=1))
    print(f"\nwrote {out}")
    return 0


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--release", required=True)
    value.add_argument("--out", required=True)
    value.add_argument("--sensors", nargs="+", default=["S2B", "S2C"])
    value.add_argument("--scenes", type=int, default=140)
    value.add_argument("--rows", type=int, default=3000)
    value.add_argument("--seed", type=int, default=11)
    return value


def main() -> None:
    sys.exit(run(parser().parse_args()))


if __name__ == "__main__":
    main()
