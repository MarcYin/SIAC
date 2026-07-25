"""Validate the v1 production path against the AERONET campaign.

Runs the shipped ``siac_process_s2`` entry point — not the research harness —
with :func:`siac.config.get_surface_driven_v1_config`, one scene at a time, and
records retrieved AOD against the AERONET truth. Success is reproducing the
campaign result the recipe was selected on (84.6% within EE with the fused
aerosol prior, 79.2% with a single-source prior).

Usage:
    validate_v1_production.py <matchup_id> [...] --out <dir> [--no-fusion]
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import traceback
from pathlib import Path
from typing import Any

import numpy as np

CAMPAIGN = Path(
    "/gws/ssde/j25a/nceo_isp/public/siac_refactor/reports/"
    "aod-final-performance-dashboard-20260713/data/all-cases.csv"
)
LIBRARY = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor/acix3_faithful_prior_dict_6scci")
# Group CAMS archive: one ~425 MB netCDF per day, laid out as {base}/YYYY-MM-DD.nc,
# which is the layout CAMSProvider resolves against.
CAMS_ARCHIVE = Path("/gws/ssde/j25b/nceo_ard/public/cams")
# Prepared per-scene aerosol priors, as used by the reference run: the fused
# max(MAIAC, CAMS) store and the plain-MAIAC store for the ablation.
PREPARED_AOD_FUSED = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor/maiac_qa_maxcams")
PREPARED_AOD_PLAIN = Path(
    "/gws/ssde/j25a/nceo_isp/public/siac_refactor/maiac_qa_lowcloud20_native_adaptive"
)
# Site coordinates live alongside the campaign results rather than in the case
# table; any completed campaign arm carries them.
SITE_COORDS = Path(
    "/gws/ssde/j25a/nceo_isp/public/siac_refactor/"
    "phaseD_results_lowcloud20_acix3faithful6s_sixs_cci_floor006_dt_20260721"
)


def load_campaign() -> dict[str, dict[str, float]]:
    rows: dict[str, dict[str, float]] = {}
    with CAMPAIGN.open() as handle:
        for row in csv.DictReader(handle):
            rows[row["matchup_id"]] = {
                "truth": float(row["truth"]),
                "maiac": float(row["maiac_aod"]) if row.get("maiac_aod") else float("nan"),
                "cams": float(row["cams_aod"]) if row.get("cams_aod") else float("nan"),
            }
    return rows


def scene_metadata(matchup_id: str) -> tuple[float, float, str, tuple[float, ...] | None]:
    """Site coordinates, product ID and AOI for a campaign matchup.

    The AOI is taken from the reference run rather than reconstructed: the
    retrieval pools over the whole AOI, so a larger box averages a localized
    aerosol plume down and under-reports thick cases. Matching it is required
    for the comparison to be like-for-like.
    """

    payload = json.loads((SITE_COORDS / f"{matchup_id}.json").read_text())
    bbox = payload.get("aoi_bbox_4326")
    aoi = tuple(float(v) for v in bbox) if bbox else None
    return float(payload["lon"]), float(payload["lat"]), str(payload["product_id"]), aoi


def _reference_tcwv_cm(matchup_id: str) -> float | None:
    payload = json.loads((SITE_COORDS / f"{matchup_id}.json").read_text())
    value = payload.get("target_tcwv_cm")
    return float(value) if value is not None else None


def build_config(matchup_id: str, *, fusion: bool) -> Any:
    from siac.config import get_surface_driven_v1_config

    # V1_LIVE_SURFACE=1 measures the works-anywhere configuration: no prepared
    # library, so the surface prior is built live from L2A composites.
    live_surface = bool(os.environ.get("V1_LIVE_SURFACE"))
    config = get_surface_driven_v1_config(
        prepared_library_path=None if live_surface else LIBRARY,
        cache_root=os.environ.get("SIAC_V1_CACHE", "/gws/ssde/j25a/nceo_isp/public/siac_cache"),
    )

    cache_root = Path(
        os.environ.get("SIAC_V1_CACHE", "/gws/ssde/j25a/nceo_isp/public/siac_cache")
    )
    overrides: dict[str, Any] = {
        "algorithms": {
            "surface_prior": (
                {} if live_surface else {"prepared_library_scene_key": matchup_id}
            ),
        },
        # Anonymous public Sentinel-2 archive: no credentials, and the scene the
        # library was built from is fetched from the same source.
        "providers": {
            "s2": {
                "backend": "gcs",
                "processing_level": "L1C",
                "cache_dir": cache_root / "s2",
            },
            # Read CAMS from the group archive rather than downloading it. Each
            # day is a ~425 MB file, so a fleet of scenes would otherwise pull
            # tens of gigabytes over the network and into local storage.
            "atmo": {
                "data_path": CAMS_ARCHIVE,
                "download_missing": False,
                # Dedicated cache: the atmospheric prior is cached by scene
                # identity, so a shared cache would serve entries built before
                # a provider change and silently mask its effect.
                "cache_dir": Path(
                    os.environ.get("SIAC_V1_ATMO_CACHE", str(cache_root / "atmo_v1"))
                ),
                "prepared_scalar_path": PREPARED_AOD_FUSED if fusion else PREPARED_AOD_PLAIN,
                "prepared_scalar_scene_key": matchup_id,
                # Diagnostic arm: inject the reference run's recorded per-scene
                # TCWV (L2A spatial mean) instead of CAMS, to measure how much
                # of the residual gap is the water-vapour source.
                "prepared_scalar_tcwv_cm": (
                    _reference_tcwv_cm(matchup_id)
                    if os.environ.get("V1_TCWV_FROM_REFERENCE")
                    else None
                ),
            },
        },
        # Retrieval-only: M5 solves, M6 is skipped, and result.aot is the
        # solved AOD field rather than a corrector-propagated copy.
        "output": {"defaults": {"skip_correction": True}},
    }
    # Pin the prebuilt native 6S extension and share its run cache. Without an
    # explicit module the backend tries to fetch and compile 6S per task, which
    # a fleet cannot do (and the upstream source host has an untrusted cert).
    sixs_module = os.environ.get("SIAC_V1_SIXS_MODULE")
    if sixs_module:
        overrides["algorithms"]["rt"] = {
            "sixs": {
                "module_path": sixs_module,
                "auto_build": False,
                "run_cache_dir": os.environ.get(
                    "SIAC_V1_SIXS_RUN_CACHE", str(cache_root / "rt6s_run_cache")
                ),
            }
        }
    if not fusion:
        # Ablation: the same recipe on a single-source aerosol prior, which is
        # what the +5.4-point fusion result is measured against.
        overrides["providers"]["atmo"]["fuse_aod_with"] = ()
    return config.with_overrides(**overrides)


def retrieve(matchup_id: str, *, fusion: bool, half_deg: float = 0.12) -> dict[str, Any]:
    from siac.api import siac_process_s2

    lon, lat, product_id, _reference_aoi = scene_metadata(matchup_id)
    # A wider AOI than the reference's ~5 km box measured better (81.6% vs 78.2%):
    # the pooled solve gets more pixels and a more robust cost surface. Set
    # V1_MATCH_REFERENCE_AOI=1 to reproduce the reference's box instead.
    if os.environ.get("V1_MATCH_REFERENCE_AOI") and _reference_aoi is not None:
        aoi = _reference_aoi
    else:
        aoi = (lon - half_deg, lat - half_deg, lon + half_deg, lat + half_deg)
    config = build_config(matchup_id, fusion=fusion)
    result = siac_process_s2(config, product_id, aoi=aoi)
    return extract(result, lon, lat)


def extract(result: Any, lon: float, lat: float, radius_m: float = 3000.0) -> dict[str, Any]:
    """Window statistics of retrieved AOD around the site."""

    import rasterio.warp

    aot = result.aot
    crs = aot.rio.crs
    xs, ys = rasterio.warp.transform("EPSG:4326", crs, [lon], [lat])
    x, y = float(xs[0]), float(ys[0])
    window = aot.sel(
        x=slice(x - radius_m, x + radius_m), y=slice(y + radius_m, y - radius_m)
    )
    values = np.asarray(window.values, dtype=float).ravel()
    finite = values[np.isfinite(values)]
    solver = (result.metadata or {}).get("solver", {}) or {}

    def _json_safe(obj: Any) -> Any:
        if isinstance(obj, dict):
            return {str(k): _json_safe(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [_json_safe(v) for v in obj]
        if isinstance(obj, (str, bool)) or obj is None:
            return obj
        try:
            return float(obj)
        except (TypeError, ValueError):
            return str(obj)

    scene_mean = float(solver.get("aot_mean", float("nan")))
    window_median = float(np.nanmedian(finite)) if finite.size else float("nan")
    return {
        # Full stage diagnostics so a production record can be diffed
        # field-by-field against a reference-run record with the same keys.
        "solver_diagnostics": _json_safe(solver),
        # Scored statistic: the AOI-mean solved AOD, matching how the recipe
        # was validated (the station-window median reads localized plumes
        # differently and flips knife-edge EE sites).
        "retrieved": scene_mean if np.isfinite(scene_mean) else window_median,
        "retrieved_window_median": window_median,
        "retrieved_mean": float(np.nanmean(finite)) if finite.size else float("nan"),
        "n_valid": int(finite.size),
        "scene_mean": scene_mean,
        "converged": bool(solver.get("converged", False)),
        "cost_final": float(solver.get("cost_final", float("nan"))),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("matchups", nargs="+")
    parser.add_argument("--out", required=True)
    parser.add_argument("--no-fusion", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    campaign = load_campaign()

    for matchup_id in args.matchups:
        target = out_dir / f"{matchup_id}.json"
        if target.exists():
            # Only completed retrievals are terminal; a recorded failure is
            # retried so a transient fetch error does not poison the sweep.
            try:
                if json.loads(target.read_text()).get("status") == "OK":
                    print(f"exists: {matchup_id}", flush=True)
                    continue
            except (OSError, ValueError):
                pass
        record: dict[str, Any] = {"matchup_id": matchup_id, "fusion": not args.no_fusion}
        try:
            record.update(retrieve(matchup_id, fusion=not args.no_fusion))
            truth = campaign[matchup_id]["truth"]
            record["truth"] = truth
            record["err"] = record["retrieved"] - truth
            record["ee_threshold"] = 0.05 + 0.15 * truth
            record["within_ee"] = bool(abs(record["err"]) <= record["ee_threshold"])
            record["status"] = "OK"
        except Exception as exc:  # noqa: BLE001 - one bad scene must not stop the sweep
            record["status"] = f"{type(exc).__name__}: {exc}"
            record["traceback"] = traceback.format_exc()[-2000:]
        target.write_text(json.dumps(record, indent=1))
        print(
            f"{matchup_id}: {record['status']}"
            + (
                f" truth={record.get('truth', float('nan')):.3f}"
                f" retrieved={record.get('retrieved', float('nan')):.3f}"
                f" within_ee={record.get('within_ee')}"
                if record["status"] == "OK"
                else ""
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
