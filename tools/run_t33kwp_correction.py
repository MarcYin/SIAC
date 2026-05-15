"""Run a real Sentinel-2 atmospheric correction against the cached T33KWP scene.

This is the "give me real code to verify it works" companion script. It
uses the same SAFE input, config TOML, and auxiliary caches that the
regression suite uses, runs the full M1-M6 pipeline through the 6S RT
backend, and prints a one-screen summary plus a goldens-diff so you
can immediately see whether the outputs match the captured baseline.

Usage::

    PYTHONPATH=python pixi run -e rt6s python tools/run_t33kwp_correction.py
    PYTHONPATH=python pixi run -e rt6s python tools/run_t33kwp_correction.py \
        --output-path /tmp/my_run

By default outputs go to a fresh ``/tmp/siac_run_t33kwp_<timestamp>``
directory so consecutive runs don't trample one another.

Wall-clock: ~6-8 minutes on an M-series Mac. The first 30 seconds is
preprocessing + atmosphere prior; the bulk is M5 solver + M6
correction.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path

# ---- inputs (matches the regression suite's defaults) --------------------
DEFAULT_SAFE = Path(
    "/Users/fengyin/Documents/SIAC/tmp/real_cdse_mcd43_t33kwp/cache/s2/"
    "S2B_MSIL1C_20260329T084559_N0512_R107_T33KWP_20260329T140503.SAFE"
)
DEFAULT_CONFIG = Path(
    "/Users/fengyin/Documents/SIAC/tmp/real_cdse_mcd43_t33kwp_sixs.toml"
)
DEFAULT_GOLDEN = (
    Path(__file__).resolve().parent.parent
    / "tests"
    / "regression"
    / "goldens"
    / "t33kwp_sixs_20260329.json"
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--safe",
        type=Path,
        default=DEFAULT_SAFE,
        help=f"Path to the .SAFE directory (default: {DEFAULT_SAFE}).",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help=f"Path to the SIAC config TOML (default: {DEFAULT_CONFIG}).",
    )
    parser.add_argument(
        "--output-path",
        type=Path,
        default=None,
        help="Output directory. Default: /tmp/siac_run_t33kwp_<UTC-timestamp>.",
    )
    parser.add_argument(
        "--compare-to-golden",
        type=Path,
        default=DEFAULT_GOLDEN,
        help=f"Captured goldens to compare against (default: {DEFAULT_GOLDEN}).",
    )
    parser.add_argument(
        "--no-golden-diff",
        action="store_true",
        help="Skip the goldens diff at the end.",
    )
    return parser.parse_args()


def _run_pipeline(safe: Path, config_path: Path, out: Path) -> object:
    """Same code path the CLI and the regression test use."""
    from siac.adapters.auth import CredentialManager
    from siac.api.requests import SceneProcessRequest
    from siac.config import SIACConfig
    from siac.workflows.scene import process_scene

    config = SIACConfig.from_file(config_path)
    request = SceneProcessRequest(
        config=config,
        input_path=safe,
        output_path=out,
        auth=CredentialManager.from_config(config),
    )
    return process_scene(request)


def _summarise_outputs(out_dir: Path) -> dict[str, dict[str, float]]:
    """Capture the same summary stats the goldens use."""
    import numpy as np
    import rasterio

    aot_files = sorted(out_dir.glob("*_AOT.tif"))
    if not aot_files:
        return {}
    prefix = aot_files[0].stem.removesuffix("_AOT")
    products = ["AOT", "TCWV", "CLOUD"]
    for band in ("B02", "B03", "B04", "B08"):
        products.append(f"BOA_{band}")

    stats: dict[str, dict[str, float]] = {}
    for label in products:
        path = out_dir / f"{prefix}_{label}.tif"
        if not path.exists():
            continue
        with rasterio.open(path) as ds:
            arr = ds.read(1, masked=True)
            finite = np.ma.compressed(arr)
            if finite.size == 0:
                stats[label] = {"all_masked": 1.0}
                continue
            stats[label] = {
                "valid_fraction": float(finite.size / arr.size),
                "mean": float(np.mean(finite)),
                "std": float(np.std(finite)),
                "p50": float(np.percentile(finite, 50)),
                "min": float(np.min(finite)),
                "max": float(np.max(finite)),
            }
    return stats


def _diff_against_golden(stats: dict[str, dict[str, float]], golden_path: Path) -> int:
    """Print a side-by-side diff vs the captured goldens.

    Returns the number of products that drifted past tolerance (1e-3 rel
    on the main stats). Zero means a clean run.
    """
    import numpy as np

    golden = json.loads(golden_path.read_text())
    print()
    print(f"  {'product':<10}  {'stat':<8}  {'golden':>14}  {'actual':>14}  delta")
    print(f"  {'-'*10}  {'-'*8}  {'-'*14}  {'-'*14}  {'-'*8}")
    drift = 0
    for product, current in stats.items():
        if product not in golden["products"]:
            print(f"  {product:<10}  (no golden)")
            continue
        gprod = golden["products"][product]
        for key in ("mean", "std", "p50"):
            if key not in current or key not in gprod:
                continue
            g, a = float(gprod[key]), float(current[key])
            close = np.isclose(a, g, rtol=1e-3, atol=1e-4)
            tag = " " if close else "!"
            if not close:
                drift += 1
            print(
                f"  {product:<10}  {key:<8}  {g:>14.6f}  {a:>14.6f}  "
                f"{a-g:+.3e} {tag}"
            )
    return drift


def main() -> int:
    args = _parse_args()

    if not args.safe.exists():
        print(f"FATAL: SAFE input does not exist: {args.safe}", file=sys.stderr)
        return 2
    if not args.config.exists():
        print(f"FATAL: config does not exist: {args.config}", file=sys.stderr)
        return 2

    out_path = args.output_path
    if out_path is None:
        stamp = datetime.utcnow().strftime("%Y%m%dT%H%M%S")
        out_path = Path(f"/tmp/siac_run_t33kwp_{stamp}")
    out_path.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print(f"  SAFE   : {args.safe}")
    print(f"  config : {args.config}")
    print(f"  output : {out_path}")
    print("=" * 72)

    t0 = time.perf_counter()
    try:
        result = _run_pipeline(args.safe, args.config, out_path)
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        return 130
    elapsed = time.perf_counter() - t0

    print()
    print(f"Pipeline completed in {elapsed:.1f}s ({elapsed/60:.1f} min).")
    print(f"AOT  mean: {float(result.aot.mean()):.5f}")
    print(f"TCWV mean: {float(result.tcwv.mean()):.5f}")

    # Files
    n = sum(1 for _ in out_path.iterdir())
    print(f"\nOutputs in {out_path} ({n} entries):")
    for path in sorted(out_path.glob("*.tif")):
        size_mb = path.stat().st_size / (1024 * 1024)
        print(f"  {path.name:<60}  {size_mb:>7.1f} MB")

    # Compare to goldens
    if not args.no_golden_diff and args.compare_to_golden.exists():
        print(f"\nDiff vs goldens at {args.compare_to_golden}:")
        stats = _summarise_outputs(out_path)
        drift = _diff_against_golden(stats, args.compare_to_golden)
        if drift == 0:
            print("  All compared stats within 1e-3 rel / 1e-4 abs tolerance ✓")
            return 0
        else:
            print(f"  {drift} stat(s) drifted past tolerance — see ! markers above.")
            return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
