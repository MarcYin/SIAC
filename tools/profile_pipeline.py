"""Profile a SIAC scene end-to-end with cProfile and emit a flamegraph-friendly summary.

This is the harness that REVIEW_FIXES.md called out as the next-step for
the deferred performance work. Run it against the same SAFE/config the
regression suite uses to find which of the REVIEW.md §1.2 / §3.x
performance hot-spots actually matter.

Usage::

    PYTHONPATH=python pixi run -e rt6s python tools/profile_pipeline.py \
        --safe tmp/.../S2B_..._T33KWP_....SAFE \
        --config tmp/real_cdse_mcd43_t33kwp_sixs.toml \
        --output-dir /tmp/profile_run \
        --report-path /tmp/profile.txt \
        --top 50

The script writes:

- ``<report-path>``: the formatted top-N by cumulative time.
- ``<report-path>.pstats``: the raw cProfile.Stats binary (for
  ``snakeviz`` or ``gprof2dot``).
- ``<report-path>.callers.txt``: top-N callers for the top hotspots
  (helpful when a hot leaf function has many call sites).

The `cProfile` overhead is non-trivial (~10-20% wall-clock) but
shape-stable enough that the relative rankings are reliable.
"""

from __future__ import annotations

import argparse
import cProfile
import pstats
import sys
from pathlib import Path


def _run_pipeline(safe: Path, config_path: Path, out: Path) -> None:
    """Run the same pipeline path the CLI / regression suite uses."""
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
    process_scene(request)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--safe", required=True, type=Path, help="Path to the .SAFE directory.")
    parser.add_argument("--config", required=True, type=Path, help="Path to the SIAC config TOML.")
    parser.add_argument(
        "--output-dir", required=True, type=Path, help="Where the pipeline writes BOA/AOT/TCWV outputs."
    )
    parser.add_argument(
        "--report-path",
        type=Path,
        default=Path("profile_report.txt"),
        help="Path to write the formatted top-N report (default: profile_report.txt).",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=40,
        help="Top-N functions to report by cumulative time (default: 40).",
    )
    parser.add_argument(
        "--sort",
        default="cumulative",
        choices=("cumulative", "tottime", "calls", "percall"),
        help="pstats sort key (default: cumulative).",
    )
    args = parser.parse_args(argv)

    if not args.safe.exists():
        parser.error(f"--safe does not exist: {args.safe}")
    if not args.config.exists():
        parser.error(f"--config does not exist: {args.config}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.report_path.parent.mkdir(parents=True, exist_ok=True)

    profiler = cProfile.Profile()
    profiler.enable()
    try:
        _run_pipeline(args.safe, args.config, args.output_dir)
    finally:
        profiler.disable()

    # Aux files share a consistent suffix-chain pattern so the trio is
    # easy to glob: ``<report>``, ``<report>.pstats``, ``<report>.callers.txt``.
    # Previously the pstats path was suffix-chained but the callers path
    # used ``with_name(stem + ".callers.txt")`` which produced
    # ``<stem>.callers.txt`` instead of ``<stem>.<ext>.callers.txt`` —
    # the verify_review_fixes.sh harness then couldn't find the callers
    # report and reported FAIL.
    pstats_path = args.report_path.with_suffix(args.report_path.suffix + ".pstats")
    profiler.dump_stats(str(pstats_path))

    # Top-N by chosen sort key.
    with args.report_path.open("w") as fh:
        stats = pstats.Stats(profiler, stream=fh)
        stats.strip_dirs().sort_stats(args.sort).print_stats(args.top)

    # Top-N callers for the same set (helps disambiguate hot leaves).
    callers_path = args.report_path.with_suffix(args.report_path.suffix + ".callers.txt")
    with callers_path.open("w") as fh:
        stats = pstats.Stats(profiler, stream=fh)
        stats.strip_dirs().sort_stats(args.sort).print_callers(args.top)

    print(f"Wrote {args.report_path}", file=sys.stderr)
    print(f"Wrote {pstats_path} (load with `python -m pstats {pstats_path}`)", file=sys.stderr)
    print(f"Wrote {callers_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
