"""CLI for the AERONET validation experiment.

Usage (from the repository root, inside the pixi environment):

    PYTHONPATH=.:python pixi run python -m tools.aeronet_validation.cli <stage> [options]

Stages: fetch-aeronet, matchup, build-manifest, run, compare, status, make-slurm.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from tools.aeronet_validation import compare, fetch_aeronet, matchup, run_retrieval
from tools.aeronet_validation.common import (
    APPROACHES,
    DEFAULT_DATA_ROOT,
    ExperimentPaths,
    setup_logging,
)

SBATCH_TEMPLATE = """#!/bin/bash
#SBATCH --job-name=siac-aeronet
#SBATCH --partition={partition}
#SBATCH --qos={qos}
{account_line}#SBATCH --time={time_limit}
#SBATCH --mem={memory}
#SBATCH --cpus-per-task={cpus}
#SBATCH --array=0-{last_index}%{max_concurrent}
#SBATCH --output={log_dir}/%A_%a.out

set -euo pipefail
cd {repo_root}
export PYTHONPATH={repo_root}:{repo_root}/python
exec pixi run python -m tools.aeronet_validation.cli run \\
    --data-root {data_root} \\
    --task-index "$SLURM_ARRAY_TASK_ID" {run_args}
"""


def _add_status_command(subparsers) -> None:
    subparsers.add_parser("status", help="Show experiment progress per stage.")


def _add_make_slurm_command(subparsers) -> None:
    parser = subparsers.add_parser(
        "make-slurm", help="Write a LOTUS sbatch array script for the run manifest."
    )
    parser.add_argument("--partition", default="standard")
    parser.add_argument(
        "--qos",
        default="high",
        help="JASMIN LOTUS QOS. standard/short/long cap jobs at cpu=1 per node, so "
        "multi-CPU retrieval tasks need 'high' (per-user cap cpu=576, mem=4500G).",
    )
    parser.add_argument("--account", default="nceo_isp", help="SLURM account (JASMIN project).")
    parser.add_argument("--time-limit", default="04:00:00")
    parser.add_argument("--memory", default="32G")
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--max-concurrent", type=int, default=20)
    parser.add_argument(
        "--run-args",
        default="",
        help="Extra arguments forwarded to the run stage (quoted string).",
    )


def _status(paths: ExperimentPaths) -> None:
    print(f"Data root: {paths.root}")
    if paths.sites_with_data_file.exists():
        sites = pd.read_csv(paths.sites_with_data_file)
        print(f"fetch-aeronet: {len(sites)} sites with data")
    else:
        print("fetch-aeronet: not run")
    if paths.matchups_file.exists():
        matchups = pd.read_csv(paths.matchups_file)
        print(f"matchup: {len(matchups)} matchups over {matchups['site'].nunique()} sites")
    else:
        print("matchup: not run")
    for approach in sorted(APPROACHES):
        approach_dir = paths.runs_dir / approach
        if not approach_dir.exists():
            print(f"run[{approach}]: not started")
            continue
        results = list(approach_dir.glob("*/result.json"))
        failed = sum(1 for r in results if '"status": "failed"' in r.read_text())
        print(f"run[{approach}]: {len(results)} done ({failed} failed)")
    summary = paths.results_dir / "summary.csv"
    if summary.exists():
        print(f"compare: {summary}")
        print(pd.read_csv(summary).round(4).to_string(index=False))
    else:
        print("compare: not run")


def _make_slurm(args: argparse.Namespace, paths: ExperimentPaths) -> None:
    if not paths.manifest_file.exists():
        raise SystemExit("No manifest; run the build-manifest stage first.")
    manifest = pd.read_csv(paths.manifest_file)
    log_dir = paths.slurm_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    repo_root = Path(__file__).resolve().parent.parent.parent
    account_line = f"#SBATCH --account={args.account}\n" if args.account else ""
    script = SBATCH_TEMPLATE.format(
        partition=args.partition,
        qos=args.qos,
        account_line=account_line,
        time_limit=args.time_limit,
        memory=args.memory,
        cpus=args.cpus,
        last_index=len(manifest) - 1,
        max_concurrent=args.max_concurrent,
        log_dir=log_dir,
        repo_root=repo_root,
        data_root=paths.root,
        run_args=args.run_args,
    )
    script_path = paths.slurm_dir / "submit_runs.sbatch"
    script_path.write_text(script)
    print(f"Wrote {script_path} ({len(manifest)} array tasks). Submit with:")
    print(f"  sbatch {script_path}")


def _register_stall_diagnostics() -> None:
    """``kill -USR1 <pid>`` dumps all thread stacks to stderr (stall debugging)."""
    import faulthandler
    import signal

    faulthandler.register(signal.SIGUSR1, all_threads=True)


def main() -> None:
    _register_stall_diagnostics()
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=DEFAULT_DATA_ROOT,
        help=f"Experiment data directory (default: {DEFAULT_DATA_ROOT}).",
    )
    parser.add_argument("--log-level", default="INFO")
    subparsers = parser.add_subparsers(dest="stage", required=True)

    fetch_parser = subparsers.add_parser(
        "fetch-aeronet", help="Download global AERONET V3 AOD measurements."
    )
    fetch_aeronet.add_arguments(fetch_parser)

    matchup_parser = subparsers.add_parser(
        "matchup", help="Search S2 acquisitions and pair with AERONET measurements."
    )
    matchup.add_arguments(matchup_parser)

    manifest_parser = subparsers.add_parser(
        "build-manifest", help="Write the run manifest (matchup x approach tasks)."
    )
    manifest_parser.add_argument(
        "--approaches", nargs="*", choices=sorted(APPROACHES), default=sorted(APPROACHES)
    )
    manifest_parser.add_argument("--limit", type=int, default=None)
    manifest_parser.add_argument(
        "--per-site-per-month",
        type=int,
        default=None,
        help="Seasonally balanced subsample: keep at most N matchups per site per month.",
    )

    run_parser = subparsers.add_parser("run", help="Run SIAC retrievals for matchups x approaches.")
    run_retrieval.add_arguments(run_parser)

    compare_parser = subparsers.add_parser(
        "compare", help="Score retrievals against AERONET and plot."
    )
    compare.add_arguments(compare_parser)

    _add_status_command(subparsers)
    _add_make_slurm_command(subparsers)

    args = parser.parse_args()
    paths = ExperimentPaths(root=args.data_root)
    setup_logging(args.log_level)

    if args.stage == "fetch-aeronet":
        fetch_aeronet.run(args, paths)
    elif args.stage == "matchup":
        matchup.run(args, paths)
    elif args.stage == "build-manifest":
        run_retrieval.build_manifest(
            paths,
            args.approaches,
            limit=args.limit,
            per_site_per_month=args.per_site_per_month,
        )
    elif args.stage == "run":
        run_retrieval.run(args, paths)
    elif args.stage == "compare":
        compare.run(args, paths)
    elif args.stage == "status":
        _status(paths)
    elif args.stage == "make-slurm":
        _make_slurm(args, paths)


if __name__ == "__main__":
    main()
