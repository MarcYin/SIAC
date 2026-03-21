"""Command-line interface for SIAC."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

from siac import __version__
from siac.config import SIACConfig
from siac.errors import SIACError


def _non_empty_text(value: str) -> str:
    text = value.strip()
    if not text:
        raise argparse.ArgumentTypeError("query must not be empty")
    return text


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="siac",
        description="SIAC command-line interface.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"siac {__version__}",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)
    process = subparsers.add_parser(
        "process-s2",
        help="Process a Sentinel-2 SAFE, product ID, or tile/date query.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    process.add_argument(
        "query",
        type=_non_empty_text,
        help="Local SAFE path, product ID, or tile/date shorthand.",
    )
    process.add_argument(
        "--config",
        type=Path,
        help="Path to a TOML config file. If omitted, built-in defaults are used.",
    )
    process.add_argument(
        "--output-path",
        type=Path,
        help="Directory where corrected products will be written.",
    )
    aoi_group = process.add_mutually_exclusive_group()
    aoi_group.add_argument(
        "--aoi",
        type=Path,
        help="Path to an AOI GeoJSON file or a raster-backed AOI source.",
    )
    aoi_group.add_argument(
        "--aoi-bbox",
        nargs=4,
        type=float,
        metavar=("MINX", "MINY", "MAXX", "MAXY"),
        help="AOI bounds in the order minx miny maxx maxy.",
    )
    process.set_defaults(handler=_run_process_s2)
    return parser


def _load_config(config_path: Path | None) -> SIACConfig:
    if config_path is None:
        try:
            return SIACConfig.load()
        except FileNotFoundError:
            return SIACConfig()
    return SIACConfig.from_file(config_path)


def _process_kwargs(args: argparse.Namespace) -> dict[str, Any]:
    kwargs: dict[str, Any] = {}
    if args.output_path is not None:
        kwargs["output_path"] = args.output_path
    if args.aoi_bbox is not None:
        kwargs["aoi"] = tuple(float(value) for value in args.aoi_bbox)
    elif args.aoi is not None:
        kwargs["aoi"] = args.aoi
    return kwargs


def _run_process_s2(args: argparse.Namespace) -> int:
    from siac.api.public import siac_process_s2

    config = _load_config(args.config)
    result = siac_process_s2(config, args.query, **_process_kwargs(args))
    print(f"Sentinel-2 processing complete. Mean AOT: {float(result.aot.mean()):.3f}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:
        code = exc.code
        return code if isinstance(code, int) else 1

    handler = getattr(args, "handler", None)
    if handler is None:
        parser.error("missing command")

    try:
        return int(handler(args))
    except KeyboardInterrupt:
        return 130
    except (FileNotFoundError, SIACError, ValueError) as exc:
        print(f"siac: error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
