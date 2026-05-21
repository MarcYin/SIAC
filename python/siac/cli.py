"""Command-line interface for SIAC."""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import Any, cast

from siac import __version__
from siac.config import SIACConfig
from siac.domain.aoi import AOI
from siac.errors import ConfigurationError, SIACError


def _non_empty_text(value: str) -> str:
    text = value.strip()
    if not text:
        raise argparse.ArgumentTypeError("query must not be empty")
    return text


def _month_number(value: str) -> int:
    month = int(value)
    if not 1 <= month <= 12:
        raise argparse.ArgumentTypeError("month must be between 1 and 12")
    return month


def _year_month(value: str) -> tuple[int, int]:
    text = value.strip()
    match = re.fullmatch(r"(\d{4})-(\d{2})", text)
    if match is None:
        raise argparse.ArgumentTypeError("year-month must be in YYYY-MM format")
    year = int(match.group(1))
    month = int(match.group(2))
    if not 1 <= month <= 12:
        raise argparse.ArgumentTypeError("month must be between 1 and 12")
    return year, month


def _add_aoi_arguments(parser: argparse.ArgumentParser, *, required: bool = False) -> None:
    aoi_group = parser.add_mutually_exclusive_group(required=required)
    aoi_group.add_argument(
        "--aoi-file",
        type=Path,
        help="Path to an AOI GeoJSON file.",
    )
    aoi_group.add_argument(
        "--aoi-wkt",
        type=_non_empty_text,
        help="AOI geometry as a WKT string.",
    )
    aoi_group.add_argument(
        "--aoi-bbox",
        nargs=4,
        type=float,
        metavar=("MINX", "MINY", "MAXX", "MAXY"),
        help="AOI bounds in the order minx miny maxx maxy.",
    )
    parser.add_argument(
        "--aoi-crs",
        default=None,
        help="CRS for AOI coordinates. If omitted, SIAC assumes WGS84 and validates degree-like coordinates.",
    )


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
    _add_aoi_arguments(process)
    process.set_defaults(handler=_run_process_s2)

    prepare = subparsers.add_parser(
        "prepare-monthly-composites",
        help="Build prepared monthly composites for a given AOI and period selection.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    prepare.add_argument(
        "--config",
        type=Path,
        help="Path to a TOML config file. If omitted, built-in defaults are used.",
    )
    prepare.add_argument(
        "--output-path",
        type=Path,
        required=True,
        help="Directory where the prepared monthly composite store will be written.",
    )
    prepare.add_argument(
        "--resolution",
        type=float,
        default=None,
        help="Target composite grid resolution. If omitted, SIAC uses the source dataset resolution.",
    )
    prepare.add_argument(
        "--year-month",
        action="append",
        type=_year_month,
        default=None,
        help="Explicit year-month selection in YYYY-MM format. Repeat to request multiple periods.",
    )
    prepare.add_argument(
        "--year",
        action="append",
        type=int,
        default=None,
        help="Year selection for cross-product period generation. Repeat to request multiple years.",
    )
    prepare.add_argument(
        "--month",
        action="append",
        type=_month_number,
        default=None,
        help="Month selection for cross-product period generation. Repeat to request multiple months.",
    )
    _add_aoi_arguments(prepare, required=True)
    prepare.set_defaults(handler=_run_prepare_monthly_composites)
    return parser


def _load_config(config_path: Path | None) -> SIACConfig:
    if config_path is None:
        try:
            return SIACConfig.load()
        except FileNotFoundError:
            return SIACConfig()
    return SIACConfig.from_file(config_path)


def _configure_logging(level_name: str) -> None:
    level = getattr(logging, str(level_name).upper(), logging.INFO)
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s:%(name)s:%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )
    # Attach the secret-redaction filter to every root handler so credentials
    # are scrubbed regardless of which submodule emitted the record. Until
    # now the filter only shadowed ``siac.adapters.auth``, leaving URLs that
    # carry tokens un-redacted in cams.py / mcd43_earthaccess.py /
    # copernicus_dataspace.py logs (REVIEW.md §2.7). A filter set on the
    # ``Logger`` only fires for records the logger emits directly, so we
    # install it on every existing handler instead.
    from siac.adapters._log_filter import SecretRedactionFilter

    for handler in logging.getLogger().handlers:
        if not any(isinstance(f, SecretRedactionFilter) for f in handler.filters):
            handler.addFilter(SecretRedactionFilter())


def _resolve_cli_aoi(args: argparse.Namespace) -> Any | None:
    if args.aoi_bbox is not None:
        return AOI.from_bounds(
            tuple(float(value) for value in args.aoi_bbox),
            crs=cast("str | None", getattr(args, "aoi_crs", None)),
        )
    if getattr(args, "aoi_file", None) is not None:
        return AOI.from_geojson(
            args.aoi_file, crs=cast("str | None", getattr(args, "aoi_crs", None))
        )
    if getattr(args, "aoi_wkt", None) is not None:
        # Parse WKT into a GeoJSON dict via siac.geo.geometry, then construct AOI.
        # Previously this passed WKT directly to AOI.from_geojson which silently
        # corrupted the AOI (REVIEW.md §1.1 #2).
        from siac.geo.geometry import _parse_wkt

        try:
            geojson = _parse_wkt(args.aoi_wkt, _target_crs=None)
        except (ValueError, AttributeError) as exc:
            raise ConfigurationError(f"Failed to parse --aoi-wkt as WKT: {exc}") from exc
        return AOI.from_geojson(geojson, crs=cast("str | None", getattr(args, "aoi_crs", None)))
    return None


def _process_kwargs(args: argparse.Namespace) -> dict[str, Any]:
    kwargs: dict[str, Any] = {}
    if args.output_path is not None:
        kwargs["output_path"] = args.output_path
    aoi = _resolve_cli_aoi(args)
    if aoi is not None:
        kwargs["aoi"] = aoi
    return kwargs


def _selected_year_months(args: argparse.Namespace) -> tuple[tuple[int, int], ...]:
    explicit = tuple(args.year_month or ())
    if explicit:
        if args.year or args.month:
            raise ValueError("Use either --year-month or --year/--month, not both")
        pairs = explicit
    else:
        years = tuple(int(year) for year in (args.year or ()))
        months = tuple(int(month) for month in (args.month or ()))
        if not years or not months:
            raise ValueError("Provide either --year-month or both --year and --month")
        pairs = tuple((year, month) for year in years for month in months)

    if len(set(pairs)) != len(pairs):
        raise ValueError("Duplicate year/month selections are not allowed")
    return tuple(sorted(pairs))


def _run_process_s2(args: argparse.Namespace) -> int:
    from siac.api.public import siac_process_s2

    config = _load_config(args.config)
    _configure_logging(config.runtime.log_level)
    result = siac_process_s2(config, args.query, **_process_kwargs(args))
    print(f"Sentinel-2 processing complete. Mean AOT: {float(result.aot.mean()):.3f}")
    return 0


def _run_prepare_monthly_composites(args: argparse.Namespace) -> int:
    from siac.api.public import prepare_monthly_composites

    config = _load_config(args.config)
    _configure_logging(config.runtime.log_level)
    aoi = _resolve_cli_aoi(args)
    if not isinstance(aoi, AOI):
        if aoi is None:
            raise ValueError("An AOI is required")
        aoi = AOI.from_geojson(aoi) if isinstance(aoi, Path) else AOI.from_geojson(str(aoi))
    result = prepare_monthly_composites(
        config,
        aoi=aoi,
        year_months=_selected_year_months(args),
        resolution=(float(args.resolution) if args.resolution is not None else None),
        output_path=args.output_path,
    )
    periods = ", ".join(result.periods)
    print(
        "Prepared monthly composites written to "
        f"{result.store_path} ({result.period_count} period(s): {periods}; representation={result.representation})."
    )
    print(
        "Use with SIAC via "
        "providers.monthly_composites.kind='prepared_store' "
        f"and providers.monthly_composites.store_path='{result.store_path}'."
    )
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
