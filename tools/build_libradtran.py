"""Build or clean the libRadtran ``uvspec`` engine from source.

Fetches + verifies + compiles libRadtran and merges the reptran/optprop data the
configured band model needs. With no arguments it builds for the default config
(``reptran medium`` base + adaptive ``reptran fine`` deep-water bands), so it
merges both the ``medium`` and ``fine`` lookup tables from the reptran archive.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from siac.algorithms.rt.direct.libradtran_build import (
    build_libradtran,
    resolve_build_paths,
)
from siac.config.algorithms import LibRadtranAlgorithmConfig


def _parse_args() -> argparse.Namespace:
    defaults = LibRadtranAlgorithmConfig()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-url", default=None, help="Override the libRadtran source archive URL."
    )
    parser.add_argument(
        "--reptran-url", default=None, help="Override the reptran data archive URL."
    )
    parser.add_argument(
        "--optprop-url", default=None, help="Override the optprop data archive URL."
    )
    parser.add_argument(
        "--build-dir", type=Path, default=None, help="Override the build/cache root."
    )
    parser.add_argument(
        "--build-profile",
        choices=("release", "parity"),
        default="release",
        help="Build profile (cache subdirectory).",
    )
    parser.add_argument(
        "--mol-abs-param",
        default=defaults.mol_abs_param,
        help="Base band model; selects which reptran tables are merged (default: %(default)r).",
    )
    parser.add_argument(
        "--no-adaptive-deep-water-fine",
        action="store_true",
        help="Do not also merge 'reptran fine' for the adaptive deep-water bands.",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove the resolved build root instead of compiling.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    base = LibRadtranAlgorithmConfig()
    config = LibRadtranAlgorithmConfig(
        source_url=args.source_url or base.source_url,
        reptran_url=args.reptran_url or base.reptran_url,
        optprop_url=args.optprop_url or base.optprop_url,
        build_dir=args.build_dir,
        build_profile=args.build_profile,
        mol_abs_param=args.mol_abs_param,
        adaptive_deep_water_fine=not args.no_adaptive_deep_water_fine,
    )
    paths = resolve_build_paths(config)

    if args.clean:
        if paths.root_dir.exists():
            shutil.rmtree(paths.root_dir)
            print(f"Removed {paths.root_dir}")
        else:
            print(f"Nothing to clean at {paths.root_dir}")
        return 0

    result = build_libradtran(config)
    print(result.uvspec)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
