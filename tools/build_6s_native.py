"""Build or clean the patched native 6SV2.1 Python extension."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from siac.algorithms.rt.direct.sixs_build import (
    build_native_sixs_module,
    resolve_build_paths,
)
from siac.config.algorithms import SixSAlgorithmConfig


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-url", default=None, help="Override the upstream 6S source tarball URL."
    )
    parser.add_argument(
        "--source-dir", type=Path, default=None, help="Use an existing unpacked 6S source tree."
    )
    parser.add_argument(
        "--build-dir", type=Path, default=None, help="Override the native build root."
    )
    parser.add_argument(
        "--module-path", type=Path, default=None, help="Override the output Python extension path."
    )
    parser.add_argument(
        "--library-path", type=Path, default=None, help="Deprecated alias for --module-path."
    )
    parser.add_argument(
        "--compiler", default="gfortran", help="Fortran compiler executable or path."
    )
    parser.add_argument(
        "--build-profile",
        choices=("release", "parity"),
        default="release",
        help="Compiler profile to use for the native 6S module.",
    )
    parser.add_argument(
        "--clean", action="store_true", help="Remove the resolved build root instead of compiling."
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    config = SixSAlgorithmConfig(
        source_url=args.source_url or SixSAlgorithmConfig().source_url,
        source_dir=args.source_dir,
        build_dir=args.build_dir,
        module_path=args.module_path or args.library_path,
        compiler=args.compiler,
        build_profile=args.build_profile,
    )
    paths = resolve_build_paths(config)

    if args.clean:
        if paths.root_dir.exists():
            shutil.rmtree(paths.root_dir)
            print(f"Removed {paths.root_dir}")
        else:
            print(f"Nothing to clean at {paths.root_dir}")
        return 0

    module_path = build_native_sixs_module(config)
    print(module_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
