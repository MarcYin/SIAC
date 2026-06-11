"""Fetch and compile libRadtran (``uvspec``) from source for the RT backend.

libRadtran is an autotools C/Fortran package; its ``uvspec`` engine is built
in-tree (``<src>/bin/uvspec``) and reads its physics tables from ``<src>/data``.
The core tarball bundles the US-Standard atmosphere and the *coarse* ``reptran``
band model, but two extras must be merged in for the remote-LUT preset: the OPAC
per-species optical properties the ``continental_average`` aerosol needs (the
``data/aerosol`` subset of the ``optprop`` archive), and — when a non-coarse band
model is configured (the production LUT used ``reptran fine``) — the matching
``reptran_solar_<resolution>.lookup.*`` tables (the ``reptran`` archive). Only
the needed subsets are extracted to keep the merge small.

This mirrors :mod:`siac.algorithms.rt.direct.sixs_build`: atomic download with
SHA-256 verification, a path-traversal-guarded unpack, an on-disk cache under
``~/.cache/siac/libradtran/<profile>/``, and an :func:`ensure_libradtran` that
returns the prebuilt engine or builds on demand. The build recipe (Apple-Silicon
triplet + conda NetCDF prefix) was validated against libRadtran 2.0.6 under the
Pixi ``libradtran`` environment.
"""

from __future__ import annotations

import contextlib
import logging
import os
import platform
import shutil
import subprocess
import sys
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from siac.algorithms.rt.direct._build_common import (
    archive_sha256,
    conda_prefix,
    fetch_archive,
    prepend_conda_lib_path,
)

if TYPE_CHECKING:
    from collections.abc import Iterator

    from siac.config.algorithms import LibRadtranAlgorithmConfig

logger = logging.getLogger(__name__)

#: Marker file identifying an unpacked libRadtran source tree.
_CONFIGURE_MARKER = "configure"
#: Subdirectory of the cache root holding the unpacked + built source tree.
_SRC_DIRNAME = "src"
#: SHA-256 of the upstream archives this harness was validated against
#: (libRadtran 2.0.6 + optprop_v2.1). Env vars override per-run; ``None`` skips.
_SOURCE_SHA256: str | None = "64930cc40b6e4a37aa220520974d330fc1563796f466a649b2238131f2d69840"
_OPTPROP_SHA256: str | None = "11daa1f1f4be0fd4ddf7e881ec2005498049674a1540d37b4b1e8f5e16052c7e"
_REPTRAN_SHA256: str | None = "55893c80bcc999651bac3bf014ee64aaf602653ba640eb5bebe787a5d8eacce7"
#: Generous safety caps (the pinned SHA-256s already gate the default archives;
#: these protect the hash-unset / env-override path from a decompression bomb or
#: a runaway download). Far above any legitimate libRadtran archive/subset.
_MAX_DOWNLOAD_BYTES = 4 * 1024**3
_MAX_EXTRACT_BYTES = 12 * 1024**3


@dataclass(frozen=True)
class LibRadtranBuildPaths:
    """Resolved paths used while fetching + building libRadtran."""

    root_dir: Path
    source_archive: Path
    optprop_archive: Path
    reptran_archive: Path
    extract_dir: Path


@dataclass(frozen=True)
class LibRadtranPaths:
    """Resolved runtime paths to a usable libRadtran engine."""

    uvspec: Path
    data_dir: Path


def default_build_root() -> Path:
    return Path.home().expanduser() / ".cache" / "siac" / "libradtran"


def _archive_name_from_url(url: str) -> str:
    """Derive a SAFE cache filename from a libRadtran download URL.

    Handles both the plain ``.../libRadtran-2.0.6.tar.gz`` form and the
    DokuWiki ``fetch.php?media=download:optprop_v2.1.tar.gz`` form. The result is
    reduced to a bare basename so a crafted URL (e.g. ``media=...:../../etc/x``)
    cannot escape the cache directory (path traversal / write-what-where).
    """
    if "media=" in url:
        token = url.split("media=", 1)[1].split("&", 1)[0]
        raw = token.split(":")[-1]
    else:
        raw = url.rsplit("/", 1)[-1].split("?", 1)[0]
    name = Path(raw.strip().replace("\\", "/").rstrip("/")).name
    if not name or name in (".", ".."):
        raise ValueError(f"Cannot derive a safe archive filename from URL: {url!r}")
    return name


def resolve_build_paths(config: LibRadtranAlgorithmConfig) -> LibRadtranBuildPaths:
    default_root = default_build_root() / str(getattr(config, "build_profile", "release"))
    root_dir = Path(config.build_dir or default_root).expanduser()
    return LibRadtranBuildPaths(
        root_dir=root_dir,
        source_archive=root_dir / _archive_name_from_url(config.source_url),
        optprop_archive=root_dir / _archive_name_from_url(config.optprop_url),
        reptran_archive=root_dir / _archive_name_from_url(config.reptran_url),
        extract_dir=root_dir / _SRC_DIRNAME,
    )


def _reptran_resolution(mol_abs_param: str) -> str | None:
    """Return the reptran resolution token that needs auxiliary data, else None.

    The core tarball bundles only the ``coarse`` reptran lookup tables, so
    ``fine``/``medium`` require the separate ``reptran`` archive. Plain
    ``reptran`` (and non-reptran band models) need nothing extra.
    """
    tokens = mol_abs_param.lower().split()
    if not tokens or tokens[0] != "reptran":
        return None
    resolution = tokens[1] if len(tokens) > 1 else "coarse"
    return resolution if resolution in ("fine", "medium") else None


def _required_reptran_resolutions(config: LibRadtranAlgorithmConfig) -> list[str]:
    """All non-coarse reptran resolutions the runtime config will request.

    The runner uses ``mol_abs_param`` as the base plus, per spectral region,
    either explicit ``mol_abs_regions`` models or - when those are unset and
    ``adaptive_deep_water_fine`` is on - ``reptran fine`` in the deep H2O bands.
    The build must merge lookup tables for ALL of these, not only the base
    (otherwise a fresh build of the adaptive default would lack the ``fine``
    tables its water-band segments need and fail at runtime).
    """
    resolutions: set[str] = set()
    base = _reptran_resolution(str(config.mol_abs_param))
    if base is not None:
        resolutions.add(base)
    regions = getattr(config, "mol_abs_regions", None)
    if regions:
        for _lo, _hi, model in regions:
            res = _reptran_resolution(str(model))
            if res is not None:
                resolutions.add(res)
    elif getattr(config, "adaptive_deep_water_fine", False):
        # Adaptive default upgrades the deep H2O bands to fine - but only require
        # the fine tables if a band actually overlaps the configured window (a
        # local import avoids a module-load cycle with the runner).
        from siac.algorithms.rt.direct.libradtran_runner import DEEP_WATER_H2O_BANDS_NM

        wmin, wmax = float(config.wavelength_min_nm), float(config.wavelength_max_nm)
        if any(min(hi, wmax) > max(lo, wmin) for lo, hi in DEEP_WATER_H2O_BANDS_NM):
            resolutions.add("fine")
    return sorted(resolutions)


def _archive_sha256(path: Path) -> str:
    return archive_sha256(path)


def _fetch_archive(
    url: str, archive_path: Path, *, expected_sha256: str | None, env_var: str
) -> None:
    """Download ``url`` to ``archive_path`` atomically and verify its SHA-256."""
    fetch_archive(
        url,
        archive_path,
        expected_sha256=expected_sha256,
        sha_env_var=env_var,
        max_bytes=_MAX_DOWNLOAD_BYTES,
        what="libRadtran",
    )


def _safe_members(
    archive: tarfile.TarFile,
    dest: Path,
    *,
    prefix: str | None = None,
    name_prefix: str | None = None,
) -> list[tarfile.TarInfo]:
    """Return archive members that stay under ``dest``.

    ``prefix`` restricts to a directory subtree; ``name_prefix`` restricts to
    members whose path begins with a literal string (for file-prefix subsets
    like ``.../reptran_solar_fine``).
    """
    root = dest.resolve()
    members: list[tarfile.TarInfo] = []
    for member in archive.getmembers():
        if prefix is not None and not (
            member.name == prefix or member.name.startswith(prefix + "/")
        ):
            continue
        if name_prefix is not None and not member.name.startswith(name_prefix):
            continue
        target = (root / member.name).resolve()
        if target != root and root not in target.parents:
            raise RuntimeError(
                f"Refusing to extract unsafe libRadtran archive member: {member.name!r}"
            )
        if member.islnk() or member.issym():
            logger.debug("Skipping linked archive member: %s", member.name)
            continue
        members.append(member)
    total = sum(max(0, int(m.size)) for m in members)
    if total > _MAX_EXTRACT_BYTES:
        raise RuntimeError(
            f"Refusing to extract libRadtran archive: selected members total "
            f"{total // 1024**2} MB exceed the {_MAX_EXTRACT_BYTES // 1024**3} GB cap "
            "(possible decompression bomb)."
        )
    return members


def _safe_extractall(
    archive: tarfile.TarFile, *, path: Path, members: list[tarfile.TarInfo]
) -> None:
    """``extractall`` with the tar ``data`` filter where available.

    Defense-in-depth on top of the ``_safe_members`` traversal/symlink guard: the
    ``data`` filter (Python 3.12+, backported to recent 3.9-3.11 patch releases)
    independently rejects absolute paths, ``..`` traversal, and escaping links.
    """
    if hasattr(tarfile, "data_filter"):
        archive.extractall(path=path, members=members, filter="data")  # noqa: S202
    else:  # pragma: no cover - only on old patch releases lacking the filter
        archive.extractall(path=path, members=members)  # noqa: S202


def _locate_source_tree(extract_dir: Path) -> Path:
    """Return the directory containing ``configure`` under ``extract_dir``."""
    if (extract_dir / _CONFIGURE_MARKER).exists():
        return extract_dir
    candidates = [
        path
        for path in extract_dir.iterdir()
        if path.is_dir() and (path / _CONFIGURE_MARKER).exists()
    ]
    if len(candidates) == 1:
        return candidates[0]
    raise RuntimeError(
        f"Unable to locate an unpacked libRadtran source tree (with a {_CONFIGURE_MARKER!r}) "
        f"under {extract_dir}."
    )


def _unpack_source(paths: LibRadtranBuildPaths) -> Path:
    """Unpack the core source tarball and return the source-tree directory."""
    if paths.extract_dir.exists():
        try:
            return _locate_source_tree(paths.extract_dir)
        except RuntimeError:
            shutil.rmtree(paths.extract_dir)
    paths.extract_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(paths.source_archive) as archive:
        _safe_extractall(
            archive, path=paths.extract_dir, members=_safe_members(archive, paths.extract_dir)
        )
    return _locate_source_tree(paths.extract_dir)


def _merge_optprop_aerosol(archive_path: Path, source_dir: Path) -> None:
    """Merge only the ``data/aerosol`` OPAC subset of optprop into the source tree.

    The optprop archive is dominated by cloud (wc/ic) optical properties that are
    irrelevant to clear-sky atmospheric correction; only the OPAC aerosol
    optical properties (~150 MB vs ~5 GB unpacked) are needed for the
    ``continental_average`` species, so just that subtree is extracted.
    """
    if (source_dir / "data" / "aerosol" / "OPAC" / "optprop").is_dir() and any(
        (source_dir / "data" / "aerosol" / "OPAC" / "optprop").glob("*.cdf")
    ):
        return  # already merged
    with tarfile.open(archive_path) as archive:
        members = _safe_members(archive, source_dir, prefix="data/aerosol")
        if not members:
            raise RuntimeError(
                f"optprop archive {archive_path} contains no 'data/aerosol' members; "
                "cannot provide OPAC optical properties for continental_average aerosol."
            )
        _safe_extractall(archive, path=source_dir, members=members)


def _merge_reptran(archive_path: Path, source_dir: Path, resolution: str) -> None:
    """Merge the ``reptran_solar_<resolution>`` lookup tables into the source tree.

    The core tarball bundles only the ``coarse`` reptran lookups; ``fine`` and
    ``medium`` (the production setting) need the separate reptran archive. Only
    the solar lookups for the requested resolution are extracted.
    """
    name_prefix = f"data/correlated_k/reptran/reptran_solar_{resolution}"
    lookups = sorted(
        (source_dir / "data" / "correlated_k" / "reptran").glob(
            f"reptran_solar_{resolution}.lookup.*.cdf"
        )
    )
    if lookups:
        return  # already merged
    with tarfile.open(archive_path) as archive:
        members = _safe_members(archive, source_dir, name_prefix=name_prefix)
        if not members:
            raise RuntimeError(
                f"reptran archive {archive_path} contains no {name_prefix!r} members; "
                f"cannot provide the 'reptran {resolution}' band model."
            )
        _safe_extractall(archive, path=source_dir, members=members)


def _conda_prefix() -> str:
    return conda_prefix()


def _configure_args() -> list[str]:
    """Configure flags: conda NetCDF prefix + Apple-Silicon triplet workaround."""
    prefix = _conda_prefix()
    args = [f"--with-netcdf4={prefix}"]
    if sys.platform == "darwin" and platform.machine() in ("arm64", "aarch64"):
        # libRadtran 2.0.6 ships a 2013 config.sub that rejects conda's
        # ``arm64-apple-darwin`` triplet but accepts ``aarch64-apple-darwin``.
        args = ["--build=aarch64-apple-darwin", "--host=aarch64-apple-darwin", *args]
    return args


def _build_environment() -> dict[str, str]:
    """Environment for configure/make surfacing the conda prefix to autotools."""
    env = dict(os.environ)
    prefix = _conda_prefix()
    inc = str(Path(prefix) / "include")
    lib = str(Path(prefix) / "lib")
    env["CPPFLAGS"] = f"-I{inc} " + env.get("CPPFLAGS", "")
    env["LDFLAGS"] = f"-L{lib} " + env.get("LDFLAGS", "")
    return prepend_conda_lib_path(env)


def _run(cmd: list[str], *, cwd: Path, env: dict[str, str], step: str) -> None:
    logger.info("libRadtran build [%s]: %s (cwd=%s)", step, " ".join(cmd), cwd)
    result = subprocess.run(  # noqa: S603 - fixed argv, no shell
        cmd, cwd=str(cwd), env=env, capture_output=True, text=True, check=False
    )
    if result.returncode != 0:
        tail = (result.stdout[-2000:] + "\n" + result.stderr[-2000:]).strip()
        raise RuntimeError(
            f"libRadtran build step {step!r} failed (exit {result.returncode}).\n{tail}"
        )


def find_uvspec(source_dir: Path) -> Path | None:
    """Return ``<src>/bin/uvspec`` if it exists."""
    candidate = source_dir / "bin" / "uvspec"
    return candidate if candidate.exists() else None


def _toolchain_available() -> bool:
    cc = os.environ.get("CC") or "cc"
    has_cc = any(shutil.which(c) for c in (cc, "clang", "gcc") if c)
    return shutil.which("make") is not None and bool(has_cc)


@contextlib.contextmanager
def _build_lock(root_dir: Path) -> Iterator[None]:
    """Best-effort exclusive lock so concurrent builds don't clobber the tree.

    Two processes building the same profile (e.g. an auto-build racing a manual
    ``build-libradtran``) would otherwise ``rmtree``/extract into the same source
    tree at once. Falls back to a no-op where ``fcntl`` is unavailable.
    """
    try:
        import fcntl
    except ImportError:  # pragma: no cover - non-POSIX
        yield
        return
    with (root_dir / ".build.lock").open("w") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def build_libradtran(config: LibRadtranAlgorithmConfig) -> LibRadtranPaths:
    """Fetch, merge OPAC aerosol (+ reptran) data, and compile; return engine paths."""
    paths = resolve_build_paths(config)
    paths.root_dir.mkdir(parents=True, exist_ok=True)
    with contextlib.suppress(OSError):
        paths.root_dir.chmod(0o700)  # keep the built engine off other users on shared hosts
    with _build_lock(paths.root_dir):
        return _build_libradtran_locked(config, paths)


def _build_libradtran_locked(
    config: LibRadtranAlgorithmConfig, paths: LibRadtranBuildPaths
) -> LibRadtranPaths:
    """Body of :func:`build_libradtran`, run while holding the cross-process lock."""
    if not _toolchain_available():
        raise RuntimeError(
            "libRadtran build requires a C/Fortran toolchain and `make`. Use the Pixi "
            "`libradtran` environment (`pixi run -e libradtran build-libradtran`)."
        )

    _fetch_archive(
        config.source_url,
        paths.source_archive,
        expected_sha256=_SOURCE_SHA256,
        env_var="SIAC_LIBRADTRAN_SOURCE_SHA256",
    )
    _fetch_archive(
        config.optprop_url,
        paths.optprop_archive,
        expected_sha256=_OPTPROP_SHA256,
        env_var="SIAC_LIBRADTRAN_OPTPROP_SHA256",
    )
    # reptran 'medium'/'fine' lookup tables are NOT bundled (only 'coarse' is).
    # Fetch + merge EVERY non-coarse resolution the config will request - the
    # base mol_abs_param AND the mol_abs_regions / adaptive deep-water-fine
    # segments - so a fresh adaptive build has the 'fine' tables its water bands
    # need. (medium and fine both come from the same reptran archive.)
    reptran_resolutions = _required_reptran_resolutions(config)
    if reptran_resolutions:
        _fetch_archive(
            config.reptran_url,
            paths.reptran_archive,
            expected_sha256=_REPTRAN_SHA256,
            env_var="SIAC_LIBRADTRAN_REPTRAN_SHA256",
        )

    source_dir = _unpack_source(paths)
    _merge_optprop_aerosol(paths.optprop_archive, source_dir)
    for resolution in reptran_resolutions:
        _merge_reptran(paths.reptran_archive, source_dir, resolution)

    uvspec = find_uvspec(source_dir)
    if uvspec is None:
        env = _build_environment()
        _run(["./configure", *_configure_args()], cwd=source_dir, env=env, step="configure")
        _run(["make"], cwd=source_dir, env=env, step="make")
        uvspec = find_uvspec(source_dir)
    if uvspec is None:
        raise RuntimeError(
            f"libRadtran build completed but {source_dir / 'bin' / 'uvspec'} was not produced."
        )
    data_dir = source_dir / "data"
    if not data_dir.is_dir():
        raise RuntimeError(f"libRadtran data directory not found at {data_dir}.")
    return LibRadtranPaths(uvspec=uvspec, data_dir=data_dir)


def _ensure_reptran_data(config: LibRadtranAlgorithmConfig, source_dir: Path) -> None:
    """Ensure every requested reptran resolution's lookups exist in a cached tree.

    Guards the case where ``uvspec`` was built before a resolution was requested
    (the binary exists, but its ``reptran_solar_<res>.lookup.*`` tables are
    absent - e.g. the base was upgraded, or mol_abs_regions / the adaptive
    deep-water-fine default added ``fine``). Only fetches when ``auto_build`` is
    enabled.
    """
    reptran_dir = source_dir / "data" / "correlated_k" / "reptran"
    missing = [
        res
        for res in _required_reptran_resolutions(config)
        if not list(reptran_dir.glob(f"reptran_solar_{res}.lookup.*.cdf"))
    ]
    if not missing:
        return
    if not config.auto_build:
        raise RuntimeError(
            f"libRadtran reptran lookup tables {missing} are missing under {reptran_dir} "
            "and auto_build is disabled. Run `pixi run -e libradtran build-libradtran`."
        )
    paths = resolve_build_paths(config)
    _fetch_archive(
        config.reptran_url,
        paths.reptran_archive,
        expected_sha256=_REPTRAN_SHA256,
        env_var="SIAC_LIBRADTRAN_REPTRAN_SHA256",
    )
    for resolution in missing:
        _merge_reptran(paths.reptran_archive, source_dir, resolution)


def ensure_libradtran(config: LibRadtranAlgorithmConfig) -> LibRadtranPaths:
    """Return a usable libRadtran engine, building it on demand if allowed."""
    # 1. Explicit prebuilt binary.
    if config.uvspec_path is not None:
        uvspec = Path(config.uvspec_path).expanduser()
        if not uvspec.exists():
            raise RuntimeError(f"Configured libradtran.uvspec_path does not exist: {uvspec}")
        data_dir = (
            Path(config.data_dir).expanduser()
            if config.data_dir is not None
            else uvspec.parent.parent / "data"
        )
        return LibRadtranPaths(uvspec=uvspec, data_dir=data_dir)

    # 2. Pre-unpacked source tree.
    if config.source_dir is not None:
        source_dir = Path(config.source_dir).expanduser()
        candidate = find_uvspec(source_dir)
        if candidate is not None:
            _ensure_reptran_data(config, source_dir)
            return LibRadtranPaths(uvspec=candidate, data_dir=source_dir / "data")

    # 3. Cached build.
    paths = resolve_build_paths(config)
    if paths.extract_dir.exists():
        try:
            cached = _locate_source_tree(paths.extract_dir)
        except RuntimeError:
            cached = None
        if cached is not None:
            candidate = find_uvspec(cached)
            if candidate is not None:
                _ensure_reptran_data(config, cached)
                return LibRadtranPaths(uvspec=candidate, data_dir=cached / "data")

    if not config.auto_build:
        raise RuntimeError(
            "libRadtran (uvspec) is not available and auto_build is disabled. Expected a built "
            f"engine under {paths.root_dir}. Run `pixi run -e libradtran build-libradtran`."
        )
    return build_libradtran(config)
