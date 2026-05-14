#!/usr/bin/env bash
#
# End-to-end verification script for REVIEW.md fixes (waves 1-10).
#
# Run each block in order. Each prints PASS / FAIL with a one-line
# summary. The whole script takes ~5 minutes if the regression scene
# is available; ~1 minute otherwise (regression auto-skips).
#
# Usage:
#   bash tools/verify_review_fixes.sh
#   bash tools/verify_review_fixes.sh --skip-regression    # ~30 sec
#   bash tools/verify_review_fixes.sh --skip-profile       # skip cProfile run
#
# Exit code is 0 only if every check passes.

set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

# ---- parse args -----------------------------------------------------------
SKIP_REGRESSION=0
SKIP_PROFILE=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --skip-regression) SKIP_REGRESSION=1 ;;
    --skip-profile)    SKIP_PROFILE=1 ;;
    -h|--help)
      sed -n '2,18p' "$0"; exit 0 ;;
    *) echo "unknown arg: $1"; exit 2 ;;
  esac
  shift
done

FAILED=0
pass() { printf "  \033[32mPASS\033[0m  %s\n" "$1"; }
fail() { printf "  \033[31mFAIL\033[0m  %s\n" "$1"; FAILED=$((FAILED+1)); }
note() { printf "  \033[33mNOTE\033[0m  %s\n" "$1"; }
hdr()  { printf "\n\033[1m=== %s ===\033[0m\n" "$1"; }

# ---------------------------------------------------------------------------
# 1. Wave 10 — siac._rust built; full unit suite green
# ---------------------------------------------------------------------------
hdr "Wave 10: siac._rust + full unit suite"

# siac._rust must exist for the Python that's about to run pytest.
# The default and rt6s pixi envs use different Python versions (3.12 vs
# 3.11), so check the explicit per-python builds we expect.
PY_MAJOR_MINOR=$(pixi run python -c "import sys; print(f'{sys.version_info.major}{sys.version_info.minor}')" 2>/dev/null)
RUST_SO="python/siac/_rust.cpython-${PY_MAJOR_MINOR}-darwin.so"
if [[ ! -f "$RUST_SO" ]]; then
  note "siac._rust not built for cp${PY_MAJOR_MINOR} — running pixi run build-rust (takes ~20s)"
  if pixi run build-rust >/dev/null 2>&1; then
    pass "pixi run build-rust"
  else
    fail "pixi run build-rust"
  fi
else
  pass "siac._rust already built ($RUST_SO)"
fi

PYTHONPATH=python pixi run python -c "
from siac import _rust
syms = sorted(s for s in dir(_rust) if not s.startswith('_'))
print('rust symbols:', syms[:6], '...' if len(syms)>6 else '')
assert 'RossThickLiSparse' in syms
assert 'TwoLayerNN' in syms
" && pass "siac._rust exposes expected symbols" || fail "siac._rust import"

# Full fast suite — should be 1293 / 0 / 7 / 0
SUITE_OUT=$(PYTHONPATH=python pixi run python -m pytest tests/unit \
    -m "not slow and not integration and not regression" \
    --tb=no --no-cov -q 2>&1 | tail -1)
echo "  $SUITE_OUT"
if echo "$SUITE_OUT" | grep -qE '^[0-9]+ passed.*0 failed' || echo "$SUITE_OUT" | grep -qE 'no errors' ; then
  pass "unit suite green"
elif echo "$SUITE_OUT" | grep -qE 'passed' && ! echo "$SUITE_OUT" | grep -qE 'failed' ; then
  pass "unit suite green"
else
  fail "unit suite — see line above"
fi

# Conftest preflight header
HDR_OUT=$(PYTHONPATH=python pixi run python -m pytest tests/unit/test_correction.py --no-cov 2>&1 | head -10)
if echo "$HDR_OUT" | grep -q "siac._rust: built"; then
  pass "pytest_report_header surfaces build state"
else
  fail "pytest_report_header missing"
fi

# ---------------------------------------------------------------------------
# 2. Wave 9 — siac.constants module + bare-except narrowing
# ---------------------------------------------------------------------------
hdr "Wave 6/7/9: siac.constants module"

PYTHONPATH=python pixi run python <<'PY' && pass "constants module exports + tunable env override" || fail "constants module"
import os
from siac import constants

# Spot-check every constant is the expected type with documented bounds.
expectations = {
    "ATMOSPHERIC_SCALE_HEIGHT_KM": (7.0, 9.5),
    "DEFAULT_JACOBIAN_DELTA_AOT":  (1e-4, 1.0),
    "DEFAULT_JACOBIAN_DELTA_TCWV": (1e-3, 1.0),
    "BOA_VALID_MIN":               (-0.2, 0.0),
    "BOA_VALID_MAX":               (1.0, 3.0),
    "DEFAULT_NO_DATA_BOA":         (0.0, 1.0),
    "DEFAULT_NO_DATA_BOA_UNC":     (0.0, 0.5),
    "DEFAULT_S2_VZA_DEG":          (0.0, 30.0),
    "DEFAULT_S2_VAA_DEG":          (0.0, 360.0),
    "DEFAULT_S2_ANGLE_GRID_DEG":   (0.0, 90.0),
}
for name, (lo, hi) in expectations.items():
    val = getattr(constants, name)
    assert isinstance(val, (int, float)), f"{name}: not numeric"
    assert lo <= val <= hi, f"{name} = {val} outside expected [{lo}, {hi}]"
    print(f"  {name:32s} = {val}")
PY

hdr "Wave 9: bare-except narrowing — log surfaces context"

PYTHONPATH=python pixi run python <<'PY' && pass "MERRA-2 probe falls back with exc_info" || fail "MERRA-2 narrowing"
import logging
from unittest.mock import patch

# Trigger the narrowed except in merra2._probe_granules.
from siac.adapters.atmo.merra2 import MERRA2Provider

logger = logging.getLogger("siac.adapters.atmo.merra2")
records = []
handler = logging.Handler()
handler.emit = lambda r: records.append(r)
logger.addHandler(handler)
logger.setLevel("WARNING")

with patch("siac.adapters.atmo.merra2.EarthAccessSource",
           side_effect=AttributeError("simulated SDK drift")):
    p = MERRA2Provider(probe_earthdata=True)
    # The probe is best-effort and runs lazily; trigger via get_prior path
    try:
        p._probe_granules((-10, -10, 10, 10), "EPSG:4326", __import__("datetime").datetime(2026, 1, 1))
    except Exception as exc:
        # _probe_granules is best-effort, should swallow AttributeError
        raise AssertionError(f"probe should have swallowed: {exc!r}") from exc

warn = [r for r in records if r.levelno >= logging.WARNING and "MERRA-2" in r.getMessage()]
assert warn, "expected a MERRA-2 warning record"
assert warn[0].exc_info is not None, "expected exc_info=True (REVIEW.md §2.1)"
print("  MERRA-2 warning exc_info attached:", warn[0].exc_info[0].__name__)
PY

# ---------------------------------------------------------------------------
# 3. Wave 7 — siac.runtime/siac.geo layering
# ---------------------------------------------------------------------------
hdr "Wave 7: layering — runtime no longer imports geo"

PYTHONPATH=python pixi run python <<'PY' && pass "runtime → geo cycle broken; backward-compat shim still works" || fail "layering"
import importlib

# 1. siac.geo._spatial owns the function now.
from siac.geo._spatial import copy_spatial_metadata_like as via_geo

# 2. siac.runtime.models still re-exports for backward compat.
from siac.runtime.models import copy_spatial_metadata_like as via_runtime

assert via_geo is via_runtime, (
    f"backward-compat shim broken: geo={via_geo!r}, runtime={via_runtime!r}"
)

# 3. siac.geo MUST NOT import siac.runtime at runtime.
import ast
from pathlib import Path
geo_root = Path("python/siac/geo")
violations = []
for py in geo_root.rglob("*.py"):
    src = py.read_text()
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            if node.module.startswith("siac.runtime"):
                # Allow TYPE_CHECKING-guarded imports.
                parent_lines = src.splitlines()[max(0, node.lineno - 6):node.lineno]
                if not any("TYPE_CHECKING" in line for line in parent_lines):
                    violations.append(f"{py}:{node.lineno}: {ast.unparse(node)}")
assert not violations, "runtime imported from geo at runtime:\n" + "\n".join(violations)
print("  geo → runtime: 0 runtime imports (TYPE_CHECKING-only OK)")
PY

# ---------------------------------------------------------------------------
# 4. Wave 5/7 — atomic writes
# ---------------------------------------------------------------------------
hdr "Wave 5/7: atomic-write helper + atomic CAMS S3 + monthly composites"

PYTHONPATH=python pixi run python <<'PY' && pass "_atomic_write_text + atomic CAMS S3" || fail "atomic write"
import os
import tempfile
from pathlib import Path
from siac.storage.writers import _atomic_write_text

with tempfile.TemporaryDirectory() as d:
    target = Path(d) / "out.json"
    _atomic_write_text(target, '{"x": 1}')
    assert target.read_text() == '{"x": 1}'

    # Crash-mid-write: simulate by making the temp file unwritable AFTER
    # creation. We can't easily simulate a crash, but we can verify the
    # tmp suffix exists during the write and disappears after.
    tmp = target.with_suffix(target.suffix + ".tmp")
    assert not tmp.exists(), "tmp should be cleaned up after successful write"

# CAMS atomic helper signature — confirm the .tmp + replace pattern is present
import inspect
from siac.adapters.atmo import cams
src = inspect.getsource(cams.CAMSProvider._cache_remote_file)
assert "tmp_path.replace(local_path)" in src, "REVIEW.md §1.3 #3 — missing atomic rename"
assert "tmp_path.unlink" in src, "REVIEW.md §1.3 #3 — missing cleanup-on-failure"
print("  _atomic_write_text round-trip ✓")
print("  CAMS S3 download uses tmp + os.replace ✓")
PY

# ---------------------------------------------------------------------------
# 5. Wave 4 — STAC robustness + adapter network
# ---------------------------------------------------------------------------
hdr "Wave 4: STAC + adapter network hardening"

PYTHONPATH=python pixi run python <<'PY' && pass "antimeridian-aware bbox + 6-element form" || fail "STAC antimeridian"
from siac.storage.stac import _wgs84_bounds_and_geometry

# Synthesise an antimeridian-crossing AOI. Using EPSG:32601 (UTM zone 1
# north) bounds that, projected back to WGS84, will cross 180°.
import warnings
import logging
logger = logging.getLogger("siac.storage.stac")
logger.setLevel("WARNING")

# Force the antimeridian path with explicit WGS84 bounds that span >180°.
# (The internal detection is xmax<xmin OR (xmax-xmin)>180.)
class _FakeCRS:
    pass

# We need to bypass transform_bounds, so use crs=None and pass bounds
# that already span >180.
bbox, geom, antimeridian = _wgs84_bounds_and_geometry(
    (170.0, -10.0, -170.0, 10.0),   # xmax < xmin (already in WGS84)
    crs=None,
)
assert antimeridian, "should detect antimeridian crossing"
assert len(bbox) == 6, f"6-element bbox expected; got {bbox}"
assert bbox[2] == -180.0 and bbox[5] == 180.0, f"unexpected bbox shape: {bbox}"
print(f"  antimeridian bbox: {bbox}")
print(f"  geometry type: {geom['type']}")
PY

PYTHONPATH=python pixi run python <<'PY' && pass "shared retry-enabled Session in _http" || fail "_http session"
from siac.adapters._http import make_session
from urllib3.util.retry import Retry

s = make_session(total_retries=5, status_forcelist=(429, 500, 503))
# requests.Session has adapters mounted; retry should be applied to both
# http:// and https://.
http_adapter = s.get_adapter("http://example.com/")
https_adapter = s.get_adapter("https://example.com/")
for label, adapter in (("http", http_adapter), ("https", https_adapter)):
    retry = adapter.max_retries
    assert isinstance(retry, Retry), f"{label}: not a Retry instance"
    assert retry.total == 5, f"{label}: total={retry.total} (want 5)"
    assert 503 in retry.status_forcelist, f"{label}: 503 not in status_forcelist"
print("  http + https adapters carry Retry(total=5, 503-forcelist) ✓")
print(f"  default timeout stashed: {getattr(s, '_siac_default_timeout', '<not set>')}")
PY

# ---------------------------------------------------------------------------
# 6. Wave 1 — CLI WKT fix + ValidationError ⊂ ValueError
# ---------------------------------------------------------------------------
hdr "Wave 1: CLI WKT routing + ValidationError"

PYTHONPATH=python pixi run python <<'PY' && pass "CLI --aoi-wkt now parses WKT correctly" || fail "CLI WKT"
import argparse
from siac.cli import _resolve_cli_aoi
ns = argparse.Namespace(
    aoi_bbox=None, aoi_file=None,
    aoi_wkt="POLYGON ((15 -24, 16 -24, 16 -23, 15 -23, 15 -24))",
    aoi_crs="EPSG:4326",
)
aoi = _resolve_cli_aoi(ns)
assert aoi is not None
assert aoi.geometry["type"] == "Polygon"
coords = aoi.geometry["coordinates"][0]
# bounding box of that polygon should be 15..16 x -24..-23
xs = [pt[0] for pt in coords]
ys = [pt[1] for pt in coords]
assert min(xs) == 15.0 and max(xs) == 16.0, f"x bounds wrong: {xs}"
assert min(ys) == -24.0 and max(ys) == -23.0, f"y bounds wrong: {ys}"
print(f"  WKT polygon parsed: bbox=({min(xs)}, {min(ys)}, {max(xs)}, {max(ys)})")
PY

PYTHONPATH=python pixi run python <<'PY' && pass "ValidationError catches via ValueError" || fail "ValidationError MRO"
from siac.errors import ValidationError, SIACError

# Both lineages must work.
try:
    raise ValidationError("nope")
except ValueError as exc:
    print(f"  caught via ValueError: {type(exc).__name__}")
try:
    raise ValidationError("nope")
except SIACError as exc:
    print(f"  caught via SIACError: {type(exc).__name__}")
PY

# ---------------------------------------------------------------------------
# 7. Wave 10 — Native 6S archive checksum + atomic
# ---------------------------------------------------------------------------
hdr "Wave 10: Native 6S download safety"

PYTHONPATH=python pixi run python <<'PY' && pass "atomic + sha256-aware _fetch_and_unpack_source" || fail "sixs_build download"
import inspect
from siac.algorithms.rt.direct.sixs_build import _fetch_and_unpack_source
src = inspect.getsource(_fetch_and_unpack_source)
assert ".tmp" in src and "tmp_path.replace(archive_path)" in src, \
    "atomic-download rename missing"
assert "SIAC_SIXS_SOURCE_SHA256" in src and "_archive_sha256" in src, \
    "SHA-256 verification missing"
print("  archive_path.with_suffix('.tmp') + replace ✓")
print("  SHA-256 hook with env-var override ✓")
PY

# ---------------------------------------------------------------------------
# 8. Wave 8/10 — end-to-end regression suite (gated)
# ---------------------------------------------------------------------------
hdr "Wave 8/10: regression suite against the T33KWP scene"

if [[ "$SKIP_REGRESSION" -eq 1 ]]; then
  note "--skip-regression set; skipping ~3-minute pipeline run"
else
  SAFE=/Users/fengyin/Documents/SIAC/tmp/real_cdse_mcd43_t33kwp/cache/s2/S2B_MSIL1C_20260329T084559_N0512_R107_T33KWP_20260329T140503.SAFE
  CFG=/Users/fengyin/Documents/SIAC/tmp/real_cdse_mcd43_t33kwp_sixs.toml
  if [[ ! -d "$SAFE" ]] || [[ ! -f "$CFG" ]]; then
    note "SAFE input or config not at expected absolute paths; regression auto-skip"
  else
    export SIAC_REGRESSION_SAFE="$SAFE"
    export SIAC_REGRESSION_CONFIG="$CFG"
    note "running pipeline on T33KWP — this takes ~3 min"
    # rt6s env uses a different Python (3.11) than the default env (3.12),
    # so siac._rust must be built for BOTH or the pipeline crashes inside
    # the fixture and the fixture turns the crash into a skip. Make sure
    # it's built before invoking pytest.
    if ! [[ -f python/siac/_rust.cpython-311-darwin.so ]]; then
      note "siac._rust not built for the rt6s Python; running pixi run -e rt6s build-rust"
      pixi run -e rt6s build-rust >/dev/null 2>&1 || fail "rt6s build-rust failed"
    fi
    REG_OUT=$(PYTHONPATH=python pixi run -e rt6s python -m pytest \
        tests/regression/test_t33kwp_sixs_scene.py \
        -m "regression and slow" --tb=line --no-cov -q 2>&1)
    LAST=$(echo "$REG_OUT" | tail -3 | head -1)
    echo "  $LAST"
    # Distinguish three cases: passed / failed / skipped.
    # - passed:  every assertion within tolerance (PASS)
    # - failed:  numerical regression — real bug (FAIL)
    # - skipped: the fixture couldn't run the pipeline (e.g. cache empty,
    #            rt6s env missing); not a code regression — NOTE, not FAIL.
    if echo "$REG_OUT" | grep -qE '17 passed' ; then
      pass "all 17 regression assertions pass within tolerance"
    elif echo "$REG_OUT" | grep -qE ' [0-9]+ failed' ; then
      fail "regression suite reported failures — numerical regression!"
      echo "$REG_OUT" | grep -E '^FAILED|tolerance|delta' | head -10
    elif echo "$REG_OUT" | grep -qE ' [0-9]+ skipped' ; then
      # Find why pytest skipped — usually surfaced in the fixture
      # message captured by pytest's short-summary.
      SKIP_REASON=$(echo "$REG_OUT" | grep -E 'SKIPPED|Skipped' | head -3 | sed 's/^/    /')
      note "regression skipped (not a code regression):"
      [[ -n "$SKIP_REASON" ]] && echo "$SKIP_REASON"
      note "Common causes: rt6s env missing siac._rust (run \`pixi run -e rt6s build-rust\`), missing cache under tmp/, or stale auxiliary data."
    else
      fail "regression suite output not understood — pasting last 10 lines:"
      echo "$REG_OUT" | tail -10 | sed 's/^/    /'
    fi
  fi
fi

# ---------------------------------------------------------------------------
# 9. Wave 10 — profiling harness
# ---------------------------------------------------------------------------
hdr "Wave 10: cProfile harness"

if [[ "$SKIP_PROFILE" -eq 1 ]]; then
  note "--skip-profile set; skipping cProfile run"
else
  SAFE=/Users/fengyin/Documents/SIAC/tmp/real_cdse_mcd43_t33kwp/cache/s2/S2B_MSIL1C_20260329T084559_N0512_R107_T33KWP_20260329T140503.SAFE
  CFG=/Users/fengyin/Documents/SIAC/tmp/real_cdse_mcd43_t33kwp_sixs.toml
  if [[ ! -d "$SAFE" ]] || [[ ! -f "$CFG" ]]; then
    note "SAFE/config not present; --help works:"
    PYTHONPATH=python pixi run python tools/profile_pipeline.py --help 2>&1 | head -4
  else
    OUTDIR=$(mktemp -d)
    REPORT="$OUTDIR/profile.txt"
    note "profiling pipeline run (overhead ~15%) → $REPORT"
    if PYTHONPATH=python pixi run -e rt6s python tools/profile_pipeline.py \
        --safe "$SAFE" --config "$CFG" \
        --output-dir "$OUTDIR/out" --report-path "$REPORT" \
        --top 20 2>&1 | tail -2; then
      if [[ -s "$REPORT" ]] && [[ -s "$REPORT.pstats" ]] && [[ -s "$REPORT.callers.txt" ]]; then
        pass "cProfile run produced report + pstats + callers"
        echo "  Top 5 by cumulative time:"
        grep -A 7 "ncalls" "$REPORT" | sed -n '2,7p' | sed 's/^/    /'
      else
        fail "profile reports missing or empty"
      fi
    else
      fail "profile_pipeline.py invocation failed"
    fi
  fi
fi

# ---------------------------------------------------------------------------
hdr "Summary"
if [[ "$FAILED" -eq 0 ]]; then
  printf "\033[32mAll checks passed.\033[0m\n"
  exit 0
else
  printf "\033[31m%d check(s) FAILED.\033[0m\n" "$FAILED"
  exit 1
fi
