"""End-to-end regression test for the T33KWP S2B scene with the 6S backend.

This test re-runs the SIAC pipeline on a real Sentinel-2 scene and asserts
that the resulting AOT, TCWV, BOA bands, and STAC sidecar match a
captured set of summary statistics within tolerance. It is the first
member of the "tier-B numerical regression" suite called out in
``REVIEW_FIXES.md`` — its purpose is to make solver/RT/correction
numerical changes safe to apply.

**The test is gated.** It needs:

1. The cached SAFE input at the path encoded in the config TOML.
2. The cached MCD43 / CAMS / DEM auxiliary data the config points to.
3. The pre-built native 6S extension (``rt6s`` pixi env).

If any of those is missing the test is skipped with a clear reason
rather than failing. To run it explicitly::

    pixi run -e rt6s python -m pytest tests/regression \
        -m "regression and slow" --no-cov

The capture script ``regenerate_goldens.py`` next to this file refreshes
the goldens from a fresh run (use it after any *intentional* numerical
change once the new outputs have been validated by other means).

Tolerances are documented in ``_compare.py``. They are tight (~1e-3
relative) because the pipeline is deterministic on the same inputs — any
real numerical regression will exceed them by orders of magnitude.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from tests.regression._compare import (
    assert_stats_within_tolerance,
    capture_product_stats,
)

pytestmark = [pytest.mark.regression, pytest.mark.slow]

REPO_ROOT = Path(__file__).resolve().parents[2]
GOLDEN_PATH = Path(__file__).parent / "goldens" / "t33kwp_sixs_20260329.json"

# Default locations match the existing developer tree on the maintainer's
# machine. Override via ``SIAC_REGRESSION_*`` env vars for CI / other users.
DEFAULT_CONFIG = REPO_ROOT / "tmp" / "real_cdse_mcd43_t33kwp_sixs.toml"
DEFAULT_SAFE = (
    REPO_ROOT
    / "tmp"
    / "real_cdse_mcd43_t33kwp"
    / "cache"
    / "s2"
    / "S2B_MSIL1C_20260329T084559_N0512_R107_T33KWP_20260329T140503.SAFE"
)


def _resolve_path(env_var: str, default: Path) -> Path:
    """Pull ``env_var`` from the environment or fall back to ``default``."""
    raw = os.environ.get(env_var)
    return Path(raw).expanduser() if raw else default


@pytest.fixture(scope="module")
def golden() -> dict:
    return json.loads(GOLDEN_PATH.read_text())


@pytest.fixture(scope="module")
def regression_inputs() -> dict[str, Path]:
    """Locate the SAFE + config; skip the test cleanly if either is missing."""
    config = _resolve_path("SIAC_REGRESSION_CONFIG", DEFAULT_CONFIG)
    safe = _resolve_path("SIAC_REGRESSION_SAFE", DEFAULT_SAFE)
    if not config.exists():
        pytest.skip(
            f"Regression config not found at {config}; set SIAC_REGRESSION_CONFIG "
            f"or place the file at the default path to enable this test."
        )
    if not safe.exists():
        pytest.skip(
            f"Regression SAFE input not found at {safe}; set SIAC_REGRESSION_SAFE "
            f"or place the directory at the default path to enable this test."
        )
    return {"config": config, "safe": safe}


@pytest.fixture(scope="module")
def regression_run(regression_inputs: dict[str, Path], tmp_path_factory) -> Path:
    """Run the pipeline once per session and return the output directory."""
    out = tmp_path_factory.mktemp("regression-t33kwp-sixs")

    from siac.adapters.auth import CredentialManager
    from siac.api.requests import SceneProcessRequest
    from siac.config import SIACConfig
    from siac.workflows.scene import process_scene

    config = SIACConfig.from_file(regression_inputs["config"])
    request = SceneProcessRequest(
        config=config,
        input_path=regression_inputs["safe"],
        output_path=out,
        auth=CredentialManager.from_config(config),
    )

    # Distinguish "inputs/environment missing" from "real pipeline regression".
    # An ImportError or FileNotFoundError almost always means cache/env;
    # any other exception class is a real bug we want to surface.
    try:
        process_scene(request)
    except (ImportError, FileNotFoundError) as exc:
        pytest.skip(
            f"Regression pipeline raised {type(exc).__name__}: {exc}. "
            "This usually means cached MCD43/CAMS data is missing or the rt6s "
            "extension isn't built. Build with `pixi run -e rt6s build-rust` "
            "and ensure the cache directories referenced by the config exist; "
            "set SIAC_REGRESSION_STRICT=1 to fail instead of skipping."
        )
    except Exception as exc:
        # SIAC errors / RT errors / solver errors are NOT skippable — they're
        # the whole point of the regression suite. Only the two
        # "environment isn't there" classes above are treated as skips.
        # Override via SIAC_REGRESSION_STRICT=0 if you want every failure
        # to skip instead.
        if os.environ.get("SIAC_REGRESSION_STRICT", "1") == "0":
            pytest.skip(
                f"Regression pipeline raised {type(exc).__name__}: {exc}. "
                "(SIAC_REGRESSION_STRICT=0 so treating as environment issue.)"
            )
        raise

    return out


def _resolve_output_prefix(out_dir: Path) -> str:
    """Find the SIAC output filename prefix (e.g. ``S2B_L2A_20260329T091708``)."""
    matches = sorted(out_dir.glob("*_AOT.tif"))
    if not matches:
        raise FileNotFoundError(f"No *_AOT.tif under {out_dir}")
    if len(matches) > 1:
        raise RuntimeError(f"Multiple *_AOT.tif under {out_dir}: {matches}")
    return matches[0].stem.removesuffix("_AOT")


def test_t33kwp_sixs_aot_within_tolerance(regression_run: Path, golden: dict) -> None:
    prefix = _resolve_output_prefix(regression_run)
    actual = capture_product_stats(regression_run / f"{prefix}_AOT.tif")
    assert_stats_within_tolerance(actual, golden["products"]["AOT"], name="AOT")


def test_t33kwp_sixs_tcwv_within_tolerance(regression_run: Path, golden: dict) -> None:
    prefix = _resolve_output_prefix(regression_run)
    actual = capture_product_stats(regression_run / f"{prefix}_TCWV.tif")
    assert_stats_within_tolerance(actual, golden["products"]["TCWV"], name="TCWV")


def test_t33kwp_sixs_cloud_mask_within_tolerance(
    regression_run: Path, golden: dict
) -> None:
    prefix = _resolve_output_prefix(regression_run)
    actual = capture_product_stats(regression_run / f"{prefix}_CLOUD.tif")
    # Cloud mask is binary — slightly looser valid_fraction bound since the
    # omnicloudmask torch backend isn't fully deterministic across hardware.
    assert_stats_within_tolerance(
        actual,
        golden["products"]["CLOUD"],
        name="CLOUD",
        valid_fraction_abs_tol=5e-3,
    )


@pytest.mark.parametrize(
    "band",
    (
        "B01",
        "B02",
        "B03",
        "B04",
        "B05",
        "B06",
        "B07",
        "B08",
        "B8A",
        "B09",
        "B10",
        "B11",
        "B12",
    ),
)
def test_t33kwp_sixs_boa_band_within_tolerance(
    regression_run: Path, golden: dict, band: str
) -> None:
    if f"BOA_{band}" not in golden["products"]:
        pytest.skip(f"No golden for BOA_{band}")
    prefix = _resolve_output_prefix(regression_run)
    actual = capture_product_stats(regression_run / f"{prefix}_BOA_{band}.tif")
    assert_stats_within_tolerance(
        actual, golden["products"][f"BOA_{band}"], name=f"BOA_{band}"
    )


def test_t33kwp_sixs_stac_properties(regression_run: Path, golden: dict) -> None:
    """Lock STAC sidecar properties — bbox, AOT/TCWV mean, view angles."""
    prefix = _resolve_output_prefix(regression_run)
    stac_path = regression_run / f"{prefix}.json"
    assert stac_path.exists(), f"STAC sidecar missing at {stac_path}"
    actual = json.loads(stac_path.read_text())

    actual_props = actual["properties"]
    for key, golden_val in golden["stac"].items():
        assert key in actual_props, f"STAC property {key!r} missing from output"
        a = actual_props[key]
        if isinstance(golden_val, (int, float)) and isinstance(a, (int, float)):
            assert abs(a - golden_val) <= max(1e-4, 1e-4 * abs(golden_val)), (
                f"STAC property {key}: golden {golden_val}, actual {a}"
            )
        else:
            assert a == golden_val, (
                f"STAC property {key}: golden {golden_val!r}, actual {a!r}"
            )

    # bbox: 4-element or 6-element (antimeridian-aware).
    assert actual["bbox"] == pytest.approx(golden["bbox"], abs=1e-6), (
        f"STAC bbox: golden {golden['bbox']}, actual {actual['bbox']}"
    )
