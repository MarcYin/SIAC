"""Helpers for comparing a pipeline output against a captured golden snapshot.

The functions here are factored out of ``test_t33kwp_sixs_scene.py`` so the
same comparison logic is reusable from the ``regenerate_goldens.py`` helper
and any future regression tests against other scenes.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from pathlib import Path

import numpy as np

# Statistic keys we lock per-product. Every product in the goldens carries
# these (except where ``all_masked`` is True). Keep this list in sync with
# ``capture_product_stats`` below and with ``regenerate_goldens.py``.
_PRODUCT_STAT_KEYS: tuple[str, ...] = (
    "valid_fraction",
    "mean",
    "std",
    "p01",
    "p50",
    "p99",
    "min",
    "max",
)

# Default tolerances for a deterministic-on-same-input pipeline. The pipeline
# IS bit-deterministic in principle, but float reductions across threads
# (``np.mean`` over many cells, OpenMP-backed Rust kernels) can shift the
# last few bits. ``rel`` covers magnitudes that matter (~1e-2 and up);
# ``abs_floor`` prevents tight relative tolerances from biting on values
# that are essentially zero.
DEFAULT_REL_TOL: float = 1e-3
DEFAULT_ABS_TOL: float = 1e-4
DEFAULT_VALID_FRACTION_ABS_TOL: float = 1e-3


def capture_product_stats(path: Path) -> dict[str, Any]:
    """Read a single-band raster and return the same stats dict the goldens use.

    Returns a dict with keys ``shape``, ``dtype``, ``valid_fraction``, and
    distribution stats. If the raster is fully masked, returns
    ``{"shape": ..., "all_masked": True}``.
    """
    import rasterio

    with rasterio.open(path) as ds:
        arr = ds.read(1, masked=True)
        finite = np.ma.compressed(arr)
        if finite.size == 0:
            return {"shape": list(arr.shape), "all_masked": True}
        return {
            "shape": list(arr.shape),
            "dtype": str(arr.dtype),
            "valid_fraction": float(finite.size / arr.size),
            "mean": float(np.mean(finite)),
            "std": float(np.std(finite)),
            "p01": float(np.percentile(finite, 1)),
            "p50": float(np.percentile(finite, 50)),
            "p99": float(np.percentile(finite, 99)),
            "min": float(np.min(finite)),
            "max": float(np.max(finite)),
        }


def assert_stats_within_tolerance(
    actual: dict[str, Any],
    golden: dict[str, Any],
    *,
    name: str,
    rel: float = DEFAULT_REL_TOL,
    abs_tol: float = DEFAULT_ABS_TOL,
    valid_fraction_abs_tol: float = DEFAULT_VALID_FRACTION_ABS_TOL,
) -> None:
    """Assert ``actual`` stats match ``golden`` within tolerance.

    Raises ``AssertionError`` with the offending key on first mismatch.
    Shape and dtype must match exactly. ``valid_fraction`` uses
    ``valid_fraction_abs_tol`` (default 1e-3 — allows a handful of pixels
    to differ from cloud-mask drift). All distribution stats use
    ``np.isclose(rtol=rel, atol=abs_tol)``.
    """
    assert actual["shape"] == golden["shape"], (
        f"{name}: shape mismatch — golden {golden['shape']}, actual {actual['shape']}"
    )
    if golden.get("all_masked") or actual.get("all_masked"):
        assert actual.get("all_masked") == golden.get("all_masked"), (
            f"{name}: one is fully masked but the other is not"
        )
        return

    assert actual["dtype"] == golden["dtype"], (
        f"{name}: dtype mismatch — golden {golden['dtype']!r}, actual {actual['dtype']!r}"
    )

    delta_vf = abs(actual["valid_fraction"] - golden["valid_fraction"])
    assert delta_vf <= valid_fraction_abs_tol, (
        f"{name}: valid_fraction drift {delta_vf:.2e} > "
        f"{valid_fraction_abs_tol:.2e} "
        f"(golden {golden['valid_fraction']:.6f}, "
        f"actual {actual['valid_fraction']:.6f})"
    )

    for key in _PRODUCT_STAT_KEYS:
        if key == "valid_fraction":
            continue  # already handled above
        g = golden[key]
        a = actual[key]
        if not np.isclose(a, g, rtol=rel, atol=abs_tol, equal_nan=True):
            raise AssertionError(
                f"{name}.{key}: golden {g:.8f}, actual {a:.8f} "
                f"(delta {abs(a - g):.3e}, rtol={rel}, atol={abs_tol})"
            )
