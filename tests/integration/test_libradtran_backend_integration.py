"""Integration test that actually runs the libRadtran ``uvspec`` engine.

Skip-guarded: requires a compiled/available libRadtran (via the cached build,
an explicit ``uvspec_path``, or ``uvspec`` on PATH). Build it with
``pixi run -e libradtran build-libradtran`` to enable this test locally.
"""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.rt.direct.libradtran import LibRadtranBackend
from siac.algorithms.rt.direct.libradtran_build import ensure_libradtran
from siac.config.algorithms import LibRadtranAlgorithmConfig
from siac.domain.sensors import SensorBand
from siac.runtime import AtmosphericState, GeometryAngles, RTCoefficients

pytestmark = pytest.mark.integration


def _require_libradtran() -> LibRadtranAlgorithmConfig:
    """Return a config if a libRadtran engine is available, else skip."""
    # Keep the test quick: visible+NIR window, a 2x2 aot/tcwv grid.
    cfg = LibRadtranAlgorithmConfig(
        wavelength_min_nm=400.0,
        wavelength_max_nm=900.0,
        scene_lut_max_nodes_per_axis=2,
        scene_lut_max_cases=8,
        native_threads=4,
        auto_build=False,  # do not trigger a multi-GB download/compile in CI
    )
    try:
        ensure_libradtran(cfg)
    except RuntimeError as exc:  # not built / unavailable
        pytest.skip(f"libRadtran engine unavailable: {exc}")
    return cfg


def _da(vals: list[list[float]]) -> xr.DataArray:
    return xr.DataArray(np.asarray(vals, dtype=np.float64), dims=("y", "x"))


def _scene() -> tuple[GeometryAngles, AtmosphericState]:
    # The reused ZarrLUTBackend requires finite geometry AND atmosphere (it
    # rejects NaNs in the interpolation inputs, unlike 6S which masks). In real
    # runs the solver always feeds finite per-pixel atmo to the RT backend, so
    # this matches production usage.
    geom = GeometryAngles.from_degrees(
        _da([[30.0, 31.0], [30.5, 30.2]]),
        _da([[150.0, 150.0], [150.0, 150.0]]),
        _da([[5.0, 5.0], [5.0, 5.0]]),
        _da([[110.0, 110.0], [110.0, 110.0]]),
    )
    zeros = _da([[0.0, 0.0], [0.0, 0.0]])
    atmo = AtmosphericState(
        aot=_da([[0.15, 0.25], [0.20, 0.30]]),
        tcwv=_da([[1.2, 1.6], [1.4, 1.5]]),
        tco3=_da([[0.30, 0.30], [0.30, 0.30]]),
        aot_unc=zeros,
        tcwv_unc=zeros,
        tco3_unc=zeros,
        elevation=_da([[0.05, 0.05], [0.05, 0.05]]),
    )
    return geom, atmo


def test_libradtran_backend_produces_physical_coefficients() -> None:
    cfg = _require_libradtran()
    backend = LibRadtranBackend(libradtran_config=cfg)
    geom, atmo = _scene()
    band = SensorBand(
        name="B02", center_wavelength=492.0, bandwidth=66.0, resolution=10.0, band_index=1
    )

    coeffs = backend.compute_coefficients(geom, atmo, band)
    assert isinstance(coeffs, RTCoefficients)
    assert coeffs.xap.shape == (2, 2)

    xap = np.asarray(coeffs.xap.values)
    xcp = np.asarray(coeffs.xcp.values)
    # Physical ranges for a visible band at modest AOT.
    finite = np.isfinite(xap)
    assert finite.any()
    assert np.all(xap[finite] > 1.0)  # xap = 1 / T_total, T_total < 1
    assert np.all(xap[finite] < 5.0)
    assert np.all(xcp[np.isfinite(xcp)] >= 0.0)  # spherical albedo non-negative

    # Multi-band reuses the cached scene LUT (no second grid build).
    bands = [
        band,
        SensorBand(
            name="B04", center_wavelength=665.0, bandwidth=31.0, resolution=10.0, band_index=3
        ),
    ]
    multi = backend.compute_coefficients_multi(geom, atmo, bands)
    assert len(multi) == 2
    assert all(isinstance(c, RTCoefficients) for c in multi)


def test_libradtran_hybrid_mol_abs_regions_runs() -> None:
    """End-to-end per-region run: coarse base + a fine 0.94um segment, stitched."""
    base = _require_libradtran()
    cfg = base.model_copy(
        update={
            "mol_abs_param": "reptran",  # coarse base (bundled, fast)
            "wavelength_min_nm": 400.0,
            "wavelength_max_nm": 1050.0,
            "mol_abs_regions": [(915.0, 985.0, "reptran fine")],
        }
    )
    backend = LibRadtranBackend(libradtran_config=cfg)
    geom, atmo = _scene()
    # B09 (945 nm) falls inside the fine segment; B08 (833 nm) in the coarse base.
    for name, center, bw, idx in [("B08", 833.0, 106.0, 7), ("B09", 945.0, 26.0, 9)]:
        band = SensorBand(
            name=name, center_wavelength=center, bandwidth=bw, resolution=20.0, band_index=idx
        )
        coeffs = backend.compute_coefficients(geom, atmo, band)
        xap = np.asarray(coeffs.xap.values)
        assert np.all(np.isfinite(xap[np.isfinite(atmo.aot.values)]))
        assert np.all(xap[np.isfinite(xap)] > 1.0)
