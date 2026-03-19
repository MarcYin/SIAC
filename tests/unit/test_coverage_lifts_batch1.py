"""Targeted coverage lifts for low-coverage core modules."""

from __future__ import annotations

from datetime import datetime
from typing import TYPE_CHECKING

import numpy as np
import pytest
import xarray as xr

from siac.adapters.satellite.base import (
    BaseSatellitePreprocessor,
    detect_sensor,
    resample_angles_to_data,
)
from siac.algorithms.solver.cost import CostFunction, CostFunctionConfig, create_sparse_laplacian
from siac.domain import (
    SENTINEL2A_CONFIG,
    AtmosphericState,
    BRDFKernelWeights,
    GeometryAngles,
    RTCoefficients,
    SensorBand,
    SurfacePrior,
)
from siac.domain.protocols import (
    AerosolSolver,
    AtmosphericPriorProvider,
    BRDFProductProvider,
    OutputWriter,
    RTModelBackend,
    SatellitePreprocessor,
    SurfacePriorDeriver,
)
from siac.domain.validation import _spatial_shape

if TYPE_CHECKING:
    from pathlib import Path

# ------------------------------
# protocol line execution
# ------------------------------

def test_protocol_method_stubs_execute_lines():
    # Property/method stubs with ellipsis bodies are executable; invoke each once.
    assert SatellitePreprocessor.sensor_config.fget(object()) is None
    assert SatellitePreprocessor.load_toa(object(), "/tmp") is None
    assert SatellitePreprocessor.extract_geometry(object(), "/tmp") is None
    assert SatellitePreprocessor.extract_cloud_mask(object(), "/tmp") is None
    assert SatellitePreprocessor.get_metadata(object(), "/tmp") is None

    assert AtmosphericPriorProvider.get_prior(object(), (0, 0, 1, 1), "EPSG:4326", datetime(2024, 1, 1), 1000.0) is None
    assert AtmosphericPriorProvider.source_name.fget(object()) is None

    assert BRDFProductProvider.get_brdf_parameters(
        object(), (0, 0, 1, 1), "EPSG:4326", datetime(2024, 1, 1), 500.0, [1, 2, 3], 16
    ) is None
    assert BRDFProductProvider.source_name.fget(object()) is None

    assert SurfacePriorDeriver.compute_surface_prior(object(), object(), object(), None) is None

    assert RTModelBackend.compute_coefficients(object(), object(), object(), object(), False) is None
    assert RTModelBackend.supports_jacobian(object()) is None
    assert RTModelBackend.backend_name.fget(object()) is None
    assert RTModelBackend.is_available_for_sensor(object(), "MSI", "S2A") is None

    assert AerosolSolver.solve(object(), object()) is None

    assert OutputWriter.write_boa(object(), xr.Dataset(), "/tmp/out", "EPSG:4326", (0, 1, 0, 0, 0, -1)) is None
    assert OutputWriter.write_auxiliary(
        object(), xr.DataArray(np.zeros((1, 1))), xr.DataArray(np.zeros((1, 1))), "/tmp/out", "EPSG:4326", (0, 1, 0, 0, 0, -1)
    ) is None


# ------------------------------
# validation fallback shape path
# ------------------------------

def test_spatial_shape_fallback_and_error():
    ds_3d = xr.Dataset({"v": xr.DataArray(np.zeros((2, 3, 4), dtype=np.float32), dims=["b", "row", "col"])})
    assert _spatial_shape(ds_3d) == (3, 4)

    ds_1d = xr.Dataset({"v": xr.DataArray(np.zeros((2,), dtype=np.float32), dims=["x"])})
    with pytest.raises(ValueError, match="fewer than 2"):
        _ = _spatial_shape(ds_1d)


# ------------------------------
# satellite.base utility branches
# ------------------------------

class _DummyPreprocessor(BaseSatellitePreprocessor):
    @property
    def sensor_config(self):
        return SENTINEL2A_CONFIG

    def load_toa(self, input_path: str | Path) -> xr.Dataset:
        del input_path
        return xr.Dataset({"B02": xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"])})

    def extract_geometry(self, input_path: str | Path) -> GeometryAngles:
        del input_path
        arr = xr.DataArray(np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"])
        return GeometryAngles(sza=arr, saa=arr, vza=arr, vaa=arr)

    def extract_cloud_mask(self, input_path: str | Path) -> xr.DataArray:
        del input_path
        return xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"])

    def get_metadata(self, input_path: str | Path) -> dict:
        del input_path
        return {"observation_time": datetime(2024, 1, 1)}


def test_base_preprocessor_and_sensor_detection_branches(tmp_path: Path, monkeypatch):
    pp = _DummyPreprocessor()
    out = pp.preprocess(tmp_path)
    assert set(out) == {"toa", "geometry", "cloud_mask", "metadata"}

    # resample_angles_to_data delegates to reproject_match
    target = xr.DataArray(np.zeros((3, 3), dtype=np.float32), dims=["y", "x"])
    called = {}

    def _fake_reproject_match(angles, tgt, resampling="linear"):
        called["resampling"] = resampling
        return xr.full_like(tgt, 42.0)

    monkeypatch.setattr("siac.io.reprojection.reproject_match", _fake_reproject_match)
    res = resample_angles_to_data(out["geometry"].sza, target, method="nearest")
    assert called["resampling"] == "nearest"
    assert float(res.mean()) == 42.0

    # Landsat 5/7 branch (previously uncovered)
    l5_dir = tmp_path / "l5"
    l5_dir.mkdir()
    (l5_dir / "foo_MTL.txt").write_text("LANDSAT_5")
    assert detect_sensor(l5_dir) == "l5"


# ------------------------------
# solver.cost uncovered helpers
# ------------------------------

def _geom(shape=(3, 3)) -> GeometryAngles:
    return GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
        saa=xr.DataArray(np.full(shape, 2.0), dims=["y", "x"]),
        vza=xr.DataArray(np.full(shape, 0.2), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
    )


def _atmo(shape=(3, 3)) -> AtmosphericState:
    def da(v):
        return xr.DataArray(np.full(shape, v, dtype=np.float32), dims=["y", "x"])
    return AtmosphericState(aot=da(0.15), tcwv=da(2.0), tco3=da(0.3), aot_unc=da(0.05), tcwv_unc=da(0.3), tco3_unc=da(0.03), elevation=da(0.1))


def _surface(shape=(3, 3)) -> SurfacePrior:
    def da(v):
        return xr.DataArray(np.full(shape, v, dtype=np.float32), dims=["y", "x"])
    kernels = BRDFKernelWeights(f0=da(0.1), f1=da(0.05), f2=da(0.02), f0_unc=da(0.01), f1_unc=da(0.01), f2_unc=da(0.01))
    return SurfacePrior(boa=da(0.2), boa_unc=da(0.1), kernels=kernels, mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]))


class _NoMultiRT:
    def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):
        del atmo_state, band, compute_jacobian
        shape = geometry.sza.shape
        xap = xr.DataArray(np.full(shape, 1.0), dims=["y", "x"])
        xbp = xr.DataArray(np.full(shape, 0.0), dims=["y", "x"])
        xcp = xr.DataArray(np.full(shape, 0.0), dims=["y", "x"])
        d = xr.concat([
            xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            xr.DataArray(np.full(shape, 0.2), dims=["y", "x"]),
        ], dim="param").assign_coords(param=["aot", "tcwv"])
        return RTCoefficients(xap=xap, xbp=xbp, xcp=xcp, d_xap=d, d_xbp=d, d_xcp=d)

    def supports_jacobian(self):
        return True

    @property
    def backend_name(self):
        return "fake"

    def is_available_for_sensor(self, sensor_id, satellite_id):
        del sensor_id, satellite_id
        return True


def test_cost_function_fallback_compute_coeffs_and_helper_methods():
    shape = (3, 3)
    toa = xr.DataArray(np.full((1, *shape), 0.3, dtype=np.float32), dims=["band", "y", "x"])
    surface = _surface(shape)
    geom = _geom(shape)
    atmo = _atmo(shape)
    mask = xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"])
    bands = [SensorBand("B02", 490.0, 65.0, 10.0, 0)]

    cfg = CostFunctionConfig()
    cf = CostFunction(toa, surface, geom, atmo, _NoMultiRT(), bands, mask, cfg)

    # Use helper diagnostic methods (were partially uncovered).
    j_obs, g_obs = cf.observation_cost_only(atmo.aot.values, atmo.tcwv.values)
    j_prior, g_prior = cf.prior_cost_only(atmo.aot.values, atmo.tcwv.values)
    j_sm, g_sm = cf.smoothness_cost_only(atmo.aot.values, atmo.tcwv.values)

    assert np.isfinite(j_obs)
    assert np.isfinite(j_prior)
    assert np.isfinite(j_sm)
    assert g_obs.size == 2 * np.prod(shape)
    assert g_prior.size == 2 * np.prod(shape)
    assert g_sm.size == 2 * np.prod(shape)

    lap = create_sparse_laplacian(4, 3)
    assert lap.shape == (12, 12)
    assert lap.nnz > 0
