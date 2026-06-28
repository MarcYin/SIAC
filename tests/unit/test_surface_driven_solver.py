"""Tests for the surface-driven aerosol solver and its Rust kernel."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

from siac._rust_compat import surface_driven_pool_argmin
from siac.app.registry import SOLVER_METHOD_REGISTRY
from siac.config.types import SolverMethod


def _aot_axis(n: int = 40) -> np.ndarray:
    return np.geomspace(0.001, 2.5, n).astype(np.float64)


# --------------------------------------------------------------------------- #
# Rust kernel: surface_driven_pool_argmin (f64; window=1 -> no pooling)
# --------------------------------------------------------------------------- #
class TestPoolArgmin:
    def test_picks_min_node_without_backstop(self) -> None:
        ny, nx, n_aot = 3, 3, 5
        axis = np.linspace(0.05, 0.45, n_aot)
        cube = np.ones((n_aot, ny, nx), dtype=np.float64)
        cube[2] = 0.1  # node 2 is the clear minimum everywhere
        prior = np.full((ny, nx), 0.2)
        prior_unc = np.full((ny, nx), 1e6)  # backstop ~off
        valid = np.ones((ny, nx), dtype=bool)

        aot, unc, jmin = surface_driven_pool_argmin(cube, axis, prior, prior_unc, valid, 1, 1)
        aot = np.asarray(aot)
        assert np.allclose(aot, axis[2])
        assert np.allclose(np.asarray(jmin), 0.1)  # obs-only min pooled cost
        assert np.all(np.asarray(unc) >= 0.02)

    def test_backstop_pulls_flat_cost_to_prior(self) -> None:
        ny, nx, n_aot = 2, 2, 7
        axis = np.linspace(0.05, 0.65, n_aot)
        cube = np.full((n_aot, ny, nx), 0.5, dtype=np.float64)  # flat surface cost
        prior_val = float(axis[4])
        prior = np.full((ny, nx), prior_val)
        prior_unc = np.full((ny, nx), 0.02)  # tight backstop
        valid = np.ones((ny, nx), dtype=bool)

        aot, _unc, _jmin = surface_driven_pool_argmin(cube, axis, prior, prior_unc, valid, 1, 1)
        assert np.allclose(np.asarray(aot), prior_val)

    def test_all_nodes_finite_required(self) -> None:
        # numpy resolve medians only pixels where every AOD node is pooled-finite
        # (isfinite(pooled).all(0)); a pixel missing any node -> NaN.
        ny, nx, n_aot = 2, 2, 4
        axis = _aot_axis(n_aot)
        cube = np.empty((n_aot, ny, nx), dtype=np.float64)
        cube[:] = 0.5
        cube[2] = [[0.5, 0.5], [0.5, 0.5]]
        cube[2, 0, 1] = np.nan  # pixel (0,1) missing node 2 -> not all-finite
        prior = np.full((ny, nx), 0.1)
        prior_unc = np.full((ny, nx), 1e6)
        valid = np.ones((ny, nx), dtype=bool)
        valid[1, 1] = False  # explicitly invalid pixel

        aot = np.asarray(surface_driven_pool_argmin(cube, axis, prior, prior_unc, valid, 1, 1)[0])
        assert np.isfinite(aot[0, 0])  # all nodes finite
        assert np.isnan(aot[0, 1])  # missing a node
        assert np.isnan(aot[1, 1])  # invalid

    def test_spatial_pooling_smooths_outlier(self) -> None:
        ny, nx, n_aot = 3, 3, 3
        axis = np.array([0.1, 0.2, 0.3], dtype=np.float64)
        # Node 0 is best everywhere except a single-pixel spike that would
        # otherwise pick node 2; median pooling over a 3x3 window heals it.
        cube = np.zeros((n_aot, ny, nx), dtype=np.float64)
        cube[0] = 0.1
        cube[1] = 0.5
        cube[2] = 0.9
        cube[2, 1, 1] = 0.0  # spike makes node 2 look best at the centre pixel
        cube[0, 1, 1] = 0.8
        prior = np.full((ny, nx), 0.2)
        prior_unc = np.full((ny, nx), 1e6)
        valid = np.ones((ny, nx), dtype=bool)

        # window=3 (lo=hi=1): 3x3 centered median, min_periods=1.
        aot_pooled = np.asarray(
            surface_driven_pool_argmin(cube, axis, prior, prior_unc, valid, 3, 1)[0]
        )
        assert np.isclose(aot_pooled[1, 1], axis[0])  # outlier healed


# --------------------------------------------------------------------------- #
# SurfaceDrivenSolver end-to-end
# --------------------------------------------------------------------------- #
class _DummyRT:
    """RT model whose BOA correction is ``boa = toa - aot`` (xbp = aot)."""

    def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
        del band, compute_jacobian
        a = float(np.asarray(atmo_state.aot.values).reshape(-1)[0])
        shp = tuple(geometry.sza.shape)
        from siac.runtime.models import RTCoefficients

        return RTCoefficients(
            xap=xr.DataArray(np.ones(shp, dtype=np.float32), dims=["y", "x"]),
            xbp=xr.DataArray(np.full(shp, a, dtype=np.float32), dims=["y", "x"]),
            xcp=xr.DataArray(np.zeros(shp, dtype=np.float32), dims=["y", "x"]),
        )

    @property
    def backend_name(self) -> str:
        return "dummy"

    def is_available_for_sensor(self, sensor_id, satellite_id):  # noqa: ANN001
        del sensor_id, satellite_id
        return True


def _make_inputs(shape=(4, 4), toa_val=0.30, prior_boa=0.10):
    from siac.domain import SensorBand
    from siac.runtime import AtmosphericState, GeometryAngles, SurfacePrior

    toa = xr.DataArray(np.full((1, *shape), toa_val, dtype=np.float32), dims=["band", "y", "x"])
    cloud_mask = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"])
    surface_prior = SurfacePrior(
        boa=xr.DataArray(np.full(shape, prior_boa, dtype=np.float32), dims=["y", "x"]),
        boa_unc=xr.DataArray(np.full(shape, 0.02, dtype=np.float32), dims=["y", "x"]),
        kernels=None,
        mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
    )
    geometry = GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
        saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
        vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
    )
    atmo_prior = AtmosphericState(
        aot=xr.DataArray(np.full(shape, 0.15, dtype=np.float32), dims=["y", "x"]),
        tcwv=xr.DataArray(np.full(shape, 2.5, dtype=np.float32), dims=["y", "x"]),
        tco3=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
        aot_unc=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=["y", "x"]),
        tcwv_unc=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
        tco3_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=["y", "x"]),
        elevation=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
    )
    bands = [SensorBand("B02", 490.0, 65.0, 10.0, 0)]
    return toa, surface_prior, geometry, atmo_prior, cloud_mask, bands


def _solver_config(**overrides):
    from siac.config.algorithms import SolverAlgorithmConfig

    base = {
        "method": "surface_driven",
        "grid_search_aot_points": 48,
        "surface_driven_backstop_calibrated": False,
        "surface_driven_pool_radius_m": 0.0,
        "aerosol_resolution": 120.0,
    }
    base.update(overrides)
    return SolverAlgorithmConfig(**base)


class TestSurfaceDrivenSolver:
    def test_recovers_surface_min_aot(self) -> None:
        from siac.algorithms.solver import SurfaceDrivenSolver

        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs(
            toa_val=0.30, prior_boa=0.10
        )
        solver = SurfaceDrivenSolver(_solver_config())
        result = solver.solve(
            toa, surface_prior, geometry, atmo_prior, _DummyRT(), cloud_mask, bands
        )
        assert result.success
        aot = result.aot.values
        assert np.all(np.isfinite(aot))
        # boa = toa - aot; mismatch minimised at aot = toa - prior = 0.20.
        assert abs(float(np.median(aot)) - 0.20) < 0.06
        assert result.tcwv.values == pytest.approx(2.5)

    def test_solve_bundle_matches_solve(self) -> None:
        from siac.algorithms.solver import SurfaceDrivenSolver
        from siac.runtime import SolverInputBundle

        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs()
        bundle = SolverInputBundle(
            toa=toa,
            geometry=geometry,
            cloud_mask=cloud_mask,
            sensor_config=None,
            bands=bands,
            atmo_prior=atmo_prior,
            surface_prior=surface_prior,
            rt_model=_DummyRT(),
            aux_resolution_m=480.0,
            aerosol_resolution_m=120.0,
        )
        solver = SurfaceDrivenSolver(_solver_config())
        result = solver.solve_bundle(bundle)
        assert result.success
        assert result.aot.shape == (4, 4)

    def test_water_mask_pixels_fall_back_to_prior(self) -> None:
        from siac.algorithms.solver import SurfaceDrivenSolver

        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs()
        water = np.zeros((4, 4), dtype=bool)
        water[0, 0] = True  # one water pixel -> excluded from the solve
        water_mask = xr.DataArray(water, dims=["y", "x"])
        solver = SurfaceDrivenSolver(_solver_config())
        result = solver.solve(
            toa,
            surface_prior,
            geometry,
            atmo_prior,
            _DummyRT(),
            cloud_mask,
            bands,
            water_mask=water_mask,
        )
        assert result.success
        # Excluded pixel falls back to the AOT prior (0.15); solved pixels ~0.20.
        assert result.aot.values[0, 0] == pytest.approx(0.15, abs=1e-4)
        assert abs(float(result.aot.values[3, 3]) - 0.20) < 0.06

    def test_calibrated_backstop_branch(self) -> None:
        from siac.algorithms.solver import SurfaceDrivenSolver

        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs()
        solver = SurfaceDrivenSolver(_solver_config(surface_driven_backstop_calibrated=True))
        result = solver.solve(
            toa, surface_prior, geometry, atmo_prior, _DummyRT(), cloud_mask, bands
        )
        assert result.success and np.all(np.isfinite(result.aot.values))


# --------------------------------------------------------------------------- #
# Registry dispatch (config.algorithms.solver.method)
# --------------------------------------------------------------------------- #
class TestSolverDispatch:
    def test_registry_has_both_methods(self) -> None:
        # Registration happens as an import side-effect of the assembly module
        # (the canonical trigger, as for every other SIAC registry).
        import siac.app._assembly_solver  # noqa: F401

        assert "surface_driven" in SOLVER_METHOD_REGISTRY.names()
        assert "multigrid" in SOLVER_METHOD_REGISTRY.names()

    def test_resolve_solver_dispatches_surface_driven(self) -> None:
        from siac.app._assembly_solver import resolve_solver

        cfg = SimpleNamespace(
            algorithms=SimpleNamespace(solver=_solver_config(method=SolverMethod.SURFACE_DRIVEN))
        )
        fn = resolve_solver(cfg)
        assert callable(fn)
