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


def test_cost_curve_diagnostics_marks_edge_minimum() -> None:
    from siac.algorithms.solver.surface_driven import _cost_curve_diagnostics

    axis = np.array([0.1, 0.2, 0.3], dtype=np.float64)
    cube = np.ones((3, 2, 2), dtype=np.float64)
    cube[0] = 0.1
    cube[1] = 0.2
    cube[2] = 0.3

    diagnostics = _cost_curve_diagnostics(cube, axis)

    assert diagnostics["surface_cost_curve_min_aot"] == pytest.approx(0.1)
    assert diagnostics["surface_cost_curve_min_at_edge"] is True
    assert diagnostics["surface_cost_curve_curvature"] is None


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


class TestIntegratedCostFieldAod:
    def _arrays(self):
        axis = np.array([0.1, 0.2, 0.3], dtype=np.float64)
        cube = np.ones((3, 5, 5), dtype=np.float64)
        cube[1] = 0.0
        prior = np.full((5, 5), 0.2, dtype=np.float64)
        prior_unc = np.full((5, 5), 1e6, dtype=np.float64)
        valid = np.ones((5, 5), dtype=bool)
        x = np.arange(5, dtype=np.float64) * 60.0
        y = np.arange(5, dtype=np.float64) * 60.0
        return axis, cube, prior, prior_unc, valid, x, y

    def test_resolves_window_median_from_cost_cube(self) -> None:
        from tools.aeronet_validation.cost_field_analysis import integrated_cost_field_aod

        axis, cube, prior, prior_unc, valid, x, y = self._arrays()
        band_cost = np.ones((2, 3, 5, 5), dtype=np.float64)
        band_cost[0, 1] = 0.0
        band_cost[1, 2] = 0.0
        band_residual = np.full_like(band_cost, 0.5)
        band_residual[0, 1] = 0.02
        band_residual[1, 2] = 0.03
        out = integrated_cost_field_aod(
            cube=cube,
            band_cost_cube=band_cost,
            band_residual_cube=band_residual,
            band_names=["B02", "B04"],
            aot_axis=axis,
            aot_prior=prior,
            aot_prior_unc=prior_unc,
            solve_valid=valid,
            x=x,
            y=y,
            center_x=120.0,
            center_y=120.0,
            radius_m=90.0,
            pool_window=1,
            min_count=1,
        )

        assert out["aod"] == pytest.approx(0.2)
        assert out["selected_pass"] == "main"
        assert out["main"]["window_shape"] == [3, 3]
        assert out["main"]["n_finite"] == 9
        diagnostics = out["diagnostics"]
        assert diagnostics["local_cost_curve_min_aot"] == pytest.approx(0.2)
        assert diagnostics["local_cost_curve_relative_second_delta"] > 0.0
        assert diagnostics["local_band_B02_argmin_aot"] == pytest.approx(0.2)
        assert diagnostics["local_band_B04_argmin_aot"] == pytest.approx(0.3)
        assert diagnostics["local_band_argmin_spread"] == pytest.approx(0.1)
        assert diagnostics["local_band_B02_residual_final_node"] == pytest.approx(0.02)

    def test_auto2_uses_site_level_abs_gate(self) -> None:
        from tools.aeronet_validation.cost_field_analysis import integrated_cost_field_aod

        axis, cube, prior, prior_unc, valid, x, y = self._arrays()
        cube_abs = np.ones_like(cube)
        cube_abs[0] = 0.0

        clean_tail = integrated_cost_field_aod(
            cube=cube,
            cube_abs=cube_abs,
            aot_axis=axis,
            aot_prior=prior,
            aot_prior_unc=prior_unc,
            solve_valid=valid,
            x=x,
            y=y,
            center_x=120.0,
            center_y=120.0,
            radius_m=90.0,
            pool_window=1,
            min_count=1,
            clean_threshold=0.15,
            high_threshold=0.6,
        )
        assert clean_tail["mode"] == "auto2"
        assert clean_tail["selected_pass"] == "abs"
        assert clean_tail["aod"] == pytest.approx(0.1)

        cube_abs[:] = 1.0
        cube_abs[2] = 0.0
        moderate = integrated_cost_field_aod(
            cube=cube,
            cube_abs=cube_abs,
            aot_axis=axis,
            aot_prior=prior,
            aot_prior_unc=prior_unc,
            solve_valid=valid,
            x=x,
            y=y,
            center_x=120.0,
            center_y=120.0,
            radius_m=90.0,
            pool_window=1,
            min_count=1,
            clean_threshold=0.15,
            high_threshold=0.6,
        )
        assert moderate["selected_pass"] == "shape"
        assert moderate["aod"] == pytest.approx(0.2)

    def test_loads_solver_dump_npz(self, tmp_path) -> None:  # noqa: ANN001
        from tools.aeronet_validation.cost_field_analysis import integrated_cost_field_aod_from_npz

        axis, cube, prior, prior_unc, valid, x, y = self._arrays()
        band_cost = np.ones((2, 3, 5, 5), dtype=np.float64)
        band_cost[0, 1] = 0.0
        band_cost[1, 2] = 0.0
        band_residual = np.full_like(band_cost, 0.5)
        band_residual[0, 1] = 0.02
        band_residual[1, 2] = 0.03
        path = tmp_path / "cost_cube.npz"
        np.savez(
            path,
            cube=cube.astype(np.float32),
            cube_abs=np.zeros(0, dtype=np.float32),
            band_cost_cube=band_cost.astype(np.float32),
            band_residual_cube=band_residual.astype(np.float32),
            band_names=np.asarray(["B02", "B04"]),
            aot_axis=axis.astype(np.float32),
            aot_prior=prior.astype(np.float32),
            aot_prior_unc=prior_unc.astype(np.float32),
            solve_valid=valid,
            x=x,
            y=y,
            pool_window=np.asarray([1]),
            min_count=np.asarray([1]),
        )

        out = integrated_cost_field_aod_from_npz(
            str(path), center_x=120.0, center_y=120.0, radius_m=90.0
        )
        assert out["aod"] == pytest.approx(0.2)
        assert out["diagnostics"]["local_band_argmin_spread"] == pytest.approx(0.1)


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

    coords = {
        "y": np.linspace(40.0, 40.3, shape[0], dtype=np.float64),
        "x": np.linspace(13.0, 13.3, shape[1], dtype=np.float64),
    }
    toa = xr.DataArray(
        np.full((1, *shape), toa_val, dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": ["B02"], **coords},
    )
    cloud_mask = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"], coords=coords)
    surface_prior = SurfacePrior(
        boa=xr.DataArray(
            np.full(shape, prior_boa, dtype=np.float32), dims=["y", "x"], coords=coords
        ),
        boa_unc=xr.DataArray(
            np.full(shape, 0.02, dtype=np.float32), dims=["y", "x"], coords=coords
        ),
        kernels=None,
        mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"], coords=coords),
    )
    geometry = GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"], coords=coords),
        saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"], coords=coords),
        vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"], coords=coords),
        vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"], coords=coords),
    )
    atmo_prior = AtmosphericState(
        aot=xr.DataArray(np.full(shape, 0.15, dtype=np.float32), dims=["y", "x"], coords=coords),
        tcwv=xr.DataArray(np.full(shape, 2.5, dtype=np.float32), dims=["y", "x"], coords=coords),
        tco3=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"], coords=coords),
        aot_unc=xr.DataArray(
            np.full(shape, 0.05, dtype=np.float32), dims=["y", "x"], coords=coords
        ),
        tcwv_unc=xr.DataArray(
            np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"], coords=coords
        ),
        tco3_unc=xr.DataArray(
            np.full(shape, 0.01, dtype=np.float32), dims=["y", "x"], coords=coords
        ),
        elevation=xr.DataArray(
            np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"], coords=coords
        ),
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
        diagnostics = result.diagnostics
        assert diagnostics["surface_cost_mode"] == "chi2"
        assert diagnostics["surface_cost_curve_valid_nodes"] > 0
        assert diagnostics["surface_cost_curve_min_aot"] == pytest.approx(0.20, abs=0.06)
        assert diagnostics["surface_final_aot_median"] == pytest.approx(0.20, abs=0.06)
        assert diagnostics["surface_band_B02_argmin_aot"] == pytest.approx(0.20, abs=0.06)
        assert diagnostics["surface_band_B02_residual_final_node"] < 0.04

    def test_cost_dump_preserves_signed_residuals(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path
    ) -> None:  # noqa: ANN001
        from siac.algorithms.solver import SurfaceDrivenSolver

        dump_path = tmp_path / "cost_cube.npz"
        monkeypatch.setenv("SIAC_DUMP_COST_CUBE", str(dump_path))
        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs(
            toa_val=0.30, prior_boa=0.10
        )
        result = SurfaceDrivenSolver(_solver_config()).solve(
            toa, surface_prior, geometry, atmo_prior, _DummyRT(), cloud_mask, bands
        )

        assert result.success
        with np.load(dump_path) as dumped:
            absolute = np.asarray(dumped["band_residual_cube"])
            signed = np.asarray(dumped["band_signed_residual_cube"])
        assert np.nanmax(signed) > 0.0
        assert np.nanmin(signed) < 0.0
        assert np.allclose(absolute, np.abs(signed), equal_nan=True)

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

    def test_tau_gate_none_is_single_pass(self) -> None:
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
        # No tau_gate configured -> plain single pass, unchanged result.
        solver = SurfaceDrivenSolver(_solver_config())
        result = solver.solve_bundle(bundle)
        assert result.success
        assert result.aot.shape == (4, 4)

    def test_tau_dependent_config_reaches_solve(self) -> None:
        """solve_bundle must not disable a configured tau-dependent prior."""
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
        solver = SurfaceDrivenSolver(_solver_config(surface_driven_tau_dependent_prior=True))
        result = solver.solve_bundle(bundle)
        direct = solver.solve(
            toa,
            surface_prior,
            geometry,
            atmo_prior,
            _DummyRT(),
            cloud_mask,
            bands,
        )

        # solve_bundle is a thin unpack of solve: no gate may reinterpret the
        # configured tau-dependent prior on the way through.
        assert result.success
        assert (
            result.diagnostics["surface_tau_prior_enabled"]
            == (direct.diagnostics["surface_tau_prior_enabled"])
        )

    def test_cloud_water_and_missing_observation_pixels_remain_unretrieved(self) -> None:
        from siac.algorithms.solver import SurfaceDrivenSolver

        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs()
        cloud_mask.values[0, 1] = True
        toa.values[0, 1, 0] = np.nan
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
        assert np.isnan(result.aot.values[0, 0])  # water
        assert np.isnan(result.aot.values[0, 1])  # cloud
        assert np.isnan(result.aot.values[1, 0])  # missing TOA
        assert np.isnan(result.aot_unc.values[0, 0])
        assert abs(float(result.aot.values[3, 3]) - 0.20) < 0.06
        assert result.diagnostics["surface_valid_observation_count"] == 13
        assert result.diagnostics["surface_solved_pixel_count"] == 13

    def test_masked_pixels_are_removed_before_spatial_pooling(self, monkeypatch) -> None:  # noqa: ANN001
        import siac.algorithms.solver.surface_driven as sd
        from siac.algorithms.solver import SurfaceDrivenSolver

        captured: dict[str, np.ndarray] = {}

        def capture_pool(cube, axis, prior, prior_unc, valid, pool_window, min_count):  # noqa: ANN001
            del axis, prior, prior_unc, pool_window, min_count
            captured["cube"] = np.asarray(cube).copy()
            captured["valid"] = np.asarray(valid).copy()
            aot = np.where(valid, 0.2, np.nan)
            unc = np.where(valid, 0.03, np.nan)
            cost = np.where(valid, 1.0, np.nan)
            return aot, unc, cost

        monkeypatch.setattr(sd, "surface_driven_pool_argmin", capture_pool)
        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs()
        cloud_mask.values[0, 1] = True
        water = xr.DataArray(np.zeros((4, 4), dtype=bool), dims=["y", "x"])
        water.values[0, 0] = True
        toa.values[0, 1, 0] = np.nan

        SurfaceDrivenSolver(_solver_config()).solve(
            toa,
            surface_prior,
            geometry,
            atmo_prior,
            _DummyRT(),
            cloud_mask,
            bands,
            water_mask=water,
        )

        invalid = ~captured["valid"]
        assert np.isnan(captured["cube"][:, invalid]).all()
        assert np.isfinite(captured["cube"][:, ~invalid]).all()

    def test_fully_cloudy_scene_returns_no_retrieval(self) -> None:
        from siac.algorithms.solver import SurfaceDrivenSolver

        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs()
        cloud_mask.values[:] = True
        result = SurfaceDrivenSolver(_solver_config()).solve(
            toa, surface_prior, geometry, atmo_prior, _DummyRT(), cloud_mask, bands
        )

        assert not result.success
        assert np.isnan(result.aot.values).all()
        assert np.isnan(result.aot_unc.values).all()
        assert result.diagnostics["surface_valid_observation_count"] == 0
        assert result.diagnostics["surface_solved_pixel_count"] == 0

    def test_cloud_rescue_keeps_water_masked(self) -> None:
        from siac.algorithms.solver import SurfaceDrivenSolver

        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs()
        cloud_mask.values[:] = True
        water = np.zeros((4, 4), dtype=bool)
        water[0, 0] = True
        water_mask = xr.DataArray(water, dims=["y", "x"])
        solver = SurfaceDrivenSolver(_solver_config(surface_driven_allow_cloud_retrieval=True))
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
        assert np.isnan(result.aot.values[0, 0])
        assert np.isfinite(result.aot.values[1:, :]).all()
        assert result.diagnostics["surface_cloud_mask_bypassed"] is True
        assert result.diagnostics["surface_water_mask_bypassed"] is False

    def test_calibrated_backstop_branch(self) -> None:
        from siac.algorithms.solver import SurfaceDrivenSolver

        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs()
        solver = SurfaceDrivenSolver(_solver_config(surface_driven_backstop_calibrated=True))
        result = solver.solve(
            toa, surface_prior, geometry, atmo_prior, _DummyRT(), cloud_mask, bands
        )
        assert result.success and np.all(np.isfinite(result.aot.values))

    def test_species_mode_requires_native_sixs_backend(self) -> None:
        from siac.algorithms.solver import SurfaceDrivenSolver

        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs()
        solver = SurfaceDrivenSolver(
            _solver_config(surface_driven_aerosol_species="cci_climatology_exact")
        )

        with pytest.raises(ValueError, match='requires algorithms.rt.backend="sixs"'):
            solver.solve(toa, surface_prior, geometry, atmo_prior, _DummyRT(), cloud_mask, bands)

    def test_reference_tcwv_routes_into_cost_cube_rt(self) -> None:
        """The configured reference TCWV feeds the cost-cube RT; None keeps the
        real scene-median (2.5 here) so existing behaviour is unchanged."""
        from siac.algorithms.solver import SurfaceDrivenSolver

        class _CapturingRT(_DummyRT):
            def __init__(self) -> None:
                self.seen_tcwv: list[float] = []

            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
                self.seen_tcwv.append(float(np.asarray(atmo_state.tcwv.values).reshape(-1)[0]))
                return super().compute_coefficients(geometry, atmo_state, band, compute_jacobian)

        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs()

        rt_default = _CapturingRT()
        SurfaceDrivenSolver(_solver_config()).solve(
            toa, surface_prior, geometry, atmo_prior, rt_default, cloud_mask, bands
        )
        assert rt_default.seen_tcwv  # provider was exercised
        assert all(v == pytest.approx(2.5) for v in rt_default.seen_tcwv)

        rt_ref = _CapturingRT()
        SurfaceDrivenSolver(_solver_config(surface_driven_reference_tcwv=2.0)).solve(
            toa, surface_prior, geometry, atmo_prior, rt_ref, cloud_mask, bands
        )
        assert rt_ref.seen_tcwv
        assert all(v == pytest.approx(2.0) for v in rt_ref.seen_tcwv)

    def test_fixed_spatial_tcwv_reaches_cost_cube_rt_without_joint_lut(self) -> None:
        """A spatial WVP field stays fixed per pixel while only AOT is swept."""
        from dataclasses import replace

        from siac.algorithms.solver import SurfaceDrivenSolver

        class _SpatialTCWVRT(_DummyRT):
            def __init__(self) -> None:
                self.tcwv_spreads: list[float] = []
                self.joint_lut_calls = 0

            def build_joint_grid_search_lut(self, **_kwargs):  # noqa: ANN003
                self.joint_lut_calls += 1
                raise AssertionError("spatial fixed TCWV must not build a scalar joint LUT")

            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
                tcwv = np.asarray(atmo_state.tcwv.values, dtype=np.float64)
                self.tcwv_spreads.append(float(np.nanmax(tcwv) - np.nanmin(tcwv)))
                return super().compute_coefficients(geometry, atmo_state, band, compute_jacobian)

        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs()
        spatial_tcwv = xr.DataArray(
            np.linspace(1.0, 3.0, 16, dtype=np.float32).reshape(4, 4),
            dims=("y", "x"),
            coords=atmo_prior.tcwv.coords,
        )
        atmo_prior = replace(atmo_prior, tcwv=spatial_tcwv)
        rt = _SpatialTCWVRT()

        result = SurfaceDrivenSolver(_solver_config()).solve(
            toa, surface_prior, geometry, atmo_prior, rt, cloud_mask, bands
        )

        assert result.success
        assert rt.joint_lut_calls == 0
        assert rt.tcwv_spreads and all(spread == pytest.approx(2.0) for spread in rt.tcwv_spreads)

    def test_scene_mean_geometry_routes_into_cost_cube_rt(self) -> None:
        """With scene_mean_geometry the cost-cube RT sees a single (mean) SZA for
        every call; default (False) passes the per-pixel field through unchanged."""
        from siac.algorithms.solver import SurfaceDrivenSolver

        class _GeomRT(_DummyRT):
            def __init__(self) -> None:
                self.seen_sza_spread: list[float] = []

            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
                sza = np.asarray(geometry.sza.values)
                self.seen_sza_spread.append(float(np.nanmax(sza) - np.nanmin(sza)))
                return super().compute_coefficients(geometry, atmo_state, band, compute_jacobian)

        toa, surface_prior, _geom, atmo_prior, cloud_mask, bands = _make_inputs()
        # Per-pixel geometry with a real SZA spread across the grid.
        ny, nx = cloud_mask.sizes["y"], cloud_mask.sizes["x"]
        sza = np.linspace(0.4, 0.6, ny * nx, dtype=np.float64).reshape(ny, nx)
        from siac.runtime import GeometryAngles

        geometry = GeometryAngles(
            sza=xr.DataArray(sza, dims=["y", "x"]),
            saa=xr.DataArray(np.full((ny, nx), 2.5), dims=["y", "x"]),
            vza=xr.DataArray(np.full((ny, nx), 0.1), dims=["y", "x"]),
            vaa=xr.DataArray(np.full((ny, nx), 1.5), dims=["y", "x"]),
        )

        rt_default = _GeomRT()
        SurfaceDrivenSolver(_solver_config()).solve(
            toa, surface_prior, geometry, atmo_prior, rt_default, cloud_mask, bands
        )
        assert rt_default.seen_sza_spread
        assert max(rt_default.seen_sza_spread) > 0.1  # per-pixel spread preserved

        rt_mean = _GeomRT()
        SurfaceDrivenSolver(_solver_config(surface_driven_scene_mean_geometry=True)).solve(
            toa, surface_prior, geometry, atmo_prior, rt_mean, cloud_mask, bands
        )
        assert rt_mean.seen_sza_spread
        assert all(s == pytest.approx(0.0) for s in rt_mean.seen_sza_spread)  # collapsed to mean


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
