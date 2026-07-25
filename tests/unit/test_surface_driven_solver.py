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


def test_profile_scale_cost_absorbs_common_brightness_mismatch() -> None:
    from siac.algorithms.solver.surface_driven import _profile_scale_cost_node

    ref = np.array([0.10, 0.20, 0.30], dtype=np.float64)[:, None, None]
    boa = 1.08 * ref
    uncertainty = np.full_like(ref, 0.01)

    cost = _profile_scale_cost_node(
        boa,
        ref,
        uncertainty,
        [0, 1, 2],
        scale_sigma=0.2,
    )
    absolute = np.sum(np.square((boa - ref) / uncertainty), axis=0)

    assert cost.shape == (1, 1)
    assert cost[0, 0] < absolute[0, 0]


def test_additive_offset_cost_absorbs_common_reflectance_mismatch() -> None:
    from siac.algorithms.solver.surface_driven import _additive_offset_cost_node

    ref = np.array([0.10, 0.20, 0.30], dtype=np.float64)[:, None, None]
    uncertainty = np.full_like(ref, 0.01)
    common_offset = ref + 0.03
    spectral_mismatch = common_offset.copy()
    spectral_mismatch[0] += 0.04

    common_cost = _additive_offset_cost_node(
        common_offset,
        ref,
        uncertainty,
        [0, 1, 2],
    )
    spectral_cost = _additive_offset_cost_node(
        spectral_mismatch,
        ref,
        uncertainty,
        [0, 1, 2],
    )

    assert common_cost[0, 0] == pytest.approx(0.0, abs=1e-12)
    assert spectral_cost[0, 0] > common_cost[0, 0]


def test_additive_offset_cost_requires_two_bands() -> None:
    from siac.algorithms.solver.surface_driven import _additive_offset_cost_node

    values = np.ones((1, 1, 1), dtype=np.float64)
    with pytest.raises(ValueError, match="at least two bands"):
        _additive_offset_cost_node(values, values, values, [0])


def test_loo_scale_cost_penalizes_an_unpredicted_band() -> None:
    from siac.algorithms.solver.surface_driven import _loo_scale_cost_node

    ref = np.array([0.10, 0.20, 0.30], dtype=np.float64)[:, None, None]
    uncertainty = np.full_like(ref, 0.01)
    consistent = 1.05 * ref
    inconsistent = consistent.copy()
    inconsistent[1] += 0.05

    consistent_cost = _loo_scale_cost_node(
        consistent,
        ref,
        uncertainty,
        [0, 1, 2],
        scale_sigma=0.2,
    )
    inconsistent_cost = _loo_scale_cost_node(
        inconsistent,
        ref,
        uncertainty,
        [0, 1, 2],
        scale_sigma=0.2,
    )

    assert inconsistent_cost[0, 0] > consistent_cost[0, 0]


def test_loo_scale_cost_requires_three_bands() -> None:
    from siac.algorithms.solver.surface_driven import _loo_scale_cost_node

    values = np.ones((2, 1, 1), dtype=np.float64)
    with pytest.raises(ValueError, match="at least three bands"):
        _loo_scale_cost_node(
            values,
            values,
            values,
            [0, 1],
            scale_sigma=0.1,
        )


def test_trimmed_chi2_drops_one_inconsistent_band() -> None:
    from siac.algorithms.solver.surface_driven import _trimmed_chi2_cost_node

    ref = np.array([0.10, 0.20, 0.30, 0.40], dtype=np.float64)[:, None, None]
    boa = ref + 0.01
    boa[2] += 0.20
    inv_var = np.full_like(ref, 10000.0)

    cost = _trimmed_chi2_cost_node(boa, ref, inv_var, [0, 1, 2, 3])

    assert cost[0, 0] == pytest.approx(3.0)


def test_trimmed_chi2_requires_three_bands() -> None:
    from siac.algorithms.solver.surface_driven import _trimmed_chi2_cost_node

    values = np.ones((2, 1, 1), dtype=np.float64)
    with pytest.raises(ValueError, match="at least three bands"):
        _trimmed_chi2_cost_node(values, values, values, [0, 1])


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
        from siac.algorithms.solver.surface_driven import integrated_cost_field_aod

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
        from siac.algorithms.solver.surface_driven import integrated_cost_field_aod

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

    def test_auto2_prefers_shape_when_abs_tail_not_materially_better(self) -> None:
        from siac.algorithms.solver.surface_driven import _select_auto2_solution

        selected_aod, selected_unc, selected_cost, diagnostics = _select_auto2_solution(
            aot_main=np.array([0.2, 0.4, 1.8], dtype=np.float64),
            unc_main=np.array([0.01, 0.01, 0.01], dtype=np.float64),
            cost_main=np.array([1.0, 1.0, 10.0], dtype=np.float64),
            aot_abs=np.array([0.05, 0.4, 2.2], dtype=np.float64),
            unc_abs=np.array([0.01, 0.01, 0.01], dtype=np.float64),
            cost_abs=np.array([0.70, 0.90, 10.0], dtype=np.float64),
            clean_threshold=0.15,
            high_threshold=0.6,
            cost_gain=0.2,
        )
        assert selected_aod[0] == pytest.approx(0.05)
        assert selected_aod[1] == pytest.approx(0.4)
        assert selected_aod[2] == pytest.approx(1.8)
        assert selected_cost[0] == pytest.approx(0.7)
        assert selected_cost[1] == pytest.approx(1.0)
        assert selected_cost[2] == pytest.approx(10.0)
        assert diagnostics["surface_auto2_abs_selected_pixels"] == 1
        assert diagnostics["surface_auto2_shape_selected_pixels"] == 2

    def test_loads_solver_dump_npz(self, tmp_path) -> None:  # noqa: ANN001
        from siac.algorithms.solver.surface_driven import integrated_cost_field_aod_from_npz

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

    def test_tau_dependent_config_reaches_solve_without_gate(self, monkeypatch) -> None:  # noqa: ANN001
        """Pure tau-dependent mode must not be disabled by solve_bundle."""
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
        original_solve = solver.solve
        force_values = []

        def wrapped_solve(*args, **kwargs):  # noqa: ANN001
            force_values.append(kwargs.get("_force_tau_prior", "missing"))
            return original_solve(*args, **kwargs)

        monkeypatch.setattr(solver, "solve", wrapped_solve)
        result = solver.solve_bundle(bundle)

        assert result.success
        assert force_values == [None]

    def test_tau_gate_forces_second_pass_when_static_cost_exceeds_gate(self, monkeypatch) -> None:  # noqa: ANN001
        from dataclasses import replace

        from siac.algorithms.solver import SurfaceDrivenSolver
        from siac.runtime import SolverInputBundle

        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs()
        surface_prior = replace(surface_prior, tau_predictor={"trees": []})
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
        solver = SurfaceDrivenSolver(_solver_config(surface_driven_tau_gate_cost=1.0e-12))
        original_solve = solver.solve
        force_values = []

        def wrapped_solve(*args, **kwargs):  # noqa: ANN001
            force_values.append(kwargs.get("_force_tau_prior", "missing"))
            return original_solve(*args, **kwargs)

        monkeypatch.setattr(solver, "solve", wrapped_solve)
        result = solver.solve_bundle(bundle)

        assert result.success
        assert force_values == [False, True]
        assert result.diagnostics["surface_tau_gate_configured"] is True
        assert result.diagnostics["surface_tau_available"] is True
        assert result.diagnostics["surface_tau_gate_fired"] is True
        assert result.diagnostics["surface_static_cost_per_band"] > 0.0
        assert result.diagnostics["surface_tau_cost_per_band"] > 0.0

    def test_tau_gate_skips_second_pass_when_static_fits(self) -> None:
        # A tau_gate is set but the static solve fits well (low cost) and there
        # is no tau_predictor payload, so the gate is a no-op and the result
        # matches the ungated solve.
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
        gated = SurfaceDrivenSolver(_solver_config(surface_driven_tau_gate_cost=5.0))
        plain = SurfaceDrivenSolver(_solver_config())
        r_gated = gated.solve_bundle(bundle)
        r_plain = plain.solve_bundle(bundle)
        assert np.allclose(r_gated.aot.values, r_plain.aot.values)
        assert r_gated.diagnostics["surface_tau_gate_configured"] is True
        assert r_gated.diagnostics["surface_tau_available"] is False
        assert r_gated.diagnostics["surface_tau_gate_fired"] is False

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

    def test_backstop_uncertainty_scale_multiplies_selected_width(self) -> None:
        from siac.algorithms.solver import SurfaceDrivenSolver

        prior = np.array([[0.4]], dtype=np.float32)
        raw_unc = np.array([[0.1]], dtype=np.float32)
        base = SurfaceDrivenSolver(
            _solver_config(
                surface_driven_backstop_calibrated=False,
                surface_driven_backstop_uncertainty_scale=1.0,
            )
        )._backstop_unc(prior, raw_unc)
        widened = SurfaceDrivenSolver(
            _solver_config(
                surface_driven_backstop_calibrated=False,
                surface_driven_backstop_uncertainty_scale=3.0,
            )
        )._backstop_unc(prior, raw_unc)

        assert widened == pytest.approx(3.0 * base)

    def test_species_mode_requires_native_sixs_backend(self) -> None:
        from siac.algorithms.solver import SurfaceDrivenSolver

        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs()
        solver = SurfaceDrivenSolver(
            _solver_config(surface_driven_aerosol_species="cci_climatology")
        )

        with pytest.raises(ValueError, match='requires algorithms.rt.backend="sixs"'):
            solver.solve(toa, surface_prior, geometry, atmo_prior, _DummyRT(), cloud_mask, bands)

    def test_species_mode_clones_three_candidates_and_selects_lowest_cost(self) -> None:
        from types import SimpleNamespace

        from siac.algorithms.solver import SurfaceDrivenSolver
        from siac.config import RTSetupConfig
        from siac.runtime.models import RTCoefficients

        class _SpeciesRT:
            def __init__(self, bias: float = 0.0, clones: list[RTSetupConfig] | None = None):
                self.bias = bias
                self.clones = clones if clones is not None else []
                self._config = SimpleNamespace(month=6)
                self.observation_time = None
                self.rt_setup = RTSetupConfig()

            @property
            def backend_name(self) -> str:
                return "sixs"

            def with_rt_setup(self, rt_setup: RTSetupConfig) -> _SpeciesRT:
                biases = [0.5, 0.0, 0.5]
                bias = biases[len(self.clones)]
                self.clones.append(rt_setup)
                return _SpeciesRT(bias=bias, clones=self.clones)

            def compute_coefficients(
                self,
                geometry,  # noqa: ANN001
                atmo_state,  # noqa: ANN001
                band,  # noqa: ANN001
                compute_jacobian: bool = False,
            ) -> RTCoefficients:
                del geometry, band, compute_jacobian
                template = atmo_state.aot
                xbp = template + np.float32(self.bias)
                return RTCoefficients(
                    xap=xr.ones_like(template),
                    xbp=xbp,
                    xcp=xr.zeros_like(template),
                )

        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs()
        rt = _SpeciesRT()
        solver = SurfaceDrivenSolver(
            _solver_config(
                surface_driven_aerosol_species="cci_climatology",
                surface_driven_aerosol_species_candidates=3,
            )
        )
        result = solver.solve(toa, surface_prior, geometry, atmo_prior, rt, cloud_mask, bands)

        assert result.success
        assert len(rt.clones) == 3
        assert all(clone.aerosol is not None for clone in rt.clones)
        assert all(clone.aerosol.profile == "multimodal_log_normal" for clone in rt.clones)
        assert abs(float(np.median(result.aot.values)) - 0.20) < 0.06
        assert result.diagnostics["surface_rt_branch"] == "tile_species_min_median_cost"
        assert result.diagnostics["surface_species_selected_candidate"] == 1

        exact_rt = _SpeciesRT()
        exact_solver = SurfaceDrivenSolver(
            _solver_config(surface_driven_aerosol_species="cci_climatology_exact")
        )
        exact_result = exact_solver.solve(
            toa,
            surface_prior,
            geometry,
            atmo_prior,
            exact_rt,
            cloud_mask,
            bands,
        )

        assert exact_result.success
        assert len(exact_rt.clones) == 1
        assert exact_rt.clones[0].aerosol is not None
        assert exact_rt.clones[0].aerosol.profile == "multimodal_log_normal"

    def test_canonical_species_mode_clones_standard_sixs_basis(self) -> None:
        from types import SimpleNamespace

        from siac.algorithms.solver.surface_driven import _species_candidate_rt_models
        from siac.config import RTSetupConfig

        class _CanonicalRT:
            def __init__(self, clones: list[RTSetupConfig] | None = None):
                self.clones = clones if clones is not None else []
                self._config = SimpleNamespace(month=6)
                self.observation_time = None
                self.rt_setup = RTSetupConfig()

            @property
            def backend_name(self) -> str:
                return "sixs"

            def with_rt_setup(self, rt_setup: RTSetupConfig) -> _CanonicalRT:
                self.clones.append(rt_setup)
                return _CanonicalRT(clones=self.clones)

        _, _, _, atmo_prior, _, _ = _make_inputs()
        rt = _CanonicalRT()
        candidates = _species_candidate_rt_models(
            rt_model=rt,
            config=_solver_config(surface_driven_aerosol_species="canonical_6s"),
            template=atmo_prior.aot,
        )

        assert len(candidates) == 5
        assert [str(clone.aerosol.profile) for clone in rt.clones] == [
            "continental",
            "maritime",
            "urban",
            "desert",
            "biomass_burning",
        ]

    def test_species_selection_uses_one_candidate_for_the_whole_tile(self, monkeypatch) -> None:  # noqa: ANN001
        import siac.algorithms.solver.surface_driven as sd
        from siac.algorithms.solver import SurfaceDrivenSolver

        monkeypatch.setattr(
            sd,
            "_species_candidate_rt_models",
            lambda **_kwargs: [_DummyRT(), _DummyRT(), _DummyRT()],
        )
        calls = 0

        def fake_pool_argmin(cube, axis, prior, prior_unc, valid, pool_window, min_count):  # noqa: ANN001
            nonlocal calls
            del cube, axis, prior, prior_unc, pool_window, min_count
            candidate_index = calls
            calls += 1
            aot = np.where(valid, 0.1 * (candidate_index + 1), np.nan)
            unc = np.where(valid, 0.03, np.nan)
            cost = np.where(valid, float(candidate_index + 1), np.nan)
            if candidate_index == 1:
                # Candidate 1 wins one pixel but loses the robust tile median.
                cost[0, 0] = 0.1
            return aot, unc, cost

        monkeypatch.setattr(sd, "surface_driven_pool_argmin", fake_pool_argmin)
        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs()
        result = SurfaceDrivenSolver(_solver_config()).solve(
            toa,
            surface_prior,
            geometry,
            atmo_prior,
            _DummyRT(),
            cloud_mask,
            bands,
        )

        assert result.success
        assert calls == 3
        assert np.allclose(result.aot.values, 0.1)
        assert result.diagnostics["surface_species_selected_candidate"] == 0
        assert result.diagnostics["surface_species_candidate_median_costs"] == [1.0, 2.0, 3.0]
        assert result.diagnostics["surface_rt_branch"] == "tile_species_min_median_cost"

    def test_backstop_escape_replaces_rail_low_aod(self, monkeypatch) -> None:  # noqa: ANN001
        import siac.algorithms.solver.surface_driven as sd
        from siac.algorithms.solver import SurfaceDrivenSolver

        def fake_cost_curve_diagnostics(cube, aot_axis):  # noqa: ANN001
            axis = np.asarray(aot_axis, dtype=np.float64)
            return {
                "surface_cost_curve_valid_nodes": int(axis.size),
                "surface_cost_curve_curvature": None,
                "surface_cost_curve_min_at_edge": True,
                "surface_cost_curve_min_aot": float(axis[0]),
                "surface_cost_curve_min_cost": 0.20,
                "surface_cost_curve_second_aot": float(axis[1]),
                "surface_cost_curve_second_delta": 0.001,
                "surface_cost_curve_relative_second_delta": 0.001,
            }

        def fake_pool_argmin(cube, axis, prior, prior_unc, valid, pool_window, min_count):  # noqa: ANN001
            del pool_window, min_count
            target_aot = 0.01
            if np.all(np.isinf(prior_unc)):
                target_aot = 3.0
            aod = np.where(valid, target_aot, np.nan)
            unc = np.where(valid, 0.03, np.nan)
            cost = np.where(valid, 10.0 if target_aot > 0.5 else 30.0, np.nan)
            return aod.astype(np.float64), unc.astype(np.float64), cost.astype(np.float64)

        monkeypatch.setattr(sd, "_cost_curve_diagnostics", fake_cost_curve_diagnostics)
        monkeypatch.setattr(sd, "surface_driven_pool_argmin", fake_pool_argmin)

        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs()
        solver = SurfaceDrivenSolver(
            _solver_config(
                surface_driven_backstop_escape_enabled=True,
                surface_driven_pool_radius_m=0.0,
                surface_driven_aot_axis="acixthree",
            )
        )
        result = solver.solve(
            toa, surface_prior, geometry, atmo_prior, _DummyRT(), cloud_mask, bands
        )
        assert result.success
        assert float(np.median(result.aot.values[~np.isnan(result.aot.values)])) == pytest.approx(
            3.0, abs=1e-6
        )

    def test_backstop_escape_respects_disabled_switch(self, monkeypatch) -> None:  # noqa: ANN001
        import siac.algorithms.solver.surface_driven as sd
        from siac.algorithms.solver import SurfaceDrivenSolver

        def fake_cost_curve_diagnostics(cube, aot_axis):  # noqa: ANN001
            axis = np.asarray(aot_axis, dtype=np.float64)
            return {
                "surface_cost_curve_valid_nodes": int(axis.size),
                "surface_cost_curve_curvature": None,
                "surface_cost_curve_min_at_edge": True,
                "surface_cost_curve_min_aot": float(axis[0]),
                "surface_cost_curve_min_cost": 0.20,
                "surface_cost_curve_second_aot": float(axis[1]),
                "surface_cost_curve_second_delta": 0.001,
                "surface_cost_curve_relative_second_delta": 0.001,
            }

        def fake_pool_argmin(cube, axis, prior, prior_unc, valid, pool_window, min_count):  # noqa: ANN001
            del pool_window, min_count, cube, axis, prior
            aod = np.where(valid, 0.01, np.nan)
            unc = np.where(valid, 0.08, np.nan)
            cost = np.where(valid, 30.0, np.nan)
            return aod.astype(np.float64), unc.astype(np.float64), cost.astype(np.float64)

        monkeypatch.setattr(sd, "_cost_curve_diagnostics", fake_cost_curve_diagnostics)
        monkeypatch.setattr(sd, "surface_driven_pool_argmin", fake_pool_argmin)

        toa, surface_prior, geometry, atmo_prior, cloud_mask, bands = _make_inputs()
        solver = SurfaceDrivenSolver(
            _solver_config(
                surface_driven_backstop_escape_enabled=False,
                surface_driven_pool_radius_m=0.0,
                surface_driven_aot_axis="acixthree",
            )
        )
        result = solver.solve(
            toa, surface_prior, geometry, atmo_prior, _DummyRT(), cloud_mask, bands
        )
        assert result.success
        assert float(np.median(result.aot.values[~np.isnan(result.aot.values)])) == pytest.approx(
            0.01, abs=1e-6
        )

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
