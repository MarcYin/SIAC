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
        out = integrated_cost_field_aod(
            cube=cube,
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

    def test_loads_solver_dump_npz(self, tmp_path) -> None:  # noqa: ANN001
        from siac.algorithms.solver.surface_driven import integrated_cost_field_aod_from_npz

        axis, cube, prior, prior_unc, valid, x, y = self._arrays()
        path = tmp_path / "cost_cube.npz"
        np.savez(
            path,
            cube=cube.astype(np.float32),
            cube_abs=np.zeros(0, dtype=np.float32),
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
