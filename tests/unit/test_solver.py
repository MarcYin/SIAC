"""
Unit tests for solver components.
"""

import numpy as np
import pytest
import xarray as xr

import siac._rust as native_rust
from siac._rust import (
    evaluate_grid_search_candidate_cost,
    evaluate_grid_search_cost_cube_with_provider,
    evaluate_grid_search_cost_cube_with_provider_qa,
    quadratic_refine_grid_search,
    quadratic_refine_grid_search_qa,
)
from siac._rust_compat import evaluate_block_grid_search_cost_cube_with_provider_qa
from siac.algorithms.solver.cost import (
    CostFunction,
    CostFunctionConfig,
    apply_smoothness_filter,
    compute_laplacian_eigenvalues,
)
from siac.algorithms.solver.multigrid import (
    MultiGridConfig,
    MultiGridSolver,
    SolverResult,
    SolverStageConfig,
    StagedMultiGridSolver,
    build_solver_valid_mask,
)
from siac.runtime import AtmosphericState, SurfacePrior

HAS_NATIVE_BLOCK_COST_CUBE = hasattr(
    native_rust, "evaluate_block_grid_search_cost_cube_with_provider_qa"
)


def _quadratic_refine_python_reference(
    costs: np.ndarray,
    aot_axis: np.ndarray,
    tcwv_axis: np.ndarray,
    valid_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Reference implementation mirroring Rust local quadratic refinement."""
    n_aot, n_tcwv, ny, nx = costs.shape
    aot_best = np.zeros((ny, nx), dtype=np.float32)
    tcwv_best = np.zeros((ny, nx), dtype=np.float32)
    da = max(float(abs(aot_axis[1] - aot_axis[0])), 1e-6) if n_aot > 1 else 0.05
    dt = max(float(abs(tcwv_axis[1] - tcwv_axis[0])), 1e-6) if n_tcwv > 1 else 0.2
    aot_unc = np.full((ny, nx), da, dtype=np.float32)
    tcwv_unc = np.full((ny, nx), dt, dtype=np.float32)

    for iy in range(ny):
        for ix in range(nx):
            if not bool(valid_mask[iy, ix]):
                continue

            best_ia = 0
            best_it = 0
            best_cost = np.inf
            for ia in range(n_aot):
                for it in range(n_tcwv):
                    c = float(costs[ia, it, iy, ix])
                    if np.isfinite(c) and c < best_cost:
                        best_cost = c
                        best_ia = ia
                        best_it = it

            aot_best[iy, ix] = float(aot_axis[best_ia])
            tcwv_best[iy, ix] = float(tcwv_axis[best_it])

            if best_ia == 0 or best_ia + 1 >= n_aot or best_it == 0 or best_it + 1 >= n_tcwv:
                continue

            ia = best_ia
            it = best_it
            f00 = float(costs[ia, it, iy, ix])
            fxm = float(costs[ia - 1, it, iy, ix])
            fxp = float(costs[ia + 1, it, iy, ix])
            fym = float(costs[ia, it - 1, iy, ix])
            fyp = float(costs[ia, it + 1, iy, ix])
            fmm = float(costs[ia - 1, it - 1, iy, ix])
            fmp = float(costs[ia - 1, it + 1, iy, ix])
            fpm = float(costs[ia + 1, it - 1, iy, ix])
            fpp = float(costs[ia + 1, it + 1, iy, ix])
            if not np.isfinite([f00, fxm, fxp, fym, fyp, fmm, fmp, fpm, fpp]).all():
                continue

            dfdx = (fxp - fxm) / (2.0 * da)
            dfdy = (fyp - fym) / (2.0 * dt)
            d2fdx2 = (fxp - 2.0 * f00 + fxm) / (da * da)
            d2fdy2 = (fyp - 2.0 * f00 + fym) / (dt * dt)
            d2fdxdy = (fpp - fpm - fmp + fmm) / (4.0 * da * dt)

            a = d2fdx2
            b = d2fdxdy
            c = d2fdy2
            det = a * c - b * b
            if not (a > 1e-12 and c > 1e-12 and det > 1e-12):
                continue

            delta_a = (c * dfdx - b * dfdy) / det
            delta_t = (-b * dfdx + a * dfdy) / det
            a_fit = float(aot_axis[ia]) - delta_a
            t_fit = float(tcwv_axis[it]) - delta_t
            a_lo = float(aot_axis[ia - 1])
            a_hi = float(aot_axis[ia + 1])
            t_lo = float(tcwv_axis[it - 1])
            t_hi = float(tcwv_axis[it + 1])
            aot_best[iy, ix] = np.float32(np.clip(a_fit, min(a_lo, a_hi), max(a_lo, a_hi)))
            tcwv_best[iy, ix] = np.float32(np.clip(t_fit, min(t_lo, t_hi), max(t_lo, t_hi)))
            aot_unc[iy, ix] = np.float32(np.sqrt(max(c / det, 1e-12)))
            tcwv_unc[iy, ix] = np.float32(np.sqrt(max(a / det, 1e-12)))

    return aot_best, tcwv_best, aot_unc, tcwv_unc


def _candidate_cost_python_reference(
    toa: np.ndarray,
    xap: np.ndarray,
    xbp: np.ndarray,
    xcp: np.ndarray,
    boa_prior: np.ndarray,
    boa_unc: np.ndarray,
    band_weights: np.ndarray,
    valid_mask: np.ndarray,
    aot_val: float,
    tcwv_val: float,
    aot_prior: np.ndarray,
    tcwv_prior: np.ndarray,
    aot_prior_unc: np.ndarray,
    tcwv_prior_unc: np.ndarray,
) -> np.ndarray:
    """Reference implementation for candidate cost map."""
    n_band, ny, nx = toa.shape
    out = np.zeros((ny, nx), dtype=np.float32)
    for ib in range(n_band):
        y = xap[ib] * toa[ib] - xbp[ib]
        denom = 1.0 + xcp[ib] * y
        boa_model = np.divide(
            y,
            denom,
            out=np.full_like(y, np.nan, dtype=np.float32),
            where=np.abs(denom) > 1e-12,
        )
        diff = boa_model - boa_prior[ib]
        weight = band_weights[ib] / np.maximum(np.square(boa_unc[ib]), 1e-12)
        valid = valid_mask & np.isfinite(diff) & np.isfinite(weight)
        out += 0.5 * np.where(valid, weight * np.square(diff), 0.0).astype(np.float32)

    prior = 0.5 * (
        np.square(aot_val - aot_prior) / np.maximum(np.square(aot_prior_unc), 1e-12)
        + np.square(tcwv_val - tcwv_prior) / np.maximum(np.square(tcwv_prior_unc), 1e-12)
    )
    out += np.where(valid_mask, prior, 0.0).astype(np.float32)
    return out


def _sum_blocks_reference(values: np.ndarray, block_size: int) -> np.ndarray:
    """Aggregate the last two dimensions by summing over non-overlapping blocks."""
    if block_size <= 1:
        return np.asarray(values)

    source = np.asarray(values)
    by = (source.shape[-2] + block_size - 1) // block_size
    bx = (source.shape[-1] + block_size - 1) // block_size
    out = np.zeros((*source.shape[:-2], by, bx), dtype=source.dtype)
    for iy in range(by):
        y0 = iy * block_size
        y1 = min(y0 + block_size, source.shape[-2])
        for ix in range(bx):
            x0 = ix * block_size
            x1 = min(x0 + block_size, source.shape[-1])
            out[..., iy, ix] = np.sum(source[..., y0:y1, x0:x1], axis=(-2, -1))
    return out


def _any_blocks_reference(mask: np.ndarray, block_size: int) -> np.ndarray:
    """Aggregate a boolean mask by non-overlapping any-valid blocks."""
    if block_size <= 1:
        return np.asarray(mask, dtype=bool)

    source = np.asarray(mask, dtype=bool)
    by = (source.shape[0] + block_size - 1) // block_size
    bx = (source.shape[1] + block_size - 1) // block_size
    out = np.zeros((by, bx), dtype=bool)
    for iy in range(by):
        y0 = iy * block_size
        y1 = min(y0 + block_size, source.shape[0])
        for ix in range(bx):
            x0 = ix * block_size
            x1 = min(x0 + block_size, source.shape[1])
            out[iy, ix] = bool(np.any(source[y0:y1, x0:x1]))
    return out


class TestLaplacianEigenvalues:
    """Tests for Laplacian eigenvalue computation."""

    def test_shape(self):
        """Eigenvalues should have correct shape."""
        lambda_vals = compute_laplacian_eigenvalues(10, 20)
        assert lambda_vals.shape == (20, 10)

    def test_origin(self):
        """Eigenvalue at origin should be zero."""
        lambda_vals = compute_laplacian_eigenvalues(10, 10)
        assert lambda_vals[0, 0] == 0.0

    def test_positive(self):
        """Eigenvalues should be non-negative."""
        lambda_vals = compute_laplacian_eigenvalues(15, 15)
        assert np.all(lambda_vals >= 0)

    def test_symmetry(self):
        """Eigenvalues should have symmetry properties."""
        n = 16
        lambda_vals = compute_laplacian_eigenvalues(n, n)

        # Check i,j symmetry with j,i (separable and symmetric in x/y)
        for i in range(n):
            for j in range(n):
                np.testing.assert_allclose(lambda_vals[i, j], lambda_vals[j, i], rtol=1e-10)


class TestSmoothnessFilter:
    """Tests for edge-preserving smoothness filter application."""

    def test_identity_gamma_zero(self):
        """With gamma=0, output should equal input."""
        x = np.random.default_rng(0).standard_normal((20, 20))

        x_smooth = apply_smoothness_filter(x, gamma=0)

        np.testing.assert_allclose(x_smooth, x, rtol=1e-10)

    def test_smoothing_effect(self):
        """With gamma>0, noise should be suppressed."""
        # Create noisy data (small random noise on a smooth field)
        rng = np.random.default_rng(0)
        x = rng.normal(0.15, 0.01, (32, 32))

        x_smooth = apply_smoothness_filter(x, gamma=5.0, delta=0.02, n_iter=50)

        # Smoothed data should have lower variance
        assert x_smooth.var() < x.var()
        # Mean should be approximately preserved
        np.testing.assert_allclose(x.mean(), x_smooth.mean(), rtol=0.05)

    def test_preserves_mean(self):
        """Mean should be approximately preserved."""
        x = np.random.default_rng(1).standard_normal((20, 20)) + 5.0

        x_smooth = apply_smoothness_filter(x, gamma=2.0)

        np.testing.assert_allclose(x.mean(), x_smooth.mean(), rtol=0.05)

    def test_preserves_hotspot(self):
        """A strong localized feature should survive with large delta."""
        rng = np.random.default_rng(42)
        x = rng.normal(0.15, 0.005, (32, 32))
        # Add a strong hotspot
        x[14:18, 14:18] = 0.8

        x_smooth = apply_smoothness_filter(x, gamma=5.0, delta=0.02)

        # The hotspot peak should be largely preserved (within 30%)
        assert x_smooth[15:17, 15:17].mean() > 0.55
        # Background should still be smooth
        bg = x_smooth[:10, :10]
        assert bg.std() < x[:10, :10].std()

    def test_grid_search_smoother_uses_smoothed_seed_values(self):
        """Trusted retrievals seed smoothing but are not restored afterward."""
        values = np.array(
            [
                [0.1, 0.1, 0.1],
                [0.1, 0.5, 0.1],
                [0.1, 0.1, 0.1],
            ],
            dtype=np.float32,
        )
        trusted = np.ones_like(values, dtype=bool)

        smoothed = MultiGridSolver._smooth_grid_search_field(
            values,
            gamma=2.0,
            delta=1.0,
            n_iter=40,
            trusted_mask=trusted,
        )

        assert smoothed[1, 1] < 0.25
        assert smoothed[1, 1] != pytest.approx(values[1, 1])


class TestCostFunctionConfig:
    """Tests for CostFunctionConfig."""

    def test_defaults(self):
        """Default config should have sensible values."""
        config = CostFunctionConfig()

        assert config.aot_gamma == 10.0
        assert config.tcwv_gamma == 5.0
        assert config.aot_min == 0.001
        assert config.aot_max == 2.5

    def test_custom_values(self):
        """Custom values should be respected."""
        config = CostFunctionConfig(
            aot_gamma=20.0,
            tcwv_gamma=10.0,
        )

        assert config.aot_gamma == 20.0
        assert config.tcwv_gamma == 10.0


class TestMultiGridConfig:
    """Tests for MultiGridConfig."""

    def test_defaults(self):
        """Default config should have sensible values."""
        config = MultiGridConfig()

        assert config.n_levels == 6
        assert config.min_grid_size == 8
        assert config.aot_bounds == (0.001, 2.5)
        assert config.fixed_atmospheric_parameter == "none"
        assert config.quadratic_block_min_valid_fraction == pytest.approx(0.5)

    def test_bounds(self):
        """Bounds should be tuple of two floats."""
        config = MultiGridConfig(
            aot_bounds=(0.01, 1.5),
            tcwv_bounds=(0.5, 5.0),
        )

        assert config.aot_bounds[0] == 0.01
        assert config.tcwv_bounds[1] == 5.0


class TestStagedMultiGridSolver:
    def test_staged_solver_chains_supported_aot_tcwv_passes(
        self,
        monkeypatch,
        mock_atmospheric_state,
        mock_geometry,
        mock_surface_prior,
        mock_rt_model,
    ):
        from siac.domain import SensorBand

        calls = []
        bands = [
            SensorBand("B02", 490.0, 65.0, 10.0, 0),
            SensorBand("B04", 665.0, 30.0, 10.0, 1),
        ]
        shape = mock_atmospheric_state.aot.shape
        toa = xr.DataArray(
            np.full((2, *shape), 0.2, dtype=np.float32),
            dims=["band", "y", "x"],
        )
        cloud_mask = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"])
        water_mask = xr.DataArray(
            np.array([[False, True], [False, False]], dtype=bool),
            dims=["y", "x"],
        )

        def _fake_solve(
            solver,
            toa,
            surface_prior,
            geometry,
            atmo_prior,
            rt_model,
            cloud_mask,
            bands,
            sharp_transition_mask=None,
            water_mask=None,
        ):
            del surface_prior, geometry, rt_model, cloud_mask, sharp_transition_mask
            call_index = len(calls)
            calls.append(
                {
                    "fixed": solver.config.fixed_atmospheric_parameter,
                    "aot_mean": float(atmo_prior.aot.mean()),
                    "bands": [band.name for band in bands],
                    "toa_shape": toa.shape,
                    "water_mask": water_mask,
                }
            )
            if call_index == 0:
                aot = xr.full_like(atmo_prior.aot, 0.31, dtype=np.float32)
                tcwv = atmo_prior.tcwv
            else:
                aot = atmo_prior.aot
                tcwv = xr.full_like(atmo_prior.tcwv, 3.2, dtype=np.float32)
            return SolverResult(
                aot=aot,
                tcwv=tcwv,
                aot_unc=atmo_prior.aot_unc,
                tcwv_unc=atmo_prior.tcwv_unc,
                n_iterations=call_index + 1,
                final_cost=float(call_index),
                success=True,
                message="ok",
                level_history=[{"level": 0}],
            )

        monkeypatch.setattr(MultiGridSolver, "solve", _fake_solve)
        config = MultiGridConfig(
            stages=(
                SolverStageConfig(
                    name="aot_pass",
                    solve=("aot",),
                    fixed=("tcwv", "tco3"),
                    bands=("B02",),
                ),
                SolverStageConfig(
                    name="tcwv_pass",
                    solve=("tcwv",),
                    fixed=("aot", "tco3"),
                    bands=("B04",),
                ),
            )
        )

        result = StagedMultiGridSolver(config).solve(
            toa=toa,
            surface_prior=mock_surface_prior,
            geometry=mock_geometry,
            atmo_prior=mock_atmospheric_state,
            rt_model=mock_rt_model,
            cloud_mask=cloud_mask,
            bands=bands,
            water_mask=water_mask,
        )

        assert calls[0]["fixed"] == "tcwv"
        assert calls[0]["bands"] == ["B02"]
        assert calls[0]["toa_shape"] == (1, *shape)
        assert calls[0]["water_mask"] is water_mask
        assert calls[1]["fixed"] == "aot"
        assert calls[1]["bands"] == ["B04"]
        assert calls[1]["aot_mean"] == pytest.approx(0.31)
        assert calls[1]["water_mask"] is water_mask
        assert result.n_iterations == 3
        assert result.level_history[0]["stage"] == "aot_pass"
        assert result.level_history[1]["stage"] == "tcwv_pass"
        assert result.atmo_state is not None
        np.testing.assert_allclose(result.atmo_state.aot, 0.31)
        np.testing.assert_allclose(result.atmo_state.tcwv, 3.2)
        np.testing.assert_allclose(result.atmo_state.tco3, mock_atmospheric_state.tco3)

    def test_staged_solver_rejects_tco3_solve_until_rt_support_exists(
        self,
        mock_atmospheric_state,
        mock_geometry,
        mock_surface_prior,
        mock_rt_model,
    ):
        from siac.domain import SensorBand

        shape = mock_atmospheric_state.aot.shape
        config = MultiGridConfig(
            stages=(
                SolverStageConfig(
                    name="ozone_pass",
                    solve=("tco3",),
                    fixed=("aot", "tcwv"),
                ),
            )
        )
        solver = StagedMultiGridSolver(config)

        with pytest.raises(NotImplementedError, match="TCO3"):
            solver.solve(
                toa=xr.DataArray(
                    np.full((1, *shape), 0.2, dtype=np.float32),
                    dims=["band", "y", "x"],
                ),
                surface_prior=mock_surface_prior,
                geometry=mock_geometry,
                atmo_prior=mock_atmospheric_state,
                rt_model=mock_rt_model,
                cloud_mask=xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"]),
                bands=[SensorBand("B09", 945.0, 20.0, 60.0, 0)],
            )


class TestSharpTransitionObservationExclusion:
    def test_build_solver_valid_mask_honors_sharp_transition_mask(self) -> None:
        cloud_mask = xr.DataArray(
            np.array([[False, False], [False, False]], dtype=bool),
            dims=["y", "x"],
        )
        toa = xr.DataArray(
            np.full((1, 2, 2), 0.2, dtype=np.float32),
            dims=["band", "y", "x"],
        )
        surface_prior = SurfacePrior(
            boa=xr.DataArray(np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"]),
            boa_unc=xr.DataArray(np.full((2, 2), 0.01, dtype=np.float32), dims=["y", "x"]),
            kernels=None,
            mask=xr.DataArray(np.ones((2, 2), dtype=bool), dims=["y", "x"]),
        )
        exclusion_mask = xr.DataArray(
            np.array([[False, True], [False, False]], dtype=bool),
            dims=["y", "x"],
        )

        valid = build_solver_valid_mask(
            cloud_mask,
            toa,
            surface_prior,
            sharp_transition_mask=exclusion_mask,
        )

        expected = np.array([[True, False], [True, True]], dtype=bool)
        np.testing.assert_array_equal(valid.values, expected)

    def test_build_solver_valid_mask_honors_water_mask(self) -> None:
        cloud_mask = xr.DataArray(
            np.array([[False, False], [False, False]], dtype=bool),
            dims=["y", "x"],
        )
        toa = xr.DataArray(
            np.full((1, 2, 2), 0.2, dtype=np.float32),
            dims=["band", "y", "x"],
        )
        surface_prior = SurfacePrior(
            boa=xr.DataArray(np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"]),
            boa_unc=xr.DataArray(np.full((2, 2), 0.01, dtype=np.float32), dims=["y", "x"]),
            kernels=None,
            mask=xr.DataArray(np.ones((2, 2), dtype=bool), dims=["y", "x"]),
        )
        water_mask = xr.DataArray(
            np.array([[False, True], [False, False]], dtype=bool),
            dims=["y", "x"],
        )

        valid = build_solver_valid_mask(
            cloud_mask,
            toa,
            surface_prior,
            water_mask=water_mask,
        )

        expected = np.array([[True, False], [True, True]], dtype=bool)
        np.testing.assert_array_equal(valid.values, expected)

    def test_lbfgsb_solver_excluded_water_pixels_are_gap_filled(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from siac.domain import SensorBand

        class _DummyRT:
            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):
                del geometry, atmo_state, band, compute_jacobian
                raise AssertionError("Optimisation is stubbed in this test")

            def supports_jacobian(self):
                return True

            @property
            def backend_name(self):
                return "dummy"

            def is_available_for_sensor(self, sensor_id, satellite_id):
                del sensor_id, satellite_id
                return True

        class _OptResult:
            nit = 1
            fun = 0.0
            success = True
            message = "ok"

        solver = MultiGridSolver(MultiGridConfig(n_levels=1, min_grid_size=2, max_iter_per_level=1))
        shape = (2, 2)
        toa = xr.DataArray(np.full((1, *shape), 0.2, dtype=np.float32), dims=["band", "y", "x"])
        cloud_mask = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"])
        water_mask = xr.DataArray(
            np.array([[False, True], [False, False]], dtype=bool), dims=["y", "x"]
        )
        surface_prior = SurfacePrior(
            boa=xr.DataArray(np.full(shape, 0.12, dtype=np.float32), dims=["y", "x"]),
            boa_unc=xr.DataArray(np.full(shape, 0.02, dtype=np.float32), dims=["y", "x"]),
            kernels=None,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )
        from siac.runtime import GeometryAngles

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

        def _fake_optimize_level(self, cost_func, aot_level, tcwv_level, level):  # noqa: ANN001
            del self, cost_func, aot_level, tcwv_level, level
            return (
                np.array([[0.4, 9.0], [0.4, 0.4]], dtype=np.float32),
                np.array([[3.0, 99.0], [3.0, 3.0]], dtype=np.float32),
                _OptResult(),
            )

        def _fake_uncertainties(self, aot, tcwv, atmo_prior_arg, cost_func):  # noqa: ANN001
            del self, aot, tcwv, atmo_prior_arg, cost_func
            return (
                np.full(shape, 0.02, dtype=np.float32),
                np.full(shape, 0.1, dtype=np.float32),
            )

        monkeypatch.setattr(MultiGridSolver, "_optimize_level", _fake_optimize_level)
        monkeypatch.setattr(MultiGridSolver, "_estimate_uncertainties", _fake_uncertainties)

        result = solver.solve(
            toa,
            surface_prior,
            geometry,
            atmo_prior,
            _DummyRT(),
            cloud_mask,
            [SensorBand("B02", 490.0, 65.0, 10.0, 0)],
            water_mask=water_mask,
        )

        assert result.qa is not None
        assert bool(result.qa["water_mask_excluded"].values[0, 1])
        assert result.aot.values[0, 1] == pytest.approx(0.4, rel=1e-6)
        assert result.tcwv.values[0, 1] == pytest.approx(3.0, rel=1e-6)

    def test_build_solver_qa_dataset_exposes_sharp_transition_layer(self) -> None:
        solver = MultiGridSolver()
        template = xr.DataArray(np.zeros((2, 2), dtype=np.float32), dims=["y", "x"])

        qa = solver._build_solver_qa_dataset(
            template=template,
            valid_mask=np.array([[True, True], [False, True]], dtype=bool),
            aot=np.full((2, 2), 0.2, dtype=np.float32),
            tcwv=np.full((2, 2), 2.0, dtype=np.float32),
            invalid_mask=np.array([[False, False], [False, True]], dtype=bool),
            zero_obs_mask=np.array([[False, False], [True, False]], dtype=bool),
            insufficient_support_mask=np.array([[False, False], [True, False]], dtype=bool),
            no_observation_mask=np.array([[True, False], [False, False]], dtype=bool),
            sharp_transition_mask=np.array([[False, True], [True, False]], dtype=bool),
            water_mask=np.array([[False, False], [False, True]], dtype=bool),
        )

        assert "sharp_transition_excluded" in qa.data_vars
        assert "water_mask_excluded" in qa.data_vars
        assert "insufficient_observation_support" in qa.data_vars
        assert "no_observation" in qa.data_vars
        np.testing.assert_array_equal(
            qa["sharp_transition_excluded"].values,
            np.array([[False, True], [True, False]], dtype=bool),
        )
        np.testing.assert_array_equal(
            qa["water_mask_excluded"].values,
            np.array([[False, False], [False, True]], dtype=bool),
        )
        np.testing.assert_array_equal(
            qa["zero_obs_support"].values,
            np.array([[False, False], [True, False]], dtype=bool),
        )
        np.testing.assert_array_equal(
            qa["insufficient_observation_support"].values,
            np.array([[False, False], [True, False]], dtype=bool),
        )
        np.testing.assert_array_equal(
            qa["no_observation"].values,
            np.array([[True, False], [False, False]], dtype=bool),
        )
        np.testing.assert_array_equal(
            qa["low_quality"].values,
            np.array([[True, True], [True, True]], dtype=bool),
        )

    def test_build_solver_qa_dataset_includes_fitting_cost_when_provided(self) -> None:
        solver = MultiGridSolver()
        template = xr.DataArray(np.zeros((2, 2), dtype=np.float32), dims=["y", "x"])
        cost_map = np.array([[0.01, 0.02], [0.03, 0.04]], dtype=np.float32)

        qa = solver._build_solver_qa_dataset(
            template=template,
            valid_mask=np.ones((2, 2), dtype=bool),
            aot=np.full((2, 2), 0.2, dtype=np.float32),
            tcwv=np.full((2, 2), 2.0, dtype=np.float32),
            invalid_mask=None,
            zero_obs_mask=None,
            insufficient_support_mask=None,
            no_observation_mask=None,
            sharp_transition_mask=None,
            water_mask=None,
            fitting_cost=cost_map,
        )

        assert "fitting_cost" in qa.data_vars
        np.testing.assert_array_equal(qa["fitting_cost"].values, cost_map)
        assert qa["fitting_cost"].dtype == np.float32

    def test_build_solver_qa_dataset_masks_unsupported_and_nonfinite_fitting_cost_pixels_to_nan(
        self,
    ) -> None:
        solver = MultiGridSolver()
        template = xr.DataArray(np.zeros((2, 2), dtype=np.float32), dims=["y", "x"])
        cost_map = np.array([[0.01, np.inf], [0.03, 0.04]], dtype=np.float32)

        qa = solver._build_solver_qa_dataset(
            template=template,
            valid_mask=np.array([[True, True], [True, False]], dtype=bool),
            aot=np.full((2, 2), 0.2, dtype=np.float32),
            tcwv=np.full((2, 2), 2.0, dtype=np.float32),
            invalid_mask=None,
            zero_obs_mask=None,
            insufficient_support_mask=np.array([[False, False], [True, False]], dtype=bool),
            no_observation_mask=None,
            sharp_transition_mask=None,
            water_mask=None,
            fitting_cost=cost_map,
        )

        fitting_cost = qa["fitting_cost"].values
        assert fitting_cost.dtype == np.float32
        assert fitting_cost[0, 0] == pytest.approx(np.float32(0.01))
        assert np.isnan(fitting_cost[0, 1])
        assert np.isnan(fitting_cost[1, 0])
        assert np.isnan(fitting_cost[1, 1])

    def test_build_solver_qa_dataset_omits_fitting_cost_when_none(self) -> None:
        solver = MultiGridSolver()
        template = xr.DataArray(np.zeros((2, 2), dtype=np.float32), dims=["y", "x"])

        qa = solver._build_solver_qa_dataset(
            template=template,
            valid_mask=np.ones((2, 2), dtype=bool),
            aot=np.full((2, 2), 0.2, dtype=np.float32),
            tcwv=np.full((2, 2), 2.0, dtype=np.float32),
            invalid_mask=None,
            zero_obs_mask=None,
            insufficient_support_mask=None,
            no_observation_mask=None,
            sharp_transition_mask=None,
            water_mask=None,
            fitting_cost=None,
        )

        assert "fitting_cost" not in qa.data_vars


class TestMultiGridSolver:
    """Tests for MultiGridSolver."""

    def test_creation(self):
        """Solver should be creatable."""
        solver = MultiGridSolver()
        assert solver.config.n_levels == 6

    def test_with_config(self):
        """Solver should accept custom config."""
        config = MultiGridConfig(n_levels=4, aot_gamma=15.0)
        solver = MultiGridSolver(config)

        assert solver.config.n_levels == 4
        assert solver.config.aot_gamma == 15.0

    def test_compute_grid_levels(self):
        """Grid levels should be computed correctly."""
        solver = MultiGridSolver(MultiGridConfig(n_levels=5, min_grid_size=8))

        shapes = solver._compute_grid_levels((256, 256))

        # Should have progressively finer grids
        assert len(shapes) <= 5
        assert shapes[0][0] <= shapes[-1][0]  # Coarse to fine

    def test_resample_field(self):
        """Field resampling should work."""
        solver = MultiGridSolver()

        field = np.random.default_rng(2).standard_normal((100, 100))
        resampled = solver._resample_field(field, (50, 50))

        assert resampled.shape == (50, 50)

    def test_resample_same_size(self):
        """Same-size resampling should return identical array."""
        solver = MultiGridSolver()

        field = np.random.default_rng(3).standard_normal((64, 64))
        resampled = solver._resample_field(field, (64, 64))

        np.testing.assert_array_equal(resampled, field)

    def test_resample_mask_to_grid_uses_center_assigned_footprint(self):
        solver = MultiGridSolver()
        mask = xr.DataArray(np.ones((5, 5), dtype=bool), dims=["y", "x"])
        mask.values[-1, -1] = False

        resampled = solver._resample_mask_to_grid(mask, (2, 2))

        expected = np.array([[True, True], [True, False]])
        np.testing.assert_array_equal(resampled.values, expected)

    def test_solve_returns_prior_when_no_valid_pixels(self):
        """Solver should short-circuit when cloud/quality masking leaves no valid pixels."""
        solver = MultiGridSolver()
        shape = (8, 8)

        toa = xr.DataArray(np.full((2, *shape), 0.2, dtype=np.float32), dims=["band", "y", "x"])
        cloud_mask = xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"])

        from siac.domain import SensorBand
        from siac.runtime import BRDFKernelWeights, GeometryAngles, SurfacePrior

        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 1.0), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.3), dims=["y", "x"]),
        )
        atmo_prior = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        )
        brdf = BRDFKernelWeights(
            f0=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            f1=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            f2=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
            f0_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            f1_unc=xr.DataArray(np.full(shape, 0.005), dims=["y", "x"]),
            f2_unc=xr.DataArray(np.full(shape, 0.002), dims=["y", "x"]),
        )
        surface_prior = SurfacePrior(
            boa=xr.DataArray(np.full(shape, 0.12), dims=["y", "x"]),
            boa_unc=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
            kernels=brdf,
            mask=xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"]),
        )

        class _DummyRT:
            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
                raise NotImplementedError

            def supports_jacobian(self) -> bool:
                return False

            @property
            def backend_name(self) -> str:
                return "dummy"

            def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool:
                _ = (sensor_id, satellite_id)
                return True

        result = solver.solve(
            toa=toa,
            surface_prior=surface_prior,
            geometry=geometry,
            atmo_prior=atmo_prior,
            rt_model=_DummyRT(),
            cloud_mask=cloud_mask,
            bands=[SensorBand("B01", 443.0, 20.0, 60.0, 0)],
        )

        assert result.n_iterations == 0
        assert not result.success
        assert "No valid pixels" in result.message
        np.testing.assert_allclose(result.aot.values, atmo_prior.aot.values)

        nodata_surface_prior = SurfacePrior(
            boa=surface_prior.boa,
            boa_unc=surface_prior.boa_unc,
            kernels=brdf,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )
        nodata_result = solver.solve(
            toa=xr.DataArray(
                np.full((2, *shape), np.nan, dtype=np.float32), dims=["band", "y", "x"]
            ),
            surface_prior=nodata_surface_prior,
            geometry=geometry,
            atmo_prior=atmo_prior,
            rt_model=_DummyRT(),
            cloud_mask=xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"]),
            bands=[SensorBand("B01", 443.0, 20.0, 60.0, 0)],
        )

        assert not np.any(np.isfinite(nodata_result.aot.values))
        assert not np.any(np.isfinite(nodata_result.tcwv.values))

    def test_solve_uses_grid_search_when_rt_has_no_jacobian(self):
        """Solver should switch to grid-search + quadratic fit for non-jacobian RT backends."""
        from siac.domain import SensorBand
        from siac.runtime import BRDFKernelWeights, GeometryAngles, RTCoefficients, SurfacePrior

        shape = (10, 10)
        config = MultiGridConfig(
            n_levels=2,
            min_grid_size=8,
            grid_search_aot_points=5,
            grid_search_tcwv_points=5,
        )
        solver = MultiGridSolver(config)

        toa = xr.DataArray(
            np.stack(
                [
                    np.full(shape, 0.25, dtype=np.float32),
                    np.full(shape, 0.28, dtype=np.float32),
                ]
            ),
            dims=["band", "y", "x"],
        )
        cloud_mask = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"])

        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
        )
        atmo_prior = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.2), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.0), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        )
        brdf = BRDFKernelWeights(
            f0=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            f1=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            f2=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
            f0_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            f1_unc=xr.DataArray(np.full(shape, 0.005), dims=["y", "x"]),
            f2_unc=xr.DataArray(np.full(shape, 0.002), dims=["y", "x"]),
        )
        surface_prior = SurfacePrior(
            boa=xr.DataArray(
                np.stack([np.full(shape, 0.2), np.full(shape, 0.22)]), dims=["band", "y", "x"]
            ),
            boa_unc=xr.DataArray(
                np.stack([np.full(shape, 0.02), np.full(shape, 0.02)]), dims=["band", "y", "x"]
            ),
            kernels=brdf,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )

        class _DummyRT:
            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
                _ = (band, compute_jacobian)
                factor = 1.0 / np.maximum(
                    0.2 + atmo_state.aot.values + 0.1 * atmo_state.tcwv.values,
                    1e-6,
                )
                xap = xr.DataArray(factor, dims=geometry.sza.dims, coords=geometry.sza.coords)
                xbp = xr.zeros_like(xap)
                xcp = xr.zeros_like(xap)
                return RTCoefficients(xap=xap, xbp=xbp, xcp=xcp)

            def supports_jacobian(self) -> bool:
                return False

            @property
            def backend_name(self) -> str:
                return "dummy"

            def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool:
                _ = (sensor_id, satellite_id)
                return True

        bands = [
            SensorBand("B02", 490.0, 65.0, 10.0, 0),
            SensorBand("B03", 560.0, 35.0, 10.0, 1),
        ]

        result = solver.solve(
            toa=toa,
            surface_prior=surface_prior,
            geometry=geometry,
            atmo_prior=atmo_prior,
            rt_model=_DummyRT(),
            cloud_mask=cloud_mask,
            bands=bands,
        )

        assert result.success
        assert len(result.level_history) == 1
        assert tuple(result.level_history[0]["shape"]) == shape
        assert all(item["method"] == "grid_search" for item in result.level_history)
        assert np.isfinite(result.aot.values).all()
        assert np.isfinite(result.tcwv.values).all()
        assert (result.aot_unc.values > 0).all()
        assert (result.tcwv_unc.values > 0).all()

    def test_solve_skips_coarse_levels_for_grid_search_backend(self):
        """Grid-search backend should solve once at full resolution, not per multigrid level."""
        from siac.domain import SensorBand
        from siac.runtime import BRDFKernelWeights, GeometryAngles, RTCoefficients, SurfacePrior

        shape = (64, 64)
        config = MultiGridConfig(
            n_levels=4,
            min_grid_size=8,
            grid_search_aot_points=5,
            grid_search_tcwv_points=5,
        )
        solver = MultiGridSolver(config)

        toa = xr.DataArray(
            np.stack(
                [
                    np.full(shape, 0.25, dtype=np.float32),
                    np.full(shape, 0.28, dtype=np.float32),
                ]
            ),
            dims=["band", "y", "x"],
        )
        cloud_mask = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"])

        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
        )
        atmo_prior = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.2), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.0), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        )
        brdf = BRDFKernelWeights(
            f0=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            f1=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            f2=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
            f0_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            f1_unc=xr.DataArray(np.full(shape, 0.005), dims=["y", "x"]),
            f2_unc=xr.DataArray(np.full(shape, 0.002), dims=["y", "x"]),
        )
        surface_prior = SurfacePrior(
            boa=xr.DataArray(
                np.stack([np.full(shape, 0.2), np.full(shape, 0.22)]),
                dims=["band", "y", "x"],
            ),
            boa_unc=xr.DataArray(
                np.stack([np.full(shape, 0.02), np.full(shape, 0.02)]),
                dims=["band", "y", "x"],
            ),
            kernels=brdf,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )

        class _DummyRT:
            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
                _ = (band, compute_jacobian)
                factor = 1.0 / np.maximum(
                    0.2 + atmo_state.aot.values + 0.1 * atmo_state.tcwv.values,
                    1e-6,
                )
                xap = xr.DataArray(factor, dims=geometry.sza.dims, coords=geometry.sza.coords)
                xbp = xr.zeros_like(xap)
                xcp = xr.zeros_like(xap)
                return RTCoefficients(xap=xap, xbp=xbp, xcp=xcp)

            def supports_jacobian(self) -> bool:
                return False

            @property
            def backend_name(self) -> str:
                return "dummy"

            def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool:
                _ = (sensor_id, satellite_id)
                return True

        bands = [
            SensorBand("B02", 490.0, 65.0, 10.0, 0),
            SensorBand("B03", 560.0, 35.0, 10.0, 1),
        ]

        result = solver.solve(
            toa=toa,
            surface_prior=surface_prior,
            geometry=geometry,
            atmo_prior=atmo_prior,
            rt_model=_DummyRT(),
            cloud_mask=cloud_mask,
            bands=bands,
        )

        assert result.success
        assert solver._compute_grid_levels(shape) == [(8, 8), (16, 16), (32, 32), (64, 64)]
        assert len(result.level_history) == 1
        assert result.n_iterations == config.grid_search_aot_points * config.grid_search_tcwv_points
        assert tuple(result.level_history[0]["shape"]) == shape
        assert result.level_history[0]["method"] == "grid_search"

    def test_solver_marks_supported_lower_bound_solution_as_low_quality(self):
        """Supported boundary solutions should be flagged as low-quality, not invalid."""
        from siac.domain import SensorBand
        from siac.runtime import BRDFKernelWeights, GeometryAngles, RTCoefficients, SurfacePrior

        shape = (2, 3)
        # Lower bound above the prior (0.3) so the supported solution is driven
        # to the AOT floor: the log-spaced grid otherwise resolves the prior.
        solver = MultiGridSolver(
            MultiGridConfig(
                grid_search_aot_points=5, grid_search_tcwv_points=3, aot_bounds=(0.4, 2.5)
            )
        )

        toa = xr.DataArray(
            np.stack([np.full(shape, 0.1, dtype=np.float32)]),
            dims=["band", "y", "x"],
        )
        cloud_mask = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"])
        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5, dtype=np.float32), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5, dtype=np.float32), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5, dtype=np.float32), dims=["y", "x"]),
        )
        atmo_prior = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.0, dtype=np.float32), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
        )
        brdf = BRDFKernelWeights(
            f0=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
            f1=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=["y", "x"]),
            f2=xr.DataArray(np.full(shape, 0.02, dtype=np.float32), dims=["y", "x"]),
            f0_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=["y", "x"]),
            f1_unc=xr.DataArray(np.full(shape, 0.005, dtype=np.float32), dims=["y", "x"]),
            f2_unc=xr.DataArray(np.full(shape, 0.002, dtype=np.float32), dims=["y", "x"]),
        )
        surface_prior = SurfacePrior(
            boa=xr.DataArray(
                np.stack([np.full(shape, 0.2, dtype=np.float32)]), dims=["band", "y", "x"]
            ),
            boa_unc=xr.DataArray(
                np.stack([np.full(shape, 0.02, dtype=np.float32)]), dims=["band", "y", "x"]
            ),
            kernels=brdf,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )

        class _DummyRT:
            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
                _ = (band, compute_jacobian)
                xap = xr.DataArray(
                    np.ones(shape, dtype=np.float32),
                    dims=geometry.sza.dims,
                    coords=geometry.sza.coords,
                )
                xbp = xr.DataArray(
                    0.05 * np.asarray(atmo_state.aot.values, dtype=np.float32),
                    dims=geometry.sza.dims,
                    coords=geometry.sza.coords,
                )
                return RTCoefficients(xap=xap, xbp=xbp, xcp=xr.zeros_like(xap))

            def supports_jacobian(self) -> bool:
                return False

            @property
            def backend_name(self) -> str:
                return "dummy"

            def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool:
                _ = (sensor_id, satellite_id)
                return True

        result = solver.solve(
            toa=toa,
            surface_prior=surface_prior,
            geometry=geometry,
            atmo_prior=atmo_prior,
            rt_model=_DummyRT(),
            cloud_mask=cloud_mask,
            bands=[SensorBand("B02", 490.0, 65.0, 10.0, 0)],
        )

        assert result.qa is not None
        np.testing.assert_allclose(
            result.aot.values, np.full(shape, solver.config.aot_bounds[0], dtype=np.float32)
        )
        assert not np.any(result.qa["invalid_retrieval"].values)
        assert not np.any(result.qa["zero_obs_support"].values)
        assert np.all(result.qa["aot_lower_boundary"].values)
        assert not np.any(result.qa["aot_upper_boundary"].values)
        assert np.all(result.qa["parameter_boundary"].values)
        assert np.all(result.qa["low_quality"].values)
        assert result.level_history[-1]["qa_final_aot_lower_boundary_pixels"] == pytest.approx(
            float(np.prod(shape))
        )
        assert result.level_history[-1]["qa_final_low_quality_pixels"] == pytest.approx(
            float(np.prod(shape))
        )

    def test_grid_search_prefers_rust_refiner_when_available(self, monkeypatch):
        """Grid-search refinement should call Rust helper."""
        from siac.algorithms.solver import multigrid as mg_mod
        from siac.domain import SensorBand
        from siac.runtime import BRDFKernelWeights, GeometryAngles, RTCoefficients, SurfacePrior

        shape = (4, 4)
        solver = MultiGridSolver(
            MultiGridConfig(grid_search_aot_points=3, grid_search_tcwv_points=3)
        )

        toa = xr.DataArray(
            np.stack([np.full(shape, 0.25, dtype=np.float32)]),
            dims=["band", "y", "x"],
        )
        mask = xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"])
        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
        )
        atmo_prior = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.2), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.0), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        )
        brdf = BRDFKernelWeights(
            f0=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            f1=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            f2=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
            f0_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            f1_unc=xr.DataArray(np.full(shape, 0.005), dims=["y", "x"]),
            f2_unc=xr.DataArray(np.full(shape, 0.002), dims=["y", "x"]),
        )
        surface_prior = SurfacePrior(
            boa=xr.DataArray(np.stack([np.full(shape, 0.2)]), dims=["band", "y", "x"]),
            boa_unc=xr.DataArray(np.stack([np.full(shape, 0.02)]), dims=["band", "y", "x"]),
            kernels=brdf,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )

        class _DummyRT:
            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
                _ = (band, compute_jacobian)
                factor = np.full(shape, 0.9, dtype=np.float32)
                xap = xr.DataArray(factor, dims=geometry.sza.dims, coords=geometry.sza.coords)
                return RTCoefficients(xap=xap, xbp=xr.zeros_like(xap), xcp=xr.zeros_like(xap))

            def supports_jacobian(self) -> bool:
                return False

            @property
            def backend_name(self) -> str:
                return "dummy"

            def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool:
                _ = (sensor_id, satellite_id)
                return True

        called = {"n": 0}

        def _fake_rust_refiner(
            costs,
            obs_counts,
            aot_axis,
            tcwv_axis,
            valid_mask,
            fixed_parameter="none",
        ):  # noqa: ANN001
            assert fixed_parameter == "none"
            called["n"] += 1
            out_aot = np.full(valid_mask.shape, 0.33, dtype=np.float32)
            out_tcwv = np.full(valid_mask.shape, 1.23, dtype=np.float32)
            out_aot_unc = np.full(valid_mask.shape, 0.07, dtype=np.float32)
            out_tcwv_unc = np.full(valid_mask.shape, 0.21, dtype=np.float32)
            assert costs.shape[:2] == (3, 3)
            assert obs_counts.shape == costs.shape
            assert aot_axis.shape[0] == 3
            assert tcwv_axis.shape[0] == 3
            return (
                out_aot,
                out_tcwv,
                out_aot_unc,
                out_tcwv_unc,
                np.zeros(valid_mask.shape, dtype=bool),
                np.zeros(valid_mask.shape, dtype=bool),
                np.zeros(valid_mask.shape, dtype=bool),
                np.zeros(valid_mask.shape, dtype=bool),
            )

        monkeypatch.setattr(mg_mod, "quadratic_refine_grid_search_qa", _fake_rust_refiner)

        out_aot, out_tcwv, out_aot_unc, out_tcwv_unc, _ = solver._solve_level_grid_search(
            toa=toa,
            surface_prior=surface_prior,
            geometry=geometry,
            atmo_prior=atmo_prior,
            rt_model=_DummyRT(),
            mask=mask,
            bands=[SensorBand("B02", 490.0, 65.0, 10.0, 0)],
            cost_config=CostFunctionConfig(),
        )

        assert called["n"] == 1
        np.testing.assert_allclose(out_aot, 0.33)
        np.testing.assert_allclose(out_tcwv, 1.23)
        np.testing.assert_allclose(out_aot_unc, 0.07)
        np.testing.assert_allclose(out_tcwv_unc, 0.21)

    def test_grid_search_can_fix_tcwv_to_prior(self, monkeypatch):
        """Fixed TCWV mode should evaluate only the AOT axis and return TCWV prior."""
        from siac.algorithms.solver import multigrid as mg_mod
        from siac.domain import SensorBand
        from siac.runtime import BRDFKernelWeights, GeometryAngles, RTCoefficients, SurfacePrior

        shape = (2, 3)
        solver = MultiGridSolver(
            MultiGridConfig(
                grid_search_aot_points=5,
                grid_search_tcwv_points=5,
                fixed_atmospheric_parameter="tcwv",
            )
        )
        toa = xr.DataArray(
            np.stack([np.full(shape, 0.25, dtype=np.float32)]), dims=["band", "y", "x"]
        )
        mask = xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"])
        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
        )
        tcwv_prior_values = np.array([[np.nan, 1.5, 2.0], [2.5, 3.0, 3.5]], dtype=np.float32)
        provider_tcwv_values = tcwv_prior_values.copy()
        provider_tcwv_values[0, 0] = np.nanmean(tcwv_prior_values)
        tcwv_unc_values = np.full(shape, 0.3, dtype=np.float32)
        atmo_prior = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.2, dtype=np.float32), dims=["y", "x"]),
            tcwv=xr.DataArray(tcwv_prior_values, dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(tcwv_unc_values, dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
        )
        brdf = BRDFKernelWeights(
            f0=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            f1=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            f2=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
            f0_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            f1_unc=xr.DataArray(np.full(shape, 0.005), dims=["y", "x"]),
            f2_unc=xr.DataArray(np.full(shape, 0.002), dims=["y", "x"]),
        )
        surface_prior = SurfacePrior(
            boa=xr.DataArray(
                np.stack([np.full(shape, 0.2, dtype=np.float32)]), dims=["band", "y", "x"]
            ),
            boa_unc=xr.DataArray(
                np.stack([np.full(shape, 0.02, dtype=np.float32)]), dims=["band", "y", "x"]
            ),
            kernels=brdf,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )

        seen_tcwv: list[np.ndarray] = []

        class _DummyRT:
            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
                _ = (band, compute_jacobian)
                seen_tcwv.append(np.asarray(atmo_state.tcwv.values, dtype=np.float32).copy())
                xap = xr.DataArray(
                    np.full(geometry.sza.shape, 0.9, dtype=np.float32), dims=["y", "x"]
                )
                return RTCoefficients(xap=xap, xbp=xr.zeros_like(xap), xcp=xr.zeros_like(xap))

            def supports_jacobian(self) -> bool:
                return False

            @property
            def backend_name(self) -> str:
                return "dummy"

            def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool:
                _ = (sensor_id, satellite_id)
                return True

        def _fake_pixel_eval(
            _provider,
            aot_axis,
            tcwv_axis,
            _toa,
            _boa_prior,
            _boa_unc,
            _band_weights,
            _valid_mask,
            _aot_prior,
            _tcwv_prior,
            _aot_prior_unc,
            _tcwv_prior_unc,
            fixed_parameter,
        ):  # noqa: ANN001
            assert fixed_parameter == "tcwv"
            assert aot_axis.shape == (5,)
            assert tcwv_axis.shape == (1,)
            _provider(float(aot_axis[2]), 999.0)
            costs = np.full((aot_axis.shape[0], 1, *shape), 10.0, dtype=np.float32)
            costs[2, 0, :, :] = 0.0
            return costs, np.ones_like(costs, dtype=np.uint16)

        def _fake_refiner(costs, obs_counts, aot_axis, tcwv_axis, valid_mask, fixed_parameter):  # noqa: ANN001
            _ = (costs, obs_counts, aot_axis, tcwv_axis)
            assert fixed_parameter == "tcwv"
            return (
                np.full(valid_mask.shape, 0.44, dtype=np.float32),
                np.full(valid_mask.shape, -999.0, dtype=np.float32),
                np.full(valid_mask.shape, 0.07, dtype=np.float32),
                np.full(valid_mask.shape, 0.21, dtype=np.float32),
                np.zeros(valid_mask.shape, dtype=bool),
                np.zeros(valid_mask.shape, dtype=bool),
                np.zeros(valid_mask.shape, dtype=bool),
                np.zeros(valid_mask.shape, dtype=bool),
            )

        monkeypatch.setattr(
            mg_mod, "evaluate_grid_search_cost_cube_with_provider_qa", _fake_pixel_eval
        )
        monkeypatch.setattr(mg_mod, "quadratic_refine_grid_search_qa", _fake_refiner)

        out_aot, out_tcwv, _out_aot_unc, out_tcwv_unc, diag = solver._solve_level_grid_search(
            toa=toa,
            surface_prior=surface_prior,
            geometry=geometry,
            atmo_prior=atmo_prior,
            rt_model=_DummyRT(),
            mask=mask,
            bands=[SensorBand("B02", 490.0, 65.0, 10.0, 0)],
            cost_config=CostFunctionConfig(aot_gamma=0.0, tcwv_gamma=0.0),
        )

        np.testing.assert_allclose(seen_tcwv[-1], provider_tcwv_values)
        np.testing.assert_allclose(out_aot, 0.44)
        np.testing.assert_allclose(out_tcwv, tcwv_prior_values, equal_nan=True)
        np.testing.assert_allclose(out_tcwv_unc, tcwv_unc_values)
        assert diag["evaluations"] == pytest.approx(5.0)

    @pytest.mark.skipif(
        not HAS_NATIVE_BLOCK_COST_CUBE,
        reason="native extension does not export block cost cube helper in this environment",
    )
    def test_rust_block_cost_cube_aggregates_pixel_costs(self):
        """Rust block cost cube should sum per-pixel candidate costs into NxN blocks."""
        aot_axis = np.array([0.2, 0.5, 0.8], dtype=np.float32)
        tcwv_axis = np.array([1.0, 2.0, 3.0], dtype=np.float32)
        shape = (4, 5)
        block_size = 3
        valid_mask = np.array(
            [
                [True, True, False, True, True],
                [True, True, True, True, False],
                [False, True, True, True, True],
                [True, False, True, True, True],
            ],
            dtype=bool,
        )

        toa = np.full((1, *shape), 0.25, dtype=np.float32)
        boa_prior = np.full((1, *shape), 0.2, dtype=np.float32)
        boa_unc = np.full((1, *shape), 0.02, dtype=np.float32)
        band_weights = np.array([1.0], dtype=np.float32)
        aot_prior = np.full(shape, 0.2, dtype=np.float32)
        tcwv_prior = np.full(shape, 2.0, dtype=np.float32)
        aot_prior_unc = np.full(shape, 0.05, dtype=np.float32)
        tcwv_prior_unc = np.full(shape, 0.3, dtype=np.float32)

        def _provider(
            _aot_val: float, _tcwv_val: float
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            xap = np.ones((1, *shape), dtype=np.float32)
            zeros = np.zeros((1, *shape), dtype=np.float32)
            return xap, zeros, zeros

        full_costs, full_counts = evaluate_grid_search_cost_cube_with_provider_qa(
            _provider,
            aot_axis,
            tcwv_axis,
            toa,
            boa_prior,
            boa_unc,
            band_weights,
            valid_mask,
            aot_prior,
            tcwv_prior,
            aot_prior_unc,
            tcwv_prior_unc,
        )
        block_costs, block_counts, block_valid = (
            evaluate_block_grid_search_cost_cube_with_provider_qa(
                _provider,
                aot_axis,
                tcwv_axis,
                toa,
                boa_prior,
                boa_unc,
                band_weights,
                valid_mask,
                aot_prior,
                tcwv_prior,
                aot_prior_unc,
                tcwv_prior_unc,
                block_size,
            )
        )

        full_costs = np.asarray(full_costs, dtype=np.float32)
        full_counts = np.asarray(full_counts, dtype=np.uint16)
        block_costs = np.asarray(block_costs, dtype=np.float32)
        block_counts = np.asarray(block_counts, dtype=np.uint16)
        block_valid = np.asarray(block_valid, dtype=bool)

        assert block_costs.shape == (aot_axis.size, tcwv_axis.size, 2, 2)
        assert block_counts.shape == block_costs.shape

        for by in range(block_costs.shape[2]):
            y0 = by * block_size
            y1 = min(y0 + block_size, shape[0])
            for bx in range(block_costs.shape[3]):
                x0 = bx * block_size
                x1 = min(x0 + block_size, shape[1])
                expected_cost = np.sum(full_costs[:, :, y0:y1, x0:x1], axis=(2, 3))
                expected_count = np.sum(full_counts[:, :, y0:y1, x0:x1], axis=(2, 3))
                np.testing.assert_allclose(
                    block_costs[:, :, by, bx],
                    expected_cost,
                    rtol=1.0e-5,
                    atol=1.0e-4,
                )
                np.testing.assert_array_equal(block_counts[:, :, by, bx], expected_count)
        np.testing.assert_array_equal(block_valid, _any_blocks_reference(valid_mask, block_size))

        def _coarse_provider(
            _aot_val: float, _tcwv_val: float
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            xap = np.ones((1, *block_valid.shape), dtype=np.float32)
            zeros = np.zeros((1, *block_valid.shape), dtype=np.float32)
            return xap, zeros, zeros

        coarse_costs, coarse_counts, coarse_valid = (
            evaluate_block_grid_search_cost_cube_with_provider_qa(
                _coarse_provider,
                aot_axis,
                tcwv_axis,
                toa,
                boa_prior,
                boa_unc,
                band_weights,
                valid_mask,
                aot_prior,
                tcwv_prior,
                aot_prior_unc,
                tcwv_prior_unc,
                block_size,
            )
        )
        np.testing.assert_allclose(np.asarray(coarse_costs, dtype=np.float32), block_costs)
        np.testing.assert_array_equal(np.asarray(coarse_counts, dtype=np.uint16), block_counts)
        np.testing.assert_array_equal(np.asarray(coarse_valid, dtype=bool), block_valid)

    def test_grid_search_block_solver_broadcasts_shared_block_solution(self, monkeypatch):
        """Block solver should return one shared solution per NxN block and broadcast it back."""
        from siac.algorithms.solver import multigrid as mg_mod
        from siac.domain import SensorBand
        from siac.runtime import BRDFKernelWeights, GeometryAngles, RTCoefficients, SurfacePrior

        shape = (5, 4)
        solver = MultiGridSolver(
            MultiGridConfig(
                grid_search_aot_points=5,
                grid_search_tcwv_points=5,
                aot_bounds=(0.2, 0.8),
                tcwv_bounds=(1.0, 3.0),
                quadratic_block_size=3,
            )
        )

        toa_values = np.full(shape, 0.25, dtype=np.float32)
        toa_values[0, 0] = np.nan
        toa_values[0, 1] = np.nan
        toa = xr.DataArray(np.stack([toa_values]), dims=["band", "y", "x"])
        mask = xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"])
        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
        )
        atmo_prior = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.2), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.0), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        )
        brdf = BRDFKernelWeights(
            f0=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            f1=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            f2=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
            f0_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            f1_unc=xr.DataArray(np.full(shape, 0.005), dims=["y", "x"]),
            f2_unc=xr.DataArray(np.full(shape, 0.002), dims=["y", "x"]),
        )
        surface_prior = SurfacePrior(
            boa=xr.DataArray(np.stack([np.full(shape, 0.2)]), dims=["band", "y", "x"]),
            boa_unc=xr.DataArray(np.stack([np.full(shape, 0.02)]), dims=["band", "y", "x"]),
            kernels=brdf,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )

        coeff_shapes: list[tuple[int, int]] = []

        class _DummyRT:
            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
                _ = (geometry, atmo_state, band, compute_jacobian)
                coeff_shapes.append(tuple(geometry.sza.shape))
                factor = np.full(geometry.sza.shape, 0.9, dtype=np.float32)
                xap = xr.DataArray(factor, dims=["y", "x"])
                return RTCoefficients(xap=xap, xbp=xr.zeros_like(xap), xcp=xr.zeros_like(xap))

            def supports_jacobian(self) -> bool:
                return False

            @property
            def backend_name(self) -> str:
                return "dummy"

            def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool:
                _ = (sensor_id, satellite_id)
                return True

        def _fail_pixel_eval(*_args, **_kwargs):  # noqa: ANN002, ANN003
            raise AssertionError(
                "Per-pixel evaluator should not be used when quadratic_block_size > 1"
            )

        def _fake_block_eval(
            _provider,
            aot_axis,
            tcwv_axis,
            _toa,
            _boa_prior,
            _boa_unc,
            _band_weights,
            _valid_mask,
            _aot_prior,
            _tcwv_prior,
            _aot_prior_unc,
            _tcwv_prior_unc,
            block_size,
            fixed_parameter="none",
        ):  # noqa: ANN001
            assert block_size == 3
            assert fixed_parameter == "none"
            xap, xbp, xcp = _provider(float(aot_axis[0]), float(tcwv_axis[0]))
            assert xap.shape == (1, 2, 2)
            assert xbp.shape == (1, 2, 2)
            assert xcp.shape == (1, 2, 2)
            costs = np.full((aot_axis.shape[0], tcwv_axis.shape[0], 2, 2), 50.0, dtype=np.float32)
            obs_counts = np.ones_like(costs, dtype=np.uint16)
            minima = {
                (0, 0): (1, 1),
                (0, 1): (2, 1),
                (1, 0): (1, 2),
                (1, 1): (2, 2),
            }
            for (by, bx), (ia, it) in minima.items():
                for da in (-1, 0, 1):
                    for dt in (-1, 0, 1):
                        ia_idx = ia + da
                        it_idx = it + dt
                        if 0 <= ia_idx < aot_axis.shape[0] and 0 <= it_idx < tcwv_axis.shape[0]:
                            costs[ia_idx, it_idx, by, bx] = float(da * da + dt * dt)
            block_valid = np.ones((2, 2), dtype=bool)
            return costs, obs_counts, block_valid

        monkeypatch.setattr(
            mg_mod, "evaluate_grid_search_cost_cube_with_provider_qa", _fail_pixel_eval
        )
        monkeypatch.setattr(
            mg_mod, "evaluate_block_grid_search_cost_cube_with_provider_qa", _fake_block_eval
        )

        out_aot, out_tcwv, out_aot_unc, out_tcwv_unc, level_diag = solver._solve_level_grid_search(
            toa=toa,
            surface_prior=surface_prior,
            geometry=geometry,
            atmo_prior=atmo_prior,
            rt_model=_DummyRT(),
            mask=mask,
            bands=[SensorBand("B02", 490.0, 65.0, 10.0, 0)],
            cost_config=CostFunctionConfig(aot_gamma=0.0, tcwv_gamma=0.0),
        )

        # Solver now log-spaces the AOT grid (see multigrid._log_aot_axis).
        aot_axis = np.geomspace(
            max(float(solver.config.aot_bounds[0]), 1e-4),
            float(solver.config.aot_bounds[1]),
            solver.config.grid_search_aot_points,
            dtype=np.float32,
        )
        tcwv_axis = np.linspace(
            solver.config.tcwv_bounds[0],
            solver.config.tcwv_bounds[1],
            solver.config.grid_search_tcwv_points,
            dtype=np.float32,
        )

        expected_aot = np.array(
            [
                [np.nan, np.nan, aot_axis[1], aot_axis[2]],
                [aot_axis[1], aot_axis[1], aot_axis[1], aot_axis[2]],
                [aot_axis[1], aot_axis[1], aot_axis[1], aot_axis[2]],
                [aot_axis[1], aot_axis[1], aot_axis[1], aot_axis[2]],
                [aot_axis[1], aot_axis[1], aot_axis[1], aot_axis[2]],
            ],
            dtype=np.float32,
        )
        expected_tcwv = np.array(
            [
                [np.nan, np.nan, tcwv_axis[1], tcwv_axis[1]],
                [tcwv_axis[1], tcwv_axis[1], tcwv_axis[1], tcwv_axis[1]],
                [tcwv_axis[1], tcwv_axis[1], tcwv_axis[1], tcwv_axis[1]],
                [tcwv_axis[2], tcwv_axis[2], tcwv_axis[2], tcwv_axis[2]],
                [tcwv_axis[2], tcwv_axis[2], tcwv_axis[2], tcwv_axis[2]],
            ],
            dtype=np.float32,
        )

        np.testing.assert_allclose(out_aot, expected_aot, equal_nan=True)
        np.testing.assert_allclose(out_tcwv, expected_tcwv, equal_nan=True)
        assert np.isnan(out_aot_unc[0, 0])
        assert np.isnan(out_tcwv_unc[0, 0])
        np.testing.assert_allclose(out_aot_unc[1:3, :3], out_aot_unc[1, 0])
        np.testing.assert_allclose(out_tcwv_unc[1:3, :3], out_tcwv_unc[1, 0])
        assert level_diag["valid_pixels"] == pytest.approx(float(np.prod(shape) - 2))
        assert level_diag["qa_zero_obs_pixels"] == pytest.approx(2.0)
        assert level_diag["qa_no_observation_pixels"] == pytest.approx(2.0)
        assert coeff_shapes == [(2, 2)]

    def test_grid_search_block_solver_requires_half_valid_block_support(self, monkeypatch):
        """Blocked retrieval should skip blocks with too few valid observation/prior pixels."""
        from siac.algorithms.solver import multigrid as mg_mod
        from siac.domain import SensorBand
        from siac.runtime import BRDFKernelWeights, GeometryAngles, RTCoefficients, SurfacePrior

        shape = (6, 6)
        solver = MultiGridSolver(
            MultiGridConfig(
                grid_search_aot_points=3,
                grid_search_tcwv_points=3,
                aot_bounds=(0.2, 0.8),
                tcwv_bounds=(1.0, 3.0),
                quadratic_block_size=3,
                quadratic_block_min_valid_fraction=0.5,
            )
        )

        toa_values = np.full(shape, 0.25, dtype=np.float32)
        toa_values[:3, :3] = -0.1
        toa_values[:2, :2] = 0.25
        toa = xr.DataArray(np.stack([toa_values]), dims=["band", "y", "x"])
        mask = xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"])
        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
        )
        atmo_prior = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.2), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.0), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        )
        brdf = BRDFKernelWeights(
            f0=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            f1=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            f2=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
            f0_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            f1_unc=xr.DataArray(np.full(shape, 0.005), dims=["y", "x"]),
            f2_unc=xr.DataArray(np.full(shape, 0.002), dims=["y", "x"]),
        )
        surface_prior = SurfacePrior(
            boa=xr.DataArray(np.stack([np.full(shape, 0.2)]), dims=["band", "y", "x"]),
            boa_unc=xr.DataArray(np.stack([np.full(shape, 0.02)]), dims=["band", "y", "x"]),
            kernels=brdf,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )

        class _DummyRT:
            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
                _ = (atmo_state, band, compute_jacobian)
                xap = xr.DataArray(
                    np.full(geometry.sza.shape, 0.9, dtype=np.float32), dims=["y", "x"]
                )
                return RTCoefficients(xap=xap, xbp=xr.zeros_like(xap), xcp=xr.zeros_like(xap))

            def supports_jacobian(self) -> bool:
                return False

            @property
            def backend_name(self) -> str:
                return "dummy"

            def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool:
                _ = (sensor_id, satellite_id)
                return True

        def _fail_pixel_eval(*_args, **_kwargs):  # noqa: ANN002, ANN003
            raise AssertionError(
                "Per-pixel evaluator should not be used when quadratic_block_size > 1"
            )

        def _fake_block_eval(
            _provider,
            aot_axis,
            tcwv_axis,
            _toa,
            _boa_prior,
            _boa_unc,
            _band_weights,
            valid_mask,
            _aot_prior,
            _tcwv_prior,
            _aot_prior_unc,
            _tcwv_prior_unc,
            block_size,
            fixed_parameter="none",
        ):  # noqa: ANN001
            assert block_size == 3
            assert fixed_parameter == "none"
            assert int(np.count_nonzero(valid_mask[:3, :3])) == 4
            costs = np.full((aot_axis.shape[0], tcwv_axis.shape[0], 2, 2), 50.0, dtype=np.float32)
            costs[1, 1, :, :] = 0.0
            obs_counts = np.ones_like(costs, dtype=np.uint16)
            return costs, obs_counts, np.ones((2, 2), dtype=bool)

        def _fake_refiner(
            costs,
            obs_counts,
            aot_axis,
            tcwv_axis,
            valid_mask,
            fixed_parameter="none",
        ):  # noqa: ANN001
            _ = (costs, obs_counts)
            assert fixed_parameter == "none"
            np.testing.assert_array_equal(
                valid_mask,
                np.array([[False, True], [True, True]], dtype=bool),
            )
            out_aot = np.full(valid_mask.shape, float(aot_axis[1]), dtype=np.float32)
            out_tcwv = np.full(valid_mask.shape, float(tcwv_axis[1]), dtype=np.float32)
            out_aot_unc = np.full(valid_mask.shape, 0.07, dtype=np.float32)
            out_tcwv_unc = np.full(valid_mask.shape, 0.21, dtype=np.float32)
            invalid = ~valid_mask
            boundary = np.zeros(valid_mask.shape, dtype=bool)
            lower_aot_boundary = np.zeros(valid_mask.shape, dtype=bool)
            zero_obs = ~valid_mask
            return (
                out_aot,
                out_tcwv,
                out_aot_unc,
                out_tcwv_unc,
                invalid,
                boundary,
                lower_aot_boundary,
                zero_obs,
            )

        monkeypatch.setattr(
            mg_mod, "evaluate_grid_search_cost_cube_with_provider_qa", _fail_pixel_eval
        )
        monkeypatch.setattr(
            mg_mod, "evaluate_block_grid_search_cost_cube_with_provider_qa", _fake_block_eval
        )
        monkeypatch.setattr(mg_mod, "quadratic_refine_grid_search_qa", _fake_refiner)

        _out_aot, _out_tcwv, _out_aot_unc, _out_tcwv_unc, level_diag = (
            solver._solve_level_grid_search(
                toa=toa,
                surface_prior=surface_prior,
                geometry=geometry,
                atmo_prior=atmo_prior,
                rt_model=_DummyRT(),
                mask=mask,
                bands=[SensorBand("B02", 490.0, 65.0, 10.0, 0)],
                cost_config=CostFunctionConfig(aot_gamma=0.0, tcwv_gamma=0.0),
            )
        )

        assert level_diag["valid_pixels"] == pytest.approx(31.0)
        assert level_diag["qa_invalid_pixels"] == pytest.approx(9.0)
        assert level_diag["qa_solve_invalid_pixels"] == pytest.approx(4.0)
        assert level_diag["qa_zero_obs_pixels"] == pytest.approx(9.0)
        assert level_diag["qa_insufficient_support_pixels"] == pytest.approx(9.0)
        assert np.all(level_diag["qa_insufficient_support_mask"][:3, :3])
        assert not np.any(level_diag["qa_insufficient_support_mask"][3:, :])
        assert not np.any(level_diag["qa_insufficient_support_mask"][:, 3:])

    def test_grid_search_uses_invalid_qa_to_restore_priors(self, monkeypatch):
        """Invalid QA pixels should be smoothed from trusted neighbouring retrievals."""
        from siac.algorithms.solver import multigrid as mg_mod
        from siac.domain import SensorBand
        from siac.runtime import BRDFKernelWeights, GeometryAngles, RTCoefficients, SurfacePrior

        shape = (3, 4)
        solver = MultiGridSolver(
            MultiGridConfig(grid_search_aot_points=3, grid_search_tcwv_points=3)
        )

        toa = xr.DataArray(
            np.stack([np.full(shape, 0.25, dtype=np.float32)]),
            dims=["band", "y", "x"],
        )
        mask = xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"])
        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
        )
        atmo_prior = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.2), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.0), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        )
        brdf = BRDFKernelWeights(
            f0=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            f1=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            f2=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
            f0_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            f1_unc=xr.DataArray(np.full(shape, 0.005), dims=["y", "x"]),
            f2_unc=xr.DataArray(np.full(shape, 0.002), dims=["y", "x"]),
        )
        surface_prior = SurfacePrior(
            boa=xr.DataArray(np.stack([np.full(shape, 0.2)]), dims=["band", "y", "x"]),
            boa_unc=xr.DataArray(np.stack([np.full(shape, 0.02)]), dims=["band", "y", "x"]),
            kernels=brdf,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )

        class _DummyRT:
            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
                _ = (band, compute_jacobian)
                factor = np.full(shape, 0.9, dtype=np.float32)
                xap = xr.DataArray(factor, dims=geometry.sza.dims, coords=geometry.sza.coords)
                return RTCoefficients(xap=xap, xbp=xr.zeros_like(xap), xcp=xr.zeros_like(xap))

            def supports_jacobian(self) -> bool:
                return False

            @property
            def backend_name(self) -> str:
                return "dummy"

            def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool:
                _ = (sensor_id, satellite_id)
                return True

        def _fake_rust_refiner(
            costs,
            obs_counts,
            aot_axis,
            tcwv_axis,
            valid_mask,
            fixed_parameter="none",
        ):  # noqa: ANN001
            _ = (costs, obs_counts, aot_axis, tcwv_axis)
            assert fixed_parameter == "none"
            out_aot = np.full(valid_mask.shape, 0.33, dtype=np.float32)
            out_tcwv = np.full(valid_mask.shape, 1.23, dtype=np.float32)
            out_aot_unc = np.full(valid_mask.shape, 0.07, dtype=np.float32)
            out_tcwv_unc = np.full(valid_mask.shape, 0.21, dtype=np.float32)
            invalid = np.zeros(valid_mask.shape, dtype=bool)
            invalid[0, 0] = True
            boundary = np.zeros(valid_mask.shape, dtype=bool)
            lower_aot_boundary = np.zeros(valid_mask.shape, dtype=bool)
            lower_aot_boundary[1, 1] = True
            zero_obs = np.zeros(valid_mask.shape, dtype=bool)
            return (
                out_aot,
                out_tcwv,
                out_aot_unc,
                out_tcwv_unc,
                invalid,
                boundary,
                lower_aot_boundary,
                zero_obs,
            )

        monkeypatch.setattr(mg_mod, "quadratic_refine_grid_search_qa", _fake_rust_refiner)

        out_aot, out_tcwv, out_aot_unc, out_tcwv_unc, diag = solver._solve_level_grid_search(
            toa=toa,
            surface_prior=surface_prior,
            geometry=geometry,
            atmo_prior=atmo_prior,
            rt_model=_DummyRT(),
            mask=mask,
            bands=[SensorBand("B02", 490.0, 65.0, 10.0, 0)],
            cost_config=CostFunctionConfig(),
        )

        # QA-invalid pixels are filled by spatial smoothing from neighbours,
        # not restored to the prior.  The smoothed value should be finite and
        # within a reasonable range, but not necessarily the exact prior.
        assert np.isfinite(out_aot[0, 0])
        assert np.isfinite(out_tcwv[0, 0])
        # Uncertainty at invalid pixels should be inflated.
        assert out_aot_unc[0, 0] >= 0.1
        assert out_tcwv_unc[0, 0] >= 0.5
        # Valid pixel values are smoothed but should remain close to original.
        assert out_aot[1, 1] == pytest.approx(0.33, abs=0.05)
        assert diag["qa_invalid_pixels"] == pytest.approx(1.0)
        assert diag["qa_lower_aot_boundary_pixels"] == pytest.approx(1.0)

    def test_rust_quadratic_refiner_matches_python_reference(self):
        """Rust quadratic refiner should match Python reference behavior."""
        rng = np.random.default_rng(123)

        n_aot, n_tcwv, ny, nx = 7, 6, 8, 9
        aot_axis = np.linspace(0.02, 1.2, n_aot, dtype=np.float32)
        tcwv_axis = np.linspace(0.1, 5.0, n_tcwv, dtype=np.float32)
        valid_mask = rng.random((ny, nx)) > 0.15

        a0 = rng.uniform(float(aot_axis[1]), float(aot_axis[-2]), size=(ny, nx))
        t0 = rng.uniform(float(tcwv_axis[1]), float(tcwv_axis[-2]), size=(ny, nx))
        wa = rng.uniform(0.5, 2.0, size=(ny, nx))
        wt = rng.uniform(0.3, 1.5, size=(ny, nx))
        wxy = rng.uniform(-0.15, 0.15, size=(ny, nx))
        bias = rng.uniform(0.0, 0.2, size=(ny, nx))

        costs = np.empty((n_aot, n_tcwv, ny, nx), dtype=np.float32)
        for ia, a in enumerate(aot_axis.astype(np.float64)):
            da = a - a0
            for it, t in enumerate(tcwv_axis.astype(np.float64)):
                dt = t - t0
                val = wa * da * da + wt * dt * dt + wxy * da * dt + bias
                costs[ia, it] = val.astype(np.float32)

        ref = _quadratic_refine_python_reference(costs, aot_axis, tcwv_axis, valid_mask)
        got = quadratic_refine_grid_search(costs, aot_axis, tcwv_axis, valid_mask)
        got = tuple(np.asarray(x, dtype=np.float32) for x in got)

        for expected, actual in zip(ref, got, strict=True):
            np.testing.assert_allclose(actual, expected, rtol=1e-5, atol=1e-6)

    def test_rust_quadratic_refiner_supports_fixed_tcwv_axis(self):
        """Fixed-TCWV mode should refine AOT using the 1-D local quadratic."""
        aot_axis = np.linspace(0.0, 1.0, 5, dtype=np.float32)
        tcwv_axis = np.array([2.0], dtype=np.float32)
        valid_mask = np.array([[True]], dtype=bool)
        target_aot = 0.4
        costs = np.empty((aot_axis.size, 1, 1, 1), dtype=np.float32)
        for ia, aot_val in enumerate(aot_axis):
            costs[ia, 0, 0, 0] = np.float32((float(aot_val) - target_aot) ** 2 + 0.25)
        obs_counts = np.ones_like(costs, dtype=np.uint16)

        (
            aot_best,
            tcwv_best,
            aot_unc,
            _tcwv_unc,
            invalid,
            boundary,
            lower_aot_boundary,
            zero_obs,
        ) = quadratic_refine_grid_search_qa(
            costs,
            obs_counts,
            aot_axis,
            tcwv_axis,
            valid_mask,
            "tcwv",
        )

        assert float(np.asarray(aot_best)[0, 0]) == pytest.approx(target_aot, abs=1.0e-6)
        assert float(np.asarray(tcwv_best)[0, 0]) == pytest.approx(2.0)
        assert np.isfinite(np.asarray(aot_unc)[0, 0])
        assert not np.asarray(invalid, dtype=bool)[0, 0]
        assert not np.asarray(boundary, dtype=bool)[0, 0]
        assert not np.asarray(lower_aot_boundary, dtype=bool)[0, 0]
        assert not np.asarray(zero_obs, dtype=bool)[0, 0]

    def test_rust_quadratic_refiner_supports_fixed_aot_axis(self):
        """Fixed-AOT mode should refine TCWV using the 1-D local quadratic."""
        aot_axis = np.array([0.2], dtype=np.float32)
        tcwv_axis = np.linspace(1.0, 5.0, 5, dtype=np.float32)
        valid_mask = np.array([[True]], dtype=bool)
        target_tcwv = 2.8
        costs = np.empty((1, tcwv_axis.size, 1, 1), dtype=np.float32)
        for it, tcwv_val in enumerate(tcwv_axis):
            costs[0, it, 0, 0] = np.float32((float(tcwv_val) - target_tcwv) ** 2 + 0.25)
        obs_counts = np.ones_like(costs, dtype=np.uint16)

        (
            aot_best,
            tcwv_best,
            _aot_unc,
            tcwv_unc,
            invalid,
            boundary,
            lower_aot_boundary,
            zero_obs,
        ) = quadratic_refine_grid_search_qa(
            costs,
            obs_counts,
            aot_axis,
            tcwv_axis,
            valid_mask,
            "aot",
        )

        assert float(np.asarray(aot_best)[0, 0]) == pytest.approx(0.2)
        assert float(np.asarray(tcwv_best)[0, 0]) == pytest.approx(target_tcwv, abs=1.0e-6)
        assert np.isfinite(np.asarray(tcwv_unc)[0, 0])
        assert not np.asarray(invalid, dtype=bool)[0, 0]
        assert not np.asarray(boundary, dtype=bool)[0, 0]
        assert not np.asarray(lower_aot_boundary, dtype=bool)[0, 0]
        assert not np.asarray(zero_obs, dtype=bool)[0, 0]

    def test_rust_grid_search_qa_distinguishes_prior_only_floor_from_true_boundary(self):
        """QA refiner should separate unsupported floor hits from supported boundary minima."""
        aot_axis = np.linspace(0.001, 2.5, 11, dtype=np.float32)
        tcwv_axis = np.array([1.0, 2.0, 3.0], dtype=np.float32)
        valid_mask = np.array([[True, True]], dtype=bool)

        costs = np.zeros((len(aot_axis), len(tcwv_axis), 1, 2), dtype=np.float32)
        obs_counts = np.zeros_like(costs, dtype=np.uint16)

        prior_aot = 0.05
        prior_tcwv = 2.0
        for ia, aot_val in enumerate(aot_axis):
            for it, tcwv_val in enumerate(tcwv_axis):
                costs[ia, it, 0, 0] = np.float32(
                    0.5 * ((float(aot_val) - prior_aot) ** 2 + (float(tcwv_val) - prior_tcwv) ** 2)
                )
                costs[ia, it, 0, 1] = np.float32(float(ia) + 0.1 * (float(tcwv_val) - 2.0) ** 2)
                obs_counts[ia, it, 0, 1] = 2

        (
            aot_best,
            tcwv_best,
            _aot_unc,
            _tcwv_unc,
            invalid_mask,
            boundary_mask,
            lower_aot_boundary_mask,
            zero_obs_mask,
        ) = quadratic_refine_grid_search_qa(
            costs,
            obs_counts,
            aot_axis,
            tcwv_axis,
            valid_mask,
        )

        aot_best = np.asarray(aot_best, dtype=np.float32)
        tcwv_best = np.asarray(tcwv_best, dtype=np.float32)
        invalid_mask = np.asarray(invalid_mask, dtype=bool)
        boundary_mask = np.asarray(boundary_mask, dtype=bool)
        lower_aot_boundary_mask = np.asarray(lower_aot_boundary_mask, dtype=bool)
        zero_obs_mask = np.asarray(zero_obs_mask, dtype=bool)

        assert aot_best[0, 0] == pytest.approx(aot_axis[0])
        assert tcwv_best[0, 0] == pytest.approx(prior_tcwv)
        assert invalid_mask[0, 0]
        assert zero_obs_mask[0, 0]
        assert not boundary_mask[0, 0]
        assert not lower_aot_boundary_mask[0, 0]

        assert aot_best[0, 1] == pytest.approx(aot_axis[0])
        assert tcwv_best[0, 1] == pytest.approx(2.0)
        assert not invalid_mask[0, 1]
        assert not zero_obs_mask[0, 1]
        assert boundary_mask[0, 1]
        assert lower_aot_boundary_mask[0, 1]

    def test_rust_cost_cube_provider_qa_reproduces_prior_only_floor_case(self):
        """Zero observation support plus a low prior should collapse to the AOT floor."""
        aot_axis = np.linspace(0.001, 2.5, 11, dtype=np.float32)
        tcwv_axis = np.array([1.0, 2.0, 3.0], dtype=np.float32)
        valid_mask = np.array([[True]], dtype=bool)

        toa = np.array([[[0.2]]], dtype=np.float32)
        boa_prior = np.array([[[0.1]]], dtype=np.float32)
        boa_unc = np.array([[[np.nan]]], dtype=np.float32)
        band_weights = np.array([1.0], dtype=np.float32)
        aot_prior = np.array([[0.05]], dtype=np.float32)
        tcwv_prior = np.array([[2.0]], dtype=np.float32)
        aot_prior_unc = np.array([[0.05]], dtype=np.float32)
        tcwv_prior_unc = np.array([[0.5]], dtype=np.float32)

        def _provider(
            _aot_val: float, _tcwv_val: float
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            xap = np.array([[[1.0]]], dtype=np.float32)
            zeros = np.array([[[0.0]]], dtype=np.float32)
            return xap, zeros, zeros

        costs, obs_counts = evaluate_grid_search_cost_cube_with_provider_qa(
            _provider,
            aot_axis,
            tcwv_axis,
            toa,
            boa_prior,
            boa_unc,
            band_weights,
            valid_mask,
            aot_prior,
            tcwv_prior,
            aot_prior_unc,
            tcwv_prior_unc,
        )
        costs = np.asarray(costs, dtype=np.float32)
        obs_counts = np.asarray(obs_counts, dtype=np.uint16)

        assert np.all(obs_counts == 0)

        old_aot, _old_tcwv, _old_aot_unc, _old_tcwv_unc = quadratic_refine_grid_search(
            costs,
            aot_axis,
            tcwv_axis,
            valid_mask,
        )
        assert np.asarray(old_aot, dtype=np.float32)[0, 0] == pytest.approx(aot_axis[0])

        (
            _aot_best,
            _tcwv_best,
            _aot_unc,
            _tcwv_unc,
            invalid_mask,
            boundary_mask,
            lower_aot_boundary_mask,
            zero_obs_mask,
        ) = quadratic_refine_grid_search_qa(
            costs,
            obs_counts,
            aot_axis,
            tcwv_axis,
            valid_mask,
        )

        invalid_mask = np.asarray(invalid_mask, dtype=bool)
        boundary_mask = np.asarray(boundary_mask, dtype=bool)
        lower_aot_boundary_mask = np.asarray(lower_aot_boundary_mask, dtype=bool)
        zero_obs_mask = np.asarray(zero_obs_mask, dtype=bool)

        assert invalid_mask[0, 0]
        assert zero_obs_mask[0, 0]
        assert not boundary_mask[0, 0]
        assert not lower_aot_boundary_mask[0, 0]

    def test_rust_cost_cube_provider_qa_does_not_count_zero_weight_support(self):
        """Zero-weight bands should not be treated as observation support."""
        aot_axis = np.linspace(0.001, 2.5, 11, dtype=np.float32)
        tcwv_axis = np.array([1.0, 2.0, 3.0], dtype=np.float32)
        valid_mask = np.array([[True]], dtype=bool)

        toa = np.array([[[0.2]]], dtype=np.float32)
        boa_prior = np.array([[[0.1]]], dtype=np.float32)
        boa_unc = np.array([[[0.02]]], dtype=np.float32)
        band_weights = np.array([0.0], dtype=np.float32)
        aot_prior = np.array([[0.05]], dtype=np.float32)
        tcwv_prior = np.array([[2.0]], dtype=np.float32)
        aot_prior_unc = np.array([[0.05]], dtype=np.float32)
        tcwv_prior_unc = np.array([[0.5]], dtype=np.float32)

        def _provider(
            _aot_val: float, _tcwv_val: float
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            xap = np.array([[[1.0]]], dtype=np.float32)
            zeros = np.array([[[0.0]]], dtype=np.float32)
            return xap, zeros, zeros

        _costs, obs_counts = evaluate_grid_search_cost_cube_with_provider_qa(
            _provider,
            aot_axis,
            tcwv_axis,
            toa,
            boa_prior,
            boa_unc,
            band_weights,
            valid_mask,
            aot_prior,
            tcwv_prior,
            aot_prior_unc,
            tcwv_prior_unc,
        )
        obs_counts = np.asarray(obs_counts, dtype=np.uint16)

        assert np.all(obs_counts == 0)

    def test_high_aot_second_pass_adopts_only_strict_lower_upward_moves(self):
        """Fine high-AOD pass should switch only where it improves and moves upward."""
        near_top = np.array(
            [
                [True, True, False],
                [True, True, True],
            ],
            dtype=bool,
        )
        pass1_costs = np.full((2, 1, 2, 3), 10.0, dtype=np.float32)
        pass2_costs = np.full_like(pass1_costs, 11.0)
        pass2_costs[0, 0, 0, 0] = 8.0  # switch: near top, lower cost, upward
        pass2_costs[0, 0, 0, 1] = 8.0  # no switch: second pass moves down
        pass2_costs[0, 0, 0, 2] = 8.0  # no switch: primary pixel was not high
        pass2_costs[0, 0, 1, 0] = 10.0  # no switch: not strictly lower
        pass2_costs[0, 0, 1, 1] = 7.0  # switch
        pass2_costs[0, 0, 1, 2] = np.float32(10.0 - 5.0e-7)  # below tolerance

        aot1 = np.full((2, 3), 0.9, dtype=np.float32)
        tcwv1 = np.full((2, 3), 2.0, dtype=np.float32)
        aot_unc1 = np.full((2, 3), 0.1, dtype=np.float32)
        tcwv_unc1 = np.full((2, 3), 0.2, dtype=np.float32)
        aot2 = np.array(
            [
                [1.2, 0.8, 1.3],
                [1.4, 1.5, 1.6],
            ],
            dtype=np.float32,
        )
        tcwv2 = np.full((2, 3), 3.0, dtype=np.float32)
        aot_unc2 = np.full((2, 3), 0.3, dtype=np.float32)
        tcwv_unc2 = np.full((2, 3), 0.4, dtype=np.float32)

        out_aot, out_tcwv, out_aot_unc, out_tcwv_unc = MultiGridSolver._adopt_lower_cost_high_aot(
            near_top=near_top,
            block_size=1,
            shape=near_top.shape,
            pass1_costs=pass1_costs,
            pass2_costs=pass2_costs,
            pass1=(aot1, tcwv1, aot_unc1, tcwv_unc1),
            pass2=(aot2, tcwv2, aot_unc2, tcwv_unc2),
        )

        expected_switch = np.array(
            [
                [True, False, False],
                [False, True, False],
            ],
            dtype=bool,
        )
        np.testing.assert_allclose(out_aot, np.where(expected_switch, aot2, aot1))
        np.testing.assert_allclose(out_tcwv, np.where(expected_switch, tcwv2, tcwv1))
        np.testing.assert_allclose(out_aot_unc, np.where(expected_switch, aot_unc2, aot_unc1))
        np.testing.assert_allclose(out_tcwv_unc, np.where(expected_switch, tcwv_unc2, tcwv_unc1))

    def test_high_aot_second_pass_broadcasts_block_cost_improvement(self):
        """Block-grid second-pass improvements should map back to pixel space."""
        near_top = np.array(
            [
                [False, False, True],
                [False, False, True],
                [True, False, True],
            ],
            dtype=bool,
        )
        pass1_costs = np.full((2, 1, 2, 2), 10.0, dtype=np.float32)
        pass2_costs = np.full_like(pass1_costs, 11.0)
        pass2_costs[0, 0, 0, 1] = 8.0
        pass2_costs[0, 0, 1, 0] = 7.0

        aot1 = np.full((3, 3), 0.9, dtype=np.float32)
        tcwv1 = np.full((3, 3), 2.0, dtype=np.float32)
        aot_unc1 = np.full((3, 3), 0.1, dtype=np.float32)
        tcwv_unc1 = np.full((3, 3), 0.2, dtype=np.float32)
        aot2 = np.full((3, 3), 1.2, dtype=np.float32)
        tcwv2 = np.full((3, 3), 3.0, dtype=np.float32)
        aot_unc2 = np.full((3, 3), 0.3, dtype=np.float32)
        tcwv_unc2 = np.full((3, 3), 0.4, dtype=np.float32)

        out_aot, out_tcwv, out_aot_unc, out_tcwv_unc = MultiGridSolver._adopt_lower_cost_high_aot(
            near_top=near_top,
            block_size=2,
            shape=near_top.shape,
            pass1_costs=pass1_costs,
            pass2_costs=pass2_costs,
            pass1=(aot1, tcwv1, aot_unc1, tcwv_unc1),
            pass2=(aot2, tcwv2, aot_unc2, tcwv_unc2),
        )

        expected_switch = np.array(
            [
                [False, False, True],
                [False, False, True],
                [True, False, False],
            ],
            dtype=bool,
        )
        np.testing.assert_allclose(out_aot, np.where(expected_switch, aot2, aot1))
        np.testing.assert_allclose(out_tcwv, np.where(expected_switch, tcwv2, tcwv1))
        np.testing.assert_allclose(out_aot_unc, np.where(expected_switch, aot_unc2, aot_unc1))
        np.testing.assert_allclose(out_tcwv_unc, np.where(expected_switch, tcwv_unc2, tcwv_unc1))

    def test_grid_search_level_diag_reports_prior_only_floor_scenario(self):
        """End-to-end grid search should report unsupported floor hits and restore those pixels to the prior."""
        from siac.domain import SensorBand
        from siac.runtime import BRDFKernelWeights, GeometryAngles, RTCoefficients, SurfacePrior

        shape = (2, 2)
        solver = MultiGridSolver(
            MultiGridConfig(grid_search_aot_points=11, grid_search_tcwv_points=3)
        )

        toa = xr.DataArray(
            np.stack([np.full(shape, 0.2, dtype=np.float32)]),
            dims=["band", "y", "x"],
        )
        mask = xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"])
        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5, dtype=np.float32), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5, dtype=np.float32), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5, dtype=np.float32), dims=["y", "x"]),
        )
        atmo_prior = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.0, dtype=np.float32), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
        )
        brdf = BRDFKernelWeights(
            f0=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
            f1=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=["y", "x"]),
            f2=xr.DataArray(np.full(shape, 0.02, dtype=np.float32), dims=["y", "x"]),
            f0_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=["y", "x"]),
            f1_unc=xr.DataArray(np.full(shape, 0.005, dtype=np.float32), dims=["y", "x"]),
            f2_unc=xr.DataArray(np.full(shape, 0.002, dtype=np.float32), dims=["y", "x"]),
        )
        surface_prior = SurfacePrior(
            boa=xr.DataArray(
                np.stack([np.full(shape, 0.2, dtype=np.float32)]), dims=["band", "y", "x"]
            ),
            boa_unc=xr.DataArray(
                np.stack([np.full(shape, np.nan, dtype=np.float32)]), dims=["band", "y", "x"]
            ),
            kernels=brdf,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )

        class _DummyRT:
            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
                _ = (geometry, atmo_state, band, compute_jacobian)
                nan = xr.DataArray(np.full(shape, np.nan, dtype=np.float32), dims=["y", "x"])
                return RTCoefficients(xap=nan, xbp=nan.copy(), xcp=nan.copy())

            def supports_jacobian(self) -> bool:
                return False

            @property
            def backend_name(self) -> str:
                return "dummy"

            def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool:
                _ = (sensor_id, satellite_id)
                return True

        out_aot, _out_tcwv, _out_aot_unc, _out_tcwv_unc, level_diag = (
            solver._solve_level_grid_search(
                toa=toa,
                surface_prior=surface_prior,
                geometry=geometry,
                atmo_prior=atmo_prior,
                rt_model=_DummyRT(),
                mask=mask,
                bands=[SensorBand("B02", 490.0, 65.0, 10.0, 0)],
                cost_config=CostFunctionConfig(),
            )
        )

        np.testing.assert_allclose(out_aot, np.full(shape, 0.05, dtype=np.float32))
        assert level_diag["qa_invalid_pixels"] == pytest.approx(float(np.prod(shape)))
        assert level_diag["qa_solve_invalid_pixels"] == pytest.approx(0.0)
        assert level_diag["qa_zero_obs_pixels"] == pytest.approx(float(np.prod(shape)))
        assert level_diag["qa_invalid_floor_pixels"] == pytest.approx(0.0)
        assert level_diag["qa_prior_floor_pixels"] == pytest.approx(0.0)
        assert level_diag["qa_boundary_pixels"] == pytest.approx(0.0)

    def test_grid_search_level_diag_distinguishes_prior_floor_passthrough(self, monkeypatch):
        """Final floor-valued AOT can come from prior restoration rather than a raw floor hit."""
        from siac.algorithms.solver import multigrid as mg_mod
        from siac.domain import SensorBand
        from siac.runtime import BRDFKernelWeights, GeometryAngles, RTCoefficients, SurfacePrior

        shape = (2, 2)
        solver = MultiGridSolver(
            MultiGridConfig(grid_search_aot_points=3, grid_search_tcwv_points=3)
        )

        toa = xr.DataArray(
            np.stack([np.full(shape, 0.25, dtype=np.float32)]),
            dims=["band", "y", "x"],
        )
        mask = xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"])
        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
        )
        atmo_prior = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.001, dtype=np.float32), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.0, dtype=np.float32), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
        )
        brdf = BRDFKernelWeights(
            f0=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            f1=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            f2=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
            f0_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            f1_unc=xr.DataArray(np.full(shape, 0.005), dims=["y", "x"]),
            f2_unc=xr.DataArray(np.full(shape, 0.002), dims=["y", "x"]),
        )
        surface_prior = SurfacePrior(
            boa=xr.DataArray(np.stack([np.full(shape, 0.2)]), dims=["band", "y", "x"]),
            boa_unc=xr.DataArray(np.stack([np.full(shape, 0.02)]), dims=["band", "y", "x"]),
            kernels=brdf,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )

        class _DummyRT:
            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
                _ = (band, compute_jacobian)
                factor = np.full(shape, 0.9, dtype=np.float32)
                xap = xr.DataArray(factor, dims=geometry.sza.dims, coords=geometry.sza.coords)
                return RTCoefficients(xap=xap, xbp=xr.zeros_like(xap), xcp=xr.zeros_like(xap))

            def supports_jacobian(self) -> bool:
                return False

            @property
            def backend_name(self) -> str:
                return "dummy"

            def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool:
                _ = (sensor_id, satellite_id)
                return True

        def _fake_rust_refiner(
            costs,
            obs_counts,
            aot_axis,
            tcwv_axis,
            valid_mask,
            fixed_parameter="none",
        ):  # noqa: ANN001
            _ = (costs, obs_counts, aot_axis, tcwv_axis)
            assert fixed_parameter == "none"
            out_aot = np.full(valid_mask.shape, 0.33, dtype=np.float32)
            out_tcwv = np.full(valid_mask.shape, 1.23, dtype=np.float32)
            out_aot_unc = np.full(valid_mask.shape, 0.07, dtype=np.float32)
            out_tcwv_unc = np.full(valid_mask.shape, 0.21, dtype=np.float32)
            invalid = np.ones(valid_mask.shape, dtype=bool)
            boundary = np.zeros(valid_mask.shape, dtype=bool)
            lower_aot_boundary = np.zeros(valid_mask.shape, dtype=bool)
            zero_obs = np.ones(valid_mask.shape, dtype=bool)
            return (
                out_aot,
                out_tcwv,
                out_aot_unc,
                out_tcwv_unc,
                invalid,
                boundary,
                lower_aot_boundary,
                zero_obs,
            )

        monkeypatch.setattr(mg_mod, "quadratic_refine_grid_search_qa", _fake_rust_refiner)

        out_aot, _out_tcwv, _out_aot_unc, _out_tcwv_unc, level_diag = (
            solver._solve_level_grid_search(
                toa=toa,
                surface_prior=surface_prior,
                geometry=geometry,
                atmo_prior=atmo_prior,
                rt_model=_DummyRT(),
                mask=mask,
                bands=[SensorBand("B02", 490.0, 65.0, 10.0, 0)],
                cost_config=CostFunctionConfig(),
            )
        )

        np.testing.assert_allclose(out_aot, np.full(shape, 0.001, dtype=np.float32))
        assert level_diag["qa_invalid_pixels"] == pytest.approx(float(np.prod(shape)))
        assert level_diag["qa_invalid_floor_pixels"] == pytest.approx(0.0)
        assert level_diag["qa_prior_floor_pixels"] == pytest.approx(float(np.prod(shape)))

    def test_rust_candidate_cost_matches_python_reference(self):
        """Rust candidate cost helper should match Python reference."""
        rng = np.random.default_rng(7)

        n_band, ny, nx = 3, 6, 7
        toa = rng.uniform(0.02, 0.4, size=(n_band, ny, nx)).astype(np.float32)
        xap = rng.uniform(0.7, 1.3, size=(n_band, ny, nx)).astype(np.float32)
        xbp = rng.uniform(0.0, 0.1, size=(n_band, ny, nx)).astype(np.float32)
        xcp = rng.uniform(-0.15, 0.15, size=(n_band, ny, nx)).astype(np.float32)
        boa_prior = rng.uniform(0.02, 0.5, size=(n_band, ny, nx)).astype(np.float32)
        boa_unc = rng.uniform(0.01, 0.08, size=(n_band, ny, nx)).astype(np.float32)
        band_weights = np.array([0.2, 0.3, 0.5], dtype=np.float32)
        valid_mask = rng.random((ny, nx)) > 0.2

        aot_val = np.float32(0.24)
        tcwv_val = np.float32(2.1)
        aot_prior = rng.uniform(0.05, 0.6, size=(ny, nx)).astype(np.float32)
        tcwv_prior = rng.uniform(0.8, 4.0, size=(ny, nx)).astype(np.float32)
        aot_prior_unc = rng.uniform(0.02, 0.15, size=(ny, nx)).astype(np.float32)
        tcwv_prior_unc = rng.uniform(0.1, 0.8, size=(ny, nx)).astype(np.float32)

        expected = _candidate_cost_python_reference(
            toa=toa,
            xap=xap,
            xbp=xbp,
            xcp=xcp,
            boa_prior=boa_prior,
            boa_unc=boa_unc,
            band_weights=band_weights,
            valid_mask=valid_mask,
            aot_val=float(aot_val),
            tcwv_val=float(tcwv_val),
            aot_prior=aot_prior,
            tcwv_prior=tcwv_prior,
            aot_prior_unc=aot_prior_unc,
            tcwv_prior_unc=tcwv_prior_unc,
        )
        got = evaluate_grid_search_candidate_cost(
            toa,
            xap,
            xbp,
            xcp,
            boa_prior,
            boa_unc,
            band_weights,
            valid_mask,
            aot_val,
            tcwv_val,
            aot_prior,
            tcwv_prior,
            aot_prior_unc,
            tcwv_prior_unc,
        )
        np.testing.assert_allclose(
            np.asarray(got, dtype=np.float32),
            expected,
            rtol=1e-5,
            atol=1e-6,
        )

    def test_rust_cost_cube_provider_matches_python_reference(self):
        """Rust cost-cube provider loop should match Python candidate-loop reference."""
        rng = np.random.default_rng(11)

        n_aot, n_tcwv, n_band, ny, nx = 4, 5, 2, 3, 4
        aot_axis = np.linspace(0.05, 0.8, n_aot, dtype=np.float32)
        tcwv_axis = np.linspace(0.2, 4.2, n_tcwv, dtype=np.float32)

        toa = rng.uniform(0.03, 0.45, size=(n_band, ny, nx)).astype(np.float32)
        boa_prior = rng.uniform(0.02, 0.5, size=(n_band, ny, nx)).astype(np.float32)
        boa_unc = rng.uniform(0.01, 0.09, size=(n_band, ny, nx)).astype(np.float32)
        band_weights = np.array([0.4, 0.6], dtype=np.float32)
        valid_mask = rng.random((ny, nx)) > 0.15
        aot_prior = rng.uniform(0.04, 0.6, size=(ny, nx)).astype(np.float32)
        tcwv_prior = rng.uniform(0.5, 3.8, size=(ny, nx)).astype(np.float32)
        aot_prior_unc = rng.uniform(0.02, 0.15, size=(ny, nx)).astype(np.float32)
        tcwv_prior_unc = rng.uniform(0.1, 0.9, size=(ny, nx)).astype(np.float32)

        xap_cube = np.empty((n_aot, n_tcwv, n_band, ny, nx), dtype=np.float32)
        xbp_cube = np.empty_like(xap_cube)
        xcp_cube = np.empty_like(xap_cube)
        base = rng.uniform(0.8, 1.2, size=(n_band, ny, nx)).astype(np.float32)
        for ia, aot_val in enumerate(aot_axis):
            for it, tcwv_val in enumerate(tcwv_axis):
                xap_cube[ia, it] = base + 0.15 * float(aot_val) + 0.03 * float(tcwv_val)
                xbp_cube[ia, it] = 0.02 + 0.01 * float(aot_val)
                xcp_cube[ia, it] = -0.05 + 0.02 * float(tcwv_val) / 4.2

        expected = np.empty((n_aot, n_tcwv, ny, nx), dtype=np.float32)
        for ia, aot_val in enumerate(aot_axis):
            for it, tcwv_val in enumerate(tcwv_axis):
                expected[ia, it] = _candidate_cost_python_reference(
                    toa=toa,
                    xap=xap_cube[ia, it],
                    xbp=xbp_cube[ia, it],
                    xcp=xcp_cube[ia, it],
                    boa_prior=boa_prior,
                    boa_unc=boa_unc,
                    band_weights=band_weights,
                    valid_mask=valid_mask,
                    aot_val=float(aot_val),
                    tcwv_val=float(tcwv_val),
                    aot_prior=aot_prior,
                    tcwv_prior=tcwv_prior,
                    aot_prior_unc=aot_prior_unc,
                    tcwv_prior_unc=tcwv_prior_unc,
                )

        def _provider(aot_val: float, tcwv_val: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            ia = int(np.argmin(np.abs(aot_axis - np.float32(aot_val))))
            it = int(np.argmin(np.abs(tcwv_axis - np.float32(tcwv_val))))
            return xap_cube[ia, it], xbp_cube[ia, it], xcp_cube[ia, it]

        got = evaluate_grid_search_cost_cube_with_provider(
            _provider,
            aot_axis,
            tcwv_axis,
            toa,
            boa_prior,
            boa_unc,
            band_weights,
            valid_mask,
            aot_prior,
            tcwv_prior,
            aot_prior_unc,
            tcwv_prior_unc,
        )
        np.testing.assert_allclose(
            np.asarray(got, dtype=np.float32),
            expected,
            rtol=1e-5,
            atol=1e-6,
        )

    def test_grid_search_accepts_surface_prior_with_extra_bands(self):
        """Grid-search should handle surface priors with more bands than solver bands."""
        from siac.domain import SensorBand
        from siac.runtime import BRDFKernelWeights, GeometryAngles, RTCoefficients, SurfacePrior

        shape = (5, 6)
        solver = MultiGridSolver(
            MultiGridConfig(grid_search_aot_points=3, grid_search_tcwv_points=3)
        )

        toa = xr.DataArray(
            np.stack([np.full(shape, 0.24, dtype=np.float32)]),
            dims=["band", "y", "x"],
        )
        mask = xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"])
        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
        )
        atmo_prior = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.2), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.0), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        )
        brdf = BRDFKernelWeights(
            f0=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            f1=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            f2=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
            f0_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            f1_unc=xr.DataArray(np.full(shape, 0.005), dims=["y", "x"]),
            f2_unc=xr.DataArray(np.full(shape, 0.002), dims=["y", "x"]),
        )
        surface_prior = SurfacePrior(
            boa=xr.DataArray(
                np.stack(
                    [
                        np.full(shape, 0.20, dtype=np.float32),
                        np.full(shape, 0.21, dtype=np.float32),
                        np.full(shape, 0.22, dtype=np.float32),
                    ]
                ),
                dims=["band", "y", "x"],
            ),
            boa_unc=xr.DataArray(
                np.stack(
                    [
                        np.full(shape, 0.02, dtype=np.float32),
                        np.full(shape, 0.02, dtype=np.float32),
                        np.full(shape, 0.02, dtype=np.float32),
                    ]
                ),
                dims=["band", "y", "x"],
            ),
            kernels=brdf,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )

        class _DummyRT:
            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
                _ = (band, compute_jacobian)
                factor = np.full(shape, 0.9, dtype=np.float32)
                xap = xr.DataArray(factor, dims=geometry.sza.dims, coords=geometry.sza.coords)
                return RTCoefficients(xap=xap, xbp=xr.zeros_like(xap), xcp=xr.zeros_like(xap))

            def supports_jacobian(self) -> bool:
                return False

            @property
            def backend_name(self) -> str:
                return "dummy"

            def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool:
                _ = (sensor_id, satellite_id)
                return True

        out_aot, out_tcwv, out_aot_unc, out_tcwv_unc, _ = solver._solve_level_grid_search(
            toa=toa,
            surface_prior=surface_prior,
            geometry=geometry,
            atmo_prior=atmo_prior,
            rt_model=_DummyRT(),
            mask=mask,
            bands=[SensorBand("B02", 490.0, 65.0, 10.0, 0)],
            cost_config=CostFunctionConfig(),
        )

        assert np.isfinite(out_aot).all()
        assert np.isfinite(out_tcwv).all()
        assert np.isfinite(out_aot_unc).all()
        assert np.isfinite(out_tcwv_unc).all()


class TestProtocolValidation:
    """Tests for protocol validation at boundaries."""

    def test_cost_function_invalid_rt_model(self, cost_function_inputs):
        """CostFunction should reject non-RTModelBackend objects."""

        with pytest.raises(TypeError, match="rt_model must implement RTModelBackend"):
            CostFunction(
                toa=cost_function_inputs["toa"],
                surface_prior=cost_function_inputs["surface_prior"],
                geometry=cost_function_inputs["geometry"],
                atmo_prior=cost_function_inputs["atmo_prior"],
                rt_model=object(),
                bands=cost_function_inputs["bands"],
                mask=cost_function_inputs["mask"],
            )

    def test_multigrid_solver_invalid_rt_model(self):
        """MultiGridSolver.solve() should reject non-RTModelBackend objects."""
        solver = MultiGridSolver()

        shape = (16, 16)
        toa = xr.DataArray(np.full((3, *shape), 0.2), dims=["band", "y", "x"])
        cloud_mask = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"])

        from siac.domain import SensorBand
        from siac.runtime import BRDFKernelWeights, GeometryAngles, SurfacePrior

        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
        )

        atmo_prior = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        )

        brdf = BRDFKernelWeights(
            f0=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            f1=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            f2=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
            f0_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            f1_unc=xr.DataArray(np.full(shape, 0.005), dims=["y", "x"]),
            f2_unc=xr.DataArray(np.full(shape, 0.002), dims=["y", "x"]),
        )

        surface_prior = SurfacePrior(
            boa=xr.DataArray(np.full(shape, 0.12), dims=["y", "x"]),
            boa_unc=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
            kernels=brdf,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )

        bands = [SensorBand("B02", 490.0, 65.0, 10.0, 0)]

        with pytest.raises(TypeError, match="rt_model must implement RTModelBackend"):
            solver.solve(toa, surface_prior, geometry, atmo_prior, object(), cloud_mask, bands)


class TestResampleFieldEquivalence:
    """Tests for resampling field consistency."""

    def test_resample_upsample(self):
        """Upsampling should produce larger array with interpolated values."""
        solver = MultiGridSolver()
        field = np.array([[1.0, 2.0], [3.0, 4.0]])
        result = solver._resample_field(field, (4, 4))

        assert result.shape == (4, 4)
        # Corners should be close to original values
        assert abs(result[0, 0] - 1.0) < 0.5
        assert abs(result[-1, -1] - 4.0) < 0.5

    def test_resample_downsample(self):
        """Downsampling should produce smaller array with averaged values."""
        solver = MultiGridSolver()
        field = np.ones((8, 8)) * 5.0
        result = solver._resample_field(field, (4, 4))

        assert result.shape == (4, 4)
        np.testing.assert_allclose(result, 5.0, rtol=1e-6)

    def test_resample_roundtrip(self):
        """Downsample then upsample should approximately preserve smooth fields."""
        solver = MultiGridSolver()
        # Create smooth field
        y, x = np.mgrid[:32, :32]
        field = 0.5 + 0.1 * np.sin(2 * np.pi * x / 32) + 0.1 * np.cos(2 * np.pi * y / 32)

        small = solver._resample_field(field, (8, 8))
        restored = solver._resample_field(small, (32, 32))

        # Mean should be preserved
        np.testing.assert_allclose(field.mean(), restored.mean(), rtol=0.1)


class TestCostFunctionGradient:
    """Tests for cost function gradient accuracy."""

    @pytest.fixture
    def simple_cost_inputs(self):
        """Create simple inputs for cost function testing."""
        shape = (16, 16)

        # Simple geometry
        sza = xr.DataArray(np.full(shape, 0.5), dims=["y", "x"])
        saa = xr.DataArray(np.full(shape, 2.5), dims=["y", "x"])
        vza = xr.DataArray(np.full(shape, 0.1), dims=["y", "x"])
        vaa = xr.DataArray(np.full(shape, 1.5), dims=["y", "x"])

        from siac.runtime import GeometryAngles

        geometry = GeometryAngles(sza=sza, saa=saa, vza=vza, vaa=vaa)

        # Atmospheric state
        atmo_prior = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        )

        return geometry, atmo_prior, shape

    def test_smoothness_cost_gradient(self, simple_cost_inputs):
        """Pseudo-Huber smoothness gradient should match numerical gradient."""
        _geometry, _atmo_prior, shape = simple_cost_inputs

        gamma = 5.0
        delta = 0.02

        aot = np.random.default_rng(4).standard_normal(shape) * 0.1 + 0.15

        # Analytical gradient via _pseudo_huber_cost_grad
        j_analytical, dj_analytical = CostFunction._pseudo_huber_cost_grad(
            aot,
            gamma,
            delta,
        )

        # Numerical gradient
        eps = 1e-6
        numerical_grad = np.zeros_like(aot)

        for i in range(4):  # Test a few interior and boundary points
            for j in range(4):
                aot_plus = aot.copy()
                aot_plus[i, j] += eps
                j_plus, _ = CostFunction._pseudo_huber_cost_grad(
                    aot_plus,
                    gamma,
                    delta,
                )

                aot_minus = aot.copy()
                aot_minus[i, j] -= eps
                j_minus, _ = CostFunction._pseudo_huber_cost_grad(
                    aot_minus,
                    gamma,
                    delta,
                )

                numerical_grad[i, j] = (j_plus - j_minus) / (2 * eps)

        # Compare at tested points
        np.testing.assert_allclose(
            dj_analytical[:4, :4], numerical_grad[:4, :4], rtol=1e-4, atol=1e-8
        )

    def test_smoothness_hotspot_preservation(self, simple_cost_inputs):
        """Pseudo-Huber penalty on large gradients should be sub-quadratic."""
        _geometry, _atmo_prior, shape = simple_cost_inputs

        gamma = 5.0
        delta = 0.02

        # Smooth field
        smooth = np.full(shape, 0.15)
        j_smooth, _ = CostFunction._pseudo_huber_cost_grad(smooth, gamma, delta)

        # Field with a hotspot (large local gradient)
        hotspot = smooth.copy()
        hotspot[shape[0] // 2, shape[1] // 2] = 0.8

        j_hotspot, _ = CostFunction._pseudo_huber_cost_grad(hotspot, gamma, delta)

        # For comparison: what L2 would give (quadratic in gradient)
        dy_h = np.diff(hotspot, axis=0)
        dx_h = np.diff(hotspot, axis=1)
        j_l2 = 0.5 * gamma**2 * (np.sum(dy_h**2) + np.sum(dx_h**2))

        # Pseudo-Huber cost should be strictly less than L2 for large gradients
        assert j_hotspot < j_l2
        # And both should be > smooth cost
        assert j_hotspot > j_smooth

    def test_prior_cost_ignores_masked_pixels(self, cost_function_inputs):
        class _DummyRT:
            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):
                del geometry, atmo_state, band, compute_jacobian
                raise AssertionError("RT coefficients are not used in prior_cost_only")

            def supports_jacobian(self):
                return True

            @property
            def backend_name(self):
                return "dummy"

            def is_available_for_sensor(self, sensor_id, satellite_id):
                del sensor_id, satellite_id
                return True

        mask = xr.DataArray(
            np.array([[True, False], [False, False]], dtype=bool),
            dims=["y", "x"],
        )
        toa = xr.DataArray(np.full((3, 2, 2), 0.2, dtype=np.float32), dims=["band", "y", "x"])
        surface_prior = SurfacePrior(
            boa=xr.DataArray(np.full((2, 2), 0.12, dtype=np.float32), dims=["y", "x"]),
            boa_unc=xr.DataArray(np.full((2, 2), 0.02, dtype=np.float32), dims=["y", "x"]),
            kernels=cost_function_inputs["surface_prior"].kernels,
            mask=xr.DataArray(np.ones((2, 2), dtype=bool), dims=["y", "x"]),
        )
        atmo_prior = AtmosphericState(
            aot=xr.DataArray(np.full((2, 2), 0.15, dtype=np.float32), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full((2, 2), 2.5, dtype=np.float32), dims=["y", "x"]),
            tco3=xr.DataArray(np.full((2, 2), 0.3, dtype=np.float32), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full((2, 2), 0.05, dtype=np.float32), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full((2, 2), 0.3, dtype=np.float32), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full((2, 2), 0.01, dtype=np.float32), dims=["y", "x"]),
            elevation=xr.DataArray(np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"]),
        )
        cf = CostFunction(
            toa,
            surface_prior,
            cost_function_inputs["geometry"],
            atmo_prior,
            _DummyRT(),
            cost_function_inputs["bands"],
            mask,
        )

        aot = np.full((2, 2), 0.25, dtype=np.float32)
        tcwv = np.full((2, 2), 3.0, dtype=np.float32)
        j_prior, grad = cf.prior_cost_only(aot, tcwv)
        n = aot.size
        grad_aot = grad[:n].reshape(aot.shape)
        grad_tcwv = grad[n:].reshape(tcwv.shape)

        expected = 0.5 * ((0.10**2) / (0.05**2) + (0.50**2) / (0.3**2))
        assert j_prior == pytest.approx(expected)
        np.testing.assert_allclose(grad_aot, np.array([[40.0, 0.0], [0.0, 0.0]], dtype=np.float32))
        np.testing.assert_allclose(
            grad_tcwv,
            np.array([[5.5555553, 0.0], [0.0, 0.0]], dtype=np.float32),
            rtol=1e-6,
            atol=1e-6,
        )

    def test_smoothness_cost_ignores_edges_touching_masked_pixels(self, cost_function_inputs):
        class _DummyRT:
            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):
                del geometry, atmo_state, band, compute_jacobian
                raise AssertionError("RT coefficients are not used in smoothness_cost_only")

            def supports_jacobian(self):
                return True

            @property
            def backend_name(self):
                return "dummy"

            def is_available_for_sensor(self, sensor_id, satellite_id):
                del sensor_id, satellite_id
                return True

        mask = xr.DataArray(
            np.array([[True, False], [False, False]], dtype=bool),
            dims=["y", "x"],
        )
        toa = xr.DataArray(np.full((3, 2, 2), 0.2, dtype=np.float32), dims=["band", "y", "x"])
        surface_prior = SurfacePrior(
            boa=xr.DataArray(np.full((2, 2), 0.12, dtype=np.float32), dims=["y", "x"]),
            boa_unc=xr.DataArray(np.full((2, 2), 0.02, dtype=np.float32), dims=["y", "x"]),
            kernels=cost_function_inputs["surface_prior"].kernels,
            mask=xr.DataArray(np.ones((2, 2), dtype=bool), dims=["y", "x"]),
        )
        atmo_prior = AtmosphericState(
            aot=xr.DataArray(np.full((2, 2), 0.15, dtype=np.float32), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full((2, 2), 2.5, dtype=np.float32), dims=["y", "x"]),
            tco3=xr.DataArray(np.full((2, 2), 0.3, dtype=np.float32), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full((2, 2), 0.05, dtype=np.float32), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full((2, 2), 0.3, dtype=np.float32), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full((2, 2), 0.01, dtype=np.float32), dims=["y", "x"]),
            elevation=xr.DataArray(np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"]),
        )
        cf = CostFunction(
            toa,
            surface_prior,
            cost_function_inputs["geometry"],
            atmo_prior,
            _DummyRT(),
            cost_function_inputs["bands"],
            mask,
        )

        aot = np.array([[0.15, 1.0], [2.0, 3.0]], dtype=np.float32)
        tcwv = np.array([[2.5, 1.0], [2.0, 3.0]], dtype=np.float32)
        j_smooth, grad = cf.smoothness_cost_only(aot, tcwv)

        assert j_smooth == pytest.approx(0.0)
        np.testing.assert_allclose(grad, 0.0, atol=1e-8)


# ── Error path / edge-case tests ────────────────────────────────────


def _make_geometry(shape: tuple[int, int]):  # noqa: ANN201
    from siac.runtime import GeometryAngles

    return GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.5, dtype=np.float32), dims=["y", "x"]),
        saa=xr.DataArray(np.full(shape, 1.0, dtype=np.float32), dims=["y", "x"]),
        vza=xr.DataArray(np.full(shape, 0.2, dtype=np.float32), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(shape, 2.0, dtype=np.float32), dims=["y", "x"]),
    )


class _StubRT:
    """Minimal stub that satisfies isinstance(rt, RTModelBackend)."""

    pass


# Register _StubRT as a virtual subclass of RTModelBackend so isinstance checks pass.
from siac.domain.protocols import RTModelBackend as _RTProtocol  # noqa: E402

_RTProtocol.register(_StubRT)


class TestSolverErrorPaths:
    """Tests for solver error conditions and edge cases."""

    def test_solver_all_cloudy_returns_prior(self):
        """When all pixels are cloudy, solver should return the prior unchanged."""
        shape = (4, 4)
        toa = xr.DataArray(
            np.full((2, *shape), 0.2, dtype=np.float32),
            dims=["band", "y", "x"],
        )
        cloud_mask = xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"])

        atmo = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.15, dtype=np.float32), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.0, dtype=np.float32), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
        )
        sp = SurfacePrior(
            boa=xr.DataArray(np.full((2, *shape), 0.1, dtype=np.float32), dims=["band", "y", "x"]),
            boa_unc=xr.DataArray(
                np.full((2, *shape), 0.02, dtype=np.float32), dims=["band", "y", "x"]
            ),
            kernels=None,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )
        from siac.domain import SensorBand

        bands = [SensorBand("B1", 490.0, 10.0, 10.0, 0), SensorBand("B2", 560.0, 10.0, 10.0, 1)]
        solver = MultiGridSolver(MultiGridConfig())
        result = solver.solve(toa, sp, _make_geometry(shape), atmo, _StubRT(), cloud_mask, bands)
        assert not result.success
        assert result.n_iterations == 0

    def test_solver_rt_type_validation(self):
        """Solver should reject invalid RT model type."""
        shape = (2, 2)
        toa = xr.DataArray(np.full((1, *shape), 0.2, dtype=np.float32), dims=["band", "y", "x"])
        cloud_mask = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"])
        atmo = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.15, dtype=np.float32), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.0, dtype=np.float32), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
        )
        sp = SurfacePrior(
            boa=xr.DataArray(np.full((1, *shape), 0.1, dtype=np.float32), dims=["band", "y", "x"]),
            boa_unc=xr.DataArray(
                np.full((1, *shape), 0.02, dtype=np.float32), dims=["band", "y", "x"]
            ),
            kernels=None,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )
        solver = MultiGridSolver()
        with pytest.raises(TypeError, match="RTModelBackend"):
            solver.solve(toa, sp, _make_geometry(shape), atmo, "not_a_model", cloud_mask, [])

    def test_build_solver_valid_mask_all_nan_toa(self):
        """build_solver_valid_mask with all-NaN TOA should produce all-False mask."""
        shape = (3, 3)
        cloud = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"])
        toa = xr.DataArray(np.full(shape, np.nan, dtype=np.float32), dims=["y", "x"])
        sp = SurfacePrior(
            boa=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
            boa_unc=xr.DataArray(np.full(shape, 0.02, dtype=np.float32), dims=["y", "x"]),
            kernels=None,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )
        result = build_solver_valid_mask(cloud, toa, sp)
        assert not result.values.any()

    def test_build_solver_valid_mask_with_water_exclusion(self):
        """Water mask should exclude pixels from valid mask."""
        shape = (4, 4)
        cloud = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"])
        toa = xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"])
        sp = SurfacePrior(
            boa=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
            boa_unc=xr.DataArray(np.full(shape, 0.02, dtype=np.float32), dims=["y", "x"]),
            kernels=None,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )
        water = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"])
        water.values[0, 0] = True  # mark one pixel as water

        result = build_solver_valid_mask(cloud, toa, sp, water_mask=water)
        assert not result.values[0, 0]  # water pixel excluded
        assert result.values[1, 1]  # non-water pixel included

    def test_build_solver_valid_mask_out_of_range_toa(self):
        """TOA values outside (0, 1) should be excluded."""
        shape = (3, 3)
        cloud = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"])
        toa = xr.DataArray(np.full(shape, 0.5, dtype=np.float32), dims=["y", "x"])
        toa.values[0, 0] = -0.1  # below range
        toa.values[0, 1] = 1.5  # above range
        sp = SurfacePrior(
            boa=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
            boa_unc=xr.DataArray(np.full(shape, 0.02, dtype=np.float32), dims=["y", "x"]),
            kernels=None,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )
        result = build_solver_valid_mask(cloud, toa, sp)
        assert not result.values[0, 0]  # negative TOA excluded
        assert not result.values[0, 1]  # TOA > 1 excluded
        assert result.values[1, 1]  # valid pixel included

    def test_solver_config_defaults(self):
        """MultiGridConfig default values should be sane."""
        config = MultiGridConfig()
        assert config.n_levels >= 1
        assert config.min_grid_size >= 2
        assert config.aot_bounds[0] < config.aot_bounds[1]
        assert config.tcwv_bounds[0] < config.tcwv_bounds[1]
        assert config.max_iter_per_level >= 1
