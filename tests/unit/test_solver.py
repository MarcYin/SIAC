"""
Unit tests for solver components.
"""

import numpy as np
import pytest
import xarray as xr

from siac._rust import (
    evaluate_grid_search_candidate_cost,
    evaluate_grid_search_cost_cube_with_provider,
    quadratic_refine_grid_search,
)
from siac.algorithms.solver.cost import (
    CostFunction,
    CostFunctionConfig,
    apply_smoothness_filter,
    compute_laplacian_eigenvalues,
)
from siac.algorithms.solver.multigrid import MultiGridConfig, MultiGridSolver
from siac.runtime import AtmosphericState


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
                np.testing.assert_allclose(
                    lambda_vals[i, j],
                    lambda_vals[j, i],
                    rtol=1e-10
                )


class TestSmoothnessFilter:
    """Tests for smoothness filter application."""

    def test_identity_gamma_zero(self):
        """With gamma=0, output should equal input."""
        x = np.random.default_rng(0).standard_normal((20, 20))
        lambda_vals = compute_laplacian_eigenvalues(20, 20)

        x_smooth = apply_smoothness_filter(x, gamma=0, lambda_vals=lambda_vals)

        np.testing.assert_allclose(x_smooth, x, rtol=1e-10)

    def test_smoothing_effect(self):
        """With gamma>0, high frequencies should be suppressed."""
        # Create noisy data
        x = np.zeros((32, 32))
        x[16, 16] = 1.0  # Single spike

        lambda_vals = compute_laplacian_eigenvalues(32, 32)

        x_smooth = apply_smoothness_filter(x, gamma=5.0, lambda_vals=lambda_vals)

        # Smoothed data should have lower variance
        assert x_smooth.var() < x.var()
        # Peak should be reduced
        assert x_smooth.max() < x.max()

    def test_preserves_mean(self):
        """Mean should be approximately preserved."""
        x = np.random.default_rng(1).standard_normal((20, 20)) + 5.0
        lambda_vals = compute_laplacian_eigenvalues(20, 20)

        x_smooth = apply_smoothness_filter(x, gamma=2.0, lambda_vals=lambda_vals)

        np.testing.assert_allclose(x.mean(), x_smooth.mean(), rtol=0.01)


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

    def test_bounds(self):
        """Bounds should be tuple of two floats."""
        config = MultiGridConfig(
            aot_bounds=(0.01, 1.5),
            tcwv_bounds=(0.5, 5.0),
        )

        assert config.aot_bounds[0] == 0.01
        assert config.tcwv_bounds[1] == 5.0


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
            boa=xr.DataArray(np.stack([np.full(shape, 0.2), np.full(shape, 0.22)]), dims=["band", "y", "x"]),
            boa_unc=xr.DataArray(np.stack([np.full(shape, 0.02), np.full(shape, 0.02)]), dims=["band", "y", "x"]),
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

    def test_grid_search_prefers_rust_refiner_when_available(self, monkeypatch):
        """Grid-search refinement should call Rust helper."""
        from siac.algorithms.solver import multigrid as mg_mod
        from siac.domain import SensorBand
        from siac.runtime import BRDFKernelWeights, GeometryAngles, RTCoefficients, SurfacePrior

        shape = (4, 4)
        solver = MultiGridSolver(MultiGridConfig(grid_search_aot_points=3, grid_search_tcwv_points=3))

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

        def _fake_rust_refiner(costs, aot_axis, tcwv_axis, valid_mask):  # noqa: ANN001
            called["n"] += 1
            out_aot = np.full(valid_mask.shape, 0.33, dtype=np.float32)
            out_tcwv = np.full(valid_mask.shape, 1.23, dtype=np.float32)
            out_aot_unc = np.full(valid_mask.shape, 0.07, dtype=np.float32)
            out_tcwv_unc = np.full(valid_mask.shape, 0.21, dtype=np.float32)
            assert costs.shape[:2] == (3, 3)
            assert aot_axis.shape[0] == 3
            assert tcwv_axis.shape[0] == 3
            return out_aot, out_tcwv, out_aot_unc, out_tcwv_unc

        monkeypatch.setattr(mg_mod, "quadratic_refine_grid_search", _fake_rust_refiner)

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
        valid_mask = (rng.random((ny, nx)) > 0.2)

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
        valid_mask = (rng.random((ny, nx)) > 0.15)
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
        solver = MultiGridSolver(MultiGridConfig(grid_search_aot_points=3, grid_search_tcwv_points=3))

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
        """Smoothness gradient should match numerical gradient."""
        geometry, atmo_prior, shape = simple_cost_inputs

        config = CostFunctionConfig(aot_gamma=5.0, tcwv_gamma=3.0)

        # Just test the smoothness part directly
        from scipy.fftpack import dct, idct

        aot = np.random.default_rng(4).standard_normal(shape) * 0.1 + 0.15
        lambda_vals = compute_laplacian_eigenvalues(shape[1], shape[0])

        # Analytical gradient
        gamma = config.aot_gamma
        aot_dct = dct(dct(aot, axis=0, norm="ortho"), axis=1, norm="ortho")
        0.5 * gamma**2 * np.sum(lambda_vals * aot_dct**2)
        dj_dct = gamma**2 * lambda_vals * aot_dct
        dj_aot = idct(idct(dj_dct, axis=1, norm="ortho"), axis=0, norm="ortho")

        # Numerical gradient
        eps = 1e-5
        numerical_grad = np.zeros_like(aot)

        for i in range(3):  # Test just a few points
            for j in range(3):
                aot_plus = aot.copy()
                aot_plus[i, j] += eps
                aot_dct_plus = dct(dct(aot_plus, axis=0, norm="ortho"), axis=1, norm="ortho")
                j_plus = 0.5 * gamma**2 * np.sum(lambda_vals * aot_dct_plus**2)

                aot_minus = aot.copy()
                aot_minus[i, j] -= eps
                aot_dct_minus = dct(dct(aot_minus, axis=0, norm="ortho"), axis=1, norm="ortho")
                j_minus = 0.5 * gamma**2 * np.sum(lambda_vals * aot_dct_minus**2)

                numerical_grad[i, j] = (j_plus - j_minus) / (2 * eps)

        # Compare at tested points
        np.testing.assert_allclose(
            dj_aot[:3, :3], numerical_grad[:3, :3], rtol=1e-4, atol=1e-6
        )
