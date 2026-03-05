"""
Unit tests for solver components.
"""

import numpy as np
import pytest
import xarray as xr

from siac.core.types import AtmosphericState
from siac.solver.cost import (
    CostFunction,
    CostFunctionConfig,
    apply_smoothness_filter,
    compute_laplacian_eigenvalues,
)
from siac.solver.multigrid import MultiGridConfig, MultiGridSolver


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

        from siac.core.types import BRDFKernelWeights, GeometryAngles, SensorBand, SurfacePrior

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

        from siac.core.types import GeometryAngles
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
