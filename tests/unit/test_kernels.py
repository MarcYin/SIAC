"""Unit tests for BRDF kernel calculations."""

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.brdf.kernels import (
    BRDFKernels,
    compute_black_sky_albedo,
    compute_kernels,
    compute_reflectance,
    compute_white_sky_albedo,
)


class TestBRDFKernels:
    """Tests for BRDFKernels class."""

    @pytest.fixture
    def kernels(self):
        """Create kernel calculator with MODIS parameters."""
        return BRDFKernels(hb=2.0, br=1.0)

    def test_creation(self, kernels):
        """Kernels should be creatable."""
        assert kernels.hb == 2.0
        assert kernels.br == 1.0

    def test_compute_nadir(self, kernels):
        """Nadir geometry should give known kernel values."""
        # Nadir: vza=0, sza=30deg, raa=0
        vza = np.array([[0.0]])
        sza = np.array([[np.deg2rad(30.0)]])
        raa = np.array([[0.0]])

        k_vol, k_geo = kernels.compute(vza, sza, raa)

        # Ross-Thick at nadir with sza=30 should be negative
        assert k_vol[0, 0] < 0
        # Li-Sparse at nadir should be around -1 to -2
        assert -3 < k_geo[0, 0] < 0

    def test_compute_hotspot(self, kernels):
        """Hotspot geometry (vza=sza, raa=0) should have specific behavior."""
        sza = np.deg2rad(30.0)
        vza = np.array([[sza]])
        sza_arr = np.array([[sza]])
        raa = np.array([[0.0]])  # Backscatter

        k_vol, k_geo = kernels.compute(vza, sza_arr, raa)

        # At hotspot, geometric kernel should be less negative
        assert k_geo[0, 0] > -2

    def test_compute_symmetry(self, kernels):
        """Kernels should be symmetric in RAA."""
        vza = np.array([[np.deg2rad(10.0)]])
        sza = np.array([[np.deg2rad(30.0)]])
        raa_pos = np.array([[np.deg2rad(45.0)]])
        raa_neg = np.array([[np.deg2rad(-45.0)]])

        k_vol_pos, k_geo_pos = kernels.compute(vza, sza, raa_pos)
        k_vol_neg, k_geo_neg = kernels.compute(vza, sza, raa_neg)

        np.testing.assert_allclose(k_vol_pos, k_vol_neg, rtol=1e-5)
        np.testing.assert_allclose(k_geo_pos, k_geo_neg, rtol=1e-5)

    def test_compute_xarray_input(self, kernels):
        """Should handle xarray inputs."""
        shape = (10, 10)
        vza = xr.DataArray(np.full(shape, 0.1), dims=["y", "x"])
        sza = xr.DataArray(np.full(shape, 0.5), dims=["y", "x"])
        raa = xr.DataArray(np.full(shape, 1.0), dims=["y", "x"])

        k_vol, k_geo = kernels.compute(vza, sza, raa)

        assert isinstance(k_vol, xr.DataArray)
        assert k_vol.shape == shape

    def test_compute_batch(self, kernels):
        """Should handle batched inputs efficiently."""
        shape = (100, 100)
        vza = np.random.default_rng().uniform(0, 0.5, shape)
        sza = np.random.default_rng().uniform(0.2, 1.0, shape)
        raa = np.random.default_rng().uniform(0, np.pi, shape)

        k_vol, k_geo = kernels.compute(vza, sza, raa)

        assert k_vol.shape == shape
        assert np.all(np.isfinite(k_vol))


class TestComputeReflectance:
    """Tests for compute_reflectance function."""

    def test_basic(self):
        """Basic reflectance computation."""
        f0 = np.array([[0.1]])
        f1 = np.array([[0.05]])
        f2 = np.array([[0.02]])
        k_vol = np.array([[-0.3]])
        k_geo = np.array([[-1.0]])

        ref = compute_reflectance(f0, f1, f2, k_vol, k_geo)

        expected = 0.1 + 0.05 * (-0.3) + 0.02 * (-1.0)
        np.testing.assert_allclose(ref, expected)

    def test_xarray(self):
        """Should work with xarray inputs."""
        shape = (5, 5)
        f0 = xr.DataArray(np.full(shape, 0.1), dims=["y", "x"])
        f1 = xr.DataArray(np.full(shape, 0.05), dims=["y", "x"])
        f2 = xr.DataArray(np.full(shape, 0.02), dims=["y", "x"])
        k_vol = xr.DataArray(np.full(shape, -0.3), dims=["y", "x"])
        k_geo = xr.DataArray(np.full(shape, -1.0), dims=["y", "x"])

        ref = compute_reflectance(f0, f1, f2, k_vol, k_geo)

        assert isinstance(ref, xr.DataArray)


class TestAlbedoFunctions:
    """Tests for albedo computation functions."""

    def test_white_sky_albedo(self):
        """White sky albedo should use integration weights."""
        f0 = np.array([[0.1]])
        f1 = np.array([[0.05]])
        f2 = np.array([[0.02]])

        wsa = compute_white_sky_albedo(f0, f1, f2)

        # WSA = f0 * 1.0 + f1 * 0.189184 + f2 * (-1.377622)
        expected = 0.1 + 0.05 * 0.189184 + 0.02 * (-1.377622)
        np.testing.assert_allclose(wsa, expected, rtol=1e-5)

    def test_black_sky_albedo_nadir(self):
        """Black sky albedo at nadir (sza=0)."""
        f0 = np.array([[0.1]])
        f1 = np.array([[0.05]])
        f2 = np.array([[0.02]])
        sza = np.array([[0.0]])

        bsa = compute_black_sky_albedo(f0, f1, f2, sza)

        # At sza=0, polynomial coefficients simplify
        assert np.isfinite(bsa[0, 0])
        # BSA should be close to WSA at nadir
        wsa = compute_white_sky_albedo(f0, f1, f2)
        assert abs(bsa[0, 0] - wsa[0, 0]) < 0.1


class TestKernelCoverageExtras:
    def test_compute_kernels_modis_and_invalid(self):
        vza = np.array([[0.1]], dtype=np.float32)
        sza = np.array([[0.2]], dtype=np.float32)
        raa = np.array([[0.3]], dtype=np.float32)

        k_vol, k_geo = compute_kernels(vza, sza, raa, kernel_type="modis")
        assert k_vol.shape == (1, 1)
        assert k_geo.shape == (1, 1)

        with pytest.raises(ValueError, match="Unknown kernel type"):
            compute_kernels(vza, sza, raa, kernel_type="roujean")

    def test_rust_kernels_preserve_non_2d_input_shape(self):
        kernels = BRDFKernels()
        vza = np.array([0.0, 0.1, 0.2], dtype=np.float32)
        sza = np.array([0.3, 0.4, 0.5], dtype=np.float32)
        raa = np.array([0.1, 0.2, 0.3], dtype=np.float32)

        k_vol, k_geo = kernels.compute(vza, sza, raa)
        assert k_vol.shape == (3,)
        assert k_geo.shape == (3,)
        assert np.all(np.isfinite(k_vol))
        assert np.all(np.isfinite(k_geo))
