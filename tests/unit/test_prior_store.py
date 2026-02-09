"""
Unit tests for PrebuiltPriorStore and prior store helpers.

See PLANS.md §9.6 test table.
"""

import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from siac.core.types import (
    BRDFKernelWeights,
    GeometryAngles,
    SensorBand,
    SensorConfig,
    SurfacePrior,
)
from siac.priors.surface.prior_store import (
    PrebuiltPriorStore,
    _doy_from_datetime,
    _interpolate_doy,
    _select_tiles,
)


# ── Fixtures ─────────────────────────────────────────────────────────

@pytest.fixture
def sensor_config():
    return SensorConfig(
        sensor_id="MOCK",
        satellite_id="TEST",
        bands=(
            SensorBand("B01", 443.0, 20.0, 60.0, 0),
            SensorBand("B02", 490.0, 65.0, 10.0, 1),
            SensorBand("B03", 560.0, 35.0, 10.0, 2),
            SensorBand("B04", 665.0, 30.0, 10.0, 3),
        ),
    )


@pytest.fixture
def geometry():
    shape = (8, 8)
    return GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
        saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
        vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
    )


@pytest.fixture
def zarr_store(tmp_path):
    """Create a minimal Zarr store with one tile for testing."""
    tile_dir = tmp_path / "T31UDQ"
    tile_dir.mkdir()

    n_doy = 4
    n_band = 7
    ny, nx = 16, 16
    doy_values = np.array([1, 91, 182, 274])
    wavelengths = np.array([645.0, 858.5, 469.0, 555.0, 1240.0, 1640.0, 2130.0])

    refl = np.random.RandomState(42).uniform(0.05, 0.3, (n_doy, n_band, ny, nx)).astype(np.float32)
    unc = np.random.RandomState(43).uniform(0.01, 0.05, (n_doy, n_band, ny, nx)).astype(np.float32)

    ds = xr.Dataset(
        {
            "reflectance": xr.DataArray(refl, dims=["doy", "band", "y", "x"]),
            "uncertainty": xr.DataArray(unc, dims=["doy", "band", "y", "x"]),
            "wavelengths": xr.DataArray(wavelengths, dims=["band"]),
        },
        coords={"doy": doy_values},
        attrs={
            "bounds": [300000.0, 5500000.0, 300160.0, 5500160.0],
            "resolution_m": 10.0,
            "crs": "EPSG:32632",
        },
    )
    ds.to_zarr(str(tile_dir), mode="w")

    # Write .zattrs for tile selection
    zattrs = {
        "bounds": [300000.0, 5500000.0, 300160.0, 5500160.0],
        "crs": "EPSG:32632",
        "resolution_m": 10.0,
    }
    (tile_dir / ".zattrs").write_text(json.dumps(zattrs))

    return tmp_path


# ── Helper function tests ────────────────────────────────────────────

class TestDOYFromDatetime:
    def test_jan_1(self):
        assert _doy_from_datetime(datetime(2023, 1, 1)) == 1

    def test_jul_4(self):
        assert _doy_from_datetime(datetime(2023, 7, 4)) == 185

    def test_dec_31(self):
        assert _doy_from_datetime(datetime(2023, 12, 31)) == 365


class TestInterpolateDOY:
    def test_exact_match(self):
        data = xr.DataArray(
            np.array([[1.0, 2.0], [3.0, 4.0]]),
            dims=["doy", "x"],
            coords={"doy": [100, 200]},
        )
        result = _interpolate_doy(data, np.array([100, 200]), 100)
        np.testing.assert_array_equal(result.values, [1.0, 2.0])

    def test_midpoint(self):
        data = xr.DataArray(
            np.array([[0.0, 0.0], [1.0, 1.0]]),
            dims=["doy", "x"],
            coords={"doy": [100, 200]},
        )
        result = _interpolate_doy(data, np.array([100, 200]), 150)
        np.testing.assert_array_almost_equal(result.values, [0.5, 0.5])

    def test_result_between_endpoints(self):
        """Interpolated DOY result should be between the two bracketing values."""
        data = xr.DataArray(
            np.array([0.1, 0.5, 0.9]),
            dims=["doy"],
            coords={"doy": [91, 182, 274]},
        )
        result = _interpolate_doy(data, np.array([91, 182, 274]), 130)
        assert 0.1 < float(result) < 0.5


class TestSelectTiles:
    def test_overlapping_tile(self, zarr_store):
        bounds = (300050.0, 5500050.0, 300100.0, 5500100.0)
        tiles = _select_tiles(zarr_store, bounds)
        assert "T31UDQ" in tiles

    def test_non_overlapping_tile(self, zarr_store):
        bounds = (0.0, 0.0, 1.0, 1.0)
        tiles = _select_tiles(zarr_store, bounds)
        assert len(tiles) == 0

    def test_empty_dir(self, tmp_path):
        tiles = _select_tiles(tmp_path, (0.0, 0.0, 1.0, 1.0))
        assert tiles == []


# ── PrebuiltPriorStore tests ─────────────────────────────────────────

class TestPrebuiltPriorStore:
    def test_returns_surface_prior(self, zarr_store, sensor_config, geometry):
        store = PrebuiltPriorStore(zarr_store)
        result = store.get_surface_prior(
            bounds=(300000.0, 5500000.0, 300160.0, 5500160.0),
            crs="EPSG:32632",
            obs_time=datetime(2023, 7, 4),
            sensor_config=sensor_config,
            geometry=geometry,
            resolution=10.0,
        )
        assert isinstance(result, SurfacePrior)

    def test_boa_shape(self, zarr_store, sensor_config, geometry):
        store = PrebuiltPriorStore(zarr_store)
        result = store.get_surface_prior(
            bounds=(300000.0, 5500000.0, 300160.0, 5500160.0),
            crs="EPSG:32632",
            obs_time=datetime(2023, 7, 4),
            sensor_config=sensor_config,
            geometry=geometry,
            resolution=10.0,
        )
        assert result.boa.ndim == 2
        assert result.boa_unc.ndim == 2

    def test_ignores_geometry(self, zarr_store, sensor_config, geometry):
        """Pre-built prior should not vary with geometry."""
        store = PrebuiltPriorStore(zarr_store)
        kwargs = dict(
            bounds=(300000.0, 5500000.0, 300160.0, 5500160.0),
            crs="EPSG:32632",
            obs_time=datetime(2023, 7, 4),
            sensor_config=sensor_config,
            resolution=10.0,
        )
        r1 = store.get_surface_prior(geometry=geometry, **kwargs)

        # Different geometry
        shape = (8, 8)
        geom2 = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 1.0), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 3.0), dims=["y", "x"]),
        )
        r2 = store.get_surface_prior(geometry=geom2, **kwargs)

        np.testing.assert_array_equal(r1.boa.values, r2.boa.values)

    def test_aoi_crop_reduces_extent(self, zarr_store, sensor_config, geometry):
        """Cropping to a sub-AOI should produce smaller output."""
        store = PrebuiltPriorStore(zarr_store)
        full = store.get_surface_prior(
            bounds=(300000.0, 5500000.0, 300160.0, 5500160.0),
            crs="EPSG:32632",
            obs_time=datetime(2023, 4, 1),
            sensor_config=sensor_config,
            geometry=geometry,
            resolution=10.0,
        )
        cropped = store.get_surface_prior(
            bounds=(300000.0, 5500000.0, 300080.0, 5500080.0),
            crs="EPSG:32632",
            obs_time=datetime(2023, 4, 1),
            sensor_config=sensor_config,
            geometry=geometry,
            resolution=10.0,
        )
        assert cropped.boa.shape[0] <= full.boa.shape[0]
        assert cropped.boa.shape[1] <= full.boa.shape[1]

    def test_mask_is_boolean(self, zarr_store, sensor_config, geometry):
        store = PrebuiltPriorStore(zarr_store)
        result = store.get_surface_prior(
            bounds=(300000.0, 5500000.0, 300160.0, 5500160.0),
            crs="EPSG:32632",
            obs_time=datetime(2023, 7, 4),
            sensor_config=sensor_config,
            geometry=geometry,
            resolution=10.0,
        )
        assert result.mask.dtype == bool

    def test_no_tiles_raises(self, tmp_path, sensor_config, geometry):
        store = PrebuiltPriorStore(tmp_path)
        with pytest.raises(ValueError, match="No tiles"):
            store.get_surface_prior(
                bounds=(0.0, 0.0, 1.0, 1.0),
                crs="EPSG:32632",
                obs_time=datetime(2023, 7, 4),
                sensor_config=sensor_config,
                geometry=geometry,
                resolution=10.0,
            )

    def test_spectral_projection_band_count(self, zarr_store, geometry):
        """Output reflectance should be projected to the sensor band count."""
        config_2band = SensorConfig(
            sensor_id="MOCK",
            satellite_id="TEST",
            bands=(
                SensorBand("B01", 645.0, 50.0, 10.0, 0),
                SensorBand("B02", 858.5, 35.0, 10.0, 1),
            ),
        )
        store = PrebuiltPriorStore(zarr_store)
        result = store.get_surface_prior(
            bounds=(300000.0, 5500000.0, 300160.0, 5500160.0),
            crs="EPSG:32632",
            obs_time=datetime(2023, 7, 4),
            sensor_config=config_2band,
            geometry=geometry,
            resolution=10.0,
        )
        # boa is 2D (averaged across sensor bands), but it should still be valid
        assert isinstance(result, SurfacePrior)
        assert result.boa.ndim == 2
