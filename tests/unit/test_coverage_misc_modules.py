"""Coverage-oriented unit tests for utility and provider modules."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
import rasterio
import xarray as xr
from rasterio.transform import from_origin

from siac.core.aoi import AOI, _detect_crs
from siac.core.exceptions import (
    ConfigurationError,
    DataNotFoundError,
    RTModelError,
    SIACError,
    SensorNotSupportedError,
    SolverConvergenceError,
    ValidationError,
)
from siac.core.types import BRDFKernelWeights, GeometryAngles
from siac.io.copernicus_dataspace import CopernicusDataspaceBackend, download_cdse, search_cdse
from siac.io.gcs_sentinel2 import GCSSentinel2Backend, download_gcs, search_gcs
from siac.io.s2_data_source import S2Product, S2Query
from siac.priors.atmospheric.cams import CAMSProvider
from siac.priors.surface.kernel_model import KernelModelDeriver, PSFConvolver
from siac.rt.direct import __all__ as direct_all


def _spatial_da(values: np.ndarray, crs: str | None = "EPSG:32632") -> xr.DataArray:
    """Create a spatial DataArray with x/y dims and optional CRS."""
    h, w = values.shape[-2], values.shape[-1]
    x = np.linspace(500000.0, 500000.0 + (w - 1) * 10.0, w)
    y = np.linspace(4500000.0 + (h - 1) * 10.0, 4500000.0, h)
    da = xr.DataArray(values, dims=["y", "x"], coords={"y": y, "x": x})
    if crs is not None:
        da = da.rio.write_crs(crs)
    return da


def _geometry(shape: tuple[int, int] = (8, 8)) -> GeometryAngles:
    """Small deterministic geometry fixture."""
    sza = _spatial_da(np.full(shape, np.deg2rad(30.0), dtype=np.float32))
    saa = _spatial_da(np.full(shape, np.deg2rad(150.0), dtype=np.float32))
    vza = _spatial_da(np.full(shape, np.deg2rad(5.0), dtype=np.float32))
    vaa = _spatial_da(np.full(shape, np.deg2rad(100.0), dtype=np.float32))
    return GeometryAngles(sza=sza, saa=saa, vza=vza, vaa=vaa)


def _brdf_weights(shape: tuple[int, int] = (8, 8)) -> BRDFKernelWeights:
    """Small deterministic BRDF weights."""
    return BRDFKernelWeights(
        f0=_spatial_da(np.full(shape, 0.10, dtype=np.float32)),
        f1=_spatial_da(np.full(shape, 0.05, dtype=np.float32)),
        f2=_spatial_da(np.full(shape, 0.02, dtype=np.float32)),
        f0_unc=_spatial_da(np.full(shape, 0.01, dtype=np.float32)),
        f1_unc=_spatial_da(np.full(shape, 0.005, dtype=np.float32)),
        f2_unc=_spatial_da(np.full(shape, 0.002, dtype=np.float32)),
    )


class TestAOI:
    def test_from_bounds_and_to_geojson(self):
        aoi = AOI.from_bounds((0.0, 1.0, 2.0, 3.0), crs="EPSG:4326")
        assert aoi.crs == "EPSG:4326"
        assert aoi.to_geojson()["type"] == "Polygon"
        assert aoi.get_bounds() == (0.0, 1.0, 2.0, 3.0)

    def test_from_geojson_and_detect_crs(self, tmp_path: Path):
        geojson = {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [[[0, 0], [2, 0], [2, 1], [0, 1], [0, 0]]],
            },
            "crs": {"type": "name", "properties": {"name": "EPSG:4326"}},
        }
        p = tmp_path / "aoi.json"
        p.write_text(json.dumps(geojson))

        assert _detect_crs(geojson) == "EPSG:4326"
        assert _detect_crs(p) == "EPSG:4326"
        assert _detect_crs("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))") is None

        aoi = AOI.from_geojson(p)
        assert aoi.crs == "EPSG:4326"
        assert aoi.get_bounds() == (0.0, 0.0, 2.0, 1.0)

    def test_from_raster_and_crs_transform(self):
        da = _spatial_da(np.ones((4, 4), dtype=np.float32), crs="EPSG:32632")
        aoi = AOI.from_raster(da)
        assert aoi.crs == "EPSG:32632"
        xmin, ymin, xmax, ymax = aoi.get_bounds(target_crs="EPSG:4326")
        assert xmin < xmax
        assert ymin < ymax

    def test_from_raster_missing_crs_raises(self):
        da = _spatial_da(np.ones((2, 2), dtype=np.float32), crs=None)
        with pytest.raises(ValueError, match="CRS is not set"):
            AOI.from_raster(da)


class TestExceptionsAndStubs:
    def test_exception_hierarchy(self):
        for exc in [
            ConfigurationError,
            DataNotFoundError,
            SensorNotSupportedError,
            RTModelError,
            SolverConvergenceError,
            ValidationError,
        ]:
            assert issubclass(exc, SIACError)

    def test_rt_direct_module(self):
        assert direct_all == []

    def test_cdse_and_gcs_stubs(self, tmp_path: Path):
        query = S2Query(mgrs_tile="31UDQ")
        prod = S2Product(
            product_id="P",
            mgrs_tile="31UDQ",
            sensing_date=datetime(2024, 1, 1),
            processing_baseline="N0500",
            cloud_cover=10.0,
            satellite="S2A",
            orbit_number=1,
            source_url="s3://dummy",
        )

        with pytest.raises(NotImplementedError):
            search_cdse(query)
        with pytest.raises(NotImplementedError):
            download_cdse(prod, tmp_path)
        with pytest.raises(NotImplementedError):
            CopernicusDataspaceBackend().search(query)
        with pytest.raises(NotImplementedError):
            CopernicusDataspaceBackend().download(prod, tmp_path)

        with pytest.raises(NotImplementedError):
            search_gcs(query)
        with pytest.raises(NotImplementedError):
            download_gcs(prod, tmp_path)
        with pytest.raises(NotImplementedError):
            GCSSentinel2Backend().search(query)
        with pytest.raises(NotImplementedError):
            GCSSentinel2Backend().download(prod, tmp_path)


class TestCAMSProvider:
    def test_default_prior_when_missing(self, tmp_path: Path):
        p = CAMSProvider(tmp_path)
        state = p.get_prior((0.0, 0.0, 10.0, 10.0), "EPSG:4326", datetime(2024, 1, 1), 1.0)
        assert state.aot.shape == (10, 10)
        assert float(state.aot.mean()) == pytest.approx(0.15)
        assert float(state.tcwv.mean()) == pytest.approx(1.5)
        assert float(state.tco3.mean()) == pytest.approx(0.3)

    def test_extract_variable_missing_var_fallback(self, tmp_path: Path):
        p = CAMSProvider(tmp_path)
        ds = xr.Dataset({"dummy": (("y", "x"), np.ones((3, 3), dtype=np.float32))})
        arr = p._extract_variable(ds, "aod550", (0.0, 0.0, 4.0, 4.0), "EPSG:4326", 1.0, datetime(2024, 1, 1))
        assert arr.shape == (4, 4)
        assert float(arr.mean()) == pytest.approx(0.15)

    def test_get_prior_from_dataset_and_no_temporal_interp(self, tmp_path: Path):
        time = np.array([np.datetime64("2024-01-01T00:00:00"), np.datetime64("2024-01-01T03:00:00")])
        lat = np.array([1.0, 0.0, -1.0], dtype=np.float32)
        lon = np.array([0.0, 1.0], dtype=np.float32)
        shape = (2, 3, 2)
        ds = xr.Dataset(
            {
                "aod550": (("time", "latitude", "longitude"), np.full(shape, 0.2, dtype=np.float32)),
                "tcwv": (("time", "latitude", "longitude"), np.full(shape, 2.2, dtype=np.float32)),
                "gtco3": (("time", "latitude", "longitude"), np.full(shape, 0.006, dtype=np.float32)),
            },
            coords={"time": time, "latitude": lat, "longitude": lon},
        )

        p = CAMSProvider(tmp_path, temporal_interp=False)
        # Exercise _load_cams_data pattern search path
        f = tmp_path / "cams_20240101_test.nc"
        f.write_text("dummy")

        # Monkeypatch loader internals using closure
        orig = xr.open_mfdataset
        try:
            xr.open_mfdataset = lambda *args, **kwargs: ds  # type: ignore[assignment]
            state = p.get_prior((0.0, -1.0, 1.0, 1.0), "EPSG:4326", datetime(2024, 1, 1, 1), 1.0)
        finally:
            xr.open_mfdataset = orig  # type: ignore[assignment]

        assert state.aot.size > 0
        assert float(state.aot.mean()) == pytest.approx(0.2)
        # Ozone conversion path executes
        assert float(state.tco3.mean()) > 0.0
        assert float(state.aot_unc.mean()) >= 0.05

    def test_load_cams_handles_open_error(self, tmp_path: Path):
        p = CAMSProvider(tmp_path)
        (tmp_path / "cams_20240101_test.nc").write_text("dummy")

        orig = xr.open_mfdataset
        try:
            def _boom(*args, **kwargs):
                raise RuntimeError("read failed")
            xr.open_mfdataset = _boom  # type: ignore[assignment]
            assert p._load_cams_data(datetime(2024, 1, 1)) is None
        finally:
            xr.open_mfdataset = orig  # type: ignore[assignment]

    def test_load_cams_from_direct_netcdf_file(self, tmp_path: Path):
        ds = xr.Dataset(
            {
                "aod550": (("time", "latitude", "longitude"), np.full((1, 2, 2), 0.2, dtype=np.float32)),
                "tcwv": (("time", "latitude", "longitude"), np.full((1, 2, 2), 2.0, dtype=np.float32)),
                "gtco3": (("time", "latitude", "longitude"), np.full((1, 2, 2), 0.006, dtype=np.float32)),
            },
            coords={
                "time": [np.datetime64("2024-01-01T00:00:00")],
                "latitude": [1.0, 0.0],
                "longitude": [0.0, 1.0],
            },
        )
        nc_path = tmp_path / "cams_direct.nc"
        ds.to_netcdf(nc_path)

        p = CAMSProvider(nc_path)
        loaded = p._load_cams_data(datetime(2024, 1, 1))
        assert loaded is not None
        assert {"aod550", "tcwv", "gtco3"}.issubset(set(loaded.data_vars))

    def test_load_cams_from_direct_tif_file(self, tmp_path: Path):
        tif_path = tmp_path / "cams_stack.tif"
        transform = from_origin(0.0, 2.0, 1.0, 1.0)
        stack = np.stack(
            [
                np.full((2, 2), 0.2, dtype=np.float32),
                np.full((2, 2), 2.5, dtype=np.float32),
                np.full((2, 2), 0.006, dtype=np.float32),
            ],
            axis=0,
        )
        with rasterio.open(
            tif_path,
            "w",
            driver="GTiff",
            height=2,
            width=2,
            count=3,
            dtype="float32",
            crs="EPSG:4326",
            transform=transform,
        ) as dst:
            dst.write(stack)
            dst.set_band_description(1, "aod550")
            dst.set_band_description(2, "tcwv")
            dst.set_band_description(3, "gtco3")

        p = CAMSProvider(tif_path)
        loaded = p._load_cams_data(datetime(2024, 1, 1))
        assert loaded is not None
        assert {"aod550", "tcwv", "gtco3"}.issubset(set(loaded.data_vars))

        state = p.get_prior((0.0, 0.0, 2.0, 2.0), "EPSG:4326", datetime(2024, 1, 1), 1.0)
        assert float(state.aot.mean()) == pytest.approx(0.2)
        assert float(state.tcwv.mean()) == pytest.approx(2.5)
        expected_tco3 = 0.006 / 2.1415e-5 / 1000
        assert float(state.tco3.mean()) == pytest.approx(expected_tco3, rel=1e-6)

    def test_single_band_tif_infers_variable_from_filename(self, tmp_path: Path):
        tif_path = tmp_path / "cams_tcwv_20240101.tif"
        transform = from_origin(0.0, 2.0, 1.0, 1.0)
        arr = np.full((2, 2), 3.0, dtype=np.float32)
        with rasterio.open(
            tif_path,
            "w",
            driver="GTiff",
            height=2,
            width=2,
            count=1,
            dtype="float32",
            crs="EPSG:4326",
            transform=transform,
        ) as dst:
            dst.write(arr, 1)

        p = CAMSProvider(tif_path)
        state = p.get_prior((0.0, 0.0, 2.0, 2.0), "EPSG:4326", datetime(2024, 1, 1), 1.0)
        assert float(state.tcwv.mean()) == pytest.approx(3.0)
        assert float(state.aot.mean()) == pytest.approx(0.15)
        assert float(state.tco3.mean()) == pytest.approx(0.3)


class TestKernelModel:
    def test_compute_surface_prior_without_psf(self):
        geom = _geometry((6, 6))
        weights = _brdf_weights((6, 6))
        kmd = KernelModelDeriver(apply_psf=False)

        prior = kmd.compute_surface_prior(weights, geom)
        assert prior.boa.shape == (6, 6)
        assert prior.boa_unc.shape == (6, 6)
        assert prior.mask.dtype == bool

    def test_compute_surface_prior_with_psf_and_override(self):
        geom = _geometry((6, 6))
        # Add band dim to hit per-band PSF path
        f0 = xr.DataArray(
            np.full((2, 6, 6), 0.1, dtype=np.float32),
            dims=["band", "y", "x"],
            coords={"band": [1, 2]},
        ).rio.write_crs("EPSG:32632")
        f1 = xr.full_like(f0, 0.05)
        f2 = xr.full_like(f0, 0.02)
        weights = BRDFKernelWeights(
            f0=f0,
            f1=f1,
            f2=f2,
            f0_unc=xr.full_like(f0, 0.01),
            f1_unc=xr.full_like(f0, 0.005),
            f2_unc=xr.full_like(f0, 0.002),
        )
        # Geometry remains 2D and broadcasts over band
        kmd = KernelModelDeriver(apply_psf=True)

        prior = kmd.compute_surface_prior(weights, geom, psf_params=(1.0, 1.0))
        assert prior.boa.shape == (2, 6, 6)
        assert prior.boa_unc.shape == (2, 6, 6)

    def test_convolution_helpers_handle_nan(self):
        kmd = KernelModelDeriver(apply_psf=True)
        arr = np.array([[1.0, np.nan], [2.0, 3.0]], dtype=np.float32)
        out = kmd._convolve_2d(arr, sigma_x=1.0, sigma_y=1.0)
        assert out.shape == arr.shape
        assert np.isfinite(out[0, 0])

    def test_psf_convolver_numpy_and_xarray(self):
        c = PSFConvolver(sigma_x=0.8, sigma_y=1.2)
        arr = np.array([[1.0, 2.0], [3.0, np.nan]], dtype=np.float32)
        out_np = c.convolve(arr)
        assert out_np.shape == arr.shape

        da = _spatial_da(arr)
        out_da = c.convolve(da)
        assert isinstance(out_da, xr.DataArray)
        assert out_da.shape == da.shape
