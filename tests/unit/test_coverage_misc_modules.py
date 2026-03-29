"""Coverage-oriented unit tests for utility and provider modules."""

from __future__ import annotations

import json
import zipfile
from datetime import date, datetime
from io import BytesIO
from typing import TYPE_CHECKING

import numpy as np
import pytest
import rasterio
import requests
import xarray as xr
from rasterio.transform import from_origin

from siac.adapters.atmo.cams import CAMSProvider
from siac.adapters.data.copernicus_dataspace import (
    CopernicusDataspaceBackend,
    download_cdse,
    search_cdse,
)
from siac.adapters.data.gcs_sentinel2 import GCSSentinel2Backend
from siac.adapters.data.s2_data_source import S2Product, S2Query
from siac.algorithms.rt.direct import __all__ as direct_all
from siac.domain.aoi import AOI, _detect_crs
from siac.errors import (
    ConfigurationError,
    DataNotFoundError,
    RTModelError,
    SensorNotSupportedError,
    SIACError,
    SolverConvergenceError,
    ValidationError,
)
from siac.runtime import BRDFKernelWeights, GeometryAngles

if TYPE_CHECKING:
    from pathlib import Path


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


def _is_cdse_targeted_cov_run(request: pytest.FixtureRequest) -> bool:
    cov_sources = getattr(request.config.option, "cov_source", None) or []
    return "siac.adapters.data.copernicus_dataspace" in cov_sources


def _skip_native_heavy_for_cdse_cov(request: pytest.FixtureRequest) -> None:
    if _is_cdse_targeted_cov_run(request):
        pytest.skip("Skip native-heavy tests for CDSE-only coverage run.")


def _kernel_model_classes(request: pytest.FixtureRequest):
    """Import kernel-model classes lazily to avoid collection-time SciPy import."""
    _skip_native_heavy_for_cdse_cov(request)

    from siac.algorithms.surface.kernel_model import KernelModelDeriver, PSFConvolver

    return KernelModelDeriver, PSFConvolver


class TestAOI:
    def test_from_bounds_and_to_geojson(self):
        aoi = AOI.from_bounds((0.0, 1.0, 2.0, 3.0), crs="EPSG:4326")
        assert aoi.crs == "EPSG:4326"
        assert aoi.to_geojson()["type"] == "Polygon"
        assert aoi.get_bounds() == (0.0, 1.0, 2.0, 3.0)

    def test_from_bounds_defaults_to_validated_wgs84(self):
        aoi = AOI.from_bounds((10.0, 20.0, 12.0, 21.0))
        assert aoi.crs == "EPSG:4326"

        with pytest.raises(ValueError, match="longitude bounds must look like degrees"):
            AOI.from_bounds((500000.0, 4100000.0, 501000.0, 4101000.0))

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

        aoi_override = AOI.from_geojson(
            "POLYGON ((500000 4100000, 501000 4100000, 501000 4101000, 500000 4101000, 500000 4100000))",
            crs="EPSG:32632",
        )
        assert aoi_override.crs == "EPSG:32632"

    def test_from_raster_and_crs_transform(self, request: pytest.FixtureRequest):
        _skip_native_heavy_for_cdse_cov(request)
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

    def test_from_raster_missing_rio_accessor_raises(self):
        with pytest.raises(ValueError, match="missing rioxarray accessor"):
            AOI.from_raster(object())  # type: ignore[arg-type]

    def test_detect_crs_additional_paths(self, tmp_path: Path):
        bad_json = tmp_path / "bad.json"
        bad_json.write_text("{not-json")
        assert _detect_crs(bad_json) is None

        assert _detect_crs({"type": "Feature", "crs": "EPSG:4326"}) == "EPSG:4326"
        assert _detect_crs({"type": "Feature", "crs": {"properties": {"other": "x"}}}) is None


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

    def test_cdse_and_gcs_backends(self, tmp_path: Path, monkeypatch):
        query = S2Query(
            mgrs_tile="31UDQ",
            date=date(2024, 1, 1),
            max_cloud_cover=30.0,
        )
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

        class _Resp:
            def __init__(self, *, status_code=200, json_data=None, raw=b""):
                self.status_code = status_code
                self._json = json_data or {}
                self._raw = raw

            def raise_for_status(self):
                if self.status_code >= 400:
                    raise requests.HTTPError(response=self)

            def json(self):
                return self._json

            def iter_content(self, chunk_size=1024 * 1024):
                for i in range(0, len(self._raw), chunk_size):
                    yield self._raw[i : i + chunk_size]

        item = {
            "id": "S2A_MSIL1C_20240101T103101_N0500_R008_T31UDQ_20240101T120000",
            "properties": {
                "datetime": "2024-01-01T10:31:01Z",
                "eo:cloud_cover": 12.5,
                "s2:mgrs_tile": "31UDQ",
            },
            "assets": {
                "Product": {
                    "href": "https://download.example/products/1/$value",
                    "file:size": 1024,
                }
            },
        }
        page = {"features": [item], "links": []}

        zip_buffer = BytesIO()
        with zipfile.ZipFile(zip_buffer, mode="w") as zf:
            zf.writestr(f"{item['id']}.SAFE/manifest.safe", "ok")
        payload = zip_buffer.getvalue()

        def _fake_post(url, json=None, timeout=60, **kwargs):  # noqa: ARG001
            if "openid-connect/token" in url:
                return _Resp(json_data={"access_token": "token-123"})
            return _Resp(json_data=page)

        def _fake_get(url, headers=None, timeout=60, stream=False, **kwargs):  # noqa: ARG001
            if "download.example" in url:
                return _Resp(raw=payload)
            return _Resp(json_data=page)

        monkeypatch.setattr("siac.adapters.data.copernicus_dataspace.requests.post", _fake_post)
        monkeypatch.setattr("siac.adapters.data.copernicus_dataspace.requests.get", _fake_get)

        found = search_cdse(query)
        assert len(found) == 1
        safe = download_cdse(found[0], tmp_path, access_key="user", secret_key="pass")
        assert safe.exists()
        assert (safe / "manifest.safe").exists()

        backend = CopernicusDataspaceBackend(access_key="user", secret_key="pass")
        found_backend = backend.search(query)
        assert len(found_backend) == 1
        safe_backend = backend.download(found_backend[0], tmp_path)
        assert safe_backend.exists()

        import siac.adapters.data.gcs_sentinel2 as gcs_mod

        gcs_safe = tmp_path / f"{prod.product_id}.SAFE"
        gcs_safe.mkdir(parents=True, exist_ok=True)

        monkeypatch.setattr(gcs_mod, "search_gcs", lambda q: [prod])  # noqa: ARG005
        monkeypatch.setattr(gcs_mod, "download_gcs", lambda p, d: gcs_safe)  # noqa: ARG005

        assert gcs_mod.search_gcs(query) == [prod]
        assert gcs_mod.download_gcs(prod, tmp_path) == gcs_safe
        assert GCSSentinel2Backend().search(query) == [prod]
        assert GCSSentinel2Backend().download(prod, tmp_path) == gcs_safe


class TestCAMSProvider:
    def test_default_prior_when_missing(self, tmp_path: Path):
        p = CAMSProvider(tmp_path)
        state = p.get_prior((0.0, 0.0, 10.0, 10.0), "EPSG:4326", datetime(2024, 1, 1), 1.0)
        assert state.aot.shape == (10, 10)
        assert float(state.aot.mean()) == pytest.approx(0.15)
        assert float(state.tcwv.mean()) == pytest.approx(1.5)
        assert float(state.tco3.mean()) == pytest.approx(0.3)

    def test_default_prior_with_reversed_bounds(self, tmp_path: Path):
        p = CAMSProvider(tmp_path)
        state = p.get_prior((10.0, 10.0, 0.0, 0.0), "EPSG:4326", datetime(2024, 1, 1), 1.0)
        assert state.aot.shape == (10, 10)
        assert float(state.aot.mean()) == pytest.approx(0.15)

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
            xr.open_mfdataset = lambda *_args, **_kwargs: ds  # type: ignore[assignment]
            state = p.get_prior((0.0, -1.0, 1.0, 1.0), "EPSG:4326", datetime(2024, 1, 1, 1), 1.0)
        finally:
            xr.open_mfdataset = orig  # type: ignore[assignment]

        assert state.aot.size > 0
        assert float(state.aot.mean()) == pytest.approx(0.2)
        assert float(state.tcwv.mean()) == pytest.approx(0.22)
        expected_tco3 = 0.006 / 2.1415e-5 / 1000
        assert float(state.tco3.mean()) == pytest.approx(expected_tco3, rel=1e-6)
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

    def test_load_cams_from_direct_netcdf_file(self, tmp_path: Path, request: pytest.FixtureRequest):
        _skip_native_heavy_for_cdse_cov(request)
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

    def test_load_cams_from_direct_tif_file(self, tmp_path: Path, request: pytest.FixtureRequest):
        _skip_native_heavy_for_cdse_cov(request)
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
        assert float(state.tcwv.mean()) == pytest.approx(0.25)
        expected_tco3 = 0.006 / 2.1415e-5 / 1000
        assert float(state.tco3.mean()) == pytest.approx(expected_tco3, rel=1e-6)

    def test_single_band_tif_infers_variable_from_filename(
        self,
        tmp_path: Path,
        request: pytest.FixtureRequest,
    ):
        _skip_native_heavy_for_cdse_cov(request)
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
        assert float(state.tcwv.mean()) == pytest.approx(0.3)
        assert float(state.aot.mean()) == pytest.approx(0.15)
        assert float(state.tco3.mean()) == pytest.approx(0.3)


class TestKernelModel:
    def test_compute_surface_prior_without_psf(self, request: pytest.FixtureRequest):
        KernelModelDeriver, _ = _kernel_model_classes(request)
        geom = _geometry((6, 6))
        weights = _brdf_weights((6, 6))
        kmd = KernelModelDeriver(apply_psf=False)

        prior = kmd.compute_surface_prior(weights, geom)
        assert prior.boa.shape == (6, 6)
        assert prior.boa_unc.shape == (6, 6)
        assert prior.mask.dtype == bool

    def test_compute_surface_prior_with_psf_and_override(self, request: pytest.FixtureRequest):
        KernelModelDeriver, _ = _kernel_model_classes(request)
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

    def test_compute_surface_prior_handles_misaligned_coords(self, request: pytest.FixtureRequest):
        KernelModelDeriver, _ = _kernel_model_classes(request)
        geom = _geometry((8, 8))

        y_ref = np.linspace(0.0, 7.0, 4, dtype=np.float32)
        x_ref = np.linspace(0.0, 7.0, 4, dtype=np.float32)
        f0 = xr.DataArray(
            np.full((4, 4), 0.1, dtype=np.float32),
            dims=["y", "x"],
            coords={"y": y_ref, "x": x_ref},
        )
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

        kmd = KernelModelDeriver(apply_psf=False)
        prior = kmd.compute_surface_prior(weights, geom)

        assert prior.boa.shape == (4, 4)
        assert prior.boa_unc.shape == (4, 4)

    def test_convolution_helpers_handle_nan(self, request: pytest.FixtureRequest):
        KernelModelDeriver, _ = _kernel_model_classes(request)
        kmd = KernelModelDeriver(apply_psf=True)
        arr = np.array([[1.0, np.nan], [2.0, 3.0]], dtype=np.float32)
        out = kmd._convolve_2d(arr, sigma_x=1.0, sigma_y=1.0)
        assert out.shape == arr.shape
        assert np.isfinite(out[0, 0])

    def test_psf_convolver_numpy_and_xarray(self, request: pytest.FixtureRequest):
        _, PSFConvolver = _kernel_model_classes(request)
        c = PSFConvolver(sigma_x=0.8, sigma_y=1.2)
        arr = np.array([[1.0, 2.0], [3.0, np.nan]], dtype=np.float32)
        out_np = c.convolve(arr)
        assert out_np.shape == arr.shape

        da = _spatial_da(arr)
        out_da = c.convolve(da)
        assert isinstance(out_da, xr.DataArray)
        assert out_da.shape == da.shape
