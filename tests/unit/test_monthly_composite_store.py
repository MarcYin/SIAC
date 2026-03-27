from __future__ import annotations

import json
from typing import cast

import numpy as np
import pytest
import rasterio
import xarray as xr

from siac.algorithms.surface.brdf_monthly_composite import (
    MonthlyBestPixelComposite,
    MonthlyCompositeCollection,
    MonthlyKernelWeightComposite,
)
from siac.algorithms.surface.monthly_composite_store import (
    MonthlyCompositeStoreGridSpec,
    read_monthly_composite_collection,
    read_monthly_composite_store_manifest,
    write_monthly_composite_collection,
)
from siac.domain import SensorBand
from siac.runtime import BRDFKernelWeights


def _source_bands() -> tuple[SensorBand, ...]:
    return (
        SensorBand(
            "B02",
            490.0,
            65.0,
            10.0,
            0,
            rsrf_wavelengths_nm=np.array([470.0, 490.0, 510.0], dtype=np.float64),
            rsrf_response=np.array([0.1, 1.0, 0.1], dtype=np.float64),
            rsrf_sensor_unit_id="MSI",
            rsrf_representation_variant="native",
            rsrf_band_id="B02",
        ),
        SensorBand("B08", 842.0, 115.0, 10.0, 1),
    )


def _kernel_collection() -> MonthlyCompositeCollection:
    bands = ["B02", "B08"]
    coords = {"band": bands, "y": [0, 1], "x": [0]}
    cube = xr.DataArray(np.full((2, 2, 1), 0.2, dtype=np.float32), dims=["band", "y", "x"], coords=coords)
    quality = xr.DataArray(np.full((2, 1), 0.03, dtype=np.float32), dims=["y", "x"], coords={"y": [0, 1], "x": [0]})
    sample_index = xr.DataArray(np.zeros((2, 1), dtype=np.int16), dims=["y", "x"], coords={"y": [0, 1], "x": [0]})
    composite = MonthlyKernelWeightComposite(
        kernels=BRDFKernelWeights(
            f0=cube,
            f1=xr.zeros_like(cube),
            f2=xr.zeros_like(cube),
            f0_unc=xr.full_like(cube, 0.01),
            f1_unc=xr.full_like(cube, 0.01),
            f2_unc=xr.full_like(cube, 0.01),
        ),
        quality=quality,
        sample_index=sample_index,
        year=2023,
        month=7,
    )
    return MonthlyCompositeCollection(
        composites=(composite,),
        source_bands=_source_bands(),
        source_name="test-kernel-store",
    )


def _nan_kernel_collection() -> MonthlyCompositeCollection:
    bands = ["B02", "B08"]
    coords = {"band": bands, "y": [0, 1], "x": [0]}
    cube = xr.DataArray(np.full((2, 2, 1), np.nan, dtype=np.float32), dims=["band", "y", "x"], coords=coords)
    quality = xr.DataArray(np.full((2, 1), np.nan, dtype=np.float32), dims=["y", "x"], coords={"y": [0, 1], "x": [0]})
    sample_index = xr.DataArray(np.full((2, 1), -1, dtype=np.int16), dims=["y", "x"], coords={"y": [0, 1], "x": [0]})
    composite = MonthlyKernelWeightComposite(
        kernels=BRDFKernelWeights(
            f0=cube,
            f1=cube.copy(deep=True),
            f2=cube.copy(deep=True),
            f0_unc=cube.copy(deep=True),
            f1_unc=cube.copy(deep=True),
            f2_unc=cube.copy(deep=True),
            reflectance_unc=cube.copy(deep=True),
        ),
        quality=quality,
        sample_index=sample_index,
        year=2023,
        month=9,
    )
    return MonthlyCompositeCollection(
        composites=(composite,),
        source_bands=_source_bands(),
        source_name="test-nan-kernel-store",
    )


def _reflectance_collection() -> MonthlyCompositeCollection:
    composite = MonthlyBestPixelComposite(
        reflectance=xr.DataArray(
            np.array(
                [
                    [[0.11], [0.12]],
                    [[0.51], [0.52]],
                ],
                dtype=np.float32,
            ),
            dims=["band", "y", "x"],
            coords={"band": ["B02", "B08"], "y": [0, 1], "x": [0]},
        ),
        quality=xr.DataArray(np.full((2, 1), 0.05, dtype=np.float32), dims=["y", "x"], coords={"y": [0, 1], "x": [0]}),
        sample_index=xr.DataArray(np.ones((2, 1), dtype=np.int16), dims=["y", "x"], coords={"y": [0, 1], "x": [0]}),
        year=2024,
        month=8,
    )
    return MonthlyCompositeCollection(
        composites=(composite,),
        source_bands=_source_bands(),
        source_name="test-reflectance-store",
    )


def _georeferenced_kernel_collection() -> MonthlyCompositeCollection:
    bands = ["B02", "B08"]
    x = [400210.0, 400710.0]
    y = [4699790.0, 4699290.0]
    transform = rasterio.transform.from_bounds(399960.0, 4699040.0, 400960.0, 4700040.0, 2, 2)
    def _with_geo(data: xr.DataArray) -> xr.DataArray:
        return data.rio.set_spatial_dims(x_dim="x", y_dim="y").rio.write_crs("EPSG:32615").rio.write_transform(transform)

    cube = _with_geo(xr.DataArray(
        np.array(
            [
                [[1.0, 2.0], [3.0, 4.0]],
                [[10.0, 20.0], [30.0, 40.0]],
            ],
            dtype=np.float32,
        ),
        dims=["band", "y", "x"],
        coords={"band": bands, "y": y, "x": x},
    ))
    quality = _with_geo(xr.DataArray(
        np.array([[0.1, 0.2], [0.3, 0.4]], dtype=np.float32),
        dims=["y", "x"],
        coords={"y": y, "x": x},
    ))
    sample_index = _with_geo(xr.DataArray(
        np.array([[0, 1], [2, 3]], dtype=np.int16),
        dims=["y", "x"],
        coords={"y": y, "x": x},
    ))
    composite = MonthlyKernelWeightComposite(
        kernels=BRDFKernelWeights(
            f0=cube,
            f1=_with_geo((cube + 100.0).astype(np.float32)),
            f2=_with_geo((cube + 200.0).astype(np.float32)),
            f0_unc=_with_geo(xr.full_like(cube, 0.01)),
            f1_unc=_with_geo(xr.full_like(cube, 0.02)),
            f2_unc=_with_geo(xr.full_like(cube, 0.03)),
        ),
        quality=quality,
        sample_index=sample_index,
        year=2023,
        month=7,
    )
    return MonthlyCompositeCollection(
        composites=(composite,),
        source_bands=_source_bands(),
        source_name="test-georeferenced-kernel-store",
    )


def test_monthly_composite_store_round_trips_kernel_collection(tmp_path) -> None:
    store_path = write_monthly_composite_collection(_kernel_collection(), tmp_path / "kernel_store")
    loaded = read_monthly_composite_collection(store_path)
    manifest = read_monthly_composite_store_manifest(store_path)

    assert loaded.source_name == "test-kernel-store"
    assert manifest.version == 3
    assert manifest.entries[0].path == "2023_07"
    assert manifest.entries[0].format == "geotiff"
    assert manifest.entries[0].assets is not None
    assert manifest.entries[0].assets["f0"] == "f0.tif"
    assert manifest.entries[0].assets["sample_index"] == "sample_index.tif"
    assert [band.name for band in loaded.source_bands] == ["B02", "B08"]
    assert loaded.source_bands[0].has_rsrf
    np.testing.assert_allclose(
        loaded.source_bands[0].rsrf_wavelengths_nm,
        np.array([470.0, 490.0, 510.0], dtype=np.float64),
    )
    np.testing.assert_allclose(
        loaded.source_bands[0].rsrf_response,
        np.array([0.1, 1.0, 0.1], dtype=np.float64),
    )
    assert loaded.source_bands[0].rsrf_sensor_unit_id == "MSI"
    assert loaded.source_bands[0].rsrf_representation_variant == "native"
    assert loaded.source_bands[0].rsrf_band_id == "B02"
    composite = loaded.composites[0]
    assert isinstance(composite, MonthlyKernelWeightComposite)
    np.testing.assert_allclose(composite.kernels.f0.values, 0.2)
    np.testing.assert_allclose(composite.quality.values, 0.03)
    with rasterio.open(store_path / manifest.entries[0].path / manifest.entries[0].assets["f0"]) as src:
        assert src.driver == "GTiff"
        assert src.is_tiled
        assert str(src.profile["compress"]).lower() == "deflate"


def test_monthly_composite_store_honors_explicit_grid_metadata(tmp_path) -> None:
    grid = MonthlyCompositeStoreGridSpec.from_bounds(
        (399960.0, 4590240.0, 509760.0, 4700040.0),
        crs="EPSG:32615",
        resolution=500.0,
    )
    store_path = write_monthly_composite_collection(
        _georeferenced_kernel_collection(),
        tmp_path / "kernel_store_on_target_grid",
        grid=grid,
    )
    manifest = read_monthly_composite_store_manifest(store_path)

    assert manifest.grid == grid
    with rasterio.open(store_path / manifest.entries[0].path / manifest.entries[0].assets["f0"]) as src:
        assert src.crs is not None
        assert src.crs.to_string() == "EPSG:32615"
        assert src.width == grid.width
        assert src.height == grid.height
        assert src.bounds == pytest.approx(grid.bounds)
        assert src.transform == rasterio.transform.from_bounds(*grid.bounds, grid.width, grid.height)
        b02 = src.read(1)
        np.testing.assert_allclose(b02[:2, :2], np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32))
        assert int(np.isfinite(b02).sum()) == 4

    loaded = read_monthly_composite_collection(store_path)
    composite = cast("MonthlyKernelWeightComposite", loaded.composites[0])
    assert composite.kernels.f0.sizes["x"] == grid.width
    assert composite.kernels.f0.sizes["y"] == grid.height
    assert composite.sample_index.rio.crs is not None
    assert composite.sample_index.rio.crs.to_string() == "EPSG:32615"
    assert composite.sample_index.rio.transform(recalc=True) == rasterio.transform.from_bounds(
        *grid.bounds,
        grid.width,
        grid.height,
    )
    np.testing.assert_allclose(
        np.asarray(composite.kernels.f0.coords["x"].values, dtype=np.float64),
        np.linspace(
            grid.bounds[0] + grid.resolution / 2.0,
            grid.bounds[2] - grid.resolution / 2.0,
            grid.width,
        ),
    )
    np.testing.assert_allclose(
        np.asarray(composite.kernels.f0.coords["y"].values, dtype=np.float64),
        np.linspace(
            grid.bounds[3] - grid.resolution / 2.0,
            grid.bounds[1] + grid.resolution / 2.0,
            grid.height,
        ),
    )
    loaded_b02 = composite.kernels.f0.sel(band="B02").values
    np.testing.assert_allclose(loaded_b02[:2, :2], np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32))
    assert int(np.isfinite(loaded_b02).sum()) == 4


def test_monthly_composite_store_round_trips_reflectance_collection(tmp_path) -> None:
    store_path = write_monthly_composite_collection(_reflectance_collection(), tmp_path / "reflectance_store")
    loaded = read_monthly_composite_collection(store_path)
    manifest = read_monthly_composite_store_manifest(store_path)

    assert loaded.source_name == "test-reflectance-store"
    assert manifest.entries[0].path == "2024_08"
    assert manifest.entries[0].assets is not None
    assert manifest.entries[0].assets["reflectance"] == "reflectance.tif"
    composite = loaded.composites[0]
    assert isinstance(composite, MonthlyBestPixelComposite)
    np.testing.assert_allclose(
        composite.reflectance.sel(band="B02").values,
        np.array([[0.11], [0.12]], dtype=np.float32),
    )


def test_monthly_composite_store_manifest_rejects_unsupported_version(tmp_path) -> None:
    root = tmp_path / "bad_store"
    root.mkdir()
    (root / "manifest.json").write_text(json.dumps({"version": 99, "entries": []}), encoding="utf-8")

    with pytest.raises(ValueError, match="Unsupported monthly composite store version"):
        read_monthly_composite_store_manifest(root)


def test_monthly_composite_store_skips_all_nan_kernel_month(tmp_path) -> None:
    store_path = write_monthly_composite_collection(_nan_kernel_collection(), tmp_path / "nan_kernel_store")

    manifest = read_monthly_composite_store_manifest(store_path)
    assert manifest.entries == ()
    assert not (store_path / "2023_09").exists()

    loaded = read_monthly_composite_collection(store_path)
    assert loaded.composites == ()


def test_monthly_composite_store_removes_stale_entry_when_new_run_skips_nan_month(tmp_path) -> None:
    store_path = tmp_path / "nan_kernel_store"
    write_monthly_composite_collection(_kernel_collection(), store_path)
    stale_month = store_path / "2023_07"
    assert stale_month.exists()

    write_monthly_composite_collection(_nan_kernel_collection(), store_path)

    manifest = read_monthly_composite_store_manifest(store_path)
    assert manifest.entries == ()
    assert not stale_month.exists()


def test_monthly_composite_store_rejects_irregular_spatial_axes(tmp_path) -> None:
    collection = _reflectance_collection()
    composite = cast("MonthlyBestPixelComposite", collection.composites[0])
    irregular = MonthlyBestPixelComposite(
        reflectance=xr.DataArray(
            np.array(
                [
                    [[0.11, 0.12, 0.13], [0.14, 0.15, 0.16]],
                    [[0.51, 0.52, 0.53], [0.54, 0.55, 0.56]],
                ],
                dtype=np.float32,
            ),
            dims=["band", "y", "x"],
            coords={"band": ["B02", "B08"], "y": [0.0, 1.0], "x": [0.0, 1.5, 3.5]},
        ),
        quality=xr.DataArray(
            np.full((2, 3), 0.05, dtype=np.float32),
            dims=["y", "x"],
            coords={"y": [0.0, 1.0], "x": [0.0, 1.5, 3.5]},
        ),
        sample_index=xr.DataArray(
            np.ones((2, 3), dtype=np.int16),
            dims=["y", "x"],
            coords={"y": [0.0, 1.0], "x": [0.0, 1.5, 3.5]},
        ),
        year=composite.year,
        month=composite.month,
    )
    collection = MonthlyCompositeCollection(
        composites=(irregular,),
        source_bands=collection.source_bands,
        source_name=collection.source_name,
    )

    with pytest.raises(ValueError, match="regularly spaced x coordinates"):
        write_monthly_composite_collection(collection, tmp_path / "irregular_store")
