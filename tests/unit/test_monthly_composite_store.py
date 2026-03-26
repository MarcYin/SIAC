from __future__ import annotations

import json

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
