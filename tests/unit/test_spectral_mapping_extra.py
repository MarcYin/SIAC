from __future__ import annotations

from types import SimpleNamespace
from typing import TYPE_CHECKING

import numpy as np
import pytest
import rasterio
import xarray as xr

if TYPE_CHECKING:
    from pathlib import Path

import siac.algorithms.surface.spectral_mapping as spectral_mapping_mod
from siac.algorithms.surface import _spectral_curve_utils as curve_utils
from siac.algorithms.surface.spectral_mapping import SpectralMapper, SpectralMappingConfig
from siac.domain import SensorBand


def _source_bands() -> tuple[SensorBand, ...]:
    return (
        SensorBand(
            "Band3", 469.0, 20.0, 500.0, 0, rsrf_sensor_unit_id="terra_modis", rsrf_band_id="B3"
        ),
        SensorBand(
            "Band4", 555.0, 20.0, 500.0, 1, rsrf_sensor_unit_id="terra_modis", rsrf_band_id="B4"
        ),
        SensorBand(
            "Band2", 858.5, 35.0, 500.0, 2, rsrf_sensor_unit_id="terra_modis", rsrf_band_id="B2"
        ),
        SensorBand(
            "Band6", 1640.0, 24.0, 500.0, 3, rsrf_sensor_unit_id="terra_modis", rsrf_band_id="B6"
        ),
        SensorBand(
            "Band7", 2130.0, 50.0, 500.0, 4, rsrf_sensor_unit_id="terra_modis", rsrf_band_id="B7"
        ),
    )


def _sampled_band(
    name: str,
    center_wavelength: float,
    bandwidth: float,
    resolution: float,
    band_index: int,
    *,
    sensor_unit_id: str,
    rsrf_band_id: str,
) -> SensorBand:
    wavelengths = np.arange(
        center_wavelength - 3.0 * bandwidth,
        center_wavelength + 3.0 * bandwidth + 1.0,
        1.0,
        dtype=np.float32,
    )
    sigma = bandwidth / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    response = np.exp(-0.5 * ((wavelengths - center_wavelength) / sigma) ** 2).astype(np.float32)
    return SensorBand(
        name,
        center_wavelength,
        bandwidth,
        resolution,
        band_index,
        rsrf_wavelengths_nm=wavelengths,
        rsrf_response=response,
        rsrf_sensor_unit_id=sensor_unit_id,
        rsrf_band_id=rsrf_band_id,
    )


def _target_bands() -> tuple[SensorBand, ...]:
    return (
        _sampled_band(
            "B02", 490.0, 65.0, 10.0, 0, sensor_unit_id="sentinel-2a_msi", rsrf_band_id="B02"
        ),
        _sampled_band(
            "B03", 560.0, 35.0, 10.0, 1, sensor_unit_id="sentinel-2a_msi", rsrf_band_id="B03"
        ),
        _sampled_band(
            "B08", 842.0, 115.0, 10.0, 2, sensor_unit_id="sentinel-2a_msi", rsrf_band_id="B08"
        ),
        _sampled_band(
            "B11", 1610.0, 90.0, 20.0, 3, sensor_unit_id="sentinel-2a_msi", rsrf_band_id="B11"
        ),
        _sampled_band(
            "B12", 2190.0, 180.0, 20.0, 4, sensor_unit_id="sentinel-2a_msi", rsrf_band_id="B12"
        ),
    )


def _gaussian_target_band() -> SensorBand:
    return SensorBand("T01", 705.0, 20.0, 10.0, 0)


def _with_geo(data: xr.DataArray) -> xr.DataArray:
    height = int(data.sizes["y"])
    width = int(data.sizes["x"])
    resolution = 10.0
    xmin = 500000.0
    ymax = 4100000.0
    x = np.linspace(
        xmin + resolution / 2.0,
        xmin + width * resolution - resolution / 2.0,
        width,
        dtype=np.float64,
    )
    y = np.linspace(
        ymax - resolution / 2.0,
        ymax - height * resolution + resolution / 2.0,
        height,
        dtype=np.float64,
    )
    transform = rasterio.transform.from_origin(xmin, ymax, resolution, resolution)
    out = data.assign_coords({"x": x, "y": y}).rio.set_spatial_dims(x_dim="x", y_dim="y")
    return out.rio.write_crs("EPSG:32632").rio.write_transform(transform)


def test_split_mapping_inputs_accepts_config_and_rejects_invalid_input() -> None:
    config = spectral_mapping_mod._split_mapping_inputs(SpectralMappingConfig(knn_backend="numpy"))
    assert isinstance(config, SpectralMappingConfig)
    assert config.knn_backend == "numpy"

    with pytest.raises(TypeError, match="SpectralMappingConfig"):
        spectral_mapping_mod._split_mapping_inputs({"extra": True})


def test_source_schema_indices_ignore_unpublished_source_bands() -> None:
    schema_band_ids = ("blue", "green", "red", "nir", "swir1", "swir2")
    source_bands = _source_bands() + (
        SensorBand(
            "Band5", 1240.0, 20.0, 500.0, 5, rsrf_sensor_unit_id="terra_modis", rsrf_band_id="B5"
        ),
    )

    indices, ignored = spectral_mapping_mod._source_schema_indices_for_bands(
        source_bands,
        sensor_id="terra_modis",
        source_schema_band_ids=schema_band_ids,
    )

    assert indices == (0, 1, 3, 4, 5, None)
    assert ignored == ("Band5",)


def test_target_output_band_ids_resolve_native_target_bands() -> None:
    target_output_band_ids = spectral_mapping_mod._target_output_band_ids_for_bands(
        _target_bands(),
        sensor_id="sentinel-2a_msi",
        target_schema_band_ids=("ultra_blue", "blue", "green", "red", "nir", "swir1", "swir2"),
    )
    assert target_output_band_ids == ("blue", "green", "nir", "swir1", "swir2")


def test_map_identity_and_nonidentity_helpers_cover_runtime_branches(tmp_path: Path) -> None:
    source = xr.DataArray(
        np.array([[[0.2]], [[0.3]]], dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": ["B02", "B03"], "y": [0], "x": [0]},
    )

    identity = object.__new__(SpectralMapper)
    identity.source_bands = (
        SensorBand("B02", 490.0, 65.0, 10.0, 0),
        SensorBand("B03", 560.0, 35.0, 10.0, 1),
    )
    identity.target_bands = identity.source_bands
    identity.k_neighbors = 1
    identity._identity = True
    identity._target_band_names = ["B02", "B03"]
    mapped, mapped_unc, source_fit = identity.map(source, source_uncertainty=None)
    np.testing.assert_allclose(mapped.values, source.values)
    np.testing.assert_allclose(mapped_unc.values, np.full_like(source.values, 0.005))
    np.testing.assert_allclose(source_fit.values, np.zeros((1, 1), dtype=np.float32))

    nonidentity = object.__new__(SpectralMapper)
    nonidentity.source_bands = _source_bands()
    nonidentity.target_bands = (_target_bands()[0],)
    nonidentity.k_neighbors = 1
    nonidentity._identity = False
    nonidentity._target_band_names = ["B02"]
    nonidentity._mapping_config = SpectralMappingConfig()
    nonidentity._runtime = SimpleNamespace(
        prepared_root=tmp_path / "runtime" / "prepared",
        source_sensor_id="terra_modis",
        source_schema_band_ids=("blue", "green", "red", "nir", "swir1", "swir2"),
        source_schema_indices=(0, 1, 3, 4, 5),
        target_sensor_id="sentinel-2a_msi",
        target_output_band_ids=("blue",),
        ignored_source_band_names=(),
    )
    nonidentity._supported_source_input_indices = (0, 1, 2, 3, 4)
    nonidentity._package_mapper = SimpleNamespace(
        map_reflectance_batch_arrays_ndarray=lambda **_kwargs: SimpleNamespace(
            reflectance=np.full((1, 1), 0.4, dtype=np.float64),
            source_fit_rmse=np.array([0.01], dtype=np.float32),
            output_columns=("blue",),
        )
    )

    single = xr.DataArray(
        np.array([[[0.25]], [[0.22]], [[0.35]], [[0.18]], [[0.11]]], dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": [band.name for band in _source_bands()], "y": [0], "x": [0]},
    )
    mapped, mapped_unc, source_fit = nonidentity.map(single)
    assert float(mapped.values[0, 0, 0]) == pytest.approx(0.4)
    assert float(mapped_unc.values[0, 0, 0]) > 0.005
    assert np.isfinite(source_fit.values[0, 0])

    uninitialized = object.__new__(SpectralMapper)
    uninitialized.source_bands = nonidentity.source_bands
    uninitialized.target_bands = nonidentity.target_bands
    uninitialized._identity = False
    uninitialized._package_mapper = None
    uninitialized._runtime = None
    with pytest.raises(RuntimeError, match="not initialized"):
        uninitialized.map(single)


def test_package_query_rows_expand_to_published_source_band_order() -> None:
    mapper = object.__new__(SpectralMapper)
    mapper._runtime = SimpleNamespace(
        source_schema_band_ids=("blue", "green", "red", "nir", "swir1", "swir2"),
        source_schema_indices=(0, 1, 3, 4, 5),
    )
    source = xr.DataArray(
        np.array([[[0.2]], [[0.3]], [[0.4]], [[0.5]], [[0.6]]], dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": [band.name for band in _source_bands()], "y": [0], "x": [0]},
    )

    flattened = SpectralMapper._flatten_source_cube(mapper, source)
    query_rows, valid_masks, valid_rows = SpectralMapper._package_query_rows(mapper, flattened)

    assert query_rows.shape == (1, 6)
    np.testing.assert_allclose(
        query_rows[0],
        np.array([0.2, 0.3, np.nan, 0.4, 0.5, 0.6], dtype=np.float64),
        equal_nan=True,
    )
    assert valid_masks[0].tolist() == [True, True, False, True, True, True]
    assert valid_rows.tolist() == [True]


def test_map_does_not_cache_distance_metrics_under_upstream_runtime(tmp_path: Path) -> None:
    mapper = object.__new__(SpectralMapper)
    mapper.source_bands = _source_bands()
    mapper.target_bands = (_target_bands()[0],)
    mapper.k_neighbors = 1
    mapper._identity = False
    mapper._target_band_names = ["B02"]
    mapper._mapping_config = SpectralMappingConfig()
    mapper._runtime = SimpleNamespace(
        source_sensor_id="terra_modis",
        ignored_source_band_names=("Band5",),
        prepared_root=tmp_path / "runtime" / "prepared",
        source_schema_band_ids=("blue", "green", "red", "nir", "swir1", "swir2"),
        source_schema_indices=(0, 1, 3, 4, 5),
        target_sensor_id="sentinel-2a_msi",
        target_output_band_ids=("blue",),
    )
    mapper._supported_source_input_indices = (0, 1, 2, 3, 4)
    mapper._package_mapper = SimpleNamespace(
        map_reflectance_batch_arrays_ndarray=lambda **_kwargs: SimpleNamespace(
            reflectance=np.full((1, 1), 0.4, dtype=np.float64),
            source_fit_rmse=np.array([0.05], dtype=np.float32),
            output_columns=("blue",),
        ),
    )

    source = _with_geo(
        xr.DataArray(
            np.array([[[0.25]], [[0.24]], [[0.23]], [[0.22]], [[0.21]]], dtype=np.float32),
            dims=["band", "y", "x"],
            coords={"band": [band.name for band in _source_bands()], "y": [0], "x": [0]},
        )
    )

    mapper.map(source)

    diagnostics_dir = tmp_path / "runtime" / "diagnostics"
    assert not diagnostics_dir.exists()


def test_curve_utils_cover_classification_and_error_paths() -> None:
    visible = SensorBand("B03", 560.0, 35.0, 10.0, 0)
    swir = SensorBand("B11", 1610.0, 90.0, 20.0, 1)
    assert curve_utils.classify_band_region(visible) == "visible"
    assert curve_utils.classify_band_region(swir) == "infrared"
    assert curve_utils.segment_for_band(swir) == "swir"
    assert curve_utils.primary_nir_band_index((visible, swir)) is None

    with pytest.raises(ValueError, match="fwhm_nm must be positive"):
        curve_utils.normalized_band_response(
            SensorBand("X", 490.0, 0.0, 10.0, 0),
            np.array([490.0, 491.0], dtype=np.float32),
        )

    with pytest.raises(ValueError, match="at least two wavelength samples"):
        curve_utils.canonicalize_curve(
            np.array([490.0], dtype=np.float32),
            np.array([1.0], dtype=np.float32),
        )

    with pytest.raises(ValueError, match="positive response sample"):
        curve_utils.canonicalize_curve(
            np.array([490.0, 491.0], dtype=np.float32),
            np.array([0.0, 0.0], dtype=np.float32),
        )
