from __future__ import annotations

import json
from types import SimpleNamespace
from typing import TYPE_CHECKING

import numpy as np
import pytest
import rasterio
import xarray as xr

import siac.algorithms.surface.spectral_mapping as spectral_mapping_mod
from siac.algorithms.surface import _spectral_curve_utils as curve_utils
from siac.algorithms.surface.spectral_mapping import HyperspectralLibrary, SpectralMapper
from siac.domain import SensorBand

if TYPE_CHECKING:
    from pathlib import Path


def _library() -> HyperspectralLibrary:
    wavelengths = np.arange(390.0, 2511.0, 10.0, dtype=np.float32)
    spectra = np.stack(
        [
            np.linspace(0.05, 0.4, wavelengths.size, dtype=np.float32),
            np.linspace(0.15, 0.3, wavelengths.size, dtype=np.float32),
        ],
        axis=0,
    )
    return HyperspectralLibrary(
        wavelengths_nm=wavelengths,
        spectra=spectra,
        sample_ids=("a", "b"),
        source_id="extra-library",
    )


def _source_bands() -> tuple[SensorBand, ...]:
    return (
        SensorBand("B03", 560.0, 35.0, 10.0, 0),
        SensorBand("B08", 842.0, 115.0, 10.0, 1),
        SensorBand("B11", 1610.0, 90.0, 20.0, 2),
        SensorBand("B12", 2190.0, 180.0, 20.0, 3),
    )


def _target_bands() -> tuple[SensorBand, ...]:
    return (
        SensorBand("B02", 490.0, 65.0, 10.0, 0),
        SensorBand("B03", 560.0, 35.0, 10.0, 1),
        SensorBand("B08", 842.0, 115.0, 10.0, 2),
        SensorBand("B11", 1610.0, 90.0, 20.0, 3),
        SensorBand("B12", 2190.0, 180.0, 20.0, 4),
    )


def _landsat_source_subset_bands() -> tuple[SensorBand, ...]:
    return (
        SensorBand(
            "B2",
            482.0,
            60.0,
            30.0,
            0,
            rsrf_sensor_unit_id="landsat-8_oli",
            rsrf_band_id="B2",
        ),
        SensorBand(
            "B4",
            654.5,
            37.0,
            30.0,
            1,
            rsrf_sensor_unit_id="landsat-8_oli",
            rsrf_band_id="B4",
        ),
        SensorBand(
            "B6",
            1608.5,
            85.0,
            30.0,
            2,
            rsrf_sensor_unit_id="landsat-8_oli",
            rsrf_band_id="B6",
        ),
    )


def _landsat_target_subset_bands() -> tuple[SensorBand, ...]:
    return (
        SensorBand(
            "B1",
            443.0,
            16.0,
            30.0,
            0,
            rsrf_sensor_unit_id="landsat-8_oli",
            rsrf_band_id="B1",
        ),
        SensorBand(
            "B2",
            482.0,
            60.0,
            30.0,
            1,
            rsrf_sensor_unit_id="landsat-8_oli",
            rsrf_band_id="B2",
        ),
        SensorBand(
            "B4",
            654.5,
            37.0,
            30.0,
            2,
            rsrf_sensor_unit_id="landsat-8_oli",
            rsrf_band_id="B4",
        ),
        SensorBand(
            "B5",
            865.0,
            28.0,
            30.0,
            3,
            rsrf_sensor_unit_id="landsat-8_oli",
            rsrf_band_id="B5",
        ),
        SensorBand(
            "B6",
            1608.5,
            85.0,
            30.0,
            4,
            rsrf_sensor_unit_id="landsat-8_oli",
            rsrf_band_id="B6",
        ),
        SensorBand(
            "B7",
            2200.5,
            187.0,
            30.0,
            5,
            rsrf_sensor_unit_id="landsat-8_oli",
            rsrf_band_id="B7",
        ),
    )


def _sentinel_noncanonical_target_bands() -> tuple[SensorBand, ...]:
    return (
        SensorBand(
            "B02",
            490.0,
            65.0,
            10.0,
            0,
            rsrf_sensor_unit_id="sentinel-2a_msi",
            rsrf_band_id="B02",
        ),
        SensorBand(
            "B08",
            842.0,
            115.0,
            10.0,
            1,
            rsrf_sensor_unit_id="sentinel-2a_msi",
            rsrf_band_id="B08",
        ),
        SensorBand(
            "B11",
            1610.0,
            90.0,
            20.0,
            2,
            rsrf_sensor_unit_id="sentinel-2a_msi",
            rsrf_band_id="B11",
        ),
    )


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


def test_runtime_input_helpers_cover_sensor_inputs_and_root_loading(tmp_path: Path) -> None:
    shared_rsrf = (
        SensorBand("B02", 490.0, 65.0, 10.0, 0, rsrf_sensor_unit_id="sentinel-2a_msi"),
        SensorBand("B03", 560.0, 35.0, 10.0, 1, rsrf_sensor_unit_id="sentinel-2a_msi"),
    )
    sensor_input, band_ids = spectral_mapping_mod._sensor_input_for_bands(shared_rsrf)
    assert isinstance(sensor_input, spectral_mapping_mod.SensorInput)
    assert band_ids == ("B02", "B03")
    assert sensor_input.bands[0].rsrf_sensor_id == "sentinel-2a_msi"
    assert sensor_input.bands[0].rsrf_band_id == "B02"

    package_library = spectral_mapping_mod._package_library_input_from_hyperspectral_library(
        _library()
    )
    assert isinstance(package_library, spectral_mapping_mod.PackageHyperspectralLibraryInput)
    assert package_library.sample_ids == ("a", "b")
    assert float(package_library.wavelengths_nm[0]) == 400.0
    assert float(package_library.wavelengths_nm[-1]) == 2500.0
    assert package_library.metadata_rows == [
        {"source_version": "1"},
        {"source_version": "1"},
    ]

    root = tmp_path / "library"
    tabular = root / "tabular"
    tabular.mkdir(parents=True)
    (tabular / "siac_spectra_metadata.csv").write_text(
        "source_id,spectrum_id,sample_name\nroot,alpha,Alpha\nroot,beta,Beta\n",
        encoding="utf-8",
    )
    (tabular / "siac_normalized_spectra.csv").write_text(
        "source_id,spectrum_id,sample_name,nm_400,nm_401\n"
        "root,alpha,Alpha,0.1,0.2\n"
        "root,beta,Beta,0.3,0.4\n",
        encoding="utf-8",
    )
    loaded = spectral_mapping_mod._package_library_input_from_root(root)
    assert loaded.sample_ids == ("alpha", "beta")
    np.testing.assert_allclose(loaded.wavelengths_nm, np.array([400.0, 401.0], dtype=np.float64))
    np.testing.assert_allclose(
        loaded.spectra,
        np.array([[0.1, 0.2], [0.3, 0.4]], dtype=np.float32),
    )
    assert loaded.metadata_rows == [
        {"source_id": "root", "spectrum_id": "alpha", "sample_name": "Alpha"},
        {"source_id": "root", "spectrum_id": "beta", "sample_name": "Beta"},
    ]


def test_package_library_input_rejects_missing_canonical_coverage(
    tmp_path: Path,
) -> None:
    bad = HyperspectralLibrary(
        wavelengths_nm=np.arange(450.0, 2501.0, 10.0, dtype=np.float32),
        spectra=np.stack(
            [
                np.linspace(0.1, 0.2, 206, dtype=np.float32),
                np.linspace(0.2, 0.3, 206, dtype=np.float32),
            ],
            axis=0,
        ),
        sample_ids=("a", "b"),
    )

    with pytest.raises(ValueError, match="canonical 400-2500 nm range"):
        spectral_mapping_mod._package_library_input_from_hyperspectral_library(bad)


def test_prepare_runtime_uses_in_memory_runtime_builder(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    calls: dict[str, object] = {}

    class _FakePreparedRuntime:
        def __init__(self) -> None:
            self.prepared_root = tmp_path / "cache" / "runtime"
            self.source_sensor_ids = ("custom_source",)
            self.target_sensor_ids = ("custom_target",)
            self.target_band_ids = {"custom_target": ("B1", "B2", "B4", "nir", "B6", "B7")}

    def _fake_build_mapping_runtime(
        *,
        library,
        source_sensors,
        target_sensors,
        cache_root,
        verify_checksums,
    ):  # noqa: ANN001
        calls["library"] = library
        calls["source_sensors"] = list(source_sensors)
        calls["target_sensors"] = list(target_sensors)
        calls["cache_root"] = cache_root
        calls["verify_checksums"] = verify_checksums
        return _FakePreparedRuntime()

    monkeypatch.setattr(spectral_mapping_mod, "build_mapping_runtime", _fake_build_mapping_runtime)

    siac_root = tmp_path / "siac-library"
    (siac_root / "tabular").mkdir(parents=True)
    (siac_root / "tabular" / "siac_spectra_metadata.csv").write_text(
        "source_id,spectrum_id,sample_name\nroot,a,A\n",
        encoding="utf-8",
    )
    (siac_root / "tabular" / "siac_normalized_spectra.csv").write_text(
        "source_id,spectrum_id,sample_name,"
        + ",".join(f"nm_{int(wl)}" for wl in spectral_mapping_mod._CANONICAL_WAVELENGTHS_NM)
        + "\nroot,a,A,"
        + ",".join("0.1" for _ in spectral_mapping_mod._CANONICAL_WAVELENGTHS_NM)
        + "\n",
        encoding="utf-8",
    )
    runtime = spectral_mapping_mod._prepare_runtime(
        _landsat_source_subset_bands(),
        _landsat_target_subset_bands(),
        library=None,
        config=spectral_mapping_mod.SpectralMappingConfig(
            cache_dir=tmp_path,
            siac_library_root=siac_root,
        ),
    )

    assert calls["cache_root"] == tmp_path
    assert calls["verify_checksums"] is False
    assert isinstance(calls["library"], spectral_mapping_mod.PackageHyperspectralLibraryInput)
    assert isinstance(calls["source_sensors"][0], spectral_mapping_mod.SensorInput)
    assert isinstance(calls["target_sensors"][0], spectral_mapping_mod.SensorInput)
    assert calls["source_sensors"][0].bands[0].rsrf_sensor_id == "landsat-8_oli"
    assert calls["target_sensors"][0].bands[3].band_id == "nir"
    assert runtime.source_sensor_id == "custom_source"
    assert runtime.target_sensor_id == "custom_target"
    assert runtime.target_band_ids == ("B1", "B2", "B4", "nir", "B6", "B7")


def test_sensor_input_uses_requested_noncanonical_target_subset() -> None:
    sensor_input, band_ids = spectral_mapping_mod._sensor_input_for_bands(
        _sentinel_noncanonical_target_bands()
    )

    assert isinstance(sensor_input, spectral_mapping_mod.SensorInput)
    assert band_ids == ("B02", "nir", "B11")
    assert sensor_input.bands[1].band_id == "nir"
    assert sensor_input.bands[1].rsrf_band_id == "B08"


def test_package_query_rows_preserve_source_band_order() -> None:
    mapper = object.__new__(SpectralMapper)
    source = xr.DataArray(
        np.array([[[0.2]], [[0.3]], [[0.4]]], dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": ["B2", "B4", "B6"], "y": [0], "x": [0]},
    )

    flattened = SpectralMapper._flatten_source_cube(mapper, source)
    query_rows, valid_masks, valid_rows = SpectralMapper._package_query_rows(mapper, flattened)

    assert query_rows.shape == (1, 3)
    np.testing.assert_allclose(query_rows[0], np.array([0.2, 0.3, 0.4], dtype=np.float64))
    assert valid_masks[0].tolist() == [True, True, True]
    assert valid_rows.tolist() == [True]


def test_map_identity_and_nonidentity_helpers_cover_runtime_branches() -> None:
    source = xr.DataArray(
        np.array([[[0.2]], [[0.3]]], dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": ["B02", "B03"], "y": [0], "x": [0]},
    )
    unc = xr.full_like(source, 0.02)

    identity = object.__new__(SpectralMapper)
    identity.source_bands = (
        SensorBand("B02", 490.0, 65.0, 10.0, 0),
        SensorBand("B03", 560.0, 35.0, 10.0, 1),
    )
    identity.target_bands = identity.source_bands
    identity._identity = True
    identity._target_band_names = ["B02", "B03"]
    mapped, mapped_unc, source_fit = identity.map(source, source_uncertainty=None)
    np.testing.assert_allclose(mapped.values, source.values)
    np.testing.assert_allclose(mapped_unc.values, np.full_like(source.values, 0.005))
    np.testing.assert_allclose(source_fit.values, np.zeros((1, 1), dtype=np.float32))

    nonidentity = object.__new__(SpectralMapper)
    nonidentity.source_bands = (SensorBand("B02", 490.0, 65.0, 10.0, 0),)
    nonidentity.target_bands = (SensorBand("T01", 500.0, 50.0, 10.0, 0),)
    nonidentity.k_neighbors = 1
    nonidentity._identity = False
    nonidentity._target_band_names = ["T01"]
    nonidentity._mapping_config = SimpleNamespace(
        min_valid_bands=1,
        neighbor_estimator="distance_weighted_mean",
        knn_backend="numpy",
        knn_eps=0.0,
    )
    nonidentity._runtime = SimpleNamespace(source_sensor_id="src", target_sensor_id="target")
    nonidentity._package_mapper = SimpleNamespace(
        map_reflectance_batch_arrays_ndarray=lambda **_kwargs: SimpleNamespace(
            reflectance=np.array([[0.4]], dtype=np.float64),
            source_fit_rmse=np.array([0.01], dtype=np.float32),
            output_columns=("T01",),
        )
    )
    nonidentity._target_internal_to_output_index = {"T01": 0}

    single = xr.DataArray(
        np.array([[[0.25]]], dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": ["B02"], "y": [0], "x": [0]},
    )
    mapped, mapped_unc, source_fit = nonidentity.map(
        single, source_uncertainty=unc.sel(band=["B02"])
    )
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


def test_map_uses_batch_array_api_and_returns_source_fit_rmse() -> None:
    mapper = object.__new__(SpectralMapper)
    mapper.source_bands = (SensorBand("B02", 490.0, 65.0, 10.0, 0),)
    mapper.target_bands = (SensorBand("T01", 500.0, 50.0, 10.0, 0),)
    mapper.k_neighbors = 1
    mapper._identity = False
    mapper._target_band_names = ["T01"]
    mapper._mapping_config = SimpleNamespace(
        min_valid_bands=1,
        neighbor_estimator="distance_weighted_mean",
        knn_backend="numpy",
        knn_eps=0.0,
    )
    mapper._runtime = SimpleNamespace(source_sensor_id="src", target_sensor_id="target")
    mapper._target_internal_to_output_index = {"T01": 0}

    calls: dict[str, int] = {"batch_arrays": 0}

    def _mock_batch_arrays(**_kwargs):  # noqa: ANN003
        calls["batch_arrays"] += 1
        return SimpleNamespace(
            reflectance=np.array([[0.4]], dtype=np.float64),
            source_fit_rmse=np.array([0.05], dtype=np.float32),
            output_columns=("T01",),
        )

    mapper._package_mapper = SimpleNamespace(
        map_reflectance_batch_arrays_ndarray=_mock_batch_arrays,
    )

    source = xr.DataArray(
        np.array([[[0.25]]], dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": ["B02"], "y": [0], "x": [0]},
    )

    mapped, mapped_unc, source_fit = mapper.map(source)

    assert calls["batch_arrays"] == 1
    assert float(mapped.values[0, 0, 0]) == pytest.approx(0.4)
    assert float(mapped_unc.values[0, 0, 0]) > 0.005
    assert np.isfinite(source_fit.values[0, 0])
    assert float(source_fit.values[0, 0]) == pytest.approx(0.05)


def test_map_caches_distance_metrics_to_disk(tmp_path: Path) -> None:
    mapper = object.__new__(SpectralMapper)
    mapper.source_bands = (SensorBand("B02", 490.0, 65.0, 10.0, 0),)
    mapper.target_bands = (SensorBand("T01", 500.0, 50.0, 10.0, 0),)
    mapper.k_neighbors = 1
    mapper._identity = False
    mapper._target_band_names = ["T01"]
    mapper._mapping_config = SimpleNamespace(
        min_valid_bands=1,
        neighbor_estimator="distance_weighted_mean",
        knn_backend="numpy",
        knn_eps=0.0,
    )
    mapper._runtime = SimpleNamespace(
        source_sensor_id="src",
        target_sensor_id="target",
        prepared_root=tmp_path / "runtime" / "prepared",
    )
    mapper._target_internal_to_output_index = {"T01": 0}

    mapper._package_mapper = SimpleNamespace(
        map_reflectance_batch_arrays_ndarray=lambda **_kwargs: SimpleNamespace(
            reflectance=np.array([[0.4]], dtype=np.float64),
            source_fit_rmse=np.array([0.05], dtype=np.float32),
            output_columns=("T01",),
        ),
    )

    source = _with_geo(
        xr.DataArray(
            np.array([[[0.25]]], dtype=np.float32),
            dims=["band", "y", "x"],
            coords={"band": ["B02"], "y": [0], "x": [0]},
        )
    )

    mapper.map(source)

    diagnostics_dir = tmp_path / "runtime" / "diagnostics"
    data_paths = sorted(diagnostics_dir.glob("spectral_mapping_distances_*.tif"))
    metadata_paths = sorted(diagnostics_dir.glob("spectral_mapping_distances_*.json"))
    assert len(data_paths) == 1
    assert len(metadata_paths) == 1

    for path in data_paths:
        with rasterio.open(path) as src:
            assert src.crs.to_string() == "EPSG:32632"
            assert src.transform == source.isel(band=0, drop=True).rio.transform()

    metadata = json.loads(metadata_paths[0].read_text(encoding="utf-8"))
    assert metadata["source_sensor_id"] == "src"
    assert metadata["target_sensor_id"] == "target"
    assert sorted(metadata["metrics"]) == ["source_fit_rmse"]


def test_wrapper_and_input_split_helpers_cover_error_and_delegation_paths(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    assert spectral_mapping_mod._split_mapping_inputs(None)[0] is None
    with pytest.raises(TypeError, match="HyperspectralLibrary"):
        spectral_mapping_mod._split_mapping_inputs(object())

    captured: dict[str, object] = {}

    class _FakeMapper:
        def __init__(self, source_bands, target_bands, *, spectral_library, k_neighbors):  # noqa: ANN001
            captured["init"] = {
                "source_bands": source_bands,
                "target_bands": target_bands,
                "spectral_library": spectral_library,
                "k_neighbors": k_neighbors,
            }

        def map(self, source_reflectance, *, source_uncertainty=None):  # noqa: ANN001
            captured["map"] = {
                "source_reflectance": source_reflectance,
                "source_uncertainty": source_uncertainty,
            }
            return "mapped", "unc", "fit"

    monkeypatch.setattr(spectral_mapping_mod, "SpectralMapper", _FakeMapper)
    source = xr.DataArray(
        np.array([[[0.2]]], dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": ["B02"], "y": [0], "x": [0]},
    )

    mapped, unc = spectral_mapping_mod.map_multispectral_reflectance(
        source,
        source_bands=(SensorBand("B02", 490.0, 65.0, 10.0, 0),),
        target_bands=(SensorBand("T01", 500.0, 50.0, 10.0, 0),),
        source_uncertainty=None,
        spectral_library=_library(),
        k_neighbors=2,
    )

    assert (mapped, unc) == ("mapped", "unc")
    assert captured["init"]["k_neighbors"] == 2
    assert captured["map"]["source_reflectance"] is source


def test_curve_utils_cover_classification_and_error_paths() -> None:
    visible = SensorBand("B03", 560.0, 35.0, 10.0, 0)
    swir = SensorBand("B11", 1610.0, 90.0, 20.0, 1)
    assert curve_utils.classify_band_region(visible) == "visible"
    assert curve_utils.classify_band_region(swir) == "infrared"
    assert curve_utils.segment_for_band(swir) == "swir"
    assert curve_utils.primary_nir_band_index((visible, swir)) is None

    with pytest.raises(ValueError, match="zero support"):
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
