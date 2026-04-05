from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import rasterio
import xarray as xr

import siac.algorithms.surface.spectral_mapping as spectral_mapping_mod
from siac.algorithms.surface import _spectral_curve_utils as curve_utils
from siac.algorithms.surface.spectral_mapping import HyperspectralLibrary, SpectralMapper
from siac.domain import SensorBand


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


def _with_geo(data: xr.DataArray) -> xr.DataArray:
    height = int(data.sizes["y"])
    width = int(data.sizes["x"])
    resolution = 10.0
    xmin = 500000.0
    ymax = 4100000.0
    x = np.linspace(xmin + resolution / 2.0, xmin + width * resolution - resolution / 2.0, width, dtype=np.float64)
    y = np.linspace(ymax - resolution / 2.0, ymax - height * resolution + resolution / 2.0, height, dtype=np.float64)
    transform = rasterio.transform.from_origin(xmin, ymax, resolution, resolution)
    out = data.assign_coords({"x": x, "y": y}).rio.set_spatial_dims(x_dim="x", y_dim="y")
    return out.rio.write_crs("EPSG:32632").rio.write_transform(transform)


def test_runtime_export_helpers_cover_hashing_and_cache_reuse(tmp_path: Path) -> None:
    shared_rsrf = (
        SensorBand("B02", 490.0, 65.0, 10.0, 0, rsrf_sensor_unit_id="sentinel-2a_msi"),
        SensorBand("B03", 560.0, 35.0, 10.0, 1, rsrf_sensor_unit_id="sentinel-2a_msi"),
    )
    assert spectral_mapping_mod._sensor_id_for_bands(shared_rsrf, prefix="source") == "sentinel-2a_msi"

    payload, schema_bands = spectral_mapping_mod._schema_payload_for_bands(
        "source",
        _source_bands(),
        rsrf_root=None,
    )
    assert payload["sensor_id"] == "source"
    assert len(payload["bands"]) == len(_source_bands())
    assert schema_bands[1].schema_band_id == "nir"
    assert all(band["support_min_nm"] <= band["support_max_nm"] for band in payload["bands"])

    export_root = tmp_path / "export"
    spectral_mapping_mod._export_hyperspectral_library_root(export_root, _library())
    assert (export_root / "tabular" / "siac_spectra_metadata.csv").exists()
    assert (export_root / "tabular" / "siac_normalized_spectra.csv").exists()

    config = spectral_mapping_mod.SpectralMappingConfig(
        cache_dir=tmp_path / "cache",
        siac_library_root=None,
    ).normalized()
    ensured = spectral_mapping_mod._ensure_siac_library_root(
        Path(config.cache_dir),
        "sig",
        _library(),
        config,
    )
    reused = spectral_mapping_mod._ensure_siac_library_root(
        Path(config.cache_dir),
        "sig",
        _library(),
        config,
    )
    assert ensured == reused

    explicit = spectral_mapping_mod.SpectralMappingConfig(siac_library_root=tmp_path / "explicit").normalized()
    assert spectral_mapping_mod._ensure_siac_library_root(tmp_path, "sig", None, explicit) == explicit.siac_library_root


def test_export_hyperspectral_library_root_rejects_missing_canonical_coverage(tmp_path: Path) -> None:
    bad = HyperspectralLibrary(
        wavelengths_nm=np.arange(450.0, 2501.0, 10.0, dtype=np.float32),
        spectra=np.stack(
            [np.linspace(0.1, 0.2, 206, dtype=np.float32), np.linspace(0.2, 0.3, 206, dtype=np.float32)],
            axis=0,
        ),
        sample_ids=("a", "b"),
    )

    with pytest.raises(ValueError, match="canonical 400-2500 nm range"):
        spectral_mapping_mod._export_hyperspectral_library_root(tmp_path / "bad", bad)


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
    mapped, mapped_unc, source_fit = nonidentity.map(single, source_uncertainty=unc.sel(band=["B02"]))
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


def test_wrapper_and_input_split_helpers_cover_error_and_delegation_paths(monkeypatch: pytest.MonkeyPatch) -> None:
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
