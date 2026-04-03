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
        map_reflectance_batch=lambda **_kwargs: SimpleNamespace(
            results=[
                SimpleNamespace(
                    target_reflectance=[0.4],
                    target_band_ids=["T01"],
                    diagnostics={"segments": {}},
                )
            ]
        )
    )
    nonidentity._target_internal_to_output_index = {"T01": 0}
    nonidentity._source_retrieval_indices_by_segment = {
        "vnir": np.array([0], dtype=np.int32),
        "swir": np.array([], dtype=np.int32),
    }
    nonidentity._target_schema_by_band_id = {}

    single = xr.DataArray(
        np.array([[[0.25]]], dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": ["B02"], "y": [0], "x": [0]},
    )
    mapped, mapped_unc, source_fit = nonidentity.map(single, source_uncertainty=unc.sel(band=["B02"]))
    assert float(mapped.values[0, 0, 0]) == pytest.approx(0.4)
    assert float(mapped_unc.values[0, 0, 0]) == pytest.approx(0.005)
    assert float(source_fit.values[0, 0]) == pytest.approx(0.0)

    uninitialized = object.__new__(SpectralMapper)
    uninitialized.source_bands = nonidentity.source_bands
    uninitialized.target_bands = nonidentity.target_bands
    uninitialized._identity = False
    uninitialized._package_mapper = None
    uninitialized._runtime = None
    with pytest.raises(RuntimeError, match="not initialized"):
        uninitialized.map(single)


def test_uncertainty_helpers_cover_normalization_loading_and_projection() -> None:
    mapper = object.__new__(SpectralMapper)
    mapper.target_bands = (SensorBand("T01", 500.0, 50.0, 10.0, 0),)
    mapper._runtime = SimpleNamespace(target_sensor_id="target")
    mapper._package_mapper = SimpleNamespace(
        _row_index_by_id={"n1": 1, "n2": 0},
        _load_hyperspectral=lambda _segment: np.array([[2.0, 4.0], [1.0, 3.0]], dtype=np.float64),
        _band_response=lambda *_args, **_kwargs: np.array([0.0, 0.0], dtype=np.float64),
    )
    mapper._target_internal_to_output_index = {"T01": 0}
    mapper._target_schema_by_band_id = {"T01": SimpleNamespace(band_id="T01", segment="vnir")}
    mapper._source_retrieval_indices_by_segment = {
        "vnir": np.array([0, 1], dtype=np.int32),
        "swir": np.array([], dtype=np.int32),
    }

    assert SpectralMapper._prepare_segment_diagnostics(None) is None
    assert SpectralMapper._prepare_segment_diagnostics({"status": "bad"}) is None
    assert SpectralMapper._prepare_segment_diagnostics({"status": "ok", "neighbor_ids": [], "neighbor_weights": [1.0]}) is None
    assert SpectralMapper._prepare_segment_diagnostics({"status": "ok", "neighbor_ids": ["n1"], "neighbor_weights": []}) is None
    assert SpectralMapper._prepare_segment_diagnostics({"status": "ok", "neighbor_ids": ["n1"], "neighbor_weights": [0.0]}) is None

    prepared = SpectralMapper._prepare_segment_diagnostics(
        {
            "status": "ok",
            "neighbor_ids": ["n1", "n2"],
            "neighbor_weights": [1.0, 3.0],
            "query_valid_mask": [True, False],
            "source_fit_rmse": 0.2,
        }
    )
    assert prepared is not None
    np.testing.assert_allclose(prepared.neighbor_weights, np.array([0.25, 0.75], dtype=np.float64))

    loaded = mapper._load_neighbor_spectra(("n1", "n2"), "vnir")
    np.testing.assert_allclose(loaded, np.array([[1.0, 3.0], [2.0, 4.0]], dtype=np.float64))
    assert mapper._segment_input_uncertainty("vnir", np.array([], dtype=bool), None) == 0.0
    assert mapper._segment_input_uncertainty("vnir", np.array([], dtype=bool), np.array([0.1, 0.2], dtype=np.float32)) == 0.0
    assert mapper._segment_input_uncertainty(
        "vnir",
        np.array([True, False], dtype=bool),
        np.array([0.3, 0.4], dtype=np.float32),
    ) == pytest.approx(0.3)
    assert mapper._neighbor_target_projection(loaded, SimpleNamespace(segment="vnir")) is None

    mapper._package_mapper = SimpleNamespace(
        _row_index_by_id={"n1": 0, "n2": 1},
        _load_hyperspectral=lambda _segment: np.array([[1.0, 3.0], [3.0, 5.0]], dtype=np.float64),
        _band_response=lambda *_args, **_kwargs: np.array([1.0, 1.0], dtype=np.float64),
    )
    projection = mapper._neighbor_target_projection(
        np.array([[1.0, 3.0], [3.0, 5.0]], dtype=np.float64),
        SimpleNamespace(segment="vnir"),
    )
    np.testing.assert_allclose(projection, np.array([2.0, 4.0], dtype=np.float64))

    result = SimpleNamespace(
        target_reflectance=None,
        target_band_ids=[],
        diagnostics={
            "segments": {
                "vnir": {
                    "status": "ok",
                    "neighbor_ids": ["n1", "n2"],
                    "neighbor_weights": [1.0, 1.0],
                    "query_valid_mask": [True, True],
                    "source_fit_rmse": 0.1,
                }
            }
        },
    )
    output = mapper._estimate_uncertainty(
        result,
        source_uncertainty=np.array([0.1, 0.3], dtype=np.float32),
    )
    assert output.shape == (1,)
    assert float(output[0]) > 0.1


def test_map_prefers_minimal_batch_path_when_package_private_hooks_exist() -> None:
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
    mapper._source_retrieval_indices_by_segment = {
        "vnir": np.array([0], dtype=np.int32),
        "swir": np.array([], dtype=np.int32),
    }
    mapper._target_schema_by_band_id = {"T01": SimpleNamespace(band_id="T01", segment="vnir")}

    calls = {"private": 0, "public": 0}

    retrieval_ok = SimpleNamespace(
        success=True,
        reconstructed=np.array([0.4, 0.6], dtype=np.float64),
        neighbor_ids=("n1",),
        neighbor_weights=np.array([1.0], dtype=np.float64),
        query_valid_mask=np.array([True], dtype=np.bool_),
        source_fit_rmse=0.1,
    )
    retrieval_empty = SimpleNamespace(
        success=False,
        reconstructed=None,
        neighbor_ids=(),
        neighbor_weights=np.array([], dtype=np.float64),
        query_valid_mask=np.array([], dtype=np.bool_),
        source_fit_rmse=0.0,
    )

    def _public_map_reflectance_batch(**_kwargs):  # noqa: ANN003
        calls["public"] += 1
        raise AssertionError("public batch path should not be used when private hooks are available")

    def _retrieve_segment_batch(*, segment, **_kwargs):  # noqa: ANN003
        calls["private"] += 1
        if segment == "vnir":
            return (retrieval_ok,)
        return (retrieval_empty,)

    mapper._package_mapper = SimpleNamespace(
        map_reflectance_batch=_public_map_reflectance_batch,
        _candidate_rows=lambda _candidate_rows: np.array([0], dtype=np.int64),
        _retrieve_segment_batch=_retrieve_segment_batch,
        _row_index_by_id={"n1": 0},
        _load_source_matrix=lambda _sensor, _segment: np.array([[0.25]], dtype=np.float64),
        _load_hyperspectral=lambda _segment: np.array([[0.4, 0.6]], dtype=np.float64),
        _band_response=lambda *_args, **_kwargs: np.array([1.0, 1.0], dtype=np.float64),
    )

    source = xr.DataArray(
        np.array([[[0.25]]], dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": ["B02"], "y": [0], "x": [0]},
    )
    source_unc = xr.full_like(source, 0.02)

    mapped, mapped_unc, source_fit = mapper.map(source, source_uncertainty=source_unc)

    assert calls["public"] == 0
    # Fast path bypasses _retrieve_segment_batch; uses _load_source_matrix directly.
    assert np.isfinite(mapped.values).any()
    assert float(mapped_unc.values[0, 0, 0]) >= 0.0


def test_map_returns_source_fit_rmse() -> None:
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
    mapper._source_retrieval_indices_by_segment = {
        "vnir": np.array([0], dtype=np.int32),
        "swir": np.array([], dtype=np.int32),
    }
    mapper._target_schema_by_band_id = {"T01": SimpleNamespace(band_id="T01", segment="vnir")}

    retrieval_ok = SimpleNamespace(
        success=True,
        reconstructed=np.array([0.4, 0.6], dtype=np.float64),
        neighbor_ids=("n1",),
        neighbor_weights=np.array([1.0], dtype=np.float64),
        query_valid_mask=np.array([True], dtype=np.bool_),
        source_fit_rmse=0.1,
    )
    retrieval_empty = SimpleNamespace(
        success=False,
        reconstructed=None,
        neighbor_ids=(),
        neighbor_weights=np.array([], dtype=np.float64),
        query_valid_mask=np.array([], dtype=np.bool_),
        source_fit_rmse=0.0,
    )

    mapper._package_mapper = SimpleNamespace(
        _candidate_rows=lambda _candidate_rows: np.array([0], dtype=np.int64),
        _retrieve_segment_batch=lambda *, segment, **_kwargs: ((retrieval_ok,) if segment == "vnir" else (retrieval_empty,)),
        _row_index_by_id={"n1": 0},
        _load_source_matrix=lambda _sensor, _segment: np.array([[0.25]], dtype=np.float64),
        _load_hyperspectral=lambda _segment: np.array([[0.4, 0.6]], dtype=np.float64),
        _band_response=lambda *_args, **_kwargs: np.array([1.0, 1.0], dtype=np.float64),
    )

    source = xr.DataArray(
        np.array([[[0.25]]], dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": ["B02"], "y": [0], "x": [0]},
    )

    mapped, mapped_unc, source_fit = mapper.map(source)

    assert np.isfinite(mapped.values[0, 0, 0])
    assert float(mapped_unc.values[0, 0, 0]) >= 0.0
    assert np.isfinite(source_fit.values[0, 0])


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
    mapper._source_retrieval_indices_by_segment = {
        "vnir": np.array([0], dtype=np.int32),
        "swir": np.array([], dtype=np.int32),
    }
    mapper._target_schema_by_band_id = {"T01": SimpleNamespace(band_id="T01", segment="vnir")}

    retrieval_ok = SimpleNamespace(
        success=True,
        reconstructed=np.array([0.4, 0.6], dtype=np.float64),
        neighbor_ids=("n1",),
        neighbor_weights=np.array([1.0], dtype=np.float64),
        query_valid_mask=np.array([True], dtype=np.bool_),
        source_fit_rmse=0.1,
    )
    retrieval_empty = SimpleNamespace(
        success=False,
        reconstructed=None,
        neighbor_ids=(),
        neighbor_weights=np.array([], dtype=np.float64),
        query_valid_mask=np.array([], dtype=np.bool_),
        source_fit_rmse=0.0,
    )

    mapper._package_mapper = SimpleNamespace(
        _candidate_rows=lambda _candidate_rows: np.array([0], dtype=np.int64),
        _retrieve_segment_batch=lambda *, segment, **_kwargs: ((retrieval_ok,) if segment == "vnir" else (retrieval_empty,)),
        _row_index_by_id={"n1": 0},
        _load_source_matrix=lambda _sensor, _segment: np.array([[0.25]], dtype=np.float64),
        _load_hyperspectral=lambda _segment: np.array([[0.4, 0.6]], dtype=np.float64),
        _band_response=lambda *_args, **_kwargs: np.array([1.0, 1.0], dtype=np.float64),
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
    assert len(data_paths) == 3
    assert len(metadata_paths) == 1

    for path in data_paths:
        with rasterio.open(path) as src:
            assert src.crs.to_string() == "EPSG:32632"
            assert src.transform == source.isel(band=0, drop=True).rio.transform()

    metadata = json.loads(metadata_paths[0].read_text(encoding="utf-8"))
    assert metadata["source_sensor_id"] == "src"
    assert metadata["target_sensor_id"] == "target"
    assert sorted(metadata["metrics"]) == sorted(["source_fit_rmse", "vnir_source_fit_rmse", "swir_source_fit_rmse"])


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
