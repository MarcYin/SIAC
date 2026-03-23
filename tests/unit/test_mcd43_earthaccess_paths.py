"""Focused coverage lifts for Earthaccess-backed BRDF providers."""

from __future__ import annotations

from datetime import date, datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

import siac.adapters.brdf.mcd43_earthaccess as mcd_mod
from siac.adapters.brdf.mcd43_earthaccess import (
    MCD19EarthAccessProvider,
    MCD43EarthAccessProvider,
)
from siac.adapters.earthdata_common import ProductBandDefinition
from siac.domain import SensorBand


class _StubSource:
    def __init__(self) -> None:
        self.search_calls: list[dict[str, object]] = []
        self.download_calls: list[tuple[list[object], Path]] = []

    def search_granules(self, **kwargs: object) -> list[object]:
        self.search_calls.append(kwargs)
        return [{"id": "g1"}]

    def download_granules(self, granules: list[object], dest_dir: Path) -> list[Path]:
        self.download_calls.append((list(granules), dest_dir))
        return [dest_dir / "granule_1.hdf"]


class _StubCatalog:
    def __init__(self) -> None:
        self.keys: list[str] = []

    def resolve_short_name(self, key: str) -> str:
        self.keys.append(key)
        return f"SN:{key}"


class _DummyProvider(mcd_mod._EarthAccessBRDFProvider):
    product_key = "dummy_brdf"
    _source_name = "DUMMY"
    _product_bands = (
        ProductBandDefinition("Band3", 469.0, 20.0, "params_3", "qa_3"),
        ProductBandDefinition("Band4", 555.0, 20.0, "params_4", "qa_4"),
    )
    _legacy_band_map = {1: "Band3", 2: "Band4"}

    def _load_native_band_stack(
        self,
        path: str | Path,
        product_band: ProductBandDefinition,
    ) -> tuple[tuple[xr.DataArray, xr.DataArray, xr.DataArray], xr.DataArray]:
        raise AssertionError(f"unexpected call for {path} / {product_band.label}")


def _install_runtime(
    monkeypatch: pytest.MonkeyPatch,
    *,
    source: _StubSource | None = None,
    catalog: _StubCatalog | None = None,
) -> tuple[_StubSource, _StubCatalog]:
    src = source or _StubSource()
    cat = catalog or _StubCatalog()
    monkeypatch.setattr(
        mcd_mod,
        "build_earthaccess_runtime",
        lambda **_kwargs: SimpleNamespace(
            cache_dir=Path("/tmp/mcd43-cache"),
            source=src,
            catalog=cat,
        ),
    )
    return src, cat


def _da(
    value: float | np.ndarray,
    *,
    dims: tuple[str, ...] = ("y", "x"),
    coords: dict[str, object] | None = None,
) -> xr.DataArray:
    values = np.asarray(value, dtype=np.float32)
    if coords is None and dims == ("y", "x"):
        coords = {"y": np.arange(values.shape[0]), "x": np.arange(values.shape[1])}
    return xr.DataArray(values, dims=dims, coords=coords)


def _fake_native_grid(
    values: np.ndarray,
    *,
    granule_path: str | Path,
    dims: tuple[str, ...] | None = None,
    coords: dict[str, object] | None = None,
) -> xr.DataArray:
    del granule_path
    arr = np.asarray(values, dtype=np.float32)
    if dims is None:
        dims = ("y", "x")
    data_coords = dict(coords or {})
    if "y" in dims and "y" not in data_coords:
        data_coords["y"] = np.arange(arr.shape[dims.index("y")], dtype=np.float32)
    if "x" in dims and "x" not in data_coords:
        data_coords["x"] = np.arange(arr.shape[dims.index("x")], dtype=np.float32)
    return xr.DataArray(arr, dims=dims, coords=data_coords)


def test_granule_helpers_cover_render_dict_parse_fallbacks_and_repr() -> None:
    render_obj = SimpleNamespace(
        render_dict={
            "umm": {
                "TemporalExtent": {
                    "RangeDateTime": {
                        "BeginningDateTime": "2024-01-02T00:00:00.000Z",
                        "EndingDateTime": "2024-01-02T23:59:59.999Z",
                    }
                },
                "GranuleUR": "unused",
            }
        }
    )
    assert mcd_mod._granule_day(render_obj) == np.datetime64("2024-01-02", "D")

    parsed_from_name = {
        "umm": {
            "TemporalExtent": {"SingleDateTime": "not-a-date"},
            "ProducerGranuleId": "MCD43A1.A2024010.h29v07.061.fake.hdf",
        }
    }
    assert mcd_mod._granule_day(parsed_from_name) == np.datetime64("2024-01-10", "D")
    assert mcd_mod._granule_day(object()) is None

    assert mcd_mod._granule_key({"meta": {"native-id": "native-1"}}) == "native-1"
    assert mcd_mod._granule_key({"meta": {"concept-id": "concept-1"}}) == "concept-1"
    assert mcd_mod._granule_key({"umm": {"GranuleUR": "", "ProducerGranuleId": "producer-1"}}) == "producer-1"
    fallback = object()
    assert mcd_mod._granule_key(fallback) == repr(fallback)


def test_provider_validation_band_resolution_and_time_helpers(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_runtime(monkeypatch)
    provider = _DummyProvider(probe_earthdata=False)

    resolved = provider._resolve_requested_bands(
        [
            1,
            SensorBand("B02", 490.0, 65.0, 10.0, 1),
        ]
    )
    assert resolved[0][0] == 1
    assert resolved[0][1].label == "Band3"
    assert resolved[1][0] == "B02"
    assert resolved[1][1].label == "Band3"

    with pytest.raises(ValueError, match="non-empty sequence"):
        provider._resolve_requested_bands([])
    with pytest.raises(KeyError, match="not available"):
        provider._resolve_requested_bands([99])

    assert provider._resolved_short_name() == "SN:dummy_brdf"
    assert mcd_mod._EarthAccessBRDFProvider._coerce_sample_day(np.datetime64("2024-01-02T12:00:00")) == np.datetime64(
        "2024-01-02",
        "D",
    )
    assert mcd_mod._EarthAccessBRDFProvider._coerce_sample_day(datetime(2024, 1, 3, 4, 5, 6)) == np.datetime64(
        "2024-01-03",
        "D",
    )
    assert list(mcd_mod._EarthAccessBRDFProvider._time_axis(datetime(2024, 1, 2, 12, 0, 0), 1)) == [
        np.datetime64("2024-01-01"),
        np.datetime64("2024-01-02"),
        np.datetime64("2024-01-03"),
    ]
    with pytest.raises(ValueError, match="must not be empty"):
        mcd_mod._EarthAccessBRDFProvider._coerce_sample_time_axis([])

    sample_axis = mcd_mod._EarthAccessBRDFProvider._coerce_sample_time_axis(
        [date(2024, 1, 5), datetime(2024, 1, 2, 10, 0, 0), np.datetime64("2024-01-05T12:00:00")]
    )
    assert list(sample_axis) == [np.datetime64("2024-01-02"), np.datetime64("2024-01-05")]

    assert mcd_mod._EarthAccessBRDFProvider._temporal_search_window(
        datetime(2024, 1, 3, 12, 0, 0),
        2,
        np.array(["2024-01-01", "2024-01-08"], dtype="datetime64[D]"),
    ) == ("2024-01-01T00:00:00Z", "2024-01-08T23:59:59Z")

    for fn in (provider.get_brdf_parameters, provider.get_temporal_brdf_parameters):
        with pytest.raises(ValueError, match="must be > 0"):
            fn(
                bounds=(0.0, 0.0, 1.0, 1.0),
                crs="EPSG:4326",
                obs_time=datetime(2024, 1, 1, 12, 0, 0),
                target_resolution=0.0,
                bands=[1],
            )
        with pytest.raises(ValueError, match="non-empty sequence"):
            fn(
                bounds=(0.0, 0.0, 1.0, 1.0),
                crs="EPSG:4326",
                obs_time=datetime(2024, 1, 1, 12, 0, 0),
                target_resolution=500.0,
                bands=[],
            )

    with pytest.raises(ValueError, match="must be > 0"):
        provider.get_temporal_brdf_parameters_batch(
            bounds=(0.0, 0.0, 1.0, 1.0),
            crs="EPSG:4326",
            obs_times=[datetime(2024, 1, 1, 12, 0, 0)],
            target_resolution=0.0,
            bands=[1],
            temporal_windows=[1],
            sample_date_sets=[[datetime(2024, 1, 1)]],
        )
    with pytest.raises(ValueError, match="non-empty sequence"):
        provider.get_temporal_brdf_parameters_batch(
            bounds=(0.0, 0.0, 1.0, 1.0),
            crs="EPSG:4326",
            obs_times=[datetime(2024, 1, 1, 12, 0, 0)],
            target_resolution=500.0,
            bands=[],
            temporal_windows=[1],
            sample_date_sets=[[datetime(2024, 1, 1)]],
        )
    with pytest.raises(ValueError, match="same length"):
        provider.get_temporal_brdf_parameters_batch(
            bounds=(0.0, 0.0, 1.0, 1.0),
            crs="EPSG:4326",
            obs_times=[datetime(2024, 1, 1, 12, 0, 0)],
            target_resolution=500.0,
            bands=[1],
            temporal_windows=[1, 2],
            sample_date_sets=[[datetime(2024, 1, 1)]],
        )


def test_download_and_filter_helpers_cover_warning_and_passthrough_branches(
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    source, catalog = _install_runtime(monkeypatch)
    provider = _DummyProvider(source=source, catalog=catalog, probe_earthdata=True)

    assert provider._search_granules(
        short_name="SN:dummy_brdf",
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        temporal=("2024-01-01T00:00:00Z", "2024-01-02T23:59:59Z"),
    ) == [{"id": "g1"}]
    assert source.search_calls[-1]["count"] == provider.max_granules

    monkeypatch.setattr(provider, "_search_granules", lambda **_kwargs: [])
    with caplog.at_level("WARNING"):
        assert provider._download_granules(
            bounds=(0.0, 0.0, 1.0, 1.0),
            crs="EPSG:4326",
            obs_time=datetime(2024, 1, 1, 12, 0, 0),
            temporal_window=1,
        ) == []
    assert "returned no results" in caplog.text

    granules = [{"meta": {"native-id": "MCD43A1.A2024001.h29v07.061.fake.hdf"}}]
    sample_dates = np.array(["2024-01-08"], dtype="datetime64[D]")
    monkeypatch.setattr(provider, "_filter_granules_to_sample_dates", lambda _g, _d: [])
    with caplog.at_level("WARNING"):
        assert provider._filter_sampled_granules(
            granules,
            sample_dates,
            empty_message="%s empty after filtering",
        ) == []
    assert "empty after filtering" in caplog.text
    assert provider._filter_sampled_granules(granules, None) == granules

    selected_calls: dict[str, object] = {}
    monkeypatch.setattr(
        mcd_mod,
        "select_candidate_paths",
        lambda paths, **kwargs: selected_calls.update({"paths": paths, **kwargs}) or paths[:1],
    )
    selected = provider._select_candidate_paths(
        [Path("/tmp/a.hdf"), Path("/tmp/b.hdf")],
        datetime(2024, 1, 1, 12, 0, 0),
        (0.0, 0.0, 1.0, 1.0),
        "EPSG:4326",
        sample_dates=sample_dates,
    )
    assert selected == [Path("/tmp/a.hdf")]
    assert selected_calls["crs"] == "EPSG:4326"

    cached = provider._download_granules_to_cache([{"id": "g1"}], short_name="SN:dummy_brdf")
    assert cached == [source.download_calls[-1][1] / "granule_1.hdf"]
    assert source.download_calls[-1][1] == Path("/tmp/mcd43-cache")

    monkeypatch.setattr(
        provider,
        "_merge_search_batches",
        lambda *_args, **_kwargs: [
            (
                np.datetime64("2024-01-01"),
                np.datetime64("2024-01-02"),
                np.array(["2024-01-01"], dtype="datetime64[D]"),
            )
        ],
    )
    monkeypatch.setattr(provider, "_search_granules", lambda **_kwargs: [])
    with caplog.at_level("WARNING"):
        assert provider._download_granules_batch(
            bounds=(0.0, 0.0, 1.0, 1.0),
            crs="EPSG:4326",
            request_specs=[(datetime(2024, 1, 1, 12, 0, 0), sample_dates)],
            temporal_windows=[1],
        ) == []
    assert "no sampled results" in caplog.text


def test_load_from_granules_fills_missing_values_and_builds_uncertainties(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_runtime(monkeypatch)
    provider = _DummyProvider(probe_earthdata=False)

    params = (
        _da([[np.nan, 0.4], [0.3, np.nan]]),
        _da([[np.nan, 0.2], [0.1, np.nan]]),
        _da([[np.nan, 0.06], [0.04, np.nan]]),
    )
    qa = _da([[0.0, 1.0], [np.nan, -1.0]])
    monkeypatch.setattr(provider, "_load_native_band_stack", lambda _path, _band: (params, qa))
    monkeypatch.setattr(provider, "_merge_reprojected_tiles", lambda arrays, **_kwargs: arrays[0])

    weights = provider._load_from_granules(
        [Path("/tmp/day1.hdf")],
        requested=[("B02", provider._product_bands[0])],
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        target_resolution=500.0,
    )

    assert weights.f0.dims == ("band", "y", "x")
    np.testing.assert_allclose(weights.f0.sel(band="B02").values[0], np.array([0.20, 0.4], dtype=np.float32))
    np.testing.assert_allclose(weights.f1.sel(band="B02").values[0], np.array([0.05, 0.2], dtype=np.float32))
    np.testing.assert_allclose(weights.f2.sel(band="B02").values[0], np.array([0.02, 0.06], dtype=np.float32))
    assert float(weights.reflectance_unc.sel(band="B02").values[0, 0]) == pytest.approx(0.015)
    assert float(weights.reflectance_unc.sel(band="B02").values[1, 0]) == pytest.approx(0.08)


def test_temporal_loading_covers_missing_days_target_grid_coercion_and_fallback(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_runtime(monkeypatch)
    provider = _DummyProvider(probe_earthdata=False)
    time_axis = np.array(["2024-01-01", "2024-01-02", "2024-01-03"], dtype="datetime64[D]")
    day_map = {
        Path("/tmp/day1.hdf"): datetime(2024, 1, 1),
        Path("/tmp/day3.hdf"): datetime(2024, 1, 3),
    }

    monkeypatch.setattr(mcd_mod, "parse_granule_date", lambda path: day_map[Path(path)])
    monkeypatch.setattr(provider, "_grid", lambda _bounds, _resolution: (np.array([10.0, 20.0]), np.array([30.0, 40.0])))

    def _stack_for_path(path: Path, _band: ProductBandDefinition):
        base = 0.1 if "day1" in path.name else 0.3
        coords = {
            "y": [0.0, 1.0],
            "x": [0.0, 1.0],
            "extra": xr.DataArray([1.0, 2.0], dims=("y",)),
        }
        params = (
            xr.DataArray(np.full((2, 2), base, dtype=np.float32), dims=("y", "x"), coords=coords),
            xr.DataArray(np.full((2, 2), base / 2.0, dtype=np.float32), dims=("y", "x"), coords=coords),
            xr.DataArray(np.full((2, 2), base / 4.0, dtype=np.float32), dims=("y", "x"), coords=coords),
        )
        qa = xr.DataArray(np.zeros((2, 2), dtype=np.float32), dims=("y", "x"), coords=coords)
        return params, qa

    monkeypatch.setattr(provider, "_load_native_band_stack", _stack_for_path)
    monkeypatch.setattr(provider, "_merge_reprojected_tiles", lambda arrays, **_kwargs: arrays[0])

    temporal = provider._load_temporal_from_granules(
        [Path("/tmp/day1.hdf"), Path("/tmp/day3.hdf")],
        requested=[("B02", provider._product_bands[0])],
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        target_resolution=500.0,
        time_axis=time_axis,
    )

    assert temporal.f0.dims == ("time", "band", "y", "x")
    assert np.isnan(temporal.f0.sel(time=np.datetime64("2024-01-02")).values).all()
    assert float(temporal.f0.sel(time=np.datetime64("2024-01-01")).mean()) == pytest.approx(0.1)
    assert float(temporal.f0.sel(time=np.datetime64("2024-01-03")).mean()) == pytest.approx(0.3)
    assert list(temporal.f0.coords["y"].values) == [10.0, 20.0]
    assert "extra" not in temporal.f0.coords

    sentinel = provider._default_temporal_weights((0.0, 0.0, 1.0, 1.0), 500.0, ["B02"], time_axis)
    monkeypatch.setattr(
        provider,
        "_load_native_band_stack",
        lambda _path, _band: (
            (
                _da(np.full((2, 2), np.nan)),
                _da(np.full((2, 2), np.nan)),
                _da(np.full((2, 2), np.nan)),
            ),
            _da(np.full((2, 2), np.nan)),
        ),
    )
    monkeypatch.setattr(provider, "_default_temporal_weights", lambda *_args, **_kwargs: sentinel)
    fallback = provider._load_temporal_from_granules(
        [Path("/tmp/day1.hdf")],
        requested=[("B02", provider._product_bands[0])],
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        target_resolution=500.0,
        time_axis=time_axis[:1],
    )
    assert fallback is sentinel


def test_stack_parameter_provider_validates_shape_and_splits_parameters(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_runtime(monkeypatch)
    provider = MCD43EarthAccessProvider(probe_earthdata=False)
    product_band = provider._product_bands[0]

    invalid_calls: list[str] = []

    def _invalid_read_dataset(_path: str | Path, dataset_name: str):
        invalid_calls.append(dataset_name)
        if dataset_name == product_band.parameter_dataset:
            return np.ones((2, 2), dtype=np.int16), {}
        return np.zeros((2, 2), dtype=np.int16), {}

    monkeypatch.setattr(MCD43EarthAccessProvider, "_read_dataset", staticmethod(_invalid_read_dataset))
    monkeypatch.setattr(mcd_mod, "apply_scale_and_mask", lambda values, _attrs: np.asarray(values, dtype=np.float32))
    monkeypatch.setattr(mcd_mod, "make_native_grid_dataarray", _fake_native_grid)

    with pytest.raises(ValueError, match="3-parameter BRDF stack"):
        provider._load_native_band_stack(Path("/tmp/granule.hdf"), product_band)
    assert invalid_calls == [product_band.parameter_dataset, product_band.qa_dataset]

    def _valid_read_dataset(_path: str | Path, dataset_name: str):
        if dataset_name == product_band.parameter_dataset:
            return np.dstack(
                [
                    np.full((2, 2), 10, dtype=np.int16),
                    np.full((2, 2), 20, dtype=np.int16),
                    np.full((2, 2), 30, dtype=np.int16),
                ]
            ), {}
        return np.full((2, 2), 4, dtype=np.int16), {}

    monkeypatch.setattr(MCD43EarthAccessProvider, "_read_dataset", staticmethod(_valid_read_dataset))
    (f0, f1, f2), qa = provider._load_native_band_stack(Path("/tmp/granule.hdf"), product_band)
    assert float(f0.mean()) == pytest.approx(10.0)
    assert float(f1.mean()) == pytest.approx(20.0)
    assert float(f2.mean()) == pytest.approx(30.0)
    assert float(qa.mean()) == pytest.approx(4.0)


def test_mcd19_native_stack_reads_expected_fields(monkeypatch: pytest.MonkeyPatch) -> None:
    _install_runtime(monkeypatch)
    provider = MCD19EarthAccessProvider(probe_earthdata=False)
    product_band = provider._product_bands[2]
    seen: list[str] = []

    def _fake_read_hdf4_dataset(_path: str | Path, dataset_name: str):
        seen.append(dataset_name)
        values = {
            "Kiso_Band3": np.full((2, 2), 10, dtype=np.int16),
            "Kvol_Band3": np.full((2, 2), 20, dtype=np.int16),
            "Kgeo_Band3": np.full((2, 2), 30, dtype=np.int16),
            "Status_QA": np.full((2, 2), 1, dtype=np.int16),
        }[dataset_name]
        return values, {}

    monkeypatch.setattr(mcd_mod, "read_hdf4_dataset", _fake_read_hdf4_dataset)
    monkeypatch.setattr(mcd_mod, "apply_scale_and_mask", lambda values, _attrs: np.asarray(values, dtype=np.float32))
    monkeypatch.setattr(mcd_mod, "make_native_grid_dataarray", _fake_native_grid)

    (f0, f1, f2), qa = provider._load_native_band_stack(Path("/tmp/granule.hdf"), product_band)
    assert seen == ["Kiso_Band3", "Kvol_Band3", "Kgeo_Band3", "Status_QA"]
    assert float(f0.mean()) == pytest.approx(10.0)
    assert float(f1.mean()) == pytest.approx(20.0)
    assert float(f2.mean()) == pytest.approx(30.0)
    assert float(qa.mean()) == pytest.approx(1.0)
