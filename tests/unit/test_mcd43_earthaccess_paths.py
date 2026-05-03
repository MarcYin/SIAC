"""Focused coverage lifts for Earthaccess-backed BRDF providers."""

from __future__ import annotations

import logging
import sys
from datetime import date, datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr
from affine import Affine

import siac.adapters.brdf.mcd43_earthaccess as mcd_mod
import siac.adapters.earthdata as earth_adapter
from siac.adapters.brdf.mcd43_earthaccess import (
    MCD19EarthAccessProvider,
    MCD43EarthAccessProvider,
)
from siac.adapters.earthdata_common import ProductBandDefinition
from siac.domain import SensorBand
from siac.runtime import BRDFKernelWeights


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
    monkeypatch.setattr(
        mcd_mod,
        "build_source_aligned_target_template",
        lambda *_args, **kwargs: (
            xr.DataArray(
                np.zeros((1, 1), dtype=np.float32),
                dims=("y", "x"),
                coords={"y": [0.0], "x": [0.0]},
            )
            .rio.set_spatial_dims(x_dim="x", y_dim="y")
            .rio.write_crs("EPSG:32615"),
            float(kwargs.get("resolution") or 500.0),
        ),
    )
    monkeypatch.setattr(
        mcd_mod,
        "resolve_gdal_subdataset_path",
        lambda path, dataset_name: f"{Path(path).name}:{dataset_name}",
    )
    return src, cat


def _resampling_name(value: object) -> object:
    return getattr(value, "name", value)


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
    filename_wins = {
        "umm": {
            "TemporalExtent": {
                "RangeDateTime": {
                    "BeginningDateTime": "2023-07-01T00:00:00.000Z",
                    "EndingDateTime": "2023-07-16T23:59:59.999Z",
                }
            },
            "ProducerGranuleId": "MCD43A1.A2023197.h10v04.061.fake.hdf",
        }
    }
    assert mcd_mod._granule_day(filename_wins) == np.datetime64("2023-07-16", "D")

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
    assert (
        mcd_mod._granule_key({"umm": {"GranuleUR": "", "ProducerGranuleId": "producer-1"}})
        == "producer-1"
    )
    fallback = object()
    assert mcd_mod._granule_key(fallback) == repr(fallback)


def test_provider_validation_band_resolution_and_time_helpers(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_runtime(monkeypatch)
    provider = _DummyProvider(probe_earthdata=False)

    resolved = provider._resolve_requested_bands([SensorBand("B02", 490.0, 65.0, 10.0, 1)])
    assert resolved[0][0] == "B02"
    assert resolved[0][1].label == "Band3"

    with pytest.raises(ValueError, match="non-empty sequence"):
        provider._resolve_requested_bands([])
    with pytest.raises(TypeError, match="SensorBand"):
        provider._resolve_requested_bands([1])

    assert provider._resolved_short_name() == "SN:dummy_brdf"
    assert mcd_mod._EarthAccessBRDFProvider._coerce_sample_day(
        np.datetime64("2024-01-02T12:00:00")
    ) == np.datetime64(
        "2024-01-02",
        "D",
    )
    assert mcd_mod._EarthAccessBRDFProvider._coerce_sample_day(
        datetime(2024, 1, 3, 4, 5, 6)
    ) == np.datetime64(
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
                bands=[provider.source_bands[0]],
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
            bands=[provider.source_bands[0]],
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
            bands=[provider.source_bands[0]],
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
        assert (
            provider._download_granules(
                bounds=(0.0, 0.0, 1.0, 1.0),
                crs="EPSG:4326",
                obs_time=datetime(2024, 1, 1, 12, 0, 0),
                temporal_window=1,
            )
            == []
        )
    assert "returned no results" in caplog.text

    granules = [{"meta": {"native-id": "MCD43A1.A2024001.h29v07.061.fake.hdf"}}]
    sample_dates = np.array(["2024-01-08"], dtype="datetime64[D]")
    monkeypatch.setattr(provider, "_filter_granules_to_sample_dates", lambda _g, _d: [])
    with caplog.at_level("WARNING"):
        assert (
            provider._filter_sampled_granules(
                granules,
                sample_dates,
                empty_message="%s empty after filtering",
            )
            == []
        )
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
    expected_cache_dir = source.download_calls[-1][1]
    assert cached == [expected_cache_dir / "granule_1.hdf"]
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
        assert (
            provider._download_granules_batch(
                bounds=(0.0, 0.0, 1.0, 1.0),
                crs="EPSG:4326",
                request_specs=[(datetime(2024, 1, 1, 12, 0, 0), sample_dates)],
                temporal_windows=[1],
            )
            == []
        )
    assert "no sampled results" in caplog.text


def test_temporal_batch_loader_reads_combined_stack_once_and_slices_per_request(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source, catalog = _install_runtime(monkeypatch)
    provider = _DummyProvider(source=source, catalog=catalog, probe_earthdata=True)
    downloaded = [Path("/tmp/january.hdf"), Path("/tmp/february.hdf")]

    monkeypatch.setattr(
        provider,
        "_download_granules_batch",
        lambda **_kwargs: downloaded,
    )

    def _select_candidate_paths(paths, obs_time, bounds, crs, sample_dates=None):  # noqa: ANN001
        del obs_time, bounds, crs
        first_day = str(np.asarray(sample_dates, dtype="datetime64[D]")[0])
        if first_day.startswith("2024-01"):
            return [paths[0]]
        return [paths[1]]

    monkeypatch.setattr(provider, "_select_candidate_paths", _select_candidate_paths)

    load_calls: list[dict[str, object]] = []

    def _load_temporal_from_granules(  # noqa: ANN001
        paths,
        *,
        requested,
        bounds,
        crs,
        target_resolution,
        time_axis,
    ):
        del requested, bounds, crs, target_resolution
        load_calls.append(
            {
                "paths": list(paths),
                "time_axis": np.asarray(time_axis, dtype="datetime64[D]").copy(),
            }
        )
        coords = {
            "time": np.asarray(time_axis, dtype="datetime64[D]"),
            "band": ["Band3"],
            "y": [0],
            "x": [0],
        }
        base = np.arange(len(time_axis), dtype=np.float32).reshape(len(time_axis), 1, 1, 1)
        unc = np.full_like(base, 0.1, dtype=np.float32)
        return BRDFKernelWeights(
            f0=xr.DataArray(base, dims=["time", "band", "y", "x"], coords=coords),
            f1=xr.DataArray(np.zeros_like(base), dims=["time", "band", "y", "x"], coords=coords),
            f2=xr.DataArray(np.zeros_like(base), dims=["time", "band", "y", "x"], coords=coords),
            f0_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
            f1_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
            f2_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
            reflectance_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
        )

    monkeypatch.setattr(provider, "_load_temporal_from_granules", _load_temporal_from_granules)

    outputs = provider.get_temporal_brdf_parameters_batch(
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        obs_times=[datetime(2024, 1, 4, 12, 0, 0), datetime(2024, 2, 4, 12, 0, 0)],
        target_resolution=500.0,
        bands=[provider.source_bands[0]],
        temporal_windows=[10, 10],
        sample_date_sets=[
            [datetime(2024, 1, 1), datetime(2024, 1, 8)],
            [datetime(2024, 2, 1), datetime(2024, 2, 8)],
        ],
    )

    assert len(load_calls) == 1
    assert load_calls[0]["paths"] == downloaded
    assert list(load_calls[0]["time_axis"]) == [
        np.datetime64("2024-01-01"),
        np.datetime64("2024-01-08"),
        np.datetime64("2024-02-01"),
        np.datetime64("2024-02-08"),
    ]
    assert list(outputs[0].f0.coords["time"].values) == [
        np.datetime64("2024-01-01"),
        np.datetime64("2024-01-08"),
    ]
    assert list(outputs[1].f0.coords["time"].values) == [
        np.datetime64("2024-02-01"),
        np.datetime64("2024-02-08"),
    ]
    np.testing.assert_allclose(
        outputs[0].f0[:, 0, 0, 0].values, np.array([0.0, 1.0], dtype=np.float32)
    )
    np.testing.assert_allclose(
        outputs[1].f0[:, 0, 0, 0].values, np.array([2.0, 3.0], dtype=np.float32)
    )


def test_provider_accepts_unbounded_max_granules(monkeypatch: pytest.MonkeyPatch) -> None:
    source, catalog = _install_runtime(monkeypatch)
    provider = _DummyProvider(source=source, catalog=catalog, probe_earthdata=True, max_granules=-1)

    provider._search_granules(
        short_name="SN:dummy_brdf",
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        temporal=("2024-01-01T00:00:00Z", "2024-12-31T23:59:59Z"),
    )

    assert provider.max_granules == -1
    assert source.search_calls[-1]["count"] == -1


def test_batched_temporal_brdf_slicing_preserves_missing_requested_days(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_runtime(monkeypatch)
    provider = _DummyProvider(probe_earthdata=False)

    coords = {
        "time": np.array(["2024-01-01", "2024-01-08"], dtype="datetime64[D]"),
        "band": ["Band3"],
        "y": [0],
        "x": [0],
    }
    base = np.array([[[[0.1]]], [[[0.2]]]], dtype=np.float32)
    unc = np.full_like(base, 0.05, dtype=np.float32)
    weights = BRDFKernelWeights(
        f0=xr.DataArray(base, dims=["time", "band", "y", "x"], coords=coords),
        f1=xr.DataArray(np.zeros_like(base), dims=["time", "band", "y", "x"], coords=coords),
        f2=xr.DataArray(np.zeros_like(base), dims=["time", "band", "y", "x"], coords=coords),
        f0_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
        f1_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
        f2_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
        reflectance_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
    )

    sliced = provider._slice_temporal_weights(
        weights,
        source_time_axis=np.array(["2024-01-01", "2024-01-08"], dtype="datetime64[D]"),
        request_time_axis=np.array(
            ["2024-01-01", "2024-01-04", "2024-01-08"], dtype="datetime64[D]"
        ),
    )

    assert list(sliced.f0.coords["time"].values) == [
        np.datetime64("2024-01-01"),
        np.datetime64("2024-01-04"),
        np.datetime64("2024-01-08"),
    ]
    np.testing.assert_allclose(
        sliced.f0[[0, 2], 0, 0, 0].values,
        np.array([0.1, 0.2], dtype=np.float32),
    )
    assert np.isnan(float(sliced.f0[1, 0, 0, 0].values))


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
    np.testing.assert_allclose(
        weights.f0.sel(band="B02").values[0], np.array([0.20, 0.4], dtype=np.float32)
    )
    np.testing.assert_allclose(
        weights.f1.sel(band="B02").values[0], np.array([0.05, 0.2], dtype=np.float32)
    )
    np.testing.assert_allclose(
        weights.f2.sel(band="B02").values[0], np.array([0.02, 0.06], dtype=np.float32)
    )
    assert float(weights.reflectance_unc.sel(band="B02").values[0, 0]) == pytest.approx(0.015)
    assert float(weights.reflectance_unc.sel(band="B02").values[1, 0]) == pytest.approx(0.08)


def test_load_from_granules_batches_requested_bands_into_one_payload_merge(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_runtime(monkeypatch)
    provider = _DummyProvider(probe_earthdata=False)

    def _stack_for_band(_path: Path, product_band: ProductBandDefinition):
        base = 0.1 if product_band.label == "Band3" else 0.2
        params = (
            _da(np.full((2, 2), base, dtype=np.float32)),
            _da(np.full((2, 2), base / 2.0, dtype=np.float32)),
            _da(np.full((2, 2), base / 4.0, dtype=np.float32)),
        )
        qa = _da(np.full((2, 2), 0.0 if product_band.label == "Band3" else 1.0, dtype=np.float32))
        return params, qa

    merge_calls: list[tuple[tuple[str, ...], int]] = []

    def _merge(arrays, **_kwargs):  # type: ignore[no-untyped-def]
        merge_calls.append((arrays[0].dims, int(arrays[0].sizes["layer"])))
        return arrays[0]

    monkeypatch.setattr(provider, "_load_native_band_stack", _stack_for_band)
    monkeypatch.setattr(provider, "_merge_reprojected_tiles", _merge)

    weights = provider._load_from_granules(
        [Path("/tmp/day1.hdf")],
        requested=[("B02", provider._product_bands[0]), ("B03", provider._product_bands[1])],
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        target_resolution=500.0,
    )

    assert merge_calls == [(("layer", "y", "x"), 8)]
    assert list(weights.f0.coords["band"].values) == ["B02", "B03"]
    np.testing.assert_allclose(
        weights.f0.sel(band="B02").values, np.full((2, 2), 0.1, dtype=np.float32)
    )
    np.testing.assert_allclose(
        weights.f0.sel(band="B03").values, np.full((2, 2), 0.2, dtype=np.float32)
    )
    assert float(weights.reflectance_unc.sel(band="B02").values[0, 0]) == pytest.approx(0.015)
    assert float(weights.reflectance_unc.sel(band="B03").values[0, 0]) == pytest.approx(
        0.015 * (2.0**1.6)
    )


def test_payload_stack_round_trips_through_real_merge_helper() -> None:
    provider = _DummyProvider(probe_earthdata=False)
    band_coords = xr.IndexVariable("band", ["B02"])
    param_coords = {
        "band": band_coords,
        "parameter": xr.IndexVariable("parameter", ["f0", "f1", "f2"]),
        "y": [1.0, 0.0],
        "x": [0.0, 1.0],
    }
    right_coords = {
        "band": band_coords,
        "parameter": xr.IndexVariable("parameter", ["f0", "f1", "f2"]),
        "y": [1.0, 0.0],
        "x": [1.0, 2.0],
    }

    def _native_with_transform(
        values: np.ndarray, *, dims: tuple[str, ...], coords: dict[str, object]
    ) -> xr.DataArray:
        x_values = np.asarray(coords["x"], dtype=np.float64)
        y_values = np.asarray(coords["y"], dtype=np.float64)
        x_step = float(x_values[1] - x_values[0]) if x_values.size > 1 else 1.0
        y_step = float(y_values[1] - y_values[0]) if y_values.size > 1 else -1.0
        return (
            xr.DataArray(values, dims=dims, coords=coords)
            .rio.set_spatial_dims(x_dim="x", y_dim="y")
            .rio.write_crs("EPSG:4326")
            .rio.write_transform(
                Affine(
                    x_step,
                    0.0,
                    float(x_values[0] - (x_step / 2.0)),
                    0.0,
                    y_step,
                    float(y_values[0] - (y_step / 2.0)),
                )
            )
        )

    left_params = _native_with_transform(
        np.array(
            [
                [
                    [[1.0, 2.0], [3.0, 4.0]],
                    [[10.0, 20.0], [30.0, 40.0]],
                    [[100.0, 200.0], [300.0, 400.0]],
                ],
            ],
            dtype=np.float32,
        ),
        dims=("band", "parameter", "y", "x"),
        coords=param_coords,
    )
    left_unc = _native_with_transform(
        np.array([[[0.1, 0.2], [0.3, 0.4]]], dtype=np.float32),
        dims=("band", "y", "x"),
        coords={"band": band_coords, "y": [1.0, 0.0], "x": [0.0, 1.0]},
    )
    right_params = _native_with_transform(
        np.array(
            [
                [
                    [[9.0, 10.0], [11.0, 12.0]],
                    [[90.0, 100.0], [110.0, 120.0]],
                    [[900.0, 1000.0], [1100.0, 1200.0]],
                ],
            ],
            dtype=np.float32,
        ),
        dims=("band", "parameter", "y", "x"),
        coords=right_coords,
    )
    right_unc = _native_with_transform(
        np.array([[[0.9, 1.0], [1.1, 1.2]]], dtype=np.float32),
        dims=("band", "y", "x"),
        coords={"band": band_coords, "y": [1.0, 0.0], "x": [1.0, 2.0]},
    )

    merged = earth_adapter.merge_reprojected_tiles(
        [
            provider._pack_payload_stack(left_params, left_unc),
            provider._pack_payload_stack(right_params, right_unc),
        ],
        bounds=(-0.5, -0.5, 2.5, 1.5),
        crs="EPSG:4326",
        resolution=1.0,
        resampling="nearest",
        nodata=np.nan,
    )
    params, unc = provider._unpack_payload_stack(
        merged, requested=[("B02", provider._product_bands[0])]
    )

    assert merged.dims == ("layer", "y", "x")
    assert params.rio.crs is not None
    assert params.rio.crs.to_string() == "EPSG:4326"
    assert unc.rio.crs is not None
    assert unc.rio.crs.to_string() == "EPSG:4326"
    np.testing.assert_allclose(
        params.sel(band="B02", parameter="f0").values,
        np.array([[1.0, 2.0, 10.0], [3.0, 4.0, 12.0]], dtype=np.float32),
    )
    np.testing.assert_allclose(
        params.sel(band="B02", parameter="f1").values,
        np.array([[10.0, 20.0, 100.0], [30.0, 40.0, 120.0]], dtype=np.float32),
    )
    np.testing.assert_allclose(
        params.sel(band="B02", parameter="f2").values,
        np.array([[100.0, 200.0, 1000.0], [300.0, 400.0, 1200.0]], dtype=np.float32),
    )
    np.testing.assert_allclose(
        unc.sel(band="B02").values, np.array([[0.1, 0.2, 1.0], [0.3, 0.4, 1.2]], dtype=np.float32)
    )


def test_stack_parameter_provider_uses_direct_vrt_payload_path(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_runtime(monkeypatch)
    provider = MCD43EarthAccessProvider(probe_earthdata=False)
    requested = [("B02", provider._product_bands[0]), ("B03", provider._product_bands[1])]
    source_groups_seen: list[list[str]] = []
    band_counts_seen: list[int] = []
    resampling_seen: list[object] = []

    monkeypatch.setattr(mcd_mod, "gdal_available", lambda: True)
    monkeypatch.setattr(
        provider,
        "_read_named_dataset_attrs",
        lambda _path, dataset_names: {dataset_name: {} for dataset_name in dataset_names},
    )
    monkeypatch.setattr(
        mcd_mod,
        "resolve_gdal_subdataset_path",
        lambda path, dataset_name: f"{Path(path).name}:{dataset_name}",
    )
    monkeypatch.setattr(
        provider,
        "_load_native_band_stack",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("array fallback should not be used")
        ),
    )
    monkeypatch.setattr(
        mcd_mod,
        "build_source_aligned_target_template",
        lambda *_args, **_kwargs: (
            xr.DataArray(
                np.zeros((1, 1), dtype=np.float32),
                dims=("y", "x"),
                coords={"y": [0.0], "x": [0.0]},
            )
            .rio.set_spatial_dims(x_dim="x", y_dim="y")
            .rio.write_crs("EPSG:32615"),
            500.0,
        ),
    )

    def _fake_vrt_reader(source_groups, *, group_band_counts, **kwargs):  # type: ignore[no-untyped-def]
        source_groups_seen[:] = [list(group) for group in source_groups]
        band_counts_seen[:] = list(group_band_counts)
        target = kwargs["target_template"]
        values = np.array(
            [
                np.full((1, 1), 10.0, dtype=np.float32),
                np.full((1, 1), 11.0, dtype=np.float32),
                np.full((1, 1), 12.0, dtype=np.float32),
                np.full((1, 1), 0.0, dtype=np.float32),
                np.full((1, 1), 20.0, dtype=np.float32),
                np.full((1, 1), 21.0, dtype=np.float32),
                np.full((1, 1), 22.0, dtype=np.float32),
                np.full((1, 1), 1.0, dtype=np.float32),
            ],
            dtype=np.float32,
        )
        return xr.DataArray(
            values,
            dims=("layer", "y", "x"),
            coords={
                "layer": np.arange(values.shape[0], dtype=np.int32),
                "y": target.coords["y"],
                "x": target.coords["x"],
            },
        )

    def _capturing_vrt_reader(source_groups, *, group_band_counts, **kwargs):  # type: ignore[no-untyped-def]
        resampling_seen.append(_resampling_name(kwargs["resampling"]))
        return _fake_vrt_reader(source_groups, group_band_counts=group_band_counts, **kwargs)

    monkeypatch.setattr(mcd_mod, "read_virtual_stack_to_target", _capturing_vrt_reader)

    weights = provider._load_from_granules(
        [Path("/tmp/tile_a.hdf"), Path("/tmp/tile_b.hdf")],
        requested=requested,
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        target_resolution=500.0,
    )

    assert band_counts_seen == [3, 1, 3, 1]
    assert source_groups_seen == [
        [
            f"tile_a.hdf:{provider._product_bands[0].parameter_dataset}",
            f"tile_b.hdf:{provider._product_bands[0].parameter_dataset}",
        ],
        [
            f"tile_a.hdf:{provider._product_bands[0].qa_dataset}",
            f"tile_b.hdf:{provider._product_bands[0].qa_dataset}",
        ],
        [
            f"tile_a.hdf:{provider._product_bands[1].parameter_dataset}",
            f"tile_b.hdf:{provider._product_bands[1].parameter_dataset}",
        ],
        [
            f"tile_a.hdf:{provider._product_bands[1].qa_dataset}",
            f"tile_b.hdf:{provider._product_bands[1].qa_dataset}",
        ],
    ]
    assert resampling_seen == ["nearest"]
    np.testing.assert_allclose(
        weights.f0.sel(band="B02").values, np.full((1, 1), 10.0, dtype=np.float32)
    )
    np.testing.assert_allclose(
        weights.f2.sel(band="B03").values, np.full((1, 1), 22.0, dtype=np.float32)
    )
    assert weights.f0.rio.crs is not None
    assert weights.f0.rio.crs.to_string() == "EPSG:32615"
    assert float(weights.reflectance_unc.sel(band="B02").values[0, 0]) == pytest.approx(0.015)
    assert float(weights.reflectance_unc.sel(band="B03").values[0, 0]) == pytest.approx(
        0.015 * (2.0**1.6)
    )


def test_mcd19_provider_uses_direct_vrt_payload_path(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_runtime(monkeypatch)
    provider = MCD19EarthAccessProvider(probe_earthdata=False)
    requested = [("B02", provider._product_bands[0]), ("B03", provider._product_bands[1])]
    source_groups_seen: list[list[str]] = []
    band_counts_seen: list[int] = []
    resampling_seen: list[object] = []

    monkeypatch.setattr(mcd_mod, "gdal_available", lambda: True)
    monkeypatch.setattr(
        mcd_mod,
        "read_hdf4_datasets_attrs",
        lambda _path, dataset_names: {dataset_name: {} for dataset_name in dataset_names},
    )
    monkeypatch.setattr(
        mcd_mod,
        "resolve_gdal_subdataset_path",
        lambda path, dataset_name: f"{Path(path).name}:{dataset_name}",
    )
    monkeypatch.setattr(
        provider,
        "_load_native_band_stack",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("array fallback should not be used")
        ),
    )
    monkeypatch.setattr(
        mcd_mod,
        "build_source_aligned_target_template",
        lambda *_args, **_kwargs: (
            xr.DataArray(
                np.zeros((1, 1), dtype=np.float32),
                dims=("y", "x"),
                coords={"y": [0.0], "x": [0.0]},
            )
            .rio.set_spatial_dims(x_dim="x", y_dim="y")
            .rio.write_crs("EPSG:32615"),
            500.0,
        ),
    )

    def _fake_vrt_reader(source_groups, *, group_band_counts, **kwargs):  # type: ignore[no-untyped-def]
        source_groups_seen[:] = [list(group) for group in source_groups]
        band_counts_seen[:] = list(group_band_counts)
        target = kwargs["target_template"]
        values = np.array(
            [
                np.full((1, 1), 10.0, dtype=np.float32),
                np.full((1, 1), 20.0, dtype=np.float32),
                np.full((1, 1), 30.0, dtype=np.float32),
                np.full((1, 1), 40.0, dtype=np.float32),
                np.full((1, 1), 50.0, dtype=np.float32),
                np.full((1, 1), 60.0, dtype=np.float32),
                np.full((1, 1), 1.0, dtype=np.float32),
            ],
            dtype=np.float32,
        )
        return xr.DataArray(
            values,
            dims=("layer", "y", "x"),
            coords={
                "layer": np.arange(values.shape[0], dtype=np.int32),
                "y": target.coords["y"],
                "x": target.coords["x"],
            },
        )

    def _capturing_vrt_reader(source_groups, *, group_band_counts, **kwargs):  # type: ignore[no-untyped-def]
        resampling_seen.append(_resampling_name(kwargs["resampling"]))
        return _fake_vrt_reader(source_groups, group_band_counts=group_band_counts, **kwargs)

    monkeypatch.setattr(mcd_mod, "read_virtual_stack_to_target", _capturing_vrt_reader)

    weights = provider._load_from_granules(
        [Path("/tmp/tile_a.hdf"), Path("/tmp/tile_b.hdf")],
        requested=requested,
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        target_resolution=500.0,
    )

    assert band_counts_seen == [1, 1, 1, 1, 1, 1, 1]
    assert source_groups_seen[-1] == [
        "tile_a.hdf:Status_QA",
        "tile_b.hdf:Status_QA",
    ]
    assert resampling_seen == ["nearest"]
    np.testing.assert_allclose(
        weights.f0.sel(band="B02").values, np.full((1, 1), 10.0, dtype=np.float32)
    )
    np.testing.assert_allclose(
        weights.f1.sel(band="B03").values, np.full((1, 1), 50.0, dtype=np.float32)
    )
    assert float(weights.reflectance_unc.sel(band="B02").values[0, 0]) == pytest.approx(
        float(weights.reflectance_unc.sel(band="B03").values[0, 0])
    )


def test_stack_parameter_provider_uses_direct_temporal_vrt_payload_path(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_runtime(monkeypatch)
    provider = MCD43EarthAccessProvider(probe_earthdata=False)
    requested = [("B02", provider._product_bands[0]), ("B03", provider._product_bands[1])]
    time_axis = np.array(["2024-01-01", "2024-01-02", "2024-01-03"], dtype="datetime64[D]")
    day_map = {
        Path("/tmp/day1_a.hdf"): datetime(2024, 1, 1),
        Path("/tmp/day1_b.hdf"): datetime(2024, 1, 1),
        Path("/tmp/day3.hdf"): datetime(2024, 1, 3),
    }
    source_groups_seen: list[list[str]] = []
    band_counts_seen: list[int] = []
    resampling_seen: list[object] = []

    monkeypatch.setattr(mcd_mod, "parse_granule_date", lambda path: day_map[Path(path)])
    monkeypatch.setattr(mcd_mod, "gdal_available", lambda: True)
    monkeypatch.setattr(
        provider,
        "_read_named_dataset_attrs",
        lambda _path, dataset_names: {dataset_name: {} for dataset_name in dataset_names},
    )
    monkeypatch.setattr(
        mcd_mod,
        "resolve_gdal_subdataset_path",
        lambda path, dataset_name: f"{Path(path).name}:{dataset_name}",
    )
    monkeypatch.setattr(
        provider,
        "_load_native_band_stack",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("array fallback should not be used")
        ),
    )
    monkeypatch.setattr(
        mcd_mod,
        "build_source_aligned_target_template",
        lambda *_args, **_kwargs: (
            xr.DataArray(
                np.zeros((1, 1), dtype=np.float32),
                dims=("y", "x"),
                coords={"y": [0.0], "x": [0.0]},
            )
            .rio.set_spatial_dims(x_dim="x", y_dim="y")
            .rio.write_crs("EPSG:32615"),
            500.0,
        ),
    )

    def _fake_vrt_reader(source_groups, *, group_band_counts, **kwargs):  # type: ignore[no-untyped-def]
        source_groups_seen[:] = [list(group) for group in source_groups]
        band_counts_seen[:] = list(group_band_counts)
        target = kwargs["target_template"]
        values = np.array(
            [
                np.full((1, 1), 10.0, dtype=np.float32),
                np.full((1, 1), 11.0, dtype=np.float32),
                np.full((1, 1), 12.0, dtype=np.float32),
                np.full((1, 1), 0.0, dtype=np.float32),
                np.full((1, 1), 20.0, dtype=np.float32),
                np.full((1, 1), 21.0, dtype=np.float32),
                np.full((1, 1), 22.0, dtype=np.float32),
                np.full((1, 1), 1.0, dtype=np.float32),
                np.full((1, 1), 30.0, dtype=np.float32),
                np.full((1, 1), 31.0, dtype=np.float32),
                np.full((1, 1), 32.0, dtype=np.float32),
                np.full((1, 1), 2.0, dtype=np.float32),
                np.full((1, 1), 40.0, dtype=np.float32),
                np.full((1, 1), 41.0, dtype=np.float32),
                np.full((1, 1), 42.0, dtype=np.float32),
                np.full((1, 1), 3.0, dtype=np.float32),
            ],
            dtype=np.float32,
        )
        return xr.DataArray(
            values,
            dims=("layer", "y", "x"),
            coords={
                "layer": np.arange(values.shape[0], dtype=np.int32),
                "y": target.coords["y"],
                "x": target.coords["x"],
            },
        )

    def _capturing_vrt_reader(source_groups, *, group_band_counts, **kwargs):  # type: ignore[no-untyped-def]
        resampling_seen.append(_resampling_name(kwargs["resampling"]))
        return _fake_vrt_reader(source_groups, group_band_counts=group_band_counts, **kwargs)

    monkeypatch.setattr(mcd_mod, "read_virtual_stack_to_target", _capturing_vrt_reader)

    temporal = provider._load_temporal_from_granules(
        [Path("/tmp/day1_a.hdf"), Path("/tmp/day1_b.hdf"), Path("/tmp/day3.hdf")],
        requested=requested,
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        target_resolution=500.0,
        time_axis=time_axis,
    )

    assert band_counts_seen == [3, 1, 3, 1, 3, 1, 3, 1]
    assert source_groups_seen[0] == [
        f"day1_a.hdf:{provider._product_bands[0].parameter_dataset}",
        f"day1_b.hdf:{provider._product_bands[0].parameter_dataset}",
    ]
    assert source_groups_seen[4] == [f"day3.hdf:{provider._product_bands[0].parameter_dataset}"]
    assert resampling_seen == ["nearest"]
    np.testing.assert_allclose(
        temporal.f0.sel(time=np.datetime64("2024-01-01"), band="B02").values,
        np.full((1, 1), 10.0, dtype=np.float32),
    )
    np.testing.assert_allclose(
        temporal.f2.sel(time=np.datetime64("2024-01-03"), band="B03").values,
        np.full((1, 1), 42.0, dtype=np.float32),
    )
    assert temporal.f0.rio.crs is not None
    assert temporal.f0.rio.crs.to_string() == "EPSG:32615"
    assert np.isnan(temporal.f0.sel(time=np.datetime64("2024-01-02")).values).all()
    assert float(
        temporal.reflectance_unc.sel(time=np.datetime64("2024-01-03"), band="B02").values[0, 0]
    ) == pytest.approx(0.015 * (3.0**1.6))


def test_stack_parameter_provider_daily_payload_fallback_uses_one_final_warp(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_runtime(monkeypatch)
    provider = MCD43EarthAccessProvider(probe_earthdata=False)
    requested = [("B02", provider._product_bands[0]), ("B03", provider._product_bands[1])]
    time_axis = np.array(["2024-01-01", "2024-01-02", "2024-01-03"], dtype="datetime64[D]")
    day_map = {
        Path("/tmp/day1_a.hdf"): datetime(2024, 1, 1),
        Path("/tmp/day1_b.hdf"): datetime(2024, 1, 1),
        Path("/tmp/day3.hdf"): datetime(2024, 1, 3),
    }
    build_calls: list[tuple[list[list[str]], list[int]]] = []
    final_calls: list[tuple[list[list[str]], list[int]]] = []
    unlink_calls: list[str] = []
    resampling_seen: list[object] = []

    monkeypatch.setattr(mcd_mod, "parse_granule_date", lambda path: day_map[Path(path)])
    monkeypatch.setattr(mcd_mod, "gdal_available", lambda: True)

    def _no_temporal_payload(_grouped_paths, **_kwargs):  # noqa: ANN001, ANN003
        return None

    monkeypatch.setattr(provider, "_load_temporal_payload_vrt", _no_temporal_payload)
    monkeypatch.setattr(
        provider,
        "_read_named_dataset_attrs",
        lambda _path, dataset_names: {dataset_name: {} for dataset_name in dataset_names},
    )
    monkeypatch.setattr(
        mcd_mod,
        "resolve_gdal_subdataset_path",
        lambda path, dataset_name: f"{Path(path).name}:{dataset_name}",
    )
    monkeypatch.setattr(
        provider,
        "_merge_requested_payload",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("per-day merges should not be used")
        ),
    )
    monkeypatch.setattr(
        mcd_mod,
        "build_source_aligned_target_template",
        lambda *_args, **_kwargs: (
            xr.DataArray(
                np.zeros((1, 1), dtype=np.float32),
                dims=("y", "x"),
                coords={"y": [0.0], "x": [0.0]},
            )
            .rio.set_spatial_dims(x_dim="x", y_dim="y")
            .rio.write_crs("EPSG:32615"),
            500.0,
        ),
    )

    monkeypatch.setitem(
        sys.modules,
        "osgeo",
        SimpleNamespace(gdal=SimpleNamespace(Unlink=lambda path: unlink_calls.append(path))),
    )

    def _fake_build_virtual_stack_vrt(source_groups, *, group_band_counts, **_kwargs):  # type: ignore[no-untyped-def]
        build_calls.append(([list(group) for group in source_groups], list(group_band_counts)))
        index = len(build_calls) - 1
        return SimpleNamespace(
            path=f"/vsimem/day_{index}.vrt",
            cleanup_paths=(f"/vsimem/day_{index}.vrt", f"/vsimem/day_{index}_group_0.vrt"),
            expected_layers=int(sum(group_band_counts)),
        )

    def _fake_vrt_reader(source_groups, *, group_band_counts, **kwargs):  # type: ignore[no-untyped-def]
        final_calls.append(([list(group) for group in source_groups], list(group_band_counts)))
        target = kwargs["target_template"]
        values = np.array(
            [
                np.full((1, 1), 10.0, dtype=np.float32),
                np.full((1, 1), 11.0, dtype=np.float32),
                np.full((1, 1), 12.0, dtype=np.float32),
                np.full((1, 1), 0.0, dtype=np.float32),
                np.full((1, 1), 20.0, dtype=np.float32),
                np.full((1, 1), 21.0, dtype=np.float32),
                np.full((1, 1), 22.0, dtype=np.float32),
                np.full((1, 1), 1.0, dtype=np.float32),
                np.full((1, 1), 30.0, dtype=np.float32),
                np.full((1, 1), 31.0, dtype=np.float32),
                np.full((1, 1), 32.0, dtype=np.float32),
                np.full((1, 1), 2.0, dtype=np.float32),
                np.full((1, 1), 40.0, dtype=np.float32),
                np.full((1, 1), 41.0, dtype=np.float32),
                np.full((1, 1), 42.0, dtype=np.float32),
                np.full((1, 1), 3.0, dtype=np.float32),
            ],
            dtype=np.float32,
        )
        return xr.DataArray(
            values,
            dims=("layer", "y", "x"),
            coords={
                "layer": np.arange(values.shape[0], dtype=np.int32),
                "y": target.coords["y"],
                "x": target.coords["x"],
            },
        )

    monkeypatch.setattr(mcd_mod, "_build_virtual_stack_vrt", _fake_build_virtual_stack_vrt)

    def _capturing_vrt_reader(source_groups, *, group_band_counts, **kwargs):  # type: ignore[no-untyped-def]
        resampling_seen.append(_resampling_name(kwargs["resampling"]))
        return _fake_vrt_reader(source_groups, group_band_counts=group_band_counts, **kwargs)

    monkeypatch.setattr(mcd_mod, "read_virtual_stack_to_target", _capturing_vrt_reader)

    temporal = provider._load_temporal_from_granules(
        [Path("/tmp/day1_a.hdf"), Path("/tmp/day1_b.hdf"), Path("/tmp/day3.hdf")],
        requested=requested,
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        target_resolution=500.0,
        time_axis=time_axis,
    )

    assert len(build_calls) == 2
    assert build_calls[0][1] == [3, 1, 3, 1]
    assert build_calls[1][1] == [3, 1, 3, 1]
    assert final_calls == [([["/vsimem/day_0.vrt"], ["/vsimem/day_1.vrt"]], [8, 8])]
    assert resampling_seen == ["nearest"]
    assert unlink_calls == [
        "/vsimem/day_0.vrt",
        "/vsimem/day_0_group_0.vrt",
        "/vsimem/day_1.vrt",
        "/vsimem/day_1_group_0.vrt",
    ]
    np.testing.assert_allclose(
        temporal.f0.sel(time=np.datetime64("2024-01-01"), band="B02").values,
        np.full((1, 1), 10.0, dtype=np.float32),
    )
    np.testing.assert_allclose(
        temporal.f2.sel(time=np.datetime64("2024-01-03"), band="B03").values,
        np.full((1, 1), 42.0, dtype=np.float32),
    )
    assert np.isnan(temporal.f0.sel(time=np.datetime64("2024-01-02")).values).all()


def test_stack_parameter_provider_skips_oversized_temporal_vrt_and_logs_per_day_progress(
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    _install_runtime(monkeypatch)
    provider = MCD43EarthAccessProvider(probe_earthdata=False)
    requested = [("B02", provider._product_bands[0]), ("B03", provider._product_bands[1])]
    time_axis = np.array(["2024-01-01", "2024-01-02", "2024-01-03"], dtype="datetime64[D]")
    day_map = {
        Path("/tmp/day1_a.hdf"): datetime(2024, 1, 1),
        Path("/tmp/day1_b.hdf"): datetime(2024, 1, 1),
        Path("/tmp/day3.hdf"): datetime(2024, 1, 3),
    }
    build_counter = {"count": 0}
    merge_calls: list[list[str]] = []
    unlink_calls: list[str] = []

    monkeypatch.setattr(mcd_mod, "parse_granule_date", lambda path: day_map[Path(path)])
    monkeypatch.setattr(mcd_mod, "gdal_available", lambda: True)
    monkeypatch.setattr(mcd_mod, "_MAX_ONE_SHOT_TEMPORAL_VRT_OUTPUT_BYTES", 0)
    monkeypatch.setattr(
        provider,
        "_read_named_dataset_attrs",
        lambda _path, dataset_names: {dataset_name: {} for dataset_name in dataset_names},
    )
    monkeypatch.setitem(
        sys.modules,
        "osgeo",
        SimpleNamespace(gdal=SimpleNamespace(Unlink=lambda path: unlink_calls.append(path))),
    )

    def _fake_build_virtual_stack_vrt(source_groups, *, group_band_counts, **_kwargs):  # type: ignore[no-untyped-def]
        index = build_counter["count"]
        build_counter["count"] += 1
        return SimpleNamespace(
            path=f"/vsimem/day_{index}.vrt",
            cleanup_paths=(f"/vsimem/day_{index}.vrt",),
            expected_layers=int(sum(group_band_counts)),
        )

    monkeypatch.setattr(mcd_mod, "_build_virtual_stack_vrt", _fake_build_virtual_stack_vrt)
    monkeypatch.setattr(
        mcd_mod,
        "read_virtual_stack_to_target",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("oversized one-shot VRT should be skipped")
        ),
    )

    def _fake_merge_requested_payload(
        day_paths,  # type: ignore[no-untyped-def]
        *,
        requested,
        target_template,
        **_kwargs,
    ):
        merge_calls.append([Path(path).name for path in day_paths])
        base = 10.0 if "day1" in Path(day_paths[0]).name else 30.0
        params_values = np.array(
            [
                [[[base]], [[base + 1.0]], [[base + 2.0]]],
                [[[base + 10.0]], [[base + 11.0]], [[base + 12.0]]],
            ],
            dtype=np.float32,
        )
        unc_values = np.full((len(requested), 1, 1), 0.05, dtype=np.float32)
        band_coords = [band_coord for band_coord, _ in requested]
        params = provider._target_array_like(
            target_template,
            params_values,
            dims=("band", "parameter", "y", "x"),
            coords={
                "band": xr.IndexVariable("band", band_coords),
                "parameter": xr.IndexVariable("parameter", ["f0", "f1", "f2"]),
                "y": target_template.coords["y"],
                "x": target_template.coords["x"],
            },
        )
        unc = provider._target_array_like(
            target_template,
            unc_values,
            dims=("band", "y", "x"),
            coords={
                "band": xr.IndexVariable("band", band_coords),
                "y": target_template.coords["y"],
                "x": target_template.coords["x"],
            },
        )
        return provider._pack_payload_stack(params, unc)

    monkeypatch.setattr(provider, "_merge_requested_payload", _fake_merge_requested_payload)

    with caplog.at_level(logging.INFO):
        temporal = provider._load_temporal_from_granules(
            [Path("/tmp/day1_a.hdf"), Path("/tmp/day1_b.hdf"), Path("/tmp/day3.hdf")],
            requested=requested,
            bounds=(0.0, 0.0, 1.0, 1.0),
            crs="EPSG:4326",
            target_resolution=500.0,
            time_axis=time_axis,
        )

    assert merge_calls == [["day1_a.hdf", "day1_b.hdf"], ["day3.hdf"]]
    assert "estimated warped payload" in caplog.text
    assert "per-day fallback 1/2 day=2024-01-01" in caplog.text
    assert "per-day fallback 2/2 day=2024-01-03" in caplog.text
    np.testing.assert_allclose(
        temporal.f0.sel(time=np.datetime64("2024-01-01"), band="B02").values,
        np.full((1, 1), 10.0, dtype=np.float32),
    )
    np.testing.assert_allclose(
        temporal.f2.sel(time=np.datetime64("2024-01-03"), band="B03").values,
        np.full((1, 1), 42.0, dtype=np.float32),
    )
    assert np.isnan(temporal.f0.sel(time=np.datetime64("2024-01-02")).values).all()


def test_mcd19_provider_uses_direct_temporal_vrt_payload_path(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_runtime(monkeypatch)
    provider = MCD19EarthAccessProvider(probe_earthdata=False)
    requested = [("B02", provider._product_bands[0]), ("B03", provider._product_bands[1])]
    time_axis = np.array(["2024-01-01", "2024-01-02", "2024-01-03"], dtype="datetime64[D]")
    day_map = {
        Path("/tmp/day1.hdf"): datetime(2024, 1, 1),
        Path("/tmp/day3.hdf"): datetime(2024, 1, 3),
    }
    source_groups_seen: list[list[str]] = []
    band_counts_seen: list[int] = []

    monkeypatch.setattr(mcd_mod, "parse_granule_date", lambda path: day_map[Path(path)])
    monkeypatch.setattr(mcd_mod, "gdal_available", lambda: True)
    monkeypatch.setattr(
        mcd_mod,
        "read_hdf4_datasets_attrs",
        lambda _path, dataset_names: {dataset_name: {} for dataset_name in dataset_names},
    )
    monkeypatch.setattr(
        mcd_mod,
        "resolve_gdal_subdataset_path",
        lambda path, dataset_name: f"{Path(path).name}:{dataset_name}",
    )
    monkeypatch.setattr(
        provider,
        "_load_native_band_stack",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("array fallback should not be used")
        ),
    )

    def _fake_vrt_reader(source_groups, *, group_band_counts, **kwargs):  # type: ignore[no-untyped-def]
        source_groups_seen[:] = [list(group) for group in source_groups]
        band_counts_seen[:] = list(group_band_counts)
        target = kwargs["target_template"]
        values = np.array(
            [
                np.full((1, 1), 10.0, dtype=np.float32),
                np.full((1, 1), 20.0, dtype=np.float32),
                np.full((1, 1), 30.0, dtype=np.float32),
                np.full((1, 1), 40.0, dtype=np.float32),
                np.full((1, 1), 50.0, dtype=np.float32),
                np.full((1, 1), 60.0, dtype=np.float32),
                np.full((1, 1), 1.0, dtype=np.float32),
                np.full((1, 1), 70.0, dtype=np.float32),
                np.full((1, 1), 80.0, dtype=np.float32),
                np.full((1, 1), 90.0, dtype=np.float32),
                np.full((1, 1), 100.0, dtype=np.float32),
                np.full((1, 1), 110.0, dtype=np.float32),
                np.full((1, 1), 120.0, dtype=np.float32),
                np.full((1, 1), 2.0, dtype=np.float32),
            ],
            dtype=np.float32,
        )
        return xr.DataArray(
            values,
            dims=("layer", "y", "x"),
            coords={
                "layer": np.arange(values.shape[0], dtype=np.int32),
                "y": target.coords["y"],
                "x": target.coords["x"],
            },
        )

    monkeypatch.setattr(mcd_mod, "read_virtual_stack_to_target", _fake_vrt_reader)

    temporal = provider._load_temporal_from_granules(
        [Path("/tmp/day1.hdf"), Path("/tmp/day3.hdf")],
        requested=requested,
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        target_resolution=500.0,
        time_axis=time_axis,
    )

    assert band_counts_seen == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    assert source_groups_seen[6] == ["day1.hdf:Status_QA"]
    assert source_groups_seen[-1] == ["day3.hdf:Status_QA"]
    np.testing.assert_allclose(
        temporal.f1.sel(time=np.datetime64("2024-01-01"), band="B03").values,
        np.full((1, 1), 50.0, dtype=np.float32),
    )
    np.testing.assert_allclose(
        temporal.f2.sel(time=np.datetime64("2024-01-03"), band="B03").values,
        np.full((1, 1), 120.0, dtype=np.float32),
    )
    assert temporal.f0.rio.crs is not None
    assert temporal.f0.rio.crs.to_string() == "EPSG:32615"
    assert np.isnan(temporal.f0.sel(time=np.datetime64("2024-01-02")).values).all()
    assert float(
        temporal.reflectance_unc.sel(time=np.datetime64("2024-01-03"), band="B03").values[0, 0]
    ) == pytest.approx(0.015 * (3.0**1.6))


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
    monkeypatch.setattr(
        provider,
        "_grid",
        lambda _bounds, _resolution: (np.array([10.0, 20.0]), np.array([30.0, 40.0])),
    )
    monkeypatch.setattr(
        mcd_mod,
        "build_source_aligned_target_template",
        lambda *_args, **kwargs: (
            xr.DataArray(
                np.zeros((2, 2), dtype=np.float32),
                dims=("y", "x"),
                coords={"y": [10.0, 20.0], "x": [30.0, 40.0]},
            )
            .rio.set_spatial_dims(x_dim="x", y_dim="y")
            .rio.write_crs("EPSG:32615"),
            float(kwargs.get("resolution") or 500.0),
        ),
    )

    def _stack_for_path(path: Path, _band: ProductBandDefinition):
        base = 0.1 if "day1" in path.name else 0.3
        coords = {
            "y": [0.0, 1.0],
            "x": [0.0, 1.0],
            "extra": xr.DataArray([1.0, 2.0], dims=("y",)),
        }
        params = (
            xr.DataArray(np.full((2, 2), base, dtype=np.float32), dims=("y", "x"), coords=coords),
            xr.DataArray(
                np.full((2, 2), base / 2.0, dtype=np.float32), dims=("y", "x"), coords=coords
            ),
            xr.DataArray(
                np.full((2, 2), base / 4.0, dtype=np.float32), dims=("y", "x"), coords=coords
            ),
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
    assert temporal.f0.rio.crs is not None
    assert temporal.f0.rio.crs.to_string() == "EPSG:32615"

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
    fallback = provider._load_temporal_from_granules(
        [Path("/tmp/day1.hdf")],
        requested=[("B02", provider._product_bands[0])],
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        target_resolution=500.0,
        time_axis=time_axis[:1],
    )
    assert np.isnan(fallback.f0.values).all()
    assert np.isnan(fallback.reflectance_unc.values).all()
    assert list(fallback.f0.coords["y"].values) == [10.0, 20.0]
    assert list(fallback.f0.coords["x"].values) == [30.0, 40.0]
    assert fallback.f0.rio.crs is not None
    assert fallback.f0.rio.crs.to_string() == "EPSG:32615"


def test_default_temporal_weights_preserve_target_crs() -> None:
    provider = _DummyProvider(probe_earthdata=False)
    time_axis = np.array(["2024-01-01", "2024-01-02"], dtype="datetime64[D]")

    weights = provider._default_temporal_weights(
        bounds=(399960.0, 4590240.0, 400960.0, 4600240.0),
        resolution=500.0,
        bands=["B02"],
        time_axis=time_axis,
        crs="EPSG:32615",
    )

    assert weights.f0.rio.crs is not None
    assert weights.f0.rio.crs.to_string() == "EPSG:32615"
    assert weights.reflectance_unc is not None
    assert weights.reflectance_unc.rio.crs is not None
    assert weights.reflectance_unc.rio.crs.to_string() == "EPSG:32615"


def test_temporal_loading_fallback_batches_requested_bands_into_one_payload_merge_per_day(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_runtime(monkeypatch)
    provider = _DummyProvider(probe_earthdata=False)
    time_axis = np.array(["2024-01-01", "2024-01-02"], dtype="datetime64[D]")
    day_map = {
        Path("/tmp/day1.hdf"): datetime(2024, 1, 1),
        Path("/tmp/day2.hdf"): datetime(2024, 1, 2),
    }

    monkeypatch.setattr(mcd_mod, "parse_granule_date", lambda path: day_map[Path(path)])
    monkeypatch.setattr(
        provider,
        "_grid",
        lambda _bounds, _resolution: (np.array([10.0, 20.0]), np.array([30.0, 40.0])),
    )
    monkeypatch.setattr(
        mcd_mod,
        "build_source_aligned_target_template",
        lambda *_args, **kwargs: (
            xr.DataArray(
                np.zeros((2, 2), dtype=np.float32),
                dims=("y", "x"),
                coords={"y": [10.0, 20.0], "x": [30.0, 40.0]},
            )
            .rio.set_spatial_dims(x_dim="x", y_dim="y")
            .rio.write_crs("EPSG:32615"),
            float(kwargs.get("resolution") or 500.0),
        ),
    )

    def _stack_for_path(path: Path, product_band: ProductBandDefinition):
        day_offset = 0.1 if "day1" in path.name else 0.3
        band_offset = 0.0 if product_band.label == "Band3" else 0.05
        base = day_offset + band_offset
        params = (
            _da(np.full((2, 2), base, dtype=np.float32)),
            _da(np.full((2, 2), base / 2.0, dtype=np.float32)),
            _da(np.full((2, 2), base / 4.0, dtype=np.float32)),
        )
        qa = _da(np.full((2, 2), 0.0 if product_band.label == "Band3" else 1.0, dtype=np.float32))
        return params, qa

    merge_calls: list[tuple[tuple[str, ...], int]] = []

    def _merge(arrays, **_kwargs):  # type: ignore[no-untyped-def]
        merge_calls.append((arrays[0].dims, int(arrays[0].sizes["layer"])))
        return arrays[0]

    monkeypatch.setattr(provider, "_load_native_band_stack", _stack_for_path)
    monkeypatch.setattr(provider, "_merge_reprojected_tiles", _merge)

    temporal = provider._load_temporal_from_granules(
        [Path("/tmp/day1.hdf"), Path("/tmp/day2.hdf")],
        requested=[("B02", provider._product_bands[0]), ("B03", provider._product_bands[1])],
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        target_resolution=500.0,
        time_axis=time_axis,
    )

    assert merge_calls == [
        (("layer", "y", "x"), 8),
        (("layer", "y", "x"), 8),
    ]
    np.testing.assert_allclose(
        temporal.f0.sel(time=np.datetime64("2024-01-01"), band="B02").values,
        np.full((2, 2), 0.1, dtype=np.float32),
    )
    np.testing.assert_allclose(
        temporal.f0.sel(time=np.datetime64("2024-01-02"), band="B03").values,
        np.full((2, 2), 0.35, dtype=np.float32),
    )


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

    monkeypatch.setattr(
        MCD43EarthAccessProvider, "_read_dataset", staticmethod(_invalid_read_dataset)
    )
    monkeypatch.setattr(
        mcd_mod, "apply_scale_and_mask", lambda values, _attrs: np.asarray(values, dtype=np.float32)
    )
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

    monkeypatch.setattr(
        MCD43EarthAccessProvider, "_read_dataset", staticmethod(_valid_read_dataset)
    )
    (f0, f1, f2), qa = provider._load_native_band_stack(Path("/tmp/granule.hdf"), product_band)
    assert float(f0.mean()) == pytest.approx(10.0)
    assert float(f1.mean()) == pytest.approx(20.0)
    assert float(f2.mean()) == pytest.approx(30.0)
    assert float(qa.mean()) == pytest.approx(4.0)


def test_stack_parameter_provider_batches_requested_native_reads(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_runtime(monkeypatch)
    provider = MCD43EarthAccessProvider(probe_earthdata=False)
    requested = [("B02", provider._product_bands[0]), ("B03", provider._product_bands[1])]
    seen: list[tuple[str, ...]] = []

    def _fake_read_named_datasets(_path: str | Path, dataset_names: tuple[str, ...]):
        seen.append(dataset_names)
        outputs: dict[str, tuple[np.ndarray, dict[str, object]]] = {}
        for dataset_name in dataset_names:
            if "Parameters" in dataset_name:
                base = 10 if "Band1" in dataset_name else 20
                outputs[dataset_name] = (
                    np.dstack(
                        [
                            np.full((2, 2), base, dtype=np.int16),
                            np.full((2, 2), base + 1, dtype=np.int16),
                            np.full((2, 2), base + 2, dtype=np.int16),
                        ]
                    ),
                    {},
                )
            else:
                qa_value = 0 if "Band1" in dataset_name else 1
                outputs[dataset_name] = (np.full((2, 2), qa_value, dtype=np.int16), {})
        return outputs

    monkeypatch.setattr(provider, "_read_named_datasets", _fake_read_named_datasets)
    monkeypatch.setattr(
        mcd_mod, "apply_scale_and_mask", lambda values, _attrs: np.asarray(values, dtype=np.float32)
    )
    monkeypatch.setattr(mcd_mod, "make_native_grid_dataarray", _fake_native_grid)

    params, qa = provider._load_native_requested_stack(Path("/tmp/granule.hdf"), requested)

    assert seen == [
        (
            provider._product_bands[0].parameter_dataset,
            provider._product_bands[0].qa_dataset,
            provider._product_bands[1].parameter_dataset,
            provider._product_bands[1].qa_dataset,
        )
    ]
    assert params.dims == ("band", "parameter", "y", "x")
    assert qa.dims == ("band", "y", "x")
    assert list(params.coords["band"].values) == ["B02", "B03"]
    np.testing.assert_allclose(
        params.sel(band="B02", parameter="f0").values, np.full((2, 2), 10.0, dtype=np.float32)
    )
    np.testing.assert_allclose(
        params.sel(band="B03", parameter="f2").values, np.full((2, 2), 22.0, dtype=np.float32)
    )


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
    monkeypatch.setattr(
        mcd_mod, "apply_scale_and_mask", lambda values, _attrs: np.asarray(values, dtype=np.float32)
    )
    monkeypatch.setattr(mcd_mod, "make_native_grid_dataarray", _fake_native_grid)

    (f0, f1, f2), qa = provider._load_native_band_stack(Path("/tmp/granule.hdf"), product_band)
    assert seen == ["Kiso_Band3", "Kvol_Band3", "Kgeo_Band3", "Status_QA"]
    assert float(f0.mean()) == pytest.approx(10.0)
    assert float(f1.mean()) == pytest.approx(20.0)
    assert float(f2.mean()) == pytest.approx(30.0)
    assert float(qa.mean()) == pytest.approx(1.0)
