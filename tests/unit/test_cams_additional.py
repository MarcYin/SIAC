from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

from siac.adapters.atmo import cams as cams_mod
from siac.adapters.atmo.cams import CAMSProvider
from siac.adapters.auth import CDSES3Credentials, CredentialManager


def _dataset(var_names: tuple[str, ...]) -> xr.Dataset:
    return xr.Dataset(
        {
            name: xr.DataArray(
                np.full((2, 2), 1.0, dtype=np.float32),
                dims=["latitude", "longitude"],
                coords={"latitude": [1.0, 0.0], "longitude": [0.0, 1.0]},
            )
            for name in var_names
        }
    )


def test_calibrated_aot_uncertainty_is_tight_at_low_aod_and_saturates() -> None:
    aot = xr.DataArray(
        np.array(
            [[-0.1, 0.0, cams_mod._CAMS_AOT_UNC_HALF, 1.5]],
            dtype=np.float32,
        ),
        dims=["y", "x"],
    )

    unc = cams_mod._calibrated_aot_uncertainty(aot)
    values = unc.values[0]

    assert unc.dims == aot.dims
    assert unc.dtype == np.float32
    assert values[0] == pytest.approx(cams_mod._CAMS_AOT_UNC_FLOOR)
    assert values[1] == pytest.approx(cams_mod._CAMS_AOT_UNC_FLOOR)
    assert values[2] == pytest.approx(
        np.sqrt(cams_mod._CAMS_AOT_UNC_FLOOR**2 + 0.5 * cams_mod._CAMS_AOT_UNC_PLATEAU**2),
        rel=1e-6,
    )
    assert cams_mod._CAMS_AOT_UNC_FLOOR < values[2] < values[3]
    assert values[3] < np.sqrt(cams_mod._CAMS_AOT_UNC_FLOOR**2 + cams_mod._CAMS_AOT_UNC_PLATEAU**2)


def test_calibrate_cams_aot_center_is_gated_default_off(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    aot = xr.DataArray(
        np.array([[0.0, 0.3, 1.2]], dtype=np.float32),
        dims=["y", "x"],
    )

    monkeypatch.setattr(cams_mod, "_APPLY_AOT_CENTER_CALIBRATION", False)
    assert cams_mod._calibrate_cams_aot_center(aot) is aot

    monkeypatch.setattr(cams_mod, "_APPLY_AOT_CENTER_CALIBRATION", True)
    corrected = cams_mod._calibrate_cams_aot_center(aot)
    powered = np.maximum(aot.values, 0.0) ** 2
    expected = aot.values + cams_mod._CAMS_CENTER_GAIN * powered / (
        powered + cams_mod._CAMS_CENTER_K**2
    )

    assert corrected.dtype == np.float32
    np.testing.assert_allclose(corrected.values, expected.astype(np.float32), rtol=1e-6)
    assert corrected.values[0, 0] == pytest.approx(0.0)
    assert corrected.values[0, 2] > aot.values[0, 2]


def test_load_cams_data_handles_explicit_missing_sources_and_local_dir_creation(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    remote_explicit = CAMSProvider("https://example.com/cams_2024-01-01.nc", download_missing=True)

    def fake_load_explicit(_path):  # type: ignore[no-untyped-def]
        return None

    def fake_complete(_dataset, _time):  # type: ignore[no-untyped-def]
        return "completed"

    monkeypatch.setattr(remote_explicit, "_load_from_explicit_path", fake_load_explicit)
    monkeypatch.setattr(remote_explicit, "_complete_cams_dataset", fake_complete)
    assert remote_explicit._load_cams_data(datetime(2024, 1, 1)) == "completed"

    local_explicit = CAMSProvider(tmp_path / "missing.nc")
    monkeypatch.setattr(local_explicit, "_load_from_explicit_path", fake_load_explicit)
    with caplog.at_level("WARNING"):
        assert local_explicit._load_cams_data(datetime(2024, 1, 1)) is None
    assert "does not exist" in caplog.text

    local_explicit_download = CAMSProvider(tmp_path / "missing-download.nc", download_missing=True)
    monkeypatch.setattr(local_explicit_download, "_load_from_explicit_path", fake_load_explicit)
    monkeypatch.setattr(local_explicit_download, "_complete_cams_dataset", fake_complete)
    with caplog.at_level("WARNING"):
        assert local_explicit_download._load_cams_data(datetime(2024, 1, 1)) == "completed"
    assert "missing-download.nc" in caplog.text

    local_dir = tmp_path / "new-cams-dir"
    provider = CAMSProvider(local_dir, download_missing=True)
    monkeypatch.setattr(provider, "_load_from_explicit_path", fake_load_explicit)

    def fake_load_cds(_time):  # type: ignore[no-untyped-def]
        return "downloaded"

    monkeypatch.setattr(provider, "_load_cds_dataset", fake_load_cds)
    assert provider._load_cams_data(datetime(2024, 1, 1)) == "downloaded"
    assert local_dir.exists()


def test_load_cams_data_handles_non_directory_and_failed_local_candidates(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    not_a_dir = tmp_path / "cams_cache"
    not_a_dir.write_text("x")
    provider = CAMSProvider(not_a_dir)
    monkeypatch.setattr(provider, "_load_from_explicit_path", lambda _path: None)
    assert provider._load_cams_data(datetime(2024, 1, 1)) is None

    local_dir = tmp_path / "local_dir"
    local_dir.mkdir()
    candidate = local_dir / "cams_20240101.nc"
    candidate.write_text("x")

    fallback = CAMSProvider(local_dir, download_missing=False)
    monkeypatch.setattr(fallback, "_load_from_explicit_path", lambda _path: None)
    monkeypatch.setattr(fallback, "_load_cams_tif_group", lambda _date, _iso: None)
    monkeypatch.setattr(
        xr,
        "open_mfdataset",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("bad nc")),
    )
    assert fallback._load_cams_data(datetime(2024, 1, 1)) is None


def test_load_cams_data_uses_open_dataset_for_single_local_netcdf(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    local_dir = tmp_path / "cams"
    local_dir.mkdir()
    candidate = local_dir / "CAMS_2024-01-01.nc"
    candidate.write_text("x")
    provider = CAMSProvider(local_dir)
    expected = _dataset(("aod550", "tcwv", "gtco3"))

    monkeypatch.setattr(xr, "open_dataset", lambda path: expected if path == candidate else None)
    monkeypatch.setattr(
        xr,
        "open_mfdataset",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("open_mfdataset should not be used")
        ),
    )

    loaded = provider._load_cams_data(datetime(2024, 1, 1))
    assert loaded is expected


def test_load_cams_data_falls_back_for_remote_and_local_sources(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    remote = CAMSProvider("https://example.com/cams/", download_missing=False)
    monkeypatch.setattr(remote, "_load_from_explicit_path", lambda _path: None)
    monkeypatch.setattr(remote, "_load_from_remote_base", lambda _base, _time: None)
    assert remote._load_cams_data(datetime(2024, 1, 1)) is None

    local_dir = tmp_path / "cams"
    local_dir.mkdir()
    provider = CAMSProvider(local_dir, download_missing=True)
    monkeypatch.setattr(provider, "_load_from_explicit_path", lambda _path: None)
    monkeypatch.setattr(provider, "_load_cams_tif_group", lambda _date, _iso: None)
    monkeypatch.setattr(provider, "_has_cams_tif_candidates", lambda _date, _iso: False)
    monkeypatch.setattr(provider, "_load_cds_dataset", lambda _time: "cds")
    assert provider._load_cams_data(datetime(2024, 1, 1)) == "cds"


def test_cams_merge_helpers_cover_none_identity_and_interp_fallback(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    variable = xr.DataArray(
        np.full((2, 2), 2.0, dtype=np.float32),
        dims=["latitude", "longitude"],
        coords={"latitude": [1.5, 0.5], "longitude": [0.5, 1.5]},
    )
    assert CAMSProvider._align_missing_variable(None, variable) is variable

    reference = xr.DataArray(
        np.full((2, 2), 1.0, dtype=np.float32),
        dims=["latitude", "longitude"],
        coords={"latitude": [1.0, 0.0], "longitude": [0.0, 1.0]},
    )
    original_interp = xr.DataArray.interp

    def fake_interp(self, coords=None, method=None, **kwargs):  # type: ignore[no-untyped-def]
        if method == "linear":
            raise ValueError("force nearest")
        return original_interp(self, coords=coords, method=method, **kwargs)

    monkeypatch.setattr(xr.DataArray, "interp", fake_interp)
    aligned = CAMSProvider._align_missing_variable(reference, variable)
    assert aligned.coords["latitude"].identical(reference.coords["latitude"])
    assert aligned.coords["longitude"].identical(reference.coords["longitude"])

    fallback = _dataset(("aod550", "tcwv"))
    merged_none, added_none = CAMSProvider._merge_missing_variables(None, fallback)
    assert merged_none.equals(fallback)
    assert added_none == ["aod550", "tcwv"]

    primary = _dataset(("aod550", "tcwv"))
    merged_same, added_same = CAMSProvider._merge_missing_variables(primary, _dataset(("aod550",)))
    assert merged_same.equals(primary)
    assert added_same == []


def test_complete_cams_dataset_handles_complete_and_empty_fallbacks(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    provider = CAMSProvider(tmp_path, download_missing=True)
    complete = _dataset(("aod550", "tcwv", "gtco3"))
    assert provider._complete_cams_dataset(complete, datetime(2024, 1, 1)).equals(complete)

    partial = _dataset(("aod550",))
    empty_fallback = xr.Dataset(
        {"other": xr.DataArray(np.ones((1, 1), dtype=np.float32), dims=["latitude", "longitude"])}
    )
    cds = _dataset(("tcwv", "gtco3"))
    calls: list[str] = []

    def fake_remote(base_url: str, _time: datetime) -> xr.Dataset | None:
        calls.append(base_url)
        return empty_fallback

    monkeypatch.setattr(provider, "_load_from_remote_base", fake_remote)
    monkeypatch.setattr(provider, "_load_cds_dataset", lambda _time: cds)

    with caplog.at_level("INFO"):
        merged = provider._complete_cams_dataset(partial, datetime(2024, 1, 1))

    assert merged is not None
    assert set(merged.data_vars) == {"aod550", "tcwv", "gtco3"}
    assert "Supplemented CAMS variables from CDS download" in caplog.text
    assert calls == [provider._JASMIN_CAMS_BASE_URL]


def test_select_longitude_window_and_standardize_temporal_dim_edge_cases(tmp_path: Path) -> None:
    provider = CAMSProvider(tmp_path)
    var = xr.DataArray(
        np.arange(5, dtype=np.float32).reshape(1, 5),
        dims=["latitude", "longitude"],
        coords={"latitude": [0.0], "longitude": [0.0, 5.0, 10.0, 15.0, 20.0]},
    )
    only_right = provider._select_longitude_window(var, xmin=25.0, xmax=10.0)
    only_left = provider._select_longitude_window(var, xmin=10.0, xmax=-5.0)
    assert np.array_equal(only_right.coords["longitude"].values, np.array([0.0, 5.0, 10.0]))
    assert np.array_equal(only_left.coords["longitude"].values, np.array([10.0, 15.0, 20.0]))

    descending = xr.DataArray(
        np.arange(5, dtype=np.float32).reshape(1, 5),
        dims=["latitude", "longitude"],
        coords={"latitude": [0.0], "longitude": [20.0, 15.0, 10.0, 5.0, 0.0]},
    )
    descending_wrap = provider._select_longitude_window(descending, xmin=25.0, xmax=10.0)
    descending_left_only = provider._select_longitude_window(descending, xmin=10.0, xmax=-5.0)
    assert np.array_equal(descending_wrap.coords["longitude"].values, np.array([10.0, 5.0, 0.0]))
    assert np.array_equal(
        descending_left_only.coords["longitude"].values, np.array([20.0, 15.0, 10.0])
    )

    valid_time_2d = xr.DataArray(
        np.array([[np.datetime64("2024-01-01T00:00:00"), np.datetime64("2024-01-01T03:00:00")]]),
        dims=["forecast_reference_time", "forecast_period"],
    )
    with_valid_time = xr.DataArray(
        np.array([[[1.0], [2.0]]], dtype=np.float32),
        dims=["forecast_reference_time", "forecast_period", "longitude"],
        coords={
            "forecast_reference_time": [np.datetime64("2024-01-01T00:00:00")],
            "forecast_period": np.array([0, 3], dtype="timedelta64[h]"),
            "valid_time": valid_time_2d,
            "longitude": [0.0],
        },
    )
    standardized = provider._standardize_temporal_dims(with_valid_time, datetime(2024, 1, 1, 2))
    assert "time" in standardized.dims
    assert "forecast_period" not in standardized.coords

    without_valid_time = xr.DataArray(
        np.ones((2, 1), dtype=np.float32),
        dims=["forecast_period", "longitude"],
        coords={"forecast_period": np.array([0, 3], dtype="timedelta64[h]"), "longitude": [0.0]},
    )
    assert provider._standardize_temporal_dims(without_valid_time, datetime(2024, 1, 1)).identical(
        without_valid_time
    )

    mismatch_valid_time = xr.DataArray(
        np.ones((2, 2, 1), dtype=np.float32),
        dims=["forecast_reference_time", "forecast_period", "longitude"],
        coords={
            "forecast_reference_time": [
                np.datetime64("2024-01-01T00:00:00"),
                np.datetime64("2024-01-01T06:00:00"),
            ],
            "forecast_period": np.array([0, 3], dtype="timedelta64[h]"),
            "valid_time": (
                ("forecast_reference_time", "forecast_period"),
                np.array(
                    [
                        [
                            np.datetime64("2024-01-01T00:00:00"),
                            np.datetime64("2024-01-01T03:00:00"),
                        ],
                        [
                            np.datetime64("2024-01-01T06:00:00"),
                            np.datetime64("2024-01-01T09:00:00"),
                        ],
                    ],
                ),
            ),
            "longitude": [0.0],
        },
    )
    selected = provider._standardize_temporal_dims(mismatch_valid_time, datetime(2024, 1, 1, 5))
    assert "time" in selected.dims

    bad_valid_time = xr.DataArray(
        np.ones((2, 1), dtype=np.float32),
        dims=["forecast_period", "longitude"],
        coords={
            "forecast_period": np.array([0, 3], dtype="timedelta64[h]"),
            "valid_time": (
                ("forecast_period", "longitude"),
                np.array(
                    [[np.datetime64("2024-01-01T00:00:00")], [np.datetime64("2024-01-01T03:00:00")]]
                ),
            ),
            "longitude": [0.0],
        },
    )
    assert provider._standardize_temporal_dims(bad_valid_time, datetime(2024, 1, 1)).identical(
        bad_valid_time
    )


def test_normalize_source_and_storage_context_helpers(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    assert CAMSProvider._normalize_data_source("  cams_local  ") == Path("cams_local").expanduser()
    assert CAMSProvider._remote_scheme(tmp_path / "x") == ""
    assert CAMSProvider._normalize_data_source("s3://eodata/CAMS").startswith("s3://")

    tif_dir = tmp_path / "tifs"
    tif_dir.mkdir()
    (tif_dir / "cams_20240101_aod550.tif").write_text("x")
    provider = CAMSProvider(tif_dir)
    assert provider._has_cams_tif_candidates("20240101", "2024-01-01") is True

    remote = CAMSProvider("s3://eodata/CAMS/GLOBAL")
    with pytest.raises(TypeError, match="remote CAMS source"):
        remote._require_local_data_dir()

    assert remote._load_from_explicit_path("https://example.com/cams.invalid") is None

    local = CAMSProvider(tmp_path)

    def fake_load_local(path: Path, source_name: str | None = None) -> Path:
        del source_name
        return path

    monkeypatch.setattr(local, "_load_from_local_explicit_path", fake_load_local)
    resolved = local._load_from_explicit_path(tmp_path / "file.nc")
    assert resolved == tmp_path / "file.nc"

    manager = CredentialManager()
    manager.set_credentials("cdse", key="user", secret="secret")
    cdse_provider = CAMSProvider("s3://eodata/CAMS/GLOBAL", auth=manager)

    @contextmanager
    def fake_temp_creds():
        yield CDSES3Credentials("AK", "SK")

    monkeypatch.setattr(manager.cdse(), "temporary_s3_credentials", fake_temp_creds)
    with cdse_provider._remote_storage_options_context("s3://eodata/CAMS/GLOBAL/file.nc") as opts:
        assert opts["key"] == "AK"

    aws_manager = CredentialManager()
    aws_manager.set_credentials("aws", key="AWS_KEY", secret="AWS_SECRET")
    aws_provider = CAMSProvider("s3://bucket/path", auth=aws_manager)
    with aws_provider._remote_storage_options_context("s3://bucket/path/file.nc") as opts:
        assert opts == {"key": "AWS_KEY", "secret": "AWS_SECRET"}

    no_auth_provider = CAMSProvider("https://example.com/cams/")
    with no_auth_provider._remote_storage_options_context(
        "https://example.com/cams/file.nc"
    ) as opts:
        assert opts == {}


def test_extract_variable_adds_halo_before_reprojection(tmp_path: Path) -> None:
    provider = CAMSProvider(tmp_path)
    data = xr.Dataset(
        {
            "tcwv": xr.DataArray(
                np.arange(25, dtype=np.float32).reshape(1, 5, 5),
                dims=["time", "latitude", "longitude"],
                coords={
                    "time": [np.datetime64("2026-04-05T03:00:00")],
                    "latitude": [35.6, 35.2, 34.8, 34.4, 34.0],
                    "longitude": [118.0, 118.4, 118.8, 119.2, 119.6],
                },
            )
        }
    )
    bounds = (600000.0, 3790200.0, 709800.0, 3900000.0)
    crs = "EPSG:32650"
    obs_time = datetime(2026, 4, 5, 3, 0, 0)

    extracted = provider._extract_variable(data, "tcwv", bounds, crs, 60.0, obs_time)

    assert extracted.shape == (5, 5)
    assert np.array_equal(
        extracted.coords["latitude"].values, np.array([35.6, 35.2, 34.8, 34.4, 34.0])
    )
    assert np.array_equal(
        extracted.coords["longitude"].values, np.array([118.0, 118.4, 118.8, 119.2, 119.6])
    )


def test_remote_cache_and_loading_helpers_cover_error_paths(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    provider = CAMSProvider("s3://eodata/CAMS/GLOBAL", cache_dir=tmp_path / "cache")

    class _FakeFS:
        def __init__(self) -> None:
            self.calls: list[tuple[str, str]] = []

        def get(self, rpath: str, lpath: str) -> None:
            self.calls.append((rpath, lpath))
            Path(lpath).write_text("cached")

    fake_fs = _FakeFS()

    def fake_filesystem(protocol: str, **storage_options: object) -> _FakeFS:
        del protocol, storage_options
        return fake_fs

    def fake_open_local(urlpath: str, **kwargs: object) -> str:
        del urlpath, kwargs
        return str(tmp_path / "http_cached.nc")

    monkeypatch.setattr(
        provider,
        "_import_fsspec",
        lambda: SimpleNamespace(
            filesystem=fake_filesystem,
            open_local=fake_open_local,
        ),
    )

    s3_cached = provider._cache_remote_file(
        "s3://eodata/CAMS/GLOBAL/2024/01/01/file.nc", storage_options={"anon": True}
    )
    assert s3_cached.exists()
    assert fake_fs.calls[0][0] == "eodata/CAMS/GLOBAL/2024/01/01/file.nc"

    http_provider = CAMSProvider("https://example.com/cams/", cache_dir=tmp_path / "cache")
    monkeypatch.setattr(
        http_provider,
        "_import_fsspec",
        lambda: SimpleNamespace(
            filesystem=fake_filesystem,
            open_local=fake_open_local,
        ),
    )
    http_cached = http_provider._cache_remote_file("https://example.com/cams/file.nc")
    assert http_cached == tmp_path / "http_cached.nc"

    def missing_cache(url: str, storage_options: dict | None = None) -> Path:
        del storage_options
        raise FileNotFoundError(url)

    monkeypatch.setattr(provider, "_cache_remote_file", missing_cache)
    with caplog.at_level("WARNING"):
        assert (
            provider._cache_remote_path_with_options(
                "s3://eodata/CAMS/GLOBAL/missing.nc",
                missing_ok=False,
                storage_options={},
            )
            is None
        )
    assert "not found" in caplog.text or "Failed to cache" in caplog.text

    assert (
        provider._cache_remote_path_with_options(
            "s3://eodata/CAMS/GLOBAL/missing.nc", missing_ok=True, storage_options={}
        )
        is None
    )

    assert provider._load_from_remote_url("s3://eodata/CAMS/GLOBAL/file.txt") is None

    def resolve_none(*args, **kwargs):  # type: ignore[no-untyped-def]
        del args, kwargs
        return None

    monkeypatch.setattr(provider, "_resolve_remote_local_path", resolve_none)
    assert (
        provider._load_from_remote_url("s3://eodata/CAMS/GLOBAL/file.nc", missing_ok=True) is None
    )

    assert provider._load_from_remote_s3_base("s3://bucket/not-cams", datetime(2024, 1, 1)) is None
    monkeypatch.setattr(provider, "_select_cdse_cams_files", lambda _base, _time, _opts: [])
    assert (
        provider._load_from_remote_s3_base("s3://eodata/CAMS/GLOBAL", datetime(2024, 1, 1)) is None
    )

    monkeypatch.setattr(
        provider,
        "_select_cdse_cams_files",
        lambda _base, _time, _opts: ["s3://eodata/CAMS/GLOBAL/file.nc"],
    )

    def fake_load_remote_url(
        url: str, missing_ok: bool = False, storage_options: dict | None = None
    ) -> None:
        del url, missing_ok, storage_options
        return None

    monkeypatch.setattr(provider, "_load_from_remote_url", fake_load_remote_url)
    assert (
        provider._load_from_remote_s3_base("s3://eodata/CAMS/GLOBAL", datetime(2024, 1, 1)) is None
    )


def test_select_cdse_files_tif_dataset_and_download_helpers(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    import fsspec

    provider = CAMSProvider("s3://eodata/CAMS/GLOBAL", cache_dir=tmp_path / "cache")

    class _MissingFS:
        def ls(self, path, detail=False):  # noqa: ARG002
            raise FileNotFoundError(path)

    monkeypatch.setattr(fsspec, "filesystem", lambda _protocol, **_kwargs: _MissingFS())
    assert (
        provider._select_cdse_cams_files("s3://eodata/CAMS/GLOBAL", datetime(2024, 1, 1), {}) == []
    )

    class _ChoiceFS:
        def ls(self, path, detail=False):  # noqa: ARG002
            return [
                "eodata/CAMS/GLOBAL/2024/01/01/not_a_match",
                "eodata/CAMS/GLOBAL/2024/01/01/z_cams_c_ecmf_20240101000000_prod_an_sfc_000_unknown",
                "eodata/CAMS/GLOBAL/2024/01/01/z_cams_c_ecmf_20240101000000_prod_fc_sfc_003_aod550",
                "eodata/CAMS/GLOBAL/2024/01/01/z_cams_c_ecmf_20240101000000_prod_an_sfc_000_aod550",
            ]

    monkeypatch.setattr(fsspec, "filesystem", lambda _protocol, **_kwargs: _ChoiceFS())
    urls = provider._select_cdse_cams_files("s3://eodata/CAMS/GLOBAL", datetime(2024, 1, 1, 1), {})
    assert len(urls) == 1
    assert "aod550" in urls[0]

    tif_provider = CAMSProvider(tmp_path)
    da = xr.DataArray(
        np.ones((3, 1, 1), dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": [1, 2, 3]},
        attrs={"long_name": ["aod550", "aod550", "gtco3"]},
    )
    ds = tif_provider._dataset_from_multiband_tif(da)
    assert set(ds.data_vars) == {"aod550", "gtco3"}

    file_target = tmp_path / "cams-file.nc"
    file_target.write_text("x")
    file_provider = CAMSProvider(file_target)

    def unreachable_client(**kwargs):  # type: ignore[no-untyped-def]
        del kwargs
        raise RuntimeError("should not be called")

    monkeypatch.setattr(
        file_provider, "_import_cdsapi", lambda: SimpleNamespace(Client=unreachable_client)
    )
    monkeypatch.setattr(file_provider, "_auth", CredentialManager())
    # REVIEW.md §2.1, §3.3 cams.py:1058-1063: cdsapi failures are now
    # raised as ``DataNotFoundError`` so auth/quota failures don't
    # silently masquerade as "data unavailable". The ``RuntimeError``
    # raised by the fake client is wrapped and re-raised.
    from siac.errors import DataNotFoundError

    with pytest.raises(DataNotFoundError, match="should not be called"):
        file_provider._download_cams_file(datetime(2024, 1, 1))

    auth = CredentialManager()
    auth.set_credentials("cds", key="api-key")
    downloader = CAMSProvider(file_target, auth=auth)
    captured: dict[str, object] = {}

    class _Request:
        def download(self, path: str) -> str:
            captured["path"] = path
            return path

    class _Client:
        def __init__(self, **kwargs):
            captured["kwargs"] = kwargs

        def retrieve(self, dataset: str, request: dict[str, object]) -> _Request:
            captured["dataset"] = dataset
            captured["request"] = request
            return _Request()

    monkeypatch.setattr(downloader, "_import_cdsapi", lambda: SimpleNamespace(Client=_Client))
    out = downloader._download_cams_file(datetime(2024, 1, 2))
    assert out == file_target.parent / "CAMS_2024-01-02.nc"
    assert captured["dataset"] == downloader._CDS_DATASET
