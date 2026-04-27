from __future__ import annotations

from datetime import datetime
from typing import TYPE_CHECKING

import numpy as np
import pytest
import xarray as xr

import siac.adapters.earthdata as earthdata_mod
import siac.adapters.output as output_mod
import siac.algorithms.cloud.providers.omnicloudmask as omnicloudmask_mod
from siac.algorithms.cloud.providers.omnicloudmask import OmniCloudMaskProvider
from siac.config.schema import OutputDefaultsConfig
from siac.runtime import CorrectionResult

if TYPE_CHECKING:
    from types import ModuleType


def _correction_result(
    *,
    include_uncertainty: bool,
    include_rgb: bool,
) -> CorrectionResult:
    coords = {"y": [0, 1], "x": [0, 1]}
    boa_vars: dict[str, xr.DataArray] = {
        "B03": xr.DataArray(
            np.full((2, 2), 0.10, dtype=np.float32), dims=["y", "x"], coords=coords
        ),
        "B02": xr.DataArray(
            np.full((2, 2), 0.08, dtype=np.float32), dims=["y", "x"], coords=coords
        ),
    }
    if include_rgb:
        boa_vars["B04"] = xr.DataArray(
            np.full((2, 2), 0.12, dtype=np.float32),
            dims=["y", "x"],
            coords=coords,
        )

    boa = xr.Dataset(boa_vars)
    boa_unc = (
        xr.Dataset(
            {
                name: xr.DataArray(
                    np.full((2, 2), 0.01, dtype=np.float32),
                    dims=["y", "x"],
                    coords=coords,
                )
                for name in boa.data_vars
            }
        )
        if include_uncertainty
        else None
    )
    return CorrectionResult(
        boa=boa,
        boa_unc=boa_unc,
        aot=xr.DataArray(np.full((2, 2), 0.2, dtype=np.float32), dims=["y", "x"], coords=coords),
        tcwv=xr.DataArray(np.full((2, 2), 2.0, dtype=np.float32), dims=["y", "x"], coords=coords),
        cloud_mask=xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"], coords=coords),
    )


def test_earthdata_helpers_cover_selection_grid_and_merge_paths(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    p0 = tmp_path / "p0.hdf"
    p1 = tmp_path / "p1.hdf"
    p2 = tmp_path / "p2.hdf"

    assert (
        earthdata_mod.select_candidate_paths(
            [],
            obs_time=datetime(2024, 1, 2),
            bounds=(0.0, 0.0, 1.0, 1.0),
            crs="EPSG:4326",
        )
        == []
    )

    monkeypatch.setattr(
        earthdata_mod,
        "parse_granule_date",
        lambda path: (
            (_ for _ in ()).throw(RuntimeError("bad metadata"))
            if path == p0
            else datetime(2024, 1, 2, 12, 0, 0)
        ),
    )
    assert earthdata_mod.select_candidate_paths(
        [p0, p1],
        obs_time=datetime(2024, 1, 2),
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
    ) == [p0, p1]

    monkeypatch.setattr(
        earthdata_mod,
        "parse_granule_date",
        lambda path: {
            p1: datetime(2024, 1, 2, 8, 0, 0),
            p2: datetime(2024, 1, 2, 10, 0, 0),
        }[path],
    )
    monkeypatch.setattr(
        earthdata_mod,
        "granule_intersects_bounds",
        lambda path, **_kwargs: path != p1,
    )
    monkeypatch.setattr(
        earthdata_mod,
        "parse_tile_indices",
        lambda path: (_ for _ in ()).throw(RuntimeError("missing tile")) if path == p2 else (1, 1),
    )

    selected = earthdata_mod.select_candidate_paths(
        [p1, p2],
        obs_time=datetime(2024, 1, 2, 9, 0, 0),
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        sample_dates=np.array(["2024-01-02"], dtype="datetime64[D]"),
    )
    assert selected == [p2]

    none_selected = earthdata_mod.select_candidate_paths(
        [p1],
        obs_time=datetime(2024, 1, 2, 9, 0, 0),
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        sample_dates=np.array(["2024-01-03"], dtype="datetime64[D]"),
    )
    assert none_selected == []

    with pytest.raises(ValueError, match="target resolution must be > 0"):
        earthdata_mod.target_grid_coords(
            (0.0, 0.0, 1.0, 1.0), 0.0, resolution_name="target resolution"
        )

    grid = earthdata_mod.constant_target_array((0.0, 0.0, 2.0, 2.0), 1.0, 3.5)
    assert grid.shape == (2, 2)
    assert float(grid.values[0, 0]) == pytest.approx(3.5)

    band_grid = earthdata_mod.constant_target_band_array(
        ["B02", "B03"], (0.0, 0.0, 2.0, 2.0), 1.0, 7.0
    )
    assert band_grid.dims == ("band", "y", "x")
    assert band_grid.coords["band"].values.tolist() == ["B02", "B03"]

    with pytest.raises(ValueError, match="at least one array"):
        earthdata_mod.merge_reprojected_tiles(
            [],
            bounds=(0.0, 0.0, 1.0, 1.0),
            crs="EPSG:4326",
            resolution=1.0,
            resampling="nearest",  # type: ignore[arg-type]
            nodata=None,
        )

    arrays = [
        xr.DataArray(np.array([[1.0, np.nan]], dtype=np.float32), dims=["y", "x"]),
        xr.DataArray(np.array([[np.nan, 2.0]], dtype=np.float32), dims=["y", "x"]),
    ]
    monkeypatch.setattr(earthdata_mod, "reproject_native_to_target", lambda arr, **_kwargs: arr)
    merged = earthdata_mod.merge_reprojected_tiles(
        arrays,
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        resolution=1.0,
        resampling="nearest",  # type: ignore[arg-type]
        nodata=None,
    )
    np.testing.assert_allclose(merged.values, np.array([[1.0, 2.0]], dtype=np.float32))
    assert merged.dtype == np.float32


def test_omnicloudmask_helpers_cover_import_and_normalization_paths(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _raise_import(_name: str) -> ModuleType:
        raise ModuleNotFoundError("missing")

    monkeypatch.setattr(omnicloudmask_mod, "import_module", _raise_import)
    with pytest.raises(ImportError, match="omnicloudmask is required"):
        OmniCloudMaskProvider()._default_predictor()

    template = xr.DataArray(
        np.zeros((2, 2), dtype=np.float32), dims=["y", "x"], coords={"y": [0, 1], "x": [0, 1]}
    )
    single_channel = OmniCloudMaskProvider._normalize_raw_output(
        np.array([[[1, 2], [3, 4]]], dtype=np.uint8),
        template,
    )
    np.testing.assert_array_equal(single_channel.values, np.array([[1, 2], [3, 4]], dtype=np.uint8))

    channel_last = OmniCloudMaskProvider._normalize_raw_output(
        np.array(
            [
                [[0.1, 0.9], [0.8, 0.2]],
                [[0.7, 0.3], [0.2, 0.8]],
            ],
            dtype=np.float32,
        ),
        template,
    )
    np.testing.assert_array_equal(channel_last.values, np.array([[1, 0], [0, 1]]))


def test_omnicloudmask_predict_honors_custom_mapping_and_missing_mask() -> None:
    red = xr.DataArray(
        np.array([[0.1, np.nan]], dtype=np.float32), dims=["y", "x"], coords={"y": [0], "x": [0, 1]}
    )
    green = xr.DataArray(
        np.array([[0.2, 0.3]], dtype=np.float32), dims=["y", "x"], coords={"y": [0], "x": [0, 1]}
    )
    nir = xr.DataArray(
        np.array([[0.4, 0.5]], dtype=np.float32), dims=["y", "x"], coords={"y": [0], "x": [0, 1]}
    )

    provider = OmniCloudMaskProvider(predictor=lambda _rgbnir: np.array([[0, 1]], dtype=np.uint8))
    out = provider.predict(
        red,
        green,
        nir,
        class_mapping={1: [0], 3: [1]},
        unmapped_to_missing=False,
    )

    assert out.name == "cloud_classes"
    assert out.dtype == np.uint8
    np.testing.assert_array_equal(out.values, np.array([[1, 0]], dtype=np.uint8))


@pytest.mark.parametrize(("fmt", "ext"), [("netcdf", ".nc"), ("zarr", ".zarr")])
def test_output_writer_array_formats_cover_uncertainty_aux_and_quicklook(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
    fmt: str,
    ext: str,
) -> None:
    calls: list[tuple[str, xr.Dataset, object]] = []

    def _record_writer(dataset: xr.Dataset, path):
        calls.append((path.name, dataset, path))
        return path

    def _record_quicklook(_dataset: xr.Dataset, path):
        return path

    monkeypatch.setattr(output_mod, "write_netcdf", _record_writer)
    monkeypatch.setattr(output_mod, "write_zarr", _record_writer)
    monkeypatch.setattr(output_mod, "write_rgb_quicklook", _record_quicklook)

    writer = output_mod.ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format=fmt,
            include_uncertainty=True,
            include_auxiliary=True,
            include_rgb=True,
        )
    )
    artifacts = writer.write(
        _correction_result(include_uncertainty=True, include_rgb=True), tmp_path
    )

    prefix = "SAT_L2A_00000000T000000"
    assert artifacts["boa"] == tmp_path / f"{prefix}_BOA{ext}"
    assert "boa_unc" in artifacts
    assert "auxiliary" in artifacts
    assert artifacts["quicklook.rgb"] == tmp_path / f"{prefix}_RGB.tif"
    written_names = {name for name, _, _ in calls}
    assert f"{prefix}_BOA{ext}" in written_names
    assert f"{prefix}_BOA_UNC{ext}" in written_names
    assert f"{prefix}_AUX{ext}" in written_names


def test_output_writer_raster_skips_optional_products_when_absent(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    raster_calls: list[str] = []

    def _fake_write_raster(data: xr.DataArray, path, **kwargs):
        raster_calls.append(path.name)
        return path

    def _unexpected_quicklook(_dataset: xr.Dataset, _path):
        raise AssertionError("should not write quicklook")

    monkeypatch.setattr(output_mod, "write_raster", _fake_write_raster)
    monkeypatch.setattr(output_mod, "write_rgb_quicklook", _unexpected_quicklook)

    writer = output_mod.ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="geotiff",
            include_uncertainty=True,
            include_auxiliary=False,
            include_rgb=False,
        )
    )
    result = _correction_result(include_uncertainty=False, include_rgb=False)
    artifacts = writer.write(result, tmp_path)

    # Should have BOA bands but no uncertainty (result has none), no quicklook
    assert any(k.startswith("boa.") for k in artifacts)
    assert "boa_unc.B03" not in artifacts
    assert "quicklook.rgb" not in artifacts


def test_output_writer_rgb_helper_returns_none_when_disabled_or_missing_bands(tmp_path) -> None:
    prefix = "SAT_L2A_00000000T000000"
    defaults = OutputDefaultsConfig(include_rgb=False)
    writer = output_mod.ConfiguredOutputWriter(defaults)
    assert (
        writer._write_rgb_if_available(
            _correction_result(include_uncertainty=False, include_rgb=True), tmp_path, prefix
        )
        is None
    )

    defaults = OutputDefaultsConfig(include_rgb=True)
    writer = output_mod.ConfiguredOutputWriter(defaults)
    assert (
        writer._write_rgb_if_available(
            _correction_result(include_uncertainty=False, include_rgb=False), tmp_path, prefix
        )
        is None
    )
