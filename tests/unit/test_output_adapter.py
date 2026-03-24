from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
import xarray as xr

import siac.adapters.output as output_module
from siac.adapters.output import ConfiguredOutputWriter
from siac.config.schema import OutputDefaultsConfig
from siac.runtime import CorrectionResult

if TYPE_CHECKING:
    from pathlib import Path


def _result(*, include_uncertainty: bool = True) -> CorrectionResult:
    coords = {"y": [0, 1], "x": [0, 1]}
    boa = xr.Dataset(
        {
            "B04": xr.DataArray(np.full((2, 2), 0.12, dtype=np.float32), dims=["y", "x"], coords=coords),
            "B03": xr.DataArray(np.full((2, 2), 0.10, dtype=np.float32), dims=["y", "x"], coords=coords),
            "B02": xr.DataArray(np.full((2, 2), 0.08, dtype=np.float32), dims=["y", "x"], coords=coords),
        }
    )
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


def test_write_raster_products_emits_boa_unc_aux_and_quicklook(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset_calls: list[tuple[xr.Dataset, Path, dict[str, object]]] = []
    raster_calls: list[tuple[Path, dict[str, object], np.dtype[np.generic]]] = []
    quicklook_calls: list[Path] = []

    def _fake_write_dataset(dataset: xr.Dataset, output_dir: Path, **kwargs: object) -> dict[str, Path]:
        dataset_calls.append((dataset, output_dir, kwargs))
        return {name: output_dir / f"{name}.tif" for name in dataset.data_vars}

    def _fake_write_cog(data: xr.DataArray, path: Path, **kwargs: object) -> Path:
        raster_calls.append((path, kwargs, data.dtype))
        return path

    def _fake_write_rgb_quicklook(dataset: xr.Dataset, path: Path) -> Path:
        quicklook_calls.append(path)
        return path

    monkeypatch.setattr(output_module, "write_dataset", _fake_write_dataset)
    monkeypatch.setattr(output_module, "write_cog", _fake_write_cog)
    monkeypatch.setattr(output_module, "write_rgb_quicklook", _fake_write_rgb_quicklook)

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="cog",
            boa_dtype="uint16",
            include_uncertainty=True,
            include_auxiliary=True,
            include_rgb=True,
        )
    )
    artifacts = writer.write(_result(), tmp_path)

    assert len(dataset_calls) == 2
    assert dataset_calls[0][1] == tmp_path / "boa"
    assert dataset_calls[0][2]["as_cog"] is True
    assert dataset_calls[0][2]["dtype"] == "uint16"
    assert dataset_calls[0][2]["nodata"] == 0
    assert all(str(dataset[name].dtype) == "uint16" for dataset, _, _ in dataset_calls for name in dataset.data_vars)

    assert len(raster_calls) == 3
    assert {path.name for path, _, _ in raster_calls} == {"aot.tif", "tcwv.tif", "cloud_mask.tif"}
    cloud_mask_call = next(call for call in raster_calls if call[0].name == "cloud_mask.tif")
    assert cloud_mask_call[1]["compression"] == "lzw"
    assert cloud_mask_call[1]["dtype"] == "uint8"
    assert cloud_mask_call[1]["nodata"] == 255

    assert quicklook_calls == [tmp_path / "quicklook.tif"]
    assert {
        "boa.B04",
        "boa.B03",
        "boa.B02",
        "boa_unc.B04",
        "boa_unc.B03",
        "boa_unc.B02",
        "auxiliary.aot",
        "auxiliary.tcwv",
        "auxiliary.cloud_mask",
        "quicklook.rgb",
    } <= set(artifacts)


def test_write_netcdf_alias_emits_cog_band_products_and_respects_toggles(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    dataset_calls: list[tuple[xr.Dataset, Path, dict[str, object]]] = []

    def _fake_write_dataset(dataset: xr.Dataset, output_dir: Path, **kwargs: object) -> dict[str, Path]:
        dataset_calls.append((dataset, output_dir, kwargs))
        return {name: output_dir / f"{name}.tif" for name in dataset.data_vars}

    monkeypatch.setattr(output_module, "write_dataset", _fake_write_dataset)

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="netcdf",
            include_uncertainty=False,
            include_auxiliary=False,
            include_rgb=False,
        )
    )
    artifacts = writer.write(_result(include_uncertainty=False), tmp_path)

    assert len(dataset_calls) == 1
    assert dataset_calls[0][1] == tmp_path / "boa"
    assert dataset_calls[0][2]["as_cog"] is True
    assert artifacts == {
        "boa.B04": tmp_path / "boa" / "B04.tif",
        "boa.B03": tmp_path / "boa" / "B03.tif",
        "boa.B02": tmp_path / "boa" / "B02.tif",
    }
    assert "treating it as 'cog'" in caplog.text


def test_write_zarr_products_route_and_respects_toggles(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[Path] = []

    def _fake_writer(dataset: xr.Dataset, path: Path) -> Path:
        calls.append(path)
        return path

    monkeypatch.setattr(output_module, "write_zarr", _fake_writer)

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="zarr",
            include_uncertainty=False,
            include_auxiliary=False,
            include_rgb=False,
        )
    )
    artifacts = writer.write(_result(include_uncertainty=False), tmp_path)

    assert calls == [tmp_path / "boa.zarr"]
    assert artifacts == {"boa": tmp_path / "boa.zarr"}


def test_write_rejects_unsupported_format(tmp_path: Path) -> None:
    defaults = OutputDefaultsConfig.model_construct(
        output_dir=None,
        format="bad-format",
        compression="deflate",
        include_rgb=True,
        include_uncertainty=True,
        include_auxiliary=True,
        boa_dtype="float32",
        boa_scale=10000.0,
        boa_nodata=0.0,
    )

    with pytest.raises(ValueError, match="Unsupported output format"):
        ConfiguredOutputWriter(defaults).write(_result(), tmp_path)
