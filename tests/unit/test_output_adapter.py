from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
import xarray as xr

import siac.adapters.output as output_module
from siac.adapters.output import ConfiguredOutputWriter
from siac.config.schema import OutputDefaultsConfig
from siac.runtime import (
    AOTScatterBandDiagnostics,
    CorrectionDiagnostics,
    CorrectionResult,
    MonthlyCompositeOutput,
)

if TYPE_CHECKING:
    from pathlib import Path


def _result(
    *,
    include_uncertainty: bool = True,
    include_surface_prior: bool = False,
    include_monthly_composites: bool = False,
    include_diagnostics: bool = False,
    include_solver_qa: bool = False,
) -> CorrectionResult:
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
    surface_prior = (
        xr.Dataset(
            {
                "B04": xr.DataArray(np.full((2, 2), 0.11, dtype=np.float32), dims=["y", "x"], coords=coords),
                "B03": xr.DataArray(np.full((2, 2), 0.09, dtype=np.float32), dims=["y", "x"], coords=coords),
            }
        )
        if include_surface_prior
        else None
    )
    surface_prior_unc = (
        xr.Dataset(
            {
                name: xr.DataArray(
                    np.full((2, 2), 0.02, dtype=np.float32),
                    dims=["y", "x"],
                    coords=coords,
                )
                for name in ("B04", "B03")
            }
        )
        if include_surface_prior and include_uncertainty
        else None
    )
    monthly_composites = (
        {
            "2023_07": MonthlyCompositeOutput(
                reflectance=xr.Dataset(
                    {
                        "B04": xr.DataArray(np.full((2, 2), 0.13, dtype=np.float32), dims=["y", "x"], coords=coords),
                        "B03": xr.DataArray(np.full((2, 2), 0.12, dtype=np.float32), dims=["y", "x"], coords=coords),
                    }
                ),
                quality=xr.DataArray(np.full((2, 2), 0.4, dtype=np.float32), dims=["y", "x"], coords=coords),
                sample_index=xr.DataArray(np.full((2, 2), 3, dtype=np.int16), dims=["y", "x"], coords=coords),
            )
        }
        if include_monthly_composites
        else None
    )
    solver_qa = (
        xr.Dataset(
            {
                "invalid_retrieval": xr.DataArray(
                    np.array([[False, True], [False, False]], dtype=bool),
                    dims=["y", "x"],
                    coords=coords,
                ),
                "low_quality": xr.DataArray(
                    np.array([[True, True], [False, False]], dtype=bool),
                    dims=["y", "x"],
                    coords=coords,
                ),
            }
        )
        if include_solver_qa
        else None
    )
    return CorrectionResult(
        boa=boa,
        boa_unc=boa_unc,
        aot=xr.DataArray(np.full((2, 2), 0.2, dtype=np.float32), dims=["y", "x"], coords=coords),
        tcwv=xr.DataArray(np.full((2, 2), 2.0, dtype=np.float32), dims=["y", "x"], coords=coords),
        cloud_mask=xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"], coords=coords),
        surface_prior=surface_prior,
        surface_prior_unc=surface_prior_unc,
        solver_qa=solver_qa,
        monthly_composites=monthly_composites,
        diagnostics=CorrectionDiagnostics(
            processing_time_s=0.25,
            aot_scatter_plots=(
                AOTScatterBandDiagnostics(
                    band_name="B02",
                    surface_reflectance=np.array([0.1, 0.2], dtype=np.float32),
                    observed_toa=np.array([0.15, 0.25], dtype=np.float32),
                    simulated_toa=np.array([0.14, 0.24], dtype=np.float32),
                    total_valid_count=2,
                ),
            ) if include_diagnostics else (),
        ),
    )


def test_write_raster_products_emits_boa_unc_aux_and_quicklook(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset_calls: list[tuple[xr.Dataset, Path, dict[str, object]]] = []
    raster_calls: list[tuple[Path, dict[str, object], np.dtype[np.generic]]] = []
    quicklook_calls: list[Path] = []
    scatter_calls: list[Path] = []

    def _fake_write_dataset(dataset: xr.Dataset, output_dir: Path, **kwargs: object) -> dict[str, Path]:
        dataset_calls.append((dataset, output_dir, kwargs))
        return {name: output_dir / f"{name}.tif" for name in dataset.data_vars}

    def _fake_write_cog(data: xr.DataArray, path: Path, **kwargs: object) -> Path:
        raster_calls.append((path, kwargs, data.dtype))
        return path

    def _fake_write_rgb_quicklook(dataset: xr.Dataset, path: Path) -> Path:
        quicklook_calls.append(path)
        return path

    def _fake_write_aot_scatter_plot(_scatter: object, path: Path) -> Path:
        scatter_calls.append(path)
        return path

    monkeypatch.setattr(output_module, "write_dataset", _fake_write_dataset)
    monkeypatch.setattr(output_module, "write_cog", _fake_write_cog)
    monkeypatch.setattr(output_module, "write_rgb_quicklook", _fake_write_rgb_quicklook)
    monkeypatch.setattr(output_module, "write_aot_scatter_plot", _fake_write_aot_scatter_plot)

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="cog",
            boa_dtype="uint16",
            include_uncertainty=True,
            include_auxiliary=True,
            include_rgb=True,
        )
    )
    artifacts = writer.write(_result(include_diagnostics=True), tmp_path)

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
    assert scatter_calls == [tmp_path / "diagnostics" / "aot_scatter_B02.png"]
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
        "diagnostics.scatter.B02",
        "quicklook.rgb",
    } <= set(artifacts)


def test_write_raster_products_emits_solver_qa_masks(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    raster_calls: list[tuple[Path, dict[str, object], np.dtype[np.generic]]] = []

    monkeypatch.setattr(
        output_module,
        "write_dataset",
        lambda dataset, output_dir, **_kwargs: {name: output_dir / f"{name}.tif" for name in dataset.data_vars},
    )

    def _fake_write_cog(data: xr.DataArray, path: Path, **kwargs: object) -> Path:
        raster_calls.append((path, kwargs, data.dtype))
        return path

    monkeypatch.setattr(output_module, "write_cog", _fake_write_cog)
    monkeypatch.setattr(output_module, "write_rgb_quicklook", lambda _dataset, path: path)
    monkeypatch.setattr(output_module, "write_aot_scatter_plot", lambda _scatter, path: path)

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="cog",
            boa_dtype="uint16",
            include_uncertainty=False,
            include_auxiliary=True,
            include_rgb=False,
        )
    )

    artifacts = writer.write(_result(include_uncertainty=False, include_solver_qa=True), tmp_path)

    qa_calls = [call for call in raster_calls if call[0].name.startswith("qa_")]
    assert {path.name for path, _, _ in qa_calls} == {"qa_invalid_retrieval.tif", "qa_low_quality.tif"}
    assert all(kwargs["dtype"] == "uint8" for _, kwargs, _ in qa_calls)
    assert all(kwargs["nodata"] == 255 for _, kwargs, _ in qa_calls)
    assert {
        "auxiliary.qa.invalid_retrieval",
        "auxiliary.qa.low_quality",
    } <= set(artifacts)


def test_write_netcdf_products_route_and_respects_toggles(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[Path] = []

    def _fake_write_netcdf(dataset: xr.Dataset, path: Path, **kwargs: object) -> Path:
        del dataset, kwargs
        calls.append(path)
        return path

    monkeypatch.setattr(output_module, "write_netcdf", _fake_write_netcdf)
    monkeypatch.setattr(output_module, "write_aot_scatter_plot", lambda _scatter, path: path)

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="netcdf",
            include_uncertainty=False,
            include_auxiliary=False,
            include_rgb=False,
        )
    )
    artifacts = writer.write(_result(include_uncertainty=False), tmp_path)

    assert calls == [tmp_path / "boa.nc"]
    assert artifacts == {"boa": tmp_path / "boa.nc"}


def test_write_netcdf_auxiliary_aligns_cloud_mask_to_atmo_grid(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, xr.Dataset] = {}

    def _fake_write_netcdf(dataset: xr.Dataset, path: Path, **kwargs: object) -> Path:
        del kwargs
        captured[str(path)] = dataset
        return path

    monkeypatch.setattr(output_module, "write_netcdf", _fake_write_netcdf)

    coarse_coords = {"y": [1500.0, 500.0], "x": [500.0, 1500.0]}
    fine_coords = {"y": [1750.0, 1250.0, 750.0, 250.0], "x": [250.0, 750.0, 1250.0, 1750.0]}
    result = CorrectionResult(
        boa=xr.Dataset(
            {
                "B04": xr.DataArray(np.full((4, 4), 0.12, dtype=np.float32), dims=["y", "x"], coords=fine_coords),
                "B03": xr.DataArray(np.full((4, 4), 0.10, dtype=np.float32), dims=["y", "x"], coords=fine_coords),
                "B02": xr.DataArray(np.full((4, 4), 0.08, dtype=np.float32), dims=["y", "x"], coords=fine_coords),
            }
        ),
        boa_unc=None,
        aot=xr.DataArray(np.full((2, 2), 0.2, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
        tcwv=xr.DataArray(np.full((2, 2), 2.0, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
        cloud_mask=xr.DataArray(np.zeros((4, 4), dtype=bool), dims=["y", "x"], coords=fine_coords),
    )
    result.cloud_mask.values[0, 0] = True

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="netcdf",
            include_uncertainty=False,
            include_auxiliary=True,
            include_rgb=False,
        )
    )

    artifacts = writer.write(result, tmp_path)

    aux_ds = captured[str(tmp_path / "auxiliary.nc")]
    assert artifacts["auxiliary"] == tmp_path / "auxiliary.nc"
    assert aux_ds.sizes == {"y": 2, "x": 2}
    assert aux_ds["cloud_mask"].coords["x"].identical(result.aot.coords["x"])
    assert aux_ds["cloud_mask"].coords["y"].identical(result.aot.coords["y"])
    assert aux_ds["cloud_mask"].dtype == np.uint8


def test_write_netcdf_auxiliary_includes_solver_qa_masks(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, xr.Dataset] = {}

    def _fake_write_netcdf(dataset: xr.Dataset, path: Path, **kwargs: object) -> Path:
        del kwargs
        captured[str(path)] = dataset
        return path

    monkeypatch.setattr(output_module, "write_netcdf", _fake_write_netcdf)

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="netcdf",
            include_uncertainty=False,
            include_auxiliary=True,
            include_rgb=False,
        )
    )

    writer.write(_result(include_uncertainty=False, include_solver_qa=True), tmp_path)

    aux_ds = captured[str(tmp_path / "auxiliary.nc")]
    assert {"invalid_retrieval", "low_quality"} <= set(aux_ds.data_vars)
    assert aux_ds["invalid_retrieval"].dtype == np.uint8
    assert aux_ds["low_quality"].dtype == np.uint8


def test_write_zarr_auxiliary_aligns_cloud_mask_to_atmo_grid(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, xr.Dataset] = {}

    def _fake_write_zarr(dataset: xr.Dataset, path: Path) -> Path:
        captured[str(path)] = dataset
        return path

    monkeypatch.setattr(output_module, "write_zarr", _fake_write_zarr)

    coarse_coords = {"y": [1500.0, 500.0], "x": [500.0, 1500.0]}
    fine_coords = {"y": [1750.0, 1250.0, 750.0, 250.0], "x": [250.0, 750.0, 1250.0, 1750.0]}
    result = CorrectionResult(
        boa=xr.Dataset(
            {
                "B04": xr.DataArray(np.full((4, 4), 0.12, dtype=np.float32), dims=["y", "x"], coords=fine_coords),
                "B03": xr.DataArray(np.full((4, 4), 0.10, dtype=np.float32), dims=["y", "x"], coords=fine_coords),
                "B02": xr.DataArray(np.full((4, 4), 0.08, dtype=np.float32), dims=["y", "x"], coords=fine_coords),
            }
        ),
        boa_unc=None,
        aot=xr.DataArray(np.full((2, 2), 0.2, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
        tcwv=xr.DataArray(np.full((2, 2), 2.0, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
        cloud_mask=xr.DataArray(np.zeros((4, 4), dtype=bool), dims=["y", "x"], coords=fine_coords),
    )
    result.cloud_mask.values[0, 0] = True

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="zarr",
            include_uncertainty=False,
            include_auxiliary=True,
            include_rgb=False,
        )
    )

    artifacts = writer.write(result, tmp_path)

    aux_ds = captured[str(tmp_path / "auxiliary.zarr")]
    assert artifacts["auxiliary"] == tmp_path / "auxiliary.zarr"
    assert aux_ds.sizes == {"y": 2, "x": 2}
    assert aux_ds["cloud_mask"].coords["x"].identical(result.aot.coords["x"])
    assert aux_ds["cloud_mask"].coords["y"].identical(result.aot.coords["y"])
    assert aux_ds["cloud_mask"].dtype == np.uint8


def test_write_netcdf_auxiliary_reorders_same_shape_cloud_mask_coords(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, xr.Dataset] = {}

    def _fake_write_netcdf(dataset: xr.Dataset, path: Path, **kwargs: object) -> Path:
        del kwargs
        captured[str(path)] = dataset
        return path

    monkeypatch.setattr(output_module, "write_netcdf", _fake_write_netcdf)

    coarse_coords = {"y": [1500.0, 500.0], "x": [500.0, 1500.0]}
    shifted_mask = xr.DataArray(
        np.array([[0, 1], [0, 0]], dtype=bool),
        dims=["y", "x"],
        coords={"y": [500.0, 1500.0], "x": [1500.0, 500.0]},
    )
    result = CorrectionResult(
        boa=xr.Dataset(
            {
                "B04": xr.DataArray(np.full((2, 2), 0.12, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
                "B03": xr.DataArray(np.full((2, 2), 0.10, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
                "B02": xr.DataArray(np.full((2, 2), 0.08, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
            }
        ),
        boa_unc=None,
        aot=xr.DataArray(np.full((2, 2), 0.2, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
        tcwv=xr.DataArray(np.full((2, 2), 2.0, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
        cloud_mask=shifted_mask,
    )

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="netcdf",
            include_uncertainty=False,
            include_auxiliary=True,
            include_rgb=False,
        )
    )

    writer.write(result, tmp_path)

    aux_ds = captured[str(tmp_path / "auxiliary.nc")]
    expected = np.array([[False, False], [True, False]], dtype=bool)
    np.testing.assert_array_equal(aux_ds["cloud_mask"].values.astype(bool), expected)


def test_write_raster_products_emits_surface_prior_and_monthly_composites(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset_calls: list[tuple[xr.Dataset, Path, dict[str, object]]] = []
    raster_calls: list[tuple[Path, dict[str, object], np.dtype[np.generic]]] = []

    def _fake_write_dataset(dataset: xr.Dataset, output_dir: Path, **kwargs: object) -> dict[str, Path]:
        dataset_calls.append((dataset, output_dir, kwargs))
        return {name: output_dir / f"{name}.tif" for name in dataset.data_vars}

    def _fake_write_cog(data: xr.DataArray, path: Path, **kwargs: object) -> Path:
        raster_calls.append((path, kwargs, data.dtype))
        return path

    monkeypatch.setattr(output_module, "write_dataset", _fake_write_dataset)
    monkeypatch.setattr(output_module, "write_cog", _fake_write_cog)
    monkeypatch.setattr(output_module, "write_rgb_quicklook", lambda _dataset, path: path)
    monkeypatch.setattr(output_module, "write_aot_scatter_plot", lambda _scatter, path: path)

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="cog",
            boa_dtype="uint16",
            include_uncertainty=True,
            include_auxiliary=False,
            include_rgb=False,
        )
    )
    artifacts = writer.write(
        _result(include_surface_prior=True, include_monthly_composites=True),
        tmp_path,
    )

    assert [output_dir for _, output_dir, _ in dataset_calls] == [
        tmp_path / "boa",
        tmp_path / "boa_unc",
        tmp_path / "surface_prior",
        tmp_path / "surface_prior_unc",
        tmp_path / "monthly_composites" / "2023_07",
    ]
    assert {path.name for path, _, _ in raster_calls} == {"quality.tif", "sample_index.tif"}
    assert {
        "surface_prior.B04",
        "surface_prior.B03",
        "surface_prior_unc.B04",
        "surface_prior_unc.B03",
        "monthly_composites.2023_07.B04",
        "monthly_composites.2023_07.B03",
        "monthly_composites.2023_07.quality",
        "monthly_composites.2023_07.sample_index",
    } <= set(artifacts)


def test_write_zarr_products_route_and_respects_toggles(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[Path] = []

    def _fake_writer(dataset: xr.Dataset, path: Path) -> Path:
        calls.append(path)
        return path

    monkeypatch.setattr(output_module, "write_zarr", _fake_writer)
    monkeypatch.setattr(output_module, "write_aot_scatter_plot", lambda _scatter, path: path)

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


def test_write_zarr_products_include_surface_prior_and_monthly_composites(
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
            include_uncertainty=True,
            include_auxiliary=False,
            include_rgb=False,
        )
    )
    artifacts = writer.write(
        _result(include_surface_prior=True, include_monthly_composites=True),
        tmp_path,
    )

    assert calls == [
        tmp_path / "boa.zarr",
        tmp_path / "boa_unc.zarr",
        tmp_path / "surface_prior.zarr",
        tmp_path / "surface_prior_unc.zarr",
        tmp_path / "monthly_composites" / "2023_07.zarr",
    ]
    assert {
        "surface_prior",
        "surface_prior_unc",
        "monthly_composites.2023_07",
    } <= set(artifacts)


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
