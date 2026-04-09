from __future__ import annotations

import json
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

import siac.adapters.output as output_module
from siac.adapters.output import ConfiguredOutputWriter, _derive_scene_prefix
from siac.config.schema import OutputDefaultsConfig
from siac.runtime import (
    AOTScatterBandDiagnostics,
    CorrectionDiagnostics,
    CorrectionResult,
    MonthlyCompositeOutput,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _result(
    *,
    include_uncertainty: bool = True,
    include_surface_prior: bool = False,
    include_monthly_composites: bool = False,
    include_diagnostics: bool = False,
    include_solver_qa: bool = False,
    metadata: dict | None = None,
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
    if metadata is None:
        metadata = {
            "satellite": "S2A",
            "observation_time": datetime(2024, 3, 15, 10, 30, 45, tzinfo=timezone.utc),
            "tile_id": "32UQD",
        }
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
        metadata=metadata,
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
            )
            if include_diagnostics
            else (),
        ),
    )


# ---------------------------------------------------------------------------
# Prefix derivation
# ---------------------------------------------------------------------------


class TestScenePrefix:
    def test_sentinel2_prefix(self):
        prefix = _derive_scene_prefix({
            "satellite": "S2A",
            "observation_time": datetime(2024, 3, 15, 10, 30, 45),
            "tile_id": "32UQD",
        })
        assert prefix == "S2A_L2A_20240315T103045"

    def test_sentinel2b_prefix(self):
        prefix = _derive_scene_prefix({
            "satellite": "S2B",
            "observation_time": datetime(2024, 1, 2, 8, 15, 0),
            "tile_id": "T10SFG",
        })
        assert prefix == "S2B_L2A_20240102T081500"

    def test_landsat_prefix(self):
        prefix = _derive_scene_prefix({
            "satellite": "L8",
            "observation_time": datetime(2024, 6, 20, 14, 0, 0),
        })
        assert prefix == "L8_L2A_20240620T140000"

    def test_missing_metadata_defaults(self):
        prefix = _derive_scene_prefix({})
        assert prefix == "SAT_L2A_00000000T000000"

    def test_tile_id_is_not_included_in_output_prefix(self):
        prefix = _derive_scene_prefix({
            "satellite": "S2A",
            "observation_time": datetime(2024, 1, 1, 0, 0, 0),
            "tile_id": "T32UQD",
        })
        assert prefix == "S2A_L2A_20240101T000000"

    def test_esa_long_tile_id_is_not_included_in_output_prefix(self):
        prefix = _derive_scene_prefix({
            "satellite": "S2C",
            "observation_time": datetime(2025, 8, 26, 17, 22, 22),
            "tile_id": "S2C_OPER_MSI_L1C_TL_2CPS_20250826T204551_A005087_T15TVG_N05.11",
        })
        assert prefix == "S2C_L2A_20250826T172222"


# ---------------------------------------------------------------------------
# Raster output (COG/GeoTIFF)
# ---------------------------------------------------------------------------


def _mock_writers(monkeypatch):
    """Patch all writer functions to record calls without writing files."""
    write_fn_calls = []
    dataset_calls = []
    quicklook_calls = []
    scatter_calls = []

    def _fake_write_fn(data, path, **kwargs):
        write_fn_calls.append((path, kwargs, data.dtype))
        return path

    def _fake_write_dataset(dataset, output_dir, **kwargs):
        dataset_calls.append((dataset, output_dir, kwargs))
        return {name: output_dir / f"{name}.tif" for name in dataset.data_vars}

    def _fake_write_rgb_quicklook(dataset, path):
        quicklook_calls.append(path)
        return path

    def _fake_write_aot_scatter_plot(_scatter, path):
        scatter_calls.append(path)
        return path

    monkeypatch.setattr(output_module, "write_cog", _fake_write_fn)
    monkeypatch.setattr(output_module, "write_raster", _fake_write_fn)
    monkeypatch.setattr(output_module, "write_dataset", _fake_write_dataset)
    monkeypatch.setattr(output_module, "write_rgb_quicklook", _fake_write_rgb_quicklook)
    monkeypatch.setattr(output_module, "write_aot_scatter_plot", _fake_write_aot_scatter_plot)

    return write_fn_calls, dataset_calls, quicklook_calls, scatter_calls


def test_raster_output_uses_satellite_prefix_for_boa_bands(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    write_fn_calls, _, _, _ = _mock_writers(monkeypatch)
    monkeypatch.setattr(output_module, "build_stac_item_from_result", lambda *_a, **_kw: {"type": "Feature"})

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

    prefix = "S2A_L2A_20240315T103045"

    # BOA bands have satellite prefix
    assert "boa.B04" in artifacts
    assert artifacts["boa.B04"].name == f"{prefix}_BOA_B04.tif"
    assert artifacts["boa.B03"].name == f"{prefix}_BOA_B03.tif"
    assert artifacts["boa.B02"].name == f"{prefix}_BOA_B02.tif"

    # BOA uncertainty bands
    assert artifacts["boa_unc.B04"].name == f"{prefix}_BOA_UNC_B04.tif"

    # Auxiliary with prefix
    assert artifacts["auxiliary.aot"].name == f"{prefix}_AOT.tif"
    assert artifacts["auxiliary.tcwv"].name == f"{prefix}_TCWV.tif"
    assert artifacts["auxiliary.cloud_mask"].name == f"{prefix}_CLOUD.tif"

    # Cloud mask encoding
    cloud_call = next(c for c in write_fn_calls if c[0].name == f"{prefix}_CLOUD.tif")
    assert cloud_call[1]["compression"] == "lzw"
    assert cloud_call[1]["dtype"] == "uint8"
    assert cloud_call[1]["nodata"] == 255

    # RGB quicklook
    assert artifacts["quicklook.rgb"].name == f"{prefix}_RGB.tif"

    # Scatter diagnostics
    assert "preview.scatter.B02" in artifacts

    # STAC item
    assert "stac_item" in artifacts
    assert artifacts["stac_item"].name == f"{prefix}.json"


def test_correction_boa_stream_writes_boa_once_and_finishes_stac(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    writes: dict[Path, xr.DataArray] = {}

    def _fake_write_fn(data, path, **_kwargs):  # noqa: ANN001
        path = Path(path)
        writes[path] = data
        return path

    def _fake_read_raster(path, **_kwargs):  # noqa: ANN001
        return writes[Path(path)]

    monkeypatch.setattr(output_module, "write_cog", _fake_write_fn)
    monkeypatch.setattr("siac.storage.readers.read_raster", _fake_read_raster)
    monkeypatch.setattr(output_module, "build_stac_item_from_result", lambda *_a, **_kw: {"type": "Feature"})
    monkeypatch.setattr(output_module, "write_rgb_quicklook", lambda *_a, **_kw: None)
    monkeypatch.setattr(output_module, "write_false_colour_preview", lambda *_a, **_kw: None)
    monkeypatch.setattr(output_module, "write_field_preview", lambda *_a, **_kw: None)
    monkeypatch.setattr(output_module, "write_cloud_mask_preview", lambda *_a, **_kw: None)

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="cog",
            boa_dtype="float32",
            include_uncertainty=False,
            include_auxiliary=False,
            include_rgb=False,
        )
    )
    result = _result(include_uncertainty=False)
    stream = writer.open_correction_boa_stream(tmp_path, metadata=result.metadata)

    assert stream is not None
    streamed_boa = xr.Dataset(
        {
            name: stream.write_boa_band(name, field)
            for name, field in result.boa.data_vars.items()
        }
    )
    streamed_result = replace(result, boa=streamed_boa)

    artifacts = stream.finish(streamed_result)

    prefix = "S2A_L2A_20240315T103045"
    assert set(artifacts) >= {"boa.B04", "boa.B03", "boa.B02", "stac_item"}
    assert {path.name for path in writes} == {
        f"{prefix}_BOA_B04.tif",
        f"{prefix}_BOA_B03.tif",
        f"{prefix}_BOA_B02.tif",
    }
    assert stream.has_written is True


def test_correction_boa_stream_can_skip_reopen_for_speed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    writes: dict[Path, xr.DataArray] = {}

    def _fake_write_fn(data, path, **_kwargs):  # noqa: ANN001
        path = Path(path)
        writes[path] = data
        return path

    def _unexpected_read_raster(*_args, **_kwargs):  # noqa: ANN001
        raise AssertionError("read_raster should not be called when reopen_streamed_boa=false")

    monkeypatch.setattr(output_module, "write_cog", _fake_write_fn)
    monkeypatch.setattr("siac.storage.readers.read_raster", _unexpected_read_raster)
    monkeypatch.setattr(output_module, "build_stac_item_from_result", lambda *_a, **_kw: {"type": "Feature"})
    monkeypatch.setattr(output_module, "write_rgb_quicklook", lambda *_a, **_kw: None)
    monkeypatch.setattr(output_module, "write_false_colour_preview", lambda *_a, **_kw: None)
    monkeypatch.setattr(output_module, "write_field_preview", lambda *_a, **_kw: None)
    monkeypatch.setattr(output_module, "write_cloud_mask_preview", lambda *_a, **_kw: None)

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="cog",
            boa_dtype="float32",
            include_uncertainty=False,
            include_auxiliary=False,
            include_rgb=False,
            reopen_streamed_boa=False,
        )
    )
    result = _result(include_uncertainty=False)
    stream = writer.open_correction_boa_stream(tmp_path, metadata=result.metadata)

    assert stream is not None
    streamed_boa = xr.Dataset(
        {
            name: stream.write_boa_band(name, field)
            for name, field in result.boa.data_vars.items()
        }
    )
    streamed_result = replace(result, boa=streamed_boa)

    artifacts = stream.finish(streamed_result)

    prefix = "S2A_L2A_20240315T103045"
    assert set(artifacts) >= {"boa.B04", "boa.B03", "boa.B02", "stac_item"}
    assert {path.name for path in writes} == {
        f"{prefix}_BOA_B04.tif",
        f"{prefix}_BOA_B03.tif",
        f"{prefix}_BOA_B02.tif",
    }
    assert stream.has_written is True


def test_raster_output_emits_solver_qa_masks(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    write_fn_calls, _, _, _ = _mock_writers(monkeypatch)
    monkeypatch.setattr(output_module, "build_stac_item_from_result", lambda *_a, **_kw: {"type": "Feature"})

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

    prefix = "S2A_L2A_20240315T103045"
    qa_calls = [c for c in write_fn_calls if "QA_" in c[0].name]
    assert {c[0].name for c in qa_calls} == {
        f"{prefix}_QA_invalid_retrieval.tif",
        f"{prefix}_QA_low_quality.tif",
    }
    assert all(c[1]["dtype"] == "uint8" for c in qa_calls)
    assert {
        "auxiliary.qa.invalid_retrieval",
        "auxiliary.qa.low_quality",
    } <= set(artifacts)


def test_raster_output_emits_fitting_cost_as_float32(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """fitting_cost is a float32 QA field, not a boolean mask."""
    write_fn_calls, _, _, _ = _mock_writers(monkeypatch)
    monkeypatch.setattr(output_module, "build_stac_item_from_result", lambda *_a, **_kw: {"type": "Feature"})

    coords = {"y": [0, 1], "x": [0, 1]}
    result = _result(include_uncertainty=False, include_solver_qa=True)
    qa_with_cost = result.solver_qa.assign(
        fitting_cost=xr.DataArray(
            np.array([[0.01, 0.02], [0.03, 0.04]], dtype=np.float32),
            dims=["y", "x"],
            coords=coords,
        ),
    )
    result = replace(result, solver_qa=qa_with_cost)

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="cog",
            boa_dtype="uint16",
            include_uncertainty=False,
            include_auxiliary=True,
            include_rgb=False,
        )
    )
    artifacts = writer.write(result, tmp_path)

    cost_calls = [c for c in write_fn_calls if "QA_fitting_cost" in c[0].name]
    assert len(cost_calls) == 1
    assert cost_calls[0][1]["dtype"] == "float32"

    # Boolean QA fields should still use uint8
    bool_calls = [c for c in write_fn_calls if "QA_" in c[0].name and "fitting_cost" not in c[0].name]
    assert all(c[1]["dtype"] == "uint8" for c in bool_calls)

    assert "auxiliary.qa.fitting_cost" in artifacts


def test_raster_output_emits_surface_prior_and_monthly_composites(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    write_fn_calls, dataset_calls, _, _ = _mock_writers(monkeypatch)
    monkeypatch.setattr(output_module, "build_stac_item_from_result", lambda *_a, **_kw: {"type": "Feature"})

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

    prefix = "S2A_L2A_20240315T103045"
    # Surface prior uses prefix
    assert artifacts["surface_prior.B04"].name == f"{prefix}_SURF_B04.tif"
    assert artifacts["surface_prior.B03"].name == f"{prefix}_SURF_B03.tif"
    assert artifacts["surface_prior_unc.B04"].name == f"{prefix}_SURF_UNC_B04.tif"

    # Monthly composites still use subdirectory (too many files)
    assert {path.name for path, _, _ in write_fn_calls if "monthly" in str(path)} == {"quality.tif", "sample_index.tif"}
    assert {
        "monthly_composites.2023_07.B04",
        "monthly_composites.2023_07.B03",
        "monthly_composites.2023_07.quality",
        "monthly_composites.2023_07.sample_index",
    } <= set(artifacts)


# ---------------------------------------------------------------------------
# NetCDF output
# ---------------------------------------------------------------------------


def test_netcdf_output_uses_satellite_prefix(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[Path] = []

    def _fake_write_netcdf(dataset, path, **kwargs):
        del dataset, kwargs
        calls.append(path)
        return path

    monkeypatch.setattr(output_module, "write_netcdf", _fake_write_netcdf)
    monkeypatch.setattr(output_module, "write_aot_scatter_plot", lambda _s, path: path)
    monkeypatch.setattr(output_module, "build_stac_item_from_result", lambda *_a, **_kw: {"type": "Feature"})

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="netcdf",
            include_uncertainty=False,
            include_auxiliary=False,
            include_rgb=False,
        )
    )
    artifacts = writer.write(_result(include_uncertainty=False), tmp_path)

    prefix = "S2A_L2A_20240315T103045"
    assert calls == [tmp_path / f"{prefix}_BOA.nc"]
    assert artifacts["boa"] == tmp_path / f"{prefix}_BOA.nc"


def test_netcdf_auxiliary_aligns_cloud_mask_to_atmo_grid(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, xr.Dataset] = {}

    def _fake_write_netcdf(dataset, path, **kwargs):
        del kwargs
        captured[str(path)] = dataset
        return path

    monkeypatch.setattr(output_module, "write_netcdf", _fake_write_netcdf)
    monkeypatch.setattr(output_module, "build_stac_item_from_result", lambda *_a, **_kw: {"type": "Feature"})

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
        metadata={"satellite": "S2A", "observation_time": datetime(2024, 1, 1, tzinfo=timezone.utc)},
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

    aux_key = next(k for k in captured if "AUX" in k)
    aux_ds = captured[aux_key]
    assert "auxiliary" in artifacts
    assert aux_ds.sizes == {"y": 2, "x": 2}
    assert aux_ds["cloud_mask"].coords["x"].identical(result.aot.coords["x"])
    assert aux_ds["cloud_mask"].coords["y"].identical(result.aot.coords["y"])
    assert aux_ds["cloud_mask"].dtype == np.uint8


def test_netcdf_auxiliary_includes_solver_qa_masks(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, xr.Dataset] = {}

    def _fake_write_netcdf(dataset, path, **kwargs):
        del kwargs
        captured[str(path)] = dataset
        return path

    monkeypatch.setattr(output_module, "write_netcdf", _fake_write_netcdf)
    monkeypatch.setattr(output_module, "build_stac_item_from_result", lambda *_a, **_kw: {"type": "Feature"})

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="netcdf",
            include_uncertainty=False,
            include_auxiliary=True,
            include_rgb=False,
        )
    )
    writer.write(_result(include_uncertainty=False, include_solver_qa=True), tmp_path)

    aux_key = next(k for k in captured if "AUX" in k)
    aux_ds = captured[aux_key]
    assert {"invalid_retrieval", "low_quality"} <= set(aux_ds.data_vars)
    assert aux_ds["invalid_retrieval"].dtype == np.uint8


# ---------------------------------------------------------------------------
# Zarr output
# ---------------------------------------------------------------------------


def test_zarr_output_uses_satellite_prefix(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[Path] = []

    def _fake_writer(dataset, path):
        calls.append(path)
        return path

    monkeypatch.setattr(output_module, "write_zarr", _fake_writer)
    monkeypatch.setattr(output_module, "write_aot_scatter_plot", lambda _s, path: path)
    monkeypatch.setattr(output_module, "build_stac_item_from_result", lambda *_a, **_kw: {"type": "Feature"})

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="zarr",
            include_uncertainty=False,
            include_auxiliary=False,
            include_rgb=False,
        )
    )
    artifacts = writer.write(_result(include_uncertainty=False), tmp_path)

    prefix = "S2A_L2A_20240315T103045"
    assert calls == [tmp_path / f"{prefix}_BOA.zarr"]
    assert artifacts["boa"] == tmp_path / f"{prefix}_BOA.zarr"


def test_zarr_auxiliary_aligns_cloud_mask_to_atmo_grid(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, xr.Dataset] = {}

    def _fake_write_zarr(dataset, path):
        captured[str(path)] = dataset
        return path

    monkeypatch.setattr(output_module, "write_zarr", _fake_write_zarr)
    monkeypatch.setattr(output_module, "build_stac_item_from_result", lambda *_a, **_kw: {"type": "Feature"})

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
        metadata={"satellite": "S2A", "observation_time": datetime(2024, 1, 1, tzinfo=timezone.utc)},
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

    aux_key = next(k for k in captured if "AUX" in k)
    aux_ds = captured[aux_key]
    assert "auxiliary" in artifacts
    assert aux_ds.sizes == {"y": 2, "x": 2}
    assert aux_ds["cloud_mask"].dtype == np.uint8


def test_zarr_surface_prior_and_monthly_composites(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[Path] = []

    def _fake_writer(dataset, path):
        calls.append(path)
        return path

    monkeypatch.setattr(output_module, "write_zarr", _fake_writer)
    monkeypatch.setattr(output_module, "build_stac_item_from_result", lambda *_a, **_kw: {"type": "Feature"})

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

    prefix = "S2A_L2A_20240315T103045"
    assert calls == [
        tmp_path / f"{prefix}_BOA.zarr",
        tmp_path / f"{prefix}_BOA_UNC.zarr",
        tmp_path / f"{prefix}_SURF.zarr",
        tmp_path / f"{prefix}_SURF_UNC.zarr",
        tmp_path / "monthly_composites" / "2023_07.zarr",
    ]
    assert {
        "surface_prior",
        "surface_prior_unc",
        "monthly_composites.2023_07",
    } <= set(artifacts)


def test_netcdf_auxiliary_reorders_same_shape_cloud_mask_coords(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, xr.Dataset] = {}

    def _fake_write_netcdf(dataset, path, **kwargs):
        del kwargs
        captured[str(path)] = dataset
        return path

    monkeypatch.setattr(output_module, "write_netcdf", _fake_write_netcdf)
    monkeypatch.setattr(output_module, "build_stac_item_from_result", lambda *_a, **_kw: {"type": "Feature"})

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
        metadata={"satellite": "S2A", "observation_time": datetime(2024, 1, 1, tzinfo=timezone.utc)},
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

    aux_key = next(k for k in captured if "AUX" in k)
    aux_ds = captured[aux_key]
    expected = np.array([[False, False], [True, False]], dtype=bool)
    np.testing.assert_array_equal(aux_ds["cloud_mask"].values.astype(bool), expected)


# ---------------------------------------------------------------------------
# STAC item generation
# ---------------------------------------------------------------------------


def test_stac_item_generated_automatically(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The writer should always produce a STAC JSON file."""
    write_fn_calls, _, _, _ = _mock_writers(monkeypatch)

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="cog",
            include_uncertainty=False,
            include_auxiliary=False,
            include_rgb=False,
        )
    )
    artifacts = writer.write(_result(include_uncertainty=False), tmp_path)

    stac_path = artifacts["stac_item"]
    assert stac_path.exists()
    item = json.loads(stac_path.read_text())

    assert item["type"] == "Feature"
    assert item["stac_version"] == "1.1.0"
    props = item["properties"]
    assert "datetime" in props
    assert props["siac:satellite"] == "S2A"
    assert "siac:aot_mean" in props
    assert "siac:tcwv_mean" in props
    assert "eo:bands" in props
    assert len(props["eo:bands"]) == 3
    assert "processing:software" in props

    # Assets should include BOA bands
    assert "B04" in item["assets"]
    assert "B03" in item["assets"]
    assert "B02" in item["assets"]
    assert item["assets"]["B04"]["roles"] == ["data", "reflectance"]


def test_stac_item_includes_view_geometry_when_available(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """View extension props present when geometry is in metadata."""
    _mock_writers(monkeypatch)

    from siac.runtime import GeometryAngles

    shape = (2, 2)
    coords = {"y": [0, 1], "x": [0, 1]}
    geometry = GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.5, dtype=np.float32), dims=["y", "x"], coords=coords),
        saa=xr.DataArray(np.full(shape, 2.5, dtype=np.float32), dims=["y", "x"], coords=coords),
        vza=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"], coords=coords),
        vaa=xr.DataArray(np.full(shape, 1.5, dtype=np.float32), dims=["y", "x"], coords=coords),
    )
    meta = {
        "satellite": "S2A",
        "observation_time": datetime(2024, 3, 15, 10, 30, 45, tzinfo=timezone.utc),
        "tile_id": "32UQD",
        "geometry": geometry,
    }

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(format="cog", include_uncertainty=False, include_auxiliary=False, include_rgb=False)
    )
    artifacts = writer.write(_result(include_uncertainty=False, metadata=meta), tmp_path)

    item = json.loads(artifacts["stac_item"].read_text())
    props = item["properties"]
    assert "view:sun_azimuth" in props
    assert "view:sun_elevation" in props
    assert "view:off_nadir" in props
    assert "view:azimuth" in props


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


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


def test_output_with_empty_metadata_uses_fallback_prefix(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _mock_writers(monkeypatch)
    monkeypatch.setattr(output_module, "build_stac_item_from_result", lambda *_a, **_kw: {"type": "Feature"})

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(format="cog", include_uncertainty=False, include_auxiliary=False, include_rgb=False)
    )
    artifacts = writer.write(_result(include_uncertainty=False, metadata={}), tmp_path)

    # Should use fallback prefix
    boa_path = artifacts["boa.B04"]
    assert boa_path.name == "SAT_L2A_00000000T000000_BOA_B04.tif"


def test_skip_correction_writes_auxiliary_without_boa(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Empty BOA (skip_correction mode): no BOA bands written, but AOT/TCWV/cloud are."""
    _mock_writers(monkeypatch)
    monkeypatch.setattr(output_module, "build_stac_item_from_result", lambda *_a, **_kw: {"type": "Feature"})

    coords = {"y": [0, 1], "x": [0, 1]}
    result = CorrectionResult(
        boa=xr.Dataset(),
        boa_unc=None,
        aot=xr.DataArray(np.full((2, 2), 0.15, dtype=np.float32), dims=["y", "x"], coords=coords),
        tcwv=xr.DataArray(np.full((2, 2), 1.5, dtype=np.float32), dims=["y", "x"], coords=coords),
        cloud_mask=xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"], coords=coords),
        metadata={"skip_correction": True},
    )

    writer = ConfiguredOutputWriter(
        OutputDefaultsConfig(
            format="cog",
            include_uncertainty=False,
            include_auxiliary=True,
            include_rgb=False,
        )
    )
    artifacts = writer.write(result, tmp_path)

    assert not any(k.startswith("boa.") for k in artifacts), "No BOA bands should be written"
    assert "auxiliary.aot" in artifacts
    assert "auxiliary.tcwv" in artifacts
    assert "auxiliary.cloud_mask" in artifacts
