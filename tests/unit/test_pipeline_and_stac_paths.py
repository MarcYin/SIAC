from __future__ import annotations

import json
from concurrent.futures import TimeoutError as FuturesTimeoutError
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

from siac.adapters.brdf.gee_stub import GEEBRDFProvider
from siac.adapters.s2_backend import build_s2_backend
from siac.algorithms.surface.brdf_monthly_composite import build_monthly_best_pixel_composite
from siac.algorithms.surface.brdf_monthly_database import (
    _dataset_to_cube,
    build_monthly_composite_database,
)
from siac.storage import stac as stac_mod
from siac.workflows import pipeline


def _composite_inputs() -> tuple[xr.DataArray, xr.DataArray]:
    reflectance = xr.DataArray(
        np.ones((2, 2, 1, 1), dtype=np.float32),
        dims=["time", "band", "y", "x"],
        coords={"time": [0, 1], "band": ["B02", "B08"], "y": [0], "x": [0]},
    )
    quality = xr.DataArray(
        np.ones((2, 1, 1), dtype=np.float32),
        dims=["time", "y", "x"],
        coords={"time": [0, 1], "y": [0], "x": [0]},
    )
    return reflectance, quality


def _make_monthly_composites() -> list:
    composites = []
    bands = ["B02", "B08", "B11", "B12"]
    for idx in range(15):
        cube = np.zeros((len(bands), 1, 1), dtype=np.float32)
        cube[bands.index("B02"), 0, 0] = 0.1 + idx * 0.01
        cube[bands.index("B08"), 0, 0] = 0.4 + idx * 0.01
        cube[bands.index("B11"), 0, 0] = 0.3 + idx * 0.01
        cube[bands.index("B12"), 0, 0] = 0.2 + idx * 0.01
        reflectance = xr.DataArray(
            cube,
            dims=["band", "y", "x"],
            coords={"band": bands, "y": [0], "x": [0]},
        )
        composites.append(
            SimpleNamespace(
                reflectance=reflectance,
                quality=xr.DataArray(np.zeros((1, 1), dtype=np.int16), dims=["y", "x"], coords={"y": [0], "x": [0]}),
                sample_index=xr.DataArray(np.zeros((1, 1), dtype=np.int16), dims=["y", "x"], coords={"y": [0], "x": [0]}),
                year=2024,
                month=idx + 1,
            )
        )
    return composites


def test_small_adapter_and_monthly_validation_paths() -> None:
    with pytest.raises(ValueError, match="Unknown S2 backend"):
        build_s2_backend(SimpleNamespace(s2_data=SimpleNamespace(backend="weird")))

    provider = GEEBRDFProvider("arg", option=True)
    with pytest.raises(NotImplementedError, match="not implemented"):
        provider()

    reflectance, quality = _composite_inputs()
    with pytest.raises(ValueError, match="reflectance must have dims"):
        build_monthly_best_pixel_composite(reflectance.isel(time=0), quality, year=2024, month=1)
    with pytest.raises(ValueError, match="quality must have dims"):
        build_monthly_best_pixel_composite(reflectance, quality.isel(time=0), year=2024, month=1)
    with pytest.raises(ValueError, match="same time axis"):
        build_monthly_best_pixel_composite(
            reflectance,
            quality.isel(time=slice(0, 1)),
            year=2024,
            month=1,
        )
    with pytest.raises(ValueError, match="same spatial shape"):
        build_monthly_best_pixel_composite(
            reflectance,
            xr.DataArray(np.ones((2, 2, 2), dtype=np.float32), dims=["time", "y", "x"]),
            year=2024,
            month=1,
        )

    composites = _make_monthly_composites()
    with pytest.raises(ValueError, match="query_bands must not be empty"):
        build_monthly_composite_database(composites, query_bands=(), visible_bands=("B02",))
    with pytest.raises(ValueError, match="visible_bands must not be empty"):
        build_monthly_composite_database(composites, query_bands=("B08",), visible_bands=())

    bad_composites = list(composites)
    bad_composites[1] = SimpleNamespace(
        reflectance=xr.DataArray(
            np.ones((4, 2, 1), dtype=np.float32),
            dims=["band", "y", "x"],
            coords={"band": ["B02", "B08", "B11", "B12"], "y": [0, 1], "x": [0]},
        ),
        quality=bad_composites[1].quality,
        sample_index=bad_composites[1].sample_index,
        year=2024,
        month=2,
    )
    with pytest.raises(ValueError, match="same spatial shape"):
        build_monthly_composite_database(bad_composites, query_bands=("B08",), visible_bands=("B02",))

    with pytest.raises(ValueError, match="band"):
        _dataset_to_cube(xr.DataArray(np.ones((1, 1), dtype=np.float32), dims=["y", "x"]), ["B08"])

    database = build_monthly_composite_database(
        composites,
        query_bands=("B08", "B11", "B12"),
        visible_bands=("B02",),
    )
    with pytest.raises(ValueError, match="k_neighbors must be >= 1"):
        database.predict_visible(
            xr.Dataset({"B08": xr.DataArray([[0.1]], dims=["y", "x"]), "B11": xr.DataArray([[0.1]], dims=["y", "x"]), "B12": xr.DataArray([[0.1]], dims=["y", "x"])}),
            k_neighbors=0,
        )
    visible, visible_unc, visible_quality = database.predict_visible(
        xr.Dataset(
            {
                "B08": xr.DataArray(np.ones((2, 1), dtype=np.float32), dims=["y", "x"]),
                "B11": xr.DataArray(np.ones((2, 1), dtype=np.float32), dims=["y", "x"]),
                "B12": xr.DataArray(np.ones((2, 1), dtype=np.float32), dims=["y", "x"]),
            }
        ),
    )
    assert visible.shape == (1, 2, 1)
    assert visible_unc.shape == visible.shape
    assert visible_quality.shape == (2, 1)


def test_stac_helper_functions_cover_fallback_and_optional_paths(tmp_path: Path) -> None:
    assert stac_mod._parse_satellite_id(None, {"satellite": "S2B"}, "fallback") == "S2B"
    assert stac_mod._parse_satellite_id(None, {"satellite": ""}, "fallback") == "fallback"
    assert stac_mod._platform_name("L8") == "landsat-8"
    assert stac_mod._platform_name("unknown") is None
    assert stac_mod._constellation_name("L9") == "landsat"
    assert stac_mod._constellation_name("unknown") is None
    assert stac_mod._safe_float("bad") is None
    assert stac_mod._safe_float(np.inf) is None
    assert stac_mod._mean_deg(xr.DataArray(np.array([], dtype=np.float32), dims=["y"])) is None
    assert stac_mod._mean_deg(xr.DataArray(np.array([np.nan], dtype=np.float32), dims=["y"])) is None
    assert stac_mod._cloud_cover_percent(xr.DataArray(np.array([np.nan], dtype=np.float32), dims=["y"])) is None

    first_band = SimpleNamespace(
        rio=SimpleNamespace(
            bounds=lambda: (_ for _ in ()).throw(RuntimeError("no bounds")),
            crs=SimpleNamespace(to_epsg=lambda: (_ for _ in ()).throw(RuntimeError("no epsg"))),
            transform=lambda: (_ for _ in ()).throw(RuntimeError("no transform")),
            resolution=lambda: (_ for _ in ()).throw(RuntimeError("no resolution")),
        ),
        sizes={"y": 2, "x": 3},
    )
    assert stac_mod._native_bounds(first_band, (1.0, 2.0, 3.0, 4.0)) == (1.0, 2.0, 3.0, 4.0)
    bbox, geometry = stac_mod._wgs84_bounds_and_geometry((1.0, 2.0, 3.0, 4.0), None)
    assert bbox == [1.0, 2.0, 3.0, 4.0]
    assert geometry["coordinates"][0][0] == [1.0, 2.0]
    assert stac_mod._proj_properties(first_band, (1.0, 2.0, 3.0, 4.0)) == {
        "proj:bbox": [1.0, 2.0, 3.0, 4.0],
        "proj:shape": [2, 3],
    }
    assert stac_mod._gsd(first_band) is None
    assert "common_name" not in stac_mod._band_metadata(SimpleNamespace(name="UNK", wavelength_um=0.5, bandwidth=20.0))
    assert stac_mod._file_size(tmp_path / "missing.bin") is None
    assert stac_mod._asset_dict("file.tif", title="T", media_type="image/tiff", roles=["data"]) == {
        "href": "file.tif",
        "type": "image/tiff",
        "title": "T",
        "roles": ["data"],
    }

    template = xr.DataArray(np.ones((1, 1), dtype=np.float32), dims=["y", "x"], coords={"y": [0], "x": [0]})
    obs = SimpleNamespace(
        bounds=(1.0, 2.0, 3.0, 4.0),
        crs="EPSG:4326",
        metadata={"observation_time": "bad"},
        cloud_mask=xr.DataArray(np.zeros((1, 1), dtype=bool), dims=["y", "x"]),
        geometry=SimpleNamespace(
            sza=template,
            saa=template,
            vza=template,
            vaa=template,
            raa=template,
        ),
        sensor_config=SimpleNamespace(
            sensor_id="msi",
            satellite_id="S2A",
            get_band=lambda name: SimpleNamespace(name=name, wavelength_um=0.5, bandwidth=20.0),
        ),
    )
    result = SimpleNamespace(
        boa=xr.Dataset({"B02": template}),
        aot=template,
        tcwv=template,
        cloud_mask=xr.DataArray(np.zeros((1, 1), dtype=bool), dims=["y", "x"]),
        diagnostics=SimpleNamespace(processing_time_s=None),
    )
    with pytest.raises(TypeError, match="datetime observation_time"):
        stac_mod.build_stac_item(obs, result, output_dir=tmp_path, boa_assets={"B02": tmp_path / "B02.tif"})

    good_obs = SimpleNamespace(**{**obs.__dict__, "metadata": {"observation_time": datetime(2024, 1, 1)}})
    path = stac_mod.write_stac_item(
        good_obs,
        result,
        output_dir=tmp_path,
        boa_assets={"B02": tmp_path / "B02.tif"},
        item_id="item-1",
        item_href=tmp_path / "custom.json",
    )
    assert path == tmp_path / "custom.json"
    written = json.loads(path.read_text(encoding="utf-8"))
    assert written["id"] == "item-1"
    assert written["links"] == [{"rel": "self", "href": "custom.json", "type": "application/geo+json"}]


def test_pipeline_helper_branches_cover_fallbacks_and_preload_failures(
    monkeypatch: pytest.MonkeyPatch,
    mock_observation_bundle,
    mock_atmospheric_state,
    mock_surface_prior,
    mock_grid_assembler,
    mock_solver_fn,
    mock_corrector_fn,
    mock_rt_model,
    caplog: pytest.LogCaptureFixture,
) -> None:
    sensor_config = SimpleNamespace(
        bands=["fallback-1", "fallback-2", "fallback-3"],
        select_bands_in_range=lambda *_args: [],
    )
    assert pipeline._select_solver_bands_for_preload(sensor_config) == ["fallback-1", "fallback-2"]

    with pytest.raises(RuntimeError, match="Unreachable retry state"):
        pipeline._call_with_retries(lambda: None, (), retries=-1, stage_name="M0")

    monkeypatch.setattr(
        pipeline,
        "_prepare_observation",
        lambda *_args, **_kwargs: (mock_observation_bundle, mock_observation_bundle.bounds, mock_observation_bundle.crs, mock_observation_bundle.metadata["observation_time"], 1000.0),
    )
    monkeypatch.setattr(pipeline, "_run_tail", lambda *_args, **_kwargs: "ok")

    class _Future:
        def __init__(self, exc: Exception | None = None) -> None:
            self.exc = exc
            self.canceled = False

        def result(self, timeout=None):  # noqa: ANN001
            if self.exc is not None:
                raise self.exc
            return None

        def cancel(self) -> None:
            self.canceled = True

    def _fake_preload_timeout(*_args, **_kwargs):
        return _Future(FuturesTimeoutError())

    with caplog.at_level("WARNING"):
        monkeypatch.setattr(pipeline, "_maybe_submit_lut_preload", _fake_preload_timeout)
        out = pipeline._run_pipeline_thread(
            Path("/fake"),
            None,
            None,
            preprocessor=lambda *_args: mock_observation_bundle,
            atmo_provider=lambda *_args: mock_atmospheric_state,
            surface_prior_provider=lambda *_args: mock_surface_prior,
            grid_assembler=mock_grid_assembler,
            solver=mock_solver_fn,
            corrector=mock_corrector_fn,
            rt_model=mock_rt_model,
            settings={"max_workers": 2, "retries": 0, "stage_timeout_s": 1.0},
        )
    assert out == "ok"
    assert "timed out" in caplog.text

    class _FailingProvider:
        requires_atmo_prior = True

        def __call__(self, *_args):  # noqa: ANN002
            return mock_surface_prior

    with caplog.at_level("WARNING"):
        monkeypatch.setattr(pipeline, "_maybe_submit_lut_preload", lambda *_args, **_kwargs: _Future(RuntimeError("boom")))
        out = pipeline._run_pipeline_thread(
            Path("/fake"),
            None,
            None,
            preprocessor=lambda *_args: mock_observation_bundle,
            atmo_provider=lambda *_args: mock_atmospheric_state,
            surface_prior_provider=_FailingProvider(),
            grid_assembler=mock_grid_assembler,
            solver=mock_solver_fn,
            corrector=mock_corrector_fn,
            rt_model=mock_rt_model,
            settings={"max_workers": 2, "retries": 0, "stage_timeout_s": 1.0},
        )
    assert out == "ok"
    assert "LUT preload failed" in caplog.text
