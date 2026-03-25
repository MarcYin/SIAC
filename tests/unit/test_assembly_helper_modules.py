from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.grid import assembler as grid_assembler_mod
from siac.algorithms.surface.brdf_monthly_composite import (
    MonthlyCompositeCollection,
    MonthlyKernelWeightComposite,
)
from siac.algorithms.surface.monthly_composite_store import (
    MonthlyCompositeStoreGridSpec,
    MonthlyCompositeStoreManifest,
    write_monthly_composite_collection,
)
from siac.app import _assembly_providers as providers_mod
from siac.app import _assembly_runtime as runtime_mod
from siac.app import _assembly_surface as surface_mod
from siac.domain import SensorBand, SensorConfig
from siac.runtime import BRDFKernelWeights, SolvedAtmosphere, SolverInputBundle


def _cloud_mask_config() -> object:
    return SimpleNamespace(model_dump=lambda **_kwargs: {"method": "mock"})


def test_build_preprocessor_runtime_supports_auto_sensor_and_typeerror_fallback(
    mock_sensor_config,
) -> None:
    toa = xr.Dataset(
        {"B02": xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"])}
    )
    raw = {
        "toa": toa,
        "geometry": "geometry",
        "cloud_mask": xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"]),
        "metadata": {"observation_time": "2024-01-01T00:00:00"},
    }
    fake_preprocessor = SimpleNamespace(
        sensor_config=mock_sensor_config,
        config={},
        preprocess=lambda _path: raw,
    )
    calls: dict[str, object] = {}

    def _detect_sensor(path: Path) -> str:
        calls["detect"] = path
        return "sentinel2"

    def _get_preprocessor(sensor_name: str, config: dict[str, object] | None = None) -> object:
        calls.setdefault("configs", []).append(config)
        if config is not None:
            raise TypeError("legacy signature")
        calls["sensor_name"] = sensor_name
        return fake_preprocessor

    config = SimpleNamespace(
        sensor="auto",
        cloud_mask=_cloud_mask_config(),
        paths=SimpleNamespace(rsrf_root="/tmp/rsrf"),
    )
    default_aoi = SimpleNamespace(crs="EPSG:4326", get_bounds=lambda: (1.0, 2.0, 3.0, 4.0))

    runtime = runtime_mod.build_preprocessor_runtime(
        config,
        input_path=Path("scene.zip"),
        default_aoi_resolver=lambda _toa: default_aoi,
        detect_sensor_fn=_detect_sensor,
        get_preprocessor_fn=_get_preprocessor,
    )
    result = runtime.preprocessor(Path("scene.zip"))

    assert runtime.sensor_config is mock_sensor_config
    assert calls["detect"] == Path("scene.zip")
    assert calls["sensor_name"] == "sentinel2"
    assert fake_preprocessor.config["cloud_mask"] == {"method": "mock"}
    assert fake_preprocessor.config["rsrf_root"] == "/tmp/rsrf"
    assert result.crs == "EPSG:4326"
    assert result.bounds == (1.0, 2.0, 3.0, 4.0)


def test_build_preprocessor_runtime_requires_input_path_for_auto_sensor() -> None:
    config = SimpleNamespace(sensor="auto", cloud_mask=_cloud_mask_config(), paths=None)

    with pytest.raises(ValueError, match="without an input path"):
        runtime_mod.build_preprocessor_runtime(config)


def test_build_preprocessor_runtime_uses_module_defaults_and_explicit_aoi(
    monkeypatch: pytest.MonkeyPatch,
    mock_sensor_config,
) -> None:
    toa = xr.Dataset(
        {"B02": xr.DataArray(np.ones((1, 1), dtype=np.float32), dims=["y", "x"])}
    )
    fake_preprocessor = SimpleNamespace(
        sensor_config=mock_sensor_config,
        preprocess=lambda _path: {
            "toa": toa,
            "geometry": "geometry",
            "cloud_mask": xr.DataArray(np.zeros((1, 1), dtype=bool), dims=["y", "x"]),
            "metadata": {"observation_time": "2024-01-02T00:00:00"},
        },
    )
    calls: dict[str, object] = {}

    def _detect_sensor(path: Path) -> str:
        calls["detect"] = path
        return "s2"

    monkeypatch.setattr(runtime_mod, "detect_sensor", _detect_sensor)

    def _get_preprocessor(sensor_name: str, config: dict[str, object] | None = None) -> object:
        calls["sensor_name"] = sensor_name
        calls["config"] = config
        return fake_preprocessor

    monkeypatch.setattr(runtime_mod, "get_preprocessor", _get_preprocessor)

    config = SimpleNamespace(sensor="auto", cloud_mask=_cloud_mask_config(), paths=None)
    explicit_aoi = SimpleNamespace(crs="EPSG:4326", get_bounds=lambda: (4.0, 5.0, 6.0, 7.0))

    runtime = runtime_mod.build_preprocessor_runtime(config, input_path=Path("scene.zip"))
    result = runtime.preprocessor(Path("scene.zip"), aoi=explicit_aoi)

    assert calls["detect"] == Path("scene.zip")
    assert calls["sensor_name"] == "s2"
    assert calls["config"] == {"cloud_mask": {"method": "mock"}}
    assert result.crs == "EPSG:4326"
    assert result.bounds == (4.0, 5.0, 6.0, 7.0)


def test_build_preprocessor_runtime_default_aoi_and_unknown_sensor(
    monkeypatch: pytest.MonkeyPatch,
    mock_sensor_config,
) -> None:
    toa = xr.Dataset(
        {"B02": xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"])}
    )
    raw = {
        "toa": toa,
        "geometry": "geometry",
        "cloud_mask": xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"]),
        "metadata": {"observation_time": "2024-01-01T00:00:00"},
    }
    def _default_aoi(_raster: object) -> object:
        return SimpleNamespace(crs="EPSG:3857", get_bounds=lambda: (0, 1, 2, 3))

    def _build_preprocessor(_sensor_name: str, config: dict[str, object] | None = None) -> object:  # noqa: ARG001
        return SimpleNamespace(
            sensor_config=mock_sensor_config,
            preprocess=lambda _path: raw,
        )

    def _raise_unknown(sensor_name: str, config: dict[str, object] | None = None) -> object:  # noqa: ARG001
        raise KeyError(sensor_name)

    config = SimpleNamespace(sensor="s2", cloud_mask=_cloud_mask_config(), paths=None)
    monkeypatch.setattr(runtime_mod.AOI, "from_raster", _default_aoi)

    runtime = runtime_mod.build_preprocessor_runtime(
        config,
        get_preprocessor_fn=_build_preprocessor,
    )
    result = runtime.preprocessor(Path("scene.zip"))
    assert result.crs == "EPSG:3857"
    assert result.bounds == (0, 1, 2, 3)

    with pytest.raises(ValueError, match="Unknown sensor"):
        runtime_mod.build_preprocessor_runtime(
            config,
            get_preprocessor_fn=_raise_unknown,
        )


def test_build_preprocessor_runtime_primes_path_sensitive_sensor_config() -> None:
    toa = xr.Dataset(
        {"B02": xr.DataArray(np.ones((1, 1), dtype=np.float32), dims=["y", "x"])}
    )
    generic_sensor_config = object()
    resolved_sensor_config = object()
    raw = {
        "toa": toa,
        "geometry": "geometry",
        "cloud_mask": xr.DataArray(np.zeros((1, 1), dtype=bool), dims=["y", "x"]),
        "metadata": {"observation_time": "2024-01-03T00:00:00"},
    }

    class _PathSensitivePreprocessor:
        def __init__(self) -> None:
            self._sensor_config = generic_sensor_config
            self.metadata_paths: list[Path] = []

        @property
        def sensor_config(self) -> object:
            return self._sensor_config

        def get_metadata(self, path: Path) -> dict[str, str]:
            self.metadata_paths.append(path)
            self._sensor_config = resolved_sensor_config
            return {"observation_time": "2024-01-03T00:00:00"}

        def preprocess(self, _path: Path) -> dict[str, object]:
            return raw

    fake_preprocessor = _PathSensitivePreprocessor()
    config = SimpleNamespace(sensor="s2", cloud_mask=_cloud_mask_config(), paths=None)

    runtime = runtime_mod.build_preprocessor_runtime(
        config,
        input_path=Path("scene.safe"),
        get_preprocessor_fn=lambda *_args, **_kwargs: fake_preprocessor,
    )
    explicit_aoi = SimpleNamespace(crs="EPSG:4326", get_bounds=lambda: (7.0, 8.0, 9.0, 10.0))
    result = runtime.preprocessor(Path("scene.safe"), aoi=explicit_aoi)

    assert fake_preprocessor.metadata_paths == [Path("scene.safe")]
    assert runtime.sensor_config is resolved_sensor_config
    assert result.sensor_config is resolved_sensor_config


def test_resolve_solver_and_rt_model_forward_expected_inputs(
    monkeypatch: pytest.MonkeyPatch,
    mock_observation_bundle,
    mock_atmospheric_state,
    mock_surface_prior,
    mock_rt_model,
) -> None:
    captured: dict[str, object] = {}

    class _FakeSolver:
        def __init__(self, config: object) -> None:
            captured["solver_config"] = config

        def solve(self, *args):  # noqa: ANN002
            captured["solve_args"] = args
            shape = mock_atmospheric_state.aot.shape
            return SimpleNamespace(
                aot=xr.DataArray(np.full(shape, 0.2, dtype=np.float32), dims=["y", "x"]),
                tcwv=xr.DataArray(np.full(shape, 3.0, dtype=np.float32), dims=["y", "x"]),
                aot_unc=xr.DataArray(np.full(shape, 0.04, dtype=np.float32), dims=["y", "x"]),
                tcwv_unc=xr.DataArray(np.full(shape, 0.2, dtype=np.float32), dims=["y", "x"]),
                final_cost=1.5,
                n_iterations=4,
                success=True,
            )

    monkeypatch.setattr(runtime_mod, "MultiGridConfig", lambda **kwargs: kwargs)
    monkeypatch.setattr(runtime_mod, "MultiGridSolver", _FakeSolver)
    monkeypatch.setattr(
        runtime_mod,
        "build_rt_model",
        lambda config, auth=None, sensor_config=None: ("rt", config, auth, sensor_config),
    )

    config = SimpleNamespace(
        solver=SimpleNamespace(
            aot_gamma=1.0,
            tcwv_gamma=2.0,
            aot_bounds=(0.0, 1.0),
            tcwv_bounds=(0.0, 5.0),
        )
    )
    toa = mock_observation_bundle.toa["B02"]
    bundle = SolverInputBundle(
        toa=toa,
        geometry=mock_observation_bundle.geometry,
        cloud_mask=mock_observation_bundle.cloud_mask,
        sensor_config=mock_observation_bundle.sensor_config,
        bands=list(mock_observation_bundle.sensor_config.bands),
        atmo_prior=mock_atmospheric_state,
        surface_prior=mock_surface_prior,
        rt_model=mock_rt_model,
        aux_resolution_m=500.0,
        aerosol_resolution_m=1000.0,
    )

    solved = runtime_mod.resolve_solver(config)(bundle, config)
    rt_model = runtime_mod.resolve_rt_model_for_pipeline(
        config,
        auth="auth",
        sensor_config=mock_observation_bundle.sensor_config,
    )

    assert isinstance(solved, SolvedAtmosphere)
    assert solved.cost_final == 1.5
    assert solved.n_iterations == 4
    assert solved.converged is True
    assert captured["solver_config"]["aot_gamma"] == 1.0
    assert captured["solve_args"][0] is toa
    assert rt_model == ("rt", config, "auth", mock_observation_bundle.sensor_config)


def test_resolve_output_writer_and_grid_assembler_use_runtime_entrypoints(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    class _FakeWriter:
        def __init__(self, defaults: object) -> None:
            captured["defaults"] = defaults

    def _fake_assemble(*args, **kwargs):  # noqa: ANN002, ANN003
        return args, kwargs

    monkeypatch.setattr(runtime_mod, "ConfiguredOutputWriter", _FakeWriter)
    monkeypatch.setattr(grid_assembler_mod, "assemble_grids", _fake_assemble)

    writer = runtime_mod.resolve_output_writer(
        SimpleNamespace(output=SimpleNamespace(defaults={"format": "netcdf"}))
    )
    assembler = runtime_mod.resolve_grid_assembler()

    assert isinstance(writer, _FakeWriter)
    assert captured["defaults"] == {"format": "netcdf"}
    assert assembler is _fake_assemble


def test_resample_helpers_and_corrector_cover_interpolation_zoom_and_passthrough(
    monkeypatch: pytest.MonkeyPatch,
    mock_observation_bundle,
    mock_atmospheric_state,
) -> None:
    template = mock_observation_bundle.toa["B02"]
    same = xr.DataArray(np.ones(template.shape, dtype=np.float32), dims=template.dims, coords=template.coords)
    assert runtime_mod._shares_template_grid(same, template) is True
    assert runtime_mod._resample_field_to_template(same, template) is same

    field = xr.DataArray(
        np.array([[0.0, 2.0], [10.0, 12.0]], dtype=np.float32),
        dims=["y", "x"],
        coords={"y": [1.0, 0.0], "x": [0.0, 2.0]},
    )
    interp_template = xr.DataArray(
        np.zeros((2, 2), dtype=np.float32),
        dims=["y", "x"],
        coords={"y": [1.0, 0.0], "x": [0.0, 1.0]},
    )
    interp = runtime_mod._resample_field_to_template(field, interp_template)
    assert float(interp.sel(y=1.0, x=1.0)) == pytest.approx(1.0)

    zoomed = runtime_mod._resample_field_to_template(
        xr.DataArray(np.array([[1.0, 2.0]], dtype=np.float32), dims=["y", "x"]),
        template,
    )
    assert zoomed.shape == template.shape
    latlon_template = xr.DataArray(
        np.zeros((2, 3), dtype=np.float32),
        dims=["latitude", "longitude"],
        coords={
            "latitude": [50.0, 49.5],
            "longitude": [-1.0, 0.0, 1.0],
        },
    )
    latlon_zoomed = runtime_mod._resample_field_to_template(
        xr.DataArray(np.array([[1.0, 2.0]], dtype=np.float32), dims=["y", "x"]),
        latlon_template,
    )
    assert latlon_zoomed.shape == latlon_template.shape
    assert latlon_zoomed.dims == latlon_template.dims
    assert latlon_zoomed.coords["latitude"].identical(latlon_template.coords["latitude"])
    assert latlon_zoomed.coords["longitude"].identical(latlon_template.coords["longitude"])
    empty = runtime_mod._resample_field_to_template(
        xr.DataArray(np.empty((0, 0), dtype=np.float32), dims=["y", "x"]),
        template,
    )
    assert np.isnan(empty.values).all()

    three_d = xr.DataArray(np.ones((1, 2, 2), dtype=np.float32), dims=["band", "y", "x"])
    assert runtime_mod._resample_field_to_template(three_d, template) is three_d

    captured: dict[str, object] = {}

    class _FakeCorrector:
        def __init__(self, rt_model: object, sensor_config: object) -> None:
            captured["init"] = (rt_model, sensor_config)

        def correct(self, toa, geometry, atmo, cloud_mask):  # noqa: ANN001
            captured["atmo_shapes"] = {
                "aot": atmo.aot.shape,
                "tcwv": atmo.tcwv.shape,
            }
            return "corrected"

    monkeypatch.setattr(runtime_mod, "AtmosphericCorrector", _FakeCorrector)

    small = xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"])
    solved = SolvedAtmosphere(
        atmo_state=type(mock_atmospheric_state)(
            aot=small,
            tcwv=small,
            tco3=small,
            aot_unc=small,
            tcwv_unc=small,
            tco3_unc=small,
            elevation=small,
        ),
        aot=small,
        tcwv=small,
        aot_unc=small,
        tcwv_unc=small,
        cost_final=0.0,
        n_iterations=1,
        converged=True,
    )

    result = runtime_mod.resolve_corrector(SimpleNamespace())(
        mock_observation_bundle,
        solved,
        "rt-model",
    )

    assert result == "corrected"
    assert captured["init"] == ("rt-model", mock_observation_bundle.sensor_config)
    assert captured["atmo_shapes"]["aot"] == template.shape


def test_resolve_corrector_preserves_coarse_atmo_outputs(
    monkeypatch: pytest.MonkeyPatch,
    mock_observation_bundle,
) -> None:
    import rioxarray  # noqa: F401

    from siac.runtime import AtmosphericState

    coarse = xr.DataArray(
        np.array([[0.2, 0.3], [0.4, 0.5]], dtype=np.float32),
        dims=["y", "x"],
        coords={"y": [1500.0, 500.0], "x": [500.0, 1500.0]},
    ).rio.write_crs("EPSG:32632")
    coarse_unc = xr.full_like(coarse, 0.05)
    solved = SolvedAtmosphere(
        atmo_state=AtmosphericState(
            aot=coarse,
            tcwv=coarse + 1.0,
            tco3=coarse + 0.1,
            aot_unc=coarse_unc,
            tcwv_unc=coarse_unc,
            tco3_unc=coarse_unc,
            elevation=xr.zeros_like(coarse),
        ),
        aot=coarse,
        tcwv=coarse + 1.0,
        aot_unc=coarse_unc,
        tcwv_unc=coarse_unc,
        cost_final=0.0,
        n_iterations=1,
        converged=True,
    )

    class _FakeCorrector:
        def __init__(self, rt_model: object, sensor_config: object) -> None:
            _ = (rt_model, sensor_config)

        def correct(self, toa, geometry, atmo, cloud_mask):  # noqa: ANN001
            from siac.runtime import CorrectionResult

            _ = (toa, geometry)
            return CorrectionResult(
                boa=mock_observation_bundle.toa,
                boa_unc=None,
                aot=atmo.aot,
                tcwv=atmo.tcwv,
                cloud_mask=cloud_mask,
            )

    monkeypatch.setattr(runtime_mod, "AtmosphericCorrector", _FakeCorrector)

    result = runtime_mod.resolve_corrector(SimpleNamespace())(
        mock_observation_bundle,
        solved,
        "rt-model",
    )

    assert result.aot.shape == (2, 2)
    assert result.tcwv.shape == (2, 2)
    assert result.aot.rio.crs is not None
    assert result.aot.rio.crs.to_string() == "EPSG:32632"
    assert result.aot.coords["x"].values.tolist() == [500.0, 1500.0]
    assert result.aot.coords["y"].values.tolist() == [1500.0, 500.0]


def test_resolve_corrector_fills_nonfinite_upsampled_atmo_for_lut_inputs(
    monkeypatch: pytest.MonkeyPatch,
    mock_observation_bundle,
) -> None:
    from siac.runtime import AtmosphericState

    coarse = xr.DataArray(
        np.array([[0.2, np.nan], [0.4, 0.5]], dtype=np.float32),
        dims=["y", "x"],
        coords={"y": [1500.0, 500.0], "x": [500.0, 1500.0]},
    )
    solved = SolvedAtmosphere(
        atmo_state=AtmosphericState(
            aot=coarse,
            tcwv=coarse + 1.0,
            tco3=xr.where(np.isfinite(coarse), coarse + 0.1, np.nan),
            aot_unc=xr.where(np.isfinite(coarse), 0.05, np.nan),
            tcwv_unc=xr.where(np.isfinite(coarse), 0.05, np.nan),
            tco3_unc=xr.where(np.isfinite(coarse), 0.05, np.nan),
            elevation=xr.zeros_like(coarse),
        ),
        aot=coarse,
        tcwv=coarse + 1.0,
        aot_unc=xr.where(np.isfinite(coarse), 0.05, np.nan),
        tcwv_unc=xr.where(np.isfinite(coarse), 0.05, np.nan),
        cost_final=0.0,
        n_iterations=1,
        converged=True,
    )
    captured: dict[str, np.ndarray] = {}

    class _FakeCorrector:
        def __init__(self, rt_model: object, sensor_config: object) -> None:
            _ = (rt_model, sensor_config)

        def correct(self, toa, geometry, atmo, cloud_mask):  # noqa: ANN001
            from siac.runtime import CorrectionResult

            _ = (toa, geometry, cloud_mask)
            captured["aot"] = np.asarray(atmo.aot.values, dtype=np.float32)
            captured["tcwv"] = np.asarray(atmo.tcwv.values, dtype=np.float32)
            captured["tco3"] = np.asarray(atmo.tco3.values, dtype=np.float32)
            return CorrectionResult(
                boa=mock_observation_bundle.toa,
                boa_unc=None,
                aot=atmo.aot,
                tcwv=atmo.tcwv,
                cloud_mask=mock_observation_bundle.cloud_mask,
            )

    monkeypatch.setattr(runtime_mod, "AtmosphericCorrector", _FakeCorrector)

    runtime_mod.resolve_corrector(SimpleNamespace())(
        mock_observation_bundle,
        solved,
        "rt-model",
    )

    assert np.isfinite(captured["aot"]).all()
    assert np.isfinite(captured["tcwv"]).all()
    assert np.isfinite(captured["tco3"]).all()


def test_shares_template_grid_false_cases_cover_mismatch_paths() -> None:
    template = xr.DataArray(
        np.ones((2, 2), dtype=np.float32),
        dims=["y", "x"],
        coords={"y": [1.0, 0.0], "x": [0.0, 1.0]},
    )

    assert runtime_mod._shares_template_grid(xr.DataArray(np.ones((1, 2), dtype=np.float32), dims=["y", "x"]), template) is False
    assert runtime_mod._shares_template_grid(xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["x", "y"]), template) is False

    missing_coords = xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"])
    assert runtime_mod._shares_template_grid(missing_coords, template) is False

    different_coords = xr.DataArray(
        np.ones((2, 2), dtype=np.float32),
        dims=["y", "x"],
        coords={"y": [0.0, 1.0], "x": [0.0, 1.0]},
    )
    assert runtime_mod._shares_template_grid(different_coords, template) is False

    coordless = xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"])
    assert runtime_mod._shares_template_grid(coordless, coordless) is True


def test_resample_field_to_template_falls_back_after_interp_error_and_pads_zoom(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    field = xr.DataArray(
        np.array([[1.0, 2.0]], dtype=np.float32),
        dims=["y", "x"],
        coords={"y": [0.0], "x": [0.0, 1.0]},
    )
    template = xr.DataArray(
        np.zeros((3, 3), dtype=np.float32),
        dims=["y", "x"],
        coords={"y": [2.0, 1.0, 0.0], "x": [0.0, 1.0, 2.0]},
    )

    def _raise_interp(_self: xr.DataArray, *args: object, **kwargs: object) -> xr.DataArray:  # noqa: ARG001
        raise RuntimeError("boom")

    def _short_zoom(src: np.ndarray, factors: tuple[float, float], order: int = 1) -> np.ndarray:  # noqa: ARG001
        return np.array([[5.0, 6.0]], dtype=np.float32)

    monkeypatch.setattr(xr.DataArray, "interp", _raise_interp)
    monkeypatch.setattr("scipy.ndimage.zoom", _short_zoom)

    out = runtime_mod._resample_field_to_template(field, template)

    assert out.shape == (3, 3)
    np.testing.assert_allclose(out.values[0, :2], np.array([5.0, 6.0], dtype=np.float32))
    assert np.isnan(out.values[0, 2])
    assert np.isnan(out.values[1:, :]).all()


def test_surface_helper_selection_and_mapping_runtime(
    monkeypatch: pytest.MonkeyPatch,
    mock_sensor_config,
) -> None:
    swir_sensor = SensorConfig(
        sensor_id="SWIR_ONLY",
        satellite_id="TEST",
        bands=(
            SensorBand("B11", 1610.0, 90.0, 20.0, 0),
            SensorBand("B12", 2190.0, 180.0, 20.0, 1),
        ),
        default_ref_scale=1.0,
        default_ref_offset=0.0,
    )
    msi_sensor = SensorConfig(
        sensor_id="MSI",
        satellite_id="S2",
        bands=(
            SensorBand("B08", 842.0, 115.0, 10.0, 0),
            SensorBand("B11", 1610.0, 90.0, 20.0, 1),
            SensorBand("B12", 2190.0, 180.0, 20.0, 2),
        ),
        default_ref_scale=1.0,
        default_ref_offset=0.0,
    )
    msi_visible_sensor = SensorConfig(
        sensor_id="MSI",
        satellite_id="S2",
        bands=(
            SensorBand("B01", 443.0, 20.0, 60.0, 0),
            SensorBand("B02", 490.0, 65.0, 10.0, 1),
            SensorBand("B03", 560.0, 35.0, 10.0, 2),
            SensorBand("B04", 665.0, 30.0, 10.0, 3),
            SensorBand("B08", 842.0, 115.0, 10.0, 4),
            SensorBand("B11", 1610.0, 90.0, 20.0, 5),
            SensorBand("B12", 2190.0, 180.0, 20.0, 6),
        ),
        default_ref_scale=1.0,
        default_ref_offset=0.0,
    )

    assert surface_mod._select_surface_prior_bands(None) == list(range(1, 8))
    assert [band.name for band in surface_mod._select_surface_prior_bands(mock_sensor_config)] == ["B01", "B02"]
    assert [band.name for band in surface_mod._select_visible_surface_prior_bands(mock_sensor_config)] == ["B01", "B02", "B03", "B04"]
    assert [band.name for band in surface_mod._select_visible_surface_prior_bands(msi_visible_sensor)] == ["B01", "B02", "B04"]
    assert [band.name for band in surface_mod._select_visible_surface_prior_bands(swir_sensor)] == ["B11", "B12"]
    assert [band.name for band in surface_mod._select_route_b_query_bands(msi_sensor)] == ["B08", "B11", "B12"]
    assert len(surface_mod._select_route_b_query_bands(mock_sensor_config)) == 3

    config = SimpleNamespace(
        surface_prior=SimpleNamespace(
            spectral_mapping=SimpleNamespace(
                enabled=False,
                k_neighbors=4,
                neighbor_estimator="distance_weighted_mean",
                knn_backend="numpy",
                knn_eps=0.0,
                min_valid_bands=1,
                siac_library_root=None,
                rsrf_root=None,
                cache_dir=None,
            )
        ),
        paths=SimpleNamespace(
            spectral_library_root="/tmp/library",
            rsrf_root="/tmp/rsrf",
            caches=SimpleNamespace(spectral_mapping="/tmp/cache"),
        ),
    )
    assert surface_mod._surface_spectral_mapping_runtime(config) == (None, 4)

    monkeypatch.setattr(
        "siac.algorithms.surface.spectral_mapping.needs_spectral_mapping",
        lambda source, target: source != target,
    )
    with pytest.raises(ValueError, match="enabled=false"):
        surface_mod._surface_spectral_mapping_runtime(
            config,
            source_bands=(mock_sensor_config.bands[0],),
            target_bands=(mock_sensor_config.bands[1],),
            context="tests",
        )

    class _RuntimeConfig:
        def __init__(self, **kwargs) -> None:
            self.kwargs = kwargs

    config.surface_prior.spectral_mapping.enabled = True
    monkeypatch.setattr(
        "siac.algorithms.surface.spectral_mapping.SpectralMappingConfig",
        _RuntimeConfig,
    )
    spectral_config, neighbors = surface_mod._surface_spectral_mapping_runtime(
        config,
        source_bands=(mock_sensor_config.bands[0],),
        target_bands=(mock_sensor_config.bands[0],),
    )
    assert neighbors == 4
    assert spectral_config.kwargs["siac_library_root"] == "/tmp/library"
    assert spectral_config.kwargs["rsrf_root"] == "/tmp/rsrf"
    assert spectral_config.kwargs["cache_dir"] == "/tmp/cache"


def test_surface_monthly_runtime_and_query_helpers(
    monkeypatch: pytest.MonkeyPatch,
    mock_observation_bundle,
    mock_atmospheric_state,
) -> None:
    captured: dict[str, object] = {}

    def _fake_resample_geometry(observation, *, resolution):  # noqa: ANN001
        captured["resample"] = (observation, resolution)
        return "geometry"

    def _fake_build_database(**kwargs):  # noqa: ANN003
        captured["database"] = kwargs
        return "db"

    def _fake_query_database(**kwargs):  # noqa: ANN003
        captured["query"] = kwargs
        return "prior"

    monkeypatch.setattr(
        "siac.algorithms.surface.swir_refine.resample_geometry_for_surface_prior",
        _fake_resample_geometry,
    )
    monkeypatch.setattr(
        "siac.algorithms.surface.swir_refine.build_monthly_surface_prior_database",
        _fake_build_database,
    )
    monkeypatch.setattr(
        "siac.algorithms.surface.swir_refine.query_surface_prior_from_monthly_database",
        _fake_query_database,
    )

    provider = SimpleNamespace(source_bands=list(mock_observation_bundle.sensor_config.bands))
    config = SimpleNamespace(
        surface_prior=SimpleNamespace(
            spectral_mapping=SimpleNamespace(
                enabled=True,
                k_neighbors=7,
                neighbor_estimator="distance_weighted_mean",
                knn_backend="numpy",
                knn_eps=0.0,
                min_valid_bands=1,
                siac_library_root=None,
                rsrf_root=None,
                cache_dir=None,
            )
        ),
        brdf=SimpleNamespace(temporal_window=16),
        paths=None,
    )

    runtime = surface_mod._prepare_monthly_surface_prior_runtime(
        config,
        provider,
        observation=mock_observation_bundle,
        resolution=500.0,
    )
    result = surface_mod._query_monthly_surface_prior(
        mock_observation_bundle,
        mock_atmospheric_state,
        "rt-model",
        runtime,
    )
    request = surface_mod._surface_prior_brdf_request(
        mock_observation_bundle,
        brdf_provider=provider,
        target_resolution=500.0,
        temporal_window=16,
    )
    def provider_fn(*args, **kwargs):  # noqa: ANN002, ANN003
        return args, kwargs

    assert runtime.database == "db"
    assert runtime.spectral_k_neighbors == 7
    assert captured["database"]["geometry"] == "geometry"
    assert captured["query"]["database"] == "db"
    assert captured["query"]["k_neighbors"] == 7
    assert result == "prior"
    assert request["obs_time"] == mock_observation_bundle.metadata["observation_time"]
    assert surface_mod._mark_surface_prior_metadata(provider_fn, requires_atmo_prior=True) is provider_fn
    assert provider_fn.requires_atmo_prior is True


def test_provider_builders_cover_registry_and_source_resolution(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    class _FakeCAMSProvider:
        def __init__(self, *args, **kwargs) -> None:  # noqa: ANN002, ANN003
            captured["cams"] = (args, kwargs)

    class _FakePriorProvider:
        def __init__(self, **kwargs) -> None:
            self.kwargs = kwargs

        def get_prior(self, *args, **kwargs):  # noqa: ANN002, ANN003
            return (args, kwargs, self.kwargs)

    monkeypatch.setattr(providers_mod, "CAMSProvider", _FakeCAMSProvider)
    monkeypatch.setattr(providers_mod, "earthaccess_source_from_auth", lambda auth: f"source:{auth}")
    monkeypatch.setattr("siac.adapters.atmo.merra2.MERRA2Provider", _FakePriorProvider)
    monkeypatch.setattr("siac.adapters.atmo.mcd19_earthaccess.MCD19AODProvider", _FakePriorProvider)
    monkeypatch.setattr("siac.adapters.atmo.mcd19_earthaccess.VNP19AODProvider", _FakePriorProvider)
    monkeypatch.setattr("siac.adapters.brdf.mcd43_earthaccess.MCD43EarthAccessProvider", _FakePriorProvider)
    monkeypatch.setattr("siac.adapters.brdf.mcd43_earthaccess.MCD19EarthAccessProvider", _FakePriorProvider)
    monkeypatch.setattr("siac.adapters.brdf.gee_stub.GEEBRDFProvider", lambda: "gee")

    config = SimpleNamespace(
        atmo_prior=SimpleNamespace(
            provider="merra2",
            data_path="/tmp/cams.nc",
            temporal_interpolation="linear",
            download_missing=True,
            cache_dir="/tmp/cache",
        ),
        brdf=SimpleNamespace(provider="gee", cache_dir="/tmp/cache"),
    )

    providers_mod._build_cams_provider(config, auth="token")
    merra = providers_mod._build_merra2_provider(config, auth="earth")
    mcd19 = providers_mod._build_mcd19_provider(config, auth="earth")
    vnp19 = providers_mod._build_vnp19_provider(config, auth="earth")
    gee = providers_mod._build_gee_brdf_provider(config)
    brdf = providers_mod.resolve_brdf_provider(config, auth="earth")
    get_prior = providers_mod.resolve_atmo_provider(config, auth="earth")

    assert captured["cams"][0][0] == "/tmp/cams.nc"
    assert captured["cams"][1]["auth"] == "token"
    assert merra.kwargs["source"] == "source:earth"
    assert mcd19.kwargs["source"] == "source:earth"
    assert vnp19.kwargs["source"] == "source:earth"
    assert gee == "gee"
    assert brdf == "gee"
    assert callable(get_prior)


def test_build_registered_component_delegates_common_signature() -> None:
    registry = SimpleNamespace(get=lambda name: lambda config, auth: (name, config, auth))
    config = SimpleNamespace()

    assert providers_mod._build_registered_component(registry, "x", config, auth="y") == ("x", config, "y")


def test_prepared_store_monthly_provider_loads_collection(tmp_path) -> None:
    bands = (
        SensorBand("B02", 490.0, 65.0, 10.0, 0),
        SensorBand("B08", 842.0, 115.0, 10.0, 1),
    )
    coords = {"band": ["B02", "B08"], "y": [0, 1], "x": [0]}
    cube = xr.DataArray(np.full((2, 2, 1), 0.2, dtype=np.float32), dims=["band", "y", "x"], coords=coords)
    collection = MonthlyCompositeCollection(
        composites=(
            MonthlyKernelWeightComposite(
                kernels=BRDFKernelWeights(
                    f0=cube,
                    f1=xr.zeros_like(cube),
                    f2=xr.zeros_like(cube),
                    f0_unc=xr.full_like(cube, 0.01),
                    f1_unc=xr.full_like(cube, 0.01),
                    f2_unc=xr.full_like(cube, 0.01),
                ),
                quality=xr.DataArray(np.full((2, 1), 0.03, dtype=np.float32), dims=["y", "x"], coords={"y": [0, 1], "x": [0]}),
                sample_index=xr.DataArray(np.zeros((2, 1), dtype=np.int16), dims=["y", "x"], coords={"y": [0, 1], "x": [0]}),
                year=2023,
                month=7,
            ),
        ),
        source_bands=bands,
        source_name="prepared-test",
    )
    store_path = write_monthly_composite_collection(collection, tmp_path / "prepared_store")
    config = SimpleNamespace(
        monthly_composites=SimpleNamespace(
            provider="prepared_store",
            store_path=store_path,
            strict_coverage=False,
        )
    )

    provider = providers_mod.resolve_monthly_composite_provider(config)
    loaded = provider.get_monthly_composites("observation", 500.0)

    assert provider.source_name == "prepared-test"
    assert [band.name for band in provider.source_bands] == ["B02", "B08"]
    assert isinstance(loaded.composites[0], MonthlyKernelWeightComposite)


def test_prepared_store_monthly_provider_is_lazy_until_first_use(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    calls = {"manifest": 0, "collection": 0}
    manifest = MonthlyCompositeStoreManifest(
        version=2,
        source_name="prepared-test",
        source_bands=(),
        entries=(),
    )
    collection = MonthlyCompositeCollection(composites=(), source_bands=(), source_name="prepared-test")

    def _fake_read_manifest(_path):  # noqa: ANN001
        calls["manifest"] += 1
        return manifest

    def _fake_read_collection(_path):  # noqa: ANN001
        calls["collection"] += 1
        return collection

    monkeypatch.setattr(
        "siac.algorithms.surface.monthly_composite_store.read_monthly_composite_store_manifest",
        _fake_read_manifest,
    )
    monkeypatch.setattr(
        "siac.algorithms.surface.monthly_composite_store.read_monthly_composite_collection",
        _fake_read_collection,
    )

    config = SimpleNamespace(
        monthly_composites=SimpleNamespace(
            provider="prepared_store",
            store_path=tmp_path / "prepared_store",
            strict_coverage=False,
        )
    )

    provider = providers_mod.resolve_monthly_composite_provider(config)

    assert calls == {"manifest": 0, "collection": 0}
    assert provider.source_name == "prepared-test"
    assert calls == {"manifest": 1, "collection": 0}
    assert provider.get_monthly_composites("observation", 500.0) is collection
    assert calls == {"manifest": 1, "collection": 1}


def test_prepared_store_monthly_provider_rejects_mismatched_grid_when_strict(tmp_path) -> None:
    bands = (
        SensorBand("B02", 490.0, 65.0, 10.0, 0),
        SensorBand("B08", 842.0, 115.0, 10.0, 1),
    )
    coords = {"band": ["B02", "B08"], "y": [750.0, 250.0], "x": [250.0, 750.0]}
    cube = xr.DataArray(np.full((2, 2, 2), 0.2, dtype=np.float32), dims=["band", "y", "x"], coords=coords)
    collection = MonthlyCompositeCollection(
        composites=(
            MonthlyKernelWeightComposite(
                kernels=BRDFKernelWeights(
                    f0=cube,
                    f1=xr.zeros_like(cube),
                    f2=xr.zeros_like(cube),
                    f0_unc=xr.full_like(cube, 0.01),
                    f1_unc=xr.full_like(cube, 0.01),
                    f2_unc=xr.full_like(cube, 0.01),
                ),
                quality=xr.DataArray(
                    np.full((2, 2), 0.03, dtype=np.float32),
                    dims=["y", "x"],
                    coords={"y": [750.0, 250.0], "x": [250.0, 750.0]},
                ),
                sample_index=xr.DataArray(
                    np.zeros((2, 2), dtype=np.int16),
                    dims=["y", "x"],
                    coords={"y": [750.0, 250.0], "x": [250.0, 750.0]},
                ),
                year=2023,
                month=7,
            ),
        ),
        source_bands=bands,
        source_name="prepared-test",
    )
    store_path = write_monthly_composite_collection(
        collection,
        tmp_path / "prepared_store",
        grid=MonthlyCompositeStoreGridSpec.from_bounds((0.0, 0.0, 1000.0, 1000.0), crs="EPSG:32632", resolution=500.0),
    )
    config = SimpleNamespace(
        monthly_composites=SimpleNamespace(
            provider="prepared_store",
            store_path=store_path,
            strict_coverage=True,
        )
    )

    provider = providers_mod.resolve_monthly_composite_provider(config)
    observation = SimpleNamespace(bounds=(0.0, 0.0, 2000.0, 1000.0), crs="EPSG:32632")

    with pytest.raises(ValueError, match="uses shape"):
        provider.get_monthly_composites(observation, 500.0)
