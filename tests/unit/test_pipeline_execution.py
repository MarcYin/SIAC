"""Unit tests for pipeline execution settings, retries, and backends."""

from __future__ import annotations

import dataclasses
import json
import sys
import time
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any

import numpy as np
import pytest
import xarray as xr

import siac.workflows.pipeline as pipeline
from siac.runtime import AtmosphericState, CorrectionResult


def _install_fake_dask(monkeypatch: pytest.MonkeyPatch, *, mode: str = "success") -> Any:
    class FakeDaskTimeoutError(Exception):
        pass

    class FakeFuture:
        def __init__(self, *, value: Any = None, error: Exception | None = None):
            self._value = value
            self._error = error
            self.canceled = False

        def result(self, timeout: float | None = None) -> Any:
            _ = timeout
            if self._error is not None:
                raise self._error
            return self._value

        def cancel(self) -> None:
            self.canceled = True

    class FakeLocalCluster:
        last_kwargs: dict[str, Any] | None = None

        def __init__(self, **kwargs: Any):
            self.kwargs = kwargs
            FakeLocalCluster.last_kwargs = kwargs

        def __enter__(self) -> FakeLocalCluster:
            return self

        def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> bool:
            _ = (exc_type, exc, tb)
            return False

    class FakePerformanceReport:
        last_filename: str | None = None
        entered = False
        exited = False

        def __init__(self, filename: str):
            self.filename = filename
            FakePerformanceReport.last_filename = filename

        def __enter__(self) -> FakePerformanceReport:
            FakePerformanceReport.entered = True
            return self

        def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> bool:
            _ = (exc_type, exc, tb)
            FakePerformanceReport.exited = True
            return False

    class FakeClient:
        last_instance: FakeClient | None = None
        mode: str = "success"

        def __init__(self, cluster: FakeLocalCluster):
            _ = cluster
            self.dashboard_link = "http://fake-dashboard:8787/status"
            self.futures: list[FakeFuture] = []
            self._submit_count = 0
            FakeClient.last_instance = self

        def __enter__(self) -> FakeClient:
            return self

        def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> bool:
            _ = (exc_type, exc, tb)
            return False

        def submit(self, fn: Any, *args: Any, **kwargs: Any) -> FakeFuture:
            if self.mode == "timeout" and self._submit_count == 0:
                fut = FakeFuture(error=FakeDaskTimeoutError("timed out"))
                self._submit_count += 1
                self.futures.append(fut)
                return fut

            if self.mode == "error" and self._submit_count == 0:
                fut = FakeFuture(error=RuntimeError("future failure"))
                self._submit_count += 1
                self.futures.append(fut)
                return fut

            try:
                value = fn(*args, **kwargs)
                fut = FakeFuture(value=value)
            except Exception as exc:  # pragma: no cover - defensive path
                fut = FakeFuture(error=exc)
            self._submit_count += 1
            self.futures.append(fut)
            return fut

    def _performance_report(*, filename: str) -> FakePerformanceReport:
        return FakePerformanceReport(filename)

    FakeClient.mode = mode

    fake_distributed = ModuleType("dask.distributed")
    fake_distributed.Client = FakeClient  # type: ignore[attr-defined]
    fake_distributed.LocalCluster = FakeLocalCluster  # type: ignore[attr-defined]
    fake_distributed.TimeoutError = FakeDaskTimeoutError  # type: ignore[attr-defined]
    fake_distributed.performance_report = _performance_report  # type: ignore[attr-defined]

    fake_dask = ModuleType("dask")
    fake_dask.distributed = fake_distributed  # type: ignore[attr-defined]

    monkeypatch.setitem(sys.modules, "dask", fake_dask)
    monkeypatch.setitem(sys.modules, "dask.distributed", fake_distributed)

    return SimpleNamespace(
        client_cls=FakeClient,
        cluster_cls=FakeLocalCluster,
        report_cls=FakePerformanceReport,
    )


def test_execution_values_handles_dict_and_object() -> None:
    values = pipeline._execution_values({"backend": "dask", "max_workers": 3, "extra": 1})
    assert values == {"backend": "dask", "max_workers": 3}

    obj = SimpleNamespace(backend="thread", retries=4, ignored=True)
    assert pipeline._execution_values(obj) == {"backend": "thread", "retries": 4}


def test_resolve_execution_settings_applies_overrides_and_coercions(tmp_path: Path) -> None:
    cfg = SimpleNamespace(
        runtime=SimpleNamespace(
            execution=SimpleNamespace(
                backend="THREAD",
                max_workers=2,
                retries=1,
                stage_timeout_s="15",
                stage_timeouts={"M2.atmospheric_prior": "8"},
                dashboard=1,
                dashboard_address=":8787",
                performance_report_path=str(tmp_path / "reports" / "perf.html"),
                show_progress=1,
                profiling_sample_interval_s="2.5",
                progress_heartbeat_s="12",
            )
        )
    )

    settings = pipeline._resolve_execution_settings(
        cfg,
        execution={"retries": 3, "stage_timeouts": {"M3.surface_prior": "7"}},
        max_workers=6,
    )

    assert settings["backend"] == "thread"
    assert settings["max_workers"] == 6
    assert settings["retries"] == 3
    assert settings["stage_timeout_s"] == 15.0
    assert settings["stage_timeouts"] == {"M3.surface_prior": 7.0}
    assert settings["dashboard"] is True
    assert settings["show_progress"] is True
    assert isinstance(settings["performance_report_path"], Path)
    assert settings["profiling_sample_interval_s"] == 2.5
    assert settings["progress_heartbeat_s"] == 12.0


@pytest.mark.parametrize(
    ("execution", "message"),
    [
        ({"backend": "unknown"}, "Unsupported execution backend"),
        ({"max_workers": 0}, "max_workers must be >= 1"),
        ({"retries": -1}, "retries must be >= 0"),
        ({"stage_timeout_s": 0}, "stage_timeout_s must be > 0"),
        ({"stage_timeouts": []}, "stage_timeouts must be a mapping"),
        ({"stage_timeouts": {"M2": 0}}, "stage_timeouts values must be > 0"),
        ({"profiling_sample_interval_s": 0}, "profiling_sample_interval_s must be > 0"),
        ({"progress_heartbeat_s": 0}, "progress_heartbeat_s must be > 0"),
    ],
)
def test_resolve_execution_settings_rejects_invalid_values(
    execution: dict[str, Any],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        pipeline._resolve_execution_settings(
            SimpleNamespace(), execution=execution, max_workers=None
        )


def test_aerosol_resolution_reads_solver_config() -> None:
    cfg = SimpleNamespace(
        algorithms=SimpleNamespace(solver=SimpleNamespace(aerosol_resolution=120.0))
    )
    assert pipeline._aerosol_resolution(cfg) == 120.0


def test_run_tail_passes_configured_aerosol_resolution_to_grid_assembler(
    mock_observation_bundle,
    mock_atmospheric_state,
    mock_surface_prior,
    mock_solver_input_bundle,
    mock_solver_fn,
    mock_corrector_fn,
    mock_rt_model,
) -> None:
    captured: dict[str, Any] = {}

    def _grid_assembler(obs, atmo, surface, rt_model, **kwargs):  # noqa: ANN001
        captured["args"] = (obs, atmo, surface, rt_model)
        captured["kwargs"] = kwargs
        return mock_solver_input_bundle

    cfg = SimpleNamespace(
        algorithms=SimpleNamespace(
            solver=SimpleNamespace(
                aerosol_resolution=120.0,
                water_mask_buffer_pixels=3,
                stages=(
                    SimpleNamespace(name="aot_pass", bands=("B01", "B02", "B04")),
                    SimpleNamespace(name="tcwv_pass", bands=("B02", "B04")),
                ),
            )
        )
    )

    result = pipeline._run_tail(
        mock_observation_bundle,
        mock_atmospheric_state,
        mock_surface_prior,
        cfg,
        grid_assembler=_grid_assembler,
        solver=mock_solver_fn,
        corrector=mock_corrector_fn,
        rt_model=mock_rt_model,
    )

    assert isinstance(result, CorrectionResult)
    assert captured["args"][0] is mock_observation_bundle
    assert captured["kwargs"]["aerosol_resolution_m"] == 120.0
    assert captured["kwargs"]["water_mask_path"] == pipeline.DEFAULT_WATER_MASK_VRT_URL
    assert captured["kwargs"]["water_mask_buffer_pixels"] == 3
    assert captured["kwargs"]["solver_band_names"] == ("B01", "B02", "B04")


def test_call_grid_assembler_requires_standardized_kwargs(
    mock_observation_bundle,
    mock_atmospheric_state,
    mock_surface_prior,
    mock_rt_model,
) -> None:
    calls: list[dict[str, Any]] = []

    class _OpaqueAssembler:
        @property
        def __signature__(self):  # noqa: ANN204
            raise ValueError("signature unavailable")

        def __call__(self, obs, atmo, surface, rt_model, *, aerosol_resolution_m):  # noqa: ANN001
            calls.append(
                {
                    "obs": obs,
                    "atmo": atmo,
                    "surface": surface,
                    "rt_model": rt_model,
                    "aerosol_resolution_m": aerosol_resolution_m,
                }
            )
            return "solver-inputs"

    with pytest.raises(TypeError, match="sharp_transition_filter"):
        pipeline._call_grid_assembler(
            _OpaqueAssembler(),
            mock_observation_bundle,
            mock_atmospheric_state,
            mock_surface_prior,
            mock_rt_model,
            aerosol_resolution_m=120.0,
            sharp_transition_filter=SimpleNamespace(enabled=True),
            water_mask_path="https://example.com/water-mask.vrt",
            water_mask_cache_dir=Path("/tmp/water-mask-cache"),
            water_mask_buffer_pixels=2,
            solver_band_names=("B01", "B02", "B04"),
        )

    assert calls == []


def test_select_solver_bands_for_preload_includes_stage_requested_bands(mock_sensor_config) -> None:
    bands = pipeline._select_solver_bands_for_preload(
        mock_sensor_config,
        requested_band_names=("B04", "B01"),
    )

    assert [band.name for band in bands] == ["B01", "B02", "B04"]


def test_run_pipeline_lut_preload_includes_stage_requested_bands(
    mock_preprocessor,
    mock_surface_prior_provider,
    mock_grid_assembler,
    mock_solver_fn,
    mock_corrector_fn,
) -> None:
    coarse_coords = {
        "y": np.arange(4, dtype=np.float32),
        "x": np.arange(4, dtype=np.float32),
    }

    def _coarse_atmo(bounds: Any, crs: str, obs_time: Any, res: float) -> AtmosphericState:
        _ = (bounds, crs, obs_time, res)
        aot = xr.DataArray(
            np.full((4, 4), 0.2, dtype=np.float32), dims=("y", "x"), coords=coarse_coords
        )
        tcwv = xr.DataArray(
            np.full((4, 4), 2.0, dtype=np.float32), dims=("y", "x"), coords=coarse_coords
        )
        tco3 = xr.DataArray(
            np.full((4, 4), 0.3, dtype=np.float32), dims=("y", "x"), coords=coarse_coords
        )
        unc = xr.DataArray(
            np.full((4, 4), 0.01, dtype=np.float32), dims=("y", "x"), coords=coarse_coords
        )
        elevation = xr.DataArray(
            np.zeros((4, 4), dtype=np.float32), dims=("y", "x"), coords=coarse_coords
        )
        return AtmosphericState(
            aot=aot,
            tcwv=tcwv,
            tco3=tco3,
            aot_unc=unc,
            tcwv_unc=unc,
            tco3_unc=unc,
            elevation=elevation,
        )

    captured: dict[str, Any] = {}

    class _PreloadRTModel:
        def preload_scene_subset(
            self, geometry: Any, atmo_state: AtmosphericState, bands: list[Any]
        ) -> None:
            captured["geometry_shape"] = tuple(geometry.sza.shape)
            captured["atmo_shape"] = tuple(atmo_state.aot.shape)
            captured["bands"] = [band.name for band in bands]

    result = pipeline.run_pipeline(
        Path("/fake"),
        None,
        SimpleNamespace(
            algorithms=SimpleNamespace(
                solver=SimpleNamespace(
                    stages=(SimpleNamespace(name="aot_pass", bands=("B01", "B02", "B04")),),
                )
            )
        ),
        preprocessor=mock_preprocessor,
        atmo_provider=_coarse_atmo,
        surface_prior_provider=mock_surface_prior_provider,
        grid_assembler=mock_grid_assembler,
        solver=mock_solver_fn,
        corrector=mock_corrector_fn,
        rt_model=_PreloadRTModel(),
        execution={"backend": "thread", "retries": 0},
    )

    assert isinstance(result, CorrectionResult)
    assert captured["geometry_shape"] == (4, 4)
    assert captured["atmo_shape"] == (4, 4)
    assert captured["bands"] == ["B01", "B02", "B04"]


def test_run_tail_attaches_surface_prior_and_monthly_composite_outputs(
    mock_observation_bundle,
    mock_atmospheric_state,
    mock_surface_prior,
    mock_solver_input_bundle,
    mock_solver_fn,
    mock_corrector_fn,
    mock_rt_model,
) -> None:
    composite = SimpleNamespace(
        year=2023,
        month=7,
        reflectance=xr.DataArray(
            np.stack([np.full(mock_surface_prior.boa.shape, 0.14, dtype=np.float32)]),
            dims=["band", "y", "x"],
            coords={"band": ["B02"]},
        ),
        quality=xr.DataArray(
            np.full(mock_surface_prior.boa.shape, 0.5, dtype=np.float32), dims=["y", "x"]
        ),
        sample_index=xr.DataArray(
            np.full(mock_surface_prior.boa.shape, 2, dtype=np.int16), dims=["y", "x"]
        ),
    )
    surface_with_monthly = dataclasses.replace(
        mock_surface_prior,
        monthly_composites=(composite,),
    )

    result = pipeline._run_tail(
        mock_observation_bundle,
        mock_atmospheric_state,
        surface_with_monthly,
        SimpleNamespace(
            algorithms=SimpleNamespace(solver=SimpleNamespace(aerosol_resolution=1000.0))
        ),
        grid_assembler=lambda *_args, **_kwargs: mock_solver_input_bundle,
        solver=mock_solver_fn,
        corrector=mock_corrector_fn,
        rt_model=mock_rt_model,
    )

    assert result.surface_prior is not None
    assert result.surface_prior_unc is not None
    assert "surface_prior" in result.surface_prior.data_vars
    assert result.monthly_composites is not None
    assert "2023_07" in result.monthly_composites
    assert "B02" in result.monthly_composites["2023_07"].reflectance.data_vars


def test_run_tail_attaches_aot_scatter_diagnostics(
    mock_observation_bundle,
    mock_atmospheric_state,
    mock_surface_prior,
    mock_solver_input_bundle,
    mock_solver_fn,
    mock_corrector_fn,
    mock_rt_model,
) -> None:
    result = pipeline._run_tail(
        mock_observation_bundle,
        mock_atmospheric_state,
        mock_surface_prior,
        SimpleNamespace(
            algorithms=SimpleNamespace(solver=SimpleNamespace(aerosol_resolution=1000.0))
        ),
        grid_assembler=lambda *_args, **_kwargs: mock_solver_input_bundle,
        solver=mock_solver_fn,
        corrector=mock_corrector_fn,
        rt_model=mock_rt_model,
    )

    assert result.diagnostics.aot_scatter_plots
    scatter = result.diagnostics.aot_scatter_plots[0]
    assert scatter.band_name == "B02"
    assert scatter.total_valid_count > 0
    assert (
        scatter.surface_reflectance.shape
        == scatter.observed_toa.shape
        == scatter.simulated_toa.shape
    )


def test_run_tail_attaches_solver_quality_metadata(
    mock_observation_bundle,
    mock_atmospheric_state,
    mock_surface_prior,
    mock_solver_input_bundle,
    mock_solver_fn,
    mock_corrector_fn,
    mock_rt_model,
) -> None:
    result = pipeline._run_tail(
        mock_observation_bundle,
        mock_atmospheric_state,
        mock_surface_prior,
        SimpleNamespace(
            algorithms=SimpleNamespace(solver=SimpleNamespace(aerosol_resolution=1000.0))
        ),
        grid_assembler=lambda *_args, **_kwargs: mock_solver_input_bundle,
        solver=mock_solver_fn,
        corrector=mock_corrector_fn,
        rt_model=mock_rt_model,
    )

    solver = result.metadata["solver"]
    assert solver["cost_final"] == pytest.approx(0.001)
    assert solver["n_iterations"] == 5
    assert solver["converged"] is True
    assert solver["solve_band_count"] == len(mock_solver_input_bundle.bands)
    assert solver["cost_final_per_band"] == pytest.approx(
        0.001 / len(mock_solver_input_bundle.bands)
    )
    assert solver["aot_finite_fraction"] == pytest.approx(1.0)
    assert solver["aot_median"] == pytest.approx(0.15)
    assert solver["aot_unc_median"] == pytest.approx(0.05)


def test_aot_scatter_diagnostics_fill_nan_atmo_for_lut_call(
    mock_solver_input_bundle,
    mock_solved_atmosphere,
) -> None:
    from siac.runtime import RTCoefficients

    nan_aot = mock_solved_atmosphere.atmo_state.aot.copy()
    nan_aot.values[0, 0] = np.nan
    atmo_state = dataclasses.replace(mock_solved_atmosphere.atmo_state, aot=nan_aot)
    solved = dataclasses.replace(mock_solved_atmosphere, atmo_state=atmo_state, aot=nan_aot)
    calls = {"n": 0}

    class _FiniteRT:
        def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
            _ = (band, compute_jacobian)
            calls["n"] += 1
            assert np.isfinite(atmo_state.aot.values).all()
            xap = xr.full_like(geometry.sza, 0.95, dtype=np.float32)
            return RTCoefficients(xap=xap, xbp=xr.zeros_like(xap), xcp=xr.zeros_like(xap))

    solver_inputs = dataclasses.replace(mock_solver_input_bundle, rt_model=_FiniteRT())

    scatter = pipeline._build_aot_scatter_diagnostics(solver_inputs, solved)

    assert calls["n"] > 0
    assert scatter
    assert scatter[0].total_valid_count == int(np.count_nonzero(np.isfinite(nan_aot.values)))


def test_call_with_retries_recovers_and_logs(caplog: pytest.LogCaptureFixture) -> None:
    calls = {"n": 0}

    def flaky() -> str:
        calls["n"] += 1
        if calls["n"] == 1:
            raise RuntimeError("temporary")
        return "ok"

    with caplog.at_level("WARNING"):
        out = pipeline._call_with_retries(flaky, (), retries=1, stage_name="M2")

    assert out == "ok"
    assert calls["n"] == 2
    assert "retrying" in caplog.text


def test_call_with_retries_raises_after_exhaustion() -> None:
    def always_fail() -> None:
        raise ValueError("permanent")

    with pytest.raises(ValueError, match="permanent"):
        pipeline._call_with_retries(always_fail, (), retries=0, stage_name="M3")


def test_call_with_retries_does_not_retry_deterministic_errors(
    caplog: pytest.LogCaptureFixture,
) -> None:
    calls = {"n": 0}

    def invalid() -> None:
        calls["n"] += 1
        raise ValueError("invalid input")

    with caplog.at_level("WARNING"), pytest.raises(ValueError, match="invalid input"):
        pipeline._call_with_retries(invalid, (), retries=3, stage_name="M3")

    assert calls["n"] == 1
    assert "retrying" not in caplog.text


def test_fetch_priors_uses_stage_specific_timeouts(
    mock_observation_bundle,
    mock_atmospheric_state,
    mock_surface_prior,
) -> None:
    class _Future:
        def __init__(self, value: Any):
            self.value = value
            self.timeouts: list[float | None] = []
            self.canceled = False

        def result(self, timeout: float | None = None) -> Any:
            self.timeouts.append(timeout)
            return self.value

        def cancel(self) -> None:
            self.canceled = True

    futures: dict[str, _Future] = {}

    def _submit(fn: Any, *args: Any, **kwargs: Any) -> _Future:  # noqa: ARG001
        stage_name = kwargs["stage_name"]
        value = {
            "M2.atmospheric_prior": mock_atmospheric_state,
            "M3.surface_prior": mock_surface_prior,
            "LUT.preload": None,
        }[stage_name]
        future = _Future(value)
        futures[stage_name] = future
        return future

    class _RTModel:
        def preload_scene_subset(
            self, geometry: Any, atmo_state: AtmosphericState, bands: list[Any]
        ) -> None:
            _ = (geometry, atmo_state, bands)

    atmo, surface = pipeline._fetch_priors(
        submit_fn=_submit,
        lut_submit_fn=_submit,
        obs=mock_observation_bundle,
        config=SimpleNamespace(),
        atmo_provider=lambda *_args: mock_atmospheric_state,
        surface_prior_provider=lambda *_args: mock_surface_prior,
        rt_model=_RTModel(),
        settings={
            "retries": 0,
            "stage_timeout_s": 99.0,
            "stage_timeouts": {
                "M2.atmospheric_prior": 1.5,
                "M3.surface_prior": 2.5,
                "LUT.preload": 3.5,
            },
        },
        backend_label="thread",
    )

    assert atmo is mock_atmospheric_state
    assert surface is mock_surface_prior
    assert futures["M2.atmospheric_prior"].timeouts == [1.5]
    assert futures["M3.surface_prior"].timeouts == [2.5]
    assert futures["LUT.preload"].timeouts == [3.5]


def test_run_pipeline_thread_timeout_raises(
    mock_preprocessor,
    mock_atmospheric_state,
    mock_surface_prior,
    mock_grid_assembler,
    mock_solver_fn,
    mock_corrector_fn,
    mock_rt_model,
) -> None:
    def slow_atmo(bounds: Any, crs: str, obs_time: Any, res: float) -> Any:
        _ = (bounds, crs, obs_time, res)
        time.sleep(0.05)
        return mock_atmospheric_state

    def slow_surface(obs: Any, atmo_prior: Any, rt_model: Any, res: float) -> Any:
        _ = (obs, atmo_prior, rt_model, res)
        time.sleep(0.05)
        return mock_surface_prior

    with pytest.raises(TimeoutError, match="thread backend"):
        pipeline._run_pipeline_thread(
            Path("/fake"),
            None,
            None,
            preprocessor=mock_preprocessor,
            atmo_provider=slow_atmo,
            surface_prior_provider=slow_surface,
            grid_assembler=mock_grid_assembler,
            solver=mock_solver_fn,
            corrector=mock_corrector_fn,
            rt_model=mock_rt_model,
            settings={"max_workers": 2, "retries": 0, "stage_timeout_s": 0.01},
        )


def test_run_pipeline_thread_propagates_provider_error(
    mock_preprocessor,
    mock_surface_prior_provider,
    mock_grid_assembler,
    mock_solver_fn,
    mock_corrector_fn,
    mock_rt_model,
) -> None:
    def bad_atmo(bounds: Any, crs: str, obs_time: Any, res: float) -> Any:
        _ = (bounds, crs, obs_time, res)
        raise RuntimeError("atmo boom")

    with pytest.raises(RuntimeError, match="atmo boom"):
        pipeline._run_pipeline_thread(
            Path("/fake"),
            None,
            None,
            preprocessor=mock_preprocessor,
            atmo_provider=bad_atmo,
            surface_prior_provider=mock_surface_prior_provider,
            grid_assembler=mock_grid_assembler,
            solver=mock_solver_fn,
            corrector=mock_corrector_fn,
            rt_model=mock_rt_model,
            settings={"max_workers": 2, "retries": 0, "stage_timeout_s": None},
        )


def test_run_pipeline_dask_success_path(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    mock_preprocessor,
    mock_atmo_provider,
    mock_surface_prior_provider,
    mock_grid_assembler,
    mock_solver_fn,
    mock_corrector_fn,
    mock_rt_model,
) -> None:
    fake = _install_fake_dask(monkeypatch, mode="success")
    report_path = tmp_path / "nested" / "perf" / "report.html"

    result = pipeline._run_pipeline_dask(
        Path("/fake"),
        None,
        None,
        preprocessor=mock_preprocessor,
        atmo_provider=mock_atmo_provider,
        surface_prior_provider=mock_surface_prior_provider,
        grid_assembler=mock_grid_assembler,
        solver=mock_solver_fn,
        corrector=mock_corrector_fn,
        rt_model=mock_rt_model,
        settings={
            "max_workers": 3,
            "retries": 0,
            "stage_timeout_s": 5.0,
            "dashboard": True,
            "dashboard_address": ":8787",
            "performance_report_path": report_path,
            "show_progress": True,
        },
    )

    assert isinstance(result, CorrectionResult)
    assert fake.cluster_cls.last_kwargs is not None
    assert fake.cluster_cls.last_kwargs["n_workers"] == 3
    assert fake.cluster_cls.last_kwargs["dashboard_address"] == ":8787"
    assert report_path.parent.exists()
    assert fake.report_cls.last_filename == str(report_path)
    assert fake.report_cls.entered is True
    assert fake.report_cls.exited is True


def test_run_pipeline_dask_timeout_cancels_futures(
    monkeypatch: pytest.MonkeyPatch,
    mock_preprocessor,
    mock_atmo_provider,
    mock_surface_prior_provider,
    mock_grid_assembler,
    mock_solver_fn,
    mock_corrector_fn,
    mock_rt_model,
) -> None:
    fake = _install_fake_dask(monkeypatch, mode="timeout")

    with pytest.raises(TimeoutError, match="dask backend"):
        pipeline._run_pipeline_dask(
            Path("/fake"),
            None,
            None,
            preprocessor=mock_preprocessor,
            atmo_provider=mock_atmo_provider,
            surface_prior_provider=mock_surface_prior_provider,
            grid_assembler=mock_grid_assembler,
            solver=mock_solver_fn,
            corrector=mock_corrector_fn,
            rt_model=mock_rt_model,
            settings={
                "max_workers": 2,
                "retries": 0,
                "stage_timeout_s": 1.0,
                "dashboard": False,
                "dashboard_address": None,
                "performance_report_path": None,
                "show_progress": False,
            },
        )

    client = fake.client_cls.last_instance
    assert client is not None
    assert len(client.futures) >= 2
    assert client.futures[0].canceled is True
    assert client.futures[1].canceled is True


def test_run_pipeline_dask_provider_error_cancels_futures(
    monkeypatch: pytest.MonkeyPatch,
    mock_preprocessor,
    mock_atmo_provider,
    mock_surface_prior_provider,
    mock_grid_assembler,
    mock_solver_fn,
    mock_corrector_fn,
    mock_rt_model,
) -> None:
    fake = _install_fake_dask(monkeypatch, mode="error")

    with pytest.raises(RuntimeError, match="future failure"):
        pipeline._run_pipeline_dask(
            Path("/fake"),
            None,
            None,
            preprocessor=mock_preprocessor,
            atmo_provider=mock_atmo_provider,
            surface_prior_provider=mock_surface_prior_provider,
            grid_assembler=mock_grid_assembler,
            solver=mock_solver_fn,
            corrector=mock_corrector_fn,
            rt_model=mock_rt_model,
            settings={
                "max_workers": 2,
                "retries": 0,
                "stage_timeout_s": 1.0,
                "dashboard": False,
                "dashboard_address": None,
                "performance_report_path": None,
                "show_progress": False,
            },
        )

    client = fake.client_cls.last_instance
    assert client is not None
    assert len(client.futures) >= 2
    assert client.futures[0].canceled is True
    assert client.futures[1].canceled is True


def test_run_pipeline_thread_writes_execution_profile(
    tmp_path: Path,
    mock_preprocessor,
    mock_atmo_provider,
    mock_surface_prior_provider,
    mock_grid_assembler,
    mock_solver_fn,
    mock_corrector_fn,
    mock_rt_model,
) -> None:
    report_path = tmp_path / "reports" / "execution_profile.json"

    result = pipeline.run_pipeline(
        Path("/fake"),
        None,
        None,
        preprocessor=mock_preprocessor,
        atmo_provider=mock_atmo_provider,
        surface_prior_provider=mock_surface_prior_provider,
        grid_assembler=mock_grid_assembler,
        solver=mock_solver_fn,
        corrector=mock_corrector_fn,
        rt_model=mock_rt_model,
        execution={
            "backend": "thread",
            "performance_report_path": report_path,
            "profiling_sample_interval_s": 0.01,
            "progress_heartbeat_s": 0.01,
        },
    )

    assert isinstance(result, CorrectionResult)
    assert report_path.exists()
    summary = json.loads(report_path.read_text(encoding="utf-8"))
    assert summary["run"]["status"] == "success"
    assert "M1.preprocessing" in summary["stages"]
    assert "M2.atmospheric_prior" in summary["stages"]
    assert "M6.correction" in summary["stages"]
    assert summary["counters"]["m2_done"] == 1
    assert summary["counters"]["m3_done"] == 1
    assert summary["resources"]["peak_threads"] >= 1
    assert len(summary["resources"]["resource_samples"]) >= 1
    assert Path(summary["run"]["events_path"]).exists()


def test_run_pipeline_thread_failure_writes_execution_profile(
    tmp_path: Path,
    mock_preprocessor,
    mock_surface_prior_provider,
    mock_grid_assembler,
    mock_solver_fn,
    mock_corrector_fn,
    mock_rt_model,
) -> None:
    def bad_atmo(bounds: Any, crs: str, obs_time: Any, res: float) -> Any:
        _ = (bounds, crs, obs_time, res)
        raise RuntimeError("atmo boom")

    report_path = tmp_path / "reports" / "execution_profile.json"

    with pytest.raises(RuntimeError, match="atmo boom"):
        pipeline.run_pipeline(
            Path("/fake"),
            None,
            None,
            preprocessor=mock_preprocessor,
            atmo_provider=bad_atmo,
            surface_prior_provider=mock_surface_prior_provider,
            grid_assembler=mock_grid_assembler,
            solver=mock_solver_fn,
            corrector=mock_corrector_fn,
            rt_model=mock_rt_model,
            execution={
                "backend": "thread",
                "performance_report_path": report_path,
                "profiling_sample_interval_s": 0.01,
            },
        )

    summary = json.loads(report_path.read_text(encoding="utf-8"))
    assert summary["run"]["status"] == "error"
    assert summary["stages"]["M2.atmospheric_prior"]["failed"] == 3
    assert any(error["error_message"] == "atmo boom" for error in summary["errors"])


def test_run_pipeline_thread_preloads_lut_on_atmo_grid(
    mock_preprocessor,
    mock_surface_prior_provider,
    mock_grid_assembler,
    mock_solver_fn,
    mock_corrector_fn,
) -> None:
    coarse_coords = {"y": np.arange(4), "x": np.arange(4)}

    def _coarse_atmo(bounds: Any, crs: str, obs_time: Any, res: float) -> AtmosphericState:
        _ = (bounds, crs, obs_time, res)
        aot = xr.DataArray(
            np.full((4, 4), 0.2, dtype=np.float32), dims=("y", "x"), coords=coarse_coords
        )
        tcwv = xr.DataArray(
            np.full((4, 4), 2.0, dtype=np.float32), dims=("y", "x"), coords=coarse_coords
        )
        tco3 = xr.DataArray(
            np.full((4, 4), 0.3, dtype=np.float32), dims=("y", "x"), coords=coarse_coords
        )
        unc = xr.DataArray(
            np.full((4, 4), 0.01, dtype=np.float32), dims=("y", "x"), coords=coarse_coords
        )
        elevation = xr.DataArray(
            np.zeros((4, 4), dtype=np.float32), dims=("y", "x"), coords=coarse_coords
        )
        return AtmosphericState(
            aot=aot,
            tcwv=tcwv,
            tco3=tco3,
            aot_unc=unc,
            tcwv_unc=unc,
            tco3_unc=unc,
            elevation=elevation,
        )

    captured: dict[str, Any] = {}

    class _PreloadRTModel:
        def preload_scene_subset(
            self, geometry: Any, atmo_state: AtmosphericState, bands: list[Any]
        ) -> None:
            captured["geometry_shape"] = tuple(geometry.sza.shape)
            captured["atmo_shape"] = tuple(atmo_state.aot.shape)
            captured["bands"] = [band.name for band in bands]

    result = pipeline.run_pipeline(
        Path("/fake"),
        None,
        None,
        preprocessor=mock_preprocessor,
        atmo_provider=_coarse_atmo,
        surface_prior_provider=mock_surface_prior_provider,
        grid_assembler=mock_grid_assembler,
        solver=mock_solver_fn,
        corrector=mock_corrector_fn,
        rt_model=_PreloadRTModel(),
        execution={"backend": "thread", "retries": 0},
    )

    assert isinstance(result, CorrectionResult)
    assert captured["geometry_shape"] == (4, 4)
    assert captured["atmo_shape"] == (4, 4)
    assert captured["bands"]


def test_run_pipeline_thread_preloads_lut_on_atmo_grid_with_nonstandard_dims(
    mock_preprocessor,
    mock_surface_prior_provider,
    mock_grid_assembler,
    mock_solver_fn,
    mock_corrector_fn,
) -> None:
    coarse_coords = {
        "latitude": np.array([51.5, 51.0, 50.5, 50.0], dtype=np.float32),
        "longitude": np.array([-1.5, -1.0, -0.5, 0.0], dtype=np.float32),
    }

    def _coarse_atmo(bounds: Any, crs: str, obs_time: Any, res: float) -> AtmosphericState:
        _ = (bounds, crs, obs_time, res)
        aot = xr.DataArray(
            np.full((4, 4), 0.2, dtype=np.float32),
            dims=("latitude", "longitude"),
            coords=coarse_coords,
        )
        tcwv = xr.DataArray(
            np.full((4, 4), 2.0, dtype=np.float32),
            dims=("latitude", "longitude"),
            coords=coarse_coords,
        )
        tco3 = xr.DataArray(
            np.full((4, 4), 0.3, dtype=np.float32),
            dims=("latitude", "longitude"),
            coords=coarse_coords,
        )
        unc = xr.DataArray(
            np.full((4, 4), 0.01, dtype=np.float32),
            dims=("latitude", "longitude"),
            coords=coarse_coords,
        )
        elevation = xr.DataArray(
            np.zeros((4, 4), dtype=np.float32),
            dims=("latitude", "longitude"),
            coords=coarse_coords,
        )
        return AtmosphericState(
            aot=aot,
            tcwv=tcwv,
            tco3=tco3,
            aot_unc=unc,
            tcwv_unc=unc,
            tco3_unc=unc,
            elevation=elevation,
        )

    captured: dict[str, Any] = {}

    class _PreloadRTModel:
        def preload_scene_subset(
            self, geometry: Any, atmo_state: AtmosphericState, bands: list[Any]
        ) -> None:
            captured["geometry_shape"] = tuple(geometry.sza.shape)
            captured["geometry_dims"] = tuple(geometry.sza.dims)
            captured["atmo_dims"] = tuple(atmo_state.aot.dims)
            captured["bands"] = [band.name for band in bands]

    result = pipeline.run_pipeline(
        Path("/fake"),
        None,
        None,
        preprocessor=mock_preprocessor,
        atmo_provider=_coarse_atmo,
        surface_prior_provider=mock_surface_prior_provider,
        grid_assembler=mock_grid_assembler,
        solver=mock_solver_fn,
        corrector=mock_corrector_fn,
        rt_model=_PreloadRTModel(),
        execution={"backend": "thread", "retries": 0},
    )

    assert isinstance(result, CorrectionResult)
    assert captured["geometry_shape"] == (4, 4)
    assert captured["geometry_dims"] == ("latitude", "longitude")
    assert captured["atmo_dims"] == ("latitude", "longitude")
    assert captured["bands"]


# ---------------------------------------------------------------------------
# skip_correction tests
# ---------------------------------------------------------------------------


def test_run_tail_skip_correction_returns_auxiliary_only(
    mock_observation_bundle,
    mock_atmospheric_state,
    mock_surface_prior,
    mock_solver_input_bundle,
    mock_solved_atmosphere,
    mock_rt_model,
) -> None:
    """When skip_correction=True, M5 still runs but M6 is skipped and BOA is empty."""
    solver_called = False
    corrector_called = False

    def _solver(inputs, config):
        nonlocal solver_called
        solver_called = True
        return mock_solved_atmosphere

    def _corrector(obs, solved, rt_model):
        nonlocal corrector_called
        corrector_called = True

    cfg = SimpleNamespace(
        solver=SimpleNamespace(aerosol_resolution=1000.0),
        output=SimpleNamespace(defaults=SimpleNamespace(skip_correction=True)),
    )

    result = pipeline._run_tail(
        mock_observation_bundle,
        mock_atmospheric_state,
        mock_surface_prior,
        cfg,
        grid_assembler=lambda *_args, **_kwargs: mock_solver_input_bundle,
        solver=_solver,
        corrector=_corrector,
        rt_model=mock_rt_model,
    )

    assert solver_called, "Solver should still run when skip_correction=True"
    assert not corrector_called, "Corrector should not be called when skip_correction=True"
    assert isinstance(result, CorrectionResult)
    assert len(result.boa.data_vars) == 0
    assert result.boa_unc is None
    assert result.aot is not None
    assert result.tcwv is not None
    assert result.cloud_mask is not None
    assert result.surface_prior is not None
    assert result.metadata.get("skip_correction") is True


def test_run_tail_skip_correction_preserves_monthly_composites(
    mock_observation_bundle,
    mock_atmospheric_state,
    mock_surface_prior,
    mock_solver_input_bundle,
    mock_solved_atmosphere,
    mock_rt_model,
) -> None:
    """Monthly composite outputs are still attached in skip_correction mode."""
    composite = SimpleNamespace(
        year=2023,
        month=7,
        reflectance=xr.DataArray(
            np.stack([np.full(mock_surface_prior.boa.shape, 0.14, dtype=np.float32)]),
            dims=["band", "y", "x"],
            coords={"band": ["B02"]},
        ),
        quality=xr.DataArray(
            np.full(mock_surface_prior.boa.shape, 0.5, dtype=np.float32), dims=["y", "x"]
        ),
        sample_index=xr.DataArray(
            np.full(mock_surface_prior.boa.shape, 2, dtype=np.int16), dims=["y", "x"]
        ),
    )
    surface_with_monthly = dataclasses.replace(mock_surface_prior, monthly_composites=(composite,))

    cfg = SimpleNamespace(
        solver=SimpleNamespace(aerosol_resolution=1000.0),
        output=SimpleNamespace(defaults=SimpleNamespace(skip_correction=True)),
    )

    result = pipeline._run_tail(
        mock_observation_bundle,
        mock_atmospheric_state,
        surface_with_monthly,
        cfg,
        grid_assembler=lambda *_args, **_kwargs: mock_solver_input_bundle,
        solver=lambda *_a, **_k: mock_solved_atmosphere,
        corrector=lambda *_a, **_k: None,
        rt_model=mock_rt_model,
    )

    assert result.monthly_composites is not None
    assert "2023_07" in result.monthly_composites
