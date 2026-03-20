"""Unit tests for pipeline execution settings, retries, and backends."""

from __future__ import annotations

import sys
import time
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any

import pytest

import siac.workflows.pipeline as pipeline
from siac.runtime import CorrectionResult


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
        execution=SimpleNamespace(
            backend="THREAD",
            max_workers=2,
            retries=1,
            stage_timeout_s="15",
            dashboard=1,
            dashboard_address=":8787",
            performance_report_path=str(tmp_path / "reports" / "perf.html"),
            show_progress=1,
        )
    )

    settings = pipeline._resolve_execution_settings(
        cfg,
        execution={"retries": 3},
        max_workers=6,
    )

    assert settings["backend"] == "thread"
    assert settings["max_workers"] == 6
    assert settings["retries"] == 3
    assert settings["stage_timeout_s"] == 15.0
    assert settings["dashboard"] is True
    assert settings["show_progress"] is True
    assert isinstance(settings["performance_report_path"], Path)


@pytest.mark.parametrize(
    ("execution", "message"),
    [
        ({"backend": "unknown"}, "Unsupported execution backend"),
        ({"max_workers": 0}, "max_workers must be >= 1"),
        ({"retries": -1}, "retries must be >= 0"),
        ({"stage_timeout_s": 0}, "stage_timeout_s must be > 0"),
    ],
)
def test_resolve_execution_settings_rejects_invalid_values(
    execution: dict[str, Any],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        pipeline._resolve_execution_settings(SimpleNamespace(), execution=execution, max_workers=None)


def test_aerosol_resolution_reads_solver_config() -> None:
    cfg = SimpleNamespace(solver=SimpleNamespace(aerosol_resolution=250.0))
    assert pipeline._aerosol_resolution(cfg) == 250.0


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
