"""Execution progress and profiling helpers."""

from __future__ import annotations

import json
import logging
import threading
import time
import uuid
from contextlib import contextmanager
from contextvars import ContextVar
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any

import psutil

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from collections.abc import Iterator

_CURRENT_OBSERVER: ContextVar[ExecutionObserver | None] = ContextVar(
    "siac_execution_observer",
    default=None,
)
_REGISTRY_LOCK = threading.Lock()
_OBSERVER_REGISTRY: dict[str, ExecutionObserver] = {}


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _json_safe(value: Any) -> Any:
    if value is None or isinstance(value, str | int | float | bool):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, datetime):
        return value.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_json_safe(item) for item in value]

    item = getattr(value, "item", None)
    if callable(item):
        try:
            return _json_safe(item())
        except Exception:
            pass
    return repr(value)


def _atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    tmp_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    tmp_path.replace(path)


def _events_path_for(summary_path: Path | None) -> Path | None:
    if summary_path is None:
        return None
    return summary_path.with_name(f"{summary_path.stem}_events.jsonl")


def derive_execution_report_path(output_path: Path | str | None) -> Path | None:
    if output_path is None:
        return None
    return Path(output_path) / "reports" / "execution_profile.json"


def derive_summary_report_path(
    configured_path: Path | str | None,
    *,
    backend: str,
) -> Path | None:
    if configured_path is None:
        return None

    path = Path(configured_path)
    if backend == "dask":
        return path.with_name(f"{path.stem}_summary.json")
    if path.suffix.lower() == ".json":
        return path
    if path.suffix:
        return path.with_suffix(".json")
    return path.with_name(f"{path.name}.json")


def current_execution_observer() -> ExecutionObserver | None:
    return _CURRENT_OBSERVER.get()


def register_execution_observer(observer: ExecutionObserver) -> str:
    with _REGISTRY_LOCK:
        _OBSERVER_REGISTRY[observer.run_id] = observer
    return observer.run_id


def unregister_execution_observer(observer_id: str) -> None:
    with _REGISTRY_LOCK:
        _OBSERVER_REGISTRY.pop(observer_id, None)


def resolve_execution_observer(observer_id: str | None) -> ExecutionObserver | None:
    if observer_id is None:
        return None
    with _REGISTRY_LOCK:
        return _OBSERVER_REGISTRY.get(observer_id)


@contextmanager
def bind_execution_observer(
    observer: ExecutionObserver | None,
) -> Iterator[ExecutionObserver | None]:
    token = _CURRENT_OBSERVER.set(observer)
    try:
        yield observer
    finally:
        _CURRENT_OBSERVER.reset(token)


@contextmanager
def observe_stage(
    stage: str,
    *,
    details: dict[str, Any] | None = None,
) -> Iterator[None]:
    observer = current_execution_observer()
    if observer is None:
        yield
        return

    with observer.stage(stage, details=details):
        yield


class ExecutionObserver:
    """Low-overhead execution tracing with optional resource sampling."""

    def __init__(
        self,
        *,
        backend: str,
        input_path: Path | str,
        summary_path: Path | None,
        show_progress: bool,
        sample_interval_s: float,
        heartbeat_interval_s: float,
        configured_report_path: Path | None = None,
    ) -> None:
        self.run_id = (
            f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"
        )
        self.backend = str(backend)
        self.input_path = str(input_path)
        self.summary_path = summary_path
        self.events_path = _events_path_for(summary_path)
        self.configured_report_path = (
            str(configured_report_path) if configured_report_path else None
        )
        self.show_progress = bool(show_progress)
        self.sample_interval_s = float(sample_interval_s)
        self.heartbeat_interval_s = float(heartbeat_interval_s)
        self.started_at = _utc_now()
        self.finished_at: str | None = None
        self.status = "running"

        self._lock = threading.RLock()
        self._started_monotonic_s = time.monotonic()
        self._finished_monotonic_s: float | None = None
        self._last_heartbeat_s = self._started_monotonic_s
        self._next_stage_token = 1
        self._active_stages: dict[int, dict[str, Any]] = {}
        self._stages: dict[str, dict[str, Any]] = {}
        self._counters: dict[str, int | float] = {}
        self._errors: list[dict[str, Any]] = []
        self._warnings: list[dict[str, Any]] = []
        self._resource_summary: dict[str, Any] = {
            "peak_rss_bytes": 0,
            "peak_threads": 0,
            "cpu_user_s": 0.0,
            "cpu_system_s": 0.0,
            "read_bytes": None,
            "write_bytes": None,
            "sample_interval_s": self.sample_interval_s,
        }
        self._resource_samples: list[dict[str, Any]] = []
        self._sample_resources = self.show_progress or self.summary_path is not None
        self._sampler_stop = threading.Event()
        self._sampler_thread: threading.Thread | None = None
        self._observer_id: str | None = None
        self._process = psutil.Process()

    def start(self) -> None:
        with self._lock:
            if self._observer_id is not None:
                return
            self._observer_id = register_execution_observer(self)
            self.emit(
                "run_start",
                backend=self.backend,
                input_path=self.input_path,
                configured_report_path=self.configured_report_path,
            )
            if self._sample_resources and self._sampler_thread is None:
                self._sampler_thread = threading.Thread(
                    target=self._sample_loop,
                    name="siac-execution-observer",
                    daemon=True,
                )
                self._sampler_thread.start()

    def finish(self, status: str, *, error: BaseException | None = None) -> None:
        with self._lock:
            if self.finished_at is not None:
                return
            self.status = status
            self.finished_at = _utc_now()
            self._finished_monotonic_s = time.monotonic()

        if self._sample_resources:
            self._capture_sample()
            self._sampler_stop.set()
            if self._sampler_thread is not None:
                self._sampler_thread.join(timeout=max(1.0, self.sample_interval_s + 1.0))

        if error is not None:
            self.record_error(
                stage="pipeline",
                error_type=type(error).__name__,
                error_message=str(error),
            )
        self.emit(
            "run_end",
            status=self.status,
            wall_time_s=self.wall_time_s,
            error_type=type(error).__name__ if error is not None else None,
            error_message=str(error) if error is not None else None,
        )
        if self._observer_id is not None:
            unregister_execution_observer(self._observer_id)
            self._observer_id = None

    @property
    def wall_time_s(self) -> float:
        finished = self._finished_monotonic_s
        if finished is None:
            finished = time.monotonic()
        return finished - self._started_monotonic_s

    @contextmanager
    def stage(
        self,
        stage: str,
        *,
        details: dict[str, Any] | None = None,
    ) -> Iterator[None]:
        token = self._begin_stage(stage, details=details)
        try:
            yield
        except Exception as exc:
            self._end_stage(
                token,
                status="error",
                error_type=type(exc).__name__,
                error_message=str(exc),
            )
            raise
        else:
            self._end_stage(token, status="ok")

    def emit(self, event: str, stage: str | None = None, **fields: Any) -> None:
        payload = {
            "ts": _utc_now(),
            "run_id": self.run_id,
            "event": event,
            "stage": stage,
            **{key: _json_safe(value) for key, value in fields.items() if value is not None},
        }
        with self._lock:
            self._append_event_locked(payload)

    def increment_counter(
        self,
        counter: str,
        amount: int | float = 1,
        *,
        stage: str | None = None,
    ) -> None:
        with self._lock:
            current = self._counters.get(counter, 0)
            if isinstance(current, (int, float)):
                self._counters[counter] = current + amount
            else:
                self._counters[counter] = amount
            self._append_event_locked(
                {
                    "ts": _utc_now(),
                    "run_id": self.run_id,
                    "event": "counter",
                    "stage": stage,
                    "counter": counter,
                    "value": self._counters[counter],
                }
            )

    def record_retry(
        self,
        *,
        stage: str,
        attempt: int,
        max_attempts: int,
        duration_s: float,
        error_type: str,
        error_message: str,
    ) -> None:
        self.emit(
            "warning",
            stage=stage,
            status="retrying",
            attempt=attempt,
            max_attempts=max_attempts,
            duration_s=duration_s,
            error_type=error_type,
            error_message=error_message,
        )

    def record_timeout(self, *, stage: str, timeout_s: float, backend: str) -> None:
        self.increment_counter("timeouts_total", stage=stage)
        self.emit(
            "warning",
            stage=stage,
            status="timeout",
            timeout_s=timeout_s,
            backend=backend,
        )

    def record_error(self, *, stage: str, error_type: str, error_message: str) -> None:
        with self._lock:
            self._errors.append(
                {
                    "stage": stage,
                    "error_type": error_type,
                    "error_message": error_message,
                    "ts": _utc_now(),
                }
            )
            self._flush_summary_locked()

    def log_progress(self, message: str, *, stage: str | None = None) -> None:
        self.emit("progress", stage=stage, message=message)

    def _begin_stage(self, stage: str, *, details: dict[str, Any] | None) -> int:
        with self._lock:
            token = self._next_stage_token
            self._next_stage_token += 1
            started = time.monotonic()
            self._active_stages[token] = {
                "name": stage,
                "started_monotonic_s": started,
            }
            stats = self._stages.setdefault(
                stage,
                {
                    "attempts": 0,
                    "completed": 0,
                    "failed": 0,
                    "total_duration_s": 0.0,
                    "max_duration_s": 0.0,
                    "last_duration_s": None,
                    "status": "running",
                },
            )
            stats["attempts"] += 1
            stats["status"] = "running"
            self._append_event_locked(
                {
                    "ts": _utc_now(),
                    "run_id": self.run_id,
                    "event": "stage_start",
                    "stage": stage,
                    "details": _json_safe(details) if details else None,
                }
            )
            return token

    def _end_stage(
        self,
        token: int,
        *,
        status: str,
        error_type: str | None = None,
        error_message: str | None = None,
    ) -> None:
        with self._lock:
            active = self._active_stages.pop(token, None)
            if active is None:
                return
            stage = str(active["name"])
            duration = time.monotonic() - float(active["started_monotonic_s"])
            stats = self._stages[stage]
            stats["last_duration_s"] = duration
            stats["total_duration_s"] += duration
            stats["max_duration_s"] = max(float(stats["max_duration_s"]), duration)
            stats["status"] = status
            if status == "ok":
                stats["completed"] += 1
            else:
                stats["failed"] += 1
                if error_type is not None and error_message is not None:
                    self._errors.append(
                        {
                            "stage": stage,
                            "error_type": error_type,
                            "error_message": error_message,
                            "ts": _utc_now(),
                        }
                    )
            self._append_event_locked(
                {
                    "ts": _utc_now(),
                    "run_id": self.run_id,
                    "event": "stage_end",
                    "stage": stage,
                    "status": status,
                    "duration_s": duration,
                    "error_type": error_type,
                    "error_message": error_message,
                }
            )

    def _append_event_locked(self, payload: dict[str, Any]) -> None:
        if self.events_path is not None:
            self.events_path.parent.mkdir(parents=True, exist_ok=True)
            with self.events_path.open("a", encoding="utf-8") as handle:
                handle.write(json.dumps(_json_safe(payload), sort_keys=True))
                handle.write("\n")

        if payload["event"] == "warning":
            self._warnings.append(payload)
        self._flush_summary_locked()

    def _flush_summary_locked(self) -> None:
        if self.summary_path is None:
            return

        payload = {
            "run": {
                "run_id": self.run_id,
                "status": self.status,
                "backend": self.backend,
                "input_path": self.input_path,
                "started_at": self.started_at,
                "finished_at": self.finished_at,
                "wall_time_s": self.wall_time_s,
                "configured_report_path": self.configured_report_path,
                "events_path": str(self.events_path) if self.events_path is not None else None,
            },
            "stages": _json_safe(self._stages),
            "active_stages": _json_safe(self._active_stage_snapshot_locked()),
            "counters": _json_safe(self._counters),
            "resources": _json_safe(
                {
                    **self._resource_summary,
                    "resource_samples": self._resource_samples,
                }
            ),
            "warnings": _json_safe(self._warnings),
            "errors": _json_safe(self._errors),
        }
        _atomic_write_json(self.summary_path, payload)

    def _sample_loop(self) -> None:
        while not self._sampler_stop.wait(self.sample_interval_s):
            self._capture_sample()

    def _capture_sample(self) -> None:
        try:
            mem = self._process.memory_info()
            cpu = self._process.cpu_times()
            num_threads = self._process.num_threads()
        except (psutil.Error, OSError):
            return
        io_counters = getattr(self._process, "io_counters", None)
        io = None
        if callable(io_counters):
            try:
                io = io_counters()
            except (psutil.Error, OSError):
                io = None

        sample = {
            "ts": _utc_now(),
            "run_elapsed_s": self.wall_time_s,
            "rss_bytes": mem.rss,
            "cpu_user_s": cpu.user,
            "cpu_system_s": cpu.system,
            "num_threads": num_threads,
            "read_bytes": getattr(io, "read_bytes", None) if io is not None else None,
            "write_bytes": getattr(io, "write_bytes", None) if io is not None else None,
            "active_stages": _json_safe(self.active_stage_snapshot()),
        }

        with self._lock:
            self._resource_summary["peak_rss_bytes"] = max(
                int(self._resource_summary["peak_rss_bytes"]),
                int(mem.rss),
            )
            self._resource_summary["peak_threads"] = max(
                int(self._resource_summary["peak_threads"]),
                int(num_threads),
            )
            self._resource_summary["cpu_user_s"] = float(cpu.user)
            self._resource_summary["cpu_system_s"] = float(cpu.system)
            self._resource_summary["read_bytes"] = (
                getattr(io, "read_bytes", None) if io is not None else None
            )
            self._resource_summary["write_bytes"] = (
                getattr(io, "write_bytes", None) if io is not None else None
            )
            if self.summary_path is not None:
                self._resource_samples.append(sample)
                self._append_event_locked(
                    {
                        "ts": sample["ts"],
                        "run_id": self.run_id,
                        "event": "resource_sample",
                        "stage": None,
                        **sample,
                    }
                )
            else:
                self._flush_summary_locked()

        if not self.show_progress:
            return

        now = time.monotonic()
        if now - self._last_heartbeat_s < self.heartbeat_interval_s:
            return

        active = self.active_stage_snapshot()
        if not active:
            return
        active_text = ", ".join(f"{item['stage']} ({item['elapsed_s']:.1f}s)" for item in active)
        logger.info(
            "Progress heartbeat: active=%s; run=%.1fs; rss=%.2f GiB; threads=%d",
            active_text,
            sample["run_elapsed_s"],
            float(mem.rss) / (1024.0**3),
            num_threads,
        )
        self._last_heartbeat_s = now

    def active_stage_snapshot(self) -> list[dict[str, Any]]:
        with self._lock:
            return self._active_stage_snapshot_locked()

    def _active_stage_snapshot_locked(self) -> list[dict[str, Any]]:
        now = time.monotonic()
        return [
            {
                "stage": str(item["name"]),
                "elapsed_s": now - float(item["started_monotonic_s"]),
            }
            for item in self._active_stages.values()
        ]
