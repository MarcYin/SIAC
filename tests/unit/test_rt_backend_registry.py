"""Tests for the RT backend registry and optional-capability accessors."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from siac.adapters import rt as rt_adapters
from siac.adapters.rt import RTBuildContext, build_rt_model, register_rt_backend
from siac.domain.protocols import (
    RTModelBackend,
    rt_optional_capability,
    rt_supports_jacobian,
)


def _config(backend: str, **rt_extra: object) -> SimpleNamespace:
    rt = SimpleNamespace(
        backend=backend,
        fallback_to_lut=False,
        lut_path=None,
        emulator_dir=None,
        lut_storage_options={},
        lut_interpolation="linear",
        **rt_extra,
    )
    return SimpleNamespace(algorithms=SimpleNamespace(rt=rt), sensor="s2a", paths=None)


def test_unknown_backend_raises_with_backend_name() -> None:
    with pytest.raises(ValueError, match="backend='nope'"):
        build_rt_model(_config("nope"))


def test_register_rt_backend_dispatches_custom_builder() -> None:
    seen: dict[str, RTBuildContext] = {}

    def _builder(ctx: RTBuildContext) -> str:
        seen["ctx"] = ctx
        return "custom-backend"

    register_rt_backend("custom-test", _builder)
    try:
        result = build_rt_model(_config("custom-test"))
    finally:
        rt_adapters._RT_BACKEND_BUILDERS.pop("custom-test", None)
    assert result == "custom-backend"
    assert seen["ctx"].sensor_id == "MSI"
    assert seen["ctx"].satellite_id == "S2A"


def test_builders_honor_module_level_test_seams(monkeypatch: pytest.MonkeyPatch) -> None:
    """The module-level backend-class globals override the lazy imports."""

    class _FakeSixS:
        def __init__(self, **kwargs: object) -> None:
            self.kwargs = kwargs

    monkeypatch.setattr(rt_adapters, "SixSBackend", _FakeSixS)
    model = build_rt_model(_config("sixs", sixs=SimpleNamespace(), setup=None))
    assert isinstance(model, _FakeSixS)
    assert "sixs_config" in model.kwargs


def test_lut_backend_without_path_is_unresolvable() -> None:
    with pytest.raises(ValueError, match="backend='lut'"):
        build_rt_model(_config("lut"))


def test_rt_optional_capability_returns_callable_or_none() -> None:
    class _Backend:
        def preload_scene_subset(self) -> str:
            return "preloaded"

        not_callable = 42

    backend = _Backend()
    hook = rt_optional_capability(backend, "preload_scene_subset")
    assert hook is not None
    assert hook() == "preloaded"
    assert rt_optional_capability(backend, "not_callable") is None
    assert rt_optional_capability(backend, "missing_hook") is None


def test_rt_supports_jacobian_false_on_missing_or_raising() -> None:
    class _NoHook:
        pass

    class _Raises:
        def supports_jacobian(self) -> bool:
            raise RuntimeError("boom")

    class _Yes:
        def supports_jacobian(self) -> bool:
            return True

    assert rt_supports_jacobian(_NoHook()) is False
    assert rt_supports_jacobian(_Raises()) is False
    assert rt_supports_jacobian(_Yes()) is True


def test_real_backends_satisfy_rt_model_protocol() -> None:
    """Engine-free constructors conform to the runtime-checkable protocol."""
    from siac.algorithms.rt.direct.libradtran import LibRadtranBackend
    from siac.algorithms.rt.lut.backend import ZarrLUTBackend
    from siac.config.algorithms import LibRadtranAlgorithmConfig

    lut = ZarrLUTBackend("dummy.zarr.zip", scene_cache_enabled=False)
    libradtran = LibRadtranBackend(libradtran_config=LibRadtranAlgorithmConfig())
    assert isinstance(lut, RTModelBackend)
    assert isinstance(libradtran, RTModelBackend)
