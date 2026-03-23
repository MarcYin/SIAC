"""Tests for the optional Rust compatibility layer."""

from __future__ import annotations

import numpy as np
import pytest

import siac._rust_compat as rust_compat
from siac.algorithms.brdf.kernels import BRDFKernels


def test_missing_rust_extension_raises_clear_function_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(rust_compat, "_NATIVE_RUST", None)
    monkeypatch.setattr(
        rust_compat,
        "_NATIVE_IMPORT_ERROR",
        ModuleNotFoundError("No module named 'siac._rust'"),
    )

    with pytest.raises(ImportError, match="siac\\._rust"):
        rust_compat.quadratic_refine_grid_search(
            np.zeros((3, 3, 1, 1), dtype=np.float32),
            np.linspace(0.0, 1.0, 3, dtype=np.float32),
            np.linspace(0.0, 1.0, 3, dtype=np.float32),
            np.ones((1, 1), dtype=bool),
        )


def test_missing_rust_extension_raises_clear_class_error(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(rust_compat, "_NATIVE_RUST", None)
    monkeypatch.setattr(
        rust_compat,
        "_NATIVE_IMPORT_ERROR",
        ModuleNotFoundError("No module named 'siac._rust'"),
    )

    with pytest.raises(ImportError, match="RossThickLiSparse"):
        BRDFKernels()


def test_load_native_rust_success_and_failure_paths(monkeypatch: pytest.MonkeyPatch) -> None:
    native_module = object()
    monkeypatch.setattr(
        rust_compat.importlib,
        "import_module",
        lambda name: native_module if name == "siac._rust" else None,
    )
    rust_compat._NATIVE_IMPORT_ERROR = None
    assert rust_compat._load_native_rust() is native_module  # noqa: SLF001
    assert rust_compat._NATIVE_IMPORT_ERROR is None

    def _boom(_name: str) -> object:
        raise RuntimeError("native boom")

    monkeypatch.setattr(rust_compat.importlib, "import_module", _boom)
    assert rust_compat._load_native_rust() is None  # noqa: SLF001
    assert isinstance(rust_compat._NATIVE_IMPORT_ERROR, RuntimeError)


def test_missing_message_and_require_native_missing_symbol(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(rust_compat, "_NATIVE_IMPORT_ERROR", RuntimeError("broken import"))
    message = rust_compat._missing_message("apply_laplacian")  # noqa: SLF001
    assert "apply_laplacian" in message
    assert "broken import" in message

    monkeypatch.setattr(rust_compat, "_NATIVE_RUST", object())
    with pytest.raises(ImportError, match="TwoLayerNN"):
        rust_compat._require_native("TwoLayerNN")  # noqa: SLF001


def test_proxy_classes_and_functions_delegate_to_native_impl(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: dict[str, object] = {}

    class _FakeRoss:
        def __init__(self, hb: float, br: float) -> None:
            calls["ross_init"] = (hb, br)
            self.scale = 3

        def compute(self, vza, sza, raa):  # noqa: ANN001
            return ("ross", vza, sza, raa)

    class _FakeNN:
        def __init__(self, *args, **kwargs) -> None:  # noqa: ANN002, ANN003
            calls["nn_init"] = (args, kwargs)
            self.kind = "nn"

        def predict(self, *args, **kwargs):  # noqa: ANN002, ANN003
            return ("predict", args, kwargs)

    class _FakePSF:
        def __init__(self, *args, **kwargs) -> None:  # noqa: ANN002, ANN003
            calls["psf_init"] = (args, kwargs)
            self.kind = "psf"

        def convolve(self, *args, **kwargs):  # noqa: ANN002, ANN003
            return ("convolve", args, kwargs)

    native = type(
        "Native",
        (),
        {
            "RossThickLiSparse": _FakeRoss,
            "TwoLayerNN": _FakeNN,
            "PSFConvolver": _FakePSF,
            "apply_laplacian": staticmethod(lambda *args, **kwargs: ("lap", args, kwargs)),
            "evaluate_grid_search_candidate_cost": staticmethod(
                lambda *args, **kwargs: ("candidate", args, kwargs)
            ),
            "evaluate_grid_search_cost_cube_with_provider": staticmethod(
                lambda *args, **kwargs: ("cube", args, kwargs)
            ),
            "interpolate_to_fine_grid": staticmethod(lambda *args, **kwargs: ("interp", args, kwargs)),
            "quadratic_refine_grid_search": staticmethod(lambda *args, **kwargs: ("quad", args, kwargs)),
            "remap_to_coarse_grid": staticmethod(lambda *args, **kwargs: ("remap", args, kwargs)),
            "whittaker_smooth_cube": staticmethod(lambda *args, **kwargs: ("smooth", args, kwargs)),
        },
    )()
    monkeypatch.setattr(rust_compat, "_NATIVE_RUST", native)
    monkeypatch.setattr(rust_compat, "_NATIVE_IMPORT_ERROR", None)

    ross = rust_compat.RossThickLiSparse(2.5, 1.5)
    nn = rust_compat.TwoLayerNN("weights", backend="cpu")
    psf = rust_compat.PSFConvolver("kernel")

    assert ross.compute("vza", "sza", "raa") == ("ross", "vza", "sza", "raa")
    assert ross.scale == 3
    assert nn.predict("x") == ("predict", ("x",), {})
    assert nn.kind == "nn"
    assert psf.convolve("x") == ("convolve", ("x",), {})
    assert psf.kind == "psf"
    assert rust_compat.apply_laplacian(1) == ("lap", (1,), {})
    assert rust_compat.evaluate_grid_search_candidate_cost(2) == ("candidate", (2,), {})
    assert rust_compat.evaluate_grid_search_cost_cube_with_provider(3) == ("cube", (3,), {})
    assert rust_compat.interpolate_to_fine_grid(4) == ("interp", (4,), {})
    assert rust_compat.quadratic_refine_grid_search(5) == ("quad", (5,), {})
    assert rust_compat.remap_to_coarse_grid(6) == ("remap", (6,), {})
    assert rust_compat.whittaker_smooth_cube(7) == ("smooth", (7,), {})
