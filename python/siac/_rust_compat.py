"""Compatibility layer for the optional ``siac._rust`` native extension."""

from __future__ import annotations

import importlib
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from types import ModuleType

_NATIVE_IMPORT_ERROR: Exception | None = None


def _load_native_rust() -> ModuleType | None:
    """Import the optional native extension if it is available.

    NOTE: ``Exception`` is broad on purpose — native extension load failures
    can surface as ``OSError``/``RuntimeError`` from the dynamic linker (e.g.
    missing libgomp), not just ``ImportError``. ``BaseException`` subclasses
    such as ``SystemExit``/``KeyboardInterrupt`` still propagate. The contract
    is exercised by ``tests/unit/test_rust_compat.py``.
    """
    global _NATIVE_IMPORT_ERROR
    try:
        return importlib.import_module("siac._rust")
    except Exception as exc:  # noqa: BLE001  - see docstring
        _NATIVE_IMPORT_ERROR = exc
        return None


_NATIVE_RUST = _load_native_rust()


def _missing_message(symbol: str) -> str:
    message = (
        "SIAC's optional Rust extension `siac._rust` is unavailable. "
        f"The symbol `{symbol}` is required for this operation. "
        "Build the extension with `maturin develop --release` or install a wheel "
        "that includes the native module."
    )
    if _NATIVE_IMPORT_ERROR is not None:
        message = f"{message} Original import failure: {_NATIVE_IMPORT_ERROR}"
    return message


def _require_native(symbol: str) -> Any:
    if _NATIVE_RUST is None:
        raise ImportError(_missing_message(symbol)) from _NATIVE_IMPORT_ERROR
    try:
        return getattr(_NATIVE_RUST, symbol)
    except AttributeError as exc:  # pragma: no cover - defensive
        raise ImportError(_missing_message(symbol)) from exc


class RossThickLiSparse:
    """Proxy for the native BRDF kernel calculator."""

    def __init__(self, hb: float = 2.0, br: float = 1.0) -> None:
        native_cls = _require_native("RossThickLiSparse")
        self._impl: Any = native_cls(hb, br)

    def compute(self, vza: Any, sza: Any, raa: Any) -> Any:
        return self._impl.compute(vza, sza, raa)

    def __getattr__(self, name: str) -> Any:
        return getattr(self._impl, name)


class TwoLayerNN:
    """Proxy for the native neural-network emulator."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        native_cls = _require_native("TwoLayerNN")
        self._impl: Any = native_cls(*args, **kwargs)

    def predict(self, *args: Any, **kwargs: Any) -> Any:
        return self._impl.predict(*args, **kwargs)

    def __getattr__(self, name: str) -> Any:
        return getattr(self._impl, name)


class PSFConvolver:
    """Proxy for the native PSF convolution helper."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        native_cls = _require_native("PSFConvolver")
        self._impl: Any = native_cls(*args, **kwargs)

    def convolve(self, *args: Any, **kwargs: Any) -> Any:
        return self._impl.convolve(*args, **kwargs)

    def __getattr__(self, name: str) -> Any:
        return getattr(self._impl, name)


def apply_laplacian(*args: Any, **kwargs: Any) -> Any:
    return _require_native("apply_laplacian")(*args, **kwargs)


def evaluate_grid_search_candidate_cost(*args: Any, **kwargs: Any) -> Any:
    return _require_native("evaluate_grid_search_candidate_cost")(*args, **kwargs)


def evaluate_grid_search_cost_cube_with_provider(*args: Any, **kwargs: Any) -> Any:
    return _require_native("evaluate_grid_search_cost_cube_with_provider")(*args, **kwargs)


def evaluate_grid_search_cost_cube_with_provider_qa(*args: Any, **kwargs: Any) -> Any:
    return _require_native("evaluate_grid_search_cost_cube_with_provider_qa")(*args, **kwargs)


def evaluate_block_grid_search_cost_cube_with_provider_qa(*args: Any, **kwargs: Any) -> Any:
    return _require_native("evaluate_block_grid_search_cost_cube_with_provider_qa")(*args, **kwargs)


def interpolate_to_fine_grid(*args: Any, **kwargs: Any) -> Any:
    return _require_native("interpolate_to_fine_grid")(*args, **kwargs)


def quadratic_refine_grid_search(*args: Any, **kwargs: Any) -> Any:
    return _require_native("quadratic_refine_grid_search")(*args, **kwargs)


def quadratic_refine_grid_search_qa(*args: Any, **kwargs: Any) -> Any:
    return _require_native("quadratic_refine_grid_search_qa")(*args, **kwargs)


def remap_to_coarse_grid(*args: Any, **kwargs: Any) -> Any:
    return _require_native("remap_to_coarse_grid")(*args, **kwargs)


def surface_driven_pool_argmin(*args: Any, **kwargs: Any) -> Any:
    return _require_native("surface_driven_pool_argmin")(*args, **kwargs)


def whittaker_smooth_cube(*args: Any, **kwargs: Any) -> Any:
    return _require_native("whittaker_smooth_cube")(*args, **kwargs)
