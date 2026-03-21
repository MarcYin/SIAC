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
