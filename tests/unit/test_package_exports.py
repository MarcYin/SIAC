"""Coverage tests for package-level lazy exports."""

from __future__ import annotations

import pytest

import siac
import siac.core as siac_core


def test_siac_module_lazy_exports_and_unknown_attr():
    assert siac.SIAC.__name__ == "SIAC"
    assert siac.SIACConfig.__name__ == "SIACConfig"
    assert callable(siac.process_sentinel2)
    assert callable(siac.process_landsat8)
    assert callable(siac.resolve_s2_input)
    assert callable(siac.siac_process_s2)
    assert callable(siac.search_sentinel2)

    with pytest.raises(AttributeError, match="has no attribute"):
        _ = siac.NOT_A_REAL_EXPORT


def test_siac_core_lazy_config_export_and_unknown_attr():
    assert siac_core.SIACConfig.__name__ == "SIACConfig"

    with pytest.raises(AttributeError, match="has no attribute"):
        _ = siac_core.NOT_A_REAL_EXPORT
