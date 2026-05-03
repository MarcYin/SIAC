"""Coverage tests for package-level lazy exports."""

from __future__ import annotations

import pytest

import siac
import siac.config as siac_config
import siac.storage as siac_storage
from siac.adapters import (
    load_band_rsrf,
    load_sensor_config_with_rsrf,
)
from siac.algorithms.rt.lut.rsrf_kernel import AlignedRSRFKernel, build_aligned_rsrf_kernel
from siac.algorithms.surface import load_reference_rsrf
from siac.domain import RelativeSpectralResponse


def test_siac_module_lazy_exports_and_unknown_attr():
    assert siac.SIAC.__name__ == "SIAC"
    assert siac.SIACConfig.__name__ == "SIACConfig"
    assert callable(siac.process_sentinel2)
    assert callable(siac.resolve_s2_input)
    assert callable(siac.siac_process_s2)
    assert callable(siac.search_sentinel2)

    with pytest.raises(AttributeError, match="has no attribute"):
        _ = siac.NOT_A_REAL_EXPORT


def test_siac_config_export_and_unknown_attr():
    assert siac_config.SIACConfig.__name__ == "SIACConfig"

    with pytest.raises(AttributeError, match="has no attribute"):
        _ = siac_config.NOT_A_REAL_EXPORT


def test_siac_storage_package_exports():
    assert callable(siac_storage.read_multiband)
    assert callable(siac_storage.write_dataset)
    assert callable(siac_storage.build_stac_item_from_result)


def test_current_spectral_exports_are_available():
    assert RelativeSpectralResponse.__name__ == "RelativeSpectralResponse"
    assert callable(load_band_rsrf)
    assert callable(load_sensor_config_with_rsrf)
    assert callable(load_reference_rsrf)
    assert AlignedRSRFKernel.__name__ == "AlignedRSRFKernel"
    assert callable(build_aligned_rsrf_kernel)
