"""Guardrails for removed compatibility import paths."""

from __future__ import annotations

import importlib

import pytest

import siac
import siac.config as siac_config
from siac.domain import SensorBand


@pytest.mark.parametrize(
    "module_name",
    [
        "siac.config.schema",
        "siac.app.assembly",
        "siac.app._assembly_runtime",
        "siac.algorithms.rt.lut.srf_kernel",
        "siac.adapters.brdf.vnp43_earthaccess",
    ],
)
def test_removed_compatibility_modules_are_not_importable(module_name: str) -> None:
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module(module_name)


@pytest.mark.parametrize(
    ("owner", "name"),
    [
        (siac, "process_landsat8"),
        (siac_config, "CredentialConfig"),
        (siac_config, "OutputSectionConfig"),
        (siac_config, "AtmoPriorConfig"),
        (siac_config, "BRDFConfig"),
        (siac_config, "S2DataAccessConfig"),
    ],
)
def test_removed_compatibility_exports_are_absent(owner: object, name: str) -> None:
    assert not hasattr(owner, name)


def test_sensor_band_rejects_srf_keyword_aliases() -> None:
    with pytest.raises(TypeError):
        SensorBand(  # type: ignore[call-arg]
            "B02",
            490.0,
            65.0,
            10.0,
            0,
            srf_wavelengths_nm=[],
        )
