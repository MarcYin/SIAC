from __future__ import annotations

import pytest

from siac.config.schema import RTAlgorithmConfig
from siac.rt_setup import (
    DEFAULT_LUT_RT_SETUP,
    DEFAULT_SIXS_RT_SETUP,
    resolve_effective_rt_setup,
    validate_lut_requested_setup,
)


def test_native_sixs_effective_setup_uses_generic_defaults() -> None:
    effective = resolve_effective_rt_setup(RTAlgorithmConfig(backend="sixs"), "sixs")

    assert effective.atmosphere is not None
    assert effective.atmosphere.profile == DEFAULT_SIXS_RT_SETUP.atmosphere.profile
    assert effective.atmosphere.columns_mode == "input_columns"
    assert effective.aerosol is not None
    assert effective.aerosol.profile == "continental"
    assert effective.surface is not None
    assert effective.surface.mode == "homogeneous_lambertian"
    assert effective.surface.target is not None
    assert effective.surface.target.constant == pytest.approx(0.0)
    assert effective.atmospheric_correction is not None
    assert effective.atmospheric_correction.mode == "lambertian_reflectance"
    assert effective.reference_reflectance == pytest.approx(0.1)


def test_requested_rt_setup_overlays_native_defaults() -> None:
    rt_config = RTAlgorithmConfig(
        backend="sixs",
        setup={
            "atmosphere": {
                "profile": "user_water_ozone",
            },
            "surface": {
                "mode": "homogeneous_brdf",
                "brdf": {
                    "model": "rahman",
                    "parameters": {
                        "intensity": 0.12,
                        "asymmetry_factor": 0.03,
                        "structural_parameter": 0.45,
                    },
                },
            },
            "atmospheric_correction": {
                "mode": "brdf_reflectance",
                "value": 0.2,
            },
            "reference_reflectance": 0.2,
        },
    )

    effective = resolve_effective_rt_setup(rt_config, "sixs")

    assert effective.atmosphere is not None
    assert effective.atmosphere.profile == "user_water_ozone"
    assert effective.atmosphere.columns_mode == "input_columns"
    assert effective.aerosol is not None
    assert effective.aerosol.profile == "continental"
    assert effective.surface is not None
    assert effective.surface.mode == "homogeneous_brdf"
    assert effective.surface.brdf is not None
    assert effective.surface.brdf.model == "rahman"
    assert effective.atmospheric_correction is not None
    assert effective.atmospheric_correction.mode == "brdf_reflectance"
    assert effective.reference_reflectance == pytest.approx(0.2)


def test_native_sixs_rejects_lut_specific_aerosol_family() -> None:
    rt_config = RTAlgorithmConfig(
        backend="sixs",
        setup={"aerosol": {"profile": "continental_average"}},
    )

    with pytest.raises(ValueError, match="not a native 6S profile"):
        resolve_effective_rt_setup(rt_config, "sixs")


def test_lut_setup_validation_rejects_incompatible_requested_structure() -> None:
    with pytest.raises(ValueError, match="fixed packaged remote libRadtran LUT preset"):
        validate_lut_requested_setup(
            resolve_effective_rt_setup(
                RTAlgorithmConfig(
                    backend="emulator",
                    setup={
                        "surface": {
                            "mode": "homogeneous_brdf",
                            "brdf": {"model": "rahman"},
                        }
                    },
                ),
                "emulator",
            )
        )


@pytest.mark.parametrize(
    ("setup", "expected"),
    [
        (
            {"surface": {"target": {"kind": "constant", "constant": 0.1}}},
            "surface.target cannot be configured",
        ),
        (
            {"reference_reflectance": 0.1},
            "reference_reflectance cannot be configured",
        ),
        (
            {"atmospheric_correction": {"mode": "lambertian_reflectance", "value": 0.1}},
            "atmospheric_correction cannot be configured",
        ),
        (
            {"atmosphere": {"columns_mode": "input_columns"}},
            "atmosphere.columns_mode cannot be configured",
        ),
    ],
)
def test_lut_setup_validation_rejects_semantic_overrides(
    setup: dict[str, object],
    expected: str,
) -> None:
    with pytest.raises(ValueError, match=expected):
        resolve_effective_rt_setup(RTAlgorithmConfig(backend="lut", setup=setup), "lut")


def test_lut_effective_setup_exposes_packaged_preset() -> None:
    rt_config = RTAlgorithmConfig(
        backend="lut",
        setup={
            "atmosphere": {"profile": "us_standard_62"},
            "aerosol": {"profile": "continental_average"},
            "surface": {"mode": "homogeneous_lambertian"},
        },
    )

    effective = resolve_effective_rt_setup(rt_config, "lut")

    assert effective.atmosphere is not None
    assert effective.atmosphere.profile == DEFAULT_LUT_RT_SETUP.atmosphere.profile
    assert effective.aerosol is not None
    assert effective.aerosol.profile == "continental_average"
    assert effective.surface is not None
    assert effective.surface.mode == "homogeneous_lambertian"
    assert effective.surface.target is None
    assert effective.atmospheric_correction is None
    assert effective.reference_reflectance is None
