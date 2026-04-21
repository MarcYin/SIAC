"""Generic radiative-transfer setup helpers."""

from __future__ import annotations

from typing import Any

from siac.config.schema import (
    RTAerosolSetupConfig,
    RTAlgorithmConfig,
    RTAtmosphereSetupConfig,
    RTAtmosphericCorrectionSetupConfig,
    RTSetupConfig,
    RTSurfaceSetupConfig,
    SixSSpectralReflectanceConfig,
)

DEFAULT_SIXS_RT_SETUP = RTSetupConfig(
    atmosphere=RTAtmosphereSetupConfig(
        profile="us_standard_62",
        columns_mode="input_columns",
    ),
    aerosol=RTAerosolSetupConfig(profile="continental"),
    surface=RTSurfaceSetupConfig(
        mode="homogeneous_lambertian",
        target=SixSSpectralReflectanceConfig(kind="constant", constant=0.0),
    ),
    atmospheric_correction=RTAtmosphericCorrectionSetupConfig(
        mode="lambertian_reflectance",
        value=0.1,
    ),
    reference_reflectance=0.1,
)

DEFAULT_LUT_RT_SETUP = RTSetupConfig(
    atmosphere=RTAtmosphereSetupConfig(profile="us_standard_62"),
    aerosol=RTAerosolSetupConfig(profile="continental_average"),
    surface=RTSurfaceSetupConfig(mode="homogeneous_lambertian"),
)

_SIXS_UNSUPPORTED_RT_AEROSOLS = {"continental_average"}


def _coerce_rt_setup(value: Any) -> RTSetupConfig:
    if isinstance(value, RTSetupConfig):
        return value
    if value is None:
        return RTSetupConfig()
    if hasattr(value, "model_dump"):
        return RTSetupConfig.model_validate(value.model_dump(mode="python"))
    if hasattr(value, "__dict__"):
        return RTSetupConfig.model_validate(vars(value))
    return RTSetupConfig.model_validate(value)


def _merge_model_payload(base: Any, override: Any) -> dict[str, Any]:
    payload = base.model_dump(mode="python") if hasattr(base, "model_dump") else dict(base or {})
    if override is None:
        return payload
    payload.update(override.model_dump(mode="python", exclude_none=True))
    return payload


def _merge_rt_setup(base: RTSetupConfig, override: RTSetupConfig | None) -> RTSetupConfig:
    if override is None or not override.has_overrides():
        return base.model_copy(deep=True)
    payload = base.model_dump(mode="python")
    if override.atmosphere is not None:
        payload["atmosphere"] = _merge_model_payload(payload.get("atmosphere") or {}, override.atmosphere)
    if override.aerosol is not None:
        payload["aerosol"] = _merge_model_payload(payload.get("aerosol") or {}, override.aerosol)
    if override.surface is not None:
        payload["surface"] = _merge_model_payload(payload.get("surface") or {}, override.surface)
    if override.atmospheric_correction is not None:
        payload["atmospheric_correction"] = _merge_model_payload(
            payload.get("atmospheric_correction") or {},
            override.atmospheric_correction,
        )
    if override.reference_reflectance is not None:
        payload["reference_reflectance"] = override.reference_reflectance
    return RTSetupConfig.model_validate(payload)


def resolve_backend_rt_setup(backend: str, setup: Any) -> RTSetupConfig:
    """Resolve an explicit RT setup payload against backend defaults."""

    requested = _coerce_rt_setup(setup)
    if backend == "sixs":
        effective = _merge_rt_setup(DEFAULT_SIXS_RT_SETUP, requested)
        aerosol = effective.aerosol
        if aerosol is not None and aerosol.profile in _SIXS_UNSUPPORTED_RT_AEROSOLS:
            raise ValueError(
                "RT setup aerosol profile "
                f"{aerosol.profile!r} is not a native 6S profile. "
                "Use a native 6S aerosol family such as 'continental', or select the LUT backend "
                "for the packaged libRadtran continental-average preset."
            )
        return effective
    if backend == "lut":
        validate_lut_requested_setup(requested)
        return DEFAULT_LUT_RT_SETUP.model_copy(deep=True)
    return requested


def resolve_effective_rt_setup(rt_config: RTAlgorithmConfig, backend: str) -> RTSetupConfig:
    """Resolve the effective generic RT setup for a backend."""

    return resolve_backend_rt_setup(backend, getattr(rt_config, "setup", None))


def validate_lut_requested_setup(requested: RTSetupConfig) -> None:
    """Validate that a requested generic setup matches the fixed packaged LUT preset."""

    mismatches: list[str] = []
    atmosphere = requested.atmosphere
    aerosol = requested.aerosol
    surface = requested.surface

    if atmosphere is not None:
        if atmosphere.profile not in {None, "us_standard_62"}:
            mismatches.append(
                f"atmosphere.profile={atmosphere.profile!r} is incompatible with the packaged LUT "
                "preset atmosphere.profile='us_standard_62'"
            )
        if atmosphere.columns_mode is not None:
            mismatches.append(
                "atmosphere.columns_mode cannot be configured for the packaged LUT preset; "
                "the remote libRadtran LUT already uses the scene atmospheric-state columns directly"
            )
        if atmosphere.profile_latitude is not None:
            mismatches.append("atmosphere.profile_latitude is not supported by the packaged LUT preset")
        if atmosphere.radiosonde_profile is not None:
            mismatches.append("atmosphere.radiosonde_profile is not supported by the packaged LUT preset")

    if aerosol is not None:
        if aerosol.profile not in {None, "continental_average"}:
            mismatches.append(
                f"aerosol.profile={aerosol.profile!r} is incompatible with the packaged LUT "
                "preset aerosol.profile='continental_average'"
            )
        for field_name in (
            "mixture",
            "distribution",
            "sun_photometer_aerosol",
            "layer_profile",
            "model_path",
        ):
            if getattr(aerosol, field_name) is not None:
                mismatches.append(f"aerosol.{field_name} is not supported by the packaged LUT preset")

    if surface is not None:
        if surface.mode not in {None, "homogeneous_lambertian"}:
            mismatches.append(
                f"surface.mode={surface.mode!r} is incompatible with the packaged LUT "
                "preset surface.mode='homogeneous_lambertian'"
            )
        if surface.target is not None:
            mismatches.append(
                "surface.target cannot be configured for the packaged LUT preset; "
                "the remote libRadtran LUT fixes the surface semantics internally"
            )
        if surface.environment is not None:
            mismatches.append("surface.environment is not supported by the packaged LUT preset")
        if surface.brdf is not None:
            mismatches.append("surface.brdf is not supported by the packaged LUT preset")

    correction = requested.atmospheric_correction
    if correction is not None:
        mismatches.append(
            "atmospheric_correction cannot be configured for the packaged LUT preset; "
            "the remote libRadtran LUT exposes fixed SIAC coefficient semantics"
        )

    if requested.reference_reflectance is not None:
        mismatches.append(
            "reference_reflectance cannot be configured for the packaged LUT preset"
        )

    if mismatches:
        detail = "; ".join(mismatches)
        raise ValueError(
            "The requested generic RT setup is incompatible with the fixed packaged remote libRadtran LUT preset: "
            + detail
        )


__all__ = [
    "DEFAULT_LUT_RT_SETUP",
    "DEFAULT_SIXS_RT_SETUP",
    "resolve_backend_rt_setup",
    "resolve_effective_rt_setup",
    "validate_lut_requested_setup",
]
