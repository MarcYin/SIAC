"""Thin adapter layer over the external RSRF package."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rsrf

from siac.catalog import get_sensor_config
from siac.domain import SensorBand, SensorConfig
from siac.domain.spectral import RelativeSpectralResponse


def _normalize_rsrf_root(rsrf_root: Path | str | None) -> Path | None:
    if rsrf_root is None:
        return None
    return Path(rsrf_root).expanduser().resolve()


def load_band_rsrf(
    band: SensorBand,
    *,
    sensor_id: str,
    satellite_id: str,
    rsrf_root: Path | str | None = None,
) -> RelativeSpectralResponse:
    """Load a canonical RSRF for a SIAC band from the external RSRF catalog."""
    if band.rsrf_sensor_unit_id is None:
        raise ValueError(f"Band {band.name!r} does not define an RSRF sensor unit id")

    root = _normalize_rsrf_root(rsrf_root)
    response_definition = rsrf.load_response_definition(
        band.rsrf_sensor_unit_id,
        band.rsrf_band_id or band.name,
        band.rsrf_representation_variant or "band_average",
        root=root,
    )
    if hasattr(response_definition, "wavelength_nm") and hasattr(response_definition, "response"):
        curve = response_definition
        source_variant = getattr(response_definition, "source_variant", None)
    else:
        curve = rsrf.realize_curve(
            response_definition,
            source_variant=band.rsrf_representation_variant or None,
        )
        source_variant = getattr(curve, "source_variant", None)

    return RelativeSpectralResponse.from_tabulated(
        sensor_id=sensor_id,
        satellite_id=satellite_id,
        band_name=band.name,
        wavelengths_nm=np.asarray(curve.wavelength_nm, dtype=np.float32),
        response=np.asarray(curve.response, dtype=np.float32),
        source_id="RSRF",
        source_version=None if source_variant is None else str(source_variant),
    )


def load_sensor_config_with_rsrf(
    sensor_id: str,
    satellite_id: str,
    *,
    rsrf_root: Path | str | None = None,
    base_config: SensorConfig | None = None,
) -> SensorConfig:
    """Attach RSRF curves to a built-in sensor config."""
    resolved_base = base_config or get_sensor_config(sensor_id, satellite_id)
    bands: list[SensorBand] = []

    for band in resolved_base.bands:
        band_rsrf = load_band_rsrf(
            band,
            sensor_id=resolved_base.sensor_id,
            satellite_id=resolved_base.satellite_id,
            rsrf_root=rsrf_root,
        )
        bands.append(
            SensorBand(
                name=band.name,
                center_wavelength=float(
                    band_rsrf.effective_wavelength_nm or band.center_wavelength
                ),
                bandwidth=float(band_rsrf.fwhm_nm or band.bandwidth),
                resolution=band.resolution,
                band_index=band.band_index,
                rsrf_wavelengths_nm=band_rsrf.wavelengths_nm.copy(),
                rsrf_response=band_rsrf.response.copy(),
                rsrf_sensor_unit_id=band.rsrf_sensor_unit_id,
                rsrf_representation_variant=band.rsrf_representation_variant,
                rsrf_band_id=band.rsrf_band_id,
            )
        )

    return SensorConfig(
        sensor_id=resolved_base.sensor_id,
        satellite_id=resolved_base.satellite_id,
        bands=tuple(bands),
        default_ref_scale=resolved_base.default_ref_scale,
        default_ref_offset=resolved_base.default_ref_offset,
    )


load_band_srf_from_rsrf = load_band_rsrf
load_sensor_config_from_rsrf = load_sensor_config_with_rsrf


__all__ = [
    "load_band_rsrf",
    "load_sensor_config_with_rsrf",
    "load_band_srf_from_rsrf",
    "load_sensor_config_from_rsrf",
]
