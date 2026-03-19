"""Helpers for building SensorConfig objects from SRF sources."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from siac.core.types import SensorBand, SensorConfig
from siac.srf.types import SpectralResponseFunction

if TYPE_CHECKING:
    from collections.abc import Mapping


@dataclass(frozen=True)
class BandCharacterization:
    """Per-band spectral characterization parsed from scene metadata."""

    band_name: str
    center_wavelength_nm: float
    fwhm_nm: float


def build_sensor_config_from_tabulated_srfs(
    base_config: SensorConfig,
    srfs: Mapping[str, SpectralResponseFunction],
) -> SensorConfig:
    """Attach canonical tabulated SRFs to a SensorConfig."""
    bands: list[SensorBand] = []
    for band in base_config.bands:
        srf = srfs.get(band.name)
        if srf is None:
            raise KeyError(
                f"Missing SRF for {base_config.satellite_id} band {band.name}"
            )
        bands.append(
            SensorBand(
                name=band.name,
                center_wavelength=float(
                    srf.effective_wavelength_nm or band.center_wavelength
                ),
                bandwidth=float(srf.fwhm_nm or band.bandwidth),
                resolution=band.resolution,
                band_index=band.band_index,
                srf_wavelengths_nm=srf.wavelengths_nm.copy(),
                srf_response=srf.response.copy(),
                rsrf_sensor_unit_id=band.rsrf_sensor_unit_id,
                rsrf_representation_variant=band.rsrf_representation_variant,
                rsrf_band_id=band.rsrf_band_id,
            )
        )
    return SensorConfig(
        sensor_id=base_config.sensor_id,
        satellite_id=base_config.satellite_id,
        bands=tuple(bands),
        default_ref_scale=base_config.default_ref_scale,
        default_ref_offset=base_config.default_ref_offset,
    )


def _gaussian_srf_from_band_characterization(
    *,
    sensor_id: str,
    satellite_id: str,
    characterization: BandCharacterization,
    spacing_nm: float,
    support_sigma_factor: float,
) -> SpectralResponseFunction:
    if spacing_nm <= 0.0:
        raise ValueError("spacing_nm must be positive")
    if support_sigma_factor <= 0.0:
        raise ValueError("support_sigma_factor must be positive")

    sigma = characterization.fwhm_nm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    half_window = max(float(spacing_nm), float(support_sigma_factor) * float(sigma))
    start = characterization.center_wavelength_nm - half_window
    stop = characterization.center_wavelength_nm + half_window
    wavelengths = np.arange(
        start,
        stop + 0.5 * spacing_nm,
        spacing_nm,
        dtype=np.float32,
    )
    response = np.exp(
        -0.5 * np.square(
            (wavelengths - characterization.center_wavelength_nm) / sigma
        )
    ).astype(np.float32)
    if response.size >= 2:
        response[0] = 0.0
        response[-1] = 0.0
    return SpectralResponseFunction.from_tabulated(
        sensor_id=sensor_id,
        satellite_id=satellite_id,
        band_name=characterization.band_name,
        wavelengths_nm=wavelengths,
        response=response,
        source_id="metadata-band-characterization",
    )


def build_sensor_config_from_band_characterization(
    base_config: SensorConfig,
    characterization: Mapping[str, BandCharacterization],
    *,
    spacing_nm: float = 1.0,
    support_sigma_factor: float = 4.0,
) -> SensorConfig:
    """Build tabulated Gaussian SRFs from per-band center wavelength and FWHM."""
    srfs = {
        band_name: _gaussian_srf_from_band_characterization(
            sensor_id=base_config.sensor_id,
            satellite_id=base_config.satellite_id,
            characterization=band_meta,
            spacing_nm=spacing_nm,
            support_sigma_factor=support_sigma_factor,
        )
        for band_name, band_meta in characterization.items()
    }
    return build_sensor_config_from_tabulated_srfs(base_config, srfs)
