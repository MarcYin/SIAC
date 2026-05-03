"""Shared helpers for built-in sensor catalog entries."""

from siac.domain.sensors import SensorBand


def rsrf_band(
    name: str,
    center_wavelength: float,
    bandwidth: float,
    resolution: float,
    band_index: int,
    *,
    sensor_unit_id: str,
    representation_variant: str = "band_average",
    rsrf_band_id: str | None = None,
) -> SensorBand:
    """Construct a band with an attached RSRF repository identity."""
    return SensorBand(
        name=name,
        center_wavelength=center_wavelength,
        bandwidth=bandwidth,
        resolution=resolution,
        band_index=band_index,
        rsrf_sensor_unit_id=sensor_unit_id,
        rsrf_representation_variant=representation_variant,
        rsrf_band_id=rsrf_band_id or name,
    )
