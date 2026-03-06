"""Local repository for canonical spectral response functions."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable

    from siac.core.srf import SpectralResponseFunction


class SRFRepository:
    """Simple in-memory SRF repository keyed by platform and band."""

    def __init__(self, srfs: Iterable[SpectralResponseFunction]):
        self._srfs = {
            (srf.sensor_id, srf.satellite_id, srf.band_name): srf
            for srf in srfs
        }

    def get_band_srf(
        self,
        sensor_id: str,
        satellite_id: str,
        band_name: str,
    ) -> SpectralResponseFunction:
        """Return the SRF for a specific platform band."""
        key = (sensor_id, satellite_id, band_name)
        try:
            return self._srfs[key]
        except KeyError as exc:
            raise KeyError(
                f"SRF not found for sensor={sensor_id!r}, satellite={satellite_id!r}, band={band_name!r}"
            ) from exc

    def get_sensor_srfs(
        self,
        sensor_id: str,
        satellite_id: str,
    ) -> dict[str, SpectralResponseFunction]:
        """Return all SRFs for a single platform keyed by band name."""
        return {
            band_name: srf
            for (sid, sat_id, band_name), srf in self._srfs.items()
            if sid == sensor_id and sat_id == satellite_id
        }
