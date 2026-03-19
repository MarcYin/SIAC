"""High-level SRF loading entry points."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from siac.srf.builders import (
    BandCharacterization,
    build_sensor_config_from_band_characterization,
    build_sensor_config_from_tabulated_srfs,
)
from siac.srf.catalog import SRFSourceAccess, list_known_srf_sources
from siac.srf.metadata.emit import load_emit_sensor_config_from_metadata
from siac.srf.metadata.enmap import load_enmap_sensor_config_from_metadata
from siac.srf.metadata.prisma import load_prisma_sensor_config_from_metadata
from siac.srf.sources.landsat import load_landsat_sensor_config
from siac.srf.sources.planet import load_planet_sensor_config
from siac.srf.sources.sentinel2 import load_sentinel2_sensor_config

if TYPE_CHECKING:
    from siac.domain import SensorConfig


def _raise_no_loader(
    *,
    sensor_id: str,
    satellite_id: str,
    access_kind: SRFSourceAccess,
) -> None:
    matches = list_known_srf_sources(sensor_id=sensor_id, satellite_id=satellite_id)
    if matches:
        hints = "; ".join(
            f"{match.source_name} [{match.source_access}, {match.implementation_status}]"
            for match in matches
        )
        raise NotImplementedError(
            f"No {access_kind} SRF loader is implemented for {sensor_id}/{satellite_id}. "
            f"Known source catalog entries: {hints}"
        )
    raise NotImplementedError(
        f"No SRF source is catalogued for {sensor_id}/{satellite_id}."
    )


def load_sensor_config_from_remote_srf(
    sensor_id: str,
    satellite_id: str,
    *,
    cache_dir: str | Path | None = None,
    refresh: bool = False,
) -> SensorConfig:
    """Load a SensorConfig from an official remote SRF source."""
    if sensor_id == "MSI" and satellite_id in {"S2A", "S2B", "S2C"}:
        return load_sentinel2_sensor_config(
            satellite_id,
            cache_dir=cache_dir,
            refresh=refresh,
        )
    if (sensor_id, satellite_id) in {("OLI", "L8"), ("OLI2", "L9")}:
        return load_landsat_sensor_config(
            sensor_id,
            satellite_id,
            cache_dir=cache_dir,
            refresh=refresh,
        )
    if sensor_id in {"DOVE", "SUPERDOVE"}:
        return load_planet_sensor_config(
            sensor_id,
            satellite_id,
            cache_dir=cache_dir,
            refresh=refresh,
        )
    _raise_no_loader(
        sensor_id=sensor_id,
        satellite_id=satellite_id,
        access_kind="remote",
    )


def load_sensor_config_from_local_srf_file(
    sensor_id: str,
    satellite_id: str,
    local_path: str | Path,
    *,
    cache_dir: str | Path | None = None,
    refresh: bool = False,
) -> SensorConfig:
    """Load a SensorConfig from a user-provided local SRF file."""
    local_path = Path(local_path)
    if sensor_id == "MSI" and satellite_id in {"S2A", "S2B", "S2C"}:
        return load_sentinel2_sensor_config(
            satellite_id,
            cache_dir=cache_dir,
            refresh=refresh,
            workbook_path=local_path,
        )
    if (sensor_id, satellite_id) in {("OLI", "L8"), ("OLI2", "L9")}:
        return load_landsat_sensor_config(
            sensor_id,
            satellite_id,
            cache_dir=cache_dir,
            refresh=refresh,
            local_path=local_path,
        )
    if sensor_id in {"DOVE", "SUPERDOVE"}:
        return load_planet_sensor_config(
            sensor_id,
            satellite_id,
            cache_dir=cache_dir,
            refresh=refresh,
            local_path=local_path,
        )
    _raise_no_loader(
        sensor_id=sensor_id,
        satellite_id=satellite_id,
        access_kind="local",
    )


def load_sensor_config_from_metadata_srf(
    sensor_id: str,
    satellite_id: str,
    metadata_path: str | Path,
    *,
    cache_dir: str | Path | None = None,
    refresh: bool = False,
) -> SensorConfig:
    """Load a SensorConfig from metadata-derived spectral characterization."""
    del cache_dir, refresh
    metadata_path = Path(metadata_path)
    if (sensor_id, satellite_id) == ("ENMAP_HSI", "EnMAP"):
        return load_enmap_sensor_config_from_metadata(metadata_path)
    if (sensor_id, satellite_id) == ("EMIT", "ISS"):
        return load_emit_sensor_config_from_metadata(metadata_path)
    if (sensor_id, satellite_id) == ("PRISMA_HSI", "PRISMA"):
        return load_prisma_sensor_config_from_metadata(metadata_path)
    _raise_no_loader(
        sensor_id=sensor_id,
        satellite_id=satellite_id,
        access_kind="metadata",
    )


def load_sensor_config_from_srf(
    sensor_id: str,
    satellite_id: str,
    *,
    cache_dir: str | Path | None = None,
    refresh: bool = False,
    metadata_path: str | Path | None = None,
    local_path: str | Path | None = None,
) -> SensorConfig:
    """Load a SensorConfig using the best available SRF source for a platform."""
    if local_path is not None:
        return load_sensor_config_from_local_srf_file(
            sensor_id,
            satellite_id,
            local_path,
            cache_dir=cache_dir,
            refresh=refresh,
        )
    if metadata_path is not None:
        return load_sensor_config_from_metadata_srf(
            sensor_id,
            satellite_id,
            metadata_path,
            cache_dir=cache_dir,
            refresh=refresh,
        )
    return load_sensor_config_from_remote_srf(
        sensor_id,
        satellite_id,
        cache_dir=cache_dir,
        refresh=refresh,
    )


__all__ = [
    "BandCharacterization",
    "build_sensor_config_from_band_characterization",
    "build_sensor_config_from_tabulated_srfs",
    "load_sensor_config_from_local_srf_file",
    "load_sensor_config_from_metadata_srf",
    "load_sensor_config_from_remote_srf",
    "load_sensor_config_from_srf",
]
