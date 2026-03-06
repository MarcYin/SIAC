"""Generic SRF source registry and load helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Literal

from siac.core.srf_builders import (
    BandCharacterization,
    build_sensor_config_from_band_characterization,
    build_sensor_config_from_tabulated_srfs,
)
from siac.core.srf_source_sentinel2 import (
    S2_DOCUMENTS_PAGE_URL,
    S2_MISSION_PAGE_URL,
    load_sentinel2_sensor_config,
)

if TYPE_CHECKING:
    from collections.abc import Iterable

    from siac.core.types import SensorConfig

SRFSourceAccess = Literal["remote", "metadata", "local"]
SRFSpectralDefinition = Literal[
    "tabulated_rsr",
    "metadata_band_characterization",
    "user_provided",
]
SRFImplementationStatus = Literal["implemented", "planned"]


@dataclass(frozen=True)
class KnownSRFSource:
    """Catalog entry describing an official or user-supplied SRF source."""

    sensor_id: str
    satellite_id: str
    source_access: SRFSourceAccess
    spectral_definition: SRFSpectralDefinition
    source_name: str
    landing_page_url: str | None
    asset_url: str | None
    implementation_status: SRFImplementationStatus
    notes: str = ""


_KNOWN_SRF_SOURCES: tuple[KnownSRFSource, ...] = (
    KnownSRFSource(
        sensor_id="MSI",
        satellite_id="S2A",
        source_access="remote",
        spectral_definition="tabulated_rsr",
        source_name="SentiWiki Sentinel-2 Spectral Response Functions",
        landing_page_url=S2_MISSION_PAGE_URL,
        asset_url=S2_DOCUMENTS_PAGE_URL,
        implementation_status="implemented",
        notes="Official ESA/Copernicus workbook linked from the Sentinel-2 mission documents page.",
    ),
    KnownSRFSource(
        sensor_id="MSI",
        satellite_id="S2B",
        source_access="remote",
        spectral_definition="tabulated_rsr",
        source_name="SentiWiki Sentinel-2 Spectral Response Functions",
        landing_page_url=S2_MISSION_PAGE_URL,
        asset_url=S2_DOCUMENTS_PAGE_URL,
        implementation_status="implemented",
        notes="Official ESA/Copernicus workbook linked from the Sentinel-2 mission documents page.",
    ),
    KnownSRFSource(
        sensor_id="MSI",
        satellite_id="S2C",
        source_access="remote",
        spectral_definition="tabulated_rsr",
        source_name="SentiWiki Sentinel-2 Spectral Response Functions",
        landing_page_url=S2_MISSION_PAGE_URL,
        asset_url=S2_DOCUMENTS_PAGE_URL,
        implementation_status="implemented",
        notes="Official ESA/Copernicus workbook linked from the Sentinel-2 mission documents page.",
    ),
    KnownSRFSource(
        sensor_id="OLI",
        satellite_id="L8",
        source_access="remote",
        spectral_definition="tabulated_rsr",
        source_name="USGS Landsat Spectral Characteristics Viewer",
        landing_page_url="https://landsat.usgs.gov/spectral-characteristics-viewer",
        asset_url="https://landsat.usgs.gov/spectral-characteristics-viewer",
        implementation_status="planned",
        notes="Official USGS viewer/export endpoint for Landsat 8 OLI spectral response functions.",
    ),
    KnownSRFSource(
        sensor_id="OLI2",
        satellite_id="L9",
        source_access="remote",
        spectral_definition="tabulated_rsr",
        source_name="USGS Landsat Spectral Characteristics Viewer",
        landing_page_url="https://landsat.usgs.gov/spectral-characteristics-viewer",
        asset_url="https://landsat.usgs.gov/spectral-characteristics-viewer",
        implementation_status="planned",
        notes="Official USGS viewer/export endpoint for Landsat 9 OLI-2 spectral response functions.",
    ),
    KnownSRFSource(
        sensor_id="MODIS",
        satellite_id="Terra",
        source_access="remote",
        spectral_definition="tabulated_rsr",
        source_name="NASA Ocean Color RSR tables",
        landing_page_url="https://oceancolor.gsfc.nasa.gov/resources/docs/rsr_tables/",
        asset_url="https://oceancolor.gsfc.nasa.gov/resources/docs/rsr_tables/",
        implementation_status="planned",
        notes="Official NASA Ocean Color RSR tables include MODIS/Terra response functions.",
    ),
    KnownSRFSource(
        sensor_id="MODIS",
        satellite_id="Aqua",
        source_access="remote",
        spectral_definition="tabulated_rsr",
        source_name="NASA Ocean Color RSR tables",
        landing_page_url="https://oceancolor.gsfc.nasa.gov/resources/docs/rsr_tables/",
        asset_url="https://oceancolor.gsfc.nasa.gov/resources/docs/rsr_tables/",
        implementation_status="planned",
        notes="Official NASA Ocean Color RSR tables include MODIS/Aqua response functions.",
    ),
    KnownSRFSource(
        sensor_id="VIIRS",
        satellite_id="SNPP",
        source_access="remote",
        spectral_definition="tabulated_rsr",
        source_name="NASA Ocean Color RSR tables",
        landing_page_url="https://oceancolor.gsfc.nasa.gov/resources/docs/rsr_tables/",
        asset_url="https://oceancolor.gsfc.nasa.gov/resources/docs/rsr_tables/",
        implementation_status="planned",
        notes="Official NASA Ocean Color RSR tables include VIIRS/SNPP response functions.",
    ),
    KnownSRFSource(
        sensor_id="VIIRS",
        satellite_id="NOAA-20",
        source_access="remote",
        spectral_definition="tabulated_rsr",
        source_name="NOAA STAR ICVS VIIRS spectral response pages",
        landing_page_url="https://www.star.nesdis.noaa.gov/icvs/status_N20_VIIRS.php",
        asset_url="https://www.star.nesdis.noaa.gov/icvs/status_N20_VIIRS.php",
        implementation_status="planned",
        notes="NOAA STAR publishes NOAA-20 VIIRS spectral response plots and supporting files.",
    ),
    KnownSRFSource(
        sensor_id="VIIRS",
        satellite_id="NOAA-21",
        source_access="remote",
        spectral_definition="tabulated_rsr",
        source_name="NOAA STAR ICVS VIIRS spectral response pages",
        landing_page_url="https://www.star.nesdis.noaa.gov/icvs/status_N21_VIIRS.php",
        asset_url="https://www.star.nesdis.noaa.gov/icvs/status_N21_VIIRS.php",
        implementation_status="planned",
        notes="NOAA STAR publishes NOAA-21 VIIRS spectral response plots and supporting files.",
    ),
    KnownSRFSource(
        sensor_id="OLCI",
        satellite_id="S3A",
        source_access="remote",
        spectral_definition="tabulated_rsr",
        source_name="Sentinel-3 OLCI instrument performance documents",
        landing_page_url="https://sentiwiki.copernicus.eu/web/s3-olci-instrument-performance",
        asset_url="https://sentiwiki.copernicus.eu/web/s3-olci-instrument-performance",
        implementation_status="planned",
        notes="Official Sentinel-3 OLCI instrument documents link pre-flight spectral response products for S3A/S3B.",
    ),
    KnownSRFSource(
        sensor_id="OLCI",
        satellite_id="S3B",
        source_access="remote",
        spectral_definition="tabulated_rsr",
        source_name="Sentinel-3 OLCI instrument performance documents",
        landing_page_url="https://sentiwiki.copernicus.eu/web/s3-olci-instrument-performance",
        asset_url="https://sentiwiki.copernicus.eu/web/s3-olci-instrument-performance",
        implementation_status="planned",
        notes="Official Sentinel-3 OLCI instrument documents link pre-flight spectral response products for S3A/S3B.",
    ),
    KnownSRFSource(
        sensor_id="SLSTR",
        satellite_id="S3A",
        source_access="remote",
        spectral_definition="tabulated_rsr",
        source_name="Sentinel-3 SLSTR instrument performance documents",
        landing_page_url="https://sentiwiki.copernicus.eu/web/s3-slstr-instrument-performance",
        asset_url="https://sentiwiki.copernicus.eu/web/s3-slstr-instrument-performance",
        implementation_status="planned",
        notes="Official Sentinel-3 SLSTR calibration and spectral response material lives under SentiWiki instrument performance pages.",
    ),
    KnownSRFSource(
        sensor_id="PRISMA_HSI",
        satellite_id="PRISMA",
        source_access="metadata",
        spectral_definition="metadata_band_characterization",
        source_name="ASI PRISMA mission products",
        landing_page_url="https://www.asi.it/en/earth-science/prisma/",
        asset_url=None,
        implementation_status="planned",
        notes="Treat PRISMA as metadata-driven unless a stable official tabulated RSR file is published; current public path is mission/product metadata.",
    ),
    KnownSRFSource(
        sensor_id="ENMAP_HSI",
        satellite_id="EnMAP",
        source_access="metadata",
        spectral_definition="metadata_band_characterization",
        source_name="EnMAP Level-1/2 product specification",
        landing_page_url="https://www.enmap.org/data/doc/EN-PCV-ICD-2009-2_HSI_Product_Specification_Level1_Level2.pdf",
        asset_url="https://www.enmap.org/data/doc/EN-PCV-ICD-2009-2_HSI_Product_Specification_Level1_Level2.pdf",
        implementation_status="planned",
        notes="Official EnMAP products expose center wavelength and FWHM per band in metadata.",
    ),
    KnownSRFSource(
        sensor_id="EMIT",
        satellite_id="ISS",
        source_access="metadata",
        spectral_definition="metadata_band_characterization",
        source_name="EMIT L2A reflectance ATBD",
        landing_page_url="https://lpdaac.usgs.gov/documents/2147/EMIT_L2A-RFL_ATBD_V1.pdf",
        asset_url="https://lpdaac.usgs.gov/documents/2147/EMIT_L2A-RFL_ATBD_V1.pdf",
        implementation_status="planned",
        notes="Official EMIT granules carry wavelengths and FWHM arrays in metadata; use metadata-driven spectral definitions.",
    ),
    KnownSRFSource(
        sensor_id="DOVE",
        satellite_id="PS2",
        source_access="remote",
        spectral_definition="tabulated_rsr",
        source_name="Planet RSR support article",
        landing_page_url="https://support.planet.com/hc/en-us/articles/4411132050451-How-Can-I-Access-Relative-Spectral-Responses-RSRs-",
        asset_url="https://support.planet.com/hc/en-us/articles/4411132050451-How-Can-I-Access-Relative-Spectral-Responses-RSRs-",
        implementation_status="planned",
        notes="Official Planet support article provides downloadable PlanetScope RSRs.",
    ),
    KnownSRFSource(
        sensor_id="SUPERDOVE",
        satellite_id="PS2.SD",
        source_access="remote",
        spectral_definition="tabulated_rsr",
        source_name="Planet RSR support article",
        landing_page_url="https://support.planet.com/hc/en-us/articles/4411132050451-How-Can-I-Access-Relative-Spectral-Responses-RSRs-",
        asset_url="https://support.planet.com/hc/en-us/articles/4411132050451-How-Can-I-Access-Relative-Spectral-Responses-RSRs-",
        implementation_status="planned",
        notes="Official Planet support article provides downloadable SuperDove / PS2.SD RSRs.",
    ),
    KnownSRFSource(
        sensor_id="SUPERDOVE",
        satellite_id="PSB.SD",
        source_access="remote",
        spectral_definition="tabulated_rsr",
        source_name="Planet RSR support article",
        landing_page_url="https://support.planet.com/hc/en-us/articles/4411132050451-How-Can-I-Access-Relative-Spectral-Responses-RSRs-",
        asset_url="https://support.planet.com/hc/en-us/articles/4411132050451-How-Can-I-Access-Relative-Spectral-Responses-RSRs-",
        implementation_status="planned",
        notes="Official Planet support article provides downloadable PlanetScope SuperDove / PSB.SD RSRs.",
    ),
)


def list_known_srf_sources(
    *,
    sensor_id: str | None = None,
    satellite_id: str | None = None,
) -> tuple[KnownSRFSource, ...]:
    """Return catalogued SRF sources known to SIAC."""
    sources: Iterable[KnownSRFSource] = _KNOWN_SRF_SOURCES
    if sensor_id is not None:
        sources = (source for source in sources if source.sensor_id == sensor_id)
    if satellite_id is not None:
        sources = (source for source in sources if source.satellite_id == satellite_id)
    return tuple(sources)


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
    del metadata_path, cache_dir, refresh
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
    "KnownSRFSource",
    "build_sensor_config_from_band_characterization",
    "build_sensor_config_from_tabulated_srfs",
    "list_known_srf_sources",
    "load_sensor_config_from_local_srf_file",
    "load_sensor_config_from_metadata_srf",
    "load_sensor_config_from_remote_srf",
    "load_sensor_config_from_srf",
]
