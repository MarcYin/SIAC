"""Catalog of known official and user-supplied SRF sources."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from siac.srf.sources.sentinel2 import S2_DOCUMENTS_PAGE_URL, S2_MISSION_PAGE_URL

SRFSourceAccess = Literal["remote", "metadata", "local"]
SRFSpectralDefinition = Literal[
    "tabulated_rsr",
    "metadata_band_characterization",
    "user_provided",
]
SRFImplementationStatus = Literal["implemented", "planned"]


@dataclass(frozen=True)
class KnownSRFSource:
    """Catalog entry describing an SRF source and access pattern."""

    sensor_id: str
    satellite_id: str
    source_access: SRFSourceAccess
    spectral_definition: SRFSpectralDefinition
    source_name: str
    landing_page_url: str | None
    asset_url: str | None
    implementation_status: SRFImplementationStatus
    notes: str = ""


KNOWN_SRF_SOURCES: tuple[KnownSRFSource, ...] = (
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
    sources = KNOWN_SRF_SOURCES
    if sensor_id is not None:
        sources = tuple(source for source in sources if source.sensor_id == sensor_id)
    if satellite_id is not None:
        sources = tuple(source for source in sources if source.satellite_id == satellite_id)
    return tuple(sources)


__all__ = [
    "KNOWN_SRF_SOURCES",
    "KnownSRFSource",
    "SRFImplementationStatus",
    "SRFSourceAccess",
    "SRFSpectralDefinition",
    "list_known_srf_sources",
]
