"""Adapter-layer helpers for external systems and backends."""

from siac.adapters.atmo import CAMSProvider, MCD19AODProvider, MERRA2Provider, VNP19AODProvider
from siac.adapters.auth import (
    AWSAuth,
    CDSAuth,
    CDSEAuth,
    CredentialManager,
    CredentialSpec,
    EarthdataAuth,
    GCSAuth,
)
from siac.adapters.brdf import (
    MCD19EarthAccessProvider,
    MCD43EarthAccessProvider,
    VNP43EarthAccessProvider,
)
from siac.adapters.earthdata import earthaccess_source_from_auth
from siac.adapters.earthdata_common import MODLAND_SINUSOIDAL_CRS, modland_tile_coords
from siac.adapters.output import ConfiguredOutputWriter
from siac.adapters.rsrf import load_band_srf_from_rsrf, load_sensor_config_from_rsrf
from siac.adapters.rt import build_rt_model
from siac.adapters.s2_backend import build_s2_backend
from siac.adapters.satellite import (
    BaseSatellitePreprocessor,
    Sentinel2Preprocessor,
    detect_sensor,
    get_preprocessor,
)

__all__ = [
    "CAMSProvider",
    "MERRA2Provider",
    "MCD19AODProvider",
    "VNP19AODProvider",
    "AWSAuth",
    "CDSAuth",
    "CDSEAuth",
    "CredentialManager",
    "CredentialSpec",
    "EarthdataAuth",
    "GCSAuth",
    "MCD43EarthAccessProvider",
    "VNP43EarthAccessProvider",
    "MCD19EarthAccessProvider",
    "MODLAND_SINUSOIDAL_CRS",
    "modland_tile_coords",
    "ConfiguredOutputWriter",
    "load_band_srf_from_rsrf",
    "load_sensor_config_from_rsrf",
    "BaseSatellitePreprocessor",
    "Sentinel2Preprocessor",
    "detect_sensor",
    "get_preprocessor",
    "earthaccess_source_from_auth",
    "build_rt_model",
    "build_s2_backend",
]
