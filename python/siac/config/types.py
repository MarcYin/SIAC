"""Shared literals and coercion helpers for SIAC configuration models."""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

DEFAULT_LUT_URL = (
    "https://gws-access.jasmin.ac.uk/public/nceo_isp/"
    "libradtran_continental_average_lut_1nm.zarr.zip"
)


class SIACEnum(str, Enum):
    """String enum with predictable string formatting on Python 3.10."""

    def __str__(self) -> str:
        return str(self.value)


class SensorName(SIACEnum):
    S2 = "s2"
    AUTO = "auto"


class LogLevel(SIACEnum):
    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"


class AtmoProviderKind(SIACEnum):
    CAMS = "cams"
    MERRA2 = "merra2"
    MCD19 = "mcd19"
    VNP19 = "vnp19"


class BRDFProviderKind(SIACEnum):
    MCD43 = "mcd43"
    MCD43_GEE = "mcd43_gee"
    VNP43 = "vnp43"
    MCD19 = "mcd19"


class MonthlyCompositeProviderKind(SIACEnum):
    GENERATED_BRDF = "generated_brdf"
    USER_CALLABLE = "user_callable"
    PREPARED_STORE = "prepared_store"
    BESTPIXEL = "bestpixel"


class AtmosphericParameterName(SIACEnum):
    AOT = "aot"
    TCWV = "tcwv"
    TCO3 = "tco3"


class SurfacePriorMethod(SIACEnum):
    KERNEL_MODEL = "kernel_model"
    WHITTAKER = "whittaker"
    MONTHLY_DATABASE = "monthly_database"
    BESTPIXEL = "bestpixel"


class SolverMethod(SIACEnum):
    """Atmospheric solver method.

    ``multigrid`` is the Bayesian multi-grid L-BFGS-B inversion (default).
    ``surface_driven`` sweeps a fixed AOT axis at the prior TCWV, scoring the
    surface-prior mismatch per node and picking the arg-min (cheap, robust;
    delegates the per-pixel cost cube to the Rust kernel).
    """

    MULTIGRID = "multigrid"
    SURFACE_DRIVEN = "surface_driven"


class MonthlyDatabaseResolutionPolicy(SIACEnum):
    PROVIDER_OR_COARSER = "provider_or_coarser"
    AEROSOL = "aerosol"


class TemporalInterpolation(SIACEnum):
    NEAREST = "nearest"
    LINEAR = "linear"


class S2Backend(SIACEnum):
    CDSE = "cdse"
    GCS = "gcs"
    LOCAL = "local"


class S2ProcessingLevel(SIACEnum):
    L1C = "L1C"
    L2A = "L2A"


class SixSMode(SIACEnum):
    DIRECT = "direct"
    SCENE_LUT = "scene_lut"
    AUTO = "auto"


class SixSParallelBackend(SIACEnum):
    OPENMP = "openmp"
    WORKER_LIBRARIES = "worker_libraries"


class SixSBuildProfile(SIACEnum):
    RELEASE = "release"
    PARITY = "parity"


class SixSAtmosphericColumnsMode(SIACEnum):
    INPUT_COLUMNS = "input_columns"
    PROFILE_DEFAULT = "profile_default"


class SixSAtmosphericProfile(SIACEnum):
    NO_GAS = "no_gas"
    TROPICAL = "tropical"
    MIDLATITUDE_SUMMER = "midlatitude_summer"
    MIDLATITUDE_WINTER = "midlatitude_winter"
    SUBARCTIC_SUMMER = "subarctic_summer"
    SUBARCTIC_WINTER = "subarctic_winter"
    US_STANDARD_62 = "us_standard_62"
    AUTO_LATITUDE_DATE = "auto_latitude_date"
    USER_WATER_OZONE = "user_water_ozone"
    USER_PROFILE = "user_profile"


class SixSAerosolProfile(SIACEnum):
    NONE = "none"
    CONTINENTAL = "continental"
    MARITIME = "maritime"
    URBAN = "urban"
    DESERT = "desert"
    BIOMASS_BURNING = "biomass_burning"
    STRATOSPHERIC = "stratospheric"
    USER_MIXTURE = "user_mixture"
    MULTIMODAL_LOG_NORMAL = "multimodal_log_normal"
    MODIFIED_GAMMA = "modified_gamma"
    JUNGE_POWER_LAW = "junge_power_law"
    SUN_PHOTOMETER = "sun_photometer"
    LAYERED_PROFILE = "layered_profile"
    USER_MODEL = "user_model"


class RTAerosolProfile(SIACEnum):
    NONE = "none"
    CONTINENTAL = "continental"
    CONTINENTAL_AVERAGE = "continental_average"
    MARITIME = "maritime"
    URBAN = "urban"
    DESERT = "desert"
    BIOMASS_BURNING = "biomass_burning"
    STRATOSPHERIC = "stratospheric"
    USER_MIXTURE = "user_mixture"
    MULTIMODAL_LOG_NORMAL = "multimodal_log_normal"
    MODIFIED_GAMMA = "modified_gamma"
    JUNGE_POWER_LAW = "junge_power_law"
    SUN_PHOTOMETER = "sun_photometer"
    LAYERED_PROFILE = "layered_profile"
    USER_MODEL = "user_model"


class SixSAerosolLayerType(SIACEnum):
    NONE = "none"
    CONTINENTAL = "continental"
    MARITIME = "maritime"
    URBAN = "urban"
    DESERT = "desert"
    BIOMASS_BURNING = "biomass_burning"
    STRATOSPHERIC = "stratospheric"


class SixSBuiltinReflectance(SIACEnum):
    GREEN_VEGETATION = "green_vegetation"
    CLEAR_WATER = "clear_water"
    SAND = "sand"
    LAKE_WATER = "lake_water"


class SixSSurfaceReflectanceKind(SIACEnum):
    CONSTANT = "constant"
    BUILT_IN = "built_in"
    SPECTRUM = "spectrum"


class SixSSurfaceMode(SIACEnum):
    HOMOGENEOUS_LAMBERTIAN = "homogeneous_lambertian"
    HETEROGENEOUS_LAMBERTIAN = "heterogeneous_lambertian"
    HOMOGENEOUS_BRDF = "homogeneous_brdf"


class SixSBRDFModel(SIACEnum):
    USER_DEFINED = "user_defined"
    HAPKE = "hapke"
    VERSTRAETE = "verstraete"
    ROUJEAN = "roujean"
    WALTHALL = "walthall"
    MINNAERT = "minnaert"
    OCEAN = "ocean"
    IAQUINTA_PINTY = "iaquinta_pinty"
    RAHMAN = "rahman"
    KUUSK = "kuusk"
    MODIS = "modis"
    ROSS_LI_MAIGNAN = "ross_li_maignan"


class SixSAtmosphericCorrectionMode(SIACEnum):
    NONE = "none"
    LAMBERTIAN_REFLECTANCE = "lambertian_reflectance"
    LAMBERTIAN_RADIANCE = "lambertian_radiance"
    BRDF_REFLECTANCE = "brdf_reflectance"
    BRDF_RADIANCE = "brdf_radiance"


class RTBackend(SIACEnum):
    EMULATOR = "emulator"
    LUT = "lut"
    SIXS = "sixs"
    LIBRADTRAN = "libradtran"


class LUTInterpolationMethod(SIACEnum):
    LINEAR = "linear"
    NEAREST = "nearest"
    CUBIC = "cubic"


class SolverStageInitialState(SIACEnum):
    PRIOR = "prior"
    PREVIOUS = "previous"


class FixedAtmosphericParameter(SIACEnum):
    NONE = "none"
    AOT = "aot"
    TCWV = "tcwv"


class CloudMaskMode(SIACEnum):
    AUTO = "auto"
    EXTERNAL_FILE = "external_file"
    USER_CALLABLE = "user_callable"
    NONE = "none"


class CloudMaskProvider(SIACEnum):
    OMNICLOUDMASK = "omnicloudmask"


class ResolutionPolicy(SIACEnum):
    AUTO = "auto"
    FORCE = "force"


class ExecutionBackend(SIACEnum):
    THREAD = "thread"
    DASK = "dask"


class OutputFormat(SIACEnum):
    GEOTIFF = "geotiff"
    COG = "cog"
    ZARR = "zarr"
    NETCDF = "netcdf"


class OutputCompression(SIACEnum):
    DEFLATE = "deflate"
    LZW = "lzw"
    ZSTD = "zstd"
    NONE = "none"


class BOADType(SIACEnum):
    FLOAT32 = "float32"
    FLOAT64 = "float64"
    UINT16 = "uint16"


DEFAULT_SIXS_SOURCE_URL = "https://salsa.umd.edu/files/6S/6sV2.1.tar"

#: libRadtran source + auxiliary-data archives compiled from source for the
#: ``backend = "libradtran"`` RT model. 2.0.6 is the only version the upstream
#: server still serves (older version URLs redirect to it); it reproduces the
#: same RT preset the remote ``libradtran_continental_average_lut_1nm.zarr.zip``
#: LUT encodes. ``reptran`` (molecular absorption band model) and the OPAC
#: optical properties (``continental_average`` aerosol) are separate downloads
#: that must be merged into the source tree's ``data/`` directory before
#: building (they are not bundled in the core tarball).
DEFAULT_LIBRADTRAN_SOURCE_URL = "https://www.libradtran.org/download/libRadtran-2.0.6.tar.gz"
DEFAULT_LIBRADTRAN_REPTRAN_URL = (
    "https://www.libradtran.org/lib/exe/fetch.php?media=download:reptran_2024_all.tar.gz"
)
DEFAULT_LIBRADTRAN_OPTPROP_URL = (
    "https://www.libradtran.org/lib/exe/fetch.php?media=download:optprop_v2.1.tar.gz"
)

_REMOTE_URI_SCHEMES = {"http", "https", "s3", "file", "gs"}


def _coerce_pathlike(value: Any) -> Path | None:
    if value is None or isinstance(value, Path):
        return value
    text = str(value).strip()
    if not text:
        return None
    return Path(text).expanduser()


def _coerce_path_or_url(value: Any) -> str | Path | None:
    if value is None or isinstance(value, Path):
        return value
    text = str(value).strip()
    if not text:
        return None
    if text.startswith("/vsi"):
        return text
    if urlparse(text).scheme.lower() in _REMOTE_URI_SCHEMES:
        return text
    return Path(text).expanduser()


def _coerce_float_tuple(value: Any) -> tuple[float, ...] | None:
    if value is None:
        return None
    if isinstance(value, tuple):
        return tuple(float(item) for item in value)
    if isinstance(value, list):
        return tuple(float(item) for item in value)
    return (float(value),)


def _coerce_float_matrix(value: Any) -> tuple[tuple[float, ...], ...] | None:
    if value is None:
        return None
    rows = []
    for row in value:
        rows.append(tuple(float(item) for item in row))
    return tuple(rows)


def _coerce_float_series_with_broadcast(
    value: Any,
    *,
    target_length: int,
) -> tuple[float, ...] | None:
    series = _coerce_float_tuple(value)
    if series is None:
        return None
    if len(series) == 1:
        return tuple(series[0] for _ in range(target_length))
    return series


def _normalize_aerosol_mixture_payload(value: Any) -> tuple[float, float, float, float] | None:
    if value is None:
        return None
    series = _coerce_float_tuple(value)
    if series is None:
        return None
    if len(series) != 4:
        raise ValueError(
            "Aerosol mixtures must contain four components: dust, water, oceanic, soot."
        )
    return (
        float(series[0]),
        float(series[1]),
        float(series[2]),
        float(series[3]),
    )
