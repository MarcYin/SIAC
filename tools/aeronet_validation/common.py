"""Shared paths, constants, and helpers for the AERONET validation experiment."""

from __future__ import annotations

import json
import logging
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

DEFAULT_DATA_ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")

#: Machine assets shared by every retrieval run. The LUT zip sits on the same
#: JASMIN GWS as the experiment data root, so prefer the local file and fall
#: back to the public HTTPS mirror elsewhere.
LOCAL_LUT_PATH = Path(
    "/gws/ssde/j25a/nceo_isp/public/libradtran_continental_average_lut_1nm.zarr.zip"
)
DEFAULT_DEM = (
    "/vsicurl/https://raw.githubusercontent.com/MarcYin/Copernicus_GLO_30_DEM_VRT/"
    "refs/heads/main/copernicus_GLO_30_dem.vrt"
)
DEFAULT_WATER_MASK = "https://zenodo.org/records/14899246/files/landWater2020.vrt?download=1"
_SLURM_ENV_VARS: Final[tuple[str, ...]] = (
    "SLURM_JOB_ID",
    "SLURM_ARRAY_JOB_ID",
    "SLURM_JOB_NAME",
)

#: Surface-prior approaches under comparison. Each entry is a nested config
#: overlay merged into the shared base config, so the solver/RT side stays
#: identical across approaches and the surface prior is the only variable.
APPROACHES: dict[str, dict[str, Any]] = {
    "kernel_model": {
        "algorithms": {"surface_prior": {"method": "kernel_model"}},
    },
    "whittaker": {
        "algorithms": {"surface_prior": {"method": "whittaker"}},
    },
    "monthly_database": {
        "algorithms": {"surface_prior": {"method": "monthly_database"}},
        "providers": {"monthly_composites": {"kind": "generated_brdf"}},
    },
    "monthly_database_spread_only": {
        "algorithms": {
            "surface_prior": {
                "method": "monthly_database",
                "monthly_database_filter": {"composite_uncertainty_scale": 0.0},
            }
        },
        "providers": {"monthly_composites": {"kind": "generated_brdf"}},
    },
    "monthly_database_spread_only_aot25": {
        "algorithms": {
            "surface_prior": {
                "method": "monthly_database",
                "monthly_database_filter": {"composite_uncertainty_scale": 0.0},
            },
            "solver": {"bounds": {"aot": (0.001, 2.5)}},
        },
        "providers": {"monthly_composites": {"kind": "generated_brdf"}},
    },
    "monthly_database_spread_only_aot25_legacy_resample": {
        "algorithms": {
            "surface_prior": {
                "method": "monthly_database",
                "monthly_database_filter": {"composite_uncertainty_scale": 0.0},
            },
            "solver": {
                "bounds": {"aot": (0.001, 2.5)},
                "water_mask_buffer_pixels": 0,
            },
        },
        "providers": {"monthly_composites": {"kind": "generated_brdf"}},
        "runtime": {"grid_resample_workers": 6},
    },
}

AERONET_WAVELENGTHS_NM = (440.0, 500.0, 675.0, 870.0)


@dataclass(frozen=True)
class ExperimentPaths:
    """Directory layout of the experiment under the data root."""

    root: Path

    @property
    def aeronet_dir(self) -> Path:
        return self.root / "aeronet"

    @property
    def aeronet_raw_dir(self) -> Path:
        return self.aeronet_dir / "raw"

    @property
    def aeronet_parsed_dir(self) -> Path:
        return self.aeronet_dir / "parsed"

    @property
    def sites_file(self) -> Path:
        return self.aeronet_dir / "sites.csv"

    @property
    def sites_with_data_file(self) -> Path:
        return self.aeronet_dir / "sites_with_data.csv"

    @property
    def matchup_dir(self) -> Path:
        return self.root / "matchups"

    @property
    def matchup_search_dir(self) -> Path:
        return self.matchup_dir / "search"

    @property
    def matchups_file(self) -> Path:
        return self.matchup_dir / "matchups.csv"

    @property
    def cache_dir(self) -> Path:
        return self.root / "cache"

    @property
    def runs_dir(self) -> Path:
        return self.root / "runs"

    @property
    def manifest_file(self) -> Path:
        return self.runs_dir / "manifest.csv"

    @property
    def results_dir(self) -> Path:
        return self.root / "results"

    @property
    def slurm_dir(self) -> Path:
        return self.root / "slurm"

    def run_dir(self, approach: str, matchup_id: str) -> Path:
        return self.runs_dir / approach / matchup_id

    def ensure(self) -> None:
        for path in (
            self.aeronet_raw_dir,
            self.aeronet_parsed_dir,
            self.matchup_search_dir,
            self.cache_dir,
            self.runs_dir,
            self.results_dir,
            self.slurm_dir,
        ):
            path.mkdir(parents=True, exist_ok=True)


def is_slurm_job() -> bool:
    """Return True when current process runs under a Slurm allocation."""
    return any(var in os.environ for var in _SLURM_ENV_VARS)


def require_slurm_execution(
    workload_name: str,
    *,
    allow_login: bool = False,
    suggestion: str = "Submit through the matching Slurm script instead of running locally.",
) -> None:
    """Guard against running expensive compute on login nodes.

    Some experiment workflows are designed to run under batch scheduling and can
    be costly in both CPU and memory. By default we fail fast when those jobs
    are started interactively, and ask for explicit override.
    """
    if allow_login or is_slurm_job():
        return
    raise SystemExit(f"Refusing to run {workload_name} on the login node.\n{suggestion}")


def setup_logging(level: str = "INFO", log_file: Path | None = None) -> None:
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s %(levelname)s:%(name)s:%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=handlers,
        force=True,
    )


def map_legacy_earthdata_env() -> None:
    """Map legacy ``Earthdata_user``/``Earthdata_pass`` env vars to the names
    earthaccess and SIAC read (``EARTHDATA_USERNAME``/``EARTHDATA_PASSWORD``)."""
    legacy_user = os.environ.get("Earthdata_user")
    legacy_pass = os.environ.get("Earthdata_pass")
    if legacy_user and "EARTHDATA_USERNAME" not in os.environ:
        os.environ["EARTHDATA_USERNAME"] = legacy_user
    if legacy_pass and "EARTHDATA_PASSWORD" not in os.environ:
        os.environ["EARTHDATA_PASSWORD"] = legacy_pass


def resolve_lut_path() -> str:
    if LOCAL_LUT_PATH.exists():
        return str(LOCAL_LUT_PATH)
    return (
        "https://gws-access.jasmin.ac.uk/public/nceo_isp/"
        "libradtran_continental_average_lut_1nm.zarr.zip"
    )


#: Daily global CAMS forecasts (`YYYY-MM-DD.nc`, 2015->present) mirrored on the
#: nceo_ard GWS. Reading these locally avoids the ADS request queue, which can
#: stall a retrieval for tens of minutes per scene.
LOCAL_CAMS_MIRROR = Path("/gws/ssde/j25b/nceo_ard/public/cams")
JASMIN_CAMS_MIRROR_URL = "https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/"


def resolve_cams_data_path(cache_dir: Path) -> str:
    if LOCAL_CAMS_MIRROR.is_dir():
        return str(LOCAL_CAMS_MIRROR)
    return str(cache_dir / "cams")


def matchup_id_for(site: str, mgrs_tile: str, sensing_time_compact: str) -> str:
    """Stable per-matchup identifier, safe as a directory name."""
    safe_site = site.replace("/", "-").replace(" ", "_")
    return f"{safe_site}__T{mgrs_tile.lstrip('T')}_{sensing_time_compact}"


def aod_550_from_angstrom(aod_reference: float, reference_nm: float, angstrom: float) -> float:
    return float(aod_reference * (550.0 / reference_nm) ** (-angstrom))


def aod_550_from_channels(aods_by_nm: dict[float, float]) -> float | None:
    """AOD at 550 nm via log-log interpolation over available AERONET channels.

    Prefers AOD_500nm + Angstrom-style two-point fits but degrades gracefully:
    a least-squares line in (ln wavelength, ln AOD) over every positive channel,
    evaluated at 550 nm. Returns None when fewer than two channels are usable.
    """
    points = [
        (math.log(wavelength), math.log(aod))
        for wavelength, aod in sorted(aods_by_nm.items())
        if aod is not None and aod > 0.0 and not math.isnan(aod)
    ]
    if len(points) < 2:
        return None
    n = float(len(points))
    sum_x = sum(x for x, _ in points)
    sum_y = sum(y for _, y in points)
    sum_xx = sum(x * x for x, _ in points)
    sum_xy = sum(x * y for x, y in points)
    denominator = n * sum_xx - sum_x * sum_x
    if denominator == 0.0:
        return None
    slope = (n * sum_xy - sum_x * sum_y) / denominator
    intercept = (sum_y - slope * sum_x) / n
    return float(math.exp(intercept + slope * math.log(550.0)))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str))


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())
