"""
Main SIAC class for atmospheric correction.

Orchestrates the full pipeline: preprocessing, prior retrieval,
aerosol solving, and atmospheric correction.
"""

from __future__ import annotations

import logging
from datetime import date, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from siac.core.aoi import AOI
from siac.core.auth import CredentialManager
from siac.core.config import SIACConfig
from siac.core.exceptions import DataNotFoundError
from siac.core.types import (
    AtmosphericState,
    GeometryAngles,
    ObservationBundle,
    SensorConfig,
    SolvedAtmosphere,
    SolverInputBundle,
)
from siac.correction import AtmosphericCorrector, CorrectionResult
from siac.io.earthaccess_source import EarthAccessSource
from siac.pipeline import (
    AtmoPriorFn,
    CorrectorFn,
    GridAssemblerFn,
    PreprocessorFn,
    SolverFn,
    SurfacePriorFn,
    run_pipeline,
)
from siac.priors.atmospheric import CAMSProvider
from siac.priors.surface.brdf_whittaker import BRDFWhittakerDeriver
from siac.priors.surface.kernel_model import KernelModelDeriver
from siac.priors.surface.swir_refine import (
    build_monthly_surface_prior_database,
    query_surface_prior_from_monthly_database,
    resample_geometry_for_surface_prior,
)
from siac.rt.emulator import TwoLayerNNEmulator
from siac.rt.lut import ZarrLUTBackend
from siac.satellite import detect_sensor, get_preprocessor
from siac.solver import MultiGridConfig, MultiGridSolver

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    import xarray as xr

    from siac.io.s2_data_source import S2Product, S2Query


def _earthaccess_source_from_auth(
    auth: CredentialManager | None,
    *,
    provider: str | None = None,
) -> EarthAccessSource:
    """Build an Earthaccess source that can authenticate non-interactively when creds exist."""
    if auth is None:
        return EarthAccessSource(provider=provider)
    return auth.earthdata().build_earthaccess_source(provider=provider)


def _select_surface_prior_bands(sensor_config: SensorConfig | None) -> list[Any]:
    """Return the sensor bands needed by the surface prior solver path."""
    if sensor_config is None:
        return list(range(1, 8))

    bands = sensor_config.select_bands_in_range(400.0, 520.0)
    if not bands:
        bands = list(sensor_config.bands[:2])
    return bands


def _select_visible_surface_prior_bands(sensor_config: SensorConfig) -> list[Any]:
    """Visible bands used as the Route-B retrieval target."""
    bands = [band for band in sensor_config.bands if 400.0 <= band.center_wavelength < 700.0]
    if bands:
        return bands
    return list(sensor_config.bands[: min(3, len(sensor_config.bands))])


def _select_route_b_query_bands(sensor_config: SensorConfig) -> list[Any]:
    """Return NIR, SWIR-1, SWIR-2 bands for Route-B database queries."""
    preferred_by_sensor = {
        "MSI": ("B08", "B11", "B12"),
        "OLI": ("B5", "B6", "B7"),
    }
    preferred_names = preferred_by_sensor.get(sensor_config.sensor_id, ())
    selected = [band for name in preferred_names for band in sensor_config.bands if band.name == name]
    if len(selected) == 3:
        return selected

    targets = (860.0, 1610.0, 2200.0)
    remaining = list(sensor_config.bands)
    resolved: list[Any] = []
    for target in targets:
        candidates = [band for band in remaining if band.center_wavelength >= 700.0] or remaining
        match = min(candidates, key=lambda band: abs(band.center_wavelength - target))
        resolved.append(match)
        remaining = [band for band in remaining if band.name != match.name]
    return resolved


class SIAC:
    """
    Sensor-Invariant Atmospheric Correction.

    Main entry point for atmospheric correction of satellite imagery.

    Example:
        >>> from siac import SIAC
        >>> siac = SIAC.from_config("config.yaml")
        >>> result = siac.process("/path/to/S2_SAFE/")
    """

    def __init__(self, config: SIACConfig):
        self.config = config
        self._auth = CredentialManager.from_config(config)
        self._preprocessor = None
        self._atmo_provider = None
        self._brdf_provider = None
        self._rt_model = None
        self._solver = None

    @classmethod
    def from_config(cls, config_path: str | Path) -> SIAC:
        """Create SIAC from configuration file."""
        config = SIACConfig.from_yaml(config_path)
        return cls(config)

    @classmethod
    def from_defaults(cls, sensor: str = "auto") -> SIAC:
        """Create SIAC with default configuration."""
        config = SIACConfig(sensor=sensor)
        return cls(config)

    def process(self, input_path: str | Path, output_path: str | Path | None = None) -> CorrectionResult:
        """
        Run full atmospheric correction pipeline.

        Delegates to :func:`run_pipeline` so that all pipeline
        validations, concurrency backends, and retry logic apply.

        Args:
            input_path: Path to satellite data
            output_path: Optional output path (default: input_path/BOA/)

        Returns:
            CorrectionResult with BOA reflectance
        """
        input_path = Path(input_path)
        logger.info(f"Processing: {input_path}")

        # Detect sensor and build a PreprocessorFn that returns ObservationBundle
        sensor_id = self.config.sensor if self.config.sensor != "auto" else detect_sensor(input_path)
        preprocessor = get_preprocessor(sensor_id)
        if isinstance(getattr(preprocessor, "config", None), dict):
            preprocessor.config = {
                **preprocessor.config,
                "cloud_mask": self.config.cloud_mask.model_dump(exclude={"user_callable"}),
            }
        sensor_config = preprocessor.sensor_config

        def _preprocessor_fn(path: Path, aoi=None) -> ObservationBundle:
            raw = preprocessor.preprocess(path)
            toa = raw["toa"]
            resolved_aoi = aoi or self._resolve_aoi(toa)
            return ObservationBundle(
                toa=toa,
                geometry=raw["geometry"],
                cloud_mask=raw["cloud_mask"],
                sensor_config=sensor_config,
                metadata=raw["metadata"],
                crs=str(resolved_aoi.crs),
                bounds=resolved_aoi.get_bounds(),
            )

        # Resolve remaining pipeline components
        atmo_provider = _resolve_atmo_provider(self.config, auth=self._auth)
        surface_prior_provider = _resolve_surface_prior_provider(self.config, auth=self._auth)
        grid_assembler = _resolve_grid_assembler()
        solver = _resolve_solver(self.config)
        corrector = _resolve_corrector(self.config)
        rt_model = self._get_rt_model(sensor_config)

        aoi = None  # let preprocessor resolve AOI internally

        result = run_pipeline(
            input_path,
            aoi,
            self.config,
            preprocessor=_preprocessor_fn,
            atmo_provider=atmo_provider,
            surface_prior_provider=surface_prior_provider,
            grid_assembler=grid_assembler,
            solver=solver,
            corrector=corrector,
            rt_model=rt_model,
        )

        # Save output if path provided
        if output_path is not None:
            self._save_output(result, output_path)

        logger.info(f"Complete. Mean AOT: {float(result.aot.mean()):.3f}")
        return result

    def _resolve_aoi(self, toa: xr.Dataset) -> AOI:
        """
        Resolve AOI from config or from the TOA extent.

        If config.aoi is set, creates AOI from GeoJSON/bounds/WKT.
        Otherwise, extracts AOI from the TOA dataset extent.
        """
        if self.config.aoi is not None:
            aoi_spec = self.config.aoi
            # Could be a file path, WKT, or bounds string
            if isinstance(aoi_spec, (list, tuple)) and len(aoi_spec) == 4:
                return AOI.from_bounds(tuple(aoi_spec))
            return AOI.from_geojson(aoi_spec)

        # Default: use TOA extent
        first_var = list(toa.data_vars)[0]
        return AOI.from_raster(toa[first_var])

    def _get_atmospheric_prior(self, aoi: AOI, metadata: dict) -> AtmosphericState:
        """Get atmospheric prior from configured provider, scoped to AOI."""
        if self._atmo_provider is None:
            provider_name = self.config.atmo_prior.provider
            if provider_name == "cams":
                self._atmo_provider = CAMSProvider(
                    self.config.atmo_prior.data_path,
                    temporal_interp=self.config.atmo_prior.temporal_interpolation == "linear",
                    download_missing=self.config.atmo_prior.download_missing,
                    auth=self._auth,
                    cache_dir=self.config.atmo_prior.cache_dir,
                )
            elif provider_name == "merra2":
                from siac.priors.atmospheric.merra2 import MERRA2Provider
                self._atmo_provider = MERRA2Provider(
                    cache_dir=self.config.atmo_prior.cache_dir,
                    source=_earthaccess_source_from_auth(self._auth),
                )
            elif provider_name == "mcd19":
                from siac.priors.atmospheric.mcd19_earthaccess import MCD19AODProvider
                self._atmo_provider = MCD19AODProvider(
                    cache_dir=self.config.atmo_prior.cache_dir,
                    source=_earthaccess_source_from_auth(self._auth),
                )
            elif provider_name == "vnp19":
                from siac.priors.atmospheric.mcd19_earthaccess import VNP19AODProvider
                self._atmo_provider = VNP19AODProvider(
                    cache_dir=self.config.atmo_prior.cache_dir,
                    source=_earthaccess_source_from_auth(self._auth),
                )
            else:
                logger.warning(
                    f"Unknown atmospheric provider '{provider_name}', "
                    f"falling back to CAMS"
                )
                self._atmo_provider = CAMSProvider(
                    self.config.atmo_prior.data_path,
                    temporal_interp=self.config.atmo_prior.temporal_interpolation == "linear",
                    download_missing=self.config.atmo_prior.download_missing,
                    auth=self._auth,
                    cache_dir=self.config.atmo_prior.cache_dir,
                )

        # Get bounds and CRS from AOI
        bounds = aoi.get_bounds()
        crs = aoi.crs

        # Try to get native CRS resolution info
        resolution = 10.0  # default
        obs_time = metadata.get("observation_time", datetime.now())

        return self._atmo_provider.get_prior(bounds, crs, obs_time, resolution)

    def _get_surface_prior(self, aoi: AOI, geometry: GeometryAngles, metadata: dict):
        """Get surface prior from BRDF, scoped to AOI."""
        # Get BRDF provider
        brdf_provider = self._get_brdf_provider()

        # Fetch BRDF parameters scoped to AOI
        bounds = aoi.get_bounds()
        obs_time = metadata.get("observation_time", datetime.now())
        resolution = 500.0  # MODIS native resolution for BRDF
        bands = _select_surface_prior_bands(metadata.get("sensor_config"))

        method = getattr(self.config.surface_prior, "method", "kernel_model")
        if method == "monthly_database":
            raise RuntimeError(
                "Route-B monthly-database surface priors require the full ObservationBundle "
                "and are only available through the main pipeline."
            )
        if method == "whittaker":
            brdf_weights = brdf_provider.get_temporal_brdf_parameters(
                bounds=bounds,
                crs=aoi.crs,
                obs_time=obs_time,
                target_resolution=resolution,
                bands=bands,
                temporal_window=self.config.brdf.temporal_window,
            )
            deriver = BRDFWhittakerDeriver(
                temporal_lambda=self.config.surface_prior.whittaker_lambda,
                psf_sigma_x=self.config.surface_prior.psf_sigma_x,
                psf_sigma_y=self.config.surface_prior.psf_sigma_y,
                apply_psf=self.config.surface_prior.apply_psf,
            )
            return deriver.compute_surface_prior(brdf_weights, geometry, obs_time=obs_time)

        brdf_weights = brdf_provider.get_brdf_parameters(
            bounds=bounds,
            crs=aoi.crs,
            obs_time=obs_time,
            target_resolution=resolution,
            bands=bands,
            temporal_window=self.config.brdf.temporal_window,
        )
        deriver = KernelModelDeriver(
            psf_sigma_x=self.config.surface_prior.psf_sigma_x,
            psf_sigma_y=self.config.surface_prior.psf_sigma_y,
            apply_psf=self.config.surface_prior.apply_psf,
        )
        return deriver.compute_surface_prior(brdf_weights, geometry)

    def _get_brdf_provider(self):
        """Get or create the BRDF product provider."""
        if self._brdf_provider is not None:
            return self._brdf_provider

        provider_name = self.config.brdf.provider
        if provider_name == "mcd43":
            from siac.priors.brdf.mcd43_earthaccess import MCD43EarthAccessProvider
            self._brdf_provider = MCD43EarthAccessProvider(
                cache_dir=self.config.brdf.cache_dir,
                source=_earthaccess_source_from_auth(self._auth),
            )
        elif provider_name == "vnp43":
            from siac.priors.brdf.vnp43_earthaccess import VNP43EarthAccessProvider
            self._brdf_provider = VNP43EarthAccessProvider(
                cache_dir=self.config.brdf.cache_dir,
                source=_earthaccess_source_from_auth(self._auth),
            )
        elif provider_name == "mcd19":
            from siac.priors.brdf.mcd43_earthaccess import MCD19EarthAccessProvider
            self._brdf_provider = MCD19EarthAccessProvider(
                cache_dir=self.config.brdf.cache_dir,
                source=_earthaccess_source_from_auth(self._auth),
            )
        elif provider_name == "gee":
            from siac.priors.brdf.gee_stub import GEEBRDFProvider
            self._brdf_provider = GEEBRDFProvider()
        else:
            raise ValueError(
                f"Unknown BRDF provider '{provider_name}'. "
                f"Available: 'mcd43', 'vnp43', 'mcd19', 'gee'"
            )

        return self._brdf_provider

    def _get_rt_model(self, sensor_config: SensorConfig):
        """Get RT model backend."""
        if self._rt_model is not None:
            return self._rt_model

        rt_config = self.config.rt_model
        if rt_config.backend == "emulator":
            try:
                self._rt_model = TwoLayerNNEmulator(
                    emulator_dir=rt_config.emulator_dir,
                    sensor_id=sensor_config.sensor_id,
                    satellite_id=sensor_config.satellite_id,
                )
            except Exception as e:
                logger.warning(f"Emulator not available: {e}, falling back to LUT")
                rt_config.backend = "lut"

        if rt_config.backend == "lut" and rt_config.lut_path:
            self._rt_model = ZarrLUTBackend(
                rt_config.lut_path,
                interpolation_method=rt_config.lut_interpolation,
                storage_options=rt_config.lut_storage_options,
            )

        return self._rt_model

    def _solve_atmosphere(self, toa, surface_prior, geometry, atmo_prior, rt_model, cloud_mask, sensor_config):
        """Solve for atmospheric parameters."""
        if self._solver is None:
            solver_config = MultiGridConfig(
                aot_gamma=self.config.solver.aot_gamma,
                tcwv_gamma=self.config.solver.tcwv_gamma,
                aot_bounds=tuple(self.config.solver.aot_bounds),
                tcwv_bounds=tuple(self.config.solver.tcwv_bounds),
            )
            self._solver = MultiGridSolver(solver_config)

        bands = [sensor_config.get_band(name) for name in sensor_config.band_names[:6]]
        return self._solver.solve(toa, surface_prior, geometry, atmo_prior, rt_model, cloud_mask, bands)

    def _save_output(self, result: CorrectionResult, output_path: Path):
        """Save correction results."""
        from siac.io import write_dataset
        output_path = Path(output_path)
        output_path.mkdir(parents=True, exist_ok=True)
        write_dataset(result.boa, output_path / "boa.nc")
        logger.info(f"Saved output to {output_path}")


def process_sentinel2(input_path: str, output_path: str | None = None, **_kwargs) -> CorrectionResult:
    """Convenience function for Sentinel-2 processing."""
    siac = SIAC.from_defaults(sensor="s2")
    return siac.process(input_path, output_path)


def process_landsat8(input_path: str, output_path: str | None = None, **_kwargs) -> CorrectionResult:
    """Convenience function for Landsat 8 processing."""
    siac = SIAC.from_defaults(sensor="l8")
    return siac.process(input_path, output_path)


def _resolve_s2_backend(config: SIACConfig, *, auth: CredentialManager | None = None):
    """Resolve configured Sentinel-2 data backend object."""
    backend_name = config.s2_data.backend
    if backend_name == "cdse":
        from siac.io.copernicus_dataspace import CopernicusDataspaceBackend
        return CopernicusDataspaceBackend(
            access_key=config.s2_data.cdse_access_key,
            secret_key=config.s2_data.cdse_secret_key,
            auth=auth,
        )
    if backend_name == "gcs":
        from siac.io.gcs_sentinel2 import GCSSentinel2Backend
        return GCSSentinel2Backend()
    if backend_name == "local":
        return None
    raise ValueError(f"Unknown S2 backend: {backend_name!r}")


def _coerce_date(value: date | str | None) -> date | None:
    if value is None:
        return None
    if isinstance(value, datetime):
        return value.date()
    if isinstance(value, date):
        return value
    return datetime.strptime(str(value), "%Y-%m-%d").date()


def _apply_s2_query_defaults(query: S2Query, *, config: SIACConfig) -> S2Query:
    if query.max_cloud_cover == 100.0:
        query.max_cloud_cover = config.s2_data.max_cloud_cover
    if query.processing_level == "L1C" and config.s2_data.processing_level != "L1C":
        query.processing_level = config.s2_data.processing_level
    return query


def _coerce_s2_query(query: S2Query | str | Path, *, config: SIACConfig):
    from siac.io.s2_data_source import S2Query

    if isinstance(query, S2Query):
        q = S2Query(**query.__dict__)
    else:
        raw = str(query).strip()
        try:
            q = S2Query.from_tile_date(raw)
        except ValueError:
            q = S2Query.from_product_id(raw)
            if "MSIL2A" in raw:
                q.processing_level = "L2A"
            elif "MSIL1C" in raw:
                q.processing_level = "L1C"
    return _apply_s2_query_defaults(q, config=config)


def resolve_s2_input(
    query: S2Query | str | Path,
    config: SIACConfig,
    *,
    auth: CredentialManager | None = None,
) -> Path:
    """Resolve S2 query/path to local SAFE directory for M1 preprocessing."""
    local_candidate = Path(query).expanduser() if isinstance(query, Path) else Path(str(query)).expanduser()
    if local_candidate.exists():
        return local_candidate

    auth_obj = auth or CredentialManager.from_config(config)
    backend = _resolve_s2_backend(config, auth=auth_obj)
    if backend is None:
        raise DataNotFoundError(
            "S2 backend is 'local', but input path does not exist. "
            "Provide a local SAFE path or switch config.s2_data.backend to 'cdse' or 'gcs'."
        )

    from siac.io.s2_data_source import S2DataAccess

    cache_dir = config.s2_data.cache_dir
    accessor = S2DataAccess(backend=backend, cache_dir=cache_dir)
    q = _coerce_s2_query(query, config=config)
    return accessor.get(q, dest_dir=cache_dir)


def search_sentinel2(
    *,
    tile: str | None = None,
    date: date | str | None = None,
    date_value: date | str | None = None,
    start_date: date | str | None = None,
    end_date: date | str | None = None,
    bbox: tuple[float, float, float, float] | None = None,
    max_cloud_cover: float = 80.0,
    backend: str = "cdse",
    config: SIACConfig | None = None,
    auth: CredentialManager | None = None,
) -> list[S2Product]:
    """Convenience API to search Sentinel-2 products without downloading."""
    from siac.io.s2_data_source import S2Query, search_s2

    cfg = config or SIACConfig(sensor="s2")
    cfg = cfg.with_overrides(s2_data={"backend": backend, "max_cloud_cover": max_cloud_cover})
    auth_obj = auth or CredentialManager.from_config(cfg)
    backend_obj = _resolve_s2_backend(cfg, auth=auth_obj)
    if backend_obj is None:
        raise ValueError("search_sentinel2 does not support backend='local'.")

    query = S2Query(
        mgrs_tile=tile,
        date=_coerce_date(date if date is not None else date_value),
        start_date=_coerce_date(start_date),
        end_date=_coerce_date(end_date),
        bbox=bbox,
        max_cloud_cover=max_cloud_cover,
        processing_level=cfg.s2_data.processing_level,
    )
    return search_s2(backend_obj, query)


def siac_process_s2(
    config: SIACConfig,
    query: S2Query | str | Path,
    *,
    output_path: str | Path | None = None,
    auth: CredentialManager | None = None,
) -> CorrectionResult:
    """Run full SIAC process from S2 query/path (product ID, tile+date, or local SAFE)."""
    auth_obj = auth or CredentialManager.from_config(config)
    input_path = resolve_s2_input(query, config, auth=auth_obj)
    siac_obj = SIAC(config)
    siac_obj._auth = auth_obj
    return siac_obj.process(input_path, output_path)

def siac_process(
    config: SIACConfig,
    input_path: Path,
    *,
    aoi: AOI | None = None,
    auth: CredentialManager | None = None,
    preprocessor: PreprocessorFn | None = None,
    atmo_provider: AtmoPriorFn | None = None,
    surface_prior_provider: SurfacePriorFn | None = None,
    grid_assembler: GridAssemblerFn | None = None,
    solver: SolverFn | None = None,
    corrector: CorrectorFn | None = None,
    rt_model: Any | None = None,
) -> CorrectionResult:
    """Public entry point for the modular pipeline.

    Resolves ``None`` arguments to config-driven defaults, then delegates
    to :func:`run_pipeline`.
    """
    if auth is None:
        auth = CredentialManager.from_config(config)

    preprocessor = preprocessor or _resolve_preprocessor(config)
    atmo_provider = atmo_provider or _resolve_atmo_provider(config, auth=auth)
    surface_prior_provider = surface_prior_provider or _resolve_surface_prior_provider(config, auth=auth)
    grid_assembler = grid_assembler or _resolve_grid_assembler()
    solver = solver or _resolve_solver(config)
    corrector = corrector or _resolve_corrector(config)
    rt_model = rt_model or _resolve_rt_model_for_pipeline(config, auth=auth)

    return run_pipeline(
        input_path,
        aoi,
        config,
        preprocessor=preprocessor,
        atmo_provider=atmo_provider,
        surface_prior_provider=surface_prior_provider,
        grid_assembler=grid_assembler,
        solver=solver,
        corrector=corrector,
        rt_model=rt_model,
    )


# ── _resolve_* helpers ────────────────────────────────────────────────

def _resolve_preprocessor(config: SIACConfig) -> PreprocessorFn:
    """Return a callable ``(path, aoi) -> ObservationBundle``."""
    sensor = config.sensor
    if sensor in ("s2", "sentinel2"):
        from inspect import signature

        from siac.satellite.sentinel2 import Sentinel2Preprocessor

        cloud_cfg = config.cloud_mask.model_dump(exclude={"user_callable"})
        params = signature(Sentinel2Preprocessor).parameters
        if "config" in params:
            pp = Sentinel2Preprocessor(config={"cloud_mask": cloud_cfg})
        else:
            pp = Sentinel2Preprocessor()
            if hasattr(pp, "config") and isinstance(pp.config, dict):
                pp.config.setdefault("cloud_mask", cloud_cfg)
        return pp.preprocess
    raise ValueError(f"Unknown sensor: {sensor!r}")


def _resolve_atmo_provider(
    config: SIACConfig,
    auth: CredentialManager | None = None,
) -> AtmoPriorFn:
    """Return a callable ``(bounds, crs, time, res) -> AtmosphericState``."""
    provider_name = config.atmo_prior.provider
    if provider_name == "cams":
        provider = CAMSProvider(
            config.atmo_prior.data_path,
            temporal_interp=config.atmo_prior.temporal_interpolation == "linear",
            download_missing=config.atmo_prior.download_missing,
            auth=auth,
            cache_dir=config.atmo_prior.cache_dir,
        )
        return provider.get_prior
    if provider_name == "merra2":
        from siac.priors.atmospheric.merra2 import MERRA2Provider
        provider = MERRA2Provider(
            cache_dir=config.atmo_prior.cache_dir,
            source=_earthaccess_source_from_auth(auth),
        )
        return provider.get_prior
    if provider_name == "mcd19":
        from siac.priors.atmospheric.mcd19_earthaccess import MCD19AODProvider
        provider = MCD19AODProvider(
            cache_dir=config.atmo_prior.cache_dir,
            source=_earthaccess_source_from_auth(auth),
        )
        return provider.get_prior
    if provider_name == "vnp19":
        from siac.priors.atmospheric.mcd19_earthaccess import VNP19AODProvider
        provider = VNP19AODProvider(
            cache_dir=config.atmo_prior.cache_dir,
            source=_earthaccess_source_from_auth(auth),
        )
        return provider.get_prior
    raise ValueError(f"Unknown atmo provider: {provider_name!r}")


def _resolve_surface_prior_provider(
    config: SIACConfig,
    auth: CredentialManager | None = None,
) -> SurfacePriorFn:
    """Return a callable matching the M3 signature -> SurfacePrior."""
    provider_name = getattr(config.brdf, "provider", "mcd43")

    if provider_name == "mcd43":
        from siac.priors.brdf.mcd43_earthaccess import MCD43EarthAccessProvider
        provider_cls = MCD43EarthAccessProvider
    elif provider_name == "vnp43":
        from siac.priors.brdf.vnp43_earthaccess import VNP43EarthAccessProvider
        provider_cls = VNP43EarthAccessProvider
    elif provider_name == "mcd19":
        from siac.priors.brdf.mcd43_earthaccess import MCD19EarthAccessProvider
        provider_cls = MCD19EarthAccessProvider
    else:
        raise ValueError(f"Unknown BRDF provider for surface prior: {provider_name!r}")

    # Create provider and deriver once, then close over them
    brdf_prov = provider_cls(
        cache_dir=config.brdf.cache_dir,
        source=_earthaccess_source_from_auth(auth),
    )
    method = getattr(config.surface_prior, "method", "kernel_model")
    fallback_deriver = KernelModelDeriver(
        psf_sigma_x=config.surface_prior.psf_sigma_x,
        psf_sigma_y=config.surface_prior.psf_sigma_y,
        apply_psf=config.surface_prior.apply_psf,
    )
    if method == "whittaker":
        deriver = BRDFWhittakerDeriver(
            temporal_lambda=config.surface_prior.whittaker_lambda,
            psf_sigma_x=config.surface_prior.psf_sigma_x,
            psf_sigma_y=config.surface_prior.psf_sigma_y,
            apply_psf=config.surface_prior.apply_psf,
        )

        def _brdf_surface_prior(observation, atmo_prior, rt_model, resolution):
            _ = (atmo_prior, rt_model)
            brdf_weights = brdf_prov.get_temporal_brdf_parameters(
                bounds=observation.bounds,
                crs=observation.crs,
                obs_time=observation.metadata["observation_time"],
                target_resolution=resolution,
                bands=_select_surface_prior_bands(observation.sensor_config),
                temporal_window=config.brdf.temporal_window,
            )
            return deriver.compute_surface_prior(
                brdf_weights,
                observation.geometry,
                obs_time=observation.metadata["observation_time"],
            )

        _brdf_surface_prior.requires_atmo_prior = False
        return _brdf_surface_prior

    if method == "monthly_database":
        def _monthly_surface_prior(observation, atmo_prior, rt_model, resolution):
            if atmo_prior is None:
                raise ValueError("Route-B monthly_database surface prior requires an atmospheric prior")

            visible_bands = _select_visible_surface_prior_bands(observation.sensor_config)
            query_bands = _select_route_b_query_bands(observation.sensor_config)
            target_geometry = resample_geometry_for_surface_prior(
                observation,
                resolution=resolution,
            )
            database = build_monthly_surface_prior_database(
                observation=observation,
                brdf_provider=brdf_prov,
                resolution=resolution,
                geometry=target_geometry,
                visible_bands=visible_bands,
                query_bands=query_bands,
            )
            prior = query_surface_prior_from_monthly_database(
                observation=observation,
                atmo_prior=atmo_prior,
                rt_model=rt_model,
                database=database,
                query_band_names=tuple(band.name for band in query_bands),
                visible_band_names=tuple(band.name for band in visible_bands),
                k_neighbors=3,
            )
            if bool(np.asarray(prior.mask.values).any()):
                return prior

            logger.warning(
                "Route-B monthly database produced no valid surface-prior pixels; "
                "falling back to kernel-model BRDF priors."
            )
            brdf_weights = brdf_prov.get_brdf_parameters(
                bounds=observation.bounds,
                crs=observation.crs,
                obs_time=observation.metadata["observation_time"],
                target_resolution=resolution,
                bands=visible_bands,
                temporal_window=config.brdf.temporal_window,
            )
            return fallback_deriver.compute_surface_prior(brdf_weights, target_geometry)

        _monthly_surface_prior.requires_atmo_prior = True
        return _monthly_surface_prior

    def _brdf_surface_prior(observation, atmo_prior, rt_model, resolution):
        _ = (atmo_prior, rt_model)
        brdf_weights = brdf_prov.get_brdf_parameters(
            bounds=observation.bounds,
            crs=observation.crs,
            obs_time=observation.metadata["observation_time"],
            target_resolution=resolution,
            bands=_select_surface_prior_bands(observation.sensor_config),
            temporal_window=config.brdf.temporal_window,
        )
        return fallback_deriver.compute_surface_prior(brdf_weights, observation.geometry)

    _brdf_surface_prior.requires_atmo_prior = False
    return _brdf_surface_prior


def _resolve_grid_assembler() -> GridAssemblerFn:
    """Return the default grid assembler function."""
    from siac.grid.assembler import assemble_grids
    return assemble_grids


def _resolve_solver(config: SIACConfig) -> SolverFn:
    """Return the default solver callable."""
    def _default_solver(inputs: SolverInputBundle, _cfg) -> SolvedAtmosphere:
        solver_config = MultiGridConfig(
            aot_gamma=config.solver.aot_gamma,
            tcwv_gamma=config.solver.tcwv_gamma,
            aot_bounds=tuple(config.solver.aot_bounds),
            tcwv_bounds=tuple(config.solver.tcwv_bounds),
        )
        mg_solver = MultiGridSolver(solver_config)
        result = mg_solver.solve(
            inputs.toa,
            inputs.surface_prior,
            inputs.geometry,
            inputs.atmo_prior,
            inputs.rt_model,
            inputs.cloud_mask,
            inputs.bands,
        )
        solved_atmo = inputs.atmo_prior.with_updated_aot_tcwv(
            aot=result.aot,
            tcwv=result.tcwv,
            aot_unc=result.aot_unc,
            tcwv_unc=result.tcwv_unc,
        )
        return SolvedAtmosphere(
            atmo_state=solved_atmo,
            aot=result.aot,
            tcwv=result.tcwv,
            aot_unc=result.aot_unc,
            tcwv_unc=result.tcwv_unc,
            cost_final=float(result.final_cost),
            n_iterations=result.n_iterations,
            converged=result.success,
        )

    return _default_solver


def _resolve_corrector(_config: SIACConfig) -> CorrectorFn:
    """Return the default corrector callable."""
    def _default_corrector(obs: ObservationBundle, solved: SolvedAtmosphere, rt_model) -> CorrectionResult:
        corrector_obj = AtmosphericCorrector(rt_model, obs.sensor_config)
        first_band = obs.toa[next(iter(obs.toa.data_vars))]
        atmo = solved.atmo_state
        matched_atmo = AtmosphericState(
            aot=_resample_field_to_template(atmo.aot, first_band),
            tcwv=_resample_field_to_template(atmo.tcwv, first_band),
            tco3=_resample_field_to_template(atmo.tco3, first_band),
            aot_unc=_resample_field_to_template(atmo.aot_unc, first_band),
            tcwv_unc=_resample_field_to_template(atmo.tcwv_unc, first_band),
            tco3_unc=_resample_field_to_template(atmo.tco3_unc, first_band),
            elevation=_resample_field_to_template(atmo.elevation, first_band),
        )
        return corrector_obj.correct(obs.toa, obs.geometry, matched_atmo, obs.cloud_mask)

    return _default_corrector


def _resample_field_to_template(field, template):
    """Resample a 2-D field to match a template grid for M6 correction."""
    if field.shape == template.shape:
        return field

    if (
        all(dim in field.dims for dim in ("y", "x"))
        and all(dim in template.dims for dim in ("y", "x"))
        and all(coord in field.coords for coord in ("y", "x"))
        and all(coord in template.coords for coord in ("y", "x"))
    ):
        try:
            return field.interp(y=template.coords["y"], x=template.coords["x"], method="linear")
        except Exception:
            pass

    src = np.asarray(field.values, dtype=np.float32)
    if src.ndim != 2:
        return field

    from scipy import ndimage
    h_out = int(template.sizes["y"])
    w_out = int(template.sizes["x"])
    if src.shape[0] == 0 or src.shape[1] == 0:
        out = np.full((h_out, w_out), np.nan, dtype=np.float32)
    else:
        out = ndimage.zoom(src, (h_out / src.shape[0], w_out / src.shape[1]), order=1)
        out = out[:h_out, :w_out]
        if out.shape != (h_out, w_out):
            padded = np.full((h_out, w_out), np.nan, dtype=np.float32)
            padded[: out.shape[0], : out.shape[1]] = out
            out = padded

    return field.__class__(
        out,
        dims=template.dims,
        coords={d: template.coords[d] for d in template.dims if d in template.coords},
    )


# Default sensor mapping for config.sensor -> (sensor_id, satellite_id)
_SENSOR_DEFAULTS: dict[str, tuple[str, str]] = {
    "s2": ("MSI", "S2A"),
    "s2a": ("MSI", "S2A"),
    "s2b": ("MSI", "S2B"),
    "s2c": ("MSI", "S2C"),
    "sentinel2": ("MSI", "S2A"),
    "l8": ("OLI", "L8"),
    "l9": ("OLI", "L9"),
    "auto": ("MSI", "S2A"),  # fallback
}


def _resolve_rt_model_for_pipeline(
    config: SIACConfig,
    auth: CredentialManager | None = None,
    *,
    sensor_config: SensorConfig | None = None,
):
    """Return an RT model instance from config.

    Parameters
    ----------
    sensor_config : SensorConfig, optional
        If provided, ``sensor_id`` and ``satellite_id`` are read from
        the config.  Otherwise they are inferred from ``config.sensor``.
    """
    rt_config = config.rt_model
    if sensor_config is not None:
        sid, satid = sensor_config.sensor_id, sensor_config.satellite_id
    else:
        sid, satid = _SENSOR_DEFAULTS.get(getattr(config, "sensor", None), ("MSI", "S2A"))
    if rt_config.backend == "emulator":
        return TwoLayerNNEmulator(
            emulator_dir=rt_config.emulator_dir,
            sensor_id=sid,
            satellite_id=satid,
        )
    if rt_config.backend == "lut" and rt_config.lut_path:
        storage_options = dict(rt_config.lut_storage_options)
        # Inject AWS credentials if not already present
        if (
            auth is not None
            and auth.aws().has_credentials()
            and "key" not in storage_options
            and str(rt_config.lut_path).startswith("s3://")
        ):
            storage_options.update(auth.aws().storage_options())
        return ZarrLUTBackend(
            rt_config.lut_path,
            interpolation_method=rt_config.lut_interpolation,
            storage_options=storage_options,
        )
    raise ValueError(f"Cannot resolve RT model from config: backend={rt_config.backend!r}")
