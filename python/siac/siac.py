"""
Main SIAC class for atmospheric correction.

Orchestrates the full pipeline: preprocessing, prior retrieval,
aerosol solving, and atmospheric correction.
"""

from __future__ import annotations
import logging
from datetime import datetime
from pathlib import Path
from typing import Any

import xarray as xr

from siac.core.aoi import AOI
from siac.core.config import SIACConfig
from siac.core.types import AtmosphericState, GeometryAngles, SensorConfig
from siac.satellite import get_preprocessor, detect_sensor
from siac.priors.atmospheric import CAMSProvider
from siac.priors.surface.kernel_model import KernelModelDeriver
from siac.rt.emulator import TwoLayerNNEmulator, EmulatorRegistry
from siac.rt.lut import ZarrLUTBackend
from siac.solver import MultiGridSolver, MultiGridConfig
from siac.correction import AtmosphericCorrector, CorrectionResult

logger = logging.getLogger(__name__)


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

        Args:
            input_path: Path to satellite data
            output_path: Optional output path (default: input_path/BOA/)

        Returns:
            CorrectionResult with BOA reflectance
        """
        input_path = Path(input_path)
        logger.info(f"Processing: {input_path}")

        # 1. Detect sensor and get preprocessor
        sensor_id = self.config.sensor if self.config.sensor != "auto" else detect_sensor(input_path)
        preprocessor = get_preprocessor(sensor_id)
        sensor_config = preprocessor.sensor_config

        # 2. Preprocess satellite data
        logger.info("Preprocessing satellite data...")
        preprocess_result = preprocessor.preprocess(input_path)
        toa = preprocess_result["toa"]
        geometry = preprocess_result["geometry"]
        cloud_mask = preprocess_result["cloud_mask"]
        metadata = preprocess_result["metadata"]

        # 3. Resolve AOI from config or TOA extent
        aoi = self._resolve_aoi(toa)
        logger.info(f"AOI resolved: bounds={aoi.get_bounds()}, crs={aoi.crs}")

        # 4. Get atmospheric priors (now AOI-scoped)
        logger.info("Fetching atmospheric priors...")
        atmo_prior = self._get_atmospheric_prior(aoi, metadata)

        # 5. Get surface prior from BRDF (now AOI-scoped)
        logger.info("Computing surface prior...")
        surface_prior = self._get_surface_prior(aoi, geometry, metadata)

        # 6. Setup RT model
        rt_model = self._get_rt_model(sensor_config)

        # 7. Solve for atmospheric parameters
        logger.info("Solving for atmospheric parameters...")
        solver_result = self._solve_atmosphere(
            toa, surface_prior, geometry, atmo_prior, rt_model, cloud_mask, sensor_config
        )

        # 8. Apply atmospheric correction
        logger.info("Applying atmospheric correction...")
        solved_atmo = atmo_prior.with_updated_aot_tcwv(
            aot=solver_result.aot,
            tcwv=solver_result.tcwv,
            aot_unc=solver_result.aot_unc,
            tcwv_unc=solver_result.tcwv_unc,
        )

        corrector = AtmosphericCorrector(rt_model, sensor_config)
        result = corrector.correct(toa, geometry, solved_atmo, cloud_mask)

        # 9. Save output if path provided
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
                self._atmo_provider = CAMSProvider(self.config.atmo_prior.data_path)
            elif provider_name == "merra2":
                from siac.priors.atmospheric.merra2 import MERRA2Provider
                self._atmo_provider = MERRA2Provider(
                    cache_dir=self.config.atmo_prior.cache_dir,
                )
            else:
                logger.warning(
                    f"Unknown atmospheric provider '{provider_name}', "
                    f"falling back to CAMS"
                )
                self._atmo_provider = CAMSProvider(self.config.atmo_prior.data_path)

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

        # Default MODIS bands for atmospheric correction (1-7)
        bands = list(range(1, 8))

        brdf_weights = brdf_provider.get_brdf_parameters(
            bounds=bounds,
            crs=aoi.crs,
            obs_time=obs_time,
            target_resolution=resolution,
            bands=bands,
            temporal_window=self.config.brdf.temporal_window,
        )

        # Derive surface prior via kernel model
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
            )
        elif provider_name == "gee":
            from siac.priors.brdf.gee_stub import GEEBRDFProvider
            self._brdf_provider = GEEBRDFProvider()
        else:
            raise ValueError(
                f"Unknown BRDF provider '{provider_name}'. "
                f"Available: 'mcd43', 'gee'"
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
        from siac.io import write_dataset, write_cog
        output_path = Path(output_path)
        output_path.mkdir(parents=True, exist_ok=True)
        write_dataset(result.boa, output_path / "boa.nc")
        logger.info(f"Saved output to {output_path}")


def process_sentinel2(input_path: str, output_path: str | None = None, **kwargs) -> CorrectionResult:
    """Convenience function for Sentinel-2 processing."""
    siac = SIAC.from_defaults(sensor="s2")
    return siac.process(input_path, output_path)


def process_landsat8(input_path: str, output_path: str | None = None, **kwargs) -> CorrectionResult:
    """Convenience function for Landsat 8 processing."""
    siac = SIAC.from_defaults(sensor="l8")
    return siac.process(input_path, output_path)


# =====================================================================
# New modular pipeline entry point (PLANS.md §7)
# =====================================================================

from siac.core.types import (
    ObservationBundle,
    SolvedAtmosphere,
    SolverInputBundle,
    SurfacePrior,
)
from siac.pipeline import (
    AtmoPriorFn,
    CorrectorFn,
    GridAssemblerFn,
    PreprocessorFn,
    SolverFn,
    SurfacePriorFn,
    run_pipeline,
)


def siac_process(
    config: SIACConfig,
    input_path: Path,
    *,
    aoi: AOI | None = None,
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
    preprocessor = preprocessor or _resolve_preprocessor(config)
    atmo_provider = atmo_provider or _resolve_atmo_provider(config)
    surface_prior_provider = surface_prior_provider or _resolve_surface_prior_provider(config)
    grid_assembler = grid_assembler or _resolve_grid_assembler()
    solver = solver or _resolve_solver(config)
    corrector = corrector or _resolve_corrector(config)
    rt_model = rt_model or _resolve_rt_model_for_pipeline(config)

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
        from siac.satellite.sentinel2 import Sentinel2Preprocessor
        pp = Sentinel2Preprocessor()
        return pp.preprocess
    raise ValueError(f"Unknown sensor: {sensor!r}")


def _resolve_atmo_provider(config: SIACConfig) -> AtmoPriorFn:
    """Return a callable ``(bounds, crs, time, res) -> AtmosphericState``."""
    provider_name = config.atmo_prior.provider
    if provider_name == "cams":
        provider = CAMSProvider(config.atmo_prior.data_path)
        return provider.get_prior
    if provider_name == "merra2":
        from siac.priors.atmospheric.merra2 import MERRA2Provider
        provider = MERRA2Provider(cache_dir=config.atmo_prior.cache_dir)
        return provider.get_prior
    raise ValueError(f"Unknown atmo provider: {provider_name!r}")


def _resolve_surface_prior_provider(config: SIACConfig) -> SurfacePriorFn:
    """Return a callable matching the M3 signature -> SurfacePrior."""
    # Default: BRDF-derived prior via kernel model
    def _brdf_surface_prior(bounds, crs, obs_time, sensor_config, geometry, resolution):
        from siac.priors.brdf.mcd43_earthaccess import MCD43EarthAccessProvider
        brdf_prov = MCD43EarthAccessProvider(
            cache_dir=config.brdf.cache_dir,
        )
        brdf_weights = brdf_prov.get_brdf_parameters(
            bounds=bounds,
            crs=crs,
            obs_time=obs_time,
            target_resolution=resolution,
            bands=list(range(1, 8)),
            temporal_window=config.brdf.temporal_window,
        )
        deriver = KernelModelDeriver(
            psf_sigma_x=config.surface_prior.psf_sigma_x,
            psf_sigma_y=config.surface_prior.psf_sigma_y,
            apply_psf=config.surface_prior.apply_psf,
        )
        return deriver.compute_surface_prior(brdf_weights, geometry)

    return _brdf_surface_prior


def _resolve_grid_assembler() -> GridAssemblerFn:
    """Return the default grid assembler function."""
    from siac.grid.assembler import assemble_grids
    return assemble_grids


def _resolve_solver(config: SIACConfig) -> SolverFn:
    """Return the default solver callable."""
    def _default_solver(inputs: SolverInputBundle, cfg) -> SolvedAtmosphere:
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


def _resolve_corrector(config: SIACConfig) -> CorrectorFn:
    """Return the default corrector callable."""
    def _default_corrector(obs: ObservationBundle, solved: SolvedAtmosphere, rt_model) -> CorrectionResult:
        corrector_obj = AtmosphericCorrector(rt_model, obs.sensor_config)
        return corrector_obj.correct(obs.toa, obs.geometry, solved.atmo_state, obs.cloud_mask)

    return _default_corrector


def _resolve_rt_model_for_pipeline(config: SIACConfig):
    """Return an RT model instance from config."""
    rt_config = config.rt_model
    if rt_config.backend == "emulator":
        return TwoLayerNNEmulator(
            emulator_dir=rt_config.emulator_dir,
            sensor_id="MSI",
            satellite_id="S2A",
        )
    if rt_config.backend == "lut" and rt_config.lut_path:
        return ZarrLUTBackend(
            rt_config.lut_path,
            interpolation_method=rt_config.lut_interpolation,
            storage_options=rt_config.lut_storage_options,
        )
    raise ValueError(f"Cannot resolve RT model from config: backend={rt_config.backend!r}")
