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

        # 3. Get atmospheric priors
        logger.info("Fetching atmospheric priors...")
        atmo_prior = self._get_atmospheric_prior(toa, metadata)

        # 4. Get surface prior from BRDF
        logger.info("Computing surface prior...")
        surface_prior = self._get_surface_prior(geometry, metadata)

        # 5. Setup RT model
        rt_model = self._get_rt_model(sensor_config)

        # 6. Solve for atmospheric parameters
        logger.info("Solving for atmospheric parameters...")
        solver_result = self._solve_atmosphere(
            toa, surface_prior, geometry, atmo_prior, rt_model, cloud_mask, sensor_config
        )

        # 7. Apply atmospheric correction
        logger.info("Applying atmospheric correction...")
        solved_atmo = atmo_prior.with_updated_aot_tcwv(
            aot=solver_result.aot,
            tcwv=solver_result.tcwv,
            aot_unc=solver_result.aot_unc,
            tcwv_unc=solver_result.tcwv_unc,
        )

        corrector = AtmosphericCorrector(rt_model, sensor_config)
        result = corrector.correct(toa, geometry, solved_atmo, cloud_mask)

        # 8. Save output if path provided
        if output_path is not None:
            self._save_output(result, output_path)

        logger.info(f"Complete. Mean AOT: {float(result.aot.mean()):.3f}")
        return result

    def _get_atmospheric_prior(self, toa: xr.Dataset, metadata: dict) -> AtmosphericState:
        """Get atmospheric prior from configured provider."""
        if self._atmo_provider is None:
            cams_dir = self.config.atmo_prior.data_path
            self._atmo_provider = CAMSProvider(cams_dir)

        # Get bounds from TOA
        first_var = list(toa.data_vars)[0]
        da = toa[first_var]
        if hasattr(da, 'rio'):
            bounds = da.rio.bounds()
            crs = str(da.rio.crs)
            resolution = abs(da.rio.resolution()[0])
        else:
            bounds = (0, 0, 1, 1)
            crs = "EPSG:4326"
            resolution = 10.0

        obs_time = metadata.get("observation_time", datetime.now())
        return self._atmo_provider.get_prior(bounds, crs, obs_time, resolution)

    def _get_surface_prior(self, geometry: GeometryAngles, metadata: dict):
        """Get surface prior from BRDF."""
        deriver = KernelModelDeriver(
            psf_sigma_x=self.config.surface_prior.psf_sigma_x,
            psf_sigma_y=self.config.surface_prior.psf_sigma_y,
            apply_psf=True,
        )
        # TODO: Load BRDF weights from MCD43 provider
        # For now return None - would need brdf_weights
        return None

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
            self._rt_model = ZarrLUTBackend(rt_config.lut_path)

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
