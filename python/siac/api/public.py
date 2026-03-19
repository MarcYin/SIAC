"""Public SIAC API."""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any

from siac.adapters.auth import CredentialManager
from siac.api.requests import (
    SceneProcessRequest,
    Sentinel2ProcessRequest,
    Sentinel2ResolveRequest,
    Sentinel2SearchRequest,
)
from siac.app.assembly import (
    resolve_atmo_provider,
    resolve_brdf_provider,
    resolve_rt_model_for_pipeline,
    resolve_solver,
    resolve_surface_prior_provider,
)
from siac.app.planning import coerce_aoi_spec, resolve_run_config
from siac.config import SIACConfig
from siac.domain.aoi import AOI
from siac.workflows.scene import process_scene, save_output
from siac.workflows.sentinel2 import (
    apply_s2_query_defaults,
    coerce_date,
    coerce_s2_query,
    process_s2,
)
from siac.workflows.sentinel2 import (
    resolve_s2_input as resolve_s2_input_workflow,
)
from siac.workflows.sentinel2 import (
    search_sentinel2 as search_sentinel2_workflow,
)

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    import xarray as xr

    from siac.adapters.data.s2_data_source import S2Product, S2Query
    from siac.algorithms.correction import CorrectionResult
    from siac.domain import AtmosphericState, GeometryAngles, SensorConfig


class SIAC:
    """User-facing SIAC session facade."""

    def __init__(self, config: SIACConfig):
        self.config = config
        self._auth = CredentialManager.from_config(config)
        self._brdf_provider = None
        self._rt_model = None
        self._solver = None

    @classmethod
    def from_config(cls, config_path: str | Path) -> SIAC:
        return cls(SIACConfig.from_file(config_path))

    @classmethod
    def from_defaults(cls, sensor: str = "auto") -> SIAC:
        return cls(SIACConfig(sensor=sensor))

    def process(
        self,
        input_path: str | Path,
        output_path: str | Path | None = None,
        *,
        aoi: AOI | Path | str | tuple[float, float, float, float] | list[float] | None = None,
    ) -> CorrectionResult:
        request = SceneProcessRequest(
            config=self.config,
            input_path=input_path,
            output_path=output_path,
            aoi=aoi,
            auth=self._auth,
        )
        result = process_scene(
            request.config,
            Path(request.input_path),
            output_path=None,
            aoi=request.aoi,
            auth=request.auth,
            aoi_resolver=self._resolve_aoi,
        )
        if output_path is not None:
            self._save_output(result, output_path)
        logger.info("Complete. Mean AOT: %.3f", float(result.aot.mean()))
        return result

    def _resolve_aoi(self, toa: xr.Dataset) -> AOI:
        if self.config.aoi is not None:
            return coerce_aoi_spec(self.config.aoi)
        first_var = list(toa.data_vars)[0]
        return AOI.from_raster(toa[first_var])

    def _get_atmospheric_prior(self, aoi: AOI, metadata: dict) -> AtmosphericState:
        provider_name = self.config.atmo_prior.provider
        config = self.config
        if provider_name not in {"cams", "merra2", "mcd19", "vnp19"}:
            logger.warning(
                "Unknown atmospheric provider '%s', falling back to CAMS",
                provider_name,
            )
            config = self.config.with_overrides(providers={"atmo": {"kind": "cams"}})

        provider_fn = resolve_atmo_provider(config, auth=self._auth)
        bounds = aoi.get_bounds()
        crs = aoi.crs
        resolution = 10.0
        obs_time = metadata.get("observation_time", datetime.now())
        return provider_fn(bounds, crs, obs_time, resolution)

    def _get_surface_prior(self, aoi: AOI, geometry: GeometryAngles, metadata: dict):
        method = getattr(self.config.surface_prior, "method", "kernel_model")
        if method == "monthly_database":
            raise RuntimeError(
                "Route-B monthly-database surface priors require the full ObservationBundle "
                "and are only available through the main pipeline."
            )

        _ = self._get_brdf_provider()
        obs = SimpleNamespace(
            bounds=aoi.get_bounds(),
            crs=aoi.crs,
            metadata=metadata,
            geometry=geometry,
            sensor_config=metadata.get("sensor_config"),
        )
        surface_fn = resolve_surface_prior_provider(self.config, auth=self._auth)
        return surface_fn(obs, None, "rt", 500.0)

    def _get_brdf_provider(self):
        if self._brdf_provider is None:
            self._brdf_provider = resolve_brdf_provider(self.config, auth=self._auth)
        return self._brdf_provider

    def _get_rt_model(self, sensor_config: SensorConfig):
        if self._rt_model is None:
            self._rt_model = resolve_rt_model_for_pipeline(
                self.config,
                auth=self._auth,
                sensor_config=sensor_config,
            )
        return self._rt_model

    def _solve_atmosphere(self, toa, surface_prior, geometry, atmo_prior, rt_model, cloud_mask, sensor_config):
        if self._solver is None:
            self._solver = resolve_solver(self.config)

        bands = [sensor_config.get_band(name) for name in sensor_config.band_names[:6]]
        inputs = SimpleNamespace(
            toa=toa,
            surface_prior=surface_prior,
            geometry=geometry,
            atmo_prior=atmo_prior,
            rt_model=rt_model,
            cloud_mask=cloud_mask,
            bands=bands,
        )
        return self._solver(inputs, self.config)

    def _save_output(self, result: CorrectionResult, output_path: Path | str) -> None:
        save_output(result, output_path)


def resolve_s2_input(
    query: S2Query | str | Path,
    config: SIACConfig,
    *,
    auth: CredentialManager | None = None,
) -> Path:
    request = Sentinel2ResolveRequest(config=config, query=query, auth=auth)
    return resolve_s2_input_workflow(
        request.query,
        request.config,
        auth=request.auth,
    )


def search_sentinel2(
    *,
    tile: str | None = None,
    date: datetime | str | None = None,
    date_value: datetime | str | None = None,
    start_date: datetime | str | None = None,
    end_date: datetime | str | None = None,
    bbox: tuple[float, float, float, float] | None = None,
    max_cloud_cover: float = 80.0,
    backend: str = "cdse",
    config: SIACConfig | None = None,
    auth: CredentialManager | None = None,
) -> list[S2Product]:
    request = Sentinel2SearchRequest(
        tile=tile,
        date=date,
        date_value=date_value,
        start_date=start_date,
        end_date=end_date,
        bbox=bbox,
        max_cloud_cover=max_cloud_cover,
        backend=backend,
        config=config,
        auth=auth,
    )
    return search_sentinel2_workflow(
        tile=request.tile,
        date=request.date,
        date_value=request.date_value,
        start_date=request.start_date,
        end_date=request.end_date,
        bbox=request.bbox,
        max_cloud_cover=request.max_cloud_cover,
        backend=request.backend,
        config=request.config,
        auth=request.auth,
    )


def process_sentinel2(input_path: str, output_path: str | None = None, **kwargs) -> CorrectionResult:
    return SIAC.from_defaults(sensor="s2").process(input_path, output_path, **kwargs)


def process_landsat8(input_path: str, output_path: str | None = None, **kwargs) -> CorrectionResult:
    return SIAC.from_defaults(sensor="l8").process(input_path, output_path, **kwargs)


def siac_process_s2(
    config: SIACConfig,
    query: S2Query | str | Path,
    *,
    output_path: str | Path | None = None,
    aoi: AOI | Path | str | tuple[float, float, float, float] | list[float] | None = None,
    auth: CredentialManager | None = None,
) -> CorrectionResult:
    request = Sentinel2ProcessRequest(
        config=config,
        query=query,
        output_path=output_path,
        aoi=aoi,
        auth=auth,
    )
    return process_s2(
        request.config,
        request.query,
        output_path=request.output_path,
        aoi=request.aoi,
        auth=request.auth,
    )


def siac_process(
    config: SIACConfig,
    input_path: Path,
    *,
    aoi: AOI | Path | str | tuple[float, float, float, float] | list[float] | None = None,
    auth: CredentialManager | None = None,
    preprocessor=None,
    atmo_provider=None,
    surface_prior_provider=None,
    grid_assembler=None,
    solver=None,
    corrector=None,
    rt_model: Any | None = None,
) -> CorrectionResult:
    request = SceneProcessRequest(
        config=config,
        input_path=input_path,
        aoi=aoi,
        auth=auth,
    )
    return process_scene(
        request.config,
        Path(request.input_path),
        aoi=request.aoi,
        auth=request.auth,
        preprocessor=preprocessor,
        atmo_provider=atmo_provider,
        surface_prior_provider=surface_prior_provider,
        grid_assembler=grid_assembler,
        solver=solver,
        corrector=corrector,
        rt_model=rt_model,
    )


__all__ = [
    "SIAC",
    "apply_s2_query_defaults",
    "coerce_aoi_spec",
    "coerce_date",
    "coerce_s2_query",
    "process_landsat8",
    "process_sentinel2",
    "resolve_run_config",
    "resolve_s2_input",
    "search_sentinel2",
    "siac_process",
    "siac_process_s2",
]
