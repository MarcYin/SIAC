"""Public SIAC API."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

from siac.adapters.auth import CredentialManager
from siac.api.requests import (
    AOISpec,
    PathLike,
    SceneProcessRequest,
    Sentinel2ProcessRequest,
    Sentinel2QuerySpec,
    Sentinel2ResolveRequest,
    Sentinel2SearchRequest,
)
from siac.app.sentinel2 import (
    apply_s2_query_defaults as apply_s2_query_defaults,
)
from siac.app.sentinel2 import (
    coerce_date as coerce_date,
)
from siac.app.sentinel2 import (
    coerce_s2_query as coerce_s2_query,
)
from siac.app.sentinel2 import (
    resolve_s2_input as app_resolve_s2_input,
)
from siac.app.sentinel2 import (
    search_sentinel2 as app_search_sentinel2,
)
from siac.config import SensorName, SIACConfig
from siac.public_models import (  # noqa: TC001
    PreparedMonthlyCompositeBuildResult as PreparedMonthlyCompositeBuildResult,
)
from siac.workflows.monthly_composites import (
    prepare_monthly_composites as workflow_prepare_monthly_composites,
)
from siac.workflows.scene import process_scene as workflow_process_scene
from siac.workflows.sentinel2 import process_s2 as workflow_process_s2

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from datetime import datetime
    from pathlib import Path

    from siac.adapters.data.s2_data_source import S2Product
    from siac.runtime import CorrectionResult


class SIAC:
    """User-facing SIAC session API."""

    def __init__(self, config: SIACConfig):
        self.config = config
        self._auth = CredentialManager.from_config(config)

    @classmethod
    def from_config(cls, config_path: str | Path) -> SIAC:
        return cls(SIACConfig.from_file(config_path))

    @classmethod
    def from_defaults(cls, sensor: SensorName = SensorName.AUTO) -> SIAC:
        return cls(SIACConfig(sensor=sensor))

    def process(
        self,
        input_path: str | Path,
        output_path: str | Path | None = None,
        *,
        aoi: AOISpec = None,
    ) -> CorrectionResult:
        request = SceneProcessRequest(
            config=self.config,
            input_path=input_path,
            output_path=output_path,
            aoi=aoi,
            auth=self._auth,
        )
        result = workflow_process_scene(request)
        logger.info("Complete. Mean AOT: %.3f", float(result.aot.mean()))
        return result


def resolve_s2_input(
    query: Sentinel2QuerySpec,
    config: SIACConfig,
    *,
    auth: CredentialManager | None = None,
) -> Path:
    request = Sentinel2ResolveRequest(config=config, query=query, auth=auth)
    return app_resolve_s2_input(request)


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
    return app_search_sentinel2(request)


def process_sentinel2(
    input_path: PathLike,
    output_path: PathLike | None = None,
    **kwargs: Any,
) -> CorrectionResult:
    return SIAC.from_defaults(sensor=SensorName.S2).process(input_path, output_path, **kwargs)


def siac_process_s2(
    config: SIACConfig,
    query: Sentinel2QuerySpec,
    *,
    output_path: str | Path | None = None,
    aoi: AOISpec = None,
    auth: CredentialManager | None = None,
) -> CorrectionResult:
    request = Sentinel2ProcessRequest(
        config=config,
        query=query,
        output_path=output_path,
        aoi=aoi,
        auth=auth,
    )
    return workflow_process_s2(request)


def siac_process(
    config: SIACConfig,
    input_path: PathLike,
    *,
    aoi: AOISpec = None,
    auth: CredentialManager | None = None,
) -> CorrectionResult:
    request = SceneProcessRequest(
        config=config,
        input_path=input_path,
        aoi=aoi,
        auth=auth,
    )
    return workflow_process_scene(request)


def prepare_monthly_composites(
    config: SIACConfig,
    *,
    aoi: AOISpec,
    year_months: list[tuple[int, int]] | tuple[tuple[int, int], ...],
    resolution: float | None,
    output_path: str | Path,
    auth: CredentialManager | None = None,
) -> PreparedMonthlyCompositeBuildResult:
    return workflow_prepare_monthly_composites(
        config,
        aoi=aoi,
        year_months=year_months,
        resolution=resolution,
        output_path=output_path,
        auth=auth,
    )
