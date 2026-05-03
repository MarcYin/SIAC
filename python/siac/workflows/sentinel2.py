"""Sentinel-2 processing workflow."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, cast

from siac.adapters.auth import CredentialManager
from siac.app.planning import resolve_run_config
from siac.app.requests import (
    SceneProcessRequest,
    Sentinel2ProcessRequest,
    Sentinel2ResolveRequest,
)
from siac.app.sentinel2 import resolve_s2_input as resolve_s2_scene_input
from siac.workflows.scene import process_scene

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path

    from siac.runtime import CorrectionResult


def process_s2(
    request: Sentinel2ProcessRequest,
    *,
    resolve_s2_input_fn: Callable[[Sentinel2ResolveRequest], Path] | None = None,
) -> CorrectionResult:
    resolved_config = resolve_run_config(
        request.config,
        sensor=request.config.sensor if request.config.sensor != "auto" else "s2",
        aoi=request.aoi if request.aoi is not None else request.config.aoi,
        s2_query=cast("Any", request.query),
    )
    auth_obj = request.auth or CredentialManager.from_config(resolved_config)
    if resolve_s2_input_fn is None:
        resolve_s2_input_fn = resolve_s2_scene_input
    input_path = resolve_s2_input_fn(
        Sentinel2ResolveRequest(
            config=request.config,
            query=request.query,
            auth=auth_obj,
        )
    )
    return process_scene(
        request=SceneProcessRequest(
            config=request.config,
            input_path=input_path,
            output_path=request.output_path,
            aoi=request.aoi,
            auth=auth_obj,
        )
    )
