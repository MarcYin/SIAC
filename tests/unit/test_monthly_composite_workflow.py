from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import xarray as xr

import siac.workflows.monthly_composites as workflow_mod
from siac.algorithms.surface.brdf_monthly_composite import (
    MonthlyCompositeCollection,
    MonthlyKernelWeightComposite,
)
from siac.algorithms.surface.monthly_composite_store import MonthlyCompositeStoreGridSpec
from siac.domain import SensorBand
from siac.domain.aoi import AOI
from siac.runtime import BRDFKernelWeights


def test_prepare_monthly_composites_orchestrates_provider_build_and_store(
    monkeypatch,
    tmp_path: Path,
) -> None:
    aoi = AOI.from_bounds((1.0, 2.0, 3.0, 4.0), crs="EPSG:4326")
    coords = {"band": ["B02"], "y": [3.5, 2.5], "x": [1.5, 2.5]}
    cube = xr.DataArray(np.full((1, 2, 2), 0.2, dtype=np.float32), dims=["band", "y", "x"], coords=coords)
    collection = MonthlyCompositeCollection(
        composites=(
            MonthlyKernelWeightComposite(
                kernels=BRDFKernelWeights(
                    f0=cube,
                    f1=xr.zeros_like(cube),
                    f2=xr.zeros_like(cube),
                    f0_unc=xr.full_like(cube, 0.01),
                    f1_unc=xr.full_like(cube, 0.01),
                    f2_unc=xr.full_like(cube, 0.01),
                ),
                quality=xr.DataArray(
                    np.full((2, 2), 0.03, dtype=np.float32),
                    dims=["y", "x"],
                    coords={"y": [3.5, 2.5], "x": [1.5, 2.5]},
                ),
                sample_index=xr.DataArray(
                    np.zeros((2, 2), dtype=np.int16),
                    dims=["y", "x"],
                    coords={"y": [3.5, 2.5], "x": [1.5, 2.5]},
                ),
                year=2023,
                month=7,
            ),
        ),
        source_bands=(SensorBand("B02", 490.0, 65.0, 10.0, 0),),
        source_name="prepared-test",
    )
    captured: dict[str, object] = {}

    def _fake_resolve_run_config(config_obj, **kwargs):  # noqa: ANN001
        captured["resolved"] = kwargs
        return config_obj

    def _fake_resolve_brdf_provider(config_obj, auth=None):  # noqa: ANN001
        captured["provider"] = (config_obj, auth)
        return "provider"

    def _fake_build_monthly_composites_from_brdf(**kwargs):  # noqa: ANN003
        captured["build_kwargs"] = kwargs
        return collection

    monkeypatch.setattr(workflow_mod, "resolve_run_config", _fake_resolve_run_config)
    monkeypatch.setattr(workflow_mod, "resolve_brdf_provider", _fake_resolve_brdf_provider)
    monkeypatch.setattr(workflow_mod, "build_monthly_composites_from_brdf", _fake_build_monthly_composites_from_brdf)

    def _fake_write(collection_obj, output_path, *, grid):  # noqa: ANN001
        captured["write"] = (collection_obj, output_path, grid)
        return Path(output_path)

    monkeypatch.setattr(workflow_mod, "write_monthly_composite_collection", _fake_write)

    auth = object()
    config = SimpleNamespace(sensor="auto")
    result = workflow_mod.prepare_monthly_composites(
        config,
        aoi=aoi,
        year_months=((2023, 7), (2022, 8)),
        resolution=500.0,
        output_path=tmp_path / "prepared_store",
        auth=auth,
    )

    assert captured["resolved"] == {"sensor": "auto", "aoi": aoi}
    build_kwargs = captured["build_kwargs"]
    assert build_kwargs["bounds"] == (1.0, 2.0, 3.0, 4.0)
    assert build_kwargs["crs"] == "EPSG:4326"
    assert build_kwargs["resolution"] == 500.0
    assert build_kwargs["year_months"] == ((2023, 7), (2022, 8))
    _, _, grid = captured["write"]
    assert grid == MonthlyCompositeStoreGridSpec.from_bounds((1.0, 2.0, 3.0, 4.0), crs="EPSG:4326", resolution=500.0)
    assert result.store_path == tmp_path / "prepared_store"
    assert result.periods == ("2022-08", "2023-07")
    assert result.period_count == 2
    assert result.source_name == "prepared-test"
    assert result.source_band_names == ("B02",)
    assert result.representation == "kernel_weights"
