"""High-value API and helper coverage for the greenfield structure."""

from __future__ import annotations

from datetime import date, datetime
from types import SimpleNamespace
from typing import TYPE_CHECKING

import numpy as np
import pytest
import xarray as xr

from siac.adapters.data.s2_data_source import S2Query
from siac.api.public import (
    SIAC,
    apply_s2_query_defaults,
    coerce_date,
    coerce_s2_query,
    process_landsat8,
    process_sentinel2,
    resolve_s2_input,
    search_sentinel2,
    siac_process,
)
from siac.catalog import SENTINEL2A_CONFIG
from siac.config import SIACConfig
from siac.domain import AtmosphericState, CorrectionResult, GeometryAngles
from siac.errors import DataNotFoundError

if TYPE_CHECKING:
    from pathlib import Path


def _da(value: float, shape: tuple[int, int] = (2, 2)) -> xr.DataArray:
    return xr.DataArray(np.full(shape, value, dtype=np.float32), dims=["y", "x"])


def _geometry(shape: tuple[int, int] = (2, 2)) -> GeometryAngles:
    return GeometryAngles(sza=_da(0.5, shape), saa=_da(1.1, shape), vza=_da(0.2, shape), vaa=_da(1.4, shape))


def _atmo_state(shape: tuple[int, int] = (2, 2)) -> AtmosphericState:
    return AtmosphericState(
        aot=_da(0.1, shape),
        tcwv=_da(2.0, shape),
        tco3=_da(0.3, shape),
        aot_unc=_da(0.01, shape),
        tcwv_unc=_da(0.2, shape),
        tco3_unc=_da(0.01, shape),
        elevation=_da(0.1, shape),
    )


def _toa_dataset(shape: tuple[int, int] = (2, 2)) -> xr.Dataset:
    return xr.Dataset(
        {
            "B02": _da(0.2, shape),
            "B03": _da(0.25, shape),
            "B04": _da(0.3, shape),
            "B08": _da(0.35, shape),
            "B11": _da(0.15, shape),
            "B12": _da(0.12, shape),
        }
    )


def _empty_auth():
    from siac.adapters.auth import CredentialManager

    return CredentialManager()


def _earthdata_auth():
    auth = _empty_auth()
    auth.set_credentials("earthdata", key="u", secret="p")
    return auth


def _aws_auth():
    auth = _empty_auth()
    auth.set_credentials("aws", key="AK", secret="SK")
    return auth


def test_siac_from_helpers(monkeypatch):
    cfg = SIACConfig(sensor="s2")
    monkeypatch.setattr("siac.api.public.CredentialManager.from_config", lambda _config: _empty_auth())
    monkeypatch.setattr("siac.api.public.SIACConfig.from_file", classmethod(lambda _cls, _path: cfg))

    siac_from_file = SIAC.from_config("dummy.toml")
    siac_default = SIAC.from_defaults(sensor="l8")

    assert isinstance(siac_from_file, SIAC)
    assert siac_from_file.config is cfg
    assert siac_default.config.sensor == "l8"


def test_process_delegates_to_scene_workflow(monkeypatch, tmp_path: Path):
    cfg = SIACConfig(sensor="auto")
    monkeypatch.setattr("siac.api.public.CredentialManager.from_config", lambda _config: _empty_auth())
    siac_obj = SIAC(cfg)

    captured: dict[str, object] = {}
    final_result = CorrectionResult(
        boa=_toa_dataset(),
        boa_unc=None,
        aot=_da(0.2),
        tcwv=_da(2.3),
        cloud_mask=xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"]),
        metadata={"ok": True},
    )

    def _fake_process_scene(config, input_path, **kwargs):  # noqa: ANN001
        captured["config"] = config
        captured["input_path"] = input_path
        captured["kwargs"] = kwargs
        return final_result

    monkeypatch.setattr("siac.api.public.process_scene", _fake_process_scene)

    result = siac_obj.process(tmp_path / "input.SAFE", output_path=tmp_path / "out")

    assert result is final_result
    assert captured["input_path"] == tmp_path / "input.SAFE"
    assert callable(captured["kwargs"]["aoi_resolver"])


def test_resolve_aoi_branches(monkeypatch):
    monkeypatch.setattr("siac.api.public.CredentialManager.from_config", lambda _config: _empty_auth())

    cfg_bounds = SIACConfig(sensor="s2")
    cfg_bounds.aoi = [1.0, 2.0, 3.0, 4.0]
    siac_bounds = SIAC(cfg_bounds)
    monkeypatch.setattr("siac.api.public.AOI.from_bounds", lambda bounds: ("bounds", bounds))
    assert siac_bounds._resolve_aoi(_toa_dataset()) == ("bounds", (1.0, 2.0, 3.0, 4.0))

    cfg_geojson = SIACConfig(sensor="s2", aoi="/tmp/aoi.geojson")
    siac_geojson = SIAC(cfg_geojson)
    monkeypatch.setattr("siac.api.public.AOI.from_geojson", lambda spec: ("geojson", spec))
    assert siac_geojson._resolve_aoi(_toa_dataset()) == ("geojson", "/tmp/aoi.geojson")

    cfg_default = SIACConfig(sensor="s2")
    siac_default = SIAC(cfg_default)
    monkeypatch.setattr("siac.api.public.AOI.from_raster", lambda band: ("raster", band.name))
    assert siac_default._resolve_aoi(_toa_dataset()) == ("raster", "B02")


def test_get_atmospheric_prior_providers_and_fallback(monkeypatch):
    monkeypatch.setattr("siac.api.public.CredentialManager.from_config", lambda _config: _earthdata_auth())
    calls: list[tuple[object, object, object, object]] = []

    def _provider_factory(config, auth=None):  # noqa: ANN001
        calls.append(("factory", config.atmo_prior.provider, auth, None))

        def _provider(bounds, crs, obs_time, resolution):
            calls.append((bounds, crs, obs_time, resolution))
            return (bounds, crs, obs_time, resolution)

        return _provider

    monkeypatch.setattr("siac.api.public.resolve_atmo_provider", _provider_factory)

    fake_aoi = SimpleNamespace(get_bounds=lambda: (1.0, 2.0, 3.0, 4.0), crs="EPSG:4326")
    metadata = {"observation_time": datetime(2026, 1, 1, 0, 0, 0)}

    cfg = SIACConfig(sensor="s2", providers={"atmo": {"provider": "cams"}})
    siac_obj = SIAC(cfg)
    out = siac_obj._get_atmospheric_prior(fake_aoi, metadata)
    assert out[0] == (1.0, 2.0, 3.0, 4.0)

    cfg_unknown = SIACConfig(sensor="s2", providers={"atmo": {"provider": "cams"}})
    cfg_unknown.atmo_prior.provider = "era5"
    siac_unknown = SIAC(cfg_unknown)
    warnings: list[str] = []
    monkeypatch.setattr("siac.api.public.logger.warning", lambda msg, *_args: warnings.append(str(msg)))
    _ = siac_unknown._get_atmospheric_prior(fake_aoi, metadata)
    assert warnings


def test_get_surface_prior_and_brdf_provider_paths(monkeypatch):
    monkeypatch.setattr("siac.api.public.CredentialManager.from_config", lambda _config: _earthdata_auth())

    cfg = SIACConfig(sensor="s2", providers={"brdf": {"provider": "mcd43"}})
    siac_obj = SIAC(cfg)

    monkeypatch.setattr("siac.api.public.resolve_brdf_provider", lambda *_args, **_kwargs: "provider")

    provider_calls: dict[str, object] = {}

    def _fake_surface_provider(_cfg, auth=None):  # noqa: ANN001
        provider_calls["auth"] = auth

        def _surface(obs, atmo_prior, rt_model, resolution):  # noqa: ANN001
            provider_calls["obs"] = obs
            provider_calls["resolution"] = resolution
            return "surface"

        return _surface

    monkeypatch.setattr("siac.api.public.resolve_surface_prior_provider", _fake_surface_provider)

    fake_aoi = SimpleNamespace(get_bounds=lambda: (0.0, 0.0, 1.0, 1.0), crs="EPSG:4326")
    geom = _geometry()
    out = siac_obj._get_surface_prior(
        fake_aoi,
        geom,
        {"observation_time": datetime(2026, 1, 2), "sensor_config": SENTINEL2A_CONFIG},
    )

    assert out == "surface"
    assert provider_calls["resolution"] == 500.0

    siac_obj._brdf_provider = "cached"
    assert siac_obj._get_brdf_provider() == "cached"


def test_get_rt_model_and_solver_delegate(monkeypatch):
    monkeypatch.setattr("siac.api.public.CredentialManager.from_config", lambda _config: _empty_auth())

    cfg = SIACConfig(sensor="s2")
    siac_obj = SIAC(cfg)

    monkeypatch.setattr("siac.api.public.resolve_rt_model_for_pipeline", lambda *_args, **_kwargs: "rt-model")
    assert siac_obj._get_rt_model(SENTINEL2A_CONFIG) == "rt-model"
    siac_obj._rt_model = "cached"
    assert siac_obj._get_rt_model(SENTINEL2A_CONFIG) == "cached"

    monkeypatch.setattr("siac.api.public.resolve_solver", lambda _cfg: lambda inputs, _config: {"bands": inputs.bands})
    out = siac_obj._solve_atmosphere(
        toa="toa",
        surface_prior="surface",
        geometry=_geometry(),
        atmo_prior=_atmo_state(),
        rt_model="rt",
        cloud_mask="mask",
        sensor_config=SENTINEL2A_CONFIG,
    )
    assert len(out["bands"]) == 6


def test_sensor_wrappers_delegate(monkeypatch):
    calls: list[object] = []

    class _Runner:
        def process(self, input_path, output_path=None):  # noqa: ANN001
            calls.append((input_path, output_path))
            return "ok"

    monkeypatch.setattr("siac.api.public.SIAC.from_defaults", lambda sensor: calls.append(sensor) or _Runner())

    assert process_sentinel2("in1", "out1") == "ok"
    assert process_landsat8("in2", "out2") == "ok"
    assert calls[0] == "s2"
    assert calls[2] == "l8"


def test_public_s2_helper_exports_and_local_errors(tmp_path: Path):
    cfg = SIACConfig(sensor="s2", providers={"s2": {"backend": "local"}})

    assert coerce_date(None) is None
    assert coerce_date(date(2026, 1, 2)) == date(2026, 1, 2)
    assert coerce_date(datetime(2026, 1, 2, 8, 0, 0)) == date(2026, 1, 2)
    assert coerce_date("2026-01-03") == date(2026, 1, 3)

    q = S2Query.from_tile_date("T50QLD_20260102")
    q.max_cloud_cover = 100.0
    q.processing_level = "L1C"
    out_q = apply_s2_query_defaults(q, config=SIACConfig(sensor="s2", providers={"s2": {"max_cloud_cover": 33.0, "processing_level": "L2A"}}))
    assert out_q.max_cloud_cover == 33.0
    assert out_q.processing_level == "L2A"
    assert coerce_s2_query(out_q, config=cfg) is not out_q

    with pytest.raises(DataNotFoundError, match="backend is 'local'"):
        resolve_s2_input("missing.SAFE", cfg)

    with pytest.raises(ValueError, match="does not support backend='local'"):
        search_sentinel2(tile="50QLD", date_value="2026-01-02", backend="local")


def test_siac_process_delegates_to_scene(monkeypatch, tmp_path: Path):
    cfg = SIACConfig(sensor="s2")
    calls: dict[str, object] = {}

    def _fake_process_scene(config, input_path, **kwargs):  # noqa: ANN001
        calls["config"] = config
        calls["input_path"] = input_path
        calls["kwargs"] = kwargs
        return "pipeline-result"

    monkeypatch.setattr("siac.api.public.process_scene", _fake_process_scene)

    out = siac_process(cfg, tmp_path / "in.SAFE")
    assert out == "pipeline-result"
    assert calls["input_path"] == tmp_path / "in.SAFE"
