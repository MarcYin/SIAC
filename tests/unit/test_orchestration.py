"""Public-orchestration tests for the greenfield API layer."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import xarray as xr

from siac.api.public import SIAC
from siac.catalog import SENTINEL2A_CONFIG
from siac.config import SIACConfig
from siac.domain import CorrectionResult, GeometryAngles


def _da(value: float, shape: tuple[int, int] = (2, 2)) -> xr.DataArray:
    return xr.DataArray(np.full(shape, value, dtype=np.float32), dims=["y", "x"])


def _geometry(shape: tuple[int, int] = (2, 2)) -> GeometryAngles:
    return GeometryAngles(sza=_da(0.5, shape), saa=_da(1.1, shape), vza=_da(0.2, shape), vaa=_da(1.4, shape))


def _toa_dataset(shape: tuple[int, int] = (2, 2)) -> xr.Dataset:
    return xr.Dataset({"B02": _da(0.2, shape), "B03": _da(0.25, shape)})


def _empty_auth():
    from siac.adapters.auth import CredentialManager

    return CredentialManager()


def _final_result() -> CorrectionResult:
    return CorrectionResult(
        boa=_toa_dataset(),
        boa_unc=None,
        aot=_da(0.2),
        tcwv=_da(2.3),
        cloud_mask=xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"]),
        metadata={"ok": True},
    )


class TestProcessDelegatesToSceneWorkflow:
    def test_delegates_to_process_scene(self, monkeypatch, tmp_path: Path):
        cfg = SIACConfig(sensor="s2")
        monkeypatch.setattr("siac.api.public.CredentialManager.from_config", lambda _config: _empty_auth())

        captured: dict[str, object] = {}

        def _fake_process_scene(config, input_path, **kwargs):  # noqa: ANN001
            captured["config"] = config
            captured["input_path"] = input_path
            captured["kwargs"] = kwargs
            return _final_result()

        monkeypatch.setattr("siac.api.public.process_scene", _fake_process_scene)

        siac_obj = SIAC(cfg)
        result = siac_obj.process(tmp_path / "in.SAFE")

        assert result.metadata["ok"] is True
        assert captured["input_path"] == tmp_path / "in.SAFE"
        assert captured["kwargs"]["aoi_resolver"] == siac_obj._resolve_aoi

    def test_output_saved_when_path_provided(self, monkeypatch, tmp_path: Path):
        cfg = SIACConfig(sensor="s2")
        monkeypatch.setattr("siac.api.public.CredentialManager.from_config", lambda _config: _empty_auth())
        monkeypatch.setattr("siac.api.public.process_scene", lambda *_args, **_kwargs: _final_result())

        siac_obj = SIAC(cfg)
        saved: dict[str, Path] = {}

        monkeypatch.setattr(siac_obj, "_save_output", lambda _res, path: saved.setdefault("path", Path(path)))

        siac_obj.process(tmp_path / "in.SAFE", output_path=tmp_path / "out")
        assert saved["path"] == tmp_path / "out"

    def test_resolve_aoi_uses_runtime_default(self, monkeypatch):
        cfg = SIACConfig(sensor="s2")
        monkeypatch.setattr("siac.api.public.CredentialManager.from_config", lambda _config: _empty_auth())
        monkeypatch.setattr("siac.api.public.AOI.from_raster", lambda band: ("raster", band.name))

        siac_obj = SIAC(cfg)
        assert siac_obj._resolve_aoi(_toa_dataset()) == ("raster", "B02")


class TestConvenienceHelpers:
    def test_get_brdf_provider_is_cached(self, monkeypatch):
        cfg = SIACConfig(sensor="s2")
        monkeypatch.setattr("siac.api.public.CredentialManager.from_config", lambda _config: _empty_auth())
        monkeypatch.setattr("siac.api.public.resolve_brdf_provider", lambda *_args, **_kwargs: "provider")

        siac_obj = SIAC(cfg)
        assert siac_obj._get_brdf_provider() == "provider"
        siac_obj._brdf_provider = "cached"
        assert siac_obj._get_brdf_provider() == "cached"

    def test_get_rt_model_is_cached(self, monkeypatch):
        cfg = SIACConfig(sensor="s2")
        monkeypatch.setattr("siac.api.public.CredentialManager.from_config", lambda _config: _empty_auth())
        monkeypatch.setattr("siac.api.public.resolve_rt_model_for_pipeline", lambda *_args, **_kwargs: "rt")

        siac_obj = SIAC(cfg)
        assert siac_obj._get_rt_model(SENTINEL2A_CONFIG) == "rt"
        siac_obj._rt_model = "cached"
        assert siac_obj._get_rt_model(SENTINEL2A_CONFIG) == "cached"

    def test_get_surface_prior_uses_surface_resolver(self, monkeypatch):
        cfg = SIACConfig(sensor="s2")
        monkeypatch.setattr("siac.api.public.CredentialManager.from_config", lambda _config: _empty_auth())
        monkeypatch.setattr("siac.api.public.resolve_brdf_provider", lambda *_args, **_kwargs: "provider")

        captured: dict[str, object] = {}

        def _fake_surface(_config, auth=None):  # noqa: ANN001
            captured["auth"] = auth

            def _surface(obs, _atmo, _rt, resolution):  # noqa: ANN001
                captured["obs"] = obs
                captured["resolution"] = resolution
                return "surface"

            return _surface

        monkeypatch.setattr("siac.api.public.resolve_surface_prior_provider", _fake_surface)

        siac_obj = SIAC(cfg)
        out = siac_obj._get_surface_prior(
            SimpleNamespace(get_bounds=lambda: (0.0, 0.0, 1.0, 1.0), crs="EPSG:4326"),
            _geometry(),
            {"observation_time": datetime(2026, 1, 1), "sensor_config": SENTINEL2A_CONFIG},
        )

        assert out == "surface"
        assert captured["resolution"] == 500.0

    def test_solve_atmosphere_delegates_to_solver(self, monkeypatch):
        cfg = SIACConfig(sensor="s2")
        monkeypatch.setattr("siac.api.public.CredentialManager.from_config", lambda _config: _empty_auth())
        monkeypatch.setattr("siac.api.public.resolve_solver", lambda _cfg: lambda inputs, _config: {"bands": inputs.bands})

        siac_obj = SIAC(cfg)
        out = siac_obj._solve_atmosphere(
            toa="toa",
            surface_prior="surface",
            geometry=_geometry(),
            atmo_prior="atmo",
            rt_model="rt",
            cloud_mask="mask",
            sensor_config=SENTINEL2A_CONFIG,
        )
        assert len(out["bands"]) == 6
