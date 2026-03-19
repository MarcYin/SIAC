"""Test that SIAC.process() correctly orchestrates the pipeline.

Verifies:
- process() delegates to run_pipeline()
- The preprocessor closure correctly wraps dict -> ObservationBundle
- Module-level resolvers are called
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import xarray as xr

from siac.config import SIACConfig
from siac.domain import (
    SENTINEL2A_CONFIG,
    CorrectionResult,
    GeometryAngles,
    ObservationBundle,
)
from siac.siac import SIAC


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


class TestProcessDelegatesToPipeline:
    """SIAC.process() must resolve components and call run_pipeline()."""

    def test_delegates_to_run_pipeline(self, monkeypatch, tmp_path: Path):
        cfg = SIACConfig(sensor="s2")
        monkeypatch.setattr("siac.siac.CredentialManager.from_config", lambda _config: _empty_auth())
        siac_obj = SIAC(cfg)

        class _PP:
            sensor_config = SENTINEL2A_CONFIG
            def preprocess(self, path):
                return {
                    "toa": _toa_dataset(),
                    "geometry": _geometry(),
                    "cloud_mask": xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"]),
                    "metadata": {"observation_time": datetime(2026, 1, 1)},
                }

        monkeypatch.setattr("siac.siac.get_preprocessor", lambda _s: _PP())
        monkeypatch.setattr(siac_obj, "_resolve_aoi", lambda _toa: SimpleNamespace(get_bounds=lambda: (0, 0, 1, 1), crs="EPSG:4326"))
        monkeypatch.setattr(siac_obj, "_get_rt_model", lambda _sc: "rt")

        def _fake_atmo(_c, auth=None):
            _ = auth
            return "atmo"

        def _fake_surface(_c, auth=None):
            _ = auth
            return "surf"

        monkeypatch.setattr("siac.siac._resolve_atmo_provider", _fake_atmo)
        monkeypatch.setattr("siac.siac._resolve_surface_prior_provider", _fake_surface)
        monkeypatch.setattr("siac.siac._resolve_grid_assembler", lambda: "grid")
        monkeypatch.setattr("siac.siac._resolve_solver", lambda _c: "solver")
        monkeypatch.setattr("siac.siac._resolve_corrector", lambda _c: "corrector")
        monkeypatch.setattr("siac.siac._resolve_rt_model_for_pipeline", lambda *_a, **_kw: "rt")

        captured = {}
        final = _final_result()

        def _mock_run_pipeline(input_path, aoi, config, **kwargs):
            captured.update(kwargs)
            captured["input_path"] = input_path
            return final

        monkeypatch.setattr("siac.siac.run_pipeline", _mock_run_pipeline)

        result = siac_obj.process(tmp_path / "in.SAFE")
        assert result is final
        assert captured["rt_model"] == "rt"
        assert callable(captured["preprocessor"])

    def test_preprocessor_closure_returns_observation_bundle(self, monkeypatch, tmp_path: Path):
        """The _preprocessor_fn closure wraps the dict into ObservationBundle."""
        cfg = SIACConfig(sensor="s2")
        monkeypatch.setattr("siac.siac.CredentialManager.from_config", lambda _config: _empty_auth())
        siac_obj = SIAC(cfg)

        obs_time = datetime(2026, 6, 15)
        class _PP:
            sensor_config = SENTINEL2A_CONFIG
            def preprocess(self, path):
                return {
                    "toa": _toa_dataset(),
                    "geometry": _geometry(),
                    "cloud_mask": xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"]),
                    "metadata": {"observation_time": obs_time},
                }

        monkeypatch.setattr("siac.siac.get_preprocessor", lambda _s: _PP())
        monkeypatch.setattr(siac_obj, "_resolve_aoi", lambda _toa: SimpleNamespace(get_bounds=lambda: (10, 20, 30, 40), crs="EPSG:32633"))
        monkeypatch.setattr(siac_obj, "_get_rt_model", lambda _sc: "rt")

        def _fake_atmo(_c, auth=None):
            _ = auth
            return "atmo"

        def _fake_surface(_c, auth=None):
            _ = auth
            return "surf"

        monkeypatch.setattr("siac.siac._resolve_atmo_provider", _fake_atmo)
        monkeypatch.setattr("siac.siac._resolve_surface_prior_provider", _fake_surface)
        monkeypatch.setattr("siac.siac._resolve_grid_assembler", lambda: "grid")
        monkeypatch.setattr("siac.siac._resolve_solver", lambda _c: "solver")
        monkeypatch.setattr("siac.siac._resolve_corrector", lambda _c: "corrector")
        monkeypatch.setattr("siac.siac._resolve_rt_model_for_pipeline", lambda *_a, **_kw: "rt")

        captured_preprocessor = {}

        def _mock_run_pipeline(input_path, aoi, config, **kwargs):
            # Exercise the preprocessor closure
            pp = kwargs["preprocessor"]
            bundle = pp(input_path, None)
            captured_preprocessor["bundle"] = bundle
            return _final_result()

        monkeypatch.setattr("siac.siac.run_pipeline", _mock_run_pipeline)

        siac_obj.process(tmp_path / "in.SAFE")
        bundle = captured_preprocessor["bundle"]
        assert isinstance(bundle, ObservationBundle)
        assert bundle.crs == "EPSG:32633"
        assert bundle.bounds == (10, 20, 30, 40)
        assert bundle.sensor_config is SENTINEL2A_CONFIG

    def test_output_saved_when_path_provided(self, monkeypatch, tmp_path: Path):
        cfg = SIACConfig(sensor="s2")
        monkeypatch.setattr("siac.siac.CredentialManager.from_config", lambda _config: _empty_auth())
        siac_obj = SIAC(cfg)

        class _PP:
            sensor_config = SENTINEL2A_CONFIG
            def preprocess(self, path):
                return {
                    "toa": _toa_dataset(),
                    "geometry": _geometry(),
                    "cloud_mask": xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"]),
                    "metadata": {"observation_time": datetime(2026, 1, 1)},
                }

        monkeypatch.setattr("siac.siac.get_preprocessor", lambda _s: _PP())
        monkeypatch.setattr(siac_obj, "_resolve_aoi", lambda _toa: SimpleNamespace(get_bounds=lambda: (0, 0, 1, 1), crs="EPSG:4326"))
        monkeypatch.setattr(siac_obj, "_get_rt_model", lambda _sc: "rt")

        def _fake_atmo(_c, auth=None):
            _ = auth
            return "atmo"

        def _fake_surface(_c, auth=None):
            _ = auth
            return "surf"

        monkeypatch.setattr("siac.siac._resolve_atmo_provider", _fake_atmo)
        monkeypatch.setattr("siac.siac._resolve_surface_prior_provider", _fake_surface)
        monkeypatch.setattr("siac.siac._resolve_grid_assembler", lambda: "grid")
        monkeypatch.setattr("siac.siac._resolve_solver", lambda _c: "solver")
        monkeypatch.setattr("siac.siac._resolve_corrector", lambda _c: "corrector")
        monkeypatch.setattr("siac.siac._resolve_rt_model_for_pipeline", lambda *_a, **_kw: "rt")

        def _fake_run_pipeline(*_args, **_kwargs):
            return _final_result()

        monkeypatch.setattr("siac.siac.run_pipeline", _fake_run_pipeline)

        saved = {}

        def _fake_save_output(_res, p):
            saved.update(path=Path(p))

        monkeypatch.setattr(siac_obj, "_save_output", _fake_save_output)

        siac_obj.process(tmp_path / "in.SAFE", output_path=tmp_path / "out")
        assert saved["path"] == tmp_path / "out"
