"""High-coverage tests for siac.siac orchestration and resolver helpers."""

from __future__ import annotations

import sys
from datetime import date, datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

from siac.core.auth import CredentialManager
from siac.core.config import SIACConfig
from siac.core.exceptions import DataNotFoundError
from siac.core.types import SENTINEL2A_CONFIG, AtmosphericState, CorrectionResult, GeometryAngles
from siac.io.s2_data_source import S2Query
from siac.siac import (
    SIAC,
    _apply_s2_query_defaults,
    _coerce_date,
    _coerce_s2_query,
    _resolve_atmo_provider,
    _resolve_corrector,
    _resolve_preprocessor,
    _resolve_rt_model_for_pipeline,
    _resolve_s2_backend,
    _resolve_solver,
    _resolve_surface_prior_provider,
    process_landsat8,
    process_sentinel2,
    search_sentinel2,
    siac_process,
)


# NOTE: _get_brdf_provider is an SIAC instance method; small wrapper keeps callsites clear.
def _resolve_brdf_provider(siac_obj: SIAC):
    return siac_obj._get_brdf_provider()


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


def _empty_auth() -> CredentialManager:
    return CredentialManager()


def _earthdata_auth() -> CredentialManager:
    auth = CredentialManager()
    auth.set_credentials("earthdata", key="u", secret="p")
    return auth


def _aws_auth() -> CredentialManager:
    auth = CredentialManager()
    auth.set_credentials("aws", key="AK", secret="SK")
    return auth


def test_siac_from_helpers(monkeypatch):
    cfg = SIACConfig(sensor="s2")
    monkeypatch.setattr("siac.siac.CredentialManager.from_config", lambda _config: _empty_auth())
    monkeypatch.setattr("siac.siac.SIACConfig.from_yaml", classmethod(lambda _cls, _path: cfg))

    siac_from_yaml = SIAC.from_config("dummy.yaml")
    siac_default = SIAC.from_defaults(sensor="l8")

    assert isinstance(siac_from_yaml, SIAC)
    assert siac_from_yaml.config is cfg
    assert isinstance(siac_default, SIAC)
    assert siac_default.config.sensor == "l8"


def test_process_happy_path_calls_all_steps(monkeypatch, tmp_path: Path):
    cfg = SIACConfig(sensor="auto")
    monkeypatch.setattr("siac.siac.CredentialManager.from_config", lambda _config: _empty_auth())
    siac_obj = SIAC(cfg)

    obs_time = datetime(2026, 1, 2, 3, 4, 5)
    preprocess_result = {
        "toa": _toa_dataset(),
        "geometry": _geometry(),
        "cloud_mask": xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"]),
        "metadata": {"observation_time": obs_time},
    }

    class _FakePreprocessor:
        sensor_config = SENTINEL2A_CONFIG

        def preprocess(self, input_path):
            return preprocess_result

    final_result = CorrectionResult(
        boa=_toa_dataset(),
        boa_unc=None,
        aot=_da(0.2),
        tcwv=_da(2.3),
        cloud_mask=xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"]),
        metadata={"ok": True},
    )

    monkeypatch.setattr("siac.siac.detect_sensor", lambda _path: "s2")
    monkeypatch.setattr("siac.siac.get_preprocessor", lambda _sensor: _FakePreprocessor())
    monkeypatch.setattr(siac_obj, "_resolve_aoi", lambda _toa: SimpleNamespace(get_bounds=lambda: (0, 0, 1, 1), crs="EPSG:4326"))
    monkeypatch.setattr(siac_obj, "_get_rt_model", lambda _sc: "rt")

    # Mock the module-level resolvers that process() now calls
    def _fake_atmo_provider(_cfg, auth=None):
        _ = auth
        return "atmo-provider"

    def _fake_surface_provider(_cfg, auth=None):
        _ = auth
        return "surface-prior-provider"

    monkeypatch.setattr("siac.siac._resolve_atmo_provider", _fake_atmo_provider)
    monkeypatch.setattr("siac.siac._resolve_surface_prior_provider", _fake_surface_provider)
    monkeypatch.setattr("siac.siac._resolve_grid_assembler", lambda: "grid-assembler")
    monkeypatch.setattr("siac.siac._resolve_solver", lambda _cfg: "solver")
    monkeypatch.setattr("siac.siac._resolve_corrector", lambda _cfg: "corrector")

    # Mock run_pipeline to verify process() delegates to it
    pipeline_calls = {}

    def _fake_run_pipeline(input_path, aoi, config, **kwargs):
        pipeline_calls["input_path"] = input_path
        pipeline_calls["kwargs"] = kwargs
        return final_result

    monkeypatch.setattr("siac.siac.run_pipeline", _fake_run_pipeline)

    saved = {}

    def _fake_save(result, output_path):
        saved["path"] = Path(output_path)

    monkeypatch.setattr(siac_obj, "_save_output", _fake_save)

    result = siac_obj.process(tmp_path / "input.SAFE", output_path=tmp_path / "out")

    assert result is final_result
    assert saved["path"] == tmp_path / "out"
    # Verify pipeline was called with resolved components
    assert pipeline_calls["input_path"] == tmp_path / "input.SAFE"
    assert pipeline_calls["kwargs"]["rt_model"] == "rt"
    assert callable(pipeline_calls["kwargs"]["preprocessor"])


def test_resolve_aoi_branches(monkeypatch):
    monkeypatch.setattr("siac.siac.CredentialManager.from_config", lambda _config: _empty_auth())

    cfg_bounds = SIACConfig(sensor="s2")
    cfg_bounds.aoi = [1.0, 2.0, 3.0, 4.0]
    siac_bounds = SIAC(cfg_bounds)
    monkeypatch.setattr("siac.siac.AOI.from_bounds", lambda bounds: ("bounds", bounds))
    assert siac_bounds._resolve_aoi(_toa_dataset()) == ("bounds", (1.0, 2.0, 3.0, 4.0))

    cfg_geojson = SIACConfig(sensor="s2", aoi="/tmp/aoi.geojson")
    siac_geojson = SIAC(cfg_geojson)
    monkeypatch.setattr("siac.siac.AOI.from_geojson", lambda spec: ("geojson", spec))
    assert siac_geojson._resolve_aoi(_toa_dataset()) == ("geojson", "/tmp/aoi.geojson")

    cfg_default = SIACConfig(sensor="s2")
    siac_default = SIAC(cfg_default)
    monkeypatch.setattr("siac.siac.AOI.from_raster", lambda band: ("raster", band.name))
    assert siac_default._resolve_aoi(_toa_dataset()) == ("raster", "B02")


def test_get_atmospheric_prior_providers_and_fallback(monkeypatch):
    monkeypatch.setattr("siac.siac.CredentialManager.from_config", lambda _config: _earthdata_auth())

    calls = []

    class _FakeProvider:
        def __init__(self, *args, **kwargs):
            calls.append((args, kwargs))

        def get_prior(self, bounds, crs, obs_time, resolution):
            return (bounds, crs, obs_time, resolution)

    fake_aoi = SimpleNamespace(get_bounds=lambda: (1.0, 2.0, 3.0, 4.0), crs="EPSG:4326")
    metadata = {"observation_time": datetime(2026, 1, 1, 0, 0, 0)}

    cfg_cams = SIACConfig(
        sensor="s2",
        atmo_prior={
            "provider": "cams",
            "data_path": "https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/",
            "cache_dir": "/tmp/cams-cache",
        },
    )
    siac_cams = SIAC(cfg_cams)
    monkeypatch.setattr("siac.siac.CAMSProvider", _FakeProvider)
    out_cams = siac_cams._get_atmospheric_prior(fake_aoi, metadata)
    assert out_cams[0] == (1.0, 2.0, 3.0, 4.0)
    assert out_cams[1] == "EPSG:4326"
    assert calls[0][1]["download_missing"] is True
    assert calls[0][0][0] == "https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/"
    assert str(calls[0][1]["cache_dir"]) == "/tmp/cams-cache"

    cfg_merra = SIACConfig(sensor="s2", atmo_prior={"provider": "merra2", "cache_dir": "/tmp/cache"})
    siac_merra = SIAC(cfg_merra)
    monkeypatch.setattr("siac.priors.atmospheric.merra2.MERRA2Provider", _FakeProvider)
    assert siac_merra._get_atmospheric_prior(fake_aoi, metadata)[3] == 10.0
    assert calls[1][1]["source"].earthdata_username == "u"

    cfg_mcd19 = SIACConfig(sensor="s2", atmo_prior={"provider": "mcd19", "cache_dir": "/tmp/cache"})
    siac_mcd19 = SIAC(cfg_mcd19)
    monkeypatch.setattr("siac.priors.atmospheric.mcd19_earthaccess.MCD19AODProvider", _FakeProvider)
    assert siac_mcd19._get_atmospheric_prior(fake_aoi, metadata)[3] == 10.0

    cfg_vnp19 = SIACConfig(sensor="s2", atmo_prior={"provider": "vnp19", "cache_dir": "/tmp/cache"})
    siac_vnp19 = SIAC(cfg_vnp19)
    monkeypatch.setattr("siac.priors.atmospheric.mcd19_earthaccess.VNP19AODProvider", _FakeProvider)
    assert siac_vnp19._get_atmospheric_prior(fake_aoi, metadata)[3] == 10.0

    cfg_unknown = SIACConfig(sensor="s2", atmo_prior={"provider": "cams"})
    cfg_unknown.atmo_prior.provider = "era5"
    siac_unknown = SIAC(cfg_unknown)
    warnings = []
    monkeypatch.setattr("siac.siac.CAMSProvider", _FakeProvider)
    monkeypatch.setattr("siac.siac.logger.warning", lambda msg: warnings.append(msg))
    _ = siac_unknown._get_atmospheric_prior(fake_aoi, metadata)
    assert warnings


def test_get_surface_prior_and_brdf_provider_paths(monkeypatch):
    monkeypatch.setattr("siac.siac.CredentialManager.from_config", lambda _config: _earthdata_auth())

    cfg = SIACConfig(sensor="s2", brdf={"provider": "mcd43"})
    siac_obj = SIAC(cfg)

    provider_calls = []

    class _FakeBRDFProvider:
        def __init__(self, cache_dir=None, source=None):
            self.cache_dir = cache_dir
            self.source = source
            self.source_bands = [SENTINEL2A_CONFIG.get_band("B01"), SENTINEL2A_CONFIG.get_band("B02")]

        def get_brdf_parameters(self, **kwargs):
            provider_calls.append(kwargs)
            return "weights"

    class _FakeDeriver:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

        def compute_surface_prior(self, weights, geometry, **kwargs):
            return (weights, geometry, kwargs)

    monkeypatch.setattr("siac.priors.brdf.mcd43_earthaccess.MCD43EarthAccessProvider", _FakeBRDFProvider)
    monkeypatch.setattr("siac.siac.KernelModelDeriver", _FakeDeriver)

    fake_aoi = SimpleNamespace(get_bounds=lambda: (0.0, 0.0, 1.0, 1.0), crs="EPSG:4326")
    geom = _geometry()
    out = siac_obj._get_surface_prior(
        fake_aoi,
        geom,
        {
            "observation_time": datetime(2026, 1, 2),
            "sensor_config": SENTINEL2A_CONFIG,
        },
    )

    assert out[0] == "weights"
    assert out[1] is geom
    assert provider_calls and [band.name for band in provider_calls[0]["bands"]] == ["B01", "B02"]
    assert siac_obj._brdf_provider.source.earthdata_username == "u"

    cached = object()
    siac_obj._brdf_provider = cached
    assert _resolve_brdf_provider(siac_obj) is cached


def test_get_brdf_provider_other_branches(monkeypatch):
    monkeypatch.setattr("siac.siac.CredentialManager.from_config", lambda _config: _empty_auth())

    cfg_vnp43 = SIACConfig(sensor="s2", brdf={"provider": "vnp43"})
    siac_vnp43 = SIAC(cfg_vnp43)
    monkeypatch.setattr(
        "siac.priors.brdf.vnp43_earthaccess.VNP43EarthAccessProvider",
        lambda **_kwargs: "vnp",
    )
    assert siac_vnp43._get_brdf_provider() == "vnp"

    cfg_mcd19 = SIACConfig(sensor="s2", brdf={"provider": "mcd19"})
    siac_mcd19 = SIAC(cfg_mcd19)
    monkeypatch.setattr(
        "siac.priors.brdf.mcd43_earthaccess.MCD19EarthAccessProvider",
        lambda **_kwargs: "mcd19",
    )
    assert siac_mcd19._get_brdf_provider() == "mcd19"

    cfg_gee = SIACConfig(sensor="s2", brdf={"provider": "mcd43"})
    cfg_gee.brdf.provider = "gee"
    siac_gee = SIAC(cfg_gee)
    fake_gee_module = SimpleNamespace(GEEBRDFProvider=lambda: "gee")
    monkeypatch.setitem(sys.modules, "siac.priors.brdf.gee_stub", fake_gee_module)
    assert siac_gee._get_brdf_provider() == "gee"

    cfg_bad = SIACConfig(sensor="s2", brdf={"provider": "mcd43"})
    cfg_bad.brdf.provider = "invalid"
    siac_bad = SIAC(cfg_bad)
    with pytest.raises(ValueError, match="Unknown BRDF provider"):
        _ = siac_bad._get_brdf_provider()


def test_get_rt_model_branches(monkeypatch):
    monkeypatch.setattr("siac.siac.CredentialManager.from_config", lambda _config: _empty_auth())

    sensor_config = SENTINEL2A_CONFIG

    cfg_emu = SIACConfig(sensor="s2", rt_model={"backend": "emulator", "emulator_dir": "/tmp/emu"})
    siac_emu = SIAC(cfg_emu)
    monkeypatch.setattr("siac.siac.TwoLayerNNEmulator", lambda **kwargs: ("emulator", kwargs))
    assert siac_emu._get_rt_model(sensor_config)[0] == "emulator"

    cfg_fallback = SIACConfig(
        sensor="s2",
        rt_model={"backend": "emulator", "lut_path": "s3://bucket/lut.zarr.zip", "lut_interpolation": "nearest"},
    )
    siac_fallback = SIAC(cfg_fallback)

    def _raise_emulator(**kwargs):
        raise RuntimeError("missing emu")

    monkeypatch.setattr("siac.siac.TwoLayerNNEmulator", _raise_emulator)
    monkeypatch.setattr("siac.siac.ZarrLUTBackend", lambda path, interpolation_method, storage_options: (path, interpolation_method, storage_options))
    out = siac_fallback._get_rt_model(sensor_config)
    assert out[0] == "s3://bucket/lut.zarr.zip"
    assert siac_fallback.config.rt_model.backend == "lut"

    cfg_lut = SIACConfig(sensor="s2", rt_model={"backend": "lut", "lut_path": "/tmp/lut.zarr"})
    siac_lut = SIAC(cfg_lut)
    monkeypatch.setattr("siac.siac.ZarrLUTBackend", lambda *_args, **_kwargs: "lut-backend")
    assert siac_lut._get_rt_model(sensor_config) == "lut-backend"

    siac_lut._rt_model = "cached"
    assert siac_lut._get_rt_model(sensor_config) == "cached"


def test_solve_atmosphere_and_save_output(monkeypatch, tmp_path: Path):
    monkeypatch.setattr("siac.siac.CredentialManager.from_config", lambda _config: _empty_auth())
    cfg = SIACConfig(sensor="s2")
    siac_obj = SIAC(cfg)

    class _FakeSolver:
        def __init__(self, config):
            self.config = config

        def solve(self, toa, surface_prior, geometry, atmo_prior, rt_model, cloud_mask, bands):
            return {"bands": bands, "toa": toa}

    monkeypatch.setattr("siac.siac.MultiGridSolver", _FakeSolver)
    out = siac_obj._solve_atmosphere(
        toa="toa",
        surface_prior="surface",
        geometry="geometry",
        atmo_prior="atmo",
        rt_model="rt",
        cloud_mask="mask",
        sensor_config=SENTINEL2A_CONFIG,
    )
    assert len(out["bands"]) == 6

    result = CorrectionResult(
        boa=_toa_dataset(),
        boa_unc=None,
        aot=_da(0.2),
        tcwv=_da(2.1),
        cloud_mask=xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"]),
        metadata={},
    )
    saved = {}

    def _fake_write_dataset(ds, path):
        saved["path"] = Path(path)

    monkeypatch.setattr("siac.io.write_dataset", _fake_write_dataset)
    siac_obj._save_output(result, tmp_path / "out")
    assert saved["path"] == tmp_path / "out" / "boa.nc"


def test_sensor_wrappers_delegate(monkeypatch):
    calls = []

    class _Runner:
        def process(self, input_path, output_path=None):
            calls.append((input_path, output_path))
            return "ok"

    def _fake_from_defaults(sensor):
        calls.append(sensor)
        return _Runner()

    monkeypatch.setattr("siac.siac.SIAC.from_defaults", _fake_from_defaults)
    assert process_sentinel2("in1", "out1") == "ok"
    assert process_landsat8("in2", "out2") == "ok"
    assert calls[0] == "s2"
    assert calls[2] == "l8"


def test_resolve_s2_backend_branches(monkeypatch):
    cfg = SIACConfig(sensor="s2", s2_data={"backend": "cdse", "cdse_access_key": "u", "cdse_secret_key": "p"})

    class _FakeCDSE:
        def __init__(self, access_key=None, secret_key=None, auth=None):
            self.access_key = access_key
            self.secret_key = secret_key
            self.auth = auth

    class _FakeGCS:
        pass

    monkeypatch.setattr("siac.io.copernicus_dataspace.CopernicusDataspaceBackend", _FakeCDSE)
    monkeypatch.setattr("siac.io.gcs_sentinel2.GCSSentinel2Backend", _FakeGCS)

    auth = SimpleNamespace(name="auth")
    cdse = _resolve_s2_backend(cfg, auth=auth)
    assert isinstance(cdse, _FakeCDSE)
    assert cdse.access_key == "u"

    cfg_gcs = SIACConfig(sensor="s2", s2_data={"backend": "gcs"})
    assert isinstance(_resolve_s2_backend(cfg_gcs), _FakeGCS)

    cfg_local = SIACConfig(sensor="s2", s2_data={"backend": "local"})
    assert _resolve_s2_backend(cfg_local) is None

    cfg_bad = SIACConfig(sensor="s2", s2_data={"backend": "local"})
    cfg_bad.s2_data.backend = "bad"
    with pytest.raises(ValueError, match="Unknown S2 backend"):
        _resolve_s2_backend(cfg_bad)


def test_coerce_date_and_query_helpers():
    cfg = SIACConfig(sensor="s2", s2_data={"max_cloud_cover": 33.0, "processing_level": "L2A"})

    assert _coerce_date(None) is None
    assert _coerce_date(date(2026, 1, 2)) == date(2026, 1, 2)
    assert _coerce_date(datetime(2026, 1, 2, 8, 0, 0)) == date(2026, 1, 2)
    assert _coerce_date("2026-01-03") == date(2026, 1, 3)

    q = S2Query.from_tile_date("T50QLD_20260102")
    q.max_cloud_cover = 100.0
    q.processing_level = "L1C"
    out_q = _apply_s2_query_defaults(q, config=cfg)
    assert out_q.max_cloud_cover == 33.0
    assert out_q.processing_level == "L2A"

    copied = _coerce_s2_query(out_q, config=cfg)
    assert copied is not out_q

    tile_q = _coerce_s2_query("T50QLD_20260102", config=cfg)
    assert tile_q.mgrs_tile == "50QLD"

    l2a_id = "S2A_MSIL2A_20260102T024121_N0511_R089_T50QLD_20260102T035433"
    assert _coerce_s2_query(l2a_id, config=cfg).processing_level == "L2A"

    cfg_l1c = SIACConfig(sensor="s2", s2_data={"max_cloud_cover": 33.0, "processing_level": "L1C"})
    l1c_id = "S2C_MSIL1C_20260102T024121_N0511_R089_T50QLD_20260102T035433.SAFE"
    assert _coerce_s2_query(l1c_id, config=cfg_l1c).processing_level == "L1C"


def test_resolve_s2_input_and_search_local_errors(tmp_path: Path):
    cfg = SIACConfig(sensor="s2", s2_data={"backend": "local"})

    with pytest.raises(DataNotFoundError, match="backend is 'local'"):
        from siac.siac import resolve_s2_input

        _ = resolve_s2_input("missing.SAFE", cfg)

    with pytest.raises(ValueError, match="does not support backend='local'"):
        _ = search_sentinel2(tile="50QLD", date_value="2026-01-02", backend="local")


def test_siac_process_uses_resolvers_and_injections(monkeypatch, tmp_path: Path):
    cfg = SIACConfig(sensor="s2")
    calls = {}

    monkeypatch.setattr("siac.siac.CredentialManager.from_config", lambda _config: "auth")
    monkeypatch.setattr("siac.siac._resolve_preprocessor", lambda _config: "pp")
    monkeypatch.setattr("siac.siac._resolve_atmo_provider", lambda _config, auth=None: ("ap", auth))
    monkeypatch.setattr("siac.siac._resolve_surface_prior_provider", lambda _config, auth=None: ("sp", auth))
    monkeypatch.setattr("siac.siac._resolve_grid_assembler", lambda: "ga")
    monkeypatch.setattr("siac.siac._resolve_solver", lambda _config: "sv")
    monkeypatch.setattr("siac.siac._resolve_corrector", lambda _config: "cr")
    monkeypatch.setattr("siac.siac._resolve_rt_model_for_pipeline", lambda _config, auth=None: ("rt", auth))

    def _fake_run_pipeline(input_path, aoi, config, **kwargs):
        calls.update(kwargs)
        return "pipeline-result"

    monkeypatch.setattr("siac.siac.run_pipeline", _fake_run_pipeline)

    out = siac_process(cfg, tmp_path / "in.SAFE")
    assert out == "pipeline-result"
    assert calls["preprocessor"] == "pp"
    assert calls["atmo_provider"] == ("ap", "auth")
    assert calls["surface_prior_provider"] == ("sp", "auth")

    # Explicit injections should bypass resolver defaults.
    out2 = siac_process(
        cfg,
        tmp_path / "in2.SAFE",
        auth="provided-auth",
        preprocessor="pp2",
        atmo_provider="ap2",
        surface_prior_provider="sp2",
        grid_assembler="ga2",
        solver="sv2",
        corrector="cr2",
        rt_model="rt2",
    )
    assert out2 == "pipeline-result"
    assert calls["preprocessor"] == "pp2"
    assert calls["atmo_provider"] == "ap2"


def test_resolve_preprocessor_and_atmo_provider_branches(monkeypatch):
    cfg = SIACConfig(sensor="s2", atmo_prior={"provider": "cams"})

    class _FakePP:
        def preprocess(self, path, aoi=None):
            return "obs"

    monkeypatch.setattr("siac.satellite.sentinel2.Sentinel2Preprocessor", _FakePP)
    pp_fn = _resolve_preprocessor(cfg)
    assert pp_fn("/tmp/scene") == "obs"

    calls = []

    class _FakeProvider:
        def __init__(self, *args, **kwargs):
            calls.append((args, kwargs))

        def get_prior(self, bounds, crs, obs_time, resolution):
            return (bounds, crs, obs_time, resolution)

    monkeypatch.setattr("siac.siac.CAMSProvider", _FakeProvider)
    atmo_fn = _resolve_atmo_provider(cfg, auth="auth")
    assert atmo_fn((0, 0, 1, 1), "EPSG:4326", datetime(2026, 1, 1), 1000.0)[3] == 1000.0
    assert calls[0][1]["download_missing"] is True
    assert calls[0][1]["cache_dir"] is None

    cfg_merra = SIACConfig(sensor="s2", atmo_prior={"provider": "merra2", "cache_dir": "/tmp/cache"})
    monkeypatch.setattr("siac.priors.atmospheric.merra2.MERRA2Provider", _FakeProvider)
    assert callable(_resolve_atmo_provider(cfg_merra))
    assert "source" in calls[1][1]


def test_resolve_surface_prior_solver_corrector_and_rt(monkeypatch):
    cfg = SIACConfig(sensor="s2", brdf={"provider": "mcd43"}, rt_model={"backend": "lut", "lut_path": "s3://bucket/lut"})

    class _FakeBRDFProvider:
        def __init__(self, cache_dir=None, source=None):
            self.cache_dir = cache_dir
            self.source = source
            self.source_bands = [SENTINEL2A_CONFIG.get_band("B02"), SENTINEL2A_CONFIG.get_band("B03")]

        def get_brdf_parameters(self, **kwargs):
            return "weights"

    class _FakeDeriver:
        def __init__(self, **kwargs):
            pass

        def compute_surface_prior(self, brdf_weights, geometry, **kwargs):
            return (brdf_weights, geometry, kwargs)

    monkeypatch.setattr("siac.priors.brdf.mcd43_earthaccess.MCD43EarthAccessProvider", _FakeBRDFProvider)
    monkeypatch.setattr("siac.siac.KernelModelDeriver", _FakeDeriver)

    surf_fn = _resolve_surface_prior_provider(cfg, auth=_earthdata_auth())
    geom = _geometry()
    obs = SimpleNamespace(
        bounds=(0, 0, 1, 1),
        crs="EPSG:4326",
        metadata={"observation_time": datetime(2026, 1, 2)},
        sensor_config=SENTINEL2A_CONFIG,
        geometry=geom,
        toa=_toa_dataset(),
        cloud_mask=xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"]),
    )
    assert getattr(surf_fn, "requires_atmo_prior", False) is False
    assert surf_fn(obs, None, "rt", 500.0)[0] == "weights"

    cfg_whittaker = SIACConfig(
        sensor="s2",
        brdf={"provider": "mcd43"},
        surface_prior={"method": "whittaker"},
        rt_model={"backend": "lut", "lut_path": "s3://bucket/lut"},
    )

    class _FakeTemporalBRDFProvider(_FakeBRDFProvider):
        def get_temporal_brdf_parameters(self, **kwargs):
            return "temporal-weights"

    class _FakeWhittakerDeriver:
        def __init__(self, **kwargs):
            pass

        def compute_surface_prior(self, brdf_weights, geometry, obs_time=None, **kwargs):
            return (brdf_weights, geometry, obs_time, kwargs)

    monkeypatch.setattr("siac.priors.brdf.mcd43_earthaccess.MCD43EarthAccessProvider", _FakeTemporalBRDFProvider)
    monkeypatch.setattr("siac.siac.BRDFWhittakerDeriver", _FakeWhittakerDeriver)

    surf_whittaker = _resolve_surface_prior_provider(cfg_whittaker, auth=_earthdata_auth())
    assert getattr(surf_whittaker, "requires_atmo_prior", False) is False
    assert surf_whittaker(obs, None, "rt", 500.0)[0] == "temporal-weights"

    cfg_monthly = SIACConfig(
        sensor="s2",
        brdf={"provider": "mcd43"},
        surface_prior={"method": "monthly_database"},
        rt_model={"backend": "lut", "lut_path": "s3://bucket/lut"},
    )

    captured: dict[str, object] = {}

    def _fake_build_monthly_database(*, observation, brdf_provider, resolution, geometry, visible_bands, query_bands):
        captured["observation"] = observation
        captured["provider"] = brdf_provider
        captured["resolution"] = resolution
        captured["geometry"] = geometry
        captured["visible_bands"] = tuple(visible_bands)
        captured["query_bands"] = tuple(query_bands)
        return "monthly-db"

    def _fake_query_surface_prior(*, observation, atmo_prior, rt_model, database, query_band_names, visible_band_names, k_neighbors):
        captured["query_observation"] = observation
        captured["atmo_prior"] = atmo_prior
        captured["rt_model"] = rt_model
        captured["database"] = database
        captured["query_band_names"] = tuple(query_band_names)
        captured["visible_band_names"] = tuple(visible_band_names)
        captured["k_neighbors"] = k_neighbors
        return SimpleNamespace(
            mask=xr.DataArray(np.ones((2, 2), dtype=bool), dims=["y", "x"]),
        )

    monkeypatch.setattr("siac.priors.brdf.mcd43_earthaccess.MCD43EarthAccessProvider", _FakeBRDFProvider)
    monkeypatch.setattr("siac.siac.build_monthly_surface_prior_database", _fake_build_monthly_database)
    monkeypatch.setattr("siac.siac.query_surface_prior_from_monthly_database", _fake_query_surface_prior)

    surf_monthly = _resolve_surface_prior_provider(cfg_monthly, auth=_earthdata_auth())
    assert getattr(surf_monthly, "requires_atmo_prior", False) is True
    monthly_out = surf_monthly(obs, _atmo_state(), "rt", 500.0)
    assert monthly_out.mask.shape == (2, 2)
    assert captured["database"] == "monthly-db"
    assert captured["atmo_prior"].aot.shape == (2, 2)

    cfg_bad_brdf = SIACConfig(sensor="s2", brdf={"provider": "mcd43"})
    cfg_bad_brdf.brdf.provider = "unknown"
    with pytest.raises(ValueError, match="Unknown BRDF provider"):
        _resolve_surface_prior_provider(cfg_bad_brdf)

    class _FakeMGSolver:
        def __init__(self, cfg):
            self.cfg = cfg

        def solve(self, toa, surface_prior, geometry, atmo_prior, rt_model, cloud_mask, bands):
            return SimpleNamespace(
                aot=_da(0.2),
                tcwv=_da(2.2),
                aot_unc=_da(0.02),
                tcwv_unc=_da(0.25),
                final_cost=0.123,
                n_iterations=5,
                success=True,
            )

    monkeypatch.setattr("siac.siac.MultiGridSolver", _FakeMGSolver)
    solve_fn = _resolve_solver(cfg)

    atmo_prior = _atmo_state()
    inputs = SimpleNamespace(
        toa=xr.DataArray(np.ones((1, 2, 2), dtype=np.float32), dims=["band", "y", "x"]),
        surface_prior="surface",
        geometry=geom,
        atmo_prior=atmo_prior,
        rt_model="rt",
        cloud_mask=xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"]),
        bands=SENTINEL2A_CONFIG.select_bands_in_range(400.0, 900.0)[:3],
    )
    solved = solve_fn(inputs, cfg)
    assert solved.converged is True

    class _FakeCorrector:
        def __init__(self, rt_model, sensor_config):
            self.rt_model = rt_model
            self.sensor_config = sensor_config

        def correct(self, toa, geometry, atmo_state, cloud_mask):
            assert atmo_state.aot.shape == toa["B02"].shape
            return "corrected"

    monkeypatch.setattr("siac.siac.AtmosphericCorrector", _FakeCorrector)
    correct_fn = _resolve_corrector(cfg)
    obs = SimpleNamespace(toa=_toa_dataset(), geometry=geom, cloud_mask=inputs.cloud_mask, sensor_config=SENTINEL2A_CONFIG)
    assert correct_fn(obs, solved, "rt") == "corrected"
    solved_mismatch = SimpleNamespace(atmo_state=_atmo_state(shape=(1, 1)))
    assert correct_fn(obs, solved_mismatch, "rt") == "corrected"

    monkeypatch.setattr("siac.siac.TwoLayerNNEmulator", lambda **_kwargs: "emu")
    cfg_emu = SIACConfig(sensor="s2", rt_model={"backend": "emulator", "emulator_dir": "/tmp/e"})
    assert _resolve_rt_model_for_pipeline(cfg_emu) == "emu"

    captured = {}

    def _fake_lut(path, interpolation_method, storage_options):
        captured["path"] = path
        captured["storage"] = dict(storage_options)
        return "lut"

    monkeypatch.setattr("siac.siac.ZarrLUTBackend", _fake_lut)
    out = _resolve_rt_model_for_pipeline(cfg, auth=_aws_auth())
    assert out == "lut"
    assert captured["path"] == "s3://bucket/lut"
    assert captured["storage"]["key"] == "AK"
