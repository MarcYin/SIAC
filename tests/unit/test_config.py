"""Unit tests for the reimplemented SIAC config system."""

from __future__ import annotations

from pathlib import Path

import pytest

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11
    import tomli as tomllib

from siac.config import (
    DEFAULT_LUT_URL,
    AtmoPriorConfig,
    BRDFConfig,
    CloudMaskConfig,
    ExecutionConfig,
    OutputConfig,
    RunRequest,
    S2DataAccessConfig,
    SIACConfig,
    SolverConfig,
    get_jasmin_config,
    get_lut_config,
    load_system_config,
    overlay_env_secrets,
)
from siac.config.schema import SharpTransitionFilterConfig


class TestSIACConfig:
    def test_default_config(self):
        config = SIACConfig()

        assert config.sensor == "auto"
        assert config.atmo_prior.provider == "cams"
        assert config.brdf.provider == "mcd43"
        assert config.s2_data.backend == "local"
        assert config.cloud_mask.mode == "auto"
        assert config.execution.backend == "thread"
        assert config.rt_model.backend == "emulator"
        assert config.paths.lut_path == DEFAULT_LUT_URL

    def test_config_from_dict(self):
        config = SIACConfig(
            sensor="s2",
            providers={
                "atmo": {"provider": "merra2"},
                "brdf": {"provider": "vnp43", "temporal_window": 8},
                "s2": {"backend": "gcs", "max_cloud_cover": 30.0},
            },
            algorithms={
                "cloud_mask": {
                    "mode": "external_file",
                    "external_mask_path": "/tmp/cloud.tif",
                    "class_mapping": {2: [8, 9, 10], 3: [3], 1: [4, 5, 6, 7]},
                }
            },
            runtime={"execution": {"backend": "thread", "max_workers": 8, "retries": 1}},
        )

        assert config.sensor == "s2"
        assert config.atmo_prior.provider == "merra2"
        assert config.brdf.provider == "vnp43"
        assert config.brdf.temporal_window == 8
        assert config.s2_data.backend == "gcs"
        assert config.s2_data.max_cloud_cover == 30.0
        assert config.cloud_mask.mode == "external_file"
        assert str(config.cloud_mask.external_mask_path) == "/tmp/cloud.tif"
        assert config.execution.max_workers == 8
        assert config.execution.retries == 1

    def test_atmo_prior_remote_url_is_preserved(self):
        config = SIACConfig(
            providers={
                "atmo": {
                    "provider": "cams",
                    "data_path": "https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/",
                }
            }
        )

        assert isinstance(config.atmo_prior.data_path, str)
        assert config.atmo_prior.data_path == "https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/"

    def test_atmo_prior_remote_s3_url_is_preserved(self):
        config = SIACConfig(
            providers={
                "atmo": {
                    "provider": "cams",
                    "data_path": "s3://eodata/CAMS/GLOBAL",
                }
            }
        )

        assert isinstance(config.atmo_prior.data_path, str)
        assert config.atmo_prior.data_path == "s3://eodata/CAMS/GLOBAL"

    def test_config_to_and_from_toml(self, tmp_path: Path):
        config = SIACConfig(
            sensor="l8",
            paths={"lut_path": "s3://bucket/lut.zarr"},
            providers={"brdf": {"provider": "vnp43", "temporal_window": 8}},
            algorithms={
                "surface_prior": {
                    "spectral_mapping": {
                        "k_neighbors": 7,
                    }
                },
                "rt": {"backend": "lut"},
            },
            auth={"earthdata": {"username": "user", "password": "secret"}},
        )

        toml_path = tmp_path / "config.toml"
        config.to_toml(toml_path)

        with toml_path.open("rb") as handle:
            data = tomllib.load(handle)

        assert "sensor" not in data
        assert data["providers"]["brdf"]["kind"] == "vnp43"
        assert data["algorithms"]["surface_prior"]["spectral_mapping"]["k_neighbors"] == 7

        loaded = SIACConfig.from_file(toml_path)
        assert loaded.sensor == "auto"
        assert loaded.brdf.temporal_window == 8
        assert loaded.surface_prior.spectral_mapping.k_neighbors == 7

    def test_load_system_config_expands_user_paths(self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
        monkeypatch.setenv("HOME", str(tmp_path))
        config_path = tmp_path / "config.toml"
        config_path.write_text('paths = { lut_path = "s3://bucket/lut.zarr" }\n', encoding="utf-8")

        loaded = load_system_config(Path("~/config.toml"))

        assert loaded.paths.lut_path == "s3://bucket/lut.zarr"

    def test_from_yaml_is_rejected(self, tmp_path: Path):
        yaml_path = tmp_path / "config.yaml"
        yaml_path.write_text("sensor: s2\n")
        with pytest.raises(ValueError, match="TOML"):
            SIACConfig.from_yaml(yaml_path)

    def test_with_overrides(self):
        config = SIACConfig(sensor="s2")
        updated = config.with_overrides(
            sensor="l8",
            providers={"s2": {"backend": "gcs"}},
            runtime={"execution": {"max_workers": 9}},
        )

        assert config.sensor == "s2"
        assert updated.sensor == "l8"
        assert updated.s2_data.backend == "gcs"
        assert updated.execution.max_workers == 9

    def test_resolve_uses_python_only_sensor_and_aoi_defaults(self):
        config = SIACConfig(sensor="s2", aoi="/tmp/aoi.geojson")

        resolved = config.resolve(RunRequest(input_path="/tmp/input.SAFE"))

        assert resolved.sensor == "s2"
        assert resolved.aoi == "/tmp/aoi.geojson"

    def test_resolve_uses_default_output_dir_when_run_output_missing(self, tmp_path: Path):
        config = SIACConfig(output={"defaults": {"output_dir": tmp_path / "products"}})

        resolved = config.resolve(RunRequest(input_path="/tmp/input.SAFE"))

        assert resolved.run.output_path == tmp_path / "products"

    def test_snapshot_redacts_auth(self):
        config = SIACConfig(
            auth={
                "earthdata": {"username": "user", "password": "secret"},
                "cds": {"api_key": "uid:key"},
                "gcs": {"credentials_file": "/tmp/creds.json"},
            }
        )

        snapshot = config.snapshot()

        assert snapshot["config"]["auth"]["earthdata"]["username"] == "<redacted>"
        assert snapshot["config"]["auth"]["earthdata"]["password"] == "<redacted>"
        assert snapshot["config"]["auth"]["cds"]["api_key"] == "<redacted>"
        assert snapshot["config"]["auth"]["gcs"]["credentials_file"] == "<redacted>"

    def test_overlay_env_secrets_expands_gcs_credentials(self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
        monkeypatch.setenv("HOME", str(tmp_path))
        monkeypatch.setenv("GOOGLE_APPLICATION_CREDENTIALS", "~/creds.json")
        config = SIACConfig()

        updated = overlay_env_secrets(config)

        assert updated.auth.gcs.credentials_file == tmp_path / "creds.json"

    def test_get_jasmin_config_uses_nested_models(self):
        cfg = get_jasmin_config()

        assert cfg.runtime.n_jobs == 8
        assert cfg.providers.atmo.kind == "cams"
        assert cfg.providers.brdf.kind == "mcd43"
        assert cfg.paths.dem.startswith("/vsicurl/")

    def test_write_default_config(self, tmp_path: Path):
        written = SIACConfig.write_default_config(tmp_path / "siac.toml")

        assert written == tmp_path / "siac.toml"
        loaded = SIACConfig.from_file(written)
        assert loaded.paths.lut_path == DEFAULT_LUT_URL


class TestAtmoPriorConfig:
    def test_valid_providers(self):
        for provider in ["cams", "merra2", "mcd19", "vnp19", "era5", "user"]:
            if provider == "user":
                config = AtmoPriorConfig(provider=provider, user_aot=Path("/tmp/aot.tif"))
            else:
                config = AtmoPriorConfig(provider=provider)
            assert config.provider == provider

    def test_user_provider_requires_data(self):
        with pytest.raises(ValueError):
            AtmoPriorConfig(provider="user")


class TestBRDFConfig:
    def test_temporal_window_bounds(self):
        config = BRDFConfig(temporal_window=16)
        assert config.temporal_window == 16

        with pytest.raises(ValueError):
            BRDFConfig(temporal_window=0)

        with pytest.raises(ValueError):
            BRDFConfig(temporal_window=50)


class TestS2DataAccessConfig:
    def test_s2_backend_choices(self):
        for backend in ["cdse", "gcs", "local"]:
            cfg = S2DataAccessConfig(backend=backend)
            assert cfg.backend == backend

    def test_s2_cloud_cover_bounds(self):
        cfg = S2DataAccessConfig(max_cloud_cover=25.0)
        assert cfg.max_cloud_cover == 25.0
        with pytest.raises(ValueError):
            S2DataAccessConfig(max_cloud_cover=-1.0)
        with pytest.raises(ValueError):
            S2DataAccessConfig(max_cloud_cover=101.0)


class TestCloudMaskConfig:
    def test_modes(self):
        for mode in ["auto", "external_file", "user_callable", "none"]:
            cfg = CloudMaskConfig(mode=mode)
            assert cfg.mode == mode

    def test_target_resolution_must_be_positive(self):
        with pytest.raises(ValueError):
            CloudMaskConfig(target_resolution_m=0.0)

    def test_resolution_policy_choices(self):
        cfg = CloudMaskConfig(resolution_policy="force")
        assert cfg.resolution_policy == "force"


class TestExecutionConfig:
    def test_backend_choices(self):
        for backend in ["thread", "dask"]:
            cfg = ExecutionConfig(backend=backend)
            assert cfg.backend == backend

    def test_max_workers_and_retries_bounds(self):
        cfg = ExecutionConfig(max_workers=2, retries=3)
        assert cfg.max_workers == 2
        assert cfg.retries == 3
        with pytest.raises(ValueError):
            ExecutionConfig(max_workers=0)
        with pytest.raises(ValueError):
            ExecutionConfig(retries=-1)

    def test_stage_timeout_optional_and_positive(self):
        assert ExecutionConfig(stage_timeout_s=None).stage_timeout_s is None
        assert ExecutionConfig(stage_timeout_s=30.0).stage_timeout_s == 30.0
        with pytest.raises(ValueError):
            ExecutionConfig(stage_timeout_s=0.0)


class TestSolverConfig:
    def test_valid_bounds(self):
        config = SolverConfig(bounds={"aot": (0.01, 2.0), "tcwv": (0.1, 6.0)})
        assert config.aot_bounds == (0.01, 2.0)

    def test_invalid_bounds(self):
        with pytest.raises(ValueError):
            SolverConfig(bounds={"aot": (2.0, 0.01)})

        with pytest.raises(ValueError):
            SolverConfig(bounds={"tcwv": (5.0, 5.0)})


class TestSharpTransitionFilterConfig:
    def test_legacy_fields_are_normalized_to_cv2_detector(self):
        cfg = SharpTransitionFilterConfig(
            enabled=True,
            context_window_pixels_native=31,
            road_std_z_threshold_native=2.0,
            road_coherence_threshold_native=0.95,
            point_range_z_threshold_native=3.0,
            point_outlier_fraction_max_native=0.20,
            outlier_sigma_native=2.5,
            dilation_pixels=5,
        )

        assert cfg.blur_kernel_pixels_native == 31
        assert cfg.residual_threshold_uint8 == 12
        dumped = cfg.model_dump()
        assert "context_window_pixels_native" not in dumped
        assert "road_std_z_threshold_native" not in dumped

    def test_even_blur_kernel_rounds_up_to_odd(self):
        cfg = SharpTransitionFilterConfig(blur_kernel_pixels_native=30)
        assert cfg.blur_kernel_pixels_native == 31


class TestOutputConfig:
    def test_valid_formats(self):
        for fmt in ["geotiff", "cog", "zarr", "netcdf"]:
            config = OutputConfig(format=fmt)
            assert config.format == fmt

    def test_valid_compression(self):
        for comp in ["deflate", "lzw", "zstd", "none"]:
            config = OutputConfig(compression=comp)
            assert config.compression == comp


class TestLUTConfigHelpers:
    def test_get_lut_config_preserves_s3_url(self):
        cfg = get_lut_config("s3://bucket/path/lut.zarr")
        assert cfg.rt_model.backend == "lut"
        assert cfg.paths.lut_path == "s3://bucket/path/lut.zarr"
