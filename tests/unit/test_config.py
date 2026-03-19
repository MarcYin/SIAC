"""
Unit tests for SIAC configuration system.
"""

from pathlib import Path

import pytest
import yaml

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11
    import tomli as tomllib

from siac.core.config import (
    DEFAULT_LUT_URL,
    AtmoPriorConfig,
    BRDFConfig,
    CloudMaskConfig,
    ExecutionConfig,
    OutputConfig,
    RTModelConfig,
    S2DataAccessConfig,
    SIACConfig,
    SolverConfig,
    get_lut_config,
)


class TestSIACConfig:
    """Tests for main SIACConfig class."""

    def test_default_config(self):
        """Default config should be valid."""
        config = SIACConfig()

        assert config.sensor == "auto"
        assert config.atmo_prior.provider == "cams"
        assert config.brdf.provider == "mcd43"
        assert config.s2_data.backend == "local"
        assert config.cloud_mask.mode == "auto"
        assert config.execution.backend == "thread"
        assert config.rt_model.backend == "emulator"
        assert config.rt_model.lut_path == DEFAULT_LUT_URL

    def test_config_from_dict(self):
        """Config should be creatable from dict."""
        data = {
            "sensor": "s2",
            "atmo_prior": {"provider": "merra2"},
            "brdf": {"provider": "vnp43", "temporal_window": 8},
            "s2_data": {"backend": "gcs", "max_cloud_cover": 30.0},
            "cloud_mask": {
                "mode": "external_file",
                "external_mask_path": "/tmp/cloud.tif",
                "class_mapping": {2: [8, 9, 10], 3: [3], 1: [4, 5, 6, 7]},
            },
            "execution": {"backend": "thread", "max_workers": 8, "retries": 1},
        }

        config = SIACConfig(**data)

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
            sensor="s2",
            atmo_prior={
                "provider": "cams",
                "data_path": "https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/",
            },
        )

        assert isinstance(config.atmo_prior.data_path, str)
        assert config.atmo_prior.data_path == "https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/"

    def test_atmo_prior_remote_s3_url_is_preserved(self):
        config = SIACConfig(
            sensor="s2",
            atmo_prior={
                "provider": "cams",
                "data_path": "s3://eodata/CAMS/GLOBAL",
            },
        )

        assert isinstance(config.atmo_prior.data_path, str)
        assert config.atmo_prior.data_path == "s3://eodata/CAMS/GLOBAL"

    def test_config_to_yaml(self, tmp_path: Path):
        """Config should serialize to YAML."""
        config = SIACConfig(
            sensor="l8",
            rt_model=RTModelConfig(backend="lut"),
        )

        yaml_path = tmp_path / "config.yaml"
        config.to_yaml(yaml_path)

        assert yaml_path.exists()

        with yaml_path.open() as f:
            data = yaml.safe_load(f)

        assert data["inputs"]["sensor"] == "l8"
        assert data["models"]["rt_model"]["backend"] == "lut"

    def test_config_from_yaml(self, tmp_path: Path):
        """Config should load from YAML."""
        yaml_content = """
inputs:
  sensor: s2
providers:
  atmo_prior:
    provider: era5
  brdf:
    provider: mcd43
    temporal_window: 16
models:
  solver:
    aot_gamma: 15.0
"""
        yaml_path = tmp_path / "config.yaml"
        yaml_path.write_text(yaml_content)

        config = SIACConfig.from_yaml(yaml_path)

        assert config.sensor == "s2"
        assert config.atmo_prior.provider == "era5"
        assert config.brdf.temporal_window == 16
        assert config.solver.aot_gamma == 15.0

    def test_config_from_legacy_flat_yaml(self, tmp_path: Path):
        yaml_content = """
sensor: s2
atmo_prior:
  provider: era5
brdf:
  provider: mcd43
  temporal_window: 8
"""
        yaml_path = tmp_path / "legacy.yaml"
        yaml_path.write_text(yaml_content)

        config = SIACConfig.from_file(yaml_path)

        assert config.sensor == "s2"
        assert config.atmo_prior.provider == "era5"
        assert config.brdf.temporal_window == 8

    def test_config_to_and_from_toml(self, tmp_path: Path):
        config = SIACConfig(
            sensor="s2",
            paths={"spectral_library_root": "/tmp/library", "lut_path": "s3://bucket/lut.zarr"},
            surface_prior={
                "spectral_mapping": {
                    "siac_library_root": "/tmp/override-library",
                    "k_neighbors": 7,
                }
            },
            credentials={"earthdata_username": "user", "earthdata_password": "secret"},
        )

        toml_path = tmp_path / "config.toml"
        config.to_toml(toml_path)

        with toml_path.open("rb") as handle:
            data = tomllib.load(handle)

        assert data["inputs"]["sensor"] == "s2"
        assert data["paths"]["spectral_library_root"] == "/tmp/library"
        assert data["models"]["surface_prior"]["spectral_mapping"]["k_neighbors"] == 7
        assert data["auth"]["earthdata_username"] == "user"

        round_tripped = SIACConfig.from_toml(toml_path)
        assert round_tripped.sensor == "s2"
        assert str(round_tripped.paths.spectral_library_root) == "/tmp/library"
        assert round_tripped.surface_prior.spectral_mapping.k_neighbors == 7

    def test_to_dict_uses_categorized_layout(self):
        config = SIACConfig(sensor="l8")

        data = config.to_dict()

        assert set(data) == {"inputs", "paths", "runtime", "providers", "processing", "models", "auth", "output"}
        assert data["inputs"]["sensor"] == "l8"
        assert "execution" in data["runtime"]

    def test_categorized_state_redacts_auth(self, tmp_path: Path):
        yaml_path = tmp_path / "config.yaml"
        yaml_path.write_text("inputs:\n  sensor: s2\n")
        config = SIACConfig.load(
            yaml_path,
            credentials={"earthdata_username": "user", "earthdata_password": "secret"},
            s2_data={"cdse_access_key": "access", "cdse_secret_key": "secret-key"},
        )

        state = config.categorized_state()

        assert state["load"]["config_file"] == str(yaml_path)
        assert state["config"]["auth"]["earthdata_username"] == "<redacted>"
        assert state["config"]["auth"]["earthdata_password"] == "<redacted>"
        assert state["config"]["inputs"]["s2_data"]["cdse_access_key"] == "<redacted>"
        assert state["resolved"]["paths"]["lut_path_resolved"] == DEFAULT_LUT_URL

    def test_write_default_config_uses_requested_format(self, tmp_path: Path):
        written = SIACConfig.write_default_config(tmp_path / "siac", format="toml")

        assert written == tmp_path / "siac.toml"
        assert written.exists()
        loaded = SIACConfig.from_file(written)
        assert loaded.sensor == "auto"

    def test_config_with_overrides(self):
        """with_overrides should create new config."""
        config1 = SIACConfig(sensor="s2")
        config2 = config1.with_overrides(sensor="l8")

        assert config1.sensor == "s2"
        assert config2.sensor == "l8"

    def test_invalid_sensor(self):
        """Invalid sensor should raise error."""
        with pytest.raises(ValueError):
            SIACConfig(sensor="invalid")


class TestAtmoPriorConfig:
    """Tests for atmospheric prior configuration."""

    def test_valid_providers(self):
        """All valid providers should work."""
        for provider in ["cams", "merra2", "mcd19", "vnp19", "era5", "user"]:
            if provider == "user":
                # User provider needs at least one path
                config = AtmoPriorConfig(
                    provider=provider,
                    user_aot=Path("/tmp/aot.tif"),
                )
            else:
                config = AtmoPriorConfig(provider=provider)
            assert config.provider == provider

    def test_user_provider_requires_data(self):
        """User provider should require data paths."""
        with pytest.raises(ValueError):
            AtmoPriorConfig(provider="user")


class TestBRDFConfig:
    """Tests for BRDF configuration."""

    def test_temporal_window_bounds(self):
        """Temporal window should be bounded."""
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
    """Tests for solver configuration."""

    def test_valid_bounds(self):
        """Valid bounds should work."""
        config = SolverConfig(
            aot_bounds=(0.01, 2.0),
            tcwv_bounds=(0.1, 6.0),
        )
        assert config.aot_bounds == (0.01, 2.0)

    def test_invalid_bounds(self):
        """Invalid bounds (min >= max) should fail."""
        with pytest.raises(ValueError):
            SolverConfig(aot_bounds=(2.0, 0.01))

        with pytest.raises(ValueError):
            SolverConfig(tcwv_bounds=(5.0, 5.0))


class TestOutputConfig:
    """Tests for output configuration."""

    def test_valid_formats(self):
        """All valid formats should work."""
        for fmt in ["geotiff", "cog", "zarr", "netcdf"]:
            config = OutputConfig(format=fmt)
            assert config.format == fmt

    def test_valid_compression(self):
        """All valid compression options should work."""
        for comp in ["deflate", "lzw", "zstd", "none"]:
            config = OutputConfig(compression=comp)
            assert config.compression == comp


class TestLUTConfigHelpers:
    def test_get_lut_config_preserves_s3_url(self):
        cfg = get_lut_config("s3://bucket/path/lut.zarr")
        assert cfg.rt_model.backend == "lut"
        assert cfg.rt_model.lut_path == "s3://bucket/path/lut.zarr"
