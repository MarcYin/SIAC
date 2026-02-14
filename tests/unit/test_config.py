"""
Unit tests for SIAC configuration system.
"""

import pytest
import yaml
from pathlib import Path

from siac.core.config import (
    DEFAULT_LUT_URL,
    SIACConfig,
    AtmoPriorConfig,
    BRDFConfig,
    RTModelConfig,
    SolverConfig,
    OutputConfig,
    get_lut_config,
    get_default_config,
)


class TestSIACConfig:
    """Tests for main SIACConfig class."""

    def test_default_config(self):
        """Default config should be valid."""
        config = SIACConfig()

        assert config.sensor == "auto"
        assert config.atmo_prior.provider == "cams"
        assert config.brdf.provider == "mcd43"
        assert config.rt_model.backend == "emulator"
        assert config.rt_model.lut_path == DEFAULT_LUT_URL

    def test_config_from_dict(self):
        """Config should be creatable from dict."""
        data = {
            "sensor": "s2",
            "atmo_prior": {"provider": "merra2"},
            "brdf": {"provider": "vnp43", "temporal_window": 8},
        }

        config = SIACConfig(**data)

        assert config.sensor == "s2"
        assert config.atmo_prior.provider == "merra2"
        assert config.brdf.provider == "vnp43"
        assert config.brdf.temporal_window == 8

    def test_config_to_yaml(self, tmp_path: Path):
        """Config should serialize to YAML."""
        config = SIACConfig(
            sensor="l8",
            rt_model=RTModelConfig(backend="lut"),
        )

        yaml_path = tmp_path / "config.yaml"
        config.to_yaml(yaml_path)

        assert yaml_path.exists()

        with open(yaml_path) as f:
            data = yaml.safe_load(f)

        assert data["sensor"] == "l8"
        assert data["rt_model"]["backend"] == "lut"

    def test_config_from_yaml(self, tmp_path: Path):
        """Config should load from YAML."""
        yaml_content = """
sensor: s2
atmo_prior:
  provider: era5
brdf:
  provider: mcd43
  temporal_window: 16
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
        for provider in ["cams", "merra2", "mcd19", "era5", "user"]:
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
