"""Unit tests for the reimplemented SIAC config system."""

from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import ValidationError

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11
    import tomli as tomllib

from siac.config import (
    DEFAULT_LUT_URL,
    AtmoProviderConfig,
    BRDFProviderConfig,
    CloudMaskAlgorithmConfig,
    ExecutionRuntimeConfig,
    OutputDefaultsConfig,
    RTAerosolSetupConfig,
    RTSetupConfig,
    RunRequest,
    S2ProviderConfig,
    SharpTransitionFilterConfig,
    SIACConfig,
    SixSAlgorithmConfig,
    SolverAlgorithmConfig,
    SurfacePriorAlgorithmConfig,
    get_jasmin_config,
    get_lut_config,
    load_system_config,
    overlay_env_secrets,
)
from siac.sixs_outputs import SIXS_DEFAULT_OUTPUT_VARIABLES


class TestSIACConfig:
    def test_default_config(self):
        config = SIACConfig()

        assert config.sensor == "auto"
        assert config.providers.atmo.kind == "cams"
        assert config.providers.brdf.kind == "mcd43"
        assert config.providers.s2.backend == "local"
        assert config.algorithms.cloud_mask.mode == "auto"
        assert config.runtime.execution.backend == "thread"
        assert config.algorithms.rt.backend == "emulator"
        assert config.paths.lut_path == DEFAULT_LUT_URL

    def test_config_from_dict(self):
        config = SIACConfig(
            sensor="s2",
            providers={
                "atmo": {"kind": "merra2"},
                "brdf": {"kind": "vnp43", "temporal_window": 8},
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
        assert config.providers.atmo.kind == "merra2"
        assert config.providers.brdf.kind == "vnp43"
        assert config.providers.brdf.temporal_window == 8
        assert config.providers.s2.backend == "gcs"
        assert config.providers.s2.max_cloud_cover == 30.0
        assert config.algorithms.cloud_mask.mode == "external_file"
        assert str(config.algorithms.cloud_mask.external_mask_path) == "/tmp/cloud.tif"
        assert config.runtime.execution.max_workers == 8
        assert config.runtime.execution.retries == 1

    def test_solver_stages_config(self):
        config = SIACConfig(
            algorithms={
                "solver": {
                    "stages": [
                        {
                            "name": "aot_pass",
                            "solve": ["aot"],
                            "fixed": ["tcwv", "tco3"],
                            "bands": ["B02", "B04"],
                        },
                        {
                            "name": "tcwv_pass",
                            "solve": "tcwv",
                            "fixed": ["aot", "tco3"],
                            "initial_state": "previous",
                        },
                    ]
                }
            }
        )

        stages = config.algorithms.solver.stages
        assert stages[0].name == "aot_pass"
        assert stages[0].solve == ("aot",)
        assert stages[0].fixed == ("tcwv", "tco3")
        assert stages[0].bands == ("B02", "B04")
        assert stages[1].solve == ("tcwv",)

    def test_surface_prior_monthly_database_resolution_policy(self):
        config = SIACConfig(
            algorithms={
                "surface_prior": {
                    "method": "monthly_database",
                    "monthly_database_resolution_policy": "aerosol",
                }
            }
        )

        assert config.algorithms.surface_prior.monthly_database_resolution_policy == "aerosol"

    def test_solver_water_mask_buffer_pixels_config(self):
        config = SIACConfig(
            algorithms={
                "solver": {
                    "water_mask_buffer_pixels": 4,
                }
            }
        )

        assert config.algorithms.solver.water_mask_buffer_pixels == 4

    def test_example_config_parses_and_matches_real_staged_setup(self):
        config = SIACConfig.from_file(Path("docs/siac-config.example.toml"))

        assert config.algorithms.surface_prior.method == "monthly_database"
        assert config.algorithms.rt.backend == "lut"
        assert config.algorithms.solver.water_mask_buffer_pixels == 32
        assert config.algorithms.solver.stages[0].bands == ("B01", "B02", "B04")
        assert config.algorithms.solver.sharp_transition_filter.enabled is True

    def test_solver_stage_rejects_solve_fixed_overlap(self):
        with pytest.raises(ValueError, match="cannot both solve and fix"):
            SIACConfig(
                algorithms={
                    "solver": {
                        "stages": [
                            {"name": "bad", "solve": ["aot"], "fixed": ["aot"]},
                        ]
                    }
                }
            )

    def test_atmo_prior_remote_url_is_preserved(self):
        config = SIACConfig(
            providers={
                "atmo": {
                    "kind": "cams",
                    "data_path": "https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/",
                }
            }
        )

        assert isinstance(config.providers.atmo.data_path, str)
        assert (
            config.providers.atmo.data_path
            == "https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/"
        )

    def test_atmo_prior_remote_s3_url_is_preserved(self):
        config = SIACConfig(
            providers={
                "atmo": {
                    "kind": "cams",
                    "data_path": "s3://eodata/CAMS/GLOBAL",
                }
            }
        )

        assert isinstance(config.providers.atmo.data_path, str)
        assert config.providers.atmo.data_path == "s3://eodata/CAMS/GLOBAL"

    def test_config_to_and_from_toml(self, tmp_path: Path):
        config = SIACConfig(
            sensor="s2",
            paths={"lut_path": "s3://bucket/lut.zarr"},
            providers={"brdf": {"kind": "vnp43", "temporal_window": 8}},
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
        assert loaded.providers.brdf.temporal_window == 8
        assert loaded.algorithms.surface_prior.spectral_mapping.k_neighbors == 7

    def test_load_system_config_expands_user_paths(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ):
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
            sensor="auto",
            providers={"s2": {"backend": "gcs"}},
            runtime={"execution": {"max_workers": 9}},
        )

        assert config.sensor == "s2"
        assert updated.sensor == "auto"
        assert updated.providers.s2.backend == "gcs"
        assert updated.runtime.execution.max_workers == 9

    def test_unsupported_processing_sensor_is_rejected(self):
        with pytest.raises(ValidationError):
            SIACConfig(sensor="l8")

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

    def test_overlay_env_secrets_expands_gcs_credentials(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ):
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


class TestAtmoProviderConfig:
    def test_valid_providers(self):
        for provider in ["cams", "merra2", "mcd19", "vnp19"]:
            config = AtmoProviderConfig(kind=provider)
            assert config.kind == provider

    @pytest.mark.parametrize("provider", ["era5", "user"])
    def test_unregistered_providers_are_rejected(self, provider: str):
        with pytest.raises(ValidationError):
            AtmoProviderConfig(kind=provider)


class TestBRDFProviderConfig:
    def test_temporal_window_bounds(self):
        config = BRDFProviderConfig(temporal_window=16)
        assert config.temporal_window == 16

        with pytest.raises(ValueError):
            BRDFProviderConfig(temporal_window=0)

        with pytest.raises(ValueError):
            BRDFProviderConfig(temporal_window=50)

    @pytest.mark.parametrize("provider", ["gee", "zarr", "user"])
    def test_unregistered_providers_are_rejected(self, provider: str):
        with pytest.raises(ValidationError):
            BRDFProviderConfig(kind=provider)


class TestSurfacePriorAlgorithmConfig:
    @pytest.mark.parametrize("method", ["kernel_model", "whittaker", "monthly_database"])
    def test_valid_methods(self, method: str):
        config = SurfacePriorAlgorithmConfig(method=method)
        assert config.method == method

    @pytest.mark.parametrize("method", ["neural", "direct"])
    def test_unregistered_methods_are_rejected(self, method: str):
        with pytest.raises(ValidationError):
            SurfacePriorAlgorithmConfig(method=method)


class TestS2ProviderConfig:
    def test_s2_backend_choices(self):
        for backend in ["cdse", "gcs", "local"]:
            cfg = S2ProviderConfig(backend=backend)
            assert cfg.backend == backend

    def test_s2_cloud_cover_bounds(self):
        cfg = S2ProviderConfig(max_cloud_cover=25.0)
        assert cfg.max_cloud_cover == 25.0
        with pytest.raises(ValueError):
            S2ProviderConfig(max_cloud_cover=-1.0)
        with pytest.raises(ValueError):
            S2ProviderConfig(max_cloud_cover=101.0)


class TestCloudMaskAlgorithmConfig:
    def test_modes(self):
        for mode in ["auto", "external_file", "user_callable", "none"]:
            cfg = CloudMaskAlgorithmConfig(mode=mode)
            assert cfg.mode == mode

    def test_target_resolution_must_be_positive(self):
        with pytest.raises(ValueError):
            CloudMaskAlgorithmConfig(target_resolution_m=0.0)

    def test_resolution_policy_choices(self):
        cfg = CloudMaskAlgorithmConfig(resolution_policy="force")
        assert cfg.resolution_policy == "force"


class TestSixSConfig:
    def test_runtime_config_defaults_cover_execution_only(self):
        cfg = SixSAlgorithmConfig()

        assert cfg.build_profile == "release"
        # Wave 18 flipped the default from "openmp" to "worker_libraries"
        # so the joint-LUT band loop can dispatch concurrently across
        # isolated 6S library copies. The OpenMP path is still selectable
        # by setting ``parallel_backend="openmp"`` explicitly.
        assert cfg.parallel_backend == "worker_libraries"
        assert cfg.output_variables == SIXS_DEFAULT_OUTPUT_VARIABLES
        assert not hasattr(cfg, "atmospheric_profile")

    def test_output_variables_normalize_and_deduplicate(self):
        cfg = SixSAlgorithmConfig(output_variables=["xap", "tgasm", "tgasm", "sutott", "rooceaw"])

        assert cfg.output_variables == ("xap", "tgasm", "sutott", "rooceaw")

    def test_surface_driven_species_knobs_are_opt_in(self):
        from siac.config.algorithms import SolverAlgorithmConfig

        default_cfg = SolverAlgorithmConfig()
        assert default_cfg.surface_driven_aerosol_species == "none"

        cfg = SolverAlgorithmConfig(surface_driven_aerosol_species="cci_climatology_exact")
        assert cfg.surface_driven_aerosol_species == "cci_climatology_exact"

        exact = SolverAlgorithmConfig(surface_driven_aerosol_species="cci_climatology_exact")
        assert exact.surface_driven_aerosol_species == "cci_climatology_exact"

    def test_rt_setup_accepts_canonical_profiles_and_components(self):
        cfg = RTSetupConfig(
            atmosphere={"profile": "auto_latitude_date", "profile_latitude": 51.5},
            aerosol={
                "profile": "user_mixture",
                "mixture": (0.55, 0.3, 0.05, 0.1),
            },
        )

        assert cfg.atmosphere is not None
        assert cfg.atmosphere.profile == "auto_latitude_date"
        assert cfg.aerosol is not None
        assert cfg.aerosol.profile == "user_mixture"
        assert cfg.aerosol.mixture == (0.55, 0.3, 0.05, 0.1)

        mie_cfg = RTAerosolSetupConfig(
            profile="user_model",
            model_path=Path("/tmp/sample.mie"),
        )

        assert mie_cfg.profile == "user_model"
        assert mie_cfg.model_path == Path("/tmp/sample.mie")

    def test_radiosonde_canonical_fields_are_accepted(self):
        levels = tuple(float(index) for index in range(34))
        cfg = RTSetupConfig(
            atmosphere={
                "profile": "user_profile",
                "radiosonde_profile": {
                    "altitude_km": levels,
                    "pressure_mb": levels,
                    "temperature_k": levels,
                    "water_g_m3": levels,
                    "ozone_g_m3": levels,
                },
            }
        )

        assert cfg.atmosphere is not None
        assert cfg.atmosphere.profile == "user_profile"
        assert cfg.atmosphere.radiosonde_profile is not None
        assert cfg.atmosphere.radiosonde_profile.altitude_km == levels

    def test_custom_aerosol_requires_mixture(self):
        with pytest.raises(ValueError, match="rt.setup.aerosol.mixture must be provided"):
            RTAerosolSetupConfig(profile="user_mixture")

    def test_user_aerosol_model_requires_path(self):
        with pytest.raises(ValueError, match="rt.setup.aerosol.model_path must be provided"):
            RTAerosolSetupConfig(profile="user_model")

    def test_siac_config_accepts_native_sixs_backend(self):
        cfg = SIACConfig(
            algorithms={
                "rt": {
                    "backend": "sixs",
                    "setup": {
                        "atmosphere": {"profile": "tropical", "columns_mode": "input_columns"},
                        "aerosol": {"profile": "continental"},
                    },
                    "sixs": {
                        "mode": "direct",
                        "output_variables": ["xap", "xbp", "xcp", "tgasm", "sutott"],
                    },
                }
            }
        )

        assert cfg.algorithms.rt.backend == "sixs"
        assert cfg.algorithms.rt.sixs.mode == "direct"
        assert cfg.algorithms.rt.setup.atmosphere is not None
        assert cfg.algorithms.rt.setup.atmosphere.profile == "tropical"
        assert cfg.algorithms.rt.sixs.output_variables == (
            "xap",
            "xbp",
            "xcp",
            "tgasm",
            "sutott",
        )

    def test_unknown_rt_backend_is_rejected(self):
        with pytest.raises(ValidationError, match="backend"):
            SIACConfig(algorithms={"rt": {"backend": "unsupported"}})


class TestExecutionRuntimeConfig:
    def test_backend_choices(self):
        for backend in ["thread", "dask"]:
            cfg = ExecutionRuntimeConfig(backend=backend)
            assert cfg.backend == backend

    def test_max_workers_and_retries_bounds(self):
        cfg = ExecutionRuntimeConfig(max_workers=2, retries=3)
        assert cfg.max_workers == 2
        assert cfg.retries == 3
        with pytest.raises(ValueError):
            ExecutionRuntimeConfig(max_workers=0)
        with pytest.raises(ValueError):
            ExecutionRuntimeConfig(retries=-1)

    def test_stage_timeout_optional_and_positive(self):
        assert ExecutionRuntimeConfig(stage_timeout_s=None).stage_timeout_s is None
        assert ExecutionRuntimeConfig(stage_timeout_s=30.0).stage_timeout_s == 30.0
        with pytest.raises(ValueError):
            ExecutionRuntimeConfig(stage_timeout_s=0.0)

    def test_stage_timeouts_values_are_positive(self):
        cfg = ExecutionRuntimeConfig(stage_timeouts={"M2.atmospheric_prior": 10.0})
        assert cfg.stage_timeouts == {"M2.atmospheric_prior": 10.0}
        with pytest.raises(ValueError, match="stage_timeouts values"):
            ExecutionRuntimeConfig(stage_timeouts={"M2.atmospheric_prior": 0.0})


class TestSolverAlgorithmConfig:
    def test_valid_bounds(self):
        config = SolverAlgorithmConfig(bounds={"aot": (0.01, 2.0), "tcwv": (0.1, 6.0)})
        assert config.aot_bounds == (0.01, 2.0)

    def test_invalid_bounds(self):
        with pytest.raises(ValueError):
            SolverAlgorithmConfig(bounds={"aot": (2.0, 0.01)})

        with pytest.raises(ValueError):
            SolverAlgorithmConfig(bounds={"tcwv": (5.0, 5.0)})


class TestSharpTransitionFilterConfig:
    def test_unknown_fields_are_rejected(self):
        with pytest.raises(ValidationError):
            SharpTransitionFilterConfig(context_window_pixels_native=31)

    def test_even_blur_kernel_rounds_up_to_odd(self):
        cfg = SharpTransitionFilterConfig(blur_kernel_pixels_native=30)
        assert cfg.blur_kernel_pixels_native == 31


class TestOutputDefaultsConfig:
    def test_valid_formats(self):
        for fmt in ["geotiff", "cog", "zarr", "netcdf"]:
            config = OutputDefaultsConfig(format=fmt)
            assert config.format == fmt

    def test_valid_compression(self):
        for comp in ["deflate", "lzw", "zstd", "none"]:
            config = OutputDefaultsConfig(compression=comp)
            assert config.compression == comp


class TestLUTConfigHelpers:
    def test_get_lut_config_preserves_s3_url(self):
        cfg = get_lut_config("s3://bucket/path/lut.zarr")
        assert cfg.algorithms.rt.backend == "lut"
        assert cfg.paths.lut_path == "s3://bucket/path/lut.zarr"
