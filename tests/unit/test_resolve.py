"""Resolver tests for the canonical assembly layer."""

from __future__ import annotations

import pytest

from siac.app.assembly import (
    resolve_atmo_provider,
    resolve_corrector,
    resolve_grid_assembler,
    resolve_preprocessor,
    resolve_rt_model_for_pipeline,
    resolve_solver,
    resolve_surface_prior_provider,
)


class TestResolveGridAssembler:
    def test_returns_callable(self):
        fn = resolve_grid_assembler()
        assert callable(fn)

    def test_is_assemble_grids(self):
        from siac.algorithms.grid.assembler import assemble_grids

        fn = resolve_grid_assembler()
        assert fn is assemble_grids


class TestResolvePreprocessor:
    def test_unknown_sensor_raises(self):
        class _FakeConfig:
            sensor = "unknown_sensor_xyz"
            cloud_mask = type("CloudMask", (), {"model_dump": lambda _self, **_kwargs: {}})()

        with pytest.raises(ValueError, match="Unknown sensor"):
            resolve_preprocessor(_FakeConfig())


class TestResolveAtmoProvider:
    def test_mcd19_provider_returns_callable(self):
        class _FakeAtmo:
            provider = "mcd19"
            cache_dir = None

        class _FakeConfig:
            atmo_prior = _FakeAtmo()

        assert callable(resolve_atmo_provider(_FakeConfig()))

    def test_vnp19_provider_returns_callable(self):
        class _FakeAtmo:
            provider = "vnp19"
            cache_dir = None

        class _FakeConfig:
            atmo_prior = _FakeAtmo()

        assert callable(resolve_atmo_provider(_FakeConfig()))

    def test_unknown_provider_raises(self):
        class _FakeAtmo:
            provider = "nonexistent_provider"

        class _FakeConfig:
            atmo_prior = _FakeAtmo()

        with pytest.raises(ValueError, match="Unknown atmo provider"):
            resolve_atmo_provider(_FakeConfig())


class TestResolveSurfacePriorProvider:
    def test_vnp43_returns_callable(self):
        class _FakeBrdf:
            provider = "vnp43"
            cache_dir = None
            temporal_window = 16

        class _FakeSurfacePrior:
            method = "kernel_model"
            psf_sigma_x = 29.75
            psf_sigma_y = 39.0
            apply_psf = True
            spectral_mapping = type(
                "SpectralMapping",
                (),
                {
                    "enabled": False,
                    "k_neighbors": 5,
                    "neighbor_estimator": "distance_weighted_mean",
                    "knn_backend": "numpy",
                    "knn_eps": 0.0,
                    "min_valid_bands": 1,
                },
            )()

        class _FakeConfig:
            brdf = _FakeBrdf()
            surface_prior = _FakeSurfacePrior()
            paths = None

        fn = resolve_surface_prior_provider(_FakeConfig())
        assert callable(fn)
        assert getattr(fn, "requires_atmo_prior", None) is False

    def test_mcd19_returns_callable(self):
        class _FakeBrdf:
            provider = "mcd19"
            cache_dir = None
            temporal_window = 16

        class _FakeSurfacePrior:
            method = "kernel_model"
            psf_sigma_x = 29.75
            psf_sigma_y = 39.0
            apply_psf = True
            spectral_mapping = type(
                "SpectralMapping",
                (),
                {
                    "enabled": False,
                    "k_neighbors": 5,
                    "neighbor_estimator": "distance_weighted_mean",
                    "knn_backend": "numpy",
                    "knn_eps": 0.0,
                    "min_valid_bands": 1,
                },
            )()

        class _FakeConfig:
            brdf = _FakeBrdf()
            surface_prior = _FakeSurfacePrior()
            paths = None

        fn = resolve_surface_prior_provider(_FakeConfig())
        assert callable(fn)
        assert getattr(fn, "requires_atmo_prior", None) is False

    def test_monthly_database_marks_atmo_dependency(self):
        class _FakeBrdf:
            provider = "mcd19"
            cache_dir = None
            temporal_window = 16

        class _FakeSurfacePrior:
            method = "monthly_database"
            psf_sigma_x = 29.75
            psf_sigma_y = 39.0
            apply_psf = True
            whittaker_lambda = 10.0
            spectral_mapping = type(
                "SpectralMapping",
                (),
                {
                    "enabled": False,
                    "k_neighbors": 5,
                    "neighbor_estimator": "distance_weighted_mean",
                    "knn_backend": "numpy",
                    "knn_eps": 0.0,
                    "min_valid_bands": 1,
                },
            )()

        class _FakeConfig:
            brdf = _FakeBrdf()
            surface_prior = _FakeSurfacePrior()
            paths = None

        fn = resolve_surface_prior_provider(_FakeConfig())
        assert callable(fn)
        assert getattr(fn, "requires_atmo_prior", None) is True


class TestResolveSolver:
    def test_returns_callable(self):
        class _FakeSolverCfg:
            aot_gamma = 10.0
            tcwv_gamma = 5.0
            aot_bounds = (0.0, 3.0)
            tcwv_bounds = (0.0, 8.0)

        class _FakeConfig:
            solver = _FakeSolverCfg()

        assert callable(resolve_solver(_FakeConfig()))


class TestResolveCorrector:
    def test_returns_callable(self):
        class _FakeConfig:
            pass

        assert callable(resolve_corrector(_FakeConfig()))


class TestResolveRTModel:
    def test_unknown_backend_raises(self):
        class _FakeRTConfig:
            backend = "unknown"
            lut_path = None

        class _FakeConfig:
            rt_model = _FakeRTConfig()
            sensor = "s2"

        with pytest.raises(ValueError, match="Cannot resolve RT model"):
            resolve_rt_model_for_pipeline(_FakeConfig())
