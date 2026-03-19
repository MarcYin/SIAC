"""
Layer 7 — Config resolution tests.

Tests the ``_resolve_*`` helpers in ``siac.siac`` that convert a
``SIACConfig`` into module callables for the pipeline.

Note: ``siac.siac`` depends on ``siac.config.SIACConfig`` which
requires ``pydantic_settings``.  If that package is not installed, these
tests are skipped.
"""

import pytest

try:
    from siac.siac import (
        _resolve_atmo_provider,
        _resolve_corrector,
        _resolve_grid_assembler,
        _resolve_preprocessor,
        _resolve_rt_model_for_pipeline,
        _resolve_solver,
        _resolve_surface_prior_provider,
    )
    _HAS_SIAC = True
except ImportError:
    _HAS_SIAC = False

pytestmark = pytest.mark.skipif(
    not _HAS_SIAC,
    reason="siac.siac requires pydantic_settings",
)


class TestResolveGridAssembler:
    def test_returns_callable(self):
        fn = _resolve_grid_assembler()
        assert callable(fn)

    def test_is_assemble_grids(self):
        from siac.algorithms.grid.assembler import assemble_grids
        fn = _resolve_grid_assembler()
        assert fn is assemble_grids


class TestResolvePreprocessor:
    def test_unknown_sensor_raises(self):
        """Unknown sensor should raise ValueError."""

        class _FakeConfig:
            sensor = "unknown_sensor_xyz"

        with pytest.raises(ValueError, match="Unknown sensor"):
            _resolve_preprocessor(_FakeConfig())


class TestResolveAtmoProvider:
    def test_mcd19_provider_returns_callable(self):
        """MCD19 atmo provider branch should resolve to a callable."""

        class _FakeAtmo:
            provider = "mcd19"
            cache_dir = None

        class _FakeConfig:
            atmo_prior = _FakeAtmo()

        fn = _resolve_atmo_provider(_FakeConfig())
        assert callable(fn)

    def test_vnp19_provider_returns_callable(self):
        """VNP19 atmo provider branch should resolve to a callable."""

        class _FakeAtmo:
            provider = "vnp19"
            cache_dir = None

        class _FakeConfig:
            atmo_prior = _FakeAtmo()

        fn = _resolve_atmo_provider(_FakeConfig())
        assert callable(fn)

    def test_unknown_provider_raises(self):
        """Unknown atmo provider should raise ValueError."""

        class _FakeAtmo:
            provider = "nonexistent_provider"

        class _FakeConfig:
            atmo_prior = _FakeAtmo()

        with pytest.raises(ValueError, match="Unknown atmo provider"):
            _resolve_atmo_provider(_FakeConfig())


class TestResolveSurfacePriorProvider:
    def test_vnp43_returns_callable(self):
        """VNP43 branch should resolve to a callable."""

        class _FakeBrdf:
            provider = "vnp43"
            cache_dir = None
            temporal_window = 16

        class _FakeSurfacePrior:
            psf_sigma_x = 29.75
            psf_sigma_y = 39.0
            apply_psf = True

        class _FakeConfig:
            brdf = _FakeBrdf()
            surface_prior = _FakeSurfacePrior()

        fn = _resolve_surface_prior_provider(_FakeConfig())
        assert callable(fn)

    def test_mcd19_returns_callable(self):
        """MCD19 BRDF branch should resolve to a callable."""

        class _FakeBrdf:
            provider = "mcd19"
            cache_dir = None
            temporal_window = 16

        class _FakeSurfacePrior:
            psf_sigma_x = 29.75
            psf_sigma_y = 39.0
            apply_psf = True

        class _FakeConfig:
            brdf = _FakeBrdf()
            surface_prior = _FakeSurfacePrior()

        fn = _resolve_surface_prior_provider(_FakeConfig())
        assert callable(fn)

    def test_returns_callable(self):
        """Default resolver should return a callable."""

        class _FakeBrdf:
            cache_dir = None
            temporal_window = 16

        class _FakeSurfacePrior:
            psf_sigma_x = 29.75
            psf_sigma_y = 39.0
            apply_psf = True

        class _FakeConfig:
            brdf = _FakeBrdf()
            surface_prior = _FakeSurfacePrior()

        fn = _resolve_surface_prior_provider(_FakeConfig())
        assert callable(fn)


class TestResolveSolver:
    def test_returns_callable(self):
        """Default solver resolver should return a callable."""

        class _FakeSolverCfg:
            aot_gamma = 10.0
            tcwv_gamma = 5.0
            aot_bounds = (0.0, 3.0)
            tcwv_bounds = (0.0, 8.0)

        class _FakeConfig:
            solver = _FakeSolverCfg()

        fn = _resolve_solver(_FakeConfig())
        assert callable(fn)


class TestResolveCorrector:
    def test_returns_callable(self):
        """Default corrector resolver should return a callable."""

        class _FakeConfig:
            pass

        fn = _resolve_corrector(_FakeConfig())
        assert callable(fn)


class TestResolveRTModel:
    def test_unknown_backend_raises(self):
        """Unknown RT backend should raise ValueError."""

        class _FakeRTConfig:
            backend = "unknown"
            lut_path = None

        class _FakeConfig:
            rt_model = _FakeRTConfig()

        with pytest.raises(ValueError, match="Cannot resolve RT model"):
            _resolve_rt_model_for_pipeline(_FakeConfig())
