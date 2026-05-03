"""Resolver tests for the canonical assembly layer."""

from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

import siac.app._assembly_surface as surface_assembly_mod
import siac.app.s2_backend as s2_backend_mod
from siac.app._assembly_correction import resolve_corrector
from siac.app._assembly_io import resolve_grid_assembler
from siac.app._assembly_preprocessor import resolve_preprocessor
from siac.app._assembly_providers import resolve_atmo_provider
from siac.app._assembly_rt import resolve_rt_model_for_pipeline
from siac.app._assembly_solver import resolve_solver
from siac.app._assembly_surface import (
    _build_kernel_surface_prior,
    _build_monthly_surface_prior,
    _build_whittaker_surface_prior,
    resolve_surface_prior_provider,
)
from siac.app.s2_backend import resolve_s2_backend
from siac.errors import SensorNotSupportedError
from siac.geo.resample import resample_field_to_template as _resample_field_to_template
from siac.observability import ExecutionObserver, bind_execution_observer


def _monthly_surface_prior_test_config():
    return SimpleNamespace(
        algorithms=SimpleNamespace(
            surface_prior=SimpleNamespace(
                psf_sigma_x=29.75,
                psf_sigma_y=39.0,
                apply_psf=True,
                spectral_mapping=_spectral_mapping(enabled=True),
                monthly_database_filter=_monthly_database_filter(),
                monthly_database_resolution_policy="provider_or_coarser",
            )
        ),
        providers=SimpleNamespace(brdf=SimpleNamespace(temporal_window=16)),
        paths=None,
    )


def _spectral_mapping(*, enabled: bool = False):
    return SimpleNamespace(
        enabled=enabled,
        k_neighbors=5,
        neighbor_estimator="distance_weighted_mean",
        knn_backend="numpy",
        knn_eps=0.0,
        min_valid_bands=1,
    )


def _monthly_database_filter():
    return SimpleNamespace(
        enabled=True,
        max_prediction_uncertainty=0.05,
        max_composite_quality=0.05,
        max_source_fit_rmse=0.05,
        max_knn_feature_distance=0.05,
    )


def _surface_prior_config(*, method: str = "kernel_model", mapping_enabled: bool = False):
    return SimpleNamespace(
        method=method,
        psf_sigma_x=29.75,
        psf_sigma_y=39.0,
        apply_psf=True,
        whittaker_lambda=10.0,
        spectral_mapping=_spectral_mapping(enabled=mapping_enabled),
        monthly_database_filter=_monthly_database_filter(),
        monthly_database_resolution_policy="provider_or_coarser",
    )


def _canonical_config(
    *,
    atmo_kind: str = "mcd19",
    brdf_kind: str = "mcd43",
    surface_method: str = "kernel_model",
    s2_backend: str = "gcs",
    rt_backend: str = "unknown",
):
    return SimpleNamespace(
        sensor="s2",
        algorithms=SimpleNamespace(
            cloud_mask=SimpleNamespace(model_dump=lambda **_kwargs: {}),
            surface_prior=_surface_prior_config(method=surface_method),
            solver=SimpleNamespace(
                aot_gamma=10.0,
                tcwv_gamma=5.0,
                aot_bounds=(0.0, 3.0),
                tcwv_bounds=(0.0, 8.0),
            ),
            rt=SimpleNamespace(backend=rt_backend, lut_path=None),
        ),
        providers=SimpleNamespace(
            atmo=SimpleNamespace(kind=atmo_kind, cache_dir=None),
            brdf=SimpleNamespace(kind=brdf_kind, cache_dir=None, temporal_window=16),
            s2=SimpleNamespace(backend=s2_backend),
            monthly_composites=SimpleNamespace(kind="generated_brdf"),
        ),
        runtime=SimpleNamespace(execution=SimpleNamespace(max_workers=1)),
        paths=None,
    )


def _surface_prior_with_mask(mock_surface_prior, *, valid: bool):
    mask = xr.ones_like(mock_surface_prior.mask, dtype=bool)
    if not valid:
        mask = xr.zeros_like(mock_surface_prior.mask, dtype=bool)
    return type(mock_surface_prior)(
        boa=mock_surface_prior.boa,
        boa_unc=mock_surface_prior.boa_unc,
        kernels=mock_surface_prior.kernels,
        mask=mask,
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
        config = _canonical_config()
        config.sensor = "unknown_sensor_xyz"
        with pytest.raises(SensorNotSupportedError, match="no registered scene preprocessor"):
            resolve_preprocessor(config)


class TestResolveAtmoProvider:
    def test_mcd19_provider_returns_callable(self):
        assert callable(resolve_atmo_provider(_canonical_config(atmo_kind="mcd19")))

    def test_vnp19_provider_returns_callable(self):
        assert callable(resolve_atmo_provider(_canonical_config(atmo_kind="vnp19")))

    def test_unknown_provider_raises(self):
        with pytest.raises(ValueError, match="Unknown atmo provider"):
            resolve_atmo_provider(_canonical_config(atmo_kind="nonexistent_provider"))


class TestResolveSurfacePriorProvider:
    def test_vnp43_returns_callable(self):
        fn = resolve_surface_prior_provider(_canonical_config(brdf_kind="vnp43"))
        assert callable(fn)
        assert getattr(fn, "requires_atmo_prior", None) is False

    def test_mcd19_returns_callable(self):
        fn = resolve_surface_prior_provider(_canonical_config(brdf_kind="mcd19"))
        assert callable(fn)
        assert getattr(fn, "requires_atmo_prior", None) is False

    def test_monthly_database_marks_atmo_dependency(self):
        fn = resolve_surface_prior_provider(
            _canonical_config(brdf_kind="mcd19", surface_method="monthly_database")
        )
        assert callable(fn)
        assert getattr(fn, "requires_atmo_prior", None) is True

    def test_unknown_method_raises(self):
        with pytest.raises(ValueError, match="Unknown surface prior method: 'nope'"):
            resolve_surface_prior_provider(_canonical_config(surface_method="nope"))

    def test_whittaker_builder_forwards_observation_time(
        self,
        monkeypatch,
        mock_observation_bundle,
        mock_surface_prior,
    ):
        source_bands = list(
            mock_observation_bundle.sensor_config.select_bands_in_range(400.0, 520.0)
        )
        captured: dict[str, object] = {}

        config = _canonical_config()

        class _FakeBRDFProvider:
            def get_brdf_parameters(self, **_kwargs):
                raise AssertionError("whittaker builder should use temporal BRDF fetches")

            def get_temporal_brdf_parameters(self, **kwargs):
                captured["brdf_kwargs"] = kwargs
                return "weights"

        class _FakeWhittakerDeriver:
            def __init__(self, *args, **kwargs):
                captured["deriver_init"] = kwargs

            def compute_surface_prior(self, brdf_weights, geometry, **kwargs):
                captured["compute"] = {
                    "brdf_weights": brdf_weights,
                    "geometry": geometry,
                    **kwargs,
                }
                return mock_surface_prior

        monkeypatch.setattr(surface_assembly_mod, "BRDFWhittakerDeriver", _FakeWhittakerDeriver)

        brdf_provider = _FakeBRDFProvider()
        brdf_provider.source_bands = source_bands
        fn = _build_whittaker_surface_prior(config, brdf_provider)
        result = fn(mock_observation_bundle, None, object(), 500.0)

        assert result is mock_surface_prior
        assert captured["brdf_kwargs"]["bands"] == source_bands
        assert (
            captured["compute"]["obs_time"] == mock_observation_bundle.metadata["observation_time"]
        )
        assert captured["compute"]["source_bands"] == source_bands

    def test_kernel_builder_uses_single_time_brdf_fetch(
        self,
        monkeypatch,
        mock_observation_bundle,
        mock_surface_prior,
    ):
        source_bands = list(
            mock_observation_bundle.sensor_config.select_bands_in_range(400.0, 520.0)
        )
        captured: dict[str, object] = {}

        config = _canonical_config()

        class _FakeBRDFProvider:
            def get_brdf_parameters(self, **kwargs):
                captured["brdf_kwargs"] = kwargs
                return "weights"

            def get_temporal_brdf_parameters(self, **_kwargs):
                raise AssertionError("kernel builder should use single-time BRDF fetches")

        class _FakeKernelDeriver:
            def __init__(self, *args, **kwargs):
                captured["deriver_init"] = kwargs

            def compute_surface_prior(self, brdf_weights, geometry, **kwargs):
                captured["compute"] = {
                    "brdf_weights": brdf_weights,
                    "geometry": geometry,
                    **kwargs,
                }
                return mock_surface_prior

        monkeypatch.setattr(surface_assembly_mod, "KernelModelDeriver", _FakeKernelDeriver)

        brdf_provider = _FakeBRDFProvider()
        brdf_provider.source_bands = source_bands
        fn = _build_kernel_surface_prior(config, brdf_provider)
        result = fn(mock_observation_bundle, None, object(), 500.0)

        target_bands = list(
            mock_observation_bundle.sensor_config.select_bands_in_range(400.0, 520.0)
        )
        assert result is mock_surface_prior
        assert captured["brdf_kwargs"]["bands"] == source_bands
        assert (
            captured["brdf_kwargs"]["obs_time"]
            == mock_observation_bundle.metadata["observation_time"]
        )
        assert captured["compute"]["source_bands"] == source_bands
        assert [band.name for band in captured["compute"]["target_bands"]] == [
            band.name for band in target_bands
        ]

    def test_monthly_database_fallback_preserves_geometry_and_visible_bands(
        self,
        monkeypatch,
        mock_observation_bundle,
        mock_atmospheric_state,
        mock_surface_prior,
    ):
        captured: dict[str, object] = {}
        fallback_geometry = object()
        config = _monthly_surface_prior_test_config()
        invalid_prior = _surface_prior_with_mask(mock_surface_prior, valid=False)

        class _FakeBRDFProvider:
            source_bands = list(mock_observation_bundle.sensor_config.bands)

            def get_monthly_composites(self, _observation, _resolution):
                return SimpleNamespace(source_bands=tuple(self.source_bands), composites=())

            def get_brdf_parameters(self, **kwargs):
                captured["fallback_brdf_kwargs"] = kwargs
                return "fallback_weights"

        class _FakeKernelDeriver:
            def __init__(self, *args, **kwargs):
                captured["fallback_deriver_init"] = kwargs

            def compute_surface_prior(self, brdf_weights, geometry, **kwargs):
                captured["fallback_compute"] = {
                    "brdf_weights": brdf_weights,
                    "geometry": geometry,
                    **kwargs,
                }
                return mock_surface_prior

        def _fake_resample_geometry(observation, *, resolution):
            captured["resample"] = {"observation": observation, "resolution": resolution}
            return fallback_geometry

        def _fake_build_database(**kwargs):
            captured["database"] = kwargs
            return "db"

        def _fake_query_database(**kwargs):
            captured["query"] = kwargs
            return invalid_prior

        real_prepare = surface_assembly_mod.prepare_monthly_surface_prior_runtime
        real_query = surface_assembly_mod.query_monthly_surface_prior

        def _fake_prepare(config, monthly_provider, *, observation, resolution):
            return real_prepare(
                config,
                monthly_provider,
                observation=observation,
                resolution=resolution,
                build_database_fn=_fake_build_database,
                resample_geometry_fn=_fake_resample_geometry,
            )

        def _fake_query(observation, atmo_prior, rt_model, runtime):
            return real_query(
                observation,
                atmo_prior,
                rt_model,
                runtime,
                query_database_fn=_fake_query_database,
            )

        monkeypatch.setattr(
            surface_assembly_mod, "prepare_monthly_surface_prior_runtime", _fake_prepare
        )
        monkeypatch.setattr(surface_assembly_mod, "query_monthly_surface_prior", _fake_query)
        monkeypatch.setattr(surface_assembly_mod, "KernelModelDeriver", _FakeKernelDeriver)

        fn = _build_monthly_surface_prior(config, _FakeBRDFProvider())
        result = fn(mock_observation_bundle, mock_atmospheric_state, object(), 500.0)

        visible_bands = [
            band
            for band in mock_observation_bundle.sensor_config.bands
            if 400.0 <= band.center_wavelength < 700.0
        ]
        query_bands = list(captured["database"]["query_bands"])

        assert result is mock_surface_prior
        assert captured["database"]["geometry"] is fallback_geometry
        assert [band.name for band in captured["database"]["visible_bands"]] == [
            band.name for band in visible_bands
        ]
        assert captured["query"]["database"] == "db"
        assert captured["query"]["visible_band_names"] == tuple(band.name for band in visible_bands)
        assert captured["query"]["query_band_names"] == tuple(band.name for band in query_bands)
        assert captured["fallback_compute"]["geometry"] is fallback_geometry
        assert [band.name for band in captured["fallback_compute"]["target_bands"]] == [
            band.name for band in visible_bands
        ]
        assert (
            captured["fallback_compute"]["spectral_library"]
            is captured["database"]["spectral_library"]
        )
        assert (
            captured["fallback_compute"]["spectral_k_neighbors"]
            == captured["database"]["spectral_k_neighbors"]
        )
        assert captured["fallback_brdf_kwargs"]["target_resolution"] == 500.0
        assert (
            captured["fallback_brdf_kwargs"]["obs_time"]
            == mock_observation_bundle.metadata["observation_time"]
        )

    def test_monthly_database_success_path_skips_fallback_deriver(
        self,
        monkeypatch,
        mock_observation_bundle,
        mock_atmospheric_state,
        mock_surface_prior,
    ):
        captured: dict[str, object] = {}
        success_geometry = object()
        config = _monthly_surface_prior_test_config()
        valid_prior = _surface_prior_with_mask(mock_surface_prior, valid=True)

        class _FakeBRDFProvider:
            source_bands = list(mock_observation_bundle.sensor_config.bands)

            def get_monthly_composites(self, _observation, _resolution):
                return SimpleNamespace(source_bands=tuple(self.source_bands), composites=())

            def get_brdf_parameters(self, **_kwargs):
                raise AssertionError(
                    "monthly success path should not fall back to kernel BRDF fetches"
                )

        class _FailKernelDeriver:
            def __init__(self, *args, **kwargs):
                _ = (args, kwargs)

            def compute_surface_prior(self, *args, **kwargs):
                raise AssertionError(
                    "monthly success path should not invoke fallback kernel deriver"
                )

        def _fake_resample_geometry(observation, *, resolution):
            captured["resample"] = {"observation": observation, "resolution": resolution}
            return success_geometry

        def _fake_build_database(**kwargs):
            captured["database"] = kwargs
            return "db"

        def _fake_query_database(**kwargs):
            captured["query"] = kwargs
            return valid_prior

        real_prepare = surface_assembly_mod.prepare_monthly_surface_prior_runtime
        real_query = surface_assembly_mod.query_monthly_surface_prior

        def _fake_prepare(config, monthly_provider, *, observation, resolution):
            return real_prepare(
                config,
                monthly_provider,
                observation=observation,
                resolution=resolution,
                build_database_fn=_fake_build_database,
                resample_geometry_fn=_fake_resample_geometry,
            )

        def _fake_query(observation, atmo_prior, rt_model, runtime):
            return real_query(
                observation,
                atmo_prior,
                rt_model,
                runtime,
                query_database_fn=_fake_query_database,
            )

        monkeypatch.setattr(
            surface_assembly_mod, "prepare_monthly_surface_prior_runtime", _fake_prepare
        )
        monkeypatch.setattr(surface_assembly_mod, "query_monthly_surface_prior", _fake_query)
        monkeypatch.setattr(surface_assembly_mod, "KernelModelDeriver", _FailKernelDeriver)

        fn = _build_monthly_surface_prior(config, _FakeBRDFProvider())
        result = fn(mock_observation_bundle, mock_atmospheric_state, object(), 500.0)

        assert result is valid_prior
        assert captured["database"]["geometry"] is success_geometry
        assert captured["query"]["database"] == "db"
        assert captured["query"]["visible_band_names"] == tuple(
            band.name for band in captured["database"]["visible_bands"]
        )
        assert captured["query"]["query_band_names"] == tuple(
            band.name for band in captured["database"]["query_bands"]
        )
        assert "fallback_brdf_kwargs" not in captured
        assert "fallback_compute" not in captured

    def test_monthly_database_requires_atmo_prior(
        self,
        mock_observation_bundle,
    ):
        class _FakeBRDFProvider:
            source_bands = list(mock_observation_bundle.sensor_config.bands)

        fn = _build_monthly_surface_prior(_monthly_surface_prior_test_config(), _FakeBRDFProvider())

        with pytest.raises(ValueError, match="requires an atmospheric prior"):
            fn(mock_observation_bundle, None, object(), 500.0)

    def test_monthly_database_records_observer_substages(
        self,
        monkeypatch,
        tmp_path,
        mock_observation_bundle,
        mock_atmospheric_state,
        mock_surface_prior,
    ):
        success_geometry = object()
        config = _monthly_surface_prior_test_config()
        valid_prior = _surface_prior_with_mask(mock_surface_prior, valid=True)

        class _FakeBRDFProvider:
            source_bands = list(mock_observation_bundle.sensor_config.bands)

            def get_monthly_composites(self, _observation, _resolution):
                return SimpleNamespace(source_bands=tuple(self.source_bands), composites=())

            def get_brdf_parameters(self, **_kwargs):
                raise AssertionError("monthly success path should not fall back")

        class _FailKernelDeriver:
            def __init__(self, *args, **kwargs):
                _ = (args, kwargs)

            def compute_surface_prior(self, *args, **kwargs):
                raise AssertionError("monthly success path should not invoke fallback deriver")

        def _fake_resample_geometry(observation, *, resolution):
            _ = (observation, resolution)
            return success_geometry

        def _fake_build_database(**_kwargs):
            return "db"

        def _fake_query_database(**_kwargs):
            return valid_prior

        real_prepare = surface_assembly_mod.prepare_monthly_surface_prior_runtime
        real_query = surface_assembly_mod.query_monthly_surface_prior

        def _fake_prepare(config, monthly_provider, *, observation, resolution):
            return real_prepare(
                config,
                monthly_provider,
                observation=observation,
                resolution=resolution,
                build_database_fn=_fake_build_database,
                resample_geometry_fn=_fake_resample_geometry,
            )

        def _fake_query(observation, atmo_prior, rt_model, runtime):
            return real_query(
                observation,
                atmo_prior,
                rt_model,
                runtime,
                query_database_fn=_fake_query_database,
            )

        monkeypatch.setattr(
            surface_assembly_mod, "prepare_monthly_surface_prior_runtime", _fake_prepare
        )
        monkeypatch.setattr(surface_assembly_mod, "query_monthly_surface_prior", _fake_query)
        monkeypatch.setattr(surface_assembly_mod, "KernelModelDeriver", _FailKernelDeriver)

        report_path = tmp_path / "execution_profile.json"
        observer = ExecutionObserver(
            backend="thread",
            input_path=tmp_path / "in.SAFE",
            summary_path=report_path,
            show_progress=False,
            sample_interval_s=0.01,
            heartbeat_interval_s=0.01,
        )
        observer.start()
        try:
            with bind_execution_observer(observer):
                fn = _build_monthly_surface_prior(config, _FakeBRDFProvider())
                result = fn(mock_observation_bundle, mock_atmospheric_state, object(), 500.0)
        finally:
            observer.finish("success")

        assert result is valid_prior
        summary = json.loads(report_path.read_text(encoding="utf-8"))
        assert "M3.monthly_database.prepare_runtime" in summary["stages"]
        assert "M3.monthly_database.query" in summary["stages"]


class TestResolveSolver:
    def test_returns_callable(self):
        assert callable(resolve_solver(_canonical_config()))


class TestResolveCorrector:
    def test_returns_callable(self):
        class _FakeConfig:
            pass

        assert callable(resolve_corrector(_FakeConfig()))

    def test_resample_field_interpolates_same_shape_coordinate_mismatch(self):
        field = xr.DataArray(
            np.array([[0.0, 2.0], [10.0, 12.0]], dtype=np.float32),
            dims=["y", "x"],
            coords={"y": [1.0, 0.0], "x": [0.0, 2.0]},
        )
        template = xr.DataArray(
            np.zeros((2, 2), dtype=np.float32),
            dims=["y", "x"],
            coords={"y": [1.0, 0.0], "x": [0.0, 1.0]},
        )

        out = _resample_field_to_template(field, template)

        assert out.coords["x"].identical(template.coords["x"])
        assert float(out.sel(y=1.0, x=1.0)) == pytest.approx(1.0)

    def test_resample_field_interpolates_non_xy_coordinate_mismatch(self):
        field = xr.DataArray(
            np.array([[0.0, 2.0], [10.0, 12.0]], dtype=np.float32),
            dims=["latitude", "longitude"],
            coords={"latitude": [1.0, 0.0], "longitude": [0.0, 2.0]},
        )
        template = xr.DataArray(
            np.zeros((2, 2), dtype=np.float32),
            dims=["latitude", "longitude"],
            coords={"latitude": [1.0, 0.0], "longitude": [0.0, 1.0]},
        )

        out = _resample_field_to_template(field, template)

        assert out.coords["longitude"].identical(template.coords["longitude"])
        assert float(out.sel(latitude=1.0, longitude=1.0)) == pytest.approx(1.0)


class TestResolveRTModel:
    def test_unknown_backend_raises(self):
        with pytest.raises(ValueError, match="Cannot resolve RT model"):
            resolve_rt_model_for_pipeline(_canonical_config(rt_backend="unknown"))


class TestResolveS2Backend:
    def test_gcs_backend_forwards_auth(self, monkeypatch):
        calls: list[object | None] = []

        def _fake_build_s2_backend(config, auth=None):
            _ = config
            calls.append(auth)
            return "backend"

        auth = object()
        monkeypatch.setattr(s2_backend_mod, "build_s2_backend", _fake_build_s2_backend)

        resolved = resolve_s2_backend(_canonical_config(s2_backend="gcs"), auth=auth)

        assert resolved == "backend"
        assert calls == [auth]

    def test_local_backend_ignores_auth(self, monkeypatch):
        calls: list[object | None] = []

        def _fake_build_s2_backend(config, auth=None):
            _ = config
            calls.append(auth)
            return "backend"

        monkeypatch.setattr(s2_backend_mod, "build_s2_backend", _fake_build_s2_backend)

        resolved = resolve_s2_backend(_canonical_config(s2_backend="local"), auth=object())

        assert resolved == "backend"
        assert calls == [None]

    def test_unknown_backend_raises(self):
        with pytest.raises(ValueError, match="Unknown S2 backend: 'nope'"):
            resolve_s2_backend(_canonical_config(s2_backend="nope"))
