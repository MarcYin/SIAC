"""
Integration tests: solver -> corrector pipeline + run_pipeline orchestration.
"""

import dataclasses
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.correction.atmospheric import AtmosphericCorrector
from siac.algorithms.solver.multigrid import MultiGridConfig, MultiGridSolver
from siac.catalog import SENTINEL2A_CONFIG
from siac.domain import SensorBand
from siac.errors import ValidationError
from siac.runtime import (
    AtmosphericState,
    BRDFKernelWeights,
    CorrectionResult,
    GeometryAngles,
    SurfacePrior,
)
from siac.workflows.pipeline import run_pipeline


@pytest.mark.integration
class TestPipelineSmoke:
    """End-to-end smoke test with synthetic data."""

    @pytest.fixture
    def synthetic_scene(self):
        """Create a synthetic 32x32 scene."""
        shape = (32, 32)

        # TOA reflectance (3 bands as Dataset for corrector, DataArray for solver)
        bands_list = [
            SensorBand("B02", 490.0, 65.0, 10.0, 0),
            SensorBand("B03", 560.0, 35.0, 10.0, 1),
            SensorBand("B04", 665.0, 30.0, 10.0, 2),
        ]

        toa_vals = np.random.RandomState(42).uniform(0.05, 0.35, (3, *shape)).astype(np.float32)
        toa_da = xr.DataArray(toa_vals, dims=["band", "y", "x"])
        toa_ds = xr.Dataset({
            b.name: xr.DataArray(toa_vals[i], dims=["y", "x"])
            for i, b in enumerate(bands_list)
        })

        # Geometry
        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
        )

        # Atmospheric state prior
        atmo_prior = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        )

        # Surface prior
        brdf_weights = BRDFKernelWeights(
            f0=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            f1=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            f2=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
            f0_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            f1_unc=xr.DataArray(np.full(shape, 0.005), dims=["y", "x"]),
            f2_unc=xr.DataArray(np.full(shape, 0.002), dims=["y", "x"]),
        )

        surface_prior = SurfacePrior(
            boa=xr.DataArray(np.full(shape, 0.12), dims=["y", "x"]),
            boa_unc=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
            kernels=brdf_weights,
            mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
        )

        cloud_mask = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"])

        return {
            "toa_da": toa_da,
            "toa_ds": toa_ds,
            "geometry": geometry,
            "atmo_prior": atmo_prior,
            "surface_prior": surface_prior,
            "cloud_mask": cloud_mask,
            "bands": bands_list,
        }

    def test_solver_then_corrector(self, synthetic_scene, mock_rt_model):
        """Solver output feeds into corrector and produces valid BOA."""
        scene = synthetic_scene

        # Run solver (minimal config for speed)
        config = MultiGridConfig(n_levels=2, min_grid_size=8, max_iter_per_level=5)
        solver = MultiGridSolver(config)

        solver_result = solver.solve(
            toa=scene["toa_da"],
            surface_prior=scene["surface_prior"],
            geometry=scene["geometry"],
            atmo_prior=scene["atmo_prior"],
            rt_model=mock_rt_model,
            cloud_mask=scene["cloud_mask"],
            bands=scene["bands"],
        )

        # Basic solver result checks
        assert solver_result.aot.shape == (32, 32)
        assert solver_result.tcwv.shape == (32, 32)
        assert solver_result.success or solver_result.n_iterations > 0

        # Build updated atmospheric state from solver output
        solved_atmo = scene["atmo_prior"].with_updated_aot_tcwv(
            aot=solver_result.aot,
            tcwv=solver_result.tcwv,
        )

        # Run corrector
        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG)
        correction_result = corrector.correct(
            scene["toa_ds"], scene["geometry"], solved_atmo
        )

        # BOA should be in physically valid range
        for band_name in ["B02", "B03", "B04"]:
            boa = correction_result.boa[band_name].values
            valid = np.isfinite(boa) & (boa > 0)
            assert valid.sum() > 0, f"No valid BOA pixels for {band_name}"
            assert np.nanmax(boa[valid]) < 1.5, f"BOA too high for {band_name}"
            assert np.nanmin(boa[valid]) > -0.1, f"BOA too low for {band_name}"

    def test_corrector_standalone(self, synthetic_scene, mock_rt_model):
        """Corrector alone should produce valid BOA from synthetic inputs."""
        scene = synthetic_scene

        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG)
        result = corrector.correct(
            scene["toa_ds"], scene["geometry"], scene["atmo_prior"]
        )

        assert len(result.boa.data_vars) == 3
        for band_name in result.boa.data_vars:
            boa = result.boa[band_name].values
            valid = np.isfinite(boa)
            assert valid.any()


# =====================================================================
# Layer 4 — Pipeline orchestration (run_pipeline)
# =====================================================================

@pytest.mark.integration
class TestRunPipeline:
    """Tests for run_pipeline() happy path and call ordering."""

    def test_pipeline_happy_path(
        self,
        mock_preprocessor,
        mock_atmo_provider,
        mock_surface_prior_provider,
        mock_grid_assembler,
        mock_solver_fn,
        mock_corrector_fn,
        mock_rt_model,
    ):
        result = run_pipeline(
            input_path=Path("/fake/path"),
            aoi=None,
            config=None,
            preprocessor=mock_preprocessor,
            atmo_provider=mock_atmo_provider,
            surface_prior_provider=mock_surface_prior_provider,
            grid_assembler=mock_grid_assembler,
            solver=mock_solver_fn,
            corrector=mock_corrector_fn,
            rt_model=mock_rt_model,
        )
        assert isinstance(result, CorrectionResult)

    def test_pipeline_calls_all_modules(
        self,
        mock_observation_bundle,
        mock_atmospheric_state,
        mock_surface_prior,
        mock_solver_input_bundle,
        mock_solved_atmosphere,
        mock_rt_model,
    ):
        """Wrap each mock with a call counter."""
        call_counts = {"m1": 0, "m2": 0, "m3": 0, "m4": 0, "m5": 0, "m6": 0}

        def pp(path, aoi=None):
            call_counts["m1"] += 1
            return mock_observation_bundle

        def atmo(bounds, crs, obs_time, res):
            call_counts["m2"] += 1
            return mock_atmospheric_state

        def surf(obs, atmo_prior, rt_model, res):
            _ = (obs, atmo_prior, rt_model, res)
            call_counts["m3"] += 1
            return mock_surface_prior

        def assemble(obs, at, sp, rt, aux_res=500.0, aero_res=1000.0):
            call_counts["m4"] += 1
            return mock_solver_input_bundle

        def solve(inputs, cfg):
            call_counts["m5"] += 1
            return mock_solved_atmosphere

        def correct(obs, solved, rt):
            call_counts["m6"] += 1
            return CorrectionResult(
                boa=obs.toa, boa_unc=None, aot=solved.aot,
                tcwv=solved.tcwv, cloud_mask=obs.cloud_mask, metadata={},
            )

        run_pipeline(
            Path("/fake"), None, None,
            preprocessor=pp, atmo_provider=atmo,
            surface_prior_provider=surf, grid_assembler=assemble,
            solver=solve, corrector=correct, rt_model=mock_rt_model,
        )
        assert all(v == 1 for v in call_counts.values()), call_counts

    def test_pipeline_module_call_order(
        self,
        mock_observation_bundle,
        mock_atmospheric_state,
        mock_surface_prior,
        mock_solver_input_bundle,
        mock_solved_atmosphere,
        mock_rt_model,
    ):
        """M1 before M4; M4 before M5; M5 before M6."""
        order = []

        def pp(path, aoi=None):
            order.append("m1")
            return mock_observation_bundle

        def atmo(bounds, crs, obs_time, res):
            order.append("m2")
            return mock_atmospheric_state

        def surf(obs, atmo_prior, rt_model, res):
            _ = (obs, atmo_prior, rt_model, res)
            order.append("m3")
            return mock_surface_prior

        def assemble(obs, at, sp, rt, aux_res=500.0, aero_res=1000.0):
            order.append("m4")
            return mock_solver_input_bundle

        def solve(inputs, cfg):
            order.append("m5")
            return mock_solved_atmosphere

        def correct(obs, solved, rt):
            order.append("m6")
            return CorrectionResult(
                boa=obs.toa, boa_unc=None, aot=solved.aot,
                tcwv=solved.tcwv, cloud_mask=obs.cloud_mask, metadata={},
            )

        run_pipeline(
            Path("/fake"), None, None,
            preprocessor=pp, atmo_provider=atmo,
            surface_prior_provider=surf, grid_assembler=assemble,
            solver=solve, corrector=correct, rt_model=mock_rt_model,
        )
        assert order.index("m1") < order.index("m4")
        assert order.index("m4") < order.index("m5")
        assert order.index("m5") < order.index("m6")

    def test_pipeline_returns_correction_result(
        self,
        mock_preprocessor,
        mock_atmo_provider,
        mock_surface_prior_provider,
        mock_grid_assembler,
        mock_solver_fn,
        mock_corrector_fn,
        mock_rt_model,
    ):
        result = run_pipeline(
            Path("/fake"), None, None,
            preprocessor=mock_preprocessor,
            atmo_provider=mock_atmo_provider,
            surface_prior_provider=mock_surface_prior_provider,
            grid_assembler=mock_grid_assembler,
            solver=mock_solver_fn,
            corrector=mock_corrector_fn,
            rt_model=mock_rt_model,
        )
        assert isinstance(result, CorrectionResult)
        assert isinstance(result.boa, xr.Dataset)


@pytest.mark.integration
class TestPipelineValidation:
    """run_pipeline catches invalid module outputs early."""

    def test_pipeline_invalid_m1_missing_time(
        self,
        mock_observation_bundle,
        mock_atmo_provider,
        mock_surface_prior_provider,
        mock_grid_assembler,
        mock_solver_fn,
        mock_corrector_fn,
        mock_rt_model,
    ):
        bad_meta = {k: v for k, v in mock_observation_bundle.metadata.items()
                    if k != "observation_time"}
        bad_obs = dataclasses.replace(mock_observation_bundle, metadata=bad_meta)

        def bad_pp(path, aoi=None):
            return bad_obs

        with pytest.raises(ValidationError, match="observation_time"):
            run_pipeline(
                Path("/fake"), None, None,
                preprocessor=bad_pp,
                atmo_provider=mock_atmo_provider,
                surface_prior_provider=mock_surface_prior_provider,
                grid_assembler=mock_grid_assembler,
                solver=mock_solver_fn,
                corrector=mock_corrector_fn,
                rt_model=mock_rt_model,
            )

    def test_pipeline_invalid_m2_negative_unc(
        self,
        mock_preprocessor,
        mock_atmospheric_state,
        mock_surface_prior_provider,
        mock_grid_assembler,
        mock_solver_fn,
        mock_corrector_fn,
        mock_rt_model,
    ):
        bad_atmo = dataclasses.replace(
            mock_atmospheric_state,
            aot_unc=xr.DataArray(
                np.full(mock_atmospheric_state.aot_unc.shape, -0.1),
                dims=["y", "x"],
            ),
        )

        def bad_m2(bounds, crs, obs_time, res):
            return bad_atmo

        with pytest.raises(ValidationError, match="non-negative"):
            run_pipeline(
                Path("/fake"), None, None,
                preprocessor=mock_preprocessor,
                atmo_provider=bad_m2,
                surface_prior_provider=mock_surface_prior_provider,
                grid_assembler=mock_grid_assembler,
                solver=mock_solver_fn,
                corrector=mock_corrector_fn,
                rt_model=mock_rt_model,
            )

    def test_pipeline_m1_exception_propagates(
        self,
        mock_atmo_provider,
        mock_surface_prior_provider,
        mock_grid_assembler,
        mock_solver_fn,
        mock_corrector_fn,
        mock_rt_model,
    ):
        def bad_pp(path, aoi=None):
            raise FileNotFoundError("SAFE dir missing")

        with pytest.raises(FileNotFoundError, match="SAFE dir missing"):
            run_pipeline(
                Path("/fake"), None, None,
                preprocessor=bad_pp,
                atmo_provider=mock_atmo_provider,
                surface_prior_provider=mock_surface_prior_provider,
                grid_assembler=mock_grid_assembler,
                solver=mock_solver_fn,
                corrector=mock_corrector_fn,
                rt_model=mock_rt_model,
            )


@pytest.mark.integration
class TestConcurrency:
    """M2 and M3 should run concurrently."""

    def test_m2_m3_run_concurrently(
        self,
        mock_observation_bundle,
        mock_atmospheric_state,
        mock_surface_prior,
        mock_solver_input_bundle,
        mock_solved_atmosphere,
        mock_rt_model,
    ):
        def slow_m2(bounds, crs, obs_time, res):
            time.sleep(0.3)
            return mock_atmospheric_state

        def slow_m3(obs, atmo_prior, rt_model, res):
            _ = (obs, atmo_prior, rt_model, res)
            time.sleep(0.3)
            return mock_surface_prior

        def pp(path, aoi=None):
            return mock_observation_bundle

        def assemble(obs, at, sp, rt, aux_res=500.0, aero_res=1000.0):
            return mock_solver_input_bundle

        def solve(inputs, cfg):
            return mock_solved_atmosphere

        def correct(obs, solved, rt):
            return CorrectionResult(
                boa=obs.toa, boa_unc=None, aot=solved.aot,
                tcwv=solved.tcwv, cloud_mask=obs.cloud_mask, metadata={},
            )

        t0 = time.monotonic()
        run_pipeline(
            Path("/fake"), None, None,
            preprocessor=pp, atmo_provider=slow_m2,
            surface_prior_provider=slow_m3, grid_assembler=assemble,
            solver=solve, corrector=correct, rt_model=mock_rt_model,
        )
        elapsed = time.monotonic() - t0
        # If truly concurrent, should take ~0.3s not ~0.6s
        assert elapsed < 0.55, f"M2+M3 took {elapsed:.2f}s — expected < 0.55s (concurrent)"

    def test_lut_preload_runs_in_parallel_with_m3(
        self,
        mock_observation_bundle,
        mock_atmospheric_state,
        mock_surface_prior,
        mock_solver_input_bundle,
        mock_solved_atmosphere,
    ):
        def slow_m2(bounds, crs, obs_time, res):
            time.sleep(0.2)
            return mock_atmospheric_state

        def slow_m3(obs, atmo_prior, rt_model, res):
            _ = (obs, atmo_prior, rt_model, res)
            time.sleep(0.4)
            return mock_surface_prior

        def pp(path, aoi=None):
            return mock_observation_bundle

        def assemble(obs, at, sp, rt, aux_res=500.0, aero_res=1000.0):
            return mock_solver_input_bundle

        def solve(inputs, cfg):
            return mock_solved_atmosphere

        def correct(obs, solved, rt):
            return CorrectionResult(
                boa=obs.toa, boa_unc=None, aot=solved.aot,
                tcwv=solved.tcwv, cloud_mask=obs.cloud_mask, metadata={},
            )

        calls: dict[str, object] = {"started": False, "n_bands": 0}

        class _RTWithPreload:
            def preload_scene_subset(self, geometry, atmo_state, bands):  # noqa: ANN001
                calls["started"] = True
                calls["n_bands"] = len(bands)
                time.sleep(0.25)

        t0 = time.monotonic()
        run_pipeline(
            Path("/fake"), None, None,
            preprocessor=pp, atmo_provider=slow_m2,
            surface_prior_provider=slow_m3, grid_assembler=assemble,
            solver=solve, corrector=correct, rt_model=_RTWithPreload(),
        )
        elapsed = time.monotonic() - t0

        assert calls["started"] is True
        assert int(calls["n_bands"]) >= 1
        # With preload started right after M2, total should stay close to max(M3, preload)
        # rather than M3 + preload.
        assert elapsed < 0.60, f"Elapsed {elapsed:.2f}s suggests LUT preload did not overlap M3."

    def test_route_b_surface_prior_waits_for_m2_but_overlaps_lut_preload(
        self,
        mock_observation_bundle,
        mock_atmospheric_state,
        mock_surface_prior,
        mock_solver_input_bundle,
        mock_solved_atmosphere,
    ):
        def slow_m2(bounds, crs, obs_time, res):
            _ = (bounds, crs, obs_time, res)
            time.sleep(0.2)
            return mock_atmospheric_state

        calls: dict[str, object] = {"atmo_seen": None, "started": False}

        def route_b_surface(obs, atmo_prior, rt_model, res):
            _ = (obs, rt_model, res)
            calls["atmo_seen"] = atmo_prior
            time.sleep(0.4)
            return mock_surface_prior

        route_b_surface.requires_atmo_prior = True

        def pp(path, aoi=None):
            _ = (path, aoi)
            return mock_observation_bundle

        def assemble(obs, at, sp, rt, aux_res=500.0, aero_res=1000.0):
            _ = (obs, at, sp, rt, aux_res, aero_res)
            return mock_solver_input_bundle

        def solve(inputs, cfg):
            _ = (inputs, cfg)
            return mock_solved_atmosphere

        def correct(obs, solved, rt):
            _ = (rt,)
            return CorrectionResult(
                boa=obs.toa, boa_unc=None, aot=solved.aot,
                tcwv=solved.tcwv, cloud_mask=obs.cloud_mask, metadata={},
            )

        class _RTWithPreload:
            def preload_scene_subset(self, geometry, atmo_state, bands):  # noqa: ANN001
                _ = (geometry, atmo_state, bands)
                calls["started"] = True
                time.sleep(0.25)

        t0 = time.monotonic()
        run_pipeline(
            Path("/fake"), None, None,
            preprocessor=pp, atmo_provider=slow_m2,
            surface_prior_provider=route_b_surface, grid_assembler=assemble,
            solver=solve, corrector=correct, rt_model=_RTWithPreload(),
        )
        elapsed = time.monotonic() - t0

        assert calls["atmo_seen"] is mock_atmospheric_state
        assert calls["started"] is True
        assert elapsed < 0.72, f"Elapsed {elapsed:.2f}s suggests Route-B M3 did not overlap LUT preload."

    def test_run_pipeline_dispatches_backend(self, monkeypatch):
        calls = {}

        def _fake_thread(*args, **kwargs):  # noqa: ANN001
            calls["backend"] = "thread"
            return "thread-result"

        def _fake_dask(*args, **kwargs):  # noqa: ANN001
            calls["backend"] = "dask"
            return "dask-result"

        monkeypatch.setattr("siac.workflows.pipeline._run_pipeline_thread", _fake_thread)
        monkeypatch.setattr("siac.workflows.pipeline._run_pipeline_dask", _fake_dask)

        cfg_thread = SimpleNamespace(
            execution=SimpleNamespace(
                backend="thread",
                max_workers=2,
                retries=0,
                stage_timeout_s=None,
                dashboard=False,
                dashboard_address=None,
                performance_report_path=None,
                show_progress=False,
            )
        )
        out_thread = run_pipeline(
            Path("/fake"),
            None,
            cfg_thread,
            preprocessor=lambda _path, _aoi=None: None,
            atmo_provider=lambda *_args, **_kwargs: None,
            surface_prior_provider=lambda *_args, **_kwargs: None,
            grid_assembler=lambda *_args, **_kwargs: None,
            solver=lambda *_args, **_kwargs: None,
            corrector=lambda *_args, **_kwargs: None,
            rt_model=None,
        )
        assert out_thread == "thread-result"
        assert calls["backend"] == "thread"

        cfg_dask = SimpleNamespace(
            execution=SimpleNamespace(
                backend="dask",
                max_workers=2,
                retries=0,
                stage_timeout_s=None,
                dashboard=False,
                dashboard_address=None,
                performance_report_path=None,
                show_progress=False,
            )
        )
        out_dask = run_pipeline(
            Path("/fake"),
            None,
            cfg_dask,
            preprocessor=lambda _path, _aoi=None: None,
            atmo_provider=lambda *_args, **_kwargs: None,
            surface_prior_provider=lambda *_args, **_kwargs: None,
            grid_assembler=lambda *_args, **_kwargs: None,
            solver=lambda *_args, **_kwargs: None,
            corrector=lambda *_args, **_kwargs: None,
            rt_model=None,
        )
        assert out_dask == "dask-result"
        assert calls["backend"] == "dask"

    def test_dask_backend_missing_dependency_raises(self, monkeypatch):
        monkeypatch.setitem(sys.modules, "dask", None)
        monkeypatch.setitem(sys.modules, "dask.distributed", None)

        with pytest.raises(RuntimeError, match="dask.distributed is not installed"):
            run_pipeline(
                Path("/fake"),
                None,
                SimpleNamespace(execution=SimpleNamespace(backend="dask")),
                preprocessor=lambda _path, _aoi=None: None,
                atmo_provider=lambda *_args, **_kwargs: None,
                surface_prior_provider=lambda *_args, **_kwargs: None,
                grid_assembler=lambda *_args, **_kwargs: None,
                solver=lambda *_args, **_kwargs: None,
                corrector=lambda *_args, **_kwargs: None,
                rt_model=None,
            )
