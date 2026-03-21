"""
Layer 5 — Custom provider injection tests.

Proves that users can replace any module with a custom callable
and the pipeline still works end-to-end.
"""

import dataclasses
from pathlib import Path

import pytest

from siac.errors import ValidationError
from siac.runtime import CorrectionResult
from siac.workflows.pipeline import run_pipeline


@pytest.mark.integration
class TestFunctionInjection:
    """Inject plain functions as module overrides."""

    def test_inject_preprocessor_function(
        self,
        mock_observation_bundle,
        mock_atmo_provider,
        mock_surface_prior_provider,
        mock_grid_assembler,
        mock_solver_fn,
        mock_corrector_fn,
        mock_rt_model,
    ):
        called = {"pp": False}

        def my_loader(path, aoi=None):
            called["pp"] = True
            return mock_observation_bundle

        result = run_pipeline(
            Path("/fake"),
            None,
            None,
            preprocessor=my_loader,
            atmo_provider=mock_atmo_provider,
            surface_prior_provider=mock_surface_prior_provider,
            grid_assembler=mock_grid_assembler,
            solver=mock_solver_fn,
            corrector=mock_corrector_fn,
            rt_model=mock_rt_model,
        )
        assert called["pp"]
        assert isinstance(result, CorrectionResult)

    def test_inject_atmo_function(
        self,
        mock_preprocessor,
        mock_atmospheric_state,
        mock_surface_prior_provider,
        mock_grid_assembler,
        mock_solver_fn,
        mock_corrector_fn,
        mock_rt_model,
    ):
        called = {"atmo": False}

        def constant_atmo(bounds, crs, obs_time, resolution):
            called["atmo"] = True
            return mock_atmospheric_state

        result = run_pipeline(
            Path("/fake"),
            None,
            None,
            preprocessor=mock_preprocessor,
            atmo_provider=constant_atmo,
            surface_prior_provider=mock_surface_prior_provider,
            grid_assembler=mock_grid_assembler,
            solver=mock_solver_fn,
            corrector=mock_corrector_fn,
            rt_model=mock_rt_model,
        )
        assert called["atmo"]
        assert isinstance(result, CorrectionResult)

    def test_inject_solver_function(
        self,
        mock_preprocessor,
        mock_atmo_provider,
        mock_surface_prior_provider,
        mock_grid_assembler,
        mock_solved_atmosphere,
        mock_corrector_fn,
        mock_rt_model,
    ):
        """Inject a passthrough solver that returns fixed AOT."""
        called = {"solver": False}

        def passthrough_solver(inputs, config):
            called["solver"] = True
            return mock_solved_atmosphere

        result = run_pipeline(
            Path("/fake"),
            None,
            None,
            preprocessor=mock_preprocessor,
            atmo_provider=mock_atmo_provider,
            surface_prior_provider=mock_surface_prior_provider,
            grid_assembler=mock_grid_assembler,
            solver=passthrough_solver,
            corrector=mock_corrector_fn,
            rt_model=mock_rt_model,
        )
        assert called["solver"]
        assert isinstance(result, CorrectionResult)

    def test_inject_corrector_function(
        self,
        mock_preprocessor,
        mock_observation_bundle,
        mock_atmo_provider,
        mock_surface_prior_provider,
        mock_grid_assembler,
        mock_solver_fn,
        mock_rt_model,
    ):
        """Inject an identity corrector (BOA = TOA)."""
        called = {"corrector": False}

        def identity_corrector(obs, solved, rt_model):
            called["corrector"] = True
            return CorrectionResult(
                boa=obs.toa,
                boa_unc=None,
                aot=solved.aot,
                tcwv=solved.tcwv,
                cloud_mask=obs.cloud_mask,
                metadata={"corrector": "identity"},
            )

        result = run_pipeline(
            Path("/fake"),
            None,
            None,
            preprocessor=mock_preprocessor,
            atmo_provider=mock_atmo_provider,
            surface_prior_provider=mock_surface_prior_provider,
            grid_assembler=mock_grid_assembler,
            solver=mock_solver_fn,
            corrector=identity_corrector,
            rt_model=mock_rt_model,
        )
        assert called["corrector"]
        assert result.boa.identical(mock_observation_bundle.toa)


@pytest.mark.integration
class TestContractViolationByCustomProvider:
    """Custom providers that violate contracts should fail with clear messages."""

    def test_inject_bad_preprocessor_missing_field(
        self,
        mock_observation_bundle,
        mock_atmo_provider,
        mock_surface_prior_provider,
        mock_grid_assembler,
        mock_solver_fn,
        mock_corrector_fn,
        mock_rt_model,
    ):
        """Custom preprocessor returns ObservationBundle missing observation_time."""
        bad_meta = {
            k: v for k, v in mock_observation_bundle.metadata.items() if k != "observation_time"
        }
        bad_obs = dataclasses.replace(mock_observation_bundle, metadata=bad_meta)

        def bad_pp(path, aoi=None):
            return bad_obs

        with pytest.raises(ValidationError, match="observation_time"):
            run_pipeline(
                Path("/fake"),
                None,
                None,
                preprocessor=bad_pp,
                atmo_provider=mock_atmo_provider,
                surface_prior_provider=mock_surface_prior_provider,
                grid_assembler=mock_grid_assembler,
                solver=mock_solver_fn,
                corrector=mock_corrector_fn,
                rt_model=mock_rt_model,
            )

    def test_inject_bad_atmo_wrong_type(
        self,
        mock_preprocessor,
        mock_surface_prior_provider,
        mock_grid_assembler,
        mock_solver_fn,
        mock_corrector_fn,
        mock_rt_model,
    ):
        """Custom atmo provider returns a dict instead of AtmosphericState."""

        def bad_atmo(bounds, crs, obs_time, resolution):
            return {"aot": 0.15}  # wrong type

        with pytest.raises((AttributeError, AssertionError)):
            run_pipeline(
                Path("/fake"),
                None,
                None,
                preprocessor=mock_preprocessor,
                atmo_provider=bad_atmo,
                surface_prior_provider=mock_surface_prior_provider,
                grid_assembler=mock_grid_assembler,
                solver=mock_solver_fn,
                corrector=mock_corrector_fn,
                rt_model=mock_rt_model,
            )
