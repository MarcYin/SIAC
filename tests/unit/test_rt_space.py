"""RT-space identity and the prior/solver consistency guard."""

from __future__ import annotations

import pytest

from siac.domain.rt_space import RTSpace, describe_aerosol_model
from siac.runtime.validation import ValidationError, validate_rt_space_consistency


class _Setup:
    def __init__(self, aerosol: object) -> None:
        self.aerosol = aerosol


class _Backend:
    def __init__(self, backend_name: str, rt_setup: object) -> None:
        self.backend_name = backend_name
        self.rt_setup = rt_setup


class _Prior:
    def __init__(self, rt_space: RTSpace | None) -> None:
        self.rt_space = rt_space


def _cci_setup(lon: float, lat: float, month: int) -> object:
    from siac.algorithms.rt.aerosol_species import climatology_cci_aerosol_setup

    return climatology_cci_aerosol_setup(lon, lat, month)


def test_describe_aerosol_model_handles_missing_aerosol() -> None:
    assert describe_aerosol_model(None) == "default"


def test_named_profile_descriptor_is_the_profile() -> None:
    class Aerosol:
        profile = "continental"
        distribution = None
        mixture = None
        model_path = None

    assert describe_aerosol_model(Aerosol()) == "continental"


def test_mixture_fractions_are_normalized_into_the_descriptor() -> None:
    class Aerosol:
        profile = "user_mixture"
        distribution = None
        mixture = (0.5, 0.5, 0.0, 0.0)
        model_path = None

    descriptor = describe_aerosol_model(Aerosol())
    assert descriptor.startswith("user_mixture[")
    # Equal halves of a two-component mixture normalize to 0.50 each.
    assert "0.50" in descriptor


def test_distinct_cci_mixtures_produce_distinct_spaces() -> None:
    desert = RTSpace.from_setup("sixs", _Setup(_cci_setup(39.1, 22.3, 9)))
    urban = RTSpace.from_setup("sixs", _Setup(_cci_setup(100.5, 13.7, 5)))

    assert desert.matches(desert)
    assert not desert.matches(urban)


def test_same_aerosol_different_backend_does_not_match() -> None:
    setup = _Setup(_cci_setup(39.1, 22.3, 9))
    assert not RTSpace.from_setup("sixs", setup).matches(RTSpace.from_setup("lut", setup))


def test_from_rt_model_returns_none_without_backend_name() -> None:
    assert RTSpace.from_rt_model(object()) is None


def test_consistency_guard_passes_for_matching_spaces() -> None:
    setup = _Setup(_cci_setup(39.1, 22.3, 9))
    space = RTSpace.from_setup("sixs", setup)

    validate_rt_space_consistency(_Prior(space), _Backend("sixs", setup))


def test_consistency_guard_rejects_cross_space_prior() -> None:
    prior_space = RTSpace("lut", "continental_average")
    solve_setup = _Setup(_cci_setup(39.1, 22.3, 9))

    with pytest.raises(ValidationError, match="corrected in RT space"):
        validate_rt_space_consistency(_Prior(prior_space), _Backend("sixs", solve_setup))


def test_consistency_guard_skips_priors_without_managed_space() -> None:
    # A BRDF-kernel prior from an externally corrected product carries no RT
    # identity, so it cannot be checked and must not be rejected.
    validate_rt_space_consistency(_Prior(None), _Backend("sixs", _Setup(None)))


def test_consistency_guard_skips_unidentifiable_backends() -> None:
    validate_rt_space_consistency(_Prior(RTSpace("sixs", "continental")), object())


class _SolverCfg:
    def __init__(self, mode: str) -> None:
        self.surface_driven_aerosol_species = mode


def test_scene_adaptive_species_identity_is_the_rule_not_the_mixture() -> None:
    # The solver builds the scene's CCI mixture internally, so the backend it is
    # handed still carries the base aerosol. A library prepared under the same
    # rule must be accepted even though the resolved fractions differ per scene.
    backend = _Backend("sixs", _Setup(None))  # base setup, no mixture resolved
    space = RTSpace.for_solver(backend, _SolverCfg("cci_climatology_exact"))

    assert space == RTSpace("sixs", "cci_climatology_exact")
    validate_rt_space_consistency(_Prior(space), backend, _SolverCfg("cci_climatology_exact"))


def test_scene_adaptive_library_rejected_against_a_different_species_rule() -> None:
    backend = _Backend("sixs", _Setup(None))
    prior = _Prior(RTSpace("sixs", "cci_climatology_exact"))

    with pytest.raises(ValidationError, match="corrected in RT space"):
        validate_rt_space_consistency(prior, backend, _SolverCfg("canonical_6s"))


def test_species_mode_none_falls_back_to_the_resolved_setup() -> None:
    setup = _Setup(_cci_setup(39.1, 22.3, 9))
    backend = _Backend("sixs", setup)

    assert RTSpace.for_solver(backend, _SolverCfg("none")) == RTSpace.from_rt_model(backend)
