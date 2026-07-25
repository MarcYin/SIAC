"""Assembly wiring for the fused aerosol prior."""

from __future__ import annotations

from typing import Any

import pytest

from siac.app import _assembly_providers as providers
from siac.app.registry import ATMO_PROVIDER_REGISTRY


class _Recorder:
    """Stand-in provider that records which registry name built it."""

    def __init__(self, name: str) -> None:
        self.source_name = name

    def get_prior(self, *args: Any, **kwargs: Any) -> Any:  # pragma: no cover - unused
        raise AssertionError("not called in wiring tests")


@pytest.fixture
def stub_registry(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        providers,
        "_build_registered_component",
        lambda registry, name, config, auth=None: _Recorder(str(name)),
    )


def _config(**atmo: Any) -> Any:
    from siac.config import SIACConfig

    return SIACConfig.model_validate({"providers": {"atmo": atmo}})


def test_single_provider_is_used_directly(stub_registry: None) -> None:
    resolved = providers.resolve_atmo_provider(_config(kind="mcd19"))

    assert resolved.__self__.source_name == "mcd19"  # type: ignore[attr-defined]


def test_fusion_wraps_the_primary_with_the_extra_sources(stub_registry: None) -> None:
    resolved = providers.resolve_atmo_provider(
        _config(kind="mcd19", fuse_aod_with=["cams"], fuse_aod_op="max")
    )

    # The validated recipe: MAIAC primary, CAMS fused in under max.
    assert resolved.__self__.source_name == "max(mcd19, cams)"  # type: ignore[attr-defined]


def test_fusion_operator_is_honoured(stub_registry: None) -> None:
    resolved = providers.resolve_atmo_provider(
        _config(kind="mcd19", fuse_aod_with=["cams"], fuse_aod_op="mean")
    )

    assert resolved.__self__.source_name.startswith("mean(")  # type: ignore[attr-defined]


def test_all_registered_atmo_kinds_are_fusable() -> None:
    # Any registered provider must be usable as a fusion source; a kind missing
    # from the registry would fail only at run time, deep inside a scene.
    from siac.config.types import AtmoProviderKind

    for kind in AtmoProviderKind:
        assert ATMO_PROVIDER_REGISTRY.get(kind.value) is not None
