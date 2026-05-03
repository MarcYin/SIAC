"""Sentinel-2 data-backend registry assembly."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from siac.adapters.s2_backend import build_s2_backend
from siac.app.registry import S2_BACKEND_REGISTRY

if TYPE_CHECKING:
    from siac.adapters.auth import CredentialManager


def _build_s2_backend_with_auth(
    config: Any,
    auth: CredentialManager | None,
    *,
    use_auth: bool,
) -> Any:
    return build_s2_backend(config, auth=auth if use_auth else None)


@S2_BACKEND_REGISTRY.register("cdse")
def _build_cdse_s2_backend(config: Any, auth: CredentialManager | None = None) -> Any:
    return _build_s2_backend_with_auth(config, auth, use_auth=True)


@S2_BACKEND_REGISTRY.register("gcs")
def _build_gcs_s2_backend(config: Any, auth: CredentialManager | None = None) -> Any:
    return _build_s2_backend_with_auth(config, auth, use_auth=True)


@S2_BACKEND_REGISTRY.register("local")
def _build_local_s2_backend(config: Any, auth: CredentialManager | None = None) -> Any:
    return _build_s2_backend_with_auth(config, auth, use_auth=False)


def resolve_s2_backend(config: Any, *, auth: CredentialManager | None = None) -> Any:
    return S2_BACKEND_REGISTRY.get(config.providers.s2.backend)(config, auth)
