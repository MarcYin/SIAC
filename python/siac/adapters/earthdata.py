"""Earthdata adapter helpers."""

from __future__ import annotations

from typing import TYPE_CHECKING

from siac.adapters.data.earthaccess_source import EarthAccessSource

if TYPE_CHECKING:
    from siac.adapters.auth import CredentialManager


def earthaccess_source_from_auth(
    auth: CredentialManager | None,
    *,
    provider: str | None = None,
) -> EarthAccessSource:
    """Build an ``EarthAccessSource`` using configured credentials when available."""
    if auth is None:
        return EarthAccessSource(provider=provider)
    return auth.earthdata().build_earthaccess_source(provider=provider)


__all__ = ["earthaccess_source_from_auth"]
