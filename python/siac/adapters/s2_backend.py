"""Sentinel-2 data-backend assembly."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from siac.adapters.auth import CredentialManager


def build_s2_backend(config, *, auth: CredentialManager | None = None):
    """Build the configured Sentinel-2 backend adapter."""
    backend_name = config.s2_data.backend
    if backend_name == "cdse":
        from siac.adapters.data.copernicus_dataspace import CopernicusDataspaceBackend

        return CopernicusDataspaceBackend(auth=auth)
    if backend_name == "gcs":
        from siac.adapters.data.gcs_sentinel2 import GCSSentinel2Backend

        return GCSSentinel2Backend()
    if backend_name == "local":
        return None
    raise ValueError(f"Unknown S2 backend: {backend_name!r}")


__all__ = ["build_s2_backend"]
