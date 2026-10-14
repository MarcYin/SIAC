"""Configured CAMS credentials reach cdsapi through provider assembly."""

from __future__ import annotations

import sys
from datetime import datetime
from types import SimpleNamespace
from typing import TYPE_CHECKING

import pytest

from siac.adapters.auth import CredentialManager
from siac.app._assembly_providers import _build_cams_provider
from siac.app._assembly_surface import _resolve_bestpixel_anchor_atmo_source
from siac.config import SIACConfig

if TYPE_CHECKING:
    from pathlib import Path


@pytest.mark.parametrize("factory", ["registered", "anchor"])
@pytest.mark.parametrize("credentials", ["config", "manager", "defaults"])
def test_cams_assembly_passes_credentials_to_cdsapi(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    factory: str,
    credentials: str,
) -> None:
    config = SIACConfig()
    if credentials != "defaults":
        config.auth.cds.api_key = "configured-test-key"
    auth = None
    if credentials == "manager":
        auth = CredentialManager()
        auth.set_credentials("cds", key="manager-test-key")

    captured: dict[str, object] = {}

    class FakeClient:
        def __init__(self, **kwargs):
            captured["kwargs"] = kwargs

        def retrieve(self, dataset, request):
            return self

        def download(self, path):
            captured["path"] = path

    monkeypatch.setitem(sys.modules, "cdsapi", SimpleNamespace(Client=FakeClient))
    if factory == "registered":
        config.providers.atmo.data_path = str(tmp_path)
        provider = _build_cams_provider(config, auth)
    else:
        anchor_config = SimpleNamespace(
            auth=config.auth,
            paths=SimpleNamespace(cache_root=tmp_path),
            algorithms=SimpleNamespace(
                surface_prior=SimpleNamespace(
                    bestpixel_predict_visible=True,
                    bestpixel_predict_visible_anchor_source="cams",
                    bestpixel_predict_visible_cams_data_path=str(tmp_path),
                )
            ),
        )
        prior = _resolve_bestpixel_anchor_atmo_source(anchor_config, auth)
        provider = prior.__self__

    output = provider._download_cams_file(datetime(2024, 1, 2))
    expected = {
        "config": {"key": "configured-test-key"},
        "manager": {"key": "manager-test-key"},
        "defaults": {},
    }
    assert captured["kwargs"] == expected[credentials]
    assert captured["path"] == str(output)
    if auth is not None:
        assert provider._auth is auth
