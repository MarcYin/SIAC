"""Unit tests for the SIAC auth layer."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

import siac.adapters.auth as auth_mod
from siac.config import (
    AuthConfig,
    AWSAuthConfig,
    CDSAuthConfig,
    CDSEAuthConfig,
    EarthdataAuthConfig,
    GCSAuthConfig,
    SIACConfig,
)
from siac.errors import AuthenticationError


class _FakeResponse:
    def __init__(self, payload: object):
        self._payload = payload

    def raise_for_status(self) -> None:
        return None

    def json(self) -> object:
        return self._payload


class _FakeSession:
    """Stand-in for the shared retry-enabled :class:`requests.Session`."""

    def __init__(
        self,
        *,
        post_handler: Any | None = None,
        delete_handler: Any | None = None,
    ) -> None:
        self._post = post_handler
        self._delete = delete_handler
        self.calls: list[tuple[str, str, dict[str, Any]]] = []

    def post(self, url: str, **kwargs: Any) -> Any:
        self.calls.append(("POST", url, kwargs))
        if self._post is None:
            raise AssertionError(f"Unexpected POST to {url}")
        return self._post(url, **kwargs)

    def delete(self, url: str, **kwargs: Any) -> Any:
        self.calls.append(("DELETE", url, kwargs))
        if self._delete is None:
            raise AssertionError(f"Unexpected DELETE to {url}")
        return self._delete(url, **kwargs)


def _patch_session(monkeypatch: pytest.MonkeyPatch, session: _FakeSession) -> None:
    monkeypatch.setattr(auth_mod, "_get_session", lambda: session)


def test_from_config_loads_nested_credentials_and_builds_earthaccess_source(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class _FakeEarthAccessSource:
        def __init__(self, **kwargs: Any) -> None:
            self.kwargs = kwargs

    monkeypatch.setattr(
        "siac.adapters.data.earthaccess_source.EarthAccessSource",
        _FakeEarthAccessSource,
    )

    config = SIACConfig(
        auth=AuthConfig(
            cdse=CDSEAuthConfig(username="cdse-user", password="cdse-pass"),
            cds=CDSAuthConfig(api_key="cds-key"),
            aws=AWSAuthConfig(access_key_id="aws-key", secret_access_key="aws-secret"),
            earthdata=EarthdataAuthConfig(username="ed-user", password="ed-pass"),
            gcs=GCSAuthConfig(credentials_file=Path("/tmp/creds.json")),
        )
    )

    manager = auth_mod.CredentialManager.from_config(config)

    assert manager.has_credentials("cdse")
    assert manager.get_credentials("cds").key == "cds-key"
    assert manager.get_credentials("gcs").key == "/tmp/creds.json"

    source = manager.earthdata().build_earthaccess_source(provider="LPDAAC_ECS")
    assert isinstance(source, _FakeEarthAccessSource)
    assert source.kwargs["provider"] == "LPDAAC_ECS"
    assert source.kwargs["earthdata_username"] == "ed-user"
    assert source.kwargs["earthdata_password"] == "ed-pass"
    assert source.kwargs["login_strategy"] == "environment"
    assert source.kwargs["persist"] is False


def test_cdse_token_exchange_rejects_non_object_response(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = _FakeSession(post_handler=lambda url, **kwargs: _FakeResponse([]))  # noqa: ARG005
    _patch_session(monkeypatch, session)

    with pytest.raises(AuthenticationError, match="JSON object"):
        auth_mod._cdse_token_exchange("user", "secret")


def test_cdse_create_temporary_s3_credentials_validates_payload(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manager = auth_mod.CredentialManager()
    manager.set_credentials("cdse", key="user", secret="secret")

    def _fake_post(url: str, **kwargs: Any) -> _FakeResponse:  # noqa: ARG001
        if url == auth_mod._CDSE_TOKEN_URL:
            return _FakeResponse({"access_token": "token", "expires_in": "60"})
        if url == auth_mod._CDSE_S3_CREDENTIALS_URL:
            return _FakeResponse({"secret": "temporary-secret"})
        raise AssertionError(f"Unexpected URL: {url}")

    _patch_session(monkeypatch, _FakeSession(post_handler=_fake_post))

    with pytest.raises(AuthenticationError, match="missing access_id"):
        manager.cdse().create_temporary_s3_credentials()


def test_cdse_wait_for_temporary_s3_credentials_requires_s3fs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        auth_mod,
        "import_module",
        lambda name: (_ for _ in ()).throw(ModuleNotFoundError(name)),  # noqa: ARG005
    )

    credentials = auth_mod.CDSES3Credentials("key", "secret")
    with pytest.raises(AuthenticationError, match="s3fs is required"):
        auth_mod.CredentialManager().cdse().wait_for_temporary_s3_credentials(
            credentials,
            activation_delays=(),
        )


def test_cds_make_client_requires_cdsapi(monkeypatch: pytest.MonkeyPatch) -> None:
    manager = auth_mod.CredentialManager()
    manager.set_credentials("cds", key="api-key")

    monkeypatch.setattr(
        auth_mod,
        "import_module",
        lambda name: (
            (_ for _ in ()).throw(ModuleNotFoundError(name)) if name == "cdsapi" else object()
        ),  # noqa: ARG005
    )

    with pytest.raises(AuthenticationError, match="cdsapi is not installed"):
        manager.cds().make_client()
