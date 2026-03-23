"""Additional unit coverage for auth helpers and provider adapters."""

from __future__ import annotations

import os
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

import siac.adapters.auth as auth_mod
from siac.errors import AuthenticationError


class _FakeResponse:
    def __init__(self, payload: object) -> None:
        self._payload = payload

    def raise_for_status(self) -> None:
        return None

    def json(self) -> object:
        return self._payload


def test_optional_dependency_and_response_object_helpers() -> None:
    json_mod = auth_mod._load_optional_dependency("json", "unused")  # noqa: SLF001
    assert json_mod.dumps({"ok": True}) == '{"ok": true}'

    with pytest.raises(AuthenticationError, match="install me"):
        auth_mod._load_optional_dependency("module_that_does_not_exist", "install me")  # noqa: SLF001

    response = _FakeResponse({"key": "value"})
    assert auth_mod._response_json_object(response, "bad") == {"key": "value"}  # noqa: SLF001

    with pytest.raises(AuthenticationError, match="bad object"):
        auth_mod._response_json_object(_FakeResponse([]), "bad object")  # noqa: SLF001


def test_cdse_helpers_cover_header_revocation_and_exchange_validation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manager = auth_mod.CredentialManager()
    manager.set_credentials("cdse", key="user", secret="secret")
    manager._oauth_tokens["cdse"] = auth_mod.OAuthToken(access_token="cached-token", expires_at=float("inf"))

    delete_calls: dict[str, Any] = {}
    monkeypatch.setattr(
        auth_mod.requests,
        "delete",
        lambda url, **kwargs: delete_calls.update({"url": url, "kwargs": kwargs}) or _FakeResponse({}),
    )

    assert manager.cdse().authorization_header() == {"Authorization": "Bearer cached-token"}
    manager.cdse().revoke_temporary_s3_credentials("ak/id")
    assert delete_calls["url"].endswith("/access_id/ak%2Fid")
    assert delete_calls["kwargs"]["headers"]["Authorization"] == "Bearer cached-token"

    monkeypatch.setattr(
        auth_mod.requests,
        "post",
        lambda *args, **kwargs: _FakeResponse({"access_token": "tok", "expires_in": "bad"}),  # noqa: ARG005
    )
    with pytest.raises(AuthenticationError, match="invalid expires_in"):
        auth_mod._cdse_token_exchange("user", "secret")


def test_cdse_wait_and_temporary_context_cover_failure_and_verify_false(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    credentials = auth_mod.CDSES3Credentials("ak", "sk")
    manager = auth_mod.CredentialManager()

    class _AlwaysFailS3FS:
        def __init__(self, **kwargs: Any) -> None:
            self.kwargs = kwargs

        def ls(self, path: str, detail: bool = False) -> list[str]:  # noqa: ARG002
            raise PermissionError(path)

    monkeypatch.setattr(auth_mod, "import_module", lambda _name: SimpleNamespace(S3FileSystem=_AlwaysFailS3FS))
    monkeypatch.setattr(auth_mod.time, "sleep", lambda _seconds: None)

    with pytest.raises(AuthenticationError, match="did not become active"):
        manager.cdse().wait_for_temporary_s3_credentials(credentials, activation_delays=(0, 1))

    seen: dict[str, Any] = {}
    monkeypatch.setattr(
        manager.cdse(),
        "create_temporary_s3_credentials",
        lambda timeout=60: seen.setdefault("created", []).append(timeout) or credentials,
    )
    monkeypatch.setattr(
        manager.cdse(),
        "revoke_temporary_s3_credentials",
        lambda access_key_id, timeout=60: seen.setdefault("revoked", []).append((access_key_id, timeout)),
    )
    monkeypatch.setattr(
        manager.cdse(),
        "wait_for_temporary_s3_credentials",
        lambda credentials, activation_delays=(0, 1, 2, 4): seen.setdefault("waited", []).append(activation_delays),  # noqa: ARG005
    )

    with manager.cdse().temporary_s3_credentials(timeout=7, verify=False) as yielded:
        assert yielded is credentials

    assert seen["created"] == [7]
    assert seen["revoked"] == [("ak", 7)]
    assert "waited" not in seen


def test_cds_aws_and_gcs_helpers_cover_empty_and_success_paths(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manager = auth_mod.CredentialManager()

    assert manager.cds().client_kwargs() == {}

    manager.set_credentials("cds", key="api-key")
    monkeypatch.setattr(auth_mod, "import_module", lambda _name: SimpleNamespace(Client=lambda **kwargs: kwargs))
    assert manager.cds().make_client(timeout=30) == {"key": "api-key", "timeout": 30}

    manager.set_credentials("aws", key="AK", secret=None)
    assert manager.aws().storage_options() == {"key": "AK"}

    manager.set_credentials("gcs", key=None, secret=None)
    assert manager.gcs().storage_options() == {}


def test_earthdata_environment_and_source_kwargs_cover_restore_and_errors(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manager = auth_mod.CredentialManager()

    monkeypatch.setenv("EARTHDATA_USERNAME", "native-user")
    monkeypatch.setenv("EARTHDATA_PASSWORD", "native-pass")
    with manager.earthdata().activate_environment():
        assert os.environ["EARTHDATA_USERNAME"] == "native-user"
        assert os.environ["EARTHDATA_PASSWORD"] == "native-pass"

    assert manager.earthdata().source_kwargs(provider="LPDAAC_ECS") == {
        "provider": "LPDAAC_ECS",
        "login_strategy": "all",
        "persist": False,
    }

    manager.set_credentials("earthdata", key="user", secret="pass")
    monkeypatch.setenv("EARTHDATA_USERNAME", "old-user")
    monkeypatch.setenv("EARTHDATA_PASSWORD", "old-pass")
    with manager.earthdata().activate_environment():
        assert os.environ["EARTHDATA_USERNAME"] == "user"
        assert os.environ["EARTHDATA_PASSWORD"] == "pass"
    assert os.environ["EARTHDATA_USERNAME"] == "old-user"
    assert os.environ["EARTHDATA_PASSWORD"] == "old-pass"

    kwargs = manager.earthdata().source_kwargs(provider=None)
    assert kwargs["earthdata_username"] == "user"
    assert kwargs["earthdata_password"] == "pass"
    assert kwargs["login_strategy"] == "environment"

    class _FakeEarthAccessSource:
        def __init__(self, **kwargs: Any) -> None:
            self.kwargs = kwargs

    monkeypatch.setattr("siac.adapters.data.earthaccess_source.EarthAccessSource", _FakeEarthAccessSource)
    source = manager.earthdata().build_earthaccess_source(provider="LPDAAC_ECS", extra="x")
    assert source.kwargs["provider"] == "LPDAAC_ECS"
    assert source.kwargs["extra"] == "x"

    broken = auth_mod.CredentialManager()
    broken.set_credentials("earthdata", key="user", secret=None)
    with pytest.raises(AuthenticationError, match="require both username"):
        broken.earthdata().source_kwargs()


def test_manager_helpers_cover_maybe_set_and_from_config_none() -> None:
    manager = auth_mod.CredentialManager.from_config(None)
    assert manager.has_credentials("cdse") is False
    assert isinstance(manager.cdse(), auth_mod.CDSEAuth)
    assert isinstance(manager.cds(), auth_mod.CDSAuth)
    assert isinstance(manager.aws(), auth_mod.AWSAuth)
    assert isinstance(manager.gcs(), auth_mod.GCSAuth)
    assert isinstance(manager.earthdata(), auth_mod.EarthdataAuth)

    auth_mod._maybe_set(manager, "custom", None, None)  # noqa: SLF001
    assert manager.has_credentials("custom") is False
    auth_mod._maybe_set(manager, "custom", "key", None)  # noqa: SLF001
    assert manager.get_credentials("custom").key == "key"

    nested = SimpleNamespace(
        cdse=SimpleNamespace(username="u", password="p"),
        cds=SimpleNamespace(api_key="api"),
        aws=SimpleNamespace(access_key_id="ak", secret_access_key="sk"),
        earthdata=SimpleNamespace(username="eu", password="ep"),
        gcs=SimpleNamespace(credentials_file=Path("/tmp/creds.json")),
    )
    loaded = auth_mod.CredentialManager()
    auth_mod._load_from_auth_config(loaded, nested)  # noqa: SLF001
    assert loaded.get_credentials("cdse").secret == "p"
    assert loaded.get_credentials("cds").key == "api"
    assert loaded.get_credentials("aws").secret == "sk"
    assert loaded.get_credentials("earthdata").key == "eu"
    assert loaded.get_credentials("gcs").key == "/tmp/creds.json"
