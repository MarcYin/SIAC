"""Tests for the centralized authentication manager."""

from __future__ import annotations

import sys
import threading
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from siac.adapters.auth import (
    AWSAuth,
    CDSEAuth,
    CredentialManager,
    CredentialSpec,
    EarthdataAuth,
    GCSAuth,
    OAuthToken,
    _cdse_token_exchange,
)
from siac.config import SIACConfig
from siac.errors import AuthenticationError

# ── 1. set/get roundtrip ────────────────────────────────────────────


class TestCredentialSetGet:
    def test_roundtrip(self):
        mgr = CredentialManager()
        mgr.set_credentials("cdse", key="user", secret="pass")
        cred = mgr.get_credentials("cdse")
        assert cred.key == "user"
        assert cred.secret == "pass"

    def test_missing_raises(self):
        mgr = CredentialManager()
        with pytest.raises(AuthenticationError, match="no_such_provider"):
            mgr.get_credentials("no_such_provider")

    def test_has_credentials_false(self):
        mgr = CredentialManager()
        assert mgr.has_credentials("cdse") is False

    def test_has_credentials_true(self):
        mgr = CredentialManager()
        mgr.set_credentials("cdse", key="u")
        assert mgr.has_credentials("cdse") is True


# ── 2. from_config with auth config ──────────────────────────────────


class TestFromConfig:
    def test_reads_credential_config_fields(self):
        """Config fields are loaded into the manager."""
        config = SIACConfig(
            auth={
                "cdse": {"username": "cdse_u", "password": "cdse_p"},
                "cds": {"api_key": "cds_k"},
                "aws": {"access_key_id": "aws_k", "secret_access_key": "aws_s"},
                "earthdata": {"username": "ed_u", "password": "ed_p"},
                "gcs": {"credentials_file": "/tmp/gcs.json"},
            }
        )

        mgr = CredentialManager.from_config(config)

        assert mgr.get_credentials("cdse") == CredentialSpec(key="cdse_u", secret="cdse_p")
        assert mgr.get_credentials("cds") == CredentialSpec(key="cds_k", secret=None)
        assert mgr.get_credentials("aws") == CredentialSpec(key="aws_k", secret="aws_s")
        assert mgr.get_credentials("earthdata") == CredentialSpec(key="ed_u", secret="ed_p")
        assert mgr.get_credentials("gcs") == CredentialSpec(key="/tmp/gcs.json", secret=None)

    def test_reads_cdse_env_vars_via_config_overlay(self, monkeypatch):
        monkeypatch.setenv("SIAC_CDSE_USERNAME", "env_u")
        monkeypatch.setenv("SIAC_CDSE_PASSWORD", "env_p")

        mgr = CredentialManager.from_config(SIACConfig().with_env_overlay())

        cred = mgr.get_credentials("cdse")
        assert cred.key == "env_u"
        assert cred.secret == "env_p"

    def test_reads_native_provider_env_vars_via_config_overlay(self, monkeypatch):
        monkeypatch.setenv("CDSAPI_KEY", "cds-env")
        monkeypatch.setenv("AWS_ACCESS_KEY_ID", "aws-env-k")
        monkeypatch.setenv("AWS_SECRET_ACCESS_KEY", "aws-env-s")
        monkeypatch.setenv("EARTHDATA_USERNAME", "earth-u")
        monkeypatch.setenv("EARTHDATA_PASSWORD", "earth-p")
        monkeypatch.setenv("GOOGLE_APPLICATION_CREDENTIALS", "/tmp/gcs-native.json")

        mgr = CredentialManager.from_config(SIACConfig().with_env_overlay())

        assert mgr.get_credentials("cds") == CredentialSpec(key="cds-env", secret=None)
        assert mgr.get_credentials("aws") == CredentialSpec(key="aws-env-k", secret="aws-env-s")
        assert mgr.get_credentials("earthdata") == CredentialSpec(key="earth-u", secret="earth-p")
        assert mgr.get_credentials("gcs") == CredentialSpec(key="/tmp/gcs-native.json", secret=None)

    def test_config_takes_priority_over_env(self, monkeypatch):
        monkeypatch.setenv("SIAC_CDSE_USERNAME", "env_u")
        monkeypatch.setenv("SIAC_CDSE_PASSWORD", "env_p")

        mgr = CredentialManager.from_config(
            SIACConfig(
                auth={"cdse": {"username": "cfg_u", "password": "cfg_p"}},
            ).with_env_overlay()
        )

        cred = mgr.get_credentials("cdse")
        assert cred.key == "cfg_u"
        assert cred.secret == "cfg_p"

    def test_config_credentials_are_completed_from_env(self, monkeypatch):
        monkeypatch.setenv("EARTHDATA_USERNAME", "env-user")
        monkeypatch.setenv("EARTHDATA_PASSWORD", "env-pass")

        mgr = CredentialManager.from_config(
            SIACConfig(
                auth={"earthdata": {"username": "cfg-user"}},
            ).with_env_overlay()
        )

        assert mgr.get_credentials("earthdata") == CredentialSpec(key="cfg-user", secret="env-pass")


# ── 3. CDSE token caching ───────────────────────────────────────────


class TestCDSEToken:
    def _make_mgr_with_cdse_creds(self):
        mgr = CredentialManager()
        mgr.set_credentials("cdse", key="user", secret="pass")
        return mgr

    @patch("siac.adapters.auth._cdse_token_exchange")
    def test_caches_token(self, mock_exchange):
        """Two calls should only trigger one HTTP exchange."""
        mock_exchange.return_value = ("tok123", 3600)
        mgr = self._make_mgr_with_cdse_creds()

        t1 = mgr.cdse().get_token()
        t2 = mgr.cdse().get_token()

        assert t1 == "tok123"
        assert t2 == "tok123"
        assert mock_exchange.call_count == 1

    @patch("siac.adapters.auth._cdse_token_exchange")
    @patch("siac.adapters.auth.time")
    def test_refreshes_expired_token(self, mock_time, mock_exchange):
        """Expired token triggers a new exchange."""
        mock_time.monotonic = MagicMock(
            side_effect=[
                0.0,  # first get_cdse_token: check cache
                0.0,  # first get_cdse_token: set expires_at = 0 + 300 = 300
                400.0,  # second get_cdse_token: check cache (expired)
                400.0,  # second get_cdse_token: set expires_at = 400 + 300 = 700
            ]
        )
        mock_exchange.side_effect = [("tok1", 300), ("tok2", 300)]

        mgr = self._make_mgr_with_cdse_creds()
        t1 = mgr.cdse().get_token()
        t2 = mgr.cdse().get_token()

        assert t1 == "tok1"
        assert t2 == "tok2"
        assert mock_exchange.call_count == 2

    @patch("siac.adapters.auth._cdse_token_exchange")
    def test_thread_safety(self, mock_exchange):
        """10 concurrent calls should not raise."""
        mock_exchange.return_value = ("tok_thread", 3600)
        mgr = self._make_mgr_with_cdse_creds()
        results: list[str] = []
        errors: list[Exception] = []

        def _get():
            try:
                results.append(mgr.cdse().get_token())
            except Exception as exc:
                errors.append(exc)

        threads = [threading.Thread(target=_get) for _ in range(10)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        assert not errors
        assert all(r == "tok_thread" for r in results)

    def test_missing_cdse_creds_raises(self):
        mgr = CredentialManager()
        with pytest.raises(AuthenticationError, match="cdse"):
            mgr.cdse().get_token()

    def test_incomplete_cdse_creds_raises(self):
        mgr = CredentialManager()
        mgr.set_credentials("cdse", key="user", secret=None)
        with pytest.raises(AuthenticationError, match="username.*password|password.*username"):
            mgr.cdse().get_token()

    @patch("siac.adapters.auth._cdse_token_exchange")
    @patch("siac.adapters.auth.time")
    def test_existing_newer_token_is_kept(self, mock_time, mock_exchange):
        mgr = self._make_mgr_with_cdse_creds()
        # Stale wrt margin (so refresh path is entered), but still newer than fetched token.
        mgr._oauth_tokens["cdse"] = OAuthToken(access_token="existing", expires_at=1000.0)

        mock_time.monotonic = MagicMock(side_effect=[950.0, 950.0])
        mock_exchange.return_value = ("fresh", 1)  # expires_at=951

        token = mgr.cdse().get_token()
        assert token == "existing"


# ── 4. storage adapter helpers ───────────────────────────────────────


class TestStorageOptions:
    def test_aws_storage_options(self):
        mgr = CredentialManager()
        mgr.set_credentials("aws", key="AK", secret="SK")
        assert isinstance(mgr.aws(), AWSAuth)
        assert mgr.aws().storage_options() == {"key": "AK", "secret": "SK"}

    def test_gcs_storage_options(self):
        mgr = CredentialManager()
        mgr.set_credentials("gcs", key="/path/to/creds.json")
        assert isinstance(mgr.gcs(), GCSAuth)
        assert mgr.gcs().storage_options() == {"token": "/path/to/creds.json"}

    def test_unsupported_provider_adapter_requires_credentials(self):
        mgr = CredentialManager()
        with pytest.raises(AuthenticationError, match="No credentials registered"):
            mgr.aws().storage_options()


# ── 5. CredentialSpec immutability ───────────────────────────────────


class TestCredentialSpec:
    def test_frozen(self):
        cs = CredentialSpec(key="k", secret="s")
        with pytest.raises(AttributeError):
            cs.key = "new"  # type: ignore[misc]


class TestProviderAdapters:
    def test_cdse_adapter_roundtrip(self):
        mgr = CredentialManager()
        assert isinstance(mgr.cdse(), CDSEAuth)

    def test_cdse_temporary_s3_credentials_roundtrip(self, monkeypatch):
        """The CDSE auth path goes through a shared retry-enabled requests.Session
        (see ``_get_session`` in ``siac.adapters.auth``), so patching
        ``requests.post``/``requests.delete`` at the module level no longer
        intercepts the call — patch the session instead.
        """
        from unittest.mock import MagicMock

        class _Resp:
            def __init__(self, body):
                self._body = body

            def raise_for_status(self):
                return None

            def json(self):
                return self._body

        session = MagicMock()
        session.post.return_value = _Resp({"access_id": "ak", "secret": "sk"})
        session.delete.return_value = _Resp({})
        monkeypatch.setattr("siac.adapters.auth._get_session", lambda: session)

        mgr = CredentialManager()
        mgr.set_credentials("cdse", key="user", secret="pass")
        mgr._oauth_tokens["cdse"] = OAuthToken(access_token="tok", expires_at=float("inf"))

        creds = mgr.cdse().create_temporary_s3_credentials()
        assert creds.access_key_id == "ak"
        assert creds.secret_access_key == "sk"
        assert creds.bucket == "eodata"
        assert creds.storage_options()["key"] == "ak"

        mgr.cdse().revoke_temporary_s3_credentials("ak")
        assert session.post.call_count == 1
        assert session.delete.call_count == 1

    def test_cdse_temporary_s3_context_verifies_with_retry(self, monkeypatch):
        """Same session-vs-requests patching fix as above (see comment)."""
        from unittest.mock import MagicMock

        class _Resp:
            def __init__(self, body):
                self._body = body

            def raise_for_status(self):
                return None

            def json(self):
                return self._body

        attempts = {"count": 0}

        class _FakeS3FS:
            def __init__(self, **kwargs):
                self.kwargs = kwargs

            def ls(self, path, detail=False):
                attempts["count"] += 1
                if attempts["count"] < 3:
                    raise PermissionError("InvalidAccessKeyId")
                return [f"{path}/ok"]

        session = MagicMock()
        session.post.return_value = _Resp({"access_id": "ak", "secret": "sk"})
        session.delete.return_value = _Resp({})
        monkeypatch.setattr("siac.adapters.auth._get_session", lambda: session)

        mgr = CredentialManager()
        mgr.set_credentials("cdse", key="user", secret="pass")
        mgr._oauth_tokens["cdse"] = OAuthToken(access_token="tok", expires_at=float("inf"))
        monkeypatch.setitem(sys.modules, "s3fs", SimpleNamespace(S3FileSystem=_FakeS3FS))
        monkeypatch.setattr("siac.adapters.auth.time.sleep", lambda _seconds: None)

        with mgr.cdse().temporary_s3_credentials() as creds:
            assert creds.access_key_id == "ak"
            assert creds.secret_access_key == "sk"

        assert attempts["count"] == 3
        assert session.delete.call_count == 1

    def test_cdse_temporary_s3_credentials_missing_fields_raise(self, monkeypatch):
        from unittest.mock import MagicMock

        class _Resp:
            def __init__(self, body):
                self._body = body

            def raise_for_status(self):
                return None

            def json(self):
                return self._body

        session = MagicMock()
        session.post.return_value = _Resp({"access_id": "ak"})
        monkeypatch.setattr("siac.adapters.auth._get_session", lambda: session)

        mgr = CredentialManager()
        mgr.set_credentials("cdse", key="user", secret="pass")
        mgr._oauth_tokens["cdse"] = OAuthToken(access_token="tok", expires_at=float("inf"))

        with pytest.raises(AuthenticationError, match="missing secret"):
            mgr.cdse().create_temporary_s3_credentials()

    def test_cds_adapter_client_kwargs_from_store(self):
        mgr = CredentialManager()
        mgr.set_credentials("cds", key="store-key")
        assert mgr.cds().client_kwargs() == {"key": "store-key"}

    def test_cds_make_client_requires_configured_credentials(self):
        mgr = CredentialManager()
        with pytest.raises(AuthenticationError, match="config.auth.cds.api_key"):
            mgr.cds().make_client()

    def test_earthdata_adapter_builds_earthaccess_source(self):
        mgr = CredentialManager()
        mgr.set_credentials("earthdata", key="user", secret="pass")
        assert isinstance(mgr.earthdata(), EarthdataAuth)

        source = mgr.earthdata().build_earthaccess_source(provider="LPDAAC_ECS")
        assert source.provider == "LPDAAC_ECS"
        assert source.earthdata_username == "user"
        assert source.earthdata_password == "pass"
        assert source.login_strategy == "environment"
        assert source.persist is False

    def test_earthdata_adapter_defaults_to_native_earthaccess_flow(self):
        mgr = CredentialManager()
        source = mgr.earthdata().build_earthaccess_source(provider="LPDAAC_ECS")
        assert source.provider == "LPDAAC_ECS"
        assert source.earthdata_username is None
        assert source.earthdata_password is None
        assert source.login_strategy == "all"
        assert source.persist is False

    def test_earthdata_adapter_incomplete_credentials_raise(self):
        mgr = CredentialManager()
        mgr.set_credentials("earthdata", key="user", secret=None)
        with pytest.raises(AuthenticationError, match="Earthdata credentials require"):
            mgr.earthdata().build_earthaccess_source()


# ── 6. CDSE backend integration (mock) ──────────────────────────────


class TestCDSEBackendWithAuth:
    @patch("siac.adapters.auth._cdse_token_exchange")
    def test_backend_uses_auth_manager_token(self, mock_exchange):
        """CopernicusDataspaceBackend uses auth manager when provided."""
        mock_exchange.return_value = ("mgr_token", 3600)

        from siac.adapters.data.copernicus_dataspace import CopernicusDataspaceBackend

        mgr = CredentialManager()
        mgr.set_credentials("cdse", key="u", secret="p")

        backend = CopernicusDataspaceBackend(auth=mgr)
        header = backend._get_auth_header()
        assert header == {"Authorization": "Bearer mgr_token"}

    def test_backend_explicit_keys_take_priority(self):
        """Explicit access_key/secret_key still take priority when supplied."""
        from siac.adapters.data.copernicus_dataspace import CopernicusDataspaceBackend

        mgr = CredentialManager()  # no CDSE creds
        backend = CopernicusDataspaceBackend(
            access_key="direct_token",
            auth=mgr,
        )
        header = backend._get_auth_header()
        assert header == {"Authorization": "Bearer direct_token"}


class TestCDSETokenExchange:
    def test_exchange_success_and_missing_token(self, monkeypatch):
        from unittest.mock import MagicMock

        class _Resp:
            def __init__(self, body):
                self._body = body

            def raise_for_status(self):
                return None

            def json(self):
                return self._body

        session = MagicMock()
        session.post.return_value = _Resp({"access_token": "abc", "expires_in": 123})
        monkeypatch.setattr("siac.adapters.auth._get_session", lambda: session)

        token, exp = _cdse_token_exchange("u", "p")
        assert token == "abc"
        assert exp == 123

        session.post.return_value = _Resp({"expires_in": 123})
        with pytest.raises(AuthenticationError, match="access_token"):
            _cdse_token_exchange("u", "p")
