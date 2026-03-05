"""Tests for the centralised authentication manager (core/auth.py)."""

from __future__ import annotations

import threading
from unittest.mock import MagicMock, patch

import pytest

from siac.core.auth import (
    CredentialManager,
    CredentialSpec,
    OAuthToken,
    _cdse_token_exchange,
)
from siac.core.exceptions import AuthenticationError

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


# ── 2. from_config with CredentialConfig ─────────────────────────────

class TestFromConfig:
    def test_reads_credential_config_fields(self):
        """Config fields are loaded into the manager."""
        # Build a minimal mock that looks like SIACConfig.credentials
        cred_cfg = MagicMock()
        cred_cfg.cdse_username = "cdse_u"
        cred_cfg.cdse_password = "cdse_p"
        cred_cfg.cds_api_key = "cds_k"
        cred_cfg.aws_access_key_id = "aws_k"
        cred_cfg.aws_secret_access_key = "aws_s"
        cred_cfg.earthdata_username = "ed_u"
        cred_cfg.earthdata_password = "ed_p"
        cred_cfg.gcs_credentials_file = "/tmp/gcs.json"

        config = MagicMock()
        config.credentials = cred_cfg

        mgr = CredentialManager.from_config(config)

        assert mgr.get_credentials("cdse") == CredentialSpec(key="cdse_u", secret="cdse_p")
        assert mgr.get_credentials("cds") == CredentialSpec(key="cds_k", secret=None)
        assert mgr.get_credentials("aws") == CredentialSpec(key="aws_k", secret="aws_s")
        assert mgr.get_credentials("earthdata") == CredentialSpec(key="ed_u", secret="ed_p")
        assert mgr.get_credentials("gcs") == CredentialSpec(key="/tmp/gcs.json", secret=None)

    def test_reads_env_vars(self, monkeypatch):
        """Env vars are picked up when no config is provided."""
        monkeypatch.setenv("SIAC_CDSE_USERNAME", "env_u")
        monkeypatch.setenv("SIAC_CDSE_PASSWORD", "env_p")

        mgr = CredentialManager.from_config(None)

        cred = mgr.get_credentials("cdse")
        assert cred.key == "env_u"
        assert cred.secret == "env_p"

    def test_config_takes_priority_over_env(self, monkeypatch):
        """Config fields win over env vars."""
        monkeypatch.setenv("SIAC_CDSE_USERNAME", "env_u")
        monkeypatch.setenv("SIAC_CDSE_PASSWORD", "env_p")

        cred_cfg = MagicMock()
        cred_cfg.cdse_username = "cfg_u"
        cred_cfg.cdse_password = "cfg_p"
        cred_cfg.cds_api_key = None
        cred_cfg.aws_access_key_id = None
        cred_cfg.aws_secret_access_key = None
        cred_cfg.earthdata_username = None
        cred_cfg.earthdata_password = None
        cred_cfg.gcs_credentials_file = None

        config = MagicMock()
        config.credentials = cred_cfg

        mgr = CredentialManager.from_config(config)

        cred = mgr.get_credentials("cdse")
        assert cred.key == "cfg_u"
        assert cred.secret == "cfg_p"

    def test_aws_standard_fallback(self, monkeypatch):
        """Standard AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY are used as fallback."""
        monkeypatch.setenv("AWS_ACCESS_KEY_ID", "std_k")
        monkeypatch.setenv("AWS_SECRET_ACCESS_KEY", "std_s")

        mgr = CredentialManager.from_config(None)

        cred = mgr.get_credentials("aws")
        assert cred.key == "std_k"
        assert cred.secret == "std_s"

    def test_cdsapirc_fallback(self, monkeypatch, tmp_path):
        """~/.cdsapirc is read when no CDS credentials are set."""
        rc_file = tmp_path / ".cdsapirc"
        rc_file.write_text("url: https://cds.climate.copernicus.eu/api\nkey: my-cds-key\n")

        monkeypatch.setattr("siac.core.auth.Path.home", staticmethod(lambda: tmp_path))

        mgr = CredentialManager.from_config(None)

        cred = mgr.get_credentials("cds")
        assert cred.key == "my-cds-key"

    def test_cdsapirc_read_error_is_ignored(self, monkeypatch, tmp_path):
        rc_file = tmp_path / ".cdsapirc"
        rc_file.write_text("key: ignored")
        monkeypatch.setattr("siac.core.auth.Path.home", staticmethod(lambda: tmp_path))
        monkeypatch.setattr("siac.core.auth.Path.read_text", lambda _self: (_ for _ in ()).throw(OSError("nope")))

        mgr = CredentialManager.from_config(None)
        assert mgr.has_credentials("cds") is False


# ── 3. CDSE token caching ───────────────────────────────────────────

class TestCDSEToken:
    def _make_mgr_with_cdse_creds(self):
        mgr = CredentialManager()
        mgr.set_credentials("cdse", key="user", secret="pass")
        return mgr

    @patch("siac.core.auth._cdse_token_exchange")
    def test_caches_token(self, mock_exchange):
        """Two calls should only trigger one HTTP exchange."""
        mock_exchange.return_value = ("tok123", 3600)
        mgr = self._make_mgr_with_cdse_creds()

        t1 = mgr.get_cdse_token()
        t2 = mgr.get_cdse_token()

        assert t1 == "tok123"
        assert t2 == "tok123"
        assert mock_exchange.call_count == 1

    @patch("siac.core.auth._cdse_token_exchange")
    @patch("siac.core.auth.time")
    def test_refreshes_expired_token(self, mock_time, mock_exchange):
        """Expired token triggers a new exchange."""
        mock_time.monotonic = MagicMock(side_effect=[
            0.0,    # first get_cdse_token: check cache
            0.0,    # first get_cdse_token: set expires_at = 0 + 300 = 300
            400.0,  # second get_cdse_token: check cache (expired)
            400.0,  # second get_cdse_token: set expires_at = 400 + 300 = 700
        ])
        mock_exchange.side_effect = [("tok1", 300), ("tok2", 300)]

        mgr = self._make_mgr_with_cdse_creds()
        t1 = mgr.get_cdse_token()
        t2 = mgr.get_cdse_token()

        assert t1 == "tok1"
        assert t2 == "tok2"
        assert mock_exchange.call_count == 2

    @patch("siac.core.auth._cdse_token_exchange")
    def test_thread_safety(self, mock_exchange):
        """10 concurrent calls should not raise."""
        mock_exchange.return_value = ("tok_thread", 3600)
        mgr = self._make_mgr_with_cdse_creds()
        results: list[str] = []
        errors: list[Exception] = []

        def _get():
            try:
                results.append(mgr.get_cdse_token())
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
            mgr.get_cdse_token()

    def test_incomplete_cdse_creds_raises(self):
        mgr = CredentialManager()
        mgr.set_credentials("cdse", key="user", secret=None)
        with pytest.raises(AuthenticationError, match="username.*password|password.*username"):
            mgr.get_cdse_token()

    @patch("siac.core.auth._cdse_token_exchange")
    @patch("siac.core.auth.time")
    def test_existing_newer_token_is_kept(self, mock_time, mock_exchange):
        mgr = self._make_mgr_with_cdse_creds()
        # Stale wrt margin (so refresh path is entered), but still newer than fetched token.
        mgr._oauth_tokens["cdse"] = OAuthToken(access_token="existing", expires_at=1000.0)

        mock_time.monotonic = MagicMock(side_effect=[950.0, 950.0])
        mock_exchange.return_value = ("fresh", 1)  # expires_at=951

        token = mgr.get_cdse_token()
        assert token == "existing"


# ── 4. get_storage_options ───────────────────────────────────────────

class TestStorageOptions:
    def test_aws_storage_options(self):
        mgr = CredentialManager()
        mgr.set_credentials("aws", key="AK", secret="SK")
        opts = mgr.get_storage_options("aws")
        assert opts == {"key": "AK", "secret": "SK"}

    def test_gcs_storage_options(self):
        mgr = CredentialManager()
        mgr.set_credentials("gcs", key="/path/to/creds.json")
        opts = mgr.get_storage_options("gcs")
        assert opts == {"token": "/path/to/creds.json"}

    def test_unsupported_provider_raises(self):
        mgr = CredentialManager()
        mgr.set_credentials("cdse", key="u", secret="p")
        with pytest.raises(AuthenticationError, match="not supported"):
            mgr.get_storage_options("cdse")


# ── 5. CredentialSpec immutability ───────────────────────────────────

class TestCredentialSpec:
    def test_frozen(self):
        cs = CredentialSpec(key="k", secret="s")
        with pytest.raises(AttributeError):
            cs.key = "new"  # type: ignore[misc]


# ── 6. CDSE backend integration (mock) ──────────────────────────────

class TestCDSEBackendWithAuth:
    @patch("siac.core.auth._cdse_token_exchange")
    def test_backend_uses_auth_manager_token(self, mock_exchange):
        """CopernicusDataspaceBackend uses auth manager when provided."""
        mock_exchange.return_value = ("mgr_token", 3600)

        from siac.io.copernicus_dataspace import CopernicusDataspaceBackend

        mgr = CredentialManager()
        mgr.set_credentials("cdse", key="u", secret="p")

        backend = CopernicusDataspaceBackend(auth=mgr)
        header = backend._get_auth_header()
        assert header == {"Authorization": "Bearer mgr_token"}

    def test_backend_explicit_keys_take_priority(self):
        """Explicit access_key/secret_key still work (backward compat)."""
        from siac.io.copernicus_dataspace import CopernicusDataspaceBackend

        mgr = CredentialManager()  # no CDSE creds
        backend = CopernicusDataspaceBackend(
            access_key="direct_token", auth=mgr,
        )
        header = backend._get_auth_header()
        assert header == {"Authorization": "Bearer direct_token"}


class TestCDSETokenExchange:
    def test_exchange_success_and_missing_token(self, monkeypatch):
        class _Resp:
            def __init__(self, body):
                self._body = body

            def raise_for_status(self):
                return None

            def json(self):
                return self._body

        monkeypatch.setattr(
            "siac.core.auth.requests.post",
            lambda *_args, **_kwargs: _Resp({"access_token": "abc", "expires_in": 123}),
        )
        token, exp = _cdse_token_exchange("u", "p")
        assert token == "abc"
        assert exp == 123

        monkeypatch.setattr(
            "siac.core.auth.requests.post",
            lambda *_args, **_kwargs: _Resp({"expires_in": 123}),
        )
        with pytest.raises(AuthenticationError, match="access_token"):
            _cdse_token_exchange("u", "p")
