"""
Centralised authentication manager for SIAC v2.

Provides a single registry for credentials used by all remote data
providers (CDSE, CDS API, AWS/S3, NASA Earthdata, GCS).  Thread-safe
token caching is included for OAuth2-based providers (CDSE).

Usage::

    auth = CredentialManager.from_config(config)
    token = auth.get_cdse_token()
    opts  = auth.get_storage_options("aws")
"""

from __future__ import annotations

import logging
import os
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import requests

from siac.core.exceptions import AuthenticationError

if TYPE_CHECKING:
    from siac.core.config import SIACConfig

logger = logging.getLogger(__name__)

# CDSE OAuth2 endpoint
_CDSE_TOKEN_URL = (
    "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/"
    "protocol/openid-connect/token"
)

# ── Env-var mapping ──────────────────────────────────────────────────

_ENV_MAP: dict[str, tuple[str, str | None]] = {
    "cdse": ("SIAC_CDSE_USERNAME", "SIAC_CDSE_PASSWORD"),
    "cds": ("SIAC_CDS_API_KEY", None),
    "aws": ("SIAC_AWS_ACCESS_KEY_ID", "SIAC_AWS_SECRET_ACCESS_KEY"),
    "earthdata": ("SIAC_EARTHDATA_USERNAME", "SIAC_EARTHDATA_PASSWORD"),
    "gcs": ("SIAC_GCS_CREDENTIALS_FILE", None),
}

# AWS also falls back to standard boto3 env vars
_AWS_FALLBACK_KEY = "AWS_ACCESS_KEY_ID"
_AWS_FALLBACK_SECRET = "AWS_SECRET_ACCESS_KEY"


# ── Data classes ─────────────────────────────────────────────────────

@dataclass(frozen=True)
class CredentialSpec:
    """A (key, secret) credential pair.

    Interpretation is provider-dependent:
    - CDSE: username / password
    - CDS:  API key / None
    - AWS:  access key id / secret access key
    - Earthdata: username / password
    - GCS:  path to credentials JSON / None
    """

    key: str | None = None
    secret: str | None = None


@dataclass
class OAuthToken:
    """Cached OAuth2 bearer token with expiry."""

    access_token: str = ""
    expires_at: float = 0.0  # time.monotonic() deadline


# ── CredentialManager ────────────────────────────────────────────────

class CredentialManager:
    """Thread-safe credential registry for all SIAC data providers."""

    def __init__(self) -> None:
        self._credentials: dict[str, CredentialSpec] = {}
        self._oauth_tokens: dict[str, OAuthToken] = {}
        self._lock = threading.Lock()

    # ── public API ───────────────────────────────────────────────────

    def set_credentials(
        self,
        provider: str,
        key: str | None = None,
        secret: str | None = None,
    ) -> None:
        """Register credentials for *provider*."""
        with self._lock:
            self._credentials[provider] = CredentialSpec(key=key, secret=secret)

    def get_credentials(self, provider: str) -> CredentialSpec:
        """Retrieve credentials.  Raises ``AuthenticationError`` if missing."""
        with self._lock:
            cred = self._credentials.get(provider)
        if cred is None:
            raise AuthenticationError(
                f"No credentials registered for provider {provider!r}.  "
                f"Call set_credentials() or set the appropriate SIAC_* env vars."
            )
        return cred

    def has_credentials(self, provider: str) -> bool:
        with self._lock:
            return provider in self._credentials

    # ── CDSE OAuth2 token (cached, auto-refresh) ────────────────────

    def get_cdse_token(self, margin_seconds: float = 60) -> str:
        """Return a valid CDSE bearer token, refreshing if necessary.

        Uses double-check locking: the lock protects dict reads/writes
        only — the HTTP call is made outside the lock so it cannot block
        other providers.
        """
        now = time.monotonic()

        # Fast path (no lock on the HTTP call)
        with self._lock:
            cached = self._oauth_tokens.get("cdse")
            if cached is not None and cached.expires_at - now > margin_seconds:
                return cached.access_token

        # Slow path — fetch a new token (outside lock)
        cred = self.get_credentials("cdse")
        if not cred.key or not cred.secret:
            raise AuthenticationError(
                "CDSE credentials require both username (key) and password (secret)."
            )

        access_token, expires_in = _cdse_token_exchange(cred.key, cred.secret)
        new_token = OAuthToken(
            access_token=access_token,
            expires_at=time.monotonic() + expires_in,
        )

        with self._lock:
            # Double-check: another thread may have refreshed in the meantime
            existing = self._oauth_tokens.get("cdse")
            if existing is not None and existing.expires_at > new_token.expires_at:
                return existing.access_token
            self._oauth_tokens["cdse"] = new_token

        return new_token.access_token

    # ── fsspec storage_options helpers ────────────────────────────────

    def get_storage_options(self, provider: str) -> dict[str, Any]:
        """Produce an fsspec-compatible ``storage_options`` dict."""
        cred = self.get_credentials(provider)

        if provider == "aws":
            opts: dict[str, Any] = {}
            if cred.key:
                opts["key"] = cred.key
            if cred.secret:
                opts["secret"] = cred.secret
            return opts

        if provider == "gcs":
            opts = {}
            if cred.key:
                opts["token"] = cred.key
            return opts

        raise AuthenticationError(
            f"get_storage_options() is not supported for provider {provider!r}."
        )

    # ── Factory ──────────────────────────────────────────────────────

    @classmethod
    def from_config(cls, config: SIACConfig | None = None) -> CredentialManager:
        """Build a ``CredentialManager`` by resolving credentials in priority order:

        1. ``config.credentials`` fields (if *config* is provided)
        2. ``SIAC_*`` environment variables
        3. External config files (``~/.cdsapirc`` for CDS)
        """
        mgr = cls()

        # 1. Config fields (highest priority)
        if config is not None:
            cred_cfg = getattr(config, "credentials", None)
            if cred_cfg is not None:
                _load_from_credential_config(mgr, cred_cfg)

        # 2. Env vars (fill in any gaps)
        _load_from_env(mgr)

        # 3. External files (lowest priority)
        _load_from_cdsapirc(mgr)

        return mgr


# ── Private helpers ──────────────────────────────────────────────────

def _cdse_token_exchange(username: str, password: str, timeout: int = 60) -> tuple[str, int]:
    """Exchange CDSE credentials for an access token.

    Returns ``(access_token, expires_in_seconds)``.
    """
    resp = requests.post(
        _CDSE_TOKEN_URL,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        data={
            "client_id": "cdse-public",
            "username": username,
            "password": password,
            "grant_type": "password",
        },
        timeout=timeout,
    )
    resp.raise_for_status()
    body = resp.json()
    token = body.get("access_token")
    if not isinstance(token, str) or not token:
        raise AuthenticationError("CDSE token response does not contain access_token.")
    expires_in = int(body.get("expires_in", 300))
    return token, expires_in


def _load_from_credential_config(mgr: CredentialManager, cred_cfg: Any) -> None:
    """Populate *mgr* from a ``CredentialConfig`` pydantic model."""
    _maybe_set(mgr, "cdse", getattr(cred_cfg, "cdse_username", None), getattr(cred_cfg, "cdse_password", None))
    _maybe_set(mgr, "cds", getattr(cred_cfg, "cds_api_key", None), None)
    _maybe_set(mgr, "aws", getattr(cred_cfg, "aws_access_key_id", None), getattr(cred_cfg, "aws_secret_access_key", None))
    _maybe_set(mgr, "earthdata", getattr(cred_cfg, "earthdata_username", None), getattr(cred_cfg, "earthdata_password", None))

    gcs_file = getattr(cred_cfg, "gcs_credentials_file", None)
    if gcs_file is not None:
        _maybe_set(mgr, "gcs", str(gcs_file), None)


def _load_from_env(mgr: CredentialManager) -> None:
    """Fill gaps from environment variables."""
    for provider, (key_var, secret_var) in _ENV_MAP.items():
        if mgr.has_credentials(provider):
            continue
        key_val = os.environ.get(key_var)
        secret_val = os.environ.get(secret_var) if secret_var else None
        if key_val:
            mgr.set_credentials(provider, key=key_val, secret=secret_val)

    # AWS standard fallback
    if not mgr.has_credentials("aws"):
        aws_key = os.environ.get(_AWS_FALLBACK_KEY)
        aws_secret = os.environ.get(_AWS_FALLBACK_SECRET)
        if aws_key:
            mgr.set_credentials("aws", key=aws_key, secret=aws_secret)


def _load_from_cdsapirc(mgr: CredentialManager) -> None:
    """Read ``~/.cdsapirc`` for CDS API key (lowest priority)."""
    if mgr.has_credentials("cds"):
        return
    rc_path = Path.home() / ".cdsapirc"
    if not rc_path.exists():
        return
    try:
        text = rc_path.read_text()
        for line in text.splitlines():
            if line.strip().startswith("key"):
                parts = line.split(":", 1)
                if len(parts) == 2:
                    api_key = parts[1].strip()
                    if api_key:
                        mgr.set_credentials("cds", key=api_key)
                        return
    except OSError:
        pass


def _maybe_set(
    mgr: CredentialManager,
    provider: str,
    key: str | None,
    secret: str | None,
) -> None:
    """Set credentials only if at least one value is non-None."""
    if key is not None or secret is not None:
        mgr.set_credentials(provider, key=key, secret=secret)
