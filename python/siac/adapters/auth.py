"""
Centralized authentication manager for SIAC.
"""

from __future__ import annotations

import logging
import os
import threading
import time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import requests

from siac.errors import AuthenticationError

if TYPE_CHECKING:
    from siac.config import SIACConfig

logger = logging.getLogger(__name__)

_CDSE_TOKEN_URL = (
    "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/"
    "protocol/openid-connect/token"
)
_CDSE_S3_CREDENTIALS_URL = "https://s3-keys-manager.cloudferro.com/api/user/credentials"
_CDSE_S3_ENDPOINT_URL = "https://eodata.dataspace.copernicus.eu"
_CDSE_S3_BUCKET = "eodata"


@dataclass(frozen=True)
class CredentialSpec:
    """A provider-specific credential pair."""

    key: str | None = None
    secret: str | None = None


@dataclass
class OAuthToken:
    """Cached OAuth2 bearer token with expiry."""

    access_token: str = ""
    expires_at: float = 0.0


@dataclass(frozen=True)
class CDSES3Credentials:
    """Temporary S3 credentials minted from a CDSE bearer token."""

    access_key_id: str
    secret_access_key: str
    endpoint_url: str = _CDSE_S3_ENDPOINT_URL
    bucket: str = _CDSE_S3_BUCKET

    def storage_options(self) -> dict[str, Any]:
        return {
            "key": self.access_key_id,
            "secret": self.secret_access_key,
            "client_kwargs": {"endpoint_url": self.endpoint_url},
        }


class _ProviderAuthBase:
    """Small adapter base for provider-specific auth helpers."""

    provider: str = ""

    def __init__(self, manager: CredentialManager) -> None:
        self._manager = manager

    def has_credentials(self) -> bool:
        return self._manager.has_credentials(self.provider)

    def get_credentials(self) -> CredentialSpec:
        return self._manager.get_credentials(self.provider)


class CDSEAuth(_ProviderAuthBase):
    """CDSE OAuth2 token cache and bearer-header helper."""

    provider = "cdse"

    def get_token(self, margin_seconds: float = 60) -> str:
        now = time.monotonic()

        with self._manager._lock:
            cached = self._manager._oauth_tokens.get("cdse")
            if cached is not None and cached.expires_at - now > margin_seconds:
                return cached.access_token

        cred = self.get_credentials()
        if not cred.key or not cred.secret:
            raise AuthenticationError(
                "CDSE credentials require both username (key) and password (secret)."
            )

        access_token, expires_in = _cdse_token_exchange(cred.key, cred.secret)
        new_token = OAuthToken(
            access_token=access_token,
            expires_at=time.monotonic() + expires_in,
        )

        with self._manager._lock:
            existing = self._manager._oauth_tokens.get("cdse")
            if existing is not None and existing.expires_at > new_token.expires_at:
                return existing.access_token
            self._manager._oauth_tokens["cdse"] = new_token

        return new_token.access_token

    def authorization_header(self, margin_seconds: float = 60) -> dict[str, str]:
        return {"Authorization": f"Bearer {self.get_token(margin_seconds=margin_seconds)}"}

    def create_temporary_s3_credentials(self, timeout: int = 60) -> CDSES3Credentials:
        resp = requests.post(
            _CDSE_S3_CREDENTIALS_URL,
            headers={
                **self.authorization_header(),
                "Accept": "application/json",
            },
            timeout=timeout,
        )
        resp.raise_for_status()
        body = resp.json()
        access_id = body.get("access_id")
        secret = body.get("secret")
        if not isinstance(access_id, str) or not access_id:
            raise AuthenticationError("CDSE temporary S3 credential response missing access_id.")
        if not isinstance(secret, str) or not secret:
            raise AuthenticationError("CDSE temporary S3 credential response missing secret.")
        return CDSES3Credentials(access_key_id=access_id, secret_access_key=secret)

    def revoke_temporary_s3_credentials(self, access_key_id: str, timeout: int = 60) -> None:
        resp = requests.delete(
            f"{_CDSE_S3_CREDENTIALS_URL}/access_id/{access_key_id}",
            headers={
                **self.authorization_header(),
                "Accept": "application/json",
            },
            timeout=timeout,
        )
        resp.raise_for_status()

    def wait_for_temporary_s3_credentials(
        self,
        credentials: CDSES3Credentials,
        activation_delays: tuple[int, ...] = (0, 1, 2, 4),
    ) -> None:
        try:
            import s3fs  # type: ignore[import-not-found]
        except Exception as exc:  # pragma: no cover
            raise AuthenticationError("s3fs is required to verify CDSE S3 credentials.") from exc

        last_exc: Exception | None = None
        for delay in activation_delays:
            if delay > 0:
                time.sleep(delay)
            try:
                fs = s3fs.S3FileSystem(**credentials.storage_options())
                fs.ls(credentials.bucket, detail=False)[:1]
                return
            except Exception as exc:  # pragma: no cover
                last_exc = exc

        raise AuthenticationError("CDSE temporary S3 credentials did not become active.") from last_exc

    @contextmanager
    def temporary_s3_credentials(
        self,
        *,
        timeout: int = 60,
        activation_delays: tuple[int, ...] = (0, 1, 2, 4),
        verify: bool = True,
    ) -> Any:
        credentials = self.create_temporary_s3_credentials(timeout=timeout)
        try:
            if verify:
                self.wait_for_temporary_s3_credentials(
                    credentials,
                    activation_delays=activation_delays,
                )
            yield credentials
        finally:
            self.revoke_temporary_s3_credentials(credentials.access_key_id, timeout=timeout)


class CDSAuth(_ProviderAuthBase):
    """CDS API auth helper for cdsapi clients."""

    provider = "cds"

    def has_external_credentials(self) -> bool:
        return bool(os.getenv("CDSAPI_KEY")) or Path.home().joinpath(".cdsapirc").exists()

    def has_any_credentials(self) -> bool:
        return self.has_credentials() or self.has_external_credentials()

    def client_kwargs(self) -> dict[str, str]:
        if not self.has_credentials():
            return {}
        cred = self.get_credentials()
        kwargs: dict[str, str] = {}
        if cred.key:
            kwargs["key"] = cred.key
        return kwargs

    def make_client(self, **kwargs: Any) -> Any:
        if not self.has_any_credentials():
            raise AuthenticationError(
                "CDS credentials are not configured; set CDSAPI_KEY or ~/.cdsapirc."
            )
        try:
            import cdsapi  # type: ignore[import-not-found]
        except Exception as exc:  # pragma: no cover
            raise AuthenticationError("cdsapi is not installed.") from exc

        client_kwargs = self.client_kwargs()
        client_kwargs.update(kwargs)
        return cdsapi.Client(**client_kwargs)


class AWSAuth(_ProviderAuthBase):
    """AWS/S3 storage-option helper."""

    provider = "aws"

    def storage_options(self) -> dict[str, Any]:
        cred = self.get_credentials()
        opts: dict[str, Any] = {}
        if cred.key:
            opts["key"] = cred.key
        if cred.secret:
            opts["secret"] = cred.secret
        return opts


class GCSAuth(_ProviderAuthBase):
    """GCS storage-option helper."""

    provider = "gcs"

    def storage_options(self) -> dict[str, Any]:
        cred = self.get_credentials()
        opts: dict[str, Any] = {}
        if cred.key:
            opts["token"] = cred.key
        return opts


class EarthdataAuth(_ProviderAuthBase):
    """Earthdata helper for environment activation and source construction."""

    provider = "earthdata"

    def _complete_credentials(self) -> CredentialSpec | None:
        if not self.has_credentials():
            return None
        cred = self.get_credentials()
        if cred.key and cred.secret:
            return cred
        raise AuthenticationError(
            "Earthdata credentials require both username (key) and password (secret)."
        )

    @contextmanager
    def activate_environment(self) -> Any:
        cred = self._complete_credentials()
        if cred is None:
            yield
            return

        previous_username = os.environ.get("EARTHDATA_USERNAME")
        previous_password = os.environ.get("EARTHDATA_PASSWORD")
        os.environ["EARTHDATA_USERNAME"] = cred.key or ""
        os.environ["EARTHDATA_PASSWORD"] = cred.secret or ""
        try:
            yield
        finally:
            if previous_username is None:
                os.environ.pop("EARTHDATA_USERNAME", None)
            else:
                os.environ["EARTHDATA_USERNAME"] = previous_username

            if previous_password is None:
                os.environ.pop("EARTHDATA_PASSWORD", None)
            else:
                os.environ["EARTHDATA_PASSWORD"] = previous_password

    def source_kwargs(self, *, provider: str | None = None) -> dict[str, Any]:
        kwargs: dict[str, Any] = {
            "provider": provider,
            "login_strategy": "all",
            "persist": False,
        }
        cred = self._complete_credentials()
        if cred is not None:
            kwargs.update(
                earthdata_username=cred.key,
                earthdata_password=cred.secret,
                login_strategy="environment",
                persist=False,
            )
        return kwargs

    def build_earthaccess_source(self, *, provider: str | None = None, **kwargs: Any) -> Any:
        from siac.adapters.data.earthaccess_source import EarthAccessSource

        source_kwargs = self.source_kwargs(provider=provider)
        source_kwargs.update(kwargs)
        return EarthAccessSource(**source_kwargs)


class CredentialManager:
    """Thread-safe secret registry and auth-adapter factory for SIAC providers."""

    def __init__(self) -> None:
        self._credentials: dict[str, CredentialSpec] = {}
        self._oauth_tokens: dict[str, OAuthToken] = {}
        self._lock = threading.Lock()
        self._cdse_auth = CDSEAuth(self)
        self._cds_auth = CDSAuth(self)
        self._aws_auth = AWSAuth(self)
        self._gcs_auth = GCSAuth(self)
        self._earthdata_auth = EarthdataAuth(self)

    def set_credentials(
        self,
        provider: str,
        key: str | None = None,
        secret: str | None = None,
    ) -> None:
        with self._lock:
            self._credentials[provider] = CredentialSpec(key=key, secret=secret)

    def get_credentials(self, provider: str) -> CredentialSpec:
        with self._lock:
            cred = self._credentials.get(provider)
        if cred is None:
            raise AuthenticationError(
                f"No credentials registered for provider {provider!r}. "
                "Call set_credentials() or populate config.auth."
            )
        return cred

    def has_credentials(self, provider: str) -> bool:
        with self._lock:
            return provider in self._credentials

    def cdse(self) -> CDSEAuth:
        return self._cdse_auth

    def cds(self) -> CDSAuth:
        return self._cds_auth

    def aws(self) -> AWSAuth:
        return self._aws_auth

    def gcs(self) -> GCSAuth:
        return self._gcs_auth

    def earthdata(self) -> EarthdataAuth:
        return self._earthdata_auth

    @classmethod
    def from_config(cls, config: SIACConfig | None = None) -> CredentialManager:
        mgr = cls()

        if config is not None:
            auth_cfg = getattr(config, "auth", None)
            if auth_cfg is not None:
                _load_from_auth_config(mgr, auth_cfg)
            else:
                cred_cfg = getattr(config, "credentials", None)
                if cred_cfg is not None:
                    _load_from_legacy_credential_config(mgr, cred_cfg)

        return mgr


def _cdse_token_exchange(username: str, password: str, timeout: int = 60) -> tuple[str, int]:
    """Exchange CDSE credentials for an access token."""
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


def _load_from_auth_config(mgr: CredentialManager, auth_cfg: Any) -> None:
    """Populate *mgr* from the nested auth config model."""
    cdse = getattr(auth_cfg, "cdse", None)
    cds = getattr(auth_cfg, "cds", None)
    aws = getattr(auth_cfg, "aws", None)
    earthdata = getattr(auth_cfg, "earthdata", None)
    gcs = getattr(auth_cfg, "gcs", None)

    if cdse is not None:
        _maybe_set(mgr, "cdse", getattr(cdse, "username", None), getattr(cdse, "password", None))
    if cds is not None:
        _maybe_set(mgr, "cds", getattr(cds, "api_key", None), None)
    if aws is not None:
        _maybe_set(
            mgr,
            "aws",
            getattr(aws, "access_key_id", None),
            getattr(aws, "secret_access_key", None),
        )
    if earthdata is not None:
        _maybe_set(
            mgr,
            "earthdata",
            getattr(earthdata, "username", None),
            getattr(earthdata, "password", None),
        )
    if gcs is not None:
        gcs_file = getattr(gcs, "credentials_file", None)
        if gcs_file is not None:
            _maybe_set(mgr, "gcs", str(gcs_file), None)


def _load_from_legacy_credential_config(mgr: CredentialManager, cred_cfg: Any) -> None:
    """Populate *mgr* from the legacy flat credential config shape."""
    _maybe_set(mgr, "cdse", getattr(cred_cfg, "cdse_username", None), getattr(cred_cfg, "cdse_password", None))
    _maybe_set(mgr, "cds", getattr(cred_cfg, "cds_api_key", None), None)
    _maybe_set(
        mgr,
        "aws",
        getattr(cred_cfg, "aws_access_key_id", None),
        getattr(cred_cfg, "aws_secret_access_key", None),
    )
    _maybe_set(
        mgr,
        "earthdata",
        getattr(cred_cfg, "earthdata_username", None),
        getattr(cred_cfg, "earthdata_password", None),
    )

    gcs_file = getattr(cred_cfg, "gcs_credentials_file", None)
    if gcs_file is not None:
        _maybe_set(mgr, "gcs", str(gcs_file), None)


def _maybe_set(
    mgr: CredentialManager,
    provider: str,
    key: str | None,
    secret: str | None,
) -> None:
    if key is not None or secret is not None:
        mgr.set_credentials(provider, key=key, secret=secret)


__all__ = [
    "AWSAuth",
    "CDSAuth",
    "CDSEAuth",
    "CDSES3Credentials",
    "CredentialManager",
    "CredentialSpec",
    "EarthdataAuth",
    "GCSAuth",
    "OAuthToken",
    "_cdse_token_exchange",
]
