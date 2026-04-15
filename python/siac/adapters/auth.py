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
from importlib import import_module
from typing import TYPE_CHECKING, Any, cast
from urllib.parse import quote

from siac.errors import AuthenticationError

if TYPE_CHECKING:
    from collections.abc import Iterator

    from siac.adapters.data.earthaccess_source import EarthAccessSource
    from siac.config import SIACConfig

requests = cast("Any", import_module("requests"))

logger = logging.getLogger(__name__)

from siac.adapters._log_filter import SecretRedactionFilter
logger.addFilter(SecretRedactionFilter())

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


def _load_optional_dependency(module_name: str, install_hint: str) -> Any:
    try:
        return import_module(module_name)
    except ModuleNotFoundError as exc:  # pragma: no cover - env-dependent
        raise AuthenticationError(install_hint) from exc


def _response_json_object(response: Any, error_message: str) -> dict[str, Any]:
    body = response.json()
    if not isinstance(body, dict):
        raise AuthenticationError(error_message)
    return cast("dict[str, Any]", body)


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
        try:
            resp.raise_for_status()
        except requests.HTTPError as exc:
            # Wrap to prevent URL tokens leaking in the error message.
            raise AuthenticationError(
                f"CDSE S3 credential request failed (HTTP {resp.status_code}). "
                "Check your CDSE credentials."
            ) from exc
        body = _response_json_object(
            resp,
            "CDSE temporary S3 credential response must be a JSON object.",
        )
        access_id = body.get("access_id")
        secret = body.get("secret")
        if not isinstance(access_id, str) or not access_id:
            raise AuthenticationError("CDSE temporary S3 credential response missing access_id.")
        if not isinstance(secret, str) or not secret:
            raise AuthenticationError("CDSE temporary S3 credential response missing secret.")
        return CDSES3Credentials(access_key_id=access_id, secret_access_key=secret)

    def revoke_temporary_s3_credentials(self, access_key_id: str, timeout: int = 60) -> None:
        resp = requests.delete(
            f"{_CDSE_S3_CREDENTIALS_URL}/access_id/{quote(access_key_id, safe='')}",
            headers={
                **self.authorization_header(),
                "Accept": "application/json",
            },
            timeout=timeout,
        )
        try:
            resp.raise_for_status()
        except requests.HTTPError as exc:
            raise AuthenticationError(
                f"CDSE S3 credential revocation failed (HTTP {resp.status_code})."
            ) from exc

    def wait_for_temporary_s3_credentials(
        self,
        credentials: CDSES3Credentials,
        activation_delays: tuple[int, ...] = (0, 1, 2, 4),
    ) -> None:
        s3fs = _load_optional_dependency(
            "s3fs",
            "s3fs is required to verify CDSE S3 credentials.",
        )

        last_exc: Exception | None = None
        for delay in activation_delays:
            if delay > 0:
                time.sleep(delay)
            try:
                fs = s3fs.S3FileSystem(**credentials.storage_options())
                fs.ls(credentials.bucket, detail=False)[:1]
                return
            except (OSError, PermissionError, ConnectionError) as exc:  # pragma: no cover
                last_exc = exc

        # Avoid leaking credential details in the exception chain.
        raise AuthenticationError(
            "CDSE temporary S3 credentials did not become active "
            f"after {len(activation_delays)} attempts."
        ) from last_exc

    @contextmanager
    def temporary_s3_credentials(
        self,
        *,
        timeout: int = 60,
        activation_delays: tuple[int, ...] = (0, 1, 2, 4),
        verify: bool = True,
    ) -> Iterator[CDSES3Credentials]:
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

    def client_kwargs(self) -> dict[str, str]:
        if not self.has_credentials():
            return {}
        cred = self.get_credentials()
        kwargs: dict[str, str] = {}
        if cred.key:
            kwargs["key"] = cred.key
        return kwargs

    def make_client(self, **kwargs: Any) -> Any:
        if not self.has_credentials():
            raise AuthenticationError(
                "CDS credentials are not configured; populate config.auth.cds.api_key."
            )
        cdsapi = _load_optional_dependency("cdsapi", "cdsapi is not installed.")

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
    def activate_environment(self) -> Iterator[None]:
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

    def build_earthaccess_source(
        self,
        *,
        provider: str | None = None,
        **kwargs: Any,
    ) -> EarthAccessSource:
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
        manager = cls()

        if config is not None:
            auth_config = getattr(config, "auth", None)
            if auth_config is not None:
                _load_from_auth_config(manager, auth_config)

        return manager


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
    try:
        resp.raise_for_status()
    except requests.HTTPError as exc:
        raise AuthenticationError(
            f"CDSE token exchange failed (HTTP {resp.status_code}). "
            "Check your CDSE username and password."
        ) from exc
    body = _response_json_object(resp, "CDSE token response must be a JSON object.")
    token = body.get("access_token")
    if not isinstance(token, str) or not token:
        raise AuthenticationError("CDSE token response does not contain access_token.")
    expires_in_raw = body.get("expires_in", 300)
    try:
        expires_in = int(expires_in_raw)
    except (TypeError, ValueError) as exc:
        raise AuthenticationError("CDSE token response contains an invalid expires_in value.") from exc
    return token, expires_in


def _load_from_auth_config(manager: CredentialManager, auth_config: Any) -> None:
    """Populate *manager* from the nested auth config model."""
    cdse = getattr(auth_config, "cdse", None)
    cds = getattr(auth_config, "cds", None)
    aws = getattr(auth_config, "aws", None)
    earthdata = getattr(auth_config, "earthdata", None)
    gcs = getattr(auth_config, "gcs", None)

    if cdse is not None:
        _maybe_set(manager, "cdse", getattr(cdse, "username", None), getattr(cdse, "password", None))
    if cds is not None:
        _maybe_set(manager, "cds", getattr(cds, "api_key", None), None)
    if aws is not None:
        _maybe_set(
            manager,
            "aws",
            getattr(aws, "access_key_id", None),
            getattr(aws, "secret_access_key", None),
        )
    if earthdata is not None:
        _maybe_set(
            manager,
            "earthdata",
            getattr(earthdata, "username", None),
            getattr(earthdata, "password", None),
        )
    if gcs is not None:
        gcs_file = getattr(gcs, "credentials_file", None)
        if gcs_file is not None:
            _maybe_set(manager, "gcs", str(gcs_file), None)


def _maybe_set(
    manager: CredentialManager,
    provider: str,
    key: str | None,
    secret: str | None,
) -> None:
    if key is not None or secret is not None:
        manager.set_credentials(provider, key=key, secret=secret)


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
