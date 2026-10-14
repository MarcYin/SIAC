"""Optional imports distinguish missing packages from broken installations."""

from __future__ import annotations

import builtins
import importlib
import runpy
from datetime import datetime
from typing import TYPE_CHECKING

import pytest

from siac.adapters import auth, earthdata_common
from siac.adapters.atmo.cams import CAMSProvider
from siac.adapters.data.earthaccess_source import EarthAccessSource
from siac.algorithms.cloud.providers.omnicloudmask import OmniCloudMaskProvider
from siac.algorithms.rt.lut.http_zip_store import _build_s3_filesystem
from siac.errors import AuthenticationError

if TYPE_CHECKING:
    from pathlib import Path


def _import_optional(package: str, tmp_path: Path):
    if package == "earthaccess":
        return EarthAccessSource()._get_module()
    if package == "omnicloudmask":
        return OmniCloudMaskProvider._default_predictor(None)
    if package == "s3fs":
        return _build_s3_filesystem({})
    if package == "cdsapi":
        return CAMSProvider(tmp_path)._download_cams_file(datetime(2024, 1, 2))
    return runpy.run_path(earthdata_common.__file__)


def _fail_import(monkeypatch: pytest.MonkeyPatch, package: str, error: ImportError) -> None:
    original_import = builtins.__import__
    original_import_module = importlib.import_module

    def import_package(name, *args, **kwargs):
        if name.split(".")[0] == package:
            raise error
        return original_import(name, *args, **kwargs)

    def import_module(name, *args, **kwargs):
        if name == package:
            raise error
        return original_import_module(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", import_package)
    monkeypatch.setattr(importlib, "import_module", import_module)
    monkeypatch.setattr(auth, "import_module", import_module)
    monkeypatch.setattr("siac.adapters.atmo.cams.import_module", import_module)
    monkeypatch.setattr(
        "siac.algorithms.cloud.providers.omnicloudmask.import_module", import_module
    )


@pytest.mark.parametrize("package", ["earthaccess", "omnicloudmask", "s3fs", "cdsapi", "pyhdf"])
@pytest.mark.parametrize("failure", ["transitive_missing", "incompatible"])
def test_optional_installed_package_import_failure(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, package: str, failure: str
) -> None:
    error = (
        ModuleNotFoundError("No module named dependency", name="dependency")
        if failure == "transitive_missing"
        else ImportError("cannot import name from dependency")
    )
    _fail_import(monkeypatch, package, error)
    with pytest.raises(ImportError, match=f"{package} is installed but failed to import") as raised:
        _import_optional(package, tmp_path)
    assert raised.value.__cause__ is error


@pytest.mark.parametrize("package", ["earthaccess", "omnicloudmask", "s3fs", "cdsapi", "pyhdf"])
def test_optional_missing_package_keeps_install_hint(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, package: str
) -> None:
    error = ModuleNotFoundError(f"No module named {package}", name=package)
    _fail_import(monkeypatch, package, error)
    if package == "cdsapi":
        assert _import_optional(package, tmp_path) is None
    elif package == "pyhdf":
        module = _import_optional(package, tmp_path)
        with pytest.raises(ModuleNotFoundError, match="pyhdf is required"):
            module["SD"]("unused.hdf")
    else:
        with pytest.raises(ImportError, match=f"{package} is required") as raised:
            _import_optional(package, tmp_path)
        assert raised.value.__cause__ is error


def test_auth_optional_import_preserves_missing_dependency(monkeypatch: pytest.MonkeyPatch) -> None:
    error = ModuleNotFoundError("No module named dependency", name="dependency")
    _fail_import(monkeypatch, "cdsapi", error)
    with pytest.raises(ModuleNotFoundError) as raised:
        auth._load_optional_dependency("cdsapi", "cdsapi is not installed.")
    assert raised.value is error


def test_auth_optional_import_keeps_missing_package_hint(monkeypatch: pytest.MonkeyPatch) -> None:
    error = ModuleNotFoundError("No module named cdsapi", name="cdsapi")
    _fail_import(monkeypatch, "cdsapi", error)
    with pytest.raises(AuthenticationError, match="cdsapi is not installed") as raised:
        auth._load_optional_dependency("cdsapi", "cdsapi is not installed.")
    assert raised.value.__cause__ is error
