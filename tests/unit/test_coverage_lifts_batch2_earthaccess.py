"""Coverage lifts for Earthaccess and atmospheric prior wrappers."""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import pytest

from siac.core.auth import CredentialManager
from siac.io.earthaccess_source import EarthAccessSource
from siac.priors.atmospheric.cams import CAMSProvider
from siac.priors.atmospheric.mcd19_earthaccess import MCD19AODProvider, VNP19AODProvider
from siac.priors.atmospheric.merra2 import MERRA2Provider
from siac.priors.brdf.mcd43_earthaccess import (
    MCD19EarthAccessProvider,
    MCD43EarthAccessProvider,
    VNP43EarthAccessProvider,
)


class _FakeEA:
    def __init__(self):
        self.login_calls = []
        self.search_dataset_calls = []
        self.search_data_calls = []

    def login(self, **kwargs):
        self.login_calls.append(kwargs)

    def search_datasets(self, **kwargs):
        self.search_dataset_calls.append(kwargs)
        out = [{"id": 1}, {"id": 2}]
        count = kwargs.get("count")
        if count == 0:
            return []
        if isinstance(count, int) and count > 0:
            return out[:count]
        return out

    def search_data(self, **kwargs):
        self.search_data_calls.append(kwargs)
        out = [{"id": "a"}, {"id": "b"}]
        count = kwargs.get("count")
        if count == 0:
            return []
        if isinstance(count, int) and count > 0:
            return out[:count]
        return out

    def open(self, granules):
        return ["opened", len(granules)]

    def download(self, granules, dest):
        if granules == ["none"]:
            return None
        if granules == ["one"]:
            return str(Path(dest) / "one.h5")
        if granules == ["many"]:
            return [str(Path(dest) / "a.h5"), str(Path(dest) / "b.h5")]
        return 123


class _FakeEATypedLogin(_FakeEA):
    def login(self, **kwargs):
        # Force the TypeError fallback branch in _ensure_auth.
        if kwargs:
            raise TypeError("unexpected kwargs")
        self.login_calls.append(kwargs)


class _StubSource:
    def __init__(self, granules):
        self.granules = granules
        self.calls = []

    def search_granules(self, **kwargs):
        self.calls.append(kwargs)
        return list(self.granules)


class _StubCatalog:
    def __init__(self, name):
        self.name = name

    def resolve_short_name(self, key):
        return f"{self.name}:{key}"


def test_earthaccess_source_auth_and_download_variants(monkeypatch, tmp_path: Path):
    fake = _FakeEA()
    monkeypatch.setitem(sys.modules, "earthaccess", fake)

    src = EarthAccessSource(provider="LP", login_strategy="interactive", persist=True)

    # Search operations should not force authentication.
    ds = src.search_datasets(keyword="k", count=None)
    assert len(ds) == 2
    assert fake.login_calls == []

    # search_datasets count truncation and provider override branch
    out = src.search_datasets(short_name="MCD43A1", provider="PO", count=1)
    assert len(out) == 1
    assert fake.search_dataset_calls[-1]["provider"] == "PO"
    assert fake.search_dataset_calls[-1]["count"] == 1

    # search_granules temporal string branch + count truncation
    gran = src.search_granules(short_name="MCD43A1", temporal="2024-01-01/2024-01-02", count=1)
    assert len(gran) == 1
    assert fake.search_data_calls[-1]["temporal"] == ("2024-01-01", "2024-01-02")
    assert fake.search_data_calls[-1]["count"] == 1

    # download output variants exercise authentication path
    assert src.download_granules(["none"], tmp_path / "d0") == []
    assert fake.login_calls[-1] == {"strategy": "interactive", "persist": True}
    out1 = src.download_granules(["one"], tmp_path / "d1")
    assert len(out1) == 1 and out1[0].name == "one.h5"
    out_many = src.download_granules(["many"], tmp_path / "d2")
    assert len(out_many) == 2
    assert src.download_granules(["weird"], tmp_path / "d3") == []


def test_earthaccess_source_typeerror_login_and_bounds_transform(monkeypatch):
    fake = _FakeEATypedLogin()
    monkeypatch.setitem(sys.modules, "earthaccess", fake)

    src = EarthAccessSource(login_kwargs={"strategy": "not-supported"})
    _ = src.open_granules([{"id": "x"}])
    assert src.is_authenticated is True
    assert fake.login_calls == [{}]

    class _FakeTransformer:
        def transform(self, x, y):
            return x + 1.0, y + 2.0

    monkeypatch.setattr(
        "siac.io.earthaccess_source.Transformer.from_crs",
        lambda *_args, **_kwargs: _FakeTransformer(),
    )
    b = EarthAccessSource.normalize_bounds_to_wgs84((0.0, 0.0, 1.0, 1.0), "EPSG:3857")
    assert b == (1.0, 2.0, 2.0, 3.0)


@pytest.mark.parametrize(
    "provider_cls,key",
    [
        (MCD19AODProvider, "mcd19_aod"),
        (MERRA2Provider, "merra2_atmo"),
    ],
)
def test_atmo_earthaccess_provider_branches(provider_cls, key, monkeypatch):
    source = _StubSource(granules=[])
    catalog = _StubCatalog("SN")
    p = provider_cls(source=source, catalog=catalog, probe_earthdata=True)

    # source_name property and probe path with no granules warning branch.
    assert isinstance(p.source_name, str)
    state = p.get_prior((0.0, 0.0, 2000.0, 2000.0), "EPSG:4326", datetime(2024, 1, 1), 1000.0)
    assert state.aot.shape == (2, 2)
    assert source.calls and source.calls[0]["short_name"].startswith("SN:")

        (VNP19AODProvider, "vnp19_aod"),
    # _grid validation branch
    with pytest.raises(ValueError, match="resolution must be > 0"):
        p._grid((0.0, 0.0, 1.0, 1.0), 0.0)


def test_brdf_earthaccess_provider_branches():
    src = _StubSource(granules=[])
    cat = _StubCatalog("SN")
    p = MCD43EarthAccessProvider(source=src, catalog=cat, probe_earthdata=True)
    assert p.source_name == "MCD43"

    weights = p.get_brdf_parameters(
        bounds=(0.0, 0.0, 2000.0, 2000.0),
        crs="EPSG:4326",
        obs_time=datetime(2024, 1, 1),
        target_resolution=1000.0,
        bands=[1, 2],
    )
    assert weights.f0.shape[0] == 2

    with pytest.raises(ValueError, match="target_resolution must be > 0"):
        p._grid((0.0, 0.0, 1.0, 1.0), 0.0)

    with pytest.raises(ValueError, match="bands must be a non-empty"):
        p._default_weights((0.0, 0.0, 1.0, 1.0), 1.0, [])

    v = VNP43EarthAccessProvider(source=src, catalog=cat, probe_earthdata=False)
    assert v.source_name == "VNP43"


def test_cams_download_and_explicit_path_branches(monkeypatch, tmp_path: Path):
    p = CAMSProvider(tmp_path)

    # Unsupported explicit file extension branch.
    bad = tmp_path / "cams.txt"
    bad.write_text("x")
    assert p._load_from_explicit_path(bad) is None

    # ImportError branch in _download_cams_file.
    monkeypatch.setitem(sys.modules, "cdsapi", None)
    if "cdsapi" in sys.modules:
        del sys.modules["cdsapi"]
    m = MCD19EarthAccessProvider(source=src, catalog=cat, probe_earthdata=False)
    assert m.source_name == "MCD19"
    assert p._download_cams_file(datetime(2024, 1, 1)) is None

    # Success branch with auth credentials path.
    class _FakeRequest:
        def __init__(self):
            self.downloaded = None

        def download(self, path):
            self.downloaded = path

    class _FakeClient:
        def __init__(self, **kwargs):
            self.kwargs = kwargs
            self.req = _FakeRequest()

        def retrieve(self, dataset, request):
            self.dataset = dataset
            self.request = request
            return self.req

    fake_mod = SimpleNamespace(Client=_FakeClient)
    monkeypatch.setitem(sys.modules, "cdsapi", fake_mod)

    auth = CredentialManager()
    auth.set_credentials("cds", key="abc:123")
    p_auth = CAMSProvider(tmp_path, auth=auth)
    out = p_auth._download_cams_file(datetime(2024, 1, 2))
    assert out is not None
    assert out.name == "CAMS_2024-01-02.nc"

    # Failure branch in retrieve/download.
    class _BoomClient:
        def __init__(self, **kwargs):
            pass

        def retrieve(self, dataset, request):
            raise RuntimeError("boom")

    monkeypatch.setitem(sys.modules, "cdsapi", SimpleNamespace(Client=_BoomClient))
    assert p_auth._download_cams_file(datetime(2024, 1, 3)) is None


def test_cams_download_missing_credentials_branch(monkeypatch, tmp_path: Path):
    monkeypatch.setitem(sys.modules, "cdsapi", SimpleNamespace(Client=lambda **_kwargs: None))
    monkeypatch.delenv("CDSAPI_KEY", raising=False)

    class _HomePath(type(Path())):
        @classmethod
        def home(cls):
            return tmp_path

    monkeypatch.setattr("siac.priors.atmospheric.cams.Path", _HomePath)

    p = CAMSProvider(tmp_path, auth=CredentialManager())
    assert p._download_cams_file(datetime(2024, 1, 4)) is None


def test_cams_download_uses_cache_dir_for_remote_source(monkeypatch, tmp_path: Path):
    captured: dict[str, object] = {}

    class _FakeReq:
        def download(self, path):
            captured["download_path"] = path

    class _FakeClient:
        def __init__(self, **kwargs):
            captured["kwargs"] = kwargs

        def retrieve(self, dataset, request):
            captured["dataset"] = dataset
            captured["request"] = request
            return _FakeReq()

    monkeypatch.setitem(sys.modules, "cdsapi", SimpleNamespace(Client=_FakeClient))

    auth = CredentialManager()
    auth.set_credentials("cds", key="abc:123")
    provider = CAMSProvider(
        "https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/",
        auth=auth,
        cache_dir=tmp_path / "cams-cache",
    )
    out = provider._download_cams_file(datetime(2024, 1, 5))
    assert out == tmp_path / "cams-cache" / "CAMS_2024-01-05.nc"
    assert captured["download_path"] == str(out)


def test_cams_load_data_download_missing_path(monkeypatch, tmp_path: Path):
    p = CAMSProvider(tmp_path / "missing_dir", download_missing=True)

    downloaded = tmp_path / "missing_dir" / "CAMS_2024-01-01.nc"

    monkeypatch.setattr(p, "_download_cams_file", lambda _t: downloaded)
    monkeypatch.setattr(p, "_load_from_explicit_path", lambda _path: "dataset")

    assert p._load_cams_data(datetime(2024, 1, 1)) == "dataset"
