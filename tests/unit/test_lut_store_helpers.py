"""Focused unit tests for LUT store helper branches."""

from __future__ import annotations

import pytest

import siac.algorithms.rt.lut.store as lut_store


def test_coerce_bool_and_split_storage_options_cover_reference_flags(tmp_path) -> None:
    assert lut_store._coerce_bool(None, default=True) is True  # noqa: SLF001
    assert lut_store._coerce_bool(True, default=False) is True  # noqa: SLF001
    assert lut_store._coerce_bool(" yes ", default=False) is True  # noqa: SLF001
    assert lut_store._coerce_bool("OFF", default=True) is False  # noqa: SLF001
    assert lut_store._coerce_bool(0, default=True) is False  # noqa: SLF001

    with pytest.raises(TypeError, match="use_reference"):
        lut_store._split_storage_options({"use_reference": True})  # noqa: SLF001

    options, reference_options = lut_store._split_storage_options(  # noqa: SLF001
        {
            "anon": True,
            "reference_refresh": "true",
            "reference_json": str(tmp_path / "refs.json"),
            "reference_cache_dir": str(tmp_path / "cache"),
        }
    )
    assert options == {"anon": True}
    assert reference_options.refresh is True
    assert reference_options.reference_json == tmp_path / "refs.json"
    assert reference_options.cache_dir == tmp_path / "cache"


def test_reference_key_and_document_helpers_cover_skip_and_error_paths(tmp_path) -> None:
    assert lut_store._normalize_mapper_root("") == ""  # noqa: SLF001
    assert lut_store._normalize_mapper_root(".") == ""  # noqa: SLF001
    assert lut_store._normalize_mapper_root("/lut.zarr/") == "lut.zarr"  # noqa: SLF001

    assert lut_store._relative_reference_key("", "lut.zarr") is None  # noqa: SLF001
    assert lut_store._relative_reference_key("lut.zarr", "lut.zarr") is None  # noqa: SLF001
    assert lut_store._relative_reference_key("other/.zgroup", "lut.zarr") is None  # noqa: SLF001
    assert lut_store._relative_reference_key("lut.zarr/.zgroup", "lut.zarr") == ".zgroup"  # noqa: SLF001
    assert lut_store._relative_reference_key("arr/0", "") == "arr/0"  # noqa: SLF001

    with pytest.raises(ValueError, match="does not expose indexed file metadata"):
        lut_store._build_reference_document("https://example.com/lut.zip", object())  # noqa: SLF001

    class _NegativeMapper:
        fs = type("FS", (), {"_files": {"lut.zarr/.zgroup": {"offset": -1, "size": 3}}})()
        root = "lut.zarr"

    with pytest.raises(ValueError, match="Invalid byte range"):
        lut_store._build_reference_document("https://example.com/lut.zip", _NegativeMapper())  # noqa: SLF001

    class _EmptyMapper:
        fs = type("FS", (), {"_files": {"lut.zarr": {"children": ["lut.zarr/.zgroup"]}}})()
        root = "lut.zarr"

    with pytest.raises(ValueError, match="No files available"):
        lut_store._build_reference_document("https://example.com/lut.zip", _EmptyMapper())  # noqa: SLF001

    assert lut_store._is_valid_reference_document({"refs": {}}) is False  # noqa: SLF001
    assert lut_store._is_valid_reference_document({"refs": {"arr/0": ["x", 1, 2]}}) is False  # noqa: SLF001
    assert lut_store._is_valid_reference_document({"refs": {".zgroup": ["x", 1, 2]}}) is True  # noqa: SLF001

    explicit = lut_store._reference_json_path(  # noqa: SLF001
        "https://example.com/lut.zarr.zip",
        lut_store._ReferenceOptions(
            refresh=False,
            reference_json=tmp_path / "explicit.json",
            cache_dir=None,
        ),
    )
    assert explicit == tmp_path / "explicit.json"


def test_reference_remote_and_open_mapper_cover_protocol_variants(monkeypatch: pytest.MonkeyPatch, tmp_path) -> None:
    with pytest.raises(TypeError, match="headers"):
        lut_store._reference_remote("https://example.com/lut.zip", {"headers": "bad"})  # noqa: SLF001

    assert lut_store._reference_remote("lut.zarr.zip", {}) == (None, {})  # noqa: SLF001

    protocol, options = lut_store._reference_remote("s3://bucket/lut.zip", {"anon": True})  # noqa: SLF001
    assert protocol == "s3"
    assert options["anon"] is True
    assert options["asynchronous"] is True

    reference_path = tmp_path / "refs.json"
    reference_path.write_text('{"version":1,"refs":{".zgroup":["x",1,2]}}', encoding="utf-8")

    captured: dict[str, object] = {}

    monkeypatch.setattr(
        "fsspec.get_mapper",
        lambda path, **kwargs: captured.update({"path": path, "kwargs": kwargs}) or {"ok": True},
    )

    out = lut_store._open_reference_mapper(  # noqa: SLF001
        "https://example.com/lut.zip",
        {"timeout": 7.0, "headers": {"Authorization": "Bearer token"}},
        reference_path,
    )
    assert out == {"ok": True}
    assert captured["path"] == "reference://"
    assert captured["kwargs"]["remote_protocol"] == "https"
    assert captured["kwargs"]["remote_options"] == {
        "timeout": 7.0,
        "asynchronous": True,
        "headers": {"Authorization": "Bearer token"},
    }

    with pytest.raises(ValueError, match="Invalid reference JSON"):
        lut_store._open_reference_mapper(  # noqa: SLF001
            "https://example.com/lut.zip",
            {},
            reference_path,
            document={"refs": []},
        )


def test_remote_reference_mapper_refresh_and_read_failure_rebuild(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    build_calls = {"n": 0}
    open_calls = {"n": 0}
    reference_path = tmp_path / "cached.refs.json"
    reference_path.write_text("not json", encoding="utf-8")

    class _Mapper:
        fs = type("FS", (), {"_files": {"lut.zarr/.zgroup": {"offset": 1, "size": 2}}})()
        root = "lut.zarr"

    def _fake_build_readonly_zip_mapper(_path: str, _options: dict[str, object]) -> _Mapper:
        build_calls["n"] += 1
        return _Mapper()

    monkeypatch.setattr(lut_store, "build_readonly_zip_mapper", _fake_build_readonly_zip_mapper)
    monkeypatch.setattr(
        lut_store,
        "_open_reference_mapper",
        lambda path, storage_options, reference_path, document=None: open_calls.__setitem__("n", open_calls["n"] + 1) or {"document": document},  # noqa: ARG005
    )

    out = lut_store._build_remote_zip_reference_mapper(  # noqa: SLF001
        "https://example.com/lut.zarr.zip",
        {"timeout": 5.0},
        lut_store._ReferenceOptions(
            refresh=False,
            reference_json=reference_path,
            cache_dir=None,
        ),
    )

    assert out["document"]["refs"][".zgroup"] == ["https://example.com/lut.zarr.zip", 1, 2]
    assert build_calls["n"] == 1
    assert open_calls["n"] == 1

    build_calls["n"] = 0
    out_refresh = lut_store._build_remote_zip_reference_mapper(  # noqa: SLF001
        "https://example.com/lut.zarr.zip",
        {"timeout": 5.0},
        lut_store._ReferenceOptions(
            refresh=True,
            reference_json=reference_path,
            cache_dir=None,
        ),
    )
    assert out_refresh["document"]["refs"][".zgroup"] == ["https://example.com/lut.zarr.zip", 1, 2]
    assert build_calls["n"] == 1
