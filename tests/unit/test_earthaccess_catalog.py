"""Unit tests for Earthaccess catalog and validation logic."""

from __future__ import annotations

import pytest

from siac.adapters.data.earthaccess_catalog import (
    EarthAccessCatalog,
    EarthaccessProduct,
    ProductValidationResult,
)
from siac.errors import DataNotFoundError


class _FakeSource:
    def __init__(self, datasets):
        self.datasets = datasets
        self.calls = []

    def search_datasets(self, **kwargs):
        self.calls.append(kwargs)
        return list(self.datasets)


class _ObjDataset:
    def __init__(self, *, short_name=None, umm=None):
        self.short_name = short_name
        self.umm = umm


class TestEarthAccessCatalog:
    def test_default_keys_include_core_products(self):
        catalog = EarthAccessCatalog(source=_FakeSource([]))
        keys = set(catalog.keys())

        assert "landsat" in keys
        assert "mcd19_aod" in keys
        assert "vnp19_aod" in keys
        assert "mcd43_brdf" in keys
        assert "mcd19_brdf" in keys
        assert "vnp43_brdf" in keys

    def test_default_vnp43_points_to_moderate_product(self):
        catalog = EarthAccessCatalog(source=_FakeSource([]))
        assert catalog.get_product("vnp43_brdf").short_name == "VNP43MA1"

    def test_validate_product_infers_short_name(self):
        product = EarthaccessProduct(
            key="test",
            description="test product",
            keyword="TEST",
            short_name=None,
        )
        source = _FakeSource([{"umm": {"ShortName": "TEST_SN"}}])
        catalog = EarthAccessCatalog(source=source, products={"test": product})

        result = catalog.validate_product("test")

        assert result.selected_short_name == "TEST_SN"
        assert catalog.get_product("test").short_name == "TEST_SN"

    def test_validate_product_uses_configured_short_name(self):
        product = EarthaccessProduct(
            key="mcd43",
            description="mcd43",
            keyword="MCD43",
            short_name="MCD43A1",
        )
        source = _FakeSource([{"short_name": "MCD43A1"}])
        catalog = EarthAccessCatalog(source=source, products={"mcd43": product})

        result = catalog.validate_product("mcd43")

        assert result.selected_short_name == "MCD43A1"

    def test_validate_product_no_results_raises(self):
        product = EarthaccessProduct(
            key="missing",
            description="missing",
            keyword="NOT_REAL",
        )
        catalog = EarthAccessCatalog(source=_FakeSource([]), products={"missing": product})

        with pytest.raises(DataNotFoundError, match="No Earthdata datasets found"):
            _ = catalog.validate_product("missing")

    def test_set_override(self):
        product = EarthaccessProduct(
            key="p",
            description="desc",
            keyword="k",
        )
        catalog = EarthAccessCatalog(source=_FakeSource([]), products={"p": product})

        catalog.set_override("p", short_name="X1", provider="LPDAAC_ECS")
        updated = catalog.get_product("p")
        assert updated.short_name == "X1"
        assert updated.provider == "LPDAAC_ECS"

    def test_get_product_unknown_key_raises(self):
        catalog = EarthAccessCatalog(source=_FakeSource([]), products={})

        with pytest.raises(KeyError, match="Unknown Earthaccess product key"):
            _ = catalog.get_product("not-a-key")

    def test_set_override_provider_only(self):
        product = EarthaccessProduct(key="p", description="desc", keyword="k", short_name="S")
        catalog = EarthAccessCatalog(source=_FakeSource([]), products={"p": product})

        catalog.set_override("p", provider="POCLOUD")
        updated = catalog.get_product("p")
        assert updated.short_name == "S"
        assert updated.provider == "POCLOUD"

    def test_validate_product_infers_short_name_from_object_attribute(self):
        product = EarthaccessProduct(key="obj", description="obj", keyword="OBJ")
        source = _FakeSource([_ObjDataset(short_name="OBJ_SN")])
        catalog = EarthAccessCatalog(source=source, products={"obj": product})

        result = catalog.validate_product("obj")
        assert result.selected_short_name == "OBJ_SN"
        assert catalog.get_product("obj").short_name == "OBJ_SN"

    def test_validate_product_infers_short_name_from_object_umm(self):
        product = EarthaccessProduct(key="obj2", description="obj2", keyword="OBJ2")
        source = _FakeSource([_ObjDataset(umm={"ShortName": "OBJ_UMM_SN"})])
        catalog = EarthAccessCatalog(source=source, products={"obj2": product})

        result = catalog.validate_product("obj2")
        assert result.selected_short_name == "OBJ_UMM_SN"

    def test_validate_product_configured_short_name_missing_raises(self):
        product = EarthaccessProduct(
            key="mismatch",
            description="mismatch",
            keyword="MIS",
            short_name="EXPECTED",
        )
        source = _FakeSource([{"short_name": "OTHER"}])
        catalog = EarthAccessCatalog(source=source, products={"mismatch": product})

        with pytest.raises(DataNotFoundError, match="Configured short_name"):
            _ = catalog.validate_product("mismatch")

    def test_validate_product_no_extractable_short_name_raises(self):
        product = EarthaccessProduct(key="none", description="none", keyword="NONE")
        source = _FakeSource([_ObjDataset()])
        catalog = EarthAccessCatalog(source=source, products={"none": product})

        with pytest.raises(DataNotFoundError, match="Could not infer dataset short_name"):
            _ = catalog.validate_product("none")

    def test_validate_all_uses_default_and_explicit_keys(self, monkeypatch):
        products = {
            "a": EarthaccessProduct(key="a", description="A", keyword="A"),
            "b": EarthaccessProduct(key="b", description="B", keyword="B"),
        }
        catalog = EarthAccessCatalog(source=_FakeSource([]), products=products)
        calls = []

        def _fake_validate(key, *, count=20):
            calls.append((key, count))
            return ProductValidationResult(
                key=key,
                dataset_count=1,
                discovered_short_names=("X",),
                selected_short_name="X",
            )

        monkeypatch.setattr(catalog, "validate_product", _fake_validate)

        all_result = catalog.validate_all(count=7)
        assert set(all_result.keys()) == {"a", "b"}
        assert ("a", 7) in calls
        assert ("b", 7) in calls

        calls.clear()
        sub_result = catalog.validate_all(keys=["b"], count=3)
        assert set(sub_result.keys()) == {"b"}
        assert calls == [("b", 3)]

    def test_resolve_short_name_uses_cached_or_validates(self, monkeypatch):
        cached_product = EarthaccessProduct(
            key="cached",
            description="cached",
            keyword="CACHED",
            short_name="READY",
        )
        cached_catalog = EarthAccessCatalog(source=_FakeSource([]), products={"cached": cached_product})

        marker = {"called": False}

        def _never_called(key, *, count=20):
            marker["called"] = True
            raise AssertionError("validate_product should not be called")

        monkeypatch.setattr(cached_catalog, "validate_product", _never_called)
        assert cached_catalog.resolve_short_name("cached") == "READY"
        assert marker["called"] is False

        uncached_product = EarthaccessProduct(
            key="uncached",
            description="uncached",
            keyword="UNCACHED",
            short_name=None,
        )
        uncached_catalog = EarthAccessCatalog(source=_FakeSource([]), products={"uncached": uncached_product})
        monkeypatch.setattr(
            uncached_catalog,
            "validate_product",
            lambda key, *, _count=20: ProductValidationResult(
                key=key,
                dataset_count=2,
                discovered_short_names=("DISCOVERED",),
                selected_short_name="DISCOVERED",
            ),
        )
        assert uncached_catalog.resolve_short_name("uncached") == "DISCOVERED"
