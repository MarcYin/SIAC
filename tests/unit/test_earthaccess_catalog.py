"""Unit tests for Earthaccess catalog and validation logic."""

from __future__ import annotations

import pytest

from siac.core.exceptions import DataNotFoundError
from siac.io.earthaccess_catalog import EarthAccessCatalog, EarthaccessProduct


class _FakeSource:
    def __init__(self, datasets):
        self.datasets = datasets
        self.calls = []

    def search_datasets(self, **kwargs):
        self.calls.append(kwargs)
        return list(self.datasets)


class TestEarthAccessCatalog:
    def test_default_keys_include_core_products(self):
        catalog = EarthAccessCatalog(source=_FakeSource([]))
        keys = set(catalog.keys())

        assert "landsat" in keys
        assert "mcd19_aod" in keys
        assert "mcd43_brdf" in keys
        assert "vnp43_brdf" in keys

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
