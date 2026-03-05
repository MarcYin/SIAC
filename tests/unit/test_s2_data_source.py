"""
Unit tests for S2 data access types and helper functions.
"""

from datetime import date, datetime

import pytest

from siac.io.s2_data_source import (
    S2Product,
    S2Query,
    deduplicate_products,
    select_best_product,
)

# ── S2Query ───────────────────────────────────────────────────────────

class TestS2Query:
    def test_from_product_id(self):
        q = S2Query.from_product_id(
            "S2B_MSIL1C_20210801T103629_N0500_R008_T31UDQ_20230731T120241"
        )
        assert q.product_id is not None
        assert "MSIL1C" in q.product_id

    def test_from_tile_date(self):
        q = S2Query.from_tile_date("T31UDQ_20210801")
        assert q.mgrs_tile == "31UDQ"
        assert q.date == date(2021, 8, 1)

    def test_from_tile_date_no_prefix(self):
        q = S2Query.from_tile_date("31UDQ_20210801")
        assert q.mgrs_tile == "31UDQ"

    def test_from_tile_date_invalid(self):
        with pytest.raises(ValueError, match="Cannot parse"):
            S2Query.from_tile_date("invalid_string")

    def test_validate_no_constraint_raises(self):
        q = S2Query()
        with pytest.raises(ValueError, match="at least one spatial constraint"):
            q.validate()

    def test_validate_with_product_id(self):
        q = S2Query(product_id="some_id")
        q.validate()  # should not raise

    def test_validate_with_bbox(self):
        q = S2Query(bbox=(1.5, 50.0, 2.5, 51.0))
        q.validate()  # should not raise

    def test_validate_with_tile(self):
        q = S2Query(mgrs_tile="31UDQ")
        q.validate()  # should not raise


# ── S2Product ─────────────────────────────────────────────────────────

class TestS2Product:
    def _make_product(self, baseline="N0500", sensing_date=None):
        return S2Product(
            product_id=f"S2B_MSIL1C_20210801__{baseline}",
            mgrs_tile="31UDQ",
            sensing_date=sensing_date or datetime(2021, 8, 1, 10, 30),
            processing_baseline=baseline,
            cloud_cover=15.0,
            satellite="S2B",
            orbit_number=8,
            source_url="https://example.com/product",
        )

    def test_baseline_number(self):
        p = self._make_product("N0500")
        assert p.baseline_number == 500

    def test_baseline_number_n0301(self):
        p = self._make_product("N0301")
        assert p.baseline_number == 301


# ── Deduplication ─────────────────────────────────────────────────────

class TestDeduplication:
    def _make_product(self, baseline, tile="31UDQ", dt=None):
        dt = dt or datetime(2021, 8, 1, 10, 30)
        return S2Product(
            product_id=f"S2B__{baseline}__{tile}",
            mgrs_tile=tile,
            sensing_date=dt,
            processing_baseline=baseline,
            cloud_cover=10.0,
            satellite="S2B",
            orbit_number=8,
            source_url="https://example.com",
        )

    def test_dedup_keeps_newest_baseline(self):
        products = [
            self._make_product("N0301"),
            self._make_product("N0400"),
            self._make_product("N0500"),
        ]
        result = deduplicate_products(products)
        assert len(result) == 1
        assert result[0].processing_baseline == "N0500"

    def test_dedup_different_tiles(self):
        products = [
            self._make_product("N0500", tile="31UDQ"),
            self._make_product("N0500", tile="32UPU"),
        ]
        result = deduplicate_products(products)
        assert len(result) == 2

    def test_dedup_different_dates(self):
        products = [
            self._make_product("N0500", dt=datetime(2021, 8, 1)),
            self._make_product("N0500", dt=datetime(2021, 8, 5)),
        ]
        result = deduplicate_products(products)
        assert len(result) == 2

    def test_select_best_product(self):
        products = [
            self._make_product("N0301"),
            self._make_product("N0500"),
        ]
        best = select_best_product(products)
        assert best.processing_baseline == "N0500"

    def test_select_best_empty_raises(self):
        with pytest.raises(ValueError, match="No products"):
            select_best_product([])
