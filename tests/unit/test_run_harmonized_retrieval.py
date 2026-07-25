import pytest
from tools.aeronet_validation.run_harmonized_retrieval import _band_names


def test_direct_history_bands_default_includes_green() -> None:
    assert _band_names(None) == ("B01", "B02", "B03", "B04")


def test_direct_history_bands_are_deduplicated_in_order() -> None:
    assert _band_names("B02, B03,B02,B04") == ("B02", "B03", "B04")


def test_direct_history_bands_reject_unknown_names() -> None:
    with pytest.raises(ValueError, match="unsupported direct-history bands"):
        _band_names("B02,B99")
