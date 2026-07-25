from __future__ import annotations

from tools.aeronet_validation.analyze_low_cloud_aod_failures import (
    _hit_value,
    _mechanism,
    _severity_bin,
    _truth_bin,
)


def _record(
    *,
    retrieved: float,
    curve_min: float,
    bands: tuple[float, float, float],
) -> dict:
    return {
        "matchup_id": "case",
        "truth": 0.5,
        "retrieved": retrieved,
        "solver": {
            "surface_cost_curve_min_aot": curve_min,
            "surface_band_B02_argmin_aot": bands[0],
            "surface_band_B03_argmin_aot": bands[1],
            "surface_band_B04_argmin_aot": bands[2],
            "surface_band_argmin_spread": max(bands) - min(bands),
        },
    }


def test_mechanism_priority_is_mutually_exclusive() -> None:
    assert _mechanism(_record(retrieved=0.2, curve_min=0.5, bands=(0.1, 0.2, 0.3)), 0.2) == "A"
    assert _mechanism(_record(retrieved=0.2, curve_min=0.2, bands=(0.1, 0.2, 0.3)), 0.2) == "B"
    assert _mechanism(_record(retrieved=0.2, curve_min=0.2, bands=(0.2, 0.5, 0.6)), 0.5) == "C"
    assert _mechanism(_record(retrieved=0.2, curve_min=0.2, bands=(0.1, 0.2, 0.3)), 0.5) == "D"
    assert _mechanism(_record(retrieved=0.8, curve_min=0.8, bands=(0.7, 0.75, 0.8)), 0.5) == "E"
    assert _mechanism(_record(retrieved=0.2, curve_min=0.2, bands=(0.40, 0.45, 0.50)), 0.5) == "F"
    assert _mechanism(_record(retrieved=0.8, curve_min=0.8, bands=(0.50, 0.55, 0.60)), 0.5) == "G"


def test_expected_error_boundary_is_inclusive() -> None:
    assert _hit_value(0.5 + 0.05 + 0.15 * 0.5, 0.5)
    assert not _hit_value(0.7, 0.5)


def test_failure_bins_have_stable_boundaries() -> None:
    assert _truth_bin(0.1) == "0.1-0.2"
    assert _truth_bin(1.0) == ">=1.0"
    assert _severity_bin(1.25) == "1.25-1.50 x EE"
    assert _severity_bin(3.0) == ">=3.00 x EE"
