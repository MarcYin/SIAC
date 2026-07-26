"""CIBR water-vapour retrieval: LUT loading, azimuth folding, inversion, masks.

The look-up table is vendored, so these tests assert against the real one: that
it is the complete regular grid the interpolation assumes, that interpolation at
a node reproduces the node exactly, and that the inversion round-trips a known
water-vapour column back out of a synthesised band ratio.
"""

from __future__ import annotations

import numpy as np
import pytest

from siac.algorithms.water_vapour import (
    DEFAULT_MIN_B8A,
    LUT_AXIS_NAMES,
    load_lut,
    relative_azimuth_deg,
    retrieve_water_vapour,
)

_BAND_RATIO = 955.19 / 813.04


def _synthesise_b09(water_vapour_cm: float, b8a: float, b0: float, b1: float) -> float:
    """TOA B09 that inverts back to ``water_vapour_cm`` for these coefficients."""
    log_cibr = (np.sqrt(water_vapour_cm) - b1) / b0
    return float(_BAND_RATIO * b8a / (10.0**log_cibr))


# --------------------------------------------------------------------------- #
# the look-up table
# --------------------------------------------------------------------------- #
def test_lut_is_the_complete_regular_grid_the_interpolation_assumes() -> None:
    axes, coefficients = load_lut()

    assert len(axes) == len(LUT_AXIS_NAMES) == 5
    assert [axis.size for axis in axes] == [7, 5, 5, 5, 7]
    assert coefficients.shape == (7, 5, 5, 5, 7, 2)
    # 7*5*5*5*7 == 6125 rows, every node filled -- no holes to interpolate over.
    assert np.isfinite(coefficients).all()
    np.testing.assert_allclose(axes[0], [0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8])
    np.testing.assert_allclose(axes[1], [0.0, 15.0, 30.0, 45.0, 60.0])
    np.testing.assert_allclose(axes[2], [0.0, 5.0, 10.0, 15.0, 20.0])
    np.testing.assert_allclose(axes[3], [0.0, 45.0, 90.0, 125.0, 180.0])
    np.testing.assert_allclose(axes[4], [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0])


def test_lut_loads_from_the_installed_package_not_the_source_tree() -> None:
    # resources.files() resolves through the import system, so this is the same
    # path a wheel install would take.
    from importlib import resources

    handle = resources.files("siac.data").joinpath("water_vapour_lut_CIBR_9_8A.txt")
    assert handle.is_file()
    text = handle.read_text(encoding="utf-8")
    assert "Brockmann Consult" in text  # provenance header preserved
    assert "surface_reflectance;sun_zenith" in text


def test_first_lut_row_is_reproduced_exactly_at_its_own_node() -> None:
    # Row 1 of the table: 0.05;0;0;0;0.0;-2.74805004884;-0.140468062016
    axes, coefficients = load_lut()
    np.testing.assert_allclose(coefficients[0, 0, 0, 0, 0], [-2.74805004884, -0.140468062016])


# --------------------------------------------------------------------------- #
# relative azimuth
# --------------------------------------------------------------------------- #
def test_relative_azimuth_folds_rather_than_wrapping() -> None:
    # The table stops at 180 because it is symmetric about the principal plane,
    # so 200 degrees is the same geometry as 160 -- NOT as 20, which is what the
    # research code's `% 180` produced.
    np.testing.assert_allclose(relative_azimuth_deg(200.0, 0.0), 160.0)
    np.testing.assert_allclose(relative_azimuth_deg(0.0, 200.0), 160.0)
    np.testing.assert_allclose(relative_azimuth_deg(190.0, 10.0), 180.0)
    np.testing.assert_allclose(relative_azimuth_deg(90.0, 0.0), 90.0)
    np.testing.assert_allclose(relative_azimuth_deg(-30.0, 0.0), 30.0)
    folded = relative_azimuth_deg(np.array([10.0, 350.0]), np.array([0.0, 0.0]))
    np.testing.assert_allclose(folded, [10.0, 10.0])
    assert np.all((folded >= 0.0) & (folded <= 180.0))


# --------------------------------------------------------------------------- #
# retrieval
# --------------------------------------------------------------------------- #
def test_retrieval_inverts_a_synthesised_column_at_a_lut_node() -> None:
    # At an exact grid node the interpolator must return that node's own
    # coefficients, so a band ratio built from them inverts back to the column.
    _axes, coefficients = load_lut()
    b0, b1 = coefficients[1, 1, 1, 1, 1]  # rho=0.1, sza=15, vza=5, raa=45, alt=0.5
    truth = 2.4
    b09 = _synthesise_b09(truth, 0.1, float(b0), float(b1))

    result = retrieve_water_vapour(
        toa_b09=np.full((3, 4), b09),
        toa_b8a=np.full((3, 4), 0.1),
        sza_deg=15.0,
        vza_deg=5.0,
        raa_deg=45.0,
        elevation_km=0.5,
    )

    assert result.valid.all()
    np.testing.assert_allclose(result.water_vapour_cm, truth, rtol=1e-4)
    np.testing.assert_allclose(result.uncertainty_cm, 0.1)
    assert result.masked_fraction == 0.0


def test_retrieval_is_vectorised_over_a_per_pixel_geometry_and_terrain() -> None:
    _axes, coefficients = load_lut()
    b0, b1 = coefficients[1, 1, 1, 1, 1]
    b09 = _synthesise_b09(2.4, 0.1, float(b0), float(b1))
    elevation = np.full((2, 2), 0.5)
    elevation[1, 1] = 1.0  # a different LUT node -> a different column

    result = retrieve_water_vapour(
        toa_b09=np.full((2, 2), b09),
        toa_b8a=np.full((2, 2), 0.1),
        sza_deg=np.full((2, 2), 15.0),
        vza_deg=np.full((2, 2), 5.0),
        raa_deg=np.full((2, 2), 45.0),
        elevation_km=elevation,
    )

    assert result.valid.all()
    np.testing.assert_allclose(result.water_vapour_cm[0, 0], 2.4, rtol=1e-4)
    assert result.water_vapour_cm[1, 1] != pytest.approx(result.water_vapour_cm[0, 0])


def test_dark_targets_and_missing_bands_are_masked_and_filled_from_the_scene() -> None:
    _axes, coefficients = load_lut()
    b0, b1 = coefficients[1, 1, 1, 1, 1]
    b09 = _synthesise_b09(2.4, 0.1, float(b0), float(b1))
    b8a = np.full((2, 3), 0.1)
    b8a[0, 0] = DEFAULT_MIN_B8A * 0.5  # too dark for a stable ratio
    band09 = np.full((2, 3), b09)
    band09[0, 1] = np.nan  # missing observation

    result = retrieve_water_vapour(
        toa_b09=band09,
        toa_b8a=b8a,
        sza_deg=15.0,
        vza_deg=5.0,
        raa_deg=45.0,
        elevation_km=0.5,
    )

    assert not result.valid[0, 0]
    assert not result.valid[0, 1]
    assert result.valid.sum() == 4
    # Masked pixels take this scene's OWN median, not an assumed climatology.
    np.testing.assert_allclose(result.water_vapour_cm, 2.4, rtol=1e-4)
    assert result.uncertainty_cm[0, 0] == pytest.approx(1.0)
    assert result.uncertainty_cm[1, 1] == pytest.approx(0.1)
    assert result.masked_fraction == pytest.approx(2 / 6)


def test_masked_pixels_stay_nan_when_fill_is_off() -> None:
    result = retrieve_water_vapour(
        toa_b09=np.array([[0.05, np.nan]]),
        toa_b8a=np.array([[0.2, 0.2]]),
        sza_deg=15.0,
        vza_deg=5.0,
        raa_deg=45.0,
        elevation_km=0.5,
        fill=False,
    )

    assert np.isnan(result.water_vapour_cm[0, 1])


def test_a_cloud_mask_is_honoured() -> None:
    _axes, coefficients = load_lut()
    b0, b1 = coefficients[1, 1, 1, 1, 1]
    b09 = _synthesise_b09(2.4, 0.1, float(b0), float(b1))
    cloud = np.array([[True, False]])

    result = retrieve_water_vapour(
        toa_b09=np.full((1, 2), b09),
        toa_b8a=np.full((1, 2), 0.1),
        sza_deg=15.0,
        vza_deg=5.0,
        raa_deg=45.0,
        elevation_km=0.5,
        cloud_mask=cloud,
    )

    assert list(result.valid.ravel()) == [False, True]
    assert result.uncertainty_cm[0, 0] == pytest.approx(1.0)


def test_an_implausible_inversion_is_masked() -> None:
    # A near-zero B09 against a bright continuum drives the column far past any
    # real atmosphere; that is a failed inversion, not an unusual sky.
    result = retrieve_water_vapour(
        toa_b09=np.full((2, 2), 0.001),
        toa_b8a=np.full((2, 2), 0.4),
        sza_deg=15.0,
        vza_deg=5.0,
        raa_deg=45.0,
        elevation_km=0.5,
    )

    assert not result.valid.any()
    assert result.masked_fraction == 1.0
    # Nothing valid to fill from, so the column stays NaN rather than invented.
    assert np.isnan(result.water_vapour_cm).all()


def test_queries_outside_the_table_are_clamped_not_extrapolated() -> None:
    _axes, coefficients = load_lut()
    b0, b1 = coefficients[-1, -1, -1, -1, -1]  # rho=0.8, sza=60, vza=20, raa=180, alt=4
    b09 = _synthesise_b09(2.0, 0.9, float(b0), float(b1))

    # rho 0.9 > 0.8, sza 75 > 60, vza 35 > 20, altitude 6 > 4: all clamped to
    # the last node, so this reproduces that node's own coefficients.
    result = retrieve_water_vapour(
        toa_b09=np.full((1, 1), b09),
        toa_b8a=np.full((1, 1), 0.9),
        sza_deg=75.0,
        vza_deg=35.0,
        raa_deg=180.0,
        elevation_km=6.0,
    )

    np.testing.assert_allclose(result.water_vapour_cm, 2.0, rtol=1e-4)


def test_mismatched_band_grids_are_rejected() -> None:
    with pytest.raises(ValueError, match="same grid"):
        retrieve_water_vapour(
            toa_b09=np.zeros((2, 2)),
            toa_b8a=np.zeros((3, 3)),
            sza_deg=15.0,
            vza_deg=5.0,
            raa_deg=45.0,
            elevation_km=0.5,
        )


def test_interpolation_keeps_the_full_reflectance_dependence() -> None:
    """Guards the choice of interpolator, measured between grid nodes.

    A 32-neighbour distance-weighted KNN on unnormalised features -- what both
    research implementations used -- puts reflectance (0.05-0.8) on the same
    axis as sun zenith (0-60), so the neighbourhood is dominated by the angles
    and retains only ~13% of this span. That flattening is what the research
    ``* 0.716`` was patching; multilinear interpolation of the grid needs no
    such patch, and measured against CAMS it is the least biased of the three.

    Queried strictly between nodes: exactly *on* a node every interpolator
    returns that node's own value, so the difference only shows where
    interpolation actually happens.
    """
    from siac.algorithms.water_vapour import _coefficients

    dark = np.array([[0.15, 22.0, 7.0, 60.0, 0.8]])
    bright = np.array([[0.55, 22.0, 7.0, 60.0, 0.8]])

    span = abs(float(_coefficients(bright)[0, 0]) - float(_coefficients(dark)[0, 0]))

    # Measured: 0.166 across this reflectance range (raw-feature KNN gave 0.021).
    assert span == pytest.approx(0.166, abs=0.01)
