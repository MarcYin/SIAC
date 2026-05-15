"""Regression tests for the S2 MTD_TL.xml angle parsers.

REVIEW: a long-standing bug — the parsers built the namespace prefix
from the root element but child elements in S2 MTD_TL are unprefixed,
so every ``findall(f".//{ns}TAG")`` returned nothing and the helpers
silently fell back to the ``DEFAULT_S2_*`` magic-number grids. Goldens
captured from any real scene were therefore characterising a run using
made-up angles. These tests pin the corrected contract using
synthesized fixtures that mirror the namespace layout of a real
Sentinel-2 SAFE.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from textwrap import dedent

import numpy as np

from siac.adapters.satellite.sentinel2 import Sentinel2Preprocessor


def _preprocessor() -> Sentinel2Preprocessor:
    """Instance without invoking __init__ (we only need the parser methods)."""
    return Sentinel2Preprocessor.__new__(Sentinel2Preprocessor)


# A minimal MTD_TL.xml that reproduces the real-world namespace layout:
# the root carries ``xmlns:n1`` but every descendant is unprefixed.
_FIXTURE_XML = dedent(
    """\
    <?xml version="1.0" encoding="UTF-8" standalone="no"?>
    <n1:Level-1C_Tile_ID xmlns:n1="https://psd-15.sentinel2.eo.esa.int/PSD/S2_PDI_Level-1C_Tile_Metadata.xsd">
      <n1:Geometric_Info>
        <Tile_Angles>
          <Sun_Angles_Grid>
            <Zenith>
              <COL_STEP unit="m">5000</COL_STEP>
              <ROW_STEP unit="m">5000</ROW_STEP>
              <Values_List>
                <VALUES>30.0 31.0 32.0</VALUES>
                <VALUES>30.5 31.5 32.5</VALUES>
                <VALUES>31.0 32.0 33.0</VALUES>
              </Values_List>
            </Zenith>
            <Azimuth>
              <COL_STEP unit="m">5000</COL_STEP>
              <ROW_STEP unit="m">5000</ROW_STEP>
              <Values_List>
                <VALUES>120.0 121.0 122.0</VALUES>
                <VALUES>120.5 121.5 122.5</VALUES>
                <VALUES>121.0 122.0 123.0</VALUES>
              </Values_List>
            </Azimuth>
          </Sun_Angles_Grid>
          <Mean_Sun_Angle>
            <ZENITH_ANGLE unit="deg">31.5</ZENITH_ANGLE>
            <AZIMUTH_ANGLE unit="deg">121.5</AZIMUTH_ANGLE>
          </Mean_Sun_Angle>
          <Mean_Viewing_Incidence_Angle_List>
            <Mean_Viewing_Incidence_Angle bandId="0">
              <ZENITH_ANGLE unit="deg">8.0</ZENITH_ANGLE>
              <AZIMUTH_ANGLE unit="deg">280.0</AZIMUTH_ANGLE>
            </Mean_Viewing_Incidence_Angle>
            <Mean_Viewing_Incidence_Angle bandId="1">
              <ZENITH_ANGLE unit="deg">9.0</ZENITH_ANGLE>
              <AZIMUTH_ANGLE unit="deg">290.0</AZIMUTH_ANGLE>
            </Mean_Viewing_Incidence_Angle>
          </Mean_Viewing_Incidence_Angle_List>
        </Tile_Angles>
      </n1:Geometric_Info>
    </n1:Level-1C_Tile_ID>
    """
)


def test_parse_view_angles_reads_unprefixed_children() -> None:
    """The view-angle parser must find unprefixed descendants under an n1: root."""
    root = ET.fromstring(_FIXTURE_XML)
    parsed = _preprocessor()._parse_view_angles(root)

    # The fixture has two bands with VZA 8.0/9.0 and VAA 280/290.
    # Arithmetic mean → VZA=8.5, VAA=285.
    assert np.allclose(parsed["zenith"], 8.5), float(np.mean(parsed["zenith"]))
    assert np.allclose(parsed["azimuth"], 285.0), float(np.mean(parsed["azimuth"]))
    # The grid is the simplified 23x23 uniform shape (full per-detector
    # parsing is a separate roadmap item).
    assert parsed["zenith"].shape == (23, 23)


def test_parse_view_angles_does_not_fall_back_when_xml_has_data() -> None:
    """The DEFAULT_S2_VZA_DEG=5.0 / VAA=100.0 fallback must NOT fire on a
    well-formed MTD_TL.xml. Pre-wave-14 it always did.
    """
    root = ET.fromstring(_FIXTURE_XML)
    parsed = _preprocessor()._parse_view_angles(root)

    # If the fallback fires, VZA would be 5.0deg and VAA 100.0deg.
    assert not np.allclose(parsed["zenith"], 5.0)
    assert not np.allclose(parsed["azimuth"], 100.0)


def test_parse_sun_angles_reads_unprefixed_values_list() -> None:
    """Sun angle grid values must be read from unprefixed Values_List/VALUES."""
    root = ET.fromstring(_FIXTURE_XML)
    parsed = _preprocessor()._parse_sun_angles(root)

    # Zenith grid from the fixture is the 3x3 we ship in the XML.
    assert "zenith" in parsed and "azimuth" in parsed
    expected_zen = np.array([[30.0, 31.0, 32.0], [30.5, 31.5, 32.5], [31.0, 32.0, 33.0]])
    expected_az = np.array(
        [[120.0, 121.0, 122.0], [120.5, 121.5, 122.5], [121.0, 122.0, 123.0]]
    )
    assert np.allclose(parsed["zenith"], expected_zen.astype(np.float32))
    assert np.allclose(parsed["azimuth"], expected_az.astype(np.float32))


def test_parse_sun_angles_does_not_fall_back_to_uniform_grid() -> None:
    """The DEFAULT_S2_ANGLE_GRID_DEG=30.0 uniform fallback must NOT fire when
    Values_List is present and populated.
    """
    root = ET.fromstring(_FIXTURE_XML)
    parsed = _preprocessor()._parse_sun_angles(root)

    # If the fallback fires, the grid is a uniform 30.0 of shape (23, 23).
    zenith = parsed["zenith"]
    assert zenith.shape != (23, 23) or not np.allclose(zenith, 30.0)
