"""Additional coverage tests for Sentinel-2 preprocessor internals."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from siac.satellite.sentinel2 import Sentinel2Preprocessor


def _da(shape: tuple[int, int] = (8, 8), value: float = 1000.0) -> xr.DataArray:
    x = np.linspace(500000.0, 500000.0 + (shape[1] - 1) * 10.0, shape[1])
    y = np.linspace(4500000.0 + (shape[0] - 1) * 10.0, 4500000.0, shape[0])
    da = xr.DataArray(np.full(shape, value, dtype=np.float32), dims=["y", "x"], coords={"y": y, "x": x})
    return da.rio.write_crs("EPSG:32632")


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def _safe_tree(tmp_path: Path, safe_name: str = "S2A_TEST.SAFE") -> Path:
    safe = tmp_path / safe_name
    img_data = safe / "GRANULE" / "L1C_TILE" / "IMG_DATA"
    img_data.mkdir(parents=True, exist_ok=True)
    for band in ["B02", "B04", "B10"]:
        (img_data / f"TILE_{band}.jp2").write_text("x")
    return safe


class TestSentinel2Internals:
    def test_resolve_paths_variants(self, tmp_path: Path):
        p = Sentinel2Preprocessor()

        safe_a = _safe_tree(tmp_path, "S2A_DEMO.SAFE")
        p._resolve_paths(safe_a)
        assert p._granule_path is not None
        assert p.sensor_config.satellite_id == "S2A"

        p2 = Sentinel2Preprocessor()
        safe_b = _safe_tree(tmp_path, "S2B_DEMO.SAFE")
        p2._resolve_paths(safe_b)
        assert p2.sensor_config.satellite_id == "S2B"

        p3 = Sentinel2Preprocessor()
        aws = tmp_path / "aws_s2"
        aws.mkdir()
        (aws / "metadata.xml").write_text("<root/>")
        p3._resolve_paths(aws)
        assert p3._granule_path == aws

        p4 = Sentinel2Preprocessor()
        bad = tmp_path / "bad"
        bad.mkdir()
        with pytest.raises(FileNotFoundError):
            p4._resolve_paths(bad)

    def test_metadata_parsers_and_finders(self, tmp_path: Path):
        p = Sentinel2Preprocessor()
        safe = _safe_tree(tmp_path)
        p._resolve_paths(safe)

        product_xml = """<root>
            <PROCESSING_BASELINE>05.00</PROCESSING_BASELINE>
            <QUANTIFICATION_VALUE>10000</QUANTIFICATION_VALUE>
            <RADIO_ADD_OFFSET band_id="1">-100</RADIO_ADD_OFFSET>
            <RADIO_ADD_OFFSET band_id="2">-50</RADIO_ADD_OFFSET>
        </root>"""
        granule_xml = """<root>
            <SENSING_TIME>2024-01-01T10:00:00.000Z</SENSING_TIME>
            <TILE_ID>T31UDQ</TILE_ID>
        </root>"""
        _write(safe / "MTD_MSIL1C.xml", product_xml)
        _write(p._granule_path / "MTD_TL.xml", granule_xml)

        assert p._find_product_xml(safe) is not None
        assert p._find_granule_xml() is not None

        md_prod = p._parse_product_metadata(ET.fromstring(product_xml))
        assert md_prod["quantification_value"] == 10000.0
        assert md_prod["radiometric_offsets"]["B01"] == -100.0

        md_gran = p._parse_granule_metadata(ET.fromstring(granule_xml))
        assert md_gran["tile_id"] == "T31UDQ"
        assert isinstance(md_gran["observation_time"], datetime)

        # fallback datetime parse branch
        granule_xml_alt = """<root><SENSING_TIME>2024-01-01T10:00:00Z</SENSING_TIME></root>"""
        md_gran_alt = p._parse_granule_metadata(ET.fromstring(granule_xml_alt))
        assert isinstance(md_gran_alt["observation_time"], datetime)

    def test_angle_parsers_and_grid_interpolation(self, tmp_path: Path, monkeypatch):
        p = Sentinel2Preprocessor()
        safe = _safe_tree(tmp_path)
        p._resolve_paths(safe)

        xml_text = """
        <root>
          <Sun_Angles_Grid>
            <Zenith>
              <Values_List>
                <VALUES>30 31</VALUES>
                <VALUES>32 33</VALUES>
              </Values_List>
            </Zenith>
            <Azimuth>
              <Values_List>
                <VALUES>120 121</VALUES>
                <VALUES>122 123</VALUES>
              </Values_List>
            </Azimuth>
          </Sun_Angles_Grid>
          <Mean_Viewing_Incidence_Angle>
            <ZENITH_ANGLE>5.0</ZENITH_ANGLE>
            <AZIMUTH_ANGLE>100.0</AZIMUTH_ANGLE>
          </Mean_Viewing_Incidence_Angle>
          <Mean_Viewing_Incidence_Angle>
            <ZENITH_ANGLE>7.0</ZENITH_ANGLE>
            <AZIMUTH_ANGLE>110.0</AZIMUTH_ANGLE>
          </Mean_Viewing_Incidence_Angle>
        </root>
        """
        root = ET.fromstring(xml_text)
        sun = p._parse_sun_angles(root)
        view = p._parse_view_angles(root)
        assert sun["zenith"].shape == (2, 2)
        assert view["zenith"].shape == (23, 23)

        # Fallbacks
        fallback = ET.fromstring("<root><Zenith/></root>")
        vals = p._parse_angle_grid(fallback.find("Zenith"), "")
        assert vals.shape == (23, 23)
        view_fallback = p._parse_view_angles(ET.fromstring("<root/>"))
        assert view_fallback["azimuth"].shape == (23, 23)

        def _fake_reproject(source, target, resampling="bilinear"):
            return xr.full_like(target, float(np.nanmean(source.values)))

        monkeypatch.setattr("siac.satellite.sentinel2.reproject_match", _fake_reproject)
        out = p._angles_to_grid(np.full((23, 23), 30.0, dtype=np.float32), _da((8, 8), 1.0))
        assert out.shape == (8, 8)

    def test_load_toa_extract_geometry_cloud_and_metadata(self, tmp_path: Path, monkeypatch):
        safe = _safe_tree(tmp_path)
        p = Sentinel2Preprocessor()

        def _fake_read_raster(path):
            name = Path(path).name
            if "B10" in name:
                return _da((6, 6), value=500.0)
            if "B02" in name:
                return _da((6, 6), value=4500.0)
            return _da((6, 6), value=3500.0)

        monkeypatch.setattr("siac.satellite.sentinel2.read_raster", _fake_read_raster)
        monkeypatch.setattr("siac.satellite.sentinel2.reproject_match", lambda src, tgt, resampling="bilinear": src)

        # product + granule metadata files
        _write(
            safe / "MTD_MSIL1C.xml",
            "<root><QUANTIFICATION_VALUE>10000</QUANTIFICATION_VALUE><RADIO_ADD_OFFSET band_id=\"4\">-100</RADIO_ADD_OFFSET></root>",
        )
        _write(
            safe / "GRANULE" / "L1C_TILE" / "MTD_TL.xml",
            "<root><SENSING_TIME>2024-01-01T10:00:00.000Z</SENSING_TIME><TILE_ID>T31UDQ</TILE_ID></root>",
        )

        toa = p.load_toa(safe)
        assert {"B02", "B04", "B10"}.issubset(toa.data_vars)
        assert float(toa["B04"].max()) <= 1.5

        cloud = p.extract_cloud_mask(safe)
        assert cloud.dtype == bool
        assert cloud.name == "cloud_mask"

        md = p.get_metadata(safe)
        assert md["sensor"] == "MSI"
        assert md["tile_id"] == "T31UDQ"

        # Extract geometry path with XML parser and resampling
        xml_text = """
        <root>
          <Sun_Angles_Grid>
            <Zenith><Values_List><VALUES>30 30</VALUES><VALUES>30 30</VALUES></Values_List></Zenith>
            <Azimuth><Values_List><VALUES>120 120</VALUES><VALUES>120 120</VALUES></Values_List></Azimuth>
          </Sun_Angles_Grid>
          <Mean_Viewing_Incidence_Angle>
            <ZENITH_ANGLE>6.0</ZENITH_ANGLE>
            <AZIMUTH_ANGLE>101.0</AZIMUTH_ANGLE>
          </Mean_Viewing_Incidence_Angle>
        </root>
        """
        _write(safe / "GRANULE" / "L1C_TILE" / "MTD_TL.xml", xml_text)
        geom = p.extract_geometry(safe)
        assert geom.sza.size > 0
        assert geom.vza.size > 0
