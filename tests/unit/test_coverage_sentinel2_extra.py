"""Additional coverage tests for Sentinel-2 preprocessor internals."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from siac.adapters.satellite.sentinel2 import Sentinel2Preprocessor
from siac.catalog import SENTINEL2A_CONFIG
from siac.domain import SensorConfig


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
    def test_resolve_paths_variants(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
        def _load_config(sensor_id: str, satellite_id: str, *, rsrf_root=None) -> SensorConfig:
            _ = rsrf_root
            return SensorConfig(
                sensor_id=sensor_id,
                satellite_id=satellite_id,
                bands=SENTINEL2A_CONFIG.bands,
                default_ref_scale=SENTINEL2A_CONFIG.default_ref_scale,
                default_ref_offset=SENTINEL2A_CONFIG.default_ref_offset,
            )

        monkeypatch.setattr(
            "siac.adapters.satellite.sentinel2.load_sensor_config_with_rsrf",
            _load_config,
        )

        p = Sentinel2Preprocessor()

        safe_a = _safe_tree(tmp_path, "S2A_DEMO.SAFE")
        p._resolve_paths(safe_a)
        assert p._granule_path is not None
        assert p.sensor_config.satellite_id == "S2A"

        p2 = Sentinel2Preprocessor()
        safe_b = _safe_tree(tmp_path, "S2B_DEMO.SAFE")
        p2._resolve_paths(safe_b)
        assert p2.sensor_config.satellite_id == "S2B"

        p25 = Sentinel2Preprocessor()
        safe_c = _safe_tree(tmp_path, "S2C_DEMO.SAFE")
        p25._resolve_paths(safe_c)
        assert p25.sensor_config.satellite_id == "S2C"

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

    def test_resolve_paths_refreshes_cached_state_for_new_input(self, tmp_path: Path):
        p = Sentinel2Preprocessor()

        safe_a = _safe_tree(tmp_path, "S2A_FIRST.SAFE")
        safe_b = _safe_tree(tmp_path, "S2B_SECOND.SAFE")

        p._resolve_paths(safe_a)
        first_granule = p._granule_path
        assert first_granule is not None
        assert p.sensor_config.satellite_id == "S2A"

        p._resolve_paths(safe_b)

        assert p._resolved_input_path == safe_b
        assert p._granule_path == safe_b / "GRANULE" / "L1C_TILE"
        assert p._granule_path != first_granule
        assert p.sensor_config.satellite_id == "S2B"

    def test_resolve_paths_rejects_missing_input_path(self, tmp_path: Path):
        p = Sentinel2Preprocessor()

        with pytest.raises(FileNotFoundError, match="does not exist"):
            p._resolve_paths(tmp_path / "missing.SAFE")

    def test_sensor_config_uses_real_srf_loader(self, monkeypatch: pytest.MonkeyPatch):
        p = Sentinel2Preprocessor(config={"rsrf_root": "/tmp/rsrf-root"})
        p._satellite_id = "S2C"

        class _Cfg:
            satellite_id = "S2C"

        seen: dict[str, object] = {}

        def _load(sensor_id: str, satellite_id: str, *, rsrf_root=None):
            seen["sensor_id"] = sensor_id
            seen["satellite_id"] = satellite_id
            seen["rsrf_root"] = rsrf_root
            return _Cfg() if satellite_id == "S2C" else None

        monkeypatch.setattr(
            "siac.adapters.satellite.sentinel2.load_sensor_config_with_rsrf",
            _load,
        )
        assert p.sensor_config.satellite_id == "S2C"
        assert seen == {
            "sensor_id": "MSI",
            "satellite_id": "S2C",
            "rsrf_root": "/tmp/rsrf-root",
        }

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

        monkeypatch.setattr("siac.adapters.satellite.sentinel2.reproject_match", _fake_reproject)
        out = p._angles_to_grid(np.full((23, 23), 30.0, dtype=np.float32), _da((8, 8), 1.0))
        assert out.shape == (8, 8)

    def test_extract_geometry_requires_granule_metadata(self, tmp_path: Path):
        p = Sentinel2Preprocessor()
        safe = _safe_tree(tmp_path)
        p._resolve_paths(safe)

        with pytest.raises(FileNotFoundError, match="No granule metadata XML found"):
            p.extract_geometry(safe)

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

        monkeypatch.setattr("siac.adapters.satellite.sentinel2.read_raster", _fake_read_raster)
        monkeypatch.setattr(
            "siac.adapters.satellite.sentinel2.reproject_match",
            lambda src, _tgt, **_kwargs: src,
        )

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

    def test_cloud_mask_error_paths_and_settings(self, tmp_path: Path, monkeypatch):
        safe = _safe_tree(tmp_path)
        p = Sentinel2Preprocessor(config={"cloud_mask": {"external_mask_path": "masks/cloud.tif"}})

        toa = _toa = xr.Dataset({"B04": _da((2, 2), 0.2)})

        calls = {"n": 0}

        def _fake_build(*args, **kwargs):  # noqa: ANN001
            calls["n"] += 1
            if kwargs.get("mode") == "none":
                return xr.full_like(_da((2, 2), 0.2), 1, dtype=np.uint8)
            raise ValueError("Could not find any red band")

        monkeypatch.setattr("siac.adapters.satellite.sentinel2.build_cloud_classes", _fake_build)
        out = p.extract_cloud_mask(safe, toa=toa)
        assert out.dtype == bool
        assert calls["n"] == 2  # auto + fallback none

        monkeypatch.setattr(
            "siac.adapters.satellite.sentinel2.build_cloud_classes",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(ValueError("other")),
        )
        with pytest.raises(ValueError, match="other"):
            p.extract_cloud_mask(safe, toa=toa)

        settings = p._cloud_mask_settings(safe)
        assert settings["external_mask_path"] == safe / "masks" / "cloud.tif"
        assert p._get_namespace(ET.fromstring("<root/>")) == ""

    def test_cloud_mask_settings_expands_user_paths(self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
        monkeypatch.setenv("HOME", str(tmp_path))
        safe = _safe_tree(tmp_path)
        p = Sentinel2Preprocessor(config={"cloud_mask": {"external_mask_path": "~/clouds/mask.tif"}})

        settings = p._cloud_mask_settings(safe)
        assert settings["external_mask_path"] == tmp_path / "clouds" / "mask.tif"

    def test_finders_and_img_data_fallback_paths(self, tmp_path: Path):
        p = Sentinel2Preprocessor()
        safe = _safe_tree(tmp_path)
        p._resolve_paths(safe)

        # Product XML second-candidate branch and no-candidate branch.
        md = safe / "metadata.xml"
        md.write_text("<root/>")
        if (safe / "MTD_MSIL1C.xml").exists():
            (safe / "MTD_MSIL1C.xml").unlink()
        assert p._find_product_xml(safe) == md
        md.unlink()
        assert p._find_product_xml(safe) is None

        # Granule XML fallback candidate and no-candidate branch.
        g = p._granule_path / "X_MTD_ALT.xml"
        g.write_text("<root/>")
        assert p._find_granule_xml() == g
        g.unlink()
        assert p._find_granule_xml() is None

        # IMG_DATA fallback for AWS-like granule layout.
        p2 = Sentinel2Preprocessor()
        aws = tmp_path / "aws_scene"
        aws.mkdir()
        (aws / "metadata.xml").write_text("<root/>")
        p2._resolve_paths(aws)
        assert p2._get_img_data_path() == aws
