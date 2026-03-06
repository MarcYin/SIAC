"""Tests for official Sentinel-2 SRF source discovery and parsing."""

from __future__ import annotations

from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pytest

from siac.core.types import SENTINEL2C_CONFIG
from siac.srf.builders import (
    BandCharacterization,
    build_sensor_config_from_band_characterization,
)
from siac.srf.catalog import list_known_srf_sources
from siac.srf.loaders import (
    load_sensor_config_from_local_srf_file,
    load_sensor_config_from_remote_srf,
    load_sensor_config_from_srf,
)
from siac.srf.sources.sentinel2 import (
    build_sentinel2_sensor_config,
    load_sentinel2_sensor_config,
    parse_sentinel2_srf_workbook,
    parse_sentiwiki_s2_srf_assets,
)
from siac.srf.types import SpectralResponseFunction


def _xlsx_sheet_xml(satellite_id: str) -> str:
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <sheetData>
    <row r="1">
      <c r="A1" t="inlineStr"><is><t>SR_WL</t></is></c>
      <c r="B1" t="inlineStr"><is><t>{satellite_id}_SR_AV_B1</t></is></c>
      <c r="C1" t="inlineStr"><is><t>{satellite_id}_SR_AV_B2</t></is></c>
      <c r="D1" t="inlineStr"><is><t>{satellite_id}_SR_AV_B8A</t></is></c>
    </row>
    <row r="2">
      <c r="A2"><v>440</v></c>
      <c r="B2"><v>0.0</v></c>
      <c r="C2"><v>0.0</v></c>
      <c r="D2"><v>0.0</v></c>
    </row>
    <row r="3">
      <c r="A3"><v>450</v></c>
      <c r="B3"><v>0.5</v></c>
      <c r="C3"><v>0.5</v></c>
      <c r="D3"><v>0.5</v></c>
    </row>
    <row r="4">
      <c r="A4"><v>460</v></c>
      <c r="B4"><v>1.0</v></c>
      <c r="C4"><v>1.0</v></c>
      <c r="D4"><v>1.0</v></c>
    </row>
    <row r="5">
      <c r="A5"><v>470</v></c>
      <c r="B5"><v>0.5</v></c>
      <c r="C5"><v>0.5</v></c>
      <c r="D5"><v>0.5</v></c>
    </row>
    <row r="6">
      <c r="A6"><v>480</v></c>
      <c r="B6"><v>0.0</v></c>
      <c r="C6"><v>0.0</v></c>
      <c r="D6"><v>0.0</v></c>
    </row>
  </sheetData>
</worksheet>
"""


def _minimal_s2_workbook(path: Path) -> Path:
    workbook_xml = """<?xml version="1.0" encoding="UTF-8"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
          xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <sheets>
    <sheet name="Spectral Responses (S2A)" sheetId="1" r:id="rId1"/>
    <sheet name="Spectral Responses (S2B)" sheetId="2" r:id="rId2"/>
    <sheet name="Spectral Responses (S2C)" sheetId="3" r:id="rId3"/>
  </sheets>
</workbook>
"""
    rels_xml = """<?xml version="1.0" encoding="UTF-8"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet2.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet3.xml"/>
</Relationships>
"""
    with ZipFile(path, "w") as zf:
        zf.writestr("xl/workbook.xml", workbook_xml)
        zf.writestr("xl/_rels/workbook.xml.rels", rels_xml)
        zf.writestr("xl/worksheets/sheet1.xml", _xlsx_sheet_xml("S2A"))
        zf.writestr("xl/worksheets/sheet2.xml", _xlsx_sheet_xml("S2B"))
        zf.writestr("xl/worksheets/sheet3.xml", _xlsx_sheet_xml("S2C"))
    return path


def _all_band_srfs(satellite_id: str) -> dict[str, SpectralResponseFunction]:
    band_names = [
        "B01", "B02", "B03", "B04", "B05", "B06", "B07",
        "B08", "B8A", "B09", "B10", "B11", "B12",
    ]
    result: dict[str, SpectralResponseFunction] = {}
    for idx, band_name in enumerate(band_names):
        wl = np.array([440.0, 450.0, 460.0, 470.0, 480.0], dtype=np.float32) + idx
        resp = np.array([0.0, 0.5, 1.0, 0.5, 0.0], dtype=np.float32)
        result[band_name] = SpectralResponseFunction.from_tabulated(
            sensor_id="MSI",
            satellite_id=satellite_id,
            band_name=band_name,
            wavelengths_nm=wl,
            response=resp,
        )
    return result


def test_parse_sentiwiki_s2_srf_assets_picks_latest_attachment():
    html = """
    <a class="filename" href="../__attachments/1/old.xlsx?inst-v=x"
       data-filename="COPE-GSEG-EOPG-TN-15-0007 - Sentinel-2 Spectral Response Functions 2022 - 3.2.xlsx"></a>
    <a class="filename" href="../__attachments/2/new.xlsx?inst-v=y"
       data-filename="COPE-GSEG-EOPG-TN-15-0007 - Sentinel-2 Spectral Response Functions 2024 - 4.0.xlsx"></a>
    """

    assets = parse_sentiwiki_s2_srf_assets(html, base_url="https://sentiwiki.copernicus.eu/web/s2-documents")

    assert len(assets) == 2
    assert assets[-1].filename.endswith("2024 - 4.0.xlsx")
    assert assets[-1].url == "https://sentiwiki.copernicus.eu/__attachments/2/new.xlsx?inst-v=y"


def test_parse_sentiwiki_s2_srf_assets_uses_anchor_text_when_data_filename_missing():
    html = """
    <a class="filename" href="../__attachments/2/new.xlsx?inst-v=y">
      COPE-GSEG-EOPG-TN-15-0007 - Sentinel-2 Spectral Response Functions 2024 - 4.0.xlsx
    </a>
    """

    assets = parse_sentiwiki_s2_srf_assets(html, base_url="https://sentiwiki.copernicus.eu/web/s2-documents")

    assert len(assets) == 1
    assert assets[0].filename.endswith("2024 - 4.0.xlsx")
    assert assets[0].version == (4, 0)


def test_parse_sentinel2_srf_workbook_extracts_platform_bands(tmp_path: Path):
    workbook = _minimal_s2_workbook(tmp_path / "s2_srf.xlsx")

    srfs = parse_sentinel2_srf_workbook(
        workbook,
        source_url="https://example.com/s2.xlsx",
        source_version="4.0",
    )

    assert set(srfs) == {"S2A", "S2B", "S2C"}
    assert set(srfs["S2C"]) == {"B01", "B02", "B8A"}
    srf = srfs["S2C"]["B02"]
    assert srf.sensor_id == "MSI"
    assert srf.satellite_id == "S2C"
    assert srf.band_name == "B02"
    assert srf.source_version == "4.0"
    assert np.trapezoid(srf.response, srf.wavelengths_nm) == pytest.approx(1.0)


def test_build_sentinel2_sensor_config_attaches_real_srfs(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        "siac.srf.sources.sentinel2.load_sentinel2_srfs",
        lambda **_kwargs: {
            "S2A": _all_band_srfs("S2A"),
            "S2B": _all_band_srfs("S2B"),
            "S2C": _all_band_srfs("S2C"),
        },
    )

    config = build_sentinel2_sensor_config("S2C")

    assert config.satellite_id == "S2C"
    assert len(config.bands) == 13
    band = config.get_band("B02")
    assert band.has_srf
    assert band.srf_wavelengths_nm is not None
    assert band.srf_response is not None
    assert band.center_wavelength == pytest.approx(461.0)


def test_load_sentinel2_sensor_config_has_no_builtin_fallback(monkeypatch: pytest.MonkeyPatch):
    def _boom(*_args, **_kwargs):
        raise RuntimeError("download failed")

    monkeypatch.setattr("siac.srf.sources.sentinel2.build_sentinel2_sensor_config", _boom)

    with pytest.raises(RuntimeError, match="download failed"):
        load_sentinel2_sensor_config("S2C")


def test_load_sensor_config_from_remote_srf_delegates_to_sentinel2_loader(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        "siac.srf.loaders.load_sentinel2_sensor_config",
        lambda satellite_id, **kwargs: (satellite_id, kwargs),
    )

    out = load_sensor_config_from_remote_srf(
        "MSI",
        "S2C",
        cache_dir="/tmp/srf-cache",
        refresh=True,
    )

    assert out == ("S2C", {"cache_dir": "/tmp/srf-cache", "refresh": True})


def test_load_sensor_config_from_local_srf_file_delegates_to_sentinel2_loader(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        "siac.srf.loaders.load_sentinel2_sensor_config",
        lambda satellite_id, **kwargs: (satellite_id, kwargs),
    )

    out = load_sensor_config_from_local_srf_file(
        "MSI",
        "S2C",
        "/tmp/local_s2.xlsx",
        cache_dir="/tmp/srf-cache",
    )

    assert out == (
        "S2C",
        {
            "cache_dir": "/tmp/srf-cache",
            "refresh": False,
            "workbook_path": Path("/tmp/local_s2.xlsx"),
        },
    )


def test_load_sensor_config_from_srf_routes_to_local_loader(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        "siac.srf.loaders.load_sensor_config_from_local_srf_file",
        lambda *args, **kwargs: (args, kwargs),
    )

    out = load_sensor_config_from_srf(
        "MSI",
        "S2C",
        local_path="/tmp/local_s2.xlsx",
    )

    assert out[0] == ("MSI", "S2C", "/tmp/local_s2.xlsx")


def test_build_sensor_config_from_band_characterization_attaches_gaussian_srfs():
    characterization = {
        band.name: BandCharacterization(
            band_name=band.name,
            center_wavelength_nm=band.center_wavelength + 1.0,
            fwhm_nm=band.bandwidth + 2.0,
        )
        for band in SENTINEL2C_CONFIG.bands
    }

    config = build_sensor_config_from_band_characterization(SENTINEL2C_CONFIG, characterization)

    band = config.get_band("B02")
    assert band.has_srf
    assert band.center_wavelength == pytest.approx(SENTINEL2C_CONFIG.get_band("B02").center_wavelength + 1.0, abs=0.5)
    assert band.bandwidth == pytest.approx(SENTINEL2C_CONFIG.get_band("B02").bandwidth + 2.0, abs=1.0)


def test_list_known_srf_sources_includes_requested_sensor_families():
    known = {(item.sensor_id, item.satellite_id) for item in list_known_srf_sources()}

    expected = {
        ("MSI", "S2A"),
        ("MSI", "S2B"),
        ("MSI", "S2C"),
        ("OLI", "L8"),
        ("OLI2", "L9"),
        ("MODIS", "Terra"),
        ("MODIS", "Aqua"),
        ("VIIRS", "SNPP"),
        ("VIIRS", "NOAA-20"),
        ("OLCI", "S3A"),
        ("OLCI", "S3B"),
        ("PRISMA_HSI", "PRISMA"),
        ("ENMAP_HSI", "EnMAP"),
        ("EMIT", "ISS"),
        ("DOVE", "PS2"),
        ("SUPERDOVE", "PS2.SD"),
        ("SUPERDOVE", "PSB.SD"),
    }

    assert expected.issubset(known)
