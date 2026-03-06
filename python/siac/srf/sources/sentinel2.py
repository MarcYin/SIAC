"""Sentinel-2 official SRF source discovery and parsing."""

from __future__ import annotations

import logging
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from functools import lru_cache
from html import unescape
from html.parser import HTMLParser
from pathlib import Path
from typing import TYPE_CHECKING
from urllib.parse import urljoin
from urllib.request import urlopen
from zipfile import ZipFile

import numpy as np

from siac.core.types import SENTINEL2A_CONFIG, SENTINEL2B_CONFIG, SENTINEL2C_CONFIG, SensorConfig
from siac.srf.builders import build_sensor_config_from_tabulated_srfs
from siac.srf.types import SpectralResponseFunction

if TYPE_CHECKING:
    from collections.abc import Callable

logger = logging.getLogger(__name__)

S2_MISSION_PAGE_URL = "https://sentiwiki.copernicus.eu/web/s2-mission"
S2_DOCUMENTS_PAGE_URL = (
    "https://sentiwiki.copernicus.eu/web/s2-documents"
    "?inheritRedirect=true#S2Documents-SPECTRALRESPONSEFUNCTIONS"
)
_DEFAULT_SRF_CACHE_DIR = Path.home() / ".cache" / "siac" / "srf"
_S2_SHEET_NAMES = {
    "S2A": "Spectral Responses (S2A)",
    "S2B": "Spectral Responses (S2B)",
    "S2C": "Spectral Responses (S2C)",
}
_S2_BASE_CONFIGS = {
    "S2A": SENTINEL2A_CONFIG,
    "S2B": SENTINEL2B_CONFIG,
    "S2C": SENTINEL2C_CONFIG,
}
_S2_FILENAME_RE = re.compile(
    r"Sentinel-2 Spectral Response Functions (?P<year>\d{4}) - (?P<version>\d+(?:\.\d+)*)\.xlsx",
    re.IGNORECASE,
)
@dataclass(frozen=True)
class RemoteSRFAsset:
    """Remote SRF asset discovered from an official source page."""

    filename: str
    url: str
    year: int
    version: tuple[int, ...]

    @property
    def version_string(self) -> str:
        return ".".join(str(part) for part in self.version)


def _parse_version_string(version: str) -> tuple[int, ...]:
    return tuple(int(part) for part in version.split("."))


class _SentiWikiSRFHTMLParser(HTMLParser):
    """Extract Sentinel-2 SRF workbook links from SentiWiki attachment anchors."""

    def __init__(self, base_url: str) -> None:
        super().__init__()
        self.base_url = base_url
        self.assets: list[RemoteSRFAsset] = []
        self._current_href: str | None = None
        self._current_filename: str | None = None
        self._text_chunks: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "a":
            return
        attr_map = dict(attrs)
        self._current_href = attr_map.get("href")
        self._current_filename = attr_map.get("data-filename")
        self._text_chunks = []

    def handle_data(self, data: str) -> None:
        if self._current_href is None:
            return
        stripped = data.strip()
        if stripped:
            self._text_chunks.append(stripped)

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() != "a" or self._current_href is None:
            return

        filename = self._current_filename or " ".join(self._text_chunks)
        href = self._current_href
        self._current_href = None
        self._current_filename = None
        self._text_chunks = []

        if not filename or not href:
            return
        filename = unescape(filename.strip())
        meta = _S2_FILENAME_RE.search(filename)
        if meta is None:
            return
        self.assets.append(
            RemoteSRFAsset(
                filename=filename,
                url=urljoin(self.base_url, unescape(href)),
                year=int(meta.group("year")),
                version=_parse_version_string(meta.group("version")),
            )
        )


def parse_sentiwiki_s2_srf_assets(html: str, *, base_url: str = S2_DOCUMENTS_PAGE_URL) -> list[RemoteSRFAsset]:
    """Parse official Sentinel-2 SRF workbook links from SentiWiki HTML."""
    parser = _SentiWikiSRFHTMLParser(base_url)
    parser.feed(html)
    assets = parser.assets
    assets.sort(key=lambda asset: (asset.year, asset.version, asset.filename))
    return assets


def _read_url_text(url: str, *, opener: Callable[..., object] = urlopen) -> str:
    with opener(url, timeout=30) as response:
        payload = response.read()
    return payload.decode("utf-8", errors="ignore")


def resolve_latest_sentinel2_srf_asset(
    *,
    opener: Callable[..., object] = urlopen,
) -> RemoteSRFAsset:
    """Resolve the latest official Sentinel-2 SRF workbook from SentiWiki."""
    html = _read_url_text(S2_DOCUMENTS_PAGE_URL, opener=opener)
    assets = parse_sentiwiki_s2_srf_assets(html, base_url=S2_DOCUMENTS_PAGE_URL)
    if not assets:
        raise RuntimeError(
            "Could not find any Sentinel-2 SRF workbook links on the SentiWiki S2 documents page."
        )
    return assets[-1]


def _latest_cached_workbook(cache_dir: Path) -> Path | None:
    candidates = sorted(cache_dir.glob("COPE-GSEG-EOPG-TN-15-0007 - Sentinel-2 Spectral Response Functions *.xlsx"))
    return candidates[-1] if candidates else None


def download_sentinel2_srf_workbook(
    *,
    cache_dir: str | Path | None = None,
    refresh: bool = False,
    opener: Callable[..., object] = urlopen,
) -> tuple[Path, RemoteSRFAsset | None]:
    """Download the latest official Sentinel-2 SRF workbook into a local cache."""
    cache_root = Path(cache_dir) if cache_dir is not None else _DEFAULT_SRF_CACHE_DIR
    cache_root.mkdir(parents=True, exist_ok=True)

    try:
        asset = resolve_latest_sentinel2_srf_asset(opener=opener)
        workbook_path = cache_root / asset.filename
        if refresh or not workbook_path.exists():
            logger.info("Downloading official Sentinel-2 SRF workbook %s", asset.filename)
            with opener(asset.url, timeout=60) as response:
                workbook_path.write_bytes(response.read())
        else:
            logger.info("Using cached Sentinel-2 SRF workbook %s", workbook_path)
        return workbook_path, asset
    except Exception:
        cached = _latest_cached_workbook(cache_root)
        if cached is None:
            raise
        logger.warning(
            "Falling back to cached Sentinel-2 SRF workbook at %s because remote resolution failed.",
            cached,
        )
        return cached, None


def _xlsx_shared_strings(zf: ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in zf.namelist():
        return []
    root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
    ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    values: list[str] = []
    for item in root.findall("x:si", ns):
        values.append(
            "".join(
                text.text or ""
                for text in item.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t")
            )
        )
    return values


def _xlsx_value(cell: ET.Element, ns: dict[str, str], shared_strings: list[str]) -> str | None:
    cell_type = cell.attrib.get("t")
    if cell_type == "s":
        value = cell.find("x:v", ns)
        if value is None or value.text is None:
            return None
        return shared_strings[int(value.text)]
    if cell_type == "inlineStr":
        is_elem = cell.find("x:is", ns)
        if is_elem is None:
            return None
        return "".join(text.text or "" for text in is_elem.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t"))
    value = cell.find("x:v", ns)
    return value.text if value is not None else None


def _xlsx_sheet_map(zf: ZipFile) -> dict[str, str]:
    workbook = ET.fromstring(zf.read("xl/workbook.xml"))
    rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
    rel_ns = {"r": "http://schemas.openxmlformats.org/package/2006/relationships"}
    rel_map = {
        rel.attrib["Id"]: rel.attrib["Target"]
        for rel in rels.findall("r:Relationship", rel_ns)
    }
    ns = {
        "x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
        "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    }
    out: dict[str, str] = {}
    for sheet in workbook.find("x:sheets", ns):
        rel_id = sheet.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
        out[sheet.attrib["name"]] = f"xl/{rel_map[rel_id]}"
    return out


def _parse_sheet_rows(zf: ZipFile, sheet_path: str, shared_strings: list[str]) -> list[dict[str, str]]:
    root = ET.fromstring(zf.read(sheet_path))
    ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    rows: list[dict[str, str]] = []
    for row in root.findall(".//x:sheetData/x:row", ns):
        values: dict[str, str] = {}
        for cell in row.findall("x:c", ns):
            value = _xlsx_value(cell, ns, shared_strings)
            if value is None:
                continue
            ref = cell.attrib.get("r", "")
            col = "".join(ch for ch in ref if ch.isalpha())
            if col:
                values[col] = value
        if values:
            rows.append(values)
    return rows


def _normalize_s2_band_name(raw_name: str) -> str:
    band = raw_name.upper()
    if band == "B8A":
        return "B8A"
    if band.startswith("B") and band[1:].isdigit() and len(band) == 2:
        return f"B0{band[1:]}"
    return band


def parse_sentinel2_srf_workbook(
    workbook_path: str | Path,
    *,
    source_url: str | None = None,
    source_version: str | None = None,
) -> dict[str, dict[str, SpectralResponseFunction]]:
    """Parse official Sentinel-2 SRFs from the official workbook."""
    result: dict[str, dict[str, SpectralResponseFunction]] = {}
    with ZipFile(workbook_path) as zf:
        shared_strings = _xlsx_shared_strings(zf)
        sheet_map = _xlsx_sheet_map(zf)
        for satellite_id, sheet_name in _S2_SHEET_NAMES.items():
            if sheet_name not in sheet_map:
                raise KeyError(f"Missing required SRF sheet {sheet_name!r} in workbook {workbook_path}")

            rows = _parse_sheet_rows(zf, sheet_map[sheet_name], shared_strings)
            if not rows:
                raise ValueError(f"SRF sheet {sheet_name!r} is empty in workbook {workbook_path}")
            header = rows[0]
            data_rows = [
                row for row in rows[1:]
                if "A" in row and row["A"] not in {"", None}
            ]
            wavelengths = np.array([float(row["A"]) for row in data_rows], dtype=np.float32)

            srfs_for_satellite: dict[str, SpectralResponseFunction] = {}
            for col, col_name in header.items():
                if col == "A" or not col_name.startswith(f"{satellite_id}_SR_AV_"):
                    continue
                band_name = _normalize_s2_band_name(col_name.rsplit("_", 1)[-1])
                response = np.array(
                    [float(row.get(col, "0") or 0.0) for row in data_rows],
                    dtype=np.float32,
                )
                srfs_for_satellite[band_name] = SpectralResponseFunction.from_tabulated(
                    sensor_id="MSI",
                    satellite_id=satellite_id,
                    band_name=band_name,
                    wavelengths_nm=wavelengths,
                    response=response,
                    source_id="sentiwiki-s2-srf",
                    source_version=source_version,
                    source_url=source_url or S2_MISSION_PAGE_URL,
                )
            result[satellite_id] = srfs_for_satellite
    return result


@lru_cache(maxsize=8)
def _cached_parse_sentinel2_srf_workbook(
    workbook_path: str,
    source_url: str | None,
    source_version: str | None,
) -> dict[str, dict[str, SpectralResponseFunction]]:
    return parse_sentinel2_srf_workbook(
        workbook_path,
        source_url=source_url,
        source_version=source_version,
    )


def load_sentinel2_srfs(
    *,
    cache_dir: str | Path | None = None,
    refresh: bool = False,
    workbook_path: str | Path | None = None,
) -> dict[str, dict[str, SpectralResponseFunction]]:
    """Load official Sentinel-2 SRFs from the local cache or official source."""
    asset = None
    if workbook_path is None:
        workbook_path, asset = download_sentinel2_srf_workbook(
            cache_dir=cache_dir,
            refresh=refresh,
        )
    else:
        workbook_path = Path(workbook_path)
        logger.info("Loading Sentinel-2 SRFs from local workbook %s", workbook_path)

    source_url = None if asset is None else asset.url
    source_version = None if asset is None else asset.version_string
    logger.info("Loading Sentinel-2 SRFs from %s", workbook_path)
    return _cached_parse_sentinel2_srf_workbook(
        str(workbook_path.resolve()),
        source_url,
        source_version,
    )


def build_sentinel2_sensor_config(
    satellite_id: str,
    *,
    cache_dir: str | Path | None = None,
    refresh: bool = False,
    workbook_path: str | Path | None = None,
) -> SensorConfig:
    """Build a Sentinel-2 SensorConfig backed by official platform SRFs."""
    if satellite_id not in _S2_BASE_CONFIGS:
        raise KeyError(f"Unsupported Sentinel-2 platform: {satellite_id!r}")

    base = _S2_BASE_CONFIGS[satellite_id]
    srfs = load_sentinel2_srfs(
        cache_dir=cache_dir,
        refresh=refresh,
        workbook_path=workbook_path,
    )[satellite_id]
    return build_sensor_config_from_tabulated_srfs(
        base,
        srfs,
    )


def load_sentinel2_sensor_config(
    satellite_id: str,
    *,
    cache_dir: str | Path | None = None,
    refresh: bool = False,
    workbook_path: str | Path | None = None,
) -> SensorConfig:
    """Public helper for loading an SRF-backed Sentinel-2 sensor config."""
    return build_sentinel2_sensor_config(
        satellite_id,
        cache_dir=cache_dir,
        refresh=refresh,
        workbook_path=workbook_path,
    )
