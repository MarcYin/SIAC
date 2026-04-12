from __future__ import annotations

from typing import TYPE_CHECKING

from siac.adapters.data import water_mask as water_mask_mod

if TYPE_CHECKING:
    from pathlib import Path

    import pytest


def test_required_water_mask_tiles_selects_expected_zenodo_tiles() -> None:
    tiles = water_mask_mod.required_water_mask_tiles(
        bounds=(29.9, -59.9, 30.2, -29.8),
        crs="EPSG:4326",
    )

    assert tiles == (
        "landWater2020_0000_-60.tif",
        "landWater2020_0000_-30.tif",
        "landWater2020_0030_-60.tif",
        "landWater2020_0030_-30.tif",
    )


def test_ensure_local_water_mask_source_downloads_missing_assets_once(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    downloads: list[tuple[str, Path]] = []

    def _fake_download(url: str, destination: Path, *, session=None):  # noqa: ANN001
        del session
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(b"placeholder")
        downloads.append((url, destination))
        return destination

    monkeypatch.setattr(water_mask_mod, "_download_file", _fake_download)

    local_vrt = water_mask_mod.ensure_local_water_mask_source(
        bounds=(30.1, -59.9, 30.2, -59.8),
        crs="EPSG:4326",
        source=water_mask_mod.DEFAULT_WATER_MASK_VRT_URL,
        cache_dir=tmp_path,
    )

    assert local_vrt == tmp_path / "landWater2020.vrt"
    assert [path.name for _url, path in downloads] == [
        "landWater2020.vrt",
        "landWater2020_0030_-60.tif",
    ]

    water_mask_mod.ensure_local_water_mask_source(
        bounds=(30.1, -59.9, 30.2, -59.8),
        crs="EPSG:4326",
        source=water_mask_mod.DEFAULT_WATER_MASK_VRT_URL,
        cache_dir=tmp_path,
    )

    assert len(downloads) == 2


def test_ensure_local_water_mask_source_uses_missing_local_path_parent_as_cache(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    downloads: list[Path] = []

    def _fake_download(url: str, destination: Path, *, session=None):  # noqa: ANN001
        del url, session
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(b"placeholder")
        downloads.append(destination)
        return destination

    monkeypatch.setattr(water_mask_mod, "_download_file", _fake_download)

    local_vrt = water_mask_mod.ensure_local_water_mask_source(
        bounds=(30.1, -59.9, 30.2, -59.8),
        crs="EPSG:4326",
        source=tmp_path / "custom-cache" / "landWater2020.vrt",
    )

    assert local_vrt == tmp_path / "custom-cache" / "landWater2020.vrt"
    assert downloads[0] == local_vrt
