"""Focused coverage lifts for STAC, validation, and backend helper modules."""

from __future__ import annotations

import dataclasses
import sys
from datetime import datetime, timedelta, timezone
from types import ModuleType, SimpleNamespace
from typing import TYPE_CHECKING

import numpy as np
import pytest
import xarray as xr

import siac.adapters.s2_backend as s2_backend_mod
import siac.storage.stac as stac_mod
from siac.domain import SensorBand
from siac.errors import ValidationError
from siac.runtime import CorrectionDiagnostics, CorrectionResult, GeometryAngles
from siac.runtime.validation import (
    spatial_shape,
    validate_correction_result,
    validate_observation_bundle,
    validate_solved_atmosphere,
    validate_solver_input_bundle,
)

if TYPE_CHECKING:
    from pathlib import Path


def _make_result(
    mock_observation_bundle,
    mock_solved_atmosphere,
    *,
    boa: xr.Dataset | None = None,
    cloud_mask: xr.DataArray | None = None,
    metadata: object | None = None,
    processing_time_s: float | None = 1.5,
) -> CorrectionResult:
    return CorrectionResult(
        boa=boa if boa is not None else mock_observation_bundle.toa[["B02"]],
        boa_unc=None,
        aot=mock_solved_atmosphere.aot,
        tcwv=mock_solved_atmosphere.tcwv,
        cloud_mask=cloud_mask
        if cloud_mask is not None
        else xr.zeros_like(mock_observation_bundle.cloud_mask, dtype=bool),
        metadata={} if metadata is None else metadata,
        diagnostics=CorrectionDiagnostics(processing_time_s=processing_time_s),
    )


class _RaisingRio:
    def bounds(self) -> tuple[float, float, float, float]:
        raise RuntimeError("no bounds")

    @property
    def crs(self) -> object:
        raise RuntimeError("no crs")

    def transform(self) -> tuple[float, ...]:
        raise RuntimeError("no transform")

    def resolution(self) -> tuple[float, float]:
        raise RuntimeError("no resolution")


class _TransformlessRio:
    def __init__(self, *, epsg: int | None = None) -> None:
        self._epsg = epsg

    @property
    def crs(self) -> object:
        if self._epsg is None:
            raise RuntimeError("no epsg")
        return SimpleNamespace(to_epsg=lambda: self._epsg)

    def transform(self) -> tuple[float, ...]:
        raise RuntimeError("no transform")


class _NaNResolutionRio:
    def resolution(self) -> tuple[float, float]:
        return (float("nan"), 10.0)


class _FiniteRio:
    def __init__(self, *, epsg: int = 32632) -> None:
        self.crs = SimpleNamespace(to_epsg=lambda: epsg)

    def transform(self) -> tuple[float, float, float, float, float, float]:
        return (10.0, 0.0, 500000.0, 0.0, -10.0, 4500000.0)

    def resolution(self) -> tuple[float, float]:
        return (10.0, -10.0)


def _install_module(
    monkeypatch: pytest.MonkeyPatch,
    module_name: str,
    **attrs: object,
) -> ModuleType:
    module = ModuleType(module_name)
    for key, value in attrs.items():
        setattr(module, key, value)
    monkeypatch.setitem(sys.modules, module_name, module)
    return module


def test_stac_scalar_helpers_cover_identifier_and_fallback_paths(tmp_path: Path) -> None:
    naive = datetime(2026, 3, 1, 12, 0, 0)
    aware = datetime(2026, 3, 1, 13, 30, tzinfo=timezone(timedelta(hours=1)))
    assert stac_mod._isoformat_utc(naive) == "2026-03-01T12:00:00Z"
    assert stac_mod._isoformat_utc(aware) == "2026-03-01T12:30:00Z"

    assert stac_mod._relative_href(tmp_path / "boa" / "B02.tif", tmp_path) == "boa/B02.tif"

    assert (
        stac_mod._parse_satellite_id(
            "S2B_MSIL1C_20260101T000000.SAFE",
            metadata={},
            fallback="TEST",
        )
        == "S2B"
    )
    assert stac_mod._parse_satellite_id(None, metadata={"satellite": "L8"}, fallback="TEST") == "L8"
    assert stac_mod._parse_satellite_id(None, metadata={}, fallback="TEST") == "TEST"
    assert (
        stac_mod._parse_satellite_id(
            "unexpected_name.zip",
            metadata={"satellite": "S2D"},
            fallback="TEST",
        )
        == "S2D"
    )

    assert stac_mod._platform_name("S2C") == "sentinel-2c"
    assert stac_mod._platform_name("L8") == "landsat-8"
    assert stac_mod._platform_name("TEST") is None
    assert stac_mod._constellation_name("S2A") == "sentinel-2"
    assert stac_mod._constellation_name("L9") == "landsat"
    assert stac_mod._constellation_name("TEST") is None

    assert stac_mod._safe_float("2.5") == pytest.approx(2.5)
    assert stac_mod._safe_float("not-a-number") is None
    assert stac_mod._safe_float(float("inf")) is None


def test_stac_numeric_helpers_handle_empty_nan_and_clipped_values() -> None:
    finite = xr.DataArray(np.array([[np.pi / 6, np.pi / 3]], dtype=np.float64), dims=["y", "x"])
    assert stac_mod._mean_deg(finite) == pytest.approx(45.0)

    empty = xr.DataArray(np.empty((0, 0), dtype=np.float64), dims=["y", "x"])
    assert stac_mod._mean_deg(empty) is None

    with pytest.warns(RuntimeWarning):
        all_nan = xr.DataArray(np.array([[np.nan]], dtype=np.float64), dims=["y", "x"])
        assert stac_mod._mean_deg(all_nan) is None

    assert (
        stac_mod._cloud_cover_percent(
            xr.DataArray(np.array([[2.0, 2.0]], dtype=np.float32), dims=["y", "x"])
        )
        == 100.0
    )
    assert (
        stac_mod._cloud_cover_percent(
            xr.DataArray(np.array([[np.nan]], dtype=np.float32), dims=["y", "x"])
        )
        is None
    )


def test_stac_asset_and_projection_helpers_cover_failure_paths(tmp_path: Path) -> None:
    first_band = SimpleNamespace(sizes={"y": 3, "x": 4}, rio=_TransformlessRio(epsg=None))
    native_bounds = stac_mod._native_bounds(SimpleNamespace(rio=_RaisingRio()), (1, 2, 3, 4))
    assert native_bounds == (1.0, 2.0, 3.0, 4.0)

    bbox, geometry, antimeridian = stac_mod._wgs84_bounds_and_geometry(
        native_bounds, crs=None
    )
    assert bbox == [1.0, 2.0, 3.0, 4.0]
    assert geometry["coordinates"][0][0] == [1.0, 2.0]
    assert antimeridian is False

    proj = stac_mod._proj_properties(first_band, native_bounds)
    assert proj["proj:bbox"] == [1.0, 2.0, 3.0, 4.0]
    assert proj["proj:shape"] == [3, 4]
    assert "proj:epsg" not in proj
    assert "proj:transform" not in proj

    assert stac_mod._gsd(SimpleNamespace(rio=_RaisingRio())) is None
    assert stac_mod._gsd(SimpleNamespace(rio=_NaNResolutionRio())) is None

    common = stac_mod._band_metadata(SensorBand("B02", 490.0, 65.0, 10.0, 1))
    uncommon = stac_mod._band_metadata(SensorBand("CUSTOM", 500.0, 20.0, 30.0, 9))
    assert common["common_name"] == "blue"
    assert "common_name" not in uncommon

    existing = tmp_path / "artifact.tif"
    existing.write_bytes(b"abc")
    assert stac_mod._file_size(existing) == 3
    assert stac_mod._file_size(tmp_path / "missing.tif") is None

    asset = stac_mod._asset_dict(
        "boa/B02.tif",
        title="BOA reflectance B02",
        media_type="image/tiff",
        roles=["data"],
        extra_fields={"custom": "value"},
        file_size=3,
    )
    assert asset["file:size"] == 3
    assert asset["custom"] == "value"


def test_stac_helper_positive_paths_cover_transform_epsg_and_resolution(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _transform_bounds(_src: object, _dst: object, *bounds: float, densify_pts: int = 21):
        _ = densify_pts
        return tuple(v + 0.5 for v in bounds)

    monkeypatch.setattr(stac_mod, "transform_bounds", _transform_bounds)

    bbox, geometry, antimeridian = stac_mod._wgs84_bounds_and_geometry(
        (1.0, 2.0, 3.0, 4.0), "EPSG:32632"
    )
    assert bbox == [1.5, 2.5, 3.5, 4.5]
    assert geometry["coordinates"][0][2] == [3.5, 4.5]
    assert antimeridian is False

    first_band = SimpleNamespace(sizes={"y": 6, "x": 7}, rio=_FiniteRio())
    proj = stac_mod._proj_properties(first_band, (1.0, 2.0, 3.0, 4.0))
    assert proj["proj:epsg"] == 32632
    assert proj["proj:transform"] == [10.0, 0.0, 500000.0, 0.0, -10.0, 4500000.0]
    assert stac_mod._gsd(first_band) == 10.0

    asset = stac_mod._asset_dict(
        "boa/B03.tif",
        title="BOA reflectance B03",
        media_type="image/tiff",
        roles=["data"],
    )
    assert "file:size" not in asset


def test_build_stac_item_uses_fallbacks_and_omits_optional_fields(
    tmp_path: Path,
    mock_observation_bundle,
    mock_solved_atmosphere,
) -> None:
    obs = dataclasses.replace(
        mock_observation_bundle,
        metadata={
            "observation_time": datetime(2026, 1, 2, 3, 4, 5),
            "tile_id": 7,
            "sensor": "MOCK",
        },
    )
    result = _make_result(
        obs,
        mock_solved_atmosphere,
        cloud_mask=xr.ones_like(obs.cloud_mask, dtype=bool),
        metadata={
            **obs.metadata,
            "satellite": "TEST",
            "sensor_config": obs.sensor_config,
            "geometry": obs.geometry,
            "crs": obs.crs,
            "bounds": obs.bounds,
        },
        processing_time_s=None,
    )

    item = stac_mod.build_stac_item_from_result(
        result,
        output_dir=tmp_path,
        artifacts={
            "boa.B02": tmp_path / "missing_b02.tif",
            "preview.false_colour": tmp_path / "preview.tif",
        },
        item_id="custom-item",
    )

    assert item["id"] == "custom-item"
    assert item["properties"]["siac:satellite"] == "TEST"
    assert item["properties"]["siac:tile_id"] == 7
    assert item["properties"]["siac:sensor"] == "MOCK"
    assert "platform" not in item["properties"]
    assert "constellation" not in item["properties"]
    assert "gsd" not in item["properties"]
    assert "siac:processing_time_s" not in item["properties"]
    assert "cloud-mask" not in item["assets"]
    assert "file:size" not in item["assets"]["B02"]
    assert item["assets"]["B02"]["href"] == "missing_b02.tif"
    # Self link points at the on-disk JSON filename; root is "./".
    # (REVIEW.md §3.7 stac.py:540-545)
    assert item["links"] == [
        {
            "rel": "self",
            "href": "custom-item.json",
            "type": "application/geo+json",
        },
        {
            "rel": "root",
            "href": "./",
            "type": "application/json",
        },
    ]


def test_build_stac_item_from_result_includes_optional_assets_and_derived_link(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mock_observation_bundle,
    mock_solved_atmosphere,
) -> None:
    obs = dataclasses.replace(
        mock_observation_bundle,
        metadata={
            "observation_time": datetime(2026, 1, 2, 3, 4, 5, tzinfo=timezone.utc),
            "tile_id": "T50QLD",
            "processing_baseline": "N0511",
            "sensor": "MSI",
            "satellite": "S2C",
        },
    )
    result = _make_result(
        obs,
        mock_solved_atmosphere,
        cloud_mask=xr.zeros_like(obs.cloud_mask, dtype=bool),
        metadata={
            **obs.metadata,
            "input_path": tmp_path
            / "S2C_MSIL1C_20260102T024121_N0511_R089_T50QLD_20260102T035433.SAFE",
            "sensor_config": obs.sensor_config,
            "geometry": obs.geometry,
            "crs": obs.crs,
            "bounds": obs.bounds,
        },
        processing_time_s=12.5,
    )

    boa_path = tmp_path / "boa" / "B02.tif"
    qa_path = tmp_path / "qa" / "cloud_mask.tif"
    boa_path.parent.mkdir()
    qa_path.parent.mkdir()
    boa_path.write_bytes(b"boa")
    qa_path.write_bytes(b"mask")

    monkeypatch.setattr(stac_mod, "_gsd", lambda _first_band: 10.0)

    item = stac_mod.build_stac_item_from_result(
        result,
        output_dir=tmp_path,
        artifacts={
            "boa.B02": boa_path,
            "auxiliary.cloud_mask": qa_path,
        },
    )

    assert item["properties"]["platform"] == "sentinel-2c"
    assert item["properties"]["constellation"] == "sentinel-2"
    assert item["properties"]["gsd"] == 10.0
    assert item["properties"]["eo:cloud_cover"] == 0.0
    assert item["properties"]["siac:processing_time_s"] == 12.5
    assert item["assets"]["cloud-mask"]["href"] == "qa/cloud_mask.tif"
    # Links order: self, root, derived_from (REVIEW.md §3.7 stac.py:540-545).
    rels = [link["rel"] for link in item["links"]]
    assert rels[:2] == ["self", "root"]
    assert "derived_from" in rels


def test_build_stac_item_from_result_rejects_non_datetime_observation_time(
    tmp_path: Path,
    mock_observation_bundle,
    mock_solved_atmosphere,
) -> None:
    """STAC requires a real datetime; string-typed times must raise.

    REVIEW.md §3.7 stac.py:341 — when ``properties['datetime']`` is null,
    STAC mandates ``start_datetime`` + ``end_datetime`` instead. Rather
    than emit an invalid item with ``datetime: null``, callers must hand
    in a real ``datetime`` instance.
    """
    obs = dataclasses.replace(
        mock_observation_bundle,
        metadata={"observation_time": "2026-01-02T03:04:05"},
    )
    result = _make_result(
        obs,
        mock_solved_atmosphere,
        metadata={
            **obs.metadata,
            "sensor_config": obs.sensor_config,
            "geometry": obs.geometry,
            "crs": obs.crs,
            "bounds": obs.bounds,
        },
    )

    with pytest.raises(ValueError, match="STAC item generation requires a datetime"):
        stac_mod.build_stac_item_from_result(
            result,
            output_dir=tmp_path,
            artifacts={"boa.B02": tmp_path / "B02.tif"},
        )


def test_wgs84_bounds_geometry_detects_antimeridian_crossing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """REVIEW.md §3.7 stac.py:145-167 — wrap-around AOI emits 6-element bbox.

    When ``transform_bounds`` returns ``xmax < xmin`` (the canonical
    antimeridian signal from rasterio.warp), we should emit the
    six-element STAC 1.1 bbox ``[xmin, ymin, -180, xmax, ymax, 180]``
    rather than synthesise a globe-spanning rectangle.
    """

    # ``transform_bounds`` for an AOI straddling 180° degrees east in
    # rasterio returns ``xmax`` < ``xmin``. Simulate that.
    def _fake_transform(_src, _dst, *_bounds, densify_pts=21):
        _ = densify_pts
        return (179.5, -1.0, -179.5, 1.0)

    monkeypatch.setattr(stac_mod, "transform_bounds", _fake_transform)
    bbox, _geometry, antimeridian = stac_mod._wgs84_bounds_and_geometry(
        (0.0, 0.0, 1.0, 1.0), "EPSG:32660"
    )
    assert antimeridian is True
    assert len(bbox) == 6
    assert bbox[2] == -180.0
    assert bbox[5] == 180.0


def test_wgs84_bounds_geometry_flags_wide_span() -> None:
    """A >180 deg span is treated as antimeridian-suspect."""
    bbox, _geometry, antimeridian = stac_mod._wgs84_bounds_and_geometry(
        (-170.0, -10.0, 170.0, 10.0), crs=None
    )
    assert antimeridian is True
    assert bbox == [-170.0, -10.0, -180.0, 170.0, 10.0, 180.0]


def test_build_stac_item_emits_root_link_and_filename_self_link(
    tmp_path: Path,
    mock_observation_bundle,
    mock_solved_atmosphere,
) -> None:
    """Self link must reference the on-disk filename; root link is './'.

    REVIEW.md §3.7 stac.py:540-545.
    """
    obs = dataclasses.replace(
        mock_observation_bundle,
        metadata={"observation_time": datetime(2026, 1, 2, tzinfo=timezone.utc)},
    )
    result = _make_result(
        obs,
        mock_solved_atmosphere,
        metadata={
            **obs.metadata,
            "sensor_config": obs.sensor_config,
            "geometry": obs.geometry,
            "crs": obs.crs,
            "bounds": obs.bounds,
            "stac_collection_id": "siac-l2a",
        },
    )

    item = stac_mod.build_stac_item_from_result(
        result,
        output_dir=tmp_path,
        artifacts={"boa.B02": tmp_path / "B02.tif"},
        item_id="MY_ITEM",
    )
    rels = {link["rel"]: link for link in item["links"]}
    assert rels["self"]["href"] == "MY_ITEM.json"
    assert rels["self"]["type"] == "application/geo+json"
    assert rels["root"] == {
        "rel": "root",
        "href": "./",
        "type": "application/json",
    }
    assert rels["collection"]["title"] == "siac-l2a"


def test_build_stac_item_sanitises_invalid_directory_name_to_uuid(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mock_observation_bundle,
    mock_solved_atmosphere,
) -> None:
    """Invalid output_dir.name -> UUID4 fallback (REVIEW.md §3.7 stac.py:316-318)."""
    weird_dir = tmp_path / "has spaces and $pecial!"
    weird_dir.mkdir()

    obs = dataclasses.replace(
        mock_observation_bundle,
        metadata={"observation_time": datetime(2026, 1, 2, tzinfo=timezone.utc)},
    )
    result = _make_result(
        obs,
        mock_solved_atmosphere,
        metadata={
            **obs.metadata,
            "sensor_config": obs.sensor_config,
            "geometry": obs.geometry,
            "crs": obs.crs,
            "bounds": obs.bounds,
        },
    )

    item = stac_mod.build_stac_item_from_result(
        result,
        output_dir=weird_dir,
        artifacts={"boa.B02": weird_dir / "B02.tif"},
    )
    assert item["id"].startswith("siac-")
    # UUID-suffix sanity-check: 8-4-4-4-12 hex segments after the prefix.
    import re as _re

    assert _re.fullmatch(
        r"siac-[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}",
        item["id"],
    )


def test_build_stac_item_requires_explicit_bounds(
    tmp_path: Path,
    mock_observation_bundle,
    mock_solved_atmosphere,
) -> None:
    """Missing bounds + no rio.bounds() -> ValueError, not silent fabrication.

    REVIEW.md §3.7 stac.py:295.
    """
    template = xr.DataArray(
        np.ones((1, 1), dtype=np.float32),
        dims=["y", "x"],
        coords={"y": [0], "x": [0]},
    )
    # `rio` accessor exists on plain DataArrays only after rioxarray import,
    # but ``rio.bounds()`` raises when no transform is set, which is what we
    # want for this test (forces the missing-bounds branch).
    result = SimpleNamespace(
        boa=xr.Dataset({"B02": template}),
        boa_unc=None,
        aot=template,
        tcwv=template,
        cloud_mask=xr.DataArray(np.zeros((1, 1), dtype=bool), dims=["y", "x"]),
        diagnostics=SimpleNamespace(processing_time_s=None),
        metadata={
            "observation_time": datetime(2026, 1, 2, tzinfo=timezone.utc),
            "sensor_config": mock_observation_bundle.sensor_config,
            "geometry": mock_observation_bundle.geometry,
            # No "bounds" and no "crs"
        },
    )

    with pytest.raises(ValueError, match="STAC item generation requires bounds"):
        stac_mod.build_stac_item_from_result(
            result,
            output_dir=tmp_path,
            artifacts={"boa.B02": tmp_path / "B02.tif"},
        )


def test_build_stac_item_drops_eo_extension_when_no_eo_fields(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mock_observation_bundle,
    mock_solved_atmosphere,
) -> None:
    """eo extension URL only present when eo:bands or eo:cloud_cover is set.

    REVIEW.md §3.7 stac.py:551-554. Achieved here by emptying the BOA
    dataset and forcing cloud_cover to be ``None``.
    """
    obs = dataclasses.replace(
        mock_observation_bundle,
        metadata={"observation_time": datetime(2026, 1, 2, tzinfo=timezone.utc)},
    )
    # Empty BOA -> first_band falls through to AOT, eo:bands is empty.
    empty_boa = xr.Dataset()
    result = CorrectionResult(
        boa=empty_boa,
        boa_unc=None,
        aot=mock_solved_atmosphere.aot,
        tcwv=mock_solved_atmosphere.tcwv,
        cloud_mask=xr.zeros_like(mock_observation_bundle.cloud_mask, dtype=bool),
        metadata={
            **obs.metadata,
            "sensor_config": obs.sensor_config,
            "geometry": obs.geometry,
            "crs": obs.crs,
            "bounds": obs.bounds,
        },
        diagnostics=CorrectionDiagnostics(processing_time_s=None),
    )
    monkeypatch.setattr(stac_mod, "_cloud_cover_percent", lambda _mask: None)

    item = stac_mod.build_stac_item_from_result(
        result,
        output_dir=tmp_path,
        artifacts={},
        item_id="test-item",
    )
    assert "eo:bands" not in item["properties"]
    assert "eo:cloud_cover" not in item["properties"]
    assert stac_mod._EO_EXTENSION not in item["stac_extensions"]


def test_build_stac_item_drops_view_fields_for_out_of_range_angles(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mock_observation_bundle,
    mock_solved_atmosphere,
) -> None:
    """SZA > 90° degrees and VZA > 90° must be skipped, not emitted as-is.

    REVIEW.md §3.7 stac.py:373, :377.
    """
    obs = dataclasses.replace(
        mock_observation_bundle,
        metadata={"observation_time": datetime(2026, 1, 2, tzinfo=timezone.utc)},
    )
    result = _make_result(
        obs,
        mock_solved_atmosphere,
        metadata={
            **obs.metadata,
            "sensor_config": obs.sensor_config,
            "geometry": obs.geometry,
            "crs": obs.crs,
            "bounds": obs.bounds,
        },
    )
    # Force out-of-range mean angles: SZA = 200 deg, VZA = 120 deg.
    counter = {"calls": 0}

    def _fake_mean_deg(da):  # noqa: ANN001
        # Order matches build_stac_item_from_result: sza, saa, vza, vaa.
        counter["calls"] += 1
        return [200.0, 12.5, 120.0, 30.0][counter["calls"] - 1]

    monkeypatch.setattr(stac_mod, "_mean_deg", _fake_mean_deg)

    item = stac_mod.build_stac_item_from_result(
        result,
        output_dir=tmp_path,
        artifacts={"boa.B02": tmp_path / "B02.tif"},
        item_id="test-item",
    )
    # SAA / VAA are azimuth fields and pass through; SZA-derived
    # sun_elevation and VZA-derived off_nadir must be dropped.
    assert "view:sun_elevation" not in item["properties"]
    assert "view:off_nadir" not in item["properties"]
    assert "view:sun_azimuth" in item["properties"]
    assert "view:azimuth" in item["properties"]


def test_band_metadata_uses_siac_bandwidth_namespace() -> None:
    """REVIEW.md §3.7 stac.py:208 — drop ``full_width_half_max`` mislabel."""
    band = SensorBand("B02", 490.0, 65.0, 10.0, 1)
    meta = stac_mod._band_metadata(band)
    assert "full_width_half_max" not in meta
    assert meta["siac:band_bandwidth_um"] == pytest.approx(0.065)


def test_spatial_shape_and_observation_validation_error_branches(
    mock_observation_bundle,
) -> None:
    spatial_ds = xr.Dataset({"band": xr.DataArray(np.ones((2, 3)), dims=["y", "x"])})
    fallback_ds = xr.Dataset(
        {"band": xr.DataArray(np.ones((2, 3, 4)), dims=["band", "row", "col"])}
    )
    one_dimensional = xr.Dataset({"band": xr.DataArray(np.ones(3), dims=["time"])})

    assert spatial_shape(spatial_ds) == (2, 3)
    assert spatial_shape(fallback_ds) == (3, 4)
    with pytest.raises(ValueError, match="fewer than 2 dimensions"):
        spatial_shape(one_dimensional)

    no_spatial = dataclasses.replace(
        mock_observation_bundle,
        toa=xr.Dataset({"B02": xr.DataArray(1.0)}),
    )
    with pytest.raises(ValidationError, match="spatial dimensions"):
        validate_observation_bundle(no_spatial)

    bad_bounds = dataclasses.replace(mock_observation_bundle, bounds=(1.0, 2.0, 3.0))
    with pytest.raises(ValidationError, match="bounds must have 4 elements"):
        validate_observation_bundle(bad_bounds)

    bad_crs = dataclasses.replace(mock_observation_bundle, crs="")
    with pytest.raises(ValidationError, match="non-empty string"):
        validate_observation_bundle(bad_crs)

    bad_geometry = dataclasses.replace(
        mock_observation_bundle,
        geometry=GeometryAngles(
            sza=xr.DataArray(np.ones((0, 3)), dims=["y", "x"]),
            saa=mock_observation_bundle.geometry.saa,
            vza=mock_observation_bundle.geometry.vza,
            vaa=mock_observation_bundle.geometry.vaa,
        ),
    )
    with pytest.raises(ValidationError, match=r"geometry\.sza"):
        validate_observation_bundle(bad_geometry)


def test_validation_numeric_guardrails(
    mock_observation_bundle,
    mock_solved_atmosphere,
    mock_solver_input_bundle,
) -> None:
    bad_converged = dataclasses.replace(mock_solved_atmosphere, converged=np.bool_(True))
    with pytest.raises(ValidationError, match="converged must be a boolean"):
        validate_solved_atmosphere(bad_converged)

    bad_iterations_type = dataclasses.replace(mock_solved_atmosphere, n_iterations=1.5)
    with pytest.raises(ValidationError, match="n_iterations must be an int"):
        validate_solved_atmosphere(bad_iterations_type)

    bad_cost = dataclasses.replace(mock_solved_atmosphere, cost_final="high")
    with pytest.raises(ValidationError, match="cost_final must be numeric"):
        validate_solved_atmosphere(bad_cost)

    bad_iterations_value = dataclasses.replace(mock_solved_atmosphere, n_iterations=-1)
    with pytest.raises(ValidationError, match="non-negative"):
        validate_solved_atmosphere(bad_iterations_value)

    wrong_metadata = _make_result(mock_observation_bundle, mock_solved_atmosphere, metadata=[])
    with pytest.raises(ValidationError, match="metadata must be a dictionary"):
        validate_correction_result(wrong_metadata)

    wrong_time_type = _make_result(
        mock_observation_bundle,
        mock_solved_atmosphere,
        processing_time_s="fast",  # type: ignore[arg-type]
    )
    with pytest.raises(ValidationError, match="must be finite when provided"):
        validate_correction_result(wrong_time_type)

    nonfinite_time = _make_result(
        mock_observation_bundle,
        mock_solved_atmosphere,
        processing_time_s=float("inf"),
    )
    with pytest.raises(ValidationError, match="must be finite when provided"):
        validate_correction_result(nonfinite_time)

    bad_aux_resolution = dataclasses.replace(mock_solver_input_bundle, aux_resolution_m=0.0)
    with pytest.raises(ValidationError, match="aux_resolution_m must be positive"):
        validate_solver_input_bundle(bad_aux_resolution)

    bad_aerosol_resolution = dataclasses.replace(
        mock_solver_input_bundle,
        aerosol_resolution_m=-1.0,
    )
    with pytest.raises(ValidationError, match="aerosol_resolution_m must be positive"):
        validate_solver_input_bundle(bad_aerosol_resolution)


def test_build_s2_backend_selects_expected_adapter(monkeypatch: pytest.MonkeyPatch) -> None:
    class _CDSEBackend:
        def __init__(self, *, auth=None) -> None:
            self.auth = auth

    class _GCSBackend:
        pass

    _install_module(
        monkeypatch,
        "siac.adapters.data.copernicus_dataspace",
        CopernicusDataspaceBackend=_CDSEBackend,
    )
    _install_module(
        monkeypatch,
        "siac.adapters.data.gcs_sentinel2",
        GCSSentinel2Backend=_GCSBackend,
    )

    auth = object()
    cdse_config = SimpleNamespace(providers=SimpleNamespace(s2=SimpleNamespace(backend="cdse")))
    gcs_config = SimpleNamespace(providers=SimpleNamespace(s2=SimpleNamespace(backend="gcs")))
    local_config = SimpleNamespace(providers=SimpleNamespace(s2=SimpleNamespace(backend="local")))
    bad_config = SimpleNamespace(providers=SimpleNamespace(s2=SimpleNamespace(backend="mystery")))

    cdse_backend = s2_backend_mod.build_s2_backend(cdse_config, auth=auth)
    gcs_backend = s2_backend_mod.build_s2_backend(gcs_config)

    assert isinstance(cdse_backend, _CDSEBackend)
    assert cdse_backend.auth is auth
    assert isinstance(gcs_backend, _GCSBackend)
    assert s2_backend_mod.build_s2_backend(local_config) is None
    with pytest.raises(ValueError, match="Unknown S2 backend"):
        s2_backend_mod.build_s2_backend(bad_config)
