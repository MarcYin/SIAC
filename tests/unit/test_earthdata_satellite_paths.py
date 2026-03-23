from __future__ import annotations

from datetime import datetime
from types import SimpleNamespace
from typing import TYPE_CHECKING

import numpy as np
import pytest
import xarray as xr

import siac.adapters.earthdata_common as earth_mod
from siac.adapters.satellite.base import (
    BaseSatellitePreprocessor,
    compute_relative_azimuth,
    detect_sensor,
)
from siac.catalog import SENTINEL2A_CONFIG

if TYPE_CHECKING:
    from pathlib import Path


class _DummyPreprocessor(BaseSatellitePreprocessor):
    @property
    def sensor_config(self):
        return SENTINEL2A_CONFIG

    def load_toa(self, input_path: str | Path) -> xr.Dataset:
        del input_path
        return xr.Dataset(
            {"B02": xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"])}
        )

    def extract_geometry(self, input_path: str | Path):
        del input_path
        from siac.runtime import GeometryAngles

        arr = xr.DataArray(np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"])
        return GeometryAngles(sza=arr, saa=arr, vza=arr, vaa=arr)

    def extract_cloud_mask(self, input_path: str | Path, toa=None):  # noqa: ANN001
        del input_path, toa
        return xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"])

    def get_metadata(self, input_path: str | Path) -> dict:
        del input_path
        return {"observation_time": datetime(2024, 1, 1)}


def test_satellite_base_covers_cloud_classes_observation_bundle_and_detection(tmp_path: Path) -> None:
    preprocessor = _DummyPreprocessor()
    preprocessor._last_cloud_classes = xr.DataArray(  # type: ignore[attr-defined]
        np.ones((2, 2), dtype=np.uint8),
        dims=["y", "x"],
    )
    out = preprocessor.preprocess(tmp_path)
    assert "cloud_classes" in out

    bundle = preprocessor.to_observation_bundle(
        tmp_path,
        bounds=(1.0, 2.0, 3.0, 4.0),
        crs="EPSG:4326",
    )
    assert bundle.bounds == (1.0, 2.0, 3.0, 4.0)
    assert bundle.crs == "EPSG:4326"

    saa = xr.DataArray(np.array([[0.2]], dtype=np.float32), dims=["y", "x"])
    vaa = xr.DataArray(np.array([[0.5]], dtype=np.float32), dims=["y", "x"])
    raa = compute_relative_azimuth(saa, vaa)
    assert float(raa.values[0, 0]) == pytest.approx(0.3)

    l8_dir = tmp_path / "l8"
    l8_dir.mkdir()
    (l8_dir / "foo_MTL.txt").write_text("LANDSAT_8", encoding="utf-8")
    assert detect_sensor(l8_dir) == "l8"

    aws_s2 = tmp_path / "aws-s2"
    aws_s2.mkdir()
    (aws_s2 / "metadata.xml").write_text("<xml/>", encoding="utf-8")
    assert detect_sensor(aws_s2) == "s2"


def test_earthdata_common_scalar_and_grid_helpers_cover_error_paths(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    assert earth_mod.decode_attr(np.array([5], dtype=np.int32)) == 5
    assert earth_mod.decode_attr(np.array([b"a", b"b"])) == ["a", "b"]
    assert earth_mod.decode_attr(np.bytes_(b"abc")) == "abc"
    assert earth_mod.attr_scalar({"k": []}, "k", default=7) == 7
    assert earth_mod.attr_scalar({"k": None}, "k", default=9) == 9

    with pytest.raises(ValueError, match="tile indices"):
        earth_mod.parse_tile_indices("bad-name")
    with pytest.raises(ValueError, match="acquisition date"):
        earth_mod.parse_granule_date("bad-name")
    assert earth_mod._decode_gctp_angle(123.0) == 123.0
    assert earth_mod._decode_gctp_angle(160000000.0) == pytest.approx(160.0)
    assert "+lon_0=10.0" in earth_mod._build_sinusoidal_crs(radius=1.0, central_meridian=10.0)

    earth_mod._parse_hdf5_grid_metadata.cache_clear()
    monkeypatch.setattr(earth_mod, "_read_hdf5_struct_metadata", lambda _p: "not matching")
    assert earth_mod._parse_hdf5_grid_metadata("x.h5") == {}

    grid = earth_mod._grid_definition_from_extent(
        width=2,
        height=2,
        upper_left=(0.0, 10.0),
        lower_right=(10.0, 0.0),
        crs="EPSG:4326",
    )
    assert grid.x == (2.5, 7.5)
    assert grid.y == (7.5, 2.5)

    assert earth_mod.granule_geographic_bounds("granule.hdf") is None
    with pytest.raises(ValueError, match="positive"):
        earth_mod.modland_tile_coords(1, 2, 0, 4)
    with pytest.raises(ValueError, match="at least 2 dimensions"):
        earth_mod.make_native_grid_dataarray(np.array([1.0], dtype=np.float32), granule_path="VNP19A2.A2024165.h25v05.h5")
    with pytest.raises(ValueError, match="resolution must be > 0"):
        earth_mod.build_target_template((0.0, 0.0, 1.0, 1.0), "EPSG:4326", 0.0)

    arr = earth_mod.apply_scale_and_mask(
        np.array([[1, -9999, 5]], dtype=np.float32),
        {"_FillValue": -9999, "valid_range": [0, 4], "scale_factor": 2.0, "add_offset": 1.0},
    )
    assert np.isnan(arr[0, 1])
    assert np.isnan(arr[0, 2])
    assert float(arr[0, 0]) == pytest.approx(3.0)

    two_d = earth_mod.reduce_orbit_stack(np.array([[1.0, 2.0]], dtype=np.float32))
    assert two_d.shape == (1, 2)
    reduced = earth_mod.reduce_orbit_stack(
        np.array(
            [
                [[1.0, np.nan]],
                [[3.0, np.nan]],
            ],
            dtype=np.float32,
        )
    )
    assert float(reduced[0, 0]) == pytest.approx(2.0)
    assert np.isnan(reduced[0, 1])

    dataset_path = tmp_path / "dataset.h5"
    with earth_mod.h5py.File(dataset_path, "w") as handle:
        ds = handle.create_dataset("field", data=np.array([[1, 2]], dtype=np.int16))
        ds.attrs["scale_factor"] = np.array([2.0], dtype=np.float64)
    values, attrs = earth_mod.read_hdf5_dataset(dataset_path, "field")
    assert values.shape == (1, 2)
    assert attrs["scale_factor"] == 2.0


def test_earthdata_common_native_grid_and_intersection_fallbacks(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    earth_mod._read_hdf5_root_bounds.cache_clear()
    earth_mod._parse_hdf5_grid_metadata.cache_clear()

    assert earth_mod._read_hdf5_root_bounds("missing.h5") is None

    class _FakeGroup:
        def __iter__(self):
            return iter(["gridA"])

    class _FakeHandle(dict):
        attrs = {}

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):  # noqa: ANN001
            return False

        def get(self, key, default=None):  # noqa: ANN001
            if key == "HDFEOS/GRIDS":
                return _FakeGroup()
            return super().get(key, default)

    monkeypatch.setattr(earth_mod, "_parse_hdf5_grid_metadata", lambda _p: {})
    monkeypatch.setattr(earth_mod.h5py, "File", lambda *_args, **_kwargs: _FakeHandle())
    assert earth_mod._read_hdf5_native_grid_definition("x.h5") is None

    grid = SimpleNamespace(x=(1.0,), y=(2.0,), crs="EPSG:4326")
    monkeypatch.setattr(earth_mod, "_read_hdf5_native_grid_definition", lambda *_args, **_kwargs: grid)
    bounds, crs = earth_mod.granule_native_bounds("x.h5", width=1, height=1)
    assert bounds == (-499.0, -498.0, 501.0, 502.0)
    assert crs == "EPSG:4326"

    monkeypatch.setattr(earth_mod, "_read_hdf5_native_grid_definition", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(earth_mod, "parse_tile_indices", lambda _p: (1, 2))
    monkeypatch.setattr(earth_mod, "modland_tile_bounds", lambda h, v: (float(h), float(v), float(h + 1), float(v + 1)))
    bounds, crs = earth_mod.granule_native_bounds("fallback.h5")
    assert bounds == (1.0, 2.0, 2.0, 3.0)
    assert crs == earth_mod.MODLAND_SINUSOIDAL_CRS

    monkeypatch.setattr(earth_mod, "granule_geographic_bounds", lambda _p: (0.0, 0.0, 2.0, 2.0))
    monkeypatch.setattr(earth_mod, "transform_bounds", lambda bounds, src, dst: bounds)  # noqa: ARG005
    assert earth_mod.granule_intersects_bounds("g.h5", bounds=(1.0, 1.0, 3.0, 3.0), crs="EPSG:4326")
    assert not earth_mod.granule_intersects_bounds("g.h5", bounds=(3.0, 3.0, 4.0, 4.0), crs="EPSG:4326")

    monkeypatch.setattr(earth_mod, "granule_geographic_bounds", lambda _p: None)
    monkeypatch.setattr(earth_mod, "granule_native_bounds", lambda _p: ((0.0, 0.0, 2.0, 2.0), "EPSG:4326"))
    assert earth_mod.granule_intersects_bounds("n.h5", bounds=(1.0, 1.0, 3.0, 3.0), crs="EPSG:4326")
