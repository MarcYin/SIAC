"""Coverage lifts for grid.assembler and priors.surface.prior_store."""

from __future__ import annotations

from datetime import datetime
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr

from siac.algorithms.grid import assembler as asm
from siac.algorithms.surface import prior_store as ps
from siac.domain import (
    AtmosphericState,
    BRDFKernelWeights,
    GeometryAngles,
    ObservationBundle,
    SensorBand,
    SensorConfig,
    SurfacePrior,
)

if TYPE_CHECKING:
    import pytest


def _make_sensor_config(*bands: SensorBand) -> SensorConfig:
    return SensorConfig(sensor_id="MOCK", satellite_id="TEST", bands=tuple(bands))


def _make_obs(sensor_config: SensorConfig, toa_vars: dict[str, xr.DataArray]) -> ObservationBundle:
    shape = next(iter(toa_vars.values())).shape
    geometry = GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.5, dtype=np.float32), dims=["y", "x"]),
        saa=xr.DataArray(np.full(shape, 2.5, dtype=np.float32), dims=["y", "x"]),
        vza=xr.DataArray(np.full(shape, 0.2, dtype=np.float32), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(shape, 1.5, dtype=np.float32), dims=["y", "x"]),
    )
    return ObservationBundle(
        toa=xr.Dataset(toa_vars),
        geometry=geometry,
        cloud_mask=xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"]),
        sensor_config=sensor_config,
        metadata={"observation_time": datetime(2024, 1, 1, 10, 0, 0)},
        crs="EPSG:32632",
        bounds=(0.0, 0.0, float(shape[1]), float(shape[0])),
    )


def _make_atmo(shape: tuple[int, int]) -> AtmosphericState:
    base = xr.DataArray(np.full(shape, 0.2, dtype=np.float32), dims=["y", "x"])
    return AtmosphericState(
        aot=base,
        tcwv=base + 1.0,
        tco3=base + 0.1,
        aot_unc=base * 0.2,
        tcwv_unc=base * 0.2,
        tco3_unc=base * 0.2,
        elevation=xr.zeros_like(base),
    )


def _make_surface(shape: tuple[int, int]) -> SurfacePrior:
    base = xr.DataArray(np.full(shape, 0.2, dtype=np.float32), dims=["y", "x"])
    kernels = BRDFKernelWeights(
        f0=base,
        f1=base,
        f2=base,
        f0_unc=base,
        f1_unc=base,
        f2_unc=base,
    )
    return SurfacePrior(
        boa=base,
        boa_unc=base * 0.1,
        kernels=kernels,
        mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
    )


def test_assembler_nearest_and_padding_paths(monkeypatch: pytest.MonkeyPatch) -> None:
    da = xr.DataArray(np.arange(4, dtype=np.float32).reshape(2, 2), dims=["y", "x"])

    nearest = asm._resample_da(da, (4, 4), method="nearest")
    assert nearest.shape == (4, 4)

    monkeypatch.setattr(
        asm,
        "zoom",
        lambda arr, _factors, **_kwargs: np.array([[arr[0, 0]]], dtype=np.float64),
    )
    padded = asm._resample_da(da, (3, 3), method="bilinear")
    assert padded.shape == (3, 3)

    mask = xr.DataArray(np.array([[True]], dtype=bool), dims=["y", "x"])
    cloud = asm._resample_cloud_mask(mask, (3, 3))
    assert cloud.shape == (3, 3)
    assert bool(cloud.values[1, 1])


def test_assemble_grids_fallback_band_selection_and_toa_fallback() -> None:
    shape = (4, 4)

    # No aerosol-sensitive bands -> fallback to first two sensor bands.
    config_no_aerosol = _make_sensor_config(
        SensorBand("B10", 700.0, 20.0, 10.0, 0),
        SensorBand("B11", 800.0, 20.0, 10.0, 1),
        SensorBand("B12", 900.0, 20.0, 10.0, 2),
    )
    toa = {
        "B10": xr.DataArray(np.ones(shape, dtype=np.float32), dims=["y", "x"]),
        "B11": xr.DataArray(np.ones(shape, dtype=np.float32) * 2, dims=["y", "x"]),
    }
    obs = _make_obs(config_no_aerosol, toa)
    sib = asm.assemble_grids(obs, _make_atmo(shape), _make_surface(shape), rt_model=object(), aux_resolution_m=10.0)
    assert [b.name for b in sib.bands] == ["B10", "B11"]

    # Selected band names missing in TOA dataset -> fallback to first available variable.
    config_with_aerosol = _make_sensor_config(
        SensorBand("B01", 443.0, 20.0, 10.0, 0),
        SensorBand("B02", 490.0, 20.0, 10.0, 1),
    )
    obs_missing = _make_obs(
        config_with_aerosol,
        {"X": xr.DataArray(np.ones(shape, dtype=np.float32), dims=["y", "x"])},
    )
    sib_missing = asm.assemble_grids(
        obs_missing,
        _make_atmo(shape),
        _make_surface(shape),
        rt_model=object(),
        aux_resolution_m=10.0,
    )
    assert sib_missing.toa.shape[0] == 1


def test_prior_store_wraparound_crop_and_select_edge_cases(tmp_path) -> None:
    data = xr.DataArray(np.array([[0.0], [1.0]], dtype=np.float32), dims=["doy", "x"], coords={"doy": [100, 200]})

    before_first = ps._interpolate_doy(data, np.array([100, 200]), 50)
    after_last = ps._interpolate_doy(data, np.array([100, 200]), 250)
    assert float(before_first.values[0]) > 0.0
    assert 0.0 < float(after_last.values[0]) < 1.0

    # Non-directory store path returns empty.
    not_dir = tmp_path / "not_dir.txt"
    not_dir.write_text("x")
    assert ps._select_tiles(not_dir, (0.0, 0.0, 1.0, 1.0)) == []

    root = tmp_path / "store"
    root.mkdir()
    (root / "README.txt").write_text("ignore")

    t1 = root / "T1"
    t1.mkdir()  # no .zattrs -> include by default

    t2 = root / "T2"
    t2.mkdir()
    (t2 / ".zattrs").write_text('{"bounds": null}')

    t3 = root / "T3"
    t3.mkdir()
    (t3 / ".zattrs").write_text("{bad-json")

    tiles = ps._select_tiles(root, (0.0, 0.0, 1.0, 1.0))
    assert tiles == ["T1", "T2", "T3"]

    da2 = xr.DataArray(np.arange(16, dtype=np.float32).reshape(4, 4), dims=["y", "x"])
    cropped2 = ps._crop_to_bounds(da2, (0.0, 0.0, 4.0, 4.0), (1.0, 1.0, 3.0, 3.0), 1.0)
    assert cropped2.ndim == 2

    da4 = xr.DataArray(np.ones((1, 1, 4, 4), dtype=np.float32), dims=["a", "b", "y", "x"])
    unchanged = ps._crop_to_bounds(da4, (0.0, 0.0, 4.0, 4.0), (1.0, 1.0, 3.0, 3.0), 1.0)
    assert unchanged is da4


def test_prior_store_projection_out_of_range_and_default_wavelengths(monkeypatch: pytest.MonkeyPatch, tmp_path) -> None:
    sensor = _make_sensor_config(SensorBand("B1", 900.0, 20.0, 10.0, 0))

    reflectance = np.ones((1, 2, 2), dtype=np.float32)
    projected = ps._project_to_sensor(reflectance, np.array([400.0, 500.0, 600.0]), sensor)
    assert np.all(projected == 0.0)

    # Build a store dataset without "wavelengths" to hit fallback wavelengths path.
    ds = xr.Dataset(
        {
            "reflectance": xr.DataArray(
                np.full((1, 1, 2, 2), 0.3, dtype=np.float32),
                dims=["doy", "band", "y", "x"],
                coords={"doy": [100]},
            ),
            "uncertainty": xr.DataArray(
                np.full((1, 1, 2, 2), 0.05, dtype=np.float32),
                dims=["doy", "band", "y", "x"],
                coords={"doy": [100]},
            ),
        },
        attrs={"bounds": [0.0, 0.0, 2.0, 2.0], "resolution_m": 1.0},
    )

    monkeypatch.setattr(ps, "_select_tiles", lambda store_path, bounds: ["TILE"])  # noqa: ARG005
    monkeypatch.setattr(xr, "open_zarr", lambda _path: ds)

    geom = GeometryAngles(
        sza=xr.DataArray(np.zeros((2, 2), dtype=np.float32), dims=["y", "x"]),
        saa=xr.DataArray(np.zeros((2, 2), dtype=np.float32), dims=["y", "x"]),
        vza=xr.DataArray(np.zeros((2, 2), dtype=np.float32), dims=["y", "x"]),
        vaa=xr.DataArray(np.zeros((2, 2), dtype=np.float32), dims=["y", "x"]),
    )

    prior = ps.PrebuiltPriorStore(tmp_path).get_surface_prior(
        bounds=(0.0, 0.0, 2.0, 2.0),
        crs="EPSG:4326",
        obs_time=datetime(2024, 4, 9),
        sensor_config=sensor,
        geometry=geom,
        resolution=1.0,
    )
    assert prior.boa.shape == (2, 2)
