"""Unit tests for Earthaccess-backed provider scaffolds (Phase M0)."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

from siac.core.types import SensorBand
from siac.priors.atmospheric.mcd19_earthaccess import MCD19AODProvider, VNP19AODProvider
from siac.priors.atmospheric.merra2 import MERRA2Provider
from siac.priors.brdf.mcd43_earthaccess import (
    MCD19EarthAccessProvider,
    MCD43EarthAccessProvider,
    VNP43EarthAccessProvider,
)
from siac.priors.earthdata_common import MODLAND_SINUSOIDAL_CRS, modland_tile_coords


class _StubEarthAccessSource:
    def __init__(self, downloaded_paths: list[Path]):
        self.downloaded_paths = downloaded_paths
        self.search_calls: list[dict] = []

    def search_granules(self, **kwargs):
        self.search_calls.append(kwargs)
        return [{"id": f"g{i}"} for i in range(len(self.downloaded_paths))]

    def download_granules(self, granules, dest_dir):
        _ = granules, dest_dir
        return list(self.downloaded_paths)


def _full_tile_bounds(h_index: int, v_index: int, shape: tuple[int, int]) -> tuple[float, float, float, float]:
    x, y = modland_tile_coords(h_index, v_index, shape[0], shape[1])
    xres = float(x[1] - x[0]) if len(x) > 1 else 1000.0
    yres = float(y[0] - y[1]) if len(y) > 1 else 1000.0
    return (
        float(x.min() - xres / 2.0),
        float(y.min() - yres / 2.0),
        float(x.max() + xres / 2.0),
        float(y.max() + yres / 2.0),
    )


def test_merra2_provider_returns_default_prior_without_probe():
    provider = MERRA2Provider(probe_earthdata=False)
    state = provider.get_prior(
        bounds=(0.0, 0.0, 1000.0, 1000.0),
        crs="EPSG:4326",
        obs_time=datetime(2024, 1, 1, 12, 0, 0),
        resolution=500.0,
    )

    assert provider.source_name == "MERRA-2"
    assert state.aot.shape == (2, 2)
    assert float(state.aot.mean()) == pytest.approx(0.15)
    assert state.aot_unc.shape == (2, 2)


def test_mcd19_provider_returns_default_prior_without_probe():
    provider = MCD19AODProvider(probe_earthdata=False)
    state = provider.get_prior(
        bounds=(0.0, 0.0, 1000.0, 1000.0),
        crs="EPSG:4326",
        obs_time=datetime(2024, 1, 1, 12, 0, 0),
        resolution=500.0,
    )

    assert provider.source_name == "MCD19"
    assert state.aot.shape == (2, 2)
    assert float(state.aot.mean()) == pytest.approx(0.12)
    assert state.tco3_unc.shape == (2, 2)


def test_mcd43_provider_returns_default_weights_without_probe():
    provider = MCD43EarthAccessProvider(probe_earthdata=False)
    weights = provider.get_brdf_parameters(
        bounds=(0.0, 0.0, 1000.0, 1000.0),
        crs="EPSG:4326",
        obs_time=datetime(2024, 1, 1, 12, 0, 0),
        target_resolution=500.0,
        bands=[1, 2],
        temporal_window=16,
    )

    assert provider.source_name == "MCD43"
    assert weights.f0.shape == (2, 2, 2)
    assert list(weights.f0.coords["band"].values) == [1, 2]
    assert float(weights.f0.mean()) == pytest.approx(0.20)


def test_vnp43_provider_returns_default_weights_without_probe():
    provider = VNP43EarthAccessProvider(probe_earthdata=False)
    weights = provider.get_brdf_parameters(
        bounds=(0.0, 0.0, 1000.0, 1000.0),
        crs="EPSG:4326",
        obs_time=datetime(2024, 1, 1, 12, 0, 0),
        target_resolution=500.0,
        bands=[3],
        temporal_window=16,
    )

    assert provider.source_name == "VNP43"
    assert weights.f0.shape == (1, 2, 2)
    assert float(weights.f2.mean()) == pytest.approx(0.02)


def test_mcd43_provider_parses_real_kernel_fields(monkeypatch):
    granule = Path("/tmp/MCD43A1.A2024001.h29v07.061.fake.hdf")
    provider = MCD43EarthAccessProvider(
        source=_StubEarthAccessSource([granule]),
        probe_earthdata=True,
    )

    def _fake_read_dataset(_path, dataset_name):
        shape = (4, 4)
        params = {
            "BRDF_Albedo_Parameters_Band3": np.dstack(
                [np.full(shape, 100), np.full(shape, 50), np.full(shape, 20)]
            ).astype(np.int16),
            "BRDF_Albedo_Parameters_Band4": np.dstack(
                [np.full(shape, 200), np.full(shape, 60), np.full(shape, 30)]
            ).astype(np.int16),
        }
        if dataset_name in params:
            return params[dataset_name], {"scale_factor": 0.001, "_FillValue": 32767, "valid_range": [0, 32766]}
        if dataset_name.startswith("BRDF_Albedo_Band_Mandatory_Quality_"):
            return np.zeros(shape, dtype=np.uint8), {"_FillValue": 255, "valid_range": [0, 1]}
        raise KeyError(dataset_name)

    monkeypatch.setattr(MCD43EarthAccessProvider, "_read_dataset", staticmethod(_fake_read_dataset))

    sensor_bands = [
        SensorBand("B02", 490.0, 65.0, 10.0, 1),
        SensorBand("B03", 560.0, 35.0, 10.0, 2),
    ]
    bounds = _full_tile_bounds(29, 7, (4, 4))
    resolution = (bounds[2] - bounds[0]) / 4.0

    weights = provider.get_brdf_parameters(
        bounds=bounds,
        crs=MODLAND_SINUSOIDAL_CRS,
        obs_time=datetime(2024, 1, 1, 12, 0, 0),
        target_resolution=resolution,
        bands=sensor_bands,
        temporal_window=2,
    )

    assert list(weights.f0.coords["band"].values) == ["B02", "B03"]
    assert float(weights.f0.sel(band="B02").mean()) == pytest.approx(0.1)
    assert float(weights.f1.sel(band="B02").mean()) == pytest.approx(0.05)
    assert float(weights.f2.sel(band="B03").mean()) == pytest.approx(0.03)


def test_mcd19_brdf_provider_parses_kernel_fields(monkeypatch):
    granule = Path("/tmp/MCD19A3.A2024001.h29v07.061.fake.hdf")
    provider = MCD19EarthAccessProvider(
        source=_StubEarthAccessSource([granule]),
        probe_earthdata=True,
    )

    def _fake_read_dataset(_path, dataset_name):
        shape = (4, 4)
        fields = {
            "Kiso_Band3": np.full(shape, 150, dtype=np.int16),
            "Kvol_Band3": np.full(shape, 70, dtype=np.int16),
            "Kgeo_Band3": np.full(shape, 40, dtype=np.int16),
            "Status_QA": np.ones(shape, dtype=np.int16),
        }
        if dataset_name not in fields:
            raise KeyError(dataset_name)
        return fields[dataset_name], {"scale_factor": 0.001, "_FillValue": -28672, "valid_range": [0, 32766]}

    monkeypatch.setattr("siac.priors.brdf.mcd43_earthaccess.read_hdf4_dataset", _fake_read_dataset)

    bounds = _full_tile_bounds(29, 7, (4, 4))
    resolution = (bounds[2] - bounds[0]) / 4.0
    weights = provider.get_brdf_parameters(
        bounds=bounds,
        crs=MODLAND_SINUSOIDAL_CRS,
        obs_time=datetime(2024, 1, 1),
        target_resolution=resolution,
        bands=[SensorBand("B02", 490.0, 65.0, 10.0, 1)],
        temporal_window=2,
    )

    assert float(weights.f0.mean()) == pytest.approx(0.15)
    assert float(weights.f1.mean()) == pytest.approx(0.07)
    assert float(weights.f2.mean()) == pytest.approx(0.04)


def test_mcd19_provider_parses_aod_and_tcwv(monkeypatch):
    granule = Path("/tmp/MCD19A2.A2024001.h29v07.061.fake.hdf")
    provider = MCD19AODProvider(
        source=_StubEarthAccessSource([granule]),
        probe_earthdata=True,
    )

    def _fake_read_dataset(_path, dataset_name):
        shape = (2, 4, 4)
        fields = {
            "Optical_Depth_055": np.stack(
                [np.full((4, 4), 100), np.full((4, 4), 300)],
                axis=0,
            ).astype(np.int16),
            "AOD_Uncertainty": np.stack(
                [np.full((4, 4), 200), np.full((4, 4), 400)],
                axis=0,
            ).astype(np.int16),
            "Column_WV": np.stack(
                [np.full((4, 4), 1500), np.full((4, 4), 2500)],
                axis=0,
            ).astype(np.int16),
            "AOD_QA": np.ones(shape, dtype=np.uint16),
        }
        if dataset_name not in fields:
            raise KeyError(dataset_name)
        scale = 0.0001 if dataset_name == "AOD_Uncertainty" else 0.001
        return fields[dataset_name], {"scale_factor": scale, "_FillValue": -28672, "valid_range": [0, 30000]}

    monkeypatch.setattr(MCD19AODProvider, "_read_dataset", staticmethod(_fake_read_dataset))

    bounds = _full_tile_bounds(29, 7, (4, 4))
    resolution = (bounds[2] - bounds[0]) / 4.0
    state = provider.get_prior(
        bounds=bounds,
        crs=MODLAND_SINUSOIDAL_CRS,
        obs_time=datetime(2024, 1, 1, 12, 0, 0),
        resolution=resolution,
    )

    assert float(state.aot.mean()) == pytest.approx(0.2)
    assert float(state.aot_unc.mean()) == pytest.approx(0.03)
    assert float(state.tcwv.mean()) == pytest.approx(2.0)


def test_vnp19_provider_parses_aod_and_defaults_tcwv(monkeypatch):
    granule = Path("/tmp/VNP19A2.A2024001.h29v07.002.fake.h5")
    provider = VNP19AODProvider(
        source=_StubEarthAccessSource([granule]),
        probe_earthdata=True,
    )

    def _fake_read_dataset(_path, dataset_name):
        shape = (2, 4, 4)
        fields = {
            "HDFEOS/GRIDS/grid750m/Data Fields/Optical_Depth_055": np.stack(
                [np.full((4, 4), 100), np.full((4, 4), 300)],
                axis=0,
            ).astype(np.int16),
            "HDFEOS/GRIDS/grid750m/Data Fields/AOD_Uncertainty": np.stack(
                [np.full((4, 4), 100), np.full((4, 4), 300)],
                axis=0,
            ).astype(np.int16),
            "HDFEOS/GRIDS/grid750m/Data Fields/AOD_QA": np.ones(shape, dtype=np.uint16),
        }
        if dataset_name not in fields:
            raise KeyError(dataset_name)
        scale = 0.0001 if dataset_name.endswith("AOD_Uncertainty") else 0.001
        return fields[dataset_name], {"scale_factor": scale, "_FillValue": -28672, "valid_range": [0, 30000]}

    monkeypatch.setattr(VNP19AODProvider, "_read_dataset", staticmethod(_fake_read_dataset))

    bounds = _full_tile_bounds(29, 7, (4, 4))
    resolution = (bounds[2] - bounds[0]) / 4.0
    state = provider.get_prior(
        bounds=bounds,
        crs=MODLAND_SINUSOIDAL_CRS,
        obs_time=datetime(2024, 1, 1, 12, 0, 0),
        resolution=resolution,
    )

    assert float(state.aot.mean()) == pytest.approx(0.2)
    assert float(state.tcwv.mean()) == pytest.approx(1.5)
