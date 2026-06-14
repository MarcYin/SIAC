"""Regression: a geographic AOI must not override a projected observation grid.

When the scene TOA is in a projected CRS (e.g. UTM) and the caller supplies a
geographic (EPSG:4326) AOI, clipping keeps the TOA in UTM. The resulting
ObservationBundle must report the TOA's UTM crs/bounds (metres), not the AOI's
degrees — otherwise the metre-based aerosol solver grid collapses to 1x1
(`ceil(0.1 degrees / 120 metres) == 1`).
"""

from __future__ import annotations

import numpy as np
import rioxarray  # noqa: F401
import xarray as xr

from siac.adapters.satellite.observation import raw_output_to_observation_bundle
from siac.domain.aoi import AOI
from siac.domain.sensors import SensorBand, SensorConfig
from siac.runtime import GeometryAngles


def _utm_dataarray(name: str, value: float = 0.1) -> xr.DataArray:
    # 200x200 grid at 10 m in UTM 31N over Palaiseau.
    res = 10.0
    x = 445000.0 + np.arange(200) * res
    y = 5402000.0 - np.arange(200) * res
    da = xr.DataArray(
        np.full((len(y), len(x)), value, dtype=np.float32),
        dims=("y", "x"),
        coords={"x": x, "y": y},
        name=name,
    )
    return da.rio.write_crs("EPSG:32631")


def _raw_utm_scene() -> dict:
    toa = xr.Dataset({b: _utm_dataarray(b) for b in ("B02", "B04", "B08")})
    angle = _utm_dataarray("sza", 0.5)
    return {
        "toa": toa,
        "geometry": GeometryAngles(sza=angle, saa=angle, vza=angle, vaa=angle),
        "cloud_mask": _utm_dataarray("cloud", 0.0).astype(bool),
        "metadata": {"observation_time": "2024-02-17T11:00:00"},
    }


def _sensor_config() -> SensorConfig:
    bands = (
        SensorBand(
            name="B02", center_wavelength=490.0, bandwidth=65.0, resolution=10.0, band_index=0
        ),
        SensorBand(
            name="B04", center_wavelength=665.0, bandwidth=30.0, resolution=10.0, band_index=1
        ),
        SensorBand(
            name="B08", center_wavelength=842.0, bandwidth=115.0, resolution=10.0, band_index=2
        ),
    )
    return SensorConfig(sensor_id="MSI", satellite_id="S2A", bands=bands)


def test_geographic_aoi_keeps_projected_observation_grid() -> None:
    aoi = AOI.from_bounds((2.18, 48.68, 2.25, 48.74), crs="EPSG:4326")
    bundle = raw_output_to_observation_bundle(
        _raw_utm_scene(), sensor_config=_sensor_config(), aoi=aoi
    )
    # Observation must report the TOA's projected grid, not the AOI's degrees.
    assert "32631" in bundle.crs
    xmin, ymin, xmax, ymax = bundle.bounds
    assert (xmax - xmin) > 1000.0  # metres, not ~0.07 degrees
    assert (ymax - ymin) > 1000.0
    # And it must be consistent with the actual TOA grid extent.
    toa_bounds = bundle.toa["B02"].rio.bounds()
    assert np.allclose(bundle.bounds, toa_bounds, rtol=0, atol=20.0)


def test_no_aoi_still_uses_scene_grid() -> None:
    bundle = raw_output_to_observation_bundle(
        _raw_utm_scene(), sensor_config=_sensor_config(), aoi=None
    )
    assert "32631" in bundle.crs
    assert (bundle.bounds[2] - bundle.bounds[0]) > 1000.0
