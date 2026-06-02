"""Tests for the thin bestpixel monthly-composite adapter.

The external ``bestpixel`` package is an optional dependency that reaches out
to STAC endpoints at call time, so every test here injects a fake
``bestpixel`` module into ``sys.modules`` (the adapter imports it lazily) and
asserts on the conversion + grid-alignment logic, never the network.
"""

from __future__ import annotations

import sys
import types
from datetime import datetime
from typing import TYPE_CHECKING, Any

import numpy as np
import pytest
import xarray as xr

from siac.adapters.bestpixel import (
    BestPixelMonthlyCompositeProvider,
    _quality_to_cost,
    bestpixel_source_bands,
)
from siac.app._assembly_providers import resolve_monthly_composite_provider
from siac.catalog import SENTINEL2C_CONFIG
from siac.config import SIACConfig
from siac.runtime.models import GeometryAngles, ObservationBundle

if TYPE_CHECKING:
    from collections.abc import Sequence

# UTM 33N grid near the T33KWP scene.
_CRS = "EPSG:32633"
_RES = 60.0
_BOUNDS = (300000.0, 8199400.0, 300600.0, 8200000.0)  # 600 m -> 10x10 px at 60 m


def _geo_layer(value: float, x: np.ndarray, y: np.ndarray) -> xr.DataArray:
    data = xr.DataArray(
        np.full((y.size, x.size), value, dtype=np.float32),
        dims=("y", "x"),
        coords={"y": y, "x": x},
    )
    return data.rio.set_spatial_dims(x_dim="x", y_dim="y").rio.write_crs(_CRS)


def _observation() -> tuple[ObservationBundle, np.ndarray, np.ndarray]:
    width = int((_BOUNDS[2] - _BOUNDS[0]) / _RES)
    height = int((_BOUNDS[3] - _BOUNDS[1]) / _RES)
    x = _BOUNDS[0] + (np.arange(width) + 0.5) * _RES
    y = _BOUNDS[3] - (np.arange(height) + 0.5) * _RES
    geometry = GeometryAngles(
        sza=_geo_layer(30.0, x, y),
        saa=_geo_layer(150.0, x, y),
        vza=_geo_layer(5.0, x, y),
        vaa=_geo_layer(110.0, x, y),
    )
    observation = ObservationBundle(
        toa=xr.Dataset(),
        geometry=geometry,
        cloud_mask=_geo_layer(0.0, x, y),
        sensor_config=SENTINEL2C_CONFIG,
        metadata={"observation_time": datetime(2026, 3, 29, 8, 45)},
        crs=_CRS,
        bounds=_BOUNDS,
    )
    return observation, x, y


def _install_fake_bestpixel(
    monkeypatch: pytest.MonkeyPatch,
    *,
    seen: dict[str, Any] | None = None,
) -> None:
    """Register a fake ``bestpixel`` whose composites overlap the observation.

    The adapter calls ``build_composite`` once per (year, month) with a
    calendar-correct ``datetime`` range, so the fake mirrors that contract and
    records every date range it was asked for under ``seen["datetimes"]``.
    """
    if seen is not None:
        seen.setdefault("datetimes", [])
        seen.setdefault("kwargs", [])

    def _build_composite(
        *,
        bbox: Sequence[float],
        datetime: str,
        resolution: float,
        bands: Sequence[str],
        **kwargs: Any,
    ) -> dict[str, Any]:
        if seen is not None:
            seen["datetimes"].append(datetime)
            seen["kwargs"].append(dict(kwargs))
            seen.setdefault("bbox", list(bbox))
            seen.setdefault("bands", list(bands))
            seen.setdefault("resolution", resolution)
        width, height = 6, 5
        xmin, ymax = 300000.0, 8200000.0
        transform = [resolution, 0.0, xmin, 0.0, -resolution, ymax]
        band_arrays = {
            name: np.full((height, width), 1000 + 200 * i, dtype=np.uint16)
            for i, name in enumerate(bands)
        }
        quality = np.zeros((height, width), dtype=np.uint16)
        quality[0, 0] = 65535  # one nodata pixel
        return {
            "bands": band_arrays,
            "quality": quality,
            "grid": {
                "bounds": [xmin, ymax - height * resolution, xmin + width * resolution, ymax],
                "epsg": 32633,
                "crs": _CRS,
                "resolution": resolution,
                "width": width,
                "height": height,
                "transform": transform,
            },
            "band_names": list(bands),
            "source_ids": ["FAKE_SCENE_1"],
        }

    fake = types.ModuleType("bestpixel")
    fake.build_composite = _build_composite  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "bestpixel", fake)


def test_bestpixel_source_bands_known_and_unknown() -> None:
    bands = bestpixel_source_bands(("coastal", "red", "swir22"))
    assert [b.name for b in bands] == ["coastal", "red", "swir22"]
    assert bands[0].center_wavelength == pytest.approx(443.0)
    assert bands[2].center_wavelength == pytest.approx(2202.0)
    assert [b.band_index for b in bands] == [0, 1, 2]

    # Tagged with the published-library S2 source basis so SIAC's spectral
    # mapper can project bestpixel reflectance onto the scene's sensor.
    assert all(b.rsrf_sensor_unit_id == "sentinel-2a_msi" for b in bands)
    assert [b.rsrf_band_id for b in bands] == ["B01", "B04", "B12"]

    with pytest.raises(ValueError, match="Unknown bestpixel band"):
        bestpixel_source_bands(("not_a_band",))


def test_bestpixel_source_bands_mcd43a4_uses_modis_basis() -> None:
    # mcd43a4 returns MODIS NBAR -> tag as terra_modis with MODIS band ids.
    bands = bestpixel_source_bands(
        ("blue", "green", "red", "nir", "swir16", "swir22"), endpoint="mcd43a4"
    )
    assert all(b.rsrf_sensor_unit_id == "terra_modis" for b in bands)
    assert [b.rsrf_band_id for b in bands] == ["B3", "B4", "B1", "B2", "B6", "B7"]

    # MODIS has no 443 nm coastal band -> requesting it for mcd43a4 must fail
    # loudly rather than silently mis-map.
    with pytest.raises(ValueError, match="no 'terra_modis' source-basis entry"):
        bestpixel_source_bands(("coastal", "red"), endpoint="mcd43a4")


def test_bestpixel_source_bands_hls_uses_s2_basis() -> None:
    # HLS is bandpass-harmonised toward S2, so it uses the S2A source basis.
    bands = bestpixel_source_bands(("blue", "nir"), endpoint="hls")
    assert all(b.rsrf_sensor_unit_id == "sentinel-2a_msi" for b in bands)
    assert [b.rsrf_band_id for b in bands] == ["B02", "B08"]


def test_resolve_years_uses_full_years_before_scene() -> None:
    provider = BestPixelMonthlyCompositeProvider(lookback_years=5)
    years = provider._resolve_years(datetime(2026, 3, 29))
    assert years == (2021, 2022, 2023, 2024, 2025)


def test_resolve_months_defaults_to_scene_month() -> None:
    # Unset -> the scene's own calendar month (seasonal prior).
    provider = BestPixelMonthlyCompositeProvider()
    assert provider._resolve_months(datetime(2026, 3, 29)) == (3,)
    assert provider._resolve_months(datetime(2026, 11, 2)) == (11,)
    # Explicit -> honoured verbatim.
    provider = BestPixelMonthlyCompositeProvider(months=[6, 7, 8])
    assert provider._resolve_months(datetime(2026, 3, 29)) == (6, 7, 8)


def test_quality_to_cost_mapping() -> None:
    quality = np.array([[0, 1, 2, 65535, 7]], dtype=np.uint16)
    cost = _quality_to_cost(quality)
    assert cost[0, 0] == pytest.approx(0.0)
    assert cost[0, 1] == pytest.approx(0.02)
    assert cost[0, 2] == pytest.approx(0.04)
    assert np.isnan(cost[0, 3])  # nodata
    assert cost[0, 4] == pytest.approx(0.04)  # unknown code -> dark tier


def test_get_monthly_composites_builds_collection_on_observation_grid(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen: dict[str, Any] = {}
    _install_fake_bestpixel(monkeypatch, seen=seen)
    observation, x, y = _observation()

    provider = BestPixelMonthlyCompositeProvider(
        endpoint="pc",
        lookback_years=5,
    )
    collection = provider.get_monthly_composites(observation, _RES)

    # Default: 5 lookback years x the scene's own month (March, from the
    # 2026-03-29 observation time) -> a seasonally matched prior.
    assert len(collection.composites) == 5
    assert collection.source_name == "bestpixel:pc monthly composites"
    assert [b.name for b in collection.source_bands] == [b.name for b in provider.source_bands]

    # bestpixel was asked for one March composite per lookback year.
    assert len(seen["datetimes"]) == 5
    assert all(kw["endpoint"] == "pc" for kw in seen["kwargs"])
    assert sorted(seen["datetimes"]) == [
        "2021-03-01/2021-03-31",
        "2022-03-01/2022-03-31",
        "2023-03-01/2023-03-31",
        "2024-03-01/2024-03-31",
        "2025-03-01/2025-03-31",
    ]
    assert {c.month for c in collection.composites} == {3}
    west, south, east, north = seen["bbox"]
    assert -180.0 <= west < east <= 180.0
    assert -90.0 <= south < north <= 90.0

    composite = collection.composites[0]
    assert composite.year == 2021
    assert composite.month == 3
    assert composite.reflectance.dims == ("band", "y", "x")
    assert composite.reflectance.shape == (len(provider.source_bands), y.size, x.size)
    assert composite.quality.shape == (y.size, x.size)
    assert composite.sample_index.shape == (y.size, x.size)

    # Grid alignment with the observation target grid (pixel-perfect).
    np.testing.assert_allclose(composite.reflectance.coords["x"].values, x)
    np.testing.assert_allclose(composite.reflectance.coords["y"].values, y)

    # DN 1000 -> 0.1 reflectance for the first band (coastal).
    coastal = np.asarray(composite.reflectance.sel(band="coastal").values)
    finite = coastal[np.isfinite(coastal)]
    assert finite.size > 0
    np.testing.assert_allclose(finite, 0.1, atol=1e-3)

    # Pixels with finite payload are flagged 0, fully-NaN pixels -1.
    sample = np.asarray(composite.sample_index.values)
    assert set(np.unique(sample)).issubset({-1, 0})
    assert (sample == 0).any()


def test_factory_builds_bestpixel_provider_from_config(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_fake_bestpixel(monkeypatch)
    config = SIACConfig.model_validate(
        {
            "providers": {
                "monthly_composites": {
                    "kind": "bestpixel",
                    "bestpixel_endpoint": "earth-search",
                    "bestpixel_lookback_years": 3,
                    "bestpixel_bands": ["blue", "green", "red"],
                    "bestpixel_months": [6, 7, 8],
                }
            }
        }
    )
    provider = resolve_monthly_composite_provider(config)
    assert isinstance(provider, BestPixelMonthlyCompositeProvider)
    assert provider.source_name == "bestpixel:earth-search monthly composites"
    assert [b.name for b in provider.source_bands] == ["blue", "green", "red"]

    observation, _x, _y = _observation()
    collection = provider.get_monthly_composites(observation, _RES)
    # 3 years x 3 months.
    assert len(collection.composites) == 9
    assert {c.month for c in collection.composites} == {6, 7, 8}


def test_fetch_resolution_fetches_fine_then_area_averages(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen: dict[str, Any] = {}
    _install_fake_bestpixel(monkeypatch, seen=seen)
    observation, x, y = _observation()  # prior grid at _RES (60 m)

    provider = BestPixelMonthlyCompositeProvider(
        endpoint="pc",
        lookback_years=1,
        resolution_m=_RES,  # prior / database grid = 60 m
        fetch_resolution_m=30.0,  # fetch finer, then average down to 60 m
    )
    collection = provider.get_monthly_composites(observation, _RES)

    # bestpixel was asked to FETCH at the fine resolution (30 m), not the prior
    # resolution — the adapter area-averages from 30 m down to the 60 m grid.
    assert seen["resolution"] == 30.0

    # The composite still lands on the observation's 60 m prior grid.
    composite = collection.composites[0]
    assert composite.reflectance.shape == (len(provider.source_bands), y.size, x.size)
    np.testing.assert_allclose(composite.reflectance.coords["x"].values, x)
    assert np.isfinite(composite.reflectance.values).any()


def test_fetch_resolution_none_fetches_at_prior_resolution(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen: dict[str, Any] = {}
    _install_fake_bestpixel(monkeypatch, seen=seen)
    observation, _x, _y = _observation()

    provider = BestPixelMonthlyCompositeProvider(endpoint="pc", lookback_years=1, resolution_m=_RES)
    provider.get_monthly_composites(observation, _RES)
    # No fetch override -> fetch directly at the prior resolution (no averaging).
    assert seen["resolution"] == _RES


def _install_flaky_bestpixel(monkeypatch: pytest.MonkeyPatch) -> None:
    """A fake bestpixel that errors on Aprils and returns nodata for Februaries."""
    width, height = 6, 5
    xmin, ymax = 300000.0, 8200000.0

    def _build_composite(
        *,
        bbox: Sequence[float],
        datetime: str,
        resolution: float,
        bands: Sequence[str],
        **kwargs: Any,
    ) -> dict[str, Any]:
        month = int(datetime[5:7])
        if month == 4:
            raise RuntimeError("STAC search non-2xx: 400 Bad Request")
        transform = [resolution, 0.0, xmin, 0.0, -resolution, ymax]
        all_nodata = month == 2
        band_arrays = {
            name: np.full(
                (height, width),
                65535 if all_nodata else 1500,
                dtype=np.uint16,
            )
            for name in bands
        }
        quality = np.full((height, width), 65535 if all_nodata else 0, dtype=np.uint16)
        return {
            "bands": band_arrays,
            "quality": quality,
            "grid": {
                "bounds": [xmin, ymax - height * resolution, xmin + width * resolution, ymax],
                "epsg": 32633,
                "crs": _CRS,
                "resolution": resolution,
                "width": width,
                "height": height,
                "transform": transform,
            },
            "band_names": list(bands),
            "source_ids": [] if all_nodata else ["FAKE"],
        }

    fake = types.ModuleType("bestpixel")
    fake.build_composite = _build_composite  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "bestpixel", fake)


def test_get_monthly_composites_skips_error_and_nodata_periods(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_flaky_bestpixel(monkeypatch)
    observation, _x, _y = _observation()
    # Explicit full-year window so the error/nodata months are actually
    # requested (the default would only fetch the scene's own month).
    provider = BestPixelMonthlyCompositeProvider(
        lookback_years=1, months=list(range(1, 13))
    )  # 2025, all 12 months

    collection = provider.get_monthly_composites(observation, _RES)

    months = sorted(c.month for c in collection.composites)
    # April (4) errored and February (2) was all-nodata -> both dropped.
    assert 4 not in months
    assert 2 not in months
    assert months == [1, 3, 5, 6, 7, 8, 9, 10, 11, 12]


def test_explicit_months_use_calendar_correct_date_ranges(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Guards the bestpixel day-31 bug workaround: end dates are real."""
    seen: dict[str, Any] = {}
    _install_fake_bestpixel(monkeypatch, seen=seen)
    observation, _x, _y = _observation()
    provider = BestPixelMonthlyCompositeProvider(
        lookback_years=2, months=[2, 4, 7]
    )  # leap + non-leap Feb, a 30-day month, a 31-day month

    provider.get_monthly_composites(observation, _RES)

    ranges = set(seen["datetimes"])
    assert "2024-02-01/2024-02-29" in ranges  # leap February
    assert "2025-02-01/2025-02-28" in ranges  # non-leap February
    assert "2024-04-01/2024-04-30" in ranges  # 30-day month, never -04-31
    assert "2024-07-01/2024-07-31" in ranges  # 31-day month
    assert not any(
        dt.endswith(("-02-30", "-02-31", "-04-31", "-06-31", "-09-31", "-11-31")) for dt in ranges
    )


def test_get_monthly_composites_raises_when_all_periods_fail(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _always_error(**kwargs: Any) -> dict[str, Any]:
        raise RuntimeError("STAC search non-2xx: 400 Bad Request")

    fake = types.ModuleType("bestpixel")
    fake.build_composite = _always_error  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "bestpixel", fake)

    observation, _x, _y = _observation()
    provider = BestPixelMonthlyCompositeProvider(lookback_years=1, months=[7])
    with pytest.raises(RuntimeError, match="no usable composites"):
        provider.get_monthly_composites(observation, _RES)
