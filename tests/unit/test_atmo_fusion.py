"""Fused aerosol-prior provider (max(MAIAC, CAMS) and friends)."""

from __future__ import annotations

from datetime import datetime

import numpy as np
import pytest
import xarray as xr

from siac.adapters.atmo.fusion import FusedAODProvider
from siac.runtime.models import AtmosphericState


def _field(values: object, shape: tuple[int, int] = (2, 2)) -> xr.DataArray:
    return xr.DataArray(np.full(shape, values, dtype=np.float32), dims=("y", "x"))


def _state(aot: xr.DataArray, tcwv: float = 2.0) -> AtmosphericState:
    return AtmosphericState(
        aot=aot,
        tcwv=_field(tcwv, aot.shape),
        tco3=_field(0.3, aot.shape),
        aot_unc=_field(0.1, aot.shape),
        tcwv_unc=_field(0.5, aot.shape),
        tco3_unc=_field(0.05, aot.shape),
        elevation=_field(0.0, aot.shape),
    )


class _Stub:
    def __init__(self, aot: xr.DataArray, name: str, tcwv: float = 2.0) -> None:
        self._aot = aot
        self._name = name
        self._tcwv = tcwv

    @property
    def source_name(self) -> str:
        return self._name

    def get_prior(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        resolution: float,
    ) -> AtmosphericState:
        return _state(self._aot, self._tcwv)


ARGS = ((0.0, 0.0, 1.0, 1.0), "EPSG:4326", datetime(2024, 5, 5, 3, 35), 60.0)


def test_max_fusion_takes_the_higher_source() -> None:
    maiac = _Stub(_field(0.30), "MAIAC")
    cams = _Stub(_field(0.45), "CAMS")

    fused = FusedAODProvider(maiac, [cams], op="max").get_prior(*ARGS)

    assert float(fused.aot.mean()) == pytest.approx(0.45)


def test_max_fusion_keeps_primary_when_it_is_higher() -> None:
    maiac = _Stub(_field(0.60), "MAIAC")
    cams = _Stub(_field(0.22), "CAMS")

    fused = FusedAODProvider(maiac, [cams], op="max").get_prior(*ARGS)

    assert float(fused.aot.mean()) == pytest.approx(0.60)


def test_fusion_debiases_two_low_sources() -> None:
    # The measured motivation: both sources under-read a thick-aerosol truth in
    # the same direction, so their max is never further from truth than either
    # source and is strictly closer than the weaker one. Averaging instead would
    # preserve the shared under-estimate.
    truth = 0.99
    maiac, cams = 0.82, 0.71
    stubs = (_Stub(_field(maiac), "MAIAC"), [_Stub(_field(cams), "CAMS")])

    fused = FusedAODProvider(*stubs, op="max").get_prior(*ARGS)
    averaged = FusedAODProvider(*stubs, op="mean").get_prior(*ARGS)

    fused_error = abs(float(fused.aot.mean()) - truth)
    assert fused_error <= abs(maiac - truth) + 1e-6
    assert fused_error < abs(cams - truth)
    assert fused_error < abs(float(averaged.aot.mean()) - truth)


def test_mean_and_min_operators_are_available() -> None:
    maiac = _Stub(_field(0.30), "MAIAC")
    cams = _Stub(_field(0.50), "CAMS")

    mean = FusedAODProvider(maiac, [cams], op="mean").get_prior(*ARGS)
    minimum = FusedAODProvider(maiac, [cams], op="min").get_prior(*ARGS)

    assert float(mean.aot.mean()) == pytest.approx(0.40)
    assert float(minimum.aot.mean()) == pytest.approx(0.30)


def test_missing_pixels_in_one_source_fall_back_to_the_other() -> None:
    partial = xr.DataArray(
        np.array([[np.nan, 0.9], [0.2, np.nan]], dtype=np.float32), dims=("y", "x")
    )
    fused = FusedAODProvider(
        _Stub(_field(0.30), "MAIAC"), [_Stub(partial, "CAMS")], op="max"
    ).get_prior(*ARGS)

    # Gap pixels keep the primary value rather than becoming NaN.
    assert np.allclose(fused.aot.values, np.array([[0.3, 0.9], [0.3, 0.3]]), atol=1e-6)


def test_non_aod_state_comes_from_the_primary_provider() -> None:
    maiac = _Stub(_field(0.30), "MAIAC", tcwv=1.5)
    cams = _Stub(_field(0.50), "CAMS", tcwv=9.9)

    fused = FusedAODProvider(maiac, [cams], op="max").get_prior(*ARGS)

    assert float(fused.tcwv.mean()) == pytest.approx(1.5)


def test_source_name_reports_the_fusion() -> None:
    provider = FusedAODProvider(_Stub(_field(0.3), "MAIAC"), [_Stub(_field(0.4), "CAMS")])

    assert provider.source_name == "max(MAIAC, CAMS)"


def test_fusion_requires_at_least_one_extra_source() -> None:
    with pytest.raises(ValueError, match="at least one additional AOD source"):
        FusedAODProvider(_Stub(_field(0.3), "MAIAC"), [])


def test_sources_with_different_dimension_names_are_aligned() -> None:
    # Regression: a satellite retrieval returns a projected (y, x) grid while a
    # model returns geographic (latitude, longitude). Concatenating those without
    # aligning first broadcasts them into a 4-D product, which fails downstream
    # resampling instead of producing a fused 2-D field.
    projected = xr.DataArray(
        np.full((3, 4), 0.20, dtype=np.float32),
        dims=("y", "x"),
        coords={"y": [30.0, 20.0, 10.0], "x": [0.0, 10.0, 20.0, 30.0]},
    )
    geographic = xr.DataArray(
        np.full((2, 2), 0.55, dtype=np.float32),
        dims=("latitude", "longitude"),
        coords={"latitude": [30.0, 10.0], "longitude": [0.0, 30.0]},
    )

    fused = FusedAODProvider(
        _Stub(projected, "MAIAC"), [_Stub(geographic, "CAMS")], op="max"
    ).get_prior(*ARGS)

    assert fused.aot.dims == ("y", "x")
    assert fused.aot.shape == projected.shape
    assert np.isfinite(fused.aot.values).all()


def test_sources_on_a_different_grid_are_aligned() -> None:
    coarse = xr.DataArray(
        np.full((1, 1), 0.7, dtype=np.float32),
        dims=("y", "x"),
        coords={"y": [0.5], "x": [0.5]},
    )
    fine = xr.DataArray(
        np.full((2, 2), 0.2, dtype=np.float32),
        dims=("y", "x"),
        coords={"y": [0.25, 0.75], "x": [0.25, 0.75]},
    )
    fused = FusedAODProvider(_Stub(fine, "MAIAC"), [_Stub(coarse, "CAMS")], op="max").get_prior(
        *ARGS
    )

    assert fused.aot.shape == fine.shape
