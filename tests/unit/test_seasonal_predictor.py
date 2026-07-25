"""Tests for seasonal surface-prior prediction helpers."""

from __future__ import annotations

import numpy as np
import pytest
import rioxarray  # noqa: F401
import xarray as xr

from siac.algorithms.surface.seasonal_predictor import (
    _anchor_match_weights,
    _field_on_template,
    _weighted_median,
    seasonal_extra_tree_prior,
)
from siac.domain import SensorBand, SensorConfig
from siac.runtime import (
    AtmosphericState,
    GeometryAngles,
    ObservationBundle,
    RTCoefficients,
    SurfacePrior,
)


def _grid(values: np.ndarray, *, x: np.ndarray, y: np.ndarray) -> xr.DataArray:
    return xr.DataArray(values, dims=("y", "x"), coords={"y": y, "x": x}).rio.write_crs(
        "EPSG:32632"
    )


class _FakeRT:
    def __init__(self) -> None:
        self.calls: list[str] = []

    def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
        _ = (atmo_state, compute_jacobian)
        self.calls.append(band.name)
        xap = xr.full_like(geometry.sza, 0.8)
        return RTCoefficients(xap=xap, xbp=xr.zeros_like(xap), xcp=xr.zeros_like(xap))


def _scene(height: int = 16, width: int = 16):
    x0 = 500000.0
    y1 = 1000.0
    res = 60.0
    x = x0 + (np.arange(width) + 0.5) * res
    y = y1 - (np.arange(height) + 0.5) * res
    xx, yy = np.meshgrid(np.linspace(0.0, 1.0, width), np.linspace(0.0, 1.0, height))

    comp = np.zeros((3, 7, height, width), dtype=np.float32)
    for index in range(comp.shape[0]):
        offset = index * 0.003
        comp[index, 4] = 0.20 + 0.04 * xx + offset
        comp[index, 5] = 0.30 + 0.03 * yy + offset
        comp[index, 6] = 0.25 + 0.02 * (xx + yy) + offset
        comp[index, 0] = 0.03 + 0.10 * xx
        comp[index, 1] = 0.02 + 0.45 * comp[index, 4]
        comp[index, 2] = 0.04 + 0.20 * yy
        comp[index, 3] = 0.03 + 0.35 * comp[index, 5]

    boa = xr.DataArray(
        np.full((2, height, width), 0.1, dtype=np.float32),
        dims=("band", "y", "x"),
        coords={"band": ["B02", "B04"], "y": y, "x": x},
    ).rio.write_crs("EPSG:32632")
    prior = SurfacePrior(
        boa=boa,
        boa_unc=xr.full_like(boa, 0.02),
        kernels=None,
        mask=_grid(np.ones((height, width), dtype=bool), x=x, y=y),
    )
    sensor = SensorConfig(
        sensor_id="MSI",
        satellite_id="S2A",
        bands=(
            SensorBand("B02", 490.0, 65.0, 10.0, 1),
            SensorBand("B04", 665.0, 30.0, 10.0, 3),
            SensorBand("B8A", 865.0, 20.0, 20.0, 4),
            SensorBand("B11", 1610.0, 90.0, 20.0, 5),
            SensorBand("B12", 2190.0, 180.0, 20.0, 6),
        ),
    )
    observation = ObservationBundle(
        toa=xr.Dataset(
            {
                "B8A": _grid(0.24 + 0.03 * xx, x=x, y=y),
                "B11": _grid(0.32 + 0.02 * yy, x=x, y=y),
                "B12": _grid(0.28 + 0.01 * (xx + yy), x=x, y=y),
            }
        ),
        geometry=GeometryAngles(
            sza=_grid(np.full((height, width), np.radians(30.0)), x=x, y=y),
            saa=_grid(np.zeros((height, width)), x=x, y=y),
            vza=_grid(np.full((height, width), np.radians(5.0)), x=x, y=y),
            vaa=_grid(np.full((height, width), np.radians(20.0)), x=x, y=y),
        ),
        cloud_mask=_grid(np.zeros((height, width), dtype=bool), x=x, y=y),
        sensor_config=sensor,
        metadata={},
        crs="EPSG:32632",
        bounds=(float(x.min()), float(y.min()), float(x.max()), float(y.max())),
    )
    atmo = AtmosphericState(
        aot=_grid(np.full((height, width), 0.2), x=x, y=y),
        tcwv=_grid(np.full((height, width), 2.0), x=x, y=y),
        tco3=_grid(np.full((height, width), 0.3), x=x, y=y),
        aot_unc=_grid(np.full((height, width), 0.1), x=x, y=y),
        tcwv_unc=_grid(np.full((height, width), 0.3), x=x, y=y),
        tco3_unc=_grid(np.full((height, width), 0.03), x=x, y=y),
        elevation=_grid(np.zeros((height, width)), x=x, y=y),
    )
    transform = [res, 0.0, x0, 0.0, -res, y1]
    return prior, observation, atmo, comp, transform


def test_anchor_match_weights_favour_the_nearest_realization() -> None:
    target = np.array([[0.2, 0.3, 0.4]], dtype=np.float64)
    historical = np.array(
        [
            [[0.201, 0.299, 0.401]],
            [[0.30, 0.40, 0.50]],
        ],
        dtype=np.float64,
    )

    weights = _anchor_match_weights(historical, target, scale=0.05)

    assert weights[0, 0] == pytest.approx(1.0)
    assert weights[0, 0] > 5.0 * weights[1, 0]


def test_weighted_median_uses_anchor_weights() -> None:
    values = np.array([[0.1], [0.3], [0.5]], dtype=np.float64)
    weights = np.array([[0.8], [0.1], [0.1]], dtype=np.float64)

    assert _weighted_median(values, weights)[0] == pytest.approx(0.1)


def test_anchor_alignment_keeps_spatial_solver_field_without_explicit_crs() -> None:
    _, observation, _, _, _ = _scene(height=8, width=8)
    template = observation.toa["B8A"]
    source = xr.DataArray(
        np.arange(16, dtype=np.float32).reshape(4, 4),
        dims=("y", "x"),
        coords={
            "y": template.y.values[::2],
            "x": template.x.values[::2],
        },
    )

    aligned = _field_on_template(source, template, fallback=-1.0, label="anchor AOD")

    assert aligned.shape == template.shape
    assert np.isfinite(aligned.values).all()
    assert float(np.ptp(aligned.values)) > 0.0


def test_seasonal_extra_tree_prior_corrects_anchor_with_rt_model() -> None:
    prior, observation, atmo, comp, transform = _scene()
    rt_model = _FakeRT()

    out = seasonal_extra_tree_prior(
        prior,
        observation,
        seasonal_composites=comp,
        epsg=32632,
        transform=transform,
        anchor_aot=0.2,
        atmo_prior=atmo,
        rt_model=rt_model,
    )

    assert rt_model.calls == ["B8A", "B11", "B12"]
    assert not out.boa.identical(prior.boa)
    assert float(out.boa.sel(band="B02").mean()) > 0.0


def test_seasonal_extra_tree_prior_preserves_spatial_anchor_state_and_geometry() -> None:
    prior, observation, atmo, comp, transform = _scene()
    yy, xx = np.indices(atmo.aot.shape, dtype=np.float32)
    atmo = AtmosphericState(
        aot=atmo.aot + 0.001 * xx,
        tcwv=atmo.tcwv + 0.01 * yy,
        tco3=atmo.tco3 + 0.0001 * xx,
        aot_unc=atmo.aot_unc + 0.001 * yy,
        tcwv_unc=atmo.tcwv_unc + 0.001 * xx,
        tco3_unc=atmo.tco3_unc + 0.0001 * yy,
        elevation=atmo.elevation + 0.02 * xx,
    )
    geometry = GeometryAngles(
        sza=observation.geometry.sza + 0.001 * yy,
        saa=observation.geometry.saa + 0.002 * xx,
        vza=observation.geometry.vza + 0.001 * xx,
        vaa=observation.geometry.vaa + 0.002 * yy,
    )
    observation = ObservationBundle(
        toa=observation.toa,
        geometry=geometry,
        cloud_mask=observation.cloud_mask,
        sensor_config=observation.sensor_config,
        metadata=observation.metadata,
        crs=observation.crs,
        bounds=observation.bounds,
    )

    class _CaptureRT:
        def __init__(self) -> None:
            self.states: list[AtmosphericState] = []
            self.geometries: list[GeometryAngles] = []

        def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
            _ = (band, compute_jacobian)
            self.states.append(atmo_state)
            self.geometries.append(geometry)
            xap = xr.full_like(geometry.sza, 0.8)
            return RTCoefficients(xap=xap, xbp=xr.zeros_like(xap), xcp=xr.zeros_like(xap))

    rt_model = _CaptureRT()
    seasonal_extra_tree_prior(
        prior,
        observation,
        seasonal_composites=comp,
        epsg=32632,
        transform=transform,
        anchor_aot=0.2,
        anchor_aot_field=atmo.aot,
        atmo_prior=atmo,
        rt_model=rt_model,
    )

    assert len(rt_model.states) == 3
    state = rt_model.states[0]
    geom = rt_model.geometries[0]
    np.testing.assert_allclose(state.aot.values, atmo.aot.values)
    np.testing.assert_allclose(state.tcwv.values, atmo.tcwv.values)
    np.testing.assert_allclose(state.tco3.values, atmo.tco3.values)
    np.testing.assert_allclose(state.elevation.values, atmo.elevation.values)
    np.testing.assert_allclose(geom.sza.values, geometry.sza.values, atol=1e-7)
    np.testing.assert_allclose(geom.saa.values, geometry.saa.values, atol=1e-7)
    np.testing.assert_allclose(geom.vza.values, geometry.vza.values, atol=1e-7)
    np.testing.assert_allclose(geom.vaa.values, geometry.vaa.values, atol=1e-7)


def test_seasonal_extra_tree_prior_scene_mean_geometry_keeps_atmosphere_spatial() -> None:
    prior, observation, atmo, comp, transform = _scene()
    yy, xx = np.indices(atmo.aot.shape, dtype=np.float32)
    atmo = AtmosphericState(
        aot=atmo.aot + 0.001 * xx,
        tcwv=atmo.tcwv + 0.01 * yy,
        tco3=atmo.tco3,
        aot_unc=atmo.aot_unc,
        tcwv_unc=atmo.tcwv_unc,
        tco3_unc=atmo.tco3_unc,
        elevation=atmo.elevation + 0.02 * xx,
    )
    geometry = GeometryAngles(
        sza=observation.geometry.sza + 0.001 * yy,
        saa=observation.geometry.saa + 0.002 * xx,
        vza=observation.geometry.vza + 0.001 * xx,
        vaa=observation.geometry.vaa + 0.002 * yy,
    )
    observation = ObservationBundle(
        toa=observation.toa,
        geometry=geometry,
        cloud_mask=observation.cloud_mask,
        sensor_config=observation.sensor_config,
        metadata=observation.metadata,
        crs=observation.crs,
        bounds=observation.bounds,
    )

    class _CaptureRT:
        def __init__(self) -> None:
            self.states: list[AtmosphericState] = []
            self.geometries: list[GeometryAngles] = []

        def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
            _ = (band, compute_jacobian)
            self.states.append(atmo_state)
            self.geometries.append(geometry)
            xap = xr.full_like(geometry.sza, 0.8)
            return RTCoefficients(xap=xap, xbp=xr.zeros_like(xap), xcp=xr.zeros_like(xap))

    rt_model = _CaptureRT()
    seasonal_extra_tree_prior(
        prior,
        observation,
        seasonal_composites=comp,
        epsg=32632,
        transform=transform,
        anchor_aot=0.2,
        anchor_aot_field=atmo.aot,
        atmo_prior=atmo,
        rt_model=rt_model,
        scene_mean_geometry=True,
    )

    state = rt_model.states[0]
    geom = rt_model.geometries[0]
    assert float(np.ptp(state.aot.values)) > 0.0
    assert float(np.ptp(state.tcwv.values)) > 0.0
    assert float(np.ptp(state.elevation.values)) > 0.0
    for field in (geom.sza, geom.saa, geom.vza, geom.vaa):
        assert float(np.ptp(field.values)) == pytest.approx(0.0)
    assert float(geom.sza.values[0, 0]) == pytest.approx(float(geometry.sza.mean()))
    assert float(geom.vza.values[0, 0]) == pytest.approx(float(geometry.vza.mean()))


def test_seasonal_extra_tree_prior_requires_rt_model_and_atmo_prior() -> None:
    prior, observation, atmo, comp, transform = _scene()

    with pytest.raises(ValueError, match="RT backend"):
        seasonal_extra_tree_prior(
            prior,
            observation,
            seasonal_composites=comp,
            epsg=32632,
            transform=transform,
            anchor_aot=0.2,
            atmo_prior=atmo,
            rt_model=None,
        )
    with pytest.raises(ValueError, match="atmospheric prior"):
        seasonal_extra_tree_prior(
            prior,
            observation,
            seasonal_composites=comp,
            epsg=32632,
            transform=transform,
            anchor_aot=0.2,
            atmo_prior=None,
            rt_model=_FakeRT(),
        )


def test_seasonal_extra_tree_prior_applies_affine_debias() -> None:
    prior, observation, atmo, comp, transform = _scene()
    common = {
        "seasonal_composites": comp,
        "epsg": 32632,
        "transform": transform,
        "anchor_aot": 0.4,
        "atmo_prior": atmo,
        "rt_model": _FakeRT(),
    }
    base = seasonal_extra_tree_prior(prior, observation, **common)
    shifted = seasonal_extra_tree_prior(
        prior,
        observation,
        debias={"B02": (0.05, 0.1)},
        **common,
    )
    delta = float((shifted.boa.sel(band="B02") - base.boa.sel(band="B02")).median())
    assert delta == pytest.approx(0.05 + 0.1 * 0.4, abs=1e-6)
    delta_b04 = float((shifted.boa.sel(band="B04") - base.boa.sel(band="B04")).median())
    assert delta_b04 == pytest.approx(0.0, abs=1e-6)


def test_seasonal_extra_tree_prior_can_blend_toward_composite_reference() -> None:
    prior, observation, atmo, comp, transform = _scene()

    blended = seasonal_extra_tree_prior(
        prior,
        observation,
        seasonal_composites=comp,
        epsg=32632,
        transform=transform,
        anchor_aot=0.2,
        atmo_prior=atmo,
        rt_model=_FakeRT(),
        composite_blend_weight=1.0,
    )

    expected = float(np.median(comp[:, 1, 8, 8]))
    actual = float(blended.boa.sel(band="B02").isel(y=8, x=8))
    assert actual == pytest.approx(expected, rel=0.05)


def test_seasonal_extra_tree_prior_attaches_tau_predictor_payload() -> None:
    prior, observation, atmo, comp, transform = _scene()
    out = seasonal_extra_tree_prior(
        prior,
        observation,
        seasonal_composites=comp,
        epsg=32632,
        transform=transform,
        anchor_aot=0.2,
        atmo_prior=atmo,
        rt_model=_FakeRT(),
        debias={"B02": (0.01, 0.02)},
        attach_tau_predictor=True,
    )
    payload = out.tau_predictor
    assert payload is not None
    assert len(payload["trees"]) == comp.shape[0]
    assert payload["anchor_bands"] == ("B8A", "B11", "B12")
    assert payload["target_bands"] == ("B02", "B04")
    assert payload["localizer"].shape == (4, 16, 16)
    assert payload["debias"]["B02"] == (0.01, 0.02)
    # default off keeps the field empty
    out_off = seasonal_extra_tree_prior(
        prior,
        observation,
        seasonal_composites=comp,
        epsg=32632,
        transform=transform,
        anchor_aot=0.2,
        atmo_prior=atmo,
        rt_model=_FakeRT(),
    )
    assert out_off.tau_predictor is None
