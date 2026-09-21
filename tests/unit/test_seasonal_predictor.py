"""Tests for seasonal surface-prior prediction helpers."""

from __future__ import annotations

import numpy as np
import pytest
import rioxarray  # noqa: F401
import xarray as xr

from siac.algorithms.surface.seasonal_predictor import (
    _anchor_match_weights,
    _field_on_template,
    _robust_clip_composites,
    _weighted_median,
    predict_visible_from_tau_payload,
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


def test_robust_clip_preserves_only_a_recurrent_snow_mode() -> None:
    comp = np.full((10, 7, 2, 2), 0.1, dtype=np.float32)
    comp[:, 5] = 0.25
    # Snow recurs in four realizations at [0,0], but only once at [0,1].
    comp[:4, 2, 0, 0] = 0.7
    comp[:4, 4, 0, 0] = 0.7
    comp[:4, 5, 0, 0] = 0.07
    comp[:4, 6, 0, 0] = 0.06
    comp[:1, 2, 0, 1] = 0.7
    comp[:1, 4, 0, 1] = 0.7
    comp[:1, 5, 0, 1] = 0.07
    comp[:1, 6, 0, 1] = 0.06

    clipped = _robust_clip_composites(
        comp,
        1.5,
        preserve_recurrent_snow_fraction=0.2,
    )

    assert np.isfinite(clipped[:4, :, 0, 0]).all()
    assert np.isnan(clipped[0, 2, 0, 1])


def test_tau_prediction_uses_archived_anchor_aggregation_weights() -> None:
    class ConstantTree:
        def __init__(self, value: float) -> None:
            self.value = value

        def predict(self, features):
            return np.full((features.shape[0], 1), self.value, dtype=np.float64)

    shape = (8, 8)
    weights = np.stack(
        [
            np.full(shape, 0.9, dtype=np.float64),
            np.full(shape, 0.1, dtype=np.float64),
        ]
    )
    predicted = predict_visible_from_tau_payload(
        np.full((1, *shape), 0.2, dtype=np.float64),
        band_names=("B02",),
        tau_payload={
            "localizer_grid": np.full((1, *shape), 0.1, dtype=np.float64),
            "trees": [ConstantTree(0.7), ConstantTree(0.1)],
            "target_bands": ("B02",),
            "debias": {},
            "debias_scale": 1.0,
            "aggregation_weights_grid": weights,
        },
        anchor_boa=np.full((3, *shape), 0.4, dtype=np.float64),
        aot=0.1,
    )

    assert np.allclose(predicted, 0.7)


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


def test_relative_uncertainty_floor_defaults_to_the_committed_absolute_floor() -> None:
    prior, observation, atmo, comp, transform = _scene()
    common = {
        "seasonal_composites": comp,
        "epsg": 32632,
        "transform": transform,
        "anchor_aot": 0.4,
        "atmo_prior": atmo,
        "rt_model": _FakeRT(),
    }
    committed = seasonal_extra_tree_prior(prior, observation, **common)
    explicit_zero = seasonal_extra_tree_prior(
        prior, observation, relative_uncertainty_floor=0.0, **common
    )
    np.testing.assert_array_equal(committed.boa_unc.values, explicit_zero.boa_unc.values)
    # Every finite value still sits at or above the committed 0.006 floor.
    finite = committed.boa_unc.values[np.isfinite(committed.boa_unc.values)]
    assert finite.size
    assert float(finite.min()) >= 0.006 - 1e-9


def test_relative_uncertainty_floor_tracks_predicted_reflectance() -> None:
    prior, observation, atmo, comp, transform = _scene()
    common = {
        "seasonal_composites": comp,
        "epsg": 32632,
        "transform": transform,
        "anchor_aot": 0.4,
        "atmo_prior": atmo,
        "rt_model": _FakeRT(),
    }
    # A small absolute guard plus a proportional term: the floor should follow
    # the prediction instead of pinning every ordinary pixel to one constant.
    relative = seasonal_extra_tree_prior(
        prior,
        observation,
        uncertainty_floor=0.001,
        relative_uncertainty_floor=0.25,
        **common,
    )
    reflectance = np.abs(relative.boa.values)
    uncertainty = relative.boa_unc.values
    finite = np.isfinite(reflectance) & np.isfinite(uncertainty)
    assert finite.any()
    # Floor is respected everywhere ...
    assert np.all(
        uncertainty[finite] >= np.minimum(0.25 * reflectance[finite], uncertainty[finite]) - 1e-9
    )
    assert np.all(uncertainty[finite] >= 0.001 - 1e-9)
    # ... and it actually binds, so sigma is no longer a single constant.
    proportional = np.isclose(uncertainty[finite], 0.25 * reflectance[finite], atol=1e-9)
    assert proportional.any(), "relative floor never bound"
    assert float(np.std(uncertainty[finite])) > 0.0


def test_pooled_fit_trains_one_model_over_every_realization() -> None:
    # The per-realization ensemble bounds each member by the values in one
    # composite, so the median across members collapses toward the seasonal
    # middle. A pooled fit sees every composite at once; the payload keeps the
    # forest's member trees so the tau path and the aggregation agree.
    prior, observation, atmo, comp, transform = _scene()
    common = {
        "seasonal_composites": comp,
        "epsg": 32632,
        "transform": transform,
        "anchor_aot": 0.4,
        "atmo_prior": atmo,
        "rt_model": _FakeRT(),
        "attach_tau_predictor": True,
    }

    pooled = seasonal_extra_tree_prior(
        prior, observation, predictor_model="extra_trees_20_pooled", **common
    )
    per_realization = seasonal_extra_tree_prior(
        prior, observation, predictor_model="extra_trees_20", **common
    )

    assert len(pooled.tau_predictor["trees"]) == 20
    assert len(per_realization.tau_predictor["trees"]) == comp.shape[0]
    assert np.isfinite(pooled.boa.sel(band="B02").values).all()
    assert not pooled.boa.identical(per_realization.boa)


def test_pooled_fit_still_produces_a_varying_uncertainty() -> None:
    # The sigma changes meaning under pooling -- spread across trees rather than
    # across realizations -- so guard that it does not degenerate to the floor.
    prior, observation, atmo, comp, transform = _scene()

    out = seasonal_extra_tree_prior(
        prior,
        observation,
        seasonal_composites=comp,
        epsg=32632,
        transform=transform,
        anchor_aot=0.4,
        atmo_prior=atmo,
        rt_model=_FakeRT(),
        predictor_model="extra_trees_20_pooled",
        uncertainty_floor=0.001,
    )

    sigma = out.boa_unc.sel(band="B02").values
    assert np.isfinite(sigma).all()
    assert float(sigma.min()) >= 0.001
    assert float(np.ptp(sigma)) > 0.0


def test_pooled_fit_refuses_anchor_weighted_aggregation() -> None:
    # Anchor weights carry one plane per realization; a pooled fit has no such
    # axis, and the tau payload's shape guard would reject the combination.
    prior, observation, atmo, comp, transform = _scene()

    with pytest.raises(ValueError, match="pooled"):
        seasonal_extra_tree_prior(
            prior,
            observation,
            seasonal_composites=comp,
            epsg=32632,
            transform=transform,
            anchor_aot=0.4,
            atmo_prior=atmo,
            rt_model=_FakeRT(),
            predictor_model="extra_trees_20_pooled",
            ensemble_aggregation="anchor_weighted",
        )


def test_pooled_row_cap_bounds_the_training_set() -> None:
    # Pooling multiplies rows by the realization count; the cap keeps fit time
    # and memory comparable to the per-realization path.
    prior, observation, atmo, comp, transform = _scene()

    capped = seasonal_extra_tree_prior(
        prior,
        observation,
        seasonal_composites=comp,
        epsg=32632,
        transform=transform,
        anchor_aot=0.4,
        atmo_prior=atmo,
        rt_model=_FakeRT(),
        predictor_model="extra_trees_20_pooled",
        pooled_max_rows=64,
    )

    assert np.isfinite(capped.boa.sel(band="B02").values).all()


def _common_kwargs(atmo, comp, transform):  # noqa: ANN001, ANN202
    return {
        "seasonal_composites": comp,
        "epsg": 32632,
        "transform": transform,
        "anchor_aot": 0.2,
        "atmo_prior": atmo,
        "rt_model": _FakeRT(),
    }


def test_wide_composite_without_band_names_is_refused() -> None:
    """Columns 4/5/6 of a wide library are red-edge bands, not the anchors.

    The old guard only rejected composites NARROWER than seven bands, so a
    twelve-band library trained on B05/B06/B07 while calling them B8A/B11/B12
    and returned a wrong prior with no error anywhere.
    """
    prior, observation, atmo, comp, transform = _scene()
    wide = np.concatenate([comp, comp[:, :5]], axis=1)

    with pytest.raises(ValueError, match="composite_band_names"):
        seasonal_extra_tree_prior(prior, observation, **_common_kwargs(atmo, wide, transform))


def test_named_wide_composite_matches_the_seven_band_layout() -> None:
    """A wide library must predict exactly what the narrow one does.

    The extra planes are inserted between the visible bands and the anchors,
    which is precisely where positional indexing goes wrong, so an identical
    prediction is what demonstrates the columns are resolved by name.
    """
    prior, observation, atmo, comp, transform = _scene()
    narrow = seasonal_extra_tree_prior(prior, observation, **_common_kwargs(atmo, comp, transform))

    filler = np.full_like(comp[:, :1], 0.5)
    wide = np.concatenate([comp[:, :4], filler, filler, filler, filler, comp[:, 4:]], axis=1)
    names = ["B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]
    widened = seasonal_extra_tree_prior(
        prior,
        observation,
        composite_band_names=names,
        **_common_kwargs(atmo, wide, transform),
    )

    np.testing.assert_allclose(
        np.asarray(widened.boa.values), np.asarray(narrow.boa.values), rtol=0, atol=0
    )


def test_composite_band_name_aliases_resolve() -> None:
    """The seven-band libraries name their planes by role, not by band."""
    prior, observation, atmo, comp, transform = _scene()
    unnamed = seasonal_extra_tree_prior(prior, observation, **_common_kwargs(atmo, comp, transform))
    aliased = seasonal_extra_tree_prior(
        prior,
        observation,
        composite_band_names=["coastal", "blue", "green", "red", "nir", "swir16", "swir22"],
        **_common_kwargs(atmo, comp, transform),
    )

    np.testing.assert_allclose(
        np.asarray(aliased.boa.values), np.asarray(unnamed.boa.values), rtol=0, atol=0
    )


def test_composite_missing_an_anchor_band_is_refused() -> None:
    prior, observation, atmo, comp, transform = _scene()

    with pytest.raises(ValueError, match="missing required bands"):
        seasonal_extra_tree_prior(
            prior,
            observation,
            composite_band_names=["B01", "B02", "B03", "B04", "B8A", "B11", "B05"],
            **_common_kwargs(atmo, comp, transform),
        )


def test_composite_band_names_must_match_the_composite_width() -> None:
    prior, observation, atmo, comp, transform = _scene()

    with pytest.raises(ValueError, match="entries for a"):
        seasonal_extra_tree_prior(
            prior,
            observation,
            composite_band_names=["B01", "B02", "B03"],
            **_common_kwargs(atmo, comp, transform),
        )


def _four_band_prior(prior: SurfacePrior) -> SurfacePrior:
    """The fixture prior widened to B01..B04 so extra bands can be targeted."""
    one = prior.boa.isel(band=0, drop=True)
    boa = xr.concat([one, one, one, one], dim="band").assign_coords(
        band=["B01", "B02", "B03", "B04"]
    )
    return SurfacePrior(boa=boa, boa_unc=xr.full_like(boa, 0.02), kernels=None, mask=prior.mask)


def _visible(result: SurfacePrior) -> np.ndarray:
    return np.asarray(result.boa.sel(band=["B02", "B04"]).values)


@pytest.mark.parametrize(
    "predictor_model", ["extra_tree", "extra_trees_20", "extra_trees_20_pooled"]
)
def test_split_fit_leaves_the_first_group_bit_identical(predictor_model: str) -> None:
    """Adding a separately fitted group must not move the bands already fitted.

    A joint multi-output fit chooses every split for all targets at once, so
    requesting more bands changes the trees of the bands already requested.
    A group fitted on its own must predict exactly what a standalone fit of
    those bands predicts.
    """
    prior, observation, atmo, comp, transform = _scene()
    wide = _four_band_prior(prior)
    kwargs = {**_common_kwargs(atmo, comp, transform), "predictor_model": predictor_model}

    standalone = seasonal_extra_tree_prior(
        wide, observation, target_band_columns={"B02": 1, "B04": 3}, **kwargs
    )
    split = seasonal_extra_tree_prior(
        wide,
        observation,
        target_band_columns={"B02": 1, "B04": 3, "B01": 0, "B03": 2},
        target_fit_groups=[["B02", "B04"], ["B01", "B03"]],
        **kwargs,
    )

    np.testing.assert_allclose(_visible(split), _visible(standalone), rtol=0, atol=0)
    np.testing.assert_allclose(
        np.asarray(split.boa_unc.sel(band=["B02", "B04"]).values),
        np.asarray(standalone.boa_unc.sel(band=["B02", "B04"]).values),
        rtol=0,
        atol=0,
    )


def test_split_fit_still_predicts_the_added_group() -> None:
    prior, observation, atmo, comp, transform = _scene()
    wide = _four_band_prior(prior)

    split = seasonal_extra_tree_prior(
        wide,
        observation,
        target_band_columns={"B02": 1, "B04": 3, "B01": 0, "B03": 2},
        target_fit_groups=[["B02", "B04"], ["B01", "B03"]],
        **_common_kwargs(atmo, comp, transform),
    )

    for band in ("B01", "B03"):
        plane = np.asarray(split.boa.sel(band=band).values)
        assert np.isfinite(plane).all()
        # The carrier plane was a flat 0.1; a fitted prediction varies.
        assert float(np.std(plane)) > 0.0


def test_split_fit_tau_replay_matches_the_standalone_fit() -> None:
    """The teacher's archived surface comes from the tau replay, not the fit."""
    prior, observation, atmo, comp, transform = _scene()
    wide = _four_band_prior(prior)
    kwargs = {**_common_kwargs(atmo, comp, transform), "attach_tau_predictor": True}

    standalone = seasonal_extra_tree_prior(
        wide, observation, target_band_columns={"B02": 1, "B04": 3}, **kwargs
    )
    split = seasonal_extra_tree_prior(
        wide,
        observation,
        target_band_columns={"B02": 1, "B04": 3, "B01": 0, "B03": 2},
        target_fit_groups=[["B02", "B04"], ["B01", "B03"]],
        **kwargs,
    )
    assert split.tau_predictor["fit_groups"] == (("B02", "B04"), ("B01", "B03"))
    assert split.tau_predictor["realizations_excluded_by_later_group"] == 0

    anchor_boa = np.stack(
        [np.asarray(observation.toa[name].values) for name in ("B8A", "B11", "B12")]
    )
    names = ["B01", "B02", "B03", "B04"]
    base = np.asarray(wide.boa.values)

    def replay(result: SurfacePrior) -> np.ndarray:
        payload = {
            **result.tau_predictor,
            "localizer_grid": np.asarray(result.tau_predictor["localizer"]),
        }
        return predict_visible_from_tau_payload(
            base, band_names=names, tau_payload=payload, anchor_boa=anchor_boa, aot=0.2
        )

    for band in ("B02", "B04"):
        index = names.index(band)
        np.testing.assert_allclose(replay(split)[index], replay(standalone)[index], rtol=0, atol=0)


def test_default_fit_groups_keep_the_joint_fit() -> None:
    prior, observation, atmo, comp, transform = _scene()
    wide = _four_band_prior(prior)
    targets = {"B02": 1, "B04": 3, "B01": 0, "B03": 2}
    kwargs = _common_kwargs(atmo, comp, transform)

    implicit = seasonal_extra_tree_prior(wide, observation, target_band_columns=targets, **kwargs)
    explicit = seasonal_extra_tree_prior(
        wide,
        observation,
        target_band_columns=targets,
        target_fit_groups=[["B02", "B04", "B01", "B03"]],
        **kwargs,
    )

    np.testing.assert_allclose(
        np.asarray(explicit.boa.values), np.asarray(implicit.boa.values), rtol=0, atol=0
    )


@pytest.mark.parametrize(
    ("groups", "message"),
    [
        ([["B02", "B04"], ["B04"]], "more than one group"),
        ([["B02", "B04", "B07"]], "not targets"),
        ([["B02"]], "unassigned"),
    ],
)
def test_malformed_fit_groups_are_refused(groups: list[list[str]], message: str) -> None:
    prior, observation, atmo, comp, transform = _scene()

    with pytest.raises(ValueError, match=message):
        seasonal_extra_tree_prior(
            prior,
            observation,
            target_band_columns={"B02": 1, "B04": 3},
            target_fit_groups=groups,
            **_common_kwargs(atmo, comp, transform),
        )
