"""Phase C — bestpixel surface prior -> surface-driven solver, end-to-end.

Proves the opt-in stack composes through the REAL pipeline pieces with only
the external ``bestpixel`` package mocked (and the MAIAC day-AOD gate injected):

    build_bestpixel_surface_prior  (M3, real)
        -> assemble_grids          (M4, real)  -> SolverInputBundle
        -> resolve_solver(...)     (config-resolved surface_driven solver, real)
        -> SolvedAtmosphere

The band set is kept identity (source_bands == target_bands == the sensor's
solver bands) so no rsrf spectral-library is needed — mirroring the unit tests.
"""

from __future__ import annotations

import dataclasses
import sys
import types
from datetime import datetime
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest
import xarray as xr

from siac.adapters.bestpixel_surface_prior import build_bestpixel_surface_prior
from siac.algorithms.grid.assembler import assemble_grids
from siac.algorithms.solver.multigrid import _log_aot_axis
from siac.app._assembly_solver import resolve_solver
from siac.app._assembly_surface import resolve_surface_prior_provider
from siac.config import SIACConfig
from siac.config.types import SolverMethod, SurfacePriorMethod
from siac.domain import SensorBand, SensorConfig
from siac.runtime import AtmosphericState, GeometryAngles, ObservationBundle, SolvedAtmosphere
from siac.runtime.validation import validate_solver_input_bundle

pytestmark = pytest.mark.integration

SHAPE = (64, 64)
_CRS = "EPSG:32632"
_BOUNDS = (300000.0, 5500000.0, 300640.0, 5500640.0)
_AEROSOL_RES = 60.0
_SOLVER_BANDS = ("B01", "B02")  # default_aerosol_solver_bands for the SYNTH sensor


# ── Synthetic inputs (style of tests/integration/test_e2e_synthetic.py) ──────
def _make_sensor_config() -> SensorConfig:
    return SensorConfig(
        sensor_id="SYNTH",
        satellite_id="TEST",
        bands=(
            SensorBand("B01", 443.0, 20.0, 60.0, 0),
            SensorBand("B02", 490.0, 65.0, 10.0, 1),
            SensorBand("B03", 560.0, 35.0, 10.0, 2),
            SensorBand("B04", 665.0, 30.0, 10.0, 3),
            SensorBand("B08", 842.0, 115.0, 10.0, 4),
            SensorBand("B11", 1610.0, 90.0, 20.0, 5),
        ),
    )


def _make_obs_bundle(*, toa_value: float | None = None) -> ObservationBundle:
    rng = np.random.RandomState(7)
    config = _make_sensor_config()
    toa = xr.Dataset(
        {
            b.name: xr.DataArray(
                (
                    np.full(SHAPE, toa_value, dtype=np.float32)
                    if toa_value is not None
                    else rng.uniform(0.05, 0.4, SHAPE).astype(np.float32)
                ),
                dims=["y", "x"],
            )
            for b in config.bands
        }
    )
    geometry = GeometryAngles(
        sza=xr.DataArray(np.full(SHAPE, 30.0), dims=["y", "x"]),
        saa=xr.DataArray(np.full(SHAPE, 150.0), dims=["y", "x"]),
        vza=xr.DataArray(np.full(SHAPE, 5.0), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(SHAPE, 110.0), dims=["y", "x"]),
    )
    cloud = xr.DataArray(np.zeros(SHAPE, dtype=bool), dims=["y", "x"])
    return ObservationBundle(
        toa=toa,
        geometry=geometry,
        cloud_mask=cloud,
        sensor_config=config,
        metadata={"observation_time": datetime(2023, 7, 15, 10, 30)},
        crs=_CRS,
        bounds=_BOUNDS,
    )


def _make_atmo(aot_value: float = 0.2) -> AtmosphericState:
    return AtmosphericState(
        aot=xr.DataArray(np.full(SHAPE, aot_value), dims=["y", "x"]),
        tcwv=xr.DataArray(np.full(SHAPE, 2.0), dims=["y", "x"]),
        tco3=xr.DataArray(np.full(SHAPE, 0.35), dims=["y", "x"]),
        aot_unc=xr.DataArray(np.full(SHAPE, 0.1), dims=["y", "x"]),
        tcwv_unc=xr.DataArray(np.full(SHAPE, 0.5), dims=["y", "x"]),
        tco3_unc=xr.DataArray(np.full(SHAPE, 0.02), dims=["y", "x"]),
        elevation=xr.DataArray(np.zeros(SHAPE), dims=["y", "x"]),
    )


# ── Fake bestpixel composites on the observation AOI grid ────────────────────
def _period_dict(band_dn: dict[str, int], *, unc_dn: int = 200) -> dict[str, Any]:
    # 11x11 at 60 m over the observation AOI (matches _build_target_template).
    width = height = 11
    xmin, ymax = _BOUNDS[0], _BOUNDS[3]
    transform = [_AEROSOL_RES, 0.0, xmin, 0.0, -_AEROSOL_RES, ymax]
    return {
        "bands": {
            name: np.full((height, width), dn, dtype=np.uint16) for name, dn in band_dn.items()
        },
        "boa_unc": {
            name: np.full((height, width), unc_dn, dtype=np.float32) for name in band_dn
        },
        "observation_count": np.full((height, width), 5, dtype=np.uint16),
        "quality": np.zeros((height, width), dtype=np.uint16),
        "grid": {
            "bounds": [xmin, ymax - height * _AEROSOL_RES, xmin + width * _AEROSOL_RES, ymax],
            "epsg": 32632,
            "crs": _CRS,
            "resolution": _AEROSOL_RES,
            "width": width,
            "height": height,
            "transform": transform,
        },
    }


def _install_fake_bestpixel(
    monkeypatch: pytest.MonkeyPatch, periods: list[dict[str, Any]]
) -> None:
    def _build_monthly_composites(**_kwargs: Any) -> list[dict[str, Any]]:
        return periods

    fake = types.ModuleType("bestpixel")
    fake.build_monthly_composites = _build_monthly_composites  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "bestpixel", fake)


def _identity_solver_bands(obs: ObservationBundle) -> tuple[SensorBand, ...]:
    # source == target == the sensor's solver bands -> SpectralMapper identity,
    # no rsrf library needed; prior band coords match the solver bands by name.
    return tuple(obs.sensor_config.get_band(name) for name in _SOLVER_BANDS)


def _prior_config(robust_clip: float = 0.0) -> SimpleNamespace:
    return SimpleNamespace(
        algorithms=SimpleNamespace(
            surface_prior=SimpleNamespace(
                bestpixel_robust_clip=robust_clip,
                bestpixel_aod_max=None,
                bestpixel_low_aod_frac=0.6,
            )
        ),
        providers=SimpleNamespace(
            monthly_composites=SimpleNamespace(
                bestpixel_endpoint="pc",
                bestpixel_bands=None,
                bestpixel_lookback_years=3,
                bestpixel_months=None,
                bestpixel_top_k=3,
                bestpixel_max_cloud_cover=90.0,
                bestpixel_disk_cache=None,
            )
        ),
    )


def _build_prior(
    obs: ObservationBundle,
    periods: list[dict[str, Any]],
    *,
    robust_clip: float = 0.0,
) -> Any:
    bands = _identity_solver_bands(obs)
    return build_bestpixel_surface_prior(
        _prior_config(robust_clip),
        obs,
        _AEROSOL_RES,
        bands=_SOLVER_BANDS,
        source_bands=bands,
        target_bands=bands,
        spectral_library=None,
        k_neighbors=5,
        maiac_day_aod=lambda *_a, **_k: {},  # empty gate -> keep all days
    )


def _surface_driven_config(**solver_overrides: Any) -> SIACConfig:
    solver = {
        "method": "surface_driven",
        "aerosol_resolution": _AEROSOL_RES,
        "grid_search_aot_points": 11,
        "surface_driven_backstop_calibrated": True,
    }
    solver.update(solver_overrides)
    return SIACConfig.model_validate({"algorithms": {"solver": solver}})


# ── Dummy linear RT for the controlled-recovery case (boa = toa - aot) ───────
class _DummyLinearRT:
    """RT whose correction gives corrected_boa = toa - aot exactly.

    With ``xap=1, xbp=aot, xcp=0`` the solver's BOA reduction
    ``boa = (xap*toa - xbp) / (1 + xcp*(...))`` collapses to ``toa - aot``.
    """

    def compute_coefficients(
        self, geometry: Any, atmo_state: Any, band: Any, compute_jacobian: bool = False
    ) -> Any:
        from siac.runtime import RTCoefficients

        shape = geometry.sza.shape
        return RTCoefficients(
            xap=xr.DataArray(np.ones(shape, dtype=np.float64), dims=["y", "x"]),
            xbp=xr.DataArray(np.asarray(atmo_state.aot.values, dtype=np.float64), dims=["y", "x"]),
            xcp=xr.DataArray(np.zeros(shape, dtype=np.float64), dims=["y", "x"]),
            d_xap=None,
            d_xbp=None,
            d_xcp=None,
        )

    def supports_jacobian(self) -> bool:
        return False

    @property
    def backend_name(self) -> str:
        return "dummy-linear"

    def is_available_for_sensor(self, sensor_id: Any, satellite_id: Any) -> bool:
        return True


# ── Tests ────────────────────────────────────────────────────────────────────
def test_bestpixel_surface_driven_smoke(
    monkeypatch: pytest.MonkeyPatch, mock_rt_model: Any
) -> None:
    """Full stack composes; SolvedAtmosphere is valid and the mask is respected."""
    periods = [
        _period_dict({"B01": 1000, "B02": 1200}),
        _period_dict({"B01": 1500, "B02": 1300}),
        _period_dict({"B01": 2000, "B02": 1400}),
    ]
    _install_fake_bestpixel(monkeypatch, periods)
    obs = _make_obs_bundle()

    # Cloud out a block so we can assert the masked-pixel prior fallback.
    cloud = np.zeros(SHAPE, dtype=bool)
    cloud[10:40, 10:40] = True
    obs = dataclasses.replace(obs, cloud_mask=xr.DataArray(cloud, dims=["y", "x"]))

    atmo = _make_atmo(aot_value=0.2)
    prior = _build_prior(obs, periods)

    # Sanity: the prior is banded on the solver bands with finite pixels.
    assert prior.boa.dims == ("band", "y", "x")
    assert list(np.asarray(prior.boa.coords["band"].values)) == list(_SOLVER_BANDS)
    assert bool(np.isfinite(prior.boa.values).all())

    sib = assemble_grids(obs, atmo, prior, mock_rt_model, aerosol_resolution_m=_AEROSOL_RES)
    validate_solver_input_bundle(sib)

    config = _surface_driven_config()
    solver_fn = resolve_solver(config)
    solved = solver_fn(sib, config)

    assert isinstance(solved, SolvedAtmosphere)
    # n_iterations == 1 is the surface-driven sweep's signature (not multigrid).
    assert solved.n_iterations == 1
    assert isinstance(solved.converged, bool)

    aot = np.asarray(solved.aot.values)
    assert aot.shape == sib.cloud_mask.shape
    assert np.isfinite(aot).all()
    lo, hi = config.algorithms.solver.aot_bounds
    assert np.all(aot >= lo - 1e-6) and np.all(aot <= hi + 1e-6)

    # Masked (cloudy) solver pixels fall back to the atmospheric prior AOT.
    solver_cloud = np.asarray(sib.cloud_mask.values, dtype=bool)
    prior_aot = np.asarray(sib.atmo_prior.aot.values)
    assert solver_cloud.any()
    np.testing.assert_allclose(aot[solver_cloud], prior_aot[solver_cloud], atol=1e-5)


def test_bestpixel_surface_driven_controlled_recovery(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With a linear dummy RT the argmin lands on the planted AOT node."""
    axis = _log_aot_axis(0.001, 2.5, 11)
    true_aot = float(axis[6])  # an interior log-spaced node
    toa_value = 0.3
    prior_boa = toa_value - true_aot
    assert prior_boa > 0.0
    prior_dn = int(round(prior_boa * 10000))

    periods = [_period_dict({"B01": prior_dn, "B02": prior_dn}, unc_dn=100)]
    _install_fake_bestpixel(monkeypatch, periods)

    obs = _make_obs_bundle(toa_value=toa_value)
    atmo = _make_atmo(aot_value=true_aot)  # backstop agrees with the planted node
    prior = _build_prior(obs, periods)

    dummy_rt = _DummyLinearRT()
    sib = assemble_grids(obs, atmo, prior, dummy_rt, aerosol_resolution_m=_AEROSOL_RES)
    validate_solver_input_bundle(sib)

    config = _surface_driven_config()
    solved = resolve_solver(config)(sib, config)

    aot = np.asarray(solved.aot.values)
    assert np.isfinite(aot).all()
    # The argmin must land on the planted node within one log-axis step.
    node_step = float(axis[7] - axis[6])
    np.testing.assert_allclose(np.median(aot), true_aot, atol=node_step)


def test_surface_driven_example_config_resolves(monkeypatch: pytest.MonkeyPatch) -> None:
    """The shipped example config resolves the bestpixel -> surface_driven stack."""
    from pathlib import Path

    config = SIACConfig.from_file(Path("docs/siac-config.surface-driven.example.toml"))

    # config -> provider wiring: enums select the opt-in stack.
    assert config.algorithms.surface_prior.method == SurfacePriorMethod.BESTPIXEL
    assert config.algorithms.solver.method == SolverMethod.SURFACE_DRIVEN

    # resolve_solver returns the surface-driven solver (n_iterations==1 signature).
    periods = [_period_dict({"B01": 1500, "B02": 1500}, unc_dn=100)]
    _install_fake_bestpixel(monkeypatch, periods)
    obs = _make_obs_bundle(toa_value=0.3)
    atmo = _make_atmo(aot_value=0.2)
    prior = _build_prior(obs, periods)
    sib = assemble_grids(obs, atmo, prior, _DummyLinearRT(), aerosol_resolution_m=_AEROSOL_RES)

    solver_fn = resolve_solver(config)
    solved = solver_fn(sib, config)
    assert isinstance(solved, SolvedAtmosphere)
    assert solved.n_iterations == 1

    # The surface-prior provider also resolves from the config (no network here:
    # the MAIAC gate source is built lazily, never invoked).
    surface_fn = resolve_surface_prior_provider(config)
    assert callable(surface_fn)
    assert surface_fn.requires_atmo_prior is False
