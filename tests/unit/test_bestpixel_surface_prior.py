"""Tests for the bestpixel-backed surface prior (surface-driven solver).

The external ``bestpixel`` and ``earthaccess`` packages are optional and reach
the network at call time, so every test injects fakes: a fake ``bestpixel``
module in ``sys.modules`` and an injected per-day MAIAC AOD gate callable. The
reduction math, gate logic, and registry dispatch are asserted, never the
network.
"""

from __future__ import annotations

import sys
import types
from datetime import datetime
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any

import numpy as np
import pytest
import rioxarray  # noqa: F401
import xarray as xr

from siac.adapters.bestpixel import bestpixel_source_bands
from siac.adapters.bestpixel_surface_prior import (
    _default_maiac_day_aod,
    _resolve_months,
    _resolve_years,
    build_aod_by_day_gate,
    build_bestpixel_surface_prior,
    direct_source_indices,
    reduce_realizations,
)
from siac.app._assembly_surface import (
    make_bestpixel_surface_prior_fn,
    resolve_surface_prior_provider,
)
from siac.catalog.sensors.sentinel2 import SENTINEL2A_CONFIG
from siac.config import SIACConfig
from siac.runtime import SurfacePrior

if TYPE_CHECKING:
    from collections.abc import Sequence

# A 10x10 grid at 60 m over a UTM 33N AOI (matches _build_target_template).
_CRS = "EPSG:32633"
_RES = 60.0
_BOUNDS = (300000.0, 8199400.0, 300600.0, 8200000.0)
_BANDS = ("coastal", "blue", "red")


# --------------------------------------------------------------------------- #
# reduce_realizations
# --------------------------------------------------------------------------- #
def test_reduce_realizations_median_and_rmse() -> None:
    # Three realizations of a single (band=1, y=1, x=1) pixel: 0.1, 0.2, 0.3.
    stack = np.array([0.1, 0.2, 0.3], dtype=np.float32).reshape(3, 1, 1, 1)
    boa, boa_unc = reduce_realizations(stack, None)
    # median = 0.2; boa_unc = RMSE about it = sqrt(mean([0.01, 0, 0.01])) = sqrt(0.02/3)
    # (matches the comp_ref harness; NOT the MAD estimator).
    np.testing.assert_allclose(boa[0, 0, 0], 0.2, atol=1e-5)
    np.testing.assert_allclose(boa_unc[0, 0, 0], np.sqrt(0.02 / 3.0), atol=1e-5)


def test_reduce_realizations_single_realization_uses_per_composite_unc() -> None:
    stack = np.array([0.25], dtype=np.float32).reshape(1, 1, 1, 1)
    unc_first = np.array([0.04], dtype=np.float32).reshape(1, 1, 1)
    boa, boa_unc = reduce_realizations(stack, unc_first)
    np.testing.assert_allclose(boa[0, 0, 0], 0.25, atol=1e-6)
    np.testing.assert_allclose(boa_unc[0, 0, 0], 0.04, atol=1e-6)


def test_reduce_realizations_single_without_unc_floors() -> None:
    stack = np.array([0.25], dtype=np.float32).reshape(1, 1, 1, 1)
    _boa, boa_unc = reduce_realizations(stack, None, floor=0.006)
    np.testing.assert_allclose(boa_unc[0, 0, 0], 0.006, atol=1e-6)


def test_reduce_realizations_floor_applied() -> None:
    # Identical realizations -> MAD 0 -> floored.
    stack = np.full((4, 1, 1, 1), 0.3, dtype=np.float32)
    _boa, boa_unc = reduce_realizations(stack, None, floor=0.006)
    np.testing.assert_allclose(boa_unc[0, 0, 0], 0.006, atol=1e-6)


def test_reduce_realizations_robust_clip_drops_outlier() -> None:
    # Tight cluster {0.10, 0.11, 0.12} plus a far outlier 1.0.
    stack = np.array([0.10, 0.11, 0.12, 1.0], dtype=np.float32).reshape(4, 1, 1, 1)
    boa_clip, _ = reduce_realizations(stack, None, robust_clip=2.0)
    boa_noclip, _ = reduce_realizations(stack, None, robust_clip=0.0)
    # With clipping the outlier is removed, pulling the median back to the cluster.
    assert boa_clip[0, 0, 0] < boa_noclip[0, 0, 0]
    assert 0.10 <= boa_clip[0, 0, 0] <= 0.12


def test_reduce_realizations_nan_boa_keeps_nan_unc() -> None:
    stack = np.full((3, 1, 1, 1), np.nan, dtype=np.float32)
    boa, boa_unc = reduce_realizations(stack, None)
    assert np.isnan(boa[0, 0, 0])
    assert np.isnan(boa_unc[0, 0, 0])


def test_reduce_realizations_rejects_bad_shape() -> None:
    with pytest.raises(ValueError, match="must be"):
        reduce_realizations(np.zeros((3, 1, 1), dtype=np.float32), None)


# --------------------------------------------------------------------------- #
# build_aod_by_day_gate
# --------------------------------------------------------------------------- #
def test_build_aod_by_day_gate_empty_map_disables_gate() -> None:
    gate, eff = build_aod_by_day_gate({}, aod_max=None, low_aod_frac=0.6)
    assert gate is None
    assert eff is None


def test_build_aod_by_day_gate_keeps_lowest_fraction_per_window() -> None:
    raw = {
        "2021-03-01": 0.1,
        "2021-03-05": 0.5,
        "2021-03-09": 0.2,
        "2021-03-15": 0.9,
        "2022-03-02": 0.3,
        "2022-03-08": 0.05,
    }
    gate, eff = build_aod_by_day_gate(raw, aod_max=None, low_aod_frac=0.5)
    assert gate is not None
    # 2021-03 has 4 days -> keep ceil(0.5*4)=2 lowest (0.1, 0.2).
    assert "2021-03-01" in gate
    assert "2021-03-09" in gate
    assert "2021-03-05" not in gate
    assert "2021-03-15" not in gate
    # 2022-03 has 2 days -> keep 1 lowest (0.05).
    assert "2022-03-08" in gate
    assert "2022-03-02" not in gate
    # effective aod_max covers the highest kept day (so bestpixel won't trim it).
    assert eff is not None and eff > max(gate.values())


def test_build_aod_by_day_gate_aod_max_threshold() -> None:
    raw = {"2021-03-01": 0.1, "2021-03-05": 0.5, "2021-03-09": 0.2}
    gate, eff = build_aod_by_day_gate(raw, aod_max=0.25, low_aod_frac=0.6)
    assert gate == {"2021-03-01": 0.1, "2021-03-09": 0.2}
    assert eff == 0.25


def test_build_aod_by_day_gate_threshold_removes_all_disables_gate() -> None:
    raw = {"2021-03-01": 0.9}
    gate, eff = build_aod_by_day_gate(raw, aod_max=0.1, low_aod_frac=0.6)
    assert gate is None
    assert eff is None


# --------------------------------------------------------------------------- #
# build_bestpixel_surface_prior (full reduction with a fake bestpixel)
# --------------------------------------------------------------------------- #
def _fake_config(*, robust_clip: float = 0.0, lookback_years: int = 2) -> SimpleNamespace:
    surface_prior = SimpleNamespace(
        bestpixel_source="l2a",
        bestpixel_aod_max=None,
        bestpixel_low_aod_frac=0.6,
        bestpixel_robust_clip=robust_clip,
        bestpixel_window_reduction="window",
    )
    monthly_composites = SimpleNamespace(
        bestpixel_endpoint="pc",
        bestpixel_bands=None,
        bestpixel_lookback_years=lookback_years,
        bestpixel_months=None,
        bestpixel_top_k=3,
        bestpixel_max_cloud_cover=90.0,
        bestpixel_disk_cache=None,
    )
    return SimpleNamespace(
        algorithms=SimpleNamespace(surface_prior=surface_prior),
        providers=SimpleNamespace(monthly_composites=monthly_composites),
    )


def _fake_observation() -> SimpleNamespace:
    return SimpleNamespace(
        bounds=_BOUNDS,
        crs=_CRS,
        metadata={"observation_time": datetime(2026, 3, 15, 9, 30)},
    )


def _period_dict(
    bands: Sequence[str],
    *,
    band_dn: dict[str, int],
    unc_dn: int = 500,
    width: int = 10,
    height: int = 10,
) -> dict[str, Any]:
    xmin, ymax = 300000.0, 8200000.0
    transform = [_RES, 0.0, xmin, 0.0, -_RES, ymax]
    return {
        "bands": {name: np.full((height, width), band_dn[name], dtype=np.uint16) for name in bands},
        "boa_unc": {name: np.full((height, width), unc_dn, dtype=np.float32) for name in bands},
        "observation_count": np.full((height, width), 5, dtype=np.uint16),
        "quality": np.zeros((height, width), dtype=np.uint16),
        "grid": {
            "bounds": [xmin, ymax - height * _RES, xmin + width * _RES, ymax],
            "epsg": 32633,
            "crs": _CRS,
            "resolution": _RES,
            "width": width,
            "height": height,
            "transform": transform,
        },
    }


def _install_fake_bestpixel(
    monkeypatch: pytest.MonkeyPatch,
    periods: list[dict[str, Any]],
    *,
    seen: dict[str, Any] | None = None,
) -> None:
    def _build_monthly_composites(**kwargs: Any) -> list[dict[str, Any]]:
        if seen is not None:
            seen["kwargs"] = dict(kwargs)
        return periods

    fake = types.ModuleType("bestpixel")
    fake.build_monthly_composites = _build_monthly_composites  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "bestpixel", fake)


def _install_fake_daily_bestpixel(
    monkeypatch: pytest.MonkeyPatch,
    *,
    seen: dict[str, Any] | None = None,
) -> None:
    def _build_composite(**kwargs: Any) -> dict[str, Any]:
        if seen is not None:
            seen.setdefault("calls", []).append(dict(kwargs))
        day = str(kwargs["datetime"]).split("/", 1)[0]
        coastal_dn = 1000 if day.endswith("-01") else 3000
        return _period_dict(
            tuple(kwargs["bands"]),
            band_dn={"coastal": coastal_dn, "blue": 800, "red": 1200},
        )

    fake = types.ModuleType("bestpixel")
    fake.build_composite = _build_composite  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "bestpixel", fake)


def _identity_bands() -> Any:
    # source == target so SpectralMapper takes the identity path (no library).
    bands = bestpixel_source_bands(_BANDS, endpoint="pc")
    return bands


def test_direct_source_indices_identity_for_s2_subset() -> None:
    # 7-band bestpixel source basis; S2 solve-band subset target (B01,B02,B04).
    source = bestpixel_source_bands(
        ("coastal", "blue", "green", "red", "nir", "swir16", "swir22"), endpoint="pc"
    )
    target = [SENTINEL2A_CONFIG.get_band(n) for n in ("B01", "B02", "B04")]
    # coastal=0, blue=1, red=3 in the source basis (matched by RSRF band id).
    assert direct_source_indices(source, target) == [0, 1, 3]


def test_direct_source_indices_none_for_cross_sensor() -> None:
    # mcd43a4 MODIS NBAR source basis has no 443 nm coastal -> B01 has no
    # direct counterpart, so identity passthrough is declined (KNN fallback).
    source = bestpixel_source_bands(("blue", "red"), endpoint="mcd43a4")
    target = [SENTINEL2A_CONFIG.get_band(n) for n in ("B01", "B02", "B04")]
    assert direct_source_indices(source, target) is None


def test_build_surface_prior_identity_passthrough_no_knn(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # S2 L2A composite -> S2 solve bands: the prior must use the composite bands
    # DIRECTLY (coastal->B01, blue->B02, red->B04), with NO KNN distortion, even
    # though source(3) and target(3) have different *names*. spectral_library is
    # None: if the KNN path ran it would need a library, so identity is required.
    periods = [_period_dict(_BANDS, band_dn={"coastal": 1000, "blue": 800, "red": 1200})]
    _install_fake_bestpixel(monkeypatch, periods)
    source = bestpixel_source_bands(_BANDS, endpoint="pc")
    target = [SENTINEL2A_CONFIG.get_band(n) for n in ("B01", "B02", "B04")]

    prior = build_bestpixel_surface_prior(
        _fake_config(),
        _fake_observation(),
        _RES,
        bands=_BANDS,
        source_bands=source,
        target_bands=target,
        spectral_library=None,
        k_neighbors=5,
        maiac_day_aod=lambda *_a, **_k: {},
    )

    assert list(prior.boa["band"].values) == ["B01", "B02", "B04"]
    # Identity: prior reflectance == composite DN * scale, no remap mixing.
    np.testing.assert_allclose(np.nanmedian(prior.boa.sel(band="B01").values), 0.1, atol=1e-5)
    np.testing.assert_allclose(np.nanmedian(prior.boa.sel(band="B02").values), 0.08, atol=1e-5)
    np.testing.assert_allclose(np.nanmedian(prior.boa.sel(band="B04").values), 0.12, atol=1e-5)


def test_build_surface_prior_predicts_visible_with_uncertainty_floor(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bands = ("coastal", "blue", "green", "red", "nir", "swir16", "swir22")
    periods = [
        _period_dict(
            bands,
            band_dn={
                "coastal": 700,
                "blue": 900,
                "green": 1000,
                "red": 1300,
                "nir": 2600,
                "swir16": 3100,
                "swir22": 2900,
            },
            width=20,
            height=20,
        )
        for _ in range(3)
    ]
    _install_fake_bestpixel(monkeypatch, periods)
    config = _fake_config(lookback_years=3)
    config.algorithms.surface_prior.bestpixel_predict_visible = True
    config.algorithms.surface_prior.bestpixel_predict_visible_bands = ("B02", "B04")
    config.algorithms.surface_prior.bestpixel_predict_visible_uncertainty_floor = 0.03

    x = 300000.0 + (np.arange(20, dtype=float) + 0.5) * _RES
    y = 8200000.0 - (np.arange(20, dtype=float) + 0.5) * _RES

    def _grid(value: float) -> xr.DataArray:
        return xr.DataArray(
            np.full((20, 20), value, dtype=np.float32),
            dims=("y", "x"),
            coords={"y": y, "x": x},
        ).rio.write_crs(_CRS)

    observation = SimpleNamespace(
        bounds=(300000.0, 8198800.0, 301200.0, 8200000.0),
        crs=_CRS,
        metadata={"observation_time": datetime(2026, 3, 15, 9, 30)},
        sensor_config=SENTINEL2A_CONFIG,
        toa=xr.Dataset({"B8A": _grid(0.26), "B11": _grid(0.31), "B12": _grid(0.29)}),
        geometry=SimpleNamespace(sza=_grid(0.4), vza=_grid(0.1), raa=_grid(2.0)),
    )
    atmo = SimpleNamespace(
        aot=_grid(0.2),
        tcwv=_grid(2.0),
        tco3=_grid(0.3),
        aot_unc=_grid(0.1),
        tcwv_unc=_grid(0.3),
        tco3_unc=_grid(0.03),
        elevation=_grid(0.0),
    )

    class _FakeRT:
        def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
            _ = (atmo_state, band, compute_jacobian)
            from siac.runtime import RTCoefficients

            xap = xr.full_like(geometry.sza, 0.8)
            return RTCoefficients(xap=xap, xbp=xr.zeros_like(xap), xcp=xr.zeros_like(xap))

    source = bestpixel_source_bands(bands, endpoint="pc")
    target = [SENTINEL2A_CONFIG.get_band(n) for n in ("B02", "B04")]

    prior = build_bestpixel_surface_prior(
        config,
        observation,
        _RES,
        bands=bands,
        source_bands=source,
        target_bands=target,
        spectral_library=None,
        k_neighbors=5,
        maiac_day_aod=lambda *_a, **_k: {},
        atmo_prior=atmo,
        rt_model=_FakeRT(),
    )

    assert list(prior.boa["band"].values) == ["B02", "B04"]
    np.testing.assert_allclose(
        np.nanmedian(prior.boa_unc.sel(band="B02").values),
        0.03,
        atol=1e-6,
    )
    np.testing.assert_allclose(
        np.nanmedian(prior.boa_unc.sel(band="B04").values),
        0.03,
        atol=1e-6,
    )


def test_build_surface_prior_median_across_realizations(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # coastal: 0.1, 0.2, 0.3 across the three realizations.
    periods = [
        _period_dict(_BANDS, band_dn={"coastal": 1000, "blue": 800, "red": 1200}),
        _period_dict(_BANDS, band_dn={"coastal": 2000, "blue": 800, "red": 1200}),
        _period_dict(_BANDS, band_dn={"coastal": 3000, "blue": 800, "red": 1200}),
    ]
    seen: dict[str, Any] = {}
    _install_fake_bestpixel(monkeypatch, periods, seen=seen)
    bands = _identity_bands()

    prior = build_bestpixel_surface_prior(
        _fake_config(lookback_years=3),
        _fake_observation(),
        _RES,
        bands=_BANDS,
        source_bands=bands,
        target_bands=bands,
        spectral_library=None,
        k_neighbors=5,
        maiac_day_aod=lambda *_a, **_k: {},  # empty gate -> keep all days
    )

    assert isinstance(prior, SurfacePrior)
    assert prior.kernels is None
    assert prior.boa.dims == ("band", "y", "x")
    assert prior.boa.shape == (3, 10, 10)
    assert prior.boa_unc.shape == (3, 10, 10)
    assert prior.mask.dims == ("y", "x")
    assert prior.mask.shape == (10, 10)
    assert bool(prior.mask.values.all())

    coastal = prior.boa.sel(band="coastal").values
    np.testing.assert_allclose(np.nanmedian(coastal), 0.2, atol=1e-4)
    coastal_unc = prior.boa_unc.sel(band="coastal").values
    # RMSE about the median (matches comp_ref): sqrt(mean([0.01, 0, 0.01])) = sqrt(0.02/3).
    np.testing.assert_allclose(np.nanmedian(coastal_unc), np.sqrt(0.02 / 3.0), atol=1e-3)

    # Gate was empty -> bestpixel called WITHOUT an aod_by_day gate.
    assert "aod_by_day" not in seen["kwargs"]
    assert seen["kwargs"]["emit_uncertainty"] is True
    assert seen["kwargs"]["output_crs"] == "utm"
    assert seen["kwargs"]["years"] == [2023, 2024, 2025]
    assert seen["kwargs"]["months"] == [3]


def test_build_surface_prior_daily_median_reduces_clean_days_first(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen: dict[str, Any] = {}
    _install_fake_daily_bestpixel(monkeypatch, seen=seen)
    config = _fake_config(lookback_years=1)
    config.algorithms.surface_prior.bestpixel_window_reduction = "daily_median"
    bands = _identity_bands()

    prior = build_bestpixel_surface_prior(
        config,
        _fake_observation(),
        _RES,
        bands=_BANDS,
        source_bands=bands,
        target_bands=bands,
        spectral_library=None,
        k_neighbors=5,
        maiac_day_aod=lambda *_a, **_k: {
            "2025-03-01": 0.1,
            "2025-03-05": 0.2,
        },
    )

    calls = seen["calls"]
    assert [call["datetime"] for call in calls] == [
        "2025-03-01/2025-03-02",
        "2025-03-05/2025-03-06",
    ]
    assert all(call["top_k"] == 1 for call in calls)
    coastal = prior.boa.sel(band="coastal").values
    np.testing.assert_allclose(np.nanmedian(coastal), 0.2, atol=1e-5)


def test_build_surface_prior_resolve_on_prior_grid_adopts_native_grid(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # The composite native grid is OFFSET from the obs-bounds grid (the real
    # ~20 m tile-grid shift). With surface_driven_resolve_on_prior_grid the prior
    # must carry the NATIVE grid as ``solver_grid`` (so M4 resolves on it and the
    # prior is never smeared onto obs bounds), not the obs-bounds grid.
    xmin_native, ymax_native = 300020.0, 8200020.0  # +20 m vs the obs-bounds origin
    period = _period_dict(_BANDS, band_dn={"coastal": 1000, "blue": 800, "red": 1200})
    period["grid"]["transform"] = [_RES, 0.0, xmin_native, 0.0, -_RES, ymax_native]
    _install_fake_bestpixel(monkeypatch, [period])
    bands = _identity_bands()

    cfg = _fake_config(lookback_years=1)
    cfg.algorithms.solver = SimpleNamespace(surface_driven_resolve_on_prior_grid=True)

    prior = build_bestpixel_surface_prior(
        cfg,
        _fake_observation(),
        _RES,
        bands=_BANDS,
        source_bands=bands,
        target_bands=bands,
        spectral_library=None,
        k_neighbors=5,
        maiac_day_aod=lambda *_a, **_k: {},
    )

    assert prior.solver_grid is not None
    assert prior.solver_grid.rio.crs is not None
    # Pixel-center origin reflects the +20 m native offset, NOT the obs-bounds grid.
    np.testing.assert_allclose(
        float(prior.solver_grid["x"].values[0]), xmin_native + 0.5 * _RES, atol=1e-6
    )
    np.testing.assert_allclose(
        float(prior.solver_grid["y"].values[0]), ymax_native - 0.5 * _RES, atol=1e-6
    )
    # The prior reflectance lives on that same native grid.
    np.testing.assert_allclose(prior.boa["x"].values, prior.solver_grid["x"].values, atol=1e-6)


def test_build_surface_prior_default_solver_grid_is_none(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Default config carries no algorithms.solver flag -> legacy obs-bounds grid
    # and no solver_grid (M4 keeps building the grid from the observation bounds).
    period = _period_dict(_BANDS, band_dn={"coastal": 1000, "blue": 800, "red": 1200})
    _install_fake_bestpixel(monkeypatch, [period])
    bands = _identity_bands()

    prior = build_bestpixel_surface_prior(
        _fake_config(lookback_years=1),
        _fake_observation(),
        _RES,
        bands=_BANDS,
        source_bands=bands,
        target_bands=bands,
        spectral_library=None,
        k_neighbors=5,
        maiac_day_aod=lambda *_a, **_k: {},
    )
    assert prior.solver_grid is None


def test_build_surface_prior_single_realization_uses_composite_unc(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    periods = [
        _period_dict(_BANDS, band_dn={"coastal": 1500, "blue": 800, "red": 1200}, unc_dn=500),
    ]
    _install_fake_bestpixel(monkeypatch, periods)
    bands = _identity_bands()

    prior = build_bestpixel_surface_prior(
        _fake_config(lookback_years=1),
        _fake_observation(),
        _RES,
        bands=_BANDS,
        source_bands=bands,
        target_bands=bands,
        spectral_library=None,
        k_neighbors=5,
        maiac_day_aod=lambda *_a, **_k: {},
    )
    # Single realization -> boa_unc = per-composite bestpixel unc (DN 500 -> 0.05).
    coastal_unc = prior.boa_unc.sel(band="coastal").values
    np.testing.assert_allclose(np.nanmedian(coastal_unc), 0.05, atol=1e-3)


def test_build_surface_prior_passes_gate_to_bestpixel(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    periods = [_period_dict(_BANDS, band_dn={"coastal": 1000, "blue": 800, "red": 1200})]
    seen: dict[str, Any] = {}
    _install_fake_bestpixel(monkeypatch, periods, seen=seen)
    bands = _identity_bands()

    def _gate(*_a: Any, **_k: Any) -> dict[str, float]:
        return {"2025-03-01": 0.1, "2025-03-08": 0.4}

    build_bestpixel_surface_prior(
        _fake_config(lookback_years=1),
        _fake_observation(),
        _RES,
        bands=_BANDS,
        source_bands=bands,
        target_bands=bands,
        spectral_library=None,
        k_neighbors=5,
        maiac_day_aod=_gate,
    )
    assert "aod_by_day" in seen["kwargs"]
    assert seen["kwargs"]["reject_unknown"] is True
    # frac default 0.6 -> ceil(0.6*2)=2 kept, both days survive.
    assert set(seen["kwargs"]["aod_by_day"]) == {"2025-03-01", "2025-03-08"}


def test_build_surface_prior_gate_failure_is_tolerated(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    periods = [_period_dict(_BANDS, band_dn={"coastal": 1000, "blue": 800, "red": 1200})]
    seen: dict[str, Any] = {}
    _install_fake_bestpixel(monkeypatch, periods, seen=seen)
    bands = _identity_bands()

    def _broken_gate(*_a: Any, **_k: Any) -> dict[str, float]:
        raise RuntimeError("earthaccess down")

    prior = build_bestpixel_surface_prior(
        _fake_config(lookback_years=1),
        _fake_observation(),
        _RES,
        bands=_BANDS,
        source_bands=bands,
        target_bands=bands,
        spectral_library=None,
        k_neighbors=5,
        maiac_day_aod=_broken_gate,
    )
    assert isinstance(prior, SurfacePrior)
    # Gate failed -> no aod_by_day, all days kept.
    assert "aod_by_day" not in seen["kwargs"]


def test_build_surface_prior_raises_when_all_periods_nodata(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    period = _period_dict(_BANDS, band_dn={"coastal": 1000, "blue": 800, "red": 1200})
    period["quality"] = np.full((10, 10), 65535, dtype=np.uint16)  # all nodata
    _install_fake_bestpixel(monkeypatch, [period])
    bands = _identity_bands()
    with pytest.raises(RuntimeError, match="no usable composites"):
        build_bestpixel_surface_prior(
            _fake_config(lookback_years=1),
            _fake_observation(),
            _RES,
            bands=_BANDS,
            source_bands=bands,
            target_bands=bands,
            spectral_library=None,
            k_neighbors=5,
            maiac_day_aod=lambda *_a, **_k: {},
        )


def test_build_surface_prior_missing_band_raises(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    period = _period_dict(("coastal", "blue"), band_dn={"coastal": 1000, "blue": 800})
    _install_fake_bestpixel(monkeypatch, [period])
    bands = _identity_bands()
    with pytest.raises(KeyError, match="missing requested band"):
        build_bestpixel_surface_prior(
            _fake_config(lookback_years=1),
            _fake_observation(),
            _RES,
            bands=_BANDS,  # asks for red, which the period omits
            source_bands=bands,
            target_bands=bands,
            spectral_library=None,
            k_neighbors=5,
            maiac_day_aod=lambda *_a, **_k: {},
        )


def test_build_surface_prior_requires_observation_time(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_fake_bestpixel(monkeypatch, [])
    bands = _identity_bands()
    observation = SimpleNamespace(bounds=_BOUNDS, crs=_CRS, metadata={})
    with pytest.raises(ValueError, match="observation_time"):
        build_bestpixel_surface_prior(
            _fake_config(),
            observation,
            _RES,
            bands=_BANDS,
            source_bands=bands,
            target_bands=bands,
            spectral_library=None,
            k_neighbors=5,
            maiac_day_aod=lambda *_a, **_k: {},
        )


# --------------------------------------------------------------------------- #
# year / month / default-gate resolution helpers
# --------------------------------------------------------------------------- #
def test_resolve_years_lookback_window() -> None:
    assert _resolve_years(datetime(2026, 3, 15), 3) == (2023, 2024, 2025)


def test_resolve_years_rejects_non_positive_lookback() -> None:
    with pytest.raises(ValueError, match="lookback_years"):
        _resolve_years(datetime(2026, 3, 15), 0)


def test_resolve_months_scene_month_or_explicit() -> None:
    assert _resolve_months(datetime(2026, 3, 15), None) == (3,)
    assert _resolve_months(datetime(2026, 3, 15), [6, 7, 8]) == (6, 7, 8)


def test_resolve_months_seasonal_window() -> None:
    # +/-1 around the scene month (matches the comp_ref harness 3-month window).
    assert _resolve_months(datetime(2026, 5, 15), None, 1) == (4, 5, 6)
    # Wraps across the year boundary (Jan +/-1 -> Dec, Jan, Feb).
    assert _resolve_months(datetime(2026, 1, 15), None, 1) == (1, 2, 12)
    # Window 0 keeps the scene month only; explicit months override the window.
    assert _resolve_months(datetime(2026, 5, 15), None, 0) == (5,)
    assert _resolve_months(datetime(2026, 5, 15), [7], 1) == (7,)


def test_default_maiac_day_aod_is_callable() -> None:
    # Constructs the earthaccess-backed provider lazily (no network at build).
    gate = _default_maiac_day_aod()
    assert callable(gate)


# --------------------------------------------------------------------------- #
# resolver dispatch + l1c guard
# --------------------------------------------------------------------------- #
def test_resolve_surface_prior_provider_handles_bestpixel() -> None:
    config = SIACConfig.model_validate({"algorithms": {"surface_prior": {"method": "bestpixel"}}})
    fn = resolve_surface_prior_provider(config)
    assert callable(fn)
    assert fn.requires_atmo_prior is False


def test_bestpixel_l1c_source_not_implemented() -> None:
    config = SIACConfig.model_validate(
        {"algorithms": {"surface_prior": {"method": "bestpixel", "bestpixel_source": "l1c"}}}
    )
    fn = make_bestpixel_surface_prior_fn(config)
    observation = SimpleNamespace(
        sensor_config=None,
        bounds=_BOUNDS,
        crs=_CRS,
        metadata={"observation_time": datetime(2026, 3, 15)},
    )
    with pytest.raises(NotImplementedError, match="l1c"):
        fn(observation, None, None, _RES)
