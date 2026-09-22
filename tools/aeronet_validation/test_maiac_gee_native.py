from __future__ import annotations

import datetime as dt

import numpy as np
import pytest
from tools.aeronet_validation import maiac_cmr_replay, maiac_gee_native
from tools.aeronet_validation.maiac_cmr_replay import GranuleId
from tools.aeronet_validation.maiac_gee_native import NativeGranule, TeacherMaiacPrior

from siac.adapters.atmo.mcd19_earthaccess import NoAtmosphericDataError

CRS = "EPSG:32652"
OBS = dt.datetime(2024, 1, 27, 2, 19, 41)
GRANULE = GranuleId(day="A2024027", tile="h28v05", production="2024030172711")
#: Best-quality AOD_QA: clear (001), normal adjacency (000), QA level 0.
BEST_QA = 1
#: Clear and normal adjacency, but QA level 1: passes only the loose ``qa > 0`` mask.
LOOSE_QA = 1 | (1 << 8)


def _aoi() -> tuple[float, float, float, float]:
    """A 128 x 60 m teacher AOI around Osan, inside MODLAND tile h28v05."""
    from pyproj import Transformer

    x, y = Transformer.from_crs("EPSG:4326", CRS, always_xy=True).transform(127.03, 37.09)
    x0, y0 = 60.0 * np.floor(x / 60.0), 60.0 * np.floor(y / 60.0)
    return (x0, y0, x0 + 128 * 60.0, y0 + 128 * 60.0)


def _native(aod: np.ndarray, qa: np.ndarray, times: tuple[dt.datetime, ...]) -> NativeGranule:
    x, y, _, _, _, window = maiac_gee_native._native_window(GRANULE.tile, _aoi(), CRS, 2)
    return NativeGranule(
        granule=GRANULE,
        times=times,
        values={
            "Optical_Depth_055": aod.astype(np.int32),
            "AOD_Uncertainty": np.full(aod.shape, 500, dtype=np.int32),
            "AOD_QA": qa.astype(np.int32),
        },
        window=window,
        x=x,
        y=y,
    )


def _block_shape() -> tuple[int, int]:
    r0, r1, c0, c1 = maiac_gee_native._native_window(GRANULE.tile, _aoi(), CRS, 2)[-1]
    return r1 - r0 + 1, c1 - c0 + 1


@pytest.fixture
def offline(monkeypatch: pytest.MonkeyPatch) -> dict[str, object]:
    """Replace every network touch point; the test sets what each returns."""
    state: dict[str, object] = {"granules": [GRANULE], "native": None, "whole": {}}

    def teacher_granules(*args: object, **kwargs: object) -> list[GranuleId]:
        del args, kwargs
        return list(state["granules"])  # type: ignore[call-overload]

    def native_granules(granules: list[GranuleId], *args: object, **kwargs: object) -> dict:
        del args, kwargs
        return dict.fromkeys(granules, state["native"])

    def whole_granule_has(granule: GranuleId, loose: bool) -> bool:
        del granule
        answers = state["whole"]
        assert isinstance(answers, dict)
        if loose not in answers:
            raise AssertionError(f"unexpected whole-granule query (loose={loose})")
        return bool(answers[loose])

    monkeypatch.setattr(maiac_gee_native, "_initialize", lambda: None)
    monkeypatch.setattr(maiac_cmr_replay, "teacher_granules", teacher_granules)
    monkeypatch.setattr(maiac_gee_native, "native_granules", native_granules)
    monkeypatch.setattr(maiac_gee_native, "_whole_granule_has", whole_granule_has)
    return state


def test_granule_id_round_trips_the_cmr_name() -> None:
    granule = GranuleId.from_granule_ur("MCD19A2.A2024027.h28v05.061.2024030172711")

    assert granule == GRANULE
    assert granule.date == dt.date(2024, 1, 27)
    assert granule.filename == "MCD19A2.A2024027.h28v05.061.2024030172711.hdf"
    with pytest.raises(ValueError, match="unexpected MCD19A2"):
        GranuleId.from_granule_ur("MOD04_L2.A2024027.0215.061.2024027150000")


def test_teacher_prior_takes_the_nearest_valid_orbit_per_pixel(offline) -> None:
    rows, cols = _block_shape()
    aod = np.stack([np.full((rows, cols), 100), np.full((rows, cols), 300)])
    qa = np.full((2, rows, cols), BEST_QA)
    # The nearer orbit (03:00, 41 min away) is missing on the upper half of the
    # block; there the earlier orbit (01:30, 49 min away) must supply the value.
    qa[1, : rows // 2] = 0
    offline["native"] = _native(
        aod, qa, (dt.datetime(2024, 1, 27, 1, 30), dt.datetime(2024, 1, 27, 3, 0))
    )

    prior = maiac_gee_native.lineage_teacher_prior(_aoi(), CRS, OBS)

    values = np.asarray(prior.aot)
    assert values.shape == (128, 128)
    assert values.dtype == np.float32
    assert set(np.unique(values)) == {np.float32(0.1), np.float32(0.3)}
    assert np.all(values[0] == np.float32(0.1)) and np.all(values[-1] == np.float32(0.3))
    # The lineage scales uncertainty in float32: 500 * 0.0001 lands one ulp below 0.05.
    np.testing.assert_allclose(np.asarray(prior.aot_unc), 0.05, rtol=1e-6)
    assert str(prior.aot.rio.crs) == CRS
    assert prior.granules == [GRANULE.filename]
    assert prior.fallback_granules == []


def test_loose_fallback_only_when_the_whole_granule_lacks_best_pixels(offline) -> None:
    rows, cols = _block_shape()
    aod = np.full((1, rows, cols), 250)
    qa = np.full((1, rows, cols), LOOSE_QA)
    offline["native"] = _native(aod, qa, (dt.datetime(2024, 1, 27, 2, 20),))

    # Best-quality pixels elsewhere in the tile: the committed provider kept
    # the strict mask, so this AOI has no MAIAC at all.
    offline["whole"] = {False: True}
    with pytest.raises(NoAtmosphericDataError, match="no QA-valid AOD"):
        maiac_gee_native.lineage_teacher_prior(_aoi(), CRS, OBS)

    # No best-quality pixel anywhere in the tile: it degraded to ``qa > 0``.
    offline["whole"] = {False: False}
    prior = maiac_gee_native.lineage_teacher_prior(_aoi(), CRS, OBS)
    np.testing.assert_array_equal(np.asarray(prior.aot), np.float32(0.25))
    assert prior.fallback_granules == [GRANULE.filename]


def test_absent_maiac_is_signalled_as_missing_data(offline) -> None:
    offline["granules"] = []
    with pytest.raises(NoAtmosphericDataError, match="no granule reaching"):
        maiac_gee_native.lineage_teacher_prior(_aoi(), CRS, OBS)


def test_atmospheric_state_carries_the_provider_defaults() -> None:
    import xarray as xr

    aot = xr.DataArray(
        np.full((2, 3), 0.2, dtype=np.float32),
        dims=("y", "x"),
        coords={"y": [1.0, 0.0], "x": [0.0, 1.0, 2.0]},
    )
    state = TeacherMaiacPrior(
        aot=aot, aot_unc=aot * 0.5, granules=[], fallback_granules=[]
    ).atmospheric_state()

    np.testing.assert_array_equal(state.aot.values, aot.values)
    np.testing.assert_array_equal(state.aot_unc.values, np.float32(0.1))
    for field, expected in (
        ("tcwv", 1.5),
        ("tco3", 0.30),
        ("tcwv_unc", 0.3),
        ("tco3_unc", 0.03),
        ("elevation", 0.0),
    ):
        value = getattr(state, field)
        assert value.dtype == np.float32
        np.testing.assert_array_equal(value.values, np.float32(expected))
        np.testing.assert_array_equal(value.x.values, aot.x.values)


def test_catalogue_replay_refuses_to_log_in() -> None:
    class Source:
        is_authenticated = False

        def search_granules(self, **kwargs: object) -> list[dict]:
            del kwargs
            self._ensure_auth()
            return []

    class Provider:
        source = Source()
        provider = "LPCLOUD"

    with pytest.raises(RuntimeError, match="never log in"):
        maiac_cmr_replay._search(Provider(), "MCD19A2", (0.0, 0.0, 1.0, 1.0), CRS, None, 8)
