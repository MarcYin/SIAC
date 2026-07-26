"""Live L1C surface library: index reading, mosaicking, correction, RT space.

Every test injects fakes. The real library reads L1C imagery from the public
Google Cloud Sentinel-2 bucket, Sentinel-2 L2A water vapour from the Planetary
Computer, ozone from CAMS and terrain from a DEM, then corrects through a
compiled 6S backend. Here fake readers supply flat planes and a fake RT backend
returns ``xap = 1 + aot`` (so a correct pipeline yields ``(1 + aod) * toa``
exactly), which pins the per-pixel AOD selection, the mosaic assembly, the
measured-state reduction and the grid handling without touching the network.

The second theme of these tests is what happens when a measured input is
*missing*: every such case must raise and name the input, never substitute a
plausible constant.
"""

from __future__ import annotations

import json
from datetime import datetime
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any

import numpy as np
import pytest
import xarray as xr

from siac.adapters.live_l1c_library import (
    LIBRARY_BAND_NAMES,
    CAMSTCO3Source,
    IndexImage,
    LiveL1CSurfaceLibrary,
    MissingLibraryInputError,
    _baseline_offset,
    _keep_realizations,
    correct_toa_realization,
    read_mosaic_index,
    resolve_mgrs_tile,
)
from siac.adapters.surface_library import PREDICTOR_BAND_ORDER
from siac.catalog.sensors.sentinel2 import SENTINEL2A_CONFIG
from siac.config.algorithms import RTSetupConfig
from siac.domain.rt_space import RTSpace

if TYPE_CHECKING:
    from collections.abc import Sequence
    from pathlib import Path

# 10x10 grid at 60 m in UTM 33N (matches _build_target_template exactly).
_CRS = "EPSG:32633"
_RES = 60.0
_BOUNDS = (300000.0, 8199400.0, 300600.0, 8200000.0)
_TRANSFORM = (60.0, 0.0, 300000.0, 0.0, -60.0, 8200000.0)
_DEM = "/vsicurl/https://example.invalid/dem.vrt"
_CAMS = "/gws/example/cams"
_ELEVATION_KM = 0.42
_TCO3 = 0.31

_IMAGE_A = {
    "idx": "20220605T101031_20220605T101025_T33WXP",
    "day": "2022-06-05",
    "ratio": 1.0,
    "sza": 40.0,
    "saa": 150.0,
    "vza": 3.0,
    "vaa": 190.0,
}
_IMAGE_B = {
    "idx": "20220615T101031_20220615T101025_T33WXP",
    "day": "2022-06-15",
    "ratio": 1.0,
    "sza": 38.0,
    "saa": 152.0,
    "vza": 3.5,
    "vaa": 191.0,
}
_IMAGE_C = {
    "idx": "20220710T101031_20220710T101025_T33WXP",
    "day": "2022-07-10",
    "ratio": 1.0,
    "sza": 35.0,
    "saa": 155.0,
    "vza": 2.0,
    "vaa": 189.0,
}
_DAY_AOD = {"2022-06-05": 0.2, "2022-06-15": 0.4, "2022-07-10": 0.3}


# --------------------------------------------------------------------------- #
# fakes
# --------------------------------------------------------------------------- #
class _FakeCoefficients:
    def __init__(self, xap: float) -> None:
        self.xap = np.array([[xap]])
        self.xbp = np.array([[0.0]])
        self.xcp = np.array([[0.0]])


class _FakeRTModel:
    """Backend whose ``xap`` is ``1 + aot``: correction is ``(1 + aod) * toa``."""

    backend_name = "sixs"

    def __init__(
        self,
        rt_setup: Any | None = None,
        observation_time: datetime | None = None,
        clones: list[_FakeRTModel] | None = None,
    ) -> None:
        self.rt_setup = rt_setup
        self.observation_time = observation_time
        self.clones: list[_FakeRTModel] = [] if clones is None else clones
        self.aot_nodes: list[float] = []
        self.states: list[tuple[float, float, float]] = []

    def compute_coefficients_multi(
        self, geometry: Any, atmo_state: Any, bands: Sequence[Any]
    ) -> list[_FakeCoefficients]:
        _ = geometry

        def _scalar(field: Any) -> float:
            return float(np.ravel(np.asarray(field.values))[0])

        aot = _scalar(atmo_state.aot)
        self.aot_nodes.append(aot)
        self.states.append(
            (_scalar(atmo_state.tcwv), _scalar(atmo_state.tco3), _scalar(atmo_state.elevation))
        )
        return [_FakeCoefficients(1.0 + aot) for _ in bands]

    def with_rt_setup(self, rt_setup: Any) -> _FakeRTModel:
        clone = _FakeRTModel(rt_setup, self.observation_time, self.clones)
        self.clones.append(clone)
        return clone


class _FakeReader:
    """Returns a flat TOA plane per sensing token; unknown tokens read as None."""

    def __init__(self, planes: dict[str, float]) -> None:
        self._planes = planes
        self.calls: list[dict[str, Any]] = []

    def read(
        self,
        *,
        mgrs_tile: str,
        sensing_token: str,
        crs: str,
        transform: tuple[float, float, float, float, float, float],
        width: int,
        height: int,
    ) -> np.ndarray | None:
        self.calls.append(
            {
                "mgrs_tile": mgrs_tile,
                "sensing_token": sensing_token,
                "crs": crs,
                "transform": transform,
                "width": width,
                "height": height,
            }
        )
        value = self._planes.get(sensing_token)
        if value is None:
            return None
        return np.full((len(LIBRARY_BAND_NAMES), height, width), value, dtype=np.float32)


class _FakeWVPReader:
    """Returns a flat L2A water-vapour plane (cm) per sensing token."""

    def __init__(self, values: dict[str, float]) -> None:
        self._values = values
        self.calls: list[dict[str, Any]] = []

    def read(
        self,
        *,
        mgrs_tile: str,
        sensing_token: str,
        day: str,
        crs: str,
        transform: tuple[float, float, float, float, float, float],
        width: int,
        height: int,
    ) -> np.ndarray | None:
        _ = (crs, transform)
        self.calls.append({"mgrs_tile": mgrs_tile, "sensing_token": sensing_token, "day": day})
        value = self._values.get(sensing_token)
        if value is None:
            return None
        return np.full((height, width), value, dtype=np.float32)


class _FakeTCO3:
    def __init__(self, values: dict[str, float] | float = _TCO3) -> None:
        self._values = values
        self.days: list[str] = []

    def __call__(self, day: str, *, lon: float, lat: float) -> float:
        _ = (lon, lat)
        self.days.append(day)
        if isinstance(self._values, dict):
            if day not in self._values:
                raise MissingLibraryInputError(f"No CAMS ozone for {day}")
            return self._values[day]
        return float(self._values)


@pytest.fixture(autouse=True)
def _fake_dem(monkeypatch: pytest.MonkeyPatch) -> None:
    """Sample a constant terrain height instead of reading the remote DEM."""
    import siac.geo.dem as dem_module

    def _read(template: xr.DataArray, dem_path: str | None) -> xr.DataArray:
        _ = dem_path
        return xr.full_like(template, _ELEVATION_KM, dtype=np.float32)

    monkeypatch.setattr(dem_module, "read_elevation_km", _read)


def _observation(scene_key: str = "SITE__T33WXP_20220801T101031") -> Any:
    return SimpleNamespace(
        sensor_config=SENTINEL2A_CONFIG,
        bounds=_BOUNDS,
        crs=_CRS,
        metadata={"scene_key": scene_key, "observation_time": datetime(2022, 8, 1, 10, 10, 31)},
    )


def _write_index(
    root: Path,
    *,
    scene_key: str = "SITE__T33WXP_20220801T101031",
    winners: np.ndarray | None = None,
    day_aod: dict[str, float] | None = None,
) -> Path:
    """Two months: June has two winning acquisitions, July has one."""
    if winners is None:
        winners = np.zeros((2, 10, 10), dtype=np.int16)
        winners[0, 5:, :] = 1  # bottom half of June won by the second acquisition
        winners[1, 0, 0] = -1  # one July pixel never cleared
    scalars = _DAY_AOD if day_aod is None else day_aod
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"{scene_key}.npz"
    np.savez(
        path,
        months=np.array(["2022-06", "2022-07"]),
        winners=winners,
        image_table=np.array([json.dumps(im) for im in (_IMAGE_A, _IMAGE_B, _IMAGE_C)]),
        day_scalars=np.array([json.dumps({"day": d, "aod": a}) for d, a in scalars.items()]),
    )
    return path


_ALL_WVP = {"20220605T101031": 1.2, "20220615T101031": 1.6, "20220710T101031": 2.4}


def _library(
    root: Path,
    reader: _FakeReader,
    rt_model: _FakeRTModel,
    *,
    wvp: dict[str, float] | _FakeWVPReader | None = None,
    tco3: Any | None = None,
    **kwargs: Any,
) -> Any:
    wvp_reader = wvp if isinstance(wvp, _FakeWVPReader) else _FakeWVPReader(wvp or _ALL_WVP)
    return LiveL1CSurfaceLibrary(
        root,
        rt_model=rt_model,
        dem_path=_DEM,
        cams_data_path=_CAMS,
        toa_reader=reader,
        wvp_reader=wvp_reader,
        tco3_source=tco3 if tco3 is not None else _FakeTCO3(),
        **kwargs,
    )


# --------------------------------------------------------------------------- #
# index reading
# --------------------------------------------------------------------------- #
def test_index_splits_the_flat_image_table_into_month_blocks(tmp_path: Path) -> None:
    path = _write_index(tmp_path)

    index = read_mosaic_index(path)

    assert index.months == ("2022-06", "2022-07")
    assert [image.day for image in index.images_by_month["2022-06"]] == [
        "2022-06-05",
        "2022-06-15",
    ]
    assert [image.day for image in index.images_by_month["2022-07"]] == ["2022-07-10"]
    assert index.day_aod == _DAY_AOD
    assert index.winners.shape == (2, 10, 10)


def test_index_month_count_mismatch_is_reported(tmp_path: Path) -> None:
    path = _write_index(tmp_path, winners=np.zeros((3, 10, 10), dtype=np.int16))

    with pytest.raises(ValueError, match="2 months"):
        read_mosaic_index(path)


def test_sensing_token_is_the_first_field_of_the_image_id() -> None:
    image = IndexImage(
        image_id="20220605T101031_20220605T101025_T33WXP",
        day="2022-06-05",
        sza=1.0,
        saa=2.0,
        vza=3.0,
        vaa=4.0,
    )
    assert image.sensing_token == "20220605T101031"


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def test_mgrs_tile_comes_from_the_scene_key_or_metadata() -> None:
    assert resolve_mgrs_tile(_observation()) == "33WXP"
    bare = SimpleNamespace(metadata={"mgrs_tile": "T31UDQ"})
    assert resolve_mgrs_tile(bare) == "31UDQ"
    tile_id = SimpleNamespace(metadata={"tile_id": "S2A_OPER_MSI_L1C_TL_x_x_T52SCG_N05.00"})
    assert resolve_mgrs_tile(tile_id) == "52SCG"
    # An explicit scene key wins over metadata.
    assert resolve_mgrs_tile(_observation(), "OTHER__T10SEG_20200101T000000") == "10SEG"


def test_mgrs_tile_missing_is_reported() -> None:
    with pytest.raises(ValueError, match="MGRS tile"):
        resolve_mgrs_tile(SimpleNamespace(metadata={"observation_time": datetime(2022, 8, 1)}))


def test_baseline_offset_follows_the_radiometric_shift() -> None:
    assert _baseline_offset("N0500") == 1000.0
    assert _baseline_offset("N0400") == 1000.0
    assert _baseline_offset("N0209") == 0.0
    assert _baseline_offset(None) == 0.0
    assert _baseline_offset("N0000") == 0.0


def test_keep_filter_yields_when_too_few_realizations_survive() -> None:
    dense = np.ones((2, 3, 3), dtype=np.float32)
    sparse = np.full((2, 3, 3), np.nan, dtype=np.float32)

    # Five realizations, four dense: the sparse one is dropped.
    assert _keep_realizations([dense, dense, sparse, dense, dense], 0.6) == [0, 1, 3, 4]
    # Dropping would leave fewer than four: keep everything rather than destroy
    # the temporal spread that becomes the prior uncertainty.
    assert _keep_realizations([dense, dense, sparse, dense, sparse], 0.6) == [0, 1, 2, 3, 4]
    assert _keep_realizations([dense, sparse], 0.6) == [0, 1]


# --------------------------------------------------------------------------- #
# correction
# --------------------------------------------------------------------------- #
def test_correction_interpolates_the_rt_curves_per_pixel() -> None:
    rt_model = _FakeRTModel()
    bands = [SENTINEL2A_CONFIG.get_band(name) for name in LIBRARY_BAND_NAMES]
    toa = np.full((len(bands), 2, 2), 0.1, dtype=np.float32)
    aod = np.array([[0.2, 0.4], [0.2, 0.4]], dtype=np.float32)

    corrected = correct_toa_realization(
        rt_model=rt_model,
        sensor_bands=bands,
        toa=toa,
        aod=aod,
        sza_deg=40.0,
        saa_deg=150.0,
        vza_deg=3.0,
        vaa_deg=190.0,
        tcwv=1.4,
        tco3=_TCO3,
        elevation_km=_ELEVATION_KM,
    )

    # xap = 1 + aot, xbp = xcp = 0 -> corrected = (1 + aod) * toa, exactly, because
    # linear interpolation of a linear curve is exact between the AOD nodes.
    np.testing.assert_allclose(corrected[:, :, 0], 1.2 * 0.1, rtol=1e-5)
    np.testing.assert_allclose(corrected[:, :, 1], 1.4 * 0.1, rtol=1e-5)
    # One 6S invocation per AOD node, spanning the mosaic's own AOD range, all at
    # the measured water vapour / ozone / terrain state.
    assert len(rt_model.aot_nodes) == 8
    assert min(rt_model.aot_nodes) <= 0.2 and max(rt_model.aot_nodes) >= 0.4
    np.testing.assert_allclose(
        np.array(rt_model.states), [[1.4, _TCO3, _ELEVATION_KM]] * 8, rtol=1e-6
    )


def test_correction_rejects_a_band_count_mismatch() -> None:
    bands = [SENTINEL2A_CONFIG.get_band("B02")]
    with pytest.raises(ValueError, match="band, y, x"):
        correct_toa_realization(
            rt_model=_FakeRTModel(),
            sensor_bands=bands,
            toa=np.zeros((3, 2, 2), dtype=np.float32),
            aod=np.full((2, 2), 0.1, dtype=np.float32),
            sza_deg=40.0,
            saa_deg=150.0,
            vza_deg=3.0,
            vaa_deg=190.0,
            tcwv=1.4,
            tco3=_TCO3,
            elevation_km=0.0,
        )


def test_correction_without_any_finite_aod_is_rejected() -> None:
    bands = [SENTINEL2A_CONFIG.get_band(name) for name in LIBRARY_BAND_NAMES]
    with pytest.raises(MissingLibraryInputError, match="no aerosol state"):
        correct_toa_realization(
            rt_model=_FakeRTModel(),
            sensor_bands=bands,
            toa=np.zeros((len(bands), 2, 2), dtype=np.float32),
            aod=np.full((2, 2), np.nan, dtype=np.float32),
            sza_deg=40.0,
            saa_deg=150.0,
            vza_deg=3.0,
            vaa_deg=190.0,
            tcwv=1.4,
            tco3=_TCO3,
            elevation_km=0.0,
        )


# --------------------------------------------------------------------------- #
# realizations
# --------------------------------------------------------------------------- #
def test_realizations_mosaic_the_winning_days_and_correct_at_their_own_aod(
    tmp_path: Path,
) -> None:
    _write_index(tmp_path)
    reader = _FakeReader({"20220605T101031": 0.1, "20220615T101031": 0.2, "20220710T101031": 0.3})
    library = _library(tmp_path, reader, _FakeRTModel())

    realizations = library.realizations(_observation(), _RES, ["blue"])

    assert [r.label for r in realizations] == ["2022-06", "2022-07"]
    for realization in realizations:
        assert realization.band_names == PREDICTOR_BAND_ORDER
        assert realization.crs == _CRS
        assert realization.transform == _TRANSFORM
        assert realization.reflectance.shape == (len(LIBRARY_BAND_NAMES), 10, 10)

    june = realizations[0].reflectance
    # Top half won by the 2022-06-05 acquisition (toa 0.1 at AOD 0.2).
    np.testing.assert_allclose(june[:, :5, :], 1.2 * 0.1, rtol=1e-5)
    # Bottom half won by 2022-06-15 (toa 0.2 at AOD 0.4).
    np.testing.assert_allclose(june[:, 5:, :], 1.4 * 0.2, rtol=1e-5)

    july = realizations[1].reflectance
    np.testing.assert_allclose(july[:, 0, 1], 1.3 * 0.3, rtol=1e-5)
    # The pixel no acquisition cleared stays unfilled.
    assert np.all(np.isnan(july[:, 0, 0]))


def test_each_month_is_corrected_at_its_own_measured_state(tmp_path: Path) -> None:
    _write_index(tmp_path)
    reader = _FakeReader({"20220605T101031": 0.1, "20220615T101031": 0.2, "20220710T101031": 0.3})
    rt_model = _FakeRTModel()
    ozone = _FakeTCO3({"2022-06-05": 0.30, "2022-06-15": 0.32, "2022-07-10": 0.28})
    library = _library(tmp_path, reader, rt_model, tco3=ozone)

    library.realizations(_observation(), _RES, ["blue"])

    states = np.unique(np.round(np.array(rt_model.states), 6), axis=0)
    # June: half the pixels at WVP 1.2 and half at 1.6 -> median 1.4; ozone
    # likewise medians the two winning days' columns (0.30, 0.32) to 0.31.
    # July: the single winning acquisition's 2.4 cm and 0.28 atm-cm.
    np.testing.assert_allclose(
        states, [[1.4, 0.31, _ELEVATION_KM], [2.4, 0.28, _ELEVATION_KM]], atol=1e-6
    )
    assert sorted(set(ozone.days)) == ["2022-06-05", "2022-06-15", "2022-07-10"]


def test_only_the_winning_acquisitions_are_read_and_on_the_solver_grid(
    tmp_path: Path,
) -> None:
    # June's second acquisition wins nothing: it must never be fetched.
    winners = np.zeros((2, 10, 10), dtype=np.int16)
    _write_index(tmp_path, winners=winners)
    reader = _FakeReader({"20220605T101031": 0.1, "20220710T101031": 0.3})
    wvp = _FakeWVPReader(_ALL_WVP)
    library = _library(tmp_path, reader, _FakeRTModel(), wvp=wvp)

    library.realizations(_observation(), _RES, ["blue"])

    tokens = sorted(call["sensing_token"] for call in reader.calls)
    assert tokens == ["20220605T101031", "20220710T101031"]
    assert sorted(call["sensing_token"] for call in wvp.calls) == tokens
    for call in reader.calls:
        assert call["mgrs_tile"] == "33WXP"
        assert call["crs"] == _CRS
        assert call["transform"] == _TRANSFORM
        assert (call["width"], call["height"]) == (10, 10)


def test_an_unreadable_acquisition_degrades_its_month(tmp_path: Path) -> None:
    _write_index(tmp_path)
    # June's winners are both unavailable on GCS; July still builds.
    reader = _FakeReader({"20220710T101031": 0.3})
    library = _library(tmp_path, reader, _FakeRTModel())

    realizations = library.realizations(_observation(), _RES, ["blue"])

    assert [r.label for r in realizations] == ["2022-07"]


def test_no_readable_acquisition_is_an_error(tmp_path: Path) -> None:
    _write_index(tmp_path)
    library = _library(tmp_path, _FakeReader({}), _FakeRTModel())

    with pytest.raises(MissingLibraryInputError, match="no usable realizations"):
        library.realizations(_observation(), _RES, ["blue"])


def test_a_missing_scene_index_reports_the_store_and_key(tmp_path: Path) -> None:
    _write_index(tmp_path, scene_key="OTHER")
    library = _library(tmp_path, _FakeReader({}), _FakeRTModel())

    with pytest.raises(FileNotFoundError, match="mosaic index"):
        library.realizations(_observation(), _RES, ["blue"])


def test_days_the_index_never_scored_fall_back_to_the_maiac_gate(tmp_path: Path) -> None:
    _write_index(tmp_path, day_aod={})
    reader = _FakeReader({"20220605T101031": 0.1, "20220710T101031": 0.3})
    calls: list[Any] = []

    def _gate(bounds: Any, crs: str, periods: Any) -> dict[str, float]:
        calls.append((bounds, crs, list(periods)))
        return {"2022-06-05": 0.5, "2022-06-15": 0.5, "2022-07-10": 0.5}

    library = _library(tmp_path, reader, _FakeRTModel(), maiac_day_aod=_gate)
    realizations = library.realizations(_observation(), _RES, ["blue"])

    assert calls == [(_BOUNDS, _CRS, [(2022, 6), (2022, 7)])]
    np.testing.assert_allclose(realizations[0].reflectance[:, 0, 0], 1.5 * 0.1, rtol=1e-5)


# --------------------------------------------------------------------------- #
# missing measured inputs must raise, never default
# --------------------------------------------------------------------------- #
def test_a_failing_maiac_gate_is_an_error_not_a_default_aod(tmp_path: Path) -> None:
    _write_index(tmp_path, day_aod={})
    reader = _FakeReader({"20220605T101031": 0.1, "20220710T101031": 0.3})

    def _gate(bounds: Any, crs: str, periods: Any) -> dict[str, float]:
        raise RuntimeError("earthaccess down")

    library = _library(tmp_path, reader, _FakeRTModel(), maiac_day_aod=_gate)

    with pytest.raises(MissingLibraryInputError, match="measured AOD"):
        library.realizations(_observation(), _RES, ["blue"])


def test_an_empty_maiac_gate_is_an_error(tmp_path: Path) -> None:
    _write_index(tmp_path, day_aod={})
    reader = _FakeReader({"20220605T101031": 0.1})
    library = _library(tmp_path, reader, _FakeRTModel(), maiac_day_aod=lambda *_a, **_k: {})

    with pytest.raises(MissingLibraryInputError, match="assumed aerosol loading"):
        library.realizations(_observation(), _RES, ["blue"])


def test_a_winning_day_without_an_aod_scalar_is_an_error(tmp_path: Path) -> None:
    # The index scores June's days but not July's winner.
    _write_index(tmp_path, day_aod={"2022-06-05": 0.2, "2022-06-15": 0.4})
    reader = _FakeReader({"20220605T101031": 0.1, "20220615T101031": 0.2, "20220710T101031": 0.3})
    library = _library(tmp_path, reader, _FakeRTModel())

    with pytest.raises(MissingLibraryInputError, match="No measured AOD for 2022-07-10"):
        library.realizations(_observation(), _RES, ["blue"])


def test_one_acquisition_without_l2a_water_vapour_leaves_its_pixels_unscored(
    tmp_path: Path,
) -> None:
    # The 2022-06-15 acquisition has no L2A coverage; June's column is then the
    # measured 1.2 cm of the acquisition that does, not an invented default.
    _write_index(tmp_path)
    reader = _FakeReader({"20220605T101031": 0.1, "20220615T101031": 0.2, "20220710T101031": 0.3})
    rt_model = _FakeRTModel()
    library = _library(
        tmp_path,
        reader,
        rt_model,
        wvp={"20220605T101031": 1.2, "20220710T101031": 2.4},
    )

    library.realizations(_observation(), _RES, ["blue"])

    assert sorted({round(state[0], 6) for state in rt_model.states}) == [1.2, 2.4]


def test_a_month_with_no_l2a_water_vapour_at_all_is_an_error(tmp_path: Path) -> None:
    _write_index(tmp_path)
    reader = _FakeReader({"20220605T101031": 0.1, "20220615T101031": 0.2, "20220710T101031": 0.3})
    library = _library(tmp_path, reader, _FakeRTModel(), wvp={"20220710T101031": 2.4})

    with pytest.raises(MissingLibraryInputError, match="No Sentinel-2 L2A water vapour for any"):
        library.realizations(_observation(), _RES, ["blue"])


def test_implausible_water_vapour_is_an_error(tmp_path: Path) -> None:
    _write_index(tmp_path)
    reader = _FakeReader({"20220605T101031": 0.1, "20220615T101031": 0.2, "20220710T101031": 0.3})
    library = _library(
        tmp_path,
        reader,
        _FakeRTModel(),
        wvp=dict.fromkeys(_ALL_WVP, 40.0),
    )

    with pytest.raises(MissingLibraryInputError, match="outside the plausible range"):
        library.realizations(_observation(), _RES, ["blue"])


def test_a_day_without_cams_ozone_is_an_error(tmp_path: Path) -> None:
    _write_index(tmp_path)
    reader = _FakeReader({"20220605T101031": 0.1, "20220615T101031": 0.2, "20220710T101031": 0.3})
    library = _library(tmp_path, reader, _FakeRTModel(), tco3=_FakeTCO3({"2022-06-05": 0.30}))

    with pytest.raises(MissingLibraryInputError, match="No CAMS ozone for 2022-06-15"):
        library.realizations(_observation(), _RES, ["blue"])


def test_a_scene_without_a_dem_is_an_error(tmp_path: Path) -> None:
    _write_index(tmp_path)
    library = LiveL1CSurfaceLibrary(
        tmp_path,
        rt_model=_FakeRTModel(),
        dem_path=None,
        cams_data_path=_CAMS,
        toa_reader=_FakeReader({"20220605T101031": 0.1}),
        wvp_reader=_FakeWVPReader(_ALL_WVP),
        tco3_source=_FakeTCO3(),
    )

    with pytest.raises(MissingLibraryInputError, match="no DEM is configured"):
        library.realizations(_observation(), _RES, ["blue"])


def test_a_sea_level_dem_sentinel_is_also_rejected(tmp_path: Path) -> None:
    _write_index(tmp_path)
    library = LiveL1CSurfaceLibrary(
        tmp_path,
        rt_model=_FakeRTModel(),
        dem_path="none",
        cams_data_path=_CAMS,
        toa_reader=_FakeReader({"20220605T101031": 0.1}),
        wvp_reader=_FakeWVPReader(_ALL_WVP),
        tco3_source=_FakeTCO3(),
    )

    with pytest.raises(MissingLibraryInputError, match="real terrain height"):
        library.realizations(_observation(), _RES, ["blue"])


def test_a_scene_without_a_cams_archive_is_an_error(tmp_path: Path) -> None:
    _write_index(tmp_path)
    library = LiveL1CSurfaceLibrary(
        tmp_path,
        rt_model=_FakeRTModel(),
        dem_path=_DEM,
        cams_data_path=None,
        toa_reader=_FakeReader({"20220605T101031": 0.1}),
        wvp_reader=_FakeWVPReader(_ALL_WVP),
    )

    with pytest.raises(MissingLibraryInputError, match="no CAMS source is configured"):
        library.realizations(_observation(), _RES, ["blue"])


def test_a_winner_index_on_a_different_grid_is_rejected(tmp_path: Path) -> None:
    # Half-size winner planes: the index carries no CRS/transform, so nothing
    # can prove it covers the same ground. It must fail, not be resampled.
    _write_index(tmp_path, winners=np.zeros((2, 5, 5), dtype=np.int16))
    library = _library(tmp_path, _FakeReader({}), _FakeRTModel())

    with pytest.raises(ValueError, match="cannot be co-registered"):
        library.realizations(_observation(), _RES, ["blue"])


# --------------------------------------------------------------------------- #
# CAMS ozone source
# --------------------------------------------------------------------------- #
#: CAMS publishes gtco3 in kg m^-2; this is the atm-cm value the scale yields.
_KG_M2_PER_DU = 2.1415e-5


def _cams_dataset(
    tco3_atm_cm: float,
    *,
    hour: int = 10,
    longitudes: Sequence[float] = (14.0, 15.0, 16.0),
) -> Any:
    """A one-day CAMS file shaped like the archive's: gtco3 in kg m^-2."""
    raw = tco3_atm_cm * _KG_M2_PER_DU * 1000.0
    values = np.zeros((24, 3, len(longitudes)), dtype=np.float32)
    # Only ``hour`` carries the value, so reading the wrong period yields 0.
    values[hour] = raw
    return xr.Dataset(
        {"gtco3": (("forecast_period", "latitude", "longitude"), values)},
        coords={
            "forecast_period": np.arange(24.0),
            "latitude": np.array([50.0, 49.5, 49.0]),
            "longitude": np.asarray(longitudes, dtype=float),
        },
    )


def _cams_source(monkeypatch: pytest.MonkeyPatch, *, data: Any, calls: list[Any] | None) -> Any:
    source = CAMSTCO3Source("/gws/example/cams")

    def _load(when: datetime) -> Any:
        if calls is not None:
            calls.append(when)
        return data

    monkeypatch.setattr(source._provider, "_load_cams_data", _load)
    return source


def test_cams_ozone_reads_the_scene_point_at_the_overpass_hour_and_caches_it(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[Any] = []
    source = _cams_source(monkeypatch, data=_cams_dataset(0.29), calls=calls)

    first = source("2022-06-05", lon=15.0, lat=49.5)
    second = source("2022-06-05", lon=15.0, lat=49.5)

    # Converted out of kg m^-2 by the provider's own published scale.
    assert first == pytest.approx(0.29, abs=1e-6)
    assert second == first
    # Cached: the second call never re-reads the archive.
    assert len(calls) == 1
    # Sampled at the ~10:30 local-solar overpass, not at midnight (lon 15E -> 10).
    assert calls[0] == datetime(2022, 6, 5, 10)


def test_cams_ozone_handles_a_0_to_360_longitude_grid(monkeypatch: pytest.MonkeyPatch) -> None:
    # A western-hemisphere scene on a 0..360 grid must not fall off the edge.
    source = _cams_source(
        monkeypatch, data=_cams_dataset(0.33, hour=14, longitudes=(0.0, 180.0, 300.0)), calls=None
    )

    assert source("2022-06-05", lon=-60.0, lat=49.5) == pytest.approx(0.33, abs=1e-6)


def test_cams_ozone_missing_for_a_day_is_an_error(monkeypatch: pytest.MonkeyPatch) -> None:
    source = _cams_source(monkeypatch, data=None, calls=None)

    with pytest.raises(MissingLibraryInputError, match="No CAMS ozone for 2022-06-05"):
        source("2022-06-05", lon=15.0, lat=49.5)


def test_cams_file_without_the_ozone_variable_is_an_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = _cams_source(monkeypatch, data=xr.Dataset({"aod550": ((), 0.1)}), calls=None)

    with pytest.raises(MissingLibraryInputError, match="no 'gtco3' variable"):
        source("2022-06-05", lon=15.0, lat=49.5)


def test_cams_ozone_outside_the_plausible_range_is_an_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = _cams_source(monkeypatch, data=_cams_dataset(0.9), calls=None)

    with pytest.raises(MissingLibraryInputError, match="outside the plausible"):
        source("2022-06-05", lon=15.0, lat=49.5)


# --------------------------------------------------------------------------- #
# RT space
# --------------------------------------------------------------------------- #
def test_rt_space_is_the_backend_default_without_a_species_mode(tmp_path: Path) -> None:
    library = _library(tmp_path, _FakeReader({}), _FakeRTModel(RTSetupConfig()))

    assert library.rt_space == RTSpace(backend="sixs", aerosol="default")


def test_rt_space_reports_the_solvers_scene_adaptive_rule(tmp_path: Path) -> None:
    solver_config = SimpleNamespace(surface_driven_aerosol_species="cci_climatology_exact")
    library = _library(
        tmp_path, _FakeReader({}), _FakeRTModel(RTSetupConfig()), solver_config=solver_config
    )

    # The rule, not one scene's resolved mixture — exactly what the solve declares.
    assert library.rt_space == RTSpace(backend="sixs", aerosol="cci_climatology_exact")


def test_species_mode_corrects_through_the_scenes_cci_mixture(tmp_path: Path) -> None:
    _write_index(tmp_path)
    reader = _FakeReader({"20220605T101031": 0.1, "20220710T101031": 0.3})
    rt_model = _FakeRTModel(RTSetupConfig())
    solver_config = SimpleNamespace(surface_driven_aerosol_species="cci_climatology_exact")
    library = _library(tmp_path, reader, rt_model, solver_config=solver_config)

    library.realizations(_observation(), _RES, ["blue"])

    # The base backend was cloned once with the scene's resolved CCI mixture, and
    # the correction ran on the clone (the base one saw no RT calls).
    assert len(rt_model.clones) == 1
    assert rt_model.aot_nodes == []
    assert rt_model.clones[0].aot_nodes
    aerosol = rt_model.clones[0].rt_setup.aerosol
    assert aerosol.profile == "multimodal_log_normal"
    assert aerosol.distribution is not None


def test_an_unknown_species_mode_is_rejected(tmp_path: Path) -> None:
    _write_index(tmp_path)
    solver_config = SimpleNamespace(surface_driven_aerosol_species="something_else")
    library = _library(
        tmp_path, _FakeReader({}), _FakeRTModel(RTSetupConfig()), solver_config=solver_config
    )

    with pytest.raises(ValueError, match="species mode"):
        library.realizations(_observation(), _RES, ["blue"])


def test_species_mode_needs_a_cloneable_rt_backend(tmp_path: Path) -> None:
    _write_index(tmp_path)
    solver_config = SimpleNamespace(surface_driven_aerosol_species="cci_climatology_exact")
    library = _library(
        tmp_path, _FakeReader({}), _FakeRTModel(rt_setup=None), solver_config=solver_config
    )

    with pytest.raises(ValueError, match="with_rt_setup"):
        library.realizations(_observation(), _RES, ["blue"])


# --------------------------------------------------------------------------- #
# wiring: config -> surface-library resolution -> l1c source guard
# --------------------------------------------------------------------------- #
def _config(**surface_prior: Any) -> Any:
    from siac.config import SIACConfig

    return SIACConfig.model_validate(
        {"algorithms": {"surface_prior": {"method": "bestpixel", **surface_prior}}}
    )


def test_config_accepts_a_live_index_path(tmp_path: Path) -> None:
    config = _config(bestpixel_source="l1c", live_l1c_index_path=str(tmp_path))

    assert config.algorithms.surface_prior.live_l1c_index_path == tmp_path


def test_resolver_builds_the_live_library_from_the_index_path(tmp_path: Path) -> None:
    from siac.adapters.bestpixel_surface_prior import _resolve_surface_library

    config = _config(bestpixel_source="l1c", live_l1c_index_path=str(tmp_path))
    config.paths.dem = _DEM
    config.providers.atmo.data_path = _CAMS

    library = _resolve_surface_library(config, rt_model=_FakeRTModel(RTSetupConfig()))

    assert isinstance(library, LiveL1CSurfaceLibrary)
    # The run's own DEM and CAMS archive supply the measured terrain and ozone.
    assert (library._dem_path, library._cams_data_path) == (_DEM, _CAMS)


def test_resolver_prefers_a_prepared_store_over_the_live_build(tmp_path: Path) -> None:
    from siac.adapters.bestpixel_surface_prior import _resolve_surface_library
    from siac.adapters.surface_library import PreparedSurfaceLibrary

    library = _resolve_surface_library(
        _config(
            bestpixel_source="l1c",
            prepared_library_path=str(tmp_path / "prepared"),
            live_l1c_index_path=str(tmp_path / "index"),
        ),
        rt_model=_FakeRTModel(RTSetupConfig()),
    )

    assert isinstance(library, PreparedSurfaceLibrary)


def test_resolver_needs_an_rt_model_for_the_live_build(tmp_path: Path) -> None:
    from siac.adapters.bestpixel_surface_prior import _resolve_surface_library

    with pytest.raises(ValueError, match="requires an RT backend"):
        _resolve_surface_library(_config(bestpixel_source="l1c", live_l1c_index_path=str(tmp_path)))


def test_l1c_source_is_accepted_once_a_live_index_is_configured(tmp_path: Path) -> None:
    from siac.app._assembly_surface import make_bestpixel_surface_prior_fn

    surface_prior_fn = make_bestpixel_surface_prior_fn(
        _config(bestpixel_source="l1c", live_l1c_index_path=str(tmp_path))
    )

    # The guard no longer fires; the run gets as far as looking for the scene's
    # index, which this empty store does not hold.
    with pytest.raises(FileNotFoundError, match="mosaic index"):
        surface_prior_fn(_observation(), None, _FakeRTModel(RTSetupConfig()), _RES)


def test_l1c_source_without_any_library_is_still_rejected() -> None:
    from siac.app._assembly_surface import make_bestpixel_surface_prior_fn

    surface_prior_fn = make_bestpixel_surface_prior_fn(_config(bestpixel_source="l1c"))

    with pytest.raises(NotImplementedError, match="live_l1c_index_path"):
        surface_prior_fn(_observation(), None, _FakeRTModel(RTSetupConfig()), _RES)
