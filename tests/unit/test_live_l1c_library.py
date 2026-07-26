"""Live L1C surface library: index reading, mosaicking, correction, RT space.

Every test injects fakes. The real library reads L1C imagery from the public
Google Cloud Sentinel-2 bucket and corrects it through a compiled 6S backend;
here a fake reader supplies flat TOA planes and a fake RT backend returns
``xap = 1 + aot`` (so a correct pipeline yields ``(1 + aod) * toa`` exactly),
which pins the per-pixel AOD selection, the mosaic assembly, and the grid
handling without touching the network or 6S.
"""

from __future__ import annotations

import json
from datetime import datetime
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any

import numpy as np
import pytest

from siac.adapters.live_l1c_library import (
    LIBRARY_BAND_NAMES,
    IndexImage,
    LiveL1CSurfaceLibrary,
    _baseline_offset,
    _keep_realizations,
    _resample_winners,
    correct_toa_realization,
    median_scalar,
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

    def compute_coefficients_multi(
        self, geometry: Any, atmo_state: Any, bands: Sequence[Any]
    ) -> list[_FakeCoefficients]:
        _ = geometry
        aot = float(np.ravel(np.asarray(atmo_state.aot.values))[0])
        self.aot_nodes.append(aot)
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
        day_scalars=np.array(
            [json.dumps({"day": day, "aod": aod}) for day, aod in scalars.items()]
        ),
    )
    return path


def _library(root: Path, reader: _FakeReader, rt_model: _FakeRTModel, **kwargs: Any) -> Any:
    return LiveL1CSurfaceLibrary(root, rt_model=rt_model, toa_reader=reader, **kwargs)


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


def test_median_scalar_ignores_non_finite_values() -> None:
    assert median_scalar(None) is None
    assert median_scalar(np.array([np.nan, np.nan])) is None
    assert median_scalar(np.array([1.0, np.nan, 3.0])) == 2.0
    assert median_scalar(SimpleNamespace(values=np.array([[2.0, 4.0]]))) == 3.0


def test_winner_resampling_only_runs_on_a_shape_mismatch() -> None:
    winners = np.arange(2 * 4 * 4, dtype=np.int16).reshape(2, 4, 4)
    assert _resample_winners(winners, 4, 4) is winners

    resampled = _resample_winners(winners, 8, 8)
    assert resampled.shape == (2, 8, 8)
    # Nearest-neighbour upsampling repeats each source cell.
    assert resampled[0, 0, 0] == winners[0, 0, 0]
    assert resampled[0, 7, 7] == winners[0, 3, 3]


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
        tcwv=2.0,
        tco3=0.3,
        elevation_km=0.0,
    )

    # xap = 1 + aot, xbp = xcp = 0 -> corrected = (1 + aod) * toa, exactly, because
    # linear interpolation of a linear curve is exact between the AOD nodes.
    np.testing.assert_allclose(corrected[:, :, 0], 1.2 * 0.1, rtol=1e-5)
    np.testing.assert_allclose(corrected[:, :, 1], 1.4 * 0.1, rtol=1e-5)
    # One 6S invocation per AOD node, spanning the mosaic's own AOD range.
    assert len(rt_model.aot_nodes) == 8
    assert min(rt_model.aot_nodes) <= 0.2 and max(rt_model.aot_nodes) >= 0.4


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
            tcwv=2.0,
            tco3=0.3,
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


def test_only_the_winning_acquisitions_are_read_and_on_the_solver_grid(
    tmp_path: Path,
) -> None:
    # June's second acquisition wins nothing: it must never be fetched.
    winners = np.zeros((2, 10, 10), dtype=np.int16)
    _write_index(tmp_path, winners=winners)
    reader = _FakeReader({"20220605T101031": 0.1, "20220710T101031": 0.3})
    library = _library(tmp_path, reader, _FakeRTModel())

    library.realizations(_observation(), _RES, ["blue"])

    tokens = sorted(call["sensing_token"] for call in reader.calls)
    assert tokens == ["20220605T101031", "20220710T101031"]
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

    with pytest.raises(RuntimeError, match="no usable realizations"):
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


def test_a_failing_maiac_gate_falls_back_to_a_default_aod(tmp_path: Path) -> None:
    _write_index(tmp_path, day_aod={})
    reader = _FakeReader({"20220605T101031": 0.1, "20220710T101031": 0.3})

    def _gate(bounds: Any, crs: str, periods: Any) -> dict[str, float]:
        raise RuntimeError("earthaccess down")

    library = _library(tmp_path, reader, _FakeRTModel(), maiac_day_aod=_gate)
    realizations = library.realizations(_observation(), _RES, ["blue"])

    np.testing.assert_allclose(realizations[0].reflectance[:, 0, 0], 1.1 * 0.1, rtol=1e-5)


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

    atmo_prior = SimpleNamespace(
        tcwv=SimpleNamespace(values=np.array([[1.5]])),
        tco3=SimpleNamespace(values=np.array([[0.32]])),
        elevation=SimpleNamespace(values=np.array([[1.2]])),
    )
    library = _resolve_surface_library(
        _config(bestpixel_source="l1c", live_l1c_index_path=str(tmp_path)),
        rt_model=_FakeRTModel(RTSetupConfig()),
        atmo_prior=atmo_prior,
    )

    assert isinstance(library, LiveL1CSurfaceLibrary)
    # The atmospheric prior sets the correction's scalar water/ozone/elevation.
    assert (library._tcwv, library._tco3, library._elevation_km) == (1.5, 0.32, 1.2)


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
