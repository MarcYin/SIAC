"""Prepared surface library: reading, band mapping and RT-space declaration."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import numpy as np
import pytest

from siac.adapters.surface_library import (
    PREDICTOR_BAND_ORDER,
    PreparedSurfaceLibrary,
    canonical_band_name,
    realization_to_period,
)
from siac.domain.rt_space import RTSpace

if TYPE_CHECKING:
    from pathlib import Path


class _Observation:
    def __init__(self, key: str | None = "SCENE_A") -> None:
        self.metadata = {"scene_key": key} if key else {}


def _write_store(
    root: Path,
    *,
    scene: str = "SCENE_A",
    realizations: int = 3,
    bands: int = 7,
    rt_space: dict[str, str] | None = None,
    band_names: list[str] | None = None,
) -> np.ndarray:
    comp = np.linspace(0.02, 0.4, realizations * bands * 6, dtype=np.float32)
    comp = comp.reshape(realizations, bands, 2, 3)
    payload: dict[str, object] = {
        "comp": comp,
        "epsg": np.array(32637),
        "transform": np.array([60.0, 0.0, 5e5, 0.0, -60.0, 2.47e6]),
        "realizations": np.array([f"2020-{i + 1:02d}" for i in range(realizations)]),
    }
    if band_names is not None:
        payload["band_names"] = np.array(band_names)
    root.mkdir(parents=True, exist_ok=True)
    np.savez(root / f"{scene}.npz", **payload)
    if rt_space is not None:
        (root / "library.json").write_text(json.dumps({"rt_space": rt_space}))
    return comp


def test_canonical_band_names_map_sensor_aliases() -> None:
    assert canonical_band_name("nir08") == "nir"
    assert canonical_band_name("B8A") == "nir"
    assert canonical_band_name("B02") == "blue"
    assert canonical_band_name("blue") == "blue"


def test_library_reads_realizations_with_default_band_order(tmp_path: Path) -> None:
    comp = _write_store(tmp_path)

    realizations = PreparedSurfaceLibrary(tmp_path).realizations(_Observation(), 60.0, [])

    assert len(realizations) == comp.shape[0]
    assert realizations[0].band_names == PREDICTOR_BAND_ORDER
    assert realizations[0].crs == "EPSG:32637"
    assert [r.label for r in realizations] == ["2020-01", "2020-02", "2020-03"]


def test_stored_band_names_are_canonicalized(tmp_path: Path) -> None:
    _write_store(
        tmp_path,
        bands=3,
        band_names=["B01", "B02", "nir08"],
    )

    realizations = PreparedSurfaceLibrary(tmp_path).realizations(_Observation(), 60.0, [])

    assert realizations[0].band_names == ("coastal", "blue", "nir")


def test_rt_space_is_read_from_the_sidecar(tmp_path: Path) -> None:
    _write_store(tmp_path, rt_space={"backend": "sixs", "aerosol": "cci_climatology_exact"})

    library = PreparedSurfaceLibrary(tmp_path)

    assert library.rt_space == RTSpace("sixs", "cci_climatology_exact")


def test_rt_space_is_none_without_a_sidecar(tmp_path: Path) -> None:
    _write_store(tmp_path)

    assert PreparedSurfaceLibrary(tmp_path).rt_space is None


def test_explicit_rt_space_overrides_the_sidecar(tmp_path: Path) -> None:
    _write_store(tmp_path, rt_space={"backend": "lut", "aerosol": "continental_average"})
    space = RTSpace("sixs", "cci_climatology_exact")

    assert PreparedSurfaceLibrary(tmp_path, rt_space=space).rt_space == space


def test_missing_scene_reports_the_store_and_key(tmp_path: Path) -> None:
    _write_store(tmp_path, scene="OTHER")

    with pytest.raises(FileNotFoundError, match="SCENE_A"):
        PreparedSurfaceLibrary(tmp_path).realizations(_Observation(), 60.0, [])


def test_observation_without_a_scene_key_is_rejected(tmp_path: Path) -> None:
    _write_store(tmp_path)

    with pytest.raises(ValueError, match="needs a scene key"):
        PreparedSurfaceLibrary(tmp_path).realizations(_Observation(None), 60.0, [])


def test_band_count_mismatch_is_reported(tmp_path: Path) -> None:
    _write_store(tmp_path, bands=4)

    with pytest.raises(ValueError, match="declares 7 band names"):
        PreparedSurfaceLibrary(tmp_path).realizations(_Observation(), 60.0, [])


def test_realization_round_trips_through_the_composite_payload(tmp_path: Path) -> None:
    comp = _write_store(tmp_path)
    realization = PreparedSurfaceLibrary(tmp_path).realizations(_Observation(), 60.0, [])[0]

    period = realization_to_period(realization)

    recovered = period["bands"]["blue"] * period["reflectance_scale"]
    assert np.allclose(recovered, comp[0, 1], atol=1e-6)
    assert period["grid"]["width"] == 3
    assert period["grid"]["height"] == 2
    assert period["grid"]["crs"] == "EPSG:32637"


def test_pixels_missing_in_every_band_become_nodata() -> None:
    reflectance = np.full((2, 2, 2), 0.2, dtype=np.float32)
    reflectance[:, 0, 0] = np.nan
    from siac.adapters.surface_library import SurfaceLibraryRealization

    period = realization_to_period(
        SurfaceLibraryRealization(
            reflectance=reflectance,
            band_names=("blue", "red"),
            crs="EPSG:32637",
            transform=(60.0, 0.0, 5e5, 0.0, -60.0, 2.47e6),
        )
    )

    assert period["quality"][0, 0] == 65535
    assert period["quality"][1, 1] == 0


def test_non_3d_reflectance_is_rejected() -> None:
    from siac.adapters.surface_library import SurfaceLibraryRealization

    flat = SurfaceLibraryRealization(
        reflectance=np.zeros((2, 2), dtype=np.float32),
        band_names=("blue",),
        crs="EPSG:32637",
        transform=(60.0, 0.0, 5e5, 0.0, -60.0, 2.47e6),
    )

    with pytest.raises(ValueError, match=r"must be \(band, y, x\)"):
        realization_to_period(flat)
