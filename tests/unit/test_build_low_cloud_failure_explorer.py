from __future__ import annotations

import numpy as np
from tools.aeronet_validation.build_low_cloud_failure_explorer import (
    _argmin_map,
    _candidate_outputs,
    _jsonable,
    _normalise_curve,
)


def test_argmin_map_preserves_missing_pixels() -> None:
    cube = np.asarray(
        [
            [[3.0, np.nan], [2.0, 4.0]],
            [[1.0, np.nan], [3.0, 2.0]],
        ]
    )

    result = _argmin_map(cube, np.asarray([0.1, 0.2]))

    np.testing.assert_allclose(result[[0, 1], [0, 1]], [0.2, 0.2])
    assert np.isnan(result[0, 1])
    assert result[1, 0] == 0.1


def test_normalised_curve_keeps_minimum_and_handles_missing_nodes() -> None:
    result = _normalise_curve(np.asarray([4.0, 2.0, np.nan, 3.0]))

    assert result[1] == 0.0
    assert np.isnan(result[2])
    assert result[0] > result[3] > result[1]


def test_jsonable_replaces_non_finite_values() -> None:
    result = _jsonable(
        {
            "values": np.asarray([1.0, np.nan, np.inf]),
            "count": np.int64(3),
            "flag": np.bool_(True),
        }
    )

    assert result == {"values": [1.0, None, None], "count": 3, "flag": True}


def test_candidate_outputs_keep_source_kind_and_expected_error_flag() -> None:
    current = {
        "truth": 0.5,
        "retrieved": 0.3,
        "retrieval_extraction": {"station": 0.4, "winmed": 0.5},
        "solver": {
            "surface_cost_curve_min_aot": 0.45,
            "surface_band_B02_argmin_aot": 0.4,
            "surface_band_B03_argmin_aot": 0.5,
            "surface_band_B04_argmin_aot": 0.6,
        },
    }
    variants = {
        "b03_shape": {
            "case": {"status": "OK", "retrieved": 0.55},
        },
        "b03_auto2": {
            "case": {"status": "FAILED", "retrieved": 0.5},
        },
    }

    rows = _candidate_outputs(
        "case",
        failure_row={"oof_aod": 0.7},
        current=current,
        baseline={"retrieved": 0.65},
        prior={"cams_aot": 0.4},
        variants=variants,
        solver_cams=0.35,
        no_backstop=0.5,
    )

    by_label = {row["label"]: row for row in rows}
    assert by_label["Current scene mean"]["within_ee"] is False
    assert by_label["1.5 km median"]["within_ee"] is True
    assert by_label["Surface-prior anchor CAMS"]["value"] == 0.4
    assert by_label["Solver CAMS backstop"]["value"] == 0.35
    assert by_label["No-backstop replay, scene mean"]["within_ee"] is True
    assert by_label["B02/B03/B04 spectral shape"] == {
        "kind": "variant",
        "label": "B02/B03/B04 spectral shape",
        "value": 0.55,
        "within_ee": True,
    }
    assert "B02/B03/B04 auto2" not in by_label
