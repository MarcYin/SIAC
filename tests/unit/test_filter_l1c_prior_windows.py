from __future__ import annotations

import json

import numpy as np
from tools.aeronet_validation.filter_l1c_prior_windows import filter_prior


def test_filter_prior_keeps_selected_composites_and_metadata(tmp_path):
    source = tmp_path / "source.npz"
    destination = tmp_path / "filtered.npz"
    quality = [
        {"window": "a", "selected_aod_median": 0.30},
        {"window": "b", "selected_aod_median": 0.10},
        {"window": "c", "selected_aod_median": 0.20},
        {"window": "d", "selected_aod_median": 0.12},
    ]
    np.savez(
        source,
        comp=np.arange(4 * 2 * 2 * 2).reshape(4, 2, 2, 2),
        realizations=np.array(["a", "b", "c", "d"]),
        clean_quality_json=np.array(json.dumps(quality)),
        epsg=np.array(32632),
    )

    result = filter_prior(source, destination, max_median_aod=0.15, min_windows=3)

    with np.load(destination, allow_pickle=False) as filtered:
        assert filtered["realizations"].tolist() == ["b", "c", "d"]
        assert filtered["comp"].shape[0] == 3
        assert int(filtered["epsg"]) == 32632
        selection = json.loads(str(filtered["clean_window_selection_json"].item()))
    assert selection["n_fallback_windows"] == 1
    assert result["kept_realizations"] == ["b", "c", "d"]
