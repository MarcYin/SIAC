from __future__ import annotations

import json

from tools.aeronet_validation.summarize_low_cloud_aod_screen import _load_records


def test_load_records_normalizes_nested_multigrid_result(tmp_path) -> None:
    result_dir = tmp_path / "runs" / "monthly_database" / "site__scene"
    result_dir.mkdir(parents=True)
    (result_dir / "result.json").write_text(
        json.dumps(
            {
                "matchup_id": "site__scene",
                "status": "ok",
                "aot_window_mean": 0.3,
            }
        ),
        encoding="utf-8",
    )

    records = _load_records(tmp_path, truth_by_id={"site__scene": 0.4})

    assert records["site__scene"]["status"] == "OK"
    assert records["site__scene"]["retrieved"] == 0.3
    assert records["site__scene"]["truth"] == 0.4
