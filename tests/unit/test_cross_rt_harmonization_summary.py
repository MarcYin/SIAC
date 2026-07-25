import hashlib
import json

from tools.aeronet_validation.summarize_cross_rt_harmonization import (
    HISTORY_TAGS,
    _best_complete_candidate,
    _cohort_definition,
    _correction_markdown,
    _correction_model_report,
    _evaluation_split,
    _history_archive_health,
    _html,
    _pair_archive_health,
    _parse_variant_tag,
    _result_commentary,
    _retrieval_metrics,
)


def _metrics(n: int, rate: float | None) -> dict:
    return {
        "n": n,
        "within_ee": round(n * rate) if rate is not None else 0,
        "within_ee_rate": rate,
        "bias": -0.01 if n else None,
        "mae": 0.08 if n else None,
        "rmse": 0.12 if n else None,
    }


def _variant(name: str, n: int, rate: float | None, development_rate: float | None) -> dict:
    return {
        "variant": name,
        "metrics": _metrics(n, rate),
        "development_metrics": _metrics(91, development_rate),
        "holdout_metrics": _metrics(61, rate - 0.02 if rate is not None else None),
        "regimes": {
            "medium_0p25_0p85": {
                **_metrics(n // 2, rate),
                "bias": -0.02 if n else None,
                "mae": 0.09 if n else None,
            }
        },
        "transitions_vs_identity": {"gain": 4, "loss": 2},
    }


def _summary(candidate_n: int) -> dict:
    return {
        "cohort_size": 152,
        "evaluation_split": {"development_n": 91, "holdout_n": 61},
        "surface": [
            {"model": "identity", "cap": "uncapped", "visible_scene_mae": 0.010},
            {"model": "cross_rt_terrain_a1", "cap": "cap_0.015", "visible_scene_mae": 0.006},
        ],
        "retrieval": [
            _variant(
                name,
                152 if name == "identity_daily" else candidate_n,
                0.74 if name == "identity_daily" else (0.88 if index == 1 else 0.80),
                0.74
                if name == "identity_daily"
                else (0.84 if index == 1 else (0.90 if index == 2 else 0.80)),
            )
            for index, name in enumerate(HISTORY_TAGS)
        ],
    }


def test_report_does_not_rank_incomplete_retrievals() -> None:
    summary = _summary(candidate_n=151)

    assert _best_complete_candidate(summary) is None
    comments = _result_commentary(summary)
    assert "40.0% below" in comments[0]
    assert "must not yet be ranked" in comments[1]


def test_retrieval_rate_keeps_non_ok_cases_in_strict_denominator() -> None:
    metrics = _retrieval_metrics(
        [
            {
                "status": "OK",
                "truth": 0.3,
                "retrieved": 0.3,
                "within_ee": True,
                "solver": {
                    "surface_cloud_mask_bypassed": True,
                    "surface_water_mask_bypassed": True,
                },
            },
            {"status": "OK", "truth": 0.4, "retrieved": 0.8, "within_ee": False},
            {"status": "MISSING", "truth": 0.5},
        ]
    )

    assert metrics["within_ee"] == 1
    assert metrics["n"] == 2
    assert metrics["expected"] == 3
    assert metrics["within_ee_rate"] == 1 / 3
    assert metrics["valid_within_ee_rate"] == 1 / 2
    assert metrics["cloud_mask_bypass_count"] == 1
    assert metrics["water_mask_bypass_count"] == 1


def test_report_comments_only_on_complete_best_candidate() -> None:
    summary = _summary(candidate_n=152)

    best = _best_complete_candidate(summary)
    assert best is not None
    assert best["variant"] == HISTORY_TAGS[2]
    comments = _result_commentary(summary)
    assert any("90.0% within EE" in comment for comment in comments)
    assert any("+16.0 percentage points" in comment for comment in comments)
    assert any("7.0 percentage points below" in comment for comment in comments)


def test_pair_health_separates_sparse_filtering_from_other_errors(tmp_path) -> None:
    matchup_id = "Example__T00XXX_20240101T000000"
    (tmp_path / f"{matchup_id}.json").write_text(
        json.dumps(
            {
                "status": "OK",
                "attempted_scenes": 4,
                "successful_scenes": 2,
                "errors": [
                    {
                        "error_type": "RuntimeError",
                        "reason": "only 12 paired clear-land pixels",
                    },
                    {"error_type": "TimeoutError", "reason": "service timeout"},
                ],
            }
        ),
        encoding="utf-8",
    )

    health = _pair_archive_health(tmp_path, [matchup_id, "Missing"])

    assert health["status_counts"] == {"OK": 1}
    assert health["missing_or_malformed"] == 1
    assert health["retained_acquisitions"] == 2
    assert health["retained_tile_scene_components"] == 2
    assert health["sparse_clear_land_rejections"] == 1
    assert health["other_rejections"] == 1


def test_evaluation_split_uses_locked_manifest(tmp_path) -> None:
    manifest = tmp_path / "summary.json"
    manifest.write_text(
        json.dumps(
            {
                "cohort": {
                    "split_seed": 42,
                    "holdout_folds": [2, 3],
                    "fold_by_matchup_id": {"development": 0, "holdout": 3},
                }
            }
        ),
        encoding="utf-8",
    )

    assignment, provenance = _evaluation_split(manifest, ["development", "holdout"])

    assert assignment == {"development": "development", "holdout": "holdout"}
    assert provenance["development_n"] == 1
    assert provenance["holdout_n"] == 1


def test_history_health_reports_same_source_fallbacks(tmp_path) -> None:
    mapped_output = tmp_path / "mapped.npz"
    fallback_output = tmp_path / "fallback.npz"
    mapped_output.touch()
    fallback_output.touch()
    (tmp_path / "mapped.json").write_text(
        json.dumps(
            {
                "status": "OK",
                "matchup_id": "mapped",
                "application": "per acquisition before tile mosaic and temporal median",
                "outputs": {"candidate": str(mapped_output)},
                "errors": [{"reason": "one tile was sparse"}],
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "fallback.json").write_text(
        json.dumps(
            {
                "status": "OK",
                "matchup_id": "fallback",
                "application": "uncorrected single-source L2A fallback",
                "mapping_applied": False,
                "skip_reason": "insufficient_harmonized_history_windows",
                "low_temporal_support": True,
                "outputs": {"candidate": str(fallback_output)},
            }
        ),
        encoding="utf-8",
    )

    health = _history_archive_health(tmp_path, ["mapped", "fallback", "missing"])

    assert health["terminal_cases"] == 2
    assert health["per_acquisition_harmonized_cases"] == 1
    assert health["uncorrected_single_source_fallback_cases"] == 1
    assert health["fallback_reason_counts"] == {"insufficient_harmonized_history_windows": 1}
    assert health["present_candidate_outputs"] == 2
    assert health["nonfatal_scene_errors"] == 1


def test_cohort_definition_uses_masked_retrieval_cloud_fraction(tmp_path) -> None:
    mids = tmp_path / "lowcloud.txt"
    serialized = "first\nsecond\n"
    mids.write_text(serialized, encoding="utf-8")
    manifest = tmp_path / "lowcloud.json"
    manifest.write_text(
        json.dumps(
            {
                "selected_mids_file": str(mids),
                "selected_mids_sha256": hashlib.sha256(serialized.encode()).hexdigest(),
                "cloud_source_directory": "masked-r2-results",
                "cloud_field": "cloud_frac",
                "comparison": "strictly_less_than",
                "threshold_fraction": 0.2,
                "campaign_count": 250,
            }
        ),
        encoding="utf-8",
    )

    definition = _cohort_definition(
        manifest,
        ["first", "second"],
        [{"cloud_cover": "5.0"}, {"cloud_cover": "30.0"}],
    )

    assert definition["cloud_field"] == "cloud_frac"
    assert definition["selected_count"] == 2
    assert definition["metadata_scene_cloud_cover_below_20_count"] == 1


def _band_error_metrics(bias: float, mae: float) -> dict:
    return {
        "n": 1000,
        "bias": bias,
        "mae": mae,
        "rmse": mae * 1.3,
        "p95_abs": mae * 2.5,
        "scene_bias": bias,
        "scene_mae": mae,
        "scene_rmse": mae * 1.2,
    }


def _correction_fixture(tmp_path) -> tuple[dict, dict, "object"]:
    bands = ("blue", "green", "red")
    band_models = {
        band: {
            "mean": [0.1, 0.0],
            "scale": [0.05, 0.1],
            "coef": [0.01, -0.02],
            "intercept": 0.001,
            "feature_names": [f"l2a_{band}", "delta_aot_maiac_minus_sen2cor"],
        }
        for band in bands
    }
    artifact = {
        "target": "same-day L1C corrected with MAIAC AOD and the current libRadtran LUT",
        "crossfit_protocol": {"holdout_model_fold": 5},
        "models": {
            "cross_rt_terrain_a1": {
                "alpha": 1.0,
                "weighting_group": "acquisition",
                "folds": {"5": band_models},
                "full": band_models,
                "oof_residual_scale": {
                    band: {"median": 0.0001, "mad_to_sigma": 0.005, "rmse": 0.007} for band in bands
                },
            }
        },
    }
    surface_metrics = {
        "identity": {band: _band_error_metrics(-0.003, 0.008) for band in bands},
        "candidates": {
            "cross_rt_terrain_a1": {
                "cap_0.030": {band: _band_error_metrics(-0.001, 0.005) for band in bands}
            }
        },
    }
    scene_csv = tmp_path / "surface_scene_metrics.csv"
    header = ["delta_aot_maiac_minus_sen2cor"]
    for band in bands:
        header.extend([f"identity_{band}_bias", f"cross_rt_terrain_a1_cap_0p030_{band}_bias"])
    rows = [
        {"delta_aot_maiac_minus_sen2cor": 0.15, "correction": 0.005},
        {"delta_aot_maiac_minus_sen2cor": 0.0, "correction": 0.0},
        {"delta_aot_maiac_minus_sen2cor": -0.15, "correction": -0.004},
    ]
    lines = [",".join(header)]
    for row in rows:
        values = [str(row["delta_aot_maiac_minus_sen2cor"])]
        for _band in bands:
            values.extend(["-0.006", str(-0.006 + row["correction"])])
        lines.append(",".join(values))
    scene_csv.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return artifact, surface_metrics, scene_csv


def test_parse_variant_tag_extracts_model_mode_and_cap() -> None:
    assert _parse_variant_tag("cross_rt_terrain_a1_solver_cap0p030") == (
        "cross_rt_terrain_a1",
        "solver",
        "cap_0.030",
        0.030,
    )
    assert _parse_variant_tag("cross_rt_terrain_a1_all_cap0p015") == (
        "cross_rt_terrain_a1",
        "all",
        "cap_0.015",
        0.015,
    )
    assert _parse_variant_tag("identity_daily") is None


def test_correction_model_report_shows_derivation_effect_and_accuracy(tmp_path) -> None:
    artifact, surface_metrics, scene_csv = _correction_fixture(tmp_path)

    report = _correction_model_report(
        artifact,
        surface_metrics,
        scene_csv,
        "cross_rt_terrain_a1_solver_cap0p030",
        selection_basis="development-selected complete variant",
    )

    assert report is not None
    assert report["corrected_bands"] == ["blue", "green", "red"]
    assert report["regression"]["feature_count"] == 2
    blue = next(row for row in report["effect"]["bands"] if row["band"] == "blue")
    assert abs(blue["scene_mae_change_pct"] - (-37.5)) < 1.0e-9
    assert blue["oof_residual_sigma_mad"] == 0.005
    features = {row["feature"] for row in report["coefficients"]["rows"]}
    assert features == {"l2a_band", "delta_aot_maiac_minus_sen2cor"}
    shared = next(row for row in report["coefficients"]["rows"] if row["feature"] == "l2a_band")
    assert set(shared["coefficients"]) == {"blue", "green", "red"}
    applied = report["scene_corrections"]["per_band"]["blue"]
    assert applied["n"] == 3
    assert applied["median"] == 0.0
    bins = {row["label"]: row for row in report["scene_corrections"]["delta_aot_bins"]}
    assert bins["> 0.10"]["n"] == 1
    assert abs(bins["> 0.10"]["mean_correction"]["blue"] - 0.005) < 1.0e-12
    assert abs(bins["> 0.10"]["corrected_blue_bias_mean"] - (-0.001)) < 1.0e-12
    assert bins["-0.02 to 0.02"]["n"] == 1
    assert bins["< -0.10"]["n"] == 1


def test_correction_model_report_rejects_unknown_variants(tmp_path) -> None:
    artifact, surface_metrics, scene_csv = _correction_fixture(tmp_path)

    assert (
        _correction_model_report(
            artifact,
            surface_metrics,
            scene_csv,
            "identity_daily",
            selection_basis="development-selected complete variant",
        )
        is None
    )
    assert (
        _correction_model_report(
            artifact,
            surface_metrics,
            scene_csv,
            "unknown_model_solver_cap0p030",
            selection_basis="development-selected complete variant",
        )
        is None
    )


def test_correction_markdown_renders_all_three_sections(tmp_path) -> None:
    artifact, surface_metrics, scene_csv = _correction_fixture(tmp_path)
    report = _correction_model_report(
        artifact,
        surface_metrics,
        scene_csv,
        "cross_rt_terrain_a1_solver_cap0p030",
        selection_basis="development-selected complete variant",
    )

    text = "\n".join(_correction_markdown(report))

    assert "## Correction model: derivation" in text
    assert "## Correction effect on the L2A history" in text
    assert "## Correction accuracy" in text
    assert "ridge regression" in text
    assert "delta_aot_maiac_minus_sen2cor" in text
    assert "-37.5%" in text


def test_html_report_is_static_and_links_machine_readable_outputs() -> None:
    document = _html(
        "# Experiment\n\n## Retrieval validation\n\n"
        "| Variant | Within EE |\n|---|---:|\n| identity | 74.0% |"
    )

    assert document.startswith("<!doctype html>")
    assert "<table>" in document
    assert 'href="summary.json"' in document
    assert 'href="retrieval_metrics.csv"' in document
    assert 'href="correction_effect.csv"' in document
    assert 'href="visual.html"' in document
    assert "<script" not in document
