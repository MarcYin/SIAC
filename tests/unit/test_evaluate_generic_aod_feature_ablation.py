from __future__ import annotations

from tools.aeronet_validation.evaluate_generic_aod_feature_ablation import (
    ablation_schemas,
)


def test_ablation_schemas_remove_only_declared_feature_groups() -> None:
    names = (
        "siac_retrieved",
        "time_utc_sin",
        "time_year_sin",
        "context_cams_total_aerosol_optical_depth_at_550nm_surface",
        "context_cams_dust_aerosol_optical_depth_at_550nm_surface",
        "consistency_cams_minus_siac",
        "atmo_aot_mean",
        "solver_cost",
    )
    schemas = ablation_schemas(names)
    assert "time_utc_sin" not in schemas["without_utc_phase"]
    assert "time_year_sin" in schemas["without_utc_phase"]
    assert (
        "context_cams_total_aerosol_optical_depth_at_550nm_surface"
        in schemas["without_cams_species_components"]
    )
    assert (
        "context_cams_dust_aerosol_optical_depth_at_550nm_surface"
        not in schemas["without_cams_species_components"]
    )
    assert schemas["retrieval_surface_solver_only"] == (
        "siac_retrieved",
        "solver_cost",
    )
