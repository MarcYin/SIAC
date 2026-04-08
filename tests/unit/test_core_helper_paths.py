from __future__ import annotations

import dataclasses
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

import siac._rust_compat as rust_compat
from siac.algorithms.rt.lut import rsrf_kernel
from siac.domain import SensorBand
from siac.domain.spectral import RelativeSpectralResponse
from siac.errors import ValidationError
from siac.runtime import CorrectionDiagnostics, CorrectionResult
from siac.runtime.validation import (
    spatial_shape,
    validate_atmospheric_state,
    validate_correction_result,
    validate_observation_bundle,
    validate_solved_atmosphere,
    validate_solver_input_bundle,
)


def _rsrf() -> RelativeSpectralResponse:
    return RelativeSpectralResponse.from_tabulated(
        sensor_id="MSI",
        satellite_id="S2A",
        band_name="B02",
        wavelengths_nm=np.array([440.0, 450.0, 460.0, 470.0, 480.0], dtype=np.float32),
        response=np.array([0.0, 0.5, 1.0, 0.5, 0.0], dtype=np.float32),
    )


def test_validation_helpers_cover_remaining_error_paths(
    mock_observation_bundle,
    mock_atmospheric_state,
    mock_solver_input_bundle,
    mock_solved_atmosphere,
) -> None:
    one_dim = xr.Dataset({"v": xr.DataArray(np.array([1.0], dtype=np.float32), dims=["only"])})
    with pytest.raises(ValueError, match="fewer than 2 dimensions"):
        spatial_shape(one_dim)

    bad_obs = dataclasses.replace(mock_observation_bundle, bounds=(1.0, 2.0, 3.0))
    with pytest.raises(ValidationError, match="bounds must have 4 elements"):
        validate_observation_bundle(bad_obs)

    bad_obs = dataclasses.replace(mock_observation_bundle, crs="")
    with pytest.raises(ValidationError, match="crs must be a non-empty string"):
        validate_observation_bundle(bad_obs)

    bad_geom = dataclasses.replace(
        mock_observation_bundle.geometry,
        sza=xr.DataArray(np.zeros(5, dtype=np.float32), dims=["y"]),
    )
    bad_obs = dataclasses.replace(mock_observation_bundle, geometry=bad_geom)
    with pytest.raises(ValidationError, match="geometry.sza"):
        validate_observation_bundle(bad_obs)

    bad_sib = dataclasses.replace(mock_solver_input_bundle, aux_resolution_m=0.0)
    with pytest.raises(ValidationError, match="aux_resolution_m must be positive"):
        validate_solver_input_bundle(bad_sib)

    bad_sib = dataclasses.replace(mock_solver_input_bundle, aerosol_resolution_m=0.0)
    with pytest.raises(ValidationError, match="aerosol_resolution_m must be positive"):
        validate_solver_input_bundle(bad_sib)

    bad_solved = dataclasses.replace(mock_solved_atmosphere, converged="yes")
    with pytest.raises(ValidationError, match="converged must be a boolean"):
        validate_solved_atmosphere(bad_solved)

    bad_solved = dataclasses.replace(mock_solved_atmosphere, n_iterations=1.5)
    with pytest.raises(ValidationError, match="n_iterations must be an int"):
        validate_solved_atmosphere(bad_solved)

    bad_solved = dataclasses.replace(mock_solved_atmosphere, cost_final="bad")
    with pytest.raises(ValidationError, match="cost_final must be numeric"):
        validate_solved_atmosphere(bad_solved)

    bad_solved = dataclasses.replace(mock_solved_atmosphere, n_iterations=-1)
    with pytest.raises(ValidationError, match="n_iterations must be non-negative"):
        validate_solved_atmosphere(bad_solved)

    result = CorrectionResult(
        boa=mock_observation_bundle.toa,
        boa_unc=None,
        aot=mock_solved_atmosphere.aot,
        tcwv=mock_solved_atmosphere.tcwv,
        cloud_mask=mock_observation_bundle.cloud_mask,
        metadata="bad",  # type: ignore[arg-type]
        diagnostics=CorrectionDiagnostics(processing_time_s=0.1),
    )
    with pytest.raises(ValidationError, match="metadata must be a dictionary"):
        validate_correction_result(result)

    result = CorrectionResult(
        boa=mock_observation_bundle.toa,
        boa_unc=None,
        aot=mock_solved_atmosphere.aot,
        tcwv=mock_solved_atmosphere.tcwv,
        cloud_mask=mock_observation_bundle.cloud_mask,
        diagnostics=CorrectionDiagnostics(processing_time_s=np.nan),
    )
    with pytest.raises(ValidationError, match="must be finite"):
        validate_correction_result(result)

    bad_atmo = dataclasses.replace(
        mock_atmospheric_state,
        aot_unc=xr.DataArray(np.array([[np.nan]], dtype=np.float32), dims=["y", "x"]),
    )
    validate_atmospheric_state(bad_atmo)


def test_relative_spectral_response_helpers_cover_validation_and_fwhm_paths(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delattr("siac.domain.spectral.np.trapezoid")
    assert rust_compat is not None  # keep import used for module coverage grouping
    assert RelativeSpectralResponse.__name__ == "RelativeSpectralResponse"

    from siac.domain import spectral as spectral_mod

    assert spectral_mod._trapezoid(np.array([1.0, 1.0]), np.array([0.0, 2.0])) == pytest.approx(2.0)
    assert spectral_mod._fwhm(np.array([500.0]), np.array([1.0])) is None
    assert spectral_mod._fwhm(np.array([500.0, 510.0]), np.array([0.0, 0.0])) is None
    assert spectral_mod._fwhm(np.array([500.0, 510.0]), np.array([-1.0, -2.0])) is None

    canonical = _rsrf()
    assert canonical.centre_wavelength_nm == canonical.center_wavelength_nm

    with pytest.raises(ValueError, match="1-D"):
        RelativeSpectralResponse(
            sensor_id="MSI",
            satellite_id="S2A",
            band_name="B02",
            wavelengths_nm=np.array([[440.0, 450.0]], dtype=np.float32),
            response=np.array([0.4, 0.6], dtype=np.float32),
        )

    with pytest.raises(ValueError, match="equal length"):
        RelativeSpectralResponse(
            sensor_id="MSI",
            satellite_id="S2A",
            band_name="B02",
            wavelengths_nm=np.array([440.0, 450.0], dtype=np.float32),
            response=np.array([1.0], dtype=np.float32),
        )

    with pytest.raises(ValueError, match="at least two samples"):
        RelativeSpectralResponse.from_tabulated(
            sensor_id="MSI",
            satellite_id="S2A",
            band_name="B02",
            wavelengths_nm=np.array([440.0], dtype=np.float32),
            response=np.array([1.0], dtype=np.float32),
        )

    with pytest.raises(ValueError, match="must be finite"):
        RelativeSpectralResponse.from_tabulated(
            sensor_id="MSI",
            satellite_id="S2A",
            band_name="B02",
            wavelengths_nm=np.array([440.0, np.nan], dtype=np.float32),
            response=np.array([1.0, 1.0], dtype=np.float32),
        )

    with pytest.raises(ValueError, match="integrate to a positive value"):
        RelativeSpectralResponse.from_tabulated(
            sensor_id="MSI",
            satellite_id="S2A",
            band_name="B02",
            wavelengths_nm=np.array([440.0, 440.0, 440.0], dtype=np.float32),
            response=np.array([1.0, 1.0, 1.0], dtype=np.float32),
        )


def test_rsrf_kernel_covers_sparse_support_and_error_paths() -> None:
    with pytest.raises(ValueError, match="at least two samples"):
        rsrf_kernel.build_aligned_rsrf_kernel(
            _rsrf(),
            lut_wavelengths_nm=np.array([440.0], dtype=np.float32),
            lut_id="lut-v1",
        )

    with pytest.raises(ValueError, match="support_padding must be non-negative"):
        rsrf_kernel.build_aligned_rsrf_kernel(
            _rsrf(),
            lut_wavelengths_nm=np.array([440.0, 450.0], dtype=np.float32),
            lut_id="lut-v1",
            support_padding=-1,
        )

    with pytest.raises(ValueError, match="zero overlap"):
        rsrf_kernel.build_aligned_rsrf_kernel(
            _rsrf(),
            lut_wavelengths_nm=np.array([100.0, 200.0, 300.0], dtype=np.float32),
            lut_id="lut-v1",
            support_padding=0,
        )

    with pytest.raises(ValueError, match="Solar irradiance must match"):
        rsrf_kernel.build_aligned_rsrf_kernel(
            _rsrf(),
            lut_wavelengths_nm=np.array([440.0, 450.0, 460.0], dtype=np.float32),
            lut_id="lut-v1",
            solar_irradiance=np.array([1.0, 2.0], dtype=np.float32),
        )


def test_rust_compat_proxies_delegate_when_native_symbols_exist(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: dict[str, object] = {}

    class _FakeKernel:
        def __init__(self, *args, **kwargs):  # noqa: ANN002, ANN003
            calls.setdefault("init", []).append((args, kwargs))
            self.extra = "value"

        def compute(self, *args):  # noqa: ANN002
            return ("compute", args)

        def predict(self, *args, **kwargs):  # noqa: ANN002, ANN003
            return ("predict", args, kwargs)

        def convolve(self, *args, **kwargs):  # noqa: ANN002, ANN003
            return ("convolve", args, kwargs)

    fake_native = SimpleNamespace(
        RossThickLiSparse=_FakeKernel,
        TwoLayerNN=_FakeKernel,
        PSFConvolver=_FakeKernel,
        apply_laplacian=lambda *args, **kwargs: ("lap", args, kwargs),
        evaluate_block_grid_search_cost_cube_with_provider_qa=lambda *args, **kwargs: (
            "block_cube_qa",
            args,
            kwargs,
        ),
        evaluate_grid_search_candidate_cost=lambda *args, **kwargs: ("cand", args, kwargs),
        evaluate_grid_search_cost_cube_with_provider=lambda *args, **kwargs: ("cube", args, kwargs),
        evaluate_grid_search_cost_cube_with_provider_qa=lambda *args, **kwargs: ("cube_qa", args, kwargs),
        interpolate_to_fine_grid=lambda *args, **kwargs: ("fine", args, kwargs),
        quadratic_refine_grid_search=lambda *args, **kwargs: ("quad", args, kwargs),
        quadratic_refine_grid_search_qa=lambda *args, **kwargs: ("quad_qa", args, kwargs),
        remap_to_coarse_grid=lambda *args, **kwargs: ("coarse", args, kwargs),
        whittaker_smooth_cube=lambda *args, **kwargs: ("smooth", args, kwargs),
    )
    monkeypatch.setattr(rust_compat, "_NATIVE_RUST", fake_native)
    monkeypatch.setattr(rust_compat, "_NATIVE_IMPORT_ERROR", None)

    brdf = rust_compat.RossThickLiSparse(2.0, 1.0)
    nn = rust_compat.TwoLayerNN("weights")
    psf = rust_compat.PSFConvolver("kernel")

    assert brdf.compute(1, 2, 3) == ("compute", (1, 2, 3))
    assert nn.predict("x", scale=2) == ("predict", ("x",), {"scale": 2})
    assert psf.convolve("img") == ("convolve", ("img",), {})
    assert brdf.extra == "value"
    assert nn.extra == "value"
    assert psf.extra == "value"

    assert rust_compat.apply_laplacian(1) == ("lap", (1,), {})
    assert rust_compat.evaluate_block_grid_search_cost_cube_with_provider_qa(2) == (
        "block_cube_qa",
        (2,),
        {},
    )
    assert rust_compat.evaluate_grid_search_candidate_cost(2) == ("cand", (2,), {})
    assert rust_compat.evaluate_grid_search_cost_cube_with_provider(3) == ("cube", (3,), {})
    assert rust_compat.evaluate_grid_search_cost_cube_with_provider_qa(3) == ("cube_qa", (3,), {})
    assert rust_compat.interpolate_to_fine_grid(4) == ("fine", (4,), {})
    assert rust_compat.quadratic_refine_grid_search(5) == ("quad", (5,), {})
    assert rust_compat.quadratic_refine_grid_search_qa(5) == ("quad_qa", (5,), {})
    assert rust_compat.remap_to_coarse_grid(6) == ("coarse", (6,), {})
    assert rust_compat.whittaker_smooth_cube(7) == ("smooth", (7,), {})


def test_rust_compat_missing_message_includes_original_import_error(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(rust_compat, "_NATIVE_IMPORT_ERROR", ModuleNotFoundError("boom"))
    message = rust_compat._missing_message("PSFConvolver")
    assert "PSFConvolver" in message
    assert "boom" in message


def test_rust_compat_missing_symbol_raises_import_error(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(rust_compat, "_NATIVE_RUST", SimpleNamespace())
    monkeypatch.setattr(rust_compat, "_NATIVE_IMPORT_ERROR", None)
    with pytest.raises(ImportError, match="apply_laplacian"):
        rust_compat.apply_laplacian()


def test_validation_solver_band_error_uses_real_sensor_band(mock_solver_input_bundle) -> None:
    bad_band = SensorBand("FAKE", 9999.0, 10.0, 10.0, 99)
    bad_sib = dataclasses.replace(mock_solver_input_bundle, bands=[bad_band])
    with pytest.raises(ValidationError, match="not in sensor_config"):
        validate_solver_input_bundle(bad_sib)
