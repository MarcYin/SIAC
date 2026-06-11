#!/usr/bin/env python3
"""Compare SIAC LUT and emulator RT backends under identical inputs."""

from __future__ import annotations

import argparse
import json
import logging
from collections import OrderedDict
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

from siac.adapters.rt import build_rt_model
from siac.algorithms.rt.emulator.two_nn import _BandEmulator
from siac.app.planning import build_execution_plan
from siac.app.requests import SceneProcessRequest
from siac.config import SIACConfig
from siac.runtime import AtmosphericState, GeometryAngles, RTCoefficients, SurfacePrior
from siac.workflows.scene_setup import aerosol_resolution, call_grid_assembler, select_band_slice

LOGGER = logging.getLogger("compare_rt_backends")


class LegacySplitEmulator:
    """Adapter for legacy master-branch emulator weights split by coefficient."""

    _SATELLITE_FALLBACKS: dict[str, str] = {"S2C": "S2A"}

    def __init__(self, emulator_dir: Path, *, sensor_id: str, satellite_id: str) -> None:
        self.emulator_dir = emulator_dir
        self.sensor_id = sensor_id
        self.satellite_id = satellite_id
        self._cache: OrderedDict[tuple[str, str], _BandEmulator] = OrderedDict()
        self._max_cache = 32

    def _candidate_satellites(self) -> tuple[str, ...]:
        fallback = self._SATELLITE_FALLBACKS.get(self.satellite_id)
        return tuple(sat for sat in (self.satellite_id, fallback) if sat is not None)

    def _find_path(self, band_name: str, coeff_name: str) -> Path:
        for sat_id in self._candidate_satellites():
            patterns = (
                f"*{self.sensor_id}*_{sat_id}_{band_name}_{coeff_name}.npz",
                f"*{sat_id}_{band_name}_{coeff_name}.npz",
                f"*{sat_id}*_{band_name}_{coeff_name}.npz",
            )
            for pattern in patterns:
                matches = sorted(self.emulator_dir.glob(pattern))
                if matches:
                    return matches[0]
        raise FileNotFoundError(
            f"No legacy emulator found for {self.sensor_id}/{self.satellite_id} "
            f"band {band_name} coeff {coeff_name} in {self.emulator_dir}"
        )

    def _load_head(self, band_name: str, coeff_name: str) -> _BandEmulator:
        key = (band_name, coeff_name)
        if key in self._cache:
            self._cache.move_to_end(key)
            return self._cache[key]

        path = self._find_path(band_name, coeff_name)
        data = np.load(path, allow_pickle=True)  # noqa: S301
        hidden_layers = data["Hidden_Layers"].tolist()
        output_layers = data["Output_Layers"].tolist()
        for layer in hidden_layers:
            layer[1] = np.atleast_1d(layer[1]).ravel().astype(np.float32)
        for layer in output_layers:
            layer[1] = np.atleast_1d(layer[1]).ravel().astype(np.float32)

        emulator = _BandEmulator(hidden_layers=hidden_layers, output_layers=output_layers)
        self._cache[key] = emulator
        if len(self._cache) > self._max_cache:
            self._cache.popitem(last=False)
        return emulator

    @staticmethod
    def _prepare_inputs(
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
    ) -> tuple[np.ndarray, tuple[int, ...], xr.DataArray]:
        original_shape = geometry.sza.shape
        inputs = np.column_stack(
            [
                np.cos(geometry.sza.values).ravel(),
                np.cos(geometry.vza.values).ravel(),
                np.cos(geometry.raa.values).ravel(),
                atmo_state.aot.values.ravel(),
                atmo_state.tcwv.values.ravel(),
                atmo_state.tco3.values.ravel(),
                atmo_state.elevation.values.ravel(),
            ]
        ).astype(np.float32)
        return inputs, original_shape, geometry.sza

    def compute_coefficients(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        band: Any,
        compute_jacobian: bool = False,
    ) -> RTCoefficients:
        if compute_jacobian:
            raise NotImplementedError("Legacy split emulator adapter does not implement Jacobians.")

        inputs, original_shape, template = self._prepare_inputs(geometry, atmo_state)
        outputs = []
        for coeff_name in ("xap", "xbp", "xcp"):
            emulator = self._load_head(band.name, coeff_name)
            coeff_output, _ = emulator.forward(inputs, compute_jacobian=False)
            outputs.append(coeff_output[:, 0].reshape(original_shape).astype(np.float32))

        return RTCoefficients(
            xap=xr.DataArray(outputs[0], dims=template.dims, coords=template.coords),
            xbp=xr.DataArray(outputs[1], dims=template.dims, coords=template.coords),
            xcp=xr.DataArray(outputs[2], dims=template.dims, coords=template.coords),
        )


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare LUT and emulator RT backends on the same scene with the "
            "same atmosphere, geometry, and surface prior."
        )
    )
    parser.add_argument("--config", type=Path, required=True, help="SIAC config TOML.")
    parser.add_argument("--input", type=Path, required=True, help="Scene SAFE path.")
    parser.add_argument(
        "--reference-output-dir",
        type=Path,
        required=True,
        help="Existing SIAC output directory providing surface_prior.nc and auxiliary.nc.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where comparison outputs will be written.",
    )
    parser.add_argument(
        "--emulator-dir",
        type=Path,
        required=True,
        help="Directory containing the RT emulator weight files.",
    )
    parser.add_argument(
        "--sensor",
        default="s2",
        help="Sensor override for config resolution. Defaults to s2.",
    )
    parser.add_argument(
        "--elevation-source",
        default="provider",
        choices=("provider", "dem"),
        help=(
            "Which elevation to feed into the RT backends. "
            "'provider' uses the atmospheric provider field as-is; "
            "'dem' reprojects config.paths.dem onto the solver grid and converts meters to km."
        ),
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        help="Logging verbosity.",
    )
    return parser.parse_args()


def _setup_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )


def _surface_prior_from_reference(path: Path) -> SurfacePrior:
    ds = xr.load_dataset(path)
    boa = ds.to_array(dim="band").astype(np.float32)
    mask = xr.DataArray(
        np.all(np.isfinite(boa.values), axis=0),
        dims=("y", "x"),
        coords={"y": boa.coords["y"], "x": boa.coords["x"]},
    )
    boa_unc = xr.zeros_like(boa, dtype=np.float32)
    return SurfacePrior(boa=boa, boa_unc=boa_unc, kernels=None, mask=mask)


def resample_field_to_template(field: xr.DataArray, template: xr.DataArray) -> xr.DataArray:
    if tuple(field.shape) == tuple(template.shape) and field.dims == template.dims:
        coords_match = all(
            dim in field.coords
            and dim in template.coords
            and np.array_equal(
                np.asarray(field.coords[dim].values), np.asarray(template.coords[dim].values)
            )
            for dim in template.dims
        )
        if coords_match:
            return field.astype(np.float32)
    if (
        len(field.dims) == 2
        and field.dims == template.dims
        and all(dim in field.coords for dim in field.dims)
        and all(dim in template.coords for dim in template.dims)
    ):
        try:
            return field.interp(
                coords={dim: template.coords[dim] for dim in template.dims},
                method="linear",
            ).astype(np.float32)
        except Exception:
            LOGGER.debug(
                "Interpolation failed for field %s; falling back to scipy zoom.",
                getattr(field, "name", "?"),
                exc_info=True,
            )

    src = np.asarray(field.values, dtype=np.float32)
    if src.ndim != 2 or len(template.dims) != 2:
        raise ValueError(
            f"Cannot resample field with dims {field.dims} to template dims {template.dims}"
        )

    from scipy import ndimage

    h_out = int(template.sizes[template.dims[0]])
    w_out = int(template.sizes[template.dims[1]])
    if src.shape[0] == 0 or src.shape[1] == 0:
        out: np.ndarray = np.full((h_out, w_out), np.nan, dtype=np.float32)
    else:
        out = ndimage.zoom(src, (h_out / src.shape[0], w_out / src.shape[1]), order=1)
        out = out[:h_out, :w_out]
        if out.shape != (h_out, w_out):
            padded = np.full((h_out, w_out), np.nan, dtype=np.float32)
            padded[: out.shape[0], : out.shape[1]] = out
            out = padded

    return xr.DataArray(
        out.astype(np.float32),
        dims=template.dims,
        coords={dim: template.coords[dim] for dim in template.dims if dim in template.coords},
        attrs=field.attrs,
    )


def _reconstruct_run_atmosphere(
    reference_output_dir: Path, atmo_prior: AtmosphericState
) -> AtmosphericState:
    aux = xr.load_dataset(reference_output_dir / "auxiliary.nc")
    template = aux["aot"].astype(np.float32)
    prior_aot = resample_field_to_template(atmo_prior.aot, template)
    prior_tcwv = resample_field_to_template(atmo_prior.tcwv, template)
    return AtmosphericState(
        aot=xr.where(np.isfinite(template), template, prior_aot).astype(np.float32),
        tcwv=xr.where(np.isfinite(aux["tcwv"]), aux["tcwv"], prior_tcwv).astype(np.float32),
        tco3=resample_field_to_template(atmo_prior.tco3, template).astype(np.float32),
        aot_unc=resample_field_to_template(atmo_prior.aot_unc, template).astype(np.float32),
        tcwv_unc=resample_field_to_template(atmo_prior.tcwv_unc, template).astype(np.float32),
        tco3_unc=resample_field_to_template(atmo_prior.tco3_unc, template).astype(np.float32),
        elevation=resample_field_to_template(atmo_prior.elevation, template).astype(np.float32),
    )


def _load_dem_elevation(template: xr.DataArray, dem_path: str | Path) -> xr.DataArray:
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.vrt import WarpedVRT

    width = int(template.sizes[template.dims[-1]])
    height = int(template.sizes[template.dims[-2]])
    target_crs = template.rio.crs
    if target_crs is None:
        raise ValueError("Template grid does not have a CRS; cannot project DEM elevation.")

    with (
        rasterio.open(str(dem_path)) as src,
        WarpedVRT(
            src,
            crs=str(target_crs),
            transform=template.rio.transform(recalc=True),
            width=width,
            height=height,
            resampling=Resampling.bilinear,
            src_nodata=src.nodata,
            nodata=np.nan,
        ) as vrt,
    ):
        values = vrt.read(1, out_dtype="float32", masked=True).filled(np.nan)

    # Match the master branch convention: DEM is converted from meters to km and gaps are filled with 0.
    values_km = np.where(np.isfinite(values), values * np.float32(0.001), np.float32(0.0))
    return template.copy(data=values_km.astype(np.float32))


def _with_elevation(atmo_state: AtmosphericState, elevation: xr.DataArray) -> AtmosphericState:
    return AtmosphericState(
        aot=atmo_state.aot.astype(np.float32),
        tcwv=atmo_state.tcwv.astype(np.float32),
        tco3=atmo_state.tco3.astype(np.float32),
        aot_unc=atmo_state.aot_unc.astype(np.float32),
        tcwv_unc=atmo_state.tcwv_unc.astype(np.float32),
        tco3_unc=atmo_state.tco3_unc.astype(np.float32),
        elevation=elevation.astype(np.float32),
    )


def _corrcoef(a: np.ndarray, b: np.ndarray) -> float:
    if a.size < 2 or b.size < 2:
        return float("nan")
    if np.allclose(a, a[0]) or np.allclose(b, b[0]):
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _pair_metrics(left: np.ndarray, right: np.ndarray) -> dict[str, Any]:
    delta = right - left
    return {
        "left_mean": float(np.mean(left)),
        "right_mean": float(np.mean(right)),
        "bias_right_minus_left": float(np.mean(delta)),
        "mae_right_minus_left": float(np.mean(np.abs(delta))),
        "rmse_right_minus_left": float(np.sqrt(np.mean(delta**2))),
        "corr": _corrcoef(left, right),
    }


def _coeff_metrics(
    left: RTCoefficients,
    right: RTCoefficients,
    valid: np.ndarray,
) -> tuple[dict[str, Any], xr.Dataset]:
    metrics: dict[str, Any] = {}
    ds = xr.Dataset()
    for name in ("xap", "xbp", "xcp"):
        left_field = getattr(left, name).astype(np.float32)
        right_field = getattr(right, name).astype(np.float32)
        band_valid = valid & np.isfinite(left_field.values) & np.isfinite(right_field.values)
        values_left = np.asarray(left_field.values[band_valid], dtype=np.float64)
        values_right = np.asarray(right_field.values[band_valid], dtype=np.float64)
        metrics[name] = {
            "valid_count": int(values_left.size),
            **(_pair_metrics(values_left, values_right) if values_left.size else {}),
        }
        ds[f"lut_{name}"] = left_field
        ds[f"emulator_{name}"] = right_field
        ds[f"emulator_minus_lut_{name}"] = (right_field - left_field).astype(np.float32)
    return metrics, ds


def _toa_metrics(
    observed: xr.DataArray,
    simulated_lut: xr.DataArray,
    simulated_emulator: xr.DataArray,
    valid: np.ndarray,
) -> dict[str, Any]:
    lut_valid = valid & np.isfinite(observed.values) & np.isfinite(simulated_lut.values)
    emu_valid = valid & np.isfinite(observed.values) & np.isfinite(simulated_emulator.values)
    shared_valid = lut_valid & emu_valid
    obs = np.asarray(observed.values[shared_valid], dtype=np.float64)
    lut_values = np.asarray(simulated_lut.values[shared_valid], dtype=np.float64)
    emu_values = np.asarray(simulated_emulator.values[shared_valid], dtype=np.float64)
    if obs.size == 0:
        return {"valid_count": 0}
    lut_resid = lut_values - obs
    emu_resid = emu_values - obs
    return {
        "valid_count": int(obs.size),
        "lut": {
            "simulated_mean": float(np.mean(lut_values)),
            "bias_simulated_minus_observed": float(np.mean(lut_resid)),
            "rmse_simulated_minus_observed": float(np.sqrt(np.mean(lut_resid**2))),
            "corr": _corrcoef(obs, lut_values),
        },
        "emulator": {
            "simulated_mean": float(np.mean(emu_values)),
            "bias_simulated_minus_observed": float(np.mean(emu_resid)),
            "rmse_simulated_minus_observed": float(np.sqrt(np.mean(emu_resid**2))),
            "corr": _corrcoef(obs, emu_values),
        },
        "emulator_minus_lut": {
            "bias": float(np.mean(emu_values - lut_values)),
            "rmse": float(np.sqrt(np.mean((emu_values - lut_values) ** 2))),
            "corr": _corrcoef(lut_values, emu_values),
        },
    }


def _geometry_on_template(geometry: GeometryAngles, template: xr.DataArray) -> GeometryAngles:
    return GeometryAngles(
        sza=resample_field_to_template(geometry.sza, template),
        saa=resample_field_to_template(geometry.saa, template),
        vza=resample_field_to_template(geometry.vza, template),
        vaa=resample_field_to_template(geometry.vaa, template),
    )


def _build_emulator_backend(config: Any, plan: Any, obs: Any, emulator_dir: Path) -> Any:
    if (
        any(emulator_dir.glob("*_xap.npz"))
        and any(emulator_dir.glob("*_xbp.npz"))
        and any(emulator_dir.glob("*_xcp.npz"))
    ):
        LOGGER.info("Using legacy split-format emulator adapter for %s", emulator_dir)
        return LegacySplitEmulator(
            emulator_dir=emulator_dir,
            sensor_id=obs.sensor_config.sensor_id,
            satellite_id=obs.sensor_config.satellite_id,
        )

    emulator_config = config.model_copy(deep=True)
    emulator_config.rt_model.backend = "emulator"
    emulator_config.rt_model.fallback_to_lut = False
    emulator_config.paths.emulator_dir = emulator_dir
    try:
        return build_rt_model(emulator_config, auth=plan.auth, sensor_config=obs.sensor_config)
    except Exception as exc:
        LOGGER.warning(
            "Native emulator loader unavailable; trying legacy split-format adapter (%s)", exc
        )
        return LegacySplitEmulator(
            emulator_dir=emulator_dir,
            sensor_id=obs.sensor_config.sensor_id,
            satellite_id=obs.sensor_config.satellite_id,
        )


def main() -> None:
    args = _parse_args()
    _setup_logging(args.log_level)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    if not args.emulator_dir.exists():
        raise FileNotFoundError(f"Emulator directory does not exist: {args.emulator_dir}")

    config = SIACConfig.from_file(args.config).with_overrides(sensor=args.sensor)
    request = SceneProcessRequest(config=config, input_path=args.input)
    plan = build_execution_plan(request)

    LOGGER.info("Preprocessing observation from %s", args.input)
    obs = plan.preprocessor(plan.input_path, plan.runtime_aoi)
    aerosol_resolution = aerosol_resolution(plan.config)

    LOGGER.info("Fetching atmospheric prior")
    atmo_prior = plan.atmo_provider(
        obs.bounds, obs.crs, obs.metadata["observation_time"], aerosol_resolution
    )
    run_atmo = _reconstruct_run_atmosphere(args.reference_output_dir, atmo_prior)
    surface_prior = _surface_prior_from_reference(args.reference_output_dir / "surface_prior.nc")

    LOGGER.info("Building LUT and emulator RT backends")
    lut_config = config.model_copy(deep=True)
    lut_config.rt_model.backend = "lut"

    lut_rt = build_rt_model(lut_config, auth=plan.auth, sensor_config=obs.sensor_config)
    emulator_rt = _build_emulator_backend(config, plan, obs, args.emulator_dir)

    summary: dict[str, Any] = {
        "config_path": str(args.config.resolve()),
        "input_path": str(args.input.resolve()),
        "reference_output_dir": str(args.reference_output_dir.resolve()),
        "emulator_dir": str(args.emulator_dir.resolve()),
        "emulator_backend": type(emulator_rt).__name__,
        "elevation_source": args.elevation_source,
        "observation_time": obs.metadata["observation_time"].isoformat(),
        "aerosol_resolution_m": float(aerosol_resolution),
        "scenarios": {},
    }

    for scenario_name, atmosphere in {"atmo_prior": atmo_prior, "run_solution": run_atmo}.items():
        LOGGER.info("Assembling solver grid for scenario %s", scenario_name)
        solver_inputs = call_grid_assembler(
            plan.grid_assembler,
            obs,
            atmosphere,
            surface_prior,
            lut_rt,
            aerosol_resolution_m=float(aerosol_resolution),
        )
        valid_mask = solver_inputs.cloud_mask.values.astype(bool) == 0
        rt_atmo = solver_inputs.atmo_prior
        if args.elevation_source == "dem":
            dem_path = getattr(config.paths, "dem", None)
            if dem_path is None:
                raise ValueError(
                    "elevation_source='dem' requires config.paths.dem to be configured."
                )
            dem_elevation = _load_dem_elevation(solver_inputs.atmo_prior.aot, dem_path)
            rt_atmo = _with_elevation(solver_inputs.atmo_prior, dem_elevation)
        geometry = _geometry_on_template(solver_inputs.geometry, rt_atmo.aot)

        elevation_values = np.asarray(rt_atmo.elevation.values, dtype=np.float64)
        finite_elevation = elevation_values[np.isfinite(elevation_values)]
        scenario_summary: dict[str, Any] = {
            "elevation_km": {
                "min": float(np.min(finite_elevation)) if finite_elevation.size else float("nan"),
                "max": float(np.max(finite_elevation)) if finite_elevation.size else float("nan"),
                "mean": float(np.mean(finite_elevation)) if finite_elevation.size else float("nan"),
            }
        }
        coeff_ds = xr.Dataset()
        toa_ds = xr.Dataset()
        coeff_ds["elevation_km"] = rt_atmo.elevation.astype(np.float32)

        for band_index, band in enumerate(solver_inputs.bands):
            band_name = band.name
            observed = select_band_slice(
                solver_inputs.toa, band_name=band_name, band_index=band_index
            )
            surface = select_band_slice(
                solver_inputs.surface_prior.boa,
                band_name=band_name,
                band_index=band_index,
            )
            if observed is None or surface is None:
                continue
            lut_coeffs = lut_rt.compute_coefficients(geometry, rt_atmo, band, False)
            emulator_coeffs = emulator_rt.compute_coefficients(geometry, rt_atmo, band, False)

            band_valid = valid_mask & np.isfinite(observed.values) & np.isfinite(surface.values)
            coeff_metrics, band_coeff_ds = _coeff_metrics(lut_coeffs, emulator_coeffs, band_valid)
            simulated_lut = lut_coeffs.simulate_toa(surface).astype(np.float32)
            simulated_emulator = emulator_coeffs.simulate_toa(surface).astype(np.float32)
            toa_metrics = _toa_metrics(observed, simulated_lut, simulated_emulator, band_valid)

            scenario_summary[band_name] = {
                "coefficients": coeff_metrics,
                "toa": toa_metrics,
            }

            coeff_ds.update(
                {f"{band_name}_{name}": band_coeff_ds[name] for name in band_coeff_ds.data_vars}
            )
            toa_ds[f"observed_toa_{band_name}"] = observed.astype(np.float32)
            toa_ds[f"simulated_toa_lut_{band_name}"] = simulated_lut
            toa_ds[f"simulated_toa_emulator_{band_name}"] = simulated_emulator
            toa_ds[f"simulated_toa_emulator_minus_lut_{band_name}"] = (
                simulated_emulator - simulated_lut
            ).astype(np.float32)
            toa_ds[f"residual_toa_lut_{band_name}"] = (simulated_lut - observed).astype(np.float32)
            toa_ds[f"residual_toa_emulator_{band_name}"] = (simulated_emulator - observed).astype(
                np.float32
            )

        summary["scenarios"][scenario_name] = scenario_summary
        coeff_ds.to_netcdf(args.output_dir / f"coefficients_{scenario_name}.nc")
        toa_ds.to_netcdf(args.output_dir / f"toa_{scenario_name}.nc")

    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    LOGGER.info("RT backend comparison written to %s", args.output_dir)


if __name__ == "__main__":
    main()
