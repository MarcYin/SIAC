"""Analytic inversion of Sen2Cor's rugged-terrain correction.

Sen2Cor normalizes slope reflectance by the local irradiance (direct beam on
the tilted surface + sky-view-weighted diffuse, with an empirical shade
limiter); our RT frame is flat-terrain. ``deterrain_factor`` converts L2A back
to the flat frame with zero fitted parameters: the per-band direct-beam
fraction comes from the libRadtran LUT's ``Eg_dir_rho1``/``Eg_rho1`` irradiance
spectra at the scene state Sen2Cor itself used, and the limiter from the
Sen2Cor/ATCOR configuration.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from pathlib import Path

#: Sen2Cor/ATCOR empirical-BRDF lower bound on the illumination ratio.
G_SHADE_LIMITER = 0.25
#: ATCOR empirical-BRDF exponent: full correction below 720 nm, 3/4 beyond.
BRDF_EXPONENT_VIS = 1.0
BRDF_EXPONENT_NIR_SWIR = 0.75
NIR_SWIR_BANDS = ("nir08", "swir16", "swir22")


class DirectFractionLUT:
    """Band-effective direct-beam irradiance fraction from the spectral LUT."""

    def __init__(self, lut_path: str | Path) -> None:
        import zarr

        self._store = zarr.storage.ZipStore(str(lut_path), mode="r")
        root = zarr.open(self._store, mode="r")
        self._eg = root["Eg_rho1"]
        self._eg_dir = root["Eg_dir_rho1"]
        self._axes = {
            name: np.asarray(root[name][:], dtype=np.float64)
            for name in ("wavelength", "sza", "aot", "tcwv", "ozone", "altitude")
        }
        sizes = {name: values.size for name, values in self._axes.items()}
        if len(set(sizes.values())) != len(sizes):
            raise ValueError("LUT axis sizes are not unique; cannot infer dimension order")
        by_size = {size: name for name, size in sizes.items()}
        self._dim_names = [by_size[dim] for dim in self._eg.shape]
        if self._dim_names[0] != "wavelength":
            raise ValueError(f"unexpected LUT layout {self._dim_names}")
        self._spectra_cache: dict[tuple[int, ...], tuple[np.ndarray, np.ndarray]] = {}

    def _spectra(self, indices: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray]:
        cached = self._spectra_cache.get(indices)
        if cached is None:
            selector = (slice(None), *indices)
            cached = (
                np.asarray(self._eg_dir[selector], dtype=np.float64),
                np.asarray(self._eg[selector], dtype=np.float64),
            )
            self._spectra_cache[indices] = cached
        return cached

    def band_direct_fractions(
        self,
        bands: list,
        *,
        aot: float,
        sza_deg: float,
        tcwv_cm: float,
        tco3_atm_cm: float,
        elevation_km: float,
    ) -> dict[str, float]:
        """RSRF-weighted E_dir/E_g per band, bilinear in (sza, aot)."""
        state = {
            "sza": float(sza_deg),
            "aot": float(np.clip(aot, self._axes["aot"][0], self._axes["aot"][-1])),
            "tcwv": float(tcwv_cm) * 10.0 if self._axes["tcwv"].max() > 20 else float(tcwv_cm),
            "ozone": float(tco3_atm_cm) * 1000.0
            if self._axes["ozone"].max() > 10
            else float(tco3_atm_cm),
            "altitude": float(elevation_km),
        }
        nearest = {
            name: int(np.argmin(np.abs(self._axes[name] - value))) for name, value in state.items()
        }

        def bracket(name: str) -> tuple[int, int, float]:
            axis = self._axes[name]
            value = np.clip(state[name], axis[0], axis[-1])
            upper = int(np.searchsorted(axis, value))
            upper = min(max(upper, 1), axis.size - 1)
            lower = upper - 1
            weight = float((value - axis[lower]) / (axis[upper] - axis[lower]))
            return lower, upper, weight

        sza_lo, sza_hi, sza_w = bracket("sza")
        aot_lo, aot_hi, aot_w = bracket("aot")
        interp_dims = list(self._dim_names[1:])
        direct = np.zeros_like(self._axes["wavelength"])
        total = np.zeros_like(self._axes["wavelength"])
        for sza_index, sza_weight in ((sza_lo, 1.0 - sza_w), (sza_hi, sza_w)):
            for aot_index, aot_weight in ((aot_lo, 1.0 - aot_w), (aot_hi, aot_w)):
                point = {**nearest, "sza": sza_index, "aot": aot_index}
                indices = tuple(point[name] for name in interp_dims)
                spectrum_dir, spectrum_total = self._spectra(indices)
                weight = sza_weight * aot_weight
                direct += weight * spectrum_dir
                total += weight * spectrum_total

        fractions: dict[str, float] = {}
        wavelength = self._axes["wavelength"]
        for band in bands:
            if band.rsrf_wavelengths_nm is not None and band.rsrf_response is not None:
                weights = np.interp(
                    wavelength,
                    np.asarray(band.rsrf_wavelengths_nm, dtype=np.float64),
                    np.asarray(band.rsrf_response, dtype=np.float64),
                    left=0.0,
                    right=0.0,
                )
            else:
                half = float(band.bandwidth) / 2.0
                weights = (
                    (wavelength >= band.center_wavelength - half)
                    & (wavelength <= band.center_wavelength + half)
                ).astype(np.float64)
            denominator = float(np.sum(weights * total))
            if denominator <= 0.0:
                raise ValueError(f"band {band.name} has no LUT spectral support")
            fractions[band.name] = float(np.sum(weights * direct) / denominator)
        return fractions


def deterrain_factor(
    c_direct: float,
    incidence_cos: np.ndarray,
    slope_deg: np.ndarray,
    cos_sza: float,
    g_limiter: float = G_SHADE_LIMITER,
    rho_background: float = 0.0,
    exponent: float = 1.0,
) -> np.ndarray:
    """Multiplier converting Sen2Cor terrain-corrected reflectance to flat frame.

    ``rho_flat = rho_L2A * factor``. On flat terrain the factor is exactly 1.
    ``rho_background`` adds the ATCOR terrain-reflected irradiance term
    ``rho * (1 - V_sky)`` — slopes are also lit by their (bright) surroundings,
    which matters most in the SWIR over deserts.
    """
    sky_view = 0.5 * (1.0 + np.cos(np.radians(np.clip(slope_deg, 0.0, 90.0))))
    effective = np.maximum(np.asarray(incidence_cos, dtype=np.float64), g_limiter * cos_sza)
    flat_irradiance = c_direct * cos_sza + (1.0 - c_direct)
    terrain_reflected = float(np.clip(rho_background, 0.0, 1.0)) * (1.0 - sky_view)
    numerator = (
        c_direct * effective + (1.0 - c_direct) * sky_view + terrain_reflected * flat_irradiance
    )
    ratio = np.clip(numerator / flat_irradiance, 0.2, 3.0)
    return np.power(ratio, float(exponent))


def deterrain_l2a_surface(
    surface: np.ndarray,
    *,
    l2a_aot: np.ndarray,
    l2a_tcwv: np.ndarray,
    terrain,
    sza_deg: float,
    saa_deg: float,
    elevation_km: float,
    satellite_id: str,
    fraction_lut: DirectFractionLUT,
    sensor_cache: dict,
    band_names: tuple[str, ...],
    s2_band_map: dict[str, str],
    tco3_atm_cm: float = 0.30,
) -> tuple[np.ndarray, dict]:
    """Convert a terrain-corrected L2A grid stack back to the flat-RT frame.

    ``surface`` is (band, y, x) in canonical band order; NaNs pass through.
    The Sen2Cor state is taken from its own per-pixel AOT/WVP rasters, per the
    mimic-the-processor rule.
    """
    from tools.aeronet_validation.terrain_features import local_solar_incidence

    from siac.adapters.rsrf import load_sensor_config_with_rsrf

    sensor = sensor_cache.get(satellite_id)
    if sensor is None:
        sensor = load_sensor_config_with_rsrf("MSI", satellite_id)
        sensor_cache[satellite_id] = sensor
    bands = [sensor.get_band(s2_band_map[name]) for name in band_names]
    state = {
        "aot": float(np.nanmedian(l2a_aot)),
        "tcwv_cm": float(np.nanmedian(l2a_tcwv)),
        "tco3_atm_cm": float(tco3_atm_cm),
    }
    fractions = fraction_lut.band_direct_fractions(
        bands,
        aot=state["aot"],
        sza_deg=float(sza_deg),
        tcwv_cm=state["tcwv_cm"],
        tco3_atm_cm=state["tco3_atm_cm"],
        elevation_km=float(elevation_km),
    )
    incidence = local_solar_incidence(terrain, sza_deg=float(sza_deg), saa_deg=float(saa_deg))
    cos_sza = float(np.cos(np.radians(float(sza_deg))))
    result = np.asarray(surface, dtype=np.float32).copy()
    c_by_band: dict[str, float] = {}
    rho_by_band: dict[str, float] = {}
    for index, name in enumerate(band_names):
        c_direct = fractions[s2_band_map[name]]
        rho_background = float(np.nanmedian(result[index]))
        if not np.isfinite(rho_background):
            rho_background = 0.0
        exponent = BRDF_EXPONENT_NIR_SWIR if name in NIR_SWIR_BANDS else BRDF_EXPONENT_VIS
        factor = deterrain_factor(
            c_direct,
            incidence,
            terrain.slope_deg,
            cos_sza,
            rho_background=rho_background,
            exponent=exponent,
        )
        result[index] = (result[index] * factor).astype(np.float32)
        c_by_band[name] = round(float(c_direct), 4)
        rho_by_band[name] = round(rho_background, 4)
    provenance = {
        "applied": True,
        "g_shade_limiter": G_SHADE_LIMITER,
        "sen2cor_state": {key: round(value, 4) for key, value in state.items()},
        "c_by_band": c_by_band,
        "rho_background_by_band": rho_by_band,
    }
    return result, provenance
