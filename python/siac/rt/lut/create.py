"""Utilities for generating LUT Zarr stores from Py6S."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import xarray as xr

logger = logging.getLogger(__name__)


def create_lut_from_py6s(
    output_path: str | Path,
    wavelengths: list[float],
    sza_range: tuple[float, float, float] = (0, 75, 5),
    vza_range: tuple[float, float, float] = (0, 45, 5),
    raa_range: tuple[float, float, float] = (0, 180, 15),
    aot_values: list[float] | None = None,
    tcwv_values: list[float] | None = None,
    aerosol_type: str = "continental",
    ozone: float = 0.3,
) -> None:
    """
    Create a Zarr LUT from Py6S simulations.

    This utility function generates the pre-computed LUT by running
    6S simulations for the specified parameter ranges.

    Args:
        output_path: Path for output Zarr store
        wavelengths: List of wavelengths to simulate (nm)
        sza_range: (min, max, step) for solar zenith angle
        vza_range: (min, max, step) for view zenith angle
        raa_range: (min, max, step) for relative azimuth angle
        aot_values: AOT values to simulate (default: log-spaced 0.01-2.5)
        tcwv_values: TCWV values to simulate (default: 0.5-6.0)
        aerosol_type: Aerosol model type
        ozone: Fixed ozone value (atm-cm)
    """
    try:
        from Py6S import (
            AeroProfile,
            AtmosProfile,
            Geometry,
            SixS,
            Wavelength,
        )
    except ImportError as exc:
        raise ImportError(
            "Py6S is required to create LUT. Install with: pip install Py6S"
        ) from exc

    output_path = Path(output_path)

    # Default parameter values
    if aot_values is None:
        aot_values = np.logspace(-2, np.log10(2.5), 15).tolist()
    if tcwv_values is None:
        tcwv_values = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]

    # Create coordinate arrays
    sza_arr = np.arange(*sza_range)
    vza_arr = np.arange(*vza_range)
    raa_arr = np.arange(*raa_range)
    aot_arr = np.array(aot_values)
    tcwv_arr = np.array(tcwv_values)
    wl_arr = np.array(wavelengths)

    # Initialize output arrays
    shape = (len(sza_arr), len(vza_arr), len(raa_arr), len(aot_arr), len(tcwv_arr), len(wl_arr))
    path_ref = np.zeros(shape, dtype=np.float32)
    trans_down = np.zeros(shape, dtype=np.float32)
    trans_up = np.zeros(shape, dtype=np.float32)
    sph_alb = np.zeros(shape, dtype=np.float32)

    # Initialize 6S
    s = SixS()
    s.atmos_profile = AtmosProfile.UserWaterAndOzone(tcwv_arr[0], ozone)

    if aerosol_type == "continental":
        s.aero_profile = AeroProfile.Continental
    elif aerosol_type == "maritime":
        s.aero_profile = AeroProfile.Maritime
    else:
        s.aero_profile = AeroProfile.Continental

    total_sims = np.prod(shape)
    logger.info(f"Running {total_sims} 6S simulations...")

    # Run simulations
    count = 0
    for i_sza, sza in enumerate(sza_arr):
        for i_vza, vza in enumerate(vza_arr):
            for i_raa, raa in enumerate(raa_arr):
                for i_aot, aot in enumerate(aot_arr):
                    for i_tcwv, tcwv in enumerate(tcwv_arr):
                        for i_wl, wl in enumerate(wl_arr):
                            # Configure 6S
                            s.geometry = Geometry.User()
                            s.geometry.solar_z = sza
                            s.geometry.solar_a = 0
                            s.geometry.view_z = vza
                            s.geometry.view_a = raa
                            s.geometry.day = 1
                            s.geometry.month = 6

                            s.atmos_profile = AtmosProfile.UserWaterAndOzone(tcwv, ozone)
                            s.aot550 = aot
                            s.wavelength = Wavelength(wl / 1000.0)  # nm to um

                            # Run simulation
                            s.run()

                            # Extract outputs
                            path_ref[i_sza, i_vza, i_raa, i_aot, i_tcwv, i_wl] = (
                                s.outputs.atmospheric_intrinsic_reflectance
                            )
                            trans_down[i_sza, i_vza, i_raa, i_aot, i_tcwv, i_wl] = (
                                s.outputs.transmittance_total_scattering.downward
                            )
                            trans_up[i_sza, i_vza, i_raa, i_aot, i_tcwv, i_wl] = (
                                s.outputs.transmittance_total_scattering.upward
                            )
                            sph_alb[i_sza, i_vza, i_raa, i_aot, i_tcwv, i_wl] = (
                                s.outputs.spherical_albedo
                            )

                            count += 1
                            if count % 1000 == 0:
                                logger.info(f"Progress: {count}/{total_sims} ({100*count/total_sims:.1f}%)")

    # Create xarray dataset
    ds = xr.Dataset(
        {
            "path_reflectance": (["sza", "vza", "raa", "aot", "tcwv", "wavelength"], path_ref),
            "transmittance_down": (["sza", "vza", "raa", "aot", "tcwv", "wavelength"], trans_down),
            "transmittance_up": (["sza", "vza", "raa", "aot", "tcwv", "wavelength"], trans_up),
            "spherical_albedo": (["sza", "vza", "raa", "aot", "tcwv", "wavelength"], sph_alb),
        },
        coords={
            "sza": sza_arr,
            "vza": vza_arr,
            "raa": raa_arr,
            "aot": aot_arr,
            "tcwv": tcwv_arr,
            "wavelength": wl_arr,
        },
    )

    ds.attrs["aerosol_type"] = aerosol_type
    ds.attrs["ozone"] = ozone
    ds.attrs["creation_date"] = str(np.datetime64("now"))

    # Save to Zarr
    logger.info(f"Saving LUT to {output_path}")
    ds.to_zarr(output_path, mode="w", consolidated=True)
    logger.info("LUT creation complete")
