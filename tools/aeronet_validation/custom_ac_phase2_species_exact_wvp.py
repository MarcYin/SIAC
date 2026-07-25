"""Build species-resolved 6S coefficients at each historical scene's WVP.

Unlike the legacy five-node 0.5-4.5 cm LUT, this computes one coefficient set
at the scene's recorded Sen2Cor L2A WVP. The output is merged with
``atmo_sidecar_merge_exact_wvp.py`` so every scene carries its own singleton
TCWV node and no interpolation or clamping is possible.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
sys.path.insert(0, str(ROOT))

import aerosol_species as asp  # noqa: E402
import custom_ac_phase2_species as legacy  # noqa: E402

from siac.algorithms.rt.direct.sixs import SixSBackend  # noqa: E402
from siac.catalog.sensors.sentinel2 import SENTINEL2A_CONFIG  # noqa: E402
from siac.config import SixSAlgorithmConfig  # noqa: E402
from siac.config.algorithms import RTAerosolSetupConfig, RTSetupConfig  # noqa: E402
from siac.runtime import AtmosphericState, GeometryAngles  # noqa: E402

BASE_SIDECAR_BANDS = ["coastal", "blue", "green", "red", "nir08", "swir16", "swir22"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--meta", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--lon", type=float, required=True)
    parser.add_argument("--lat", type=float, required=True)
    parser.add_argument("--month", type=int, required=True)
    args = parser.parse_args()

    bands = legacy.build_bands()
    sidecar_bands = list(BASE_SIDECAR_BANDS)
    if os.environ.get("EXTENDED_ACIX_ANCHORS", "0") == "1":
        bands.extend(
            [
                SENTINEL2A_CONFIG.get_band("B06"),
                SENTINEL2A_CONFIG.get_band("B08"),
            ]
        )
        sidecar_bands.extend(["rededge2", "nirwide"])
    fraction = asp.candidate_fraction_sets(args.lon, args.lat, args.month, n=1)[0]
    distribution = asp.cci_distribution(fraction)
    elevation = legacy.site_elevation_km(args.lon, args.lat)
    setup = RTSetupConfig(
        aerosol=RTAerosolSetupConfig(profile="multimodal_log_normal", distribution=distribution)
    )
    config = SixSAlgorithmConfig(
        build_profile="release",
        output_variables=("xap", "xbp", "xcp"),
        native_threads=2,
    ).model_copy(update={"module_path": legacy.MOD, "auto_build": False})
    backend = SixSBackend(sixs_config=config, rt_setup=setup)
    backend.set_observation_time(datetime(2022, 7, 15, 8, 45))

    metadata = json.loads(Path(args.meta).read_text(encoding="utf-8"))["selected"]
    coefficients = np.full((len(metadata), 2, 1, 3, len(bands)), np.nan, dtype=np.float64)
    tcwv_by_day = np.full(len(metadata), np.nan, dtype=np.float64)
    for day_index, day in enumerate(metadata):
        tcwv = float(day["wvp"])
        if not math.isfinite(tcwv) or not 0.05 <= tcwv <= 15.0:
            raise ValueError(f"invalid historical WVP {tcwv!r} for {day.get('day')}")
        tcwv_by_day[day_index] = tcwv
        geometry = GeometryAngles.from_degrees(
            legacy.field(day["sza"]),
            legacy.field(day["saa"]),
            legacy.field(day["vza"]),
            legacy.field(day["vaa"]),
        )
        aods = [
            float(day["maiac"]),
            float(day["aeronet_op"] if day["aeronet_op"] is not None else day["maiac"]),
        ]
        computed: dict[float, int] = {}
        for aod_index, aod in enumerate(aods):
            if aod in computed:
                coefficients[day_index, aod_index] = coefficients[day_index, computed[aod]]
                continue
            computed[aod] = aod_index
            atmosphere = AtmosphericState(
                aot=legacy.field(aod),
                tcwv=legacy.field(tcwv),
                tco3=legacy.field(0.30),
                aot_unc=legacy.field(0.02),
                tcwv_unc=legacy.field(0.05),
                tco3_unc=legacy.field(0.01),
                elevation=legacy.field(elevation),
            )
            band_coefficients = backend.compute_coefficients_multi(geometry, atmosphere, bands)
            for variable_index, variable in enumerate(("xap", "xbp", "xcp")):
                coefficients[day_index, aod_index, 0, variable_index] = [
                    float(result.get_output(variable).values.ravel()[0])
                    for result in band_coefficients
                ]
        print(
            f"{day['day']}: exact-WVP species 6S done "
            f"(aot {aods[0]:.3f}/{aods[1]:.3f}, tcwv {tcwv:.3f} cm)",
            flush=True,
        )

    np.savez(
        args.out,
        coeffs=coefficients,
        bands=np.asarray(sidecar_bands),
        tcwv_by_day=tcwv_by_day,
    )
    print(
        f"saved {args.out} {coefficients.shape} "
        f"(exact WVP {tcwv_by_day.min():.3f}-{tcwv_by_day.max():.3f} cm)",
        flush=True,
    )


if __name__ == "__main__":
    main()
