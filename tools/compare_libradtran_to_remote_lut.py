"""Compare the direct libRadtran (uvspec) backend to the remote libRadtran LUT.

Both paths ultimately run the same SIAC spectral math (RSRF convolution +
``derive_standard_rt_coefficients``); the only difference is that this backend
computes the ``TOA_rho``/``Eg_rho`` terms on the fly with ``uvspec`` rather than
reading them from the precomputed remote LUT. Tight agreement is therefore the
real acceptance test of the reconstructed uvspec deck (the preset, the band
model, and the two-albedo normalization).

Usage (requires a built uvspec — ``pixi run -e libradtran build-libradtran`` —
and network access to the remote LUT)::

    pixi run -e libradtran python tools/compare_libradtran_to_remote_lut.py
    pixi run -e libradtran python tools/compare_libradtran_to_remote_lut.py --bands B02,B04,B8A --json report.json
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass

import numpy as np
import xarray as xr

from siac.algorithms.rt.direct.libradtran import LibRadtranBackend
from siac.algorithms.rt.lut import DEFAULT_LUT_URL, ZarrLUTBackend
from siac.config.algorithms import LibRadtranAlgorithmConfig
from siac.domain.sensors import SensorBand
from siac.rt_setup import DEFAULT_LIBRADTRAN_RT_SETUP, DEFAULT_LUT_RT_SETUP
from siac.runtime import AtmosphericState, GeometryAngles, RTCoefficients

# Sentinel-2 visible/NIR/SWIR bands (nominal centre/bandwidth in nm).
_BAND_CATALOG: dict[str, SensorBand] = {
    "B01": SensorBand(name="B01", center_wavelength=443.0, bandwidth=20.0, resolution=60.0, band_index=0),
    "B02": SensorBand(name="B02", center_wavelength=492.0, bandwidth=66.0, resolution=10.0, band_index=1),
    "B03": SensorBand(name="B03", center_wavelength=560.0, bandwidth=36.0, resolution=10.0, band_index=2),
    "B04": SensorBand(name="B04", center_wavelength=665.0, bandwidth=31.0, resolution=10.0, band_index=3),
    "B8A": SensorBand(name="B8A", center_wavelength=865.0, bandwidth=21.0, resolution=20.0, band_index=8),
}


@dataclass(frozen=True)
class PointCase:
    name: str
    sza: float
    vza: float
    raa: float
    aot: float
    tcwv: float  # cm
    tco3: float  # atm-cm
    elevation: float  # km


def _cases() -> list[PointCase]:
    return [
        PointCase("clear_nadir", 30.0, 5.0, 60.0, 0.10, 1.5, 0.30, 0.05),
        PointCase("hazy_oblique", 45.0, 20.0, 120.0, 0.40, 2.5, 0.32, 0.20),
        PointCase("dry_lowsun", 60.0, 10.0, 90.0, 0.20, 0.8, 0.28, 0.50),
    ]


def _scene(case: PointCase) -> tuple[GeometryAngles, AtmosphericState]:
    def _u(value: float) -> xr.DataArray:
        return xr.DataArray(np.full((2, 2), value, dtype=np.float64), dims=("y", "x"))

    geom = GeometryAngles.from_degrees(_u(case.sza), _u(150.0), _u(case.vza), _u(150.0 - case.raa))
    zeros = _u(0.0)
    atmo = AtmosphericState(
        aot=_u(case.aot), tcwv=_u(case.tcwv), tco3=_u(case.tco3),
        aot_unc=zeros, tcwv_unc=zeros, tco3_unc=zeros, elevation=_u(case.elevation),
    )
    return geom, atmo


def _scalars(coeffs: RTCoefficients) -> dict[str, float]:
    return {
        "xap": float(np.nanmean(coeffs.xap.values)),
        "xbp": float(np.nanmean(coeffs.xbp.values)),
        "xcp": float(np.nanmean(coeffs.xcp.values)),
    }


def _rel(ref: float, cand: float) -> float:
    denom = max(abs(ref), 1e-6)
    return abs(cand - ref) / denom


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bands", default="B02,B04,B8A", help="Comma-separated band names.")
    parser.add_argument("--lut-path", default=DEFAULT_LUT_URL, help="Remote/local LUT path.")
    parser.add_argument("--wavelength-max", type=float, default=2500.0)
    parser.add_argument("--nodes", type=int, default=3, help="Grid nodes per aot/tcwv axis.")
    parser.add_argument("--json", type=str, default=None, help="Optional path for a JSON report.")
    args = parser.parse_args()

    bands = [_BAND_CATALOG[name] for name in args.bands.split(",")]
    lr_cfg = LibRadtranAlgorithmConfig(
        wavelength_min_nm=400.0,
        wavelength_max_nm=float(args.wavelength_max),
        scene_lut_max_nodes_per_axis=int(args.nodes),
        scene_lut_max_cases=2 * int(args.nodes) ** 2,
    )
    libradtran = LibRadtranBackend(libradtran_config=lr_cfg, rt_setup=DEFAULT_LIBRADTRAN_RT_SETUP)
    remote = ZarrLUTBackend(args.lut_path, rt_setup=DEFAULT_LUT_RT_SETUP)

    report: list[dict] = []
    print(f"{'case':16s} {'band':5s} {'coeff':4s} {'libradtran':>12s} {'remoteLUT':>12s} {'rel_diff':>10s}")
    for case in _cases():
        geom, atmo = _scene(case)
        for band in bands:
            lr = _scalars(libradtran.compute_coefficients(geom, atmo, band))
            rm = _scalars(remote.compute_coefficients(geom, atmo, band))
            for coeff in ("xap", "xbp", "xcp"):
                rel = _rel(rm[coeff], lr[coeff])
                report.append(
                    {"case": case.name, "band": band.name, "coeff": coeff,
                     "libradtran": lr[coeff], "remote_lut": rm[coeff], "rel_diff": rel}
                )
                print(f"{case.name:16s} {band.name:5s} {coeff:4s} {lr[coeff]:12.5f} {rm[coeff]:12.5f} {rel:10.4f}")

    rels = [r["rel_diff"] for r in report]
    print(f"\nrel_diff: mean={np.mean(rels):.4f} median={np.median(rels):.4f} max={np.max(rels):.4f}")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(
                {"identity": {"generator": "libRadtran", "aerosol_profile": "continental_average",
                              "atmospheric_profile_shape": "us_standard"},
                 "entries": report,
                 "summary": {"mean_rel": float(np.mean(rels)), "median_rel": float(np.median(rels)),
                             "max_rel": float(np.max(rels))}},
                fh, indent=2,
            )
        print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
