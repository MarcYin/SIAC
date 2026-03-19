"""
Core data contracts for SIAC atmospheric correction.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import xarray as xr


@dataclass(frozen=True)
class GeometryAngles:
    """View and sun geometry for atmospheric correction."""

    sza: xr.DataArray
    saa: xr.DataArray
    vza: xr.DataArray
    vaa: xr.DataArray

    @property
    def raa(self) -> xr.DataArray:
        return self.vaa - self.saa

    @property
    def cos_sza(self) -> xr.DataArray:
        return np.cos(self.sza)

    @property
    def cos_vza(self) -> xr.DataArray:
        return np.cos(self.vza)

    @property
    def cos_raa(self) -> xr.DataArray:
        return np.cos(self.raa)

    def to_emulator_input(self) -> xr.Dataset:
        return xr.Dataset(
            {
                "cos_sza": self.cos_sza,
                "cos_vza": self.cos_vza,
                "cos_raa": self.cos_raa,
            }
        )

    @classmethod
    def from_degrees(
        cls,
        sza_deg: xr.DataArray,
        saa_deg: xr.DataArray,
        vza_deg: xr.DataArray,
        vaa_deg: xr.DataArray,
    ) -> GeometryAngles:
        deg_to_rad = np.pi / 180.0
        return cls(
            sza=sza_deg * deg_to_rad,
            saa=saa_deg * deg_to_rad,
            vza=vza_deg * deg_to_rad,
            vaa=vaa_deg * deg_to_rad,
        )


@dataclass(frozen=True)
class AtmosphericState:
    """Atmospheric parameters for radiative transfer calculations."""

    aot: xr.DataArray
    tcwv: xr.DataArray
    tco3: xr.DataArray
    aot_unc: xr.DataArray
    tcwv_unc: xr.DataArray
    tco3_unc: xr.DataArray
    elevation: xr.DataArray

    def to_emulator_input(self, geometry: GeometryAngles) -> xr.Dataset:
        return xr.Dataset(
            {
                "cos_sza": geometry.cos_sza,
                "cos_vza": geometry.cos_vza,
                "cos_raa": geometry.cos_raa,
                "aot": self.aot,
                "tcwv": self.tcwv,
                "tco3": self.tco3,
                "elevation": self.elevation,
            }
        )

    def with_updated_aot_tcwv(
        self,
        aot: xr.DataArray,
        tcwv: xr.DataArray,
        aot_unc: xr.DataArray | None = None,
        tcwv_unc: xr.DataArray | None = None,
    ) -> AtmosphericState:
        return AtmosphericState(
            aot=aot,
            tcwv=tcwv,
            tco3=self.tco3,
            aot_unc=aot_unc if aot_unc is not None else self.aot_unc,
            tcwv_unc=tcwv_unc if tcwv_unc is not None else self.tcwv_unc,
            tco3_unc=self.tco3_unc,
            elevation=self.elevation,
        )


@dataclass(frozen=True)
class RTCoefficients:
    """Radiative transfer coefficients from emulator or LUT."""

    xap: xr.DataArray
    xbp: xr.DataArray
    xcp: xr.DataArray
    d_xap: xr.DataArray | None = None
    d_xbp: xr.DataArray | None = None
    d_xcp: xr.DataArray | None = None

    @property
    def has_jacobian(self) -> bool:
        return self.d_xap is not None

    def apply_correction(self, toa: xr.DataArray) -> xr.DataArray:
        y = self.xap * toa - self.xbp
        return y / (1.0 + self.xcp * y)

    def compute_boa_jacobian(
        self, toa: xr.DataArray
    ) -> tuple[xr.DataArray, xr.DataArray]:
        if not self.has_jacobian:
            raise ValueError("RTCoefficients does not have Jacobian information")

        y = self.xap * toa - self.xbp
        denom = 1.0 + self.xcp * y

        d_y_aot = self.d_xap.sel(param="aot") * toa - self.d_xbp.sel(param="aot")
        d_y_tcwv = self.d_xap.sel(param="tcwv") * toa - self.d_xbp.sel(param="tcwv")

        d_denom_aot = self.d_xcp.sel(param="aot") * y
        d_denom_tcwv = self.d_xcp.sel(param="tcwv") * y

        d_boa_aot = (d_y_aot * denom - y * d_denom_aot) / (denom**2)
        d_boa_tcwv = (d_y_tcwv * denom - y * d_denom_tcwv) / (denom**2)

        return d_boa_aot, d_boa_tcwv


@dataclass(frozen=True)
class BRDFKernelWeights:
    """Ross-Thick Li-Sparse BRDF kernel weights."""

    f0: xr.DataArray
    f1: xr.DataArray
    f2: xr.DataArray
    f0_unc: xr.DataArray
    f1_unc: xr.DataArray
    f2_unc: xr.DataArray
    reflectance_unc: xr.DataArray | None = None

    def compute_reflectance(
        self, k_vol: xr.DataArray, k_geo: xr.DataArray
    ) -> xr.DataArray:
        return self.f0 + self.f1 * k_vol + self.f2 * k_geo

    def compute_reflectance_uncertainty(
        self, k_vol: xr.DataArray, k_geo: xr.DataArray
    ) -> xr.DataArray:
        if self.reflectance_unc is not None:
            return self.reflectance_unc
        var = self.f0_unc**2 + (k_vol * self.f1_unc) ** 2 + (k_geo * self.f2_unc) ** 2
        return np.sqrt(var)


@dataclass(frozen=True)
class SurfacePrior:
    """Surface reflectance prior derived from BRDF model."""

    boa: xr.DataArray
    boa_unc: xr.DataArray
    kernels: BRDFKernelWeights | None
    mask: xr.DataArray


@dataclass(frozen=True)
class ObservationBundle:
    """Complete observation data from satellite preprocessing."""

    toa: xr.Dataset
    geometry: GeometryAngles
    cloud_mask: xr.DataArray
    sensor_config: Any
    metadata: dict[str, Any]
    crs: str
    bounds: tuple[float, float, float, float]


@dataclass(frozen=True)
class SolverInputBundle:
    """All inputs to the aerosol solver, resampled to solver grids."""

    toa: xr.DataArray
    geometry: GeometryAngles
    cloud_mask: xr.DataArray
    sensor_config: Any
    bands: list[Any]
    atmo_prior: AtmosphericState
    surface_prior: SurfacePrior
    rt_model: Any
    aux_resolution_m: float
    aerosol_resolution_m: float


@dataclass(frozen=True)
class SolvedAtmosphere:
    """Solver output: retrieved atmospheric parameters + diagnostics."""

    atmo_state: AtmosphericState
    aot: xr.DataArray
    tcwv: xr.DataArray
    aot_unc: xr.DataArray
    tcwv_unc: xr.DataArray
    cost_final: float
    n_iterations: int
    converged: bool


@dataclass(frozen=True)
class CorrectionResult:
    """Final output of atmospheric correction."""

    boa: xr.Dataset
    boa_unc: xr.Dataset | None
    aot: xr.DataArray
    tcwv: xr.DataArray
    cloud_mask: xr.DataArray
    metadata: dict[str, Any]
