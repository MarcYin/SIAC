"""Xarray-backed runtime payload models used during SIAC execution."""

from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING, Any

import numpy as np
import xarray as xr

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping, Sequence

    from siac.domain.sensors import SensorBand, SensorConfig


@dataclass
class BandLoaderContext:
    """Wraps a TOA dataset with an optional lazy band loader and cache.

    Instead of storing callables in xarray Dataset attrs (which is fragile
    and non-discoverable), pass this context object alongside the dataset.
    """

    toa: xr.Dataset
    band_loader: Callable[[str], xr.DataArray] | None = None
    _cache: dict[str, xr.DataArray] = field(default_factory=dict)

    def get_band(self, band_name: str) -> xr.DataArray:
        """Get a band from the TOA dataset, loading lazily if needed."""
        band = self.toa.data_vars.get(band_name)
        if band is not None:
            return band

        cached = self._cache.get(band_name)
        if cached is not None:
            return cached

        if self.band_loader is None:
            raise KeyError(band_name)

        loaded = self.band_loader(band_name)
        if not isinstance(loaded, xr.DataArray):
            raise TypeError(
                f"Band loader must return xr.DataArray for {band_name!r}, "
                f"got {type(loaded).__name__}"
            )
        self._cache[band_name] = loaded
        return loaded


def _as_data_array(value: object) -> xr.DataArray:
    if not isinstance(value, xr.DataArray):
        if isinstance(value, (np.ndarray, np.generic, int, float)):
            return xr.DataArray(np.asarray(value))
        raise TypeError(f"Expected xr.DataArray, got {type(value).__name__}")
    return value


# ``copy_spatial_metadata_like`` lives in ``siac.geo._spatial`` so that
# ``siac.geo.resample`` and the surface-prior pipelines can use it without
# the previous ``geo -> runtime`` layering inversion (REVIEW.md §1.4). The
# re-export here keeps backward compatibility for any caller still doing
# ``from siac.runtime.models import copy_spatial_metadata_like``. ``E402``
# is suppressed because the import is intentionally below the docstring
# preamble — the comment block explains why.
from siac.geo._spatial import copy_spatial_metadata_like as copy_spatial_metadata_like  # noqa: E402


@dataclass(frozen=True)
class GeometryAngles:
    """View and sun geometry for atmospheric correction.

    Internal SIAC convention:

    - ``sza``, ``saa``, ``vza``, ``vaa`` are stored in radians.
    - ``raa`` is derived as ``vaa - saa`` in radians.
    - neural-network emulators consume ``cos(sza)``, ``cos(vza)``, and
      ``cos(raa)``.
    - the LUT backend converts the radian fields back to degrees before LUT
      interpolation.
    """

    sza: xr.DataArray
    saa: xr.DataArray
    vza: xr.DataArray
    vaa: xr.DataArray

    def __post_init__(self) -> None:
        for name in ("sza", "saa", "vza", "vaa"):
            val = getattr(self, name)
            if not isinstance(val, xr.DataArray):
                raise TypeError(
                    f"GeometryAngles.{name} must be xr.DataArray, got {type(val).__name__}"
                )

    @property
    def raa(self) -> xr.DataArray:
        """Relative azimuth angle in radians, wrapped to [0, 2*pi)."""
        return _as_data_array(np.abs(self.vaa - self.saa) % (2.0 * np.pi))

    @property
    def cos_sza(self) -> xr.DataArray:
        return _as_data_array(np.cos(self.sza))

    @property
    def cos_vza(self) -> xr.DataArray:
        return _as_data_array(np.cos(self.vza))

    @property
    def cos_raa(self) -> xr.DataArray:
        return _as_data_array(np.cos(self.raa))

    def to_emulator_input(self) -> xr.Dataset:
        """Return emulator-ready geometry features from the radian angle fields."""
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
            sza=_as_data_array(sza_deg * deg_to_rad),
            saa=_as_data_array(saa_deg * deg_to_rad),
            vza=_as_data_array(vza_deg * deg_to_rad),
            vaa=_as_data_array(vaa_deg * deg_to_rad),
        )


@dataclass(frozen=True)
class AtmosphericState:
    """Atmospheric parameters for radiative transfer calculations.

    Shared SIAC RT units used by the solver, LUT backend, and emulator backend:

    - ``aot``: aerosol optical depth at 550 nm, unitless
    - ``tcwv``: total column water vapour in cm precipitable water
      (numerically equal to g cm^-2)
    - ``tco3``: total column ozone in atm-cm (DU / 1000)
    - ``elevation``: terrain altitude in km above mean sea level
    """

    aot: xr.DataArray
    tcwv: xr.DataArray
    tco3: xr.DataArray
    aot_unc: xr.DataArray
    tcwv_unc: xr.DataArray
    tco3_unc: xr.DataArray
    elevation: xr.DataArray

    def __post_init__(self) -> None:
        for name in ("aot", "tcwv", "tco3", "aot_unc", "tcwv_unc", "tco3_unc", "elevation"):
            val = getattr(self, name)
            if not isinstance(val, xr.DataArray):
                raise TypeError(
                    f"AtmosphericState.{name} must be xr.DataArray, got {type(val).__name__}"
                )

    def to_emulator_input(self, geometry: GeometryAngles) -> xr.Dataset:
        """Return the shared RT state in the emulator feature convention.

        Geometry is emitted as cosines of the radian angle fields. Atmospheric
        variables are passed through unchanged in the shared SIAC RT units:
        ``aot`` unitless, ``tcwv`` in cm, ``tco3`` in atm-cm, and elevation in
        km.
        """
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

    def get_parameter(self, name: str) -> xr.DataArray:
        """Return one atmospheric parameter field by canonical parameter name."""
        if name not in {"aot", "tcwv", "tco3"}:
            raise KeyError(f"Unknown atmospheric parameter {name!r}")
        return getattr(self, name)

    def get_uncertainty(self, name: str) -> xr.DataArray:
        """Return one atmospheric parameter uncertainty field by parameter name."""
        if name not in {"aot", "tcwv", "tco3"}:
            raise KeyError(f"Unknown atmospheric parameter {name!r}")
        return getattr(self, f"{name}_unc")

    def with_updated_parameters(
        self,
        values: Mapping[str, xr.DataArray],
        uncertainties: Mapping[str, xr.DataArray] | None = None,
    ) -> AtmosphericState:
        """Return a state with selected atmospheric parameters replaced.

        ``values`` and ``uncertainties`` keys are canonical parameter names:
        ``aot``, ``tcwv``, and ``tco3``.
        """
        allowed = {"aot", "tcwv", "tco3"}
        fields = {name: getattr(self, name) for name in allowed}
        uncs = {name: getattr(self, f"{name}_unc") for name in allowed}

        for name, value in values.items():
            if name not in allowed:
                raise KeyError(f"Unknown atmospheric parameter {name!r}")
            fields[name] = copy_spatial_metadata_like(value, getattr(self, name))

        if uncertainties is not None:
            for name, value in uncertainties.items():
                if name not in allowed:
                    raise KeyError(f"Unknown atmospheric uncertainty {name!r}")
                uncs[name] = copy_spatial_metadata_like(
                    value,
                    getattr(self, f"{name}_unc"),
                )

        return AtmosphericState(
            aot=fields["aot"],
            tcwv=fields["tcwv"],
            tco3=fields["tco3"],
            aot_unc=uncs["aot"],
            tcwv_unc=uncs["tcwv"],
            tco3_unc=uncs["tco3"],
            elevation=self.elevation,
        )

    def with_updated_aot_tcwv(
        self,
        aot: xr.DataArray,
        tcwv: xr.DataArray,
        aot_unc: xr.DataArray | None = None,
        tcwv_unc: xr.DataArray | None = None,
    ) -> AtmosphericState:
        uncertainties: dict[str, xr.DataArray] = {}
        if aot_unc is not None:
            uncertainties["aot"] = aot_unc
        if tcwv_unc is not None:
            uncertainties["tcwv"] = tcwv_unc
        return self.with_updated_parameters(
            {"aot": aot, "tcwv": tcwv},
            uncertainties or None,
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
    extras: Mapping[str, xr.DataArray] = field(default_factory=dict)

    def __post_init__(self) -> None:
        reserved = {"xap", "xbp", "xcp", "d_xap", "d_xbp", "d_xcp", "extras"}
        normalized_extras: dict[str, xr.DataArray] = {}
        for name, value in self.extras.items():
            key = str(name).strip()
            if not key:
                raise ValueError("RTCoefficients extras keys must be non-empty strings")
            if key in reserved:
                raise ValueError(
                    f"RTCoefficients extras key {key!r} conflicts with a reserved field"
                )
            normalized_extras[key] = _as_data_array(value)
        object.__setattr__(self, "extras", MappingProxyType(normalized_extras))

    @property
    def has_jacobian(self) -> bool:
        return self.d_xap is not None

    @property
    def output_names(self) -> tuple[str, ...]:
        return ("xap", "xbp", "xcp", *self.extras.keys())

    def get_output(self, name: str) -> xr.DataArray:
        if name == "xap":
            return self.xap
        if name == "xbp":
            return self.xbp
        if name == "xcp":
            return self.xcp
        try:
            return self.extras[name]
        except KeyError as exc:
            raise KeyError(
                f"Unknown RT output {name!r}. Available outputs: {', '.join(self.output_names)}"
            ) from exc

    def to_output_arrays(
        self,
        names: Sequence[str] | None = None,
    ) -> dict[str, xr.DataArray]:
        selected = self.output_names if names is None else tuple(names)
        return {name: self.get_output(name) for name in selected}

    def resample_to_template(self, template: xr.DataArray) -> RTCoefficients:
        """Resample every field of this coefficient bundle to *template* grid.

        The function used to live in ``siac.geo.resample`` and reach back
        into ``siac.runtime`` to reconstruct ``RTCoefficients``, creating
        the runtime↔geo cycle that REVIEW.md §1.4 flagged. Moving the
        construction onto the class itself lets ``siac.geo.resample``
        keep its public wrapper as a one-liner that just delegates here.
        """
        from siac.geo.resample import (
            resample_field_or_param_stack_to_template,
            resample_field_to_template,
        )

        resampled_extras: dict[str, xr.DataArray] = {}
        # ``extra_field`` (not ``field``) — the latter is the
        # ``dataclasses.field`` imported above; reusing the name as a loop
        # variable triggers F402.
        for name, extra_field in self.extras.items():
            resampled = resample_field_or_param_stack_to_template(extra_field, template)
            resampled_extras[name] = extra_field if resampled is None else resampled

        return RTCoefficients(
            xap=resample_field_to_template(self.xap, template),
            xbp=resample_field_to_template(self.xbp, template),
            xcp=resample_field_to_template(self.xcp, template),
            d_xap=resample_field_or_param_stack_to_template(self.d_xap, template),
            d_xbp=resample_field_or_param_stack_to_template(self.d_xbp, template),
            d_xcp=resample_field_or_param_stack_to_template(self.d_xcp, template),
            extras=resampled_extras,
        )

    def apply_correction(self, toa: xr.DataArray) -> xr.DataArray:
        y = _as_data_array(self.xap * toa - self.xbp)
        denom = _as_data_array(1.0 + self.xcp * y)
        stable = _as_data_array(np.isfinite(denom) & (np.abs(denom) > 1.0e-10))
        return _as_data_array((y / denom).where(stable))

    def simulate_toa(self, boa: xr.DataArray) -> xr.DataArray:
        """Forward-simulate dimensionless TOA reflectance from BOA reflectance."""
        denom = _as_data_array(1.0 - self.xcp * boa)
        stable = _as_data_array(
            np.isfinite(denom) & (np.abs(denom) > 1.0e-6) & (np.abs(self.xap) > 1.0e-12)
        )
        y = _as_data_array(boa / denom)
        toa = _as_data_array((y + self.xbp) / self.xap)
        return _as_data_array(toa.where(stable))

    def compute_boa_jacobian(
        self,
        toa: xr.DataArray,
    ) -> tuple[xr.DataArray, xr.DataArray]:
        if not self.has_jacobian:
            raise ValueError("RTCoefficients does not have Jacobian information")
        if self.d_xap is None or self.d_xbp is None or self.d_xcp is None:
            raise ValueError("Jacobian arrays (d_xap, d_xbp, d_xcp) must all be non-None")

        y = _as_data_array(self.xap * toa - self.xbp)
        denom = _as_data_array(1.0 + self.xcp * y)

        d_y_aot = _as_data_array(self.d_xap.sel(param="aot") * toa - self.d_xbp.sel(param="aot"))
        d_y_tcwv = _as_data_array(self.d_xap.sel(param="tcwv") * toa - self.d_xbp.sel(param="tcwv"))

        d_denom_aot = _as_data_array(self.d_xcp.sel(param="aot") * y + self.xcp * d_y_aot)
        d_denom_tcwv = _as_data_array(self.d_xcp.sel(param="tcwv") * y + self.xcp * d_y_tcwv)

        d_boa_aot = _as_data_array((d_y_aot * denom - y * d_denom_aot) / (denom**2))
        d_boa_tcwv = _as_data_array((d_y_tcwv * denom - y * d_denom_tcwv) / (denom**2))

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
        self,
        k_vol: xr.DataArray,
        k_geo: xr.DataArray,
    ) -> xr.DataArray:
        return _as_data_array(self.f0 + self.f1 * k_vol + self.f2 * k_geo)

    def compute_reflectance_uncertainty(
        self,
        k_vol: xr.DataArray,
        k_geo: xr.DataArray,
    ) -> xr.DataArray:
        if self.reflectance_unc is not None:
            return self.reflectance_unc
        var = self.f0_unc**2 + (k_vol * self.f1_unc) ** 2 + (k_geo * self.f2_unc) ** 2
        return _as_data_array(np.sqrt(var))


@dataclass(frozen=True)
class SurfacePrior:
    """Surface reflectance prior derived from BRDF model."""

    boa: xr.DataArray
    boa_unc: xr.DataArray
    kernels: BRDFKernelWeights | None
    mask: xr.DataArray
    monthly_composites: tuple[Any, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class MonthlyCompositeOutput:
    """Persistable monthly-composite product bundle."""

    reflectance: xr.Dataset
    quality: xr.DataArray
    sample_index: xr.DataArray


@dataclass(frozen=True)
class ObservationBundle:
    """Complete observation data from satellite preprocessing."""

    toa: xr.Dataset
    geometry: GeometryAngles
    cloud_mask: xr.DataArray
    sensor_config: SensorConfig
    metadata: dict[str, Any]
    crs: str
    bounds: tuple[float, float, float, float]


@dataclass(frozen=True)
class SolverInputBundle:
    """All inputs to the aerosol solver, resampled to solver grids."""

    toa: xr.DataArray
    geometry: GeometryAngles
    cloud_mask: xr.DataArray
    sensor_config: SensorConfig
    bands: list[SensorBand]
    atmo_prior: AtmosphericState
    surface_prior: SurfacePrior
    rt_model: Any
    aux_resolution_m: float
    aerosol_resolution_m: float
    sharp_transition_mask: xr.DataArray | None = None
    water_mask: xr.DataArray | None = None


@dataclass(frozen=True)
class AOTScatterBandDiagnostics:
    """Sampled scatter diagnostics for one aerosol-solver band."""

    band_name: str
    surface_reflectance: np.ndarray
    observed_toa: np.ndarray
    simulated_toa: np.ndarray
    total_valid_count: int


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
    qa: xr.Dataset | None = None
    level_history: tuple[dict[str, Any], ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class CorrectionDiagnostics:
    """Structured diagnostics produced by the correction stage."""

    processing_time_s: float | None = None
    aot_scatter_plots: tuple[AOTScatterBandDiagnostics, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class CorrectionResult:
    """Final output of atmospheric correction."""

    boa: xr.Dataset
    boa_unc: xr.Dataset | None
    aot: xr.DataArray
    tcwv: xr.DataArray
    cloud_mask: xr.DataArray
    #: Raw input cloud mask from M1 (OmniCloudMask), before the correction
    #: stage ORs in per-band invalid-BOA pixels. ``cloud_mask`` is the union
    #: (True = cloudy OR invalid-BOA) for output masking; this field isolates
    #: actual cloud so downstream consumers can distinguish "cloudy" from
    #: "correction flagged the pixel unreliable".
    cloud_mask_m1: xr.DataArray | None = None
    surface_prior: xr.Dataset | None = None
    surface_prior_unc: xr.Dataset | None = None
    solver_qa: xr.Dataset | None = None
    monthly_composites: dict[str, MonthlyCompositeOutput] | None = None
    metadata: dict[str, Any] = field(default_factory=dict)
    diagnostics: CorrectionDiagnostics = field(default_factory=CorrectionDiagnostics)
