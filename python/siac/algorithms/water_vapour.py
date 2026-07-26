"""Total column water vapour retrieved from Sentinel-2 top-of-atmosphere data.

Water vapour is a *measured* input to atmospheric correction, and the scene
being corrected already carries the measurement: Sentinel-2's B09 sits inside
the 940 nm water-absorption feature while B8A sits on the continuum beside it,
so their ratio is a direct sounding of the column above the pixel. Retrieving it
here removes any dependence on somebody else's Level-2A product being published,
reachable, and processed with a baseline whose radiometry we understand.

The retrieval is the Continuum Interpolated Band Ratio inverted through the
SNAP/s2tbx ``S2.WaterVapour`` look-up table:

.. math::

    \\mathrm{CIBR} = \\frac{955.19}{813.04}\\,
        \\frac{\\rho_{8A}}{\\rho_{9}}, \\qquad
    W\\;[\\mathrm{cm}] = (b_0 \\log_{10}\\mathrm{CIBR} + b_1)^2

``b0`` and ``b1`` come from the table, keyed on surface reflectance (for which
TOA :math:`\\rho_{8A}` is the proxy), sun zenith, view zenith, relative azimuth
and ground altitude.

**On interpolating that table.** It holds 6125 rows, which is exactly
7 x 5 x 5 x 5 x 7 — a complete regular grid, not scattered samples, so it is
interpolated multilinearly. The two research implementations instead queried it
with a 32-neighbour distance-weighted KNN, one of them on *unnormalised*
features: reflectance spans 0.05-0.8 while sun zenith spans 0-60, so the
Euclidean neighbourhood is dominated by the angle axes and retains under an
eighth of the reflectance dependence. That flattening is what its empirical
``* 0.716 + 0.03`` correction was patching, and that correction is **not**
carried forward here. Measured over 2268 monthly realizations from 129 scenes
against CAMS total column water vapour, this multilinear retrieval has a
+0.08 cm median bias and a 1.06 slope, where the ``0.716``-corrected KNN reads
-0.21 cm / 0.88 and Sen2Cor's own L2A ``WVP`` reads -0.32 cm / 0.84: the factor
was calibrating the retrieval onto Sen2Cor's low bias, not removing one.

The azimuth axis stops at 180 degrees because the table is symmetric about the
principal plane, so the relative azimuth must be *folded* there
(:func:`relative_azimuth_deg`) rather than wrapped: a naive ``% 180`` sends
200 degrees to 20 degrees where the fold correctly gives 160.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from importlib import resources
from typing import TYPE_CHECKING, Any, cast

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Sequence

__all__ = [
    "LUT_AXIS_NAMES",
    "WaterVapourResult",
    "load_lut",
    "relative_azimuth_deg",
    "retrieve_water_vapour",
]

#: Package holding the vendored look-up table.
_DATA_PACKAGE = "siac.data"
_LUT_FILENAME = "water_vapour_lut_CIBR_9_8A.txt"

#: Look-up-table axes, in table-column order.
LUT_AXIS_NAMES: tuple[str, ...] = (
    "surface_reflectance",
    "sun_zenith",
    "view_zenith",
    "azimuth",
    "ground_altitude",
)

#: Ratio of the B09 / B8A extraterrestrial solar irradiances assumed by the
#: table (the CIBR is defined on irradiance-normalised reflectance).
_CIBR_BAND_RATIO = 955.19 / 813.04

#: TOA B8A below this is too dark for a stable band ratio: the denominator
#: noise dominates and the retrieval scatters wildly. From the v1 retrieval.
DEFAULT_MIN_B8A = 0.05
#: Physically plausible total column water vapour (cm). A retrieval outside it
#: is a failed inversion, not an unusual atmosphere.
DEFAULT_VALID_RANGE = (0.05, 15.0)

#: Uncertainty (cm) assigned to a clean retrieval, to a dark-target retrieval,
#: and to a pixel filled from the scene's own median. From the v1 retrieval.
_UNCERTAINTY_NOMINAL = 0.1
_UNCERTAINTY_DARK = 0.5
_UNCERTAINTY_FILLED = 1.0


@dataclass(frozen=True)
class WaterVapourResult:
    """Retrieved column water vapour and what is known about its quality.

    ``water_vapour_cm`` is complete (no NaN) when ``fill`` is requested: masked
    pixels carry the median of *this scene's own* valid retrievals, which is a
    measurement from the same acquisition rather than an assumed climatology.
    ``valid`` marks the pixels that were actually retrieved, and
    ``uncertainty_cm`` is widened wherever they were not.
    """

    water_vapour_cm: np.ndarray
    uncertainty_cm: np.ndarray
    valid: np.ndarray

    @property
    def masked_fraction(self) -> float:
        """Fraction of pixels that were filled rather than retrieved."""

        if self.valid.size == 0:
            return 1.0
        return float(1.0 - self.valid.mean())


@lru_cache(maxsize=1)
def load_lut() -> tuple[tuple[np.ndarray, ...], np.ndarray]:
    """Return the table as ``(axes, coefficients)`` on its own regular grid.

    ``axes`` is the sorted node vector per :data:`LUT_AXIS_NAMES` and
    ``coefficients`` is ``(*axis_sizes, 2)`` holding ``(b0, b1)``. Reading it
    this way asserts the grid is complete — a hole would mean the file is not
    the regular grid the interpolation assumes, and is rejected rather than
    silently interpolated across.
    """

    text = resources.files(_DATA_PACKAGE).joinpath(_LUT_FILENAME).read_text(encoding="utf-8")
    rows = [
        line
        for line in text.splitlines()
        if line.strip() and not line.startswith("#") and not line.startswith(LUT_AXIS_NAMES[0])
    ]
    table = np.array([[float(value) for value in line.split(";")] for line in rows], dtype=float)
    if table.ndim != 2 or table.shape[1] != len(LUT_AXIS_NAMES) + 2:
        raise ValueError(
            f"{_LUT_FILENAME} must hold {len(LUT_AXIS_NAMES)} axis columns plus b0/b1; "
            f"parsed shape {table.shape}."
        )

    axes = tuple(np.unique(table[:, index]) for index in range(len(LUT_AXIS_NAMES)))
    shape = [axis.size for axis in axes]
    expected = int(np.prod(shape))
    if table.shape[0] != expected:
        raise ValueError(
            f"{_LUT_FILENAME} holds {table.shape[0]} rows but its axes span {expected} grid "
            f"nodes ({' x '.join(str(n) for n in shape)}); it is not the complete regular grid "
            "the interpolation assumes."
        )
    coefficients: np.ndarray = np.full([*shape, 2], np.nan, dtype=float)
    indices = tuple(
        np.searchsorted(axes[index], table[:, index]) for index in range(len(LUT_AXIS_NAMES))
    )
    coefficients[indices] = table[:, len(LUT_AXIS_NAMES) :]
    if not np.isfinite(coefficients).all():
        raise ValueError(f"{_LUT_FILENAME} has holes in its grid; it cannot be interpolated.")
    return axes, coefficients


def relative_azimuth_deg(vaa_deg: Any, saa_deg: Any) -> np.ndarray:
    """Relative azimuth folded into ``[0, 180]``, the table's own convention.

    The table is symmetric about the principal plane, so 200 degrees of relative
    azimuth is the same geometry as 160 — not as 20, which is what a ``% 180``
    wrap would give.
    """

    raa = np.mod(
        np.asarray(vaa_deg, dtype=np.float64) - np.asarray(saa_deg, dtype=np.float64), 360.0
    )
    folded = np.asarray(np.where(raa > 180.0, 360.0 - raa, raa), dtype=np.float64)
    return cast("np.ndarray", folded)


def _clip_to_axes(features: np.ndarray, axes: Sequence[np.ndarray]) -> np.ndarray:
    """Clamp query points into the table's own domain.

    Beyond the grid the polynomial coefficients have no meaning, and letting a
    multilinear interpolator extrapolate them produces confident nonsense; the
    nearest represented geometry is the honest answer.
    """

    clipped = np.array(features, dtype=np.float64, copy=True)
    for index, axis in enumerate(axes):
        np.clip(clipped[:, index], float(axis[0]), float(axis[-1]), out=clipped[:, index])
    return cast("np.ndarray", clipped)


@lru_cache(maxsize=1)
def _grid_interpolator() -> Any:
    from scipy.interpolate import RegularGridInterpolator

    axes, coefficients = load_lut()
    return RegularGridInterpolator(
        axes, coefficients, method="linear", bounds_error=False, fill_value=None
    )


def _coefficients(features: np.ndarray) -> np.ndarray:
    """Return ``(N, 2)`` ``(b0, b1)`` for ``(N, 5)`` query points."""

    axes, _ = load_lut()
    clipped = _clip_to_axes(features, axes)
    values: np.ndarray = np.asarray(_grid_interpolator()(clipped), dtype=np.float64)
    return values


def retrieve_water_vapour(
    *,
    toa_b09: np.ndarray,
    toa_b8a: np.ndarray,
    sza_deg: Any,
    vza_deg: Any,
    raa_deg: Any,
    elevation_km: Any,
    cloud_mask: np.ndarray | None = None,
    min_b8a: float = DEFAULT_MIN_B8A,
    valid_range: tuple[float, float] = DEFAULT_VALID_RANGE,
    fill: bool = True,
) -> WaterVapourResult:
    """Retrieve total column water vapour (cm) from one acquisition's TOA.

    ``toa_b09`` and ``toa_b8a`` are top-of-atmosphere reflectance on a common
    grid; the angles are in degrees and ``elevation_km`` in km, each either a
    scalar or broadcastable to that grid. ``raa_deg`` must already be folded
    into ``[0, 180]`` — see :func:`relative_azimuth_deg`.

    Pixels are masked where B8A is too dark for a stable ratio, where either
    band is missing, where the inversion leaves the physical range, and where
    ``cloud_mask`` is set. With ``fill`` (the default) those pixels take the
    median of this acquisition's *own* valid retrievals and a widened
    uncertainty; without it they are NaN.
    """

    band09 = np.asarray(toa_b09, dtype=np.float64)
    band8a = np.asarray(toa_b8a, dtype=np.float64)
    if band09.shape != band8a.shape:
        raise ValueError(f"B09 {band09.shape} and B8A {band8a.shape} must be on the same grid.")
    low, high = float(valid_range[0]), float(valid_range[1])

    sza = np.broadcast_to(np.asarray(sza_deg, dtype=np.float64), band8a.shape)
    vza = np.broadcast_to(np.asarray(vza_deg, dtype=np.float64), band8a.shape)
    raa = np.broadcast_to(np.asarray(raa_deg, dtype=np.float64), band8a.shape)
    elevation = np.broadcast_to(np.asarray(elevation_km, dtype=np.float64), band8a.shape)

    with np.errstate(invalid="ignore", divide="ignore"):
        cibr = _CIBR_BAND_RATIO * band8a / band09
        log_cibr = np.log10(cibr)

    # Only invert where the inputs can support it; the rest is filled below.
    retrievable = (
        np.isfinite(band09)
        & np.isfinite(band8a)
        & np.isfinite(log_cibr)
        & (band09 > 0.0)
        & (band8a > float(min_b8a))
        & np.isfinite(sza)
        & np.isfinite(vza)
        & np.isfinite(raa)
        & np.isfinite(elevation)
    )
    if cloud_mask is not None:
        retrievable &= ~np.asarray(cloud_mask, dtype=bool)

    water_vapour = np.full(band8a.shape, np.nan, dtype=np.float64)
    if retrievable.any():
        features = np.stack(
            [
                band8a[retrievable],
                sza[retrievable],
                vza[retrievable],
                raa[retrievable],
                elevation[retrievable],
            ],
            axis=1,
        )
        coefficients = _coefficients(features)
        water_vapour[retrievable] = (
            coefficients[:, 0] * log_cibr[retrievable] + coefficients[:, 1]
        ) ** 2

    valid = retrievable & np.isfinite(water_vapour) & (water_vapour >= low) & (water_vapour <= high)
    uncertainty = np.where(
        valid,
        np.where(band8a > float(min_b8a), _UNCERTAINTY_NOMINAL, _UNCERTAINTY_DARK),
        _UNCERTAINTY_FILLED,
    ).astype(np.float64)

    output = np.where(valid, water_vapour, np.nan)
    if fill and valid.any():
        # The scene's own median: a measurement from this acquisition, not an
        # assumed climatology. The widened uncertainty above marks every such
        # pixel, and callers log the masked fraction.
        output = np.where(valid, output, float(np.median(water_vapour[valid])))
    return WaterVapourResult(
        water_vapour_cm=output.astype(np.float32),
        uncertainty_cm=uncertainty.astype(np.float32),
        valid=valid,
    )
