"""Shared physical and numerical constants used across SIAC.

REVIEW.md §2.4 flagged a sprawl of magic numbers across the RT, solver, and
correction layers. This module collects the most-cited ones in a single
place so that:

- Each constant has a name, a unit, a citation (where one exists), and a
  short note about what changes if you tune it.
- A reader looking at e.g. ``ATMOSPHERIC_SCALE_HEIGHT_KM`` doesn't have to
  guess whether ``8.5`` in ``backend.py`` is the same physical quantity as
  ``8.5`` in some other module.
- A grep for ``ATMOSPHERIC_SCALE_HEIGHT_KM`` finds every callsite, making
  future re-tuning a single PR.

Constants live here only when they are referenced from more than one
module, or when they correspond to a physical quantity that deserves a
named name. Module-local thresholds (e.g. test cosmetic tolerances) stay
where they are.
"""

from __future__ import annotations

from typing import Final

# ---------------------------------------------------------------------------
# Atmospheric / RT physical constants
# ---------------------------------------------------------------------------

#: Atmospheric Rayleigh-equivalent scale height in km.
#:
#: Used by the LUT-backend altitude correction: ``correction = exp(-h / H)``.
#: 8.5 km is the standard value for dry-air Rayleigh scattering and matches
#: most introductory atmospheric-optics texts. Empirical fits give values
#: between 8.0 and 8.5 km depending on temperature profile assumptions.
#:
#: Source: e.g. Bodhaine et al. 1999, "On Rayleigh optical depth calculations".
ATMOSPHERIC_SCALE_HEIGHT_KM: Final[float] = 8.5

# ---------------------------------------------------------------------------
# Numerical solver constants
# ---------------------------------------------------------------------------

#: Default forward-difference step for Jacobian estimation w.r.t. AOT.
#:
#: Used in ``ZarrLUTBackend._compute_jacobian_numerical``. Smaller values
#: give a more accurate first-order approximation but amplify floating-point
#: noise; larger values average over more LUT-grid curvature. 0.01 strikes
#: a workable balance for AOT in [0.001, 2.5].
DEFAULT_JACOBIAN_DELTA_AOT: Final[float] = 0.01

#: Default forward-difference step for Jacobian estimation w.r.t. TCWV (cm).
#:
#: Used in ``ZarrLUTBackend._compute_jacobian_numerical``. TCWV's natural
#: scale is roughly 10x AOT's, so the step is 10x larger. Same trade-off as
#: ``DEFAULT_JACOBIAN_DELTA_AOT``.
DEFAULT_JACOBIAN_DELTA_TCWV: Final[float] = 0.1

# ---------------------------------------------------------------------------
# Reflectance validity ranges
# ---------------------------------------------------------------------------

#: Lower acceptance bound for valid bottom-of-atmosphere reflectance.
#:
#: Slightly negative values are accepted to allow for residual atmospheric
#: correction noise after the inversion (BOA can dip below zero in deep
#: shadows or over very dark targets such as still water).
BOA_VALID_MIN: Final[float] = -0.05

#: Upper acceptance bound for valid bottom-of-atmosphere reflectance.
#:
#: Conservative ceiling (well above any realistic land cover) used to flag
#: numerical blow-ups in the correction. Real BOA never exceeds 1.0 except
#: in very rare specular cases (active fires, glints).
BOA_VALID_MAX: Final[float] = 1.5

# ---------------------------------------------------------------------------
# Surface-prior fallbacks
# ---------------------------------------------------------------------------

#: Default reflectance to use for pixels with no temporal BRDF data.
#:
#: Conservative mid-range value used by ``brdf_whittaker`` when a pixel
#: has no usable observations across the temporal window. Picked to lie
#: between dark soil (~0.05) and snow (~0.8) so the prior pulls weakly
#: in either direction. The companion uncertainty
#: ``DEFAULT_NO_DATA_BOA_UNC`` is large enough that the data term in the
#: cost function dominates whenever real observations are present.
DEFAULT_NO_DATA_BOA: Final[float] = 0.20

#: Default reflectance uncertainty for pixels with no temporal BRDF data.
#:
#: Wide enough (1-σ ≈ 8% reflectance) that the prior contribution to
#: the cost is small whenever real observations exist; tight enough
#: that the prior still acts as a regulariser for fully-cloudy pixels.
DEFAULT_NO_DATA_BOA_UNC: Final[float] = 0.08

# ---------------------------------------------------------------------------
# Sentinel-2 angle parsing fallbacks
# ---------------------------------------------------------------------------

#: Fallback Sentinel-2 view zenith angle (degrees) when MTD_TL.xml is missing
#: ``Mean_Viewing_Incidence_Angle``. Roughly the median across-swath VZA for
#: an MSI scene; the warning emitted alongside this value tells operators
#: their angle parsing fell through.
DEFAULT_S2_VZA_DEG: Final[float] = 5.0

#: Fallback Sentinel-2 view azimuth angle (degrees) — see
#: ``DEFAULT_S2_VZA_DEG``. ~100° is a plausible descending-orbit value.
DEFAULT_S2_VAA_DEG: Final[float] = 100.0

#: Fallback Sentinel-2 angle-grid value for sun_zenith / sun_azimuth when
#: the per-detector ``Values_List`` is empty. 30° is a defensive
#: mid-latitude / mid-day estimate.
DEFAULT_S2_ANGLE_GRID_DEG: Final[float] = 30.0
