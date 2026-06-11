"""Shared helpers for the direct (engine-driven) RT backends."""

from __future__ import annotations

from typing import Any

from siac.runtime import RTCoefficients

_REQUIRED_COEFFICIENTS = ("xap", "xbp", "xcp")


def coerce_rt_coefficients(outputs: Any, *, runner_name: str) -> RTCoefficients:
    """Coerce a runner result into :class:`RTCoefficients`.

    Runners may return a ready :class:`RTCoefficients` or a mapping with the
    required ``xap``/``xbp``/``xcp`` arrays plus optional extras. This is the
    single definition of that contract for every direct backend (6S,
    libRadtran); ``runner_name`` only flavours the error messages.
    """
    if isinstance(outputs, RTCoefficients):
        return outputs
    if not isinstance(outputs, dict):
        raise TypeError(
            f"{runner_name} runner must return RTCoefficients or dict[str, xr.DataArray], "
            f"got {type(outputs).__name__}"
        )
    missing = set(_REQUIRED_COEFFICIENTS) - set(outputs)
    if missing:
        raise KeyError(
            f"{runner_name} runner did not return required coefficients: "
            + ", ".join(sorted(missing))
        )
    extras = {
        name: value for name, value in outputs.items() if name not in set(_REQUIRED_COEFFICIENTS)
    }
    return RTCoefficients(
        xap=outputs["xap"],
        xbp=outputs["xbp"],
        xcp=outputs["xcp"],
        extras=extras,
    )
