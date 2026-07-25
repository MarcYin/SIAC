"""Radiative-transfer space identity for atmospherically corrected products.

A surface-reflectance product that was produced by atmospheric correction only
means something relative to the radiative-transfer model used to correct it: the
RT backend and, critically, its aerosol model. Correcting a library with one RT
model and then solving against it with another injects a systematic offset —
measured at roughly 0.008 reflectance in the blue between a continental-average
LUT and native 6S with the exact CCI climatology mixture, which is large enough
to bias retrieved AOT by a wide margin because surface error amplifies into AOT
error by one to two orders of magnitude.

:class:`RTSpace` makes that identity explicit so a mismatch between the space a
surface prior was corrected in and the space the solver runs in is rejected when
the execution plan is built, rather than silently degrading the retrieval.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

__all__ = ["RTSpace", "describe_aerosol_model"]

#: Fraction rounding used when canonicalising mixture components. Two decimals
#: keeps genuinely different mixtures distinct while tolerating float noise from
#: repeated normalisation.
_FRACTION_PRECISION = 2


def _normalized_fractions(values: dict[str, float]) -> str:
    total = sum(float(v) for v in values.values())
    if total <= 0.0:
        return "empty"
    parts = [
        f"{name}={float(value) / total:.{_FRACTION_PRECISION}f}"
        for name, value in sorted(values.items())
        if float(value) > 0.0
    ]
    return ",".join(parts) if parts else "empty"


def describe_aerosol_model(aerosol: Any) -> str:
    """Return a canonical, comparable descriptor for an aerosol setup.

    ``aerosol`` is an ``RTAerosolSetupConfig`` (or ``None``). The descriptor
    captures everything that changes the RT coefficients: the named profile and,
    for mixtures, the component fractions. Two setups that produce the same
    coefficients produce the same descriptor.
    """

    if aerosol is None:
        return "default"

    profile = getattr(aerosol, "profile", None) or "default"

    distribution = getattr(aerosol, "distribution", None)
    if distribution is not None:
        components = getattr(distribution, "components", None)
        if components:
            # Components carry no name; their log-normal mode radius identifies
            # the species (CCI dust/sea-salt/soot/water-soluble each have their
            # own rmean), so key the fractions by it.
            fractions = {
                f"r{float(getattr(component, 'rmean', 0.0)):.3f}": float(
                    getattr(component, "percentage_density", 0.0)
                )
                for component in components
            }
            return f"{profile}[{_normalized_fractions(fractions)}]"

    mixture = getattr(aerosol, "mixture", None)
    if mixture is not None:
        names = ("dust", "water_soluble", "oceanic", "soot")
        fractions = {name: float(value) for name, value in zip(names, mixture, strict=False)}
        return f"{profile}[{_normalized_fractions(fractions)}]"

    model_path = getattr(aerosol, "model_path", None)
    if model_path is not None:
        return f"{profile}[{model_path}]"

    return str(profile)


@dataclass(frozen=True)
class RTSpace:
    """Identity of the RT model a reflectance product was corrected in.

    ``backend`` is the RT backend name (``sixs``, ``lut``, ``libradtran``,
    ``emulator``) and ``aerosol`` is the canonical aerosol-model descriptor from
    :func:`describe_aerosol_model`.
    """

    backend: str
    aerosol: str

    def __str__(self) -> str:
        return f"{self.backend}/{self.aerosol}"

    @classmethod
    def from_setup(cls, backend: str, setup: Any) -> RTSpace:
        """Build an :class:`RTSpace` from a backend name and an ``RTSetupConfig``."""

        return cls(
            backend=str(backend),
            aerosol=describe_aerosol_model(getattr(setup, "aerosol", None)),
        )

    @classmethod
    def for_solver(cls, rt_model: Any, solver_config: Any) -> RTSpace | None:
        """Return the RT space the solver *effectively* runs in.

        Scene-adaptive aerosol modes resolve their mixture from the scene's
        location and month inside the solver, so the backend handed to the prior
        still carries the base aerosol setup. For those modes the identity is the
        rule (``cci_climatology_exact``) rather than one scene's resolved
        fractions — a library built under the same rule is consistent with the
        solve even though the numbers differ from scene to scene.
        """

        backend = getattr(rt_model, "backend_name", None)
        if backend is None:
            return None
        mode = str(getattr(solver_config, "surface_driven_aerosol_species", "none") or "none")
        if mode != "none":
            return cls(backend=str(backend), aerosol=mode)
        return cls.from_rt_model(rt_model)

    @classmethod
    def from_rt_model(cls, rt_model: Any) -> RTSpace | None:
        """Build an :class:`RTSpace` from a live RT backend instance.

        Returns ``None`` when the backend does not expose enough information to
        identify its space (in which case consistency cannot be checked).
        """

        backend = getattr(rt_model, "backend_name", None)
        if backend is None:
            return None
        setup = getattr(rt_model, "rt_setup", None)
        return cls.from_setup(str(backend), setup)

    def matches(self, other: RTSpace) -> bool:
        """Return whether two spaces are the same RT model."""

        return self.backend == other.backend and self.aerosol == other.aerosol
