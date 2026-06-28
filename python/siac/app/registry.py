"""Small registries used by SIAC assembly."""

from __future__ import annotations

from typing import TYPE_CHECKING, Generic, TypeVar

if TYPE_CHECKING:
    from collections.abc import Callable

T = TypeVar("T")


class NamedRegistry(Generic[T]):
    """Dictionary-backed named factory registry."""

    def __init__(self, label: str):
        self._label = label
        self._factories: dict[str, Callable[..., T]] = {}

    def register(self, name: str) -> Callable[[Callable[..., T]], Callable[..., T]]:
        key = name.lower()

        def _decorator(factory: Callable[..., T]) -> Callable[..., T]:
            self._factories[key] = factory
            return factory

        return _decorator

    def get(self, name: str) -> Callable[..., T]:
        key = name.lower()
        if key not in self._factories:
            raise ValueError(f"Unknown {self._label}: {name!r}")
        return self._factories[key]

    def names(self) -> tuple[str, ...]:
        return tuple(sorted(self._factories))


ATMO_PROVIDER_REGISTRY: NamedRegistry[object] = NamedRegistry("atmo provider")
BRDF_PROVIDER_REGISTRY: NamedRegistry[object] = NamedRegistry("BRDF provider")
MONTHLY_COMPOSITE_PROVIDER_REGISTRY: NamedRegistry[object] = NamedRegistry(
    "monthly composite provider"
)
S2_BACKEND_REGISTRY: NamedRegistry[object] = NamedRegistry("S2 backend")
SURFACE_PRIOR_METHOD_REGISTRY: NamedRegistry[object] = NamedRegistry("surface prior method")
SOLVER_METHOD_REGISTRY: NamedRegistry[object] = NamedRegistry("solver method")
