"""Google Earth Engine BRDF adapter placeholder."""

from __future__ import annotations

from typing import Any


class GEEBRDFProvider:
    """Placeholder provider until a real GEE adapter is implemented."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        _ = (args, kwargs)

    def __call__(self, *args: Any, **kwargs: Any) -> None:
        raise NotImplementedError("GEE BRDF provider is not implemented in this package.")
