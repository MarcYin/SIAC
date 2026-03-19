"""Google Earth Engine BRDF adapter placeholder."""


class GEEBRDFProvider:
    """Placeholder provider until a real GEE adapter is implemented."""

    def __init__(self, *args, **kwargs):  # noqa: D401, ANN002, ANN003
        _ = (args, kwargs)

    def __call__(self, *args, **kwargs):  # noqa: ANN002, ANN003
        raise NotImplementedError("GEE BRDF provider is not implemented in this package.")
