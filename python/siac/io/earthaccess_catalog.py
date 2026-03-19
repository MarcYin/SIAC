"""Earthaccess dataset catalog and validation utilities."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from siac.errors import DataNotFoundError
from siac.io.earthaccess_source import EarthAccessSource


@dataclass
class EarthaccessProduct:
    """Catalog entry describing one logical Earthdata product family."""

    key: str
    description: str
    keyword: str
    short_name: str | None = None
    provider: str | None = None


@dataclass
class ProductValidationResult:
    """Validation result for one Earthaccess catalog entry."""

    key: str
    dataset_count: int
    discovered_short_names: tuple[str, ...]
    selected_short_name: str


def _extract_short_name(dataset: Any) -> str | None:
    """Best-effort short-name extraction from earthaccess dataset objects."""
    if isinstance(dataset, dict):
        if isinstance(dataset.get("short_name"), str):
            return dataset["short_name"]
        umm = dataset.get("umm")
        if isinstance(umm, dict) and isinstance(umm.get("ShortName"), str):
            return umm["ShortName"]

    short_name = getattr(dataset, "short_name", None)
    if isinstance(short_name, str):
        return short_name

    umm_obj = getattr(dataset, "umm", None)
    if isinstance(umm_obj, dict):
        umm_short = umm_obj.get("ShortName")
        if isinstance(umm_short, str):
            return umm_short

    return None


def default_products() -> dict[str, EarthaccessProduct]:
    """Return the default logical product registry for SIAC Earthaccess use."""
    return {
        "landsat": EarthaccessProduct(
            key="landsat",
            description="Landsat Collection products available through Earthdata",
            keyword="Landsat Collection",
            short_name=None,
        ),
        "merra2_atmo": EarthaccessProduct(
            key="merra2_atmo",
            description="MERRA-2 atmospheric products",
            keyword="MERRA-2",
            short_name=None,
        ),
        "mcd19_aod": EarthaccessProduct(
            key="mcd19_aod",
            description="MODIS MCD19 atmospheric aerosol optical depth",
            keyword="MCD19",
            short_name="MCD19A2",
        ),
        "vnp19_aod": EarthaccessProduct(
            key="vnp19_aod",
            description="VIIRS VNP19 atmospheric aerosol optical depth",
            keyword="VNP19",
            short_name="VNP19A2",
        ),
        "vj119_aod": EarthaccessProduct(
            key="vj119_aod",
            description="VIIRS NOAA-20 VJ119 atmospheric aerosol optical depth",
            keyword="VJ119",
            short_name="VJ119A2",
        ),
        "vj219_aod": EarthaccessProduct(
            key="vj219_aod",
            description="VIIRS NOAA-21 VJ219 atmospheric aerosol optical depth",
            keyword="VJ219",
            short_name="VJ219A2",
        ),
        "mcd43_brdf": EarthaccessProduct(
            key="mcd43_brdf",
            description="MODIS MCD43 BRDF/albedo products",
            keyword="MCD43",
            short_name="MCD43A1",
        ),
        "mcd19_brdf": EarthaccessProduct(
            key="mcd19_brdf",
            description="MODIS MCD19 MAIAC BRDF kernel products",
            keyword="MCD19",
            short_name="MCD19A3",
        ),
        "vnp43_brdf": EarthaccessProduct(
            key="vnp43_brdf",
            description="VIIRS VNP43 BRDF/albedo products",
            keyword="VNP43",
            short_name="VNP43MA1",
        ),
    }


class EarthAccessCatalog:
    """Registry and validation for Earthaccess-backed logical product keys."""

    def __init__(
        self,
        source: EarthAccessSource | None = None,
        products: dict[str, EarthaccessProduct] | None = None,
    ) -> None:
        self.source = source or EarthAccessSource()
        self._products = dict(products or default_products())

    def keys(self) -> tuple[str, ...]:
        """Return available logical product keys in deterministic order."""
        return tuple(sorted(self._products))

    def get_product(self, key: str) -> EarthaccessProduct:
        """Return catalog entry by key."""
        try:
            return self._products[key]
        except KeyError as exc:
            raise KeyError(f"Unknown Earthaccess product key: {key!r}") from exc

    def set_override(
        self,
        key: str,
        *,
        short_name: str | None = None,
        provider: str | None = None,
    ) -> None:
        """Override entry fields from runtime/config inputs."""
        product = self.get_product(key)
        if short_name is not None:
            product.short_name = short_name
        if provider is not None:
            product.provider = provider

    def validate_product(self, key: str, *, count: int = 20) -> ProductValidationResult:
        """Validate one product key against Earthdata dataset discovery."""
        product = self.get_product(key)

        # First attempt: resolve by explicit short_name if provided.
        datasets: list[Any] = []
        if product.short_name is not None:
            datasets = self.source.search_datasets(
                short_name=product.short_name,
                provider=product.provider,
                count=count,
            )

        # Fallback: keyword search.
        if not datasets:
            datasets = self.source.search_datasets(
                keyword=product.keyword,
                provider=product.provider,
                count=count,
            )

        if not datasets:
            raise DataNotFoundError(
                f"No Earthdata datasets found for key={key!r} keyword={product.keyword!r}"
            )

        short_names = sorted(
            {
                short_name
                for short_name in (_extract_short_name(item) for item in datasets)
                if short_name is not None
            }
        )

        if product.short_name is not None:
            selected = product.short_name
            if short_names and selected not in short_names:
                raise DataNotFoundError(
                    f"Configured short_name={selected!r} not present for key={key!r}. "
                    f"Discovered={short_names}"
                )
        else:
            if not short_names:
                raise DataNotFoundError(
                    f"Could not infer dataset short_name for key={key!r}; "
                    "set an explicit short_name override."
                )
            selected = short_names[0]
            product.short_name = selected

        return ProductValidationResult(
            key=key,
            dataset_count=len(datasets),
            discovered_short_names=tuple(short_names),
            selected_short_name=selected,
        )

    def validate_all(
        self,
        keys: list[str] | None = None,
        *,
        count: int = 20,
    ) -> dict[str, ProductValidationResult]:
        """Validate all (or selected) product keys."""
        selected_keys = keys if keys is not None else list(self.keys())
        return {k: self.validate_product(k, count=count) for k in selected_keys}

    def resolve_short_name(self, key: str) -> str:
        """Return short_name for key, validating/discovering it when needed."""
        product = self.get_product(key)
        if product.short_name is not None:
            return product.short_name
        return self.validate_product(key).selected_short_name
