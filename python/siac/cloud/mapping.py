"""Class mapping utilities for cloud/cloud-shadow masks."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import xarray as xr

EXPECTED_CLASSES = (0, 1, 2, 3)


def validate_class_mapping(mapping: dict[int, Iterable[int]] | None) -> dict[int, list[int]]:
    """Validate class mapping and ensure no source class maps to multiple targets."""
    if mapping is None:
        return {}
    if not isinstance(mapping, dict):
        raise TypeError("class_mapping must be a dictionary if provided")

    normalized: dict[int, list[int]] = {}
    reverse: dict[int, int] = {}

    for target_raw, source_values in mapping.items():
        target = int(target_raw)
        if target not in EXPECTED_CLASSES:
            raise ValueError(
                f"class_mapping target class must be one of {EXPECTED_CLASSES}, got {target}"
            )

        if isinstance(source_values, (int, np.integer)):
            sources = [int(source_values)]
        elif isinstance(source_values, Iterable):
            sources = [int(v) for v in source_values]
        else:
            raise TypeError(
                f"class_mapping[{target}] must be an int or iterable of ints"
            )

        uniq = sorted(set(sources))
        normalized[target] = uniq

        for source in uniq:
            if source in reverse and reverse[source] != target:
                raise ValueError(
                    "Invalid class_mapping: one source class maps to multiple "
                    f"target classes ({source} -> {reverse[source]} and {target})"
                )
            reverse[source] = target

    return normalized


def _ensure_expected_classes(out: xr.DataArray) -> xr.DataArray:
    values = out.values
    if values.size == 0:
        return out
    uniques = np.unique(values)
    if not set(int(v) for v in uniques).issubset(set(EXPECTED_CLASSES)):
        raise ValueError(
            f"Mapped cloud class values must be in {EXPECTED_CLASSES}, got {uniques.tolist()}"
        )
    return out


def apply_class_mapping(
    data: xr.DataArray,
    mapping: dict[int, Iterable[int]] | None = None,
    *,
    unmapped_to_missing: bool = True,
    output_name: str = "cloud_classes",
) -> xr.DataArray:
    """Map arbitrary class labels to SIAC standardized classes {0,1,2,3}."""
    normalized = validate_class_mapping(mapping)

    src = np.asarray(data.values)
    finite = np.isfinite(src)
    out = np.zeros(src.shape, dtype=np.uint8)

    if not normalized:
        # Identity mapping for already-standard classes.
        for klass in EXPECTED_CLASSES:
            out[(src == klass) & finite] = np.uint8(klass)
        known = np.isin(src, np.array(EXPECTED_CLASSES, dtype=src.dtype)) & finite
        if not unmapped_to_missing and not bool(np.all(known | ~finite)):
            unknown = np.unique(src[(~known) & finite])
            raise ValueError(
                f"Input cloud classes include unmapped values: {unknown.tolist()}"
            )
    else:
        all_sources: set[int] = set()
        for target, sources in normalized.items():
            if not sources:
                continue
            all_sources.update(sources)
            out[np.isin(src, sources) & finite] = np.uint8(target)

        if not unmapped_to_missing:
            known = np.isin(src, np.array(sorted(all_sources), dtype=src.dtype)) & finite
            if not bool(np.all(known | ~finite)):
                unknown = np.unique(src[(~known) & finite])
                raise ValueError(
                    f"Input cloud classes include unmapped values: {unknown.tolist()}"
                )

    result = xr.DataArray(
        out,
        dims=data.dims,
        coords=data.coords,
        attrs=dict(data.attrs),
        name=output_name,
    )
    result.attrs["class_values"] = {0: "missing", 1: "clear", 2: "cloud", 3: "shadow"}
    return _ensure_expected_classes(result)
