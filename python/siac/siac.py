"""Compatibility shim for the public SIAC API."""

from siac.api.public import (
    SIAC,
    process_landsat8,
    process_sentinel2,
    resolve_s2_input,
    search_sentinel2,
    siac_process,
    siac_process_s2,
)

__all__ = [
    "SIAC",
    "process_sentinel2",
    "process_landsat8",
    "resolve_s2_input",
    "search_sentinel2",
    "siac_process_s2",
    "siac_process",
]
