"""Tests for the network-free staged MCD19 MAIAC provider."""

from __future__ import annotations

import logging
from datetime import datetime
from typing import TYPE_CHECKING

from siac.adapters.atmo.mcd19_earthaccess import CachedMCD19AODProvider

if TYPE_CHECKING:
    from pathlib import Path


def test_cached_provider_uses_only_local_same_window_granules(tmp_path: Path, monkeypatch) -> None:
    for stamp in ("2024001", "2024002", "2024003", "2024005"):
        (tmp_path / f"MCD19A2.A{stamp}.h18v08.061.fake.hdf").touch()

    provider = CachedMCD19AODProvider(tmp_path, temporal_window_days=1)
    observed: dict[str, object] = {}
    sentinel = object()

    def select(paths, *_args):  # noqa: ANN001
        observed["paths"] = [path.name for path in paths]
        return list(paths)

    def load(paths, **kwargs):  # noqa: ANN001
        observed["selected"] = [path.name for path in paths]
        observed["kwargs"] = kwargs
        return sentinel

    monkeypatch.setattr(provider, "_select_candidate_paths", select)
    monkeypatch.setattr(provider, "_load_from_granules", load)

    state = provider.get_cached_prior(
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        obs_time=datetime(2024, 1, 2, 10, 0),
        resolution=60.0,
    )

    assert state is sentinel
    assert observed["paths"] == [
        "MCD19A2.A2024001.h18v08.061.fake.hdf",
        "MCD19A2.A2024002.h18v08.061.fake.hdf",
        "MCD19A2.A2024003.h18v08.061.fake.hdf",
    ]
    assert observed["selected"] == observed["paths"]
    assert observed["kwargs"] == {
        "bounds": (0.0, 0.0, 1.0, 1.0),
        "crs": "EPSG:4326",
        "resolution": 60.0,
        "short_name": "MCD19A2",
        "obs_time": datetime(2024, 1, 2, 10, 0),
    }


def test_cached_provider_returns_none_without_matching_granules(tmp_path: Path) -> None:
    provider = CachedMCD19AODProvider(tmp_path)

    state = provider.get_cached_prior(
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        obs_time=datetime(2024, 1, 2, 10, 0),
        resolution=60.0,
    )

    assert state is None


def test_cached_provider_treats_empty_qa_coverage_as_expected(tmp_path: Path, monkeypatch, caplog) -> None:
    """A staged fallback should not report valid-but-empty MAIAC data as a parser failure."""
    path = tmp_path / "MCD19A2.A2024002.h18v08.061.fake.hdf"
    path.touch()
    provider = CachedMCD19AODProvider(tmp_path)

    monkeypatch.setattr(provider, "_select_candidate_paths", lambda *_args: [path])

    def no_qa(*_args, **_kwargs):  # noqa: ANN001
        raise ValueError("MCD19 has no QA-valid AOD after reprojection to the requested AOI")

    monkeypatch.setattr(provider, "_load_from_granules", no_qa)
    with caplog.at_level(logging.INFO):
        state = provider.get_cached_prior(
            bounds=(0.0, 0.0, 1.0, 1.0),
            crs="EPSG:4326",
            obs_time=datetime(2024, 1, 2, 10, 0),
            resolution=60.0,
        )

    assert state is None
    assert "no QA-valid AOD" in caplog.text
    assert "granule parsing failed" not in caplog.text
