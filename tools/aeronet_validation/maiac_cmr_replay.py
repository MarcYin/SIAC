#!/usr/bin/env python3
"""Which MCD19A2 granules the lineage used, replayed from the public catalogue.

The AERONET-lineage priors were built from whatever NASA's CMR catalogue
returned to the lineage's own search: its temporal semantics decide which days
are in (a +/-2 day teacher window around 02:19 UTC returned four days, not the
five a calendar overlap suggests), and its ``count`` caps silently truncate
(64 granules a month cut a polar AOI's month at day 16). Deciding the granule set
any other way reproduces neither.

CMR search is public -- only data downloads authenticate -- so this replays the
lineage's exact search through its own ``EarthAccessSource.search_granules``
without ever logging in to Earthdata, and refuses to return if a login happened.
Earth Engine then supplies the data for exactly these granules, matched down to
the production version CMR named.
"""

from __future__ import annotations

import datetime as dt
from dataclasses import dataclass


@dataclass(frozen=True, order=True)
class GranuleId:
    """One MCD19A2 tile-day, as CMR names it: ``MCD19A2.A2024027.h28v05.061.<production>``."""

    day: str  # e.g. "A2024027"
    tile: str  # e.g. "h28v05"
    production: str  # 13-digit production time, e.g. "2024030172711"

    @classmethod
    def from_granule_ur(cls, granule_ur: str) -> GranuleId:
        parts = str(granule_ur).split(".")
        if len(parts) < 5 or not parts[0].startswith("MCD19A2"):
            raise ValueError(f"unexpected MCD19A2 granule id {granule_ur!r}")
        return cls(day=parts[1], tile=parts[2], production=parts[4])

    @property
    def date(self) -> dt.date:
        return dt.datetime.strptime(self.day[1:], "%Y%j").date()

    @property
    def filename(self) -> str:
        return f"MCD19A2.{self.day}.{self.tile}.061.{self.production}.hdf"


#: The providers' search sources are built only to reuse their catalogue, short
#: name and provider settings; nothing is downloaded, so the cache path is inert.
_CACHE = "/gws/ssde/j25a/nceo_isp/public/siac_refactor/cache/maiac_day_aod"


def _refuse_login() -> None:
    raise RuntimeError("the catalogue replay must never log in to Earthdata")


def _search(provider, short_name, bounds, crs, temporal, count) -> list[GranuleId]:
    # earthaccess searches CMR anonymously; the provider's source would only log
    # in on a download. Make any login attempt through it fail outright.
    provider.source._ensure_auth = _refuse_login
    granules = provider.source.search_granules(
        short_name=short_name,
        bounds=bounds,
        crs=crs,
        temporal=temporal,
        provider=provider.provider,
        count=count,
    )
    if provider.source.is_authenticated:
        raise RuntimeError("the catalogue replay authenticated to Earthdata; it must never log in")
    return [GranuleId.from_granule_ur(g["umm"]["GranuleUR"]) for g in granules]


def teacher_granules(
    bounds: tuple[float, float, float, float],
    crs: str,
    obs_time: dt.datetime,
    *,
    window_days: int = 2,
    max_granules: int = 8,
) -> list[GranuleId]:
    """The granules ``MCD19AODProvider.get_prior`` searched for, in CMR order."""
    from siac.adapters.atmo.mcd19_earthaccess import MCD19AODProvider
    from siac.adapters.data.earthaccess_source import EarthAccessSource

    provider = MCD19AODProvider(cache_dir=_CACHE, best_quality_qa=True, allow_default_prior=False)
    short_name = provider.short_name or provider.catalog.resolve_short_name(
        provider.product_keys[0]
    )
    temporal = EarthAccessSource.temporal_window(obs_time, window_days)
    return _search(provider, short_name, bounds, crs, temporal, max_granules)


def month_granules(
    bounds: tuple[float, float, float, float],
    crs: str,
    year: int,
    month: int,
    *,
    max_granules: int = 64,
) -> list[GranuleId]:
    """The granules ``MAIACDayAODProvider`` searched for one calendar month."""
    from siac.adapters.atmo import maiac_day_aod

    provider = maiac_day_aod.MAIACDayAODProvider(cache_dir=_CACHE)
    short_name = provider.catalog.resolve_short_name(maiac_day_aod._MCD19_PRODUCT_KEY)
    temporal = maiac_day_aod._month_temporal_range(int(year), int(month))
    return _search(provider, short_name, bounds, crs, temporal, max_granules)
