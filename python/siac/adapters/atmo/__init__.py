"""Atmospheric prior providers (CAMS, MERRA-2, ERA5)."""

from siac.adapters.atmo.cams import CAMSProvider as CAMSProvider
from siac.adapters.atmo.mcd19_earthaccess import CachedMCD19AODProvider as CachedMCD19AODProvider
from siac.adapters.atmo.mcd19_earthaccess import MCD19AODProvider as MCD19AODProvider
from siac.adapters.atmo.mcd19_earthaccess import VNP19AODProvider as VNP19AODProvider
from siac.adapters.atmo.merra2 import MERRA2Provider as MERRA2Provider
