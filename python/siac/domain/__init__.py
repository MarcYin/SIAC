"""Domain-layer types and protocols.

Most types here (sensors, spectral, protocols) are dependency-free pure
Python. The one exception is :class:`siac.domain.aoi.AOI`, which depends
on :mod:`siac.geo` because Area-of-Interest handling fundamentally needs
spatial libraries (rasterio, pyproj). The dependency is one-way —
``siac.geo`` does not import ``siac.domain`` — so the conceptual layering
is intact: ``domain → geo`` (and never the reverse). REVIEW.md §1.4
originally flagged this as a layering inversion; treating AOI as a
spatial type rather than a pure domain type clarifies that it isn't.
"""

from siac.domain.aoi import AOI as AOI
from siac.domain.protocols import (
    AerosolSolver as AerosolSolver,
)
from siac.domain.protocols import (
    AtmosphericPriorProvider as AtmosphericPriorProvider,
)
from siac.domain.protocols import (
    BRDFProductProvider as BRDFProductProvider,
)
from siac.domain.protocols import (
    MonthlyCompositeProvider as MonthlyCompositeProvider,
)
from siac.domain.protocols import (
    OutputWriter as OutputWriter,
)
from siac.domain.protocols import (
    RTModelBackend as RTModelBackend,
)
from siac.domain.protocols import (
    SatellitePreprocessor as SatellitePreprocessor,
)
from siac.domain.protocols import (
    SurfacePriorDeriver as SurfacePriorDeriver,
)
from siac.domain.sensors import SensorBand as SensorBand
from siac.domain.sensors import SensorConfig as SensorConfig
from siac.domain.spectral import RelativeSpectralResponse as RelativeSpectralResponse
