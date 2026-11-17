"""amchi."""

__version__ = "0.0.1"

from ._amchi import from_geometry, from_ts_geometries, string
from ._candidates import add_bond_candidate, candidates, remove_candidate
from ._connect import connectivity, ts_connectivity
from ._element import AmchiError, UnsupportedElementError
from ._parity import ParityMargins, parity_margins
from ._types import Bond, Candidates, Connectivity, GeometryLike, IllDefinedParityError

__all__ = [
    "AmchiError",
    "Bond",
    "Candidates",
    "Connectivity",
    "GeometryLike",
    "IllDefinedParityError",
    "ParityMargins",
    "UnsupportedElementError",
    "add_bond_candidate",
    "candidates",
    "connectivity",
    "from_geometry",
    "from_ts_geometries",
    "parity_margins",
    "remove_candidate",
    "string",
    "ts_connectivity",
]
