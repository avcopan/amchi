"""amchi."""

__version__ = "0.0.0"

from ._connect import connectivity
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
    "connectivity",
    "parity_margins",
]
