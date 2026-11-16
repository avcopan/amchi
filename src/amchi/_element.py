"""Valences, lone pairs and atomic radii of the supported elements.

The values are the old AutoMol's: a valence and a lone-pair count for each main group
(neutral atoms with standard valences), the covalent radii of Alvarez (2008), and the
van der Waals radii of Mantina et al. (2009), as qcelemental gives them.
"""

from typing import NamedTuple


class AmchiError(Exception):
    """Base class for the errors that amchi raises."""


class UnsupportedElementError(AmchiError):
    """An element that amchi has no data for."""


class _Element(NamedTuple):
    valence: int
    lone_pairs: int
    covalent_radius: float  # Å
    van_der_waals_radius: float  # Å


# The main-group elements from H to Xe
_ELEMENTS: dict[str, _Element] = {
    "H": _Element(1, 0, 0.31, 1.10),
    "He": _Element(0, 4, 0.28, 1.40),
    "Li": _Element(1, 0, 1.28, 1.81),
    "Be": _Element(2, 0, 0.96, 1.53),
    "B": _Element(3, 0, 0.84, 1.92),
    "C": _Element(4, 0, 0.76, 1.70),
    "N": _Element(3, 1, 0.71, 1.55),
    "O": _Element(2, 2, 0.66, 1.52),
    "F": _Element(1, 3, 0.57, 1.47),
    "Ne": _Element(0, 4, 0.58, 1.54),
    "Na": _Element(1, 0, 1.66, 2.27),
    "Mg": _Element(2, 0, 1.41, 1.73),
    "Al": _Element(3, 0, 1.21, 1.84),
    "Si": _Element(4, 0, 1.11, 2.10),
    "P": _Element(3, 1, 1.07, 1.80),
    "S": _Element(2, 2, 1.05, 1.80),
    "Cl": _Element(1, 3, 1.02, 1.75),
    "Ar": _Element(0, 4, 1.06, 1.88),
    "K": _Element(1, 0, 2.03, 2.75),
    "Ca": _Element(2, 0, 1.76, 2.31),
    "Ga": _Element(3, 0, 1.22, 1.87),
    "Ge": _Element(4, 0, 1.20, 2.11),
    "As": _Element(3, 1, 1.19, 1.85),
    "Se": _Element(2, 2, 1.20, 1.90),
    "Br": _Element(1, 3, 1.20, 1.83),
    "Kr": _Element(0, 4, 1.16, 2.02),
    "Rb": _Element(1, 0, 2.20, 3.03),
    "Sr": _Element(2, 0, 1.95, 2.49),
    "In": _Element(3, 0, 1.42, 1.93),
    "Sn": _Element(4, 0, 1.39, 2.17),
    "Sb": _Element(3, 1, 1.39, 2.06),
    "Te": _Element(2, 2, 1.38, 2.06),
    "I": _Element(1, 3, 1.39, 1.98),
    "Xe": _Element(0, 4, 1.40, 2.16),
}


def _element(symbol: str) -> _Element:
    try:
        return _ELEMENTS[symbol]
    except KeyError:
        msg = f"Unsupported element: {symbol!r}"
        raise UnsupportedElementError(msg) from None


def valence(symbol: str) -> int:
    """The number of bonds that a neutral atom of this element forms."""
    return _element(symbol).valence


def lone_pairs(symbol: str) -> int:
    """The number of lone pairs of a neutral atom of this element."""
    return _element(symbol).lone_pairs


def covalent_radius(symbol: str) -> float:
    """The covalent radius in Å."""
    return _element(symbol).covalent_radius


def van_der_waals_radius(symbol: str) -> float:
    """The van der Waals radius in Å."""
    return _element(symbol).van_der_waals_radius
