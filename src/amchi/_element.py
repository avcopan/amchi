"""Valences, lone pairs and covalent radii of the supported elements.

The values are the old AutoMol's: a valence and a lone-pair count for each main group
(neutral atoms with standard valences), and the covalent radii of Alvarez (2008).
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


# The main-group elements from H to Xe
_ELEMENTS: dict[str, _Element] = {
    "H": _Element(1, 0, 0.31),
    "He": _Element(0, 4, 0.28),
    "Li": _Element(1, 0, 1.28),
    "Be": _Element(2, 0, 0.96),
    "B": _Element(3, 0, 0.84),
    "C": _Element(4, 0, 0.76),
    "N": _Element(3, 1, 0.71),
    "O": _Element(2, 2, 0.66),
    "F": _Element(1, 3, 0.57),
    "Ne": _Element(0, 4, 0.58),
    "Na": _Element(1, 0, 1.66),
    "Mg": _Element(2, 0, 1.41),
    "Al": _Element(3, 0, 1.21),
    "Si": _Element(4, 0, 1.11),
    "P": _Element(3, 1, 1.07),
    "S": _Element(2, 2, 1.05),
    "Cl": _Element(1, 3, 1.02),
    "Ar": _Element(0, 4, 1.06),
    "K": _Element(1, 0, 2.03),
    "Ca": _Element(2, 0, 1.76),
    "Ga": _Element(3, 0, 1.22),
    "Ge": _Element(4, 0, 1.20),
    "As": _Element(3, 1, 1.19),
    "Se": _Element(2, 2, 1.20),
    "Br": _Element(1, 3, 1.20),
    "Kr": _Element(0, 4, 1.16),
    "Rb": _Element(1, 0, 2.20),
    "Sr": _Element(2, 0, 1.95),
    "In": _Element(3, 0, 1.42),
    "Sn": _Element(4, 0, 1.39),
    "Sb": _Element(3, 1, 1.39),
    "Te": _Element(2, 2, 1.38),
    "I": _Element(1, 3, 1.39),
    "Xe": _Element(0, 4, 1.40),
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
