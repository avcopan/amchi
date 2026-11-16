"""Tests for the element data."""

import pytest

from amchi import AmchiError, UnsupportedElementError
from amchi._element import covalent_radius, lone_pairs, valence


@pytest.mark.parametrize(
    ("symbol", "expected_valence", "expected_lone_pairs", "expected_radius"),
    [
        ("H", 1, 0, 0.31),
        ("C", 4, 0, 0.76),
        ("N", 3, 1, 0.71),
        ("O", 2, 2, 0.66),
        ("S", 2, 2, 1.05),
        ("Cl", 1, 3, 1.02),
    ],
)
def test__element(
    symbol: str,
    expected_valence: int,
    expected_lone_pairs: int,
    expected_radius: float,
) -> None:
    """The old AutoMol's values."""
    assert valence(symbol) == expected_valence
    assert lone_pairs(symbol) == expected_lone_pairs
    assert covalent_radius(symbol) == expected_radius


@pytest.mark.parametrize("symbol", ["Fe", "X", "c", ""])
def test__unsupported_element(symbol: str) -> None:
    """A transition metal, a dummy atom or a misspelled symbol raises a named error."""
    for function in (valence, lone_pairs, covalent_radius):
        with pytest.raises(UnsupportedElementError, match="Unsupported element"):
            function(symbol)
    assert issubclass(UnsupportedElementError, AmchiError)
