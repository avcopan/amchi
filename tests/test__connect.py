"""Tests for the connectivity from a geometry."""

from dataclasses import dataclass

import pytest
from numpy.typing import ArrayLike

import amchi


@dataclass
class Geometry:
    """A minimal GeometryLike."""

    symbols: list[str]
    coordinates: ArrayLike


def _bonds(pairs: list) -> set[tuple[int, int]]:
    return {(min(pair), max(pair)) for pair in pairs}


def test__connectivity(records: list[dict]) -> None:
    """The old bonds of all 900 records: for a TS, of the reactant and the product."""
    for record in records:
        old = record["old"]
        if record["kind"] == "ts":
            sides = [
                ("xyz_reactant", old["bonds"] + old["breaking"]),
                ("xyz_product", old["bonds"] + old["forming"]),
            ]
        else:
            sides = [("xyz", old["bonds"])]
        for key, old_bonds in sides:
            geometry = Geometry(record["symbols"], record[key])
            connectivity = amchi.connectivity(geometry)
            assert set(connectivity.bonds) == _bonds(old_bonds), (record["id"], key)
            assert connectivity.breaking == connectivity.forming == ()


def test__connectivity__hypervalent() -> None:
    """A carbon with five hydrogens in range loses the longest bond."""
    coordinates = [
        [0, 0, 0],
        [1.09, 0, 0],
        [-1.09, 0, 0],
        [0, 1.09, 0],
        [0, -1.09, 0],
        [0, 0, 1.2],
    ]
    geometry = Geometry(["C", *"HHHHH"], coordinates)
    assert amchi.connectivity(geometry).bonds == ((0, 1), (0, 2), (0, 3), (0, 4))


def test__connectivity__errors() -> None:
    """Coordinates of the wrong shape, and an unsupported element."""
    with pytest.raises(ValueError, match="coordinates"):
        amchi.connectivity(Geometry(["H", "H"], [[0, 0, 0]]))
    with pytest.raises(amchi.UnsupportedElementError):
        amchi.connectivity(Geometry(["Fe", "H"], [[0, 0, 0], [1.6, 0, 0]]))
