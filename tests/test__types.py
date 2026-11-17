"""Tests for the input models."""

import pytest
from pydantic import ValidationError

from amchi import AmchiError, Candidates, Connectivity, IllDefinedParityError


def test__connectivity() -> None:
    """Pairs and lists are sorted, and lists of lists are accepted."""
    connectivity = Connectivity(bonds=[[2, 0], [1, 0]], breaking=[[4, 3]])
    assert connectivity.bonds == ((0, 1), (0, 2))
    assert connectivity.breaking == ((3, 4),)
    assert connectivity.forming == ()


@pytest.mark.parametrize(
    "lists",
    [
        {"bonds": [(0, 1)], "breaking": [(1, 0)]},  # a pair in two lists
        {"bonds": [(0, 1), (1, 0)]},  # a repeated pair
        {"bonds": [(1, 1)]},  # a self-bond
        {"bonds": [(-1, 1)]},  # a negative index
        {"bonds": [(0, 1, 2)]},  # not a pair
        {"bonds": [("a", 1)]},  # not an index
        {"bonds": 1},  # not a list
    ],
)
def test__connectivity__invalid(lists: dict) -> None:
    """Invalid bonds are rejected."""
    with pytest.raises(ValidationError):
        Connectivity(**lists)


def test__candidates() -> None:
    """The neighbor order is kept, and a reversed bond key swaps the sides."""
    candidates = Candidates(
        atoms={0: [4, 3, 2, 1], 5: [6, 7, 8]},
        bonds={(2, 1): ([9, 10, 11], [12])},
    )
    assert candidates.atoms == {0: (4, 3, 2, 1), 5: (6, 7, 8)}
    assert candidates.bonds == {(1, 2): ((12,), (9, 10, 11))}
    assert Candidates().atoms == Candidates().bonds == {}
    assert Candidates().reverse_atoms == {}


def test__candidates__reverse_atoms() -> None:
    """A TS atom can have product-side neighbors, in the order given."""
    candidates = Candidates(atoms={0: [1, 2, 3, 4]}, reverse_atoms={0: [5, 2, 3, 4]})
    assert candidates.reverse_atoms == {0: (5, 2, 3, 4)}


@pytest.mark.parametrize(
    "candidates",
    [
        {"atoms": {0: [1, 2]}},  # 2 neighbors
        {"atoms": {0: [1, 2, 3, 4, 5]}},  # 5 neighbors
        {"atoms": {0: [1, 2, 2]}},  # a repeated neighbor
        {"atoms": {0: [0, 1, 2]}},  # the atom itself
        {"bonds": {(0, 1): ([], [2])}},  # an empty side
        {"bonds": {(0, 1): ([2, 3, 4, 5], [6])}},  # 4 on a side
        {"bonds": {(0, 1): ([1], [2])}},  # the other bond atom
        {"bonds": {(0, 0): ([2], [3])}},  # a self-bond
        {"bonds": {(0, 1, 2): ([3], [4])}},  # not a pair
        {"bonds": 1},  # not a dict
        {"atoms": {0: [1, 2, 3]}, "reverse_atoms": {0: [1, 2]}},  # 2 neighbors
        {"atoms": {0: [1, 2, 3]}, "reverse_atoms": {0: [0, 1, 2]}},  # the atom
        {"atoms": {0: [1, 2, 3]}, "reverse_atoms": {5: [1, 2, 3]}},  # not in atoms
    ],
)
def test__candidates__invalid(candidates: dict) -> None:
    """Invalid candidates are rejected."""
    with pytest.raises(ValidationError):
        Candidates(**candidates)


def test__errors() -> None:
    """The parity error is an amchi error."""
    assert issubclass(IllDefinedParityError, AmchiError)
