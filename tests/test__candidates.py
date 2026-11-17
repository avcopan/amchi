"""Tests for the stereocenter candidates and their edits."""

from collections.abc import Sequence
from dataclasses import dataclass

import pytest
from numpy.typing import ArrayLike

import amchi
from amchi import Candidates, Connectivity


@dataclass
class Geometry:
    """A minimal GeometryLike."""

    symbols: list[str]
    coordinates: ArrayLike


def _split(symbols: list[str], neighbors: Sequence[int]) -> tuple[frozenset, tuple]:
    """The hydrogens as a set, and the other neighbors in order."""
    hydrogens = frozenset(n for n in neighbors if symbols[n] == "H")
    return hydrogens, tuple(n for n in neighbors if symbols[n] != "H")


def test__candidates(records: list[dict]) -> None:
    """The old candidates of the 600 species, with their neighbors."""
    for record in records:
        if record["kind"] == "ts":
            continue
        symbols, old = record["symbols"], record["old"]
        result = amchi.candidates(symbols, Connectivity(bonds=old["bonds"]))
        old_atoms = {
            atom: _split(symbols, neighbors)
            for atom, neighbors in old["candidates"]["atoms"]
        }
        old_bonds = {
            tuple(bond): tuple(_split(symbols, side) for side in sides)
            for bond, sides in old["candidates"]["bonds"]
        }
        atoms = {atom: _split(symbols, n) for atom, n in result.atoms.items()}
        bonds = {
            bond: tuple(_split(symbols, side) for side in sides)
            for bond, sides in result.bonds.items()
        }
        assert atoms == old_atoms, record["id"]
        assert bonds == old_bonds, record["id"]
        # Hydrogens first (D21)
        for neighbors in [*result.atoms.values(), *sum(result.bonds.values(), ())]:
            is_h = [symbols[n] == "H" for n in neighbors]
            assert is_h == sorted(is_h, reverse=True), record["id"]


def test__candidates__stjohn_279708(records: list[dict]) -> None:
    """C#C[C@](C)(O)CC([CH2])=O: the bond (6, 7) is the only candidate.

    Of the rigid-planar bonds (6, 7) and (6, 8), (6, 8) fails, since O8 has no other
    neighbor. The triple bond (0, 1) isn't rigid and planar.
    """
    (record,) = [r for r in records if r["id"] == "stjohn:279708_c2065f"]
    result = amchi.candidates(
        record["symbols"], Connectivity(bonds=record["old"]["bonds"])
    )
    assert result.bonds == {(6, 7): ((5, 8), (16, 17))}


def test__candidates__hnno(hand_cases: dict[str, list[dict]]) -> None:
    """The documented failure, and its override with the edits."""
    for case in hand_cases["hnno"]:
        connectivity = amchi.connectivity(Geometry(case["symbols"], case["xyz"]))
        result = amchi.candidates(case["symbols"], connectivity)
        assert list(result.atoms) == [1]
        assert result.bonds == {}
        override = case["override"]
        for atom in override["remove_atoms"]:
            result = amchi.remove_candidate(result, atom)
        for bond in override["add_bonds"]:
            result = amchi.add_bond_candidate(result, connectivity, bond)
        assert result.atoms == {}
        assert list(result.bonds) == [(0, 1)]


def _as_sets(atoms: dict, bonds: dict) -> tuple[dict, dict]:
    """Each neighbor list as a set."""
    return (
        {atom: set(neighbors) for atom, neighbors in atoms.items()},
        {bond: tuple(map(set, sides)) for bond, sides in bonds.items()},
    )


def test__candidates__ts(records: list[dict]) -> None:
    """The old candidates of the 300 TSs, in both directions (W4, step 7).

    Read in the given direction, `atoms` are the dumped ones; read reversed, `atoms`
    with `reverse_atoms` in their place are. Swapping breaking and forming swaps the
    two dumps.
    """
    for record in records:
        if record["kind"] != "ts":
            continue
        old = record["old"]
        given, reverse = old["candidates"], old["reverse_candidates"]
        for breaking, forming, first, second in (
            (old["breaking"], old["forming"], given, reverse),
            (old["forming"], old["breaking"], reverse, given),
        ):
            connectivity = Connectivity(
                bonds=old["bonds"], breaking=breaking, forming=forming
            )
            result = amchi.candidates(record["symbols"], connectivity)
            old_bonds = {tuple(bond): sides for bond, sides in first["bonds"]}
            assert _as_sets(result.atoms, result.bonds) == _as_sets(
                dict(first["atoms"]), old_bonds
            ), record["id"]
            assert (
                _as_sets(result.atoms | result.reverse_atoms, {})[0]
                == _as_sets(dict(second["atoms"]), {})[0]
            ), record["id"]
            # Only atoms with different neighbors on the two sides
            for atom, neighbors in result.reverse_atoms.items():
                assert set(neighbors) != set(result.atoms[atom]), record["id"]


def test__candidates__substitution() -> None:
    """F moves from one C of CH3F to the other: two neighbor lists for the C."""
    connectivity = Connectivity(
        bonds=[(0, 1), (0, 2), (0, 3)], breaking=[(0, 4)], forming=[(0, 5)]
    )
    result = amchi.candidates(["C", "H", "H", "H", "F", "F"], connectivity)
    assert result.atoms == {0: (1, 2, 3, 4)}
    assert result.reverse_atoms == {0: (1, 2, 3, 5)}
    # Tetrahedral on the product side only: the product-side neighbors, no reverse
    connectivity = Connectivity(bonds=[(0, 1), (0, 2), (0, 3)], forming=[(0, 4)])
    result = amchi.candidates(["C", "H", "H", "H", "F"], connectivity)
    assert result.atoms == {0: (1, 2, 3, 4)}
    assert result.reverse_atoms == {}


def test__candidates__ts_ring() -> None:
    """The ring test leaves and reaches a bond through its listed neighbors only.

    A ring closure: C3 of CH2=CH-CH2-CH2 adds to C0, so C0=C1 is rigid and planar on
    the reactant side only, where C0 lists H4 and H5. The 4-ring C0-C1-C2-C3 of the
    TS graph leaves C0 through C3, which isn't listed, so the bond stays a candidate.
    """
    symbols = ["C", "C", "C", "C", "H", "H", "H", "H", "H", "H", "H"]
    bonds = [(0, 1), (1, 2), (2, 3), (0, 4), (0, 5), (1, 6)]
    bonds += [(2, 7), (2, 8), (3, 9), (3, 10)]
    connectivity = Connectivity(bonds=bonds, forming=[(0, 3)])
    result = amchi.candidates(symbols, connectivity)
    assert result.bonds == {(0, 1): ((4, 5), (6, 2))}
    # C0 and C3 are tetrahedral on the product side only
    assert result.atoms == {0: (4, 5, 1, 3), 2: (7, 8, 1, 3), 3: (9, 10, 0, 2)}
    assert result.reverse_atoms == {}


def test__add_bond_candidate__ts() -> None:
    """A TS bond candidate lists the neighbors over all three bond lists."""
    connectivity = Connectivity(
        bonds=[(0, 1), (0, 2), (1, 3), (5, 6)], breaking=[(1, 4)], forming=[(0, 5)]
    )
    result = amchi.add_bond_candidate(Candidates(), connectivity, (0, 1))
    assert result.bonds == {(0, 1): ((2, 5), (3, 4))}
    # The bond itself may be a reacting one
    result = amchi.add_bond_candidate(Candidates(), connectivity, (5, 0))
    assert result.bonds == {(0, 5): ((1, 2), (6,))}
    with pytest.raises(ValueError, match="Not a bond"):
        amchi.add_bond_candidate(Candidates(), connectivity, (2, 3))


def test__candidates__small_rings() -> None:
    """A double bond in a 7-ring is no candidate; in an 8-ring it is."""
    for size, expected in [(7, False), (8, True)]:
        ring = [(k, (k + 1) % size) for k in range(size)]
        # One hydrogen on atoms 0 and 1, two on the others
        carbons = [0, 1] + [k for k in range(2, size) for _ in range(2)]
        hydrogens = [(k, size + n) for n, k in enumerate(carbons)]
        symbols = ["C"] * size + ["H"] * len(carbons)
        result = amchi.candidates(symbols, Connectivity(bonds=ring + hydrogens))
        assert ((0, 1) in result.bonds) is expected


def test__add_bond_candidate__errors() -> None:
    """Not a bond, or a side without neighbors."""
    connectivity = Connectivity(bonds=[(0, 1), (1, 2)])
    with pytest.raises(ValueError, match="Not a bond"):
        amchi.add_bond_candidate(Candidates(), connectivity, (0, 2))
    with pytest.raises(ValueError, match="neighbor on each side"):
        amchi.add_bond_candidate(Candidates(), connectivity, (0, 1))


def test__remove_candidate__errors() -> None:
    """Not a candidate; the input is left unchanged."""
    candidates = Candidates(atoms={0: (1, 2, 3)}, bonds={(4, 5): ((6,), (7,))})
    with pytest.raises(KeyError):
        amchi.remove_candidate(candidates, 1)
    with pytest.raises(KeyError):
        amchi.remove_candidate(candidates, (0, 1))
    # The input is unchanged
    amchi.remove_candidate(candidates, 0)
    assert candidates.atoms == {0: (1, 2, 3)}


def test__remove_candidate__reverse_atoms() -> None:
    """An atom goes with its reverse neighbors; the others keep theirs."""
    candidates = Candidates(
        atoms={0: (1, 2, 3), 4: (5, 6, 7)},
        reverse_atoms={0: (1, 2, 8), 4: (5, 6, 9)},
        bonds={(10, 11): ((12,), (13,))},
    )
    result = amchi.remove_candidate(candidates, 0)
    assert result.reverse_atoms == {4: (5, 6, 9)}
    result = amchi.remove_candidate(candidates, (10, 11))
    assert result.reverse_atoms == candidates.reverse_atoms
