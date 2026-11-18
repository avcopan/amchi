"""Tests for the stereocenter candidates and their edits."""

from collections.abc import Sequence
from dataclasses import dataclass

import networkx as nx
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


# The valences of the elements of the sample, for the D41 and D42 checks
VALENCE = {"H": 1, "C": 4, "N": 3, "O": 2}


def _graph(symbols: list[str], bonds: Sequence[Sequence[int]]) -> nx.Graph:
    graph = nx.Graph(list(map(tuple, bonds)))
    graph.add_nodes_from(range(len(symbols)))
    return graph


def _planar_nitrogen(
    symbols: list[str], graph: nx.Graph, atom: int, fixed: nx.Graph | None = None
) -> bool:
    """A nitrogen with 3 neighbors, one of them unsaturated, not in a 3-ring (D41).

    For a TS side, the ring must be in `fixed`, the graph without reacting bonds.
    """
    neighbors = list(graph[atom])
    return (
        symbols[atom] == "N"
        and len(neighbors) == 3
        and any(VALENCE[symbols[n]] > graph.degree(n) for n in neighbors)
        and not nx.triangles(graph if fixed is None else fixed, atom)
    )


def _ring_radical(symbols: list[str], graph: nx.Graph, atom: int) -> bool:
    """A radical with 3 saturated neighbors in a 3-membered ring (D42)."""
    neighbors = list(graph[atom])
    return (
        len(neighbors) == 3
        and VALENCE[symbols[atom]] == 4
        and all(VALENCE[symbols[n]] == graph.degree(n) for n in neighbors)
        and bool(nx.triangles(graph, atom))
    )


def test__candidates(records: list[dict]) -> None:
    """The old candidates of the 600 species, with their neighbors, after D41 and D42.

    The old candidates lose the nitrogens with an unsaturated neighbor (D41), and gain
    the radicals in 3-membered rings (D42). The bond candidates don't change.
    """
    dropped, added = 0, 0
    for record in records:
        if record["kind"] == "ts":
            continue
        symbols, old = record["symbols"], record["old"]
        result = amchi.candidates(symbols, Connectivity(bonds=old["bonds"]))
        graph = _graph(symbols, old["bonds"])
        old_atoms = {
            atom: _split(symbols, neighbors)
            for atom, neighbors in old["candidates"]["atoms"]
            if not _planar_nitrogen(symbols, graph, atom)
        }
        dropped += len(old["candidates"]["atoms"]) - len(old_atoms)
        for atom in set(result.atoms) - set(old_atoms):
            assert _ring_radical(symbols, graph, atom), record["id"]
            old_atoms[atom] = _split(symbols, result.atoms[atom])
            added += 1
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
    assert (dropped, added) == (226, 3)


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
    """The documented failure, and its override with the edits.

    The central N of HN=N(O)OH is conjugated, so it is no candidate (D41); the N=N bond
    is still missed.
    """
    for case in hand_cases["hnno"]:
        connectivity = amchi.connectivity(Geometry(case["symbols"], case["xyz"]))
        result = amchi.candidates(case["symbols"], connectivity)
        assert result.atoms == {}
        assert result.bonds == {}
        override = case["override"]
        assert override["remove_atoms"] == [1]
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


def _planar_ts_nitrogens(record: dict) -> set[int]:
    """The nitrogens planar on either side of a TS (D41)."""
    symbols, old = record["symbols"], record["old"]
    fixed = _graph(symbols, old["bonds"])
    return {
        atom
        for side in (old["bonds"] + old["breaking"], old["bonds"] + old["forming"])
        for graph in [_graph(symbols, side)]
        for atom in graph
        if _planar_nitrogen(symbols, graph, atom, fixed)
    }


def test__candidates__ts(records: list[dict]) -> None:
    """The old candidates of the 300 TSs, in both directions (W4, step 7), after D41.

    Read in the given direction, `atoms` are the dumped ones; read reversed, `atoms`
    with `reverse_atoms` in their place are. Swapping breaking and forming swaps the
    two dumps. A nitrogen planar on either side is dropped from both (D41), and a
    rigid planar bond that joined two atom candidates can become a bond candidate
    when one of them is dropped. No TS of the sample has a radical in a 3-ring (D42).
    """
    dropped, added = set(), set()
    for record in records:
        if record["kind"] != "ts":
            continue
        old = record["old"]
        planar = _planar_ts_nitrogens(record)
        given, reverse = old["candidates"], old["reverse_candidates"]
        for breaking, forming, first, second in (
            (old["breaking"], old["forming"], given, reverse),
            (old["forming"], old["breaking"], reverse, given),
        ):
            connectivity = Connectivity(
                bonds=old["bonds"], breaking=breaking, forming=forming
            )
            result = amchi.candidates(record["symbols"], connectivity)
            first_atoms = {a: n for a, n in first["atoms"] if a not in planar}
            second_atoms = {a: n for a, n in second["atoms"] if a not in planar}
            dropped |= {(record["id"], a) for a, _ in first["atoms"] if a in planar}
            old_bonds = {tuple(bond): sides for bond, sides in first["bonds"]}
            new_bonds = {
                bond: sides
                for bond, sides in result.bonds.items()
                if bond not in old_bonds
            }
            for bond in new_bonds:
                assert set(bond) & planar, record["id"]
                added.add((record["id"], bond))
            assert _as_sets(result.atoms, result.bonds) == _as_sets(
                first_atoms, old_bonds | new_bonds
            ), record["id"]
            assert (
                _as_sets(result.atoms | result.reverse_atoms, {})[0]
                == _as_sets(second_atoms, {})[0]
            ), record["id"]
            # Only atoms with different neighbors on the two sides
            for atom, neighbors in result.reverse_atoms.items():
                assert set(neighbors) != set(result.atoms[atom]), record["id"]
    assert (len(dropped), len(added)) == (242, 6)


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


def _molecule(
    heavy: list[str], hydrogens: list[int], bonds: list[tuple[int, int]]
) -> tuple[list[str], Connectivity]:
    """The heavy atoms, then their hydrogens, with all bonds."""
    symbols, all_bonds = list(heavy), list(bonds)
    for atom, count in enumerate(hydrogens):
        for _ in range(count):
            all_bonds.append((atom, len(symbols)))
            symbols.append("H")
    return symbols, Connectivity(bonds=all_bonds)


BENZENE_RING = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5)]


@pytest.mark.parametrize(
    ("heavy", "hydrogens", "bonds", "atom", "is_candidate"),
    [
        # D41: a nitrogen with an unsaturated neighbor is planar
        (["C", "O", "N", "C"], [1, 0, 1, 3], [(0, 1), (0, 2), (2, 3)], 2, False),
        (
            ["C"] * 6 + ["N", "C"],
            [0, 1, 1, 1, 1, 1, 1, 3],
            [*BENZENE_RING, (0, 6), (6, 7)],
            6,
            False,
        ),  # N-methylaniline
        (["C", "C", "N", "C"], [2, 1, 1, 3], [(0, 1), (1, 2), (2, 3)], 2, False),
        (["C", "N", "C"], [2, 1, 3], [(0, 1), (1, 2)], 1, False),  # •CH2NHCH3
        (["C", "N", "C"], [1, 1, 3], [(0, 1), (1, 2)], 1, False),  # HC:NHCH3
        (["C", "N", "O", "O"], [3, 0, 0, 0], [(0, 1), (1, 2), (1, 3)], 1, False),
        # A saturated amine stays
        (["C", "N", "C", "C"], [3, 1, 2, 3], [(0, 1), (1, 2), (2, 3)], 1, True),
        # D42: a 3-ring atom is pyramidal, conjugated or not (N-formylaziridine)
        (
            ["N", "C", "C", "C", "O"],
            [0, 2, 2, 1, 0],
            [(0, 1), (1, 2), (0, 2), (0, 3), (3, 4)],
            0,
            True,
        ),
        # D42: 1-fluoro-2-methylcyclopropyl radical; not a 4-ring or allylic one
        (
            ["C", "C", "C", "C", "F"],
            [0, 1, 2, 3, 0],
            [(0, 1), (1, 2), (0, 2), (0, 4), (1, 3)],
            0,
            True,
        ),
        (
            ["C", "C", "C", "C", "C", "F"],
            [0, 1, 2, 2, 3, 0],
            [(0, 1), (1, 2), (2, 3), (0, 3), (0, 5), (1, 4)],
            0,
            False,
        ),
        (
            ["C", "C", "C", "C", "C"],
            [0, 2, 2, 1, 2],
            [(0, 1), (1, 2), (0, 2), (0, 3), (3, 4)],
            0,
            False,
        ),  # 1-vinylcyclopropyl
    ],
    ids=[
        "amide",
        "aniline",
        "enamine",
        "alpha-amino radical",
        "aminocarbene",
        "nitro",
        "saturated amine",
        "N-acyl aziridine",
        "cyclopropyl radical",
        "cyclobutyl radical",
        "vinylcyclopropyl radical",
    ],
)
def test__candidates__low_temperature(
    heavy: list[str],
    hydrogens: list[int],
    bonds: list[tuple[int, int]],
    atom: int,
    is_candidate: bool,  # noqa: FBT001
) -> None:
    """Conjugated nitrogens are planar (D41), atoms in 3-membered rings not (D42)."""
    symbols, connectivity = _molecule(heavy, hydrogens, bonds)
    assert (atom in amchi.candidates(symbols, connectivity).atoms) == is_candidate


def test__candidates__low_temperature_ts() -> None:
    """H abstraction from N-methylethylamine (D41 in a TS).

    Next to the N, the product is an alpha-amino radical, so the N is planar there and
    is a candidate on neither side. At the far carbon, the N stays a candidate.
    """
    symbols, connectivity = _molecule(
        ["C", "N", "C", "C", "H"], [3, 1, 2, 3, 0], [(0, 1), (1, 2), (2, 3)]
    )
    for carbon, expected in ((2, False), (3, True)):
        hydrogen = next(j for i, j in connectivity.bonds if i == carbon)
        ts = Connectivity(
            bonds=[b for b in connectivity.bonds if b != (carbon, hydrogen)],
            breaking=[(carbon, hydrogen)],
            forming=[(4, hydrogen)],
        )
        result = amchi.candidates(symbols, ts)
        assert (1 in result.atoms) == expected, carbon
        assert 1 not in result.reverse_atoms


def test__candidates__ts_forming_3_ring() -> None:
    """A 3-ring closed by a forming bond doesn't keep a conjugated N (D41, D42).

    The diradical •N(CHO)CH2CH2• closes to N-formylaziridine, whose N is a candidate.
    In the TS, the N has 3 neighbors on the product side only, and is conjugated there;
    its 3-ring has the forming bond, so it doesn't count.
    """
    heavy, hydrogens = ["N", "C", "C", "C", "O"], [0, 2, 2, 1, 0]
    bonds = [(0, 1), (1, 2), (0, 3), (3, 4)]
    symbols, product = _molecule(heavy, hydrogens, [*bonds, (0, 2)])
    assert 0 in amchi.candidates(symbols, product).atoms
    _, reactant = _molecule(heavy, hydrogens, bonds)
    ts = Connectivity(bonds=reactant.bonds, forming=[(0, 2)])
    assert 0 not in amchi.candidates(symbols, ts).atoms
