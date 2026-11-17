"""Tests for the local parities."""

from dataclasses import dataclass

import numpy as np
import pytest
from numpy.typing import ArrayLike

import amchi
from amchi import Candidates, Connectivity, parity_margins
from amchi._graph import MolGraph, explicit_graph, implicit_graph
from amchi._parity import stereo_graph


@dataclass
class Geometry:
    """A minimal GeometryLike."""

    symbols: list[str]
    coordinates: ArrayLike


def _old_candidates(old_candidates: dict) -> Candidates:
    return Candidates(
        atoms=dict(old_candidates["atoms"]),
        bonds={tuple(bond): tuple(sides) for bond, sides in old_candidates["bonds"]},
    )


def _parities(graph: MolGraph) -> tuple[dict, dict]:
    """The atom parities, and each bond's table entry of its last two neighbors."""
    atoms = {a: s.parity for a, s in graph.nodes(data="stereo") if s is not None}
    bonds = {(a, b): s.parities[-1][-1] for a, b, s in graph.edges(data="stereo") if s}
    return atoms, bonds


def _stereo_graph(symbols: list[str], xyz: list, candidates: Candidates) -> MolGraph:
    bonds = [*candidates.bonds]
    bonds += [(a, n) for a, ns in candidates.atoms.items() for n in ns]
    bonds += [
        (a, n)
        for (a, b), sides in candidates.bonds.items()
        for a, ns in zip((a, b), sides, strict=True)
        for n in ns
    ]
    graph = explicit_graph(symbols, Connectivity(bonds=set(map(frozenset, bonds))))
    return stereo_graph(graph, Geometry(symbols, xyz), candidates)


def test__stereo_graph(records: list[dict]) -> None:
    """With the old neighbor order, the parities are the old ones, for all 900."""
    for record in records:
        old = record["old"]
        connectivity = Connectivity(
            bonds=old["bonds"], breaking=old["breaking"], forming=old["forming"]
        )
        graph = implicit_graph(explicit_graph(record["symbols"], connectivity))
        geometry = Geometry(record["symbols"], record["xyz"])
        graph = stereo_graph(graph, geometry, _old_candidates(old["candidates"]))
        graph.validate()
        atoms, bonds = _parities(graph)
        old_atoms = dict(old["local_parities"]["atoms"])
        old_bonds = {tuple(b): p for b, p in old["local_parities"]["bonds"]}
        assert atoms == old_atoms, record["id"]
        assert bonds == old_bonds, record["id"]
        assert None not in atoms.values()
        assert None not in bonds.values()


def test__stereo_graph__record(records: list[dict]) -> None:
    """The manual check of the walkthrough."""
    record = next(r for r in records if r["id"] == "stjohn:443143_dfe711")
    candidates = Candidates(
        atoms={2: (10, 11, 9, 1), 3: (12, 1, 4, 5)}, bonds={(0, 1): ((8, 7), (2, 3))}
    )
    graph = _stereo_graph(record["symbols"], record["xyz"], candidates)
    assert _parities(graph) == ({2: True, 3: True}, {(0, 1): False})


def test__stereo_graph__ts(records: list[dict]) -> None:
    """The 300 TSs: each neighbor list of an atom gets its own parity (D26).

    With the derived candidates, the reverse parity of an atom is the parity measured
    with the reverse neighbors as its only list, and every one is well-defined.
    """
    count = 0
    for record in records:
        if record["kind"] != "ts":
            continue
        old, symbols = record["old"], record["symbols"]
        connectivity = Connectivity(
            bonds=old["bonds"], breaking=old["breaking"], forming=old["forming"]
        )
        candidates = amchi.candidates(symbols, connectivity)
        geometry = Geometry(symbols, record["xyz"])
        graph = stereo_graph(
            explicit_graph(symbols, connectivity), geometry, candidates
        )
        graph.validate()
        reverse = stereo_graph(
            explicit_graph(symbols, connectivity),
            geometry,
            Candidates(atoms=candidates.reverse_atoms),
        )
        for atom, stereo in graph.nodes(data="stereo"):
            if atom in candidates.reverse_atoms:
                count += 1
                expected = reverse.nodes[atom]["stereo"]
                assert stereo.reverse_neighbors == expected.neighbors, record["id"]
                assert stereo.reverse_parity == expected.parity, record["id"]
                assert stereo.reverse_parity is not None, record["id"]
            elif stereo is not None:
                assert stereo.reverse_neighbors is None, record["id"]
                assert stereo.reverse_parity is None, record["id"]
    assert count == 253  # the substitution sites of the sample (W4)


# C with four H neighbors at the corners of a tetrahedron
METHANE = [[0, 0, 0], [1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]
# trans-HC(Cl)=C(F)H: C0=C1 along z, Cl2 and H3 on C0, F4 and H5 on C1
ETHENE = [[0, 0, 0], [0, 0, 1.3], [1, 0, -0.7], [-1, 0, -0.7], [-1, 0, 2], [1, 0, 2]]
ETHENE_SYMBOLS = ["C", "C", "Cl", "H", "F", "H"]
ETHENE_BOND = Candidates(bonds={(0, 1): ((3, 2), (5, 4))})


def test__stereo_graph__table() -> None:
    """A bond has a parity for each pair of neighbors across it: trans is True."""
    graph = _stereo_graph(ETHENE_SYMBOLS, ETHENE, ETHENE_BOND)
    stereo = graph.edges[0, 1]["stereo"]
    assert stereo.neighbors == ((3, 2), (5, 4))
    # H3-H5 and Cl2-F4 are trans, H3-F4 and Cl2-H5 cis
    assert stereo.parities == ((True, False), (False, True))


def test__stereo_graph__reflection() -> None:
    """Reflecting the coordinates inverts an atom parity and keeps a bond parity."""
    mirror = (np.array(METHANE) * [-1, 1, 1]).tolist()
    candidates = Candidates(atoms={0: (1, 2, 3, 4)})
    atoms, _ = _parities(_stereo_graph(["C", *"HHHH"], METHANE, candidates))
    mirror_atoms, _ = _parities(_stereo_graph(["C", *"HHHH"], mirror, candidates))
    assert atoms == {0: False}
    assert mirror_atoms == {0: True}

    mirror = (np.array(ETHENE) * [-1, 1, 1]).tolist()
    _, bonds = _parities(_stereo_graph(ETHENE_SYMBOLS, ETHENE, ETHENE_BOND))
    _, mirror_bonds = _parities(_stereo_graph(ETHENE_SYMBOLS, mirror, ETHENE_BOND))
    assert bonds == mirror_bonds == {(0, 1): True}


def test__stereo_graph__reverse() -> None:
    """Exchanging one neighbor for another in the same place keeps the parity.

    Methane's H4 replaced by a point on the same corner, and by one on the opposite
    side of the plane of H1, H2 and H3, which inverts the parity.
    """
    coordinates = [*METHANE, [-2, -2, 2], [0.2, 0.2, -1.6]]
    symbols = ["C", *"HHHHHH"]
    candidates = Candidates(atoms={0: (1, 2, 3, 4)}, reverse_atoms={0: (1, 2, 3, 5)})
    stereo = _stereo_graph(symbols, coordinates, candidates).nodes[0]["stereo"]
    assert (stereo.parity, stereo.reverse_parity) == (False, False)
    candidates = Candidates(atoms={0: (1, 2, 3, 4)}, reverse_atoms={0: (1, 2, 3, 6)})
    stereo = _stereo_graph(symbols, coordinates, candidates).nodes[0]["stereo"]
    assert (stereo.reverse_neighbors, stereo.reverse_parity) == ((1, 2, 3, 6), True)


@pytest.mark.parametrize(
    ("neighbors", "parity"),
    [((1, 2, 3, 4), False), ((2, 1, 3, 4), True), ((1, 2, 4, 3), True)],
)
def test__stereo_graph__swap(neighbors: tuple[int, ...], parity: bool) -> None:  # noqa: FBT001
    """Swapping two neighbors of an atom inverts its parity."""
    candidates = Candidates(atoms={0: neighbors})
    atoms, _ = _parities(_stereo_graph(["C", *"HHHH"], METHANE, candidates))
    assert atoms == {0: parity}


def test__stereo_graph__ill_defined() -> None:
    """Four points in a plane, or reference neighbors at 90°, give no parity."""
    planar = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [-1, -1, 0]]
    candidates = Candidates(atoms={0: (1, 2, 3, 4)})
    atoms, _ = _parities(_stereo_graph(["C", *"HHHH"], planar, candidates))
    assert atoms == {0: None}

    twisted = [
        [0, 0, 0],
        [0, 0, 1.3],
        [1, 0, -0.7],
        [-1, 0, -0.7],
        [0, 1, 2],
        [0, -1, 2],
    ]
    _, bonds = _parities(_stereo_graph(ETHENE_SYMBOLS, twisted, ETHENE_BOND))
    assert bonds == {(0, 1): None}


def test__parity_margins() -> None:
    """A tetrahedral atom, a planar atom, a planar bond and a twisted bond."""
    margins = parity_margins(
        Geometry(["C", *"HHHH"], METHANE), Candidates(atoms={0: (1, 2, 3, 4)})
    )
    assert margins.atoms[0] == pytest.approx(3.08, abs=0.005)

    planar = [[0, 0, 0], [1, 0, 0], [-0.5, 0.8, 0], [-0.5, -0.8, 0]]
    margins = parity_margins(
        Geometry(["N", *"HHH"], planar), Candidates(atoms={0: (1, 2, 3)})
    )
    assert margins.atoms == {0: 0.0}

    margins = parity_margins(Geometry(ETHENE_SYMBOLS, ETHENE), ETHENE_BOND)
    assert margins.bonds[0, 1] == pytest.approx(1.0)

    twisted = [
        [0, 0, 0],
        [0, 0, 1.3],
        [1, 0, -0.7],
        [-1, 0, -0.7],
        [0, 1, 2],
        [0, -1, 2],
    ]
    margins = parity_margins(Geometry(ETHENE_SYMBOLS, twisted), ETHENE_BOND)
    assert margins.bonds == {(0, 1): 0.0}

    # A reference neighbor on the bond axis
    linear = [[0, 0, 0], [0, 0, 1.3], [0, 0, -1], [-1, 0, -0.7], [-1, 0, 2], [1, 0, 2]]
    margins = parity_margins(Geometry(ETHENE_SYMBOLS, linear), ETHENE_BOND)
    assert margins.bonds == {(0, 1): 0.0}


def test__parity_margins__ts() -> None:
    """A bond's margin is the smallest over its table; reverse atoms get their own."""
    # trans-HC(Cl)=C(F)H, with a third neighbor on C1 at 90° to the plane
    coordinates = [*ETHENE, [0, 1, 1.9]]
    symbols = [*ETHENE_SYMBOLS, "H"]
    candidates = Candidates(bonds={(0, 1): ((3, 2), (5, 4, 6))})
    margins = parity_margins(Geometry(symbols, coordinates), candidates)
    assert margins.bonds == {(0, 1): pytest.approx(0.0)}
    candidates = Candidates(bonds={(0, 1): ((3, 2), (5, 4))})
    assert parity_margins(Geometry(symbols, coordinates), candidates).bonds == {
        (0, 1): pytest.approx(1.0)
    }

    # Reverse neighbors H1, H2 and H5 in one plane with the C: margin 0
    coordinates = [*METHANE, [-1, 0, 0]]
    candidates = Candidates(atoms={0: (1, 2, 3)}, reverse_atoms={0: (1, 2, 5)})
    margins = parity_margins(Geometry(["C", *"HHHHH"], coordinates), candidates)
    assert margins.atoms[0] > 0.5
    assert margins.reverse_atoms == {0: pytest.approx(0.0)}
