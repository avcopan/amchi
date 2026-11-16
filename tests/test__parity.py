"""Tests for the local parities."""

from dataclasses import dataclass

import numpy as np
import pytest
from numpy.typing import ArrayLike

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
    atoms = {a: s.parity for a, s in graph.nodes(data="stereo") if s is not None}
    bonds = {(a, b): s.parity for a, b, s in graph.edges(data="stereo") if s}
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


# C with four H neighbors at the corners of a tetrahedron
METHANE = [[0, 0, 0], [1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]
# trans-HC(Cl)=C(F)H: C0=C1 along z, Cl2 and H3 on C0, F4 and H5 on C1
ETHENE = [[0, 0, 0], [0, 0, 1.3], [1, 0, -0.7], [-1, 0, -0.7], [-1, 0, 2], [1, 0, 2]]
ETHENE_SYMBOLS = ["C", "C", "Cl", "H", "F", "H"]
ETHENE_BOND = Candidates(bonds={(0, 1): ((3, 2), (5, 4))})


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
