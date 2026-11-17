"""Tests for the rigid-planar bonds."""

import time
from dataclasses import dataclass

import pytest
from numpy.typing import ArrayLike

import amchi
from amchi import Connectivity
from amchi._graph import Atom, MolGraph, explicit_graph, implicit_graph
from amchi._resonance import pi_systems, rigid_planar_bonds, unsaturations

MAX_SECONDS = 0.05


@dataclass
class Geometry:
    """A minimal GeometryLike."""

    symbols: list[str]
    coordinates: ArrayLike


def _bonds(pairs: list) -> frozenset[tuple[int, int]]:
    return frozenset((min(pair), max(pair)) for pair in pairs)


def _heavy_graph(
    symbols: str, bonds: list[tuple[int, int]], hydrogens: list[int]
) -> MolGraph:
    """A graph of heavy atoms with implicit hydrogens."""
    graph = MolGraph()
    for atom, (symbol, count) in enumerate(zip(symbols, hydrogens, strict=True)):
        graph.add_node(atom, **dict(Atom(symbol=symbol, hydrogens=count)))
    graph.add_edges_from(bonds)
    return graph


def _ring(start: int, size: int = 6) -> list[tuple[int, int]]:
    atoms = list(range(start, start + size))
    return [(atoms[k], atoms[(k + 1) % size]) for k in range(size)]


def _coronene() -> MolGraph:
    """Coronene: inner ring 0-5, spokes 6-11, and outer atoms 12-23 with one H each."""
    bonds = _ring(0)
    for k in range(6):
        spoke, next_spoke = 6 + k, 6 + (k + 1) % 6
        outer1, outer2 = 12 + 2 * k, 13 + 2 * k
        bonds += [(k, spoke), (spoke, outer1), (outer1, outer2), (outer2, next_spoke)]
    return _heavy_graph("C" * 24, bonds, [0] * 12 + [1] * 12)


def test__rigid_planar_bonds(records: list[dict]) -> None:
    """The old rigid-planar bonds of all 900 records: for a TS, of both sides."""
    for record in records:
        old = record["old"]
        if record["kind"] == "ts":
            sides = [
                (old["bonds"] + old["breaking"], old["rigid_planar"]["reactant"]),
                (old["bonds"] + old["forming"], old["rigid_planar"]["product"]),
            ]
        else:
            sides = [(old["bonds"], old["rigid_planar"])]
        for bonds, old_rigid_planar in sides:
            graph = explicit_graph(record["symbols"], Connectivity(bonds=bonds))
            assert rigid_planar_bonds(graph) == _bonds(old_rigid_planar), record["id"]
            # The implicit graph gives the same bonds
            implicit = implicit_graph(graph)
            assert rigid_planar_bonds(implicit) == _bonds(old_rigid_planar), record[
                "id"
            ]


@pytest.mark.parametrize(
    ("name", "symbols", "bonds", "hydrogens", "expected"),
    [
        ("allyl", "CCC", [(0, 1), (1, 2)], [2, 1, 2], {(0, 1), (1, 2)}),
        ("butadiene", "CCCC", [(0, 1), (1, 2), (2, 3)], [2, 1, 1, 2], {(0, 1), (2, 3)}),
        ("propargyl", "CCC", [(0, 1), (1, 2)], [1, 0, 2], set()),
        ("allene", "CCC", [(0, 1), (1, 2)], [2, 0, 2], set()),
        ("acetonitrile", "CCN", [(0, 1), (1, 2)], [3, 0, 0], set()),
        ("vinyl", "CC", [(0, 1)], [2, 1], {(0, 1)}),
        ("benzene", "C" * 6, _ring(0), [1] * 6, set(_bonds(_ring(0)))),
        (
            "biphenyl",
            "C" * 12,
            [*_ring(0), *_ring(6), (0, 6)],
            [0, 1, 1, 1, 1, 1] * 2,
            set(_bonds(_ring(0) + _ring(6))),
        ),
    ],
)
def test__rigid_planar_bonds__by_hand(
    name: str,
    symbols: str,
    bonds: list[tuple[int, int]],
    hydrogens: list[int],
    expected: set[tuple[int, int]],
) -> None:
    """Small molecules, as graphs with implicit hydrogens."""
    graph = _heavy_graph(symbols, bonds, hydrogens)
    assert rigid_planar_bonds(graph) == expected, name


def test__rigid_planar_bonds__subset() -> None:
    """Only the bonds asked about, in either atom order of the graph."""
    graph = _heavy_graph("CCCC", [(0, 1), (1, 2), (2, 3)], [2, 1, 1, 2])
    assert rigid_planar_bonds(graph, bonds=[(1, 2), (2, 3)]) == {(2, 3)}
    assert rigid_planar_bonds(graph, bonds=[]) == frozenset()


def test__unsaturations_and_pi_systems() -> None:
    """A hypervalent atom takes no part in a π system (D23)."""
    # CH2=CH-CH2-CH=O and a five-bonded carbon next to an unsaturated one
    graph = _heavy_graph("CCCCO", [(0, 1), (1, 2), (2, 3), (3, 4)], [2, 1, 2, 1, 0])
    assert unsaturations(graph) == {0: 1, 1: 1, 2: 0, 3: 1, 4: 1}
    assert pi_systems(graph) == [frozenset({0, 1}), frozenset({3, 4})]
    hypervalent = _heavy_graph("CC", [(0, 1)], [5, 2])
    assert unsaturations(hypervalent) == {0: 0, 1: 1}
    assert pi_systems(hypervalent) == []
    assert rigid_planar_bonds(hypervalent) == frozenset()


def test__rigid_planar_bonds__coronene() -> None:
    """Fully aromatic: every bond, which the old code never finished."""
    graph = _coronene()
    assert rigid_planar_bonds(graph) == _bonds(graph.edges)


def _best_time(graph: MolGraph, repeat: int = 3) -> float:
    times = []
    for _ in range(repeat):
        start = time.perf_counter()
        rigid_planar_bonds(graph)
        times.append(time.perf_counter() - start)
    return min(times)


def test__rigid_planar_bonds__speed(hand_cases: dict[str, list[dict]]) -> None:
    """The 17 scaling molecules and coronene, each under 50 ms for all bonds."""
    graphs = {"coronene": _coronene()}
    for case in hand_cases["scaling"]:
        geometry = Geometry(case["symbols"], case["xyz"])
        graphs[case["id"]] = explicit_graph(
            case["symbols"], amchi.connectivity(geometry)
        )
    assert len(graphs) == 18
    for name, graph in graphs.items():
        assert _best_time(graph) < MAX_SECONDS, name
