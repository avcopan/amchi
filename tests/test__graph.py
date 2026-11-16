"""Tests for the molecular graph."""

import pytest
from pydantic import ValidationError

from amchi import Connectivity
from amchi._graph import (
    AtomStereo,
    MolGraph,
    connected_components,
    explicit_graph,
    implicit_graph,
)


def _graph(record: dict) -> MolGraph:
    old = record["old"]
    connectivity = Connectivity(
        bonds=old["bonds"], breaking=old["breaking"], forming=old["forming"]
    )
    return implicit_graph(explicit_graph(record["symbols"], connectivity))


def test__implicit_graph(records: list[dict]) -> None:
    """The nodes and components are the old ones, for all 900 records."""
    for record in records:
        graph = _graph(record)
        graph.validate()
        canonical = record["old"]["canonical"]
        old_nodes = {atom for form in canonical for atom, _ in form["priorities"]}
        assert set(graph) == old_nodes, record["id"]
        assert len(connected_components(graph)) == len(canonical), record["id"]
        if record["kind"] != "ts":
            hydrogens = sum(h for _, h in graph.nodes(data="hydrogens"))
            assert hydrogens == record["symbols"].count("H"), record["id"]


@pytest.mark.parametrize(
    ("symbols", "lists", "nodes", "hydrogens"),
    [
        (["H", "H"], {"bonds": [(0, 1)]}, {0}, {0: 1}),  # H2
        (["H"], {"bonds": []}, {0}, {0: 0}),  # H atom
        (["O", "H", "H"], {"bonds": [(0, 1), (0, 2)]}, {0}, {0: 2}),  # water
        # H transfer between two O atoms: the moving H is a backbone hydrogen
        (
            ["O", "H", "O", "H"],
            {"bonds": [(0, 3)], "breaking": [(0, 1)], "forming": [(1, 2)]},
            {0, 1, 2},
            {0: 1, 1: 0, 2: 0},
        ),
        # H abstraction from H2: the H of the breaking bond stays, the other goes
        (
            ["H", "H", "O"],
            {"bonds": [], "breaking": [(0, 1)]},
            {0, 1, 2},
            {0: 0, 1: 0, 2: 0},
        ),
    ],
)
def test__implicit_graph__hydrogens(
    symbols: list[str], lists: dict, nodes: set[int], hydrogens: dict[int, int]
) -> None:
    """Backbone hydrogens stay nodes, and the rest are counted on their neighbor."""
    graph = implicit_graph(explicit_graph(symbols, Connectivity(**lists)))
    assert set(graph) == nodes
    assert dict(graph.nodes(data="hydrogens")) == hydrogens


def test__explicit_graph() -> None:
    """Each bond has its change, and the fields validate."""
    connectivity = Connectivity(bonds=[(0, 1)], breaking=[(1, 2)], forming=[(2, 3)])
    graph = explicit_graph(["C", "C", "H", "O"], connectivity)
    assert {(i, j): c for i, j, c in graph.edges(data="change")} == {
        (0, 1): None,
        (1, 2): "breaking",
        (2, 3): "forming",
    }
    assert connected_components(graph) == [frozenset({0, 1, 2, 3})]
    graph.nodes[0]["stereo"] = AtomStereo(neighbors=(1, 2, 3), parity=None)
    graph.validate()
    graph.edges[0, 1]["change"] = "moving"
    with pytest.raises(ValidationError):
        graph.validate()
