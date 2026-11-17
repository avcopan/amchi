"""Tests for the layer strings.

The formula and h layers are tested against the old strings, and the c layer, which is
an intended difference (D24), against InChI's c layers (`tests/data/inchi.json`).
"""

import json
import re
from dataclasses import dataclass
from pathlib import Path

import networkx as nx
import pytest
from numpy.typing import ArrayLike

import amchi
from amchi import Candidates, Connectivity
from amchi._canon import canonical_forms
from amchi._graph import MolGraph, explicit_graph, implicit_graph
from amchi._layers import c_layer, canonical_graph, formula_layer, h_layer
from amchi._parity import stereo_graph


@dataclass
class Geometry:
    """A minimal GeometryLike."""

    symbols: list[str]
    coordinates: ArrayLike


def _graph(
    symbols: list[str], xyz: ArrayLike, connectivity: Connectivity, cands: Candidates
) -> MolGraph:
    graph = implicit_graph(explicit_graph(symbols, connectivity))
    return stereo_graph(graph, Geometry(symbols, xyz), cands)


def _old_candidates(old_candidates: dict) -> Candidates:
    return Candidates(
        atoms=dict(old_candidates["atoms"]),
        bonds={tuple(bond): tuple(sides) for bond, sides in old_candidates["bonds"]},
    )


def _layers(graph: MolGraph) -> dict[str, str]:
    """The formula and h layers of a one-component graph."""
    (form,) = canonical_forms(graph)
    canonical = canonical_graph(graph, form)
    return {"formula": formula_layer(canonical), "h": h_layer(canonical)}


def _old_layers(string: str) -> dict[str, str]:
    """The formula and h layers of an old one-component string ("" if absent)."""
    _, formula, *layers = string.split("/")
    prefixed = {layer[0]: layer[1:] for layer in layers}
    return {"formula": formula, "h": prefixed.get("h", "")}


def _inchi_graph(layer: str) -> nx.Graph:
    """The graph of an InChI c layer, keyed by InChI's numbers minus one."""
    graph = nx.Graph()
    previous = None
    branch_points = []
    for part in re.findall(r"\d+|[-(),]", layer):
        if part == "(":
            branch_points.append(previous)
        elif part == ",":
            previous = branch_points[-1]
        elif part == ")":
            previous = branch_points.pop()
        elif part != "-":
            atom = int(part) - 1
            graph.add_node(atom)
            if previous is not None:
                graph.add_edge(previous, atom)
            previous = atom
    return graph


def test__layers(records: list[dict]) -> None:
    """The formula and h layers of the 600 species equal the old ones."""
    species = [r for r in records if r["kind"] != "ts"]
    assert len(species) == 600
    for record in species:
        old = record["old"]
        connectivity = Connectivity(bonds=old["bonds"])
        candidates = _old_candidates(old["candidates"])
        graph = _graph(record["symbols"], record["xyz"], connectivity, candidates)
        assert _layers(graph) == _old_layers(old["amchi"]), record["id"]


def test__c_layer__inchi() -> None:
    """The c layer of a graph in InChI's numbering is InChI's c layer.

    The InChIs of the 600 species and the hand cases were made with RDKit
    (`baseline/inchi_c_layers.py` in the planning repo).
    """
    path = Path(__file__).parent / "data" / "inchi.json"
    rows = json.loads(path.read_text())
    assert len(rows) == 625
    rings = {0: 0, 1: 0, 2: 0}
    for row in rows:
        (layer,) = (x[1:] for x in row["inchi"].split("/") if x.startswith("c"))
        graph = _inchi_graph(layer)
        assert sorted(graph) == list(range(len(graph))), row["id"]
        assert c_layer(graph) == layer, row["id"]
        count = graph.number_of_edges() - graph.number_of_nodes() + 1
        rings[min(count, 2)] += 1
    # Chains, single rings and polycycles are all covered
    assert all(rings.values()), rings


def test__layers__hand_cases(hand_cases: dict[str, list[dict]]) -> None:
    """The layers of the hand cases, from the helpers, equal the old ones.

    The scaling group has the large polycycles (pyrene, perylene).
    """
    for group in ("hnno", "ring", "ring_control", "scaling"):
        for case in hand_cases[group]:
            geometry = Geometry(case["symbols"], case["xyz"])
            connectivity = amchi.connectivity(geometry)
            candidates = amchi.candidates(case["symbols"], connectivity)
            graph = _graph(case["symbols"], case["xyz"], connectivity, candidates)
            assert _layers(graph) == _old_layers(case["amchi"]), case["id"]


@pytest.mark.parametrize(
    ("symbols", "bonds", "expected"),
    [
        (
            ["C", "H", "H", "H", "H"],
            [(0, 1), (0, 2), (0, 3), (0, 4)],
            ("CH4", "", "1H4"),
        ),
        (["O", "H", "H"], [(0, 1), (0, 2)], ("H2O", "", "1H2")),
        (["N"], [], ("N", "", "")),
        (["H", "H"], [(0, 1)], ("H2", "", "1H")),
    ],
)
def test__layers__small(
    symbols: list[str], bonds: list[tuple[int, int]], expected: tuple[str, str, str]
) -> None:
    """A lone heavy atom has no c layer; H2 keeps one hydrogen as a node, as `h1H`."""
    xyz = [[float(i), 0.0, 0.0] for i in range(len(symbols))]
    graph = _graph(symbols, xyz, Connectivity(bonds=bonds), Candidates())
    (form,) = canonical_forms(graph)
    canonical = canonical_graph(graph, form)
    layers = (formula_layer(canonical), c_layer(canonical), h_layer(canonical))
    assert layers == expected


@pytest.mark.parametrize(
    ("bonds", "expected"),
    [
        ([(0, 1), (1, 2), (2, 3)], "1-2-3-4"),
        ([(0, 2), (1, 2), (2, 3), (3, 4)], "1-3(2)4-5"),
        ([(0, 1), (1, 3), (3, 5), (5, 4), (4, 2), (2, 0)], "1-2-4-6-5-3-1"),
        ([(0, 3), (1, 3), (2, 3), (3, 4), (4, 5)], "1-4(2,3)5-6"),
    ],
)
def test__c_layer__small(bonds: list[tuple[int, int]], expected: str) -> None:
    """Branches are written smallest first, and a ring closes at its later atom."""
    assert c_layer(nx.Graph(bonds)) == expected
