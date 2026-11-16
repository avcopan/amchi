"""The NetworkX graph with validated atom and bond fields.

The engine (`Node`, `Edge`, `Graph`) is copied from picbh. Node keys are the atom
indices of the geometry.
"""

from collections.abc import Sequence
from typing import Literal

import networkx as nx
from pydantic import BaseModel

from ._types import Connectivity


# The engine
class Node(BaseModel):
    """Base class for node data. Subclass to add validated fields."""


class Edge(BaseModel):
    """Base class for edge data. Subclass to add validated fields."""


class Graph[NodeT: Node, EdgeT: Edge](nx.Graph):
    """Generic graph with typed node and edge data."""

    node_type: type[NodeT]
    edge_type: type[EdgeT]

    def validate(self) -> None:
        """Validate node and edge data against the node and edge classes."""
        for *_, data in self.nodes(data=True):
            self.node_type.model_validate(data)

        for *_, data in self.edges(data=True):
            self.edge_type.model_validate(data)


# The fields of a molecular graph
class AtomStereo(BaseModel, frozen=True):
    """The stereo of an atom candidate: its neighbors and its local parity."""

    neighbors: tuple[int, ...]
    parity: bool | None  # None: ill-defined


class BondStereo(BaseModel, frozen=True):
    """The stereo of a bond candidate: the neighbors of each end and its parity."""

    neighbors: tuple[tuple[int, ...], tuple[int, ...]]
    parity: bool | None  # None: ill-defined


class Atom(Node):
    """An atom, with its implicit hydrogens and, for a candidate, its stereo."""

    symbol: str
    hydrogens: int = 0
    stereo: AtomStereo | None = None


class Bond(Edge):
    """A bond, with its change in a TS and, for a candidate, its stereo."""

    change: Literal["breaking", "forming"] | None = None
    stereo: BondStereo | None = None


class MolGraph(Graph[Atom, Bond]):
    """A molecular graph, keyed by atom index."""

    node_type = Atom
    edge_type = Bond


def explicit_graph(symbols: Sequence[str], connectivity: Connectivity) -> MolGraph:
    """The graph with every atom as a node, including the hydrogens."""
    graph = MolGraph()
    for atom, symbol in enumerate(symbols):
        graph.add_node(atom, **dict(Atom(symbol=symbol)))
    for change, pairs in (
        (None, connectivity.bonds),
        ("breaking", connectivity.breaking),
        ("forming", connectivity.forming),
    ):
        for i, j in pairs:
            graph.add_edge(i, j, **dict(Bond(change=change)))
    return graph


def implicit_graph(graph: MolGraph) -> MolGraph:
    """The graph without its non-backbone hydrogens, counted on their neighbors.

    A hydrogen is a backbone hydrogen, and stays a node, if
    - all its neighbors are hydrogens and its index is lower than theirs (one H of
      H2, or a lone H atom),
    - it has more than one neighbor, counting breaking and forming bonds, or
    - it is an atom of a breaking or forming bond.
    """
    implicit = [atom for atom in graph if not _is_backbone_hydrogen(graph, atom)]
    implicit = [atom for atom in implicit if graph.nodes[atom]["symbol"] == "H"]
    result = graph.subgraph(set(graph) - set(implicit)).copy()
    for hydrogen in implicit:
        (neighbor,) = graph[hydrogen]
        result.nodes[neighbor]["hydrogens"] += 1
    return result


def _is_backbone_hydrogen(graph: MolGraph, atom: int) -> bool:
    neighbors = graph[atom]
    if all(graph.nodes[n]["symbol"] == "H" and atom < n for n in neighbors):
        return True
    return len(neighbors) > 1 or any(
        data["change"] is not None for data in neighbors.values()
    )


def connected_components(graph: MolGraph) -> list[frozenset[int]]:
    """The connected components, over all bonds, ordered by their lowest atom."""
    return sorted(map(frozenset, nx.connected_components(graph)), key=min)
