"""The π-system model: which bonds are rigid and planar.

The model is the old AutoMol's: neutral atoms with standard valences, and dominant
resonance structures that pair as many unsaturated neighbors as possible. Instead of
listing the structures, each question about a bond or an atom is a search for one
maximum pairing that contains given π pairs (W3; decisions D12 and D23).

A *pairing* puts π pairs on the bonds of a π system: at most `u` at an atom of
unsaturation `u`, and at most `min(u_a, u_b, 2)` on a bond. A *maximum pairing* has
as many pairs as possible, and its *deficiency* is the number of electrons it leaves
unpaired.
"""

from collections.abc import Iterable
from itertools import combinations_with_replacement

import networkx as nx

from ._element import lone_pairs, valence
from ._graph import MolGraph
from ._types import Bond

# valence + lone pairs of an atom that can be sp2 (C, N, O, F and the elements below)
_OCTET = 4


def unsaturations(graph: MolGraph) -> dict[int, int]:
    """The unsaturation of each atom: its valence minus its sigma bonds, at least 0.

    The sigma bonds are the graph's bonds and the implicit hydrogens, so the explicit
    and the implicit graph give the same values for their common atoms.
    """
    return {
        atom: max(
            0, valence(symbol) - graph.degree(atom) - graph.nodes[atom]["hydrogens"]
        )
        for atom, symbol in graph.nodes(data="symbol")
    }


def pi_systems(graph: MolGraph) -> list[frozenset[int]]:
    """The π systems: connected unsaturated atoms, more than one, by lowest atom."""
    unsaturated = [atom for atom, u in unsaturations(graph).items() if u > 0]
    components = nx.connected_components(graph.subgraph(unsaturated))
    return sorted((frozenset(c) for c in components if len(c) > 1), key=min)


def rigid_planar_bonds(
    graph: MolGraph, bonds: Iterable[Bond] | None = None
) -> frozenset[Bond]:
    """The rigid and planar bonds among `bonds` (default: all bonds of the graph).

    A bond is rigid and planar when it can carry a π pair in some maximum pairing, and
    neither of its atoms can carry two π pairs in any (D23). Only bonds between two
    atoms of the same π system, whose elements have valence + lone pairs = 4, qualify.
    """
    if bonds is None:
        bonds = (tuple(sorted(bond)) for bond in graph.edges)
    unsaturation = unsaturations(graph)
    system_of = {atom: atoms for atoms in pi_systems(graph) for atom in atoms}
    octet = {
        atom: valence(symbol) + lone_pairs(symbol) == _OCTET
        for atom, symbol in graph.nodes(data="symbol")
    }
    searches: dict[frozenset[int], _PiSystem] = {}
    can_have_two: dict[int, bool] = {}
    result = set()
    for i, j in bonds:
        atoms = system_of.get(i)
        if atoms is None or j not in atoms or not (octet[i] and octet[j]):
            continue
        if atoms not in searches:
            searches[atoms] = _PiSystem(graph, atoms, unsaturation)
        search = searches[atoms]
        if not search.can_be_double(i, j):
            continue
        for atom in (i, j):
            if atom not in can_have_two:
                can_have_two[atom] = search.can_have_two(atom)
        if not (can_have_two[i] or can_have_two[j]):
            result.add((i, j))
    return frozenset(result)


class _PiSystem:
    """The pairing questions about one π system."""

    def __init__(
        self, graph: MolGraph, atoms: frozenset[int], unsaturation: dict[int, int]
    ) -> None:
        self.capacity = {atom: unsaturation[atom] for atom in atoms}
        self.neighbors = {
            atom: sorted(n for n in graph[atom] if n in atoms) for atom in atoms
        }
        self.multiplicity = {
            _key(a, b): min(unsaturation[a], unsaturation[b], 2)
            for a in atoms
            for b in self.neighbors[a]
        }
        deficiency = sum(self.capacity.values()) % 2
        while not _has_pairing(
            self.capacity, self.multiplicity, self.neighbors, deficiency
        ):
            deficiency += 2
        self.deficiency = deficiency

    def can_be_double(self, a: int, b: int) -> bool:
        """Can the bond carry a π pair in some maximum pairing?"""
        return self._with([(a, b)])

    def can_have_two(self, atom: int) -> bool:
        """Can the atom carry two π pairs in some maximum pairing?"""
        if self.capacity[atom] < 2:  # noqa: PLR2004
            return False
        return any(
            self._with([(atom, b), (atom, c)])
            for b, c in combinations_with_replacement(self.neighbors[atom], 2)
        )

    def _with(self, pairs: list[tuple[int, int]]) -> bool:
        """Can these π pairs be part of a maximum pairing?"""
        capacity, multiplicity = dict(self.capacity), dict(self.multiplicity)
        for a, b in pairs:
            key = _key(a, b)
            capacity[a] -= 1
            capacity[b] -= 1
            multiplicity[key] -= 1
            if capacity[a] < 0 or capacity[b] < 0 or multiplicity[key] < 0:
                return False
        return _has_pairing(capacity, multiplicity, self.neighbors, self.deficiency)


def _key(a: int, b: int) -> Bond:
    return (a, b) if a < b else (b, a)


def _has_pairing(
    capacity: dict[int, int],
    multiplicity: dict[Bond, int],
    neighbors: dict[int, list[int]],
    budget: int,
) -> bool:
    """Is there a pairing of the remaining capacities with at most `budget` unpaired?

    Pairs the most constrained atom first (the fewest available partners), and stops
    at the first success. Mutates the two dicts and restores them before returning.
    """
    atom, partners = None, None
    for a, c in capacity.items():
        if c > 0:
            options = [
                b
                for b in neighbors[a]
                if capacity[b] > 0 and multiplicity[_key(a, b)] > 0
            ]
            if partners is None or len(options) < len(partners):
                atom, partners = a, options
                if not options:
                    break
    if atom is None or partners is None:
        return True
    for b in partners:
        key = _key(atom, b)
        capacity[atom] -= 1
        capacity[b] -= 1
        multiplicity[key] -= 1
        found = _has_pairing(capacity, multiplicity, neighbors, budget)
        capacity[atom] += 1
        capacity[b] += 1
        multiplicity[key] += 1
        if found:
            return True
    if budget > 0:
        capacity[atom] -= 1
        found = _has_pairing(capacity, multiplicity, neighbors, budget - 1)
        capacity[atom] += 1
        return found
    return False
