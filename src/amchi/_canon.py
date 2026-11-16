"""The canonicalization: priorities, canonical parities, enantiomer and numbers.

The algorithm is the paper's (Schemes 1 and 3; SI section 1), in the old AutoMol's
order, since the priority values depend on it (W2; decision D22). The data is pulled
out of the graph once per component, and the loops run on plain dicts and tuples.
"""

import functools
import itertools
import math
from collections.abc import Iterable, Mapping, Sequence

from pydantic import BaseModel

from ._graph import MolGraph, connected_components
from ._types import Bond, IllDefinedParityError

type _Priorities = dict[int, int]
type _AtomParities = dict[int, bool]
type _BondParities = dict[Bond, bool]
type _Sides = tuple[tuple[int, ...], tuple[int, ...]]


class CanonicalForm(BaseModel, frozen=True):
    """The canonical form of one connected component.

    `priorities` and the parities are those of the chosen enantiomer. `numbers` are
    the canonical numbers, 0-based. `is_canonical_enantiomer` is `None` when the
    component has no `m` layer, and `is_canonical_direction` is `None` for a species.
    """

    priorities: dict[int, int]
    numbers: dict[int, int]
    atom_parities: dict[int, bool]
    bond_parities: dict[Bond, bool]
    is_canonical_enantiomer: bool | None
    is_canonical_direction: bool | None


def canonical_forms(graph: MolGraph) -> list[CanonicalForm]:
    """The canonical form of each connected component, by its lowest atom.

    The graph is the implicit graph with the candidates' local parities on it.
    """
    return [
        _Component(graph, atoms).canonical_form()
        for atoms in connected_components(graph)
    ]


def _first_pass_parities(graph: MolGraph) -> list[tuple[_AtomParities, _BondParities]]:
    """The canonical parities of each component before the enantiomer choice.

    For the tests, which compare them with the old code's.
    """
    result = []
    for atoms in connected_components(graph):
        component = _Component(graph, atoms)
        _, atom_parities, bond_parities = component.stereo_loop(
            component.local_atom_parities
        )
        result.append((atom_parities, bond_parities))
    return result


class _Component:
    """One connected component, pulled out of the graph into plain data."""

    def __init__(self, graph: MolGraph, atoms: frozenset[int]) -> None:
        self.keys = sorted(atoms)
        self.neighbors = {k: tuple(sorted(graph[k])) for k in self.keys}
        self.is_hydrogen = {k: graph.nodes[k]["symbol"] == "H" for k in self.keys}
        # The bonds of each atom that count for the degree and the bond parities.
        # W4: leave out the reacting bonds of a TS here.
        self.bonds = {
            k: tuple(_bond(k, n) for n in self.neighbors[k]) for k in self.keys
        }
        self.fixed = {
            k: (
                _symbol_key(graph.nodes[k]["symbol"]),
                len(self.bonds[k]),
                graph.nodes[k]["hydrogens"],
            )
            for k in self.keys
        }
        self.atom_candidates: dict[int, tuple[int, ...]] = {}
        self.local_atom_parities: dict[int, bool | None] = {}
        for k in self.keys:
            stereo = graph.nodes[k]["stereo"]
            if stereo is not None:
                self.atom_candidates[k] = stereo.neighbors
                self.local_atom_parities[k] = stereo.parity
        self.bond_candidates: dict[Bond, _Sides] = {}
        self.local_bond_parities: dict[Bond, bool | None] = {}
        for i, j, stereo in graph.subgraph(atoms).edges(data="stereo"):
            if stereo is not None:
                bond = _bond(i, j)
                self.bond_candidates[bond] = stereo.neighbors
                self.local_bond_parities[bond] = stereo.parity

    def canonical_form(self) -> CanonicalForm:
        """Steps 4 to 6 of W2: the stereo loop, the enantiomer, the numbers."""
        local = self.local_atom_parities
        priorities, atom_parities, bond_parities = self.stereo_loop(local)
        # W4: the direction step of a TS goes here.
        is_canonical_enantiomer = None
        if atom_parities:
            mirror = {k: None if p is None else not p for k, p in local.items()}
            reflected = self.stereo_loop(mirror)
            original = self.representation(priorities, atom_parities, bond_parities)
            other = self.representation(*reflected)
            if original < other:
                is_canonical_enantiomer = True
            elif original > other:
                is_canonical_enantiomer = False
                priorities, atom_parities, bond_parities = reflected
        numbers = self.break_ties(priorities, atom_parities, bond_parities)
        return CanonicalForm(
            priorities=priorities,
            numbers=numbers,
            atom_parities=atom_parities,
            bond_parities=bond_parities,
            is_canonical_enantiomer=is_canonical_enantiomer,
            is_canonical_direction=None,
        )

    def sort_value(
        self,
        atom: int,
        priorities: _Priorities,
        atom_parities: _AtomParities,
        bond_parities: _BondParities,
    ) -> tuple:
        """The atom invariants, compared left to right (W2, step 2).

        A missing parity is `+inf`, so that `False < True < inf`.
        """
        return (
            *self.fixed[atom],
            atom_parities.get(atom, math.inf),
            tuple(sorted(bond_parities.get(b, math.inf) for b in self.bonds[atom])),
            # W4: the bond types of a TS go here.
            tuple(sorted(priorities[n] for n in self.neighbors[atom])),
        )

    def refine(
        self,
        priorities: _Priorities,
        atom_parities: _AtomParities,
        bond_parities: _BondParities,
    ) -> _Priorities:
        """Stable partition refinement, with the old queue order (W2, step 3)."""
        # W4: backbone hydrogens of a TS get their sign back (-|p|) at the end.
        priorities = {k: abs(p) for k, p in priorities.items()}
        classes = _classes(self.keys, priorities)
        queue = sorted(
            ((i, atoms) for i, atoms in classes.items() if len(atoms) > 1),
            reverse=True,
        )
        while queue:
            index, atoms = queue.pop(0)
            values = {
                k: self.sort_value(k, priorities, atom_parities, bond_parities)
                for k in atoms
            }
            atoms = sorted(atoms, key=values.__getitem__)
            parts = [
                tuple(part)
                for _, part in itertools.groupby(atoms, key=values.__getitem__)
            ]
            if len(parts) == 1:
                continue
            new = []
            i = index
            for part in parts:
                classes[i] = part
                for k in part:
                    priorities[k] = i
                new.append(i)
                i += len(part)
            touched = {
                priorities[n]
                for i in new
                for k in classes[i]
                for n in self.neighbors[k]
            }
            touched -= {i for i, _ in queue}
            for i in sorted(touched):
                if len(classes[i]) > 1:
                    queue.insert(0, (i, classes[i]))
        return priorities

    def stereo_priority(self, priorities: _Priorities, atom: int) -> float:
        """The priority of a stereo-determining neighbor; hydrogens rank lowest."""
        if atom not in priorities:
            return -math.inf  # an implicit hydrogen
        if self.is_hydrogen[atom]:
            return -abs(priorities[atom])  # a backbone hydrogen
        return priorities[atom]

    def stereo_loop(
        self, local_atom_parities: dict[int, bool | None]
    ) -> tuple[_Priorities, _AtomParities, _BondParities]:
        """Refine, assign the parities of new stereocenters, and repeat (W2, step 4)."""
        priorities: _Priorities = dict.fromkeys(self.keys, 0)
        atom_parities: _AtomParities = {}
        bond_parities: _BondParities = {}
        previous = None
        while priorities != previous:
            previous = priorities
            priorities = self.refine(priorities, atom_parities, bond_parities)

            rank = functools.partial(self.stereo_priority, priorities)
            assigned = False
            for atom, neighbors in self.atom_candidates.items():
                if atom in atom_parities or not _all_differ(map(rank, neighbors)):
                    continue
                local = _defined(local_atom_parities, atom, "atom")
                ordered = sorted(neighbors, key=rank)
                atom_parities[atom] = local ^ _is_odd_permutation(neighbors, ordered)
                assigned = True
            for bond, sides in self.bond_candidates.items():
                if bond in bond_parities or not all(
                    _all_differ(map(rank, side)) for side in sides
                ):
                    continue
                local = _defined(self.local_bond_parities, bond, "bond")
                for side in sides:
                    local ^= max(side, key=rank) != side[-1]
                bond_parities[bond] = local
                assigned = True
            if not assigned:
                break
        return priorities, atom_parities, bond_parities

    def representation(
        self,
        priorities: _Priorities,
        atom_parities: _AtomParities,
        bond_parities: _BondParities,
    ) -> tuple:
        """The value that picks the canonical enantiomer, lowest first (W2, step 5)."""
        atoms = sorted(atom_parities, key=priorities.__getitem__)
        bonds = sorted(bond_parities, key=lambda b: sorted(priorities[k] for k in b))

        def ranks(neighbors: Sequence[int]) -> tuple[float, ...]:
            return tuple(sorted(priorities.get(n, -math.inf) for n in neighbors))

        return (
            tuple(atom_parities[a] for a in atoms),
            tuple(bond_parities[b] for b in bonds),
            tuple(priorities[a] for a in atoms),
            tuple(ranks(self.atom_candidates[a]) for a in atoms),
            tuple(tuple(map(ranks, self.bond_candidates[b])) for b in bonds),
        )

    def break_ties(
        self,
        priorities: _Priorities,
        atom_parities: _AtomParities,
        bond_parities: _BondParities,
    ) -> dict[int, int]:
        """The canonical numbers: promote the last atom of the highest tied class."""
        priorities = {k: abs(p) for k, p in priorities.items()}
        while True:
            tied = {
                i: atoms
                for i, atoms in _classes(self.keys, priorities).items()
                if len(atoms) > 1
            }
            if not tied:
                return priorities
            index = max(tied)
            priorities[tied[index][-1]] = index + len(tied[index]) - 1
            priorities = self.refine(priorities, atom_parities, bond_parities)


def _bond(i: int, j: int) -> Bond:
    return (i, j) if i < j else (j, i)


def _symbol_key(symbol: str) -> str:
    """C first, then alphabetical, and a backbone H last."""
    return "" if symbol == "C" else "zz" if symbol == "H" else symbol


def _classes(
    keys: Sequence[int], priorities: _Priorities
) -> dict[int, tuple[int, ...]]:
    """The atoms of each priority, by ascending atom index."""
    classes: dict[int, list[int]] = {}
    for k in keys:
        classes.setdefault(priorities[k], []).append(k)
    return {i: tuple(atoms) for i, atoms in classes.items()}


def _all_differ(values: Iterable[float]) -> bool:
    values = list(values)
    return len(set(values)) == len(values)


def _defined[K](parities: Mapping[K, bool | None], key: K, kind: str) -> bool:
    """The local parity of a stereocenter, which must be defined."""
    parity = parities[key]
    if parity is None:
        msg = f"The parity of the stereocenter at {kind} {key} is ill-defined"
        raise IllDefinedParityError(msg)
    return parity


def _is_odd_permutation(sequence: Sequence[int], target: Sequence[int]) -> bool:
    """Does it take an odd number of swaps to turn `sequence` into `target`?"""
    permutation = [target.index(k) for k in sequence]
    odd = False
    for i in range(len(permutation)):
        while permutation[i] != i:
            j = permutation[i]
            permutation[i], permutation[j] = permutation[j], permutation[i]
            odd = not odd
    return odd
