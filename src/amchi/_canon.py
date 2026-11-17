"""The canonicalization: priorities, canonical parities, enantiomer and numbers.

The algorithm is the paper's (Schemes 1 and 3; SI section 1), in the old AutoMol's
order, since the priority values depend on it (W2; decision D22). The data is pulled
out of the graph once per component, and the loops run on plain dicts and tuples.

A reacting component of a TS is read in each direction, as it is and reflected, and
the form with the lowest key is canonical (W4; D28, D29).
"""

import functools
import itertools
import math
from collections.abc import Iterable, Sequence

from pydantic import BaseModel

from ._graph import MolGraph, connected_components
from ._types import Bond, IllDefinedParityError

type _Priorities = dict[int, int]
type _AtomParities = dict[int, bool]
type _BondParities = dict[Bond, bool]
type _Sides = tuple[tuple[int, ...], tuple[int, ...]]
type _Table = tuple[tuple[bool | None, ...], ...]
type _Form = tuple[_Priorities, _AtomParities, _BondParities]  # one stereo loop

# The bond types of the old code: present throughout, forming, breaking
_BOND_TYPE = {None: 1, "forming": 10, "breaking": 90}
_SWAPPED = {None: None, "forming": "breaking", "breaking": "forming"}


class CanonicalForm(BaseModel, frozen=True):
    """The canonical form of one connected component.

    `priorities` and the parities are those of the chosen enantiomer. `numbers` are
    the canonical numbers, 0-based. `is_canonical_enantiomer` is `None` when the
    component has no `m` layer, and `is_canonical_direction` is `None` when it has no
    `r` layer (a component without reacting bonds, or one whose directions are equal).
    In a reacting component, `priorities` are those of the canonical direction, and a
    hydrogen's is negative (`-|p|`).
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
    result = []
    for atoms in connected_components(graph):
        component = _Component(graph, atoms)
        if component.is_reacting:
            reverse = _Component(graph, atoms, reverse=True)
            result.append(_ts_canonical_form(component, reverse))
        else:
            result.append(component.canonical_form())
    return result


def _ts_canonical_form(given: "_Component", reverse: "_Component") -> CanonicalForm:
    """The lowest of up to four forms: each direction, as it is and reflected (D29)."""
    directions = (given, reverse)
    forms = {
        (d, False): c.stereo_loop(c.local_atom_parities)
        for d, c in enumerate(directions)
    }
    if any(atom_parities for _, atom_parities, _ in forms.values()):
        for d, c in enumerate(directions):
            forms[d, True] = c.stereo_loop(_inverted(c.local_atom_parities))
    keys = {(d, m): directions[d].key(*form) for (d, m), form in forms.items()}

    # The best direction as it is and reflected; the given one if they are equal
    best = {m: min((0, 1), key=lambda d, m=m: keys[d, m]) for _, m in keys}
    is_canonical_enantiomer = None
    if True in best:
        ours, theirs = keys[best[False], False], keys[best[True], True]
        if ours != theirs:
            is_canonical_enantiomer = ours < theirs
    mirrored = is_canonical_enantiomer is False
    direction = best[mirrored]
    is_canonical_direction = None
    if keys[0, mirrored] != keys[1, mirrored]:
        is_canonical_direction = direction == 0

    priorities, atom_parities, bond_parities = forms[direction, mirrored]
    return CanonicalForm(
        priorities=priorities,
        numbers=directions[direction].break_ties(
            priorities, atom_parities, bond_parities
        ),
        atom_parities=atom_parities,
        bond_parities=bond_parities,
        is_canonical_enantiomer=is_canonical_enantiomer,
        is_canonical_direction=is_canonical_direction,
    )


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
    """One connected component, pulled out of the graph into plain data.

    A TS component is read in one direction: `reverse` swaps its breaking and forming
    bonds and uses the reverse neighbors of an atom candidate where it has them (D26).
    """

    def __init__(
        self, graph: MolGraph, atoms: frozenset[int], *, reverse: bool = False
    ) -> None:
        self.keys = sorted(atoms)
        self.neighbors = {k: tuple(sorted(graph[k])) for k in self.keys}
        self.is_hydrogen = {k: graph.nodes[k]["symbol"] == "H" for k in self.keys}
        self.change = {
            _bond(i, j): _SWAPPED[change] if reverse else change
            for i, j, change in graph.subgraph(atoms).edges(data="change")
        }
        self.forming = sorted(b for b, c in self.change.items() if c == "forming")
        self.breaking = sorted(b for b, c in self.change.items() if c == "breaking")
        self.is_reacting = bool(self.forming or self.breaking)
        # The bonds of each atom that count for the degree and the bond parities:
        # not the reacting ones
        self.bonds = {
            k: tuple(
                _bond(k, n)
                for n in self.neighbors[k]
                if self.change[_bond(k, n)] is None
            )
            for k in self.keys
        }
        # The type of the bond to each neighbor, in the order of `neighbors`
        self.bond_types = {
            k: tuple(_BOND_TYPE[self.change[_bond(k, n)]] for n in self.neighbors[k])
            for k in self.keys
        }
        self.sorted_bond_types = {
            k: tuple(sorted(t)) for k, t in self.bond_types.items()
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
            if stereo is None:
                continue
            if reverse and stereo.reverse_neighbors is not None:
                self.atom_candidates[k] = stereo.reverse_neighbors
                self.local_atom_parities[k] = stereo.reverse_parity
            else:
                self.atom_candidates[k] = stereo.neighbors
                self.local_atom_parities[k] = stereo.parity
        self.bond_candidates: dict[Bond, _Sides] = {}
        self.local_bond_parities: dict[Bond, _Table] = {}
        for i, j, stereo in graph.subgraph(atoms).edges(data="stereo"):
            if stereo is not None:
                bond = _bond(i, j)
                self.bond_candidates[bond] = stereo.neighbors
                self.local_bond_parities[bond] = stereo.parities

    def canonical_form(self) -> CanonicalForm:
        """Steps 4 to 6 of W2: the stereo loop, the enantiomer, the numbers."""
        local = self.local_atom_parities
        priorities, atom_parities, bond_parities = self.stereo_loop(local)
        is_canonical_enantiomer = None
        if atom_parities:
            reflected = self.stereo_loop(_inverted(local))
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
        """The atom invariants, compared left to right (W2, step 2; W4, step 4).

        A missing parity is `+inf`, so that `False < True < inf`. A reacting component
        adds the bond types and, last, the pairs (bond type, neighbor priority) (D28).
        """
        neighbors = self.neighbors[atom]
        value = (
            *self.fixed[atom],
            atom_parities.get(atom, math.inf),
            tuple(sorted(bond_parities.get(b, math.inf) for b in self.bonds[atom])),
        )
        if not self.is_reacting:
            return (*value, tuple(sorted(priorities[n] for n in neighbors)))
        return (
            *value,
            self.sorted_bond_types[atom],
            tuple(sorted(priorities[n] for n in neighbors)),
            tuple(
                sorted(
                    zip(
                        self.bond_types[atom],
                        (priorities[n] for n in neighbors),
                        strict=True,
                    )
                )
            ),
        )

    def refine(
        self,
        priorities: _Priorities,
        atom_parities: _AtomParities,
        bond_parities: _BondParities,
    ) -> _Priorities:
        """Stable partition refinement, with the old queue order (W2, step 3)."""
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

    def stereo_loop(self, local_atom_parities: dict[int, bool | None]) -> _Form:
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
                local = _defined(local_atom_parities[atom], "atom", atom)
                ordered = sorted(neighbors, key=rank)
                atom_parities[atom] = local ^ _is_odd_permutation(neighbors, ordered)
                assigned = True
            for bond, sides in self.bond_candidates.items():
                if bond in bond_parities or not all(
                    _all_differ(map(rank, side)) for side in sides
                ):
                    continue
                # D27: the table entry of the highest-priority neighbor of each side
                i, j = (side.index(max(side, key=rank)) for side in sides)
                table = self.local_bond_parities[bond]
                bond_parities[bond] = _defined(table[i][j], "bond", bond)
                assigned = True
            if not assigned:
                break
        if self.is_reacting:  # a hydrogen of a TS carries a negative priority
            priorities = {
                k: -abs(p) if self.is_hydrogen[k] else p for k, p in priorities.items()
            }
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

    def key(
        self,
        priorities: _Priorities,
        atom_parities: _AtomParities,
        bond_parities: _BondParities,
    ) -> tuple:
        """The value that picks the canonical form of a TS, lowest first (D29).

        The old direction representation without its parts that depend on the atom
        order, then everything the string says, in canonical numbers (W4, step 5).
        """
        numbers = self.break_ties(priorities, atom_parities, bond_parities)

        def pairs(bonds: Iterable[Bond]) -> list[tuple[int, ...]]:
            return sorted(tuple(sorted(priorities[k] for k in b)) for b in bonds)

        def numbered(bond: Bond) -> tuple[int, ...]:
            return tuple(sorted(numbers[k] for k in bond))

        atoms = sorted(atom_parities, key=priorities.__getitem__)
        bonds = sorted(bond_parities, key=lambda b: sorted(priorities[k] for k in b))
        return (
            (len(self.forming), len(self.breaking)),
            (pairs(self.forming), pairs(self.breaking)),
            tuple(atom_parities[a] for a in atoms),
            tuple(bond_parities[b] for b in bonds),
            tuple(priorities[a] for a in atoms),
            sorted((numbered(b), _BOND_TYPE[c]) for b, c in self.change.items()),
            sorted((numbers[k], self.fixed[k]) for k in self.keys),
            sorted((numbered(b), p) for b, p in bond_parities.items()),
            sorted((numbers[a], p) for a, p in atom_parities.items()),
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


def _inverted(parities: dict[int, bool | None]) -> dict[int, bool | None]:
    """The local atom parities of the mirror image."""
    return {k: None if p is None else not p for k, p in parities.items()}


def _all_differ(values: Iterable[float]) -> bool:
    values = list(values)
    return len(set(values)) == len(values)


def _defined(parity: bool | None, kind: str, key: int | Bond) -> bool:  # noqa: FBT001
    """The local parity of a stereocenter, which must be defined."""
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
