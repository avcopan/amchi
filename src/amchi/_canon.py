"""The canonicalization: priorities, canonical parities, enantiomer and numbers.

The algorithm is the paper's (Schemes 1 and 3; SI section 1), in the old AutoMol's
order, since the priority values depend on it (W2; decision D22). The data is pulled
out of the graph once per component, and the loops run on plain dicts and tuples.

A reacting component of a TS is read in each direction, as it is and reflected, and
the form with the lowest key is canonical (W4; D28, D29).

A para-stereocenter is a candidate whose neighbors keep one tied pair joined in a
ring through it, as in cis- and trans-1,3-dimethylcyclobutane. Those of a ring system
are kept if flipping one changes the form, and add Schneider's dependent-chirality
invariant to the refinement; their parities come from the tie-breaking if they are
still tied (W6; D32, D33). A component without one takes the old path.
"""

import functools
import itertools
import math
from collections.abc import Callable, Iterable, Sequence
from typing import NamedTuple, cast

import networkx as nx
from pydantic import BaseModel

from ._graph import MolGraph, connected_components
from ._types import Bond, IllDefinedParityError

type _Priorities = dict[int, int]
type _AtomParities = dict[int, bool]
type _BondParities = dict[Bond, bool]
type _Sides = tuple[tuple[int, ...], tuple[int, ...]]
type _Table = tuple[tuple[bool | None, ...], ...]
type _Form = tuple[_Priorities, _AtomParities, _BondParities]  # one stereo loop
type _Pair = tuple[int, int]
type _Member = int | Bond  # a para atom or a para bond

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
    """The lowest of up to four forms: each direction, as it is and reflected (D29).

    Each direction perceives its para-stereocenters as it is, and keeps them for its
    reflected form (W6, rule 5).
    """
    directions = (given, reverse)
    results = {}
    for d, c in enumerate(directions):
        local = c.local_atom_parities
        form = c.stereo_loop(local)
        results[d, False] = c.full_form(
            local, None if c.perceive(local, form) else form
        )
    if any(atom_parities for _, _, (_, atom_parities, _) in results.values()):
        for d, c in enumerate(directions):
            results[d, True] = c.full_form(_inverted(c.local_atom_parities))
    keys = {dm: key for dm, (key, _, _) in results.items()}

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

    _, (priorities, _, _), (numbers, atom_parities, bond_parities) = results[
        direction, mirrored
    ]
    return CanonicalForm(
        priorities=priorities,
        numbers=numbers,
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
        # The para-stereocenters (W6): each atom's tied pair, and each bond's tied
        # pair per side (None where its neighbors all differ). `viewers` gives the
        # atoms and bonds that each member of a pair reads a view of.
        self.para_atoms: dict[int, _Pair] = {}
        self.para_bonds: dict[Bond, tuple[_Pair | None, _Pair | None]] = {}
        self.viewers: dict[int, list[tuple[_Member, int]]] = {}
        # The viewers whose views read each atom's priority, which refinement must
        # requeue when that atom's class splits
        self.dependents: dict[int, set[int]] = {}
        # The local atom parities of the current pass, which `stereo_loop` sets
        self.pass_atom_parities = self.local_atom_parities

    def canonical_form(self) -> CanonicalForm:
        """Steps 4 to 6 of W2: the stereo loop, the enantiomer, the numbers."""
        local = self.local_atom_parities
        form = self.stereo_loop(local)
        if self.perceive(local, form):
            return self.para_canonical_form()
        priorities, atom_parities, bond_parities = form
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
        numbers, _, _ = self.break_ties(priorities, atom_parities, bond_parities)
        return CanonicalForm(
            priorities=priorities,
            numbers=numbers,
            atom_parities=atom_parities,
            bond_parities=bond_parities,
            is_canonical_enantiomer=is_canonical_enantiomer,
            is_canonical_direction=None,
        )

    def para_canonical_form(self) -> CanonicalForm:
        """The canonical form with para-stereocenters (W6, rule 5).

        The enantiomer is chosen on the full key, since a para parity may exist only
        in canonical numbers.
        """
        local = self.local_atom_parities
        key, (priorities, _, _), numbered = self.full_form(local)
        is_canonical_enantiomer = None
        if numbered[1]:
            other, (other_priorities, _, _), other_numbered = self.full_form(
                _inverted(local)
            )
            if key < other:
                is_canonical_enantiomer = True
            elif key > other:
                is_canonical_enantiomer = False
                priorities, numbered = other_priorities, other_numbered
        numbers, atom_parities, bond_parities = numbered
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
        A component with para-stereocenters adds, last, the views of the atom (W6).
        """
        neighbors = self.neighbors[atom]
        value = (
            *self.fixed[atom],
            atom_parities.get(atom, math.inf),
            tuple(sorted(bond_parities.get(b, math.inf) for b in self.bonds[atom])),
        )
        if not self.is_reacting:
            value = (*value, tuple(sorted(priorities[n] for n in neighbors)))
        else:
            value = (
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
        if self.viewers:
            return (*value, self.views(atom, priorities))
        return value

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
            if self.dependents:
                touched |= {
                    priorities[v]
                    for i in new
                    for k in classes[i]
                    for v in self.dependents.get(k, ())
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
        self.pass_atom_parities = local_atom_parities
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
        """The value that picks the canonical enantiomer, lowest first (W2, step 5).

        Nothing in it depends on the atom order: tied stereocenters are ordered by
        parity, and the two sides of a bond by their ranks.
        """
        atoms, bonds = _stereocenters(priorities, atom_parities, bond_parities)

        def ranks(neighbors: Sequence[int]) -> tuple[float, ...]:
            return tuple(sorted(priorities.get(n, -math.inf) for n in neighbors))

        return (
            tuple(atom_parities[a] for a in atoms),
            tuple(bond_parities[b] for b in bonds),
            tuple(priorities[a] for a in atoms),
            tuple(ranks(self.atom_candidates[a]) for a in atoms),
            tuple(tuple(sorted(map(ranks, self.bond_candidates[b]))) for b in bonds),
        )

    def key(self, form: _Form, numbered: _Form) -> tuple:
        """The value that picks the canonical form of a TS, lowest first (D29).

        The old direction representation without its parts that depend on the atom
        order, then everything the string says, in canonical numbers (W4, step 5).
        `form` is the stereo loop's, and `numbered` the tie-breaking's, whose parities
        include those of para-stereocenters.
        """
        priorities, atom_parities, bond_parities = form
        numbers, numbered_atom_parities, numbered_bond_parities = numbered

        def pairs(bonds: Iterable[Bond]) -> list[tuple[int, ...]]:
            return sorted(tuple(sorted(priorities[k] for k in b)) for b in bonds)

        def numbered(bond: Bond) -> tuple[int, ...]:
            return tuple(sorted(numbers[k] for k in bond))

        atoms, bonds = _stereocenters(priorities, atom_parities, bond_parities)
        return (
            (len(self.forming), len(self.breaking)),
            (pairs(self.forming), pairs(self.breaking)),
            tuple(atom_parities[a] for a in atoms),
            tuple(bond_parities[b] for b in bonds),
            tuple(priorities[a] for a in atoms),
            sorted((numbered(b), _BOND_TYPE[c]) for b, c in self.change.items()),
            sorted((numbers[k], self.fixed[k]) for k in self.keys),
            sorted((numbered(b), p) for b, p in numbered_bond_parities.items()),
            sorted((numbers[a], p) for a, p in numbered_atom_parities.items()),
        )

    def break_ties(
        self,
        priorities: _Priorities,
        atom_parities: _AtomParities,
        bond_parities: _BondParities,
    ) -> _Form:
        """The canonical numbers: the lowest leaf of a search over the tie-breaks.

        Returns the numbers and the parities. Each step promotes one atom of the
        highest tied class and refines. Promoting the last one, as the old code does,
        is canonical only if the tied atoms are symmetric, which they needn't be in a
        cage (`stjohn:138545_0`). So each atom is tried in turn, and the leaf with the
        lowest `leaf_key` is kept (`_TieSearch`). The old choice is tried first and
        kept on a tie, so the numbers change only where the old ones depended on the
        atom order. A para-stereocenter gets its parity as soon as its neighbors
        differ, and the priorities are refined again (W6, rule 3).
        """
        priorities = {k: abs(p) for k, p in priorities.items()}
        search = _TieSearch(self)
        search.node((), priorities, dict(atom_parities), dict(bond_parities))
        return cast("_Leaf", search.best).numbered

    def settle(
        self,
        priorities: _Priorities,
        atom_parities: _AtomParities,
        bond_parities: _BondParities,
    ) -> _Priorities:
        """Assign the para parities whose neighbors differ, and refine, until none."""
        while self.viewers and self.assign_para(
            priorities, atom_parities, bond_parities
        ):
            priorities = self.refine(priorities, atom_parities, bond_parities)
        return priorities

    def leaf_key(self, numbered: _Form) -> tuple:
        """Everything the string says about the component, in canonical numbers."""
        numbers, atom_parities, bond_parities = numbered

        def bond(b: Bond) -> tuple[int, ...]:
            return tuple(sorted(numbers[k] for k in b))

        return (
            sorted((bond(b), _BOND_TYPE[c]) for b, c in self.change.items()),
            sorted((numbers[k], self.fixed[k]) for k in self.keys),
            sorted((bond(b), p) for b, p in bond_parities.items()),
            sorted((numbers[a], p) for a, p in atom_parities.items()),
        )

    # Para-stereocenters (W6) ---------------------------------------------------------
    def full_form(
        self, local_atom_parities: dict[int, bool | None], form: _Form | None = None
    ) -> tuple[tuple, _Form, _Form]:
        """One pass: the key, the stereo loop's form and the tie-breaking's (rule 4).

        `form` is the stereo loop's result, if it is known already.
        """
        if form is None:
            form = self.stereo_loop(local_atom_parities)
        self.pass_atom_parities = local_atom_parities
        numbered = self.break_ties(*form)
        if self.is_reacting:
            return self.key(form, numbered), form, numbered
        numbers, atom_parities, bond_parities = numbered
        key = (
            *self.representation(*form),
            sorted(tuple(sorted(numbers[k] for k in b)) for b in self.change),
            sorted((numbers[k], self.fixed[k]) for k in self.keys),
            sorted(
                (tuple(sorted(numbers[k] for k in b)), p)
                for b, p in bond_parities.items()
            ),
            sorted((numbers[a], p) for a, p in atom_parities.items()),
        )
        return key, form, numbered

    def perceive(
        self, local_atom_parities: dict[int, bool | None], form: _Form
    ) -> bool:
        """Find the para-stereocenters from the first stereo loop's `form` (rules 1, 4).

        The candidates of a ring system are kept together if flipping any one alone
        changes the form as it is, and dropped together otherwise. An ill-defined
        parity is tested as `True`; if its candidate is kept, the final pass raises.
        Returns whether any are kept.
        """
        atoms, bonds = self.para_candidates(form)
        if not (atoms or bonds):
            return False
        local = dict(local_atom_parities)
        for atom in atoms:
            if local[atom] is None:
                local[atom] = True
        tables = dict(self.local_bond_parities)
        for bond in bonds:
            self.local_bond_parities[bond] = _replaced_none(tables[bond])
        try:
            self.set_para(atoms, bonds)
            reference = self.full_form(local)[0]
            for group in self.groups():
                if all(self.flipped_key(local, m) == reference for m in group):
                    for member in group:
                        if isinstance(member, int):
                            del atoms[member]
                        else:
                            del bonds[member]
        finally:
            self.local_bond_parities = tables
        self.set_para(atoms, bonds)
        return bool(self.viewers)

    def para_candidates(
        self, form: _Form
    ) -> tuple[dict[int, _Pair], dict[Bond, tuple[_Pair | None, _Pair | None]]]:
        """The candidates without a parity that keep a tied pair in a ring (rule 1)."""
        priorities, atom_parities, bond_parities = form
        rank = functools.partial(self.stereo_priority, priorities)
        atoms = {}
        for atom, neighbors in self.atom_candidates.items():
            if atom not in atom_parities:
                pair = self.tied_pair(rank, neighbors, atom)
                if pair:
                    atoms[atom] = pair
        bonds = {}
        for bond, sides in self.bond_candidates.items():
            if bond not in bond_parities:
                pairs = [
                    self.tied_pair(rank, side, end)
                    for side, end in zip(sides, bond, strict=True)
                ]
                if None not in pairs and any(pairs) and self.splits(bond, pairs):
                    bonds[bond] = (pairs[0] or None, pairs[1] or None)
        return atoms, bonds

    def splits(self, bond: Bond, pairs: list[_Pair | tuple[()] | None]) -> bool:
        """Does each tied pair of a bond lie on both sides of it?

        Its two atoms must read opposite entries against every neighbor across the
        bond. Otherwise (a ring at 90° to the plane) the geometry can't tell them
        apart, and the bond is no para-stereocenter.
        """
        sides, table = self.bond_candidates[bond], self.local_bond_parities[bond]
        rows = (table, tuple(zip(*table, strict=True)))
        for side, pair in enumerate(pairs):
            if pair:
                a, b = (rows[side][sides[side].index(k)] for k in pair)
                if any(x == y for x, y in zip(a, b, strict=True)):
                    return False
        return True

    def tied_pair(
        self, rank: Callable[[int], float], neighbors: Sequence[int], center: int
    ) -> _Pair | tuple[()] | None:
        """The one tied pair of `neighbors`, if it is joined in a ring through `center`.

        `()` if all differ, and `None` otherwise: two tied pairs, three tied, or a
        pair joined only through `center`, whose swap is a symmetry.
        """
        ranks = [rank(n) for n in neighbors]
        if _all_differ(ranks):
            return ()
        if len(set(ranks)) != len(ranks) - 1:
            return None
        a, b = (n for n, r in zip(neighbors, ranks, strict=True) if ranks.count(r) > 1)
        return (a, b) if self.joined(a, b, center) else None

    def joined(self, a: int, b: int, center: int) -> bool:
        """Are `a` and `b` connected without passing through `center`?"""
        if a not in self.neighbors or b not in self.neighbors:
            return False  # an implicit hydrogen
        seen, stack = {a, center}, [a]
        while stack:
            for n in self.neighbors[stack.pop()]:
                if n == b:
                    return True
                if n not in seen:
                    seen.add(n)
                    stack.append(n)
        return False

    def set_para(
        self,
        atoms: dict[int, _Pair],
        bonds: dict[Bond, tuple[_Pair | None, _Pair | None]],
    ) -> None:
        """Set the para-stereocenters and the atoms that read their views."""
        self.para_atoms, self.para_bonds = atoms, bonds
        self.viewers, self.dependents = {}, {}
        for atom, pair in atoms.items():
            for viewer in pair:
                self.viewers.setdefault(viewer, []).append((atom, -1))
                for k in (atom, *self.atom_candidates[atom]):
                    self.dependents.setdefault(k, set()).add(viewer)
        for bond, pairs in bonds.items():
            for side, pair in enumerate(pairs):
                for viewer in pair or ():
                    self.viewers.setdefault(viewer, []).append((bond, side))
                    for k in (*bond, *self.bond_candidates[bond][1 - side]):
                        self.dependents.setdefault(k, set()).add(viewer)

    def groups(self) -> list[set[_Member]]:
        """The para-stereocenters of each ring system (block), merged where they meet.

        A bond with a pair on both sides joins the blocks of its two sides.
        """
        graph = nx.Graph((k, n) for k in self.keys for n in self.neighbors[k])
        block = {}
        for i, edges in enumerate(nx.biconnected_component_edges(graph)):
            for u, v in edges:
                block[u, v] = block[v, u] = i
        members: list[tuple[_Member, set[int]]] = [
            (atom, {block[atom, pair[0]]}) for atom, pair in self.para_atoms.items()
        ]
        members += [
            (
                bond,
                {
                    block[end, pair[0]]
                    for end, pair in zip(bond, pairs, strict=True)
                    if pair
                },
            )
            for bond, pairs in self.para_bonds.items()
        ]
        groups: list[tuple[set[_Member], set[int]]] = []
        for member, member_blocks in members:
            group, blocks = {member}, set(member_blocks)
            for other in [g for g in groups if g[1] & blocks]:
                groups.remove(other)
                group |= other[0]
                blocks |= other[1]
            groups.append((group, blocks))
        return [group for group, _ in groups]

    def flipped_key(
        self, local_atom_parities: dict[int, bool | None], member: _Member
    ) -> tuple:
        """The key with one para-stereocenter's local stereo flipped."""
        if isinstance(member, int):
            flipped = {**local_atom_parities, member: not local_atom_parities[member]}
            return self.full_form(flipped)[0]
        table = self.local_bond_parities[member]
        self.local_bond_parities[member] = _flipped(table)
        try:
            return self.full_form(local_atom_parities)[0]
        finally:
            self.local_bond_parities[member] = table

    def views(self, atom: int, priorities: _Priorities) -> tuple:
        """The para parities read with `atom` first (Schneider's invariant, rule 2).

        A view is left out while the other neighbors don't all differ, since the
        stored neighbor order would decide it.
        """
        rank = functools.partial(self.stereo_priority, priorities)
        views = []
        for member, side in self.viewers.get(atom, ()):
            if isinstance(member, int):
                view = self.atom_view(member, atom, rank)
                if view is not None:
                    views.append((0, (priorities[member],), view))
            else:
                view = self.bond_view(member, side, atom, rank)
                if view is not None:
                    priority = tuple(sorted(priorities[k] for k in member))
                    views.append((1, priority, view))
        return tuple(sorted(views))

    def atom_view(
        self, center: int, viewer: int, rank: Callable[[int], float]
    ) -> bool | None:
        """The local parity of `center` against (`viewer`, the rest by priority)."""
        local = self.pass_atom_parities[center]
        neighbors = self.atom_candidates[center]
        rest = sorted((n for n in neighbors if n != viewer), key=rank)
        if local is None or not _all_differ(map(rank, rest)):
            return None
        return local ^ _is_odd_permutation(neighbors, (viewer, *rest))

    def bond_view(
        self, bond: Bond, side: int, viewer: int, rank: Callable[[int], float]
    ) -> bool | None:
        """The table entry of `viewer` and the top neighbor of the other side (D27)."""
        sides = self.bond_candidates[bond]
        other = sides[1 - side]
        if not _all_differ(map(rank, other)):
            return None
        i = sides[side].index(viewer)
        j = other.index(max(other, key=rank))
        table = self.local_bond_parities[bond]
        return table[i][j] if side == 0 else table[j][i]

    def assign_para(
        self,
        priorities: _Priorities,
        atom_parities: _AtomParities,
        bond_parities: _BondParities,
    ) -> bool:
        """Give a parity to each para-stereocenter whose neighbors now all differ."""
        rank = functools.partial(self.stereo_priority, priorities)
        assigned = False
        for atom in self.para_atoms:
            neighbors = self.atom_candidates[atom]
            if atom in atom_parities or not _all_differ(map(rank, neighbors)):
                continue
            local = _defined(self.pass_atom_parities[atom], "atom", atom)
            ordered = sorted(neighbors, key=rank)
            atom_parities[atom] = local ^ _is_odd_permutation(neighbors, ordered)
            assigned = True
        for bond in self.para_bonds:
            sides = self.bond_candidates[bond]
            if bond in bond_parities or not all(
                _all_differ(map(rank, side)) for side in sides
            ):
                continue
            i, j = (side.index(max(side, key=rank)) for side in sides)
            table = self.local_bond_parities[bond]
            bond_parities[bond] = _defined(table[i][j], "bond", bond)
            assigned = True
        return assigned


class _Leaf(NamedTuple):
    key: tuple
    path: tuple[int, ...]  # the promoted atoms, in order
    numbered: _Form


class _TieSearch:
    """The search tree of the tie-breaking, pruned by the automorphisms it finds.

    As in nauty: two leaves with equal keys give an automorphism, which maps the
    first one's path to the other's. The subtree where they diverge is then a copy of
    one already searched, and is left. At each node, an atom in the orbit of one
    already tried, under the automorphisms that fix the node's path, is skipped.
    """

    def __init__(self, component: "_Component") -> None:
        self.component = component
        self.first: _Leaf | None = None
        self.best: _Leaf | None = None
        self.automorphisms: list[dict[int, int]] = []

    def node(
        self,
        path: tuple[int, ...],
        priorities: _Priorities,
        atom_parities: _AtomParities,
        bond_parities: _BondParities,
    ) -> int | None:
        """Search below a node; return the depth to go back to, if a copy was found."""
        c = self.component
        priorities = c.settle(priorities, atom_parities, bond_parities)
        tied = {
            i: atoms
            for i, atoms in _classes(c.keys, priorities).items()
            if len(atoms) > 1
        }
        if not tied:
            return self.leaf(path, (priorities, atom_parities, bond_parities))
        index = max(tied)
        tried: list[int] = []
        for atom in reversed(tied[index]):  # the old choice first
            if tried and self.is_equivalent(atom, tried, path):
                continue
            tried.append(atom)
            child = dict(priorities)
            child[atom] = index + len(tied[index]) - 1
            child_atom_parities, child_bond_parities = (
                dict(atom_parities),
                dict(bond_parities),
            )
            child = c.refine(child, child_atom_parities, child_bond_parities)
            depth = self.node(
                (*path, atom), child, child_atom_parities, child_bond_parities
            )
            if depth is not None and depth < len(path):
                return depth
        return None

    def leaf(self, path: tuple[int, ...], numbered: _Form) -> int | None:
        leaf = _Leaf(self.component.leaf_key(numbered), path, numbered)
        if self.first is None or self.best is None:
            self.first = self.best = leaf
            return None
        for other in (self.first, self.best):
            if leaf.key == other.key:
                atom = {n: k for k, n in leaf.numbered[0].items()}
                self.automorphisms.append(
                    {k: atom[n] for k, n in other.numbered[0].items()}
                )
                return next(
                    (
                        d
                        for d, (a, b) in enumerate(zip(path, other.path, strict=False))
                        if a != b
                    ),
                    min(len(path), len(other.path)),
                )
        if leaf.key < self.best.key:
            self.best = leaf
        return None

    def is_equivalent(self, atom: int, tried: list[int], path: tuple[int, ...]) -> bool:
        """Is `atom` in the orbit of a tried atom, under automorphisms fixing `path`?"""
        parent = {k: k for k in self.component.keys}

        def find(k: int) -> int:
            while parent[k] != k:
                parent[k] = parent[parent[k]]
                k = parent[k]
            return k

        for g in self.automorphisms:
            if all(g[k] == k for k in path):
                for k, image in g.items():
                    parent[find(k)] = find(image)
        return find(atom) in {find(k) for k in tried}


def _stereocenters(
    priorities: _Priorities, atom_parities: _AtomParities, bond_parities: _BondParities
) -> tuple[list[int], list[Bond]]:
    """The stereo atoms and bonds by priority, ties by parity, never by atom index."""
    atoms = sorted(atom_parities, key=lambda a: (priorities[a], atom_parities[a]))
    bonds = sorted(
        bond_parities,
        key=lambda b: (sorted(priorities[k] for k in b), bond_parities[b]),
    )
    return atoms, bonds


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


def _flipped(table: _Table) -> _Table:
    """A bond's local stereo with every entry flipped."""
    return tuple(tuple(None if p is None else not p for p in row) for row in table)


def _replaced_none(table: _Table) -> _Table:
    """A bond's local stereo with an ill-defined entry read as `True`, for the test."""
    return tuple(tuple(True if p is None else p for p in row) for row in table)


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
