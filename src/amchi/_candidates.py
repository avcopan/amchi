"""The stereocenter candidates of a species or a TS, and the edits that correct them.

The rules are the old AutoMol's (`tetrahedral_atoms`, `rigid_planar_bonds`,
`stereocenter_candidates`), on the graph alone: the candidates never change with the
conformer (W3; decision D21). A TS applies them to each side, reactant and product (W4).
"""

from collections import deque
from collections.abc import Collection, Iterable, Sequence

from ._element import lone_pairs
from ._graph import MolGraph, explicit_graph, product_graph, reactant_graph
from ._resonance import rigid_planar_bonds, unsaturations
from ._types import Bond, Candidates, Connectivity

# Valence + lone pairs of a tetrahedral atom, and the largest ring that a stereo bond
# can't be in
_OCTET = 4
_MAX_RING_SIZE = 7


def candidates(symbols: Sequence[str], connectivity: Connectivity) -> Candidates:
    """The stereocenter candidates of a species or a TS, with their neighbors.

    An atom candidate has at least 3 neighbors and neighbors + lone pairs = 4. A bond
    candidate is rigid and planar, has a neighbor on each side, is in no ring of fewer
    than 8 atoms, and doesn't join two atom candidates. The neighbors are listed
    hydrogens first, then the other atoms, each by ascending index (D21).

    For a TS, the rules apply to each side. An atom's neighbors are those of the
    reactant side if it is a candidate there, else of the product side; an atom that is
    a candidate on both sides with different neighbors also gets the product-side ones
    in `reverse_atoms` (D26). A bond's neighbors are united over the sides where it is
    rigid and planar, and the ring test follows the bonds of the whole TS, leaving and
    reaching the bond through its listed neighbors only.

    Examples:
        >>> import amchi
        >>> # CHFClBr and FN=NF
        >>> connectivity = amchi.Connectivity(bonds=[(0, 1), (0, 2), (0, 3), (0, 4)])
        >>> amchi.candidates(["C", "F", "Cl", "Br", "H"], connectivity).atoms
        {0: (4, 1, 2, 3)}
        >>> connectivity = amchi.Connectivity(bonds=[(0, 1), (1, 2), (2, 3)])
        >>> amchi.candidates(["F", "N", "N", "F"], connectivity).bonds
        {(1, 2): ((0,), (3,))}
        >>> # An F atom moves from one C of CH3F to the other in the TS F...CH3...F
        >>> connectivity = amchi.Connectivity(
        ...     bonds=[(0, 1), (0, 2), (0, 3)], breaking=[(0, 4)], forming=[(0, 5)]
        ... )
        >>> result = amchi.candidates(["C", "H", "H", "H", "F", "F"], connectivity)
        >>> result.atoms, result.reverse_atoms
        ({0: (1, 2, 3, 4)}, {0: (1, 2, 3, 5)})
    """
    graph = explicit_graph(symbols, connectivity)
    is_ts = bool(connectivity.breaking or connectivity.forming)
    sides = [reactant_graph(graph), product_graph(graph)] if is_ts else [graph]

    def ordered(neighbors: Iterable[int]) -> tuple[int, ...]:
        return tuple(sorted(neighbors, key=lambda n: (symbols[n] != "H", n)))

    def tetrahedral(side: MolGraph) -> dict[int, tuple[int, ...]]:
        return {
            atom: ordered(side[atom])
            for atom in side
            if side.degree(atom) >= 3  # noqa: PLR2004
            and side.degree(atom) + lone_pairs(symbols[atom]) == _OCTET
        }

    on_side = [tetrahedral(side) for side in sides]
    atoms = on_side[-1] | on_side[0]  # the reactant side wins
    reverse_atoms = {
        atom: on_side[-1][atom]
        for atom, neighbors in on_side[0].items()
        if atom in on_side[-1] and set(on_side[-1][atom]) != set(neighbors)
    }

    # Each bond's neighbors, united over the sides where it is rigid and planar. The
    # unsaturation test is a cheap filter before the resonance search.
    neighbors: dict[Bond, tuple[set[int], set[int]]] = {}
    for side in sides:
        unsaturation = unsaturations(side)
        unsaturated = [
            (i, j)
            for i, j in (_bond_key(bond) for bond in side.edges)
            if unsaturation[i] > 0 and unsaturation[j] > 0
        ]
        for i, j in rigid_planar_bonds(side, unsaturated):
            side_i, side_j = neighbors.setdefault((i, j), (set(), set()))
            side_i |= set(side[i]) - {j}
            side_j |= set(side[j]) - {i}
    bonds = {
        (i, j): (ordered(side_i), ordered(side_j))
        for (i, j), (side_i, side_j) in sorted(neighbors.items())
        if side_i
        and side_j
        and not _in_small_ring(graph, i, j, side_i, side_j)
        and not (i in atoms and j in atoms)
    }
    return Candidates(atoms=atoms, reverse_atoms=reverse_atoms, bonds=bonds)


def add_bond_candidate(
    candidates: Candidates, connectivity: Connectivity, bond: Collection[int]
) -> Candidates:
    """Add a bond candidate, with all its neighbors by ascending index.

    No filter is applied: a caller overrides the perception on purpose. Remove an atom
    candidate of the bond yourself where it applies. For a TS, the bond and its
    neighbors may be breaking or forming bonds.

    Raises:
        ValueError: The pair isn't a bond, or one of its atoms has no other neighbor.

    Examples:
        >>> import amchi
        >>> connectivity = amchi.Connectivity(bonds=[(0, 1), (1, 2), (2, 3)])
        >>> amchi.add_bond_candidate(amchi.Candidates(), connectivity, (2, 1)).bonds
        {(1, 2): ((0,), (3,))}
    """
    i, j = key = _bond_key(bond)
    if key not in _all_bonds(connectivity):
        msg = f"Not a bond: {key}"
        raise ValueError(msg)
    side_i = tuple(sorted(n for n in _neighbors(connectivity, i) if n != j))
    side_j = tuple(sorted(n for n in _neighbors(connectivity, j) if n != i))
    if not (side_i and side_j):
        msg = f"A bond candidate needs a neighbor on each side: {key}"
        raise ValueError(msg)
    bonds = {**candidates.bonds, key: (side_i, side_j)}
    return candidates.model_copy(update={"bonds": dict(sorted(bonds.items()))})


def remove_candidate(candidates: Candidates, key: int | Collection[int]) -> Candidates:
    """Remove an atom candidate (an index) or a bond candidate (a pair).

    An atom candidate is removed with its reverse neighbors, if it has them.

    Raises:
        KeyError: It isn't a candidate.

    Examples:
        >>> import amchi
        >>> candidates = amchi.Candidates(
        ...     atoms={0: (1, 2, 3)}, bonds={(4, 5): ([6], [7])}
        ... )
        >>> amchi.remove_candidate(candidates, 0)
        Candidates(atoms={}, reverse_atoms={}, bonds={(4, 5): ((6,), (7,))})
        >>> amchi.remove_candidate(candidates, (5, 4))
        Candidates(atoms={0: (1, 2, 3)}, reverse_atoms={}, bonds={})
    """
    atoms, bonds = dict(candidates.atoms), dict(candidates.bonds)
    reverse_atoms = dict(candidates.reverse_atoms)
    if isinstance(key, int):
        if key not in atoms:
            msg = f"Not an atom candidate: {key}"
            raise KeyError(msg)
        del atoms[key]
        reverse_atoms.pop(key, None)
    else:
        bond = _bond_key(key)
        if bond not in bonds:
            msg = f"Not a bond candidate: {bond}"
            raise KeyError(msg)
        del bonds[bond]
    return Candidates(atoms=atoms, reverse_atoms=reverse_atoms, bonds=bonds)


def _bond_key(bond: Collection[int]) -> Bond:
    i, j = sorted(bond)
    return (i, j)


def _all_bonds(connectivity: Connectivity) -> tuple[Bond, ...]:
    return (*connectivity.bonds, *connectivity.breaking, *connectivity.forming)


def _neighbors(connectivity: Connectivity, atom: int) -> list[int]:
    return [
        j if i == atom else i for i, j in _all_bonds(connectivity) if atom in (i, j)
    ]


def _in_small_ring(
    graph: MolGraph, i: int, j: int, side_i: Collection[int], side_j: Collection[int]
) -> bool:
    """Is the bond in a ring of at most 7 atoms through its listed neighbors?

    That is, a path of 6 or fewer bonds from `i` to `j` that avoids the bond itself,
    leaves `i` through `side_i` and reaches `j` from `side_j`. The graph has all the
    bonds of a TS, reacting ones included.
    """
    distance = {i: 0}
    queue = deque([i])
    while queue:
        atom = queue.popleft()
        if distance[atom] >= _MAX_RING_SIZE - 1:
            continue
        for neighbor in graph[atom]:
            if atom == i and (neighbor == j or neighbor not in side_i):
                continue
            if neighbor == j:
                if atom in side_j:
                    return True
                continue
            if neighbor not in distance:
                distance[neighbor] = distance[atom] + 1
                queue.append(neighbor)
    return False
