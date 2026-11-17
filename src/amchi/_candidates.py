"""The stereocenter candidates of a species, and the edits that correct them.

The rules are the old AutoMol's (`tetrahedral_atoms`, `stereocenter_candidates`), on
the graph alone: the candidates never change with the conformer (W3; decision D21).
"""

from collections import deque
from collections.abc import Collection, Iterable, Sequence

from ._element import lone_pairs
from ._graph import MolGraph, explicit_graph
from ._resonance import rigid_planar_bonds, unsaturations
from ._types import Bond, Candidates, Connectivity

# Valence + lone pairs of a tetrahedral atom, and the largest ring that a stereo bond
# can't be in
_OCTET = 4
_MAX_RING_SIZE = 7


def candidates(symbols: Sequence[str], connectivity: Connectivity) -> Candidates:
    """The stereocenter candidates of a species, with their neighbors.

    An atom candidate has at least 3 neighbors and neighbors + lone pairs = 4. A bond
    candidate is rigid and planar, has a neighbor on each side, is in no ring of fewer
    than 8 atoms, and doesn't join two atom candidates. The neighbors are listed
    hydrogens first, then the other atoms, each by ascending index (D21).

    Examples:
        >>> import amchi
        >>> # CHFClBr and FN=NF
        >>> connectivity = amchi.Connectivity(bonds=[(0, 1), (0, 2), (0, 3), (0, 4)])
        >>> amchi.candidates(["C", "F", "Cl", "Br", "H"], connectivity).atoms
        {0: (4, 1, 2, 3)}
        >>> connectivity = amchi.Connectivity(bonds=[(0, 1), (1, 2), (2, 3)])
        >>> amchi.candidates(["F", "N", "N", "F"], connectivity).bonds
        {(1, 2): ((0,), (3,))}
    """
    if connectivity.breaking or connectivity.forming:
        msg = "The candidates of a transition state are not supported yet"
        raise NotImplementedError(msg)
    graph = explicit_graph(symbols, connectivity)

    def ordered(neighbors: Iterable[int]) -> tuple[int, ...]:
        return tuple(sorted(neighbors, key=lambda n: (symbols[n] != "H", n)))

    atoms = {
        atom: ordered(graph[atom])
        for atom in graph
        if graph.degree(atom) >= 3  # noqa: PLR2004
        and graph.degree(atom) + lone_pairs(symbols[atom]) == _OCTET
    }

    # The cheap filters first, then the resonance search on the bonds that are left
    unsaturation = unsaturations(graph)
    sides = {}
    for i, j in connectivity.bonds:
        side_i = [n for n in graph[i] if n != j]
        side_j = [n for n in graph[j] if n != i]
        if (
            unsaturation[i] > 0
            and unsaturation[j] > 0
            and side_i
            and side_j
            and not _in_small_ring(graph, i, j)
        ):
            sides[(i, j)] = (ordered(side_i), ordered(side_j))
    bonds = {
        bond: sides[bond]
        for bond in sorted(rigid_planar_bonds(graph, sides))
        if not (bond[0] in atoms and bond[1] in atoms)
    }
    return Candidates(atoms=atoms, bonds=bonds)


def add_bond_candidate(
    candidates: Candidates, connectivity: Connectivity, bond: Collection[int]
) -> Candidates:
    """Add a bond candidate, with all its neighbors by ascending index.

    No filter is applied: a caller overrides the perception on purpose. Remove an atom
    candidate of the bond yourself where it applies.

    Raises:
        ValueError: The pair isn't a bond, or one of its atoms has no other neighbor.

    Examples:
        >>> import amchi
        >>> connectivity = amchi.Connectivity(bonds=[(0, 1), (1, 2), (2, 3)])
        >>> amchi.add_bond_candidate(amchi.Candidates(), connectivity, (2, 1)).bonds
        {(1, 2): ((0,), (3,))}
    """
    if connectivity.breaking or connectivity.forming:
        msg = "The candidates of a transition state are not supported yet"
        raise NotImplementedError(msg)
    i, j = key = _bond_key(bond)
    if key not in connectivity.bonds:
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


def _neighbors(connectivity: Connectivity, atom: int) -> list[int]:
    return [j if i == atom else i for i, j in connectivity.bonds if atom in (i, j)]


def _in_small_ring(graph: MolGraph, i: int, j: int) -> bool:
    """Is the bond in a ring of at most 7 atoms?

    That is, a path of 6 or fewer bonds from `i` to `j` that avoids the bond itself.
    """
    distance = {i: 0}
    queue = deque([i])
    while queue:
        atom = queue.popleft()
        if distance[atom] >= _MAX_RING_SIZE - 1:
            continue
        for neighbor in graph[atom]:
            if atom == i and neighbor == j:
                continue
            if neighbor == j:
                return True
            if neighbor not in distance:
                distance[neighbor] = distance[atom] + 1
                queue.append(neighbor)
    return False
