"""Connectivity from a geometry.

The rule is the old AutoMol's (`automol.geom.graph_without_stereo`), since the old
strings depend on it.
"""

from itertools import combinations

import numpy as np

from ._element import covalent_radius, valence, van_der_waals_radius
from ._types import Bond, Connectivity, GeometryLike

DISTANCE_FACTOR = 1.05


def connectivity(geometry: GeometryLike) -> Connectivity:
    """The bonds of a species, from its interatomic distances.

    Two atoms are bonded within 1.05 times the larger of the sum of their covalent
    radii and the mean of their van der Waals radii. Then an atom with more bonds
    than its valence loses its longest ones. The excess and the neighbors are counted
    before any bond is removed, as in the old code.

    Examples:
        >>> import amchi
        >>> from types import SimpleNamespace
        >>> water = SimpleNamespace(
        ...     symbols=["O", "H", "H"],
        ...     coordinates=[[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]],
        ... )
        >>> amchi.connectivity(water).bonds
        ((0, 1), (0, 2))
    """
    symbols = list(geometry.symbols)
    coordinates = np.asarray(geometry.coordinates, dtype=np.float64)
    if coordinates.shape != (len(symbols), 3):
        msg = f"Expected {len(symbols)} x 3 coordinates, got {coordinates.shape}"
        raise ValueError(msg)

    def distance(i: int, j: int) -> float:
        return float(np.linalg.norm(coordinates[i] - coordinates[j]))

    bonds: set[Bond] = {
        (i, j)
        for i, j in combinations(range(len(symbols)), 2)
        if distance(i, j) <= _distance_limit(symbols[i], symbols[j])
    }

    neighbors: dict[int, list[int]] = {i: [] for i in range(len(symbols))}
    for i, j in bonds:
        neighbors[i].append(j)
        neighbors[j].append(i)

    for atom, atom_neighbors in neighbors.items():
        excess = len(atom_neighbors) - valence(symbols[atom])
        if excess > 0:
            longest = sorted(sorted(atom_neighbors), key=lambda n: distance(atom, n))
            bonds -= {(min(atom, n), max(atom, n)) for n in longest[-excess:]}

    return Connectivity(bonds=bonds)


def _distance_limit(symbol1: str, symbol2: str) -> float:
    """The longest distance at which two atoms are bonded, in Å."""
    covalent = covalent_radius(symbol1) + covalent_radius(symbol2)
    van_der_waals = (van_der_waals_radius(symbol1) + van_der_waals_radius(symbol2)) / 2
    return max(covalent, van_der_waals) * DISTANCE_FACTOR
