"""Local parities of the candidates, measured from the coordinates.

A local parity is defined against the neighbor order stored with its candidate. A
bond has one for each pair of neighbors across it, and a TS atom with reverse
neighbors has one for each neighbor list (D26, D27). A local parity is ill-defined
(`None`) only when its determinant or dot product is exactly zero, as in the old
AutoMol. How nearly planar a candidate is, is a diagnostic: `parity_margins`.
"""

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray
from pydantic import BaseModel, Field

from ._graph import AtomStereo, BondStereo, MolGraph
from ._types import Bond, Candidates, GeometryLike

type _Coordinates = NDArray[np.float64]


class ParityMargins(BaseModel, frozen=True):
    """How far each candidate is from an ill-defined parity, with the candidates' keys.

    For an atom, `|det| / d**3`, with `d` the mean distance to its neighbors: 3.08 for
    an ideal tetrahedron, 0.77 for an ideal pyramid with a lone pair, and 0 for a
    planar atom. `reverse_atoms` gives the same for a TS atom's reverse neighbors. For
    a bond, the smallest `|cos|` of the dihedral angle over the pairs of neighbors
    across it: 1 when planar, 0 at 90°.

    Examples:
        >>> import amchi
        >>> amchi.ParityMargins(atoms={0: 3.08}, bonds={}).atoms
        {0: 3.08}
    """

    atoms: dict[int, float]
    bonds: dict[Bond, float]
    reverse_atoms: dict[int, float] = Field(default_factory=dict)


def stereo_graph(
    graph: MolGraph, geometry: GeometryLike, candidates: Candidates
) -> MolGraph:
    """Write each candidate's neighbors and local parities on the graph, in place.

    Each candidate atom must be a node of the graph, and each candidate bond an edge.
    """
    coordinates = _coordinates(geometry)
    for atom, neighbors in candidates.atoms.items():
        reverse_neighbors = candidates.reverse_atoms.get(atom)
        graph.nodes[atom]["stereo"] = AtomStereo(
            neighbors=neighbors,
            parity=_atom_parity(coordinates, atom, neighbors),
            reverse_neighbors=reverse_neighbors,
            reverse_parity=None
            if reverse_neighbors is None
            else _atom_parity(coordinates, atom, reverse_neighbors),
        )
    for bond, sides in candidates.bonds.items():
        parities = tuple(
            tuple(_bond_parity(coordinates, bond, n0, n1) for n1 in sides[1])
            for n0 in sides[0]
        )
        graph.edges[bond]["stereo"] = BondStereo(neighbors=sides, parities=parities)
    return graph


def parity_margins(geometry: GeometryLike, candidates: Candidates) -> ParityMargins:
    """How nearly planar each candidate is, as a diagnostic.

    Examples:
        >>> import amchi
        >>> from types import SimpleNamespace
        >>> methane = SimpleNamespace(
        ...     symbols=["C", "H", "H", "H", "H"],
        ...     coordinates=[
        ...         [0, 0, 0], [1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]
        ...     ],
        ... )
        >>> candidates = amchi.Candidates(atoms={0: (1, 2, 3, 4)})
        >>> round(amchi.parity_margins(methane, candidates).atoms[0], 2)
        3.08
    """
    coordinates = _coordinates(geometry)
    bonds = {
        bond: min(
            _bond_margin(coordinates, bond, n0, n1)
            for n0 in sides[0]
            for n1 in sides[1]
        )
        for bond, sides in candidates.bonds.items()
    }
    return ParityMargins(
        atoms=_atom_margins(coordinates, candidates.atoms),
        bonds=bonds,
        reverse_atoms=_atom_margins(coordinates, candidates.reverse_atoms),
    )


def _coordinates(geometry: GeometryLike) -> _Coordinates:
    return np.asarray(geometry.coordinates, dtype=np.float64)


def _atom_margins(
    coordinates: _Coordinates, atoms: dict[int, tuple[int, ...]]
) -> dict[int, float]:
    margins = {}
    for atom, neighbors in atoms.items():
        det = _atom_determinant(coordinates, atom, neighbors)
        distances = np.linalg.norm(
            coordinates[list(neighbors)] - coordinates[atom], axis=1
        )
        margins[atom] = float(abs(det) / np.mean(distances) ** 3)
    return margins


def _atom_parity(
    coordinates: _Coordinates, atom: int, neighbors: Sequence[int]
) -> bool | None:
    det = _atom_determinant(coordinates, atom, neighbors)
    return None if det == 0.0 else bool(det > 0.0)


def _atom_determinant(
    coordinates: _Coordinates, atom: int, neighbors: Sequence[int]
) -> float:
    """The determinant of the rows `(1, x, y, z)` of the four points.

    With three neighbors, the atom itself is the first point.
    """
    points = [atom, *neighbors] if len(neighbors) == 3 else list(neighbors)  # noqa: PLR2004
    matrix = np.ones((4, 4))
    matrix[:, 1:] = coordinates[points]
    return float(np.linalg.det(matrix))


def _bond_parity(
    coordinates: _Coordinates, bond: Bond, n0: int, n1: int
) -> bool | None:
    """The local parity of a bond with `n0` and `n1` as reference neighbors: trans."""
    dot, _ = _bond_dot(coordinates, bond, ((n0,), (n1,)))
    return None if dot == 0.0 else bool(dot < 0.0)


def _bond_margin(coordinates: _Coordinates, bond: Bond, n0: int, n1: int) -> float:
    """The `|cos|` of the dihedral angle with `n0` and `n1` as reference neighbors."""
    dot, norms = _bond_dot(coordinates, bond, ((n0,), (n1,)))
    return abs(dot) / norms if norms else 0.0


def _bond_dot(
    coordinates: _Coordinates,
    bond: Bond,
    sides: tuple[Sequence[int], Sequence[int]],
) -> tuple[float, float]:
    """The dot product of the reference neighbors' vectors, normal to the bond axis.

    The reference neighbor of each end is the last one in the stored order. Also
    returns the product of the two vectors' lengths.
    """
    a, b = bond
    axis = coordinates[b] - coordinates[a]
    axis /= np.linalg.norm(axis)
    va = coordinates[sides[0][-1]] - coordinates[a]
    vb = coordinates[sides[1][-1]] - coordinates[b]
    va -= axis * (va @ axis)
    vb -= axis * (vb @ axis)
    return float(va @ vb), float(np.linalg.norm(va) * np.linalg.norm(vb))
