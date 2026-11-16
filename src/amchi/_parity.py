"""Local parities of the candidates, measured from the coordinates.

A local parity is defined against the neighbor order stored with its candidate. It is
ill-defined (`None`) only when its determinant or dot product is exactly zero, as in
the old AutoMol. How nearly planar a candidate is, is a diagnostic: `parity_margins`.
"""

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray
from pydantic import BaseModel

from ._graph import AtomStereo, BondStereo, MolGraph
from ._types import Bond, Candidates, GeometryLike

type _Coordinates = NDArray[np.float64]


class ParityMargins(BaseModel, frozen=True):
    """How far each candidate is from an ill-defined parity, with the candidates' keys.

    For an atom, `|det| / d**3`, with `d` the mean distance to its neighbors: 3.08 for
    an ideal tetrahedron, 0.77 for an ideal pyramid with a lone pair, and 0 for a
    planar atom. For a bond, the `|cos|` of the dihedral angle between its reference
    neighbors: 1 when planar, 0 at 90°.
    """

    atoms: dict[int, float]
    bonds: dict[Bond, float]


def stereo_graph(
    graph: MolGraph, geometry: GeometryLike, candidates: Candidates
) -> MolGraph:
    """Write each candidate's neighbors and local parity on the graph, in place.

    Each candidate atom must be a node of the graph, and each candidate bond an edge.
    """
    coordinates = _coordinates(geometry)
    for atom, neighbors in candidates.atoms.items():
        det = _atom_determinant(coordinates, atom, neighbors)
        parity = None if det == 0.0 else bool(det > 0.0)
        graph.nodes[atom]["stereo"] = AtomStereo(neighbors=neighbors, parity=parity)
    for bond, sides in candidates.bonds.items():
        dot, _ = _bond_dot(coordinates, bond, sides)
        parity = None if dot == 0.0 else bool(dot < 0.0)
        graph.edges[bond]["stereo"] = BondStereo(neighbors=sides, parity=parity)
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
    atoms = {}
    for atom, neighbors in candidates.atoms.items():
        det = _atom_determinant(coordinates, atom, neighbors)
        distances = np.linalg.norm(
            coordinates[list(neighbors)] - coordinates[atom], axis=1
        )
        atoms[atom] = float(abs(det) / np.mean(distances) ** 3)
    bonds = {}
    for bond, sides in candidates.bonds.items():
        dot, norms = _bond_dot(coordinates, bond, sides)
        bonds[bond] = float(abs(dot) / norms) if norms else 0.0
    return ParityMargins(atoms=atoms, bonds=bonds)


def _coordinates(geometry: GeometryLike) -> _Coordinates:
    return np.asarray(geometry.coordinates, dtype=np.float64)


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
