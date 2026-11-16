"""The input models and the GeometryLike protocol."""

from collections.abc import Sequence
from typing import Any, Protocol, Self

from numpy.typing import ArrayLike
from pydantic import BaseModel, Field, field_validator, model_validator

from ._element import AmchiError

type Bond = tuple[int, int]
"""A bond, as a pair of atom indices `(i, j)` with `i < j`."""


class IllDefinedParityError(AmchiError):
    """A stereocenter whose parity the coordinates don't define."""


class GeometryLike(Protocol):
    """Anything with atomic symbols and Cartesian coordinates (N x 3, Å)."""

    @property
    def symbols(self) -> Sequence[str]:
        """The atomic symbols."""
        ...

    @property
    def coordinates(self) -> ArrayLike:
        """The coordinates, N x 3 in Å."""
        ...


def _sorted_pair(pair: Any) -> Any:  # noqa: ANN401
    try:
        return tuple(sorted(pair))
    except TypeError:
        return pair  # left for pydantic to reject


class Connectivity(BaseModel, frozen=True):
    """The bonds of a species or a transition state.

    `bonds` are the bonds present throughout. A TS also has `breaking` and `forming`
    bonds, which are not repeated in `bonds`: the reactant has `bonds + breaking` and
    the product `bonds + forming`. Each pair is sorted to `(i, j)` with `i < j`.

    Examples:
        >>> import amchi
        >>> amchi.Connectivity(bonds=[[1, 0], [0, 2]]).bonds
        ((0, 1), (0, 2))
    """

    bonds: tuple[Bond, ...]
    breaking: tuple[Bond, ...] = ()
    forming: tuple[Bond, ...] = ()

    @field_validator("bonds", "breaking", "forming", mode="before")
    @classmethod
    def _sort(cls, pairs: Any) -> Any:  # noqa: ANN401
        try:
            return tuple(sorted(_sorted_pair(pair) for pair in pairs))
        except TypeError:
            return pairs  # left for pydantic to reject

    @model_validator(mode="after")
    def _check(self) -> Self:
        pairs = [*self.bonds, *self.breaking, *self.forming]
        for i, j in pairs:
            if i < 0 or i == j:
                msg = f"Invalid bond: {(i, j)}"
                raise ValueError(msg)
        if len(set(pairs)) != len(pairs):
            repeated = sorted({p for p in pairs if pairs.count(p) > 1})
            msg = f"Bonds listed more than once: {repeated}"
            raise ValueError(msg)
        return self


class Candidates(BaseModel, frozen=True):
    """The stereocenter candidates, each with its stereo-determining neighbors.

    `atoms` maps an atom to its 3 or 4 neighbors. `bonds` maps a bond `(i, j)` to
    the neighbors of `i`, then of `j`, without the other bond atom (1 to 3 on each
    side). Any neighbor order is valid: a local parity is measured against the order
    given. A bond given as `(j, i)` is sorted, with its two sides swapped.

    Examples:
        >>> import amchi
        >>> candidates = amchi.Candidates(bonds={(3, 1): ([4], [0, 5])})
        >>> candidates.bonds
        {(1, 3): ((0, 5), (4,))}
    """

    atoms: dict[int, tuple[int, ...]] = Field(default_factory=dict)
    bonds: dict[Bond, tuple[tuple[int, ...], tuple[int, ...]]] = Field(
        default_factory=dict
    )

    @field_validator("bonds", mode="before")
    @classmethod
    def _sort_bonds(cls, bonds: Any) -> Any:  # noqa: ANN401
        if not isinstance(bonds, dict):
            return bonds
        sorted_bonds = {}
        for bond, sides in bonds.items():
            if isinstance(bond, Sequence) and len(bond) == 2 and bond[0] > bond[1]:  # noqa: PLR2004
                sorted_bonds[(bond[1], bond[0])] = tuple(reversed(sides))
            else:
                sorted_bonds[bond] = sides
        return sorted_bonds

    @model_validator(mode="after")
    def _check(self) -> Self:
        for atom, neighbors in self.atoms.items():
            if len(neighbors) not in (3, 4) or not _distinct(atom, *neighbors):
                msg = f"An atom candidate needs 3 or 4 distinct neighbors: {atom}"
                raise ValueError(msg)
        for bond, sides in self.bonds.items():
            ok = _distinct(*bond) and all(
                1 <= len(side) <= 3 and _distinct(*bond, *side)  # noqa: PLR2004
                for side in sides
            )
            if not ok:
                msg = f"A bond candidate needs 1 to 3 distinct neighbors a side: {bond}"
                raise ValueError(msg)
        return self


def _distinct(*atoms: int) -> bool:
    return len(set(atoms)) == len(atoms) and min(atoms) >= 0
