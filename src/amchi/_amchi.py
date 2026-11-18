"""The public pipeline: the AMChI string of a species or a TS from its three inputs.

The components are written one by one and joined as the old AutoMol joins them
(`amchi.base.sorted_join`): sorted by the old sort key, highest first, with repeated
layers written once behind a multiplier (`2*`).
"""

import itertools
from collections import Counter
from collections.abc import Sequence

import numpy as np

from ._candidates import candidates as _candidates
from ._canon import CanonicalForm, canonical_forms
from ._connect import connectivity as _connectivity
from ._connect import ts_connectivity as _ts_connectivity
from ._element import valence
from ._graph import MolGraph, explicit_graph, implicit_graph
from ._layers import component_layers, component_string
from ._parity import stereo_graph
from ._types import Candidates, Connectivity, GeometryLike

_LETTERS = ("c", "h", "b", "t", "m", "s", "k", "f", "r")
"""The layers after the formula, in the order they are written."""


def string(
    geometry: GeometryLike,
    connectivity: Connectivity,
    candidates: Candidates,
    *,
    stereo: bool = True,
) -> str:
    """The AMChI string of a species or a transition state (TS).

    A TS has breaking or forming bonds in its connectivity, and its string has the
    `k` (breaking) and `f` (forming) layers, as read in its canonical direction, and
    an `r` layer if the given direction is the reverse of that (`r1`) or not (`r0`).
    With `stereo=False`, the candidates are ignored, and the string has no `b`, `t`,
    `m` or `s` layer.

    Raises:
        UnsupportedElementError: A symbol has no element data.
        IllDefinedParityError: A stereocenter's parity is exactly zero.
        ValueError: The inputs don't fit each other: the coordinates aren't N x 3,
            a bond names a missing atom, or a candidate's neighbors aren't the
            neighbors of its atom or bond (for a TS: aren't among them, over the
            breaking and forming bonds too).

    Examples:
        >>> import amchi
        >>> from types import SimpleNamespace
        >>> water = SimpleNamespace(
        ...     symbols=["O", "H", "H"],
        ...     coordinates=[[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]],
        ... )
        >>> connectivity = amchi.Connectivity(bonds=[(0, 1), (0, 2)])
        >>> amchi.string(water, connectivity, amchi.Candidates())
        'AMChI=2S/H2O/h1H2'
    """
    symbols = list(geometry.symbols)
    if not stereo:
        candidates = Candidates()
    _check(symbols, geometry, connectivity, candidates)
    graph = implicit_graph(explicit_graph(symbols, connectivity))
    graph = stereo_graph(graph, geometry, candidates)
    forms = canonical_forms(graph)
    return _join(graph, forms)


def from_geometry(geometry: GeometryLike, *, stereo: bool = True) -> str:
    """The AMChI string of a species, with its connectivity and candidates derived.

    Raises the errors of `string`.

    Examples:
        >>> import amchi
        >>> from types import SimpleNamespace
        >>> # CHFClBr and its mirror image
        >>> chfclbr = SimpleNamespace(
        ...     symbols=["C", "H", "F", "Cl", "Br"],
        ...     coordinates=[
        ...         [0, 0, 0], [0, 0, 1.09], [1.31, 0, -0.46],
        ...         [-0.83, 1.44, -0.59], [-0.91, -1.58, -0.65],
        ...     ],
        ... )
        >>> amchi.from_geometry(chfclbr)
        'AMChI=2S/CHBrClF/c2-1(3)4/h1H/t1-/m1/s1'
        >>> chfclbr.coordinates = [[-x, y, z] for x, y, z in chfclbr.coordinates]
        >>> amchi.from_geometry(chfclbr)
        'AMChI=2S/CHBrClF/c2-1(3)4/h1H/t1-/m0/s1'
        >>> amchi.from_geometry(chfclbr, stereo=False)
        'AMChI=2S/CHBrClF/c2-1(3)4/h1H'
    """
    connectivity = _connectivity(geometry)
    if not stereo:
        return string(geometry, connectivity, Candidates(), stereo=False)
    candidates = _candidates(geometry.symbols, connectivity)
    return string(geometry, connectivity, candidates)


def from_ts_geometries(
    geometry: GeometryLike,
    reactant_geometry: GeometryLike,
    product_geometry: GeometryLike,
    *,
    stereo: bool = True,
) -> str:
    """The AMChI string of a TS, with its connectivity and candidates derived.

    The three geometries are in the same atom order: the TS (for the parities) and
    its reactant side and product side (for the bonds and the candidates). With
    `stereo=False` the string has no `b`, `t`, `m` or `s` layer.

    Raises the errors of `string`, and `ValueError` if the atoms differ.

    Examples:
        >>> import amchi
        >>> from types import SimpleNamespace
        >>> symbols = ["C", "H", "H", "H", "H", "O", "H"]
        >>> ring = [[0, 0, 0.3], [1, 0, 0.3], [-0.5, 0.87, 0.3], [-0.5, -0.87, 0.3]]
        >>> # CH3 + H2O: H4 moves from O5 to C0
        >>> ts = SimpleNamespace(
        ...     symbols=symbols,
        ...     coordinates=[[0, 0, 0], *ring[1:], [0, 0, -1.3], [0, 0, -2.5],
        ...                  [0.9, 0, -2.8]],
        ... )
        >>> reactant = SimpleNamespace(
        ...     symbols=symbols,
        ...     coordinates=[[0, 0, 0], *ring[1:], [-0.5, 0.5, -3.2], [0, 0, -2.5],
        ...                  [0.9, 0, -2.8]],
        ... )
        >>> product = SimpleNamespace(
        ...     symbols=symbols,
        ...     coordinates=[[0, 0, 0], *ring[1:], [0, 0, -1.09], [0, 0, -3.2],
        ...                  [0.9, 0, -3.5]],
        ... )
        >>> amchi.from_ts_geometries(ts, reactant, product)
        'AMChI=2S/CH5O/c1-3-2/h2H,1H3/k3-2/f3-1/r0'
        >>> amchi.from_ts_geometries(ts, product, reactant)
        'AMChI=2S/CH5O/c1-3-2/h2H,1H3/k3-2/f3-1/r1'
    """
    connectivity = _ts_connectivity(reactant_geometry, product_geometry)
    if not stereo:
        return string(geometry, connectivity, Candidates(), stereo=False)
    candidates = _candidates(geometry.symbols, connectivity)
    return string(geometry, connectivity, candidates)


def _check(
    symbols: Sequence[str],
    geometry: GeometryLike,
    connectivity: Connectivity,
    candidates: Candidates,
) -> None:
    """Check the three inputs against each other."""
    for symbol in set(symbols):
        valence(symbol)  # raises UnsupportedElementError
    shape = np.shape(geometry.coordinates)
    if shape != (len(symbols), 3):
        msg = f"Expected {len(symbols)} x 3 coordinates, got {shape}"
        raise ValueError(msg)
    # A species' candidate lists exactly its atom's or bond's neighbors. A TS's
    # lists those of one side or of several, so only some of all its neighbors.
    is_ts = bool(connectivity.breaking or connectivity.forming)
    fits = set.issubset if is_ts else set.__eq__
    neighbors: list[set[int]] = [set() for _ in symbols]
    for i, j in (*connectivity.bonds, *connectivity.breaking, *connectivity.forming):
        if j >= len(symbols):
            msg = f"A bond names a missing atom: {(i, j)}"
            raise ValueError(msg)
        neighbors[i].add(j)
        neighbors[j].add(i)
    for atoms in (candidates.atoms, candidates.reverse_atoms):
        for atom, atom_neighbors in atoms.items():
            if atom >= len(symbols) or not fits(set(atom_neighbors), neighbors[atom]):
                msg = f"An atom candidate's neighbors aren't its atom's: {atom}"
                raise ValueError(msg)
    for (i, j), (side_i, side_j) in candidates.bonds.items():
        if (
            j >= len(symbols)
            or i not in neighbors[j]
            or not fits(set(side_i), neighbors[i] - {j})
            or not fits(set(side_j), neighbors[j] - {i})
        ):
            msg = f"A bond candidate's neighbors aren't its bond's: {(i, j)}"
            raise ValueError(msg)


def _join(graph: MolGraph, forms: Sequence[CanonicalForm]) -> str:
    """The string of all components, sorted and joined as the old code does."""
    layers = [component_layers(graph, form) for form in forms]
    symbols = _sorted_symbols({graph.nodes[k]["symbol"] for k in graph} | {"H"})
    keys = [_sort_key(graph, form, symbols) for form in forms]
    order = sorted(range(len(forms)), key=keys.__getitem__, reverse=True)
    layers = [layers[i] for i in order]

    joined = {"formula": _join_layer([x["formula"] for x in layers], "", ".")}
    for letter in _LETTERS:
        values = [x.get(letter, "") for x in layers]
        if not any(values):
            continue
        if letter == "m":
            joined[letter] = "".join(value or "." for value in values)
        elif letter == "s":
            joined[letter] = max((value for value in values if value), key=int)
        else:
            joined[letter] = _join_layer(values, "*", ";")
    return component_string(joined)


def _join_layer(values: Sequence[str], count_separator: str, separator: str) -> str:
    """Join one layer of each component, with a count before a repeated value.

    As the old `join_layer_strings`: equal neighbors are counted (`2*1H2`), and an
    empty value is kept as an empty slot between separators.
    """
    parts = []
    for value, group in itertools.groupby(values):
        count = len(list(group))
        if value:
            parts.append(f"{count}{count_separator}{value}" if count > 1 else value)
        else:
            parts.append(separator * (count - 1))
    return separator.join(parts)


def _sorted_symbols(symbols: set[str]) -> list[str]:
    """C first, then alphabetical: the old `argsort` doesn't put H second."""
    return sorted(symbols, key=lambda s: (s != "C", s))


def _sort_key(graph: MolGraph, form: CanonicalForm, symbols: Sequence[str]) -> tuple:
    """The old sort key of a component (`amchi.base.argsort`), from its layers.

    The formula counts, the neighbors and the hydrogen count of each atom by
    canonical number, the bond and atom parities in the order of their layers, and
    whether it is the inverted enantiomer. The charge, always 0 here, is left out.
    Before the last, unlike the old key: the breaking and forming bonds as read in the
    canonical direction, so that two components equal but for them are ordered by
    them, and not by their `m` (which reflecting flips; `rgd1:MR_242866_0`).
    """
    atoms = sorted(form.numbers, key=form.numbers.__getitem__)
    counts = Counter(graph.nodes[k]["symbol"] for k in atoms)
    counts["H"] += sum(graph.nodes[k]["hydrogens"] for k in atoms)
    formula = tuple(counts[s] for s in symbols)
    neighbors = [sorted(form.numbers[n] for n in graph[k]) for k in atoms]
    hydrogens = [graph.nodes[k]["hydrogens"] for k in atoms]
    bond_parities = [
        parity
        for _, parity in sorted(
            (sorted((form.numbers[i], form.numbers[j]), reverse=True), parity)
            for (i, j), parity in form.bond_parities.items()
        )
    ]
    atom_parities = [form.atom_parities[k] for k in atoms if k in form.atom_parities]
    changes = ("breaking", "forming")
    if form.is_canonical_direction is False:
        changes = changes[::-1]
    reacting = [
        sorted(
            sorted((form.numbers[i], form.numbers[j]), reverse=True)
            for i, j, change in graph.edges(atoms, data="change")
            if change == value
        )
        for value in changes
    ]
    # None (no m layer) sorts below both values, where the old code would fail
    inverted = (
        -1 if form.is_canonical_enantiomer is None else not form.is_canonical_enantiomer
    )
    return (
        formula,
        neighbors,
        hydrogens,
        bond_parities,
        atom_parities,
        reacting,
        inverted,
    )
