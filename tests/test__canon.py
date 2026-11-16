"""Tests for the canonicalization, against the old code's dumped intermediates."""

from dataclasses import dataclass

import numpy as np
import pytest
from numpy.typing import ArrayLike

from amchi import Candidates, Connectivity, IllDefinedParityError
from amchi._canon import CanonicalForm, _first_pass_parities, canonical_forms
from amchi._graph import MolGraph, explicit_graph, implicit_graph
from amchi._parity import stereo_graph


@dataclass
class Geometry:
    """A minimal GeometryLike."""

    symbols: list[str]
    coordinates: ArrayLike


def _old_candidates(old_candidates: dict) -> Candidates:
    return Candidates(
        atoms=dict(old_candidates["atoms"]),
        bonds={tuple(bond): tuple(sides) for bond, sides in old_candidates["bonds"]},
    )


def _graph(
    symbols: list[str], xyz: ArrayLike, bonds: list, candidates: Candidates
) -> MolGraph:
    graph = implicit_graph(explicit_graph(symbols, Connectivity(bonds=bonds)))
    return stereo_graph(graph, Geometry(symbols, xyz), candidates)


def _species(records: list[dict]) -> list[dict]:
    return [r for r in records if r["kind"] != "ts"]


def _record_graph(record: dict, candidates: Candidates | None = None) -> MolGraph:
    old = record["old"]
    candidates = candidates or _old_candidates(old["candidates"])
    return _graph(record["symbols"], record["xyz"], old["bonds"], candidates)


def _mirror(graph: MolGraph) -> MolGraph:
    """The same graph with every local atom parity inverted."""
    for atom, stereo in graph.nodes(data="stereo"):
        if stereo is not None:
            parity = None if stereo.parity is None else not stereo.parity
            graph.nodes[atom]["stereo"] = stereo.model_copy(update={"parity": parity})
    return graph


def test__canonical_forms(records: list[dict]) -> None:
    """The priorities, numbers and enantiomer of the 600 species are the old ones."""
    species = _species(records)
    assert len(species) == 600
    for record in species:
        old = record["old"]
        (form,) = canonical_forms(_record_graph(record))
        (ref,) = old["canonical"]
        priorities = {k: abs(p) for k, p in ref["priorities"]}
        assert form.priorities == priorities, record["id"]
        assert form.numbers == dict(ref["numbers"]), record["id"]
        assert form.is_canonical_enantiomer == ref["is_canonical_enantiomer"]
        assert form.is_canonical_direction is None


def test__first_pass_parities(records: list[dict]) -> None:
    """The canonical parities of the unreflected pass are the old ones."""
    for record in _species(records):
        old = record["old"]
        ((atoms, bonds),) = _first_pass_parities(_record_graph(record))
        assert atoms == dict(old["canonical_parities"]["atoms"]), record["id"]
        old_bonds = {tuple(b): p for b, p in old["canonical_parities"]["bonds"]}
        assert bonds == old_bonds, record["id"]


def test__canonical_forms__enantiomer_counts(records: list[dict]) -> None:
    """157 closed shells and 187 radicals have an m layer, so the test covers both."""
    counts = {"closed": 0, "radical": 0}
    for record in _species(records):
        (form,) = canonical_forms(_record_graph(record))
        counts[record["kind"]] += form.is_canonical_enantiomer is not None
    assert counts == {"closed": 157, "radical": 187}


def test__canonical_forms__mirror(records: list[dict]) -> None:
    """Inverting the local atom parities flips the enantiomer choice.

    For a chiral component the numbers and parities stay the same.
    """
    flip = {True: False, False: True, None: None}
    for record in _species(records):
        (form,) = canonical_forms(_record_graph(record))
        (mirror,) = canonical_forms(_mirror(_record_graph(record)))
        assert mirror.is_canonical_enantiomer == flip[form.is_canonical_enantiomer]
        if form.is_canonical_enantiomer is not None:
            assert _same(mirror, form), record["id"]


def test__canonical_forms__neighbor_order(records: list[dict]) -> None:
    """Listing the neighbors of each candidate in another order changes nothing."""
    rng = np.random.default_rng(0)
    for record in _species(records):
        candidates = _old_candidates(record["old"]["candidates"])
        shuffled = Candidates(
            atoms={a: _shuffled(ns, rng) for a, ns in candidates.atoms.items()},
            bonds={
                b: tuple(_shuffled(side, rng) for side in sides)
                for b, sides in candidates.bonds.items()
            },
        )
        (form,) = canonical_forms(_record_graph(record))
        (other,) = canonical_forms(_record_graph(record, shuffled))
        assert other == form, record["id"]


def _shuffled(items: tuple[int, ...], rng: np.random.Generator) -> tuple[int, ...]:
    return tuple(int(k) for k in rng.permutation(items))


def _same(form: CanonicalForm, other: CanonicalForm) -> bool:
    return (
        form.numbers == other.numbers
        and form.atom_parities == other.atom_parities
        and form.bond_parities == other.bond_parities
    )


# CHFClBr: C0 with H1, F2, Cl3, Br4 at the corners of a tetrahedron
CHFCLBR = ["C", "H", "F", "Cl", "Br"]
TETRAHEDRON = [[0, 0, 0], [1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]
PLANAR = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [-1, -1, 0]]
STAR = [(0, 1), (0, 2), (0, 3), (0, 4)]


def test__canonical_forms__chfclbr() -> None:
    """A single stereocenter: the kept enantiomer has parity False."""
    candidates = Candidates(atoms={0: (1, 2, 3, 4)})
    (form,) = canonical_forms(_graph(CHFCLBR, TETRAHEDRON, STAR, candidates))
    mirror_xyz = (np.array(TETRAHEDRON) * [-1, 1, 1]).tolist()
    (mirror,) = canonical_forms(_graph(CHFCLBR, mirror_xyz, STAR, candidates))
    assert form.atom_parities == mirror.atom_parities == {0: False}
    assert {form.is_canonical_enantiomer, mirror.is_canonical_enantiomer} == {
        True,
        False,
    }
    assert sorted(form.numbers.values()) == [0, 1, 2, 3]


def test__canonical_forms__ill_defined() -> None:
    """An ill-defined parity raises only when its candidate is a stereocenter."""
    candidates = Candidates(atoms={0: (1, 2, 3, 4)})
    with pytest.raises(IllDefinedParityError, match="atom 0"):
        canonical_forms(_graph(CHFCLBR, PLANAR, STAR, candidates))
    # CH2FCl: two implicit hydrogens tie, so C is not a stereocenter
    (form,) = canonical_forms(
        _graph(["C", "H", "H", "F", "Cl"], PLANAR, STAR, candidates)
    )
    assert form.atom_parities == {}
    assert form.is_canonical_enantiomer is None


def test__canonical_forms__components() -> None:
    """One form per component, by the lowest atom; no candidates, one refinement."""
    graph = _graph(
        ["O", "H", "H", "C", "H", "H", "H", "H"],
        np.zeros((8, 3)),
        [(0, 1), (0, 2), (3, 4), (3, 5), (3, 6), (3, 7)],
        Candidates(),
    )
    forms = canonical_forms(graph)
    assert [f.numbers for f in forms] == [{0: 0}, {3: 0}]
