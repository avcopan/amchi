"""Tests for the canonicalization, against the old code's dumped intermediates."""

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest
from numpy.typing import ArrayLike

from amchi import Candidates, Connectivity, IllDefinedParityError, candidates
from amchi._canon import CanonicalForm, _first_pass_parities, canonical_forms
from amchi._graph import MolGraph, explicit_graph, implicit_graph
from amchi._parity import stereo_graph

DATA = Path(__file__).parent / "data"


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


def test__canonical_forms(records: list[dict], ring_ids: set[str]) -> None:
    """The priorities, numbers and enantiomer of the 600 species are the old ones.

    The 12 with a para-stereocenter (D32, D33) are left out: they refine further.
    """
    species = _species(records)
    assert len(species) == 600
    for record in species:
        if record["id"] in ring_ids:
            continue
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


def test__canonical_forms__bond_table() -> None:
    """A bond's parity is the table entry of its highest-priority neighbors (D27).

    HClC=CFH, with Cl perpendicular to F but H not perpendicular to H. The stored
    order puts the hydrogens last, so only the Cl-F entry is ill-defined.
    """
    symbols = ["C", "C", "Cl", "H", "F", "H"]
    bonds = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)]
    xyz = [[0, 0, 0], [0, 0, 1.3], [1, 0, -0.7], [-1, 0.5, -0.7], [0, 1, 2], [0, -1, 2]]
    candidates = Candidates(bonds={(0, 1): ((2, 3), (4, 5))})
    graph = _graph(symbols, xyz, bonds, candidates)
    assert graph.edges[0, 1]["stereo"].parities == ((None, None), (False, True))
    with pytest.raises(IllDefinedParityError, match=r"bond \(0, 1\)"):
        canonical_forms(graph)

    # Turned so that Cl and F are trans: the Cl-F entry is the parity
    xyz[4], xyz[5] = [-1, 0, 2], [1, 0.3, 2]
    graph = _graph(symbols, xyz, bonds, candidates)
    assert graph.edges[0, 1]["stereo"].parities[0][0] is True
    (form,) = canonical_forms(graph)
    assert form.bond_parities == {(0, 1): True}


# Transition states (W4) ---------------------------------------------------------------
# The TSs whose mirror image is their own reverse (D29): reflecting one reverses it
SELF_REVERSE_MIRROR = {
    "rgd1:MR_35814_1",
    "rgd1:MR_448099_1",
    "rgd1:MR_103507_0",
    "rgd1:MR_323646_1",
}
FLIP = {True: False, False: True, None: None}


def _ts_forms(
    record: dict, *, reverse: bool = False, mirror: bool = False
) -> list[CanonicalForm]:
    """The canonical forms of a TS record, read reversed or reflected."""
    old = record["old"]
    breaking, forming = old["breaking"], old["forming"]
    if reverse:
        breaking, forming = forming, breaking
    connectivity = Connectivity(bonds=old["bonds"], breaking=breaking, forming=forming)
    xyz = np.array(record["xyz"]) * ([-1, 1, 1] if mirror else 1)
    graph = implicit_graph(explicit_graph(record["symbols"], connectivity))
    geometry = Geometry(record["symbols"], xyz)
    return canonical_forms(
        stereo_graph(graph, geometry, candidates(record["symbols"], connectivity))
    )


def _ts(records: list[dict]) -> list[dict]:
    return [r for r in records if r["kind"] == "ts"]


def test__canonical_forms__ts(records: list[dict]) -> None:
    """The 289 TSs whose string is the old one give the old canonical forms.

    The priorities with their signs, the numbers and both flags, component by
    component. The other 11 are the intended differences of D28 and D29.
    """
    differences = {
        d["id"] for d in json.loads((DATA / "ts_differences.json").read_text())
    }
    checked = [r for r in _ts(records) if r["id"] not in differences]
    assert len(checked) == 289
    for record in checked:
        dumped = {
            frozenset(k for k, _ in c["priorities"]): c
            for c in record["old"]["canonical"]
        }
        for form in _ts_forms(record):
            ref = dumped[frozenset(form.priorities)]
            assert form.priorities == dict(ref["priorities"]), record["id"]
            assert form.numbers == dict(ref["numbers"]), record["id"]
            assert form.is_canonical_direction == ref["is_canonical_direction"]
            assert form.is_canonical_enantiomer == ref["is_canonical_enantiomer"]


def test__canonical_forms__ts_reversed(records: list[dict]) -> None:
    """Swapping breaking and forming flips the direction flag only (D29).

    Where the two directions are equal (no flag), each read keeps its own given
    direction, so the numbers may differ by a symmetry of the TS.
    """
    for record in _ts(records):
        for form, other in zip(
            _ts_forms(record), _ts_forms(record, reverse=True), strict=True
        ):
            assert other.is_canonical_direction == FLIP[form.is_canonical_direction]
            assert other.is_canonical_enantiomer == form.is_canonical_enantiomer
            assert _numbered(other) == _numbered(form), record["id"]
            if form.is_canonical_direction is not None:
                assert _same(other, form), record["id"]


def test__canonical_forms__ts_mirrored(records: list[dict]) -> None:
    """Reflecting flips the enantiomer flag only (D29).

    Except for the 4 TSs whose mirror image is their reverse: they have no enantiomer
    flag, and reflecting flips the direction flag instead.
    """
    for record in _ts(records):
        exception = record["id"] in SELF_REVERSE_MIRROR
        for form, other in zip(
            _ts_forms(record), _ts_forms(record, mirror=True), strict=True
        ):
            if exception and form.is_canonical_direction is not None:
                # The mirror image is read as the reversed TS
                assert form.is_canonical_enantiomer is None
                assert other.is_canonical_enantiomer is None
                assert other.is_canonical_direction is not form.is_canonical_direction
                assert _numbered(other) == _numbered(form), record["id"]
            else:
                assert _same(other, form), record["id"]
                assert (
                    other.is_canonical_enantiomer == FLIP[form.is_canonical_enantiomer]
                )
                assert other.is_canonical_direction == form.is_canonical_direction


def _numbered(form: CanonicalForm) -> tuple[list, list]:
    """The parities by canonical number, as the string gives them."""
    numbers = form.numbers
    return (
        sorted((numbers[a], p) for a, p in form.atom_parities.items()),
        sorted(
            (sorted(numbers[k] for k in b), p) for b, p in form.bond_parities.items()
        ),
    )
