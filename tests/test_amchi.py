"""Tests for the public pipeline, against the old strings.

A string equals the old one except for the version prefix and the c layer (D24). The
c layer must describe the same graph in the same canonical numbers as the old one;
the writer itself is tested against InChI in `test__layers.py`.
"""

import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest
from numpy.typing import ArrayLike

import amchi
from amchi import Candidates, Connectivity

from .test__layers import _inchi_graph

DATA = Path(__file__).parent / "data"


@dataclass
class Geometry:
    """A minimal GeometryLike."""

    symbols: list[str]
    coordinates: ArrayLike


def _c_graphs(layer: str) -> list[set[tuple[int, int]]]:
    """The bonds of each component of a c layer, with multipliers expanded."""
    graphs = []
    for part in layer.split(";"):
        count, _, value = part.rpartition("*")
        graph = _inchi_graph(value)
        graphs += [{tuple(sorted(e)) for e in graph.edges}] * int(count or 1)
    return graphs


def _assert_old(string: str, old: str, name: str) -> None:
    """The string equals the old one but for the prefix and the writing of c."""
    prefix, *layers = string.split("/")
    old_prefix, *old_layers = old.split("/")
    assert (prefix, old_prefix) == ("AMChI=2S", "AMChI=1"), name
    c_layers = [x[1:] for x in layers if x.startswith("c")]
    old_c_layers = [x[1:] for x in old_layers if x.startswith("c")]
    assert list(map(_c_graphs, c_layers)) == list(map(_c_graphs, old_c_layers)), name
    rest = [x for x in layers if not x.startswith("c")]
    old_rest = [x for x in old_layers if not x.startswith("c")]
    assert rest == old_rest, name


def _without_stereo(string: str) -> str:
    return "/".join(x for x in string.split("/") if x[0] not in "btms")


def _shuffled(geometry: Geometry, rng: np.random.Generator) -> Geometry:
    order = rng.permutation(len(geometry.symbols))
    symbols = [geometry.symbols[k] for k in order]
    return Geometry(symbols, np.asarray(geometry.coordinates)[order])


@pytest.fixture(scope="module")
def species(records: list[dict]) -> list[dict]:
    """The 600 species, each with its geometry and its string from `from_geometry`."""
    result = []
    for record in records:
        if record["kind"] != "ts":
            geometry = Geometry(record["symbols"], record["xyz"])
            string = amchi.from_geometry(geometry)
            result.append(record | {"geometry": geometry, "amchi": string})
    assert len(result) == 600
    return result


@pytest.fixture(scope="module")
def components() -> list[dict]:
    """Species of the sample put together, written by `baseline/dump_components.py`."""
    return json.loads((DATA / "components.json").read_text())


def test__from_geometry(species: list[dict]) -> None:
    """The 600 species give the old strings (D13), with the c layer of D24."""
    for record in species:
        _assert_old(record["amchi"], record["old"]["amchi"], record["id"])


def test__from_geometry__components(components: list[dict]) -> None:
    """Components are sorted, counted (`2*`) and joined as by the old code."""
    for case in components:
        string = amchi.from_geometry(Geometry(case["symbols"], case["xyz"]))
        _assert_old(string, case["amchi"], case["id"])
    by_id = {case["id"]: case["amchi"] for case in components}
    assert "/m10/" in by_id["enantiomers"]
    assert "2*1H2" in by_id["water+chiral+water"]


def test__from_geometry__stereo_false(species: list[dict]) -> None:
    """Without stereo, the string drops the b, t, m and s layers and nothing else.

    The old code does the same on all 600 species (`geom.amchi(geo, stereo=False)`).
    """
    for record in species:
        string = amchi.from_geometry(record["geometry"], stereo=False)
        assert string == _without_stereo(record["amchi"]), record["id"]


def test__from_geometry__shuffled(species: list[dict], components: list[dict]) -> None:
    """Shuffling the atoms leaves the string unchanged."""
    rng = np.random.default_rng(0)
    geometries = [(r["id"], r["geometry"], r["amchi"]) for r in species]
    geometries += [
        (c["id"], geometry, amchi.from_geometry(geometry))
        for c in components
        for geometry in [Geometry(c["symbols"], c["xyz"])]
    ]
    for name, geometry, string in geometries:
        for _ in range(2):
            shuffled = _shuffled(geometry, rng)
            assert amchi.from_geometry(shuffled) == string, name


def test__from_geometry__hand_cases(hand_cases: dict[str, list[dict]]) -> None:
    """The ring pairs, the larger molecules and HNNO give the old strings."""
    for group in ("hnno", "ring", "ring_control", "scaling"):
        for case in hand_cases[group]:
            geometry = Geometry(case["symbols"], case["xyz"])
            _assert_old(amchi.from_geometry(geometry), case["amchi"], case["id"])


def test__string__override(hand_cases: dict[str, list[dict]]) -> None:
    """The HNNO strings with corrected candidates give the old overridden strings."""
    for case in hand_cases["hnno"]:
        geometry = Geometry(case["symbols"], case["xyz"])
        connectivity = amchi.connectivity(geometry)
        candidates = amchi.candidates(case["symbols"], connectivity)
        for atom in case["override"]["remove_atoms"]:
            candidates = amchi.remove_candidate(candidates, atom)
        for bond in case["override"]["add_bonds"]:
            candidates = amchi.add_bond_candidate(candidates, connectivity, bond)
        string = amchi.string(geometry, connectivity, candidates)
        _assert_old(string, case["amchi_overridden"], case["id"])


# CHFClBr: C0 with H1, F2, Cl3, Br4 at the corners of a tetrahedron
CHFCLBR = ["C", "H", "F", "Cl", "Br"]
TETRAHEDRON = [[0, 0, 0], [1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]
STAR = Connectivity(bonds=[(0, 1), (0, 2), (0, 3), (0, 4)])
CENTER = Candidates(atoms={0: (1, 2, 3, 4)})


def test__string__mirror_and_shuffle() -> None:
    """Reflecting flips m and nothing else; any neighbor order of a candidate works."""
    string = amchi.string(Geometry(CHFCLBR, TETRAHEDRON), STAR, CENTER)
    mirror = (np.array(TETRAHEDRON) * [-1, 1, 1]).tolist()
    mirror_string = amchi.string(Geometry(CHFCLBR, mirror), STAR, CENTER)
    assert re.findall("/m(.)/", string + mirror_string) in (["0", "1"], ["1", "0"])
    assert re.sub("/m./", "", string) == re.sub("/m./", "", mirror_string)
    reordered = Candidates(atoms={0: (2, 1, 3, 4)})
    assert amchi.string(Geometry(CHFCLBR, mirror), STAR, reordered) == mirror_string
    stereo_free = amchi.string(Geometry(CHFCLBR, mirror), STAR, CENTER, stereo=False)
    assert stereo_free == _without_stereo(string)


CHFCLBR_GEOMETRY = Geometry(CHFCLBR, TETRAHEDRON)


@pytest.mark.parametrize(
    ("geometry", "connectivity", "candidates", "error"),
    [
        (Geometry(["C", "Xx"], np.eye(2, 3)), Connectivity(bonds=[(0, 1)]),
         Candidates(), (amchi.UnsupportedElementError, "Xx")),
        (Geometry(CHFCLBR, np.zeros((4, 3))), STAR, CENTER,
         (ValueError, "coordinates")),
        (CHFCLBR_GEOMETRY, Connectivity(bonds=[(0, 5)]), Candidates(),
         (ValueError, "missing atom")),
        (CHFCLBR_GEOMETRY, Connectivity(bonds=[(0, 1)], forming=[(0, 2)]),
         Candidates(), (NotImplementedError, "transition state")),
        (CHFCLBR_GEOMETRY, STAR, Candidates(atoms={0: (1, 2, 3)}),
         (ValueError, "atom candidate")),
        (CHFCLBR_GEOMETRY, STAR, Candidates(atoms={1: (0, 2, 3)}),
         (ValueError, "atom candidate")),
        (CHFCLBR_GEOMETRY, STAR, Candidates(bonds={(0, 1): ((2, 3, 4), (2,))}),
         (ValueError, "bond candidate")),
        (CHFCLBR_GEOMETRY, STAR, Candidates(bonds={(1, 2): ((0,), (0,))}),
         (ValueError, "bond candidate")),
        (Geometry(CHFCLBR, np.zeros((5, 3))), STAR, CENTER,
         (amchi.IllDefinedParityError, "atom 0")),
    ],
)  # fmt: skip
def test__string__errors(
    geometry: Geometry,
    connectivity: Connectivity,
    candidates: Candidates,
    error: tuple[type[Exception], str],
) -> None:
    """Inputs that don't fit each other, and the named errors."""
    error_type, match = error
    with pytest.raises(error_type, match=match):
        amchi.string(geometry, connectivity, candidates)
