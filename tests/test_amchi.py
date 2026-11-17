"""Tests for the public pipeline, against the old strings.

A string equals the old one except for the version prefix and the c layer (D24). The
c layer must describe the same graph in the same canonical numbers as the old one;
the writer itself is tested against InChI in `test__layers.py`.
"""

import json
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest
from numpy.typing import ArrayLike

import amchi
from amchi import Bond, Candidates, Connectivity

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
def transition_states(records: list[dict]) -> list[dict]:
    """The 300 TSs as sampled, each with its connectivity and its string."""
    result = []
    for record in records:
        if record["kind"] == "ts":
            old = record["old"]
            geometry = Geometry(record["symbols"], record["xyz"])
            connectivity = Connectivity(
                bonds=old["bonds"], breaking=old["breaking"], forming=old["forming"]
            )
            candidates = amchi.candidates(record["symbols"], connectivity)
            string = amchi.string(geometry, connectivity, candidates)
            result.append(
                record
                | {"geometry": geometry, "connectivity": connectivity, "amchi": string}
            )
    assert len(result) == 300
    return result


@pytest.fixture(scope="module")
def ts_differences() -> dict[str, dict]:
    """The TSs whose new strings differ from the old ones (D28, D29), by id."""
    cases = json.loads((DATA / "ts_differences.json").read_text())
    return {case["id"]: case for case in cases}


@pytest.fixture(scope="module")
def components() -> list[dict]:
    """Species of the sample put together, written by `baseline/dump_components.py`."""
    return json.loads((DATA / "components.json").read_text())


def test__from_geometry(species: list[dict]) -> None:
    """The 600 species give the old strings (D13), with the c layer of D24."""
    for record in species:
        _assert_old(record["amchi"], record["old"]["amchi"], record["id"])


VARIANTS = ("forward", "reversed", "mirrored", "reversed_mirrored")


def _variant(
    geometry: Geometry, connectivity: Connectivity, variant: str
) -> tuple[Geometry, Connectivity]:
    """A TS reversed (breaking and forming swapped), mirrored (x negated), or both."""
    if "reversed" in variant:
        connectivity = Connectivity(
            bonds=connectivity.bonds,
            breaking=connectivity.forming,
            forming=connectivity.breaking,
        )
    if "mirrored" in variant:
        coordinates = np.asarray(geometry.coordinates) * [-1, 1, 1]
        geometry = Geometry(geometry.symbols, coordinates)
    return geometry, connectivity


def _ts_string(geometry: Geometry, connectivity: Connectivity) -> str:
    candidates = amchi.candidates(geometry.symbols, connectivity)
    return amchi.string(geometry, connectivity, candidates)


@pytest.fixture(scope="module")
def ts_variants(transition_states: list[dict]) -> dict[str, dict[str, str]]:
    """The strings of the four variants of the 300 TSs, by id and variant."""
    return {
        record["id"]: {
            variant: _ts_string(
                *_variant(record["geometry"], record["connectivity"], variant)
            )
            for variant in VARIANTS
        }
        for record in transition_states
    }


def _assert_variants(
    strings: dict[str, str],
    old: dict[str, str],
    difference: dict | None,
    name: str,
) -> int:
    """Each variant gives the old string, or the new one of the differences file.

    Returns the number of new strings.
    """
    differing, new = (
        (difference["differing"], difference["new"]) if difference else ([], {})
    )
    for variant in VARIANTS:
        if variant in differing:
            assert strings[variant] == new[variant], (name, variant)
        else:
            _assert_old(strings[variant], old[variant], f"{name} {variant}")
    return len(differing)


def test__string__ts(
    transition_states: list[dict],
    ts_variants: dict[str, dict[str, str]],
    ts_differences: dict[str, dict],
) -> None:
    """The four variants of the 300 TSs give the old strings, but for 34 (D28, D29).

    That is 1,166 old strings and the 34 new ones of the differences file.
    These include TSs with two components (`/k;3-2,5-1/f;5-3/r;0`) and with the
    canonical direction the reverse of the given one (`/r1`).
    """
    old_variants = {
        x["id"]: x for x in json.loads((DATA / "ts_variants.json").read_text())
    }
    differing = 0
    for record in transition_states:
        name = record["id"]
        assert old_variants[name]["forward"] == record["old"]["amchi"], name
        differing += _assert_variants(
            ts_variants[name], old_variants[name], ts_differences.get(name), name
        )
    assert differing == 34
    assert len(ts_differences) == 11


def test__string__ts__old_tests() -> None:
    """The four variants of the nine TSs of the old tests give the old strings.

    These include radical TSs and halogens, which the sample lacks.
    """
    cases = json.loads((DATA / "ts_hand_cases.json").read_text())
    assert len(cases) == 9
    for case in cases:
        geometry = Geometry(case["symbols"], case["xyz"])
        connectivity = Connectivity(
            bonds=case["bonds"], breaking=case["breaking"], forming=case["forming"]
        )
        strings = {
            variant: _ts_string(*_variant(geometry, connectivity, variant))
            for variant in VARIANTS
        }
        assert _assert_variants(strings, case, None, case["id"]) == 0


def _layer_digits(string: str, layer: str) -> tuple[str, str]:
    """The string with the digits of a one-digit layer removed, and those digits."""
    pattern = rf"/{layer}([0-9;*.]*)(?=/|$)"
    digits = "".join(re.findall(pattern, string))
    return re.sub(pattern, f"/{layer}", string), digits


SWAP_01 = str.maketrans("01", "10")

# The TSs whose mirror image is their own reverse (D29)
MIRROR_IS_REVERSE = {
    "rgd1:MR_35814_1",
    "rgd1:MR_448099_1",
    "rgd1:MR_103507_0",
    "rgd1:MR_323646_1",
}


def test__string__ts__reversed_and_mirrored(
    ts_variants: dict[str, dict[str, str]],
) -> None:
    """Reversing changes the digits of r only; reflecting those of m only.

    A TS whose two directions are equal has no r layer, so reversing changes nothing.
    The exception: a TS whose mirror image is its reverse has no m layer, and
    reflecting it gives the reversed string.
    """
    for name, strings in ts_variants.items():
        forward = strings["forward"]
        rest, r = _layer_digits(forward, "r")
        reversed_rest, reversed_r = _layer_digits(strings["reversed"], "r")
        assert reversed_rest == rest, name
        assert r.translate(SWAP_01) == reversed_r, name
        if name in MIRROR_IS_REVERSE:
            assert "/m" not in forward, name
            assert strings["mirrored"] == strings["reversed"], name
            assert strings["reversed_mirrored"] == forward, name
            continue
        rest, m = _layer_digits(forward, "m")
        mirrored_rest, mirrored_m = _layer_digits(strings["mirrored"], "m")
        assert mirrored_rest == rest, name
        assert m.translate(SWAP_01) == mirrored_m, name
    without_r = {name for name, x in ts_variants.items() if "/r" not in x["forward"]}
    assert without_r == {"rgd1:MR_348096_1", "rgd1:MR_573487_1"}


def _shuffled_ts(
    geometry: Geometry, connectivity: Connectivity, rng: np.random.Generator
) -> tuple[Geometry, Connectivity]:
    order = rng.permutation(len(geometry.symbols))
    index = {int(old): new for new, old in enumerate(order)}

    def renumber(bonds: Sequence[Bond]) -> list[Bond]:
        pairs = ((index[i], index[j]) for i, j in bonds)
        return [(min(pair), max(pair)) for pair in pairs]

    shuffled = Connectivity(
        bonds=renumber(connectivity.bonds),
        breaking=renumber(connectivity.breaking),
        forming=renumber(connectivity.forming),
    )
    coordinates = np.asarray(geometry.coordinates)[order]
    return Geometry([geometry.symbols[k] for k in order], coordinates), shuffled


def test__string__ts__shuffled(transition_states: list[dict]) -> None:
    """Shuffling the atoms of a TS leaves its string unchanged."""
    rng = np.random.default_rng(0)
    for record in transition_states:
        for _ in range(2):
            shuffled = _shuffled_ts(record["geometry"], record["connectivity"], rng)
            assert _ts_string(*shuffled) == record["amchi"], record["id"]


def test__string__ts__stereo_false(transition_states: list[dict]) -> None:
    """Without stereo, a TS string has no b, t, m or s layer.

    It is the stereo string without them, except where the stereo decides the
    direction: then the direction is chosen without it, which can drop `r` or change
    the numbering.
    """
    decided_by_stereo = []
    for record in transition_states:
        string = amchi.string(
            record["geometry"], record["connectivity"], Candidates(), stereo=False
        )
        assert string == _without_stereo(string), record["id"]
        if string != _without_stereo(record["amchi"]):
            decided_by_stereo.append(record["id"])
    assert len(decided_by_stereo) == 8


def test__from_ts_geometries(transition_states: list[dict]) -> None:
    """The geometries of the 300 TSs give the string of their dumped connectivity."""
    for record in transition_states:
        ts, reactant, product = (
            Geometry(record["symbols"], record[key])
            for key in ("xyz", "xyz_reactant", "xyz_product")
        )
        string = amchi.from_ts_geometries(ts, reactant, product)
        assert string == record["amchi"], record["id"]
        string = amchi.from_ts_geometries(ts, reactant, product, stereo=False)
        expected = amchi.string(ts, record["connectivity"], Candidates(), stereo=False)
        assert string == expected, record["id"]


# CH3 + H2O: H4 moves from O5 to C0
H_TRANSFER = Geometry(
    ["C", "H", "H", "H", "H", "O", "H"],
    [[0, 0, 0], [1, 0, 0.3], [-0.5, 0.87, 0.3], [-0.5, -0.87, 0.3],
     [0, 0, -1.3], [0, 0, -2.5], [0.9, 0, -2.8]],
)  # fmt: skip
H_TRANSFER_BONDS = [(0, 1), (0, 2), (0, 3), (5, 6)]


@pytest.mark.parametrize(
    ("breaking", "forming", "expected"),
    [
        ((4, 5), (0, 4), "AMChI=2S/CH5O/c1-3-2/h2H,1H3/k3-2/f3-1/r0"),
        ((0, 4), (4, 5), "AMChI=2S/CH5O/c1-3-2/h2H,1H3/k3-2/f3-1/r1"),
    ],
)
def test__string__ts__h_transfer(
    breaking: tuple[int, int], forming: tuple[int, int], expected: str
) -> None:
    """The transferring H is a node, in the c, k and f layers; reversing changes r."""
    connectivity = Connectivity(
        bonds=H_TRANSFER_BONDS, breaking=[breaking], forming=[forming]
    )
    candidates = amchi.candidates(H_TRANSFER.symbols, connectivity)
    assert amchi.string(H_TRANSFER, connectivity, candidates) == expected


def test__string__ts__candidates() -> None:
    """A TS candidate may list some of its atom's neighbors, but no other atom."""
    connectivity = Connectivity(bonds=H_TRANSFER_BONDS, breaking=[(4, 5)])
    candidates = Candidates(atoms={0: (1, 2, 3)})
    string = amchi.string(H_TRANSFER, connectivity, candidates)
    assert string == amchi.string(H_TRANSFER, connectivity, Candidates())
    with pytest.raises(ValueError, match="atom candidate"):
        amchi.string(H_TRANSFER, connectivity, Candidates(atoms={0: (1, 2, 5)}))
    reverse = Candidates(atoms={0: (1, 2, 3)}, reverse_atoms={0: (1, 2, 6)})
    with pytest.raises(ValueError, match="atom candidate"):
        amchi.string(H_TRANSFER, connectivity, reverse)
    with pytest.raises(ValueError, match="bond candidate"):
        amchi.string(H_TRANSFER, connectivity, Candidates(bonds={(0, 5): ((1,), (6,))}))


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


def test__from_geometry__achiral_bond_stereo() -> None:
    """An achiral species with bond stereocenters gets no m layer in any atom order.

    The enantiomers of bicyclo[4.2.2]deca-2,4,7-triene compare equal; the sides of a
    bond candidate once compared in their stored order (stjohn:122910_00879f).
    """
    rng = np.random.default_rng(0)
    for case in json.loads((DATA / "achiral_cases.json").read_text()):
        geometry = Geometry(case["symbols"], case["xyz"])
        mirror = Geometry(case["symbols"], np.asarray(case["xyz"]) * [-1, 1, 1])
        for variant in (geometry, mirror):
            _assert_old(amchi.from_geometry(variant), case["amchi"], case["id"])
            for _ in range(8):
                string = amchi.from_geometry(_shuffled(variant, rng))
                _assert_old(string, case["amchi"], case["id"])


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
        (CHFCLBR_GEOMETRY, Connectivity(bonds=[(0, 1)], forming=[(0, 5)]),
         Candidates(), (ValueError, "missing atom")),
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
