"""Fixtures: the records that the old AutoMol produced (see tests/data/)."""

import gzip
import json
from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def records() -> list[dict]:
    """The 900 baseline records: geometry, old string, and the old intermediates.

    Each has `id`, `kind` (`closed`, `radical`, `ts`), `symbols`, `xyz` (Å) and `old`,
    a dict of the old intermediates (`amchi`, `bonds`, `candidates`, ...). A TS also has
    `xyz_reactant` and `xyz_product`.
    """
    with gzip.open(DATA / "records.json.gz", "rt") as f:
        return json.load(f)


@pytest.fixture(scope="session")
def hand_cases() -> dict[str, list[dict]]:
    """The hand cases by group: `hnno`, `ring`, `ring_control` and `scaling`.

    Each has `id`, `symbols`, `xyz` (Å) and `amchi`, the old string. The HNNO cases also
    have `amchi_overridden` and the `override` that produces it.
    """
    groups: dict[str, list[dict]] = {}
    for case in json.loads((DATA / "hand_cases.json").read_text()):
        groups.setdefault(case["group"], []).append(case)
    return groups


@pytest.fixture(scope="session")
def ring_ids() -> set[str]:
    """The ids of the 12 species with a para-stereocenter (D32, D33).

    Their old intermediates (priorities, numbers, layers) no longer apply.
    """
    cases = json.loads((DATA / "ring_stereo_differences.json").read_text())
    return {case["id"] for case in cases}


@pytest.fixture(scope="session")
def low_t() -> dict[str, dict[str, dict]]:
    """The records whose strings D41 and D42 change, by fixture and id.

    Each has `decisions` and `new`: the string, or for a TS of `records` or
    `ts_hand_cases` the strings of its four variants. `fixture` is the file name, or
    `hand_cases:<group>`. The old strings stay in their own fixtures.
    """
    result: dict[str, dict[str, dict]] = {}
    for entry in json.loads((DATA / "low_t_differences.json").read_text()):
        result.setdefault(entry["fixture"], {})[entry["id"]] = entry
    return result
