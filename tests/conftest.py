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
