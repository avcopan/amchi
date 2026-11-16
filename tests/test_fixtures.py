"""Check that the test fixtures load."""

from collections import Counter


def test__records(records: list[dict]) -> None:
    """The baseline has 300 records of each kind, with consistent geometries."""
    assert Counter(r["kind"] for r in records) == {
        "closed": 300,
        "radical": 300,
        "ts": 300,
    }
    assert len({r["id"] for r in records}) == 900
    for r in records:
        assert r["old"]["amchi"].startswith("AMChI=1/")
        assert len(r["symbols"]) == len(r["xyz"])
        if r["kind"] == "ts":
            assert len(r["xyz_reactant"]) == len(r["xyz_product"]) == len(r["xyz"])
            assert "reverse_candidates" in r["old"]


def test__hand_cases(hand_cases: dict[str, list[dict]]) -> None:
    """Each hand-case group is present with geometries and strings."""
    assert {k: len(v) for k, v in hand_cases.items()} == {
        "hnno": 2,
        "ring": 4,
        "ring_control": 2,
        "scaling": 17,
    }
    for case in (c for group in hand_cases.values() for c in group):
        assert len(case["symbols"]) == len(case["xyz"])
        assert case["amchi"].startswith("AMChI=1/")
    assert all("amchi_overridden" in c for c in hand_cases["hnno"])
