# Test data

Written by the old AutoMol (never imported here). Both files are read by the fixtures in `tests/conftest.py`.

- `records.json.gz`: 900 records (300 closed shells, 300 radicals, 300 TSs) with geometries, the old strings and the old intermediates of each module.
- `hand_cases.json`: HNNO (default and overridden strings), ring-stereo pairs (the kept limitation), and 17 larger molecules.

They are regenerated from the planning repo's `baseline/` (`make_test_records.py`, `dump_hand_cases.py`).
