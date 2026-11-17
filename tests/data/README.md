# Test data

Written by the old AutoMol (never imported here), except `inchi.json`. The first two files are read by the fixtures in `tests/conftest.py`.

- `records.json.gz`: 900 records (300 closed shells, 300 radicals, 300 TSs) with geometries, the old strings and the old intermediates of each module.
- `hand_cases.json`: HNNO (default and overridden strings), ring-stereo pairs (the kept limitation), and 17 larger molecules.
- `inchi.json`: the standard InChI of the 600 species and the 25 hand cases, made with RDKit from their connectivity. `tests/test__layers.py` checks that the `c` layer writer reproduces each InChI `c` layer (D24).

They are regenerated from the planning repo's `baseline/` (`make_test_records.py`, `dump_hand_cases.py`, `inchi_c_layers.py`).
