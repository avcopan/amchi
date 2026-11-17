# Test data

Written by the old AutoMol (never imported here), except `inchi.json`. The first two files are read by the fixtures in `tests/conftest.py`.

- `records.json.gz`: 900 records (300 closed shells, 300 radicals, 300 TSs) with geometries, the old strings and the old intermediates of each module.
- `hand_cases.json`: HNNO (default and overridden strings), ring-stereo pairs (one old string for both isomers, which now differ), and 17 larger molecules.
- `components.json`: four multi-component cases (species of the sample and water, 20 Å apart), for the sorting and joining of components; none of the 600 species has two.
- `inchi.json`: the standard InChI of the 600 species and the 25 hand cases, made with RDKit from their connectivity. `tests/test__layers.py` checks that the `c` layer writer reproduces each InChI `c` layer (D24).
- `ts_variants.json`: the old strings of the four variants (as sampled, reversed, mirrored, both) of the 300 TSs of `records.json.gz`.
- `ts_hand_cases.json`: the nine TSs of the old AutoMol's own tests (radical TSs and halogens included), with their connectivity and the old strings of their four variants.
- `ts_differences.json`: the 11 TSs whose new string differs from the old one in at least one of four variants (as sampled, reversed, mirrored, both), with the old and the new strings (D28, D29). Written by the planning repo's TS prototype (`walkthroughs/2026-11-16-ts-prototype.py --write`), not by the old AutoMol.

- `achiral_cases.json`: `stjohn:122910_00879f` of the large sample (bicyclo[4.2.2]deca-2,4,7-triene, achiral, with two bond stereocenters), whose `m` layer once depended on the atom order. The geometry is from the planning repo's `baseline/oracle_read.py --ids`, the old string from `baseline/sample_large_ids_mismatches.csv`.
- `large_differences.json.gz`: the 190 records of the large sample (21,948 records) whose new string differs from the old one or where a code fails, with the decision that explains each (D34–D36 of the planning repo, or a record that can't be read), the old and new strings and errors, and the geometries. Written by the planning repo's `baseline/dump_large_differences.py`; its `compare_large.py` reads it back.
- `ring_stereo.json`: the 20 geometries of the planning repo's ring-stereo prototype (`walkthroughs/2026-11-17-ring-stereo-geometries.json`, embedded by RDKit, with their SMILES and InChI), and `main`, the string of amchi 0.0.1, before para-stereocenters. Written by the planning repo's `baseline/ring_stereo_data.py`, not by the old AutoMol.

They are regenerated from the planning repo's `baseline/` (`make_test_records.py`, `dump_hand_cases.py`, `dump_components.py`, `inchi_c_layers.py`, `ts_variants.py`, `dump_ts_hand_cases.py`, `dump_large_differences.py`).
