# Validation report (Phase 3)

How `amchi` was checked against the old AutoMol AMChI, against itself and against InChI, at the state of `main` on 2026-11-18 (after D41 and D42).
The measurements, scripts and raw tables live in the planning repo [`avcopan/ac`](https://github.com/avcopan/ac): the sections of its [`baseline/README.md`](https://github.com/avcopan/ac/blob/main/baseline/README.md) linked below, and the decision log in [`VISION.md`](https://github.com/avcopan/ac/blob/main/VISION.md#decision-log) (D-numbers).
The old AutoMol is never imported here: every comparison reads strings that it wrote to files.

**Result.** No mismatch is unexplained. Every string that differs from the old AMChI is listed as intended in a decision, and every metamorphic violation found was fixed.

## Summary

| check | records | result | section | decisions behind the exceptions |
|---|---|---|---|---|
| **Oracle, baseline sample** (old code with D41, D42) | 900 (300 closed shells, 300 radicals, 300 TSs) | 881 agree; 19 differ, all explained | [The oracle with D41 and D42](https://github.com/avcopan/ac/blob/main/baseline/README.md#the-oracle-with-d41-and-d42) | D37 ring stereo (12), D28/D29 TSs that break the old code's own rules (7) |
| **Oracle, large sample** (old code with D41, D42) | 21,948 (10,000 stjohn, 5,000 RGD1 TSs, 6,948 reactants and products) | 21,621 agree; 305 differ; 21 fail only in the old code; 1 fails in both; **327 of 327 explained** | [The larger sample](https://github.com/avcopan/ac/blob/main/baseline/README.md#the-larger-sample), [Triage](https://github.com/avcopan/ac/blob/main/baseline/README.md#triage-of-the-larger-sample), [With D41 and D42](https://github.com/avcopan/ac/blob/main/baseline/README.md#the-oracle-with-d41-and-d42) | D34 (160: the old string breaks its own rules), D37 (145: ring stereo), D35 (21: only the old code fails), 1 reader failure |
| Same, old code unchanged (before D41, D42) | 21,948 | 21,597 agree; 324 differ; 23 old-only failures; 4 both fail; 351 of 351 explained | [The oracle with D41 and D42](https://github.com/avcopan/ac/blob/main/baseline/README.md#the-oracle-with-d41-and-d42) | D34, D35, D36, D37 |
| **Test fixtures against the oracle** (`tests/data/low_t_differences.json`) | 340 entries | 296 agree up to D24; 44 differ, all known | [The oracle with D41 and D42](https://github.com/avcopan/ac/blob/main/baseline/README.md#the-oracle-with-d41-and-d42) | D28/D29, D34, D35, D37 |
| **Metamorphic tests**: atom shuffles, reflection, reversal of a TS | 591,865 (full datasets: 289,639 stjohn, 125,234 RGD1 species, 176,992 RGD1 TSs) | 0 violations, and 1 error (the unreadable record); first run: 2 violations (both fixed) and 36 errors (8 `ValidationError` fixed, 27 `IllDefinedParityError` now get strings). Rerun on all 591,865 records after D41 and D42 (2026-11-18, `d5bb51e`) | [Metamorphic tests on the full datasets](https://github.com/avcopan/ac/blob/main/baseline/README.md#metamorphic-tests-on-the-full-datasets) | D38 (cage tie-breaking), D39 (equal TS components), D40 (4 neighbors on a bond side), D41/D42 for the 27 |
| **Diagnostics** of the inputs | 591,865 | After D41 and D42 (rerun 2026-11-18, `d5bb51e`): 1 error (the `ReadError`, 0 `IllDefinedParityError`); 29 TSs with an exactly planar candidate, none a stereocenter; 1.9% with a nearly planar candidate (margin < 0.1); 9 HNNO-pattern records. Before: 28 errors, 62 exactly planar, 13.5% nearly planar, 551 HNNO | [Diagnostics on the full datasets](https://github.com/avcopan/ac/blob/main/baseline/README.md#diagnostics-on-the-full-datasets) | D36 (superseded by D41/D42), D20 |
| **What D41 and D42 change** | 591,865 | 165,480 (28%) change string, 162,207 through D41 and 3,752 through D42; all 27 error TSs get strings | [What D41 and D42 change](https://github.com/avcopan/ac/blob/main/baseline/README.md#what-d41-and-d42-change-over-the-datasets) | D41, D42 |
| **Partition against InChI** (minima) | 414,873 species: 403,227 AMChI groups, 395,937 InChI groups | InChI splits more in 323 groups (751 records), all listed gaps; AMChI splits more in 6,185 groups, all expected | [The partition against InChI](https://github.com/avcopan/ac/blob/main/baseline/README.md#the-partition-against-inchi-minima) | D17/Q7 (0 groups from symmetric rings), D41, D42, the valence cap |
| **Speed** | 600 species, 300 TSs | median 2.6–2.7 ms per species (target 10 ms), 4.1–4.4 ms per TS (target 30 ms); about 130 and 300 times faster than the old code | [Profile of the port](https://github.com/avcopan/ac/blob/main/baseline/README.md#profile-of-the-port) | none; Phase 6 (a Rust port) is not needed |

## What the exceptions are

A string that differs from the old AMChI is a bug unless a decision lists it as intended:

- **D28, D29, D34, D35:** the old code is wrong or ill-defined for the record. Its string changes under an atom shuffle, or on reversing a TS by more than `r`, or it raises. The new string stands.
- **D37:** ring stereo. Ring nitrogens and bridgeheads count as stereocenters, so the new string has a `t` entry that the old one lacks (the prefix stays `AMChI=2S`).
- **D38–D40:** found by the metamorphic tests on the full datasets, so intended differences from the old numbering of tied atoms, the order of equal TS components, and the input check.
- **D41, D42:** the low-temperature principle. A conjugated nitrogen is not a stereocenter, and a 3-ring radical atom is one. They change 28% of the records, so the comparison with the old code applies the same rules to it (`baseline/low_t.py`, without changing the old code).

## The InChI gaps

InChI splits more than AMChI in 323 of 403,227 AMChI groups:

| cause | groups |
|---|---|
| charge-separated stereo bonds (no formal charges in scope) | 191 |
| the valence cap of `amchi.connectivity` (a fourth bond on N is dropped, as in the old code) | 123 |
| parities that InChI reads as undefined | 6 |
| allenes | 3 |
| symmetric rings (Q7) | 0 |

The separate allene and cumulene gap, measured by a parity from the geometry, is 43 groups.

## Limits of the report

- **The full-dataset runs were repeated after D41 and D42** (2026-11-18, `amchi` `d5bb51e`): the metamorphic test finds 0 violations in all 591,865 records and the diagnostics 0 `IllDefinedParityError`, so the limit that they predate those decisions is lifted. The earlier outputs are kept in `baseline/` under `..._2026-11-17-pre-d41`.
- **The oracle ran on a sample.** The old code takes about 0.35 s per species and 1.3 s per TS, so it ran on the 21,948 records of the large sample and the 900 of the baseline sample (D30). The fallback of running it over the full datasets was not needed.
- **Strings are compared, not chemistry.** InChI is a separate check of the partition, and it needs RDKit to assign bond orders, so some of its disagreements come from how its input is made.
- **Known limitations** that no test can remove: no formal charges, the N=N bond of HNNO-type structures is missed by the connectivity (the `add_bond_candidate` correction handles it), and N with 4 neighbors is not reached.
