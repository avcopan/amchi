# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- `AmchiError`, the base of the errors that amchi raises, and `UnsupportedElementError` for an element without data (amchi supports the main-group elements from H to Xe).
- The input models `Connectivity` (the bonds, and the breaking and forming bonds of a TS) and `Candidates` (each stereocenter candidate with its stereo-determining neighbors), the `GeometryLike` protocol, and `IllDefinedParityError`.
- `parity_margins`, a diagnostic of how nearly planar each candidate is in a geometry. A parity is ill-defined only when it is exactly zero, as in the old AutoMol, so a nearly planar nitrogen still gets one.
- `candidates`, the stereocenter candidates of a species with their stereo-determining neighbors (hydrogens first, then by ascending index), and the edits `add_bond_candidate` and `remove_candidate`, to correct them by hand. Every candidate is kept, also an atom with two or more hydrogens.
- The canonicalization of a species (internal, for the layers to come): the priorities, the canonical parities including higher-order stereocenters, the canonical enantiomer, and the canonical numbers, in the old AutoMol's order, so they equal its values.
- `string(geometry, connectivity, candidates)`, the AMChI string of a species, and `from_geometry(geometry)`, which derives the connectivity and the candidates first. Both take `stereo=False` to leave out the `b`, `t`, `m` and `s` layers. Several components are sorted and joined as in the old AutoMol (`C4H6N.C2H5NO2.H2O`, `2*` for repeated layers). `string` raises `ValueError` when its inputs don't fit each other, and `NotImplementedError` for a transition state, for now.

### Intended differences from the old AMChI

- The strings start with `AMChI=2S`, and the `c` layer is written the way InChI writes its connections layer (one depth-first writer). The old code switched between two writers, one of which depends on Python's set hash order. The `c` layer describes the same graph in the same canonical numbers, and every other layer is unchanged.
