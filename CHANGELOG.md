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
