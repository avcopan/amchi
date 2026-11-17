# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- A user guide page (`user_guide/using-amchi.qmd`), with a section on transition states, and doctests for every public name.
- `AmchiError`, the base of the errors that amchi raises, and `UnsupportedElementError` for an element without data (amchi supports the main-group elements from H to Xe).
- The input models `Connectivity` (the bonds, and the breaking and forming bonds of a TS) and `Candidates` (each stereocenter candidate with its stereo-determining neighbors), the `GeometryLike` protocol, and `IllDefinedParityError`.
- `parity_margins`, a diagnostic of how nearly planar each candidate is in a geometry. A parity is ill-defined only when it is exactly zero, as in the old AutoMol, so a nearly planar nitrogen still gets one.
- `candidates`, the stereocenter candidates of a species with their stereo-determining neighbors (hydrogens first, then by ascending index), and the edits `add_bond_candidate` and `remove_candidate`, to correct them by hand. Every candidate is kept, also an atom with two or more hydrogens.
- The canonicalization of a species (internal, for the layers to come): the priorities, the canonical parities including higher-order stereocenters, the canonical enantiomer, and the canonical numbers, in the old AutoMol's order, so they equal its values.
- `string(geometry, connectivity, candidates)`, the AMChI string of a species, and `from_geometry(geometry)`, which derives the connectivity and the candidates first. Both take `stereo=False` to leave out the `b`, `t`, `m` and `s` layers. Several components are sorted and joined as in the old AutoMol (`C4H6N.C2H5NO2.H2O`, `2*` for repeated layers). `string` raises `ValueError` when its inputs don't fit each other.
- `ts_connectivity(reactant_geometry, product_geometry)`, the bonds of a transition state from reactant-side and product-side geometries in the TS atom order: the bonds of both sides, and the breaking and forming bonds as their difference.
- `Candidates.reverse_atoms`, the product-side neighbors of a TS atom that is tetrahedral on both sides with different neighbors (D26). `remove_candidate` removes an atom's reverse neighbors with it.
- `candidates` and `add_bond_candidate` for a transition state: the rules apply to the reactant and the product side; an atom gets the neighbors of the side where it is tetrahedral (the reactant side first), and also `reverse_atoms`; a bond gets its neighbors from every side where it is rigid and planar, and the small-ring test runs over all the TS's bonds through those neighbors only.
- `ParityMargins.reverse_atoms`, the margins of a TS atom's reverse neighbors, which get their own local parity. A bond's margin is now the smallest over every pair of neighbors across it, since a TS bond's neighbors needn't lie in one plane.
- `string` for a transition state: a connectivity with breaking or forming bonds gives the `k` (breaking) and `f` (forming) layers as read in the canonical direction, and an `r` layer (`0` if the given direction is the canonical one, `1` if it is the reverse, none if the two are equal). The layers are written in the order `c`, `h`, `b`, `t`, `m`, `s`, `k`, `f`, `r`. A TS candidate may list only some of its atom's neighbors, since its neighbors come from one side or several.
- `from_ts_geometries(geometry, reactant_geometry, product_geometry)`, the AMChI string of a transition state with its connectivity (`ts_connectivity`) and candidates derived; `stereo=False` leaves out the `b`, `t`, `m` and `s` layers.

### Intended differences from the old AMChI

- The strings start with `AMChI=2S` (D24), and the `c` layer is written the way InChI writes its connections layer (one depth-first writer). The old code switched between two writers, one of which depends on Python's set hash order. The `c` layer describes the same graph in the same canonical numbers, and every other layer is unchanged.
- A transition state is canonicalized as the smallest of four forms (each direction, as it is and reflected; D29), with one more atom invariant, the sorted pairs (bond type, neighbor priority) (D28). Where the old code was consistent, the strings are unchanged; 11 of the 300 sampled TSs differ in at least one of their four variants, where the old string depended on the atom order or on which of its two passes decided. A TS whose mirror image is its own reverse has no `m` layer, and reflecting it changes `r`.
