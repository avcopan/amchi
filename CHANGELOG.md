# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Intended differences from the old AMChI

- Para-stereocenters get `t` and `b` entries: a ring stereocenter whose two ring neighbors are equivalent, as in cis- and trans-1,3-dimethylcyclobutane, 1,4-dimethylcyclohexane and bicyclo[3.3.0]octane, which the old AMChI gives one string. A candidate counts if its neighbors keep exactly one tied pair joined in a ring through it, and the candidates of a ring system are kept if flipping any one of them changes the form. They add Schneider's dependent-chirality invariant to the refinement, and those still tied get their parity during the tie-breaking (D32, D33). The `c` layer can differ too, since the stereo splits the ring atoms before numbering, so the string with `stereo=False` is no longer always the stereo string without its stereo layers. A component without a para-stereocenter keeps its string; 12 of the 600 sampled species change, among them ring NH groups (cis or trans to a substituent, under the old convention that a pyramidal N is a stereocenter) and bicyclic bridgeheads. A pyramidal ring N and every bridgehead count, including the cages whose in/out isomer is physically impossible: the test is topological. The prefix stays `AMChI=2S`.

### Fixed

- A bond candidate of a TS may have 4 neighbors on a side, where its atom swaps both its other neighbors (two breaking and two forming bonds). `Candidates` rejected the candidates that `candidates` found, so 8 TSs of RGD1 raised a `ValidationError`; their strings now equal the old ones, one with a bond stereocenter on that bond (`rgd1:MR_447071_2`).
- The components of a TS are sorted by their breaking and forming bonds before their `m` layer, unlike the old AMChI. Two components equal but for the reacting bonds, as in a TS of RGD1 (`rgd1:MR_242866_0`, `2CH4N2`), were ordered by their `m`, so reflecting the TS swapped them in the `k` and `r` layers. No string of the large sample changes.
- The canonical numbers no longer depend on the atom order where the refinement leaves atoms tied that are not symmetric, as in a C10H10 cage of St. John's dataset (`stjohn:138545_0`), whose `c` layer had three forms. Promoting the last atom of the highest tied class, as the old AMChI does, is replaced by a search over every choice, which keeps the lowest numbered graph and is pruned by the automorphisms it finds. The old choice is tried first and kept on a tie, so no string changes where the old one was consistent (none of the large sample's 21,948 records changes).
- The canonical enantiomer no longer depends on the atom order: the two sides of a bond stereocenter, and stereocenters of equal priority, are compared without reference to the atom indices. An achiral species with bond stereocenters (bicyclo[4.2.2]deca-2,4,7-triene) got an `m` layer for some atom orders; it now has none, as in the old AMChI.

## [0.0.1] - 2026-11-17

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
