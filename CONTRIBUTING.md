# Contributing to amchi

Contributions are welcome: bug reports, documentation, and code.
For larger changes, please open an issue first to discuss the idea.

## Setup

1. Install [Pixi](https://pixi.prefix.dev/latest/installation/).
2. Fork and clone the repository.
3. In the clone, install the environments and the Git hooks:
   ```bash
   pixi install --all
   pixi run lefthook install
   ```

The project follows the [copier-pixi-python](https://github.com/avcopan/copier-pixi-python) template, whose README explains the tooling in more detail.

## Workflow

| task | what it does |
|---|---|
| `pixi run test` | Run the tests, including the doctests in `src/` |
| `pixi run test-lowest` | Run the tests with the lowest dependency versions allowed by `pyproject.toml` |
| `pixi run pre-commit` | Run every check that runs on commit: format, lint, types, import contracts, both test runs, and a clean-tree check |
| `pixi run pre-commit-check` | The same, without the tests |
| `pixi run imports` | Check the import layers |
| `pixi run docs-preview` | Build and preview the documentation site |
| `pixi run add <pkg>` | Add a run dependency to `pyproject.toml` |

To run a single test, call pytest in the `test` environment:
```bash
pixi run -e test pytest tests/test__element.py::test__element
```

Coverage must stay at or above 80%.
Docstrings use the [Google style](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings), with types taken from the annotations rather than repeated in `Args:` and `Returns:`.
Docstring examples run as doctests, so keep them working.

Run dependency bounds follow [SPEC 0](https://scientific-python.org/specs/spec-0000/), with no upper bounds (see `.claude/rules/dependencies.md`).
`pixi run add` writes `>=<latest>`, so lower the bound afterwards and run `pixi run test-lowest`.

## Design

### What amchi does

amchi turns a 3D structure into an AMChI string: the formula, the `c`, `h`, `b`, `t`, `m` and `s` layers, and, for a transition state, the `k`, `f` and `r` layers.
It reproduces the strings of the AMChI in the old AutoMol (Copan et al., *J. Phys. Chem. A* 2024), written fresh for a smaller scope.
Not in scope: formal charges, stereoisomer enumeration, SMILES, parsing AMChI strings, and the reactant and product stereo of a TS.

### Explicit inputs, separate helpers

The core takes three inputs, and separate helpers derive two of them, so that a caller can correct a perception step instead of working around it:

```python
import amchi

conn = amchi.connectivity(geo)
cands = amchi.candidates(geo.symbols, conn)
# Correct the result by hand: here atom 1 is a planar N=N, not a stereocenter
cands = amchi.remove_candidate(cands, 1)
cands = amchi.add_bond_candidate(cands, conn, (0, 1))
amchi.string(geo, conn, cands)
```

- **Core:** `string(geometry, connectivity, candidates)`.
- **Helpers:** `connectivity(geometry)`, `ts_connectivity(reactant_geometry, product_geometry)` (the breaking and forming bonds, from geometries in the TS atom order), and `candidates(symbols, connectivity)`.
- **Corrections:** `add_bond_candidate` and `remove_candidate`, or a changed bond list.
- **Convenience:** `from_geometry(geometry)` and `from_ts_geometries(ts, reactant, product)` chain the helpers and the core.
- **Inputs** are frozen Pydantic models over plain atom indices, so no graph library is needed to build them:
  `Connectivity` (`bonds`, `breaking`, `forming`), `Candidates` (the stereo-determining neighbors of each atom and bond candidate), and the `GeometryLike` protocol (anything with `symbols` and `coordinates`, N×3 in Å).
  amchi doesn't depend on autogeom.

The candidates depend on the graph alone, so an identity never changes with the conformer, apart from the parities.

### Functions at the root

Every public name is exported from the package root, and the modules are private (`_canon.py`).
Add a new function to a private module and export it in `amchi/__init__.py`.

Naming rules:
- **Functions are named for what they return** (`connectivity`, `candidates`, `string`). A conversion names its source (`from_geometry`, `from_ts_geometries`).
- **Use the paper's vocabulary:** *candidate*, *stereo-determining neighbors*, *priority* (a class number), *canonical number* (after tie-breaking), *local* and *canonical parity*, *canonical direction*, *canonical enantiomer*.
- **Use full words:** `graph`, `priorities`, `candidates`, `atom`, `bond`. The old abbreviations (`gra`, `pri_dct`, `cand_dct`, `bnd_key`, `tsg`) and closures with a trailing underscore (`par_eval_`) are not used.
- **Atoms are indices into the geometry.** A bond is an ordered pair `(i, j)` with `i < j`.
- **The layers keep their AMChI letters.**
- **Avoid the old top-level module names** (`automol`, `autochem`, …).

### Layers

Modules are layered, and a module may only import from the layers below it.
`pixi run imports` enforces this with the contract in `pyproject.toml`:

```
_amchi         the public pipeline functions        (highest)
_layers        the layer strings and the final assembly
_canon         priority refinement, canonical parities and enantiomer, TS direction, canonical numbers
_parity        local parities from the coordinates, and the check that each is well defined
_candidates    tetrahedral atoms and rigid-planar bonds with their neighbors, and the edits
_resonance     the π-system model: which bonds are rigid and planar
_connect       connectivity from a geometry; the reacting bonds of a TS
_graph         the NetworkX graph with validated atom and bond fields
_types         the input models and the GeometryLike protocol
_element       valences, lone pairs and covalent radii       (lowest)
```

Only `_connect` and `_parity` read coordinates.
Each layer is wrapped in parentheses in the contract until its module exists. Remove the parentheses when you add the module.

### How the string is made

The canonicalization is the algorithm of the paper (Schemes 1 and 3; SI section 1), so that the strings stay the same:
1. Refine the priorities with the atom invariants (stable partition refinement, Schneider et al. 2015).
2. Candidates whose stereo-determining neighbors now all differ are stereocenters. Derive their canonical parities, add them to the invariants, and repeat until nothing changes.
3. For a TS, do this for both directions and keep the canonical one.
4. For a chiral molecule, do this for the mirror image and keep the canonical enantiomer.
5. Break the remaining ties to get canonical numbers, relabel, and write the layers.

The mirror image is the same graph with its atom parities inverted.
Inside `_canon`, the neighbor lists and the fixed invariants are pulled out of the graph once, into plain lists and tuples, and the loops run on those.

Rigid planar bonds (the `b` layer) come from the old electron model: neutral atoms with standard valences, and a dominant resonance structure that pairs as many unsaturated neighbors as possible.
It is one yes-or-no question for each candidate bond, answered by a pruned search, and not an enumeration of resonance structures.

### Tests

The regression reference is the old AutoMol.
`tests/data/` holds records (geometries, and for a TS the reactant and product geometries) with the strings and the intermediates that the old code produced: the bonds, the rigid planar bonds, the candidates, the local parities, the priorities, and the string.
Each module is tested alone against its intermediates, and the full pipeline against the strings.
Never import the old code in the tests: they read files.

Also test properties that need no reference: shuffling the atoms leaves the string unchanged, reflecting the coordinates flips `m` and nothing else, and swapping reactant and product flips `r` and nothing else.

A string that differs from the old one is a bug unless the changelog lists it as an intended difference.
Known limitations of the old code are kept, with tests that record today's behavior: a species with formal charges or hypervalent atoms can get the wrong candidates, which a caller can correct.
One is fixed, as an intended difference: the old code gives cis- and trans-1,3-dimethylcyclobutane (and 1,4-dimethylcyclohexane) one string, and amchi finds their para-stereocenters.

### Errors

Report, don't guess: an unsupported element, and a stereocenter whose parity the coordinates don't define, raise a named error.
A parity is ill-defined only when its determinant or dot product is exactly zero, as in the old code: the old strings give nearly planar nitrogens a parity.
How nearly planar each candidate is (a tetrahedral atom, a bond twisted near 90°) is a diagnostic, `parity_margins`, not an error.
Validate at the boundary, in the input models, and not inside loops.
