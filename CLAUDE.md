# CLAUDE.md

amchi computes the AMChI string of a species (minimum) or a transition state from a 3D structure.
It is a small, standalone rewrite of the AMChI of the old AutoMol, and it is temporary: the Rust StereoMolGraph is meant to replace it.
Read `CONTRIBUTING.md` for the workflow and the design rules (the flat API, the naming rules, the import layers, how the tests use the old AutoMol).
`.claude/skills/` and `.claude/rules/` come from the copier template.

## Commands

- Run everything through Pixi (`pixi run <task>`, or `pixi run -e test python ...`).
  The dependencies are only installed in the Pixi environments, not in the system Python.
- `pixi run pre-commit` runs all checks. A single test: `pixi run -e test pytest tests/test_amchi.py::test__greet`.

## Things that aren't obvious from the code

- **The strings must equal the old `AMChI=1` strings.** A difference is a bug unless it is listed as intended in the changelog.
  The tests read strings that the old AutoMol wrote to files in `tests/data/`. Never import `automol` or `autochem`, here or in the tests.
- **The core takes three explicit inputs:** `string(geometry, connectivity, candidates)`.
  `connectivity()`, `ts_connectivity()` and `candidates()` are helpers that derive the second and third, so a caller can correct them (`add_bond_candidate`, `remove_candidate`) before calling the core.
- **Only `_connect` and `_parity` read coordinates.** The candidates depend on the graph alone.
- **Parities are stored local**, against the neighbor order stored with them. Canonical parities exist only inside `_canon`, derived from the local ones and the current priorities.
- **Hydrogens are explicit in the inputs and implicit in the canonicalized graph**, except backbone hydrogens such as a transferring H in a TS.
- **Validate at the boundary** (the input models and the tests), not inside loops. Build the graph once per record and never copy it in a loop.
- **A parity that isn't well defined, or an unsupported element, raises a named error.** Don't guess.
- Atoms are indices into the geometry, and a bond is an ordered pair `(i, j)` with `i < j`.
- The public API is flat: users write `amchi.string(...)`. A new function goes in a private module and must be added to `amchi/__init__.py`.
- Docstring examples run as doctests and start with `import amchi`, so they read as user code.

## Style

Keep docstrings and comments short: one-line summaries where possible, without restating what a name or type hint already says.
Use the paper's vocabulary (*candidate*, *stereo-determining neighbors*, *priority*, *canonical number*, *local* and *canonical parity*) and full words (`graph`, `priorities`, `candidates`), not the old abbreviations (`gra`, `pri_dct`, `cand_dct`, `bnd_key`).
