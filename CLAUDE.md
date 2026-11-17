# CLAUDE.md

amchi computes the AMChI string of a species (minimum) or a transition state from a 3D structure.
It is a small, standalone rewrite of the AMChI of the old AutoMol, and it is temporary: the Rust StereoMolGraph is meant to replace it.
Read `CONTRIBUTING.md` for the workflow and the design rules (the flat API, the naming rules, the import layers, how the tests use the old AutoMol).
`.claude/skills/` and `.claude/rules/` come from the copier template.

## Commands

- Run everything through Pixi (`pixi run <task>`, or `pixi run -e test python ...`).
  The dependencies are only installed in the Pixi environments, not in the system Python.
- `pixi run pre-commit` runs all checks. A single test: `pixi run -e test pytest tests/test__element.py::test__element`.

## Things that aren't obvious from the code

- **The strings must equal the old `AMChI=1` strings.** A difference is a bug unless it is listed as intended in the changelog.
  The tests read strings that the old AutoMol wrote to files in `tests/data/`. Never import `automol` or `autochem`, here or in the tests.
- **The core takes three explicit inputs:** `string(geometry, connectivity, candidates)`.
  `connectivity()`, `ts_connectivity()` and `candidates()` are helpers that derive the second and third, so a caller can correct them (`add_bond_candidate`, `remove_candidate`) before calling the core.
- **Only `_connect` and `_parity` read coordinates.** The candidates depend on the graph alone.
- **Parities are stored local**, against the neighbor order stored with them. Canonical parities exist only inside `_canon`, derived from the local ones and the current priorities.
- **Hydrogens are explicit in the inputs and implicit in the canonicalized graph**, except backbone hydrogens such as a transferring H in a TS.
- **A bond's local stereo is a table:** one parity for each pair of neighbors across the bond (`BondStereo.parities`), and `_canon` looks up the entry of the highest-priority neighbor at each end.
  Don't derive it by flipping one parity: a TS bond can have three neighbors at one end and needn't be planar, and the flip is then wrong.
- **A TS is canonicalized as the smallest of up to four forms:** each direction, as it is and reflected (`_canon`, a 9-part key).
  The `r` and `m` layers say how the given TS relates to the kept form, so reversing changes only `r` and reflecting only `m`.
  The exception is a TS whose mirror image is its own reverse: it has no `m`, and reflecting it changes `r`. Don't choose the direction first and the enantiomer afterwards, as the old code did: the two passes disagree for some TSs.
- **A TS atom can have two neighbor lists:** `Candidates.atoms` for one side and `reverse_atoms` for the other, each with its own measured parity.
  Reading the TS reversed uses `reverse_atoms`. Don't derive one parity from the other: a substitution site inverts only about half the time.
- **Validate at the boundary** (the input models and the tests), not inside loops. Build the graph once per record and never copy it in a loop.
- **An unsupported element raises a named error.** Don't guess.
- **A local parity is ill-defined (`None`) only when its determinant or dot product is exactly zero**, as in the old code, so a nearly planar nitrogen still gets a parity.
  `_canon` raises `IllDefinedParityError` only when a candidate with parity `None` turns out to be a stereocenter. How nearly planar a candidate is, is a diagnostic: `parity_margins`.
- Atoms are indices into the geometry, and a bond is an ordered pair `(i, j)` with `i < j`.
- The public API is flat: users write `amchi.string(...)`. A new function goes in a private module and must be added to `amchi/__init__.py`.
- Docstring examples run as doctests and start with `import amchi`, so they read as user code.

## Style

Keep docstrings and comments short: one-line summaries where possible, without restating what a name or type hint already says.
Use the paper's vocabulary (*candidate*, *stereo-determining neighbors*, *priority*, *canonical number*, *local* and *canonical parity*) and full words (`graph`, `priorities`, `candidates`), not the old abbreviations (`gra`, `pri_dct`, `cand_dct`, `bnd_key`).
