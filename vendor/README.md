# Vendored third-party source

## `pyodine`

An open implementation of the Butler et al. (1996) iodine-cell forward model,
by Paul Heeren, René Tronsgaard-Rasmussen and Frank Grundahl (LSW Heidelberg /
DTU Space / SAC Aarhus), MIT-licensed.  Phase 3 of this project builds its
absolute velocities on it.

| | |
|---|---|
| upstream | `https://github.com/pepeheeren/pyodine` (the docs also point at `https://gitlab.com/Heeren/pyodine`) |
| commit | `4488b0914fe5b272b787982647691045bff2604a` |
| date | 2022-11-04, "Adding pyodine to GitHub" |
| branch | `main` — and this is the repository's **only** commit |
| licence | MIT, preserved at `pyodine/LICENSE` |
| vendored | 2026-09-23, phase 3 prompt 1 |

### This is a fork, not a dependency

From the moment the code landed here it is **our source**.  Not a submodule,
not a `pip install`, not a runtime patch.  Every change to `vendor/pyodine/` is
an ordinary commit with an ordinary justification and a test that fails before
it and passes after.  Where a change is generally correct rather than
HIRES-specific it should also go back to the author.

The reasoning is in `claude_prompts/data_phase3_prompt.md`: `pyodine` is a
published code that demonstrates 0.69 m/s on SONG solar data, and it is also a
one-commit repository with no issues, no releases and no test suite, in which
phase 2 found four things that are wrong for HIRES, two of them silently.  It
is a reference implementation to adapt, not a dependency to trust.

As of this commit the tree is **byte-identical to upstream**.  Verify that, or
see the accumulated diff, with:

```
conda run -n pypeit14 python -m first_hires_exoplanet.vendor_pyodine --verify
```

### What was and was not copied

`first_hires_exoplanet/vendor_pyodine.py` clones upstream, checks that HEAD is
the pinned commit, and copies every git-tracked file, with two adjustments.

**Droppings are dropped.** Upstream commits 133 `__pycache__/*.pyc` files
(bytecode for Python 3.6/3.8/3.9) and five `.ipynb_checkpoints/*.ipynb`
(Jupyter autosaves of the tutorial notebooks).  `.gitignore` would refuse to
commit them anyway, and leaving them on disk would mean the tree here differs
from the tree in the repository — which matters, because that tree is about to
be edited and the diff against upstream is the record of what we changed.  No
source file is affected: 149 of upstream's 287 tracked files are carried.

**Bulk binaries live outside the repository.**  Thirteen files are 303 MB of
the 311 MB upstream tracks: two iodine atlases, the two tutorial tarballs, the
Arcturus reference spectrum, and the CARMENES/HITRAN telluric data.  They are
written to `../first-hires-exoplanet-data/pyodine_assets/` and symlinked back
into `vendor/pyodine/` at their upstream relative paths, because
`utilities_*/conf.py` resolves `../iodine_atlas` relative to itself and the
layout has to survive.  The symlinks are gitignored; `pyodine_assets.csv`
records path, size and SHA-256 for each.  What is committed here is 4.9 MB.

One of those assets is `iodine_atlas/Fischer_Cell_May2022_downsampled3.h5`,
which phase 2 downloaded separately to `../first-hires-exoplanet-data/atlas/`.
The two copies are byte-identical (SHA-256
`3ae788e7…a55a429`), which independently confirms the provenance of the atlas
phase 2 measured `alpha = 2.59` against.

Re-create the whole arrangement from upstream with:

```
conda run -n pypeit14 python -m first_hires_exoplanet.vendor_pyodine
```

### Does it run as shipped?

`first_hires_exoplanet/pyodine_smoke.py` answers this by running upstream's own
tutorial — the σ Dra SONG data, template → observations → velocities — because
upstream ships no test suite at all: no `tests/`, no `conftest.py`, no CI, no
`pytest` configuration.  The `utilities_lick/__pycache__/` directory contains a
compiled `pyodine_parameters_tests` whose source is not in the repository.

Findings are written up in the phase-3 report; the short version is that all 43
modules import cleanly on Python 3.14, and that the science code needs two
NumPy 2 fixes before it will run — `np.float` in `pyodine/lib/misc.py` (three
sites) and `np.NaN` in `pyodine/fitters/lmfit_wrapper.py` (five).  Those are
upstream bugs against a modern NumPy, not HIRES-specific ones, and they belong
in the report that goes back to the author.

### Two things to know before running it

**It swallows its own exceptions.**  `create_template` and
`model_single_observation` wrap their whole body in `except Exception`, write
the traceback to an error log, and return normally.  A stage that failed
outright looks, from the caller, exactly like one that worked.  Read the error
logs; never infer success from the absence of an exception.

**Do not let it plot and then fork.**  The template stage writes diagnostics;
on macOS that initialises Matplotlib's `macosx` backend, and `pathos` then
forks workers into a process holding CoreFoundation.  The workers finish and
write correct output, and the pool then hangs forever — it looks like slow
work, not like a deadlock.  Set `MPLBACKEND=Agg` before anything imports
`matplotlib`.

### Environment

Upstream's declared dependencies are `astropy`, `h5py`, `numpy`, `lmfit`,
`barycorrpy`, `pathos`, `argparse`, `matplotlib`, `progressbar2`, `dill` and
`scipy`.  All are now installed in `pypeit14`.
