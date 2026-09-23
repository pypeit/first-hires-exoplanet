# Data, Phase 3: the iodine forward model

## Goals

Item 3 of the three-step scope in `context_prompts.md` Q3: **absolute
velocities by forward-modelling the iodine cell**, in the Butler et al. (1996)
tradition, using the inputs phase 2 built.

Phase 2 ended with a measured statement of what this has to beat. Its
cross-correlation reached **702 m/s** epoch to epoch against a **72 m/s**
semiamplitude, and localised the entire shortfall in the wavelength zero point.
Phase 3 replaces the ThAr wavelength solution with one printed onto the science
photons themselves.

## The decision this document makes

Two, and they should be settled here rather than argued in every prompt.

### 1. We fork `pyodine`. It is ours from the moment we take it.

`pyodine` is a published, MIT-licensed implementation of exactly this machinery
and it is the right starting point. It is also a sparse repository — one
commit, no issue traffic, no releases — and phase 2 found four things in it that
are wrong for us, two of which fail silently:

- the iodine depth model is linear where the physics is Beer-Lambert, and at
  the exponent HIRES needs it produces negative transmission;
- the atlas is loaded on the air wavelength grid where PypeIt gives vacuum,
  an 83 km/s error;
- `bary_date` and `bary_vel_corr` are documented in the wrong units and the
  wrong time system, and both would fail silently;
- `compute_weight` knows nothing of a propagated inverse variance.

None of that makes `pyodine` a bad code — its paper demonstrates 0.69 m/s on
SONG solar data. It makes it a **reference implementation to adapt, not a
dependency to trust**, which is what the phase-2 document already said.

So: **vendor it into this repository at a recorded commit, and treat the fork as
our own source from that moment.** Not a submodule, not a `pip install`, not a
patch file applied at runtime. Every change is a normal commit with a normal
justification, and every change carries a test that fails before it and passes
after. Where a change is generally correct rather than HIRES-specific — the
depth model is — it should also go back to the author.

### 2. The target is a detection, not a match to Butler.

It would be easy to set the bar at the 3 m/s Butler et al. achieved, or the
1.2 m/s the modern catalogue reaches on these same frames, and then call
anything less a failure. That would be the wrong bar. Those numbers come from
decades of refinement of a pipeline built around this instrument, with a
deconvolved template, a measured instrumental profile, and calibration
machinery we do not have and are not going to rebuild.

The honest target is **a convincing, independent detection of the published
Keplerian from a modern open reduction**. With 31 cell-in epochs and the period
known, fitting a circular Keplerian gives an amplitude uncertainty of roughly
`sigma_epoch * sqrt(2/N)`:

| per-epoch precision | sigma_K | significance of K = 72 m/s |
|---|---|---|
| 200 m/s | 50.8 m/s | 1.4 sigma — not a detection |
| 100 m/s | 25.4 m/s | 2.8 sigma — suggestive |
| **50 m/s** | **12.7 m/s** | **5.7 sigma — a detection** |
| 30 m/s | 7.6 m/s | 9.4 sigma |
| 10 m/s | 2.5 m/s | 28.3 sigma |

**Phase 3 succeeds at 50 m/s per epoch.** That is a 14x improvement on phase 2
and roughly 17x short of Butler; both of those facts should be stated plainly
in anything we publish. Below 30 m/s we are doing well. Below 10 m/s we should
be suspicious and go looking for the mistake.

Phase 2's per-order scatter *within* one exposure was 47 m/s on the blue
orders. The forward model works on the iodine orders instead, 5000-6200 A,
which are the high-S/N end — template S/N 320-346 against 104-193 in the blue.
So the photon budget supports the target; the question is whether the model
does.

## What we already know

### Carried forward from phase 2

**The data.** 20 nights reduced, 36 science frames, 37 echelle orders (57-93)
covering 3806-6262 A. 31 frames are cell-in and are the velocity epochs; five
are cell-out. Nightly wavelength RMS 0.108-0.146 px across nine months. Summary
in `redux/run_summary.json`; reductions in
`../first-hires-exoplanet-data/redux/reduce_YYYYMMDD/`.

**The template.** `redux/template/hd187123_template_19980826.fits`, the three
1998-08-26 cell-out exposures co-added: **S/N 310 median**, 320-346 in the
iodine region, gain 1.72 against the ideal sqrt(3). Observed frame, with
`VBARY` and `MJDREF` in the primary header. **It is not deconvolved.**

**The atlas.** `../first-hires-exoplanet-data/atlas/Fischer_Cell_May2022_downsampled3.h5`,
4980-6250 A vacuum, FTS grid, R ~ 600,000-680,000, no metadata of any kind. Its
line positions match the HIRES cell at correlation **0.928**. Its depths do not:
the HIRES cell is **2.6x optically thicker**, and a single Beer-Lambert exponent
`alpha = 2.59` describes the difference.

**The adapter.** `first_hires_exoplanet/utilities_hires/`, laid out like
`pyodine`'s `utilities_lick/` so it drops into the fork unchanged. All 36
epochs load; every contract check passes.

**The quality filter.** `data/order_quality.csv`, 1329 order-spectra, 11.3%
rejected, 460 of 539 usable in the iodine region. Note that the rejection rate
in the iodine region (15%) is worse than in the blue (9%) — phase 3 loses more
than phase 2 did.

**The reference frame.** Established in phase 2 prompt 1 and implemented in the
adapter: work in the **observed** frame, divide out `VEL_CORR` unconditionally,
`bary_date` a full JD(UTC) at the exposure midpoint, `bary_vel_corr` in m/s from
`barycorrpy`.

**The precision to beat and the diagnosis.** 702 m/s epoch to epoch, 47 m/s
between the orders of one exposure. The gap is the wavelength zero point.

### The traps that are already known

Each of these has already caught this project once. None should cost time twice.

1. **Air versus vacuum.** PypeIt reports vacuum. The atlas carries both grids
   and `pyodine` reads the air one. 1.56 A at 5615 A is 83 km/s. Phase 1 was
   caught by the same distinction in a different place.
2. **`pyodine`'s linear depth scaling.** `models/spectrum.py:93`. At
   `alpha = 2.6`, 5.8% of the atlas goes negative.
3. **`bary_date` is a full JD in UTC**, not the reduced BJD `components.py:315`
   claims, and **`bary_vel_corr` is m/s**, not the km/s of line 316.
4. **Masked pixels arrive as exact zeros**, not flagged gaps. An all-zero flux
   vector makes `components.Spectrum` raise, so bad orders must be
   zero-*weighted*, never zero-*fluxed*.
5. **The modern catalogue's "BJD" column is JD(UTC)**, not a barycentric date,
   despite its description. Up to 8 minutes out.
6. **PypeIt's heliocentric correction has a sign error** worth +13 m/s and
   drifting 5 m/s over this baseline. Never reuse it; compute the correction
   fresh.
7. **Two nights are calibrated from borrowed arcs** — 1998-07-19 and 1998-09-17
   — and phase 2 measured them 1.5-2.4 km/s off. **Two nights have no pixel
   flat**: 1997-12-23 and 12-24, traced from an iodine-in flat.
8. **A plausible-looking wrong answer is the normal failure mode here**, not a
   crash. Phase 2's list: a +37 km/s annual curve from a sign error, a perfect
   +0.00 km/s agreement from quantisation, a 111 km/s line width from an
   estimator measuring its own window, blank calibration frames that were clean
   by every archive column. Every one was caught by asking whether the
   magnitude was physically possible, or by a check sharing no machinery with
   the first. A forward model has far more ways to look right.

### Open questions worth carrying

- **The instrumental profile.** `pyodine` fits an LSF per chunk; HIRES at
  R ~ 42,000 (measured from the template's weak lines) has an asymmetric
  profile that varies along the order. Which of `pyodine`'s LSF models suits
  this data is not known and is probably the single biggest lever on the
  result.
- **Gain and read noise.** Headers give `CCDGN01 = 4.8` e-/ADU and
  `CCDRN01 = 6.0` e- against PypeIt's hard-coded 1.9 and 2.8. The inverse
  variance in the spec1d files is on PypeIt's scale, so the weights are
  internally consistent, but the *relative* weighting of bright and faint pixels
  is not. Phases 1 and 2 both deliberately left this alone. Phase 3 is the first
  phase with a reason to care.
- **Whether a single `alpha` is enough.** Phase 2 measured 2.59 with a spread
  of 2.28-2.86 from *two 60 s frames*. That is the weakest measurement phase 3
  depends on.
- **The upstream report to Ryan Cooke is still unsent.** Phase 2 prompt 9
  decided it should go and listed eleven items. Deciding is not sending.
- **1998-09-13 gave 34 orders, not 37**, with a worst-order RMS of 0.477 px.
  Never diagnosed.

## Code

See the guidelines in the prompt docs in
`Projects/PypeIt/PypeIt-development-suite/claude_prompts` for how to code in
Python.

Use the `pypeit14` environment. It will need `h5py` and `barycorrpy`, neither of
which is currently installed — that is prompt 1. `first_hires_exoplanet/atlas_check.py`
currently runs in `astro` for exactly this reason and should move once the
install is done.

The vendored `pyodine` goes in `vendor/pyodine/` at the repository root, with
its upstream commit recorded in `vendor/README.md`. It is MIT-licensed;
preserve the licence file.

Existing modules to build on, not duplicate:
`koa_survey.py`, `koa_download.py`, `run_setup.py`, `reduce_run.py`,
`calib_inventory.py`, `order_shift.py`, `quality_filter.py`, `refframe_audit.py`,
`build_template.py`, `measure_velocities.py`, `atlas_check.py`, `figs_phase2.py`,
`check_env.py`, `verify_orig_fixes.py`, and `utilities_hires/`.

Keep raw frames, reductions, the atlas and the template in
`../first-hires-exoplanet-data/`, outside the repository. Small derived
products — velocity tables, diagnostics — belong in
`first_hires_exoplanet/data/` and should be committed.

Record the exact PypeIt commit for any reduction whose results we quote, and the
exact `pyodine` commit anything is derived from.

## Prompts

1. Read this file. Vendor `pyodine` into `vendor/pyodine/` at a recorded
   commit, preserving its licence, and install `h5py` and `barycorrpy` into
   `pypeit14`. **Change nothing in the vendored code yet.** Establish that it
   works as shipped: run whatever tests or examples it has, and report what
   passes, what fails, and what it needs that we do not have. Move
   `atlas_check.py` to `pypeit14`. Use Opus 5. Log your work.

2. Read this file. Make the four HIRES changes to the fork, each as a separate
   commit with a test that fails before and passes after: the Beer-Lambert
   depth model, the vacuum atlas grid, the two `components.py` docstrings, and
   whatever `compute_weight` needs to accept a propagated inverse variance.
   Report which of these are HIRES-specific and which should go back to the
   author. Use Opus 5. Log your work.

3. Read this file. Deconvolve the stellar template against the instrumental
   profile, using `pyodine`'s own machinery (`template/deconvolve.py`). Report
   the S/N cost of deconvolution, how the result depends on the assumed LSF,
   and whether the deconvolved template reproduces the observed one when
   re-convolved. Use Opus 5. Log your work.

4. Read this file. Write the remaining `utilities_hires` files — `conf.py`,
   `pyodine_parameters.py`, `timeseries_parameters.py`, `logging.json` —
   modelled on `utilities_lick`. Then fit **one chunk of one order of one
   epoch** and report every fitted parameter with its uncertainty. Do not fit
   more than one chunk. Use Opus 5. Log your work.

5. Read this file. Fit **one epoch** end to end, across all usable iodine
   orders. Report the per-chunk velocity scatter, how many chunks fail, what
   the LSF looks like as a function of position, and the runtime. Compare the
   fitted continuum and wavelength solution against PypeIt's. Use Opus 5. Log
   your work.

6. Read this file. Settle the instrumental profile. Compare `pyodine`'s
   available LSF models on the single epoch from prompt 5, and choose one on
   evidence. Report how much the per-chunk scatter depends on the choice —
   this is expected to be the largest single lever on the final precision.
   Use Opus 5. Log your work.

7. Read this file. Fit all 31 cell-in epochs. Produce a velocity table with
   honest uncertainties in `first_hires_exoplanet/data/`, reporting the
   per-chunk, per-order and epoch-to-epoch scatter separately, as phase 2
   prompt 5 did. Use Opus 5. Log your work.

8. Read this file. Assess what prompt 7 produced against the published
   3.097-day Keplerian and the modern catalogue, exactly as phase 2 prompt 6
   did so the two are comparable. State the achieved precision, the
   significance of the detection, and what the residuals contain. Generate
   figures with a script on disk. Use Opus 5. Log your work.

9. Read this file. Diagnose what is left. Do the two borrowed-arc nights
   (1998-07-19, 1998-09-17) come back into line, as they should if the forward
   model derives its own wavelength solution? Do the intra-night drifts on
   08-25 and 08-26 disappear? What happens to the two December nights with no
   pixel flat, and to 1998-09-13? Report what the forward model fixed and what
   it did not. Use Opus 5. Log your work.

10. Read this file. Send the upstream report to Ryan Cooke: the eleven items
    listed in the phase-2 assessment plus anything phase 3 has added, with
    `refframe_audit.py` and `verify_orig_fixes.py` as the reproducers. Draft it
    first and show it to me before anything is sent. Use Opus 5. Log your work.

11. Read this file. Write the phase-3 assessment, and update
    `docs/public_HD187123b.md` with what the three phases actually found —
    including, plainly, where we fell short of the 1998 result and why.
    Use Opus 5. Log your work.

## Logging

After finishing a task, append a dated entry under `## Logs` below, in the
format phases 1 and 2 established:

```
### YYYY-MM-DD (Prompt N: one-line summary of what changed)

**Task.** What the prompt asked for, in a sentence or two, and a pointer to the
Report section for the findings.

**What was done.** The concrete actions: files written, commands run, data
produced.

**Headline results.** The numbers, in a short paragraph.

**What this taught us about the repository and the data.** The part that
matters most. Not what was done, but what is now known that was not known
before — including things that turned out to be wrong, estimators that failed,
and traps that would catch the next person. Phase 2's logs are the model.

**Files added / modified.**

**Whether any git command changed state.** (It should not; git is the user's.)
```

The findings themselves go in `## Report`, under a `### Prompt N: title`
heading. The log records the work and the lessons; the Report records the
result.

## Q&A

## Report

### Prompt 1: the fork is in, and it does not run as shipped

**The fork.** `pyodine` is vendored at `vendor/pyodine/`, from
`https://github.com/pepeheeren/pyodine` at commit
`4488b0914fe5b272b787982647691045bff2604a` — `main`, 2022-11-04, "Adding
pyodine to GitHub", still the repository's only commit, and the same SHA phase
2 read through the GitHub API. The MIT licence is preserved at
`vendor/pyodine/LICENSE`. `first_hires_exoplanet/vendor_pyodine.py` performs
the vendoring and re-runs as a verifier; as of this prompt it reports **149
files identical to upstream, none changed, none missing, none extra**.

Of upstream's 287 tracked files, 149 are carried. The 138 dropped are 133
`__pycache__/*.pyc` (bytecode for Python 3.6/3.8/3.9) and five
`.ipynb_checkpoints/*.ipynb` (Jupyter autosaves of the tutorial notebooks).
No source file is affected.

Of those 149, thirteen are 303 MB of the 311 MB upstream tracks — two iodine
atlases, two tutorial tarballs, the Arcturus reference spectrum, and the
CARMENES/HITRAN telluric data. They are written to
`../first-hires-exoplanet-data/pyodine_assets/` and symlinked back into
`vendor/pyodine/` at their upstream relative paths, because
`utilities_*/conf.py` resolves `../iodine_atlas` relative to itself. The
symlinks are gitignored; `vendor/pyodine_assets.csv` records path, byte count
and SHA-256 for each. **4.9 MB is what the repository actually gains.**

A free check fell out of that: upstream's
`iodine_atlas/Fischer_Cell_May2022_downsampled3.h5` is byte-identical to the
copy phase 2 downloaded separately (SHA-256 `3ae788e7…a55a429`). The atlas
phase 2 measured `alpha = 2.59` against is the atlas `pyodine` ships.

**The environment.** `pypeit14` is Python 3.14.6, NumPy 2.5.3, SciPy 1.18.1,
Astropy 8.0.1, Matplotlib 3.11.2. Installed here: `h5py` 3.16.0 and
`barycorrpy` 0.4.4 as the prompt asked, which pulled in `astroquery` 0.4.11,
`jplephem`, `pyvo`, `keyring`, `beautifulsoup4` and `html5lib`.

Four more of upstream's declared dependencies were missing and had to be
installed too, because without them `pyodine` does not import at all, let alone
run: **`lmfit` 1.3.4, `pathos` 0.3.5, `progressbar2` 4.6.0 and `dill` 0.4.1**
(plus `asteval`, `uncertainties`, `multiprocess`, `pox`, `ppft`,
`python-utils`). `lmfit` is the fitter, `pathos` the multiprocessing layer;
neither is optional. The prompt named only two packages, and that turned out to
be an underestimate of what "install the dependencies" means here.

**Does it work as shipped? Partly, and it says nothing when it does not.**

Upstream ships **no test suite**: no `tests/`, no `conftest.py`, no CI
configuration, no `pytest` settings. The only executable check in the
repository is its own tutorial — SONG spectra of σ Draconis, reduced in three
stages. `first_hires_exoplanet/pyodine_smoke.py` is that tutorial turned into
something re-runnable, and it is the baseline any change in prompt 2 has to
keep passing.

*What passes.* All **43 modules import cleanly** on Python 3.14 — 33 `pyodine`
submodules, four `utilities_*` packages, six top-level driver scripts. The only
complaints are six `SyntaxWarning`s for invalid escape sequences in
`plot_lib.py` (`'$\AA$'`, `'$\pm$'`).

*What fails.* The science code does not run at all against NumPy 2, for two
reasons and, as far as a tree-wide grep can tell, only two:

| removed name | file | sites |
|---|---|---|
| `np.float` | `pyodine/lib/misc.py` | 226, 227, 235 |
| `np.NaN` | `pyodine/fitters/lmfit_wrapper.py` | 155, 163, 172, 189, 206 |

`np.float` is hit immediately: `rebin` is called from `guess_velocity` in the
first stage of template creation. Both aliases were removed in NumPy 2.0.
Neither is HIRES-specific; both are upstream bugs against any modern NumPy and
belong in the prompt-10 report to the author.

*What runs once those two names are restored.* Prompt 1 is forbidden from
changing the vendored code, so `pyodine_smoke.py --numpy-shim` restores
`np.float` and `np.NaN` at runtime, from outside the fork, purely to find out
what else is broken. Nothing else is: the whole tutorial then runs end to end.

| stage | result | time |
|---|---|---|
| imports | 43 modules | 1.1 s |
| template | 528 deconvolved chunks | 118 s |
| observations | 2 epochs × 528 chunks | 114 s (4 cores) |
| velocities | 2-epoch timeseries | 1.5 s |

The first fitted numbers, for the record: chunk velocities on the first σ Dra
epoch have median 1173.9 m/s and scatter 910.8 m/s over 528 chunks, with none
non-finite, and a median reduced χ² of **150516**. That last number is not a
fit-quality statement about the tutorial data; it is `pyodine` dividing by
weights that are not a variance. It is the same hole phase 2 found in
`compute_weight`, seen from the other end, and it means **reduced χ² cannot be
used as a fit diagnostic until prompt 2 gives the weights a real inverse
variance**.

**Phase 2's four findings, confirmed in the vendored tree.** All four are
present at the expected places, so prompt 2 has its targets:

| finding | location in `vendor/pyodine/` |
|---|---|
| linear, not Beer-Lambert, iodine depth | `pyodine/models/spectrum.py:93` |
| atlas read on the air grid | `utilities_lick/load_pyodine.py:35` (`h['wavelength_air']`) |
| `bary_date` / `bary_vel_corr` docstrings | `pyodine/components.py:315-316` |
| `compute_weight` knows only 'flat' and 'inverse' | `pyodine/components.py:155`, and the order loop at `:277` |

**`atlas_check.py` now runs in `pypeit14`.** It was in `astro` only because
`h5py` was missing. It reproduces phase 2's numbers unchanged: `alpha` median
2.59 with spread 2.28-2.86 over seven windows, and order 65 correlating +0.903
against the vacuum grid versus +0.078 against air.

**Four traps this prompt hit, none of them obvious.**

1. **`pyodine` swallows every exception.** `create_template`
   (`pyodine_create_templates.py:620`) and `model_single_observation` each wrap
   their entire body in `except Exception`, write "Something went wrong!" and
   the traceback to an error log, and **return normally**. The caller is told
   nothing. The `np.float` failure surfaced only as a missing output file two
   function calls later, and a stage that has failed outright is otherwise
   indistinguishable from one that worked. Every phase-3 driver must read the
   error logs; `pyodine_smoke.py` does, and reports the exception line rather
   than the last line of the log, because NumPy appends several lines of
   migration advice after the exception.
2. **Plot, then fork, and the run hangs forever.** The template stage writes
   diagnostics, which on macOS initialises Matplotlib's `macosx` backend and
   therefore CoreFoundation; the observation stage then forks four `pathos`
   workers off that process. The workers complete their fits and write correct
   result files, and the pool then never joins — 18 minutes and still waiting,
   with 222 "the process has forked and you cannot use this CoreFoundation
   functionality safely" warnings. It looks exactly like slow work.
   `MPLBACKEND=Agg` before any import removes it. It appears only when
   plotting and forking happen in the same process, which is why running the
   two stages separately the first time hid it completely.
3. **`conda run` buffers.** Without `--no-capture-output` nothing is printed
   until the process exits, so a hang and a long computation look identical.
   Both of the above were diagnosed only after that was fixed.
4. **A stock Python `.gitignore` eats vendored source.** The `lib/` rule —
   meant for build output — silently excluded `vendor/pyodine/pyodine/lib/`,
   five real source files including `misc.py`, the very file with the
   `np.float` bug. `.gitignore` now negates it. Anything vendored into a
   repository with a boilerplate ignore file should be checked with
   `git check-ignore` before it is trusted to be committed.

One smaller thing worth knowing for prompts 5 and 7: `load_results` on an
`.h5` file returns a **dict of arrays**, one entry per fitted parameter, each
of length `n_chunks` — not the list of result objects the tutorial's `dill`
path produces. The tutorial only ever shows the `dill` form.

## Logs

### 2026-09-23 (Prompt 1: `pyodine` vendored at a recorded commit; it needs two NumPy 2 fixes before it runs)

**Task.** Vendor `pyodine` into `vendor/pyodine/` at a recorded commit with its
licence, install `h5py` and `barycorrpy` into `pypeit14`, change nothing in the
vendored code, establish whether it works as shipped, and move `atlas_check.py`
to `pypeit14`. Findings in `## Report`, "Prompt 1: the fork is in, and it does
not run as shipped".

**What was done.** Wrote `first_hires_exoplanet/vendor_pyodine.py`, which
clones upstream, checks HEAD against the pinned SHA, copies every tracked file
except committed bytecode and notebook checkpoints, externalises the thirteen
files over 1 MiB to `../first-hires-exoplanet-data/pyodine_assets/` with
symlinks back into place, and writes `vendor/pyodine_assets.csv`; it also runs
as `--verify`, which currently reports 149 files identical to upstream. Wrote
`vendor/README.md` (provenance, what was and was not copied, how to verify) and
added `.gitignore` rules for the thirteen symlinks plus a negation for the
`lib/` rule that was swallowing `vendor/pyodine/pyodine/lib/`. Installed
`h5py` and `barycorrpy`, and then `lmfit`, `pathos`, `progressbar2` and `dill`,
without which nothing imports. Wrote
`first_hires_exoplanet/pyodine_smoke.py`, which runs upstream's own σ Dra
tutorial — imports, template, observations, velocities — because upstream ships
no test suite; ran it both as shipped and with a runtime NumPy shim. Updated
the environment note in `atlas_check.py` and ran it in `pypeit14`.

**Headline results.** The fork is at `4488b0914fe5b272b787982647691045bff2604a`
and is byte-identical to upstream. 4.9 MB is committed; 303 MB sits outside the
repository with a SHA-256 manifest. All 43 modules import on Python 3.14. The
science code does **not** run as shipped: `np.float` (`lib/misc.py`, three
sites) and `np.NaN` (`fitters/lmfit_wrapper.py`, five) were removed in NumPy 2,
and a tree-wide grep finds no others. With those two names restored at runtime
the tutorial runs end to end — 528 template chunks in 118 s, two epochs modelled
in 114 s on four cores, velocities in 1.5 s — giving a first fitted number of
910.8 m/s chunk scatter on a σ Dra epoch, at a median reduced χ² of 150516.
`atlas_check.py` reproduces phase 2's `alpha = 2.59` and +0.903/+0.078
vacuum-versus-air correlations from `pypeit14`.

**What this taught us about the repository and the data.**

*The code fails silently by construction.* `create_template` and
`model_single_observation` catch `Exception` around their whole body, log it,
and return normally. The `np.float` crash reached us as a `FileNotFoundError`
on the output file, two calls downstream — a textbook instance of the
phase-3 warning that a plausible-looking wrong answer is the normal failure
mode here. Nothing in phase 3 may infer success from the absence of an
exception; the error logs are the only signal, and the smoke script now reads
them.

*Reduced χ² is currently meaningless.* 150516 on a fit that visibly worked is
not a fit-quality statement, it is `compute_weight` dividing by weights that
are not a variance. This is phase 2's fourth finding, seen from the results
end, and it means the obvious fit diagnostic is unavailable until prompt 2
fixes the weights. Worth remembering before anyone reads a χ² in prompt 5.

*Plotting and forking do not mix on macOS, and the symptom is a hang, not a
crash.* Diagnostic plots initialise CoreFoundation; `pathos` then forks into
it; the workers finish and write correct output and the pool never joins. It
cost roughly forty minutes here, twice, and it will recur in prompts 5 and 7
where the same driver plots and forks. `MPLBACKEND=Agg` at the top of any
driver, before `matplotlib` is imported anywhere.

*Diagnosing any of this needs `conda run --no-capture-output`.* Otherwise
output appears only at exit and a deadlock is indistinguishable from slow work.
This applies to every long-running command in this project, not just `pyodine`.

*Vendoring into a repository with a boilerplate `.gitignore` is not safe by
default.* The stock `lib/` rule excluded five source files of the vendored
tree, including the file holding the `np.float` bug. Had that gone uncaught,
the committed fork would have been missing code that the working tree had, and
prompt 2's "test that fails before and passes after" would have been written
against a file nobody else could see. `git check-ignore` over the whole
vendored tree is now part of the recipe.

*Upstream's own asset is upstream's own atlas.* The Fischer cell file `pyodine`
ships is byte-identical to the one phase 2 downloaded and measured against, so
`alpha = 2.59` and the 0.928 line-position correlation are statements about the
same bytes the forward model will use.

*Upstream is sparser than its paper suggests.* One commit, no releases, no
issues, no tests, and a `utilities_lick/__pycache__/pyodine_parameters_tests`
whose source was never committed. The tutorial is the entire executable
specification of the code's behaviour, which is why it is now a script in this
repository rather than five notebooks in a `docs/` directory.

**Files added / modified.** Added `first_hires_exoplanet/vendor_pyodine.py`,
`first_hires_exoplanet/pyodine_smoke.py`, `vendor/README.md`,
`vendor/pyodine_assets.csv`, and `vendor/pyodine/` (136 committed files, 13
gitignored symlinks). Modified `.gitignore` (thirteen asset paths, the `lib/`
negation) and the environment note in `first_hires_exoplanet/atlas_check.py`;
`docs/figs/fig_p2_atlas.png` was rewritten by re-running `atlas_check.py` in
the new environment, with the same content. Outside the repository:
`../first-hires-exoplanet-data/pyodine_assets/` (303 MB of upstream binaries)
and `../first-hires-exoplanet-data/pyodine_tutorial/` (unpacked tutorial data
and its outputs, scratch).

**Whether any git command changed state.** No. `git clone` of upstream, into a
scratch directory outside the repository, and read-only `git ls-files`,
`git rev-parse`, `git status` and `git check-ignore`. Nothing was staged,
committed or branched in this repository.
