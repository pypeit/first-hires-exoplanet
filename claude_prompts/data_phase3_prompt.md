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

## Q&A

### After prompt 3: decisions before prompt 4

The evidence for each is in `## Report`, "Prompt 3".

1. **Which deconvolved template goes forward?** Prompt 3 wrote two, both in
   `../first-hires-exoplanet-data/redux/template/`:
   - `…_deconv_os10.h5`, at `pyodine`'s default oversampling. It does not
     reproduce the observation (χ² 5.2 against the measured noise; 390 of 700
     chunks above 4). The spline in `misc.rebin` is the cause.
   - `…_deconv_os1.h5`, with no oversampling. It reproduces the observation
     (χ² 1.1). But `pyodine` interpolates the template again when it
     evaluates the model, and could bring the same error back there.

   The options:
   - (a) carry both into prompts 4–5 and let the one-epoch fit decide;
   - (b) change the fork now, so the deconvolver fits only at detector pixels
     instead of the splined fine grid. This is a fork change, so it comes with
     its own test, and it should also go upstream;
   - (c) use oversampling 1 and change nothing.

   My recommendation: (a) now, and (b) if prompt 5 shows the model-side
   interpolation hurts.

   Answer: (a)

2. **Should prompt 4 fix the adapter's zero wavelengths?** Every order of
   every spec1d reaches `pyodine` with wavelength exactly 0 at masked pixels
   (950 pixels in one epoch). `pyodine`'s velocity guess divides by the first
   and last wavelength of each order, so prompt 5 will fail at its first step.
   The fix already exists in `deconvolve_template.fill_wavelengths`: fill
   those pixels from a smooth polynomial and keep them at zero weight. It
   would move into `utilities_hires/load_pyodine.py`, with a test. Prompt 4
   rewrites the `utilities_hires` files anyway.

   My recommendation: yes, in prompt 4.

   Answer: yes, in prompt 4.

3. **When should the B-star and iodine-flat frames be downloaded?** The
   stellar LSF, as opposed to the ThAr upper bound, should be measured from
   B stars observed through the cell. The template night has them:
   - HR 7236 at 08:21 UT and HR 8634 at 10:20 UT, bracketing the template
     exposures (09:39–09:59 UT). There are more HR 8634 and HR 838 frames
     later that night.
   - An iodine-in quartz flat at 15:32 UT.

   None is on disk. Downloading and reducing them with `koa_download.py` /
   `reduce_run.py` is new data acquisition. Should it happen now, as its own
   step, or wait for prompt 6? And is the 1998-08-26 night enough, or do you
   want the B stars from every night?

   My recommendation: 1998-08-26 only, downloaded before prompt 5, so that
   prompt 6 has them ready.

   Answer: follow your recommendation

4. **How should the noise scale be handled?** PypeIt's inverse variance is
   2.35× too large in σ; the exposure-to-exposure scatter is 0.426 of the
   propagated σ. The header gain (4.8 against PypeIt's 1.9 e⁻/ADU) explains
   1.59× of that, and the remaining 1.5× is unexplained. Every reduced χ²,
   fitter weight and per-chunk error in prompt 5 inherits it. The options:
   - (a) rescale the inverse variance by 1/0.426² in the adapter;
   - (b) leave the weights alone, since relative weights are unaffected if
     the factor is uniform, and correct only χ² and the error bars;
   - (c) re-reduce with the header gain and read noise, and find the other
     1.5× first.

   My recommendation: (b) for phase 3, with the factor recorded. (c) is a
   PypeIt question for the prompt-10 upstream report, not something to block
   on.

   Answer: (b)

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

### Prompt 2: five changes to the fork, five tests, and one that was not asked for

**What changed, and what each is worth.** Six units of work, each with its own
test file, each test written and run *before* the change and shown to fail for
the right reason.

| # | change | files in the fork | test | HIRES-specific? |
|---|---|---|---|---|
| 0 | NumPy 2 compatibility | `lib/misc.py`, `fitters/lmfit_wrapper.py` | `test_fork_numpy2.py` | **no** — upstream |
| 1 | Beer-Lambert iodine depth | `models/spectrum.py` | `test_fork_iodine_depth.py` | **no** — upstream |
| 2 | recorded atlas wavelength frame | `components.py` | `test_fork_atlas_frame.py` | mechanism no, choice yes |
| 3 | `bary_*` units, documented and guarded | `components.py` | `test_fork_bary_units.py` | **no** — upstream |
| 4 | `compute_weight` takes an inverse variance | `components.py` | `test_fork_weights.py` | **no** — upstream |
| 5 | adapter contract repair | — (ours) | `test_fork_adapter_contract.py` | **yes** — ours alone |

34 tests pass. `vendor_pyodine.py --verify` reports **145 files identical to
upstream and four changed**, which is the whole diff and is exactly the four
intended.

**Change 0 was not on the list and had to come first.** Prompt 1 established
that `np.float` (`lib/misc.py`, three sites) and `np.NaN`
(`fitters/lmfit_wrapper.py`, five) stop the code dead against NumPy 2, and
prompt 1 was forbidden from fixing them. Nothing else in this prompt could be
tested end to end until they were gone, so they are change 0. The diagnostic
shim `pyodine_smoke.py` carried for prompt 1 is deleted with them.

**Change 1: `iod_depth` is a column density, not a multiplier.** Upstream
applied `flux = iod_depth * (flux - 1) + 1` at two places — `eval` and
`clean_of_I2`, which must agree or the I2-free reconstruction contradicts the
model. Both now call one helper, `scale_iodine_depth`, which applies
`clip(flux, 0, None) ** depth`.

The test that matters is the end-to-end one: build a synthetic chunk whose
atlas has a black line core, evaluate `SimpleModel.eval` at
`iod_depth = 2.59`, and look at the model flux. Before the change it is
**−1.012**; after, it cannot go below zero by construction. That is phase 2's
prediction reproduced inside the fork rather than argued from the source.

Three further properties are pinned because phase 3 depends on them: unit depth
is exactly the identity; depths compose multiplicatively, so the adapter may
pre-scale the atlas by `alpha = 2.59` and still let the fit vary `iod_depth`
around 1 with the two composing exactly; and for weak lines the new form agrees
with the old to a relative 1e-4, so nothing that was fitting sensibly moves.

**What the change does to a real fit.** The SONG tutorial, same template, only
the depth model differing:

| | Beer-Lambert | upstream linear |
|---|---|---|
| robust per-chunk scatter, epoch 1 | **159.4 m/s** | 163.7 m/s |
| robust per-chunk scatter, epoch 2 | **159.1 m/s** | 165.7 m/s |
| median chunk velocity | 1169.6, 1169.7 m/s | 1220.4, 1218.4 m/s |
| median fitted `iod_depth` | 0.988, 0.985 | 0.975, 0.979 |

So the change is a small improvement on data that did not need it — 3% in
scatter, and the two epochs agree with each other to 0.3 m/s instead of 2.0 —
and it moves the velocity zero point by about 50 m/s, common to both epochs and
therefore largely differential-cancelling. On HIRES, where the cell is 2.6x
thicker, the old form was not merely suboptimal but unphysical.

**A statistic to distrust.** The *standard deviation* of the per-chunk
velocities moved 911 → 1483 m/s across this change, which looks like a
catastrophe and is not one. A handful of chunks fail outright and land
thousands of m/s away — the worst is −27 km/s — and they dominate the standard
deviation completely: the same two epochs give 1483 and 961 m/s of standard
deviation while agreeing on 159.4 and 159.1 m/s of robust scatter. Prompts 5
and 7 are asked for "per-chunk velocity scatter"; it must be the robust one, or
the number will describe the failures rather than the fit.
`first_hires_exoplanet/pyodine_chunk_stats.py` now reports both, with outlier
counts, for exactly this reason.

**Change 2: an atlas must say which grid it is on.** The fix is *not* to read
the vacuum grid in the fork — Lick and SONG report air wavelengths and reading
vacuum would be as wrong for them as air is for us. What was missing is that
the convention was never recorded, so a mismatch could not be detected. So:
`IodineAtlas` gains `wave_frame` ('air', 'vacuum', or None), a constructor
`from_h5(filename, wave_frame)` **with no default** — choosing for the caller
is the mistake being prevented — and `require_wave_frame`, which raises
`DataMismatchError` on a mismatch and warns when an atlas does not know. The
HIRES adapter declares `vacuum`, and the test confirms on the real file that
the two grids differ by 1.56 A at 5615 A, which is 83 km/s.

**Change 3: the two comments, plus the one check a machine can make.** They now
read `full Julian Date in UTC (JD(UTC), not BJD)` and `m/s`, which is what
`timeseries/bary_vel_corr.py` and `chunks.py` actually require. A comment is a
weak thing to test, so `Observation.check_bary_units()` was added: a reduced
Julian Date is three orders of magnitude from a full one and is caught, as is a
barycentric velocity beyond 40 km/s. A correction given in km/s is *not*
mechanically distinguishable from a small one in m/s — which is precisely why
the documentation had to be right.

**Change 4: the weights are the inverse variance.** `Spectrum` now carries an
optional `ivar`, `__getitem__` slices it — a chunk is built as
`observation[order][pixels]`, so dropping it there would lose it exactly where
the fit needs it — and `compute_weight` gains `weight_type='ivar'`. Masked,
negative and non-finite entries return zero weight rather than passing through
the `sqrt(abs(weight))` in `lmfit_wrapper.py:96` as large ones. Asking for
`'ivar'` when none was supplied raises rather than silently returning ones.
`'flat'` and `'inverse'` are untouched.

This is the change that should make reduced chi-square mean something. With
flat weights the residual is in flux units and prompt 1 measured a median
reduced chi-square of 150516; with a propagated inverse variance the residual
is a proper chi. **That is not yet demonstrated on HIRES** — the SONG tutorial
has no inverse variance to feed it, so the number above is unchanged there, and
prompt 5 is where it gets tested on our own data.

**Change 5 was not asked for and could not be avoided.** Putting the real fork
on the path — phase 2's adapter had been subclassing stand-ins — immediately
raised `AttributeError: property 'orders' has no setter`.
`components.MultiOrderSpectrum` defines `orders` as a read-only property
giving *positions* (0..nord-1), which is what `pyodine`'s loops index with,
while phase 2's adapter had assigned the echelle order *numbers* (57-93) to the
same name. The physical numbers now live under `ech_orders`; `orders` is
inherited and means what the fork means by it. Each order also carries its
inverse variance down to the `Spectrum`, so a chunk sliced out of it is
weightable on its own.

The adapter self-check passes unchanged against the real fork: full JD,
bary_vel_corr −12677 to +11581 m/s, 150 of 1329 order-spectra zero-weighted,
every order still loadable, atlas on vacuum at alpha = 2.59.

**Which of these should go back to the author.** Five of the six, which is more
than expected:

* **NumPy 2 (change 0)** — urgent and unconditional. `pyodine` does not run at
  all against any NumPy from 2.0 onward, and fails silently because the
  pipeline swallows its own exceptions.
* **Beer-Lambert depth (change 1)** — generally correct physics, backwards
  compatible near depth 1, and a measured small improvement on the author's own
  SONG data.
* **The recorded wavelength frame (change 2)** — the mechanism, not the choice.
  Upstream's four adapters can keep reading air; what they gain is that a
  mismatch becomes an exception instead of 83 km/s of nonsense.
* **The `bary_*` documentation and guard (change 3)** — a pure documentation
  bug in upstream, and the guard costs nothing.
* **`compute_weight(weight_type='ivar')` (change 4)** — purely additive.
  Instruments with no propagated variance leave `ivar` at None and nothing
  changes for them.

Only the *choice* of the vacuum grid, and change 5, are ours and stay ours.

**On "each as a separate commit".** The work is in six separable units with
independent tests, but three of them land in `pyodine/components.py` and three
touch `utilities_hires/load_pyodine.py`, so they cannot be separated with
`git add <file>` alone. The hunks are contiguous and disjoint, so `git add -p`
splits them cleanly; the map is in the log entry below. Git remains the user's:
nothing here was staged or committed.

### Prompt 3: the template deconvolves, but not at `pyodine`'s default sampling

**The three answers, in brief.**

1. **S/N cost:** real per pixel, and none in velocity. Deconvolution cuts
   the per-pixel S/N from 793 to 564 at `pyodine`'s default oversampling, and
   to 393 at oversampling 1. But the forward model re-convolves the template.
   Once re-convolved, the template carries the same velocity noise per chunk
   as the observation it came from, 11–12 m/s. That is about 0.5 m/s per
   epoch, negligible against the 50 m/s target.
2. **Dependence on the assumed LSF:** large, and in one direction. The
   phase-2 resolution of R = 42,000 is not the instrument's. A template
   deconvolved with it misfits the data at χ² ≈ 264 per pixel once blurred
   by the true profile. A ±20% error in the width costs χ² 17–35 and moves
   chunk velocities by ±15–40 m/s. The data cannot choose the width
   themselves: every narrower kernel fits better, down to fitting the noise.
3. **Does it reproduce the observation when re-convolved? Not at the default
   settings.** χ² is 5.2 per pixel against the measured noise, with 390 of 700
   chunks above 4. **At an oversampling of 1 it does:** χ² 1.1, with 19 of 700
   above 4. The cause is `pyodine`'s spline, not the LSF.

`first_hires_exoplanet/deconvolve_template.py` produces every number below.
Nothing in `vendor/pyodine/` was changed. The deconvolution is
`ChunkedDeconvolver` → `jansson` with the LSF passed through its own
`lsf_fixed`. Chunking (40 px, 12 px padding) and the Jansson settings are
Lick's `Template_Parameters`, copied so prompt 4 cannot move them underneath
this result.

**The instrumental profile is R ≈ 66,000, not 42,000.** Gaussian fits to
1780 isolated ThAr lines from the 1998-08-26 reduction give:

- FWHM **2.18 px (4.55 km/s), R = 65,900**. PypeIt's own per-order estimate
  agrees (2.21 px).
- In the 14 iodine orders (58–71), 2.15–2.59 px, widening toward the red end.

That is the slit-limited resolution of the 0.574″ B1 decker. Phase 2's "45,000
nominal" belongs to the 0.861″ slit. Phase 2's R ≈ 42,000 came from the
narrowest weak *stellar* lines (7.45 km/s), so it is an **upper bound on the
instrument's width**. Taken in quadrature, the star contributes about
5.9 km/s of its own. That is an inference, but a plausible one for a
slowly rotating G dwarf.

The arc fills the slit, so a star that under-fills it in good seeing is
narrower still. The arc width is therefore itself an upper bound, only a much
tighter one. The measured stellar LSF belongs to prompt 6, and the data for
it exist on the template night itself:
- B stars through the iodine cell bracket the three template exposures
  (09:39–09:59 UT): HR 7236 at 08:21 UT and HR 8634 at 10:20 UT.
- There is an iodine-in quartz flat at 15:32 UT.
- **None of these is downloaded yet.**

**The propagated noise is 2.35 times too large, and it hides the failure.**
All three template exposures share one arc and one wavelength solution. They
can therefore be differenced pixel by pixel, with no resampling:

| exposure pair | orders | width of difference / propagated σ | detector shift |
|---|---|---|---|
| 0 − 1 | 14 | **0.426** [0.400, 0.438] | +44.5 m/s |
| 1 − 2 | 11 | 1.03 [0.76, 1.59] | −33.0 m/s |
| 0 − 2 | 11 | 1.05 [0.84, 1.49] | +6.6 m/s |

The clean pair (0 − 1) is consistent in every order:

| | fractional scatter of the ratio a/b |
|---|---|
| observed | 0.30% |
| photon noise at the header gain of 4.8 e⁻/ADU | 0.24% |
| photon noise at PypeIt's 1.9 e⁻/ADU | 0.38% |
| PypeIt's inverse variance | 0.72% |

The gain explains 1.59× of the 2.35×. The rest is undiagnosed. It is not the
`noise_floor` (0.0 in this reduction), nor the `adderr` of the local sky
subtraction (the extraction passes `noise_floor` there too).

The co-added template's true S/N in the iodine band is therefore **≈ 790 per
pixel**, not the 320–346 phase 2 quoted from the inverse variance. That is a
second answer to the gain question the phase document carried open.

Exposure 2 disagrees with the other two by more than noise, but not by a
shift. The detector shifts close to within 5 m/s and amount to 0.01–0.02 px
beyond the barycentric motion. Exposure 2 has 22% fewer counts and no usable
pixels in orders 58–60, which fits poorer seeing and a different line profile.
So the co-add mixes two slightly different LSFs. That is an inference, but a
measured difference.

**The trap:** against the propagated noise, the default template scores χ²
= **0.94** and looks perfect. Only against the measured noise does it score
5.2. Any χ² quoted from PypeIt's inverse variance in this project is too
small by about 5.5×.

**Why the default template fails: the spline, proved three ways.**
`ChunkedDeconvolver` puts the observation on its 10×-oversampled grid with
`misc.rebin`, a cubic spline through the detector pixels. At 2.2 px FWHM the
lines are barely Nyquist-sampled. Between pixel centres near a deep line, the
spline is wrong by 6–18σ at S/N 790. Jansson then fits those wiggles, and no
blurred spectrum can make them.

1. **On real data, χ² against line depth and oversampling** (150 chunks,
   orders 58/63/69):

   | oversampling | deepest line < 0.6 | 0.6–0.9 | > 0.9 | all |
   |---|---|---|---|---|
   | **1** | **0.7** | **1.2** | **0.9** | **1.1** |
   | 2 | 14.7 | 3.0 | 1.1 | 3.8 |
   | 4 | 21.7 | 3.7 | 1.1 | 5.1 |
   | 10 (Lick) | 23.8 | 3.9 | 1.2 | 5.5 |

2. **Synthetic, exact kernel, known truth.** Noise added on the fine grid:
   χ² 0.12–0.16 and bias 0.1–0.2σ, a perfect recovery. The same truth
   sampled at the pixels and splined as `pyodine` does: χ² 2–95 and bias
   1.3–5.6σ, scaling with line depth.
3. **Neither the solver nor the LSF moves it.** On the worst chunks, 5000
   iterations, a 5× larger step, or a continuum ceiling of 1.10 change χ² by
   less than 5%. Halving the kernel width leaves χ² at 30–50.

At oversampling 1, `rebin` is the identity. The same deconvolver then
reproduces the data *and* conserves equivalent width (1.000, against 0.982 at
10). **Both templates are written.** Whether the forward model can use the
1-px one is prompt 5's question: `pyodine` interpolates the template again
when it evaluates the model, and could reintroduce the same error there. The
lasting fix is probably in the fork, making the deconvolver fit only at
detector pixels. That is a fork change, and it is left for a decision.

**The S/N cost, by Monte Carlo** (20 realisations of the measured noise, 150
chunks, ThAr LSF; medians):

| | oversampling 10 | oversampling 1 |
|---|---|---|
| S/N per pixel, observed | 793 | 793 |
| S/N, deconvolved (per detector pixel) | 564 (×0.71) | 393 (×0.49) |
| S/N, deconvolved then re-convolved | 1174 (×1.48) | 1066 (×1.34) |
| velocity noise per chunk, observed | 12.5 m/s | 12.5 m/s |
| — deconvolved | 14.8 m/s (×1.10) | 16.6 m/s (×1.29) |
| — re-convolved (what the model uses) | 11.3 m/s (×0.91) | 12.4 m/s (×1.00) |

Deconvolution buys resolution with per-pixel noise, and re-convolution hands
the noise back. In velocity it is close to a wash. The default setting looks
cheaper only because the spline floor stops it deconvolving the deep lines.

**The dependence on the assumed LSF** (150 chunks; χ² against the measured
noise; velocities relative to the ThAr-LSF template; median [16th, 84th]):

| template built with | χ², own LSF (os 1) | χ², blurred by ThAr LSF (os 1) | Δv (os 1) | χ², blurred by ThAr LSF (os 10) |
|---|---|---|---|---|
| ThAr FWHM | 1.1 | 1.1 | — | 5.5 |
| 0.8 × ThAr | 0.2 | 16.6 | 0 [−15, +14] m/s | 31 |
| 1.2 × ThAr | 1.9 | 34.6 | +3 [−21, +40] m/s | 13 |
| R = 42,000 | **7.8** | **264** | +12 [−62, +83] m/s | 176 |
| ThAr, skewed (25% red satellite) | 1.5 | 10.0 | **+61** [+29, +96] m/s | 2.7 |

Three things follow:
- **R = 42,000 cannot be right.** It fails even against itself: undoing that
  much blur would need negative flux.
- **"Reproduces the data with its own LSF" cannot choose the width.** Every
  narrower kernel does better, down to χ² 0.2, which is fitting noise. The
  stellar LSF has to come from an independent measurement, the B stars, in
  prompt 6.
- **An asymmetric LSF moves the template's velocity zero point by ~60 m/s.**
  That is mostly common to all chunks, so it cancels between epochs if the
  same LSF is used throughout, but it scatters by ±30 m/s from chunk to chunk.

**`jansson` uses the wrong adjoint, and it barely matters here.** The update
is `x += r(K y − K K x)` where least squares needs `Kᵀ`. The fixed point is
the same; the path to it is not.
- **Symmetric kernel:** the fork and a true-adjoint copy agree to machine
  precision. This is checked in the run and in a test.
- **Skewed kernel:** the fork is biased by −5.7 [−24.4, +0.3] m/s per chunk
  relative to the true adjoint, with identical χ² and RMS.

This is small against the target, but a correctness point for upstream. It
will matter more if prompt 6 chooses an asymmetric LSF model.

**Six other things this prompt tripped over.**

1. **Masked pixels have a wavelength of exactly zero, not just a flux of
   zero.** The template has zeros at the first and last pixel of every order,
   and at 878 pixels of order 57. `pyodine`'s `get_velocity_offset` divides by
   `wave[0]` and `wave[-1]`. The loader here extrapolates a polynomial over
   those pixels and keeps them at zero weight. **The phase-2 adapter has the
   same hole:** a real spec1d epoch reaches `pyodine` with 950 zero-wavelength
   pixels across all 37 orders, and its self-check never noticed. It will
   break prompt 5 at the velocity-guess step, so prompt 4 should fix it.
2. **`pyodine`'s reference atlas is on air wavelengths.** The Arcturus
   cross-correlation returns **+74.9 km/s** for a star whose observed
   velocity is about −10 km/s. The difference is the ~84 km/s air-to-vacuum
   offset in the blue orders used for the guess. It is harmless as a
   continuum aid and in differences, since template and observation carry
   it alike, but the value must never be read as a velocity.
3. **`NormalizedObservation` drops the inverse variance.** Its `__getitem__`
   rebuilds each `Spectrum` without `ivar`, undoing prompt 2's change at the
   first normalisation. Here the deconvolver does not weight anyway, so the
   noise is carried alongside. It will matter wherever a normalised
   observation is fitted.
4. **The continuum mostly does not come from the solar comparison.** It came
   from the atlas comparison in only 1 of 14 iodine orders; the other 13 fell
   back to `top`, a quadratic upper envelope. `normalize_single` swallows the
   `ValueError`, so nothing says so.
5. **`LinearWaveModel` misplaces PypeIt's wavelengths.** Its straight line
   through PypeIt's polynomial is off by up to 64 m/s (median of the
   per-chunk maxima; 162 m/s worst) across a padded chunk. That goes into the
   template's wavelength scale near chunk edges.
6. **Jansson's stopping rule is not a fit-quality test.** It compares the
   observation with the *deconvolved* estimate, not its re-convolution. Most
   chunks run all 1202 iterations; some stop early (one at 394). The
   iteration count, not the rule, sets how far it deconvolves.

**Products.**
- **Templates**, read back identically by `StellarTemplate_Chunked`, with no
  negative flux in either:
  `redux/template/hd187123_template_19980826_deconv_os10.h5` (Lick settings,
  700 chunks, 23 s) and `…_deconv_os1.h5` (700 chunks, 6 s).
- **Tables** in `first_hires_exoplanet/data/`: `deconv_arc_widths.csv`,
  `deconv_noise_check.csv`, `deconv_lsf_cases.csv`, `deconv_osample.csv`,
  `deconv_snr_mc.csv`, `deconv_adjoint.csv`, `deconv_chunks.csv`.
- **Figure:** `docs/figs/fig_p3_deconv.png`.

### Prompt 4: the machinery fits a chunk, and the model does not fit the data

**The Q&A decisions, applied.**

| question | decision | what was done |
|---|---|---|
| 1. which template | (a) both | the chunk is fitted against each |
| 2. adapter zero wavelengths | fix now | done, and the same bug in the *flux* too (below) |
| 3. B stars | 1998-08-26, before prompt 5 | 12 frames downloaded; not reduced |
| 4. noise scale | (b) weights untouched | χ² reported both ways; errors checked by Monte Carlo |

**The four `utilities_hires` files.** They are modelled on `utilities_lick`
and read by the same drivers through the same attributes. Every departure
from Lick is commented with its source.

- **`conf.py`:** a Keck `Instrument` and one atlas (index 1, the Fischer
  cell). Unlike Lick, it also records the wavelength grid the atlas must be
  read on (`vacuum`). `IodineTemplate(1)` refuses anything else.
- **`pyodine_parameters.py`:**
  - `weight_type = 'ivar'`.
  - `noise_scale = 0.426`, for reporting only.
  - `wavelength_scale = 'vacuum'`, `i2_to_use = 1`.
  - `velgues_order_range = (15, 37)`: orders 72–93. Lick's (2, 15) would sit
    inside the iodine band here.
  - `maxlag = 500`, `number_cores = 4`.
  - Run 0 is a single Gaussian started at the ThAr 2.2 px; run 1 is Lick's
    multi-Gaussian, pending prompt 6.
  - `Template_Parameters` reproduces prompt 3's deconvolution exactly (700
    chunks over orders 58–71).
- **`timeseries_parameters.py`:** Lick's settings, except that the
  barycentric correction comes from coordinates, not a Hipparcos number parsed
  from the star name. The 1998 headers name the star three different ways.
- **`logging.json`:** no file paths. Lick's hard-codes `/home/paul/…`.
  `setup_logging` deletes both file handlers unless a driver passes a path,
  so the file only needs the console.

The package `__init__` now exposes `load_pyodine`, `conf`,
`pyodine_parameters` and `timeseries_parameters`, as `utilities_lick` does.

**The adapter: three fixes, six tests, each test shown failing first.**

1. **Zero wavelengths**, as decided. Every order of every spec1d has zeros at
   masked pixels. They are now filled by a polynomial at zero weight.
   `fill_wavelengths` moved here from `deconvolve_template.py`, which imports
   it.
2. **Zero-flux runs. This was not in the decision, but the fit could not run
   without it.** PypeIt zero-fluxes masked runs *inside* orders too:
   1998-08-26 19326 has 113 consecutive zeros in order 70 (pixels 954–1066).
   One 40-pixel chunk inside such a run makes `Spectrum` raise
   `NoDataError` for the whole epoch, because `auto_wave_comoving` builds
   every chunk before any is fitted. The phase-2 guard (never zero the flux
   of a rejected order) did not cover PypeIt's own zeros. They are now
   interpolated, at zero weight.
3. **`IodineTemplate` takes the integer atlas index** the drivers pass.

The adapter self-check still passes. The test suite is at 43.

**The B stars.** Twelve frames were downloaded, all FITS-verified, 57 MB:
- 11 B-star exposures: HR 7236 ×2 (08:21 UT), HR 8634 ×6 (10:20, 11:23,
  13:16 UT) and HR 838 ×3 (15:25 UT);
- the iodine-in quartz flat (15:32 UT).

They sit in `raw/1998aug26_lsf/`, *beside* the night's directory, because
`reduce_run.frame_roles` globs the night directory and would otherwise reduce
them as HD 187123 science frames. `koa_download.py` gained the two selections.
They are not reduced; that belongs to prompt 6.

**The chunk.**
- **Epoch:** `HI.19980826.19326` (cell in, 400 s), the template's own night.
- **Order:** 63, mid-band.
- **Chunk:** 276, pixels 1064–1103, 5661.45–5663.00 Å. The rule was stated
  before fitting: the chunk nearest the order's middle with no masked pixel.
  The literal middle chunk, 275, has 14 of its 40 pixels at zero weight, a
  PypeIt-masked run near detector column 1050.
- **How it is built:** everything up to the fitting loop is
  `model_single_observation`'s own code path, and all 700 chunks are built.
  Only this one is fitted, with run 0's model: single-Gaussian LSF, linear
  wavelength, linear continuum.

**Every fitted parameter.** The "±" column is `lmfit`'s standard error. It is
already rescaled by √χ²_ν, so it absorbs both the 2.35× noise mis-scale and
the model's misfit. The "noise-only" column is the scatter of 30 refits of the
same chunk, started at the best fit, on the best-fit model plus noise at the
*measured* level.

| parameter | os 10 | ± | noise-only | os 1 | ± | noise-only |
|---|---|---|---|---|---|---|
| velocity (m/s) | −417.5 | 267.7 | 27.4 | −323.1 | 240.3 | 25.7 |
| tem_depth | 0.961 | 0.054 | 0.007 | 0.946 | 0.049 | 0.008 |
| iod_depth | 1.121 | 0.077 | 0.008 | 1.118 | 0.071 | 0.008 |
| lsf_fwhm (px) | 2.080 | 0.113 | 0.013 | 2.067 | 0.103 | 0.013 |
| wave_intercept (Å) | 5662.22920 | 0.00207 | 0.00009 | 5662.22970 | 0.00191 | 0.00009 |
| wave_slope (Å/px) | 0.0393331 | 0.0000446 | 0.0000063 | 0.0393243 | 0.0000502 | 0.0000087 |
| cont_intercept (counts) | 72236 | 550 | 51 | 72224 | 486 | 46 |
| cont_slope (counts/px) | 0.18 | 18.7 | 3.1 | 10.4 | 17.3 | 3.1 |
| χ²_ν, PypeIt σ | 8.13 | | | 7.18 | | |
| **χ²_ν, measured σ** | **44.8** | | | **39.6** | | |

**The checks against independent quantities.**

- **velocity: −417 m/s (os 10) against −410 expected**, with os 1 at −323.
  The expectation is the barycentric difference (−432.8 m/s) plus the
  planet's motion between the two times. The planet term is +23.3 m/s, from
  a circular fit to the modern catalogue that returns K = 68.3 m/s.
  Agreement to 8 m/s is within the noise-only 27 m/s, and os 1's 87 m/s gap
  is 3.4 of it. One chunk cannot choose between the templates.
- **lsf_fwhm: 2.07–2.08 px against the ThAr 2.17 px** of order 63. That is 4%
  narrower, the direction prompt 3 predicted for a star that under-fills the
  slit, though not significant against the misfit-scaled error.
- **iod_depth: 1.12.** The adapter's α = 2.59 is 12% short for this chunk,
  giving an effective exponent of 2.90, just above phase 2's 2.28–2.86 range
  from two 60 s frames.
- **tem_depth: 0.95–0.96.** The template's lines are 4–5% too deep for the
  observation.
- **The wavelength solution is −566 m/s (os 10) / −539 m/s (os 1) from
  PypeIt's ThAr solution**, with dispersion equal to within 0.01–0.03%. The
  noise-only error on the intercept is 5 m/s.

  That is the size of zero-point error phase 2 blamed for its 702 m/s. But
  one chunk cannot tell a PypeIt error, which should vary from epoch to epoch,
  from an offset of the Fischer atlas, which would be constant. Prompt 5
  compares across a whole epoch.
- **Continuum:** 72 236 counts against a chunk median of 70 874, as expected
  for a chunk whose median sits below the continuum.

**What the fit says about the model: it does not fit the data.**
χ²_ν = 45 against the measured noise. The residuals are not noise
(`docs/figs/fig_p4_one_chunk.png`): ±10–15σ, with the data *deeper* than the
model at every minimum and *higher* at every maximum. The observation has
more contrast than the model can make. Both templates show it (45 and 40), so
the template's sampling is not the main cause. Three candidates remain, and
none can be separated with one chunk:
- the atlas is a different cell, and phase 2 measured its line pattern at
  only 0.928 correlation with the HIRES cell's;
- a single Gaussian is not the HIRES profile;
- the template is under-deconvolved.

The first two are prompts 5 and 6.

**What the uncertainties mean.**

- **The error machinery is internally consistent.** Refits from the best fit
  land at χ²_ν = 0.20–0.21 against PypeIt's σ, where 0.426² = 0.18 is
  expected. Their parameter scatter matches the unscaled curvature errors.
- **So the reported ±268 m/s is about 90% model misfit.** The noise alone
  would give ±27 m/s for this chunk. With decision (b), `scale_covar` makes
  the error bars honest about the misfit, not about the photons. Prompt 7's
  velocity errors will mean "how badly the model fits", unless the model
  improves.
- **From `pyodine`'s own starting guess, only 30–37% of refits reach the
  minimum.** The rest stop at χ²_ν 0.56–2.1 against the expected 0.18. The
  start is 640 m/s from the answer in velocity and 510 m/s in wavelength
  intercept, and the two are strongly degenerate within one chunk. Run 0
  exists to get close before the cross-order smoothing, but prompt 5 should
  expect a large fraction of run-0 chunks not to have converged.

**Four things in `pyodine` this prompt tripped over.**

1. **The Chauvenet re-fit is a no-op.** `pipe_lib.model_all_chunks` re-fits
   only `if any(mask)==False`, i.e. when *every* pixel fails the criterion.
   It evidently means "if any pixel fails". Worse, the criterion runs on
   unweighted residuals, so zero-weight pixels count: on the literal middle
   chunk, 5 of its 6 flags were masked pixels. `use_chauvenet_pixels = True`
   therefore does nothing, and if repaired it would reject the wrong pixels.
2. **`red_chi_sq` is not χ²_ν.** `model_all_chunks` stores its square root
   under that name, so anything reading the results must square it first.
3. **The iodine model divides the atlas by its mean before exponentiating.**
   (T/⟨T⟩)^d = T^d·⟨T⟩^(−d): still Beer–Lambert, but the constant is absorbed
   by the continuum. `iod_depth` and `cont_intercept` are correlated by
   construction.
4. **The velocity guess is 640 m/s from the answer** (+0.226 km/s against
   −0.41). That is fine for a guess at 1 km/s steps, but it is where the
   convergence problem above starts.

**Products.**
- **Code:**
  - `first_hires_exoplanet/utilities_hires/`: `conf.py`,
    `pyodine_parameters.py`, `timeseries_parameters.py`, `logging.json`;
    `load_pyodine.py` and `__init__.py` modified;
  - `first_hires_exoplanet/fit_one_chunk.py`.
- **Tables:** `first_hires_exoplanet/data/one_chunk_fit.csv` and
  `one_chunk_mc.csv`.
- **Figure:** `docs/figs/fig_p4_one_chunk.png`.
- **Tests:** `first_hires_exoplanet/tests/test_adapter_wavelengths.py` (six).
- **Raw data:** `../first-hires-exoplanet-data/raw/1998aug26_lsf/` (12
  frames).

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

### 2026-09-23 (Prompt 2: the four HIRES changes to the fork, plus the NumPy 2 fix they all depended on)

**Task.** Make the four phase-2 changes to the vendored fork — Beer-Lambert
depth, vacuum atlas grid, the two `components.py` docstrings, and an inverse
variance in `compute_weight` — each as a separate commit with a test that fails
before and passes after, and report which are HIRES-specific. Findings in
`## Report`, "Prompt 2: five changes to the fork, five tests, and one that was
not asked for".

**What was done.** Six units, each test written and run against the unchanged
tree first and shown to fail for the right reason before the change was made.
In the fork: `lib/misc.py` and `fitters/lmfit_wrapper.py` for NumPy 2;
`models/spectrum.py` gained `scale_iodine_depth` and both iodine sites now call
it; `components.py` gained `IodineAtlas.wave_frame` / `from_h5` /
`require_wave_frame`, corrected `bary_date` and `bary_vel_corr` with a new
`Observation.check_bary_units`, and an optional `Spectrum.ivar` with
`compute_weight(weight_type='ivar')`. In the repository: six test files and a
`conftest.py` under `first_hires_exoplanet/tests/`, a new
`pyodine_chunk_stats.py`, the removal of `pyodine_smoke.py`'s NumPy shim, and
repairs to `utilities_hires/load_pyodine.py` (fork on the path, `ech_orders`,
per-order `ivar`, `'ivar'` accepted alongside `'inverse-variance'`). Ran the
SONG tutorial end to end after the changes, and again with the depth model
reverted, to measure what the change actually did.

**Headline results.** 34 tests pass; `--verify` shows four changed files and
145 identical to upstream. The end-to-end test of the depth model produced
**−1.012** flux from `SimpleModel.eval` at `iod_depth = 2.59` before the change
and cannot go negative after. On the SONG tutorial, same template, the new
depth model gives a robust per-chunk scatter of 159.4 and 159.1 m/s against
163.7 and 165.7 for upstream's linear form, and shifts the velocity zero point
by about 50 m/s. The adapter self-check passes against the real fork with
phase 2's numbers unchanged. Five of the six changes are not HIRES-specific and
should go upstream.

**What this taught us about the repository and the data.**

*The obvious scatter statistic is the wrong one, and it nearly produced a false
alarm.* The per-chunk velocity standard deviation moved 911 → 1483 m/s across
the depth change and looked like a 63% regression. It is not: a few chunks fail
completely — the worst lands at −27 km/s — and the standard deviation is
theirs, not the fit's. The same two epochs give 1483 and 961 m/s of standard
deviation while agreeing to 0.3 m/s on the robust width. Prompts 5 and 7 ask
for per-chunk scatter by name; reported as a standard deviation the number will
describe the failures. `pyodine_chunk_stats.py` exists so that both are always
printed together with an outlier count.

*Phase 2's adapter had a latent contract violation that only the real fork
could reveal.* `orders` means positions to `pyodine` and echelle numbers to
phase 2's adapter, and because phase 2 subclassed stand-ins, nothing ever
compared the two. It raised the moment the fork was on the path — loudly, which
was luck: had `MultiOrderSpectrum.orders` been a plain attribute rather than a
property, the assignment would have succeeded and every `pyodine` loop over
`obs.orders` would have indexed with 57-93 into a 37-order array. The general
lesson for prompt 4, which writes four more `utilities_hires` files: a
"drop-in" module that has never been dropped in is a claim, not a fact.

*Fixing `pyodine` for HIRES mostly means fixing `pyodine`.* Only the choice of
the vacuum grid turned out to be genuinely ours. The depth model, the NumPy 2
breakage, the mis-documented barycentric units and the absence of any way to
pass a propagated variance are all defects against any instrument; HIRES merely
made them visible, because its cell is thick, its pipeline is modern and its
data is old enough to need a modern NumPy. That is worth stating plainly in the
prompt-10 report: we are not asking the author to accommodate a special case.

*Reduced chi-square is still not usable, and will not be until prompt 5.* The
weights change is the one that should fix it, but the SONG tutorial has no
propagated variance to feed it, so 150516 is still what the tutorial reports.
The first honest chi-square in this project will come from HIRES data, and if
it is not of order unity the weights are still wrong.

*A depth model is not a free parameter change.* Making `iod_depth` an exponent
moved SONG's velocity zero point by 50 m/s. It is common to both epochs and so
mostly cancels differentially, but it is a reminder that the fork's changes
have to be finished before any velocity is quoted, not adjusted afterwards.

**Files added / modified.** In the fork: `vendor/pyodine/pyodine/lib/misc.py`,
`vendor/pyodine/pyodine/fitters/lmfit_wrapper.py`,
`vendor/pyodine/pyodine/models/spectrum.py`,
`vendor/pyodine/pyodine/components.py`. Added
`first_hires_exoplanet/tests/conftest.py`,
`first_hires_exoplanet/tests/test_fork_numpy2.py`,
`test_fork_iodine_depth.py`, `test_fork_atlas_frame.py`,
`test_fork_bary_units.py`, `test_fork_weights.py`,
`test_fork_adapter_contract.py`, and
`first_hires_exoplanet/pyodine_chunk_stats.py`. Modified
`first_hires_exoplanet/pyodine_smoke.py` (shim removed) and
`first_hires_exoplanet/utilities_hires/load_pyodine.py`.

Suggested commit split, since three changes share `components.py` and three
share `load_pyodine.py` (`git add -p` separates them; the hunks are disjoint):

1. *NumPy 2*: `lib/misc.py`, `fitters/lmfit_wrapper.py`, `pyodine_smoke.py`,
   `tests/conftest.py`, `tests/test_fork_numpy2.py`.
2. *Beer-Lambert depth*: `models/spectrum.py`, `tests/test_fork_iodine_depth.py`,
   `pyodine_chunk_stats.py`.
3. *Atlas wavelength frame*: `components.py` (the `IodineAtlas` block,
   lines ~314-400), `load_pyodine.py` (the fork-on-path and `IodineTemplate`
   hunks), `tests/test_fork_atlas_frame.py`.
4. *Barycentric units*: `components.py` (the two declarations and
   `check_bary_units`), `tests/test_fork_bary_units.py`.
5. *Inverse-variance weights*: `components.py` (`Spectrum.__init__`,
   `__getitem__`, `compute_weight`), `load_pyodine.py` (the `'ivar'` alias and
   per-order `ivar`), `tests/test_fork_weights.py`.
6. *Adapter contract*: `load_pyodine.py` (`ech_orders`),
   `tests/test_fork_adapter_contract.py`.

**Whether any git command changed state.** No. Read-only `git ls-files` and
`git rev-parse` inside `vendor_pyodine.py --verify`, which clones upstream into
a scratch directory outside the repository. Nothing was staged or committed
here; the commit split above is a suggestion for the user to apply.

### 2026-09-28 (Prompt 3: template deconvolved with `pyodine`'s own deconvolver; it reproduces the data only without the spline)

**Task.** Deconvolve the stellar template against the instrumental profile with
`pyodine`'s `template/deconvolve.py`. Report the S/N cost, the dependence on
the assumed LSF, and whether the template reproduces the observation when
re-convolved. Findings in `## Report`, "Prompt 3: the template deconvolves,
but not at `pyodine`'s default sampling".

**What was done.**

- Wrote `first_hires_exoplanet/deconvolve_template.py`, which:
  - measures the instrumental width from 1780 ThAr lines of the 1998-08-26
    reduction;
  - measures the true noise by differencing the three template exposures
    pixel by pixel;
  - presents phase 2's co-added template to `pyodine` as an Observation;
  - normalises it with `SimpleNormalizer`, chunks it with `auto_equal_width`
    (Lick settings), and takes per-chunk wavelength parameters from
    `SimpleModel.guess_params`;
  - deconvolves with `ChunkedDeconvolver` through `lsf_fixed`.
- Five experiments on 150 chunks of orders 58/63/69:
  1. five assumed LSFs, each at oversampling 10 and 1;
  2. an oversampling scan (1, 2, 4, 10);
  3. a 20-realisation Monte Carlo at oversampling 10 and 1;
  4. a synthetic adjoint test (fork `jansson` against a true-adjoint copy,
     symmetric and skewed kernels);
  5. full 700-chunk templates at oversampling 10 and 1.
- Three exploratory scripts in the scratchpad (iteration counts, χ² against
  line depth and solver settings, spline against kernel width) guided the
  design. Every number quoted in the Report comes from the module on disk.
- Added `first_hires_exoplanet/tests/test_deconvolve_template.py` (three
  tests). 37 tests pass.

**Headline results.**

- **Instrumental width:** ThAr gives FWHM 2.18 px (4.55 km/s), **R ≈ 65,900**.
  Phase 2's R = 42,000 was the stellar line width, an upper bound on the
  instrument's.
- **Noise:** PypeIt's propagated noise is **2.35× too large**. The clean
  exposure pair scatters at 0.426σ in all 14 iodine orders, so the template's
  true S/N is ≈ 790 per pixel, not ≈ 330.
- **Reproduction at Lick's oversampling of 10:** the deconvolved template does
  **not** reproduce the observation when re-convolved: χ² 5.2, 390 of 700
  chunks above 4, and 23.8 on chunks with lines deeper than 0.6.
- **Reproduction at oversampling 1:** it does (χ² 1.1, 19 of 700 above 4), and
  it conserves equivalent width.
- **S/N cost:** the per-pixel S/N falls to 0.71× (oversampling 10) or 0.49×
  (oversampling 1). The re-convolved template's velocity noise per chunk
  equals the observation's, 11–12 m/s.
- **LSF dependence:**
  - R = 42,000 is ruled out (χ² 7.8 against itself, 264 against the ThAr LSF).
  - A ±20% width error costs χ² 17–35 and ±15–40 m/s per chunk.
  - A skewed LSF shifts chunk velocities by ~60 m/s.
- **Adjoint:** `jansson`'s wrong adjoint biases a skewed-kernel deconvolution by
  −5.7 m/s per chunk (median), and does nothing for a symmetric kernel.

**What this taught us about the repository and the data.**

*`pyodine`'s default template sampling does not suit HIRES, and the reason
is the sampling, not the physics.* `ChunkedDeconvolver` splines the
observation onto a 10× grid with `misc.rebin`. For lines 2.2 pixels wide
that spline is wrong between pixel centres by up to 18σ at S/N 790, and
Jansson faithfully deconvolves the error. Three independent checks agree:
- the χ² floor tracks line depth;
- it vanishes on synthetic data when the spline step is removed;
- it vanishes on real data at oversampling 1.

More iterations, a larger step, a higher continuum ceiling or a narrower
kernel do not remove it. This is the most important thing prompt 3 found for
prompts 4–6. Whatever `pyodine` does to the template when it evaluates the
model may carry the same error, and a better-sampled instrument (the Lick
and SONG data the defaults were tuned on) would have hidden it.

*The inverse variance cannot be trusted to judge a fit.* Against it, the
broken default template scores χ² = 0.94 and looks perfect. Only the noise
measured from exposure differences exposes the misfit. The measurement needs
consecutive exposures on one wavelength solution, differenced **in pixel
space**. Resampling onto a common wavelength grid smooths the noise and would
bias the test.

The factor is 2.35. The header gain (4.8 against PypeIt's hard-coded 1.9)
explains 1.59× of it; the remaining 1.5× is not the `noise_floor` and not
the local sky subtraction's `adderr`, and is still unexplained. This
answers the phase document's open gain question with a measurement, and it
matters for prompt 5: reduced χ², the fitter's weights and any per-chunk
error bar all inherit it.

*Estimators that failed:*
- **DER_SNR** reported noise at three times the propagated σ, because it
  assumes a spectrum smooth over five pixels and HIRES lines are 2.2 pixels
  wide. It was measuring line curvature.
- **An LSF's own χ² cannot choose the LSF.** Every narrower kernel
  "reproduces the data" better, down to χ² 0.2. It can only rule out kernels
  that are too wide to be undone by non-negative flux, as R = 42,000 is.
- **Jansson's stopping rule** compares the data with the deconvolved
  estimate rather than its re-convolution, so it is not a convergence test.
  The iteration budget is the regulariser.

*The resolving power was wrong since phase 2, and in the direction that
hurts.* The stellar weak-line width (7.45 km/s) includes about 5.9 km/s of
the star's own broadening. The instrument is R ≈ 66,000, as the 0.574″ B1
slit implies; the "45,000 nominal" in phase 2 was the 0.861″ slit's.
Deconvolving with the phase-2 value would strip real line width from the
template and cost χ² ~260 against the data. The B stars that measure the
real stellar LSF are on the template night itself: HR 7236 and HR 8634,
bracketing the template exposures, plus an iodine-in quartz flat. None is
downloaded yet; that is prompt 6's input.

*Masked pixels carry a wavelength of exactly zero, not only a flux of zero.*
This is trap 4 again, one column over. `pyodine` divides by the first and
last wavelength of every order in its velocity guess, and every order of
every HIRES spec1d has a zero there. The template loader here extrapolates
over them. **The phase-2 adapter does not, and its self-check passed**
because it never called the velocity guess. Prompt 4 must fix
`utilities_hires/load_pyodine.py` before prompt 5 can run.

*Three more silent behaviours in the fork:*
- `NormalizedObservation` drops the `ivar` prompt 2 added.
- `normalize_single` swallows the failure of the solar-atlas continuum and
  falls back to an envelope fit, which happened in 13 of 14 iodine orders.
- The velocity guess against `pyodine`'s air-wavelength reference comes back
  as +74.9 km/s for a star at about −10 km/s. That is harmless in
  differences, but a trap for anyone who reads it as a velocity.

*The template night is not quite one LSF.* Exposures 0 and 1 agree to photon
noise. Exposure 2 differs from both by more than noise and not by a shift (it
has 22% less flux and three orders lost), so the co-add mixes two slightly
different profiles. That belongs to prompt 9's drift question.

**Files added / modified.**
- Added:
  - `first_hires_exoplanet/deconvolve_template.py`;
  - `first_hires_exoplanet/tests/test_deconvolve_template.py`;
  - in `first_hires_exoplanet/data/`: `deconv_arc_widths.csv`,
    `deconv_noise_check.csv`, `deconv_lsf_cases.csv`, `deconv_osample.csv`,
    `deconv_snr_mc.csv`, `deconv_adjoint.csv`, `deconv_chunks.csv`;
  - `docs/figs/fig_p3_deconv.png`.
- Outside the repository:
  `../first-hires-exoplanet-data/redux/template/hd187123_template_19980826_deconv_os10.h5`
  and `…_deconv_os1.h5`.
- Nothing in `vendor/pyodine/` or `utilities_hires/` was modified.
- Provenance: PypeIt 2.0.2.dev1217+g017bece06; `pyodine`
  4488b0914fe5b272b787982647691045bff2604a plus the prompt-2 changes.

**Whether any git command changed state.** No. Only read-only `git status`
was run. Nothing was staged or committed.

### 2026-09-29 (Prompt 4: `utilities_hires` completed, the adapter repaired, B stars fetched, one chunk fitted)

**Task.** Write `conf.py`, `pyodine_parameters.py`, `timeseries_parameters.py`
and `logging.json` modelled on `utilities_lick`. Fit one chunk of one order of
one epoch, and report every parameter with its uncertainty. First apply the
four Q&A decisions taken after prompt 3. Findings in `## Report`, "Prompt 4:
the machinery fits a chunk, and the model does not fit the data".

**What was done.**

- Wrote the four files, each departure from Lick commented with its source,
  and exposed them from the package `__init__`.
- Repaired `utilities_hires/load_pyodine.py`. Each test was written and shown
  failing before its fix:
  - zero wavelengths filled (decision 2);
  - zero-flux runs filled, at zero weight (found here: without it the epoch
    cannot be chunked);
  - `IodineTemplate` accepts the integer atlas index.
- Added `tests/test_adapter_wavelengths.py` (six tests). Moved
  `fill_wavelengths` out of `deconvolve_template.py` into the adapter.
- Added two 1998-08-26 selections to `koa_download.py` and downloaded the 11
  B-star exposures and the iodine-in flat into `raw/1998aug26_lsf/`, dry run
  first (decision 3).
- Wrote `first_hires_exoplanet/fit_one_chunk.py`. It follows
  `model_single_observation` to the fitting loop and fits one chunk (order 63,
  chunk 276, epoch 19326) against both prompt-3 templates (decision 1). It
  checks each parameter against an independent quantity, and runs two
  30-refit Monte Carlos of that chunk, one from the best fit and one from the
  driver's start.
- A scratch script confirmed the Monte Carlo's injected model is exact and
  that the Chauvenet flags fall on masked pixels.

**Headline results.**

- **The fit runs and converges.**
  - velocity −417 ± 268 m/s against −410 expected (os 10); os 1 gives −323;
  - lsf_fwhm 2.08 px against the ThAr 2.17;
  - iod_depth 1.12, tem_depth 0.96;
  - the fitted wavelength solution sits 0.54–0.57 km/s from PypeIt's.
- **But χ²_ν is 45 against the measured noise**, with residuals of ±10–15σ
  structured at every line extremum: the data have more contrast than the
  model.
- **The error machinery is internally consistent** (refits from the best fit
  give χ²_ν 0.20 against the expected 0.18). The noise-only velocity error is
  27 m/s, so ~90% of the reported ±268 m/s is misfit.
- **Only 30–37% of refits from `pyodine`'s own starting guess reach the
  minimum.**
- 43 tests pass, and the adapter self-check passes.

**What this taught us about the repository and the data.**

*A drop-in adapter has to survive the whole driver, not one call.* Prompt 3
found the zero wavelengths from the template side. This prompt found the
zero-flux runs from the chunking side: PypeIt zero-fluxes a 113-pixel run
inside order 70 of 19326, and one chunk inside it takes down the whole epoch,
because `auto_wave_comoving` builds every chunk before fitting any. Phase 2's
rule, "zero-weight, never zero-flux", was right but incomplete: it guarded
against *us* zeroing flux, not against PypeIt having done it. The new end-to-end
test (build every chunk of a real epoch against the real template) is the
one that would have caught both.

*The standard errors are honest about the misfit and silent about the
photons.* With the weights left as PypeIt made them (decision 4b),
`lmfit`'s `scale_covar` rescales every error by √χ²_ν. At χ²_ν = 45 that turns
a 27 m/s photon-limited chunk into a ±268 m/s one. That is defensible, since
the model really does not fit, but it means the velocity errors of prompts
5–7 will measure model inadequacy until the model improves. A per-chunk
scatter compared against these errors will look fine while being dominated
by the misfit. The Monte Carlo from the best fit is the check that separates
the two, and it costs 0.3 s per chunk.

*`pyodine`'s starting point is not inside the basin of its own minimum, for
one chunk.* Within 40 pixels the velocity and the wavelength zero point are
nearly degenerate. The iodine lines are what separate them, and they are
weak. From the driver's start (640 m/s and 510 m/s off respectively) two
refits in three stop short. Lick's recipe of cross-order smoothing after run
0 exists for this reason, but prompt 5 should count unconverged chunks and
not only failed ones. A converged-looking `success=True` is not convergence
here.

*The model's misfit has a shape, and it is the same for both templates.* The
data are deeper at every minimum and higher at every maximum. The template's
sampling therefore matters less than prompt 3 made it look for the forward
model (χ²_ν 45 against 40). The candidates are the atlas cell (a different
cell, 0.928 correlation), the LSF shape (single Gaussian) and
under-deconvolution. One chunk cannot separate them; prompts 5 and 6 can.

*The iodine wavelength solution disagrees with PypeIt's by half a km/s.*
This is the first direct measurement of the zero point phase 2 blamed, and
it is the right size. But an absolute offset of the Fischer atlas would look
identical in one epoch. It becomes a finding only if it varies from epoch to
epoch or order to order, which is prompt 5's comparison and then prompt 7's.

*Two more `pyodine` behaviours that fail silently:*
- **The Chauvenet re-fit never runs.** Its condition is inverted, and its
  criterion runs on unweighted residuals, so it would flag masked pixels.
- **`red_chi_sq` holds √χ²_ν**, not χ²_ν.

Both belong in the prompt-10 upstream report.

*The raw directory is an input, not storage.* `reduce_run.frame_roles` globs
the night's raw directory, so anything downloaded into it becomes part of the
next reduction of that night. Calibration-only frames for other purposes need
their own directory: `1998aug26_lsf/`, like phase 1's `inspect_only/`.

**Files added / modified.**
- Added:
  - `first_hires_exoplanet/utilities_hires/conf.py`, `pyodine_parameters.py`,
    `timeseries_parameters.py` and `logging.json`;
  - `first_hires_exoplanet/fit_one_chunk.py`;
  - `first_hires_exoplanet/tests/test_adapter_wavelengths.py`;
  - `first_hires_exoplanet/data/one_chunk_fit.csv` and `one_chunk_mc.csv`;
  - `docs/figs/fig_p4_one_chunk.png`.
- Modified:
  - `first_hires_exoplanet/utilities_hires/load_pyodine.py`
    (`fill_wavelengths`, `fill_masked_flux`, `IodineTemplate` index, and
    docstrings);
  - `first_hires_exoplanet/utilities_hires/__init__.py`;
  - `first_hires_exoplanet/deconvolve_template.py` (imports
    `fill_wavelengths`);
  - `first_hires_exoplanet/koa_download.py` (two selections).
- Outside the repository: `../first-hires-exoplanet-data/raw/1998aug26_lsf/`
  (12 frames, 57 MB).
- Nothing in `vendor/pyodine/` was changed.

**Whether any git command changed state.** No. Only read-only `git status`,
`git log` and `git check-ignore`. The prompt-3 files you had staged are as you
left them.
