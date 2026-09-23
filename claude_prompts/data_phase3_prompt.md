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

## Logs
