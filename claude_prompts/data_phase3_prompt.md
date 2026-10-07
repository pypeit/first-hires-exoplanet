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

12. Generate a technical slide deck in PowerPoint describing what you have accomplished in phases 1-3.  Place it in `docs/slides/technical_summary.pptx`.  I will then upload to a GoogleDrive.  Include figures (generate those from Python scripts, either existing or new ones).
    Use Opus 5. Log your work.

13. This is very nice.  Can you modify the text to have a font size no smaller than 20pt?  
    Use Opus 5. Log your work.

14. Can you add a slide to the deck that shows the results from the original 1998 paper?  Include the figure from the original paper that shows the velocity curve.  Use Opus 5. Log your work.

15.  Have you drafted the email for Ryan Cooke?  If so, please send it to me for review.  If not, please draft it now.  Use Opus 5. Log your work.

16. Ok, do your best job at making 2 slides to 
describe all that we have done for the general public.
Empahize: (1) they could to most of this themselves, 
eg. show the ~10 key prompts that I had to make;
(2) emphasize that I could transfer my knowledge to a high school teacher who could lead a series of lessons for their students to do this and understand it.
Use Opus 5. Log your work.

17. I have pulled this branch into the main branch.  Please generate a QR code that points at the public facing report: `https://github.com/pypeit/first-hires-exoplanet/blob/main/docs/public_HD187123b.md` that I will put on a slide.
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

### After prompt 5: decisions before prompt 6

The evidence for each is in `## Report`, "Prompt 5".

1. **The misfit follows the iodine, not the LSF. Should prompt 6 still be
   about the LSF alone?** χ²_ν tracks iodine contrast at partial ρ ≈ +0.7,
   against +0.1 for stellar contrast. The Fischer atlas is a different cell.
   The B-star frames downloaded for prompt 6 can measure both things at
   once: the stellar LSF, and the transmission of the *HIRES* cell at HIRES
   resolution (the B-star spectrum is almost pure iodine). The options:
   - (a) prompt 6 as written: compare LSF models on this epoch, atlas
     unchanged;
   - (b) prompt 6 first fits the B stars with the atlas free to rescale per
     chunk (a per-chunk `iod_depth`, i.e. α), and then compares LSF models;
   - (c) a new step before prompt 6: build an empirical HIRES-cell
     correction from the B stars.

   My recommendation: (b). It stays inside `pyodine`'s own O-star machinery
   and answers whether a per-chunk α removes the correlation, before choosing
   an LSF on top of a known-wrong atlas.

   *Added after prompt 6, which ran as (a) because this was unanswered.* The
   LSF turned out to be a ≤5% lever. The misfit's correlation with iodine
   contrast is the same under every LSF model (partial ρ +0.63 to +0.73).
   That strengthens (b), or (c), as the next step before prompt 7.

   Answer: Ok, please do (b)

2. **Should the index trap be fixed in the fork?** The driver fits
   observation chunk *i* against template chunk *i*, and a subset of orders
   breaks that silently. 11 of the 31 cell-in epochs have at least one
   rejected iodine order. The options:
   - (a) fit every epoch over all 14 orders, with rejected orders at zero
     weight, and drop them afterwards. There is no fork change, but the
     zero-weight chunks enter run 1's medians;
   - (b) fix the fork: store each observation chunk's template index and fit
     against that. It gets its own test, and it goes upstream in the prompt-10
     report.

   My recommendation: (b), before prompt 7.

   Answer: (b)

3. **Should the B stars be reduced as part of prompt 6?** The 12 frames are
   on disk in `raw/1998aug26_lsf/`, unreduced. They need a PypeIt run that
   borrows 1998-08-26's calibrations. My recommendation: yes, as the first
   step of prompt 6.

   Answer: use your recommendation

### After the prompt-5 decisions were carried out: before prompt 7

The evidence is in `## Report`, "Prompt 6, continued".

1. **The atlas's line pattern, not its depth, limits the fit. Attempt an
   empirical HIRES-cell spectrum before prompt 7?**

   Fixing each chunk's α to its B-star value changed nothing significant,
   and the pure-iodine B-star fits misfit where the iodine is strongest,
   just as HD 187123 does. The one remaining fix is option (c) from the
   round before: build the HIRES cell's own transmission from the B stars
   and use it in place of the Fischer atlas. The B stars give that at the
   stellar LSF, so it would need deconvolving, or an atlas-plus-correction
   form. A real FTS scan of the Keck cell would be better still, but none
   has been found.

   The options:
   - (a) prompt 7 now, with the current model. The epoch velocity already
     agrees with the prediction to ~1–2σ at ±12 m/s. Prompt 7's
     epoch-to-epoch scatter measures what the misfit actually costs. Revisit
     (c) if it exceeds the target.
   - (b) (c) first: an empirical cell correction from the 10 B-star frames,
     tested on the prompt-5 epoch, then prompt 7.

   My recommendation: (a). The question that decides whether (c) is worth
   the effort is "what does the misfit cost epoch to epoch", and prompt 7
   answers it. The misfit is common to every epoch's chunks, so it may cancel
   largely in relative velocities.

   Answer: (a)

### After prompt 9: before prompts 10–11

The evidence is in `## Report`, "Prompt 9".

1. **Which velocities are phase 3's result?** Prompt 7's combination uses
   all 31 epochs and every usable order. It gives 25.8 m/s per epoch and
   K = 72.4 ± 7.1 m/s (10σ). Prompt 9 found two weak groups, and neither
   criterion uses the Keplerian:
   - the red orders 58–61, where the iodine is weak and whole orders shift
     together;
   - six epochs with ≤ 8 usable orders.

   Leaving both out gives 11.5 m/s and K = 71.8 ± 3.5 (20σ) on 25 epochs.
   But the decision to try it came after the residuals were seen. The
   options:
   - (a) keep prompt 7 as the headline, and report the prompt-9 variants as
     "what a better combination would reach";
   - (b) make it a method change, stated before any Keplerian is looked at:
     weight each *order* by the scatter of its own time series (the
     chunk-level weighting prompt 7 already uses, one level up), recompute
     the errors, and re-run prompt 8. Dropping orders by hand would then be
     unnecessary;
   - (c) adopt "without 58–61 and the weak epochs" as the result.

   My recommendation: (a) for prompt 11's headline, because it is the result
   made without hindsight. (b) could be a documented follow-up; it is
   principled, but anything done now has been seen. Not (c).

>A. (a)

### After prompt 10: sending the upstream reports

The drafts are `docs/upstream/pypeit_report_draft.md` and
`docs/upstream/pyodine_report_draft.md`. Nothing has been sent.

1. **Are the drafts right to send?** Edit them in place or say what to
   change. The PypeIt draft carries a placeholder, `<repository link>`, which
   needs this repository's public URL.
>A. I will send them
2. **How should the PypeIt items go?** The options:
   - (a) email Ryan as drafted, and you send it;
   - (b) GitHub issues on `pypeit/PypeIt`, one per item or group, plus a PR
     from `orig-hires-fixes` for items 1–5;
   - (c) both: the email as the overview, then the issues and the PR.

   My recommendation: (c). Items 7, 8, 10 and 13–15 affect other users and
   want tracking; the email gives Ryan the context. I will not open issues
   or a PR myself unless you ask, because that is outward-facing.
>A. (c)
3. **Should the `pyodine` report go at all, and how?** By email to the
   authors, or as GitHub issues on `pepeheeren/pyodine` (a one-commit
   repository with no issue traffic, so email is more likely to be read). My
   recommendation: email first, offering PRs.
>A. As Issue(s) on GitHub

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

### Prompt 5: one epoch, end to end — right velocity, wrong iodine

**In one paragraph.** `pyodine`'s own driver fits all 700 chunks of an epoch
in under a minute. The epoch's median velocity agrees with the prediction to
12–18 m/s, with a statistical error of about 13 m/s. But each chunk scatters
6–7 times more than its photon noise, the model misfits at χ²_ν ≈ 23, and the
misfit follows the iodine, not the star. The fitted wavelength solution sits
about 1 km/s from PypeIt's, with structure. The two prompt-3 templates are
nearly indistinguishable at this level.

`first_hires_exoplanet/fit_one_epoch.py` produces everything below. It runs
`pyodine_model_observations.model_single_observation` unchanged through
`utilities_hires`: run 0 is a single Gaussian and run 1 is Lick's
multi-Gaussian. It does so for both templates, then analyses the saved
results.

**The epoch, and a trap that decided it.** The epoch is `HI.19980825.19425`:
cell in, 430 s, S/N 167. It is the best cell-in frame whose 14 iodine orders
all pass the quality filter, and that condition is necessary:
- The driver fits observation chunk *i* against template chunk *i*.
- But `auto_wave_comoving` builds observation chunks only for the orders it
  is passed.
- So any subset of orders silently pairs every chunk after the gap with a
  template chunk from a different order. The spline extrapolates, and there
  is no error.
- Prompt 4's epoch (19326) has orders 58–59 rejected, so it cannot be run
  correctly as the code stands.

For prompt 7 this means one of two things: epochs must be fitted over all
template orders with rejected orders at zero weight, or the driver must be
fixed.

The epoch's raw OBJECT reads `Gl 83.1`, a stale keyword. The pointing
(19:46:58 +34:25:10) and TARGNAME are HD 187123's.

**Getting it to run: one more adapter repair.** The first full run fitted all
700 chunks and then died writing the results: `results_io` reads
`obs.star.proper_motion`, which the adapter's `Star` lacked. The driver
reported this with `logging.info`, so **its error log was empty**. The adapter
now carries `proper_motion`, with a test that fits a real chunk and
round-trips it through `save_results`/`load_results`. This driver checks the
info log and the existence of every result file, never the error log.

**Runtime.**

| template | total | run 0 | run 1 |
|---|---|---|---|
| os 10 | 58.8 s | 16.6 s | 39.1 s |
| os 1 | 47.6 s | 13.7 s | 30.7 s |

That is serial, on one core, and includes the driver's plots. All 31 epochs
would be about half an hour, before any parallelism.

**How many chunks fail** (700 in each case):

| | os 10, run 0 | os 10, run 1 | os 1, run 0 | os 1, run 1 |
|---|---|---|---|---|
| lmfit failed | 0 | 0 | 0 | 0 |
| no uncertainties | 59 | 60 | 59 | 60 |
| a parameter at a bound | 2 | 0 | 0 | 0 |
| misfit (χ² > 3× order median) | 64 | 66 | 59 | 71 |
| velocity outlier (> 5 robust σ) | 35 | 24 | 32 | 23 |
| **good** | **552** | **562** | **560** | **556** |

Nothing fails loudly: every failure is a fit that returned. The
no-uncertainty chunks are spread over every order, are not tied to masked
pixels, and fit worse (median χ²_ν 34 against 24). The Chauvenet step
re-fitted nothing, as prompt 4 predicted.

**The per-chunk velocity scatter.** Good chunks, m/s:

| | os 10, run 0 | os 10, run 1 | os 1, run 0 | os 1, run 1 |
|---|---|---|---|---|
| median velocity | −684 | −677 | −694 | −671 |
| expected (barycentric + planet) | −689 | −689 | −689 | −689 |
| robust σ per chunk | 340 | 301 | 339 | 287 |
| standard deviation, good chunks | 483 | 434 | 494 | 441 |
| standard deviation, all chunks | 2245 | 808 | 2029 | 1012 |
| median lmfit error (misfit-scaled) | 179 | 216 | 174 | 224 |
| median photon-noise error | 40 | 48 | 39 | 48 |
| epoch error, robust σ / √n | 14.5 | 12.7 | 14.3 | 12.1 |
| order-to-order scatter of medians | 154 | 106 | 91 | 64 |
| median within-order robust σ | 292 | 240 | 364 | 273 |
| median χ²_ν, measured noise | 25.7 | 23.6 | 25.9 | 22.7 |

Five things follow:
- **The epoch is where it should be.** The epoch velocity agrees with the
  prediction to 5–18 m/s, and its statistical error is ~13 m/s. The median is
  robust to the selection: −677 good, −679 all finite, −682 non-outlier (os
  10, run 1). *If* the chunk errors are independent, one epoch already beats
  the 50 m/s target. Prompt 7 is where that "if" gets tested.
- **The per-chunk scatter is 6–7× the photon noise, and 1.4× even the
  misfit-inflated lmfit error.** The error bars understate the scatter even
  after absorbing the misfit.
- **The standard deviation is the wrong statistic here too.** It is 808–2245
  m/s over all chunks, set by a few dozen outliers.
- **Run 1 (multi-Gaussian) improves on run 0** in every row: −12% robust
  scatter, −30% order-to-order.
- **Red orders scatter most.** By order (os 10, run 1), the robust σ is
  410–600 m/s in orders 58–61, where the photon error is also 90–140 m/s, and
  130–250 m/s in orders 63–71, where the photon error is 23–40 m/s.

**What the misfit follows: the iodine, not the star.** For each chunk, the
contrast of the α-scaled atlas and of the template (RMS about the mean) was
set against χ²_ν. Spearman correlations over ~560 good chunks:

| | iodine | stellar | partial, iodine | partial, stellar |
|---|---|---|---|---|
| os 10, run 0 | +0.76 | +0.39 | **+0.72** | +0.17 |
| os 10, run 1 | +0.73 | +0.34 | **+0.69** | +0.11 |
| os 1, run 0 | +0.77 | +0.34 | **+0.73** | +0.11 |
| os 1, run 1 | +0.71 | +0.33 | **+0.67** | +0.11 |

χ²_ν by order peaks at ~45 in orders 65–68 and falls to ~10 at both ends,
tracking the iodine band's strength.
- An LSF error would misfit stellar and iodine structure alike.
- A template error would follow the stellar contrast.
- What the data show is misfit tracking the iodine with the star held fixed.

That is the signature of **the iodine model**: the Fischer atlas is a
different cell (phase 2: 0.928 line correlation), with one Beer–Lambert α.
It is the most important thing prompt 5 found for prompt 6. **Choosing an
LSF model cannot fix an atlas that is the wrong cell.**

**The LSF as a function of position.**

- **Run 0 (single Gaussian), os 10:**
  - median FWHM **2.00 px**; per-order, per-quarter medians range 1.76–2.60
    px;
  - the red orders (58–60) are broadest, at 2.2–2.6 px;
  - in the blue orders the width rises along the order, from ~1.8 px at pixel
    250 to ~2.05 px at pixel 1300 (fig. panel d).
- **Run 1 (multi-Gaussian):**
  - FWHM 1.99 px;
  - half-maximum asymmetry −0.01 to −0.02, essentially symmetric;
  - but a **centroid offset of +0.03 px, about 60 m/s**.

  Lick's multi-Gaussian is deliberately not re-centred. So that offset is a
  velocity zero point the LSF absorbs, and it is shared, degenerately, with
  the velocity and the wavelength intercept.
- **The same night's ThAr is 2.51 px** (1754 lines), against 2.18 px on
  08-26. So star/ThAr is 0.78–0.80, and the arc width moved 15% between
  nights while the star's stayed at ~2.0 px (prompt 4 found 2.07 on 08-26).
  **The ThAr is not a usable proxy for the stellar LSF in either direction.**
  That is one more reason prompt 6 needs the B stars.

**The wavelength solution against PypeIt's.** For each good chunk, the
fitted line (`wave_intercept`, `wave_slope`) was compared with PypeIt's
`OPT_WAVE` over the same pixels:

| | os 10, run 1 | os 1, run 1 |
|---|---|---|
| median offset over orders | **−954 m/s** | **−941 m/s** |
| order-to-order scatter | 151 m/s | 129 m/s |
| within-order robust σ | 178 m/s | 191 m/s |
| trend along an order | +15 m/s per 1000 px | +9 m/s per 1000 px |
| dispersion ratio | 0.99986 | 0.99968 |

Run 0 gives −1020 m/s and a +82–96 m/s/kpx trend, because its dispersion is
smoothed over the order before the fit.

Across the band the offset runs from about −0.7 km/s at 5000 Å to −1.0 km/s,
with sharp excursions of 1–2 km/s at the red ends of the reddest orders (fig.
panel e). That is the shape a per-order polynomial takes where it is poorly
constrained.

The offset is not constant between nights: prompt 4 measured −0.54 to −0.57
km/s on 08-26 for one chunk. A fixed atlas offset cannot do that; PypeIt's
per-epoch ThAr zero point can. This is the phase-2 diagnosis appearing
directly, but it rests on two epochs and one chunk of one of them, and
prompts 7 and 9 will test it properly.

**The continuum against PypeIt's.** The model normalises the iodine and the
template by their chunk means, so the raw `cont_intercept` is not a
continuum. Once de-normalised (`cont × ⟨T⟩^(−d) × (d_t(1/⟨S⟩−1)+1)`), it is
the flux the model puts at a line-free pixel:

| | os 10, run 1 | os 1, run 1 |
|---|---|---|
| chunk-to-chunk scatter vs OPT_FLAT, de-normalised | 2.40% | 2.40% |
| the same, raw `cont_intercept` | 5.68% | 6.42% |
| de-normalised / the adapter's 95th-percentile envelope | 1.099 ± 0.037 | 1.101 ± 0.036 |

The ratio to PypeIt's blaze has the **same U shape in every order**: +8% at
both ends of an order, −5% in the middle (panel f). That is a smooth
difference between the star's and the flat lamp's illumination of the blaze.
The flats are 1998-08-12 donors through a different decker, and scattered
light is a plausible cause. It does not matter to a chunk-by-chunk continuum,
but it would to anyone using `OPT_FLAT` as the blaze. The adapter's envelope
sits 10% below the continuum the model finds, as a running percentile does in
a line-dense band.

**The two templates.** At the epoch level they are close to
indistinguishable. Oversampling 1 is marginally better in run 1: robust σ
287 against 301, order-to-order 64 against 106 m/s, χ²_ν 22.7 against 23.6.
Its median velocity is 18 m/s further from the prediction (os 10 is 12 m/s
off). Prompt 3's re-convolution failure of the os-10 template is real, but
it is not what limits the fit now. The iodine model is.

**Products.**
- **Code:** `first_hires_exoplanet/fit_one_epoch.py`.
- **Tables:** `first_hires_exoplanet/data/one_epoch_chunks.csv` (every chunk,
  both templates and runs, 1.2 MB), `one_epoch_summary.csv`,
  `one_epoch_lsf.csv`.
- **Figure:** `docs/figs/fig_p5_one_epoch.png`.
- **Test:** `first_hires_exoplanet/tests/test_adapter_results_io.py`.
- **`pyodine`'s own results and plots:**
  `../first-hires-exoplanet-data/pyodine_runs/HI.19980825.19425_os{10,1}/`.

### Prompt 6: the LSF is settled, and it is not the lever

**Scope.** The three Q&A questions after prompt 5 were unanswered when this
prompt ran. It therefore ran **as written**, which is option (a) of question
1: LSF models compared on the prompt-5 epoch, atlas unchanged. The B stars
were not reduced (question 3) and the fork was not changed (question 2).
All three remain yours.

**The answer.** The prompt expected the LSF to be "the largest single lever
on the final precision". **It is not.**
- Of the five models `pyodine` offers, the four that work give the same
  per-chunk velocity scatter to within ±5%, and no paired difference is
  significant.
- The fifth (Hermite) is worse and biased.
- Holding a smoothed LSF fixed (run 2) never helps.

The chosen model is the **super-Gaussian**, with the **oversampling-1
template**, and it is now the `utilities_hires` default. The lever that
remains is the iodine model: prompt 5's diagnosis is unchanged by any LSF.

**What was compared.** `first_hires_exoplanet/compare_lsf_models.py` ran
prompt 5's driver on `HI.19980825.19425` ten times: five run-1 models × two
templates. Each run had the common run 0 (single Gaussian) first and a run 2
after (run 1's LSF smoothed ±160 px and ±3 orders, Lick's radii, then held
fixed as `FixedLSF`).

| name | model | LSF parameters | configuration |
|---|---|---|---|
| single | `SingleGaussian` | 1 | FWHM 0.5–4 px |
| super | `SuperGaussian` | 4 | σ 0.2–3 px, exponent 1–4, satellites 0–1 |
| multi_lick | `MultiGaussian_Lick` | 10 | Lick's layout, Lick's bounds (prompt 5's run 1) |
| multi_song | `MultiGaussian` | 10 | SONG's layout, Lick's bounds |
| hermite | `HermiteGaussian` | 7 | SONG's enabled weights 3–8, each −0.5 to 0.5 |

Two configurations depart from upstream on purpose:
- **SONG's Hermite bounds** pin every weight within ±2×10⁻¹³ of ~10⁻⁹,
  which makes the model a single Gaussian.
- **Lick's generic bounds** would allow the super-Gaussian a negative
  exponent.

The smoothing's order separation is HIRES's: 71 physical pixels (36 binned),
measured from the spec1d traces, where Lick uses 15.

**A trap the first attempt fell into.** The first pass gave the
super-Gaussian a per-chunk σ of exactly 0 and every Hermite chunk no
uncertainties. Neither fit had moved from its start.
- Lick's run-1 scheme starts from `fit_lsfs`, which fits the new model to
  run 0's Gaussian. That sends every parameter a Gaussian does not need
  (satellites, Hermite weights) to about 10⁻¹².
- `lmfit`'s `leastsq` takes finite-difference steps *relative to the value*,
  which is about 10⁻²⁰ there. The derivatives are zero, and MINPACK returns
  the start after one Jacobian: 16 evaluations for 11 parameters, every
  value at `init`.
- It is probably why SONG's own Hermite configuration pins the weights.

Fixed by starting each shape parameter at least 10⁻³ from zero and 2%
inside its bounds (`pyodine_parameters.start_inside`, tested). It belongs in
the prompt-10 report. Everything below is from the fixed run.

**The result.** Good chunks, run 1. Robust σ and order-to-order scatter are
in m/s; the epoch median is against −689 m/s expected (error ≈ ±12–17):

| model | os | good | robust σ | order-to-order | χ²_ν (measured) | epoch median | no-unc | run-1 time |
|---|---|---|---|---|---|---|---|---|
| single | 10 | 565 | 319 | 80 | 24.2 | −684 | 59 | 12 s |
| super | 10 | 554 | **295** | 109 | **23.3** | −681 | 59 | 21 s |
| multi_lick | 10 | 562 | 301 | 106 | 23.6 | −677 | 60 | 37 s |
| multi_song | 10 | 544 | 304 | 72 | 28.5 | −669 | 58 | 31 s |
| hermite | 10 | 544 | 397 | 153 | 32.7 | **−551** | 60 | 39 s |
| single | 1 | 569 | 289 | 52 | 25.0 | −669 | 59 | 11 s |
| super | 1 | 564 | 288 | **50** | **22.6** | −665 | 59 | 18 s |
| multi_lick | 1 | 556 | 287 | 64 | 22.7 | −671 | 60 | 29 s |
| multi_song | 1 | 541 | **278** | 77 | 27.3 | −658 | 56 | 28 s |
| hermite | 1 | 552 | 375 | 162 | 34.9 | **−551** | 60 | 36 s |

**Is any difference real?** Paired bootstrap over the chunks both models
call good, against `multi_lick` run 1 on the same template. The change in
robust σ, with its 68% interval:

| model | os 10 | os 1 |
|---|---|---|
| single | +13 [−11, +30] | +3 [−16, +30] |
| super | −3 [−25, +15] | −0 [−17, +21] |
| multi_song | +4 [−12, +28] | −11 [−26, +13] |
| **hermite** | **+92 [+70, +120]** | **+108 [+91, +140]** |

The four workable models are statistically indistinguishable, and Hermite is
clearly worse. Switching among the workable four moves individual chunk
velocities by 150–180 m/s (robust) but the epoch by only 5–21 m/s. Hermite
moves the epoch by +92 to +117 m/s: it finds a spurious −0.08 asymmetry
where every other model finds the profile symmetric (≤0.02).

**Run 2, the smoothed LSF held fixed, does not help.**
- Robust σ changes by −13 to +63 m/s against run 1, never significantly
  better.
- It loses 30–40 good chunks for the multi-Gaussians and raises χ²_ν by 1–4.
- It moves the epoch by up to 3.3σ (multi_lick os 1: −644) and 5σ
  (multi_song os 1: −604).

The per-chunk LSF is doing work, and smoothing it away undoes that work. The
most likely work is absorbing some of the iodine misfit chunk by chunk, since
the misfit's correlation with iodine contrast is unchanged by every model
(partial ρ +0.63 to +0.73) while χ² rises when the LSF is smoothed.

**What the profile is.** Every workable model converges on the same shape:
- FWHM **2.00–2.05 px**;
- symmetric, with half-maximum asymmetry ≤ 0.02 in magnitude;
- centred, except Lick's non-re-centred multi-Gaussian, which carries a
  +0.03–0.04 px centroid (~60–75 m/s) that is degenerate with the velocity
  and wavelength zero point.

The per-order, per-position variation of prompt 5 (1.6–2.6 px, broadest in
the red orders) is common to all models. It is the profile, not the model.

**The choice, and why.**
1. **The super-Gaussian (run 1, no run 2).** On velocity scatter it ties for
   best with both templates; that is the deciding criterion, and the tie is
   the point. The tie-breaks all favour it:
   - the lowest χ²_ν with both templates;
   - the lowest order-to-order scatter with the oversampling-1 template;
   - centred, unlike Lick's multi-Gaussian;
   - 4 parameters where the multi-Gaussians have 10, and no parameter at a
     bound;
   - run 1 in about 60% of the multi-Gaussians' time.

   The single Gaussian is the statistically equivalent fallback. BIC
   prefers `multi_lick`, but at χ²_ν ≈ 23 BIC is scoring systematic misfit
   and rewards flexibility without improving a velocity, so it was not used
   to decide.
2. **The oversampling-1 template**, which settles the Q&A decision after
   prompt 3 to let the fit decide. Paired, same model, oversampling 1
   against 10:

   | model | change in robust σ |
   |---|---|
   | single | −28 [−40, −10] |
   | super | −11 [−29, +2] |
   | multi_lick | −23 [−43, −4] |
   | multi_song | −28 [−53, −15] |

   It is also consistent with prompt 5. The two templates' epochs differ by
   +6 to +8 m/s (a zero point, irrelevant to relative velocities).

`utilities_hires/pyodine_parameters.py` now sets run 1 to the
super-Gaussian, with physical bounds and `start_inside`, in both
`Parameters` and `Template_Parameters`, and `osample_temp = 1`. Four new
tests pin this.

**The levers, ranked** (per-chunk robust σ on this epoch):

| lever | effect |
|---|---|
| per-chunk scatter against photon noise | **×6–7** (≈290 against ≈45 m/s): what is left to win |
| the iodine model (prompt 5) | the dominant term, and not yet varied: partial ρ ≈ +0.7 with iodine contrast, under every LSF |
| a broken configuration (Hermite as tried) | +30% |
| run 2, smoothed and fixed LSF | 0 to +25% |
| run 0 → run 1 (prompt 5) | −12% |
| template oversampling 10 → 1 | −4% to −9% |
| **LSF model, among the four that work** | **≤ ±5%, not significant** |

**Products.**
- **Code:** `first_hires_exoplanet/compare_lsf_models.py`;
  `utilities_hires/pyodine_parameters.py` (settled).
- **Tables:** `first_hires_exoplanet/data/lsf_models_summary.csv` (40 rows),
  `lsf_models_chunks.csv` (7000 rows, 1.9 MB) and
  `lsf_models_template_pairs.csv`.
- **Figure:** `docs/figs/fig_p6_lsf_models.png`.
- **Test:** `first_hires_exoplanet/tests/test_lsf_settled.py`.
- **`pyodine`'s results:** `../first-hires-exoplanet-data/pyodine_runs/lsf_*`
  (10 directories).

### Prompt 6, continued: the Q&A decisions after prompt 5, carried out

**The three decisions.** Taken after prompt 6 had run, and acted on here:

| question | decision |
|---|---|
| 1 | (b): fit the B stars with a per-chunk α, then compare LSF models |
| 2 | (b): fix the index trap in the fork |
| 3 | reduce the B stars first |

**The answer, in brief.**
- **The index trap is fixed in the fork**, with a test.
- **The B stars are reduced.** 10 of 11 frames are usable, after one PypeIt
  parameter and a switch to the boxcar extraction.
- **The iodine depths are not the problem.** The HIRES cell is α ≈ 3.33 with
  a gentle wavelength trend.
- **The iodine line pattern is.** Fixing each chunk's α to its B-star value
  leaves the epoch's scatter, χ² and iodine-correlated misfit unchanged, and
  the pure-iodine B-star fits show the same correlation.
- **On pure iodine the LSF models are again near-equivalent.** The
  super-Gaussian choice stands.

**1. The index trap, fixed in the fork (change 5).**
- `chunks.auto_wave_comoving` now records each chunk's template index
  (`chunk.template_index`).
- `SimpleModel.eval`, `SimpleModel.clean_of_I2` and the four
  template-plotting sites in `plot_lib` resolve the template through a new
  `models.spectrum.template_chunk_index(chunk, chunk_ind)`. It falls back to
  the list index for chunks built any other way.
- `tests/test_fork_chunk_index.py` builds one epoch's chunks for all orders
  and for a subset. It failed before the change, when the same chunk gave a
  different model from each list, and passes after.
- `vendor_pyodine --verify` shows exactly the intended diff.
- `vendor/README.md` now lists every change to the fork. It had still been
  claiming "byte-identical to upstream" since prompt 2.
- It should go upstream in the prompt-10 report. Any epoch may now be fitted
  over any subset of orders.

**2. The B stars, reduced.**
`first_hires_exoplanet/reduce_bstars.py` stages the night's own ThAr arc and
the 1998-08-12 donor flats with the 11 B-star frames. It reuses 1998-08-26's
`reduce_run.PARAM_BLOCK` and copies the night's calibrations rather than
rebuilding them. Output goes to `redux/reduce_19980826_lsf/`, and the HD
187123 reduction is untouched. Two things were needed:

- **`force_center_obj = True`**, the remedy built for Shane/Hamspec's
  slit-filling stars. PypeIt at `017bece06` has it.
  - Without it, a V ≈ 3–4 star in a 3.5″ slit defeated peak-finding: HR 8634
    (37323) was found in 1 of 37 orders, and the echelle extraction stopped
    the run.
  - With it, one object is placed at each order's centre and extracted with
    a boxcar across the order. Sky subtraction was already off.
  - It is applied to the B stars only.
- **The boxcar extraction, not the optimal one.** On these slit-filling stars
  the optimal extraction's mask collapses: 8–38% of iodine-order pixels kept
  in most frames. The boxcar's keeps 99%, at 90–95% of the S/N.
  - The adapter gained `extraction='OPT'|'BOX'` (default `OPT`, unchanged),
    with a test.
  - The boxcar does not mask the detector's bad columns (runs of 3–4 pixels
    at ~15% of the continuum), so `pyodine`'s own `BadPixelMask` is applied.
    `Template_Parameters` asks for it, and it took the first B-star fit from
    χ²_ν ≈ 1800 to ≈ 44.
  - `create_template` computes that mask on the *template* observation and
    applies it to the hot star's weights. Here it is computed on each B-star
    frame.
- **HR 838 at 20 s (55734) is unusable.** The raw frame peaks at 773 ADU:
  the star missed the slit.

**3. The iodine cell through the B stars** (`fit_bstars.py`: `pyodine`'s
hot-star path, the template's chunk grid, velocity and template depth fixed,
`iod_depth` free per chunk, run 1 the super-Gaussian).

| | |
|---|---|
| noise, HR 8634 pair differenced pixel by pixel | 1.61 × propagated (boxcar); HD 187123's optimal extraction: 0.43× |
| χ²_ν, measured noise, per frame | HR 838 1.0–1.2; HR 8634 2.0–4.0; HR 7236 3.9–5.2 |
| α, map over 10 frames | median **3.33**; between chunks 0.34; frame to frame per chunk 0.31; one chunk's error 0.44 |
| α per frame | 3.29–3.45 |
| α against wavelength / iodine contrast | Spearman +0.43 / −0.22 |
| χ² against iodine contrast | Spearman **+0.40** |
| LSF FWHM | 2.10 px (per frame 2.01–2.17), symmetric; rising along every order from ~2.0 px to ~2.2 px |
| fitted wavelengths against PypeIt | −143 m/s median; **−392 to +95 m/s from frame to frame** |

Five things follow:
- **The cell is thicker than the adapter assumes.** α ≈ 3.33, not 2.59.
  Phase 2's 2.28–2.86 came from two 60 s frames. Between chunks α varies
  less than one chunk's error: a band-wide α with a gentle red-ward rise
  (panel a) is all the data support.
- **The misfit that follows the iodine is there in pure iodine too**
  (panel b). χ² climbs steeply above an iodine contrast of ~0.3, as it does
  for HD 187123. It is not the star, the template or anything stellar.
- **The noise scale is extraction-specific.** The boxcar's inverse variance
  *understates* the noise (1.61) where the optimal one overstates it (0.43).
  A χ² is comparable only within one extraction. In PypeIt's own units the
  B-star χ²_ν (~3–13) is not clearly better than HD 187123's (~4), so "the B
  stars fit to the noise" would overstate it.
- **The iodine wavelength scale moves against PypeIt's through the night**,
  by half a km/s between B-star frames an hour apart. PypeIt's single
  end-of-night arc cannot follow that, and slit-filling B stars may also
  shift with their position in the slit. This is prompt 9's intra-night
  question, seen from the calibration side.
- **The B-star LSF (2.1 px) matches the star's (~2.0 px)** and this night's
  ThAr (2.18). Its trend along the order is the one prompt 5 saw in the star.

**The LSF models on pure iodine** (HR 8634 pair, run 1; median χ²_ν,
measured noise):

| model | 37213 | 37323 | α | FWHM |
|---|---|---|---|---|
| single | 3.14 | 3.00 | 3.19–3.21 | 2.09–2.11 |
| super | 3.12 | 3.01 | 3.30–3.36 | 2.11–2.14 |
| multi_lick | 3.41 | 3.19 | 3.20–3.24 | 2.05–2.07 |
| multi_song | 4.04 | 3.63 | 3.22–3.29 | 2.11–2.13 |
| hermite | 4.62 | **29.6** | 1.84–3.05 | 1.43–2.09 |

The single and super-Gaussians fit pure iodine best by median χ² (Lick's
multi-Gaussian has the lowest total, fewer bad tails), and α and the FWHM
barely depend on the model. Hermite is unstable. The prompt-6 choice stands.

**4. The payoff: the prompt-5 epoch with each chunk's α fixed to its B-star
value** (super-Gaussian, oversampling-1 template):

| | good | robust σ | order-to-order | χ²_ν (measured) | epoch − expected | partial ρ, iodine |
|---|---|---|---|---|---|---|
| `iod_depth` free | 564 | 288 m/s | 50 m/s | 22.6 | +24 ± 12 m/s | +0.69 |
| `iod_depth` fixed, B stars | 560 | 276 m/s | 77 m/s | 23.4 | +28 ± 12 m/s | +0.67 |

Paired, fixed against free: −7 [−26, +7] m/s, not significant. **A better
depth model buys nothing.** What limits the fit is the line pattern of an
atlas scanned from a different cell. Neither an α nor an LSF can supply it;
only the HIRES cell's own spectrum can. That is option (c) of question 1.
The new Q&A below asks whether to attempt it before prompt 7.

**Products.**
- **Fork:** `pyodine/chunks.py`, `pyodine/models/spectrum.py`,
  `pyodine/plot_lib.py`; `vendor/README.md`.
- **Code:** `first_hires_exoplanet/reduce_bstars.py`,
  `first_hires_exoplanet/fit_bstars.py`; `utilities_hires/load_pyodine.py`
  (`extraction`).
- **Tests:** `tests/test_fork_chunk_index.py`, `tests/test_adapter_boxcar.py`.
- **Tables:** `first_hires_exoplanet/data/bstar_chunks.csv` (2.0 MB),
  `bstar_alpha_map.csv`, `bstar_lsf_models.csv`, `bstar_alpha_payoff.csv`.
- **Figure:** `docs/figs/fig_p6b_bstars.png`.
- **Outside the repository:** `redux/reduce_19980826_lsf/` (11 spec1d),
  `stage_`/`setup_19980826_lsf/`, `pyodine_runs/bstar_alpha_*`.

### Prompt 7: 31 epochs, and the misfit cancels between them

**In one paragraph.** All 31 cell-in epochs fitted, with nothing changed
since prompt 6 (Q&A answer (a)). The epoch velocities scatter by **55 m/s
rms**, against 733 m/s for phase 2's cross-correlation on the same frames,
and most of that 55 m/s is the planet. The adopted per-epoch error is **19
m/s** (median), and two checks that share no machinery with it say it is
honest.
- The per-chunk scatter within an epoch is still 6–7× the photon noise, and
  still misfit-dominated.
- But the misfit is the same in every epoch. `pyodine`'s own time-series
  combination measures each chunk's offset and removes it, and that is what
  turns a 290 m/s chunk scatter into a 19 m/s epoch.

**What ran.** `first_hires_exoplanet/fit_all_epochs.py` has two modes:
- `--run`: `pyodine`'s driver, one epoch at a time, in three shards of 4
  cores each. That is 35–46 s per full-format epoch, about 20 minutes of fitting
  in all, and 8 minutes of wall time.
- the analysis: per-chunk tables, the combination, the tables and the
  figure.

The configuration is prompt 6's: run 0 a single Gaussian, run 1 the
super-Gaussian, the oversampling-1 template, PypeIt inverse-variance
weights, α = 2.59 with a free per-chunk `iod_depth`. Run 1's velocities are
used throughout.

**Which orders.** Each epoch is fitted over the iodine orders phase 2's
quality filter passes, handed to the driver as its `orders` argument. That is
the first production use of fork change 5.

| orders fitted | epochs |
|---|---|
| 14 | 20 |
| 11–13 | 5 |
| 6–8 | 3 |
| 2 | 3 (07-17, 07-18, 07-19, the last frame of each night) |

18,350 chunks were fitted and 13,950 are good in prompt 5's sense. The
median χ²_ν, against the measured noise, is 21.4, as in prompt 5.

**Two traps found on the way.** Both were silent.
1. **An order missing from the extraction is bridged.** `pyodine` maps
   template order *o* to observation position *o* + `order_correction`, a
   single shift. 1998-09-13 was extracted with 34 orders and has no order 71
   (the gap phase 2 left undiagnosed). So the shift sent template order 71
   onto echelle order 72. Fifty chunks were fitted against the wrong
   wavelengths and returned velocities of up to 7.7 × 10⁶ m/s. There was no
   error, and the epoch's `pyodine` error came out at 26 km/s.
   - `fit_all_epochs.usable_orders` now keeps a template order only if its
     position holds the same echelle order. `tests/test_epoch_orders.py`
     tests it.
   - It is fixed in the caller, not the fork. It belongs in the prompt-10
     report beside the index trap, since it is the same assumption one level
     up.
   - Only 09-13 is affected, and it now fits 13 orders.
2. **`pyodine`'s example-chunk plots index the chunk list by fixed numbers**
   (`plot_chunks = [150, 250, 400]`). An epoch with two usable orders has 100
   chunks, and the plot step raises `IndexError` after the results are
   written. My driver treats any traceback as failure (prompt 5's rule), so
   this was caught, not missed. The list is now trimmed per epoch.

Two smaller `pyodine` robustness problems in the combination:
- `robust.mean` divides by zero on a column with one finite value, so chunks
  seen in fewer than 3 epochs are blanked (29 of 700).
- The module's default weighting dictionary is mutated by the first call
  that uses it implicitly, so a copy is always passed.

**The combination.** `timeseries.combine_vels.combine_chunk_velocities` is
the iSONG algorithm, run with `pyodine`'s default (SONG) weighting parameters
on an (epoch × template chunk) array. NaN marks any chunk not fitted or with
no velocity. The steps:
1. Each chunk's velocity is referred to its epoch's robust mean.
2. Its **offset** is the robust mean of that over the time series.
3. It is weighted by the robust scatter of its own time series.
4. It is re-weighted per epoch by its deviation.

In the barycentre, `v_bary = v + BVC_epoch − BVC_template`: additive, as in
phase 2 (the multiplicative term is < 2 m/s here).

**The three scatters, separately.**

| | m/s |
|---|---|
| **per chunk, within an epoch**: good chunks, robust, as prompt 5 measured it | **290** (median over epochs) |
| the same, each chunk's offset removed | 230 |
| all finite chunks, offset removed (`pyodine`'s `c2c_scatter`) | 383, range 303–683 |
| chunk offsets: spread between chunks, the part common to every epoch | 201 |
| each chunk's own time-series scatter | 311 (median) |
| **per order, within an epoch**: robust scatter of the order means | **55** (median), range 5–117 |
| **epoch to epoch**: rms over 31 epochs, planet included | **55.3** (robust 66.2) |
| within a night (7 nights with more than one frame) | median rms 27, worst spread 66 |
| between nights (19 nightly means) | 52.2 |
| phase 2's cross-correlation, the same 31 epochs | 733 |

Three things follow:
- **The misfit that has dominated since prompt 4 is mostly a fixed pattern.**
  It is a 201 m/s spread of chunk offsets that is the same in every epoch.
  The iodine pattern of a different cell misfits each chunk the same way
  every time, so relative velocities lose it. That is the case the Q&A
  answer (a) bet on.
- **The per-order scatter is where the systematics show.** Order by order,
  over the time series (order velocity − epoch velocity):
  - orders 64–70: **31–42 m/s** rms;
  - orders 62–63 and 71: 56–65 m/s;
  - orders 58–61: **97–198 m/s**.

  The red end of the iodine band is again the worst, as in prompt 5. The
  combination weights it down by its time-series scatter.
- **The estimator matters.** Against the catalogue orbit (the magnitude check
  below), the residual rms is:

  | estimator | residual rms |
  |---|---|
  | plain median of good chunks (prompt 5's estimator) | 48.5 m/s |
  | median after removing the chunk offsets (`mdvel`) | 38.9 m/s |
  | the weighted combination | **26.2 m/s** |
  | the weighted combination, good chunks only | 26.3 m/s |

  Removing the offsets and weighting by each chunk's time series are both
  worth having. Pre-selecting "good" chunks adds nothing, because the
  re-weighting already does it.

**Honest uncertainties.** For each epoch the table carries:

| column | what it is | median |
|---|---|---|
| `sig_photon` | photon noise of the good chunks, combined | 1.4 m/s |
| `sig_pyodine` | `pyodine`'s `rv_err`: the weighted chunk scatter / √N | 10.6 m/s |
| `sig_order` | std of the order means / √(orders) | 17.3 m/s |
| **`sig_epoch`** | **the adopted error**: max(`sig_pyodine`, `sig_order`), in quadrature with any excess from night pairs | **19.1 m/s** |

- The chunk-based error assumes independent chunks, and the order-to-order
  scatter says they are not. So the adopted error takes the larger of the
  two.
- **The within-night test.** Four pairs of frames lie closer than 2 h on
  one night, where the planet moves at most 6 m/s per hour. Their rms
  difference is 35 m/s, against 44 m/s expected from `sig_epoch`. There is
  no excess, so σ_extra = 0.
- **The orbit test (magnitude only).** Minus the modern catalogue's circular
  orbit (K = 68.3 m/s) with only an offset fitted, the residual rms is
  **26.2 m/s** and χ²/dof with `sig_epoch` is **1.15**. The errors describe
  the residuals.
- The epochs fitted over 2–8 orders have errors built from very few order
  means. They are the least certain rows of the table.

This is a magnitude check, taken because trap 8 says a plausible wrong answer
is the normal failure. It is not the assessment, which is prompt 8's. It says
the result is in the "doing well" range the goals set (below 30 m/s per epoch)
and not in the suspicious one (below 10).

**A check that shares no machinery with the velocities.**
- Phase 2's cross-correlation velocity minus the iodine velocity, epoch by
  epoch, correlates with the difference between PypeIt's ThAr zero point and
  the iodine-fitted one: r = **+0.982**, slope 0.91, residual 138 m/s
  (panel c).
- The two methods disagree by exactly what the ThAr wavelength zero point is
  wrong by: up to 2.55 km/s, on the borrowed-arc night 1998-07-19.
- This is phase 2's diagnosis confirmed from the other side. The iodine
  solution sits −2.6 to +0.95 km/s from PypeIt's.
- 07-19's three frames, which cross-correlation put at +1.5 to +2.4 km/s,
  come back at +36, +42 and +96 m/s. Prompt 9 examines this properly.

**Products.**
- **Code:** `first_hires_exoplanet/fit_all_epochs.py`. `fit_one_epoch.py`
  gained `run_driver(orders=...)` and an explicit template index in
  `chunk_table`, which would otherwise have repeated the index trap in the
  analysis.
- **Test:** `first_hires_exoplanet/tests/test_epoch_orders.py` (3 tests; 59
  pass in all).
- **Tables**, in `first_hires_exoplanet/data/`:
  - `iodine_velocities.csv`: 31 epochs, with `v_bary`, `sig_epoch` and
    every component, the three alternative estimators, BVC, counts, χ², the
    wavelength zero point against PypeIt, and phase 2's velocity;
  - `iodine_velocities_per_order.csv`: 367 epoch-order velocities;
  - `iodine_chunks.csv`: 18,350 run-1 chunks, with each chunk's time-series
    σ and offset (2.3 MB).
- **Figure:** `docs/figs/fig_p7_epochs.png`. Panel (d)'s ceiling at
  1000 m/s is `pyodine`'s own clamp on chunk σ (`sig_limit_up`).
- **Outside the repository:** `pyodine_runs/epochs/`, with each epoch's
  `pyodine` results and plots, the full per-chunk table
  (`all_epochs_chunks_run1.fits`) and the analysis log.

### Prompt 8: HD 187123 b detected at 10σ, from an open reduction

**In one paragraph.** On the same frames, the same catalogue and the same
fit that phase 2 used, the iodine forward model recovers the planet.
- **K = 72.4 ± 7.1 m/s**, a **10.2σ** detection, against the published
  72 m/s and the catalogue's 69.2 m/s. Phase 2 got 408 ± 188.
- The residual is **25.8 m/s per epoch**, and the adopted errors describe
  it: χ²/dof 1.05.
- A period search that is not told the period finds **3.0954 d**, with a
  false-alarm probability of 2 × 10⁻¹¹.
- Epoch by epoch the velocities follow Teklu et al.'s at r = +0.89, with
  slope 1.05.

This is the goal the document set: a convincing, independent detection of
the published Keplerian from a modern open reduction. It is **28× better
than phase 2** and **22× short of the modern catalogue** on the same photons.

`first_hires_exoplanet/figs_phase3.py` does the assessment:

```
conda run -n pypeit14 python -m first_hires_exoplanet.figs_phase3
```

**Comparable by construction.** The script imports `figs_phase2` and reuses
its pieces:
- the catalogue loader, the discovery-era cut and the time matching;
- the forced-period circular fit `fit_circular` (P = 3.0965828 d, linear, so
  there is no minimiser);
- the unweighted K error rescaled by the scatter;
- the signal-plus-noise correlation test, and the house style.

It feeds these prompt 7's `v_bary` and `sig_epoch` where phase 2 used
`v_rel` and `sigma_empirical`. Every phase-3 addition below is labelled as
one. Figure names are `fig_ph3_*`, because `fig_p3_*` belongs to prompt 3.

| figure | what it shows |
|---|---|
| `fig_ph3_timeseries.png` | ours and the catalogue over the nine months, **on the same scale** this time |
| `fig_ph3_phasefold.png` | both folded on 3.0966 d, with residuals: both now trace the Keplerian |
| `fig_ph3_compare.png` | epoch by epoch, equal aspect; phase 2 in grey, 22 of its 30 points off the frame |
| `fig_ph3_precision.png` | each stage against K, the 50 m/s target, phase 2 and the catalogue |
| `fig_ph3_residuals.png` | residuals in time, against the ThAr zero point and the orders fitted, and the period search |

Table: `first_hires_exoplanet/data/iodine_keplerian_residuals.csv`. It has
each epoch's velocity, error, the catalogue Keplerian, the residual, and the
catalogue's own velocity where one exists.

**1. The epochs match one to one.** All 30 discovery-era catalogue rows pair
with one of our epochs, worst time difference 0.1 min, as in phase 2. The
31st epoch, `HI.19980812.29160`, is the 60 s cell-in frame the catalogue
never had. We measure it at a residual of −6 ± 12 m/s. The iodine method
has made one velocity the modern catalogue does not contain.

**2. Against the modern catalogue, epoch by epoch** (30 matched):

| | phase 2 | **phase 3** |
|---|---|---|
| catalogue spread (rms) | 47.4 m/s | 47.4 m/s |
| ours (rms) | 723.5 m/s | **55.8 m/s** |
| difference (rms) | 706.7 m/s | **25.5 m/s** |
| correlation r | +0.38 | **+0.89** |
| r expected for the signal plus our noise | +0.07 | +0.88 |
| slope, ours on catalogue | — | 1.05 |

- Phase 2's r was "not evidence of anything". Phase 3's r = +0.89 rejects
  r = 0 at 7.4σ, and it is exactly what the signal plus our own measured
  noise predicts (+0.88).
- The difference, 25.5 m/s rms, is our error: at χ²/dof 1.03 against
  `sig_epoch`, the catalogue's 1.2 m/s being effectively truth. The honest
  error bar of prompt 7 is honest against an external standard.

**3. Against the published 3.097-day Keplerian.**

| fit | K (m/s) | residual rms | significance |
|---|---|---|---|
| Teklu et al. 2025 (phase 2's number) | 69.2 ± 0.3 | 2.2 | 211σ |
| phase 2, cross-correlation | 408 ± 188 | 668 | 2.2σ, not a detection |
| **phase 3, exactly as phase 2** (unweighted, error rescaled) | **72.4 ± 7.1** | **25.5** | **10.2σ** |
| the same, the 30 matched epochs only | 72.2 ± 7.3 | | |
| *phase 3 addition:* weighted with `sig_epoch` | 78.0 ± 4.9 (χ²/dof 0.99) | 25.8 | 15.9σ |
| *addition:* orbital phase fixed to the catalogue's | 77.8 ± 4.9 | | 15.8σ |

The phase-3 additions add up to three more results:
- **Constant against Keplerian** (weighted): χ² 279.2 → 27.6 for two extra
  parameters. F = 127.6, p = 8.5 × 10⁻¹⁵, **7.8σ**. The planet is required
  by the data, not just consistent with them.
- **The orbital phase agrees** with the catalogue's to −3.1° (± ~4°).
  Leaving out any one night moves K by 74–80 m/s (weighted); 1998-08-25's four
  frames move it most.
- **The period, not told:** Lomb–Scargle over 0.5–30 d peaks at
  **3.0954 d** (the true period is 3.0966 d), with power 0.90 and a
  false-alarm probability of **2.3 × 10⁻¹¹**. The next peak, 2.80 d (power
  0.72), is the ~30-day alias of the observing-run sampling:
  1/2.80 − 1/3.10 ≈ 1/30 d⁻¹.

**What the K numbers say:**
- Phase 2's method gives 72.4 ± 7.1 m/s: the published 72 m/s almost exactly,
  and +0.5σ from the catalogue's 69.2.
- The weighted fit, 78.0 ± 4.9, is +1.8σ above the catalogue. It leans on the
  best-measured epochs, and a 1.8σ excursion in one of two estimators of the
  same quantity is not a tension.
- The answer to report is phase 2's estimator, since that is what makes the
  two phases comparable. The weighted one is the cross-check.

**4. What the residuals contain.** About the catalogue's Keplerian, offset
only:

| | |
|---|---|
| rms | **25.8 m/s** (robust 22.3) |
| χ²/dof against `sig_epoch` | **1.05** |
| largest residuals | −74 (08-26, 6 orders, −1.2σ); +53 (07-19, 2 orders, +2.0σ); −50 (09-12, 7 orders, −1.2σ) |
| within a night (7 nights, 19 frames), about the nightly mean | 23.2 m/s |
| between nights (19 nightly means) | 19.8 m/s |
| period search on the residuals | best 0.94 d, FAP 0.62: nothing |

Six things follow:
- **Noise, as far as 31 epochs can tell.** χ²/dof is 1.05, no residual
  periodicity, no excess within a night over between nights. The largest
  residuals belong to the epochs with the fewest orders, and their errors say
  so.
- **The ThAr zero point is gone.** Spearman ρ against PypeIt's wavelength
  zero-point error is −0.07 (panel b). That error was 2.55 km/s on 07-19 and
  spans 3.5 km/s across the run, and it was the whole of phase 2's 700 m/s.
- **The borrowed-arc nights are fixed.**

  | night | phase 3 residuals | phase 2 velocities |
  |---|---|---|
  | 1998-07-19 | +5, +3, +53 m/s | +2438, +1609, +1554 m/s |
  | 1998-09-17 | +1 m/s | +516 m/s |

- **The intra-night drift is gone.** 1998-08-25 spanned 1604 m/s under
  cross-correlation; its four residuals are +6, +29, +22 and +25 m/s.
- **The model misfit does not reach the velocities.** Against χ²_ν, ρ =
  −0.01. So the misfit prompts 4–6 fought is harmless to relative velocities,
  as prompt 7 found.
- **One marginal signal: season.** The residuals correlate with the
  barycentric correction at ρ = +0.40 (p = 0.03), and with time (−0.34) and
  airmass (+0.34) at p ≈ 0.06.
  - The three quantities are close to one variable over June–September.
  - It is one of eight tests, and one p ≈ 0.03 among eight is what chance
    gives.
  - It is not a detection of anything. But it is the shape a
    barycentric-correction or template-epoch error would take, so prompt 9
    should look.
  - The multiplicative BVC term that prompt 7 neglected is ≤ 2 m/s here and
    cannot produce it.

Exposure time, photon noise and the number of orders show nothing:
|ρ| ≤ 0.07 on the residual, ≤ 0.26 on its magnitude.

**5. The achieved precision, and what it means.**

| | per epoch |
|---|---|
| photon noise | 1.4 m/s |
| adopted error (prompt 7) | 19 m/s |
| **achieved: residual about the Keplerian** | **25.8 m/s** |
| the goal: "phase 3 succeeds" | 50 m/s |
| phase 2, cross-correlation | 724 m/s |
| Teklu et al. 2025, the same frames | 1.2 m/s |
| Butler et al. 1998 | 3 m/s (as quoted in the goals) |

At 25.8 m/s with 31 epochs the goals' formula predicts σ_K = 6.6 m/s, an
11σ detection. The measured result is 7.1 m/s and 10.2σ, so the table in
`## Goals` was right.

**What this does and does not demonstrate.**

It **does** demonstrate an independent detection of HD 187123 b from the
1998 discovery frames, with open software end to end:
- PypeIt's reduction;
- `pyodine`'s forward model, forked here, with its changes recorded and
  tested;
- an iodine atlas from a different cell;
- no calibration machinery from the discovery pipeline.

The amplitude, the phase and the period all come out right, and the period
comes out without being given. Phase 2 said "nothing in this result would
have looked different if HD 187123 b did not exist". That is no longer true:
without the planet the velocities would scatter 26 m/s about a constant, and
they scatter 56.

It **does not** match the 1998 result or the modern one.
- **Against Butler's 3 m/s: about 9× short.**
- **Against the 1.2 m/s the modern pipeline gets from these frames: 22×
  short.**
- The photon noise is 1.4 m/s, so the photons are not the limit. The limits
  are known from prompts 5–7:
  - the Fischer atlas is a different cell: its line pattern misfits every
    chunk at χ²_ν ≈ 20, and only the chunk-offset combination removes most
    of that;
  - the red iodine orders (58–61) scatter 100–200 m/s epoch to epoch;
  - there is no measured instrumental profile beyond a per-chunk
    super-Gaussian.

Of those, only the first has an identified remedy (the HIRES cell's own
spectrum, option (c) of the Q&A), and it has not been tried. That should be
said wherever this result is.

The honest summary for the public document: **the planet is recovered at
10σ with per-epoch precision of 26 m/s. That is 28× better than a ThAr
wavelength solution allows, and still an order of magnitude short of what the
iodine method achieves with the machinery built for this instrument.**

### Prompt 9: what the forward model fixed, and what it did not

**In one paragraph.**
- **Everything the prompt named as a calibration defect is fixed.** The
  borrowed-arc nights, the intra-night drifts, the December nights without a
  pixel flat and 1998-09-13's broken arc all come back into line: none has a
  residual beyond 2σ, and each is at the level of an ordinary epoch.
- **The forward model also measures what was wrong.** PypeIt's zero point is
  off by up to 2.6 km/s, it drifts by up to 256 m/s per hour within a night,
  and the drift grows with distance in time from the arc.
- **What is left is not calibration.** It is the four red iodine orders
  (58–61), where the iodine is weak, and the six epochs with only 2–8 usable
  orders.
- **Leaving both out halves the per-epoch residual,** 25.8 → 11.5 m/s, and K
  stays at 71.8 ± 3.5 m/s.
- **Tellurics are a small contribution.** The seasonal correlation prompt 8
  flagged remains marginal and cannot be tested within a night.

`first_hires_exoplanet/diagnose_phase3.py`:

```
conda run -n pypeit14 python -m first_hires_exoplanet.diagnose_phase3
```

It reads prompt 7's velocities and chunks, prompt 8's residuals, phase 2's
cross-correlation and the reductions' own logs and calibrations. It writes
`data/diagnose_nights.csv` (one row per night: arc, flat, arc quality,
PypeIt's zero-point error, phase 2 and 3 residuals),
`data/diagnose_variants.csv` and `docs/figs/fig_p9_diagnosis.png`.

**A measurement phase 2 could not make.** The forward model fits each
chunk's wavelength solution from the iodine, so `zp = −dv_wave` is PypeIt's
wavelength zero-point error, measured on the science photons. Across the 31
epochs it spans **−944 to +2554 m/s** (median |zp| 308 m/s on own-arc
nights). Phase 2's residual follows it at Spearman **+0.94** (p = 8 × 10⁻¹⁵)
and phase 3's at +0.07 (panel a). **The whole of phase 2's error was PypeIt's
zero point, and phase 3 has none of it.** |zp| also grows with the time
between the science frame and the night's arc (Spearman +0.61, p < 0.001),
as flexure would.

**1. The borrowed-arc nights: back in line.**

| epoch | arc | PypeIt zp error | phase 2 residual | phase 3 residual |
|---|---|---|---|---|
| 07-19 24897 | 07-18's | +2554 m/s | +2132 m/s | +5 ± 13 (+0.4σ) |
| 07-19 35708 | 07-18's | +1878 | +1296 | +3 ± 15 (+0.2σ) |
| 07-19 49996 | 07-18's | +1373 | +1237 | +53 ± 27 (+2.0σ; 2 orders) |
| 09-17 29496 | 09-15's B2 | +620 | +290 | +1 ± 12 (+0.1σ) |

The own-arc epochs' residual rms is 25.3 m/s. The borrowed-arc nights are
indistinguishable from them, and 07-19's worst point is its 2-order last
frame, not its arc. Phase 1 concluded "same-night arcs are not required",
and phase 2 corrected that: "the RMS measures the scatter of the fit, not its
zero point". Phase 3 makes the point moot. With the iodine, the arc is only
a starting guess.

**2. The intra-night drifts: gone.** Slopes against UT, fitted within each
night:

| night | frames, span | phase 2 residual | PypeIt zp error | **phase 3 residual** |
|---|---|---|---|---|
| 08-25 | 4, 6.7 h | −240 ± 2 m/s/h (span 1630) | −256 ± 1 m/s/h | **+1 ± 3 m/s/h** (span 23) |
| 08-26 | 3, 6.9 h | −77 ± 6 (span 538) | −59 ± 16 | **−2 ± 4** (span 65) |
| 07-19 | 3, 7.0 h | −123 ± 74 (span 896) | −167 ± 27 | **+5 ± 4** (span 50) |
| 07-18 | 3, 5.6 h | +4 ± 2 (span 23) | −78 ± 37 | +4 ± 6 (span 38) |

- Phase 2's drift on 08-25 and 08-26 was PypeIt's wavelength zero point
  drifting through the night. The slopes agree to within the errors on 08-25
  (−240 against −256).
- Phase 3's residuals about each night's mean give χ²/dof **0.69** over 19
  frames on 7 nights: no drift, and the within-night errors are if anything
  generous.
- The 08-26 B stars show PypeIt's zero point moving by −32 ± 22 m/s/h,
  consistent with HD 187123's −59 ± 16 on the same night (panel b). But
  back-to-back B-star pairs 110 s apart differ by 18–121 m/s. With a
  slit-filling boxcar the B stars are a noisy probe of the wavelength zero
  point, and they are not needed for this answer.

**3. The December nights: nothing visible.**
- **1997-12-24's only frame is cell-out.** It can never have an iodine
  velocity.
- **1997-12-23** (trace from an iodine-in flat, no pixel or illumination
  flat) is ordinary on every quantity a flat-field error would move. Against
  the other 19 epochs with all 14 orders:

  | | 12-23 | others, median | percentile |
  |---|---|---|---|
  | χ²_ν (measured noise) | 23.9 | 22.9 | 95 |
  | per-chunk scatter | 370 m/s | 381 | 32 |
  | per-order scatter | 60 m/s | 55 | 74 |
  | good chunks | 561 | 543 | 89 |
  | **residual** | **−18 ± 23 m/s (−0.8σ)** | | |

- The χ² is 4% above the median. By order, 12-23's chunk χ² runs 0.87–1.39×
  the other epochs'; the largest excess is order 70 (1.39). That is the only
  trace of the missing pixel flat, and the model's per-chunk continuum
  absorbs it.
- The epoch is 180 days before every other one. It anchors the only long
  baseline, and it sits on the Keplerian.

**4. 1998-09-13: a bad arc, repaired in zero point but not in noise.**
- **Why 34 orders.** The night's arc (`HI.19980913.55970`, 15:32 UT)
  matches PypeIt's ThAr archive at median cc **0.65**, below 0.8 in every
  order, against 0.91 on every other night.
  - Orders 92 and 71 failed to reidentify, and order 86 fitted at 0.477 px.
  - PypeIt flagged all three `BADWVCALIB`, `BADTILTCALIB` and
    `BADFLATCALIB` and did not extract them: 37 traced, 35 solutions, 34
    extracted.
  - The flats are the same 1998-08-12 donors as every other night, and the
    edges agree with 09-14's to 1 pixel. Their wider traced widths (9.6–10.3
    px against 8.0–8.3) are untweaked edges on already-masked orders, a
    consequence and not the cause.
  - Why the arc itself is poor (lamp, exposure, a changed configuration)
    would need its raw header and counts; it is left for the PypeIt side of
    the prompt-10 report.
- **What the iodine sees.** PypeIt's wavelength solution within each order,
  as the chunk-to-chunk scatter of `zp`, is **1.3–3.7× the typical epoch's**
  in every order (panel c): the bad arc degrades the whole solution, not just
  three orders.
- **What the forward model makes of it:**
  - the zero point is repaired: the residual is **−35 ± 46 m/s (−0.8σ)**;
  - the noise is not: the epoch has the worst per-chunk (511 m/s) and
    per-order (114 m/s) scatter, higher than every one of the 20 epochs
    with all 14 orders, so its error is
    2.8× typical;
  - the red orders 58–61 deviate by −1.4 to −2.8 of their usual scatter.

  The iodine fit starts from PypeIt's wavelengths and smooths its run-0
  solution over each order, so a poor starting solution is not fully
  forgotten. The prompt-7 error bar reflects this honestly, since the epoch
  gets small weight.
- **Its chunk χ² is *lower* than typical** (0.6–0.9×), while its velocities
  scatter more. A loose wavelength solution lets the model fit the misfit
  away.

**5. What is left.**

*The red orders.*

| order | Å | iodine contrast | stellar contrast | photon σ per chunk | chunk time-series σ | tellurics in chunk |
|---|---|---|---|---|---|---|
| 58 | 6143 | 0.07 | 0.07 | 88 m/s | 402 | 0% |
| 59 | 6040 | 0.11 | 0.03 | 136 | 513 | 14% |
| 60 | 5940 | 0.16 | 0.04 | 96 | **819** | **97%** |
| 61 | 5845 | 0.19 | 0.04 | 101 | 606 | 33% |
| 62–70 | 5095–5750 | 0.25–0.36 | 0.06–0.18 | 25–60 | 169–350 | 0–37% |
| 71 | 5021 | 0.12 | 0.19 | 55 | 195 | 10% |

- The red orders have **weak iodine** (little to fix the wavelength with) and
  **weak stellar lines** (little to measure the velocity with). They
  repeat worst, and over the time series their order velocities scatter by
  97–198 m/s about the epoch velocity (prompt 7).
- The combination weights chunks by their own scatter, as if they were
  independent. But these errors are shared by every chunk in an order, and
  per-chunk weighting does not down-weight a shift common to a whole order.
- **Tellurics are a minor part.**
  - `pyodine`'s CARMENES mask (Wallace 2011, vacuum) puts telluric lines in
    97% of order 60's chunk-epochs and 33–37% of 61–62's. No telluric mask
    was used.
  - Chunks crossed by a line repeat slightly worse overall (269 against 243
    m/s, Mann–Whitney p = 0.05), but **no worse within orders 60–62**
    (322 against 311, p = 0.48).
  - Masking every crossed chunk in every epoch lowers the residual only from
    25.8 to 24.2 m/s.
  - Orders 58–59 have almost no telluric lines and are still bad.

  **The red orders are bad because the iodine is weak, not mainly because of
  the atmosphere.**

*The weak epochs.* Six epochs have 2–8 usable iodine orders:
- the last frames of 07-17, 07-18 and 07-19 (2 orders each);
- 07-18 38646 (8), 08-26 19932 (6) and 09-12 (7).

They are the ones phase 2's quality filter cut hardest, and their residuals
are the largest in absolute terms. Their errors are correspondingly large.

*What each costs.* `pyodine`'s combination re-run on subsets (re-combining
prompt 7 unchanged reproduces its velocities to 0.2 m/s). The error column
is prompt 7's `sig_epoch`, so a χ²/dof below 1 means those errors are now too
large.

| variant | epochs | residual rms | χ²/dof | K (phase-2 method) | significance |
|---|---|---|---|---|---|
| prompt 7, everything | 31 | 25.8 m/s | 1.07 | 72.4 ± 7.1 | 10.2σ |
| without the borrowed-arc, December and 09-13 epochs | 25 | 24.6 | 1.07 | 69.5 ± 8.0 | 8.7σ |
| telluric-crossed chunks masked | 31 | 24.2 | 0.94 | 73.0 ± 6.7 | 10.9σ |
| without the 6 epochs of ≤ 8 orders | 25 | 18.2 | 0.84 | 70.8 ± 5.6 | 12.6σ |
| tellurics masked, orders 58–59 dropped | 31 | 19.4 | 0.74 | 73.8 ± 5.3 | 13.9σ |
| **without orders 58–61** | 31 | **17.2** | 0.60 | 72.4 ± 4.7 | 15.4σ |
| **without orders 58–61 and the 6 weak epochs** | 25 | **11.5** | 0.43 | 71.8 ± 3.5 | 20.3σ |

- **Removing the nights the prompt asked about changes nothing** (24.6 m/s).
  They are fixed.
- **The red orders are the largest single item.** Dropping them alone takes
  a third off the residual.
- **K does not move** (69.5–73.8 m/s) under any variant. The detection does
  not depend on any of these choices.

A caution on the variants. Neither criterion used the Keplerian: the red
orders were singled out by their own time-series scatter (prompt 7) and the
weak epochs by phase 2's quality filter. But the decision to try them was
made after prompt 8's residuals were seen. The headline stays prompt 7's; see
the Q&A.

*The seasonal correlation* (prompt 8: ρ = +0.40 with BVC):
- The residual slope is +0.75 ± 0.45 m/s per km/s of BVC (1.7σ). Over the
  1998 season BVC and time are the same variable (Pearson r = −1.00), so no
  test across nights can tell them apart.
- **Within nights**, BVC changes by Earth's rotation (up to 0.66 km/s) while
  the season does not. Over 18 pairs the slope is −26 ± 15 m/s per km/s:
  consistent with zero, and with +0.75 at 1.8σ. That lever is 30× weaker
  than the seasonal one, so it can neither confirm nor exclude a
  BVC-proportional error.
- Without the 1997 epoch the slope is +0.70 ± 0.46.
- **Verdict: marginal, and not attributable.** The candidate mechanisms
  (the ≤ 2 m/s multiplicative BVC term; the template's wavelength scale;
  tellurics sweeping through chunks) are each too small on the numbers
  above. With 31 epochs this is not worth more.

**What the forward model fixed, and what it did not.**

| | fixed? | evidence |
|---|---|---|
| the ThAr wavelength zero point (phase 2's whole 700 m/s) | **yes** | residual vs PypeIt zp: ρ +0.94 → +0.07 |
| borrowed arcs (07-19, 09-17) | **yes** | 1.2–2.1 km/s → +0.1 to +2.0σ |
| intra-night drift (08-25, 08-26, 07-19) | **yes** | −240 m/s/h → +1 ± 3 m/s/h; within-night χ²/dof 0.69 |
| no pixel flat (1997-12-23) | **yes**, to 4% in χ² | −0.8σ, ordinary in every metric |
| a bad arc (1998-09-13) | **zero point yes, noise no** | −0.8σ, but the epoch's error is 2.8× typical |
| a missing extracted order (09-13) | by the caller, not the model | `usable_orders` (prompt 7) |
| the atlas's line pattern (a different cell) | **mostly**, by the chunk offsets | 290 → 230 m/s per chunk; 48.5 → 26 m/s per epoch |
| weak iodine in orders 58–61 | **no** | 97–198 m/s per order; worth 25.8 → 17.2 m/s |
| epochs with few usable orders | **no**: extraction, not model | worth 25.8 → 18.2 m/s |
| tellurics | not modelled; minor | worth 25.8 → 24.2 m/s |
| a seasonal/BVC term | undetermined | +0.75 ± 0.45 m/s per km/s |

**Products.**
- **Code:** `first_hires_exoplanet/diagnose_phase3.py`.
- **Tables:** `first_hires_exoplanet/data/diagnose_nights.csv` (20 nights)
  and `diagnose_variants.csv`.
- **Figure:** `docs/figs/fig_p9_diagnosis.png`:
  - (a) residual against PypeIt's zero-point error, phases 2 and 3;
  - (b) the 08-25 and 08-26 drifts;
  - (c) PypeIt's within-order solution as the iodine sees it, with 09-13 and
    12-23;
  - (d) residual against BVC.

### Prompt 10: the upstream reports, drafted and not sent

**Two drafts, not one.** The prompt names Ryan Cooke, and the PypeIt items
go to him. But half of what phase 3 found is in `pyodine`, whose authors are
Paul Heeren, René Tronsgaard and Frank Grundahl. Sending those items to
PypeIt would reach nobody who can act on them. So there are two drafts, both
in `docs/upstream/`, and **neither has been sent**:

| draft | to | items |
|---|---|---|
| `pypeit_report_draft.md` | Ryan Cooke | 15, plus an information section and three usability notes |
| `pyodine_report_draft.md` | the `pyodine` authors | 6 fork changes and 15 found-and-worked-around, plus one note on the combination |

**The PypeIt report.**
- **Items 1–9** are phase 1's. The fixes for 1–5 and the docstring in 9 are
  on `origin/orig-hires-fixes` at `017bece06`, which is `develop` at
  `f3a1f1d27` plus one commit touching only `keck_hires.py`. No PR has been
  opened.
- **Items 10–11** are phase 2's `core/wave.py` findings.
  - The solar-term sign error is **still in `develop`**, at `wave.py:128`
    (phase 2 cited the function at line 89).
  - Re-run over all 36 epochs, `refframe_audit.py` gives PypeIt's
    heliocentric correction +9.1 to +14.2 m/s off, drifting 5.16 m/s, and
    its barycentric correction −4.6 m/s off. These are phase 2's numbers.
- **Items 12–15 are new in phase 3:**
  - **12:** the inverse variance is 2.35× too large for the optimal
    extraction and 1.61× too small for the boxcar. The gain explains 1.59×
    of the first.
  - **13:** `OPT_WAVE` is exactly 0 at masked pixels. Checked here: 950 of
    950 zero wavelengths in one epoch fall on masked pixels.
  - **14:** boxcar flux in a partially masked aperture is biased low and not
    flagged. Checked here: 1,730 pixels (6%) have `BOX_NPIX` below 75% of
    the aperture, their flux sits at 0.51 of the continuum (r = 0.84 with the
    fraction summed), and `BOX_MASK` is True.
  - **15:** orders that fail calibration vanish from the `spec1d` with only
    an INFO line (1998-09-13).
- **Section E** gives, for information, what the iodine measures of PypeIt's
  ThAr zero point. **Section F** carries phase 2's three usability notes.

**Checks made while drafting.**
- `verify_orig_fixes.py` runs cleanly.
- `refframe_audit.py` had been globbing `reduce_1998*`. That now also caught
  the B-star reduction `reduce_19980826_lsf/` (11 frames of other stars), and
  it has always missed December 1997. It now globs `reduce_19` followed by
  exactly six digits: all 36 HD 187123 frames and nothing else.
  - `data/refframe_audit.csv` grows from the 10 July rows committed in phase 2
    to 36.
  - The drift numbers quoted in the report come back exactly (5.16 m/s
    peak-to-peak; −4.6 m/s barycentric).
  - The one-epoch breakdown in section 5 now uses 1997-12-23 as its first
    epoch (+13.33 m/s where the July epoch gave +13.09).

**The `pyodine` report.**
- **Changes 0–5** of the fork (`vendor/README.md`), each with its test.
- **Fifteen items worked around rather than changed:**
  - the `order_correction` gap;
  - the inverted Chauvenet condition;
  - `red_chi_sq` holding √χ²_ν;
  - the `jansson` adjoint;
  - the `misc.rebin` spline at HIRES sampling;
  - the `fit_lsfs` 10⁻¹² starts;
  - `NormalizedObservation` dropping `ivar`;
  - `create_template`'s bad-pixel mask;
  - the air-wavelength Arcturus reference;
  - the atlas mean normalisation;
  - swallowed exceptions, re-read in the source for the draft:
    `create_template` logs with `logging.error`, and only
    `model_single_observation` uses `logging.info`;
  - the plot-then-fork hang;
  - fixed `plot_chunks`;
  - `robust.mean` dividing by zero;
  - the mutated weighting default.
- It closes with prompt 9's point about order-level weighting.

### Prompt 11: the phase-3 assessment

**Phase 3 succeeded by the measure it set itself, and fell an order of
magnitude short of the 1998 result.**
- **The planet is detected** at 10.2σ from an independent, open reduction:
  K = 72.4 ± 7.1 m/s, against 72 m/s published (Butler et al. 1998) and
  69.2 m/s from the modern catalogue. The goal was a detection at 50 m/s per
  epoch.
- **The per-epoch precision is 25.8 m/s**, inside the "below 30 m/s we are
  doing well" line and outside the "below 10, be suspicious" one.
- **It is 28× better than phase 2** (716 m/s about the same orbit, same
  frames).
- **It is about 9× short of the 3 m/s the method was built for, and 12× short
  of the modern pipeline** on these photons (2.2 m/s about the orbit). The
  photon budget is 1.4 m/s.

Following the Q&A after prompt 9 (answer (a)), prompt 7's velocities are
phase 3's result. The prompt-9 variants (11.5 m/s) are reported as what a
better combination would reach, not as the result.

All numbers below are from committed products:
- `data/iodine_velocities.csv`;
- `data/iodine_keplerian_residuals.csv`;
- `data/diagnose_variants.csv` and `data/diagnose_nights.csv`;
- `data/lsf_models_summary.csv` and `data/bstar_alpha_payoff.csv`;
- the `figs.py` printout for the public figures.

#### What was achieved

| | |
|---|---|
| epochs | 31 cell-in frames, 19 nights, 1997-12-23 to 1998-09-18 |
| chunks fitted | 18,350 (13,950 good) |
| K, period fixed (phase 2's method) | **72.4 ± 7.1 m/s, 10.2σ** |
| K, weighted | 78.0 ± 4.9 m/s |
| period, not told | 3.0954 d (true 3.0966), FAP 2 × 10⁻¹¹ |
| orbital phase against the catalogue | −3° (± ~4°) |
| residual about the catalogue orbit | **25.8 m/s**, χ²/dof 1.05 against the adopted errors |
| ours − catalogue, epoch by epoch | 25.5 m/s rms, r = +0.89 |
| adopted per-epoch error | 19.1 m/s (median); honest against the catalogue (χ²/dof 1.03) |
| photon noise per epoch | 1.4 m/s |
| per-chunk scatter within an epoch | 290 m/s; 230 m/s with the chunk offsets removed |
| per-order scatter within an epoch | 55 m/s |

The chain the result rests on:
- **PypeIt** 2.0.2.dev1217+g017bece06 (branch `orig-hires-fixes`, `develop`
  `f3a1f1d27` plus one commit to `keck_hires.py`);
- **`pyodine`** at `4488b0914fe5b272b787982647691045bff2604a`, plus fork
  changes 0–5 (`vendor/README.md`);
- the **Fischer May 2022 atlas** (SHA-256 `3ae788e7…a55a429`) at α = 2.59
  with a free per-chunk depth;
- the **oversampling-1 deconvolved template** from the three 1998-08-26
  cell-out exposures;
- the **super-Gaussian** instrumental profile.

#### How it got there: what each prompt settled

| prompt | settled | what it turned out to matter |
|---|---|---|
| 1–2 | fork vendored; NumPy 2 + four HIRES changes, each tested | the air/vacuum and depth-model changes were prerequisites: 83 km/s and negative transmission otherwise |
| 3 | template deconvolved; the `misc.rebin` spline breaks it at HIRES sampling | oversampling 1; −4 to −9% on the scatter |
| 4 | one chunk fits; noise scale 0.426 measured | the model misfit (χ²_ν ≈ 20) is not noise |
| 5 | one epoch: right velocity, wrong iodine | misfit follows iodine contrast (partial ρ ≈ +0.7) |
| 6 | LSF: super-Gaussian, a ≤ 5% lever | the prompt expected the largest lever; it was not |
| 6 cont. | index trap fixed; B stars: α = 3.33; per-chunk α buys nothing | the atlas's *line pattern*, not its depth, is the misfit |
| 7 | all epochs; the misfit is common to every epoch | chunk offsets cancel it: 48.5 → 26 m/s |
| 8 | detection, 10σ, comparable to phase 2 | errors honest against the catalogue |
| 9 | every named calibration defect fixed | what is left is weak red iodine and weak epochs |
| 10 | upstream reports drafted | 15 PypeIt items, 21 `pyodine` items |

**The decision that mattered most** was the Q&A before prompt 7: run the
full time series before trying to repair the atlas. A misfit that looked
fatal on one epoch (χ²_ν ≈ 22, per-chunk scatter 6–7× photon noise) was
mostly a fixed pattern that `pyodine`'s own combination removes. Building an
empirical HIRES-cell spectrum first would have cost the most effort for the
part of the error that cancels anyway.

**The prediction that failed** was the document's own. It said the
instrumental profile was "probably the single biggest lever on the result".
Four LSF models tie to ±5%. The iodine atlas was the lever, and even that
mostly cancels.

#### Where phase 3 fell short of 1998, and why

In order of what each costs, from prompts 7 and 9:
1. **The atlas is a different cell.** It misfits every chunk at χ²_ν ≈ 20.
   The chunk offsets remove the fixed part (201 m/s spread). What varies
   between epochs leaves 230 m/s per chunk against ~45 m/s photon noise. Only
   the HIRES cell's own spectrum would remove this, and no FTS scan of it has
   been found.
2. **The weak-iodine red orders (58–61).** Their order velocities scatter by
   97–198 m/s over the time series. Chunk weights cannot see an error shared
   by a whole order. They are worth 25.8 → 17.2 m/s.
3. **Six epochs with 2–8 usable orders.** This is extraction quality inherited
   from phase 2's filter. They are worth 25.8 → 18.2 m/s.
4. **A simple instrumental profile** (a four-parameter super-Gaussian per
   chunk) against the multi-Gaussian IP the Lick/Keck pipeline refined over
   years. Among the models tried it is ≤ 5%. A better one has not been tried.
5. **Tellurics** (not modelled): 25.8 → 24.2 m/s.
6. **A marginal seasonal term:** +0.75 ± 0.45 m/s per km/s of BVC. It cannot
   be separated from time in this season.

Items 2 and 3 together: 11.5 m/s. Everything named in prompt 9's question
(borrowed arcs, drifts, flats, 09-13's arc) is fixed and costs nothing.

#### What phase 3 did not do, and what a phase 4 would need

- **A HIRES-cell spectrum.** An FTS scan of the Keck cell, if one exists, or
  the empirical B-star cell spectrum the Q&A after prompt 6 deferred.
  Expected to be the largest remaining gain.
- **Order-level weighting in the combination** (Q&A after prompt 9, option
  (b)). It is principled and Keplerian-free, and would recover most of item 2
  without hand selection. Any version run now has been informed by seeing
  the residuals.
- **The PypeIt noise model** (upstream item 12). The inverse variance is 2.35×
  too large (optimal) and 1.61× too small (boxcar). It leaves χ² values only
  comparable within one extraction.
- **The 1998-09-13 arc.** Why it is poor (cc 0.65) is not diagnosed.
- **The second planet.** HD 187123 c's ten-year orbit is invisible in nine
  months and was never in scope.

#### The fork, the tests and the upstream reports

- **The fork:** six changes (`vendor/README.md`), each with a test that fails
  without it.
- **Caller-side fixes** in `utilities_hires` and the phase-3 drivers: the
  zero-wavelength and zero-flux repairs, `proper_motion`, the boxcar
  extraction, `start_inside`, the echelle-matched order selection, and the
  plot-chunk trim.
- **Tests:** 59 pass.
- **Upstream reports** are in `docs/upstream/`. By the Q&A after prompt 10:
  - you send them;
  - the PypeIt items go as an email to Ryan Cooke, then GitHub issues and a
    PR from `orig-hires-fixes`;
  - the `pyodine` items go as GitHub issues on `pepeheeren/pyodine`. The
    draft is written as a letter and will need splitting into issues.

  Nothing has been sent from this side.

#### Did phase 3 meet its own success criterion?

"Phase 3 succeeds at 50 m/s per epoch": a convincing, independent detection
of the published Keplerian from a modern open reduction.
- **Precision:** 25.8 m/s per epoch, against the 50 m/s target. ✓
- **Detection:** 10.2σ at the published period, and 3.0954 d found without
  being told. ✓
- **Independent:** no code, velocities or calibration from the discovery
  pipeline. The catalogue enters only in the assessment. ✓
- **Open:** PypeIt and `pyodine`, both public, with every change recorded
  and tested. ✓
- **"Both of those facts should be stated plainly":** the result is 28× better
  than phase 2, and ~9× short of Butler (12× of the modern pipeline). It is
  stated so here and in `docs/public_HD187123b.md`. ✓

**Phase 3 is complete.**

#### The public document

`docs/public_HD187123b.md` gains two sections before the References:
- **"Going back to the 1998 spectra"** covers the three stages: the PypeIt
  reduction and its three detector bugs; the lamp-calibrated attempt and why
  it failed; the iodine forward model and the recovered planet;
- **"Where we fell short, and why"** says it plainly: 26 against 3 and
  2.2 m/s, with the reasons in plain language and the 12 m/s variant
  explicitly not claimed.

The closing paragraph of "What happened next" no longer says "we have not
done it yet". Two figures are added to `figs.py`:
- `fig7_three_ways.png`: the same frames folded three ways, the phase-2 panel
  on a 16× larger scale;
- `fig8_precision.png`: scatter about the orbit, 716 / 26 / 2.2 m/s against
  K.

Figures 1–6 were restored from `HEAD` after regeneration, because a
different Matplotlib renders them byte-differently with no change of
content.

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

### 2026-09-29 (Prompt 5: one epoch end to end through `pyodine`'s own driver)

**Task.** Fit one epoch end to end across all usable iodine orders. Report
the per-chunk velocity scatter, the chunk failures, the LSF against position
and the runtime. Compare the fitted continuum and wavelength solution with
PypeIt's. Findings in `## Report`, "Prompt 5: one epoch, end to end — right
velocity, wrong iodine". There were no new Q&A answers since prompt 4, and
the four decisions from after prompt 3 still apply: both templates were
fitted, and the weights stayed PypeIt's.

**What was done.**

- Chose `HI.19980825.19425`, the best cell-in epoch with all 14 iodine orders
  usable. This is forced by an index trap found while reading
  `auto_wave_comoving`.
- Ran `pyodine_model_observations.model_single_observation` unchanged
  through `utilities_hires`, against both prompt-3 templates.
- The first run fitted every chunk and then failed to save: the adapter's
  `Star` lacked `proper_motion`, and the driver said so only in its info log.
  Fixed test-first (`tests/test_adapter_results_io.py`, two tests).
- Wrote `first_hires_exoplanet/fit_one_epoch.py`, which:
  - runs the driver and treats an empty error log as proof of nothing;
  - builds a per-chunk table: every parameter and error, χ² on both noise
    scales, photon-noise velocity error, the wavelength offset against
    PypeIt, the de-normalised continuum against `OPT_FLAT`, the evaluated
    LSF's FWHM/centroid/asymmetry, and the iodine and stellar contrast;
  - flags five kinds of failure;
  - bins the LSF against the same night's ThAr;
  - correlates the misfit with iodine and stellar structure.
- Ran it twice in full. The numbers reproduce exactly.
- Added a second Q&A round for prompt 6.

**Headline results.**
- 700 chunks per run, in 48–59 s per template for both runs.
- No lmfit failures, but ~8.5% of chunks without uncertainties. About 80%
  are good by all five tests.
- Epoch median velocity −677 m/s (os 10, run 1) against −689 expected, with
  an error of ~13 m/s.
- Per-chunk robust σ 287–340 m/s, 6–7× the ~40–48 m/s photon noise.
  χ²_ν ≈ 23 against the measured noise.
- The misfit tracks iodine contrast (partial ρ ≈ +0.7), not stellar (+0.1).
- Stellar LSF ≈ 2.0 px, varying 1.6–2.6 px with order and position. The same
  night's ThAr is 2.51 px, where 08-26's was 2.18.
- Fitted wavelengths sit −0.95 km/s from PypeIt's, with 130–150 m/s
  order-to-order scatter and order-edge structure.
- The continuum follows PypeIt's blaze to 2.4% chunk to chunk, after
  de-normalisation.

**What this taught us about the repository and the data.**

*The driver's chunk indexing is a silent trap for this data set.*
`auto_wave_comoving` builds chunks only for the orders it is given, while
the model fetches template chunk *i* for observation chunk *i*. On any subset
of orders every later chunk is modelled against the wrong stretch of
spectrum, with no error, because the template spline extrapolates. 11 of the
31 cell-in epochs have a rejected iodine order. This must be settled before
prompt 7 (Q&A after prompt 5, question 2).

*`pyodine`'s error log is not where its errors go.*
`model_single_observation` reports its own crash with `logging.info`. The
first run crashed after fitting every chunk and left an empty `errors.log`.
Any batch driver in prompt 7 must check for result files and scan the info
log.

*The adapter's contract gaps surface one driver step at a time.* Zero
wavelengths broke the velocity guess (prompt 3), zero-flux runs broke the
chunking (prompt 4), and a missing `proper_motion` broke the saving (here).
Each was invisible to the step before. The adapter now has an end-to-end
test for every stage the driver reaches, which is the only kind that has
caught any of them.

*The fit is limited by the iodine model, and one statistic shows it.* The
per-chunk misfit correlates with how much iodine structure a chunk carries,
with the stellar structure held fixed, and hardly at all the other way
round. That separates the three candidates prompt 4 could not:
- an LSF error would misfit both components;
- a template error would follow the star;
- an atlas that is the wrong cell follows the iodine, and that is what the
  data show.

It changes what prompt 6 is for.

*The ThAr cannot stand in for the stellar profile, in either direction.* On
08-26 the arc (2.18 px) was broader than the star (~2.05). On 08-25 the arc
is 2.51 px and the star is still ~2.0. The star's LSF is steadier than the
arc's, which makes sense if the star under-fills the slit and the arc does
not. The B stars on disk are the right measurement.

*The half-km/s wavelength offset from prompt 4 is not a constant.* It is
−0.95 km/s here against −0.55 on 08-26. It is structured along the band and
at the order ends, where a per-order polynomial is weakest. That fits
phase 2's diagnosis, that PypeIt's ThAr zero point is what failed, and it
does not fit a fixed error in the atlas. It is two epochs, so it is a lead,
not a result.

*The epoch velocity is already good; the chunk velocities are not.* A 13 m/s
epoch error from chunks that each scatter by 300 m/s holds only if the chunk
errors are independent. The order-to-order scatter (64–154 m/s) is larger
than within-order scatter would predict, so they are not fully independent.
Prompt 7's epoch-to-epoch scatter is the test that matters.

*The template's oversampling barely matters to the forward model.* The
re-convolution failure that dominated prompt 3 moves the epoch statistics by
a few percent. It matters less than the iodine model by an order of
magnitude.

**Files added / modified.**
- Added:
  - `first_hires_exoplanet/fit_one_epoch.py`;
  - `first_hires_exoplanet/tests/test_adapter_results_io.py`;
  - `first_hires_exoplanet/data/one_epoch_chunks.csv`, `one_epoch_summary.csv`
    and `one_epoch_lsf.csv`;
  - `docs/figs/fig_p5_one_epoch.png`.
- Modified: `first_hires_exoplanet/utilities_hires/load_pyodine.py`
  (`Star.proper_motion`).
- Outside the repository:
  `../first-hires-exoplanet-data/pyodine_runs/HI.19980825.19425_os10/` and
  `_os1/` (`pyodine`'s result files, logs and diagnostic plots).
- Nothing in `vendor/pyodine/` was changed.

**Whether any git command changed state.** No. Only read-only `git status`
and `git log`.

### 2026-09-29 (Prompt 6: five LSF models compared; the super-Gaussian chosen; the LSF is a ≤5% lever)

**Task.** Settle the instrumental profile: compare `pyodine`'s LSF models on
the prompt-5 epoch, choose one on evidence, and report how much the
per-chunk scatter depends on the choice. Findings in `## Report`, "Prompt 6:
the LSF is settled, and it is not the lever". The three Q&A questions after
prompt 5 were unanswered, so this ran as written (question 1, option a).
Nothing was reduced and the fork was not changed.

**What was done.**
- Wrote `first_hires_exoplanet/compare_lsf_models.py`, which runs prompt 5's
  driver (`fit_one_epoch.run_driver`) once per run-1 model and template, 10
  runs:
  - run 0 is a single Gaussian;
  - run 1 is one of single, super, Lick multi-Gaussian, SONG multi-Gaussian
    or Hermite;
  - run 2 smooths run 1's LSF (±160 px, ±3 orders, order separation 71
    physical px measured from the traces) and holds it fixed.
- The first pass returned the super-Gaussian and Hermite fits unmoved from
  their start. Diagnosed from the saved fit reports, fixed with a start rule,
  and rerun.
- Added paired bootstraps against the prompt-5 configuration and between
  templates.
- Settled the choice in `utilities_hires/pyodine_parameters.py`:
  - run 1 is the super-Gaussian, with physical bounds and `start_inside`;
  - the same in `Template_Parameters`;
  - `osample_temp = 1`.
- Added `tests/test_lsf_settled.py` (seven cases). Made `misfit_drivers`
  tolerate an empty set. Made a `--skip-run` pass keep the last runtimes.
- Updated `fit_one_epoch.py`'s docstring. Added a note under the open
  question 1.

**Headline results.**
- **The four workable models are indistinguishable:** robust per-chunk σ
  278–319 m/s, and every paired change against the Lick multi-Gaussian
  within ±13 m/s with a 68% interval including zero.
- **Hermite is worse (+92 to +108 m/s) and biases the epoch by +130 m/s.**
- **Run 2's smoothed fixed LSF never helps.** It costs good chunks and can
  bias the epoch by 3–5σ.
- **The oversampling-1 template lowers scatter** by 11–28 m/s for the same
  four models.
- **Chosen:** the super-Gaussian (tied-best scatter, lowest χ²_ν, centred, 4
  parameters) with the oversampling-1 template.
- **Every model finds the same profile:** FWHM 2.0–2.05 px, symmetric.
- **The iodine correlation of the misfit is unchanged by any of them.**

**What this taught us about the repository and the data.**

*The prompt's premise was wrong, and the data say so cleanly.* The LSF was
expected to be the largest lever; among the models that work it is under
5%. That is not because the profile does not matter. The models all find
the same profile (FWHM ~2.0 px, symmetric), and the residual that limits
the chunks is elsewhere. With prompt 5's partial correlation, the ranking
now reads: the iodine model first, then everything else together at the
~10% level.

*A flexible model started at zero is a rigid one.* `fit_lsfs` hands run 1 a
start in which every non-Gaussian shape parameter is ~10⁻¹², and `lmfit`'s
leastsq steps relative to the value, so those parameters cannot move. When
they cannot, nothing moves: the super-Gaussian's first attempt returned
every chunk at run 0's median velocity, a robust σ of *exactly zero*, which
looks like the best result in the table unless one asks whether it is
physically possible. Hermite returned "uncertainties could not be
estimated" 700 times. SONG's own configuration, which bounds its Hermite
weights at ±2×10⁻¹³, looks like someone hitting this and pinning the
parameters rather than fixing the start. The start rule now lives in
`pyodine_parameters` and is tested.

*Flexibility the data cannot support buys bias, not precision.* The Hermite
model is the most flexible that converges. It fits worse (χ²_ν 33–35), finds
an asymmetry no other model sees, and moves the epoch 130 m/s. An LSF shape
term is degenerate with velocity, and at χ²_ν ≈ 23 the fit uses every
degree of freedom it is given to chase the iodine misfit.

*Per-chunk LSF freedom is absorbing something, and smoothing it away shows
this.* Holding a smoothed LSF fixed is the Butler/Lick practice, and here it
raises χ², loses chunks and biases the epoch. On well-modelled data,
smoothing an LSF only removes noise. Here it removes a per-chunk degree of
freedom that was compensating for the atlas. That is one more sign of an
iodine-model problem, and a reason not to adopt run 2 until the atlas is
addressed.

*Model-selection statistics are not evidence at χ²_ν ≈ 23.* BIC ranks the
10-parameter Lick multi-Gaussian first, but its velocities are no better.
With systematic misfit, a likelihood criterion measures flexibility. The
decision was taken on paired velocity scatter, which is what an RV depends
on.

*The template question is closed.* Across four models, paired, oversampling
1 beats 10 by 11–28 m/s of per-chunk scatter. Prompt 3's spline finding
shows up in the velocities, modestly.

*Analysis-only reruns must not overwrite what they cannot recompute.* The
first `--skip-run` pass wrote a summary with blank runtimes over the one
with real ones. It now carries runtimes forward.

**Files added / modified.**
- Added:
  - `first_hires_exoplanet/compare_lsf_models.py`;
  - `first_hires_exoplanet/tests/test_lsf_settled.py`;
  - `first_hires_exoplanet/data/lsf_models_summary.csv`,
    `lsf_models_chunks.csv` (1.9 MB) and `lsf_models_template_pairs.csv`;
  - `docs/figs/fig_p6_lsf_models.png`.
- Modified:
  - `first_hires_exoplanet/utilities_hires/pyodine_parameters.py` (settled
    LSF, `start_inside`, `lsf_bounds`, `SUPER_BOUNDS`, `osample_temp = 1`);
  - `first_hires_exoplanet/fit_one_epoch.py` (empty-set guard, docstring);
  - this document's Q&A (a note under question 1 after prompt 5).
- Outside the repository: `../first-hires-exoplanet-data/pyodine_runs/lsf_*`
  (10 run directories).
- Nothing in `vendor/pyodine/` was changed.

**Whether any git command changed state.** No. Only read-only `git status`,
`git log` and `git diff --stat`.

### 2026-09-29 (Prompt 6, continued: index trap fixed in the fork; B stars reduced and fitted; per-chunk α buys nothing)

**Task.** Act on the three Q&A decisions after prompt 5:
1. fit the B stars with a per-chunk α, then compare LSF models;
2. fix the index trap in the fork;
3. reduce the B stars.

Findings in `## Report`, "Prompt 6, continued".

**What was done.**
- **Fork change 5, test first.** `auto_wave_comoving` records
  `chunk.template_index`, and the model and plots resolve the template
  through `template_chunk_index`. `tests/test_fork_chunk_index.py`;
  `vendor/README.md` now lists all six fork changes.
- **Wrote `reduce_bstars.py`.** It stages the night's arc and donor flats
  with the B-star frames, reuses the night's parameters and calibrations,
  and adds `force_center_obj = True`.
  - The first run stopped at HR 8634 (37323), found in 1 of 37 orders. I
    began a frame-by-frame workaround; you pointed to the Shane/Hamspec fix
    (`force_center_obj`); with it all 11 frames extract.
  - Added `extraction='BOX'` to the adapter (`tests/test_adapter_boxcar.py`),
    because the optimal mask collapses on slit-filling stars.
- **Wrote `fit_bstars.py`:** `pyodine`'s hot-star path on the template's
  chunk grid, with `BadPixelMask` on each frame.
  - Fitted 10 frames (one skipped: the star missed the slit).
  - Measured the frames' noise from the HR 8634 pair.
  - Built the α map.
  - Compared the five LSF models on the HR 8634 pair.
  - Refitted the prompt-5 epoch with `iod_depth` fixed per chunk to the
    B-star map, against a free control.
- Added a Q&A round before prompt 7. 56 tests pass.

**Headline results.**
- **The fork fix** passes its before/after test; any subset of orders is now
  modelled against the right template chunks.
- **The HIRES cell is α ≈ 3.33** (per frame 3.29–3.45), with a gentle red-ward
  rise and no chunk-to-chunk structure beyond the errors.
- **Pure-iodine fits** reach χ²_ν 1.0–5 against the measured (boxcar) noise,
  but still misfit more where the iodine is stronger (+0.40).
- **B-star LSF** 2.1 px. The LSF models are near-equivalent on pure iodine,
  Hermite unstable.
- **The payoff is negative.** Fixing α per chunk on the epoch: robust σ −7
  [−26, +7] m/s (n.s.), χ² 22.6 → 23.4, iodine partial ρ +0.69 → +0.67.
- **Iodine against PypeIt wavelengths** moves by half a km/s between B-star
  frames through the night.

**What this taught us about the repository and the data.**

*The limiting error is the atlas's line pattern, and three tests now agree.*
- Prompt 5: the misfit follows iodine structure, with the star held fixed.
- Prompt 6: no LSF model changes that.
- Here: no depth model does either, and pure iodine shows it too.

What is left is that the Fischer atlas is a different cell (phase 2: 0.928
line correlation). No parameter of this model can make up for that, which is
why the new Q&A asks whether to build the cell's own spectrum before prompt 7
or let prompt 7 measure what it costs.

*The cell is thicker than phase 2 measured.* α ≈ 3.33 from ten B-star frames
against 2.59 from two 60 s frames. The per-chunk `iod_depth` has absorbed the
difference in every fit so far (it is free), which is why fixing it did not
matter; the adapter's 2.59 stays as the starting point.

*A bright star in a short slit needs `force_center_obj`, and then a boxcar.*
Peak-finding fails on a slit-filling star, which is the same failure mode as
Shane/Hamspec, and the parameter built for it works unchanged on
`keck_hires_orig`. Then the optimal extraction's mask collapses on these
stars, as it did on HD 187123 in poor seeing (phase 2), and the boxcar is the
product to use. The boxcar does not mask the bad columns, so a bad-pixel mask
has to. Without that the fit was ruined: χ²_ν ≈ 1800, falling to 44 with the
mask.

*Every noise scale belongs to one extraction.* Pair differences put the
boxcar's inverse variance 1.61× too *small* and the optimal one's 2.35× too
large. Any χ² comparison across extractions, B stars against HD 187123 for
instance, has to be made in one set of units or not at all.

*Workarounds should wait for the known fix.* I had started reducing the
frames one by one to route around the extraction failure. The remedy already
existed in PypeIt for exactly this case; it cost one parameter and gave all
11 frames.

*The wavelength zero point moves within a night.* The iodine-derived zero
point against PypeIt's single-arc solution varies by ~0.5 km/s between B-star
frames an hour apart. That is prompt 9's intra-night drift question, and one
more reason the forward model's own wavelength solution is the point of
phase 3.

**Files added / modified.**
- Fork: `vendor/pyodine/pyodine/chunks.py`, `models/spectrum.py`,
  `plot_lib.py`; `vendor/README.md`.
- Added:
  - `first_hires_exoplanet/reduce_bstars.py` and `fit_bstars.py`;
  - `first_hires_exoplanet/tests/test_fork_chunk_index.py` and
    `test_adapter_boxcar.py`;
  - `first_hires_exoplanet/data/bstar_chunks.csv` (2.0 MB),
    `bstar_alpha_map.csv`, `bstar_lsf_models.csv` and
    `bstar_alpha_payoff.csv`;
  - `docs/figs/fig_p6b_bstars.png`.
- Modified: `first_hires_exoplanet/utilities_hires/load_pyodine.py`
  (`extraction`); this document (Report, Q&A).
- Outside the repository: `redux/stage_`, `setup_` and
  `reduce_19980826_lsf/`; `pyodine_runs/bstar_alpha_{free,fixed}_*`.

**Whether any git command changed state.** No. Only read-only `git status`
and `git log`.

### 2026-09-30 (Prompt 7: all 31 cell-in epochs fitted; 55 m/s epoch to epoch, 19 m/s honest errors)

**Task.** Fit all 31 cell-in epochs and produce a velocity table with honest
uncertainties, reporting the per-chunk, per-order and epoch-to-epoch scatter
separately, as phase 2 prompt 5 did. The Q&A before this prompt was answered
(a): run it with the current model. Findings are in `## Report`, "Prompt 7".

**What was done.**
- Wrote `first_hires_exoplanet/fit_all_epochs.py`:
  - an epoch list (phase 2's reduced frames with `IODIN = T`);
  - per-epoch usable orders, matched to the template by echelle number;
  - `pyodine`'s driver run in three shards;
  - per-chunk tables mapped to template chunks;
  - `pyodine`'s iSONG chunk combination;
  - four error estimates, a within-night pair test, a magnitude check
    against the catalogue orbit, and a cross-check against phase 2;
  - the tables and `fig_p7_epochs.png`.
- `fit_one_epoch.py`: `run_driver` passes `orders` through, and
  `chunk_table` takes an explicit template index.
- Added `tests/test_epoch_orders.py`.
- Fitted every epoch. Two epochs were refitted after the traps below.
  Ran the analysis, and all 59 tests.

**Headline results.**
- **Per chunk:** 290 m/s robust within an epoch, 230 m/s after the chunk
  offsets are removed. The offsets spread 201 m/s and are common to every
  epoch.
- **Per order:** 55 m/s within an epoch. Over the time series, 31–42 m/s in
  orders 64–70 and 97–198 m/s in orders 58–61.
- **Epoch to epoch:** 55.3 m/s rms, planet included, against 733 m/s from
  cross-correlation.
- **Adopted error:** 19.1 m/s median. Close night pairs show no excess, and
  the residual about the catalogue orbit is 26.2 m/s at χ²/dof 1.15. This is
  a magnitude check only; prompt 8 is the assessment.
- **Cross-check:** cross-correlation minus iodine tracks the ThAr-against-
  iodine zero point at r = 0.982.

**What this taught us about the repository and the data.**

*The misfit was never the whole story; its variation is.* Prompts 4–6 fought
a χ²_ν ≈ 22 misfit that no LSF, α or template sampling could remove. Most of
it is a fixed pattern per chunk, the same in every epoch, and `pyodine`'s
chunk-offset combination removes it. Judging per-chunk scatter on a single
epoch overstated what the misfit costs. The median-to-weighted step (48.5 →
26.2 m/s about the orbit) is the size of that overstatement.

*`pyodine` assumes contiguous orders, a second time.* The driver maps
template to observation orders with one constant shift. An extraction with a
missing order (1998-09-13: 34 orders, no 71) is silently bridged, and fifty
chunks are fitted against the wrong order with velocities of up to 7.7 × 10⁶
m/s. It is the index trap one level up, and it goes in the prompt-10 report.

*The analysis code can repeat a trap the fork has fixed.* `chunk_table`
looked up the template by list position. It was harmless while every epoch
had all 14 orders, and wrong from the first subset. Fixing a trap in the
fork does not fix every copy of the assumption.

*`pyodine`'s combination is fragile at the edges.*
- `robust.mean` divides by zero on a one-epoch column.
- The default weighting dictionary is mutated by its first implicit use.
- The example-chunk plots index by fixed chunk numbers and raise
  `IndexError` on a small epoch.

All three are upstream items.

*Waiting on a process by `pgrep -f` of its own command line matches the
waiter.* A background wait loop never ended because its own shell contained
the pattern. Check the actual worker processes.

**Files added / modified.**
- Added:
  - `first_hires_exoplanet/fit_all_epochs.py`;
  - `first_hires_exoplanet/tests/test_epoch_orders.py`;
  - `first_hires_exoplanet/data/iodine_velocities.csv`,
    `iodine_velocities_per_order.csv` and `iodine_chunks.csv` (2.3 MB);
  - `docs/figs/fig_p7_epochs.png`.
- Modified: `first_hires_exoplanet/fit_one_epoch.py`; this document
  (Report, Logs).
- Outside the repository: `../first-hires-exoplanet-data/pyodine_runs/epochs/`
  (31 epoch directories, the full chunk table, logs).

**Whether any git command changed state.** No. None was run.

### 2026-09-30 (Prompt 8: the planet is recovered, K = 72.4 ± 7.1 m/s at 10σ, 26 m/s per epoch)

**Task.** Assess prompt 7's velocities against the published 3.097-day
Keplerian and the modern catalogue, exactly as phase 2 prompt 6 did, and
state the precision, the significance and what the residuals contain, with
figures from a script. No new Q&A answers since prompt 7. Findings are in
`## Report`, "Prompt 8".

**What was done.**
- Wrote `first_hires_exoplanet/figs_phase3.py`. It imports `figs_phase2`
  and reuses its catalogue loader, era cut, matching, `fit_circular`,
  correlation test and style, so the two phases are compared by the same
  code.
- It adds a weighted fit, a phase-fixed amplitude, an F-test, a
  Lomb–Scargle period search, leave-one-night-out, and residual
  correlations with eight epoch quantities.
- It writes five figures (`docs/figs/fig_ph3_*.png`) and
  `data/iodine_keplerian_residuals.csv`.

**Headline results.**
- **Phase 2's method:** K = 72.4 ± 7.1 m/s (10.2σ), against 408 ± 188 in
  phase 2.
- **Weighted:** 78.0 ± 4.9 m/s.
- **Against the catalogue:** r = +0.89, and the difference is 25.5 m/s at
  χ²/dof 1.03.
- **Period, not told:** 3.0954 d, FAP 2 × 10⁻¹¹.
- **Residuals:** 25.8 m/s, χ²/dof 1.05, with no ThAr zero-point, misfit or
  order-count dependence. There is a marginal seasonal correlation
  (ρ = +0.40 with BVC, p = 0.03, one of eight tests).
- **Precision:** 28× phase 2, 1.9× inside the 50 m/s goal, 22× short of the
  catalogue.

**What this taught us about the repository and the data.**

*Prompt 7's errors are honest against an external standard.* Its adopted
error was built from internal scatter and a four-pair night test. Against
the catalogue's 1.2 m/s velocities it gives χ²/dof 1.03. So the error model
needs no fudge factor, and the result can be quoted with its own error bars.

*The goals' precision table predicted the outcome.* σ_K = σ√(2/N) gave
6.6 m/s at 25.8 m/s per epoch, and the fit returned 7.1. The planning
arithmetic in `## Goals` can be trusted for the "what would it take"
questions of prompt 11.

*Reusing phase 2's code is what makes the comparison valid.* Importing
`figs_phase2` rather than copying it guarantees the same matching, era cut
and K-error convention. A rewritten version would have had to be checked
line by line against the old one. The one place the phases must differ, the
error columns, is a two-line adapter (`load_ours`).

*Two estimators of K differ by 1.3σ of their own error* (72.4 unweighted,
78.0 weighted). Report the one that makes the phases comparable, and show
the other.

*Residuals against time, BVC and airmass are nearly one variable over a
June–September season.* A correlation with any of them cannot say which one
is the cause. Prompt 9 needs a test that separates them, for example the
intra-night BVC change on the long nights.

**Files added / modified.**
- Added: `first_hires_exoplanet/figs_phase3.py`;
  `first_hires_exoplanet/data/iodine_keplerian_residuals.csv`; `docs/figs/`
  `fig_ph3_timeseries.png`, `fig_ph3_phasefold.png`, `fig_ph3_compare.png`,
  `fig_ph3_precision.png` and `fig_ph3_residuals.png`.
- Modified: this document (Report, Logs).

**Whether any git command changed state.** No. Only read-only `git status`.

### 2026-09-30 (Prompt 9: every named calibration defect fixed; what is left is weak red iodine and weak epochs)

**Task.** Diagnose what is left:
- the borrowed-arc nights;
- the 08-25 and 08-26 drifts;
- the December nights without a pixel flat;
- 1998-09-13;
- what the forward model fixed and what it did not.

No new Q&A answers since prompt 8. Findings are in `## Report`, "Prompt 9".
A new Q&A question is posed ("After prompt 9").

**What was done.**
- Wrote `first_hires_exoplanet/diagnose_phase3.py`. It covers:
  - per-night arc quality, parsed from `run_pypeit.log` (the ThAr archive
    cross-correlation) and `run_summary.json`;
  - PypeIt's wavelength zero-point error per epoch and per order, as the
    iodine fit measures it;
  - within-night slopes of phase 2's residual, PypeIt's error and phase 3's
    residual, with the 08-26 B stars;
  - December and 09-13 against the full-format epochs;
  - 09-13's slit and wavelength calibrations (`SlitTraceSet`, `WaveCalib`,
    the slit bitmask);
  - a telluric test with `pyodine`'s CARMENES mask;
  - the seasonal term within and between nights;
  - `pyodine`'s combination re-run on seven subsets.
- It writes two tables and one figure.

**Headline results.**
- **Borrowed arcs:** 07-19 and 09-17 come back to +0.1 to +2.0σ from
  +0.3 to +2.1 km/s.
- **Drifts:** 08-25's −240 m/s/h becomes +1 ± 3 m/s/h.
- **December:** 1997-12-23 is ordinary (−0.8σ; χ² +4%); 12-24 is cell-out.
- **1998-09-13:** a bad arc (cc 0.65 against 0.91) cost three orders in
  PypeIt. The iodine repairs its zero point (−0.8σ) but not its noise (2.8×).
- **Phase 2's error:** residual against PypeIt's zero-point error is
  ρ = +0.94 in phase 2, +0.07 in phase 3.
- **What is left:** orders 58–61 (weak iodine), worth 25.8 → 17.2 m/s, and
  the six ≤ 8-order epochs, worth 25.8 → 18.2 m/s. Together 11.5 m/s, with
  K unchanged at 71.8 ± 3.5. Tellurics are worth 25.8 → 24.2. The seasonal
  term is +0.75 ± 0.45 m/s per km/s, undetermined.

**What this taught us about the repository and the data.**

*The forward model is also a calibration diagnostic.* Its per-chunk
wavelength fit measures PypeIt's zero-point error on the science photons:
2.6 km/s on a borrowed-arc night, 256 m/s per hour of drift, and growth with
time from the arc (ρ = 0.61). That is the most direct statement yet of what
phase 2 could only infer.

*1998-09-13's anomaly is the arc, not the flats or the tracing.* The
PypeIt log has the evidence. Every order cross-correlates against the archive
at cc 0.4–0.7 against ~0.9 elsewhere, and three orders were flagged
`BADWVCALIB|BADTILTCALIB|BADFLATCALIB` and dropped from extraction. Phase 2
carried it undiagnosed for four prompts; one grep of `cc =` in the log
settles it.

*Chunk weighting does not see order-level errors.* `pyodine`'s combination
treats chunks as independent. Where a whole order shifts coherently, as in
the red orders, it keeps weight it should not have. Order-level weighting is
the natural next step (Q&A option (b)).

*Tellurics are a real but small effect here.* No telluric mask was set,
and the iodine band is mostly clear of lines. Where lines exist (order 60)
the order is bad for other reasons too. Masking them buys 1.6 m/s.

*Season and BVC are one variable for this data set.* A within-night test
has 30× less lever. Do not chase the prompt-8 correlation with 31 epochs.

*Selection after seeing the answer is still selection.* Both subsets
chosen here have Keplerian-free criteria. But they were tried because of the
residuals, so they are reported as variants, and the headline stays
prompt 7's until the Q&A says otherwise.

**Files added / modified.**
- Added: `first_hires_exoplanet/diagnose_phase3.py`;
  `first_hires_exoplanet/data/diagnose_nights.csv` and
  `diagnose_variants.csv`; `docs/figs/fig_p9_diagnosis.png`.
- Modified: this document (Report, Q&A, Logs).

**Whether any git command changed state.** No. None was run.

### 2026-09-30 (Prompt 10: upstream reports drafted for PypeIt and for `pyodine`; nothing sent)

**Task.** Draft the upstream report to Ryan Cooke (phase 2's eleven items
plus phase 3's), with `refframe_audit.py` and `verify_orig_fixes.py` as
reproducers, and show it before anything is sent. The new Q&A answer ((a):
prompt 7 stays the headline) bears on prompt 11, not this one. Findings are
in `## Report`, "Prompt 10"; the sending decisions are in the Q&A.

**What was done.**
- Re-read phase 1's nine items and phase 2's two.
- Checked each against the PypeIt checkout, read-only:
  - `orig-hires-fixes` is on `origin` at `017bece06`, one commit on
    `develop` `f3a1f1d27`, touching `keck_hires.py` only;
  - no PR is open;
  - the `wave.py` sign error is still in `develop` (line 128);
  - upstream `develop` was checked through the GitHub API, not
    `git fetch`.
- Ran both reproducers. Fixed `refframe_audit.py`'s glob so it covers
  exactly the 36 HD 187123 frames.
- Verified the two new spec1d claims directly: zero wavelengths at masked
  pixels, and biased boxcar flux under partially masked apertures.
- Wrote `docs/upstream/pypeit_report_draft.md` and
  `docs/upstream/pyodine_report_draft.md`.

**Headline results.**
- **To PypeIt:** 15 items. Four are new from phase 3: the ivar scale, zero
  wavelengths at masked pixels, the unflagged boxcar bias, and silently
  dropped orders.
- **To `pyodine`:** 6 fork changes and 15 found-and-worked-around items.
- Nothing sent.

**What this taught us about the repository and the data.**

*A reproducer rots when the data around it grows.* `refframe_audit.py`'s
`reduce_1998*` glob was right when written. Since then it has silently
taken in a B-star reduction of other stars, and it had always skipped
December. Globs over a data directory other prompts write into need to name
exactly what they mean.

*An upstream item should be re-checked against upstream before it is
sent.* The sign error's line number had moved (89 → 128). Its status ("still
in `develop`") is what makes it worth sending, so it was checked, not
assumed.

*Claims about data files get verified on the files.* Two new PypeIt items
(13, 14) came from memory of earlier prompts. Both were re-measured before
drafting, and item 14 turned out to be sharper than remembered: the flux
bias follows the fraction of the aperture summed, not only bad columns.

*The findings split by audience.* Phase 3's `pyodine` findings would have
been lost in a PypeIt report.

**Files added / modified.**
- Added: `docs/upstream/pypeit_report_draft.md` and
  `docs/upstream/pyodine_report_draft.md`.
- Modified:
  - `first_hires_exoplanet/refframe_audit.py` (glob);
  - `first_hires_exoplanet/data/refframe_audit.csv` (regenerated: 36 rows,
    was 10);
  - this document (Report, Q&A, Logs).

**Whether any git command changed state.** No. Only read-only `git log`,
`git branch -r --contains`, `git show`, `git diff` and `git status`, plus
`gh api` GET requests against the public repository.

### 2026-09-30 (Prompt 11: phase-3 assessment written; the public document tells what the three phases found)

**Task.** Write the phase-3 assessment, and update `docs/public_HD187123b.md`
with what the three phases actually found, including where we fell short of
1998 and why. The new Q&A answers (after prompt 10): you send the reports;
PypeIt as email plus issues and a PR; `pyodine` as GitHub issues. Nothing in
them changes this prompt. Nothing was opened or sent from here. Findings are
in `## Report`, "Prompt 11".

**What was done.**
- Wrote the assessment from the committed tables.
- Added `fig_three_ways` and `fig_precision_ladder` to
  `first_hires_exoplanet/figs.py`, which read `iodine_velocities.csv` and
  `xcorr_velocities.csv`, and print the scatter numbers the document quotes.
- Regenerated the figures. Restored `fig1`–`fig6` from `HEAD` with
  `git show HEAD:… > file` (read-only git), since they re-render
  byte-differently.
- Replaced the public document's closing "we have not done it yet" paragraph
  and added the two new sections, the `pyodine` reference, and data notes
  for figures 7–8.

**Headline results.** Phase 3 met its own criterion:
- 25.8 m/s per epoch against a 50 m/s target;
- 10.2σ, with the period found blind;
- independent and open.

It is 28× phase 2 and ~9–12× short of the 1998 and modern results. The
public document says so in those numbers.

**What this taught us about the repository and the data.**

*The public document is the place hindsight is most tempting.* The 11.5 m/s
variant would make a better headline. Both documents keep 26 m/s, as the
Q&A decided, and say why.

*Facts written for a general reader need the same checks as the Report.*
Three statements in the first draft of the public section were wrong or
unsupported, and were corrected before finishing:
- the false-alarm odds (1 in 50 billion; the right figure is 1 in 40
  billion);
- the `pyodine` item count;
- a reference title that was not in any project record.

*`figs.py` is environment-sensitive in its bytes, not its content.*
Regenerating it in either conda environment rewrites figures 1–6. Restore
the committed ones, or accept the churn knowingly.

*The prompt doc's own prediction about the LSF was wrong, and the record
should say so.* The assessment does.

**Files added / modified.**
- Added: `docs/figs/fig7_three_ways.png` and `fig8_precision.png`.
- Modified: `first_hires_exoplanet/figs.py`; `docs/public_HD187123b.md`;
  this document (Report, Logs).

**Whether any git command changed state.** No. `git show` (read-only) wrote
the committed figures 1–6 back to the working tree, and `git status` was
used for checking.

### 2026-09-30 (Prompt 12: technical slide deck for phases 1–3, 28 slides, built by a script)

**Task.** Generate a technical PowerPoint deck on what phases 1–3
accomplished, at `docs/slides/technical_summary.pptx`, with figures from
Python scripts, for upload to Google Drive. No new Q&A answers since
prompt 11.

**What was done.**
- Installed `python-pptx` 1.0.2 into `pypeit14` with pip; it was in neither
  environment. It brought `XlsxWriter` and `lxml`.
- Wrote `first_hires_exoplanet/make_slides.py`. It:
  - takes the headline numbers from the committed tables
    (`iodine_velocities.csv`, `iodine_keplerian_residuals.csv`,
    `diagnose_variants.csv`, `xcorr_velocities.csv`) and recomputes K with
    `figs_phase2.fit_circular`;
  - draws one new schematic, `docs/figs/fig_slides_pipeline.png`;
  - assembles a 16:9 deck from the figures the phase scripts already
    generate; the docstring maps each figure to its script.
- Rendered the deck through LibreOffice to PDF and PNG (scratch only) to
  check the layout.

**The deck (28 slides).**
- **Opening:** title; the bottom line; the three-phase schematic.
- **Phase 1 (4):** the source fixes and parameters as a table;
  wavelength-solution quality; order tracing and the flat.
- **Phase 2 (5):** what was built; the atlas; the non-detection; the
  47-against-702 m/s diagnosis.
- **Phase 3 (13):** the fork's six changes; deconvolution; one epoch; LSF
  models; B stars; all epochs; the detection; the catalogue comparison; the
  residuals and the blind period search; the calibration defects fixed; the
  prompt-9 variants (with the Q&A's "reported, not claimed"); the precision
  ladder.
- **Summary (3):** why 26 m/s and not 3; the upstream reports and next steps;
  reproducibility (commits and scripts).
- Every content slide carries a one-line takeaway.

**What this taught us about the repository and the data.**

*Every figure the deck needed already existed and came from a script.* The
discipline of "calculations become scripts" paid off here: the deck is a
layout job over committed products, plus one schematic. It re-runs in
seconds if a number changes.

*Headline numbers in a deck should be computed, not typed.* K, its error, the
residual and the variant rows are read from the tables at build time. Only
numbers that live in a Report and nowhere in a table are typed: phase 1's
fixes, phase 2's diagnostics, the per-prompt findings.

*Check the layout by rendering, not by reading the code.* LibreOffice's
headless conversion plus a contact sheet showed every slide at once. It
caught one syntax error before that, and nothing in the layout needed
changing.

**Files added / modified.**
- Added: `first_hires_exoplanet/make_slides.py`;
  `docs/slides/technical_summary.pptx`;
  `docs/figs/fig_slides_pipeline.png`.
- Modified: this document (Logs). `pypeit14` gained `python-pptx`.

**Whether any git command changed state.** No. None was run.

### 2026-09-30 (Prompt 14: the 1998 paper's result and velocity figures added to the deck as slide 2)

**Task.** Add a slide showing the results of the original 1998 paper, with
its velocity-curve figure. No new Q&A answers.

**What was done.**
- Obtained the publisher's PDF of Butler, Marcy, Vogt & Apps 1998, PASP,
  110, 1389 (DOI 10.1086/316287). It is 5 pages; the paper had not been read
  end to end before (see `context_prompts.md`, 2026-09-19).
- Wrote `first_hires_exoplanet/fetch_butler1998.py`. It downloads the PDF to
  `../first-hires-exoplanet-data/literature/` unless it is already there. It
  then extracts the two page-2 images with poppler's `pdfimages`, as
  `butler1998_fig1_july.png` and `butler1998_fig2_phased.png`.
- Added `Deck.two_figures` and a new slide 2 to `make_slides.py`. The
  slide shows both figures, the paper's Table 1 (P, K, e) and M sin i, the
  data description, and a takeaway that sets our K and rms, computed from
  the tables, against the target. If the figures are missing,
  `make_slides.py` calls the fetch script.
- Rebuilt the deck (31 slides) and rendered it through LibreOffice to check
  the layout. Every text run is still ≥ 20 pt (prompt 13).

**Headline results.** The paper gives two velocity figures, and the slide
uses both:
- **Fig. 1:** the 10 velocities from the 1998 July 15–19 run, against time.
- **Fig. 2:** all 20 velocities, phased.

Its Table 1 gives P = 3.097 ± 0.003 d, K = 72.0 ± 2.0 m/s, e = 0.03 ± 0.03
and T₀ = JD 2,451,010.982 ± 0.01. The text gives more:
- the Keplerian fit has K = 73.0 m/s and rms 7.50 m/s;
- the sinusoid has K = 71.6 m/s and rms 7.62 m/s;
- internal errors are ~6 m/s: ~3 m/s from photons and ~4 m/s systematic;
- R = 87,000, over 3900–6200 Å.

**What this taught us about the repository and the data.**

*The paper answers questions phase 0 left open.*
- The discovery set is **20 velocities**, from 1997 Dec 23 to 1998 Aug.
  Our 31 cell-in epochs run to 1998 December, so they are not the same set.
- The first two velocities (1997 Dec 23 and 1998 Jun 18) differed by
  84 m/s. That difference is what triggered the five-night July run.
- The paper prints no velocity table, only the figures. So a direct
  epoch-by-epoch comparison with the 1998 numbers is still not possible
  (phase 0, Q13).

*IOPscience's bot protection is intermittent, not absolute.*
- One `curl` through the ADS `PUB_PDF` link gateway returned the PDF.
- Minutes later, the same URL and the direct IOP URL both redirected to the
  Radware validation page, from `curl`-like and from Python requests alike.
- `fetch_butler1998.py` therefore falls back to a PDF saved by hand at a
  fixed path. It does not try to get around the check.

*The figures are embedded as 300-dpi 1-bit images with their axes.* So they
can be extracted losslessly, and no page cropping is needed.

*`conda run` does not pass stdin through.* A `python - <<EOF` heredoc under
`conda run -n pypeit14` silently does nothing. Write the script to a file.

*The figures are the ASP's copyright.* They are kept outside the repository
with the PDF, and they are not committed. The `.pptx` embeds them. That is
fine for a technical talk with the citation on the slide, but it matters if
the deck is posted publicly.

*Prompt 13 has no log entry here.* Its change is in the code
(`MIN_PT = 20` and `_pt`), and the deck already had 30 slides before this
prompt, not prompt 12's 28. The prompt-13 session evidently did not log.

**Files added / modified.**
- Added: `first_hires_exoplanet/fetch_butler1998.py`.
- Modified: `first_hires_exoplanet/make_slides.py`;
  `docs/slides/technical_summary.pptx`; `docs/figs/fig_slides_pipeline.png`
  (regenerated, unchanged in content); this document (Logs).
- Outside the repository: `../first-hires-exoplanet-data/literature/`
  (the PDF and the two figure PNGs).

**Whether any git command changed state.** No. None was run.

### 2026-09-30 (Prompt 15: the Ryan Cooke draft already existed; emailed to the user for review, not to Ryan)

**Task.** Say whether the email for Ryan Cooke had been drafted. If it had,
send it to the user for review; if not, draft it. No new Q&A answers.

**What was done.**
- Confirmed that the draft exists: `docs/upstream/pypeit_report_draft.md`,
  written in prompt 10 (1,625 words, 15 items in sections A–F). It was left
  unchanged.
- Stripped the Markdown (bold, backticks, `##` headings) for a plain-text
  email, and sent it through the Gmail connector to the user only
  (jxp@ucsc.edu). The subject is prefixed "[DRAFT for review]". Gmail
  message id `1a0f4d40776d110c`.
- Put a cover note at the top of the email. It says the email has not gone
  to Ryan, and it raises the two points below.

**Headline results.** Nothing went to Ryan. Two things are open before it
does:
- `<repository link>` is still a placeholder. The git remote is
  `https://github.com/pypeit/first-hires-exoplanet`.
- The closing paragraph asks Ryan how he wants issues and PRs. The Q&A
  answer after prompt 10 was (c): the email, then the issues and the PR. The
  close could say so.

**What this taught us about the repository and the data.**
- *The upstream drafts are the source; the email is a rendering of one.*
  Edits belong in `docs/upstream/*.md`, so that the committed text matches
  what is sent.
- *Gmail's `send_message` wants plain text, not Markdown.* A `sed` pass that
  strips bold, backticks and heading marks is enough. The draft's lists and
  numbering read fine as plain text.
- *The prompt-10 Q&A has moved past the draft's ending.* The answer chose
  issues and a PR alongside the email, but the draft still asks Ryan which
  he prefers.

**Files added / modified.**
- Modified: this document (Logs).
- Outside the repository: one email sent, to the user only.

**Whether any git command changed state.** No. `git remote -v` (read-only)
was run.

### 2026-10-01 (Prompt 16: two public slides — the ten prompts, and the project as six lessons)

**Task.** Make two slides on the whole project for the general public. They
should stress (1) that people could do most of it themselves, showing the
~10 key prompts, and (2) that the knowledge could pass to a high-school
teacher who leads students through it.

**What was done.**
- Read the prompt lists and Q&A sections of all five prompt docs, to find
  which prompts the user actually wrote.
- Wrote `first_hires_exoplanet/make_public_slides.py`, which builds
  `docs/slides/public_summary.pptx` (2 slides, 16:9). It reuses `Deck` and
  the 20-pt floor from `make_slides.py`.
- The slide-1 counts are computed from the prompt docs at build time:
  - **steps:** numbered items under each doc's `## Prompts`;
  - **decisions:** `>A.` or `Answer:` lines in the Q&A sections;
  - **days:** the first to the last dated log entry.
- Rendered the slides through LibreOffice and revised once. Two tiles
  wrapped and one lesson title ran to two lines. "Mostly (a) or (b)"
  overstated it, so it became "often just (a)".

**The slides.**
- **Slide 1, "You could do most of this yourself."**
  - Ten prompts, each shortened from the original. `PROMPTS` records where
    each one is.
  - The two workhorses are highlighted: "Reread the file and execute
    prompt #N" and "Answer: (a)".
  - Four tiles: 12 days, 49 steps, 26 decisions, and "all free" (public data
    and software).
  - Takeaway: the human set the goal and made the calls; the AI wrote the
    plans, code and reports.
- **Slide 2, "A teacher could lead a class through it."** Six lessons
  follow the project's arc:
  - the wobble; the archive; light into spectra;
  - lamps fall short (phase 2's non-detection, used as a lesson);
  - a ruler of iodine; find the planet.

  Below them, what the teacher would be handed: the repository, the ten
  prompts, the public write-up, and the judgment calls. Takeaway: the scarce
  skill is knowing what to ask and when to doubt an answer.

**Headline results.** The counts are 49 numbered steps across five prompt
docs, 26 recorded decisions, and 12 days (2026-09-19 to 2026-09-30).

**What this taught us about the repository and the data.**
- *The user wrote few of the detailed prompts.*
  - Each phase's numbered prompts were written by Claude, on request (for
    example, context prompt 7: "generate a prompt doc").
  - The user's own words were the goal, a handful of steering prompts, the
    Q&A answers, and "execute prompt #N".
  - That is the point slide 1 makes, so it is stated as such and not hidden.
- *The prompt docs work as a dataset.* Their regular structure (`## Prompts`
  as a numbered list, `>A.` / `Answer:` in Q&A, `### YYYY-MM-DD` logs) means
  the project's own history can be counted by a script. Keep that format.
- *Phases 1 and 2 recorded no Q&A answers.* All 26 decisions are in
  `context_prompts.md` (14) and phase 3 (12). The phase-1 and phase-2
  decisions were steering prompts, such as phase 1 #6, "Fix the issues you
  have identified."

**Files added / modified.**
- Added: `first_hires_exoplanet/make_public_slides.py`;
  `docs/slides/public_summary.pptx`.
- Modified: this document (Logs).

**Whether any git command changed state.** No. None was run.

### 2026-10-06 (Prompt 17: QR code for the public report, written by a script and verified by decoding)

**Task.** Generate a QR code for a slide, pointing at
`https://github.com/pypeit/first-hires-exoplanet/blob/main/docs/public_HD187123b.md`.
The branch is now merged into `main`.

**What was done.**
- Checked that the URL resolves: HTTP 200.
- Installed `segno` 1.6.6 into `pypeit14` with pip. It is pure Python with
  no dependencies. Neither `segno` nor `qrcode` was installed before.
- Wrote `first_hires_exoplanet/make_qr.py`. It writes
  `docs/figs/qr_public_report.png` (1197 px square) and
  `docs/figs/qr_public_report.svg`. Both are black on white, at
  error-correction level H, with the standard 4-module quiet zone.
- Verified both images by decoding them with OpenCV's `QRCodeDetector`.
  OpenCV was installed into a scratchpad `--target` directory, not into
  `pypeit14`.

**Headline results.**
- The code is QR version 8, 57 modules across with the border.
- It decodes to the exact URL at full size, and still decodes when shrunk
  to 150 px.

**What this taught us about the repository and the data.**
- *Level H costs little here.* The URL is 84 characters. At level H it
  needs version 8, which is still coarse enough to scan from across a room.
  It also tolerates a logo or a partial occlusion, at 30% loss.
- *Keep the white border on the slide.* Without the quiet zone against a
  dark background, many phone cameras fail to read the code.
- *A QR code is only as durable as its URL.* This one points at `main`. If
  `docs/public_HD187123b.md` is renamed, or the repository moves, printed
  codes break. A short redirect, or a GitHub Pages URL, would decouple them.
  `make_qr.py --url` regenerates the code for a new address.
- *Verify a code by decoding it, not by looking at it.* No decoder was on
  the machine (no `zbarimg`, no `cv2`). A throwaway OpenCV install in the
  scratchpad did the check without touching the environment.

**Files added / modified.**
- Added: `first_hires_exoplanet/make_qr.py`;
  `docs/figs/qr_public_report.png`; `docs/figs/qr_public_report.svg`.
- Modified: this document (Logs). `pypeit14` gained `segno`.

**Whether any git command changed state.** No. `git status` and `git log`
(read-only) were run.
