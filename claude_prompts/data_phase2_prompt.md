# Data, Phase 2: relative velocities, and the inputs `pyodine` will need

## Goals

Item 2 of the three-step scope in `context_prompts.md` Q3: **relative radial
velocities by cross-correlation**, from the extracted spectra phase 1 produced.

## The decision this document makes

Prompt 8 asked whether phase 2's cross-correlation is worth doing on its own
terms, or is better framed as the step that produces the inputs and diagnostics
`pyodine` will need.

**It is framed as the latter. Cross-correlation is the instrument here, not the
deliverable.** The reasoning, from evidence already in hand:

- The I2 forest **dominates 5000-6200 A and does not move with the star**, so a
  stellar cross-correlation is pushed onto the ~20 orders blueward of 5000 A,
  which are also the low-S/N end of the reduction (33-96 against a median of
  117-146).
- Context Q3 and the phase-1 assessment independently put cross-correlation at
  **tens to hundreds of m/s** against a **72 m/s** semiamplitude. A result at
  that precision cannot confirm the planet, and presenting it as though it might
  would misrepresent what we did.
- `pyodine` exists, is published, and does the real job. Since prompt 8 was
  written, two things that were assumed to be obstacles have turned out not to
  be (see below): the FTS atlas ships with the code, and the iodine-free stellar
  template exists in the archive. Phase 3 is now a realistic target, so phase 2
  should be built to feed it.

**But phase 2 still measures velocities, and reports honestly what precision it
reached.** Not as a planet detection — as four things we actually need:

1. The **per-epoch initial velocity guess** `pyodine` requires.
2. An **end-to-end consistency check** across the nine-month baseline: the same
   star, reduced independently on many nights, should give the same velocity to
   within the method's precision. If it does not, something is wrong upstream.
3. The **number we will quote** when explaining why a forward model is necessary
   at all. "Cross-correlation gave X m/s scatter; the signal is 72 m/s" is the
   argument for phase 3, and it is worth having X measured rather than asserted.
4. A **filter** that catches bad epochs and bad orders before they reach
   `pyodine`.

Phase 2 succeeds when we have velocities for every discovery-era epoch with
honest uncertainties, a quality-filtered spectrum set, and a written statement of
exactly what `pyodine` still needs.

## What we already know

### Carried forward from phase 1

- **10 epochs reduced**, the 1998 July 15-19 run, each into **37 echelle orders
  (57-93) covering 3806-6262 A**; median S/N 117-146; all wavelength solutions
  under 0.24 px RMS. Reductions live in
  `../first-hires-exoplanet-data/redux/reduce_1998071[5-9]/`, summarised in
  `redux/run_summary.json`.
- **The reduction is externally correct**, not merely self-consistent: seven
  photospheric lines against vacuum rest wavelengths give **-14.7 km/s** against
  HD 187123's **-17 km/s** systemic velocity.
- **PypeIt reports vacuum wavelengths.** Comparing against air rest wavelengths
  produces a spurious +70 km/s. This has already caught us once.
- **Masked pixels arrive as exact zeros**, not flagged gaps.
- **One known-bad order-spectrum**: order 60 of `HI.19980719.49996`, 122 counts
  where the same order holds 51,000-110,000 elsewhere. Adjacent-order S/N
  scatter separates it cleanly (6.2% against 1.8-2.0%) and is the basis for the
  quality filter phase 2 needs.
- **Environment**: `pypeit14`, with PypeIt installed `-e` from
  `Projects/PypeIt/PypeIt`, currently on branch **`orig-hires-fixes`** (six
  fixes on top of `develop` commit `f3a1f1d27`, PypeIt `2.0.2.dev1216`). Quote
  the commit for any reduction whose results we publish.

### New, found while writing this document

Three findings that change what phase 2 should do. Each was checked, not assumed.

**1. The FTS iodine atlas is not a gate any more — it ships with `pyodine`.**
The repository has an `iodine_atlas/` directory containing
`Fischer_Cell_May2022_downsampled3.h5` (54.5 MB) and
`song_iodine_cell_01_65C.h5` (38.6 MB), MIT-licensed with the rest. The first is
almost certainly **Debra Fischer's** cell — the Lick/Keck lineage, the same
tradition as the HIRES cell. This retires the concern recorded in context Q9 and
repeated in the phase-1 assessment that a public FTS atlas could not be found.

It does **not** retire the question. An FTS atlas is a scan of *one physical
cell*. The I2 line *positions* are molecular constants and transfer to any cell;
the line *depths* depend on column density and cell temperature, which do not.
Whether the Fischer atlas is an adequate description of the HIRES cell is a real
question with a real answer, and phase 2 should establish it rather than hope.

**2. The iodine-free stellar template exists in the archive.** The Butler method
needs a spectrum of the star taken with the cell **out**. The KOA survey already
in `first_hires_exoplanet/data/koa_hd187123_science.csv` records an `iodin` /
`iodout` column, and of 37 science frames **five are cell-out**:

| Date | Exposure | KOAID |
|---|---|---|
| 1997-12-24 | 500 s | `HI.19971224.16259` |
| 1998-08-12 | 60 s | `HI.19980812.29316` |
| 1998-08-26 | 500 s | `HI.19980826.34749` |
| 1998-08-26 | 500 s | `HI.19980826.35345` |
| 1998-08-26 | 500 s | `HI.19980826.35939` |

The **three consecutive 500 s cell-out exposures on 1998-08-26** are exactly the
deliberate template observation the technique calls for; the 1998-08-12 pair
(60 s cell-in at `.29160`, 60 s cell-out at `.29316`) is a ready-made diagnostic
pair for isolating the cell's transmission empirically.

**None of them is in the July run we reduced.** So phase 2 must extend the
reduction to at least 1998-08-26, and doing so is on the critical path for phase
3, not an optional extra.

**3. `pyodine`'s instrument adapter is small and its contract is legible.**
Adapting an instrument means adding a `utilities_<instrument>/` directory
alongside `utilities_lick/`, `utilities_song/`, `utilities_mtkent/` and
`utilities_waltz/`. The Lick one is five files: `__init__.py`, `conf.py`,
`load_pyodine.py`, `pyodine_parameters.py`, `timeseries_parameters.py`. In
`load_pyodine.py`, `ObservationWrapper(components.Observation)` takes a filename
and exposes `flux` shaped `(nord, npix)`, a weight array, `nord`, `npix`,
`instrument`, `star`, `exp_time`, and — importantly — `bary_date` and
`bary_vel_corr` read from header cards.

### The trap to settle first

**`pyodine` expects to apply the barycentric correction itself; PypeIt has
already applied a heliocentric one.** Phase 1 recorded `VEL_CORR` showing a
+4.2 km/s heliocentric correction in the reduced products. Heliocentric and
barycentric reference frames differ by up to roughly **15 m/s** — about 20% of
the 72 m/s semiamplitude, and far above anything phase 3 aims at. Handing
`pyodine` heliocentric-corrected wavelengths while also giving it a
`bary_vel_corr` would double-correct.

The likely resolution is to work in the **observed** frame (PypeIt's `refframe`
parameter) and let `pyodine` apply a proper barycentric correction computed with
`astropy`, but this must be established from the code and the FITS headers, not
assumed. It is prompt 1.

### Open questions worth carrying

- The nine upstream items for Ryan Cooke from phase 1 are still unreported.
  Phase 2 is a reasonable moment to send them.
- Gain and read noise: headers give `CCDGN01 = 4.8`, `CCDRN01 = 6.0` against
  PypeIt's hard-coded 1.9 e-/ADU and 2.8 e-. Phase 1 deliberately did not change
  these. They affect the *weights* a cross-correlation and a forward model use,
  so phase 2 has more reason to care than phase 1 did.
- The KOA survey lists 32 cell-in science frames against the 30 discovery-era
  velocities in the modern catalogue. The two short 60 s frames (1998-07-14,
  1998-08-12) are the likely explanation, but it is worth confirming rather than
  assuming, since it tells us whether we can expect a one-to-one comparison.

## Code

See the guidelines in the prompt docs in
`Projects/PypeIt/PypeIt-development-suite/claude_prompts` for how to code in
Python.

Use the `pypeit14` environment. `pyodine` is a **reference implementation to
adapt, not a dependency to trust** — the repository is sparse (few commits, no
issue traffic) even though the paper is solid. Vendor or fork it deliberately
rather than `pip install` from a moving target, and record the commit.

Existing modules to build on, not duplicate:
`koa_survey.py`, `koa_download.py`, `run_setup.py`, `reduce_run.py`,
`figs_phase1.py`, `check_env.py`, `verify_orig_fixes.py`.

Keep raw frames and reductions in `../first-hires-exoplanet-data/`, outside the
repository. Small derived products — velocity tables, quality tables — belong in
`first_hires_exoplanet/data/` and should be committed.

Record the exact PypeIt commit for any reduction whose results we quote, and the
exact `pyodine` commit for anything derived from it.

## Prompts

1. Read this file. Settle the reference-frame question before anything else.
   Determine from the FITS headers and the PypeIt source exactly what correction
   was applied to the phase-1 reductions, what `refframe` was set to, and what
   `pyodine`'s `ObservationWrapper` expects in `bary_date` and `bary_vel_corr`.
   Recommend how we should carry velocities through phase 2 and 3 without
   double-correcting. Write it into the Report section. **Change nothing yet.**
   Use Opus 5. Log your work.

2. Read this file. Write the quality filter as a script on disk: read every
   order-spectrum from the five reduced nights, compute the adjacent-order S/N
   scatter diagnostic that isolated order 60 of `HI.19980719.49996`, and emit a
   per-epoch, per-order quality table to `first_hires_exoplanet/data/`. Report
   how many orders it flags and whether any flag is a false positive. Use Opus 5.
   Log your work.

3. Read this file. Extend the reduction to the rest of the discovery era: the
   1998-08-26 night (which carries the three cell-out template exposures **and**
   three cell-in science frames), then the remaining cell-in nights — 1998-08-17,
   08-18, 08-25, and 1998-09-12 through 09-18, plus 1997-12-23/24 and
   1998-06-18. Reuse `reduce_run.py`. Report which nights reduce cleanly, which
   need intervention, and whether the recipe settled in phase 1 holds across a
   nine-month baseline and a change of decker (1998-09-17 is B2, not B1).
   Use Opus 5. Log your work.

4. Read this file. Build the stellar template: co-add the three cell-out
   exposures of 1998-08-26 into a single high-S/N, iodine-free spectrum of
   HD 187123, order by order. Check it against the 1997-12-24 cell-out frame for
   consistency, and against a G2V expectation. Report the S/N achieved per order.
   Use Opus 5. Log your work.

5. Read this file. Measure relative velocities by cross-correlation against that
   template, restricted to the orders blueward of 5000 A where the stellar lines
   are not swamped by the I2 forest. Produce a velocity per epoch with an honest
   uncertainty, weighted per order, and write the table to
   `first_hires_exoplanet/data/`. Report the per-order scatter and the
   epoch-to-epoch scatter separately — they answer different questions.
   Use Opus 5. Log your work.

6. Read this file. Assess what prompt 5 produced. Compare the velocities against
   the published 3.097-day Keplerian and against the modern catalogue values in
   `hd187123_hires_rv.tsv`. State plainly what the comparison does and does not
   demonstrate, and what the achieved precision implies for phase 3. Generate
   figures with a script on disk, as in phase 1. Use Opus 5. Log your work.

7. Read this file. Settle the atlas question. Obtain
   `iodine_atlas/Fischer_Cell_May2022_downsampled3.h5` from `pyodine`, establish
   its provenance (which cell, what temperature, what resolution, what
   wavelength coverage), and compare it against what is known of the HIRES
   iodine cell. Report whether it can serve as the HIRES cell's atlas, whether
   the 1998-08-12 cell-in/cell-out pair can test that empirically, and if not
   what we would need instead. Use Opus 5. Log your work.

8. Read this file. Write the `pyodine` input adapter: a module that presents a
   PypeIt `spec1d` as the `(nord, npix)` flux, weights, wavelengths, `bary_date`
   and `bary_vel_corr` that `ObservationWrapper` expects, following the
   resolution from prompt 1 on reference frames and the quality flags from
   prompt 2. Model it on `utilities_lick/load_pyodine.py`. Do not yet attempt a
   forward model. Use Opus 5. Log your work.

9. Read this file. Write the phase-2 assessment: what the cross-correlation
   achieved and what it cost, whether the extended reduction held up, the state
   of the template and the atlas, and a specific list of what phase 3 needs —
   including anything in `pyodine` that will have to change for HIRES. Also
   decide whether the nine upstream items from phase 1 should go to Ryan Cooke
   now. Use Opus 5. Log your work.

## Q&A

## Report

## Logs
