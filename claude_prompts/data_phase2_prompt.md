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

### Prompt 1: the reference frame

Settled from the products, the PypeIt source and the `pyodine` source. Nothing
was changed. Everything below is reproduced by
`first_hires_exoplanet/refframe_audit.py`, which writes
`first_hires_exoplanet/data/refframe_audit.csv`:

```
conda run -n pypeit14 python -m first_hires_exoplanet.refframe_audit
```

Provenance: PypeIt `2.0.2.dev1217+g017bece06` (branch `orig-hires-fixes`);
`pyodine` commit `4488b0914fe5b272b787982647691045bff2604a` (`main`, 2022-11-04,
"Adding pyodine to GitHub" — the repository's only commit).

#### 1. What PypeIt actually did

`refframe` was never set in the pypeit file, so the phase-1 reductions ran on
the PypeIt **default**, `[calibrations][wavelengths] refframe = heliocentric`
(`pypeit/par/pypeitpar.py:3475`). The correction is applied at the end of
reduction in `pypeit_steps.refframe_correct`, which calls
`wave.geomotion_correct` and then `SpecObj.apply_helio` — the latter multiplies
`OPT_WAVE`/`BOX_WAVE` in place by `VEL_CORR` and records `VEL_TYPE`.

Checked on all ten epochs:

- `VEL_TYPE = 'heliocentric'` in every order of every epoch.
- `VEL_CORR` is **identical across all 37 orders** (spread exactly 0), and
  reproduces `wave.geomotion_correct` to 4e-8 m/s.
- The applied correction runs **+3.434 to +4.888 km/s** over the five nights —
  a 1.45 km/s swing, twenty thousand times the 72 m/s semiamplitude.
- Dividing the stored wavelengths by `VEL_CORR` recovers the observed frame
  **exactly** (round-trip error 0.000 m/s). The correction is a single scalar
  multiply and is losslessly invertible. No re-reduction is required to undo it.

#### 2. PypeIt's `heliocentric` has a sign error on the solar term

This was not expected and is the most consequential finding.

`pypeit/core/wave.py:89` builds the observer's barycentric velocity as
`ev + ov` (Earth's barycentric velocity plus the observatory's GCRS velocity),
then for the heliocentric frame does:

```python
if frame == "heliocentric":
    sp, sv = solar_system.get_body_barycentric_posvel('sun', time)
    velocity += sv
```

The observer's velocity *relative to the Sun* is `ev + ov - sv`. Recomputing
both ways and comparing against `astropy`'s `radial_velocity_correction`:

| projected onto the line of sight | minus astropy heliocentric |
|---|---|
| `ev + ov + sv` — what PypeIt does | **+13.09 m/s** |
| `ev + ov - sv` — the definition | **−0.00 m/s** |
| `ev + ov` — no solar term at all | +6.55 m/s |

The `-` form lands on `astropy` exactly. It is a sign error, not a convention
difference, and the error is twice the solar term: `n̂ · v_sun` is +6.5 m/s for
HD 187123, so PypeIt's heliocentric velocities are **+13.1 m/s** off.

It does not cancel. Scanned across the discovery era (1997-12 to 1998-09) the
offset moves over **+9.06 to +14.22 m/s, 5.16 m/s peak to peak** — 7% of the
semiamplitude, injected as a slow drift across exactly the baseline we care
about. **This belongs on the list for Ryan Cooke as item 10**, and it is an
upstream bug affecting every PypeIt user who takes the default and wants better
than ~15 m/s.

#### 3. PypeIt's `barycentric` is sound, to a near-constant −4.6 m/s

`refframe = barycentric` has no sign error. It differs from `astropy` by
**−4.60 to −4.57 m/s** across the ten epochs, and only **0.22 m/s peak to peak**
across the whole nine-month baseline, while the barycentric velocity itself
swings −13.4 to +16.9 km/s. It is an offset, not a drift.

The cause is that PypeIt's calculation is purely kinematic. `astropy`
additionally carries the solar gravitational redshift at the Earth
(**+2.91 m/s**) and the observer's special-relativistic time dilation
(**+1.43 m/s**), totalling **+4.34 m/s** against the measured −4.58 m/s. That
accounts for the offset to within 0.24 m/s and explains why it barely drifts.

A constant offset is harmless for relative velocities. It is still the wrong
number to hand a code that will add it to a fitted Doppler shift.

#### 4. PypeIt corrected at the *start* of the exposure

`construct_obstime` (`pypeit/metadata.py:471`) returns `Time(meta['mjd'])`, and
for `keck_hires` `meta['mjd']` is the raw header `MJD` card. Nothing in the raw
header says whether that is the start, middle or end of the exposure — there is
one time card (`UT = 10:34:14.18`), and `DATE` duplicates it.

The header's hour angle settles it. `HA + RA = LST`, so the LST implied by the
recorded `HA` can be compared against the LST computed from the recorded time:

| epoch tested | LST(HA) − LST(epoch) |
|---|---|
| header MJD | **+1.4 ± 0.5 s** |
| MJD + ½ exptime | −173.8 ± 38.2 s |
| MJD + exptime | −348.8 ± 76.7 s |

The header MJD is the **exposure start**, to within 1.4 s. Mid-exposure is
`MJD + EXPTIME/2`; exposure times run 215–427 s.

Using the start instead of the midpoint costs **+3.0 to +4.4 m/s**, varying
epoch to epoch because it is driven by the diurnal term.

#### 5. Heliocentric versus barycentric, in this dataset

`astropy` heliocentric − barycentric is **−11.16 to −11.10 m/s** over the July
run, and −11.70 to −9.24 m/s (2.47 m/s peak to peak) over the full discovery
era. The ~15 m/s figure quoted above in "The trap to settle first" is the right
order; the actual number for this target is ~11 m/s with a 2.5 m/s drift.

Total damage from taking PypeIt's heliocentric value as if it were a
barycentric one at mid-exposure: ~11 m/s frame + 13 m/s sign error + 4 m/s
epoch ≈ **28 m/s**, of which ~8 m/s drifts across the baseline. Against a 72 m/s
semiamplitude.

#### 6. What `pyodine` expects

Read from the source, not inferred.

**It wants the observed (topocentric) frame.** `pyodine/chunks.py:172` sets the
initial chunk velocity to `init_dv = temp.bary_vel_corr - obs.bary_vel_corr`,
`init_z = init_dv / c`, and shifts the template by
`wave_shifted = temp.w0 * (1 + init_z)` to line it up with the observation.
That identity only holds if both spectra carry observed-frame wavelengths:
λ_obs = λ_bary·(1 − BVC/c), so λ_sci/λ_temp = 1 + (BVC_temp − BVC_obs)/c, which
is exactly `1 + init_z`. Feeding it frame-corrected wavelengths would cancel the
very shift the model is built to track.

**The correction is applied at the very end**, in
`pyodine/timeseries/base.py:102`: `timeseries['rv_bc'] = timeseries['rv'] + bvcs`,
with `bvcs` computed fresh by `bvc_wrapper` → `barycorrpy.get_BC_vel`.

**`ObservationWrapper` does not compute either quantity** — it reads both from
header cards (`utilities_lick/load_pyodine.py:114-115`: `LICKJD`, `LICKBVC`,
with `# TODO: Re-calculate BVC` alongside). Whatever we put in those slots is
taken on trust.

Two docstrings in `pyodine/components.py` are **stale, and both are traps**:

- Line 315 says `bary_date` is a "Barycentric Reduced Julian Date (BJD −
  2400000.0)". It is not. `bvc_wrapper` passes it straight to
  `barycorrpy.get_BC_vel(JDUTC=...)` and to `utc_tdb.JDUTC_to_BJDTDB`, both of
  which want a **full JD in UTC**. The Lick helper `get_barytime` confirms it,
  returning `Time(...).jd`. Passing a reduced JD would put the observation in
  1858.
- Line 316 says `bary_vel_corr` is in km/s. It is in **m/s**: `chunks.py`
  imports `astropy.constants.c` and divides by `c.value` (2.998e8 m/s), and
  `barycorrpy.get_BC_vel` returns m/s.

#### 7. Recommendation

1. **Work in the observed frame throughout.** Set
   `[calibrations][wavelengths] refframe = observed` in the pypeit file for
   every reduction from prompt 3 onward, and re-reduce the ten July epochs
   under the same setting in that same pass, so the whole discovery-era set is
   homogeneous. Mixing settings across nine months of data is precisely the
   trap this document was written to avoid.

2. **Until then the existing ten are exactly recoverable** as
   `λ_obs = λ_stored / VEL_CORR`. Round-trip error is 0.000 m/s and `VEL_CORR`
   is one scalar per epoch. Nothing needs re-reducing to start phase 2.

3. **The prompt-8 adapter divides by `VEL_CORR` unconditionally**, reading it
   from the spec1d and treating an absent value as 1.0. That is correct under
   either `refframe` setting and is the single structural guard against
   double-correction — it cannot silently do the wrong thing if someone later
   reduces a night with a different default.

4. **Never reuse PypeIt's number as `bary_vel_corr`.** Compute it with
   `barycorrpy` — the same library `pyodine`'s own timeseries stage uses, so the
   initial guess and the final correction come from one source — at the
   flux-weighted midpoint, with HD 187123's proper motion and `rv0 = −17 km/s`.
   This avoids all three errors at once: the +13.1 m/s sign error (drifting
   5.2 m/s over the baseline), the −4.6 m/s missing relativistic terms, and the
   3–4.4 m/s start-of-exposure epoch.

5. **`bary_date` = full JD in UTC at the flux-weighted midpoint**, i.e.
   `MJD_header + EXPTIME/2/86400 + 2400000.5`. Not BJD, not reduced. And
   `bary_vel_corr` in **m/s**.

6. **Record the irreducible floor.** 1998 HIRES had no exposure meter, so the
   geometric midpoint is the best midpoint available. Under variable
   transparency the true flux-weighted midpoint can sit tens of seconds away;
   at the diurnal rate of ~0.03 m/s/s that is ~1 m/s of unremovable systematic.
   Phase 3 should carry it as a known floor rather than discover it later.

7. **Phase 2 is not threatened by any of this** — every effect here is under
   30 m/s and cross-correlation will land at tens to hundreds. But phase 2
   should adopt the observed frame and `barycorrpy` anyway, for two reasons:
   the applied heliocentric correction swings 1.45 km/s across the five July
   nights, so the epoch-to-epoch consistency check in goal 2 is meaningless
   unless the frame is handled consistently; and velocities produced this way
   are directly comparable against phase 3's rather than needing to be
   reconciled with them afterwards.


### Prompt 2: the quality filter

Written as `first_hires_exoplanet/quality_filter.py`, which emits
`first_hires_exoplanet/data/order_quality.csv` — one row per (epoch, order),
370 rows:

```
conda run -n pypeit14 python -m first_hires_exoplanet.quality_filter
```

Provenance: PypeIt `2.0.2.dev1217+g017bece06` (branch `orig-hires-fixes`).

#### Headline

**68 of 370 order-spectra are rejected (18.4%), not one.** Phase 1's "one
known-bad order-spectrum" is a substantial understatement, and the reason it
was missed is instructive: order 60 of `HI.19980719.49996` was found by chasing
an anomalous S/N *minimum*, and the other failures do not depress the S/N
minimum at all.

**All of the damage is in the iodine region.** Orders 78–89 (3977–4603 Å) are
clean in all ten epochs; every rejection falls in orders 57–77, 4597–6262 Å.

| window | usable | rejected |
|---|---|---|
| blue, < 5000 Å — phase 2 cross-correlation | 204 / 220 | 7% |
| iodine, > 5000 Å — phase 3 forward model | 98 / 150 | **35%** |

Phase 2 is barely affected. Phase 3 loses a third of the region it cannot do
without. That inversion was not visible before this ran.

#### Two failure modes, and why the tests must be sequential

Phase 1 framed this as one statistic. It is two, and they are independent.

**Mode 1 — PypeIt rejected the pixels.** `OPT_MASK` collapses to near-zero
while the flux and inverse variance are healthy. `HI.19980718.49487` order 59
has `OPT_MASK` true for **10 of 2048** pixels, yet 2040 pixels have positive
inverse variance and a median of 18,208 counts. The driver is `OPT_FRAC_USE`,
the fraction of the object profile landing on usable slit pixels: 0.428 there
against 0.997 in a clean frame. Across the dataset the masked orders have
median FWHM 2.28 px and FRAC_USE 0.850; the unmasked ones 1.90 px and 0.988.
Poor seeing widens the profile, it spills past the aperture, and PypeIt throws
the order away. **57 order-spectra**, in six of the ten frames.

Phase 1 wrote `redux/checks/frac_use_stats.py` and never drew a conclusion from
it. The conclusion is that this is the *dominant* quality problem in the
dataset, four times more common than the extraction faults.

**Mode 2 — the extraction silently lost the object.** Order 60 of
`HI.19980719.49996` has `OPT_MASK` true for 2026 of 2048 pixels and
`OPT_FRAC_USE` exactly 1.000. Nothing internal to that order says it is wrong.
Only its neighbours do. **11 order-spectra** flagged, of which 8 confirmed.

Running one test explains why the order matters: computing the smoothness
statistic on all valid pixels — including the ones PypeIt rejected — measures
how much of each order the mask discarded, not whether the extraction worked.
Done that way it flags 78 order-spectra and buries the real faults. The tests
are therefore applied in sequence: reject the mask-collapsed orders first, then
test the smoothness of what survives, using only surviving orders as
neighbours.

#### The adjacent-order diagnostic, and the false positive it would have caused

For each frame, a leave-one-out local linear fit to the six neighbouring orders
(±3) predicts each order's log S/N; the residual measures the departure from
its own frame's blaze trend. The order is excluded from its own fit, so a bad
order cannot drag its prediction towards itself and hide.

**That residual is not noise, and taking it at face value produces a textbook
false positive.** Order 89 sits at +0.062 dex in *every one of the ten epochs* —
the free spectral range steps there, and a local linear fit cannot follow it.
Against the raw scatter that is +60σ, so a naive threshold flags order 89 ten
times out of ten. Orders 88, 87 and 92 do the same at +0.042, +0.023, +0.021
dex. These are properties of the echelle format, not faults, so the per-order
median residual across epochs is subtracted before anything is flagged.

**The threshold is not a sigma cut, because a sigma cut is meaningless here.**
After de-trending, the residuals are sharply bimodal: 95% of order-spectra lie
below 0.005 dex, the worst healthy one reaches 0.015, and then the distribution
jumps — 0.035, 0.060, 0.070, 0.078, 0.089, 0.096, 0.113, 0.180, 0.198, 0.342,
0.741. The robust σ is 0.0010 dex, which would put the worst order at 700σ and
every order above the noise floor "significant". **0.02 dex** — a 5% departure
from the local blaze trend — sits in the middle of the gap and is what the
filter uses.

#### The eleven anomalies, adjudicated

The verification test shares no machinery with the flagging test: it never
looks at a neighbouring order, only at the same order in the other nine frames,
after dividing out each frame's overall S/N level measured on a set of
reference orders common to all epochs.

| epoch | order | S/N | residual | cross-epoch | verdict |
|---|---|---|---|---|---|
| `HI.19980719.49996` | 60 | 10.7 | −0.741 | 0.062 | confirmed |
| `HI.19980719.49996` | 59 | 67.6 | +0.342 | 0.387 | confirmed |
| `HI.19980718.49487` | 93 | 27.3 | −0.198 | 0.612 | confirmed |
| `HI.19980719.49996` | 61 | 68.6 | +0.180 | 0.396 | confirmed |
| `HI.19980718.49487` | 92 | 52.0 | +0.112 | 0.972 | **not confirmed** |
| `HI.19980719.49996` | 63 | 71.5 | +0.096 | 0.426 | confirmed |
| `HI.19980719.49996` | 65 | 76.7 | −0.089 | 0.467 | confirmed |
| `HI.19980719.49996` | 62 | 61.8 | +0.078 | 0.362 | confirmed |
| `HI.19980719.49996` | 58 | 42.3 | +0.070 | 0.241 | confirmed |
| `HI.19980718.49487` | 91 | 57.4 | +0.060 | 0.978 | **not confirmed** |
| `HI.19980718.49487` | 90 | 65.3 | +0.035 | 0.986 | **not confirmed** |

**Yes, there are false positives: three of eleven, and they have a single
explanation.** All three are orders 90, 91 and 92 of `HI.19980718.49487`, and
order 93 of that same frame is a confirmed fault. They sit 1, 2 and 3 orders
away — inside the ±3 neighbourhood. A broken order is a neighbour of the orders
either side of it and drags their predictions with it. The script checks this
explicitly and reports that every unconfirmed flag lies within the
neighbourhood of a confirmed one.

Iterating — barring the first pass's outliers from the second pass's
neighbourhoods — was tried and is **wrong here**, and the failure is worth
recording. In `HI.19980719.49996` the entire red block is damaged, so the
second pass strips away every neighbour order 60 has, its residual becomes
undefined, and *the one fault phase 1 found by hand stops being flagged at
all*. The filter therefore makes a single pass, over-flags in the safe
direction, and lets the independent cross-epoch test adjudicate. Three
conservative false positives on three good orders of an already-damaged frame
is a much better trade than losing the worst order in the dataset.

#### What the confirmed faults actually are

`HI.19980719.49996` is not "one bad order". Its red block reads 42, 68, **11**,
69, 62, 72, —, 77 across orders 58–65 where the other frames run 150–170 and
vary by 2% between neighbours. Seven of its orders are confirmed anomalous and
nine more are mask-collapsed: **16 of 37 rejected**. It is the last frame of
the night, 13:53 UT at airmass 1.52. Atmospheric extinction is smooth in
wavelength and cannot produce order-to-order jaggedness of this size, so this
is the extraction failing, not the sky — consistent with phase 1's guess that
the trace lost the star, but across a whole block rather than one order.

`HI.19980718.49487` order 93 — S/N 27.3 against 40–52 for that order elsewhere
— is a **second, previously unknown** extraction fault, at the far blue end,
inside phase 2's cross-correlation window.

#### Per-epoch summary

| epoch | adjacent-order scatter | rejected |
|---|---|---|
| `HI.19980715.38054` | 1.81% | 0 |
| `HI.19980716.52985` | 2.28% | 4 |
| `HI.19980717.32194` | 2.54% | 3 |
| `HI.19980717.48747` | 5.14% | 16 |
| `HI.19980718.29255` | 1.76% | 0 |
| `HI.19980718.38646` | 2.13% | 7 |
| `HI.19980718.49487` | 5.60% | 22 |
| `HI.19980719.24897` | 1.49% | 0 |
| `HI.19980719.35708` | 1.38% | 0 |
| `HI.19980719.49996` | 7.37% | 16 |

Four epochs are perfect. The frame-level scatter reproduces phase 1's
diagnostic closely — 7.37% for `HI.19980719.49996` against 1.49% and 1.38% for
the other two frames of that night, where phase 1 quoted 6.2% against 1.8–2.0%.
Phase 1's calculation was never written to disk, so the small difference cannot
be traced; the separation is the same and the conclusion is identical.

Note that `HI.19980717.48747` and `HI.19980718.49487` sit at 5.1% and 5.6% —
between the clean frames and the known-bad one. Phase 1 saw only the 6.2%
outlier and reported the run as having a single blemish. The frame-level
statistic did carry the warning; it was not read.

#### What phase 2 should do with this

1. **`order_quality.csv` is the filter.** Booleans are written as 0/1 rather
   than `True`/`False`, because `ascii.csv` round-trips Python booleans as
   strings that evaluate truthy on read — the prompt-8 adapter can consume the
   `use` column directly.
2. **Prompt 5's cross-correlation loses almost nothing**: 204 of 220 blue
   order-spectra survive, and orders 78–89 are clean in every epoch.
3. **Prompt 3 should expect this to recur.** Six of ten frames have
   mask-collapsed orders and the driver is seeing, not the recipe. Over a
   nine-month baseline and a change of decker it will be worse, not better.
   Worth testing during prompt 3 whether relaxing PypeIt's `FRAC_USE`
   threshold, or widening the extraction aperture, recovers the red orders —
   the flux is demonstrably there.
4. **Phase 3's real constraint is now visible.** 35% of the iodine region is
   unusable as reduced, and three epochs retain 2, 2 and 3 usable orders of 15.
   A forward model has nothing to work with on those nights. This should go
   into the prompt-9 assessment as a phase-3 requirement, and it strengthens
   the case for reducing the additional nights in prompt 3.


## Logs

### 2026-09-22 (Prompt 1: settled the reference frame — and found a sign error in PypeIt's heliocentric correction)

**Task.** Prompt 1 of this document: determine from the FITS headers and the
PypeIt source what correction was applied to the phase-1 reductions, what
`refframe` was set to, and what `pyodine`'s `ObservationWrapper` expects, then
recommend how to carry velocities through phases 2 and 3 without
double-correcting. Nothing was changed. Full findings are in the Report section
above; this entry records what was done and what it taught us about the
repository.

**What was done.**

- Traced `refframe` through the PypeIt source: default set at
  `pypeit/par/pypeitpar.py:3475`, applied in `pypeit_steps.refframe_correct`
  (:984), which calls `wave.geomotion_correct` and `SpecObj.apply_helio`
  (`pypeit/specobj.py:647`) — the latter multiplies the wavelength arrays in
  place and writes `VEL_TYPE`/`VEL_CORR`.
- Wrote `first_hires_exoplanet/refframe_audit.py`, which reads all ten reduced
  spec1d files and compares PypeIt's heliocentric and barycentric options
  against `astropy`'s `radial_velocity_correction` at the same site and epoch,
  tests the solar term's sign directly, scans the difference across the
  nine-month discovery baseline, accounts for the residual offset with the
  relativistic terms, and determines from the header hour angle whether the
  header MJD is the start or the end of the exposure. It writes
  `first_hires_exoplanet/data/refframe_audit.csv`.
- Read `pyodine` at commit `4488b0914fe5b272b787982647691045bff2604a` (GitHub
  API, no clone): `utilities_lick/load_pyodine.py`, `components.py`,
  `chunks.py`, `timeseries/base.py`, `timeseries/bary_vel_corr.py`.

**Headline results.** The phase-1 reductions ran on PypeIt's *default*
`refframe = heliocentric` (it was never set in the pypeit file); the applied
correction is +3.434 to +4.888 km/s and is exactly invertible. PypeIt's
heliocentric calculation has a **sign error on the solar term** — `wave.py:89`
does `velocity += sv` where the definition wants `-= sv` — worth **+13.1 m/s**
for this target and drifting **5.2 m/s** across the discovery baseline.
PypeIt's *barycentric* option is sound to a near-constant −4.6 m/s, explained
by the relativistic terms it omits. The header MJD is the exposure **start**
(proved to 1.4 s from `HA + RA = LST`), costing another 3–4.4 m/s.
`pyodine` wants **observed-frame** wavelengths, a **full JD in UTC** at the
flux-weighted midpoint, and `bary_vel_corr` in **m/s**.

**What this taught us about the repository and the data.**

- **PypeIt defaults are load-bearing and invisible.** `refframe = heliocentric`
  never appeared in our pypeit file, so nothing in the reduction inputs records
  that a 4 km/s shift was applied. The only trace is `VEL_TYPE`/`VEL_CORR` in
  the spec1d extensions. Any parameter we do not set is still a choice, and for
  an RV project the defaults must be audited, not inherited.
- **`VEL_CORR` is a clean escape hatch.** One scalar per epoch, identical across
  all 37 orders, round-trip exact to 0.000 m/s. The phase-1 products are not
  spoiled by the default — we can undo it losslessly and start phase 2 today.
  This is worth remembering as a pattern: PypeIt's frame correction is applied
  last and stored, not baked into the wavelength solution.
- **The raw HIRES headers carry exactly one time card**, and it does not say
  what it is. `UT`, `DATE` and `MJD` all give the same instant and none is
  labelled start or end. The hour angle is the only independent timing
  information in the file, and it settled the question. Worth remembering for
  the nights still to be reduced: `HA + RA = LST` is a general-purpose check on
  any epoch we are handed.
- **1998 HIRES has no exposure meter**, so the flux-weighted midpoint is not
  recoverable. The geometric midpoint is the best available and leaves ~1 m/s
  of irreducible systematic. Phase 3 should carry it as a known floor.
- **`pyodine`'s docstrings cannot be trusted where its code can.** Two comments
  in `components.py` (lines 315-316) contradict the code that consumes those
  attributes — the "Reduced Julian Date" would put observations in 1858, and
  the "km/s" is m/s. This is consistent with the document's judgement that
  `pyodine` is a reference implementation to adapt rather than a dependency to
  trust, and it raises the priority of the prompt-8 adapter being written
  against the code paths rather than the documentation.
- **A tenth upstream item for Ryan Cooke.** The heliocentric sign error is not
  specific to HIRES or to our branch — it is in `pypeit/core/wave.py` on
  `develop` and affects any PypeIt user who takes the default and wants better
  than ~15 m/s. It should go with the nine items phase 1 accumulated.

**Files added.**

- `first_hires_exoplanet/refframe_audit.py`
- `first_hires_exoplanet/data/refframe_audit.csv`

**Nothing was changed** in the reductions, the pypeit files, or PypeIt itself,
as the prompt directed.

### 2026-09-22 (Prompt 2: the quality filter — 68 bad order-spectra, not one, and all of them in the iodine region)

**Task.** Prompt 2 of this document: write the quality filter as a script on
disk, compute the adjacent-order S/N scatter diagnostic that isolated order 60
of `HI.19980719.49996` over every order of the five reduced nights, emit a
per-epoch per-order quality table, and report how many orders it flags and
whether any flag is a false positive. Full findings are in the Report section
above.

**What was done.**

- Wrote `first_hires_exoplanet/quality_filter.py`, which reads all 370
  order-spectra, applies three tests in sequence (mask survival, absolute S/N
  floor, adjacent-order smoothness), de-trends the smoothness residual against
  the echelle format, adjudicates every flag with an independent cross-epoch
  test, and writes `first_hires_exoplanet/data/order_quality.csv`.
- Reproduced phase 1's frame-level statistic for continuity: 7.37% for
  `HI.19980719.49996` against 1.49%/1.38% for the other two frames of that
  night (phase 1 quoted 6.2% against 1.8–2.0%).

**Headline results.** 68 of 370 order-spectra rejected, 18.4%. Every rejection
falls in orders 57–77 (4597–6262 Å); orders 78–89 are clean in all ten epochs.
The blue window phase 2 cross-correlates over keeps 204 of 220; the iodine
region phase 3 needs keeps only 98 of 150. Eleven orders flagged by the
adjacent-order statistic, **eight confirmed and three false positives**, all
three of them neighbours of a genuinely bad order in the same frame.

**What this taught us about the repository and the data.**

- **Phase 1's "one known-bad order-spectrum" is wrong, and the error was
  methodological.** Order 60 was found by chasing an anomalous S/N *minimum*.
  The dominant failure mode — `OPT_MASK` collapsing while the flux is healthy —
  does not move the S/N minimum, so it was invisible to that search. Anything
  found by following up a single summary statistic should be assumed to be the
  tip of something.
- **The dominant defect is `OPT_FRAC_USE`, and phase 1 had the tool but not the
  conclusion.** `redux/checks/frac_use_stats.py` exists and was never turned
  into a claim. 57 of the 68 rejections are mask collapse driven by the object
  profile spilling past the extraction aperture in poor seeing (masked orders:
  median FWHM 2.28 px, FRAC_USE 0.850; unmasked: 1.90 px, 0.988). This is a
  *recipe* question, not a data question — the flux is still in the frame, with
  18,208 counts sitting under a mask that kept 10 pixels of 2048. Prompt 3
  should test whether the aperture or the threshold can be widened.
- **The order numbering runs blue-to-red backwards from intuition.** Order 93
  is 3806 Å and order 57 is 6262 Å. Getting this the wrong way round inverts
  the central conclusion — "the damage is in the blue" versus "the damage is in
  the iodine region" — and the two lead to opposite decisions about phase 3.
  Worth stating explicitly in anything that quotes order numbers.
- **The echelle format masquerades as a defect.** Order 89 departs from a local
  linear fit by +0.062 dex in every single epoch, which is +60σ against the raw
  scatter. Any per-order outlier statistic on these data must be de-trended
  against the format before it is thresholded, or it will confidently report
  the same ten false positives every time.
- **A sigma cut is the wrong tool for this distribution.** After de-trending,
  the residuals are bimodal with a clean gap between 0.015 and 0.035 dex, and
  the robust σ is 0.0010 dex. Thresholding on the physical quantity — a 5%
  departure from the local blaze trend — is defensible; thresholding on σ would
  have flagged everything above the noise floor.
- **Iterative outlier rejection can destroy the signal it is meant to find.**
  Barring pass-one outliers from pass-two neighbourhoods is the standard move
  and it silently un-flagged order 60 of `HI.19980719.49996`, because in that
  frame the whole red block is damaged and there were no clean neighbours left.
  The single pass over-flags by three and keeps the one fault we already knew
  was real. Recorded in the code so it is not "fixed" later.
- **Two more compromised frames were hiding in plain sight.**
  `HI.19980717.48747` and `HI.19980718.49487` sit at 5.1% and 5.6%
  frame-level scatter, between the clean frames (1.4–2.5%) and the known-bad
  one (7.4%). Phase 1 computed the same statistic, saw the largest value, and
  did not look at the ranking.
- **`ascii.csv` does not round-trip booleans.** They come back as the strings
  `'True'`/`'False'`, both truthy. The quality table writes flags as 0/1. Worth
  applying to every table this project commits.

**Files added.**

- `first_hires_exoplanet/quality_filter.py`
- `first_hires_exoplanet/data/order_quality.csv`

**Nothing was changed** in the reductions or in PypeIt.
