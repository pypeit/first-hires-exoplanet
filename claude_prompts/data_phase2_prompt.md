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


### Prompt 3: extending the reduction across the discovery era

**20 of 20 nights reduced, 36 science frames**, every one into 37 echelle
orders (57-93) covering 3806-6262 A. Reproduced by:

```
conda run -n pypeit14 python -m first_hires_exoplanet.calib_inventory
conda run -n pypeit14 python -m first_hires_exoplanet.order_shift
conda run -n pypeit14 python -m first_hires_exoplanet.koa_download --night ...
conda run -n pypeit14 python -m first_hires_exoplanet.reduce_run --era
```

696 MB of raw frames and 9.8 GB of reductions, all in the sibling data tree.
Summary in `redux/run_summary.json`.

#### Does the phase-1 recipe hold over nine months? Yes.

`find_trim_edge = 1,1`, `skip_skysub = True`, `no_local_sky = True` were not
changed for any night. Every night returned **37 orders** except one, and the
wavelength-solution RMS stayed inside the same 0.05-0.24 px band it occupied in
July. Median S/N per night runs 76-182, tracking exposure time and conditions.

#### Does it survive the change of decker? Yes -- but the decker was never the problem.

**1998-09-17 reduced through B2 with the phase-1 parameters untouched**:
37 orders, median S/N 136.9, RMS 0.062 / 0.117 / 0.229, indistinguishable from
the B1 nights. `find_trim_edge = 1,1` was set in phase 1 because the B1 orders
are only ~10.3 binned pixels wide; on the wider B2 slit it is simply permissive
and does no harm.

What 09-17 actually lacked was an **arc**. It took no B2 ThAr, and only three
B2 arcs exist in the entire discovery era (1998-07-14, 09-14, 09-15). It
borrows 1998-09-15's -- 10 s, d echangle 0.0004 -- exactly as 1998-07-19
borrows 07-18's.

#### The calibration problem, which is the real content of this prompt

Phase 1 borrowed flats from 1998-07-14 because the July run had none of its own.
That was not a quirk of July. **Clean B1 flats exist on only five of the
twenty-one discovery-era nights** -- 07-14 (4), 08-12 (18), 08-17 (2), 09-17
(10), 09-18 (8) -- because this program flat-fielded through the wider **B2**
decker as a matter of course while observing through B1. Sixteen nights have
clean B2 flats; most have no usable B1 flat at all.

PypeIt only accepts a flat whose echelle angle is within ~0.01 of the science
frame, and the era spans 0.021, so there is no single flat set for the whole
run. **1998-08-12, with eighteen B1 flats near the middle of the era, is the
calibration keystone for the entire August-September block.** Until this prompt
it was noted only as the cell-in/cell-out diagnostic pair for prompt 7.

The rule adopted, applied uniformly: a night uses its own clean flats of its own
decker when it has at least four; otherwise it borrows four from 1998-08-12.

#### Four traps, each of which broke the reduction before it was found

**1. Five of 1998-08-12's eighteen B1 flats are empty frames.** KOA records
them as `flatlamp`, iodine out, hatch closed -- a perfect clean flat by every
column of the survey. The cross-disperser cover (`XCOVOPEN`) was shut: they hold
bias and read noise, **mean 768 counts with a standard deviation of 1.4**,
against 17,400 and 12,000 for a real one. The "four flats nearest in time to
the science frames" rule selected exactly those four, PypeIt silently declined
to frametype them, and three nights died several minutes later with `No frames
of type=trace provided`, which says nothing whatever about the cause.
**Calibrations cannot be chosen from archive metadata alone.**

**2. The signal test has to be the level, not the scatter.** The obvious
discriminator is spatial structure -- a real B1 flat has bright orders against
dark gaps and reaches `IM01SD01` = 4,300-12,400 where an unlit frame gives 1.4.
But it is decker-dependent: the wider B2 slit fills the frame far more evenly
and its perfectly good flats sit at **149-171**. A scatter threshold calibrated
on B1 throws every B2 flat away. The test is now on `IM01MN01` -- 6,700-20,400
for a lit flat against 763-768 for the bias.

**3. `TARGNAME = 'H187123'` on 1998-08-12.** A one-character typo by the
observer in 1998. KOA's own `targname` column normalises it to `187123`, so the
archive survey and the downloader both find the frames, and only a reader of the
raw FITS header sees it. An exact-equality test gave that night zero science
frames and it was skipped in silence -- the one night that is both the flat
donor for half the era and prompt 7's diagnostic pair. Matched as a substring
now. (The same header shows the cell-out frame labelled `OBJECT = 'Template'`:
the observer was deliberately taking a template.)

**4. December 1997 is a different cross-disperser, as far as PypeIt is
concerned.** `keck_hires.py:245` overrides the header's `XDISPERS = 'RED'` to
`RED97` for every frame before 1997-12-31, on the MAKEE DRP's claim that a
different cross-disperser was fitted. `dispname` is a hard configuration key
with **no tolerance**, so no 1998 frame can ever join a 1997 configuration. The
echelle-angle gap was never the binding constraint.

#### December 1997, and how far the orders actually move

1997-12-23 and 12-24 failed outright at first. The question of whether a flat
from another night could be forced in was settled by measurement rather than by
PypeIt's tolerance, in `order_shift.py`: collapse each night's B1 ThAr arc along
dispersion and cross-correlate the spatial order profile against 1998-07-14's.

| night | d echangle | shift | correlation | PypeIt accepts? |
|---|---|---|---|---|
| **1997-12-23** | +0.01387 | **-1 px** | 0.986 | **no** |
| **1997-12-24** | +0.01488 | **-1 px** | 0.994 | **no** |
| 1998-07-17 | -0.00021 | -3 px | 0.989 | yes |
| 1998-07-18 | -0.00025 | -3 px | 0.985 | yes |
| 1998-08-25 | -0.00512 | -2 px | 0.990 | yes |

**The December nights are better aligned with the July flats than several
nights PypeIt accepts without complaint.** The cross-disperser angle, which is
what sets spatial position, barely moved (d xdangle 0.0025 and 0.0010). The
0.01 echelle tolerance is a proxy for "same instrument setup" and is
conservative here.

That established the geometry, but the `RED97` barrier still rules out any 1998
donor, so December is calibrated **entirely from its own night**. The only B1
quartz flat it has is hatch-closed, cover-open and **iodine-in**, at exactly
the science echelle angle. That is the right frame to trace with: the cell
modulates the spectrum along dispersion without moving the orders.

What an iodine-in flat must never do is flat-field the science, because
dividing it into iodine-in science would partly cancel the 5000-6200 A
absorption phase 3 measures velocities from. The intention was to keep the
illumination correction and drop only the pixel flat; **PypeIt forbids that
combination** (`pypeitpar.py:575` validates that a slit-illumination or spectral
flat-field correction is only applied alongside the pixel flat). So both are
off, and the flat serves **tracing only**.

It works. 1997-12-23 and 12-24 both reduce to 37 orders, median S/N 138.2 and
161.6, RMS 0.067 / 0.121 / 0.208 and 0.070 / 0.121 / 0.221 -- indistinguishable
from the 1998 nights despite `RED97` and trace-only flat-fielding.

**These two nights are nonetheless flat-fielded differently from the other
eighteen**, with no pixel-flat (phase 1 measured it at +/-3%) and no
slit-illumination correction. Both are smooth and multiplicative, so neither
moves a line centre and neither biases a velocity, but anything comparing
December against the rest must know it. The forcing is expressed in
`FORCED_CALIBS` and `DECEMBER_PARAMS` in `reduce_run.py`, not by hand-editing a
`.pypeit` file.

#### Which nights reduce cleanly, and which needed intervention

**Clean, no intervention (16 nights, 31 frames):** 1998-06-18, 07-15, 07-16,
07-17, 07-18, 08-12, 08-17, 08-18, 08-25, 08-26, 09-12, 09-14, 09-15, 09-16,
09-18 -- and of these, eleven borrowed 1998-08-12's flats, which is routine
rather than intervention.

**Needed intervention (4 nights, 5 frames):**

| night | intervention |
|---|---|
| 1998-07-19 | borrows the 07-18 arc (phase 1) |
| 1998-09-17 | borrows the 1998-09-15 **B2** arc; only three B2 arcs exist in the era |
| 1997-12-23 | own iodine-in B1 flat forced as trace; pixel and illumination flats off |
| 1997-12-24 | as 12-23 |

#### Two anomalies for the quality filter

- **1998-09-13** is the one night that did not return 37 orders: **34 orders,
  35 wavelength solutions, worst-order RMS 0.477 px** against ~0.22 everywhere
  else. It should be looked at before it feeds prompt 5.
- **1998-09-17** has an S/N minimum of **13.7**, the signature of a single
  collapsed order -- exactly what prompt 2's filter is built to catch.
- 1998-08-12's low median S/N (76.4, 91.3) is not an anomaly: those are 60 s
  exposures against 200-500 s elsewhere.

`quality_filter.py` was written in prompt 2 against the five July nights and
should now be re-run over all twenty; it reads whatever reductions it finds.

#### What this gives the later prompts

- **Prompt 4** has its template: the three cell-out exposures of 1998-08-26 came
  out at median S/N **180.6, 181.4, 178.4**, the highest of the whole era, plus
  the 1997-12-24 cell-out frame at 161.6 for the consistency check.
- **Prompt 5** has 36 epochs over nine months instead of 10 over five nights.
- **Prompt 7** has the 1998-08-12 cell-in/cell-out pair, both reduced.


### Prompt 4: the iodine-free stellar template

Built by `first_hires_exoplanet/build_template.py`:

```
conda run -n pypeit14 python -m first_hires_exoplanet.build_template
```

Products: `redux/template/hd187123_template_19980826.fits` (37 order
extensions, in the data tree) and `first_hires_exoplanet/data/template_snr.csv`.

#### What was co-added, and in which frame

The three 1998-08-26 cell-out exposures, 500 s each, `IODIN = False` checked in
the raw headers rather than taken from the survey:

| KOAID | MJD | v_bary |
|---|---|---|
| `HI.19980826.34749` | 51051.402204 | −7.2433 km/s |
| `HI.19980826.35345` | 51051.409098 | −7.2588 km/s |
| `HI.19980826.35939` | 51051.415978 | −7.2740 km/s |

Following prompt 1, every wavelength was first divided by `VEL_CORR` to undo
PypeIt's heliocentric correction and recover the **observed** frame, which is
what `pyodine` expects. The three exposures span 20 minutes, over which the
barycentric velocity moves **−15.5 to −30.7 m/s**; each was shifted onto the
first exposure's frame before co-adding. That is 0.5% of a resolution element
and cosmetic, but this is the one spectrum every later velocity is measured
against.

`barycorrpy` is not installed in `pypeit14`, so `astropy` was used. Prompt 1
established the two agree to a near-constant 4.6 m/s, which is irrelevant for
the *relative* alignment of three exposures twenty minutes apart. Phase 3
should still install `barycorrpy`, as prompt 1 recommended.

#### Signal-to-noise achieved

**Median S/N 309.5 per pixel, range 103.8 to 345.6**, against 62.1–198.1 in a
single exposure. The gain is **1.72** against the ideal √3 = 1.73 — the
co-addition is doing essentially everything it can.

| order | λ (Å) | single | template | gain |
|---|---|---|---|---|
| 57 | 6210–6262 | 197.0 | 342.1 | 1.74 |
| 60 | 5899–5984 | 197.8 | 345.3 | 1.75 |
| 65 | 5445–5523 | 194.9 | 337.8 | 1.73 |
| 70 | 5056–5129 | 188.2 | 324.6 | 1.72 |
| 75 | 4719–4786 | 180.6 | 309.5 | 1.71 |
| 80 | 4424–4488 | 162.1 | 275.9 | 1.70 |
| 85 | 4164–4224 | 138.2 | 232.4 | 1.68 |
| 89 | 3977–4034 | 115.5 | 193.8 | 1.68 |
| 90 | 3933–3989 | 88.5 | 143.1 | 1.62 |
| 93 | 3806–3861 | 62.1 | 103.8 | 1.67 |

All 37 orders have all three exposures contributing at every good pixel. The
full table is in `data/template_snr.csv`. Order 90 is the weakest gain, 1.62 —
that is the Ca II H&K order, where the deep cores leave few high-flux pixels.

The iodine region that phase 3 needs, orders 57–71, sits at **S/N 320–346**.

#### Check 1: the cell really was out

Absorption minima per 100 Å, prominence > 4% of the continuum, against a
cell-in frame of the **same night** (`HI.19980826.44043`):

| window | template | cell-in |
|---|---|---|
| iodine, 5000–6200 Å | **91.7** | **313.0** |
| control, 4000–4800 Å | 229.5 | 231.9 |

Decisive. In the control window, where the cell has no lines, the two agree to
1%. In the iodine window the cell-in frame has **3.4×** the line density. The
template is genuinely iodine-free, and the comparison is internally controlled
rather than resting on a header keyword.

#### Check 2: a G2V expectation

All ten canonical G-dwarf features are present at sensible depths:

| line | order | depth |
|---|---|---|
| Ca II K | 91 | 0.801 |
| Ca II H | 90 | 0.730 |
| H-delta | 87 | 0.745 |
| H-gamma | 82 | 0.778 |
| H-beta | 73 | 0.823 |
| Mg b1/b2/b3 | 69 | 0.848 / 0.876 / 0.897 |

**Velocity zero point.** These wavelengths are in the observed frame, so the
lines are *not* at the −17.0 km/s systemic velocity — they sit at
v_sys − BVC = **−9.76 km/s**. Measured median: **−8.43 km/s**, residual
**+1.33 km/s**. Phase 1 found the same sign and size (−14.7 against −17.0, i.e.
+2.3 km/s) on a July frame, so this is the known zero-point offset rather than
a new one. It is ~0.2 of a resolution element and comes from taking the deepest
pixel of deep, asymmetric line cores.

*Comparing against −17 km/s would have manufactured an 8 km/s error that is not
there.* The first version of this check did exactly that.

**Resolution.** The named Fe I lines give FWHM 13.2–22.3 km/s, but they are all
strong and saturated with damping wings, so that is a stellar width, not an
instrumental one. Measuring the **weak** lines instead — 91 features of depth
0.10–0.40 in orders 78–81, which prompt 2 found clean in every epoch:

| percentile | FWHM |
|---|---|
| 10th | 7.45 km/s |
| 25th | 8.01 km/s |
| 50th | 9.48 km/s |

The narrowest features bound the resolution element from above, giving
**R ≳ 40,200** against the ~45,000 expected for the B1 decker. Since weak lines
still carry intrinsic width, that is a lower bound and it is consistent. It
also rules out a giant or a fast rotator.

#### Check 3: against the 1997-12-24 cell-out frame

1997-12-24's only science frame happens to be cell-out, giving an independent
comparison **eight months earlier**, on the other side of the `RED97`
cross-disperser boundary and reduced with the trace-only flat-fielding of
prompt 3.

- 31 orders compared
- correlation: **median 0.992**, minimum 0.882
- velocity offset: **+0.16 km/s, scatter 0.07 km/s** across orders

The two epochs are the same star, reduced consistently, over an eight-month
baseline and across two different calibration routes. That 160 ± 70 m/s is
**not** a planet measurement — it is far above the 72 m/s semiamplitude and is
dominated by method systematics. It is, however, a useful early number: it is
the first direct measurement in this project of the floor a cross-correlation
reaches, and it lands exactly where context Q3 predicted, in the "tens to
hundreds of m/s" band.

#### Three corrections made along the way

Recorded because each was wrong in a way that looked plausible:

1. **The velocity expectation.** Comparing observed-frame lines against the
   systemic velocity, forgetting the −7.24 km/s barycentric term the frame
   still carries.
2. **The line-width estimator.** Taking min-to-max of everything above half
   depth measures the *window* rather than the line as soon as a blend enters
   it: Fe I 4383 returned 111 km/s, which is 1.6 Å, the full window. Fixed by
   walking out from the minimum along the contiguous run above half depth and
   rejecting lines that never come back down.
3. **The cross-correlation.** Every order returned exactly +0.00 km/s with zero
   scatter — a suspiciously perfect result that was pure quantisation: the grid
   sampled ~2 km/s per lag and the true offset is smaller. Parabolic refinement
   of the peak turns it into +0.16 ± 0.07 km/s.

#### What phase 3 still needs from the template

- `pyodine` wants the template as a `StellarTemplate`, deconvolved against the
  instrumental profile. This product is the co-added observation, not yet
  deconvolved; that is prompt 8's and phase 3's work.
- The template's own `bary_vel_corr` and `bary_date` are stored in the primary
  header (`VBARY`, `MJDREF`) in the conventions prompt 1 settled: observed
  frame, astropy barycentric at the exposure midpoint.


### Prompt 5: relative velocities by cross-correlation

Measured by `first_hires_exoplanet/measure_velocities.py`:

```
conda run -n pypeit14 python -m first_hires_exoplanet.measure_velocities
```

Tables: `data/xcorr_velocities.csv` (36 epochs) and
`data/xcorr_velocities_per_order.csv` (the 717 individual order velocities
behind them).

Prompt 2's quality filter was first re-run over all twenty nights. Its file
glob had been written as `reduce_1998*` and silently excluded the two 1997
December nights; corrected to `reduce_19*`, it now covers all 36 frames and
1329 order-spectra, of which **719 of 790 blue order-spectra are usable**.

#### Method

Each epoch is cross-correlated against the prompt-4 template over the 22 orders
lying entirely blueward of 5000 Å (orders 72–93), on a common log-wavelength
grid at 1.5 km/s per pixel, searching ±60 km/s, with parabolic refinement of
the peak. Only order-spectra prompt 2 marks usable take part.

Both spectra are in the **observed** frame, `VEL_CORR` divided out per prompt 1,
so PypeIt's heliocentric correction — with the +13 m/s sign error and 5 m/s
drift prompt 1 found — never enters. In that frame a line sits at v_star − BVC,
so the correlation returns dv_obs = (v_e − v_t) − (BVC_e − BVC_t), and the
relative barycentric velocity is dv_obs + (BVC_e − BVC_t), with BVC from astropy
at the exposure midpoint.

**Internal check:** the three exposures that *are* the template return
**+2.1, +17.0, −14.2 m/s**. They should return zero, and do, to 17 m/s.

#### The two scatters

They answer different questions and conflating them would be the whole error of
this prompt.

**PER-ORDER scatter, within one exposure** — *do the 22 blue orders agree with
each other?*

> **median 47 m/s**, range 13–156 m/s

This is the precision of a single order and the internal consistency of the
wavelength solution across the format. It is dominated by photon noise and by
how well each order's dispersion solution agrees with its neighbours'.

**EPOCH-TO-EPOCH scatter, across the baseline** — *does the same star give the
same velocity on different nights?*

> **rms 702 m/s** over 36 epochs and 269 days
> (488 m/s excluding 1998-07-19, the one night with no arc of its own)

This is the end-to-end precision of the method, and it includes everything the
per-order scatter cannot see.

**The gap between them is the result.** 47 m/s within an exposure against
702 m/s between exposures means the limiting error is *not* photon noise and
*not* the quality of individual orders — those are fine. It is something common
to every order of a given exposure: the **wavelength zero point**. All 22 orders
of one exposure agree with each other and then move together, by up to a
kilometre per second, from one exposure to the next.

That is precisely the error an iodine cell removes, by printing a wavelength
fiducial onto the same photons through the same optics at the same instant.
**This is the argument for phase 3, now measured rather than asserted.**

#### Where the 702 m/s comes from

| | rms |
|---|---|
| within a night (8 nights with >1 frame) | median **132 m/s**, worst spread **1604 m/s** |
| between nights (20 nightly means) | **605 m/s** |

So it is mostly a night-to-night effect, with a substantial intra-night
component. Two nights show large drifts across a single night: 1998-08-25
spans **1604 m/s over 6.7 hours** and 1998-08-26 **627 m/s over 6.9 hours** —
consistent with instrument flexure and with a wavelength solution anchored to an
arc taken at one end of the night.

#### A phase-1 conclusion that needs correcting

**1998-07-19 took no calibration frames at all and borrows the 07-18 arc.** Its
three frames come out at **+2438, +1609, +1554 m/s**, against 0–400 m/s for
every other July night. Excluding it drops the overall rms from 702 to 488 m/s.

Phase 1 concluded that *"same-night arcs are not required"* because 07-19's
wavelength-solution RMS was indistinguishable from the other nights'. That is
true and beside the point: **the RMS measures the scatter of the fit, not its
zero point.** A solution can be beautifully tight and sit a kilometre per
second off. Any future night calibrated from a borrowed arc should be treated
as having an unknown velocity offset until it is measured.

#### Uncertainties, and why the small ones are not the honest ones

| | median |
|---|---|
| formal (Zucker 2003, combined over orders) | **2.9 m/s** |
| empirical (per-order scatter / √N) | **10.4 m/s** |
| actual epoch-to-epoch rms | **702 m/s** |

Both per-epoch error bars are optimistic by a factor of 70–240. Both are
reported per epoch in the table, and neither should be used as the uncertainty
on a velocity. The honest statement is that a single epoch's velocity is good to
roughly **0.5–0.7 km/s**, set by the epoch-to-epoch scatter, and the formal
errors describe only the internal consistency of one exposure.

This is exactly why the prompt asked for the two scatters separately.

#### The sign error, recorded because it looked like a result

The first run produced velocities running smoothly from **+37.4 km/s in June to
−10.9 km/s in September** — a clean annual curve. It was a sign error in the
cross-correlation lag: a positive lag aligns the epoch feature at the *shorter*
wavelength, so a positive lag means the epoch is blueshifted and dv_obs is
negative. With the sign reversed the barycentric term was doubled instead of
cancelled, and the output was the Earth's orbit rather than the star's velocity.

It did not look like a bug. It looked like a smooth, physically shaped,
high-precision measurement. The tell was magnitude: a 72 m/s planet cannot
produce 37 km/s. The same error was present in prompt 4's
`compare_against` and has been fixed there too; it changed only the sign of an
already-small residual, since both spectra were put in a common frame first.

#### What this gives the later prompts

- **The per-epoch initial guess `pyodine` needs**, good to ~0.5 km/s, which is
  far better than that code requires.
- **The number for the argument**: cross-correlation on this data reaches
  **~700 m/s**, or ~490 m/s discounting the arc-less night, against a **72 m/s**
  semiamplitude. The method is an order of magnitude too coarse to see the
  planet, as context Q3 and the phase-1 assessment both predicted.
- **A quality signal prompt 2 could not give**: 1998-07-19's borrowed arc, and
  the intra-night drifts on 08-25 and 08-26, are velocity-level problems
  invisible to an S/N-based filter.
- Prompt 6 should compare these against the published Keplerian and the modern
  catalogue, and state plainly what 700 m/s against 72 m/s does and does not
  demonstrate.


### Prompt 6: assessment against the Keplerian and the modern catalogue

`first_hires_exoplanet/figs_phase2.py` does the assessment and writes four
figures into `docs/figs/`:

```
conda run -n pypeit14 python -m first_hires_exoplanet.figs_phase2
```

| figure | what it shows |
|---|---|
| `fig_p2_timeseries.png` | our velocities and Teklu et al.'s over the same nine months, each on its own scale |
| `fig_p2_phasefold.png` | both folded on 3.0966 d: the catalogue traces the Keplerian, ours does not |
| `fig_p2_compare.png` | epoch by epoch on one scale, equal aspect |
| `fig_p2_precision.png` | what each stage reaches against what the planet requires |

#### The epochs match one to one

All **30** discovery-era catalogue rows pair with one of our epochs, worst time
difference **0.1 minutes**. Six of our 36 have no counterpart, and they are
exactly the ones that should not:

- the five **cell-out** frames (1997-12-24, 1998-08-12's template frame, and
  the three 1998-08-26 template exposures) — Butler's method needs the cell, so
  these were never velocity epochs;
- the one 60 s **cell-in** frame of 1998-08-12.

**That settles an open question this document has carried since it was
written**: the survey's 32 cell-in frames against 30 catalogue velocities were
indeed explained by the two short 60 s exposures. One of them is 1998-08-12's,
which we reduced and which has no catalogue velocity; the other is
1998-07-14's, which we never reduced.

#### Against the modern catalogue

| | rms |
|---|---|
| Teklu et al. 2025, discovery era | **47.4 m/s** |
| this reduction, same epochs | **723.5 m/s** |
| difference | **706.7 m/s** |

Our scatter is **15×** the catalogue's, so the difference is essentially our own
noise: their signal is buried inside it. The correlation is r = **+0.38**. If we
were measuring their signal plus independent noise of our own size we would
expect r = +0.07, and the standard error on r with 30 points is 0.19, so +0.38
sits about 1.5 standard errors above that expectation. **That is not evidence of
anything.** With 30 points it is the kind of correlation noise produces
routinely, and it should not be presented as partial recovery of the signal.

#### Against the published 3.097-day Keplerian

Fitting a circular Keplerian with the period forced to 3.0965828 d — linear in
its three parameters, so there is no minimiser to get stuck:

| | K | residual rms | significance |
|---|---|---|---|
| Teklu et al. 2025 | **69.2 ± 0.3 m/s** | 2.2 m/s | 211σ |
| this reduction | **408 ± 188 m/s** | 668 m/s | 2.2σ |

The catalogue recovers the published 72 m/s cleanly, which confirms the fit is
working and the period is right.

**Our 408 m/s is not a detection and must not be read as one.** It is 5.7× the
published value and differs from it by 1.8σ. Forcing the fit removes almost
nothing — the residual rms falls only from 724 to 668 m/s — so the fit is
absorbing noise, not signal. The honest statement is a **95% upper limit of
785 m/s**, which sits **11× above** the semiamplitude we are trying to see. Our
data are consistent with the published planet and equally consistent with no
planet at all.

#### What this does and does not demonstrate

**It does demonstrate that the reduction is sound end to end.** 36 epochs over
269 days, every one yielding a velocity. The three exposures that are themselves
the template return zero to 17 m/s. The 22 blue orders of a single exposure
agree to 47 m/s. The epochs match the archive's own record of the same
observations to a tenth of a minute. Nothing in the chain from raw frame to
velocity is broken.

**It does not detect the planet, and it does not confirm the published orbit.**
At 724 m/s against a 72 m/s semiamplitude the measurement has no power to do
either. A non-detection was the predicted outcome — context Q3 and the phase-1
assessment both put cross-correlation in the tens-to-hundreds of m/s band — and
a non-detection is what we got. **Nothing in this result would have looked
different if HD 187123 b did not exist.** Any presentation of phase 2 that
implies otherwise would misrepresent it.

It is worth being equally clear about what the *catalogue* comparison does not
show. Teklu et al. reduced the same photons with the iodine method and reach
1.2 m/s. That is not a statement about PypeIt against their pipeline in general;
it is a statement about a wavelength solution from ThAr arcs against one printed
onto the science photons by an iodine cell.

#### What the achieved precision implies for phase 3

To detect K = 72 m/s at 3σ per epoch, precision must improve by about **30×**,
from 724 m/s to roughly 24 m/s.

That sounds forbidding until the decomposition from prompt 5 is applied:

- the **per-order scatter within one exposure is 47 m/s**, already below the
  semiamplitude;
- the **epoch-to-epoch scatter is 724 m/s**, fifteen times larger.

The photons are not the limit. The orders of a single exposure agree with each
other and then move together, which localises the entire deficit in the
**wavelength zero point** — the one quantity an iodine cell fixes, by printing a
fiducial onto the same photons through the same optics at the same instant.
Teklu et al. reach 1.2 m/s on these very frames, so the data carry the signal;
only our wavelength calibration does not.

**This is the strongest possible argument for phase 3, and it is now measured
rather than asserted.** Phase 2 has done what the document said it was for.

Two specific consequences for phase 3:

1. **The initial guess is already good enough.** `pyodine` needs a per-epoch
   starting velocity; ours are good to ~0.7 km/s, far inside its tolerance.
2. **Two nights carry known velocity-level defects** that an S/N filter cannot
   see: 1998-07-19, calibrated from a borrowed arc, sits 1.5–2.4 km/s off, and
   1998-09-17 is on a borrowed arc too. A forward model should recover both,
   since it derives its own wavelength solution from the I2 lines — which makes
   them a useful test of whether phase 3 is working.


### Prompt 7: the atlas question, settled

`first_hires_exoplanet/atlas_check.py`, with `docs/figs/fig_p2_atlas.png`:

```
conda run -n astro python -m first_hires_exoplanet.atlas_check
```

**Environment note.** `h5py` is not installed in `pypeit14`. Rather than change
that environment, this module avoids importing `pypeit` — it reads the spec1d
files with plain `astropy.io.fits`, whose layout phase 1 and prompt 4 already
documented — so it runs in `astro`, which the repository already uses for
`figs.py`. Phase 3 will need `h5py` in `pypeit14`, since `pyodine` requires it.

The atlas is 55.8 MB and lives in `../first-hires-exoplanet-data/atlas/`, not
the repository. Obtained from
`raw.githubusercontent.com/pepeheeren/pyodine/main/iodine_atlas/` at commit
`4488b0914fe5b272b787982647691045bff2604a`.

#### Provenance: what the file actually is

**The file carries no metadata whatsoever** — no HDF5 attributes at root or on
any dataset. Everything below is measured from the data.

| | |
|---|---|
| datasets | `flux`, `flux_normalized`, `wavelength`, `wavelength_air`, `wavenumber` |
| points | 1,551,029 |
| coverage | **4980.01–6250.00 Å vacuum** (4978.62–6248.27 air) |
| wavenumber | 16000.01–20080.29 cm⁻¹, step 0.0025970 cm⁻¹ |
| grid | **uniform in wavenumber** — the signature of an FTS, not a grating |
| resolution | **R ≈ 600,000–680,000** |

The coverage is exactly the I2 band and nothing else, which is consistent with
a dedicated cell scan. The resolution is measured from 4653 isolated I2 lines
in 5400–5600 Å: FWHM 0.440 km/s at the 10th percentile, 0.500 km/s median.

Against that, the **thermal Doppler width** of I2 is 0.242 km/s at 50 °C,
0.248 at 65 °C, 0.251 at 75 °C. The narrowest atlas lines are within a factor
of two of the thermal limit, so the atlas genuinely resolves the forest. The
temperature cannot be pinned from the widths — the FTS profile dominates, and
50 °C and 75 °C are indistinguishable at 0.44 km/s.

**The atlas resolution is not a limitation for us.** HIRES runs at R ≈ 42,000
(measured in prompt 4); the atlas is fifteen times finer.

#### The empirical test, which is the point of this prompt

The 1998-08-12 pair — `HI.19980812.29160` cell in, `HI.19980812.29316` cell
out, 156 seconds apart — divides out the star, the blaze and the detector and
leaves the **transmission of the HIRES cell itself**, at HIRES resolution. Both
frames are reduced (prompt 3). 14 orders in the I2 band give a usable ratio.

Comparing against the atlas convolved to R = 42,000, with an *identical*
running-95th-percentile continuum applied to both sides:

| | |
|---|---|
| median correlation, measured vs atlas | **0.928** |
| median measured depth | 0.122 |
| median atlas depth | 0.044 |

**The line positions transfer. The depths do not.** This is exactly the split
the document predicted: positions are molecular constants, depths are column
density and temperature.

#### How much deeper, and is it one number?

A mean-depth ratio is the wrong statistic because it ignores saturation. The
physical comparison is Beer–Lambert: if the cells differ only in column
density then T_HIRES = T_atlas^α, and α is the column-density ratio. Fitting α
on the logarithms is a one-parameter linear fit.

| fit range | α | spread | n |
|---|---|---|---|
| full (0.02–0.97) | 2.79 | 2.26–3.53 | 13 |
| restricted (0.50–0.90) | **2.59** | **2.28–2.86** | 7 |

Over the full range α trends with wavelength — 3.5 at 6050 Å down to 2.3 at
5320 Å — and is largest exactly where the absorption is weakest. That is the
continuum, not the cell: where the forest is dense a running 95th percentile
cannot find true continuum. Restricting to transmissions of 0.50–0.90, where
neither saturation nor continuum placement bites, the spread collapses.

**A single exponent α ≈ 2.6 describes the difference.** The two cells differ in
column density, not in a way that reshuffles relative line strengths — so
temperature is not importantly different either.

*Noise control:* the measured ratio comes from two 60 s frames, and a running
upper-percentile continuum rides up on noise, deepening lines. Injecting the
same 1.7% fractional noise into the atlas and pushing it through the identical
normalisation returns α = 1.07. The noise accounts for 7%, not a factor of 2.6.
**The depth difference is real.**

`fig_p2_atlas.png` shows it directly: as shipped, the atlas is line-for-line
aligned with the measured HIRES cell and far too shallow; raised to the power
2.6 it tracks it almost exactly.

#### Can it serve as the HIRES cell's atlas?

**As shipped, no. Transformed, yes — but `pyodine` cannot do the transformation
it offers.**

`pyodine` carries a free `iod_depth` parameter (`models/spectrum.py:17`), which
looks like it settles the matter. It does not, because of how it is applied
(`models/spectrum.py:93`):

```python
flux_iod = params['iod_depth'] * (flux_iod - 1.0) + 1.0
```

That is a **linear** depth scaling, T → 1 + d(T − 1). The physical law is
Beer–Lambert, T → T^α. The two agree only for weak lines. This atlas reaches
`flux_normalized` = 0.000, and a saturated line at T = 0 scaled linearly by
d = 2.6 returns **T = −1.6**: negative transmission. **5.8% of the atlas would
go negative** under the scaling we need.

So three routes, in order of preference:

1. **Change `pyodine`'s depth model from linear to power-law** — a one-line
   change, `flux_iod ** params['iod_depth']`, which is physically correct,
   reduces to the current behaviour for weak lines, and lets the code fit α
   itself. This is the recommendation.
2. **Pre-transform the atlas** with T → T^2.6 before handing it over, and leave
   `iod_depth` free to take up the residual. Works with `pyodine` unmodified,
   but bakes in a number measured from two 60 s frames.
3. **Obtain an FTS scan of the actual HIRES cell.** The real answer, and worth
   asking for — the Keck/HIRES cell was scanned for the Butler et al. (1996)
   method, but no public copy has been found. This is the one item here that
   needs someone outside the project.

#### The air/vacuum trap, which would have been fatal

`pyodine`'s `IodineTemplate` reads **`wavelength_air`**
(`utilities_lick/load_pyodine.py`). PypeIt reports **vacuum**. The two grids
differ by 1.56 Å at 5615 Å, which is **83 km/s**.

Measured on order 65, against the HIRES cell transmission:

| atlas grid | correlation |
|---|---|
| vacuum | **+0.903** |
| air | **+0.078** |

The signal disappears entirely. This project has already been caught once by
air versus vacuum — phase 1's spurious +70 km/s — and this is the same trap in
a new place. The prompt-8 adapter must read `wavelength`, not `wavelength_air`.

#### Answers to the three questions asked

1. **Can the Fischer atlas serve as the HIRES cell's atlas?** Its line
   positions can, at correlation 0.93. Its depths cannot: the HIRES cell is
   2.6× optically thicker. The fix is a power-law rescaling, which `pyodine`'s
   linear `iod_depth` cannot express.
2. **Can the 1998-08-12 pair test that empirically?** Yes, and it has. 14
   orders, a decisive answer on both positions and depths, with a noise control
   that rules out the artefact. Nothing else in this dataset could have settled
   it, and no header keyword could.
3. **If not, what would we need instead?** An FTS scan of the HIRES cell. Not
   required to proceed — routes 1 and 2 above are sound — but it would remove
   the one remaining assumption, that a single Beer–Lambert exponent captures
   the whole difference between two cells whose temperatures we cannot measure.


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

### 2026-09-22 (Prompt 3: the whole discovery era reduced -- 20 nights, 36 frames -- after four calibration traps)

**Task.** Prompt 3 of this document: extend the reduction to the rest of the
discovery era, reusing `reduce_run.py`; report which nights reduce cleanly,
which need intervention, and whether the phase-1 recipe holds across nine months
and a change of decker. Full findings are in the Report section above.

**What was done.**

- Wrote `first_hires_exoplanet/calib_inventory.py`: per-night inventory of what
  calibrations exist in the planet-search configuration, decker-aware, plus a
  check of the flats actually on disk against what the archive claims.
- Wrote `first_hires_exoplanet/order_shift.py`: measures the spatial order shift
  between nights by cross-correlating collapsed B1 arc profiles, which is what
  decided the December question.
- Extended `koa_download.py` with the fifteen new nights, and `reduce_run.py`
  with per-night flat selection, forced frametypes and per-night parameters.
- Downloaded 696 MB of raw frames; produced 9.8 GB of reductions.

**Headline results.** 20 of 20 nights, 36 science frames, 37 orders each,
wavelength RMS in the same 0.05-0.24 px band as July. The recipe held unchanged
across nine months and across the B1-to-B2 decker change. Four nights needed
intervention: 07-19 and 09-17 borrow arcs, 12-23 and 12-24 are traced from their
own iodine-in flat with pixel and illumination flatting switched off.

**What this taught us about the repository and the data.**

- **Phase 1's flat problem was the general case, not a quirk of July.** Clean B1
  flats exist on only five of twenty-one nights; this program flat-fielded
  through B2 while observing through B1. 1998-08-12, previously noted only as
  the cell-in/cell-out pair for prompt 7, turns out to be the calibration
  keystone for the whole August-September block.
- **Archive metadata is not enough to choose calibrations.** Five of 08-12's
  eighteen B1 "flats" are clean by every KOA column and contain nothing but
  bias -- the cross-disperser cover was shut. They are distinguishable only from
  the pixel statistics, and PypeIt's downstream error (`No frames of type=trace
  provided`) points nowhere near the cause. Any future calibration selection in
  this project should verify the frames, not just the table.
- **A discriminator calibrated on one decker can be wrong on another.** The
  first version of that check tested spatial scatter, which separates lit from
  unlit B1 flats by four orders of magnitude -- and throws away every B2 flat,
  because the wide slit fills the frame evenly (sigma 149-171). Testing the
  level instead works for both.
- **Exact string matches on observer-typed header cards are fragile.**
  `TARGNAME = 'H187123'` on 1998-08-12 -- one stray character in 1998 -- silently
  removed the most important night in the era from the reduction. KOA's
  normalised `targname` hides it, so the survey and the downloader both work and
  only the raw header reveals it.
- **`dispname` is a hard configuration key and PypeIt rewrites it by date.**
  `keck_hires.py:245` turns `XDISPERS = 'RED'` into `RED97` for everything before
  1997-12-31 on the MAKEE DRP's authority. No 1998 calibration can ever reach a
  1997 night. This was the binding constraint on December, not the echelle angle
  we spent time measuring -- though that measurement was still worth having.
- **PypeIt's configuration tolerances are conservative, and that is now
  quantified.** The December nights sit 1 binned pixel from the July flats at
  correlation 0.99, closer than 1998-07-17 and 07-18 which PypeIt accepts at
  -3 px. `order_shift.py` is the general tool for asking this question of any
  future night.
- **PypeIt will not let you apply an illumination flat without a pixel flat**
  (`pypeitpar.py:575`). The plan for December was to keep the illumination
  correction and drop only the pixel flat, to keep the iodine-in flat's I2 out
  of the science frames; that combination is rejected outright, so both had to
  go and the flat serves tracing only.
- **A per-night try/except is essential in a twenty-night driver.** The first
  run aborted the whole loop on the first failure. It also exposed a formatting
  bug -- `'{:s}'.format(exc)` raises on an exception object -- which hid the real
  error behind a `TypeError`.
- **Two anomalies to carry forward:** 1998-09-13 returned 34 orders and a
  worst-order RMS of 0.477 px, and 1998-09-17 has an S/N minimum of 13.7.
  `quality_filter.py` from prompt 2 should now be re-run over all twenty nights.

**Files added.**

- `first_hires_exoplanet/calib_inventory.py`
- `first_hires_exoplanet/order_shift.py`

**Files modified.** `first_hires_exoplanet/koa_download.py` (fifteen new nights,
the 09-17 B2 arc, the December trace flats, `xcovopen` in the flat selection);
`first_hires_exoplanet/reduce_run.py` (the era nights, header-driven frame
classification, per-night flat selection with the donor rule, forced frametypes
and parameters for December, per-night error handling, decker in the summary).

**No git command changed state, and no raw frame or reduction product is in the
repository.**

### 2026-09-22 (Prompt 4: the stellar template -- S/N 310 median, iodine-free, and three broken checks fixed)

**Task.** Prompt 4 of this document: co-add the three cell-out exposures of
1998-08-26 into a single high-S/N iodine-free spectrum of HD 187123, order by
order; check it against the 1997-12-24 cell-out frame and against a G2V
expectation; report the S/N achieved per order. Full findings are in the Report
section above.

**What was done.**

- Wrote `first_hires_exoplanet/build_template.py`: loads the three exposures in
  the observed frame (dividing out `VEL_CORR` per prompt 1), shifts each onto
  the first exposure's barycentric frame, inverse-variance co-adds per order,
  and runs three independent checks.
- Wrote `redux/template/hd187123_template_19980826.fits` (37 order extensions,
  in the data tree) and `first_hires_exoplanet/data/template_snr.csv`.

**Headline results.** Median S/N 309.5 per pixel, range 103.8-345.6, gain 1.72
against the ideal sqrt(3) = 1.73. The iodine region phase 3 needs sits at
320-346. The cell-out claim is confirmed internally: 91.7 absorption minima per
100 A in 5000-6200 A against 313.0 for a cell-in frame of the same night, with
a control window at 4000-4800 A agreeing to 1%. Against the 1997-12-24 cell-out
frame, eight months earlier, 31 orders correlate at median 0.992 with a
velocity offset of +0.16 +/- 0.07 km/s.

**What this taught us about the repository and the data.**

- **The observed frame moves the expected line velocity, and it is easy to
  forget.** Having removed PypeIt's heliocentric correction, the lines sit at
  v_sys - BVC = -9.76 km/s, not at the -17.0 km/s systemic velocity. The first
  version of the check compared against -17 and reported an 8 km/s error that
  did not exist. Any check of a wavelength zero point in this project has to
  state which frame it is in first. The real residual, +1.33 km/s, matches the
  +2.3 km/s phase 1 found on a July frame -- a known offset, not a new one.
- **A "perfect" result is a bug report.** The cross-correlation against the
  December frame returned exactly +0.00 km/s with 0.00 scatter for all 31
  orders. That was not agreement, it was quantisation: the grid sampled ~2 km/s
  per lag and the true offset is smaller than one lag. Sub-pixel refinement
  turns it into +0.16 +/- 0.07 km/s. Zero scatter across 31 independent
  measurements should never be believed.
- **Strong lines do not measure resolution.** The Fe I lines chosen for a width
  check are all saturated with damping wings and gave 13-22 km/s, which says
  nothing about the instrument. Weak unsaturated lines (depth 0.10-0.40) give a
  10th-percentile FWHM of 7.45 km/s, R >~ 40,200 against the ~45,000 expected.
  Worth remembering when any line-profile diagnostic is written later.
- **An envelope estimator needs to know where the line is.** Taking min-to-max
  of everything above half depth measures the window, not the line, as soon as
  a blend enters: Fe I 4383 returned 111 km/s, exactly the 1.6 A window width.
  Walking outward from the minimum along the contiguous run above half depth,
  and rejecting profiles that never come back down, fixes it and also rejects
  the blends.
- **The cell-out check should be internally controlled.** Counting I2 lines in
  the template alone proves nothing -- a low count could mean a poor spectrum.
  Counting them against a cell-in frame of the same night, with a blue control
  window where the cell has no lines, turns it into a real test: 1% agreement
  in the control, 3.4x difference in the iodine window.
- **The eight-month cell-out comparison is the first real precision number this
  project has.** 160 +/- 70 m/s between two independent epochs, across the
  RED97 cross-disperser boundary and two different flat-fielding routes, lands
  exactly in the "tens to hundreds of m/s" band context Q3 predicted for
  cross-correlation, and is well above the 72 m/s semiamplitude. It is the
  argument for phase 3, measured rather than asserted, and prompt 5 should be
  expected to land near it.
- **1997-12-24's only science frame is cell-out**, which is why the consistency
  check was available at all. Worth noting that the December nights contribute
  a template check rather than a velocity epoch.

**Files added.**

- `first_hires_exoplanet/build_template.py`
- `first_hires_exoplanet/data/template_snr.csv`
- `redux/template/hd187123_template_19980826.fits` (data tree, not the repo)

**No git command changed state.**

### 2026-09-22 (Prompt 5: cross-correlation velocities -- 47 m/s within an exposure, 702 m/s between them)

**Task.** Prompt 5 of this document: measure relative velocities by
cross-correlation against the prompt-4 template, restricted to the orders
blueward of 5000 A, with an honest uncertainty, weighted per order; write the
table to `first_hires_exoplanet/data/`; report the per-order scatter and the
epoch-to-epoch scatter separately. Full findings are in the Report section
above.

**What was done.**

- Re-ran `quality_filter.py` over all twenty nights after fixing its file glob,
  which had been written as `reduce_1998*` and silently excluded the two 1997
  December nights. Now 36 epochs, 1329 order-spectra, 719 of 790 blue ones
  usable.
- Wrote `first_hires_exoplanet/measure_velocities.py`: log-wavelength
  cross-correlation per order against the template, both sides in the observed
  frame, barycentric correction applied explicitly, Zucker (2003) formal errors
  alongside empirical ones.
- Wrote `data/xcorr_velocities.csv` (36 epochs) and
  `data/xcorr_velocities_per_order.csv` (717 order velocities).

**Headline results.** Per-order scatter within an exposure: median 47 m/s.
Epoch-to-epoch scatter across 269 days: 702 m/s, or 488 m/s excluding the one
night with no arc of its own. The three exposures that are themselves the
template return +2.1, +17.0 and -14.2 m/s, which is the internal check passing.

**What this taught us about the repository and the data.**

- **The gap between the two scatters is the scientific result, not a
  disappointment.** 47 m/s within an exposure and 702 m/s between exposures
  means the orders agree with each other and then move together. The limiting
  error is not photon noise and not per-order wavelength quality -- it is the
  wavelength zero point, common to all 22 orders of an exposure. That is
  exactly the error an iodine cell removes, so phase 2 has now produced the
  quantitative argument for phase 3 that the document asked for.
- **A wavelength-solution RMS says nothing about its zero point, and phase 1's
  conclusion needs correcting.** Phase 1 wrote that same-night arcs are not
  required, because 1998-07-19's RMS was indistinguishable from nights with
  their own arcs. Its three frames sit 1.5-2.4 km/s off every other July night.
  The fit is tight and in the wrong place. Any night calibrated from a borrowed
  arc -- 07-19, and 09-17 from prompt 3 -- carries an unknown velocity offset
  until it is measured.
- **A sign error can look exactly like a result.** The first run returned
  velocities sweeping smoothly from +37.4 km/s in June to -10.9 km/s in
  September: a clean, physically shaped annual curve with small error bars. The
  cross-correlation lag had the wrong sign, so the barycentric term was doubled
  rather than cancelled and the output was the Earth's orbit. Nothing about the
  shape gave it away; only the magnitude did, because a 72 m/s planet cannot
  produce 37 km/s. Worth a standing habit: check the amplitude of any new
  measurement against what the physics allows before looking at its shape.
- **Formal errors on a cross-correlation are not the precision.** Zucker (2003)
  gives 2.9 m/s per epoch and the per-order scatter gives 10.4 m/s; the actual
  repeatability is 702 m/s. Both are internal-consistency measures of a single
  exposure and neither sees the zero point. Any error bar quoted in this project
  should be checked against a repeatability measurement before it is believed.
- **There is real intra-night structure.** 1998-08-25 drifts 1604 m/s across
  6.7 hours and 08-26 627 m/s across 6.9 hours, with all orders moving together.
  This is flexure plus a wavelength solution anchored to an arc at one end of
  the night. It is invisible to prompt 2's S/N-based filter, which is a useful
  reminder that quality has more than one axis.
- **The template is internally consistent to 17 m/s**, which bounds the
  co-addition and interpolation machinery well below anything that matters here.

**Files added.**

- `first_hires_exoplanet/measure_velocities.py`
- `first_hires_exoplanet/data/xcorr_velocities.csv`
- `first_hires_exoplanet/data/xcorr_velocities_per_order.csv`

**Files modified.** `first_hires_exoplanet/quality_filter.py` (glob widened to
`reduce_19*`); `first_hires_exoplanet/build_template.py` (the same
cross-correlation sign error, which there only flipped the sign of an
already-small residual).

**No git command changed state.**

### 2026-09-22 (Prompt 6: assessment -- a clean non-detection, and the measured case for phase 3)

**Task.** Prompt 6 of this document: assess what prompt 5 produced, compare
against the published 3.097-day Keplerian and the modern catalogue in
`hd187123_hires_rv.tsv`, state plainly what the comparison does and does not
demonstrate and what the precision implies for phase 3, with figures from a
script on disk. Full findings are in the Report section above.

**What was done.**

- Wrote `first_hires_exoplanet/figs_phase2.py`: matches our epochs to the
  catalogue by time, fits a circular Keplerian at the forced period to both,
  and writes four figures to `docs/figs/`.
- Figures: `fig_p2_timeseries.png`, `fig_p2_phasefold.png`,
  `fig_p2_compare.png`, `fig_p2_precision.png`.

**Headline results.** All 30 discovery-era catalogue rows match one of our
epochs to 0.1 minutes. Catalogue rms 47.4 m/s against our 723.5. Fitting the
published period: the catalogue gives K = 69.2 +/- 0.3 m/s at 211 sigma,
residual 2.2 m/s; ours gives K = 408 +/- 188 m/s, residual 668 m/s, a 95% upper
limit of 785 m/s. A clean, expected non-detection.

**What this taught us about the repository and the data.**

- **An open question in this document is now closed.** The survey's 32 cell-in
  frames against 30 catalogue velocities were indeed the two 60 s exposures.
  Six of our 36 epochs have no catalogue match and they are exactly the five
  cell-out frames plus 1998-08-12's 60 s cell-in frame; the other short frame,
  1998-07-14's, we never reduced. The one-to-one match of the remaining 30, to
  a tenth of a minute, also independently validates the timing chain from raw
  header through mid-exposure to BJD.
- **A fitted amplitude is not a detection.** Forcing the published period on
  our velocities returns K = 408 +/- 188 m/s at 2.2 sigma, which reads like a
  marginal signal until three things are checked: it is 5.7x the published
  value, it differs from the published value by 1.8 sigma, and forcing the fit
  drops the residual rms only from 724 to 668 m/s. It is the amplitude our noise
  happens to place at that period. The reportable number is the upper limit.
- **A correlation coefficient needs its null stated alongside it.** Ours against
  the catalogue gives r = +0.38, which looks encouraging. If we were measuring
  their signal plus independent noise of our own size, the expected r is +0.07,
  and the standard error with 30 points is 0.19. So +0.38 is about 1.5 standard
  errors from the expectation -- unremarkable. Quoting r alone would have
  invited a claim of partial recovery that the data do not support.
- **The non-detection is the deliverable, and it is a strong one.** Per-order
  scatter within an exposure is 47 m/s, already below the 72 m/s semiamplitude;
  epoch-to-epoch scatter is 724 m/s. The photons are not the limit. That
  localises the entire 30x shortfall in the wavelength zero point, which is
  precisely what an iodine cell removes, and Teklu et al. reach 1.2 m/s on these
  same frames. Phase 2 was framed as the step that produces the argument for
  phase 3; the argument is now a measurement.
- **Equal aspect was the right choice for the comparison figure even though it
  looks bad.** Plotting ours against the catalogue on a square where one axis is
  the other's units collapses 30 points into a vertical stripe. That is not a
  failed figure, it is the result: the catalogue's entire range is +/-142 m/s
  and ours is +/-2 km/s. Rescaling the axes independently would have hidden
  exactly what the figure exists to show.
- **VizieR's declared units are wrong and the repository already knew.**
  `figs.py` records that the RV columns are labelled km/s and are actually m/s.
  Worth having read that note before trusting the column.

**Files added.**

- `first_hires_exoplanet/figs_phase2.py`
- `docs/figs/fig_p2_timeseries.png`, `fig_p2_phasefold.png`,
  `fig_p2_compare.png`, `fig_p2_precision.png`

**No git command changed state.**

### 2026-09-23 (Prompt 7: the atlas settled -- positions transfer, depths are 2.6x off, and an 83 km/s trap)

**Task.** Prompt 7 of this document: obtain
`iodine_atlas/Fischer_Cell_May2022_downsampled3.h5` from `pyodine`, establish
its provenance, compare it against what is known of the HIRES cell, and report
whether it can serve as that cell's atlas, whether the 1998-08-12 pair can test
that empirically, and if not what we would need instead. Full findings are in
the Report section above.

**What was done.**

- Downloaded the atlas (55.8 MB) to `../first-hires-exoplanet-data/atlas/`.
- Wrote `first_hires_exoplanet/atlas_check.py`, which characterises the file,
  measures its resolution from 4653 isolated I2 lines, derives the HIRES cell's
  own transmission from the 1998-08-12 cell-in/cell-out pair, and fits the
  Beer-Lambert exponent between the two. Figure: `docs/figs/fig_p2_atlas.png`.

**Headline results.** The atlas covers 4980-6250 A vacuum on a
wavenumber-uniform FTS grid at R ~ 600,000-680,000, with no metadata of any
kind. Its line positions match the HIRES cell at correlation 0.93. Its depths
do not: the HIRES cell is 2.6x optically thicker, a single Beer-Lambert
exponent. `pyodine`'s linear `iod_depth` cannot express that. And `pyodine`
loads the atlas on the **air** wavelength grid while PypeIt reports vacuum, an
83 km/s offset that destroys the correlation entirely.

**What this taught us about the repository and the data.**

- **The 1998-08-12 pair is worth more than the document credited it with.** It
  was recorded as "a ready-made diagnostic pair for isolating the cell's
  transmission empirically", which turned out to be exactly right, and it is
  the only way this question could have been settled from our own data. It also
  justifies having reduced that night in prompt 3, where it was nearly lost to
  the `TARGNAME = 'H187123'` typo.
- **The air/vacuum trap has now caught this project twice in different
  places.** Phase 1 found a spurious +70 km/s comparing vacuum observations
  against air rest wavelengths; here `pyodine` loads `wavelength_air` while our
  spectra are vacuum, an 83 km/s offset that drops the correlation from +0.90
  to +0.08. Any external wavelength reference entering this project should have
  its frame checked before it is used, as a standing rule.
- **A free parameter is not the same as the right functional form.** `pyodine`
  has an `iod_depth` parameter, which on first reading settles the depth
  mismatch. It applies it linearly, T -> 1 + d(T-1), where the physics is
  T -> T**alpha. For the d = 2.6 we need, 5.8% of the atlas goes to negative
  transmission. Reading the parameter list was not enough; the line that uses
  it had to be read too.
- **Choose the statistic that matches the physics.** A mean-depth ratio gave
  2.50 with no way to tell whether one number could describe the difference.
  Fitting the Beer-Lambert exponent gave 2.79, and restricting to the
  transmission band where continuum placement and saturation do not bite gave
  2.59 with the spread collapsing from 2.26-3.53 to 2.28-2.86. That collapse is
  the actual result: one exponent is enough, so the cells differ in column
  density and not in temperature.
- **A trend in a fitted parameter is often the method, not the object.** Alpha
  ran 3.5 to 2.3 across the band and was largest where the absorption was
  weakest, which is where a running upper-percentile continuum works best. The
  physical reading -- a temperature difference reshuffling line strengths --
  would have been wrong.
- **Normalise both sides the same way.** The first version put the measured
  spectrum through a running percentile filter and the atlas through a single
  scalar. Any continuum difference then lands in the depth comparison.
- **Noise biases an upper-envelope continuum.** Two 60 s frames give a 1.7%
  ratio noise, and a 95th-percentile continuum rides up on it, deepening every
  line. Injecting matched noise into the atlas and re-running the identical
  normalisation returned alpha = 1.07, which bounds the artefact at 7%. Worth
  doing whenever a depth or equivalent width is measured off a noisy spectrum.
- **`h5py` is not in `pypeit14`.** `pyodine` requires it, so phase 3 will need
  it installed there. This module sidesteps the issue by reading spec1d files
  with plain `astropy.io.fits` rather than `pypeit.specobjs`, which also
  demonstrates that the spec1d layout is simple enough not to need PypeIt to
  read it -- useful for the prompt-8 adapter.

**Files added.**

- `first_hires_exoplanet/atlas_check.py`
- `docs/figs/fig_p2_atlas.png`
- `../first-hires-exoplanet-data/atlas/Fischer_Cell_May2022_downsampled3.h5`
  (data tree, not the repo)

**No git command changed state.**
