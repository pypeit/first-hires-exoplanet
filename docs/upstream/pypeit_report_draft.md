# DRAFT — not sent

**To:** Ryan Cooke
**From:** J. Xavier Prochaska
**Subject:** PypeIt: fifteen items from reducing the 1998 HIRES discovery data for HD 187123 b

---

Hi Ryan,

We have spent the past few weeks reducing the original-CCD Keck/HIRES frames
behind the 1998 discovery of HD 187123 b with PypeIt: 36 frames, 20 nights,
1997-12 to 1998-09. The iodine-cell forward model on top of that reduction
recovers the planet at 10σ: K = 72.4 ± 7.1 m/s, 26 m/s per epoch. So the
reduction is sound end to end. Along the way we collected a list of PypeIt
items. Most are small, a few matter to anyone measuring velocities, and two
are data-integrity issues. The fixes for items 1–5 and part of 9 are already
on `orig-hires-fixes` (`017bece06`, `keck_hires.py` only), which is
`develop` at `f3a1f1d27` plus that one commit.

Two scripts reproduce everything that can be reproduced without our data
directory:

- `first_hires_exoplanet/verify_orig_fixes.py` runs `keck_hires_orig` over
  real 1998 frames and prints everything items 1–6 touch (metadata, raw
  sections, the bad-pixel mask, `check_spectrograph`, frame typing), plus a
  regression check on the post-2004 `keck_hires` class. Run it before and
  after the branch and diff the output.
- `first_hires_exoplanet/refframe_audit.py` recomputes PypeIt's stored
  `VEL_CORR` for all 36 epochs against `astropy` and breaks the difference
  into its terms (items 10–11).

Both live at <repository link>. The raw frames are public in KOA; the
KOAIDs are in the scripts.

## A. Bugs in `keck_hires_orig`, fixed on `orig-hires-fixes`

1. **`get_rawimage` offsets single-amplifier frames by 21 columns.** The
   data section starts at `prepix * 2`, which is right only for
   `NUMAMPS = 2`. On a 1998 single-amp frame the data landed 21 columns late
   and the overscan was truncated to 213 of 234 columns. Fix:
   `prepix * namp` (and the same in the overscan offset).
2. **The MAKEE bad-pixel mask lands ~20 px from the defects.** MAKEE column
   numbers count raw columns including the prescan, but the mask was applied
   at trimmed column `C − 1`. Fix: `C − PREPIX`, read from the example file,
   with a warning if it cannot be.
3. **The binning string is inverted for the original CCD.** `'1,2'` was
   reported as `'2,1'`, which made `order_platescale` 0.216 instead of
   0.432″ per binned pixel. On this detector `BINNING` is
   `'<spectral>,<spatial>'`. Fix: a `compound_meta` override in
   `KECKHIRESOrigSpectrograph`.
4. **`check_spectrograph` builds a `PypeItError` and never raises it.** Fix:
   `raise`.

## B. Defaults and behaviour

5. **`scienceframe.exprng = [601, None]` excludes every bright-target
   original-CCD exposure** (the planet-search frames are 200–500 s). We set
   `[1, None]` for the Orig class only; science and standards are still
   separated by `vet_assigned_ftypes`. This is a policy change you may want
   to review.
6. **`find_trim_edge = [3, 3]` is untenable at 2× spatial binning,** where
   the B1 orders are ~10 px wide: it masks six of ten columns. A
   binning-aware default would help. We override it in the reduction; it is
   not in the branch.
7. **`run_pypeit` exits 0 having written no `spec1d`.** This was the most
   expensive behaviour we met: a silent failure that reads as success. A
   non-zero exit, or at least a closing warning, when no science output was
   written would have saved hours.
8. **`run_pypeit -o` appends to an existing `spec1d` rather than replacing
   it.** Iterating on parameters silently superimposes reductions: we found
   74 extensions for 37 orders, with the *bad* pass listed first. This is a
   data-integrity issue.
9. **Docstring and header questions.** `config_independent_frames` describes
   a DATE-OBS rule that does not exist (corrected on the branch). The headers
   give `CCDGN01 = 4.8` e⁻/ADU and `CCDRN01 = 6.0` e⁻ against the hard-coded
   1.9 and 2.8. We did not change these, but item 12 now measures what the
   difference does.

## C. `pypeit/core/wave.py` — affects every user, not just HIRES

10. **`geomotion_velocity` has a sign error on the solar term**
    (`wave.py:128`, `velocity += sv`, still present on `develop` at
    `f3a1f1d27`). For the heliocentric frame the observer's velocity
    relative to the Sun is `ev + ov − sv`.
    - Against `astropy`'s `radial_velocity_correction` over our 36 epochs,
      the `−` form agrees to 0.00 m/s and PypeIt's `+` form is **+9.1 to
      +14.2 m/s** off.
    - The error drifts by **5.2 m/s** over the nine months, so it does not
      cancel in relative velocities.
    - `refframe = heliocentric` is the **default**, so every user who does
      not override it gets this.
11. **`geomotion_correct` omits the relativistic terms.** Barycentric
    velocities come out **−4.6 m/s** from `astropy`'s: the solar
    gravitational redshift (+3.0 m/s) and the observer's time dilation
    (+1.5 m/s). It is nearly constant (0.2 m/s peak-to-peak here), so it is
    harmless for relative velocities but wrong for an absolute one.

    `refframe_audit.py` sections 5–6 reproduce both items. Separately, the
    header MJD of these frames marks the exposure *start*, to within ~2 s
    (section 7, against the hour angle). PypeIt evaluates the correction at
    that time. At mid-exposure (`MJD + EXPTIME/2`) it differs by
    +3.0 to +4.4 m/s, varying from epoch to epoch. Whether PypeIt should use
    the midpoint is your call.

## D. New since phase 2 (found by the iodine forward model)

12. **The inverse variance is mis-scaled, and differently for optimal and
    boxcar extraction.**
    - Three consecutive exposures of HD 187123 share one arc and one
      wavelength solution. Differencing them pixel by pixel, the true
      scatter is **0.426×** the σ propagated in `OPT_COUNTS_IVAR` (clean
      pair, all 14 orders consistent).
    - The header gain (4.8 against 1.9 e⁻/ADU) explains a factor of 1.59
      of that; a factor of ~1.5 remains unexplained. It is not
      `noise_floor` (0.0 here), nor the local-sky `adderr`.
    - On a boxcar extraction of a pair of B-star frames, `BOX_COUNTS_IVAR`
      goes the *other* way: the true scatter is **1.61×** the propagated σ.
    - An optimal and a boxcar extraction of the same kind of data
      disagreeing by 3.8× in their noise model looks worth a look. Every χ²
      computed from PypeIt's ivar on this data is off by these factors.
13. **Masked pixels get a wavelength of exactly 0 in the `spec1d`.** In one
    epoch, 950 pixels have `OPT_WAVE == 0.0`, every one of them where
    `OPT_MASK` is False. A pixel's wavelength is known whether or not its
    flux is usable. Code that divides by wavelength, or takes the first and
    last wavelength of an order, breaks: `pyodine` did, at its first step.
    Suggestion: write the wavelength solution everywhere and let the mask
    say what is usable, or write NaN rather than 0.
14. **Boxcar flux in partially masked apertures is biased low and not
    flagged.** In a boxcar extraction (`force_center_obj`, a slit-filling
    B star), 1,730 iodine-order pixels (6%) have `BOX_NPIX` below 75% of the
    full aperture. Their `BOX_COUNTS` sit at 0.51 of the local continuum,
    scaling with the fraction summed (r = 0.84), while `BOX_MASK` is True.
    Either renormalising by the profile fraction or flagging these pixels
    would stop them reading as absorption lines.
15. **Orders that fail calibration vanish from the `spec1d` without a
    warning.**
    - On 1998-09-13 the night's ThAr arc cross-correlates against the
      archive at a median cc of 0.65 (below 0.8 in every order), against
      0.91 on every other night.
    - Orders 92 and 71 failed to reidentify and order 86 fitted at
      0.477 px. All three were flagged
      `BADWVCALIB | BADTILTCALIB | BADFLATCALIB` and left out: 37 orders
      traced, 34 extracted.
    - The only trace is an INFO line, "Skipping bad slit".
    - Downstream code that assumes the full set of orders is silently
      wrong: `pyodine` fitted order 71's template against order 72.
    - A closing WARNING listing the orders dropped, and why, would make this
      visible.

## E. For information: what the iodine says about the ThAr wavelengths

Not a bug report. The forward model fits its own wavelength solution to the
iodine lines on every chunk, so it measures PypeIt's wavelength zero point
on the science photons:
- **The zero-point error is −0.9 to +2.6 km/s across the 36 epochs.** The
  large end is 1998-07-19, which borrows the previous night's arc.
- **It drifts by up to 256 m/s per hour within a night** (1998-08-25: 1.7
  km/s over 6.7 h).
- **Its size grows with the time between the science frame and the arc**
  (Spearman +0.61).

That is flexure against an end-of-night arc: expected, but now measured.
Anyone using PypeIt's HIRES wavelengths for velocities better than ~1 km/s
needs a same-time calibration.

Two things worked out of the box and deserve saying:
- **`force_center_obj`** (from the Hamspec work) rescued the B stars, V ≈ 3–4
  stars that overfill the 3.5″ slit. Without it, peak-finding found HR 8634
  in 1 order of 37. With it, all 11 frames extracted.
- **Cross-night arcs** reduce as well as same-night ones, as far as the fit
  RMS can tell. Item E is the caveat: the RMS is not the zero point.

## F. Three usability notes

- **"No frames of type=trace provided" names the symptom, not the cause.**
  In all three cases we hit, the cause was flats PypeIt had (rightly)
  declined to frametype. Listing the rejected candidates would have pointed
  straight at them.
- **The `RED97` rule deserves a reference in the code.** `keck_hires.py:245`
  overrides a header value on the MAKEE DRP's authority, with a comment but
  no citation. It is a hard configuration key and it partitions the archive
  at 1997-12-31. We measured December's orders one binned pixel from July's
  at correlation 0.99. That does not contradict a cross-disperser swap, but
  it is worth someone knowing.
- **`refframe` is load-bearing and invisible.** Nothing in the reduction
  inputs records that a ~4 km/s shift was applied; the only trace is
  `VEL_TYPE` and `VEL_CORR` in the spec1d extensions.

Happy to open issues or PRs for any of these. The branch could become a PR
for items 1–5 as it stands. Let me know how you would like them.

Best,
X.
