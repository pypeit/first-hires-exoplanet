# Data, Phase 1: reduce the 1998 HIRES spectra of HD 187123

## Goals

Take the original Keck/HIRES observations behind
[Butler, Marcy, Vogt & Apps (1998)](https://doi.org/10.1086/316287) and reduce
them from scratch with `PypeIt`. This is item 1 of the three-step scope agreed in
`context_prompts.md`:

1. **Reduce and inspect** — clean, wavelength-calibrated, extracted 1D spectra of
   HD 187123 with the iodine lines visible, and evidence that the reduction is
   sound. **This document.**
2. Relative radial velocities by cross-correlation.
3. Full iodine forward-model velocities, aiming at the published 72 m/s.

Phase 1 succeeds when we can look at an extracted spectrum from a 1998 night and
say, with reasons, that it is good. It does **not** attempt a velocity.

## What we already know

Established while writing `docs/public_HD187123b.md`. Check anything that looks
stale before relying on it.

- **The target.** HD 187123, G2V, Teff 5790 K, [M/H] +0.09, 46 pc, in Cygnus.
  Planet b: P = 3.0965828 d, K = 72 m/s, M sin i = 0.52 M_Jup, a = 0.042 AU.
  There is also a long-period planet c (P ~ 10.4 yr), which matters for velocity
  work later but not for phase 1.
- **The discovery-era observations.** From the modern catalogue there are **30
  epochs between 1997 Dec 23 and 1998 Sep 18**, clustered into runs: one epoch
  in 1997 Dec, one in 1998 Jun, ten across five nights 1998 Jul 15-19, four on
  Aug 17-18, seven on Aug 25-26, and seven over Sep 12-18. The July run is the
  "five-night string" the literature credits with revealing the 3-day period,
  and is the obvious first target. Epoch list:
  `first_hires_exoplanet/data/hd187123_hires_rv.tsv`.
- **The detector.** 1998 data predates the August 2004 upgrade, so it is the
  original single Tektronix CCD with 24-micron pixels, *not* the three-CCD
  mosaic.
- **PypeIt now supports that detector — but not in our checkout.** As of
  September 2026 `origin/develop` has a `KeckHIRESOrigSpectrograph`, spectrograph
  name **`keck_hires_orig`**, `ndet = 1`, described as "Pre detector upgrade
  (~August 2004) ... the original Tektronix CCD". It was added by Ryan Cooke in
  commits `a7ce0fab6` and `aed8b71b6` (1-2 Sep 2026), followed by
  `4e6bb723f` (support for the very old RED97 cross-disperser) and `b76bccacf`
  (a hard-coded Tektronix bad-pixel mask). The docs at
  `doc/spectrographs/keck_hires.rst` on develop describe both eras.

  **None of these commits is in our local `PypeIt` checkout**, which sits on
  branch `hamspec`, and the installed `pypeit` in the `pypeit14` environment is
  `2.0.2.dev891+gf8a757720` from that checkout. So the very first obstacle is an
  environment problem, not a data problem. Ryan Cooke is a collaborator (he is a
  co-author on `the-holy-grail`) and is the person to ask if the original-CCD
  support behaves oddly.
- **HIRES is still `supported = False`** in PypeIt generally. We are early
  adopters; expect rough edges and be willing to report them upstream.
- **The development suite has no HIRES raw data.** `RAW_DATA/` contains no HIRES
  directory, so there is no worked example to copy and no regression test to
  lean on. Everything comes from the Keck Observatory Archive.

### Risks worth probing early

- **Calibrations may be thin.** The planet search used the iodine cell as its
  wavelength fiducial, so those nights may carry few or no ThAr arcs, and
  whatever flats exist were taken for a programme with different priorities than
  ours. If PypeIt needs arcs that were never taken, that shapes the whole phase.
- **The iodine cell is in the science beam.** Every science frame of HD 187123
  is imprinted with thousands of I2 lines between roughly 500 and 620 nm. That
  is the point of the observation, but it will complicate anything that assumes
  a clean stellar spectrum — continuum fitting, order tracing off the science
  frames, automated wavelength cross-correlation.
- **`pypeit14` has no `setuptools`.** Noted during start-up; may bite on an
  editable install.

## Code

See the guidelines in the prompt docs in
`Projects/PypeIt/PypeIt-development-suite/claude_prompts` for how to code in
Python.

If you need to run Python, use the `pypeit14` environment — unless prompt 1
concludes we need a separate environment, in which case record the new name
here.

The development suite carries skills that apply directly to this work, under
`PypeIt-development-suite/.claude/skills/`: `diagnose-reduction` for triaging a
failed or suspect `run_pypeit` run, `wavelength-calibration` for arcs, line
lists and templates, `run-dev-suite` for the test harness, and
`add-devsuite-setup` if we end up contributing a HIRES setup back upstream.
Use them.

Record the exact PypeIt commit used for any reduction whose results we quote.
The installed version is a development checkout, so "PypeIt 2.0.2" is not a
reproducible statement on its own.

## Prompts

1. Read this file. We need a PypeIt that provides `keck_hires_orig` without
   disturbing the `hamspec` branch currently checked out in
   `Projects/PypeIt/PypeIt`. Explore the options — a separate clone, a git
   worktree, a second conda environment, or simply switching branches — and
   write a recommendation into the Report section below. **Do no installing or
   branch-switching yet.** Use Opus 5. Log your work.

2. Read this file. Carry out the recommendation from prompt 1 and verify that
   `keck_hires_orig` is importable and appears in `pypeit_show_spectrographs` (or
   the equivalent). Record the commit hash. Use Opus 5. Log your work.

3. Read this file. Survey the Keck Observatory Archive
   (`koa.ipac.caltech.edu`) for HIRES observations of HD 187123 between 1997
   December and 1998 October, and for the calibration frames taken on those same
   nights. Report what exists: frame types, how many of each, binning,
   cross-disperser and echelle angles, deckers, exposure times, and whether
   everything is public. Pay particular attention to whether ThAr arcs exist for
   the 1998 July 15-19 run. Write this into the Report section. **Do not
   download anything yet.** Use Opus 5. Log your work.

4. Read this file. Read `keck_hires_orig` in the PypeIt source and report what it
   expects: the header cards used by `init_meta`, the configuration keys, the
   detector parameters, the frame-typing rules, the required calibrations, and
   the default `PypeItPar`. Compare that against what prompt 3 found in the
   archive, and list every gap. **Do no editing.** Use Opus 5. Log your work.

5. Read this file. Download one night from the 1998 July 15-19 run, with its
   calibrations, into a directory outside the repository (the raw frames must
   not be committed; `.gitignore` excludes a top-level `data/`). Report what
   arrived and confirm the headers match what prompt 4 expected. Use Opus 5. Log
   your work.

6. Read this file. Run `pypeit_setup` on that night. Report the configurations
   it identifies, the frame types it assigns, and any frame it fails to type or
   assigns wrongly. Do not hand-edit the `.pypeit` file yet — first report what
   the automatic pass produced. Use Opus 5. Log your work.

7. Read this file. Reduce that night with `run_pypeit`. If it fails, use the
   `diagnose-reduction` skill to work out why, fix what can be fixed in the
   input file or parameters, and report anything that looks like a genuine
   PypeIt bug rather than a configuration mistake. Use Opus 5. Log your work.

8. Read this file. Inspect the reduction. Check the order tracing, the flat
   field, and the wavelength solution against the QA products, and look at the
   extracted 1D spectrum: are the iodine lines present between roughly 500 and
   620 nm, are the stellar lines where they should be for a G2V star, and what
   signal-to-noise did we get? Write a figure-generating script to disk, as with
   the public document. Use Opus 5. Log your work.

9. Read this file. Reduce the remaining nights of the 1998 July run the same
   way, and report whether the reduction is stable from night to night. Use
   Opus 5. Log your work.

10. Read this file. Write a short assessment of phase 1: whether the reduction is
    sound and how we know, what PypeIt handled well, what needed intervention,
    anything worth reporting upstream, and what phase 2 would require. Use
    Opus 5. Log your work.

## Q&A

## Report

## Logs
