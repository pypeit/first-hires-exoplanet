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
- **PypeIt now supports that detector, and our checkout has it.** As of
  September 2026 `origin/develop` has a `KeckHIRESOrigSpectrograph`, spectrograph
  name **`keck_hires_orig`**, `ndet = 1`, described as "Pre detector upgrade
  (~August 2004) ... the original Tektronix CCD". It was added by Ryan Cooke in
  commits `a7ce0fab6` and `aed8b71b6` (1-2 Sep 2026), followed by
  `4e6bb723f` (support for the very old RED97 cross-disperser) and `b76bccacf`
  (a hard-coded Tektronix bad-pixel mask). The docs at
  `doc/spectrographs/keck_hires.rst` on develop describe both eras.

  **Resolved 2026-09-19 (prompt 1).** Our `PypeIt` checkout is now on `develop`
  at `f3a1f1d27` (merge of PR #2197, `keck_hires_tektronix`), clean and in sync
  with `origin/develop`, and all four commits above are ancestors of HEAD. The
  install is editable and maps to this checkout, so `keck_hires_orig` loads and
  runs. The only residue is a stale `pypeit.__version__`
  (`2.0.2.dev891+gf8a757720`, 325 commits behind HEAD) — cosmetic for running,
  but it is the string PypeIt stamps into outputs, so it must be refreshed
  before we quote any reduction. Ryan Cooke is a collaborator (he is a co-author
  on `the-holy-grail`) and is the person to ask if the original-CCD support
  behaves oddly.
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
- **`pypeit14` has no `setuptools`.** Noted during start-up. Checked in prompt
  1: it does *not* break the existing editable install (the PEP 660 finder in
  site-packages is self-contained and no PypeIt module imports `pkg_resources`),
  and a plain `pip install -e .` still works because build isolation fetches
  setuptools into a temporary build environment. It would only bite under
  `--no-build-isolation`, or offline.

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

1. Read this file. We need to check that the `PypeIt` software is ready to work with 
   the original HIRES detector.  I have switched us onto the `develop` branch of `PypeIt`.  Let's check that it is ready to work with the original HIRES detector.  If we find that we need to do some additional development, I will switch us onto a new branch and work on that.
   Use Fable if you can. Log your work.

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

### 2026-09-19 (Prompt 1: confirmed PypeIt `develop` is ready for the original HIRES detector)

**Task.** Prompt 1 of this document: check whether the `PypeIt` software, now on
the `develop` branch, is ready to reduce original (pre-2004, single Tektronix
CCD) HIRES data, and say whether any development is needed. Investigation was
delegated to a Fable subagent as the prompt asked; every load-bearing claim was
then re-verified directly in this session.

**Verdict: ready as-is. No development and no new branch are needed for the
reduction itself.**

Evidence:

- **Branch and commits.** `/Users/xavier/Projects/PypeIt/PypeIt` is on `develop`
  at `f3a1f1d274b15ee1358f167819d77f1948fce1bd` ("Merge pull request #2197 from
  pypeit/keck_hires_tektronix", 2026-09-14). Working tree is clean and level
  with `origin/develop` (no ahead/behind). All four commits this document names
  are ancestors of HEAD, confirmed with `git merge-base --is-ancestor`:
  `a7ce0fab6`, `aed8b71b6` (original-CCD functionality), `4e6bb723f` (RED97
  cross-disperser), `b76bccacf` (hard-coded Tektronix BPM). The earlier note in
  "What we already know" — that none of these were in our checkout, which sat on
  `hamspec` — is now stale and has been updated in place.
- **Class.** `KeckHIRESOrigSpectrograph` is at
  `pypeit/spectrographs/keck_hires.py:928`, with `name = 'keck_hires_orig'`
  (:931 area), `ndet = 1`, and the "Pre detector upgrade (~August 2004)" comment.
  It subclasses `KECKHIRESBaseSpectrograph`, which carries `camera = 'HIRES'`,
  `header_name = 'HIRES'`, `pypeline = 'Echelle'`, `ech_fixed_format = False`.
- **Registered and loads.** At runtime in `pypeit14`:
  `available_spectrographs` contains both `keck_hires` and `keck_hires_orig`;
  `load_spectrograph('keck_hires_orig')` returns `name/ndet/supported =
  keck_hires_orig 1 False`; `get_detector_par(1)` builds without a traceback
  (platescale 0.216, ronoise [2.8], saturation 65535); `default_pypeit_par()`
  sets `rdx.detnum = [1]` and `calibrations.wavelengths.ech_separate_2d = False`;
  `configuration_keys()` = `['dispname', 'decker', 'filter1', 'echangle',
  'xdangle', 'binning']`.
- **Docs.** `doc/spectrographs/keck_hires.rst` documents both eras, including a
  dedicated "Original Tektronix Detector" section covering the hard-coded BPM and
  the three cross-dispersers.

**What I learned about the repository / environment.**

1. *The install is editable and already points at the checkout.* `pip show
   pypeit` reports `Editable project location:
   /Users/xavier/Projects/PypeIt/PypeIt`, and `pypeit.__file__` resolves inside
   it. So switching branches changes the code that runs, with no reinstall
   needed. **But `pypeit.__version__` is stale**: it still reads
   `2.0.2.dev891+gf8a757720`, a commit 325 behind HEAD. The version comes from
   `setuptools_scm` writing the gitignored `pypeit/pkg/version.py` at install
   time; it is only regenerated by re-running the install. This is cosmetic for
   running, but it is the string PypeIt stamps into output headers and logs — so
   under this document's own rule ("Record the exact PypeIt commit used for any
   reduction whose results we quote") it must be refreshed *before* we produce
   any reduction we intend to publish. **Recommendation carried to prompt 2:**
   `cd /Users/xavier/Projects/PypeIt/PypeIt && conda run -n pypeit14 pip install
   -e . --no-deps` (needs network, because build isolation fetches setuptools).
2. *Two claims in this document about tooling are wrong and would have cost time
   later.* There is **no `pypeit_show_spectrographs` script** — nothing matching
   in `pypeit/scripts/`. The equivalent check is the spectrograph choice list in
   `pypeit_setup -h` or `pypeit_obslog -h`, or the runtime import. And the
   import path is `from pypeit.spectrographs.util import
   available_spectrographs`, *not* `from pypeit.spectrographs import ...`, which
   raises `ImportError`. Prompt 2 should use these forms.
3. *The `setuptools` worry is mostly a non-issue.* `pypeit14` genuinely has no
   `setuptools`, but the existing editable install does not depend on it at
   runtime, and a normal `pip install -e .` works because build isolation
   supplies it. Only `--no-build-isolation` (or an offline machine) would fail.
   Noted in the risks list.
4. *`supported = False` is a project-level flag, not a blocker.* The base class
   keeps it False behind a TODO list; `keck_hires_orig` does not override it. It
   does not gate any code path we need — it is PypeIt's statement about how
   well-exercised HIRES is, which is exactly the early-adopter posture this
   document already assumes.
5. *RED97 handling is date-driven, which matters for our 1997 Dec epoch.* In
   `KECKHIRESBaseSpectrograph.compound_meta`, `dispname` is forced to `"RED97"`
   when `mjd < 50814.0` (midnight 31 Dec 1997), and for RED97 data the frame
   `idname` is derived from `IMAGETYP`, with `Bias` vs `Dark` split on
   `ELAPTIME < 0.001`. Our 1998 July target run is *after* that cutoff so it
   takes the normal RED/UV path from the header — but the single 1997 Dec 23
   epoch would be classified as RED97, i.e. a different configuration. Worth
   remembering in prompt 3 when we survey the archive, and in prompt 6 when
   `pypeit_setup` groups configurations.
6. *Minor upstream nit, not worth a report on its own:*
   `doc/spectrographs/keck_hires.rst:74` has a missing space before ``RED97``,
   an RST cosmetic. Bank it in case we open an issue for something real.

**Nothing was edited in the PypeIt source, and no git state was changed.** The
only edits were to this document: the stale checkout/branch note and the
`setuptools` risk bullet.
