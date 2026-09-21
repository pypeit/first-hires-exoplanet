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
  runs. The stale `pypeit.__version__` noted in prompt 1 has since been
  refreshed (prompt 2): it now reads **`2.0.2.dev1216+gf3a1f1d27`**, which
  embeds HEAD. **That is the version/commit string to quote for any reduction we
  publish.** Ryan Cooke is a collaborator (he is a co-author
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
   the equivalent). Record the commit hash. Use Fable if you can. Log your work.

3. Read this file. Survey the Keck Observatory Archive
   (`koa.ipac.caltech.edu`) for HIRES observations of HD 187123 between 1997
   December and 1998 October, and for the calibration frames taken on those same
   nights. Report what exists: frame types, how many of each, binning,
   cross-disperser and echelle angles, deckers, exposure times, and whether
   everything is public. Pay particular attention to whether ThAr arcs exist for
   the 1998 July 15-19 run. Write this into the Report section. **Do not
   download anything yet.** Use Fable if you can. Log your work.

4. Read this file. Read `keck_hires_orig` in the PypeIt source and report what it
   expects: the header cards used by `init_meta`, the configuration keys, the
   detector parameters, the frame-typing rules, the required calibrations, and
   the default `PypeItPar`. Compare that against what prompt 3 found in the
   archive, and list every gap. **Do no editing.** Use Fable if you can. Log your work.

5. Read this file. Download one night from the 1998 July 15-19 run, with its
   calibrations, into a directory outside the repository (the raw frames must
   not be committed; `.gitignore` excludes a top-level `data/`). Report what
   arrived and confirm the headers match what prompt 4 expected. Use Fable if you can. Log
   your work.

6. I have put `PypeIt` on a new branch called `orig-hires-fixes`.  Please fix the 
   issues you have identified in prompt 4 and 5.  I will then push the branch to the remote repository. Use Fable if you can. Log your work.

7. Read this file. Run `pypeit_setup` on that night. Report the configurations
   it identifies, the frame types it assigns, and any frame it fails to type or
   assigns wrongly. Do not hand-edit the `.pypeit` file yet — first report what
   the automatic pass produced. Use Fable if you can. Log your work.

8. Read this file. Reduce that night with `run_pypeit`. If it fails, use the
   `diagnose-reduction` skill to work out why, fix what can be fixed in the
   input file or parameters, and report anything that looks like a genuine
   PypeIt bug rather than a configuration mistake. Use Fable if you can. Log your work.

9. Read this file. Inspect the reduction. Check the order tracing, the flat
   field, and the wavelength solution against the QA products, and look at the
   extracted 1D spectrum: are the iodine lines present between roughly 500 and
   620 nm, are the stellar lines where they should be for a G2V star, and what
   signal-to-noise did we get? Write a figure-generating script to disk, as with
   the public document. Use Fable if you can. Log your work.

10. Read this file. Reduce the remaining nights of the 1998 July run the same
   way, and report whether the reduction is stable from night to night. Use
   Fable if you can. Log your work.

11. Read this file. Write a short assessment of phase 1: whether the reduction is
    sound and how we know, what PypeIt handled well, what needed intervention,
    anything worth reporting upstream, and what phase 2 would require. Use
    Fable if you can. Log your work.

## Q&A

## Report

### Prompt 3 (2026-09-20): KOA survey of HIRES frames on HD 187123, 1997 Dec - 1998 Oct

Metadata only; **no FITS files were downloaded**. Script:
`first_hires_exoplanet/koa_survey.py`. Raw results:
`first_hires_exoplanet/data/koa_hd187123_science.csv` (37 rows) and
`koa_hd187123_nights_all.csv` (3124 rows -- every frame on the 21 nights).
The script re-runs offline from those CSVs; `--refresh` re-queries KOA.

**Query.** POST to `https://koa.ipac.caltech.edu/TAP/sync` with
`REQUEST=doQuery, LANG=ADQL, FORMAT=csv, MAXREC=100000`:

    SELECT <54 columns> FROM koa_hires
     WHERE mjd >= 50783.0 AND mjd < 51118.0
       AND (upper(targname) LIKE '%187123%' OR upper(object) LIKE '%187123%')
     ORDER BY mjd

then one query per UT night, `... WHERE koaid LIKE 'HI.YYYYMMDD.%'`, with the
night list derived from the science `koaid` prefixes. Gotchas: the table is
`koa_hires` (`koa.koa_hires` returns 403); filter on `mjd` (double), because
`date_obs` is a char column and filtering it raises `ORA-01861`; it is `xdangl`,
not `xdisangl`; `filename` is not a column (use `filehand` / `ofname`).

**Public status: everything below is public.** KOA silently appends a
proprietary-period clause to every query, so it only ever returns rows whose
`propint` has elapsed. Our frames carry `propint` 115-202 months across semids
`1997B_U07H` ... `1998B_U05H`; the longest (202 months from Dec 1997) expired in
2014.

**Census.** 37 `object` frames of HD 187123 on 21 UT nights. Exactly 30 match the
30 catalogue epochs in `hd187123_hires_rv.tsv` (exposure midpoints agree to
< 0.05 min): 1 + 1 + 10 + 4 + 7 + 7 -- **the document's epoch count is
confirmed against the archive.** The 7 extra frames are 5 iodine-free
(`iodin = F`) exposures -- 1997-12-24 (500 s), 1998-08-12 (60 s), and three
consecutive 500 s frames on 1998-08-26, which look like a deliberate template
set and matter for phase 3 -- plus 2 short Latham frames (1998-07-14, 1998-08-12).
All 37 frames: `xdispers = RED`, `binning = '1,2'`, `fil1name = kv370`,
`fil2name = clear`, 1 amplifier, `ccdgain = F`; decker B1 except 1998-09-17 (B2).

#### The 1998 July 15-19 run

Program `N14H`, PI G. Marcy, "True Jupiter Analogs: A Survey within 50 parsecs",
semid `1998A_N14H`, `propint` 117. **Every night is shared with program `N01H`
(T. Bida, Mercury), whose daytime frames (02:55-04:30 UT) are in an unrelated
setup** -- decker B2, binning 1,1, filter bg24a, echangl +0.266, xdangl +1.32,
2 amplifiers. Those are excluded below and should be excluded from any download.

Counts verified directly against KOA (558 frames total across the five nights):

| UT night | HD 187123 (UT, s) | echangl | xdangl | ThAr arcs, B1 | clean flats B2 2 s | I2 flats B1 3 s | bias |
|---|---|---|---|---|---|---|---|
| 07-15 | 1 (10:34, 260) | 0.0019450 | -0.5405 | **2** (07:29, 15:07) | 16 | 2 | 0 |
| 07-16 | 1 (14:43, 400) | 0.0009302 | -0.5455 | **2** (07:39, 15:05) | 16 | 2 | 0 |
| 07-17 | 2 (08:56, 300; 13:32, 400) | 0.0009302 | -0.5445..-0.5475 | **1** (15:06) | 16 | 1 | 0 |
| 07-18 | 3 (08:07, 250; 10:44, 215; 13:44, 400) | 0.0008880 | -0.5495..-0.5525 | **1** (15:09) | 16 | 1 | 0 |
| 07-19 | 3 (06:54, 427; 09:55, 300; 13:53, 427) | -0.0000423 | -0.5495..-0.5535 | **0** | **0** | **0** | 0 |

**ThAr arcs exist, and this is the key result of the survey.** They are
`koaimtyp = arclamp`, `lampname = ThAr1`, `lampcat1 = T` (hollow cathode),
`lmirrin = T`, lamp filter ng3, 10 s, decker B1, binning 1,2, filter kv370 --
i.e. *the science configuration*. Angle offsets from the same night's science
frames are 0.0000 in echangl and <= 0.004 deg in xdangl, against PypeIt
`np.isclose` tolerances of atol 0.01 and 0.1 respectively: a margin of more than
20x. Caveat: each arc is 10 s and carries 365-1926 saturated pixels (`npixsat`),
so the strongest ThAr lines are saturated.

**July 19 has no calibrations of any kind** -- no arc, no flat, no bias. The
Marcy program's last frame that night is the 13:53 UT science exposure. The
nearest usable arc is `HI.19980718.54587` (15:09 UT on Jul 18): within tolerance
on both angles, but taken ~15.8 h before the first July 19 science frame.

**Flats come in two families, and neither is ideal.** The 16 clean quartz flats
per night (`lampqtz2 = T`, lamp filter bg14, 2 s, `iodin = F`, unsaturated) use
decker **B2** (7.0"), not the science **B1** (3.5"). The only B1 flats are the
1-2 per night taken with the **iodine cell in the beam** (`iodin = T`, 3 s) --
the planet-search "iodine flats", which carry the I2 forest. Since `decker` is
one of PypeIt's `configuration_keys`, `pypeit_setup` will put the clean flats in
a *different* configuration from the science frames.

**Hatch open on nearly all flats.** 62 of the 70 July flats have `hatopen = T`;
only the 8 B2 flats at 15:08-15:20 UT on 07-16 have `hatopen = F`. PypeIt's
HIRES `idname` rule types a frame as an internal flat only when the hatch is
closed, so most flats may come back untyped from the automatic pass. Prompt 4
should confirm against the code and prompt 6 against real headers.

**No biases or darks on any of the five nights.** The nearest same-configuration
ones are on 1998-07-14 (10 biases at 0 s, plus 39 darks of 60-180 s) and
1998-08-12. PypeIt keys bias/dark only on `dispname` and `binning`
(`config_independent_frames`), so the July 14 biases are the natural candidates.

**Angles will not split the run.** echangl is constant to 8 decimals within each
night and spans 0.0020 deg across the five nights; xdangl drifts up to 0.004 deg
within a night and 0.013 deg across the run. Both are comfortably inside PypeIt's
tolerances, so the five nights form one configuration on angle grounds. What
*will* split configurations is `decker` (B1 science vs. B2 flats).

**Binning is uniform** at `'1,2'` with `window = '0,0,0,2048,1024'` for every
frame in the science configuration, so the second value bins the FITS row axis.
Whether that axis is spatial or spectral is *inconsistent inside PypeIt* for this
class -- see risk 7.

**Iodine.** All 10 July science frames have `iodin = T`. No iodine-free
HD 187123 frame exists in the July run; the candidate templates are the three
500 s frames on 1998-08-26.

#### Other runs (summary)

1997-12-23/24 (U07H): 1 science frame each, with 2 matching B1 arcs and 2
matching B1 flats per night, no biases. PypeIt will type these as `RED97`
(`mjd < 50814`), i.e. a separate configuration. 1998-06-18 (N12H): 1 matching
arc, 2 matching flats. 1998-08-17..26 (U05H): 1-2 matching arcs and 1-3 matching
flats per night, with biases/darks present. 1998-09-12..18: 1 matching arc per
night except 09-17, which has none (and whose science frame is B2).

#### Risks carried into prompts 5-9

1. **July 19 has no calibrations at all** and needs a cross-night assignment from
   July 18 in the `.pypeit` file.
2. **Clean flats are B2, science is B1**; `decker` is a configuration key, so
   expect to hand-assign flats or drop `decker` from the configuration keys.
3. **The only B1 flats have the iodine cell in the beam** -- do not let PypeIt
   use them as pixel flats without thought.
4. **62 of 70 flats have `hatopen = T`**, which may leave them untyped.
5. **No biases on the run nights**; the 1998-07-14 set is the fallback.
6. **Arcs are 10 s with hundreds of saturated pixels**, and there are only 1-2
   per night.
7. **The BINNING axis convention is internally inconsistent** for
   `keck_hires_orig`: `get_detector_par` sets `specaxis = 1`, while the inherited
   `compound_meta('binning')` does `binspatial, binspec = parse_binning(...)`, a
   convention written for the post-2004 mosaic. This feeds
   `wavelengths.fwhm = 8.0/bin_spec` and `order_spat_range`. A real 2048x1024
   frame settles it (prompt 5). Related: `config_specific_par` scales
   `order_spat_range` by 6200 px, a mosaic-sized number, on a 2048 px chip.
8. **Each July night is shared with Bida's program** (N01H) and with 68-111
   frames of other planet-search targets in the *same* configuration. A download
   must be filtered to HD 187123 plus its calibrations, or `pypeit_setup` will be
   swamped.

### Prompt 4 (2026-09-20): what `keck_hires_orig` expects, and every gap against the archive

Source read at `develop` `f3a1f1d274b15ee1358f167819d77f1948fce1bd`,
`pypeit/spectrographs/keck_hires.py`. **Nothing was edited.** Frame-typing
predictions below were *simulated* by replaying `compound_meta` and
`check_frame_type` over the real KOA rows in `koa_hd187123_nights_all.csv`, not
guessed.

**Inheritance, settled at runtime.** `KeckHIRESOrigSpectrograph` inherits from
`KECKHIRESBaseSpectrograph`, **not** from the post-2004 `KECKHIRESSpectrograph`.
MRO: `['KeckHIRESOrigSpectrograph', 'KECKHIRESBaseSpectrograph', 'Spectrograph',
'object']`. So `config_specific_par` resolves to the generic `Spectrograph` one,
and the mosaic machinery (`get_mosaic_par`, `hires_read_1chip`, `indexing`) is
never reached.

#### Header cards `init_meta` uses

Base `init_meta` L188-218, Orig override L937-947. Key entries: `ra`/`dec` from
`RA`/`DEC` (the only ones `required` for science/standard); **`target` from
`OBJECT`** (Orig L947 overriding the base `TARGNAME`); `decker` from `DECKNAME`;
`exptime` from `ELAPTIME`; `dateobs` from `DATE-OBS`; `hatch` from `HATOPEN`;
`filter1` from `FIL1NAME`; `echangle` from `ECHANGL` (rtol 1e-3, atol 1e-2);
`xdangle` from `XDANGL` (rtol 1e-2, atol 1e-1); `frameno` from `FRAMENO`;
`instrument` from `INSTRUME`. Compound: `binning`, `mjd`, `dispname`, `idname`,
`lampstat01`. `raw_header_cards()` = `['FIL1NAME','ECHANGL','XDANGL']`.

The commented-out `idname`/`IMAGETYP` line at L945 is dead: **`idname` is derived
entirely from the lamp/cover/hatch/shutter logicals** at L259-297, not from a
card.

#### Frame typing — the decisive logic

`idname` (L275-297, non-RED97 branch), paraphrasing the actual code:

    collcoveropen = (dispname=='RED' and RCCVOPEN) or (dispname=='UV' and BCCVOPEN)
    if XCOVOPEN and collcoveropen and no LAMPCAT1/2 and no LAMPQTZ2 and LAMPNAME!='quartz1':
        if HATOPEN and AUTOSHUT:      -> 'Object'
        elif not HATOPEN:             -> 'Bias' if ELAPTIME < 0.001 else 'Dark'
    elif XCOVOPEN and collcoveropen and AUTOSHUT and (LAMPCAT1 or LAMPCAT2):
                                      -> 'Line'
    elif collcoveropen and AUTOSHUT and (LAMPQTZ2 or LAMPNAME=='quartz1') and not HATOPEN:
                                      -> 'slitlessFlat' if not XCOVOPEN else 'IntFlat'
    (otherwise falls off the end -> None)

`check_frame_type` (L367-414) then maps `Object`->science/standard,
`Line`->arc+tilt, `IntFlat`->pixelflat+illumflat+trace, `Bias`->bias,
`Dark`->dark, each gated by an exposure range.

**The hatch question from prompt 3 is answered: `not HATOPEN` is mandatory for
every flat ftype, and is NOT tested for arcs.**

Simulated typing of the real July 15-19 rows:

| frame class | `idname` | automatic frametype |
|---|---|---|
| HD 187123 science, 215-427 s | `Object` | **NONE** (see the exposure-range gap) |
| ThAr1 B1 10 s arcs | `Line` | `arc,tilt` |
| B2 2 s quartz flats, hatch open (56 of 64) | `None` | **NONE** |
| B2 2 s quartz flats, hatch closed (8, on 07-16) | `IntFlat` | `pixelflat,illumflat,trace` (B2 setup) |
| B1 3 s iodine flats (all hatch open) | `None` | **NONE** |
| D5 focus frames (ThAr) | `Line` | `arc,tilt` (own D5 setup) |
| Bida's N01H frames (83) | `None` | **NONE** (`RCCVOPEN=F`) |
| 1998-07-14 biases/darks | `None` | **NONE** (covers closed) |

#### Detector parameters (Orig `get_detector_par`, L968-1012)

`dataext=0`, `specaxis=1`, `specflip=False`, `spatflip=False`,
`platescale=0.216` "/unbinned pix, `darkcurr=0.0`, `saturation=65535`,
`nonlinear=0.7`, `mincounts=-1e10`, `numamplifiers=1`, `ronoise=[2.8]`.
**`gain` is header-dependent**: `CCDGAIN == False` -> 1.9 e-/ADU, `== True` ->
0.78, anything else -> `PypeItError("Bad CCDGAIN mode for HIRES")`. Our frames
have `ccdgain = F`. `datasec`/`oscansec` are `None` here and supplied per-file by
`get_rawimage`. `nonlinear_counts` works out to ~87,200 e-, which is what
wavecalib uses to reject saturated arc lines.

#### Default `PypeItPar` and required calibrations

Orig `default_pypeit_par` (L950-966) adds only `rdx.detnum = [1]` and
`wavelengths.ech_separate_2d = False` on top of the base. Resolved at runtime:

- `wavelengths`: `method = echelle`, `lamps = ['ThAr']` (matching our ThAr1
  arcs), `fwhm = 4.0` with `fwhm_fromlines = True`, `ech_2dfit = True`,
  `sigdetect = 5`, `rms_thresh_frac_fwhm = 0.1`, `refframe = heliocentric`.
- `use_biasimage = False` and `use_darkimage = False` for every frame type;
  `use_overscan = True` (median); science/standard have `use_pixelflat = True`
  and `use_illumflat = True`; `raise_chk_error = True`.
- **Therefore the mandatory calibrations per group are `arc`, `tilt`, `trace`,
  `pixelflat` and `illumflat`. Bias and dark are NOT required.**
- `slitedges`: `order_spat_range = None`, which resolves to `[0, nspat]` =
  `[0, 1024]` for our chip — the correct chip-sized value.

#### `get_echelle_angle_files` — not a blocker

Returns `keck_hires_orig_angle_fits.fits` and `keck_hires_composite_arc.fits`.
**Both are present on disk** in `pypeit/data/arc_lines/reid_arxiv/` (25,920 B and
15,577,920 B) and git-tracked. No `pypeit_cache_github_data` download needed. The
orig angle file covers `RED` with our echangle/xdangle inside its range.

#### Updated gap list

1. **`scienceframe.exprng = [601, None]` vs. our 215-427 s science frames** —
   **the single most certain blocker for the automatic pass.** All 10 HD 187123
   frames fail the exposure floor and come out untyped, so `pypeit_setup` will
   comment them out. Remedy: set `frametype = science` for the 10 frames in the
   `.pypeit` data block (`run_pypeit` honours the user column and does not
   re-apply `exprng`), optionally with `[scienceframe] exprng = 100, None` for
   the record. The 601 s floor is a post-2004 faint-target assumption that does
   not fit bright planet-search stars — worth reporting upstream.
2. **Flats need the hatch closed; 62 of 70 July flats have it open** — hand-type
   in the `.pypeit` file. Confirms risk 4.
3. **No clean B1 flat in the July 15-19 run — but 1998-07-14 has four.**
   `HI.19980714.07711`, `.07810`, `.07910` (02:08-02:11 UT) and `.55408`
   (15:23 UT) are 3 s `Narrowflat` frames: decker B1, `iodin = F`, **hatch
   closed**, `echangl` 0.00106-0.00114, `xdangl` -0.5420 — inside the angle
   tolerances of the whole July 15-19 run, and typed `IntFlat` automatically.
   That night also carries 8 hatch-closed B1 ThAr arcs. **This is the cleanest
   fix for the flat problem and was not visible in the prompt-3 survey.**
   (1998-08-12 has ~18 more clean B1 flats if wanted.)
4. **The B1 iodine flats carry the I2 forest and PypeIt has no `IODIN` logic** —
   nothing stops a user typing them `pixelflat`; use them as `trace` only, if at
   all. Confirms risk 3.
5. **Biases/darks are not required at all** (`use_biasimage = False`), so risk 5
   is largely void. If wanted anyway, the July-14 biases are untyped because the
   Bias/Dark branch requires the covers *open* and theirs were shut — they would
   need hand-typing.
6. **July 19 still needs cross-night calibrations**; the code imposes no date
   constraint on calibration groups, so a `calib` column assignment works.
   Confirms risk 1.
7. **Saturated arcs are handled**: wavecalib rejects lines above
   `nonlinear_counts` (~87,200 e-), and the echelle method cross-correlates
   against the archived composite arc, so one arc per night suffices. Refines
   risk 6 downward.
8. **The binning inversion is real but narrow.** Header `BINNING = '1,2'` becomes
   the PypeIt string `'2,1'` (binspec 2, binspat 1), while `specaxis = 1` means
   the *column* axis (binned x1) is spectral and the *row* axis (binned x2) is
   spatial — i.e. physically binspec 1, binspat 2, the opposite of the label.
   Consequences are limited: the echelle wavelength solution does not use
   `binspectral` at all, and `bpm` is correct *because* it inverts the naming a
   second time (L1126). What is affected is `order_platescale`
   (0.216 instead of 0.432 "/binned pixel, so anything in arcsec, including the
   1.5" boxcar radius, is off by 2x) and the `detector.binning` string written
   into spec1d/spec2d headers. **Not fixable from a `.pypeit` file** — it is an
   upstream fix, one for Ryan Cooke.
   **The `config_specific_par` / 6200-px half of prompt-3 risk 7 is void**: that
   method belongs to the post-2004 class and is unreachable here.
9. **`get_rawimage` reads `PREPIX`, not `PRECOL`** (L1039) — and the KOA metadata
   only exposes `PRECOL = 21`. If the 1998 headers lack a `PREPIX` card this is a
   hard `KeyError` on every frame. **Undecidable without a real file; prompt 5
   must check.** Expected layout: `NAXIS1 = 2*PREPIX + 2048 + POSTPIX`.
10. **Header logicals must be genuine FITS booleans.** All the typing logic and
    the gain selection do truthiness tests on raw cards. If `HATOPEN`, `AUTOSHUT`,
    `XCOVOPEN`, `RCCVOPEN`, `LAMPCAT1/2`, `LAMPQTZ2` or `CCDGAIN` are stored as
    the *strings* `'T'`/`'F'`, science frames would type as `Line` and
    `get_detector_par` would raise "Bad CCDGAIN mode". **Prompt 5 must print the
    card types.**
11. **`mjd` needs an `MJD` card, or `DATE-OBS` + `UTC`.** If neither exists the
    error is swallowed and `mjd = None`, which cascades into `dispname` and
    `idname` failing and the frame going untyped. Prompt 5 to confirm.
12. **`target` will read `Star+Iodine`**, not `187123`, because the Orig class
    takes `target` from `OBJECT`. Output basenames inherit it. Cosmetic;
    editable in the `.pypeit` data block.
13. **Filter the download.** Bida's 83 frames and 300+ other planet-search frames
    per night all arrive untyped, and the D5 focus frames would form a spurious
    arc-only setup. Confirms risk 8.
14. **Upstream nits banked**: `check_spectrograph` L567 constructs a
    `PypeItError` without `raise` (a no-op); `config_independent_frames`'
    docstring claims a `DATE-OBS` rule the code does not implement; `spectrim`,
    `PRELINE` and `POSTLINE` are read but unused in the Orig reader.

### Prompt 5 (2026-09-20): downloaded 1998-07-16 (+ 1998-07-14 flats) and checked the headers

**31 frames, 152.8 MB, all verified**, outside the repository in
`/Users/xavier/Projects/PypeIt/first-hires-exoplanet-data/raw/`. No FITS file is
inside the git repo (`find ... -name '*.fits'` returns 0). Download script:
`first_hires_exoplanet/koa_download.py`, idempotent and refusing any output root
under the repo.

**Night chosen: 1998-07-16**, because it is the only night of the July 15-19 run
with *hatch-closed* flats (8 B2 quartz frames at 15:08-15:20 UT), which PypeIt
types automatically, and it has two ThAr arcs bracketing the night. The four
clean B1 flats and two B1 arcs from 1998-07-14 identified in prompt 4 were
downloaded alongside.

**Retrieval recipe.** Query KOA for the `filehand` column (it is *not* in the
saved survey CSVs), then
`GET https://koa.ipac.caltech.edu/cgi-bin/getKOA/nph-getKOA?filehand=<filehand>`,
which returns `Content-Type: image/x-fits`, 4,743,360 B per 1-amp frame. A wrong
`filehand` silently returns an HTML error page with **status 200**, so every
download must be verified (content type, size, `SIMPLE`, FITS-openable). Note
`astropy.io.fits` refuses `memmap=True` on these files ("BZERO/BSCALE/BLANK
header keywords present") — use `memmap=False`.

#### Manifest

`raw/1998jul16/` (21 files, each 4,743,360 B, raw shape (1024, 2303)):

| koaid | role | UT | elap | decker | hatch | iodin |
|---|---|---|---|---|---|---|
| HI.19980716.52985 | **science** (HD 187123) | 14:43:05 | 400 | B1 | open | T |
| HI.19980716.27559 | ThAr arc | 07:39:19 | 10 | B1 | **open** | F |
| HI.19980716.54321 | ThAr arc | 15:05:21 | 10 | B1 | closed | F |
| .26537 .26659 .26754 .26850 .26945 .27041 .27137 .27233 | quartz flat | 07:22-07:34 | 2 | B2 | open | F |
| .54482 .54599 .54696 .54791 .54887 .54982 .55079 .55175 | quartz flat | 15:08-15:20 | 2 | B2 | **closed** | F |
| HI.19980716.27384, .54172 | iodine flat | 07:36, 15:03 | 3 | B1 | open | T |

`raw/1998jul14/` (9 files): the 4 clean B1 `Narrowflat` frames
(`.07711 .07810 .07910` at 02:08-02:12, `.55408` at 15:23), 2 hatch-closed B1
ThAr arcs (`.07383`, `.55506`), and 3 zero-second biases (`.55653 .55747 .55840`).

`raw/inspect_only/` holds the 07-16 `bias_lamp_on` frame, which turned out to be
Bida's 2-amplifier N01H frame; it is kept out of the reduction directories so
`pypeit_setup` never sees it, but it proved decisive (below).

#### Prompt 4's open header questions, answered

**Gap 10 (header logicals) — PASS.** They are genuine FITS booleans, not strings:
on the science frame `HATOPEN = True`, `AUTOSHUT = True`, `XCOVOPEN = True`,
`RCCVOPEN = True`, `LAMPCAT1 = False`, `LAMPQTZ2 = False`, `CCDGAIN = False`,
each `<class 'bool'>`. So the typing logic and the gain lookup behave as written,
and gain 1.9 e-/ADU is selected without error.

**Gap 11 (`MJD`) — PASS with a nuance.** `MJD = '51010.613267'` exists but as a
**string**; PypeIt casts it, and `get_meta_value('mjd')` returns a float.
`UTC` is **absent** (the card is `UT`), and `DATE-OBS = '1998-07-16'` is
date-only. The `DATE-OBS + UTC` fallback would raise `KeyError('UTC')` if it were
ever reached — it is not, because `MJD` is present. Anything keyed on `dateobs`
rather than `mjd` will see midnight.

**Gap 9 (`PREPIX`) — the card exists, but it exposed a genuine PypeIt bug.**
`PREPIX = 21`, `PRECOL = 21`, `POSTPIX = 234`, `NUMAMPS = 1`,
`WINDOW = '0,0,0,2048,1024'`, `NAXIS1 = 2303`, `NAXIS2 = 1024`. The arithmetic
prompt 4 predicted does **not** hold:

    2*PREPIX + 2048 + POSTPIX = 2324   != NAXIS1 = 2303
      PREPIX + 2048 + POSTPIX = 2303   == NAXIS1

Column medians of a hatch-closed flat confirm the real layout: columns 0-20 sit
at bias level (~778), column 21 jumps to 3533 and stays illuminated through
column 2068, and column 2069 drops back to bias (~784). **Data = columns
21:2069.** But `get_rawimage` (keck_hires.py:1050-1055) computes
`col0 = prepix * 2`, and at runtime returns `datasec cols 42..2089`,
`oscansec cols 2090..2302`. So **the data section is offset by exactly 21
columns**: it discards the 21 bluest real columns and includes 21 columns of
overscan at the other end of every order.

The 2-amplifier Bida frame explains the mistake: there,
`NAXIS1 = 2558 = 2*21 + 2048 + 2*234`, and `prepix * 2` is correct. **The code
should read `prepix * namp`, not `prepix * 2`** — a one-line upstream fix, and
one that only shows up on single-amplifier original-HIRES data. Verified
independently at runtime, not inferred.

**Gap 8 (binning axis) — settled empirically. `specaxis = 1` is correct.**
The raw array is (1024, 2303): axis 0 is 1024 rows (binned x2), axis 1 is 2048
chip columns plus pre/overscan (unbinned). Collapsing a B2 flat along columns
gives **25 peaks along the row axis**, median spacing 29.5 px, with 39% of the
profile dark between them — the echelle order stack. Collapsing along rows gives
a smooth blaze with no gaps. An arc shows 49 emission peaks along one order's
row. So **dispersion runs along NAXIS1 (2048, unbinned) and the x2-binned axis
is spatial**: physically binspec = 1, binspat = 2, while PypeIt reports
`binning = '2,1'`. Prompt 4's inversion is confirmed on real data.

#### PypeIt over the real files — everything ran

- `get_rawimage(science, 1)`: succeeds. Image (1024, 2303) float64, gain `[1.9]`,
  ronoise `[2.8]`, `binning '2,1'`, `specaxis 1`, `platescale 0.216`.
- `get_detector_par(1, hdu)`: gain `[1.9]`, `numamplifiers 1`.
- `bpm(science, 1)`: shape (2048, 1024), 0.67% masked.
- `check_spectrograph(science)`: does not raise.
- **Frame typing matches prompt 4's predictions 30/30.** Science -> `Object`;
  both arcs -> `Line` (*including the hatch-open 07:39 one*, confirming on real
  data that the hatch is not tested for arcs); hatch-closed B2 flats and the
  07-14 clean B1 flats -> `IntFlat`; hatch-open B2 flats and the iodine flats ->
  `None`; 07-14 biases -> `None` (covers closed).

#### Two further findings

**The hard-coded MAKEE bad-pixel mask does not land on the bad columns.** The
real bad columns are visible in the raw frames at 0-indexed NAXIS1 columns
106-107, 1148 and 2027-2028 — i.e. MAKEE's 1-indexed coordinates count *raw*
columns, including the 21 prescan columns. PypeIt applies them in the *trimmed*
frame, so the masked columns land roughly 41 px away with the current reader, and
would still be ~20 px off even after the `prepix * namp` fix. The 2027-2028 pair
is a strong hot pair (median 1078 vs 807 ADU in the 400 s science frame), so
watch for a spurious feature there in prompt 8.

**Gain and read noise disagree with the headers.** `DETECTOR = 'Tek 2048E LRIS
Eng grade device'`, with `CCDGN01 = 4.8` and `CCDRN01 = 6.0`, while PypeIt
hard-codes 1.9 e-/ADU and 2.8 e-. These KOA DQA cards may be generic rather than
measured for 1998, so this is **a question for Ryan Cooke, not a bug claim** —
but it is worth asking before we quote a signal-to-noise in prompt 8.

#### Carried into prompts 6-7

1. The 400 s science frame is `Object` but fails `exprng = [601, None]`, so it
   will be untyped — hand-type it (prompt 4 gap 1).
2. Configurations will split on `decker`: B1 (science, arcs, iodine flats, 07-14
   clean flats) vs B2 (the 07-16 quartz flats). The 07-14 angles are inside
   tolerance of 07-16, so the 07-14 B1 frames should merge into the B1
   configuration — verify in prompt 6.
3. 13 of the 30 reduction frames arrive untyped (8 hatch-open flats, 2 iodine
   flats, 3 biases). Expected.
4. **The 21-column reader offset will affect every processed frame.** A 21 px
   shift should be tolerable for flat-fielding and for the composite-arc
   cross-correlation, but the dead strip may trigger low-count masking at the
   order ends and slightly shrink the usable wavelength range. The fix is
   upstream, not in the `.pypeit` file.
5. The bad columns are effectively unmasked.

### Prompt 7 (2026-09-20): what `pypeit_setup` produced automatically

Run with PypeIt `2.0.2.dev1217+g017bece06` (branch `orig-hires-fixes`, commit
`017bece06`, i.e. **including the prompt-6 fixes**). **No `.pypeit` file was
hand-edited and no parameter overrides were passed** — everything below is the
honest automatic result. All products are outside the repository, under
`first-hires-exoplanet-data/redux/`.

**Invocation.** `-r/--root` accepts multiple paths, so no symlink tree was
needed. One quirk worth recording: `pypeit_setup` writes **either** the
per-configuration `.pypeit`/`.calib` files (with `-c`) **or**
`setup_files/*.sorted`/`.obslog` (without `-c`), never both, so each run needs
two passes:

    pypeit_setup -s keck_hires_orig -r <raw...> -c all -d <outdir>   # .pypeit files
    pypeit_setup -s keck_hires_orig -r <raw...>       -d <outdir>   # .sorted/.obslog

Both runs exited 0 with `Found 2 unique configuration(s).` and
`PypeIt file successfully vetted.` The only warnings were a benign
`No parameters necessary for median overscan method` and the expected
`Couldn't identify the following files:` list.

#### Run A — the night alone (21 frames)

Two configurations, split on `decker` alone:

| Setup | dispname | decker | filter1 | echangle | xdangle | binning |
|---|---|---|---|---|---|---|
| A | RED | **B2** | kv370 | 0.00093024 | -0.54400003 | 1,2 |
| B | RED | **B1** | kv370 | 0.00093024 | -0.54400003 | 1,2 |

Setup B (the science setup) contains the science frame and the two B1 ThAr arcs
— and **nothing else**. Its `.calib` lists only `arc`, `tilt` and `science` for
calibration group 1: **no `trace`, no `pixelflat`, no `illumflat`.** Setup A
holds the 8 hatch-closed B2 flats as `pixelflat,illumflat,trace` (and the 8
hatch-open ones commented out), but has no science frame.

**So the night on its own cannot be reduced**, and `pypeit_setup` does not warn
about it: `inputfiles.vet` checks file syntax, not calibration completeness, and
still reported "successfully vetted". `run_pypeit` on this file would fail at
order tracing. That is the single most important result of this prompt.

#### Run B — the night plus the 1998-07-14 calibrations (30 frames)

| Setup | dispname | decker | filter1 | echangle | xdangle | binning |
|---|---|---|---|---|---|---|
| **A** | RED | **B1** | kv370 | 0.00114165 | -0.542 | 1,2 |
| B | RED | B2 | kv370 | 0.00093024 | -0.54400003 | 1,2 |

(The setup letters swap relative to Run A, and Setup A's recorded angles are
those of the *first frame encountered* — the 07-14 arc — not the science
frame's. Informational only; matching is by tolerance.)

**The prediction from prompts 4 and 5 holds: the 1998-07-14 clean B1 flats and
arcs land in the same configuration as the 07-16 science frame.** The angle
differences are |dechangle| = 2.1e-4 and |dxdangle| = 0.0035, far inside the
`atol` of 0.01 and 0.1. Setup A's data block, as generated:

    # HI.19980714.55653/.55747/.55840  None                       Bias        (3 biases)
    # HI.19980716.27384/.54172         None                       iodine      (2 iodine flats)
      HI.19980714.07383/.55506         arc,tilt                   Th-Ar
      HI.19980716.27559/.54321         arc,tilt                   Th-Ar
      HI.19980714.07711/.07810/.07910/.55408
                                       pixelflat,illumflat,trace  Narrowflat
      HI.19980716.52985                science                    Star+Iodine

Calibration group 0 is **complete**: 4 arcs (both nights), 4 flats serving
`pixelflat`/`illumflat`/`trace`, and the science frame. **Nothing needs
hand-typing for this to run.**

A physical check confirmed the cross-night merge is safe: cross-correlating the
spatial order profile between 07-14 and 07-16 frames gives shifts of 0 to -1
binned spatial pixels (the two 07-16 iodine flats, 7.5 h apart, differ by -1 px
themselves), and the arc dispersion shift between nights is 0 px.

#### Frame typing: 17 typed, 13 untyped, 0 mis-typed

| frames | why untyped | verdict |
|---|---|---|
| 8 hatch-open B2 flats (2 s) | `idname` flat branch requires `not HATOPEN` (keck_hires.py:290-293); with the hatch open the frame falls off the end of every branch -> `None` | correctly left alone; B2 is not the science slit |
| 2 B1 iodine flats (3 s) | identical mechanism | hand-typeable as `trace` only if wanted; **not needed**, the 07-14 flats supply it |
| 3 biases (0 s) | `collcoveropen` is False (covers shut), so the Bias/Dark branch at L283-286 is unreachable | correctly left alone; `use_biasimage = False`, biases are not required |

**No frame is mis-typed.** Every automatic assignment matches prompt 4's
simulation and prompt 5's 30/30 verification.

#### The prompt-6 fixes in a real pipeline run

- **Fix 5 (science `exprng`) works**: `HI.19980716.52985` (400 s) is typed
  `science` — and *only* `science`, since the archive-standard positional test
  returns False, so `vet_assigned_ftypes` had nothing to arbitrate. Before the
  fix this frame would have been commented out.
- **Fix 3 (binning) works**: the `binning` column reads `1,2` for all 30 frames,
  and participates in configuration matching without incident.
- Fixes 1, 2, 4 and 6 are not exercised by `pypeit_setup`, which reads only
  headers. They come into play in prompt 8.

**No new bugs were introduced by the prompt-6 fixes.** Two PypeIt-wide
observations worth banking: the `-c`/no-`-c` either/or described above, and the
fact that `vet` will pass a `.pypeit` file whose science calibration group has no
flats at all.

#### What prompt 8 should use

`redux/setup_jul16_jul14/keck_hires_orig_A/keck_hires_orig_A.pypeit`, copied into
a fresh reduction directory so the automatic result stays on disk untouched.
**No hand-typing is required.** Optional edits only: set `target` to `HD187123`
for sane output basenames; optionally drop the two 07-14 arcs if single-night
arcs are preferred (the pixel check says stacking is harmless). Do **not**
reduce Setup B — it has flats but no science.

### Prompt 8 (2026-09-20): the reduction runs — 37 orders extracted

PypeIt `2.0.2.dev1217+g017bece06` (branch `orig-hires-fixes`). Working directory
`first-hires-exoplanet-data/redux/reduce_jul16/`, outside the repository. Every
step below was reproduced from the pristine prompt-7 output; each parameter was
added one at a time and verified.

#### The baseline fails silently

`run_pypeit` on the unmodified prompt-7 file **exits 0 and prints "Data reduction
complete"**, but writes only a `spec2d` — **no `spec1d` at all**. The log carries
`No objects found automatically` 28 times, and
`findobj_skymask.py:2144` explains itself: *"No objects were found because the
image was heavily masked, not because no source was detected."* Runtime 5 min.

**A reduction that produces no 1D spectrum should not exit 0.** That is the
first genuine issue this prompt found, and it is PypeIt-wide rather than
HIRES-specific.

#### Mechanism: the orders are too narrow for the default edge trim

Measured from `Slits_A_0_DET01`: the B1 (3.5") orders are **10.32 binned pixels
wide** (min 9.47, max 14.21). The inherited default `find_trim_edge = [3, 3]`
masks six of those ten columns, leaving **4.3 px** for object finding.

This is a direct consequence of the 2x spatial binning: at 1x1 the orders would
be ~21 px and `[3, 3]` would be unremarkable. It surfaces now precisely *because*
prompt 6 fixed the binning inversion, so `order_platescale` finally returns the
correct 0.432"/binned pixel. The default is inherited from the post-2004 mosaic
class, where the slits are much wider in pixels.

#### Three parameters, each justified by its own run

| run | parameters | result |
|---|---|---|
| baseline | (automatic file) | 28/37 orders find nothing; **no spec1d** |
| v1 | `find_trim_edge = 1,1` | 25/37 still fail; **still no spec1d** — necessary but **not sufficient** |
| v2 | + `skip_skysub = True` | **spec1d written**, 37 orders; 34 good, but orders 84, 81, 77 at S/N 17.9, **-1.9**, 17.5 |
| **v3** | + `no_local_sky = True` | **37/37 orders clean**; those three recover to S/N 94, 108, 125 |

`skip_skysub` and `no_local_sky` have the same physical cause: HD 187123 (V=7.9)
in a 400 s exposure **fills the ~10 px slit**, so neither the global nor the
local sky fit has object-free pixels to work with. The sky is well under 1% of
the stellar peak here, so skipping it costs nothing.

Notably, `min_frac_prof = 0.5` was **not** needed — the default 0.9 gives clean
extractions once the sky handling is right.

#### Final result (v3)

- **37 echelle orders extracted, 93 down to 57**, one object (`OBJ0525-DET01`).
- **S/N 33 at the blue end rising to 164 in the red**, median **124.5**; no order
  below 20.
- **Wavelength solutions for all 37 orders**, RMS **0.069-0.195 px**, median
  0.133; **none above 0.2 px**. The 10 s arcs with saturated lines were not a
  problem — `nonlinear_counts` rejection plus the archived composite arc did the
  job, exactly as prompt 4 predicted.
- Calibrations built without incident: `Edges`, `Slits`, `Flat`, `Arc`,
  `Tiltimg`, `Tilts`, `WaveCalib`. QA HTML written for the calibration group and
  the science frame.
- Runtime: 5 min for the full pass including calibrations; ~40 s for a
  science-only re-run.

The final input file is `keck_hires_orig_A.pypeit`, with the intermediate
versions kept as `.v0_auto`, `.v1_trimedge`, `.v2_skipsky`, `.v3_nolocalsky` so
the path from the automatic output is auditable.

#### Genuine issues for upstream (items 5-7 for Ryan Cooke)

5. **`run_pypeit` exits 0 having written no `spec1d`.** Object finding failed in
   every order and the run still reported success. A non-zero exit, or at least a
   prominent warning, would have saved the whole diagnosis.
6. **`find_trim_edge = [3, 3]` is untenable for `keck_hires_orig` at 2x spatial
   binning**, where the B1 orders are ~10 px wide. Worth a binning-aware default,
   or an Orig-class override, rather than leaving every user to rediscover it.
7. **The `spec1d_*.txt` summary is appended to, not overwritten, across
   `run_pypeit -o` runs**, silently accumulating duplicate tables (74 rows for 37
   orders after two runs). Minor, but it will mislead anyone parsing it.

Nothing here contradicts the prompt-6 fixes; they behaved as intended throughout.

### Prompt 9 (2026-09-21): inspection of the 1998-07-16 reduction

Figures: `first_hires_exoplanet/figs_phase1.py` ->
`docs/figs/fig_p1_{orders,wavecal,iodine,spectrum}.png`. The script reads the
reduction from the sibling data tree (`--redux` to point elsewhere) and refuses
to run against a directory holding more than one accumulated pass.

**Verdict: the reduction is sound.** Every check below passes, and two of them
are independent confirmations rather than internal consistency.

#### Order tracing and flat field

**37 orders traced, none masked.** Order widths are tight and vary smoothly
with position: median **10.32 binned px**, 9.47 to 10.5 across the detector,
with a single 12.8 px outlier at the extreme edge. That smooth run is itself
evidence the trace is following real orders rather than noise.

The pixel flat over the **0.1M illuminated pixels** has 1-99 percentile
**0.959-1.024**, i.e. pixel-to-pixel response flat to about +/-3%, single-peaked
and near-symmetric about 1.0. (Taken over the whole frame the spread looks like
0.986-1.003, but that is meaningless: pixels outside the traced orders are left
at exactly 1.0 and swamp the histogram.)

#### Wavelength solution

**All 37 orders solved, none above 0.2 px RMS**: 0.069 to 0.195, median
**0.133 px**, with no trend from blue to red. The 10 s arcs with hundreds of
saturated pixels were not a problem.

**The zero point is independently correct.** Measuring seven strong photospheric
lines against their *vacuum* rest wavelengths gives **-12 to -16 km/s,
mean -14.7 km/s**, consistent across 3935-5185 A. HD 187123's systemic velocity
is about **-17 km/s**, and `VEL_CORR` shows the +4.2 km/s heliocentric
correction was applied. So the wavelength scale is right to a few km/s against
an external truth, not merely internally self-consistent.

**A trap worth recording: PypeIt reports vacuum wavelengths.** Comparing them
against the air rest wavelengths conventionally tabulated for the optical
produces a spurious **+70 km/s** redshift (+1.2 A at 5170 A). This is noted in
the figure script's docstring so nobody repeats it in phase 2.

#### The iodine forest

Present and unmistakable. Counting absorption features per 100 A order by order:

- **3833-4951 A**: density falls smoothly from 356 to 162 lines/100 A -- the
  normal stellar metal-line crowding, thinning toward the red.
- **5021 A**: jumps to 282, then **5093-5750 A** sits on a plateau of
  **386-432 lines/100 A** -- roughly **double** the extrapolated stellar trend.
- **5844-6236 A**: falls back through 374, 298, 294, 248 to 180 as the band ends.

That is the I2 band at **~500-620 nm**, exactly where the cell puts it. The
side-by-side order panels make the same point directly: order 73 (4880 A) shows
resolved, well-separated stellar lines on a flat continuum; order 65 (5484 A)
shows the dense quasi-continuous forest.

#### Stellar spectrum and signal-to-noise

Coverage is continuous from **3806 to 6262 A** (381-626 nm) over 37 orders. Every
strong G2V marker in range is present and deep: **Ca II H and K, H-gamma,
H-beta, and the Mg b triplet**, each at its expected (vacuum) wavelength.
Na D at 5890/5896 A falls in an inter-order gap -- the free spectral range no
longer overlaps at the red end, which is expected for this setup, not a defect.

**S/N per pixel rises monotonically from 33 in the bluest order to 168 in the
reddest, median 124.** The rise is smooth with no discontinuities, as expected
for a G2V star through this cross-disperser.

#### Residual blemishes

Bad-pixel-masked columns come through the extraction as exact zeros rather than
as flagged gaps. The figure script converts non-positive flux to NaN so they
plot as gaps; anyone consuming `spec1d` numerically should mask `OPT_COUNTS <= 0`
as well.

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

### 2026-09-19 (Prompt 2: refreshed the install, verified `keck_hires_orig`, recorded the commit)

**Task.** Carry out prompt 1's recommendation, verify `keck_hires_orig` is
importable and appears in the spectrograph listing, and record the commit hash.
Delegated to a Fable subagent as the prompt asked; re-verified here.

**The commit to quote.** Every reduction we publish from this environment should
be cited as:

> PypeIt `2.0.2.dev1216+gf3a1f1d27` -- branch `develop`, commit
> `f3a1f1d274b15ee1358f167819d77f1948fce1bd`
> ("Merge pull request #2197 from pypeit/keck_hires_tektronix", 2026-09-14),
> editable install from `/Users/xavier/Projects/PypeIt/PypeIt`.

**Version metadata: fixed and confirmed.** The recommended
`conda run -n pypeit14 pip install -e . --no-deps` was run in the PypeIt
checkout. It succeeded, rebuilt the editable wheel, and rewrote the gitignored
`pypeit/pkg/version.py` to `__version__ = '2.0.2.dev1216+gf3a1f1d27'` -- HEAD's
short sha, with no `.dirty` suffix. The install is still editable and still
resolves to the checkout, and the tracked tree stayed clean. Network was needed
(build isolation fetches `setuptools`, absent from the env) and was available.

One honest wrinkle: by the time the subagent took its "before" reading, the
version was *already* `...+gf3a1f1d27`. The prompt-1 measurement of
`2.0.2.dev891+gf8a757720` was real, and `version.py` was regenerated between the
two tasks -- so the fix had landed before this task re-applied it. The reinstall
is idempotent, so re-running it changed nothing and confirmed the state. The net
effect is what matters: the version now matches HEAD, verified.

**Verification.** `keck_hires_orig` is importable and listed:

- `conda run -n pypeit14 python -m first_hires_exoplanet.check_env` -- all checks
  pass, exit 0 (see below).
- CLI: `conda run -n pypeit14 pypeit_setup -h | grep keck_hires_orig` prints
  `keck_deimos, keck_esi, keck_hires, keck_hires_orig,` in the
  `-s/--spectrograph` choices. 87 spectrographs are registered.
- Runtime: `load_spectrograph('keck_hires_orig')` gives `ndet = 1`,
  `supported = False`, `configuration_keys()` =
  `['dispname', 'decker', 'filter1', 'echangle', 'xdangle', 'binning']`,
  `get_detector_par(1)` -> platescale 0.216, ronoise [2.8], saturation 65535,
  and `default_pypeit_par()['rdx']['detnum']` -> `[1]`.

**New file: `first_hires_exoplanet/check_env.py`.** Rather than leave this
verification as a one-off, it is now a committable, re-runnable script, run from
the repo root with:

    conda run -n pypeit14 python -m first_hires_exoplanet.check_env

It prints the PypeIt version and import path, infers the checkout from
`pypeit.__file__`, reports HEAD sha / branch / commit subject / tree
cleanliness, checks that the version string actually embeds HEAD (warning loudly
and exiting 1 if not, which is precisely the trap we fell into), and confirms
`keck_hires_orig` loads and builds its detector and default parameters. It exits
non-zero on any failure, so it can gate a reduction. Style follows `figs.py`;
standard library plus `pypeit` only.

**What I learned about the repository / environment.**

1. *The stale-version trap is a standing hazard, not a one-off.* Because
   `pypeit/pkg/version.py` is written by `setuptools_scm` **at install time** and
   is gitignored, every future branch switch in the PypeIt checkout will silently
   desynchronise `pypeit.__version__` from the code actually running. Since that
   string is stamped into PypeIt output headers, a reduction can be mislabelled
   with no visible symptom. `check_env.py` exists to catch this; **run it before
   any reduction whose results we quote** (prompts 7, 8, 9).
2. *There is no spectrograph-listing script in PypeIt at all.* Re-confirmed:
   `pypeit_show_spectrographs` (named in this document's prompt 2) does not
   exist, and the only `pypeit_show_*` entry points are `show_arxiv`,
   `show_wvcalib`, `show_1dspec`, `show_2dspec`, `show_pixflat`. The spectrograph
   list surfaces *only* as the `-s/--spectrograph` choices text in
   `pypeit_setup`, `pypeit_obslog`, `pypeit_chk_for_calibs`, `pypeit_ql`,
   `run_pypeit`, `pypeit_trace_edges`, `pypeit_cache_github_data` and
   `pypeit_view_fits`, all of which import `available_spectrographs` from
   `pypeit.spectrographs.util`. The canonical check is the runtime import or
   `pypeit_setup -h`.
3. *`pip show pypeit` dumps the entire BSD licence text into its `License:`
   field.* Harmless, but it makes that command a poor way to eyeball the version;
   `python -c "import pypeit; print(pypeit.__version__)"` is the clean one.
4. *The repo now has a second module alongside `figs.py`.* `check_env.py` is the
   first non-figure code in `first_hires_exoplanet/`, and it runs as
   `python -m first_hires_exoplanet.check_env` because `__init__.py` already
   exists. It is untracked (`??`) -- staging and committing is yours.

**No PypeIt source file was edited and no git command changed state.** Changes
in this repo: the new `first_hires_exoplanet/check_env.py`, and this document
(the "Resolved" note updated with the new version string, plus this entry).

### 2026-09-20 (Prompt 3: surveyed KOA for the 1997-98 HD 187123 frames and their calibrations)

**Task.** Survey the Keck Observatory Archive for HIRES observations of
HD 187123 between 1997 December and 1998 October and for the calibrations taken
on those nights; report frame types, counts, binning, angles, deckers, exposure
times and public status; establish whether ThAr arcs exist for the 1998 July
15-19 run; download nothing. Delegated to a Fable subagent as the prompt asked.
The findings are written into the `## Report` section above. **No FITS files
were downloaded.**

**Headline results.**

- **ThAr arcs DO exist for the July run, in the science configuration** -- on
  four of the five nights. This was the biggest open risk in "Risks worth probing
  early" ("Calibrations may be thin ... If PypeIt needs arcs that were never
  taken, that shapes the whole phase"), and it has landed the good way.
- **July 19 is the exception: no arc, no flat, no bias.** It will need the
  July 18 arc, assigned across nights by hand.
- **The flats are awkward.** The clean quartz flats are decker B2 while the
  science is B1, and the only B1 flats have the iodine cell in the beam.
- **The archive confirms the document's own epoch count**: exactly 30 of the 37
  HD 187123 frames match the 30 catalogue epochs to better than 0.05 min.

I verified the phase-defining claims independently of the subagent, with a
direct query over all 558 frames on the five July nights: the per-night arc
counts (2, 2, 1, 1, 0), the ThAr1/B1/10 s configuration, the 16-per-night B2
`iodin = F` flats against the 1-2 B1 `iodin = T` flats, the near-universal
`hatopen = T`, and the absence of biases all reproduce exactly.

**New files (untracked -- staging is yours).**

- `first_hires_exoplanet/koa_survey.py` -- the survey, re-runnable, in `figs.py`
  house style. Runs **offline** from the saved CSVs by default and only touches
  the network with `--refresh`, so the committed CSVs are the reproducible
  record of what the archive said today.
- `first_hires_exoplanet/data/koa_hd187123_science.csv` (37 rows, 15 KB)
- `first_hires_exoplanet/data/koa_hd187123_nights_all.csv` (3124 rows, **1.25 MB**,
  54 columns) -- this is by far the largest file in the repository. It is
  committable, but you may want to trim columns first; flagging it rather than
  deciding for you.

**What I learned about the repository / archive.**

1. *The KOA TAP service is directly usable from `pypeit14` with nothing but
   `requests`* -- no `pykoa`, `astroquery` or `pyvo` needed, none of which are
   installed. The recipe is recorded in the Report section so we never have to
   rediscover it. Four traps cost real time: the table is `koa_hires` and not
   `koa.koa_hires`; `date_obs` is a *char* column and filtering on it raises
   `ORA-01861`, so filter on `mjd`; the column is `xdangl`, not `xdisangl`; and
   `filename` does not exist (`filehand` / `ofname` do).
2. *Public status is answered by the service's own behaviour, not by a column.*
   KOA appends a proprietary-period clause to every query, so anything returned
   is already public. That is a cleaner answer than reading `propint` and doing
   the arithmetic, and it is worth remembering for any future archive question.
3. *`.gitignore` already anticipates this.* The top-level `data/` exclusion is
   explicitly negated by `!first_hires_exoplanet/data/` (line 234), so metadata
   files placed there are tracked normally while raw frames stay out. The
   existing `hd187123_hires_rv.tsv` set that precedent.
4. *Nights are shared between programs, heavily.* Each July night carries
   T. Bida's Mercury program in a completely different instrument setup, plus
   68-111 frames of other planet-search targets in the *same* configuration as
   ours. Prompt 5 must filter the download rather than pulling whole nights, or
   `pypeit_setup` will be handed hundreds of irrelevant frames.
5. *KOA's `imtype` column does not mean what its description says* for these
   1998 headers -- it holds the amplifier mode (`ONEAMP` / `TWOAMPTOP`).
   `obstype` and `ccdspeed` are empty and `xdname` is `undefined`. Use
   `koaimtyp` plus the lamp columns for frame typing. The boolean columns
   (`iodin`, `hatopen`, `lampcat1`, ...) are `'T'`/`'F'` strings.
6. *The survey surfaced three things in the PypeIt source for prompt 4 to chase*,
   all recorded as risks in the Report: the `not HATOPEN` requirement for flat
   typing, the binning-axis inconsistency between `specaxis = 1` and the
   inherited `compound_meta('binning')`, and `config_specific_par` scaling
   `order_spat_range` by a mosaic-sized 6200 px on a 2048 px chip.

**No FITS data was downloaded, no PypeIt source was edited, and no git command
changed state.**

### 2026-09-20 (Prompt 4: read `keck_hires_orig`, compared it against the archive, listed the gaps)

**Task.** Read `keck_hires_orig` in the PypeIt source, report what it expects
(header cards, configuration keys, detector parameters, frame-typing rules,
required calibrations, default `PypeItPar`), compare that against the prompt-3
archive survey, and list every gap. No editing. Delegated to a Fable subagent as
the prompt asked; the load-bearing claims were then re-verified here against the
source and against the saved KOA CSVs. Findings are in the `## Report` section.
**Nothing was edited in either repository.**

**The one result that changes the plan:** `par['scienceframe']['exprng']` is
`[601, None]`, and every HD 187123 frame in the July run is 215-427 s. So **the
science frames will not be typed as science by the automatic pass** — they will
be commented out of the `.pypeit` file. This is now the most certain obstacle
for prompts 6-7, and it is a one-line fix in the data block rather than anything
deep. I verified the parameter and the `check_frame_type` logic directly in the
source.

**The best new finding:** the flat problem has a clean solution one night
earlier. 1998-07-14 (D. Latham, N10H) carries **four iodine-free, hatch-closed,
decker-B1 3 s "Narrowflat" frames** — `HI.19980714.07711`, `.07810`, `.07910`,
`.55408` — with `echangl` 0.00106-0.00114 and `xdangl` -0.5420, inside the angle
tolerances of the entire July 15-19 run, and typed `IntFlat` automatically. That
night also has 8 hatch-closed B1 ThAr arcs. Prompt 3 had concluded the only B1
flats were the iodine ones; that is true *within* the run, but not for the night
before it. I verified this against the saved CSV rather than taking it on trust.
**Prompt 5 should download the July 14 flats along with whichever night it
picks.**

**What I learned about the repository / PypeIt.**

1. *Inheritance matters more than it looks.* `KeckHIRESOrigSpectrograph` derives
   from `KECKHIRESBaseSpectrograph`, **not** from the post-2004
   `KECKHIRESSpectrograph`. Confirmed at runtime from the MRO. This voids half of
   prompt-3's risk 7 — `config_specific_par`, with its mosaic-sized 6200 px
   `order_spat_range` and its `fwhm = 8.0/bin_spec`, is simply unreachable for
   our class, and `order_spat_range` correctly defaults to `[0, 1024]`. A good
   reminder to check reachability before recording a risk.
2. *Frame typing hangs entirely on instrument logicals, not on an image-type
   card.* The commented-out `IMAGETYP` line in the Orig `init_meta` is dead code;
   `idname` is computed from `HATOPEN`, `AUTOSHUT`, `XCOVOPEN`, `RCCVOPEN`,
   `LAMPCAT1/2`, `LAMPQTZ2` and `LAMPNAME`. Two consequences worth remembering:
   flats require the hatch **closed** (arcs do not care), and biases require the
   covers **open**, which is why the July-14 biases will also come out untyped.
3. *Bias and dark are not required calibrations here.* `use_biasimage` and
   `use_darkimage` are False for every frame type, with overscan used instead.
   Prompt-3's "no biases on the run nights" risk is therefore mostly void — a
   useful reduction in scope for prompt 5's download.
4. *The echelle angle files are already on disk and git-tracked*
   (`keck_hires_orig_angle_fits.fits`, `keck_hires_composite_arc.fits` in
   `pypeit/data/arc_lines/reid_arxiv/`). I had expected a
   `pypeit_cache_github_data` step to be needed before any reduction; it is not.
5. *Three things are genuinely undecidable without a real frame*, and prompt 5
   should print them before anything else: whether the header has a **`PREPIX`**
   card (the Orig reader uses `PREPIX`, while KOA only exposes `PRECOL = 21`, and
   its absence is a hard `KeyError`); whether the logicals are FITS booleans
   rather than `'T'`/`'F'` strings (as strings, science frames would mistype as
   arcs and the gain lookup would raise); and whether `MJD` or `UTC` is present.
6. *The binning inversion is real but narrower than feared.* The PypeIt binning
   string comes out `'2,1'` for a `'1,2'` header while `specaxis = 1` makes the
   opposite true physically. But the echelle wavelength path never uses
   `binspectral`, and `bpm` is correct precisely because it inverts the naming
   again. The live consequences are `order_platescale` (2x, so arcsec quantities
   like the 1.5" boxcar radius) and mislabelled output headers. It cannot be
   fixed from a `.pypeit` file; it is an upstream report, and Ryan Cooke is the
   contact. Three smaller upstream nits are banked in the Report.

**No PypeIt source was edited, no files were created, and no git command changed
state.** The only change is this document.

### 2026-09-20 (Prompt 5: downloaded 1998-07-16 and settled every open header question)

**Task.** Download one night of the 1998 July run with its calibrations, outside
the repository, report what arrived, and confirm the headers match prompt 4's
expectations. Delegated to a Fable subagent; the load-bearing findings were
re-verified here directly against the downloaded FITS files. Details are in the
`## Report` section.

**What arrived.** 31 frames, 152.8 MB, in
`/Users/xavier/Projects/PypeIt/first-hires-exoplanet-data/raw/` — outside the
repo, confirmed by `find`. I chose **1998-07-16** because it is the only night of
the run with hatch-closed flats (which PypeIt types automatically) and it has two
ThAr arcs bracketing the night; the four clean B1 flats from 1998-07-14 came
too, per prompt 4's finding.

**Prompt 4's three undecidable questions are now decided.**

- *Header logicals are genuine FITS booleans* (`<class 'bool'>`), so the typing
  logic and the `CCDGAIN` lookup work as written. **Gap 10 PASS** — the worst
  case (science frames mistyping as arcs) does not happen.
- *`MJD` exists* (as a string, which PypeIt casts). **Gap 11 PASS.** Worth
  knowing: `UTC` is absent and `DATE-OBS` is date-only, so the documented
  fallback would in fact raise — it is simply never reached.
- *`PREPIX` exists*, so no `KeyError`. **But checking it turned up a real
  PypeIt bug** — see below.

**The main result: a genuine upstream bug in the original-CCD raw reader.**
`NAXIS1 = 2303 = PREPIX + 2048 + POSTPIX`, not `2*PREPIX + ...`, and the column
medians show the data occupying columns 21:2069. But `get_rawimage`
(keck_hires.py:1050-1055) computes `col0 = prepix * 2`, and at runtime hands back
`datasec cols 42..2089`. **Every frame's data section is offset by exactly 21
columns**, dropping 21 real columns at one end of each order and including 21
overscan columns at the other. The two-amplifier frame explains the slip:
`prepix * 2` is right for `NUMAMPS = 2` and wrong for `NUMAMPS = 1`, so the code
should say `prepix * namp`. This only bites single-amplifier original-HIRES data,
which is exactly our case and is presumably why it survived. I verified the
arithmetic, the column medians and the runtime `datasec` myself rather than
taking the subagent's word.

**What I learned about the data and the repository.**

1. *`specaxis = 1` is correct, and the binning inversion is confirmed on real
   data.* Collapsing a B2 flat along columns gives 25 order peaks along the row
   axis at 29.5 px spacing; collapsing along rows gives a smooth blaze. So
   dispersion runs along the unbinned 2048-column axis and the x2-binned row axis
   is spatial — physically binspec 1, binspat 2, while PypeIt reports `'2,1'`.
   Prompt 3's risk 7 is now fully resolved: half void, half a confirmed
   upstream naming bug with limited consequences.
2. *Frame typing matched prompt 4's predictions 30 out of 30.* The simulation
   against KOA metadata was reliable. The one detail worth recording: the 07:39
   UT arc has the hatch **open** and still types `Line`, confirming on a real file
   that the hatch is not tested for arcs.
3. *The hard-coded MAKEE bad-pixel mask misses the bad columns.* The real bad
   columns sit at raw NAXIS1 106-107, 1148, 2027-2028 — MAKEE's coordinates count
   raw columns *including* the 21 prescan columns, but PypeIt applies them to the
   trimmed frame. They land ~41 px off now, and would still be ~20 px off after
   the reader fix. 2027-2028 is a strong hot pair, so prompt 8 should look for a
   spurious feature there.
4. *Gain and read noise are worth querying.* The header says
   `DETECTOR = 'Tek 2048E LRIS Eng grade device'`, `CCDGN01 = 4.8`,
   `CCDRN01 = 6.0`; PypeIt hard-codes 1.9 e-/ADU and 2.8 e-. The KOA DQA cards
   may be generic rather than measured for 1998, so this is a **question for Ryan
   Cooke**, not a bug claim — but it must be settled before we quote a
   signal-to-noise.
5. *KOA retrieval traps.* `filehand` is the only retrieval key and is **not** in
   the survey CSVs, so it must be re-queried. A wrong `filehand` returns an HTML
   error page with **HTTP 200**, so downloads must be content-verified rather
   than status-checked. And `astropy.io.fits` refuses `memmap=True` on these
   files because of the `BZERO`/`BSCALE` cards.

**Three items now stand for an upstream report to Ryan Cooke**: the
`prepix * 2` reader offset (concrete, one line), the BPM column origin, and the
binning-string inversion — plus the gain/read-noise question.

**New file (untracked):** `first_hires_exoplanet/koa_download.py`. No raw data
is in the repository, no PypeIt source was edited, and no git command changed
state.

### 2026-09-20 (Prompt 6: fixed the `keck_hires_orig` bugs on PypeIt branch `orig-hires-fixes`)

**Task.** Fix the issues identified in prompts 4 and 5, on the new PypeIt branch
`orig-hires-fixes`. Delegated to a Fable subagent; I reviewed the complete diff
and re-ran every verification myself. **One file changed**,
`pypeit/spectrographs/keck_hires.py` (+131/-33). No git command changed state --
the work is in the working tree for you to commit and push.

**Six fixes, all verified against the real 1998 frames.**

1. **`get_rawimage` prescan offset** (the most important). `col0 = prepix * 2`
   is correct only for `NUMAMPS = 2`; our 1-amp frames have
   `NAXIS1 = PREPIX + 2048 + POSTPIX`. Changed to `prepix * namp` in both the
   data and overscan expressions. Verified: the science frame now gives
   `data 21..2068 (n=2048)` and `oscan 2069..2302 (n=234)` instead of
   `42..2089` / `2090..2302`. Note the overscan was previously **truncated to
   213 columns** because it ran off the array -- that is fixed too. The
   2-amplifier frame is byte-for-byte unchanged (`data 42..1065`,
   `oscan 2090..2323`).
2. **BPM column origin.** I measured the real bad columns on the science frame:
   raw 0-indexed NAXIS1 columns **106, 107, 1148, 2027, 2028**, which are
   exactly MAKEE's numbers read as raw 0-indexed columns. Since the trimmed
   frame starts at raw column `PREPIX`, the correct mapping is
   `trimmed = C - PREPIX`, not `C - 1`. `bpm` now reads `PREPIX` from the
   example file and applies that origin, falling back to the old behaviour with
   a warning when no file is available. Verified: the fully-masked spectral
   columns are now **[85, 86, 1127, 2006, 2007]**, precisely the measured bad
   columns minus 21; previously they were [105, 106, 1147, 2026, 2027], masking
   clean pixels while leaving the real defects exposed. The subagent's
   end-to-end check through `RawImage.process` found the newly masked pixels
   carry +515, +171, +40 and +15 e- of excess while the old positions showed
   exactly 0.0 -- good confirmation.
3. **Binning inversion, Orig class only.** Added a
   `KeckHIRESOrigSpectrograph.compound_meta` that overrides **`binning` alone**
   and delegates everything else to the base class, so header `'1,2'` now yields
   PypeIt `'1,2'` (binspec 1, binspat 2) rather than the inverted `'2,1'`. The
   paired change in `bpm` (`xbin, ybin = binspec, binspat`) went in with it --
   that pairing was the risky part, since the old BPM was correct *only* because
   it inverted the buggy string a second time. Verified both: binning `'1,2'`
   and the ink spot still at rows 518-565. `order_platescale` now returns the
   physically correct **0.432**"/binned spatial pixel instead of 0.216.
4. **Missing `raise`** in `check_spectrograph` (L567): post-2004 data handed to
   `keck_hires_orig` silently passed. Now raises. Our 1998 frames still pass.
5. **Science exposure floor**, Orig class only: `[601, None]` -> `[1, None]`.
   **This is the one policy choice rather than a bug fix, so flagging it.** The
   reasoning is sound and I checked it: the floor is not what separates science
   from calibrations -- `check_frame_type` requires `idname == 'Object'`, and
   every calibration has a different `idname`, so no exposure range can turn a
   2 s flat into a science frame. The only overlap is with `standard`, which
   additionally requires an archive-standard positional match (HD 187123 returns
   `False`). Verified over all 30 frames: the 400 s science frame now types as
   **`science`**, and every other frame types exactly as before. If you would
   rather this stayed a `.pypeit`-file setting, revert that one hunk -- nothing
   else depends on it.
6. **Docstrings**: the `config_independent_frames` claim about a `DATE-OBS` rule
   (which the code does not implement) corrected; the `bpm` binning convention
   and column origin documented; `get_rawimage` given a real docstring noting
   `spectrim`, `PRELINE` and `POSTLINE` are unused.

**Deliberately not changed: gain and read noise.** The headers carry
`CCDGN01 = 4.8` and `CCDRN01 = 6.0` against PypeIt's hard-coded 1.9 e-/ADU and
2.8 e-. That is an open **question** for Ryan Cooke, not a diagnosed bug, and
changing it unilaterally would silently rescale every count we later quote.

**Regression checks.** The post-2004 `keck_hires` class is untouched: it still
builds with `ndet = 3`, its `compound_meta('binning')` still inverts as before,
and its resolved parameter tree is byte-identical (5331 lines). The
2-amplifier frame reads identically. `pytest test_spectrographs.py
test_inputfiles.py test_trace.py` gives **51 passed**. No existing test
exercises `keck_hires_orig` at all, which is worth remembering.

**What I learned.**

1. *The two riskiest fixes were coupled, and only correct together.* The BPM was
   right by accident: `xbin, ybin = binspat, binspec` inverted the buggy binning
   string a second time. Fixing `compound_meta` without flipping that line would
   have silently transposed the mask -- a failure with no traceback and no
   obvious symptom. Worth remembering that a "harmless naming artefact" can have
   a compensating bug downstream holding it up.
2. *Fixing the reader changed the BPM's ground truth.* Fix 2 only makes sense
   *after* fix 1, because the trimmed frame's origin moves from raw column 42 to
   raw column 21. The two had to be reasoned about as a pair, not independently.
3. *MAKEE's mask coordinates are raw, prescan-inclusive, and effectively
   0-indexed* relative to the numbers in `MaskHIRES_1x1.dat` -- established by
   measuring the real defects rather than by reading the file's documentation.
   The hot-corner entry at column 2059 is corroborating evidence: it is off the
   end of a 2048-column trimmed detector under the old reading, but sits on real
   excess flux under the new one.
4. *`first_hires_exoplanet/verify_orig_fixes.py`* is the new re-runnable driver
   that produced these numbers; run it before and after any future change to
   this class and diff the output.

**New file (untracked):** `first_hires_exoplanet/verify_orig_fixes.py`. The
PypeIt working tree holds the six fixes on `orig-hires-fixes`, uncommitted, for
you to review and push. **Four items remain for an upstream conversation with
Ryan Cooke**: the four code bugs above (as a PR), the science-floor policy
choice, and the gain/read-noise question.

### 2026-09-20 (Prompt 7: ran `pypeit_setup`, reported the automatic pass)

**Task.** Run `pypeit_setup` on the 1998-07-16 night and report the
configurations, the frame types, and anything untyped or mistyped — without
hand-editing the `.pypeit` file. Delegated to a Fable subagent; I inspected the
generated `.pypeit`, `.calib` and `.sorted` files directly. Full findings are in
the `## Report` section. No `.pypeit` file was edited and no overrides were
passed.

**The headline: the night on its own cannot be reduced.** `pypeit_setup` splits
1998-07-16 into two configurations on `decker`, and the B1 setup that holds the
science frame and the arcs has **no `trace`, `pixelflat` or `illumflat` frame at
all** — the only flats that night are the B2 ones (different decker, hence a
different configuration) and the hatch-open/iodine ones (untyped). Adding the
four 1998-07-14 clean B1 flats fixes it completely: they merge into the science
configuration, and calibration group 0 then has arcs, tilts, all three flat
roles and the science frame. **Prompt 8 needs no hand-typing at all**, which is a
better position than prompts 3-5 suggested we would be in.

**The prompt-6 fixes hold up in a real pipeline run.** The 400 s science frame is
now typed `science` (fix 5) — and only `science`, since the archive-standard
lookup returns False — and the `binning` column reads `1,2` throughout (fix 3).
17 of 30 frames typed, 13 untyped, **0 mistyped**, matching prompt 4's
prediction exactly. No new bugs.

**What I learned about the repository / PypeIt.**

1. *`check_env.py` earned its keep on its first real outing.* Switching PypeIt to
   `orig-hires-fixes` left `pypeit.__version__` stale at `f3a1f1d27`, and the
   script caught it before any reduction product was written. I refreshed the
   install, and every file generated in this prompt is stamped
   `2.0.2.dev1217+g017bece06`. **This needs doing after every branch move** —
   the trap recurs exactly as predicted in the prompt-2 log.
2. *`pypeit_setup` writes either `.pypeit` files or `.sorted`, never both.* With
   `-c` you get the per-configuration `.pypeit`/`.calib`; without it you get
   `setup_files/*.sorted`/`.obslog`. Getting both means running it twice. Easy to
   trip over, and `run_setup.py` now encodes it.
3. *`vet` does not check calibration completeness.* Run A's science `.pypeit`
   file has no flats whatsoever and still reported "PypeIt file successfully
   vetted". A vetted file is not a reducible one — worth remembering in prompt 8,
   where a failure could otherwise look mysterious.
4. *Setup letters are not stable.* B1 is Setup B in Run A and Setup A in Run B;
   the angles recorded in the `.sorted` header are those of the first frame
   encountered, not the science frame. Never refer to a setup by letter across
   runs.
5. *The cross-night merge is physically sound, not just within tolerance.*
   Cross-correlating order profiles between 07-14 and 07-16 gives 0 to -1 binned
   spatial pixels of shift — the same as between two 07-16 flats taken 7.5 h
   apart — and 0 px in dispersion. So borrowing the 07-14 flats is not merely
   permitted by `configuration_keys`, it is justified by the data.

**New file (untracked):** `first_hires_exoplanet/run_setup.py`, which encodes the
two-pass invocation, refuses to write inside the repository, and reproduces both
runs (verified byte-identical apart from the UTC stamp). Two small diagnostic
scripts live in the data tree at `redux/checks/`; bring them into the repo if you
want them committed. No data artefact of any kind is inside the repository, and
no git command changed state.

### 2026-09-20 (Prompt 8: reduced 1998-07-16 — 37 orders, S/N 33-164)

**Task.** Reduce the night with `run_pypeit`, diagnose any failure with the
`diagnose-reduction` skill, fix what belongs in the input file, and separate
genuine PypeIt bugs from configuration mistakes. Results are in the `## Report`
section. **The reduction works**: 37 orders, S/N 33-164, all wavelength
solutions under 0.2 px RMS.

**An honest account of how this went.** I attempted to delegate this to a Fable
subagent; the launch was rejected, and I said I would do it myself. In fact the
subagent had already started, and by the time I looked at the directory it held
three runs and three edited `.pypeit` files. I misread that as activity from an
unknown source and said so to the user, when checking the timestamps against my
own actions (staging at 06:14:43, runs beginning 06:16:10) would have identified
it immediately. Two lessons: **reconcile unexpected state against my own recent
actions before attributing it elsewhere**, and note that **a rejected tool call
may still have executed** — do not assume the rejection means nothing ran.

That work was not wasted, but it could not be adopted as it stood. Its `.pypeit`
comment asserted that the baseline failed because "every order ended flagged
BADSKYSUB"; that string appears **zero** times in the run log. The actual symptom
was `No objects found automatically`. Since these comments are destined for a
public reduction recipe, I moved the directory aside as
`reduce_jul16_prior_unverified/` and rebuilt from the pristine prompt-7 file,
adding one parameter at a time. The conclusion held — but one of the three
parameters that run used (`min_frac_prof = 0.5`) turned out to be **unnecessary**,
and I would have carried it forward as folklore had I simply adopted the file.

**What I learned.**

1. *The real mechanism is geometric, and it is a consequence of our own fix.*
   The B1 orders are 10.32 binned pixels wide; `find_trim_edge = [3, 3]` leaves
   4.3 px, and object finding drowns. It only became visible now because prompt 6
   corrected the binning inversion so `order_platescale` returns the true
   0.432"/binned pixel. A fix in one place surfaced a latent default problem in
   another.
2. *Incremental verification paid for itself immediately.* `find_trim_edge` alone
   is necessary but **not** sufficient (25 of 37 orders still failed). Had I
   applied all three parameters at once I would have learned nothing about which
   mattered, and would have shipped a fourth that does nothing.
3. *`run_pypeit` exits 0 with no 1D spectrum.* The most costly property of this
   whole task: a silent failure that looks like success. Banked as an upstream
   item.
4. *The `diagnose-reduction` skill's warning about stale cached calibrations is
   the real thing to remember here.* `-o` overwrites science products but reuses
   `Calibrations/`. That is legitimate when only `reduce` parameters change (as
   in v1-v3, which is why re-runs took 40 s rather than 5 min), and a trap the
   moment frame typing or calibration grouping changes.
5. *The star fills the slit.* Both sky-related parameters trace to one physical
   fact: a V=7.9 star in a 3.5" slit at 400 s leaves no sky pixels. This will be
   true of **every** frame in this programme, so the same three parameters should
   carry to the other nights in prompt 10.

**Upstream tally is now seven items** for Ryan Cooke: the four code bugs from
prompt 6, the science-floor policy choice, the gain/read-noise question, and now
the silent-success exit code, the `find_trim_edge` default at 2x binning, and the
appending `spec1d_*.txt`.

No repository file changed; all products are in the data tree. No git command
changed state.

### 2026-09-21 (Prompt 9: inspected the reduction — it is sound)

**Task.** Check order tracing, flat field and wavelength solution against the QA
products; examine the extracted 1D spectrum for the iodine forest, for G2V
stellar lines, and for signal-to-noise; write a figure script to disk. Findings
are in the `## Report` section. **The reduction passes every check.**

Headline numbers: 37 orders traced with none masked; pixel flat good to +/-3%
over the illuminated pixels; all 37 wavelength solutions under 0.2 px RMS
(median 0.133); continuous coverage 3806-6262 A; S/N 33 rising to 168 (median
124); the I2 forest doubling the line density between 5000 and 6200 A; and every
strong G2V feature in range at its expected wavelength.

**The strongest single result is the velocity zero point.** Seven photospheric
lines give a mean **-14.7 km/s** against HD 187123's systemic **-17 km/s**. That
is an *external* check — it tests the wavelength solution against the sky rather
than against its own arc — and it is the first evidence in this phase that the
reduction is not merely self-consistent but correct.

**A mistake I made and caught, worth recording because it would have derailed
phase 2.** My first pass compared observed wavelengths against *air* rest
values and found a consistent **+70 km/s** offset across every line. I very
nearly reported it as a wavelength zero-point defect. The cause was entirely
mine: **PypeIt reports vacuum wavelengths**, and air-to-vacuum at 5000 A is
+1.4 A, almost exactly the offset I measured. The tell was that the "error" was
suspiciously constant in velocity and matched the refractive index of air. Two
lessons: an anomaly that is constant in *velocity* across a wide wavelength
range is far more likely to be a reference-frame mistake than an instrumental
one, and a 70 km/s systematic in data whose whole purpose is 72 m/s precision
deserves several minutes of scepticism before it is written down. The trap is
now recorded in the figure script's docstring.

**A second self-correction.** I first reported the pixel flat as "1-99% =
0.986-1.003", which is wrong in substance: pixels outside the traced orders are
left at exactly 1.0 and dominate the histogram. Restricted to the 0.1M
illuminated pixels the true spread is **0.959-1.024**. The flat is still
excellent, but the number I would have published was meaningless.

**What I learned.**

1. *`run_pypeit -o` accumulates into `spec1d` rather than replacing it.* The
   file held **74 entries for 37 orders** — the v2 and v3 passes superimposed,
   with the *bad* v2 extraction first in the list. Anyone iterating on
   parameters and then analysing the output would silently mix reductions. I
   cleared `Science/` and re-ran to get an unambiguous product, and
   `load_orders()` now refuses to plot a directory in that state. This upgrades
   the prompt-8 note about the `.txt` summary: the FITS product does it too, and
   that is a genuine data-integrity issue for upstream.
2. *Bad columns arrive as exact zeros, not as flagged gaps*, so naive continuum
   normalisation turns them into full-depth spurious "lines". Mask
   `OPT_COUNTS <= 0`.
3. *Na D is not covered*, falling in an inter-order gap — the free spectral range
   stops overlapping at the red end of this setup. Expected, but worth knowing
   before someone goes looking for it in phase 2.
4. *The iodine band is a clean, quantitative diagnostic.* Line density per order
   is a blunt statistic, but the factor-of-two step at 5000 A against a smoothly
   falling stellar trend is unambiguous and needs no template.

**New file (untracked):** `first_hires_exoplanet/figs_phase1.py`, plus four PNGs
in `docs/figs/`. No raw data or reduction product is in the repository; no git
command changed state.
