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

6. Read this file. Run `pypeit_setup` on that night. Report the configurations
   it identifies, the frame types it assigns, and any frame it fails to type or
   assigns wrongly. Do not hand-edit the `.pypeit` file yet — first report what
   the automatic pass produced. Use Fable if you can. Log your work.

7. Read this file. Reduce that night with `run_pypeit`. If it fails, use the
   `diagnose-reduction` skill to work out why, fix what can be fixed in the
   input file or parameters, and report anything that looks like a genuine
   PypeIt bug rather than a configuration mistake. Use Fable if you can. Log your work.

8. Read this file. Inspect the reduction. Check the order tracing, the flat
   field, and the wavelength solution against the QA products, and look at the
   extracted 1D spectrum: are the iodine lines present between roughly 500 and
   620 nm, are the stellar lines where they should be for a G2V star, and what
   signal-to-noise did we get? Write a figure-generating script to disk, as with
   the public document. Use Fable if you can. Log your work.

9. Read this file. Reduce the remaining nights of the 1998 July run the same
   way, and report whether the reduction is stable from night to night. Use
   Fable if you can. Log your work.

10. Read this file. Write a short assessment of phase 1: whether the reduction is
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
