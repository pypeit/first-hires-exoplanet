# DRAFT — not sent

**To:** Paul Heeren (and René Tronsgaard, Frank Grundahl), `pyodine`
**From:** J. Xavier Prochaska
**Subject:** `pyodine` on 1998 Keck/HIRES data: what we changed, and what we found

---

Dear Paul,

We have used `pyodine` (github.com/pepeheeren/pyodine, commit `4488b09`) to
measure velocities from the original 1998 Keck/HIRES iodine-cell spectra
behind the discovery of HD 187123 b. The spectra were reduced with PypeIt,
and the atlas is the Fischer May 2022 FTS scan, which is a different cell
from the HIRES one. It worked, and we recover the planet at 10σ
(K = 72.4 ± 7.1 m/s, 26 m/s per epoch). Thank you for making the code open;
this would not have been possible otherwise.

We vendored it as a fork (`vendor/pyodine/` in <repository link>). Every
change carries a test that fails before and passes after. Below are the
changes we think are generally correct and not HIRES-specific, then the
things we found and worked around without changing your code. Every item
has file and line references; the tests reproduce the first group.

## A. Changes in our fork that we think belong upstream

| # | change | files | test |
|---|---|---|---|
| 0 | NumPy 2: `np.float` (3 sites) and `np.NaN` (5 sites) no longer exist | `lib/misc.py`, `fitters/lmfit_wrapper.py` | `test_fork_numpy2.py` |
| 1 | **Iodine depth as Beer–Lambert**, T^α, instead of the linear scaling at `models/spectrum.py:93`. At α ≈ 2.6 (the HIRES cell against this atlas) the linear form makes 5.8% of the atlas transmission negative | `models/spectrum.py` | `test_fork_iodine_depth.py` |
| 2 | **Atlases record their wavelength frame** (`IodineAtlas.wave_frame`). Your loader reads the air grid; a vacuum-calibrated reduction is then 83 km/s out, and nothing says so | `components.py` | `test_fork_atlas_frame.py` |
| 3 | **`bary_date` and `bary_vel_corr` docstrings.** `timeseries/bary_vel_corr.py` passes `bary_date` straight to `barycorrpy` as a *full JD(UTC)*, not the reduced BJD the docstring states, and `chunks.py` treats `bary_vel_corr` as m/s, not km/s. We corrected the docs and added a guard | `components.py` | `test_fork_bary_units.py` |
| 4 | **`compute_weight(weight_type='ivar')`** and `Spectrum.ivar`, to weight by a propagated inverse variance | `components.py` | `test_fork_weights.py` |
| 5 | **The template-index trap.** The model and the plots fetch template chunk *i* for observation chunk *i*, but `auto_wave_comoving` builds chunks only for the orders passed. Any subset of orders silently pairs every chunk after the gap with the wrong template chunk. Chunks now record `template_index`, and `models.spectrum.template_chunk_index` resolves it | `chunks.py`, `models/spectrum.py`, `plot_lib.py` | `test_fork_chunk_index.py` |

## B. Found and worked around; no change made to `pyodine`

Correctness:

6. **`order_correction` is one constant shift.** The driver maps template
   order *o* to observation position *o* + `order_correction`, which assumes
   the observation's orders are contiguous. One of our nights was extracted
   without one order, and every template order past the gap was fitted
   against the neighbouring order. That gave velocities of up to 7.7 × 10⁶
   m/s and no error. Matching by physical order number would fix it; we
   filter the orders before calling the driver.
7. **The Chauvenet re-fit never runs.** `pipe_lib.model_all_chunks` re-fits
   only `if any(mask)==False`, i.e. when *every* pixel fails the criterion;
   presumably "if any pixel fails" was meant. The criterion also runs on
   unweighted residuals, so zero-weight pixels are counted: on one chunk,
   5 of its 6 flags were masked pixels. If the condition were repaired, it
   would reject the wrong pixels.
8. **`red_chi_sq` holds √χ²_ν**, not χ²_ν (`model_all_chunks`).
9. **`jansson` uses the wrong adjoint.** The update is `x += r(Ky − KKx)`
   where least squares needs `Kᵀ`. The fixed point is the same, and for a
   symmetric kernel so is everything else. For a skewed kernel we measure a
   −5.7 m/s bias per chunk against a true-adjoint copy.
10. **The template spline breaks re-convolution at HIRES sampling.**
    `ChunkedDeconvolver` puts the observation on its 10×-oversampled grid
    with `misc.rebin`, a cubic spline through the detector pixels. With lines
    2.2 px wide it overshoots near deep lines by 6–18σ at S/N ~800. Jansson
    fits the wiggles, and the re-convolved template misses the data
    (χ²_ν 5.2 against 1.1 with no oversampling). `osample_temp = 1` avoids
    it. For well-sampled SONG data it will not show.
11. **`fit_lsfs` starts shape parameters at ~10⁻¹².** Lick's run-1 scheme
    starts from `fit_lsfs`, which sends every parameter a Gaussian does not
    need (satellite amplitudes, Hermite weights) to ~10⁻¹². `lmfit`'s
    `leastsq` then takes finite-difference steps *relative to the value*,
    the derivatives are zero, and the fit returns its start after one
    Jacobian. This is probably why SONG's Hermite configuration pins its
    weights. We start each shape parameter at least 10⁻³ from zero and 2%
    inside its bounds.
12. **`NormalizedObservation` drops `ivar`.** Its `__getitem__` rebuilds each
    `Spectrum` without it.
13. **`create_template` applies the template's bad-pixel mask to the hot
    star's weights**, not a mask computed on the hot-star frame.
14. **The Arcturus reference is on air wavelengths.** The velocity guess
    against a vacuum-calibrated observation comes out at ~+75 km/s for a
    −10 km/s star. It is harmless in differences, but the value must not be
    read as a velocity.
15. **The model normalises the atlas by its chunk mean before
    exponentiating**, so `iod_depth` and `cont_intercept` are correlated by
    construction. It is not wrong, but it is worth documenting.

Robustness and usability:

16. **Exceptions are swallowed.** `create_template`
    (`pyodine_create_templates.py:620`) and `model_single_observation`
    (`pyodine_model_observations.py:521`) wrap their whole body in
    `except Exception` and return normally. The first logs with
    `logging.error`. The second prints the traceback and logs with
    `logging.info`, so its error log stays empty. Either way, a failed stage
    is indistinguishable from a successful one to the caller.
    `normalize_single` likewise swallows the solar-atlas `ValueError` and
    falls back to `top` silently: in 13 of 14 of our orders.
17. **Plot, then fork, and the pool hangs forever on macOS.** Once the
    template stage has plotted with the `macosx` backend, `pathos` forks
    workers into a process holding CoreFoundation. They finish and write
    correct results, and the pool never joins. `MPLBACKEND=Agg` avoids it;
    setting a non-interactive backend inside `pyodine` would too.
18. **`plot_chunks` indexes the chunk list by fixed numbers**
    (`[150, 250, 400]` in the Lick parameters). An observation with fewer
    chunks raises `IndexError` in `create_analysis_plots` after the results
    are written.
19. **`timeseries/misc.robust_mean` divides by zero** (`lib/robust.py:155`,
    `sqrt(len(good) − 1)`) on a column with one finite value.
20. **`combine_chunk_velocities` mutates its module-level default**
    `_weighting_pars` (`del pars['good_orders']`) when called without
    `weighting_pars`, so a second call raises `KeyError`.

## C. One thought on the combination

`combine_chunk_velocities` weights chunks by their own time-series scatter,
as if independent. On our data whole orders shift coherently: the weak-iodine
red orders scatter by 100–200 m/s as orders. Per-chunk weights do not see
that. Dropping those orders halves our residual, and an order-level weight
would do the same without hand selection. It is not a bug, but it may
interest you.

We are happy to send any of the fork changes as pull requests, one per
change with its test. Would that be welcome, and against which repository
(GitHub or GitLab)?

Best wishes,
Xavier Prochaska
