""" Deconvolve the HD 187123 template against the instrumental profile
(phase 3, prompt 3).

`pyodine` models an iodine-cell observation as

    (deconvolved stellar template  x  iodine atlas)  convolved with the LSF

so the template it wants is the star *before* the spectrograph blurred it.
Phase 2 built the co-added cell-out observation
(`redux/template/hd187123_template_19980826.fits`, S/N 310); this module runs
it through `pyodine`'s own deconvolver, `template/deconvolve.py`
(`ChunkedDeconvolver` -> Jansson), and measures three things:

  1. The S/N cost.  Deconvolution amplifies noise.  Measured by Monte Carlo
     (the Jansson iteration is non-linear, so no closed form applies):
     per-pixel S/N of the deconvolved template, of the template once it is
     re-convolved as the forward model will re-convolve it, and the velocity
     noise the template itself injects into every chunk.
  2. The dependence on the assumed LSF.  Several LSFs, symmetric and not,
     each judged by how well its template reproduces the observation and how
     far it moves the chunk velocities.
  3. Whether the deconvolved template, re-convolved with the same LSF,
     reproduces the observed one.  Chi-square per pixel, against both the
     propagated inverse variance and the noise actually measured from the
     differences between the three co-added exposures (:func:`noise_check`),
     which is 0.43 of the propagated value.

It answered a fourth question it was not asked.  At `pyodine`'s default
template oversampling of 10, the answer to (3) is no, and the reason is not
the LSF: the deconvolver splines the observation onto its fine grid
(`misc.rebin`), HIRES lines are barely Nyquist-sampled at 2.2 pixels FWHM,
and the spline's errors between pixel centres are many sigma at this S/N.
At an oversampling of 1 there is no spline, and the same deconvolver
reproduces the observation to the noise.  Both templates are written.

WHERE THE LSF COMES FROM.  `pyodine` normally takes the LSF from a fitted
hot-star-plus-iodine observation.  That fit needs the model configuration
prompt 4 writes and the LSF model prompt 6 chooses, so here the LSF is
*assumed*, and is passed through the deconvolver's own `lsf_fixed` argument.
The fiducial width is measured from the ThAr lines of the same night's
reduction, which share no machinery with the stellar spectrum.  Phase 2's
R = 42,000 is carried as the alternative: it came from the narrowest *stellar*
lines, which carry the star's own broadening, so it is an upper bound on the
instrumental width, not a measurement of it.

WHAT IS NOT CHANGED.  Nothing in `vendor/pyodine/`.  The one piece of code
here that re-implements part of the fork, :func:`jansson_true_adjoint`, exists
only to test the fork's `jansson` against an asymmetric kernel; it is not used
to build the template.

Reduction: PypeIt 2.0.2.dev1217+g017bece06.  Fork: `pyodine`
4488b0914fe5b272b787982647691045bff2604a plus the phase-3 prompt-2 changes
(see `vendor/README.md`).  Deconvolution and chunking parameters are
`utilities_lick/pyodine_parameters.py`'s `Template_Parameters` at that
commit, copied below rather than imported so that prompt 4's HIRES parameter
file cannot change them underneath this result.

Run with:

    conda run --no-capture-output -n pypeit14 \\
        python -m first_hires_exoplanet.deconvolve_template

"""

# Must precede any import of matplotlib, including pyodine's own: plotting
# followed by a fork hangs on macOS (phase 3 prompt 1).
import os
os.environ.setdefault('MPLBACKEND', 'Agg')

# Standard imports
import time
import argparse
import warnings

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table
from scipy.optimize import curve_fit

# Puts the vendored fork on sys.path, so `pyodine` below is ours
from first_hires_exoplanet.utilities_hires import load_pyodine as lp

from pyodine import components
from pyodine import chunks as pchunks
from pyodine.models import lsf as plsf
from pyodine.models import wave as pwave
from pyodine.models import cont as pcont
from pyodine.models.base import ParameterSet
from pyodine.models.spectrum import SimpleModel
from pyodine.template import deconvolve as pdeconv
from pyodine.template import normalize as pnorm
from pyodine.template.base import StellarTemplate_Chunked


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(_HERE, 'data')
FIG_DIR = os.path.join(os.path.dirname(_HERE), 'docs', 'figs')

TEMPLATE_FITS = os.path.join(lp.DEFAULT_REDUX, 'template',
                             'hd187123_template_19980826.fits')
WAVECALIB = os.path.join(lp.DEFAULT_REDUX, 'reduce_19980826', 'Calibrations',
                         'WaveCalib_A_0_DET01.fits')
#: The deconvolved templates: Lick's settings (oversampling 10), and the
#: same with no oversampling, which prompt 3 found is the only setting that
#: reproduces the observation at deep lines
TEMPLATE_H5 = {os_: os.path.join(lp.DEFAULT_REDUX, 'template',
                                 'hd187123_template_19980826_deconv_os{:d}.h5'.format(os_))
               for os_ in (10, 1)}

CKMS = 299792.458

#: The band the iodine cell absorbs in, where the forward model will work
I2_LO, I2_HI = 5000., 6200.

#: Phase 2's resolving power, from the narrowest weak stellar lines
#: (`build_template.weak_line_widths`, 10th percentile FWHM 7.45 km/s)
R_STELLAR = 42000.

#: `utilities_lick/pyodine_parameters.py`, `Template_Parameters`, at the
#: vendored commit
CHUNK_WIDTH = 40
CHUNK_PADDING = 12
DECONV_PARS = {'osample_temp': 10,
               'jansson_niter': 1200,
               'jansson_zerolevel': 0.00,
               'jansson_contlevel': 1.02,
               'jansson_conver': 0.2,
               'jansson_chi_change': 1e-6,
               'lsf_conv_width': 6.}
REF_SPECTRUM = 'arcturus'
DELTA_V, MAXLAG = 1000., 500

#: Orders used for the Monte Carlo and the LSF comparison: the two ends and
#: the middle of the iodine band.  The full template uses every iodine order.
SUBSET_ORDERS = (58, 63, 69)


# ---------------------------------------------------------------------------
# 1. The instrumental width, from ThAr
# ---------------------------------------------------------------------------

def _gauss_lin(x, amp, mu, sig, c0, c1):
    return amp * np.exp(-0.5 * ((x - mu) / sig) ** 2) + c0 + c1 * (x - mu)


def arc_line_widths(wavecalib_file=WAVECALIB, half=6, min_sep=6.):
    """ FWHM of the identified ThAr lines, per line, in pixels and km/s.

    Each line PypeIt identified is fitted with a Gaussian on a linear
    background over +/- ``half`` pixels of the extracted arc spectrum.  A line
    is kept only if no other identified line lies within ``min_sep`` pixels
    (blends widen the fit), the fit converged with its centre within one
    pixel of PypeIt's, the width is between 0.4 and 3 pixels, and the
    residual is below 5% of the amplitude.

    The arc fills the 0.574 arcsec slit, so its width is the slit-limited
    profile.  A star that under-fills the slit in good seeing is narrower
    still, so this is an upper bound on the stellar LSF, but a much tighter
    one than stellar lines give.

    Args:
        wavecalib_file (str): PypeIt WaveCalib file.
        half (int): half-width of the fitting window, pixels.
        min_sep (float): isolation criterion, pixels.

    Returns:
        `astropy.table.Table`_: one row per accepted line.
    """
    from pypeit import wavecalib
    wv = wavecalib.WaveCalib.from_file(wavecalib_file, chk_version=False)
    rows = []
    for fit in wv.wv_fits:
        if fit is None or fit.pixel_fit is None or fit.spec is None:
            continue
        spec = np.asarray(fit.spec, dtype=float)
        wave = np.asarray(fit.wave_soln, dtype=float)
        disp = np.gradient(wave)
        x = np.arange(len(spec), dtype=float)
        centres = np.sort(np.asarray(fit.pixel_fit, dtype=float))
        for p in centres:
            others = centres[centres != p]
            if len(others) and np.min(np.abs(others - p)) < min_sep:
                continue
            i0, i1 = int(round(p)) - half, int(round(p)) + half + 1
            if i0 < 0 or i1 > len(spec):
                continue
            xx, yy = x[i0:i1], spec[i0:i1]
            amp0 = yy.max() - np.median(yy)
            if amp0 <= 0:
                continue
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    popt, _ = curve_fit(_gauss_lin, xx, yy,
                                        p0=[amp0, p, 1., np.median(yy), 0.],
                                        maxfev=4000)
            except (RuntimeError, ValueError):
                continue
            amp, mu, sig = popt[0], popt[1], abs(popt[2])
            resid = yy - _gauss_lin(xx, *popt)
            if (amp <= 0 or abs(mu - p) > 1. or not 0.4 < sig < 3.
                    or np.std(resid) > 0.05 * amp):
                continue
            k = int(round(mu))
            fwhm_px = 2.354820045 * sig
            rows.append(dict(order=int(fit.ech_order), pixel=float(mu),
                             wave=float(np.interp(mu, x, wave)),
                             amp=float(amp), fwhm_px=float(fwhm_px),
                             fwhm_kms=float(fwhm_px * disp[k] / wave[k] * CKMS),
                             pypeit_fwhm_px=float(fit.fwhm)))
    return Table(rows)


def order_fwhm(arc, orders):
    """ Median arc FWHM (pixels) of each echelle order, for the orders given.

    Orders with fewer than five clean lines fall back to the median over all
    orders, and are flagged.
    """
    allmed = float(np.median(arc['fwhm_px']))
    out = {}
    for o in orders:
        sel = arc['order'] == o
        out[o] = (float(np.median(arc['fwhm_px'][sel])), False) \
            if sel.sum() >= 5 else (allmed, True)
    return out


# ---------------------------------------------------------------------------
# 2. The co-added template, as a pyodine Observation
# ---------------------------------------------------------------------------

#: Moved to the adapter in prompt 4, which needed the same fix for every
#: spec1d; kept under this name for prompt 3's test and callers.
fill_wavelengths = lp.fill_wavelengths


class TemplateObservation(lp._ObsBase):
    """ Phase 2's co-added cell-out template, in `pyodine`'s Observation shape.

    Laid out like `utilities_hires.load_pyodine.ObservationWrapper`: arrays in
    ascending echelle order, indexed by *position*, the physical numbers in
    ``ech_orders``.  Wavelengths are already in the observed frame (the
    file's ``REFFRAME``), vacuum.
    """

    def __init__(self, filename=TEMPLATE_FITS):
        with fits.open(filename) as hdul:
            header = hdul[0].header
            rows = []
            for hdu in hdul[1:]:
                order = int(hdu.name.replace('ORDER', ''))
                d = hdu.data
                rows.append((order, np.asarray(d['WAVE'], dtype=float),
                             np.asarray(d['FLUX'], dtype=float),
                             np.asarray(d['IVAR'], dtype=float),
                             np.asarray(d['NCOADD'], dtype=int)))
        if header.get('REFFRAME') != 'observed':
            raise ValueError('Template is not in the observed frame')
        rows.sort(key=lambda r: r[0])
        self.ech_orders = np.array([r[0] for r in rows], dtype=int)
        self._wave = np.array([fill_wavelengths(r[1]) for r in rows])
        self._flux = np.array([r[2] for r in rows])
        self._ivar = np.array([r[3] for r in rows])
        self.ncoadd = np.array([r[4] for r in rows])
        self.nord, self.npix = self._flux.shape

        good = (self._ivar > 0) & (self.ncoadd > 0) & (self._flux > 0)
        self._weight = np.where(good, self._ivar, 0.)
        self._cont = np.ones_like(self._flux)
        for i in range(self.nord):
            if good[i].sum() > 50:
                f = np.where(good[i], self._flux[i],
                             np.median(self._flux[i][good[i]]))
                self._cont[i] = lp.continuum(f)

        self.orig_header = header
        self.orig_filename = os.path.abspath(filename)
        self.star = lp.Star(str(header.get('OBJECT', 'HD 187123')))
        self.time_start = Time(float(header['MJDREF']), format='mjd',
                               scale='utc')
        # VBARY is km/s in the file; pyodine wants m/s (prompt 2, change 3)
        self.bary_vel_corr = float(header['VBARY']) * 1e3
        self.iodine_in_spectrum = False

    def __len__(self):
        return self.nord

    def __getitem__(self, order):
        if isinstance(order, (int, np.integer)) or hasattr(order, '__int__'):
            i = int(order)
            return components.Spectrum(self._flux[i], wave=self._wave[i],
                                       cont=self._cont[i],
                                       ivar=self._weight[i])
        if isinstance(order, (list, np.ndarray)):
            return [self.__getitem__(int(i)) for i in order]
        if isinstance(order, slice):
            return self.__getitem__([int(i) for i in np.arange(self.nord)[order]])
        raise IndexError(type(order))

    def index_of(self, ech_order):
        return int(np.where(self.ech_orders == int(ech_order))[0][0])

    def iodine_positions(self, lo=I2_LO, hi=I2_HI):
        """ Positions of the orders whose median wavelength is in the band. """
        med = np.median(self._wave, axis=1)
        return [int(i) for i in np.where((med > lo) & (med < hi))[0]]


def normalize(obs, positions):
    """ `pyodine`'s own normalisation, as `create_template` does it.

    Returns:
        tuple: (NormalizedObservation, velocity guess in m/s, reference-atlas
        wavelength offset diagnostics, list of positions whose continuum came
        from the solar comparison rather than the envelope fallback)
    """
    normalizer = pnorm.SimpleNormalizer(reference=REF_SPECTRUM)
    # The velocity guess, from orders outside the iodine band (as Lick does)
    blue = [i for i in range(obs.nord)
            if np.median(obs._wave[i]) < I2_LO and i not in positions]
    v_guess = float(normalizer.guess_velocity(obs[blue], delta_v=DELTA_V,
                                              maxlag=MAXLAG))
    norm = normalizer.normalize_obs(obs, v_guess, orders=positions)

    # Which route did each order's continuum take?  normalize_single swallows
    # the ValueError, so ask compare_sun directly.
    solar = []
    for i in positions:
        spec = obs[i]
        try:
            normalizer.compare_sun(spec.wave, spec.flux / spec.cont, v_guess)
            solar.append(i)
        except ValueError:
            pass
    return norm, v_guess, solar


def norm_sigma(obs, norm, i):
    """ 1-sigma of the normalised flux of order position ``i``; inf if unusable.

    `NormalizedObservation` drops the inverse variance (it rebuilds each
    `Spectrum` without it), so the noise is propagated here by the ratio of
    the normalised to the raw flux.
    """
    raw = obs._flux[i]
    nf = norm._flux[i]
    ivar = obs._weight[i]
    with np.errstate(divide='ignore', invalid='ignore'):
        sig = np.where(ivar > 0, 1. / np.sqrt(ivar), np.inf) * np.abs(nf / raw)
    return np.where(np.isfinite(sig) & (raw > 0), sig, np.inf)


def noise_check(obs, positions, smooth=101):
    """ Is the propagated inverse variance the right size?

    Every chi-square in this module divides by it, and the phase-3 document
    flags it as open: PypeIt reduced these frames with a gain of 1.9 e-/ADU
    and read noise 2.8 e-, against the headers' 4.8 and 6.0.  The template is
    three exposures, so the answer can be measured rather than argued.

    The three share one arc and one wavelength solution and were taken within
    twenty minutes, over which the star moves ~15 m/s on the detector, 0.007
    pixel.  So they are differenced **pixel by pixel**, with no interpolation:
    resampling smooths the noise and would bias the test low.  One exposure
    is scaled to the other's flux level by a running median of their ratio,
    and the difference is divided by its propagated 1-sigma.  If the inverse
    variance is right, that has unit width.

    The width is taken as 1.4826 x MAD, so the cosmic rays a raw difference
    always carries, and the steepest line flanks where 0.007 pixel still
    shows, do not set it.  DER_SNR was tried first and rejected: it assumes
    the spectrum is smooth over five pixels, which lines 2.2 pixels wide are
    not, and it returned "noise" three times the propagated sigma that was
    line curvature.

    Every pair with usable pixels in both members is used: the third
    exposure has no usable pixels in orders 58-60.

    Returns:
        `astropy.table.Table`_: per order and pair, the robust width of the
        normalised difference.
    """
    from scipy.ndimage import median_filter
    frames = [obs.orig_header['FRAME{:d}'.format(k)] for k in range(3)]
    sci = os.path.join(lp.DEFAULT_REDUX, 'reduce_19980826', 'Science')
    exps = []
    for f in frames:
        hit = [p for p in os.listdir(sci) if p.startswith('spec1d_' + f) and p.endswith('.fits')]
        exps.append(lp.ObservationWrapper(os.path.join(sci, hit[0])))
    rows = []
    for i in positions:
        ech = obs.ech_orders[i]
        for a, b in ((0, 1), (1, 2), (0, 2)):
            ea, eb = exps[a], exps[b]
            ja, jb = ea.index_of(ech), eb.index_of(ech)
            f1, iv1 = ea._flux[ja], ea._weight[ja]
            f2, iv2 = eb._flux[jb], eb._weight[jb]
            ok = (iv1 > 0) & (iv2 > 0)
            if ok.sum() < 500:
                continue
            ratio = np.where(ok & (f2 != 0), f1 / np.where(f2 != 0, f2, 1.), np.nan)
            scale = median_filter(np.nan_to_num(ratio, nan=np.nanmedian(ratio)), smooth)
            d = f1 - scale * f2
            with np.errstate(divide='ignore'):
                sd = np.sqrt(1. / iv1 + scale ** 2 / iv2)
            z = (d / sd)[ok]
            mad = 1.4826 * np.median(np.abs(z - np.median(z)))
            # Does one exposure sit shifted on the detector relative to the
            # other?  Same pixels, same wavelength solution, so any shift is
            # the star moving, and the co-add smeared it in.
            w = ea._wave[ja]
            core = ok & (w > 0)
            dv = velocity_shift(f1[core] / scale[core], w[core], f2[core]) \
                if core.sum() > 100 else np.nan
            rows.append(dict(order=int(ech), pair='{:d}{:d}'.format(a, b),
                             npix=int(ok.sum()), width=float(mad),
                             dv_ms=float(dv)))
    return Table(rows)


# ---------------------------------------------------------------------------
# 3. Chunks, wavelength parameters and LSFs
# ---------------------------------------------------------------------------

def build_chunks(obs, positions):
    """ Chunk the template exactly as the Lick parameters chunk an O-star. """
    return pchunks.auto_equal_width(obs, width=CHUNK_WIDTH,
                                    padding=CHUNK_PADDING, orders=positions)


def build_model():
    """ The `SimpleModel` the deconvolver asks for.

    With ``lsf_fixed`` given, `ChunkedDeconvolver` uses only its wavelength
    model; the LSF model and the atlas are needed to construct it.
    """
    return SimpleModel(plsf.SingleGaussian, pwave.LinearWaveModel,
                       pcont.LinearContinuumModel, lp.IodineTemplate(),
                       osample_factor=4, conv_width=DECONV_PARS['lsf_conv_width'])


def wave_params(model, chunks):
    """ Per-chunk wavelength parameters, from the model's own guess.

    `LinearWaveModel.guess_params` is a straight-line fit to the PypeIt
    wavelengths of the chunk.  In a normal `pyodine` run these would come
    from the O-star fit; here PypeIt's solution stands in for it.  The
    residual of that straight line against PypeIt's polynomial, across the
    *padded* chunk the deconvolver evaluates, is returned so its size is on
    record.

    Returns:
        tuple: (list of `ParameterSet`, ndarray of max residual in m/s)
    """
    params, resid = [], []
    for ch in chunks:
        p = model.guess_params(ch)
        params.append(p)
        pad = ch.padded
        w_lin = p['wave_intercept'] + p['wave_slope'] * (
            pad.abspix - ch.abspix[0] + ch.pix[0])
        resid.append(np.max(np.abs(w_lin / pad.wave - 1.)) * CKMS * 1e3)
    return params, np.array(resid)


def deconv_pars(osample=None):
    """ The Lick deconvolution parameters, optionally at another oversampling. """
    pars = dict(DECONV_PARS)
    if osample is not None:
        pars['osample_temp'] = osample
    return pars


def lsf_grid(osample=None):
    """ The pixel vector `ChunkedDeconvolver` evaluates the LSF on. """
    pars = deconv_pars(osample)
    nlsf = pars['lsf_conv_width']
    nfine = int(nlsf * pars['osample_temp'])
    return np.linspace(-int(nlsf), int(nlsf), 2 * nfine + 1)


def gaussian_lsf(fwhm_px, osample=None):
    """ `pyodine`'s SingleGaussian on the deconvolver's grid (unit sum). """
    return plsf.SingleGaussian.eval(lsf_grid(osample), ParameterSet(fwhm=fwhm_px))


def skewed_lsf(fwhm_px, osample=None, frac=0.25, offset=1.2):
    """ An asymmetric LSF: a Gaussian plus a satellite on its red side.

    The main component has the given FWHM; the satellite has the same width,
    ``frac`` of its height, ``offset`` sigma to the right.  The sum is
    re-centred so its centroid is at zero: what is left is pure skewness,
    which is what an asymmetric HIRES profile looks like to a deconvolver.
    """
    x = lsf_grid(osample)
    sig = fwhm_px / 2.354820045

    def prof(xx):
        return np.exp(-0.5 * (xx / sig) ** 2) + \
            frac * np.exp(-0.5 * ((xx - offset * sig) / sig) ** 2)
    y = prof(x)
    c0 = np.sum(x * y) / np.sum(y)
    y = prof(x + c0)
    return y / y.sum()


def lsf_array(chunks, fwhm_of_chunk, kind='gauss', osample=None):
    """ One LSF per chunk, stacked as ``lsf_fixed`` expects. """
    f = gaussian_lsf if kind == 'gauss' else skewed_lsf
    return np.array([f(fwhm_of_chunk(ch), osample) for ch in chunks])


def kms_per_pixel(ch):
    """ Velocity width of one pixel at the centre of a chunk. """
    return float(np.median(np.abs(np.diff(ch.wave)) / ch.wave[1:]) * CKMS)


# ---------------------------------------------------------------------------
# 4. Measuring a deconvolved chunk
# ---------------------------------------------------------------------------

def fine_indices(tchunk, pixels):
    """ Indices of the fine grid that sit exactly on detector pixels. """
    idx = np.searchsorted(tchunk.pixel, pixels - 1e-6)
    ok = (idx < len(tchunk.pixel)) & \
        (np.abs(tchunk.pixel[np.minimum(idx, len(tchunk.pixel) - 1)] - pixels) < 1e-6)
    return idx, ok


def reconvolve(flux_fine, lsf):
    """ The forward model's operation: the deconvolved template, blurred. """
    return np.convolve(flux_fine, lsf, 'same')


def chunk_metrics(tchunk, ch, lsf, nflux, sig):
    """ How well one deconvolved chunk reproduces the observation.

    Args:
        tchunk (`TemplateChunk`): the deconvolved chunk (fine grid, padded).
        ch (`Chunk`): the corresponding observed chunk.
        lsf (ndarray): the LSF to re-convolve with.
        nflux, sig (ndarray): normalised flux and its 1-sigma for the order.

    Returns:
        dict
    """
    pix = ch.abspix
    idx, ok = fine_indices(tchunk, pix)
    model = reconvolve(tchunk.flux, lsf)[idx[ok]]
    obs = nflux[pix[ok]]
    s = sig[pix[ok]]
    use = np.isfinite(s)
    chi = (obs[use] - model[use]) / s[use]
    # Equivalent width over the unpadded chunk, observed vs deconvolved;
    # a deconvolution should move absorption around, not create or destroy it
    step = float(np.median(np.diff(tchunk.pixel)))
    core = (tchunk.pixel >= pix[0]) & (tchunk.pixel <= pix[-1])
    ew_obs = float(np.sum(1. - obs))
    ew_dec = float(np.sum(1. - tchunk.flux[core]) * step)
    return dict(chi2=float(np.mean(chi ** 2)) if use.any() else np.nan,
                npix=int(use.sum()),
                resid_rms=float(np.std(obs[use] - model[use])) if use.any() else np.nan,
                sig_med=float(np.median(s[use])) if use.any() else np.nan,
                ew_ratio=ew_dec / ew_obs if ew_obs > 0.5 else np.nan,
                min_obs=float(obs.min()),
                min_dec=float(tchunk.flux[core].min()))


def velocity_shift(ref_flux, ref_wave, flux):
    """ Linearised velocity of ``flux`` relative to ``ref_flux``, in m/s.

    Both on the same wavelength grid.  If flux(lambda) = ref(lambda (1-v/c)),
    i.e. ``flux`` is ref redshifted by v, then to first order
    flux - ref = -(v/c) dref/dln(lambda), and v follows by least squares.
    The sign convention is checked by :func:`_check_velocity_sign`.
    """
    g = np.gradient(ref_flux, np.log(ref_wave))
    den = np.sum(g * g)
    if den <= 0:
        return np.nan
    return float(-CKMS * 1e3 * np.sum((flux - ref_flux) * g) / den)


def _check_velocity_sign():
    """ Inject a known shift and recover it; raises if the sign is wrong. """
    w = np.linspace(5500., 5501., 2001)
    ref = 1. - 0.5 * np.exp(-0.5 * ((w - 5500.5) / 0.03) ** 2)
    v = 100.                                            # m/s
    shifted = 1. - 0.5 * np.exp(-0.5 * ((w * (1 - v / CKMS / 1e3) - 5500.5)
                                        / 0.03) ** 2)
    got = velocity_shift(ref, w, shifted)
    if not abs(got - v) < 2.:
        raise RuntimeError('velocity_shift sign/scale check failed: '
                           '{:.2f} m/s for an injected {:.1f}'.format(got, v))


# ---------------------------------------------------------------------------
# 5. The reference Jansson with a true adjoint (test only)
# ---------------------------------------------------------------------------

def jansson_true_adjoint(observed, lsf, niter, a=0.0, b=1.0, delta=0.1,
                         chi_change=1e-8):
    """ `pyodine.template.deconvolve.jansson`, with the adjoint made explicit.

    The fork's iteration is ``x += r * (K y - K K x)``, applying the kernel
    where the least-squares gradient needs its transpose, ``K^T y - K^T K x``.
    For a symmetric kernel the two are identical, and this function returns
    exactly what the fork's does (checked in :func:`adjoint_test`).  For an
    asymmetric one they differ.  Used only to measure that difference.
    """
    kernel = lsf
    kernel_t = lsf[::-1]
    old = observed
    guess = np.convolve(old, kernel, 'same')
    old = guess + delta * (1.0 - np.abs(guess - (a + b) / 2.) * 2. / (b - a)) * \
        (observed - np.convolve(guess, kernel, 'same'))
    chi = 1.
    k = -1
    while k <= niter:
        k += 1
        old_chi = chi
        relax = delta * (1.0 - np.abs(old - (a + b) / 2.0) * 2.0 / (b - a))
        guess = np.convolve(old, kernel, 'same')
        convol1 = np.convolve(observed, kernel_t, 'same')
        convol2 = np.convolve(guess, kernel_t, 'same')
        old += relax * (convol1 - convol2)
        chi = np.std((observed - old) ** 2)
        if abs((old_chi / chi) - 1.) < chi_change:
            k = niter + 1
    return old


# ---------------------------------------------------------------------------
# 6. The experiments
# ---------------------------------------------------------------------------

class Setup:
    """ Everything the experiments share, built once. """

    def __init__(self, arc, noise_scale=1.):
        t0 = time.time()
        self.obs = TemplateObservation()
        self.positions = self.obs.iodine_positions()
        #: measured / propagated noise (:func:`noise_check`); chi-square is
        #: reported against both, and the Monte Carlo injects the measured
        self.noise_scale = noise_scale
        self.norm, self.v_guess, self.solar = normalize(self.obs, self.positions)
        self.chunks = build_chunks(self.obs, self.positions)
        self.model = build_model()
        self.params, self.wave_resid = wave_params(self.model, self.chunks)
        self.deconvolver = pdeconv.ChunkedDeconvolver(self.chunks, self.model,
                                                      self.params)
        self.sig = {i: norm_sigma(self.obs, self.norm, i) for i in self.positions}
        self.nflux = {i: self.norm._flux[i] for i in self.positions}

        ech = [int(self.obs.ech_orders[i]) for i in self.positions]
        self.arc_fwhm = order_fwhm(arc, ech)
        self.t_setup = time.time() - t0

    def ech(self, ch):
        return int(self.obs.ech_orders[ch.order])

    def fwhm_arc(self, ch):
        return self.arc_fwhm[self.ech(ch)][0]

    def fwhm_stellar(self, ch):
        return CKMS / R_STELLAR / kms_per_pixel(ch)

    def indices(self, ech_orders=None):
        if ech_orders is None:
            return list(range(len(self.chunks)))
        return [i for i, ch in enumerate(self.chunks) if self.ech(ch) in ech_orders]

    def deconvolve(self, i, lsfs, norm=None, osample=None):
        return self.deconvolver.deconvolve_single_chunk(
            self.norm if norm is None else norm, i, deconv_pars(osample),
            lsf_fixed=lsfs)


def lsf_cases(S, osample=None):
    """ The LSFs compared, as (name, description, lsf array). """
    return [
        ('arc', 'Gaussian, ThAr FWHM (fiducial)',
         lsf_array(S.chunks, S.fwhm_arc, osample=osample)),
        ('arc_x0.8', 'Gaussian, 0.8 x ThAr FWHM',
         lsf_array(S.chunks, lambda ch: 0.8 * S.fwhm_arc(ch), osample=osample)),
        ('arc_x1.2', 'Gaussian, 1.2 x ThAr FWHM',
         lsf_array(S.chunks, lambda ch: 1.2 * S.fwhm_arc(ch), osample=osample)),
        ('stellar_R42k', 'Gaussian, phase-2 R = 42,000',
         lsf_array(S.chunks, S.fwhm_stellar, osample=osample)),
        ('arc_skew', 'Skewed, ThAr FWHM core + 25% red satellite',
         lsf_array(S.chunks, S.fwhm_arc, kind='skew', osample=osample)),
    ]


def run_lsf_comparison(S, subset, osample=None):
    """ Deconvolve the subset with every LSF; score each against the data.

    Every template is scored twice: re-convolved with its *own* LSF (does the
    deconvolution reproduce the data it was given?) and with the *fiducial*
    one (if the true LSF is the ThAr one, how wrong is a template built under
    another assumption?).  Velocities are relative to the fiducial template.
    """
    cases = lsf_cases(S, osample)
    fid = cases[0][2]
    results, fid_tchunks = [], {}
    for name, desc, lsfs in cases:
        t0 = time.time()
        for i in subset:
            ch = S.chunks[i]
            tc = S.deconvolve(i, lsfs, osample=osample)
            if name == 'arc':
                fid_tchunks[i] = tc
            own = chunk_metrics(tc, ch, lsfs[i], S.nflux[ch.order], S.sig[ch.order])
            cross = chunk_metrics(tc, ch, fid[i], S.nflux[ch.order], S.sig[ch.order])
            ref = fid_tchunks.get(i)
            core = (tc.pixel >= ch.abspix[0]) & (tc.pixel <= ch.abspix[-1])
            dv = velocity_shift(ref.flux[core], ref.wave[core], tc.flux[core]) \
                if ref is not None else np.nan
            k2 = S.noise_scale ** 2
            results.append(dict(case=name, osample=int(deconv_pars(osample)['osample_temp']),
                                chunk=i, order=S.ech(ch),
                                pix0=int(ch.abspix[0]),
                                fwhm_px=float(np.sqrt(np.sum(lsf_grid(osample) ** 2 * lsfs[i]))
                                              * 2.354820045),
                                chi2_own=own['chi2'], chi2_fid=cross['chi2'],
                                chi2_own_meas=own['chi2'] / k2,
                                chi2_fid_meas=cross['chi2'] / k2,
                                ew_ratio=own['ew_ratio'], min_obs=own['min_obs'],
                                min_dec=own['min_dec'], dv_vs_fid=dv))
        print('  {:14s} osample {:2d}  {:4d} chunks  {:6.1f} s   ({:s})'.format(
            name, int(deconv_pars(osample)['osample_temp']), len(subset),
            time.time() - t0, desc))
    return Table(results), fid_tchunks


def run_monte_carlo(S, subset, fid_tchunks, lsfs, nmc, seed=1998, osample=None):
    """ Noise propagation through the deconvolution, by Monte Carlo.

    Each realisation adds Gaussian noise at the propagated 1-sigma to the
    normalised observation and deconvolves again.  The spread across
    realisations is the deconvolution's response to noise of the size the
    observation actually carries.  Reported per chunk:

      * S/N per detector pixel of the observation, the deconvolved template
        (fine samples averaged back to detector pixels), and the template
        re-convolved with the LSF, which is what the forward model uses;
      * the velocity noise each of those three carries, m/s.
    """
    rng = np.random.default_rng(seed)
    step = int(deconv_pars(osample)['osample_temp'])
    store = {i: dict(dec=[], rec=[], obs=[]) for i in subset}
    for r in range(nmc):
        t0 = time.time()
        pert = components.NormalizedObservation(S.obs)
        noise = {}
        for p in S.positions:
            s = S.sig[p] * S.noise_scale
            n = np.where(np.isfinite(s), rng.normal(size=len(s)) * np.where(
                np.isfinite(s), s, 0.), 0.)
            noise[p] = n
            pert[p] = S.nflux[p] + n
        for i in subset:
            ch = S.chunks[i]
            tc = S.deconvolve(i, lsfs, norm=pert, osample=osample)
            store[i]['dec'].append(tc.flux.copy())
            store[i]['rec'].append(reconvolve(tc.flux, lsfs[i]))
            store[i]['obs'].append(noise[ch.order][ch.abspix])
        if r == 0 or r == nmc - 1:
            print('  realisation {:2d}/{:d}  {:5.1f} s'.format(r + 1, nmc, time.time() - t0))

    rows = []
    for i in subset:
        ch = S.chunks[i]
        ref = fid_tchunks[i]
        core = (ref.pixel >= ch.abspix[0]) & (ref.pixel <= ch.abspix[-1])
        idx, ok = fine_indices(ref, ch.abspix)
        dec = np.array(store[i]['dec'])
        rec = np.array(store[i]['rec'])
        nob = np.array(store[i]['obs'])
        nflux = S.nflux[ch.order][ch.abspix]
        sig = S.sig[ch.order][ch.abspix] * S.noise_scale
        use = np.isfinite(sig)

        # Deconvolved template, averaged back onto detector pixels: each
        # detector pixel gets the mean of the osample fine samples centred on
        # it.  The fine samples themselves are correlated and their own
        # scatter would overstate the per-pixel noise.
        dec_core = dec[:, core]
        npx = dec_core.shape[1] // step
        dec_px = dec_core[:, :npx * step].reshape(len(dec), npx, step).mean(axis=2)
        rec_px = rec[:, idx[ok]]

        snr_obs = np.median(nflux[use] / sig[use])
        snr_dec_fine = np.median(np.mean(dec_core, axis=0) / np.std(dec_core, axis=0))
        snr_dec = np.median(np.mean(dec_px, axis=0) / np.std(dec_px, axis=0))
        snr_rec = np.median(np.mean(rec_px, axis=0) / np.std(rec_px, axis=0))

        # Velocity noise: each realisation against the realisation mean
        w_core = ref.wave[core]
        v_dec = [velocity_shift(dec_core.mean(axis=0), w_core, d) for d in dec_core]
        rec_core = rec[:, core]
        v_rec = [velocity_shift(rec_core.mean(axis=0), w_core, d) for d in rec_core]
        w_px = ch.wave
        v_obs = [velocity_shift(nflux, w_px, nflux + n) for n in nob]
        rows.append(dict(chunk=i, order=S.ech(ch), pix0=int(ch.abspix[0]),
                         snr_obs=float(snr_obs), snr_dec=float(snr_dec),
                         snr_dec_fine=float(snr_dec_fine), snr_rec=float(snr_rec),
                         sigv_obs=float(np.std(v_obs)), sigv_dec=float(np.std(v_dec)),
                         sigv_rec=float(np.std(v_rec))))
    return Table(rows)


def run_adjoint_test(S, subset, fid_tchunks, snr, seed=7):
    """ Does the fork's Jansson handle an asymmetric kernel?

    A synthetic truth is taken from the fiducial deconvolved chunks, blurred
    with a known kernel, given noise at S/N ``snr``, and deconvolved twice:
    by the fork's `jansson` and by :func:`jansson_true_adjoint`.  With a
    symmetric kernel the two must agree exactly (a check of this harness);
    with the skewed one they need not.  Scored by the velocity shift of the
    recovered template against the truth and by its RMS error.
    """
    from pyodine.lib import misc
    rng = np.random.default_rng(seed)
    p = DECONV_PARS
    rows = []
    for i in subset:
        tc = fid_tchunks[i]
        truth = tc.flux.copy()
        wave = tc.wave
        ch = S.chunks[i]
        core = (tc.pixel >= ch.abspix[0]) & (tc.pixel <= ch.abspix[-1])
        # Detector pixels of the padded chunk, and where they sit on the grid
        pad = ch.padded.abspix
        idx, ok = fine_indices(tc, pad)
        cidx, cok = fine_indices(tc, ch.abspix)
        fw = S.fwhm_arc(ch)
        for kname, k in (('gauss', gaussian_lsf(fw)), ('skew', skewed_lsf(fw))):
            # Observe the truth as the detector would: blur, sample at the
            # pixels, add noise per pixel, then spline back onto the fine
            # grid exactly as ChunkedDeconvolver does with the real data.
            blurred = reconvolve(truth, k)
            obs_px = blurred[idx[ok]]
            obs_px = obs_px + rng.normal(size=len(obs_px)) * obs_px / snr
            obs_fine = misc.rebin(pad[ok].astype(float), obs_px, tc.pixel)
            obs_c = blurred[cidx[cok]]
            noisy_c = np.interp(ch.abspix[cok], pad[ok], obs_px)
            out = {}
            for aname, fn in (('fork', pdeconv.jansson), ('adjoint', jansson_true_adjoint)):
                rec = fn(obs_fine.copy(), k, p['jansson_niter'], a=p['jansson_zerolevel'],
                         b=p['jansson_contlevel'], delta=p['jansson_conver'],
                         chi_change=p['jansson_chi_change'])
                out[aname] = rec
                model_c = reconvolve(rec, k)[cidx[cok]]
                rows.append(dict(chunk=i, order=S.ech(ch), kernel=kname, solver=aname,
                                 dv_vs_truth=velocity_shift(truth[core], wave[core], rec[core]),
                                 rms_vs_truth=float(np.std(rec[core] - truth[core])),
                                 recon_chi2=float(np.mean(((model_c - noisy_c)
                                                           / (noisy_c / snr)) ** 2)),
                                 recon_bias=float(np.std(model_c - obs_c) / np.median(obs_c / snr))))
            if kname == 'gauss':
                diff = float(np.max(np.abs(out['fork'] - out['adjoint'])))
                if diff > 1e-10:
                    raise RuntimeError('Symmetric kernel: fork and reference '
                                       'Jansson differ by {:.2e}'.format(diff))
    return Table(rows)


def run_osample_scan(S, subset, osamples=(1, 2, 4, 10)):
    """ Re-convolution chi-square against the template oversampling factor.

    `ChunkedDeconvolver` puts the observation on its fine grid with
    `misc.rebin`, a cubic spline through the detector pixels.  HIRES lines
    are 2.2 pixels wide -- barely Nyquist-sampled -- and between pixel
    centres near a deep line the spline is wrong by many sigma at this S/N.
    Jansson then fits those wiggles, which no blurred spectrum can make.  At
    ``osample_temp = 1`` `rebin` is the identity and there is no spline.

    Returns:
        tuple: (`astropy.table.Table`_, dict osample -> {chunk: TemplateChunk})
    """
    rows, kept = [], {}
    for os_ in osamples:
        lsfs = lsf_array(S.chunks, S.fwhm_arc, osample=os_)
        kept[os_] = {}
        t0 = time.time()
        for i in subset:
            ch = S.chunks[i]
            tc = S.deconvolve(i, lsfs, osample=os_)
            kept[os_][i] = tc
            m = chunk_metrics(tc, ch, lsfs[i], S.nflux[ch.order], S.sig[ch.order])
            rows.append(dict(osample=os_, chunk=i, order=S.ech(ch),
                             chi2_meas=m['chi2'] / S.noise_scale ** 2,
                             ew_ratio=m['ew_ratio'], min_obs=m['min_obs'],
                             min_dec=m['min_dec']))
        print('  osample {:2d}: {:d} chunks in {:.1f} s'.format(os_, len(subset), time.time() - t0))
    return Table(rows), kept


def run_full(S, lsfs, outfile, osample=None):
    """ A full template: every iodine chunk, through `deconvolve_obs`. """
    t0 = time.time()
    template = S.deconvolver.deconvolve_obs(S.norm, S.v_guess, S.obs.bary_vel_corr,
                                            lsf_fixed=lsfs, deconv_pars=deconv_pars(osample))
    runtime = time.time() - t0
    template.save(outfile)
    rows = []
    for i, (tc, ch) in enumerate(zip(template.chunks, S.chunks)):
        m = chunk_metrics(tc, ch, lsfs[i], S.nflux[ch.order], S.sig[ch.order])
        rows.append(dict(chunk=i, order=S.ech(ch), pix0=int(ch.abspix[0]),
                         wave=float(np.median(ch.wave)),
                         chi2=m['chi2'], chi2_meas=m['chi2'] / S.noise_scale ** 2,
                         npix=m['npix'], resid_rms=m['resid_rms'],
                         sig_med=m['sig_med'], ew_ratio=m['ew_ratio'],
                         min_obs=m['min_obs'], min_dec=m['min_dec'],
                         weight=float(tc.weight),
                         wave_resid_ms=float(S.wave_resid[i])))
    return template, Table(rows), runtime


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_figure(S, arc, lsfcmp, osc, mc, fid_tchunks, kept, outfile):
    """ Six panels: the instrumental width; one deep-line chunk and its
    residuals at oversampling 10 and 1; the S/N cost; the LSF assumption;
    and the re-convolution chi-square against line depth. """
    import matplotlib.pyplot as plt

    INK, C1, C2, C3 = '#0b0b0b', '#2a78d6', '#eb6834', '#1baf7a'
    fig, axs = plt.subplots(2, 3, figsize=(16, 9))
    k = S.noise_scale

    # (a) instrumental vs stellar width
    ax = axs[0, 0]
    ax.plot(arc['order'], arc['fwhm_kms'], '.', color=C1, alpha=0.25, ms=3,
            label='ThAr lines, 1998-08-26')
    orders = np.unique(arc['order'])
    ax.plot(orders, [np.median(arc['fwhm_kms'][arc['order'] == o]) for o in orders],
            '-', color=C1, lw=2, label='ThAr, median per order')
    ax.axhline(CKMS / R_STELLAR, color=C2, ls='--',
               label='phase 2, weak stellar lines (R = 42,000)')
    ax.axvspan(min(S.ech(c) for c in S.chunks) - 0.5,
               max(S.ech(c) for c in S.chunks) + 0.5, color='0.9', zorder=0,
               label='iodine orders')
    ax.set_xlabel('echelle order')
    ax.set_ylabel('FWHM (km/s)')
    ax.set_ylim(2, 10)
    ax.legend(fontsize=8, loc='upper right')
    ax.set_title('(a) instrumental width: R = {:.0f}'.format(
        CKMS / np.median(arc['fwhm_kms'])), loc='left', fontsize=10)

    # (b) + (c) a deep-line chunk, re-convolved at oversampling 10 and 1
    t10 = osc[osc['osample'] == 10]
    i = int(t10['chunk'][np.argmax(t10['chi2_meas'])])
    ch = S.chunks[i]
    nflux = S.nflux[ch.order]
    sig = S.sig[ch.order] * k
    ax, axr = axs[0, 1], axs[0, 2]
    ax.plot(ch.wave, nflux[ch.abspix], 'o', ms=4, color=INK, label='observed', zorder=3)
    for os_, col in ((10, C2), (1, C1)):
        tc = kept[os_][i]
        lsf = gaussian_lsf(S.fwhm_arc(ch), os_)
        ax.plot(tc.wave, tc.flux, '-', color=col, lw=0.8, alpha=0.6,
                label='deconvolved, oversampling {:d}'.format(os_))
        idx, ok = fine_indices(tc, ch.abspix)
        model = reconvolve(tc.flux, lsf)[idx[ok]]
        r = (nflux[ch.abspix[ok]] - model) / sig[ch.abspix[ok]]
        axr.plot(ch.wave[ok], r, 'o-', ms=3, lw=0.8, color=col,
                 label='oversampling {:d}: chi2 = {:.1f}'.format(os_, np.mean(r ** 2)))
    ax.set_xlim(ch.wave[0], ch.wave[-1])
    ax.set_xlabel('vacuum wavelength, observed frame (A)')
    ax.set_ylabel('normalised flux')
    ax.legend(fontsize=8, loc='lower left')
    ax.set_title('(b) order {:d}, chunk {:d}: deconvolved'.format(S.ech(ch), i),
                 loc='left', fontsize=10)
    axr.axhline(0, color='0.5', lw=0.8)
    axr.set_xlim(ch.wave[0], ch.wave[-1])
    axr.set_xlabel('vacuum wavelength, observed frame (A)')
    axr.set_ylabel('(observed - re-convolved) / measured sigma')
    axr.legend(fontsize=8)
    axr.set_title('(c) same chunk: does it reproduce the observation?', loc='left',
                  fontsize=10)

    # (d) S/N cost
    ax = axs[1, 0]
    for os_, mk in ((10, 'o'), (1, 's')):
        t = mc[mc['osample'] == os_]
        x = np.arange(len(t))
        if os_ == 10:
            ax.plot(x, t['snr_obs'], mk, color=INK, ms=3, label='observed')
        ax.plot(x, t['snr_dec'], mk, color=C2, ms=3, alpha=0.7,
                label='deconvolved, osample {:d}'.format(os_))
        ax.plot(x, t['snr_rec'], mk, color=C1, ms=3, alpha=0.7,
                label='deconvolved + re-convolved, osample {:d}'.format(os_))
    ax.set_yscale('log')
    ax.set_xlabel('chunk (orders {:s})'.format(', '.join(str(o) for o in SUBSET_ORDERS)))
    ax.set_ylabel('median S/N per detector pixel')
    ax.legend(fontsize=7)
    ax.set_title('(d) S/N cost of deconvolution (Monte Carlo)', loc='left', fontsize=10)

    # (e) the LSF assumption: template built with LSF X, blurred with ThAr
    ax = axs[1, 1]
    lsf1 = lsfcmp[lsfcmp['osample'] == 1]
    names = list(dict.fromkeys(lsf1['case']))
    for j, name in enumerate(names):
        t = lsf1[lsf1['case'] == name]
        ax.plot(np.full(len(t), j - 0.12), t['chi2_own_meas'], '.', color=C1, alpha=0.35, ms=4)
        ax.plot(np.full(len(t), j + 0.12), t['chi2_fid_meas'], '.', color=C2, alpha=0.35, ms=4)
        ax.plot([j - 0.25, j], [np.median(t['chi2_own_meas'])] * 2, '-', color=C1, lw=2)
        ax.plot([j, j + 0.25], [np.median(t['chi2_fid_meas'])] * 2, '-', color=C2, lw=2)
    ax.plot([], [], '-', color=C1, lw=2, label='re-convolved with its own LSF')
    ax.plot([], [], '-', color=C2, lw=2, label='re-convolved with the ThAr LSF')
    ax.axhline(1., color='0.5', lw=0.8)
    ax.set_yscale('log')
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel('chi-square per pixel (measured noise)')
    ax.legend(fontsize=8, loc='upper left')
    ax.set_title('(e) the assumed LSF (oversampling 1)', loc='left', fontsize=10)

    # (f) chi2 against line depth, oversampling 10 vs 1
    ax = axs[1, 2]
    for os_, col in ((10, C2), (4, C3), (1, C1)):
        t = osc[osc['osample'] == os_]
        ax.plot(t['min_obs'], t['chi2_meas'], '.', color=col, ms=5, alpha=0.6,
                label='oversampling {:d}: median {:.1f}'.format(os_, np.median(t['chi2_meas'])))
    ax.axhline(1., color='0.5', lw=0.8)
    ax.set_yscale('log')
    ax.set_xlabel('deepest observed line in the chunk (normalised flux)')
    ax.set_ylabel('chi-square per pixel (measured noise)')
    ax.legend(fontsize=8)
    ax.set_title('(f) the spline sets the floor', loc='left', fontsize=10)

    fig.tight_layout()
    fig.savefig(outfile, dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def summarize(tab, col, fmt='{:.3f}'):
    v = np.asarray(tab[col], dtype=float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return 'n/a'
    return (fmt + ' [' + fmt + ', ' + fmt + ']').format(
        np.median(v), np.percentile(v, 16), np.percentile(v, 84))


def parser(options=None):
    p = argparse.ArgumentParser(
        description='Deconvolve the HD 187123 template (phase 3, prompt 3)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--nmc', type=int, default=12,
                   help='Monte Carlo realisations for the S/N cost')
    p.add_argument('--subset', type=int, nargs='+', default=list(SUBSET_ORDERS),
                   help='Echelle orders for the MC and the LSF comparison')
    p.add_argument('--no-full', action='store_true',
                   help='Skip deconvolving the full iodine template')
    p.add_argument('--datadir', type=str, default=DATA_DIR)
    p.add_argument('--figdir', type=str, default=FIG_DIR)
    return p.parse_args() if options is None else p.parse_args(options)


def main(pargs):
    _check_velocity_sign()

    print('1. ThAr line widths, 1998-08-26')
    arc = arc_line_widths()
    arc.write(os.path.join(pargs.datadir, 'deconv_arc_widths.csv'), overwrite=True)
    print('  {:d} clean lines in {:d} orders'.format(len(arc), len(np.unique(arc['order']))))
    print('  FWHM {} px, {} km/s  ->  R = {:.0f}'.format(
        summarize(arc, 'fwhm_px'), summarize(arc, 'fwhm_kms', '{:.2f}'),
        CKMS / np.median(arc['fwhm_kms'])))
    print('  PypeIt\'s own per-order FWHM: median {:.3f} px'.format(
        np.median(np.unique(arc['pypeit_fwhm_px']))))

    print('2. Is the propagated noise the right size?  Exposure differences')
    obs = TemplateObservation()
    noise = noise_check(obs, obs.iodine_positions())
    noise.write(os.path.join(pargs.datadir, 'deconv_noise_check.csv'), overwrite=True)
    for pair in ('01', '12', '02'):
        t = noise[noise['pair'] == pair]
        print('  exposures {:s}: {:2d} orders  width {:s}  shift {:s} m/s'.format(
            pair, len(t), summarize(t, 'width'), summarize(t, 'dv_ms', '{:.1f}')))
    # The clean pair sets the scale; the third exposure differs from both
    # others by more than noise and is not used for it
    k = float(np.median(noise['width'][noise['pair'] == '01']))
    print('  measured / propagated noise = {:.3f}  (propagated S/N is low by {:.2f}x)'.format(
        k, 1. / k))

    print('3. Template, normalisation, chunks')
    S = Setup(arc, noise_scale=k)
    print('  {:d} iodine orders ({:d}-{:d}), {:d} chunks, setup {:.1f} s'.format(
        len(S.positions), S.obs.ech_orders[S.positions].min(),
        S.obs.ech_orders[S.positions].max(), len(S.chunks), S.t_setup))
    print('  velocity guess vs {:s}: {:+.3f} km/s'.format(REF_SPECTRUM, S.v_guess / 1e3))
    print('  continuum from solar comparison in {:d} of {:d} orders, '
          'envelope fallback in the rest'.format(len(S.solar), len(S.positions)))
    print('  linear wave model, padded chunk: max residual {} m/s'.format(
        summarize(Table(dict(r=S.wave_resid)), 'r', '{:.1f}')))
    fb = [o for o, (f, flag) in S.arc_fwhm.items() if flag]
    print('  arc FWHM per iodine order: {:.2f}-{:.2f} px{}'.format(
        min(f for f, _ in S.arc_fwhm.values()), max(f for f, _ in S.arc_fwhm.values()),
        '' if not fb else ' (fallback for {})'.format(fb)))

    subset = S.indices(pargs.subset)
    print('4. LSF comparison: {:d} chunks in orders {}'.format(len(subset), pargs.subset))
    from astropy.table import vstack
    lsf10, fid_tchunks = run_lsf_comparison(S, subset, osample=10)
    lsf1, _ = run_lsf_comparison(S, subset, osample=1)
    lsfcmp = vstack([lsf10, lsf1])
    lsfcmp.write(os.path.join(pargs.datadir, 'deconv_lsf_cases.csv'), overwrite=True)
    print('  chi-square per pixel against the MEASURED noise; median [16th, 84th]')
    for os_, tab in ((10, lsf10), (1, lsf1)):
        print('  -- oversampling {:d}'.format(os_))
        for name in dict.fromkeys(tab['case']):
            t = tab[tab['case'] == name]
            print('  {:14s} own LSF {:s}  ThAr LSF {:s}  EW ratio {:s}  '
                  'dv {:s} m/s'.format(name, summarize(t, 'chi2_own_meas', '{:.1f}'),
                                       summarize(t, 'chi2_fid_meas', '{:.1f}'),
                                       summarize(t, 'ew_ratio'),
                                       summarize(t, 'dv_vs_fid', '{:.1f}')))

    print('5. Template oversampling: does the spline set the floor?')
    osc, kept = run_osample_scan(S, subset)
    osc.write(os.path.join(pargs.datadir, 'deconv_osample.csv'), overwrite=True)
    bins = ((0., 0.6), (0.6, 0.9), (0.9, 2.))
    print('  median chi2 (measured noise) by the deepest observed line in the chunk')
    print('  osample  ' + '  '.join('min in [{:.1f},{:.1f}) '.format(*b) for b in bins) + '  all')
    for os_ in dict.fromkeys(osc['osample']):
        t = osc[osc['osample'] == os_]
        cells = []
        for lo, hi in bins:
            u = t[(t['min_obs'] >= lo) & (t['min_obs'] < hi)]
            cells.append('{:6.1f} (n={:3d})   '.format(np.median(u['chi2_meas']), len(u))
                         if len(u) else '   n/a            ')
        print('  {:6d}   '.format(int(os_)) + ''.join(cells) +
              '{:.1f}'.format(np.median(t['chi2_meas'])))

    fid_lsf = lsf_cases(S)[0][2]
    mcs = []
    for os_, ref in ((10, fid_tchunks), (1, kept[1])):
        print('6. Monte Carlo S/N cost, oversampling {:d}, {:d} realisations, '
              'measured noise'.format(os_, pargs.nmc))
        lsfs = fid_lsf if os_ == 10 else lsf_array(S.chunks, S.fwhm_arc, osample=os_)
        mc = run_monte_carlo(S, subset, ref, lsfs, pargs.nmc, osample=os_)
        mc['osample'] = os_
        mcs.append(mc)
        for col in ('snr_obs', 'snr_dec', 'snr_dec_fine', 'snr_rec',
                    'sigv_obs', 'sigv_dec', 'sigv_rec'):
            print('  {:13s} {:s}'.format(col, summarize(mc, col, '{:.1f}')))
    mc = vstack(mcs)
    mc.write(os.path.join(pargs.datadir, 'deconv_snr_mc.csv'), overwrite=True)

    snr = float(np.median(mc['snr_obs']))
    print('7. Adjoint test (synthetic, S/N {:.0f}, oversampling 10)'.format(snr))
    adj = run_adjoint_test(S, subset[::3], fid_tchunks, snr=snr)
    adj.write(os.path.join(pargs.datadir, 'deconv_adjoint.csv'), overwrite=True)
    for kname in ('gauss', 'skew'):
        for sv in ('fork', 'adjoint'):
            t = adj[(adj['kernel'] == kname) & (adj['solver'] == sv)]
            print('  {:5s} {:7s} dv vs truth {:s} m/s  rms vs truth {:s}  '
                  'recon chi2 {:s}  recon bias {:s} sigma'.format(
                      kname, sv, summarize(t, 'dv_vs_truth', '{:.1f}'),
                      summarize(t, 'rms_vs_truth', '{:.4f}'),
                      summarize(t, 'recon_chi2', '{:.2f}'),
                      summarize(t, 'recon_bias', '{:.2f}')))
    f = adj[(adj['kernel'] == 'skew') & (adj['solver'] == 'fork')]
    a = adj[(adj['kernel'] == 'skew') & (adj['solver'] == 'adjoint')]
    diff = np.asarray(f['dv_vs_truth']) - np.asarray(a['dv_vs_truth'])
    print('  skewed kernel, fork minus true adjoint, per chunk: {:.1f} [{:.1f}, {:.1f}] m/s'.format(
        np.median(diff), *np.percentile(diff, [16, 84])))

    fulls = {}
    if not pargs.no_full:
        for os_ in (10, 1):
            print('8. Full template, ThAr LSF, every iodine chunk, oversampling {:d}'.format(os_))
            lsfs = fid_lsf if os_ == 10 else lsf_array(S.chunks, S.fwhm_arc, osample=os_)
            template, full, runtime = run_full(S, lsfs, TEMPLATE_H5[os_], osample=os_)
            full['osample'] = os_
            fulls[os_] = full
            print('  {:d} chunks in {:.0f} s -> {:s}'.format(len(template), runtime,
                                                              TEMPLATE_H5[os_]))
            print('  chi2 vs propagated {:s}   vs measured {:s}'.format(
                summarize(full, 'chi2', '{:.2f}'), summarize(full, 'chi2_meas', '{:.1f}')))
            print('  EW ratio {:s}   min deconvolved flux {:s}'.format(
                summarize(full, 'ew_ratio'), summarize(full, 'min_dec')))
            print('  chunks with chi2 (measured) > 4: {:d} of {:d}; any flux < 0: {:d}'.format(
                int(np.sum(full['chi2_meas'] > 4)), len(full), int(np.sum(full['min_dec'] < 0))))
            back = StellarTemplate_Chunked(TEMPLATE_H5[os_])
            same = all(np.array_equal(x.flux, y.flux) for x, y in zip(back.chunks, template.chunks))
            print('  read back by StellarTemplate_Chunked: {:d} chunks, osample {}, '
                  'identical: {}'.format(len(back), back.osample, same))
        full = vstack([fulls[10], fulls[1]])
        full.write(os.path.join(pargs.datadir, 'deconv_chunks.csv'), overwrite=True)
    else:
        full = None

    make_figure(S, arc, lsfcmp, osc, mc, fid_tchunks, kept,
                os.path.join(pargs.figdir, 'fig_p3_deconv.png'))
    print('Figure: {:s}'.format(os.path.join(pargs.figdir, 'fig_p3_deconv.png')))


if __name__ == '__main__':
    main(parser())
