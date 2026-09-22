""" Figures for the public-facing document docs/public_HD187123b.md

Six figures on the discovery of HD 187123 b, the first exoplanet found with
Keck/HIRES (Butler, Marcy, Vogt & Apps 1998, PASP, 110, 1389).

Figures 1, 2, 3 and 6 are schematics generated entirely from this script.
Figures 4 and 5 plot real Keck/HIRES velocities of HD 187123 taken from:

    Teklu et al. 2025, A&A, 702, A68 -- "An updated catalog of HIRES/Keck
    radial velocity measurements", VizieR catalog J/A+A/702/A68, table A.1.

The velocity table was retrieved 2026-09-19 with:

    curl "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?\
-source=J/A%2BA/702/A68/tablea1&Name=HD187123\
&-out=Name,BJD,RV,e_RV,RVcor,e_RVcor&-out.max=unlimited&-sort=BJD"

and is stored at first_hires_exoplanet/data/hd187123_hires_rv.tsv.

NOTE ON UNITS: VizieR declares the RV columns as km/s.  That is wrong -- the
velocities are relative and in m/s, as is obvious from their range (about
-100 to +90) against the published 72 m/s semiamplitude.  We read them as m/s.

Run with:

    conda run -n astro python first_hires_exoplanet/figs.py

"""

# Standard imports
import os

import numpy as np
from astropy.time import Time

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Orbital period of HD 187123 b (d).  Butler et al. (1998) quote 3.097 d; this
#: is the refined modern value.
PERIOD = 3.0965828

#: Published semiamplitude (m/s) from Butler et al. (1998)
K_PUBLISHED = 72.

#: Semimajor axis of HD 187123 b (AU) and of Mercury, for scale
A_PLANET = 0.042
A_MERCURY = 0.387

#: BJD of the paper's publication date, 1998 Oct 22 -- the discovery-era cut
BJD_PUBLICATION = 2451110.

#: Colours.  Slots 1-3 of the dataviz skill's validated categorical palette,
#: taken in fixed order, plus its light-mode ink and surface tokens.
C_DATA = '#2a78d6'      # slot 1, blue
C_MODEL = '#eb6834'     # slot 2, orange
C_THIRD = '#1baf7a'     # slot 3, aqua
INK = '#0b0b0b'
INK_2 = '#52514e'
SURFACE = '#fcfcfb'

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA_FILE = os.path.join(_HERE, 'data', 'hd187123_hires_rv.tsv')
FIG_DIR = os.path.join(os.path.dirname(_HERE), 'docs', 'figs')


# ---------------------------------------------------------------------------
# Plot style
# ---------------------------------------------------------------------------

def set_style():
    """ Apply a common, deliberately recessive style to all figures. """
    plt.rcParams.update({
        'figure.facecolor': SURFACE,
        'axes.facecolor': SURFACE,
        'savefig.facecolor': SURFACE,
        'font.size': 11,
        'axes.labelsize': 11,
        'axes.titlesize': 12,
        'axes.labelcolor': INK,
        'axes.edgecolor': '#c9c8c3',
        'axes.linewidth': 0.8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'text.color': INK,
        'xtick.color': INK_2,
        'ytick.color': INK_2,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'grid.color': '#e6e5e0',
        'grid.linewidth': 0.8,
        'legend.frameon': False,
        'legend.fontsize': 10,
    })


def _save(fig, name):
    """ Write a figure to docs/figs and report it. """
    os.makedirs(FIG_DIR, exist_ok=True)
    outfile = os.path.join(FIG_DIR, name)
    fig.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('Wrote {:s}'.format(outfile))


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def load_velocities(discovery_era=False):
    """ Read the Keck/HIRES velocities of HD 187123.

    Args:
        discovery_era (bool):
            If True, keep only epochs before the 1998 Oct 22 publication of
            Butler et al. (1998).

    Returns:
        tuple: (bjd, rv, erv) as `numpy.ndarray`.  Velocities are relative,
        in m/s, using the NZP-corrected columns.
    """
    bjd, rv, erv = [], [], []
    with open(DATA_FILE, 'r') as ff:
        for line in ff:
            if not line.startswith('HD187123'):
                continue
            items = line.split('\t')
            # Name, BJD, RV, e_RV, RVcor, e_RVcor
            bjd.append(float(items[1]))
            rv.append(float(items[4]))      # NZP corrected
            erv.append(float(items[5]))
    bjd, rv, erv = np.array(bjd), np.array(rv), np.array(erv)

    if discovery_era:
        keep = bjd < BJD_PUBLICATION
        bjd, rv, erv = bjd[keep], rv[keep], erv[keep]

    return bjd, rv, erv


def fit_circular(bjd, rv, erv, period=PERIOD):
    """ Least-squares fit of a circular orbit at fixed period.

    The model is v = gamma + Kc*cos(2 pi t/P) + Ks*sin(2 pi t/P), which is
    linear in its three parameters, so no starting guess is needed.  We fit the
    phase rather than adopting one, because the discovery paper's reference
    epoch is not in the machine-readable record.

    Args:
        bjd, rv, erv (`numpy.ndarray`): epochs (d) and velocities (m/s).
        period (float): orbital period held fixed (d).

    Returns:
        tuple: (K, phase0, gamma) -- semiamplitude (m/s), phase of maximum
        velocity, and systemic offset (m/s).
    """
    ang = 2 * np.pi * bjd / period
    design = np.column_stack([np.ones_like(bjd), np.cos(ang), np.sin(ang)])
    weight = 1. / erv
    coeff, *_ = np.linalg.lstsq(design * weight[:, None], rv * weight, rcond=None)
    gamma, kcos, ksin = coeff
    semi = np.hypot(kcos, ksin)
    # v = gamma + K cos(2 pi (t/P - phase0)) expands to
    #     gamma + K cos(2 pi phase0) cos(ang) + K sin(2 pi phase0) sin(ang),
    # so kcos = K cos(2 pi phase0) and ksin = K sin(2 pi phase0).
    phase0 = np.arctan2(ksin, kcos) / (2 * np.pi)
    return semi, phase0, gamma


# ---------------------------------------------------------------------------
# Figure 1 -- the wobble
# ---------------------------------------------------------------------------

def fig_wobble():
    """ Star and planet about their common centre of mass, and the velocity
    curve an observer sees. """
    fig, axes = plt.subplots(1, 2, figsize=(10., 4.), width_ratios=[1., 1.3])

    # -- Left: the orbit, seen from above, star wobble exaggerated
    ax = axes[0]
    theta = np.linspace(0, 2 * np.pi, 400)
    r_planet, r_star = 1.0, 0.28        # not to scale; see caption
    ax.plot(r_planet * np.cos(theta), r_planet * np.sin(theta),
            color=C_MODEL, lw=1.2, ls='--', alpha=0.7)
    ax.plot(-r_star * np.cos(theta), -r_star * np.sin(theta),
            color=C_DATA, lw=1.2, ls='--', alpha=0.7)

    # Phase 0 is the moment the star moves straight towards Earth
    phase = 0.
    x_star, y_star = -r_star * np.cos(phase), -r_star * np.sin(phase)
    ax.plot(r_planet * np.cos(phase), r_planet * np.sin(phase), 'o',
            color=C_MODEL, ms=11, alpha=0.45)
    ax.plot(x_star, y_star, 'o', color=C_DATA, ms=22)
    ax.plot(0, 0, '+', color=INK_2, ms=11, mew=1.5)

    ax.text(r_planet * np.cos(phase) + 0.12, r_planet * np.sin(phase) - 0.02,
            'planet\n(never seen directly)', color=INK_2, fontsize=9,
            va='center')
    ax.text(-1.55, 1.02, 'star (all we can see)', color=INK, fontsize=10,
            ha='left')
    # A ghost of the star half an orbit later
    ax.plot(-x_star, -y_star, 'o', color=C_DATA, ms=22, alpha=0.25)
    ax.text(-0.05, 0.36, 'centre of mass', color=INK_2, fontsize=9, ha='right')

    # The star's motion now, and half an orbit later
    ax.annotate('', xy=(x_star, y_star - 0.62), xytext=(x_star, y_star - 0.14),
                arrowprops=dict(arrowstyle='-|>', color=C_DATA, lw=2.4))
    ax.text(x_star, y_star - 0.72, 'moving\ntowards us', color=C_DATA,
            fontsize=9.5, ha='center', va='top')
    ax.annotate('', xy=(-x_star, 0.62), xytext=(-x_star, 0.14),
                arrowprops=dict(arrowstyle='-|>', color=C_DATA, lw=2.4,
                                alpha=0.45))
    ax.text(-x_star, 0.72, 'moving away', color=C_DATA, fontsize=9.5,
            ha='center', va='bottom', alpha=0.8)

    ax.annotate('', xy=(0., -1.55), xytext=(0., -1.15),
                arrowprops=dict(arrowstyle='-|>', color=INK_2, lw=1.2))
    ax.text(0.08, -1.42, 'to Earth', color=INK_2, fontsize=9)

    ax.set_xlim(-1.6, 1.7)
    ax.set_ylim(-1.7, 1.35)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Only the star is visible', color=INK, loc='left')

    # -- Right: the velocity curve.  Phase 0 matches the left-hand panel.
    ax = axes[1]
    tt = np.linspace(0, 2, 500)
    ax.axhline(0., color='#c9c8c3', lw=0.8)
    ax.plot(tt, K_PUBLISHED * np.sin(2 * np.pi * tt), color=C_DATA, lw=2.)
    ax.annotate('', xy=(0.25, 0.), xytext=(0.25, K_PUBLISHED),
                arrowprops=dict(arrowstyle='<->', color=INK_2, lw=1.))
    ax.text(0.30, K_PUBLISHED / 2, 'K = 72 m s$^{-1}$', color=INK, fontsize=10,
            va='center')
    ax.text(1.62, 40, 'moving towards us', color=C_DATA, fontsize=9,
            ha='center')
    ax.text(1.62, -46, 'moving away from us', color=C_DATA, fontsize=9,
            ha='center')

    ax.set_xlabel('Orbits completed')
    ax.set_ylabel('Velocity of the star\ntowards or away from Earth (m s$^{-1}$)')
    ax.set_ylim(-95, 95)
    ax.grid(axis='y', alpha=0.6)
    ax.set_title('What a spectrograph records', color=INK, loc='left')

    fig.tight_layout()
    _save(fig, 'fig1_wobble.png')


# ---------------------------------------------------------------------------
# Figure 2 -- the Doppler shift, to scale
# ---------------------------------------------------------------------------

def _line_profile(wave, centre, depth=0.65, fwhm=0.15):
    """ A single Gaussian absorption line on a unit continuum. """
    sigma = fwhm / 2.3548
    return 1. - depth * np.exp(-0.5 * ((wave - centre) / sigma) ** 2)


def fig_doppler():
    """ How big the 72 m/s shift actually is next to a stellar line.

    Wavelengths are quoted in nanometres for a general audience.  Internally the
    line profile is still built in Angstroms, so the plotted axes convert:
    1 nm = 10 A.
    """
    c_light = 2.99792458e8              # m/s
    lam0 = 5500.                        # A  (= 550 nm)
    dlam = lam0 * K_PUBLISHED / c_light  # A
    dlam_nm = dlam / 10.                 # nm

    fig, axes = plt.subplots(1, 2, figsize=(10., 4.))

    wave = np.linspace(lam0 - 0.6, lam0 + 0.6, 4000)
    rest = _line_profile(wave, lam0)
    shifted = _line_profile(wave, lam0 + dlam)

    # -- Left: the whole line.  The two curves are indistinguishable.
    ax = axes[0]
    ax.plot((wave - lam0) / 10., rest, color=C_DATA, lw=2., label='at rest')
    ax.plot((wave - lam0) / 10., shifted, color=C_MODEL, lw=2., ls='--',
            label='moving at 72 m s$^{-1}$')
    ax.set_xlabel('Wavelength offset from 550 nm (nm)')
    ax.set_ylabel('Fraction of light transmitted')
    ax.set_ylim(0.3, 1.05)
    ax.legend(loc='lower right')
    ax.grid(alpha=0.6)
    ax.set_title('One absorption line', color=INK, loc='left')
    ax.text(-0.055, 0.36,
            'shift = 0.00013 nm\nline width $\\approx$ 0.015 nm',
            color=INK_2, fontsize=9, va='bottom')

    # -- Right: zoom onto the steep flank, where the shift is finally visible
    ax = axes[1]
    lo, hi = -0.080, -0.070          # A
    fine = np.linspace(lam0 + lo, lam0 + hi, 2000)
    xx = (fine - lam0) / 10. * 1e4   # units of 1e-4 nm
    ax.plot(xx, _line_profile(fine, lam0), color=C_DATA, lw=2., label='at rest')
    ax.plot(xx, _line_profile(fine, lam0 + dlam), color=C_MODEL, lw=2.,
            ls='--', label='moving at 72 m s$^{-1}$')

    x_arrow = -75.                   # in 1e-4 nm
    y_arrow = _line_profile(np.array([lam0 - 0.075]), lam0)[0]
    ax.annotate('', xy=(x_arrow + dlam_nm * 1e4, y_arrow),
                xytext=(x_arrow, y_arrow),
                arrowprops=dict(arrowstyle='<->', color=INK, lw=1.3))
    ax.text(x_arrow + 0.2, y_arrow + 0.008,
            'the entire signal: 0.00013 nm', color=INK, fontsize=10)
    ax.set_xlabel(r'Wavelength offset (10$^{-4}$ nm)')
    ax.set_ylabel('Fraction of light transmitted')
    ax.grid(alpha=0.6)
    ax.legend(loc='upper right')
    ax.set_title(r'The same line, magnified 60$\times$', color=INK, loc='left')

    fig.tight_layout()
    _save(fig, 'fig2_doppler.png')


# ---------------------------------------------------------------------------
# Figure 3 -- the iodine cell
# ---------------------------------------------------------------------------

def _synthetic_iodine(wave, seed=1998, nlines=420):
    """ An illustrative iodine absorption forest.

    This is *not* a real FTS scan -- no public machine-readable atlas was
    available.  It reproduces the character of the I2 spectrum (a dense thicket
    of narrow lines) for teaching purposes only.
    """
    rand = np.random.default_rng(seed)
    flux = np.ones_like(wave)
    centres = rand.uniform(wave.min(), wave.max(), nlines)
    depths = rand.uniform(0.1, 0.75, nlines)
    widths = rand.uniform(0.02, 0.05, nlines)
    for cen, dep, wid in zip(centres, depths, widths):
        flux *= _line_profile(wave, cen, depth=dep, fwhm=wid)
    return flux


def fig_iodine():
    """ Stellar spectrum x iodine transmission = what the detector records. """
    wave = np.linspace(5200., 5210., 6000)

    star = np.ones_like(wave)
    for cen, dep in [(5201.4, 0.55), (5203.1, 0.35), (5204.8, 0.7),
                     (5206.6, 0.3), (5208.2, 0.6), (5209.4, 0.4)]:
        star *= _line_profile(wave, cen, depth=dep, fwhm=0.16)
    iodine = _synthetic_iodine(wave)

    fig, axes = plt.subplots(3, 1, figsize=(10., 6.5), sharex=True)

    axes[0].plot(wave, star, color=C_DATA, lw=1.6)
    axes[0].set_ylabel('Star')
    axes[0].set_title('The starlight we want to measure', color=INK, loc='left')

    axes[1].plot(wave, iodine, color=C_THIRD, lw=0.9)
    axes[1].set_ylabel('Iodine')
    axes[1].set_title('The iodine cell: thousands of lines at known wavelengths, '
                      'fixed forever', color=INK, loc='left')

    axes[2].plot(wave, star * iodine, color=INK, lw=0.9)
    axes[2].set_ylabel('Recorded')
    axes[2].set_xlabel(r'Wavelength ($\AA$)')
    axes[2].set_title('What HIRES actually records: the two, multiplied',
                      color=INK, loc='left')

    for ax in axes:
        ax.set_ylim(-0.05, 1.12)
        ax.grid(alpha=0.5)

    fig.tight_layout()
    _save(fig, 'fig3_iodine.png')


# ---------------------------------------------------------------------------
# Figures 4 and 5 -- the real velocities
# ---------------------------------------------------------------------------

def fig_july1998():
    """ The July 1998 run: ten velocities over five nights. """
    bjd, rv, erv = load_velocities(discovery_era=True)
    semi, phase0, gamma = fit_circular(bjd, rv, erv)

    # The run spans 1998 Jul 15-19 UT.  Note the lower bound must sit below
    # 2451009.94 or the first night is silently dropped.
    run = (bjd > 2451009.5) & (bjd < 2451015.)
    tt = np.linspace(bjd[run].min() - 0.7, bjd[run].max() + 0.7, 800)

    def model(times):
        return gamma + semi * np.cos(2 * np.pi * (times / PERIOD - phase0))

    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    ax.plot(tt, model(tt), color=C_MODEL, lw=2.,
            label='3.097-day orbit (K = {:.0f} m s$^{{-1}}$)'.format(semi))
    ax.errorbar(bjd[run], rv[run], yerr=erv[run], fmt='o', color=C_DATA,
                ms=8, lw=1.2, capsize=0, zorder=5, mec=SURFACE, mew=1.2,
                label='Keck/HIRES, {:d} observations'.format(run.sum()))
    ax.axhline(gamma, color='#c9c8c3', lw=0.8)

    # Label the nights themselves rather than days since some arbitrary zero
    nights = Time(['1998-07-{:02d}'.format(day) for day in range(15, 20)])
    ax.set_xticks(nights.jd + 0.5)
    ax.set_xticklabels(['15 July', '16 July', '17 July', '18 July', '19 July'])
    for edge in nights.jd:
        ax.axvline(edge, color='#e6e5e0', lw=0.8, zorder=0)

    ax.text(tt[-1], gamma + 46, 'moving towards us', color=C_DATA, fontsize=9,
            ha='right')
    ax.text(tt[-1], gamma - 54, 'moving away from us', color=C_DATA,
            fontsize=9, ha='right')

    ax.set_xlabel('Night of observation, 1998')
    ax.set_ylabel('Velocity of the star (m s$^{-1}$)')
    ax.set_xlim(tt[0], tt[-1])
    ax.grid(axis='y', alpha=0.6)
    ax.legend(loc='lower left')
    fig.tight_layout()
    _save(fig, 'fig4_july1998.png')


def fig_phasefold():
    """ All 30 discovery-era velocities, folded on the 3.097-day period. """
    bjd, rv, erv = load_velocities(discovery_era=True)
    semi, phase0, gamma = fit_circular(bjd, rv, erv)

    phase = np.mod(bjd / PERIOD - phase0, 1.)
    pp = np.linspace(0, 1, 400)

    fig, ax = plt.subplots(figsize=(9., 4.6))
    ax.plot(pp, gamma + semi * np.cos(2 * np.pi * pp), color=C_MODEL, lw=2.,
            label='circular orbit, K = {:.0f} m s$^{{-1}}$'.format(semi))
    for offset in (-1., 0., 1.):
        ax.errorbar(phase + offset, rv, yerr=erv, fmt='o', color=C_DATA, ms=8,
                    lw=1.2, capsize=0, zorder=5, mec=SURFACE, mew=1.2,
                    label='Keck/HIRES, 1997 Dec - 1998 Sep' if offset == 0 else None)
    ax.plot(pp - 1., gamma + semi * np.cos(2 * np.pi * pp), color=C_MODEL,
            lw=2., alpha=0.35)
    ax.plot(pp + 1., gamma + semi * np.cos(2 * np.pi * pp), color=C_MODEL,
            lw=2., alpha=0.35)
    ax.axhline(gamma, color='#c9c8c3', lw=0.8)

    ax.set_xlabel('Orbital phase')
    ax.set_ylabel('Velocity (m s$^{-1}$)')
    ax.set_xlim(-0.25, 1.25)
    ax.grid(alpha=0.6)
    ax.legend(loc='lower right')
    fig.tight_layout()
    _save(fig, 'fig5_phasefold.png')

    return semi, len(bjd)


# ---------------------------------------------------------------------------
# Figure 6 -- how close is 0.042 AU?
# ---------------------------------------------------------------------------

def fig_orbit_scale():
    """ The orbit of HD 187123 b against Mercury's. """
    fig, ax = plt.subplots(figsize=(7., 5.))
    theta = np.linspace(0, 2 * np.pi, 400)

    ax.plot(A_MERCURY * np.cos(theta), A_MERCURY * np.sin(theta),
            color=INK_2, lw=1.4, ls='--')
    ax.plot(A_PLANET * np.cos(theta), A_PLANET * np.sin(theta),
            color=C_MODEL, lw=2.)

    r_star = 0.00465   # one solar radius, in AU
    star = plt.Circle((0, 0), r_star * 3, color=C_DATA, zorder=5)
    ax.add_patch(star)
    ax.plot(A_PLANET, 0., 'o', color=C_MODEL, ms=9, zorder=6)

    ax.text(A_MERCURY * 0.72, A_MERCURY * 0.78, "Mercury's orbit", color=INK_2,
            fontsize=10)
    ax.text(0.06, 0.05, 'HD 187123 b\n0.042 AU', color=INK, fontsize=10)
    ax.text(0., -0.10, 'the star (drawn $\\times$3)', color=INK_2, fontsize=9,
            ha='center', va='top')

    ax.set_xlim(-0.45, 0.45)
    ax.set_ylim(-0.45, 0.45)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('A year that lasts three days', color=INK, loc='left')

    fig.tight_layout()
    _save(fig, 'fig6_orbit_scale.png')


# ---------------------------------------------------------------------------

def main():
    """ Generate every figure and print the numbers quoted in the document. """
    set_style()

    fig_wobble()
    fig_doppler()
    fig_iodine()
    fig_july1998()
    semi, nobs = fig_phasefold()
    fig_orbit_scale()

    bjd, rv, erv = load_velocities(discovery_era=True)
    tt = Time(bjd, format='jd')
    print('\n--- numbers quoted in docs/public_HD187123b.md ---')
    print('discovery-era epochs : {:d}'.format(nobs))
    print('first / last         : {:s} / {:s}'.format(tt[0].iso[:10],
                                                      tt[-1].iso[:10]))
    print('baseline             : {:.0f} d'.format(bjd[-1] - bjd[0]))
    print('fitted K             : {:.1f} m/s  (published 72)'.format(semi))
    print('median error         : {:.2f} m/s'.format(np.median(erv)))
    print('all epochs in file   : {:d}'.format(len(load_velocities()[0])))


if __name__ == '__main__':
    main()
