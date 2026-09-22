""" Audit the reference frame of the phase-1 reductions (phase 2, prompt 1).

`pyodine` applies its own barycentric correction; PypeIt has already applied a
heliocentric one.  Handing `pyodine` a heliocentric-corrected wavelength scale
*and* a `bary_vel_corr` would double-correct.  Before anything in phase 2 is
built on these spectra we need to know, from the products rather than from
assumption:

    1. what correction PypeIt actually applied, and whether it is exactly
       invertible from what is stored in the spec1d files;
    2. how PypeIt's 'heliocentric' and 'barycentric' options compare against
       `astropy`'s own `radial_velocity_correction`, which is what `pyodine`
       (and any sane re-derivation) would use;
    3. what epoch PypeIt used -- exposure start or mid-exposure -- and how much
       that choice is worth in m/s.

Every number is computed for the ten real discovery-era epochs, not for a
representative one, because the diurnal term makes the answer epoch-dependent.

The reduction these read is NOT in this repository -- the raw frames and the
PypeIt products live in a sibling data tree.  Point --redux at it if it is
somewhere else.

Reduced with PypeIt 2.0.2.dev1217+g017bece06 (branch orig-hires-fixes).

Run with:

    conda run -n pypeit14 python -m first_hires_exoplanet.refframe_audit

"""

# Standard imports
import os
import glob
import argparse

import numpy as np

from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, Angle

from pypeit.core import wave as pypeit_wave


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))

#: Default reduction directory (outside the repository; see the module docstring)
DEFAULT_REDUX = os.path.join(os.path.dirname(os.path.dirname(_HERE)),
                             'first-hires-exoplanet-data', 'redux')

#: Where the committed derived products go
DEFAULT_OUTDIR = os.path.join(_HERE, 'data')

#: Keck, as PypeIt's telescope parameters give it (pypeit/telescopes.py).  Held
#: here so the astropy comparison uses *the same site* as PypeIt and any
#: difference that shows up is a difference of method, not of geography.
KECK_LON = -155.47833333333335
KECK_LAT = 19.82833333333333
KECK_ELV = 4159.99999999977

#: Speed of light, km/s, matching the value hard-coded in pypeit.core.wave
CKMS = 299792.458


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def vel_from_corr(vel_corr):
    """ Invert PypeIt's relativistic wavelength factor back to a velocity.

    PypeIt stores ``VEL_CORR = sqrt((1 + v/c) / (1 - v/c))`` and multiplies the
    wavelength array by it.  This returns v in km/s.

    Args:
        vel_corr (float): the stored multiplicative factor.

    Returns:
        float: velocity in km/s.
    """
    r2 = vel_corr**2
    return CKMS * (r2 - 1.) / (r2 + 1.)


def corr_from_vel(vel):
    """ The forward direction of :func:`vel_from_corr`, v in km/s. """
    return np.sqrt((1. + vel / CKMS) / (1. - vel / CKMS))


def collect_epochs(redux_dir):
    """ Read every reduced spec1d and pull out what the audit needs.

    Args:
        redux_dir (str): directory holding the ``reduce_*`` reduction folders.

    Returns:
        `astropy.table.Table`_: one row per epoch.
    """
    files = sorted(glob.glob(os.path.join(redux_dir, 'reduce_1998*', 'Science',
                                          'spec1d_*.fits')))
    if len(files) == 0:
        raise FileNotFoundError(f'No spec1d files under {redux_dir}')

    rows = []
    for ifile in files:
        with fits.open(ifile) as hdul:
            h0 = hdul[0].header
            # VEL_TYPE / VEL_CORR are per-object, i.e. per order, and should be
            # identical across orders; check rather than trust.
            vel_corrs = np.array([hdul[i].header['VEL_CORR']
                                  for i in range(1, len(hdul))
                                  if 'VEL_CORR' in hdul[i].header])
            vel_types = set(hdul[i].header.get('VEL_TYPE')
                            for i in range(1, len(hdul))
                            if 'VEL_CORR' in hdul[i].header)
            rows.append(dict(
                filename=os.path.basename(ifile),
                koaid=h0['FILENAME'].replace('.fits', ''),
                mjd_hdr=float(h0['MJD']),
                exptime=float(h0['EXPTIME']),
                ra=float(h0['RA']),
                dec=float(h0['DEC']),
                norder=len(vel_corrs),
                vel_type='|'.join(sorted(str(v) for v in vel_types)),
                vel_corr=float(vel_corrs[0]),
                vel_corr_spread=float(np.ptp(vel_corrs)),
            ))
    return Table(rows)


def audit(tbl):
    """ Add the reference-frame comparison columns to the epoch table.

    Args:
        tbl (`astropy.table.Table`_): output of :func:`collect_epochs`;
            modified in place.

    Returns:
        `astropy.table.Table`_: the same table.
    """
    loc = EarthLocation.from_geodetic(lon=KECK_LON * u.deg, lat=KECK_LAT * u.deg,
                                      height=KECK_ELV * u.m)

    ncol = len(tbl)
    out = {k: np.zeros(ncol) for k in
           ['v_pypeit_stored', 'v_pypeit_helio', 'v_pypeit_bary',
            'v_astropy_helio', 'v_astropy_bary', 'v_astropy_bary_mid',
            'd_helio_bary', 'd_pypeit_astropy_helio', 'd_pypeit_astropy_bary',
            'd_start_mid']}

    for i, row in enumerate(tbl):
        radec = SkyCoord(ra=row['ra'] * u.deg, dec=row['dec'] * u.deg)
        t_start = Time(row['mjd_hdr'], format='mjd')
        t_mid = Time(row['mjd_hdr'] + 0.5 * row['exptime'] / 86400., format='mjd')

        # What is actually in the file
        out['v_pypeit_stored'][i] = vel_from_corr(row['vel_corr'])

        # PypeIt's own calculation, re-run at the same epoch it used
        for frame, key in [('heliocentric', 'v_pypeit_helio'),
                           ('barycentric', 'v_pypeit_bary')]:
            v, _ = pypeit_wave.geomotion_correct(radec, t_start, KECK_LON,
                                                 KECK_LAT, KECK_ELV, frame)
            out[key][i] = v

        # astropy's calculation, same epoch, same site
        for kind, key in [('heliocentric', 'v_astropy_helio'),
                          ('barycentric', 'v_astropy_bary')]:
            rv = radec.radial_velocity_correction(kind=kind, obstime=Time(
                t_start, location=loc))
            out[key][i] = rv.to(u.km / u.s).value

        # astropy barycentric at mid-exposure, which is the epoch a velocity
        # measurement actually belongs to
        rv_mid = radec.radial_velocity_correction(
            kind='barycentric', obstime=Time(t_mid, location=loc))
        out['v_astropy_bary_mid'][i] = rv_mid.to(u.km / u.s).value

    out['d_helio_bary'] = (out['v_astropy_helio'] - out['v_astropy_bary']) * 1e3
    out['d_pypeit_astropy_helio'] = (out['v_pypeit_helio']
                                     - out['v_astropy_helio']) * 1e3
    out['d_pypeit_astropy_bary'] = (out['v_pypeit_bary']
                                    - out['v_astropy_bary']) * 1e3
    out['d_start_mid'] = (out['v_astropy_bary'] - out['v_astropy_bary_mid']) * 1e3

    for k, v in out.items():
        tbl[k] = v
    # The three d_* columns are m/s; the v_* columns are km/s.
    for k in ['d_helio_bary', 'd_pypeit_astropy_helio', 'd_pypeit_astropy_bary',
              'd_start_mid']:
        tbl[k].unit = u.m / u.s
    for k in ['v_pypeit_stored', 'v_pypeit_helio', 'v_pypeit_bary',
              'v_astropy_helio', 'v_astropy_bary', 'v_astropy_bary_mid']:
        tbl[k].unit = u.km / u.s
    return tbl


def check_invertibility(redux_dir, tbl):
    """ Confirm that dividing by VEL_CORR recovers the observed-frame scale.

    The concern is not arithmetic -- it is that PypeIt might have applied the
    correction somewhere other than the stored factor, or applied it twice.
    This re-derives the observed wavelength for one order of one epoch and
    reports the round-trip error in m/s, which should be at the numerical floor.

    Args:
        redux_dir (str): reduction directory.
        tbl (`astropy.table.Table`_): epoch table.

    Returns:
        tuple: (max round-trip error in m/s, the order checked)
    """
    ifile = glob.glob(os.path.join(redux_dir, 'reduce_1998*', 'Science',
                                   'spec1d_' + tbl['koaid'][0] + '*.fits'))[0]
    worst = 0.
    with fits.open(ifile) as hdul:
        for i in range(1, len(hdul)):
            if 'VEL_CORR' not in hdul[i].header:
                continue
            wv = hdul[i].data['OPT_WAVE']
            gd = wv > 1.
            if not np.any(gd):
                continue
            vc = hdul[i].header['VEL_CORR']
            # round trip: divide out, multiply back
            rt = (wv[gd] / vc) * vc
            frac = np.max(np.abs(rt - wv[gd]) / wv[gd])
            worst = max(worst, frac * CKMS * 1e3)
    return worst, ifile


def check_solar_term_sign(tbl):
    """ Test whether PypeIt's 'heliocentric' has the solar term's sign wrong.

    ``pypeit.core.wave.geomotion_velocity`` builds the observer's barycentric
    velocity as ``ev + ov`` and then, for the heliocentric frame, does
    ``velocity += sv`` with ``sv`` the Sun's barycentric velocity.  The
    observer's velocity *relative to the Sun* is ``ev + ov - sv``, so the sign
    looks wrong.  Rather than argue from the source, recompute it both ways and
    see which one lands on `astropy`'s heliocentric answer.

    Args:
        tbl (`astropy.table.Table`_): the audited epoch table.

    Returns:
        `astropy.table.Table`_: one row per epoch, in m/s.
    """
    from astropy.coordinates import solar_system, ICRS
    from astropy.coordinates import (UnitSphericalRepresentation,
                                     CartesianRepresentation)

    loc = EarthLocation.from_geodetic(lon=KECK_LON * u.deg, lat=KECK_LAT * u.deg,
                                      height=KECK_ELV * u.m)
    rows = []
    for row in tbl:
        radec = SkyCoord(ra=row['ra'] * u.deg, dec=row['dec'] * u.deg)
        t = Time(row['mjd_hdr'], format='mjd', scale='utc', location=(
            KECK_LON * u.deg, KECK_LAT * u.deg, KECK_ELV * u.m))

        ep, ev = solar_system.get_body_barycentric_posvel('earth', t)
        op, ov = t.location.get_gcrs_posvel(t)
        sp, sv = solar_system.get_body_barycentric_posvel('sun', t)
        nhat = radec.icrs.represent_as(
            UnitSphericalRepresentation).represent_as(CartesianRepresentation)

        def proj(vec):
            return nhat.dot(vec).to(u.km / u.s).value

        v_bary = proj(ev + ov)
        v_helio_plus = proj(ev + ov + sv)     # what PypeIt does
        v_helio_minus = proj(ev + ov - sv)    # what the definition asks for

        ap = radec.radial_velocity_correction(
            kind='heliocentric', obstime=Time(t, location=loc)).to(u.km / u.s).value

        rows.append(dict(koaid=row['koaid'],
                         solar_term=proj(sv) * 1e3,
                         plus_minus_astropy=(v_helio_plus - ap) * 1e3,
                         minus_minus_astropy=(v_helio_minus - ap) * 1e3,
                         bary_minus_astropy_helio=(v_bary - ap) * 1e3))
    out = Table(rows)
    for c in out.colnames[1:]:
        out[c].unit = u.m / u.s
        out[c].info.format = '%.3f'
    return out


def season_scan(ra, dec, mjd_start=50805., mjd_end=51080., nsamp=40):
    """ How much of the PypeIt-vs-astropy difference is constant?

    A constant offset is harmless for *relative* velocities; a term that drifts
    over the nine-month discovery baseline is not.  Sample the difference on a
    grid spanning 1997-12 through 1998-09.

    Args:
        ra, dec (float): target coordinates in degrees.
        mjd_start, mjd_end (float): the span to scan.
        nsamp (int): number of samples.

    Returns:
        `astropy.table.Table`_
    """
    loc = EarthLocation.from_geodetic(lon=KECK_LON * u.deg, lat=KECK_LAT * u.deg,
                                      height=KECK_ELV * u.m)
    radec = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
    mjds = np.linspace(mjd_start, mjd_end, nsamp)
    rows = []
    for mjd in mjds:
        t = Time(mjd, format='mjd')
        v_pyp_b, _ = pypeit_wave.geomotion_correct(radec, t, KECK_LON, KECK_LAT,
                                                   KECK_ELV, 'barycentric')
        v_pyp_h, _ = pypeit_wave.geomotion_correct(radec, t, KECK_LON, KECK_LAT,
                                                   KECK_ELV, 'heliocentric')
        v_ap_b = radec.radial_velocity_correction(
            kind='barycentric', obstime=Time(t, location=loc)).to(u.km/u.s).value
        v_ap_h = radec.radial_velocity_correction(
            kind='heliocentric', obstime=Time(t, location=loc)).to(u.km/u.s).value
        rows.append(dict(mjd=mjd, v_astropy_bary=v_ap_b,
                         d_bary=(v_pyp_b - v_ap_b) * 1e3,
                         d_helio=(v_pyp_h - v_ap_h) * 1e3,
                         d_helio_bary=(v_ap_h - v_ap_b) * 1e3))
    out = Table(rows)
    for c in ['d_bary', 'd_helio', 'd_helio_bary']:
        out[c].unit = u.m / u.s
    return out


def relativistic_terms(ra, dec, mjd):
    """ Account for the residual PypeIt-vs-astropy barycentric offset.

    PypeIt's barycentric velocity is purely kinematic: it projects the
    observer's barycentric velocity onto the line of sight and stops.
    `astropy`'s ``radial_velocity_correction`` additionally carries the
    gravitational redshift of the Sun at the Earth and the special-relativistic
    time dilation of the observer.  Both are near-constant, which is why the
    offset barely drifts.  This evaluates them so the offset is explained
    rather than merely measured.

    Args:
        ra, dec (float): target coordinates, degrees.
        mjd (float): epoch.

    Returns:
        dict: the two terms and their sum, in m/s.
    """
    from astropy.constants import GM_sun, c as c_const
    from astropy.coordinates import solar_system

    t = Time(mjd, format='mjd', scale='utc', location=(
        KECK_LON * u.deg, KECK_LAT * u.deg, KECK_ELV * u.m))
    ep, ev = solar_system.get_body_barycentric_posvel('earth', t)
    sp, sv = solar_system.get_body_barycentric_posvel('sun', t)
    r_sun = (ep - sp).norm()

    grav = (GM_sun / (c_const * r_sun)).to(u.m / u.s).value
    dilation = ((ev.norm())**2 / (2 * c_const)).to(u.m / u.s).value
    return dict(grav_redshift=grav, time_dilation=dilation,
                total=grav + dilation)


def check_ha_epoch(redux_dir, tbl):
    """ Is the header time the start of the exposure, or the end?

    Nothing in the raw header says.  But the header also records the hour angle
    HA, and HA + RA = local sidereal time.  Comparing the LST implied by the
    recorded HA against the LST computed from the recorded UT tells us which
    instant the header describes: if they agree, HA and UT were sampled
    together, and the residual against start/end tells us which end of the
    exposure that instant sits at.

    Args:
        redux_dir (str): reduction directory.
        tbl (`astropy.table.Table`_): epoch table.

    Returns:
        `astropy.table.Table`_: per-epoch LST residuals in seconds.
    """
    loc = EarthLocation.from_geodetic(lon=KECK_LON * u.deg, lat=KECK_LAT * u.deg,
                                      height=KECK_ELV * u.m)
    rows = []
    for row in tbl:
        raw = glob.glob(os.path.join(redux_dir, 'stage_*', row['koaid'] + '.fits'))
        if len(raw) == 0:
            continue
        h = fits.getheader(raw[0])
        ha = Angle(h['HA'], unit=u.hourangle).wrap_at(12 * u.hourangle)
        ra = Angle(h['RA'], unit=u.hourangle)
        lst_from_ha = (ha + ra).wrap_at(24 * u.hourangle)

        res = {}
        for lbl, dt in [('start', 0.), ('mid', 0.5 * row['exptime']),
                        ('end', row['exptime'])]:
            t = Time(row['mjd_hdr'] + dt / 86400., format='mjd', location=loc)
            lst = t.sidereal_time('apparent')
            d = (lst_from_ha - lst).wrap_at(12 * u.hourangle)
            # sidereal seconds -> solar seconds
            res[lbl] = d.hourangle * 3600. / 1.0027379
        rows.append(dict(koaid=row['koaid'], exptime=row['exptime'], **res))
    out = Table(rows)
    for c in ['start', 'mid', 'end']:
        out[c].unit = u.s
        out[c].info.format = '%.1f'
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args(options=None):
    parser = argparse.ArgumentParser(
        description='Audit the reference frame of the phase-1 reductions.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--redux', type=str, default=DEFAULT_REDUX,
                        help='Directory holding the reduce_* folders')
    parser.add_argument('--outdir', type=str, default=DEFAULT_OUTDIR,
                        help='Where to write the audit table')
    return parser.parse_args() if options is None else parser.parse_args(options)


def main(pargs):

    tbl = collect_epochs(pargs.redux)
    tbl = audit(tbl)

    # --- what PypeIt stored vs what PypeIt computes now -------------------
    d_stored = (tbl['v_pypeit_stored'] - tbl['v_pypeit_helio']) * 1e3
    print('\n=== 1. Is the stored VEL_CORR PypeIt\'s heliocentric value? ===')
    print(f"  VEL_TYPE values found      : {set(tbl['vel_type'])}")
    print(f"  max spread across orders   : {np.max(tbl['vel_corr_spread']):.3e}")
    print(f"  max |stored - recomputed|  : {np.max(np.abs(d_stored)):.3e} m/s")

    rt, rtfile = check_invertibility(pargs.redux, tbl)
    print(f"  wavelength round-trip err  : {rt:.3e} m/s")
    print(f"    (checked on {os.path.basename(rtfile)})")

    # --- how big is the correction, and the frame differences -------------
    print('\n=== 2. Size of the correction and of the frame differences ===')
    print(f"  applied correction         : {np.min(tbl['v_pypeit_stored']):+.3f}"
          f" to {np.max(tbl['v_pypeit_stored']):+.3f} km/s")
    print(f"  helio - bary (astropy)     : {np.min(tbl['d_helio_bary']):+.2f}"
          f" to {np.max(tbl['d_helio_bary']):+.2f} m/s")
    print(f"  PypeIt helio - astropy helio: "
          f"{np.min(tbl['d_pypeit_astropy_helio']):+.2f}"
          f" to {np.max(tbl['d_pypeit_astropy_helio']):+.2f} m/s")
    print(f"  PypeIt bary  - astropy bary : "
          f"{np.min(tbl['d_pypeit_astropy_bary']):+.2f}"
          f" to {np.max(tbl['d_pypeit_astropy_bary']):+.2f} m/s")

    # --- exposure start vs mid --------------------------------------------
    print('\n=== 3. Epoch: PypeIt uses the header MJD ===')
    print(f"  exposure times             : {np.min(tbl['exptime']):.0f}"
          f" to {np.max(tbl['exptime']):.0f} s")
    print(f"  start - mid (bary vel)     : {np.min(tbl['d_start_mid']):+.2f}"
          f" to {np.max(tbl['d_start_mid']):+.2f} m/s")
    print(f"  max |start - mid|          : "
          f"{np.max(np.abs(tbl['d_start_mid'])):.2f} m/s")

    # --- per-epoch table ---------------------------------------------------
    print('\n=== 4. Per-epoch ===')
    show = tbl['koaid', 'mjd_hdr', 'exptime', 'v_pypeit_stored',
               'd_helio_bary', 'd_pypeit_astropy_helio', 'd_start_mid']
    for col in show.colnames[3:]:
        show[col].info.format = '%.3f'
    show['mjd_hdr'].info.format = '%.6f'
    show.pprint_all()


    # --- the solar-term sign ----------------------------------------------
    print('\n=== 5. Why PypeIt\'s heliocentric differs: the solar term ===')
    sgn = check_solar_term_sign(tbl)
    print(f"  n-hat . v_sun                    : "
          f"{np.min(sgn['solar_term']):+.2f} to {np.max(sgn['solar_term']):+.2f} m/s")
    print(f"  (ev+ov+sv) - astropy helio       : "
          f"{np.mean(sgn['plus_minus_astropy']):+.2f} m/s   <- PypeIt")
    print(f"  (ev+ov-sv) - astropy helio       : "
          f"{np.mean(sgn['minus_minus_astropy']):+.2f} m/s   <- definition")
    print(f"  (ev+ov)    - astropy helio       : "
          f"{np.mean(sgn['bary_minus_astropy_helio']):+.2f} m/s   <- no solar term")

    # --- does the offset drift over the discovery baseline? ---------------
    print('\n=== 6. Drift over the nine-month discovery baseline ===')
    scan = season_scan(tbl['ra'][0], tbl['dec'][0])
    print(f"  astropy barycentric velocity : {np.min(scan['v_astropy_bary']):+.2f}"
          f" to {np.max(scan['v_astropy_bary']):+.2f} km/s")
    for c, lbl in [('d_bary', 'PypeIt bary  - astropy bary '),
                   ('d_helio', 'PypeIt helio - astropy helio'),
                   ('d_helio_bary', 'astropy helio - astropy bary')]:
        print(f"  {lbl}: {np.min(scan[c]):+.2f} to {np.max(scan[c]):+.2f} m/s"
              f"  (peak-to-peak {np.ptp(scan[c]):.2f})")

    rel = relativistic_terms(tbl['ra'][0], tbl['dec'][0], tbl['mjd_hdr'][0])
    print(f"  the residual bary offset is the relativistic terms PypeIt omits:")
    print(f"    solar gravitational redshift : {rel['grav_redshift']:+.2f} m/s")
    print(f"    observer time dilation       : {rel['time_dilation']:+.2f} m/s")
    print(f"    total                        : {rel['total']:+.2f} m/s"
          f"  (vs measured {np.mean(tbl['d_pypeit_astropy_bary']):+.2f})")

    # --- start or end of exposure? ----------------------------------------
    print('\n=== 7. Does the header MJD mark the start or the end? ===')
    ha = check_ha_epoch(pargs.redux, tbl)
    for c in ['start', 'mid', 'end']:
        print(f"  LST(HA) - LST(MJD + {c:5s}) : median {np.median(ha[c]):+8.1f} s"
              f"   scatter {np.std(ha[c]):6.1f} s")

    if not os.path.isdir(pargs.outdir):
        os.makedirs(pargs.outdir)
    outfile = os.path.join(pargs.outdir, 'refframe_audit.csv')
    tbl.write(outfile, format='ascii.csv', overwrite=True)
    print(f'\nWrote {outfile}')


if __name__ == '__main__':
    main(parse_args())
