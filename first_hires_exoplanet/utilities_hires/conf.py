""" Instruments and iodine atlases for Keck/HIRES (phase 3, prompt 4).

Modelled on `utilities_lick/conf.py` of `pyodine`
4488b0914fe5b272b787982647691045bff2604a, which the drivers read for two
things: an `Instrument` per observatory, and ``my_iodine_atlases``, an integer
index -> atlas path that ``Parameters.i2_to_use`` selects from.

Two differences from Lick, both deliberate:

* Upstream resolves ``../iodine_atlas`` relative to its own directory.  This
  package does not live inside the vendored tree, and the atlas is the one
  phase 2 measured against, so the path is the adapter's ``DEFAULT_ATLAS``.
  It is byte-identical to the copy `pyodine` ships (phase 3 prompt 1,
  SHA-256 ``3ae788e7...a55a429``).
* Each atlas records the wavelength grid it is to be read on.  Upstream's
  atlases are read on air; PypeIt reports vacuum, and the two differ by
  83 km/s at 5615 A (phase 2 prompt 7).  ``my_iodine_atlas_frames`` is
  consumed by ``load_pyodine.IodineTemplate`` and checked against the fork's
  ``IodineAtlas.require_wave_frame`` (phase 3 prompt 2, change 2).
"""

from pyodine.components import Instrument

from .load_pyodine import DEFAULT_ATLAS, KECK_FALLBACK


#: Keck I, from PypeIt's telescope parameters (also LON-OBS / LAT-OBS /
#: ALT-OBS in every spec1d, which is what the adapter actually reads).
#: `pyodine` wants east longitude; Mauna Kea is west, so this is negative.
my_instruments = {
    'keck_hires': Instrument(
            'HIRES (Keck I, Mauna Kea)',
            latitude=KECK_FALLBACK['lat'],
            longitude=KECK_FALLBACK['lon'],
            altitude=KECK_FALLBACK['elv']
    ),
}

#: Iodine atlases, by the index ``Parameters.i2_to_use`` selects.  There is
#: one: the Fischer May 2022 scan of *a* cell, 4980-6250 A.  No FTS scan of
#: the Keck/HIRES cell itself has been found; the adapter turns this one into
#: the HIRES cell with a Beer-Lambert exponent (``load_pyodine.ATLAS_ALPHA``).
my_iodine_atlases = {
    1: DEFAULT_ATLAS,
}

#: The wavelength grid each atlas above is to be read on.  Vacuum, because
#: PypeIt's wavelengths are vacuum.
my_iodine_atlas_frames = {
    1: 'vacuum',
}
