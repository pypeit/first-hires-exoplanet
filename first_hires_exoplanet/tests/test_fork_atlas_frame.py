""" An iodine atlas must say which wavelength convention it is on.

The FTS atlases `pyodine` ships carry *both* grids: `wavelength` (vacuum) and
`wavelength_air`.  Every one of upstream's four instrument adapters silently
reads the air grid, which is right for Lick and SONG and wrong for us: PypeIt
reports vacuum.  The two differ by 1.56 A at 5615 A, which is **83 km/s**, and
phase 2 measured the consequence — the correlation between the atlas and the
HIRES cell falls from +0.90 on the vacuum grid to +0.08 on the air one.

Nothing failed when that happened.  The atlas loaded, the model ran, and the
answer was wrong, which is the failure mode phase 3 was told to expect.  The
fix is not to change which grid is read — that would break Lick and SONG — but
to make the choice explicit and recorded, so that an atlas can be asked what
convention it is on and a mismatch can be caught rather than absorbed.

    conda run -n pypeit14 pytest first_hires_exoplanet/tests/test_fork_atlas_frame.py

"""

# Standard imports
import os
from pathlib import Path

import numpy as np
import pytest

#: The real atlas, if the external data tree is present
ATLAS = (Path(__file__).resolve().parents[2].parent / 'first-hires-exoplanet-data'
         / 'atlas' / 'Fischer_Cell_May2022_downsampled3.h5')


@pytest.fixture
def two_grid_atlas(tmp_path):
    """ A miniature FTS file with both wavelength grids, as the real ones have.

    Returns
    -------
    :obj:`pathlib.Path`
        Path to the written HDF5 file.
    """
    import h5py

    vacuum = np.linspace(5000., 5010., 1001)
    # The air grid is shorter by n-1, about 1.4 A here; the exact refractive
    # index does not matter, only that the two are distinguishable.
    air = vacuum / 1.000277
    flux = 1. - 0.5 * np.exp(-0.5 * ((vacuum - 5005.) / 0.05) ** 2)

    path = tmp_path / 'mini_atlas.h5'
    with h5py.File(path, 'w') as h:
        h.create_dataset('flux_normalized', data=flux)
        h.create_dataset('wavelength', data=vacuum)
        h.create_dataset('wavelength_air', data=air)
    return path


def test_the_frame_must_be_asked_for(two_grid_atlas):
    """ There is no default: a caller states the convention it works in. """
    from pyodine.components import IodineAtlas

    with pytest.raises(TypeError):
        IodineAtlas.from_h5(two_grid_atlas)


def test_vacuum_and_air_load_different_grids(two_grid_atlas):
    """ And the difference is the one that cost 83 km/s. """
    from pyodine.components import IodineAtlas

    vac = IodineAtlas.from_h5(two_grid_atlas, wave_frame='vacuum')
    air = IodineAtlas.from_h5(two_grid_atlas, wave_frame='air')

    assert vac.wave_frame == 'vacuum'
    assert air.wave_frame == 'air'
    assert np.allclose(vac.flux, air.flux)
    # Same lines, different labels: a shift of order an Angstrom at 5000 A.
    shift = np.median(vac.wave - air.wave)
    assert 1.0 < shift < 2.0


def test_an_unknown_frame_is_refused(two_grid_atlas):
    """ Not silently defaulted to one of them. """
    from pyodine.components import IodineAtlas

    with pytest.raises(ValueError):
        IodineAtlas.from_h5(two_grid_atlas, wave_frame='barycentric')


def test_a_missing_grid_is_refused(tmp_path):
    """ An atlas with only one grid cannot be read as the other. """
    import h5py
    from pyodine.components import IodineAtlas

    path = tmp_path / 'air_only.h5'
    with h5py.File(path, 'w') as h:
        h.create_dataset('flux_normalized', data=np.ones(10))
        h.create_dataset('wavelength_air', data=np.linspace(5000., 5001., 10))

    assert IodineAtlas.from_h5(path, wave_frame='air').wave_frame == 'air'
    with pytest.raises(KeyError):
        IodineAtlas.from_h5(path, wave_frame='vacuum')


def test_a_hand_built_atlas_admits_it_does_not_know():
    """ Constructing one directly leaves the convention unrecorded, not wrong.

    That is the honest default: it is what upstream's adapters effectively do,
    and it is what `require_wave_frame` refuses to accept.
    """
    from pyodine.components import IodineAtlas

    atlas = IodineAtlas(np.ones(10), wave=np.linspace(5000., 5001., 10))
    assert atlas.wave_frame is None

    atlas.require_wave_frame('vacuum')  # unknown: warns, does not raise


def test_a_frame_mismatch_is_caught(two_grid_atlas):
    """ The check that would have caught phase 2's 83 km/s. """
    from pyodine.components import DataMismatchError, IodineAtlas

    air = IodineAtlas.from_h5(two_grid_atlas, wave_frame='air')
    air.require_wave_frame('air')
    with pytest.raises(DataMismatchError):
        air.require_wave_frame('vacuum')


@pytest.mark.skipif(not ATLAS.exists(), reason='atlas not present on disk')
def test_the_real_atlas_shift_is_the_measured_one():
    """ 1.56 A at 5615 A, which is 83 km/s. """
    from pyodine.components import IodineAtlas

    vac = IodineAtlas.from_h5(ATLAS, wave_frame='vacuum')
    air = IodineAtlas.from_h5(ATLAS, wave_frame='air')

    here = np.argmin(np.abs(vac.wave - 5615.))
    shift = vac.wave[here] - air.wave[here]
    assert np.isclose(shift, 1.56, atol=0.02)
    assert np.isclose(2.998e5 * shift / 5615., 83., atol=2.)
