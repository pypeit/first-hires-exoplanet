from pyodine.components import Instrument

import os

i2_dir_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../iodine_atlas')

# Used by get_instrument below...  # TODO: Make a settings file
my_instruments = {
    'song_1': Instrument(
            'SONG Hertzsprung spectrograph (Tenerife)',
            latitude=28.2983,
            longitude=-16.5094,  # East longitude
            altitude=2400.0
    ),
    'song_2': Instrument(
            'SONG China spectrograph (Delingha)',
            latitude=37.378001,
            longitude=97.73167,  # East longitude
            altitude=3200.    # From https://arxiv.org/pdf/1602.00838
    ),
    'lick': Instrument(
            'Hamilton spectrograph (Lick observatory)',
            latitude=37.34139,
            longitude=238.35722,
            altitude=1283.
    ),
    'waltz': Instrument(
            'Waltz spectrograph (LSW Heidelberg)',
            latitude=49.398611,
            longitude=8.720833,
            altitude=560.
    )
}

# List of iodine atlas locations
my_iodine_atlases = {
    1: os.path.join(i2_dir_path, 'song_iodine_cell_01_65C.h5'),
    2: os.path.join(i2_dir_path, 'ftslick05_norm_new.h5')
}
