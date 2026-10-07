.. _overview_utilities:

Overview of the ``utilities`` modules
=====================================

Within the ``utilities`` modules we have defined all instrument-specific code and parameters, such as routines on how to load observation spectra of a specific instrument, or the relevant model parameters (LSF shape, chunk sizes etc.) employed in the fitting process. So far, three such modules should be present in the top level of the **pyodine** repository: ``utilities_song``, ``utilities_lick`` and ``utilities_waltz`` for the SONG, Lick and Waltz instruments, respectively (with the first two having been tested successfully). All these modules contain a few must-have routines, which we will explain in the following.

.. _overview_utilities_conf:

``conf.py``
-----------

The module ``conf.py`` mainly consists of two dictionaries with required meta-information on the instruments:

1. Geographical data, i.e. altitude, longitude and latitude of the observatories, which is used to compute the barycentric correction for the radial velocities - here an exemplary entry for SONG on Tenerife:
::

    my_instruments = {
            #other instruments,
            'song_1': Instrument(
                    'SONG Hertzsprung spectrograph (Tenerife)',
                    latitude=28.2983,
                    longitude=-16.5094,  # East longitude
                    altitude=2400.0
            )
    }

with the class :class:`pyodine.components.Instrument`.

2. Paths to the I2 atlases of the instruments (with simple indices as keys) - again, exemplary for SONG Tenerife:
::

    my_iodine_atlases = {
            #other I2 atlases,
            1: os.path.join(i2_dir_path, 'song_iodine_cell_01_65C.h5')
    }

where ``i2_dir_path`` points to the directory ``iodine_atlas`` within the **pyodine** repository.

.. _overview_utilities_load_pyodine:

``load_pyodine.py``
-------------------

The file ``load_pyodine.py`` contains all the necessary code to correctly read data products of the instrument in question:

1. The class :class:`load_pyodine.IodineTemplate` specifies how to load the I2 atlas data, which depends on the file format that the data is stored in (more on that in :ref:`new_i2atlas_format`):

.. autoclass:: utilities_song.load_pyodine.IodineTemplate

2. The class :class:`load_pyodine.ObservationWrapper` is used to store the observations acquired with your instrument in memory and make them accessible for **pyodine** - i.e. both the actual spectral data of the (already extracted and wavelength calibrated) Echelle orders as well as supplementary information such as date and object name of an observation:

.. autoclass:: utilities_song.load_pyodine.ObservationWrapper

.. _overview_utilities_pyodine_parameters:

``pyodine_parameters.py``
-------------------------

All instrument-specific parameters for the modelling of the observations, such as which Echelle orders to fit, LSF models to use, or in how many runs the fitting should be performed and what the starting values (and bounds) of the fit parameters should be, are specified in the :class:`pyodine_parameters.Parameters` input object for the instrument - this is then handed to the main routine :func:`pyodine_model_observations.model_single_observation()`:

.. autoclass:: utilities_song.pyodine_parameters.Parameters
   :inherited-members:

Very similarly, a :class:`pyodine_parameters.Template_Parameters` class exists for the creation of the deconvolved stellar templates (main routine :func:`pyodine_create_templates.create_template`):

.. autoclass:: utilities_song.pyodine_parameters.Template_Parameters

.. _overview_utilities_timeseries_parameters:

``timeseries_parameters.py``
----------------------------

Finally, parameters required in the combination algorithm for the individual chunk velocities :func:`pyodine_combine_vels.combine_velocity_results` (to arrive at a RV timeseries for a star) are defined in the class :class:`timeseries_parameters.Timeseries_Parameters`:

.. autoclass:: utilities_song.timeseries_parameters.Timeseries_Parameters
