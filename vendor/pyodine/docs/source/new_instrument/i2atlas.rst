.. _new_i2atlas:

Add a new I2 atlas
==================

The basic principle of the I2 cell method is to use an extremely high resolved atlas of the I2 cell absorption features as reference spectrum in the modelling of observations. To reach high RV precision, you want the atlas to be a measurement of exactly that I2 cell which is also used in your observations, as the exact shapes of the absorption features depend on parameters such as the temperature and pressure of the I2 gas in your cell. With the **pyodine** package we are distributing the I2 cell atlases used in the SONG project and the Hamilton spectrograph at Lick observatory, but if you want to tailor the software to a new instrument, you should add the I2 atlas of your own cell!

What features should the I2 atlas have?
---------------------------------------

The I2 atlas should be measured at very high S/N and high resolving power - e.g., both the Lick and SONG atlas have a resolving power around :math:`R = \frac{\lambda}{\Delta\lambda} \sim 10^6`, which was achieved by recording their transmission spectra with a Fourier Transform Spectrometer (FTS).

.. _new_i2atlas_format:

Which format is required by **pyodine**?
----------------------------------------

The I2 atlases in **pyodine** come as HDF5 files, which can be accessed through Python using the package `h5py <http://www.h5py.org/>`_. Check out our *Tutorial* chapter on how to :doc:`../tutorial/preparation` to get a glance into the structure of the SONG I2 atlas. But basically you need to store at least two data arrays in the file:

* A wavelength vector, provided in Angstrom. Both the SONG and Lick atlas come with vacuum AND air wavelengths (denoted with keys *wavelength* and *wavelength_air* in the files), but for both instruments only the latter is used in the modelling (as the wavelength calibration in the extraction of the Echelle spectra uses air wavelengths). For your own instrument, you can include whichever you want to work in (and also use `common transformations between both <https://pyastronomy.readthedocs.io/en/latest/pyaslDoc/aslDoc/pyasl_wvlconv.html>`_ if you want to change it later on).

* A flux vector, preferrably already normalized, i.e. with values between 0 (complete absorption) and 1 (no absorption). Both the SONG and Lick atlas come with the original measured flux values (key *flux*) and the normalized flux (key *flux_normalized*).

Of course, you can also use different key names for the data arrays, or even save the atlas in a different file format altogether (e.g. FITS), but then you will have to change how to read the information from it (see :ref:`new_utilities_load_pyodine`).


