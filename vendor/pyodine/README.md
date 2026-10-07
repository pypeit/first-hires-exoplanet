# pyodine

[![License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](https://gitlab.com/Heeren/pyodine/-/blob/master/LICENSE) [![Documentation Status](https://readthedocs.org/projects/pyodine/badge/?version=latest)](https://pyodine.readthedocs.io/en/latest/?badge=latest)

An I2 analysis code for the determination of precise radial velocities from extracted input spectra, using the I2 cell method developed by [Butler et al. (1996)](https://ui.adsabs.harvard.edu/abs/1996PASP..108..500B/abstract). The code is open-source and build in a flexible and modular approach to allow easy adaptation to different instruments; it has been tested successfully on Lick and SONG spectra thus far.

## Documentation

Please refer to [pyodine.readthedocs.io](http://pyodine.readthedocs.io/).

## Developer

Paul Heeren, LSW Heidelberg; René Tronsgaard-Rasmussen, DTU Space; Frank Grundahl, SAC Aarhus

## Contributer

People who have actively contributed to the code:

Ayk Jessen, Master student @ LSW Heidelberg

## Attribution

Tbd - we are in the process of publishing first results from the code.

## Acknowledgements

We wish to thank Paul Butler, Debra Fischer, Geoff Marcy and Sharon Wang for many useful conversations and inputs on radial velocity extraction from iodine based data.

In our software, we use spectral atlases of the Sun and Arcturus as reference spectra for a first velocity guess. This data has been collected and published by [Hinkle et al. (2000)](https://ui.adsabs.harvard.edu/abs/2000vnia.book.....H/abstract), and can be downloaded for free from the [Astro Data Lab at NSF’s NOIRLab](https://noirlab.edu/science/data-services/other), which is operated by the Association of Universities for Research in Astronomy (AURA), Inc. under a cooperative agreement with the National Science Foundation.

## License

Copyright 2021 Paul Heeren, René Tronsgaard-Rasmussen, Frank Grundahl

pyodine is free software made available under the MIT License. For details see the LICENSE file.
