
# Standard imports
import glob, os
from setuptools import setup, find_packages


# Begin setup
setup_keywords = dict()
setup_keywords['name'] = 'first_hires_exoplanet'
setup_keywords['description'] = 'Reducing and analyzing the first Keck/HIRES exoplanet spectra with PypeIt'
setup_keywords['author'] = 'J. Xavier Prochaska'
setup_keywords['author_email'] = 'jxp@ucsc.edu'
setup_keywords['license'] = 'BSD'
setup_keywords['url'] = 'https://github.com/pypeit/first-hires-exoplanet'
setup_keywords['version'] = '0.0.dev0'
# Use README.md as long_description.
setup_keywords['long_description'] = ''
if os.path.exists('README.md'):
    with open('README.md') as readme:
        setup_keywords['long_description'] = readme.read()
setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['python_requires'] = '>=3.11'
setup_keywords['install_requires'] = [
    'numpy', 'scipy', 'matplotlib', 'astropy', 'pandas',
    'tqdm', 'IPython', 'pytest',
    # The reduction engine; see requirements.txt on installing the local checkout
    'pypeit']
setup_keywords['zip_safe'] = False
setup_keywords['packages'] = find_packages()

if os.path.isdir('bin'):
    setup_keywords['scripts'] = [fname for fname in glob.glob(os.path.join('bin', '*'))
                                 if not os.path.basename(fname).endswith('.rst')]

setup(**setup_keywords)
