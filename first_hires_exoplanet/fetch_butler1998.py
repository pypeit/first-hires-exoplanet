""" Fetch the 1998 discovery paper and extract its velocity figures (phase 3, prompt 14).

Butler, Marcy, Vogt & Apps 1998, PASP, 110, 1389, "A Planet with a 3.1 Day
Period around a Solar Twin" (DOI 10.1086/316287).  Its page 2 carries the two
velocity figures, each embedded as a 300-dpi 1-bit image with its axes:

    Fig. 1  the 1998 July 15-19 run against time (JD - 2451000)
    Fig. 2  all 20 Keck velocities against orbital phase

Both, and the PDF, are written to ``../first-hires-exoplanet-data/literature/``,
outside the repository: the figures are the ASP's copyright and are not
committed.  ``make_slides.py`` reads them from there.

IOPscience has blocked automated downloads before (phase 0, context_prompts.md).
If the download fails, save the PDF by hand as
``../first-hires-exoplanet-data/literature/butler1998_pasp110_1389.pdf`` and
re-run; the extraction step then uses it.

Needs poppler's ``pdfimages`` on the PATH (``brew install poppler``).

Run with:

    conda run -n pypeit14 python -m first_hires_exoplanet.fetch_butler1998

"""

# Standard imports
import os
import glob
import shutil
import argparse
import subprocess
import urllib.request

from PIL import Image


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(os.path.dirname(os.path.dirname(_HERE)),
                         'first-hires-exoplanet-data')
LIT = os.path.join(DATA_ROOT, 'literature')
PDF = os.path.join(LIT, 'butler1998_pasp110_1389.pdf')

URL = 'https://iopscience.iop.org/article/10.1086/316287/pdf'

#: Page 2 of the PDF holds Figs. 1 and 2, in that order
FIG_PAGE = 2
#: Output name for each image on that page, by order of appearance
FIG_NAMES = ['butler1998_fig1_july.png', 'butler1998_fig2_phased.png']


def fig_paths():
    """ Where the extracted figures live (used by make_slides.py). """
    return [os.path.join(LIT, f) for f in FIG_NAMES]


def download(url=URL, out=PDF, force=False):
    """ Download the PDF unless it is already on disk. """
    if os.path.isfile(out) and not force:
        print('  have {:s}'.format(out))
        return out
    os.makedirs(os.path.dirname(out), exist_ok=True)
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=60) as resp:
        body = resp.read()
    if not body.startswith(b'%PDF'):
        raise IOError('{:s} did not return a PDF (bot protection?). Save the PDF by '
                      'hand as {:s} and re-run.'.format(url, out))
    with open(out, 'wb') as f:
        f.write(body)
    print('  wrote {:s} ({:d} bytes)'.format(out, len(body)))
    return out


def extract(pdf=PDF):
    """ Pull the page-2 images out of the PDF as PNGs, as black on white. """
    if shutil.which('pdfimages') is None:
        raise RuntimeError('pdfimages (poppler) is not on the PATH')
    tmp = os.path.join(LIT, '_tmp_butler1998')
    os.makedirs(tmp, exist_ok=True)
    subprocess.run(['pdfimages', '-f', str(FIG_PAGE), '-l', str(FIG_PAGE), '-png',
                    pdf, os.path.join(tmp, 'im')], check=True)
    ims = sorted(glob.glob(os.path.join(tmp, 'im-*.png')))
    if len(ims) != len(FIG_NAMES):
        raise ValueError('Expected {:d} images on page {:d}, found {:d}'.format(
            len(FIG_NAMES), FIG_PAGE, len(ims)))
    outs = fig_paths()
    for src, dst in zip(ims, outs):
        # 1-bit images: store as 8-bit greyscale so every viewer renders them alike
        with Image.open(src) as im:
            im.convert('L').save(dst)
        print('  wrote {:s}'.format(dst))
    shutil.rmtree(tmp)
    return outs


def main():
    p = argparse.ArgumentParser(description='Fetch Butler et al. 1998 and extract its velocity figures')
    p.add_argument('--force', action='store_true', help='Re-download the PDF')
    args = p.parse_args()
    extract(download(force=args.force))


if __name__ == '__main__':
    main()
