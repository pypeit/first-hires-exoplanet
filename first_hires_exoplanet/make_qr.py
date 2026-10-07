""" QR code for the public report, for slides (phase 3, prompt 17).

Writes, to ``docs/figs/``:

    qr_public_report.png   1200 px square, black on white, for slides
    qr_public_report.svg   vector, scales to any size

Both encode ``URL`` at error-correction level H (30% of the code can be
lost or covered and it still scans), with the standard 4-module quiet zone.
Keep the white border when placing it on a slide: a dark background up to
the edge stops many phone cameras reading it.

Needs ``segno`` (``pip install segno``; pure Python, no dependencies).

Run with:

    conda run -n pypeit14 python -m first_hires_exoplanet.make_qr

"""

# Standard imports
import os
import argparse

import segno


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(_HERE)
FIGS = os.path.join(ROOT, 'docs', 'figs')

URL = 'https://github.com/pypeit/first-hires-exoplanet/blob/main/docs/public_HD187123b.md'

#: Target width of the PNG, in pixels
PNG_PX = 1200


def make(url=URL, outroot=os.path.join(FIGS, 'qr_public_report')):
    qr = segno.make(url, error='h', micro=False, boost_error=False)
    width = qr.symbol_size(scale=1, border=4)[0]
    scale = max(1, PNG_PX // width)
    qr.save(outroot + '.png', scale=scale, border=4, dark='black', light='white')
    qr.save(outroot + '.svg', scale=10, border=4, dark='black', light='white')
    print('  {:s}'.format(url))
    print('  version {:d}, error level {:s}, {:d} modules with border, PNG {:d} px'.format(
        qr.version, qr.error, width, width * scale))
    print('  wrote {:s}.png and .svg'.format(outroot))
    return outroot + '.png'


def main():
    p = argparse.ArgumentParser(description='QR code for the public report')
    p.add_argument('--url', type=str, default=URL)
    args = p.parse_args()
    make(args.url)


if __name__ == '__main__':
    main()
