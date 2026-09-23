""" Summarise the fitted chunk parameters in a `pyodine` result file.

Phase 3 prompt 2 needed this to answer a question its own change raised: after
the iodine depth model was made Beer-Lambert, the per-chunk velocity standard
deviation on the SONG tutorial moved from 911 to 1483 m/s, which looks like a
serious regression and is not one.  The standard deviation of a chunk velocity
distribution is the wrong statistic: a handful of chunks fail outright and land
thousands of m/s away, and they dominate it.  The robust width barely moves.

So this prints both, plus the tails, for every fitted parameter.  It will be
wanted again in prompts 5 and 7, which ask for per-chunk scatter by name.

Run from the repository root with:

    conda run -n pypeit14 python -m first_hires_exoplanet.pyodine_chunk_stats \\
        <result.h5> [<result.h5> ...]

"""

# Standard imports
import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault('MPLBACKEND', 'Agg')

import numpy as np

#: The vendored fork
VENDOR_DIR = Path(__file__).resolve().parent.parent / 'vendor' / 'pyodine'

if str(VENDOR_DIR) not in sys.path:
    sys.path.insert(0, str(VENDOR_DIR))


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def robust_sigma(values):
    """ Scatter from the median absolute deviation, scaled to a Gaussian sigma.

    Parameters
    ----------
    values : :obj:`numpy.ndarray`
        The sample.

    Returns
    -------
    :obj:`float`
        1.4826 * MAD, or NaN for an empty sample.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    return 1.4826 * np.median(np.abs(values - np.median(values)))


def summarise(values, n_sigma=3.):
    """ Robust and non-robust descriptions of one fitted parameter.

    Parameters
    ----------
    values : :obj:`numpy.ndarray`
        One parameter across all chunks.
    n_sigma : :obj:`float`
        How far from the median a chunk must sit to count as an outlier.

    Returns
    -------
    :obj:`dict`
        median, std, robust sigma, extremes, and the outlier count.
    """
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    sigma = robust_sigma(finite)
    deviant = (np.abs(finite - np.median(finite)) > n_sigma * sigma
               if np.isfinite(sigma) and sigma > 0 else np.zeros(0, dtype=bool))
    return {'n': values.size,
            'n_bad': int(values.size - finite.size),
            'median': float(np.median(finite)) if finite.size else np.nan,
            'std': float(np.std(finite)) if finite.size else np.nan,
            'robust': float(sigma),
            'min': float(finite.min()) if finite.size else np.nan,
            'max': float(finite.max()) if finite.size else np.nan,
            'n_outliers': int(np.sum(deviant))}


def report(path, n_sigma=3.):
    """ Print a table of every fitted parameter in one result file.

    Parameters
    ----------
    path : :obj:`str`
        A `pyodine` `.h5` result file.
    n_sigma : :obj:`float`
        Outlier threshold, in robust sigmas.

    Returns
    -------
    :obj:`dict`
        Parameter name -> the dict from :func:`summarise`.
    """
    from pyodine.fitters.results_io import load_results

    results = load_results(path)
    print(f'\n{path}')
    print(f'{"parameter":18s} {"median":>12s} {"robust sig":>11s} '
          f'{"std":>12s} {"min":>12s} {"max":>12s} {"outliers":>9s}')

    stats = {}
    for name, values in sorted(results['params'].items()):
        row = summarise(values, n_sigma=n_sigma)
        stats[name] = row
        print(f'{name:18s} {row["median"]:12.4f} {row["robust"]:11.4f} '
              f'{row["std"]:12.4f} {row["min"]:12.4f} {row["max"]:12.4f} '
              f'{row["n_outliers"]:6d}/{row["n"]:d}')

    chi2 = np.asarray(results['redchi2'], dtype=float)
    print(f'{"redchi2":18s} {np.nanmedian(chi2):12.4g}')
    return stats


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_args(options=None):
    """ Command line.

    Parameters
    ----------
    options : :obj:`list`, optional
        Arguments to parse instead of `sys.argv`.

    Returns
    -------
    :obj:`argparse.Namespace`
    """
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[1].strip())
    parser.add_argument('results', nargs='+', help='pyodine .h5 result files')
    parser.add_argument('--n-sigma', type=float, default=3.,
                        help='outlier threshold in robust sigmas; default 3')
    return parser.parse_args() if options is None else parser.parse_args(options)


def main(options=None):
    """ Report on each result file given.

    Parameters
    ----------
    options : :obj:`list`, optional
        Arguments to parse instead of `sys.argv`.

    Returns
    -------
    :obj:`int`
        Process exit status.
    """
    args = parse_args(options)
    for path in args.results:
        report(path, n_sigma=args.n_sigma)
    return 0


if __name__ == '__main__':
    sys.exit(main())
