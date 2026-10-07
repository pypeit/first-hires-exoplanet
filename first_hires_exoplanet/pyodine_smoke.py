""" Run the vendored `pyodine` as shipped, on its own tutorial data.

Phase 3 forks `pyodine` (see `claude_prompts/data_phase3_prompt.md`).  Before
changing anything in the fork we need a baseline: does the code, unmodified,
run to completion in `pypeit14` on Python 3.14, four years after it was written
for Python 3.6-3.9?  Upstream ships no test suite -- there is no `tests/`, no
`conftest.py`, no CI -- so its own tutorial is the only executable check it
has, and this script is that tutorial turned into something re-runnable.

The tutorial reduces SONG spectra of sigma Draconis in three stages, which are
exactly the three stages phase 3 will need for HIRES:

* **template** -- deconvolve five I2-free stellar exposures against the LSF
  measured from five B-star exposures taken through the cell (prompt 3);
* **observations** -- forward-model I2 observations chunk by chunk against that
  template (prompts 4 and 5);
* **velocities** -- weight and combine the chunk velocities into a timeseries
  (prompt 7).

Each stage is timed and reported separately, and a failure in one does not stop
the others from being attempted, because "which stage breaks" is the useful
output.  The data are the tutorial tarballs that upstream tracks; they are
unpacked once into the external data directory (they are 100 MB compressed and
do not belong in the repository).

The observation stage is the slow one, so `--n-obs` limits how many of the 20
tutorial observations are modelled; two is enough to show the machinery works
and keeps a re-run to a few minutes.

As shipped, the code did not get past the template stage in this environment:
it called `np.float` and `np.NaN`, aliases NumPy removed in 2.0, against the
NumPy 2.5 in `pypeit14`.  Prompt 1 established that with a runtime shim and
prompt 2 fixed it in the fork, so this script now runs it straight.  Its job
from here is regression: the SONG tutorial is the only end-to-end check the
fork has, and every HIRES change has to leave it passing.

A word on how failure looks here.  `create_template` and
`model_single_observation` each wrap their whole body in `except Exception`,
write "Something went wrong!" and the traceback to their error log, and then
return normally.  Nothing propagates to the caller.  A stage that has failed
outright is therefore indistinguishable, from the outside, from a stage that
worked -- which is exactly the failure mode phase 3 was warned about.  So each
stage here reads its own error log and fails loudly if anything was written to
it, rather than trusting that no exception means no error.

Two environmental traps, neither of them a code bug, both of which cost time
here and will cost it again on HIRES:

* **Matplotlib's backend must not be interactive.**  The template stage writes
  diagnostic plots; on macOS that initialises `macosx`, which is CoreFoundation,
  and the observation stage then forks four `pathos` workers off a process that
  has CoreFoundation loaded.  The workers finish their fits, write their result
  files, and then hang forever in the pool teardown, so the run looks like a
  slow success rather than a deadlock.  Forcing `Agg` before anything imports
  `matplotlib` removes it.  It appears only when template and observations run
  in the same process, which is why it did not show up the first time the two
  stages were run separately.
* **`conda run` buffers.**  Without `--no-capture-output` nothing is printed
  until the process exits, which makes a hang indistinguishable from slow work.

Note that the velocity stage asks `barycorrpy` to resolve HIP96100 through
SIMBAD, so that stage needs a network connection.

Run from the repository root with:

    conda run -n pypeit14 python -m first_hires_exoplanet.pyodine_smoke
    conda run -n pypeit14 python -m first_hires_exoplanet.pyodine_smoke --n-obs 5

"""

# Standard imports
import argparse
import os
import re

# Before anything can import matplotlib: see the module docstring.  `pyodine`
# imports pyplot at module scope, so this has to happen at the top of the file,
# not inside a function.
os.environ.setdefault('MPLBACKEND', 'Agg')
import shutil
import sys
import tarfile
import time
import traceback
from pathlib import Path


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Repository root (this file lives in <root>/first_hires_exoplanet/)
REPO_ROOT = Path(__file__).resolve().parent.parent

#: The fork
VENDOR_DIR = REPO_ROOT / 'vendor' / 'pyodine'

#: Upstream's bulk binaries, externalised by `vendor_pyodine.py`
ASSET_ROOT = REPO_ROOT.parent / 'first-hires-exoplanet-data' / 'pyodine_assets'

#: Scratch space for the tutorial: unpacked inputs and all outputs
WORK_ROOT = REPO_ROOT.parent / 'first-hires-exoplanet-data' / 'pyodine_tutorial'

#: The tutorial tarballs, and the directory each unpacks into
TARBALLS = {
    'sigdra_template.tar.gz': 'sigdra_template',
    'sigdra_obs.tar.gz': 'sigdra_obs_tutorial',
}

#: SIMBAD identifier of the tutorial star, for the barycentric correction
TUTORIAL_STAR = 'HIP96100'

#: Width of the label column in the printed report
_LABEL = 14


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

def add_vendor_to_path():
    """ Put the fork on `sys.path`.

    `pyodine` is not installed into the environment on purpose: it is vendored
    source, and the tutorial itself recommends importing it by path.  The
    driver scripts (`pyodine_create_templates` and friends) live at the top of
    the tree rather than inside the package, so the tree root is what goes on
    the path.
    """
    if str(VENDOR_DIR) not in sys.path:
        sys.path.insert(0, str(VENDOR_DIR))


def unpack_tutorial_data():
    """ Unpack the tutorial tarballs into `WORK_ROOT`, once.

    Returns
    -------
    :obj:`dict`
        Directory name -> :obj:`pathlib.Path` of the unpacked tree.
    """
    WORK_ROOT.mkdir(parents=True, exist_ok=True)
    unpacked = {}
    for tarball, dirname in TARBALLS.items():
        dest = WORK_ROOT / dirname
        if not dest.exists():
            src = ASSET_ROOT / 'tutorial_data' / tarball
            print(f'  unpacking {tarball}')
            with tarfile.open(src) as tar:
                tar.extractall(WORK_ROOT, filter='data')
        unpacked[dirname] = dest
    return unpacked


def check_error_log(path):
    """ Raise if `pyodine` logged anything to an error log.

    See the module docstring: the pipeline swallows its own exceptions, so an
    error log is the only signal that a stage failed.

    Parameters
    ----------
    path : :obj:`pathlib.Path` or :obj:`str`
        The error log the stage was given.

    Raises
    ------
    RuntimeError
        If the log exists and is not empty.
    """
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return
    lines = path.read_text().splitlines()
    print(f'  {path} is not empty ({len(lines)} lines):')
    for line in lines[:3]:
        print(f'    {line}')
    # The exception line is the one worth putting in the summary table, and it
    # is not necessarily the last line: NumPy's removal messages append several
    # lines of migration advice after it.
    culprit = next((l for l in reversed(lines)
                    if re.match(r'^[\w.]*(Error|Exception|Warning)\b', l)),
                   lines[-1])
    raise RuntimeError(f'pyodine logged an error: {culprit[:160]}')


def fits_in(directory):
    """ Sorted list of FITS files in a directory.

    Parameters
    ----------
    directory : :obj:`pathlib.Path`
        Directory to list.

    Returns
    -------
    :obj:`list`
        Absolute paths, as :obj:`str`.
    """
    return sorted(str(p) for p in Path(directory).glob('*.fits'))


# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------

def stage_imports():
    """ Import every module of the fork.

    A 2022 code base on Python 3.14 can fail at import for reasons that have
    nothing to do with the science -- removed standard-library modules,
    tightened syntax -- and that failure is worth separating from a modelling
    failure.

    Returns
    -------
    :obj:`str`
        A one-line summary.
    """
    import importlib
    import pkgutil

    import pyodine

    modules = [m.name for m in pkgutil.walk_packages(pyodine.__path__,
                                                     'pyodine.')]
    modules += ['utilities_lick', 'utilities_song', 'utilities_waltz',
                'utilities_mtkent', 'pipe_lib', 'pyodine_create_templates',
                'pyodine_model_observations', 'pyodine_combine_vels',
                'pyodine_model_i2flat', 'pyodine_clean_i2']

    failed = []
    for name in modules:
        try:
            importlib.import_module(name)
        except BaseException as exc:
            failed.append(f'{name}: {type(exc).__name__}: {exc}')
    for line in failed:
        print(f'  FAILED {line}')
    if failed:
        raise ImportError(f'{len(failed)} of {len(modules)} modules')
    return f'{len(modules)} modules import'


def stage_template(data):
    """ Create the deconvolved sigma Dra template.

    Parameters
    ----------
    data : :obj:`dict`
        Unpacked tutorial directories, from :func:`unpack_tutorial_data`.

    Returns
    -------
    :obj:`str`
        A one-line summary.
    """
    import pyodine
    import pyodine_create_templates
    import utilities_song as utilities

    night = data['sigdra_template'] / '2018-05-16'
    out_dir = WORK_ROOT / 'temp_results'
    out_dir.mkdir(exist_ok=True)
    temp_outname = str(out_dir / 'temp_sigdra_2018-05-16.h5')

    ostar_files = fits_in(night / 'obs_ostar')
    temp_files = fits_in(night / 'obs_temp')
    print(f'  {len(ostar_files)} O-star, {len(temp_files)} stellar exposures')

    pyodine_create_templates.create_template(
        utilities, utilities.pyodine_parameters.Template_Parameters(),
        ostar_files, temp_files, temp_outname,
        plot_dir=str(out_dir / 'plots'),
        obs_sum_outname=str(out_dir / 'sigdra_2018-05-16_summed.fits'),
        error_log=str(out_dir / 'error.log'),
        info_log=str(out_dir / 'info.log'), quiet=True)

    check_error_log(out_dir / 'error.log')
    template = pyodine.template.base.StellarTemplate_Chunked(temp_outname)
    return f'{len(template)} template chunks written'


def stage_observations(data, n_obs):
    """ Forward-model the first `n_obs` tutorial observations.

    Parameters
    ----------
    data : :obj:`dict`
        Unpacked tutorial directories.
    n_obs : :obj:`int`
        How many observations to model.

    Returns
    -------
    :obj:`str`
        A one-line summary.
    """
    import numpy as np
    import pyodine
    import pyodine_model_observations
    import utilities_song as utilities

    temp_file = str(WORK_ROOT / 'temp_results' / 'temp_sigdra_2018-05-16.h5')
    obs_files = fits_in(data['sigdra_obs_tutorial'])[:n_obs]
    out_root = WORK_ROOT / 'obs_results'
    out_root.mkdir(exist_ok=True)

    plot_dirs, res_files, error_files, info_files = [], [], [], []
    for obs_file in obs_files:
        base = Path(obs_file).stem
        out_dir = out_root / base
        out_dir.mkdir(exist_ok=True)
        plot_dirs.append(str(out_dir))
        res_files.append([str(out_dir / f'{base}_res1.h5')])
        error_files.append(str(out_dir / 'error.log'))
        info_files.append(str(out_dir / 'info.log'))

    print(f'  modelling {len(obs_files)} observations')
    pyodine_model_observations.model_multi_observations(
        utilities, utilities.pyodine_parameters.Parameters(),
        obs_files, [temp_file] * len(obs_files), plot_dirs=plot_dirs,
        res_files=res_files, error_files=error_files, info_files=info_files,
        quiet=True)

    for error_file in error_files:
        check_error_log(error_file)

    # A chunk velocity scatter is the first number that says the model
    # actually fitted something, rather than merely running to the end.
    # `load_results` on an .h5 returns a dict of arrays, one entry per fitted
    # parameter, each of length n_chunks -- not the list of result objects the
    # tutorial's dill path gives.
    results = pyodine.fitters.results_io.load_results(res_files[0][0])
    velocity = np.asarray(results['params']['velocity'])
    print(f'  chunk velocities: median {np.nanmedian(velocity):.1f} m/s, '
          f'scatter {np.nanstd(velocity):.1f} m/s, '
          f'{np.sum(~np.isfinite(velocity))} not finite')
    print(f'  median reduced chi2: {np.nanmedian(results["redchi2"]):.1f}')
    return f'{len(obs_files)} observations, {len(velocity)} chunks each'


def stage_velocities(n_obs):
    """ Weight and combine the chunk velocities into a timeseries.

    Parameters
    ----------
    n_obs : :obj:`int`
        How many observations were modelled.

    Returns
    -------
    :obj:`str`
        A one-line summary.
    """
    import numpy as np
    import pyodine_combine_vels
    import utilities_song as utilities

    out_root = WORK_ROOT / 'obs_results'
    res_files = sorted(str(p) for p in out_root.glob('*/*_res1.h5'))[:n_obs]
    plot_dir = WORK_ROOT / 'vel_results'
    plot_dir.mkdir(exist_ok=True)
    vels_out = str(plot_dir / 'sigdra_tutorial.vels')

    pyodine_combine_vels.combine_velocity_results(
        utilities.timeseries_parameters.Timeseries_Parameters(),
        res_files=res_files, plot_dir=str(plot_dir),
        comb_res_out=str(plot_dir / 'sigdra_tutorial_comb.h5'),
        vels_out=vels_out, bary_dict={'star_name': TUTORIAL_STAR},
        error_log=str(plot_dir / 'error.log'),
        info_log=str(plot_dir / 'info.log'), quiet=True)

    check_error_log(plot_dir / 'error.log')
    rows = [l for l in open(vels_out) if l.strip() and not l.startswith('#')]
    print(f'  {vels_out}')
    for row in rows:
        print(f'    {row.rstrip()}')
    return f'{len(rows)} epochs in the velocity file'


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
    parser.add_argument('--n-obs', type=int, default=2,
                        help='observations to model (of 20); default 2')
    parser.add_argument('--stages', type=str, default='all',
                        help='comma-separated subset of imports,template,'
                             'observations,velocities; default all')
    parser.add_argument('--clean', action='store_true',
                        help='delete previous tutorial outputs first')
    return parser.parse_args() if options is None else parser.parse_args(options)


def main(options=None):
    """ Run the requested stages and report.

    Parameters
    ----------
    options : :obj:`list`, optional
        Arguments to parse instead of `sys.argv`.

    Returns
    -------
    :obj:`int`
        Process exit status: 0 if every attempted stage passed.
    """
    args = parse_args(options)
    add_vendor_to_path()

    if args.clean:
        for sub in ('temp_results', 'obs_results', 'vel_results'):
            shutil.rmtree(WORK_ROOT / sub, ignore_errors=True)

    import matplotlib

    print(f'pyodine tree : {VENDOR_DIR}')
    print(f'mpl backend  : {matplotlib.get_backend()}')
    print(f'work dir     : {WORK_ROOT}')
    data = unpack_tutorial_data()

    wanted = ('imports', 'template', 'observations', 'velocities')
    if args.stages != 'all':
        wanted = tuple(s.strip() for s in args.stages.split(','))
    runners = {
        'imports': stage_imports,
        'template': lambda: stage_template(data),
        'observations': lambda: stage_observations(data, args.n_obs),
        'velocities': lambda: stage_velocities(args.n_obs),
    }

    report = []
    for name in wanted:
        print(f'\n--- {name} ---')
        start = time.time()
        try:
            summary = runners[name]()
            report.append((name, 'PASS', time.time() - start, summary))
        except BaseException as exc:
            traceback.print_exc()
            report.append((name, 'FAIL', time.time() - start,
                           f'{type(exc).__name__}: {exc}'))

    print('\n' + '=' * 72)
    for name, status, elapsed, summary in report:
        print(f'{name:{_LABEL}s} {status}  {elapsed:7.1f} s  {summary}')
    print('=' * 72)
    return 0 if all(r[1] == 'PASS' for r in report) else 1


if __name__ == '__main__':
    sys.exit(main())
