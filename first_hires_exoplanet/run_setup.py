""" Run ``pypeit_setup`` over the 1998 HIRES frames of HD 187123 (prompt 7).

Reproduces the automatic ``pypeit_setup`` pass reported in
claude_prompts/data_phase1_prompt.md, prompt 7, so that anyone can regenerate
the same ``.pypeit``, ``.sorted``, ``.calib`` and ``.obslog`` files from the
raw frames fetched by ``koa_download.py``.  Two runs are made, each into its
own directory under the output root:

  setup_jul16         the 1998-07-16 night alone (21 frames) -- the literal
                      "run pypeit_setup on that night";
  setup_jul16_jul14   1998-07-16 plus the 1998-07-14 clean B1 flats, arcs and
                      biases (30 frames) -- the input prompt 8 reduces, since
                      07-16 has no hatch-closed flat in the science decker.

``pypeit_setup`` writes different products depending on ``-c``: with
``-c all`` it writes one ``<spectrograph>_<setup>.pypeit`` (plus ``.calib``)
per configuration and *no* ``.sorted`` file; without ``-c`` it writes
``setup_files/<spectrograph>.sorted``, ``.calib`` and ``.obslog`` and no
``.pypeit`` file (pypeit/scripts/setup.py, ``if args.cfg_split is None``).
So each run invokes ``pypeit_setup`` twice, once each way, and keeps the
full stdout/stderr of both passes next to the products:

  <run>/keck_hires_orig_A/keck_hires_orig_A.pypeit   (and _B, ...)
  <run>/keck_hires_orig_A/keck_hires_orig_A.calib
  <run>/setup_files/keck_hires_orig.{sorted,calib,obslog}
  <run>/pypeit_setup_stdout.txt          (the ``-c all`` pass)
  <run>/pypeit_setup_sorted_stdout.txt   (the setup-only pass)
  <run>/pypeit_setup_*.log               (pypeit's own log files)

Nothing here edits the generated files: the point is to record what the
*automatic* pass produces.  Hand edits belong to prompt 8, in copies.

Generated files are large-ish and tied to absolute raw-data paths, so they
live outside the git repository.  The default output root is a sibling of
the repository, and the script refuses to write anywhere under the repo.

The ``pypeit_setup`` executable is taken from the same environment as the
running interpreter, so run this as

    conda run -n pypeit14 python -m first_hires_exoplanet.run_setup

Options: ``--run setup_jul16`` (repeatable) restricts to one run;
``--rawroot DIR`` / ``--outroot DIR`` move the raw-data and output roots;
``--overwrite`` passes ``-o`` to ``pypeit_setup`` and re-runs into existing
directories; ``--dry-run`` prints the commands without running them.
"""

# Standard imports
import argparse
import os
import shutil
import subprocess
import sys


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: PypeIt spectrograph name for the pre-2004 single-CCD HIRES
SPECTROGRAPH = 'keck_hires_orig'

#: The repository root (parent of this package).  Nothing is written below it.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: Sibling of the repository holding raw and reduced data, outside git.
DATA_ROOT = os.path.join(os.path.dirname(REPO_ROOT), 'first-hires-exoplanet-data')

#: Default raw-data root: where ``koa_download.py`` puts the frames.
DEFAULT_RAWROOT = os.path.join(DATA_ROOT, 'raw')

#: Default output root for the setup products.
DEFAULT_OUTROOT = os.path.join(DATA_ROOT, 'redux')

#: Runs: output sub-directory -> list of raw sub-directories (under the raw
#: root) to hand to ``pypeit_setup -r``.  ``inspect_only/`` (T. Bida's
#: two-amplifier frame) is deliberately absent from both.
RUNS = {
    'setup_jul16': ['1998jul16'],
    'setup_jul16_jul14': ['1998jul16', '1998jul14'],
}

#: stdout/stderr capture files for the two passes
STDOUT_CFG = 'pypeit_setup_stdout.txt'
STDOUT_SORTED = 'pypeit_setup_sorted_stdout.txt'


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def check_outroot(outroot):
    """ Refuse an output root inside the git repository. """
    real_out = os.path.realpath(outroot)
    real_repo = os.path.realpath(REPO_ROOT)
    if real_out == real_repo or real_out.startswith(real_repo + os.sep):
        raise SystemExit('Refusing to write setup products into the repository: {:s}'.format(outroot))


def find_pypeit_setup():
    """ Locate ``pypeit_setup`` in the interpreter's own environment. """
    candidate = os.path.join(os.path.dirname(sys.executable), 'pypeit_setup')
    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
        return candidate
    found = shutil.which('pypeit_setup')
    if found is None:
        raise SystemExit('pypeit_setup not found next to {:s} or on PATH; '
                         'run under `conda run -n pypeit14`'.format(sys.executable))
    return found


def setup_commands(exe, raw_dirs, outdir, overwrite):
    """ Return the two ``pypeit_setup`` command lines for one run.

    Args:
        exe (str): path to ``pypeit_setup``
        raw_dirs (list of str): raw directories for ``-r``
        outdir (str): output directory for ``-d``
        overwrite (bool): pass ``-o``

    Returns:
        list of (label, argv, capture_file)
    """
    base = [exe, '-s', SPECTROGRAPH, '-r'] + raw_dirs + ['-d', outdir]
    if overwrite:
        base.append('-o')
    return [('cfg_split', base + ['-c', 'all'], STDOUT_CFG),
            ('setup_only', base, STDOUT_SORTED)]


def run_one(exe, name, raw_dirs, outroot, overwrite, dry_run):
    """ Run both passes for one entry of RUNS; return True on success. """
    outdir = os.path.join(outroot, name)
    for raw in raw_dirs:
        if not os.path.isdir(raw):
            print('  MISSING raw directory {:s}; skipping run {:s}'.format(raw, name))
            return False
    if os.path.isdir(outdir) and os.listdir(outdir) and not overwrite and not dry_run:
        print('  {:s} exists and is not empty; use --overwrite to re-run into it'.format(outdir))
        return False

    ok = True
    for label, argv, capture in setup_commands(exe, raw_dirs, outdir, overwrite):
        print('  [{:s}] {:s}'.format(label, ' '.join(argv)))
        if dry_run:
            continue
        os.makedirs(outdir, exist_ok=True)
        # pypeit_setup writes its own log file into the current directory
        with open(os.path.join(outdir, capture), 'w', encoding='utf-8') as fobj:
            proc = subprocess.run(argv, cwd=outdir, stdout=fobj, stderr=subprocess.STDOUT)
        print('           exit {:d}; stdout+stderr -> {:s}'.format(
            proc.returncode, os.path.join(outdir, capture)))
        ok &= proc.returncode == 0
    return ok


def summarize(outdir):
    """ Print the frame-type table from the ``.sorted`` file of one run. """
    sorted_file = os.path.join(outdir, 'setup_files', SPECTROGRAPH + '.sorted')
    if not os.path.isfile(sorted_file):
        print('  no .sorted file at {:s}'.format(sorted_file))
        return
    print('\n  ' + sorted_file)
    n_typed = n_untyped = 0
    with open(sorted_file, encoding='utf-8') as fobj:
        for line in fobj:
            line = line.rstrip('\n')
            if line.startswith('Setup ') or (':' in line and '|' not in line
                                             and not line.startswith('#')):
                # the setup heading and its indented configuration keys
                print('  ' + line)
            elif line.startswith('# HI.'):
                n_untyped += 1
                print('    UNTYPED  ' + line[2:].split('|')[0].strip())
            elif line.startswith('  HI.'):
                n_typed += 1
                cols = [c.strip() for c in line.split('|')]
                print('    {:22s} {:s}'.format(cols[0], cols[1]))
    print('  {:d} frames typed, {:d} commented out'.format(n_typed, n_untyped))
    products = []
    for root, _, files in os.walk(outdir):
        products += [os.path.join(root, f) for f in files
                     if f.endswith(('.pypeit', '.sorted', '.calib', '.obslog'))]
    for path in sorted(products):
        print('    ' + path)


# ---------------------------------------------------------------------------

def main():
    """ Run pypeit_setup for each entry of RUNS and summarize the products. """
    parser = argparse.ArgumentParser(
        description='Run pypeit_setup ({:s}) over the 1998 HIRES frames of HD 187123, '
                    'into a directory outside the repository.'.format(SPECTROGRAPH))
    parser.add_argument('--run', action='append', choices=sorted(RUNS),
                        help='which run to make (repeatable; default: all)')
    parser.add_argument('--rawroot', default=DEFAULT_RAWROOT,
                        help='raw-data root holding 1998jul16/ and 1998jul14/ '
                             '(default: {:s})'.format(DEFAULT_RAWROOT))
    parser.add_argument('--outroot', default=DEFAULT_OUTROOT,
                        help='output root; products go in <outroot>/<run>/ '
                             '(default: {:s})'.format(DEFAULT_OUTROOT))
    parser.add_argument('--overwrite', action='store_true',
                        help='pass -o to pypeit_setup and re-run into existing directories')
    parser.add_argument('--dry-run', action='store_true',
                        help='print the pypeit_setup commands; run nothing')
    args = parser.parse_args()

    check_outroot(args.outroot)
    exe = find_pypeit_setup()
    runs = args.run or sorted(RUNS)

    n_fail = 0
    for name in runs:
        raw_dirs = [os.path.join(args.rawroot, d) for d in RUNS[name]]
        print('\nRun {:s}: {:d} raw director{:s}'.format(name, len(raw_dirs),
                                                        'y' if len(raw_dirs) == 1 else 'ies'))
        if not run_one(exe, name, raw_dirs, args.outroot, args.overwrite, args.dry_run):
            n_fail += 1
            continue
        if not args.dry_run:
            summarize(os.path.join(args.outroot, name))

    if n_fail:
        print('\n*** {:d} run(s) did not complete ***'.format(n_fail))
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
