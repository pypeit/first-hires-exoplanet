""" Environment check for the HD 187123 reduction recipe.

Confirms that the PypeIt installation in front of us is the one we intend to
quote: an editable checkout of the PypeIt source, at a known git commit, whose
reported version string actually embeds that commit, and which knows about the
original (pre-2004, single Tektronix CCD) HIRES detector as `keck_hires_orig`.

The checks are deliberately simple and print everything they look at, so the
output can be pasted into a log.  Any failure is collected and reported at the
end, and the script then exits with a non-zero status.

Background: `pypeit.__version__` is written by setuptools_scm into the
gitignored `pypeit/pkg/version.py` *at install time*.  Switching branches in an
editable install changes the code that runs but not that file, so the version
can silently lag the checkout by hundreds of commits.  We hit exactly this
(see claude_prompts/data_phase1_prompt.md, 2026-09-19 log).  The fix is
`pip install -e . --no-deps` in the PypeIt checkout; this script tells you
whether it is needed.

Run from the repository root with:

    conda run -n pypeit14 python -m first_hires_exoplanet.check_env

"""

# Standard imports
import os
import subprocess
import sys

import pypeit
from pypeit.spectrographs.util import available_spectrographs, load_spectrograph


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: The PypeIt spectrograph for original-detector HIRES data (before ~Aug 2004)
SPECTROGRAPH = 'keck_hires_orig'

#: Width of the label column in the printed report
_LABEL = 22


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _show(label, value):
    """ Print one aligned `label : value` line. """
    print('{:<{w}s}: {}'.format(label, value, w=_LABEL))


def _git(checkout, *args):
    """ Run a git command in `checkout` and return its stripped stdout.

    Args:
        checkout (str): directory holding the git working tree.
        *args (str): arguments to pass after `git`.

    Returns:
        str or None: command output, or None if git failed or is missing.
    """
    try:
        result = subprocess.run(['git', '-C', checkout] + list(args),
                                capture_output=True, text=True, check=True)
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def _embedded_sha(version):
    """ Extract the commit hash that setuptools_scm embedded in a version.

    setuptools_scm writes the local segment as `+g<sha>` optionally followed by
    `.d<date>` or `.dirty`, e.g. `2.0.2.dev1216+gf3a1f1d27`.

    Args:
        version (str): the `pypeit.__version__` string.

    Returns:
        tuple: (sha, dirty) -- the embedded hex string (or None if there is
        none) and whether the version flags an unclean tree.
    """
    if '+' not in version:
        return None, False
    local = version.split('+', 1)[1]
    parts = local.split('.')
    sha = parts[0][1:] if parts[0].startswith('g') else None
    dirty = any(p == 'dirty' or p.startswith('d') for p in parts[1:])
    return sha, dirty


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

def check_install():
    """ Report the PypeIt version and where it is imported from.

    Returns:
        str: the directory containing the `pypeit` package, i.e. the checkout
        root for an editable install.
    """
    print('\n--- PypeIt install ---')
    _show('pypeit.__version__', pypeit.__version__)
    _show('pypeit.__file__', pypeit.__file__)
    package_dir = os.path.dirname(os.path.abspath(pypeit.__file__))
    checkout = os.path.dirname(package_dir)
    _show('checkout (inferred)', checkout)
    return checkout


def check_git(checkout, failures):
    """ Report the git state of the PypeIt checkout.

    Returns:
        str or None: the full HEAD sha, or None if this is not a git checkout.
    """
    print('\n--- PypeIt git checkout ---')
    if _git(checkout, 'rev-parse', '--is-inside-work-tree') != 'true':
        print('Not a git checkout (or git is not installed).  Cannot record a '
              'commit hash; the reduction is not reproducible from this '
              'install alone.')
        failures.append('PypeIt is not installed from a git checkout')
        return None

    sha = _git(checkout, 'rev-parse', 'HEAD')
    short = _git(checkout, 'rev-parse', '--short', 'HEAD')
    branch = _git(checkout, 'rev-parse', '--abbrev-ref', 'HEAD')
    subject = _git(checkout, 'log', '-1', '--format=%s (%ci)')
    status = _git(checkout, 'status', '--porcelain')

    _show('HEAD', sha)
    _show('HEAD (short)', short)
    _show('branch', branch)
    _show('commit', subject)
    if status:
        _show('working tree', 'DIRTY')
        for line in status.splitlines():
            print(' ' * (_LABEL + 2) + line)
        failures.append('PypeIt working tree is not clean')
    else:
        _show('working tree', 'clean')
    return sha


def check_version_matches(sha, failures):
    """ Warn loudly if `pypeit.__version__` does not embed the current HEAD. """
    print('\n--- Version string vs. HEAD ---')
    embedded, dirty = _embedded_sha(pypeit.__version__)
    _show('embedded commit', embedded)
    if sha is None:
        print('No git HEAD to compare against; skipping.')
        return
    if embedded is None or not sha.startswith(embedded) or len(embedded) < 7:
        print('*** WARNING: pypeit.__version__ does NOT embed the current HEAD.')
        print('*** The version is stamped into every PypeIt output header, so')
        print('*** anything reduced now would be mislabelled.  Refresh it with:')
        print('***     cd <checkout> && conda run -n pypeit14 pip install -e . '
              '--no-deps')
        failures.append('pypeit.__version__ is stale relative to git HEAD')
    else:
        _show('matches HEAD', 'yes')
    if dirty:
        print('*** WARNING: version string carries a dirty/date suffix, i.e. the')
        print('*** tree was modified when the version was generated.')
        failures.append('pypeit.__version__ was generated from a dirty tree')


def check_spectrograph(failures):
    """ Confirm `keck_hires_orig` is registered and fully constructs. """
    print('\n--- Spectrograph: {:s} ---'.format(SPECTROGRAPH))
    _show('n_spectrographs', len(available_spectrographs))
    listed = SPECTROGRAPH in available_spectrographs
    _show('listed', listed)
    if not listed:
        failures.append('{:s} not in available_spectrographs'.format(SPECTROGRAPH))
        return

    try:
        spec = load_spectrograph(SPECTROGRAPH)
    except Exception as err:      # report, do not crash
        print('load_spectrograph failed: {!r}'.format(err))
        failures.append('load_spectrograph({!r}) raised'.format(SPECTROGRAPH))
        return
    _show('name', spec.name)
    _show('ndet', spec.ndet)
    _show('supported', spec.supported)
    _show('configuration_keys', spec.configuration_keys())

    try:
        det = spec.get_detector_par(1)
        _show('platescale', det['platescale'])
        _show('ronoise', det['ronoise'])
        _show('saturation', det['saturation'])
    except Exception as err:
        print('get_detector_par(1) failed: {!r}'.format(err))
        failures.append('get_detector_par(1) raised')

    try:
        par = spec.default_pypeit_par()
        _show('rdx.detnum', par['rdx']['detnum'])
    except Exception as err:
        print('default_pypeit_par() failed: {!r}'.format(err))
        failures.append('default_pypeit_par() raised')


# ---------------------------------------------------------------------------

def main():
    """ Run every check; exit non-zero if any failed. """
    failures = []

    checkout = check_install()
    sha = check_git(checkout, failures)
    check_version_matches(sha, failures)
    check_spectrograph(failures)

    print('\n--- Summary ---')
    if failures:
        for item in failures:
            print('FAIL: {:s}'.format(item))
        sys.exit(1)
    print('All checks passed.')


if __name__ == '__main__':
    main()
