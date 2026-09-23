""" Vendor `pyodine` into this repository at a recorded commit.

Phase 3 forks `pyodine` rather than depending on it: from the moment the code
lands in `vendor/pyodine/` it is our source, changed by normal commits with
normal justifications (see `claude_prompts/data_phase3_prompt.md`, "We fork
pyodine").  This script performs that vendoring, and re-runs as a verifier.

Two things make a plain `cp -r` insufficient.

1. **Provenance.**  The upstream repository has exactly one commit, no tags and
   no releases, so "the version we forked" is a bare SHA that has to be pinned
   in code and checked, not remembered.  `UPSTREAM_COMMIT` below is that pin.

2. **Size.**  The upstream tree is 311 MB tracked, of which 303 MB is thirteen
   binary assets: two iodine atlases, two tutorial tarballs, the Arcturus
   reference spectrum, and the CARMENES/HITRAN telluric data.  Committing those
   into this repository would violate the project rule that bulk data lives in
   `../first-hires-exoplanet-data/`, and would bloat the history permanently.
   So files above `ASSET_THRESHOLD` are written to the external data directory
   and symlinked back into `vendor/pyodine/` at their upstream relative paths —
   `utilities_*/conf.py` resolves `../iodine_atlas` relative to itself, so the
   layout has to be preserved exactly.  The symlinks are gitignored; the
   manifest of what they point at, with sizes and SHA-256, is committed.

Nothing in the vendored tree is modified here.  Only git-tracked files are
taken, minus two exclusions: upstream commits 133 `__pycache__/*.pyc` files
(bytecode for Python 3.6/3.8/3.9, which this environment at 3.14 cannot use)
and five `.ipynb_checkpoints/*.ipynb` (Jupyter's autosaves of the tutorial
notebooks, 1.2 MB of duplicates).  Both are editor and build droppings that
`.gitignore` would refuse to commit anyway.  Dropping them keeps the tree
on disk identical to the tree in the repository, which matters because that
tree is about to be edited and the diff against upstream is the record of what
we changed.  They are build droppings, not code; no source file is touched.

Run from the repository root:

    conda run -n pypeit14 python -m first_hires_exoplanet.vendor_pyodine
    conda run -n pypeit14 python -m first_hires_exoplanet.vendor_pyodine --verify

"""

# Standard imports
import argparse
import csv
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Upstream repository.  The paper and docs point at gitlab.com/Heeren/pyodine;
#: the GitHub mirror is what is publicly reachable and is what phase 2 read.
UPSTREAM_URL = 'https://github.com/pepeheeren/pyodine.git'

#: The commit we fork.  `main`, 2022-11-04, "Adding pyodine to GitHub" -- the
#: repository's only commit.  Phase 2 read the code at this SHA via the GitHub
#: API; pinning it here means the fork's diff is meaningful.
UPSTREAM_COMMIT = '4488b0914fe5b272b787982647691045bff2604a'

#: Repository root (this file lives in <root>/first_hires_exoplanet/)
REPO_ROOT = Path(__file__).resolve().parent.parent

#: Where the fork lives
VENDOR_DIR = REPO_ROOT / 'vendor' / 'pyodine'

#: Where upstream's bulk binaries live, outside the repository
ASSET_ROOT = REPO_ROOT.parent / 'first-hires-exoplanet-data' / 'pyodine_assets'

#: Committed record of the externalised files
MANIFEST = REPO_ROOT / 'vendor' / 'pyodine_assets.csv'

#: Files at or above this size are externalised rather than committed (1 MiB)
ASSET_THRESHOLD = 1024 * 1024

#: Upstream-tracked droppings not carried into the fork (see the module
#: docstring).  A path containing any of these components is skipped.
DROPPINGS = ('__pycache__', '.ipynb_checkpoints')


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def sha256(path):
    """ SHA-256 of a file, read in chunks.

    Parameters
    ----------
    path : :obj:`pathlib.Path`
        File to digest.  Symlinks are followed.

    Returns
    -------
    :obj:`str`
        Hex digest.
    """
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def git(*args, cwd=None):
    """ Run a git command and return its stdout, raising on failure.

    Parameters
    ----------
    *args : :obj:`str`
        Arguments after ``git``.
    cwd : :obj:`str`, optional
        Working directory.

    Returns
    -------
    :obj:`str`
        Stripped stdout.
    """
    out = subprocess.run(('git',) + args, cwd=cwd, check=True,
                         capture_output=True, text=True)
    return out.stdout.strip()


def fetch_upstream(workdir):
    """ Clone upstream and check out the pinned commit.

    Parameters
    ----------
    workdir : :obj:`pathlib.Path`
        Directory to clone into (must not exist).

    Returns
    -------
    :obj:`pathlib.Path`
        The checkout.
    """
    print(f'Cloning {UPSTREAM_URL}')
    git('clone', '--quiet', UPSTREAM_URL, str(workdir))
    git('checkout', '--quiet', UPSTREAM_COMMIT, cwd=str(workdir))
    head = git('rev-parse', 'HEAD', cwd=str(workdir))
    if head != UPSTREAM_COMMIT:
        raise RuntimeError(f'checkout is {head}, expected {UPSTREAM_COMMIT}')
    print(f'  HEAD {head}')
    return workdir


def tracked_files(checkout):
    """ Git-tracked paths in the checkout, sorted, minus committed bytecode.

    Parameters
    ----------
    checkout : :obj:`pathlib.Path`
        An upstream checkout.

    Returns
    -------
    :obj:`list`
        Relative paths, as :obj:`str`.
    """
    return sorted(rel for rel in git('ls-files', cwd=str(checkout)).splitlines()
                  if not any(d in Path(rel).parts for d in DROPPINGS))


# ---------------------------------------------------------------------------
# Vendoring
# ---------------------------------------------------------------------------

def vendor(checkout):
    """ Copy the checkout into `VENDOR_DIR`, externalising the bulk binaries.

    Parameters
    ----------
    checkout : :obj:`pathlib.Path`
        Upstream checkout at `UPSTREAM_COMMIT`.

    Returns
    -------
    :obj:`list`
        One dict per externalised asset: path, bytes, sha256.
    """
    if VENDOR_DIR.exists():
        shutil.rmtree(VENDOR_DIR)
    VENDOR_DIR.mkdir(parents=True)
    ASSET_ROOT.mkdir(parents=True, exist_ok=True)

    assets, n_small, small_bytes = [], 0, 0
    for rel in tracked_files(checkout):
        src = checkout / rel
        dest = VENDOR_DIR / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        size = src.stat().st_size

        if size < ASSET_THRESHOLD:
            shutil.copy2(src, dest)
            n_small += 1
            small_bytes += size
            continue

        # Externalise: the bytes go to the data directory, and a relative
        # symlink puts them back where upstream's relative paths expect them.
        target = ASSET_ROOT / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, target)
        dest.symlink_to(os.path.relpath(target, dest.parent))
        assets.append({'path': rel, 'bytes': size, 'sha256': sha256(target)})
        print(f'  externalised {rel}  ({size / 1e6:.1f} MB)')

    print(f'{n_small} files committed ({small_bytes / 1e6:.1f} MB), '
          f'{len(assets)} externalised '
          f'({sum(a["bytes"] for a in assets) / 1e6:.1f} MB)')
    return assets


def write_manifest(assets):
    """ Write the committed record of the externalised assets.

    Parameters
    ----------
    assets : :obj:`list`
        As returned by :func:`vendor`.
    """
    with open(MANIFEST, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=('path', 'bytes', 'sha256'))
        writer.writeheader()
        writer.writerows(assets)
    print(f'Wrote {MANIFEST.relative_to(REPO_ROOT)}')


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

def verify(checkout):
    """ Compare the vendored tree against upstream, file by file.

    Reports *every* difference rather than the first, because after prompt 2
    the expected output is a short list of deliberate HIRES changes and the
    point of the check is to see that list whole.

    Parameters
    ----------
    checkout : :obj:`pathlib.Path`
        Upstream checkout at `UPSTREAM_COMMIT`.

    Returns
    -------
    :obj:`bool`
        True if the tree is byte-identical to upstream.
    """
    same, changed, missing, extra = 0, [], [], []

    expected = set()
    for rel in tracked_files(checkout):
        expected.add(rel)
        dest = VENDOR_DIR / rel
        if not dest.exists():
            missing.append(rel)
        elif sha256(dest) == sha256(checkout / rel):
            same += 1
        else:
            changed.append(rel)

    for path in VENDOR_DIR.rglob('*'):
        rel = path.relative_to(VENDOR_DIR)
        # Running the code writes 3.14 bytecode back into the tree; that is
        # not a change to the fork, and `.gitignore` drops it, so the verifier
        # ignores the same droppings the vendoring does.
        if any(d in rel.parts for d in DROPPINGS):
            continue
        if path.is_file() and str(rel) not in expected:
            extra.append(str(rel))

    print(f'{same} files identical to upstream')
    for label, paths in (('CHANGED', changed), ('MISSING', missing),
                         ('EXTRA', extra)):
        for rel in sorted(paths):
            print(f'  {label}: {rel}')
    return not (changed or missing or extra)


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
    parser.add_argument('--verify', action='store_true',
                        help='compare the vendored tree against upstream '
                             'instead of rewriting it')
    parser.add_argument('--checkout', type=str, default=None,
                        help='use this existing upstream checkout instead of '
                             'cloning (must be at the pinned commit)')
    return parser.parse_args() if options is None else parser.parse_args(options)


def main(options=None):
    """ Vendor or verify.

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
    print(f'pyodine {UPSTREAM_COMMIT}')

    tmpdir = None
    try:
        if args.checkout is not None:
            checkout = Path(args.checkout).resolve()
            head = git('rev-parse', 'HEAD', cwd=str(checkout))
            if head != UPSTREAM_COMMIT:
                print(f'ERROR: {checkout} is at {head}')
                return 1
        else:
            tmpdir = tempfile.mkdtemp(prefix='pyodine_upstream_')
            checkout = fetch_upstream(Path(tmpdir) / 'pyodine')

        if args.verify:
            return 0 if verify(checkout) else 1

        write_manifest(vendor(checkout))
        return 0
    finally:
        if tmpdir is not None:
            shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == '__main__':
    sys.exit(main())
