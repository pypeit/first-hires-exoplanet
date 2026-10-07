""" Shared setup for the tests of the vendored `pyodine` fork.

The fork is source, not an installed package: it lives at `vendor/pyodine/`
and is imported by path, exactly as `pyodine`'s own tutorial recommends.  This
file puts the tree on `sys.path` so that `import pyodine` inside a test means
*our* fork and nothing else, and pins the Matplotlib backend before anything
can import `pyplot`.

That second point is not cosmetic.  `pyodine` imports `pyplot` at module
scope; on macOS that initialises the `macosx` backend, and any later `fork`
(the fitting runs four `pathos` workers) then hangs in CoreFoundation after
writing correct output.  Phase 3 prompt 1 lost the better part of an hour to
it twice.
"""

# Standard imports
import os
import sys
from pathlib import Path

os.environ.setdefault('MPLBACKEND', 'Agg')

#: The vendored fork, at the repository root
VENDOR_DIR = Path(__file__).resolve().parents[2] / 'vendor' / 'pyodine'

if str(VENDOR_DIR) not in sys.path:
    sys.path.insert(0, str(VENDOR_DIR))
