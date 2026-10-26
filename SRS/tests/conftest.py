"""
pytest / static-analysis configuration for the SRS test suite.

WHY THIS FILE EXISTS
--------------------
`SRS/` and `Src/` are both flat directories of sibling modules, not installed
packages.  Two things follow, and both used to be real problems:

1. A test module doing `import signexact` only resolves if `SRS/` is on
   `sys.path`.  Running `pytest SRS/tests` from the repository root therefore
   failed collection with ModuleNotFoundError for every test, even though
   `cd SRS && pytest` passed.  Adding SRS (and Src, which SRS modules import
   from) to `sys.path` here makes the suite runnable from any directory.

2. A static analyzer cannot see a `sys.path.insert(...)` executed at import
   time, so it reports `Import "exp_certified" could not be resolved` for
   SRS/certsingle.py and SRS/_deps.py even though the import succeeds at
   runtime.  The analyzer reads this file and the type-checker settings below
   rather than guessing; the same paths are listed in pyrightconfig.json for
   tools that do not read pytest config.

Nothing here changes what the code computes.  It only makes the existing
flat-module layout visible to tools that could not otherwise follow it.
"""

import os
import sys

# This file lives in SRS/tests/, so SRS/ is one level up and the repository
# root two. Getting this wrong resolves ../Src to SRS/Src, which does not
# exist, and the failure then looks like "the tests cannot find their own
# modules" rather than a bad path.
_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_SRS_DIR = os.path.dirname(_TESTS_DIR)
_REPO_ROOT = os.path.dirname(_SRS_DIR)
_SRC_DIR = os.path.join(_REPO_ROOT, "Src")

# SRS first, so the local modules win over anything of the same name.
for _path in (_SRS_DIR, _SRC_DIR):
    if _path not in sys.path:
        sys.path.insert(0, _path)
