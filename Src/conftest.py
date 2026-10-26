"""
pytest / static-analysis configuration for the Src test suite.

WHY THIS FILE EXISTS
--------------------
`Src/` is a flat directory of sibling modules, not an installed package, so a
test doing `import run_all` only resolves when `Src/` is on `sys.path`.  One
test also checks `os.path.exists(mod + ".py")` for every registered harness,
which is relative to the current directory and therefore only passed when
pytest happened to be invoked from inside `Src/`.

Running the suite from the repository root exposed both:

    $ python3 -m pytest Src
    FAILED Src/test_opcount.py::TestRunAllRegistry::test_registered_modules_all_exist
    AssertionError: False is not true : exp_precision.py missing

`exp_precision.py` was never missing; the test was looking in the wrong
directory.  Putting `Src/` on `sys.path` here fixes the import half, and
`chdir`-ing to `Src/` fixes the existence check, so the suite behaves the
same from any working directory.  Nothing here changes what the code
computes.
"""

import os
import sys

_SRC_DIR = os.path.dirname(os.path.abspath(__file__))

if _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

# test_opcount asserts os.path.exists(mod + ".py") for each harness in
# run_all.HARNESSES, using a path relative to the cwd.  Anchor the cwd to Src/
# so that assertion means what it says regardless of where pytest was invoked.
os.chdir(_SRC_DIR)
