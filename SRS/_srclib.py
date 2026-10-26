"""
Load a module from the sibling `Src/` tree by explicit file path.

WHY NOT `sys.path.insert`
-------------------------
Both `certsingle.py` and `_deps.py` used to do this:

    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'Src'))
    import exp_certified

which has three problems:

1. It mutates global interpreter state as a side effect of importing a
   module.  Any later `import` anywhere in the process can now pick up a
   module from `Src/` by accident.
2. It depends on the *current* working directory in a subtle way whenever a
   relative `__file__` is involved, and on import order otherwise.
3. A static analyzer cannot follow it.  The path is computed at run time, so
   the analyzer has no way to know what `exp_certified` resolves to and
   reports `Import "exp_certified" could not be resolved` on a line that
   works fine at run time.  The message is not wrong about the code being
   unclear; it is wrong only about the outcome.

Loading by absolute file path with `importlib` fixes all three: no global
state is touched, the path is computed from `__file__` alone so the working
directory is irrelevant, and the filename is a literal in the source that an
analyzer can resolve.

WHAT IS GIVEN UP
----------------
Nothing that was wanted, as it happens.  The earlier version of this
docstring claimed the opposite: that the module is loaded under a private
name and NOT registered in `sys.modules`, so a bare `import exp_certified`
elsewhere would build a second copy.  That was never true of this code, which
registers under the plain name *before* executing the module body, and a bare
import afterwards gets the very same object.  `_deps.py` described that
behaviour correctly while this file contradicted it.

Early registration is not incidental, it is load-bearing:

* Registering first means a module that imports itself, or a sibling in the
  same tree, during its own body execution resolves to the object being
  built rather than triggering a second execution of the same file.
* Registering first means a later bare import cannot produce a *second live
  copy*.  `Src/exp_certified.py` carries Decimal context state and calls
  `sys.set_int_max_str_digits(0)` at import time, so two copies would mean
  two independent pieces of interpreter-global state that could disagree.
* Returning the cached object when the name is already present is what makes
  repeated loads idempotent.
"""

import importlib.util
import os
import sys
from types import ModuleType

_SRS_DIR = os.path.dirname(os.path.abspath(__file__))
_SRC_DIR = os.path.abspath(os.path.join(_SRS_DIR, os.pardir, "Src"))


def src_path(module_name: str) -> str:
    """Absolute path of `module_name` inside the sibling `Src/` tree."""
    return os.path.join(_SRC_DIR, f"{module_name}.py")


def load_src_module(module_name: str) -> ModuleType:
    """
    Import `module_name` from the sibling `Src/` directory.

    `module_name` must be a bare module name.  It is interpolated into a path
    here, so without a check a value like "../../etc/passwd" would resolve
    outside `Src/` and execute an arbitrary file.  In-tree callers only ever
    pass literal names, but a loader is the wrong place to rely on that.

    Raises ImportError with the concrete path that was tried, so a failure
    says which file was missing rather than just that a name did not
    resolve.
    """
    if not module_name or os.path.basename(module_name) != module_name:
        raise ImportError(
            f"cannot load {module_name!r}: expected a bare module name with no "
            f"path separators, so the result stays inside {_SRC_DIR!r}"
        )

    path = src_path(module_name)
    if not os.path.isfile(path):
        raise ImportError(
            f"cannot load {module_name!r}: expected it at {path!r}. "
            f"This tree expects SRS/ and Src/ to be siblings under the "
            f"repository root; _SRC_DIR resolved to {_SRC_DIR!r}."
        )

    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"no import spec could be built for {path!r}")

    # If this module is already loaded, hand back the existing object.  Without
    # this, a second caller gets a freshly executed copy while sys.modules keeps
    # pointing at the first, so the two are different objects carrying different
    # Decimal context state -- the exact duplication this module exists to avoid.
    if module_name in sys.modules:
        return sys.modules[module_name]

    module = importlib.util.module_from_spec(spec)

    # Register before exec_module.  This is the documented importlib recipe: it
    # makes a circular import resolve to the partially-built module instead of
    # recursing, and it guarantees the object being executed is the one the rest
    # of the process will see.
    sys.modules[module_name] = module

    # exec_module runs the module body, which may do `from exp_precision
    # import ...`.  Those are sibling imports inside Src/, so Src/ must be
    # importable for the duration.  It is added and removed around the call
    # rather than left in place, keeping the effect scoped to this one load
    # instead of permanently altering global interpreter state.
    added = _SRC_DIR not in sys.path
    if added:
        sys.path.insert(0, _SRC_DIR)
    try:
        spec.loader.exec_module(module)
    except Exception as exc:
        # Do not leave a half-built module behind: a later `import` would find
        # it in sys.modules and inherit a broken object.
        sys.modules.pop(module_name, None)
        raise ImportError(f"loading {path!r} failed: {exc}") from exc
    finally:
        if added:
            try:
                sys.path.remove(_SRC_DIR)
            except ValueError:  # pragma: no cover - defensive
                pass
    return module


if __name__ == "__main__":  # pragma: no cover
    import os as _os

    print(f"_SRC_DIR = {_SRC_DIR}")
    print(f"exists   = {_os.path.isdir(_SRC_DIR)}")
    for name in ("exp_certified", "exp_precision", "exact_decide"):
        print(f"  {name:<16} -> {src_path(name)}")
