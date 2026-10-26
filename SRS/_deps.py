"""
SRS internal dependency shim.

Loads the sibling `Src/` modules by explicit file path and re-exports the
handful of functions the SRS package uses.

WHY THIS EXISTS
---------------
`SRS/` and `Src/` are flat directories of siblings under the repository root,
not installed packages.  Two mechanisms were used to bridge them and both
were wrong:

* `sys.path.insert(...)` at import time.  It mutated global interpreter
  state, made the result depend on import order, and was invisible to static
  analyzers, which reported `Import "exp_certified" could not be resolved`
  for a line that worked at run time.
* A bare `import exp_certified`, relying on something else having already
  put `Src/` on the path.  This one is worse: it works or not depending on
  which module got imported first.

Both are replaced by `_srclib.load_src_module`, which loads by absolute path
computed from `__file__`.  No global state, no import-order dependency, and
the filename appears literally in the source so a resolver can follow it.

NOTE ON DUPLICATE COPIES
------------------------
`load_src_module` registers each module in `sys.modules` under its plain name
before executing it, so a bare `import exp_certified` anywhere later in the
process finds this same object rather than building a second one.  That
matters because these modules carry Decimal context state: two live copies
would mean two independent pieces of global state that could disagree.
"""

from _srclib import load_src_module, src_path

_cert = load_src_module("exp_certified")
_decide = load_src_module("exact_decide")

# Certified-root core: exact-integer certification of a sqrt bracket.
certify = _cert.certify
zoom = _cert.zoom
# exp_certified names it `bisection`, not `bisect`; the old alias raised
# AttributeError the moment anything imported this module.
bisection = _cert.bisection
interval_newton = _cert.interval_newton
newton_certified = _cert.newton_certified

# Exact decisions about the published / u-aware error bounds.
decide_u_aware_bound = _decide.decide_u_aware_bound
decide_published_bound = _decide.decide_published_bound
decide_error_identity = _decide.decide_error_identity

__all__ = [
    "certify",
    "zoom",
    "bisection",
    "interval_newton",
    "newton_certified",
    "decide_u_aware_bound",
    "decide_published_bound",
    "decide_error_identity",
    "src_path",
]
