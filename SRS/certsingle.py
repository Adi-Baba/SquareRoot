"""
Certified single-root bracket adapter.

Wraps the exact-integer certified square root from ../Src/exp_certified.py.

WHAT "CERTIFIED" MEANS HERE
---------------------------
`exp_certified.certify(lo, hi, alpha, digits)` decides
lo <= sqrt(alpha) <= hi and hi - lo <= 10**-digits in EXACT INTEGER
ARITHMETIC: it decomposes each Decimal endpoint as (coeff, exponent) and
compares squares, so no rounding can decide a bound in its own favour.
That is the whole reason this module exists.

`certified_sqrt` below does NOT run that check, because the check costs more
than the bracket does.  `certified_sqrt_certified` and
`get_certified_bracket` do run it and report the result.  Callers that need
a proof should use one of those two; callers that only need digits should
use `certified_sqrt` and must not describe the result as certified.
"""

from decimal import Decimal as D
from typing import Dict, Tuple

from _srclib import load_src_module

# Load the certified core by explicit path rather than mutating sys.path.
# See _srclib for why: the sys.path version worked at run time but was
# invisible to static analyzers, which reported
#   Import "exp_certified" could not be resolved
# for a line that always succeeded.  Loading by path also makes this module
# independent of the current working directory.
EC = load_src_module("exp_certified")

# Re-export the certified core.  `ALPHAS` is re-exported rather than restated
# as a literal: a second copy at the bottom of this file shadowed this one, and
# the two happened to agree only until one of them was edited.
zoom = EC.zoom
bisection = EC.bisection
interval_newton = EC.interval_newton
newton_certified = EC.newton_certified
certify = EC.certify
ALPHAS = EC.ALPHAS


def _check_alpha(alpha: int) -> int:
    if not isinstance(alpha, int) or isinstance(alpha, bool):
        raise TypeError(f"alpha must be an int, got {alpha!r}")
    if alpha < 0:
        raise ValueError(f"alpha must be >= 0, got {alpha}")
    return alpha


def _check_digits(digits: int) -> int:
    if not isinstance(digits, int) or isinstance(digits, bool):
        raise TypeError(f"digits must be an int, got {digits!r}")
    if digits <= 0:
        raise ValueError(f"digits must be positive, got {digits}")
    return digits


def certified_sqrt(alpha: int, digits: int) -> Tuple[D, D, int]:
    """
    Certified bracket for sqrt(alpha) to `digits` decimal places.

    Returns (lo, hi, depth) with lo <= sqrt(alpha) <= hi and
    hi - lo <= 10**-digits, as Decimals.  `depth` is the number of zoom
    levels used.

    The bracket comes from the cubic zoom in exp_certified, whose outward
    rounding is designed so that the result passes certify().  This function
    does not itself call certify(); see the module docstring.
    """
    _check_alpha(alpha)
    _check_digits(digits)
    return EC.zoom(alpha, digits)


def certified_sqrt_certified(alpha: int, digits: int) -> Tuple[D, D, int]:
    """
    Certified bracket for sqrt(alpha), with the proof actually run.

    Raises ValueError if certify() rejects the bracket, so a returned value
    is one that has been verified rather than one that is expected to be.
    """
    _check_alpha(alpha)
    _check_digits(digits)
    lo, hi, depth = EC.zoom(alpha, digits)
    msg = EC.certify(lo, hi, alpha, digits)
    if msg is not None:
        raise ValueError(f"certification failed for alpha={alpha} digits={digits}: {msg}")
    return lo, hi, depth


def get_certified_bracket(alpha: int, digits: int) -> Dict:
    """
    Full bracket result including the outcome of the exact proof.

    Keys: alpha, digits, lo, hi, depth, certified, cert_msg.  `certified` is
    True only when certify() returned no message.
    """
    _check_alpha(alpha)
    _check_digits(digits)
    lo, hi, depth = EC.zoom(alpha, digits)
    msg = EC.certify(lo, hi, alpha, digits)
    return {
        'alpha': alpha,
        'digits': digits,
        'lo': lo,
        'hi': hi,
        'depth': depth,
        'certified': msg is None,
        'cert_msg': msg,
    }


def certified_bracket_to_rational(lo: D, hi: D, alpha: int,
                                  digits: int) -> Dict:
    """
    Exact rational bounds for a Decimal bracket, as integers.

    Decimal.as_integer_ratio() is exact, so lo_num/lo_den really is the
    endpoint and not a rounded print of it.  No proof is performed here; the
    returned dict is a re-encoding, and inherits whatever standing the
    original bracket had.
    """
    lo_num, lo_den = lo.as_integer_ratio()
    hi_num, hi_den = hi.as_integer_ratio()
    return {
        'lo_num': lo_num, 'lo_den': lo_den,
        'hi_num': hi_num, 'hi_den': hi_den,
        'alpha': alpha,
        'digits': digits,
    }


if __name__ == "__main__":
    print("=" * 64)
    print("Certified single-root bracket")
    print("=" * 64)
    print("each bracket is re-proved by certify() in exact integer arithmetic\n")

    bad = 0
    for alpha in (2, 3, 12, 17, 999983):
        lo, hi, depth = certified_sqrt(alpha, 100)
        msg = certify(lo, hi, alpha, 100)
        status = "ok" if msg is None else f"FAIL {msg}"
        bad += msg is not None
        print(f"  alpha={alpha:<8} depth={depth:<3} width={hi - lo:.3E}  {status}")

    print(f"\n{'all brackets certified' if not bad else f'{bad} CERT FAILURES'}")
    raise SystemExit(1 if bad else 0)
