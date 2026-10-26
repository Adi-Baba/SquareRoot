"""
Bracket-Zoom Square Root -- reference implementation of Draft1.md.

Section references are to Docs/Draft1.md.

Sign conventions
----------------
p = alpha - a**2  is the "deficit" (Section 2).  With a bracket a < r < b and
w = b - a, one always has p > 0 and p = (r - a)(r + a), which is the
factorisation the whole method rests on (Section 5).

The routines are numeric-type agnostic: pass int, Fraction, or Decimal and
they work unchanged.  Use Fraction when you want the rigorous claims checked
with no floating point in the way; use Decimal when the target precision is
high enough that exact rationals would be too slow (they blow up badly past
a few hundred digits).
"""

from __future__ import annotations

__all__ = [
    "plain",
    "corrected",
    "error_bound",
    "error_bound_tight",
    "closed_form_error_from_r",
    "zoom_decimal",
    "zoom_recursive",
    "ZOOM_OK",
    "ZOOM_DIVERGES",
]

# eps = w**3 / (16 * a**2) is < w  exactly when  w < 4a.  Outside that range
# the "zoom" *widens* the bracket and the invariant a < r < b is destroyed.
ZOOM_OK = "ok"
ZOOM_DIVERGES = "diverges"


def plain(a, b, alpha):
    """Level-0 / plain estimate M(a, b) = a + p/(a + b)   (Sections 2 and 4).

    This is the same object as the paper's "a + w*p/d" since
    d = b**2 - a**2 = w*(a + b), so p/(a+b) = w*p/d.

    Theorem 1 gives r - M(a,b) = (r-a)(b-r)/(a+b) >= 0, so this estimate is
    always a *under*-estimate of sqrt(alpha).
    """
    p = alpha - a * a
    return a + p / (a + b)


def corrected(a, b, alpha):
    """Algorithm 2: est = a + p/(a + m) with m = M(a, b)   (Section 5).

    Theorem 1 with z = m < r gives est > r always: a one-sided over-estimate.
    """
    p = alpha - a * a
    m = a + p / (a + b)
    return a + p / (a + m)


def error_bound(a, w):
    """Theorem 2 as printed: 0 <= est - r <= w**3/(16*a**2).

    Valid, but loose: the true worst case is only 16/27 of this.
    """
    return w ** 3 / (16 * a * a)


def error_bound_tight(a, w):
    """Corrected Theorem 2: 0 <= est - r <= w**3/(27*a**2).

    u**2*(w-u) is maximised at u = 2w/3 with value 4*w**3/27, and the
    denominator 2a(2a+w) + u(2a+u) is >= 4a**2, giving 4*w**3/(27*4a**2).
    Sharpens the paper's 1/16 by a factor of 27/16 = 1.6875.
    """
    return w ** 3 / (27 * a * a)


def closed_form_error_from_r(a, w, r):
    """Exact error  est - r  in closed form.

        est - r = u**2 * (w - u) / ( 2a(2a+w) + u(2a+u) ),   u = r - a, w = b - a

    Obtained by substituting p = u(2a+u) into corrected().  This closed form is
    what Theorem 2's bound is bounding, and it makes the 4/27 constant and the
    u = 2w/3 maximiser visible.  Takes r directly so that u is known exactly.
    """
    u = r - a
    denom = 2 * a * (2 * a + w) + u * (2 * a + u)
    return u * u * (w - u) / denom


def _scan_digit(a, w, alpha):
    """Algorithm 1 step 1: largest k in {0..9} with f(a + k*w) keeping the sign
    of f(a), where f(x) = x**2 - alpha.  Scans downward, so it exits early.
    """
    for k in range(9, -1, -1):
        if (a + k * w) ** 2 <= alpha:
            return k
    return 0


def zoom_decimal(alpha, a, b, levels):
    """Algorithm 1: fixed-factor decimal zoom.  Returns the plain estimates.

    Each level shrinks the bracket by 10 regardless of how good the estimate
    already is -- this is the bottleneck Section 6 set out to fix.
    """
    a = a + 0
    out = []
    for _ in range(levels):
        w = b - a
        w = w / 10
        a = a + _scan_digit(a, w, alpha) * w
        b = a + w
        out.append(plain(a, b, alpha))
    return out


def zoom_recursive(alpha, a, b, tol_digits, tight=False, max_depth=64):
    """Algorithm 3: recursive cubic zoom.  Returns (a, b, depth).

    Contract: caller must supply a < sqrt(alpha) < b AND w = b - a <= 4a.
    The second condition is NOT implied by the first and is not enforced in
    the paper's listing -- without it eps >= w, the bracket widens, and the
    invariant is lost within two steps.  `a >= 1` (as in the paper) makes this
    automatic only while w <= 1.

    `tight=True` swaps Theorem 2 for the 4/27 version.
    """
    a, b = a + 0, b + 0
    tol = type(a)(10) ** -tol_digits
    depth = 0
    while b - a > tol and depth < max_depth:
        w = b - a
        est = corrected(a, b, alpha)
        eps = error_bound_tight(a, w) if tight else error_bound(a, w)
        a, b = est - eps, est
        depth += 1
    return a, b, depth


def zoom_status(a, b):
    """Classify whether a bracket will make progress under Algorithm 3."""
    a, b = a + 0, b + 0
    if a <= 0:
        return ZOOM_DIVERGES  # bound is undefined / infinite
    w = b - a
    if w <= 0:
        return ZOOM_DIVERGES
    return ZOOM_OK if w < 4 * a else ZOOM_DIVERGES
