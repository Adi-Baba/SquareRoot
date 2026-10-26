"""The n-th root zoom map: exact error identity, multiplier, and the obstruction.

Supports Section 10 of `Docs/Paper2.md`.  Everything here is `fractions.Fraction`
and `int`; there is no floating point, no sampling, and no search.  A claim
either holds identically in Q, or the test fails.

Setup.  For alpha > 0, r = alpha^(1/n), and an anchor a with 0 < a < r, put

    p_n   = r^n - a^n
    S_n(z) = sum_{k=0}^{n-1} z^(n-1-k) a^k          (degree n-1 in z)
    phi_n(z) = a + p_n / S_n(z)                       (degree n-1 rational)

so r = a + p_n / S_n(r): the n-th root analogue of r = a + p/(a+r).

The two facts this module establishes exactly:

  * the error identity  r - phi_n(z) = (r-a)(z-r) Q_n(z,r,a) / S_n(z)  with
    Q_n the complete homogeneous symmetric polynomial h_{n-2}(z, r, a); and
  * the multiplier  lambda_n = phi_n'(r) = -(r-a) S_n'(r) / S_n(r)
                           = 1 - n / sum_{k<n} (a/r)^k.

Q_2 = 1, which is precisely why n = 2 is the only case in which the error is a
ratio of *linear* forms and the map is Moebius.  For n >= 3, Q_n has degree
n-2 >= 1, and no rational change of coordinates can turn a degree-(n-1) map
into a linear scaling (see the degree argument in the paper).

A caveat that is easy to get wrong and is pinned by a test here: the scalar
map is NOT attracting for every bracket.  |lambda_n| < 1 holds only while
sum_{k<n} t^k > n/2 with t = a/r, i.e. only for a sufficiently tight anchor;
as the anchor weakens lambda_n -> -(n-1), so the fixed point becomes repelling
for n >= 3.  See `attracting_threshold`.
"""

from fractions import Fraction as F


# ------------------------------------------------------------- the map itself


def S(z, a, n):
    """S_n(z, a) = sum_{k=0}^{n-1} z^(n-1-k) a^k."""
    z, a = F(z), F(a)
    return sum(z ** (n - 1 - k) * a ** k for k in range(n))


def dS(z, a, n):
    """dS_n/dz at z = z (same expression, one degree lower)."""
    z, a = F(z), F(a)
    return sum((n - 1 - k) * z ** (n - 2 - k) * a ** k for k in range(n - 1))


def phi(z, a, r, n):
    """phi_n(z) = a + (r^n - a^n)/S_n(z,a)."""
    z, a, r = F(z), F(a), F(r)
    return a + (r ** n - a ** n) / S(z, a, n)


def Q(z, r, a, n):
    """Q_n(z,r,a) = sum_{k=0}^{n-2} sum_{j=0}^{n-2-k} z^(n-2-k-j) r^j a^k.

    This is the complete homogeneous symmetric polynomial h_{n-2}(z,r,a);
    `test_mobius_n.py` proves the identification by expansion.
    """
    z, r, a = F(z), F(r), F(a)
    return sum(z ** (n - 2 - k - j) * r ** j * a ** k
               for k in range(n - 1) for j in range(n - 1 - k))


def h_complete(vars_, degree):
    """The complete homogeneous symmetric polynomial h_degree in `vars_`."""
    out = F(0)
    # enumerate all exponent tuples summing to `degree`
    for e0 in range(degree + 1):
        for e1 in range(degree - e0 + 1):
            rest = degree - e0 - e1
            term = F(1)
            for v, e in zip(vars_, (e0, e1, rest)):
                term *= v ** e
            out += term
    return out


# ---------------------------------------------------------- the two theorems


def error_identity(z, a, r, n):
    """r - phi_n(z) == (r-a)(z-r) Q_n(z,r,a) / S_n(z,a)  -- exact, returns bool.

    The identity is decidable in Q, so this is a check on exact values, not a
    sample.  `test_mobius_n.py` drives it over a spread of (z,a,r,n).
    """
    z, a, r = F(z), F(a), F(r)
    lhs = r - phi(z, a, r, n)
    rhs = (r - a) * (z - r) * Q(z, r, a, n) / S(z, a, n)
    return lhs == rhs


def multiplier(a, r, n):
    """lambda_n = phi_n'(r), in the form that is useful for the threshold."""
    a, r = F(a), F(r)
    t = a / r
    return 1 - F(n) / sum(t ** k for k in range(n))


def multiplier_alt(a, r, n):
    """The same number as -(r-a) S_n'(r) / S_n(r); equal, proved in tests."""
    a, r = F(a), F(r)
    return -(r - a) * dS(r, a, n) / S(r, a, n)


# --------------------------------------------------------------- the threshold


def attracting(a, r, n):
    """Is the fixed point r attracting for the scalar map phi_n?

    Equivalently |lambda_n| < 1, i.e. sum_{k=0}^{n-1} t^k > n/2 with t = a/r.
    For n = 2 this holds for every t in (0,1); for n >= 3 it fails once the
    anchor is too far from the root.
    """
    a, r = F(a), F(r)
    t = a / r
    return sum(t ** k for k in range(n)) > F(n, 2)


def threshold_gap(n):
    """1 - t_crit: the relative accuracy u/r you must already have.

    Returns None for n = 2 (no requirement: the map always contracts).  For
    n >= 3, t_crit solves sum_{k<n} t^k = n/2, and the returned value is
    1 - t_crit, which tends to 1.5936.../n as n grows.

    Returned as a float: the boundary is located by bisection, and only its
    asymptotic behaviour is used.  Use `attracting` for exact decisions.
    """
    if n == 2:
        return None
    # Locate the boundary to ~50 bits in float, then return the exact rational
    # bracket.  Doing the whole search in Fraction is exact but hopeless: 200
    # halvings on exact rationals makes every power t^k a 2^200-bit numerator,
    # and threshold_gap(500) took minutes.
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2.0
        if sum(mid ** k for k in range(n)) < n / 2.0:
            lo = mid
        else:
            hi = mid
    return 1 - (lo + hi) / 2.0


# ------------------------------------------------- the n = 2 recovery checks

T = lambda z, r: (F(z) - F(r)) / (F(z) + F(r))     # noqa: E731  conjugacy


def mobius_2(z, a, r):
    """phi_2 written as (az + alpha)/(z + a) -- the Moebius form of Paper2 §4."""
    z, a, r = F(z), F(a), F(r)
    return (a * z + r * r) / (z + a)
