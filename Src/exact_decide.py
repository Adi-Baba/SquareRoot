"""Exact decision procedures for the paper's algebraic claims.

Why this file exists
--------------------
Several claims in `Docs/Paper1.md` were originally supported by *sampling* --
run N random brackets, observe zero violations, conclude the bound holds.  That
is not a proof: it shows no counterexample was found in the sample, which is a
statement about the sample and not about the theorem.  The u-aware bound
(Theorem 8) is one such claim, and "0/30000" was doing work it cannot do.

The good news is that these claims are *decidable*.  Every one of them reduces
to a polynomial inequality on a closed interval with rational coefficients:

    f(u) <= 0   for all u in [lo, hi],   f cubic, lo/hi/f rational

and the truth of that is decidable by exact algebra with no search, no
floating point, and no sampling.  A cubic attains its maximum on an interval at
an endpoint or at a critical point, and its critical points are the roots of a
quadratic.  So the maximum is attained either at rational points (evaluate
exactly) or at quadratic irrationals of the form (p + q*sqrt(D))/r (evaluate
exactly in Q(sqrt(D)), which is a two-dimensional vector space over Q and
therefore closed under arithmetic).

That is the entire argument.  `poly_nonpos_on_interval` implements it, and
every "verified over N random cases" claim in the paper can be re-expressed as
"decided exactly", with the decider itself pinned by this file's tests.

Everything here is `fractions.Fraction` and `math.isqrt`.  There is no
floating point and no randomness anywhere in this module, by construction.
"""

from fractions import Fraction as F
from math import isqrt

__all__ = [
    "poly_at",
    "quad_roots_exact",
    "poly_nonpos_on_interval",
    "decide_u_aware_bound",
    "decide_published_bound",
    "decide_error_identity",
    "sqrt_field",
]


# ---------------------------------------------------------------- polynomials
#
# Polynomials are coefficient lists in DESCENDING order: [c_n, ...,c_1, c_0]
# so that p(u) = sum c_i u^(n-i).  All coefficients are Fractions.


def poly_at(p, u):
    """Evaluate descending-coefficient list p = [c_n,...,c_0] at rational u.

    Horner must walk the list in its stored (descending) order, starting from
    the leading coefficient.  Iterating `reversed(p)` instead evaluates
    sum c_i u^i with the exponents attached backwards -- a different
    polynomial.  That bug shipped in the first draft of this file and was
    caught by the brute-force cross-check in test_exact_decide.py; both this
    and `_eval_at` are now pinned against direct expansion.
    """
    s = F(0)
    for c in p:
        s = s * u + c
    return s


def _derivative(p):
    """d/du of a descending list, dropping the (zero) leading term."""
    n = len(p) - 1
    if n <= 0:
        return [F(0)]
    return [p[i] * (n - i) for i in range(n)]


def _is_nonneg_sqrt(d):
    """True if rational d >= 0 is a square of a rational, i.e. sqrt(d) in Q."""
    if d < 0:
        return False
    dn, dd = isqrt(d.numerator), isqrt(d.denominator)
    return dn * dn == d.numerator and dd * dd == d.denominator


def quad_roots_exact(b, c):
    """Real roots of x^2 + b x + c, exactly.

    Returns a list of 0, 1 or 2 entries.  Rational roots come back as
    Fractions.  An irrational root comes back as a `sqrt_field` element, which
    is an exact representation (p + q*sqrt(D))/r with D a positive rational
    that is not a square.  Nothing is approximated.
    """
    disc = b * b - 4 * c
    if disc < 0:
        return []
    if _is_nonneg_sqrt(disc):
        s = F(isqrt(disc.numerator), isqrt(disc.denominator))
        r1, r2 = (-b + s) / 2, (-b - s) / 2
        return sorted([r1, r2])
    # irrational: represent as (p + q*sqrt(D))/r with D > 0 non-square
    return sorted([sqrt_field(-b / 2, F(1) / 2, disc),
                   sqrt_field(-b / 2, F(-1) / 2, disc)],
                  key=lambda z: z)


# ------------------------------------------------------- exact sqrt(D) field
#
# A number of the form (p + q*sqrt(D))/r  with p, q, r rational and D a
# positive rational that is not a perfect square.  Arithmetic is componentwise
# in the basis {1, sqrt(D)}; comparison to zero is exact because
# p + q*sqrt(D) has the sign of p unless p and q disagree, in which case it is
# the sign of (p^2 - q^2 D) with p's sign.  No floating point is involved.


class sqrt_field:
    """An exact element of Q(sqrt(D)).  Arithmetic and sign are exact."""

    __slots__ = ("p", "q", "D", "_r")

    def __init__(self, p, q, D):
        # Store as (p + q*sqrt(D)) with a common denominator folded into p, q.
        self.p = F(p)
        self.q = F(q)
        self.D = F(D)
        self._r = F(1)

    # --- helpers
    def _narrow(self):
        """Absorb the denominator: (p + q*sqrt D)/r == (p/r) + (q/r)*sqrt D."""
        return self

    def as_pair(self):
        """(p, q) with value p + q*sqrt(D).  D stripped if q == 0."""
        if self.q == 0:
            return (self.p, F(0), F(1))
        return (self.p, self.q, self.D)

    def is_rational(self):
        return self.q == 0

    def to_fraction(self):
        if self.q != 0:
            raise ValueError("irrational element has no Fraction value")
        return self.p

    # --- arithmetic
    def __add__(self, other):
        return self._combine(other, 1)

    def __sub__(self, other):
        return self._combine(other, -1)

    def _combine(self, other, sign):
        if isinstance(other, F):
            return sqrt_field(self.p + sign * other, self.q, self.D)
        if isinstance(other, int):
            return sqrt_field(self.p + sign * F(other), self.q, self.D)
        if isinstance(other, sqrt_field):
            # both must share a radicand to add directly; if not, square one
            if self.D == other.D:
                return sqrt_field(self.p + sign * other.p, self.q + sign * other.q,
                                  self.D)
            raise ValueError("mismatched radicands; square first")
        return NotImplemented

    def __mul__(self, other):
        if isinstance(other, (F, int)):
            return sqrt_field(self.p * F(other), self.q * F(other), self.D)
        if isinstance(other, sqrt_field):
            # (p1+q1 s)(p2+q2 s) = (p1 p2 + q1 q2 D) + (p1 q2 + q2 p1) s
            p = self.p * other.p + self.q * other.q * other.D
            q = self.p * other.q + self.q * other.p
            return sqrt_field(p, q, self.D)
        return NotImplemented

    __rmul__ = __mul__

    def __pow__(self, n):
        if not isinstance(n, int) or n < 0:
            raise ValueError("only non-negative integer powers supported")
        result = sqrt_field(F(1), F(0), self.D)
        base = self
        while n:
            if n & 1:
                result = result * base
            base = base * base
            n >>= 1
        return result

    def __eq__(self, other):
        return self._cmp(other) == 0

    def __le__(self, other):
        return self._cmp(other) <= 0

    def __ge__(self, other):
        return self._cmp(other) >= 0

    def __lt__(self, other):
        return self._cmp(other) < 0

    def __gt__(self, other):
        return self._cmp(other) > 0

    def _diff_with(self, other):
        """Return (p, q, D) for (self - other), as p + q*sqrt(D).

        Handles the two comparisons that actually arise in this file:
        against a Fraction, and against another surd with the same radicand.
        Different radicands raise -- they cannot arise, because every critical
        point handed to us comes from the single quadratic derivative p', so
        they all share one D.
        """
        if isinstance(other, (F, int)):
            return (self.p - F(other), self.q, self.D)
        if isinstance(other, sqrt_field):
            if self.D != other.D:
                raise ValueError(
                    "cannot order surds with different radicands; this cannot "
                    "occur in poly_nonpos_on_interval (all critical points "
                    "come from one quadratic)")
            return (self.p - other.p, self.q - other.q, self.D)
        return NotImplemented

    def _cmp(self, other):
        d = self._diff_with(other)
        if d is NotImplemented:
            return NotImplemented
        return sqrt_field(d[0], d[1], d[2]).sign()

    def sign(self):
        """Exact sign of p + q*sqrt(D).  Returns -1, 0 or 1."""
        if self.q == 0:
            return (self.p > 0) - (self.p < 0)
        # sign of p + q*sqrt(D): same as sign of p + q*sqrt(D) squared-ish
        # Standard trick: compare |p| vs |q*sqrt D| using p^2 vs q^2 D.
        lhs = self.p * self.p
        rhs = self.q * self.q * self.D
        if self.p == 0:
            return (self.q > 0) - (self.q < 0)
        if lhs == rhs:
            return 0
        if lhs > rhs:            # |p| dominates
            return (self.p > 0) - (self.p < 0)
        # |q sqrt D| dominates -> sign of q*sqrt(D) dominates -> sign of q
        return (self.q > 0) - (self.q < 0)

    def __repr__(self):
        if self.q == 0:
            return f"Q({self.p})"
        return f"Q({self.p} + {self.q}*sqrt({self.D}))"

    def value(self):
        """For tests/printing only: exact-ish float, never used in proofs."""
        from decimal import Decimal, getcontext
        getcontext().prec = 60
        d = Decimal(self.D.numerator) / Decimal(self.D.denominator)
        s = d.sqrt()
        return float(Decimal(self.p.numerator) / Decimal(self.p.denominator)
                     + Decimal(self.q.numerator) / Decimal(self.q.denominator) * s)


# ------------------------------------------------------------ the decider


def poly_nonpos_on_interval(p, lo, hi):
    """Decide exactly: is p(u) <= 0 for ALL u in [lo, hi]?

    p is a descending coefficient list over Fractions.  lo, hi are Fractions.

    Method: a cubic attains its maximum on a closed interval at an endpoint or
    at a critical point (Fermat).  Critical points are the real roots of p',
    which is quadratic, so they are rational or quadratic irrational and are
    found exactly by `quad_roots_exact`.  We evaluate p at the two endpoints
    and at every critical point that lies inside; the answer is "yes" iff all
    of those values are <= 0.  Exact, finite, no search, no sampling.

    Returns (bool, witness) where witness is a point (Fraction or sqrt_field)
    at which p > 0 if the answer is False, else None.  The witness makes a
    negative answer falsifiable: a reader can check p(witness) > 0 by hand.
    """
    lo, hi = F(lo), F(hi)
    if lo > hi:
        lo, hi = hi, lo
    for u in (lo, hi):
        if poly_at(p, u) > 0:
            return False, u
    for c in quad_roots_exact(*_quad_form(_derivative(p))):
        if lo <= c <= hi:
            val = _eval_at(p, c)
            if val > 0:
                return False, c
    return True, None


def _quad_form(q):
    """Normalise a descending quadratic [c2, c1, c0] to the form x^2 + b x + c.

    Dividing by the leading coefficient c2 gives
        x^2 + (c1/c2) x + c0/c2,
    so b = c1/c2 and c = c0/c2 -- both WITHOUT a sign change.  Negating b
    (as an earlier version did) reflects the roots through the origin and
    made the decider miss the true maximiser entirely: it accepted a bound
    half the true maximum.  Pinned by test_exact_decide.py.
    """
    if len(q) != 3:
        raise ValueError("expected quadratic")
    c2, c1, c0 = q
    if c2 == 0:
        raise ValueError("degenerate")
    return (c1 / c2, c0 / c2)


def _eval_at(p, c):
    """Evaluate descending-coefficient list p at a Fraction or sqrt_field.

    Same forward Horner as `poly_at`; kept separate only so the type of `c`
    is documented at the call site.
    """
    s = F(0)
    for coef in p:
        s = s * c + coef
    return s


# --------------------------------------------- the paper's claims, decided


def decide_u_aware_bound(a, w, lo, hi, B):
    """Decide exactly: is n(u)/D(u) <= B for all u in [lo, hi]?

    n(u) = u^2 (w - u),  D(u) = 2a(2a+w) + u(2a+u),  D(u) > 0 on [lo,hi]
    because a > 0 and u >= 0.  So the claim is equivalent to

        n(u) - B*D(u) <= 0   for all u in [lo, hi],

    a cubic inequality, which `poly_nonpos_on_interval` decides exactly.

    This is the exact content of Theorem 8's validity claim for the given
    inputs.  Pass the real B (the closed form from `exp_u_aware`) and this
    returns whether B is genuinely a bound, with a witness if it is not.

    Preconditions (asserted, not assumed silently):
      a > 0, w > 0, 0 <= lo <= hi <= w
    """
    a, w, lo, hi, B = F(a), F(w), F(lo), F(hi), F(B)
    assert a > 0 and w > 0, "need a > 0, w > 0"
    assert 0 <= lo <= hi <= w, f"need 0 <= lo <= hi <= w, got {lo},{hi},{w}"
    # p(u) = n(u) - B*D(u) = u^2(w-u) - B*(2a(2a+w) + u(2a+u))
    #     = -u^3 + (w - B) u^2 - 2aB u - 2aB(2a+w)
    # (Verified against direct evaluation in test_exact_decide.py -- an
    # earlier version of this line carried a spurious B^2 u term.)
    p = [F(-1),
         w - B,
         -2 * a * B,
         -2 * a * B * (2 * a + w)]
    return poly_nonpos_on_interval(p, lo, hi)


def decide_published_bound(a, w, B):
    """Decide exactly: is n(u)/D(u) <= B for all u in [0, w]?

    Same polynomial, wider interval.  With B = w^3/(27 a^2) this is Theorem 4.
    """
    return decide_u_aware_bound(a, w, F(0), F(w), B)


def decide_error_identity(a, w, u):
    """Decide exactly: is est - r == u^2(w-u) / (2a(2a+w) + u(2a+u))?

    Both sides are evaluated as exact Fractions -- this identity is a closed
    algebraic form, so it is decided by evaluation, not by search.  Returns
    (True, None) or (False, (lhs, rhs)) with both sides exact so the reader
    can compare them by hand.

    Setup, all in exact rationals:  r = a + u,  so alpha = r^2 = (a+u)^2 and
        p = alpha - a^2 = (a+u)^2 - a^2 = u(2a + u).
    Then  m = a + p/(2a+w)  (the plain estimator, z = b = a+w)  and
          est = a + p/(a + m)  (the corrected estimator).
    """
    a, w, u = F(a), F(w), F(u)
    p = u * (2 * a + u)                 # = alpha - a^2, with r = a + u
    m = a + p / (2 * a + w)             # a + p/(a+b), b = a+w
    est = a + p / (a + m)               # corrected estimate
    lhs = est - (a + u)                 # est - r
    rhs = u * u * (w - u) / (2 * a * (2 * a + w) + u * (2 * a + u))
    return (lhs == rhs, (lhs, rhs) if lhs != rhs else None)


if __name__ == "__main__":
    # small self-demo
    print("poly_nonpos_on_interval decides cubics exactly; see test_certified.py")
