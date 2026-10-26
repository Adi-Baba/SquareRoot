"""
Certified rational interval arithmetic for square-root sums.

Gives a certified bracket [lo, hi] for S = SUM_i c_i sqrt(a_i) and a
certified sign for it.

SOUNDNESS
---------
The interval kernel is the integer one in signexact: for each radicand it
computes t = isqrt(a << 2b), which certifies t^2 <= a*2^(2b) < (t+1)^2,
hence t/2^b <= sqrt(a) <= (t+1)/2^b.  Endpoints are combined with
sign-aware selection, so the enclosure is valid term by term and no
cancellation is assumed.  Every value here is an exact Fraction; nothing
passes through a binary or decimal float, and in particular the earlier
ROUND_FLOOR-with-a-fixed-precision-context approach is gone -- flooring
*both* endpoints is not an enclosure, and a context that is too narrow
silently rounds the answer rather than admitting it cannot be trusted.

WHAT "SLACK" CAN AND CANNOT MEAN HERE
------------------------------------
Slack is normally (bracket width) / (error of the value you meant to
compute).  For a square-root sum there is no separate "intended value" to be
in error against: the exact value is the thing we cannot write down, and the
bracket *is* our description of it.  So the honest quantity is

    relative width = (hi - lo) / max(|lo|, |hi|)

which says how many digits of S the bracket actually pins down, and it is
only meaningful once the bracket has cleared zero.  That is the whole of it:
there is no `slack_ratio` here, and an earlier version of this docstring
referred to one as though callers could still ask for it.  A width/error
ratio needs a true value to be in error against, and this problem has none,
so the ratio was removed rather than filled in with a number that would look
like an error bound and not be one.  Returning `inf` -- which is what the
removed function did -- would have been the dishonest version: it reads as
"infinitely uncertain" and silently satisfies any tolerance check handed to
it.
"""

from fractions import Fraction as F
from typing import Dict, Iterable, NamedTuple, Optional, Tuple

from signexact import (
    DEFAULT_MAX_BITS,
    SignUnknownError,
    _sum_bracket_int,
    certified_bracket,
    exact_sign,
)
from squarefree import decompose_and_aggregate

__all__ = [
    "CertifiedInterval",
    "certified_interval",
    "certified_interval_for_S",
    "certified_interval_for_term",
    "relative_width",
    "sign_with_bracket",
    "sign_with_bracket_for_S",
    "digit_count",
]


class CertifiedInterval(NamedTuple):
    """A certified rational enclosure of S, with the bits it cost."""

    lo: F
    hi: F
    bits: int
    sign: int

    @property
    def width(self) -> F:
        return self.hi - self.lo

    @property
    def mid(self) -> F:
        return (self.lo + self.hi) / 2

    def contains_zero(self) -> bool:
        return self.lo <= 0 <= self.hi

    def relative_width(self) -> F:
        return relative_width(self.lo, self.hi)


def relative_width(lo: F, hi: F) -> F:
    """
    (hi - lo) / max(|lo|, |hi|), exact.

    0 means a degenerate (exact) bracket.  Larger means fewer digits pinned,
    and the value is only meaningful once the bracket has cleared zero.

    The degenerate and straddling cases are spelled out because they are not
    what intuition suggests, and the previous version of this docstring got
    both wrong:

      * [0, 0] returns 0, not infinity.  A zero-width bracket is an exactly
        determined value (S = 0), which is the *best* possible outcome, so it
        is the one case that must not be reported as total uncertainty.
      * A bracket straddling zero symmetrically, [-1, 1], returns 2, not 1.
        There is no reason for it to be exactly 1; the value is simply large.

    An empty interval (lo > hi) is a caller bug and raises.
    """
    if lo > hi:
        raise ValueError(f"empty interval: lo={lo} > hi={hi}")
    scale = max(abs(lo), abs(hi))
    if scale == 0:
        return F(0)          # lo == hi == 0, a degenerate exact bracket
    return (hi - lo) / scale


def digit_count(interval: CertifiedInterval) -> Optional[int]:
    """
    Decimal digits of S pinned by the bracket, or None if it straddles zero.

    None is not an error: while the bracket contains 0 no digit of S is
    determined, and reporting 0 would falsely suggest "exactly 0".
    """
    if interval.contains_zero():
        return None
    rel = interval.relative_width()
    if rel <= 0:
        return None
    # -log10(rel), computed on integers to stay exact.
    num, den = rel.numerator, rel.denominator
    # rel > 1/10^k  <=>  den < num * 10^k
    digits = 0
    scale = num
    while den >= scale:
        scale *= 10
        digits += 1
    return digits


def certified_interval(
    agg: Dict[int, int],
    bits: int,
    max_bits: int = DEFAULT_MAX_BITS,
) -> CertifiedInterval:
    """
    Certified bracket for an already-aggregated sum, refining until it clears 0.

    `agg` maps squarefree part -> integer coefficient.  Returns a
    CertifiedInterval whose lo <= S <= hi exactly, with the sign of S when the
    bracket has cleared zero.  Raises SignUnknownError if the cap is hit.
    """
    if bits <= 0:
        raise ValueError(f"bits must be positive, got {bits}")

    work = {s: c for s, c in agg.items() if c != 0}
    if not work:
        return CertifiedInterval(F(0), F(0), 0, 0)
    if len(work) == 1 and 1 in work:
        total = F(work[1])
        return CertifiedInterval(total, total, 0, (total > 0) - (total < 0))

    b = bits
    while b <= max_bits:
        lo_s, hi_s = _sum_bracket_int(work, b)
        denom = 1 << b
        lo, hi = F(lo_s, denom), F(hi_s, denom)
        if lo > 0:
            return CertifiedInterval(lo, hi, b, 1)
        if hi < 0:
            return CertifiedInterval(lo, hi, b, -1)
        b *= 2

    raise SignUnknownError(max_bits, sum(abs(c) for c in work.values()).bit_length() - b,
                           "raise max_bits to continue")


def certified_interval_for_S(
    terms: Iterable[Tuple[int, int]],
    digits: int = 50,
    max_bits: Optional[int] = None,
) -> CertifiedInterval:
    """
    Certified bracket for S = SUM (coeff, radicand) terms.

    `digits` is a hint for the starting precision; the interval is refined
    until it clears zero regardless, so a too-small hint costs a little time
    and never a wrong answer.
    """
    materialised = list(terms)
    for coeff, rad in materialised:
        if not isinstance(coeff, int) or not isinstance(rad, int):
            raise TypeError(f"terms must be (int, int), got ({coeff!r}, {rad!r})")
        if rad <= 0:
            raise ValueError(f"radicand must be positive, got {rad}")

    agg, _ = decompose_and_aggregate(materialised)
    start = 16 + 2 * max((s.bit_length() for s in agg), default=0) + 3 * digits
    cap = DEFAULT_MAX_BITS if max_bits is None else max_bits
    return certified_interval(agg, start, cap)


def certified_interval_for_term(coeff: int, radicand: int, digits: int = 50) -> Tuple[F, F]:
    """
    Certified bracket for coeff * sqrt(radicand), as exact Fractions.

    Sign-aware: a negative coefficient swaps the endpoints.
    """
    lo, hi = certified_bracket({radicand: 1}, 16 + 2 * radicand.bit_length() + 3 * digits)
    c = F(coeff)
    if coeff >= 0:
        return (c * lo, c * hi)
    return (c * hi, c * lo)


def sign_with_bracket(agg: Dict[int, int], **kwargs) -> CertifiedInterval:
    """Certified bracket for an aggregated sum; see `certified_interval`."""
    return certified_interval(agg, **kwargs)


def sign_with_bracket_for_S(
    terms: Iterable[Tuple[int, int]], digits: int = 50, max_bits: Optional[int] = None
) -> CertifiedInterval:
    """Certified bracket and sign for a term list; see `certified_interval_for_S`."""
    return certified_interval_for_S(terms, digits, max_bits)


# Backwards-compatible alias: the sign alone, as a SignResult.
def certified_sign(terms: Iterable[Tuple[int, int]], **kwargs) -> int:
    return exact_sign(terms, **kwargs).sign


if __name__ == "__main__":
    print("=" * 64)
    print("Certified interval arithmetic for S")
    print("=" * 64)

    cases = [
        ([(1, 12), (-1, 3)], "sqrt(12) - sqrt(3) = sqrt(3)"),
        ([(1, 2), (-1, 3)], "sqrt(2) - sqrt(3)"),
        ([(1, 2), (1, 3), (-1, 5)], "sqrt(2)+sqrt(3)-sqrt(5)"),
        ([(1, 12), (-2, 3)], "sqrt(12) - 2sqrt(3) = 0 exactly"),
    ]
    for terms, why in cases:
        iv = certified_interval_for_S(terms, 20)
        pinned = digit_count(iv)
        print(f"\n{why}")
        print(f"  sign={iv.sign:+d} bits={iv.bits} digits_pinned={pinned}")
        print(f"  lo={float(iv.lo):.17g}")
        print(f"  hi={float(iv.hi):.17g}")
        if iv.contains_zero():
            print("  (bracket contains zero: no digits determined, sign from the exact zero test)")
