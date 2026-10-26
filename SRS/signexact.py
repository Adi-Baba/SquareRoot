"""
Exact sign decision for square-root sums.

Given  S = SUM_i c_i * sqrt(a_i)  with integer c_i and integer a_i > 0,
decides sign(S) exactly.

WHY THIS IS EXACT
-----------------
The decision path uses *only* Python integers.  No float, no Decimal, no
approximation of any kind enters it.  The kernel is:

    For s >= 0 and b >= 0,  t = isqrt(s << 2b)  satisfies
        t^2 <= s * 2^{2b} < (t + 1)^2
    and t >= 0, so with  lo = t/2^b,  hi = (t+1)/2^b  (both >= 0)
        lo^2 <= s <= hi^2   =>   lo <= sqrt(s) <= hi
        hi - lo = 2^{-b}

Aggregating with sign-aware endpoint selection gives integers L, H with
S in [L/2^b, H/2^b], so sign(S) = sign(L) whenever L > 0 and sign(S) = sign(H)
whenever H < 0.  Both tests are exact integer comparisons.

TERMINATION
-----------
S = 0 if and only if every aggregated coefficient vanishes, because
{ sqrt(d) : d distinct squarefree } is linearly independent over Q.  The
zero test is therefore exact, and once S != 0 the interval width
C * 2^{-b} (C = sum |c_i| over the aggregated coefficients) shrinks below
|S| for all sufficiently large b, so the doubling loop always decides.
If the precision cap is hit first the function raises SignUnknownError
rather than guessing.

THE COST
--------
Deciding the sign is the easy half; the interesting bit is *how many bits*
are needed.  It is governed by how close S sits to 0, which is the
square-root-sum separation problem, and it is exponential in the input
size for adversarial inputs.  This module reports the bit count it used so
callers can see the cost instead of being surprised by it.
"""

from math import isqrt
from typing import Dict, Iterable, Optional, Tuple

from squarefree import decompose_and_aggregate

__all__ = [
    "SignResult",
    "SignUnknownError",
    "exact_sign",
    "sign_of_sum",
    "certified_bracket",
    "certified_bracket_int",
    "DEFAULT_MAX_BITS",
]


# Refuse to spin forever on adversarial inputs.  2**22 fractional bits is far
# past anything that occurs for sums with <= 2**22 total input bits, so hitting
# this cap means the instance is pathological -- and saying so beats returning
# a sign we cannot justify.
DEFAULT_MAX_BITS = 1 << 22


class SignUnknownError(RuntimeError):
    """
    Raised when the sign could not be decided within the precision cap.

    This is deliberately an exception and not a `SignResult`: a caller that
    receives a result object believes it.  Refusing to answer must not look
    like an answer.
    """

    def __init__(self, bits: int, width_log2: int, msg: str = ""):
        self.bits = bits
        self.width_log2 = width_log2
        detail = f": {msg}" if msg else ""
        super().__init__(
            f"sign not decided within {bits} fractional bits "
            f"(bracket still ~2^{width_log2} wide){detail}"
        )


class SignResult:
    """Result of an exact sign determination.  `sign` is -1, 0 or +1."""

    __slots__ = ("sign", "method", "witness")

    def __init__(self, sign: int, method: str, witness: Optional[dict] = None):
        if sign not in (-1, 0, 1):
            raise ValueError(f"sign must be -1, 0 or 1, got {sign}")
        self.sign = sign
        self.method = method
        self.witness = witness if witness is not None else {}

    def __repr__(self) -> str:
        return f"SignResult(sign={self.sign:+d}, method={self.method})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SignResult):
            return NotImplemented
        return self.sign == other.sign and self.method == other.method

    def __hash__(self) -> int:
        return hash((self.sign, self.method))


# ------------------------------------------------------------
# Certified sqrt bracket (integer kernel)
# ------------------------------------------------------------


def certified_bracket_int(s: int, bits: int) -> Tuple[int, int]:
    """
    Exact integer bracket for sqrt(s), scaled by 2**bits.

    Returns (lo_scaled, hi_scaled) with

        lo_scaled / 2**bits  <=  sqrt(s)  <=  hi_scaled / 2**bits
        hi_scaled - lo_scaled == 1        (width exactly 2**-bits)

    Proof: t = isqrt(s << 2*bits) gives t**2 <= s*2**(2*bits) < (t+1)**2.
    Both t and t+1 are >= 0, so dividing by 2**(2*bits) and taking square
    roots gives the claim.
    """
    if s < 0:
        raise ValueError(f"sqrt bracket needs s >= 0, got {s}")
    if bits < 0:
        raise ValueError(f"bits must be >= 0, got {bits}")
    if s == 0:
        return (0, 0)
    t = isqrt(s << (2 * bits))
    return (t, t + 1)


def _sum_bracket_int(agg: Dict[int, int], bits: int) -> Tuple[int, int]:
    """
    Certified bracket for S = SUM_s agg[s] * sqrt(s), scaled by 2**bits.

    Returns (lo_scaled, hi_scaled) integers with S in
    [lo_scaled/2**bits, hi_scaled/2**bits].

    Endpoint selection is sign-aware: a positive coefficient multiplies the
    [lo, hi] bracket of sqrt(s) into [c*lo, c*hi], a negative one into
    [c*hi, c*lo].  Summing lower endpoints and upper endpoints therefore
    keeps the enclosure valid term by term -- no cancellation is assumed and
    none is lost silently.
    """
    lo = 0
    hi = 0
    for s, c in agg.items():
        if c == 0:
            continue
        if s == 1:
            # sqrt(1) = 1 exactly; no rounding, no width.
            exact = c << bits
            lo += exact
            hi += exact
            continue
        t, t1 = certified_bracket_int(s, bits)
        if c > 0:
            lo += c * t
            hi += c * t1
        else:
            lo += c * t1
            hi += c * t
    return (lo, hi)


def certified_bracket(agg: Dict[int, int], bits: int) -> Tuple["object", "object"]:
    """
    Certified rational bracket [lo, hi] for S = SUM_s agg[s] * sqrt(s).

    `agg` maps squarefree part -> integer coefficient.  Returns exact
    Fractions; the enclosure is proved by the integer kernel above.
    """
    from fractions import Fraction as F

    lo_scaled, hi_scaled = _sum_bracket_int(agg, bits)
    denom = 1 << bits
    return (F(lo_scaled, denom), F(hi_scaled, denom))


# ------------------------------------------------------------
# Bracketing driver
# ------------------------------------------------------------


def _start_bits(agg: Dict[int, int], max_input_bits: int) -> int:
    """
    First precision to try.

    The width of the bracket is C * 2**-bits with C = sum |agg[s]|, so the
    useful scale is log2(C) plus the scale of the radicands.  Starting near
    that avoids the wasted halvings of a naive start, without assuming
    anything about how close S is to 0 -- the doubling loop below still has
    to carry the real work.
    """
    c_sum = sum(abs(c) for c in agg.values())
    c_bits = c_sum.bit_length()
    # 2*max_input_bits: the naive worst case for a 2-term sum is that S sits
    # one bit above the top of the radicands, needing ~2x that in fractions.
    return 16 + c_bits + 2 * max_input_bits


def _aggregate(terms: Iterable[Tuple[int, int]]):
    """Validate and aggregate terms.  Returns (agg, max_input_bits)."""
    materialised = list(terms)
    for coeff, rad in materialised:
        if not isinstance(coeff, int) or not isinstance(rad, int):
            raise TypeError(f"terms must be (int, int), got ({coeff!r}, {rad!r})")
        if rad <= 0:
            raise ValueError(f"radicand must be positive, got {rad}")
    return decompose_and_aggregate(materialised)


def _sign_of_aggregated(
    agg: Dict[int, int],
    max_input_bits: int,
    max_bits: int,
    start_bits: Optional[int],
) -> SignResult:
    """
    Sign decision on an aggregated sum, where every key is a squarefree part.

    The zero test is the first thing that happens and it is exact, because
    {sqrt(d) : d distinct squarefree} is linearly independent over Q.
    """
    agg = {s: c for s, c in agg.items() if c != 0}

    if not agg:
        return SignResult(0, "zero-test", {"aggregated": {}})

    # A sum with no irrational part is a plain integer comparison.
    if len(agg) == 1 and 1 in agg:
        total = agg[1]
        sign = (total > 0) - (total < 0)
        return SignResult(sign, "rational", {"sum": str(total)})

    bits = start_bits if start_bits is not None else _start_bits(agg, max_input_bits)
    c_sum = sum(abs(c) for c in agg.values())

    while True:
        if bits > max_bits:
            # Width of the bracket we would have tried next, as a log2 exponent.
            raise SignUnknownError(
                max_bits,
                c_sum.bit_length() - bits,
                "raise max_bits to continue",
            )
        lo, hi = _sum_bracket_int(agg, bits)
        witness = {"bits": bits, "terms": len(agg)}
        if lo > 0:
            witness["bracket_num"] = (lo, 1 << bits)
            return SignResult(1, "certified-interval", witness)
        if hi < 0:
            witness["bracket_num"] = (hi, 1 << bits)
            return SignResult(-1, "certified-interval", witness)
        bits *= 2


def exact_sign(
    terms: Iterable[Tuple[int, int]],
    max_bits: Optional[int] = DEFAULT_MAX_BITS,
    start_bits: Optional[int] = None,
) -> SignResult:
    """
    Exact sign of S = SUM_i c_i * sqrt(a_i).

    `terms` is an iterable of (coeff, radicand) integer pairs.

    Returns a SignResult whose sign is -1, 0 or +1.  Raises SignUnknownError
    if the sign is not decided within `max_bits` fractional bits; it never
    returns an undecided guess.

    `max_bits` bounds work on adversarial inputs (where S is exponentially
    close to 0).  `start_bits` overrides the initial precision, which is
    only useful for tests.

    `max_bits=None` means `DEFAULT_MAX_BITS`, matching
    `intervalsum.certified_interval_for_S`, so a caller holding an optional
    cap can forward it to both deciders without translating it first.  This
    function used to reject `None` outright with a TypeError from the
    comparison below, which made the two deciders disagree about what an
    absent cap meant.
    """
    if max_bits is None:
        max_bits = DEFAULT_MAX_BITS
    if max_bits <= 0:
        raise ValueError(f"max_bits must be positive, got {max_bits}")

    if start_bits is not None and start_bits <= 0:
        # Checked here rather than letting it surface as a ValueError from
        # inside `certified_bracket_int`, where the message says "bits" and
        # names no parameter.
        raise ValueError(f"start_bits must be positive, got {start_bits}")

    agg, max_input_bits = _aggregate(terms)
    witness_agg = {str(s): c for s, c in agg.items()}

    if all(c == 0 for c in agg.values()):
        return SignResult(0, "zero-test", {"aggregated": witness_agg})

    result = _sign_of_aggregated(agg, max_input_bits, max_bits, start_bits)
    result.witness.setdefault("aggregated", witness_agg)
    return result


def sign_of_sum(terms: Iterable[Tuple[int, int]], **kwargs) -> SignResult:
    """Public alias for exact_sign."""
    return exact_sign(terms, **kwargs)


def exact_sign_from_aggregated(agg: Dict[int, int], **kwargs) -> SignResult:
    """
    Sign of an already-aggregated sum (squarefree part -> integer coeff).

    Skips decomposition and works directly on `agg`, so no radicand is
    squared back into a term -- a round trip through (c, s*s) would hand the
    decomposition routine a perfect square and silently turn sqrt(s) into s.
    The zero test still runs; it just runs on agg.
    """
    for s in agg:
        if not isinstance(s, int) or s < 1:
            raise ValueError(f"aggregated keys must be squarefree parts >= 1, got {s!r}")
    for c in agg.values():
        if not isinstance(c, int):
            raise TypeError(f"aggregated coefficients must be int, got {c!r}")

    max_input_bits = max((s.bit_length() for s in agg), default=0)
    return _sign_of_aggregated(
        dict(agg), max_input_bits, kwargs.pop("max_bits", DEFAULT_MAX_BITS),
        kwargs.pop("start_bits", None),
    )


# ------------------------------------------------------------
# Self-test
# ------------------------------------------------------------

if __name__ == "__main__":
    cases = [
        ([(1, 12), (-2, 3)], 0, "sqrt(12) - 2sqrt(3) = 0"),
        ([(1, 12), (-1, 3)], 1, "sqrt(12) - sqrt(3) = sqrt(3) > 0"),
        ([(1, 2), (-1, 3)], -1, "sqrt(2) - sqrt(3) < 0"),
        ([(1, 2), (1, 3), (-1, 5)], 1, "sqrt(2)+sqrt(3)-sqrt(5) > 0"),
        ([(1, 4), (-1, 9)], -1, "2 - 3 = -1 (rational)"),
        ([(3, 2), (-4, 3)], -1, "3sqrt(2) - 4sqrt(3) < 0"),
        ([(1, 1)], 1, "1 > 0"),
    ]
    bad = 0
    for terms, want, why in cases:
        got = exact_sign(terms)
        ok = got.sign == want
        bad += not ok
        print(f"{'ok ' if ok else 'FAIL'} {why:<34} -> {got.sign:+d} ({got.method})")
    print("all passed" if not bad else f"{bad} FAILED")
    raise SystemExit(1 if bad else 0)
