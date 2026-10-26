"""
Instance families for the SRS probe.

Each generator returns a dict with

    terms    [(coeff, radicand), ...]   the sum S = SUM coeff*sqrt(radicand)
    N_bits   input bit size
    k        number of terms
    family   the family tag
    metadata family-specific extras

N_bits is computed with int.bit_length(), never log2(): log2 on a float
loses the top bits of a 200-bit radicand and raises OverflowError past
2**1024, which is exactly the regime the F1 probe exists to explore.

Families:
- F1: Pell construction, sign and separation both proved in f1_theorem
- F2: random radicands
- F3: greedy near-cancel search (searches on a float surrogate, so it is a
       heuristic for finding hard instances, never a correctness claim)
- F4: product composition, (sqrt(a) - sqrt(b)) * (sqrt(c) - sqrt(d))
- F5: squaring, (sqrt(a) + sqrt(b))^2 collapsed to two radicals
- F6: rank-3 pigeonhole gap, the sorted-gap construction behind Thm 2.5 at
       r=3.  This is the family that attains the conjectured sharp rate
       E/N = (r-1)/r = 2/3, so it is the only family here that probes the
       CEILING rather than looking for something worse.
"""

import random
from fractions import Fraction
from typing import Dict, List, Tuple

from squarefree import is_squarefree as _is_squarefree_exact
from squarefree import random_squarefree

__all__ = [
    "SRSInstance",
    "generate_f1",
    "generate_f1_sequence",
    "generate_f2",
    "generate_f3",
    "generate_f4",
    "generate_f5",
    "generate_f6",
    "generate_f6_sequence",
    "generate_family",
    "generate_family_sequence",
    "FAMILY_GENERATORS",
    "FAMILY_SIGNATURES",
    "input_bits",
    "evaluate_float_surrogate",
    "is_squarefree",
]

DEFAULT_SEED = 42

# Largest radicand width any generator will produce.  This is a *cost* bound,
# not a soundness one: `squarefree.random_squarefree` builds radicands as
# products of distinct primes each below 2^64, where Miller-Rabin is proven
# correct, so a generated radicand is a certified prime product at every width
# and its squarefreeness is a theorem rather than a probable property.  Past
# this size the sampler simply needs more certified factors and gets slower.
MAX_GENERATED_BITS = 512


class SRSInstance(dict):
    """
    An SRS instance.

    A dict subclass so the existing `inst['terms']` access keeps working,
    with the same fields available as attributes.
    """

    __slots__ = ()

    @property
    def terms(self) -> List[Tuple[int, int]]:
        return self["terms"]

    @property
    def N_bits(self) -> int:
        return self["N_bits"]

    @property
    def k(self) -> int:
        return self["k"]

    @property
    def family(self) -> str:
        return self["family"]

    @property
    def metadata(self) -> dict:
        return self["metadata"]


def input_bits(terms: List[Tuple[int, int]]) -> int:
    """
    Input size N = sum over terms of the bit length of |coeff| and radicand.

    Exact: bit_length() on Python ints, so no rounding and no overflow.
    """
    total = 0
    for coeff, rad in terms:
        total += abs(coeff).bit_length() + rad.bit_length()
    return total


def is_squarefree(n: int) -> bool:
    """
    Exact squarefree test.

    Delegates to squarefree.is_squarefree, which factors whatever trial
    division leaves behind.  The local version this replaces looped to
    isqrt(n) and -- worse -- returned True unconditionally for n >= 10**12,
    so F2 would happily emit radicands with square factors and quietly
    change the instance it claimed to be generating.
    """
    return _is_squarefree_exact(n)


def _validate_terms(terms: List[Tuple[int, int]]) -> None:
    for coeff, rad in terms:
        if not isinstance(coeff, int) or not isinstance(rad, int):
            raise TypeError(f"terms must be (int, int), got ({coeff!r}, {rad!r})")
        if rad <= 0:
            raise ValueError(f"radicand must be positive, got {rad}")


def _make(terms, family: str, metadata: dict) -> SRSInstance:
    _validate_terms(terms)
    return SRSInstance(
        terms=terms,
        N_bits=input_bits(terms),
        k=len(terms),
        family=family,
        metadata=metadata,
    )


# ------------------------------------------------------------
# F1: Pell construction (sign and separation both proved)
# ------------------------------------------------------------

def generate_f1(a: int, n: int) -> SRSInstance:
    """
    F1: S_n = sqrt(p_n^2) - sqrt(alpha q_n^2) = p_n - q_n sqrt(alpha),
    alpha = a(a+1), p_n/q_n the CF convergents of sqrt(alpha).

    f1_theorem proves p_n^2 - alpha q_n^2 = +1 on odd n and -a on even n, so
    the sign is +1 / -1 respectively and needs no computation to know.
    """
    if not isinstance(a, int) or isinstance(a, bool):
        raise TypeError(f"a must be an int, got {a!r}")
    if a < 2:
        raise ValueError(f"a must be >= 2, got {a}")
    if not isinstance(n, int) or n < 0:
        raise ValueError(f"n must be >= 0, got {n!r}")

    from f1_theorem import expected_sign, exact_pell_identity

    p, q, residue = exact_pell_identity(a, n)
    alpha = a * (a + 1)
    return _make(
        [(1, p * p), (-1, alpha * q * q)],
        "F1",
        {
            "a": a, "n": n, "p": p, "q": q, "alpha": alpha,
            "residue": residue,
            "known_sign": expected_sign(a, n),
        },
    )


def generate_f1_sequence(a: int, n_start: int, n_end: int) -> List[SRSInstance]:
    if n_end < n_start:
        raise ValueError(f"n_end={n_end} < n_start={n_start}")
    return [generate_f1(a, n) for n in range(n_start, n_end + 1)]


# ------------------------------------------------------------
# F2: random radicands
# ------------------------------------------------------------

def generate_f2(k: int, n_bits: int, seed: int = None) -> SRSInstance:
    """
    F2: k terms with random squarefree radicands of about n_bits bits.

    Radicands are verified squarefree exactly, so the family is what it
    claims.  Rejection sampling can be slow for large n_bits, so the bound
    is honest: it tries, and if it cannot fill k parts in max_tries it
    raises rather than returning a short instance.
    """
    if k < 1:
        raise ValueError(f"k must be >= 1, got {k}")
    if n_bits < 2:
        raise ValueError(f"n_bits must be >= 2, got {n_bits}")
    if n_bits > MAX_GENERATED_BITS:
        raise ValueError(
            f"n_bits={n_bits} exceeds MAX_GENERATED_BITS={MAX_GENERATED_BITS}. "
            f"Radicands are built as a product of distinct primes each below "
            f"2^64, so squarefreeness stays certified above this size; the "
            f"bound is on sampler cost, which grows as the number of certified "
            f"factors. Raise MAX_GENERATED_BITS if that trade is wanted."
        )

    rng = random.Random(DEFAULT_SEED if seed is None else seed)

    parts: List[int] = []
    seen = set()
    max_tries = max(1000, 20 * k)
    tries = 0
    while len(parts) < k and tries < max_tries:
        tries += 1
        n = random_squarefree(n_bits, rng)
        if n in seen:
            continue
        seen.add(n)
        parts.append(n)
    if len(parts) < k:
        # Distinguish "the window is too small to hold k distinct values"
        # from "the sampler got unlucky".  At 2 bits the window [2,4) holds
        # exactly one prime, so k=3 there is unsatisfiable rather than
        # unlucky, and saying so saves a reader from hunting for a
        # seed-dependent explanation that does not exist.
        window = 1 << n_bits
        raise RuntimeError(
            f"generate_f2: produced only {len(parts)} distinct radicands in "
            f"{max_tries} tries. The {n_bits}-bit window [2^{n_bits - 1}, "
            f"2^{n_bits}) has {window // 2} members, so k={k} distinct "
            f"squarefree radicands may not fit; raise n_bits or lower k."
        )

    coeffs = [rng.choice([-1, 1]) for _ in parts]
    terms = [(coeffs[i], parts[i]) for i in range(k)]
    return _make(terms, "F2", {"k": k, "n_bits": n_bits, "seed": seed,
                               "parts": parts, "coeffs": coeffs})


# ------------------------------------------------------------
# F3: greedy near-cancel search
# ------------------------------------------------------------

def evaluate_float_surrogate(terms: List[Tuple[int, int]]) -> float:
    """
    Float approximation of S.  MEASUREMENT ONLY.

    Named for what it is.  The previous `evaluate_exact` returned a float
    while claiming exactness, and its callers (F3's search, the tests) then
    treated the result as if it were the value of S.  Nothing in the
    decision path may call this; use signexact.exact_sign instead.
    """
    from math import sqrt
    return sum(coeff * sqrt(rad) for coeff, rad in terms)


def generate_f3(k: int, n_bits: int, target_bits: int = 20,
                max_iter: int = 10000, seed: int = None) -> SRSInstance:
    """
    F3: hill-climb toward a small |S| by mutating one term at a time.

    The objective is the float surrogate, so this is a *search heuristic*:
    it finds instances that look near-cancelling to ~16 digits and no more.
    Whether the winner really is near-cancelling is decided exactly
    afterwards by signexact, and the reported best_value is labelled a
    surrogate.  On large instances the surrogate is in fact useless -- its
    error exceeds |S| -- which is the whole reason the exact decider exists.
    """
    inst = generate_f2(k, n_bits, seed)
    rng = random.Random(DEFAULT_SEED if seed is None else seed)

    best_terms = list(inst["terms"])
    best_val = abs(evaluate_float_surrogate(best_terms))

    rad_hi = (1 << min(n_bits, 20)) - 1
    for _ in range(max_iter):
        idx = rng.randrange(len(best_terms))
        new_coeff = rng.choice([-1, 1])
        new_rad = rng.randint(2, rad_hi)
        if not is_squarefree(new_rad):
            continue

        test_terms = list(best_terms)
        test_terms[idx] = (new_coeff, new_rad)

        val = abs(evaluate_float_surrogate(test_terms))
        if val < best_val:
            best_val = val
            best_terms = test_terms
            if best_val < 2.0 ** (-target_bits):
                break

    from signexact import exact_sign

    decided = exact_sign(best_terms)
    return _make(
        best_terms,
        "F3",
        {
            "k": k,
            "n_bits": n_bits,
            "seed": seed,
            "parts": [r for _, r in best_terms],
            "surrogate_abs_value": best_val,
            "surrogate_is_approximate": True,
            "target_bits": target_bits,
            "exact_sign": decided.sign,
            "exact_method": decided.method,
            "exact_bits": decided.witness.get("bits"),
        },
    )


# ------------------------------------------------------------
# F4: product composition
# ------------------------------------------------------------

def generate_f4(a: int, b: int, c: int, d: int) -> SRSInstance:
    """
    F4: (sqrt(a) - sqrt(b)) * (sqrt(c) - sqrt(d))
       = sqrt(ac) - sqrt(ad) - sqrt(bc) + sqrt(bd).

    Terms are combined by exact radicand equality, so a factorisation that
    collapses (a=2, b=8, c=3, d=12 -> sqrt(24) twice) shows up as one
    coefficient of -2.
    """
    for name, v in (("a", a), ("b", b), ("c", c), ("d", d)):
        if not isinstance(v, int) or isinstance(v, bool):
            raise TypeError(f"{name} must be an int, got {v!r}")
        if v < 1:
            raise ValueError(f"{name} must be >= 1, got {v}")

    combined: Dict[int, int] = {}
    for coeff, rad in ((1, a * c), (-1, a * d), (-1, b * c), (1, b * d)):
        combined[rad] = combined.get(rad, 0) + coeff
    terms = [(v, k) for k, v in sorted(combined.items()) if v != 0]

    if not terms:
        raise ValueError("F4 instance collapsed to the empty sum")
    return _make(terms, "F4", {"a": a, "b": b, "c": c, "d": d})


# ------------------------------------------------------------
# F5: squaring
# ------------------------------------------------------------

def generate_f5(a: int, b: int, s: int = -1) -> SRSInstance:
    """
    F5: squaring, as a two-term root sum.

        (sqrt(a) + sqrt(b))^2 = (a+b) + 2 sqrt(ab) = (a+b) + sqrt(4ab)
        (sqrt(a) - sqrt(b))^2 = (a+b) - 2 sqrt(ab) = (a+b) - sqrt(4ab)

    so both are instances of (a+b) + s * sqrt(4ab) for s = +1 or s = -1, and
    `s` selects which.  The emitted terms are [(1, (a+b)^2), (s, 4ab)].

    An earlier version documented this `s` parameter but did not have it: the
    function took only (a, b) and always emitted s = -1, so a caller following
    the docstring would look for a way to pick the other sign and find none.
    The default is -1, so existing callers are unaffected.

    WORTH KNOWING BEFORE USING THIS FAMILY FOR SEPARATION
    ----------------------------------------------------
    Both forms are squares of real numbers, so S >= 0 always, and S = 0
    exactly when a == b with s = -1.  F5 therefore cannot produce a negative
    instance, and its non-degenerate instances are bounded away from zero by
    the gap between distinct square roots.  That makes it a consistency check
    on the sign decider rather than a source of hard instances; the negative
    sign in SRS comes from cancellation between separate families of terms,
    not from a single square.
    """
    for name, v in (("a", a), ("b", b)):
        if not isinstance(v, int) or isinstance(v, bool):
            raise TypeError(f"{name} must be an int, got {v!r}")
        if v < 1:
            raise ValueError(f"{name} must be >= 1, got {v}")
    if not isinstance(s, int) or isinstance(s, bool):
        raise TypeError(f"s must be an int, got {s!r}")
    if s not in (1, -1):
        raise ValueError(f"s must be +1 or -1, got {s}")

    terms = [(1, (a + b) * (a + b)), (s, 4 * a * b)]
    return _make(terms, "F5", {"a": a, "b": b, "s": s})


# ------------------------------------------------------------
# Registry
# ------------------------------------------------------------

# ------------------------------------------------------------
# F6: rank-3 pigeonhole gap (the Thm 2.5 ceiling at r = 3)
# ------------------------------------------------------------

def generate_f6(B: int, radicands: Tuple[int, ...] = (2, 3, 5)) -> SRSInstance:
    """
    F6: the sorted-gap construction of Thm 2.5 at r = 3.

    Consider all (B+1)^3 sums  v(C) = C_1 sqrt(d_1) + C_2 sqrt(d_2) +
    C_3 sqrt(d_3)  with C_i in {0..B}.  They are pairwise distinct (Cor 1.3)
    and lie in an interval of length B*(sum sqrt(d_i)), so two consecutive
    ones differ by at most

        B * sum_i sqrt(d_i) / ((B+1)^3 - 1)  ~  sum_i sqrt(d_i) / B^2,

    and their difference is a nonzero sum with |C_i| <= B.  That difference is
    the instance returned here.

    WHY THIS FAMILY IS DIFFERENT FROM F1-F5.  Every other family searches for
    instances that are *worse* than expected.  This one provably attains the
    rate Thm 2.5 shows is the ceiling, E/N -> (r-1)/r = 2/3, and Thm 2.9
    shows that rate is also the *guaranteed* one.  So F6 is the family that
    tests whether the ceiling is where it is claimed to be, and it is the
    only family in this repo capable of refuting the conjecture by exhibiting
    E/N above 2/3 -- though as a corollary of Thm 2.5 that would require
    breaking Thm 2.5, not merely finding a lucky instance.

    The gap is located by sorting on a float surrogate, exactly as the
    existing F3 documents, so the SEARCH is a heuristic.  The RESULT is not:
    the returned C is re-verified exactly by the caller (see
    tests/test_families.py and tests/test_proofs.py), and Thm 2.5 bounds what
    any such C can achieve.
    """
    B = int(B)
    if B < 1:
        raise ValueError(f"F6 needs B >= 1, got {B}")
    d = tuple(int(x) for x in radicands)
    if len(d) != 3:
        raise ValueError(f"F6 is the rank-3 construction, got {len(d)} radicands")
    for x in d:
        if not _is_squarefree_exact(x) or x <= 0:
            raise ValueError(f"F6 radicands must be positive squarefree, got {x}")

    roots = [float(x) ** 0.5 for x in d]

    # Pass 1: sort the (B+1)^3 values and locate the smallest GAP between
    # consecutive ones.  Thm 2.5 bounds this by  B*sum_i sqrt(d_i) / ((B+1)^3-1).
    # Floats are used only to RANK candidates, exactly as F3 documents.
    grid = [
        (c1 * roots[0] + c2 * roots[1] + c3 * roots[2], (c1, c2, c3))
        for c1 in range(B + 1)
        for c2 in range(B + 1)
        for c3 in range(B + 1)
    ]
    grid.sort(key=lambda t: t[0])
    best_gap = None
    best_diff = None
    for i in range(len(grid) - 1):
        g = grid[i + 1][0] - grid[i][0]
        if best_gap is None or g < best_gap:
            best_gap = g
            best_diff = (grid[i][1], grid[i + 1][1])

    # Pass 2: the gap is a DIFFERENCE of two grid points, i.e. coefficients
    # in {-(B+1)..(B+1)}, not {0..B}.  Certify that difference exactly, so
    # the instance returned is verified rather than merely searched.
    (a1, a2, a3), (b1, b2, b3) = best_diff
    terms = [
        (a1 - b1, d[0]),
        (a2 - b2, d[1]),
        (a3 - b3, d[2]),
    ]
    terms = [(c, r) for c, r in terms if c]      # drop exact zeros
    from signexact import exact_sign
    res = exact_sign(terms)
    if res.sign == 0:
        raise ValueError("F6 gap coefficients cancelled exactly; cannot certify")
    num, den = res.witness["bracket_num"]
    mag = Fraction(abs(num), den)

    inst = _make(terms, "F6", {
        "B": B,
        "radicands": d,
        "exact_sign": res.sign,
        "exact_mag": mag,
        # PROOFS.md Thm 2.5's guarantee for |C_i| <= B at r = 3:
        #     |S| <= (sum_i sqrt(d_i)) / B^(r-1)  =  (sum sqrt(d)) / B^2.
        # Reported as a float for display; the exact |S| above is a Fraction
        # and the test compares them exactly.
        "thm25_bound": sum(r ** 0.5 for r in d) / B ** (len(d) - 1),
        "grid_gap": best_gap,        # float ranking of the gap, display only
        "note": "Thm 2.5 ceiling construction at r=3; E/N -> 2/3",
    })
    return inst


def generate_f6_sequence(Bs, radicands: Tuple[int, ...] = (2, 3, 5)):
    """F6 at increasing B, so the E/N -> 2/3 trend is visible."""
    return [generate_f6(B, radicands) for B in Bs]


FAMILY_GENERATORS = {
    "F1": generate_f1,
    "F2": generate_f2,
    "F3": generate_f3,
    "F4": generate_f4,
    "F5": generate_f5,
    "F6": generate_f6,
}

# Required keyword arguments per family, for iteration helpers.
FAMILY_SIGNATURES = {
    "F1": ("a", "n"),
    "F2": ("k", "n_bits"),
    "F3": ("k", "n_bits"),
    "F4": ("a", "b", "c", "d"),
    "F5": ("a", "b", "s"),
    "F6": ("B",),
}


def generate_family(family: str, **kwargs) -> SRSInstance:
    if family not in FAMILY_GENERATORS:
        raise ValueError(
            f"unknown family: {family!r}; expected one of {sorted(FAMILY_GENERATORS)}"
        )
    return FAMILY_GENERATORS[family](**kwargs)


def generate_family_sequence(family: str, param_ranges: Dict) -> List[SRSInstance]:
    """
    Cartesian product of per-parameter ranges.

    `param_ranges` maps parameter name -> iterable of values, e.g.
    {"a": [2, 3], "n": [1, 3, 5]}.  Parameters absent from the map are not
    passed.  Previously this ignored `param_ranges` entirely and returned a
    single instance, so any sweep built on it silently sampled once.
    """
    import itertools

    if family not in FAMILY_SIGNATURES:
        raise ValueError(
            f"unknown family: {family!r}; expected one of {sorted(FAMILY_SIGNATURES)}"
        )
    unknown = set(param_ranges) - set(FAMILY_SIGNATURES[family])
    if unknown:
        raise ValueError(
            f"family {family} has no parameters {sorted(unknown)}; "
            f"expected a subset of {list(FAMILY_SIGNATURES[family])}"
        )

    names = list(param_ranges)
    combos = itertools.product(*(list(param_ranges[n]) for n in names))
    out = []
    for combo in combos:
        out.append(generate_family(family, **dict(zip(names, combo))))
    return out


# ------------------------------------------------------------
# Self-test
# ------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 64)
    print("SRS instance families")
    print("=" * 64)

    print("\nF1: Pell (sign known in closed form from the residue)")
    for a in (2, 3, 4):
        for n in (1, 2, 3):
            inst = generate_f1(a, n)
            print(f"  a={a} n={n}: k={inst['k']} N={inst['N_bits']:>4} "
                  f"residue={inst['metadata']['residue']:>3} "
                  f"sign={inst['metadata']['known_sign']:+d}")

    print("\nF2: random squarefree")
    inst = generate_f2(3, 20, seed=42)
    print(f"  terms={inst['terms']} N={inst['N_bits']}")
    assert all(is_squarefree(r) for _, r in inst["terms"])

    print("\nF3: greedy near-cancel (surrogate objective)")
    inst = generate_f3(3, 20, target_bits=10, max_iter=200, seed=42)
    print(f"  terms={inst['terms']}")
    print(f"  surrogate |S| ~ {inst['metadata']['surrogate_abs_value']:.3e} "
          f"(approximate), exact sign {inst['metadata']['exact_sign']:+d}")

    print("\nF4: product")
    inst = generate_f4(2, 3, 5, 7)
    print(f"  terms={inst['terms']} N={inst['N_bits']}")
    inst = generate_f4(2, 8, 3, 12)
    print(f"  (2,8,3,12) -> {inst['terms']}   (sqrt(24) combines to -2)")

    print("\nF5: squaring")
    inst = generate_f5(3, 2)
    print(f"  terms={inst['terms']} N={inst['N_bits']}")

    print("\nF1 sequence and family sweep:")
    print(f"  generate_f1_sequence(3, 1, 3) -> {len(generate_f1_sequence(3, 1, 3))} instances")
    swept = generate_family_sequence("F1", {"a": [2, 3], "n": [1, 3]})
    print(f"  generate_family_sequence('F1', a=[2,3], n=[1,3]) -> {len(swept)} instances")

    print("\nall family generators ok")
