"""
Exact squarefree decomposition and aggregation for the Square-Root Sum Problem.

Zero test: S = 0 iff after squarefree decomposition and coefficient aggregation,
all net coefficients vanish. This is O(N) for N input integers, no approximation.
"""

from math import isqrt
from collections import Counter
from typing import List, Tuple, Dict


# ------------------------------------------------------------
# Squarefree decomposition (exact, trial division + isqrt)

_SMALL_PRIMES: List[int] = []

# Upper bound of the trial-division sieve.  Any prime factor of a radicand at
# or below this is resolved exactly and cheaply; anything larger is left to
# `isqrt`, the primality test, and Pollard's rho.  The bound is named here and
# used in `_build_small_primes` so a reader does not have to go looking for a
# second copy of the number.
_PRIME_SIEVE_BOUND = 10**6

def _build_small_primes():
    """Sieve up to 10^6 — once, on first use. Bound is documented for honesty.

    Exact throughout: the inner bound uses isqrt, not bound**0.5, because
    10**6**0.5 is a float and this module is on the decision path.
    """
    global _SMALL_PRIMES
    if _SMALL_PRIMES:
        return
    bound = _PRIME_SIEVE_BOUND
    sieve = bytearray(b'\x01') * (bound + 1)
    sieve[0:2] = b'\x00\x00'
    root = isqrt(bound)
    for p in range(2, root + 1):
        if sieve[p]:
            step = p
            start = p * p
            sieve[start:bound+1:step] = b'\x00' * ((bound - start) // step + 1)
    _SMALL_PRIMES = [i for i, is_p in enumerate(sieve) if is_p]


def _trial_division(n: int) -> Tuple[int, int, int]:
    """
    Trial-divide n by sieved primes up to `_PRIME_SIEVE_BOUND`.
    Returns (remaining, square_part, squarefree_part_so_far).

    The bound is not a parameter.  It used to be, with `bound: int = None`,
    but no caller ever passed it, so every code path and every docstring that
    mentioned a caller-chosen bound was describing a capability that did not
    exist.
    """
    _build_small_primes()
    bound = _PRIME_SIEVE_BOUND

    square_part = 1
    squarefree_part = 1
    remaining = n

    for p in _SMALL_PRIMES:
        if p > bound:
            break
        if p * p > remaining:
            break
        exp = 0
        while remaining % p == 0:
            remaining //= p
            exp += 1
        if exp > 0:
            square_part *= p ** (exp // 2)
            if exp % 2 == 1:
                squarefree_part *= p

    return remaining, square_part, squarefree_part


def _is_square(n: int) -> bool:
    """Exact perfect-square test using isqrt."""
    r = isqrt(n)
    return r * r == n


def _is_probable_prime(n: int) -> bool:
    """
    Deterministic Miller–Rabin for n < 2^64.
    For n >= 2^64, uses probabilistic bases (labelled probabilistic).
    """
    if n < 2:
        return False
    # Small prime checks
    small_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small_primes:
        if n % p == 0:
            return n == p

    # Deterministic bases for n < 2^64
    if n < 2**64:
        bases = (2, 325, 9375, 28178, 450775, 9780504, 1795265022)
    else:
        # Probabilistic for larger n — callers must handle labelling
        bases = (2, 3, 5, 7, 11, 13)

    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1

    for a in bases:
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _pollard_rho(n: int) -> int:
    """
    Find a non-trivial factor of composite n, by Brent's variant of Pollard's rho.

    Deterministic: the sequence is seeded from n itself and `c` is stepped
    deterministically, so the same n always yields the same factor.  There is
    no randomness and hence nothing to reproduce.

    Only ever called on n already known composite (Miller-Rabin said so), so
    a returned factor is always a genuine divisor.

    Raises RhoFailure if no factor is found within the effort budget.

    BUDGET, AND WHAT IT ACTUALLY BUYS
    ----------------------------------
    Pollard's rho needs on the order of the square root of the *smallest*
    prime factor to find it, so the worst case for a balanced semiprime is
    about n^(1/4): a 122-bit product of two 61-bit primes needs ~2^30
    iterations, not the ~2^17 this budget allows.

    The formula is therefore a deliberate *cheap* cap, not an estimate of
    n^(1/4):

        budget = max(2^16, (bit_length + 3) // 4 * 2^12)

    The comment here used to read "~n**(1/4) expected, with a 2**12 safety
    factor", which is wrong twice over: it treats bit_length/4 as though it
    were n^(1/4) (a linear quantity for an exponential one), and the 2^12
    factor cannot cover a gap that grows with the input.  Measured against the
    2^(bit_length/4) it should be scaling towards, the budget is 2^16 at 64
    bits -- correct there, by coincidence of the floor -- and about 2^17 at
    200 bits where 2^50 would be needed.

    That shortfall is safe, not silent: the outcome is a raised RhoFailure
    naming the radicand, never a wrong factorisation.  What it does mean is
    that a balanced semiprime above roughly 70 bits is out of reach for this
    implementation, which is why `squarefree_decompose` refuses on it and why
    `random_squarefree` never builds one.  Past that, ECM or a pre-factored
    input is the way forward.
    """
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3

    # Cheap cap, deliberately short of the n^(1/4) worst case; see above.
    budget = max(1 << 16, (n.bit_length() + 3) // 4 * (1 << 12))

    # Deterministic seed offsets, tried in turn.  A single unlucky (y, c)
    # can cycle with period dividing n without finding anything; moving to
    # the next seed is cheaper than restarting from scratch.
    for attempt in range(_RHO_RESTARTS):
        y, c, m = 1 + attempt, 1 + attempt, 1
        g = q = r = 1
        cycles = 0
        while g == 1:
            if cycles > budget:
                break
            x = y
            for _ in range(r):
                y = (y * y + c) % n
                cycles += 1
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = _gcd(q, n)
                k += m
            r *= 2
        if g == n:
            # The batch accumulated n as a whole: back up over the last block,
            # testing one difference at a time, which recovers a proper factor.
            g = 1
            while g == 1:
                ys = (ys * ys + c) % n
                g = _gcd(abs(x - ys), n)
        if 1 < g < n:
            return g
        # else: this seed failed, try the next one

    raise RhoFailure(
        f"Pollard's rho found no factor of a {n.bit_length()}-bit cofactor "
        f"within {budget} iterations x {_RHO_RESTARTS} restarts. Splitting "
        f"numbers this large needs a stronger algorithm (Pollard p-1, ECM, "
        f"SymPy, or a pre-factored input)."
    )


class RhoFailure(RuntimeError):
    """Pollard's rho could not split a cofactor within its effort budget."""


def _gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return a


# How many (y, c) seeds to try before giving up. Each is a fresh deterministic
# start; rho fails outright only when every seed degenerates.
_RHO_RESTARTS = 6


def _factor(n: int, out: Dict[int, int], _depth: int = 0) -> None:
    """
    Accumulate the full prime factorisation of n into `out` as p -> exponent.

    Exact: Miller-Rabin certifies primality for n < 2^64 and is
    probabilistically correct above that (the caller is told so), and rho
    returns a verified divisor, so the product of the collected factors
    always equals n exactly regardless of which primality branch was taken.

    Raises RhoFailure if a cofactor is too large to split.  That is a real
    limit of this implementation and is surfaced as such rather than
    returning a wrong factorisation.
    """
    if n == 1:
        return
    if _depth > _MAX_FACTOR_DEPTH:
        raise RhoFailure(
            f"factorisation recursion exceeded {_MAX_FACTOR_DEPTH} levels; "
            f"n has more prime factors than that, or rho is cycling"
        )
    if _is_probable_prime(n):
        out[n] = out.get(n, 0) + 1
        return
    d = _pollard_rho(n)
    if d == n or d == 1:
        raise RhoFailure(f"Pollard's rho returned the useless factor {d} for n={n}")
    _factor(d, out, _depth + 1)
    _factor(n // d, out, _depth + 1)


# Recursion depth bound for _factor.  Each level removes at least one bit of
# the cofactor, so ~log2(n) levels suffice; 1024 is far beyond any input whose
# bit length fits in the effort budget.
_MAX_FACTOR_DEPTH = 1024


def _finish_decompose(n: int, c: int, s: int) -> Tuple[int, int]:
    """
    The single exit for `squarefree_decompose`: check the identity, then return.

    Every branch used to `return` directly, so the `c * c * s == n` check ran
    on one of six paths and was documented as if it guarded all of them.  It
    now guards all of them, and `c * c * s == n` is a claim worth checking on
    every exit: it is what makes a misfactorisation an error instead of a
    radicand that silently corrupts a later sign decision.
    """
    if c * c * s != n:
        raise RuntimeError(
            f"squarefree_decompose({n}) is internally inconsistent: "
            f"c={c}, s={s} give c*c*s={c * c * s}"
        )
    return (c, s)


def squarefree_decompose(n: int) -> Tuple[int, int]:
    """
    Exact decomposition n = c² · s with s squarefree.
    Returns (c, s).

    Algorithm:
    1. If n is a perfect square, return (isqrt(n), 1) immediately.
    2. Trial-divide by primes up to 10^6, extracting square factors.
    3. For the remaining cofactor r:
       - r == 1: done
       - r a perfect square: absorb sqrt(r) into c
       - r probable prime (deterministic for r < 2^64): r is squarefree
       - otherwise factor r with Pollard's rho and read off the exponents,
         taking floor(e/2) into c and e mod 2 into s.

    WHAT IS PROVEN, AND WHAT IS NOT
    -------------------------------
    `c * c * s == n` is checked on every return path (`_finish_decompose`) and
    holds unconditionally, whichever primality branch was taken.

    Squarefreeness of `s` is a *consequence of how `s` was built*, not of a
    test run afterwards, and the strength of that argument depends on one
    thing: whether every primality judgement in the chain was deterministic.

      * cofactor r < 2^64 -- Miller-Rabin with the 7-base set is proven
        correct, so r is prime, so `s * r` is a product of distinct primes.
        Squarefreeness is PROVEN.
      * cofactor r >= 2^64 -- the same test is only probable, so a false
        positive would yield a composite r and a non-squarefree `s`.  The
        identity would still hold; only squarefreeness would be unsound.
        The probability is the standard Miller-Rabin bound and is negligible,
        but it is not zero, so this is stated rather than implied.

    An earlier version of this docstring claimed `s` was "verified squarefree
    rather than assumed", which no code in the function did: there was no
    squarefreeness check anywhere on the path.  The honest description is the
    one above -- sound below 2^64, and resting on a probable test above it.

    The one thing that *would* settle it above 2^64 is a primality proof
    (Pocklington, ECPP, or a Pratt certificate), which this module does not
    implement.  `random_squarefree` sidesteps the question instead: it builds
    radicands out of primes that are all below 2^64, so instances generated
    in-repo are in the proven case.
    """
    if n <= 0:
        raise ValueError("n must be positive")

    if n == 1:
        return _finish_decompose(n, 1, 1)

    # Fast path for a whole number that is itself a perfect square: the answer
    # is (isqrt(n), 1) and no factoring is needed.  This is the common case in
    # the F1 family, where radicands are p^2 and alpha*q^2, and it also keeps
    # the trial-division loop below from being handed a 2^40-bit cofactor.
    root = isqrt(n)
    if root * root == n:
        return _finish_decompose(n, root, 1)

    r, c, s = _trial_division(n)

    if r == 1:
        return _finish_decompose(n, c, s)

    # Check perfect square
    if _is_square(r):
        c *= isqrt(r)
        return _finish_decompose(n, c, s)

    # Check prime (deterministic < 2^64, probabilistic above: see docstring)
    if _is_probable_prime(r):
        return _finish_decompose(n, c, s * r)

    # Composite cofactor: factor it exactly and read off the exponents.
    exps: Dict[int, int] = {}
    try:
        _factor(r, exps)
    except RhoFailure as exc:
        # Re-raise with the original radicand attached: the caller passed n,
        # and "a 122-bit cofactor" is much less actionable than "this radicand".
        raise RhoFailure(f"cannot squarefree-decompose {n}: {exc}") from exc
    for p, e in exps.items():
        c *= p ** (e // 2)
        if e % 2:
            s *= p

    return _finish_decompose(n, c, s)


def is_squarefree(n: int) -> bool:
    """
    Exact squarefree test: does any prime square divide n?

    This needs much less than a full factorisation, and the previous version
    asked for one.  It trial-divided by the sieved small primes, then called
    `_factor` on whatever was left, which meant a 200-bit candidate went
    through rho all the way down to primes.  That is why F2 with 200-bit
    radicands raised RhoFailure: deciding squarefreeness of a semiprime does
    not require resolving it into primes, it requires noticing whether any
    prime repeats.

    The cheaper exact recursion is

        n is squarefree  <=>  d squarefree
                            and n/d squarefree
                            and gcd(d, n/d) == 1

    for any proper divisor d, since a prime that divides both d and n/d is
    exactly a prime that divides n at least twice.  So one rho split plus a
    gcd suffices, and each half is strictly smaller, so it terminates.  When
    rho genuinely cannot split a very large cofactor, RhoFailure propagates
    rather than returning a guess.
    """
    if n < 0:
        raise ValueError(f"n must be non-negative, got {n}")
    if n < 4:
        return True              # 0, 1, 2, 3
    if n % 4 == 0:
        return False
    root = isqrt(n)
    if root * root == n:
        return False

    _build_small_primes()
    remaining = n
    for p in _SMALL_PRIMES:
        if p * p > remaining:
            break
        if remaining % (p * p) == 0:
            return False
        while remaining % p == 0:
            remaining //= p

    return _is_squarefree_core(remaining)


def _is_squarefree_core(n: int) -> bool:
    """
    Squarefreeness of a number with no prime factor below 10^6.

    Cheaper than a full factorisation, as documented on `is_squarefree`.
    """
    if n == 1:
        return True
    if _is_square(n):
        return False
    if _is_probable_prime(n):
        return True

    d = _pollard_rho(n)
    other = n // d
    if other == d:
        # n = d^2 exactly, and _is_square already caught that above; reaching
        # here means rho returned a factor that is not what we think, which
        # would be a bug rather than an input problem.
        return False
    if _gcd(d, other) != 1:
        return False
    return _is_squarefree_core(d) and _is_squarefree_core(other)


def _certified_prime(floor: int, limit: int, rng, exclude) -> int:
    """
    A prime p with floor <= p < limit and p not dividing anything in `exclude`.

    Every prime this module *constructs* is below `DETERMINISTIC_MR_CEILING`,
    so the 7-base Miller-Rabin base set applies and the answer is a **proven**
    prime, not a probable one.  That is the whole reason the construction in
    `random_squarefree` uses these factors: it is what turns squarefreeness of
    a generated radicand into a theorem rather than a guess.
    """
    for _ in range(_PRIME_SEARCH_TRIES):
        cand = rng.randrange(floor, limit) | 1
        # `| 1` can round an even value up to `limit`, and the ceiling is
        # exclusive, so the bound is rechecked rather than assumed.
        if cand < limit and cand not in exclude and _is_probable_prime(cand):
            return cand
    raise RuntimeError(
        f"no certified prime in [{floor}, {limit}) after "
        f"{_PRIME_SEARCH_TRIES} draws; retry with a different seed"
    )


def random_squarefree(n_bits: int, rng) -> int:
    """
    A random squarefree integer with exactly `n_bits` bits.

    CONSTRUCTION, NOT REJECTION SAMPLING
    ------------------------------------
    The obvious approach -- draw a random n and test it -- costs about 5
    seconds per unusable candidate at 200 bits, because a random even-ish
    large number almost always has a large prime cofactor that pure-Python
    rho will not split, so the test ends in RhoFailure and the draw is
    wasted.  Measured: generating three 200-bit squarefree numbers that way
    took 21.8 s, essentially all of it in failures.

    Instead build one, as a product of **distinct primes each below
    `DETERMINISTIC_MR_CEILING`**, where Miller-Rabin is proven correct:

        lo = 2^(n_bits-1) <= Q * p < 2^n_bits = hi

    * Q accumulates distinct certified filler primes until
      `Q >= 2^(n_bits - 64)`.  That threshold is load-bearing: it forces
        p <= (hi - 1) // Q <= 2^64 - 1,
      so the final factor is certified by the deterministic base set too, and
      no factor anywhere is ever tested probabilistically.
    * p is then one more certified prime in the window forced by
      `p * Q in [lo, hi)`, namely `[ceil(lo/Q), (hi-1)//Q]`.  A product of
      distinct primes is squarefree, so the result is squarefree *by
      construction* -- there is no test of squarefreeness to pass, and none
      could fail.

    A previous version searched for a single n_bits-bit prime and returned
    it.  Above 64 bits that prime's primality is only **probable**, so the
    squarefreeness it implied was only probable as well, while the docstring
    claimed a verified result.  Splitting the same bit budget across several
    certified factors removes the gap rather than describing it away.

    The post-conditions the construction actually guarantees are re-checked
    before returning, so a bug in the arithmetic above is reported instead of
    shipped.
    """
    if n_bits < 2:
        raise ValueError(f"n_bits must be >= 2, got {n_bits}")
    if n_bits > _MAX_CONSTRUCT_BITS:
        raise ValueError(
            f"n_bits={n_bits} exceeds {_MAX_CONSTRUCT_BITS}. Every factor stays "
            f"below 2^64 so the result is a certified prime product, but the "
            f"count of certified factors grows as n_bits/{_FILLER_BITS}, so "
            f"past this size the sampler is slow rather than unsound. Raise "
            f"_MAX_CONSTRUCT_BITS if that trade is wanted."
        )

    lo = 1 << (n_bits - 1)
    hi = 1 << n_bits                      # exclusive

    # Phase 1: distinct certified fillers until the last factor is forced small.
    used: set = set()
    q = 1
    primes: List[int] = []
    while q.bit_length() < n_bits - 63:
        filler = _certified_prime(_FILLER_FLOOR, _FILLER_CEILING, rng, used)
        used.add(filler)
        primes.append(filler)
        q *= filler

    # Phase 2: one certified prime lands the product inside the target window.
    # p*q in [lo, hi) is equivalent to p in [ceil(lo/q), (hi-1)//q].  The
    # window is never narrow: it is at least (2^(n_bits-1) - 1)/q - 2 wide,
    # and q < 2^(n_bits-32), so that is comfortably above 2^30.  The only
    # clamp needed is p >= 2, since 1 is not a prime to search for.
    p_lo = max(-(-lo // q), 2)            # ceil(lo/q), and at least 2
    p_hi = (hi - 1) // q
    top = _certified_prime(p_lo, p_hi + 1, rng, used)
    used.add(top)
    primes.append(top)

    result = q * top

    # Re-check the post-conditions the argument above rests on.  These are
    # cheap (all factors are < 2^64) and each one is load-bearing:
    _verify_constructed(n_bits, lo, hi, primes, result)
    return result


def _verify_constructed(
    n_bits: int, lo: int, hi: int, primes: List[int], result: int
) -> None:
    """
    Re-check the claims `random_squarefree` makes, so a construction bug is
    reported rather than returned as a radicand.

    Pairwise-distinct certified primes prove squarefreeness, so verifying those
    two properties *is* the squarefreeness proof -- there is nothing further
    to test.  (Running `is_squarefree` here would not add anything: it reaches
    the same Miller-Rabin branch on the same factors, so it is not an
    independent witness, and it would need rho to re-split the product.)
    """
    if len(set(primes)) != len(primes):
        raise RuntimeError(
            f"random_squarefree construction repeated a prime factor: {primes}"
        )
    product = 1
    for p in primes:
        if p >= DETERMINISTIC_MR_CEILING:
            raise RuntimeError(
                f"random_squarefree construction produced a {p.bit_length()}-bit "
                f"factor {p} at or above 2^64, where Miller-Rabin is only "
                f"probabilistic; the result would not be a certified prime "
                f"product"
            )
        if not _is_probable_prime(p):
            raise RuntimeError(
                f"random_squarefree construction produced {p}, which is not "
                f"prime; the result would not be squarefree"
            )
        product *= p
    if product != result:
        raise RuntimeError(
            f"random_squarefree construction: product of factors {product} != "
            f"returned radicand {result}"
        )
    if not (lo <= result < hi):
        raise RuntimeError(
            f"random_squarefree produced {result} "
            f"({result.bit_length()} bits), outside the requested window "
            f"[2^{n_bits - 1}, 2^{n_bits})"
        )


# Below this, the 7-base Miller-Rabin set is a *proof* of primality; at or
# above it the test is only probabilistic.  Every factor this module
# constructs is kept strictly below this line, which is what makes a generated
# radicand a certified prime product.
DETERMINISTIC_MR_CEILING = 2**64

# Width of the filler primes in the construction above.  Small enough that
# finding one costs a handful of Miller-Rabin rounds, large enough that a
# 512-bit radicand needs about 14 of them.
_FILLER_BITS = 32
_FILLER_FLOOR = 1 << (_FILLER_BITS - 1)
_FILLER_CEILING = 1 << _FILLER_BITS

# Bit ceiling for the construction.  This is a *cost* bound, not a soundness
# one: nothing above it fails to be certified, it just needs more certified
# factors (n_bits/_FILLER_BITS of them) and so takes proportionally longer.
_MAX_CONSTRUCT_BITS = 512

# How many candidates to try when searching for a prime in a window.
_PRIME_SEARCH_TRIES = 400


def decompose_and_aggregate(
    terms: List[Tuple[int, int]]
) -> Tuple[Dict[int, int], int]:
    """
    Given terms [(coeff, radicand), ...], perform squarefree decomposition
    and aggregate by squarefree part.

    Each term: coeff * sqrt(radicand) = coeff * c * sqrt(s)
    Aggregates to: sum[sqrt(s)] = total_coeff * sqrt(s)

    Returns (aggregated_dict, max_bitlength) where
      aggregated_dict: {squarefree_part: total_coefficient}
      max_bitlength: maximum bit-length of any input radicand

    Zero test: S = 0 iff all total_coefficients are 0.
    """
    aggregated = Counter()
    max_bits = 0

    for coeff, rad in terms:
        if rad <= 0:
            raise ValueError(f"Radicand must be positive: {rad}")
        if coeff == 0:
            continue

        max_bits = max(max_bits, rad.bit_length())
        c, s = squarefree_decompose(rad)
        aggregated[s] += coeff * c

    return dict(aggregated), max_bits


def is_zero_sum(terms: List[Tuple[int, int]]) -> Tuple[bool, Dict[int, int]]:
    """
    Exact zero test for Σ coeff_i * sqrt(radicand_i).
    Returns (is_zero, aggregated_dict).

    Zero iff all aggregated coefficients are zero.
    """
    agg, _ = decompose_and_aggregate(terms)
    return all(v == 0 for v in agg.values()), agg


def canonical_form(terms: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Convert sum of sqrt terms to canonical list:
    [(coeff, squarefree_part), ...] sorted by squarefree part,
    with zero coefficients removed.
    """
    agg, _ = decompose_and_aggregate(terms)
    return [(v, s) for s, v in sorted(agg.items()) if v != 0]


# ------------------------------------------------------------
# Demo / self-test

if __name__ == "__main__":
    # Quick sanity checks
    assert squarefree_decompose(12) == (2, 3)      # 12 = 2^2 * 3
    assert squarefree_decompose(18) == (3, 2)      # 18 = 3^2 * 2
    assert squarefree_decompose(100) == (10, 1)    # 100 = 10^2
    assert squarefree_decompose(1) == (1, 1)

    # Aggregation
    terms = [(1, 12), (2, 3)]   # sqrt(12) + 2*sqrt(3) = 2√3 + 2√3 = 4√3
    agg, _ = decompose_and_aggregate(terms)
    assert agg == {3: 4}

    # Zero test
    terms = [(1, 12), (-2, 3)]   # 2√3 - 2√3 = 0
    is_zero, agg = is_zero_sum(terms)
    assert is_zero and all(v == 0 for v in agg.values())

    terms = [(1, 12), (-1, 3)]   # 2√3 - √3 = √3 ≠ 0
    is_zero, agg = is_zero_sum(terms)
    assert not is_zero

    print("All basic tests passed.")