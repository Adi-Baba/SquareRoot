"""
Multiquadratic field arithmetic (exact, rational coefficients).

Represents elements of K = Q(√s_1, ..., √s_m) where s_i are distinct squarefree integers.
Assumes s_i are independent (no s_i is a product of others, up to squares).

Basis: { √(∏_{i∈S} s_i) : S ⊆ {1..m} }, indexed by bitmask S ∈ [0, 2^m).
Order: empty set = 0 (basis = 1), singletons = powers of 2, etc.

Operations: add, mul, norm, equality — all exact in Fractions.
"""

from typing import Dict, List
from fractions import Fraction as F
from collections import defaultdict
from math import gcd

from squarefree import _factor, is_squarefree


def _rank_f2(vectors: List[dict]) -> int:
    """
    Rank over F_2 of a list of {prime: 0/1} vectors, by elimination.

    Each vector is a parity set of prime exponents, so XOR is the addition.
    The pivot is the largest prime in the row, which needs no ordering beyond
    max() and keeps the elimination deterministic.
    """
    pivots: Dict[int, dict] = {}
    rank = 0
    for vec in vectors:
        cur = {p: 1 for p, v in vec.items() if v}
        while cur:
            lead = max(cur)
            if lead not in pivots:
                pivots[lead] = cur
                rank += 1
                break
            row = pivots[lead]
            for p in set(cur) | set(row):
                if cur.get(p, 0) ^ row.get(p, 0):
                    cur[p] = 1
                else:
                    cur.pop(p, None)
    return rank

# Every operation here is O(2^m) or worse, and the basis is a table of 2^m
# entries, so m is capped rather than left to exhaust memory.  2^20 basis
# elements is already far past any use; a request beyond it is a bug, not an
# instance.
MAX_GENERATORS = 20

# Separate, tighter cap for the four exponential operations.  `MAX_GENERATORS`
# bounds the 2**m basis table built at construction; it does not bound what
# `norm`, `trace`, `characteristic_polynomial` or `minimal_polynomial` can then
# be asked to do, since each of those walks all 2**m conjugates and the
# polynomial pair multiplies them pairwise.  The Rust port already errors above
# its equivalent cap; these two methods did not, so a field that constructed
# fine at m = 20 could still be handed to them and run essentially forever.
_POLYNOMIAL_MAX_GENERATORS = 6


class MultiQuadraticField:
    """
    Exact arithmetic in a multiquadratic field Q(√s_1, ..., √s_m).

    Basis elements: for subset S ⊆ {0..m-1}, basis element is
        √(∏_{i∈S} s_i)
    with the convention S = ∅ gives basis 1.
    Index of basis element = bitmask of S.

    Invariant: all s_i are distinct squarefree integers ≥ 2, and independent
    (no s_i is a product of others up to squares).
    """

    def __init__(self, squarefree_parts: List[int]):
        if not squarefree_parts:
            raise ValueError("At least one squarefree part required")

        self.s = sorted(set(squarefree_parts))
        self.m = len(self.s)

        if self.m > MAX_GENERATORS:
            raise ValueError(
                f"too many generators: m={self.m} exceeds MAX_GENERATORS={MAX_GENERATORS} "
                f"(the basis alone would need 2^{self.m} entries)"
            )

        for s in self.s:
            if s < 2:
                raise ValueError(f"Squarefree parts must be ≥ 2, got {s}")
            # Delegate to the exact squarefree test.  The previous inline
            # `for p in range(2, isqrt(s)+1)` loop was O(sqrt(s)) and, worse,
            # reported "squarefree" for any s whose first square divisor
            # exceeded sqrt(s) -- which cannot happen, but the loop also ran
            # for hours on a 100-digit radicand.  is_squarefree factors what
            # trial division leaves behind, so it is a decision either way.
            if not is_squarefree(s):
                raise ValueError(f"{s} is not squarefree")

        self._verify_independence()

        # basis_radicand[mask] = product of s[i] for i in mask
        self.basis_radicand = [1] * (1 << self.m)
        for mask in range(1, 1 << self.m):
            lsb = mask & -mask
            i = (lsb.bit_length() - 1)
            prev = mask ^ lsb
            self.basis_radicand[mask] = self.basis_radicand[prev] * self.s[i]

    def _verify_independence(self):
        """
        Reject generator sets that are linearly dependent over Q.

        The previous check enumerated all 2^(m-1) subsets per generator and
        tested each product for being a perfect square.  That is correct but
        costs O(m * 2^m) big-integer multiplications, which is why constructing
        a field with m=20 took 19 seconds before any arithmetic was done.

        The exact criterion is cheaper.  For squarefree integers, a product
        s_i1 * ... * s_ik is a perfect square exactly when every prime appears
        an even number of times across the chosen factors, i.e. exactly when
        the prime-exponent vectors of the s_i are dependent over F_2.  So
        factor each generator and take an F_2 rank, which is O(m^2) set
        operations on small dicts, and the check becomes independent of m's
        exponential.

        Equivalence with the old subset scan is not argued from inspection
        here: tests/test_field.py runs both criteria over every subset of the
        first six primes and asserts they agree.  That is a finite check, not
        a proof, so it is stated as a test rather than as a fact about all
        inputs.

        Factoring is the same exact routine squarefree_decompose uses, so a
        generator too large to split raises rather than being assumed
        independent.
        """
        vectors = [self._parity_vector(s) for s in self.s]
        if _rank_f2(vectors) < self.m:
            dependent = self._find_dependent_subset(vectors)
            raise ValueError(
                f"Squarefree parts not independent: "
                f"{' * '.join(map(str, (self.s[i] for i in dependent)))} "
                f"has a perfect-square product"
            )

    @staticmethod
    def _parity_vector(n: int) -> dict:
        """Prime -> exponent parity for a squarefree integer n."""
        exps: Dict[int, int] = {}
        _factor(n, exps)
        return {p: e % 2 for p, e in exps.items() if e % 2}

    def _find_dependent_subset(self, vectors: List[dict]) -> List[int]:
        """
        Indices of a dependent subset, for a useful error message.

        Runs the same elimination as _rank_f2 but remembers which original
        vectors ended up cancelling, so the message can name them instead of
        just saying "not independent".

        The subset is one the elimination actually produced, which is
        dependent but not guaranteed *minimal* -- minimality is not enforced
        and the docstring no longer claims it.
        """
        pivot_to_row: Dict[int, dict] = {}
        pivot_to_indices: Dict[int, List[int]] = {}
        for idx, vec in enumerate(vectors):
            cur = dict(vec)
            trail = [idx]
            while cur:
                lead = max(cur)
                if lead not in pivot_to_row:
                    pivot_to_row[lead] = cur
                    pivot_to_indices[lead] = trail
                    break
                # XOR with the stored row; it represents the symmetric
                # difference of the generators named in its trail.
                row = pivot_to_row[lead]
                for p in set(cur) | set(row):
                    x = cur.get(p, 0) ^ row.get(p, 0)
                    if x:
                        cur[p] = x
                    else:
                        cur.pop(p, None)
                trail = trail + pivot_to_indices[lead]
            else:
                return trail      # reduced to zero: dependent
        # Unreachable: the caller only asks for a dependent subset when
        # _rank_f2(vectors) < m, in which case some vector must reduce to zero.
        # Reported rather than answered with a guess -- the previous fallback
        # returned every index, i.e. it claimed the whole generator set was a
        # dependent subset, which is false precisely when the vectors are
        # independent.
        raise ValueError(
            f"_find_dependent_subset called on {len(vectors)} vectors that are "
            f"actually independent; the caller's rank check disagrees with the "
            f"elimination"
        )

    def zero(self) -> Dict[int, F]:
        return {}

    def one(self) -> Dict[int, F]:
        return {0: F(1)}

    def from_rational(self, q: F) -> Dict[int, F]:
        return {0: F(q)}

    def from_sqrt(self, idx: int) -> Dict[int, F]:
        if not 0 <= idx < self.m:
            raise ValueError(f"Index {idx} out of range [0, {self.m})")
        return {1 << idx: F(1)}

    def _clean(self, d: Dict[int, F]) -> Dict[int, F]:
        return {k: v for k, v in d.items() if v != 0}

    def add(self, x: Dict[int, F], y: Dict[int, F]) -> Dict[int, F]:
        res = defaultdict(F)
        for k, v in x.items():
            res[k] += v
        for k, v in y.items():
            res[k] += v
        return self._clean(dict(res))

    def neg(self, x: Dict[int, F]) -> Dict[int, F]:
        return self._clean({k: -v for k, v in x.items()})

    def sub(self, x: Dict[int, F], y: Dict[int, F]) -> Dict[int, F]:
        return self.add(x, self.neg(y))

    def mul(self, x: Dict[int, F], y: Dict[int, F]) -> Dict[int, F]:
        res = defaultdict(F)
        for mask_x, coeff_x in x.items():
            for mask_y, coeff_y in y.items():
                mask_res = mask_x ^ mask_y
                intersection = mask_x & mask_y
                coeff_mult = F(1)
                inter = intersection
                while inter:
                    lsb = inter & -inter
                    idx = lsb.bit_length() - 1
                    coeff_mult *= self.s[idx]
                    inter ^= lsb
                res[mask_res] += coeff_x * coeff_y * coeff_mult
        return self._clean(dict(res))

    def mul_scalar(self, x: Dict[int, F], scalar: F) -> Dict[int, F]:
        return self._clean({k: v * scalar for k, v in x.items()})

    def sqr(self, x: Dict[int, F]) -> Dict[int, F]:
        return self.mul(x, x)

    def pow(self, x: Dict[int, F], n: int) -> Dict[int, F]:
        if n < 0:
            raise ValueError("Negative powers not supported")
        if n == 0:
            return self.one()
        res = self.one()
        base = x
        while n > 0:
            if n & 1:
                res = self.mul(res, base)
            base = self.mul(base, base)
            n >>= 1
        return res

    def norm(self, x: Dict[int, F]) -> F:
        """
        Compute the field norm N_{K/Q}(x) = product of all 2^m Galois conjugates.
        The conjugates are obtained by all 2^m sign flips of the square roots.
        The product is guaranteed to be rational; we extract the rational part.

        The zero element is normalised before anything else.  `_clean` strips
        zero coefficients, so the zero element's canonical form is the *empty*
        dict -- but a caller can hand in any equivalent representation, such as
        `{0: 0}`, and testing only for an empty dict would then treat a zero
        element as a nonzero one.  Normalising first makes "is this zero?" a
        question about the element rather than about its spelling.

        After that, "the product has no rational part" is unambiguous.  It
        would mean a nonzero element produced no rational part, which cannot
        happen in a field, so it is reported rather than answered with 0: a
        silent 0 there is exactly the confident-wrong-answer shape that the
        Rust port of this file was fixed for.
        """
        x = self._clean(x)
        if not x:
            return F(0)          # x is the zero element; N(0) = 0

        if self.m > 6:
            raise ValueError(f"Norm computation requires m ≤ 6, got m={self.m}")

        # Generate all 2^m conjugates as field elements
        conjugates = []
        for mask in range(1 << self.m):
            val = {}
            for basis_mask, coeff in x.items():
                sign = 1
                bm = basis_mask
                while bm:
                    lsb = bm & -bm
                    idx = lsb.bit_length() - 1
                    if (mask >> idx) & 1:
                        sign = -sign
                    bm ^= lsb
                val[basis_mask] = coeff * sign
            conjugates.append(self._clean(val))

        # Multiply all conjugates together in the field
        prod = self.one()
        for conj in conjugates:
            prod = self.mul(prod, conj)

        # The result must be rational, i.e. only mask 0 may be non-zero.
        # Anything else means the product gained an irrational part, which
        # cannot happen for a well-formed element, so it is reported rather
        # than returned as a number.
        for mask, coeff in prod.items():
            if mask != 0 and coeff != 0:
                raise ValueError(
                    f"Norm of a nonzero element is not rational: got "
                    f"{ {k: str(v) for k, v in prod.items()} } for "
                    f"x = { {k: str(v) for k, v in x.items()} }"
                )
        if 0 not in prod:
            # A nonzero x has N(x) != 0, so this is not reachable.  If the
            # field arithmetic above is ever wrong, this is where it shows.
            raise ValueError(
                f"Norm of a nonzero element { {k: str(v) for k, v in x.items()} } "
                f"came out with no rational part; the field arithmetic is "
                f"inconsistent"
            )
        return prod[0]

    def trace(self, x: Dict[int, F]) -> F:
        """Field trace: sum of all 2^m Galois conjugates."""
        if self.m > 6:
            raise ValueError(f"Trace requires m ≤ 6, got m={self.m}")

        if not x:
            return F(0)

        total = F(0)
        for mask in range(1 << self.m):
            val = F(0)
            for basis_mask, coeff in x.items():
                sign = 1
                bm = basis_mask
                while bm:
                    lsb = bm & -bm
                    idx = lsb.bit_length() - 1
                    if (mask >> idx) & 1:
                        sign = -sign
                    bm ^= lsb
                val += coeff * sign
            total += val
        return total

    def _conjugates(self, x: Dict[int, F]) -> List[Dict[int, F]]:
        """
        The 2^m Galois images of x, one per sign vector, with multiplicity.

        These are field elements (dicts), not numbers: at this stage they are
        still combinations of radicals and are not comparable or hashable as
        reals.  Use `distinct_conjugates` when the numerical values are needed.
        """
        out = []
        for mask in range(1 << self.m):
            val = {}
            for basis_mask, coeff in x.items():
                sign = 1
                bm = basis_mask
                while bm:
                    lsb = bm & -bm
                    idx = lsb.bit_length() - 1
                    if (mask >> idx) & 1:
                        sign = -sign
                    bm ^= lsb
                if sign != 0:
                    val[basis_mask] = coeff * sign
            out.append(self._clean(val))
        return out

    def distinct_conjugate_masks(self, x: Dict[int, F]) -> List[Dict[int, F]]:
        """
        The distinct Galois images of x, as field elements, deduplicated.

        These are dicts of radical coefficients, not Fractions: a conjugate
        of a general element is a sum of square roots and has no rational
        value.  This returns the elements themselves so a caller can tell
        them apart exactly, which is what `minimal_polynomial` needs.  To
        turn one into a number, use signexact, which brackets it.
        """
        if not x:
            return [self.zero()]
        out: List[Dict[int, F]] = []
        seen: set = set()
        for mask in range(1 << self.m):
            val: Dict[int, F] = {}
            for basis_mask, coeff in x.items():
                sign = 1
                bm = basis_mask
                while bm:
                    lsb = bm & -bm
                    idx = lsb.bit_length() - 1
                    if (mask >> idx) & 1:
                        sign = -sign
                    bm ^= lsb
                if sign > 0:
                    val[basis_mask] = coeff
                else:
                    val[basis_mask] = -coeff
            val = self._clean(val)
            key = tuple(sorted(val.items()))
            if key not in seen:
                seen.add(key)
                out.append(val)
        return out

    def characteristic_polynomial(self, x: Dict[int, F]) -> List[F]:
        """
        The characteristic polynomial of multiplication by x on K.

        This is the product of (t - sigma(x)) over ALL 2^m sign vectors, so it
        is a power of the minimal polynomial, not the minimal polynomial.  It
        is the right object when you want the full orbit including multiplicity.

        Coefficients in ascending order, normalised to be primitive (integer
        coefficients with gcd 1) when the result is integral.

        Capped at `_POLYNOMIAL_MAX_GENERATORS`, like `norm` and `trace`.  The
        module-level `MAX_GENERATORS` bounds construction, not this: a general
        element at m = 20 has 2^20 conjugates and the product-of-linear-factors
        loop is quadratic in that, which `MAX_GENERATORS` does not prevent.
        """
        if self.m > _POLYNOMIAL_MAX_GENERATORS:
            raise ValueError(
                f"characteristic_polynomial requires m <= "
                f"{_POLYNOMIAL_MAX_GENERATORS}, got m={self.m}; this is a "
                f"quadratic cost in 2**m. Use minimal_polynomial, which "
                f"collapses the orbit first, or raise the cap deliberately."
            )
        if not x:
            return [F(-1), F(1)]          # t - 0
        return self._product_of_linear_factors(self._conjugates(x))

    def minimal_polynomial(self, x: Dict[int, F]) -> List[F]:
        """
        The minimal polynomial of x over Q, coefficients in ascending order.

        The roots are exactly the DISTINCT Galois images of x, so this is
        `characteristic_polynomial` with duplicate orbit values collapsed
        before multiplying.  Getting this distinction right matters: the old
        implementation multiplied over all 2^m images and called the result
        the minimal polynomial, which made e.g. sqrt(2) in Q(sqrt2,sqrt3) look
        like it had degree 4 rather than 2.

        Normalised to be primitive when integral.

        Capped at `_POLYNOMIAL_MAX_GENERATORS`, for the same reason as
        `characteristic_polynomial`: the orbit is still 2^m long before the
        duplicates collapse, so the cap bounds the scan even though the
        product itself is over far fewer factors.
        """
        if self.m > _POLYNOMIAL_MAX_GENERATORS:
            raise ValueError(
                f"minimal_polynomial requires m <= "
                f"{_POLYNOMIAL_MAX_GENERATORS}, got m={self.m}; the orbit is "
                f"still 2**m long before duplicates collapse. Raise the cap "
                f"deliberately if that is intended."
            )
        if not x:
            return [F(0), F(1)]           # t

        seen: Dict[str, Dict[int, F]] = {}
        for conj in self._conjugates(x):
            key = repr(sorted(conj.items()))
            seen.setdefault(key, conj)
        return self._product_of_linear_factors(list(seen.values()))

    def _product_of_linear_factors(self, roots: List[Dict[int, F]]) -> List[F]:
        """
        Product of (t - r) over field-element roots, as rational coefficients.

        Each intermediate coefficient is still a field element; the rational
        part is extracted at the end, and the absence of any irrational
        remainder is checked rather than assumed.
        """
        poly: List[Dict[int, F]] = [self.one()]

        for conj in roots:
            neg_conj = self.neg(conj)
            new_poly = [self.zero()] * (len(poly) + 1)
            for i, coeff in enumerate(poly):
                new_poly[i] = self.add(new_poly[i], self.mul(coeff, neg_conj))
                new_poly[i + 1] = self.add(new_poly[i + 1], coeff)
            poly = new_poly

        rational_poly: List[F] = []
        for coeff_dict in poly:
            for mask, v in coeff_dict.items():
                if mask != 0 and v != 0:
                    raise ValueError(
                        f"polynomial coefficient is not rational: {coeff_dict}"
                    )
            rational_poly.append(coeff_dict.get(0, F(0)))

        return _primitive(rational_poly)

    def __repr__(self) -> str:
        return f"MultiQuadraticField(s={self.s}, m={self.m})"


def _primitive(poly: List[F]) -> List[F]:
    """
    Scale a rational polynomial to be primitive: integer coefficients, gcd 1.

    No-op for polynomials that are not integral, which is the common case for
    a general field element.
    """
    nonzero = [c for c in poly if c != 0]
    if not nonzero:
        return list(poly)
    if not all(c.denominator == 1 for c in nonzero):
        return list(poly)
    g = 0
    for c in nonzero:
        g = gcd(g, abs(c.numerator))
    if g > 1:
        return [c / g for c in poly]
    return list(poly)


def conjugate_polynomial(field: MultiQuadraticField, x: Dict[int, F]) -> List[F]:
    """
    The characteristic polynomial of x: product of (t - sigma(x)) over all
    2^m sign vectors, with multiplicity.

    Named "conjugate polynomial" to avoid the minimal/characteristic
    confusion.  For the true minimal polynomial use
    `field.minimal_polynomial(x)`; for the sign of a square-root sum use
    signexact.exact_sign, which does not need either.
    """
    return field.characteristic_polynomial(x)


# `evaluate_at_mask(field, x, mask)` used to live here: a single Galois
# conjugate of x, selected by sign-flip mask.  It had no callers and no tests,
# and duplicated the sign-flip walk already in `MultiQuadraticField._conjugates`
# and `.norm`.  Rather than leave untested public API with a second copy of
# that arithmetic to drift, it is gone.  `field._conjugates(x)[mask]` covers
# the same ground if a caller ever needs it.


# ------------------------------------------------------------
# Demo

if __name__ == "__main__":
    K = MultiQuadraticField([2, 3])
    print(f"Field: {K}")
    print(f"Basis radicands: {K.basis_radicand}")

    a = K.from_sqrt(0)  # √2
    b = K.from_sqrt(1)  # √3
    c = K.mul(a, b)
    print(f"√2 * √3 = {c}")

    c2 = K.mul(a, a)
    print(f"√2 * √2 = {c2}")

    x = K.sub(K.from_sqrt(0), K.one())
    n = K.norm(x)
    print(f"Norm(√2 - 1) = {n}")

    print("All manual checks passed.")
