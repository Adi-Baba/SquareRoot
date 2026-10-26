"""Tests for field.py — multiquadratic field arithmetic."""

import unittest
from fractions import Fraction as F

import field


def eval_poly(coeffs, x):
    """Evaluate an ascending-order coefficient list at x, in Horner form."""
    acc = F(0)
    for c in reversed(coeffs):
        acc = acc * x + c
    return acc


class TestConstruction(unittest.TestCase):
    def setUp(self):
        self.K = field.MultiQuadraticField([2, 3])
        self.K3 = field.MultiQuadraticField([2, 3, 5])

    def test_basis_order(self):
        self.assertEqual(self.K.basis_radicand, [1, 2, 3, 6])

    def test_sorts_and_dedups(self):
        K = field.MultiQuadraticField([3, 2, 3, 2])
        self.assertEqual(K.s, [2, 3])
        self.assertEqual(K.m, 2)

    def test_rejects_empty(self):
        with self.assertRaises(ValueError):
            field.MultiQuadraticField([])

    def test_rejects_part_below_two(self):
        with self.assertRaises(ValueError):
            field.MultiQuadraticField([1, 2])
        with self.assertRaises(ValueError):
            field.MultiQuadraticField([-2, 3])

    def test_rejects_non_squarefree(self):
        for bad in ([4], [4, 5], [8, 11], [12, 13]):
            with self.assertRaises(ValueError, msg=f"{bad} should be rejected"):
                field.MultiQuadraticField(bad)

    def test_rejects_dependent_generators(self):
        # 2*3 = 6, so {2,3,6} is not independent
        with self.assertRaises(ValueError):
            field.MultiQuadraticField([2, 3, 6])
        with self.assertRaises(ValueError):
            field.MultiQuadraticField([2, 8])

    def test_rejects_too_many_generators(self):
        """The basis is a 2**m table, so m must be capped before allocating."""
        primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                  31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73]
        with self.assertRaises(ValueError):
            field.MultiQuadraticField(primes)

    def test_f2_rank_agrees_with_the_exhaustive_subset_scan(self):
        """
        The O(m^2) F_2 rank replaced an O(m * 2^m) scan that tested every
        subset's product for being a perfect square.  The docstring used to
        assert they agree "on every subset of the first six primes" while no
        such check existed anywhere in the repo.  This is that check.

        The exhaustive criterion: a set of squarefree generators is
        independent iff NO nonempty subset has a perfect-square product.  The
        F_2 criterion: the prime-exponent parity vectors are linearly
        independent.  Both are evaluated here, from scratch, and compared on
        every subset of the first six primes.
        """
        import itertools
        from math import isqrt

        primes = [2, 3, 5, 7, 11, 13]

        def exhaustive_dependent(parts):
            """True iff some nonempty subset has a perfect-square product."""
            for r in range(1, len(parts) + 1):
                for combo in itertools.combinations(parts, r):
                    prod = 1
                    for v in combo:
                        prod *= v
                    root = isqrt(prod)
                    if root * root == prod:
                        return True
            return False

        checked = 0
        for r in range(1, len(primes) + 1):
            for combo in itertools.combinations(primes, r):
                checked += 1
                want = exhaustive_dependent(list(combo))
                vectors = [
                    field.MultiQuadraticField._parity_vector(v) for v in combo
                ]
                got = field._rank_f2(vectors) < len(vectors)
                self.assertEqual(
                    got, want,
                    f"{combo}: F_2 rank says dependent={got}, exhaustive scan "
                    f"says dependent={want}"
                )
        # Every subset of six primes, empty one included: 2**6.
        self.assertEqual(checked, 2 ** len(primes) - 1)

    def test_f2_rank_matches_the_constructor_acceptance(self):
        """
        The same equivalence, checked through the public constructor rather
        than the helpers: a set is accepted iff the two criteria agree.
        """
        import itertools

        primes = [2, 3, 5, 7, 11, 13]
        for r in range(2, 5):
            for combo in itertools.combinations(primes, r):
                parts = list(combo)
                dependent = field._rank_f2(
                    [field.MultiQuadraticField._parity_vector(v) for v in parts]
                ) < len(parts)
                if dependent:
                    with self.assertRaises(ValueError, msg=f"{parts} accepted"):
                        field.MultiQuadraticField(parts)
                else:
                    field.MultiQuadraticField(parts)   # must not raise

    def test_large_squarefree_radicand_is_accepted(self):
        """The old O(sqrt(s)) squarefree loop would have hung here."""
        big = 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31 * 37   # squarefree
        K = field.MultiQuadraticField([big, 41])
        self.assertEqual(K.m, 2)


class TestArithmetic(unittest.TestCase):
    def setUp(self):
        self.K = field.MultiQuadraticField([2, 3])
        self.K3 = field.MultiQuadraticField([2, 3, 5])

    def test_basic_arithmetic(self):
        a = self.K.from_sqrt(0)
        b = self.K.from_sqrt(1)
        self.assertEqual(self.K.mul(a, b), {3: F(1)})
        self.assertEqual(self.K.mul(a, a), {0: F(2)})
        self.assertEqual(self.K.mul(b, b), {0: F(3)})

    def test_add_sub(self):
        a = self.K.from_sqrt(0)
        b = self.K.from_sqrt(1)
        self.assertEqual(self.K.add(a, b), {1: F(1), 2: F(1)})
        self.assertEqual(self.K.sub(a, a), {})
        two_a = self.K.mul_scalar(a, F(2))
        self.assertEqual(self.K.sub(two_a, a), a)

    def test_sqr_pow(self):
        a = self.K.from_sqrt(0)
        self.assertEqual(self.K.sqr(a), {0: F(2)})
        self.assertEqual(self.K.pow(a, 4), {0: F(4)})
        self.assertEqual(self.K.pow(a, 3), {1: F(2)})
        self.assertEqual(self.K.pow(a, 0), self.K.one())

    def test_pow_rejects_negative(self):
        with self.assertRaises(ValueError):
            self.K.pow(self.K.from_sqrt(0), -1)

    def test_norm(self):
        # Norm of sqrt2 in Q(sqrt2,sqrt3) is (sqrt2)^2 * (-sqrt2)^2 = 4
        self.assertEqual(self.K.norm(self.K.from_sqrt(0)), F(4))
        self.assertEqual(self.K.norm(self.K.one()), F(1))
        self.assertEqual(self.K.norm({}), F(0))

    def test_norm_of_non_rational_element(self):
        """Norm(sqrt2 - 1) = ((sqrt2-1)(-sqrt2-1))^2 = (-1)^2 = 1."""
        x = self.K.sub(self.K.from_sqrt(0), self.K.one())
        self.assertEqual(self.K.norm(x), F(1))

    def test_norm_of_zero_element_is_zero_not_an_error(self):
        """
        `_clean` strips zero coefficients, so the zero element is the empty
        dict and the norm's product comes out empty too.  That case must keep
        answering 0, since N(0) = 0 is correct -- it is only the *nonzero*
        element with an empty product that indicates a bug.
        """
        self.assertEqual(self.K.norm({}), F(0))
        self.assertEqual(self.K.norm(self.K.zero()), F(0))
        self.assertEqual(self.K.norm({0: F(0)}), F(0))

    def test_norm_raises_rather_than_reporting_zero_for_a_broken_product(self):
        """
        The guard added to `norm` must be load-bearing, so break the field
        arithmetic it guards and confirm it reports instead of answering 0.

        A nonzero element cannot have a norm with no rational part, so the
        raise is unreachable through correct arithmetic.  It is reachable if
        `mul` is wrong, which is exactly the case worth catching: the earlier
        code returned 0 there, which is a confident wrong answer of the shape
        this repo elsewhere treats as a bug.
        """
        original = self.K.mul
        try:
            def broken_mul(x, y):
                # Return a product that keeps an irrational part, which no
                # norm of a field element may have.
                return {1: F(1)}
            self.K.mul = broken_mul
            with self.assertRaises(ValueError) as ctx:
                self.K.norm(self.K.from_sqrt(0))
            self.assertIn("not rational", str(ctx.exception))
        finally:
            self.K.mul = original
        # And the correct arithmetic still works afterwards.
        self.assertEqual(self.K.norm(self.K.from_sqrt(0)), F(4))

    def test_norm_raises_when_a_nonzero_element_yields_no_rational_part(self):
        """The second guard: an empty product from a nonzero element."""
        original = self.K.mul
        try:
            self.K.mul = lambda x, y: {}      # annihilate everything
            with self.assertRaises(ValueError) as ctx:
                self.K.norm(self.K.from_sqrt(0))
            self.assertIn("no rational part", str(ctx.exception))
        finally:
            self.K.mul = original
        self.assertEqual(self.K.norm(self.K.from_sqrt(0)), F(4))

    def test_trace(self):
        self.assertEqual(self.K.trace(self.K.one()), F(4))
        self.assertEqual(self.K.trace(self.K.from_sqrt(0)), F(0))
        # trace of (sqrt2 - 1) is 4 * (coefficient of 1) = -4
        x = self.K.sub(self.K.from_sqrt(0), self.K.one())
        self.assertEqual(self.K.trace(x), F(-4))

    def test_neg_and_zero_one(self):
        a = self.K.from_sqrt(0)
        self.assertEqual(self.K.neg(a), {1: F(-1)})
        self.assertEqual(self.K.zero(), {})
        self.assertEqual(self.K.one(), {0: F(1)})

    def test_mul_scalar(self):
        a = self.K.from_sqrt(0)
        self.assertEqual(self.K.mul_scalar(a, F(5)), {1: F(5)})

    def test_three_generators(self):
        a = self.K3.from_sqrt(0)
        b = self.K3.from_sqrt(1)
        c = self.K3.from_sqrt(2)
        ab = self.K3.mul(a, b)
        self.assertEqual(ab, {3: F(1)})
        self.assertEqual(self.K3.mul(ab, c), {7: F(1)})
        self.assertEqual(self.K3.norm(self.K3.from_sqrt(0)), F(16))

    def test_from_sqrt_range_check(self):
        with self.assertRaises(ValueError):
            self.K.from_sqrt(5)
        with self.assertRaises(ValueError):
            self.K.from_sqrt(-1)

    def test_norm_and_trace_cap_m(self):
        """Large m must raise, not silently return a wrong answer."""
        K = field.MultiQuadraticField([2, 3, 5, 7, 11, 13, 17])
        with self.assertRaises(ValueError):
            K.norm(K.one())
        with self.assertRaises(ValueError):
            K.trace(K.one())

    def test_all_four_exponential_ops_are_capped(self):
        """
        `MAX_GENERATORS` bounds the basis table, not the work the exponential
        operations can be asked to do.  All four must refuse above the tighter
        cap, including the two polynomial methods, which previously had no
        cap at all and so could be run effectively forever.
        """
        K = field.MultiQuadraticField([2, 3, 5, 7, 11, 13, 17])
        x = K.one()
        for name, call in (
            ("norm", lambda: K.norm(x)),
            ("trace", lambda: K.trace(x)),
            ("characteristic_polynomial", lambda: K.characteristic_polynomial(x)),
            ("minimal_polynomial", lambda: K.minimal_polynomial(x)),
        ):
            with self.subTest(op=name):
                with self.assertRaises(ValueError, msg=f"{name} did not refuse"):
                    call()

    def test_the_cap_is_a_separate_constant_from_max_generators(self):
        """
        The two caps are different on purpose, and the polynomial cap is the
        tighter one: a field that constructs at m = 20 must still refuse the
        exponential operations rather than attempt 2**20 conjugates.
        """
        self.assertLess(field._POLYNOMIAL_MAX_GENERATORS, field.MAX_GENERATORS)
        big = field.MultiQuadraticField([2, 3, 5, 7, 11, 13, 17, 19, 23, 29])
        self.assertEqual(big.m, 10)          # construction is fine
        with self.assertRaises(ValueError):
            big.characteristic_polynomial(big.one())

    def test_dependent_subset_helper_does_not_guess_when_independent(self):
        """
        `_find_dependent_subset` used to fall through to
        `return list(range(len(vectors)))`, naming *every* generator as a
        dependent subset.  That is false exactly when the vectors are
        independent, so the unreachable path reported rather than guessed.
        """
        K = field.MultiQuadraticField([2, 3])
        independent = [K._parity_vector(v) for v in (2, 3)]
        with self.assertRaises(ValueError) as ctx:
            K._find_dependent_subset(independent)
        self.assertIn("independent", str(ctx.exception))


class TestPolynomials(unittest.TestCase):
    def setUp(self):
        self.K = field.MultiQuadraticField([2, 3])

    def test_minimal_polynomial_of_sqrt2_is_quadratic(self):
        """
        The key fix: sqrt2 in Q(sqrt2,sqrt3) has degree 2, not 4.

        The old code multiplied over all 2^m conjugates and called the
        result the minimal polynomial, producing (t^2-2)^2.
        """
        s2 = self.K.from_sqrt(0)
        self.assertEqual(self.K.minimal_polynomial(s2), [F(-2), F(0), F(1)])

    def test_characteristic_polynomial_has_multiplicity(self):
        s2 = self.K.from_sqrt(0)
        char = self.K.characteristic_polynomial(s2)
        # (t^2 - 2)^2 = t^4 - 4t^2 + 4
        self.assertEqual(char, [F(4), F(0), F(-4), F(0), F(1)])

    def test_minimal_of_sqrt2_plus_sqrt3(self):
        x = self.K.add(self.K.from_sqrt(0), self.K.from_sqrt(1))
        # t^4 - 10t^2 + 1
        self.assertEqual(self.K.minimal_polynomial(x),
                         [F(1), F(0), F(-10), F(0), F(1)])

    def test_minimal_of_rational_element_is_linear(self):
        self.assertEqual(self.K.minimal_polynomial(self.K.one()), [F(-1), F(1)])
        self.assertEqual(self.K.minimal_polynomial({0: F(7)}), [F(-7), F(1)])

    def test_minimal_of_sqrt2_minus_1(self):
        x = self.K.sub(self.K.from_sqrt(0), self.K.one())
        # t^2 + 2t - 1
        self.assertEqual(self.K.minimal_polynomial(x), [F(-1), F(2), F(1)])

    def test_minimal_polynomial_vanishes_at_the_element(self):
        """
        The defining property, checked without any float: for each of the
        distinct conjugates r, evaluating the polynomial at r must give 0.
        Conjugates are sums of radicals, so this is checked symbolically by
        multiplying the polynomial coefficient-wise in the field.
        """
        K = self.K
        for x in (K.from_sqrt(0),
                  K.from_sqrt(1),
                  K.sub(K.from_sqrt(0), K.one()),
                  K.add(K.from_sqrt(0), K.from_sqrt(1)),
                  K.one(),
                  K.mul_scalar(K.from_sqrt(0), F(3))):
            poly = K.minimal_polynomial(x)
            for conj in K._conjugates(x):
                # Horner in the field, then require the result to be the
                # rational zero element.
                acc = K.zero()
                for coeff in reversed(poly):
                    acc = K.add(K.mul(acc, conj), {0: coeff} if coeff else {})
                self.assertEqual(acc, {}, f"poly does not vanish on a conjugate of {x}")

    def test_conjugate_polynomial_is_the_characteristic_one(self):
        s2 = self.K.from_sqrt(0)
        self.assertEqual(field.conjugate_polynomial(self.K, s2),
                         self.K.characteristic_polynomial(s2))

    def test_distinct_conjugate_masks_deduplicates(self):
        """
        sqrt2 in Q(sqrt2,sqrt3) has orbit {+sqrt2, -sqrt2}: size 2, not 4.

        The count is what separates the minimal polynomial's degree from the
        characteristic polynomial's.
        """
        orbit = self.K.distinct_conjugate_masks(self.K.from_sqrt(0))
        self.assertEqual(len(orbit), 2, f"expected 2 distinct images, got {orbit}")
        self.assertEqual(len(orbit), len(self.K.minimal_polynomial(self.K.from_sqrt(0))) - 1)

    def test_distinct_conjugate_masks_for_full_degree_element(self):
        """sqrt2 + sqrt3 has orbit of size 4, matching its degree."""
        x = self.K.add(self.K.from_sqrt(0), self.K.from_sqrt(1))
        self.assertEqual(len(self.K.distinct_conjugate_masks(x)), 4)
        self.assertEqual(len(self.K.minimal_polynomial(x)) - 1, 4)

    def test_distinct_conjugate_masks_of_zero(self):
        self.assertEqual(self.K.distinct_conjugate_masks({}), [{}])


if __name__ == "__main__":
    unittest.main(verbosity=2)
