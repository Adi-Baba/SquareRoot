"""Tests for squarefree.py — exact decomposition, aggregation, zero test."""

import random
import unittest
from math import isqrt

import squarefree as SF


def brute_force_squarefree(n):
    """Reference squarefree test by trial division, for cross-checking."""
    for p in range(2, isqrt(n) + 1):
        if n % (p * p) == 0:
            return False
    return True


class TestSquarefreeDecompose(unittest.TestCase):
    def test_perfect_squares(self):
        for n in [1, 4, 9, 16, 25, 36, 49, 64, 100, 10000]:
            c, s = SF.squarefree_decompose(n)
            self.assertEqual(s, 1, f"{n}: s={s} expected 1")
            self.assertEqual(c * c, n, f"{n}: c={c}, c^2={c * c}")

    def test_known_values(self):
        cases = {
            2: (1, 2), 3: (1, 3), 5: (1, 5),
            12: (2, 3), 18: (3, 2), 20: (2, 5),
            28: (2, 7), 45: (3, 5), 50: (5, 2),
            72: (6, 2), 98: (7, 2), 200: (10, 2),
        }
        for n, (want_c, want_s) in cases.items():
            got_c, got_s = SF.squarefree_decompose(n)
            self.assertEqual((got_c, got_s), (want_c, want_s), f"n={n}")
            self.assertEqual(got_c * got_c * got_s, n)
            self.assertTrue(SF.is_squarefree(got_s), f"{got_s} not squarefree")

    def test_decompose_identity(self):
        """c^2 * s == n exactly, for a spread of inputs."""
        for n in [1, 2, 3, 5, 6, 7, 10, 11, 12, 13, 14, 15, 17, 19,
                  20, 21, 22, 23, 24, 26, 28, 29, 30, 97, 121, 1000, 999983]:
            c, s = SF.squarefree_decompose(n)
            self.assertTrue(SF.is_squarefree(s), f"{n}: s={s} not squarefree")
            self.assertEqual(c * c * s, n, f"{n}: {c}^2 * {s} != {n}")

    def test_matches_brute_force_on_sampled_values(self):
        """Cross-check is_squarefree against an independent trial-division loop."""
        rng = random.Random(101)
        values = list(range(1, 3000))
        values += [rng.randrange(1, 10 ** rng.randint(2, 6)) for _ in range(1500)]
        for n in values:
            self.assertEqual(
                SF.is_squarefree(n), brute_force_squarefree(n),
                f"is_squarefree({n}) disagrees with trial division",
            )

    def test_hard_cofactor_beyond_trial_division_bound(self):
        """
        p*q with both primes above 10^6: the old code raised here.

        Pollard's rho now factors it, so the decomposition succeeds and the
        identity c^2*s == n holds.
        """
        p, q = 1000003, 1000033
        self.assertTrue(p * q == 1000036000099)
        c, s = SF.squarefree_decompose(p * q)
        self.assertEqual(c, 1)
        self.assertEqual(c * c * s, p * q)
        self.assertTrue(SF.is_squarefree(s))

    def test_hard_cofactor_with_square_factor(self):
        """p^2*q with p,q above 10^6 exercises the exponent bookkeeping."""
        p, q = 1000003, 1000033
        n = p * p * q
        c, s = SF.squarefree_decompose(n)
        self.assertEqual(c, p)
        self.assertEqual(s, q)
        self.assertEqual(c * c * s, n)

    def test_large_perfect_square_fast_path(self):
        """A huge perfect square must not go anywhere near trial division."""
        n = 10 ** 40
        c, s = SF.squarefree_decompose(n * n)
        self.assertEqual((c, s), (n, 1))
        self.assertEqual(c * c * s, n * n)

    def test_near_square_from_pell_identity(self):
        """
        N^2 - 1 is not a square and factors into two large pieces.

        This is the radicand shape the sign decider meets in the F1 family,
        and it is why rho exists.
        """
        n = 10 ** 40
        c, s = SF.squarefree_decompose(n * n - 1)
        self.assertEqual(c * c * s, n * n - 1)
        self.assertTrue(SF.is_squarefree(s))

    def test_rejects_non_positive(self):
        for bad in (0, -1, -100):
            with self.assertRaises(ValueError):
                SF.squarefree_decompose(bad)

    def test_rho_is_deterministic(self):
        """Same input, same factor, every time -- no hidden randomness."""
        n = 1000003 * 1000033
        first = SF.squarefree_decompose(n)
        for _ in range(5):
            self.assertEqual(SF.squarefree_decompose(n), first)


class TestIsSquarefree(unittest.TestCase):
    def test_small_values(self):
        for n in range(0, 200):
            self.assertEqual(SF.is_squarefree(n), brute_force_squarefree(n), f"n={n}")

    def test_squares_are_not_squarefree(self):
        # base=1 gives n=1, which is squarefree by convention; start at 2.
        for base in range(2, 40):
            n = base * base
            self.assertFalse(SF.is_squarefree(n), f"{n} should not be squarefree")

    def test_one_is_squarefree(self):
        self.assertTrue(SF.is_squarefree(1))
        self.assertEqual(SF.squarefree_decompose(1), (1, 1))

    def test_large_squarefree_primes(self):
        for p in (1000003, 1000033, 999983, 104729):
            self.assertTrue(SF.is_squarefree(p), f"prime {p}")

    def test_large_non_squarefree(self):
        p = 1000003
        self.assertFalse(SF.is_squarefree(p * p))
        self.assertFalse(SF.is_squarefree(4 * 1000003))

    def test_rejects_negative(self):
        with self.assertRaises(ValueError):
            SF.is_squarefree(-1)


class TestAggregationAndZeroTest(unittest.TestCase):
    def test_basic_aggregation(self):
        # sqrt(12) + 2sqrt(3) = 2sqrt(3) + 2sqrt(3) = 4sqrt(3)
        agg, _ = SF.decompose_and_aggregate([(1, 12), (2, 3)])
        self.assertEqual(agg, {3: 4})

    def test_zero_sum(self):
        # sqrt(12) - 2sqrt(3) = 0
        is_zero, agg = SF.is_zero_sum([(1, 12), (-2, 3)])
        self.assertTrue(is_zero)
        self.assertTrue(all(v == 0 for v in agg.values()))

    def test_nonzero_sum(self):
        # sqrt(12) - sqrt(3) = sqrt(3)
        is_zero, agg = SF.is_zero_sum([(1, 12), (-1, 3)])
        self.assertFalse(is_zero)
        self.assertEqual(agg, {3: 1})

    def test_canonical_form(self):
        # sqrt(12) + 2sqrt(3) - sqrt(3) = 3sqrt(3)
        self.assertEqual(SF.canonical_form([(1, 12), (2, 3), (-1, 3)]), [(3, 3)])

    def test_zero_test_is_sound_against_decimal_oracle(self):
        """
        The zero test claims S = 0 iff all aggregated coefficients vanish.

        That rests on {sqrt(d) : d distinct squarefree} being linearly
        independent over Q.  Check the non-trivial direction here: a sum
        that the oracle calls non-zero must not be reported as zero.
        """
        from decimal import Decimal as D, getcontext
        getcontext().prec = 200
        rng = random.Random(42)
        for _ in range(400):
            k = rng.randint(2, 5)
            terms = [
                (rng.choice([-5, -2, -1, 1, 2, 5]), rng.randint(1, 60))
                for _ in range(k)
            ]
            is_zero, _ = SF.is_zero_sum(terms)
            val = sum(D(c) * D(a).sqrt() for c, a in terms)
            oracle_zero = val == 0
            self.assertEqual(is_zero, oracle_zero, f"mismatch on {terms}")


class TestEdgeCases(unittest.TestCase):
    def test_single_term(self):
        agg, _ = SF.decompose_and_aggregate([(5, 12)])
        self.assertEqual(agg, {3: 10})

    def test_zero_coefficient(self):
        agg, _ = SF.decompose_and_aggregate([(0, 12), (1, 3)])
        self.assertEqual(agg, {3: 1})

    def test_large_squarefree(self):
        n = 2 * 3 * 5 * 7 * 11
        c, s = SF.squarefree_decompose(n)
        self.assertEqual((c, s), (1, 2310))

    def test_max_bitlength(self):
        # 2^100 is a perfect square, so it aggregates into the rational part
        # as coefficient 2^50 with squarefree part 1, not as a radical.
        agg, max_bits = SF.decompose_and_aggregate([(1, 2 ** 100), (1, 3)])
        self.assertEqual(max_bits, 101)
        self.assertEqual(agg, {1: 2 ** 50, 3: 1})

    def test_rejects_bad_radicand(self):
        with self.assertRaises(ValueError):
            SF.decompose_and_aggregate([(1, 0)])
        with self.assertRaises(ValueError):
            SF.decompose_and_aggregate([(1, -3)])

    def test_empty_terms(self):
        agg, max_bits = SF.decompose_and_aggregate([])
        self.assertEqual(agg, {})
        self.assertEqual(max_bits, 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
