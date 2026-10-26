"""
Tests for signexact.py — exact sign determination.

The oracle here is 200-digit Decimal arithmetic, computed independently of
signexact's integer kernel.  Agreement is checked on random sums, on the
Pell family where the sign is known in closed form from the residue
p^2 - alpha q^2, and on exact-zero cases.
"""

import random
import unittest
from decimal import Decimal as D, getcontext
from fractions import Fraction as F

import signexact as SE

getcontext().prec = 200


def decimal_oracle(terms):
    """Sign of SUM (c, a) computed at 200 significant digits."""
    val = sum(D(c) * D(a).sqrt() for c, a in terms)
    return 0 if val == 0 else (1 if val > 0 else -1)


class TestExactSign(unittest.TestCase):
    def test_zero_sum(self):
        # sqrt(12) - 2sqrt(3) = 2sqrt(3) - 2sqrt(3) = 0
        result = SE.sign_of_sum([(1, 12), (-2, 3)])
        self.assertEqual(result.sign, 0)
        self.assertEqual(result.method, "zero-test")

    def test_positive_sum(self):
        # sqrt(12) - sqrt(3) = sqrt(3) > 0
        result = SE.sign_of_sum([(1, 12), (-1, 3)])
        self.assertEqual(result.sign, 1)
        self.assertIn(result.method, ("certified-interval", "rational"))

    def test_negative_sum(self):
        # 2sqrt(2) - 3sqrt(3) = 2.828 - 5.196 < 0
        result = SE.sign_of_sum([(2, 2), (-3, 3)])
        self.assertEqual(result.sign, -1)

    def test_rational_sum(self):
        # 2 - 3 = -1
        result = SE.sign_of_sum([(1, 4), (-1, 9)])
        self.assertEqual(result.sign, -1)
        self.assertEqual(result.method, "rational")

    def test_empty_sum(self):
        self.assertEqual(SE.sign_of_sum([]).sign, 0)

    def test_zero_coefficient(self):
        self.assertEqual(SE.sign_of_sum([(0, 2), (1, 3)]).sign, 1)

    def test_square_perfect(self):
        # 2 - 2*3 = -4
        self.assertEqual(SE.sign_of_sum([(1, 4), (-2, 9)]).sign, -1)

    def test_deep_cancellation_is_exact(self):
        """S = sqrt(N^2-1) - N is negative but tiny; must not round to 0."""
        n = 10 ** 40
        result = SE.exact_sign([(1, n * n - 1), (-1, n * n)])
        self.assertEqual(result.sign, -1)
        self.assertEqual(result.method, "certified-interval")
        # It genuinely needed many bits, which is the point.
        self.assertGreater(result.witness["bits"], 100)

    def test_witness_stores_ints_not_strings(self):
        """
        Big instances must not blow up on int->str conversion.

        Python refuses int("...") past 4300 digits, so a witness that
        stringified its numerator raised ValueError on inputs of a few
        thousand bits -- exactly the regime this crate exists to explore.
        """
        result = SE.exact_sign([(1, 2), (-1, 3)])
        num, den = result.witness["bracket_num"]
        self.assertIsInstance(num, int)
        self.assertIsInstance(den, int)
        self.assertLess(num, 0)

    def test_scales_to_large_f1_instances(self):
        """
        The F1 family at tens of thousands of input bits.

        The sign is known in closed form from the residue, so this pins the
        decider against ground truth while exercising the refinement loop at
        a size that a Decimal surrogate cannot resolve.
        """
        import f1_theorem as FT
        for n in (200, 1000, 5000):
            a = 2
            p, q = FT.pell_convergents(a, n)[n]
            alpha = a * (a + 1)
            result = SE.exact_sign([(1, p * p), (-1, alpha * q * q)])
            self.assertEqual(result.sign, FT.expected_sign(a, n), f"n={n}")
            # bits needed must scale with the input, which is the separation law
            input_bits = (p * p).bit_length() + (alpha * q * q).bit_length()
            self.assertGreater(result.witness["bits"], input_bits // 2, f"n={n}")


class TestSignAgainstDecimalOracle(unittest.TestCase):
    def test_random_sums_agree_with_oracle(self):
        parts = [2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29, 30, 31]
        rng = random.Random(7)
        for _ in range(600):
            k = rng.randint(1, 4)
            chosen = rng.sample(parts, k)
            coeffs = [rng.choice([-3, -2, -1, 1, 2, 3]) for _ in chosen]
            terms = list(zip(coeffs, chosen))
            self.assertEqual(
                SE.exact_sign(terms).sign, decimal_oracle(terms),
                f"disagreement on {terms}",
            )

    def test_random_large_radicands_agree_with_oracle(self):
        rng = random.Random(11)
        for _ in range(300):
            k = rng.randint(2, 5)
            nums = [rng.randrange(2, 10 ** rng.randint(3, 9)) for _ in range(k)]
            coeffs = [rng.choice([-2, -1, 1, 2]) for _ in nums]
            terms = list(zip(coeffs, nums))
            self.assertEqual(
                SE.exact_sign(terms).sign, decimal_oracle(terms),
                f"disagreement on {terms}",
            )


class TestPellFamilySign(unittest.TestCase):
    """The F1 sign is known in closed form from the exact Pell residue."""

    def test_sign_matches_pell_residue(self):
        import f1_theorem as FT
        for a in (2, 3, 4, 5, 6):
            for n in range(1, 30):
                p, q = FT.pell_convergents(a, n)[n]
                alpha = a * (a + 1)
                residue = p * p - alpha * q * q
                want = 1 if residue > 0 else -1
                got = SE.exact_sign([(1, p * p), (-1, alpha * q * q)]).sign
                self.assertEqual(got, want, f"a={a}, n={n}, residue={residue}")

    def test_residue_sign_agrees_with_parity(self):
        import f1_theorem as FT
        for a in (2, 3, 4, 5):
            for n in range(0, 20):
                self.assertEqual(
                    FT.expected_sign(a, n),
                    SE.exact_sign(
                        [(1, FT.pell_convergents(a, n)[n][0] ** 2),
                         (-1, a * (a + 1) * FT.pell_convergents(a, n)[n][1] ** 2)]
                    ).sign,
                    f"a={a}, n={n}",
                )


class TestKernel(unittest.TestCase):
    def test_bracket_is_exact_for_known_sqrt(self):
        """The integer bracket must satisfy lo^2 <= s <= hi^2 after scaling."""
        for s in (2, 3, 5, 12, 17, 999, 1000003):
            for bits in (1, 4, 17, 64):
                lo, hi = SE.certified_bracket_int(s, bits)
                scale = 1 << bits
                self.assertLessEqual(lo * lo, s * scale * scale)
                self.assertLessEqual(s * scale * scale, hi * hi)
                self.assertEqual(hi - lo, 1)

    def test_bracket_width_is_exactly_one_ulp(self):
        for s in (2, 7, 12345):
            lo, hi = SE.certified_bracket_int(s, 32)
            self.assertEqual(hi - lo, 1)

    def test_zero_and_unit_sqrt(self):
        self.assertEqual(SE.certified_bracket_int(0, 10), (0, 0))
        lo, hi = SE.certified_bracket_int(1, 10)
        scale = 1 << 10
        self.assertEqual(lo, scale)
        self.assertEqual(hi, scale + 1)

    def test_rejects_negative_radicand(self):
        with self.assertRaises(ValueError):
            SE.certified_bracket_int(-1, 10)

    def test_certified_bracket_rational(self):
        lo, hi = SE.certified_bracket({2: 1, 3: -1}, 128)
        self.assertIsInstance(lo, F)
        self.assertIsInstance(hi, F)
        # sqrt(2) - sqrt(3) < 0
        self.assertLess(hi, 0)


class TestValidation(unittest.TestCase):
    def test_rejects_non_positive_radicand(self):
        for bad in (0, -1, -100):
            with self.assertRaises(ValueError):
                SE.exact_sign([(1, bad)])

    def test_rejects_non_integer_terms(self):
        with self.assertRaises(TypeError):
            SE.exact_sign([(1.0, 2)])
        with self.assertRaises(TypeError):
            SE.exact_sign([(1, 2.0)])

    def test_rejects_bad_max_bits(self):
        with self.assertRaises(ValueError):
            SE.exact_sign([(1, 2)], max_bits=0)

    def test_sign_result_rejects_bad_sign(self):
        with self.assertRaises(ValueError):
            SE.SignResult(2, "exact")


class TestSignUnknownIsAnError(unittest.TestCase):
    def test_precision_cap_raises_rather_than_guessing(self):
        """A cap too small to decide must raise, not return a guess."""
        n = 10 ** 40
        with self.assertRaises(SE.SignUnknownError):
            SE.exact_sign([(1, n * n - 1), (-1, n * n)], max_bits=8)

    def test_sign_result_sign_is_always_in_range(self):
        rng = random.Random(3)
        for _ in range(200):
            terms = [(rng.choice([-1, 1]), rng.randint(2, 500)) for _ in range(3)]
            self.assertIn(SE.exact_sign(terms).sign, (-1, 0, 1))


class TestAggregatedEntryPoint(unittest.TestCase):
    def test_matches_term_entry_point(self):
        rng = random.Random(5)
        for _ in range(200):
            agg = {s: rng.choice([-3, -1, 1, 3]) for s in rng.sample(range(2, 40), 3)}
            agg = {s: c for s, c in agg.items() if c != 0}
            if not agg:
                continue
            terms = [(c, s) for s, c in agg.items()]
            self.assertEqual(
                SE.exact_sign_from_aggregated(agg).sign,
                SE.exact_sign(terms).sign,
                f"disagreement on {agg}",
            )

    def test_rejects_bad_aggregated_keys(self):
        with self.assertRaises(ValueError):
            SE.exact_sign_from_aggregated({0: 1})
        with self.assertRaises(TypeError):
            SE.exact_sign_from_aggregated({2: 1.5})


if __name__ == "__main__":
    unittest.main(verbosity=2)
