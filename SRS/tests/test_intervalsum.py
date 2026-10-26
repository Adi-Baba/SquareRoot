"""
Tests for intervalsum.py — certified rational brackets for square-root sums.

Every bracket is checked against the defining property lo <= S <= hi, verified
by squaring rational endpoints against exact radicands rather than by
comparing to a reference float.
"""

import random
import unittest
from fractions import Fraction as F

import intervalsum as IV
from signexact import exact_sign


def _sign_of_scaled(terms, offset):
    """
    Exact sign of den * (S - offset) as a plain term list, or None if offset
    is not a rational with an integer numerator/denominator.

    S - offset = (1/den) * SUM (c_i*den) sqrt(a_i) - num, and the subtracted
    num is folded in as -num * sqrt(1).  The 1/den is positive, so the sign
    of the scaled expression is the sign of S - offset.
    """
    from signexact import exact_sign
    if offset.denominator == 0:
        return None
    den, num = offset.denominator, offset.numerator
    scaled = [(c * den, a) for c, a in terms]
    if num:
        scaled.append((-num, 1))
    if not scaled:
        return 0
    return exact_sign(scaled).sign


def bracket_contains(terms, lo, hi):
    """
    Independently verify lo <= SUM c_i sqrt(a_i) <= hi.

    Done by asking the exact sign decider about S - lo and S - hi, which is
    an exact integer decision and does not reuse the interval code under
    test.  (An earlier version rebuilt its own sqrt brackets and could not
    confirm a bracket that was tight on an exact value, since any
    interval-based check has strictly positive width of its own.)
    """
    lo = F(lo)
    hi = F(hi)
    if lo > hi:
        return False
    below = _sign_of_scaled(terms, lo)     # sign(S - lo)  must be >= 0
    above = _sign_of_scaled(terms, hi)     # sign(S - hi)  must be <= 0
    return below >= 0 and above <= 0


class TestCertifiedInterval(unittest.TestCase):
    def test_known_signs_and_containment(self):
        cases = [
            ([(1, 12), (-1, 3)], 1),
            ([(1, 2), (-1, 3)], -1),
            ([(1, 2), (1, 3), (-1, 5)], 1),
            ([(2, 2), (-3, 3)], -1),
        ]
        for terms, want in cases:
            iv = IV.certified_interval_for_S(terms, 20)
            self.assertEqual(iv.sign, want, f"wrong sign for {terms}")
            self.assertNotIn(0, (iv.lo, iv.hi))
            # The bracket must not straddle zero once the sign is decided.
            self.assertFalse(iv.contains_zero())
            # And the endpoints must be correctly ordered and tight.
            self.assertLess(iv.lo, iv.hi)
            self.assertTrue(bracket_contains(terms, iv.lo, iv.hi))

    def test_exact_zero_gives_degenerate_bracket(self):
        # sqrt(12) - 2sqrt(3) = 0
        iv = IV.certified_interval_for_S([(1, 12), (-2, 3)], 20)
        self.assertEqual(iv.sign, 0)
        self.assertEqual(iv.lo, F(0))
        self.assertEqual(iv.hi, F(0))
        self.assertTrue(iv.contains_zero())
        self.assertIsNone(IV.digit_count(iv))

    def test_rational_only_sum(self):
        # 2 - 3 = -1, no irrational part
        iv = IV.certified_interval_for_S([(1, 4), (-1, 9)], 20)
        self.assertEqual(iv.sign, -1)
        self.assertEqual(iv.lo, F(-1))
        self.assertEqual(iv.hi, F(-1))
        self.assertEqual(iv.relative_width(), F(0))

    def test_deep_cancellation_bracket_is_tight(self):
        n = 10 ** 30
        iv = IV.certified_interval_for_S([(1, n * n - 1), (-1, n * n)], 20)
        self.assertEqual(iv.sign, -1)
        # Relative width must be small, i.e. the refinement really happened.
        self.assertLess(iv.relative_width(), F(1, 10 ** 20))

    def test_agreement_with_exact_sign(self):
        rng = random.Random(19)
        for _ in range(300):
            k = rng.randint(1, 4)
            terms = [
                (rng.choice([-3, -2, -1, 1, 2, 3]), rng.randint(2, 400))
                for _ in range(k)
            ]
            iv = IV.certified_interval_for_S(terms, 15)
            self.assertEqual(iv.sign, exact_sign(terms).sign, f"on {terms}")
            self.assertTrue(bracket_contains(terms, iv.lo, iv.hi))

    def test_pell_family_agrees_with_residue(self):
        import f1_theorem as FT
        for a in (2, 3, 4):
            for n in range(1, 12):
                p, q = FT.pell_convergents(a, n)[n]
                alpha = a * (a + 1)
                terms = [(1, p * p), (-1, alpha * q * q)]
                iv = IV.certified_interval_for_S(terms, 10)
                self.assertEqual(iv.sign, FT.expected_sign(a, n), f"a={a}, n={n}")
                self.assertTrue(bracket_contains(terms, iv.lo, iv.hi))


class TestRelativeWidthAndDigits(unittest.TestCase):
    def test_relative_width_of_exact_point_is_zero(self):
        self.assertEqual(IV.relative_width(F(3), F(3)), F(0))
        self.assertEqual(IV.relative_width(F(0), F(0)), F(0))

    def test_relative_width_grows_with_bracket(self):
        narrow = IV.relative_width(F(1), F(1) + F(1, 10 ** 6))
        wide = IV.relative_width(F(1), F(2))
        self.assertLess(narrow, wide)

    def test_relative_width_rejects_empty_interval(self):
        with self.assertRaises(ValueError):
            IV.relative_width(F(2), F(1))

    def test_digit_count_none_when_straddling_zero(self):
        iv = IV.CertifiedInterval(F(-1), F(1), 10, 0)
        self.assertIsNone(IV.digit_count(iv))

    def test_digit_count_counts_pinned_digits(self):
        # bracket pins ~30 digits of 1
        iv = IV.CertifiedInterval(F(1), F(1) + F(1, 10 ** 30), 128, 1)
        self.assertGreaterEqual(IV.digit_count(iv), 29)

    def test_relative_width_no_longer_invents_infinity(self):
        """The old slack_ratio returned inf; this is a real number now."""
        lo, hi = IV.certified_interval_for_S([(1, 2), (-1, 3)], 20).lo, \
            IV.certified_interval_for_S([(1, 2), (-1, 3)], 20).hi
        self.assertNotEqual(IV.relative_width(lo, hi), float('inf'))

    def test_relative_width_edge_cases_match_the_documented_values(self):
        """
        The docstring spells out two cases that are not what intuition
        suggests, and got both wrong before.  Asserted here so the
        documentation cannot drift away from the code again.
        """
        # [0,0] is 0, not infinity: a zero-width bracket is an exactly
        # determined value, which is the best outcome, not total uncertainty.
        self.assertEqual(IV.relative_width(F(0), F(0)), F(0))
        # A symmetric straddle is 2, not 1; there is no reason for it to be 1.
        self.assertEqual(IV.relative_width(F(-1), F(1)), F(2))
        # One-sided straddles come out at 1.
        self.assertEqual(IV.relative_width(F(0), F(1)), F(1))
        self.assertEqual(IV.relative_width(F(-1), F(0)), F(1))
        # The scale uses max(|lo|,|hi|), so an asymmetric bracket is not 2.
        self.assertEqual(IV.relative_width(F(1), F(2)), F(1, 2))

    def test_documented_relative_width_returns_are_exact_rationals(self):
        """
        `relative_width` is on the decision path, so its return type is part
        of the contract: an exact Fraction, never a float, including on the
        edge cases where an implementation might be tempted to shortcut.
        """
        for lo, hi in ((F(0), F(0)), (F(-1), F(1)), (F(0), F(1)),
                       (F(-1), F(0)), (F(1), F(2)), (F(3), F(3)),
                       (F(1), F(1) + F(1, 10 ** 6))):
            with self.subTest(lo=lo, hi=hi):
                val = IV.relative_width(lo, hi)
                self.assertIsInstance(
                    val, F,
                    f"relative_width({lo},{hi}) returned {type(val).__name__}, "
                    f"which would be a float in the decision path"
                )

    def test_slack_ratio_is_not_claimed_to_exist(self):
        """
        An earlier docstring told callers `slack_ratio` was "kept for callers
        that want width/error", but no such function existed anywhere in the
        module.  It must not be reintroduced as a stub, and the docstring must
        not advertise it as available.
        """
        self.assertFalse(
            hasattr(IV, "slack_ratio"),
            "slack_ratio is documented as unanswerable; do not add a stub"
        )
        doc = IV.__doc__ or ""
        for claim in ("`slack_ratio` is kept",
                      "slack_ratio is kept",
                      "kept for callers",
                      "slack_ratio(lo"):
            self.assertNotIn(claim, doc, f"docstring still advertises {claim!r}")


class TestTermInterval(unittest.TestCase):
    def test_negative_coefficient_swaps_endpoints(self):
        lo, hi = IV.certified_interval_for_term(-1, 2, 20)
        self.assertLess(lo, hi)          # -sqrt(2) < 0 < upper
        self.assertLess(hi, 0)

    def test_positive_coefficient_keeps_order(self):
        lo, hi = IV.certified_interval_for_term(3, 2, 20)
        self.assertLess(lo, hi)
        self.assertGreater(lo, 0)


class TestValidation(unittest.TestCase):
    def test_rejects_bad_terms(self):
        with self.assertRaises(ValueError):
            IV.certified_interval_for_S([(1, 0)], 10)
        with self.assertRaises(ValueError):
            IV.certified_interval_for_S([(1, -5)], 10)
        with self.assertRaises(TypeError):
            IV.certified_interval_for_S([(1.5, 2)], 10)

    def test_precision_cap_raises(self):
        n = 10 ** 40
        with self.assertRaises(Exception):
            IV.certified_interval_for_S([(1, n * n - 1), (-1, n * n)], 10, max_bits=8)


if __name__ == "__main__":
    unittest.main(verbosity=2)
