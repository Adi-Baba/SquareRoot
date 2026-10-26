"""Tests for certsingle.py — certified single-root bracket adapter."""

import unittest
from decimal import Decimal as D

import certsingle


class TestCertifiedSqrt(unittest.TestCase):
    def test_basic_certified_sqrt(self):
        for alpha in (2, 3, 5, 7, 12, 17):
            lo, hi, depth = certsingle.certified_sqrt(alpha, 100)
            self.assertIsInstance(lo, D)
            self.assertIsInstance(hi, D)
            self.assertIsInstance(depth, int)
            self.assertGreaterEqual(lo, 0)
            self.assertLessEqual(lo, hi)
            self.assertLessEqual(hi - lo, D(10) ** -100)

    def test_certification_passes(self):
        """The exact-integer proof must accept the produced bracket."""
        for alpha in (2, 3, 5, 7, 12, 17, 999983):
            lo, hi, _ = certsingle.certified_sqrt(alpha, 50)
            msg = certsingle.certify(lo, hi, alpha, 50)
            self.assertIsNone(msg, f"alpha={alpha}: {msg}")

    def test_certified_sqrt_certified_raises_on_failure(self):
        """
        The wrapper must actually run the proof, not just call zoom.

        Patch zoom to hand back a bracket nowhere near sqrt(12) and confirm
        the wrapper raises instead of passing the bad bracket through.  This
        is what makes a returned value mean "verified".
        """
        from unittest import mock
        with mock.patch.object(certsingle.EC, 'zoom',
                               return_value=(D(5), D(6), 1)):
            with self.assertRaises(ValueError) as ctx:
                certsingle.certified_sqrt_certified(12, 50)
        self.assertIn("certification failed", str(ctx.exception))

    def test_deep_precision(self):
        """A demanding target must still certify rather than stall."""
        lo, hi, depth = certsingle.certified_sqrt_certified(12, 2000)
        self.assertIsNone(certsingle.certify(lo, hi, 12, 2000))
        self.assertLessEqual(hi - lo, D(10) ** -2000)

    def test_get_certified_bracket(self):
        result = certsingle.get_certified_bracket(12, 50)
        self.assertEqual(result['alpha'], 12)
        self.assertEqual(result['digits'], 50)
        self.assertTrue(result['certified'])
        self.assertIsNone(result['cert_msg'])

    def test_get_certified_bracket_reports_failure(self):
        """A rejected bracket must be reported, not raised away."""
        lo, hi, _ = certsingle.certified_sqrt(2, 50)
        msg = certsingle.certify(D(100), D(101), 2, 50)
        self.assertIsNotNone(msg, "certify accepted a bracket nowhere near sqrt(2)")

    def test_bracket_to_rational(self):
        lo, hi, _ = certsingle.certified_sqrt(12, 50)
        rational = certsingle.certified_bracket_to_rational(lo, hi, 12, 50)
        for key in ('lo_num', 'lo_den', 'hi_num', 'hi_den'):
            self.assertIn(key, rational)
            self.assertIsInstance(rational[key], int)
        # as_integer_ratio is exact, so the ratio must rebuild the endpoint
        # bit for bit.  Compare as Fractions, not by dividing Decimals: the
        # division would round at the ambient context precision.
        from fractions import Fraction
        self.assertEqual(Fraction(rational['lo_num'], rational['lo_den']),
                         Fraction(*lo.as_integer_ratio()))
        self.assertEqual(Fraction(rational['hi_num'], rational['hi_den']),
                         Fraction(*hi.as_integer_ratio()))

    def test_different_precisions(self):
        for digits in (10, 50, 100, 200):
            lo, hi, depth = certsingle.certified_sqrt(2, digits)
            self.assertLessEqual(hi - lo, D(10) ** -digits)
            self.assertGreater(depth, 0)
            self.assertLess(depth, 100)

    def test_zero_and_perfect_square(self):
        for alpha in (0, 1, 4, 9, 100):
            lo, hi, _ = certsingle.certified_sqrt(alpha, 30)
            self.assertIsNone(certsingle.certify(lo, hi, alpha, 30), f"alpha={alpha}")
            self.assertLessEqual(lo, D(alpha).sqrt() + D(10) ** -25)
            self.assertGreaterEqual(hi, D(alpha).sqrt() - D(10) ** -25)

    def test_large_alpha(self):
        """
        Radicands with many digits before the decimal point.

        The zoom used to cap working precision at digits+GUARD, which cannot
        hold the integer part, and raised "zoom stalled" for anything past
        ~2**100.  The F1 family's radicands are exactly this shape.
        """
        for k in (64, 200, 800):
            alpha = (1 << k) + 1
            lo, hi, _ = certsingle.certified_sqrt(alpha, 40)
            self.assertIsNone(certsingle.certify(lo, hi, alpha, 40), f"2**{k}+1")

    def test_zero_alpha(self):
        """sqrt(0) used to divide by zero inside the interval evaluation."""
        lo, hi, depth = certsingle.certified_sqrt(0, 30)
        self.assertEqual((lo, hi), (D(0), D(0)))
        self.assertIsNone(certsingle.certify(lo, hi, 0, 30))

    def test_rejects_negative_alpha(self):
        with self.assertRaises(ValueError):
            certsingle.certified_sqrt(-1, 20)

    def test_rejects_bad_params(self):
        with self.assertRaises(ValueError):
            certsingle.certified_sqrt(-1, 20)
        with self.assertRaises(ValueError):
            certsingle.certified_sqrt(2, 0)
        with self.assertRaises(ValueError):
            certsingle.certified_sqrt(2, -5)
        with self.assertRaises(TypeError):
            certsingle.certified_sqrt(2.0, 20)
        with self.assertRaises(TypeError):
            certsingle.certified_sqrt(2, 20.0)
        with self.assertRaises(TypeError):
            certsingle.certified_sqrt(True, 20)

    def test_bisection_alias_exists(self):
        """_deps re-exports `bisection`; exp_certified has no `bisect`."""
        import _deps
        self.assertIs(_deps.bisection, certsingle.bisection)


if __name__ == "__main__":
    unittest.main(verbosity=2)
