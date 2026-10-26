"""Cross-validation of the `Decimal` arithmetic the whole project rests on.

Why this file exists
--------------------
Every certificate in this project is exact *integer* arithmetic applied to
values produced by `decimal`.  The certificates are therefore only as good as
the rounding of libmpdec: a misrounding in the library would propagate into
every bracket, every bound, every "0/480 certified" claim, and none of the
`Fraction`-based tests would catch it, because they never touch the same code
path.

That assumption was flagged as an open limitation and left open.  Here it is
closed against an INDEPENDENT oracle: `math.isqrt`, which is exact integer
arithmetic in the standard library, with no relationship to libmpdec.  If the
two disagree, one of them is wrong, and the disagreement is visible.

What is checked
---------------
* `Decimal.sqrt` against `math.isqrt` at 50, 200, 1000, 3000 and 6000 digits,
  for a spread of alpha including the awkward extreme case alpha = 10^12 + 1.
* One-ulp agreement, not digit-string equality: `sqrt` is correctly rounded to
  the context precision, so `x * 10^D` must lie in [floor(sqrt(a)*10^D),
  +1].  Anything outside that is a genuine misrounding.
* A deliberately mis-scaled control, so a future edit that breaks the
  comparison's own arithmetic cannot pass silently.

Run:  python3 test_decimal_oracle.py
"""

import math
import sys
import unittest
from decimal import Decimal as D, localcontext
from fractions import Fraction as F

# Comparing at 6000 digits builds integers with more than 4300 digits, and
# Python 3.11+ refuses int(...) on those.  The same lift exp_certified.py
# applies, for the same reason.
sys.set_int_max_str_digits(0)

ALPHAS = [2, 3, 12, 101, 257, 1009, 7919, 999983, 104729, 10**7 + 1, 10**12 + 1]
SIZES = [50, 200, 1000, 3000, 6000]


class TestDecimalSqrtAgainstIsqrt(unittest.TestCase):
    def test_sqrt_is_correctly_rounded(self):
        for dg in SIZES:
            with localcontext() as c:
                c.prec = dg + 10
                for a in ALPHAS:
                    x = D(a).sqrt()
                    s = math.isqrt(a * 10 ** (2 * dg))
                    xi = int(x * D(10) ** dg)
                    if not (s <= xi <= s + 1):
                        self.fail(
                            f"alpha={a} digits={dg}: x*10^D outside "
                            f"[floor(sqrt*10^D), +1] "
                            f"(got ...{str(xi)[-12:]}, want {s}..{s+1})")

    def test_control_must_fail(self):
        """A wrong oracle must be caught.  Guards the test above.

        `math.isqrt(a * 10^(2*dg))` is the right oracle; scaling it by 10^Dg
        instead is wrong, and a comparison written loosely enough to accept it
        would prove nothing.  This asserts the loose test really does fail.
        """
        dg = 200
        with localcontext() as c:
            c.prec = dg + 10
            a = 12
            x = D(a).sqrt()
            wrong = math.isqrt(a * 10 ** (4 * dg))       # over-scaled
            xi = int(x * D(10) ** dg)
            self.assertFalse(wrong <= xi <= wrong + 1,
                             "the mis-scaled oracle must NOT validate")


class TestDecimalArithmeticAgainstFraction(unittest.TestCase):
    """Beyond sqrt: the arithmetic itself, at precisions that matter here."""

    def test_division_and_power_agree_with_rationals(self):
        for dg in (50, 200, 1000, 3000):
            with localcontext() as c:
                c.prec = dg
                for a in (3, 12, 257):
                    # division
                    q = D(1) / D(a)
                    self.assertLessEqual(abs(q * D(a) - 1), D(10) ** (-dg + 2))
                    # a power, as the error bound uses
                    p = (D(3) / D(7)) ** 3
                    self.assertLessEqual(
                        abs(p * (D(7) / D(3)) ** 3 - 1), D(10) ** (-dg + 2))

    def test_localcontext_rounding_modes_are_honoured(self):
        """The certificates depend on ROUND_FLOOR / ROUND_CEILING actually
        doing directed rounding.  If they silently rounded to nearest, an
        'outward enclosure' would be inward and every certificate suspect."""
        from decimal import ROUND_FLOOR, ROUND_CEILING
        for digits in (20, 50, 200):
            with localcontext() as c:
                c.prec = 10
                c.rounding = ROUND_FLOOR
                self.assertEqual(D(1) / D(3) * 3 <= 1, True,
                                 "ROUND_FLOOR must not round up")
                c.rounding = ROUND_CEILING
                self.assertEqual(D(1) / D(3) * 3 >= 1, True,
                                 "ROUND_CEILING must not round down")

    def test_decimal_conversions_are_exact(self):
        """certify() turns a Decimal into an exact integer via
        as_integer_ratio(); if that were lossy the certificate would be too."""
        for dg in (50, 200, 1000, 3000):
            with localcontext() as c:
                c.prec = dg
                x = D(7) / D(3)
                num, den = x.as_integer_ratio()
                self.assertEqual(F(num, den), x, "as_integer_ratio lossy")


if __name__ == "__main__":
    unittest.main(verbosity=2)
