"""Regression tests for exp_mirror.py.

Each class pins a bug that was actually found while building the mirror
experiment, so a future edit that reintroduces it fails loudly instead of
quietly hanging or silently discarding the root.
"""

import unittest
from decimal import Decimal as D, localcontext
from fractions import Fraction as F

import exp_mirror as E

SQUARES = [1, 4, 9, 16, 25, 100, 121, 10000, 1048576, 99980001, 10 ** 14]
NONSQUARES = [2, 3, 10, 12, 999983, 99998016]   # 99998016 is NOT square: 9999^2=99980001
MODES = ("plain", "mirror", "narrow")
DIGITS = 60


def _ok(est, alpha, digits=DIGITS):
    with localcontext() as ctx:
        ctx.prec = digits + 40
        return abs(est - D(alpha).sqrt()) < D(10) ** (-digits)


class TestPerfectSquaresTerminate(unittest.TestCase):
    """alpha a perfect square => alpha/a == a => mirror bracket width 0.

    Then p = 0, est = a, eps = 0, and Decimal(0).adjusted() is 0, which is
    never <= -digits, so the old loop spun forever.  plain mode was immune
    because it starts at w = 1.
    """

    def test_all_mirror_modes_return_immediately(self):
        for alpha in SQUARES:
            for mode in MODES:
                with self.subTest(alpha=alpha, mode=mode):
                    est, _depth = E.zoom_mirror(alpha, DIGITS, mode)
                    self.assertTrue(_ok(est, alpha), f"{mode} alpha={alpha}")

    def test_freehalf_returns_immediately(self):
        for alpha in SQUARES:
            with self.subTest(alpha=alpha):
                est, _depth = E.zoom_freehalf(alpha, DIGITS)
                self.assertTrue(_ok(est, alpha), f"freehalf alpha={alpha}")

    def test_exact_root_is_returned_unchanged(self):
        for alpha in SQUARES:
            with self.subTest(alpha=alpha):
                est, depth = E.zoom_mirror(alpha, DIGITS, "mirror")
                with localcontext() as ctx:
                    ctx.prec = 200
                    self.assertEqual(+est, D(alpha).sqrt())
                self.assertEqual(depth, 0)


class TestSeedingProducesGenuineBrackets(unittest.TestCase):
    """Every seeding mode must return a bracket that actually contains r."""

    def test_mirror_bracket_contains_root(self):
        for alpha in NONSQUARES + SQUARES:
            with localcontext() as ctx:
                ctx.prec = 300
                A = D(alpha)
                a = E._seed(alpha)
                if a * a == A:
                    continue
                lo, hi = a, A / a
                r = A.sqrt()
                self.assertLessEqual(lo, r, alpha)
                self.assertLessEqual(r, hi, alpha)

    def test_freehalf_bracket_contains_root(self):
        # The free half is r = GM(a,b) < AM(a,b), so r lies in [a, AM].
        for alpha in NONSQUARES + SQUARES:
            with localcontext() as ctx:
                ctx.prec = 300
                A = D(alpha)
                a = E._seed(alpha)
                if a * a == A:
                    continue
                AM = (a + A / a) / 2
                r = A.sqrt()
                self.assertLess(r, AM, alpha)          # strict: GM < AM
                self.assertLessEqual(a, r, alpha)

    def test_narrow_picks_a_bracket_that_contains_root(self):
        for alpha in NONSQUARES + SQUARES:
            with localcontext() as ctx:
                ctx.prec = 300
                A = D(alpha)
                a = E._seed(alpha)
                if a * a == A:
                    continue
                w = min(A / a - a, D(1))
                r = A.sqrt()
                self.assertLessEqual(a, r, alpha)
                self.assertLessEqual(r, a + w, alpha)


class TestPerIterationHalvingIsUnsound(unittest.TestCase):
    """The original freehalf bug: halving every iteration, not just at seed.

    Halving keeps [a, a + w/2], which is valid only while r <= a + w/2, i.e.
    only while est - r >= eps/2.  The bound gives est - r <= eps but allows it
    to be arbitrarily close to 0 when u -> 0, so the left half can miss r.
    Checked in exact rationals -- no rounding decides this.
    """

    def _state(self, a, w, u):
        r = a + u
        alpha = r * r
        p = alpha - a * a
        m = a + p / (2 * a + w)
        est = a + p / (a + m)
        eps = w ** 3 / (27 * a * a)
        return r, est, eps, est - r

    def test_root_can_sit_in_the_right_half(self):
        a, w = F(3), F(1)
        u = F(1, 100)                    # r very close to the right end
        r, est, eps, err = self._state(a, w, u)
        self.assertLess(err, eps / 2)     # r > midpoint: left half is wrong
        a_new, w_new = est - eps, eps
        self.assertGreater(r, a_new + w_new / 2)
        self.assertFalse(a_new <= r <= a_new + w_new / 2)

    def test_root_can_sit_in_the_left_half(self):
        a, w = F(3), F(1)
        u = F(2, 3) * w                   # the maximiser: err saturates eps
        r, est, eps, err = self._state(a, w, u)
        self.assertGreater(err, eps / 2)  # here the left half WOULD contain r
        a_new, w_new = est - eps, eps
        self.assertTrue(a_new <= r <= a_new + w_new / 2)

    def test_freehalf_only_halves_at_seeding(self):
        # If the halving were inside the loop the depth would drop below the
        # plain seeding's depth by more than one level; taking it once must
        # save at most roughly a single level.
        with localcontext() as ctx:
            ctx.prec = 400
            for alpha in NONSQUARES:
                plain = E.zoom_mirror(alpha, DIGITS, "plain")[1]
                free = E.zoom_freehalf(alpha, DIGITS)[1]
                self.assertGreaterEqual(free, plain - 1, alpha)
                self.assertLessEqual(free, plain, alpha)


class TestAllModesAreAccurate(unittest.TestCase):
    def test_squares_and_nonsquares(self):
        for alpha in SQUARES + NONSQUARES:
            with localcontext() as ctx:
                ctx.prec = 400
                for mode in MODES:
                    with self.subTest(alpha=alpha, mode=mode):
                        est, _ = E.zoom_mirror(alpha, DIGITS, mode)
                        self.assertTrue(_ok(est, alpha), f"{mode} {alpha}")
                with self.subTest(alpha=alpha, kind="freehalf"):
                    est, _ = E.zoom_freehalf(alpha, DIGITS)
                    self.assertTrue(_ok(est, alpha), alpha)

    def test_mirror_plus_average_is_newton(self):
        # The punchline of the experiment: (x + alpha/x)/2 is Newton's step,
        # so the mirror idea, made rigorous, is the competitor itself.
        with localcontext() as ctx:
            ctx.prec = 60
            alpha, x = 12, D(3)
            for _ in range(6):
                x = (x + D(alpha) / x) / 2
            self.assertTrue(abs(x - D(alpha).sqrt()) < D(10) ** -40)


if __name__ == "__main__":
    unittest.main(verbosity=2)
