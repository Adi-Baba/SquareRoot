"""Regression tests for exp_u_aware.py and exp_depth.py.

Each class here pins a bug that was actually found, so that a future edit which
reintroduces it fails loudly rather than quietly inflating a headline number.
"""

import random
import unittest
from fractions import Fraction as F

import exp_u_aware as U
from exp_depth import (
    _reached,
    _step,
    bound_published as d_bound_published,
    depth_adaptive,
    depth_published,
    m_le_r_le_est,
    trace_adaptive,
    trace_full,
)
from mpmath import mp, mpf


# ============================================== exp_u_aware: bound validity

class TestBoundValidity(unittest.TestCase):
    """Bug 3: the 'tight' path used to locate the maximiser by ternary search
    and return the value at the search endpoint -- strictly below the true
    maximum, so it certified nothing.  A bound must never under-shoot."""

    def test_bound_dominates_error_at_every_piece_count(self):
        rng = random.Random(7)
        for _ in range(400):
            a, w, r = U.rand_bracket(rng, w_over_a_max=100)
            m, est = U.plain(a, w, r), U.corrected(a, w, r)
            e = est - r
            for pieces in (1, 2, 7, 64):
                b = U.bound_known_interval(a, w, m, est, pieces=pieces)
                self.assertGreaterEqual(
                    b, e,
                    f"bound under-shot the error at pieces={pieces}: "
                    f"{b} < {e} (a={a}, w={w})",
                )

    def test_more_pieces_only_tighten_never_undershoot(self):
        """A refinement must move toward the truth from ABOVE, never below it."""
        rng = random.Random(8)
        for _ in range(200):
            a, w, r = U.rand_bracket(rng, w_over_a_max=100)
            m, est = U.plain(a, w, r), U.corrected(a, w, r)
            prev = None
            for pieces in (1, 2, 4, 8, 16, 64):
                b = U.bound_known_interval(a, w, m, est, pieces=pieces)
                if prev is not None:
                    self.assertLessEqual(
                        b, prev,
                        f"refining pieces {pieces} made the bound worse",
                    )
                prev = b



class TestSameFooting(unittest.TestCase):
    """Bug: the baseline returned bare 4w^3/27 with no denominator while the
    candidate divided by 2a(2a+w), inflating the gain to '100% tighter'."""

    def test_baseline_depends_on_a(self):
        w = F(1)
        self.assertNotEqual(U.bound_published(F(1), w),
                            U.bound_published(F(100), w))

    def test_baseline_is_4w3_over_27_times_denominator(self):
        a, w = F(3), F(1)
        expect = F(4) * w ** 3 / F(27) / (2 * a * (2 * a + w))
        self.assertEqual(U.bound_published(a, w), expect)

    def test_both_sides_call_the_same_routine(self):
        """Same footings means same function, different interval."""
        rng = random.Random(10)
        a, w, r = U.rand_bracket(rng)
        m, est = U.plain(a, w, r), U.corrected(a, w, r)
        # candidate restricted to the full interval must equal the baseline
        self.assertEqual(
            U._bound(a, w, F(0), w, pieces=1),
            U.bound_published(a, w, pieces=1),
        )
        self.assertEqual(
            U._bound(a, w, m - a, est - a, pieces=1),
            U.bound_known_interval(a, w, m, est, pieces=1),
        )

    def test_candidate_never_exceeds_baseline(self):
        rng = random.Random(11)
        for _ in range(300):
            a, w, r = U.rand_bracket(rng, w_over_a_max=100)
            m, est = U.plain(a, w, r), U.corrected(a, w, r)
            self.assertLessEqual(
                U.bound_known_interval(a, w, m, est),
                U.bound_published(a, w),
            )


class TestOptimalCertificate(unittest.TestCase):
    """The position-aware OPTIMAL certificate C(t;L,U) of Paper2 5.4.

    `C(t;L,U) = max_{eta in [L,U]} F(t,eta)` is the supremum of the exact
    error over the information the method possesses, so it is the best any
    position-aware certificate can be.  Theorem 8's B is an upper bound on it,
    and the two coincide iff B is already tight.

    What these tests pin is NOT that C is better -- it is measurably not,
    0.0% in the realistic regime (`exp_u_aware.measure_optimality_gap`) -- but
    that C is a legitimate certificate at all: it must never under-shoot the
    true error, and it must never exceed B.  A routine that under-shot would
    certify nothing, which is the failure mode this project keeps hitting.
    """

    PIECES = (1, 4, 64)

    def test_optimal_is_always_a_valid_bound(self):
        rng = random.Random(21)
        for _ in range(400):
            a, w, r = U.rand_bracket(rng, w_over_a_max=F(1, 100))
            m, est = U.plain(a, w, r), U.corrected(a, w, r)
            e = est - r
            for p in self.PIECES:
                self.assertGreaterEqual(U.bound_optimal(a, w, m, est, pieces=p),
                                        e, f"a={a} w={w} pieces={p}")

    def test_optimal_never_exceeds_theorem8(self):
        # Refining the piece count only ever returns a value of the bound at a
        # real point, and we take a max, so more pieces can only tighten.
        rng = random.Random(22)
        prev_pieces = self.PIECES[0]
        for _ in range(400):
            a, w, r = U.rand_bracket(rng, w_over_a_max=F(1, 100))
            m, est = U.plain(a, w, r), U.corrected(a, w, r)
            b = U.bound_known_interval(a, w, m, est, pieces=prev_pieces)
            prev = None
            for p in self.PIECES:
                c = U.bound_optimal(a, w, m, est, pieces=p)
                self.assertLessEqual(c, b, f"pieces={p} must not loosen")
                if prev is not None:
                    self.assertLessEqual(c, prev, "more pieces must not loosen")
                prev = c

    def test_optimal_reduces_to_theorem8_at_one_piece(self):
        # Same routine and same interval: at pieces=1 they are the same call.
        rng = random.Random(23)
        for _ in range(100):
            a, w, r = U.rand_bracket(rng, w_over_a_max=F(1, 100))
            m, est = U.plain(a, w, r), U.corrected(a, w, r)
            self.assertEqual(U.bound_optimal(a, w, m, est, pieces=1),
                             U.bound_known_interval(a, w, m, est, pieces=1))

    def test_theorem8_is_already_tight_on_the_recursion_interval(self):
        # The reason optimality buys nothing.  The interval [m-a, est-a] IS the
        # current error, so refining 64 pieces finds almost nothing.
        #
        # Note the quantity that matters is the MAGNITUDE of the improvement,
        # not how many states improve.  Measured in exact rationals, most
        # visited states DO improve by a hair -- 95 of 120 -- which reads like
        # a result until you ask by how much: at most 5e-6 relative, median
        # exactly 0.  Asserting on a count would have called this a win; the
        # measurement below is why it is a null result.  (In floats the
        # improvement rounds away entirely, which is the only reason the
        # printed ratio in exp_u_aware reads 1.0000.)
        rng = random.Random(24)
        checked, tighter = 0, 0
        worst = F(0)      # max over states of (B - C) / B, exactly
        for _ in range(60):
            a, w, r = U.rand_bracket(rng, w_over_a_max=F(1, 100))
            for st_a, st_w, m, est in U.recursion_states(a, w, r, levels=2):
                b = U.bound_known_interval(st_a, st_w, m, est)
                c = U.bound_optimal(st_a, st_w, m, est)
                self.assertGreaterEqual(c, est - r)
                checked += 1
                if b > 0:
                    gain = (b - c) / b
                    if gain > worst:
                        worst = gain
                    if c < b:
                        tighter += 1
        self.assertGreater(checked, 0, "no recursion states produced")
        # The load-bearing assertion: optimality recovers under 1 part in 10^4
        # even at its worst on the intervals the method really visits.
        self.assertLess(float(worst), F(1, 10 ** 4),
                        f"worst relative gain {float(worst):.3e} -- if this grew, "
                        "the 'Theorem 8 is already tight' claim needs revisiting")
        # And most states are untouched, so the mean gain is smaller still.
        self.assertLess(float(worst) * tighter / checked, F(1, 100000))


class TestExactIdentities(unittest.TestCase):
    def test_three_identities_hold_exactly(self):
        rng = random.Random(12)
        for _ in range(500):
            a, w, r = U.rand_bracket(rng, w_over_a_max=100)
            m, est = U.plain(a, w, r), U.corrected(a, w, r)
            self.assertEqual(U.plain_error(a, w, r), r - m)
            self.assertEqual(U.corrected_error(a, w, r), est - r)

    def test_known_interval_brackets_u(self):
        rng = random.Random(13)
        for _ in range(500):
            a, w, r = U.rand_bracket(rng, w_over_a_max=100)
            u = r - a
            m, est = U.plain(a, w, r), U.corrected(a, w, r)
            self.assertTrue(m - a < u < est - a)


# ================================================= exp_depth: depth counting

class TestDepthDriver(unittest.TestCase):
    """Bug 1 (baseline constant inverted, `a` ignored) and Bug 2 (adaptive
    side never updated `a`, so level 2+ ran on a bracket without r in it)."""

    def test_both_schemes_share_everything_except_the_bound(self):
        """The experiment isolates ONE variable, so `est` must be identical."""
        with mp.workdps(60):
            alpha, a = mpf(17), mpf(4)
            est_pub, b_pub = _step(a, mpf(1), alpha, u_aware=False)
            est_awa, b_awa = _step(a, mpf(1), alpha, u_aware=True)
            self.assertEqual(est_pub, est_awa, "schemes must compute the same est")
            self.assertNotEqual(b_pub, b_awa, "the bounds must actually differ")

    def test_baseline_constant_is_not_inverted(self):
        """w^3/((4/27)a^2) = 27w^3/(4a^2) > w for a=1 -- the old sign error made
        the bound GROW.  The baseline must shrink as w shrinks."""
        with mp.workdps(60):
            a = mpf(1)
            b1 = d_bound_published(a, mpf(1))
            b2 = d_bound_published(a, mpf(1) / 10)
            self.assertLess(b2, b1)
            self.assertLess(b1, mpf(1), "bound must be < w to make progress")

    def test_baseline_depends_on_a(self):
        with mp.workdps(60):
            self.assertNotEqual(d_bound_published(mpf(1), mpf(1)),
                                d_bound_published(mpf(50), mpf(1)))

    def test_bracket_invariant_holds_through_the_recursion(self):
        """If either bound ever under-shot, r would leave the bracket."""
        with mp.workdps(360):
            for al in (2, 17, 1009):
                alpha, r = mpf(al), mp.sqrt(mpf(al))
                for u_aware in (False, True):
                    a, w = mpf(int(al ** 0.5)), mpf(1)
                    for _ in range(12):
                        est, b = _step(a, w, alpha, u_aware)
                        # stop once est-r would fall under the working
                        # precision: past that the invariant is decided by
                        # rounding noise, not by the bound
                        if b < mpf(10) ** (-(mp.dps - 25)):
                            break
                        self.assertTrue(
                            m_le_r_le_est(est - b, est, r),
                            f"bracket lost r for alpha={al}, u_aware={u_aware}",
                        )
                        a, w = est - b, b

    def test_u_aware_never_needs_more_levels_than_baseline(self):
        # working precision must exceed the target: if it does not, est - b
        # rounds b away, the new left endpoint overshoots r, and the invariant
        # check reports FAILURE.  That is correct behaviour, so the test has to
        # give the recursion enough precision rather than blame the bound.
        for al in (2, 3, 17, 1009, 999983):
          for digits in (50, 500, 3000):
            with mp.workdps(digits + 60):
                alpha, r = mpf(al), mp.sqrt(mpf(al))
                a = mpf(int(al ** 0.5))
                p = depth_published(a, alpha, r, digits)
                q = depth_adaptive(a, alpha, r, digits)
                self.assertNotEqual(p, -1, f"baseline failed {al}/{digits}")
                self.assertNotEqual(q, -1, f"u-aware failed {al}/{digits}")
                self.assertLessEqual(
                    q, p,
                    f"u-aware used more levels for alpha={al}, {digits} digits",
                )

    def test_u_aware_reaches_working_precision_not_failure(self):
        """b == 0 means p underflowed: a is the correctly-rounded root.  An
        earlier version returned -1 there, reporting success as failure."""
        with mp.workdps(200):
            alpha, r = mpf(17), mp.sqrt(mpf(17))
            d = depth_adaptive(mpf(4), alpha, r, 20)
            self.assertGreater(d, 0)
            self.assertNotEqual(d, -1)

    def test_traces_agree_with_depth_counters(self):
        with mp.workdps(400):
            alpha, r = mpf(1009), mp.sqrt(mpf(1009))
            a = mpf(31)
            self.assertEqual(len(trace_full(a, alpha, r, 500)),
                             depth_published(a, alpha, r, 500))
            self.assertEqual(len(trace_adaptive(a, alpha, r, 500)),
                             depth_adaptive(a, alpha, r, 500))


class TestReachedPredicate(unittest.TestCase):
    """_reached replaces a full mp.log10 at 30000-digit working precision,
    which dominated the script.  It must agree with the thing it replaced."""

    def test_agrees_with_mp_log10_away_from_the_threshold(self):
        # Exactly at b = 10^-digits both sides are decided by the last bit of
        # an exp/log round-trip, so the oracle is as unreliable as the
        # predicate; test the unambiguous offsets instead.
        with mp.workdps(80):
            for digits in (10, 100, 1000):
                for k in (-4, -3, -2, -1, 1, 2, 3, 4):
                    b = mpf(10) ** (-digits + k)
                    self.assertEqual(
                        _reached(b, digits),
                        -mp.log10(b) >= digits,
                        f"disagreement at digits={digits}, offset={k}",
                    )

    def test_threshold_sides(self):
        with mp.workdps(80):
            self.assertTrue(_reached(mpf(10) ** -11, 10))   # 11 digits > 10
            self.assertFalse(_reached(mpf(10) ** -9, 10))   # 9 digits < 10

    def test_zero_counts_as_reached(self):
        self.assertTrue(_reached(mpf(0), 100))


class TestNoResurrectedBrokenCode(unittest.TestCase):
    """The two depth functions had an inverted constant and never updated `a`;
    decimal_log10 existed only to serve them.  They were deleted, not fixed.
    If they come back, these fail."""

    def test_fraction_depth_functions_are_gone(self):
        for name in ("depth_published", "depth_adaptive", "decimal_log10"):
            self.assertFalse(
                hasattr(U, name),
                f"exp_u_aware.{name} has come back -- it was wrong; "
                f"use exp_depth instead",
            )

    def test_recursion_states_keeps_the_bracket(self):
        rng = random.Random(14)
        for _ in range(30):
            a = F(rng.randint(1, 10 ** 3))
            w = a * F(rng.randint(1, 1000), 1000)
            r = a + w * F(rng.randint(1, 999), 1000)
            if not (a < r < a + w):
                continue
            states = U.recursion_states(a, w, r, levels=4)
            self.assertTrue(states)
            for st_a, st_w, m, est in states:
                self.assertLessEqual(m, est)


if __name__ == "__main__":
    unittest.main(verbosity=2)
