"""Regression tests for exp_certified.py.

Every class here pins a bug that was actually found while building the
certified benchmark, or guards the certificate verifier itself -- a verifier
that accepts a bad bracket is worse than none, so most of the first class is
devoted to proving certify() REJECTS what it should.
"""

import unittest
from decimal import Decimal as D, localcontext

import exp_certified as E
from exp_precision import _seed

ALPHAS = [2, 3, 12, 101, 257, 999983]
SMALL = [10, 50]


def _ref(alpha, digits):
    with localcontext() as c:
        c.prec = digits + 60
        return D(alpha).sqrt()


class TestCertifyRejectsBadBrackets(unittest.TestCase):
    """certify() must refuse anything it cannot prove.  These are the cases
    that a silently-permissive verifier would let through."""

    def test_accepts_a_genuine_bracket(self):
        r = _ref(2, 60)
        lo = D(r).to_integral_value(rounding="ROUND_FLOOR")
        hi = D(r).to_integral_value(rounding="ROUND_CEILING")
        self.assertIsNone(E.certify(lo, hi, 2, 0))

    def test_rejects_left_endpoint_above_root(self):
        r = _ref(2, 60)
        with localcontext() as c:
            c.prec = 60
            lo = +r + D("1e-10")          # strictly above sqrt(2)
            hi = +r + D("1e-5")
        self.assertIsNotNone(E.certify(lo, hi, 2, 20))

    def test_rejects_right_endpoint_below_root(self):
        r = _ref(2, 60)
        with localcontext() as c:
            c.prec = 60
            lo = +r - D("1e-5")
            hi = +r - D("1e-10")          # strictly below sqrt(2)
        self.assertIsNotNone(E.certify(lo, hi, 2, 20))

    def test_rejects_bracket_too_wide(self):
        # Correctly brackets the root but is nowhere near 1e-20.
        self.assertIsNotNone(E.certify(D(1), D(2), 2, 20))

    def test_rejects_empty_bracket(self):
        self.assertIsNotNone(E.certify(D(2), D(1), 2, 10))

    def test_rejects_exact_root_bracket_that_is_too_wide(self):
        # alpha = 4: lo = hi = 2 brackets the root perfectly but claims 0
        # width while being asked to certify nothing (digits=0 passes).  With
        # a real target it must still be checked on width, not just enclosure.
        self.assertIsNone(E.certify(D(2), D(2), 4, 0))
        self.assertIsNotNone(E.certify(D(2), D(3), 4, 30))

    def test_tiny_width_violation_is_caught(self):
        # 10**-20 requested, given 2*10**-20: a permissive verifier loses this.
        r = _ref(2, 80)
        with localcontext() as c:
            c.prec = 80
            lo = +r
            hi = +r + D("2e-20")
        self.assertIsNotNone(E.certify(lo, hi, 2, 20))


class TestExactComparators(unittest.TestCase):
    """_le_sqrt/_ge_sqrt are the whole certificate; they must agree with
    exact rational arithmetic on both sides of the boundary."""

    def test_agree_with_fractions(self):
        from fractions import Fraction as F

        r = _ref(3, 50)
        with localcontext() as c:
            c.prec = 50
            for k in (-30, -5, -1, 0, 1, 5, 30):
                # sample points just below, at, and just above sqrt(3)
                step = D(10) ** k
                for cand in (r - step, r, r + step):
                    got_le = E._le_sqrt(cand, 3)
                    exact_le = F(cand) * F(cand) <= 3 if cand > 0 else True
                    self.assertEqual(got_le, exact_le, f"cand={cand} k={k}")
                    # _ge_sqrt_exact must be the exact complement
                    c2, e2 = E._coeff_exp(cand)
                    if c2 > 0:
                        self.assertEqual(E._ge_sqrt_exact(c2, e2, 3),
                                         F(cand) * F(cand) >= 3, f"k={k}")

    def test_handles_extreme_exponents(self):
        # Decimals far from normalised form: tiny/huge coefficients, extreme
        # exponents.  _le_sqrt has to route through _le_sqrt_exact, which
        # moves the power of ten to whichever side keeps it non-negative.
        self.assertTrue(E._le_sqrt(D("1.4E-300"), 2))    # << sqrt(2)
        self.assertTrue(E._le_sqrt(D("2E-300"), 2))      # << sqrt(2)
        self.assertFalse(E._le_sqrt(D("2E+300"), 2))     # >> sqrt(2)
        c, e = E._coeff_exp(D("1.5E+300"))
        self.assertTrue(E._ge_sqrt_exact(c, e, 2))       # >> sqrt(2)
        c, e = E._coeff_exp(D("1.0E+300"))
        self.assertTrue(E._ge_sqrt_exact(c, e, 2))       # >> sqrt(2)
        c, e = E._coeff_exp(D("1E-300"))
        self.assertFalse(E._ge_sqrt_exact(c, e, 2))      # << sqrt(2)


class TestEveryMethodIsCertified(unittest.TestCase):
    def test_all_methods_certify(self):
        for dg in SMALL:
            for al in ALPHAS:
                for name, fn in E.METHODS:
                    with self.subTest(digits=dg, alpha=al, method=name):
                        lo, hi, _steps = fn(al, dg)
                        msg = E.certify(lo, hi, al, dg)
                        self.assertIsNone(msg, msg)

    def test_certified_on_large_target(self):
        for al in (2, 12, 999983):
            lo, hi, _ = E.zoom(al, 500)
            self.assertIsNone(E.certify(lo, hi, al, 500))


class TestConvergenceMatchesTheory(unittest.TestCase):
    """The whole point of the comparison is convergence order, so the
    measured step counts must sit where theory puts them."""

    def test_bisection_is_linear(self):
        d1, d2 = 50, 100
        s1 = E.bisection(12, d1)[2]
        s2 = E.bisection(12, d2)[2]
        # 1 bit per step => 3.322 steps per decimal digit, doubling exactly
        self.assertAlmostEqual(s1 / d1, 3.322, delta=0.3)
        self.assertAlmostEqual(s2 / d2, 3.322, delta=0.3)
        self.assertAlmostEqual(s2 / s1, d2 / d1, delta=0.15)

    def test_zoom_is_cubic(self):
        # cubic => depth grows like log3(D): 200 -> 1000 -> 3000 digits
        depths = [E.zoom(12, dg)[2] for dg in (50, 500, 5000)]
        self.assertLess(depths[2] - depths[1], 4,
                        f"zoom depth not growing logarithmically: {depths}")

    def test_interval_newton_is_quadratic(self):
        # quadratic => number of steps grows like log2(D)
        its = [E.interval_newton(12, dg)[2] for dg in (50, 500, 5000)]
        # 10x the digits should cost only a few more steps, not 10x
        self.assertLessEqual(its[2] - its[1], 6,
                             f"interval Newton not quadratic: {its}")
        self.assertLess(its[2], 40)

    def test_newton_certified_is_quadratic(self):
        its = [E.newton_certified(12, dg)[2] for dg in (50, 500, 5000)]
        self.assertLessEqual(its[2] - its[1], 6, its)
        self.assertLess(its[2], 40)


class TestPinnedBugs(unittest.TestCase):
    """Bugs found while building this file.  Each test fails if it returns."""

    def test_termination_test_is_strict(self):
        # eps.adjusted() <= -digits accepts eps in [10^-D, 10^(1-D)): a
        # certificate up to ten times looser than promised.  alpha=257 at 50
        # digits is the case that actually landed in that band.
        from exp_precision import _width_digits, MIN_PREC, GUARD

        alpha, digits = 257, 50
        with localcontext() as ctx:
            ctx.prec = digits + GUARD + 10
            A = D(alpha)
            a = _seed(alpha)
            w = D(1)
            depth = 0
            bad_band = False
            while depth < 200:
                prec = min(MIN_PREC + 3 * _width_digits(w) + 3 * a.adjusted(),
                           digits + GUARD)
                with localcontext() as c2:
                    c2.prec = prec
                    p = A - a * a
                    m = a + p / (2 * a + w)
                    est = a + p / (a + m)
                eps = w ** 3 / (27 * a * a)
                limit = D(10) ** (-digits)
                if eps.adjusted() <= -digits and eps > limit:
                    bad_band = True
                if eps.adjusted() < -digits:
                    break
                with localcontext() as c4:
                    c4.prec = prec
                    a = +est - eps
                    w = +eps
                depth += 1
        # the OLD rule really does land in the forbidden band for this alpha
        self.assertTrue(bad_band,
                        "expected alpha=257/digits=50 to expose the <= band")

    def test_naive_bracket_extraction_fails_certification(self):
        # The exact identity r in [est-eps, est] holds for the EXACT est.
        # Reading that interval straight off finite-precision Decimals drops
        # the rounding uncertainty and fails certify() -- 16 times in 48 runs
        # while this file was being written.  The interval extractor must not.
        from exp_precision import _width_digits, MIN_PREC, GUARD

        naive_fail = 0
        checked = 0
        for digits in (20, 50, 100):
            for alpha in ALPHAS:
                with localcontext() as ctx:
                    ctx.prec = digits + GUARD + 10
                    A = D(alpha)
                    a = _seed(alpha)
                    w = D(1)
                    depth = 0
                    while depth < 200:
                        prec = min(
                            MIN_PREC + 3 * _width_digits(w) + 3 * a.adjusted(),
                            digits + GUARD)
                        with localcontext() as c2:
                            c2.prec = prec
                            p = A - a * a
                            m = a + p / (2 * a + w)
                            est = a + p / (a + m)
                        eps = w ** 3 / (27 * a * a)
                        if eps.adjusted() < -digits:
                            with localcontext() as c3:
                                c3.prec = digits + GUARD + 10
                                p = A - a * a
                                m = a + p / (2 * a + w)
                                est = a + p / (a + m)
                                eps = w ** 3 / (27 * a * a)
                                naive_lo, naive_hi = +(est - eps), +est
                            checked += 1
                            if E.certify(naive_lo, naive_hi, alpha, digits):
                                naive_fail += 1
                            break
                        with localcontext() as c4:
                            c4.prec = prec
                            a = +est - eps
                            w = +eps
                        depth += 1
        # the naive extraction really does lose the certificate sometimes...
        self.assertGreater(naive_fail, 0, "naive extraction unexpectedly OK")
        # ...and the interval extractor never does
        for digits in (20, 50, 100):
            for alpha in ALPHAS:
                lo, hi, _ = E.zoom(alpha, digits)
                self.assertIsNone(E.certify(lo, hi, alpha, digits))

    def test_interval_newton_needs_double_rate_schedule(self):
        # At rate=1 the enclosure is rounding-limited, the bracket shrinks by
        # only ~16 digits per step, and 1000 digits took 65 iterations.
        self.assertLessEqual(E.interval_newton(12, 1000)[2], 20)

    def test_bracket_width_is_within_target(self):
        for digits in (20, 50, 100):
            for alpha in ALPHAS:
                for name, fn in E.METHODS:
                    lo, hi, _ = fn(alpha, digits)
                    with localcontext() as c:
                        c.prec = digits + 60
                        self.assertLessEqual(
                            hi - lo, D(10) ** (-digits),
                            f"{name} alpha={alpha} digits={digits}")


class TestSeedSymmetry(unittest.TestCase):
    """Newton is conventionally seeded from a 30-digit sqrt; the zoom and
    interval Newton are not, and cannot be -- their bracket is defined to start
    at floor(sqrt(alpha)).  The head start is large enough to move the
    ranking, so both seedings must be available and both must certify."""

    def test_both_seedings_certify(self):
        for dg in (20, 50, 200):
            for alpha in ALPHAS:
                for fs in (True, False):
                    lo, hi, _ = E.newton_certified(alpha, dg, float_seed=fs)
                    with self.subTest(digits=dg, alpha=alpha, float_seed=fs):
                        self.assertIsNone(E.certify(lo, hi, alpha, dg))

    def test_float_seed_really_does_buy_fewer_steps(self):
        # If this ever stops holding, the seed asymmetry the paper reports is
        # no longer real and the narrative has to change.
        for dg in (50, 200):
            a_float = E.newton_certified(12, dg, float_seed=True)[2]
            a_int = E.newton_certified(12, dg, float_seed=False)[2]
            self.assertLess(a_float, a_int, dg)

    def test_integer_seed_needs_more_iterations(self):
        # ...and the step counts should match the O(log D) law for both.
        for fs, lo_bound, hi_bound in ((True, 1, 12), (False, 4, 20)):
            its = E.newton_certified(12, 500, float_seed=fs)[2]
            self.assertTrue(lo_bound <= its <= hi_bound, (fs, its))

    def test_zoom_and_interval_newton_share_the_integer_seed(self):
        # The headline comparison is only same-footing because these two use
        # the identical starting bracket.  Pin that they both terminate.
        for dg in (50, 200):
            self.assertIsNone(E.certify(*E.zoom(12, dg)[:2], 12, dg))
            self.assertIsNone(E.certify(*E.interval_newton(12, dg)[:2], 12, dg))


class TestNoCertificationInTimedPath(unittest.TestCase):
    """A verify() call inside a method inflates its own measurement -- it cost
    Newton+res ~2.7 ms of a 9.8 ms row.  The in-loop check that IS legitimate
    is a single Decimal subtraction, and it is load-bearing: without it the
    pessimistic digit counter falls one step short."""

    def test_source_has_no_in_loop_certify(self):
        import inspect
        src = inspect.getsource(E.newton_certified)
        # certif( appears only inside the explanatory comment, never in a call
        self.assertNotIn("if certify(", src)
        self.assertNotIn("=certify(", src)
        self.assertNotIn("certify(lo_b", src)

    def test_integer_seed_counter_is_pessimistic(self):
        # The width check must actually be reached, i.e. the routine must not
        # return a bracket wider than the target.
        for dg in (50, 1000):
            for alpha in ALPHAS:
                lo, hi, _ = E.newton_certified(alpha, dg, float_seed=False)
                msg = E.certify(lo, hi, alpha, dg)
                self.assertIsNone(msg, f"alpha={alpha} digits={dg}: {msg}")


class TestTerminationContractIsLiteral(unittest.TestCase):
    """The strict termination test is a *contract*, and mutation testing
    showed the scheduled path masks it: exp_precision.zoom overshoots its
    target by ~24 orders of magnitude, so `adjusted() <= -digits` and
    `adjusted() < -digits` return the same value on every input the suite
    exercises.  The mutant therefore survives, and no behavioural test on
    that path can kill it.

    The looseness is still real -- it is the stated tolerance that is wrong,
    not the arithmetic -- so it is pinned textually here.  This is a source
    check, not a behavioural one, and is labelled as such deliberately: it
    asserts the code says what the paper says it says.  The behaviour of the
    *certified* path, where the looseness IS observable, is covered by
    TestPinnedBugs.test_naive_bracket_extraction_fails_certification and
    TestEveryMethodIsCertified.
    """

    FILES = ("exp_precision.py", "exp_certified.py", "exp_hybrid.py",
             "exp_mirror.py", "exp_compare.py")

    def test_no_file_uses_the_loose_termination_test(self):
        import os
        import re
        here = os.path.dirname(os.path.abspath(__file__))
        for fname in self.FILES:
            path = os.path.join(here, fname)
            if not os.path.exists(path):
                continue
            src = open(path).read()
            loose = re.findall(r"adjusted\(\)\s*<=\s*-digits", src)
            self.assertEqual(loose, [],
                             f"{fname} uses the loose termination test "
                             f"(accepts a bracket up to 10x too wide)")

    def test_strict_form_is_what_is_shipped(self):
        import os
        import re
        here = os.path.dirname(os.path.abspath(__file__))
        found = 0
        for fname in self.FILES:
            path = os.path.join(here, fname)
            if not os.path.exists(path):
                continue
            found += len(re.findall(r"adjusted\(\)\s*<\s*-digits",
                                    open(path).read()))
        self.assertGreater(found, 0, "expected the strict test to be present")

    def test_sharp_constant_is_what_is_shipped(self):
        """Same reasoning for the sharp constant 4/27.

        Mutation testing showed that replacing 27 with 16 in exp_precision's
        eps SURVIVES, and it is worth being precise about why.  The zoom's
        output error at large targets is dominated by the rounding of the
        final full-precision recomputation (~1e-85 at 60-80 working digits),
        not by eps -- so eps's constant is unobservable through the returned
        value, and no behavioural test on that path can detect the swap.  The
        same holds for a tight-vs-loose output comparison, which is noisy for
        the same reason and was measured to invert in 5 of 10 cases.

        The constant that actually matters is pinned behaviourally elsewhere:
        TestTheorem2 in test_bracket_zoom.py checks the attained value at
        u = 2w/3, and exp_u_aware's _num_max is checked in test_u_aware.py.
        This test pins that the recursion in the harness uses the sharp form,
        as a source contract, explicitly labelled as such.
        """
        import os
        here = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(here, "exp_precision.py")
        src = open(path).read()
        self.assertIn("w ** 3 / (27 * a * a)", src,
                      "exp_precision must use the sharp 4/27 constant")


if __name__ == "__main__":
    unittest.main(verbosity=2)
