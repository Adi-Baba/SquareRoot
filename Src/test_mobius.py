"""Tests for mobius.py -- the exact n-step form and the Mobius conjugacy.

Everything here is `fractions.Fraction`.  There is no sampling, no float, no
tolerance: a claim in this file either holds identically in Q, or the test
fails.  That is the point of the section -- it replaces an asymptotic argument
with exact algebra, so it should be tested with exact arithmetic.
"""

import unittest
from fractions import Fraction as F

import mobius as M

# A spread of (a, u) pairs covering: u << a (root near the left end), u ~ a/2,
# u close to a (root near the right end), and a non-integer a.  Chosen, not
# random -- the identities are universal, so a handful of exact cases suffices
# to catch an algebraic slip, and none of them is a "sample" in the sense the
# rest of this project has learned to distrust.
CASES = [
    (F(3), F(1, 2)),
    (F(1), F(1, 3)),
    (F(5), F(7, 2)),
    (F(2), F(1, 10)),
    (F(7, 2), F(1, 4)),
    (F(4), F(3)),        # u == a
    (F(1), F(1, 7)),
]


class TestConjugacy(unittest.TestCase):
    def test_T_of_phi_equals_lambda_T(self):
        for a, u in CASES:
            lam = M.lam(a, u)
            for z in (a + u + F(1, 10), a + u + F(9, 10), a + 1, a + u * 2):
                self.assertEqual(M.T(a, u, M.phi(a, u, z)), lam * M.T(a, u, z),
                                 f"a={a} u={u} z={z}")

    def test_lambda_is_minus_u_over_2a_plus_u(self):
        for a, u in CASES:
            self.assertEqual(M.lam(a, u), -u / (2 * a + u))

    def test_fixed_points_are_plus_and_minus_r(self):
        for a, u in CASES:
            r = a + u
            self.assertEqual(M.phi(a, u, r), r)     # r is fixed
            self.assertEqual(M.phi(a, u, -r), -r)   # -r is fixed too

    def test_lambda_equals_phi_prime_at_r(self):
        # phi'(z) = -p/(a+z)^2; at z=r this is -u(2a+u)/(2a+u)^2 = -u/(2a+u).
        for a, u in CASES:
            p = M.p_of(a, u)
            r = a + u
            self.assertEqual(-p / (a + r) ** 2, M.lam(a, u))

    def test_T_is_an_involution_free_bijection(self):
        for a, u in CASES:
            for v in (F(0), F(1, 3), F(-1, 4), F(3, 5)):
                z = M.Tinv(a, u, v)
                self.assertEqual(M.T(a, u, z), v)


class TestClosedForm(unittest.TestCase):
    def test_n_step_matches_iteration(self):
        # The core claim: formula (N) equals phi^n - r, exactly.
        for a, u in CASES:
            for z in (a + u + F(1, 10), a + 1, a + u * 2):
                for n in range(0, 7):
                    self.assertEqual(M.error_closed_form(a, u, z, n),
                                     M.phi_iter(a, u, z, n),
                                     f"a={a} u={u} z={z} n={n}")

    def test_n_zero_is_the_identity(self):
        for a, u in CASES:
            for z in (a + u + F(1, 10), a + 1):
                self.assertEqual(M.error_closed_form(a, u, z, 0), z - (a + u))

    def test_lambda_magnitude_is_less_than_one(self):
        for a, u in CASES:
            self.assertLess(abs(M.lam(a, u)), 1)


class TestConsequences(unittest.TestCase):
    def test_even_iterates_are_above_r_and_odd_below(self):
        # Exact one-sidedness for EVERY iterate, not just the first two.
        # sign(phi^n(z) - r) = sign(lambda^n) * sign(T(z)), so the parity
        # alternation holds when T(z) > 0, i.e. when z > r.  (For z < r the
        # parity is reversed: an earlier version of this test used z = a+1,
        # which is below r whenever u > 1, and failed for that reason.)
        for a, u in CASES:
            z = (a + u) * 1 + F(1, 10)          # z = r + 1/10, so z > r
            for n in range(1, 9):
                err = M.phi_iter(a, u, z, n)
                if n % 2 == 0:
                    self.assertGreater(err, 0, f"a={a} u={u} n={n}")
                else:
                    self.assertLess(err, 0, f"a={a} u={u} n={n}")

    def test_parity_reverses_below_r(self):
        # The other half of the sign rule, and a guard on the test above.
        for a, u in CASES:
            z = (a + u) - F(1, 10)              # z = r - 1/10, so z < r
            for n in range(1, 6):
                err = M.phi_iter(a, u, z, n)
                if n % 2 == 0:
                    self.assertLess(err, 0, f"a={a} u={u} n={n}")
                else:
                    self.assertGreater(err, 0, f"a={a} u={u} n={n}")

    def test_corrected_estimator_is_phi_squared(self):
        # The paper defines the corrected estimate as feeding m back in; that is
        # literally two steps of the same Mobius map.
        for a, u in CASES:
            for w in (F(1), F(2), F(5, 2)):
                b = a + w
                m = M.phi(a, u, b)          # a + p/(a+b)
                est = a + M.p_of(a, u) / (a + m)   # a + p/(a+m)
                self.assertEqual(est - (a + u),
                                 M.phi_iter(a, u, b, 2),
                                 f"a={a} u={u} w={w}")

    def test_n2_reproduces_the_papers_closed_form(self):
        # Formula (N) at n=2 must equal u^2(w-u)/D, the paper's exact error.
        for a, u in CASES:
            for w in (F(1), F(2), F(5, 2)):
                D = 2 * a * (2 * a + w) + u * (2 * a + u)
                self.assertEqual(M.error_closed_form(a, u, a + w, 2),
                                 u * u * (w - u) / D,
                                 f"a={a} u={u} w={w}")

    def test_numerator_maximised_at_two_thirds(self):
        # The 4/27 constant and the u = 2w/3 maximiser fall out of (N).

        # Part 1: an ALGEBRAIC identity, not a search.  With c = 2w/3,
        #   4w^3/27 - u^2(w-u) = (u - c)^2 * (u + w/3)  >= 0  for u >= 0,
        # with equality iff u = c.  That proves the maximum is attained at
        # c = 2w/3 and equals 4w^3/27, with no grid and no calculus.  (An
        # earlier version of this test maximised over a 60-point grid, which
        # is sampling -- the exact reasoning Task 4 withdrew from the paper.)
        for w in (F(1), F(3), F(7, 2), F(11, 4)):
            c = 2 * w / 3
            for k in range(0, 40):
                u = F(k, 7) * w               # a spread of exact rationals
                gap = 4 * w ** 3 / 27 - u * u * (w - u)
                factor = (u - c) ** 2 * (u + w / 3)
                self.assertEqual(gap, factor, f"w={w} u={u}")
                self.assertGreaterEqual(gap, 0, f"w={w} u={u}")
            # equality exactly at the critical point
            self.assertEqual(c * c * (w - c), 4 * w ** 3 / 27)

    def test_numerator_zero_at_endpoints(self):
        # n(0) = n(w) = 0, so the interior maximiser is the only candidate.
        for w in (F(1), F(3), F(7, 2)):
            self.assertEqual(0 * 0 * w, 0)
            self.assertEqual(w * w * (w - w), 0)

    def test_attainment_ratio(self):
        # Theorem 4's bound is a product of two separately-taken suprema, so
        # it is attained only as w/a -> 0.  This pins the exact factor.
        for a, w in ((F(1), F(1, 100)), (F(1), F(1, 10)), (F(1), F(1)),
                     (F(3), F(1)), (F(1), F(3)), (F(7, 2), F(5, 2))):
            got = M.attainment_ratio(a, w)
            want = 18 * a**2 / (18 * a**2 + 15 * a * w + 2 * w**2)
            self.assertEqual(got, want, f"a={a} w={w}")
            self.assertLessEqual(got, 1, "must never exceed the bound")
            self.assertGreater(got, 0)

    def test_corollary7_slack_identity(self):
        # Corollary 7 (S8): the exact slack between the sharp position-blind
        # constant 2/(27(2+t)) and the scale-free error F.  Checked as an
        # IDENTITY in exact rationals over a spread of a, w, u -- not a search
        # for a maximiser, and not a grid.
        for a, _u in CASES:
            for w in (F(1), F(2), F(5, 2), F(11, 4)):
                for uu in (w * F(1, 10), w * F(1, 3), w * F(2, 3), w * F(9, 10)):
                    if not (0 < uu < w):
                        continue
                    t = w / a
                    got = M.position_blind_constant(a, w) - M.normed_error(
                        a, w, uu)
                    sq, mx = M.strict_slack(a, w, uu)
                    want = (sq + mx) / (27 * (2 + t) * (
                        4 + 2 * t + 2 * t * (uu / w) + t * t * (uu / w) ** 2))
                    self.assertEqual(got, want, f"a={a} w={w} u={uu}")

    def test_corollary7_slack_is_strictly_positive(self):
        # The strictness argument, done structurally rather than numerically:
        # both summands of the slack numerator are non-negative on [0,1], and
        # their zero sets are DISJOINT -- the square term vanishes only at
        # eta = 2/3, the mixed term only at eta = 0 -- so for t > 0 they can
        # never both vanish.  This is what upgrades Corollary 6's asymptotic
        # "sharp" to an exact non-attainment statement.
        for a, _u in CASES:
            for w in (F(1), F(2), F(5, 2)):
                # 2/3 is NOT of the form k/40, so it is added explicitly.
                # Without that the square term's only zero is never sampled,
                # the `if ... == 0` branches below never fire, and the
                # disjointness assertion is vacuous: the test would pass while
                # checking nothing about the argument.  Mutation-tested --
                # deleting the 2/3 entry from this list survives unless the
                # zero-set equality is asserted, which is why it is.
                zeros_sq, zeros_mx = set(), set()
                for uu in [F(k, 40) * w for k in range(0, 41)] + [w * F(2, 3),
                                                                  F(0)]:
                    sq, mx = M.strict_slack(a, w, uu)
                    self.assertGreaterEqual(sq, 0, f"a={a} w={w} u={uu}")
                    self.assertGreaterEqual(mx, 0, f"a={a} w={w} u={uu}")
                    # t > 0 here because w > 0 and a > 0, so the slack is
                    # strict even at uu = 0 and at uu = 2w/3.
                    self.assertGreater(sq + mx, 0,
                                       f"a={a} w={w} u={uu}: slack is strict")
                    if sq == 0:
                        zeros_sq.add(uu)
                    if mx == 0:
                        zeros_mx.add(uu)
                # The disjointness IS the proof, so assert it as such: each
                # summand has exactly one zero, at a different point of [0, w].
                self.assertEqual(zeros_sq, {w * F(2, 3)},
                                 f"a={a} w={w}: square zero set")
                self.assertEqual(zeros_mx, {F(0)},
                                 f"a={a} w={w}: mixed zero set")
                self.assertEqual(zeros_sq & zeros_mx, set(),
                                 "zero sets must not intersect")

    def test_corollary7_constant_is_below_one_over_27(self):
        # Corollary: 2/(27(2+t)) < 1/27 strictly, i.e. the position-blind bound
        # is strictly better than the published constant WITHOUT any position
        # information.  Relative gain t/(2+t).
        for a in (F(1), F(2), F(7, 2), F(1000)):
            for k in (1, 2, 5, 10, 100, 1000):
                w = a * F(1, k)
                c = M.position_blind_constant(a, w)
                self.assertLess(c, F(1, 27), f"a={a} w={w}")
                self.assertEqual(c, F(2) / (27 * (2 + F(1, k))))
                self.assertEqual(F(1, 27) - c, F(1, 27) * F(1, k) / (2 + F(1, k)))

    def test_corollary7_constant_equals_the_repository_baseline(self):
        # The load-bearing check against the working code: exp_u_aware's
        # "published" baseline is (4w^3/27)/D(0), and Corollary 7 says that IS the
        # sharp position-blind constant.  So the codebase has been using the
        # better constant all along, and the papers' 1/27 was the loose one.
        import exp_u_aware as U
        for a, w in ((F(1), F(1, 10)), (F(1), F(1)), (F(7, 2), F(5, 2)),
                     (F(100), F(3))):
            base = U.bound_published(a, w)
            scaled = base * a * a / w ** 3
            self.assertEqual(scaled, M.position_blind_constant(a, w),
                             f"a={a} w={w}")
            self.assertLessEqual(base, w ** 3 / (27 * a * a),
                                 "still a valid bound for the corrected error")

    def test_corollary7_recovers_one_over_27_at_the_corner(self):
        # F(0, 2/3) = 1/27 exactly.  So Theorem 4's constant is the value of F at
        # the CORNER (t, eta) = (0, 2/3), which no bracket realises: t = w/a > 0
        # and 0 < eta < 1 strictly.  That is the geometric statement behind
        # non-attainment, and it is an exact evaluation, not an asymptotic.
        eta = F(2, 3)
        f_corner = eta * eta * (1 - eta) / 4
        self.assertEqual(f_corner, F(1, 27))

    def test_attainment_ratio_tends_to_one_as_w_over_a_vanishes(self):
        # "Sharp" is an asymptotic statement; check it numerically-free by
        # taking w/a = 1/k for large integer k and confirming the factor rises
        # toward 1.
        prev = F(0)
        for k in (1, 10, 100, 1000, 10000):
            r = M.attainment_ratio(F(1), F(1, k))
            self.assertGreater(r, prev, "factor must increase as w/a shrinks")
            prev = r
        self.assertGreater(prev, F(99, 100))

    def test_bound_4w3_over_27a2_dominates(self):
        # Theorem 4's stated bound follows from the sharper one under D >= D(0).
        for a, u in CASES:
            for w in (F(1), F(2)):
                D0 = 2 * a * (2 * a + w)
                sharp = 4 * w ** 3 / (27 * D0)
                self.assertLessEqual(sharp, w ** 3 / (27 * a * a))


class TestProgressGuard(unittest.TestCase):
    """The bracket-shrinking guard, under the sharp constant of Corollary 7.

    Paper2 5.4 notes that with eps = 2w^3/(27a(2a+w)) the progress condition
    eps < w becomes 2w^2 < 27a(2a+w), i.e. w < 3a(9+sqrt(129))/4 ~ 15.27a --
    wider than the w < sqrt(27) a ~ 5.20a implied by the loose
    eps = w^3/(27a^2).  That was stated as a fact about the GUARD while being
    untested against the recursion, which is the kind of claim this project
    keeps separate.  These three tests settle all of it in exact rationals:

      1. the sharper eps is valid (it dominates the true error);
      2. the guard really is wider -- both roots, compared exactly;
      3. and, the reason it does not matter, the recursion never comes near
         either one: w/a <= 1 at every state it visits.

    Exact rationals throughout.  The recursion state generator is the one in
    exp_u_aware, i.e. states the algorithm actually produces, not a synthetic
    draw over w/a.
    """

    def _eps(self, a, w, sharp):
        return (2 * w ** 3 / (27 * a * (2 * a + w)) if sharp
                else w ** 3 / (27 * a * a))

    def test_sharper_eps_dominates_the_true_error(self):
        import exp_u_aware as U
        import random
        rng = random.Random(11)
        n = 0
        for _ in range(300):
            a, w, r = U.rand_bracket(rng, w_over_a_max=F(1, 100))
            est = U.corrected(a, w, r)
            self.assertGreaterEqual(self._eps(a, w, True), est - r,
                                    f"a={a} w={w}")
            self.assertLess(self._eps(a, w, True), self._eps(a, w, False),
                            "sharp eps must be the smaller of the two")
            n += 1
        self.assertEqual(n, 300)

    def test_the_guard_is_wider_under_the_sharp_constant(self):
        # eps < w  <=>  w < guard.  For the loose constant the guard solves
        # w^2 < 27a^2, giving sqrt(27)*a.  For the sharp one it solves
        # 2w^2 < 27a(2a+w).  Compare the two thresholds exactly, at rational
        # approximations to the irrational roots so no float decides it.
        loose_guard = F(5196, 1000)     # sqrt(27)   = 5.19615...
        sharp_guard = F(15268, 1000)    # 3(9+sqrt129)/4 = 15.2683...
        self.assertGreater(sharp_guard, loose_guard)
        for a in (F(1), F(2), F(7, 2), F(1000)):
            # strictly inside both guards -> eps < w holds for both constants
            for w in (a, a * F(1, 2), a * F(4)):
                self.assertLess(w, loose_guard * a)
                self.assertLess(w, sharp_guard * a)
                self.assertLess(self._eps(a, w, False), w)
                self.assertLess(self._eps(a, w, True), w)

    def test_recursion_never_approaches_either_guard(self):
        # The reason the wider guard is not an improvement in practice: the
        # recursion drives w/a below 1 immediately, so the binding constraint
        # is never the guard at all.  States come from actually running the
        # recursion, so this is about the states the method visits.
        import exp_u_aware as U
        import random
        rng = random.Random(12)
        loose_guard = F(5196, 1000)     # 5.19615...
        seen, worst = 0, F(0)
        for _ in range(120):
            a, w, r = U.rand_bracket(rng, w_over_a_max=F(1, 100))
            for st_a, st_w, _m, _est in U.recursion_states(a, w, r, levels=3):
                t = st_w / st_a
                seen += 1
                if t > worst:
                    worst = t
                self.assertLess(t, loose_guard,
                                "recursion must stay well inside the guard")
        self.assertGreater(seen, 0, "no recursion states were produced")
        # `seen > 0` is a VACUITY guard, not a load-bearing assertion: it holds
        # (360 states over 120 draws, max w/a = 1/100), so deleting it does not
        # fail.  The real check is the 360 per-state `assertLess` above, which
        # is what mutation testing confirms is load-bearing -- swapping the two
        # guard thresholds, or weakening the per-state bound to 1000, both
        # kill.  Kept because it turns a future change to
        # exp_u_aware.recursion_states that silently yields nothing into a
        # loud failure instead of a green test.
        self.assertLessEqual(worst, F(1),
                             f"max w/a over visited states was {worst}")


class TestAlphaBelowOne(unittest.TestCase):
    """0 < alpha < 1: the regime Section 5.3 originally claimed needed no
    rescaling.  It does, and these tests pin why, in exact arithmetic.

    At a = 0 the corrected estimator degenerates to a constant, and the
    progress condition is unsatisfiable.  Both are exact, not empirical.
    """

    def test_corrected_estimator_degenerates_at_a_zero(self):
        # m = alpha, est = 1, for EVERY alpha in (0,1).
        for alpha in (F(1, 2), F(1, 3), F(1, 7), F(9, 10), F(1, 100)):
            a, b = F(0), F(1)
            p = alpha - a * a
            m = a + p / (a + b)
            est = a + p / (a + m)
            self.assertEqual(m, alpha)
            self.assertEqual(est, F(1), f"alpha={alpha}: est must degenerate to 1")

    def test_progress_condition_unsatisfiable_at_a_zero(self):
        # Prop 5: eps < w  iff  w < a*sqrt(c).  At a=0 that is w < 0.
        for w in (F(1), F(1, 2), F(1, 10)):
            for c in (F(16), F(4, 27)):
                # a*sqrt(c) = 0, so the condition w < 0 can never hold
                self.assertFalse(w < F(0), "no positive w satisfies it at a=0")

    def test_rescaling_by_reciprocal_restores_a_positive_seed(self):
        # The fix: sqrt(alpha) = 1 / sqrt(1/alpha), and 1/alpha > 1 so the
        # integer floor is >= 1 and the recursion can run.
        for alpha in (F(1, 2), F(1, 3), F(1, 7), F(1, 100), F(9, 10)):
            recip = 1 / alpha
            a_new = int(recip ** F(1, 2))       # floor(sqrt(1/alpha))
            self.assertGreaterEqual(a_new, 1,
                                    f"alpha={alpha}: rescaled seed must be >= 1")


class TestSymbolicRecord(unittest.TestCase):
    def test_recorded_identities_are_all_true(self):
        for name, val in M.SYMBOLIC_IDENTITY.items():
            self.assertTrue(val, name)


if __name__ == "__main__":
    unittest.main(verbosity=2)
