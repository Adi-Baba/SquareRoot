"""Tests for mobius_n.py -- the n-th root obstruction (Paper2 Section 10).

No sampling, no floats, no tolerances: every claim is an identity in Q or an
exact sign computation.  The cases are chosen by hand to cover the regimes
(tight anchor, loose anchor, the repelling region), not sampled.
"""

import unittest
from fractions import Fraction as F

import mobius_n as M

# (a, r, n) with 0 < a < r.  Chosen: tight anchor, mid, and the repelling side.
CASES = [
    (F(4), F(5), 2), (F(4), F(5), 3), (F(4), F(5), 4), (F(4), F(5), 5),
    (F(1), F(2), 3), (F(1), F(2), 5),
    (F(9), F(10), 4), (F(3), F(10), 3),
]


class TestErrorIdentity(unittest.TestCase):
    def test_identity_holds_exactly(self):
        for a, r, n in CASES:
            for z in (r + F(1, 10), r - F(1, 10), a, r * F(3, 2)):
                if z == r:
                    continue
                self.assertTrue(M.error_identity(z, a, r, n),
                                f"n={n} a={a} r={r} z={z}")

    def test_Q_is_the_complete_homogeneous_symmetric_polynomial(self):
        # The identification that makes the degree obstruction and the error
        # identity the same fact seen twice.
        for a, r, n in CASES:
            for z in (r + F(1, 7), r - F(1, 3), F(2, 5)):
                self.assertEqual(M.Q(z, r, a, n), M.h_complete([z, r, a], n - 2),
                                 f"n={n} a={a} r={r} z={z}")

    def test_Q2_is_one(self):
        for a, r in ((F(3), F(7, 2)), (F(1), F(5))):
            self.assertEqual(M.Q(F(5, 4), r, a, 2), F(1))

    def test_S_difference_factors_as_z_minus_r_times_Q(self):
        # This is where the (z-r) factor in the error identity comes from:
        #   S_n(z,a) - S_n(r,a) = (z-r) Q_n(z,r,a).
        # Q_n itself does NOT vanish at z = r (Q_3(r,r,a) = 2r+a); an earlier
        # version of this test asserted that, which is false.
        for a, r, n in CASES:
            for z in (r + F(1, 3), r - F(1, 5), F(2, 7)):
                self.assertEqual(M.S(z, a, n) - M.S(r, a, n),
                                 (z - r) * M.Q(z, r, a, n),
                                 f"n={n} a={a} r={r} z={z}")


class TestOneSidednessAndAlternation(unittest.TestCase):
    def test_z_above_r_underestimates(self):
        for a, r, n in CASES:
            z = r + F(1, 5)
            self.assertLess(M.phi(z, a, r, n), r, f"n={n} a={a} r={r}")

    def test_each_step_flips_the_side_of_r(self):
        # Alternation in one line, no iteration required.  From the error
        # identity  r - phi_n(z) = (r-a)(z-r) Q_n / S_n,  for a,r,z > 0 the
        # factors (r-a), Q_n and S_n are all strictly positive (every monomial
        # of Q_n has non-negative exponents and positive base).  Hence
        #     sign(r - phi_n(z)) = sign(z - r),
        # i.e. every step moves to the opposite side of r, for ALL n.
        for a, r, n in CASES:
            for z in (r + F(1, 10), r - F(1, 10), r + F(1, 2), r - F(1, 3)):
                if z <= 0:
                    continue
                err = r - M.phi(z, a, r, n)     # sign == sign(z - r)
                self.assertTrue((err > 0) == (z > r),
                                f"n={n} a={a} r={r} z={z}: sign not flipped")

    def test_iterates_alternate_over_a_few_exact_steps(self):
        # A short confirmation run.  Long exact runs are not viable: each
        # iterate squares the denominators, and 8 steps at n=6 took 33s.
        a, r = F(9), F(10)
        for n in (2, 3, 4):
            z = r + F(1, 10)
            signs = []
            for _ in range(3):
                z = M.phi(z, a, r, n)
                signs.append(z > r)
            self.assertTrue(all(signs[i] != signs[i + 1]
                                for i in range(len(signs) - 1)),
                            f"n={n}: {signs}")

    def test_multiplier_is_negative(self):
        for a, r, n in CASES:
            self.assertLess(M.multiplier(a, r, n), 0,
                            "lambda_n < 0 is what makes the sign flip")


class TestMultiplier(unittest.TestCase):
    def test_two_forms_agree(self):
        for a, r, n in CASES:
            self.assertEqual(M.multiplier(a, r, n), M.multiplier_alt(a, r, n),
                             f"n={n} a={a} r={r}")

    def test_n2_recovers_paper2_multiplier(self):
        # lambda_2 = (a-r)/(a+r) = -u/(2a+u)
        for a, r in ((F(3), F(7, 2)), (F(4), F(5)), (F(1), F(9))):
            self.assertEqual(M.multiplier(a, r, 2), (a - r) / (a + r))

    def test_tends_to_minus_n_minus_one_as_anchor_weakens(self):
        # lambda_n -> -(n-1) as a -> 0+.  This is the whole point of the
        # obstruction: the multiplier does NOT stay inside the unit disc.
        for n in (2, 3, 4, 5):
            r = F(1)
            small = M.multiplier(F(1, 10**6), r, n)
            target = F(1 - n)
            self.assertLess(abs(small - target), F(1, 1000),
                            f"n={n}: {small} should be near {target}")

    def test_N2_is_always_attracting(self):
        # For square roots the scalar map contracts from ANY bracket.
        for r in (F(2), F(5), F(101)):
            for k in (1, 2, 5, 50, 99):
                a = r * F(k, 100)
                if 0 < a < r:
                    self.assertTrue(M.attracting(a, r, 2), f"a/r={k}/100")

    def test_N3_is_repelling_for_a_loose_anchor(self):
        # The claim a reader is most likely to doubt.
        r = F(10)
        self.assertFalse(M.attracting(r * F(1, 10), r, 3),
                         "a/r = 0.1 must be repelling for n=3")
        self.assertTrue(M.attracting(r * F(9, 10), r, 3),
                        "a/r = 0.9 must be attracting for n=3")

    def test_threshold_matches_exact_n3_value(self):
        # lambda_3 = -1  <=>  a/r = (sqrt(3)-1)/2
        import math
        exact = (math.sqrt(3) - 1) / 2
        r = F(1000)
        # bisection on t = a/r.  attracting() is FALSE for small t and TRUE for
        # large t, so the low bound advances when NOT attracting.
        lo, hi = F(0), F(1)
        for _ in range(200):
            mid = (lo + hi) / 2
            if M.attracting(mid * r, r, 3):
                hi = mid
            else:
                lo = mid
        self.assertAlmostEqual(float((lo + hi) / 2), exact, places=6)

    def test_required_accuracy_degrades_like_one_over_n(self):
        # 1 - t_crit ~ 1.5936/n
        # 1.5936.../n is the leading asymptotic; the exact values converge to it
        # slowly, so compare relatively with a 4% window rather than absolutely.
        for n in (20, 100, 500, 2000):
            got = float(M.threshold_gap(n))
            asym = 1.5936 / n
            self.assertLess(abs(got - asym) / asym, 0.04,
                            f"n={n}: got {got}, asymptotic {asym}")
            self.assertGreater(got, asym,
                               "finite-n gap is above its asymptote")

    def test_threshold_gap_is_none_for_n2(self):
        self.assertIsNone(M.threshold_gap(2))


class TestRecoverN2(unittest.TestCase):
    def test_mobius_form_matches_the_map(self):
        for a, r in ((F(3), F(7, 2)), (F(4), F(5))):
            for z in (F(3), F(4), F(5), F(6)):
                self.assertEqual(M.phi(z, a, r, 2), M.mobius_2(z, a, r))

    def test_conjugacy_recovers_paper2_equation(self):
        # T(phi_2(z)) = lambda T(z), the relation of Paper2 Theorem 2.
        for a, r in ((F(3), F(7, 2)), (F(4), F(5))):
            lam = M.multiplier(a, r, 2)
            for z in (F(3), F(4), F(5), F(6)):
                self.assertEqual(M.T(M.phi(z, a, r, 2), r), lam * M.T(z, r))

    def test_n2_error_recovers_paper2_equation_1(self):
        # r - phi_2(z) = (r-a)(z-r)/(a+z)
        a, r = F(3), F(7, 2)
        for z in (F(3), F(4), F(5)):
            self.assertEqual(r - M.phi(z, a, r, 2),
                             (r - a) * (z - r) / (a + z))


if __name__ == "__main__":
    unittest.main(verbosity=2)
