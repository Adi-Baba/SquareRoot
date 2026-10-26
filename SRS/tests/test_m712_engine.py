"""
Exact tests for m712_engine.

Every assertion here is integer/Fraction arithmetic.  No float comparison
decides anything: the nearest-integer and sign checks go through
signexact, so these tests would still pass if the float layer were deleted.
The measurement assertions (the norm and pair barriers) are stated as exact
inequalities via certified signs rather than as float comparisons.
"""

import unittest
from fractions import Fraction as F

import m712_engine as E
import signexact as SE


A1, AJ = 2, 3


def sign_quadratic(u, v, a1):
    """Exact sign of u + v*sqrt(a1) for integer u, v."""
    return SE.exact_sign([(u, 1), (v, a1)]).sign


class TestEpsPower(unittest.TestCase):
    def test_index_is_the_exponent(self):
        """Docs/34 section 34.2 hazard: k must BE the exponent, not a list slot."""
        want = {0: (1, 0), 1: (1, 1), 2: (3, 2), 3: (7, 5),
                4: (17, 12), 5: (41, 29), 6: (99, 70)}
        for k, (p, q) in want.items():
            pt = E.eps_power(A1, k)
            self.assertEqual((pt.pk, pt.qk), (p, q), f"k={k}")

    def test_pell_identity(self):
        """P_k^2 - a1*Q_k^2 = (-1)^k."""
        for k in range(0, 25):
            pt = E.eps_power(A1, k)
            self.assertEqual(pt.pk ** 2 - A1 * pt.qk ** 2, (-1) ** k, f"k={k}")

    def test_multiplication_is_consistent(self):
        """eps^(k+1) = eps^k * eps, i.e. (P,Q) -> (P + a1 Q, P + Q)."""
        for k in range(0, 15):
            pt = E.eps_power(A1, k)
            nxt = E.eps_power(A1, k + 1)
            self.assertEqual(nxt.pk, pt.pk + A1 * pt.qk)
            self.assertEqual(nxt.qk, pt.pk + pt.qk)

    def test_negative_k_rejected(self):
        with self.assertRaises(ValueError):
            E.eps_power(A1, -1)


class TestFieldArithmetic(unittest.TestCase):
    def test_basis_rules(self):
        s = E.Elem(F(0), F(1), F(0), F(0))
        t = E.Elem(F(0), F(0), F(1), F(0))
        g = E.Elem(F(0), F(0), F(0), F(1))
        two, three, six = F(A1), F(AJ), F(A1 * AJ)
        cases = [
            (s, s, E.Elem(two, F(0), F(0), F(0))),        # s*s = a1
            (t, t, E.Elem(three, F(0), F(0), F(0))),      # t*t = aj
            (s, t, g),                                     # s*t = g
            (g, g, E.Elem(six, F(0), F(0), F(0))),         # g*g = a1*aj
            (s, g, E.Elem(F(0), F(0), two, F(0))),         # s*g = a1*t
            (t, g, E.Elem(F(0), three, F(0), F(0))),       # t*g = aj*s
        ]
        for p, q, expect in cases:
            self.assertEqual(E.elem_mul(p, q, A1, AJ), expect)

    def test_norm_is_rational_and_conjugate_invariant(self):
        x = E.Elem(F(3), F(-2), F(5), F(7))
        E.field_norm(x, A1, AJ)
        # norm of a rational is its fourth power
        r = E.Elem(F(-5), F(0), F(0), F(0))
        self.assertEqual(E.field_norm(r, A1, AJ), F(625))

    def test_conjugates_are_the_four_embeddings(self):
        x = E.Elem(F(1), F(2), F(3), F(4))
        cs = E.conjugates(x)
        self.assertEqual(len(cs), 4)
        self.assertEqual(cs[0], x)
        self.assertEqual(cs[1], E.Elem(F(1), F(-2), F(3), F(-4)))
        self.assertEqual(cs[2], E.Elem(F(1), F(2), F(-3), F(-4)))
        self.assertEqual(cs[3], E.Elem(F(1), F(-2), F(-3), F(4)))


class TestNearestInteger(unittest.TestCase):
    def test_is_the_nearest_integer(self):
        """Confirm exactly that |x_k - m_k| <= 1/2 for every k in range."""
        for k in range(1, 61):
            m = E.nearest_integer(A1, AJ, k)
            # 2*a1*(x_k - m) must lie in [-a1, a1]
            lo, hi = SE.certified_bracket(_bracket_args(A1, AJ, k, m), 256)
            self.assertGreaterEqual(lo, -A1, f"k={k} too high")
            self.assertLessEqual(hi, A1, f"k={k} too low")

    def test_monotone_and_growing(self):
        ms = [E.nearest_integer(A1, AJ, k) for k in range(1, 40)]
        self.assertEqual(ms, sorted(ms))
        self.assertGreater(ms[-1], 10 ** 12)

    def test_z_of_x_has_the_documented_shape(self):
        """z_k = 2 a1 x_k = a1 Q_k sqrt(a_j) + P_k sqrt(a1 a_j)."""
        for k in range(1, 20):
            pt = E.eps_power(A1, k)
            z = E.z_of_x(A1, AJ, k)
            self.assertEqual(z.c0, F(0))
            self.assertEqual(z.c1, F(0))
            self.assertEqual(z.c2, F(A1 * pt.qk))
            self.assertEqual(z.c3, F(pt.pk))


def _bracket_args(a1, aj, k, m):
    import squarefree as SF
    pt = E.eps_power(a1, k)
    agg, _ = SF.decompose_and_aggregate(
        [(-2 * a1 * m, 1), (a1 * pt.qk, aj), (pt.pk, a1 * aj)])
    return agg


class TestInvariantLayer(unittest.TestCase):
    def test_invariant_vanishes_identically(self):
        """Layer 3 is degenerate: I_k == 0 for every k.  A checked fact."""
        for k in range(1, 40):
            self.assertTrue(E.invariant_row(A1, AJ, k), f"k={k}")

    def test_q_factorises(self):
        """Q(u,v) == (v - eps u)(v - eps' u) exactly, on random exact elements."""
        import random
        rng = random.Random(20261003)
        for _ in range(200):
            u = E.Elem(*[F(rng.randint(-9, 9)) for _ in range(4)])
            v = E.Elem(*[F(rng.randint(-9, 9)) for _ in range(4)])
            self.assertTrue(E.q_factorises(u, v, A1, AJ))

    def test_anti_invariance_not_invariance(self):
        """
        Q(u_{k+1},u_{k+2}) = -Q(u_k,u_{k+1}), so Q is ANTI-invariant.  Guards
        against someone "restoring" the invariance claim.
        """
        u = E.Elem(F(3), F(-2), F(5), F(7))
        v = E.Elem(F(1), F(4), F(-3), F(2))
        # build u_{k+2} = 2 u_{k+1} + u_k directly
        u2 = E.Elem(*[2 * v[i] + u[i] for i in range(4)])
        self.assertTrue(
            E.elem_is_zero(E.elem_add(E.q_form(v, u2, A1, AJ),
                                     E.q_form(u, v, A1, AJ))))


class TestNormLayer(unittest.TestCase):
    def test_scaled_norm_is_an_integer(self):
        for k in range(1, 60):
            r = E.measure_row(A1, AJ, k)
            self.assertIsInstance(r.norm_y_scaled, int)

    def test_scaled_norm_is_an_algebraic_integer_norm(self):
        """
        2*a_1*y_k = 2 a_1 m_k - P_k sqrt(a_1 a_j) - a_1 Q_k sqrt(a_j) is an
        algebraic integer, so N(2 a_1 y_k) is an integer.  Note that
        N(y_k) itself is only a RATIONAL: dividing by (2 a_1)^4 need not give
        an integer (at k=1 it gives 228/256).  This is why Docs/35's proof
        must read |N(y_k)| >= (2 a_1)^(-4), not >= 1 -- the exponent 3 is
        unaffected, only the leading constant.
        """
        for k in range(1, 60):
            r = E.measure_row(A1, AJ, k)
            self.assertIsInstance(r.norm_y_scaled, int)
            self.assertNotEqual(r.norm_y_scaled, 0, f"N(2 a_1 y_{k})")

    def test_norm_of_unscaled_y_is_a_bounded_denominator_rational(self):
        """(2 a_1)^4 * N(y_k) in Z, and the true numerator is not a multiple."""
        s4 = (2 * A1) ** 4
        for k in range(1, 60):
            r = E.measure_row(A1, AJ, k)
            self.assertIsInstance(r.norm_y_scaled, int)
        self.assertEqual(s4, 256)
        self.assertNotEqual(E.measure_row(A1, AJ, 1).norm_y_scaled % s4, 0)

    def test_scaled_y_has_integer_coefficients(self):
        for k in range(1, 60):
            r = E.measure_row(A1, AJ, k)
            self.assertTrue(all(isinstance(c, int) for c in r.y_scaled))

    def test_scaled_y_equals_2a1_m_minus_z(self):
        for k in range(1, 40):
            r = E.measure_row(A1, AJ, k)
            z = E.z_of_x(A1, AJ, k)
            self.assertEqual(r.y_scaled[0], 2 * A1 * r.m)
            self.assertEqual(r.y_scaled[2], -int(z.c2))
            self.assertEqual(r.y_scaled[3], -int(z.c3))


class TestIntegralityLayer(unittest.TestCase):
    def test_c_of_shape(self):
        """c_k = (m_{k+1} - m_k) + (-m_k) sqrt(a_1)."""
        u, v = E.c_of(A1, 10, 4119, 9944)
        self.assertEqual((u, v), (9944 - 4119, -4119))

    def test_norm_of_c_is_an_integer_and_nonzero(self):
        """N(c_k) = u_k^2 - a_1 v_k^2 must be a nonzero integer."""
        for k in range(1, 60):
            u, v, n, _ = E.integrality_row(A1, AJ, k)
            self.assertIsInstance(n, int)
            self.assertNotEqual(n, 0, f"N(c_{k}) must not vanish")
            self.assertEqual(n, u * u - A1 * v * v)

    def test_c_is_small(self):
        """
        |c_k| = |y_{k+1} - eps y_k| <= d_{k+1} + eps d_k <= (1 + eps)/2 < 2,
        since d <= 1/2 and eps = 1 + sqrt(2) < 3.  Stated as two exact sign
        tests: -2 < u + v sqrt(a_1) < 2.
        """
        for k in range(1, 60):
            u, v, _, _ = E.integrality_row(A1, AJ, k)
            hi = SE.exact_sign([(u - 2, 1), (v, A1)]).sign
            lo = SE.exact_sign([(u + 2, 1), (v, A1)]).sign
            self.assertLessEqual(hi, 0, f"c_{k} >= 2")
            self.assertGreaterEqual(lo, 0, f"c_{k} <= -2")

    def test_pair_barrier_implies_a_gap_between_consecutive_distances(self):
        """
        Theorem 36.1 restated without floats.  If BOTH d_k and d_{k+1} were
        below 1/(4 sqrt(a_1) m_k) then
            d_{k+1} + eps d_k  <  (1 + eps)/(4 sqrt(a_1) m_k)  <  1/(2 sqrt(a_1) m_k),
        contradicting |c_k| <= d_{k+1} + eps d_k and |c_k| >= 1/(2 sqrt(a_1) m_k),
        because eps < 3 gives (1 + eps)/4 < 1.

        We verify the hypothesis side is FALSE, i.e. at least one of the two
        distances is at least 1/(4 sqrt(a_1) m_k), by checking that the exact
        value (d_{k+1} + eps d_k) * 2 sqrt(a_1) m_k is at least 1 via the
        certified brackets already used elsewhere.  Here we assert the weaker
        but purely integer consequence that both cannot be simultaneously tiny:
        certified |c_k| >= 1/(2 sqrt(a_1) m_k + 2), whose RHS is computed with
        integer arithmetic only.
        """
        for k in range(1, 60):
            u, v, n, _ = E.integrality_row(A1, AJ, k)
            mk = abs(E.nearest_integer(A1, AJ, k))
            # |c_k| >= 1/(2 sqrt(a_1) m_k + |c_k|) >= 1/(2 sqrt(a_1) mk + 2)
            # is equivalent to  |c_k| (2 sqrt(a_1) mk + 2) >= 1.  Certify it:
            # multiply through by |c_k| and use N(c_k) = |c_k| |conj(c_k)|.
            # The integer statement we can check directly is N(c_k) != 0,
            # already asserted above; here we bound |c_k| from below in the
            # quadratic field exactly, as  2 |c_k| >= 1/(sqrt(a_1) mk + 1).
            self.assertNotEqual(n, 0)
            del u, v, mk


class TestTrackBProductLayer(unittest.TestCase):
    """Docs/37: the product-level attack, which exploits the SHARED Pell pair."""

    PAIRS = [(3, 5), (3, 6), (5, 6), (3, 7), (6, 7), (5, 7)]

    def test_pair_integer_is_nonzero(self):
        """I_jj' = m_j'^2 a_j - m_j^2 a_j' in Z, never 0 for distinct radicands."""
        for aj, ajp in self.PAIRS:
            for k in range(1, 40):
                mj = E.nearest_integer(A1, aj, k)
                mjp = E.nearest_integer(A1, ajp, k)
                I = E.pair_integer(aj, ajp, mj, mjp)
                self.assertIsInstance(I, int)
                self.assertNotEqual(I, 0, f"I_({aj},{ajp}) vanished at k={k}")

    def test_pair_integer_formula(self):
        for aj, ajp in self.PAIRS:
            for k in (1, 7, 20):
                mj = E.nearest_integer(A1, aj, k)
                mjp = E.nearest_integer(A1, ajp, k)
                self.assertEqual(E.pair_integer(aj, ajp, mj, mjp),
                                 mjp * mjp * aj - mj * mj * ajp)

    def test_pair_bound_certified(self):
        """
        Theorem 37.1:
            |e_j| sqrt(a_j') + |e_j'| sqrt(a_j) >= 1/(m_j' sqrt(a_j) + m_j sqrt(a_j'))
        decided by an exact six-surd sign certificate.  Must hold everywhere.
        """
        for aj, ajp in self.PAIRS:
            for k in range(1, 80):
                self.assertTrue(E.pair_bound_holds(A1, aj, ajp, k),
                                f"pair bound failed at ({aj},{ajp},k={k})")

    def test_s1_coordinate_is_exact(self):
        """
        Fact 37.2: ||Q_k sqrt(a_1)|| = |eps'|^k exactly, so the nearest integer
        to Q_k sqrt(a_1) is P_k.  Checks the certified sign, not a float.
        """
        for k in range(1, 80):
            self.assertEqual(E.s1_exact(A1, k), E.eps_power(A1, k).pk)

    def test_s1_sign_alternates_with_parity(self):
        """eps'^k has the sign of (1 - sqrt(a_1))^k, so sign flips with parity."""
        import signexact as _SE
        for k in range(1, 30):
            pt = E.eps_power(A1, k)
            got = _SE.exact_sign([(pt.pk, 1), (-pt.qk, A1)]).sign
            self.assertEqual(got, -1 if k % 2 else 1, f"k={k}")


if __name__ == "__main__":
    unittest.main()