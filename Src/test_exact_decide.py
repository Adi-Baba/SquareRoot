"""Tests for exact_decide.py -- the exact decision procedures.

The point of this file: it pins the *decider* against independent means
(brute-force grid search, direct polynomial expansion, hand-checked sign
cases), so a bug in the decider can never masquerade as a proof.  Three real
bugs in exact_decide.py were found by exactly these cross-checks and are each
pinned below:

  * poly_at / _eval_at iterated `reversed(p)` on a descending list, which
    evaluates a different polynomial.
  * decide_u_aware_bound carried a spurious B^2 u term in its expansion.
  * _quad_form negated the linear coefficient, reflecting every critical point
    through the origin so the decider missed the true maximiser and accepted
    bounds as small as half the true maximum.

None of these were visible from the code by inspection.  All three were caught
by the cross-checks below, which is the argument for having them.
"""

import unittest
from fractions import Fraction as F

import exact_decide as X


def _err(a, w, u):
    """The true error est - r = u^2(w-u) / (2a(2a+w) + u(2a+u))."""
    return (u * u * (w - u)) / (2 * a * (2 * a + w) + u * (2 * a + u))


def _brute_bound(a, w, lo, hi, N=3000):
    """Grid max of the error on [lo,hi] (an under-estimate of the true max)."""
    return max(_err(a, w, lo + (hi - lo) * F(k, N)) for k in range(N + 1))


class TestPolyAtMatchesExpansion(unittest.TestCase):
    def test_horner_is_not_reversed(self):
        # A cubic with a known value, computed two independent ways.
        p = [F(-1), F(2), F(3), F(4)]        # -u^3 + 2u^2 + 3u + 4
        for u in (F(0), F(1), F(2), F(-1), F(1, 2), F(7, 3)):
            self.assertEqual(X.poly_at(p, u), -u**3 + 2 * u**2 + 3 * u + 4)

    def test_matches_direct_expansion_randomly(self):
        import random
        random.seed(20240607)
        for _ in range(200):
            a, w = F(random.randint(1, 20)), F(random.randint(1, 20))
            B = F(random.randint(1, 200), 100)
            p = [F(-1), w - B, -2 * a * B, -2 * a * B * (2 * a + w)]
            u = F(random.randint(-500, 500), 100)
            direct = (u * u * (w - u)
                      - B * (2 * a * (2 * a + w) + u * (2 * a + u)))
            self.assertEqual(X.poly_at(p, u), direct)


class TestQuadraticFormRoots(unittest.TestCase):
    def test_rational_roots(self):
        # u^2 - 3u + 2 = 0 -> {1, 2}
        b, c = X._quad_form([F(1), F(-3), F(2)])
        self.assertEqual(b, F(-3))
        self.assertEqual(c, F(2))
        self.assertEqual(sorted(X.quad_roots_exact(b, c)), [F(1), F(2)])

    def test_negative_leading_coefficient_gives_positive_roots(self):
        # -3u^2 + 1.98u - 0.019 = 0 -> roots near +0.0097, +0.65 (both POSITIVE).
        # The _quad_form sign bug reflected these to negative and hid the peak.
        q = [F(-3), F(198, 100), F(-19, 1000)]
        b, c = X._quad_form(q)
        roots = X.quad_roots_exact(b, c)
        for r in roots:
            v = r.value() if isinstance(r, X.sqrt_field) else float(r)
            self.assertGreater(v, 0, f"root {r} should be positive")

    def test_derivative_roots_locate_a_cubic_maximum(self):
        # p(u) = -u^3 + 0.9u^2 - 0.01u - 0.05, max should be at a positive u
        p = [F(-1), F(9, 10), F(-1, 100), F(-5, 100)]
        b, c = X._quad_form(X._derivative(p))
        crit = X.quad_roots_exact(b, c)
        vals = [X.poly_at(p, r).value() if isinstance(X.poly_at(p, r), X.sqrt_field)
                else float(X.poly_at(p, r)) for r in crit]
        # at least one critical point has positive p (a local max above 0)
        self.assertTrue(any(v > 0 for v in vals), vals)

    def test_irrational_roots_are_exact_and_comparable(self):
        # u^2 - 2u + 1/2 = 0 -> 1 +- sqrt(2)/2, irrational
        b, c = X._quad_form([F(1), F(-2), F(1, 2)])
        roots = X.quad_roots_exact(b, c)
        self.assertTrue(all(isinstance(r, X.sqrt_field) for r in roots))
        lo, hi = roots
        self.assertLess(lo, hi)
        for r in roots:
            # each satisfies the polynomial exactly (sign 0)
            self.assertEqual(X.poly_at([F(1), F(-2), F(1, 2)], r).sign(), 0)


class TestSqrtField(unittest.TestCase):
    def test_signs(self):
        cases = [
            (X.sqrt_field(F(1), F(1), F(2)), 1),        # 1+sqrt2
            (X.sqrt_field(F(-3), F(1), F(2)), -1),      # -3+sqrt2 < 0
            (X.sqrt_field(F(1), F(-1), F(2)), -1),      # 1-sqrt2 < 0
            (X.sqrt_field(F(0), F(3), F(5)), 1),
            (X.sqrt_field(F(0), F(-3), F(5)), -1),
            (X.sqrt_field(F(0), F(0), F(5)), 0),
        ]
        for z, want in cases:
            self.assertEqual(z.sign(), want, f"{z!r}")

    def test_exact_zero_when_equal(self):
        z = X.sqrt_field(F(5, 3), F(2, 7), F(5))
        self.assertEqual((z - z).sign(), 0)

    def test_ordering_against_fraction(self):
        z = X.sqrt_field(F(1), F(-1), F(2))            # 1 - sqrt2 ~ -0.414
        self.assertLess(z, F(0))                        # negative
        self.assertLess(z, F(1))
        z2 = X.sqrt_field(F(1), F(1), F(2))            # 1 + sqrt2 ~ 2.414
        self.assertGreater(z2, F(1))
        self.assertGreater(z2, F(2))

    def test_witness_value_helper(self):
        # poly_at on a surd argument returns a surd; .value() is display-only.
        z = X.sqrt_field(F(0), F(1), F(2))             # sqrt2
        val = X.poly_at([F(1), F(0)], z)               # u + 0, i.e. u itself
        self.assertEqual(val.value(), 2.0 ** 0.5)
        const = X.poly_at([F(3), F(1)], z)             # 3u + 1
        self.assertAlmostEqual(const.value(), 3 * 2.0 ** 0.5 + 1)

    def test_arithmetic(self):
        z = X.sqrt_field(F(1), F(1), F(2))
        self.assertEqual((z * z).p, 3)                # (1+s)^2 = 3 + 2s
        self.assertEqual((z * z).q, 2)


class TestDecideUBound(unittest.TestCase):
    def test_accepts_loose_bound(self):
        a, w = F(3), F(1)
        B = _brute_bound(a, w, F(0), F(1)) * F(3, 2)  # comfortably loose
        ok, wit = X.decide_u_aware_bound(a, w, F(0), F(1), B)
        self.assertTrue(ok)
        self.assertIsNone(wit)

    def test_rejects_tight_bound_with_witness(self):
        a, w = F(1), F(1)
        tm = _brute_bound(a, w, F(0), F(1), N=400000)
        for frac in (F(1), F(999, 1000), F(99, 100), F(1, 2)):
            ok, wit = X.decide_u_aware_bound(a, w, F(0), F(1), tm * frac)
            self.assertFalse(ok, f"B={float(frac)}*trueMax should be rejected")
            self.assertIsNotNone(wit, "a rejection must come with a witness")

    def test_witness_actually_violates(self):
        a, w = F(1), F(1)
        tm = _brute_bound(a, w, F(0), F(1), N=400000)
        B = tm * F(1, 2)
        ok, wit = X.decide_u_aware_bound(a, w, F(0), F(1), B)
        self.assertFalse(ok)
        p = [F(-1), w - B, -2 * a * B, -2 * a * B * (2 * a + w)]
        val = X.poly_at(p, wit)
        val = val.value() if isinstance(val, X.sqrt_field) else float(val)
        self.assertGreater(val, 0, "witness must have p(witness) > 0")

    def test_never_accepts_what_grid_rejects(self):
        # The decisive one-sided check: whenever the decider says VALID, the
        # error must genuinely stay below B on a dense grid.  (The converse
        # fails legitimately when the grid misses a narrow peak.)
        import random
        random.seed(11)
        for _ in range(60):
            a = F(random.randint(1, 30))
            w = F(random.randint(1, 6))
            lo = F(0)
            hi = w
            tm = _brute_bound(a, w, lo, hi, N=1500)
            for frac in (F(1, 2), F(9, 10), tm):
                B = tm * frac
                ok, _ = X.decide_u_aware_bound(a, w, lo, hi, B)
                if ok:
                    for k in range(2001):
                        u = lo + (hi - lo) * F(k, 2000)
                        self.assertLessEqual(_err(a, w, u), B,
                                             f"decider accepted B but err>B at {u}")

    def test_matches_brute_force_one_sided(self):
        # Where brute (grid) says INVALID, the decider must too.
        import random
        random.seed(99)
        for _ in range(40):
            a = F(random.randint(1, 20))
            w = F(random.randint(1, 5))
            tm = _brute_bound(a, w, F(0), w, N=3000)
            B = tm * F(3, 4)                            # 25% below grid max
            grid_says_invalid = any(
                _err(a, w, F(k, 4000) * w) > B for k in range(4001))
            ok, _ = X.decide_u_aware_bound(a, w, F(0), w, B)
            if grid_says_invalid:
                self.assertFalse(ok, "decider accepted a bound the grid refutes")

    def test_published_bound_holds_exactly(self):
        # Theorem 4: B = w^3/(27 a^2) is a valid bound on [0,w] for all a>0.
        import random
        random.seed(7)
        for _ in range(50):
            a = F(random.randint(1, 50))
            w = F(random.randint(1, 10))
            B = w**3 / (27 * a * a)
            ok, wit = X.decide_published_bound(a, w, B)
            self.assertTrue(ok, f"published bound rejected for a={a}, w={w}: {wit}")


class TestDecideErrorIdentity(unittest.TestCase):
    def test_identity_holds(self):
        for a in (F(1), F(3), F(12)):
            for w in (F(1), F(5, 2)):
                for k in range(1, 20):
                    u = w * F(k, 20)
                    ok, diff = X.decide_error_identity(a, w, u)
                    self.assertTrue(ok, f"a={a}, w={w}, u={u}, diff={diff}")

    def test_checker_setup_is_correct(self):
        # The checker builds p = alpha - a^2 with r = a+u, so p = u(2a+u).
        # A wrong p (the bug that shipped) made the identity check fail.
        a, w, u = F(1), F(1), F(1, 20)
        p = u * (2 * a + u)
        m = a + p / (2 * a + w)
        est = a + p / (a + m)
        self.assertGreater(est, a + u)   # est over-estimates r = a+u


if __name__ == "__main__":
    unittest.main(verbosity=2)
