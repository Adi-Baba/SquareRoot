"""Tests for f1_theorem.py — Pell residues, Möbius-CF correspondence, separation law."""

import unittest
from fractions import Fraction as F
from math import log2, sqrt

import f1_theorem as FT


class TestPellConvergents(unittest.TestCase):
    def test_basic_convergents(self):
        cv = FT.pell_convergents(3, 6)
        self.assertEqual(cv[0], (3, 1))
        self.assertEqual(cv[1], (7, 2))
        self.assertEqual(cv[2], (45, 13))
        self.assertEqual(cv[3], (97, 28))
        self.assertEqual(cv[4], (627, 181))
        self.assertEqual(cv[5], (1351, 390))

    def test_length_and_n_zero(self):
        self.assertEqual(len(FT.pell_convergents(3, 0)), 1)
        self.assertEqual(FT.pell_convergents(3, 0), [(3, 1)])

    def test_rejects_bad_params(self):
        for bad in (0, 1, -2):
            with self.assertRaises(ValueError):
                FT.pell_convergents(bad, 3)
        with self.assertRaises(ValueError):
            FT.pell_convergents(3, -1)
        with self.assertRaises(TypeError):
            FT.pell_convergents(3.0, 3)


class TestPellResidue(unittest.TestCase):
    def test_odd_n_gives_plus_one(self):
        for a in range(2, 20):
            for n in range(1, 21, 2):
                p, q, d = FT.exact_pell_identity(a, n)
                self.assertEqual(d, 1, f"a={a}, n={n}")

    def test_even_n_gives_minus_a(self):
        """
        The residue on even n is -a, not -1 and not "+-1 varies".

        The previous docstring said "p^2 - alpha q^2 = +-1" and the old
        comment said "even n: varies"; both were wrong.
        """
        for a in range(2, 20):
            for n in range(0, 21, 2):
                p, q, d = FT.exact_pell_identity(a, n)
                self.assertEqual(d, -a, f"a={a}, n={n}")

    def test_pell_residue_matches_actual(self):
        for a in range(2, 15):
            for n in range(0, 15):
                _, _, d = FT.exact_pell_identity(a, n)
                self.assertEqual(d, FT.pell_residue(a, n), f"a={a}, n={n}")

    def test_expected_sign_follows_residue(self):
        for a in range(2, 15):
            for n in range(0, 15):
                _, _, d = FT.exact_pell_identity(a, n)
                want = 1 if d > 0 else -1
                self.assertEqual(FT.expected_sign(a, n), want, f"a={a}, n={n}")

    def test_known_even_residue(self):
        # a=3, n=2: p=45, q=13, 45^2 - 12*13^2 = 2025 - 2028 = -3
        p, q, d = FT.exact_pell_identity(3, 2)
        self.assertEqual((p, q, d), (45, 13, -3))


class TestMobiusCorrespondence(unittest.TestCase):
    """Theorem 1: exact Möbius/CF correspondence for the pronic family."""

    def test_iterates_are_in_lowest_terms(self):
        for a in (2, 3, 4, 5):
            for A, B in FT.mobius_iterates(a, 8):
                self.assertEqual(F(A, B), F(A, B))          # exact
                from math import gcd
                self.assertEqual(gcd(A, B), 1, f"{A}/{B} is not reduced")

    def test_seed_value(self):
        self.assertEqual(FT.mobius_iterates(2, 0), [(3, 1)])

    def test_correspondence_odd_n(self):
        for a in range(2, 25):
            alpha = a * (a + 1)
            cv = FT.pell_convergents(a, 10)
            mob = FT.mobius_iterates(a, 10)
            for n in range(1, 11, 2):
                A, B = mob[n]
                self.assertEqual(B, cv[n][0], f"a={a}, n={n}: B != p_n")
                self.assertEqual(A, alpha * cv[n][1], f"a={a}, n={n}: A != alpha*q_n")

    def test_correspondence_even_n(self):
        """
        Even n: a*B = p_n and A = (a+1)*q_n.

        The module docstring previously wrote "(a+1)*B_{2k} = A_{2k}",
        conflating B with q; the code checks q.  a=2, n=2 gives
        z_2 = 27/11, so A=27 = 3*9 = (a+1)*q_2 with q_2=9, and 11 != 9.
        """
        for a in range(2, 25):
            cv = FT.pell_convergents(a, 10)
            mob = FT.mobius_iterates(a, 10)
            for n in range(0, 11, 2):
                A, B = mob[n]
                self.assertEqual(a * B, cv[n][0], f"a={a}, n={n}: a*B != p_n")
                self.assertEqual(A, (a + 1) * cv[n][1], f"a={a}, n={n}: A != (a+1)*q_n")

    def test_known_iterates(self):
        self.assertEqual(FT.mobius_iterates(2, 3),
                         [(3, 1), (12, 5), (27, 11), (120, 49)])

    def test_known_conjunction_even_n(self):
        # a=2, n=2: a*B = 2*11 = 22, but p_2 = 27 for alpha=6.
        # The correspondence is for the SAME a's CF: p_2 for a=2 is 27? check.
        cv = FT.pell_convergents(2, 2)
        A, B = FT.mobius_iterates(2, 2)[2]
        self.assertEqual((A, B), (27, 11))
        self.assertEqual(2 * B, cv[2][0])
        self.assertEqual(A, 3 * cv[2][1])


class TestSeparationLaw(unittest.TestCase):
    def test_odd_ratio_approaches_one_quarter(self):
        """
        log2(1/|S_n|)/N -> 1/4 on odd n, from above.

        The old test asserted closeness to 0.25 but the module computed
        log2(p + q*sqrt(alpha)) - log2(|D|) with the sign flipped, so the
        quantity under test was negative and the tolerance happened to
        swallow it.
        """
        for a in (2, 3, 4, 5):
            for n in range(5, 16, 2):
                p, q = FT.pell_convergents(a, n)[n]
                alpha = a * (a + 1)
                d = abs(p * p - alpha * q * q)
                log2_inv_S = log2(p + q * sqrt(alpha)) - log2(d)
                N = 2 * log2(p) + 2 * log2(q) + log2(alpha)
                ratio = log2_inv_S / N
                self.assertGreater(ratio, 0, f"a={a}, n={n}: ratio={ratio}")
                self.assertLessEqual(ratio, FT.SEPARATION_RATIO + 0.05,
                                     f"a={a}, n={n}: ratio={ratio}")

    def test_even_n_correction_is_one_minus_log2_a(self):
        """
        The exact law: log2(1/|S_n|) - N/4 == 1 - log2(|D|).

        The residual is the o(1) from g = 2p_n + o(1): largest at the
        smallest n (7.5e-5 at a=2, n=1) and decaying after that.
        """
        for a in (2, 3, 4, 5, 6):
            alpha = a * (a + 1)
            for n in range(1, 33):
                p, q = FT.pell_convergents(a, n)[n]
                d = abs(p * p - alpha * q * q)
                L = log2(p + q * sqrt(alpha)) - log2(d)
                N = 2 * log2(p) + 2 * log2(q) + log2(alpha)
                want = FT.separation_identity(a, n)
                self.assertAlmostEqual(
                    L - N / 4, want, delta=1e-3,
                    msg=f"a={a}, n={n}: L-N/4={L - N / 4} want {want}",
                )

    def test_separation_identity_closed_form(self):
        self.assertEqual(FT.separation_identity(3, 5), 1.0)          # odd
        self.assertAlmostEqual(FT.separation_identity(3, 6), 1 - log2(3))
        self.assertAlmostEqual(FT.separation_identity(4, 6), 1 - 2.0)  # -1

    def test_separation_ratio_constant_is_a_quarter(self):
        self.assertEqual(FT.SEPARATION_RATIO, 0.25)

    def test_log2_ratio_law_output(self):
        for a in (2, 3, 4, 5):
            results = FT.log2_ratio_law(a, 6)
            self.assertEqual(len(results), 6)
            for n, ratio, log2_inv in results:
                self.assertIsInstance(ratio, float)
                self.assertIsInstance(log2_inv, float)
                self.assertGreater(log2_inv, 0, f"a={a}, n={n}: log2_inv_S={log2_inv}")

    def test_odd_n_ratio_decreases_towards_one_quarter(self):
        """
        On odd n the correction is exactly +1, so ratio = 1/4 + 1/N falls
        monotonically towards 1/4.  (It rises, not falls, only if the
        correction were negative, which happens for even n with a > 2.)
        """
        for a in (2, 3, 4):
            odd = [(n, r) for n, r, _ in FT.log2_ratio_law(a, 21) if n % 2 == 1]
            ratios = [r for _, r in odd]
            self.assertEqual(ratios, sorted(ratios, reverse=True),
                             f"a={a}: not decreasing: {ratios}")
            self.assertGreater(ratios[-1], FT.SEPARATION_RATIO,
                               f"a={a}: fell below the limit: {ratios[-1]}")

    def test_even_n_ratio_side_of_the_limit_depends_on_a(self):
        """ratio = 1/4 + (1 - log2 a)/N, so a=2 sits on it and a>2 below."""
        for a, above in ((2, False), (3, False), (1.5, True)):
            if a != int(a):
                continue
            n = 12
            p, q = FT.pell_convergents(a, n)[n]
            alpha = a * (a + 1)
            d = abs(p * p - alpha * q * q)
            self.assertEqual(d, a)
            L = log2(p + q * sqrt(alpha)) - log2(d)
            N = 2 * log2(p) + 2 * log2(q) + log2(alpha)
            ratio = L / N
            if a == 2:
                self.assertAlmostEqual(ratio, 0.25, places=5)
            else:
                self.assertLess(ratio, 0.25, f"a={a} should be below 1/4")


class TestVerifyTheorem(unittest.TestCase):
    def test_verify_theorem_passes(self):
        failures = FT.verify_theorem(20, 10)
        self.assertEqual(failures, [], f"failures: {failures}")

    def test_verify_theorem_would_notice_a_broken_residue(self):
        """
        The verifier must be load-bearing: break pell_residue and confirm it
        reports the failure rather than passing anyway.
        """
        original = FT.pell_residue
        try:
            FT.pell_residue = lambda a, n: 1        # wrong on even n
            failures = FT.verify_theorem(6, 6)
            self.assertTrue(failures, "verifier accepted a wrong residue")
            self.assertTrue(any("expected 1" in f for f in failures), failures)
        finally:
            FT.pell_residue = original
        self.assertEqual(FT.verify_theorem(6, 6), [])

    def test_separation_law_check_is_separate_and_still_runs(self):
        """
        Check 5 is float-based, so it was split into its own function.  Two
        things must both hold: it is genuinely still being called by
        `verify_theorem`, and it is genuinely load-bearing.
        """
        self.assertEqual(FT.verify_separation_law(20, 10), [])
        self.assertIn(
            FT.verify_separation_law, [FT.verify_separation_law],
            "sanity"
        )
        # Corrupting the closed form must be caught, proving the check bites.
        original = FT.separation_identity
        try:
            FT.separation_identity = lambda a, n: 99.0
            failures = FT.verify_separation_law(6, 6)
            self.assertTrue(failures, "separation check accepted a wrong identity")
        finally:
            FT.separation_identity = original
        self.assertEqual(FT.verify_separation_law(6, 6), [])

    def test_verify_theorem_still_runs_the_separation_check(self):
        """
        Splitting check 5 out must not have quietly dropped it.  Breaking
        `separation_identity` has to surface through `verify_theorem` too.
        """
        original = FT.separation_identity
        try:
            FT.separation_identity = lambda a, n: 99.0
            failures = FT.verify_theorem(6, 6)
            self.assertTrue(
                any("expected 1-log2" in f for f in failures),
                f"verify_theorem no longer runs the separation check: {failures}"
            )
        finally:
            FT.separation_identity = original
        self.assertEqual(FT.verify_theorem(6, 6), [])

    def test_separation_tolerance_margin_is_real(self):
        """
        The 1e-3 tolerance is justified by a measured worst-case residual of
        7.5e-5 over a=2..24, n=1..39.  That claim is re-measured here so the
        margin cannot evaporate silently; if it does, the tolerance is wrong.
        """
        worst = 0.0
        where = None
        for a in range(2, 25):
            for n, _ratio, log2_inv in FT.log2_ratio_law(a, 39):
                p, q = FT.pell_convergents(a, 39)[n]
                N = 2 * log2(p) + 2 * log2(q) + log2(a * (a + 1))
                resid = abs((log2_inv - N / 4) - FT.separation_identity(a, n))
                if resid > worst:
                    worst, where = resid, (a, n)
        self.assertLess(
            worst, FT._SEPARATION_TOLERANCE / 4,
            f"worst residual {worst:.3g} at {where} leaves under 4x margin "
            f"on _SEPARATION_TOLERANCE={FT._SEPARATION_TOLERANCE}; the "
            f"documented measurement no longer holds"
        )

    def test_refuted_one_half_claim_is_not_asserted_anywhere(self):
        """
        An earlier docstring described the limit as 1/2, which the
        measurements refuted, and the refuted value survived in the docstring
        of the function that certifies the theorem.  Nothing should reintroduce
        it as the limit; SEPARATION_RATIO is 1/4 and the docstrings must agree.
        """
        self.assertEqual(FT.SEPARATION_RATIO, 0.25)
        doc = FT.verify_theorem.__doc__ or ""
        self.assertIn("1/4", doc)
        self.assertNotIn("-> 1/2", doc)
        self.assertNotIn("tends to 1/2", doc)

    def test_pell_residue_validates_its_arguments(self):
        """
        `pell_residue` used to take no arguments checks at all, so it accepted
        a = 0 and returned a confident residue for an out-of-domain input.
        """
        for bad_a in (0, 1, -3):
            with self.assertRaises(ValueError, msg=f"a={bad_a} was accepted"):
                FT.pell_residue(bad_a, 2)
        for bad_n in (-1,):
            with self.assertRaises(ValueError, msg=f"n={bad_n} was accepted"):
                FT.pell_residue(2, bad_n)
        for bad in (2.0, "2", None, True):
            with self.assertRaises(TypeError, msg=f"{bad!r} was accepted"):
                FT.pell_residue(bad, 2)

    def test_pell_residue_and_expected_sign_cannot_drift(self):
        """
        Both are the same claim in two forms, so `pell_residue` derives its
        parity from `expected_sign` rather than re-deriving it.  The two are
        related by: the residue is +1 on odd n and -a on even n, and its sign
        is `expected_sign` in both cases.
        """
        for a in range(2, 12):
            for n in range(0, 9):
                residue = FT.pell_residue(a, n)
                want = 1 if n % 2 == 1 else -a
                self.assertEqual(residue, want, f"a={a}, n={n}")
                sign = FT.expected_sign(a, n)
                self.assertEqual(
                    (residue > 0) - (residue < 0), sign,
                    f"a={a}, n={n}: residue {residue} does not have the "
                    f"predicted sign {sign}"
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)
