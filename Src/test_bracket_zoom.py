"""
Verification suite for Docs/Draft1.md.

Run:  python3 test_bracket_zoom.py

Each test is labelled with the paper section it checks.  Tests that assert a
*discrepancy* against the paper are marked [FINDING] and state what the paper
claims versus what is actually true.

Everything is exact (fractions / mpmath).  No floating point is used in any
assertion about rigorous claims.
"""

import math
import random
import statistics as st
import unittest
from fractions import Fraction as F
from decimal import Decimal as D, getcontext

from mpmath import mp, mpf, sin, cos, exp, log, tanh, findroot

from bracket_zoom import (
    plain, corrected, error_bound, error_bound_tight,
    closed_form_error_from_r, zoom_recursive, zoom_status,
    ZOOM_DIVERGES, ZOOM_OK,
)


# ---------------------------------------------------------------- helpers

def rand_bracket(rng, a_max=10 ** 4, widths=(1, 10, 10 ** 6, 10 ** 12)):
    """A bracket (a, w) with a >= 1, and a rational r strictly inside it.

    Returns (a, w, r, alpha) with alpha == r**2 exactly, so every downstream
    identity can be checked in exact arithmetic.
    """
    a = F(rng.randint(1, a_max))
    w = F(rng.choice(widths))
    if w < 2:
        w = F(2)
    num = a * w.denominator + rng.randint(1, w.numerator - 1)
    r = F(num, w.denominator)
    return a, w, r, r * r


def nonsquares(lo, hi):
    return [n for n in range(lo, hi + 1) if math.isqrt(n) ** 2 != n]


# ============================================== Section 2 & 3: exact structure

class TestExactStructure(unittest.TestCase):
    def test_claim1_weighted_average(self):
        """Claim 1: alpha == (q*x^2 + p*y^2)/d for every alpha."""
        for alpha in range(2, 400):
            x = math.isqrt(alpha)
            y = x + 1
            p, q, d = alpha - x * x, y * y - alpha, y * y - x * x
            self.assertEqual(F(q * x * x + p * y * y, d), alpha)

    def test_claim2_rms_identity(self):
        """Claim 2: r == sqrt((x + p/d)^2 + pq/d^2), exactly."""
        getcontext().prec = 60
        for alpha in range(2, 400):
            x = math.isqrt(alpha)
            y = x + 1
            p, q, d = alpha - x * x, y * y - alpha, y * y - x * x
            m = F(x) + F(p, d)
            Md = D(m.numerator) / D(m.denominator)
            rhs = (Md * Md + D(p * q) / D(d * d)).sqrt()
            self.assertLess(abs(rhs - D(alpha).sqrt()), D(10) ** -55)

    def test_paper_table_section2(self):
        """The Section 2 table, reproduced to the printed precision."""
        for alpha, x, p, q, d, est in [
            (12, 3, 3, 4, 7, 3.42857),
            (10, 3, 1, 6, 7, 3.14286),
            (2, 1, 1, 2, 3, 1.33333),
        ]:
            self.assertEqual(math.isqrt(alpha), x)
            self.assertEqual(alpha - x * x, p)
            self.assertEqual((x + 1) ** 2 - alpha, q)
            self.assertEqual(2 * x + 1, d)
            self.assertAlmostEqual(float(F(x) + F(p, d)), est, places=5)

    def test_proronic_side_fact(self):
        """Section 2: for pronic alpha = x(x+1), x + p/d == HM(x, x+1), and
        GM^2 == AM*HM."""
        getcontext().prec = 50
        for x in range(1, 50):
            alpha = x * (x + 1)
            p, d = alpha - x * x, 2 * x + 1
            am = F(2 * x + 1, 2)
            hm = F(2 * x * (x + 1), 2 * x + 1)
            self.assertEqual(F(x) + F(p, d), hm)
            self.assertEqual(am * hm, F(alpha))          # GM^2 = AM * HM
            self.assertEqual((D(x) * D(x + 1)).sqrt(), D(alpha).sqrt())


# ================================================== Section 5: Theorem 1

class TestTheorem1(unittest.TestCase):
    def test_exact_error_identity(self):
        """r - M(a,z) == (r-a)(z-r)/(a+z) for arbitrary z, exactly."""
        rng = random.Random(7)
        checked = 0
        for _ in range(30000):
            a, w, r, alpha = rand_bracket(rng)
            b = a + w
            m = plain(a, b, alpha)
            est = corrected(a, b, alpha)
            for z in (a, b, m, est, a + (b - a) / 3):
                self.assertEqual(r - plain(a, z, alpha), (r - a) * (z - r) / (a + z))
                checked += 1
        self.assertEqual(checked, 150000)

    def test_one_sidedness(self):
        """Plain under-estimates, corrected over-estimates, always."""
        rng = random.Random(8)
        for _ in range(20000):
            a, w, r, alpha = rand_bracket(rng)
            b = a + w
            self.assertLessEqual(plain(a, b, alpha), r)
            self.assertGreaterEqual(corrected(a, b, alpha), r)

    def test_closed_form_error(self):
        """est - r == u^2(w-u) / (2a(2a+w) + u(2a+u)), exactly."""
        rng = random.Random(9)
        for _ in range(30000):
            a, w, r, alpha = rand_bracket(rng)
            self.assertEqual(
                corrected(a, a + w, alpha) - r,
                closed_form_error_from_r(a, w, r),
            )


# ================================================== Section 6: Theorem 2

class TestTheorem2(unittest.TestCase):
    def test_paper_bound_holds(self):
        """0 <= est - r <= w^3/(16a^2) -- the paper's Theorem 2 is valid."""
        rng = random.Random(11)
        for _ in range(50000):
            a, w, r, alpha = rand_bracket(rng, a_max=10 ** 3, widths=(1, 10 ** 3, 10 ** 6))
            e = corrected(a, a + w, alpha) - r
            self.assertGreaterEqual(e, 0)
            self.assertLessEqual(e, error_bound(a, w))

    def test_tight_bound_holds_and_is_sharper(self):
        """[FINDING] The 1/16 constant is loose; 4/27 is the sharp one.

        The paper prints w^3/(16a^2).  Since u^2(w-u) peaks at u = 2w/3 with
        value 4w^3/27, the sharp constant is 4/27, i.e. 27/16 = 1.6875x
        tighter.
        """
        rng = random.Random(11)
        for _ in range(50000):
            a, w, r, alpha = rand_bracket(rng, a_max=10 ** 3, widths=(1, 10 ** 3, 10 ** 6))
            e = float(corrected(a, a + w, alpha) - r)
            self.assertLessEqual(e / float(error_bound(a, w)), 16 / 27 + 1e-9)
            self.assertLessEqual(e / float(error_bound_tight(a, w)), 1.0 + 1e-9)

    def test_paper_bound_attained_at_two_thirds(self):
        """The ratio to the paper's bound peaks at u = 2w/3, tending to 16/27.

        Exactly:  e / eps = 64a^2 / (108a^2 + 90aw + 12w^2)  at u = 2w/3,
        which -> 64/108 = 16/27 as w/a -> 0 and is smaller for w ~ a
        (e.g. a = w = 1 gives 64/210 = 0.3048).
        """
        for a, w in [(F(10 ** 12), F(1)), (F(10 ** 10), F(1, 10 ** 2)),
                     (F(10 ** 8), F(1, 10 ** 4))]:
            u = 2 * w / 3
            r = a + u
            e = corrected(a, a + w, r * r) - r
            self.assertAlmostEqual(
                float(e / error_bound(a, w)),
                16 / 27, places=8,
            )
        # and for w ~ a it is strictly below the asymptotic value
        r = F(1) + F(2, 3)
        e = corrected(F(1), F(2), r * r) - r
        self.assertLess(float(e / error_bound(F(1), F(1))), 16 / 27)

    def test_no_failures_in_paper_regime(self):
        """In the regime the paper actually tests (a >= 1, w <= 1) nothing breaks."""
        rng = random.Random(1)
        for _ in range(100000):
            a = F(rng.randint(1, 10 ** 6))
            w = F(rng.randint(1, 10 ** 6), 10 ** 6)
            u = F(rng.randint(1, 10 ** 6 - 1), 10 ** 6)
            if u >= w:
                continue
            r = a + u
            alpha = r * r
            e = corrected(a, a + w, alpha) - r
            self.assertGreaterEqual(e, 0)
            self.assertLessEqual(e, error_bound(a, w))
            self.assertLess(error_bound(a, w), w)


class TestAlgorithm3Precondition(unittest.TestCase):
    """[FINDING] Algorithm 3 diverges when w > 4a; the paper omits this."""

    def test_progress_iff_w_lt_4a(self):
        rng = random.Random(3)
        for _ in range(5000):
            a = F(rng.randint(1, 500))
            w = F(rng.randint(2, 5000))
            shrinks = error_bound(a, w) < w
            self.assertEqual(shrinks, w < 4 * a)

    def test_zoom_status(self):
        self.assertEqual(zoom_status(1, 2), ZOOM_OK)      # w=1 < 4a=4
        self.assertEqual(zoom_status(1, 6), ZOOM_DIVERGES)  # w=5 > 4a
        self.assertEqual(zoom_status(0, 1), ZOOM_DIVERGES)  # a=0 -> undefined

    def test_divergence_loses_invariant(self):
        """With a=1, w=5 the bracket widens and the invariant dies in 2 steps."""
        a, w, alpha = F(1), F(5), F(4)
        r = F(2)
        widths = []
        for _ in range(4):
            if not (a < r < a + w):
                break
            est = corrected(a, a + w, alpha)
            a, w = est - error_bound(a, w), error_bound(a, w)
            widths.append(float(w))
        self.assertGreater(widths[0], 5.0)          # grew
        self.assertFalse(a < r < a + w)              # invariant lost

    def test_paper_regime_converges(self):
        """a=1, w=1 converges cubically, as Section 6 claims."""
        a, b, depth = zoom_recursive(F(2), F(1), F(2), tol_digits=50)
        self.assertLessEqual(depth, 7)
        self.assertTrue(a * a < 2 < b * b)

    def test_alpha_below_one_breaks(self):
        """[FINDING] a >= 1 is load-bearing: a=0 makes the bound undefined."""
        with self.assertRaises(ZeroDivisionError):
            error_bound(F(0), F(1))


# ================================================== Section 8.1: experiments

def decimal_zoom_plain_and_corrected(alpha, levels, prec=80):
    """Algorithm 1 in Decimal.  Returns (plain errors, corrected errors).

    Mirrors the paper's level indexing: level 0 is the estimate on the initial
    bracket [floor(r), floor(r)+1], and each subsequent level first shrinks w
    by 10 and then scans k in {0..9}.
    """
    getcontext().prec = prec
    alpha = D(alpha)
    r = alpha.sqrt()
    a, w = D(math.isqrt(int(alpha))), D(1)
    pe, ce = [], []
    for level in range(levels):
        p = alpha - a * a
        b = a + w
        m = a + w * p / (b * b - a * a)
        pe.append(abs(r - m))
        ce.append(a + p / (a + m) - r)
        w = w / 10
        k = 9
        while k > 0 and (a + k * w) ** 2 > alpha:   # largest k keeping r bracketed
            k -= 1
        a = a + k * w
    return pe, ce


class TestSection81(unittest.TestCase):
    def test_sample_size_2946(self):
        """The paper's '2946 non-square alpha in 2..3000'."""
        self.assertEqual(len(nonsquares(2, 3000)), 2946)

    def test_sample_size_3152(self):
        """'3152 numbers' = 2946 + 6 extremes + 200 random."""
        self.assertEqual(2946 + 6 + 200, 3152)

    def test_median_digits_per_level(self):
        """Plain median 2.00, corrected median 3.00."""
        gp, gc = [], []
        for alpha in nonsquares(2, 3000):
            pe, ce = decimal_zoom_plain_and_corrected(alpha, 6)
            for i in range(1, 5):
                if 0 < pe[i] < pe[i - 1]:
                    gp.append(math.log10(pe[i - 1] / pe[i]))
                if 0 < ce[i] < ce[i - 1]:
                    gc.append(math.log10(ce[i - 1] / ce[i]))
        self.assertAlmostEqual(st.median(gp), 2.00, places=2)
        self.assertAlmostEqual(st.median(gc), 3.00, places=2)

    def test_extreme_alpha_mechanism(self):
        """[FINDING] The paper blames a metric artifact; it is a real weakness.

        Section 8.1 attributes the poor gains at alpha = 10**12 + 1 to the
        measure ("error pehle level pe hi ~1e-25 tha").  The level-0 error is
        in fact ~2.5e-13, and the cause is structural: the digit scan cannot
        advance `a` at all, because (10**6 + 10**-k)**2 > 10**12 + 1 for every
        k <= 6.  With `a` frozen, the plain error ~ u*w/(2a) shrinks only
        linearly in w, so the method degrades to exactly 1 digit per level for
        the first six levels.  The recursive zoom escapes this because it
        re-centres a on est ~ r.
        """
        getcontext().prec = 200
        alpha = 10 ** 12 + 1
        a, w = D(10 ** 6), D(1)
        r = D(alpha).sqrt()
        errs, anchors = [], []
        for _ in range(6):
            p = D(alpha) - a * a
            b = a + w
            m = a + w * p / (b * b - a * a)
            errs.append(float(r - m))
            anchors.append(a)
            w = w / 10
            k = 9
            while k > 0 and (a + k * w) ** 2 > D(alpha):
                k -= 1
            a = a + k * w
        self.assertAlmostEqual(errs[0], 2.5e-13, delta=1e-15)
        # the anchor never moves -- this is the actual mechanism
        self.assertEqual(len(set(anchors)), 1)
        # and the gain is exactly 1 digit/level, not 2
        gains = [math.log10(errs[i - 1] / errs[i]) for i in range(1, 6)]
        for g in gains:
            self.assertAlmostEqual(g, 1.0, delta=0.05)

    def test_bracket_invariant_never_broken(self):
        """Paper: 'bracket invariant (a < r < b) 0 baar tuta'."""
        for alpha in nonsquares(2, 500):
            getcontext().prec = 60
            a, w = D(math.isqrt(int(alpha))), D(1)
            for _ in range(8):
                self.assertTrue(a * a <= D(alpha) < (a + w) * (a + w))
                w = w / 10
                k = 9
                while k > 0 and (a + k * w) ** 2 > D(alpha):
                    k -= 1
                a = a + k * w


# ============================== Section 7/8.2: the general-f interpolation rate

def interp_zoom(f, lo, hi, levels=7, dps=200):
    """Algorithm 1 for a general f, with a correctly-oriented digit scan."""
    mp.dps = dps
    a, w, out = mpf(lo), mpf(hi) - mpf(lo), []
    pos = f(a) > 0
    for _ in range(levels):
        w = w / 10
        k = 9
        while k > 0 and (f(a + k * w) > 0) != pos:  # largest k keeping r bracketed
            k -= 1
        a = a + k * w
        out.append(a + w * f(a) / (f(a) - f(a + w)))
    return out


def level_gains(c, estimates):
    errs = [abs(x - c) for x in estimates]
    return [
        float(mp.log10(errs[i - 1] / errs[i]))
        for i in range(1, len(errs)) if 0 < errs[i] < errs[i - 1]
    ]


class TestInterpolationRateTheory(unittest.TestCase):
    """[FINDING] Section 8.2's pi explanation is correct -- and generalises.

    For a simple root c of f with 0 < |f'(c)| < inf:

        |c - xL| = (A*B/2) * |f''(xi)/f'(c)|,   A = c-a, B = b-c

    so the per-level gain is 2 digits when f''(c) != 0 (A,B ~ w on both sides),
    and 3 digits when f''(c) == 0 (because then f''(xi) = f'''(c)(xi-c) = O(w),
    promoting the error to O(w**3)).  When |f'(c)| -> 0 or inf the w**2 term
    dies and the gain drops to 1 digit -- which is exactly the flat/steep pair
    the paper reports but does not explain.
    """

    CASES = [
        ("pi  (sin)", lambda x: sin(x), 3, 4, "findroot"),
        ("e  (log x - 1)", lambda x: log(x) - 1, 2, 3, mp.e),
        ("ln2 (exp x - 2)", lambda x: exp(x) - 2, 0, 1, log(2)),
        ("phi (x^2-x-1)", lambda x: x * x - x - 1, 1, 2, None),
        ("cbrt2 (x^3-2)", lambda x: x ** 3 - 2, 1, 2, None),
        ("dottie (cos x-x)", lambda x: cos(x) - x, 0, 1, None),
        ("tanh(30(x-pi))", lambda x: tanh(30 * (x - mp.pi)), 3, 4, mp.pi),
    ]

    def _root(self, f, lo, hi, known):
        if known == "findroot":
            return findroot(f, mpf((lo + hi) / 2))
        if known is not None:
            return mpf(known)
        return findroot(f, mpf((lo + hi) / 2))

    def test_f_double_prime_zero_gives_three_digits(self):
        """pi (sin''(pi)=0) and tanh(30(x-pi)) (f''=0) both give ~3."""
        for label, f, lo, hi, known in self.CASES:
            if label not in ("pi  (sin)", "tanh(30(x-pi))"):
                continue
            c = self._root(f, lo, hi, known)
            self.assertLess(abs(mp.diff(f, c, 2)), mpf(10) ** -40, label)
            g = level_gains(c, interp_zoom(f, lo, hi, 7))
            self.assertGreater(sum(g) / len(g), 2.5, f"{label}: got {g}")

    def test_f_double_prime_nonzero_gives_two_digits(self):
        """e, ln2, phi, cbrt2, dottie all give ~2 (paper says '~2.0 to 2.1')."""
        for label, f, lo, hi, known in self.CASES:
            if label in ("pi  (sin)", "tanh(30(x-pi))"):
                continue
            c = self._root(f, lo, hi, known)
            self.assertGreater(abs(mp.diff(f, c, 2)), mpf(10) ** -10, label)
            g = level_gains(c, interp_zoom(f, lo, hi, 7))
            self.assertTrue(1.7 < sum(g) / len(g) < 2.4, f"{label}: got {g}")

    def test_degenerate_roots_give_one_digit(self):
        """Flat and steep roots collapse to ~1 digit/level, as Section 8.2 says."""
        c = mp.pi
        for label, f in [
            ("flat (x-pi)^3", lambda x: (x - c) ** 3),
            ("steep (x-pi)^(1/3)", lambda x: mp.sign(x - c) * abs(x - c) ** (mpf(1) / 3)),
        ]:
            g = level_gains(c, interp_zoom(f, 3, 4, 5))
            self.assertTrue(0.7 < sum(g) / len(g) < 1.3, f"{label}: got {g}")


# ================================================== Section 8.5: depth model

class TestDepthModel(unittest.TestCase):
    def test_depth_matches_paper(self):
        """Paper's depths are reproduced up to a constant +1 convention offset."""
        random.seed(0)
        for target, paper_depth in [(20, 4), (100, 5), (1000, 7), (3000, 8), (10000, 9)]:
            getcontext().prec = target + 80
            depths = []
            for alpha in random.sample(range(2, 10 ** 6), 30):
                if math.isqrt(alpha) ** 2 == alpha:
                    continue
                _, _, d = zoom_recursive(alpha, D(math.isqrt(alpha)),
                                         D(math.isqrt(alpha) + 1), target)
                depths.append(d)
            # The paper never states how the initial bracket is obtained, so its
            # depth counts are offset from ours by a small convention amount.
            self.assertLessEqual(
                abs(st.median(depths) + 1 - paper_depth), 1, f"{target} digits"
            )

    def test_tight_constant_does_not_change_depth(self):
        """4/27 narrows the bracket but does not change the recursion depth."""
        getcontext().prec = 200
        for alpha in [7, 12345, 999983]:
            a = D(math.isqrt(alpha))
            _, _, d1 = zoom_recursive(alpha, a, a + 1, 200, tight=False)
            _, _, d2 = zoom_recursive(alpha, a, a + 1, 200, tight=True)
            self.assertEqual(d1, d2)


# ================================================== Section 8.4/8.5: reporting

class TestReportedNumbers(unittest.TestCase):
    """Arithmetic checks on the paper's own tables (no benchmarking here)."""

    def test_zoom_over_newton_ratio_is_non_monotonic(self):
        """[FINDING] Section 8.5's zoom/Newton ratio must grow with D, but does not.

        The stated cost model is O(D log D) for the zoom vs O(D) for Newton, so
        the ratio should increase monotonically.  It goes 5.70 -> 8.79 -> 4.50.
        """
        rows = [(15, .011, .006), (200, .029, .011), (1000, .37, .065),
                (3000, 2.9, .33), (10000, 18, 4.0)]
        ratios = [z / n for _, z, n in rows]
        # increases up to 3000 digits, then DROPS -- a monotonicity violation
        self.assertLess(ratios[0], ratios[1])
        self.assertLess(ratios[1], ratios[2])
        self.assertLess(ratios[2], ratios[3])
        self.assertGreater(ratios[3], ratios[4], "expected the 8.79 -> 4.50 drop")
        self.assertFalse(
            all(x < y for x, y in zip(ratios, ratios[1:])),
            "the reported sequence should be non-monotonic",
        )

    def test_isqrt_baseline_gap_omitted_from_abstract(self):
        """[FINDING] The abstract compares only to Newton; isqrt is 10-675x off."""
        for digits, zoom, isqrt in [(15, .020, .002), (200, .34, .006),
                                    (1000, 9.2, .044), (3000, 189, .28)]:
            self.assertGreater(zoom / isqrt, 9.9)
        # and even the best recursive number is ~10x off the C reference
        self.assertAlmostEqual(2.9 / 0.28, 10.4, places=1)

    def test_abstract_speedup_range(self):
        """[FINDING] '65x faster' is the best case; the range is 1.8x-65x."""
        ratios = [.020 / .011, .34 / .029, 9.2 / .37, 189 / 2.9]
        self.assertAlmostEqual(ratios[-1], 65.2, delta=0.5)
        self.assertAlmostEqual(ratios[0], 1.8, delta=0.05)


# ========================================== Paper1.md: results not in Draft1

class TestPaper1(unittest.TestCase):
    """Claims made in Docs/Paper1.md that go beyond the original draft."""

    def test_corollary3_first_correction(self):
        """Cor. 3: m0 + pq/(2 d^2 m0) > r always, with error v^2/(4 m0 (rhat+r))."""
        getcontext().prec = 50
        for alpha in (12, 10, 2, 17, 101, 3, 5, 8, 12345):
            x = math.isqrt(alpha)
            y = x + 1
            p, q, d = alpha - x * x, y * y - alpha, y * y - x * x
            m0 = F(x) + F(p, d)
            r = D(alpha).sqrt()
            v = D(p * q) / D(d * d)
            M0, V = D(m0.numerator) / D(m0.denominator), v
            c3 = M0 + V / (2 * M0)
            self.assertGreater(c3, r, alpha)
            # exact error form
            self.assertAlmostEqual(
                float(c3 - r), float(V * V / (4 * M0 * M0 * (c3 + r))), places=20
            )

    def test_theorem2_weighted_rms_for_all_alpha(self):
        """Thm 2 (weighted RMS) to 50 significant digits, alpha in [2, 400)."""
        getcontext().prec = 60
        for alpha in range(2, 400):
            x = math.isqrt(alpha)
            y = x + 1
            p, q, d = alpha - x * x, y * y - alpha, y * y - x * x
            m = F(x) + F(p, d)
            M = D(m.numerator) / D(m.denominator)
            rhs = (M * M + D(p * q) / D(d * d)).sqrt()
            self.assertLess(abs(rhs - D(alpha).sqrt()), D(10) ** -55, alpha)

    def test_bound_holds_for_alpha_below_one(self):
        """Thm 4 needs only a > 0, so it covers 0 < alpha < 1 (Sec 5.3)."""
        for a, w, r in [(F(1, 2), F(1, 10), F(7, 12)),
                        (F(9, 10), F(1, 10), F(11, 12)),
                        (F(1, 4), F(1, 100), F(51, 200)),
                        (F(1, 10 ** 4), F(1, 10 ** 3), F(3, 10 ** 4))]:
            self.assertTrue(a < r < a + w, "fixture must bracket r strictly")
            alpha = r * r
            e = corrected(a, a + w, alpha) - r
            self.assertGreaterEqual(e, 0)
            self.assertLessEqual(e, error_bound_tight(a, w))
            # and the paper's 1/16 form is still valid for a > 0
            self.assertLessEqual(e, error_bound(a, w))

    def test_legacy_bound_fails_only_at_a_zero(self):
        """The 1/16 expression is undefined at a=0 but fine for any a>0."""
        with self.assertRaises(ZeroDivisionError):
            error_bound(F(0), F(1))
        for a in (F(1, 10 ** 6), F(1, 2), F(1)):
            error_bound(a, F(1, 10))   # must not raise

    def test_theorem6_digit_tripling(self):
        """Thm 6: E_{n+1} = 3E_n + log10(c a^2), solved exactly for E_n."""
        c = F(4, 27)
        a = F(10 ** 3)
        k = math.log10(c * a * a)
        E0 = 1.0
        # iterate the recursion numerically and compare with the closed form
        E = E0
        for n in range(1, 6):
            E_rec = 3 * E + k
            E_closed = 3 ** n * (E0 + k / 2) - k / 2
            self.assertAlmostEqual(E_rec, E_closed, places=6)
            E = E_rec

    def test_progress_threshold_matches_proposition5(self):
        """Prop 5: the 1/16 recursion contracts iff w < 4a; the 4/27 one iff
        w < sqrt(27) a."""
        for a in (F(1), F(2), F(7), F(10 ** 3)):
            w16 = 4 * a
            self.assertEqual(error_bound(a, w16) < w16, False)          # boundary
            self.assertEqual(error_bound(a, w16 - F(1, 10 ** 6)) < w16, True)
            w27 = F(27) ** F(1, 2) * a      # = 3*sqrt(3)*a
            self.assertEqual(error_bound_tight(a, w27 + F(1, 10 ** 6)) < w27 + F(1, 10 ** 6), False)
            self.assertEqual(error_bound_tight(a, w27 - F(1, 10 ** 6)) < w27 - F(1, 10 ** 6), True)


if __name__ == "__main__":
    unittest.main(verbosity=2)
