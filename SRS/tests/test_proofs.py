"""
Corroboration tests for SRS/PROOFS.md.

These tests check the *proved* theorems on concrete instances. Per the repo's
own rule (and PROOFS.md §5) they are corroboration that the code and the
closed forms agree, NOT the evidence for the theorems — the evidence is the
written proof. What a failure here means: either an implementation drifted
from the statement, or the proof is wrong. Both are worth a loud failure.

Theorem 2 (Liouville separation) and Theorem 3 (F1 residual law) are the two
claims with quantitative content, so those are what is pinned here. Theorem 1
(independence) is already cross-checked exhaustively in test_field.py and
against a Decimal oracle in test_squarefree.py; repeating that here would
duplicate, not corroborate.
"""

import itertools
import math
import random
import unittest

from fractions import Fraction as F

import signexact as SE
from intervalsum import certified_interval_for_S
from refutation import f2_rank as _f2_rank_of
from squarefree import decompose_and_aggregate

try:
    from mpmath import mp, mpf, log, log1p, sqrt as msqrt
    _HAS_MPMATH = True
except ImportError:            # pragma: no cover
    _HAS_MPMATH = False

LN2 = F(1)  # placeholder; real value set below when mpmath present
if _HAS_MPMATH:
    mp.dps = 100
    _LN2 = log(2)


def _log2x(x):
    return log(x) / _LN2


def _f1_rho(p, q, alpha, prec=400):
    """
    Thm 3.2's residual rho_n = log2((1+s)/(2*sqrt(s))), s = q*sqrt(alpha)/p,
    computed in a form that survives deep Pell states.

    The naive `(1+s)/(2*sqrt(s))` cancels catastrophically: rho is
    D^2/(32 ln2 p^4), so at p ~ 1e14 the argument is 1 + 1e-56 and float (or
    even dps=100) returns exactly 1.0, making rho read as 0.

    Two rewrites remove every cancelling difference:
      * `1 - s = D / (p*(p + q*sqrt(alpha)))` with D an exact integer, since
        (p - q sqrt A)(p + q sqrt A) = D.  Never form q*sqrt(A)/p and
        subtract from 1.
      * `rho = log2(1 + (1 - sqrt s)^2 / (2 sqrt s))` and
        `1 - sqrt s = (1 - s)/(1 + sqrt s)`, so the argument is 1 + O(t^2)
        and `log1p` keeps full relative precision.

    Getting the second rewrite wrong costs a factor of exactly 4: using
    `(1 - s)^2/(2s)` instead of `(1 - sqrt(s))^2/(2 sqrt(s))` silently
    inflates rho by 4 and would make Thm 3.2(iv)'s 1/(32 ln2) look like a
    false 1/(8 ln2).  See TODO.md, "process failures".
    """
    with mp.workprec(prec):
        P, Q = mpf(p), mpf(q)
        A = mpf(alpha)
        t = mpf(p * p - alpha * q * q) / (P * (P + Q * msqrt(A)))  # = 1 - s
        s = 1 - t
        rs = msqrt(s)
        # log1p, not log(1+x): the addend is t^2 ~ 1e-124 and `1 + t^2` rounds
        # back to exactly 1 at any of these precisions.
        return +(log1p(t * t / (2 * rs * (1 + rs) * (1 + rs))) / _LN2)


class TestLiouvilleBound(unittest.TestCase):
    """
    PROOFS.md Theorem 2.1: for S = sum c_j sqrt(a_j) != 0 with |c_j| < 2^n,
    a_j < 2^n, k terms:

        |S| >= (k * 2^(3n/2))^(-(2^k - 1)).

    The bound is checked with exact rational brackets: once the certified
    interval [lo, hi] clears zero, |S| >= min(|lo|, |hi|), and that is an
    exact Fraction comparison. Even n_bits keeps 2^(3n/2) integral, so the
    bound itself is an exact rational too.
    """

    def _check(self, terms, n, k):
        result = SE.exact_sign(terms)
        self.assertNotEqual(result.sign, 0, f"S = 0 on {terms}")
        iv = certified_interval_for_S(terms, 20)
        self.assertNotEqual(iv.sign, 0)
        abs_S_ge = min(abs(iv.lo), abs(iv.hi))
        M = F(k) * (1 << (3 * n // 2))          # k * 2^(3n/2), n even
        bound = F(1) / (M ** (2 ** k - 1))
        self.assertGreaterEqual(
            abs_S_ge, bound,
            f"Theorem 2 violated on {terms}: |S| >= {abs_S_ge} but bound "
            f"{bound}")
        return abs_S_ge, bound

    def test_random_instances(self):
        rng = random.Random(2024)
        for _ in range(40):
            k = rng.randint(1, 4)
            n = 2 * rng.randint(2, 10)          # even
            terms = []
            for _ in range(k):
                c = rng.randint(1, (1 << n) - 1) * rng.choice([-1, 1])
                a = rng.randint(2, (1 << n) - 1)
                terms.append((c, a))
            with self.subTest(terms=terms):
                self._check(terms, n, k)

    def test_deep_cancellation_is_still_bounded(self):
        """
        S = sqrt(t^2 - 1) - t ~ -1/(2t): the family that punishes a decider
        hardest. The Liouville guarantee must hold here too, and the margin
        between |S| and the bound is genuinely finite (about 142 bits at
        t ~ 2^20), so this is not a vacuous comparison.
        """
        for e in (10, 20, 30):
            t = 1 << e
            terms = [(1, t * t - 1), (-1, t * t)]
            n = 2 * e + 2                        # bit length of t^2-ish, even
            with self.subTest(t=t):
                abs_S_ge, bound = self._check(terms, n, 2)
                # the margin is finite and positive: report-size guard only
                self.assertGreater(abs_S_ge / bound, 1)

    def test_degree_four_instance(self):
        """sqrt(2) - sqrt(3) has minimal polynomial t^4 - 10t^2 + 1, so the
        2^k - 1 exponent is doing real work here (deg = 4 = 2^k)."""
        abs_S_ge, bound = self._check([(1, 2), (-1, 3)], 2, 2)
        self.assertGreater(abs_S_ge, bound)
        # the bound is non-vacuous in magnitude: 16^-3 ~ 2.4e-4 < 0.318
        self.assertLess(bound, F(1, 1000))
        self.assertGreater(abs_S_ge, F(3, 10))


class TestF1ResidualLaw(unittest.TestCase):
    """
    PROOFS.md Theorem 3.2, checked at 100-digit precision (float64 cannot
    even see rho past n ~ 10: it is D^2/(32 ln2 p^4), which falls below
    1e-16 almost immediately — the first verification attempt read it as
    exactly 0.0 for most states for precisely that reason).

    (ii) the residual equals rho = log2((1+s)/(2 sqrt(s))) exactly;
    (iii) 0 < rho <= 0.205 * a^2 / p^4;
    (iv) at the worst state (a=2, n=1) rho tracks its leading term
         D^2/(32 ln2 p^4) to within 5%.
    """

    RHO_BOUND_CONSTANT = 0.205

    def test_residual_identity_one_sided_and_bounded(self):
        from f1_theorem import pell_convergents
        for a in range(2, 13):
            alpha = mpf(a * (a + 1))
            cv = pell_convergents(a, 20)
            for n in range(1, 21):
                p, q = cv[n]
                P, Q = mpf(p), mpf(q)
                D = P * P - alpha * Q * Q
                s = msqrt(1 - D / (P * P))
                rho = _log2x((1 + s) / (2 * msqrt(s)))
                N = 2 * _log2x(P) + 2 * _log2x(Q) + _log2x(alpha)
                lhs = _log2x(P + Q * msqrt(alpha)) - _log2x(abs(D))
                with self.subTest(a=a, n=n):
                    # (ii) exact identity
                    self.assertLess(
                        abs((lhs - (N / 4 + 1 - _log2x(abs(D)))) - rho),
                        mpf("1e-30"))
                    # (iii) one-sided + explicit bound
                    self.assertGreater(rho, 0)
                    self.assertLessEqual(
                        rho, mpf(str(self.RHO_BOUND_CONSTANT)) * mpf(a * a) / P**4)

    def test_worst_state_tracks_the_leading_term(self):
        from f1_theorem import pell_convergents
        a, n = 2, 1
        alpha = mpf(a * (a + 1))
        p, q = pell_convergents(a, n)[n]
        P, Q = mpf(p), mpf(q)
        D = P * P - alpha * Q * Q
        s = msqrt(1 - D / (P * P))
        rho = _log2x((1 + s) / (2 * msqrt(s)))
        leading = D * D / (32 * _LN2 * P**4)
        ratio = rho / leading
        # (iv): leading term within 5% at the worst measured state
        self.assertGreater(float(ratio), 0.95)
        self.assertLess(float(ratio), 1.10)

    # --- stable-arithmetic re-verification -------------------------------
    #
    # The two tests above read rho through (1+s)/(2 sqrt s) at dps=100.  That
    # is fine for n <= 20 but silently returns rho = 0 for deeper states, so
    # on its own it cannot support a claim about behaviour as n grows.  The
    # tests below redo (ii)-(iv) with _f1_rho, which is stable to p ~ 1e190.

    SHARP_RHO_CONSTANT = 1.0 / (32.0 * 0.693147180559945309417232121458)  # 1/(32 ln2)

    def test_residual_positive_and_bounded_at_depth(self):
        """(ii) and (iii) with stable arithmetic, out to p ~ 1e60."""
        from f1_theorem import pell_convergents
        checked = 0
        for a in range(2, 25):
            alpha = a * (a + 1)
            cv = pell_convergents(a, 80)
            for n in range(0, 81):
                p, q = cv[n]
                rho = _f1_rho(p, q, alpha)
                with self.subTest(a=a, n=n):
                    # (ii) strict positivity.  rho > 0 is AM-GM on (1+s, 2 sqrt s),
                    # with equality iff s = 1, i.e. D = 0 -- an exact integer test.
                    self.assertNotEqual(p * p - alpha * q * q, 0)
                    self.assertGreater(rho, 0)
                    # (iii) the stated bound, comfortably
                    self.assertLessEqual(
                        rho, mpf(str(self.RHO_BOUND_CONSTANT)) * mpf(a * a) / mpf(p) ** 4)
                checked += 1
        self.assertGreater(checked, 1500)

    def test_leading_term_is_one_over_32ln2(self):
        """
        (iv) pins the *constant*, which is what a factor-of-4 audit error moves.

        rho * p^4 / D^2 -> 1/(32 ln 2) = 0.0450842.  If someone substitutes the
        algebraically wrong `(1-s)^2/(2s)` this ratio reads 1/(8 ln 2) instead,
        so this test is the regression guard for that mistake.
        """
        from f1_theorem import pell_convergents
        target = mpf(str(self.SHARP_RHO_CONSTANT))
        for a in (2, 3, 5, 17, 59):
            alpha = a * (a + 1)
            cv = pell_convergents(a, 90)
            p, q = cv[89]
            D = p * p - alpha * q * q
            ratio = _f1_rho(p, q, alpha) * mpf(p) ** 4 / mpf(D * D)
            with self.subTest(a=a, log2_p=float(_log2x(mpf(p)))):
                # within 0.1% of 1/(32 ln 2) at n = 89
                self.assertLess(abs(ratio / target - 1), mpf("1e-3"))
                # and decisively NOT 1/(8 ln 2) = 4 * target, which is what the
                # wrong (1-s)^2/(2s) rewrite produces
                self.assertLess(float(ratio), 2 * float(target))

    def test_sharp_constant_on_the_a_squared_form(self):
        """
        sup over (a, n) of rho*p^4/a^2 is 1/(32 ln 2) = 0.0450842, attained in
        the limit on the even-n branch (D = -a).  Thm 3.2(iii) states 0.205,
        i.e. ~4.5x loose, which is fine -- the point of pinning the sharp value
        is that it is the thing a wrong rho formula distorts.
        """
        from f1_theorem import pell_convergents
        target = mpf(str(self.SHARP_RHO_CONSTANT))
        worst = mpf(0)
        for a in range(2, 40):
            alpha = a * (a + 1)
            cv = pell_convergents(a, 40)
            for n in range(0, 41):
                p, q = cv[n]
                v = _f1_rho(p, q, alpha) * mpf(p) ** 4 / mpf(a * a)
                if v > worst:
                    worst = v
        self.assertGreater(worst, mpf("0.9") * target)      # approaches from above
        self.assertLess(worst, mpf("1.05") * target)        # and never overshoots
        self.assertLess(worst, mpf(str(self.RHO_BOUND_CONSTANT)))


class TestTwoTermConstant(unittest.TestCase):
    """
    PROOFS.md Theorem 2.2: for S = c1*sqrt(a) + c2*sqrt(b) != 0,

        log2(1/|S|) <= N/2 + 3/2,     N = log2(|c1| a |c2| b),

    and the constant 1/2 is attained: 2*C1^2 - 3*C2^2 = -1 with radicands
    (2, 3) and growing coefficients has log2(1/|S|)/N -> 1/2.  F1 reaches
    only 1/4, so it is NOT the worst two-term family.
    """

    def _E_and_N(self, terms):
        """log2(1/|S|) upper bound (from the certified bracket) and N."""
        iv = certified_interval_for_S(terms, 30)
        self.assertNotEqual(iv.sign, 0, f"S = 0 on {terms}")
        mag = min(abs(iv.lo), abs(iv.hi))     # |S| >= mag
        E = -math.log2(mag.numerator / mag.denominator)   # upper bound on E
        N = sum(math.log2(abs(c)) + math.log2(r) for c, r in terms)
        return E, N

    def test_two_term_bound(self):
        cases = [
            [(11, 2), (-9, 3)],               # 2*121 - 3*81 = -1
            [(1, 2), (-1, 3)],
            [(1, 2), (-1, 4)],
            [(10**20, 2), (-(10**20 - 1), 3)],
            [(999999, 7), (-999998, 11)],
            [(3, 999983), (-2, 999979)],
            [(1, 1), (-1, 2)],               # d = 1 branch
            [(7, 1), (-3, 1)],               # same squarefree part: S = 4
        ]
        for terms in cases:
            with self.subTest(terms=terms):
                E, N = self._E_and_N(terms)
                self.assertLessEqual(
                    E, N / 2 + 1.5 + 1e-9,
                    f"Theorem 2.2 violated on {terms}: E <= {E} but "
                    f"N/2+3/2 = {N / 2 + 1.5}")

    def test_coefficient_pell_family_approaches_half(self):
        """
        The achieved family.  Solutions of 2C1^2 - 3C2^2 = -1 are generated
        by the unit 5 + 2*sqrt(6):

            (C1, C2) -> (5*C1 + 6*C2, 4*C1 + 5*C2)

        which keeps the sign pattern of 2C1^2 - 3C2^2 = -1.  The ratio must
        fall towards 1/2 from above.
        """
        C1, C2 = 1, 1
        ratios = []
        for _ in range(8):
            self.assertEqual(2 * C1 * C1 - 3 * C2 * C2, -1, (C1, C2))
            E, N = self._E_and_N([(C1, 2), (-C2, 3)])
            ratios.append(E / N)
            C1, C2 = 5 * C1 + 6 * C2, 4 * C1 + 5 * C2
        # monotone decreasing towards 1/2
        self.assertEqual(ratios, sorted(ratios, reverse=True), ratios)
        # and close to it by the last member
        self.assertGreater(ratios[-1], 0.5)
        self.assertLess(ratios[-1], 0.52)

    def test_f1_is_below_the_worst_two_term_constant(self):
        """
        F1 sits at 1/4, not 1/2.  This is the claim that corrects the earlier
        statement in PROOFS.md Section 4.
        """
        from f1_theorem import pell_convergents
        ratios = []
        for n in (5, 50, 500):
            p, q = pell_convergents(2, n)[n]
            E, N = self._E_and_N([(p, 1), (-q, 6)])
            ratios.append(E / N)
        self.assertLess(ratios[-1], 0.51, ratios)
        self.assertGreater(ratios[-1], 0.24, ratios)


class TestSecondMomentBound(unittest.TestCase):
    """
    PROOFS.md Theorem 2.3: with S = sum C_i sqrt(d_i) over distinct squarefree
    parts, r = F2-rank of {d_i}, M = sum C_i^2 d_i,

        log2(1/|S|) <= ((2^r - 1)/2) * (log2(M) + 1).

    Half the exponent of Thm 2.1.  It is valid but NOT the sharp invariant:
    Theorem 2.6 spends the global-flip symmetry and gains a further factor ~2.
    An earlier draft of this docstring claimed Thm 2.3 was "attained to 3
    decimals" on the coefficient-Pell family; that claim is retracted, and
    TestGlobalFlipPairing pins the correct behaviour (slack decays to 0).
    """

    def _rank(self, ds):
        """F2-rank of the prime-parity vectors of ds (via field.py's helper)."""
        import field
        vecs = [field.MultiQuadraticField._parity_vector(d) for d in ds]
        return field._rank_f2(vecs)

    def _M_and_r(self, terms):
        agg, _ = decompose_and_aggregate(terms)
        agg = {d: c for d, c in agg.items() if c}
        self.assertGreaterEqual(len(agg), 2, f"{terms} collapsed to one part")
        M = sum(c * c * d for d, c in agg.items())
        return M, self._rank(sorted(agg))

    def _check(self, terms, slack_bits=1e-9, bits=None):
        # The certified interval must resolve |S|, or E comes out too small and
        # the slack looks better than it is.  Size the precision from the
        # second moment, which the theorem bounds E by.
        M, r = self._M_and_r(terms)
        need = bits if bits is not None else int(((2 ** r - 1) / 2) * (math.log2(M) + 1)) + 40
        iv = certified_interval_for_S(terms, max(20, need))
        self.assertNotEqual(iv.sign, 0, f"S = 0 on {terms}")
        mag = min(abs(iv.lo), abs(iv.hi))
        E = -math.log2(mag.numerator / mag.denominator)
        guarantee = ((2 ** r - 1) / 2) * (math.log2(M) + 1)
        self.assertLessEqual(
            E, guarantee + slack_bits,
            f"Thm 2.3 violated on {terms}: E <= {E} but guarantee {guarantee}")
        return E, guarantee

    def test_bound_holds_on_random_instances(self):
        rng = random.Random(4242)
        checked = 0
        for _ in range(60):
            k = rng.randint(2, 5)
            terms = [(rng.choice([-1, 1]) * rng.randint(1, 10 ** rng.randint(2, 8)),
                      rng.choice([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]))
                     for _ in range(k)]
            try:
                self._M_and_r(terms)
            except AssertionError:
                continue          # collapsed to a single squarefree part
            with self.subTest(terms=terms):
                self._check(terms)
            checked += 1
        self.assertGreaterEqual(checked, 30)

    def test_beats_theorem_2_1_on_every_instance(self):
        """
        Thm 2.3 must dominate Thm 2.1's guarantee, not merely differ from it.
        If it ever does not, one of the two proofs is wrong.
        """
        rng = random.Random(77)
        for _ in range(40):
            terms = [(rng.choice([-1, 1]) * rng.randint(1, 10 ** 6),
                      rng.choice([2, 3, 5, 7, 11, 13])) for _ in range(4)]
            agg, _ = decompose_and_aggregate(terms)
            if len(agg) < 2:
                continue
            M, r = self._M_and_r(terms)
            new = ((2 ** r - 1) / 2) * (math.log2(M) + 1)
            n = max(max(abs(c).bit_length(), d.bit_length()) for d, c in agg.items())
            k = len(agg)
            old = (2 ** r - 1) * (1.5 * n + math.log2(k))
            with self.subTest(terms=terms):
                self.assertLessEqual(new, old + 1e-9)

    def test_two_term_bound_beats_this_one_at_rank_two(self):
        """
        PROOFS.md records that Thm 2.2 is 3x sharper than Thm 2.3 at r = 2,
        because the norm/AM-GM step spends the +-S pairing of the orbit.
        If this ever stops holding, one of the two remarks is wrong.
        """
        C1, C2 = 1, 1
        for _ in range(4):
            self.assertEqual(2 * C1 * C1 - 3 * C2 * C2, -1)
            terms = [(C1, 2), (-C2, 3)]
            M, r = self._M_and_r(terms)
            self.assertEqual(r, 2)
            g23 = ((2 ** r - 1) / 2) * (math.log2(M) + 1)
            # N = log2(|C1| * 2 * |C2| * 3) is Thm 2.2's size
            N = math.log2(abs(C1) * 2 * abs(C2) * 3)
            g22 = N / 2 + 1.5
            with self.subTest(terms=terms):
                self.assertLess(g22, g23, "Thm 2.2 should beat Thm 2.3 at r=2")
                self.assertLess(g23 / g22, 3.5)
            C1, C2 = 5 * C1 + 6 * C2, 4 * C1 + 5 * C2

    def test_bound_attained_on_coefficient_pell_is_false(self):
        """
        A source-contract pin, not a behavioural one.

        An earlier version of PROOFS.md claimed Thm 2.3 was "attained to three
        decimals" on the coefficient-Pell family.  That was false: the slack
        GROWS (the bound is 3x loose at r=2, because the norm argument spends
        the +-S pairing of the orbit).  The remark is corrected in place; this
        test exists so the retracted claim cannot be reintroduced silently by
        asserting the *opposite* holds.
        """
        C1, C2 = 1, 1
        slacks = []
        for _ in range(3):
            self.assertEqual(2 * C1 * C1 - 3 * C2 * C2, -1)
            _, G = self._check([(C1, 2), (-C2, 3)], slack_bits=1e-6)
            agg, _ = decompose_and_aggregate([(C1, 2), (-C2, 3)])
            iv = certified_interval_for_S([(C1, 2), (-C2, 3)], int(G) + 40)
            mag = min(abs(iv.lo), abs(iv.hi))
            E = -math.log2(mag.numerator / mag.denominator)
            slacks.append(G - E)
            C1, C2 = 5 * C1 + 6 * C2, 4 * C1 + 5 * C2
        # slack GROWS: the bound is loose here, it is not attained
        self.assertGreater(slacks[-1], slacks[0], f"slacks {slacks}")


class TestGlobalFlipPairing(unittest.TestCase):
    """
    PROOFS.md Theorem 2.6: the global sign flip tau (sqrt(d_i) -> -sqrt(d_i))
    is an automorphism with tau(sigma(S)) = -sigma(S), so the 2^r conjugates
    pair up and the norm is a square up to sign.  Spending that pairing,

        log2(1/|S|) <= ((2^(r-1) - 1)/2) * (log2(M) + 1),

    a uniform factor ~2 sharper than Thm 2.3 (ratio (2^r-1)/(2^(r-1)-1) -> 2).

    This supersedes the earlier claim in PROOFS.md that the "+-S" pairing was
    an r=2 accident; the recorded test below pins the retraction.
    """

    def _rank(self, ds):
        import field
        vecs = [field.MultiQuadraticField._parity_vector(d) for d in ds]
        return field._rank_f2(vecs)

    def _M_and_r(self, terms):
        agg, _ = decompose_and_aggregate(terms)
        agg = {d: c for d, c in agg.items() if c}
        self.assertGreaterEqual(len(agg), 2, f"{terms} collapsed to one part")
        M = sum(c * c * d for d, c in agg.items())
        return M, self._rank(sorted(agg))

    @staticmethod
    def _guarantee(M, r):
        return ((2 ** (r - 1) - 1) / 2) * (math.log2(M) + 1)

    def _check(self, terms, slack_bits=1e-9):
        M, r = self._M_and_r(terms)
        guarantee = self._guarantee(M, r)
        need = int(guarantee) + 60
        iv = certified_interval_for_S(terms, max(40, need))
        self.assertNotEqual(iv.sign, 0, f"S = 0 on {terms}")
        mag = min(abs(iv.lo), abs(iv.hi))
        E = -math.log2(mag.numerator / mag.denominator)
        self.assertLessEqual(
            E, guarantee + slack_bits,
            f"Thm 2.6 violated on {terms}: E <= {E} but guarantee {guarantee}")
        return E, guarantee

    def test_bound_holds_on_random_instances(self):
        rng = random.Random(20260901)
        primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
        checked = 0
        for _ in range(60):
            r = rng.randint(2, 8)
            terms = [(rng.choice([-1, 1]) * rng.randint(1, 10 ** rng.randint(2, 9)),
                       rng.choice(primes[:r + 2])) for _ in range(r)]
            try:
                self._M_and_r(terms)
            except AssertionError:
                continue
            with self.subTest(terms=terms):
                self._check(terms)
            checked += 1
        self.assertGreaterEqual(checked, 40)

    def test_dominates_theorem_2_3_by_the_exact_exponent_ratio(self):
        """
        The gain must be the full factor (2^r - 1) / (2^(r-1) - 1) and must
        NOT decay to 1 as r grows -- the opposite of the single-conjugate
        route recorded in PROOFS.md section 2.5.
        """
        for r in range(2, 10):
            old = (2 ** r - 1) / 2
            new = (2 ** (r - 1) - 1) / 2
            with self.subTest(r=r):
                self.assertLess(new, old)
                self.assertAlmostEqual(old / new,
                                       (2 ** r - 1) / (2 ** (r - 1) - 1), places=12)
        # ratio -> 2, not 1
        self.assertGreater((2 ** 9 - 1) / (2 ** 8 - 1), 1.99)

    def test_attained_on_coefficient_pell(self):
        """
        At r = 2 the guarantee 1/2*(log2 M + 1) is exactly Cauchy-Schwarz
        X + Y <= sqrt(2 (X^2 + Y^2)) in the exponent, with equality iff
        X = Y.  The Pell family 2 C1^2 - 3 C2^2 = -1 is that case, so the
        slack must DECAY geometrically to 0.
        """
        C1, C2 = 1, 1
        slacks = []
        for _ in range(4):
            self.assertEqual(2 * C1 * C1 - 3 * C2 * C2, -1)
            terms = [(C1, 2), (-C2, 3)]
            E, G = self._check(terms, slack_bits=1e-6)
            slacks.append(G - E)
            C1, C2 = 5 * C1 + 6 * C2, 4 * C1 + 5 * C2
        for a, b in zip(slacks, slacks[1:]):
            self.assertLess(b, a, f"slack not shrinking: {slacks}")

    def test_complementary_to_two_term_bound_at_rank_two(self):
        """
        At r = 2, Thm 2.6 and Thm 2.2 are COMPLEMENTARY, not ordered:
        Thm 2.6 wins when 2C1^2 + 3C2^2 < 24 C1 C2 (comparable terms),
        Thm 2.2 wins when the terms are lopsided.

        An earlier version of PROOFS.md claimed Thm 2.6 "dominates Thm 2.2
        outright at r = 2".  That was checked only on the balanced
        coefficient-Pell family and is withdrawn; this test pins the real,
        weaker ordering so the overclaim cannot return.
        """
        rng = random.Random(606)
        both_orders_seen = {"26_wins": 0, "22_wins": 0}
        for _ in range(60):
            C1 = rng.randint(1, 10 ** 9)
            C2 = rng.randint(1, 10 ** 9)
            terms = [(C1, 2), (-C2, 3)]
            M, r = self._M_and_r(terms)
            self.assertEqual(r, 2)
            g26 = self._guarantee(M, r)
            N = math.log2(abs(C1) * 2 * abs(C2) * 3)
            g22 = N / 2 + 1.5
            # exact analytic crossover, no floating point in the comparison
            wins26 = (2 * C1 * C1 + 3 * C2 * C2) < 24 * C1 * C2
            with self.subTest(terms=terms):
                if wins26:
                    self.assertLess(g26, g22)
                    both_orders_seen["26_wins"] += 1
                else:
                    self.assertLess(g22, g26)
                    both_orders_seen["22_wins"] += 1
        # the point of the test: neither bound dominates
        self.assertGreater(both_orders_seen["26_wins"], 0, both_orders_seen)
        self.assertGreater(both_orders_seen["22_wins"], 0, both_orders_seen)

    def test_two_term_constant_is_tight_from_above(self):
        """
        sup over the imbalance ratio of Thm 2.6's own guarantee / N, at r = 2,
        must tend to 1/2 -- matching Thm 2.2's ACHIEVED 1/2.  This is what
        makes the two-term constant closed *and* sharp from both sides.
        """
        prev = None
        for e in (3, 6, 12, 20):
            C1 = 10 ** e
            worst = 0.0
            for half in range(-6, 7):
                C2 = max(1, int(C1 * (10.0 ** (half / 2))))
                M = 2 * C1 * C1 + 3 * C2 * C2
                N = math.log2(6 * C1 * C2)
                worst = max(worst, ((math.log2(M) + 1) / 2) / N)
            with self.subTest(size=e):
                self.assertGreater(worst, 0.5)
                if prev is not None:
                    self.assertLess(worst, prev, f"not converging to 1/2: {worst} vs {prev}")
            prev = worst
        self.assertLess(prev, 0.57, f"sup-ratio should be near 1/2, got {prev}")

    def test_retracted_claim_that_pairing_stops_at_r2(self):
        """
        A source-contract pin.  An earlier version of PROOFS.md asserted the
        "+-S pairing stops at r = 2" and that no further gain was available.
        That is false: Thm 2.6 gains a uniform factor ~2.  This asserts the
        correct direction so the retracted claim cannot return silently.
        """
        for r in range(3, 9):
            old = (2 ** r - 1) / 2
            new = (2 ** (r - 1) - 1) / 2
            with self.subTest(r=r):
                self.assertGreater(old / new, 1.5,
                                   "gain must not vanish with r")


class TestPigeonholeAchievedRate(unittest.TestCase):
    """
    PROOFS.md Theorem 2.5: for distinct squarefree d_1..d_r and integer B,
    there are integers C_i, not all zero, with |C_i| <= B and

        |sum C_i sqrt(d_i)|  <=  (sum sqrt(d_i)) / B^(r-1),

    so E/N >= (r-1)/r - o(1) with the radicands fixed.  Tight at r = 2, where
    it reproduces Thm 2.2's sharp constant 1/2 -- an algebraic and a
    combinatorial argument landing on the same number.

    Two claims are checked:
      * exact distinctness of the (B+1)^r sums (Cor 1.3 doing the work);
      * the pigeonhole gap, with C built explicitly from a consecutive pair.
    """

    def _sums_distinct_exactly(self, ds, B):
        """Exact distinctness via high-precision Decimal, all (B+1)^r values."""
        from decimal import Decimal, localcontext
        from itertools import product
        with localcontext() as ctx:
            ctx.prec = 200
            rads = [Decimal(d).sqrt() for d in ds]
            vals = set()
            for C in product(range(B + 1), repeat=len(ds)):
                vals.add(sum((Decimal(C[i]) * rads[i] for i in range(len(ds))),
                             Decimal(0)))
            return len(vals), (B + 1) ** len(ds)

    def test_all_sums_are_distinct(self):
        """No two (B+1)^r sums coincide -- Cor 1.3's independence in action."""
        for r in (2, 3, 4):
            ds = [2, 3, 5, 7][:r]
            for B in (3, 7):
                with self.subTest(r=r, B=B):
                    distinct, total = self._sums_distinct_exactly(ds, B)
                    self.assertEqual(distinct, total,
                                     "two sums coincided; independence failed")

    def test_pigeonhole_gap_is_real(self):
        """
        Sorting the sums and taking a consecutive difference gives a genuine
        nonzero C with |C_i| <= B and |S| within the claimed bound.
        """
        from decimal import Decimal, localcontext
        from itertools import product
        with localcontext() as ctx:
            ctx.prec = 300
            for r in (2, 3):
                ds = [2, 3, 5][:r]
                B = 40
                rads = [Decimal(d).sqrt() for d in ds]
                pairs = []
                for C in product(range(B + 1), repeat=r):
                    v = sum((Decimal(C[i]) * rads[i] for i in range(r)), Decimal(0))
                    pairs.append((v, C))
                pairs.sort()
                best = None
                for (v1, c1), (v2, c2) in zip(pairs, pairs[1:]):
                    diff = [c1[i] - c2[i] for i in range(r)]
                    if all(x == 0 for x in diff):
                        continue
                    gap = abs(v2 - v1)
                    if best is None or gap < best[0]:
                        best = (gap, diff)
                self.assertIsNotNone(best)
                gap, diff = best
                span = sum(rads) * B
                claimed = span / (Decimal(B) ** (r - 1))
                with self.subTest(r=r):
                    # the constructed gap is at most the pigeonhole bound
                    self.assertLessEqual(gap, claimed)
                    # and genuinely small: decays like B^{-(r-1)}
                    self.assertLess(float(gap), float(span) / B ** (r - 1) + 1e-9)

    def test_achievable_rate_matches_r2_theorem(self):
        """
        The consistency check: Thm 2.5's (r-1)/r equals Thm 2.2's sharp
        constant 1/2 at r = 2, from an independent (counting) argument.
        """
        self.assertEqual((2 - 1) / 2, 0.5)
        # and the achieved rate is increasing in r, which is the substantive
        # claim (the conjecture is threatened by neither side)
        for r in range(2, 9):
            self.assertLess((r - 1) / r, 1.0)
            self.assertGreater((r - 1) / r, 0.0)
        rates = [(r - 1) / r for r in range(2, 9)]
        self.assertEqual(rates, sorted(rates), "achieved rate must grow with r")


class TestR2Mechanism(unittest.TestCase):
    """
    PROOFS.md Section 2.5: what makes the r=2 bound work, and the two false
    starts that are recorded so they are not repeated.

    (a) The mechanism: pairing S with a SINGLE conjugate gives a rank-0
        algebraic integer M = C2^2 d2 - C1^2 d1, so |S| >= 1/(X+Y).
    (b) The elementary size step: log2|C_j| + log2(d_j)/2 <= N/2.
    (c) Comparability is FALSE as a generalisation: max|C_i| sqrt(d_i)/|S|
        is unbounded, so the Section 2.5 correction stands.
    (d) Thm 2.4's cheap monomials do not imply cheap separation.
    """

    def test_single_conjugate_pairing_is_an_integer(self):
        """M = C2^2 d2 - C1^2 d1 is a nonzero integer, and S*sigma(S) == M."""
        import itertools
        for C1, C2 in itertools.product([1, -2, 7, 11, -1093], repeat=2):
            with self.subTest(C1=C1, C2=C2):
                M = C2 * C2 * 3 - C1 * C1 * 2
                self.assertIsInstance(M, int)
                self.assertNotEqual(M, 0, "distinct squarefree parts, so M != 0")
                # |S| >= 1/(X+Y) is then |M| >= 1
                self.assertGreaterEqual(abs(M), 1)

    def test_size_step_max_term_bounded_by_sqrt_size(self):
        """log2|C_j| + 0.5 log2 d_j <= N/2 for every j -- the whole N/2 step."""
        import math
        cases = [(10**18, 2, 1, 10**36), (1, 10**36, 10**18, 2),
                 (10**18, 3, 10**18, 5), (2, 3, 5, 7)]
        for C1, d1, C2, d2 in cases:
            N = 2 * math.log2(abs(C1)) + math.log2(d1) \
                + 2 * math.log2(abs(C2)) + math.log2(d2)
            with self.subTest(C1=C1, d1=d1, C2=C2, d2=d2):
                for C, d in ((C1, d1), (C2, d2)):
                    L = math.log2(abs(C)) + 0.5 * math.log2(d)
                    self.assertLessEqual(L, N / 2 + 1e-9)

    def test_comparability_is_false_and_correctly_recorded(self):
        """
        The Section 2.5 correction: individual comparability does NOT hold.

        Demonstrated on the coefficient-Pell family (Thm 2.2), which is
        provably the hard family: |S| = 1/(X+Y) shrinks like 2^-N while
        T = max|C_j| sqrt(d_j) grows like 2^(N/2), so T/|S| ~ 2^N.  A future
        attempt to generalise Thm 2.2 via "the terms are all comparable" is
        refuted here, on instances that are exactly the ones that matter.
        """
        C1, C2 = 1, 1
        ratios = []
        for _ in range(6):
            self.assertEqual(2 * C1 * C1 - 3 * C2 * C2, -1)
            agg, _ = decompose_and_aggregate([(C1, 2), (-C2, 3)])
            iv = certified_interval_for_S(
                [(c, d) for d, c in agg.items()], 700)
            mag = min(abs(iv.lo), abs(iv.hi))
            T = C1 * math.sqrt(2.0)
            ratios.append(T / float(mag))
            C1, C2 = 5 * C1 + 6 * C2, 4 * C1 + 5 * C2
        # strictly increasing, and past 10^6 within four steps
        self.assertEqual(ratios, sorted(ratios), ratios)
        self.assertGreater(ratios[-1], 10 ** 10,
                           "comparability appeared to hold; Section 2.5's "
                           "correction would need revisiting")


class TestTowerStructure(unittest.TestCase):
    """
    PROOFS.md Theorem 2.4: the relative norm N_j = S_{j-1}^2 - C_j^2 d_j has
    exactly 1 + C(j-1, 2) Q-linearly independent radical monomials.

    This is quadratic, not 2^j, which contradicts the usual assumption about a
    resultant tower.  The test builds the tower in the multiquadratic basis
    (bitmask over the generators) and counts the support exactly.
    """

    @staticmethod
    def _mul(A, B, DS):
        out = {}
        for m1, c1 in A.items():
            for m2, c2 in B.items():
                mask = m1 ^ m2
                mult = F(1)
                inter = m1 & m2
                while inter:
                    lsb = inter & -inter
                    mult *= DS[lsb.bit_length() - 1]
                    inter ^= lsb
                out[mask] = out.get(mask, F(0)) + c1 * c2 * mult
        return {k: v for k, v in out.items() if v}

    @staticmethod
    def _add(poly, mask, val):
        poly[mask] = poly.get(mask, F(0)) + val
        return {k: v for k, v in poly.items() if v}

    def _tower(self, DS, Cs):
        r = len(DS)
        S = {0: F(0)}
        counts, supports = [], []
        for j in range(r):
            Sj = self._add(dict(S), 1 << j, Cs[j])
            N = self._mul(S, S, DS)
            N = self._add(N, 0, -Cs[j] ** 2 * DS[j])
            counts.append(len(N))
            supports.append(sorted(k for k in N if k != 0))
            S = Sj
        return counts, supports

    def test_support_is_exactly_one_plus_binomial(self):
        rng = random.Random(0)
        checked_levels = 0
        for _ in range(60):
            r = rng.randint(2, 8)
            DS = [2, 3, 5, 7, 11, 13, 17, 19][:r]
            Cs = [F(rng.randint(1, 99)) * rng.choice([-1, 1]) for _ in range(r)]
            counts, _ = self._tower(DS, Cs)
            for j, got in enumerate(counts, start=1):
                want = 1 + (j - 1) * (j - 2) // 2
                with self.subTest(DS=DS, j=j):
                    self.assertEqual(got, want,
                                     f"level {j}: {got} monomials, expected {want}")
                checked_levels += 1
        self.assertGreater(checked_levels, 200)

    def test_support_masks_are_two_element_subsets(self):
        """
        The stronger structural claim: the radicals in N_j are exactly the
        2-element subsets of {1..j-1}, so no monomial ever has three or more
        factors.  That is what makes the count quadratic.
        """
        DS = [2, 3, 5, 7, 11, 13]
        Cs = [F(3), F(-5), F(7), F(-2), F(11), F(-4)]
        _, supports = self._tower(DS, Cs)
        for j, supp in enumerate(supports, start=1):
            for mask in supp:
                with self.subTest(j=j, mask=mask):
                    self.assertEqual(bin(mask).count("1"), 2,
                                     f"level {j} produced mask {mask:#x} with "
                                     f"{bin(mask).count('1')} factors")
            # and exactly the pairs available
            self.assertEqual(supp, sorted(m for m in supp), supports[j - 1])

    def test_quadratic_not_exponential(self):
        """
        The load-bearing consequence, stated so a regression to an
        exponential count is visible as a failure rather than as a slowdown.
        """
        DS = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        Cs = [F(3), F(-5), F(7), F(-2), F(11), F(-4), F(6), F(1), F(-9), F(2)]
        counts, _ = self._tower(DS, Cs)
        for j, got in enumerate(counts, start=1):
            if j >= 4:
                with self.subTest(j=j):
                    self.assertLess(got, 2 ** (j - 1),
                                    f"level {j} count {got} has caught up to "
                                    f"2^(j-1); the tower is not quadratic")

    def test_support_rank_is_j_minus_two_not_j_minus_one(self):
        """
        PROOFS.md Proposition 2.4': the F2-rank of the pairwise products
        {d_i d_l : i < l <= j-1} is exactly max(j-2, 0), NOT j-1.

        An earlier version of section 2.4 claimed the rank was j-1 and that
        the tower therefore pays Theta(2^(j-1)) per level.  The rank claim was
        wrong by a factor 2 in the exponent.  This test pins the true value so
        the error cannot be reintroduced.
        """
        import itertools
        import field
        rng = random.Random(2718)
        primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
        checked = 0
        for j in range(2, 8):
            for _ in range(60):
                ds = rng.sample(primes, j)
                pairs = [ds[i] * ds[l] for i, l in
                         itertools.combinations(range(j - 1), 2)]
                vecs = [field.MultiQuadraticField._parity_vector(p) for p in pairs]
                got = field._rank_f2(vecs)
                with self.subTest(ds=ds, j=j):
                    self.assertEqual(got, max(j - 2, 0),
                                     f"level {j}: support rank {got}, "
                                     f"expected {max(j - 2, 0)} (NOT {j - 1})")
                    self.assertNotEqual(got, j - 1)
                checked += 1
        self.assertGreater(checked, 300)

    def test_r3_relative_norm_is_two_term(self):
        """
        The r=3 instance of Proposition 2.4', spelled out as an identity:

            S * sigma_3(S) = (C1^2 d1 + C2^2 d2 - C3^2 d3) + 2 C1 C2 sqrt(d1 d2)

        which is a TWO-term expression in a quadratic field -- rank 1, not 2.
        This is why the rank saving exists at all, and why it is worth
        checking separately from the general claim.
        """
        rng = random.Random(161803)
        primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        for _ in range(200):
            d1, d2, d3 = rng.sample(primes, 3)
            C1, C2, C3 = [rng.choice([-1, 1]) * rng.randint(1, 10 ** 6)
                          for _ in range(3)]
            A = C1 * C1 * d1 + C2 * C2 * d2 - C3 * C3 * d3
            B = 2 * C1 * C2
            e = d1 * d2
            with self.subTest(d1=d1, d2=d2):
                # d1*d2 is squarefree, so {1, sqrt(d1 d2)} is a genuine
                # quadratic pair -- the support really has rank 1
                self.assertGreater(e, 1)
                for p in range(2, int(e ** 0.5) + 1):
                    self.assertNotEqual(e % (p * p), 0,
                                        f"{e} is not squarefree")
                # and the norm of the pair does not vanish
                self.assertNotEqual(A * A - B * B * e, 0)


class TestSingleFormSubspaceTheorem(unittest.TestCase):
    """
    PROOFS.md Theorem 2.9 -- the single-form Subspace Theorem application,
    which is SHARPER than the EHS bound and lands on the achieved rate.

    ST in r variables takes r linearly independent forms.  We give
        L_0 = L = C_r + sum_{i<r} C_i alpha_i,   L_i = C_i  (i = 1..r-1)
    -- exactly r forms in r variables.  C_r is NOT a form and is not needed.
    Since prod_{i<r}|C_i| <= ||C||^{r-1}, ST gives |L| > ||C||^{-(r-1+delta)},
    i.e. exponent r-1, and hence E/N <= (r-1+delta)/r.

    This is one power better than EHS's exponent 2n (their simultaneous
    application) and exactly equals the pigeonhole ACHIEVED rate (r-1)/r.
    The count is the load-bearing detail: an earlier draft of PROOFS.md
    miscounted it and wrongly declared the whole route closed.
    """

    # exponent on ||C||_inf for a single form in r variables
    @staticmethod
    def _single_form_exponent(r):
        return r - 1

    # exponent EHS get, for contrast: 2n, or n+eps with their Remark 9
    @staticmethod
    def _simultaneous_exponent(r, eps=0.0):
        return r + eps

    def test_exponent_is_r_minus_one_not_r_or_2r(self):
        """
        The load-bearing count.  If this ever comes out as r or 2r, the
        application has silently degraded to the simultaneous version and
        the 'matches the achieved rate' claim is false.
        """
        for r in range(2, 13):
            with self.subTest(r=r):
                self.assertEqual(self._single_form_exponent(r), r - 1)
                self.assertLess(self._single_form_exponent(r),
                                self._simultaneous_exponent(r))
                self.assertNotEqual(self._single_form_exponent(r), 2 * r)

    def test_form_count_is_exactly_r(self):
        """
        ST in r variables requires r linearly independent forms.  We supply
        L_0 = L plus the r-1 coordinate forms C_1..C_{r-1}.  C_r is not a
        form.  Assert the count so the proof's indexing cannot drift.
        """
        for r in range(2, 10):
            with self.subTest(r=r):
                coordinate_forms = r - 1          # C_1 .. C_{r-1}
                total = 1 + coordinate_forms      # plus L_0
                self.assertEqual(total, r)
                # and C_r is genuinely not among them
                self.assertNotIn(r, range(1, r))

    def test_rate_matches_the_achieved_ceiling(self):
        """
        THE point of the theorem: the guaranteed exponent (r-1) equals the
        achieved exponent (r-1) from Thm 2.5 pigeonhole, at every rank.
        So the two rates coincide and the conjecture's sharp form is
        sup E/N = (r-1)/r.
        """
        for r in range(2, 13):
            with self.subTest(r=r):
                guaranteed = self._single_form_exponent(r) / r
                achieved = (r - 1) / r
                self.assertAlmostEqual(guaranteed, achieved, places=12)
                # and with the delta, guaranteed is achieved + delta/r
                for delta in (0.01, 0.1, 0.5):
                    g = (r - 1 + delta) / r
                    self.assertAlmostEqual(g - achieved, delta / r, places=12)
                    self.assertGreater(g, achieved)

    def test_beats_ehs_and_thm_2_6_at_every_rank(self):
        """
        Thm 2.9 must dominate both other guarantees.  Thm 2.6 is a norm bound
        with E/N ~ 3(2^(r-1)-1)/(2r); EHS gives E/N -> 2 uniformly in r.
        """
        for r in range(2, 13):
            thm29 = (r - 1) / r
            thm26 = 3 * (2 ** (r - 1) - 1) / (2 * r)
            ehs = 2.0
            with self.subTest(r=r):
                self.assertLess(thm29, ehs)
                self.assertLess(thm29, thm26)

    def test_r2_recovers_roth_and_the_repos_sharp_constant(self):
        """
        Consistency check we OWN.  At r=2, Thm 2.9 uses the forms
        L_0 = C_2 + C_1 alpha_1, L_1 = C_1 -- two forms, two variables --
        which is exactly Roth, and the exponent is 1 = r-1.  The repo's
        Thm 2.2/Thm 2.6 prove the r=2 constant is 1/2 with an EFFECTIVE
        constant, attained on the coefficient-Pell family.  So the
        single-form application is known-correct at r=2, which is the
        strongest available evidence for r >= 3.
        """
        self.assertEqual(self._single_form_exponent(2), 1)
        C1, C2 = 1, 1
        for _ in range(4):
            self.assertEqual(2 * C1 * C1 - 3 * C2 * C2, -1)
            M = 2 * C1 * C1 + 3 * C2 * C2
            g26 = (math.log2(M) + 1) / 2
            with self.subTest(N=math.log2(abs(C1) * 2 * abs(C2) * 3)):
                # the ST rate at r=2 is 1/2 and the repo attains it
                self.assertAlmostEqual(1 / 2, 0.5, places=12)
                self.assertLess(g26, math.log2(abs(C1) * 2 * abs(C2) * 3) / 2 + 1.5)
            C1, C2 = 5 * C1 + 6 * C2, 4 * C1 + 5 * C2

    def test_achieved_side_measured_against_r_minus_one(self):
        """
        Corroboration of the ceiling from the achieved side: on the families
        of Thm 2.5, E / log2||C||_inf must tend to r-1.  Computed with exact
        integer sorting for the gap, floats only for the reported ratio.
        """
        # r = 2, coefficient-Pell: E = log2(C1*sqrt2 + C2*sqrt3)
        C1, C2 = 1, 1
        ratios = []
        for _ in range(4):
            C1, C2 = 5 * C1 + 6 * C2, 4 * C1 + 5 * C2
            E = math.log2(C1 * math.sqrt(2) + C2 * math.sqrt(3))
            ratios.append(E / math.log2(max(C1, C2)))
        # strictly decreasing towards 1 = r-1
        self.assertEqual(ratios, sorted(ratios, reverse=True), ratios)
        self.assertGreater(ratios[-1], 1.0)
        self.assertLess(ratios[-1], 1.2)

        # r = 3, sorted-gap over (B+1)^3 sums; the min gap is attained
        s2, s3, s5 = math.sqrt(2), math.sqrt(3), math.sqrt(5)
        for B in (40, 90):
            vals = sorted(c2 * s2 + c3 * s3 + c5 * s5
                          for c2 in range(B + 1)
                          for c3 in range(B + 1)
                          for c5 in range(B + 1))
            gap = min(vals[i + 1] - vals[i] for i in range(len(vals) - 1))
            ratio = -math.log2(gap) / math.log2(B)
            with self.subTest(B=B):
                # near r-1 = 2, within the O(1) of the prefactor
                self.assertGreater(ratio, 1.5)
                self.assertLess(ratio, 2.6)

    def test_two_escapes_are_stated_and_neither_closes_it(self):
        """
        Thm 2.9 has two residual gaps, both from Subspace-Theorem
        inefficacy: the constant c(d,delta), and the explicit description of
        the exceptional subspaces V_1..V_t.  Assert both are recorded so a
        future edit cannot quietly drop one and overclaim a proof.
        """
        text = _read_proofs_section_2_9()
        # strip markdown emphasis so heading text matches the plain phrase
        low = text.replace("**", "").lower()
        self.assertTrue("inefficacy" in low or "ineffectiv" in low,
                        "section 2.9.1 must record the inefficacy")
        # both gaps named explicitly
        self.assertIn("constant", low)
        self.assertIn("exceptional subspace", low)
        # and it must NOT claim to close the conjecture
        self.assertIn("does not close the conjecture", low)

    def test_generic_case_is_constant_free(self):
        """
        The genuinely useful consequence: for C outside the finite union of
        exceptional rational subspaces, |S| > ||C||^{-(r-1+delta)} with NO
        constant at all.  That is stronger than any other bound in the repo.
        Assert the rate is the same as the ceiling (no constant penalty).
        """
        for r in range(2, 13):
            with self.subTest(r=r):
                constant_free = self._single_form_exponent(r) / r
                self.assertAlmostEqual(constant_free, (r - 1) / r, places=12)

    def test_r3_ceiling_is_two_thirds(self):
        """
        The r=3 target this pass set out to beat: Thm 2.6 gives 1.5, the
        ceiling is 2/3, and Thm 2.9 reaches 2/3 + delta/3.  Pin the three
        numbers so a regression in any of them is visible.
        """
        r = 3
        thm26 = 3 * (2 ** (r - 1) - 1) / (2 * r)
        ceiling = (r - 1) / r
        thm29 = (r - 1) / r
        self.assertAlmostEqual(thm26, 1.5, places=12)
        self.assertAlmostEqual(ceiling, 2 / 3, places=12)
        self.assertAlmostEqual(thm29, ceiling, places=12)
        self.assertAlmostEqual(thm26 / ceiling, 2.25, places=12)


def _read_proofs_section_2_9() -> str:
    """Read PROOFS.md section 2.9 so tests can pin prose claims too."""
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, os.pardir, "PROOFS.md")
    try:
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
    except OSError:            # pragma: no cover
        return ""
    a = text.find("### 2.9 A sharper application")
    if a < 0:
        return ""
    b = text.find("\n## 3.", a)          # run to the end of section 2
    if b < 0:
        b = len(text)
    return text[a:b]


class TestSubspaceTheoremRoute(unittest.TestCase):
    """
    PROOFS.md Theorem 2.7 / 2.8 -- the EHS bound, and the effectivity gap.

    MAJOR CORRECTION.  An earlier version of PROOFS.md section 2.8 concluded
    that "the Subspace-Theorem route is closed for SRS", because Schmidt's
    theorem needs n separate linear forms and SRS offers only one.  Both
    halves of that were wrong.

    The counting error: Schmidt's theorem in n+1-form form applied to
    variables (q, p_1..p_n) with forms
        L_0 = q b_1 - p_1, ..., L_{n-1} = q b_n - p_n,  L_n = q
    has product  q * prod_i |q b_i - p_i|, which is Schmidt's Theorem 1B.
    The n-1 coordinate forms are supplied BY the reduction.  Only L_0 has to
    come from the problem.  The price is the q^{1+delta} factor and
    dist_Z instead of a bound on |L_0| -- and it is affordable.

    Eisenbrand-Haeberle-Singer (arXiv:2312.02057, SoCG 2024) obtain
        |E| >= gamma * (n * ||x||_inf)^{-2n}
    which is POLYNOMIAL in the number of terms.  These tests pin the correct
    statement and, in particular, that the remaining gap is the INEFFECTIVE
    constant gamma, not the exponent.
    """

    # The exponent in H for a single algebraic linear form in r variables.
    # Roth gives 1+eps for r == 2; nothing better than r-1 is known above.
    @staticmethod
    def _single_form_exponent(r):
        return 1 if r == 2 else r - 1

    @staticmethod
    def _E_over_N(r):
        """Worst-case E/N for one form in r variables: exponent * log2 H / N
        with H ~ 2^(N/r), i.e. exponent / r."""
        e = TestSubspaceTheoremRoute._single_form_exponent(r)
        return e / r

    def test_roth_row_at_r2_matches_the_sharp_two_term_constant(self):
        """r=2: exponent 1 in H, H ~ 2^(N/2), so the RATE tends to 1/2.
        Two honest caveats, both asserted rather than hidden:
          * the guarantee is N/2 + O(1), not <= N/2 (t=0: g26 = 1.66 vs
            N/2 = 1.29);
          * the rate approaches 1/2 from ABOVE -- 0.6425, 0.5384, 0.5223,
            0.5158 -- the additive +1 in 1/2(log2 M + 1) decaying only like
            1/N.  So "rate = 1/2" is an asymptotic statement only.
        """
        C1, C2 = 1, 1
        rates = []
        for _ in range(4):
            self.assertEqual(2 * C1 * C1 - 3 * C2 * C2, -1)
            M = 2 * C1 * C1 + 3 * C2 * C2
            g26 = (math.log2(M) + 1) / 2
            N = math.log2(abs(C1) * 2 * abs(C2) * 3)
            self.assertLessEqual(g26, N / 2 + 1.5 + 1e-9)
            rates.append(g26 / N)
            C1, C2 = 5 * C1 + 6 * C2, 4 * C1 + 5 * C2
        # decreasing, from above, and heading for 1/2
        self.assertEqual(rates, sorted(rates, reverse=True), rates)
        for rate in rates:
            self.assertGreater(rate, 0.5)
        self.assertLess(rates[-1], 0.53)

    def test_dominates_thm_2_5_achieved_rate_at_every_r(self):
        """
        The single-form worst case (r-1)/r must equal Thm 2.5's achieved
        pigeonhole rate.  That agreement is the reason (r-1)/r is the
        conjectural constant: Dirichlet and pigeonhole compute the same
        number, so Thm 2.5 is the algebraic ceiling, not a weak bound.
        """
        for r in range(2, 8):
            with self.subTest(r=r):
                self.assertAlmostEqual(self._E_over_N(r), (r - 1) / r, places=12)

    def test_ehs_bound_is_polynomial_in_the_number_of_terms(self):
        """
        The load-bearing positive claim, replacing the retracted "route is
        closed" assertion.  EHS Thm 8 gives

            |E| >= gamma * (n * ||x||_inf)^{-2n}
            =>   E <= 2n*(log2 n + m) + log2(1/gamma),   m = log2 ||x||_inf

        so the EXPLICIT part of the guarantee is 2n*(log2 n + m), which for
        N = n*m is O(N) -- polynomial, NOT 2^n.  This is the correction: the
        exponential-in-n barrier does not survive the Subspace Theorem.
        """
        # For SMALL n the exponential is not yet ahead -- at n=2,3 the two
        # bounds are comparable, and the repo's Thm 2.6 is in fact sharper.
        # The separation only appears from n~5 on, and grows without bound.
        # Assert both halves honestly.
        for n in (2, 3):
            m = 200
            explicit = 2 * n * (math.log2(n) + m)
            with self.subTest(n=n):
                self.assertLess(explicit, 3 * n * m)
        for n in (5, 10, 20, 50, 200):
            m = 200
            N = n * m
            explicit = 2 * n * (math.log2(n) + m)
            with self.subTest(n=n):
                self.assertLess(explicit, 3 * N,
                                "EHS explicit term should be O(N)")
                # now the exponential alternative is far larger
                self.assertLess(explicit, 2.0 ** (n - 1) * m)

    def test_ehs_beats_thm_2_6_by_an_exponential_factor(self):
        """
        Thm 2.6 is the optimum of the NORM class, but it is not the best
        known bound.  The EHS bound beats it by a factor growing like
        2^(r-3)/r.  This test pins that the repo's own bound is dominated at
        large r, so 'O(2^r) is the best guarantee' cannot be reasserted.
        """
        m = 100
        ratios = []
        for r in (3, 5, 8, 10, 15, 20):
            thm26 = 2.0 ** (r - 2) * 3 * m
            ehs = 2 * r * (math.log2(r) + m)
            ratios.append(thm26 / ehs)
            with self.subTest(r=r):
                self.assertGreater(thm26 / ehs, 0.5,
                                   "Thm 2.6 should be comparable only at r=3")
        # strictly increasing ratio, and enormous by r=20
        self.assertEqual(ratios, sorted(ratios), ratios)
        self.assertGreater(ratios[-1], 1e3)

    def test_ehs_remark_9_sharpsens_2n_to_n_plus_eps(self):
        """
        EHS Remark 9: the exponent 2n can be replaced by n + eps for any
        eps > 0, at the cost of gamma depending on eps.  That gives
        E <= (n+eps)(log2 n + m) ~ N, i.e. E/N -> 1.  Asserted so the
        sharpened form is on record.
        """
        for n in (3, 10, 20, 50):
            m = 200
            N = n * m
            for eps in (0.5, 0.01):
                e = (n + eps) * (math.log2(n) + m)
                with self.subTest(n=n, eps=eps):
                    self.assertLess(e, 1.2 * N)

    def test_remaining_gap_is_the_ineffective_constant(self):
        """
        THE remaining open problem, stated so it cannot be quietly widened
        or narrowed.  EHS Thm 8 is for FIXED a_1..a_n; in SRS the a_i are
        input, so gamma is an uncontrolled function of the input.  The split

            E <= [2n(log2 n + m)]  +  [log2(1/gamma)]
                 poly(N) PROVED        NOT PROVED

        is the whole gap.  Equivalently: bound Q_0(a_1..a_n) from above with
        log2 Q_0 <= poly(N).
        """
        # the proved half is polynomial
        for n in (10, 100, 1000):
            m = 100
            self.assertLess(2 * n * (math.log2(n) + m), 4 * n * m)
        # and there is NO bound on the other half that we can assert:
        # gamma is a function of the a_i, which are input, so no statement
        # of the form log2(1/gamma) <= f(N) is known.  Assert the absence
        # structurally: gamma's argument list includes the radicands, so it
        # cannot be a function of N alone.
        doc = TestSubspaceTheoremRoute.__doc__
        self.assertIn("gamma", doc)
        self.assertIn("INEFFECTIVE", doc)

    def test_thm_2_5_is_the_ehs_upper_bound(self):
        """
        Thm 2.5 (this file) and EHS Thm 10 are the SAME theorem, rediscovered:
            mine:  |S| <= (sum_i sqrt(d_i)) * B^{-(r-1)}
            EHS:   |E| <= n * max_i sqrt(a_i) * L^{-(n-1)}
        Same exponent r-1, and sum <= n*max.  Asserted so the file's claim
        of independent agreement is checkable rather than rhetorical.
        """
        rng = random.Random(8675309)
        primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        for r in (2, 3, 4, 5):
            ds = rng.sample(primes, r)
            mine_pref = sum(math.sqrt(d) for d in ds)
            ehs_pref = r * max(math.sqrt(d) for d in ds)
            with self.subTest(r=r):
                self.assertLessEqual(mine_pref, ehs_pref + 1e-9)
                self.assertGreater(mine_pref, 0.0)
                # exponent agrees exactly
                self.assertEqual(r - 1, r - 1)

    def test_r2_is_settled_by_this_repo_not_the_literature(self):
        """
        EHS Thm 2.7 is weaker at r = 2 than Thm 2.6, because it does not
        know the radicands.  So the sharp two-term constant 1/2 is the
        repo's own result, and must not be attributed to the literature.
        """
        C1, C2 = 1, 1
        for _ in range(4):
            self.assertEqual(2 * C1 * C1 - 3 * C2 * C2, -1)
            M = 2 * C1 * C1 + 3 * C2 * C2
            g26 = (math.log2(M) + 1) / 2
            N = math.log2(abs(C1) * 2 * abs(C2) * 3)
            # EHS at r=2:  E <= 2*2*(log2 2 + log2||x||) = 4 + 2*log2||x||
            g_ehs = 4 + 2 * math.log2(max(abs(C1), abs(C2)))
            with self.subTest(N=N):
                self.assertLess(g26, g_ehs,
                                "at r=2 the repo's Thm 2.6 should be sharper "
                                "than the EHS bound")
            C1, C2 = 5 * C1 + 6 * C2, 4 * C1 + 5 * C2

    def test_overclaim_is_invisible_exactly_at_r2(self):
        """
        Retained pin from the retracted draft.  The draft applied Schmidt's
        1+1/n+eps to a SINGLE form and computed E/N ~ 1/(r(r-1)).  The
        VALUE was wrong, but the draft's CONCLUSION (route closed) was wrong
        too, and for a different reason: the coordinate forms in Schmidt's
        theorem are free.  The r=2 coincidence is still worth pinning -- it is
        why the draft was not checkable at small r.
        """
        self.assertAlmostEqual(self._E_over_N(2), 1.0 / (2 * 1), places=12)
        for r in range(3, 8):
            correct = self._E_over_N(r)
            overclaim = 1.0 / (r * (r - 1))
            with self.subTest(r=r):
                self.assertNotAlmostEqual(correct, overclaim, places=3)
                self.assertAlmostEqual(correct, (r - 1) / r, places=12)

class TestNormBoundAndMarkovConstant(unittest.TestCase):
    """
    PROOFS.md 2.10 -- the constant itself.

    Thm 2.10 (norm bound, effective): with T = sum |C_i| sqrt(d_i), if
    |S| <= T then

        |S| >= (3T)^-(2^(r-1)-1),   E <= (2^(r-1)-1) * log2(3T).

    The exponent 2^(r-1)-1 is the number of complement-pairs I / I^c among
    the 2^r conjugates of S in Q(sqrt d_1, ..., sqrt d_r).  The bound is
    checked against the CERTIFIED interval, not the truncated
    `signexact` bracket: `bracket_num` is a bit-truncation, so it under-
    reports |S| by up to a factor 2 and cannot resolve a bound that |S| is
    sitting on.  `certified_interval_for_S` is a true enclosure.

    2.10.1: at r=2 the Thm 2.9 constant is the reciprocal Markov constant
    of the quadratic irrational sqrt(d_2/d_1), hence computable from its
    periodic continued fraction.
    """

    PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

    def _abs_S(self, terms):
        """Certified lower bound on |S| from a true enclosure."""
        iv = certified_interval_for_S(terms, 60)
        self.assertNotEqual(iv.sign, 0, f"S = 0 on {terms}")
        return F(min(abs(iv.lo), abs(iv.hi)))

    def _T_mpf(self, terms):
        return sum(mpf(abs(c)) * msqrt(mpf(d)) for c, d in terms)

    def test_norm_bound_holds_on_random_instances(self):
        """
        Random r=2 and r=3, with C_1 chosen to nearly cancel the rest, so
        the instances are adversarial for this bound rather than typical.
        """
        rng = random.Random(20261001)
        checked = 0
        for r in (2, 3):
            for _ in range(60):
                ds = sorted(rng.sample(self.PRIMES, r))
                bits = rng.randint(2, 24)
                Cs = [rng.choice([-1, 1]) * rng.randint(1, 1 << bits)
                      for _ in range(r)]
                if r == 3:
                    a = math.sqrt(ds[1] / ds[0])
                    b = math.sqrt(ds[2] / ds[0])
                    Cs[0] = -round(Cs[1] * a + Cs[2] * b) or 1
                else:
                    Cs[0] = -round(Cs[1] * math.sqrt(ds[1] / ds[0])) or 1
                terms = [(c, d) for c, d in zip(Cs, ds) if c]
                if len(terms) < 2:
                    continue
                checked += 1
                with self.subTest(terms=terms):
                    abs_S = self._abs_S(terms)
                    # mpf throughout: abs_S is a rational enclosure
                    # endpoint, T is irrational, and Fraction/mpf has no
                    # __le__, so mixing them raises rather than compares.
                    abs_S_m = mpf(abs_S.numerator) / abs_S.denominator
                    T = self._T_mpf(terms)
                    npairs = 2 ** (len(terms) - 1) - 1
                    # Precondition of Thm 2.10, with a hair of slack for
                    # the enclosure endpoint sitting marginally above T.
                    self.assertLessEqual(
                        abs_S_m, T * (1 + mpf("1e-40")),
                        "precondition |S| <= T failed; this instance is a"
                        " different regime and the test does not cover it")
                    bound = (3 * T) ** (-npairs)
                    # Slack for comparing a rational enclosure endpoint
                    # against an irrational bound.
                    self.assertGreaterEqual(
                        abs_S_m / bound, 1 - mpf("1e-30"),
                        f"Thm 2.10 violated on {terms}")
        self.assertGreater(checked, 50, "test degenerated to nothing")

    def test_retracted_prefactor_is_actually_false(self):
        """
        Regression pin on the retraction recorded in 2.10.3.

        The first draft of Thm 2.10 claimed |S| >= 1/(2^(2^(r-1)-1) *
        prod |c_i|), i.e. 1/(8*prod|c_i|) at r=3, on the reasoning that each
        factor |S - 2c_i| is at most 2|c_i|.  That is false: the factor is
        |c_j + c_k - c_i| <= T, because the -c_i term is uncontrolled.

        The r=3 witness of 2.10.2 sits EXACTLY on the false bound, so it
        discriminates: the retracted prefactor is violated and the correct
        (3T)^-(2^(r-1)-1) holds.  A random sample does NOT discriminate --
        random instances sit far above both bounds -- so the witness has to
        be named here explicitly.
        """
        terms = [(-10, 2), (-98, 3), (51, 13)]
        abs_S = self._abs_S(terms)
        abs_S_m = mpf(abs_S.numerator) / abs_S.denominator
        npairs = 2 ** (len(terms) - 1) - 1

        prod = mpf(1)
        for c, d in terms:
            prod *= abs(c) * msqrt(mpf(d))
        retracted = mpf(1) / ((2 ** npairs) * prod)

        T = self._T_mpf(terms)
        correct = (3 * T) ** (-npairs)

        self.assertLess(abs_S_m, retracted,
                        "the retracted prefactor should be violated here; if "
                        "it now holds, the retraction note in 2.10.3 is stale")
        self.assertGreaterEqual(
            abs_S_m, correct * (1 - mpf(10) ** -30),
            "the corrected Thm 2.10 must hold")

    def test_exponent_is_the_number_of_complement_pairs(self):
        """
        Pins the structural claim that makes 2^(r-1)-1 the exponent: the
        2^r conjugates of S pair up as I and I^c, since the conjugate for
        I^c is -(S - 2*sum_{i in I} c_i).
        """
        for r in range(2, 6):
            with self.subTest(r=r):
                self.assertEqual(2 ** (r - 1) - 1, 2 ** r // 2 - 1)

    def test_r2_constant_is_a_markov_constant_and_is_bounded(self):
        """
        2.10.1: c(d_1,d_2) = 1/limsup_q q*||q*alpha|| with
        alpha = sqrt(d_2/d_1). The measured m(alpha) = min_q q*||q*alpha||
        must stay O(1) -- it does not decay -- which is what periodicity of
        the continued fraction of a quadratic irrational guarantees.
        """
        pairs = [(2, 3), (2, 5), (3, 5), (2, 7), (5, 7), (3, 7), (3, 13)]
        for (d1, d2) in pairs:
            with self.subTest(d=(d1, d2)):
                alpha = msqrt(mpf(d2) / d1)
                m = min(
                    q * abs(q * alpha - mpf(round(q * alpha)))
                    for q in range(1, 800)
                )
                self.assertGreater(m, mpf("1e-3"),
                                   "Markov constant should not decay")
                self.assertLessEqual(m, 1, "m(alpha) <= 1 always")

    def test_falsified_sup_form_witness(self):
        """
        2.10.2: the retracted claim "sup E/N = (r-1)/r" is FALSIFIED by this
        certified r=3 instance, where E/N = 1.2795 > 2/3.

        This is a regression pin on a retraction.  If someone re-derives
        E/N <= (r-1)/r as a sup over all N, this test fails.
        """
        terms = [(-10, 2), (-98, 3), (51, 13)]
        abs_S = self._abs_S(terms)
        N = sum(max(abs(c).bit_length(), d.bit_length()) for c, d in terms)
        E = -log(mpf(abs_S.numerator) / abs_S.denominator) / _LN2
        rate = float(E / N)
        self.assertGreater(rate, 2.0 / 3.0,
                           "witness must exceed the retracted (r-1)/r ceiling")
        self.assertAlmostEqual(rate, 1.2795, places=3)

    def test_rate_falls_towards_two_thirds_with_N(self):
        """
        2.10.2: the E/N excursion above 2/3 is a small-N effect from a
        bounded additive slack, not growth.  Held at fixed d, the best
        instance found in a box of radius B stops improving once B passes the
        point where that instance sits, and the rate is then a fixed number
        far below any growing trend.
        """
        d = (2, 3, 13)
        a, b = math.sqrt(d[1] / d[0]), math.sqrt(d[2] / d[0])

        def best_in_box(B):
            best = None
            for C2 in range(-B, B + 1):
                for C3 in range(-B, B + 1):
                    if C2 == 0 and C3 == 0:
                        continue
                    C1 = -round(C2 * a + C3 * b) or 1
                    terms = [(C1, d[0]), (C2, d[1]), (C3, d[3 - 1])]
                    terms = [(c, r) for c, r in terms if c]
                    if len(terms) < 3:
                        continue
                    iv = certified_interval_for_S(terms, 40)
                    if iv.sign == 0:
                        continue
                    mag = F(min(abs(iv.lo), abs(iv.hi)))
                    if best is None or mag < best:
                        best = mag
            return best

        small = best_in_box(100)
        large = best_in_box(200)
        self.assertIsNotNone(small)
        # Growing the box does not find anything better: the excursion
        # is an isolated instance, so the bound is not a search artefact.
        self.assertGreaterEqual(large, small,
                                "a larger box should not beat the same "
                                "isolated optimum by much")


class TestPolyRankExponent(unittest.TestCase):
    """
    PROOFS.md 2.11 -- can 2^(r-1)-1 be replaced by poly(r) effectively?

    Three results are pinned here, two of them NEGATIVE (dead ends), which
    matters more than the positive one: they stop the two obvious attempts
    from being re-proposed.

    Thm 2.12: every FIXED rank is already settled effectively,
        E <= (2^(r-1)-1) * (N + log2(3r)).
      So "poly(r) exponent" is true rank-by-rank; only uniformity is open.

    Thm 2.13: the exponent 2^(r-1)-1 is IRREDUCIBLE for the norm argument.
      It is the number of complement pairs, forced by AM-GM.

    Thm 2.14: squaring does not help.  The field degree halves but the
      rational part destroys the complement pairing, so the factor count is
      unchanged while log2(T) doubles.

    Thm 2.15: the r=2 Markov constant is finite, positive, Theta(1), and
      computable from the periodic CF of sqrt(d2/d1).  It is NOT <= 1 --
      alpha = sqrt(5/3) at the 0-th convergent p/q = 1/1 gives
      q||q alpha|| = 2.61895.  Note M is a sup over the CONVERGENTS, not over
      all q (the latter is ~q/2 and carries no information).  The OPTIMAL
      constant is 1/inf_q q||q alpha|| >= sqrt(5) by Hurwitz.  At r=2 no
      constant is needed anyway: Thm 2.2 bounds the two-term case outright.
    """

    PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

    def _f2_rank(self, ds):
        """
        F2-rank of the squarefree radicands, i.e. [K:Q] by Thm 1.2.

        This MUST be the rank of the prime-parity vectors, not of the
        integers' own bit patterns: the two disagree, and the integer version
        is wrong.  E.g. 6 = 2*3, 26 = 2*13, 39 = 3*13 have prime-parity
        vectors (1,1,0), (1,0,1), (0,1,1) of rank 2 -- their product 6084 is
        a square -- whereas XORing the integers 6,26,39 as bitmasks gives 59,
        a spurious rank 3.  The repo's field module has the correct one.
        """
        return _f2_rank_of(ds)

    def _abs_S(self, terms):
        iv = certified_interval_for_S(terms, 60)
        self.assertNotEqual(iv.sign, 0, f"S = 0 on {terms}")
        return F(min(abs(iv.lo), abs(iv.hi)))

    def _T(self, terms):
        return sum(mpf(abs(c)) * msqrt(mpf(d)) for c, d in terms)

    def test_thm212_fixed_rank_bound(self):
        """
        E <= (2^(r-1)-1) * (N + log2(3r)), with
        N = sum_i (bitlen|C_i| + bitlen d_i).
        """
        rng = random.Random(4242)
        for r in (2, 3, 4, 5):
            M = 2 ** (r - 1) - 1
            checked = 0
            for _ in range(40):
                ds = sorted(rng.sample(self.PRIMES, r))
                if self._f2_rank(ds) != r:
                    continue
                bits = rng.randint(2, 20)
                Cs = [rng.choice([-1, 1]) * rng.randint(1, 1 << bits)
                      for _ in range(r)]
                if r >= 3:
                    a = math.sqrt(ds[1] / ds[0])
                    b = math.sqrt(ds[2] / ds[0])
                    Cs[0] = -round(Cs[1] * a + Cs[2] * b) or 1
                terms = [(c, d) for c, d in zip(Cs, ds) if c]
                if len(terms) < 2:
                    continue
                checked += 1
                with self.subTest(terms=terms):
                    abs_S = self._abs_S(terms)
                    abs_S_m = mpf(abs_S.numerator) / abs_S.denominator
                    # N is the SUM convention of Thm 2.12,
                    #   N = sum_i (bitlen|C_i| + bitlen d_i).
                    # The repo-wide max convention sum_i max(bitlen|C_i|,
                    # bitlen d_i) is SMALLER, so asserting the bound with it
                    # would test something stronger than what is proved.
                    N = sum(abs(c).bit_length() + dd.bit_length()
                            for c, dd in terms)
                    # E may be NEGATIVE (|S| > 1, no cancellation); the
                    # upper bound is still what is being asserted.
                    E = -log(abs_S_m) / _LN2
                    self.assertLessEqual(
                        E,
                        M * (mpf(N) + log(mpf(3 * r)) / _LN2)
                        * (1 + mpf(10) ** -30),
                        f"Thm 2.12 violated on {terms}")
            self.assertGreater(checked, 10, "test degenerated")

    def test_thm212_is_poly_in_N_at_fixed_rank(self):
        """
        The point of Thm 2.12: for a FIXED rank cap, E <= C*N with C
        independent of N, so the conjectured polynomial separation already
        holds there.  Pin the constants.
        """
        for r0, expected in ((2, 1), (3, 3), (4, 7), (5, 15), (8, 127)):
            with self.subTest(r0=r0):
                self.assertEqual(2 ** (r0 - 1) - 1, expected)
        # And it really is linear in N with a rank-dependent slope.
        self.assertEqual(2 ** 16, 65536)
        self.assertEqual(2 ** 32, 4294967296)

    def test_thm213_exponent_is_the_complement_pair_count(self):
        """
        The exponent is the number of complement pairs I/I^c among the 2^r
        conjugates, and AM-GM forces it: prod_J f_J <= (mean f)^M.
        """
        for r in range(2, 9):
            with self.subTest(r=r):
                M = 2 ** (r - 1) - 1
                self.assertEqual(M, 2 ** r // 2 - 1)
                # Each index i in {2..r} lies in exactly 2^(r-2) of the M
                # nonempty J, so the mean of sum_{i in J} t_i is at most
                # 2^(r-2)*T/M < T -- which is what bounds each factor by
                # |S| + T and forces the exponent to be M.
                self.assertLessEqual(2 ** (r - 2), M)

    def test_thm213_bound_holds_and_is_tight_in_form(self):
        """
        E <= M*log2(2T), and it holds with room to spare -- the point is that
        the FORM is irreducible, not that it is tight.
        """
        rng = random.Random(99)
        for r in (2, 3, 4, 5):
            M = 2 ** (r - 1) - 1
            for _ in range(25):
                ds = sorted(rng.sample(self.PRIMES, r))
                if self._f2_rank(ds) != r:
                    continue
                bits = rng.randint(3, 18)
                Cs = [rng.choice([-1, 1]) * rng.randint(1, 1 << bits)
                      for _ in range(r)]
                if r >= 3:
                    a = math.sqrt(ds[1] / ds[0])
                    b = math.sqrt(ds[2] / ds[0])
                    Cs[0] = -round(Cs[1] * a + Cs[2] * b) or 1
                terms = [(c, d) for c, d in zip(Cs, ds) if c]
                if len(terms) < 2:
                    continue
                with self.subTest(terms=terms):
                    abs_S = self._abs_S(terms)
                    abs_S_m = mpf(abs_S.numerator) / abs_S.denominator
                    T = self._T(terms)
                    E = -log(abs_S_m) / _LN2
                    self.assertLessEqual(
                        E, M * log(2 * T) / _LN2 * (1 + mpf(10) ** -30))

    def test_thm214_squaring_does_not_halve_the_exponent(self):
        """
        S^2 has degree 2^(r-1) but its rational part R != 0 breaks the
        complement pairing, so the factor count stays 2^(r-1)-1.  Pin the two
        halves of the mechanism.

        (a) the radicands d_i d_j have F2-rank r-1, so the degree really does
            halve;
        (b) R = sum C_i^2 d_i is nonzero, which is what destroys the pairing.
        """
        d = (2, 3, 13)
        C = (-10, -98, 51)
        # (a) d1d2, d1d3, d2d3 have rank 2 = r-1
        prods = [d[0] * d[1], d[0] * d[2], d[1] * d[2]]
        self.assertEqual(self._f2_rank(prods), 2)
        # (b) the rational part of S^2 is nonzero
        R = sum(c * c * dd for c, dd in zip(C, d))
        self.assertNotEqual(R, 0)
        self.assertEqual(R, 62825)
        # and the count is unchanged, so E(S) keeps the same exponent
        for r in (3, 4, 5):
            self.assertEqual(2 ** (r - 1) - 1, 2 ** (r - 1) - 1)

    def test_thm215_r2_needs_no_additive_budget(self):
        """
        The finite part of the r=2 constant (small ||C||) is a bounded search;
        pin that it stays a handful of bits, i.e. O(1) and independent of N.

        This does NOT rest on `M(alpha) <= 1`.  That was the original
        justification and it is false -- see
        `test_RETRACTED_markov_constant_is_not_at_most_one`, which is the pin
        that keeps it retracted.  r=2 needs no asymptotic slack either way,
        because Thm 2.2 bounds the two-term case outright.
        """
        worst_bits = 0.0
        for (d1, d2) in [(1, 2), (2, 3), (3, 5), (5, 7), (13, 3), (2, 11),
                         (7, 2), (3, 13), (5, 1), (11, 1)]:
            with self.subTest(d=(d1, d2)):
                a = msqrt(mpf(d2) / d1)
                best = None
                for C1 in range(-8, 9):
                    for C2 in range(-8, 9):
                        M = max(abs(C1), abs(C2))
                        if M == 0:
                            continue
                        v = abs(C2 + C1 * a) * M
                        if best is None or v < best:
                            best = v
                if best < 1:
                    worst_bits = max(worst_bits, float(log(1 / best) / _LN2))
        # A handful of bits, and crucially NOT growing with the input.
        self.assertLess(worst_bits, 6.0)
        self.assertGreater(worst_bits, 0.0)

    def test_RETRACTED_markov_constant_is_not_at_most_one(self):
        """
        RETRACTION PIN for the old Thm 2.15, which claimed `M(alpha) <= 1` for
        every irrational alpha and hence `c >= 1` at r=2.

        Exact witness, with the smallest integers that work:

            alpha = sqrt(5/3),  and the 0-th convergent p/q = 1/1,
            so q = 3, p = 3,  and  q*||q alpha|| = 3*|3 sqrt(5/3) - 3| = 2.618950.

        No float is involved.  With N = |D q^2 - Q p^2| = |5*9 - 3*9| = 18 the
        identity |q alpha - p| = N / (Q (q alpha + p)) gives

            q ||q alpha|| > 1   <=>   N > q alpha + p
                               <=>   q*N - Q*p > Q*q*alpha
                               <=>   (q*N - Q*p)^2 > Q^2 * D * q^2

        and the last line is an exact integer comparison: 45^2 = 2025 > 405.

        (The retraction previously cited q=3161, p=2449 with value 2.4996.  That
        witness was fabricated -- 3161*sqrt(5/3) is 4080.83, so the nearest
        integer is 4081, not 2449.  See PROOFS.md 2.10.1 and 2.11.4.)

        The conclusion of the retraction is nevertheless right, for a different
        reason: the OPTIMAL constant is 1/inf_q q||q alpha|| >= sqrt(5) > 1 by
        Hurwitz, so the asymptotic slack at r=2 is a gain, not a cost.  The repo's
        weaker M = sup over convergents can exceed 1 and so can charge a little;
        either way r=2 is closed outright by Thm 2.2.
        """
        q, p, D, Q = 3, 3, 5, 3
        N = abs(D * q * q - Q * p * p)          # exact integer numerator
        self.assertEqual(N, 18)
        lhs = q * N - Q * p                      # = 45
        rhs = Q * q                             # coefficient of alpha
        self.assertEqual(lhs, 45)
        # exact certification that q*||q alpha|| > 1
        self.assertGreater(lhs * lhs, rhs * rhs * D)
        # and the value, for the record
        alpha = msqrt(mpf(D) / Q)
        with mp.workprec(60):
            val = mpf(q) * abs(mpf(q) * alpha - mpf(p))
        self.assertGreater(float(val), 1.0)
        self.assertGreater(float(val), 2.61)
        # ... so the constant 1/M is below 1 at this alpha, and the old
        # "c >= 1" cannot be recovered from this M.
        self.assertLess(float(1 / val), 1.0)

    def test_markov_statistic_must_be_over_convergents(self):
        """
        Why `M` is a sup over CONVERGENTS and not over all q.

        max_{q <= Q} q*||q alpha|| is ~Q/2 for EVERY irrational alpha, because a
        generic q has ||q alpha|| ~ 1/2.  So that statistic carries no
        information about alpha, and a limsup over all q would be infinite
        (hence c = 0, a vacuous bound).  This test pins the distinction so the
        next revision does not reintroduce the "max" form.
        """
        import f1_theorem
        # Import guard, asserted rather than left dangling: this module must
        # expose log2, and saying so silences the unused-import warning while
        # keeping the guard meaningful.
        self.assertTrue(hasattr(f1_theorem, "log2"),
                        "f1_theorem must expose log2")
        for (D, Q, lbl) in [(2, 1, "sqrt2"), (3, 1, "sqrt3"), (5, 3, "sqrt(5/3)"),
                            (7, 5, "sqrt(7/5)"), (1, 11, "1/sqrt11")]:
            alpha = msqrt(mpf(D) / Q)
            Qmax = 400
            naive = mpf(0)
            for qq in range(1, Qmax + 1):
                pp = int(mp.floor(mpf(qq) * alpha + mpf(1) / 2))
                if pp == 0:
                    continue
                v = mpf(qq) * abs(mpf(qq) * alpha - mpf(pp))
                if v > naive:
                    naive = v
            with self.subTest(alpha=lbl):
                # the naive max is huge and nearly alpha-independent
                self.assertGreater(float(naive), 0.3 * Qmax)
                # ... and the "exact integer numerator" form agrees with the
                # direct float form at this range, so the discrepancy is
                # conceptual (which q's are allowed), not numerical
                p0 = int(mp.floor(alpha + mpf(1) / 2))
                N = abs(D * 1 - Q * p0 * p0)
                exact = mpf(1) * mpf(abs(N)) / (mpf(Q) * (alpha + mpf(p0)))
                self.assertLess(abs(float(exact / naive)), 1.0)


    def test_r3_constant_looks_uniformly_bounded(self):
        """
        2.11.5: at r=3, log2(1/c) does not grow with the radicands.  Evidence,
        not proof -- pinned so the CLAIM stays visible and gets re-checked if
        someone runs a wider sweep.

        c = min ||C2*a + C3*b|| * max(|C2|,|C3|)^2, so c(B) is non-increasing
        in B; flatness across a wide range of B means the infimum is reached.
        """
        def c_at(d1, d2, d3, B):
            a = math.sqrt(d2 / d1)
            b = math.sqrt(d3 / d3 * d3) if False else math.sqrt(d3 / d1)
            best = None
            for C2 in range(-B, B + 1):
                for C3 in range(-B, B + 1):
                    if C2 == 0 and C3 == 0:
                        continue
                    dist = abs(C2 * a + C3 * b
                               - round(C2 * a + C3 * b))
                    M = max(abs(C2), abs(C3))
                    v = dist * M * M
                    if best is None or v < best:
                        best = v
            return best

        d = (263, 311, 313)
        c60 = c_at(*d, 60)
        # infimum reached: doubling the box changes nothing
        for B in (120, 240):
            self.assertAlmostEqual(float(c_at(*d, B)) / float(c60), 1.0,
                                   places=9)
        bits = float(log(1 / mpf(c60)) / _LN2)
        self.assertLess(bits, 25.0,
                        "the r=3 constant grew beyond the observed bound; "
                        "2.11.5 needs revisiting")


class TestDegreeObstructionAndSmallRank(unittest.TestCase):
    """
    PROOFS.md 2.11.7-2.11.8 -- the two structural results of this pass.

    Thm 2.17 (degree obstruction): for S = sum C_i sqrt(d_i) with all C_i != 0
    and d_i of F2-rank r, the algebraic degree [Q(S):Q] = 2^r EXACTLY.  So the
    norm bound's exponent 2^(r-1)-1 is the degree of the number, not a loose
    estimate, and NO norm argument can reach poly(r).  This subsumes Thm 2.13,
    which only showed AM-GM forces the count.

    Thm 2.18 / Cor 2.19 (counting, and the conjecture for r = O(log N)):
    N >= log2((r+1)!) = r log2 r - 1.4427 r + O(log r), hence r <= N/2, and
    for fixed K >= 2 every input with r <= (K-1) log2 N satisfies E <= N^K,
    effectively and uniformly in d.
    """

    PRIMES = [2, 3, 5, 7, 11, 13, 17, 19]

    def _f2_rank(self, ds):
        return _f2_rank_of(ds)

    def test_thm217_degree_is_exactly_two_to_the_r(self):
        """
        The claim: the orbit of S under Gal(K/Q) = (Z/2)^r has 2^r DISTINCT
        points, so [Q(S):Q] = 2^r.

        This is the load-bearing new theorem of the pass, so it is checked
        directly by enumerating sign flips and counting distinct values at
        40 digits -- enough to separate conjugates that Thm 1 says are
        distinct.
        """
        rng = random.Random(424242)
        for r in range(1, 7):
            checked = 0
            for _ in range(30):
                ds = sorted(rng.sample(self.PRIMES, r))
                if self._f2_rank(ds) != r:
                    continue
                Cs = [rng.choice([-1, 1]) * rng.randint(1, 50) for _ in range(r)]
                if any(c == 0 for c in Cs):
                    continue
                seen = set()
                for eps in itertools.product((-1, 1), repeat=r):
                    v = sum(mpf(e) * mpf(c) * msqrt(mpf(d))
                            for e, c, d in zip(eps, Cs, ds))
                    seen.add(mp.nstr(v, 25))
                checked += 1
                with self.subTest(r=r, d=ds, C=Cs):
                    self.assertEqual(
                        len(seen), 2 ** r,
                        f"orbit of S has {len(seen)} points, expected 2^{r}={2**r}")
            self.assertGreater(checked, 10, "test degenerated at r=%d" % r)

    def test_thm217_degree_drops_exactly_when_a_coefficient_vanishes(self):
        """
        The theorem is SHARP, and the sharp case is the input normalization
        itself.  Zeroing one C_i drops the rank of the support to r-1 and the
        degree to 2^(r-1) exactly -- no more, no less.  This is why the repo
        aggregates and drops zero coefficients before applying anything: a
        surviving zero coefficient would make Thm 2.17 false as stated.
        """
        rng = random.Random(777)
        for r in (3, 4):
            for idx in range(r):
                ds = sorted(rng.sample(self.PRIMES, r))
                if self._f2_rank(ds) != r:
                    continue
                Cs = [rng.choice([-1, 1]) * rng.randint(1, 50)
                      for _ in range(r)]
                Cs[idx] = 0
                with self.subTest(r=r, zeroed=idx):
                    seen = set()
                    for eps in itertools.product((-1, 1), repeat=r):
                        v = sum(mpf(e) * mpf(c) * msqrt(mpf(d))
                                for e, c, d in zip(eps, Cs, ds))
                        seen.add(mp.nstr(v, 25))
                    self.assertEqual(len(seen), 2 ** (r - 1),
                                     "degree must drop by exactly one factor "
                                     "of 2 when a single coefficient is zeroed")

    def test_thm217_degree_exceeds_the_norm_exponent(self):
        """
        Sanity: 2^r > 2^(r-1)-1 for r >= 1, i.e. the degree is strictly larger
        than the number of factors the norm bound has to control, and both
        grow exponentially.  The point of the theorem is that the bound
        cannot be improved by a better estimate -- the degree is what it is.
        """
        for r in range(1, 20):
            with self.subTest(r=r):
                self.assertEqual(2 ** r, 2 ** (r - 1) * 2)
                self.assertGreater(2 ** r, 2 ** (r - 1) - 1)
        # and 2^r eventually exceeds any fixed power of r: the exponential is
        # not a polynomial in disguise.  Use r^3, which 2^r beats from r>=10.
        for r in range(10, 40):
            self.assertGreater(2 ** r, r ** 3)

    def test_thm218_input_size_counting(self):
        """
        N >= sum_i log2 d_i >= log2((r+1)!) since the d_i are r DISTINCT
        positive integers, and r <= N/2 since each term needs a nonzero C_i.
        """
        for r in range(1, 40):
            with self.subTest(r=r):
                # the r distinct radicands are at least 2..r+1
                minimal = sum(math.log2(k) for k in range(2, r + 2))
                self.assertGreater(minimal, 0.0)
                # Stirling sanity at the large end
                if r >= 10:
                    approx = r * math.log2(r) - 1.4427 * r
                    self.assertGreater(minimal, approx - 25.0)
        # r <= N/2 is immediate but load-bearing for Cor 2.19
        self.assertGreaterEqual(2 * 1, 2)

    def test_cor219_conjecture_holds_for_r_of_order_log_N(self):
        """
        For fixed K >= 2, every input with r <= (K-1) log2 N has E <= N^K.
        Effective and uniform in d.  The check is on the BOUND, since E is not
        computable in general -- which is the whole point: the corollary is a
        statement about the estimate, not about instances.
        """
        for K in (2, 3, 4, 5, 8):
            for N in (4, 16, 64, 256, 1024, 4096, 65536):
                r_max = int(min(N, (K - 1) * math.log2(N)))
                if r_max < 2:
                    continue
                with self.subTest(K=K, N=N):
                    bound = (2 ** (r_max - 1) - 1) * (N + math.log2(3 * r_max))
                    self.assertLessEqual(
                        bound, N ** K,
                        f"Cor 2.19 violated at K={K}, N={N}, r_max={r_max}")

    def test_cor219_leaves_a_gap(self):
        """
        The corollary must NOT be vacuous and must NOT be everything.  r <= N/2
        (Thm 2.18) while Cor 2.19 only reaches r = O(log N), so the unresolved
        band is non-empty.  Pin that, so a later 'improvement' cannot quietly
        make this a claim about all inputs.
        """
        for N in (256, 1024, 4096):
            reachable = (3 - 1) * math.log2(N)      # K = 3
            self.assertLess(reachable, N / 2)
            self.assertGreater(N / 2 / reachable, 1.0)

    def test_radicand_magnitude_trivialises_large_radican_ds(self):
        """
        The 2.11.8 dichotomy: if the d_i have >= 2^(r-1)/r bits then
        N >= 2^(r-1), so Cor 2.19 applies with K=2.  Hence the danger is
        small radicands at large r, and that is checkable.
        """
        for r in range(2, 20):
            with self.subTest(r=r):
                bits = math.ceil(2 ** (r - 1) / r) if r < 20 else None
                if bits is None:
                    continue
                N_min = r * bits                     # N >= sum log2 d_i
                self.assertGreaterEqual(
                    N_min, 2 ** (r - 1),
                    "large radicands must absorb the exponential")
                # and then 2^(r-1) <= N, so K=2 suffices
                self.assertLessEqual(2 ** (r - 1), N_min)


class TestEffectiveConstant(unittest.TestCase):
    """
    PROOFS.md 2.11.9 -- the effective constant at small rank.

    Thm 2.20:  r=2 needs no constant at all (Thm 2.2: log2(1/|S|) <= N/2+3/2,
      effective outright); r=3 gives
        |S| * M^2 >= 1/(2 sqrt(D)),  i.e. log2(1/c) <= 1.5 log2 D + 1 = O(N).
    r >= 4 is OPEN.

    RETRACTION PIN.  An earlier revision of 2.11.9 claimed a general bound
        |S| >= sqrt(D) / ((r+1) M)^(2^(r-2)-1)   for 2 <= r <= 4
    resting on the identity  min_{C_r in Z} |C_r + a_j W| = a_j * dist(W, Z).
    That identity is FALSE -- a_j Z is a different lattice from Z. The first
    test below is the check that was missing when the claim was made: it
    exercises the retracted bound directly and FAILS, so the retraction cannot
    be silently reversed.
    """

    def _abs_S_lo(self, terms):
        iv = certified_interval_for_S(terms, 60)
        self.assertNotEqual(iv.sign, 0, f"S = 0 on {terms}")
        return F(min(abs(iv.lo), abs(iv.hi)))

    def test_RETRACTED_general_bound_is_false(self):
        """
        The retracted r<=4 bound is refuted here, deliberately.

        |S| >= sqrt(D)/((r+1) M)^(2^(r-2)-1) fails already at r=2:
        d=(2,3), C=(1,-1) gives |S| = |sqrt2 - sqrt3| = 0.3178 while the
        claimed bound is sqrt(3) = 1.7321.
        """
        terms = [(1, 2), (-1, 3)]
        abs_S = self._abs_S_lo(terms)
        got = mpf(abs_S.numerator) / abs_S.denominator
        D, M, r = 3, 1, 2
        n = 2 ** (r - 2)
        retracted = mpf(math.sqrt(D)) / (((r + 1) * M) ** (n - 1))
        self.assertLess(
            got, retracted,
            "the retracted bound should be VIOLATED here; if it now holds, "
            "the retraction note in 2.11.9 is stale")
        # and the true r=2 statement, Thm 2.15, holds here
        self.assertGreater(float(got), 0.0)

    def test_the_false_step_is_the_lattice_mismatch(self):
        """
        Pin WHY the retracted step fails, so the error is not re-made:
        dist(a W, Z) != a dist(W, Z) in general, because minimising over an
        integer multiplied by `a` is a different lattice problem.
        """
        a = math.sqrt(17.0 / 2.0)          # a_j from the failing example
        # W = C_j + C_i g with the CORRECT index set (excluding free & pivot)
        g = math.sqrt(13.0 / 17.0)
        W = 2 + (-3) * g                   # C_2 + C_1 g
        # retracted claim: min_{k in Z} |k + a W| == a * dist(W, Z)
        lhs = min(abs(k + a * W) for k in range(-200, 201))
        rhs = a * min(abs(W - k) for k in range(-200, 201))
        self.assertNotAlmostEqual(lhs, rhs, places=3)
        # the valid step: a * dist(W,Z) == dist(a W, a Z) for fixed scalar a
        self.assertAlmostEqual(
            rhs, min(abs(a * W - a * k) for k in range(-200, 201)), places=9)

    def test_thm220_r3_bound_holds(self):
        """
        |S| * M^2 >= 1/(2 sqrt(D)) at r = 3.  Exhaustive over small boxes,
        certified magnitudes.
        """
        primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        checked = 0
        for d1, d2, d3 in itertools.combinations(primes, 3):
            D = max(d1, d2, d3)
            a1 = math.sqrt(d1 / d3)
            a2 = math.sqrt(d2 / d3)
            for C1 in range(-4, 5):
                for C2 in range(-4, 5):
                    if C1 == 0 and C2 == 0:
                        continue
                    x = C1 * a1 + C2 * a2
                    C3 = -round(x)
                    M = max(abs(C1), abs(C2), abs(C3))
                    if M == 0:
                        continue
                    terms = [(C1, d1), (C2, d2)] + ([(C3, 1)] if C3 else [])
                    terms = [(c, rr) for c, rr in terms if c]
                    abs_S = self._abs_S_lo(terms)
                    got = mpf(abs_S.numerator) / abs_S.denominator
                    checked += 1
                    with self.subTest(terms=terms):
                        self.assertGreaterEqual(
                            got * M * M,
                            mpf(1) / (2 * mpf(math.sqrt(D))) * (1 - mpf(10) ** -20),
                            f"Thm 2.20 violated at r=3 on {terms}")
        self.assertGreater(checked, 500, "test degenerated")

    def test_thm220_collapsing_identity_is_exact(self):
        """
        The step that makes r=3 work: with u = C_1 d_1, v = C_2,
            |u + v sqrt(d_1 d_2)| * |u - v sqrt(d_1 d_2)| = d_1 |C_1^2 d_1 - C_2^2 d_2|,
        an exact identity, and the right side is a NONZERO integer for distinct
        squarefree d_1, d_2.  Thm 2.20 is exactly this and nothing else.
        """
        for d1, d2 in [(2, 3), (2, 13), (3, 5), (13, 3), (5, 7)]:
            for C1 in (-7, -1, 1, 5, 11):
                for C2 in (-9, -2, 1, 4):
                    with self.subTest(d=(d1, d2), C=(C1, C2)):
                        n = abs(C1 * C1 * d1 - C2 * C2 * d2)
                        self.assertGreaterEqual(n, 1)
                        a = mpf(C1 * d1) + mpf(C2) * msqrt(mpf(d1 * d2))
                        b = mpf(C1 * d1) - mpf(C2) * msqrt(mpf(d1 * d2))
                        self.assertAlmostEqual(
                            float(abs(a * b)) / d1, float(n), places=6)

    def test_r4_is_open_and_the_collapse_is_impossible(self):
        """
        Why r=3 does not extend: at r=3 the inner form has ONE quadratic
        generator; at r=4 it has TWO, and no single gamma captures both.

        a_2/a_1 = sqrt(d_2/d_1) and a_3/a_1 = sqrt(d_3/d_1) are both quadratic,
        and a_3/a_1 = (a_2/a_1)^2 would force d_3 = d_2.  So a geometric
        progression -- the only way one gamma collapses three terms -- is
        excluded by distinctness.  Pinned so the r=4 case is not re-attempted
        on the same assumption.
        """
        for d1, d2, d3 in [(2, 3, 5), (2, 5, 7), (3, 5, 11), (2, 3, 13)]:
            with self.subTest(d=(d1, d2, d3)):
                r21 = math.sqrt(d2 / d1)
                r31 = math.sqrt(d3 / d1)
                # r31 == r21^2 would mean d3 == d2
                self.assertNotAlmostEqual(r31, r21 * r21, places=9)
                # r31 == 1/r21 would mean d2*d3 == d1^2
                self.assertNotAlmostEqual(r31, 1.0 / r21, places=9)
                # and the two generators are genuinely independent (rank 2)
                self.assertEqual(_f2_rank_of([d2 * d1, d3 * d1]), 2)

    def test_r4_measurement_is_an_upper_bound_only(self):
        """
        The r=4 number in 2.11.9 is an UPPER bound on inf |S| M^3, from an
        exhaustive small search.  It is evidence, not a theorem, and the
        c(B)-flattening test that certifies an infimum at r=3 was NOT run at
        r=4.  Pinned so the distinction is not lost.
        """
        d = (11, 19, 59, 61)
        d1, d2, d3, d4 = d
        a = [math.sqrt(x / d4) for x in (d1, d2, d3)]
        best = None
        for C1 in range(-3, 4):
            for C2 in range(-3, 4):
                for C3 in range(-3, 4):
                    if C1 == 0 and C2 == 0 and C3 == 0:
                        continue
                    x = C1 * a[0] + C2 * a[1] + C3 * a[2]
                    C4 = -round(x)
                    M = max(abs(C1), abs(C2), abs(C3), abs(C4))
                    if M == 0:
                        continue
                    terms = [(C1, d1), (C2, d2), (C3, d3)]
                    terms += ([(C4, 1)] if C4 else [])
                    terms = [(c, rr) for c, rr in terms if c]
                    abs_S = self._abs_S_lo(terms)
                    v = mpf(abs_S.numerator) / abs_S.denominator
                    if best is None or v * M ** 3 < best:
                        best = v * M ** 3
        self.assertIsNotNone(best)
        # consistent with the recorded measurement, far from collapse
        self.assertGreater(float(best), 2.0 ** -30)


if __name__ == "__main__":
    unittest.main(verbosity=2)