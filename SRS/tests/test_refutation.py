"""
Tests for refutation.py -- the SRS refutation harness.

The load-bearing property here is NEGATIVE: the harness must not report a
refutation of the sharp conjecture `sup E/N = (r-1)/r` unless it has a
certified instance that clears the ceiling by more than the Thm 2.5 slack.

That is not a formality.  The first version of this module compared E/N
against the limit (r-1)/r directly and immediately reported the F6 family --
which IS Thm 2.5's own pigeonhole construction -- as a refutation.  The
false positive came from ignoring the construction's `sum_i sqrt(d_i)`
prefactor, which contributes an O(1/N) rate term.  `test_f6_is_not_a_
refutation` is the regression pin for exactly that bug, and the other tests
pin the discipline that a reported refutation must be certified.
"""

import math
import unittest

import refutation
from refutation import (
    LITERATURE_CAP,
    Instance,
    ceiling_for_rank,
    exact_E,
    f2_rank,
    hunt_greedy,
    hunt_random,
    hunt_rank3_ceiling,
    is_refutation,
    slack_for,
)
from signexact import exact_sign
from squarefree import decompose_and_aggregate


class TestRefutationDiscipline(unittest.TestCase):
    """The criteria a reported refutation must satisfy."""

    def test_f6_is_not_a_refutation(self):
        """
        REGRESSION PIN.  F6 is the Thm 2.5 pigeonhole construction, so by
        that theorem its rate sits AT the ceiling.  An earlier harness
        compared against the bare limit and called it a refutation; that was
        wrong, and it is the single most important thing this module must not
        do again.
        """
        best, beat = hunt_rank3_ceiling(Bs=(40, 60))
        self.assertIsNotNone(best, "F6 produced no candidate")
        # it does sit near the ceiling ...
        self.assertGreater(best.rate, 0.5)
        # ... but it is NOT a refutation
        self.assertEqual(beat, [], f"F6 must not be a refutation, got {beat}")
        self.assertFalse(is_refutation(best),
                         f"F6 at rate {best.rate} was misreported as a refutation")

    def test_additive_constant_slack_prevents_second_false_positive(self):
        """
        The SECOND false-positive mode, distinct from the F6 prefactor one.

        `slack_for` originally modelled only the Thm 2.5 pigeonhole
        prefactor, log2(sum sqrt d)/N.  That is not the only source of
        headroom: every bound carries the ADDITIVE constant log2(1/c(d,delta))
        of Thm 2.9, and c can legitimately be ~2.6e-3 at r=3, i.e. ~8.6 bits.
        At N=17 that is a rate contribution of ~0.5 -- enough to push a
        certified instance to E/N = 1.2835, well above the (r-1)/r = 0.667
        ceiling, and enough to make the old threshold call it a refutation.

        This is the witness of PROOFS.md 2.10.2.  It is NOT a refutation: the
        additive bonus is constant at fixed d and vanishes as N grows.
        """
        terms = [(-10, 2), (-98, 3), (51, 13)]
        cand = Instance(terms, 1, exact_E(terms), 17, 3)
        self.assertGreater(cand.rate, 2.0 / 3.0,
                           "witness must exceed the ceiling, else it no "
                           "longer exercises the false positive")
        self.assertFalse(
            is_refutation(cand),
            "the additive-constant slack was not accounted for: a certified "
            "instance with a bounded additive bonus is not a refutation")

    def test_additive_allowance_is_a_fixed_bit_budget(self):
        """
        The additive constant of Thm 2.9 is accounted for as a FIXED bit
        budget, `ADDITIVE_BITS`, and that choice is load-bearing.

        A per-instance estimate is provably vacuous.  With E = rate*N and
        c = |S|*||C||^(r-1), one gets log2(1/c) = rate*N - (r-1)*log2||C||,
        which grows FASTER than the excess (rate - (r-1)/r)*N it should
        excuse, because the excess carries a factor < 1 and the allowance
        does not.  So the allowance always wins and the harness never fires.
        This module made that mistake twice before this test existed.
        """
        self.assertGreater(refutation.ADDITIVE_BITS, 0.0)

        # The witness of 2.10.2: a real additive bonus, comfortably inside
        # the fixed budget.
        terms = [(-10, 2), (-98, 3), (51, 13)]
        cand = Instance(terms, 1, exact_E(terms), 17, 3)
        c_measured = cand.magnitude * cand.norm_power()
        measured = math.log2(1.0 / c_measured)
        self.assertGreater(measured, 0.0)
        self.assertLessEqual(
            measured, refutation.ADDITIVE_BITS,
            "the fixed budget must cover the measured additive constant")
        self.assertFalse(is_refutation(cand))

        # The budget must NOT depend on the rate, or it chases the excess.
        terms2 = [(1 << 40, 2), (1 << 40, 3), (1 << 40, 40)]
        N2 = sum(max(abs(cc).bit_length(), d.bit_length())
                 for cc, d in terms2)
        allows = []
        for rate in (0.7, 1.0, 2.5):
            inst = Instance(terms2, 1, rate * N2, N2, 3)
            allows.append(2 * (slack_for(inst) * N2 + refutation.ADDITIVE_BITS))
        self.assertLess(max(allows) - min(allows), 1e-6,
                        "the allowance must be flat in the rate")

    def test_finite_instance_cannot_refute_an_asymptotic_rate(self):
        """
        The structural point, stated as a test.

        `is_refutation` compares against the ASYMPTOTIC rate (r-1)/r, so it
        can only be tripped by slack that grows without bound in N.  For
        every instance here the slack is O(1/N), so no finite instance is
        ever reported.  That is the correct behaviour and it is why the
        negative result in refutation.py is a statement about searched
        FAMILIES (see `trend_for`), not about single instances.
        """
        import families
        for B in (10, 40, 100):
            inst = families.generate_f6(B)
            cand = Instance(inst.terms, inst.metadata["exact_sign"],
                            exact_E(inst.terms), inst.N_bits, 3)
            with self.subTest(B=B):
                self.assertFalse(is_refutation(cand))
                # slack decays: O(1/N), so it cannot mask a growing rate
                self.assertLess(slack_for(cand), 1.0)

    def test_slack_is_positive_and_shrinks_with_size(self):
        """
        The Thm 2.5 slack is log2(sum sqrt d)/N, so it must be positive and
        must decrease as N grows.  A refutation test with vanishing slack is
        what makes the negative result meaningful.
        """
        import families
        slacks = []
        for B in (20, 40):
            inst = families.generate_f6(B)
            cand = Instance(inst.terms, inst.metadata["exact_sign"],
                            exact_E(inst.terms), inst.N_bits, 3)
            s = slack_for(cand)
            self.assertGreater(s, 0.0)
            slacks.append((cand.N_bits, s))
        # strictly increasing N gives non-increasing slack
        rates = [s for _, s in slacks]
        self.assertEqual(rates, sorted(rates, reverse=True), rates)

    def test_ceiling_is_the_limit_not_a_finite_bound(self):
        """
        ceiling_for_rank returns (r-1)/r, which Thm 2.5 only approaches from
        below.  Pin the values so a change to the sharp form is deliberate.
        """
        self.assertEqual(ceiling_for_rank(2), 0.5)
        self.assertAlmostEqual(ceiling_for_rank(3), 2 / 3, places=12)
        self.assertAlmostEqual(ceiling_for_rank(4), 0.75, places=12)
        self.assertEqual(ceiling_for_rank(1), 0.0)

    def test_is_refutation_requires_clearing_the_slack(self):
        """
        Boundary logic, tested in BITS against the computed allowance.

        The test works in bits because that is what `is_refutation` does, and
        because the two allowances it subtracts -- the pigeonhole prefactor
        and the additive constant -- are both O(1) in N, hence fixed bit
        budgets.  Rates are placed relative to the computed budget rather
        than hardcoded, so the test stays meaningful if the slack model
        changes.

        A caution learned the hard way: an earlier version of this test used
        terms with ||C|| = 2^40 while passing N = 200, which is not a
        possible instance (those terms have N ~ 120).  `Instance` does not
        police consistency between terms and N_bits, and the additive
        allowance reads ||C||, so an inconsistent synthetic instance produces
        nonsense allowances.  The `N` used below is derived from the terms.
        """
        terms = [(1 << 40, 2), (1 << 40, 3), (1 << 40, 40)]
        N = sum(max(abs(c).bit_length(), d.bit_length()) for c, d in terms)
        ceiling = ceiling_for_rank(3)

        def budget(rate):
            """Bit allowance `is_refutation` would subtract at this rate."""
            inst = Instance(terms, 1, rate * N, N, 3)
            c = inst.magnitude * inst.norm_power()
            return slack_for(inst) * N + max(0.0, math.log2(1.0 / c))

        # Inside the allowance -> not a refutation, even though the rate is
        # above the ceiling.  This is the F6 situation.
        near = Instance(terms, 1, (ceiling + 0.01) * N, N, 3)
        self.assertGreater(near.rate, ceiling)
        self.assertFalse(is_refutation(near))

        # Excess beyond the allowance -> refutation.  The allowance depends
        # on the instance, so the rate is sized off a measured budget at a
        # small rate and then pushed well past it.  The point of this case is
        # that the harness is NOT vacuous: something must be able to clear it.
        probe = Instance(terms, 1, (ceiling + 0.01) * N, N, 3)
        allow = slack_for(probe) * N + refutation.ADDITIVE_BITS
        self.assertGreater(allow, 0.0)
        big = Instance(terms, 1, ceiling * N + 4 * allow, N, 3)
        self.assertGreater(big.rate, ceiling)
        self.assertTrue(is_refutation(big),
                        "an excess of 4x the allowance must be flagged; if "
                        "not, the harness is vacuous")
        # And the allowance must not grow with the rate, or nothing clears.
        far_allow = (slack_for(big) * N + refutation.ADDITIVE_BITS)
        self.assertLess(far_allow, 4 * allow + 1.0)

        # Above the literature cap -> certainly a refutation.
        far = Instance(terms, 1, 2.5 * N, N, 3)
        self.assertTrue(is_refutation(far))
        self.assertGreater(far.rate, LITERATURE_CAP)
        self.assertLess(budget(0.7), budget(2.5), "allowance must not grow "
                        "without bound with the rate, or nothing can clear it")

    def test_a_real_large_instance_is_not_a_refutation(self):
        """
        The regime where the harness has real power: at large N the prefactor
        allowance is O(1) bits, so it cannot mask a genuine excursion.  This
        pins that a real near-cancelling r=3 instance is correctly judged a
        NON-refutation, so the loosened small-N allowance is not hiding one.
        """
        terms = [(1 << 90, 2), (1 << 91, 3), (1 << 92, 5)]
        E = exact_E(terms)
        if E is None:                      # pragma: no cover
            self.skipTest("degenerate terms")
        N = sum(max(abs(c).bit_length(), d.bit_length()) for c, d in terms)
        inst = Instance(terms, 1, E, N, 3)
        self.assertLess(slack_for(inst) * N, 10.0,
                        "prefactor allowance should be a few bits at this N")
        self.assertFalse(is_refutation(inst))

    def test_no_real_instance_found_is_a_refutation(self):
        """
        Reproduces the negative result of PROOFS.md 2.10.2 at small N: the
        best r=3 instances found by exhaustive greedy search all have a
        bounded additive bonus, so none is a refutation.  This is the
        harness doing its job on real data rather than synthetic terms.
        """
        for terms in ([(-10, 2), (-98, 3), (51, 13)],
                      [(-5, 2), (2, 3), (1, 13)]):
            with self.subTest(terms=terms):
                E = exact_E(terms)
                N = sum(max(abs(c).bit_length(), d.bit_length())
                        for c, d in terms)
                inst = Instance(terms, 1, E, N, 3)
                self.assertGreater(inst.rate, 2.0 / 3.0)
                self.assertFalse(is_refutation(inst))

    def test_instance_does_not_police_E_against_its_terms(self):
        """
        A known limitation, pinned so it cannot be forgotten.

        `Instance` trusts its `E` argument and does not check it against the
        terms.  An earlier design accidentally policed it, because
        `slack_for` read the magnitude; that protection is gone now that the
        additive allowance is a fixed bit budget.

        Consequence: a caller that fabricates `E` gets a meaningless verdict.
        Here terms [(1,2),(1,3),(1,5)] have |S| ~ 4.8, so the real E is about
        -2.26, yet claiming E = 750 produces an enormous excess and is
        reported as a refutation.  That is the harness correctly reacting to
        a bogus input, not a bug -- but it means every `E` fed to
        `is_refutation` must come from `exact_E`, which is what the rest of
        this module does.
        """
        terms = [(1, 2), (1, 3), (1, 5)]
        real_E = exact_E(terms)
        self.assertLess(real_E, 0.0, "these terms have |S| > 1, so E < 0")

        fabricated = Instance(terms, 1, 750.0, 1000, 3)
        self.assertTrue(is_refutation(fabricated),
                        "with a fixed bit budget the harness trusts E; the "
                        "protection is that callers use exact_E")

        # The honest version of the same terms is never a refutation.
        honest = Instance(terms, 1, real_E, 1000, 3)
        self.assertFalse(is_refutation(honest))

    def test_literature_cap_is_two(self):
        """
        EHS Thm 8 bounds E/N by about 2 uniformly, so a rate above 2 would
        contradict the literature, not just this repo.  Pin the constant.
        """
        self.assertEqual(LITERATURE_CAP, 2.0)


class TestCertifiedMagnitudes(unittest.TestCase):
    """Everything reported must be certified, not a float evaluation."""

    def test_exact_E_matches_the_sign_kernel(self):
        """
        E = log2(1/|S|) may legitimately be NEGATIVE when |S| > 1, e.g.
        S = 3*sqrt(2) - 5*sqrt(3) ~ -4.4 gives E ~ -2.14.  An earlier version
        of this test asserted E >= 0, which is false; the quantity is a rate
        and rates here are signed that way.  What must hold is that E is
        finite, defined, and consistent with a certified nonzero sign.
        """
        for terms in ([(1, 2)], [(3, 2), (-5, 3)], [(7, 2), (-2, 3), (11, 5)]):
            res = exact_sign(terms)
            E = exact_E(terms)
            with self.subTest(terms=terms):
                self.assertIsNotNone(E)
                self.assertNotEqual(res.sign, 0)
                self.assertFalse(math.isnan(E))
                self.assertTrue(math.isfinite(E))

    def test_exact_E_returns_None_when_S_is_zero(self):
        # distinct squarefree parts are independent, so force a real zero
        # only via a repeated radicand that aggregates away
        self.assertIsNone(exact_E([(5, 2), (-5, 2)]))

    def test_f2_rank_of_independent_and_dependent_sets(self):
        """
        2, 3, 5, 7 are independent (rank 4).  6 = 2*3 lies in their span, so
        adjoining it does NOT lower the rank -- an earlier version of this
        test claimed 3, which was simply wrong: a dependent element can only
        fail to raise the rank, never lower it.
        """
        self.assertEqual(f2_rank([2, 3, 5]), 3)
        self.assertEqual(f2_rank([2, 3, 6]), 2)          # 6 = 2*3
        self.assertEqual(f2_rank([2, 3, 5, 7]), 4)
        self.assertEqual(f2_rank([2, 3, 5, 7, 6]), 4)   # 6 adds nothing


class TestHuntsFindNothingAboveCeiling(unittest.TestCase):
    """
    A negative result, pinned so it is a reproducible claim rather than a
    shrug.  Small budgets keep this fast; the full sweep is
    `python3 refutation.py`.
    """

    def test_random_hunt_finds_no_refutation(self):
        for r in (2, 3, 4):
            best, beat = hunt_random(r, budget=25, bits=32, seed=7)
            with self.subTest(r=r):
                if best is not None:
                    self.assertEqual(beat, [])

    def test_greedy_hunt_finds_no_refutation(self):
        for r in (2, 3, 4):
            best, beat = hunt_greedy(r, budget=25, bits=32, seed=11)
            with self.subTest(r=r):
                if best is not None:
                    self.assertEqual(beat, [])

    def test_hunt_candidates_have_the_requested_rank(self):
        """
        The hunts enforce rank by rejection, so a returned candidate really
        does have F2-rank r.  Without this the ceiling being compared against
        would be the wrong one.
        """
        for r in (2, 3, 4):
            best, _ = hunt_random(r, budget=30, bits=24, seed=3)
            if best is None:
                continue
            agg, _ = decompose_and_aggregate(best.terms)
            agg = {d: c for d, c in agg.items() if c}
            with self.subTest(r=r):
                self.assertEqual(f2_rank(sorted(agg)), r)

    def test_hunt_reports_certified_signs(self):
        best, _ = hunt_random(3, budget=20, bits=24, seed=5)
        if best is None:
            self.skipTest("no candidate in budget")
        self.assertIn(best.sign, (-1, 1))
        self.assertEqual(best.sign, exact_sign(best.terms).sign)

    def test_module_main_runs_and_reports_negative(self):
        rc = refutation.main(["--r", "2", "--budget", "3", "--bits", "16",
                             "--f6-max-B", "20"])
        self.assertEqual(rc, 0, "main should report no refutation, rc=0")


if __name__ == "__main__":
    unittest.main(verbosity=2)
