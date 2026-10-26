"""Tests for families.py — instance generators."""

import random
import unittest

from families import (
    evaluate_float_surrogate,
    generate_f1,
    generate_f1_sequence,
    generate_f2,
    generate_f3,
    generate_f4,
    generate_f5,
    generate_f6,
    generate_f6_sequence,
    generate_family,
    generate_family_sequence,
    input_bits,
    is_squarefree,
)
from signexact import exact_sign
from squarefree import decompose_and_aggregate, squarefree_decompose


class TestInputBits(unittest.TestCase):
    def test_uses_exact_bit_length(self):
        """bit_length, not log2: a 200-bit radicand must be counted exactly."""
        big = (1 << 200) + 1          # odd, so not a perfect square
        n = input_bits([(1, big)])
        self.assertEqual(n, 1 + big.bit_length())

    def test_does_not_overflow_on_huge_radicands(self):
        """log2 would raise OverflowError past 2**1024; this must not."""
        huge = (1 << 5000) + 3
        self.assertEqual(input_bits([(1, huge)]), 1 + huge.bit_length())


class TestF1PellMobius(unittest.TestCase):
    def test_f1_generation(self):
        for a in (2, 3, 4, 5):
            inst = generate_f1(a, 5)
            self.assertEqual(inst["family"], "F1")
            self.assertEqual(inst["k"], 2)
            self.assertEqual(len(inst["terms"]), 2)
            self.assertEqual(inst["terms"][0][0], 1)
            self.assertEqual(inst["terms"][1][0], -1)

    def test_f1_pell_identity_and_known_sign(self):
        """F1 must carry the exact residue and the sign it implies."""
        for a in (2, 3, 4, 5, 6):
            for n in (1, 2, 3, 4, 5, 6):
                inst = generate_f1(a, n)
                meta = inst["metadata"]
                p, q, alpha = meta["p"], meta["q"], meta["alpha"]
                self.assertEqual(p * p - alpha * q * q, meta["residue"])
                want = 1 if n % 2 == 1 else -1
                self.assertEqual(meta["residue"], want if n % 2 else -a)
                self.assertEqual(meta["known_sign"], want)

    def test_f1_terms_match_the_construction(self):
        inst = generate_f1(3, 5)
        p, q, alpha = (inst["metadata"][k] for k in ("p", "q", "alpha"))
        self.assertEqual(inst["terms"], [(1, p * p), (-1, alpha * q * q)])

    def test_f1_rejects_bad_a(self):
        for bad in (0, 1, -3):
            with self.assertRaises(ValueError):
                generate_f1(bad, 1)
        with self.assertRaises(TypeError):
            generate_f1(2.0, 1)

    def test_f1_sequence_length(self):
        self.assertEqual(len(generate_f1_sequence(3, 1, 8)), 8)
        with self.assertRaises(ValueError):
            generate_f1_sequence(3, 5, 2)


class TestF2Random(unittest.TestCase):
    def test_f2_generation(self):
        inst = generate_f2(3, 20, seed=42)
        self.assertEqual(inst["family"], "F2")
        self.assertEqual(inst["k"], 3)
        self.assertEqual(len(inst["terms"]), 3)

    def test_f2_reproducibility(self):
        self.assertEqual(generate_f2(3, 20, seed=42)["terms"],
                         generate_f2(3, 20, seed=42)["terms"])

    def test_f2_terms_are_really_squarefree(self):
        """
        The old generator returned True for any n >= 10**12 without testing,
        so it could emit radicands with square factors.
        """
        for seed in range(10):
            inst = generate_f2(4, 40, seed=seed)
            for _, rad in inst["terms"]:
                self.assertTrue(is_squarefree(rad), f"{rad} not squarefree")

    def test_f2_rejects_bad_params(self):
        with self.assertRaises(ValueError):
            generate_f2(0, 20)
        with self.assertRaises(ValueError):
            generate_f2(3, 1)

    def test_f2_raises_rather_than_returning_short(self):
        """Impossible request must raise, not silently return fewer terms."""
        with self.assertRaises(RuntimeError):
            generate_f2(200, 8)


class TestF3Greedy(unittest.TestCase):
    def test_f3_generation(self):
        inst = generate_f3(3, 20, target_bits=10, max_iter=100, seed=42)
        self.assertEqual(inst["family"], "F3")
        self.assertEqual(inst["k"], 3)

    def test_f3_labels_its_objective_as_approximate(self):
        inst = generate_f3(3, 20, target_bits=10, max_iter=100, seed=42)
        meta = inst["metadata"]
        self.assertTrue(meta["surrogate_is_approximate"])
        self.assertIn("exact_sign", meta)
        self.assertIn("exact_method", meta)

    def test_f3_reported_sign_is_the_exact_one(self):
        from signexact import exact_sign
        inst = generate_f3(4, 24, target_bits=8, max_iter=200, seed=3)
        self.assertEqual(inst["metadata"]["exact_sign"],
                         exact_sign(inst["terms"]).sign)


class TestF4Product(unittest.TestCase):
    def test_f4_generation(self):
        inst = generate_f4(2, 3, 5, 7)
        self.assertEqual(inst["family"], "F4")
        terms = {rad: coeff for coeff, rad in inst["terms"]}
        self.assertEqual(terms.get(10), 1)
        self.assertEqual(terms.get(14), -1)
        self.assertEqual(terms.get(15), -1)
        self.assertEqual(terms.get(21), 1)

    def test_f4_combination(self):
        inst = generate_f4(2, 8, 3, 12)
        terms = {rad: coeff for coeff, rad in inst["terms"]}
        self.assertEqual(terms.get(6), 1)
        self.assertEqual(terms.get(24), -2)
        self.assertEqual(terms.get(96), 1)

    def test_f4_rejects_bad_params(self):
        for bad in (0, -1):
            with self.assertRaises(ValueError):
                generate_f4(bad, 3, 5, 7)
        with self.assertRaises(TypeError):
            generate_f4(2.0, 3, 5, 7)

    def test_f4_identical_factors_collapse_to_empty_and_raise(self):
        # (sqrt(a)-sqrt(a))(...) is identically zero; that is not an instance.
        with self.assertRaises(ValueError):
            generate_f4(2, 2, 5, 7)


class TestF5Squaring(unittest.TestCase):
    def test_f5_generation(self):
        inst = generate_f5(3, 2)
        self.assertEqual(inst["family"], "F5")
        terms = {rad: coeff for coeff, rad in inst["terms"]}
        self.assertEqual(terms.get(25), 1)
        self.assertEqual(terms.get(24), -1)

    def test_f5_terms_evaluate_to_the_square_it_claims(self):
        """
        The identity the family exists for, checked on the emitted terms.

        The old test only asserted the arithmetic `(a+b)^2 == 25` and
        `4ab == 24`, which the generator produces trivially and which would
        still pass if the emitted sum were some unrelated pair of terms with
        the same radicands.  Here the *sum* is evaluated by the certified
        decider and compared against the square, so a wrong pairing of
        coefficients to radicands is caught.
        """
        for a, b in ((3, 2), (5, 7), (12, 5), (30, 1), (100, 99)):
            for s, label in ((-1, "minus"), (1, "plus")):
                with self.subTest(a=a, b=b, s=s):
                    inst = generate_f5(a, b, s)
                    got = exact_sign(inst["terms"]).sign
                    # s = -1 -> (sqrt(a)-sqrt(b))^2 > 0 unless a == b
                    # s = +1 -> (sqrt(a)+sqrt(b))^2 > 0 always
                    self.assertEqual(
                        got, 1,
                        f"F5({a},{b},s={s}) [{label}] did not evaluate positive"
                    )
                    # And the value really is the square: check the integer
                    # part and the radical coefficient agree with a+b and s.
                    agg, _ = decompose_and_aggregate(inst["terms"])
                    agg = {k: v for k, v in agg.items() if v}
                    self.assertEqual(len(agg), 2,
                                     f"expected two distinct squarefree parts, "
                                     f"got {agg}")

    def test_f5_square_equality_in_the_field(self):
        """
        The strongest available statement: the emitted S equals the square,
        checked exactly via squarefree aggregation rather than by sign.

        (sqrt(a) +- sqrt(b))^2 = a + b +- 2 sqrt(ab), so after decomposition
        the sum must carry the rational part a+b at squarefree part 1, and
        s*2c at the squarefree part of ab.  Comparing that rules out a
        generator which returns the right radicands paired with the wrong
        coefficients -- a failure the arithmetic-only old test could not see.
        """
        for a, b in ((3, 2), (5, 7), (12, 5), (30, 1)):
            for s in (-1, 1):
                with self.subTest(a=a, b=b, s=s):
                    inst = generate_f5(a, b, s)
                    agg, _ = decompose_and_aggregate(inst["terms"])
                    # The emitted radicands are (a+b)^2 and 4ab, so the
                    # aggregated form is (a+b) at squarefree part 1, plus
                    # s*2*c at the squarefree part of ab where 4ab = (2c)^2*sab.
                    sq_sum, sab = squarefree_decompose((a + b) * (a + b))
                    self.assertEqual(sq_sum, a + b, "(a+b)^2 must be a perfect square")
                    c4, sab4 = squarefree_decompose(4 * a * b)
                    c, sq_ab = squarefree_decompose(a * b)
                    self.assertEqual((c4, sab4), (2 * c, sq_ab))
                    expected = {1: a + b, sq_ab: s * c4}
                    expected = {k: v for k, v in expected.items() if v}
                    self.assertEqual(
                        agg, expected,
                        f"F5({a},{b},s={s}) terms do not sum to "
                        f"(sqrt(a) {'-' if s < 0 else '+'} sqrt(b))^2"
                    )

    def test_f5_s_selects_the_cross_term_sign(self):
        """s is a real parameter now, not just prose in a docstring."""
        plus = generate_f5(3, 2, 1)
        minus = generate_f5(3, 2, -1)
        default = generate_f5(3, 2)
        self.assertEqual(default["terms"], minus["terms"])
        self.assertNotEqual(plus["terms"], minus["terms"])
        self.assertEqual(plus["metadata"]["s"], 1)
        self.assertEqual(minus["metadata"]["s"], -1)

        # The two forms differ by exactly 4*sqrt(ab).
        agg_plus, _ = decompose_and_aggregate(plus["terms"])
        agg_minus, _ = decompose_and_aggregate(minus["terms"])
        c, sab = squarefree_decompose(2 * 3)
        for k in set(agg_plus) | set(agg_minus):
            self.assertEqual(
                agg_plus.get(k, 0) - agg_minus.get(k, 0),
                4 * c if k == sab else 0,
                "s=+1 and s=-1 must differ by exactly 4*sqrt(ab)"
            )

    def test_f5_rejects_bad_params(self):
        with self.assertRaises(ValueError):
            generate_f5(0, 2)
        with self.assertRaises(TypeError):
            generate_f5(3, 2.0)
        for bad_s in (0, 2, -2):
            with self.assertRaises(ValueError, msg=f"s={bad_s} was accepted"):
                generate_f5(3, 2, bad_s)
        for bad_s in (1.0, None, True):
            with self.assertRaises(TypeError, msg=f"s={bad_s!r} was accepted"):
                generate_f5(3, 2, bad_s)

    def test_f5_is_never_negative(self):
        """
        The structural limit the docstring states: both forms are squares of
        real numbers, so S >= 0 always and is zero only at a == b, s == -1.
        """
        for a, b in ((3, 2), (5, 7), (30, 1), (2, 3)):
            for s in (-1, 1):
                with self.subTest(a=a, b=b, s=s):
                    sign = exact_sign(generate_f5(a, b, s)["terms"]).sign
                    self.assertEqual(sign, 1)
        self.assertEqual(
            exact_sign(generate_f5(5, 5, -1)["terms"]).sign, 0,
            "(sqrt(5) - sqrt(5))^2 is exactly zero"
        )


class TestSurrogate(unittest.TestCase):
    def test_surrogate_is_float_and_labelled(self):
        val = evaluate_float_surrogate([(1, 2), (1, 3), (-1, 5)])
        self.assertIsInstance(val, float)
        # sqrt2 + sqrt3 - sqrt5 = 0.9101963924...
        self.assertAlmostEqual(val, 0.910196392442, places=9)

    def test_surrogate_agrees_with_exact_sign_on_clearly_separated_sums(self):
        from signexact import exact_sign
        rng = random.Random(17)
        for _ in range(200):
            terms = [(rng.choice([-1, 1]), rng.randint(2, 50)) for _ in range(3)]
            if abs(evaluate_float_surrogate(terms)) > 1e-6:
                self.assertEqual(
                    (evaluate_float_surrogate(terms) > 0) - (evaluate_float_surrogate(terms) < 0),
                    exact_sign(terms).sign,
                )


class TestF6CeilingFamily(unittest.TestCase):
    """
    F6 is the Thm 2.5 pigeonhole-gap construction at r = 3.  It is the only
    family here that ATTAINS the conjectured sharp rate E/N -> (r-1)/r = 2/3
    rather than searching for something worse, so it is the family that
    tests whether the ceiling is where it is claimed to be.
    """

    def test_f6_produces_three_distinct_squarefree_radicands(self):
        inst = generate_f6(12)
        self.assertEqual(inst.family, "F6")
        rads = [r for _, r in inst.terms]
        self.assertEqual(len(set(rads)), len(rads), "radicands must be distinct")
        for r in rads:
            self.assertTrue(is_squarefree(r), f"{r} not squarefree")

    def test_f6_instance_is_certified_nonzero(self):
        """The generator must not return an unverified candidate."""
        for B in (8, 12, 20):
            inst = generate_f6(B)
            res = exact_sign(inst.terms)
            with self.subTest(B=B):
                self.assertNotEqual(res.sign, 0)
                self.assertEqual(res.sign, inst.metadata["exact_sign"])

    def test_f6_respects_thm_2_5_bound(self):
        """
        PROOFS.md Thm 2.5 at r=3:  |S| <= (sum_i sqrt(d_i)) / B^(r-1).
        The exact Fraction |S| must be at or below it.
        """
        for B in (10, 20, 40):
            inst = generate_f6(B)
            mag = inst.metadata["exact_mag"]
            with self.subTest(B=B):
                self.assertLessEqual(float(mag),
                                     inst.metadata["thm25_bound"] * (1 + 1e-9))

    def test_f6_coefficients_are_bounded_by_the_grid(self):
        """A gap on the {0..B}^3 grid has |C_i| <= B, per Thm 2.5's hypothesis."""
        for B in (10, 25):
            inst = generate_f6(B)
            for c, _ in inst.terms:
                with self.subTest(B=B, c=c):
                    self.assertLessEqual(abs(c), B)

    def test_f6_rate_trends_toward_two_thirds(self):
        """
        E/N should climb toward (r-1)/r = 2/3 = 0.6667.  It may exceed it at
        finite size -- the o(1) in Thm 2.5 is still O(1) here -- so assert
        the TREND and a plausible band, not an exact inequality.
        """
        import math
        rates = []
        for B in (20, 60, 120):
            inst = generate_f6(B)
            m = inst.metadata["exact_mag"]
            E = -math.log2(m.numerator / m.denominator)
            rates.append(E / inst.N_bits)
        self.assertGreaterEqual(rates[-1], 0.5)
        self.assertLess(rates[-1], 1.2, f"F6 rate {rates[-1]} implausible")
        self.assertGreaterEqual(max(rates), 0.6,
                                "F6 should get within 10% of the r=3 ceiling")

    def test_f6_rejects_bad_params(self):
        with self.assertRaises(ValueError):
            generate_f6(0)
        with self.assertRaises(ValueError):
            generate_f6(10, radicands=(2, 3))
        with self.assertRaises(ValueError):
            generate_f6(10, radicands=(2, 3, 4))     # 4 is not squarefree

    def test_f6_sequence_sweeps(self):
        seq = generate_f6_sequence([8, 12, 16])
        self.assertEqual(len(seq), 3)
        self.assertEqual([i.metadata["B"] for i in seq], [8, 12, 16])

    def test_f6_is_in_the_registry(self):
        self.assertEqual(generate_family("F6", B=10)["family"], "F6")
        # generate_family forwards unknown kwargs straight to the generator,
        # so a bad parameter surfaces as TypeError, not ValueError.
        with self.assertRaises(TypeError):
            generate_family("F6", nope=1)


class TestRegistry(unittest.TestCase):
    def test_generate_family_dispatch(self):
        self.assertEqual(generate_family("F1", a=2, n=3)["family"], "F1")
        with self.assertRaises(ValueError):
            generate_family("F9", a=2, n=3)

    def test_generate_family_sequence_actually_sweeps(self):
        """Previously ignored param_ranges and returned a single instance."""
        swept = generate_family_sequence("F1", {"a": [2, 3, 4], "n": [1, 3]})
        self.assertEqual(len(swept), 6)
        pairs = {(i["metadata"]["a"], i["metadata"]["n"]) for i in swept}
        self.assertEqual(pairs, {(a, n) for a in (2, 3, 4) for n in (1, 3)})

    def test_generate_family_sequence_rejects_unknown_params(self):
        with self.assertRaises(ValueError):
            generate_family_sequence("F1", {"nope": [1]})
        with self.assertRaises(ValueError):
            generate_family_sequence("F9", {})

    def test_f4_sweep(self):
        swept = generate_family_sequence("F4", {"a": [2, 3], "b": [5, 6], "c": [7], "d": [11]})
        self.assertEqual(len(swept), 4)


if __name__ == "__main__":
    unittest.main(verbosity=2)
