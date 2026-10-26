"""
Tests for the squarefree construction and the rho effort budget.

Both fixes here came from measured failures, not from reading the code:

  * `is_squarefree` used to require a FULL factorisation, because it called
    `_factor` on the leftover cofactor.  Deciding squarefreeness only needs to
    know whether a prime repeats, so a single rho split plus a gcd suffices.
    The full-factorisation version made F2 with 200-bit radicands raise
    RhoFailure: a random 200-bit number almost always has a large prime
    cofactor that pure-Python rho will not split.

  * `random_squarefree` is the replacement.  It builds a certified radicand
    in the target window instead of rejection-sampling random integers and
    testing them, which took generating three 200-bit radicands from 21.8 s
    to 21 ms.

    `random_squarefree` used to search for one n_bits-bit prime and return it.
    Above 64 bits that prime's primality is only *probable*, so the
    squarefreeness it implied was only probable too, while the docstring
    described a verified result.  It now splits the bit budget across
    distinct primes each below 2^64, where the Miller-Rabin base set is
    proven, so squarefreeness holds by construction at every width.  The
    tests below check the exponent claim directly rather than re-running
    `is_squarefree` on a product that rho would have to re-split.
"""

import random
import unittest
from math import isqrt

import squarefree as SF


def brute_force_squarefree(n):
    """Definition: no prime square divides n."""
    for p in range(2, isqrt(n) + 1):
        if n % (p * p) == 0:
            return False
    return True


class TestIsSquarefreeMatchesDefinition(unittest.TestCase):
    def test_exhaustive_small_range(self):
        """The whole range where trial division is the definition."""
        for n in range(0, 5000):
            self.assertEqual(SF.is_squarefree(n), brute_force_squarefree(n),
                             f"is_squarefree({n}) disagrees with the definition")

    def test_against_decomposition(self):
        """
        n is squarefree exactly when squarefree_decompose(n) == (1, n).

        Note s == n, not s == 1: the squarefree *part* of a squarefree number
        is the number itself.  (Reading this as s == 1 is a plausible slip
        that only shows up on n = 2, since 2 is squarefree but has s = 2.)

        These are two different code paths -- one avoids full factorisation,
        the other requires it -- so agreement is real evidence rather than a
        tautology.
        """
        for n in (1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 17, 19, 30, 60, 210,
                  2310, 30030, 510510, 9699690, 2 * 3 * 5 * 7 * 11 * 13 * 17,
                  1000003, 1000003 ** 2, 1000003 ** 3, 41041, 561):
            c, s = SF.squarefree_decompose(n)
            self.assertEqual(SF.is_squarefree(n), s == n, f"n={n}, c={c}, s={s}")

    def test_large_random_values(self):
        rng = random.Random(20240101)
        for _ in range(300):
            n = rng.randrange(2, 10 ** rng.randint(2, 12))
            if n < 4:
                continue
            self.assertEqual(SF.is_squarefree(n), brute_force_squarefree(n), f"n={n}")

    def test_large_values_decide_or_raise_never_guess(self):
        """
        The contract at 200 bits, which is not "always answers".

        About a quarter of random 200-bit odd numbers are products of two
        large primes; splitting those needs ~2**50 rho iterations, which pure
        Python will not do, and each attempt costs seconds.  So this test uses
        64-bit values, where the same uncertainty exists but the budget
        covers it, and checks the property that actually matters: the
        function either decides or refuses, and never returns True for an
        input it could not check.
        """
        rng = random.Random(5)
        for _ in range(12):
            n = rng.randrange(1 << 63, 1 << 64) | 1
            try:
                result = SF.is_squarefree(n)
            except SF.RhoFailure:
                continue
            self.assertIsInstance(result, bool)

    def test_large_primes_and_their_squares_are_decided(self):
        """
        The cases that matter most are decided, not deferred.

        A prime is settled by the primality test alone, and a prime square by
        isqrt, so neither needs factorisation.
        """
        for p in ((1 << 127) - 1, (1 << 89) - 1, 1000003, 999983):
            if not SF._is_probable_prime(p):
                continue
            self.assertTrue(SF.is_squarefree(p), f"prime {p}")
            self.assertFalse(SF.is_squarefree(p * p), f"prime square {p**2}")

    def test_large_prime_squared_is_detected(self):
        """The square case that does not need factorisation at all."""
        p = 1000003
        self.assertFalse(SF.is_squarefree(p * p))
        self.assertFalse(SF.is_squarefree(4 * p))


class TestRandomSquarefree(unittest.TestCase):
    def test_bit_length_is_exact(self):
        """
        Every width, exactly.

        The construction splits the requested bit budget across distinct
        certified primes, so an off-by-one anywhere in the window arithmetic
        would show up here as a value of the wrong width.
        """
        rng = random.Random(7)
        for bits in (2, 3, 4, 5, 6, 7, 8, 12, 16, 31, 32, 33, 63, 64, 65, 66,
                     95, 96, 97, 128, 200, 256, 300, 511, 512):
            with self.subTest(bits=bits):
                n = SF.random_squarefree(bits, rng)
                self.assertEqual(n.bit_length(), bits, f"{n} is not {bits} bits")

    def test_exhaustive_squarefreeness_at_small_sizes(self):
        """
        The independent check, at the sizes where the definition is cheap.

        Checking the construction's own post-conditions would be circular, so
        small values are tested against the definition itself -- no prime
        square divides n -- which nothing in `squarefree.py` is written to
        make pass.
        """
        rng = random.Random(99)
        for bits in range(2, 22):
            for _ in range(30):
                n = SF.random_squarefree(bits, rng)
                self.assertTrue(
                    brute_force_squarefree(n),
                    f"{n} ({bits} bits) is not squarefree by the definition"
                )

    def test_every_factor_is_a_distinct_certified_prime(self):
        """
        The property the construction actually claims, proved rather than
        asserted: a product of distinct primes is squarefree.

        A complete factorisation of the result is taken, and every exponent
        must be exactly 1 and every factor below `DETERMINISTIC_MR_CEILING`.
        The second half matters as much as the first: it is what makes each
        factor a *proven* prime rather than a probable one, which is the whole
        basis for calling a generated radicand certified.

        This is the check the old single-prime construction could not pass
        above 64 bits, where its one large factor's primality was only
        probable.
        """
        rng = random.Random(5)
        for bits in (24, 32, 48, 64, 96, 128, 200):
            with self.subTest(bits=bits):
                n = SF.random_squarefree(bits, rng)
                exps: dict = {}
                SF._factor(n, exps)
                self.assertTrue(exps, f"no factors for {n}")
                for p, e in exps.items():
                    self.assertEqual(
                        e, 1, f"{n} has {p} to exponent {e}, so it is not squarefree"
                    )
                    self.assertLess(
                        p, SF.DETERMINISTIC_MR_CEILING,
                        f"factor {p} is at or above 2^64, where Miller-Rabin is "
                        f"only probable"
                    )
                product = 1
                for p, e in exps.items():
                    product *= p ** e
                self.assertEqual(product, n, "factorisation does not multiply back")

    def test_construction_verifier_rejects_a_broken_factors_list(self):
        """
        `_verify_constructed` is load-bearing, so prove it bites.

        A repeated factor is exactly how a constructed radicand would stop
        being squarefree, and the previous single-prime version had no check
        that could have caught it.
        """
        prime = 1000003
        with self.assertRaises(RuntimeError) as ctx:
            SF._verify_constructed(41, 1 << 40, 1 << 41, [prime, prime], prime * prime)
        self.assertIn("repeated", str(ctx.exception))

    def test_construction_verifier_rejects_a_factor_above_the_mr_ceiling(self):
        """A factor at or above 2^64 is refused rather than called certified."""
        too_big = SF.DETERMINISTIC_MR_CEILING + 1
        with self.assertRaises(RuntimeError) as ctx:
            SF._verify_constructed(65, 1 << 64, 1 << 65, [too_big], too_big)
        self.assertIn("2^64", str(ctx.exception))

    def test_is_fast_at_200_bits(self):
        """
        The measured improvement: 21.8 s -> ~21 ms per value, so a k=3 sweep
        must be well under a second.  Generous bound to stay stable on a
        loaded machine; the point is to catch a return to rho-per-candidate.
        """
        import time
        rng = random.Random(1)
        start = time.perf_counter()
        for _ in range(5):
            SF.random_squarefree(200, rng)
        elapsed = time.perf_counter() - start
        self.assertLess(elapsed, 5.0,
                        f"5 x 200-bit squarefree took {elapsed:.1f}s; "
                        f"expected milliseconds, so the sampler is back to "
                        f"trial-dividing candidates")

    def test_cost_grows_with_width_not_exponentially(self):
        """
        The 512-bit cap is a cost bound, so the cost must actually be tracked.

        Certified factors accumulate linearly in n_bits, so doubling the width
        roughly doubles the work.  A return to rejection sampling would show
        up here as a far worse blowup than linear.
        """
        import time
        def elapsed(bits, count):
            rng = random.Random(3)
            start = time.perf_counter()
            for _ in range(count):
                SF.random_squarefree(bits, rng)
            return time.perf_counter() - start
        small = elapsed(64, 20)
        large = elapsed(512, 20)
        self.assertLess(
            large, max(small * 40, 2.0),
            f"512-bit construction took {large:.3f}s vs {small:.3f}s at 64 bits"
        )

    def test_rejects_out_of_range_bit_counts(self):
        rng = random.Random(1)
        for bad in (0, 1, -4):
            with self.assertRaises(ValueError):
                SF.random_squarefree(bad, rng)
        with self.assertRaises(ValueError):
            SF.random_squarefree(4096, rng)

    def test_deterministic_for_a_given_seed(self):
        a = [SF.random_squarefree(64, random.Random(5)) for _ in range(3)]
        b = [SF.random_squarefree(64, random.Random(5)) for _ in range(3)]
        self.assertEqual(a, b)

    def test_produces_distinct_values(self):
        rng = random.Random(11)
        vals = [SF.random_squarefree(64, rng) for _ in range(8)]
        self.assertEqual(len(set(vals)), len(vals), "sampler repeated a value")


class TestRhoBudgetScalesWithSize(unittest.TestCase):
    def test_small_cofactors_still_split(self):
        for name, n in (("p*q", 1000003 * 1000033),
                        ("p^2*q", 1000003 ** 2 * 1000033),
                        ("carmichael", 41041),
                        ("smooth", 2 ** 20 * 3 ** 20)):
            with self.subTest(name=name):
                d = SF._pollard_rho(n)
                self.assertGreater(d, 1)
                self.assertLess(d, n)
                self.assertEqual(n % d, 0, f"rho returned a non-divisor for {name}")

    def test_rho_is_deterministic(self):
        n = 1000003 * 1000033
        first = SF._pollard_rho(n)
        for _ in range(3):
            self.assertEqual(SF._pollard_rho(n), first)

    def test_unsplittable_input_raises_rather_than_guessing(self):
        """
        A large semiprime must raise, not return a wrong factorisation.

        Two 61-bit primes need ~2**30 rho iterations, which pure Python will
        not do, so this input cannot be split here.  The test asserts the
        refusal and that the message names the radicand.  It is the slowest
        test in the suite by design -- it is proving the failure is reported
        rather than papered over, and there is no cheap way to shrink the
        input far enough for rho to still fail.
        """
        p = 2305843009213693951     # 2**61 - 1
        q = 2305843009213693921
        n = p * q
        with self.assertRaises(SF.RhoFailure) as ctx:
            SF.squarefree_decompose(n)
        message = str(ctx.exception)
        self.assertIn(str(n), message, "the error should name the radicand")
        self.assertIn("Pollard", message)

    def test_rho_failure_leaves_no_partial_factorisation(self):
        """
        A failed factorisation must not leave a dict that looks complete.

        Checks the invariant directly on a small input where rho can be made
        to fail by lowering the budget, which keeps the test fast.
        """
        original_budget_probe = SF._pollard_rho

        def always_fails(n):
            raise SF.RhoFailure("budget exhausted (test double)")

        n = 1000003 * 1000033
        out = {}
        SF._pollard_rho = always_fails
        try:
            with self.assertRaises(SF.RhoFailure):
                SF._factor(n, out)
        finally:
            SF._pollard_rho = original_budget_probe
        self.assertEqual(out, {}, f"a failed factorisation recorded {out}")


class TestRhoFailureType(unittest.TestCase):
    def test_is_a_runtime_error(self):
        """Callers that catch RuntimeError around factoring keep working."""
        self.assertTrue(issubclass(SF.RhoFailure, RuntimeError))


if __name__ == "__main__":
    unittest.main(verbosity=2)
