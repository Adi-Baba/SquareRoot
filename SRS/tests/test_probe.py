"""
Tests for probe.py — the experiment harness.

This file did not exist.  probe.py is about 40% of the SRS code by line count
and had no runtime test at all, which left its most important invariant
unexercised: `probe_instance` asserts that the algebraic decider and the
certified bracket agree, and raises if they do not.  That assertion is the
harness's whole reason for existing, and nothing checked it fired.

Also covered here, because each was a flag that did not do what it said:

  * `--max-bits` was parsed, documented, and then never passed to either
    decider, so the cap was always `signexact.DEFAULT_MAX_BITS`.
  * All four families shared one try block, so the first instance that hit the
    precision cap aborted every family after it.
  * `main` computed an `undecided` count, discarded it, and returned 0 either
    way, so an undecided run was indistinguishable from a clean one.
  * F2 and F3 produced identical `instance_id` strings on every default run.
"""

import csv
import io
import os
import tempfile
import unittest
from contextlib import redirect_stdout

import probe
from families import generate_f1, generate_f2, generate_f3, generate_f4, generate_f5
from signexact import SignUnknownError


def quiet(fn, *args, **kwargs):
    """Run something that prints, and return its value with stdout captured."""
    buf = io.StringIO()
    with redirect_stdout(buf):
        result = fn(*args, **kwargs)
    return result, buf.getvalue()


class TestProbeInstance(unittest.TestCase):
    def test_each_family_probes(self):
        """One instance per family, decided and bracketed."""
        from fractions import Fraction as F
        runner = probe.ProbeRunner()
        for inst in (generate_f1(3, 5), generate_f2(3, 16, seed=7),
                     generate_f3(3, 16, seed=7), generate_f4(2, 3, 5, 7),
                     generate_f5(3, 2), generate_f5(3, 2, 1)):
            with self.subTest(family=inst["family"]):
                r = runner.probe_instance(inst)
                self.assertIn(r.exact_sign, (-1, 0, 1))
                # The bracket endpoints are exact rationals serialised as
                # strings, so they are compared as Fractions -- comparing the
                # strings would be a meaningless lexicographic test.
                self.assertLessEqual(F(r.certified_lo), F(r.certified_hi))
                self.assertGreater(r.N_bits, 0)

    def test_f1_sign_matches_the_closed_form_residue_law(self):
        """
        F1's sign is known in closed form, so the harness can be checked
        against a source that does not use the harness at all.
        """
        import f1_theorem as FT
        runner = probe.ProbeRunner()
        for a in (2, 3, 4):
            for n in (1, 2, 3, 6, 7):
                with self.subTest(a=a, n=n):
                    r = runner.probe_instance(generate_f1(a, n))
                    self.assertEqual(
                        r.exact_sign, FT.expected_sign(a, n),
                        f"probe disagreed with the closed form for a={a}, n={n}"
                    )

    def test_bracket_clears_zero_exactly_when_the_sign_is_nonzero(self):
        """
        `bracket_clears_zero` means "the bracket excludes zero", so it is
        False for an exact zero.  The sign itself lives in `exact_sign`, and
        conflating the two would report an exact answer as an unknown one.
        """
        runner = probe.ProbeRunner()
        # (sqrt(12) - 2 sqrt(3)) == 0 exactly, by exact cancellation.
        inst = {"terms": [(1, 12), (-2, 3)], "family": "F1", "k": 2,
                "N_bits": 4, "metadata": {"a": 0, "n": 0}}
        r = runner.probe_instance(inst)
        self.assertEqual(r.exact_sign, 0)
        self.assertFalse(r.bracket_clears_zero)

    def test_sign_disagreement_is_a_loud_failure(self):
        """
        The invariant: the algebraic decider and the certified bracket must
        agree, and disagreement is a bug rather than a table row.

        `certified_interval_for_S` is stubbed to return the opposite sign,
        which is the one thing that must never pass silently.
        """
        runner = probe.ProbeRunner()
        original = probe.certified_interval_for_S
        try:
            def wrong_sign(terms, digits=50, max_bits=None):
                iv = original(terms, digits, max_bits)
                return iv._replace(sign=-iv.sign if iv.sign else 1)
            probe.certified_interval_for_S = wrong_sign
            with self.assertRaises(AssertionError) as ctx:
                runner.probe_instance(generate_f4(2, 3, 5, 7))
            self.assertIn("sign disagreement", str(ctx.exception))
        finally:
            probe.certified_interval_for_S = original
        # And the real path still works.
        runner.probe_instance(generate_f4(2, 3, 5, 7))

    def test_certificate_width_is_a_plausible_bit_estimate(self):
        """
        `certificate_width` is a bit-length *estimate* of log2(1/width), not
        an exact count.  It must at least be a sane non-negative number that
        tracks the actual width.
        """
        from fractions import Fraction as F
        import math
        runner = probe.ProbeRunner()
        r = runner.probe_instance(generate_f1(3, 8))
        self.assertIsNotNone(r.certificate_width)
        self.assertGreater(r.certificate_width, 0)
        self.assertGreaterEqual(r.certificate_width_log2, 0.0)
        # The estimate should be within about one bit of the true value.
        true_bits = -math.log2(float(F(r.certified_hi) - F(r.certified_lo)))
        self.assertLess(abs(r.certificate_width - true_bits), 1.5)


class TestMaxBitsIsThreaded(unittest.TestCase):
    def test_runner_forwards_max_bits_to_the_deciders(self):
        """A tiny cap must actually cause a refusal, not be ignored."""
        tight = probe.ProbeRunner(max_bits=8)
        with self.assertRaises(SignUnknownError):
            tight.probe_instance(generate_f1(3, 6))

    def test_default_runner_is_unaffected(self):
        # n = 6 is even, so the residue law gives S < 0.
        loose = probe.ProbeRunner()
        self.assertEqual(loose.probe_instance(generate_f1(3, 6)).exact_sign, -1)

    def test_max_bits_none_means_the_default(self):
        """
        `--max-bits` absent passes None through to both deciders, so both
        have to read None as "use the default".  `exact_sign` used to reject
        None with a TypeError while `certified_interval_for_S` accepted it.
        """
        runner = probe.ProbeRunner(max_bits=None)
        self.assertIsNone(runner.max_bits)
        self.assertIn(runner.probe_instance(generate_f1(3, 6)).exact_sign,
                      (-1, 0, 1))

    def test_exact_sign_accepts_none_for_max_bits(self):
        from signexact import exact_sign
        self.assertEqual(exact_sign([(1, 12), (-2, 3)], max_bits=None).sign, 0)


class TestMainReporting(unittest.TestCase):
    def test_clean_run_exits_zero(self):
        rc, out = quiet(probe.main, ["--family", "F1"])
        self.assertEqual(rc, 0, out)
        self.assertNotIn("stopped on the precision cap", out)

    def test_capped_run_exits_nonzero(self):
        """An undecided run must not look like a clean one."""
        rc, out = quiet(probe.main, ["--max-bits", "8"])
        self.assertEqual(rc, 1, out)
        self.assertIn("stopped on the precision cap", out)

    def test_one_capped_family_does_not_abort_the_others(self):
        """
        The families used to share a single try block, so a refusal in the
        first one silently skipped every family after it.  With a cap of 8
        bits all four refuse, and each must be reported on its own.
        """
        rc, out = quiet(probe.main, ["--max-bits", "8"])
        for family in ("F1", "F2", "F4", "F5"):
            self.assertIn(
                f"{family}: instance beyond the precision cap", out,
                f"{family} was not reported as capped"
            )
        self.assertIn("4 suite(s) stopped", out)

    def test_family_filter_selects_one_suite(self):
        rc, out = quiet(probe.main, ["--family", "F4"])
        self.assertEqual(rc, 0, out)
        self.assertIn("F4 (product composition)", out)
        self.assertNotIn("F1 (Pell", out)

    def test_csv_output_is_written_and_reparses(self):
        from fractions import Fraction as F
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "out.csv")
            rc, out = quiet(probe.main, ["--family", "F1", "--output", path])
            self.assertEqual(rc, 0, out)
            self.assertTrue(os.path.isfile(path))
            with open(path, newline="") as fh:
                rows = list(csv.DictReader(fh))
            self.assertTrue(rows)
            for row in rows:
                self.assertIn(int(row["exact_sign"]), (-1, 0, 1))
                # Endpoints survive the round trip as exact rationals, which
                # is the point of serialising them as text rather than float.
                self.assertLessEqual(F(row["certified_lo"]), F(row["certified_hi"]))


class TestInstanceIds(unittest.TestCase):
    def test_f2_and_f3_ids_do_not_collide(self):
        """
        Both families describe themselves by (k, seed) and the default suite
        runs both over the same grid, so the ids were identical on every
        default run.  Uniqueness is asserted on a real run.
        """
        runner = probe.ProbeRunner()
        quiet(runner.run_random_suite, [2], [16], num_instances=3)
        ids = [r.instance_id for r in runner.results]
        self.assertEqual(len(ids), len(set(ids)),
                         f"duplicate instance_id in a default run: {ids}")
        self.assertTrue(any(i.startswith("F2_") for i in ids), ids)
        self.assertTrue(any(i.startswith("F3_") for i in ids), ids)

    def test_f5_id_records_the_sign(self):
        """Two F5 instances differing only in s must not share an id."""
        a = probe.ProbeRunner()._instance_id("F5", {"a": 3, "b": 2, "s": 1})
        b = probe.ProbeRunner()._instance_id("F5", {"a": 3, "b": 2, "s": -1})
        self.assertNotEqual(a, b)

    def test_unknown_family_is_marked(self):
        self.assertEqual(probe.ProbeRunner()._instance_id("F9", {}), "?")


class TestSuitesRunForReal(unittest.TestCase):
    def test_run_family_returns_one_result_per_combination(self):
        """
        `run_family` used to be a stub returning [].  It must return a real
        row per point of the Cartesian product.
        """
        runner = probe.ProbeRunner()
        results, _ = quiet(runner.run_family, "F1", {"a": [2, 3], "n": [1, 2, 3]})
        self.assertEqual(len(results), 6)
        self.assertEqual(len(runner.results), 6)

    def test_run_f1_sequence_honours_its_range(self):
        """The range argument used to be ignored in favour of range(1, 10)."""
        runner = probe.ProbeRunner()
        results, _ = quiet(runner.run_f1_sequence, [2], (3, 5))
        self.assertEqual(len(results), 3)
        self.assertEqual([r.instance_id for r in results],
                         ["a2_n3", "a2_n4", "a2_n5"])

    def test_print_table_survives_an_empty_run(self):
        out, text = quiet(probe.ProbeRunner().print_results_table)
        self.assertIn("no results", text)


if __name__ == "__main__":
    unittest.main(verbosity=2)
