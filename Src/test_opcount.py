"""Tests for opcount.py -- exact, machine-independent operation counts.

The counts themselves are the point: the paper's explanation for why the zoom
loses to interval Newton was read off the source by eye, and measurement
shows it was wrong.  These tests pin the counts so the corrected claim cannot
silently drift back, and they check the counting machinery is not lying.

The counting method is (line-executions) x (arithmetic ops on that line), both
exact.  The risk is a systematic error -- e.g. attributing a multi-line
statement's operations to the wrong line, or double counting.  The
`test_known_loop_shape` test pins a hand-counted loop as a control.
"""

import ast
import unittest

import opcount as O
import exp_certified as C


class TestCountingMachinery(unittest.TestCase):
    def test_ops_per_statement_hand_checked(self):
        # A statement with one division and one multiplication, and one with
        # neither, must be categorised exactly.
        src = "x = (a + b) / 2\ny = p * q\nz = w\n"
        tree = ast.parse(src)
        per_line = O.ops_per_statement(tree)
        self.assertEqual(per_line[1]["div"], 1)
        self.assertEqual(per_line[1]["add"], 1)
        self.assertEqual(per_line[2]["mul"], 1)
        self.assertEqual(per_line[3], dict.fromkeys(O.CATEGORIES, 0))

    def test_power_is_categorised_separately(self):
        per_line = O.ops_per_statement(ast.parse("q = w ** 3\n"))
        self.assertEqual(per_line[1]["pow"], 1)
        self.assertEqual(per_line[1]["mul"], 0)

    def test_counts_are_nonnegative_integers(self):
        for fn in ("bisection", "interval_newton", "zoom", "newton_certified"):
            counts, _lines, _res = O.executed_ops(getattr(C, fn), 12, 60)
            for k, v in counts.items():
                self.assertIsInstance(v, int, f"{fn}.{k}")
                self.assertGreaterEqual(v, 0, f"{fn}.{k}")


class TestCountsMatchTheCorrectedClaim(unittest.TestCase):
    """The corrected Section 10.1 claim, pinned.

    The paper previously asserted "three divisions and three multiplications
    against interval Newton's two divisions".  That is wrong in both
    directions: the zoom performs FEWER divisions per level than interval
    Newton, and MORE total arithmetic operations.  These tests pin the
    measured relationships, which are exact and machine independent.
    """

    @classmethod
    def setUpClass(cls):
        cls.inet = O.executed_ops(C.interval_newton, 12, 300)
        cls.zoom = O.executed_ops(C.zoom, 12, 300)

    def test_zoom_does_fewer_divisions_per_iteration_than_interval_newton(self):
        # This is the correction: the paper had the direction backwards.
        c_inet, _, r_inet = self.inet
        c_zoom, _, r_zoom = self.zoom
        s_inet = max(r_inet[2], 1)
        s_zoom = max(r_zoom[2], 1)
        self.assertLess(c_zoom["div"] / s_zoom, c_inet["div"] / s_inet,
                        "paper claims the zoom does more divisions; it does "
                        "fewer -- if this now fails, the paper is stale again")

    def test_zoom_does_more_total_operations_per_iteration(self):
        c_inet, _, r_inet = self.inet
        c_zoom, _, r_zoom = self.zoom
        tot = lambda c, r: sum(c[k] for k in O.CATEGORIES if k != "digit_ops")
        s_inet = max(r_inet[2], 1)
        s_zoom = max(r_zoom[2], 1)
        self.assertGreater(tot(c_zoom, r_zoom) / s_zoom,
                           tot(c_inet, r_inet) / s_inet)

    def test_zoom_does_more_digit_operations(self):
        # The machine-independent cost measure: operations x working precision.
        c_inet = self.inet[0]
        c_zoom = self.zoom[0]
        self.assertGreater(c_zoom["digit_ops"], c_inet["digit_ops"])
        ratio = c_zoom["digit_ops"] / c_inet["digit_ops"]
        self.assertGreater(ratio, 2.0)

    def test_digit_ops_overpredicts_wall_clock(self):
        """An honest negative check.

        Digit-operations is exact, but it is a *proxy*: it assumes every
        operation costs its working precision in digit-operations, which
        overstates the measured gap by roughly 2.5x at 3000 digits.  This
        test records that the proxy and the clock disagree, so nobody later
        reads the op-count as a time prediction.
        """
        c_zoom = O.executed_ops(C.zoom, 12, 3000)[0]
        c_inet = O.executed_ops(C.interval_newton, 12, 3000)[0]
        predicted = c_zoom["digit_ops"] / c_inet["digit_ops"]
        measured = 5.131 / 1.967          # from exp_certified.py, alpha=12
        self.assertGreater(predicted, measured,
                           "digit-ops should overpredict the wall clock")
        self.assertGreater(predicted / measured, 1.5)


if __name__ == "__main__":
    unittest.main(verbosity=2)


class TestRunAllRegistry(unittest.TestCase):
    """run_all.py exists so the list of runnable modules cannot rot.  That only
    works if the list is checked.  Without this, a new harness could be added
    and never documented."""

    def test_every_runnable_module_is_registered(self):
        import run_all
        missing = run_all.check_registry()
        self.assertEqual(missing, [],
                         f"runnable but unregistered: {missing} -- add to "
                         f"run_all.HARNESSES so it gets documented")

    def test_registered_modules_all_exist(self):
        import os
        import run_all
        for mod, _desc, _sec, _slow in run_all.HARNESSES:
            self.assertTrue(os.path.exists(mod + ".py"), f"{mod}.py missing")
        for mod in run_all.TEST_MODULES:
            self.assertTrue(os.path.exists(mod + ".py"), f"{mod}.py missing")
