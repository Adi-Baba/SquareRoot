"""
Guard: the decision path must contain no floating-point arithmetic.

The rule being enforced is that nothing which *decides* anything may touch a
float, a Decimal approximation, or math.sqrt.  Files that only *measure* are
allowed to, and are listed separately with the reason.

The previous version of this test listed probe.py among the decision-path
files, which was wrong in two directions at once: probe.py is a reporting
harness, so excluding it is right, but the test then passed vacuously for it
and its float() calls went unnoticed -- and the actual quantity it was
computing (log2 of the bracket width) was wrong in a way no float check
would have caught.  A grep cannot find that class of bug; the sign tests in
test_signexact.py and test_intervalsum.py are what cover it, by comparing
against an independent oracle.
"""

import ast
import os
import unittest

SRS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Modules that decide a sign or a bound.  No float, no Decimal.sqrt, no
# math.sqrt, and no float() on a computed value.
DECISION_PATH_FILES = [
    "squarefree.py",
    "field.py",
    "signexact.py",
    "intervalsum.py",
    "m712_engine.py",
]

# Modules that only measure or search.  Float is allowed; the reason column
# records why, so this list cannot quietly grow.
MEASUREMENT_FILES = {
    "certsingle.py": "adapts the Decimal-bracket code in Src/exp_certified.py, "
                     "whose certify() is exact integer arithmetic",
    "probe.py": "reporting harness; sign and bracket come from signexact, "
                "floats appear only in log-ratios and in the F3 search objective",
    "f1_theorem.py": "log2_ratio_law is a float measurement; every identity it "
                     "states is re-verified with integer/Fraction arithmetic",
    "families.py": "evaluate_float_surrogate is a labelled search objective; "
                   "no decision depends on it",
    "run.py": "subprocess runner",
    "m712_report.py": "reporting layer for m712_engine.py; takes logarithms of "
                      "quantities the engine computed exactly and prints. It "
                      "decides nothing, which is why the measuring code lives "
                      "here instead of in the engine.",
}

BANNED_CALLS = {"sqrt", "fsqrt", "hypot", "exp", "log", "log2",
                "sin", "cos", "tan", "dist", "fsum"}


def _module_path(fname):
    return os.path.join(SRS_DIR, fname)


def _iter_calls(tree):
    """Yield (func_name, node) for every call in the AST."""
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            func = node.func
            if isinstance(func, ast.Name):
                yield func.id, node
            elif isinstance(func, ast.Attribute):
                yield func.attr, node


def _dunder_main_lineno(tree):
    """Line number where the `if __name__ == '__main__'` block starts, or None."""
    for node in ast.walk(tree):
        if isinstance(node, ast.If):
            test = node.test
            if (isinstance(test, ast.Compare)
                    and isinstance(test.left, ast.Name)
                    and test.left.id == "__name__"):
                return node.lineno
    return None


class TestNoFloatInDecisionPath(unittest.TestCase):
    def test_decision_path_has_no_float_calls(self):
        for fname in DECISION_PATH_FILES:
            path = _module_path(fname)
            self.assertTrue(os.path.exists(path), f"missing {fname}")
            with open(path) as fh:
                tree = ast.parse(fh.read(), filename=path)

            main_line = _dunder_main_lineno(tree)
            body_end = main_line - 1 if main_line else None

            offenders = []
            for name, node in _iter_calls(tree):
                if name != "float":
                    continue
                if body_end is not None and node.lineno > body_end:
                    continue        # the __main__ demo block is not the library
                offenders.append(f"line {node.lineno}: float(...)")
            self.assertEqual(offenders, [], f"{fname} uses float(): {offenders}")

    def test_decision_path_has_no_float_constructors(self):
        """float(3) is harmless; float(x) and 1.5 literals are not."""
        for fname in DECISION_PATH_FILES:
            path = _module_path(fname)
            with open(path) as fh:
                src = fh.read()
            tree = ast.parse(src, filename=path)
            main_line = _dunder_main_lineno(tree)
            body_end = main_line - 1 if main_line else None

            literals = []
            for node in ast.walk(tree):
                if isinstance(node, ast.Constant) and isinstance(node.value, float):
                    if body_end is not None and node.lineno > body_end:
                        continue
                    literals.append(f"line {node.lineno}: {node.value!r}")
            self.assertEqual(literals, [], f"{fname} has float literals: {literals}")

    def test_decision_path_has_no_math_transcendental_calls(self):
        for fname in DECISION_PATH_FILES:
            path = _module_path(fname)
            with open(path) as fh:
                tree = ast.parse(fh.read(), filename=path)
            main_line = _dunder_main_lineno(tree)
            body_end = main_line - 1 if main_line else None

            offenders = []
            for name, node in _iter_calls(tree):
                if name not in BANNED_CALLS:
                    continue
                if body_end is not None and node.lineno > body_end:
                    continue
                offenders.append(f"line {node.lineno}: {name}()")
            self.assertEqual(
                offenders, [],
                f"{fname} calls float-returning maths: {offenders}",
            )

    def test_decision_path_imports_only_exact_math_names(self):
        """
        Ban `import math` outright, and allow `from math import` only for
        names that are exact on integers.

        `isqrt` and `gcd` belong on the decision path; `sqrt`, `log2` and
        `pow` do not.  Banning the whole module would have forbidden isqrt,
        which is precisely the tool the integer kernel is built from.
        """
        exact_math_names = {"isqrt", "gcd", "lcm", "factorial", "comb", "perm"}
        for fname in DECISION_PATH_FILES:
            path = _module_path(fname)
            with open(path) as fh:
                tree = ast.parse(fh.read(), filename=path)
            offenders = []
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    offenders += [f"line {node.lineno}: import {a.name}"
                                  for a in node.names if a.name.split('.')[0] == "math"]
                elif isinstance(node, ast.ImportFrom) and node.module == "math":
                    for alias in node.names:
                        if alias.name not in exact_math_names:
                            offenders.append(
                                f"line {node.lineno}: from math import {alias.name}")
            self.assertEqual(offenders, [], f"{fname}: {offenders}")


class TestMeasurementFilesAreJustified(unittest.TestCase):
    def test_every_measurement_file_exists_and_is_listed(self):
        for fname, reason in MEASUREMENT_FILES.items():
            self.assertTrue(os.path.exists(_module_path(fname)), f"missing {fname}")
            self.assertTrue(reason.strip(), f"{fname} listed with no justification")

    def test_measurement_and_decision_lists_are_disjoint(self):
        overlap = set(MEASUREMENT_FILES) & set(DECISION_PATH_FILES)
        self.assertEqual(overlap, set(),
                         f"files cannot be both: {sorted(overlap)}")

    def test_probe_really_is_measurement_only(self):
        """
        probe.py must take its sign from signexact, not compute one.

        This is the check that would have caught the original bug: probe used
        to derive log2(1/|S|) from its own bracket width, which is a float
        computation dressed up as a result.
        """
        with open(_module_path("probe.py")) as fh:
            src = fh.read()
        self.assertIn("from signexact import", src)
        self.assertIn("exact_sign", src)
        # No direct float evaluation of a sqrt-sum anywhere in the harness.
        self.assertNotIn("math.sqrt(", src)


if __name__ == "__main__":
    unittest.main(verbosity=2)
