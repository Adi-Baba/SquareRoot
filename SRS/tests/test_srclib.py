"""
Tests for _srclib — loading sibling Src/ modules by explicit path.

The bugs these guard against were not crashes. Both were cases where the
code imported successfully and produced a subtly wrong result:

  * `sys.path.insert(...)` at import time, which a static analyzer cannot
    follow. It reported `Import "exp_certified" could not be resolved` for
    a line that worked at run time, and it mutated global interpreter state
    as a side effect of importing a module.

  * Loading a module without registering it in `sys.modules`, so a second
    caller got a second, independently executed copy carrying its own
    Decimal context state. Two copies of exp_certified in one process is a
    correctness hazard, not a tidiness complaint.

These run in a subprocess where isolation matters, because the point is
about what lands in `sys.modules` and `sys.path` and that is global state.
"""

import os
import subprocess
import sys
import unittest

SRS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(SRS_DIR)


def run_python(code):
    """Run a snippet in a fresh interpreter with SRS/ importable."""
    env = dict(os.environ)
    env["PYTHONPATH"] = SRS_DIR + os.pathsep + env.get("PYTHONPATH", "")
    proc = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True, text=True, cwd=REPO_ROOT, env=env, timeout=120,
    )
    if proc.returncode != 0:
        raise AssertionError(f"snippet failed:\n{proc.stdout}\n{proc.stderr}")
    return proc.stdout.strip()


class TestSrcPath(unittest.TestCase):
    def test_src_path_points_at_real_files(self):
        import _srclib
        for name in ("exp_certified", "exp_precision", "exact_decide"):
            path = _srclib.src_path(name)
            self.assertTrue(os.path.isfile(path), f"{path} does not exist")

    def test_src_dir_is_the_real_sibling_not_srs_src(self):
        """
        The earlier conftest bug: one level too shallow resolved ../Src to
        SRS/Src, which does not exist.
        """
        import _srclib
        self.assertTrue(os.path.isdir(_srclib._SRC_DIR), _srclib._SRC_DIR)
        self.assertEqual(os.path.basename(_srclib._SRC_DIR), "Src")
        self.assertEqual(os.path.basename(os.path.dirname(_srclib._SRC_DIR)), "Irrational")


class TestLoadSrcModule(unittest.TestCase):
    def test_missing_module_names_the_path_it_tried(self):
        import _srclib
        with self.assertRaises(ImportError) as ctx:
            _srclib.load_src_module("definitely_not_a_module_xyz")
        message = str(ctx.exception)
        self.assertIn("definitely_not_a_module_xyz", message)
        # The point of the message is to say which file was looked for.
        self.assertIn("_SRC_DIR", message)

    def test_loads_and_exposes_the_expected_functions(self):
        import _srclib
        mod = _srclib.load_src_module("exp_certified")
        for attr in ("certify", "zoom", "bisection", "interval_newton",
                     "newton_certified"):
            self.assertTrue(callable(getattr(mod, attr)), f"{attr} missing")
        # The old alias was `bisect`, which does not exist in exp_certified.
        self.assertFalse(hasattr(mod, "bisect"))

    def test_bare_module_name_is_required(self):
        """
        `module_name` is interpolated into a path, so a value containing
        separators would resolve outside `Src/` and execute that file.
        """
        import _srclib
        for bad in ("../Src/exp_certified", "../../etc/passwd", "a/b", ""):
            with self.subTest(bad=bad):
                with self.assertRaises(ImportError, msg=f"{bad!r} was accepted"):
                    _srclib.load_src_module(bad)

    def test_docstring_agrees_with_the_sys_modules_behaviour(self):
        """
        The docstring used to say the module was NOT registered in
        `sys.modules` under its plain name, which is the opposite of what the
        code does, and `_deps.py` said the correct thing.  Assert the
        documented claim matches the observed behaviour so the two files
        cannot disagree again.
        """
        out = run_python(
            "import sys, _srclib;"
            "m = _srclib.load_src_module('exp_certified');"
            "import exp_certified;"
            "print('exp_certified' in sys.modules, exp_certified is m)"
        )
        self.assertEqual(out.split(), ["True", "True"],
                         "the code must register under the plain name and "
                         "hand the same object to a later bare import")
        import _srclib
        doc = _srclib.__doc__ or ""
        self.assertNotIn(
            "rather than being\nregistered in `sys.modules`", doc,
            "docstring still claims the module is not registered in sys.modules"
        )
        self.assertNotIn("will NOT pick up this instance", doc)


class TestNoGlobalStateDamage(unittest.TestCase):
    def test_src_not_left_on_sys_path(self):
        out = run_python(
            "import sys; sys.path.insert(0, %r);"
            "import _srclib; _srclib.load_src_module('exp_certified');"
            "print(any(p.endswith('/Src') for p in sys.path))" % SRS_DIR
        )
        self.assertEqual(out, "False", "load_src_module leaked Src onto sys.path")

    def test_single_module_object_across_import_orders(self):
        """Both orders must yield one shared object, not two copies."""
        for first, second in (("_deps", "certsingle"), ("certsingle", "_deps")):
            with self.subTest(first=first):
                out = run_python(
                    f"import {first}, {second}, exp_certified;"
                    "print(certsingle.EC is exp_certified,"
                    "      _deps._cert is exp_certified)"
                )
                self.assertEqual(out.split(), ["True", "True"],
                                 f"duplicate copy via {first}")

    def test_repeat_load_returns_identical_object(self):
        out = run_python(
            "import sys; sys.path.insert(0, %r); import _srclib;"
            "a = _srclib.load_src_module('exp_certified');"
            "b = _srclib.load_src_module('exp_certified');"
            "print(a is b)" % SRS_DIR
        )
        self.assertEqual(out, "True", "second load built a duplicate module")

    def test_import_order_independent_for_deps_exports(self):
        out = run_python(
            "import sys; sys.path.insert(0, %r); import _deps;"
            "print(callable(_deps.bisection), callable(_deps.zoom),"
            "      callable(_deps.decide_u_aware_bound))" % SRS_DIR
        )
        self.assertEqual(out, "True True True")


class TestDepsAndCertsingle(unittest.TestCase):
    def test_certsingle_does_not_mutate_sys_path_on_import(self):
        out = run_python(
            "import sys; sys.path.insert(0, %r);"
            "before = list(sys.path); import certsingle;"
            "print(before == sys.path)" % SRS_DIR
        )
        self.assertEqual(out, "True", "importing certsingle changed sys.path")

    def test_no_bare_sys_path_insert_remains_in_srs(self):
        """
        A structural check, so the hack cannot come back unnoticed.

        `_srclib` and `conftest` are the only places allowed to touch
        `sys.path`; both do so in a scoped, documented way.
        """
        allowed = {"_srclib.py"}
        offenders = []
        for name in sorted(os.listdir(SRS_DIR)):
            if not name.endswith(".py") or name in allowed:
                continue
            with open(os.path.join(SRS_DIR, name)) as fh:
                for lineno, line in enumerate(fh, 1):
                    stripped = line.strip()
                    if stripped.startswith("#") or stripped.startswith("*"):
                        continue
                    if "sys.path.insert" in stripped or "sys.path.append" in stripped:
                        offenders.append(f"{name}:{lineno}: {stripped}")
        self.assertEqual(offenders, [], f"sys.path mutation outside _srclib: {offenders}")

    def test_deps_exports_are_callable(self):
        import _deps
        for name in ("certify", "zoom", "bisection", "interval_newton",
                     "newton_certified", "decide_u_aware_bound",
                     "decide_published_bound", "decide_error_identity"):
            self.assertTrue(callable(getattr(_deps, name)), f"_deps.{name} not callable")

    def test_bisection_alias_points_at_the_real_function(self):
        """Regression: the old alias was `bisect`, which never existed."""
        import _deps
        import certsingle
        self.assertIs(_deps.bisection, certsingle.bisection)


if __name__ == "__main__":
    unittest.main(verbosity=2)
