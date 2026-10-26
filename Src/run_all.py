"""One registry of everything runnable in this package, and a single runner.

Why this file exists
--------------------
`Src/` grew to ten runnable entry points: seven experiment harnesses (one per
paper section) and four verification tools.  Every one of them had to be
documented in `Readme.md`, `ToDo.md`, `Src/README.md` and `Paper1.md` -- 35
invocations across four files, all of which had to be kept in sync by hand,
and which had already drifted (an early count said 66 tests when there were 91).

So: describe them once, here, and let `run_all.py` enumerate them.  The docs
can then point at `python3 run_all.py --list` instead of repeating a list that
rots.

What counts as a harness
------------------------
Anything with a `main()` that prints a report -- i.e. anything you would run to
reproduce a number in the paper.  The `test_*.py` modules are excluded: they
run through `unittest`, not individually.

Run:
    python3 run_all.py --list      # what exists, and which paper section
    python3 run_all.py             # the fast checks (default)
    python3 run_all.py --all       # everything, including the slow sweeps
    python3 run_all.py certified   # just the harnesses matching a substring
"""

import argparse
import subprocess
import sys
import time

# (module, what it produces, paper section, slow?)
HARNESSES = [
    ("exp_precision", "Sec 10: naive vs scheduled vs Newton, certified", "10", False),
    ("exp_hybrid", "Sec 10: float-prefix hybrid (a pessimization)", "10", False),
    ("exp_compare", "Sec 10: cubic vs plain vs Newton, cost breakdown", "10", False),
    ("exp_certified", "Sec 10.1: bisection / interval Newton / zoom / "
                      "Newton+res, all exact-integer certified", "10.1", True),
    ("exp_u_aware", "Sec 11: u-aware bound, cross-check (a sample, not the proof)",
     "11", False),
    ("exp_depth", "Sec 11: recursion depth and cost table", "11", False),
    ("exp_mirror", "the mirror hypothesis: reflection across r", "--", False),
    ("verify_theorems", "the proved / decided / measured ledger", "11", False),
    ("opcount", "exact executed-arithmetic counts", "10, 10.1", False),
    ("mutation_test", "14 injected faults; 12 caught, 2 expected survivors",
     "8.1", True),
    ("exp_halley", "certified Halley (order-4 bracket); optimised intervals",
     "--", True),
    ("exp_nth_root", "certified Halley bracket for alpha^(1/n), n=2..8",
     "Paper2 §9", True),
    ("mobius_cf_theorem", "Möbius-CF correspondence for pronic numbers "
     "(Paper3)", "--", False),
    ("mobius_tree", "rational approximant tree from the Möbius closed form",
     "--", False),
]

TEST_MODULES = [
    "test_bracket_zoom", "test_u_aware", "test_exp_mirror", "test_certified",
    "test_exact_decide", "test_opcount", "test_mobius", "test_mobius_n",
    "test_decimal_oracle",
]


def _run(mod, timeout=None):
    t0 = time.monotonic()
    try:
        p = subprocess.run([sys.executable, mod + ".py"],
                           capture_output=True, text=True,
                           timeout=timeout)
        rc, out = p.returncode, p.stdout + p.stderr
    except subprocess.TimeoutExpired:
        rc, out = 124, "<timed out>"
    return rc, out, time.monotonic() - t0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("pattern", nargs="?", default=None,
                    help="only harnesses whose name contains this")
    ap.add_argument("--list", action="store_true", help="list and exit")
    ap.add_argument("--all", action="store_true",
                    help="include the slow sweeps")
    ap.add_argument("--no-tests", action="store_true",
                    help="skip the unittest suite")
    args = ap.parse_args(argv)

    if args.list:
        print(f"{'harness':<18} {'sec':<8} {'slow':<6} produces")
        print("-" * 92)
        for mod, desc, sec, slow in HARNESSES:
            print(f"{mod:<18} {sec:<8} {'yes' if slow else '-':<6} {desc}")
        print(f"\nplus {len(TEST_MODULES)} unittest modules: "
              + " ".join(TEST_MODULES))
        return 0

    failures = []

    if not args.no_tests:
        print("=" * 92)
        print("TESTS")
        print("=" * 92)
        rc, out, dt = _unittest()
        print(out.rstrip().splitlines()[-1] if out else "")
        print(f"  {'ok' if rc == 0 else 'FAILED'} in {dt:.1f}s\n")
        if rc:
            failures.append("tests")

    todo = [h for h in HARNESSES
            if args.all or not h[3]]
    if args.pattern:
        todo = [h for h in todo if args.pattern in h[0]]

    for mod, desc, sec, slow in todo:
        print("=" * 92)
        print(f"{mod}.py  --  {desc}   [Sec {sec}]")
        print("=" * 92)
        rc, out, dt = _run(mod)
        tail = out.strip().splitlines()
        for line in tail[-14:]:
            print("  " + line)
        print(f"  [{'ok' if rc == 0 else 'FAILED rc=' + str(rc)}  {dt:.1f}s]\n")
        if rc:
            failures.append(mod)

    if failures:
        print(f"FAILED: {', '.join(failures)}")
        return 1
    print("everything ran clean")
    return 0


def _unittest():
    t0 = time.monotonic()
    p = subprocess.run([sys.executable, "-m", "unittest", *TEST_MODULES],
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr, time.monotonic() - t0


def check_registry():
    """Every runnable module must be in HARNESSES or TEST_MODULES.

    Without this the list rots -- which is the whole reason this file exists.
    Run as part of the test suite so a new harness cannot be added without
    being registered (and therefore documented).
    """
    import ast
    import os
    registered = {m for m, _, _, _ in HARNESSES} | set(TEST_MODULES)
    registered.add("run_all")
    missing = []
    for fn in sorted(os.listdir(".")):
        if not fn.endswith(".py") or fn.startswith("_"):
            continue
        mod = fn[:-3]
        if mod in registered:
            continue
        tree = ast.parse(open(fn).read())
        has_main = any(
            isinstance(n, ast.If) and "main()" in ast.dump(n.test)
            for n in tree.body)
        if has_main:
            missing.append(mod)
    return missing


if __name__ == "__main__":
    raise SystemExit(main())
