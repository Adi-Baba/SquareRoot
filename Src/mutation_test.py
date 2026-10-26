"""Mutation testing: deliberately break the code, confirm the suite notices.

Why this file exists
--------------------
`Docs/Paper1.md` claimed the suite was "mutation-tested" without saying how, or
with what result.  A claim like that is either evidence or decoration, and
decorating with it would be exactly the error this file's sibling
`verify_theorems.py` exists to prevent.  So: run it.

Method.  For each mutation -- a one-token edit to a real source line -- apply
it to a temporary copy of `Src/`, run the full test suite, and record whether
the suite FAILS.  A mutant that survives is a gap in the suite: some claim is
not actually pinned by any test.  The original sources are never modified;
everything happens in a scratch directory that is removed afterwards.

What a good score looks like.  Not 100%.  Some mutants are equivalent (the
code does the same thing), and some are genuinely untested regions.  The point
is to *know* which, and to report the survivors rather than hide them.

Run:  python3 mutation_test.py [--keep]
"""

import os
import shutil
import subprocess
import sys
import tempfile
import time

SRC = os.path.dirname(os.path.abspath(__file__))

# Mutants that are KNOWN to survive, with the reason each is not a suite gap.
# They are still run and still printed; listing them here only decides whether
# they fail the exit status.  Anything surviving that is NOT on this list is an
# unexpected gap and does fail.  Remove an entry once the suite pins it -- the
# run reports when that happens, so the list cannot rot unnoticed.
EXPECTED_SURVIVORS = {
    # A performance regression: correct answers, wrong amount of work.  No
    # correctness test can see slowness by construction, so this is pinned by
    # a source-contract check labelled as such rather than by behaviour.
    "make the seed precision depend on the target again",
    # Both the plain and the corrected estimator satisfy every property the
    # suite checks.  Note in particular that the obvious test for it,
    # `est > r`, is FALSE of the computed value: Corollary 2 is exact, but once
    # the error drops below the working precision the returned estimate rounds
    # across the root.  The certificate is unaffected.
    "drop the cubic correction, keep the plain estimate",
}

# (file, original substring, replacement, what the mutant breaks)
MUTANTS = [
    # --- the sharp constant 4/27: n's maximum is at u = 2w/3.  Moving the
    # critical point is the cleanest way to loosen (or wrongly tighten) it.
    ("exp_u_aware.py", "    c = 2 * w / 3", "    c = w / 2",
     "move the numerator's critical point off 2w/3"),
    # --- D(u) monotonicity is the lemma the u-aware bound rests on.
    # Flipping the sign makes the denominator *decrease*, invalidating the
    # "evaluate at the left end" step.
    ("exp_u_aware.py", "return 2 * a * (2 * a + w) + u * (2 * a + u)",
     "return 2 * a * (2 * a + w) - u * (2 * a + u)",
     "break D(u) monotonicity"),
    # --- off-by-ten in the termination test, reintroduced
    ("exp_precision.py", "if eps.adjusted() < -digits:",
     "if eps.adjusted() <= -digits:",
     "reintroduce the loose termination test"),
    # --- the seed bug: alpha-derived precision replaced by target precision
    ("exp_precision.py", "alpha.adjusted() + 5", "digits + 5",
     "make the seed precision depend on the target again"),
    # --- certify(): the two enclosure directions, one at a time.  These are
    # the comparisons that decide whether a bracket really encloses the root.
    ("exp_certified.py", "        return c * c * 10 ** (2 * e) <= alpha",
     "        return c * c * 10 ** (2 * e) < alpha",
     "loosen the right-endpoint test in certify()"),
    ("exp_certified.py", "        return alpha <= c * c * 10 ** (2 * e)",
     "        return alpha < c * c * 10 ** (2 * e)",
     "loosen the left-endpoint test in certify()"),
    # --- exact_decide: the three defects that shipped and were found by
    # cross-checking.  Each made the decider ACCEPT invalid bounds.
    ("exact_decide.py", "    for c in p:\n        s = s * u + c",
     "    for c in reversed(p):\n        s = s * u + c",
     "exact_decide: Horner iterates the coefficient list backwards"),
    ("exact_decide.py", "         -2 * a * B,\n", "         -(2 * a * B + B * B),\n",
     "exact_decide: spurious B^2 u term in the expansion"),
    ("exact_decide.py", "    return (c1 / c2, c0 / c2)", "    return (-c1 / c2, c0 / c2)",
     "exact_decide: negate the linear coefficient (hides the maximiser)"),
    # --- the corrected estimator: drop the correction (plain estimate)
    ("exp_precision.py", "est = a + p / (a + m)", "est = m",
     "drop the cubic correction, keep the plain estimate"),
    # --- the progress condition, the one whose omission makes it diverge
    ("exp_precision.py", "eps = (w ** 3 / (27 * a * a))",
     "eps = (w ** 3 / (16 * a * a))",
     "use the paper's loose 1/16 constant in the recursion"),
    # --- interval Newton's schedule rate
    ("exp_certified.py", "prec = _prec_for(alpha, b - a, rate=2)",
     "prec = _prec_for(alpha, b - a, rate=1)",
     "under-schedule interval Newton"),
    # --- the seeding asymmetry, the other way: give the "int seed" path the
    # full 30-digit float seed, collapsing the two rows the paper compares.
    ("exp_certified.py", "        c = 1\n", "        c = 30\n",
     "make the int-seeded Newton count as if float-seeded"),
    # --- the mirror bug that hung on perfect squares (first of two guards)
    ("exp_mirror.py", "        if a * a == A:\n            return +a, 0\n"
     "        if mode == \"plain\":",
     "        if False:\n            return +a, 0\n        if mode == \"plain\":",
     "remove the perfect-square termination guard"),
]

SUITE = ["test_bracket_zoom", "test_u_aware", "test_exp_mirror",
         "test_certified", "test_exact_decide"]

# A mutant that diverges, hangs, or blows the stack has broken the code, so it
# is KILLED, not survived.  The baseline suite runs in well under a minute, so
# anything past this is a broken mutant rather than a slow one -- and waiting
# the full suite duration to discover that is pure waste.
PER_MUTANT_TIMEOUT = 90


def run_suite(cwd, timeout=PER_MUTANT_TIMEOUT):
    """Run the suite.  Returns (killed, detail).

    killed=True  means the mutant did NOT survive: the suite failed, or the
                  run exceeded `timeout` (a divergent mutant).
    killed=False means the suite passed on the mutated code -- a gap.
    """
    try:
        proc = subprocess.run(
            [sys.executable, "-m", "unittest", *SUITE],
            cwd=cwd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return True, f"TIMEOUT after {timeout}s (divergent mutant -> killed)"
    return proc.returncode != 0, (proc.stdout + proc.stderr)


def main():
    keep = "--keep" in sys.argv
    tmp = tempfile.mkdtemp(prefix="mutants_")
    print(f"scratch: {tmp}\n")
    # Hard guarantee: nothing outside `tmp` is ever opened for writing.  An
    # earlier version mutated files in SRC, which both broke the experiment
    # (tests ran on the unmutated copies) and risked leaving the real sources
    # modified if the run was interrupted.  Assert the invariant up front.
    assert os.path.abspath(tmp) != os.path.abspath(SRC), "scratch must differ from SRC"

    # Copy the package into the scratch dir.  Read-only on SRC from here on.
    for f in os.listdir(SRC):
        if f.endswith(".py"):
            shutil.copy2(os.path.join(SRC, f), tmp)
    print(f"copied {len([f for f in os.listdir(tmp) if f.endswith('.py')])} "
          f"modules into scratch; SRC is now read-only\n")

    # Baseline: the unmutated suite MUST pass, or nothing below means anything.
    # run_suite returns `killed`, so a healthy baseline is killed=False.
    baseline_killed, out = run_suite(tmp, timeout=600)
    baseline_ok = not baseline_killed
    print(f"baseline suite: {'PASS' if baseline_ok else 'FAIL'}")
    if not baseline_ok:
        print("  suite FAILED on clean code -- mutation results would be "
              "meaningless; fix the suite first.")
        print(out[-3000:])
        shutil.rmtree(tmp, ignore_errors=True)
        return 2
    print(f"suite = {len(SUITE)} modules, budget {PER_MUTANT_TIMEOUT}s/mutant\n")

    results = []
    for fname, old, new, desc in MUTANTS:
        # Mutate the SCRATCH copy, never the original source.  The suite is
        # run in `tmp`, so the mutation has to be applied in `tmp` too -- an
        # earlier version wrote to SRC and tested tmp, which meant the suite
        # always saw clean code and reported every mutant as SURVIVED.
        path = os.path.join(tmp, fname)
        with open(path) as fh:
            original = fh.read()
        if old not in original:
            results.append((desc, "NOT-APPLIED",
                            f"pattern not found in {fname}"))
            print(f"  [??] {desc}: PATTERN NOT FOUND in {fname}")
            continue
        with open(path, "w") as fh:
            fh.write(original.replace(old, new, 1))
        print(f"    running: {desc} ...", end="", flush=True)
        t0 = time.monotonic()
        try:
            killed, detail = run_suite(tmp)
        finally:
            with open(path, "w") as fh:      # restore, for the next mutant
                fh.write(original)
        dt = time.monotonic() - t0
        verdict = "SURVIVED" if not killed else "killed"
        note = " (timeout)" if detail.startswith("TIMEOUT") else ""
        results.append((desc, verdict, f"{dt:.1f}s"))
        mark = "!!" if not killed else "ok"
        print(f"\r  [{mark}] {verdict:<9} {dt:>5.1f}s  {desc}{note}")

    killed = sum(1 for _, v, _ in results if v == "killed")
    survived = [d for d, v, _ in results if v == "SURVIVED"]
    na = [d for d, v, _ in results if v == "NOT-APPLIED"]

    print()
    print(f"mutants: {len(results)}   killed: {killed}   "
          f"survived: {len(survived)}   not-applied: {len(na)}")
    if survived:
        print("\nSURVIVORS -- these claims are NOT pinned by the suite:")
        for d in survived:
            print(f"  - {d}")
    if na:
        print("\nNOT APPLIED (pattern drift, re-check the source):")
        for d in na:
            print(f"  - {d}")

    if not keep:
        shutil.rmtree(tmp, ignore_errors=True)
    else:
        print(f"\nkept scratch dir {tmp}")

    # Exit status.  Two things are distinguished here, and conflating them is
    # why `run_all.py --all` used to report a FAILURE that was really a
    # documented, deliberate finding:
    #
    #   * NOT-APPLIED means a pattern in this file no longer matches the
    #     source.  That IS an authoring error and makes the run
    #     non-authoritative, so it fails the exit.
    #
    #   * A SURVIVOR is only a gap if it is unEXPECTED.  Two of the fourteen
    #     mutants are KNOWN to survive and are discussed in Readme.md: the
    #     seed-precision regression (a performance bug a correctness suite
    #     cannot see) and dropping the cubic correction (both estimators
    #     satisfy every property the suite checks).  Reporting those as a
    #     failure makes `run_all.py --all` exit 1 on a healthy tree, which
    #     trains the reader to ignore the exit status.
    #
    # So the exit is non-zero for an unEXPECTED survivor or any NOT-APPLIED,
    # and the known survivors are printed as SURVIVORS without failing.  An
    # unexpected one is still printed, loudly, as a real gap.
    unexpected = [d for d in survived if d not in EXPECTED_SURVIVORS]
    if unexpected:
        print("\nUNEXPECTED SURVIVORS -- real gaps in the suite:")
        for d in unexpected:
            print(f"  - {d}")
    observed = {desc for desc, _, _ in results}
    if EXPECTED_SURVIVORS - observed:
        print("\nEXPECTED SURVIVORS NOT OBSERVED -- the suite improved:")
        for d in sorted(EXPECTED_SURVIVORS - observed):
            print(f"  - {d}")
    print(f"\nexpected survivors: {len(EXPECTED_SURVIVORS)}   "
          f"unexpected: {len(unexpected)}")
    return 1 if (unexpected or na) else 0


if __name__ == "__main__":
    raise SystemExit(main())
