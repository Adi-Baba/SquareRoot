#!/usr/bin/env python3
"""
SRS Runner - single entry point for all SRS experiments.

Usage:
    python3 run.py --list
    python3 run.py                    # tests + fast harnesses
    python3 run.py --all              # everything including slow sweeps
    python3 run.py probe              # just the probe
    python3 run.py probe --family F1  # just F1 family
"""

import argparse
import sys
import os
import subprocess
import time

SRS_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(os.path.dirname(SRS_DIR), 'Src')

# Registry of all runnable modules in SRS
HARNESSES = [
    ("probe", "SRS probe: family x size experiments"),
    ("certsingle", "Certified single-root brackets"),
    ("signexact", "Exact sign determination"),
    ("intervalsum", "Certified interval arithmetic"),
    ("families", "Instance generators"),
    ("f1_theorem", "F1 theorem verification"),
    ("squarefree", "Squarefree decomposition"),
    ("field", "Multiquadratic field"),
    ("refutation", "Refutation hunt for sup E/N = (r-1)/r (certified)"),
]

# Also include the original Src harnesses
SRC_HARNESSES = [
    ("exp_precision", "Sec 10: naive vs scheduled vs Newton"),
    ("exp_hybrid", "Sec 10: float-prefix hybrid"),
    ("exp_compare", "Sec 10: cubic vs plain vs Newton"),
    ("exp_certified", "Sec 10.1: certified benchmark"),
    ("exp_u_aware", "Sec 11: u-aware bound"),
    ("exp_depth", "Sec 11: recursion depth"),
    ("exp_mirror", "Mirror hypothesis"),
    ("verify_theorems", "Proved/decided/measured ledger"),
    ("opcount", "Exact executed-arithmetic counts"),
    ("mutation_test", "14 injected faults"),
    ("mobius_cf_theorem", "Möbius-CF correspondence for pronic numbers"),
    ("exp_halley", "Certified Halley (order-4 bracket)"),
    ("exp_nth_root", "Certified Halley bracket for alpha^(1/n)"),
    ("mobius_tree", "Rational approximant tree from Möbius closed form"),
]


def _run_module(module: str, args: list = None, cwd: str = None, timeout: int = None) -> tuple:
    """Run a Python module and return (rc, stdout, stderr, elapsed)."""
    if cwd is None:
        cwd = SRS_DIR
    cmd = [sys.executable, '-m', module] + (args or [])
    t0 = time.monotonic()
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
        rc = p.returncode
        out = p.stdout
        err = p.stderr
    except subprocess.TimeoutExpired as e:
        rc = 124
        out = e.stdout or ''
        err = str(e)
    elapsed = time.monotonic() - t0
    return rc, out, err, elapsed


def _run_src_module(module: str, args: list = None, timeout: int = None) -> tuple:
    """Run a module from the Src directory."""
    return _run_module(module, args, cwd=SRC_DIR, timeout=timeout)


def _run_srs_module(module: str, args: list = None, timeout: int = None) -> tuple:
    """Run a module from the SRS directory."""
    return _run_module(module, args, cwd=SRS_DIR, timeout=timeout)


def _list_harnesses():
    """Print all available harnesses."""
    print(f"{'harness':<18} {'description':<60}")
    print("-" * 80)
    for name, desc in HARNESSES:
        print(f"{name:<18} {desc}")
    print("\nSrc harnesses (run from Src/):")
    for name, desc in SRC_HARNESSES:
        print(f"  {name:<18} {desc}")


def _run_tests(verbose: bool = False) -> tuple:
    """Run the SRS test suite via unittest discovery from SRS_DIR."""
    extra = ['-v'] if verbose else []
    return _run_srs_module('unittest', ['discover', '-s', 'tests'] + extra)


def _run_src_tests(verbose: bool = False) -> tuple:
    """Run the original Src test suite."""
    args = ['-m', 'unittest', 'test_bracket_zoom', 'test_u_aware', 'test_exp_mirror',
            'test_certified', 'test_exact_decide', 'test_opcount', 'test_mobius',
            'test_mobius_n', 'test_decimal_oracle']
    if verbose:
        args.append('-v')
    return _run_src_module('unittest', args)


def _run_probe(args: list = None) -> tuple:
    """
    Run the SRS probe.

    Kept as a named wrapper for callers outside this module; `main` invokes
    `_run_srs_module('probe', ...)` directly on its own branch, so this was
    unreferenced internally.
    """
    return _run_srs_module('probe', args or [])


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--list', action='store_true', help='list all harnesses and exit')
    parser.add_argument('--all', action='store_true',
                        help='run everything including slow sweeps. All 14 Src '
                             'harnesses run, so budget ~15-20 min: exp_u_aware '
                             'alone is ~5.5 min and exp_depth ~45 s on this '
                             'machine, and the rest are seconds. This flag used '
                             'to run only 2 of the Src harnesses while still printing '
                             '"SRS FULL SUITE".')
    parser.add_argument('--no-tests', action='store_true', help='skip the test suite')
    parser.add_argument('command', nargs='?', help='command to run (probe, certsingle, etc.)')
    parser.add_argument('args', nargs=argparse.REMAINDER, help='arguments passed to the command')

    args = parser.parse_args()

    if args.list:
        _list_harnesses()
        return 0

    failures = []

    if not args.no_tests:
        print("=" * 80)
        print("SRS TESTS")
        print("=" * 80)
        rc, out, err, dt = _run_tests(verbose=True)
        print(out)
        if err:
            print("STDERR:", err)
        print(f"  {'ok' if rc == 0 else 'FAILED rc=' + str(rc)}  {dt:.1f}s\n")
        if rc != 0:
            failures.append('tests')

    if args.command:
        cmd = args.command
        cmd_args = args.args

        if cmd == 'probe':
            print("=" * 80)
            print("SRS PROBE")
            print("=" * 80)
            rc, out, err, dt = _run_srs_module('probe', args.args)
        elif cmd in dict(HARNESSES):
            print("=" * 80)
            print(f"SRS {cmd.upper()}")
            print("=" * 80)
            rc, out, err, dt = _run_srs_module(cmd, cmd_args)
        elif cmd in dict(SRC_HARNESSES):
            print("=" * 80)
            print(f"SRC {cmd.upper()}")
            print("=" * 80)
            rc, out, err, dt = _run_src_module(cmd, cmd_args)
        else:
            print(f"Unknown command: {cmd}")
            print(f"Available: {[h[0] for h in HARNESSES] + [h[0] for h in SRC_HARNESSES]}")
            return 1

        print(out)
        if err:
            print("STDERR:", err)
        print(f"  {'ok' if rc == 0 else 'FAILED rc=' + str(rc)}  {dt:.1f}s\n")
        if rc != 0:
            failures.append(cmd)

    elif args.all:
        # Run everything
        print("=" * 80)
        print("SRS FULL SUITE (--all)")
        print("=" * 80)

        # SRS tests
        rc, out, err, dt = _run_tests(verbose=True)
        print(out)
        if rc != 0:
            failures.append('tests')

        # Src tests
        rc, out, err, dt = _run_src_tests(verbose=True)
        print(out)
        if rc != 0:
            failures.append('src_tests')

        # All SRS harnesses
        for name, _ in HARNESSES:
            print("=" * 80)
            print(f"SRS {name.upper()}")
            print("=" * 80)
            rc, out, err, dt = _run_srs_module(name)
            print(out)
            if err:
                print("STDERR:", err)
            print(f"  {'ok' if rc == 0 else 'FAILED rc=' + str(rc)}  {dt:.1f}s\n")
            if rc != 0:
                failures.append(name)

        # ALL Src harnesses, not a two-name subset.
        #
        # This used to run only `exp_certified` and `mutation_test` while
        # `--list` advertised all eleven and the banner above claimed a full
        # suite.  A reader had no way to tell that nine of them had been
        # silently skipped.  If any of these is too slow to belong in `--all`,
        # the fix is to say so by name here, not to drop it quietly.
        #
        # Two of them are genuinely slow -- `exp_u_aware` (~5.5 min) and
        # `exp_depth` (~45 s) -- which is almost certainly why they were
        # dropped.  They run here, and the `--all` help text says what to
        # expect, so the cost is visible instead of the coverage being not.
        for name, _ in SRC_HARNESSES:
            print("=" * 80)
            print(f"SRC {name.upper()}")
            print("=" * 80)
            rc, out, err, dt = _run_src_module(name)
            print(out)
            if err:
                print("STDERR:", err)
            print(f"  {'ok' if rc == 0 else 'FAILED rc=' + str(rc)}  {dt:.1f}s\n")
            if rc != 0:
                failures.append(f'src:{name}')

    if failures:
        print(f"FAILED: {', '.join(failures)}")
        return 1
    print("everything ran clean")
    return 0


if __name__ == "__main__":
    sys.exit(main())