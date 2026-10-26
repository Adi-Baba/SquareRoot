"""
Probe harness for the Square-Root Sum problem.

Runs family x size experiments, decides each instance exactly, brackets it
with a certified rational interval, and tabulates the measured separation.

THE BUG THIS FILE USED TO HAVE
------------------------------
`log2_inv_S` was computed as log2(1/(hi - lo)) -- the width of the
*certification bracket*.  That is not a property of S at all: it is a
property of how many digits the harness happened to ask for, so it grew
whenever the bracket was refined and reported a "separation ratio" that was
an artefact of the harness settings.  The width of a bracket is upper-bound
information about |S|; it is not log2(1/|S|).

What is reported now: the exact sign, and a certified bracket [lo, hi] for
S itself.  From that bracket, log2(1/|S|) is bracketed too -- once the sign
is known, |S| lies between min(|lo|,|hi|) and max(|lo|,|hi|), so
log2(1/|S|) lies between -log2(max) and -log2(min).  The table prints that
interval, and a separate column carries the certificate width, clearly
labelled as a property of the certificate.

Everything in the decision path is exact.  Float appears only in
log-ratios for reporting, and the surrogate used to *search* for instances
(F3) is labelled approximate.
"""

import argparse
import csv
import math
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional, Tuple

# Plain sibling imports.  This file used to prepend its own directory to
# sys.path before importing them, which mutated global interpreter state and
# forced the imports below an E402.  SRS/ is on sys.path whenever this module
# is legitimately importable (tests/conftest.py arranges it, and running the
# file directly puts its directory there already), so the mutation bought
# nothing.
from families import (
    FAMILY_SIGNATURES,
    generate_f1,
    generate_f2,
    generate_f3,
    generate_f4,
    generate_f5,
    generate_family_sequence,
)
from intervalsum import CertifiedInterval, certified_interval_for_S
from signexact import SignUnknownError, exact_sign

__all__ = ["ProbeResult", "ProbeRunner", "main"]


@dataclass
class ProbeResult:
    """One probed instance.  All sign/interval fields are exact."""

    family: str
    instance_id: str
    k: int
    N_bits: int

    # exact decision
    exact_sign: int
    exact_method: str
    exact_bits: Optional[int]          # fractional bits the decision needed

    # certified bracket for S (decimal strings; they are long)
    certified_lo: str
    certified_hi: str
    #: True when the bracket determines the sign, i.e. excludes zero.
    #: For an exact zero the bracket is [0, 0], which contains zero, so this
    #: is False there -- the sign is known to be 0, not unknown, which is why
    #: `exact_sign` is the field to read for that case.
    bracket_clears_zero: bool
    digits_pinned: Optional[int]

    # measured (float), reporting only
    log2_inv_S_lo: Optional[float]     # -log2(max(|lo|,|hi|)), a lower bound
    log2_inv_S_hi: Optional[float]     # -log2(min(|lo|,|hi|)), an upper bound
    ratio_lo: Optional[float]          # log2_inv_S_lo / N
    ratio_hi: Optional[float]
    certificate_width_log2: Optional[float]   # property of the certificate
    #: Bit-length *estimate* of log2(1/width), not an exact bit count.  It is
    #: accurate to about 1 bit, and reading it as a precise cost is wrong.
    certificate_width: Optional[int]

    metadata: dict = field(default_factory=dict)


class ProbeRunner:
    def __init__(self, seed: int = 42, digits: int = 50,
                 max_bits: Optional[int] = None):
        self.seed = seed
        self.digits = digits
        #: Precision cap handed to both deciders.  `None` means each decider's
        #: own default.  This exists because `--max-bits` used to be parsed,
        #: documented, and then dropped on the floor: the cap was always
        #: `signexact.DEFAULT_MAX_BITS`, so the flag had no effect at all.
        self.max_bits = max_bits
        self.results: List[ProbeResult] = []

    # -- core ---------------------------------------------------------------

    def probe_instance(self, inst: dict) -> ProbeResult:
        """Decide one instance exactly and bracket it."""
        terms = inst["terms"]
        family = inst["family"]
        n_bits = inst["N_bits"]
        meta = inst.get("metadata", {})

        instance_id = self._instance_id(family, meta)

        sign_result = exact_sign(terms, max_bits=self.max_bits)
        iv = certified_interval_for_S(terms, self.digits, max_bits=self.max_bits)

        # The two must agree; if they ever do not, the harness is broken and
        # that is worth a loud failure rather than a row in a table.
        if sign_result.sign != iv.sign:
            raise AssertionError(
                f"sign disagreement on {terms}: "
                f"exact_sign={sign_result.sign} ({sign_result.method}) "
                f"vs interval={iv.sign} ({iv.bits} bits)"
            )

        lo_log2 = hi_log2 = None
        ratio_lo = ratio_hi = None
        if not iv.contains_zero():
            abs_lo = min(abs(iv.lo), abs(iv.hi))
            abs_hi = max(abs(iv.lo), abs(iv.hi))
            if abs_hi > 0:
                lo_log2 = -math.log2(abs_hi)
            if abs_lo > 0:
                hi_log2 = -math.log2(abs_lo)
            if n_bits > 0:
                ratio_lo = lo_log2 / n_bits if lo_log2 is not None else None
                ratio_hi = hi_log2 / n_bits if hi_log2 is not None else None

        width = iv.hi - iv.lo
        if width > 0:
            # `width` is a Fraction, so it is finite by construction and
            # converting it is total: there is no inf case to guard against,
            # and the earlier `width < math.inf` test could never be false.
            cert_log2 = -math.log2(float(width))
            # Report 2^(-cert_bits): a positive count of how narrow the
            # certificate is.  A signed numerator/denominator bit difference
            # would read as a negative "width", which is nonsense.  This is a
            # bit-length *estimate* of log2(1/width), not an exact bit count;
            # it is off by at most about 1 for a Fraction in lowest terms.
            cert_bits = width.denominator.bit_length() - width.numerator.bit_length()
        else:
            cert_log2 = None
            cert_bits = 0

        return ProbeResult(
            family=family,
            instance_id=instance_id,
            k=inst["k"],
            N_bits=n_bits,
            exact_sign=iv.sign,
            exact_method=sign_result.method,
            exact_bits=sign_result.witness.get("bits"),
            certified_lo=str(iv.lo),
            certified_hi=str(iv.hi),
            bracket_clears_zero=not iv.contains_zero(),
            digits_pinned=self._digits_pinned(iv),
            log2_inv_S_lo=lo_log2,
            log2_inv_S_hi=hi_log2,
            ratio_lo=ratio_lo,
            ratio_hi=ratio_hi,
            certificate_width_log2=cert_log2,
            certificate_width=cert_bits,
            metadata=dict(meta),
        )

    @staticmethod
    def _digits_pinned(iv: CertifiedInterval) -> Optional[int]:
        from intervalsum import digit_count
        return digit_count(iv)

    @staticmethod
    def _instance_id(family: str, meta: dict) -> str:
        # F2 and F3 both describe themselves by (k, seed), and the default
        # random suite runs both over the same grid, so an id built from those
        # two fields alone collides on every default run.  The family is
        # included so an id is unique within a table without relying on the
        # reader also cross-referencing the family column.
        if family == "F1":
            return f"a{meta.get('a')}_n{meta.get('n')}"
        if family in ("F2", "F3"):
            return f"{family}_k{meta.get('k', '?')}_s{meta.get('seed')}"
        if family == "F4":
            return f"a{meta.get('a')}_{meta.get('b')}_c{meta.get('c')}_{meta.get('d')}"
        if family == "F5":
            return f"a{meta.get('a')}_b{meta.get('b')}_s{meta.get('s', -1)}"
        return "?"

    # -- suites -------------------------------------------------------------

    def run_family(self, family: str, param_ranges: Dict) -> List[ProbeResult]:
        """
        Sweep one family over a Cartesian product of parameter ranges.

        `param_ranges` maps parameter name -> iterable of values; see
        families.generate_family_sequence.  Implemented for real now, where
        it previously returned an empty list and every caller got nothing.
        """
        out = []
        for inst in generate_family_sequence(family, param_ranges):
            result = self.probe_instance(inst)
            self.results.append(result)
            out.append(result)
        return out

    def run_f1_sequence(self, a_values: List[int],
                        n_range: Tuple[int, int]) -> List[ProbeResult]:
        """
        F1 sweep over `a_values` and the inclusive index range `n_range`.

        The range is used; it used to be ignored in favour of a hardcoded
        range(1, 10).
        """
        lo, hi = n_range
        out = []
        for a in a_values:
            for n in range(lo, hi + 1):
                result = self.probe_instance(generate_f1(a, n))
                self.results.append(result)
                out.append(result)
        return out

    def run_random_suite(self, k_values: List[int], n_bits_values: List[int],
                         num_instances: int = 10) -> List[ProbeResult]:
        """F2 and F3 sweeps, seeded per instance so runs are reproducible."""
        out = []
        for k in k_values:
            for n_bits in n_bits_values:
                for i in range(num_instances):
                    seed = self.seed + i
                    for gen in (generate_f2, generate_f3):
                        try:
                            inst = gen(k, n_bits, seed=seed)
                        except (RuntimeError, ValueError) as exc:
                            # A generator that cannot build the instance says
                            # so; skipping silently would hide it.
                            print(f"  skip {gen.__name__}(k={k}, n_bits={n_bits}, "
                                  f"seed={seed}): {exc}")
                            continue
                        result = self.probe_instance(inst)
                        self.results.append(result)
                        out.append(result)
        return out

    def run_f4_suite(self, params: List[Tuple[int, int, int, int]]) -> List[ProbeResult]:
        out = []
        for a, b, c, d in params:
            result = self.probe_instance(generate_f4(a, b, c, d))
            self.results.append(result)
            out.append(result)
        return out

    def run_f5_suite(self, params: List[Tuple[int, int]]) -> List[ProbeResult]:
        out = []
        for a, b in params:
            result = self.probe_instance(generate_f5(a, b))
            self.results.append(result)
            out.append(result)
        return out

    # -- reporting ----------------------------------------------------------

    def _fmt(self, lo: Optional[float], hi: Optional[float], width: int = 21) -> str:
        if lo is None and hi is None:
            return "-".rjust(width)
        lo_s = f"{lo:.4f}" if lo is not None else "?"
        hi_s = f"{hi:.4f}" if hi is not None else "?"
        return f"[{lo_s},{hi_s}]".rjust(width)

    def print_results_table(self) -> None:
        if not self.results:
            print("(no results)")
            return
        print()
        print("log2(1/|S|) is bracketed, not a point value: the exact sign plus a")
        print("certified interval bound it on both sides.  'cert bits' is the width")
        print("of the certificate itself and says nothing about |S|.")
        print()
        print(f"{'fam':<4} {'instance':<14} {'k':>2} {'N':>6} {'sign':>5} "
              f"{'method':<18} {'log2(1/|S|)':>21} {'ratio':>13} {'cert bits':>10}")
        print("-" * 100)
        for r in self.results:
            print(f"{r.family:<4} {r.instance_id:<14} {r.k:>2} {r.N_bits:>6} "
                  f"{r.exact_sign:>+5} {r.exact_method:<18} "
                  f"{self._fmt(r.log2_inv_S_lo, r.log2_inv_S_hi):>21} "
                  f"{self._fmt(r.ratio_lo, r.ratio_hi, 13):>13} "
                  f"{str(r.certificate_width):>10}")

    def save_results_csv(self, filepath: str) -> None:
        if not self.results:
            return
        rows = [asdict(r) for r in self.results]
        with open(filepath, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="SRS probe (exact decisions)")
    parser.add_argument("--family", choices=list(FAMILY_SIGNATURES) + ["all"],
                        default="all")
    parser.add_argument("--digits", type=int, default=50,
                        help="bracket precision hint; refinement is automatic")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-bits", type=int, default=None,
                        help="cap on fractional bits used by both deciders; "
                             "an instance needing more raises SignUnknownError "
                             "for that instance only")
    parser.add_argument("--output", help="write results to this CSV")
    args = parser.parse_args(argv)

    runner = ProbeRunner(seed=args.seed, digits=args.digits,
                         max_bits=args.max_bits)

    # Each family runs under its own guard.  They used to share one try block,
    # so the first instance that hit the precision cap aborted every family
    # after it and the run looked like a short but complete sweep.
    suites = (
        ("F1", "F1 (Pell, sign proved in closed form):",
         lambda: runner.run_f1_sequence([2, 3, 4, 5, 6], (1, 10))),
        ("F2", "F2/F3 (random / greedy):",
         lambda: runner.run_random_suite([2, 3], [20, 30], num_instances=3)),
        ("F4", "F4 (product composition):",
         lambda: runner.run_f4_suite([(2, 3, 5, 7), (2, 8, 3, 12), (6, 10, 7, 11)])),
        ("F5", "F5 (squaring):",
         lambda: runner.run_f5_suite([(3, 2), (5, 3), (7, 5)])),
    )

    undecided = 0
    for family, header, run_suite in suites:
        if args.family not in (family, "all"):
            continue
        print(header)
        try:
            run_suite()
        except SignUnknownError as exc:
            # The cap was reached.  That is a real limit being reported, so the
            # run is marked as not fully decided and the exit code says so --
            # previously the count was computed and then discarded, and `main`
            # returned 0 either way, making an undecided run indistinguishable
            # from a clean one.
            print(f"  {family}: instance beyond the precision cap: {exc}")
            undecided += 1

    runner.print_results_table()

    if args.output:
        runner.save_results_csv(args.output)
        print(f"\nwrote {len(runner.results)} rows to {args.output}")

    if undecided:
        print(f"\n{undecided} suite(s) stopped on the precision cap; "
              f"results above are the instances that did decide, not a "
              f"complete sweep.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
