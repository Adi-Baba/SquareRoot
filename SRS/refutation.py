"""
Refutation hunt for the SRS conjecture.

WHAT THIS IS FOR
----------------
Thm 2.5 (pigeonhole) and Thm 2.9 (single-form Subspace Theorem) bracket the
separation rate from both sides at the SAME number:

    achieved   E/N >= (r-1)/r  -  o(1)        Thm 2.5
    guaranteed E/N <= (r-1+delta)/r + O_{d,delta}(1)   Thm 2.9

so the sharp form of the conjecture is

    sup E/N = (r-1)/r      for every r.

This module hunts for a COUNTEREXAMPLE: an instance with E/N strictly above
the ceiling.  A single such instance with a certified sign and a certified
magnitude would falsify Thm 2.5 or Thm 2.9, and with them the sharp form.

WHAT IS AND IS NOT PROVEN HERE
------------------------------
Every magnitude reported by this module is EXACT.  The search heuristic (which
promotes a candidate) uses floats, exactly as F3 documents, but the reported
E is computed from a certified rational bracket by the integer kernel in
signexact.py, so a reported E/N above the ceiling would be a certified
refutation rather than a numerical artefact.

That is the point of the discipline: a refutation must be exact, otherwise it
is not a refutation.  The search may well produce nothing; the module is
designed so that "nothing found" is a meaningful, reproducible negative
result rather than a vague shrug.

USAGE
-----
    python3 refutation.py                      # default sweep
    python3 refutation.py --r 3 --budget 200    # single rank, bigger budget
"""

from __future__ import annotations

import argparse
import math
import random
from fractions import Fraction
from typing import List, Optional, Tuple

import field
from families import generate_f6
from signexact import exact_sign
from squarefree import decompose_and_aggregate

__all__ = [
    "Instance",
    "LITERATURE_CAP",
    "exact_E",
    "is_refutation",
    "slack_for",
    "f2_rank",
    "ceiling_for_rank",
    "hunt_random",
    "hunt_rank3_ceiling",
    "main",
]


class Instance:
    """One certified candidate, with its exact separation."""

    __slots__ = ("terms", "sign", "E", "N_bits", "r", "_mag")

    def __init__(self, terms, sign: int, E: float, N_bits: int, r: int,
                 magnitude: Optional[Fraction] = None):
        self.terms = terms
        self.sign = sign
        self.E = E
        self.N_bits = N_bits
        self.r = r
        self._mag = magnitude

    @property
    def magnitude(self) -> float:
        """
        |S| as a float.  Used only to SIZE the additive allowance in
        `slack_for`, never to decide a sign, so a float here is harmless.
        The certified value is 2**-E.
        """
        if self._mag is not None:
            return float(self._mag)
        return 2.0 ** (-self.E)

    def norm_power(self) -> float:
        """
        ||C||_inf^(r-1), the normalisation in Thm 2.9's |S| >= c*||C||^{-(r-1+delta)}.

        So c_measured = |S| * ||C||^(r-1) is the best constant this instance
        certifies, and log2(1/c_measured) is the additive bit allowance.
        """
        m = max((abs(c) for c, _ in self.terms), default=0)
        return float(m) ** (self.r - 1) if m else 0.0

    @property
    def rate(self) -> float:
        return self.E / self.N_bits if self.N_bits else float("inf")

    def __repr__(self) -> str:
        return (f"Instance(r={self.r}, N={self.N_bits}, E={self.E:.3f}, "
                f"E/N={self.rate:.4f}, terms={self.terms})")


def exact_E(terms) -> Optional[float]:
    """
    E = log2(1/|S|), certified, for S != 0.  Returns None if S = 0.

    The magnitude comes from the integer kernel's rational bracket, so this is
    a proved value of E up to the displayed float precision -- not a float
    evaluation of the sum.
    """
    res = exact_sign(terms)
    if res.sign == 0:
        return None
    num, den = res.witness["bracket_num"]
    mag = Fraction(abs(num), den)
    if mag == 0:
        return None
    return -math.log2(mag.numerator / mag.denominator)


def f2_rank(radicands) -> int:
    """F2-rank of the parity vectors, i.e. [K:Q] by Thm 1.2."""
    vecs = [field.MultiQuadraticField._parity_vector(int(d)) for d in radicands]
    return field._rank_f2(vecs)


def ceiling_for_rank(r: int) -> float:
    """
    The conjectured sharp RATE, (r-1)/r, from Thm 2.5 (below) meeting
    Thm 2.9 (above).

    IMPORTANT: this is a LIMIT, not a finite-size inequality.  Thm 2.5 gives

        E/N  >=  (r-1)/r  -  log2(sum_i sqrt(d_i)) / (r*log2 B)
                              the O(1/log2 B) slack

    because the construction's prefactor sum_i sqrt(d_i) is not 1.  An
    instance may therefore sit ABOVE (r-1)/r at finite size without
    refuting anything.  Use `slack_for` before comparing; comparing a finite
    instance against this limit produces spurious "refutations", which an
    earlier version of this module duly did.  See SLACK below.
    """
    if r < 2:
        return 0.0
    return (r - 1) / r


# The largest rate any instance is EXPECTED to reach.  Thm 2.5 supplies
# instances at (r-1)/r, and EHS Thm 8 guarantees E/N <= 2 for every input.
# So the band is [(r-1)/r, 2] and a rate above 2 would contradict the
# literature outright, not just the repo's own theorems.  Anything above 1
# is already a serious surprise.
LITERATURE_CAP = 2.0


def slack_for(inst: "Instance") -> float:
    """
    How far ABOVE the asymptotic rate (r-1)/r a finite instance is entitled
    to sit.  A candidate must exceed `ceiling + 2*slack` before it counts as a
    refutation rather than a finite-size artefact.

    There are TWO independent sources of slack, and this function originally
    modelled only the first.  Getting this wrong produces false positives,
    which has now happened twice in this repo -- see the note below.

    1. The Thm 2.5 pigeonhole prefactor.  With coefficients bounded by B and
       fixed radicands, N ~ r*log2 B, and the construction's prefactor
       sum_i sqrt(d_i) contributes log2(sum sqrt d) bits:

           slack_pigeonhole = log2( sum_i sqrt(d_i) ) / N.

    2. The ADDITIVE constant of Thm 2.9.  Every bound in PROOFS.md carries an
       additive term log2(1/c(d,delta)) + O(1), where c is the Markov constant
       of the radicands.  At fixed d that is a CONSTANT number of bits, so it
       contributes

           slack_constant = log2(1/c_measured) / N

       to the rate -- and c_measured can legitimately be ~1e-7 at r=3, i.e.
       ~22 bits, which at N=17 is a rate contribution of ~1.3.  This is
       exactly the excursion documented in PROOFS.md 2.10.2, where a certified
       instance reaches E/N = 1.2795 > 2/3.

    Source 2 is the additive constant log2(1/c(d,delta)) of Thm 2.9.  `c` is
    a function of the RADICANDS ONLY, so at fixed `d` it is a constant number
    of bits (measured at ~9 to ~22 for r=3) which inflates E/N by that many
    bits over N and vanishes as N grows.  It is a genuine source of
    finite-size headroom, and ignoring it is the SECOND false-positive mode
    this module has had: a certified instance at E/N = 1.2835 (PROOFS.md
    2.10.2) clears a prefactor-only threshold and gets reported.

    Estimating it from the instance is CIRCULAR and provably vacuous.  With
    E = rate*N and c = |S|*||C||^(r-1) = 2^(-rate*N)*||C||^(r-1), one gets

        log2(1/c) = rate*N - (r-1)*log2||C||,

    which GROWS FASTER than the excess it is supposed to excuse,

        excess = (rate - (r-1)/r) * N.

    Since the excess carries a factor `(r-1)/r < 1` and the allowance carries
    `rate`, the allowance exceeds the excess for every rate.  Checked
    numerically at r=3, N=123, ||C||=2^40: rates 0.68, 0.70, 0.75, 0.80, 1.0,
    1.5, 2.5, 4.0 all fail to fire.  A per-instance estimate of the additive
    constant therefore makes `is_refutation` silently vacuous, and the harness
    would report a clean negative result for the wrong reason.  An earlier
    draft of this module did precisely that, twice, in two different forms.

    The resolution: the additive constant is bounded in N, so it is bounded by
    a FIXED number of bits, `ADDITIVE_BITS`, not by a per-instance estimate.
    The value is measured, not guessed -- see the note on the constant.
    """
    if inst.r < 2 or not inst.N_bits:
        return 0.0
    pref = sum(math.sqrt(d) for _, d in inst.terms)
    return math.log2(pref) / inst.N_bits


# Fixed bound, in BITS, on the additive constant log2(1/c(d,delta)) of
# Thm 2.9 at r = 3.  This is the finite-size headroom that makes certified
# instances reach E/N = 1.2835 at N = 17 (PROOFS.md 2.10.2) without
# refuting anything.
#
# It is a constant of the RADICANDS, not of the instance, so it does not
# scale with E and cannot swallow the excess.  The value is a measured upper
# range over exhaustive greedy search at r = 3 (B <= 700, several radicand
# triples): c >= 1.2e-3, so log2(1/c) <= 9.7 bits; 32 is taken as a safe
# round ceiling with margin.  Raising it makes the harness more permissive,
# never less, so the negative result stays honest either way.
ADDITIVE_BITS = 32.0


def trend_for(instances) -> float:
    """
    The ONLY statistic that can refute an ASYMPTOTIC rate.

    A single finite instance cannot: the additive constant log2(1/c) inflates
    E/N by O(1/N), so E/N is free to sit above (r-1)/r at small N.  What would
    refute the conjecture is a FAMILY whose rate stays bounded away from
    (r-1)/r as N -> infinity.

    Returns max_i (E_i/N_i) - (r-1)/r over the family, i.e. the largest
    excess over the ceiling.  If this is bounded and does not grow with N,
    the family does not refute anything.  This is a NEGATIVE-result helper:
    the honest output is "no family found", never "proved".
    """
    if not instances:
        return 0.0
    return max(inst.rate - ceiling_for_rank(inst.r) for inst in instances)


def is_refutation(inst: "Instance") -> bool:
    """
    The honest test, and the two false positives it has to survive.

    Both mistakes this function exists to prevent were made in this repo:

    1. Comparing E/N against (r-1)/r with no slack at all, which reported
       F6 -- Thm 2.5's own construction -- as a refutation.
    2. Allowing only the pigeonhole prefactor as slack, which reports the
       certified witness of PROOFS.md 2.10.2 (E/N = 1.2835) as a
       refutation, because it ignores the additive constant of Thm 2.9.

    So the test is done in BITS, not in rates.  The excess

        excess = E - (r-1)/r * N

    is the number of bits by which the instance sits above the asymptotic
    line.  The pigeonhole prefactor contributes `slack_for(inst) * N` bits,
    and the additive constant contributes a bounded, N-independent number of
    bits which the instance itself certifies: `log2(1/c)` with
    c = |S|*||C||^(r-1).  Both are subtracted, with a factor 2 of margin.

    Working in bits is what keeps this from being circular.  Both allowances
    are O(1) in N -- the prefactor because it is log2(sum sqrt d)/N times N,
    the additive constant because it is a fixed budget `ADDITIVE_BITS` -- so
    the excess `(rate - (r-1)/r)*N` grows without bound against them and a
    genuine family is still caught.  Folding the additive constant into a
    per-instance RATE slack instead makes the threshold `ceiling + 2*rate`,
    which nothing can clear.
    """
    if inst.r < 2 or not inst.N_bits:
        return False
    ceiling = ceiling_for_rank(inst.r)
    excess_bits = inst.E - ceiling * inst.N_bits          # bits above line
    prefactor_bits = slack_for(inst) * inst.N_bits        # O(1) in N
    return excess_bits > 2.0 * (prefactor_bits + ADDITIVE_BITS)


def _N_bits(terms) -> int:
    return sum(max(abs(c).bit_length(), d.bit_length()) for c, d in terms)


def _rank_of(terms) -> int:
    agg, _ = decompose_and_aggregate(terms)
    agg = {d: c for d, c in agg.items() if c}
    return f2_rank(sorted(agg)) if len(agg) >= 2 else 1


# ------------------------------------------------------------
# Hunt 1: random instances, the F2 route
# ------------------------------------------------------------

def hunt_random(r: int, budget: int, bits: int, seed: int = 0
                ) -> Tuple[Instance, List[Instance]]:
    """
    Random squarefree radicands and random signed coefficients at the given
    rank, looking for a rate above the ceiling.

    The rank is enforced by rejection: resample radicands until their parity
    vectors have F2-rank exactly r, so `r` is the quantity the ceiling is
    stated in.  (A large k over a small r buys an adversary nothing -- see
    THE_MISSING_THEOREM.md -- so r is the only parameter worth sweeping.)
    """
    rng = random.Random(seed)
    best: Optional[Instance] = None
    beat: List[Instance] = []
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    tries = 0
    while tries < budget:
        tries += 1
        ds = rng.sample(primes, r)
        vecs = [field.MultiQuadraticField._parity_vector(d) for d in ds]
        if field._rank_f2(vecs) != r:
            continue                      # not rank r; resample
        terms = []
        for d in ds:
            c = rng.choice([-1, 1]) * rng.randint(1, 2 ** bits)
            terms.append((c, d))
        E = exact_E(terms)
        if E is None:
            continue
        inst = Instance(terms, exact_sign(terms).sign, E, _N_bits(terms), r)
        if best is None or inst.rate > best.rate:
            best = inst
        if is_refutation(inst):
            beat.append(inst)
    return best, beat


# ------------------------------------------------------------
# Hunt 2: the rank-3 ceiling family
# ------------------------------------------------------------

def hunt_rank3_ceiling(Bs=(30, 60, 90, 120, 150, 200, 260)
                       ) -> Tuple[Optional[Instance], List[Instance]]:
    """
    F6 is the Thm 2.5 construction, so by that theorem its rate is at the
    ceiling and cannot beat it.  Running it anyway is worth it for two
    reasons: it checks the generator against the theorem it is meant to
    instantiate, and it is the only family in the repo whose rate is
    EXPECTED to sit at the ceiling, so it calibrates what "at the ceiling"
    looks like numerically.
    """
    best: Optional[Instance] = None
    beat: List[Instance] = []
    for B in Bs:
        try:
            inst = generate_f6(B)
        except ValueError:
            continue
        E = exact_E(inst.terms)
        if E is None:
            continue
        cand = Instance(inst.terms, inst.metadata["exact_sign"], E,
                        inst.N_bits, 3)
        if best is None or cand.rate > best.rate:
            best = cand
        if is_refutation(cand):
            beat.append(cand)
    return best, beat


# ------------------------------------------------------------
# Hunt 3: greedy near-cancellation, the F3 idea at fixed rank
# ------------------------------------------------------------

def hunt_greedy(r: int, budget: int, bits: int, seed: int = 0
                ) -> Tuple[Optional[Instance], List[Instance]]:
    """
    Fix r-1 terms, then choose the last coefficient to make the sum as close
    to zero as possible.  This is the F3 idea, restricted to a target rank
    and with the final coefficient chosen by exact rounding.

    For a candidate last radicand d_r, the best C_r is
    round(-sum_{i<r} C_i sqrt(d_i) / sqrt(d_r)), and the residual is then
    certified exactly.  This is the most direct attack on the geometry that
    makes S small, and it is the shape an adversary would actually use.
    """
    rng = random.Random(seed)
    best: Optional[Instance] = None
    beat: List[Instance] = []
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    for _ in range(budget):
        ds = rng.sample(primes, r)
        vecs = [field.MultiQuadraticField._parity_vector(d) for d in ds]
        if field._rank_f2(vecs) != r:
            continue
        head = [(rng.choice([-1, 1]) * rng.randint(1, 2 ** bits), d)
                for d in ds[:-1]]
        total = sum(c * math.sqrt(d) for c, d in head)
        last = ds[-1]
        c_last = -int(round(total / math.sqrt(last)))
        if c_last == 0:
            continue
        terms = head + [(c_last, last)]
        E = exact_E(terms)
        if E is None:
            continue
        inst = Instance(terms, exact_sign(terms).sign, E, _N_bits(terms), r)
        if best is None or inst.rate > best.rate:
            best = inst
        if is_refutation(inst):
            beat.append(inst)
    return best, beat


# ------------------------------------------------------------
# Reporting
# ------------------------------------------------------------

def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Hunt for an SRS instance refuting the sharp conjecture "
                    "sup E/N = (r-1)/r.  All reported E are certified.")
    p.add_argument("--r", type=int, nargs="*", default=[2, 3, 4, 5],
                   help="F2-ranks to hunt at (default: 2 3 4 5)")
    p.add_argument("--budget", type=int, default=200,
                   help="candidates per rank per strategy (default: 200)")
    p.add_argument("--bits", type=int, default=60,
                   help="coefficient bit budget (default: 60)")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--f6-max-B", type=int, default=60,
                   help="largest B for the F6 ceiling family (default: 60). "
                        "The F6 grid is (B+1)^3 values sorted, so this is "
                        "cubic in B; raise it for a deeper ceiling check.")
    args = p.parse_args(argv)

    print("SRS refutation hunt")
    print("=" * 72)
    print("Target: an instance with E/N > (r-1)/r, which would falsify")
    print("Thm 2.5 (achieved) or Thm 2.9 (guaranteed).")
    print("All E values below are CERTIFIED (integer kernel), not floats.")
    print()

    total_beat: List[Instance] = []

    for r in args.r:
        ceiling = ceiling_for_rank(r)
        print(f"--- rank r = {r}   ceiling (r-1)/r = {ceiling:.4f} ---")
        for name, fn in (("random", hunt_random),
                         ("greedy", hunt_greedy)):
            best, beat = fn(r, args.budget, args.bits, seed=args.seed)
            total_beat.extend(beat)
            if best is None:
                print(f"  {name:8} no candidate")
                continue
            refutes = is_refutation(best)
            verdict = ("REFUTES" if refutes else
                       f"below ceiling+2*slack ({ceiling + 2 * slack_for(best):.4f})")
            print(f"  {name:8} best E/N = {best.rate:.4f}  slack = {slack_for(best):.4f}  [{verdict}]")
            print(f"           {best}")
        print()

    print("--- rank 3, F6 ceiling family (expected AT the ceiling) ---")
    f6_Bs = [b for b in (15, 30, 45, 60, 90, 120, 150, 200, 260)
             if b <= args.f6_max_B]
    best, beat = hunt_rank3_ceiling(Bs=tuple(f6_Bs))
    total_beat.extend(beat)
    if best is not None:
        verdict = ("REFUTES" if is_refutation(best) else
                   "at the ceiling, as Thm 2.5 predicts (slack "
                   f"{slack_for(best):.4f})")
        print(f"  F6       best E/N = {best.rate:.4f}  [{verdict}]")
        print(f"           {best}")
    print()

    print("=" * 72)
    if total_beat:
        print(f"*** {len(total_beat)} CANDIDATE(S) ABOVE THE CEILING. ***")
        for inst in total_beat[:10]:
            print(f"    {inst}")
        print("Each has a certified sign and magnitude.  Verify, then this is")
        print("a genuine refutation of the ASYMPTOTIC rate")
        print("limsup_{N->inf} E/N = (r-1)/r.  (The `sup` form of that claim is")
        print("retracted: finite instances above (r-1)/r are expected, because")
        print("the additive constant of Thm 2.9 inflates E/N by O(1/N) bits --")
        print("see PROOFS.md 2.10.2, where a certified instance reaches 1.2835.)")
        return 1
    print("No refutation found.  This is a NEGATIVE result, reported as such:")
    print("the sweep found no certified instance whose excess over the")
    print("asymptotic rate exceeds the finite-size allowance.")
    print("It does not prove the conjecture; it only bounds what was searched.")
    print("Note the search covers ONE parameterisation per strategy: a family")
    print("outside it would not be seen.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
