# Irrational

Two bodies of work: a bracketed zoom method for square roots (three papers,
`Docs/` + `Src/`), and an exact-arithmetic toolkit for the Square-Root Sum
Problem (`SRS/`, with a Rust port in `Rust/`).

**`Docs/Paper1.md` (audit)** takes the original draft, finds and fixes its
errors, and measures the method honestly against certified peers — reporting
that it loses.
**`Docs/Paper2.md` (theory)** proves the mathematics: the estimator family is
Möbius, so every iterate has a closed form, and the classical sharp bound
falls out of it.
**`Docs/Paper3.md` (Möbius–CF)** relates the iterates to the continued fraction
of `sqrt(a(a+1))`. Its proof is now written out (§2.4–2.5: the residue, cross-
term and determinant lemmas by a joint induction, then the unified identity
`zₙ = αqₙ/pₙ` with the 1/a reduction by parity); the computational check for
`a <= 100` is retained as corroboration. An earlier version rested on a false
matrix identity (`M² = aC + (a−1)I`, false at a = 2) and a misstated even-n
numerator item; both are corrected in place.

**`SRS/`** decides the sign of a square-root sum `S = SUM sqrt(a_i) - SUM
sqrt(b_j)` in exact integer arithmetic, with no float anywhere on a decision
path. The zero test is in P; the sign for `S != 0` is an **open** problem, and
this repo is a labelled negative measurement of it plus proved separation
results — **not a solution**.

The **exponent** of the separation bound is settled here: applying Schmidt's
Subspace Theorem to the *single* linear form the problem is gives
`E/N <= (r-1+delta)/r`, which is exactly the rate the pigeonhole construction
achieves from below. So the sharp conjecture is `sup E/N = (r-1)/r` and the
two sides meet; what is open is the effectivity of one input-dependent
constant. The repo's `O(2^r)` bound is the sharp statement about what a
*norm* argument can do, and remains the best fully effective bound here.
`refutation.py` hunts actively for a counterexample to the sharp form and
finds none, with certified signs and magnitudes. Start at
`SRS/THE_MISSING_THEOREM.md`.
**`Rust/srs/`** ports the exact core: `BigInt` field arithmetic, the Pell
family with its residue law, and a certified sign decision cross-checked
against the closed form.

## Layout

| Path | Contents |
|---|---|
| `Docs/Draft1.md` | The original draft (Hinglish). Source of the method and of the claims audited here. |
| `Docs/Paper1.md` | The audit paper: 13 sections, 2 appendices, Theorems 1–8. Corrections, certified benchmarking, performance. |
| `Docs/Paper2.md` | The theory paper: Möbius conjugacy, exact closed form, sharp-bound attainment, the $u$-aware certificate, general-$f$ rates, and §9 on where the structure stops for $n$-th roots. |
| `Docs/Paper3.md` | Möbius iterates vs the continued fraction of $\sqrt{a(a+1)}$. Proof written out in §2 (unified identity $z_n = \alpha q_n/p_n$); computational check retained as corroboration. |
| `SRS/` | **Square-Root Sum Problem** toolkit: exact sign, certified brackets, F1–F6 instance families, measurement harness, refutation hunt. 315 tests. Docs: `SRS/THE_MISSING_THEOREM.md` (roadmap + the open constant), `SRS/STATE_OF_THE_ART.md` (literature check + corrections), `SRS/PROOFS.md` (proofs), `SRS/README.md`, `SRS/TODO.md`. |
| `Rust/srs/` | Rust port of the exact core: `BigInt` multiquadratic field, Pell family + residue law, certified sign, `is_squarefree`, CLI probe. 75 tests. |
| `Src/bracket_zoom.py` | Algorithms 1–3. Type-agnostic: accepts `int`, `Fraction`, or `Decimal`. |
| `Src/test_bracket_zoom.py` | 35 tests covering §3–§6. |
| `Src/test_u_aware.py` | 24 tests covering §11, including two rigour regressions and the optimal-certificate null result. |
| `Src/exp_precision.py` | §10 — naive vs per-level scheduled vs Newton. |
| `Src/exp_hybrid.py` | §10 — float-prefix hybrid. |
| `Src/exp_u_aware.py` | §11 — the $u$-aware bound: derivation, exact verification, tightness. |
| `Src/exp_depth.py` | §11 — depth and cost table. |
| `Src/exp_compare.py` | cubic vs plain vs Newton: where the time actually goes. |
| `Src/exp_mirror.py` | the mirror hypothesis: reflection across r, and what it buys. |
| `Src/test_exp_mirror.py` | 11 tests: perfect-square termination, bracket validity, the unsound halving. |
| `Src/exp_certified.py` | §10.1 — bisection, interval Newton, zoom, Newton+res, all exact-integer certified. |
| `Src/exact_decide.py` | Exact decision procedures: decides cubic inequalities on rational intervals, with witnesses. |
| `Src/verify_theorems.py` | The proved / decided / measured ledger — which claims are which. |
| `Src/mutation_test.py` | Injects 14 faults into the working code and reports which the tests catch. |
| `Src/test_certified.py` | 28 tests: `certify()` must reject bad brackets; convergence orders; the seed asymmetry; every §10.1 bug pinned. |
| `Src/test_exact_decide.py` | 19 tests: the decision procedures, cross-checked against brute force and direct expansion. |
| `Src/opcount.py` | Exact executed-arithmetic counts — machine independent. |
| `Src/test_opcount.py` | 9 tests pinning the corrected op-count claims. |
| `Src/mobius.py` | The Möbius structure behind Paper2's Proposition 1 — exact rationals only. Also `normed_error` / `position_blind_constant` / `strict_slack` for Corollary 7. |
| `Src/run_all.py` | **Start here.** Registry of every runnable harness + a single runner. |
| `Src/test_mobius.py` | 29 tests: closed form vs iteration, no sampling; sharpness attainment; the α<1 degeneracy; Corollary 7 (slack identity, disjointness of zero sets, the `2/(27(2+t))` constant). |
| `Src/mobius_n.py` | The $n$-th root map: error identity, multiplier, and the $1/n$ degradation (Paper2 §9). |
| `Src/test_mobius_n.py` | 19 tests: Theorem 9 exactly for $n$=2–5; sign flip; repelling region. |
| `Src/test_decimal_oracle.py` | 5 tests: `Decimal` vs `math.isqrt` — closes the libmpdc assumption. |
| `Src/README.md` | Findings index and test-to-section map. |

## Run

```bash
# Square-root sum problem toolkit (SRS/)
cd SRS
python3 run.py --list
python3 run.py            # tests + fast harnesses
python3 run.py --all      # everything, incl. all 14 Src harnesses (~15-20 min)
python3 -m pytest tests -q     # 347 tests, 12093 subtests

# Square-root zoom method (Src/)
cd ../Src

# Everything, with one command.  --list shows what exists and which paper
# section each harness reproduces; --all includes the slow sweeps.
python3 run_all.py --list
python3 run_all.py            # tests + the fast harnesses
python3 run_all.py --all      # adds exp_certified and mutation_test
python3 run_all.py certified  # just the harnesses matching a substring

# Or drive them individually, which is what the paper's Appendix A does:
python3 -m unittest test_bracket_zoom test_u_aware test_exp_mirror \
        test_certified test_exact_decide test_opcount test_mobius \
        test_mobius_n test_decimal_oracle -v      # 179 pass
python3 verify_theorems.py                       # proved/decided/measured

# Rust port of the SRS exact core
cd ../Rust/srs && cargo test && cargo clippy --all-targets
```

## Tooling rules

These were learned the hard way during this work and are worth preserving:

* **Exact rationals for proving, `mpmath` for counting.** Use
  `fractions.Fraction` whenever a claim is being *checked* — small quantities
  where rounding decides the outcome. Use `mpmath` with a fixed
  `with mp.workdps(...)` when counting iterations; `Fraction` cubes every
  level and blows up. Never mutate `mp.dps` globally.
* **Same foot, same shoe.** Any comparison of two bounds must pass both
  through the *same* evaluation routine with the same denominator. An earlier
  version divided one side by $2a(2a+w)$ and left the other bare, inflating an
  apparent gain to "100% tighter".
* **Never search for a maximiser.** A ternary search returns a value *below*
  the true maximum and therefore yields something that is not a bound. Use
  monotonicity: the denominator $D(u)=2a(2a+w)+u(2a+u)$ is strictly
  increasing, so it is lower-bounded at the left endpoint and the numerator is
  maximised in closed form.
* **Cap the precision schedule.** A per-level budget must be
  `min(need, target)`; uncapped, one level ran at 7992 digits for a
  3000-digit target — more than the naive path — which inverted the
  comparison being tested.
* **`unittest.main()` belongs at the end of the file.** Placed mid-file it
  silently skips everything defined below it; six tests were not running.
* **Derive a seed's precision from the input, never from the target.** The
  bracket came from `D(alpha).sqrt()` at full target precision, floored to an
  integer — a 3025-digit square root producing a four-digit value, 58% of
  every call, worth 2.4x of the gap to Newton. `alpha.adjusted() + 5` suffices,
  because `sqrt(alpha) - n >= 1/(2*sqrt(alpha)+1)` unless alpha is a perfect
  square.
* **Sample the points your argument names, or the test is vacuous.** The
  strictness in Corollary 7 is a *disjointness* claim — the two summands of the
  slack vanish at `eta = 2/3` and `eta = 0` respectively — but `2/3` is not of
  the form `k/40`, so on a `k/40` grid the square term's only zero is never
  hit, the zero-set collection stays empty, and the disjointness assertion
  passes on nothing. Adding the point explicitly is what makes it bite; the
  bug was found by mutation-testing the new test, which is the argument for
  mutation-testing a new test.
* **Certify the output, not just the invariant.** The interval
  `[est-eps, est]` is exact algebra and the invariant held in 48 of 48 runs —
  but reading those endpoints straight off finite-precision decimals dropped
  the rounding uncertainty of `est` and **16 of 48 certificates failed**. The
  proof and the printed number are different objects; certify the number.
* **Tolerance tests are strict inequalities.** `eps.adjusted() <= -digits`
  accepts every value in `[10^-D, 10^(1-D))`, a bracket up to 10x wider than
  promised; `alpha=257` at 50 digits landed there. Use `<`.
* **Give every method a schedule matched to its own order.** Interval Newton
  run with the `rate=1` schedule shrank its enclosure by ~16 digits per step
  and needed 65 iterations for 1000 digits; at `rate=2` it needs 9.7. A bad
  schedule makes a sound method look absurd, and the first benchmark run did
  exactly that.
* **Give every method the same seed, or report the asymmetry.** Newton was
  seeded from a 30-digit `sqrt` while the zoom and interval Newton started at
  `floor(√α)` — worth 1.66× at 3000 digits, enough to turn a 1.71× gap into
  parity. The zoom cannot take a float seed without changing what it
  computes, so it is the one method that must report the disadvantage.
* **Interleave timed methods, and keep verification out of the timed path.**
  Machine load moved every measurement by ~1.8× during this work; timing
  methods one after another charges that drift to whichever ran last. A
  `certify()` call inside a method's convergence test cost Newton+res 2.7 ms
  of a 9.8 ms row. Cite ratios, not milliseconds.

## Status

The draft's method is sound and its exact algebra is clean. Of the audit's
ten findings, several draft claims were wrong (a missing divergence guard, a
loose constant, a misattributed benchmark). The draft's own performance
hypothesis — that precision scheduling would close the gap to Newton — was
implemented and **refuted**: $2.90\times$ from scheduling at 3000 digits, with
Newton still $3.64\times$ faster.

Chasing that refutation found a larger defect: the initial bracket was
computed as a full-target-precision `sqrt` and then floored to an integer —
**58% of every zoom call**, worth $2.4\times$ of the gap. Fixed; the gap went
from $6.97\times$ to $3.64\times$ at 3000 digits.

Those numbers pit certified work against an uncertified baseline, so §10.1
re-ran the comparison against certified peers — bisection, interval Newton,
Newton with a residual certificate — every bracket proved by exact integer
arithmetic with no reference square root (320 of 320 runs, 480 of 480 on a
wider sweep). **The result refutes the claim we set out to make**: on
identical seeds the method is $417\times$ faster than bisection, but interval
Newton is $2.61\times$ faster than it at 3000 digits, and the zoom ranks third
of four. The apparent $1.71\times$ gap to Newton turned out to be **entirely a
seeding artefact** — give Newton the same integer seed the other three use and
the two are at parity. Building the benchmark surfaced a termination test 10x
looser than its stated tolerance, a bracket that certified before being printed
and failed after, an under-scheduled competitor, and that seeding asymmetry.

All 179 `Src/` tests pass, and all 315 `SRS/` tests pass.

## Rigor: what is proved, decided, and only measured

The three are kept distinct throughout, and `Src/verify_theorems.py` runs the
ledger. This matters because the obvious move — "run 30 000 cases, report zero
violations" — is a statement about the *sample*, not about the universally
quantified claim it is offered to support. That framing has been withdrawn.

* **Proved** — universally quantified, by argument. The error identities, the
  sharp constant $4/27$, Theorem 4's bound, and the $u$-aware bound's validity
  and dominance. Short written proofs; no sampling stands in for them.
* **Decided** — a claim about a specific artifact (a bracket, a run, a bug),
  settled exactly in integer/rational arithmetic with a witness on failure.
  `Src/exact_decide.py` decides "is $n(u)/D(u)\le B$ on $[\ell,h]$?" as a cubic
  inequality on a closed interval, evaluating at endpoints and the real roots of
  the derivative (rational, or quadratic irrational, both exact). The
  certified benchmark's brackets are decided in integer arithmetic by
  `certify()`. When a bound is too tight the decider returns a point where the
  inequality actually fails, so a negative answer is hand-checkable.
* **Measured** — wall-clock performance. These are measurements of specific
  implementations on one machine; they are *not* theorems and are not used to
  establish any claim about the method. Only ratios are quoted.

Building the decision procedures to replace the sampling was not a
rephrasing: `exact_decide.py` carried three defects, each caught only by
cross-checking against brute force, and each of which had made the decider
*accept invalid bounds*. A test that only ever fed it valid bounds would have
passed all three.

**Mutation testing** (`Src/mutation_test.py`) asks whether the tests check
anything at all: it injects 14 single-edit faults and reports which the suite
catches. **12 of 14 are caught**, including all three defects that shipped
inside `exact_decide.py` — a test that only ever fed the decision procedure
valid bounds would have passed every one of them. The 2 survivors are declared
in `EXPECTED_SURVIVORS` with their reasons, and printed rather than hidden,
because both are informative. (They used to make the harness exit 1 on a
healthy tree, which trained the reader to ignore the exit status. It now exits
non-zero only for an *unexpected* survivor or a mutation pattern that no longer
applies, and it prints when an expected survivor stops surviving, so the list
cannot rot unnoticed.)

* reverting the seed to depend on the target is a *performance* regression, and
  a correctness suite cannot see slowness by construction;
* replacing the corrected estimate with the plain one leaves every *property*
  the suite checks satisfied — and the obvious test for it, `est > r`, is
  actually **false** of the computed value: Corollary 2 is exact, and once the
  error drops below the working precision the returned estimate rounds across
  the root. The certificate is unaffected; the point estimate alone does not
  carry the one-sidedness.

Known limits, stated rather than papered over: nothing is machine-checked
(no proof assistant); `libmpdec`'s correct rounding is trusted rather than
verified against a second implementation; and the performance ratios are as
strong as one Python process on one machine.
