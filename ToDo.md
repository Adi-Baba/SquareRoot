# To-Do

A changelog of finished work, and what remains. Full reasoning for any item
lives in `Docs/Paper1.md` (§ references below) or `Readme.md`; this file
records *what was decided*, not the argument for it.

Done

### Draft audit and paper

- [X] Audit `Docs/Draft1.md` — 10 findings (missing divergence guard, loose
  `1/16` constant, misattributed benchmarks, non-monotonic ratios,
  omitted `isqrt` baseline, …). Corrected statements in paper §9.
- [X] Write `Docs/Paper1.md` — 13 sections, 2 appendices, Theorems 1–8.
- [X] **Split into two papers.** `Paper1.md` is the audit (corrections,
  certified benchmarking, performance). `Paper2.md` is the theory (Möbius
  conjugacy and its corollaries, sharp-bound attainment, the $u$-aware
  certificate, general-$f$ rates). The Möbius proof is stated once, in
  `Paper2`; `Paper1` §5.4 carries the consequence and a pointer, so the
  argument is not duplicated. This also retired "Proposition 9" and its
  three corollaries from `Paper1`, which are now Theorem 2 and Corollaries
  4-6 of `Paper2`.

### Task 1 — precision scheduling (§10)

- [X] Per-level precision scheduling in `Src/exp_precision.py`: **2.90×** at
  3000 digits, not the predicted ~7.6×. The draft's hypothesis is
  **refuted**; §10 rewritten against measurement.
- [X] Capped the schedule (`min(need, target)`) — uncapped, a level ran at
  7992 digits for a 3000-digit target and *inverted* the comparison.
- [X] **Fixed the full-precision seed** (the largest defect found): `zoom`
  computed `D(alpha).sqrt()` at the *target* precision then floored it —
  58% of every call. Now `alpha.adjusted() + 5`. Gap 6.97× → **3.64×**.
- [X] Float-prefix hybrid re-measured after the seed fix: now a
  **pessimization** (7.37 vs 5.71 ms at 3000 digits).

### Task 2 — the `u`-aware estimator (§11, Theorem 8)

- [X] Bound that uses the root's position inside the bracket — free, since
  `m` and `est` are already computed. Median `0.27`–`0.62×` the published
  bound; never worse. `Src/exp_u_aware.py`.
- [X] Same-footing (both sides through one bound routine), closed-form
  numerator maximiser (never a search — a search returns an under-estimate
  and certifies nothing), exact verification, depth/cost table.
- [X] `Src/test_u_aware.py` — 20 regression tests.

### Task 3 — certified benchmark (§10.1)

- [X] `Src/exp_certified.py`: four bracket-carrying methods on identical
  footing. `certify()` proves `lo² ≤ α ≤ hi²` and `hi − lo ≤ 10^-D` in
  exact integer arithmetic, **no reference sqrt**. 320/320 certified;
  480/480 on a wider sweep.
- [X] **Refuted the claim we set out to make.** On identical seeds the zoom is
  **2.61× slower than interval Newton** at 3000 digits despite fewer
  steps (7.6 vs 11.2); **417× faster than bisection**, which is not a
  serious competitor at 3.322 steps/digit.
- [X] Fixed: termination test 10× too lax (`<= -digits` admits `[10^-D, 10^(1-D))`; α=257 @ 50 digits lands there) — 1/64 → 0/64.
- [X] Fixed: bracket certified before printing and failed after — reading
  `[est-eps, est]` off finite-precision decimals dropped `est`'s rounding
  uncertainty. 16/48 certificates failed; `_est_bracket` (outward
  enclosure) → 48/48.
- [X] Fixed: certification inside the timed region (+2.7 ms of a 9.8 ms row).
  The legitimate in-loop check is one subtraction `hi_b - lo_b > 10^-D`;
  without it the digit counter falls a step short (30/480 fail).
- [X] Found: interval Newton under-scheduled at `rate=1` (65 iterations vs
  9.7). The first benchmark run was unfair to it.
- [X] Found: **seeding asymmetry worth 1.66×**, which moved the ranking. On
  the shared integer seed the zoom is at *parity* with Newton, so the
  headline gap to Newton was a seeding artefact, not an algorithmic one.
- [X] Found: machine load (~1.8×) was being charged to whichever method ran
  last. Harness now interleaves methods per α with best-of-N.
- [X] `Src/test_certified.py` — 28 tests, including that `certify()` must
  *reject* malformed brackets.

### Task 4 — sampling replaced by exact decision (§11)

- [X] **Withdrew "0 violations in 30 000 brackets" as a proof.** It described a
  sample, not a universally quantified claim. Retained only as a labelled
  cross-check that the code matches the closed form the proof uses.
- [X] `Src/exact_decide.py`: decides `n(u)/D(u) ≤ B` on `[ℓ,h]` as a cubic
  inequality on a closed interval, via Horner at the endpoints and the real
  roots of the derivative (rational or quadratic irrational, both exact).
  Rejections carry a **witness**.
- [X] **Three defects in the decider itself**, each found only by
  cross-checking, each of which made it *accept invalid bounds*: Horner
  iterating the coefficient list backwards; a spurious `B²u` term; and a
  sign error in the quadratic normalisation that reflected every critical
  point through the origin, so the true maximiser was never checked (it
  accepted bounds as small as half the true maximum).
- [X] `Src/verify_theorems.py` — the ledger separating **proved / decided /
  measured**. Paper §11 rewritten around that distinction.
- [X] `Src/test_exact_decide.py` — 19 tests pinning the deciders against brute
  force, direct expansion, and hand-checked signs.

### Task 5 — substantiating "mutation-tested" (§8.1)

- [X] `Src/mutation_test.py`: 14 single-edit faults into a scratch copy,
  **12 killed**. Includes all three `exact_decide` defects, both
  `certify()` directions, `D` monotonicity, the `2w/3` critical point, the
  interval-Newton schedule rate, the perfect-square guard, and the seed
  asymmetry.
- [X] **Found a bug I had claimed to have pinned and had not.** Reverting the
  termination test to `<=` *survives*, because the scheduled zoom
  overshoots its target by ~24 orders of magnitude, so both comparisons
  agree on every exercised input. Pinned now by an explicit source-contract
  check, labelled as such. Same for the `4/27 → 1/16` swap: output error
  there is rounding-limited (~1e-85), not bound-limited, so the constant
  is unobservable through that path.
- [X] Two survivors reported rather than hidden: the seed regression (a
  *performance* bug — a correctness suite cannot see slowness) and dropping
  the cubic correction (both estimators satisfy every property checked).
- [X] Found while trying to kill the second: **`est > r` is false of the
  computed value.** Corollary 2 holds exactly, but at α=2, 30 digits the
  returned estimate lands *below* r — once the error drops below the
  working precision the recomputed estimate rounds across the root. The
  certificate is unaffected (it is an outward enclosure, not a point).
  Added as a remark at Corollary 2.

### Task 6 — the op-count claim, measured (§10, §10.1)

- [X] `Src/opcount.py`: counts **executed** arithmetic as (line executions) ×
  (per-statement ops), so it follows helper calls and separates divisions
  from multiplications. Exact and machine independent — unlike wall clock,
  an op count belongs in the *decided* column.
- [X] **§10.1's stated cause was backwards.** The paper said "three divisions
  and three multiplications against interval Newton's two divisions".
  Measured: the zoom does **fewer** divisions per level (11.2 vs 19.5) and
  **more** total ops (88.8 vs 72.0). What costs it is precision-weighted
  volume — **6.5×** the digit-operations at 3000 digits.
- [X] **§10's claim was wrong too.** "Three divisions per level against
  Newton's one division per step" counted only the inline level body and
  omitted `_prec_for` / `_width_digits` / `_seed`. Measured: zoom 99.8 ops
  and 11.7 divisions per level against Newton's 25.0 ops and **9.0**
  divisions per step.
- [X] Recorded that the proxy is **not** a time model: digit-operations
  overpredicts the measured gap by ~2.5×. It explains the direction and
  rules out the original explanation; the magnitude stays a measurement.
- [X] `Src/test_opcount.py` — 7 tests, including one that fails if the stale
  "more divisions" claim becomes true again, and one pinning the
  overprediction.

### Task 7 — closing the stated assumptions (§12)

- [X] **`libmpdec`'s rounding cross-validated** against `math.isqrt` — an
  independent exact-integer oracle — at 50/200/1000/3000/6000 digits over
  11 values of α (`Src/test_decimal_oracle.py`). Also checks that
  ROUND_FLOOR/ROUND_CEILING are genuinely directed, and that
  `as_integer_ratio()` is lossless. **This was listed as an open
  assumption for three tasks; it is now measured, not trusted.**
- [X] **The $0<\alpha<1$ claim was wrong and is withdrawn.** §5.3 said the
  method extends there "with no rescaling". In exact arithmetic it does
  not: at $a=0$ the corrected estimator degenerates to $\hat r=1$ for
  every $\alpha$, and the progress condition becomes $w<0$. The *bound* is
  valid at $a=0$; the *recursion* is not. §5.3 corrected in place,
  §12 updated, and `TestAlphaBelowOne` pins both facts.
- [X] **Removed a sampling check I had introduced.** `test_mobius` verified
  the $4w^3/27$ maximiser on a 60-point grid — the exact reasoning Task 4
  withdrew from the paper. Replaced with the algebraic identity
  $\frac{4w^3}{27}-u^2(w-u)=(u-\frac{2w}{3})^2(u+\frac{w}{3})\ge0$.

### Task 8 — non-attainment of the sharp constant (Paper2 §5.4, Corollary 7)

- [X] **A proposed "position-aware optimal certification" paper was declined,
  and the reason is recorded because it is the useful part.** The
  construction is Theorem 8 verbatim in both papers (interval `[m−a, r̂−a]`,
  maximise `u²(w−u)` over it), and the sharp certificate `max_{L≤η≤U} F(t,η)`
  beats it by a factor `1.000000` at median *and* worst, because `h−ℓ` is
  already the current error — median `3.7e-4·w`. Sharpening a bound below
  the resolution of the interval it is computed on is not a result.
- [X] The theorem asked for — "`C < 1/27` except near `eta = 2/3`" — is
  **false as stated**. There is no exception: the inequality is strict for
  every `eta ∈ [0,1]` and every `t > 0`. `1/27` is a supremum over a
  *corner* of the state space (`t=0`, `eta=2/3`), which no bracket reaches
  since `w/a > 0` and `0 < u < w` strictly.
- [X] **Corollary 7**, from the exact factorisation
  `27η²(1−η) − (4+2t+2tη+t²η²) = −27(η−2/3)²(η+1/3) − 2t(1+η) − t²η²`:
  the slack against `2/(27(2+t))` is a sum of two non-negative terms whose
  zero sets (`eta=2/3`, `eta=0`) are **disjoint**, hence strictly positive.
  Replaces Corollary 6's asymptotic attainment factor with an exact
  non-attainment statement.
- [X] **The sharp position-blind constant is `2/(27(2+t))`, not `1/27`** — and
  `exp_u_aware.bound_published` already computes it (it is `4w³/27 / D(0)`).
  The papers were carrying the loose constant while the code used the sharp
  one; neither recorded the gap. Relative gain `t/(2+t)`: 33% at `w/a=1`,
  4.8% at `1/10`, 0.05% at `1/1000`.
- [X] **The "unwieldy quartic maximiser" of §8.3 is `2/3 − t/27 + t²/243`**, and
  `sup F = 1/27 − 5t/162 + 16t²/729`. So the ratio maximiser and the
  numerator maximiser agree to zeroth order and part company at `O(t)`;
  the residual uncertainty is `O(t)` about a known point, not an unlocated
  quartic root.
- [X] `Src/mobius.py` (`normed_error`, `position_blind_constant`,
  `strict_slack`), 5 new tests in `Src/test_mobius.py`, a `PROVED` row in
  `Src/verify_theorems.py`. **179 pass.**
- [X] **"Quantify the average and worst-case reduction" — done, and the answer
  is zero** (Paper2 §5.5, `Src/exp_u_aware.py` §(2b),
  `measure_optimality_gap`). `C(t;L,U) = max_{L<=eta<=U} F(t,eta)` is the
  optimal position-aware certificate; evaluated by 64-fold interval
  refinement against Theorem 8's closed form on the same interval through
  the same routine, `C/B` = **1.0000 median, 1.0000 worst** at `w<=a/100`
  and 0.9979 median at `w<=a`.
- [X] **The instructive part is a trap, and it is now pinned.** In exact
  rationals **95 of 120** visited states *do* improve — so a count reads
  like a result — but by at most `5e-6` relative, median exactly zero, and
  it vanishes entirely in floats. Asserting on frequency rather than
  magnitude would have "confirmed" a null result as a win; my first test
  did exactly that and failed. `TestOptimalCertificate` now asserts on the
  magnitude (`< 1e-4`), not the count.
- [X] **Structural reason, not just numerical.** The bracket handed to the
  certificate is `[m-a, rhat-a]` and its width `h-l = rhat-m` *is the
  current error* — the interval over which `C` maximises is the quantity
  being bounded, so refinement buys resolution far below the quantity's own
  size.
- [X] **The §5.4 progress-guard corollary was asserted, then decided.** With
  the sharp `ε` the bracket-shrinking guard widens from `w < √27·a ≈ 5.20a`
  to `w < 3a(9+√129)/4 ≈ 15.27a`. Stated in §5.4 as a fact about the
  guard while nothing tested it against the recursion — the exact failure
  mode this repo keeps warning about. Now `TestProgressGuard` (3 tests,
  rationals): `ε` is valid, both guard roots compared exactly, and running
  the recursion over its 16 radicands gives max `w/a = 1` at level 0 and
  `≤ 1e-2` after — so the guard is not the binding constraint anywhere the
  method reaches. 3 of 5 mutations killed; the 2 survivors are provably
  implied by the 360 per-state assertions, the third (`seen > 0`) is a
  vacuity guard that holds and so cannot fail.
- [X] **Mutation-tested the new tests, and the new tests had a bug.** Deleting
  `2/3` from the sample grid leaves the disjointness assertion vacuous (the
  square term's only zero is never hit) and *survives* in `test_mobius.py`
  while being caught in `verify_theorems.py`. Fixed by asserting the
  zero-set equalities, not just collecting them. 8/8 algebraic faults
  killed (wrong factor `η+1/3→η+1/2`, mixed term dropped, summands
  swapped, `27→26`, constant reverted to `1/27`, three slips inside `F`).

### Also investigated

- [X] **Cost analysis** (`Src/exp_compare.py`): the plain (uncorrected)
  estimate measures *faster* than the cubic one at 1000 and 3000 digits,
  and needs no divergence guard (plain always contracts). See Open.
- [X] **Mirror hypothesis** (`Src/exp_mirror.py`): `x ↦ α/x` is a reflection,
  and `r = GM < AM` gives one free halving. Result: 1.11× at 1000 digits,
  depth 7.6 → 7.2 — not the missing 3.6×. Punchline: `(x + α/x)/2` *is*
  Newton's step, so the idea converges to the competitor. Two bugs fixed
  (perfect squares hung; the free half is only valid at seeding).
- [X] Housekeeping: `unittest.main()` mid-file (6 tests never ran), corollary
  renumbering, `1/16` misattribution, cross-reference audit.

## Open

Ordered roughly by how much they weaken the work.

- [X] **Where the structure stops: the $n$-th root** (Paper2 §9). Error
  polynomial is the complete homogeneous symmetric polynomial
  $h_{n-2}(z,r,a)$, constant only at $n=2$; hence degree $n-1$ and **no**
  rational linearisation for $n\ge3$ (degree argument). Multiplier
  $\lambda_n\to-(n-1)$, so the scalar iteration is repelling unless the root
  is known to within a relative $\sim1.59/n$. Alternation survives for
  all $n$ (Cor 12) — convergence does not. `Src/mobius_n.py` +
  `Src/test_mobius_n.py` (19 tests, exact rationals). The companion-pencil
  construction is noted as classical and explicitly *not* claimed as a
  remedy.
- [X] **Exact closed form for the iteration** (Paper2 Prop. 1, Thms 2-3; this
      was the retired "Prop. 9" of an earlier draft). The
  estimator family is Möbius; $T(z)=\frac{z-r}{z+r}$ conjugates it to a
  scaling with $\lambda=-u/(2a+u)$, so the $n$-step iterate is exact.
  Subsumes the error identity, the $4/27$ constant and the $2w/3$
  maximiser; makes one-sidedness exact for *every* iterate.
  `Src/mobius.py` + `Src/test_mobius.py` (29 tests, exact rationals, no
  sampling; five of them are Task 8's Corollary 7). Two of my own test bugs
  found and fixed on the way.
- [ ] **Derive and test a 4th-order estimator** (digit quadrupling against
  tripling). Never attempted — nothing in `Src/` or the paper touches it.
  Read the caveat first: §10.1 measured the zoom losing on
  *precision-weighted volume* (6.5× the digit-operations of interval
  Newton at 3000 digits), not on divisions — it already does fewer
  divisions per level (11.2 vs 19.5). A 4th-order variant buys *levels*,
  not cheaper levels, so on the metric now in use it would likely lose by
  more, not less. The per-level cost has to come down first, and
  `_est_bracket` (two directed-rounding quotients per level) is the
  obvious target. Test the cheap version before the general one: does the
  degree-4 closed form even exist, and what is its error constant?
- [ ] **Nothing is machine-checked.** The universal claims rest on written
  argument; the finite ones on `exact_decide.py`. Neither is verified by a
  proof assistant. A Lean/Coq or SymPy encoding of the two-line Theorem 4
  argument and the `u`-aware bound would close this. Brute-force
  cross-checking caught the decider's three defects, but brute force is
  not a proof assistant.
- [ ] **Why does digit-operations overpredict the measured gap by ~2.5×?**
  The proxy charges every op its full working precision; `Decimal` does
  not. A real cost model (`M(D)` for mul/div at `D` digits) would turn
  the op-count from a direction claim into a quantitative one, and would
  settle the next item.
- [ ] **Attack the op count, or keep the negative finding.** §10.1 now
  attributes the loss to precision-weighted volume (6.5× the
  digit-operations at 3000 digits) rather than to divisions. Either reduce
  that — the certificate extraction `_est_bracket` builds two
  directed-rounding quotients per level and is the obvious target — or
  leave it as the stated negative result.
- [ ] **Decide the fate of the cubic correction.** `exp_compare.py` measures
  the plain estimate as *faster* at 1000 and 3000 digits, with no
  divergence guard needed. If that holds, the paper's headline feature is
  a pessimization and the framing has to change.
- [ ] **Does the `u`-aware bound generalise to §7's arbitrary `f`?** It rests
  on `D'(u) > 0`, an artefact of `x² − α`. For general `f` the denominator
  may not be monotone and only a piecewise bound survives. The decider
  extends to any polynomial error form, so this becomes decidable once `f`
  is fixed.
- [ ] **§11 depth table is one trace per digit target**, not an ensemble. The
  "saves 0–2 levels" figure should be measured across many α before it is
  quoted as anything but illustrative.
- [ ] **Re-measure the draft's §8.5 tables.** Still the draft's own timings,
  arithmetic-audited but not re-run. Either re-run or keep labelling them
  quoted.
- [ ] **Hybrid's weakened guarantee.** The float prefix is certified only from
  the first `Decimal` level. Acceptable for publication, or should a
  float-rounding error term be carried explicitly?

## Verification

```bash
cd Src

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
```
