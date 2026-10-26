# Paper 7 — Error Log

Version 1.0. Every error found and retracted, in one place, including the
ones made during the repair pass. Each entry records: the false claim, why it
was false, the actual truth, and where the correction now stands.

The point of this file is that the same mistakes recur; naming the exact step
at which each argument dies is the only way to stop reinventing them.

---

## E1. "`M(α) ≤ 1` implies the `r = 2` constant is `≥ 1` and free"

**False claim.** In an earlier revision of this repository, the closing of
`r = 2` was justified by "`M(α) ≤ 1` for every irrational `α`, so the constant
at `r = 2` is `≥ 1` and costs no bits."

**Why false.** It conflates two different constants. `M(α)` (sup over
convergents) is indeed `< 1`, but that does *not* bound the separation
constant. The correct separation constant is `λ(α) = inf_q q‖qα‖`, and Hurwitz
gives `λ(α) ≤ 1/√5`, i.e. `log₂(1/λ) ≥ log₂√5`. The constant is a real,
positive, computable slack — it is never zero.

**Current status.** Paper D states `λ` correctly and cites Hurwitz. The `r = 2`
case is still closed, by Paper D and Paper E, but it is **not** free.

**Related audit errors.** Two measurement errors compounded this:

- The supporting "witness" `q = 3161, p = 2449, value 2.4996` for
  `M(√(5/3)) > 1` is arithmetically impossible (`3161·√(5/3) ≈ 4080.83`).
  The clean witness formula `M(√D) = √D − ⌊√D⌋` at the `0`-th convergent
  gives `M(α) < 1` trivially, so the "counterexample to `M ≤ 1`" cannot be
  repaired — the statement `M ≤ 1` is actually *true*; the error was never in
  `M`, it was in inferring the separation constant from it.
- The supporting "exact" argument cited `max_{q≤400} q‖qα‖` as the statistic.
  That max is `~Q/2` for *every* `α` (a generic `q` has `‖qα‖ ≈ 1/2`), so it
  carries no information about `α` and is not a Markov constant. The relevant
  statistic is the infimum over reduced denominators, not the max over all
  denominators.

---

## E2. The `dist(aW, Z) = a·dist(W, Z)` lattice identity

**False claim.** A general effective constant for `2 ≤ r ≤ 4` rested on

```
min_{C_r ∈ Z} |C_r + a_j·W| = a_j · dist(W, Z),
```

where `W` is a point in a number field and `a_j ∈ R_{>0}`.

**Why false.** `dist(a·W, Z) = a·dist(W, Z/a)`, i.e. the lattice on the right
is rescaled to `Z/a`, which is a different lattice from `Z` unless `a ∈ {±1}`.
So the "`= a·dist(W, Z)`" step is invalid for every irrational `a`, which is
the generic case.

**Witness.** `d = (2,3)`, `C = (1,−1)`, `S = √2 − √3`. The claimed bound
would give `|S| ≥ √3 = 1.7321`; the true value is `|√2 − √3| = 0.3178`. The
inequality fails, and the failure is caused by the lattice rescaling.

**Current status.** The claim is retracted in Paper F and pinned by a test
(`test_RETRACTED_general_bound_is_false`) that exercises the inequality
directly and *fails* it, so the retraction cannot be silently reversed.

---

## E3. `sup E/N = (r−1)/r` stated as a supremum over all `N`

**False claim.** The files stated the sharp form as
`sup_{S ≠ 0} E/N = (r−1)/r` over all instances.

**Why false.** That is a *pointwise-over-all-N* statement, and it admits small
excursions. The certified counterexample is

```
S = −10√2 − 98√3 + 51√13,   |S| = 2.8318252869537366×10⁻⁷,
```

with `E ≈ 21.75`, `N = 17`, `E/N ≈ 1.2795 > 2/3`. The excess is the additive
slack `log₂(1/c)` divided by a small `N`, not a rate phenomenon.

**Correct statement.** The sharp form is the asymptotic one:

```
limsup_{N→∞} E/N = (r−1)/r,
```

which is what the Subspace Theorem (Paper C) plus the pigeonhole construction
(Paper B) together establish. The finite-`N` excess is `O(log(1/c)/N)`.

**Current status.** `sup` replaced with `limsup` in the current files; the
counterexample instance is recorded in Paper F and is the anchor of the error
log.

---

## E4. `bracket_num` described as a "bit-truncation"

**False claim.** An earlier version of the docs called `signexact`'s
`bracket_num` "a bit-truncation" and implied it was unsound.

**Why false.** `signexact._sum_bracket_int` computes a true enclosure
`S ∈ [lo/2^bits, hi/2^bits]`, and `bracket_num` is the lower endpoint of that
same enclosure. It is a sound *lower* bound, only coarser than
`intervalsum.certified_interval_for_S`'s enclosure.

**Current status.** §2.11.10's footnote records the correction.

---

## E5. Thm 2.12's rank condition, stated too strongly

**False claim.** The old statement used `r` without the support-rank
distinction, so the claimed `log₂(1/c)` budget was too small by a factor
related to the support functions.

**Why false.** The exponent that survives is `r − 1`, but the intermediate
monomial counting involves supports of rank `j − 2`, not `j − 1`; the
effective constant came out off by one in the small-support case.

**Current status.** Corrected to `2^{r−1} − 1` in Paper B and Cor B.2.

---

## E6. EHS Thm 4 stated as `E/N = O(N)`

**False claim.** The repo recorded EHS Thm 4 as giving `E/N = O(N)`.

**Why false.** That is the *shape* of a non-uniform bound, not the EHS
statement. EHS Theorem 8 gives `E ≤ 2n·N`, i.e. `E = O(nN)` with an
**ineffective** constant hidden in the Subspace-Theorem step. Writing
`E/N = O(N)` would be an overstatement of effectivity.

**Current status.** Paper C states EHS as `E ≤ 2n·N`, ineffective in the
constant.

---

## E7. My own repair-pass errors (the ones just made)

These are the errors introduced *while repairing the others*. They are the
most likely to recur.

1. **Conflating `M` and `λ` again.** My first repair of the `r = 2` constant
   treated the repo's `M(α) = limsup q_n‖q_nα‖` as if its reciprocal were the
   separation constant. It is not; Paper D defines `λ` separately and uses it.

2. **A non-reduced convergent "witness".** I briefly offered `q = 3, p = 3`
   for `α = √(5/3)` as a counterexample to `M(α) ≤ 1`. But `3/3` reduces to
   `1/1`; the pair is not a convergent denominator pair, and
   `3·|3α − 3| = 9·|α − 1|` is just the reduced convergent's value inflated by
   `k²`. So it does not count. The lesson: a witness pair must be in lowest
   terms.

3. **An unstable formula for `ρ_n`.** An early version of the F1 asymptotic
   used `(1−s_n)²/(2s_n)` where the correct identity is
   `(1−√s_n)²/(2√s_n)`. The two differ by exactly a factor of `4`, so the
   "`1/(8 ln 2)`" constant I derived was a phantom of the algebra. The correct
   identity is `ρ_n = log₂(1 + (1−√s_n)²/(2√s_n))`, and Paper E re-derives it
   correctly.

4. **`max_{q ≤ 400}` as a Markov statistic.** As above in E1, this is not the
   object; the infimum over reduced denominators is.

5. **Panic-editing a correct theorem.** During the `ρ_n` audit I almost
   "fixed" Theorem 3.2(iv), believing it was off by a factor of 4. The
   discrepancy was entirely my own wrong formula, not the paper's. The
   invariant: **before editing a proof, re-derive the quantity by hand and
   check the constant it converges to against the constant the proof predicts.**

---

## E8. The float trap behind every numeric claim in this log

The common thread: several "violations" were just float cancellation. The two
forms that cannot be evaluated in double precision are:

```
q·α − p        (at p, q ~ 1e65, this is ~1e−130 and reads as exactly 0)
s_n − 1        (Pell residue, same issue past n ≈ 10)
```

The numerically stable forms are:

```
q‖qα‖   via  q·|q²α² − p²·α| / (q·α + p)   ... but even that loses the
          subtraction; correct is |D·q² − Q·p²| / (Q(qα + p))  with
          D, Q the exact integers defining α = √(D/Q).
```

```
1 − √(s)  via  (1 − s)/(1 + √s)
```

All certified numeric statements in this repository must go through an exact
integer kernel; any "measurement" that does not is a candidate for the error
log, not for the theorems.

---

## E9. The process lesson

Each error above was found by an adversarial test, not by re-reading the
proof. The standing rule for this repository is therefore:

> A bound is not claimed until a test tries to refute it and fails, and even
> then it is labelled **proved** only after a hand derivation.

The retractions are pinned by regression tests so the corrected statements
cannot quietly revert.