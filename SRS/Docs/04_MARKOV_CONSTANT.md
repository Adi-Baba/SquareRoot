# Paper D — The `r = 2` Constant, Done Right

Version 1.0. This paper states the `r = 2` constant with the correct
definition, proves the needed lemmas, and explains why the earlier "M(α) ≤ 1"
argument (kept in this repo for some time) was wrong. It does **not** re-prove
Hurwitz's theorem; that is quoted.

---

## D.0 The two definitions, told apart

Two Markov-type constants live in this area and are easy to conflate:

```
λ(α)  :=  inf_{q ≥ 1}  q·‖qα‖,                (the separation constant)
M(α)  :=  sup_{n ≥ 0}  q_n·‖q_nα‖,  convergents  (a limsup)
```

They serve different purposes.

- `λ(α)` is the *best constant in the uniform lower bound*
  `|qα + p| ≥ c/q`. The largest such `c` is exactly `λ(α)`.
- `M(α)` is a *measure of the supremum of the convergent errors*. It is what
  one reads off a continued fraction directly. It is **always** `< 1` for
  non-trivial reasons, and it is what the "F1 achieves rate 1/2" computation
  actually uses.

**They are not equal.** In particular `M(α) ≤ 1` does **not** imply that
`1/M(α)` is a lower-bound constant. The repo's earlier claim "`M(α) ≤ 1`, so
the `r = 2` constant is `≥ 1`" committed exactly this conflation. The correct
constant is `λ(α)`, and computing it requires the infimum, not the supremum.

---

## D.1 The separation constant `λ`

**Lemma D.1.** For any irrational `α` and any integers `p, q` with `q ≥ 1`,

```
|qα + p| ≥ λ(α) / q.
```

*Proof.* Choose `p` to minimize `|qα + p|`; then `|qα + p| = ‖qα‖`, so
`q|qα + p| = q‖qα‖ ≥ inf_{q'} q'‖q'α‖ = λ(α)`. ∎

Conversely any `c < λ(α)` fails for sufficiently good convergents, so `λ(α)`
is the best possible such constant.

**Lemma D.2 (reduction of general pairs to lowest terms).** Write `p/q = P/Q`
in lowest terms, so `p = kP`, `q = kQ`. Then

```
|qα + p| = k·|Qα + P| ≥ k·λ(α)/Q = k²·λ(α)/q ≥ λ(α)/q.
```

So Lemma D.1 already covers every integer pair. There is no extra factor of
`k`: the worst case is `k = 1`, i.e. a reduced denominator. ∎

---

## D.2 Hurwitz, and what it costs

**Hurwitz's theorem (quoted).** For every irrational `α`,

```
λ(α) ≤ 1/√5,
```

and `1/√5` is attained exactly when `α` is equivalent to the golden ratio
`(1+√5)/2`.

Combining with Lemma D.1: the constant `c₀ = λ(α)` in `|qα + p| ≥ c₀/q`
satisfies

```
c₀ ≤ 1/√5.
```

So the *optimal* `c₀` is at most `1/√5`; equivalently the *slack*
`log₂(1/c₀) = −log₂ λ(α)` is always `≥ log₂√5 ≈ 1.16` bits. In words:
**the `r = 2` separation bound costs at least 1.16 bits of slack**, and this is
sharp for the worst-case `α`. It is not free — and it is not zero.

However `c₀ > 0` for every quadratic irrational `α` (bounded partial
quotients, i.e. badly-approximable), so `log₂(1/c₀)` is a *finite* quantity
that can be computed from the period of `α`'s continued fraction. Hence:

**Theorem D.3 (`r = 2` is effective, and the constant needs no extra
hypothesis).** For two squarefree radicands `d₁ ≠ d₂`, `α = √(d₂/d₁)` has
bounded partial quotients, so `λ(α) > 0` is computable from the period, and

```
log₂(1/|S|) ≤ N/2 + log₂(1/λ(α)) + O(1),
```

with `0 < λ(α) ≤ 1/√5`. The additive cost is a computable constant
independent of `N`. So `r = 2` needs no separate Markov constant assumed *in
addition* to the algorithm — only that the CF of `α` be periodic, which
Lagrange's theorem guarantees for every quadratic irrational.

---

## D.3 Why the old `M(α) ≤ 1` claim was wrong

The former statement was: "`M(α) ≤ 1` for every irrational `α`, so the `r=2`
constant is `≥ 1` and costs no bits." Three distinct errors are packed in that
sentence.

1. **Wrong direction of the inequality.** `λ(α) ≤ 1/√5 < 1` is the true
   half of Hurwitz. The claim "`c₀ ≥ 1`" is the opposite.

2. **Wrong witness topology.** `M(α)` is a supremum *over convergents*. The
   convergent denominators are exactly the reduced denominators `q`, so by
   D.2 every convergent pair is already covered by `λ`. Hence `M(α)` and
   `λ(α)` are tied to the same arithmetic, and one cannot conclude `c₀` from
   `M(α)`.

3. **A fabricated exact witness.** The earlier "`q = 3161, p = 2449`, value
   `2.4996`" used to prop up `M(√(5/3)) > 1` is arithmetically impossible
   (`3161·√(5/3) ≈ 4080.83`, so the nearest integer is `4081`, not `2449`).
   No such witness exists. The clean fact is the other direction: the
   convergent at the `0`-th convergent `p/q = ⌊√D⌋/1` already gives
   `q·‖qα‖ = ‖√D‖ = √D − ⌊√D⌋ < 1`, i.e. the *sup over convergents* is
   always `< 1`, and the *inf* tells us `λ` can be much smaller.

The valid takeaway is precisely: at `r = 2` the constant is computable and the
case is closed; it is never "free" in the sense of costing zero bits.

---

## D.4 What the r=3 bound would need

The `r = 2` success rests on `α` being one quadratic irrational, so `|qα+p|`
has a one-variable lower bound. At `r = 3` the form is
`C₁ + C₂√(d₂/d₁) + C₃√(d₃/d₁)`, with *two* independent generators. The
analogue of `λ` would be an infimum over two integer parameters, and no
Hurwitz-scale universal lower bound is available. That is the single reason
`r = 3` is harder, and why Theorem D.3 does not generalize naively.

---

## D.5 Summary

- Use `λ(α) = inf_q q‖qα‖`, never `1/M(α)`.
- `0 < λ(α) ≤ 1/√5`; the `r = 2` slack is always at least `log₂√5`.
- The cost is a computable constant for quadratic `α` (Lagrange), so `r = 2`
  is closed but **not** for free.
- Nothing in this paper is a uniform bound across radicands; it is per-`d`.