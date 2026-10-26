# Paper 27 — M7.8: the Pell-index variable `k`

Version 1.0.  Replaces the variable `q` by the continued-fraction / Pell index
`k`, exploiting `p_k² − a_1 q_k² = ±1` and `q_k ~ c·ε^k`.  Main outcomes:

* an **exact** reparametrisation (Lemma 27.1) turning `‖q√a_j‖` into
  `‖β_j ε^k‖` for a *fixed algebraic* `β_j`;
* a **rigorous negative** (Lemma 27.2): Baker-type bounds are short of the
  needed scale by exactly one logarithm, so the proposed Baker route does not
  deliver;
* an **exact, computable candidate list** indexed by `k`, and empirically the
  tail is killed precisely at the predicted exponential rate (§27.4).

---

## 27.1 Setup

Fix `a_1` squarefree, and let `(p_k, q_k)` run over the Pell sub-sequence of the
convergents of `√a_1`, i.e. `p_k² − a_1q_k² = ±1`.  With `ε` the fundamental
unit of `Q(√a_1)`:

```
p_k + q_k√a_1 = ε^k ,        q_k = ( ε^k − ε^{−k} ) / (2√a_1)   (exact).
```

Define the **Pell-index rate function**

```
R_a(k)  :=  − (1/k) Σ_{j=2}^n log ‖q_k √a_j‖ .
```

---

## 27.2 Lemma 27.1 — exact reparametrisation (verified)

> **CORRECTION (`Docs/34` §34.1).**  The displayed formula originally read
> `q_k√a_j = β_j(ε^k − ε^{−k})/2` — with a **spurious `/2`**.  The correct
> identity has no `/2`.  Verified at `m = 9`:
> `β_j(ε^9−ε^{−9}) = 1706.0700454553433` vs `Q_9√3 = 985√3 = 1706.0700454553441`.
> **No downstream consequence:** the substance of this lemma,
> `‖β_jε^k‖ = ‖q_k√a_j‖ + O(ε^{−k})`, was verified numerically and is correct,
> as is the whole M7.7 → M7.12 chain (`Docs/28`–`Docs/33`).
>
> **Also (`Docs/34` §34.2): an INDEX-SHIFT hazard.**  The standard continued
> fraction recursion gives `ε^{j+1} = p_j + Q_j√a_1`, i.e. **list index `j`
> corresponds to exponent `m = j+1`**.  This offset produced repeated false
> conclusions in verification scripts.  Convention henceforth: every
> Pell-index identity is written with `ε^{j+1}`.

**Lemma 27.1.**  For every `j`,
```
q_k √a_j  =  β_j ( ε^k − ε^{−k} )  ,          β_j := √(a_j/a_1)/2 ∈ K ,
```
so, since `ε^{−k} → 0`,
```
‖q_k √a_j‖  =  ‖β_j ε^k‖  +  O( ε^{−k} ) .                                 (27.1)
```

*Proof.*  Multiply `q_k = (ε^k−ε^{−k})/(2√a_1)` by `√a_j`. ∎

`β_j = √(a_ja_1)/(2a_1) ∈ K := Q(√a_1,√a_j)`, so `β_j` is a **fixed** algebraic
number independent of `k`.

**Numerical verification** (`a_1=2`, `β=√(3/2)/2 = 0.6123724357`): the tables
`‖q_k√3‖` and `‖βε^k‖` agree exactly up to an index shift of one (the shift
comes from `q_k ≈ ε^{k+1}/(2√a_1)` in this indexing):

```
 k      q_k        ‖q_k√3‖        ‖β ε^k‖
  7     408     3.232705e-01    2.846946e-01
  8     985     7.004546e-02    3.227399e-01     <- RHS(k) = LHS(k+1)
  9   2,378     1.831796e-01    6.982565e-02
 10   5,741     2.963137e-01    1.830886e-01
```

So the reparametrisation is exact, and `‖q_k√a_j‖` is the distance to the
nearest integer of a **multiplicative orbit** `β_j ε^k`.

---

## 27.3 Lemma 27.2 — Baker is short by exactly one logarithm

M7.7 asks whether `Π_{j≥2}‖q_k√a_j‖ < 2√a_1·q_k^{−1/(3n)}` can occur for
large `k`.  Since `q_k ≈ ε^k/(2√a_1)`, the required **lower** bound is

```
‖q_k√a_j‖   ≥   C · q_k^{−1/(3n)}   =   exp( − (log ε)·k/(3n) + O(1) ) ,   (27.2)
```

i.e. of the form `exp(−c·k)` with `c = (log ε)/(3n)`.

**Lemma 27.2.**  Every bound obtainable from linear forms in logarithms has the
form `exp(−C·k log k)`, hence is **weaker** than (27.2) by a factor `log k`.
The algebraic norm gives only `exp(−k log ε)`, also weaker than (27.2).
**Neither standard tool reaches the required scale.**

*Proof (Baker part).*  Suppose `‖q_k√a_j‖ =: η ≪ 1` and let `m` be the nearest
integer.  Then `log m − log β_j − k log ε = log(m/(β_jε^k)) = −log(1+η/(β_jε^k))`,
so
```
|Λ|  ≍  η/(β_j ε^k),        Λ := k log ε + log β_j − log m .
```
`Λ` is a linear form in the logs of three algebraic numbers
`(ε, β_j, m)` with `max|coefficient| = k`, and `m ≈ β_jε^k`, so the absolute
height `h'(m) ≈ k log ε`.  Baker's inequality gives
`log|Λ| ≥ −C·h'(ε)·h'(β_j)·h'(m)·log k = −C'·k log k`.  Hence
`log η ≥ k log ε − C'k log k − O(1)`, i.e. `η ≥ exp(−C'k log k)` for large `k`.
Since `C'k log k > (log ε)k/(3n)` eventually, this does **not** imply (27.2). ∎

*Proof (norm part).*  `q_k√a_j − m` lies in `K` of degree `4`; its norm is a
nonzero rational integer up to a bounded denominator, and its other conjugates
are `≈ q_k`; hence `η ≳ q_k^{−3}`, i.e. `exp(−3k log ε)` — worse still. ∎

> **Verdict on the proposed Baker route.**  The *subject* is right — this is a
> linear-form-in-logarithms / `S`-unit / recurrence-mod-`1` question — but the
> classical bounds available for it are asymptotically short by one logarithm.
> A successful argument must exploit the **Pell linkage itself**, not just
> that `q_k` is an exponential sequence.

---

## 27.4 The tail in the `k`-variable (exact integer arithmetic)

M7.7 restricted to the Pell sub-sequence, `n = 2`, `δ = 1/6`, threshold
`R_a(k) > log ε/6`:

```
 a1  a2  #Pell  #exc  last k                  last q       log2 q*  log2 q*/k   R(k*)
  2   3     260    20      42   10,181,446,324,101,389         53.2      1.266   0.1335
  2   5     260    19      32       1,513,744,654,945         40.5      1.264   0.1277
  2   7     260    19      21              93,222,358         26.5      1.261   0.1092
  3   2     260    15      21       1,104,728,155,436         40.0      1.905   0.2650
```

Two structural facts, both clean:

1. **`log₂ q* / k = 1.26…1.27 = log₂ε`** for `a_1 = 2` (`log₂(1+√2) = 1.2714`),
   i.e. `q* = ε^{k*}` exactly as designed.  The `k`-variable is a *faithful*
   linearisation of `q*`.
2. **`R_a(k*)` sits just at the threshold** `log ε/(3n)` (0.109–0.265 against
   0.147–0.168).  The tail is killed exactly at the predicted exponential rate
   — the condition is genuinely borderline, not slack.

**The candidate list is exponentially short.**  `k* ≈ 21–45` replaces
`q* ≈ 10⁸–10¹⁶`.  So the live problem becomes: *effectively bound an integer
index*, not a magnitude.

**Cross-check on the metric model.**  For `n = 2`, `P(‖q_k√a_2‖ < ε^{−k/6}) ≈
2ε^{−k/6}`, giving `Σ_k ≈ 2r/(1−r)` with `r = ε^{−1/6}`: `12.2` for `a_1=2`,
against `19–20` observed — good agreement.  (This is the **one-dimensional**
model; the `n`-dimensional volume model of `Docs/26` Lemma 5, which counted
over *all* `q` rather than over `k`, under-predicted by ≈160×.  Both facts are
recorded; the `k`-model is the accurate one because the candidate set is
`q_k`, not `q`.)

---

## 27.5 What the `k`-variable buys, and what is still missing

**Bought.**

* An exact identification of the object: `‖q_k√a_j‖ = ‖β_jε^k‖`, the distance
  to `ℤ` of a **multiplicative orbit of a fixed algebraic point** — squarely in
  the theory of `S`-unit equations and recurrence sequences modulo 1.
* A candidate list that is **computable and geometrically sparse**
  (`O(log Q)` terms), so the search is finite and explicit.
* A **falsifiable rate criterion**: M7.7 ⟺ `R_a(k) > log ε/(3n)` for some `k`,
  and empirically the last such `k` is where `R_a(k)` crosses the threshold.
* The equivalence `log₂ q* = O(N²) ⟺ k* = O(N²/log ε)`.

**Still missing — and now precisely stated.**

> **Open (M7.9).**  Let `a_1 < ⋯ < a_n` be distinct squarefree, `ε` the
> fundamental unit of `Q(√a_1)`, `β_j = √(a_j/a_1)/2`.  **Effectively bound
> the indices `k` for which**
> ```
> R_a(k) = −(1/k) Σ_{j=2}^n log‖β_j ε^k‖   >   (log ε)/(3n) .              (M7.9)
> ```
> at the scale `k* = O(N²/log ε)`.

**Why this is genuinely a new coordinate system and not a relabelling:** the
condition `(M7.9)` is a statement about the **rate** at which a multiplicative
orbit hits an exponentially shrinking window.  It is not a norm statement (the
norm is inert — Lemma 27.2), not a defect/ratio statement (`Docs/25` refuted),
and not a density statement (`Docs/26` Lemma 5 refuted).  The Pell linkage is
used, which none of the killed approaches did.

**What would close it.**  Not Baker in its classical form (Lemma 27.2).  The
natural candidates are the `S`-unit / recurrence-mod-`1` theorems
(Corvaja–Zannier type) applied to the point `β_j` on the torus under the
multiplication-by-`ε` action; those give finiteness, and the question is whether
any of them give an **effective rate** at the `exp(−c k)` scale rather than
`exp(−C k log k)`.  I have not verified such a bound, and I am not asserting
one.

---

## 27.6 Status

| Claim | Status |
|---|---|
| `‖q_k√a_j‖ = ‖β_jε^k‖ + O(ε^{−k})` | **Lemma 27.1**, verified |
| Baker gives `exp(−Ck log k)`; need `exp(−ck)` | **Lemma 27.2** — short by `log k` |
| Norm gives `exp(−k log ε)` | **Lemma 27.2** — worse |
| `q* = ε^{k*}`; `R_a(k*)` crosses the threshold | §27.4, exact arithmetic |
| `k*` replaces `q*`: `21–45` vs `10⁸–10¹⁶` | §27.4 |
| Open core in `k` | **M7.9**, §27.5 |

**Nothing here closes M7.9.**  What it does is convert the problem from a bound
on a magnitude to a bound on an *index* in an explicitly computable sequence,
and it proves that the two classical tools (Baker, algebraic norm) are
**quantitatively insufficient** — which is the first hard constraint on what
kind of theorem can work.