# Paper 35 — M7.12 sharpened: the exact exponential gap, quantified

Version 1.0.  `Docs/34` left M7.12 qualitative ("prove `‖β_jε^k‖ ≥ k^{-C}`").
This paper replaces that with a **quantified exponential gap of a factor
`9n(n−1)`**, identifies the one mechanism that would close it, and shows that
the natural competitor (Ridout / Subspace-Theorem power laws) is *strictly
worse* and therefore useless here.

---

## 35.1 Setup and notation

Distinct squarefree `a₁ ≠ a_j`, `ε = 1+√a₁`, `ε' = 1−√a₁`, `ρ = ε'/ε`,
`β_j = √(a_j/a_1)/2`, and
```
L  =  Q(√a₁, √a_j)  =  Q(√a₁, √(a₁a_j)),        [L : Q] = 4 .
```
The object of study is the distance to the nearest integer of
```
x_k  :=  β_j ε^k .
```
By `Docs/27` Lemma 27.1 (corrected, `Docs/34` §34.1), `‖x_k‖ = ‖q_k√a_j‖ + O(ε^{-k})`
with `q_k ≍ ε^k`, so `x_k` **is** the SRS quantity.

## 35.2 Theorem 35.1 — the norm (Pisot-type) lower bound

> **Theorem 35.1.**  For every `m ∈ Z` and every `k ≥ 1`,
> ```
> | m − β_j ε^k |   ≥   ε^{−3k} / ( 2 β_j³ )   =   |ε'|^{3k} / (2 β_j³)
>                                          =   exp( −3k·log ε + O(1) ) .
> ```
>
> **CORRECTION (`Docs/36` §36.1).**  This was first printed as `|ρ|^{3k}/(2β_j³)`.
> That is **wrong**: `|ρ| = |ε'/ε| = 0.1716`, whereas `ε^{−1} = |ε'| = 0.4142`, so
> `|ρ|^{3k} ≠ ε^{−3k}`.  The norm argument gives `ε^{−3k}` (the surviving
> conjugate product is `~ 2β_j³ε^{3k}`).  The prose, the gap analysis of
> Thm 35.4, and the numerical verification all used `log ε` and are unaffected;
> only the displayed `|ρ|` was wrong.  It mattered in the *optimistic*
> direction — the true barrier `ε^{−3k}` is **stronger** than the erroneous
> `|ρ|^{3k}`.

*Proof.*  `y := m − β_jε^k ∈ L`.  **Correction (`Docs/36` §36.1b):** `y` is
*not* an algebraic integer, so `N_{L/Q}(y)` is a **rational**, not an integer —
at `k=1, a₁=2, a_j=3` it is `228/256`.  What *is* an integer is
`N_{L/Q}(2a₁y) = 4a₁m − P_k√(a₁a_j) − a₁Q_k√a_j` being an algebraic integer.
Hence `|N(y)| ≥ (2a₁)^{−4}` rather than `≥ 1`, which changes only the leading
constant:  The conjugates of `y`
are `m ∓ β_jε^{±k}` (§34.3), i.e. of moduli `|m−β_jε^k|`, `|m+β_jε^{-k}|`,
`|m+β_jε^k|`, `|m−β_jε^{-k}|`.  For fixed `m,k` at most one factor is `<1`, and
since `y ≠ 0` (β_jε^k irrational) the norm is a nonzero integer, so `|N| ≥ 1`.
Dividing by the other three factors — each `≤ 2β_jε^k + |m|`, and for the
minimising `m = ‖·‖` we have `|m| ≤ β_jε^k` — gives
`1 ≤ |y|·(2β_jε^k)·(β_jε^k)·(2β_jε^k)`, hence the claim. ∎

**Exact verification** (integer Pell pairs, 60-digit magnitudes, `a₁=2, a_j=3`,
`k = 2…41`, including `q = 1.7469e15`): the quantity
`‖x_k‖·2β_j³·ε^{3k} ≥ 1` holds at **every** point, minimum `39.18`.  **Zero
violations.**

## 35.3 Theorem 35.2 — the automatic-Pisot mechanism is absent

For a Pisot number `θ` and `λ ∈ Q(θ)`, one gets **for free**
`‖λθ^k‖ ≈ |λ'|^k → 0`, because `λθ^k + λ'θ'^k ∈ Z`.  This is the entire
Pisot mechanism, and it is exactly what would make M7.12 *false*.

> **Theorem 35.2.**  `β_j ∉ Q(√a₁)`  (for distinct squarefree `a₁ ≠ a_j`).

*Proof.*  `β_j = √(a₁a_j)/(2a₁)`, so `β_j ∈ Q(√a₁)` ⟺ `√(a₁a_j) ∈ Q(√a₁)` ⟺
`a₁a_j` is a rational square ⟺ `a₁ = a_j` (distinct squarefree). ∎

Checked directly for `a₁=2, a_j=3`: `β_j = √6/4`; if `√6/4 = r+s√2` then
`rs = 0` and `r² = s² = 3/8`, neither a rational square. ∎

**Consequence.**  `x_k` is a **trace-zero** element whose conjugate `−β_jε^k`
has *equal* modulus (`Tr_L(x_k) = 0`, Thm 35.3).  There is no Pisot-type
"automatic integrality": smallness of `‖x_k‖` must be earned, not automatic.
This is the clean algebraic reason `Docs/30`'s pseudo-Pisot hypothesis fails.

## 35.4 Theorem 35.3 — the recurrence is order-2 and trace-zero

> **Theorem 35.3.**  `x_{k+2} = 2x_{k+1} + x_k` (since `ε² = 2ε+1`), and
> `Tr_L(x_k) = 0` for all `k`.

*Proof.*  The recurrence is immediate.  For the trace, the conjugates of
`β_jε^k` are `{β_jε^k, −β_jε^{-k}, −β_jε^k, β_jε^{-k}}` (each `±` twice), which
sum to `0`. ∎

So M7.12 is precisely: **how small can the fractional part of a trace-zero
order-2 linear recurrence get?**  This is a classical, well-posed question, and
it is the honest final form of the problem.

## 35.5 Theorem 35.4 — THE SHARP GAP

> **Theorem 35.4.**  `Docs/29`'s threshold to exclude is `‖β_jε^k‖ > θ₀^k` with
> `θ₀ = ε^{-1/(3n(n−1))}`, i.e. `exp( −k·log ε / (3n(n−1)) )`.  Theorem 35.1
> delivers only `exp( −3k·log ε )`.  **The proven bound is weaker than the
> needed one by a factor `9n(n−1)` in the exponent**, so it excludes nothing:
>
> ```
> 3 log ε   /   ( log ε / (3n(n−1)) )   =   9n(n−1) .
> ```
>
> For `n = 2` this is a factor **18**, i.e. proven `exp(−2.6442k)` versus needed
> `exp(−0.1469k)`.

This is exactly the "short by `9n`" of `Docs/28` §28.3, now made exact: the
factor is `9n(n−1)`, not `9n`.

## 35.6 Corollary 35.5 — what would close M7.12, and what will not

> **Corollary 35.5.**  It suffices to prove, for each `j`, a bound
> ```
> ‖β_j ε^k‖   ≥   exp( −k·log ε / (3n(n−1)) + O(1) )   for all large k,
> ```
> i.e. to improve Theorem 35.1's exponent by the factor `9n(n−1)`.  **No**
> polynomial/height argument can do this (Theorem 34.1).

**And the obvious competitor fails, provably.**  One would expect
`‖β_jε^k‖` — a 3-term linear form `q_j√a_j + p_j√(a₁a_j)/a₁ − m` with
coefficients `≍ ε^k` — to be attacked by Ridout / Subspace-Theorem power laws,
which give `≥ q^{-1-ε} = exp(−(1+ε)k·log ε)`.  But
```
(1+ε)·log ε   <   3·log ε        (ε < 2),
```
so **Ridout is strictly weaker than the elementary norm argument here**, by a
factor `≈3`.  The norm bound wins because it exploits the *Pell relation*
coupling `p_k` and `q_k`; a generic 3-term-form treatment throws that away.

This is a real and useful conclusion: **the standard Diophantine-approximation
toolkit is exhausted, and the one remaining algebraic structure — the Pell
coupling — is already fully used by Theorem 35.1.**

## 35.7 Empirical note (not a proof)

For `a₁=2, a_j=3`, `k = 2…41`: `‖x_k‖ ∈ [0.0040, 0.4881]`, and
`min_k ‖x_k‖ = 0.00404` at `k=41` versus `θ₀^{41} = 0.00242`.  **No exception
occurs anywhere in the tested range**, and the norm bound of Thm 35.1 is never
close to tight (minimum slack `39.18`).  So the exception set is plausibly
*empty*; what is missing is a proof that it is.

## 35.8 Final M7.12 statement

> **M7.12 (final, quantified).**  For distinct squarefree `a₁ ≠ a_j`, the set
> ```
> { k ≥ 1  :  ‖ (√a_j / 2√a₁)·(1+√a₁)^k ‖   ≤   (1+√a₁)^{−k/(3n(n−1))} }
> ```
> is finite, with the largest element `k*` satisfying `k*·log ε ≤ poly(N)`.
>
> **Best available:** `‖·‖ ≥ exp(−3k log ε)` (Thm 35.1) — short by a factor
> `9n(n−1)`.  **Blocked on:** an argument improving the exponent by that factor.
> Heights, Liouville, and standard Subspace-Theorem power laws are each proved
> insufficient (`Docs/34` Thm 34.1, `Docs/35` §35.6).

**Assessment.**  The problem is now a sharply quantified open question about
the fractional parts of a trace-zero order-2 recurrence, with a clean
`9n(n−1)`-factor gap.  I am not claiming a solution; the record now states
precisely what is missing and why the standard methods do not reach it.