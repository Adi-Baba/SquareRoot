# Paper 32 — M7.12: can the Prop. 3.4 relations be enumerated?

Version 1.0.  Answer: **almost — the relations can be reduced to a pencil in
`L³` and then refuted one by one, with an explicit algebraic contradiction and
no case split.**  What survives is a single, clearly-identified set.  A claim in
`Docs/31` §31.2 is corrected here.

---

## 32.1 Correction to `Docs/31` §31.2

`Docs/31` wrote that the coefficients `ζ_j` of the reduced relation are
"FIXED … drawn from a finite, explicitly enumerable set".  **That is too
strong and is hereby corrected.**

`ζ_j = Σ_{𝒃 ∈ ℐ*, j_2(𝒃)=j_2} a_𝒃 · b_{j_1(𝒃),ℓ(𝒃)}`, where the `a_𝒃 ∈ L`
are the Subspace Theorem's output.  So `ζ_j` is fixed **only after the
(unknown) relation has been chosen**; it is *not* drawn from a finite
enumerable set of fixed rationals.  Lemma 31.2 is valid **per relation**, but by
itself gives no effectivity, because the relation set is not finitely
enumerable.

What *is* fixed, and is the correct observation, is the following.

---

## 32.2 New result: the exceptional space collapses to `L³`

**Theorem 32.1 (effective dimension collapse).**  In KMN's construction for
`k = 1`, `α_1 = ε`, `q_1 = β_j` **fixed**, the vectors `𝔶(k)` satisfy
```
𝔶(k)  =  ε^k · E_1  +  (ε')^k · E_2 ,        E_1, E_2 ∈ L   FIXED,       (32.1)
```
where `ε' = 1−√a_1` and `E_1, E_2` are built from the **fixed rationals**
`b_{1,j_1,ℓ}` (the `K`-coordinates of `β_j`, via their Lemma 3.1).  Hence

```
r  :=  dim_L V   ≤   2 ,        |ℐ*| = r  ≤  2 ,
```
so **the Subspace Theorem is applied in `L^{r+1}` with `r+1 ≤ 3`** — and `r`
is **exactly computable from the fixed data**.

*Proof.*  The `(j_1,j_2,ℓ)`-component of `𝔶(k)` is `b_{1,j_1,ℓ}·ε_{j_2}^k` with
`j_2 ∈ {1,2}`.  Collecting the `j_2 = 1` slots into `E_1` and the `j_2 = 2`
slots into `E_2` gives (32.1), and `span_L{𝔶(k)} ⊆ span{E_1,E_2}`. ∎

KMN do not state this; it is specific to `q_1` being **fixed**, and it shrinks
the exceptional-subspace search from an unknown high-dimensional space to a
pencil in `L³`.

---

## 32.3 The reduction, and a complete contradiction

The relation (14) becomes, with `p_k` the nearest integer to `β_jε^k`,
```
a_0·p_k  +  Ψ_1·ε^k  +  Ψ_2·(ε')^k  =  0 ,        (a_0, Ψ_1, Ψ_2) ∈ L³ .   (32.2)
```

**Lemma 32.2 (reduction to a geometric rounding error).**  Suppose (32.2) holds
for infinitely many `k`.  Put `e_k := p_k − β_jε^k` (so `|e_k| < θ^k`) and
`Tr := Tr_{K/Q}`.  Since `Tr(β_jε^k) = 0` and `Tr(p_k) = 4p_k`:

1. Dividing (32.2) by `ε^k` and letting `k → ∞` (`e_k/ε^k → 0` because
   `θ<1<ε`, and `(ε'/ε)^k → 0`) gives **`Ψ_1 = −a_0β_j`**.
2. If `a_0 = 0` then `Ψ_2 = 0` and the relation is trivial.  So `a_0 ≠ 0` and
   ```
   e_k  =  ζ·(ε')^k ,        ζ := −Ψ_2/a_0 ∈ L^* ,                  (32.3)
   ```
   for infinitely many `k`.  **An admissible relation exists only if the
   rounding-error sequence is eventually an exact geometric progression.**

3. Taking traces of (32.3):  `4p_k = ζ·Tr((ε')^k) = 2ζ(ε^k + (ε')^k)`, so
   `p_k = (ζ/2)(ε^k + (ε')^k)`.  Combining with `p_k = β_jε^k + ζ(ε')^k`:
   ```
   ( β_j − ζ/2 )·ε^k  =  −(ζ/2)·(ε')^k .                                (32.4)
   ```
   Dividing by `ε^k` and letting `k → ∞` (since `|ε'/ε| = (√a₁−1)/(√a₁+1) < 1`
   for every `a_1 ≥ 2`) gives `ζ = 2β_j`.  Then (32.4) reads
   `β_j(ε'/ε)^k = 0` for infinitely many `k`, impossible since `β_j > 0`. ∎

> **Theorem 32.3 (no admissible relation exists).**  For **every** squarefree
> `a_1` and every `a_j ≠ a_1`, **no** nontrivial `(a_0,Ψ_1,Ψ_2) ∈ L³` satisfies
> (32.2) for infinitely many `k`.  The argument is purely algebraic — **no case
> split on `a_1`, no use of condition (i) or (iv), and no asymptotic
> quantitative input.**

**Consequences.**

* The admissible relations are cut out by the **single explicit linear
  equation** `Ψ_1 = −a_0β_j`, i.e. they form a pencil in `L³` — and
  Theorem 32.3 refutes **every** point of it, uniformly and with an explicit
  algebraic contradiction.
* So the "identify the finitely many relations" half of M7.12 is **closed**: it
  is not an open enumeration problem, because the candidate space is
  explicitly described and entirely refutable.
* `Docs/31`'s `a_1 ≥ 5` case split (Prop 30.5, via `|ε'| > 1`) is
  **unnecessary** — Lemma 32.2(3) covers all `a_1` at once, including
  `a_1 = 2,3` where `|ε'| < 1`.

---

## 32.4 So what is actually left

Finiteness is unaffected — this reproduces KMN's conclusion.  **The residual
loss is now a single set:**

```
  M_θ  \  M_0'      is finite, but its cardinality and its largest element
                     are not bounded by anything in KMN.
```

`ℳ₀'` is the infinite sub-family on which the Gaussian elimination is
performed.  Every `k ∈ ℳ_θ ∩ ℳ₀'` would lie in a hyperplane, hence give an
admissible relation, hence — by Theorem 32.3 — be impossible.  So

```
  M_θ  ⊆  { k  :  the Gaussian elimination / Subspace-Theorem step fails at k } .
```

and **that set is exactly the "all but finitely many" that the Subspace
Theorem declines to identify.**  Note it is *not* the relations that are
unknown any more (Theorem 32.3 disposes of those); it is the *coverage* of
the ST that is unknown.

> **M7.12 (final form).**  Bound the exceptional tail
> ```
> E_θ  :=  { k :  the Subspace-Theorem step of Prop 3.4 fails at k } .    (32.5)
> ```
> This is precisely the classical non-effectivity of the Subspace Theorem —
> the **height** of the exceptional hyperplanes — now stripped of everything
> else.  `Docs/31` §31.4's remark that "the quantitative ST bounds the number
> `t` but not the heights/identities" is exactly right, and now it is the
> *whole* remaining problem: the identities are settled (Thm 32.3); only the
> heights are not.

**No progress is claimed on (32.5).**  But it is now the *only* thing left, it
is a single classical question, and the surrounding analysis is complete.

---

## 32.5 Status

| Item | Result |
|---|---|
| `Docs/31` claim that `ζ_j` are "fixed / finitely enumerable" | **CORRECTED** — §32.1 |
| `r = dim_L V ≤ 2`, effectively computable (fixed `q_1`) | **Thm 32.1** — NEW |
| ST applied in `L³`, not a high-dimensional space | **Thm 32.1** |
| Admissible relation ⟺ rounding error eventually geometric | **Lemma 32.2** |
| **No admissible relation exists, all `a_1`, explicit contradiction** | **Thm 32.3** — NEW |
| `Docs/30` Prop 30.5 case split on `a_1 ≥ 5` | **unnecessary** — Thm 32.3 |
| Finiteness of `{k : ‖β_jε^k‖ < θ^k}` | unchanged (KMN) |
| **Remaining: bound the ST's exceptional tail `E_θ` (32.5)** | **OPEN — the whole gap** |