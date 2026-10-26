# Paper 12 — M4.1: A Two-Embedding Lower Bound (one direction of the dichotomy)

Version 1.0.  We prove that two *distinct* sign-embeddings of the *same*
inner form cannot both be extremely small unless the coefficient vector is
extremely large, with an explicit bound.  The content uses only the same
principle as the fixed-rank norm bound, applied to the *difference* of the
two embeddings.

---

## 12.1 Statement

**Setup.**  Let `d₁,…,d_m` be distinct squarefree, `r = log₂[K:Q]`, and
write the quantity of interest simply as a signed radical sum:

```
T_σ(C) := Σ_{i=1}^{r} C_i·σ_i·√d_i ,     σ ∈ {±1}^r,
```

i.e. `T_σ` is the value of `S` under the embedding `√d_i ↦ σ_i√d_i`.  Let
`C = (C_1,…,C_r)` be a coefficient vector with `S ≠ 0`.  When two embeddings
of the *same* vector are simultaneously small, we get a lower bound on how
small they can both be.

**Lemma M4.1a (two small embeddings cannot both beat the norm).**
Suppose `σ ≠ τ` in `{±1}^r` satisfy

```
|T_σ(C)| ≤ ε      and        |T_τ(C)| ≤ ε .
```

Let `J = { i : σ_i ≠ τ_i }`, and put `M := max(|C_i|, max_i d_i)`.
Assume the genuine-difference hypothesis that there is at least one
`i ∈ J` with `C_i ≠ 0`.  Then

```
ε ≥ ( 2 · ( 2·r·M^{3/2} )^{2^{r_J} − 1} )^{-1},
```

where `r_J = log₂[Q(√d_i : i ∈ J, C_i≠0):Q] ≤ |J| ≤ r`.  When the two
embeddings genuinely differ on the support, one can take `r_J ≤ r−1`, and the
worst case `r_J = r−1` recovers the full exponent `2^{r−1} − 1`.

---

## 12.2 Proof

**Step 1 (the difference is an algebraic integer).**  Set

```
X := T_σ(C) − T_τ(C) = Σ_{i=1}^{r} (σ_i − τ_i)·C_i·√d_i .
```

Since `σ_i − τ_i ∈ {−2,0,2}`, we have `(σ_i − τ_i)·C_i ∈ 2Z`.  Each `√d_i`
is an algebraic integer (root of the monic integer polynomial `x² − d_i`),
so `X` is an algebraic integer.

Let `J' = { i ∈ J : C_i ≠ 0 }`, which is nonempty by hypothesis, and set

```
a_i := (σ_i − τ_i)·C_i     (a_i ∈ 2Z, a_i ≠ 0 for i ∈ J').
```

Then `X = Σ_{i∈J'} a_i √d_i`.

**Step 2 (`X ≠ 0`).**  By the character sum: a `Q`-linear combination of
`√d_i` over distinct squarefree `d_i` that is zero must have all rational
coefficients zero (Paper A, Thm A.3 — the distinct radicals are
`Q`-linearly independent).  Since `a_i ≠ 0` for `i ∈ J'` and `J' ≠ ∅`, we
get `X ≠ 0`.

**Step 3 (bound all its conjugates).**  In the field
`K_{J'} := Q(√d_i : i ∈ J')`, of degree `2^{r_J}` with `r_J = r_J(C,σ,τ)`
the rank of the subset `J'`, every embedding sends each `√d_i` to `±√d_i`.
So every conjugate of `X` has the form
`Σ_{i∈J'} a_i·(±√d_i)` and therefore

```
|σ(X)| ≤ Σ_{i∈J'} |a_i|·√d_i ≤ 2·Σ_{i∈J'} |C_i|·√d_i
       ≤ 2·r·M·√M = 2·r·M^{3/2}.
```

**Step 4 (norm).**  `N_{K_{J'}/Q}(X) = Π_{σ} σ(X)` is a non-zero integer
(because `X` is a non-zero algebraic integer, Step 2), so its absolute value
is `≥ 1`.  There are `2^{r_J}` factors, one of which is `|X|` itself; the
remaining `2^{r_J} − 1` factors are each `≤ 2·r·M^{3/2}`.  Hence

```
1 ≤ |X| · ( 2·r·M^{3/2} )^{2^{r_J} − 1},
```

so

```
|X| ≥ ( 2·r·M^{3/2} )^{-(2^{r_J} − 1)}.
```

**Step 5 (convert to `ε`).**  But

```
|X| = |T_σ(C) − T_τ(C)| ≤ |T_σ| + |T_τ| ≤ 2ε.
```

Combine:

```
2ε ≥ ( 2·r·M^{3/2} )^{-(2^{r_J} − 1)},
```

```
ε ≥ ( 2·( 2·r·M^{3/2} )^{2^{r_J} − 1} )^{-1}.
```

∎∎

---

## 12.3 What this says (and what it does not)

*If two different embeddings of the same `C` are simultaneously bounded by
`ε`, then `ε` is bounded below by the norm-bound threshold of the
sub-form supported on their difference set `J'`.*  Taking `r_J = r−1` gives
the absolute statement

```
two simultaneous embeddings both < ε  ⇒  ε ≥ (2·(2r M^{3/2})^{2^{r-1}−1})^{-1},
```

i.e. *you cannot have two embeddings simultaneously below the rank-`r−1`
norm threshold of a linear form in `r−1` variables.*  This rigorously
delimits the "thin simultaneous regime" that M2 could not see: the first
place two embeddings can coincide is exactly at a rank-`r−1` norm-form lower
bound.  For `r = 2`, `r_J = 1`, and the bound becomes the Pell/`λ` scale,
consistent with the effective rank-2 treatment in Paper D.

**What it does not say.**  This is a *lower* bound on the threshold for two
small embeddings; it does not assert that the simultaneous regime is actually
populated by interior points, and it does not give the polynomial-in-`N`
scale the uniform conjecture needs.  It is the *sharp local obstruction*:
below this threshold, "two small embeddings" cannot happen, so the only way
the Subspace-sharp exponent `r−1` can be saturated is by a single small
embedding — the generic picture measured in Paper M2/11.

---

## 12.4 Consequence for M4's dual branch

The dichotomy is now quantitative:

* **One small embedding**, all others bounded away from the threshold →
  apply the norm bound, already effective.
* **Two small embeddings** → they are both at or above the rank-`r_J`
  norm threshold; their difference yields a new integer radical form whose
  smallness is explained by its own norm.

The second branch reduces, formally, to a **rank-`r_J` instance of the same
problem**.  M4.2 should test computationally whether, in a putative
two-small-embedding point, the extracted difference form can itself be
certified small (which, by the present Lemma, forces it to sit at the
rank-`r_J` lower-bound scale).  That is the next observable.