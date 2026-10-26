# Paper 10 — Progressive Relative Norms Do Not Beat `2^{r-1}-1`

Version 1.0.  A negative statement, proved.  It explains *why* the tower-of-quadratics
idea cannot by itself reduce the effective exponent, and therefore what kind of
input is required to escape the conjugate explosion.

---

## 10.1 Setup

Let `K = Q(√d₁, …, √d_r)` with `F₂`-rank `r`, `[K:Q]=2^r`.  Take the Galois
group elements as sign choices `σ = (ε₁,…,ε_r) ∈ {±1}^r`, acting by
`σ(√d_i) = ε_i·√d_i`.  For a generic `S = Σ_i C_i√d_i`, every conjugate is
`σ(S) = Σ_i C_i·ε_i·√d_i`.

Order the radicands by a tower

```
Q = K₀ ⊊ K₁ ⊊ … ⊊ K_r = K,    K_i = K_{i-1}(√d_i),   [K_i : K_{i-1}] = 2.
```

Let `τ_i` denote the automorphism of `K_i` fixing `K_{i-1}` and sending
`√d_i ↦ −√d_i`.

---

## 10.2 The progressive norm is the full norm

For an intermediate form `X_i ∈ K_i`, the relative norm is

```
N_{K_i / K_{i-1}}(X_i) = X_i · τ_i(X_i)  ∈ K_{i-1}.
```

Iterating from `X = S ∈ K_r = K`:

```
N_{K_r/K_{r-1}}(S) = S·τ_r(S),
```

and this element of `K_{r-1}` has the form `Σ C(ε_1,…,ε_{r-1})√d_1·…`
summed over the `2^{r-1}` choices of the first `r-1` signs.  Repeating:
`N_{K_i/K_{i-1}}` of the product over the first `i`-sign choices equals the
product over all `2^i` sign choices on the first `i` coordinates:

```
N_{K/Q}^{(r)}(S) := Π_{σ ∈ {±1}^r} σ(S) = N_{K₁/K₀}( · ) ∘ N_{K₂/K₁}( · ) ∘ …
```

i.e. the telescoped relative-norm process **is exactly the full field norm**:

```
Π_{i=1..r} N_{K_i/K_{i-1}} = Π_{σ} σ.
```

In particular its modulus is `|N_{K/Q}(S)| ≥ 1`, and the bound derived from
it is forced to be the one used in Paper B:

```
|S| ≥ |N_{K/Q}(S)| / Π_{σ ≠ id} |σ(S)|
    ≥ 1 / (m·M^{3/2})^{2^r - 1}.
```

The exponent `2^r - 1` counts the conjugates; any purely relative-norm
rewriting has the **same** exponent, because it computes the same norm.

**Theorem 10.1.**  Any argument that isolates `|S|` from the integrality of
`N_{K/Q}(S)` alone, including iterated relative norms through any quadratic
tower, yields at most the bound `|S| ≥ (m M^{3/2})^{-(2^r-1)}`.  It cannot
recover an exponent `poly(r)`.

*Proof.* The quantity being lower-bounded by the integrality hypothesis is
exactly `1/|S| · Π_{σ≠id}|σ(S)|`.  The product over all `σ ≠ id` is determined
by the multiset of `2^r` conjugates, which is independent of the order in
which the relative norms are taken.  Each factor satisfies
`|σ(S)| ≤ m M^{3/2}`, and equality `|σ(S)| = m M^{3/2}` is attained when all
terms in the sum have the same size and sign structure (which occurs for a
nonempty open set of `C` and at least one choice of signs/embeddings, e.g.
`σ` real with each `σ(√d_i)` having the same sign as the corresponding
term).  So the inequality `Π_{σ≠id}|σ(S)| ≤ (m M^{3/2})^{2^r-1}` is sharp up
to constants, and the exponent is a structural function of the degree, not of
the tower order. ∎

---

## 10.3 Why pairing only helps by a factor of two

The all-ones sign flip (when it exists, Paper B) groups the `2^r` conjugates
into `2^{r-1}` pairs `{T, -T}`, so the norm is a product of
`2^{r-1}` squares.  The exponent of the resulting effective bound is
`2^{r-1}-1`.  This is a genuine factor-of-two saving but does **not**
collapse the count; for dependent support it fails outright.  There is no
order-2 automorphism that fixes more than two embeddings.

---

## 10.4 What additional input would beat the count

To do better than `2^{r-1}-1`, one must exploit a property of the small
conjugate **beyond** its existence.  The only such structure present is:

- **Smallness of one `L_σ`** constrains the coefficient vector `C` to lie
  near the hyperplane `{L_σ = 0}`.
- `L_σ = 0` is a `Q`-linear relation among `√d_i`, which is exactly rank-`r`
  independence disallowing.  The deviation from that hyperplane is what is
  measured.

In other words, the right object is a **linear form in `r−1` algebraic
genunates inside `K`**, and a lower bound on it is precisely the
Subspace-Theorem-style statement of degree of approximation `r−1` with a good
constant.  No relative-norm rewriting supplies that constant; it must be
*assumed or proved separately*.

---

## 10.5 What this rules in (and out)

- **Out.** Any "improved norm" via towers, squaring, or sign pairing beyond
  the factor of two.
- **Out.** Any argument based only on the integrality of the full norm.
- **In.** Simultaneous Diophantine approximation with explicit constants:
  this is the content of the open constant, and it is exactly where the work
  must go.  Concretely, the required theorem is: there exist effective
  constants `A, K` such that the approximation `degree of approximation of the
  vector `(√(d₂/d_r),…,√(d_{r-1}/d_r))` by a rational point has
  `|L| ≥ 2^{−AN^K}`.  That is equivalent to the missing constant, and is the
  content of the next milestone.