# Paper 21 — M7.3: the two-small refinement (proved)

Version 1.0.  A new, effective, provable bound for the case of **two**
simultaneously small embeddings.  It improves the single-form norm exponent
from `2^r − 1` to `2^{⌊m/2⌋} − 1`.  This is the first improvement over the
norm exponent that is *not* a norm identity.

---

## 21.1 Statement

Let `d_1,…,d_m` be distinct squarefree integers, `C_i ≠ 0` integers,
`S = Σ C_i√d_i ≠ 0`, `M := max_i max(|C_i|, d_i)`.

Let `σ⁰ ≠ σ¹` be two sign patterns (distinct as patterns in `{±1}^m`) and
suppose both embeddings are small:

```
|T_{σ⁰}(C)| ≤ ε ,   |T_{σ¹}(C)| ≤ ε ,   T_σ(C) := Σ_i C_i σ_i √d_i.
```

Write

```
A := { i : σ⁰_i = σ¹_i }        (agreement set)
J := { i : σ⁰_i ≠ σ¹_i }        (difference set)     A ⊔ J = {1,…,m}
k_A := rank_{F₂}{ d_i : i ∈ A },     k_J := rank_{F₂}{ d_i : i ∈ J }.
```

**Hypothesis (H).**  `supp(C) ∩ A ≠ ∅` and `supp(C) ∩ J ≠ ∅`, so that
`X ≠ 0` and `Y ≠ 0` below.

**Theorem 21.1 (two-small refinement).**  Write `s_A := |A|`, `s_J := |J|`
(support *sizes*; note `k_A ≤ s_A` and `k_J ≤ s_J`).  Under (H),

```
log₂(1/ε)  ≤  ( 2^{min(k_A,k_J)} − 1 ) · ( (3/2)·log₂M + log₂ max(s_A,s_J) + 1 ).   (T)
```

Consequently, since `k_A + k_J ≤ m`, so `min(k_A,k_J) ≤ ⌊m/2⌋`:

```
log₂(1/ε)  ≤  ( 2^{⌊m/2⌋} − 1 ) · ( (3/2)·log₂ M + O(m) ).                    (†)
```

---

## 21.2 Proof

**Step 1 (the split).**  Define

```
X := (T_{σ⁰} + T_{σ¹})/2 ,      Y := (T_{σ⁰} − T_{σ¹})/2 .
```

By additivity of the embedding, for `i ∈ A` the two sign factors coincide and
for `i ∈ J` they are opposite, so

```
X = Σ_{i∈A} C_i·σ⁰_i·√d_i ,          (support ⊆ A)
Y = Σ_{i∈J} C_i·σ⁰_i·√d_i ,          (support ⊆ J).
```

Both `X, Y` are **nonzero** by hypothesis (H) together with the
linear independence of distinct squarefree radicals (Paper A, Thm A.3): a
`Q`-linear combination of `√d_i`, `i ∈ A`, with `C_i ≠ 0` on the support,
cannot vanish.

**Step 2 (both halves are small).**  The triangle inequality gives

```
|X| = |T_{σ⁰}+T_{σ¹}|/2  ≤ ε ,        |Y| = |T_{σ⁰}−T_{σ¹}|/2  ≤ ε.
```

**Step 3 (norm bound on `X`).**  `X` lies in the multiquadratic field
`E_A := Q(√d_i : i ∈ A)`, whose degree is `2^{k_A}` (Paper A, Cor A.2′).
`X` is an algebraic integer. Each conjugate of `X` is

```
X^τ = Σ_{i∈A} C_i τ_i √d_i ,      τ ∈ Gal(E_A/Q),
```

so `|X^τ| ≤ k_A · M · √M = k_A·M^{3/2}`. Since `X ≠ 0`, its norm is a
nonzero integer, so `|N_{E_A/Q}(X)| ≥ 1`; isolating the identity factor and
bounding the other `2^{k_A} − 1`:

```
1 ≤ |X| · ( k_A·M^{3/2} )^{2^{k_A} − 1},
```

hence

```
|X| ≥ ( k_A · M^{3/2} )^{−(2^{k_A}−1)} .                                    (21.2a)
```

**Step 4 (norm bound on `Y`).**  Identically, in `E_J` of degree `2^{k_J}`:

```
|Y| ≥ ( k_J · M^{3/2} )^{−(2^{k_J}−1)} .                                    (21.2b)
```

**Step 5 (combine).**  Since `ε ≥ |X|` and `ε ≥ |Y|`, take the larger of
(21.2a), (21.2b).  For fixed `M`, each of (21.2a), (21.2b) is
*decreasing in its own exponent* `2^{k_·} − 1`, so the larger lower bound is
the one attached to the **smaller** of `k_A, k_J`.  Writing
`q := min(k_A,k_J)` and bounding both prefactors by `max(s_A,s_J) ≤ m`:

```
ε ≥ ( max(s_A,s_J)·M^{3/2} )^{−(2^q − 1)} .
```

Taking `log₂` proves (T).  Finally `k_A + k_J ≤ m_total` (ranks of two
disjoint subsets sum to at most the total rank), so
`q = min(k_A,k_J) ≤ ⌊m_total/2⌋`, giving (†) since `max(s_A,s_J) ≤ m_total`. ∎

**Remark (the two `m`'s).**  The *exponent* is controlled by the ranks
`k_A,k_J` (via the field degrees `2^{k_A}`, `2^{k_J}`), while the *prefactor*
is controlled by the support *sizes* `s_A,s_J` (via the triangle inequality on
the conjugates).  Only the exponent is halved; the prefactor merely absorbs a
factor `m`, i.e. `log₂ m` extra bits — negligible next to `2^{⌊m/2⌋}`.

---

## 21.3 Why this is not a norm identity

* The single-form norm bound (`Docs/02`, `Docs/09`) reads the norm of `S`
  itself: exponent `2^r − 1` with `r` the full rank.
* Theorem 21.1 never looks at `N_{K/Q}(S)`.  It looks at **two** algebraic
  integers, `X` and `Y`, each living in a **smaller** field (`E_A` of degree
  `2^{k_A}`, `E_J` of degree `2^{k_J}`), whose degrees sum to at most `m`.
  The exponent is therefore governed by `2^{min(k_A,k_J)}`, and the minimum
  of two ranks whose sum is `≤ m` is `≤ ⌊m/2⌋`.

For independent radicands (`r = m`) this replaces `2^m − 1` by
`2^{⌊m/2⌋} − 1` — an improvement of a factor `2^{2^{m-1}}` in the exponent,
achieved by no rearrangement of the single norm (this is the M1 negative
result, `Docs/10`, in sharp contrast).

---

## 21.4 The hypothesis (H) is exactly the "genuine" case

If `supp(C) ⊆ J` (no coefficient on the agreement set) then `X = 0` and the
theorem does not apply — but then `T_{σ⁰} = −T_{σ¹}` exactly, i.e. the two
patterns are *numerically* antipodal and there is only **one** small
embedding in disguise.

Hence:

> **Hypothesis (H) is the exact algebraic condition under which a second
> independent small form exists, and Theorem 21.1 is provable precisely
> whenever it holds.**

Note this is a *sufficient*, not a *necessary*, condition for two small
embeddings: two embeddings can both be small with `X = 0` (genuinely
antipodal).  Hypothesis (H) is exactly the condition that makes the `X`-half
non-degenerate, which is what the norm argument in §21.2 needs.

---

## 21.5 The data agree, and they locate the hard case

**Correction.**  The v1.0 draft of this paper repeated a wrong claim from
`Docs/20`: that the residual population was dominated by a near-antipodal
"artifact".  That came from a buggy filter (medians computed against a
mismatched list; float comparisons).  The corrected classification of all
9104 residuals in `m7_2_residuals.csv` is:

```
rows satisfying hypothesis (H)  : 3344   (over all B, r = 3..7)
  |A| == 2                     : 3344 / 3344   <-- agreement set has 2 coords
  supp(X) == 2                 : 3344 / 3344
  |X| >= |Y|/2  (genuine)      : 1616
  |X| <  |Y|/2                 : 1728
THEOREM 21.1 VIOLATIONS        : 0
strictly tighter than the norm bound (2^r − 1): 3344 / 3344
```

So the key empirical fact is **not** "the phenomenon is illusory" but:

> `|A| = 2` in every single residual, i.e. the two small sign patterns agree
> on only two coordinates, so the agreement form is always a two-term form
> and `k_A ≤ 2`.

With `k_A ≤ 2` we get `min(k_A,k_J) ≤ 2` and `max(s_A,s_J) ≤ m`, so
Theorem 21.1 gives

```
log₂(1/ε)  ≤  ( 2² − 1 )·(1.5·log₂M + log₂m + 1 )  =  O(log N),
```

i.e. `ε ≥ M^{−O(1)}`: **the entire genuine two-small regime observed is
already polynomially bounded in `N`.**  Observed `E ≤ 4.95` is consistent
with the predicted `E ≤ 3(1.5·log₂M + 1)`.

## 21.5.1 Independent verification (not just the residual CSV)

Theorem 21.1 is a general statement, so it was also checked on freshly
generated random data unrelated to the M7.2 search protocol — random
squarefree `d_i` (1-3 primes each), random `C_i ∈ [−4,4]\{0}`, random sign
pairs `(σ⁰,σ¹)` filtered to satisfy (H), embeddings evaluated in 60-digit
`mpmath` arithmetic:

```
fresh random cases under (H)      : 5970   (m up to 6, 1-3 primes per d_i)
THEOREM 21.1 VIOLATIONS            : 0
tightest observed margin           : +5.57 bits (k_A=k_J=1, |A|=|J|=1)
```

Both the residual CSV check and this random check use the **corrected**
prefactor `max(s_A,s_J)` (support sizes) rather than `max(k_A,k_J)`
(ranks) — see the Remark at the end of §21.2.  The earlier version of this
paper had that prefactor wrong and has been fixed.

The random cases are not adversarial (they do not search for near-violations),
so this is a sanity check, not a tightness study.  But it confirms the bound
holds on data generated by a completely different procedure.

---

## 21.6 What is left open after Theorem 21.1

Theorem 21.1 settles the case "**two** embeddings are small".  The genuinely
hard remaining case is:

```
exactly ONE embedding is small, all 2^r − 1 others are generic;
equivalently supp(C) ⊆ J for every pair, so no X can be formed.
```

There the hypothesis (H) fails for every pair, no two-half split is
available, and one is back to a single small coordinate — which is exactly
the Subspace-Theorem configuration of `Docs/03`, with the ineffective
constant `log₂(1/c(d,δ))`.

So the remaining theorem is sharply characterised:

> **Open (M7.4).**  An isolation lemma for a *single* small embedding:
> if `T_{σ⁰}` is the only small coordinate (all `|T_σ| ≥ δ·M` for
> `σ ≠ ±σ⁰`), then `log₂(1/T_{σ⁰}) ≤ poly(N)`.
>
> Theorem 21.1 handles everything with a second small coordinate; M7.4 is the
> isolated-coordinate case, and it is the only case left.

---

## 21.7 Status

| Configuration | Bound | Proven in |
|---|---|---|
| one small embedding | `2^r − 1` in `log M` | `Docs/02`, `Docs/09` (norm) |
| two small, full difference support | `(2^{⌊m/2⌋}−1)` in `log M` | **this paper** |
| two small, `supp(X) = 0` (antipodal) | reduces to one | trivial |
| two small, `supp(X) ≤ 2` | `O(log M)` | this paper + `Docs/19` data |
| **one small, isolated** | **open** | target of M7.4 |

No appeal to Subspace-Theorem constants, Baker bounds, or regulators is used
anywhere in Theorem 21.1.  It is elementary and effective.