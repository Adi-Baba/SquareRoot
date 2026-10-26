# Paper 37 — Tracks A and B: the literature audit, and the product-level attack

Version 1.0.  Executes the two-track plan.  Norm/invariant/height
experimentation is **frozen** from here on.

* **Track A (literature).**  Searched for an *effective* theorem on
  `‖βε^k‖` with `β` outside the quadratic field of `ε`.  Nothing located
  supplies the missing effectivity.  §37.1.
* **Track B (product level).**  Two new theorems and one exact new fact, all
  certified — including the first genuine use of the *shared* Pell pair, which
  the factorwise route of `Docs/29` discards.  But the attack **does not close
  the gap**, for a reason that is now precise.  §37.2–37.5.

---

## 37.1 Track A — literature audit

The question was re-aimed per plan: not "better norm bounds" but *dynamics /
recurrence / distribution mod 1*, and specifically an **effective bound on the
largest exceptional `k`** (a bare finiteness theorem is worthless now).

**Closest hit — the Loxton–Van der Poorten family.**
Noubissie, *Quantitative growth of linear recurrences*, **arXiv:2504.09519**
(Apr 2025), proves an **effective** upper bound for the **number** of solutions
of

```
| u_n |  <  ( max_i |alpha_i| )^{ n (1 - eps) }
```

for a non-degenerate integer linear recurrence `u_n = Σ P_i(n) alpha_i^n`,
completing a programme of **Schmidt** and **Van der Poorten** whose qualitative
(non-effective) forms were proved by **Fuchs–Heintze (2021)** and, independently,
**Karimov et al. (2023)**.

Our condition *is* of this shape: `θ₀^k = ε^{−k/(3n(n−1))}` is
`ρ^{n(1−ε)}` with `ε = 1 − 1/(3n(n−1))`.  So this is the right family.  It
nonetheless does **not** supply what M7.12 needs, for two reasons:

1. **It bounds the COUNT, not the largest index.**  A count bound is what
   `Docs/33` Thm 33.3–33.4 already delivers (`|ℳ_θ| ≤ n₀ + t`), and a count
   bound does not bound `k*`.
2. **It is for integer LRS.**  `u_n` must be an integral recurrence.
   `x_k = β_jε^k` is **not** one: `β_jε^k = Q_k√a_j/2 + P_k√(a₁a_j)/(2a₁)`
   (§37.2), a genuine two-term form.  The integer recurrences here are the Pell
   sequences `P_k, Q_k`, which satisfy no such smallness condition.

**Other families checked and set aside.**  *Multiplicative dependence in linear
recurrence sequences* (arXiv:2501.17365) — bounds tuples whose **values** are
multiplicatively dependent, a different relation.  *Effective results on the
Skolem Problem* (arXiv:1505.07147) and *On large zeros of LRS* — both about
**zero terms**, i.e. `|u_n| = 0`, the opposite extreme from `‖·‖ ≪ 1`.

> **Track A verdict.**  No located theorem gives an effective bound on the
> largest `k` with `‖βε^k‖ < θ^k` for `β ∉ Q(ε)`.  The nearest neighbour counts
> solutions instead of locating the last one.  I am recording this as a search
> result, **not** as a claim that no such theorem exists.

## 37.2 The object, stated correctly

A preliminary correction, because everything below depends on it.  With
`ε^k = P_k + Q_k√a₁` (index = exponent, per `Docs/34` §34.2) and
`β_j = √(a_j/a₁)/2`:

```
beta_j eps^k   =   Q_k sqrt(a_j) / 2   +   P_k sqrt(a_1 a_j) / (2 a_1) .        (37.1)
```

Verified exactly at `k = 5, 10, 20` to 20 digits.  So `‖β_jε^k‖` is a genuine
**two-term form**, not `‖Q_k√a_j‖`.

> **Correction to `Docs/35` §35.1.**  The corrected identity there reads
> `Q_k√a_j = β_j(ε^k − ε^{−k})`.  It must be `(ε')^k`, not `ε^{−k}`: since
> `ε' = −ε^{−1}` we have `ε'^k = (−1)^k ε^{−k}`, and the two differ at **odd**
> `k`.  Confirmed: the identity as printed in `Docs/35` matches at `k = 2, 10`
> and **fails at `k = 1, 3, 5`**.  The *conclusion* of `Docs/27` — that
> `‖β_jε^k‖ = ‖Q_k√a_j‖ + O(ε^{−k})` — is unaffected, since the two differ by
> exactly `β_jε^{−k}` (measured: `1.354e−8` at `k = 20`, matching
> `β_jε^{−20} = 1.354e−8`).

## 37.3 Fact 37.2 — the `a₁` coordinate is exact

From `ε^k = P_k + Q_k√a₁` and `ε'^k = P_k − Q_k√a₁`:

```
Q_k sqrt(a_1)   =   P_k − eps'^k        =>        || Q_k sqrt(a_1) ||   =   |eps'|^k   =   eps^(-k) .
```

Certified for `k = 1…79` (`s1_exact`: the nearest integer to `Q_k√a₁` is `P_k`,
with the sign of `ε'^k`, alternating in parity).

> This is **not** a generic coordinate.  The `a₁` factor in the Markov product
> is *exactly* `q^{−1}` up to `2√a₁`.  It therefore consumes the entire `q^{−1}`
> in the EHS budget `∏_i ‖q√a_i‖ < q^{−1−δ}`, leaving only
> ```
> prod_{i=2..n} || q sqrt(a_i) ||   <~   q^{-delta} .                           (37.2)
> ```

That is a genuine sharpening of the reduction and it sharpens the target to a
*small* exponent.

## 37.4 Theorem 37.1 — the pair bound (first use of the shared Pell pair)

The factorwise route never uses the fact that every coordinate shares the same
`(P_k, Q_k)`.  Multiplying (37.1)-based equations for two coordinates by
`√a_{j'}` and `√a_j` respectively makes the `P_k` and `Q_k` terms **cancel**:

```
m_{j'} sqrt(a_j) − m_j sqrt(a_{j'})   =   e_j sqrt(a_{j'}) − e_{j'} sqrt(a_j) ,
e_j := beta_j eps^k − m_j ,   m_j := nint(beta_j eps^k) .                       (37.3)
```

Since `(m_{j'}√a_j − m_j√a_{j'})(m_{j'}√a_j + m_j√a_{j'}) = m_{j'}²a_j − m_j²a_{j'}`
is a **nonzero integer** (zero would make `a_j/a_{j'}` a rational square), we get

> **Theorem 37.1.**  For `j ≠ j' ≥ 2`,
> ```
> |e_j| sqrt(a_j') + |e_j'| sqrt(a_j)   >=   1 / ( m_{j'} sqrt(a_j) + m_j sqrt(a_j') )
>                                          ~   eps^{-k} / ( 2 beta_{j'} sqrt(a_j) ) .
> ```

**Certified exactly.**  Clearing the radicals turns it into a sign test on a
six-surd linear form in `{1, √a_j, √a_{j'}, √(a₁a_j), √(a₁a_{j'}), √(a_ja_{j'})}`
— `pair_bound_holds`.  Result: **474/474** certified over 6 radicand pairs ×
`k = 1…79`; `I_{j,j'} ≠ 0` in 354/354.

This is the exponent-**1** bound, and it is genuinely new relative to `Docs/36`
(it lives in the degree-2 field `Q(√a_j, √(a_{j'}))` *after* the cancellation).

## 37.5 Theorem 37.3 — the product bound, and why it still does not close

Order `s_2 ≤ s_3 ≤ … ≤ s_n`, `s_i := |e_i|`.  Theorem 37.1 with the pair
`(2, j')` and `s_2 ≤ s_{j'}` gives `s_{j'} ≳ ε^{−k}` for every `j' ≥ 3`, while
`s_2` is stuck at the single-index exponent `3` (`Docs/35` Thm 35.1, which
`Docs/36` Thm 36.2 proves is structural).  Hence

> **Theorem 37.3.**  `prod_{i=2..n} s_i  >=  C_n · q^{−(n+1)}`, with `C_n > 0`
> explicit.  Verified numerically for `n = 3, 4, 5` (`a₁=2`), e.g.
> `C_3 = 3.27e−16`, `C_4 = 5.95e−30`.

The exponent is `3 + (n−2)·1 = n+1`: **one** coordinate stuck at `3`, the rest
pushed from `3` to `1`.

**Now combine with the target.**  By (37.2) the requirement is
`∏_{i=2..n} s_i ≲ q^{−δ}`, so Theorem 37.3 contradicts it **iff `n+1 < δ`**:

| `n` | bound exponent | closes when |
|---|---|---|
| 3 | 4 | `δ > 4` |
| 4 | 5 | `δ > 5` |
| 5 | 6 | `δ > 6` |

> **Track B verdict — a negative result, stated plainly.**  For the small `δ`
> that the EHS exponent is used with, Theorem 37.3 gives **no contradiction**.
> Multi-coordinate structure buys exactly this much and no more: it moves `n−2`
> coordinates from exponent `3` to exponent `1`, and leaves one coordinate at
> `3`.  Since the requirement is an exponent `δ < 1`, **`Docs/36` Thm 36.2's wall
> survives the product attack unchanged.**  Track B does not rescue M7.12.

**But it does relocate the question decisively.**  Theorem 37.3 is
non-vacuous precisely when `δ > n+1`.  So the live question is no longer
Diophantine — it is:

> **Does EHS admit `δ > n+1` in its hypothesis?**

That is checkable against the EHS statement, and it is the only remaining lever
Track B found.  I have **not** checked it here, and I am not claiming either
answer.

## 37.6 Corrections and bugs found in this pass

* **`Docs/35` §35.1:** `ε^{−k}` should be `(ε')^k`; they differ at odd `k`
  (§37.2).  `Docs/35` corrected in place.
* **Engine bug, found by the exact certificate and fixed.**  The local helper
  computing `sign(L_{j'})` closed over the radicand `a_j` instead of taking
  `a_{j'}` as a parameter, so `j'`'s sign was decided in the wrong field.  This
  made the certificate report **68 spurious failures out of 295**.  Two earlier
  versions of `pair_bound_holds` were also wrong — the first discarded the
  `|e_{j'}|` term (a weakening that fails where the truth holds), and the
  assembly of the six-surd form was checked term-by-term against
  `squarefree.decompose_and_aggregate` to confirm each coefficient.  Recorded
  because the failure mode — a silently wrong radicand in a sign test — is
  exactly the class of defect that a green test suite can hide.

## 37.7 State of M7.12

| route | outcome |
|---|---|
| heights / Liouville | refuted (`Docs/34` Thm 34.1) |
| Ridout / Subspace-Theorem power laws | strictly worse (`Docs/35` §35.6) |
| norm / lattice / product-of-conjugates | structurally exhausted (`Docs/36` Thm 36.2) |
| invariant layer | degenerate (`Docs/36` §36.3) |
| conjugate amplification | measured flat (`Docs/36` §36.6) |
| **Track A (literature)** | **no effective largest-`k` theorem located** (§37.1) |
| **Track B (product)** | **new Thm 37.1, Fact 37.2, Thm 37.3; does not close, wall unchanged** (§37.5) |
| remaining lever | whether EHS admits `δ > n+1` — a question about EHS |

**Recommendation, as planned.**  Norm/invariant/height experimentation is frozen.
The one cheap, decisive check left is the EHS parameter range; if `δ ≤ n+1` is
forced, M7.12 should be **parked** with `Docs/36` Thm 36.2 as its stopping
theorem and `Docs/37` as the record of why the product route cannot help either.