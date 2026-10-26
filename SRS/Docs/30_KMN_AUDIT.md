# Paper 30 — M7.11a: audit of the KMN non-degeneracy and pseudo-Pisot conditions

Version 1.0.  Audits, against the **verbatim definitions** of
Kulkarni–Mavraki–Nguyen (arXiv:1511.08525v1, Trans. AMS 2018), whether our pair
`(β_j, ε)` satisfies the hypotheses of their Theorem 1.4, and what their
theorem actually delivers.

**Outcome.**

* **M7.11a is resolved — and it is *vacuous*.**  Their Definition 1.2 makes
  non-degeneracy vacuous at `k = 1`.
* **The pseudo-Pisot condition (i) fails *explicitly*, by an archimedean size
  comparison**, and so does (iv) whenever `a_1 ≥ 5`.
* **Consequence:** `‖β_jε^k‖ < θ^k` has only **finitely many** solutions for
  **every** `θ ∈ (0,1)`, uniformly, with **no quantitative constant entering**.
* **But this is finiteness, not effectivity** — and the gap is now precisely
  measurable: KMN delivers `η_k ≥ e^{-o(k)}`; `Docs/28` Thm 28.3 needs
  `η_k ≥ k^{-c} = e^{-O(log k)}`.  **One step remains, and it is now measured
  against a published non-effective bound rather than against nothing.**

---

## 30.1 The verbatim hypotheses

**Definition 1.2 (non-degenerate), verbatim:**
> "The tuple of non-zero algebraic numbers `(α_1,…,α_k)` is called
> **non-degenerate** if `α_i/α_j ∉ μ` for `1 ≤ i ≠ j ≤ k`."

**Definition 1.3 (pseudo-Pisot), verbatim:**
> "Let `(β_1,…,β_k)` be a tuple of distinct non-zero algebraic numbers. Write
> ```
> B := { β ∈ Q̄* \ {β_1,…,β_k}  :  β = σ(β_i) for some σ ∈ G_Q, 1 ≤ i ≤ k } .
> ```
> Then the tuple `(β_1,…,β_k)` is called **pseudo-Pisot** if
> `Σ_i β_i + Σ_{β∈B} β ∈ Z` and `|β| < 1` for every `β ∈ B`."

**Theorem 1.4, verbatim (the four properties).**  Let `(α_1,…,α_k)` be
non-degenerate with `|α_i| ≥ 1`, `f` sublinear, and let `M` (the set of
`(n,q_1,…,q_k)` with `‖Σ q_iα_i^n‖_Z < θ^n` and `max h(q_i) < f(n)`) be
infinite.  Then for all but finitely many elements of `M`:

```
(i)   the tuple (q_1α_1^n, …, q_kα_k^n)  is pseudo-Pisot;
(ii)  α_i  is an algebraic integer  for all i;
(iii) σ(q_iα_i^n) = q_jα_j^n   ⟺   σ(α_i)/α_j ∈ μ;
(iv)  for σ ∈ G_Q, 1≤i,j≤k:  if σ(α_i)/α_j ∉ μ  then  |σ(α_i)| < 1.
```

And their own remark on why these are the *right* conditions, verbatim:
> "In fact `q_1α_1^n + ⋯ + q_kα_k^n` is a **subsum of an integer-valued trace**,
> and, according to properties (iii) and (iv), any term that does not appear in
> this subsum must have the form `σ(q_iα_i^n)` with `|σ(α_i)| < 1`."

---

## 30.2 M7.11a resolved: non-degeneracy is vacuous at `k = 1`

Our problem is `k = 1`, `α_1 = ε`, `q_1 = β_j` (fixed, so `h(q_1) = O(1) < f(k)`
for any sublinear `f`; and `|ε| = 1+√a_1 > 1`).

> **Corollary 30.1.**  For `k = 1` the non-degeneracy condition of Definition 1.2
> is **vacuous**: it quantifies over pairs `i ≠ j`, of which there are none.
> **M7.11a requires no work and no verification — it is automatic.**

This is the first of the user's two proposed sublemmas, and it is closed.

---

## 30.3 Condition (i) fails *explicitly* — the pseudo-Pisot audit

The tuple in (i) is `(β_jε^k)`, a single element.  So `B` is the set of its
Galois conjugates in `K = Q(√a_1,√a_j)` other than itself.  By `Docs/28`
Thm 28.1(4) the conjugates of `β_jε^k` are `{ ±β_jε^k, ±β_jε^{-k} }`, hence

```
B = { −β_jε^k ,  β_jε^{-k}·(±1) } .
```

Pseudo-Pisotness therefore requires **`|β| < 1` for every `β ∈ B`**, i.e. in
particular

```
|−β_jε^k|  <  1 ,     i.e.   β_jε^k < 1 .
```

But `β_j > 0` and `ε > 1`, so `β_jε^k → ∞`.  Therefore:

> **Proposition 30.2 (explicit failure of (i)).**  For all sufficiently large
> `k`, the tuple `(β_jε^k)` is **not** pseudo-Pisot, because `B` contains
> `−β_jε^k` with `|−β_jε^k| = β_jε^k > 1`.  This is a **plain archimedean size
> comparison**, with no Subspace-Theorem content.

(The integrality half of Def 1.3 *does* hold: `Σβ_1 + Σ_{β∈B}β` is the trace
of `β_jε^k`, which is `0 ∈ Z`.  It is the `|β|<1` clause that fails.)

**Consequence.**  Theorem 1.4 asserts (i) for all but finitely many elements of
an infinite `M`.  Proposition 30.2 says (i) fails for all large `k`.  Hence:

> **Theorem 30.3 (finiteness, uniform in `θ`).**  Let `a_1 ≠ a_j` be distinct
> squarefree integers, `ε = 1+√a_1`, `β_j = √(a_j/a_1)/2`.  Then for **every**
> `θ ∈ (0,1)`,
> ```
> { k ∈ ℕ  :  ‖β_jε^k‖ < θ^k }      is finite.                       (30.1)
> ```
> No constant depending on the data appears in the argument.

Equivalently, with `η_k := ‖β_jε^k‖`:

> **Corollary 30.4 (a new quantitative statement).**  `η_k ≥ θ^k` for all
> sufficiently large `k`, for **every** `θ < 1`.  Hence
> ```
> liminf_k η_k^{1/k} ≥ 1 ,      i.e.      η_k = e^{−o(k)} ,
> ```
> `η_k` decays **sub-exponentially**.

---

## 30.4 The independent second obstruction: (iv) fails for `a_1 ≥ 5`

Condition (iv) requires `|σ(α_1)| < 1` whenever `σ(α_1)/α_1 ∉ μ`.  For
`α_1 = ε = 1+√a_1` the nontrivial conjugate is `1−√a_1`, and

```
|1−√a_1| = √a_1 − 1  <  1   ⟺   a_1 < 4   ⟺   a_1 ∈ {2,3} .
```

> **Proposition 30.5.**  For `a_1 ≥ 5`, `ε = 1+√a_1` is **not a Pisot
> number** and condition (iv) of Theorem 1.4 **fails**.  This gives a second,
> independent route to (30.1), and it coincides with the classical fact that
> `1+√a_1` is Pisot only for `a_1 = 2,3`.

So: for `a_1 = 2,3` the finiteness is forced by (i); for `a_1 ≥ 5` it is forced
by (iv).  Either way (30.1) holds.  **This is a clean dichotomy and it is a
useful sanity check on the whole reduction.**

---

## 30.5 Why KMN and not Corvaja–Zannier — and a warning

KMN state, verbatim: *"The case `k = 1` and `q_1` is essentially a rational
integer was obtained by Corvaja and Zannier."*  Our `q_1 = β_j` is **not**
rational, so the Corvaja–Zannier Theorem 1.1 does **not** apply; **KMN is
required.**  `Docs/29` §29.2 should be read accordingly.

This distinction is not cosmetic.  If `q_1 = 1`, i.e. if we were studying
`‖ε^k‖`, then for `a_1 = 2,3` we would have `ε` Pisot and
`‖ε^k‖ ≈ |1−√a_1|^k`, so `‖ε^k‖ < θ^k` for **every** `k` once
`θ > √a_1−1`.  Infinitely many exceptions.

> **So the condition `β_j ∉ Q` is exactly what removes the problem from the
> degenerate regime**, and it is why the `q_1 ∈ K*` extension in KMN — not
> Corvaja–Zannier — is the operative theorem.  `Docs/29`'s citation of "the
> `k=1` case is CZ's" is correct only for rational `q_1`; corrected here.

---

## 30.6 What is still missing, measured precisely

**Finiteness is not effectivity.**  `Docs/28` Thm 28.3 requires
`η_k ≥ k^{-c}` eventually (equivalently `η_k ≥ e^{-O(log k)}`), not merely that
each `θ^k` is exceeded finitely often.  The three statements now sit at:

```
  elementary (Liouville + norm, Docs/28 Thm 28.2) :  eta_k >= e^{−3k log eps} = e^{−Θ(k)}
  KMN 2021/2018, non-effective (Thm 30.3)          :  eta_k >= e^{−o(k)}
  REQUIRED by Docs/28 Thm 28.3                      :  eta_k >= e^{−O(log k)} = k^{−c}
```

`e^{-o(k)}` does **not** imply `k^{-c}`: a rate `e^{-√k}` is `e^{-o(k)}` yet is
`≪ k^{-c}`.  And KMN's proof rests on the Subspace Theorem, so it inherits the
ineffectivity — no explicit `K_ST` or `c` is produced.

> **Open (M7.11b′, restated and now sharply bounded on both sides).**  Improve
> ```
>      eta_k = ‖β_j ε^k‖  ≥  e^{−o(k)}          (KMN, non-effective)
> ```
> to
> ```
>      eta_k  ≥  e^{−C log k}   for an explicit C = C(β_j, ε).         (M7.11b′)
> ```
> Note `C` may depend on `β_j, ε` — i.e. on `a_1, a_j` — **without any need to
> control that dependence in `N`**, because `Docs/28` Thm 28.3 absorbs it into
> `K_0 = O(n c/ log ε)`, which is already independent of `N` (`Docs/29` §29.5).

**This is the sharpest statement of the remaining gap in the project.**  It is a
*quantitative* improvement of a *published non-effective* bound, by one
exponential step — not an appeal to an unknown theorem.

---

## 30.7 Status

| Item | Result |
|---|---|
| KMN Definition 1.2 (non-degeneracy), verbatim | obtained |
| KMN Definition 1.3 (pseudo-Pisot), verbatim | obtained |
| KMN Theorem 1.4 (i)–(iv), verbatim | obtained |
| **M7.11a** | **RESOLVED — vacuous at `k=1`** (Cor 30.1) |
| (i) fails explicitly, archimedean | **Prop 30.2** |
| (iv) fails for `a_1 ≥ 5`; `ε` Pisot iff `a_1∈{2,3}` | **Prop 30.5** |
| `‖β_jε^k‖ < θ^k` finite for **every** `θ<1`, uniform | **Thm 30.3** |
| `η_k = e^{-o(k)}` | **Cor 30.4** — new |
| CZ vs KMN: `q_1` must be in `K*` | §30.5; `Docs/29` citation corrected |
| **M7.11b′: `e^{-o(k)} → e^{-O(log k)}` with explicit `C`** | **OPEN — the whole remaining gap** |