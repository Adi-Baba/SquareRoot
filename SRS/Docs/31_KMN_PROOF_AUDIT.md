# Paper 31 — M7.12: does KMN's proof contain the needed rate?

Version 1.0.  Direct answer to the question: **no** — and the loss is located to
a *single sentence*.  A structural rigidity observation is also obtained which
gives genuine leverage on exactly that step.

All quotations are verbatim from arXiv:1511.08525v1.

---

## 31.1 Where the Subspace Theorem is used at all

The proof of KMN Theorem 1.4 invokes the Subspace Theorem (their Theorem 2.1)
in exactly **two** places:

* **Proposition 2.2** (finiteness of exact zeros) — by minimality in `k`.
* **Proposition 3.4** (the key technical result).

Everything else is Northcott, Gaussian elimination, or elementary.  Crucially:

> **Proposition 2.3, verbatim:** "When `k = 1`, the given inequality together
> with the **Northcott** property imply that the set `𝒮` must be finite.  Hence
> we may assume `k ≥ 2`.  We now apply the Subspace Theorem as in the proof of
> [CZ04, Lemma 1]."

**So for `k = 1` the "CZ Lemma 1" step is elementary.**  The *only*
ineffective step in the whole `k = 1` argument is inside **Proposition 3.4**.

There, the Subspace Theorem is applied to `𝐱 = (x_0, 𝑥_𝒃) = (p, 𝑦_𝒃)` and
produces, verbatim, "a non-trivial linear relation … satisfied by infinitely
many elements of `ℳ₀′`", i.e.

```
a_0·p  +  Σ_𝒃 a_𝒃·𝑦_𝒃  =  0 ,        a_0, a_𝒃 ∈ L   (independent of n).
```

From this the proof derives (15): `p = Σ_{i,j} η_{i,j}α_{i,j}^n` with
`η_{i,j} ∈ L` and `h(η_{i,j}) = o(n)`, then the Galois compatibility (ii), and
finally the properties (i)–(iv).

> **The entire loss of effectivity is the phrase "for all but finitely many
> elements of `ℳ`".**  Everything downstream is elementary.  KMN therefore
> proves `ℳ_θ` is finite but gives **no control on `max ℳ_θ`** — which is
> exactly what `Docs/28` Thm 28.3 and `Docs/30` M7.11b′ require.

**Can the quantitative Subspace Theorem repair it?**  *(Corrected by `Docs/33`.)*
The original claim here was "No — it bounds the **number** `t`, not the
heights/identities."  **That is right for `k*` but wrong as stated for the
cardinality.**  Because `Docs/32` Thm 32.1 gives `r ≤ 2`, the ST is applied in
`L³` and Schlickewei's bound `t ≤ c(3)δ^{-3}(\log(1/\delta)-1)^2` has an
**absolute** `c(3)`.  Combined with `Docs/33` Thm 33.2 (each exceptional
hyperplane carries at most one exception, computable) this yields the **explicit
cardinality bound** `|ℳ_θ| ≤ n_0 + t` (`Docs/33` Thm 33.3).  What counting still
cannot do is *locate* the hyperplanes, hence bound `k*`.

**Original text follows (the height half of the claim is unaffected).**  No.
Schmidt (1989),
Evertse (1995), Evertse–Schlickewei (2002) and Evertse–Ferretti (2013) bound
the **number** `t` of exceptional subspaces (e.g. Schlickewei's
`t ≤ c(n)δ^{-3}(\log(1/\delta)-1)^2`).  **The heights and identities of those
subspaces remain unbounded.**  Counting is not locating.

---

## 31.2 The structural rigidity observation (the promising part)

Here the Pell/biquadratic structure *does* bite — on precisely the step KMN
loses.

> **CORRECTION (`Docs/32` §32.1).**  The sentence below originally concluded
> that the `ζ_j` are "drawn from a finite, explicitly enumerable set".  **That
> is too strong and is withdrawn.**  In fact
> `ζ_j = Σ_{𝒃 : j_2(𝒃)=j_2} a_𝒃·b_{j_1(𝒃),ℓ(𝒃)}`, and the `a_𝒃 ∈ L` are the
> **Subspace Theorem's unknown output**.  So `ζ_j` is fixed only *after* the
> unknown relation is chosen, and Lemma 31.2 is valid **per relation** but
> yields no effectivity on its own.  What is genuinely fixed — and what turns
> out to be decisive — is only the `b_{1,j_1,ℓ}`, i.e. the `K`-coordinates of
> the fixed `q_1 = β_j`.

**Observation 31.1 (corrected).**  In our case `k = 1`, `α_1 = ε`, and **`q_1 = β_j` is
FIXED.**  In KMN's construction the coordinates are

```
𝑦_{(j_1,j_2,ℓ)}  =  b_{1,j_1,ℓ} · ε_{j_2}^k ,
```

where `b_{1,j_1,ℓ}` are the **Q-coordinates of the fixed number `q_1 = β_j`** in
a basis of `K` (Lemma 3.1), and `ε_{j_2}` are the four conjugates of `ε`.
**Hence every `b_{1,j_1,ℓ}` is a fixed rational number, independent of `k`.**

Therefore the relation `a_0p + Σ_𝒃 a_𝒃𝑦_𝒃 = 0` collapses to

```
a_0·p_k  +  Σ_{j=1}^{4}  ζ_j · ε_j^k  =  0 ,        ζ_j ∈ K  FIXED,          (31.1)
```

with coefficients `ζ_j` that are **fixed linear combinations of the fixed
rationals `b_ℓ`** — i.e. drawn from a **finite, explicitly enumerable set**
once the relation `(a_0,a_𝒃)` is fixed.

Now write `p_k = β_jε^k + O(θ^k)` (the rounding error).  Substituting into
(31.1):

```
( a_0β_j + ζ_1 )·ε^k  +  ζ_2 ε_2^k + ζ_3 ε_3^k + ζ_4 ε_4^k  +  a_0·O(θ^k)  =  0 .
```

Because `|ε^k| → ∞`, `|ε_j^k| → 0` (`j ≥ 2`) and `θ^k → 0`:

> **Lemma 31.2 (effective failure of each fixed relation).**  Fix `(a_0,ζ_1,ζ_2,
> ζ_3,ζ_4) ≠ 0` arising as above.  Then
> ```
> a_0·p_k + Σ_j ζ_j ε_j^k = 0
> ```
> holds for at most finitely many `k`, and the largest such `k` is bounded by an
> **explicit function of `(a_0,ζ_1,…,ζ_4)` and `θ`** — obtained by comparing
> successive terms in the display above.

*Proof.*  If `a_0β_j + ζ_1 ≠ 0` the LHS grows like `ε^k` and cannot be
cancelled; if `a_0β_j+ζ_1 = 0` then `ζ_2 ε_2^k` dominates the `O(θ^k)` error
unless `ζ_2 = 0`, and so on.  Each step is an explicit size comparison. ∎

**Why this matters.**  KMN's loss is that it cannot say *which* `k` are
covered.  Lemma 31.2 says that **for each candidate relation, the excluded `k`
are explicitly bounded.**  So the effectivity gap is reduced to a finite,
bounded enumeration problem:

> **M7.12 (restated precisely).**  Obtain an effective bound on, or an
> effective identification of, the exceptional relations produced by
> Proposition 3.4 for `k = 1`, `α_1 = ε = 1+√a_1`, `q_1 = β_j` fixed.  Each
> individual relation is harmless (Lemma 31.2); the difficulty is that the
> quantitative Subspace Theorem bounds only how many there are, not which.

This is a genuine narrowing: the problem is no longer about rates, nor about
the Subspace Theorem's height bound, but about **identifying finitely many
fixed algebraic relations.**

---

## 31.3 And the honest reason the obvious fix fails

One might hope that the rigidity in Lemma 31.2, iterated, closes everything.
It does not, and the reason is worth recording.

From (15), KMN has `p = Σ_j η_jε_j^k` with `η_j ∈ K` and `h(η_j) = o(k)`.
Since `p_k ≈ β_jε^k` and `Σ_jε_j^k = Tr(ε^k) ≈ ε^k`, this forces

```
η_j/β_j  =  1 + O(e^{−ck})        for some c > 0 ,
```

i.e. **`η_j` is an element of the fixed field `K` agreeing with the fixed
quadratic irrational `β_j` to exponential precision.**  Liouville then says
only

```
|η_j − β_j|  ≥  c / den(η_j)²  =  e^{−h(η_j)}  =  e^{−o(k)} .
```

which is **fully compatible** with `e^{−ck}`.

> **The `o(k)` in the Subspace-Theorem output is exactly calibrated to the
> Liouville denominator scale.**  So the biquadratic/Pell structure supplies no
> purchase on that particular gap: any improvement of `h(η) = o(k)` to
> `h(η) = O(log k)` would be needed, and nothing in the fixed-data structure
> delivers it.

Note where the `o(k)` comes from: it is **inherited from the hypothesis**
`h(q_i) < f(n)` with `f` sublinear.  With `q_1 = β_j` fixed the input is far
better than `o(k)`, so there is a plausible route to an improved output — but
establishing it requires a strengthened Subspace-Theorem statement, not a
calculation.

---

## 31.4 Status

| Question | Answer |
|---|---|
| Does KMN's proof contain the needed rate? | **No** |
| How many Subspace-Theorem applications? | **Two** (Prop 2.2, Prop 3.4) |
| Is Prop 2.3's `k=1` case effective? | **Yes — Northcott, no ST** (verbatim, §31.1) |
| The single ineffective step | Prop 3.4: "for all but finitely many elements of `ℳ`" |
| Can the quantitative ST repair it? | **No** — it bounds the **number** `t` of subspaces, not their heights/identities |
| Does the fixed-`q_1` structure help? | **Yes, on that step** — Lemma 31.2: each fixed relation fails for an **explicitly bounded** set of `k` |
| Does it close the gap? | **No** — §31.3: `h(η)=o(k)` is calibrated to Liouville; `O(log k)` is needed and is not obtained |
| **M7.12, precisely** | identify the finitely many relations of Prop 3.4 effectively (count is bounded; identity is not) |

**Nothing here closes M7.12.**  What it does is (a) locate the entire
ineffectivity of the `k = 1` argument in one sentence, (b) show that the
quantitative Subspace Theorem is the wrong repair, and (c) reduce the problem
to identifying finitely many *fixed algebraic relations*, each of which is
already provably harmless with an explicit threshold.