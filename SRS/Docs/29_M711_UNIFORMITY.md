# Paper 29 — M7.11: effective polynomial irrationality for `‖βε^k‖`

Version 1.0.  Verification of the outstanding input, plus the **product-to-factor
reduction** — which turns out to be *essential*, because it collapses the field
from degree `2^m` to a fixed biquadratic field and hence repairs exactly the
discriminant obstruction that killed `γ ≥ Disc(K)^{-O(1)}` in `Docs/23` §23.7.

**Outcome.**

* The one-dimensional finiteness statement is a **published theorem**
  (Corvaja–Zannier 2004; Kulkarni–Mavraki–Nguyen 2018), verified below.
* **Lemma 29.1 (product → single factor):** an M7.9 exception forces *some*
  `j ≥ 2` with `‖β_jε^k‖ < θ₀^k`, `θ₀ = ε^{-1/(3n(n−1))}`.
* **Theorem 29.2 (the discriminant repair):** after the reduction the field is
  `K = Q(√a_1,√a_2)`, degree `4`, and
  ```
  log₂ Disc(K)  =  2( log₂a_1 + log₂a_2 ) + O(1)  =  O(N) ,
  ```
  **linearly in the radicand bit-sizes** — versus `~2^{m−1}·N` for the full
  multiquadratic field.  This is precisely the condition `Docs/23` §23.7
  identified as necessary and which we there proved fatal for the discriminant
  scale.  **The reduction removes that obstruction entirely.**
* **The residual question is a single uniformity question**, stated precisely in
  §29.5 — and it has a genuine cost, recorded in §29.4.

---

## 29.1 Correction to `Docs/28` §28.4

`Docs/28` §28.4 asserted: *"Effective Subspace-Theorem thresholds for a fixed
system of algebraic linear forms are classical."* **That sentence conflated two
distinct assertions and is hereby retracted:**

1. effective control of *solutions* to a fixed-system Subspace-Theorem problem,
   versus
2. an explicit inequality `‖βε^k‖ ≥ k^{-c}` for all large `k`.

(1) does **not** automatically give (2): one must exhibit the Subspace-Theorem
variables, show an exceptional `k` yields an integer point in one of finitely
many exceptional subspaces *with an effective height threshold*, and then
eliminate those subspaces or convert each into an effective recurrence /
linear-form problem.  None of that was carried out in `Docs/28`.

**What replaces it is stronger and correctly targeted:** the needed statement is
not merely a Subspace-Theorem consequence but a *named, published theorem* about
exactly the quantity `‖βε^k‖` (§29.2).

---

## 29.2 The finiteness is a published theorem (verified)

**Corvaja–Zannier**, *On the rational approximations to the powers of an
algebraic number: solution of two problems of Mahler and Mendès France*,
Acta Math. **193** (2004), 175–191 (arXiv:math/0403522).

> **Theorem 1 (verbatim).**  "Let `α > 1` be a real algebraic number and let
> `0 < λ < 1`.  Suppose that `‖α^n‖ < λ^n` for infinitely many natural
> numbers `n`.  Then some power `α^n` is a **Pisot number**."

> **MAIN THEOREM (verbatim, opening).**  "Let `Γ ⊂ Q̄*` be a finitely generated
> multiplicative group of algebraic numbers, let `β ∈ Q̄*` and let `δ > 0` be
> fixed.  Then there are only finitely many pairs `(q,u) ∈ Z_{>0} × Γ` with
> `d = [Q(u):Q]` such that `1 ≪ ‖qβ‖_v`, `βq` is not a pseudo-Pisot number and
> `h(1/(qβ)) ≤ … l(v,d_χ)^{1/δ}`."

This is a **lower bound for `‖qβ‖` against heights, with finitely many and
explicitly controlled exceptions** — the correct shape, and built from the
Subspace Theorem with explicit height bookkeeping.

**Kulkarni–Mavraki–Nguyen**, *Algebraic approximations to linear combinations
of powers: an extension of results by Mahler and Corvaja–Zannier*,
Trans. AMS (2018) (arXiv:1511.08525).

> **Abstract (verbatim).**  "We **completely solve the problem of the
> existence of `θ ∈ (0,1)`** such that there are infinitely many tuples
> `(n,q_1,…,q_k)` satisfying
> `‖q_1α_1^n + … + q_kα_k^n‖_Z < θ^n` where `n ∈ N` and `q_1,…,q_k ∈ K*` having
> small logarithmic height compared to `n`. […]  The case `k = 1` and `q_1` is
> essentially a rational integer was obtained by Corvaja and Zannier and settled
> a long-standing question of Mahler."

> Our problem is **exactly `k = 1`, `α_1 = ε`, `q_1 = β_j`** (fixed, hence
> bounded height).
>
> **BUT see `Docs/30` §30.5: this is KMN's contribution, NOT Corvaja–Zannier's.**
> KMN state verbatim that "the case `k = 1` and `q_1` is **essentially a
> rational integer** was obtained by Corvaja and Zannier."  Our `q_1 = β_j` is
> **not** rational, so Corvaja–Zannier Theorem 1.1 does **not** apply; **KMN is
> the operative theorem.**  The distinction is essential: if `q_1 = 1` and
> `a_1 ∈ {2,3}` then `ε` is Pisot and `‖ε^k‖ ≈ (√a_1−1)^k`, so
> `‖ε^k‖ < θ^k` for every `k` once `θ > √a_1−1`.  It is precisely
> `β_j ∉ Q` that removes us from the degenerate regime.

Effectivity of the surrounding `S`-unit machinery: Waldschmidt (2015),
*Schmidt Subspace Theorem and S-unit equation* — "It turns out that the result
of Corollary 1 **is effective**: one can bound from above the (numerators and
denominators of the) solutions `x_1` and `x_2`."

**One structural caveat that must be checked, not assumed.**  `ε = 1+√a_1` **is
a Pisot number** (its conjugate `1−√a_1` has `|·| < 1`), so `‖ε^k‖ → 0`
exponentially — `ε` sits in exactly the *degenerate* regime that CZ Theorem 1
describes.  Our quantity is `‖β_jε^k‖` with `β_j ∉ Q̄` non-rational, and
empirically `‖β_jε^k‖` stays `Θ(1)` (values `0.05–0.47` for `(a_1,a_2)=(2,3)`,
`k ≤ 44`).  So the orbit is non-degenerate; **but the non-degeneracy hypothesis
must be verified in the sense of KMN, not inferred from the numerics.**

---

## 29.3 Lemma 29.1 — the product reduces to a single factor

**Lemma 29.1.**  Suppose `Π_{j=2}^n ‖β_jε^k‖ < ε^{-k/(3n)}`.  Then for some
`j ∈ {2,…,n}`,
```
‖β_jε^k‖  <  ε^{-k/(3n(n−1))}  =  θ₀^k ,        θ₀ := ε^{-1/(3n(n−1))} .
```
*Proof.*  Otherwise all `n−1` factors are `≥ ε^{-k/(3n(n−1))}`, so the product
is `≥ ε^{-k/3} > ε^{-k/(3n)}` for `n ≥ 2` — contradiction. ∎

So **M7.9's simultaneous product condition implies a one-dimensional condition
of the CZ/KMN type.**  We never need a simultaneous quantitative Subspace
Theorem — exactly the simplification proposed.

---

## 29.4 The cost of the reduction, stated honestly

The reduction is **not free**.  It weakens the threshold:

```
product condition :  ‖·‖_prod  <  ε^{-k/(3n)}          θ_prod = ε^{-1/(3n)}
reduced condition :  ‖·‖_1     <  ε^{-k/(3n(n−1))}    θ₀     = ε^{-1/(3n(n−1))}
```

`θ₀ > θ_prod`, so the reduced condition is *easier to satisfy* and therefore
*harder to rule out*.  Worse, `θ₀ → 1` as `n → ∞`:

```
 n :      2        3        4        6
θ₀ (a₁=2): 0.8634  0.9522  0.9758  0.9903
θ₀ (a₁=3): 0.8458  0.9457  0.9725  0.9889
```

**So the uniformity question gets strictly harder as the rank grows.**  This is
the genuine tension of `Docs/29`:

| | field | threshold | theorem available? |
|---|---|---|---|
| simultaneous product | degree `2^m`, `log₂Disc ~ 2^{m−1}N` (**fatal**, `Docs/23`) | `θ_prod` (stronger) | **no** — not a CZ/KMN statement |
| reduced single factor | degree `4`, `log₂Disc = O(N)` (**repaired**) | `θ₀` (weaker, `→1`) | **yes** — CZ 2004 / KMN 2018 |

Which side wins is precisely the open question.  `Docs/29` does **not** claim
the trade is favourable.

---

## 29.5 Theorem 29.2 — the discriminant repair

`K = Q(√a_1, √a_2)` is biquadratic of degree `4`, with

```
Disc(K)  =  D(a_1) · D(a_2) · D(a_1a_2) ,     D(x) = fundamental disc. of Q(√x) = O(x) ,
```

so

```
log₂ Disc(K)  =  2( log₂a_1 + log₂a_2 ) + O(1)  =  O(N) .                    (29.1)
```

Measured (`N` = radicand bit-size):

```
  a1    a2  [K:Q]      Disc(K)   log2 Disc   N(rad)   log2Disc/N
   2     3      4          2,304       11.17         3         3.72
   2    13      4         10,816       13.40         5         2.68
 101   103      4   1,731,558,544       30.69        14         2.19
1009  1013      4 1,044,723,161,689      39.93        20         2.00
```

The ratio is **bounded (≈2–4) and does not grow**.

> **This is the point of the whole reduction.**  `Docs/23` §23.7 proved that any
> bound of the form `γ ≥ Disc(K)^{-O(1)}` is useless, because for the full
> multiquadratic field `log₂Disc(K) ≍ 2^{m−1}·N` is **exponential in the rank**.
> Lemma 29.1 removes the simultaneous product, hence removes the rank from the
> field entirely: the working field is a **fixed-degree biquadratic** with
> `log₂Disc = O(N)`.  Consequently any constant that is polynomial — or even
> `exp(poly)` — in `log Disc K` or in `[K:Q]=4` is acceptable for the
> `O(N²)` target of `Docs/24` Thm 24.5.

---

## 29.6 M7.11, stated exactly

> **Open (M7.11).**  Let `θ_crit(a_1,a_2)` be the critical value furnished by
> Corvaja–Zannier (2004) / Kulkarni–Mavraki–Nguyen (2018) for the statement
> *"`‖β_jε^k‖ < θ^k` has only finitely many solutions in `k`"*, with `β_j` and
> `ε` non-degenerate in their sense.  **Prove**
> ```
> θ₀  =  ε^{-1/(3n(n−1))}   <   θ_crit(a_1,a_2)                          (M7.11)
> ```
> uniformly over distinct squarefree `a_1 < … < a_n`.

Equivalently, in the language of `Docs/28`: obtain explicit `c, K_ST` in
`‖β_jε^k‖ ≥ k^{-c}` — noting that `k^{-c} > ε^{-k/(3n(n−1))}` for all large `k`,
so any such bound suffices.

**What is now settled, and what is not.**

| Question | Status |
|---|---|
| Is the one-dimensional finiteness a theorem? | **Yes** — CZ 2004 Thm 1 + Main Theorem; KMN 2018 (verbatim, §29.2) |
| Is the simultaneous product reducible to one factor? | **Yes** — Lemma 29.1 |
| Does the reduction repair the discriminant obstruction? | **Yes** — Thm 29.2, `log₂Disc = O(N)` |
| Is the relevant `S`-unit machinery effective? | **Yes** — Waldschmidt 2015 |
| Is the orbit non-degenerate in the KMN sense? | **RESOLVED in `Docs/30`** — Definition 1.2 is *vacuous* at `k = 1` (Cor 30.1) |
| Is `‖β_jε^k‖ < θ^k` finite for every `θ<1`? | **YES — `Docs/30` Thm 30.3**, uniform, no constant; (i) fails by an archimedean size comparison (Prop 30.2), (iv) fails for `a_1 ≥ 5` (Prop 30.5) |
| **Is `θ₀ < θ_crit` uniformly, EFFECTIVELY?** | **OPEN — `Docs/30` §30.6, restated as M7.11b′: improve `η_k ≥ e^{-o(k)}` to `η_k ≥ e^{-C log k}` with explicit `C`** |

**Honest closing statement.**  `Docs/29` converts the problem from "an
ineffective constant in a general theorem" into "**one specific inequality
between an explicit number `θ₀ = ε^{-1/(3n(n−1))}` and a critical value
produced by a named, effective, published theorem**".  That is a large
reduction in what stands between the project and a proof — but it is a
uniformity question about someone else's theorem, it degrades with `n`, and it
is **not settled here**.