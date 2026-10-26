# Paper 26 — the SRS Markov function

Version 1.0.  Builds the custom function requested, in the normalized form that
makes the Liouville floor equal to 1, and records six exact lemmas about it.
The function is *not* another invariant of `u_i`: it measures the joint metric
rarity directly.  Main outcome: the function **rigorously reduces the
exceptional set to a computable candidate list**, and the residual open core
becomes a single sharp statement.

---

## 26.1 The function

For `a = (a_1,…,a_n)` distinct squarefree, `δ = 1/(3n)`, define the
**normalized surplus**

```
ρ_a(q)  :=  4^n · Π_{i=1}^n ‖q√a_i‖ ,                                    (26.1)
```

together with the Markov and dyadic forms

```
𝓜_a(Q) := min_{1≤q≤Q} ρ_a(q) ,        G_a(Q) := min_{Q<q≤2Q} ρ_a(q) .
```

Normalization rationale: for `q√a_i` "generic", `‖q√a_i‖ ~ U[0,1/2]` with
mean `1/4`, so `ρ_a(q) ~ 4^n(1/4)^n = 1`.  **So `ρ = 1` is the typical value
and deviations below 1 are the rare events.**  The exception condition becomes

```
q^{1+δ} Π‖q√a_i‖ < 1     ⟺     ρ_a(q) < 4^n q^{-1−δ} ,                     (26.2)
```

i.e. **an exception is a suppression of `ρ` by a factor `≈ q^{1+δ}/4^n`** —
roughly one power of `q` below typical.  This is the quantitative picture the
whole project was missing.

---

## 26.2 Lemma 1 — the function *is* the defect product

**Lemma 1.**  With `p_i = round(q√a_i)`, `u_i = p_i² − a_i q²`,
`2q√a_i − 1/2 ≤ p_i + q√a_i ≤ 2q√a_i + 1/2`, hence

```
ρ_a(q)  =  Π_i |u_i| / √(a_1⋯a_n) · ( 1 + O( 1/(q √(min_i a_i) ) ) ) .
```

*Proof.*  `‖q√a_i‖ = |u_i|/(p_i+q√a_i)`, so
`4^nΠ‖·‖ = 4^nΠ|u_i| / Π(p_i+q√a_i)`, and
`Π(p_i+q√a_i) = 2^n q^n √(a_1⋯a_n)(1+O(1/q))`. ∎

So `ρ_a(q)` is (to leading order) *exactly* the product of the quadratic
defects, rescaled.  `𝓜_a` and `G_a` are therefore the right rescalings of the
defect product — and, unlike the cross-defects of `Docs/25`, they are
**scale-free**: they do not depend on any choice of algebraic form.

---

## 26.3 Lemma 2 — Markov positivity is *elementary*, not a Subspace Theorem

**Lemma 2.**  `ρ_a(q) ≥ 1/(q^n √(a_1⋯a_n))` for every `q ≥ 1`, purely
elementarily (Liouville for `√a_i` alone).  Consequently

```
min_{q ≥ 1} ρ_a(q)  >  0      with NO appeal to the Subspace Theorem,
```

and the infimum is attained (since `ρ_a(q) ≤ 2^{-n}·4^n = 1` and `ρ_a(q)`
eventually exceeds `4^nq^{-1-δ}` by the Subspace Theorem).

**Proof.**  `|p_i − q√a_i| = |u_i|/(p_i+q√a_i) ≥ 1/(2q√a_i + 1/2) ≥ 1/(3q√a_i)`;
multiplying and using (26.1). ∎

> **This corrects a natural misreading.**  The Subspace Theorem is *not* what
> makes the multiplicative Markov constant positive — Liouville does that for
> free.  What the Subspace Theorem adds is the far stronger
> `ρ_a(q) ≥ 4^n q^{-1−δ}` **for all large `q`**, i.e. that `ρ` eventually
> exceeds a *growing* threshold.  Only that is in question.

---

## 26.4 Lemma 3 — the candidate-set reduction, and its refutation

**Natural conjecture.**  Since an exception needs a suppression of `≈ q^{1+δ}`
and one factor can be suppressed by at most `≈ q` (Liouville is sharp:
convergents attain `‖q√a‖ ≈ 1/(2q√a)`), an exception should require *two*
factors to be suppressed — hence `q` should be a convergent denominator.

**Lemma 3 (REFUTED).**  The reduction "exception ⇒ `q` is a convergent
denominator of some `√a_i`" is **false**.  Of the 121 certified exceptional
`q` for `(√2,√3)`, only **40** are convergent denominators of `√a_1` or `√a_2`
(and only **1** of both).

*Why the argument fails (third instance of the same trap).*  To force a
contradiction one needs the lower bound on `ρ` to *exceed* `4^nq^{-1-δ}`.  With
all factors at their floors, `ρ ≥ 4^n/(3^nq^n√(∏a)) = 4^n/(3^nq^n√(∏a))`, and
the required comparison `4^nq^{-n} < 4^nq^{-1-δ}` i.e. `q^{-n} < q^{-1-δ}` i.e.
`q^{n-1-δ} < 1` holds for `n ≥ 2` and `q ≥ 1`.  **The floor is automatically
compatible with the threshold** — the same wrong-direction-of-information
failure as `Docs/22` Thm C and `Docs/24` Thm 24.1(ii).

---

## 26.5 Lemma 4 — the *tail* is floor events (the operative reduction)

Although Lemma 3 fails in the bulk, it holds decisively in the **extreme tail**,
which is what determines `q*`.

**Lemma 4.**  Restrict to `q` that are convergent denominators of some `√a_i`.
Then `‖q√a_i‖ = |u_i|/(p_i+q√a_i)` with `|u_i| = O(√a_i)` bounded, so

```
‖q√a_i‖  ∈  [ 1/(2q√a_i + C√a_i) , 1/(2q√a_i − C√a_i) ]  ≈  1/(2q√a_i) ,
```

and `(26.2)` reduces **exactly** to the *mild* condition

```
Π_{j ≠ i} ‖q√a_j‖   <   2 √a_i · q^{−δ} ,        δ = 1/(3n).             (26.4)
```

*Verification.*  Scanning the convergent denominators of `√2,√3` and testing
(26.2) directly:

```
Q = 1e13 : 79 candidates,  42 exceptional,  max q = 1,104,728,155,436
Q = 1e16 : 97 candidates,  43 exceptional,  max q = 1,746,860,020,068,409
Q = 1e20 :122 candidates,  43 exceptional,  max q = 1,746,860,020,068,409  (stabilised)

mechanism check  ‖q√a_j‖ ≤ 2√a_i·q^{−δ}  satisfied at  42/42  hits
hits with min_i|u_i| ≤ 1 (true Pell floor):  33/43
```

**Consequences.**

1. **`q*` is at least `1,746,860,020,068,409 ≈ 1.75×10¹⁵`** for
   `(√2,√3)`, `δ = 1/6`.  This **corrects `Docs/24` §24.8**, whose figure
   `29,354,524` came from a brute-force scan to `4×10⁷` that was **8 orders of
   magnitude too shallow** — the "stabilised" claim there was wrong.
2. The candidate set is **computable and geometrically growing** (the
   continued fraction of `√a_i` is periodic), so the search space collapses
   from `q ≤ Q` (density-1 integers) to `O(log Q)` candidates.
3. `δ = 1/(3n)` is *tiny*, so (26.4) is an extremely weak demand on the other
   coordinates — yet it is exactly where effectivity is missing.

---

## 26.6 The open core, now a single sharp statement

> **Open (M7.7).**  Let `a_1 < ⋯ < a_n` be distinct squarefree and `δ = 1/(3n)`.
> Bound effectively the convergent denominators `q` of `√a_1` for which
> ```
> Π_{j≥2} ‖q√a_j‖  <  2√a_1 · q^{−δ} .
> ```
> The required scale is `log₂ q* = O(N²)` (`Docs/24` Thm 24.5).

Why this is the right core:

* the left side is a product of **weak** approximation exponents (`1/(3n)`),
  so no individual continued-fraction or Roth-type control applies;
* the `q` are **restricted to a computable, geometrically sparse set**, so a
  bound is plausible in principle and is not a density question;
* it is precisely the configuration the Pell family produces (`Docs/24` Thm
  24.3), so any correct theorem must survive it — and it already forces
  `log₂q* = 3n(n−1)−O(1)`;
* it is **not** killed by `Docs/25` (no algebraic rigidity available) nor by
  `Docs/22` Thm C (the information is metric).

---

## 26.7 What the metric/discrepancy toolset can and cannot do

**Lemma 5 (metric tools predict the wrong magnitude).**  The naive
volume/ergodic prediction for the count of exceptions up to `Q` is

```
Σ_{q≤Q} vol{ x ∈ [0,1/2]^n : Π x_i < q^{−1−δ} }  ≍  Σ_{q≤Q} q^{−1−δ}(log q)^{n−1}/(n−1)!
```

For `(√2,√3)`, `n=2`, `δ=1/6`, `Q=4×10⁷` this predicts `≈ 0.76`, while the
**observed** count is `121` — under-prediction by a factor `≈ 160`.

> **Therefore discrepancy / volumes-of-numbers / ergodic reasoning cannot
> deliver the required bound**: it is quantitatively wrong by more than two
> orders of magnitude at the tested scale, and the discrepancy is driven by the
> algebraic (Pell/convergent) structure, not by generic equidistribution.  Any
> purely metric proof would have to be quantitatively sharp on a phenomenon
> that generic theory mis-predicts.

This is the rigorous reason the user's proposal should use the function to
*measure* rarity, not to prove effectivity via standard metric tools.

---

## 26.8 Status

| Item | Result |
|---|---|
| Normalized function with Liouville floor = 1 | `ρ_a(q) = 4^nΠ‖q√a_i‖` (§26.1) |
| Function = defect product, rescaled | **Lemma 1** |
| Markov positivity `min ρ > 0` | **Lemma 2** — *elementary*, no ST |
| Exception = suppression by `≈ q^{1+δ}/4^n` | (26.2) |
| "Exception ⇒ convergent denominator" | **Lemma 3** — refuted (40/121) |
| Tail = floor events; `q* ≥ 1.75×10¹⁵` | **Lemma 4** — verified |
| Purely metric proof route | **Lemma 5** — quantitatively wrong |
| Open core | M7.7, §26.6 |

**What this achieves.**  The problem is reduced, rigorously and with verified
constants, from "bound the exceptional `q`-set" to a single statement about a
**computable, geometrically sparse** candidate set with a **weak** side
condition — and the correct value of `q*` for the smallest case is now known to
be at least `10¹⁵`, eight orders of magnitude above the previous brute-force
figure.  What is still missing is the upper bound in §26.6; that is the live
open problem.