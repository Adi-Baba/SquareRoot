# Paper 24 — M7.5: attacking the EHS exceptional `q`-set

Version 1.0.  A dedicated attack on the single remaining obstruction identified
in `Docs/23`: the ineffective constant `γ` in Eisenbrand–Haeberle–Singer,
equivalently the largest exceptional `q` in

```
q^{1+δ} · Π_{i=1}^n dist_ℤ(q √a_i)  <  1 ,        δ = 1/(3n).              (EHS)
```

**Main results of this paper:**

1. An exact integer reformulation (E) of (EHS).
2. **Theorem 24.1** — a sharp boundary: for `n = 1` the exceptional set is
   elementary and effective; for `n ≥ 2` every per-coordinate Liouville bound
   is *provably useless*.  (This is where I first went wrong by a sign error,
   caught numerically; §24.4 records that.)
3. **Theorem 24.3 (new, verified)** — a **Pell lower bound**: the exceptional
   set provably contains convergent denominators up to `2^{Ω(n²)}`.
4. **Theorem 24.4** — **the target of `Docs/23` §23.7 is REFUTED.**  No
   absolute `C` exists with `γ ≥ (d_1⋯d_m)^{-C}`.
5. **Theorem 24.5** — the corrected target, which survives: `log₂ q* = O(N²)`.

---

## 24.1 Setup

`a_1 < ⋯ < a_n` distinct squarefree `>1`, `δ = 1/(3n)` (EHS's choice), and

```
E(δ) := { q ∈ ℕ⁺ : q^{1+δ} · Π_i dist_ℤ(q√a_i) < 1 } .
```

By the Subspace Theorem `E(δ)` is finite.  Write `q* := max E(δ)`.

The chain from `q*` to `γ`: EHS Thm 4 + Thm 8 give
`|E_sum| ≥ Q^{-(n+1)}` for `Q = Q₀·(2n‖x‖∞)^{3/2}` provided `Q₀ > q*`, hence
`γ = Q₀^{-(n+1)-O(n)}` and therefore

```
log₂(1/γ)  =  (n+1)·log₂Q₀  +  O(n)  ≥  (n+1)·log₂q*  +  O(n).           (24.1)
```

So **`log₂ q*` is exactly the quantity to bound**, and any *lower* bound on
`q*` is a *lower* bound on how strong the sought theorem must be.

---

## 24.2 Exact integer reformulation

Let `p_i := round(q√a_i)` and `n_i := p_i² - a_i q² ∈ ℤ`.  Since `a_i` is
nonsquare, `n_i ≠ 0`, so `|n_i| ≥ 1`.  And

```
dist_ℤ(q√a_i)  =  |n_i| / (p_i + q√a_i),        p_i + q√a_i ≥ 2q√a_i .
```

Hence (EHS) is *equivalent* to

```
(E)   Π_{i=1}^n |p_i² − a_i q²|   ≲   2^n √(a_1⋯a_n) · q^{n−1−δ} .
```

> **Interpretation.**  `E(δ)` is the set of common denominators `q` at which
> **all `n` quadratic forms `p² − a_i q²` are simultaneously small**, in the
> product sense.  This is the honest arithmetic content of the EHS constant.

---

## 24.3 The sharp `n = 1` boundary (Theorem 24.1)

**Theorem 24.1.**  Let `δ > 0` and `a > 1` squarefree.

* **(i)** `E₁(δ) := { q : q^{1+δ} dist_ℤ(q√a) < 1 }` is finite, with the
  **explicit** bound `q < (2√a + 1)^{1/δ}`.
* **(ii)** For `n ≥ 2`, if `dist_ℤ(q√a_i) ≥ c/(q√a_i)` for all `i`, then
  (EHS) implies only `q^{n−1−δ} > (1/c)^n (a_1⋯a_n)^{−1/2}`, a **lower**
  bound on `q`, hence no constraint whatsoever.

*Proof of (i).*  `dist_ℤ(q√a) = |p−q√a| = |n|/(p+q√a) ≥ 1/(p+q√a)`, and
`p ≤ q√a + 1`, so `p+q√a ≤ 2q√a+1`.  Hence `dist ≥ 1/(2q√a+1)`, so
`q^{1+δ}/(2q√a+1) < 1`, i.e. `q^δ < 2q√a+1 ≤ 3q√a` for `q ≥ 1`, giving
`q^{1−δ} < 3√a`; using `q^δ ≤ q` for `δ ≤ 1` gives `q < (3√a)^{1/(1−δ)}`.
To get the clean form `q < (2√a+1)^{1/δ}`: `dist ≥ 1/(2q√a+1)` and
`q^{1+δ}/(2q√a+1) < 1` ⟹ `q^{1+δ} < 2q√a+1 ≤ (2√a+1)q` ⟹ `q^δ < 2√a+1`. ∎

*Proof of (ii).*  `Π_i dist ≥ c^n q^{-n}(a_1⋯a_n)^{-1/2}`.  (EHS) then forces
`q^{1+δ−n} < c^{-n}(a_1⋯a_n)^{1/2}`.  For `n ≥ 2` the exponent `1+δ−n` is
**negative**, so this is automatically satisfied for large `q`: it is a lower
bound, not an upper bound. ∎

**This is the wall.**  A single quadratic irrational is separated from the
integers at rate `~1/q`, which already beats the `q^{-1−δ}` threshold of
(EHS).  But `n` of them only give `~q^{-n}`, which for `n ≥ 2` is *far worse*
than `q^{-1−δ}`.  The threshold is therefore only ever attained by **rare
simultaneous** near-convergents — and no per-coordinate estimate can detect
rarity.

---

## 24.4 A sign error I made, and why the numerics caught it

I first "proved" `q* < (4^n√(a_1⋯a_n))^{1/(n−1−δ)}`, i.e. an explicit bound
finishing the whole problem.  The step was

```
Π_i dist ≥ 4^{-n} q^{-n} (∏a_i)^{-1/2}   (a LOWER bound)
(EHS)  q^{1+δ} Π_i dist < 1
  ⟹  q^{1+δ} · 4^{-n} q^{-n} (∏a_i)^{-1/2} < 1
  ⟹  q^{1+δ−n} < 4^n √(∏a_i)
  ⟹  "q^{n−1−δ} < 4^n √(∏a_i)"          ← WRONG: inverting flips the sign
```

The correct conclusion is `q^{n−1−δ} > 1/(4^n√(∏a_i))`, a **lower** bound.
An exhaustive integer scan immediately refuted the claim: for `(√2,√3)`,
`δ=1/6`, the predicted bound was `q < 81.6` while exceptional `q` reach
`29,354,524`.

**This is the same structural obstruction as `Docs/22` Theorem C, one level up:**
a lower bound on a product, combined with a requirement that the product be
*small*, is arithmetically inert.  Recorded because the error was tempting and
because it is exactly the mistake any future attempt will make.

---

## 24.5 The Pell lower bound (Theorem 24.3)

The exceptional set is *large*, and provably so.

**Theorem 24.3.**  Let `n ≥ 2`, `δ = 1/(3n)`, `a_1<⋯<a_n` distinct squarefree
`>1`.

* **(a)** Every **convergent denominator** `q` of `√a_i` with
  `q < 2^{3n(n−1)}` lies in `E(δ)`.
* **(b)** Every **Pell denominator** `q` of `√a_i` (i.e. `p²−a_iq² = ±1`) with
  `q < (2^n √a_i)^{3n}` lies in `E(δ)`.

*Proof of (a).*  For a convergent, `dist_ℤ(q√a_i) < 1/q`.  For `j ≠ i`,
`dist ≤ 1/2`.  Hence
`q^{1+δ}Π_j dist_j < q^{1+δ}·(1/q)·2^{-(n−1)} = q^{δ}2^{-(n−1)}`,
which is `< 1` iff `q < 2^{(n−1)/δ} = 2^{3n(n−1)}`. ∎

*Proof of (b).*  For a Pell denominator, `dist = 1/(p+q√a_i) ≤ 1/(2q√a_i)`.
Then `q^{1+δ}Π_j dist_j ≤ q^{1+δ}(2q√a_i)^{-1}2^{-(n−1)} = q^{δ}/(2^n√a_i)`,
which is `< 1` iff `q < (2^n √a_i)^{1/δ} = (2^n√a_i)^{3n}`. ∎

**Verification (exhaustive, exact integer arithmetic).**  For
`n = 2,3,4` and radicand sets `(2,3)`, `(2,3,13)`, `(2,3,5,7)`, **both (a) and
(b) hold with zero exceptions**:

```
a=(2,3)      n=2 : (a)  5+8=13 tested, 13/13 OK      (b) 7+5=12 tested, 12/12 OK
a=(2,3,13)   n=3 : (a) 15+20+20=55,  55/55 OK        (b) 16+11+4=31,   31/31  OK
a=(2,3,5,7)  n=4 : (a) 29+39+18+37=123, 123/123 OK   (b) 31+20+19+10=80, 80/80 OK
```

And the predicted scales match the observed behaviour:

```
n=3, (2,3,13), (b) for sqrt2 : threshold 3.04e9,  largest Pell q below = 1,311,738,121
n=4, (2,3,5,7), (b) for sqrt2: threshold 1.80e16, largest Pell q below = 1.02e16
```

which explains why the brute-force scan to `4×10⁷` was *still finding new
exceptional `q`* at `n = 3` — the true `q*` there is far beyond `10⁹`.

---

## 24.6 Consequence: the `Docs/23` target is refuted (Theorem 24.4)

`Docs/23` §23.7 proposed

```
γ ≥ (d_1⋯d_m)^{-C}   for an absolute C,      i.e.   log₂(1/γ) ≤ C·Σ log₂ d_i .
```

**This is FALSE.**

**Theorem 24.4.**  No absolute constant `C` satisfies `γ ≥ (d_1⋯d_m)^{-C}`.

*Proof.*  By Theorem 24.3(a) the exceptional set contains convergent
denominators up to `2^{3n(n−1)}`, so `log₂ q* ≥ 3n(n−1) − O(1)`.  By (24.1),
`log₂(1/γ) ≥ (n+1)(3n(n−1)−O(1)) = Ω(n³)`.  On the other hand, `d_i` distinct
squarefree means `Σ log₂ d_i ≥ log₂(n!) = Θ(n log n)` (the `d_i` are at least the
first `n` primes).  Hence

```
log₂(1/γ) / Σ log₂ d_i   ≥   Ω( n³ / (n log n) )  =  Ω( n² / log n )  →  ∞ .
```

So the ratio is unbounded and no absolute `C` exists. ∎

**Numerically** (using the weaker but safer version `(a)`):

```
 n   log2 q* >=      log2(1/gamma) >=      N >= log2(n!)    ratio
 2          6                  18               1.0          18.0
 3         18                  72               2.6          27.9
 4         36                 180               4.6          39.3
 5         60                 360               6.9          52.1
 6         90                 630               9.5          66.4
 8        168                1512              15.3          98.8
10        270                2970              21.8         136.3
```

The ratio grows without bound, as the proof requires.

**So `Docs/23` §23.7 is hereby retracted.**  It was a *necessary*-looking
plausible guess, and it is false — falsified by an elementary Pell argument
that took one page.

---

## 24.7 What survives: the corrected target (Theorem 24.5)

SRS only needs `log₂(1/γ) ≤ poly(N)`, not a linear bound.  With
`N ≥ Σ log₂ d_i ≥ log₂(n!) = Θ(n log n)`:

**Theorem 24.5 (corrected target).**  It suffices to prove

```
log₂ q*  =  O( n² + n·log₂(max_i d_i) + n·log n ) .                        (T24)
```

Indeed (24.1) then gives `log₂(1/γ) = O(n³ + n² log max d + n² log n)`, and
since `n ≤ N/log n` and `log max d ≤ N`,

```
log₂(1/γ)  =  O(N³) ,      and      E  ≤  O(N²) + O(N³)  =  O(N³) .
```

> **The corrected falsifiable target is therefore: an effective bound
> `q* ≤ exp(O(N²))`, i.e. SRS with exponent 3.**
>
> The scale `n²` in (T24) is **forced**: Theorem 24.3 shows `log₂ q* ≥
> 3n(n−1)`.  So the target is tight up to the value of the constant and the
> `n·log max d` term.  Any proposed theorem with `log₂ q* = o(n²)` is
> **already refuted** by Theorem 24.3.

This is the sharpest statement currently available, and it is a genuine
advance over `Docs/23`: it supplies a **lower bound** on the unknown constant,
a **refutation** of the natural guess, and a **tight** reformulated target.

---

## 24.8 Empirical structure of the exceptional set

For `(√2,√3)`, `δ = 1/6`, exhaustive integer scan to `4×10⁷`:

```
|E| = 121 ,   q* = 29,354,524   (interior — the set has stabilised)
defects n_1 = p_1^2−2q^2 ,  n_2 = p_2^2−3q^2   at the top of the list:
      195,025        n_1 =        −1     n_2 =   −141,026
      470,832        n_1 =         1     n_2 =    88,353
      564,719        n_1 =   −429,233     n_2 =         1
    2,107,560        n_1 =   384,400     n_2 =         1
   15,994,428        n_1 =         1     n_2 = 3,656,692
   29,354,524        n_1 = 3,503,204     n_2 =         1
```

Structure (after fixing a bug in my continued-fraction routine — the first
version was wrong and gave a false "not convergents" conclusion):

* convergent denominator of `√2`: 18/121;  of `√3`: 23/121;
* of **at least one**: 40/121;  of **both**: **1/121**.

> **So ~2/3 of the exceptional `q` are genuine simultaneous near-convergents,
> explained by neither individual continued fraction.**  The phenomenon of
> §24.2 is real and is not a continued-fraction artefact of one coordinate.
> The `n_i` are `~q^{0.4}`–`q^{0.9}`: far above the Pell floor `±1`, far below
> `q` — the "soft" regime where no sharp theory exists.

---

## 24.9 Status

> **CORRECTION (`Docs/26` Lemma 4).**  The `q* = 29,354,524` figure in
> §24.8 was **wrong**: the brute-force scan to `4×10⁷` was ~8 orders of
> magnitude too shallow, and the set had **not** stabilised.  Scanning the
> convergent denominators (a computable, geometrically sparse candidate set)
> to `10²⁰` gives **43** hits with
> ```
> q*  ≥  1,746,860,020,068,409  ≈  1.75×10^15   for (sqrt2,sqrt3), delta=1/6.
> ```
> The structure is: `q` a convergent denominator of one `√a_i` (so
> `‖q√a_i‖ = |u_i|/(p_i+q√a_i)` with `|u_i| = O(√a_i)` bounded), and the
> remaining coordinates satisfy only the **mild** condition
> `Π_{j≠i}‖q√a_j‖ < 2√a_i·q^{−δ}` — verified at 42/42 hits.  See `Docs/26`.

| Question | Answer |
|---|---|
| Is `q*` effectively bounded at all? | **Open** (this is the ST effectivity problem) |
| Can per-coordinate Liouville bounds do it? | **No** — Theorem 24.1(ii), proven |
| Is `γ ≥ (∏d)^{-C}` for absolute `C`? | **No** — Theorem 24.4, proven |
| How large must any bound be? | `log₂ q* ≥ 3n(n−1)` — Theorem 24.3, proven |
| What target survives? | `log₂ q* = O(N²)` ⟹ SRS with exponent 3 — Theorem 24.5 |

**This paper supplies lower bounds and a refutation.  It supplies no upper
bound.**  The upper bound — the actual effectivity problem — is untouched and
remains open, exactly as `Docs/23` §23.5 documents.