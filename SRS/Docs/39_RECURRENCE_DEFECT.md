# Paper 39 — the recurrence-defect route, assessed

Version 1.0.  Tests the proposed architecture: recurrence defect `R_k` →
integer-recurrence shadowing → coefficient separation against `β_j ∉ Q(√a₁)`.

**Verdict: the architecture fails, and not only at Step 1.** Step 3–4 is
**false**, by a four-line argument (§39.3). Step 1 is separately **empirically
void** (§39.2). What survives is a genuine *dynamical reformulation* (§39.5),
which is a reformulation and not a proof.

---

## 39.1 The mechanism, and Theorem 39.1

With `η_k := x_k − m_k` (`m_k = nint(x_k)`, so `|η_k| = d_k ≤ 1/2`) and
`x_{k+2} = 2x_{k+1} + x_k` exact:

> **Theorem 39.1.**  `R_k := η_{k+2} − 2η_{k+1} − η_k = −( m_{k+2} − 2m_{k+1} − m_k ) ∈ ℤ`,
> and `|R_k| ≤ d_{k+2} + 2d_{k+1} + d_k ≤ 2`.  Moreover `R_k = 0` forces
> `m_{k+2} = 2m_{k+1} + m_k`, i.e. the nearest-integer sequence shadows the exact
> Pell recurrence at that step.

This is exact and needs **no floating point** — `R_k` is computed from the
integers `m` alone. (A float computation of `R_k` *appears* non-integral for
`k ≳ 40`; that is catastrophic cancellation, not mathematics. The engine's
`m_k` was verified against an independent 400-digit reference at
`k = 40, 100, 200, 296, 500, 1000` — all exact matches.)

## 39.2 Step 1 is empirically void: the shadowing is ubiquitous

Measured over `k = 1…298`, `a₁ = 2, a_j = 3`:

```
distribution of R_k :  { -1 : 76 ,  0 : 159 ,  1 : 62 }        (|R_k| <= 2, so no +-2)
R_k = 0 in 159 of 298 cases
R_k = R_{k+1} = 0 in 69 cases, the last at k = 296
```

So the exact recurrence `m_{k+2} = 2m_{k+1} + m_k` holds **more than half the
time**, at arbitrarily large `k`, and two consecutive instances occur
constantly.

**This kills Step 1 empirically.** Step 1 asks whether an *exceptionally small*
`d_k` forces local shadowing. But shadowing is the **generic** behaviour, not an
exception-triggered one. A condition satisfied ~53% of the time carries
essentially no information about whether `d_k` is tiny. So even if Step 1 were
proved, the shadowing it produces would not be a *rare, structured* event that
could contradict anything.

## 39.3 Step 3–4 is FALSE: there is no coefficient separation

The proposed Step 3–4 needs `|β_j − A| ≥ exp(−poly(N,k))` for the shadowing
coefficient `A`. It does not exist.

If `R_k = R_{k+1} = 0` then `m` satisfies the exact recurrence on that stretch,
so with `A` the dominant coefficient,

```
A = ( m_{k+1} − eps' m_k ) / ( eps^k (eps − eps') ) ,
A − beta_j = ( eps' eta_k − eta_{k+1} ) / ( eps^k (eps − eps') ) .           (39.1)
```

Numerator `m_{k+1} − ε'm_k ∈ Z[√a₁]` and denominator `ε^k(ε−ε') ∈ Z[√a₁]`, whose
norm is `±4a₁` (because `P_k² − a₁Q_k² = (−1)^k`). Multiplying by the conjugate
gives

> **`A ∈ Λ := (1/(4a₁)) · Z[√a₁]`.**

> **Lemma 39.2.**  `Λ` is **dense in R**, so `dist(β_j, Λ) = 0`.
>
> *Proof.*  Fix `x ∈ R`, write `4a₁x = n + θ` with `n ∈ Z`. For any `N ≥ 1`,
> density of `{q√a₁ mod 1}` — immediate from pigeonhole on
> `0, √a₁, …, N√a₁` mod 1 — gives `q` with `‖q√a₁ − θ‖ < 1/N`. Then
> `p = n − round(q√a₁ − θ)` is an integer and
> `|(p + q√a₁)/(4a₁) − x| < 1/(4a₁N)`. Let `N → ∞. ∎

`Z[√a₁]` is dense because `1` and `√a₁` generate a rank-2 subgroup of `R` with
an irrational ratio. Hence:

> **Theorem 39.3.**  No positive lower bound on `|β_j − A|` exists for
> `A ∈ Λ`.  The "algebraic separation bound" of Step 4 is **false**, not merely
> unproved. ∎

**Two false numerics, recorded as a methodology warning.**  A brute-force search
over `|p|,|q| ≤ 400` reported `δ = 6.44e−5`; a second search restricted to Pell
denominators `q = Q_k` plateaued at `0.1010`. **Both were artifacts** — a bounded
search box, and a measure-zero subset of `q`. The limit argument settles it in
four lines. A `δ` that is *small* is not a `δ` that is *positive*.

## 39.4 The quantitative version also fails — it reproduces `Docs/36` Thm 36.1

`dist(β_j, Λ) = 0` is for *unrestricted* `A`. `A` is not unrestricted: from
(39.1), `4a₁A` has integer coefficients of size `≲ ε^{2k}`, so `A` is confined to
a **finite** piece of `Λ`, and the best separation available is of order
`ε^{−2k}`.

But `ε^{−2k} ≪ ε^{−k}`, and the shadowing error in (39.1) is already
`ε^{−k}`.  So the separation bound and the shadowing bound are **the same order
and compatible** — no contradiction.  In the best case the mechanism yields

```
d_{k+1} + |eps'| d_k   >=   eps^{-k} · (1/(eps−eps')) ,
```

which is **exponent 1 — exactly `Docs/36` Theorem 36.1.**

> **Theorem 39.4.**  The recurrence-defect / shadowing / coefficient-separation
> architecture yields, at best, `Docs/36` Theorem 36.1.  It provides **no
> evasion** of `Docs/36` Theorem 36.2. ∎

## 39.5 What genuinely survives: a dynamical reformulation

`R_k ∈ {−1, 0, 1}` is a **finite-alphabet digit sequence**, and it is the
honest new object here.  Given `η_k, η_{k+1}`, the digit `ℓ_k := −R_k` is forced:
it is the unique integer bringing `2η_{k+1} + η_k + ℓ_k` back into `[−1/2, 1/2]`.
So `(η_k, η_{k+1})` is the orbit of an explicit **piecewise-linear expanding map**

```
T(u, v) = ( v , 2v + u + l(u,v) ) ,     l(u,v) = the digit keeping 2v+u+l in [-1/2,1/2]
```

on `[-1/2,1/2]²`, and `η_k = {β_j ε^k}` is a specific orbit of it.  M7.12 is
then exactly:

> Does this particular orbit enter `{ |u| < ε^{−ck} }` infinitely often?

That is a **dynamical / shrinking-target question for a piecewise-linear
expanding map**, which is a genuinely different kind of input from everything
`Docs/34`–`Docs/38` excluded — it uses the digit sequence and the map, not a
norm, a height, or a conjugate product.  It is also a well-studied subject.

**It is a reformulation, not a proof**, and no answer is claimed.

## 39.6 Status

The M7.12 attack is **closed**, now with one more mechanism excluded:

| route | excluded by | how |
|---|---|---|
| … the eleven in `Docs/38` §38.4 | `Docs/38` | — |
| recurrence defect `R_k` | `Docs/39` §39.2 | shadowing is generic (53%), not exception-triggered |
| coefficient separation at the `A`-level | `Docs/39` Thm 39.3 | **false**: `dist(β_j, Λ) = 0` |
| best case of that architecture | `Docs/39` Thm 39.4 | reduces to `Docs/36` Thm 36.1 |

**Surviving reformulation for the literature** (§39.5): the piecewise-linear
expanding map `T` on `[-1/2,1/2]²` driven by the digit `ℓ_k`, and whether the
orbit of `(η_0, η_1)` enters an exponentially shrinking target infinitely often.
That is the concrete object to raise with someone working on shrinking targets
for PL maps or on Diophantine approximation for quadratic recurrences.

**A methodological note worth keeping.**  I twice reported a small positive `δ`
that was an artifact of a bounded search, and I also predicted `ℓ_k = 0 ⟹
k = O(1)` — which is false, because a single `ℓ_k = 0` does not pin `A` and, even
with two consecutive, `δ = 0`.  The failures were all *numerical or scope*
errors, not algebra; each was caught only by an independent check.  This is the
concrete content of `Docs/37` §37.6's warning, now with a second instance.