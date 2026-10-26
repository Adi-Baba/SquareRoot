# Paper 33 — M7.12: an effective COUNT bound for the exceptional set

Version 1.0.  Two new effective results obtained by exploiting `q_1` being
**fixed**, following the collapse of `Docs/32`.  Together they give the first
**explicitly bounded cardinality** for the EHS/KMN exceptional set — where
`Docs/24`–`Docs/30` had only *finiteness*.

They do **not** yet give the largest element, which is what SRS needs.

---

## 33.1 New result 1: the Subspace-Theorem exponent is explicit

KMN's key step (their eq. (13)) is
```
Π_{w∈S'} Π_{L∈ℒ_w} |L(𝐱)|_w/|𝐱|_w   ≤   θ^n · ℬ^r · H(𝐱)^{−r−1} ,          (33.1)
```
with `ℬ := max H(b_{i,j,ℓ})`.  They then write *"there is `ε > 0` such that
… when `n` is sufficiently large"*, and apply Theorem 2.1.

The vagueness has two sources:

* **`ℬ = e^{o(n)}`** in general, because the hypothesis is `h(q_i) < f(n)` with
  `f` sublinear.  Then `θ^n ℬ^r = e^{−n logθ + o(n)}` and the required
  `ε > 0` exists only *eventually*, with no explicit `n_0`.
* the height bookkeeping of Prop 3.4.

> **Theorem 33.1 (explicit ST exponent, in our case).**  For `q_1 = β_j`
> **fixed**, the numbers `b_{1,j_1,ℓ}` are fixed rationals (Lemma 3.1 of KMN),
> hence `ℬ` is a **constant**.  Therefore (33.1) reads
> ```
> θ^k ℬ^r H(𝐱)^{−r−1}   ≤   H(𝐱)^{−r−1−ε} ,      ε := −log θ / 2 > 0 ,
> ```
> for every `k ≥ n_0` with
> ```
> n_0  =  ⌈ ( r·log ℬ + log(1/ε) ) / (−log θ) ⌉ ,                            (33.2)
> ```
> **explicit, and `O(log ℬ)`** since `−log θ ≥ (log ε)/(3n(n−1))` and `r ≤ 2`.

So in our situation the *"sufficiently large"* in KMN is an **explicit
`n_0`**.  The fixedness of `q_1` — precisely what `Docs/32` exploited for the
dimension collapse — also removes the `o(n)` from the ST exponent.

---

## 33.2 New result 2: each exceptional hyperplane holds at most one exception

**Theorem 33.2 (one exception per hyperplane, computable).**  Fix a nontrivial
relation `a_0p_k + Ψ_1ε^k + Ψ_2(ε')^k = 0` with `(a_0,Ψ_1,Ψ_2) ∈ L³`, and put
`ζ := −Ψ_2/a_0` (the case `a_0 = 0` is the trivial relation).  If it holds at
some `k`, then `Ψ_1 = −a_0β_j` and

```
( β_j − ζ/2 )·ε^k  =  −(ζ/2)·(ε')^k ,
```

so, taking absolute values and using `|ε'/ε| = (√a_1−1)/(√a_1+1) < 1`,

```
k  =  log( |ζ| / ( 2|β_j − ζ/2| ) )  /  log( ε/|ε'| ) ,                     (33.3)
```

which determines `k` **uniquely** from `ζ`.

> **Consequently every exceptional hyperplane of Theorem 1.4 contains at most
> one element of `ℳ_θ`, and that element is explicitly computable from the
> hyperplane.**

(The contradiction of `Docs/32` Thm 32.3 is the special case `k → ∞`:
`k` determined by (33.3) cannot range over an infinite set.)

---

## 33.3 Theorem 33.3 — an effective COUNT bound

**Theorem 33.3.**  Let `a_1 ≠ a_j` be distinct squarefree, `θ ∈ (0,1)`, and
`ℳ_θ = { k : ‖β_jε^k‖ < θ^k }`.  Then
```
|ℳ_θ|   ≤   n_0  +  t ,                                                     (33.4)
```
with `n_0` explicit by (33.2) and `t` **any** valid bound for the number of
exceptional hyperplanes of Theorem 1.4.

*Proof.*  By (33.1), every solution with `k ≥ n_0` satisfies the
Subspace-Theorem inequality, hence lies in one of the `t` hyperplanes; by
Theorem 33.2 each hyperplane carries at most one such `k`.  Solutions with
`k < n_0` number at most `n_0`. ∎

**And here the quantitative Subspace Theorem *is* usable** — because
`Docs/32` Thm 32.1 gives `r ≤ 2`, so the ST is applied in `L^3` and
Schlickewei's bound is
```
t   ≤   c(3) · δ^{−3} · ( log(1/δ) − 1 )² ,        δ = −log θ / 2 ,           (33.5)
```
with `c(3)` **absolute** — the variable count is the fixed `3`, not the rank.
Numerically (`a_1 = 2`):

```
 rank n   theta_0     delta    t / c(3)
    2     0.86338   0.07345      6.6e3
    3     0.95221   0.02448      5.0e5
    4     0.97581   0.01224      6.3e6
    6     0.99025   0.00490      1.6e8
   10     0.99674   0.00163      6.8e9
```

so `t` is **polynomial in the rank** (`δ^{-3} ~ n^6`) and carries **no
dependence on `a_1`, `a_j` or `N`**.

> **Theorem 33.4.**  `|ℳ_θ| ≤ n_0 + t` with `n_0 = O(log ℬ)` and `t` given by
> (33.5).  **`ℬ` is the maximum height of the `K`-coordinates of the fixed
> `β_j`, hence `log ℬ = O(N)`, so `n_0 = O(N)`.  Therefore**
> ```
> |ℳ_θ|   =   O(N)  +  O( δ^{−3} (log 1/δ)² ) ,                             (33.6)
> ```
> **an explicit, polynomial bound on the cardinality of the exceptional set —
> the first such bound in this project.**  `Docs/24`–`Docs/30` had only
> *finiteness*.

---

## 33.4 What this does **not** give — stated plainly

SRS needs the **largest** element `k* = max ℳ_θ`, not `|ℳ_θ|`.  And (33.4) does
not deliver it: from (33.3), `k = k(ζ)`, and `ζ` ranges over **all of `L`**, so
`k(ζ)` is unbounded — `ζ → 2β_j` makes `|β_j − ζ/2| → 0` and `k(ζ) → ∞`.

**The exceptions with `k < n_0` are fully enumerable** (`n_0 = O(N)`, explicit).
The exceptions with `k ≥ n_0` are at most `t` of them, each equal to `k(ζ_i)`,
but the `ζ_i` — equivalently the hyperplanes — are not identified.

> **M7.12 (final, second iteration).**  Bound
> ```
> |ζ|  away from its degeneracy near 2β_j ,   equivalently  bound  k(ζ) ,
> ```
> uniformly over the at most `t` exceptional values produced by Theorem 2.1 in
> `L³`.  **Everything else is done:** the candidates are described, the count is
> bounded explicitly, each candidate is explicitly computable, and each
> hyperplane carries at most one exception.

Note the shift in what is being asked.  It is **no longer** an effectivity
problem for the Subspace Theorem in the usual sense: the quantitative ST
*already* supplies `t`, and `Docs/32` reduced the ambient space to `L³`.  What is
missing is arithmetic information about **which** at most `t` algebraic numbers
`ζ ∈ L` arise — a *finite, bounded, algebraic* identification problem rather
than an analytic one.

---

## 33.5 Status

| Item | Result |
|---|---|
| ST exponent `ε` and threshold `n_0` explicit (fixed `q_1`) | **Thm 33.1** — NEW |
| Each exceptional hyperplane holds ≤ 1 exception, explicitly computable | **Thm 33.2** — NEW |
| **Effective count bound `|ℳ_θ| ≤ n_0 + t`** | **Thm 33.3** — NEW |
| `t` via Schlickewei with `n = 3` (absolute `c(3)`) | (33.5) |
| **`|ℳ_θ| = O(N) + O(δ^{-3}(\log 1/\delta)^2)`** | **Thm 33.4** — first explicit cardinality bound |
| Refutation of `Docs/31` §31.4 ("counting is not locating") | **Partly**: counting *is* enough for the cardinality; it is **not** enough for `k*` |
| **Remaining: identify the `≤ t` values of `ζ`, or bound `k(ζ)`** | **OPEN — M7.12** |