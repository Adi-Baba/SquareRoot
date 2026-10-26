# Paper B — The Effective Fixed-Rank Norm Bound

Version 1.0. The effective bound that closes every *fixed* rank. The proof
is a classical norm/trace argument, done here in full, and we identify the
two places where it cannot be pushed to a uniform poly(r) bound — that is
Paper F.

---

## B.0 Setup and notation

`S = Σ_{i=1}^m C_i·sqrt(d_i)` with `d_i` distinct squarefree, `C_i`
non-zero integers, `r = log₂[K:Q]`, `K = Q(sqrt(d₁),…,sqrt(d_m))`. Put

```Markdown
M := max( ‖C‖_∞ , max_i d_i ),      m is the number of terms,
N := Σ_i (log₂|C_i| + log₂ d_i)     (real size convention).
```

---

## B.1 The norm argument, exactly

**Theorem B.1.** For `S ≠ 0`, with **no** extra hypothesis,

```
log₂(1/|S|) ≤ (2^r − 1)·( log₂ m + (3/2)·log₂ M ).
```

Moreover, the exponent drops to `2^{r−1} − 1` exactly when the radicands
admit the **all-ones sign automorphism** — equivalently, when every
`F₂`-dependency among `{v(d_i)}` has **even support**. In particular this
always holds when the `d_i` are `F₂`-independent, and then

```
log₂(1/|S|) ≤ (2^{r−1} − 1)·( log₂ m + (3/2)·log₂ M ).
```

*Proof.* Every conjugate of `S` over `Q` has the form
`S^σ = Σ_i C_i ε_i sqrt(d_i)` with `ε_i ∈ {±1}`, and there are at most `2^r`
such conjugates (exactly `2^r` for F₂-independent radicands). A naive
estimate is

```
|S^σ| ≤ Σ_i |C_i| sqrt(d_i) ≤ m · M · sqrt(M) = m·M^{3/2}.
```

`S` is an algebraic integer, so `N_{K/Q}(S) = Π_σ S^σ` is a non-zero
integer, hence `|N_{K/Q}(S)| ≥ 1`. One factor is `S` itself; the other
`2^r − 1` factors are each `≤ m·M^{3/2}`. Therefore

```
1 ≤ |S| · (m·M^{3/2})^{2^r − 1},
```

so `|S| ≥ (m·M^{3/2})^{-(2^r − 1)}`, i.e.

```
log₂(1/|S|) ≤ (2^r − 1)·(log₂ m + (3/2)·log₂ M).
```

That is the unconditional bound. For the conditional improvement: the all-ones
sign map `τ : sqrt(d_i) ↦ −sqrt(d_i)` is an automorphism of `K` iff assigning
`χ(e_i)=1` is a well-defined character of the radicand lattice, which holds
iff every `F₂`-dependency of the `v(d_i)` has even support. (If a dependency
has odd support — e.g. `v(2)+v(3)+v(6)=0` on three radicands — then
`τ(√6)=τ(√2)τ(√3)=√6 ≠ −√6`, so no such automorphism exists.) When it does
exist, `τ(S^σ) = −S^σ`, the conjugates pair as `{T, −T}`, and
`N_{K/Q}(S) = ±Π_{pairs} T²`. Isolating the pair `{S, −S}`:

```
1 ≤ |N| = |S|² · Π_{other pairs} |T|² ≤ |S|²·(m·M^{3/2})^{2(2^{r−1}−1)},
```

so `|S| ≥ (m·M^{3/2})^{-(2^{r−1}−1)}`, i.e. the conditional bound. ∎

---

## B.2 Why `N` linearizes

`N ≥ m` (each `d_i ≥ 2` contributes at least one bit and each `|C_i| ≥ 1`
at least one, so `N ≥ 2m ≥ m`), hence `log₂ m ≤ log₂ N`. Also `M ≤ 2^N`
because `M` is one of the numbers `|C_i|, d_i`, each `≤ 2^N`. And the rank
satisfies `r ≤ m`. Collecting:

Consequently the unconditional form is

```
log₂(1/|S|) ≤ (2^r − 1)·( log₂ N + (3/2)·N ),
```

and, whenever the all-ones sign automorphism exists (in particular for
F₂-independent radicands), the improved form is

```
log₂(1/|S|) ≤ (2^{r−1} − 1)·( log₂ N + (3/2)·N ).
```

So for fixed `r` the bound is `E ≤ (2^{r−1} − 1)·((3/2)N + log₂N)` (resp.
`2^r − 1` without the pairing), i.e. `E/N ≤ (3/2)(2^{r−1}−1) + o(1)` in the
paired case. (An earlier version of this file wrote `(2^{r−1} − 1)(N + O(log r))`;
that dropped a factor `3/2` on the main term, because `M^{3/2}`, not `M`,
bounds each non-trivial conjugate.) Either way the constant is exponential in
`r`; the exponent is a **constant** at fixed rank, so every fixed rank is
settled.

---

## B.3 From fixed rank to `r = O(log N)`

The bound in B.1 is the only effective poly-in-`N` statement we have. To reach
`E/N ≤ (r−1)/r` we need the Subspace Theorem (`Paper C`). But we can combine:
if `2^r ≤ N^K` for some fixed `K` (i.e. `r ≤ K·log₂N`), then the norm bound
itself gives `E ≤ (2^r)((3/2)N + log N) = O(N^{K+1})`, which is polynomial.
Hence:

**Corollary B.2.** For `r ≤ c·log₂ N` with `c` fixed, `E ≤ poly(N)`
effectively. So the conjecture is **proved** in the band `r = O(log N)`.

The open band is therefore `O(log N) ≤ r ≤ N/2` — and there the norm bound
is exponential in `r` while the Subspace bound is ineffective. Paper C
explains the complement.

---

## B.4 Why this cannot be improved to poly(r)

Two hard facts, proved in Paper F:

1. `[Q(S):Q] = 2^r` (Paper A, Prop A.4). The norm has `2^r` conjugates. Any
   argument that bounds `|S|` by comparing it to its full norm necessarily
   multiplies `2^r` factors, each `≤ poly(M)`. Squeezing to `2^{r−1}` is
   already done by the sign pairing; no further pairing collapses distinct
   conjugates (the only fixed value of the global sign flip on `Q(S)` is `0`,
   and no other automorphism has fewer than `2^r` distinct images on generic
   `S`).
2. Squaring `S` does **not** reduce the count: `S²` lies in the same degree-`2^r`
   field (it is not rational — the rational part `Σ C_i²d_i ≠ 0`), and its
   conjugates also number `2^r`. So degree-halving tricks fail.

These are not gaps in the proof; they are the reason the exponent in this
argument is at best `2^{r−1} − 1`, and is genuinely `2^r − 1` whenever the
radicands lack the all-ones automorphism. The only routes left are
(a) use the Subspace Theorem (Paper C), which is ineffective, or (b) find a
different structure. This is the wall.

---

## B.5 The achieved rate (pigeonhole), for the matching lower bound

**Theorem B.5.** Let `d₁, …, d_r` be distinct squarefree integers and `B ≥ 1`.
Then there are integers `C₁, …, C_r`, not all zero, `|C_i| ≤ B`, with

```Markdown
|Σ_i C_i sqrt(d_i)| ≤ (Σ_i sqrt(d_i)) / B^{r−1}.
```

*Proof.* The `(B+1)^r` values `v(C) = Σ C_i sqrt(d_i)`, `C ∈ {0,…,B}^r`, are
pairwise distinct (an equality is a rational dependence among the
`sqrt(d_i)`, forbidden by Paper A, Thm A.3). They all lie in
`[0, B·Σ_i sqrt(d_i)]`. Sorting them, two consecutive values differ by at
most `B·Σ sqrt(d_i) / ((B+1)^r − 1) ≤ 2(Σ sqrt(d_i))·B^{-(r−1)}`, and their
difference is the desired sum. ∎

So at fixed radicands with `|C_i| ≤ B`, `N ≍ r·log₂B` and the construction
gives `E ≥ (r−1)log₂B − O(1)`, i.e. `E/N ≥ (r−1)/r − o(1)`. Together with
the Subspace-Theorem upper bound (Paper C) this makes the limsup sharp:
`limsup E/N = (r−1)/r`.
