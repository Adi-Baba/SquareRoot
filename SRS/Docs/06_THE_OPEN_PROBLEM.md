# Paper F — The Open Problem

Version 1.0. What is actually unconjectured, unproved, and unproven. No new
theorem is proved here; the content is a precise statement of the missing
result and of the exact point where every known route fails.

---

## F.0 The conjecture

**Conjecture (SRS).** For every `S = Σ c_i √(a_i) ≠ 0` with
`|c_i| < 2^n`, `0 < a_i < 2^n`, the separation satisfies

```
log₂(1/|S|) ≤ poly(N),
```

polynomially in the total input size, with exponent independent of the number
of terms. Equivalently: there is a fixed `K` and a constant `A` (depending on
the radicand field rank at most through `r`) such that
`log₂(1/|S|) ≤ A·N^K`.

This is a uniform statement: one polynomial must work for *all* radicand
tuples `d` and all coefficient vectors `C`.

---

## F.1 What is known

**Theorem (effective, fixed rank).** For fixed `r`, `E ≤ (2^{r−1}−1)(N +
O(log r))` (Paper B). Hence every fixed rank satisfies the conjecture, and in
the band `r = O(log N)` the conjecture holds effectively (Paper B, Cor B.2).

**Theorem (Subspace Theorem).** For every `δ > 0` and every fixed radicand
tuple,

```
E ≤ (r−1+δ)/r · N + O_{d,δ}(1).
```

This is the sharp *rate*; its constant is ineffective (Paper C).

**Theorem (achieved rate).** There are instances with `E/N ≥ (r−1)/r − o(1)`
(Paper B, Thm B.5). So the rate `(r−1)/r` is sharp, and the limsup is
`limsup E/N = (r−1)/r` once the Subspace side is included.

**But the uniform conjecture is OPEN.** The Subspace-Theorem constant
`c(d,δ)` is not known to be `≥ 2^{−poly(N)}`.

---

## F.2 The exact missing statement

**Missing Theorem (FMT).** There exists `K ≥ 1` and `A > 0` such that for
every `δ ∈ (0,1)`, every rank `r`, every squarefree tuple `d =
(d_1,…,d_r)`, and every `C ∈ Z^r`,

```
|L(C)| ≥ 2^{-A·N(d,C)^{K}},
```

where `L` is the single linear form of Paper C and `N(d,C)` its real input
size.

Equivalently, the Subspace-Theorem constant satisfies
`log₂(1/c(d,δ)) ≤ A·N^K` uniformly in `d` and `C`. This is EHS Problem 33.

---

## F.3 Why `r = 2` is done but `r = 3` is only measured

At `r = 2` the inner form `C₁ + C₂√(d₂/d₁)` has *one* quadratic generator.
The separation `|C₁ + C₂α|` is controlled by the one-variable Markov constant
`λ(α)`, which is computable and positive (Paper D). Hence `r = 2` is closed.

At `r = 3` the inner form is `C₁ + C₂g₂ + C₃g₃` with two independent square
roots `g₂, g₃`. One can still write `‖C₁ + C₂g₂ + C₃g₃‖`, but the natural
analogue of `λ` is now an infimum over *two* integer variables, and there is
no Hurwitz-type universal lower bound. The best known result is the measured
finite one: `min c_S ≈ 1.7×10^{−6}` at `r = 3` over an exhaustive search, an
upper bound on the infimum only. No effective lower bound is known. The
`r = 3` case is therefore **open in the uniform sense**, even though the rate
is settled.

---

## F.4 Why `r ≥ 4` resists every known route

Every attempt either is ineffective or does not survive scrutiny:

1. **Norm / minimal polynomial (Paper B).** Gives exponent `2^{r−1} − 1`,
   exponential in `r`. Pushed to its limit it is still exponential because
   `[Q(S):Q] = 2^r` and the sign pairing is the only compression available.

2. **Squaring (`S²`).** Does not reduce the degree: the rational part
   `Σ C_i²d_i` is non-zero, so `S²` is still generic of degree `2^r`. The
   complement pairing is destroyed, and the bound does not improve.

3. **Geometric-progression hope.** One could try to pick a single
   `γ` capturing all three generators at `r = 3`. That requires
   `a₃/a₁ = (a₂/a₁)²`, forcing `d₃ = d₂` — a contradiction — or
   `a₃/a₁ = (a₂/a₁)^{−1}`, forcing `d₂d₃ = d₁²`. So the one algebraic trick
   that would have collapsed `r = 3` to `r = 2` is excluded.

4. **Subspace Theorem (Paper C).** Always available, always ineffective in the
   constant.

The pattern is not a coincidence. The obstruction is genuine: either the
constant is non-uniform in `d`, or the argument is exponential in `r`. No
third option is known.

---

## F.5 The retracted `r ≤ 4` claim

An earlier revision of this repository claimed an effective general bound for
`2 ≤ r ≤ 4`,

```
|S| ≥ sqrt(D) / ((r+1)M)^{2^{r−2}−1},
```

resting on the identity `dist(aW, ℤ) = a·dist(W, ℤ)` for a lattice `W`. That
identity is false unless `aℤ = ℤ`, i.e. unless `a = ±1`. The claim is
therefore **retracted**. The retraction is pinned by a regression test that
exhibits a counterexample to the claimed inequality at `r = 2`
(`d = (2,3)`, `C = (1,−1)` gives `|S| = |√2 − √3| = 0.3178`, while the claimed
bound is `sqrt(3) = 1.7321`).

---

## F.6 What would settle it

A proof of FMT, i.e. an *effective* lower bound

```
|L(C)| ≥ 2^{−poly(N)}
```

uniform in `d` and `C`. Anything short of that is either an ineffective
existence theorem or a non-uniform bound on a measured set, and does not count.
The repository deliberately states this as the only open statement; nothing
else is claimed.