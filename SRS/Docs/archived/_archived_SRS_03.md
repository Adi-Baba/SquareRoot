# SRS — Theorem 3: The F1 Separation Law, Exactly

Version 1.0. This document states and proves Theorem 3 — the exact separation
law for the F1 (Pell) family — from first principles, including the residue
law (Fact 3.1) on which the whole statement rests. It replaces the float-based
tolerance in `f1_theorem.py` (`_SEPARATION_TOLERANCE = 1e-3`) with a proved
statement, and proves that tolerance is a theorem's slack, not an empirical
margin.

All bounds use the bit-length-free real size `N_n` of §0 in
`SRS_00_Definitions_and_Notation.md`.

---

## 3.1 The F1 family

Fix an integer `a ≥ 2` and set

```
α = a(a+1),        S_n = p_n − q_n·sqrt(α),
```

where `p_n/q_n` are the continued-fraction convergents of `sqrt(α)`. The SRS
instance at index `n` is the two-term sum `p_n·1 − q_n·sqrt(α)` (coefficients
`p_n, −q_n`, radicands `1, α`); it has `m = 2`, `r = 1`, and input size

```
N_n = log₂(p_n²) + log₂(α q_n²) = 2log₂p_n + 2log₂q_n + log₂α.
```

Define the residue

```
D_n = p_n² − α q_n².
```

---

## 3.2 Fact 3.1 (residue law), proved

**Claim.** `D_n = +1` when `n` is odd, and `D_n = −a` when `n` is even. In
particular `|D_n| ≤ a` and `D_n ≠ 0` for every `n ≥ 0`.

*Proof.* Write the continued fraction `sqrt(α) = [a_0; a_1, a_2, …]` and the
associated invariants `(m_k, d_k, a_k)` with

```
m_0 = 0, d_0 = 1, a_0 = floor(sqrt(α)) = a,
```

`m_{k+1} = d_k a_k − m_k`, `d_{k+1} = (α − m_{k+1}²)/d_k`,
`a_{k+1} = floor((a + m_{k+1})/d_{k+1})`.

A direct computation of the first few complete quotients gives

```
m_1 = a,   d_1 = (α − a²)/1 = a,   a_1 = floor((a + a)/a) = 2,

m_2 = a·2 − a = a,   d_2 = (α − a²)/a = 1,   a_2 = floor((a + a)/1) = 2a,

m_3 = 1·2a − a = a,   d_3 = (α − a²)/1 = a,   a_3 = floor((a + a)/a) = 2.
```

The cycle then repeats, so `(m_k, d_k)` alternates: `d_k = a, 1, a, 1, …`
with `d_{2j+1} = a`, `d_{2j} = 1`, and `(a_k)` alternates `a_1 = 2,
a_2 = 2a, a_3 = 2, …`. Hence `sqrt(α) = [a; over{2, 2a}]`. The point for
us is the alternation of `d_k`.

The Pell identity is `p_n² − α q_n² = (−1)^{n+1} d_{n+1}`. With `d_{n+1}`
alternating `a, 1, a, 1, …` (starting at `d_1 = a`), we get

```
D_0 = (−1)^1 d_1 = −a,
D_1 = (−1)^2 d_2 = +1,
D_2 = (−1)^3 d_3 = −a,
D_3 = (−1)^4 d_4 = +1,  …
```

i.e. `D_n = −a` for even `n`, `D_n = +1` for odd `n`. ∎

*Remark on the period.* For `a = 2`, `α = 6`, this gives
`sqrt(6) = [2; overline{2,4}]`, matching the classical expansion. The residue
table `D_n = −2, +1, −2, +1, …` follows.

---

## 3.3 Theorem 3.2 — the exact separation law

Let

```
s_n = q_n·sqrt(α)/p_n  =  sqrt(1 − D_n / p_n²)   (a positive real),
```

and define the residual

```
ρ_n := log₂(1/|S_n|) − (N_n/4 + 1 − log₂|D_n|).
```

**Theorem 3.2.** For every `a ≥ 2` and `n ≥ 0`:

**(i)** *(exact)* `N_n / 4 = log₂ p_n + (1/2)·log₂ s_n`.

**(ii)** *(exact, one-sided)*  
`ρ_n = log₂((1 + s_n) / (2·sqrt(s_n))) > 0`.

**(iii)** *(bound)*  
`ρ_n ≤ 0.205 · a² / p_n⁴`, and uniformly `ρ_n ≤ 1.31×10⁻³` for `n ≥ 1`.

**(iv)** *(asymptotic)*  
`ρ_n = D_n² / (32·ln 2·p_n⁴) · (1 + O(a / p_n²))`.

*Proof.*

**(i).** `2 log₂ q_n + log₂ α = 2 log₂(q_n·sqrt(α)) = 2 log₂(p_n·s_n)
= 2 log₂ p_n + 2 log₂ s_n`. Hence

```
N_n = 2 log₂ p_n + (2 log₂ p_n + 2 log₂ s_n) = 4 log₂ p_n + 2 log₂ s_n,
```

and dividing by `4` gives `(i)`.

**(ii).** Rationalise using the residue:

```
|S_n|·(p_n + q_n·sqrt(α)) = |D_n|,
```

so

```
log₂(1/|S_n|) = log₂(p_n + q_n·sqrt(α)) − log₂|D_n|
             = log₂ p_n + log₂(1 + s_n) − log₂|D_n|.
```

Subtract `(N_n/4 + 1 − log₂|D_n|)` and use `(i)` and `1 = log₂ 2`. The
`log₂ p_n`, the `−log₂|D_n|`, and the `2log₂s_n`-part of `N_n/4` reduce to

```
ρ_n = log₂(1 + s_n) − (1/2) log₂ s_n − 1 = log₂((1 + s_n)/(2·sqrt(s_n))).
```

For the sign: `(1+s)/(2√s) = (√s + 1/√s)/2 ≥ 1` by AM–GM, equality iff
`s = 1`, i.e. `D_n = 0` — impossible by Fact 3.1. Hence `ρ_n > 0`.  Note
`f(s) = ln((1+s)/2) − (1/2) ln s` satisfies `f(s) = f(1/s)`, so the sign
statement is parity-independent.

**(iii).** Let `f(s) = ln((1+s)/2) − (1/2) ln s`, so `ρ_n = f(s_n)/ln 2`.
Then `f(1) = 0`, `f'(1) = 1/2 − 1/2 = 0`, and

```
f''(s) = 1/(2s²) − 1/(1+s)²,   f''(1) = 1/4.
```

Since `|D_n| ≤ a`, `s_n = sqrt(1 − D_n/p_n²)` lies in
`[sqrt(1 − a/p_n²), sqrt(1 + a/p_n²)]`. For `n ≥ 1`, `p_n ≥ p_1 = 2a+1`, and
`a/(2a+1)²` is decreasing in `a ≥ 2`, so `a/p_n² ≤ a/(2a+1)² ≤ 2/25` (at
`a = 2`). Hence `s_n ∈ [sqrt(23/25), sqrt(27/25)] ⊂ [0.959, 1.041]`.

On this interval `f''` is decreasing (its derivative is negative for
`s < 1/(2^{1/3}−1) ≈ 3.85`), so `f''(s) ≤ f''(sqrt(23/25)) ≈ 0.283`. Taylor:

```
f(s) = (1/2) f''(ξ)(s−1)² ≤ 0.1415·(s−1)².
```

The mean-value step `|sqrt(1+x) − 1| ≤ |x|` for `x > −1` gives
`(s_n − 1)² ≤ (D_n/p_n²)² ≤ a²/p_n⁴`. In `log₂` units,

```
ρ_n ≤ (0.1415 / ln 2)·(a²/p_n⁴) ≈ 0.204·a²/p_n⁴ ≤ 0.205·a²/p_n⁴.
```

For the uniform bound, `a²/p_n⁴ ≤ a²/(2a+1)⁴`, and `a²/(2a+1)⁴` is decreasing
in `a` (derivative sign `(2a+1)^{-5}·2a(1−2a) < 0`), so its maximum over
`a ≥ 2` is at `a = 2`: `4/625`. Hence

```
ρ_n ≤ 0.205·(4/625) = 1.312×10⁻³.
```

**(iv).** Since `f(1) = f'(1) = 0` and `f''(1)/2 = 1/8`,

```
f(s) = (s − 1)²/8 + O((s−1)³).
```

Also `s_n − 1 = sqrt(1 − D_n/p_n²) − 1 = −D_n/(2p_n²) + O(a²/p_n⁴)`.
Squaring, `(s_n−1)² = D_n²/(4p_n⁴) + O(a³/p_n⁶)`, and dividing by `8·ln 2`:

```
ρ_n = D_n²/(32·ln 2·p_n⁴)·(1 + O(a/p_n²)).   ∎
```

---

## 3.4 Corollary 3.3 (the Θ(N) law, with constants)

For `a ≥ 2`, `n ≥ 1`:

```
N_n/4 + 1 − log₂|D_n|  <  log₂(1/|S_n|)  ≤  N_n/4 + 1 − log₂|D_n| + 1.31×10⁻³,
```

so, with `θ_n := log₂(1/|S_n|) − N_n/4` = `1 − log₂|D_n| + ρ_n`,

```
1 − log₂ a  <  θ_n  ≤  1 + 1.31×10⁻³   (even n sharper: θ_n ≤ 1 − log₂a + 1.31×10⁻³).
```

In particular `log₂(1/|S_n|) = N_n/4 + O(log a)`, i.e. the F1 family achieves
`E = Θ(N)` with an explicit constant `1/4` in front of `N`.  The one-sidedness
`ρ_n > 0` says the guaranteed rate is never *exceeded* in the other direction.

---

## 3.5 Why the residual is second-order

Both sides of the law carry a first-order correction `−D_n/(4p_n²·ln 2)`:
`log₂(1+s_n)` on the left and `(1/2)log₂ s_n` inside `N_n/4` on the right.
Theorem 3.2(ii) shows these cancel **exactly**, leaving the second-order term
`D_n²/(32·ln 2·p_n⁴)`.  This is why the measured worst residual (`7.5×10⁻⁵`
at `a = 2, n = 1`, over `a ≤ 24`, `n ≤ 39`) is roughly 400× smaller than
either correction alone (`~0.03`), and why part (iv)'s leading term predicts
`1/(32·ln 2·5⁴) = 7.2×10⁻⁵` there — within 5% of the measurement.  A float
that forms `q_n·sqrt(α)/p_n` and then subtracts from `1` reads the residual as
exactly `0` past `n ≈ 10`; that is a numerical artefact, and is why the proof
in this file never forms `s_n − 1` directly — it factors through
`(1−sqrt(s_n)) = (1−s_n)/(1+sqrt(s_n))` with `1 − s_n = D_n/(p_n(p_n +
q_n sqrt α))`.

---

## 3.6 `N` conventions

`N_n` is the real quantity `2log₂p_n + 2log₂q_n + log₂α`.  The instance
generator's `N_bits = bit_length(p_n²) + bit_length(α q_n²)` satisfies
`N_bits ∈ (N_n, N_n + 2]`.  Every ratio statement transfers with that `O(1)`
slack.

---

## 3.7 What is not claimed

- No uniform bound over all radicands — this is a one-family (two radicands,
  `r = 1`) statement.
- The constant `1.31×10⁻³` is proved for the `N_n` quantity; rederiving it for
  `N_bits` changes it by the additive slack `O(1)/N`.