# Paper E — The F1/Pell Family

Version 1.0. The exact separation law for the two-term Pell family, stated
and proved from the CF recurrence. This is the one family where the
separation is fully explicit.

---

## E.0 Setup

Fix `a ≥ 2`, set `α = a(a+1)`, and let `p_n/q_n` be the convergents of
`√α`. The F1 instance is

```
S_n = p_n − q_n√α,        C = (p_n, −q_n),  radicands (1, α).
```

So `m = 2`, `r = 1`, and the real input size is

```
N_n = log₂(p_n²) + log₂(α q_n²) = 2 log₂ p_n + 2 log₂ q_n + log₂ α.
```

Define the residue `D_n = p_n² − α q_n²`.

---

## E.1 Fact E.1 — the residue law

**Fact E.1.** `D_n = +1` for odd `n`, and `D_n = −a` for even `n`. Hence
`|D_n| ∈ {1, a}` and `D_n ≠ 0`.

*Proof.* The invariants of the CF satisfy

```
m_0 = 0, d_0 = 1, a_0 = a,
m_{k+1} = d_k a_k − m_k,  d_{k+1} = (α − m_{k+1}²)/d_k,
a_{k+1} = floor((a + m_{k+1})/d_{k+1}).
```

Direct two-step computation:

```
m_1 = a, d_1 = (α − a²)/1 = a, a_1 = floor(2a/a) = 2;
m_2 = a·2 − a = a, d_2 = (α − a²)/a = 1, a_2 = floor(2a/1) = 2a;
m_3 = 1·2a − a = a, d_3 = (α − a²)/1 = a, a_3 = 2.
```

Thereafter the cycle repeats: `d_k` alternates `a, 1, a, 1, …` with
`d_{2j+1} = a`, `d_{2j} = 1`, and `√α = [a; \overline{2, 2a}]`. The Pell
identity `p_n² − α q_n² = (−1)^{n+1} d_{n+1}` then gives

```
D_0 = −d_1 = −a,  D_1 = d_2 = +1,  D_2 = −d_3 = −a,  D_3 = +1, …
```

i.e. `D_n = −a` (even `n`), `+1` (odd `n`). ∎

---

## E.2 The exact law

Let `s_n = q_n√α/p_n = √(1 − D_n/p_n²)`. From `|S_n|·(p_n + q_n√α) = |D_n|`,

```
log₂(1/|S_n|) = log₂(p_n + q_n√α) − log₂|D_n|
             = log₂ p_n + log₂(1 + s_n) − log₂|D_n|.
```

Theorem of Paper E: with `N_n = 4 log₂ p_n + 2 log₂ s_n` (from `N_n =
2log₂p_n + 2log₂q_n + log₂α` and `q_n√α = p_ns_n`), the residual

```
ρ_n := log₂(1/|S_n|) − (N_n/4 + 1 − log₂|D_n|)
```

takes the **exact** form

```
ρ_n = log₂((1 + s_n)/(2√s_n)).
```

*Proof.* Substitute the two expressions and use `1 = log₂2`: the `log₂p_n`,
the `−log₂|D_n|`, and the `log₂s_n` of `N_n/4` cancel, leaving
`log₂(1+s_n) − ½log₂s_n − 1`. ∎

**Positivity.** `(1+s)/(2√s) = (√s + 1/√s)/2 ≥ 1` by AM–GM, equality iff
`s = 1`, i.e. `D_n = 0` — excluded by E.1. So `ρ_n > 0` strictly.

---

## E.3 The bound and the asymptotic

**Bound.**

```
ρ_n ≤ 0.205·a²/p_n⁴,   uniformly  ρ_n ≤ 1.31×10⁻³  (n ≥ 1).
```

*Sketch.* `f(s) := ln((1+s)/2) − ½ln s` has `f(1) = f'(1) = 0`,
`f''(1) = 1/4`, decreasing `f''` near 1, and Taylor gives
`f(s) ≤ 0.1415(s−1)²` there. Since `|1 − √(1+x)| ≤ |x|` for `x > −1` and
`s_n − 1 = O(a/p_n²)`, one gets `f(s_n) ≤ 0.1415·a²/p_n⁴`, hence
`ρ_n = f(s_n)/ln2 ≤ 0.204·a²/p_n⁴`. For `n ≥ 1`, `p_n ≥ 2a+1` and
`a²/(2a+1)⁴ ≤ 4/625`, giving `ρ_n ≤ 1.31×10⁻³`. ∎

**Asymptotic.**

```
ρ_n = D_n²/(32·ln2·p_n⁴)·(1 + O(a/p_n²)).
```

So the leading term of the residual is *second order*; the first-order pieces
in `log₂(1+s_n)` and in `N_n/4` cancel exactly. That cancellation is the whole
content of the law, and it is why the measured residual is `~10⁻⁵` rather than
`~10⁻²`.

---

## E.4 The effective rate

Substituting into the law,

```
log₂(1/|S_n|) = N_n/4 + 1 − log₂|D_n| + ρ_n,
               ∈ ( N_n/4 − log₂ a,  N_n/4 + 1 + 1.31×10⁻³ ].
```

Hence

```
E = N_n/4 + O(log a),   i.e.  E/N_n → 1/4   as n → ∞ (a fixed).
```

Under the bit-length convention `N_bits ∈ (N_n, N_n + 2]`, the same limit
holds: `E/N_bits → 1/4`. (Earlier summaries quoted "rate `1/2` under
bit-length"; that was the wrong cell. The F1 family saturates `1/4` of `N`,
and the family achieving `(r−1)/r` in general is the pigeonhole construction
of Paper B, not F1.)

---

## E.5 What this buys and what it doesn't

F1 is the unique family where the separation `E` is pinned to `N/4 + O(log a)`
with a proof and a certified residue law. It also tells us the discriminating
role of `D_n`: the sign of `S_n` flips with `n` because `D_n` alternates, and
`|S_n|` is always explained by `|D_n| ∈ {1, a}`. What F1 does **not** do is
address `r ≥ 2`: it is a one-radicand (`r = 1`) probe. The general rate
`(r−1)/r` needs the multi-radicand argument of Papers B/C, and the constant
there is open.