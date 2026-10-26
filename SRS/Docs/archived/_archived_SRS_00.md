# SRS — Definitions and Notation

Version 1.0. This document establishes the exact objects, conventions, and
standing assumptions used throughout the SRS proof papers. Every theorem in
`SRS_02_Main_Theorems.md` and `SRS_03_F1_Separation_Law.md` refers to this
notation.

---

## 0.1 The Square-Root Sum

We are given `k ≥ 1` pairs of integers `(c_j, a_j)` with `c_j ≠ 0`,
`a_j ≥ 1`. The **square-root sum** is

```
S = Σ_{j=1..k} c_j · sqrt(a_j).
```

The sum is over a multiset of radicands; we allow `a_i = a_j` and coefficients
with any sign. The standard convention is `sqrt(1) = 1`, so `c_j sqrt(1)` is
an integer term. We **exclude** the trivial case `S = 0` from all separation
bounds — the zero test is handled separately (see `SRS_01_Zero_Test_Correctness.md`).

---

## 0.2 Squarefree Decomposition and the Aggregated Form

Every positive integer `a` factors uniquely as `a = t² · d` with `d`
squarefree (`d = 1` is squarefree). We call `d` the **squarefree part** of `a`
and `t` its **coefficient multiplier**.

Applying this to each radicand `a_j = t_j² · d_j` gives

```
S = Σ_{j=1..k} (c_j · t_j) · sqrt(d_j).
```

Since `d_j` are squarefree, we can **aggregate** by collecting terms with the
same squarefree part. Let `D = {d_1, …, d_m}` be the set of **distinct**
squarefree parts appearing, with `d_i ≠ 1` (integer terms with `d=1` are kept
separate). Then

```
S = C_0 + Σ_{i=1..m} C_i · sqrt(d_i),
```

where `C_i = Σ_{j: d_j = d_i} c_j · t_j` are non-zero integers, and `C_0` is the
sum of all `c_j` with `d_j = 1`. For the separation bounds the integer term
`C_0` can be absorbed into the smallest `C_i` or handled as an extra term
(`d = 1` behaves like `sqrt(1) = 1`), so we usually state results with
`C_i ≠ 0` and `d_i` distinct squarefree integers `≥ 2`, and note that `C_0`
only adds `O(1)` to `N`.

**Convention.** In the theorems below, `m` denotes the number of non-rational
squarefree parts, `d_1, …, d_m` are distinct squarefree integers `≥ 2`,
and `C_1, …, C_m` are non-zero integers. We write `S = Σ_{i=1..m} C_i sqrt(d_i)`
and treat the integer part as a constant additive shift to `N`.

---

## 0.3 Input Size `N`

There are two natural conventions for the input size. We use **both** and
explicitly state which one a theorem uses.

### 0.3.1 The Bit-Length Convention (`N_bits`)

```
N_bits = Σ_{i=1..m} ( bit_length(|C_i|) + bit_length(d_i) ),
```

where `bit_length(x)` is the number of bits in the binary representation of the
positive integer `x` (i.e. `floor(log₂ x) + 1` for `x ≥ 1`). This is the
**exact integer** input size used by the Python implementation. For a single
term `C sqrt(d)` we have `bit_length(C) ∈ (log₂|C|, log₂|C| + 1]` and
`bit_length(d) ∈ (log₂ d, log₂ d + 1]`, so

```
N_bits ∈ ( N_real, N_real + 2m ],
```

where `N_real` is the real-valued sum defined next. Every ratio statement
transfers between `N_bits` and `N_real` with additive slack `O(m)`; no theorem
depends on which is used beyond that slack.

### 0.3.2 The Real-Valued Convention (`N_real`)

```
N_real = Σ_{i=1..m} ( log₂|C_i| + log₂ d_i ).
```

This is the quantity that appears naturally in the analytic proofs. It differs
from `N_bits` by at most `2m` and is sometimes more convenient because it
linearises products. Unless explicitly noted, `N` in theorems means `N_real`.

---

## 0.4 The Exponent `E` and the Separation Rate `E/N`

```
E = log₂( 1 / |S| ),
```

so `E > 0` exactly when `|S| < 1`. The **separation rate** is `E / N`. The
central conjecture of SRS is that `E/N` is bounded by a polynomial in `N`
(and hence in `r` and `log N`).

---

## 0.5 F₂-Rank `r`

For a squarefree integer `d ≥ 2`, its **parity vector** is the map
`p → exponent of p in d (mod 2)`, i.e. the vector of prime exponents modulo 2.
The **F₂-rank** of a set of distinct squarefree integers `{d_1, …, d_m}` is the
dimension over `F₂` of the span of their parity vectors. Equivalently:

```
r = log₂ [ Q(sqrt(d_1), …, sqrt(d_m)) : Q ].
```

The field extension degree is exactly `2^r`. The rank satisfies `1 ≤ r ≤ m`,
and `r = m` iff the `d_i` are F₂-independent (no nonempty subproduct is a
perfect square).

**Important.** The rank `r` is the number of *independent* quadratic generators
in the field. It is not the number of terms `m`. A term with `d = 1` (rational
part) does not increase `r`.

---

## 0.6 Distance and Norm Notation

For a real number `x`, `‖x‖ = min_{z ∈ Z} |x − z|` is the distance to the
nearest integer. Note `‖x‖ ∈ [0, 1/2]`.

For a vector `v ∈ R^m`, `‖v‖_∞ = max_i |v_i|`.

The **Markov constant** of an irrational `α` is the badly-approximable constant

```
λ(α) = inf_{q ≥ 1} q · ‖qα‖.
```

`λ(α) > 0` iff `α` is badly approximable (bounded partial quotients). For
quadratic irrationals this holds by Lagrange's theorem.

---

## 0.7 The Two Conventions for `r=2` Markov Constants

Two different Markov-type constants appear in the literature and in this
repository. We distinguish them carefully:

- **Optimal constant (infimum):**  
  `λ(α) = inf_{q ≥ 1} q ‖qα‖`.  
  This is the *tightest* constant `c₀` such that `|qα − p| ≥ c₀/q` for all
  integers `p, q (q ≠ 0)`. By Hurwitz, `λ(α) ≤ 1/√5` for every irrational `α`.
  The corresponding separation constant is `1/λ(α) ≥ √5`.

- **Repo's constant (supremum over convergents):**  
  `M(α) = sup_{n ≥ 0} q_n ‖q_n α‖`, where `p_n/q_n` are the **reduced**
  continued-fraction convergents of `α`.  
  It is always true that `M(α) ≤ 1` (because `|q_n α − p_n| < 1/q_{n+1}`),
  so `1/M(α) ≥ 1`. However `1/M(α)` is **not** a valid uniform lower bound for
  `|qα − p|` over all `p, q` — it can exceed the optimal constant and is not
  admissible for the separation theorem.

The theorems below state their bounds in terms of `λ(α)` (the infimum). The
convention `c = 1/M(α)` is recorded only for historical reference and is
**not** used in the proofs.

---

## 0.8 `r=2` Convergent Witnesses (Hand-Checkable)

For a quadratic irrational `α = sqrt(D)` with `D` squarefree and not a square,
the `0`-th convergent is `p_0/q_0 = floor(sqrt(D))/1`. Its value is

```
q_0 ‖q_0 α‖ = |sqrt(D) − floor(sqrt(D))| = {sqrt(D)} < 1,
```

which is a trivial and hand-checkable witness that `M(α) ≤ 1`. No deeper
convergents are needed to see that the repo's `M` is bounded by `1`.

---

## 0.9 Assumptions on the Radicands

Unless stated otherwise, the `d_i` are assumed to be **distinct squarefree
integers `≥ 2`**, and the `C_i` are **non-zero integers**. The integer term
`C_0` is omitted from the main formulas but its presence is noted in the
bounds as an additive `O(1)` to `N` and `E`.

---

## 0.10 Quick Reference

| Symbol | Meaning |
|--------|---------|
| `S` | `Σ C_i sqrt(d_i)`, square-root sum |
| `d_i` | distinct squarefree `≥ 2` |
| `C_i` | non-zero integer coefficients |
| `m` | number of non-rational squarefree parts |
| `r` | F₂-rank of `{d_i}` (`= log₂ [Q(sqrt(d_i)):Q]`) |
| `N` | input size (`N_real` or `N_bits`, stated per theorem) |
| `E` | `log₂(1/|S|)` |
| `‖x‖` | distance from `x` to nearest integer |
| `λ(α)` | `inf_{q≥1} q ‖qα‖` (Markov/badly-approximable constant) |
| `M(α)` | `sup_n q_n ‖q_n α‖` (convergent-supremum, `≤ 1`) |
| `p_n/q_n` | reduced convergents of `α` |