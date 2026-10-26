# SRS — Zero Test Correctness

Version 1.0. This document proves the correctness of the zero test: given
integer pairs `(c_j, a_j)`, the algorithm `decompose_and_aggregate` correctly
decides whether `Σ c_j sqrt(a_j) = 0`.

---

## 1.1 Statement

Let `k ≥ 1`, `c_j ∈ Z`, `a_j ∈ Z_{≥1}`. Define `S = Σ_{j=1..k} c_j sqrt(a_j)`.
The algorithm:

1. Decomposes each `a_j` as `a_j = t_j² · d_j` with `d_j` squarefree.
2. Aggregates coefficients by `d_j`: for each distinct squarefree `d`,
   `C_d = Σ_{j: d_j = d} c_j t_j`.
3. Returns `S = 0` iff `C_d = 0` for all `d`.

**Theorem 1.1 (Zero Test Correctness).** The algorithm returns `S = 0`
iff the true mathematical sum is zero.

---

## 1.2 Proof

The proof relies on the linear independence of square roots of distinct
squarefree integers over `Q`.

**Lemma 1.2 (Irrationality Step).** Let `t_1, …, t_r` be positive integers with
linearly independent parity vectors over `F₂`, and let `s ≥ 1` be an integer
such that `{t_1, …, t_r, s}` also has independent parity vectors. Then

```
sqrt(s) ∉ Q(sqrt(t_1), …, sqrt(t_r)).
```

*Proof.* By induction on `r`.

*Base `r = 0`.* The field is `Q`. `sqrt(s) ∉ Q` for non-square `s`. Since the
parity vector of `s` is non-zero (otherwise it would be a square and dependent
with the empty set), `s` is not a square. So `sqrt(s) ∉ Q`.

*Inductive step.* Assume the claim holds for `r − 1`. Suppose for contradiction
`sqrt(s) ∈ Q(sqrt(t_1), …, sqrt(t_r))`. By the induction hypothesis,
`sqrt(t_r) ∉ Q(sqrt(t_1), …, sqrt(t_{r-1}))`, so the degree of the extension is
`[Q(sqrt(t_1), …, sqrt(t_r)) : Q(sqrt(t_1), …, sqrt(t_{r-1}))] = 2`. Therefore
every element of the larger field can be written uniquely as
`a + b sqrt(t_r)` with `a, b ∈ Q(sqrt(t_1), …, sqrt(t_{r-1}))`. Write

```
sqrt(s) = a + b sqrt(t_r).
```

Squaring gives `s = a² + b² t_r + 2ab sqrt(t_r)`. Since `sqrt(t_r)` is
irrational over the smaller field, the coefficient of `sqrt(t_r)` must vanish:
`2ab = 0`. So either `a = 0` or `b = 0`.

- If `a = 0`, then `sqrt(s) = b sqrt(t_r)`, so `s/b² = t_r`. But `b` is in the
  smaller field, so `s/b²` has parity vector in the span of `t_1, …, t_{r-1}`.
  Then `t_r` would be dependent on `t_1, …, t_{r-1}` — contradiction.
- If `b = 0`, then `sqrt(s) = a ∈ Q(sqrt(t_1), …, sqrt(t_{r-1}))`, so by the
  induction hypothesis `s` is dependent on `t_1, …, t_{r-1}` — contradiction.

Thus `sqrt(s) ∉ Q(sqrt(t_1), …, sqrt(t_r))`. ∎

**Corollary 1.3 (Linear Independence).** If `d_1, …, d_m` are distinct
squarefree integers `≥ 2`, then `{sqrt(d_1), …, sqrt(d_m)}` is linearly
independent over `Q`.

*Proof.* The squarefree parts correspond to parity vectors of prime
exponents mod 2. Distinct squarefree integers have distinct parity vectors,
and a set of parity vectors is linearly independent over `F₂` iff no nonempty
subproduct of the corresponding squarefree integers is a perfect square.
By Lemma 1.2, the square roots are linearly independent over `Q`. ∎

**Proof of Theorem 1.1.** Decompose `a_j = t_j² · d_j` with `d_j` squarefree.
Then

```
S = Σ_j c_j t_j sqrt(d_j) = Σ_{d ∈ D} ( Σ_{j: d_j = d} c_j t_j ) sqrt(d),
```

where `D` is the set of distinct squarefree parts. By Corollary 1.3, the
`sqrt(d)` are linearly independent over `Q`. Since the coefficients `C_d` are
integers, the sum is zero iff every `C_d = 0`. The algorithm checks exactly
this condition. ∎

---

## 1.4 Remarks

- The algorithm runs in polynomial time: squarefree decomposition is
  polynomial-time (trial division up to `√a_j`, or using the exact routine in
  `squarefree.py`), and aggregation is linear in `k`.
- The `r = 0` case (all `a_j` perfect squares) is handled by the same logic:
  then `d_j = 1` for all `j`, and `S = Σ c_j t_j` is an integer; the algorithm
  correctly checks whether that integer is zero.

---

## 1.5 What Is Not Claimed

- The algorithm does not *certify* the sign of `S` when `S ≠ 0`. Sign
  certification is the subject of `SRS_02_Main_Theorems.md`.
- The algorithm does not give a separation bound. It only decides exact
  equality.