# Paper 9 — The Rank-3 Quartic Norm Form (exact) and the Effective Bound

Version 1.0.  We derive the exact quartic norm form for the inner form at
rank 3 and record the *unconditional* effective separation it already gives.
The conclusion is: at fixed `r = 3` the SRS conjecture `E ≤ poly(N)` is
**already a theorem**, with an explicit exponent.  The open question is only
uniformity in `r`.

---

## 9.1 Setup

Let `d_1, d_2, d_3` be distinct squarefree, and put

```
g₂ = √(d₁/d₃),   g₃ = √(d₂/d₃),   u = C₃ + C₁ g₂ + C₂ g₃.
```

Then `S = C₁√d₁ + C₂√d₂ + C₃√d₃ = √d₃ · u`.  The four conjugates of `u`
under `g₂ ↦ ±g₂`, `g₃ ↦ ±g₃` are

```
v_{εη} = C₃ + ε C₁ g₂ + η C₂ g₃,   ε, η ∈ {±1}.
```

Note `√d₃·v_{εη}` are four of the (up to) eight conjugates of `S`; their
product is `N_{Q(g₂,g₃)/Q}(u)`.

---

## 9.2 Theorem 9.1 (exact quartic form)

Set `x = C₁ g₂`, `y = C₂ g₃`.  Then

```
F(C) := v_{++}v_{+-}v_{-+}v_{--}
      = [(C₃+x)² − y²][(C₃−x)² − y²]
      = (C₃² + x² − y²)² − 4 C₃² x².
```

Expanding,

```
F(C) = C₃⁴ − 2 C₃²(x² + y²) + (x² − y²)².      (★)
```

Because `x² = C₁² d₁/d₃` and `y² = C₂² d₂/d₃`, the right side is free of
irrationals.  Clearing the common denominator by `d₃²`:

```
d₃²·F(C) = C₃⁴ d₃² − 2 C₃² d₃·(C₁²d₁ + C₂²d₂) + (C₁²d₁ − C₂²d₂)².
                                               (integer, all C_i, d_i ∈ Z)
```

So `H(C) := d₃²·F(C) ∈ Z`, an integral quartic form.

---

## 9.3 Theorem 9.2 (no interior zeros)

`F(C) = 0` iff some `v_{εη} = 0`, i.e. `C₃ = −ε C₁ g₂ − η C₂ g₃`.  That is a
`Q`-relation among `1, g₂, g₃`.  But `g₂, g₃` have radicands
`d₁/d₃, d₂/d₃` whose mod-square classes are two independent elements of
`Q^×/(Q^×)²` whenever `d_1,d_2,d_3` have full `F₂`-rank 3, so
`1, √(d₁/d₃), √(d₂/d₃)` are linearly independent over `Q` (Paper A, Thm A.3
applied to the independent subset of `{d₁d₃, d₂d₃, d₃²} = {d₁d₃, d₂d₃, 1}`).
Hence the relation is only possible if `C₁ = C₂ = C₃ = 0`.

**Therefore `F(C) ≠ 0` for every nonzero integer `C`, and in particular**
`H(C) = d₃²F(C)` is a **non-zero integer**.  Equivalently,

```
|F(C)| ≥ 1/d₃²       (all C ∈ Z³ \ {0}).
```

---

## 9.4 The effective rank-3 bound

Take the conjugate `v_{++} = u` to be the one whose reciprocal we bound.
Using Theorem 9.2,

```
|u| = |v_{++}| ≥ 1 / (d₃² · |v_{+-}| |v_{-+}| |v_{--}|).
```

Each `v_{εη}` is bounded above by the triangle inequality.  With
`M = max(|C_i|, max d_i)`,

```
|v_{εη}| ≤ |C₃| + |C₁|√(d₁/d₃) + |C₂|√(d₂/d₃)
         ≤ M + 2·M·√(M/2)
         ≤ 3·M^{3/2}       (for M ≥ 2).
```

So

```
|u| ≥ 1/(d₃² · (3 M^{3/2})³) = 1/(27·d₃²·M^{9/2}),
```

and, since `S = √d₃·u`,

```
|S| ≥ 1/(27·d₃^{3/2}·M^{9/2}).
```

**Theorem 9.3 (effective rank-3 separation).** For squarefree distinct
`d₁,d₂,d₃`, every integer vector `C ≠ 0` with `S ≠ 0` satisfies

```
log₂(1/|S|) ≤ (9/2)·log₂ M + (3/2)·log₂ d₃ + log₂ 27.
```

Because `M` and `d₃` are each `≤ 2^N`, this is

```
E ≤ (9/2)·N + (3/2)·N + 5 = 6·N + 5,
```

a purely effective, polynomial bound in `N`.  **Consequently the SRS
conjecture `E ≤ poly(N)` is TRUE and effective at `r = 3`.**  The same
argument with `2^{r−1}−1` complementary factors reproduces the general
fixed-rank effective bound of Paper B.

---

## 9.5 Where the sharpness is lost (a comment, not a theorem)

The bound `|v_{εη}| ≤ 3M^{3/2}` is generic.  The Subspace-Theorem predicts
the sharper exponent `2` (i.e. `|L| ≳ c/‖C‖²`).  The gap to `9/2` is exactly
`|v_{+-}v_{-+}v_{--}| ≲ M^{9/2}` versus an expected `|v_{+-}v_{-+}v_{--}| ≲ M²`
on the thin set of small `u`.  Proving that extra room is, once again, an
*effective Subspace-Theorem* statement at `r=3`; it is not available from the
norm argument, and not available from Baker theory at polynomial scale.
So the effective rank-3 theorem stands at exponent `9/2`, while the sharp
rate `2/3` is known to hold asymptotically (Paper C, Thm C.1) via the
ineffective constant.

---

## 9.6 Corollary for the program

At every fixed rank `r`, SRS is effective.  There is no uniform conjecture gap
at fixed `r`.  The open problem is strictly: the effective uniformity of the
Subspace-Theorem constant `c(d,δ)` across **growing** `r`.  Any proof must
use that `r` is not fixed — the norm form has degree `2^r` and its denominator
clearing gives `E ≤ (2^{r−1}−1)·(3/2)N`, which is polynomial in `N` iff `r`
is bounded.  Hence no amount of purely algebraic manipulation of the norm
form at fixed `r` can prove the uniform conjecture; one needs a new input
that beats the degree growth.  This is recorded as the residual open target.