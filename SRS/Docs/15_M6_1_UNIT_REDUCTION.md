# Paper 15 — M6.1: Exact unit reduction at `r = 3`

Version 1.0.  This is the first M6 deliverable, and it corrects a frequent
slip.  The branches on `F(C)=±1` do NOT, by themselves, put `u` in the unit
group of the order.  The correct integral unit is `ρ := d₃·u²`, and the
small-`L` branch is exactly the condition that this `ρ` be a unit.

---

## 15.1 The two fields, kept separate

Throughout, `d_1,d_2,d_3` are positive distinct squarefree integers.

```
K  :=  Q(√d₁, √d₂, √d₃),       full field, degree ≤ 8.
E  :=  Q(√(d₁d₃), √(d₂d₃)),   biquadratic, degree 4 (generically).
```

We have `u := C₃ + C₁√(d₁/d₃) + C₂√(d₂/d₃)` and `S = √d₃·u`.  Observe

```
u = C₃ + C₁·(√d₁/√d₃) + C₂·(√d₂/√d₃)
  = C₃ + C₁·(√(d₁d₃)/d₃) + C₂·(√(d₂d₃)/d₃),
```

so `u ∈ E`, and

```
V := d₃·u ∈ O_E.                                             (1)
```

`E` is a real biquadratic field `Q(√α,√β)` with `α = sqfree(d₁d₃)`,
`β = sqfree(d₂d₃)`.  Its Galois group is `(Z/2Z)²` (real, four embeddings).
**Unit rank.**  All four embeddings are real, so by Dirichlet,

```
O_E^× = {±1} × ⟨ε₁,ε₂,ε₃⟩,     rank = r₁ + r₂ − 1 = 4 + 0 − 1 = 3.                        (2)
```

So there are **three** independent units, not two.  This number 3 is what
forces the linear form in three logarithms below.

---

## 15.2 The four conjugates and the norm (Paper 9 restated)

The four real embeddings of `E` send

```
g₂ := √(d₁/d₃) ↦ ε·g₂,  g₃ := √(d₂/d₃) ↦ η·g₃,   ε,η ∈ {±1},
```

because they fix `g₂² = d₁/d₃` and `g₃² = d₂/d₃`.  Call the resulting four
images of `u`

```
v_{ε,η} := C₃ + ε C₁ g₂ + η C₂ g₃.
```

Then

```
F(C) := ∏_{ε,η} v_{ε,η} = N_{E/Q}(u).                                          (3)
```

We derived an exact expansion and an integral quartic form:

```
F(C) = C₃⁴ − 2C₃²(x²+y²) + (x²−y²)²,  x² = C₁²d₁/d₃,  y²=C₂²d₂/d₃,
d₃²F(C) = H(C)  ∈ Z.                                                        (4)
```

For `(V)` from (1), the norm is

```
N_{E/Q}(V) = d₃⁴F(C) = d₃²·H(C).                                             (5)
```

---

## 15.3 F(C) = ±1 does NOT make `u` a unit – the correct unit is `ρ`

The condition `F(C)=±1` is a rational-number condition on `N(u)`, with `u`
generically having a denominator.  The genuine integral object is `V` (1),
and (5) shows already that in the dangerous branch `H = ±1` one gets only

```
N_{E/Q}(V) = ±d₃²,          NOT ±1.                                          (6)
```

So `V` has norm `±d₃²`; in particular `V` is **not** a unit of `O_E`
(`E` is not the unit ring when `|N(V)| = d₃² > 1`).  The "unit" appearing
after the correct clearing is the quadratic element

```
┌────────────────────────────────────────────────┐
│   ρ := d₃·u²  =  V² / d₃  ∈ E                 │                              (7)
└────────────────────────────────────────────────┘
```

and its norm and integrality follow from (5):

```
N(ρ) = d₃⁴·F(C)² = d₃⁴·(H²/d₃⁴) = H(C)²,                                      (8)
```

so `ρ` is a unit of `O_E` **exactly in the branch `H(C) = ±1`**, i.e.

```
F-branch    H(C) = d₃²F(C) = ± 1   ⇔   ρ ∈ O_E^×,   N(ρ)=1.                   (9)
```

Moreover (7) expands with no irrationals remaining, because

```
x² = C₁² d₁/d₃ ⇒ d₃·u² = d₃C₃² + 2C₃(C₁ dg₂ + C₂ dg₃) + (C₁²g₂²+C₂²g₃²)d₃ + 2C₁C₂ d₃ g₂g₃,
```

and using `d₃g₂ = √(d₁d₃)`, `d₃g₃ = √(d₂d₃)`, `d₃g₂g₃ = √(d₁d₂)` (the last
is an identity in `E`: `g₂g₃ = √(d₁d₂)/d₃`), we obtain the **explicit
biquadratic quadratic form**

```
ρ = (d₁C₁² + d₂C₂² + d₃C₃²)
    + 2 C₁C₃·√(d₁d₃)
    + 2 C₂C₃·√(d₂d₃)
    + 2 C₁C₂·√(d₁d₂).                                                          (10)
```

Every coefficient is an integer.  In the branch (9), `ρ` is automatically
integral (it is a sum of integer-multiples of `√(integer)`), and its norm is
`+1`: it is a unit with `N(ρ)=1` (a *norm-one unit*).

**This is the concrete unit attached to the small-`L` branch.**

---

## 15.4 Where the small `L` lives inside `ρ`

Choose the distinguished embedding `σ₀` (the all-principal-roots one) for
which we assume `L = u_{++}` is extremely small.  Then from (7),

```
ρ_{σ₀} = d₃ · u_{++}² =  d₃·L² .                                              (11)
```

Let `|S|` be the actual SRS sum; `S_{σ₀} = √d₃·L`, so

```
|S_{σ₀}|² = d₃·|L|² = |ρ_{σ₀}|.                                             (12)
```

For the *other three* embeddings of `E`, the crude bound from Paper 9
applies via (3): on the support where each `v_{εη}` is one of the four
conjugates, `|v_{εη}| ≤ 3M^{3/2}`, `M = max(|C_i|, max d_i)`, so

```
|ρ_σ| = d₃·|v_{εη}|² ≤ 9 · d₃·M³   for σ ≠ σ₀.                              (13)
```

Taking the product of the four real absolute values and using (8)/(9),

```
|ρ_{σ₀}| · ∏_{σ≠σ₀} |ρ_σ| = |N(ρ)| = 1.                                       (14)
```

Hence the *small-`L` condition is quantitatively:*

```
log(1/|S|) = ½·log(1/|ρ_{σ₀}|) = −½·log(d₃|L|²)
           ≤  ½·( log(d₃) + 3log(9 d₃ M³) )· ...
```

More precisely, (14) plus (13) gives the trivial-of-P9 bound

```
2 log(1/|S|) = −log|ρ_{σ₀}| = Σ_{σ≠σ₀} log|ρ_σ|
             ≤ 3·(log 9 + log d₃ + 3 log M).                                (15)
```

This is the point at which the norm argument says: the lack of information
beyond the norm forces the `9/2` exponent.  To improve it, we must say more
than `|ρ_σ| ≤ 9 d₃ M³`: we must exploit that `ρ` is not an arbitrary
bounded element but a specific **unit**, i.e. a value of a multiplicative
group.

---

## 15.5 The exact logarithmic reduction

By (2), there are fundamental units `ε₁,ε₂,ε₃` and `a₁,a₂,a₃ ∈ Z` such that,
up to a root of unity,

```
ρ = ± ε₁^{a₁} ε₂^{a₂} ε₃^{a₃}.                                              (16)
```

Define the standard logarithmic embedding of `K^×` into `R^{4}` (all real),
using absolute values:

```
λ(β)  :=  ( log|σ_1(β)|, log|σ_2(β)|, log|σ_3(β)|, log|σ_4(β)| )ᵀ ∈ R⁴,
```

numbered so that `σ₀` is one of `σ₁,…,σ₄`, say `σ₁ = σ₀`.  Then from (16),

```
λ(ρ) = Σ_{i=1}^{3} a_i · λ(ε_i)  ∈ R⁴,                                       (17)
```

and from (14), the coordinates satisfy

```
Σ_{t=1}^{4} λ(ρ)_t = log|N(ρ)| = 0.                                         (18)
```

Likewise each `λ(ε_i)` satisfies the same zero-sum relation (norm-one
units: `Σ_t λ(ε_i)_t = log 1 = 0`), so (17) lies in the three-dimensional
zero-sum hyperplane

```
H := { x ∈ R⁴ : x₁+x₂+x₃+x₄ = 0 }.
```

The regulator argument: the four embeddings give a `3×4` real matrix
`M_{i,t} := log|σ_t(ε_i)|`, of rank exactly `3` (a classical fact; the
`3×3` minor over the last three coordinates is the real regulator
`R_E` of `E`, which is non-zero).  Hence the three rows `λ(ε_i)` span `H`,
and (17) gives an **explicit rational inversion formula** for the unit
exponents:

```
a_i  =  bnd of (Σ_t λ(ρ)_t·θ_{i,t}) / R_E ,                                  (19)
```

where `θ` is the adjugate/(3×3)-minor inverse; equivalently, in the usual
basis,

```
R_E·a_i = det(  [λ(ρ), λ(ε_1),…, λ(ε_3)] with row i replaced  ).            (19′)
```

We do *not* need a closed formula here — only the structural statement
that the integer vector `a := (a₁,a₂,a₃)` is determined by `λ(ρ)`, which has
first coordinate

```
λ(ρ)_{σ₀} = log|ρ_{σ₀}| = log(d₃|L|²),                                       (20)
```

and the remaining three coordinates `λ(ρ)_σ` for `σ ≠ σ₀`, each of which is
bounded above by (13).

---

## 15.6 Bound for the height of `ρ` and of `a`

From (10), `ρ` is an explicit integral quadratic form in `C₁,C₂,C₃` with
coefficients bounded by `max(d₁,d₂,d₃) := d_max`.  The raw multiplicative
height satisfies

```
h(ρ) ≤ log( d₁C₁² + d₂C₂² + d₃C₃² + 2|C₁C₃|√(d₁d₃)
          + 2|C₂C₃|√(d₂d₃) + 2|C₁C₂|√(d₁d₂) )
     ≤ log d_max + 2 log M + log 6 ,   M = max(|C_i|, d_max).                (21)
```

From (19), the integer vector `a` is controlled by `λ(ρ)` and the regulator
of `E`.  Because the four `λ(ρ)_σ` sum to `0` and the `σ₀` one is the most
negative, the three others are positive-summing:

```
Σ_{σ≠σ₀} λ(ρ)_σ = − λ(ρ)_{σ₀} = log(1/(d₃|L|²)).
```

Each `λ(ρ)_σ (σ≠σ₀)` is in turn bounded above by (13) in absolute value.
So every entry of `λ(ρ)` is at most `log(1/(d₃|L|²)) + log d_max + 3log M +
const`, and so

```
|a_i|  ≤  μ(E)·( log(1/(d₃|L|²)) + log d_max + 3 log M + O(1) ),            (22)
```

where `μ(E)` is an algebraic constant of `E` (a bound for the inverse of the
logarithmic-embedding residue; it depends on the fundamental units but is
fixed once `d₁,d₂,d₃` are fixed).

This gives the **height/exponent control** the uniform conjecture would
need at `r=3`; the dependence on `log(1/|L|)` is logarithmic, as it should
be.

---

## 15.7 The exact reduction statement

Putting (16), (17), (20) together, in the branch `H(C) = ±1` the smallness
of `L` is *by definition* a smallness of one coordinate of a rank-3 lattice
in the unit-log hyperplane.  Specifically, with `Λ := λ(ρ) ∈ H`, we have
produced an explicit integer vector `a=(a₁,a₂,a₃)`, determined uniquely by
(19), satisfying

```
Σ_{i=1}^{3} a_i · λ(ε_i)_{σ₀} = log(d₃|L|²).                              (★)
```

The other three coordinates of (17) are the same linear form applied to the
other three embeddings, and (14)–(15) show they cannot cancel the `σ₀`-side
excess by more than `3 log(9 d₃ M³)` without violating norm-one.  Hence
(★) is the only equation that involves the genuinely small quantity `|L|`;
the remaining three coordinates are constrained to the single real interval
`0 ≤ λ(ρ)_σ ≤ log(9 d₃ M³)` for `σ ≠ σ₀`.

**Conclusion.**  The SRS separation `E ≤ poly(N)` at fixed `r=3` is now
equivalent to the *effective inequality*

> `` Given a fixed real biquadratic field `E`, three fundamental units
> `ε₁,ε₂,ε₃`, and an integer `a ∈ Z³`, if the first coordinate of
> `Σ a_i λ(ε_i)` is at most `log(d₃|L|²)` while the other three coordinates
> of the same sum are all at most `log(9 d₃ M³)`, then `log(1/|L|)` is
> bounded above by a polynomial in the input size. ''

In **M6.1** we have reduced rank-3 SRS to exactly that statement.  It is a
single effective lower-bound problem on a **three-term linear form in
logarithms** of the fundamental units of a real biquadratic field, with the
coefficient-vector `a` controlled by (22).  Applying any effective estimate
to *this* quantity is M6.3.

---

## 15.8 Reference: what is *not* claimed

- We do **not** claim that `F(C)=±1` puts `u` in `O_E^×`.  It does not (6).
  The norm-`1` unit is `ρ = d₃u²` (7)–(8).
- We do **not** claim `H(C) = ±1` automatically solves the problem.  This
  section merely exhibits an exact reduction; M6.3 evaluates the exponent.
- The branch `|H(C)| ≥ 2` is already effective by Paper 9 together with the
  norm argument: `|L| ≥ 2/(d₃²·(3M^{3/2})³)`.  Only `H(C)=±1` is contentful.