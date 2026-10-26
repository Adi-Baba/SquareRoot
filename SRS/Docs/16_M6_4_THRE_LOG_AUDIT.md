# Paper 16 — M6.4: quantitative audit of the trilinear (three-log) bound

Version 1.1 (surgically revised).  We pick one effective theorem of
Matveev/Baker–Wüstholz type, write it with the `M6.1` data, and propagate
every parameter dependence.  Two consequences are now kept strictly in view:

- A three-log lower bound `|Λ| ≥ e^{−𝔇}` is a **lower bound on `E_L`**
  (in the small-`L` branch), not an upper bound; it therefore cannot
  certify `T(3,δ)`.  The only effective upper bound on `E_L` in this folder
  is the elementary `E_L ≤ 6N+O(1)` of Paper 9.
- The failure of the direct route is numerical, not logical: the heights
  `Π h(ε_i)` and the inverse-regulator coefficient `μ(E)` feed the Matveev
  constant, and no known bound for them is polynomial in the SRS input
  length `N`; the crudest generic figures give `𝔇 ≤ 2^{O(N)}`.

We keep parameters separate (D, h_i, regulator, B, N) and only then
combine.  All logs are base 2 unless marked ln.

---

## 16.1 Parameters inherited from M6.1

```
E = Q(√a, √b),    a = sqfree(d₁d₃),  b = sqfree(d₂d₃),   [E:Q] = 4 =: D.
ρ = d₃·u² ∈ O_E,   N_{E/Q}(ρ) = H(C)²,   H(C) = d₃²F(C) ∈ Z.   (M6.1 (8))
u = C₃ + C₁√(d₁/d₃) + C₂√(d₂/d₃),   S = √d₃·u.
σ₀ := the principal-roots real embedding, and the branch H(C) = ±1.
a := (a₁,a₂,a₃) ∈ Z³ in ρ = ±ε₁^{a₁}ε₂^{a₂}ε₃^{a₃}.    (M6.1 (16))
λ := ( log|σ_t(ρ)| )_{t=1..4} ∈ R⁴ = Σ_i a_i·(log|σ_t(ε_i)|)_{t}.  (M6.1 (17))
Height bound of ρ:  h(ρ) ≤ log d_max + 2 log M + log 6.           (M6.1 (21))
Exponent bound: |a_i| ≤ μ(E)·(log(1/(d₃|L|²)) + log d_max + 3 log M + O(1)).  (M6.1 (22))
```

The branch of interest is

```
log₂(1/|L|)  is as large as possible.  (M6.1 (★))
```

---

## 16.2 The trilinear theorem we use and its parameters

We use a standard effective lower bound for a linear form in three logarithms
of algebraic numbers (Matveev, *An explicit lower bound for a homogeneous
form in three logarithms...*, 2000; the placeholder constant is denoted `C₃`;
the shape is exact, the numerical constant is not the object of the audit):

> There are absolute `C₃` such that for nonnegative integers `a_i` and
> nonzero algebraic `α_i` with `Λ := a_1 log α_1 + a_2 log α_2 + a_3 log α_3
> ≠ 0`, defined with real logarithms, and for `B ≥ max(3,|a_i|)`,
>
> ```
> |Λ| ≥ exp( − C₃ · D^{5} · (1+log B)·Λlow ),   where   Λlow := (∏_i h(α_i)).
> ```

The abbreviative reading of the constant is

```
𝔇  :=  C₃ · D^{5} · (1+log B) · (Π_i h(α_i)) · (poly-log of algebraic heights),
```

with `log B` measured with the natural embedding.  The point of M6.4 is to
compute what 𝔇 becomes as a function of the SRS input size `N`.

*Direction of inequality (fixed early, so it is not lost).*  In the small-`S`
branch, `Λ = log(d₃L²) = log d₃ − 2E_L` is negative, `|Λ| = 2E_L − O(1)`.  A
three-log lower bound `|Λ| ≥ e^{−𝔇}` therefore gives

```
E_L ≥ e^{−𝔇}/2 ,
```

a **lower bound on `E_L`**, which is the *opposite* of the upper bound
`E_L ≤ poly(N)` that SRS needs.  The only effective upper bound on `E_L`
available in this folder is the elementary one of Paper 9,
`E_L ≤ 6N+O(1)`; a lower-bound triangle from Matveev cannot tighten it.

---

## 16.3 What `𝔇` becomes at fixed `r=3` (`D=4`, `n=3`)

```
𝔇  =  C₃ · 4⁵ · (1+log B) · ( h(ε₁)·h(ε₂)·h(ε₃) ) · (polylog of heights).
```

The coefficient factor `B = max(3, |a_i|)`, by M6.1 (22),

```
log₂ B  ≤  log₂ μ(E) + log₂( log(1/(d₃|L|²)) + log d_max + 3 log M + O(1) ).
```

For the branch `H=±1` we may use the norm-ceiling on the *content-free*
size: `log(1/(d₃|L|²)) ≤ log₂(1/|S|)·2 + O(1)`... but the only *a priori*
upper bound on `E = log(1/|S|)` that we currently own is the norm bound
`E ≤ (9/2)·log M + (3/2)log d₃ + O(1)`.  Substituting it into (22) gives

```
log₂ B ≤ log₂ μ(E) + log₂( E·N^{O(1)} + N^{O(1)} ) ≤ log₂ μ(E) + O(log N).   (26)
```

So, **under the elementary norm bound**, `B = O(μ(E)·poly(N))`.

---

## 16.4 The two irreducible N-dependent factors (the loss)

Two factors survive in `𝔇` and both are invisible to `N`:

### (i) The fundamental-unit heights `h(ε_i)` (already logarithmic scale —
 in unit size)

For a real quadratic field `Q(√a)` the fundamental unit `ε₀` is determined by
a Pell solution, and the real regulator means
`h(ε₀) = ½·log₂ ε₀ ˜  ½·Reg(E)`.  Known *general-discriminant* bounds do not
express `ε₀` as a polynomial function of the input height; the crude
Minkowski/Ankeny-style bound is exponential in the discriminant, i.e.

```
h(ε_i)  ≤  O( √d_max  )   (crude worst case of the associated quadratic
                           regulator; its `bit`-length is O(N), but the height
                           *itself* is an absolute logarithmic quantity that the
                           input alone does not control polynomially).
```

What is *meant* by the claimed failure of polynomiality is therefore at the
level of the conversion from the height-product in the Matveev constant to
`N`: no known upper bound `log₂(Π h(ε_i)) ≤ poly(N)` is available without an
extra regulation-style estimate.  We **record this coarse scale
(`Π h(ε_i) ≤ d_max^{3/2}`) as the worst-case bound that simpler estimates
produce, not as a proved property of the fields that arise from the `H(C)=±1`
branch.**

**Terminological precision.**  `h(ε_i)` is already a *logarithmic* invariant
(the unit size is `~ 2^{h(ε_i)}`).  Hence the "driver" of the blow-up is the
*logarithm-of-the-height-product* `log₂(Π h(ε_i))`, and this is the thing
that an effective theorem would have to bound by `poly(N)`.  The data available
in the literature bound it only by the associated block quadratic
field discriminants, and that is too coarse for the SRS input length.

### (ii) The inverse-regulator constant `μ(E)`

From M6.1 (19′), the exponent vector `a` is a fixed adjugate of the four
log-embeddings `log|σ(ε_i)|`.  Dividing by the regulator determinant gives

```
μ(E)  ≍  1/|det(M)|     (M := the 3×4 log-minor matrix),
```

which is fixed once the basis `{ε_i}` is fixed.  In `N` this can be as large
as `Reg(E)^{-1}`, a quantity that, in the worst-case regime, is again bounded
only by the surviving quadratic regulators.  It would be `poly(N)` exactly
where (i) is also `poly(N)`; the two are WFM-related via
`log₂ Π h(ε_i) ≈ log₁₀(1/R(E))^{-1}` up to the factor of degree.  So the
genuinely unmet ingredient is one scaler that makes the product
`Π h(ε_i)`-times-`μ(E)` polynomial in `N`.

### Direction of inequality (a checkpoint that later text must not lose)

For the dangerous small-`L` branch, `Λ = log(d₃L²) = log(d₃) − 2E_L` is
*negative*, and `|Λ| = 2E_L − O(1)`.  A three-log lower bound
`|Λ| ≥ e^{−𝔇}` therefore yields

```
2 E_L − O(1) ≥ e^{−𝔇}, i.e.  E_L ≥ e^{−𝔇}/2,
```

a **lower bound on `E_L`**, NOT an upper bound on it.  Plane arithmetic:
the norm-conclusion `|ρ_{σ₀}|·Π_{σ≠σ₀}|ρ_σ| = 1` plus the ceiling
`|ρ_σ| ≤ 9 d₃ M³` (`σ ≠ σ₀`) already gives the *effective upper* bound
`E_L ≤ 6N+O(1)`.  A lower-bound estimate in the direction Matveev provides
cannot improve that upper bound.  Hence the evaluation of 𝔇 below is a
*quantitative record of non-improvement*, not a candidate proof of
`T(3,δ)`.  The only honest conclusion is: the three-log route does not
recover `T(3,δ)`; the norm bound remains the effective rank-3 bound.

A technical hypothesis correction: the displayed Matveev form has
*non-negative* coefficients, but (M6.1 (16)) has arbitrary `a_i ∈ Z`.  To
use the displayed statement literally, absorb the minus one by replacing
`ε_i` by `ε_i^{−1}` whenever the corresponding `a_i < 0`, or quote a
standard integer-coefficient version as in Matveev 2000 directly.  From
here on `h(α_i)` denotes the height taken with this sign convention
(`h(ε_i^{−1}) = h(ε_i)`).

---

## 16.5 The M6.4 figure of loss, and why it is not `T(3,δ)`

We now propagate into (26) with `D=4`:

```
𝔇  =  C₃·4⁵·(1+log B)·(Π h(ε_i))·polylog(heights)
   ≤  C'·(1+log μ(E)+log N)·d_max^{O(1)}·polylog(N).
```

Since `d_max ≤ 2^N`,

```
𝔇  ≤  2^{O(N)}.
```

Baker would give the one-coordinate bound `|Λ| ≥ e^{−𝔇}` around (26) with
`D=4`.  Reading it correctly, this yields only `E_L ≥ e^{−𝔇}/2` — a lower
bound on `E_L` — and *not* an estimate of the requested kind
`E_L ≤ poly(N)`.  The provably valid upper bound on `E_L` remains the
elementary one carried through the document, `E_L ≤ 6N+O(1)`.

With that correction recorded, the audit provides the honest measurement of
*where* the trilinear formula would *deviate from being polynomial* if it
did succeed: it is the heights `Π h(ε_i)` and the inverse regulator
coefficient `μ(E)` in `𝔇`, both of which are bounded by the coarsest
estimates only with an `d_max^{O(1)}` factor.  Since a necessary step
would be `𝔇 ≤ poly(N)` and the current available bound only gives
`𝔇 ≤ 2^{O(N)}`, we have a recorded *exponential loss* in the three-log
route — even before converting the sign to the upper bound on `E_L`, which
would require the extra simultaneous input not present in a single
coordinate.  For `d_max ≥ 100`, the direct norm bound `E ≤ 6N+O(1)` is
therefore preferred; the three-log track is strictly worse and cannot
recover `T(3,δ)`.

**Consequence of the audit (not a proof of T(3,δ)):**  The three-log form
of the unit invariant `ρ = d₃u²` does not yield a polynomial-in-`N` bound
through the available effective-exponential estimates; it yields
`E_L ≥ e^{−𝔇}/2`, which is the opposite inequality to `T(3,δ)`, and where
the crude bound for `𝔇` is exponential (`𝔇 ≤ 2^{O(N)}`).  Consequently
`T(3,δ)` is *not established* by that method.  The sole effective rank-3
bound remains the elementary `E ≤ 6N+O(1)` provided by the full norm
(Paper 9).

---

## 16.6 The loss table requested by the user

| Source of loss            | Quantity(ies)                        | Desired scale         | M6.4 bound (after audit) |
|---------------------------|--------------------------------------|------------------------|--------------------------|
| field degree              | `D = [E:Q]`                          | `O(1)` (fixed `r==3`) | `4` (constant) ✓          |
| inverse-reg / exponents   | `μ(E) ` (M6.1 22), from `1/R_E`      | `poly(N)`             | not reducible to `poly(N)` from what is known (tied to the same regulator scale as (i)) |
| unit heights              | `Π h(ε_i)`                           | `poly(N)`             | not reducible to `poly(N)` from what is known; crudest generic bound leaves `≤ d_max^{O(1)}` |
| trilinear constant        | `C₃·D^5·(1+log B)`                   | `poly(N)`             | `poly(N)` once `B=poly` (holds by 26) |
| coefficient-size `B`      | `log B`                              | `poly(N)`             | `≤ log₂μ(E) + O(log N)`     |
| density of the ε_i        | a subtle half-rank issue?           | —                     | full rank, regulator≠0 (fixed)      |
| **resulting E**            | the desired `log₂(1/|S|)` upper bound | `N^{O(1)}`            | NOT obtained: the three-log route yields `E_L ≥ e^{−𝔇}/2` (lower) and 𝔇 ≤ 2^{O(N)}; the only known *upper* bound remains Paper 9: `E ≤ 6N+O(1)` |

---

## 16.7 M6.5: the exact source of the blow-up

The single source of the loss over the required scale is:

> **The fundamental-unit heights of a real biquadratic field are not bounded
> by the input data in any way that the currently available estimates
> express.**  Coarse discriminant-based estimates leave
> `h(ε_i) ≤ O(√d_max)` as the largest general figure, and since the
> three-log constant multiplies by `Π h(ε_i)`, the resulting `𝔇` is bounded
> only by `2^{O(N)}`.  This is not a proved property of the specific
> `H(C)=±1` fields; it is the *limit of what standard regulator bounds
> certify*.  The proved necessary reduction is: we would need an upper bound
> `log₂(Π h(ε_i)) ≤ O(log N)` (i.e. the fundamental units, in bit-size, are
> at most polynomial in the input height) in order to get `𝔇 ≤ poly(N)` from
> this route.

in the EU is actually what we need: the *heights* `h(ε_i)` themselves are at
most polynomial in `N`:

```
max_i h(ε_i)  ≤  N^{O(1)}   (♦)
```

(where (♦) is the precise regulatory inequality — the fundamental units of
the fields produced by the `H=±1` branch must have *logarithmic* height
polynomial in `N`, not merely bounded by the discriminant).  Under (♦) one
gets `Π h(ε_i) ≤ N^{O(1)}`, `(1+log B) ≤ O(log N)` (26), so
`𝔇 ≤ poly(N)`, and — after crossing to the correct simultaneous
multiembedding formulation — `E_L ≤ poly(N)`.  The point of M6.5 is that no
general bound gives (♦); current partial results for fundamental units only
give it at the *logarithm* of discriminant through `log_h`-level estimates.

---

## 16.8 Why this is OK for the overall program

The useful M6.4 output is: **for `r = 3`, a one-coordinate three-log estimate
cannot upper-bound `E`.  The only effective polynomial-in-`N` upper bound on
`E` that we possess is the direct algebraic one (Paper 9):
`E ≤ (9/2)·log M + (3/2)·log d₃ + O(1) ≤ 6N+O(1)`.  The three-log route
instead supplies the lower-bound `E_L ≥ e^{−𝔇}/2` with `𝔇 ≤
2^{O(N)}` under the crudest available regulator bounds; starting from
the weaker, standard estimates, it is strictly dominated by the norm bound.
So it cannot recover `T(3,δ)`.  The original open problem remains a
*uniform-in-`r`* statement `c · ‖C‖^{−(r-1+δ)}` with `log(1/c) ≤ poly(N)`,
and Paper 9/10/13/14 collectively show it is not produced by any full norm,
bimodular norm tower, or sign-pattern-linear-algebra manipulation.  The
precise sufficient condition for the *Baker-form* of the argument to close it
is the regulator inequality (♦) of §16.7.

See `TODO.md`, item M6.4/M6.5 for the small table in-place, and item M5
for the uniformity target.