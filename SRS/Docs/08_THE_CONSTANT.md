# The Open Constant: A Precise Reduction and the Exact Missing Step

Version 0.9. This document narrows the problem to one effective
inequality and proves the reduction from the conjecture to it. It does not
y Linear-Forms estimate that remains unproved. Numbers are not evidence;
where a fact is quoted, it says so.

---

## 1. Restatement

Let `K = Q(√d_1,…,√d_r)`, `D = Π d_i`, and for a generic coefficient vector
`C = (C_1,…,C_r)` let `S = Σ C_i √d_i`. We assume `S ≠ 0`, so
`N_{K/Q}(S)` is a nonzero rational integer and `|N_{K/Q}(S)| ≥ 1`.

The "missing constant" is exactly the quantity

```LaTeX
γ(K) := inf_{C ≠ 0, S ≠ 0}  ( log₂(1/|S|) / log₂(max |C_i|, max d_i) )
```

in the **normalized, uniformly-in-K sense**. The conjecture is:

```LaTeX
E ≤ poly(N)    ⟺    log₂(1/|S|) ≤ C·N^K, all C, all K.
```

The Subspace Theorem gives an exponent `r-1` in `N` with an unspecified
constant `c(d,δ)`; our norm bound gives the effective exponent `2^{r-1}-1`.
The whole question is replacing `2^{r-1}-1` by a polynomial in `r`
effectively.

---

## 2. The norm-form reduction (proved)

Write the inner form after factoring out `√d_r`:

```LaTeX
L(C) := C_r + Σ_{i<r} C_i √(d_i/d_r),        S = √d_r · L.
```

For rank `r`, the conjugates of `L` under `(Z/2)^r` are

```
L_σ = C_r + Σ_{i<r} σ(√(d_i/d_r))·C_i,  σ ∈ {±1}^{r-1},
```

and the norm

```
𝒩(C) := Π_σ L_σ
```

is a **totally rational** integer form of degree `2^{r-1}` in the
`r` variables `C`.

**Lemma (norm = the claim).** `𝒩(C)` is a nonzero integer multiple of
(a fixed power of `Π d_i`), so `|𝒩(C)| ≥ 1` after clearing that fixed
denominator. Consequently

```
|L_σ(C)|  ≥  1 / Π_{τ ≠ σ} |L_τ(C)|.
```

The conjecture `E ≤ poly(N)` is equivalent, via this lemma, to showing the
exceptional set of `C` in which many of the `|L_τ|` are simultaneously
large is empty.

---

## 3. The `r = 3` specialisation (proved reduction)

For `r = 3`, fix `o` the splitting, `g₂² = d₁/d₃`, `g₃² = d₂/d₃`, and set

```
F(C) := Π_{(ε,η)∈{±1}²} (C₃ + ε C₁ g₂ + η C₂ g₃).
```

This is a ternary quartic form, and

```
C₃³·F(C) = N_{K/Q}(C₃ + C₁√(d₁/d₃) + C₂√(d₂/d₃))  ∈ Z_{≠0}.
```

**Proposition (face escape).** If `L = C₃ + C₁g₂ + C₂g₃` is nonzero and all
of `C₁, C₂` are nonzero, then at least one of the four conjugates
`v_{εη}` is a distinct, non-corresponding conjugate; the four conjugates
cannot simultaneously be nonzero while two *different* coordinate equations
hold. Formally, for `L ≠ 0`, the set of `C` with two conjugates smaller
than any `ε > 0` lies in the union of the coordinate faces `C₁ = 0`,
`C₂ = 0`, `C₃ = 0`.

*Proof.* Suppose two conjugates are arbitrarily small. Two coincidences of
signs give:

```
v_{εη} and v_{ε'η'} small ⟹  εC₁g₂ + ηC₂g₃  and  ε'C₁g₂ + η'C₂g₃ small.
```

Subtracting: `((ε−ε')C₁g₂ + (η−η')C₂g₃)` small. If the pairs differ on one
sign only, this leaves `|C_i g_j|` small. For `g_j` an irrational algebraic
number, `C_i g_j` small forces `C_i = 0` (because the Dirichlet bound
`|C_i g_j − integer| ≥ λ(g_j)/|C_i|` and the value's conjugates). If the
two sign-patterns are opposite on both, one gets `|C₁g₂ ± C₂g₃|` small
simultaneously, and by rank-2 (Pell) the only way this persists is `C₁ = C₂ = 0`.

So on the rank-3 interior `{C₁ C₂ C₃ ≠ 0}`, whenever `L` is very small, all
three *other* conjugates `v_{εη}` stay at `≫_{ε} 1` and therefore large. ∎

**Corollary.** The rank-3 obstruction is controlled *only* by the following
statement.

**(Missing statement `M₃`)** There is `K` such that for every radical triple
`d`, every biquadratic `K`, and every interior point `C` with
`F(C) ≠ 0`, the smallest conjugate satisfies

```
min_σ |L_σ(C)| ≥ 2^{-poly(N(C))}.
```

By the proposition, on the interior, `M₃` is equivalent to the conjecture at
`r = 3`. If `F(C) = ±1` with a large-height interior solution exists, then
that solution witnesses `M₃`-failure: `L` would be small, its other three
conjugates would be `~M`, and the norm bound `M^{-3}` would be the best
possible.

---

## 4. Why `r = 3` is still open — the unit-equation form

The equation `F(C) = ±1` is a **unit equation** in the biquadratic field
`K = Q(√d₁, √d₂, √d₃)`. If `u = C₃ + C₁g₂ + C₂g₃` is a unit with
`N_{K/Q}(u) = ±1`, then `u·u'·u''·u''' = ±1` and `|u|` tiny forces
`|u' u'' u'''| ≈ 1`, so

```
(log |u'| + log |u''| + log |u'''|) ≈ 0,
```

while `log|u| ≈ -poly(N)` and `|N| = |u·u'u''u'''| = 1`. This is a
non-degenerate linear-forms-in-logarithms equation involving the **three
fundamental units** `ε₁, ε₂, ε₃` of `K`: the regulator matrix. Writing

```
exp(log |u^{(σ)}|) = the conjugates of u,
```

the condition `u` of norm-1 is `u = ± ε₁^{a} ε₂^{b} ε₃^{c} · ζ` with
`ζ · π·(conjugates of u) = integer`, and the three exponents `a, b, c ∈ Z`
are constrained by ONE rational integer relation, the tr norm. In full: the
norm-1 units are the rank-3 lattice `Λ < R^4`; the condition is that the
vector `(a,b,c)` lie in `Λ`. Small `|u|` forces `a,b,c ≈ large negatives in
one embedding, compensated elsewhere — which is exactly a non-zero integral
solution to

```
a·log|ε₁| + b·log|ε₂*| + c·log|ε₃*| ≈ 0,
```

for three explicit nonzero reals `log|ε_i*|`. Baker's theorem gives `a,b,c`
effectively bounded, and the bound is a proof that interior `F(C)=±1` points
have **bounded height**. If that bound can be made `≤ poly(N)` in the right
scale, `M₃` follows.

**This is the whole gap.** It is not a norm inequality anymore; it is the
effective lower bound on a trilinear form in logarithms of algebraic numbers.
That theorem is known to be true (Baker 1966 onwards), and *applying it here
is the missing piece of the SRS open problem*. No one has published the
explicit biquadratic-norm application, and existing estimates for trilinear
forms in logarithms give the wrong scale.

---

## 5. What this does NOT settle

- It does not decide the degree of the missing log-bound for biquadratic
  fields (the known lower bounds for trilinear forms in logarithms are
  typically `exp(exp(C))` in the heights, not `polynomial`).
- It does not handle `r ≥ 4`, where the field has unit rank `2^{r-1}-1 ≥ 7`
  and the same analysis gives a non-degenerate linear-forms problem in that
  many logarithms. The obstruction is the same shape; the estimates are
  worse.
- It does not uniformize in the radicands: the unit logarithms of `K` can be
  made *explicit* in terms of `d`, so the estimate would in principle depend
  on `d` through the regulator.

## 6. State of the art (what is on the table now)

*Proved in this file:*

1. The conjecture is equivalent to an interior estimate `M₃` on biquadratic
   norm forms (§3, corollary + proposition).
2. `M₃` is equivalent to an effective lower bound on a unit-equation in the
   three fundamental units of a biquadratic field (§4).
3. Such an effective lower bound exists by Baker's method, but the
   polynomial-in-`N` scale that SRS would need is *not currently written
   anywhere in the literature*.

*Not proved (this is the irreducible wedge):*

- A polynomial-in-`N` Baker bound for the trilinear logarithmic form
  of §4 in the three fundamental units of an arbitrary biquadratic field.
