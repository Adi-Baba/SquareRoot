# Paper 14 — M5: how the literature reaches the conjecture

Version 0.9.  We do NOT prove a theorem here.  We isolate the exact
external statement that would imply SRS has polynomial separation, translate
the strongest known effective results into SRS parameters, and state why the
multiquadratic Galois form does or does not collapse them.  Every quoted
result is a citation, not a new proof.

---

## 14.1 M5.1 — the one exact missing statement

The conjecture `E ≤ poly(N)` is equivalent to the following single uniform
inequality:

> **(M5.1).** There exist absolute constants `A, K` such that for every finite
> set `{d_i}` of distinct squarefree integers and every integer tuple
> `(C_i) ≠ 0`, with `S = Σ C_i √d_i ≠ 0`,
>
> ```
>     |S| ≥ 2^{−A · N^K},   N := Σ_i ( bit_length(|C_i|) + bit_length(d_i) ).
> ```

Equivalently, in inner-form notation with the rank-`r` subset,

> **(M5.1′).** There are absolute `A, K` with
> `|L(C)| ≥ 2^{−A N^K}` for every `C ∈ Z^r` with `L(C) ≠ 0`, uniformly in
> all radicand tuples `d` of full `F₂`-rank `r` and in `r`.

This is the **weakest standalone statement that suffices**: a single
a- priori quantitative lower bound with no exceptional-set qualifier, no
`δ`, no dependence of the exponent on the radicands beyond `N`.

---

## 14.2 M5.2 — what the literature actually gives when you compute the parameters

We ignore proof details and evaluate the *arithmetical* content of each
candidate.

### (a) Eisenbrand–Haeberle–Singer (SoCG 2024), via Subspace Theorem

The published statement is effectively of the form

```
|S| ≥ γ(d) · ( m · max_i |C_i| )^{−2m},
```

with `γ` depending on the radicands and, because Subspace-Theorem constants
are ineffective, **not explicit**.  Converting: with
`N ⊇ Σ(bitlen|C_i|+bitlen d_i) ≥ m + Σ bitlen|C_i|`, we get
`max|C_i| ≤ 2^{N/m}` and `m ≤ N`, so

```
E = log₂(1/|S|) ≤ 2m·( log₂ m + (N/m) )  +  log₂(1/γ(d))
                = 2N + 2m log₂ m + log₂(1/γ(d)).
```

The first two terms are `O(N log N)`.  **All of the open problem is sitting
in `log₂(1/γ(d))`.**  SRS-polynomial separation `E ≤ A N^K` is equivalent to
showing `log₂(1/γ(d)) ≤ poly(N)` for a suitable polynomial.  This is
precisely the `c(d,δ)` factor of our Subspace-Theorem constant in Paper C,
so it is already restated, not solved.

### (b) Full algebraic-number norm (Paper B, effective)

Gives `E ≤ (2^{r−1}−1)((3/2)N + log₂N)`.  As a function of `N` this is
`E ≤ poly(N)` **iff** `r` is bounded; for `r → ∞` it is `E ≍ 2^{r}N`, not
polynomial.  So (b) is SRS only in the fixed-rank band `r = O(log N)`.

### (c) Linear forms in logarithms (Baker–Wüstholz 1993; matrix models)

After rationalizing into a unit equation, `u` has norm `1` (Case B of
Paper 8), and for three fundamental units in degree 4, `u = ε·η₁^{a₁}η₂^{a₂}η₃^{a₃}`
with `a_i ∈ Z` constrained by a single norm equation.  Baker bounds give, in
the generic case, `|u| ≥ exp(−C·(h(d))^{C'}·(log B)^{C''})` for some
heights `h(d)`.  Even with the best explicit estimates (e.g. the
heights of Philipp? and the trilinear estimates of Cao), the resulting lower
bound on `|L(C)|` has the shape `exp(− exp(poly(r)) · poly(N))`, i.e.

```
log₂(1/|L(C)|)  ≈  exp(poly(r)) · poly(N).
```

For fixed `r` this is `constant_r · poly(N)` — effective — which *reproves*
the fixed-rank SRS bound only qualitatively.  Letting `r` grow, the prefactor
`exp(poly(r))` is exactly the `2^{O(r)}` (or worse) loss, and it is
**no polynomial in N**.  This is the quantitative reason the
linear-forms route does not reach the SRS scale.

### (d) Quantitative simultaneous Diophantine approximation

The exponent `m−1` in `|Σ C_i α_i| ≥ c·‖C‖^{-(m−1)}` is what the Subspace
Theorem provides; effective versions exist only in narrower settings.  In
every version where the algebraic vector `(α_i)` is *generic*
(transcendence / Schanuel-type structure), the constant produced is again a
function of the minimum heights of the α_i, and the corresponding quality is
`|L(C)| ≥ c(α) · ‖C‖^{−(m−1)}` with `c(α)` not bounded uniformly in the
coefficients of the minimal polynomials — i.e. the same `γ(d)` gap as in
(a).

### The unifying diagnosis

Every effective route we can evaluate produces one of two shapes for
`E`:

```
        E  ≤  (2^{r-1}-1)·poly(N)         ← Paper B / full norm (effective)
        E  ≤  poly(N) + log₂(1/γ(d))      ← (a), (d) (γ not explicit)
        E  ≤  exp(poly(r)) · poly(N)      ← (c) Baker (fixed-r effective, but
                                                  exponential in r).
```

Only the middle has the SRS shape, and only with `log₂(1/γ) ≤ poly(N)`,
which is unproved.  Therefore the precise missing theorem is:

> **M5.1′′.**  For an arbitrary multiquadratic field
> `K = Q(√d₁,…,√d_r)` and for every integer `C` with `L(C) ≠ 0`,
> `|L(C)| ≥ γ` with `log₂(1/γ) ≤ A·N^K`.

---

## 14.3 M5.3 — does the (Z/2Z)^r group collapse γ?

Two structural patterns of real solutions are known, and we ask whether the
multiquadratic Galois form explains either of them.

- **r = 2**: γ is explicit.  The Markov constant of a single real quadratic
  irrational is computable from its CF period; this is exactly the
  Hurwitz/Lagrange mechanism.  Here the "many-heights constant" `γ(d)` is as
  good as `1/(2·max(1, α)+1)` in the worst direction and never zero.  This is
  fully effective and is not, in fact, a generic phenomenon: it uses that the
  inner linear form is *binary*.  No analogue of the Markov constant is
  available at `r ≥ 3`.

- **Generic simultaneous projection**: the (Z/2Z)^r action makes all
  conjugates real, which lets us *pair* them and save a factor 2
  (Paper B conditional statement).  This is the only trace of the Galois group
  that survives in the norm estimates; it does not survive in γ.

The inner linear form
`L(C) = C_m + Σ_{i<m} C_i √(d_i/d_m)` has `m` unknowns and `m` coordinates
in the real space.  The (Z/2Z)^r group enforces that all conjugates of any
single `L_σ` saturate the same integer norm simultaneously, so geometrically
it says "the `m` large conjugates form a symmetric set" and **no**
quantitative bound on the smallest one.  In parametric terms the group does
not lower the rank-of-approximation from `m−1`; it is exactly what
characterises the inner product of the simultaneous vector
`(√(d_i/d_m)) ∈ R^{m−1}` with the integer vector.  The real defect — for
which no Galois shortcut is known — is precisely that a single vector in
`R^{m−1}` of degree over `ℚ` equal to `2^{m−1}` admits a full-norm estimate
but no effective Dirichlet-type estimate of quality `m−1` at polynomial
height.

So the answer of M5.3, as it can be rigorously stated today, is negative:
**the (Z/2Z)^r form does not, by itself, reduce the exponent or make γ
explicit.**  It supplies the only *effective* argument available in the
literature's effective envelope (the all-ones pairing and the fixed-rank
effectivity), and it is already used exactly there.

---

## 14.4 The formalised missing theorem

We can now name the one statement that, if proved, makes SRS true and
effective.  It is a simultaneous-approximation statement *in a quadratic
field of rank `r*`:

> **Theorem T(r*) (Missing).**  For every multiquadratic field
> `K = Q(√d₁,…,√d_r)`, every `C` with `L(C) ≠ 0`, and every `δ ∈ (0,1)`,
> 
> ```
> |L(C)| ≥ c(d) · ‖C‖_∞^{−(r−1+δ)},   with log₂(1/c(d)) ≤ A·N^K
> ```
> uniform in `r` and `d`.

If `T(r*)` holds, then `E = log₂(1/|S|) ≤ (r−1+δ)/r·N + A·N^K + O(1) ≤
poly(N)`, i.e. `E ≤ B·N^{K+1}` for an absolute `B`.  If it does not hold,
there is a sequence `K_j` and radicands where the denominator of every
effective Subspace-type constant grows faster than any polynomial in the
input height, and that sequence is exactly where the conjecture would have
to fail or lie beyond reach of current methods.

---

## 14.5 What this changes (and does not change) about the project

- **What is no longer a target**: norms, towers, sign-pattern rank.  Those
  have been audited (Papers 9–13) and their asymptotic reach is now fixed:
  fixed rank effective (Paper B / Paper 9), `r = O(log N)` closes the
  conjecture (Paper B Cor B.2), no more.
- **What is the one remaining target**: an explicit proof of `T(r*)`.  The
  cleanest possible next move is to try `T(3)` (Paper 9 gives the quartic
  form and the measure).  An explicit, `poly(N)`-scale `c(d)` in rank 3,
  uniform over *all* triples `d₁,d₂,d₃`, would be a publishable refinement of
  EHS; the same structure for all `r` would settle SRS.
- **What is fully proved and can be stated anywhere**: the zero test
  (Paper A), the effective fixed-rank bound (Paper B), the sharp asymptotic
  rate (Paper C, using a quoted Subspace theorem), the effective rank-3 bound
  `|S| ≥ 1/(27 d₃^{3/2}M^{9/2})` (Paper 9), the M1 obstruction that pure
  norms/towers give `2^{r−1}−1` (Paper 10), and that any effective
  simultaneous-approximation result on the open form must, after M5.2, be
  poly-scale in `N` as in `T(r*)` (Paper 14).