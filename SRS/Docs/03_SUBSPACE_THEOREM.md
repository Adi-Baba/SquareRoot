# Paper C — The Subspace Theorem Bound

Version 1.0. States what the Subspace Theorem gives for SRS and exactly why
that constant is ineffective. We do not re-prove the Subspace Theorem (it is
quoted from the literature); we do the specialization in full and we say
plainly what it does and does not buy.

---

## C.0 The single form that SRS actually is

Fix radicands whose `F₂`-rank is `r`, and let `d_1, …, d_r ≥ 2` be a
spanning independent set. Put

```
α_i = sqrt(d_i / d_r)   (i < r).
```

Then `√d_1, …, √d_r` are `Q`-independent (Paper A, Cor A.2 with independent
parity vectors), and from `Σ_{i<r} c_i·α_i = q ∈ Q` we get
`Σ c_i √d_i = q·√d_r`, forcing `c_i = q = 0`. Hence
`1, α_1, …, α_{r−1}` are `Q`-linearly independent.
Define the single linear form

```
L(C₁,…,C_r) = C_r + Σ_{i<r} C_i·α_i.
```

Then `S = √d_r · L`, since
`√d_r·(C_r + Σ_{i<r} C_i √(d_i/d_r)) = C_r√d_r + Σ_{i<r} C_i √d_i`.
So `L` is exactly `S` up to the non-zero factor `√d_r`, and `|S|` and `|L|`
have the same rate.

---

## C.1 Theorem C.1 (the exponent `r − 1`)

**Theorem C.1 (Subspace Theorem, specialized).** For every `δ > 0` there are
a constant `c = c(d₁,…,d_r, δ) > 0` and finitely many proper rational
subspaces `V₁, …, V_t ⊊ Q^r` such that for every integer vector `C` not in
any `V_j`,

```
|L(C)| > c · ‖C‖^{-(r−1+δ)},
```

and consequently, with `m = log₂‖C‖_∞` and fixed radicands (`N = r·m`),

```
log₂(1/|S|) ≤ (r−1+δ)/r · N + log₂(1/c) + O(1).
```

*Proof.* Apply W. M. Schmidt's Subspace Theorem in `r` variables to the `r`
forms

```
L_0 = L,    L_1 = C₁,    L_2 = C_2,   …,   L_{r−1} = C_{r−1}.
```

These are `r` forms in `r` variables, linearly independent over `ℚ̄` because
`1, α_1, …, α_{r−1}` are `ℚ`-independent (Paper A). The Subspace Theorem
says: for every `δ > 0`, every integer point `C` with

```
|L(C)| · Π_{i<r}|C_i| ≤ ‖C‖^{-δ}
```

lies in a finite union of proper rational subspaces `V₁,…,V_t ⊊ Q^r`. Since
`Π_{i<r}|C_i| ≤ ‖C‖^{r−1}`, every `C` **outside** that union satisfies

```
|L(C)| > ‖C‖^{-δ}/‖C‖^{r−1} = ‖C‖^{-(r−1+δ)}.
```

It remains only to identify what the points *inside* some `V_j` contribute.
Each `V_j ∩ Z^r` is a sublattice of rank `≤ r−1`. Parametrise it by a basis
`b_1,…,b_s` (`s ≤ r−1`, integer matrix with rational entries): `C = Σ a_k b_k`.
Restricted to `V_j`, the form `L` becomes a linear form
`L^{(j)}(a) = Σ_k (Σ_i b_{k,i}·c_i)·a_k` in `s` integer variables. Its
coefficients are `ℚ`-linear combinations of `1, α_1,…,α_{r−1}`. If `L^{(j)}`
vanishes identically on `V_j`, then `V_j ⊆ ker L`, i.e. `L` is identically
zero on a rank-`≥1` lattice, which is a `ℚ`-relation among the coefficients of
`L` — impossible unless `C` is the zero vector, which is excluded. So on each
`V_j`, `L^{(j)}` is a non-zero linear form, of the same multiquadratic shape
(its non-constant coefficients live in the same degree-`2^r` field), in
`s ≤ r−1` integer variables, with coefficients effectively bounded by
`‖C‖^{s}` and with input size `N_j ≤ N`.

Apply the same argument to `L^{(j)}` on `V_j`. Either the resulting finite
exception set handles `V_j`-points of `L^{(j)}` outside a smaller union, or
`L^{(j)}`'s exceptional subspaces cut `V_j` further. Each induction step
strictly lowers the variable count, so after at most `r` rounds every integer
point is either in a terminal proper subspace where `L` is identically zero
(hence `S = 0`, excluded) or satisfies the bound. Collecting the (finite,
`d`-dependent) constants from each round into one `c(d,δ)` gives

```
|L(C)| > c · ‖C‖^{-(r−1+δ)}.
```

This recursion explains *why* `c` is dynamically ineffective: at each stage we
must know the exceptional subspaces, and the Subspace Theorem provides only
their existence, not a computable description. ∎

---

## C.2 The exact deficit: why `c` is ineffective

The constant `c` produced by the Subspace Theorem is **not computable from**
`d, δ` alone. The proof is a contradiction argument: it assumes an infinite
sequence of solutions and derives that it is impossible; it does not locate
the solutions. Hence no a priori upper bound on `log₂(1/c)` comes out.

So the correct reading of Theorem C.1 is the **asymptotic** statement: for
every `ε > 0`,

```
E/N ≤ (r−1)/r + ε    (as N → ∞, radicands fixed).
```

Combined with the pigeonhole lower bound (Paper B, Thm B.5), this gives the
sharp rate

```
limsup_{N→∞} E/N = (r−1)/r.
```

What it does **not** give is a uniform-in-`d` polynomial `log₂(1/c)`. That is
the single missing ingredient, and it is exactly EHS Problem 33.

---

## C.3 Comparison with Eisenbrand–Haeberle–Singer

EHS (arXiv:2312.02057, SoCG 2024) apply the Subspace Theorem to the
simultaneous-approximation form of the problem and obtain

```
E ≤ (2n)·N        (EHS Thm 8),
```

i.e. their exponent is `2n` rather than `r−1`. Our Theorem C.1 is **sharper by
an exponential factor in `r`** because it applies the Subspace Theorem to the
single linear form that SRS actually is, not to a simultaneous system. The two
estimates agree on effectivity: both are ineffective in the constant, so
neither settles the conjecture.

---

## C.4 What would make this effective

A uniform bound `log₂(1/c) ≤ poly(N)` over all `d, δ`. Equivalent statements:

```
|L(C)| ≥ 2^{−poly(N)}   uniformly over all integer C and all squarefree d,
```

or `c(d,δ) ≥ 2^{−poly(N)}`. This is the exact open statement of Paper F.
Everything in this file is conditional on it, which is why we state it last.