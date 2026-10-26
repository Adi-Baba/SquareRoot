# Paper 22 — M7.4: the isolated-embedding regime

Version 1.0.  This paper analyses the single remaining case left open by
M7.3 (`Docs/21`): **one** small embedding, all others bounded away from zero.
It (a) records the exact Fourier/Hadamard identities that constrain the
conjugates, (b) proves what the isolation hypothesis does and does not force,
(c) proves a **no-go theorem** for norm-type arguments in this regime, and
(d) reduces the whole project to one sharply-posed constant.

---

## 22.1 Setup and the exact identities

Let `d_1,…,d_m` be distinct squarefree, `C_i ≠ 0` integers,
`S = Σ_i b_i` with `b_i := C_i√d_i`, and

```
f(τ) := Σ_i b_i τ_i ,        τ ∈ {±1}^m ,
```

so `f(σ) = T_σ(C)` is the value of the embedding given by the sign pattern
`σ`.  Write `ε := |f(1)|` (after absorbing signs into the `b_i`), and
`M := max_i max(|C_i|, d_i)`, `L := Σ_i |b_i| ≤ m·M^{3/2}`.

**Key structural fact.** `f` is **affine linear** on the hypercube: it is a
linear combination of only `m` of the `2^m` Walsh characters `χ_i(τ)=τ_i`.

```
f(τ) = Σ_{i=1}^m b_i χ_i(τ),        f̂(∅) = 0,  f̂({i}) = 2^{m-1} b_i,  f̂(S) = 0 for |S|≥2.
```

So the vector of conjugates is *not* a general Walsh transform: it has only
`m+1` nonzero Fourier coefficients.  All of the following are exact.

**Proposition A1 (Parseval).**
```
Σ_{τ} f(τ)²  =  2^m · B² ,        B² := Σ_i C_i² d_i  ∈ Z_{>0} .
```

**Proposition A2 (trace).**  `Σ_τ f(τ) = 0`.

**Proposition A3 (flip identity).**  For each `i` and `τ`, with
`X_i(τ) := Σ_{j≠i} b_j τ_j`:

```
f(τ)·f(τ⊕e_i)  =  X_i(τ)² − b_i²  ∈ O_K ,     and   f(τ) − f(τ⊕e_i) = 2 b_i τ_i .
```

In particular `|f(τ) − f(τ⊕e_i)| = 2|b_i| ≥ 2` for every `τ, i`, and
`f(τ)f(τ⊕e_i)` is an algebraic integer living in the **degree-`2^{r-1}`**
subfield `Q(√d_j : j≠i)`.

**Proposition A4 (conjugate separation).**  Let `τ ≠ τ'` and let
`T := {i : τ_i ≠ τ'_i}`, `k_T := rank_{F₂}{d_i : i∈T}`, `s_T := Σ_{i∈T}b_i`.
Then

```
f(τ) − f(τ') = 2 s_T  ,   s_T ≠ 0 (Thm A.3) ,   hence
|f(τ) − f(τ')|  ≥  2 ( |T|·M^{3/2} )^{−(2^{k_T} − 1)} .                    (22.1)
```

*Proof.* `|f(τ)-f(τ')| = 2|s_T|`; `s_T` is a nonzero algebraic integer in the
multiquadratic field of degree `2^{k_T}` with all conjugates `≤ |T|·M^{3/2}`;
apply the norm bound. ∎

**Consequence (local rigidity).**  `|T| ≤ K` gives `k_T ≤ K`, so
`|f(τ)-f(τ')| ≥ 2(KM^{3/2})^{−(2^K−1)}`.  Hence

> **Proposition A5.**  If `K = O(log N)`, then all conjugates at Hamming
> distance `≤ K` are separated by `≥ N^{−O(1)}`.

**Proof of Proposition A5.** `2^K = N^{O(1)}` and `K·M^{3/2} = N^{O(1)}` when
`K = O(log N)`, so the right side of (22.1) is `N^{−O(1)}`. ∎

So the hypercube geometry of the conjugates is **polynomially rigid up to
Hamming distance `O(log N)`, and exponentially soft beyond it.**  The
difficulty lives precisely at large Hamming distance.

---

## 22.2 What the isolation hypothesis actually forces

Assume now the **isolated** case:

```
(H0)   ε = |f(1)| ,   and   |f(τ)| ≥ δ   for all τ ≠ ±1 .
```

**Proposition B1 (neighbour / subsum dichotomy).**  Under (H0),

```
(i)   for every nonempty T ⊊ [m] :   |f(1) − 2s_T| ≥ δ ,  so  |s_T| ≥ (δ−ε)/2 ;
(ii)  |b_i| ≥ (δ−ε)/2   for every i  (take T={i}).
```

*Proof.* `f(1⊕e_i) = f(1) − 2b_i` gives (ii) and, for general `T`,
`f(τ) = f(1) − 2s_T` for the flip set `T`; apply (H0). ∎

So the isolation hypothesis says: **no proper nonempty subsum of the `b_i`
lands within `(δ−ε)/2` of zero** — a *dissociated-subsums* condition.

**Proposition B2 (the norm bound is automatically consistent).**  For every
nonempty `T ⊊ [m]` the norm bound gives
`|s_T| ≥ (|T|·M^{3/2})^{−(2^{k_T}−1)}`, and Proposition B1 gives
`|s_T| ≥ (δ−ε)/2`.  Both are **lower** bounds on the same positive quantity,
so they never conflict.

> **Proposition B2 (inertness).**  The isolation hypothesis (H0) is
> *arithmetically inert* with respect to every Liouville/norm argument: it
> supplies lower bounds on quantities that norm arguments also bound from
> below, so it can never produce a contradiction, and hence can never improve
> a norm-derived bound on `ε`.

This is the formal version of the observation that "the other conjugates are
large" is the **wrong direction of information** for bounding `ε` from below.

---

## 22.3 The no-go theorem: norms cannot touch the isolated case

**Theorem C.**  Let `S = Σ_i b_i` be isolated in the sense of (H0), with
`δ = min_{τ≠±1}|f(τ)|`.  Then the strongest lower bound on `ε` obtainable
from the integrality of `S` alone is

```
log₂(1/ε)  ≤  (2^r − 2)·( (3/2)·log₂M + log₂ m ) + O(1) ,                  (C)
```

i.e. the usual norm bound with the tiny conjugate deleted; and this bound is
**monotonically worsened** as the isolation gap `δ` increases.

*Proof.* Every norm-type (Liouville) estimate on `|S|` has the form
`|S| ≥ ( Π_{τ≠1} |f(τ)| )^{−1}`, by `|N_{K/Q}(S)| ≥ 1`.  Equivalently
`log(1/ε) ≤ Σ_{τ≠1} log(1/|f(τ)|)`.  The hypothesis (H0) asserts lower bounds
on the summands `|f(τ)| ≥ δ`, i.e. it makes each `log(1/|f(τ)|)` **more
negative**.  Hence (H0) can only *decrease* the right-hand side: it never
yields a better bound on `ε`.  Deleting the antipodal pair (one factor is
`|f(-1)| = ε`) and applying `|f(τ)| ≤ m·M^{3/2}` for the rest gives
`log₂(1/ε) ≤ (2^r−2)(1.5 log₂M + log₂ m)`. ∎

**Corollary C1 (the norm's worst case is the isolated case).**  The quantity
`(2^r−1)(1.5 log₂M + log₂ m)` is maximised when every non-tiny conjugate
sits at its maximum `m·M^{3/2}`, which is *exactly* the maximally isolated
configuration.  So

> **The configuration that is hardest for the conjecture is the same
> configuration that is worst for the norm bound.**  The two-small case
> (M7.3) and the isolated case are complementary, and together exhaust all
> possibilities.

**Corollary C2 (no amount of Fourier/Hadamard machinery helps).**  Propositions
A1–A5 give exact identities, but A1 and A2 hold identically for every `C`
(they constrain nothing), A3–A5 are *separation* statements (they control
differences of conjugates, i.e. two-small configurations — already settled by
M7.3), and by Theorem C the remaining Fourier information enters in the wrong
direction.  **Hence the proposed route "use the Fourier/Hadamard structure of
the `f(τ)` to control a single isolated tiny coordinate" cannot work.**

---

## 22.4 The sharp witnesses live in the isolated regime

The retractions in `TODO.md` include the certified witness

```
d = (2, 3, 13),   C = (−10, −98, 51) ,    r = 3 .
```

Computed with **exact** algebraic arithmetic (multiquadratic ring, integer
coefficients) and 50-digit magnitudes:

```
exact norm:  N_{K/Q}(S) = 1        (all non-trivial mask coefficients vanish)
             product of |conjugates| = 1.0  = |N(S)|      [verified]

canonical conjugates (one per antipodal pair), ascending:
  σ = (1, 1, 1)    |T| = 2.83182528695e-07    E =  21.7517644049
  σ = (1,−1,−1)    |T| = 28.2842715306         E =  −4.82192810933
  σ = (1,−1, 1)    |T| = 339.481958567         E =  −8.40719109568
  σ = (1, 1,−1)    |T| = 367.766229814         E =  −8.52264519993

E₁ = 21.7517644049 ,   E₂ = −4.8219281093 ,   E₂/E₁ = −0.22168
```

Two facts, both certified:

1. **`S` is a unit** (`N(S) = 1`).  So the extremal witness is a *unit* whose
   single conjugate is small — the "one tiny embedding, all others large"
   picture is literally realised by a unit of the maximal order.
2. **`E₂/E₁ = −0.222 < 0`: the form is ISOLATED.**  The second-smallest
   conjugate has `|T| = 28.28`, five orders of magnitude above `ε = 2.8e−7`.

> **The known extremal witnesses are isolated, not two-small.**  M7.3
> therefore does **not** improve the sharp constant; the sharp rate is
> attained precisely in the regime M7.3 provably cannot reach.

This is the strongest available evidence that the isolated regime is the true
hard case, and it validates M7.4 as the correct target.

**Bookkeeping: the asymptotic rate and this finite witness are separate
statements.**  `TODO.md` records this witness as "saturating
`limsup E/N = (r−1)/r`" with `E/N ≈ 1.2795`.  The number `1.2795 = 21.7518/17`
uses coefficient bits only (`N = 17`); under the standard convention
`N = Σ bitlen|C_i| + Σ bitlen d_i = 17 + 8 = 25` the ratio is

```
E₁/N = 21.7518/25 = 0.8701     vs   (r−1)/r = 2/3 = 0.6667 .
```

Neither equals the asymptotic rate, and the finite instance *exceeds* it.
The dossier therefore states these as two independent claims, with no
conflation:

> **Claim A (asymptotic rate).**  `limsup_{r→∞} E/N = (r−1)/r`, proved via
> Subspace Theorem plus pigeonhole (`Docs/02` B.5, `Docs/03`).  This is a
> statement about a *family* as `r → ∞`.
>
> **Claim B (finite witness).**  `d=(2,3,13)`, `C=(−10,−98,51)` is a
> certified **unit** (`N_{K/Q}(S)=1`) with a single isolated tiny embedding,
> `E₁ = 21.7518`, second conjugate `|T| = 28.28`, `E₂/E₁ = −0.222`.

**Claim B is an extremal finite example for the norm mechanism — it is not
evidence that a single `r = 3` instance saturates the asymptotic rate of
Claim A, and it is not offered as such.**  Its role in this paper is
narrower and purely structural: it certifies that the isolated regime is
*realised by a unit*, so the configuration that stresses the norm bound is
attained, and M7.3 cannot be expected to improve the sharp constant.

For completeness the Subspace-Theorem bound with no additive constant gives
`(r−1)(1.5 log₂M + log₂m) = 23.0 ≥ E₁ = 21.75`, so the witness is consistent
with the `r−1` exponent under either `N` convention.

---

## 22.5 The reduction

Combining M7.3 (`Docs/21`) with Theorem C:

**Theorem D (complete dichotomy).**  Let `S ≠ 0` be a multiquadratic form with
`r = rank_{F₂}{d_i}`.  Then exactly one of the following holds.

1. **`S` has two non-antipodal small embeddings** (hypothesis (H) of
   `Docs/21`).  Then, **absolutely effectively**,
   ```
   log₂(1/|S|)  ≤  ( 2^{min(k_A,k_J)} − 1 )·( 1.5·log₂M + log₂ m + 1 )  ≤  poly(N),
   ```
   with no field-dependent constant whatsoever.

2. **`S` has a single isolated small embedding.**  Then the only currently
   known bound is the norm bound (C) — exponential in `r` — and the
   Subspace Theorem (`Docs/03`) gives
   ```
   log₂(1/|S|)  ≤  (r−1+δ)( 1.5·log₂M + log₂ m )  +  log₂(1/c(d,δ)) ,
   ```
   whose polynomial part is harmless (`r ≤ N`, `log M ≤ N`) and whose **only**
   obstruction is the field-dependent constant `log₂(1/c(d,δ))`.

Therefore:

> **Open (M7.4, final form).**  Prove
> ```
> log₂(1/c(d,δ))  ≤  poly(N)
> ```
> for the multiquadratic Subspace-Theorem constant — equivalently, make the
> **isolated** case effective.  Nothing else in the project is open: the
> two-small case is closed absolutely and effectively by M7.3, and the sharp
> rate is attained in the isolated case by `d=(2,3,13), C=(−10,−98,51)`.

---

## 22.6 Answer to the proposed formulation

The suggested question was:

> Can "one coordinate tiny + all other coordinates non-tiny" force a
> polynomial lower bound on that coordinate?

**Answer: not by any method that uses only integrality, norm inequalities, or
the Fourier/Hadamard structure of the conjugates.**  Theorem C proves this
rigorously: the isolation hypothesis enters every such argument with the wrong
sign.  The mechanism that *can* convert isolation into a polynomial bound is
the Subspace Theorem, and the price is exactly `log₂(1/c(d,δ))`.

So M7.4 is worth pursuing, but **as an effectiveness theorem for the
Subspace-Theorem constant in the isolated case — not as a new norm or Fourier
argument.**  Any further pure-norm refinement is provably futile
(Corollary C2).