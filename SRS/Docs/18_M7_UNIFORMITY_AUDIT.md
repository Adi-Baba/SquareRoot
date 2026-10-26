# Paper 18 — M7: Uniformity Audit — what information is actually available?

Specification of the next honest objective.  This note *states* the mandate;
it does not yet prove anything new.

---

## 18.1 The target statement

For every rank `r`, every squarefree tuple `d=(d_1,…,d_m)`, every `C ≠ 0`
with `S = Σ C_i√d_i ≠ 0`, we want an effective lower bound containing a
Subspace-typed exponent along with an effective constant:

```
log₂(1/c(d,δ)) ≤ poly(r, N)        for all δ ∈ (0,1),
```

in

```
|L(C)| > c(d,δ)·‖C‖^{−(r−1+δ)},   L(C) = C_r + Σ_{i<r} C_i √(d_i/d_r).
```

This is exactly the statement that, on the full input `N`, the exponent and
the constant are *both* polynomial.

---

## 18.2 Everything presently exhausted (M1–M6 virtual repository)

| Route we were led from | Why it cannot close uniform SRS | Where documented |
|---|---|---|
| Full algebraic norm `N_{K/Q}(S)` | exponential `2^{r−1}−1` in `r` | `Docs/02_NORM_BOUND.md` |
| Relative-norm quadratic towers | telescopes to same full norm (M1 negative) | `Docs/10` |
| Simultaneous-approximation through sign-pattern matrices | M4-theory: does not compress the conjugate count | `Docs/13` |
| Effective fixed-rank | effective but the norm bound is exponential in `r` | `Docs/02`, `Docs/09` |
| Complete norm lower bound at r=3 | already strong (fixed-rank), but only for fixed `r` | `Docs/09_RANK3_NORM_FORM.md` |
| Unit-equation on `ρ = d₃u²` at r=3 | the Matveev constant `𝔇 ≤ 2^{O(N)}`, no polynomial scale; `T(3,δ)` unproved | `Docs/15`, `Docs/16` |
| The sufficient inequality `max_i h(ε_i) ≤ N^{O(1)}` (♦) | not a theorem for the `H=±1` fields   and of comparable difficulty | `Docs/16` M6.5 |

---

## 18.3 Three open information buckets

Every future mechanistic route must obtain at least one of the following
from SRS data alone — none is implied by the items already closed:

### Bucket 1: Special-field structure beyond `(Z/2Z)^r`

Concretely: a relationship among `{d_i}` beyond independence, between `C`
and the radicands, or in the *order* of `K`, that the full norm does not see.
Example targets to test:

* the shape of the norm-one units `ε_i` of the *rank-r* real subfields of
*specific* `H=±1` slices (the "M6-surviving" family from `Docs/16`);
* whether the discriminant `disc(K)` is constrained by `S ≠ 0` on the
small- conjugate hitting set;
* relative integrality/variants of `d`'s parity patterns beyond the square of
the GCD cube.

### Bucket 2: Simultaneous approximation that does not follow one linear form

The two-small-embedding lemma is a sharp one-shot: if any *two* real
embeddings of `L(C)` are `≤ ε`, then `ε` is itself bounded below by the
`max(2^{r_J}−1)`-scale of the difference form.  We must analyse the *class*
of coefficient vectors producing *many* simultaneous small embeddings — a
lattice condition.  The sign-difference matrix captures the F₂ semivariance of
the small set, and what we need is a lower bound on the minimum determinant
of the induced rational approximation type.

### Bucket 3: A different algebraic object

Any entire revision in which `S` is built around a generator whose
fine arithmetic (a unit, a discriminant, or a resultant of the `d_i`
and `C_i`) is not dictated by the sign-embeddings of `S` itself.  This
bucket is
intentionally empty and is where a genuine surprise would have to live.

---

## 18.4 First deliverable when M7 is resumed

Create the adversarial enemy for each bucket separately: given a route in
rows of §18.2, construct a *small* coefficient vector that realises its
maximum blurring: e.g. for full norm, choose `C_i sqrt(d_i)` phases so that
exactly one conjugate vanishes; for the two-small-embedding constraint,
choose a sign matrix with a two-dimensional support to saturate M4.
The aim is to make precise which *feature* of those enemies a hypothetical
uniform theorem would have to reject — i.e. the next "significant bit" of
the open SRS statement that current estimates all miss.

This is the job of M7.  It is explicitly **not** a search for another
generic Baker/Subspace application.  It is a classification: which
sub-object of the multiquadratic tower `Q ⊂ Q(√d_{i₁}) ⊂ … ⊂ K` can the
remaining magnitudes that current bounds (norm, tower, two-embedding, sign-pattern,
generic subspace, generic Baker) all, separately, regard as generic and
therefore cannot individually collapse.