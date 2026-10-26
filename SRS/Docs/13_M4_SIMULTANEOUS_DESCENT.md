# Paper 13 — Simultaneous Descent: exact sign-pattern/rank analysis (M4-theory)

Version 1.0.  This is the theoretical branch requested after Paper 12: if
`k` sign-embeddings of the same `S` are simultaneously small, what does the
sign-difference matrix actually buy, combinatorially and over a field?  The
answer is: a clean column-pattern decomposition, an exact rank count, and an
**obstruction** — sign-difference rank never compresses the exponent.  We
prove that, and record the one full-rank case, and then state exactly why the
open constant must be attacked through *simultaneous Diophantine approximation*
rather than through a descent of radical supports.

---

## 13.1 Setup with correct indexing

Let

```
S(C) = Σ_{i=1}^{m} C_i·√d_i,   d_i distinct squarefree ≥ 2,
```

let `r = log₂[K:Q] ≤ m` be the rank (the field degree), and let `σ` range
over the `2^r` real embeddings (the sign vectors that act on the independent
generators; dependent radicands have their signs forced).  Record the
quantity measured by the embedding profile:

```
T_σ := Σ_{i=1}^{m} C_i·σ(√d_i) ,        σ ∈ G_K.
```

Here `G_K` is the group of real embeddings of `K` (in bijection with sign
classes, and for generic `S` all give distinct values).  `T_σ` is one real
number per embedding.

For one fixed coefficient vector `C`, `T_σ` is one real number per choice of
sign pattern `σ`.  The event "many embeddings are small" is

```
|T_{σ₁}(C)| ≤ ε,   |T_{σ₂}(C)| ≤ ε,   …,   |T_{σ_k}(C)| ≤ ε.
```

---

## 13.2 From many small embeddings to a sign-difference matrix

Pick one small embedding as the base, say `σ⁰`.  For each of the other
`j = 1,…,k-1` small embeddings `σ^j`, form the difference row

```
D^{(j)} := σ^j − σ⁰   ∈ {−2,0,2}^m,
```

where we read the sign of `σ^j(√d_i)` relative to `σ⁰(√d_i)`.  Set

```
w := ( C₁√d₁, C₂√d₂, …, C_m√d_m )ᵀ  ∈ R^m.
```

Then the smallness bounds give, for every `j`,

```
|⟨D^{(j)}, w⟩| = |T_{σ^j} − T_{σ⁰}| ≤ 2ε .                              (∗)
```

Let `D` be the `(k-1)×m` real matrix with these rows.  Its only structural
property is in the columns.

---

## 13.3 The column-pattern decomposition (exact)

At each coordinate `i ∈ {1,…,m}`, the `(i)`-th column of `D` over the
`k-1` difference rows is the vector of relative signs of embedding
`σ^j(√d_i)` versus `σ⁰(√d_i)`:

```
col_i(D) := ( (σ^j(√d_i))/(σ⁰(√d_i)) − 1 )_{j=1..k-1}  ∈ {−2,0,2}^{k-1}.
```

So `col_i(D)` depends only on the **column-pattern** `c_i := col_i(D)`, and
coordinates whose columns are equal form a block

```
V_c := { i : col_i(D) = c } ⊆ {1,…,m}.
```

There are at most `2^{k-1}` distinct column-patterns, and for each
non-zero pattern `c` the contribution of `V_c` to the row constraints is
the same: in row `j`, the block `V_c` contributes `c_j · Σ_{i∈V_c} w_i`.
Define

```
B_c := Σ_{i ∈ V_c} C_i·√d_i .
```

Then the vector of differences `Dw` across the `k-1` rows is exactly

```
D w  =  Σ_{c ≠ 0}  B_c · c    (a linear combination of the column-pattern
                               vectors c ∈ {−2,0,2}^{k-1}).
```

**Lemma 13.1 (rank of the sign matrix).**

```
rank_R D  =  rank_R( { c ∈ {−2,0,2}^{k-1} : c ≠ 0 } ) ≤ min(k−1, m).
```

*Proof.* The columns of `D` are exactly the vectors `c_i` grouped into
constant blocks `V_c`, so the column space of `D` equals the span of the
distinct nonzero column-pattern vectors `c`.  There are at most
`min(k−1, 2^{k−1}−1, m)` of them, and they live in `R^{k-1}`, so the rank is
at most `k−1`. ∎

---

## 13.4 The exact content of many simultaneous small embeddings

The `k-1` inequalities in (∗) say, in the coordinates just introduced,

```
| Σ_{c ≠ 0} B_c · c_j | ≤ 2ε    for j = 1,…,k-1.
```

If the distinct nonzero patterns `c` are `c⁽¹⁾,…,c⁽q⁾` with
`q = rank_R D`, and we form the `q×(k-1)` matrix `P` whose columns are the
`c⁽ⁿ⁾`, then (∗) is the single statement that

```
|Pᵀ · B| ≤ 2ε   ⇒   B_c satisfies q independent small equations.        (♣)
```

Crucially, (♣) does **not** say each `B_c` is small.  Each `B_c` is an
integer radical form supported on the block `V_c`:

```
B_c = Σ_{i∈V_c} C_i √d_i .
```

It is a nonzero algebraic integer as soon as the signs on `V_c` are not all
zero in any radical dependency; its absolute value jumps in `[1/d, ∞)` only
after normalizing by the denominator of the radicand product, and
`|B_c| ≥ 1/(denom)` once we clear radicals.  But the individual `B_c` can
**stay at scale `M^{1/2}` while their *combination* through `Pᵀ` is small**:
that is precisely the extra condition a multi-embedding configuration imposes.
The equations (♣) are the content of "many small embeddings"; reducing them
to a *single* short radical relation would require inverting `Pᵀ` on the
vector `B`, i.e. knowing all the `B_c` independently — which is the thing we
do not know.

---

## 13.5 From a *short support* relation to M4.1 (the single clean branch)

If a single column pattern `c` is *nonzero*, the row on `V_c` becomes
(for k=2, or by isolating one row):

```
X_c := Σ_{i∈V_c} (c_i/2)·C_i·√d_i,   |X_c| ≤ 2ε .
```

Each coefficient `(c_i/2)·C_i` is an even integer; `X_c` is a nonzero
radical form on `V_c`.  Its effective norm cost is the rank of the radicand
subset on `V_c`, call it `r_c = log₂[Q(√d_i:i∈V_c):Q]`, and M4.1 applied to
`X_c` gives

```
ε ≥ (2·( 2·r·M^{3/2} )^{2^{r_c} − 1})^{-1}.
```

So a **short-support** difference (one block with a single radical) is the
clean descent: it forces `varepsilon` down to a rank-`r_c` norm threshold with
`r_c` equal to the number of *independent radicals on that block*.

---

## 13.6 The obstruction: rank does not imply short support

Two sign patterns `σ⁰, σ` can differ on a *large* block (every coordinate
object if the radicands carry full rank), in which case `r_c` can be as large
as `r`.  In that generic case the norm cost of the corresponding `X_c` is the
**full** one – `(2^{r}−1)` – and no support-based descent occurs: the
coefficient vector does not concentrate on fewer generators.

Even collecting all `k-1` difference rows does not structurally reduce the
support: by Lemma 13.1, `rank D ≤ k−1`, so to cut the real vector
`(C_i√d_i)` down to the origin by *these* equations one needs

```
rank D = m   ⇒   k−1 ≥ m ≥ r.
```

Only when `k−1 ≥ m` and the column-pattern vectors are independent can
`D` have maximal rank.  But then `D w ≈ 0` is an `m×m` integer system
inverting `w = (C_i√d_i)`, which collapses – via `M^{-1}` on the
minors – to `|C_i√d_i| ≤ (measure)·ε` for every i.  In other words: *the
only way many small embeddings of the same `C` can *fully* determine `C` is
when the number of distinct small embeddings `k` is at least `m+1`, and then
the bound on `C` is as exponential in `m` as the original norm bound.*

---

## 13.7 The two clean cases and their consequences

**Case A (rank exploratory, `q < m`).**  If only `q ≤ k−1 < m` independent
column patterns occur, then the smallness of `B` enters only the finite
projection `PᵀB`.  The complement (`m − q ≥ 1` free directions) means there
exists an integer vector `x ∈ R^m` with `x` in the kernel of `Pᵀ` for which
`C`-combinations are undisturbed; formally, nothing single-handedly improves
the norm exponent, and M4.1 applied to each individual relation gives the
**same** exponent `2^{r_c}-1` as for that one relation.  So less than `m`
embeddings cannot by themselves trigger a genuine rank drop.

**Case B (rank saturation, `q = m`).**  If `rank D = m`, then `D` is a
full-column-rank integer `(k-1)×m` system, and `min ‖D w‖_∞` over nonzero
`C` with `S≠0` is achieved on an open (discrete) set.  It cannot be zero
(`D w = 0` implies `w = 0` since `rank D = m`, hence all `C_i = 0`).  So the
smallest value is strictly positive, and whenever `2ε` is below that
minimum, no nonzero integer `C` can have all `k` embeddings simultaneously
`2ε`-small.  The search for that minimum involves precisely the same
`2^r`-factor norm-type product as the full norm (`m` independent radical
constraints), so the certified threshold cannot be better than the full norm
bound.  This case therefore reproduces the known exponent.

**Conclusion 13.2.**  In every case accessible to this sign-difference linear
algebra, the effective lower bound on `ε` it can certify is dominated by the
norm cost of one radical form; those costs are
`(2^{r_c}−1)` or `2^{r}−1` with `r_c ≤ r`.  Sign matrices can mix the
differences, but they cannot, by themselves, produce a form whose *norm cost*
is smaller than the norm cost of at least one difference form.  Therefore
*support-based descent does not beat the conjugate-count exponent.*

---

## 13.8 Exact equivalent statement (the actual open content)

The obstruction can be written as one clean algebraic demand.  Given any `C`
and any embedding `σ⁰`, define the **neighbourhood**

```
N_C(σ⁰,δ) := { σ ∈ {±1}^m : |T_σ(C) − T_{σ⁰}(C)| ≤ δ } .
```

If `|N_C(σ⁰, δ)| = k`, then the k points of the neighbourhood produce the
`k×m` difference matrix above, of rank `q ≤ k−1`.  There is **no** algebraic
identity that, given only these `k` small differences, lowers the effective
norm cost: the induced constraints live in `ker(Pᵀ)` directions, each
grounded in a rank-`≤r` sub-norm.*

  ────────────────────────────────────────────────────────────
  * precise means: for each of the independent column patterns, one has a
  small integer radical integer form whose own norm threshold is
  `2^{r_c}−1`.  A *smaller* threshold would require a shorter column pattern,
  which does not come for free from more embeddings: the columns *are* the
  independent data.  ∎

Consequently the missing object is not a bigger sign matrix: it is
*external simultaneous Diophantine information about the real vector
`(√d_i)_i` of degree `r−1`*.  That is the exact statement we cannot derive
from the sign geometry and which carries the full content of the open
uniform constant.

---

## 13.9 What this doc freezes

- Rank of the difference matrix equals the rank of the nonzero column-pattern
  vectors, `≤ min(k−1, m)` (Lemma 13.1).
- The `k` smallness conditions are equivalent to a smallness of a
  `(rank)`-dimensional projection of the vector `(B_c)_c`; they do **not**
  bound a single `B_c` (§13.4).
- M4.1 already gives `ε ≥ threshold` pedometrical for each individual
  constituent form on its support, with exponent `2^{r_c}−1` — no improvement
  with `k` (§13.5, §13.7 Case A/B).
- Therefore: multi-embedding smallness **by itself** never yields a lower
  effective exponent than one small embedding on the worst support.  Any
  genuine exponent reduction must import an effective simultaneous
  approximation result on `(√d_i)`.  This is the precise, proved boundary
  of the embedding-geometry attack.

*(No numeric claim in this file is asserted; all statements are formal.)*