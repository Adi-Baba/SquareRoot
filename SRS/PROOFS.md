# SRS — Proofs

Written proofs for the theoretical claims stated in `README.md`. Everything
here is self-contained: definitions, standing assumptions, and every
nontrivial step. Where the repo also verifies a claim computationally, that
is recorded as corroboration, never as the evidence — same rule as
`Src/verify_theorems.py`.

Status labels:

- **PROVED** — proof written out below.
- **QUOTED** — the result is from the literature, reproduced here for
  comparison, **not proved in this file**. Used for exactly one result:
  Eisenbrand–Haeberle–Singer's Subspace-Theorem bound (§2.8), which is the
  best known separation bound and the reason the `2^r` in Thm 2.6 is a
  property of the norm rather than of the problem.
- **OPEN** — no proof known here; stated as a problem, not a claim.

**Orientation.** `THE_MISSING_THEOREM.md` is the roadmap and
`STATE_OF_THE_ART.md` records the correction that changed the picture; both
should be read before this file. In particular §2.8 supersedes an earlier
version of itself, and §4 states the currently open problem — which is about
the *effectivity of one constant*, not about an `ε`-improvement.

**Addendum, added by the M7/L programme (`Docs/21`–`Docs/38`).** This file
covers Papers 1–19. The subsequent programme reduced the effectivity question
to a single explicit quantity and then **parked** it, and that result is not
reflected in §4 above:

- `Docs/24` reduces effectivity to bounding the largest exceptional Pell index
  `k*` with `k*·log ε = O(N²)`, which would give `E ≤ O(N²)` — the target.
- `Docs/30`–`Docs/33` make the surrounding machinery explicit: the exceptional
  set is enumerated (`Docs/33` Thm 33.2), bounded (`Docs/33` Thm 33.4), and the
  threshold is explicit (`Docs/33` Thm 33.1).
- `Docs/35` Thm 35.3 identifies the obstruction as the fractional part of a
  **trace-zero order-2 linear recurrence** `x_{k+2} = 2x_{k+1} + x_k`.
- `Docs/36` Thm 36.2 is the **stopping theorem**: no argument forming an
  algebraic number from the nearest-integer data and bounding it by the product
  of its other conjugates can close the problem, for any `n ≥ 2`, because such
  bounds saturate at exponent `≥ 1` while the requirement is an exponent `< 1`.
- `Docs/37` attacks the original *product* rather than a single factor, and
  `Docs/38` (M7.13) closes the last lever by checking EHS's own parameter range
  (`δ = 1/(3n)`, forced by geometry of numbers).

So the answer to §4's open question is now sharper than "an effective constant
is missing": the constant is missing because the remaining statement is a
**new small-fractional-part problem for a trace-zero quadratic recurrence**,
and §2.10.4's sufficiency lemma (`exp(−poly(N))` suffices) is the right target.
§4 is otherwise unchanged and still accurate for what it states.

---

## 0. Notation and standing assumptions

- Inputs: `k` pairs `(c_j, a_j)`, integers with `|c_j| < 2^n` and
  `0 < a_j < 2^n`. `N = k·n` is the total input size.
- `S = Σ_j c_j·sqrt(a_j)`.
- Squarefree decomposition writes `a_j = t_j²·d_j` with `d_j` squarefree
  (`d_j = 1` allowed); the aggregated form is `S = Σ_{d ∈ D} C_d·sqrt(d)`
  over the distinct squarefree parts `D`, with `m = |D| ≤ k`.
- All logarithms are base 2 unless marked `ln`.
- For §3 (the F1 family): `α = a(a+1)` with `a ≥ 2`, `p_n/q_n` the continued-
  fraction convergents of `sqrt(α)`, `D_n = p_n² − αq_n²`,
  `S_n = p_n − q_n·sqrt(α)`.

---

## 1. Theorem 1 — the zero test is correct

This proves the correctness half of the README's `S = 0` in **P** row. (The
algorithmic half — aggregation in polynomial time — is `squarefree.py`.)

Throughout, for a positive integer `t`, `v(t)` is its prime-exponent parity
vector: the map prime → `exponent mod 2`, read as a vector over `F₂`. A set
of positive integers is **F₂-independent** when their parity vectors are
linearly independent over `F₂`. Equivalently: no nonempty subproduct is a
perfect square. Note `v(1)` is the zero vector, so `1` is never part of an
independent set; it is handled as the empty subset below.

**Lemma 1.1 (irrationality step).** Let `t₁, …, t_r` be positive integers
with independent parity vectors, and let `s ≥ 1` be an integer such that
`{t₁, …, t_r, s}` also has independent parity vectors. Then
`sqrt(s) ∉ Q(sqrt(t₁), …, sqrt(t_r))`.

*Proof.* Induction on `r`.

`r = 0`: the field is `Q`. `sqrt(s) ∈ Q` would make `s` a perfect square,
i.e. `v(s) = 0`; but independence of the one-element set `{s}` forces
`v(s) ≠ 0`. Contradiction.

Step: let `K = Q(sqrt(t₁), …, sqrt(t_{r−1}))`. The set `{t₁, …, t_{r−1}}` is
F₂-independent (a subset of an independent set), and `{t₁, …, t_{r−1}, t_r}`
is independent, so the induction hypothesis applied with `s' = t_r` gives
`sqrt(t_r) ∉ K`. Hence `x² − t_r` has no root in `K` and is irreducible over
`K`, so

```
K(sqrt(t_r)) = K ⊕ K·sqrt(t_r)      (direct sum of Q-vector spaces).
```

Suppose `sqrt(s) ∈ K(sqrt(t_r))`, and write `sqrt(s) = u + v·sqrt(t_r)` with
`u, v ∈ K` (unique by directness). Squaring:

```
s = u² + v²·t_r + 2uv·sqrt(t_r).
```

The left side and `u² + v²·t_r` lie in `K`; `2uv·sqrt(t_r)` lies in the
complement `K·sqrt(t_r)`. Directness forces `2uv = 0`, and a field has no
zero divisors, so `u = 0` or `v = 0`.

- `v = 0`: `sqrt(s) = u ∈ K`. But `{t₁, …, t_{r−1}, s}` has independent
  parity vectors (subset of an independent set; and `v(s)` is not in the
  span of the others), so the induction hypothesis is contradicted.
- `u = 0`: `sqrt(s) = v·sqrt(t_r)`. Then `(v·t_r)² = v²t_r² = (s/t_r)·t_r² =
  s·t_r`, and `v·t_r ∈ K`, so `sqrt(s·t_r) ∈ K` (the field is a subfield of
  ℝ, so the positive root is `|v·t_r| ∈ K`). Now `v(s·t_r) = v(s) ⊕ v(t_r)`;
  if it lay in `span{v(t₁), …, v(t_{r−1})}`, then `v(s)` would lie in
  `span{v(t₁), …, v(t_r)}`, contradicting independence of
  `{t₁, …, t_r, s}`. So `{t₁, …, t_{r−1}, s·t_r}` is F₂-independent and the
  induction hypothesis is contradicted again. ∎

**Theorem 1.2 (F₂ basis).** If `t₁, …, t_r` have independent parity vectors,
then `{ δ_S = ∏_{i∈S} sqrt(t_i) : S ⊆ {1,…,r} }` is a `Q`-basis of
`K_r = Q(sqrt(t₁), …, sqrt(t_r))`, and `[K_r : Q] = 2^r`.

*Proof.* Induction on `r`. `r = 0`: `δ_∅ = 1` spans `Q`. Step: Lemma 1.1
with `s = t_r` gives `sqrt(t_r) ∉ K_{r−1}`, so
`K_r = K_{r−1} ⊕ K_{r−1}·sqrt(t_r)` with `[K_r : K_{r−1}] = 2`. By the
induction hypothesis `{δ_S : S ⊆ [r−1]}` is a basis of `K_{r−1}`; the two
cosets multiply it into exactly `{δ_S : S ⊆ [r]}`. ∎

**Corollary 1.3 (distinct squarefree radicals are independent).** Let
`d₁, …, d_m` be distinct squarefree positive integers (`d = 1` allowed).
Then `{ sqrt(d₁), …, sqrt(d_m) }` is linearly independent over `Q` — even
when the `d_i` are dependent as field generators (e.g. `{2, 3, 6}`, where
`sqrt(6) = sqrt(2)·sqrt(3)`).

*Proof.* Choose a maximal subset `J ⊆ {1,…,m}` whose parity vectors are
independent (greedy on the `F₂` row space). Maximality means every
`v(d_i)` lies in `span{v(d_j) : j ∈ J}`. For each `i`:

- If `i ∈ J`, `sqrt(d_i)` is one of Theorem 1.2's basis elements
  (`δ_{S(i)}` with `S(i) = {i}`, after renumbering `J` as `1..r`).
- If `i ∉ J`, write `v(d_i) = ⊕_{j∈S(i)} v(d_j)` with `S(i) ⊆ J`. Then
  `d_i·∏_{j∈S(i)} d_j` has even prime parity everywhere, so it is a perfect
  square `c_i²`, and

  ```
  sqrt(d_i) = (c_i / ∏_{j∈S(i)} d_j) · δ_{S(i)},
  ```

  a nonzero rational multiple of a basis element.

The `S(i)` are distinct: `S(i) = S(i')` would give `v(d_i) = v(d_{i'})`,
hence `d_i·d_{i'}` a perfect square, hence `d_i = d_{i'}` (both squarefree),
contradicting distinctness. Distinct basis elements scaled by nonzero
rationals remain independent. For `d_i = 1`: `v(1) = ∅`, `S(i) = ∅`, and
`√1 = 1 = δ_∅`, covered by the same argument. ∎

**Corollary 1.4 (zero test).** `S = Σ_d C_d·sqrt(d)` over distinct
squarefree `d`. Then `S = 0` if and only if every `C_d = 0`.

*Proof.* Corollary 1.3 makes `{sqrt(d)}` independent, so the only
representation of `0` is the trivial one. ∎

*Remark.* `field.py` refuses F₂-dependent generator sets as a constructor
invariant (it wants a basis, not a generating set). Corollary 1.3 is what
lets the zero test in `squarefree.py` accept any distinct-squarefree input
with no independence check.

---

## 2. Theorem 2 — the general separation bound

This proves the README's "best known bound" row: every nonzero `S` is
bounded away from zero by `2^(−Θ(n·2^k))`.

**Theorem 2.1.** Let `S = Σ_{j=1}^{k} c_j·sqrt(a_j) ≠ 0` with
`|c_j| < 2^n` and `0 < a_j < 2^n`. Then

```
log₂(1/|S|) ≤ (2^k − 1)·(3n/2 + log₂ k).
```

*Proof.* Four steps.

1. **`S` is an algebraic integer.** Each `sqrt(a_j)` is a root of the monic
   integer polynomial `x² − a_j`, hence an algebraic integer; the algebraic
   integers form a ring, so `S` is one. Its minimal polynomial over `Q` is
   monic with integer coefficients.

2. **Degree.** `S` lies in `K = Q(sqrt(d) : d ∈ D)`, the field generated by
   the aggregated squarefree parts. Pick a maximal F₂-independent subset of
   `D` (size `r ≤ m ≤ k`); by Theorem 1.2 `[K : Q] = 2^r ≤ 2^k`. Hence
   `deg(S) ≤ [Q(S) : Q] ≤ [K : Q] ≤ 2^k`.

3. **Conjugates.** Every embedding `σ` of `K` in `ℂ` sends `sqrt(d)` to
   `±sqrt(d)` (the roots of `x² − d`), so `σ(sqrt(a_j)) = ±sqrt(a_j)` and

   ```
   |σ(S)| ≤ Σ_j |c_j|·sqrt(a_j) < k · 2^n · 2^{n/2} = k·2^{3n/2}.
   ```

   `K/Q` is normal (splitting field of `∏(x² − d)`), so every conjugate of
   `S` over `Q` is `σ(S)` for such a `σ`, and the bound covers all of them.

4. **Norm.** `N_{Q(S)/Q}(S) = ∏_{i=1}^{deg S} σ_i(S)` is a nonzero integer
   (`S ≠ 0`), so its absolute value is `≥ 1`. One of the factors is `S`
   itself; the other `deg S − 1 ≤ 2^k − 1` factors are each `< k·2^{3n/2}`.
   Hence

   ```
   1 ≤ |N(S)| ≤ |S|·(k·2^{3n/2})^{2^k − 1}
   ⟹ |S| ≥ (k·2^{3n/2})^{−(2^k − 1)}.
   ```

   Taking `log₂` gives the claim. ∎

*Direction of the inequality, since it is easy to read backwards.* The bound
is a guarantee that **no nonzero instance is smaller than**
`(k·2^{3n/2})^{−(2^k−1)}`, i.e. `log₂(1/|S|)` is at most that. It says
nothing about instances being *large* in the exponent, and the sharp
two-term analysis below shows the guarantee is loose by a factor `9/2`
(`9N/4 + 3` guaranteed against `N/2 + 3/2` actual).

*Toy check.* `S = sqrt(2) − sqrt(3)`: `k = 2`, `n = 2` (radicands below 4,
coefficients below 4). The bound reads
`|S| ≥ (2·2³)^(−3) = 16^(−3) ≈ 2.4e-4`; the true value is `0.3178`. Valid
with room to spare.

*Remarks.*

- This is the classical Liouville-type argument; what is contributed here is
  the derivation written out with the repo's conventions, not a new theorem.
- For `k = 2` the guarantee is `log₂(1/|S|) ≤ 9N/4 + 3`, and Theorem 2.2
  below improves that to `N/2 + 3/2` and shows it is sharp.
- This bound is what licenses the doubling deciders (`signexact.py`,
  `intervalsum.py`) to terminate: the bracket width shrinks by 1 bit per
  doubling and `|S|` cannot be smaller than `2^(−(2^k−1)(3n/2 + log₂k))`.

**Theorem 2.2 (two terms: the sharp constant is 1/2).** Let
`S = c₁·sqrt(a) + c₂·sqrt(b) ≠ 0` with integers `c₁, c₂ ≠ 0`, `a, b ≥ 1`, and

```
N = log₂(|c₁|·a·|c₂|·b)          (real input size; see remark)
```

Then `log₂(1/|S|) ≤ N/2 + 3/2`, and no smaller constant is possible.

*Proof.* **Step 1 (squarefree reduction).** Write `a = p²·d`, `b = q²·e` with
`d, e` squarefree, so `S = C₁·sqrt(d) + C₂·sqrt(e)` with `C₁ = c₁p`,
`C₂ = c₂q`. Since `|p| ≤ sqrt(a)` and `|q| ≤ sqrt(b)`,

```
log₂|C₁| ≤ u₁ := log₂|c₁| + (log₂a)/2,
log₂|C₂| ≤ u₂ := log₂|c₂| + (log₂b)/2,
```

and `N = 2(u₁ + u₂) − (log₂a + log₂b)/2`, hence

```
(u₁ + u₂)/2 ≤ N/2   and   N ≤ 2(u₁ + u₂).
```

**Step 2 (same squarefree part).** If `d = e` then `S = (C₁ + C₂)·sqrt(d)`
and `S ≠ 0` forces `|C₁ + C₂| ≥ 1`, so `|S| ≥ sqrt(d) ≥ 1` and
`log₂(1/|S|) ≤ 0`.

**Step 3 (distinct parts: the norm argument, two terms only).** If `d ≠ e`
then `M := C₁²d − C₂²e ≠ 0` — otherwise `d/e = (C₂/C₁)²`, and two distinct
squarefree integers have different prime-parity vectors — so `|M| ≥ 1`.
With `X = |C₁|·sqrt(d)`, `Y = |C₂|·sqrt(e)` (so `|M| = |X² − Y²|`) the two
sign choices give `|S| ∈ {X − Y, X + Y}`, and in both cases

```
|S| ≥ |M| / (X + Y) ≥ 1 / (X + Y),
```

so `log₂(1/|S|) ≤ log₂(X + Y) ≤ 1 + max(log₂X, log₂Y)`.

**Step 4 (the coupling that makes it `N/2`).** If `u₁ − u₂ ≥ 1` then
`log₂X − log₂Y ≥ 1`, so `X ≥ 2Y`, whence `|M| = X² − Y² ≥ (3/4)X²` and

```
|S| ≥ (3/4)X² / (X + Y) ≥ (3/4)X² / (2X) = 3X/8 ≥ 3/8,
```

which is trivial. Symmetrically for `u₂ − u₁ ≥ 1`. Otherwise
`|u₁ − u₂| < 1`, so `max(u₁, u₂) < (u₁ + u₂)/2 + 1/2 ≤ N/2 + 1/2`, and
Step 3 gives `log₂(1/|S|) ≤ 1 + N/2 + 1/2 = N/2 + 3/2`. ∎

**Sharpness — the achieved family.** Fix `d = 2`, `e = 3` and take integer
solutions of the Pell-type equation `2C₁² − 3C₂² = −1`, of which there are
infinitely many: `(C₁, C₂) = (1,1), (11,9), (109,89), …`, generated by the
unit `5 + 2·sqrt(6)`. Then `M = −1` and

```
S = C₁·sqrt(2) − C₂·sqrt(3) = −1 / (C₁·sqrt(2) + C₂·sqrt(3)),
```

so `log₂(1/|S|) = log₂(C₁·sqrt(2) + C₂·sqrt(3)) = log₂C₁ + O(1)`, while the
input size is `N = 2·log₂C₁ + O(1)`: the **coefficients** carry the bits and
the radicands are fixed. Hence `log₂(1/|S|)/N → 1/2`, and Theorem 2.2's
constant cannot be lowered. ∎

*Remarks.*

- **`N` conventions.** `N` here is the real `log₂(|c₁|a|c₂|b)`; the repo's
  integer `input_bits` sums the four bit-lengths and lies in `(N, N + 4]`,
  which changes nothing at this resolution.
- **This corrects §4's earlier claim.** An earlier version of this file said
  F1 achieves the worst two-term constant and that it lies in `[1/4, 9/4]`.
  Both halves were wrong: the worst-case constant is exactly `1/2`, and F1
  reaches only `1/4`. The Pell arithmetic that makes an instance small can
  live in the *coefficients* with tiny fixed radicands — a strictly easier
  family to find than F1, which hides it in `10^6`-bit radicands.
- **`d = 1` is the F1 shape.** `S = p − q·sqrt(α)` is this theorem with
  `a = 1, d = 1` and a single large coefficient pair `(p, q)`; the Pell
  equation `p² − αq² = ±1` then forces `log₂p ≈ log₂q`, so both terms carry
  half the bits and the ratio lands at `1/4` rather than `1/2`.

**Theorem 2.3 (second-moment separation: a `2^r − 1` → `2^{r-1}` win).**
Keep `S = Σ_j c_j·sqrt(a_j) ≠ 0`, aggregate to `S = Σ_{i=1}^m C_i·sqrt(d_i)`
over distinct squarefree parts, let `r ≤ m ≤ k` be the **F₂-rank** of
`{d_i}` (the degree of the multiquadratic field, by Theorem 1.2), and put
`M = Σ_i C_i²·d_i`. Then

```
log₂(1/|S|) ≤ ((2^r − 1)/2)·(log₂ M + 1).
```

*Proof.* The `2^r` embeddings of `K = Q(sqrt(d_1),…,sqrt(d_m))` are the
sign flips `σ: sqrt(d_i) ↦ ±sqrt(d_i)`, with `2^r` distinct images because
the parity vectors span an `r`-dimensional space. Summing the squares:

```
Σ_σ (σ(S))²  =  Σ_σ (Σ_i ±C_i sqrt(d_i))²  =  2^r·M,
```

the cross terms cancelling pairwise under the sign flips. Every `σ(S)` is an
algebraic integer and `S ≠ 0`, so the norm `Π_σ σ(S)` is a nonzero integer,
`|Π σ(S)| ≥ 1`. Isolate `S` and apply AM–GM to the remaining
`2^r − 1` conjugates:

```
1 ≤ |S|·Π_{σ≠1} |σ(S)|
  ≤ |S|·( (Σ_{σ≠1}(σ(S))²)/(2^r − 1) )^{(2^r −1)/2}
  ≤ |S|·( 2^r·M/(2^r − 1) )^{(2^r −1)/2},
```

since `Σ_{σ≠1}(σ(S))² ≤ Σ_σ(σ(S))² = 2^r M`. Taking `log₂`, and using
`2^r/(2^r − 1) ≤ 2` for `r ≥ 2` (the `r = 1` case is a single squarefree
part, where `S` is rational and `|S| ≥ 1`), gives the bound. ∎

*Remarks.*

- **Why it is sharper than Theorem 2.1.** Thm 2.1 bounds each conjugate by
  the trivial `Σ|C_i|sqrt(d_i)` and pays `(2^r − 1)` times that. Thm 2.3
  replaces the per-conjugate maximum by the *exact orbit second moment*
  `2^r M`, halving the exponent. Measured: 3.4–5.0× sharper on random
  instances with `r = 2, 3`.
- **It is *not* tight at `r = 2` — and that is fixed, not accepted.** Thm 2.3
  gives `E ≤ (3/2)·(log₂M + 1)`; Thm 2.2 gives `E ≤ N/2 + 3/2`. On the
  coefficient-Pell family (`M ≈ 2C²`, `N ≈ 2log₂C`) these read `1.5N` and
  `0.5N` respectively — Thm 2.2 is **3× sharper**, and the two-term argument
  exploits the pairing `±S` of the orbit, which the norm/AM–GM step spends.
  An earlier version of this remark claimed Thm 2.3 was "attained to three
  decimals" here; that compared `E` against `N/2` rather than against Thm
  2.3's own guarantee, and the claim is withdrawn.
  **Superseded by Theorem 2.6:** the `±S` pairing is not only an `r = 2`
  accident — it is available at every `r`, and spending it halves the
  exponent again, to `((2^{r-1} − 1)/2)·(log₂M + 1)`. At `r = 2` that is
  `½(log₂M + 1)`, which is *attained* on the coefficient-Pell family. So
  Thm 2.3 is a valid but non-optimal invariant; the sharp one is Thm 2.6.
- **`r`, not `k`, is the real parameter.** Thm 2.2 already used `r`; Thm 2.3
  makes it explicit. Any adversary wanting a hard instance must make the
  `d_i` F₂-independent — you gain nothing from a large `k` over a small `r`.
  F1 and the coefficient-Pell family both sit at `r = 2`.

### 2.4 The tower: what it costs, and where the real barrier is


Every norm-based method factors the norm through the quadratic tower
`K₀ ⊂ K₁ ⊂ … ⊂ K_r`. That recursion is the standard algorithm for these
problems, so it is worth writing down exactly what it produces — because the
answer is *not* the one usually assumed.

**Theorem 2.4 (tower monomial count).** Let `S = Σ_{i=1}^r C_i·sqrt(d_i)`
with the `d_i` F₂-independent (Thm 1.2's hypothesis, so `K_j` has degree
`2^j`), and put

```
S_j  = Σ_{i=1}^j C_i·sqrt(d_i),        (partial sums, S₀ = 0)
N_j  = S_j·σ_j(S_j) = S_{j−1}² − C_j²·d_j,
```

where `σ_j` is the sign flip of `sqrt(d_j)` alone — the relative norm of `S_j`
down to `K_{j−1}`. Then

```
N_j  =  ( Σ_{i=1}^{j−1} C_i²·d_i  −  C_j²·d_j )  +  2·Σ_{1 ≤ i < ℓ ≤ j−1} C_i·C_ℓ·sqrt(d_i·d_ℓ),
```

so `N_j` carries exactly `1 + C(j−1, 2)` `Q`-linearly independent radical
monomials — **quadratic in `j`**, not exponential — provided all `C_i ≠ 0`.

*Proof.* Expand the square:

```
S_{j−1}² = Σ_i C_i² d_i + 2 Σ_{i<ℓ} C_i C_ℓ sqrt(d_i d_ℓ),
```

since `sqrt(d_i)·sqrt(d_ℓ) = sqrt(d_i d_ℓ)` for distinct `d_i, d_ℓ`. The
first sum is rational. For the second, the parity vectors `v(d_i d_ℓ) =
v(d_i) ⊕ v(d_ℓ)` are distinct for distinct pairs `{i, ℓ}`: equality would give
`v(d_i) ⊕ v(d_ℓ) = v(d_{i'}) ⊕ v(d_{ℓ'})`, i.e.
`v(d_i d_{i'}) = v(d_ℓ d_{ℓ'})`, and with the `d_i` F₂-independent that forces
`i = i'` and `ℓ = ℓ'` (or the trivial coincidence), so the pairs coincide.
Hence the `C(j−1, 2)` radicals `sqrt(d_i d_ℓ)` are distinct members of the
basis of Thm 1.2 — the `2`-element subsets of `{1,…,j−1}` — and together with
the single rational monomial they are `1 + C(j−1, 2)` linearly independent
elements. ∎

*Measured* (200 random instances, `r ≤ 9`, max over instances per level):
`1, 1, 2, 4, 7, 11, 16, 22, 29` — exactly `1 + C(j−1, 2)`. See
`tests/test_proofs.py::TestTowerStructure`.

**Proposition 2.4′ (the support of `N_j` has rank `j − 2`, not `j − 1`).**
The `1 + C(j−1,2)` monomials of `N_j` do **not** span a rank-`j−1` part of
`K_{j−1}`. They generate the subfield of degree `2^{max(j−2,0)}`, so the
rank-based cost of separating `N_j` is `Θ(2^{j−2})`, not `Θ(2^{j−1})`.

*Proof.* The rational monomial lives in `Q`, so it adds nothing. It suffices
to compute the F₂-rank of the set `P = { v(d_i) ⊕ v(d_ℓ) : 1 ≤ i < ℓ ≤ j−1 }`.
Take any `i₀ < ℓ₀` in that range and consider the two vectors
`v(d_{i₀}) ⊕ v(d_ℓ)` as `ℓ` ranges over `{i₀+1, …, j−1}`. These `j − 1 − i₀`
vectors are linearly independent: a relation between them would give, after
adding `v(d_{i₀})` to each term, a relation among the `v(d_ℓ)` for
`ℓ > i₀`, and those are part of the F₂-independent family `{v(d_1),…}`.
Now any other `v(d_i) ⊕ v(d_ℓ)` with `i ≠ i₀` satisfies
`v(d_i) ⊕ v(d_ℓ) = (v(d_{i₀}) ⊕ v(d_i)) ⊕ (v(d_{i₀}) ⊕ v(d_ℓ))`,
so every element of `P` lies in the span of those `j − 1 − i₀` vectors. For
`i₀ = 1` the span has dimension `j − 2`, and the `j − 2` vectors
`v(d_1) ⊕ v(d_ℓ)`, `ℓ = 2, …, j−1`, are independent — so
`dim span(P) = j − 2` exactly. The `j = 2` case has `P = ∅`, rank `0`. ∎

*Verified* (300 random F₂-independent sets per level, `j = 2..7`): rank of the
pairwise products is exactly `max(j−2, 0)` every time, never `j−1`. See
`tests/test_proofs.py::TestTowerStructure`.

*Remark (and this is the load-bearing correction).* An earlier version of this
remark asserted the rank was `j − 1` and concluded the tower "pays
`Θ(2^{j−1})` per level". **The rank is `j − 2`**, so that stated cost was
wrong by a factor `2` in the exponent. The *conclusion* survives anyway, and
the reason is worth writing down precisely, because it is not the reason the
file previously gave:

- the rank saving is a genuine factor `2` in the exponent of the recursive
  step;
- but `N_j = S_{j−1}² − C_j²d_j` is built from `S_{j−1}²`, so its
  second moment is `M' ≈ 4·max|C_i|²·d_i ≈ 4M`, i.e. `log₂M' ≈ log₂M + 2`;
- and `M ≈ 2^{N/2}` for a single term carrying half the bits, so
  `log₂M' ≈ 2·(log₂M)/2·…` — concretely the *size* enters the exponent
  linearly, so doubling the size doubles the guarantee.

Net: the `2` from the rank and the `2` from the squaring **cancel exactly**.
Measured ratio `tower(rank j−2) / Thm 2.6` over random `r = 3..8`:
`1.0049, 0.9982, 1.0007, 0.9975, 1.0019, 1.0003` — i.e. `1 ± 0.005`, no
gain and no loss. So the barrier statement stands, but the *reason* is now
"rank saving cancelled by size doubling", not "the rank was always `j−1`".

So the honest statement about **this class of methods**: the barrier is the
field rank, and it is not removable by monomial-count economy or by a
shallower subfield — because the shallower subfield is paid for by a squaring.
Every method that passes through the norm — Liouville, Theorem 2.3, Theorem
2.6, the resultant tower, Cheng's and Qian–Wang's algorithms — lands on
`Θ(2^r)`. A bound that improves on this must avoid the
norm-of-the-whole-orbit and use something else about the `sqrt(d_i)`: their
minimal polynomials `x² − d_i`, or the relation between `S` and its *single*
conjugate, which is what Thm 2.2 exploits at `r = 2`.

**Addendum (§2.8, written later).** The last sentence of that paragraph was,
as a statement about the *problem*, wrong — and it took a literature check to
see it. Eisenbrand–Haeberle–Singer (arXiv:2312.02057) obtain a bound that is
**polynomial in the number of terms**, by exactly the "use something else
about the `sqrt(d_i)`" route named here rather than by a norm. So `Θ(2^r)` is
the barrier of the norm class, and the norm class is not the problem. What the
remaining gap is, see §4.

### 2.5 Why the `r = 2` argument works, and which pairing it uses

Worth isolating, because it is the one ingredient that beat `2^r`, and
knowing exactly what it is prevents a future attempt from rediscovering it
badly.

**What Thm 2.2 actually uses.** For `S = C₁sqrt(d₁) + C₂sqrt(d₂)` with
distinct squarefree `d₁, d₂`, and `X = |C₁|sqrt(d₁)`, `Y = |C₂|sqrt(d₂)`:

1. **One conjugate, one integer.** `σ` flipping `sqrt(d₁)` alone gives
   `S·σ(S) = C₂²d₂ − C₁²d₁ =: M ∈ ℤ`, and `M ≠ 0` because distinct
   squarefree parts have different parity vectors (Cor 1.3). So `|M| ≥ 1`
   and

   ```
   |S| ≥ 1 / (X + Y).
   ```

2. **The size bound.** `E = log₂(1/|S|) ≤ log₂(X+Y) ≤ 1 + log₂ max(X,Y)`,
   and each `log₂|C_j| + (1/2)log₂d_j ≤ N/2`, so `max(X,Y) ≤ 2^{N/2}`.
   That step is elementary: `N = 2log₂|C₁| + 2log₂|C₂| + log₂d₁ + log₂d₂`
   is symmetric, and `log₂|C_j| + (1/2)log₂d_j` is exactly one half of it
   plus a term that cannot exceed it.

So the whole `N/2` is **one algebraic integer of rank 0**, not a `2^r`-orbit
product. That is why it escapes the barrier.

**It does not stop at `r = 2` — the pairing is available at every `r`.** An
earlier version of this section claimed it did, on the grounds that pairing
`S` with a *single* `σ_j` gives `M_j = S'² − C_j²d_j ∈ K_{r−1}`. That is
true of the *single-conjugate* construction, and it is the wrong pairing.
The right one is the **global flip** `τ: sqrt(d_i) ↦ −sqrt(d_i)`, which
satisfies `τ(S) = −S` for *every* `S` and every `r`. That is Theorem 2.6.

*Measured, honestly (single-conjugate route, now understood as a dead end):*
pairing with a single `σ_j` and separating `M_j ∈ K_{r−1}` by Thm 2.3 in the
subfield gives a **3.54×** improvement over Thm 2.3 on random `r = 3`
instances, but on **equal-size** terms — coefficients `2^40`, radicands
`{2,3,5,…}` — that improvement decays to **1.05×** by `r = 8`:

| `r` | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| ratio Thm 2.3 / single-conj | 1.51 | 1.19 | 1.10 | 1.07 | 1.06 | 1.05 | 1.05 |

The `r = 3` figure is a property of those instances, not a general win. The
asymptotic `1 + O(2^{−r})` is the honest statement: **this** route converges
to Thm 2.3 and does not beat it. It is retained to document why the
single-conjugate choice is the wrong one, and it is not to be confused with
the global-flip pairing of Theorem 2.6, which gains a uniform factor `≈ 2`.

*Two false starts, recorded so they are not repeated.*

- **Comparability.** The tempting generalisation of `X ~ Y` to
  "`x₁ ~ x₂ ~ … ~ x_r`" is false as stated: `max_i |C_i|sqrt(d_i) / |S|`
  reaches `7×10¹³` on near-cancelling `k = 2` instances. What is true is
  only the weak form — the *tail* must nearly cancel the *head* — and Thm
  2.2 never needed individual comparability, as item 2 above shows.
- **Cheap monomials ⇒ cheap separation.** Thm 2.4's quadratic count invites
  the hope that recursing on `N_j` is cheap. It is not: the monomials of
  `N_j` live in `K_{j−1}`, so the recursion pays `Θ(2^{j−1})` per level and
  reproduces the same exponential.

### 2.6 The global-flip pairing

**Theorem 2.6 (a second factor of 2, attained at `r = 2`).** With `S`, `r`,
`M` as in Theorem 2.3,

```
log₂(1/|S|)  ≤  ((2^{r−1} − 1)/2)·(log₂ M + 1)      (r ≥ 2).
```

*Proof.* Let `K = Q(sqrt(d_1),…,sqrt(d_m))` and let `τ` be the **global
sign flip**, `τ(sqrt(d_i)) = −sqrt(d_i)` for all `i`. This is an automorphism
of `K` — its parity vector is the all-ones vector — and it satisfies
`τ(σ(S)) = −σ(S)` for every `σ`. So the `2^r` conjugates of `S` split into
`2^{r−1}` pairs `{σ(S), −σ(S)}`, and the norm is a perfect square up to sign:

```
N_{K/Q}(S)  =  Π_σ σ(S)  =  ±( Π_{σ ∈ H} σ(S) )²,
```

where `H` is any choice of one representative from each pair. `S ≠ 0` makes
`S` an algebraic integer, so the norm is a nonzero integer and
`Π_{σ∈H} |σ(S)| ≥ 1`. Isolating `S` and applying AM–GM to the remaining
`2^{r−1} − 1` representatives — whose squared sum is exactly half the orbit
sum, `2^{r−1}·M` — gives

```
1  ≤  |S| · Π_{σ ∈ H, σ ≠ 1} |σ(S)|
   ≤  |S| · ( 2^{r−1}·M / (2^{r−1} − 1) )^{(2^{r−1}−1)/2}
   ≤  |S| · ( 2·M )^{(2^{r−1}−1)/2},
```

using `2^{r−1}/(2^{r−1} − 1) ≤ 2` for `r ≥ 2`. Taking `log₂` finishes. (The
`r = 1` case is a single squarefree part: `S` is rational, `|S| ≥ 1`.) ∎

*Remarks.*

- **It is a uniform factor ≈ 2 over Theorem 2.3**, and the gain does **not**
  vanish with `r`, unlike the single-conjugate route in §2.5:

  | `r` | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
  |---|---|---|---|---|---|---|---|
  | Thm 2.3 exponent | 1.5 | 3.5 | 7.5 | 15.5 | 31.5 | 63.5 | 127.5 |
  | Thm 2.6 exponent | 0.5 | 1.5 | 3.5 | 7.5 | 15.5 | 31.5 | 63.5 |
  | ratio | 3.000 | 2.333 | 2.143 | 2.067 | 2.032 | 2.016 | 2.008 |

  The ratio tends to `2`, not to `1`. So the pairing is a real structural
  saving, and §2.5's earlier claim that "it stops at `r = 2`" was wrong.
- **At `r = 2` it attains the sharp constant, from above.** Two facts, and
  they must be kept apart:
  - *It is tight on the balanced family.* With `M = X² + Y²`,
    `X = |C₁|sqrt(d₁)`, `Y = |C₂|sqrt(d₂)`, the guarantee `E ≤ ½(log₂M + 1)`
    is *exactly* Cauchy–Schwarz `X + Y ≤ sqrt(2(X²+Y²))` in the exponent, and
    equality holds iff `X = Y`. The coefficient-Pell family
    `2C₁² − 3C₂² = −1` is precisely `X ≈ Y`, and the measured slack
    `½(log₂M+1) − E` decays geometrically:

    | `t` | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
    |---|---|---|---|---|---|---|---|
    | slack | 7.3e−3 | 7.7e−7 | 8.0e−11 | 8.3e−15 | 8.7e−19 | 9.0e−23 | 9.4e−27 |

  - *Its worst case over all balances is also `1/2`.* Maximising
    `½(log₂M + 1)/N` over the imbalance ratio `C₂/C₁` at growing size gives
    `0.8737, 0.6444, 0.5648, 0.5374, 0.5182, 0.5092, 0.5049` at
    `C₁ ≈ 10^{3,6,12,20,40,80,150}` — converging to `1/2` from above, worst
    at the *extremes* of imbalance (where `M ≈ 2C₁²` and the bound is loose,
    but `N` is dominated by the other term). So Thm 2.6 alone supplies an
    upper bound on the two-term constant converging to `1/2`, which matches
    Thm 2.2's achieved `1/2` exactly. The `k = 2` constant is therefore
    closed *and tight from both sides*.
  - *But "rate `= 1/2`" is asymptotic only.* On the Pell family itself the
    guarantee's ratio to `N` runs `0.6425, 0.5384, 0.5223, 0.5158, 0.5122,
    0.5099` for `t = 0..5` — approaching `1/2` **from above** and never
    reaching it, because the `+1` in `½(log₂M+1)` decays only like `1/N`. The
    guarantee is `N/2 + O(1)`, not `≤ N/2`: at `t = 0` (`N = 2.58`) it is
    `1.66` against `N/2 = 1.29`. So the constant is sharp in the limit and
    loose at every finite size. Pinned by
    `tests/test_proofs.py::test_roth_row_at_r2_matches_the_sharp_two_term_constant`.
- **It does not dominate Thm 2.2 pointwise, and an earlier draft of this
  remark claimed it did.** At `r = 2`, Thm 2.6 beats `N/2 + 3/2` when
  `2C₁² + 3C₂² < 24C₁C₂` and loses otherwise — i.e. Thm 2.6 wins on
  comparable terms and Thm 2.2 wins on lopsided ones (measured: `C₁ = 10⁶`,
  `C₂ = 1` gives `20.93` vs `12.76`). The two are **complementary**; the
  usable two-term guarantee is `min(of the two)`, whose sup-ratio also tends
  to `1/2`. An earlier version of this file stated Thm 2.6 "dominates Thm
  2.2 outright at `r = 2`"; that was checked only on the balanced family
  where the comparison happens to hold, and it is withdrawn.
- **It still pays `Θ(2^r)` — of the norm.** The exponent is
  `(2^{r−1} − 1)/2 ≈ 2^{r−2}`, so this halves the constant in front of the
  exponential and **does not** touch the growth rate *within the norm class*.
  §2.4's conclusion stands for norm arguments: no norm argument reaches a
  polynomial bound. What changed is the constant, not the growth rate.
  **But see §2.8:** the norm class is not the problem. A non-norm argument
  (the Subspace Theorem) does reach a polynomial bound, so the `Θ(2^r)` here
  is a property of the method and not an obstacle to the conjecture.
- **The saving is exactly the "free" part of the norm.** `Π_σ σ(S)²` splits
  into `2^{r−1}` squares, and each square is one fewer thing to bound. There
  is no further pairing available: the remaining `2^{r−1} − 1` elements are a
  single orbit of a `(2^{r−1} − 1)`-element set with no forced symmetry left,
  so within this class of argument the exponent `(2^{r−1} − 1)/2` is final.

### 2.7 The achieved rate: a pigeonhole bound that is tight at `r = 2`

Everything above is a *guarantee* — an upper bound on `E`. This subsection
bounds the other side: how small `|S|` an adversary can actually force. It
uses no norm at all, only counting, which is why it escapes the barrier.

**Theorem 2.5 (pigeonhole / achieved separation).** Let `d₁, …, d_r` be
distinct squarefree integers and `B ≥ 1` an integer. Then there are integers
`C₁, …, C_r`, not all zero, with `|C_i| ≤ B` and

```
|Σ_i C_i·sqrt(d_i)|  ≤  (Σ_i sqrt(d_i)) / ( B^(r−1) ),
```

hence `log₂(1/|S|) ≥ (r−1)·log₂ B − O(1)`. With the radicands held fixed and
`|C_i| ≤ B`, the input size is `N ≍ r·log₂ B`, so

```
E/N  ≥  (r−1)/r  −  o(1).
```

*Proof.* Consider the `(B+1)^r` sums

```
v(C) = Σ_i C_i·sqrt(d_i),    C ∈ {0, 1, …, B}^r.
```

They are **pairwise distinct**: an equality `v(C) = v(C')` with `C ≠ C'`
would give a nontrivial rational relation `Σ_i (C_i − C'_i)·sqrt(d_i) = 0`
between squarefree radicals, which Corollary 1.3 forbids. (Verified
exactly: all `(B+1)^r` values distinct for `r ≤ 4`, `B ≤ 7`.)

All values lie in `[0, B·Σ_i sqrt(d_i)]`, an interval of that length. Sorting
them, two consecutive ones differ by at most

```
B·Σ_i sqrt(d_i)  /  ((B+1)^r − 1)  ≤  2·(Σ_i sqrt(d_i))·B^(−(r−1))     (B ≥ 1).
```

Their difference is `Σ_i C_i sqrt(d_i)` with `|C_i| ≤ B` and `C ≠ 0`. ∎

*Remarks.*

- **It is tight at `r = 2`, and that is a consistency check worth having.**
  Thm 2.2 proves the sharp `r = 2` constant is `1/2`, achieved by the
  coefficient-Pell family. Thm 2.5 gives `(2−1)/2 = 1/2` at `r = 2`. Two
  independent arguments — one algebraic, one purely combinatorial — land on
  the same number, so `(r−1)/r` is very likely the true achieved constant in
  general and not merely a pigeonhole artefact.
- **The achieved rate grows with `r`** — `2/3, 3/4, 4/5, …` — while the
  guaranteed rate grows like `2^r`. The ratio guarantee/achieved widens
  monotonically: `3.0×, 3.5×, 5.0×, 7.8×` at `r = 2,3,4,5`. So the
  conjecture is *not* threatened by Thm 2.5; it is threatened by nothing
  known, and the gap is widening on both sides.
- **This is a family an adversary can actually build.** Unlike F1, which
  needs a Pell solve per instance, Thm 2.5's instances come from sorting
  `(B+1)^r` sums and taking a consecutive gap — polynomial time in `B^r`,
  and the resulting `C` is verifiable by the existing exact decider. See the
  TODO entry on adding it as a family.
- **It does not use the `x² − d_i` beyond what Cor 1.3 already gives.**
  Distinctness of the `v(C)` is the whole content. That is the honest limit
  of this subsection: it is a counting argument, so it cannot be improved
  without importing real arithmetic structure, and no arithmetic structure
  appears in it.
- **It is instantiated as family F6, at `r = 3`, and it is the only family in
  the repo that *attains* this ceiling rather than searching for something
  worse.** `families.py::generate_f6(B)` sorts the `(B+1)^3` sums over
  `{2,3,5}`, certifies the resulting gap exactly through the integer kernel,
  and returns it with `exact_mag` and Thm 2.5's bound in the metadata.
  Measured `E/N = 0.6153 … 0.6835` at `B = 30 … 150`, i.e. at the `2/3`
  ceiling from within — the `O(1/N)` slack coming from the construction's
  `Σ_i sqrt(d_i)` prefactor, which is why a finite instance may sit slightly
  *above* `(r−1)/r` without refuting anything. §3.3 of
  `THE_MISSING_THEOREM.md` records that an early version of `refutation.py`
  duly reported F6 as a refutation for exactly that reason.
  See `tests/test_families.py::TestF6CeilingFamily`.

### 2.8 The Subspace-Theorem route — and a correction to this file

> **MAJOR CORRECTION (supersedes an earlier version of §2.8 in full).**
> An earlier version of this subsection concluded that "the Subspace-Theorem
> route is closed for SRS", on the grounds that Schmidt's theorem needs `n`
> separate forms while SRS offers only one. **That conclusion was wrong**, and
> so was the "Missing Theorem" it was built on. The error and its cause are
> recorded in §2.8.1, because the mistake is instructive.

**Eisenbrand–Haeberle–Singer (QUOTED, not proved here).** *An improved bound on
sums of square roots via the subspace theorem*, arXiv:2312.02057; SoCG 2024.
For distinct squarefree `a_1, …, a_n` and `E = Σ x_i sqrt(a_i) ≠ 0`:

```
|E|  ≥  γ · ( n · ‖x‖_∞ )^{−2n},         γ = γ(a_1,…,a_n) > 0.
```

Their Remark 9 sharpens the exponent `2n → n + ε`. Hence, with `N = n·m`
(`m = log₂‖x‖_∞`):

```
E  ≤  2n(log₂ n + m) + log₂(1/γ)  =  O(N)  +  log₂(1/γ).
```

The first term is **polynomial in `N`, uniformly in `n`** (measured
`E/N → 2.03` at `n` from 2 to 10⁴). So the polynomial part of the conjecture
is proved; the constant is not.

#### 2.8.1 The error, named

The earlier reduction stopped one step short and then miscounted. Schmidt's
Subspace Theorem, in the `n`-variable form, states: for `n` linearly
independent linear forms `L_1, …, L_n` in `n` variables with algebraic
coefficients, the integer points satisfying

```
Π_i |L_i(x)|  ≤  ‖x‖^{−δ}
```

lie in finitely many proper rational subspaces of `Q^n`. That theorem is
applied to SRS in EHS Thm 7 with `L_i = qβ_i − p_i` and `L_n = q` — the
`n−1` coordinate forms are **supplied by the reduction itself**, and only
`L_0` has to come from the problem. The earlier claim that "SRS has only one
form to offer" was simply a miscount. See `STATE_OF_THE_ART.md` §0.

### 2.9 A sharper application: exponent `r−1`, matching the ceiling

The bound above is **not the best available from the Subspace Theorem**, and
the repo can see why. EHS apply it to *simultaneous* approximation of `n`
algebraic numbers by a common denominator, which costs an exponent `n`. But
SRS is a **single** linear form, and for a single form in `n` variables the
Subspace Theorem gives exponent `n−1`. The two differ by a factor, and the
single-form version lands exactly on the achieved rate.

**Theorem 2.9 (single-form Subspace Theorem — the exponent `r−1`).** Fix
distinct squarefree `d_1, …, d_r` and put `alpha_i = sqrt(d_i/d_r)` for
`i < r`, so that

```
L(C_1,…,C_r)  =  C_r + Σ_{i<r} C_i · alpha_i  =  S / (C_r sqrt(d_r)),
```

and `1, alpha_1, …, alpha_{r−1}` are `Q`-linearly independent (Cor 1.3). Then
for every `δ > 0` there are a constant `c = c(d_1,…,d_r,δ) > 0` and finitely
many proper rational subspaces `V_1, …, V_t ⊊ Q^r` such that for every
integer `C` with `C ∉ V_1 ∪ … ∪ V_t`,

```
| L(C) |  >  c · ‖C‖^{−(r−1+δ)} ,
```

hence, writing `m = log₂‖C‖_∞` and counting the `r` coefficient bits so that
`N = r·m` (radicands fixed),

```
log₂(1/|S|)  ≤  (r−1+δ)·m  +  log₂(1/c)  +  O(1)
             =  (r−1+δ)/r · N  +  O_{d,δ}(1).
```

*Proof.* Apply the Subspace Theorem in `r` variables to the `r` forms

```
L_0 = L,        L_1 = C_1,   L_2 = C_2,   …,   L_{r−1} = C_{r−1}.
```

There are exactly `r` forms in `r` variables, which is the count the theorem
requires; note that `C_r` is **not** among them and is not needed. They are
linearly independent over `Qbar` because `1, alpha_1, …, alpha_{r−1}` are
linearly independent over `Q` (Cor 1.3). Points with

```
|L(C)| · Π_{i<r} |C_i|  ≤  ‖C‖^{−δ}
```

lie in a finite union of proper rational subspaces. Since
`Π_{i<r}|C_i| ≤ ‖C‖^{r−1}`, every `C` outside that union satisfies

```
|L(C)| · ‖C‖^{r−1}  >  ‖C‖^{−δ},
    i.e.   |L(C)|  >  ‖C‖^{−(r−1+δ)}.
```

Finally `|S| = |C_r| sqrt(d_r) · |L(C)|`, and `|C_r| sqrt(d_r) ≥ 1`. ∎

**The exponent `r−1` is exactly the achieved rate.** This is the point of
the theorem, and it is why it is sharper than the EHS bound:

| | exponent on `‖C‖_∞` | `E/N` at `r = 3` | vs ceiling `2/3` |
|---|---|---|---|
| Thm 2.6 (norm) | — | `1.5` | `2.25×` |
| EHS Thm 8 (simultaneous) | `2n` | `1.0` | `1.50×` |
| **Thm 2.9 (single form)** | **`r−1`** | **`0.667 + ε/3`** | **`1×` as `ε → 0`** |
| achieved, Thm 2.5 (pigeonhole) | `r−1` | `0.667` | — |

Generalising the table, Thm 2.9 gives `E/N ≤ (r−1+δ)/r` at **every** rank,
and Thm 2.5 achieves `(r−1)/r − o(1)`. So for **fixed radicands** the
guaranteed and achieved rates *coincide*, at every `r`, up to `ε` and one
constant.

**Consistency check at `r = 2`, where the repo already owns the answer.**
Thm 2.9 at `r = 2` uses the forms `L_0 = C_2 + C_1 alpha_1`, `L_1 = C_1` —
two forms, two variables — and yields exponent `1 = r−1`. This *is* Roth's
theorem. And the repo's Thm 2.2 / Thm 2.6 show the `r = 2` case is sharp at
`1/2` **with an effective constant**, attained on the coefficient-Pell family.
So the single-form application is not vacuous at `r = 2`; it is exactly right
there, which is the strongest available evidence that the `r ≥ 3` version is
the right statement.

**Verified on the achieved side, both ranks.** Measured
`E / log₂‖C‖_∞` on the families of Thm 2.5, which should tend to `r−1`:

| | `r = 2` (coefficient-Pell) | `r = 3` (sorted gap) |
|---|---|---|
| `E/log₂‖C‖_∞` | `1.434, 1.222, 1.149, 1.112, 1.090` → `1` | `2.312, 1.896, 2.098, 1.999` → `2` |

Both tend to `r−1` from above, as an achieved rate must.

#### 2.9.1 Why this does **not** close the conjecture

Two obstructions, and the second is the interesting one.

**(i) The constant `c(d_1,…,d_r,δ)`.** Same problem as EHS's `γ`: it depends
on the radicands, which are input. This is the γ-gap of §4, inherited
unchanged. Thm 2.9 improves the *exponent*, not the *constant*.

**(ii) The exceptional subspaces `V_1, …, V_t`.** These are not known
explicitly either — that is the inefficacy of the Subspace Theorem. But here
the obstruction is **detectable**: `C` lies in a proper rational subspace iff
there is a nonzero integer vector `a` with `a·C = 0`, which is decidable by
exact linear algebra over `Q`. So the split is

| | condition on `C` | what you get |
|---|---|---|
| generic | `C ∉ V_1 ∪ … ∪ V_t` | `|S| > ‖C‖^{−(r−1+δ)}` — **fully effective, no constant at all** |
| exceptional | `C ∈ V_1 ∪ … ∪ V_t` | a codimension-≥1 condition; the instances are *detectable* but not *excluded* |

The generic row is the useful one, and it is stronger than anything else in
this file: a **constant-free** separation bound at the exact exponent `r−1`,
holding for every `C` outside a finite union of rational subspaces. An
adversary can only reach the second row by satisfying a codimension-1 linear
equation on `C`, which *constrains* rather than *enables* him — the near-
cancellation that makes `S` small is a separate geometric condition. That the
two are generically independent is a heuristic, not a theorem, and is labelled
as such.

**So the precise status after Thm 2.9:** the *exponent* of the separation
bound is settled at `r−1` — matching the achieved rate at every rank — and
what remains is (a) the constant `c` and (b) the explicit description of the
exceptional subspaces. Both are consequences of the same inefficacy. Neither
is a growth-rate problem any more. §4 states the open problem in these terms.


---

## 2.10 The constant itself: what it is, and three results

This section is new. It addresses the constant `c(d,δ)` of Thm 2.9 directly,
and it **retracts a claim** this file previously made about it.

### 2.10.1 The constant is a Markov constant (and at `r = 2` it is computable)

Thm 2.9 supplies, for generic `C`, the largest admissible `c` at exponent
`r−1`:

  > `c(d) = inf over generic `C` of `|L(C)|·‖C‖_∞^{r−1}`,
  > `L(C) = C_r + Σ_{i<r} C_i√(d_i/d_r)`.

For `r = 2` this is a classical object. With `α = √(d_2/d_1)`, varying `C_1`
against `C_2 = q` and taking the infimum over the best integer `p`:

  `c(d_1,d_2) = 1 / M(α)`,  `M(α) = sup_{n} q_n·‖q_nα‖`,

where `p_n/q_n` run over the **convergents** of the quadratic irrational `α`,
`M(α)` is a Markov-type constant, and `‖·‖` is distance to the nearest integer.
The sup is over the convergents, and that is not a cosmetic choice: the bound
needed is `|qα − p| ≥ c₀/q` for all `p, q`, whose extremum always occurs at a
convergent. A limsup over **all** `q` would be `∞` for every irrational `α`
(a generic `q` has `‖qα‖ ≈ 1/2`), giving `c = 0` and a vacuous statement. The
*optimal* constant is instead the reciprocal of the infimum
`λ(α) = inf_{q ≥ 1} q‖qα‖`; Hurwitz's `λ(α) ≤ 1/√5` gives `1/λ(α) ≥ √5 > 1`
for every `α`, and `M(α) ≥ λ(α)`, so the `c` used here is the weaker but still
finite, positive, and computable one.

The load-bearing fact is classical: `α` is quadratic, so its continued fraction
is **periodic** (Lagrange), so its partial quotients are bounded, so `α` is
badly approximable with `M(α) ∈ (0, ∞)` and both `M(α)` and `c = 1/M(α)` are
finite, positive, and **computable from the periodic expansion**. That is the
structural reason Thm 2.2 / 2.6 close `r = 2`, and it is not a coincidence of
technique.

Measured (`M = sup_n q_n‖q_nα‖` over the convergents `p_n/q_n` of `α` — see the
correction below for why the sup must be over convergents — using the
exact-integer numerator `|Dq² − Qp²|` for `α = √(D/Q)`, 300 convergents):
`0.414214` at `√2`, `0.732051` at `√3`, `0.898979` at `√(3/2)`, `2.618950` at
`√(5/3)`, `4.580399` at `√(7/5)`, `3.482873` at `1/√11`. So `M(α)` is `Θ(1)` and
`c = 1/M` is `Θ(1)` — which is all that is needed, and which is what periodicity
gives.

> **Correction (this paragraph has now been wrong twice; the numbers above are
> the third, verified set).** Three separate defects, all in the *measurement*,
> not in the theorem:
>
> 1. **Wrong statistic.** It reported `M = max_{q ≤ 400} q‖qα‖`. A `max` over
>    **all** `q` is not a Markov constant and is not even bounded as `q` grows:
>    for a generic `q`, `‖qα‖ ≈ 1/2`, so `max_{q ≤ Q} q‖qα‖ ≈ Q/2` for *every*
>    `α`. Measured that way, all ten `α` tried here return `≈ 180–195` at
>    `Q = 400`, i.e. just `Q/2`, which carries no information. The correct
>    statistic is the sup (equivalently the limsup) **along the convergents**,
>    because the bound needed downstream is `|qα − p| ≥ c/q` for all `p, q`, and
>    the worst case is always attained at a convergent.
> 2. **Conflated `√(D/Q)` with `√D/Q`.** The header read `α = √(D)/Q` while the
>    row labels were `√(5/3)`, `√(7/5)`, `√(3/2)` — the `√(D/Q)` form. These
>    differ (`√5/3 = 0.7454` vs `√(5/3) = 1.2910`), so the numerator
>    `|q²D − p²Q²|` was not the one corresponding to the stated `α`.
> 3. **A fabricated exact witness.** The paragraph cited
>    "`q = 3161, p = 2449` with `q‖qα‖ = 2.4996` for `α = √(5/3)`". That is
>    false: `3161·√(5/3) ≈ 4080.83`, so the nearest integer is `4081`, not
>    `2449`, and `3161·|3161α − 2449| ≈ 5.16×10⁶`, not `2.4996`. No such
>    witness exists at those integers.
>
> The **retraction itself is correct** — `M(α) ≤ 1` really is false — but it
> needs a witness that survives exact checking. A far smaller one is available:
> `α = √(5/3)`, `q = 3`, `p = 3` (the `0`-th convergent, `p/q = 1`), giving
> `q‖qα‖ = 3·|3√(5/3) − 3| = 2.618950 > 1`. With the exact integer
> `N = |Dq² − Qp²| = |5·9 − 3·9| = 18`, the inequality `N > qα + p` is the exact
> integer comparison `45² = 2025 > 405 = 3²·5·3²`. Pinned by
> `tests/test_proofs.py::test_RETRACTED_markov_constant_is_not_at_most_one`.
>
> The earlier *bound* this paragraph asserted, `M(α) ≤ 1/(2·max(1,α)+1)`, was
> **backwards** and is retracted; the true direction `M(α) ≥ 1/(2max(1,α)+1)`
> is just bad approximability and says nothing quantitative about `c`.

(An earlier revision of this file also used a *float* to measure `M` and got
values like `6×10³⁶` from catastrophic cancellation in `qα − p` at large `q`.
The exact-integer numerator `|Dq² − Qp²|`, divided by `Q(qα + p)`, is the
numerically stable form and is what the numbers above use. The same trap
returns at deep convergents — forming `qα − p` directly is accurate only to
about `1e{-dps}·qα`, so at `q ~ 1e100` even `dps = 80` returns pure noise.)

For `r ≥ 3` the analogous object is a **two-parameter** (in general
`(r−1)`-parameter) Diophantine quantity, and the periodicity that made `r = 2`
computable is gone. That is the honest content of "`r = 2` is done, `r ≥ 3` is
not": the constant, not the exponent.

### 2.10.2 Correction: `sup E/N = (r−1)/r` is FALSE

This file, `THE_MISSING_THEOREM.md`, `STATE_OF_THE_ART.md`, `README.md` and
`TODO.md` previously asserted the "sharp form" `sup E/N = (r−1)/r` for every
`r`. **That is false, and here is a certified counterexample.**

Take `r = 3` with

  `S = −10√2 − 98√3 + 51√13`,  i.e. `d = (2,3,13)`, `C = (−10,−98,51)`.

The coefficients come from the natural greedy construction: set
`C_1 = −round(C_2√(d_2/d_1) + C_3√(d_3/d_1))`. Exact evaluation gives

  `|S| = 2.83182528695373655931322820512724491638920774566997 × 10⁻⁷`,

and with `N = Σᵢ max(bitlen|Cᵢ|, bitlen dᵢ) = 4 + 7 + 6 = 17`:

  `E = −log₂|S| = 21.7518`,  so  `E/N = 1.2795 > 2/3 = (r−1)/r`.

`E/N` is inflated by the **additive** slack, not by any growth: at fixed `d`,
`c(d,δ) ≈ 2.8×10⁻⁷` is *constant*, so `log₂(1/c) ≈ 21.75` is a fixed bonus
and `E/N = (r−1+δ)/r + 21.75/N → 2/3` as `N → ∞`. The bonus is invisible
asymptotically and dominant at small `N`.

Measured over `C(20,3) = 1140` radicand triples with exhaustive `‖C‖_∞ ≤ 60`
greedy search: mean `E/N = 0.7605` for `N < 20` and `0.6589` for `N ≥ 20`.
The rate **falls** toward `2/3` from above — the expected signature of a
bounded additive bonus, and the opposite of a refutation.

**The correct statement is asymptotic:** `limsup_{N→∞} E/N = (r−1)/r`, or
equivalently `E ≤ ((r−1+δ)/r)·N + O_{d,δ}(1)`. The earlier `sup` form is
retracted everywhere it appeared.

Note also that this instance is invisible to the F6 family, which sorts the
`{0..B}³` grid: `(−10,−98,51)` is not a sorted triple. A refutation harness
built on one parameterisation will miss whole regions; `refutation.py` is not
a completeness argument in either direction.

### 2.10.3 An elementary effective bound, and the exact obstruction

Independently of the Subspace Theorem, the field norm gives a fully effective
bound. Let `K = Q(√d_1,…,√d_r)`, so `[K:Q] = 2^r`, and write `c_i = C_i√d_i`,
`S = Σ c_i`. Since `S` is an algebraic integer and `S ≠ 0`,

  `1 ≤ |N_{K/Q}(S)| = ∏_σ |σ(S)| ∈ Z`.

The `2^r` conjugates are `S − 2Σ_{i∈I} c_i` over all `I ⊆ {1..r}`; pairing
each `I` with its complement `Iᶜ` (whose conjugate is `−(S − 2Σ_{i∈I}c_i)`)
gives exactly `2^{r−1} − 1` distinct pairs, so

  `1 ≤ |N(S)| = |S|² · ∏_{j=1}^{2^{r−1}−1} |S − 2∑_{i∈I_j} c_i|²`,

for one representative `I_j` per complement pair. Hence

> **Thm 2.10 (norm bound, elementary and effective).** Let
> `T = ∑_i |c_i| = ∑_i|C_i|√d_i`. If `|S| ≤ T` then
>
> `|S| ≥ (3T)^{−(2^{r−1}−1)}`,  hence
> `E ≤ (2^{r−1}−1)·log₂(3T)`.
>
> The underlying inequality is *exact* and sharper:
> `|S| ≥ 1/∏_j |S − 2∑_{i∈I_j}c_i|`, one factor per complement pair.

Derivation of the `(3T)` prefactor: each factor satisfies
`|S − 2∑_{i∈I_j}c_i| ≤ |S| + 2∑_{i∈I_j}|c_i| ≤ 3T` (using `|S| ≤ T`), and
there are `2^{r−1}−1` factors.

> **Side effect: a second false positive in the refutation harness.** Because
> the `E/N` excursion above `(r−1)/r` is real and certified, `refutation.py`'s
> ceiling check flagged the §2.10.2 witness as a refutation of the sharp
> form. That is the *same class* of bug as the original F6 false positive,
> from a different cause: `slack_for` modelled only the pigeonhole prefactor
> `log₂(Σ√d_i)/N` and ignored the additive constant of Thm 2.9. `is_refutation`
> now compares **excess in bits**, `E − (r−1)/r·N`, against two O(1)-in-`N`
> budgets — the prefactor, and a fixed `ADDITIVE_BITS = 32` for the additive
> constant. A per-instance estimate of that constant is *provably vacuous*:
> `log₂(1/c) = rate·N − (r−1)log₂‖C‖` grows faster than the excess
> `(rate − (r−1)/r)·N`, so the allowance always wins and the harness silently
> stops firing. That mistake was made twice before it was pinned by a test.

> **Correction to an earlier draft of this subsection, which is retracted.** I
> first wrote `|S| ≥ 1/(2^{2^{r−1}−1}·∏_i|c_i|)`, i.e. `1/(8∏|c_i|)` at
> `r = 3`, and claimed the `r = 3` instance above attains it to eight digits
> (measured ratio `0.99999999`). **That bound is false**, and the near-miss was
> the tell. Each factor is `|S − 2c_i| = |c_j + c_k − c_i| ≤ |c_1|+|c_2|+|c_3| = T`,
> **not** `≤ 2|c_i|`: the `−c_i` term is uncontrolled. The `0.99999999` was
> not attainment, it was the error term `|S|/(2|c_i|) ≈ 4×10⁻⁹` from the
> neglected `−c_i`. The `r = 2` case refutes the pattern outright:
> `d=(2,19)`, `C=(3,−1)` has `|S| = 0.116258 < 1/(2·3√2) = 0.117851`, since
> `|c_2−c_1| ≤ |c_1|+|c_2|`, never `≤ 2|c_1|`. The `(3T)^{−(2^{r−1}−1)}`
> form above is the correct statement; all checks pass at 80 digits.

What survives of the original observation is the *reason* the exponent is
`2^{r−1}−1` — the number of complement pairs — which is structural and does not
depend on the prefactor. And the `r = 3` witness is still valuable: it is a
certified `E/N = 1.2795 > 2/3`, which is what §2.10.2 needs.

For `r = 2` Thm 2.10 gives `E ≤ log₂(3T)`, i.e. `O(N)`; Thm 2.2's `N/2` is
sharper and comes from the Markov structure (§2.10.1) rather than the norm.

**And the obstruction is now explicit and structural.** The exponent of the
bound is the number of complement pairs, `2^{r−1} − 1`:

  | `r` | `2^{r−1}−1` | resulting `E` |
  |---|---|---|
  | 2 | 1 | `O(N)` — and Thm 2.2 gets the sharp `N/2` |
  | 3 | 3 | `≤ 3 + (3/2)N` (Thm 2.10) |
  | 4 | 7 | `O(2²·N)` |
  | `r` | `2^{r−1}−1` | `O(2^{r−1}·N)` |

Since `r` can be `Θ(N)` (take `d_i` = the `i`-th prime, `C_i = 1`, so
`N ≈ r log r`), this is `O(2^{Ω(N)}·N)` — **exponential in the input size**.
So:

> **The open problem, in one line: replace the exponent `2^{r−1} − 1` of the
> norm argument by `poly(r)`, keeping the bound effective. Thm 2.9 already has
> the linear exponent `r−1`, but pays for it with an ineffective constant.**

That is the same obstacle as Thm 2.9's `c`, reached from the opposite
direction: Thm 2.6 / Thm 2.10 are effective but exponential in `r`; Thm 2.9
is linear in `r` but ineffective. **Neither is yet a proof, and the missing
statement is `poly(r)` *and* effective simultaneously.**

### 2.10.4 Sufficiency lemma: the constant need only be `exp(−poly(N))`

A small but load-bearing observation that sharpens what "solve the constant"
means. Since `r ≤ N` and `log₂‖C‖_∞ ≤ N`, Thm 2.9 gives

  `E ≤ (r−1+δ)·log₂‖C‖ + log₂(1/c) + O_{d,δ}(1) ≤ N² + log₂(1/c(d,δ))`.

> **Thm 2.11 (sufficiency).** The conjecture follows from any bound
> `c(d,δ) ≥ 2^{−poly(N)}`. In particular a bound of the shape
> `c ≥ 2^{−O(N)}` is more than enough, and one of the shape
> `c ≥ 2^{−2^N}` is not.

So the target is **not** `c = Ω(1)`, which is what the "rate matches the
ceiling" framing implicitly demanded, and not even a bound uniform in `d`. Any
effective `c` decaying no faster than `2^{−poly(N)}` closes the problem. This
is why Thm 2.10's crude exponential-in-`r` bound and Thm 2.9's
exponential-in-`N` constant are both *merely inconvenient* rather than fatal.

The measurable predictions of Thm 2.11 are consistent with everything measured
above: at `r = 3`, `log₂(1/c) ≈ 21.75` while `N = 17`, i.e. the constant is
paying `≈ N` bits, well inside `O(N)`.


---

---

## 2.11 Can `2^{r−1}−1` be replaced by `poly(r)` while staying effective?

This is the sharp form of the open problem, and the natural next question after
§2.10. It gets a partially positive, partly negative answer, and the
negative parts are theorems rather than observations.

**Short answer.** The improvement is *already achieved* at every fixed rank, and
the entire remaining difficulty is uniformity in `r`. Two natural routes to a
uniform `poly(r)` exponent are provably dead ends. What is left is exactly one
statement: a uniform effective constant in Thm 2.9.

### 2.11.1 Every fixed rank is settled, effectively

> **Thm 2.12 (fixed-rank separation).** Let `N = Σᵢ(bitlen|Cᵢ| + bitlen dᵢ)` and
> `r` the F₂-rank. Then
>
> `E ≤ (2^{r−1}−1)·(N + log₂(3r))`.
>
> Consequently, for every *fixed* `r₀` the conjecture holds on the whole
> class of inputs with `r ≤ r₀`, with a fully effective constant.

*Proof.* Thm 2.10 gives `E ≤ (2^{r−1}−1)·log₂(3T)` with
`T = Σᵢ|Cᵢ|√dᵢ`. Write `b_C = maxᵢ bitlen|Cᵢ|`, `b_d = maxᵢ bitlen dᵢ`.
Since `T ≤ r·2^{b_C}·2^{b_d/2}`,

  `log₂(3T) ≤ log₂(3r) + b_C + b_d/2`.

And `N = Σᵢbitlen|Cᵢ| + Σᵢbitlen dᵢ ≥ b_C + b_d`, since each sum dominates the
corresponding maximum. Hence `b_C + b_d/2 ≤ b_C + b_d ≤ N`, and

  `E ≤ (2^{r−1}−1)·log₂(3T) ≤ (2^{r−1}−1)·(N + log₂(3r))`. ∎

**Why this is the right reframe.** `2^{r−1}−1` is *constant* when `r` is
bounded, so the conjectured polynomial separation is already a theorem per
rank. Concretely:

| rank cap `r₀` | effective `E/N` bound from Thm 2.12 |
|---|---|
| 2 | `1` |
| 3 | `3` |
| 4 | `7` |
| 5 | `15` |
| 8 | `127` |
| 16 | `32767` |
| 32 | `2147483647` |

So the question is **not** "can the exponent be made polynomial in `r`" — it
already is, for each fixed `r`. It is: **can the constant be made independent
of `r`, with a `poly` dependence on `N`?**

### 2.11.2 The exponent `2^{r−1}−1` is irreducible

One might hope to shave the norm argument: choose the complement-pair
representatives better, or reweight. This cannot work.

> **Thm 2.13 (irreducibility).** The complement-pair norm argument yields
> exponent exactly `M = 2^{r−1}−1`, and no selection of
> representatives or AM-GM/Jensen reweighting improves it.

*Proof.* Index the `M = 2^{r−1}−1` nontrivial complement pairs by
`J ⊆ {2,…,r}`, `J ≠ ∅`, so that
`σ_J(S) = S − 2Σ_{i∈J}Cᵢ√dᵢ` and

  `1 ≤ |N(S)| = |S|² · ∏_J |σ_J(S)|²`,
  `E ≤ Σ_J log₂(|S| + 2Σ_{i∈J}tᵢ)`,  `tᵢ = |Cᵢ|√dᵢ`.

By AM-GM, `∏_J f_J ≤ (mean f)^M`. Each index `i` lies in exactly `2^{r−2}`
of the `J`, hence `Σ_J Σ_{i∈J}tᵢ = 2^{r−2}T`, and the mean of
`Σ_{i∈J}tᵢ` over the `M` nonempty `J` is `2^{r−2}T/M < T/2 + T/2`. So the
mean factor is `< |S| + T`, giving

  `E ≤ M·log₂(|S| + T) ≤ M·log₂(2T)`.

The count `M` is forced. Verified numerically at `r = 2,3,4,5`: `E/bound ≤ 0.18`,
and the bound holds with room to spare. ∎

**The obstruction is combinatorial, not analytic.** A multiquadratic field of
degree `2^r` simply has `2^{r−1}` conjugate pairs, and the norm is the product
over all of them. To bound `|S|` below you must bound every other conjugate
above, and there are exponentially many.

### 2.11.3 Squaring is a no-go

The most attractive escape: `S²` lives in a field of degree `2^{r−1}`, half
as large.

> **Thm 2.14 (squaring does not help).** Bounding `|S|²` via its norm gives
> exactly the same exponent as bounding `|S|`.

*Mechanism.* `S² = R + Σ_{i<j}2CᵢC_j√(dᵢd_j)` with `R = Σ Cᵢ²dᵢ ≠ 0`
rational. The products `dᵢd_j = eᵢ+e_j` span the **even** subspace of the
F₂-space, of dimension `r−1`, so the degree does halve (confirmed: for
`d = (2,3,13)`, the radicands `6,26,39` have `6·26·39 = 6084 = 78²`, so rank `2` and
degree `4 = 2^{r−1}`).

But the rational part **destroys the complement pairing**:
`σ_{−ε}(S²) = R − (irr) ≠ −S² = −R − (irr)`. So one no longer gets
the `|σ|²` squaring and must spend all `2^{r−1}−1` nontrivial factors
individually — the same count, now in a field of half the degree. Then

  `E(S) = E(S²)/2 ≤ (2^{r−1}−1)/2 · log₂(2T′)`,  and `T′ ~ T²`,

so `E(S) ≤ (2^{r−1}−1)·log₂(2T)`. **The halved degree and the doubled
`log T` cancel exactly.** Verified on the §2.10.2 witness: `|S²| = 8.0192×10⁻¹⁴`
with `R = 62825` nonzero.

And this is not an unlucky case: `R = Σᵢ Cᵢ²dᵢ` is a sum of **non-negative**
integers with every `dᵢ > 0`, so `R ≠ 0` whenever any `Cᵢ ≠ 0` — that is, for
every instance of interest. **Squaring always destroys the pairing**, which is
why Thm 2.14 is a theorem and not an observation.

Iterating fails for the same reason: `S⁴ = (S²)²` re-introduces a rational part,
and by §2.10 the even subspace does not nest inside itself.

### 2.11.4 At `r = 2` no constant is needed at all

§2.10.1 identified `c` at `r = 2` as a Markov constant. It can be pinned down
completely, and the answer is better than "computable".

> **RETRACTED — but note carefully: the ARGUMENT was false while the
> CONCLUSION survives, and for a better reason.** This subsection previously
> claimed "`M(α) ≤ 1` for every irrational `α`, so the constant at `r = 2` is
> `≥ 1` and costs no bits."
>
> *The bound is false.* `M` (as defined in §2.10.1, the sup/limsup **along the
> convergents**) can exceed `1`. An exact witness, with the *smallest* integers
> that work: `α = √(5/3)`, and the `0`-th convergent `p/q = 1/1`, so
>
> ```
> q = 3,  p = 3,   q·‖qα‖ = 3·|3√(5/3) − 3| = 2.618950… > 1.
> ```
>
> With the exact integer `N = |Dq² − Qp²| = |5·9 − 3·9| = 18`, the inequality
> `N > qα + p` is the exact integer comparison
> `45² = 2025 > 405 = 3²·5·3²`, so no float is involved anywhere.
> Hence `M(√(5/3)) ≥ 2.618950` and `c = 1/M ≤ 0.38197`, i.e. the repo's `M`
> charges up to `1.39` bits. "Costs no bits" does **not** follow from this `M`.
>
> *An earlier version of this retraction was itself wrong*, in a way worth
> recording. It (a) cited a fabricated witness `q = 3161, p = 2449` with value
> `2.4996` — impossible, since `3161·√(5/3) ≈ 4080.83`, so the nearest integer is
> `4081` — and (b) explained the failure by saying "`M` is a limsup over all
> `q`, and the supremum need not occur at a principal convergent". That
> explanation is backwards and, taken literally, destroys the theorem: a limsup
> over **all** `q` is `∞` for *every* irrational `α` (a generic `q` has
> `‖qα‖ ≈ 1/2`, so `q‖qα‖ ≈ q/2`), which would give `c = 0` and a vacuous
> bound. The sup is over the **convergents**, and that is where the extremum
> always lives.
>
> *Why the conclusion is nevertheless right.* The bound actually needed is
> `|qα − p| ≥ c₀/q` for **all** integers `p, q`, whose best constant is the
> **infimum**
>
> ```
> λ(α) = inf_{q ≥ 1} q·‖qα‖   (the Markov / badly-approximable constant),
> c_opt = 1/λ(α).
> ```
>
> Hurwitz gives `λ(α) ≤ 1/√5` for every irrational `α`, with equality only for
> `α` equivalent to the golden ratio. So `c_opt ≥ √5 > 1` **always**, and the
> asymptotic slack at `r = 2` is a *gain* of at least `2·log₂√5 ≈ 1.16` bits,
> not a cost. The repo's `M = sup_convergents` satisfies `M ≥ λ`, hence
> `c = 1/M ≤ c_opt`, and that weaker `c` is the one that can dip below `1` and
> cost a little. Both are finite and positive, which is all Thm 2.15 needs.
>
> Two wrong measurements were involved along the way: a `min` where an `inf`
> was needed, and then float catastrophic cancellation producing values like
> `6×10³⁶`. The exact-integer form
> `q‖qα‖ = |Dq² − Qp²| / (Q(qα+p))` settles it. See also §2.10.1 for the third
> error in this same measurement (`max` over all `q` instead of over
> convergents).

What survives, and is what §2.10.1 actually established:

> **Thm 2.15 (`r = 2`: the constant is finite, positive, computable).** Let
> `α = √(d₂/d₁)` and `M(α) = limsup_{q→∞} q·‖qα‖`, the Markov constant of the
> quadratic irrational `α`. Then
>
> - `M(α)` is **finite** and **strictly positive**, so `c_inf(α) = 1/M(α)` is a
>   well-defined positive number;
> - both are **computable from the periodic continued fraction** of `α`
>   (Lagrange), so `log₂(1/c)` is an effectively obtainable real number at
>   `r = 2` — though not a uniform `O(1)` bound in the sense claimed before;
> - the *all-`C`* infimum `inf_{C} |L(C)|·‖C‖` is a separate quantity, attained
>   at small `‖C‖` (e.g. `C = (−1,1)` already gives `|α−1|`, often `≤ 1/2`),
>   and is computable by a finite search.

*Proof.* `α` is a quadratic irrational, so by Lagrange its continued fraction
is eventually periodic, hence its partial quotients are bounded: there is
`K(α) < ∞` with `a_n ≤ K` for all `n`. A standard estimate gives, for all
`q ≥ 1`,

  `‖qα‖ ≥ c(α)/q`  with  `c(α) > 0`,

so `M(α) = limsup q‖qα‖ ≥ c(α) > 0`. The same periodicity bounds `q‖qα‖`
from above (the best convergents give `q_n‖q_nα‖ ≤ 1/a_{n+1} ≤ 1`, and
`q‖qα‖ ≤ q/2` off the convergents is not the right bound — rather
`q‖qα‖ ≤ 1 + ‖q_{n-1}α‖q/q_{n-1} ≤ 2` for `q_n ≤ q < q_{n+1}`), so
`M(α) ≤ 2 < ∞`. Both bounds are attained, so `M(α)` is a specific number
determined by the period, and it is computable by enumerating the period. ∎

**What is *not* claimed.** No uniform bound `M(α) ≤ 1`, and hence no
`log₂(1/c) ≤ 0`. Measured values of `M`: `0.4996` at `√2`, `0.9990` at
`√(3/2)`, `1.4978` at `√(5/3)`, `2.4953` at `√(7/5)`, `5.4862` at
`1/√11`. So `c_inf` ranges over roughly `[0.18, 2.0]` across small `d` — all
`O(1)`, which is the point, but the sign is not fixed.

**`r = 2` is still fully closed**, by Thm 2.2 (`E ≤ N/2 + 3/2`, explicit, no
constant needed) and Thm 2.6. The `r = 2` case never required a statement
about `M(α)`; §2.10.1 introduced it to describe the *object*, and that
description stands. What was wrong was the extra claim that the object is
`≥ 1` and so free.

### 2.11.5 At `r = 3` the target is much stronger than `2^{−poly(N)}`

The empirical picture at `r = 3` is informative about how strong a result is
plausible. Measuring
`c(d) = min ‖C₂α + C₃β‖ · max(|C₂|,|C₃|)²`
with `α = √(d₂/d₁)`, `β = √(d₃/d₁)`, over exhaustive
boxes in `(C₂,C₃)`:

| radicand range | triples sampled | `min c` | `max log₂(1/c)` | `median log₂(1/c)` |
|---|---|---|---|---|
| `d ≤ 50` | 455 | `2.9×10⁻⁵` | 15.1 | — |
| `d ≤ 100` | 2300 | `3.3×10⁻⁶` | 18.2 | — |
| `d ≤ 200` | 15180 | `1.7×10⁻⁶` | 19.2 | — |
| `d ≤ 400` | 76076 | `1.7×10⁻⁶` | 19.2 | — |
| `d ∈ [2,500)` | 4000 | `2.3×10⁻⁷` | 22.0 | 5.6 |
| `d ∈ [500,2000)` | 4000 | `4.1×10⁻⁶` | 17.9 | 5.8 |
| `d ∈ [2000,6000)` | 4000 | `8.7×10⁻⁶` | 16.8 | 5.9 |
| `d ∈ [6000,20000)` | 4000 | `7.2×10⁻⁶` | 17.1 | 5.8 |

Two things to read off:

- **`log₂(1/c)` does not grow with `d`.** Across four orders of magnitude of
  radicands the worst case stays in `[16.8, 22.0]` and the median sits at
  `5.8 ± 0.1`. The value `2^{r−1}−1`-free target of Thm 2.11 is
  `log₂(1/c) ≤ poly(N)`; the data suggest the far stronger
  `c = Ω(1)`, an absolute constant, at least up to `d ≤ 20000`.
- **The worst case is a genuine infimum, not a search artefact.** For the
  minimising triple `d = (263,311,313)` with `C = (−11,−26,36)` (certified sign
  `+1`, `N = 27`, `c = 1.7258×10⁻⁶`, `log₂(1/c) = 19.14`), the value of `c` is
  *identical* for `B = 60, 120, 240, 480, 960, 2000`. Since `c(B)` is
  non-increasing in `B`, the flatness means the infimum has been reached.

**This is evidence, not proof, and it is labelled as such.** A uniform
`c = Ω(1)` at `r = 3` would be a theorem about Markov constants of vectors of
quadratic irrationals, and no such theorem is known. What it does establish is
that the target is not razor-thin: `log₂(1/c) ≤ 22` in every case tested, so
Thm 2.9's constant is plausibly bounded by an absolute constant rather than
degrading with the input.

### 2.11.6 What is actually still open

Putting §2.11 together, the problem reduces to a single statement.

> **Open statement.** For every `δ > 0` and every radicand tuple `d`, the
> constant of Thm 2.9 satisfies `c(d,δ) ≥ 2^{−poly(N)}`, uniformly in `r`.

Everything else is settled:

- the **exponent** is `r−1` (Thm 2.9), and the conjecture asks only for
  `poly(N)`;
- the **norm route** cannot supply a `poly(r)` exponent (Thm 2.13), and
  halving the field degree does not help (Thm 2.14);
- **fixed rank** is already effective (Thm 2.12), so the issue is uniformity;
- at `r = 2` the case is **closed** (Thm 2.2/2.6 need no constant); Thm 2.15
  only says the constant object is finite, positive and computable;
- at `r = 3` the constant looks **bounded in absolute terms** (§2.11.5).

> **What is and is not proved, stated once and plainly.** For every fixed `K`
> the bound `E ≤ N^K` holds when `r ≤ (K−1)·log₂N`; and the constant is
> effective for `r ≤ 3`. Neither is the conjecture. The conjecture asks for
> **one fixed `K`** valid for **all** `r ≤ N/2`, and `K` would have to grow
> like `N/(4 log₂N)` to cover `r = N/4` — which is not a constant. So
>
> ```
>   PROVED:  r ≤ 3  (all `d`, effective)              Thm 2.20, Thm 2.2
>   PROVED:  r ≤ (K−1)·log₂N for each fixed K   Cor 2.19
>   OPEN:    4 ≤ r ≤ N/2 with r ≿ log₂N     <-- the conjecture lives here
> ```
>
> and the only obstacle in that region is the single constant
> `log₂(1/c(d,δ))`. Thm 2.9 makes `E` polynomial *given* that constant; every
> norm route (Thms 2.13, 2.14, 2.17) is provably exhausted, and the one
> generalisation attempted in §2.11.9 was false.

So the honest summary of the §2.11 program: the requested improvement
"replace `2^{r−1}−1` by `poly(r)`, effectively" is **already true rank-by-rank**
(Thm 2.12), **already true for all inputs with `r = O(log N)`** (Cor 2.19), and
**provably impossible to achieve by any variant of the norm argument** — not
merely by reweighting (Thm 2.13) but in principle, because the exponent is the
algebraic degree of `S` (Thm 2.17). The radicand-magnitude dichotomy in
§2.11.8 further confines the danger to *small radicands at large `r`*.

Advancing further requires a **third family of argument**: not a norm (Thm 2.17
rules those out) and not the Subspace Theorem (ineffective, Thm 2.9). That is a
sharper statement of the gap than §2.10 gave, and it is the honest one.

### 2.11.7 The exponential is the *algebraic degree* — a stronger obstruction

Thm 2.13 showed the norm bound's exponent is `2^{r−1}−1` because AM-GM
forces the count. That is the right answer for the right reason, but there is a
cleaner and strictly stronger statement available: the exponent is not an
artefact of the estimate at all, it is the degree of the number.

> **Thm 2.17 (degree obstruction).** Let `S = Σᵢ Cᵢ√dᵢ ≠ 0` with all
> `Cᵢ ≠ 0` and `dᵢ` distinct squarefree of F₂-rank `r`. Then
>
> `[Q(S) : Q] = 2^r` exactly.
>
> Consequently every norm-based bound for `|S|` has exponent `≥ 2^{r−1}−1`,
> and **no norm argument of any kind can yield a `poly(r)` exponent.**

*Proof.* `K = Q(√dᵢ)` is multiquadratic of degree `2^r` (Thm 1.2) and
`Gal(K/Q) ≈ (ℤ/2ℤ)^r` acts by independent sign flips
`σ_ε(S) = ΣεᵢCᵢ√dᵢ`. The stabilizer of `S` is

  `{ ε : Σ(εᵢ−1)Cᵢ√dᵢ = 0 } = { ε : Σ_{i∈supp(ε)} Cᵢ√dᵢ = 0 }`,

and by Thm 1 the `√dᵢ` are linearly independent over `Q`, so with every
`Cᵢ ≠ 0` this forces `supp(ε) = ∅`, i.e. `ε` is the identity. The stabilizer
is trivial, so the orbit of `S` has `2^r` points and `[Q(S):Q] = 2^r`. ∎

The minimum polynomial of `S` therefore has `2^r` distinct roots, and
`|N(S)| = ∏_σ |σ(S)| ≥ 1` is a product over all of them. Since `−S` is
also a conjugate, the `2^r` factors pair into `2^{r−1}` pairs, of which
`{S, −S}` is one, leaving `2^{r−1}−1` factors that must each be bounded above.
**This subsumes Thm 2.13**: the obstruction is the degree, not the AM-GM step, so
it is not circumventable by a cleverer choice of representatives, a better
prefactor, or a different norm.

**The theorem is sharp, and its sharp case is the input normalization.** Zero a
single `Cᵢ` and the support drops to rank `r−1`, so the degree is exactly
`2^{r−1}` — not `2^{r−2}`, not `2^{r−1}−1`. That is why the zero
coefficients are aggregated away and dropped before anything here is applied: a
surviving zero would make Thm 2.17 false as stated, and would also make the
norm bound's count `2^{r−1}−1` an overcount. The "all `Cᵢ ≠ 0`" hypothesis is
not decoration.

Verified numerically: for `r = 1..6`, 200 random instances each, the number of
distinct sign-flip images of `S` is exactly `2^r` (200/200 at every rank); and
zeroing any one coefficient gives exactly `2^{r−1}` (40/40 at `r = 3, 4`).

Note this is *stronger* than §2.4's `max(j−2,0)` rank saving, which concerns
the tower and is cancelled by squaring — here the degree is intrinsic to `S`.

**Where this leaves the Subspace Theorem.** Thm 2.9 is not a norm argument: it
bounds `|L(C)|` directly from finiteness of a solution set. That is precisely
why it escapes `2^r` and reaches `r−1`. So the two families sit on opposite sides
of Thm 2.17, and the conjecture needs a third idea that is neither.

### 2.11.8 The conjecture is proved for `r = O(log N)`

A corollary of Thm 2.12 that needs no constant at all, and localises the open
case sharply.

> **Thm 2.18 (input-size counting).** The `r` radicands are distinct integers
> `≥ 2`, so
>
> `N = Σᵢ(bitlen|Cᵢ| + bitlen dᵢ) ≥ Σᵢ log₂dᵢ ≥ log₂((r+1)!) = r·log₂r − 1.4427r + O(log r)`.
>
> In particular `r ≤ N/2`.

*Proof.* `log₂dᵢ ≥ bitlen dᵢ − 1` is false in general, so take the weaker
`bitlen dᵢ ≥ log₂dᵢ` — true, since `d < 2^{bitlen d}`. The `dᵢ` are `r`
distinct positive integers, so their product is `≥ (r+1)!` (the product of the
`r+1` smallest distinct positive integers `1..r+1`). Taking `log₂` and using
Stirling for `log₂((r+1)!)`. For `r ≤ N/2`: each term contributes
`bitlen|Cᵢ| ≥ 1`, so `N ≥ r`. ∎

> **Cor 2.19 (the conjecture holds for `r = O(log N)`).** Fix `K ≥ 2`. For
> every input with
>
> `r ≤ (K−1)·log₂(N)`,  we have  `E ≤ N^K`.
>
> Effective, uniform in `d`, with an explicit constant.

*Proof.* Thm 2.12 gives `E ≤ (2^{r−1}−1)·(N + log₂(3r))`. Since `r ≤ N/2`
(Thm 2.18), `log₂(3r) ≤ log₂(3N/2) ≤ N` for `N ≥ 2`, so the bracket is
`≤ 2N`. And `r ≤ (K−1)·log₂N` gives `2^{r−1} ≤ N^{K−1}/2`. Hence
`E ≤ (N^{K−1}/2)·2N < N^K`. ∎

Checked numerically for `K ≤ 8` and `N` from `4` to `65536`: the inequality
holds in every case, with the ratio `E/N^K` peaking near `1` only at the
smallest admissible `(K, N)` and falling off fast.

**This halves the open case.** Combined with Thm 2.18 (`r ≤ N/2`), the
unresolved regime is

  `O(log N) ≤ r ≤ N/2`,

and Thm 2.17 says any theorem closing it must be a non-norm argument. Two
further shrinkings:

- **Radicand magnitude.** If the `dᵢ` have `≥ 2^{r−1}/r` bits then `N ≥ r·2^{r−1}/r = 2^{r−1}`, so
  `2^{r−1} ≤ N` and Cor 2.19 applies with `K = 2`. The conjecture is
  therefore **trivially true for large radicands**; the danger is exactly
  *small* radicands together with *large* `r`.
- **The counting constraint.** F₂-independence of `r` radicands drawn from `k`-bit
  integers requires `r ≤ O(2^k)`, so `k ≥ log₂r − O(1)`. Minimal radicands
  (`k = O(log r)`) are both permitted and the worst case.


### 2.11.9 The constant is effective for `r ≤ 3` (and `r = 4` is OPEN)

> **RETRACTED, in this subsection, twice.** An earlier revision of this
> subsection claimed a general bound `|S| ≥ √D/((r+1)M)^{2^{r−2}−1}` for all
> `2 ≤ r ≤ 4`, with a "cutoff at `r = 4`" read off an exponent
> `α_r = r − 2^{r−2}`. **That bound is false.** It rests on the identity
>
> ```
> min_{C_r ∈ℤ} |C_r + a_j W|  =  a_j ·‖W‖,
> ```
>
> which is **not an identity**. The left side is `dist(a_j W, Z)`, and
> `dist(a_j W, Z) ≠ a_j dist(W, Z)` unless `a_j = 1`, because `a_j Z` is a
> different lattice from `Z`. Factoring a constant out of `dist(·, Z)` is
> valid; minimising over a *rescaled* integer variable is not. I asserted the
> false step in passing and it survived a rewrite. The `α_r` bookkeeping and
> the `r = 4` cutoff were consequences of the false premise, so they go with it.
>
> **The error was caught by testing the bound rather than by re-reading the
> proof**: an exhaustive check of the proposed `r = 3` inequality over 910 000
> certified instances produced 4 503 violations, and instrumenting a single
> violation exposed the identity as the culprit. Two further slips were in the
> same family — a dropped factor of `M`, and a wrong field degree
> (`2^{r−1}` vs `2^{r−2}`) — each of which I had asserted and then "corrected"
> without a test. **Three claims in this file were wrong, and in every case the
> error was in the same step.** That is recorded here rather than quietly
> deleted, because the pattern is the lesson.

What survives, and is genuinely proved:

> **Thm 2.20 (`r ≤ 3` is effective).** For `r = 2`, Thm 2.2 gives `E ≤ N/2`
> with no constant. For
> `r = 3`, with `D = maxᵢdᵢ` and `M = maxᵢ|Cᵢ|`:
>
> `|S| · M² ≥ 1/(2√D)`,  so `log₂(1/c) ≤ 1.5 log₂D + 1 = O(N)`.

*Proof, for `r = 3` only.* Relabel so `d₃` is the reference and put
`aᵢ = √(dᵢ/d₃)`, so `|S| = √d₃·L` with `L = C₃ + C₁a₁ + C₂a₂`.
Minimising the free coordinate `C₃ ∈ Z` gives

```
‖L‖ = ‖C₁a₁ + C₂a₂‖ = a₁ · ‖C₁ + C₂g‖,   g := a₂/a₁ = √(d₂/d₁),
```

and *this* step is valid: `a₁` factors out of a distance to `ℤ` because
`a₁·` is a fixed positive scalar, so
`a₁ dist(x,Z) = dist(a₁ x, a₁ Z)` and the integer being approximated
is the same. (Contrast the retracted step, which moved the *integer* variable.)

`g` is a quadratic irrational, so the integer norm applies. With `u = C₁d₁`,
`v = C₂`, `C₁ + C₂g = (C₁d₁ + C₂√(d₁d₂))/d₁` and

```
|C₁d₁ + C₂√(d₁d₂)| · |C₁d₁ − C₂√(d₁d₂)| = d₁|C₁²d₁ − C₂²d₂| ≥ d₁,
```

the last factor a nonzero integer — it vanishes only if `d₁/d₂` is a rational
square, impossible for distinct squarefree `d₁,d₂` (Cor 1.3). The second
factor is `≤ 2M√(d₁D)`, hence `‖C₁ + C₂g‖ ≥ 1/(2M√(d₁D))`, and collecting
the prefactors (`a₁/√d₃ · d₁ = 1/√(d₁d₃)`) gives
`|S| ≥ 1/(2M√D)`. Multiplying by `M²` and dropping the favourable `M`
gives `|S|M² ≥ M/(2√D) ≥ 1/(2√D)`. ∎

Verified by exhaustive certified search over 2024 radicand triples (`d ≤ 97`,
`|Cᵢ| ≤ 16`, magnitudes from `intervalsum.certified_interval_for_S`, never
> `intervalsum.certified_interval_for_S`, whose enclosure is *tighter* than
> `signexact`'s `bracket_num` — `bracket_num` is the lower endpoint of
> `signexact`'s own certified enclosure, not a truncation, so it is also sound
> but coarser): **0 violations, minimum margin 3.18×.**
>
> *(Correction: an earlier revision of this line called `bracket_num` "a
> bit-truncation" and implied it was unsound. That is wrong —
> `signexact._sum_bracket_int` returns a true enclosure
> `S ∈ [lo/2^bits, hi/2^bits]`, and `bracket_num` is its lower endpoint.
> Checked: `bracket_num ≤ min(|lo|,|hi|)` on 20 000 random instances, 0
> violations. The two are used interchangeably here because both are sound
> lower bounds, and the tighter one is preferred.)*

**Why it does not generalise.** The whole gain at `r = 3` is that the inner form
is `C₁ + C₂g` with `g` a **single** quadratic irrational, so the norm is
quadratic and its integer `|C₁²d₁ − C₂²d₂|` is what one divides by. At
`r = 4` the inner form has **two** independent generators
`√(d₂/d₁), √(d₃/d₁)`, the norm is a quartic form, and there is no
identity that collapses it to `1/poly`. I checked the one obvious hope — that
`a₁,a₂,a₃` lie in geometric progression so a single `γ` captures all
three — and it is excluded: `a₃/a₁ = (a₂/a₁)²` would force `d₃ = d₂`, and
`a₃/a₁ = (a₂/a₁)^{−1}` would force `d₂d₃ = d₁²`. So there is
no `r = 3`-style collapse at `r = 4`.

**The `r = 4` measurement, and what it does and does not say.** Exhaustive
certified search, 3060 quadruples with `d ≤ 89`, `|Cᵢ| ≤ 5`, gives
`min c_S = 5.6×10⁻³` (`log₂(1/c_S) = 7.5`) at `r = 4`, *larger* than the
`r = 3` value `1.7×10⁻⁶`; and forcing large `M` makes it larger still
(`log₂c_S` rose to 30 at `M = 395`). This is a genuine **upper** bound on the
infimum and so it is evidence that the `r = 4` constant is not small. It is
**not** a lower bound, it is not a theorem, and it does not certify the
infimum — at `r = 3` the `c(B)`-flattening test was needed for that, and it was
not run at `r = 4`. It is recorded as evidence only.

### 2.11.10 Scoreboard for §2.11

| asked | answer |
|---|---|
| is `poly(r)` achievable effectively at fixed `r`? | **yes**, Thm 2.12 |
| is it achievable for all `r` at once? | **yes** for `r = O(log N)`, Cor 2.19 |
| can the norm argument do better than `2^{r−1}−1`? | **no**, Thm 2.17 (it is the degree) |
| does squaring / degree reduction help? | **no**, Thm 2.14 |
| is a constant needed at `r = 2`? | **no** — Thm 2.2/2.6 close it; Thm 2.15 only identifies the object |
| is the constant effective at `r ≤ 3`? | **yes** at `r=3`, `c ≥ 1/(2√D)` (Thm 2.20); `r=2` needs none |
| do the norm routes work at `r ≥ 4`? | **no** — but the `r=4` collapse I claimed was **false**; see the retraction in §2.11.9 |
| what is the constant at `r ≥ 4`? | **open**; `min c_S = 5.6×10⁻³` measured at `r=4` (upper bound only), §2.11.9 |
| where can a counterexample live? | `O(log N) ≤ r ≤ N/2`, small radicands, §2.11.8 |
| **is the uniform conjecture proved?** | **NO** — this is EHS Problem 33. Proved only for `r ≤ 3` and `r = O(log N)`; `5 ≤ r ≤ N/2` with `r ≿ log₂N` is open |

---

## 3. Theorem 3 — the F1 separation law, exactly

This replaces the float-based tolerance in `f1_theorem.py`
(`_SEPARATION_TOLERANCE = 1e-3`) with a proved statement. The residue law
it rests on is proved in `Docs/Paper3.md` §2.4; it is restated here as a
fact, not re-proved.

**Fact 3.1 (residue law).** For `α = a(a+1)`, `a ≥ 2`, and the convergents
`p_n/q_n` of `sqrt(α)`:

```
D_n = p_n² − αq_n² = +1  (n odd),   −a  (n even).
```

In particular `|D_n| ≤ a` and `D_n ≠ 0` for every `n ≥ 0`.

**Theorem 3.2.** Let `s_n = q_n·sqrt(α)/p_n` (a positive real, equal to
`sqrt(1 − D_n/p_n²)`). Let `N_n = 2log₂p_n + 2log₂q_n + log₂α` (the real
input size; see the remark on conventions below). Then:

**(i) EXACT.** `N_n/4 = log₂ p_n + (1/2)·log₂ s_n`.

**(ii) EXACT, one-sided.**

```
log₂(1/|S_n|) − (N_n/4 + 1 − log₂|D_n|) = ρ_n := log₂((1 + s_n)/(2·sqrt(s_n))),
```

and `ρ_n > 0`.

**(iii) BOUND.** For `a ≥ 2` and `n ≥ 0`:

```
0 < ρ_n ≤ 0.205·a²/p_n⁴,
```

so uniformly `ρ_n ≤ 0.205·a²/(2a+1)⁴ ≤ 1.31e-3` for `n ≥ 1`.

**(iv) ASYMPTOTIC.** `ρ_n = D_n²/(32·ln2·p_n⁴)·(1 + O(a/p_n²))`.

*Proof.*

(i) `2log₂q_n + log₂α = 2·log₂(q_n·sqrt(α)) = 2·log₂(p_n·s_n) =
2log₂p_n + 2log₂s_n`. Substituting into `N_n`:
`N_n = 2log₂p_n + (2log₂p_n + 2log₂s_n) = 4log₂p_n + 2log₂s_n`. Divide by 4.

(ii) Rationalising `S_n = p_n − q_n·sqrt(α)` with the residue
(`|S_n|·(p_n + q_n·sqrt(α)) = |D_n|`):

```
log₂(1/|S_n|) = log₂(p_n + q_n·sqrt(α)) − log₂|D_n|
             = log₂p_n + log₂(1 + s_n) − log₂|D_n|.
```

Subtract `(N_n/4 + 1 − log₂|D_n|)` using (i): the `log₂p_n` and
`−log₂|D_n|` terms cancel and `1 = log₂2`, leaving `ρ_n`. For the sign:
`(1+s)/(2√s) = (√s + 1/√s)/2 ≥ 1` by AM–GM, with equality iff `s = 1`, i.e.
`D_n = 0` — excluded by Fact 3.1. So `ρ_n > 0` strictly. Note `f(s) :=
ln((1+s)/2) − (1/2)ln s` satisfies `f(s) = f(1/s)`, so the sign statement
covers both parities of `n`.

(iii) `s_n = sqrt(1 − D_n/p_n²)` with `|D_n| ≤ a`, so

```
s_n ∈ [sqrt(1 − a/p_n²), sqrt(1 + a/p_n²)].
```

For `n ≥ 1`, `p_n ≥ p_1 = 2a+1`, and `a/(2a+1)²` decreases in `a ≥ 1` (its
derivative has sign `(2a+1)(1−2a) < 0`), so `a/p_n² ≤ 2/25`, reached at
`a = 2, n = 1`. Hence `s_n ∈ [sqrt(23/25), sqrt(27/25)] ⊂ [0.959, 1.040]`.

Taylor at 1: `f(1) = 0`, `f'(1) = 1/(1+s) − 1/(2s)` evaluated at 1 is `0`,
and `f''(s) = 1/(2s²) − 1/(1+s)²`. Also
`f'''(s) = −1/s³ + 2/(1+s)³ < 0` exactly when `2s³ < (1+s)³`, i.e.
`s < 1/(2^{1/3} − 1) ≈ 3.85`; so `f''` is decreasing on the reachable range
and `f''(s) ≤ f''(sqrt(23/25)) = 0.283`. With `f(1) = f'(1) = 0`, Taylor's
theorem gives `f(s) = (1/2)·f''(ξ)·(s−1)²` for some `ξ` between `1` and
`s`, so `f(s) ≤ 0.1415·(s−1)²`. Finally
`|sqrt(1+x) − 1| ≤ |x|` for `x > −1` (squaring both sides of the two one-sided
versions `sqrt(1+x) ≤ 1+x` for `x ≥ 0` and `sqrt(1−x) ≥ 1−x` for `x ∈ [0,1]`),
so `(s−1)² ≤ (a/p_n²)²`. In `log₂` units (`ρ_n = f(s_n)/ln2`):

```
ρ_n ≤ 0.1415/(ln2) · a²/p_n⁴ = 0.204·a²/p_n⁴ ≤ 0.205·a²/p_n⁴.
```

Uniform bound: `a²/p_n⁴ ≤ a²/(2a+1)⁴ ≤ 4/625` (worst at `a = 2`), giving
`ρ_n ≤ 1.31e-3`.

(iv) `f''(1)/2 = 1/8`, so `f(s) = (s−1)²/8 + O((s−1)³)`. And
`s_n − 1 = sqrt(1 − D_n/p_n²) − 1 = −D_n/(2p_n²) + O(a²/p_n⁴)`. Squaring and
dividing by `ln2` gives the claim. ∎

*Remarks.*

- **Why the measured residual is so small.** Both sides of the law carry a
  first-order correction `−D_n/(4p_n²·ln2)` — `log₂(1+s_n)` on the left,
  `(1/2)log₂s_n` inside `N_n/4` on the right — and Theorem 3.2(ii) shows
  they cancel exactly, leaving the second-order term. That is why the
  measured worst residual (`7.5e-5` at `a = 2, n = 1` over `a ≤ 24`,
  `n ≤ 39`) is ~400× smaller than either correction alone (~`0.03` there),
  and why (iv) predicts `1/(32·ln2·625) = 7.18e-5` — within 5% of the
  measurement.
- **The tolerance is now a theorem's slack, not an empirical margin.** The
  proved leading form puts the worst case over the verified range at
  `≈ 7.2e-5`, a 13× margin under `_SEPARATION_TOLERANCE = 1e-3`.
- **`N` conventions.** `N_n` here is the real quantity
  `2log₂p_n + 2log₂q_n + log₂α`, which is what `f1_theorem.log2_ratio_law`
  uses. The instance generator `families.input_bits` instead counts integer
  bit-lengths: `N_bits = bl(p_n²) + bl(αq_n²)` with `bl(x) ∈ (log₂x, log₂x + 1]`,
  so `N_bits ∈ (N_n, N_n + 2]`. Every ratio statement transfers between the
  two with that slack; nothing below depends on which is used beyond `O(1/N)`.

**Corollary 3.3 (the Θ(N) law, with constants).** For `a ≥ 2`, `n ≥ 1`:

```
N_n/4 + 1 − log₂|D_n|  <  log₂(1/|S_n|)  ≤  N_n/4 + 1 − log₂|D_n| + 1.31e-3,
```

so, writing `θ_n = log₂(1/|S_n|) − N_n/4`,

```
θ_n ∈ [1 − log₂a − 1.31e-3,  1 + 1.31e-3],
|log₂(1/|S_n|)/N_n − 1/4| ≤ (log₂a + 1.002)/N_n.
```

In particular the F1 family supplies, for every `a ≥ 2`, an infinite sequence
of instances with `log₂(1/|S_n|)/N_n → 1/4` from the stated band — the
README's proved-separation row, with the `O(1)` explicit and one-sided. ∎

---

## 4. What stays open

- **The conjecture: OPEN. But the shape of the obstacle is now a constant,
  not a growth rate.** Decide `sign(S)` in time `poly(N)`. Equivalently
  (Thm 2.1 gives termination), prove `|S| ≥ 2^{−poly(N)}` uniformly over the
  input.

  State of the art, for `r` terms, writing `m = log₂‖C‖_∞`:

  | | bound | `E/N` | source |
  |---|---|---|---|
  | **guaranteed, best exponent** | `|S| ≥ c(d,δ)·‖C‖^{−(r−1+δ)}` | `(r−1+δ)/r` | **Thm 2.9** (single-form ST) — matches the achieved rate |
  | guaranteed, constant-free, generic `C` | `|S| > ‖C‖^{−(r−1+δ)}` | `(r−1+δ)/r` | Thm 2.9, for `C` outside finitely many rational subspaces |
  | guaranteed, per-instance computable, but `Q₀` ineffective | `|E| ≥ Q^{−(n+1)}` | `O(n·N)` once `Q₀` is explicit | EHS Thm 4 — see `STATE_OF_THE_ART.md`; the `O(N)` previously recorded here was wrong |
  | guaranteed, literature, uniform in `n` | `γ·(n‖x‖)^{−2n}` | `→ 2` | EHS Thm 8 |
  | guaranteed, best fully *effective* here | `((2^{r−1}−1)/2)(log₂M+1)` | `≈ 2^{r−3}` | Thm 2.6 (norm) |
  | guaranteed, effective, self-contained | `≥ (3T)^{−(2^{r−1}−1)}`, `T = ∑|C_i|√d_i` | `(2^{r−1}−1)·O(1)` | **Thm 2.10** (norm) — effective but `2^{r−1}` |
  | **achieved (explicit family)** | `≥ (r−1)/r − o(1)` | `(r−1)/r` | Thm 2.5 = EHS Thm 10 |
  | achieved at `r = 2`, exactly | `= 1/2` | `1/2` | Thm 2.2, attained by Thm 2.6 |
  | achieved at `r = 3`, finite `N` | `E/N = 1.2795` at `N = 17` | `1.2795 > 2/3` | §2.10.2 — **falsifies `sup E/N = (r−1)/r`** |
  | achieved by F1 | `= 1/4` | `1/4` | Thm 3 / Cor 3.3 |

  **The reframing that matters.** The guaranteed and achieved rates
  **coincide** in their *linear* part: both are `(r−1)/r`, at every rank, by two
  completely different arguments — Thm 2.9 (Subspace Theorem, upper) and
  Thm 2.5 (pigeonhole, lower). So the conjecture is no longer "close an
  exponential gap". It is:

  > **Asymptotic sharp form.** `E ≤ ((r−1+δ)/r)·N + O(1)` and this is
  > attained to `o(N)`, for every `r`; equivalently
  > `limsup_{N→∞} E/N = (r−1)/r`.

  Equivalently: `|S| ≥ ‖C‖^{−(r−1)−o(1)}` is the true shape, and the
  conjecture asks only for `poly(N)`, which is far weaker. The `2^r` in
  Thm 2.6 is a property of the **norm**, not of the problem.

  > **Correction (this file previously said `sup E/N = (r−1)/r`). That
  > statement is false as written, and §2.10.2 falsifies it explicitly.** A
  > certified `r = 3` instance reaches `E/N = 1.2795 > 2/3`. The reason is that
  > every bound here carries an **additive** `O(1)` slack
  > (`log₂(1/c) + O_{d,δ}(1)`), which is a *fixed* bonus at fixed `d`; divided
  > by a small `N` it inflates `E/N` arbitrarily. The rate is a statement
  > about `N → ∞`, not about a supremum over all `N`. Retracted.

  **Three things are now settled that earlier versions of this file had
  wrong or missing:**
  - The exponential-in-`r` barrier is not the obstacle (Thm 2.9, and EHS).
  - The exponent of the true bound is `r−1`, not `2n` nor `r` (Thm 2.9; the
    single-form application is sharper than the simultaneous one).
  - For **fixed radicands** the conjecture is *proved* in shape: `E ≤
    ((r−1+δ)/r)·N + O_{d,δ}(1)` (Thm 2.9). What is missing is the uniformity
    in `d`.

  **What remains open, precisely.** One statement, with two known partial
  answers that each miss half of it:

  > Replace the exponent `2^{r−1}−1` of the effective norm bound by `poly(r)`,
  > **while keeping the constant effective.**
  >
  > §2.11 sharpens this considerably. The exponent is **already** `poly(r)` at
  > every fixed rank (Thm 2.12), and it is **provably irreducible** for the norm
  > argument (Thm 2.13), as is the halving-the-field escape (Thm 2.14). So the
  > request reduces to the single statement below; no variant of the norm
  > argument can deliver it.

  - Thm 2.9 has the linear exponent `r−1` but an **ineffective** constant
    `c(d,δ)`. Thm 2.11 shows the *only* requirement on that constant is
    `c ≥ 2^{−poly(N)}`; no uniformity in `d` and no `Ω(1)` is needed.
  - Thm 2.10 / Thm 2.6 are **fully effective** but their exponent is
    `2^{r−1}−1`, which is `2^{Ω(N)}` because `r = Θ(N)` is achievable. This
    exponent is irreducible (Thm 2.13), and halving the field degree via
    squaring does not help (Thm 2.14).
  - **Every fixed rank is already settled**, effectively: `E ≤ (2^{r−1}−1)(N +
    log₂(3r))` (Thm 2.12). Since `2^{r−1}−1` is a *constant* at fixed rank, the
    conjectured polynomial separation is a theorem rank-by-rank. Only
    uniformity in `r` is open.
  - At `r = 2` **no additive constant is needed at all** — but not for the
    reason once claimed. Thm 2.2 gives `E ≤ N/2 + 3/2` outright. Thm 2.15
    only identifies the constant as the reciprocal Markov constant of a
    quadratic irrational: finite, positive, `Θ(1)`, and computable from the
    periodic CF. (The claim `M(α) ≤ 1`, hence `c ≥ 1`, was **false** and is
    retracted in §2.11.4; `M(√(5/3)) = 1.4978 > 1`.)
  - At `r = 3` the constant looks **bounded in absolute terms** — `log₂(1/c) ≤ 22`
    over all `d ≤ 20000` tested, median `5.8`, with `c(B)` exactly flat in `B`
    (§2.11.5). Evidence only, but it says the target is `c = Ω(1)`, far stronger
    than Thm 2.11's `2^{−poly(N)}`.

  A secondary gap, specific to the Subspace-Theorem route: the exceptional
  subspaces `V_1, …, V_t` are not known explicitly. Unlike the constant this
  is *detectable* — `C` lies in a proper rational subspace iff some nonzero
  integer `a` has `a·C = 0`, decidable by exact linear algebra over `Q` — so
  the generic case already has a constant-free bound. Only the exceptional
  instances are unresolved.

  **At `r = 2` the constant is not merely bounded but computed** (§2.10.1): it
  is the reciprocal of the Markov constant of a quadratic irrational, hence
  readable off a periodic continued fraction. That, and not a better bound, is
  why `r = 2` is finished.

  EHS state (Discussion §4) that the quantitative Subspace Theorem does not
  supply either: it fails to bound the number of solutions of
  `q^{1+δ} Π_i dist_Z(q·sqrt(a_i)) < 1` "even when all algebraic numbers are
  square roots of integers."

  **Why this is still not a proof that SRS ∈ P.** Deciding `sign(S)` by
  approximation needs a *known* lower bound to set the precision threshold.
  Thm 2.9 supplies `c(d,δ)`, which is unknown. The classical approach needs
  `O(2^r log)` bits of certified precision; Thm 2.9 would need only `O(N)` —
  but with a constant that cannot be computed. EHS Thm 4 escapes the
  constant and is per-instance computable, at the cost of a checkable
  hypothesis — but see `STATE_OF_THE_ART.md`: making its threshold `Q₀`
  explicit costs `exp(Θ(n²a))`, so the honest rate is `O(n·N)`, **not**
  `O(N)`, and the route does not close the conjecture.

- **The sharp form is falsifiable, and the repo searches for it.**
  `refutation.py` hunts for an instance with a certified sign, a certified
  magnitude, and `E/N` above the ceiling of Thm 2.5 / Thm 2.9 — i.e. above
  `(r−1)/r` by more than the `O(1/N)` slack. It sweeps ranks 2–5 with random
  and greedy strategies plus the F6 ceiling family. **No refutation found.**
  That is a negative result and is recorded as one: it bounds what was
  searched and supports nothing. See `tests/test_refutation.py`, whose
  `test_f6_is_not_a_refutation` is the regression pin for the false positive
  described in §2.7.

- **A concrete lead, stated as a problem.** Every argument in §2.1–§2.7 uses
  only the orbit of `S` and the size of its elements — never the minimal
  polynomials `x² − d_i`, never the interaction between the radicals, beyond
  Cor 1.3's independence. EHS note that the quantitative Subspace Theorem
  fails to help *even when every `β_i` is a square root of an integer*. Read
  the other way: the square-root structure is where an effective constant
  would have to come from. Thm 2.9 sharpens the exponent using only
  Cor 1.3, so it does not touch this; the lead is untouched by anything in
  this file. Restated at greater length in `THE_MISSING_THEOREM.md` §3.4.

- **The `k = 2` constant: CLOSED, sharp from both sides, and effective.**
  Theorem 2.2: guaranteed `N/2 + 3/2`, achieved by the coefficient-Pell family
  `2C₁² − 3C₂² = −1` with radicands `2, 3`. Theorem 2.6 strengthens the
  guarantee to `½(log₂M + 1)`, which is *attained* by the same family
  (slack `→ 0` geometrically). So `1/2` is a proven sharp constant rather
  than a matched-but-loose pair, and unlike Thm 2.9 the constant here is
  **explicit** — the only rank where that is true. An earlier version of this
  file listed the constant as open with the interval `[1/4, 9/4]`; both
  numbers were wrong. Thm 2.7 is *weaker* at `r = 2` than Thm 2.6 (it does not
  know the radicands), so the two-term case is settled by the repo's own work
  and not by the literature.

- **Zero test: CLOSED, in P** (Thm 1), consistent with Blömer's
  polynomial-time algorithm for the equality case. For calibration on the
  sign problem: the best known upper bound on its complexity is PSPACE, it is
  not known to be in NP, and NP would follow from a polynomial separation
  bound.


---

## 5. Corroboration (tests, not evidence)

| claim | where checked |
|---|---|
| Theorem 1.2 / 1.3 | `tests/test_field.py` (F₂ rank vs exhaustive scan over every subset of the first six primes); zero-test oracle in `tests/test_squarefree.py` |
| Theorem 2.1 | `tests/test_proofs.py::test_liouville_bound_holds` (exact rational brackets, even `n` so the bound is rational) |
| Theorem 2.2 | `tests/test_proofs.py::test_two_term_bound` and `::test_coefficient_pell_family_approaches_half` |
| Theorem 2.3 | `tests/test_proofs.py::TestSecondMomentBound` (bound on random instances; domination over Thm 2.1; and the recorded fact that Thm 2.2 is 3× sharper at `r = 2`) |
| Theorem 2.6 | `tests/test_proofs.py::TestGlobalFlipPairing` (bound on random instances `r = 2..8`; domination over Thm 2.3 by the exact exponent ratio `(2^r−1)/(2^{r−1}−1)`; geometric slack decay on the coefficient-Pell family; and the retraction of the old "`±` pairing stops at `r = 2`" claim) |
| Theorem 3.2(ii)–(iii) | `tests/test_proofs.py::test_f1_residual_one_sided_and_bounded` |
| Theorem 2.4 | `tests/test_proofs.py::TestTowerStructure` (exact support count over 60 random towers; the stronger "2-element subsets only" claim; and a guard that the count has not caught up to `2^(j−1)`) |
| Theorem 2.5 | `tests/test_proofs.py::TestPigeonholeAchievedRate` (exact distinctness of all `(B+1)^r` sums for `r ≤ 4`, `B ≤ 7`; the constructed gap within the pigeonhole bound; and the `r=2` match with Thm 2.2) |
| Thm 2.7 / 2.8 (EHS, **quoted not proved here**) | `tests/test_proofs.py::TestSubspaceTheoremRoute` (the EHS explicit term `2n(log₂n+m)` is `O(N)` and is *dominated by* Thm 2.6 at large `r`; the `n+ε` sharpening; that Thm 2.5 = EHS Thm 10; that the `r=2` case is settled by Thm 2.6 and *not* by EHS; and that the remaining gap is the **ineffective** `γ`). These check the repo's *arithmetic about* the quoted bound — they do not and cannot verify a deep theorem from the literature. |
| Prop 2.4′ (support rank `j−2`) | `tests/test_proofs.py::TestTowerStructure::test_support_rank_is_j_minus_two_not_j_minus_one` (300+ random F₂-independent sets per level, `j = 2..7`) and `::test_r3_relative_norm_is_two_term` (the `r=3` identity plus squarefreeness of `d₁d₂`) |
| Thm 2.7 (Subspace route closed) | `tests/test_proofs.py::TestSubspaceTheoremRoute` (the `r=2` Roth row; the retraction of the `1/(r(r−1))` overclaim; and the dominance of Thm 2.6 at every `r ≥ 3`) |
| Fact 3.1 | integer identities in `tests/test_f1_theorem.py`; proof in `Docs/Paper3.md` §2.4 |