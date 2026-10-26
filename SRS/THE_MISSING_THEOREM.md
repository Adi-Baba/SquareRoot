# SRS — The Exact Missing Theorem

The roadmap for the Square-Root Sum Problem. One file that states, precisely:

1. what the problem is,
2. what is proved on both sides, and by whom,
3. **the exact missing theorem** that would settle it,
4. every route tried, and which ones are closed,
5. every error made and retracted along the way,
6. what a proof would have to contain.

| file | role |
|---|---|
| **`STATE_OF_THE_ART.md`** | the literature check, and the correction that changed the answer. **Read first.** |
| `PROOFS.md` | the actual proofs. Where this file and that one disagree, `PROOFS.md` wins. |
| `README.md` | what the code does. |

---

## 0. Status

`sign(S)` for `S ≠ 0` is **open**. The obstacle is **not** the exponential
barrier that earlier versions of this file described, and it is not a growth
rate at all.

**The exponent is settled.** `PROOFS.md` Thm 2.9 applies Schmidt's Subspace
Theorem to the **single** linear form that SRS actually is, and gets exponent
`r−1` on `‖C‖_∞`:

```
|S|  ≥  c(d,δ) · ‖C‖_∞^{−(r−1+δ)}        ⇒   E/N  ≤  (r−1+δ)/r
```

This is sharper than Eisenbrand–Haeberle–Singer (arXiv:2312.02057; SoCG
2024), whose simultaneous-approximation application gives exponent `2n`, and
it **coincides with the achieved rate**. Thm 2.5 (pigeonhole) supplies
instances at `E/N ≥ (r−1)/r − o(1)`. So:

> **The sharp form of the conjecture is asymptotic:**
> `limsup_{N→∞} E/N = (r−1)/r`, at every rank. The guaranteed and achieved
> rates meet. There is no exponential gap left to close.
>
> **Retracted:** these files previously stated this as `sup E/N = (r−1)/r`,
> a supremum over *all* `N`. That is **false** and is falsified explicitly in
> `PROOFS.md` §2.10.2 by the certified `r = 3` instance `d=(2,3,13)`,
> `C=(−10,−98,51)`, for which `E/N = 1.2795 > 2/3`. The excursion is an
> additive slack `log₂(1/c) ≈ 21.75` divided by a small `N = 17`; it is
> invisible as `N → ∞`.

**What remains is one constant.** Thm 2.9's `c(d_1,…,d_r,δ)` depends on the
radicands, which are input — the same γ-gap as EHS, inherited unchanged.
Thm 2.9 improves the exponent and does not touch the constant.
**§2.11 narrows this further:** every *fixed* rank is already settled
effectively (Thm 2.12), the norm exponent is provably irreducible
(Thm 2.13), `r = 2` needs no constant at all (Thm 2.15), and at `r = 3` the
constant measures `≤ 22` bits over all `d ≤ 20000`. The single open statement
is a **uniform** `c(d,δ) ≥ 2^{−poly(N)}` across `r = Θ(N)`.

**Two things follow that are worth stating plainly.**

- The `2^r` in this repo's Thm 2.6 is a property of the **norm**, not of the
  problem. Thm 2.6 remains the sharp statement about what a norm argument can
  do, and the best fully *effective* bound here.
- For **fixed radicands** the conjecture is *proved in shape*:
  `E ≤ ((r−1+δ)/r)·N + O_{d,δ}(1)`. What is missing is uniformity in `d`.

**No proof of the conjecture is claimed or implied anywhere in this
repository.** The refutation harness in `refutation.py` searches actively for
a counterexample to the sharp form and has found none.

---

## 1. The problem

Input: integers `c_1, …, c_k` and squarefree `a_1, …, a_k`. Let

```
S  =  Σ_j  c_j · sqrt(a_j)
```

Reduce by squarefree decomposition and coefficient aggregation
(`squarefree.py`) to

```
S  =  Σ_{i=1}^{m}  C_i · sqrt(d_i),        d_i  distinct squarefree, C_i ≠ 0.
```

Let `r ≤ m ≤ k` be the **F₂-rank** of `{d_i}` — the F₂-dimension of the prime
parity vectors `v(d_i)`. By Thm 1.2 this equals `[K : Q]` for
`K = Q(sqrt(d_1),…,sqrt(d_m))`.

**The precondition that makes it non-vacuous.** With unbounded input,
`sqrt(M²+1) − M ≈ 1/(2M)`, so "small separation" is trivially available. The
repo therefore fixes `a_i, b_j < 2^n`, `k` terms, total bit size `N = k·n`. All
statements below are in these terms.

**The quantity of interest.** With `E = log₂(1/|S|)` for `S ≠ 0`:

| question | status |
|---|---|
| `S = 0`? | **solved, in P** — Thm 1: distinct squarefree radicals are ℚ-linearly independent. Matches Blömer's polynomial-time algorithm for the equality case. |
| `sign(S)`, `S ≠ 0`? | **open.** Best known complexity class containing it is PSPACE; not known to be in NP. |
| separation `\|S\| ≥ ?` | polynomial in `n` with an ineffective constant (§3); conjectured `2^{−poly(N)}` uniformly. |

The conjecture is equivalent to a uniform polynomial separation bound: given
Thm 2.1's termination guarantee, `sign(S)` is decidable in `poly(N)` **iff**
`|S| ≥ 2^{−poly(N)}` for all `S ≠ 0` of input size `N` — with the constant
depending on nothing, in particular not on the radicands.

---

## 2. What is proved, on both sides

### Guaranteed (upper bounds on the achievable `E/N`)

| bound | statement | source |
|---|---|---|
| Thm 2.1 | `E ≤ (2^k − 1)(3n/2 + log k)` | Liouville on the norm |
| Thm 2.3 | `E ≤ ((2^r − 1)/2)(log₂M + 1)`, `M = Σ C_i²d_i` | exact orbit second moment `Σ_σ σ(S)² = 2^r M` |
| Thm 2.6 | `E ≤ ((2^{r−1} − 1)/2)(log₂M + 1)` | global-flip pairing — **the optimum of the norm class**, and the best fully *effective* bound here |
| Thm 2.7 = EHS Thm 8 | `E ≤ 2n(log₂ n + m) + log₂(1/γ)` | EHS, **quoted**: simultaneous-approximation application, exponent `2n`; `γ` ineffective and a function of the input |
| **Thm 2.9** | **`E ≤ ((r−1+δ)/r)·N + O_{d,δ}(1)`** | **single-form Subspace Theorem, exponent `r−1` — matches the achieved rate** |

**Read the last three rows together.** Thm 2.9 is the repo's own application
and it is the sharpest of the three, because SRS is a *single* linear form and
a single form in `r` variables has Subspace-Theorem exponent `r−1`, whereas
EHS's simultaneous application pays exponent `2n`. At `r = 3`:

| | exponent on `‖C‖_∞` | `E/N` | vs ceiling `2/3` |
|---|---|---|---|
| Thm 2.6 (norm) | — | `1.5` | `2.25×` |
| EHS Thm 8 | `2n` | `1.0` | `1.50×` |
| **Thm 2.9** | **`r−1`** | **`0.667 + δ/3`** | **`1×` as `δ → 0`** |
| achieved, Thm 2.5 | `r−1` | `0.667` | — |

Thm 2.7 still beats Thm 2.6 by an exponential factor at large `r` — measured
ratio of Thm 2.6 to Thm 2.7 at `r = 3, 5, 8, 10, 15, 20` is
`0.98, 2.4, 12, 37, 788, 1.9e4`, growing like `2^{r−3}/r` — so `O(2^r)` is
not the state of the art, and Thm 2.9 improves on Thm 2.7 as well.

**Why Thm 2.6 is the norm optimum** (exact, not heuristic). The `2^r`
conjugates pair as `{σ(S), −σ(S)}` under the global flip
`τ: sqrt(d_i) ↦ −sqrt(d_i)`, so `Π_σ σ(S) = ±(Π_half σ(S))²` and
`Π_half |σ(S)| ≥ 1` is *free*. AM–GM over the `2^{r−1}−1` remaining half-orbit
elements halves the exponent by a uniform factor

```
(2^r − 1) / (2^{r−1} − 1)  →  2        (measured 3.000, 2.333, 2.143, 2.067, 2.032, 2.016, 2.008 at r = 2..8)
```

The gain does **not** decay with `r`, and after the flip no further pairing is
available: the remaining `2^{r−1}−1` elements carry no forced symmetry. So
`((2^{r−1}−1)/2)(log₂M+1)` is final *for this class*.

### Achieved (explicit hard instances)

| rate | family | source |
|---|---|---|
| `1/4` | F1: `p_n − q_n sqrt(α_n)`, `α = a(a+1)` | Thm 3 / Cor 3.3, exact with constants |
| `1/2` | coefficient-Pell: `2C₁² − 3C₂² = −1`, radicands `2, 3` | Thm 2.2 |
| `(r−1)/r − o(1)` | pigeonhole on `(B+1)^r` sums | Thm 2.5 |

### The bracket

```
achieved:    E/N  =  (r−1)/r  +  o(1)        linear in r, polynomial in N
guaranteed:  E/N  =  O(1)  +  log2(1/gamma)/N      polynomial, but gamma unknown
best effective in this repo:  E/N  =  O(2^r)
```

Both ends are load-bearing:

- **The achieved side is optimal, not weak.** `(r−1)/r` is the **Dirichlet
  exponent**, and it is attained by pigeonhole — Thm 2.5 is *the same theorem*
  as EHS Thm 10 (their `n·max√a_i` prefactor is dominated by our
  `Σ√d_i`, same exponent `n−1`). So no adversary will find something worse.
  The conjecture is not threatened from below.
- **The guaranteed side is polynomial in the literature.** An earlier version
  of this file put `O(2^r)` here and called the norm optimum the state of the
  art. That was wrong. The residual gap is the constant `γ`, not the growth
  rate.

The conjecture asks only for `E/N ≤ poly(1)`. Both sides are already
polynomial. The entire remaining difficulty is one unknown constant.

### The `k = 2` case: closed and sharp, and it is this repo's own result

Worth singling out — it is the only rank where everything meets, and it is
where the repo is *stronger* than the literature.

- **Achieved `= 1/2` exactly** — coefficient-Pell `2C₁² − 3C₂² = −1` (Thm 2.2).
- **Guaranteed, tight from above** — Thm 2.6 at `r = 2` is `½(log₂M + 1)` with
  `M = X² + Y²`, which *is* Cauchy–Schwarz `X + Y ≤ sqrt(2(X²+Y²))` in the
  exponent, with equality iff `X = Y`. The Pell family is precisely
  `X ≈ Y`, and the measured slack decays geometrically:
  `7.3e−3, 7.7e−7, 8.0e−11, 8.3e−15, 8.7e−19, 9.0e−23, 9.4e−27`.
- **Stronger than EHS at this rank.** At `r = 2, 3` the two bounds are
  comparable and Thm 2.6 is in fact sharper, because the EHS bound does not
  know the radicands. The constant `1/2` should be attributed to this repo,
  not to the literature.
- **Caveat, stated because it is easy to overclaim.** The *rate* reaches `1/2`
  only asymptotically and from above: `0.6425, 0.5384, 0.5223, 0.5158, 0.5122,
  0.5099` for `t = 0..5`. The guarantee is `N/2 + O(1)`, not `≤ N/2` — at
  `t = 0` (`N = 2.58`) it is `1.66` against `N/2 = 1.29`. Sharp in the limit,
  loose at every finite size.
- **Thm 2.6 and Thm 2.2 are complementary, not ordered.** Thm 2.6 wins iff
  `2C₁² + 3C₂² < 24C₁C₂` (comparable terms); Thm 2.2 wins on lopsided ones
  (`C₁=10⁶, C₂=1`: `20.93` vs `12.76`). Use the min.

---

## 3. THE MISSING THEOREM

**Short version: it is a constant, not a growth rate.** The exponent of the
separation bound is settled at `r−1`, matching the achieved rate. What is left
is the effectiveness of one input-dependent constant.

### 3.1 What is already in hand

**Theorem 2.9 (repo's own; `PROOFS.md` §2.9).** Fix distinct squarefree
`d_1, …, d_r` and set `alpha_i = sqrt(d_i/d_r)`, so

```
L(C_1,…,C_r)  =  C_r + Σ_{i<r} C_i·alpha_i  =  S / (C_r·sqrt(d_r)).
```

By Cor 1.3, `1, alpha_1, …, alpha_{r−1}` are ℚ-linearly independent. Apply
the Subspace Theorem in `r` variables to the `r` forms

```
L_0 = L,        L_1 = C_1,   L_2 = C_2,   …,   L_{r−1} = C_{r−1}.
```

Exactly `r` forms in `r` variables — the count the theorem requires; note
`C_r` is **not** among them and is not needed. Since
`Π_{i<r}|C_i| ≤ ‖C‖^{r−1}`, every `C` outside a finite union of proper rational
subspaces satisfies `|L(C)| > ‖C‖^{−(r−1+δ)}`, hence

```
E  ≤  (r−1+δ)·m  +  log₂(1/c)  +  O(1)  =  ((r−1+δ)/r)·N  +  O_{d,δ}(1),
```

with `m = log₂‖C‖_∞` and `N = r·m` for fixed radicands.

**The exponent `r−1` is the whole content, and it is sharp.** A single linear
form in `r` variables has Subspace-Theorem exponent `r−1`; EHS's
simultaneous application pays `2n` instead, which is why Thm 2.9 is
`1.5×` better than theirs at `r = 3` and better at every rank. Measured
corroboration from the achieved side, `E / log₂‖C‖_∞`:

| | `r = 2` (coefficient-Pell) | `r = 3` (sorted gap) |
|---|---|---|
| measured | `1.434, 1.222, 1.149, 1.112, 1.090` → `1` | `2.312, 1.896, 2.098, 1.999` → `2` |
| `r−1` | `1` | `2` |

**Consistency check the repo owns outright.** At `r = 2` the forms are
`L_0 = C_2 + C_1·alpha_1`, `L_1 = C_1` — two forms, two variables — which
*is* Roth, and the exponent is `1 = r−1`. Thm 2.2/2.6 prove the `r = 2` case
sharp at `1/2` **with an explicit constant**, attained on the coefficient-Pell
family. So the single-form application is known-correct at `r = 2`, which is
the strongest available evidence that the `r ≥ 3` version is the right
statement rather than a formal artefact.

### 3.2 The missing theorem, stated

> **Open Problem 1 — the constant.** Bound
> `c = c(d_1,…,d_r,δ)` in Thm 2.9 from below so that
>
> ```
> log₂(1/c)  ≤  poly(N),        N = total input bit size.
> ```
>
> This is the same object as EHS's `γ` (arXiv:2312.02057), reached by a
> different route. Thm 2.9 improves the exponent and leaves the constant
> untouched.

> **Open Problem 2 — the exceptional subspaces.** The finitely many proper
> rational subspaces `V_1, …, V_t` of Thm 2.9 are not known explicitly. This
> is the same inefficacy, in a second guise.
>
> Unlike Open Problem 1, this one is **detectable**: `C` lies in a proper
> rational subspace iff some nonzero integer vector `a` satisfies `a·C = 0`,
> decidable by exact linear algebra over ℚ. So the generic case already has a
> **constant-free** bound:

> | condition on `C` | what you get |
> |---|---|
> | `C ∉ V_1 ∪ … ∪ V_t` (generic) | `|S| > ‖C‖^{−(r−1+δ)}` — **fully effective, no constant** |
> | `C ∈ V_1 ∪ … ∪ V_t` (exceptional) | codimension ≥ 1; detectable but not excluded |

> An adversary reaches the second row only by satisfying a codimension-1 linear
> equation on `C`, which *constrains* rather than *enables* him — the
> near-cancellation that makes `S` small is a separate geometric condition.
> That the two are generically independent is a **heuristic, not a theorem**,
> and is labelled as such.

**Why this is still not a proof that SRS ∈ P.** Deciding `sign(S)` by
approximation needs a *known* lower bound to set the precision threshold.
Thm 2.9 supplies `c(d,δ)`, which is unknown. The classical norm approach
needs `O(2^r·log)` bits of certified precision and Thm 2.9 would need only
`O(N)` — but with a constant that cannot be computed. EHS Thm 4 escapes the
constant and is per-instance computable, at the cost of a checkable
hypothesis: if `Q ≥ (2n‖x‖_∞)^{3/2}` and the first successive minimum of an
explicit lattice satisfies `λ₁ ≥ Q^{−(1+1/3n)}`, then `|E| ≥ Q^{−(n+1)}`. `λ₁`
is computable by lattice reduction, so this is an *algorithm* rather than a
hope — it is the most concrete implementable target in the repo.

### 3.3 What would falsify the sharp form

The asymptotic conjecture `limsup_{N→∞} E/N = (r−1)/r` is falsifiable, and
**`refutation.py` hunts for it actively.** A single instance with a certified
sign, a certified magnitude, and `E/N` above the ceiling — clearing the
`O(1/N)` slack that the pigeonhole prefactor contributes — would refute Thm
2.5 or Thm 2.9.

Results so far, `r = 2..5`, random and greedy strategies plus the F6 ceiling
family: **no refutation found.** Reported as a negative result. It bounds what
was searched and proves nothing.

The correction that made the harness honest is worth recording, because the
same mistake is easy to make in writing: an early version compared `E/N`
against the bare limit `(r−1)/r` and immediately reported **F6 — which *is*
Thm 2.5's own construction — as a refutation**. The false positive came from
ignoring the construction's `Σ√d_i` prefactor, which contributes an `O(1/N)`
rate term. The test `test_f6_is_not_a_refutation` is the regression pin.

### 3.4 The most concrete lead

Every argument in `PROOFS.md` §2.1–§2.7 uses only the orbit of `S` and the
size of its elements — never the minimal polynomials `x² − d_i`, never the
interaction between the radicals, beyond Cor 1.3's independence. EHS note
(§4) that the quantitative Subspace Theorem fails to help *even when every
`β_i` is a square root of an integer*. Read the other way: the square-root
structure is where an effective constant would have to come from.

Thm 2.9 sharpens the exponent using **only** Cor 1.3, so it does not touch
this. The lead is untouched by anything in the repo, and it is the honest
place to put effort.

## 4. Every route tried, and which are closed

Routes 1–5 are closed **as norm methods**. Routes 6 and 7 are **not** closed;
route 7 is the sharpest bound available anywhere.

| # | route | outcome | why |
|---|---|---|---|
| 1 | **Liouville / norm** | `E ≤ (2^k−1)(3n/2 + log k)` | pays the `2^k`-orbit |
| 2 | **Orbit second moment** | Thm 2.3, `≈2^{r−2}` | exact `2^rM`, halves the exponent once |
| 3 | **Single-conjugate pairing** | dead end | `M_j ∈ K_{r−1}`; gain `3.54×` at `r=3` decaying to `1.05×` at `r=8`, i.e. `1+O(2^{−r})` — converges *to* Thm 2.3, never beats it |
| 4 | **Global-flip pairing** | **Thm 2.6, `≈2^{r−3}`** | the norm optimum; no further pairing exists in this class |
| 5 | **Quadratic tower recursion** | dead end | monomial count is quadratic, but the rank saving is cancelled by the squaring that creates the shallower subfield (§6.2) |
| 6 | **Subspace Theorem, simultaneous** | works, but not the best | EHS Thm 8: exponent `2n`, so `E/N → 2`. **An earlier version of this file wrongly declared this route closed — see §6.1 row 5.** |
| 7 | **Subspace Theorem, single form** | **BEST KNOWN — matches the ceiling** | **Thm 2.9**: exponent `r−1`, so `E/N ≤ (r−1+δ)/r`, equal to the achieved rate. The repo's own application. |

**The single most important distinction in this file: route 3 vs route 4.**
An earlier version of `PROOFS.md` §2.5 asserted the `±S` pairing "stops at
`r = 2`". **That was false**, and the reason it was believed is instructive:
the *single-conjugate* pairing does stop (route 3), and that is a different,
worse pairing. The **global flip** `τ: sqrt(d_i) ↦ −sqrt(d_i)` is not a single
conjugate — it is a symmetry of the whole orbit — and it works at every `r`.
Conflating the two cost a factor of `2` in the exponent.

**The honest summary of routes 1–5:** *within the norm*, the barrier is the
field rank `2^r`, and it is not removable by monomial-count economy (route 5),
nor by a shallower subfield (route 5 again), nor by a cleverer pairing (route 4
is the optimum of that class). Every method passing through the norm of the
whole orbit — Liouville, Thm 2.3, Thm 2.6, resultants, the quadratic tower,
and the published algorithms of Cheng and Qian–Wang — lands on `Θ(2^r)`.

**And route 6 leaves that class entirely.** So `Θ(2^r)` is a barrier of the
*method*, not of the *problem*. That is the correction
`STATE_OF_THE_ART.md` exists to record, and the reason §3 is about a constant
rather than about a growth rate.

---

## 5. What a proof would have to contain

Stated as a target, so it can be recognised when it appears:

1. **Not** a norm over the whole orbit. That is `Θ(2^r)` and is final as a
   statement about the norm (Thm 2.6 + §4) — but it is not required, and Thm
   2.9 already reaches the ceiling without ever taking it.
2. **Not** an `ε`-improvement for one algebraic linear form. That exists and
   is applied — it is Thm 2.9, exponent `r−1`, and at `r = 2` it is literally
   Roth. See §6.1 row 5 for the error of believing otherwise.
3. It **would** be an **effective constant**: `log₂(1/c(d,δ)) ≤ poly(N)` in
   Thm 2.9, or the equivalent `γ` in EHS Thm 8, i.e. a bound on the number
   of `q` with `q^{1+δ} Π_i dist_Z(q√a_i) < 1`. §3.2, Open Problem 1.
4. **And** the explicit description of Thm 2.9's exceptional rational
   subspaces — detectable but not yet known. §3.2, Open Problem 2.
5. **Or** a per-instance effective route: EHS Thm 4's lattice certificate
   (`λ₁ ≥ Q^{−(1+1/3n)}` on an explicit lattice, computable by reduction).
   **This was the most promising lead and it does not survive contact with the
   paper.** `Q` is a free parameter, and the object that is *not* computable is
   the threshold `Q₀`; making it explicit via
   `dist_Z(q√a) ≥ 1/(2q√a+1)` costs `exp(Θ(n²a))`, so the honest rate is
   `E = O(n·N)`, not `O(N)`. The `λ₁`-computable/threshold-ineffective
   conflation is the trap. See `STATE_OF_THE_ART.md`.
6. **Plausibly**, something exploiting that every `β_i` is a **square root of a
   rational** — §3.4. This is the most concrete lead available, and the one the
   literature has flagged without solving.
7. **Or** a counterexample. `refutation.py` hunts for one; an instance above
   the ceiling would refute the sharp form outright. §3.3.

---

## 6. Errors made and retracted

Recorded because each was plausible, and each is now pinned by a test asserting
the *correct* direction, so it cannot return silently.

**6.1 Eight substantive retractions**

| # | the error | the truth | pin |
|---|---|---|---|
| 1 | Thm 2.3 "attained to 3 decimals" on the Pell family | slack *grows*; the norm step spends the `±S` pairing. Thm 2.6 spends it properly and *is* tight | `test_bound_attained_on_coefficient_pell_is_false`, `TestGlobalFlipPairing::test_attained_on_coefficient_pell` |
| 2 | The `±S` pairing "stops at `r = 2`" | the global flip works at every `r`, uniform factor `≈2` | `TestGlobalFlipPairing::test_retracted_claim_that_pairing_stops_at_r2` |
| 3 | Thm 2.6 "dominates Thm 2.2 at `r = 2`" | complementary, not ordered; crossover `2C₁²+3C₂² = 24C₁C₂`. Checked only on the balanced family | `TestGlobalFlipPairing::test_complementary_to_two_term_bound_at_rank_two` |
| 4 | Tower support rank is `j−1` | rank is `max(j−2, 0)`. The norm-barrier conclusion survives, for a *different* reason (§6.2) | `TestTowerStructure::test_support_rank_is_j_minus_two_not_j_minus_one`, `::test_r3_relative_norm_is_two_term` |
| 5 | The Subspace route is closed; the missing theorem is an `ε`-improvement for one algebraic form in `r ≥ 3` variables | **The route works.** Schmidt's Thm 1B applies — the `n−1` coordinate forms are supplied by the reduction, only `L_0` must come from the problem — and EHS get a **polynomial** bound. The error was a form count, invisible at `r = 2` | `TestSubspaceTheoremRoute::test_ehs_bound_is_polynomial_in_the_number_of_terms`, `::test_ehs_beats_thm_2_6_by_an_exponential_factor`, `::test_overclaim_is_invisible_exactly_at_r2` |
| 6 | "`O(2^r)` is the best known guarantee" | false. `2^r` is the optimum **of the norm**; the literature's bound is polynomial in `n`. The confusion was treating a property of the *method* as a property of the *problem* | `TestSubspaceTheoremRoute::test_ehs_beats_thm_2_6_by_an_exponential_factor` |
| 7 | The EHS bound `γ·(n‖x‖)^{−2n}` is the sharpest available | **false, and the repo was leaving a factor on the table.** A *single-form* application of the same theorem gives exponent `r−1`, matching the achieved rate — Thm 2.9. EHS's simultaneous application pays `2n` | `TestSingleFormSubspaceTheorem` (8 tests) |
| 8 | `refutation.py` reported F6 as a refutation of the sharp conjecture | **false positive, and the worst kind.** F6 *is* Thm 2.5's own pigeonhole construction, so it must sit at the ceiling. The harness compared a finite instance against the *limit* `(r−1)/r` and ignored the `Σ√d_i` prefactor's `O(1/N)` term | `tests/test_refutation.py::test_f6_is_not_a_refutation` |

**Two of these are my own errors found by tests I wrote in the same pass**
(rows 7 and 8). Row 7 is a missed improvement, not a false claim; row 8 was a
false "REFUTES" verdict printed by a harness whose entire purpose is to be
trusted on that word.

**6.2 Why retraction 4 does not break the norm barrier** — the subtlety worth
keeping

`N_j = S_{j−1}² − C_j²d_j` genuinely lives in a *shallower* subfield than
`K_{j−1}`: its support rank is `max(j−2, 0)`, verified exactly over 300+
random F₂-independent sets per level for `j = 2..7`. The `r = 3` instance is
the identity

```
S · σ₃(S)  =  (C₁²d₁ + C₂²d₂ − C₃²d₃)  +  2C₁C₂ · sqrt(d₁d₂)
```

— a **two-term** expression in a quadratic field, rank `1` not `2`. So the
rank saving is real, and it is a factor `2` in the exponent.

But `N_j` is built from `S_{j−1}²`, so its second moment is `M' ≈ 4M`, i.e.
`log₂M' ≈ 2log₂M` instead of `≈ log₂M`. And the size enters the guarantee
linearly. So the factor `2` from the rank and the factor `2` from the squaring
**cancel exactly**:

```
tower(rank j−2) / Thm 2.6   at r = 3..8:
1.0049, 0.9982, 1.0007, 0.9975, 1.0019, 1.0003        →  1 ± 0.005
```

No gain, no loss **within the norm class**. The old *reason* was wrong; the
*conclusion* holds for the norm, and the correct reason is now written down.
(It does **not** survive as a statement about the problem — route 6 bypasses
the norm entirely.)

**6.3 Two smaller corrections**

- The "worst two-term constant" was long recorded as the open interval
  `[1/4, 9/4]`. Both numbers were wrong. The constant is exactly `1/2`, proved
  from above (Thm 2.6) and achieved from below (Thm 2.2) — and it is this
  repo's result, not the literature's.
- F1 reaches `1/4`, not `1/2`. The old claim that F1 was the worst two-term
  family was wrong: the hard arithmetic can live in the *coefficients* with
  tiny fixed radicands (`2, 3`), a strictly easier family to construct than F1,
  which hides it in `10^6`-bit radicands.

**6.4 The methodological lessons**

*Lesson 1 (retraction 5): check negative claims about the literature as
carefully as positive ones.* I asserted "this route is closed" and "that case
is open" without checking either, while verifying positive claims
attentively. To say a theorem does not apply, one must read its actual
hypotheses — and here they *were* satisfied, because the reduction supplies
forms I had not counted. The correction came from searching the literature,
not from reasoning harder about the same reduction.

*Lesson 2 (retractions 7 and 8): a result you have quoted is not automatically
the best result, and a threshold you have written is not automatically
correct.* Row 7: having found the Subspace Theorem applied, I stopped there
and did not ask whether a different choice of forms does better. The
single-form application is `1.5×` sharper at `r = 3` and better at every
rank, and it is *the same theorem*. Row 8: I wrote a harness whose job is to
say "REFUTES" or not, and its very first output was a false "REFUTES", because
I compared a finite instance against an asymptotic limit and forgot the
prefactor. Both errors were found by tests written after the fact — which is
the only reason they are retractions rather than published claims.

*Lesson 3 (general).* Every one of the eight retractions in this section was
plausible on first reading, and five of them were caught only by an explicit
check — a test, or an arithmetic verification before writing the claim down.
The pattern that produced all eight is the same: stating a rate, an exponent,
or a comparison without first computing it at the size in question. The
discipline that catches them is dull and non-negotiable: **compute the number
before you write the sentence that contains it.**

---

## 7. Provenance and honesty notes

- Every theorem above is proved in `PROOFS.md`, except the EHS results, which
  are **quoted** from arXiv:2312.02057 and labelled as such wherever they
  appear. **Thm 2.9 is this repo's own** and is proved in `PROOFS.md` §2.9.
  Every retraction above is pinned by a named test.
- Computation in this repo is **corroboration, not evidence**. The measured
  tables (`3.000, 2.333, …`; the slack decay; the `1 ± 0.005` ratio; the
  `E/log₂‖C‖ → r−1` measurements; the `0.68` F6 rate against its `0.11`
  slack) confirm that code and closed forms agree. The evidence for the
  theorems is the written proof.
- **Literature check performed.** Roth (1955), Schmidt's Subspace Theorem and
  its Theorem 1B corollary were checked. The decisive reference —
  Eisenbrand–Haeberle–Singer, arXiv:2312.02057 (SoCG 2024) — was read in full.
  Davenport–Schmidt (1967) was checked and **excluded**: it concerns
  approximating a real number by quadratic irrationals, not the
  `ε`-improvement at issue. See `STATE_OF_THE_ART.md` §6 for exactly what was
  and was not verified.
- **A live negative result.** `refutation.py` reports no counterexample to
  the *asymptotic* rate over the ranks and strategies searched. This bounds
  what was searched and **proves nothing**. It is recorded because a
  reproducible negative is worth more than an unrecorded shrug, not because it
  supports the conjecture.
- **Attribution.**
  - Thm 2.5 is EHS Thm 10, rediscovered independently.
  - Thm 2.9 is the repo's own single-form application. It is *not* in the EHS
    paper, which uses the simultaneous form; at `r = 2` it reduces to Roth,
    where the repo's Thm 2.2/2.6 independently confirm it with an explicit
    constant.
  - The sharp `k = 2` constant `1/2` (Thm 2.2 + Thm 2.6) is this repo's, and
    is **stronger** than EHS at that rank, which does not use the radicands.
  - Thm 2.6, Prop 2.4′, and F6 have no counterpart in the reference.
- `2^r` is the F₂-rank parameter, **not** the number of terms `k`. A large `k`
  over a small `r` buys an adversary nothing: F1, the coefficient-Pell family,
  and F6 all sit at small `r`. Any hard instance must make the `d_i`
  F₂-independent, and `refutation.py` enforces that by rejection.
- Complexity context, for calibration: the zero test is in P (Thm 1, matching
  Blömer); the sign problem's best known upper bound is PSPACE; it is not
  known to be in NP; and NP would follow from a polynomial separation bound.

---

## 8. Summary

**The exponent of the separation bound is settled; the constant is not.**

Schmidt's Subspace Theorem applies to SRS as a *single* linear form, and
`PROOFS.md` Thm 2.9 does exactly that: exponent `r−1` on `‖C‖_∞`, hence
`E/N ≤ (r−1+δ)/r`. The pigeonhole construction of Thm 2.5 supplies instances
at `E/N ≥ (r−1)/r − o(1)`. **The guaranteed and achieved rates coincide in
their linear part, at every rank**, so the sharp form is the *asymptotic*
statement `limsup_{N→∞} E/N = (r−1)/r` and there is no growth-rate gap left to
close. (The `sup` version of this claim is retracted — see `PROOFS.md`
§2.10.2, where a certified `r = 3` instance gives `E/N = 1.2795 > 2/3`.) At `r = 2` this is Roth, and the
repo's Thm 2.2/2.6 confirm it with an explicit constant attained on the
coefficient-Pell family — so the application is known-correct at the one rank
where the answer is owned.

This supersedes the earlier framing in two ways. The `2^r` in Thm 2.6 is a
property of the **norm**, not of the problem; and EHS's own bound, which this
file previously treated as the best available, is `1.5×` looser at `r = 3`
because it applies the same theorem in a more expensive shape.

Two open problems remain, both consequences of the Subspace Theorem's
inefficacy: **bound the constant** `c(d,δ) ≤ 2^{−poly(N)}`, and **describe the
exceptional rational subspaces** — the second being detectable by exact linear
algebra, so the generic case already carries a constant-free bound. The most
concrete implementable target is EHS Thm 4's per-instance lattice certificate,
which escapes the constant at the cost of a checkable hypothesis. And
`refutation.py` hunts actively for a counterexample to the sharp form, with
certified signs and magnitudes; it has found none, which bounds what was
searched and proves nothing.

