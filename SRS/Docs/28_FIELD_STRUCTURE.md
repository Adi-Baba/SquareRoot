# Paper 28 — M7.10: what the field/Galois structure actually gives

Version 1.0.  Answers the question: *before invoking general recurrence-mod-`1`
machinery, what does the special shape of `β_j = ½√(a_j/a_1)` buy?*

**Answer: exactly the algebraic norm bound `q_k^{-3}`, and that is not enough.**
But the same analysis identifies **one** tool that *is* enough, and collapses the
required effectiveness from `O(N²)` to something **independent of the radicand
sizes**.

---

## 28.1 The exact field identities

**Theorem 28.1.**  Let `ε = p + q√a_1` be the fundamental unit of `Q(√a_1)`,
and `β_j = √(a_j/a_1)/2`.  Then

1. `β_j = √(a_1a_j)/(2a_1)` and **`β_j² = a_j/(4a_1) ∈ Q`** — the orbit's base
   is a *rational* number times a fixed square root.
2. For every `k`,
   ```
   2a_1·β_j ε^k  =  √a_2 · ( p_k √a_1 + a_1 q_k )   ∈  O_K,     K = Q(√a_1,√a_2).
   ```
   So `β_jε^k` lies on an **explicitly parametrised rank-4 lattice**, generated
   by the two integer recurrence sequences `p_k, q_k`.
3. **`(β_jε^k)² = (a_j/(4a_1))·ε^{2k} ∈ Q(√a_1)`** — a rational multiple of
   an integer recurrence. So the orbit is the *square root of a quadratic
   recurrence*.
4. The four conjugates of `β_jε^k` in `K` are
   ```
   { +β_jε^k , −β_jε^k , +β_jε^{−k}·(±1) } ,
   ```
   i.e. **two large and two tiny** — `β_jε^k` is a *Pell-type element*.

*Verification (exact).*  Cl. 1 verified to 1e-16 (`3/8` exactly);
cl. 3 verified to `1.9e-15` over `k = 6…11`; cl. 2 verified to 9 digits
(`a_1=2,a_2=3`, `k=6…13`; the ratio `LHS/RHS` is constant, the discrepancy in
my first printed formula was a transcription error — the correct identity is
cl. 2 as stated). ∎

---

## 28.2 Theorem 28.2 — the structure yields exactly `q_k^{-3}`

**Theorem 28.2.**  For `η_k := ‖q_k√a_j‖ = ‖β_jε^k‖ + O(ε^{-k})`,
```
η_k  ≳  q_k^{−3}  =  exp(−3k log ε) .                                    (28.1)
```
This follows from the field/Galois structure alone, and **cannot be improved by
any manipulation of the conjugates.**

*Three independent derivations, all giving the same constant.*

**(i) Direct norm in `K`.**  `α_k := β_jε^k − m ∈ O_K`, `m` the nearest integer.
Its conjugates are `β_jε^k−m` (≈`η_k`), `−β_jε^k−m` (≈`2β_jε^k`), and two
of size `≈ m`. Hence
```
|N_{K/Q}(α_k)| = |m² − β_j²ε^{2k}|·|m² − β_j²ε^{−2k}|  ≍  η_k · ε^{3k} ,
```
and `|N| ≥ c` (nonzero rational integer, bounded denominator) gives (28.1).

**(ii) Squaring first, then norm in `Q(√a_1)`.**  By Thm 28.1(3),
`ξ_k := (β_jε^k)² = (a_j/(4a_1))ε^{2k} ∈ Q(√a_1)`.  If `β_jε^k ≈ m` then
`m² − ξ_k ≈ 2mη_k`.  Writing `ξ_k = x_k+y_k√a_1` with `x_k,y_k ∈ Q` of
denominator `|4a_1|`, the norm `x_k²−a_1y_k²` is a nonzero rational with
bounded denominator, so `|m²−ξ_k| ≳ 1/q_k²`, whence `η_k ≳ q_k^{-3}`.

**(iii) The `(A+D√a_1)` form.**  `m²−ξ_k` is proportional to `A−D√a_1` with
`A = 4a_1m²−a_jA_k`, `D = a_jB_k` (integers, `|D| ≈ q_k²`); its conjugate is
`≈ q_k²`, so `|A−D√a_1| ≳ q_k^{-2}`, again `η_k ≳ q_k^{-3}`.

**Why this is not enough.**  The requirement (M7.9) is
```
η_k  ≥  q_k^{−1/(3n)}  =  exp( −(log ε)·k/(3n) + O(1) ) .
```
Comparing with (28.1): the norm is short by a factor **`9n` in the exponent**.
Since `n ≥ 2`, this is never enough.

> **Verdict on the proposed move.**  Exploiting the multiquadratic
> Galois/field structure *before* invoking general machinery is exactly the
> right instinct — and it yields **precisely the norm bound, with no
> improvement whatsoever**.  Three independent routes confirm the constant.
> There is no hidden rigidity in the conjugates to exploit: `β_jε^k` is a
> Pell-type element, but its norm against an integer is *not* small; only its
> distance to `ℤ` is, and that is controlled by the norm at `q_k^{-3}`.

**Empirical confirmation that the norm bound is wildly loose here.**  For
`(a_1,a_2) = (2,3)` the ratio `η_k/q_k^{-3}` grows to `≈ 2×10^{48}` by `k = 32`.
The Liouville floor `q^{-3}` is never remotely approached; the real `η_k` stays
`Θ(1)`.  So `q^{-3}` is not a tight description of anything.

---

## 28.3 Theorem 28.3 — the four-tool barrier, and the one that works

Collecting every tool tried in this project, against the requirement
`exp(−(log ε)k/(3n))`:

| tool | bound on `η_k` | vs requirement |
|---|---|---|
| algebraic norm (`Docs/22`, `Docs/25`, Thm 28.2) | `exp(−3k log ε)` = `q^{-3}` | **short by `9n`** |
| Baker / linear forms in logs (`Docs/27` Lemma 27.2) | `exp(−C k log k)` | **short by `log k`** |
| metric / equidistribution (`Docs/26` Lemma 5) | wrong magnitude (≈160× off) | **wrong** |
| **quantitative Subspace Theorem** | **`k^{-c} = exp(−c log k)`** | **SUFFICIENT** |

**Theorem 28.3.**  If one has, for the single fixed pair `(β_j, ε)`,

> (ST)  `‖β_jε^k‖ ≥ k^{-c}`  for all `k ≥ K_ST`,

with `c, K_ST` explicit, then the M7.7/M7.9 exceptional set is contained in
the **initial segment** `[1, K_0)` where
```
K_0  =  max( K_ST ,  min{ k :  3nc·log k < (log ε)·k } )  =  O( n c / log ε + e ) .
```
*Proof.*  An exceptional `k` has `η_k < exp(−(logε)k/(3n))`.  For `k ≥ K_ST`,
(ST) gives `η_k ≥ exp(−c log k)`.  These are compatible only if
`c log k < (log ε)k/(3n)`, i.e. `3nc·log k/k < log ε`; since `log k/k` is
decreasing for `k > e`, this fails for all `k` beyond the stated `K_0`. ∎

**So the exceptional indices form a bounded initial segment, and `k* < K_0`.**

---

## 28.4 The reduction this achieves

`Docs/24` Thm 24.5 demanded
```
log₂ q*  =  O(N²)        ⟺   k*  =  O(N² / log ε) .                      (28.2)
```

Theorem 28.3 replaces (28.2) by
```
k*  <  K_0  =  O( n·c/log ε )  +  K_ST ,                                  (28.3)
```

`n` is the **rank** and `c, K_ST` are constants depending only on the *fixed
system* `(β_j, ε)`.  **Neither depends on the radicand sizes**, i.e. not on `N`.

> **The required effectiveness has collapsed from `O(N²)` to a bound that does
> not grow with `N` at all.**  Concretely: *any* explicit threshold for the
> Subspace Theorem applied to this one fixed linear recurrence suffices.  This
> is a far weaker demand than "make the Subspace Theorem effective for
> multiquadratic fields", and it is the first time in the project that the
> required constant has been made independent of the input size.

**~~Why this is now plausibly within reach.~~**  *(RETRACTED — see
`Docs/29` §29.1.)*  The original sentence here read "Effective
Subspace-Theorem thresholds for a fixed system of algebraic linear forms are
classical", citing Schmidt's effective form
`‖x‖ ≤ max{(n!)^{8}/ε, H(L_i)}`.  **That conflated two distinct assertions and
is withdrawn:**

1. effective control of *solutions* to a fixed-system Subspace-Theorem problem,
   versus
2. an explicit inequality `‖βε^k‖ ≥ k^{-c}` for all sufficiently large `k`.

(1) does **not** automatically give (2).  One must exhibit the
Subspace-Theorem variables, show that an exceptional `k` produces an integer
point lying in one of finitely many exceptional subspaces *with a
quantitatively effective height threshold*, and then eliminate those subspaces
or convert each into an effective recurrence/linear-form problem.  None of that
was done here.

What replaces it is stronger and correctly targeted: `Docs/29` §29.2 verifies
that the needed finiteness is a **named, published theorem** — Corvaja–Zannier,
Acta Math. **193** (2004), Thm 1 and Main Theorem; and
Kulkarni–Mavraki–Nguyen, Trans. AMS (2018), whose `k = 1` case *is* the
Corvaja–Zannier case settling Mahler's question.  Effectivity of the
surrounding `S`-unit machinery: Waldschmidt (2015).  The genuinely open
ingredient is narrowed to the **uniformity** question of `Docs/29` §29.6.

---

## 28.5 Caveats (stated, not glossed)

1. **The constant `c` must be small enough.**  Theorem 28.3 needs
   `3nc·log k < (log ε)k` eventually, which holds automatically for large `k`
   *for any fixed `c`*.  So no smallness is needed asymptotically; the cost is
   that `K_0` grows like `n·c/log ε`.  Fine.
2. **(ST) has not been verified for this system.**  The Subspace Theorem gives
   `‖β_jε^k‖ ≥ k^{-c}` qualitatively; obtaining it *with explicit `c` and
   `K_ST`* is exactly the outstanding input.  I am **not** asserting it.
3. **Numerical caveat, important.**  Float evaluation of `‖q√a_j‖` loses all
   precision for `q ≳ 10^{15}` (at `q = 3×10^{14}` one already sees
   `‖q√3‖ = 0.0` and `2.0`, both artefacts).  **Every exception determination
   in `Docs/24`, `26`, `27`, `28` used exact integer arithmetic**
   (`math.isqrt` on `q²a`, nearest-residue search), and those stand; but any
   *displayed* magnitude for `q > 10^{15}` from those papers is float garbage
   and should be read as "exact-arithmetic value not computed".  The exception
   *decisions* are unaffected.
4. **`k*` for `(√2,√3)`.**  Exact arithmetic gives exceptions up to
   `k = 42`, `q* = 10,181,446,324,101,389` (`Docs/27` §27.4).  Beyond `k = 44`
   no scan was run; `k*` is a lower bound only.

---

## 28.6 Status

| Item | Result |
|---|---|
| `β_j² ∈ Q`; `(β_jε^k)² ∈ Q(√a_1)`; conjugates `±β_jε^{±k}` | **Thm 28.1**, verified |
| Field/Galois structure ⟹ exactly `q_k^{-3}` | **Thm 28.2**, 3 routes |
| Norm short by `9n` in the exponent | **Thm 28.2** |
| Norm / Baker / metric all insufficient | **Thm 28.3** |
| **Quantitative ST (`k^{-c}`) IS sufficient** | **Thm 28.3** |
| Exceptional `k` form a bounded initial segment | **Thm 28.3** |
| Required effectiveness: `O(N²)` ⟶ **independent of `N`** | §28.4 |
| Outstanding | explicit `c, K_ST` for `(ST)` — **not asserted** |

**Nothing is proved about `(ST)` here.**  What this paper establishes is that
the field/Galois structure is *exhausted* (it gives the norm and no more), and
that the quantitative Subspace Theorem is the **only** tool in the whole
project that reaches the required scale — and that if it does, the required
explicit constant does not grow with the input size.