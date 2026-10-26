# Paper 36 — M7.12: the engine, and a structural no-go

Version 1.0.  Builds `m712_engine.py` (four layers), records the layer-3
degeneracy as a checked fact, corrects a formula in `Docs/35`, and replaces the
"factor `9n(n−1)` deficit" of `Docs/35` Thm 35.4 with a **structural** statement:
the required improvement is impossible for *every* integrality or
product-of-conjugates argument, for **every** `n ≥ 2`.

---

## 36.1 Correction to `Docs/35` Thm 35.1

`Docs/35` printed the barrier as `|ρ|^{3k}/(2β_j³)`.  That is **wrong**:

```
|ρ| = |ε'/ε| = 0.1715728753        ε^{−1} = |ε'| = 0.4142135624
```

so `|ρ|^{3k} ≠ ε^{−3k}`.  The norm argument divides by the surviving conjugate
product `~ 2β_j³ε^{3k}`, so the correct barrier is

> **`|m − β_jε^k|  ≥  ε^{−3k} / (2β_j³)  =  |ε'|^{3k}/(2β_j³)  =  exp(−3k·log ε + O(1))`.**

The prose of `Docs/35`, its gap arithmetic (`3 log ε = 2.6441`), and its
numerical verification all used `log ε` and are unaffected.  Only the displayed
`|ρ|` was wrong, and it erred in the **optimistic** direction — `ε^{−3k}` is a
*stronger* bound than the erroneous `|ρ|^{3k}`.  `Docs/35` corrected in place.

Exact re-verification of the corrected bound, `k = 1…60`: slack
`d_k·2β_j³·ε^{3k} ≥ 1` everywhere, minimum `3.0916702` at `k=1`.

## 36.1b A second correction: `N_L(y_k)` is a rational, not an integer

`Docs/35`'s proof opened "`y := m − β_jε^k ∈ L`, so `N_{L/Q}(y) ∈ Z`".  **False.**
`y_k` is not an algebraic integer; the integer is the norm of the *scaled*
element `2a₁y_k = 2a₁m_k − P_k√(a₁a_j) − a₁Q_k√a_j`.  Concretely
`N(y_1) = 228/256` for `a₁=2, a_j=3`.

So `|N(y_k)| ≥ (2a₁)^{−4}`, not `≥ 1`.  This changes **only the leading
constant** of Theorem 35.1 (by `(2a₁)^4 = 256` at `a₁=2`), never the exponent
`3`, so every conclusion of `Docs/35` and of §36.5 below stands unchanged.
`Docs/35` corrected in place.  (Theorem 36.1 is unaffected: `c_k ∈ Z[√a₁]` *is*
an algebraic integer, so `N(c_k) ∈ Z` genuinely.)

## 36.2 The engine

`m712_engine.py`.  Layer 1 recurrence, layer 2 norm, layer 3 invariant, layer 4
near-integrality.  Design points:

* **Index safety.**  `eps_power` builds `ε^k = P_k + Q_k√a₁` by the exact
  recurrence `(P,Q) ↦ (P + a₁Q, P + Q)`, so `k` *is* the exponent by
  construction.  This retires the `Docs/34` §34.2 hazard by construction rather
  than by vigilance.  Checked: `P_k² − 2Q_k² = (−1)^k` for `k = 0…6`.
* **Exactness.**  `nearest_integer` brackets `z_k = 2a₁x_k` with
  `signexact.certified_bracket` (integer arithmetic, doubling loop) and confirms
  the candidate with a second bracket; undecided ⇒ `UndecidedError`, never a
  guess.  `2a₁x_k = a₁Q_k√a_j + P_k√(a₁a_j)`, derived in the docstring.
* **Split enforced by test.**  All log/float arithmetic lives in `__main__`;
  `tests/test_no_float.py` gains `m712_engine.py` as a decision-path module, so
  the library core cannot regress into float.

## 36.3 Layer 3 is degenerate — recorded, not silently dropped

`Q(u,v) = v² − 2uv − u²` satisfies `Q(u_{k+1},u_{k+2}) = −Q(u_k,u_{k+1})`: it is
**anti**-invariant, not invariant.  More decisively it factorises exactly,

```
Q(u,v)  =  (v − εu)(v − ε'u)          (ε + ε' = 2,  εε' = −1)
```

and for a pure eigen-solution `x_k = β_jε^k` the first factor vanishes
identically.  So

> **`I_k = x_{k+1}² − 2x_kx_{k+1} − x_k²  =  0`  exactly, for every `k`.**

Verified by the engine (`invariant_row` true for `k = 1…39`) and the
factorisation by 200/200 randomised exact `Elem` pairs.  The proposed
invariant layer therefore carries **no** information here; the surviving content
is layer 4.  Sign correction found while testing: with `y_k := m_k − x_k`,
`c_k = y_{k+1} − εy_k` (**not** `+`); verified at `k=10`,
`0.296351 − 2.414214·0.183089 = −0.145663 = −c_{10}`.

## 36.4 Layer 4, and the best bound it gives

`c_k = m_{k+1} − εm_k = (m_{k+1} − m_k) − m_k√a₁ ∈ Z[√a₁]`, so

> **Theorem 36.1 (pair / lattice bound).**
> `d_{k+1} + ε·d_k  ≥  |c_k|  ≥  1/(2√a₁·m_k)  ~  ε^{−k}/(2√2·β_j)`.

*Proof.*  `c_k` is a nonzero element of `Z[√a₁]`; its conjugate is
`c_k − 2√a₁m_k`, of modulus `≤ 2√a₁m_k + |c_k|`.  `N(c_k) = c_k·conj(c_k)` is a
nonzero integer, so `|c_k| ≥ 1/(2√a₁m_k + |c_k|)`, and `|c_k| = |y_{k+1} − εy_k|
≤ d_{k+1} + εd_k`. ∎

`N(c_k) ≠ 0` verified for `k = 1…59`.  The bound is exact: slack
`(d_{k+1}+εd_k)·2√a₁m_k ≥ 1` for `k = 1…59`, minimum `4.4852814` at `k=1`.

Note this is **stronger** than `Docs/35` Thm 35.1 (exponent `1` vs `3`), because
it lives in the degree-2 subfield rather than in `L`.

## 36.5 Theorem 36.2 — the STRUCTURAL no-go

> **Theorem 36.2.**  *No* argument that (i) forms an algebraic number out of the
> integers `m_k` (and the fixed data), and (ii) bounds it by the product of its
> other conjugates, can establish M7.12 — **for any `n ≥ 2`**.

*Proof.*  Such a bound has the shape `|α| ≥ 1/∏_{σ≠σ₀}|σ(α)|`, and its exponent
is governed by the size of the *surviving* conjugates.  There are exactly two
cases.

* **Single index** (`k` alone).  The only non-integer available at one index is
  `x_k = β_jε^k`, and by `Docs/35` Thm 35.2 (`β_j ∉ Q(√a₁)`, equivalently
  `a₁a_j` not a rational square) it has degree 4, not 2.  Its three surviving
  conjugates are *forced*: `m_k + β_jε^k ≈ 2m_k`, `m_k ± β_jε^{−k} ≈ m_k`, all
  `~ ε^k`.  Hence the exponent is **3** — this is `Docs/35` Thm 35.1, and no
  single-index object can do better, because degree 2 is unavailable at one
  index.
* **Two consecutive indices.**  Degree 2 *is* available, but only as
  `c_k = m_{k+1} − εm_k`, and its surviving conjugate is
  `≈ 2√a₁m_k ≈ ε^k`.  Exponent **1** — Theorem 36.1.  But the conclusion is a
  bound on the **sum** `d_{k+1} + εd_k`, not on `d_k`.

So the best exponent obtainable on `d_k` *alone* is `3`, and the best obtainable
anywhere is `1`, attained only on a sum.  M7.12 asks to exclude a **single**
`d_k < θ₀^k` with `θ₀ = ε^{−1/(3n(n−1))}`, i.e. it needs exponent

```
1/(3n(n−1))   <   1      for every n ≥ 2        ( = 1/6, 1/18, 1/30, … )
```

A *smaller* exponent is a *stronger* bound, so the requirement is strictly
beyond every value the toolkit produces. ∎

> **Corollary 36.3.**  M7.12 is not reachable by norms, lattices, heights,
> Liouville, Ridout, or the Subspace Theorem.  `Docs/34` Thm 34.1 killed the
> height/Liouville family and `Docs/35` §35.6 killed the power-law family; this
> theorem kills the *remaining* family — pure product-of-conjugates — with a
> single argument valid for all `n`.  What M7.12 needs is not a better
> inequality but a genuinely different kind of input.

## 36.6 Measured confirmation, and a negative result on the diagnostic

The engine's purpose was to find the quantity that grows when `d_k` falls.  It
does not exist in this setting, and the reason is quantitative:

> **Conjugate amplification is exactly flat.**
> `A_k = |N(2a₁y_k)|/(|y_k|ε^{3k}) → 2β_j³(2a₁)^4`, measured ratio `1.00000`
> (`117.57551` vs predicted `117.57551`).

So the three crude conjugate estimates in the norm proof are asymptotically
**tight**: there is no slack in the norm estimate to convert into a better
exponent.  This is the direct measurement refuting the conjugate-amplification
strategy.  Consistently, no layer quantity correlates with `−log d_k`:

```
corr( −log d_k , G_norm      ) =  0.147
corr( −log d_k , log A_conj  ) =  0.078
corr( −log d_k , log|N(c_k)|/k) = −0.205
```

and empirically `d_k` never decays: `min_k d_k = 0.0036793` at `k=43`, with no
exception below `θ₀^k` beyond `k=16` (the `k ≤ 16` "exceptions" are the
artifact that `θ₀^k > 1/2` while `k` is small).

**A correction to my own speculation.**  I suggested Theorem 36.1 might isolate
exceptions.  It does not, and the engine says so: adjacent exception pairs
occur (`(1,2), (2,3), …, (6,7), (9,10)`, minimum gap `1`).  The reason is that
Theorem 36.1's `ε^{−k}` sits *below* the threshold `ε^{−k/(3n(n−1))}` for `n ≥ 2`,
so it constrains a regime deeper than the one M7.12 asks about.

## 36.7 Where M7.12 now stands

| ingredient | status |
|---|---|
| best single-index bound | `d_k ≥ ε^{−3k}/(2β_j³(2a₁)^4)` (`Docs/35` Thm 35.1, corrected twice) |
| best two-index bound | `d_{k+1}+εd_k ≥ ε^{−k}/(2√2β_j)` (Thm 36.1, **new**) |
| required exponent | `1/(3n(n−1)) < 1` for all `n ≥ 2` |
| gap | **structural**, not a tuning deficit (Thm 36.2) |
| heights / Liouville | refuted (`Docs/34` Thm 34.1) |
| Ridout / Subspace-Theorem power laws | strictly worse (`Docs/35` §35.6) |
| conjugate amplification | measured flat (§36.6) |
| invariant layer | degenerate (§36.3) |
| empirical behaviour | no exception for `k > 16`; exception set plausibly empty |

**Final M7.12.**  Prove that for distinct squarefree `a₁ ≠ a_j` the set
`{ k : ‖(√a_j/2√a₁)(1+√a₁)^k‖ ≤ (1+√a₁)^{−k/(3n(n−1))} }` is finite with
`k*·log ε ≤ poly(N)`.  **Blocked**, and now for a stated structural reason: every
integrality- and norm-based method saturates at exponent `≥ 1`, while the
requirement is an exponent `< 1`.  I am not claiming a solution, and I am not
asserting that no theorem reaches it — only that the algebraic toolkit
represented here provably does not, and that the missing input is of a different
kind.