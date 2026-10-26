# Paper 23 — Effectiveness phase: what the literature actually delivers

Version 1.0.  A literature-level audit of the effectiveness question left
open by `Docs/22`.  Result: **the polynomial part of SRS is already solved
in the literature; the entire remaining obstruction is a single ineffective
Subspace-Theorem constant — exactly the reduction `Docs/22` reached by
internal reasoning, now confirmed from the outside.**  We also identify a
sharper form of the target than the natural one, and we record a quantitative
cap on the unknown constant derived from our own certified witnesses.

---

## 23.1 Nomenclature map

Our problem is the classical **sum-of-square-roots separation bound**.  With
`n = m`, `a_i = d_i`, `x_i = C_i`, `E = S`:

| this project | literature |
|---|---|
| `E = log₂(1/\|S\|)` | separation bound for `E_sum = Σ x_i√a_i` |
| `N = Σ(bitlen\|C_i\| + bitlen d_i)` | `O(n log n + Σ log a_i + n log M)` |
| `\|C\|∞ = max\|C_i\|` | `‖x‖∞` |
| `c(d,δ)` | `γ` (Eisenbrand–Haeberle–Singer) |

Linear independence of distinct squarefree radicals is **Besicovitch (1940)**
(`Docs/01` Thm A.3); the polynomial-time zero-test is **Blömer (1991, 1998)**
(`Docs/00`, aggregation test).  The norm bound is the **Burnikel–Fleischer–
Mehlhorn–Schirra (2000)** separation bound.

---

## 23.2 The classical bound = our norm bound

Burnikel et al. (Algorithmica 2000), as stated in EHS §1:

```
|E| ≥ ( n · max_i (|x_i|·√a_i) )^{−(2^n − 1)} .
```

This is **identically** `Docs/02`: the `2^r − 1` exponent is the norm of `S`
over the degree-`2^r` field.  So `Docs/02` is the known state of the art
among *explicit* bounds before 2023, and `Docs/22` Theorem C explains exactly
why: this bound is worst-case optimal for isolated forms.

---

## 23.3 The state of the art: the polynomial exponent is SOLVED

**Eisenbrand–Haeberle–Singer**, arXiv:2312.02057 / SoCG 2024 (LIPIcs
SoCG 2024:54):

> **Theorem 8.**  Let `a_1,…,a_n` be distinct squarefree and `x ≠ 0 ∈ Z^n`.
> Then `E ≠ 0` implies
> ```
> |E| ≥ γ · ( n·‖x‖∞ )^{−2n} ,
> ```
> where `γ > 0` depends only on `√a_1,…,√a_n`.
>
> **"The constant `γ` is not explicit and stems from the subspace theorem."**

Translating to our variables (`n = m ≤ N`, `max|C_i| ≤ M`, `log₂M ≤ N`):

```
E  =  log₂(1/|S|)  ≤  2n·log₂(n·max|C_i|)  +  log₂(1/γ)
                  ≤  2N·(log₂N + N)  +  log₂(1/γ)
                  =  O(N²)  +  log₂(1/γ) .                            (23.1)
```

**Remark 9 of EHS** sharpens this further: the exponent `2n` can be replaced
by `n+ε` for any `ε > 0` (at the price of making `γ` depend on `ε`), giving

```
E  ≤  (n+ε)·(log₂n + log₂M)  +  log₂(1/γ_ε)  =  O(N²)  +  log₂(1/γ_ε) .
```

> **Consequence (major, and it validates `Docs/22`).**  The *exponent* half of
> the SRS conjecture is a **solved, published theorem**.  What remains is
> **exactly** `log₂(1/γ) ≤ poly(N)` — the single ineffective constant.  Our
> `Docs/03` reduction and `Docs/22` Theorem D are therefore not merely
> plausible: they match the published frontier precisely.

Their bound is also asymptotically almost tight: pigeonhole (their Thm 10)
gives `x` with `‖x‖∞ ≤ L` and `0 < |E| ≤ n·max√a_i / L^{n−1}`, matching
`L^{−(n−1)}` against the lower bound `L^{−2n}` up to the constant.

---

## 23.4 The best *explicit* constant, and where our M7.3 beats the literature

**Aymone–Figueredo–Iyer–Táfula**, arXiv:2609.14161 (Sept 2026) — the most
recent result, and the closest to our own work:

> For `K ≥ 8` and nonzero integers `m_1,…,m_K`,
> ```
> |Σ_{n≤K} m_n√a_n|  >  e^{ (2^{K−1}−1)/K − 1/2 } · ( max_n|m_n|√a_n · √K )^{−(2^{K−1}−1)} .
> ```
> "This inequality improves the dependence on `K` in the classical product
> bound."

Translating (with `M = max(|m_i|√a_i)`):

```
E  ≤  (2^{K−1}−1)·( log₂M + ½log₂K )  −  (2^{K−1}−1)·log₂(e)/K  +  ½log₂(e) .   (23.2)
```

This is **fully explicit** and beats the norm exponent `2^K − 1`, with a
large explicit prefactor `e^{(2^{K−1}−1)/K}`.  But the exponent `2^{K−1} − 1`
is still **exponential in the rank**, so (23.2) does *not* prove SRS.

**The frontier is therefore a clean trade-off:**

| | exponent in rank | constant |
|---|---|---|
| norm bound (`Docs/02`) = Burnikel et al. | `2^r − 1` | explicit |
| **Aymone et al. 2026** | `2^{r−1} − 1` | **explicit**, prefactor `e^{(2^{r−1}−1)/r}` |
| **Eisenbrand–Haeberle–Singer 2023/24** | `2r` (or `r+ε`) | **INEFFECTIVE** `γ` |
| `Docs/21` (**M7.3, this project**) | `2^{min(k_A,k_J)} − 1 ≤ 2^{⌊m/2⌋} − 1` | **explicit** |

> **Our M7.3 strictly dominates Aymone et al. in its regime.**  M7.3 gives
> exponent `2^{⌊m/2⌋} − 1` where Aymone et al. give `2^{m−1} − 1` — and M7.3
> does so with *no* prefactor loss, by an elementary norm argument on the two
> agreement/difference subfields.  M7.3 only applies to two-small forms
> (`Docs/21` (H)), so this is a genuine but partial improvement over the
> literature, obtained by an argument no one in that literature uses.

**Cross-check on our certified data** (40-digit arithmetic):

```
r=3 witness d=(2,3,13), C=(-10,-98,51),  E=21.7517644:
   norm  bound  (exponent 7)      ->  6.43e-20   OK
   AFIT bound    (exponent 3)      ->  5.10e-08   OK   (actual |E|=2.83e-07)
   EHS forces    gamma <= |E|(n*max|x|)^(2n)  =  1.83e+08
r=4:  norm 1.19e-27 | AFIT 1.22e-10 | EHS gamma <= 5.83e+12
r=5:  norm 3.67e-57 | AFIT 1.05e-21 | EHS gamma <= 1.85e+14
```

All three published bounds hold on our data, and the `r=3` witness caps
`γ ≤ 1.83e8` for `d=(2,3,13)`.

---

## 23.5 Why the constant is hard (the effectivity literature)

The obstruction is the classical non-effectivity of the Thue–Siegel–Roth–
Schmidt method.  Precise state of play:

* **Schmidt (1972)** (`Norm form equations`, Ann. Math. 96:526–551), ST
  itself: the union of subspaces containing the solutions is **not**
  effectively determined.
* **Schmidt (1989)** and **Evertse–Ferretti (2013)** give **explicit bounds
  on the NUMBER `t` of exceptional subspaces**, e.g.
  `t ≤ (2d)^{26n}·ε^{−2}` for forms over a degree-`d` field.  This does *not*
  bound their **heights**, nor the solutions inside them.
* **Schmidt (1989) / Vojta (1989) / Faltings–Wüstholz (1994)**: there is a
  finite, **effectively determinable** collection `U_1,…,U_r` of proper
  subspaces, independent of `δ`, such that all solutions outside their union
  form a **finite but non-effectively bounded** set.  This is precisely the
  state of play: the descent target is explicit, the terminal set is not.
* **Koleda**, *Probabilistic effectivity in the Subspace Theorem*,
  arXiv:2411.01247 / Res. Number Theory 2026: "One of the most fundamental
  open problems in Diophantine Analysis is to establish an effective version
  of this statement; namely, to determine the maximal height attained by the
  subspaces in the finite union."
* **EHS §4**, our exact question, verbatim: *"The subspace theorem does not
  provide explicit bounds on the number of solutions `q ∈ N⁺`.  The existing
  quantitative versions of the subspace theorem … do not provide such bounds
  either. **This is still the case when all algebraic numbers are square roots
  of integers.**  An explicit bound on the number of solutions would
  immediately apply to a separation bound for the sum of square roots."*

So: **our open target is a long-standing, explicitly stated open problem in
the literature, and it is open for exactly our special case.**

---

## 23.6 The mechanism EHS use (and why it is not a norm argument)

EHS do **not** improve the norm.  They route through simultaneous
approximation:

1. Build the lattice `Λ(B)` with `B` dual to the "fractional-part" lattice,
   dimension `n+1`, `det = Q^{−(n+1)}`.
2. Minkowski's second theorem gives a lattice vector `v = B(q, −p)ᵀ` with
   `‖v‖∞ ≤ 1/Q`, i.e. `dist_ℤ(qβ_i) ≤ 1/Q`.
3. If additionally `λ₁ ≥ Q^{−1−1/(3n)}`, then *all* successive minima are
   good, and one of the `n+1` lattice vectors is of the form `B(q,−p)ᵀ`, so
   `|βᵀx| ≥ 1/(2q) ≥ Q^{−(n+1)}` — **the separation bound, with exponent
   `n+1`, not `2^n`**.
4. The condition `λ₁ ≥ Q^{−1−δ}` with `δ = 1/(3n)` is exactly the
   Subspace-Theorem inequality `q^{1+δ}·Π_i dist_ℤ(qβ_i) < 1`.  Schmidt's
   Theorem 1B says only finitely many `q` satisfy it.  **`γ` is determined by
   the largest such `q`.**

So the entire constant is *"how far can `q` go before the Subspace Theorem
finiteness kicks in"*, and that is the non-effective step.  Note this is
**consistent with `Docs/22` Theorem C**: EHS obtain a polynomial exponent not
from the norm but from Minkowski + the Subspace Theorem, i.e. from a genuinely
different mechanism — and they are explicit that the residual constant cannot
be extracted.

---

## 23.7 A sharper target than "polynomial in the discriminant"

The natural guess is `γ ≥ Disc(K)^{−O(1)}`.  **That guess is wrong**, and it is
worth recording why, because it is a trap.

For a multiquadratic `K = Q(√d_1,…,d_r)`,

```
Disc(K) = ∏_{v ∈ Λ \ {0}} D(v),        D(v) = fundamental disc. of Q(√v),
```
and each `d_i` occurs in `2^{r−1}` of the `2^r − 1` nontrivial products, so

```
log₂ Disc(K)  ≍  2^{r−1} · Σ_i log₂ d_i  +  2·2^r ,                  (23.3)
```

which is **exponential in the rank**.  Measured:

```
r=3:  log2 Disc = 37.14   vs   log2(prod d_i) =  6.29   ratio  5.9
r=4:  log2 Disc = 85.71   vs   log2(prod d_i) =  7.71   ratio 11.1
r=5:  log2 Disc = 226.78  vs   log2(prod d_i) = 11.17   ratio 20.3
```

So `γ ≥ Disc(K)^{−O(1)}` would only give `log₂(1/γ) = O(2^{r}·N)` — still
exponential in `r`, hence useless.

> **Sharp target (M7.4, restated).**  ~~Prove~~
> ```
> gamma  >=  ( d_1*...*d_m )^(-C)      for some absolute C,
> ```
> equivalently `log2(1/gamma) <= C * sum_i log2 d_i <= C * N`.

> ## RETRACTED BY `Docs/24` (Theorem 24.4)
>
> **This target is FALSE.**  No absolute `C` exists.  The exceptional `q`-set
> provably contains convergent denominators up to `2^{3n(n-1)}`
> (`Docs/24` Thm 24.3), forcing `log2(1/gamma) = Omega(n^3)`, while distinct
> squarefree radicands give `sum log2 d_i >= log2(n!) = Theta(n log n)`.  The
> ratio diverges.  Refuted by an elementary Pell argument.
>
> **Corrected target (`Docs/24` Thm 24.5):**  it suffices to prove
> ```
> log2 q*  =  O( n^2 + n*log2(max_i d_i) + n*log n ) ,
> ```
> equivalently `q* <= exp(O(N^2))`, which yields **SRS with exponent 3**.
> The `n^2` scale is forced, so the target is tight up to constants.

Equivalently, in EHS's own terms (their closing sentence): produce an
effective upper bound on the `q ∈ N⁺` satisfying
`q^{1+δ}·Π_i dist_ℤ(q√a_i) < 1` with `δ = 1/poly(n)`.

*(That last formulation is the one `Docs/24` attacks directly; it obtains a
lower bound and a refutation of the scale above, but no upper bound.)*

---

## 23.8 Assessment

**Solved / settled by the literature:**
* the exponent `2^r − 1` is replaceable by `2r` or `r+ε` — EHS;
* the classical norm bound is optimal in the isolated regime — `Docs/22` Thm C;
* the `2^{r−1}−1` explicit improvement — Aymone et al.;
* the zero test and linear independence — Besicovitch, Blömer.

**Ours, and apparently not in the literature:**
* **M7.3** (`Docs/21`): elementary exponent `2^{min(k_A,k_J)} − 1 ≤ 2^{⌊m/2⌋}−1`
  for two-small forms, dominating Aymone et al.'s universal `2^{r−1}−1` in
  that regime, with an elementary proof and no lost prefactor.
* **M7.4** (`Docs/22`): the no-go theorem that isolation is inert for
  norm/Fourier arguments, and hence that the isolated case *must* be attacked
  through effectivity rather than algebraic refinement.

**Still open, and open in the literature:**
```
γ ≥ (d_1···d_m)^{−C} ,   equivalently   log₂(1/c(d,δ)) ≤ poly(N).
```

**What would count as progress (in decreasing strength):**
1. Any explicit bound on the EHS `q`-set for `β_i = √a_i` — immediately yields
   an explicit polynomial-exponent separation bound, settling SRS.
2. A bound `log₂(1/γ) = O(r·N)` instead of `O(N)` — still polynomial in `N`
   for fixed `r`, and would extend the two-small/isolated dichotomy to a full
   polynomial-in-`N` statement at bounded rank.
3. Making the *descent* explicit: since `U_0` is effectively determinable and
   the rank strictly drops at each descent (≤ `r` steps), the only
   non-effective residue is the terminal finite set.  A height bound on that
   terminal set is exactly item 1.

**Do not re-propose** (all closed, with citations): pure-norm refinements
(`Docs/22` Thm C), norm towers (`Docs/10`), sign-pattern rank descent
(`Docs/13`), generic Baker/Matveev at polynomial scale (`Docs/16`), and
any attempt to use the isolation hypothesis or the Fourier/Hadamard structure
of the conjugates (`Docs/22` Prop. B2, Cor. C2).