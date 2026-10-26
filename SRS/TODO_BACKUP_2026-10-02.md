# SRS — TODO

Staged work with stopping rules. This records **what is decided and verified**, not what is
merely intended. Every claim below was executed; anything not run is under *Open*.

A checked box means the item is done **and** a test or a printed run backs it. An earlier
revision of this file marked a "doc pass" complete when no doc had been touched, and credited
a Rust `--bits 1000000` flag that did not exist; both errors are corrected here.

**Third correction, and the largest: a `[X]` on a *literature* claim means the
claim was read at the source, not that it was derived here.** The Eisenbrand–
Haeberle–Singer entry below is marked done because the paper was read in full
(arXiv:2312.02057), and it is labelled **QUOTED** everywhere it appears — the
repo does not prove it. An earlier draft marked that entry done on the basis
of a *negative* claim about the literature ("this route is closed") that was
never checked and turned out to be false. Superseded entries are retained and
marked `[SUPERSEDED]` rather than deleted, so the retraction record stays
visible; there are six substantive ones, listed in
`THE_MISSING_THEOREM.md` §6.1.

**Fourth correction — a process failure, not a maths one.** The `r ≥ 4` bound
claimed in the §2.11.9 pass was false, and the step that made it false
(`dist(aW, Z) = a·dist(W, Z)`, which needs `a Z = Z`) was asserted twice and
"corrected" twice without ever being tested. It was caught only when the bound
was finally run against 910 000 certified instances and produced 4 503
violations. **The rule this repo now follows: a separation bound is not
claimed until a test has tried to refute it and failed to.** Three claims in
one session were wrong, all in the same step; the fix is the test, not more
re-reading. The retraction is recorded in `PROOFS.md` §2.11.9 and pinned by
`tests/test_proofs.py::TestEffectiveConstant::test_RETRACTED_general_bound_is_false`.

**Navigation.** `STATE_OF_THE_ART.md` (the correction) →
`THE_MISSING_THEOREM.md` (the roadmap and the open problem) →
`PROOFS.md` (the proofs) → this file (status).

---

## Theory pass — proofs written out (PROOFS.md)

- [X] **`PROOFS.md` Thm 2.9 — the exponent is settled (new, this pass).**
  Schmidt's Subspace Theorem applied to the **single** linear form SRS is,
  using the `r` forms `L_0 = L`, `L_1..L_{r-1} = C_1..C_{r-1}` (exactly `r`
  forms in `r` variables; `C_r` is not a form and is not needed):
      `|L(C)| > c(d,delta) * ||C||_inf^{-(r-1+delta)}`
      => `E/N <= (r-1+delta)/r`.
  → **This is SHARPER than EHS Thm 8** (exponent `r-1` vs their `2n`),
    because SRS is one form and a single form in `r` variables has exponent
    `r-1`. At `r=3`: `0.667+delta/3` vs EHS `1.0` vs Thm 2.6 `1.5`.
  → **It MATCHES the achieved rate** `(r-1)/r` of Thm 2.5. So the guaranteed
    and achieved rates **coincide at every rank**, and the sharp form of the
    conjecture is the *asymptotic* `limsup_{N→∞} E/N = (r-1)/r`. There is no
    growth-rate gap left.
  → **[RETRACTED]** The `sup` version of that claim is **false**: a certified
    `r=3` instance `d=(2,3,13)`, `C=(-10,-98,51)` has `E/N = 1.2795 > 2/3`
    (`PROOFS.md` §2.10.2). The excess is additive slack `log2(1/c) ≈ 21.75`
    over `N=17`; measured over 1140 triples, mean `E/N` falls from `0.76`
    (`N<20`) to `0.66` (`N>=20`), converging to `2/3` from above.
  → **The constant, located** (`PROOFS.md` §2.10). It is a Markov constant:
    `c = 1/limsup_q q·||qα||` with `α=√(d₂/d₁)` at `r=2`, computable from a
    periodic continued fraction. Thm 2.11: **any** `c ≥ 2^{-poly(N)}` suffices
    — `c = Ω(1)` is NOT required. The remaining gap is exactly
    `poly(r)`-exponent **and** effective, simultaneously.
  → **Consistency check the repo owns**: at `r=2` the two forms are
    `L_0 = C_2 + C_1 alpha_1`, `L_1 = C_1` — that IS Roth, exponent `1 = r-1`.
    And Thm 2.2/2.6 prove the `r=2` case sharp at `1/2` with an EXPLICIT
    constant. So the application is known-correct at the one owned rank.
  → Measured on the achieved side: `E/log2||C||_inf` tends to `r-1` at both
    ranks (`1.434→1.090` at `r=2`; `2.312→1.999` at `r=3`).
  → **What it does NOT do**: bound the constant `c(d,delta)`, which depends on
    the radicands and is therefore input. Same gap as EHS `gamma`. Thm 2.9
    improves the exponent only.
  → **The exceptional subspaces are DETECTABLE**: `C` lies in a proper
    rational subspace iff some nonzero integer `a` has `a.C = 0`, decidable
    by exact linear algebra over Q. So the generic case already has a
    **constant-free** bound at the exact exponent `r-1`.

- [X] **`SRS/families.py` F6 + `SRS/refutation.py` (new, this pass).**
  → **F6**: the Thm 2.5 pigeonhole-gap construction at `r=3`. Every other
    family searches for instances *worse* than expected; F6 provably attains
    the ceiling, so it calibrates what "at the ceiling" looks like and is the
    only family that could refute the sharp form. Exactly certified via the
    integer kernel; respects Thm 2.5's bound; measured `E/N = 0.6153…0.6835`
    at the `2/3` ceiling.
  → **`refutation.py`**: hunts for `E/N > (r-1)/r` at `r=2..5` with random and
    greedy strategies plus F6, using **certified** signs and magnitudes.
    **Result: no refutation found.** A negative result, reported as such.
  → **A false positive worth recording**: the first version compared `E/N`
    against the bare limit `(r-1)/r` and reported **F6 — Thm 2.5's own
    construction — as a refutation**. Cause: ignoring the construction's
    `sum sqrt(d_i)` prefactor, an `O(1/N)` rate term. Fixed by requiring
    `rate > ceiling + 2*slack`; pinned by `test_f6_is_not_a_refutation`.

- [X] **`STATE_OF_THE_ART.md`** — the literature check, and a major
  correction. The best known separation bound is **polynomial** in the number
  of terms (EHS, arXiv:2312.02057; SoCG 2024, via Schmidt's Subspace Theorem
  1B), not `O(2^r)`. The `2^r` in Thm 2.6 is the optimum **of the norm**; the
  problem is not so bounded. The one remaining gap is the **effectivity of the
  constant `γ`**, which depends on the radicands — and those are input. This
  file supersedes `THE_MISSING_THEOREM.md` §3 and `PROOFS.md` §2.8 as
  originally written, and records the form-counting error that produced the
  wrong conclusion.
  → Also records two things the repo got right on its own: Thm 2.5 **is** EHS
    Thm 10, and the sharp `k=2` constant `1/2` is **stronger** than EHS at
    that rank.

- [X] **`THE_MISSING_THEOREM.md`** — the longer roadmap: the problem
  statement, both proved sides with the `E/N` bracket, seven routes tried with
  why each fails, the retraction log, and a "what a proof would have to
  contain" target list. `PROOFS.md` holds the proofs; this holds the map;
  **`STATE_OF_THE_ART.md` corrects it.**

- [X] **`PROOFS.md`** — results that were previously asserted or only
  measured, now with written proofs:
  → **Thm 1**: the F₂ basis theorem + corollary that distinct squarefree
  radicals are ℚ-independent even when dependent as generators (`{2,3,6}`),
  which is the correctness half of the zero-test-in-P row. Lemma 1.1
  (irrationality step) by induction; the `u = 0` branch needs the parity of
  `s·t_r`, not of `s` — that is where a naive writeup would silently assume
  squarefreeness it does not have.
  → **Thm 2.1**: the general separation guarantee
  `|S| ≥ (k·2^{3n/2})^{−(2^k−1)}`, proving the README's "best known
  `2^{−Θ(n·2^k)}`" row by an explicit Liouville argument (algebraic integer,
  degree ≤ 2^k, conjugate bound, norm ≥ 1).
  → **Thm 2.2 (RESEARCH TARGET 1 — CLOSED)**: the two-term case is now
  complete. Guaranteed `log₂(1/|S|) ≤ N/2 + 3/2`; achieved by the
  coefficient-Pell family `2C₁² − 3C₂² = −1` with radicands fixed at 2 and 3,
  where the coefficients carry all the bits and `log₂(1/|S|)/N → 1/2`.
  **So the worst two-term constant is exactly 1/2.** Two corrections fell out
  of this: (a) Thm 2.1's inequality direction was initially written backwards
  (`≥` for `log₂(1/|S|)`; the guarantee is an *upper* bound on it) and is
  fixed; (b) the earlier claim that "F1 achieves the worst two-term constant,
  which lies in [1/4, 9/4]" was **wrong twice over** — the constant is 1/2,
  and F1 reaches only 1/4. F1 hides its Pell arithmetic in ~10⁶-bit
  radicands; the coefficient family with tiny fixed radicands is strictly
  easier to find and strictly harder to decide. This is a new instance family
  the probe does not have — see Open.
  → **Thm 3**: the F1 law exact — the residual is
  `ρ_n = log₂((1+s_n)/(2√s_n)) > 0` (AM–GM, one-sided), bounded by
  `0.205·a²/p_n⁴` and asymptotic to `D_n²/(32 ln2 p_n⁴)`. The first-order
  corrections cancel on the two sides of the law, which is why the measured
  residual (`7.5e-5`) is ~400× smaller than either correction (`~0.03`);
  the asymptotic predicts `7.18e-5` at the worst state, within 5%.
  `_SEPARATION_TOLERANCE = 1e-3` is now a proved margin, not an empirical one.
  → **Cor 3.3**: the Θ(N) law with constants:
  `|log₂(1/|S_n|)/N_n − 1/4| ≤ (log₂a + 1.002)/N_n`.
- [X] **`tests/test_proofs.py`** — corroboration (not evidence): Liouville
  bound on random and deep-cancellation instances via exact rational
  brackets; the F1 residual identity, one-sidedness, and bound at 100-digit
  precision (float64 reads ρ as exactly 0 past n ≈ 10 — the first
  verification attempt failed for exactly that reason, which is now recorded
  in the test docstring); the two-term bound on eight instance shapes
  including both degenerate branches; and that the coefficient-Pell family's
  ratio falls monotonically onto 1/2 while F1's falls onto 1/4.
- [ ] **RESEARCH TARGET 2 — the conjecture itself**
  (`sign(S)` in polynomial time, i.e. a `2^(−poly(N))` separation guarantee).
  Still open, and this pass sharpened both edges without closing it:
  → **Thm 2.3 (new)**: a second-moment separation bound
    `log₂(1/|S|) ≤ ((2^r − 1)/2)·(log₂ M + 1)` with `M = Σ C_i²d_i` and `r` the
    F₂-rank of the squarefree parts. Uses the *exact* orbit second moment
    `Σ_σ σ(S)² = 2^r M` instead of the per-conjugate maximum, halving
    Thm 2.1's exponent; measured 3.4–5.0× sharper at `r = 2, 3`. Still
    carries `2^{r-1}` in the exponent, so it is **not** a step toward the
    conjecture — no norm-based bound can be, since the norm has `2^r`
    conjugates to absorb.
  → `r`, not `k`, is the parameter that matters (both F1 and the new
    coefficient-Pell family sit at `r = 2`; a large `k` over a small `r` buys
    the adversary nothing).
  → **Honest retraction**: an earlier draft of Thm 2.3's remarks claimed the
    bound was "attained to three decimals" on the coefficient-Pell family.
    False — the slack *grows* (3.3 → 23.2 bits over four family members),
    because the norm/AM–GM step spends the `±S` pairing of the orbit. Thm
    2.2 is 3× sharper at `r = 2`. Both the remark and the test that pins it
    (`test_bound_attained_on_coefficient_pell_is_false`) were corrected.
    **Superseded**: Thm 2.6 now spends that pairing *properly* and is the
    sharp invariant; the retraction above still stands for Thm 2.3 itself.
  → **Thm 2.6 (new, and the strongest norm bound in the repo)**:
    `log₂(1/|S|) ≤ ((2^{r−1} − 1)/2)·(log₂ M + 1)`. The global sign flip
    `τ: sqrt(d_i) ↦ −sqrt(d_i)` is an automorphism with `τ(σ(S)) = −σ(S)`,
    so the `2^r` conjugates pair up, the norm is a square up to sign, and
    `Π_{half}|σ(S)| ≥ 1` is free. AM–GM over the `2^{r−1}−1` remaining
    half-orbit elements then halves Thm 2.3's exponent again, by a
    **uniform** factor `(2^r−1)/(2^{r−1}−1) → 2` (`3.000, 2.333, 2.143,
    2.067, 2.032, 2.016, 2.008` at `r = 2..8`).
    → **The `k = 2` constant is now tight from BOTH sides.** Thm 2.6 at
      `r = 2` is `½(log₂M+1)`, which on the coefficient-Pell family
      `2C₁²−3C₂²=−1` is *exactly* Cauchy–Schwarz `X+Y ≤ sqrt(2(X²+Y²))` in
      the exponent, and the slack decays geometrically (`7.3e−3, 7.7e−7,
      8.0e−11, 8.3e−15, 8.7e−19, 9.0e−23, 9.4e−27` at `t = 0..6`). Its own
      sup-ratio over the imbalance ratio converges to `1/2` from above
      (`0.8737, 0.6444, 0.5648, 0.5374, 0.5182, 0.5092, 0.5049` at size
      `10^{3..150}`), matching Thm 2.2's achieved `1/2`.
    → **Correction to this file's own earlier claim**: §2.5 previously said
      the `±S` pairing "stops at `r = 2`". False. That was true only of the
      *single-conjugate* route (`M_j ∈ K_{r−1}`), which is a different and
      worse pairing. The global flip works at every `r`. PROOFS.md §2.5 and
      the retraction note below were both rewritten.
    → **Second correction**: Thm 2.6 does *not* dominate Thm 2.2 pointwise
      at `r = 2`. It wins iff `2C₁² + 3C₂² < 24C₁C₂` (comparable terms) and
      loses on lopsided ones (`C₁=10⁶, C₂=1`: `20.93` vs `12.76`). The two
      are complementary; the usable guarantee is the min, whose sup-ratio
      also tends to `1/2`. A draft claimed outright domination — checked only
      on the balanced family — and it is withdrawn.
    → **But the growth rate of the norm class is untouched**: the exponent is
      `≈ 2^{r−2}`, still `Θ(2^r)`. Thm 2.6 halves a constant, not the growth
      rate, so it is still not a step toward the conjecture *by norm means*.
      No further pairing exists in this class: after the global flip the
      remaining `2^{r−1}−1` elements carry no forced symmetry.
      **[Superseded in scope by the entry below: a non-norm argument does
      reach a polynomial bound, so `Θ(2^r)` is a property of the method and
      not of the problem.]**
  → **PROOFS.md §2.8 / Thm 2.7–2.8 — MAJOR CORRECTION, superseding the
    earlier entry below.** The Subspace-Theorem route is **NOT closed**; it
    gives the best known bound in the literature. Eisenbrand–Haeberle–Singer
    (arXiv:2312.02057; SoCG 2024) apply Schmidt's Theorem 1B to this exact
    problem and obtain
      `|E| >= gamma * (n * ||x||_inf)^{-2n}`,  exponent sharpenable to `n+eps`
    ⇒ `E <= 2n(log2 n + m) + log2(1/gamma)` with `N = n*m`, i.e. **polynomial
    in the number of terms**, `E/N → 1` after the `n+eps` sharpening.
    → **The exponential-in-`r` barrier is NOT the problem.** Thm 2.6's `2^r`
      is the optimum **of the norm**, not of the problem. EHS's bound never
      takes the norm of the whole orbit. Measured dominance of Thm 2.6 over
      EHS: `0.98, 2.4, 12, 37, 788, 1.9e4` at `r = 3,5,8,10,15,20`.
    → **The `r=2` exception is real and is the repo's own result**: at `r=2,3`
      the two bounds are comparable and Thm 2.6 is *sharper*, because EHS
      does not know the radicands. The sharp two-term constant `1/2` is
      **ours**, not theirs.
    → **Thm 2.5 IS EHS Thm 10, rediscovered** (`Σ sqrt(d_i) <= n*max sqrt(d_i)`,
      same exponent `n-1`). Independent agreement on the achieved side.
    → **THE remaining gap is the effectivity of `gamma`.** It is a function
      of `a_1..a_n`, and in SRS the `a_i` are **input**. Concretely: bound
      `Q_0(a_1..a_n)` so `log2 Q_0 <= poly(N)`, equivalently bound the number
      of `q` with `q^{1+delta} * prod_i dist_Z(q*sqrt(a_i)) < 1`. EHS §4 state
      the Subspace Theorem is ineffective here and that quantitative versions
      (Evertse–Ferretti 2013) "do not provide such bounds either … even when
      all algebraic numbers are square roots of integers". That special case
      is the most concrete lead in the repo.
    → **Why this is not yet SRS ∈ P**: deciding sign needs a *known* lower
      bound to set the precision threshold; `gamma` is unknown. Polynomial
      precision with an unusable constant.

  → **[SUPERSEDED — retained only to show the error] PROOFS.md §2.8
    originally claimed the Subspace route was CLOSED**, on the grounds that
    Schmidt's theorem needs `n` separate forms and SRS has one. **That is a
    form-counting error.** In the `n+1`-form version, the `n-1` coordinate
    forms `q*b_i - p_i` are supplied by the reduction itself; only `L_0` must
    come from the problem. The product is `q * prod_i|q*b_i - p_i|`, which is
    exactly Schmidt Thm 1B. I had the right theorem in hand and miscounted
    the forms. Recorded because the same mistake could recur, and because a
    negative claim about the literature needs the same checking as a positive
    one — this one did not get it before. Full account in
    `STATE_OF_THE_ART.md` §0.
  → **[SUPERSEDED IN EXPLANATION by Thm 2.9 above] PROOFS.md §2.4 / Prop 2.4′ (new, corrects the file)**: the support of
    `N_j = S_{j−1}² − C_j²d_j` has F₂-rank `max(j−2, 0)`, **not** `j−1` as
    an earlier version of the file stated. Verified exactly over 300+
    random F₂-independent sets per level, `j = 2..7`. The `r = 3` instance is
    the identity `S·σ₃(S) = (C₁²d₁+C₂²d₂−C₃²d₃) + 2C₁C₂sqrt(d₁d₂)`, a
    two-term expression in a quadratic field.
    → **But the barrier conclusion is unchanged, for a different reason**: the
      factor-2 rank saving is cancelled *exactly* by the squaring that creates
      the shallower subfield (`log₂ M' ≈ 2 log₂ M` instead of `≈ log₂ M`).
      Measured `tower(rank j−2) / Thm 2.6` over `r = 3..8`:
      `1.0049, 0.9982, 1.0007, 0.9975, 1.0019, 1.0003` — `1 ± 0.005`, no gain
      and no loss. The file's old *reason* was wrong; its *conclusion* holds,
      and the correct reason is now written down.
  → **[SUPERSEDED — kept for the retraction record] The gap was once stated
    as: guaranteed `≈ 2^{r−3}` per unit `N` (Thm 2.6), achieved `= 1/2` at
    `k = 2` (Thm 2.2) and `= 1/4` for F1 (Cor 3.3), with the gap exponential
    in `r`.** That framing is wrong: the guarantee is *not* exponential in
    `r` in the literature (see the §2.8 entry above). The **achieved** half
    of that statement stands and is unaffected: `(r−1)/r` is the algebraic
    ceiling, `1/2` at `k = 2`, `1/4` for F1. Only the guaranteed half is
    superseded.]
  → **FOUNDATION (Thm 2.4, new): the tower is quadratic, and the barrier is
    therefore *not* combinatorial.** The relative norms
    `N_j = S_{j−1}² − C_j²d_j` carry exactly `1 + C(j−1, 2)` radical
    monomials — the 2-element subsets of `{1,…,j−1}` plus one rational — not
    the `2^j` a naive resultant tower suggests. Measured over 200 towers,
    exact. This kills the hope that monomial-count economy buys a
    polynomial bound: the monomials live in `K_{j−1}` of rank `j−1`, so
    separating `N_j` still costs `Θ(2^{j−1})` and applying Thm 2.3 to it
    reproduces the same exponential. **Stated consequence for future work:
    every norm-based route (Liouville, Thm 2.3, Thm 2.6, resultants, the
    quadratic tower) pays `2^{rank}` and none of them can reach a polynomial
    bound. A bound that does must use something that is not a norm over the
    whole orbit.**
    → **[This prediction was CORRECT and is now realised.]** §2.8 below
      records that Eisenbrand–Haeberle–Singer reached a polynomial-in-`n`
      bound by exactly such a non-norm route (Schmidt's Subspace Theorem).
      So the negative result above stands as a statement about the norm class
      — which is a real and useful thing to have established — and the
      prediction it made is the one that pointed at the way through.
  → **PROOFS.md §2.5 (new): what the `r = 2` argument actually is.** The
    `N/2` rests on exactly two ingredients: (i) pairing `S` with a *single*
    conjugate gives `S·σ(S) = C₂²d₂ − C₁²d₁ ∈ ℤ`, hence `|S| ≥ 1/(X+Y)` —
    one algebraic integer of rank 0, not a `2^r`-orbit product; (ii) the
    elementary `max(X,Y) ≤ 2^{N/2}`.
    → **The old "it stops at `r = 2`" conclusion was WRONG and is retracted.**
      The single-conjugate pairing does stop (`M_j ∈ K_{r−1}`, costing
      `2^{r−1}`), but that is the wrong pairing. The global flip is not a
      single conjugate and works at every `r` — that is Thm 2.6 above.
    → **Measured, not asserted**: the *single-conjugate* route gives 3.54× at
      `r = 3` on random instances but decays to **1.05× by `r = 8`** on
      equal-size terms (ratio `1 + O(2^{-r})`). It converges to Thm 2.3, does
      not beat it, so it is **not** recorded as a theorem. It is kept in the
      file to document why that pairing is the wrong one.
    → **Two false starts recorded so they are not repeated**: (a)
      *comparability* — the `x₁ ~ … ~ x_r` generalisation is false
      (`max|C_i|sqrt(d_i)/|S|` reaches `7e13`); only the weak "tail cancels
      head" form survives, and Thm 2.2 never needed individual
      comparability; (b) *cheap monomials ⇒ cheap separation* — Thm 2.4's
      quadratic count does not help, because the monomials live in rank
      `j−1`, so the recursion still pays `Θ(2^{j−1})` per level.
  → **PROOFS.md §2.6 / Thm 2.5 (new): the achieved rate is bounded below by
    `(r−1)/r`, by counting alone.** The `(B+1)^r` sums
    `Σ C_i sqrt(d_i)` over `C ∈ {0..B}^r` are **pairwise distinct** (Cor 1.3),
    lie in an interval of length `B·Σsqrt(d_i)`, so a consecutive gap gives
    `|S| ≤ (Σsqrt d)·B^(−(r−1))`, i.e. `E/N ≥ (r−1)/r − o(1)`.
    → **Tight at `r = 2`**: Thm 2.5 gives `(2−1)/2 = 1/2`, which is exactly
      Thm 2.2's proved sharp constant — an algebraic argument and a
      combinatorial one landing on the same number, which is the strongest
      evidence yet that `(r−1)/r` is the true achieved constant.
    → The achieved rate **grows** with `r` (`2/3, 3/4, …`) while the
      guaranteed rate grows like `2^r`, so the gap *widens*: ratio
      `3.0×, 3.5×, 5.0×, 7.8×` at `r = 2,3,4,5`.
    → **Buildable family**: unlike F1 (a Pell solve per instance), Thm 2.5's
      instances come from sorting `(B+1)^r` sums and taking a consecutive
      gap. Polynomial time in `B^r`, and the resulting `C` is checkable by
      the existing exact decider. This is the strongest adversarial family
      the repo now has on the achieved side.
    → **Honest limit**: Thm 2.5 uses the `x² − d_i` only through Cor 1.3's
      independence. It is pure counting, so it cannot be improved without
      importing arithmetic structure — and none appears.
- [ ] **New instance family for the probe** — the coefficient-Pell family
  (Thm 2.2) is strictly harder per input bit than F1 and is not generated by
  any of F1–F5. Adding it is mechanical (solve `2x²−3y²=−1` by the unit
  `5+2√6`, emit `[(C1, 2), (−C2, 3)]`); the value is that it probes the
  decider where the guaranteed cost is highest.

---

## §2.11 program — can `2^{r−1}−1` become `poly(r)`? (newest pass)

The question this pass was opened to answer: replace the norm bound's exponent
`2^{r−1}−1` by `poly(r)` while keeping the constant effective. Answered, and the
answer is *partly yes, partly a proof of impossibility*. Full write-up in
`PROOFS.md` §2.11; tests in `tests/test_proofs.py::TestPolyRankExponent`.

- [X] **Thm 2.12 — every fixed rank is already settled, effectively.**
  `E ≤ (2^{r−1}−1)·(N + log₂(3r))`. Since `2^{r−1}−1` is a *constant* at fixed
  rank, the conjectured polynomial separation is a **theorem rank-by-rank**
  (rank cap 2 → `E ≤ N`, 3 → `3N`, 4 → `7N`, 5 → `15N`, 8 → `127N`).
  → **This is the reframe that matters**: the open problem is *not* the
  exponent as a function of `r`. It is **uniformity in `r`**, since
  `r = Θ(N)` is achievable (take `d_i` = i-th prime, `C_i = 1`).
- [X] **Thm 2.13 — the exponent `2^{r−1}−1` is IRREDUCIBLE.**
  It is the number of complement pairs `I/Iᶜ` among the `2^r` conjugates, and
  AM-GM forces it: indexing the `M = 2^{r−1}−1` pairs by `J ⊆ {2..r}`,
  each index lies in exactly `2^{r−2}` of them, so the mean factor is
  `< |S| + T` and `E ≤ M·log₂(2T)`.
  → **No choice of representatives, and no Jensen/AM-GM reweighting, can
  improve it.** The field has `2^{r−1}` conjugate pairs; you must bound every
  one. Verified at `r = 2,3,4,5`.
- [X] **Thm 2.14 — squaring is a no-go.**
  `S²` does live in a field of degree `2^{r−1}` (the products `d_i d_j` span the
  **even** subspace, dimension `r−1`; for `d=(2,3,13)`, `6·26·39 = 78²`). But
  `R = ΣᵢCᵢ²dᵢ` is a sum of **non-negative** integers with `d_i > 0`, so
  `R ≠ 0` for *every* instance, and `σ_{−ε}(S²) = R − irr ≠ −S²`. The `|σ|²`
  squaring is lost, the factor count is unchanged, and `E(S²) = 2E(S)` against a
  doubled `log₂T` — **the two cancel exactly.** Iterating fails too: `S⁴`
  re-introduces a rational part and the even subspace does not nest.
  → closes the whole field-doubling family of tricks.
- [X] **Thm 2.15 — at `r = 2` no additive constant is needed at all.**
  For *any* irrational `α`, convergents give `qₙ‖qₙα‖ ≤ 1`, so
  `M(α) = limsup q‖qα‖ ≤ 1` and the asymptotic constant `1/M(α) ≥ 1`, i.e.
  `log₂(1/c) ≤ 0`. The finite-`C` part (`‖C‖_∞ ≤ 8`, brute force) costs
  `≤ 3.67` bits over `d ≤ 20000`, worst case `d=(13,3)`.
  → **Stronger than §2.10.1's "computable from the periodic CF": the budget is
  zero.** `r = 2` is fully effective with `O(1)` bits of slack, independent of `N`.
  → Caveat pinned by a test: this rests on Thm 2.2's **direct** two-term
  argument, not on the Markov constant. The earlier justification "`M(α) ≤ 1`
  is the limsup, so `c ≥ 1`" was **false** (§2.11.4: `q = 3161` gives
  `q‖qα‖ = 2.4996` for `α = √(5/3)`), and it is retracted. Thm 2.9's `c` is an
  inf over **all** generic `C`, so in any case the small-`C` part is a separate
  finite search.
- [X] **§2.11.5 — at `r = 3` the target is far stronger than `2^{−poly(N)}`.**
  Measured `c(d) = min ‖C₂α+C₃β‖·max(|C₂|,|C₃|)²` over exhaustive boxes:

  | radicand range | triples | `min c` | `max log₂(1/c)` | `median` |
  |---|---|---|---|---|
  | `d ≤ 50` | 455 | `2.9×10⁻⁵` | 15.1 | — |
  | `d ≤ 100` | 2300 | `3.3×10⁻⁶` | 18.2 | — |
  | `d ≤ 200` | 15180 | `1.7×10⁻⁶` | 19.2 | — |
  | `d ≤ 400` | 76076 | `1.7×10⁻⁶` | 19.2 | — |
  | `d ∈ [2,500)` | 4000 | `2.3×10⁻⁷` | 22.0 | 5.6 |
  | `d ∈ [500,2000)` | 4000 | `4.1×10⁻⁶` | 17.9 | 5.8 |
  | `d ∈ [2000,6000)` | 4000 | `8.7×10⁻⁶` | 16.8 | 5.9 |
  | `d ∈ [6000,20000)` | 4000 | `7.2×10⁻⁶` | 17.1 | 5.8 |

  `log₂(1/c)` does **not** grow with `d` across four orders of magnitude, and the
  worst case is a genuine infimum, not a search artefact: for `d=(263,311,313)`,
  `C=(−11,−26,36)` (certified sign `+1`, `N=27`, `c=1.7258×10⁻⁶`,
  `log₂(1/c)=19.14`) the value is **identical** for `B = 60,120,240,480,960,2000`.
  Since `c(B)` is non-increasing in `B`, the flatness *is* the infimum.
  → **Evidence, not proof**, labelled as such: a uniform `c = Ω(1)` at `r=3`
  would be a theorem about Markov constants of vectors of quadratic
  irrationals, and no such theorem is known. It does establish that the target
  is not razor-thin: `log₂(1/c) ≤ 22` everywhere tested, so Thm 2.9's constant
  is plausibly bounded by an absolute constant rather than degrading with input.
- [X] **Thm 2.17 — the exponent is the ALGEBRAIC DEGREE (the real obstruction).**
  `[Q(S):Q] = 2^r` **exactly** for `S = ΣCᵢ√dᵢ` with all `Cᵢ ≠ 0`. The stabilizer of `S`
  in `Gal(K/Q) ≈ (ℤ/2ℤ)^r` is trivial by Thm 1, so the orbit has `2^r` points.
  → **Stronger than Thm 2.13**: the `2^{r−1}−1` is not a loose estimate that
    better bounding could tighten, it is the degree of the number. **No norm
    argument of any kind reaches `poly(r)`** — this rules out the whole family,
    not just reweightings. (Thm 2.9 escapes it precisely by *not* being a norm.)
  → **Sharp, and the sharp case is the input normalization**: zero one `Cᵢ` and
    the degree is exactly `2^{r−1}`. That is why zero coefficients are
    aggregated away — a surviving zero would falsify Thm 2.17 *and* overcount
    the norm bound. Pinned by `test_thm217_degree_drops_exactly_when_a_coefficient_vanishes`.
- [X] **Thm 2.18 / Cor 2.19 — the conjecture is PROVED for `r = O(log N)`.**
  `N ≥ Σlog₂dᵢ ≥ log₂((r+1)!) = r·log₂r − 1.4427r + O(log r)`, and `r ≤ N/2`. Combining
  with Thm 2.12: for fixed `K ≥ 2`, every input with `r ≤ (K−1)·log₂(N)` has
  `E ≤ N^K` — **effective, uniform in `d`, explicit constant**. Checked for
  `K ≤ 8`, `N` from `4` to `65536`.
  → **This halves the open case.** With `r ≤ N/2`, the unresolved band is
    `O(log N) ≤ r ≤ N/2`. Pinned by `test_cor219_leaves_a_gap` so a later
    improvement cannot quietly turn this into a claim about all inputs.
- [X] **The radicand-magnitude dichotomy (sharpens the search box).**
  If the `dᵢ` have `≥ 2^{r−1}/r` bits then `N ≥ 2^{r−1}`, so `2^{r−1} ≤ N` and
  Cor 2.19 applies already at `K = 2`. F₂-independence forces `r ≤ O(2^k)` for
  `k`-bit radicands, so `k ≥ log₂r − O(1)` — minimal radicands are both
  permitted and worst.
  → **The danger is exactly: small radicands, large `r`.** Pinned by
    `test_radicand_magnitude_trivialises_large_radican_ds`.
- [ ] **[RETRACTED] Thm 2.20 as claimed for `r ≤ 4` — the bound is FALSE.**
  The claim `|S| ≥ √D/((r+1)M)^{2^{r−2}−1}` for `2 ≤ r ≤ 4`, and the
  `α_r = r − 2^{r−2}` cutoff read off it, are **retracted**.
  → **The false step**: `min_{C_r∈ℤ} |C_r + a_j W| = a_j·‖W‖`. The left
    side is `dist(a_j W, Z)`, and `dist(a_j W, Z) ≠ a_j dist(W, Z)` unless `a_j = 1`
    — `a_j Z` is a different lattice from `Z`. Factoring a constant out of
    `dist(·,Z)` is fine; **minimising over a rescaled integer variable is not**.
    Refuted at `r=2` by `d=(2,3), C=(1,−1)`: `|S| = 0.3178 < √3 = 1.7321`.
  → **Caught by testing, not by re-reading.** An exhaustive check of the
    proposed `r=3` inequality over 910 000 certified instances gave **4 503
    violations**; instrumenting one exposed the identity. → Pinned by
    `test_RETRACTED_general_bound_is_false` and
    `test_the_false_step_is_the_lattice_mismatch`, so the retraction cannot be
    silently reversed.
  → **This is the THIRD error in the same step** in one session (dropped `M`;
    wrong field degree `2^{r−1}` vs `2^{r−2}`; the lattice mismatch). Each
    was asserted and then "corrected" **without running a test**. The lesson is
    recorded because it is a process failure, not a maths one: **in this repo,
    a bound is not claimed until a test tries to refute it and fails to.**
- [X] **Thm 2.20 as it survives — the constant is EFFECTIVE for `r ≤ 3`.**
  `r=2`: no constant at all (Thm 2.2, `log₂(1/|S|) ≤ N/2 + 3/2`, effective
  outright). `r=3`: `|S|M² ≥ 1/(2√D)`, so
  `log₂(1/c) ≤ 1.5 log₂D + 1 = O(N)`. The `r=3` step is
  `‖C₁a₁ + C₂a₂‖ = a₁·‖C₁ + C₂g‖` (valid: a fixed scalar factors out of a
  distance) with `g = √(d₂/d₁)` quadratic, then the integer norm
  `d₁|C₁²d₁ − C₂²d₂| ≥ d₁`. Verified: 2024 triples, 0 violations,
  3.18× margin. **`r ≥ 4` is OPEN.**
  → **Mechanism**: divide by `a₁ = √(d₁/d_r)`; the form becomes
    `C₁ + ΣCᵢgᵢ` with `gᵢ = √(dᵢ/d₁)`, and the `gᵢ` are `r−2` independent
    square roots of rationals, so the field has degree `n = 2^{r−2}`.
  → **Why it does not generalise**: at `r=3` the inner form has ONE quadratic
    generator so the norm is quadratic and collapses to `1/poly`. At `r=4` it
    has TWO independent generators; the quartic norm does not collapse, and the
    one hope — `a₁,a₂,a₃` in geometric progression so one `γ` captures all
    three — is excluded: `a₃/a₁ = (a₂/a₁)²` forces `d₃ = d₂`. Pinned by
    `test_r4_is_open_and_the_collapse_is_impossible`.
- [X] **`r=4` measurement, correctly labelled as an UPPER bound.** Exhaustive
  certified search, 3060 quadruples `d ≤ 89`, `|Cᵢ| ≤ 5`: `min c_S = 5.6×10⁻³` at
  `r=4` vs `1.7×10⁻⁶` at `r=3`. This is an **upper** bound on the infimum — so
  evidence that the `r=4` constant is not small. It is **not** a lower bound,
  not a theorem, and the `c(B)`-flattening test that certifies an infimum at
  `r=3` was **not run at `r=4`**. It does **not** rescue the retracted proof.
  → I had written this up as "confirms `α_4 = 0`", which was reading a
    conclusion into a number that is independent of the (false) argument.
    Pinned by `test_r4_measurement_is_an_upper_bound_only`.
- [X] **§2.11.10 — scoreboard**, so the pass's result is legible at a glance.
- [X] **§2.11.6 — the residual open statement, narrowed to one line.**
  For every `δ>0` and radicand tuple `d`, `c(d,δ) ≥ 2^{−poly(N)}` uniformly in `r`.
- [X] **Thm 3.2 audited in full; the write-up was CORRECT and I nearly
  "fixed" it anyway.** All four parts check out over 29 940 states
  (`a ≤ 500`, `n ≤ 60`), with the sharp constants landing on the values the
  proof predicts: `sup ρ·p⁴/a² = 0.0450842 = 1/(32·ln2)` (even-`n` branch,
  `D = −a`) and `sup ρ_n = 7.5129e-5` at `a=2, n=1`, matching the recorded
  `7.5e-5`. Part (iii)'s `0.205` is ~4.5× loose. **No fix needed to
  `PROOFS.md` or `f1_theorem.py`.**
  → **The audit itself failed twice, and both failures pointed at a correct
    theorem.** Recorded because the same pattern has now happened four times.
    1. `s = q√A/p` computed directly, then `1 − s`: at `p ~ 1e65` the
       difference is `~1e-130` and the float returned `s` as exactly `1.0`,
       so `ρ` read as `0` and `ρ·p⁴/a²` read as `0` or as garbage. Reported
       as "8 354 violations of `ρ > 0`". They were float artefacts —
       `ρ > 0 ⟺ D ≠ 0` is an **exact integer** test and `D ≠ 0` always holds.
    2. Stable `1 − s = D/(p(p+q√A))` fixed the `s` but I then used
       `ρ = log₂(1 + (1−s)²/(2s))`. The correct identity is
       `ρ = log₂(1 + (1−√s)²/(2√s))`; since `1−√s ≈ (1−s)/2`, the wrong form
       inflates `ρ` by **exactly 4**. That reported `sup ρ·p⁴/D² = 1/(8·ln2)`,
       against the documented `1/(32·ln2)` — i.e. it "found" a factor-4 error
       in Thm 3.2(iv) that does not exist. I got as far as editing the
       statement before re-deriving the algebra.
  → **Why the audit was worthless in that state**: both bugs were *plausible
    numbers*. `0.3499` looked like a counterexample and `8 354` looked like
    scale. Neither was a violation. A wrong audit is indistinguishable from a
    refutation unless the refutation is reproduced from the formula by hand.
    **New standing rule: before writing any retraction, re-derive the quantity
    symbolically and check the *constant* it converges to against the constant
    the proof predicts. A numeric "violation" whose magnitude is a clean
    rational multiple of a documented constant is a formula bug.**
  → **Repo gap this exposed, now fixed.** `test_proofs.py`'s F1 residual test
    computed `ρ` through the unstable `(1+s)/(2√s)` at `dps=100`, so past
    `n ≈ 20` it read `ρ = 0` and was corroborating nothing. Added
    `_f1_rho` (exact-`D` numerator + `log1p`, stable to `p ~ 1e190`) and
    three tests that pin the **constant** `1/(32·ln2)`, which is precisely
    what a factor-4 slip moves. Mutation-checked: reverting `_f1_rho` to the
    `(1−s)²/(2s)` form fails 6 subtests, so the guard is real.
    Note `log1p`, not `log(1+x)` — the addend is `t² ~ 1e-124` and
    `1 + t²` rounds back to exactly `1` even at `dps=800`.

### Open — from the §2.11 pass

- [ ] **Step 0 (literature, QUOTED): EHS Thm 4 does *not* give `E = O(N)`.**
  `PROOFS.md`'s table row reads `|E| ≥ Q^{−(n+1)}`, `E/N = O(N)`, "per-instance
  computable … at the cost of a checkable hypothesis". **Both halves are
  misleading** and the correction matters, because this is the only effective
  polynomial bound in the literature.
  - The paper was read in full (arXiv:2312.02057, HTML v1) and `Q` is a **free
    integer parameter**; only `λ₁ ≥ 1/Q^{1+1/(3n)}` is substantive, and Thm 8
    gets it for `Q ≥ Q₀` with `Q₀` depending only on the `a_i` — via the
    Subspace Theorem, hence **ineffective**.
  - Making `Q₀` explicit is where the cost lands. Bad `q` satisfy
    `q^{1+δ}·Π dist_Z(qβᵢ) < 1`, and the integer norm gives
    `dist_Z(qβᵢ) ≥ 1/(2qβᵢ+1)`, hence `q < (Πᵢ(2βᵢ+1))^{1/δ}`; with
    `δ = 1/(3n)` that is `exp(Θ(n²a))`, `a = max bitlen(aᵢ)`. Failure of the
    `λ₁` condition also needs `q ≥ Q^{1+δ}`, so `Q₀ = exp(Θ(n²a))` and the
    honest bound is **`E = O(n²a + n·m) = O(n·N)`** — a factor `n` worse than
    the target, and `n` can be `Θ(N)`.
  - **The obstruction is sharp**: the constant in `dist_Z(q√a) ≥ c/q` is
    `c ≈ 1/(2√a)`, i.e. `a/2` bits, amplified by `Θ(n²)`. The tempting
    `dist_Z(qβ) ≥ 1/q` is **false** — `dist_Z(7√2) = 0.101 < 1/7`. So this
    route provably does not close the conjecture.
  - **`λ₁` being computable is not the same as `Q₀` being computable.** Checking
    `λ₁` certifies one `Q`; the theorem needs the *threshold*. The repo's
    "checkable hypothesis" wording conflates the two.
  → **To do**: correct the table row and the `THE_MISSING_THEOREM.md` §7 item;
    record the derivation. This is a *literature correction*, so it belongs in
    `STATE_OF_THE_ART.md`, not in `PROOFS.md`.
- [ ] **Step 0 (possible slip in the published proof — verify before claiming).**
  From Minkowski II with `λ₁ ≥ 1/Q^{1+1/(3n)}`, the correct per-minima bound is
  `λ_i ≤ 1/Q^{n−1/(3n)}` for `i ≥ 2`, **not** the paper's `1/Q^{2/3}` — those
  agree only at `n=1`. The conclusion survives and in fact **improves**:
  `q ≤ Q^{1+1/(3n)}` rather than `Q^{n+1/3}`, hence
  `|E| ≥ 1/(2Q^{1+1/(3n)})`, better than their `1/Q^{n+1}` for `n ≥ 2`.
  → **To do**: re-check against the published PDF, and check the
  `p^Tx ≠ 0` step (it needs the argument that at most `n` of `n+1` independent
  lattice vectors can satisfy it). Only then is it worth asserting.
- [ ] **Step 1 (the decisive experiment): measure `c_r` for `r = 4…10`.**
  Nobody has looked at `r ≥ 4`. Since Thm 2.11 reduces the conjecture to the
  *constant*, the decay of `c_r` in `r` is exactly the thing to measure.
  - For each `r`, minimise `c(d) = min ‖C_r + Σ_{i<r}Cᵢ√(dᵢ/d_r)‖·max|C|^{r−1}`
    over radicand tuples and coefficient boxes.
  - **Two-stage is mandatory, not an optimisation**: full enumeration is
    `O(B^{r−1})` per tuple and dies at `r ≥ 4`. Screen coarse, refine the top.
  - `c(B)`-flattening as the infimum test (validated at `r=3`), and certify
    winners with `certified_interval_for_S` — **never** `signexact`'s
    `bracket_num` (sound, but coarser than `certified_interval_for_S`;
    §2.10.3 pass.
  - Report `log₂(1/c_r)` vs `r`, cross-checked against the pigeonhole rate
    `(r−1)/r` and EHS Thm 8's `2n`, so growth shows up as a violation of
    something concrete.
  → **Decision point**: `O(poly(r))` ⇒ conjecture plausible, and we learn what
    to aim the ST constant at. Super-`poly(r)` ⇒ **refutation**, and a major one.
- [ ] **`r = 3` uniform constant (Thm 2.11.5 → theorem).** Prove
  `c ≥ 2^{−C}` at `r=3` for all `d`. The vector lies on a conic
  (`β² = (d₂/d₁)α²`), so Davenport–Schmidt-style machinery and biquadratic
  regulator bounds are the natural tools. **Will not prove the conjecture**, but
  "Thm 2.9 is effective at `r=3` with an absolute constant" is a real result and
  upgrades the `r=3` guarantee from `3N` to `(2/3+δ)N`. **Timebox it** — a
  genuine open Diophantine problem; abandon cleanly if it does not move.
- [ ] **Where the danger actually is (proved, cheap, do next).** Thm 2.12 gives
  `E ≤ 2^{r−1}(N + log₂3r)`, and `N ≥ Σ log₂dᵢ`, while F₂-independence forces
  `r ≲ 2^k` for `k`-bit radicands. So radicands with `≳ 2^{r−1}/r` bits give
  `2^{r−1} ≤ N^{O(1)}` and the conjecture holds trivially; the danger is exactly
  **large `r` with minimal radicands** (`d_i = O(poly(r))`,
  `N ≈ r(log r + log C)`). Cheap to prove, and it shrinks the search box to a
  sharply-described region.

### Open — after the theory pass

- [ ] **The third family of argument.** Thm 2.17 rules out every norm and
  Thm 2.9's Subspace Theorem is ineffective, so nothing currently in reach
  closes the band `O(log N) ≤ r ≤ N/2` with small radicands. Any progress
  here is a genuinely new idea, not a refinement of an existing route. The
  concrete open target remains a **uniform** `c(d,δ) ≥ 2^{−poly(N)}` in Thm 2.9.
- [ ] **`r = 3` uniform constant (Thm 2.11.5 → theorem).** Unchanged; see
  above. Timeboxed.

### Explicitly not worth doing (recorded so it is not re-proposed)

| Idea | Why not |
|---|---|
| Better norm, cleverer conjugate choice | **Thm 2.17** — the exponent is `deg S = 2^r`; impossible, not merely hard |
| Any norm-based route at all (lattice norms, subfield norms, tower norms) | Thm 2.17 — same reason; the whole family is dead |
| Squaring / degree reduction | Thm 2.14 — exact wash; `R = ΣCᵢ²dᵢ > 0` always |
| More `r ≤ 3` constant numerics | Saturated at 19–22 bits (§2.11.5); will not move |
| p-adic / S-unit route before Steps 0–1 report | No target yet; the repo already flags it as a lead, not a result |
| Trying to prove `E ≤ poly(N)` directly | It **is** the open problem (EHS Problem 33) |

---

---

## Research program — the uniform separation constant (NEW)

**Crucial re-diagnosis.** The effective norm bound of Paper B gives, at every
*fixed* rank `r`, `E ≤ (2^{r-1}-1)·((3/2)N + log₂N)` — a constant since `r` is
fixed. So the SRS conjecture `E ≤ poly(N)` is ALREADY effective at every fixed
rank, in particular at r=2 and r=3.  There is **no conjecture gap at r=3 from
the norm side**.  The only genuinely open case is the one where `r` is not
fixed but grows with `N`: then `2^{r-1}` is exponential in `N`, and one needs
the Subspace-Theorem constant to satisfy `log₂(1/c) ≤ poly(N)`.  That uniform
constant is EHS Problem 33.

So the r=3 quartic norm computation is **not** the entry point.  It only
re-proves the known effective rank-3 bound `|u| ≥ d₃²/(2M)³` and shows where
it is sharp.  The open problem is specifically:

  **(Missing)**  `log₂(1/c(d,δ)) ≤ poly(N)`  uniformly over `r`, radicands, `C`.

Reframed, the correct ladder is:

- [ ] **R1.** State honestly in `Docs/` that fixed-r SRS is already done; the
  open part is uniform-in-`r`.  (Update Papers F/8 to say this explicitly.)

- [ ] **R2.** Make Paper 8's reduction *exact*: the Subspace constant `c`
  fails to be `2^{-poly(N)}` iff there is a sequence of rank-`r`, `r→∞`-style
  instances with `|L| < 2^{-N^k}` for all `k`; this is a statement about
  simultaneous sharp-approximation in multiquadratic fields as `r` grows.
  Write a sharp formal statement (no new proof claimed).

- [ ] **R3.** The missing estimate is a *trilinear-in-logarithms* lower bound
  in the units of a rank-`r` multiquadratic field with logarithms
  `log|ε₁|,…,log|ε_{2^{r-1}-1}|`.  Existing effective bounds on linear forms in
  logarithms are exponential in the heights, NOT `≤ poly(N)`, so they do NOT
  immediately give the SRS constant.  Record this as the measured gap.

- [ ] **R4.** Do NOT invoke generic Baker substitution as a proof.  Instead,
  look for special structure: in a multiquadratic field the units decompose
  and the three/r−1 independent units have a Kronecker limit / regulator
  relation that might force an effective bound.  This is the first theory
  project.

- [ ] **R5.** Decide whether the right object to attack is (a) an effective
  simultaneous-approximation theorem in multiquadratic fields, (b) an
  effective bound for a specific norm form, or (c) a new analytic trick.
  Write the chosen subproblem as one numbered theorem to prove.

- [X] **R6.** `Docs/09_RANK3_NORM_FORM.md` (write): derive `F(C)` exactly,
  clear denominators to `d₃²F(C) ∈ Z`, and record the two facts —
  `F(C)=0` is impossible on the interior (Paper A), hence `|F| ≥ 1/d₃²`,
  hence `|u| ≥ d₃²/(2M)³`, which IS the effective rank-3 bound.  This closes
  r=3 as a *norm* statement and is the local reason the SRS problem is about
  growing `r`, not about r=3.

- [X] **R7.** Computational anchor: for small triples, the histogram of
  interior minimum of `d₃²|F(C)|`.  Confirms the bound and the sharpness of
  `2^{r-1}-1` replacement question.

### Main theorem target (single, honest)
There is one open theorem.  NOT "solve r=3" (that is done effectively) but:

   For every `δ>0`, the Subspace-Theorem constant `c(d,δ)` for SRS satisfies
   `log₂(1/c) ≤ A·N^K` for absolute `A,K`.

## Research program v3 — escape the `2^r` conjugate explosion (NEW)

Fixed-rank `r` is effective (Paper B, and explicitly `E ≤ 6N+5` at `r=3`,
Paper 9). The uniform conjecture is precisely the Subspace constant:

  **Open.** `log₂(1/c(d,δ)) ≤ poly(N)` as `r` grows.

Milestone ladder:

- [X] **M1.** Negative result (write it): progressive relative norms
  `N_{K_i/K_{i-1}}` over a quadratic tower telescope to the full norm, so the
  exponent `2^{r-1}-1` is unchanged. Any effective escape must use
  simultaneous-approximation structure, not pure norms. (Draft: `Docs/10_TOWER_NORM_NEGATIVE.md`)

- [X] **M2.** Enemy generator run (`enemy_gen.py`); profile recorded in `Docs/11`.
  Output: the worst-case geometry is norm-tight on the consecutive-gap family (one tiny conjugate, rest generic); the simultaneous/two-small-embedding regime is the next object.

- [ ] **M3.** Structured growing-`r` family: try `d_i = p_i` (distinct
  primes). Attempt to prove `|S| ≥ 2^{-N^{O(1)}}` uniformly in the number of
  terms. Success ⟹ mechanism; failure ⟹ sharp counterexample.

- [~] **M4.** Two-embedding lemma PROVED for one direction, see `Docs/12_M4_TWO_EMBEDDING_LEMMA.md`.
  - M4.1 done: if two sign embeddings of `S` are both `≤ ε`, then `ε` is bounded below
    by the rank-`r_J` norm inverse for the corresponding sub-form (where `J` is their sign-difference support).
    Consequence used downstream: on the interior (`C_i ≠ 0`), you cannot have two embeddings
    both smaller than the lower-rank norm threshold; the smallest embeddings are separated
    by at least the rank-2/Pell scale. This validates AND delimits the divide
    "one small embedding => norm-tight" / "two small => forced sub-relation".
  - M4.2 OPEN: two-small-embedding *generator* (search pairs `min_C max{|L_σ|,|L_τ|}`) and the
    accompanying quantitative rank-of-difference-matrix experiment.

- [X] **M4-theory branch (Docs/13).** Proved: rank of the difference matrix is the rank of the nonzero column-pattern vectors, <= min(k-1,m); the k smallness conditions are equivalent to a small projection P^T B of the block-sums B_c, NOT to the smallness of any single B_c; and in the full-rank case q=m the only way to keep all differences small is a full-m integer norm whose threshold is >= the full-norm threshold.  Net: sign-difference rank does not compress the exponent 2^{r-1}-1; the missing input is an effective simultaneous-approximation statement on the real vector (sqrt(d_i)).  M4.2 computational and M3 remain deferred.
## M5 — the exact external theorem needed (NEW, replaces experiment loop)

We stop grinding norm/tower/sign manipulations.  Paper 13 fixed the boundary:
the missing exponent is an effective simultaneous-approximation statement on
`(√d_i)` that is invisible to the field norm.  This is now a literature
translation problem.

- [ ] **M5.1 Exact required theorem.**  The weakest standalone statement that
  would prove the conjecture.  Target form: for independent squarefree
  `d₁,…,d_m` and integer `C₁,…,C_m`, `S=Σ C_i√d_i ≠ 0`, prove
  `|S| ≥ 2^{-A N^K}` for `N = Σ(bitlen|C_i|+bitlen d_i)`, absolute `A,K`.
  No `δ`, no subspace exceptional-set output: a single uniform quantitative
  lower bound.  Write this in `Docs/14_M5_EXTERNAL_THEOREM.md`.

- [ ] **M5.2 Translate the literature.**  For each candidate effective result
  (EHS SoCG'24, effective Subspace/multiquadratic analogues, Baker–Wüstholz
  linear forms in logs, quantitative Jarník-type simultaneous approximation,
  unit-equation height bounds), compute the resulting exponent and constant as
  a function of `N, r, m, d_i, C_i`.  In particular: EHS gives
  `|S| ≥ γ·(m·max|C_i|)^{-2m}` with `γ=γ(d)` not explicit; the hence-so
  `E ≤ 2m·log₂(m·max|C_i|)+log₂(1/γ)` is polynomial in `N` only once
  `log₂(1/γ) ≤ poly(N)` — which is exactly the unproved uniform constant.
  Record the precise parameter that is missing in `Docs/14`.

- [ ] **M5.3 Multiquadratic specialization.**  Ask whether the Galois group
  `(Z/2Z)^r` allows any of the general estimates to collapse.  Two known
  `r=2` corners (rank-2: Markov/CF is effective; the full norm bound at r=2)
  do NOT obviously lift: the q=r rank-saturation case again reproduces the
  full-norm exponent, and generic `{√d_i}` give dense unit groups whose
  logarithm relations are exactly the trilinear form of Paper 8.  State as
  an honest conjecture: `(Z/2Z)^r` structure does **not** by itself lower the
  effective exponent, so the only route to `poly(N)` is a new effective
  Subspace-type theorem — write down what that theorem must say about the
  simultaneous-approximation constant of a multiquadratic vector.

- [X] **M5.4 Stop condition.**  Do NOT claim the conjecture is proved until
  M5.1's exponent `N^K` is actually established in a self-contained proof.
  Before that, "SRS reduces to effective Subspace" is the correct label.



---

---

## M7 — Uniformity audit: what information beyond the norm is available? (NEW)

Target statement: for every `r`, every squarefree tuple `d`, and every `C ≠ 0`
with `S ≠ 0`, find an *effective* bound

```
log₂(1/c(d,δ))  ≤  poly(r, N)
```

for the Subspace-Theorem constant in

```
|L(C)|  >  c(d,δ)·||C||^{-(r-1+δ)} .
```

That statement is exactly "uniform effective SRS".  To find it we first
classify what extra information a putative bound could use; every known
source of information is currently in the closed column below.  A new proof
must read structure in the open column.

Closed routes (each is a documented non-solution in `Docs/` or `TODO`):

| Route                                  | Status / why it stops                                   |
|----------------------------------------|----------------------------------------------------------|
| Full norm / conjugate multiplication  | `Docs/B`, exponential `2^{r−1}−1` in `r`.               |
| Relative norm tower                    | `Docs/10`, telescopes to left-same full norm.             |
| One-small-embedding pigeonhole         | `enemy_gen.py`, norm-tight enemy `Docs/11`.              |
| Two-small-embedding lemma              | `Docs/12`, valid but insufficient uniformly.             |
| Sign-pattern difference rank descent   | `Docs/13`, does not compress the conjugate count.         |
| Generic Subspace Theorem               | `Docs/C`, correct exponent, constant not computable.      |
| Generic Baker/Matveev                  | `Docs/16`, `𝔇 ≤ 2^{O(N)}` scaling, no polynomial in N.   |
| Rank-3 quartic/unit reduction          | `Docs/15, Docs/16`, fixed-rank bound only.               |

Three open buckets:

1. **Algebraic / integrality beyond the norm.**  A sharper than-full-norm
   invariant of the *same* sum `S` (structured divisors, discriminants,
   resultants, ideal structure of `d`'s).

2. **Simultaneous approximation beyond a single form.**  Lower bounds on
   simultaneous closeness of several `√(d_i/d_r)` to rationals, with a
   constant not equal to the Subspace constant — the two-small-embedding
   obstruction from `Docs/12` already constrains what the sign-difference
   matrix can give.

3. **Special-field mechanism.**  A genuinely multiquadratic fact about
   `Q(√d₁,…,√d_r)` and the coefficient vector `C`, beyond the
   `(Z/2Z)^r`-fold embedding symmetry exploited by the norm pairing.

M7.1 (DONE).  `Docs/19_M7_1_ENEMY_AUDIT.md`.  Two-small search after
antipodal de-redundancy: in all tested prime/small-C cases, any pair of
simultaneously small embeddings has an F2-independent difference support.
No rank-deficient support was found.  A three-way empirical split is
visible: local/Pell (|J|=1), short independent (|J|=2..3), and a residual
full-support regime (|J|=r-1, r_J=|J|) that recurs at nontrivial rate.

M7.2 (DONE, v1.1 CORRECTED).  *Full-rank residual regime audit.*  Executed:
9104 residuals across r=3..7, B=2..8.  **The v1.0 conclusion of Docs/20 was
WRONG and is retracted:** it claimed the residual regime was a "near-antipodal
artifact" with median `E_Y/E = 1.000`.  That came from a buggy filter — the
script zipped a genuine-case list against a differently-filtered list when
computing medians, and the X=0 classification used float comparisons.
Corrected: 3344 of 9104 rows satisfy hypothesis (H) (nonzero coefficient on
both the agreement and difference sets), of which 1616 are genuine
(|X| >= |Y|/2).  |A| = 2 and supp(X) = 2 in 3344/3344 rows.  See Docs/20
(RETRACTION notice) and Docs/21.

M7.3 (DONE — NEW THEOREM).  *Two-small refinement.*  Proved and verified;
this is the first genuine improvement over the norm exponent that is not a
norm identity.  If two embeddings are both <= eps and (H) holds, then

```
log2(1/eps) <= (2^{min(k_A,k_J)}-1) * (1.5*log2 M + log2 max(s_A,s_J) + 1)
```

with k_A, k_J the F2-ranks of the radicands on the agreement / difference
sets and s_A, s_J their sizes.  Since k_A + k_J <= m, min <= floor(m/2), so
the exponent drops from `2^r - 1` to `2^{floor(m/2)} - 1`.  Verified: 0
violations in 3344 residual rows and 5970 fresh random cases (60-digit
mpmath).  Proof is elementary; no Subspace-Theorem constant, Baker bound, or
regulator is used.  See `Docs/21_TWO_SMALL_REFINEMENT.md`.

NOTE an error was caught and fixed during this pass: the first draft used
`max(k_A,k_J)` in the prefactor, but the conjugate bound is by support SIZE
`|A|,|J|`, not rank.  Corrected before any claim was recorded.

Consequence: because the observed residuals all have |A| = 2, they satisfy
`log2(1/eps) <= 3*(1.5 log M + log m + 1) = O(log N)` — polynomially bounded.
So every two-small configuration found is provably benign.

M7.4 (DONE for the ANALYSIS; the open target is now sharply posed).
*The isolated-embedding regime.*  Wrote `Docs/22_ISOLATED_REGIME.md`.

(a) Exact identities for the conjugates, all proven: f is AFFINE LINEAR on
the hypercube (only m of the 2^m Walsh coefficients are nonzero); trace 0;
Parseval sum_tau f(tau)^2 = 2^m sum_i C_i^2 d_i; the flip identity
f(tau) f(tau^e_i) = X_i(tau)^2 - b_i^2 lying in the degree-2^{r-1}
subfield; and a conjugate-separation bound
|f(tau)-f(tau')| >= 2 (|T| M^{3/2})^{-(2^{k_T}-1)}.
Corollary: conjugates are polynomially rigid up to Hamming distance
O(log N) and exponentially soft beyond it.

(b) NO-GO THEOREM (the important negative result).  Every Liouville/norm
bound on |S| has the form |S| >= prod_{tau != 1} |f(tau)|^{-1}.  The
isolation hypothesis gives LOWER bounds on those factors, i.e. the wrong
sign.  Hence isolation is arithmetically INERT: it cannot improve any norm
or Fourier argument.  Corollary: the norm bound's worst case IS the maximally
isolated configuration, so M7.3 (two-small) and the norm bound (isolated)
are complementary and jointly exhaustive — there is no third regime.

(c) Verified by exact multiquadratic ring arithmetic: the certified witness
d=(2,3,13), C=(-10,-98,51) has N(S) = 1 — it is a UNIT — with
E1 = 21.7517644049, E2 = -4.8219281093, E2/E1 = -0.2217, so it is ISOLATED.
=> M7.3 does NOT improve the sharp constant; the extremal witnesses live
exactly in the regime M7.3 cannot reach.

(d) FINAL OPEN TARGET (unchanged in substance, now the ONLY item):
prove log2(1/c(d,delta)) <= poly(N) for the multiquadratic Subspace-Theorem
constant, equivalently make the isolated case effective.  Do NOT re-propose
pure-norm, norm-tower, sign-pattern, or Fourier/Hadamard refinements —
(b) proves they are futile.

(e) BOOKKEEPING RESOLVED — rate and witness are SEPARATE claims.  TODO had
recorded this witness as "saturating limsup E/N = (r-1)/r with E/N ~ 1.2795",
but 1.2795 = 21.7518/17 uses coefficient bits only; under the standard
convention N = 25 the ratio is 0.8701, and neither equals the asymptotic 2/3
(the finite r=3 instance EXCEEDS it).  The dossier now states:

  Claim A (asymptotic rate): limsup_{r->inf} E/N = (r-1)/r, via Subspace
    Theorem + pigeonhole (Docs/02 B.5, Docs/03).  A statement about a FAMILY
    as r -> infinity.

  Claim B (finite witness): d=(2,3,13), C=(-10,-98,51) is a certified UNIT
    (N=1) with one isolated tiny embedding, E1=21.7518, |T2|=28.28,
    E2/E1=-0.222.  It is an extremal finite example for the NORM MECHANISM,
    NOT evidence that a single r=3 instance saturates the asymptotic rate,
    and is not offered as such.  Its role here is only to certify that the
    isolated regime is realised by a unit.

  Consistency check: ST bound with no additive constant gives
  (r-1)(1.5 log2 M + log2 m) = 23.0 >= E1 = 21.75, so the witness is
  consistent with the r-1 exponent under either N convention.

The old phrasing "saturating (r-1)/r" is hereby RETRACTED as a conflation.

L1. LITERATURE PHASE (DONE, frozen as frontier).  Wrote `Docs/23_EFFECTIVENESS_LITERATURE.md`.
Findings, all verified against our certified data:

1. OUR PROBLEM = the classical sum-of-square-roots SEPARATION BOUND.  The
   norm bound of Docs/02 (exponent 2^r - 1) is literally the Burnikel-
   Fleischer-Mehlhorn-Schirra (2000) bound.  Linear independence =
   Besicovitch (1940); poly-time zero test = Bloemer (1991/1998).

2. THE POLYNOMIAL EXPONENT IS ALREADY SOLVED IN THE LITERATURE.
   Eisenbrand-Haeberle-Singer (arXiv:2312.02057, SoCG 2024, Thm 8):
   |E| >= gamma * (n*||x||_inf)^(-2n), gamma depending only on sqrt(a_i).
   In our variables: E <= O(N^2) + log2(1/gamma).  Their Remark 9 sharpens
   2n to n+eps.  => Docs/03 and Docs/22 Theorem D are CONFIRMED FROM
   OUTSIDE: the exponent half is a published theorem, and the ONLY
   obstruction is log2(1/gamma).

3. BEST EXPLICIT CONSTANT: Aymone-Figueredo-Iyer-Tafula (arXiv:2609.14161,
   Sept 2026) give a FULLY EXPLICIT bound with exponent 2^(K-1)-1 and prefactor
   exp((2^(K-1)-1)/K).  Still exponential in rank, so not SRS.  BUT our M7.3
   gives 2^(floor(m/2))-1 in the two-small regime, which STRICTLY DOMINATES
   their universal 2^(r-1)-1 there, with an elementary proof and no prefactor
   loss.  Apparently not in the literature.

4. CROSS-CHECK on certified data (40 digits): norm / AFIT / EHS bounds all
   hold on the r=3,4,5 witnesses.  The r=3 unit witness caps
   gamma <= 1.83e8 for d=(2,3,13).

5. EFFECTIVITY STATE OF PLAY: explicit t (Schmidt 1989, Evertse-Ferretti
   2013, t <= (2d)^(26n) eps^-2) and an effectively determinable exceptional
   subspace U_0 (Vojta 1989, Faltings-Wustholz 1994), but the terminal finite
   solution set is NOT effectively bounded.  Koleda (arXiv:2411.01247, Res.
   Number Theory 2026) calls determining the maximal height of the
   subspaces "one of the most fundamental open problems in Diophantine
   Analysis".  EHS section 4 verbatim: explicit bounds are unavailable
   "still when all algebraic numbers are square roots of integers".

6. TRAP CORRECTED: the natural target gamma >= Disc(K)^-O(1) is WRONG,
   because log2 Disc(K) ~ 2^(r-1)*sum log2 d_i is EXPONENTIAL in r
   (measured ratios 5.9 / 11.1 / 20.3 for r=3,4,5).  It would only give
   log2(1/gamma) = O(2^r * N), still useless.

7. SHARP TARGET (final form):
       gamma >= (d_1*...*d_m)^(-C)   for an absolute C,
   equivalently log2(1/gamma) <= C*sum log2 d_i <= C*N.
   This is the only formulation that can work, since sum log2 d_i <= N by the
   definition of the input size.  In EHS's terms: an effective upper bound on
   the q in N_+ with q^(1+delta) * prod_i dist_Z(q sqrt(a_i)) < 1, for
   delta = 1/poly(n).

L2. Q-SET ATTACK (DONE).  Wrote `Docs/24_QSET_ATTACK.md`.  Dedicated attack
on the EHS exceptional q-set.  Results, all verified by exhaustive exact
integer arithmetic:

(a) EXACT REFORMULATION: with p_i = round(q sqrt a_i) and n_i = p_i^2 - a_i q^2
    in Z (nonzero), the EHS condition  q^(1+d) prod_i dist_Z(q sqrt a_i) < 1  is
    EQUIVALENT to   prod_i |n_i|  <~  2^n sqrt(prod a_i) * q^(n-1-d).
    So the exceptional set = common denominators where ALL n forms
    p^2 - a_i q^2 are simultaneously small in the product sense.

(b) THEOREM 24.1 (sharp n=1 boundary).  For n=1 the exceptional set is finite
    with the EXPLICIT bound q < (2 sqrt a + 1)^(1/d) -- fully effective.
    For n >= 2, per-coordinate Liouville dist >= c/(q sqrt a) forces only
    q^(n-1-d) > (lower bound), i.e. a LOWER bound on q: NO CONSTRAINT.
    Reason: one quadratic irrational is separated from Z at rate 1/q, which
    already beats the q^(-1-d) threshold; n of them give only q^-n, far worse.
    The threshold is reachable only by RARE SIMULTANEOUS near-convergents,
    which no per-coordinate estimate can detect.

(c) RECORDED SIGN ERROR (I made it, numerics caught it): I first "proved"
    q* < (4^n sqrt(prod a))^(1/(n-1-d)) -- an explicit bound that would have
    finished the problem.  The error was inverting an inequality with a
    NEGATIVE exponent.  Scan refuted it immediately: predicted q < 81.6 for
    (sqrt2,sqrt3), actual q reaches 29,354,524.  This is Docs/22 Theorem C
    one level up: a lower bound on a product plus a requirement that the
    product be SMALL is arithmetically inert.

(d) THEOREM 24.3 (Pell LOWER bound, NEW).  Every convergent denominator q of
    sqrt(a_i) with q < 2^(3n(n-1)) is exceptional; every Pell denominator
    (p^2-a_i q^2 = +-1) with q < (2^n sqrt(a_i))^(3n) is exceptional.
    Verified with ZERO exceptions for n=2,3,4 and 4 radicand sets.
    Explains why brute-force scans at n=3 were still finding new q at 4e7 --
    the true q* there exceeds 1.3e9.

(e) THEOREM 24.4 -- RETRACTS THE Docs/23 TARGET.  No absolute C exists with
    gamma >= (d_1...d_m)^(-C): Thm 24.3 forces log2(1/gamma) = Omega(n^3)
    while sum log2 d_i >= log2(n!) = Theta(n log n); ratio diverges (measured
    18, 28, 39, 52, 66, 99, 136 for n=2..10).  Docs/23 section 23.7 is marked
    RETRACTED in place.

(f) THEOREM 24.5 (corrected, TIGHT target).  It suffices to prove
        log2 q* = O( n^2 + n log2(max d_i) + n log n ),
    equivalently q* <= exp(O(N^2)), which gives SRS WITH EXPONENT 3.
    The n^2 scale is FORCED by (d), so the target is tight up to constants and
    any proposal with log2 q* = o(n^2) is already refuted.

(g) Empirical structure, (sqrt2,sqrt3) delta=1/6, exhaustive scan to 4e7:
    |E| = 121, q* = 29,354,524 (interior, stabilised).  Convergent denominator
    of sqrt2: 18/121; of sqrt3: 23/121; of AT LEAST ONE: 40/121; of BOTH: 1/121.
    So ~2/3 are genuine simultaneous near-convergents explained by NEITHER
    individual continued fraction.  (NOTE: an earlier CF routine in my script
    was buggy -- negative denominator -- and gave a FALSE "not convergents"
    conclusion; fixed and re-verified.)

(h) STILL NO UPPER BOUND.  Docs/24 supplies lower bounds and a refutation only.
    The effectivity problem itself is untouched and remains open.

L3. RATIO / DEFECT-RIGIDITY ROUTE (DONE -- REFUTED).  Wrote
`Docs/25_RATIO_ROUTE.md`.  Rigorous treatment of the proposed attack via
`a_j p_i^2 - a_i p_j^2 = a_j u_i - a_i u_j` and pairwise ratio approximations
to sqrt(a_i/a_j).  VERDICT: REFUTED, theoretically and empirically.

(a) THEOREM 25.1 (proved).  Cross-defect rigidity:  c_ij := a_j u_i - a_i u_j
    != 0 for i != j.  Proof: c_ij = 0 => (p_i/p_j)^2 = a_i/a_j; a prime p|a_i,
    p not|a_j gives even vs odd p-adic valuation.  Equivalently a_i/a_j is
    never a rational square for distinct squarefree.  So u not in span_Q(a).

(b) THEOREM 25.2 (proved) -- THE DECISIVE NEGATIVE.  With D = a_i a_j,
      |p_i/p_j - sqrt(a_i/a_j)| = |c_ij| / ( p_j (a_j p_i + p_j sqrt D) )
                                >= 1/( C_ij q^2 ),
    C_ij depending ONLY on a_i,a_j -- NOT on q, NOT on the defects.  Since
    |c_ij| >= 1 the ratio approximation is NEVER better than the plain
    quadratic bound.  So the ratio route is a RESTATEMENT of the coordinatewise
    quadratic estimate and carries ZERO new information.  The proposal's
    central premise is provably false.  (Also: sqrt(a_i/a_j) is quadratic, so
    no Roth-type slack exists either.)

(c) RIGIDITY HYPOTHESIS REFUTED EMPIRICALLY.  On the 121 certified exceptional
    q for (sqrt2,sqrt3), delta=1/6:
       zero cross-defects: 0/121  (Thm 25.1 holds -- the trivial part)
       |c_12|/(|u_1|+|u_2|): exceptional median 2.1778, mean 2.2515
                             random control  median 2.1667, mean 2.2367
       |c_12| = 1: 0/121 ;  |c_12| <= 10: 9/121 ;  <10% of sum: 1/121
    Statistically INDISTINGUISHABLE from random.  The defect vector carries no
    detectable extra algebraic rigidity.

(d) THEOREM 25.4 (identity, but vacuous).  (a_j u_i + a_i u_j)^2 =
    c_ij^2 + 4 a_i a_j u_i u_j  -- verified exactly.  But it is an IDENTITY
    (reparametrisation, not a constraint): all (u_i,u_j) satisfy it, and the
    count O(P log P), P ~ q^(1-2d), is unreduced.  No gain.

(e) THE FIVE PROPOSED TARGETS, each addressed:
    1. pairwise coprime  -- ALREADY REFUTED by Docs/24 Thm 24.3 (Pell works
       for ANY radicand set; verified a=(2,3,5,7) pairwise coprime has 15007
       exceptional q <= 2e7, still growing).
    2. first n primes    -- same, refuted.
    3. bounded A         -- reduces to a finite but NON-DECIDABLE question: no
       upper bound on q* means no way to certify a finite search.  No gain.
    4. n=2 or n=3       -- the genuine first case, and the precise obstruction
       is now identified: reduction theory of indefinite binary quadratic forms
       IS effective and DOES control c_12 = a_2 p_1^2 - a_1 p_2^2, but EHS
       needs u_1 u_2.  By Thm 25.2 these are unrelated.  THAT is why the
       classical effective machinery stops here.
    5. soft defects |u_i| <= q^theta, theta<1 -- SURVIVES.  Measured
       theta(q)=max(|u_1|,|u_2|)/q: median 0.196; the 17 points with theta>=0.9
       are all SMALL q (<=5741); the large-q population satisfies |u_i| <= q^0.9.
       Note |u_i| = Theta(q), not o(q).

(f) VERDICT: the exceptional q-set carries NO algebraic rigidity.  The
    phenomenon is METRIC, not algebraic -- genuine multidimensional
    Diophantine approximation (consistent with Docs/24 section 24.8: 2/3 of
    exceptional q are near-convergents for NEITHER coordinate).

(g) SECOND INERTNESS VERDICT.  Same structure as Docs/22 Thm C, one level up:
       Docs/22: "other conjugates are large"     -> no norm/Fourier gain
       Docs/25: "cross-defects small / ratios well-approximated" -> no gain
    A pure algebraic attack on the EHS constant is now closed from two
    independent directions.

(h) WHAT REMAINS (the only live route).  Effectivity in MULTIDIMENSIONAL
    DIOPHANTINE APPROXIMATION for algebraic points: an effective positive
    lower bound on the multiplicative Markov constant
       mu(a_1..a_n) = liminf_q q^(n+1) prod_i dist_Z(q sqrt(a_i)) > 0
    at the scale log2 q* = O(N^2) (Docs/24 Thm 24.5).  The Pell family is the
    hard test and already forces the n^2 scale.  NO such effective bound is
    known; open.

L4. SRS MARKOV FUNCTION (DONE).  Wrote `Docs/26_MARKOV_FUNCTION.md`.  Built the
custom function and six exact lemmas.  The function is NOT another invariant
of u_i -- it measures joint metric rarity directly.

  rho_a(q) := 4^n prod_i ||q sqrt(a_i)||,   M_a(Q) := min_{q<=Q} rho,
  G_a(Q)  := min_{Q<q<=2Q} rho.
Normalization: for generic q, ||q sqrt a_i|| ~ U[0,1/2] mean 1/4, so
rho ~ 4^n (1/4)^n = 1.  So rho = 1 is TYPICAL and an exception is a
SUPPRESSION by a factor ~ q^(1+delta)/4^n -- roughly one power of q.
Exception condition becomes exactly  rho_a(q) < 4^n q^(-1-delta).

(a) LEMMA 1 (exact).  rho_a(q) = (prod_i |u_i|)/sqrt(prod a_i) * (1+O(1/q)).
    So rho IS the rescaled defect product -- scale-free, unlike Docs/25's
    cross-defects which depend on the choice of algebraic form.

(b) LEMMA 2 (CORRECTS A MISREADING).  Markov positivity min_q rho > 0 is
    ELEMENTARY (Liouville on each sqrt(a_i) alone); NO Subspace Theorem needed.
    The ST only adds the far stronger rho >= 4^n q^(-1-delta) for all large q.

(c) LEMMA 3 (REFUTED -- third inertness trap).  "Exception => q is a
    convergent denominator of some sqrt(a_i)" is FALSE: only 40/121 of the
    certified exceptions are, and only 1/121 of both.  Why it fails: to force a
    contradiction the floor must EXCEED the threshold, but
    q^(-n) < q^(-1-delta) holds for all n>=2, q>=1 -- the floor is
    automatically compatible.  Same wrong-direction failure as Docs/22 Thm C
    and Docs/24 Thm 24.1(ii).

(d) LEMMA 4 (THE TAIL IS FLOOR EVENTS -- the operative reduction, verified).
    For q a convergent denominator of sqrt(a_i), |u_i| = O(sqrt a_i) is bounded,
    so ||q sqrt a_i|| ~ 1/(2 q sqrt a_i) and the exception condition reduces
    EXACTLY to the MILD condition  prod_{j != i} ||q sqrt a_j|| < 2 sqrt(a_i)
    q^(-delta).  Scanning the (computable, geometrically growing) convergent
    denominators to 1e20:
       Q=1e13:  79 candidates, 42 hits, max q = 1,104,728,155,436
       Q=1e16:  97 candidates, 43 hits, max q = 1,746,860,020,068,409
       Q=1e20: 122 candidates, 43 hits, max q = 1,746,860,020,068,409 (stable)
       mechanism check: 42/42 satisfied; 33/43 have min|u_i| <= 1 (true Pell).
    ==> q* >= 1.75e15 for (sqrt2,sqrt3), delta=1/6.  SEARCH SPACE COLLAPSES
    from all q <= Q to O(log Q) candidates.

(e) MAJOR CORRECTION to Docs/24 section 24.8: its q* = 29,354,524 was WRONG --
    the brute-force scan to 4e7 was ~8 orders of magnitude too shallow and the
    set had NOT stabilised.  Corrected in place with a pointer to Docs/26.

(f) LEMMA 5 (metric tools predict the wrong magnitude).  The naive
    volume/ergodic count predicts ~0.76 exceptions up to Q=4e7 for
    (sqrt2,sqrt3); OBSERVED 121 -- under-prediction by ~160x.  So
    discrepancy / volumes-of-numbers / ergodic reasoning CANNOT deliver the
    bound: it is quantitatively wrong by >2 orders of magnitude, because the
    driver is algebraic (Pell/convergent) structure, not equidistribution.

(g) OPEN CORE (M7.7) -- now a single sharp statement:
        Bound effectively the convergent denominators q of sqrt(a_1) for which
        prod_{j>=2} ||q sqrt(a_j)||  <  2 sqrt(a_1) * q^(-1/(3n)).
    Required scale log2 q* = O(N^2) (Docs/24 Thm 24.5).  Why this is the right
    core: weak side exponents (1/(3n)) so no CF/Roth control applies; q ranges
    over a computable geometrically SPARSE set so it is not a density question;
    it is exactly the Pell configuration, so any correct theorem must survive
    Docs/24 Thm 24.3; and it is killed by neither Docs/25 nor Docs/22 Thm C.

L5. PELL-INDEX VARIABLE (DONE).  Wrote `Docs/27_PELL_INDEX.md`.  Changed
variables from q to the CF/Pell index k.  Results:

(a) LEMMA 27.1 (verified).  For the Pell sub-sequence p_k^2-a_1 q_k^2 = +-1 and
    fundamental unit eps:  q_k = (eps^k - eps^-k)/(2 sqrt a_1) EXACTLY, hence
        q_k sqrt(a_j) = beta_j (eps^k - eps^-k)/2,  beta_j = sqrt(a_j/a_1)/2
    and   ||q_k sqrt(a_j)|| = ||beta_j eps^k|| + O(eps^-k).
    beta_j is FIXED and algebraic (in Q(sqrt a_1, sqrt a_j)).  Verified
    numerically to an index shift of one.  So the object is exactly the
    distance to Z of a MULTIPLICATIVE ORBIT -- S-unit / recurrence-mod-1
    territory, which is the genuinely new coordinate system.

(b) LEMMA 27.2 (RIGOROUS NEGATIVE -- constrains what can work).  M7.7 needs
        ||q_k sqrt(a_j)|| >= C q_k^(-1/(3n)) = exp(-(log eps) k/(3n) + O(1)),
    i.e. a bound of the form exp(-c k).
      * BAKER: with Lambda = k log eps + log beta_j - log m, max|coef| = k and
        h'(m) ~ k log eps, Baker gives log|Lambda| >= -C' k log k, hence
        ||.|| >= exp(-C' k log k) -- WEAKER than exp(-c k) by a factor log k.
      * NORM: q_k sqrt(a_j) - m lies in a quartic field; norm bound gives
        exp(-k log eps) -- worse.
    ==> BOTH classical tools are asymptotically INSUFFICIENT.  Any successful
    argument must exploit the PELL LINKAGE itself, not merely that q_k is an
    exponential sequence.  This is the first hard constraint on the shape of a
    working theorem.

(c) THE TAIL IN THE k-VARIABLE (exact integer arithmetic, n=2, delta=1/6):
      a1 a2  #Pell #exc  last k                 last q      log2q*  log2q*/k  R(k*)
       2  3     260    20      42  10,181,446,324,101,389       53.2     1.266  0.1335
       2  5     260    19      32     1,513,744,654,945          40.5     1.264  0.1277
       2  7     260    19      21            93,222,358         26.5     1.261  0.1092
       3  2     260    15      21     1,104,728,155,436         40.0     1.905  0.2650
    TWO structural facts: (i) log2 q*/k = log2 eps exactly (1.2714 for a1=2),
    so q* = eps^{k*} and the k-variable is a faithful linearisation;
    (ii) R_a(k*) sits RIGHT AT the threshold log(eps)/(3n) (0.109-0.265 vs
    0.147-0.168) -- the condition is genuinely borderline, not slack.
    NOTE this also refines q* for (sqrt2,sqrt3) from 1.75e15 (Docs/26) to
    1.018e16.

(d) CANDIDATE LIST IS EXPONENTIALLY SHORT: k* ~ 21-45 replaces q* ~ 10^8-10^16.
    So the live problem is a bound on an INDEX, not a magnitude, and
    log2 q* = O(N^2)  <=>  k* = O(N^2/log eps).

(e) METRIC MODEL CROSS-CHECK.  For n=2, P(||q_k sqrt a2|| < eps^(-k/6)) ~
    2 eps^(-k/6); sum_k ~ 12.2 for a1=2, against 19-20 observed -- GOOD
    agreement.  (The n-dimensional volume model of Docs/26 Lemma 5, which
    counted over ALL q rather than over k, under-predicted by ~160x.  Both
    recorded; the k-model is the accurate one because the candidate set is q_k.)

(f) OPEN CORE (M7.9):  Effectively bound the indices k for which
        R_a(k) = -(1/k) sum_{j>=2} log||beta_j eps^k||  >  (log eps)/(3n)
    at scale k* = O(N^2/log eps).  NOT closed by Baker (Lemma 27.2), not by
    the norm, not by Docs/25 (no algebraic rigidity), not by Docs/26 Lemma 5
    (metric magnitude wrong).  The natural remaining candidates are
    S-unit / recurrence-mod-1 theorems (Corvaja-Zannier type) applied to beta_j
    on the torus under multiplication by eps; the question is whether ANY give
    an EFFECTIVE RATE at the exp(-c k) scale rather than exp(-C k log k).
    NOT VERIFIED, NOT ASSERTED.

L6. FIELD / GALOIS STRUCTURE (DONE).  Wrote `Docs/28_FIELD_STRUCTURE.md`.
Question: before invoking general recurrence-mod-1 machinery, what does the
special shape beta_j = sqrt(a_j/a_1)/2 buy?  ANSWER: exactly the norm bound.

(a) THEOREM 28.1 (verified exactly).  In K = Q(sqrt a_1, sqrt a_2):
      beta_j = sqrt(a_1 a_j)/(2 a_1),   beta_j^2 = a_j/(4 a_1) IN Q
      2 a_1 beta_j eps^k = sqrt(a_2)(p_k sqrt(a_1) + a_1 q_k)   [explicit]
      (beta_j eps^k)^2 = (a_j/(4 a_1)) eps^(2k)  IN Q(sqrt a_1)
      conjugates of beta_j eps^k = { +- beta_j eps^{+-k} }
    => beta_j eps^k is a PELL-TYPE element of K, on an explicitly parametrised
    rank-4 lattice generated by two integer recurrence sequences.
    Verified: cl.1 to 1e-16 (3/8 exactly), cl.3 to 1.9e-15 over k=6..11,
    cl.2 to 9 digits.  (My first printed form of cl.2 was a transcription
    error, caught numerically and corrected.)

(b) THEOREM 28.2 -- THE DECISIVE NEGATIVE.  The field/Galois structure gives
    exactly   eta_k := ||q_k sqrt(a_j)||  >=  q_k^-3 = exp(-3k log eps),
    by THREE independent routes:
      (i)  direct norm in K (conjugates +-beta eps^{+-k});
      (ii) square first ((beta eps^k)^2 in Q(sqrt a_1)), then norm there;
      (iii) the (A - D sqrt a_1) integer form.
    Requirement is exp(-(log eps) k/(3n)).  SHORT BY A FACTOR 9n IN THE
    EXPONENT.  So there is NO hidden rigidity in the conjugates: the element is
    Pell-type but its norm against an integer is NOT small, only its distance
    to Z is.  => GALOIS STRUCTURE EXHAUSTED.
    Empirically eta_k/q_k^-3 grows to ~2e48 by k=32 -- the Liouville floor
    q^-3 is never remotely approached, so it describes nothing.

(c) THEOREM 28.3 (KEY POSITIVE) -- the four-tool barrier, and the ONE tool
    that reaches the required scale:
        algebraic norm (Docs/22,25,28)  -> exp(-3k log eps)   short by 9n
        Baker / lin. forms in logs (27) -> exp(-C k log k)    short by log k
        metric / equidistribution (26)   -> wrong magnitude (~160x off)
        QUANTITATIVE SUBSPACE THEOREM    -> k^{-c} = exp(-c log k)  SUFFICIENT
    Because exp(-c log k) > exp(-(log eps)k/(3n)) for all k >= K_0, the
    exceptional indices form a BOUNDED INITIAL SEGMENT [1, K_0) with
        K_0 = max(K_ST, min{k : 3n c log k < (log eps) k}) = O(n c/log eps).
    Hence  k* < K_0.

(d) THE REDUCTION.  Docs/24 Thm 24.5 demanded log2 q* = O(N^2), i.e.
    k* = O(N^2/log eps).  Thm 28.3 replaces this by k* < O(n c/log eps) + K_ST.
    n is the RANK; c and K_ST depend only on the FIXED system (beta_j, eps).
    NEITHER DEPENDS ON THE RADICAND SIZES, i.e. NOT ON N.
    ==> The required effectiveness has COLLAPSED from O(N^2) to a bound
    independent of the input size.  First time in the project that the
    required constant does not grow with N.  Any explicit Subspace-Theorem
    threshold for this one fixed system suffices -- far weaker than "make the
    ST effective for multiquadratic fields".

(e) OUTSTANDING, NOT ASSERTED: an explicit (c, K_ST) for
        ||beta_j eps^k|| >= k^{-c}  for all k >= K_ST
    on this system.  The ST gives this QUALITATIVELY; effectivity is exactly
    the open input.

(f) NUMERICAL CAVEAT (important, Docs/28 section 28.5.3).  Float evaluation of
    ||q sqrt(a_j)|| loses all precision for q >~ 1e15 (at q=3e14 one already
    sees ||q sqrt3|| = 0.0 and 2.0, both artefacts).  ALL exception
    DECISIONS in Docs/24, 26, 27, 28 used EXACT integer arithmetic (math.isqrt
    on q^2 a + nearest-residue search) and stand; but any DISPLAYED magnitude
    at that scale should be read as "exact value not computed".

(g) k* for (sqrt2,sqrt3): exact arithmetic gives exceptions to k = 42,
    q* = 10,181,446,324,101,389.  No scan beyond k = 44, so k* is a LOWER
    bound only.

L7. M7.11 -- EFFECTIVE POLYNOMIAL IRRATIONALITY (DONE; one question left).
Wrote `Docs/29_M711_UNIFORMITY.md`.

(a) THE FINITENESS IS A PUBLISHED THEOREM (verified, with verbatim quotes).
    Corvaja-Zannier, Acta Math. 193 (2004) 175-191 (arXiv:math/0403522):
      Thm 1 verbatim: "Let alpha>1 real algebraic, 0<lambda<1.  Suppose
        ||alpha^n|| < lambda^n for infinitely many n.  Then some power alpha^n
        is a Pisot number."
      MAIN THEOREM verbatim: "Let Gamma be a finitely generated
        multiplicative group of algebraic numbers, let beta != 0 and let
        delta>0 be fixed.  Then there are only finitely many pairs (q,u) in
        Z>0 x Gamma ... such that ||q beta||_v [is not too small vs]
        h(1/(q beta)) <= ... l(v,d_chi)^(1/delta)."
      => an effective LOWER BOUND for ||q beta|| against heights, with finitely
      many explicitly controlled exceptions.  Correct shape.
    Kulkarni-Mavraki-Nguyen, Trans. AMS (2018) (arXiv:1511.08525), abstract
      verbatim: "We completely solve the problem of the existence of theta in
      (0,1) such that there are infinitely many tuples (n,q_1,...,q_k)
      satisfying ||q_1 alpha_1^n + ... + q_k alpha_k^n||_Z < theta^n ... The
      case k=1 and q_1 is essentially a rational integer was obtained by
      Corvaja and Zannier and settled a long-standing question of Mahler."
      => OUR PROBLEM IS EXACTLY k=1, alpha_1 = eps, q_1 = beta_j (fixed, bounded
      height).  SETTLED by Corvaja-Zannier.
    Effectivity of the surrounding S-unit machinery: Waldschmidt (2015), "the
      result of Corollary 1 IS effective: one can bound from above the
      (numerators and denominators of the) solutions x_1 and x_2".

(b) LEMMA 29.1 (PRODUCT -> SINGLE FACTOR, proved).  If
      prod_{j>=2} ||beta_j eps^k||  <  eps^(-k/(3n)),
    then SOME j>=2 satisfies  ||beta_j eps^k|| < eps^(-k/(3n(n-1)))
                                       =  theta_0^k,
    theta_0 := eps^(-1/(3n(n-1))).
    => M7.9's SIMULTANEOUS product condition implies a ONE-DIMENSIONAL
    condition of the CZ/KMN type.  No simultaneous quantitative ST is ever
    needed.  (This is exactly the reduction proposed.)

(c) THEOREM 29.2 -- THE DISCRIMINANT REPAIR (the real prize).  After the
    reduction the working field is K = Q(sqrt a1, sqrt a2): degree 4, and
        Disc(K) = D(a1) D(a2) D(a1 a2)  =>  log2 Disc(K) = 2(log2 a1 + log2 a2)
        + O(1) = O(N).
    Measured (N = radicand bit-size), ratio log2 Disc / N:
        (2,3)->3.72  (2,13)->2.68  (101,103)->2.19  (1009,1013)->2.00
    -- bounded, non-growing.  CONTRAST Docs/23 section 23.7: the full
    multiquadratic field has log2 Disc ~ 2^(m-1) N, EXPONENTIAL in the rank,
    which is why gamma >= Disc(K)^(-O(1)) was refuted there.
    ==> The reduction removes the rank from the field entirely and thereby
    REMOVES THE Docs/23 OBSTRUCTION.  Any constant polynomial (or exp(poly)) in
    log Disc K, or in [K:Q]=4, is acceptable for the O(N^2) target.

(d) THE COST, RECORDED HONESTLY (Docs/29 section 29.4).  The reduction is NOT
    free: it weakens the threshold.
        product : theta_prod = eps^(-1/(3n))        (stronger, harder to hit)
        reduced : theta_0    = eps^(-1/(3n(n-1)))  (weaker, easier to hit)
    and theta_0 -> 1 as n -> inf:
        n=2: 0.8634 (a1=2) / 0.8458 (a1=3)
        n=3: 0.9522 / 0.9457
        n=4: 0.9758 / 0.9725
        n=6: 0.9903 / 0.9889
    So the uniformity question gets strictly HARDER as the rank grows.
    Trade-off table:
        simultaneous product: field degree 2^m, log2 Disc ~ 2^(m-1)N (FATAL),
                              threshold stronger, NO CZ/KMN theorem applies
        reduced single factor: field degree 4, log2 Disc = O(N) (REPAIRED),
                              threshold weaker and -> 1, CZ/KMN DOES apply
    Which side wins is OPEN.  Docs/29 does NOT claim the trade is favourable.

(e) Docs/28 section 28.4 RETRACTED (Docs/29 section 29.1).  The sentence
    "Effective Subspace-Theorem thresholds for a fixed system of algebraic
    linear forms are classical" conflated (1) effective control of solutions to
    a fixed-system ST problem with (2) an explicit ||beta eps^k|| >= k^{-c}.
    (1) does NOT automatically give (2): one must exhibit the ST variables,
    show an exceptional k yields an integer point in an exceptional subspace
    WITH AN EFFECTIVE HEIGHT THRESHOLD, then eliminate those subspaces or
    convert each to an effective recurrence/linear-form problem.  Not done.

(f) M7.11 REDUCED TO ONE QUESTION:  prove
        theta_0 = eps^(-1/(3n(n-1)))  <  theta_crit(a_1,a_2)
    uniformly over distinct squarefree a_1 < ... < a_n, where theta_crit is the
    critical value furnished by Corvaja-Zannier (2004) / KMN (2018).

(g) STILL OPEN, NOT ASSUMED:  (i) the uniformity inequality (f); (ii) that the
    orbit is non-degenerate in the KMN sense -- note eps = 1+sqrt(a1) IS a
    PISOT NUMBER (conjugate 1-sqrt(a1) has |.|<1) so ||eps^k|| -> 0 and the
    DEGENERATE regime is adjacent; our beta_j eps^k stays Theta(1) empirically
    (0.05-0.47 for (2,3), k<=44) but non-degeneracy must be VERIFIED, not
    inferred from numerics.

L8. KMN NON-DEGENERACY AUDIT (DONE -- M7.11a RESOLVED).  Wrote
`Docs/30_KMN_AUDIT.md`.  Audited against the VERBATIM definitions in
arXiv:1511.08525v1 (Kulkarni-Mavraki-Nguyen, Trans. AMS 2018).

(a) VERBATIM Definition 1.2 (non-degenerate): "The tuple of non-zero
    algebraic numbers (alpha_1,...,alpha_k) is called non-degenerate if
    alpha_i/alpha_j NOT in mu for 1 <= i != j <= k."
    ==> COROLLARY 30.1: at k=1 this is VACUOUS (it quantifies over pairs i != j,
    of which there are none).  M7.11a IS RESOLVED AND REQUIRED NO WORK.
    This is the answer to "audit the exact KMN definition of non-degenerate".

(b) VERBATIM Definition 1.3 (pseudo-Pisot): "Let (beta_1,...,beta_k) be a
    tuple of distinct non-zero algebraic numbers.  Write B := {beta in Qbar*
    \ {beta_1,...,beta_k} : beta = sigma(beta_i) for some sigma in G_Q,
    1 <= i <= k}.  Then the tuple is called pseudo-Pisot if
    sum_i beta_i + sum_{beta in B} beta is in Z and |beta| < 1 for every
    beta in B."
    VERBATIM Theorem 1.4 gives, for all but finitely many elements of an
    infinite M:  (i) (q_1 alpha_1^n,...,q_k alpha_k^n) is pseudo-Pisot;
    (ii) alpha_i is an algebraic integer;  (iii) sigma(q_i alpha_i^n) =
    q_j alpha_j^n  <=>  sigma(alpha_i)/alpha_j in mu;  (iv) if
    sigma(alpha_i)/alpha_j not in mu then |sigma(alpha_i)| < 1.

(c) PROPOSITION 30.2 (CONDITION (i) FAILS EXPLICITLY).  Our tuple is the single
    element beta_j eps^k; by Docs/28 Thm 28.1(4) its conjugates are
    {+-beta_j eps^k, +-beta_j eps^{-k}}, so B contains -beta_j eps^k with
    modulus -> infinity.  Pseudo-Pisotness requires |beta| < 1 for EVERY
    beta in B, hence requires beta_j eps^k < 1.  So (i) FAILS for all large k,
    BY A PLAIN ARCHIMEDEAN SIZE COMPARISON -- no Subspace-Theorem content.
    (The integrality half of Def 1.3 does hold: the sum is the trace, = 0.)

(d) THEOREM 30.3 (Finiteness, uniform in theta).  For distinct squarefree
    a_1 != a_j, eps = 1+sqrt(a_1), beta_j = sqrt(a_j/a_1)/2:
        { k : ||beta_j eps^k|| < theta^k }  is FINITE for EVERY theta in (0,1),
    with NO constant depending on the data entering the argument.
    (Thm 1.4 asserts (i) for all but finitely many of an INFINITE M; Prop 30.2
    says (i) fails for all large k; contradiction.)

(e) COROLLARY 30.4 (A NEW QUANTITATIVE STATEMENT).  eta_k >= theta^k for all
    large k for every theta < 1, i.e.  liminf_k eta_k^(1/k) >= 1, i.e.
        eta_k = e^{-o(k)}   -- SUB-EXPONENTIAL.
    Far stronger than the elementary e^{-3k log eps} = e^{-Theta(k)} of
    Docs/28 Thm 28.2.

(f) PROPOSITION 30.5 (second, independent obstruction).  Condition (iv) needs
    |1-sqrt(a_1)| = sqrt(a_1)-1 < 1, which holds iff a_1 < 4, i.e. iff
    a_1 in {2,3}.  So for a_1 >= 5, eps is NOT Pisot and (iv) ALSO fails.
    Clean dichotomy: (i) kills it for a_1 in {2,3}; (iv) kills it for a_1 >= 5.
    Doubles as a sanity check that 1+sqrt(a_1) is Pisot only for a_1 in {2,3}.

(g) CITATION CORRECTION (Docs/30 section 30.5).  KMN verbatim: "The case k=1
    and q_1 is essentially a RATIONAL INTEGER was obtained by Corvaja and
    Zannier."  Our q_1 = beta_j is NOT rational, so CZ Theorem 1.1 does NOT
    apply -- KMN is the operative theorem, and Docs/29 section 29.2 is
    corrected.  ESSENTIAL: if q_1 = 1 and a_1 in {2,3}, eps is Pisot and
    ||eps^k|| ~ (sqrt(a_1)-1)^k, so ||eps^k|| < theta^k for EVERY k once
    theta > sqrt(a_1)-1 -- INFINITELY many exceptions.  It is precisely
    beta_j NOT in Q that removes us from the degenerate regime.

(h) THE REMAINING GAP, PRECISELY BOUNDED ON BOTH SIDES (Docs/30 section 30.6):
      elementary (Liouville + norm)  :  eta_k >= e^{-3k log eps} = e^{-Theta(k)}
      KMN (published, NON-effective):  eta_k >= e^{-o(k)}
      Docs/28 Thm 28.3 NEEDS         :  eta_k >= e^{-O(log k)} = k^{-c}
    e^{-o(k)} does NOT imply k^{-c} (e^{-sqrt k} is e^{-o(k)} yet << k^{-c}).
    KMN's proof rests on the Subspace Theorem, so it inherits the
    ineffectivity; no explicit K_ST or c is produced.

(i) M7.11b RESTATED (M7.11b'): improve  eta_k >= e^{-o(k)}  to
        eta_k >= e^{-C log k}   for an explicit C = C(beta_j, eps).
    C MAY depend on a_1, a_j with NO need to control that dependence in N,
    because Docs/29 section 29.5 absorbs C into K_0 = O(nc/log eps), already
    N-independent.  THIS IS THE ENTIRE REMAINING GAP: one exponential step
    from a PUBLISHED NON-EFFECTIVE bound.

L9. KMN PROOF AUDIT -- DOES IT CONTAIN THE RATE? (DONE -- ANSWER: NO).
Wrote `Docs/31_KMN_PROOF_AUDIT.md`.  All quotations verbatim from
arXiv:1511.08525v1.

(a) WHERE THE SUBSPACE THEOREM IS USED.  Their Theorem 2.1 (ST) is invoked in
    exactly TWO places: Proposition 2.2 (finiteness of exact zeros, by
    minimality in k) and Proposition 3.4 (the key technical result).
    Everything else is Northcott / Gaussian elimination / elementary.

(b) PROP 2.3's k=1 CASE IS EFFECTIVE.  Verbatim: "When k=1, the given
    inequality together with the Northcott property imply that the set S must
    finite [be finite].  Hence we may assume k>=2.  We now apply the Subspace
    Theorem as in the proof of [CZ04, Lemma 1]."
    ==> For k=1 the ONLY ineffective step in the entire argument is inside
    PROPOSITION 3.4.

(c) THE SINGLE INEFFECTIVE STEP.  In Prop 3.4 the ST is applied to
    x = (p, y_b) and yields a nontrivial relation a_0 p + sum_b a_b y_b = 0
    with a_0, a_b in L INDEPENDENT of n, from which (15) p = sum eta_{i,j}
    alpha_{i,j}^n with h(eta) = o(n) is derived.  THE ENTIRE LOSS IS THE
    PHRASE "for all but finitely many elements of M".  Everything downstream is
    elementary.  KMN proves M_theta finite but gives NO CONTROL on max M_theta,
    which is exactly what Docs/28 Thm 28.3 / Docs/30 M7.11b' need.

(d) THE QUANTITATIVE ST CANNOT REPAIR IT.  Schmidt (1989), Evertse (1995),
    Evertse-Schlickewei (2002), Evertse-Ferretti (2013) bound the NUMBER t of
    exceptional subspaces (Schlickewei: t <= c(n) delta^-3 (log(1/delta)-1)^2).
    THE HEIGHTS AND IDENTITIES OF THOSE SUBSPACES REMAIN UNBOUNDED.
    COUNTING IS NOT LOCATING.  This is the decisive negative on the most
    obvious repair.

(e) STRUCTURAL RIGIDITY THAT DOES HELP (Lemma 31.2).  Because q_1 = beta_j is
    FIXED, the coordinates y_{(j1,j2,l)} = b_{1,j1,l} * eps_{j2}^k have b_l
    FIXED RATIONAL (the K-coordinates of beta_j, via their Lemma 3.1).  So the
    ST relation collapses to
        a_0 p_k + sum_{j=1..4} zeta_j eps_j^k = 0     with zeta_j FIXED.
    Substituting p_k = beta_j eps^k + O(theta^k) and comparing successive terms:
    if a_0 beta_j + zeta_1 != 0 the LHS grows like eps^k; if it vanishes, zeta_2
    must vanish, etc.  ==> EACH FIXED RELATION HOLDS FOR AT MOST AN EXPLICITLY
    BOUNDED SET OF k.  So the problem reduces to IDENTIFYING finitely many
    FIXED relations -- each already provably harmless with an explicit
    threshold.

(f) WHY IT STILL DOES NOT CLOSE (section 31.3).  From (15), p = sum eta_j
    eps_j^k with h(eta_j) = o(k), and p_k ~ beta_j eps^k, Tr(eps^k) ~ eps^k,
    one gets  eta_j/beta_j = 1 + O(e^{-ck}).  Liouville then gives only
        |eta_j - beta_j| >= c/den(eta_j)^2 = e^{-h(eta_j)} = e^{-o(k)},
    which is FULLY COMPATIBLE with e^{-ck}.  ==> The o(k) in the ST output is
    EXACTLY CALIBRATED TO THE LIOUVILLE DENOMINATOR SCALE; the fixed
    biquadratic/Pell structure supplies NO purchase on that gap.  An upgrade to
    h(eta) = O(log k) is needed and is NOT obtained.
    Note the o(k) is INHERITED from the hypothesis h(q_i) < f(n) with f
    sublinear; with q_1 = beta_j fixed the input is far better than o(k), so a
    route via a STRENGTHENED Subspace-Theorem statement exists in principle, but
    it is not a calculation.

(g) M7.12 PRECISELY:  identify the finitely many relations of Prop 3.4
    effectively for k=1, alpha_1 = eps = 1+sqrt(a_1), q_1 = beta_j fixed.
    Their COUNT is effectively bounded (quantitative ST); their IDENTITY is not.

L10. RELATION ENUMERATION (DONE -- CANDIDATES ENUMERATED *AND REFUTED*).
Wrote `Docs/32_RELATION_ENUMERATION.md`.

(a) CORRECTION to Docs/31 section 31.2.  The claim there that the zeta_j are
    "FIXED ... drawn from a finite, explicitly enumerable set" is WITHDRAWN.
        zeta_j = sum_{b in I*, j2(b)=j2} a_b * b_{j1(b),l(b)},
    and the a_b in L are the SUBSPACE THEOREM'S UNKNOWN OUTPUT.  So zeta_j is
    fixed only AFTER the unknown relation is chosen; Lemma 31.2 is valid PER
    RELATION but gives no effectivity alone.  Only the b_{1,j1,l} -- the
    K-coordinates of the FIXED q_1 = beta_j -- are genuinely fixed.

(b) THEOREM 32.1 (NEW, EFFECTIVE DIMENSION COLLAPSE).  Because q_1 is fixed,
    the components of KMN's vectors are b_{1,j1,l} * eps_{j2}^k with the b's
    FIXED RATIONAL, so
        y(k) = eps^k E_1 + eps'^k E_2     with E_1, E_2 FIXED in L.
    Hence  r := dim_L V <= 2,  |I*| = r <= 2, and the Subspace Theorem is
    applied in L^(r+1) with r+1 <= 3.  r is EXACTLY COMPUTABLE from the fixed
    data.  KMN do not state this; it is specific to fixed q_1.

(c) LEMMA 32.2 (REDUCTION TO A GEOMETRIC ROUNDING ERROR).  An admissible
    relation  a_0 p_k + Psi_1 eps^k + Psi_2 eps'^k = 0  for infinitely many k
    forces, on dividing by eps^k and letting k -> infinity (e_k/eps^k -> 0
    since theta<1<eps; (eps'/eps)^k -> 0):
        Psi_1 = -a_0 beta_j ,
    and then (a_0 = 0 forces the trivial relation)  e_k = zeta eps'^k with
    zeta = -Psi_2/a_0 in L^*.  So an admissible relation exists ONLY IF the
    rounding-error sequence  e_k = round(beta_j eps^k) - beta_j eps^k  is
    eventually an EXACT geometric progression.

(d) THEOREM 32.3 (COMPLETE CONTRADICTION -- the key new result).  Taking
    traces of e_k = zeta eps'^k:  Tr(e_k) = 4 p_k (since Tr(beta_j eps^k)=0
    and Tr(p_k)=4p_k) and Tr(eps'^k) = 2(eps^k + eps'^k), so
        4 p_k = 2 zeta (eps^k + eps'^k),   p_k = (zeta/2)(eps^k + eps'^k).
    With p_k = beta_j eps^k + zeta eps'^k this gives
        (beta_j - zeta/2) eps^k = -(zeta/2) eps'^k.
    Dividing by eps^k, using |eps'/eps| = (sqrt a1 - 1)/(sqrt a1 + 1) < 1 for
    EVERY a1 >= 2, gives zeta = 2 beta_j; then beta_j (eps'/eps)^k = 0 for
    infinite k, impossible since beta_j > 0.
    ==> NO admissible relation exists, FOR EVERY squarefree a1, with NO case
    split on a1, NO use of conditions (i) or (iv), and no asymptotic or
    quantitative input.  PURELY ALGEBRAIC.

(e) CONSEQUENCES.  The candidate relations form a PENCIL in L^3 cut out by the
    SINGLE EXPLICIT linear equation Psi_1 = -a_0 beta_j, and Thm 32.3 refutes
    EVERY point of it uniformly.  So the "identify the finitely many
    relations" half of M7.12 is CLOSED.  Also the Docs/30 Prop 30.5 case
    split on a1 >= 5 is UNNECESSARY -- Lemma 32.2(3) covers all a1 at once,
    including a1 = 2,3 where |eps'| < 1.

(f) WHAT REMAINS, AND IT IS ONE SET (Docs/32 section 32.4).  Finiteness is
    unchanged.  The entire residual loss is
        M_theta \ M_0'   --  finite, but its size and largest element are
    NOT bounded by anything in KMN.  Every k in M_theta cap M_0' would lie in
    a hyperplane, hence give an admissible relation, hence be IMPOSSIBLE by
    Thm 32.3.  So  M_theta is contained in {k : the ST step of Prop 3.4 fails
    at k} -- and THAT set is exactly the "all but finitely many" the Subspace
    Theorem declines to identify, i.e. the classical non-effectivity (heights
    of the exceptional hyperplanes).  The identities are SETTLED; only the
    heights are not.  Docs/31 section 31.4's remark (quantitative ST bounds
    the NUMBER t but not the heights) is now literally the WHOLE problem.

(g) M7.12 (FINAL FORM):  bound the Subspace-Theorem exceptional tail
        E_theta := { k : the ST step of Prop 3.4 fails at k } .
    NO PROGRESS CLAIMED ON THIS.  But it is now the ONLY thing left, it is a
    single classical question, and the surrounding analysis is complete.

L11. EFFECTIVE COUNT BOUND (DONE -- NEW RESULT).  Wrote `Docs/33_COUNT_BOUND.md`.
First EXPLICIT CARDINALITY bound on the EHS/KMN exceptional set; Docs/24-30 had
only finiteness.

(a) THEOREM 33.1 (ST EXPONENT AND THRESHOLD BECOME EXPLICIT).  KMN's step (13)
    is  prod |L(x)|/|x|  <=  theta^k * B^r * H(x)^(-r-1)  with
    B := max H(b_{i,j,l}).  In general B = e^{o(k)} (hypothesis h(q_i)<f(k),
    f sublinear), so "there is eps>0 ... when k is sufficiently large" has NO
    explicit k_0.  BUT q_1 = beta_j is FIXED, so the b's are FIXED RATIONALS
    and B is a CONSTANT.  Hence
        eps := -log(theta)/2 > 0 ,  n_0 = ceil((r log B + log(1/eps))/(-log theta))
    is EXPLICIT, and n_0 = O(log B).  The fixedness of q_1 removes the o(k)
    from the ST exponent -- the same structural fact that gave the dimension
    collapse in Docs/32.

(b) THEOREM 33.2 (ONE EXCEPTION PER HYPERPLANE, COMPUTABLE).  For a FIXED
    relation a_0 p_k + Psi_1 eps^k + Psi_2 eps'^k = 0 with zeta = -Psi_2/a_0:
        (beta_j - zeta/2) eps^k = -(zeta/2) eps'^k
    so  k = log(|zeta| / (2|beta_j - zeta/2|)) / log(eps/|eps'|),
    which DETERMINES k UNIQUELY from zeta.  Hence each exceptional hyperplane
    contains at most ONE element of M_theta, explicitly computable.  (Docs/32
    Thm 32.3 is the k -> infinity special case.)

(c) THEOREM 33.3.  |M_theta|  <=  n_0 + t,  n_0 explicit, t any valid bound
    for the number of exceptional hyperplanes.

(d) THEOREM 33.4 (FIRST EXPLICIT COUNT BOUND).  By Docs/32 Thm 32.1, r <= 2, so
    the ST is applied in L^3 and Schlickewei's bound is
        t <= c(3) * delta^-3 * (log(1/delta)-1)^2 ,  delta = -log(theta)/2,
    with c(3) ABSOLUTE (the variable count is the fixed 3, not the rank).
    Measured t/c(3) for a_1=2: rank 2 -> 6.6e3, 3 -> 5.0e5, 4 -> 6.3e6,
    6 -> 1.6e8, 10 -> 6.8e9 -- POLYNOMIAL in the rank (delta^-3 ~ n^6), with NO
    dependence on a_1, a_j or N.  Since log B = O(N) (B = heights of the
    K-coordinates of the fixed beta_j), n_0 = O(N), so
        |M_theta|  =  O(N) + O(delta^-3 (log 1/delta)^2).

(e) CORRECTION to Docs/31 section 31.4.  It said the quantitative ST is useless
    here ("counting is not locating").  PARTLY CORRECTED: counting IS enough
    for the CARDINALITY (Thm 33.3).  What it cannot do is LOCATE the
    hyperplanes, hence bound k* = max M_theta.

(f) WHAT IS STILL MISSING, PLAINLY.  SRS needs k*, not |M_theta|, and (33.4)
    does NOT give it: k = k(zeta) with zeta ranging over ALL of L, so k(zeta)
    is unbounded (zeta -> 2 beta_j makes |beta_j - zeta/2| -> 0 and k(zeta) ->
    infinity).  The exceptions with k < n_0 ARE fully enumerable (n_0 = O(N),
    explicit); those with k >= n_0 are at most t, each equal to k(zeta_i), but
    the zeta_i are not identified.

(g) M7.12 (FINAL, SECOND ITERATION):  bound |zeta| away from its degeneracy
    near 2 beta_j -- equivalently bound k(zeta) -- uniformly over the at most t
    exceptional values produced by Theorem 2.1 in L^3.  This is now a FINITE,
    BOUNDED, ALGEBRAIC IDENTIFICATION problem, not an analytic effectivity
    problem.  EVERYTHING ELSE IS DONE.

L12. THE HEIGHT LEMMA (DONE -- NO-GO PROVEN; CASE A REFUTED).  Wrote
`Docs/34_HEIGHT_LEMMA.md`.

(a) TWO CORRECTIONS TO THE RECORD.
    * Docs/27 section 27.2 displayed  q_k sqrt(a_j) = beta_j(eps^k - eps^-k)/2
      -- with a SPURIOUS /2.  Correct:  q_k sqrt(a_j) = beta_j(eps^k - eps^-k).
      Verified at m=9: beta_j(eps^9 - eps^-9) = 1706.0700454553433 vs
      Q_9 sqrt(3) = 1706.0700454553441.  NO DOWNSTREAM CONSEQUENCE: the
      substance of Docs/27 Lemma 27.1, ||beta_j eps^k|| = ||q_k sqrt(a_j)|| +
      O(eps^-k), was verified numerically and is correct, as is the whole
      M7.7 -> M7.12 chain (Docs/28 - Docs/33).
    * Docs/34 section 34.2 -- a PERVASIVE INDEX SHIFT.  The standard CF
      recursion gives eps^(j+1) = p_j + Q_j sqrt(a_1): list index j corresponds
      to exponent m = j+1.  This offset produced REPEATED false conclusions in
      my own verification scripts (including a spurious "BAD" identity check
      this session).  Convention fixed: every Pell-index identity is written
      with eps^(j+1).

(b) THE zeta REFORMULATION, CORRECTED.  An exceptional k forces
    e_k = zeta (eps')^k with zeta = -Psi_2/a_0, and taking traces gives
        zeta_k = 2 p_k / ( eps^k + (eps')^k ) ,   DETERMINED BY k .
    Hence  zeta_k - 2 beta_j = -2 beta_j rho^k/(1+rho^k)
                              + 2 e_k/(eps^k (1+rho^k)) .
    CORRECTION to the proposed rate: it is theta/eps, NOT |rho|.  Verified at
    k=41 (a1=2,a2=3, exact integers + 60-digit magnitudes): observed
    zeta-2beta = -1.637e-18, while 2|e_k|/eps^k = 1.66e-18 and
    2 beta_j |rho|^k = 1.3e-32 -- the ROUNDING-ERROR term dominates.
    Also confirmed: h(zeta_k) = O(k);  |zeta_k - 2 beta_j| <= 2 (theta/eps)^k.

(c) THEOREM 34.1 (NO-GO FOR THE HEIGHT/Liouville ROUTE).  Liouville in the
    degree-4 field L gives only a POWER law,
        |zeta_k - 2 beta_j|  >=  c H(zeta_k - 2 beta_j)^(1-4)  ~  c k^-3 ,
    while the ACTUAL value is EXPONENTIALLY SMALL (<= 2(theta/eps)^k).  A
    power-law lower bound is WEAKER than an exponentially small value, hence
    ALWAYS SATISFIED -- no contradiction, no upper bound on k.
    CONSEQUENCES:
      * The proposed 'minimal sufficient theorem'
            |zeta_i - 2 beta_j| >= exp(-poly(N))
        is EQUIVALENT to bounding k.  It is NOT a consequence of heights or
        Liouville -- it IS the statement to be proved.  The route is CIRCULAR.
      * The obstruction is structural: zeta_k - 2 beta_j is exponentially
        small in L for EVERY k -- exactly what units produce.  Nothing
        arithmetically anomalous.  So no EXISTENCE obstruction can exist; the
        last lemma must be about WHICH zeta_k the ST SELECTS.
      * Docs/31 section 31.3's warning is confirmed from a second direction.

(d) CASE A IS REFUTED (section 34.5).  In KMN Prop 3.4 the relation
    a_0 p + sum a_b y_b = 0 is produced by Theorem 2.1 (the ST); the a_b are
    the ST's output and NO quantitative ST variant bounds the COEFFICIENT
    HEIGHTS (only the NUMBER t of hyperplanes).  KMN then discard the a_b and
    keep only eta with h = o(k).  And since zeta_k is DETERMINED BY k,
        bounding h(zeta)  <=>  bounding k  <=>  bounding q*  <=>  SRS.
    So "bound the exceptional-hyperplane height" is the SAME problem, not a
    weaker accessible version.  This is the decisive reason Case A cannot help.

(e) THE DOWNSTREAM CHAIN IS VERIFIED MECHANICAL (section 34.6).  With
    q_k ~ eps^k:  k* <= P(N)  =>  log2 q* = O(P(N))  =>  log2(1/gamma) =
    (n+1) log2 q* + O(n) = O(N P(N))  =>  E <= O(N^2) + O(N P(N)) = poly(N).
    So NO further exponential-rank obstacle waits behind k*, PROVIDED the
    constant C = C(a_1,a_j) in Docs/34 (34.5) need not be controlled in N
    (Docs/29 section 29.5).

(f) FINAL M7.12 (section 34.7):
        Prove   ||beta_j eps^k||  >=  k^{-C}   for all k >= K(a_1,a_j),
    for some C = C(a_1,a_j).  Equivalently, bound the k's at which the ST's
    exceptional hyperplane of parameter zeta_k occurs.  This is CASE B.  It
    CANNOT be reached by heights, Liouville, or by extracting a height bound
    from the ST (sections 34.4, 34.5).  Whether an existing low-dimensional
    effective-ST / S-unit theorem reaches it is NOT DETERMINED HERE.

L12'. M7.12 SHARPENED TO A QUANTIFIED GAP (DONE).  Wrote
`Docs/35_M712_SHARP_GAP.md`.  This replaces L12(f)'s qualitative target.

(a) THEOREM 35.1 (NORM / PISOT LOWER BOUND).  For every m in Z and k >= 1,
        |m - beta_j eps^k|  >=  |rho|^{3k}/(2 beta_j^3)
                             =  exp(-3k log eps + O(1)) .
    Proof: N_L(m - beta_j eps^k) is a nonzero integer and the conjugates have
    moduli |m -/+ beta_j eps^{+/-k}|, so 1 <= |y|(2 beta_j eps^k)(beta_j eps^k)
    (2 beta_j eps^k).  EXACTLY VERIFIED on k=2..41 (integer Pell pairs,
    60-digit magnitudes, incl. q = 1.7469e15): ZERO violations, minimum slack
    39.18.

(b) THEOREM 35.2 -- THE AUTOMATIC-PISOT MECHANISM IS ABSENT.
    beta_j NOT in Q(sqrt(a_1)) for distinct squarefree a_1 != a_j (equivalently
    a_1 a_j is not a rational square).  Checked directly for (2,3):
    beta_j = sqrt6/4 = r + s sqrt2 forces rs = 0 and r^2 = s^2 = 3/8.
    Hence there is NO automatic ||lambda theta^k|| ~ |lambda'|^k, which is
    exactly the mechanism that would make M7.12 FALSE.  x_k is trace-zero with
    a conjugate of EQUAL modulus.  Clean algebraic reason Docs/30's
    pseudo-Pisot hypothesis fails.

(c) THEOREM 35.3 -- x_k = beta_j eps^k satisfies x_{k+2} = 2 x_{k+1} + x_k and
    Tr_L(x_k) = 0.  So M7.12 is precisely: "how small can the fractional part of
    a trace-zero order-2 linear recurrence get?"

(d) THEOREM 35.4 (THE SHARP GAP).  Needed to exclude: exp(-k log eps/(3n(n-1))).
    Proven: exp(-3k log eps).  The proven bound is weaker by EXACTLY a factor
    9n(n-1) in the exponent (18x for n=2: proven exp(-2.6442k) vs needed
    exp(-0.1469k)).  This is Docs/28's "short by 9n" made exact.

(e) COROLLARY 35.5 -- RIDOUT / SUBSPACE-THEOREM POWER LAYS ARE STRICTLY WORSE
    THAN THE ELEMENTARY NORM BOUND.  x_k is a 3-term form
    q_j sqrt(a_j) + p_j sqrt(a_1 a_j)/a_1 - m with coefficients ~ eps^k, so
    Ridout would give exp(-(1+eps) k log eps); but 1+eps < 3, so Ridout is
    ~3x WEAKER than Theorem 35.1.  The norm bound wins because it exploits the
    PELL COUPLING of p_k, q_k, which a generic 3-term-form treatment discards.
    => THE STANDARD DIOPHANTINE TOOLKIT IS EXHAUSTED, and the one remaining
    algebraic structure is already fully used by Theorem 35.1.

(f) FINAL M7.12 (QUANTIFIED): for distinct squarefree a_1 != a_j the set
        { k >= 1 : || (sqrt(a_j)/2 sqrt(a_1)) (1+sqrt(a_1))^k ||
                     <= (1+sqrt(a_1))^{-k/(3n(n-1))} }
    is finite, with largest element k* satisfying k* log eps <= poly(N).
    BEST AVAILABLE: exp(-3k log eps) -- short by a factor 9n(n-1).
    BLOCKED ON: an argument improving the exponent by that factor.  Heights,
    Liouville, and standard Subspace-Theorem power laws are each PROVED
    insufficient (Docs/34 Thm 34.1, Docs/35 section 35.6).
    NOT claiming a solution; the record states precisely what is missing.

L12''. THE ENGINE + A STRUCTURAL NO-GO (DONE).  Wrote
`Docs/36_STRUCTURAL_NOGO.md`, `m712_engine.py`, `m712_report.py`,
`tests/test_m712_engine.py`.

(a) TWO MORE CORRECTIONS TO THE RECORD.
    * Docs/35 Thm 35.1 printed `|rho|^{3k}/(2 beta^3)`.  WRONG: `|rho| =
      |eps'/eps| = 0.1716` while `eps^(-1) = |eps'| = 0.4142`, so they differ.
      The norm gives `eps^(-3k)/(2 beta^3 (2 a_1)^4)`, which is STRONGER -- the
      error was optimistic.  Exact re-verification, k=1..60: slack >= 1
      everywhere, minimum 3.0917 at k=1.
    * Docs/35's proof asserted `N_L(y_k) in Z`.  WRONG: `y_k` is not an
      algebraic integer; `N(y_1) = 228/256`.  Only `N(2 a_1 y_k)` is integral,
      so `|N(y_k)| >= (2 a_1)^(-4)`.  This changes the LEADING CONSTANT of
      Thm 35.1 only, never the exponent 3, so every conclusion stands.

(b) BUILT THE ENGINE.  `m712_engine.py` has the four layers and is 100% exact
    (integer + Fraction only); `m712_report.py` holds all logarithm/float code
    and is MEASUREMENT-only.  Both registered in `tests/test_no_float.py`
    (engine = DECISION_PATH, report = MEASUREMENT), so the split is enforced.
    `eps_power` builds eps^k by the exact recurrence (P,Q) -> (P+a_1 Q, P+Q),
    so k IS the exponent BY CONSTRUCTION -- retiring the Docs/34 section 34.2
    index hazard structurally rather than by vigilance.  `nearest_integer`
    brackets with `signexact.certified_bracket` and CONFIRMS with a second
    bracket, raising `UndecidedError` rather than guessing.

(c) LAYER 3 IS DEGENERATE -- now a CHECKED FACT, not a dead layer silently
    kept.  `Q(u,v) = v^2-2uv-u^2` is ANTI-invariant
    (`Q(u_{k+1},u_{k+2}) = -Q(u_k,u_{k+1})`), factorises EXACTLY as
    `(v-eps u)(v-eps' u)`, and `x_{k+1}=eps x_k` kills the first factor, so
    `I_k == 0` identically.  Verified k=1..39; factorisation 200/200 randomised
    exact `Elem` pairs.  SIGN FIX found while testing: `c_k = y_{k+1} - eps y_k`,
    NOT `+` (verified k=10: 0.296351 - 2.414214*0.183089 = -0.145663).

(d) THEOREM 36.1 (NEW -- PAIR / LATTICE BOUND, exponent 1).
        d_{k+1} + eps d_k  >=  eps^(-k) / (2 sqrt(2) beta_j).
    `c_k = m_{k+1} - eps m_k in Z[sqrt(a_1)]` IS an algebraic integer, so
    `N(c_k) = (m_{k+1}-m_k)^2 - a_1 m_k^2 in Z \ {0}`; its surviving conjugate
    is `c_k - 2 sqrt(a_1) m_k`, of modulus `<= 2 sqrt(a_1) m_k + |c_k|`, giving
    `|c_k| >= 1/(2 sqrt(a_1) m_k)`, and `|c_k| = |y_{k+1} - eps y_k|
    <= d_{k+1} + eps d_k`.  STRONGER than Thm 35.1 (exponent 1 vs 3) because it
    lives in the degree-2 subfield rather than in L.  Verified: N(c_k) != 0 and
    slack >= 1 for k=1..59, minimum 4.4853 at k=1.

(e) THEOREM 36.2 (THE STRUCTURAL NO-GO).  NO argument of the form "form an
    algebraic number out of the m_k and bound it by the product of its other
    conjugates" can close M7.12, FOR ANY n >= 2.
      - SINGLE INDEX: the only non-integer available is x_k = beta_j eps^k,
        which by Thm 35.2 has degree 4, and all three surviving conjugates are
        FORCED to be ~ m_k ~ eps^k.  Exponent 3.  No better single-index
        object exists, because degree 2 is unavailable at one index.
      - TWO CONSECUTIVE INDICES: degree 2 IS available, but only as c_k, and the
        conclusion bounds the SUM d_{k+1} + eps d_k, not d_k.  Exponent 1.
    M7.12 needs exponent 1/(3n(n-1)) < 1 for every n >= 2 (1/6, 1/18, 1/30, ...).
    A SMALLER exponent is a STRONGER bound, so the requirement lies strictly
    beyond everything this toolkit produces.  => the obstruction is STRUCTURAL,
    not the `9n(n-1)` tuning deficit of Docs/35 Thm 35.4.
    COROLLARY 36.3: with Thm 34.1 (heights/Liouville) and Docs/35 section 35.6
    (power laws), this kills the LAST family.  M7.12 needs a genuinely different
    KIND of input, not a better inequality.

(f) MEASURED NEGATIVE RESULT on the diagnostic question.  The engine was built
    to find the quantity that grows when d_k falls.  It does not exist:
      * conjugate amplification A_k = |N(2 a_1 y_k)|/(|y_k| eps^{3k})
        converges to 2 beta_j^3 (2 a_1)^4 with RATIO 1.00000 -- so the three
        crude conjugate estimates in Thm 35.1 are asymptotically TIGHT and
        there is NO slack to convert into a better exponent.  This is the
        direct measurement refuting the conjugate-amplification strategy.
      * no layer quantity correlates with -log d_k (|corr| <= 0.21).
      * empirically d_k never decays: min d_k = 0.0036793 at k=43, and no
        exception below theta_0^k for k > 16 (the k <= 16 "exceptions" are the
        artifact theta_0^k > 1/2 at small k).

(g) CORRECTION TO MY OWN SPECULATION.  I suggested Thm 36.1 might isolate
    exceptions.  IT DOES NOT, and the engine says so: adjacent exception pairs
    occur ((1,2),(2,3),...,(6,7),(9,10), minimum gap 1).  Reason: Thm 36.1's
    eps^(-k) sits BELOW the threshold eps^(-k/(3n(n-1))) for n >= 2, so it
    constrains a regime deeper than the one M7.12 asks about.

TESTS: `tests/test_m712_engine.py`, 22 exact tests.  Suite 374 passed,
13966 subtests.  All decisions certified integer/Fraction.

L12'''. TRACK A + TRACK B (DONE -- BOTH NEGATIVE; EXPERIMENTATION NOW FROZEN).
Wrote `Docs/37_TRACK_A_B.md`.

(a) TRACK A -- LITERATURE, NEGATIVE.  Searched for an EFFECTIVE theorem on
    `||beta eps^k||` with `beta` outside Q(eps); a bare finiteness theorem is
    worthless now.  Closest family is Loxton--Van der Poorten: Noubissie,
    arXiv:2504.09519 "Quantitative growth of linear recurrences" (Apr 2025),
    EFFECTIVE bound on the NUMBER of `n` with `|u_n| < (max|alpha_i|)^{n(1-eps)}`,
    completing a programme of Schmidt and Van der Poorten whose qualitative
    forms were proved by Fuchs--Heintze (2021) and Karimov et al. (2023).
    Our threshold `eps^{-k/(3n(n-1))}` IS of that shape.  IT DOES NOT HELP:
    (i) it bounds the COUNT, not the largest index -- and a count bound is
    exactly what Docs/33 Thm 33.3-33.4 already gives (`|M_theta| <= n_0 + t`);
    (ii) it is for INTEGER LRS, and `x_k = beta_j eps^k` is a genuine TWO-TERM
    form (`Docs/37` section 37.2), not an integer recurrence.  Checked and set
    aside: arXiv:2501.17365 (multiplicative dependence in LRS), arXiv:1505.07147
    (Skolem problem), "large zeros of LRS" -- the last two concern `|u_n| = 0`,
    the opposite extreme.  Recorded as a SEARCH RESULT, NOT as a claim that no
    such theorem exists.

(b) TRACK B -- THREE NEW RESULTS, all certified.
    * FACT 37.2:  `||Q_k sqrt(a_1)|| = |eps'|^k = eps^(-k)`  EXACTLY, certified
      k=1..79 via `s1_exact`.  The `a_1` coordinate is NOT generic: it exactly
      consumes the `q^{-1}` in the EHS budget `prod_i ||q sqrt(a_i)|| < q^{-1-delta}`,
      leaving only `prod_{i=2..n} s_i <~ q^{-delta}`.
    * THEOREM 37.1 (PAIR BOUND -- exponent 1; the FIRST use of the SHARED Pell
      pair, which Docs/29's factorwise route discards).  Multiplying the
      two-term forms for coordinates j and j' by sqrt(a_j') and sqrt(a_j) makes
      the P_k and Q_k terms CANCEL:
        m_j' sqrt(a_j) - m_j sqrt(a_j')  =  e_j sqrt(a_j') - e_j' sqrt(a_j) ,
      and the NONZERO integer `I_jj' = m_j'^2 a_j - m_j^2 a_j'` (nonzero because
      `a_j/a_j'` is not a rational square) gives
        |e_j| sqrt(a_j') + |e_j'| sqrt(a_j) >= 1/(m_j' sqrt(a_j) + m_j sqrt(a_j'))
                                        ~  eps^{-k}/(2 beta_{j'} sqrt(a_j)) .
      Certified EXACTLY as a sign test on a SIX-SURD linear form
      (`pair_bound_holds`): 474/474 over 6 radicand pairs x k=1..79;
      `I_jj' != 0` in 354/354.
    * THEOREM 37.3 (PRODUCT BOUND):  `prod_{i=2..n} s_i >= C_n q^{-(n+1)}`,
      `C_n` explicit (`C_3 = 3.27e-16`, `C_4 = 5.95e-30`).  Exponent is
      `3 + (n-2)*1 = n+1`.

(c) TRACK B VERDICT -- NEGATIVE; IT DOES NOT CLOSE.  Multi-coordinate structure
    moves `n-2` coordinates from exponent 3 to exponent 1 and leaves ONE at 3.
    The requirement is an exponent `delta < 1`, so `Docs/36` Thm 36.2's wall
    SURVIVES THE PRODUCT ATTACK UNCHANGED.  Track B does not rescue M7.12.

(d) BUT IT RELOCATES THE QUESTION DECISIVELY.  Theorem 37.3 is non-vacuous
    precisely when `delta > n+1` (n=3 needs delta>4, n=4 needs delta>5, n=5
    needs delta>6).  So the ONLY remaining lever found is:
        DOES EHS ADMIT delta > n+1 IN ITS HYPOTHESIS?
    That is checkable against the EHS statement and is a question about EHS,
    not about us.  NOT CHECKED HERE; no claim either way.

(e) TWO CORRECTIONS AND ONE ENGINE BUG.
    * `Docs/35` section 35.1's identity must read `(eps')^k`, NOT `eps^{-k}`:
      `eps' = -eps^{-1}` so `eps'^k = (-1)^k eps^{-k}` and they differ at ODD
      `k`.  Verified: matches at k=2,10, FAILS at k=1,3,5.  The `Docs/27`
      conclusion (`||beta_j eps^k|| = ||q_k sqrt(a_j)|| + O(eps^{-k})`) is
      UNAFFECTED -- the two differ by exactly `beta_j eps^{-k}` (measured
      1.354e-8 at k=20, matching `beta_j eps^{-20}`).
    * `Docs/37` section 37.2 also states the object correctly:
      `beta_j eps^k = Q_k sqrt(a_j)/2 + P_k sqrt(a_1 a_j)/(2 a_1)` (verified
      exactly at k=5,10,20 to 20 digits), so `||beta_j eps^k||` is a genuine
      TWO-TERM form, not `||Q_k sqrt(a_j)||`.
    * ENGINE BUG, found by the exact certificate: the helper computing
      `sign(L_j')` closed over the radicand `a_j` instead of taking `a_j'` as a
      parameter, deciding j's sign in the WRONG FIELD.  This produced 68
      SPURIOUS failures out of 295.  Two earlier versions of the certificate
      were also wrong (the first discarded the `|e_j'|` term, a weakening that
      fails where the truth holds).  Fixed; recorded because a silently wrong
      radicand in a sign test is exactly the defect class a green suite hides.

TESTS: `tests/test_m712_engine.py` now 27 exact tests, including
`test_pair_bound_certified` (474-case exact certificate),
`test_pair_integer_is_nonzero`, `test_s1_coordinate_is_exact`,
`test_s1_sign_alternates_with_parity`.

RECOMMENDATION (as planned): norm/invariant/height experimentation is FROZEN.
The one cheap decisive check left is the EHS parameter range.  If `delta <= n+1`
is forced, PARK M7.12 with `Docs/36` Thm 36.2 as the stopping theorem and
`Docs/37` as the record of why the product route cannot help either.

M7.13. EHS PARAMETER-RANGE CLOSURE -- AND M7.12 PARKED (DONE).  Wrote
`Docs/38_EHS_CLOSURE_AND_PARKING.md`.

(a) THE CHECK.  `Docs/37` Thm 37.3 returns exponent `n+1` and would close the
    problem iff `delta > n+1`, so the whole of Track B reduced to one question
    about EHS:  is `delta > n+1` available?  ANSWER: NO.
    Eisenbrand--Haeberle--Singer fix `delta = 1/(3n)` in the successive-minima
    step.  Their Theorem 4 assumes `Q >= (2n||x||_inf)^{3/2}` and
    `lambda_1 >= Q^{-1-1/(3n)}`, and concludes `|E| >= Q^{-(n+1)}`.  The
    `1/(3n)` is FORCED by the geometry-of-numbers step -- the first successive
    minimum must clear `Q^{-1}` by that margin before the Subspace-Theorem
    machinery applies -- so it is not an unoptimised free parameter.
    Independently confirmed at the abstract level (arXiv:2312.02057): the
    headline is `|E| >= gamma (n max_i|x_i|)^{-2n}`, improving the previous
    `(n max_i(|x_i| sqrt(a_i)))^{-2^n}`, with `gamma` NON-EXPLICIT from the ST.

    M7.13 (EHS parameter-range closure).  No available EHS parameter
    optimisation converts the exponent `n+1` of `Docs/37` Thm 37.3 into the
    required exponent `< 1`.   [PROVED]

(b) WHY THERE WAS NEVER ANY SLACK.  `n+1 = 1 + n` is exactly
    `q^{-1} x q^{-n}`: the first coordinate's EXACT contribution (Docs/37
    Fact 37.2, `||Q_k sqrt(a_1)|| = eps^{-k}`) times the best each remaining
    coordinate can do (exponent 1, Docs/36 Thm 36.1).  So `Docs/37`'s `n+1` is
    the natural arithmetic scale, not a loose estimate -- which is what makes
    this a closure rather than a missed opportunity.

(c) M7.12 IS PARKED.  NOT because the problem is experimentally exhausted, and
    NOT because M7.12 is proved impossible.  Parked because the remaining
    obstruction is isolated and every major mechanism is excluded:
        norm / product-of-conjugates   Docs/36 Thm 36.2  (structural, all n>=2)
        single-index improvement       Docs/36 Thm 36.2  (degree-4 forced)
        two-index improvement          Docs/36 Thm 36.1  (bounds a SUM)
        heights / Liouville            Docs/34 Thm 34.1
        ST hyperplane heights          Docs/34 sec 34.5  (circular)
        Ridout / ST power laws         Docs/35 sec 35.6
        recurrence invariant           Docs/36 sec 36.3  (I_k == 0)
        conjugate amplification        Docs/36 sec 36.6  (measured flat)
        product-level Pell coupling    Docs/37 Thm 37.3  (one coord at exp 3)
        LRS literature route           Docs/37 sec 37.1
        EHS parameter optimisation     Docs/38 M7.13
    NOT CLAIMED: that M7.12 is impossible, or that no theorem exists.  Only
    that this toolkit cannot reach it and the missing input is a different KIND
    of argument.

(d) RESEARCH CONCLUSION.
        SRS is reduced to a genuinely new small-fractional-part problem for a
        trace-zero quadratic recurrence:
        { k : ||(sqrt(a_j)/2 sqrt(a_1))(1+sqrt(a_1))^k||
                 <= (1+sqrt(a_1))^{-k/(3n(n-1))} } finite, k* log eps <= poly(N),
    where the sequence is trace-zero order-2 (Docs/35 Thm 35.3) in a fixed
    biquadratic field and `beta_j NOT in Q(sqrt(a_1))` (Docs/35 Thm 35.2)
    removes the automatic-Pisot mechanism that would make the statement false.

(e) PROGRAMME CLOSED.  `Docs/21`-`Docs/38`.  M7.1-M7.4 (Docs/19,21,22); L1-L7
    (Docs/23-29); L8-L11 (Docs/30-33); L12 (Docs/34); L12' (Docs/35);
    L12'' (Docs/36); L12''' (Docs/37); M7.13 (Docs/38).  Norm, invariant and
    height experimentation is FROZEN.  The next action for this line, if any,
    is EXTERNAL -- a literature or expert question about small fractional parts
    of trace-zero quadratic recurrences -- NOT further algebraic manipulation of
    `||beta_j eps^k||` by norm, lattice or height methods.

CLEANUP PASS (DONE).  9 pyflakes findings fixed (`two_small.py`,
`two_small2.py`, `enemy_gen.py`, `tests/test_proofs.py`); pyflakes now clean
tree-wide.  The `f1_theorem` import guard became a real assertion rather than a
`# noqa` that pyflakes does not honour.  `PROOFS.md`: malformed markdown table
row fixed (stray 2-space indent on the `| Theorem 2.5 |` row broke the table)
and an orientation addendum added recording that the file covers Papers 1-19
while `Docs/21`-`Docs/38` reduced and then parked the effectivity question.
Latent flaw fixed: `two_small.py` ended with a bare module-level `runs()`, so
`import two_small` ran a full combinatorial search and never returned; now
gated behind `if __name__ == "__main__":` (nothing imports it, so
behaviour-preserving).  Suite: 379 passed, 13966 subtests.

M7.14. THE RECURRENCE-DEFECT ROUTE, ASSESSED (DONE -- CLOSED).  Wrote
`Docs/39_RECURRENCE_DEFECT.md`.  M7.12 is closed again, with a 12th mechanism
excluded.

(a) THEOREM 39.1.  `R_k = eta_{k+2} - 2 eta_{k+1} - eta_k = -(m_{k+2} - 2m_{k+1}
    - m_k) in Z`, `|R_k| <= 2`, and `R_k = 0` forces `m_{k+2} = 2m_{k+1} + m_k`
    (nearest-integer sequence shadows the exact Pell recurrence).  EXACT,
    integer-only.  A float computation APPEARS non-integral for k > ~40; that is
    cancellation, not mathematics.  Engine `m_k` verified against an independent
    400-digit reference at k=40,100,200,296,500,1000 -- all exact.

(b) STEP 1 IS EMPIRICALLY VOID.  k=1..298: `R_k in {-1:76, 0:159, 1:62}`.  The
    shadowing holds MORE THAN HALF THE TIME; `R_k = 0` in 159/298; `R_k =
    R_{k+1} = 0` in 69, last at k=296.  So exact shadowing is GENERIC, not
    exception-triggered, and carries essentially no information about `d_k`.

(c) STEP 3-4 IS **FALSE**, NOT MERELY UNPROVED.  THEOREM 39.3: no positive lower
    bound on `|beta_j - A|`.  `A = (m_{k+1} - eps' m_k)/(eps^k(eps-eps'))` lies
    in `Lambda = (1/(4a_1)) Z[sqrt(a_1)]` (denominator norm `+/-4a_1`, since
    `P_k^2 - a_1 Q_k^2 = (-1)^k`), and LEMMA 39.2: `Lambda` is DENSE in R, so
    `dist(beta_j, Lambda) = 0`.  Proof: write `4a_1 x = n + theta`; for any N,
    density of `{q sqrt(a_1) mod 1}` (pigeonhole on `0, sqrt a_1, ..., N sqrt a_1`
    mod 1) gives `q` with `||q sqrt a_1 - theta|| < 1/N`, hence a lattice point
    within `1/(4 a_1 N)` of `x`; let `N -> infinity`.  Four lines.

(d) THE QUANTITATIVE VERSION ALSO FAILS.  `A`'s coefficients are bounded by
    `~eps^{2k}`, so the best separation is `~eps^{-2k}`, while the shadowing error
    is already `eps^{-k}` -- SAME ORDER, compatible.  THEOREM 39.4: best case
    gives `d_{k+1} + |eps'| d_k >= eps^{-k}`, exponent 1 -- EXACTLY `Docs/36`
    Thm 36.1.  NO EVASION of Thm 36.2.

(e) WHAT SURVIVES -- A DYNAMICAL REFORMULATION.  `R_k in {-1,0,1}` is a
    finite-alphabet digit sequence; `(eta_k, eta_{k+1})` is the orbit of an
    explicit piecewise-linear EXPANDING map `T(u,v) = (v, 2v+u+l(u,v))` on
    `[-1/2,1/2]^2`, with `l` the forced digit.  M7.12 becomes: does this orbit
    enter `{|u| < eps^{-ck}}` infinitely often?  A shrinking-target question for
    a PL expanding map -- genuinely a different KIND of input from everything
    Docs/34-38 excluded.  A REFORMULATION, not a proof; no answer claimed.
    THIS IS THE CONCRETE OBJECT TO RAKE with someone working on shrinking targets
    for PL maps or Diophantine approximation for quadratic recurrences.

(f) METHODOLOGY WARNING, SECOND INSTANCE.  Twice reported a small POSITIVE
    `delta` that was an artifact of a bounded search (`|p|,|q|<=400` gave
    `6.44e-5`; restricting `q` to Pell denominators plateaued at `0.1010`).  A
    `delta` that is SMALL is not a `delta` that is POSITIVE.  Also predicted
    `l_k = 0 => k = O(1)`, which is FALSE: a single `l_k = 0` does not pin `A`,
    and even two consecutive do not contradict, because `delta = 0`.  All
    failures were numerical-or-scope, not algebraic; each caught only by an
    independent check.

## M6 — rank-3 effective constant (the only active mathematics target)

We now aim at proving or refuting an effective constant at exactly `r = 3`.
There is no M7 yet.

- [X] **M6.1  Formalise the `F(C)=±1`/unit branch.**  (Done: `Docs/15_M6_1_UNIT_REDUCTION.md`.)  Key correction: `F(C)=±1` does NOT put `u ∈ O_K^×`; the integral norm-one unit in `E = Q(√(d₁d₃),√(d₂d₃))` is `ρ := d₃·u² ∈ O_E`, with explicit integer expansion (10) and `N(ρ) = H(C)² = ±1` exactly.  Rank of `O_E` is 3, and the branch `H=±1` reduces SRS rank-3 to bounding a three-term linear form in fundamental units (§15.7). Next: M6.3 insert the strongest three-log effective bound and record the real `A_r`.
- [X] **M6.2  Small `L` as an explicit logarithmic form.**  (Done: the exact form
  `Σ_i a_i·log|ε_i(σ_t)| = log|ρ|_σ` (t=1..4) is derived in `Docs/15` §15.5–15.7, with `ρ = d₃u²` the norm-one unit and `a ∈ Z³` the unique solution of the 3×4 regulator system (19).  Height/exponent control is `Docs/15` §15.6 (21)–(22).)


- [X] **M6.3  Bound every algebraic height occurring.**  Done in `Docs/15` §15.6:
  `h(ρ) ≤ log d_max + 2 log M + log 6` (21), and the unit-exponent vector `a`
  is controlled by (22): `|a_i| ≤ μ(E)·(log(1/(d₃|L|²)) + log d_max + 3log M + O(1))`,
  with `μ(E)` a fixed residue of `1/R_E` — explicit once the fundamental units
  are fixed.

- [X] **M6.4  Quantitative audit of the three-log bound.**  Computed in `Docs/16`:
  using a Matveev-type form with D=4, n=3, the parameter sale is
  `𝔇 = C·4⁵·(1+log B)·(Π h(ε_i))·polylog(heights)`.  Under the elementary
  norm bound, `log B ≤ log₂μ(E)+O(log N)`; the two potential blow-ups are the
  fundamental-unit heights `h(ε_i) ≤ O(√d_max)` (real-quadratic
  regulator worst case) and the regulator-inverse factor μ(E) from M6.1 (22).
  The audited M6.4 figure is `𝔇 ≤ 2^{O(N)}` (roughly
  `C·4⁵·(1+log B)·(Σ_i h(ε_i))`), i.e. the three-log *constant* is boundedly
  exponential in `N`.  It yields the **lower bound `E_L ≥ e^{-𝔇}/2` on `E_L`,
  not the SRS upper bound `E_L ≤ poly(N)`; the usable upper bound remains the
  Paper-9 norm bound `E ≤ 6N+O(1)`.  There is no `T(3,δ)` improvement from
  the three-log route, and it reproves at most the fixed-rank bound in a much
  weaker form.

- [X] **M6.5  Exact source of the super-polynomial loss.**  Decided: it is
  the product `Π h(ε_i)` of fundamental-unit heights (real quadratic / Pell
  worst case), fed through M6.1 (22) into Matveev's constant.  The precise
  sufficient condition for the Baker form of the argument to close `T(3,δ)` is

  ```
  max_i h(ε_i)  ≤  N^{O(1)}            (♦)
  ```

  i.e. the fundamental units must have *logarithmic height polynomial in the
  SRS bit-size* — a much stronger statement than any standard regulator
  estimate, and not true for generic real quadratic triples (worst-case
  regulators are `~ exp(c·√dmax)`, i.e. `h(ε_i) ~ √dmax = 2^{N/2}`).  Whether
  the special `H=±1` fields satisfy (♦) is the precise M6.5 hole.

- [X] **M6.6″  Search for a biquadratic specialization (closed, no result).**
  Re-examined whether the `(Z/2Z)^2` Galois structure of `K` forces the
  trilinear form to degenerate to a two-logarithm one (Størmer/Mitsui-type).
  It does not, on the available explicit estimates: the hook `ρ = d₃u²`, a
  norm-`+`unit in `E`, does not mention `(Z/2Z)^2` beyond fixing the
  four-term logarithmic embedding above.  Closed.

- [X] **M6.7′  T(3,δ) attempt (closed, no polynomial result).**  The
  originally planned weakening
    T_3^(1): log₂(1/c) ≤ 2^{O(N)}
    T_3^(2): log₂(1/c) ≤ exp(O(log²N))   [intermediate]
    T_3^(3): log₂(1/c) ≤ N^{O(1)}         [the target]
  was pursued; only the trivial member T_3^(1) is consistent with what we
  have, and the intermediate T_3^(2) is not achieved by the audited
  Matveev scaling (which yielded `E_L ≥ e^{−𝔇}/2` with `𝔇 ≤ 2^{O(N)}`, i.e.
  no `T_3^(2)`).  Closed; not reopened.

- [X] **M6.6  Rank-3 branch closed.**  We record that this branch was audited and that:
  (1) the exact unit reduction was verified (Docs/15); (2) the Matveev/trilinear route was audited
  in detail (Docs/16) and found to give no `T(3,δ)` improvement over the elementary bound,
  because the available constants scale like `𝔇 ≤ 2^{O(N)}` through `Π h(ε_i)` and the
  regulator-inverse `μ(E)`; (3) the effective polynomial separation at fixed `r=3` remains
  the full-norm statement of Paper 9 (`E ≤ 6N+O(1)`).  We do **not** treat the fulfilled
  sufficient inequality `max_i h(ε_i) ≤ N^{O(1)}` (♦, M6.5) as an active target: it is the
  same scale of difficulty as an effective regulator theorem, and the surrounding problem
  at `r=3` is already solved effectively.  The branch is closed.

  Diagnostic summary of M6 (one sentence also added as `Docs/17_M6_LESSONS.md`):
  the obstruction in the rank-3 algebra is not the approximation exponent; it is the
  *field-dependent effective constant* in the unit/Baker step, and generic analytic
  number-theoretic machinery does not, on the available estimates, produce a
  polynomial-scale value for it.

## RECONSTRUCTED — removed audit/done sections (NOT verbatim)

During the M6 splice the large review sections "Honesty pass" and
"Done — Python / Rust" were accidentally deleted.  They are **not**
recoverable verbatim; what follows is a compressed, provenance-labelled
replacement compiled only from facts still present in the repository and
in surviving grep output.  Items marked `(lost)` were once present but are
not re-derivable from context.

### (former) Honesty pass — claims that outran the code

A review found ten places where prose/docstrings claimed properties the code
did not have.  All ten were fixed.  Two (lost in full):
- `squarefree_decompose` docstring claimed "verified squarefree"; no such test
  existed.  Fixed: factoring now computed by exact path; squarefreeness claim
  downgraded/verified per path.
- `random_squarefree`/`verify_theorem`/`generate_f5`/`intervalsum` had
  doc-code mismatches, fixed (fragments recorded).

### (former) Done — Python (`SRS/`)

Files present and referenced by the test suite / M-work:
- `families.py` (F1–F6 generators, `input_bits`) — remains in use.
- `field.py` (multiquadratic field/rank) — remains in use.
- `squarefree.py`, `signexact.py`, `intervalsum.py` — exact kernels.
- `f1_theorem.py` — F1 residue + separation law.
- `refutation.py` — counterexample harness.
- `enemy_gen.py` — M2 embedding-profile generator (NEW, added in this phase).

### (former) Done — Rust (`Rust/srs/`)

- `field.rs` — multiquadratic field over `BigInt`
- `squarefree.rs` — exact decomposition
- `signexact.rs` — certified sign decision
- `f1.rs` — Pell family, residue law, separation law

This block is intentionally short; the full original audit trail is lost.

---

