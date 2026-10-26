# SRS — The State of the Art, Corrected

**Read this first.** It records the literature check that changed the answer,
and the error that check exposed.

| | |
|---|---|
| **Written** | after checking the literature (arXiv:2312.02057; SoCG 2024) |
| **Corrects** | `PROOFS.md` §2.8 as originally written; the first version of `THE_MISSING_THEOREM.md` §3 |
| **Bottom line** | the conjecture is **not** blocked on a missing `ε`-improvement. That improvement exists, it is Schmidt's Theorem 1B, and it has already been applied to this exact problem. What remains open is the **effectivity of one constant**. |
| **Where the current account lives** | `THE_MISSING_THEOREM.md` §3 (the open problem, restated correctly), `PROOFS.md` §2.8, **§2.9** and §4 |
| **Superseded again by Thm 2.9** | `PROOFS.md` §2.9 applies the same theorem to the *single* linear form SRS is, getting exponent `r−1` instead of EHS's `2n` — so the EHS bound recorded in §1 below is **not** the sharpest available. Thm 2.9 matches the achieved rate. This file remains the record of the correction; the sharp bound is Thm 2.9. |

Scope note: this file is the *record of the correction* — what was claimed,
what is true, and what was verified. It is deliberately narrow. For the
current statement of the problem, the bounds, and the remaining gap, read
`THE_MISSING_THEOREM.md`; for the proofs, read `PROOFS.md`.

---

## 0. What I got wrong, in one paragraph

My previous turn stated, as the "exact missing theorem":

> an `ε`-improvement for a **single** algebraic linear form in `r ≥ 3`
> variables — Roth is two-variable, Schmidt's `n`-variable exponent needs
> `n` **separate** forms, and SRS has only **one** form to offer.

**Every clause of that is wrong except the arithmetic.** The counting error:
Schmidt's Subspace Theorem, in its `n+1`-form form, is a statement about a
**product** of forms. Applied to variables `(q, p_1, …, p_n)` with

```
L_0 = q·β_1 − p_1,   L_1 = q·β_2 − p_2,   …,   L_{n−1} = q·β_n − p_n,   L_n = q,
```

the product is `q · Π_i |qβ_i − p_i|`, and with `p_i = O(q)` this is exactly
Schmidt's Theorem 1B. The `n−1` coordinate forms are **supplied by the
reduction itself** — they are free. Only `L_0` has to come from the problem.
The price is the `q^{1+δ}` factor and `dist_Z` in place of a bound on `|L_0|`,
and that price is affordable.

I had the right theorem in hand and miscounted the forms it needs. That is
the whole error.

---

## 1. The result that settles the exponent

**Eisenbrand, Haeberle, Singer — "An improved bound on sums of square roots
via the subspace theorem", arXiv:2312.02057, SoCG 2024.**

**Theorem (their Thm 8).** Let `a_1, …, a_n` be distinct squarefree positive
integers and `E = Σ_i x_i sqrt(a_i) ≠ 0`, `x ∈ Z^n`. There is a constant
`γ > 0` depending only on `a_1, …, a_n` such that

```
|E|  ≥  γ · ( n · ‖x‖_∞ )^{−2n}.
```

With `m = log₂‖x‖_∞` and `N = n·m` (total input size), this is

```
E  ≤  2n·(log₂ n + m)  +  log₂(1/γ)  =  O(N)  +  log₂(1/γ).
```

**Remark 9 of the paper** sharpens the exponent `2n → n + ε` for any `ε > 0`
(at the cost of `γ` depending on `ε`), giving `E ≤ N + O(n log n)`, i.e.
`E/N → 1`.

**And the bound above is not the sharpest available.** `PROOFS.md` Thm 2.9
applies the same Subspace Theorem to the **single** linear form that SRS is,
rather than to `n` simultaneous approximations, and gets exponent `r−1`
instead of `2n`. At `r = 3` that is `E/N ≤ (r−1+δ)/r = 0.667 + δ/3` against
this bound's `1.0` and Thm 2.6's `1.5` — and it **equals the achieved rate**,
so the sharp conjecture is the *asymptotic* `limsup_{N→∞} E/N = (r−1)/r`. See
§2 below.

> **Retracted (`sup` form).** This file previously wrote the sharp conjecture
> as `sup E/N = (r−1)/r`, a supremum over all `N`. That is **false**: a
> certified `r = 3` instance gives `E/N = 1.2795 > 2/3`
> (`PROOFS.md` §2.10.2, `d=(2,3,13)`, `C=(−10,−98,51)`). The excess is a
> bounded additive slack over a small `N`, and vanishes as `N → ∞`.

**Correction — EHS Thm 4 does not give `E = O(N)`.** Their Thm 4 (and the
`PROOFS.md` table row derived from it) is recorded as "per-instance
computable, `E/N = O(N)`, at the cost of a checkable hypothesis". That is
misleading in two ways, and the correction matters because this is the **only
effective polynomial bound in the literature**. Read from the paper
(arXiv:2312.02057, HTML v1) rather than from a summary:

- `Q` is a **free integer parameter**, not an input. The lattice is generated
  by `B` with `det = 1/Q^{n+1}`; Thm 4's hypothesis is
  `Q ≥ (2n‖x‖_∞)^{3/2}` *and* `λ₁ ≥ 1/Q^{1+1/(3n)}`. Only the second is
  substantive, and Thm 8 obtains it for all `Q ≥ Q₀` with `Q₀` depending only
  on the `a_i` — **via the Subspace Theorem, hence ineffective**.
- So the cost hides in `Q₀`. The "bad" `q` are those with
  `q^{1+δ}·Πᵢ dist_Z(qβᵢ) < 1`. Since `|q²a − p²|` is a nonzero integer,
  `dist_Z(qβᵢ) ≥ 1/(2qβᵢ+1)`, giving
  `q < (Πᵢ(2βᵢ+1))^{1/δ}`; with `δ = 1/(3n)` that is `exp(Θ(n²a))` for
  `a = max bitlen(aᵢ)`. Failure of the `λ₁` condition *also* requires
  `q ≥ Q^{1+δ}`, so `Q₀ = exp(Θ(n²a))` and the honest bound is

  `E = O(n²a + n·m) = O(n·N)`,  `m = log₂‖x‖_∞`,

  a factor `n` worse than the target — and `n` can be `Θ(N)`.

- **The obstruction is sharp.** The constant in `dist_Z(q√a) ≥ c/q` is
  `c ≈ 1/(2√a)`, i.e. `a/2` bits, amplified by `Θ(n²)`. The tempting
  `dist_Z(qβ) ≥ 1/q` is **false**: `dist_Z(7√2) = 0.101 < 1/7`. So this route
  provably does not close the conjecture.
- **"Computable hypothesis" conflates two objects.** `λ₁` is computable, but
  checking it certifies a *single* `Q`; the theorem needs the *threshold*
  `Q₀`, which is the ineffective part.

| | this file previously said | corrected |
|---|---|---|
| Thm 4 bound | `E/N = O(N)` | `E = O(n·N)` once `Q₀` is made explicit |

**Unverified, flagged rather than claimed.** From Minkowski II with
`λ₁ ≥ 1/Q^{1+1/(3n)}` the per-minima bound is `λ_i ≤ 1/Q^{n−1/(3n)}` for
`i ≥ 2`, not the paper's `1/Q^{2/3}` (these agree only at `n=1`). The
conclusion survives and improves: `q ≤ Q^{1+1/(3n)}` rather than `Q^{n+1/3}`,
hence `|E| ≥ 1/(2Q^{1+1/(3n)})`. This has **not** been checked against the
published PDF and is not asserted here.

**So the exponential-in-`n` barrier is gone.** Measured against this repo's
Thm 2.6 (bits, `m = 100`):

| `r` | 3 | 5 | 8 | 10 | 15 | 20 |
|---|---|---|---|---|---|---|
| Thm 2.6 (norm) | 600 | 2400 | 1.9e4 | 7.7e4 | 2.5e6 | 7.9e7 |
| EHS, exponent `2n` | 610 | 1023 | 1648 | 2066 | 3117 | 4173 |
| ratio | 0.98 | 2.4 | 12 | 37 | 788 | 1.9e4 |

The ratio grows like `2^{r−3}/r`. **Thm 2.6 is no longer the best known bound
in this area** — it remains the sharp statement about what the *norm* can do,
which is a different and still worthwhile claim.

Caveat, stated because it matters at small `r`: at `r = 2, 3` the two bounds
are comparable and Thm 2.6 is in fact **sharper**, because EHS does not know
the radicands. The `r = 2` constant `1/2` is this repo's own result and should
not be attributed to the literature.

---

## 2. A sharper application than EHS: the single form

Worth flagging separately, because it means the EHS bound quoted above is not
the best bound obtainable from the same theorem. SRS is a **single** linear
form

```
L(C_1,…,C_r)  =  C_r + Σ_{i<r} C_i·sqrt(d_i/d_r)
```

in `r` variables, and a single form in `r` variables has Subspace-Theorem
exponent `r−1` — the `r` forms being `L_0 = L` and `L_1..L_{r−1} = C_1..C_{r−1}`
(`C_r` is not a form and is not needed). EHS instead approximate the `r−1`
numbers `sqrt(a_i/d_r)` *simultaneously* with a common denominator, which is a
strictly more expensive shape and costs exponent `2n`.

| | exponent on `‖C‖_∞` | `E/N` at `r = 3` | vs ceiling `2/3` |
|---|---|---|---|
| Thm 2.6 (norm) | — | `1.5` | `2.25×` |
| EHS Thm 8 (simultaneous) | `2n` | `1.0` | `1.50×` |
| **Thm 2.9 (single form)** | **`r−1`** | **`0.667 + δ/3`** | **`1×` as `δ → 0`** |
| achieved, Thm 2.5 | `r−1` | `0.667` | — |

Proved in `PROOFS.md` §2.9, corroborated by
`tests/test_proofs.py::TestSingleFormSubspaceTheorem`. At `r = 2` the same
application is literally Roth, and the repo's Thm 2.2/2.6 confirm it sharp at
`1/2` with an explicit constant — so it is known-correct at the one rank where
the answer is owned.

## 3. The matching upper bound — and a rediscovery

**Theorem (their Thm 10).** For every `L ≥ 2` there is `x ∈ Z^n`, `x ≠ 0`,
`‖x‖_∞ ≤ L`, with

```
0 < |E|  ≤  n · max_i sqrt(a_i) / L^{n−1}.
```

Compare this repo's **Thm 2.5**:

```
|S|  ≤  ( Σ_i sqrt(d_i) ) / B^{r−1}.
```

`Σ_i sqrt(d_i) ≤ n · max_i sqrt(d_i)`, and the exponent is `n−1` in both.
**These are the same theorem.** Thm 2.5 was rediscovered here independently,
by sorting `(B+1)^r` sums instead of applying the pigeonhole principle to a
lattice. That agreement is the substance of §2.7's "tight at `r = 2`" remark,
and it is a genuinely good sign for the repo's own proofs.

Their §3 also gives the near-tightness statement: the exponent `n−1` is
asymptotically optimal for fixed `a_1, …, a_n`.

---

## 4. The actual remaining gap: effectivity of the constant

This is the whole problem now, and it is a much narrower target than the one
I gave you last turn.

The bound is stated for **fixed** `a_1, …, a_n`, with `γ` a function of them.
**In SRS the `a_i` are part of the input.** The split is therefore

```
E  ≤  2n·(log₂ n + m)  +  log₂(1/γ(a_1,…,a_n))
   \___________________/     \_________________________/
    polynomial in N — PROVED     NOT PROVED   ← the entire gap
```

> ### The open problem, precisely
>
> Bound `Q_0 = Q_0(a_1, …, a_n)` — the threshold in EHS's proof — from
> above so that `log₂ Q_0 ≤ poly(N)`.
>
> Equivalently: bound the number of positive integers `q` satisfying
>
> ```
> q^{1+δ} · Π_i dist_Z( q · sqrt(a_i) )  <  1.
> ```
>
> The Subspace Theorem is **ineffective** here. EHS state in their Discussion
> (§4) that the quantitative versions — including Evertse–Ferretti (2013) —
> "do not provide such bounds either … [t]his is still the case when all
> algebraic numbers are square roots of integers." They name it as the open
> problem, and note that an explicit bound "would immediately apply to a
> separation bound for the sum of square roots."

**Why this is not already a proof that SRS ∈ P.** Deciding `sign(E)` by
approximation requires a *known* lower bound on `|E|` to set the precision
threshold. `γ` is unknown and depends on the input, so no threshold can be
set. EHS is a **separation bound**, not a complexity-class result. The
exponential-in-`n` precision of the classical norm bound is replaced by a
polynomial-in-`n` precision *with an unusable constant*.

---

## 5. Where this leaves the repo's own results

Nothing in §2.1–§2.6 of `PROOFS.md` is invalidated, but its role changes.

| result | status after the correction |
|---|---|
| Thm 2.1 (Liouville) | correct, and now visibly the classical `2^n` bound |
| Thm 2.3 (orbit second moment) | correct; best *norm* bound before Thm 2.6 |
| Thm 2.4 / Prop 2.4′ (tower, support rank `j−2`) | correct, including the retracted rank claim |
| Thm 2.6 (global-flip pairing) | **still the sharp statement about the norm**, and the best fully *effective* bound here — but not the best bound overall: see Thm 2.9 in the next row |
| **Thm 2.9 (single-form Subspace Theorem, this repo)** | **the sharpest bound available**: exponent `r−1`, `E/N ≤ (r−1+δ)/r`, equal to the achieved rate at every rank |
| Thm 2.5 (pigeonhole) | **= EHS Thm 10**, rediscovered independently |
| Thm 2.2 + Thm 2.6 at `r=2` (`= 1/2`) | still the repo's own sharp result; **stronger than EHS at that rank**, which does not know the radicands |
| Thm 3 / Cor 3.3 (F1) | unaffected |
| §2.8 as originally written | **withdrawn** — see §0 |

The `2^r` in Thm 2.6 is a property of the **norm**, not of the problem. That
distinction is the single most important thing this correction fixes: the
previous version of the file treated the norm's barrier as the problem's
barrier, and it is not. Nothing in §2.1–§2.7 of `PROOFS.md` is invalidated —
only its role changes, from "the obstacle" to "the sharp effective bound".

---

## 6. What a proof would have to contain, revised

1. **Not** a better norm. Thm 2.6 is its optimum; the `2^r` is intrinsic to
   taking the norm of the whole `2^r`-element orbit.
2. **Not** an `ε`-improvement for one algebraic linear form. That exists, it
   is Schmidt's Subspace Theorem, and it is already applied — Thm 2.9 reaches
   the ceiling with it (§2).
3. It **would** be an effective constant: `log₂(1/c(d,δ)) ≤ poly(N)` in
   Thm 2.9, or equivalently `log₂(1/γ(a_1,…,a_n)) ≤ poly(N)` in EHS Thm 8
   — a bound on the number of `q` with `q^{1+δ} Π_i dist_Z(q sqrt(a_i)) < 1`.
   Plus the explicit description of Thm 2.9's exceptional rational
   subspaces, which are *detectable* by exact linear algebra, so the generic
   case is already constant-free. See §4 and
   `THE_MISSING_THEOREM.md` §3.2 for the full statement.
4. Plausibly something exploiting that every `β_i` is a **square root of a
   rational** — EHS note the quantitative Subspace Theorem fails even in this
   special case, so the square-root structure is exactly where the missing
   ingredient would have to live. This is now the most concrete lead in the
   repo, and it is a lead the literature has flagged and not solved.

Restated at greater length, with the precision argument, in
`THE_MISSING_THEOREM.md` §3 and §5.

---

## 7. Provenance

- The EHS result was found by searching the literature during this session and
  reading the paper (arXiv:2312.02057, HTML full text). The bound, the
  `n+ε` sharpening, the sharpness statement, and the inefficacy remark are
  quoted from that paper; the tables comparing it to Thm 2.6 and the
  identification of Thm 2.5 with their Thm 10 are this repo's own arithmetic.
- Roth (1955) and Schmidt's Subspace Theorem / Theorem 1B were checked.
  Davenport–Schmidt (1967) was checked and **excluded** — it concerns
  approximating a real number by quadratic irrationals, not the
  `ε`-improvement at issue.
- Nothing about the *current* best quantitative bound on the number of
  exceptional `q` was verified beyond what the EHS paper itself reports. If
  this file is used to plan work, that specific question should be re-checked
  against the literature as of the date of use.
- **Thm 2.9 (`PROOFS.md` §2.9) is this repo's own**, not quoted, and is
  sharper than the EHS bound recorded in §1. It follows from the same
  Subspace Theorem applied to a single linear form rather than to `n`
  simultaneous approximations. It is corroborated by
  `tests/test_proofs.py::TestSingleFormSubspaceTheorem` and, on the achieved
  side, by family F6 in `families.py` with `refutation.py` searching for a
  counterexample to the *asymptotic* rate and finding none. (The `sup` form of
  the conjecture is retracted — see `PROOFS.md` §2.10.2, where a certified
  `r = 3` instance attains `E/N = 1.2795 > 2/3`.)
- The retraction in §0 is a real error in this repo's prior documentation, not
  a hypothetical. `tests/test_proofs.py::TestSubspaceTheoremRoute` now pins
  the corrected statement, including a test that Thm 2.6 is *dominated* at
  large `r` so the false claim cannot be reintroduced.
