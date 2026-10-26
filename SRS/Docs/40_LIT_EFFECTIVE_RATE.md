# Paper 40 — L13/T1: is there a published polynomial-rate theorem for `‖βε^k‖`?

Version 1.0 (first pass). Question (TODO L13 T1, from `Docs/34` §34.7 "NOT DETERMINED"):
does ANY published effective low-dimensional theorem give

```
‖β_j ε^k‖  ≥  k^{-C}     for explicit C = C(β_j, ε), all k ≥ K
```

for fixed algebraic β_j, ε (non-Pisot base)? This single theorem closes M7.12
(`Docs/28` Thm 28.3), hence SRS. This file records one row per candidate source.
Verdict rule: a "no" is a SEARCH RESULT, never a proof of nonexistence.

Empirical calibration from T3 (`m712_decay_probe.py`, this pass): on (2,3),
`η_k ≥ k^{-1.535}` over all k ≤ 10^5, record minima consistent (within fluctuations)
with iid scale 1/(2K). Truth-shape is C ≈ 1.5–2 if the metric picture is right.
Any theorem with explicit C > 2 and matching this scale would be sharp enough.

---

## Rows

**R1. Corvaja–Zannier, Acta Math. 193 (2004).** Qualitative: `{n : ‖α^n‖ < λ^n}`
finite for non-Pisot α (Thm 1); Main Theorem gives finiteness with height
conditions. **No rate.** Status: QUOTED via `Docs/29` §29.2 (verbatim).

**R2. Kulkarni–Mavraki–Nguyen, Trans. AMS (2018) (arXiv:1511.08525).** Operative
theorem for our tuple (q_1 = β_j fixed, non-rational). Delivers, via our audit
(`Docs/30` Thm 30.3/Cor 30.4): `{k : ‖β_jε^k‖ < θ^k}` finite for EVERY θ < 1,
hence `η_k = e^{−o(k)}`. **Non-effective; no k_0, no C.** Status: QUOTED verbatim
(`Docs/30` §30.1); proof audited (`Docs/31`): single ineffective sentence
("for all but finitely many elements of M"), Prop 3.4.

**R3. Quantitative Subspace Theorem line.** Schmidt 1989; Evertse–Schlickewei
2002; Evertse–Ferretti 2010/2013 (arXiv:1008.2340). Bounds the NUMBER t of
exceptional subspaces; heights/identities unbounded (`Docs/31` §31.1(d)).
**New confirmation this pass:** Evertse's survey arXiv:1008.2268 (abstract,
verbatim): for the small solutions outside the effectively determinable
subspace (Faltings–Wüstholz refinement), "Estimating the number of these
solutions is still an open problem. We give some motivation that this problem
is very hard." So even the COUNT of the terminal solution set is open in this
line — a height/rate version is further still. Status: CONFIRMED count-only.

**R4. Loxton–Van der Poorten programme.** Noubissie, arXiv:2504.09519 (Apr
2025): EFFECTIVE bound on the NUMBER of n with `|u_n| < (max|α_i|)^{n(1−ε)}`
for integer linear recurrence sequences. Count, not largest index — exactly
what `Docs/33` already gives for our set. Also integer-LRS only; our x_k is a
genuine two-term form (`Docs/37` §37.2). Status: checked and set aside
(`Docs/37` §37.1).

**R5. Chen–Ye–Zheng, arXiv:2604.14036 (Aug 2026), "Distribution modulo one of
linear recurrent sequences".** Criteria for finiteness of the set of limit
values of fractional parts of linear recurrences of REAL numbers; lower bounds
for the maximal distance between limit values. Generalizes Flatto–Lagarias–
Pollington and Dubickas. **Qualitative distribution theory, no rate, no
distance-to-Z lower bound.** Expanded version of §4 of arXiv:2511.21324v3.
Status: NEW this pass (arXiv API). Adjacent, not the target.

**R6. Accossato–Sanna, arXiv:2401.07661 (2024).** Number of distinct residues
mod m of integer second-order recurrences; corollary on the number of limit
points of `frac(ξα^n)`. Different problem (residue counts / limit-point
counts). Status: NEW this pass (arXiv API). Not the target.

**R7. Evertse–Győry, *Unit Equations in Diophantine Numbers* (effective unit
equations, rank ≤ 2–3).** NOT YET CHECKED against our shape. Note the shape
mismatch to resolve when checking: our object is `m ≈ βε^k` (approximation of
an orbit point by integers), not a unit equation x+y=1; CZ Main Theorem is the
right wrapper, and it is R2. Status: OPEN row.

**R8. Waldschmidt 2015 (effectivity of CZ Corollary 1).** Quoted in `Docs/29`
§29.2 as "one can bound from above the (numerators and denominators of the)
solutions x_1 and x_2" — that is effectivity for S-UNIT EQUATION SOLUTIONS,
not for the `‖qβ‖` threshold. Whether the quoted effectiveness transfers to a
`k^{-C}` rate is NOT DETERMINED. Status: OPEN row — must read the source.

---

## First-pass verdict

Four arXiv API queries ("fractional parts"+"linear recurrence"; Pisot+effective
+Subspace Theorem; "quantitative subspace theorem"; "powers of algebraic
numbers"+"nearest integer"), plus the R1–R4 rows already audited in `Docs/29`–
`Docs/37`:

> **NO published polynomial-rate theorem found.** The published frontier for
> `‖βε^k‖` with fixed algebraic data is exactly R2: `e^{−o(k)}`, non-effective.
> The quantitative ST cannot supply rates even at count level (R3, Evertse
> survey verbatim). Recorded as a SEARCH RESULT; R7/R8 remain open rows.

Consequence for the programme: the required input is a NEW theorem — a
polynomial-rate lower bound for small fractional parts of two-term
algebraic recurrences (or, in the M7.14(e) frame, a shrinking-target rate for
the PL expanding map T(u,v)). The T3 measurement says the true constant is
small (C ≈ 1.5–2, iid-like dips at scale 1/(2K)), so such a theorem is
plausibly within reach of a metric argument made effective for this algebraic
setting — which is precisely the "effective multidimensional Diophantine
approximation" gap named in `Docs/25` §25.7.

## Limits

- arXiv API search, not a MathSciNet/zbl sweep; older results (pre-1990
  Vijayaraghavan-type limsup statements) are covered qualitatively by the
  e^{−o(k)} discussion but were not row-checked.
- R5/R6 read at abstract level only; full texts not inspected.
- R7/R8 open.