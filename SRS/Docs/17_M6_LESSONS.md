# M6 Lessons — the rank-3 audit, in one paragraph

Version 1.0.

The full audit of the rank-3 branch (Docs/15 for the reduction,
Docs/16 for the quantitative Matveev walk-through) yields one sentence:

> **The bottleneck at `r = 3` is not the approximation exponent and not the
> algebraic structure of `S`; it is the effective size of the field-dependent
> constant that the Subspace/Subspace-type route produces.**

Two exact consequences follow and should be carried into the rest of the
dossier:

1. **A polynomial bound at fixed `r = 3` already exists** as a pure
   algebraic-integrality (norm) statement, written down in
   `Docs/9_RANK3_NORM_FORM.md`: `E ≤ (9/2)·log M + (3/2)·log d₃ + O(1)`,
   i.e. `E ≤ 6N + O(1)`.  This is the reason the uniform SRS conjecture
   keeps returning to the *growing-`r`* regime: for fixed `r`, `2^{r-1}−1`
   is a constant and the norm argument already proves effectiveness.
2. **No argument that goes through a single logarithmic coordinate of
   `ρ = d₃u²` can improve this**.  As established in `Docs/16`, a
   three-log lower bound `|Λ| ≥ e^{−𝔇}` is in the direction
   `E_L ≥ e^{−𝔇}/2`, i.e. it does not give an upper bound on `E_L` at all,
   and the estimate for its effective constant reduces to the heights
   `Π h(ε_i)` and the inverse regulator residue `μ(E)`.  A sufficient
   closing hypothesis is (♦), and verifying that short of an effective
   regulator theorem is exactly the unsolved uniform target.

**Non-reopening clause.**  We do **not** treat (♦) as a new main
conjecture.  The fixed-`r=3` SRS theorem is available (Paper 9), so the
unit/Baker branch serves only as a diagnostic: generic multiplicative-unit
machinery does not carry the uniform-in-`r` constant.  That diagnostic is
what feeds the next module, M7.

---

*Reduced idea bank for M7:* any successful uniform proof must distinguish
itself from the closed list `full norm / norm tower / one-small-embedding
pigeonhole / two-small-embedding lemma / sign-pattern descent / generic
Subspace Theorem / generic Baker-Matveev`, all of which are documented in
`Docs/10`–`Docs/16` and in the `M7` closed-routes table in `TODO.md`.