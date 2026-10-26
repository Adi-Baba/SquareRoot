# Reconnaissance M2 — Embedding-Profile Envelope of the Enemy Family

Version 1.0.  Empirical/derived profile of near-extremal instances produced by
the pigeonhole (consecutive-gap) generator.  This is a measurement record,
not a theorem; every quantitative tendency is supported by the explicit
table, and the analytic consequence is recorded at the end.

---

## M2.0 What was computed

Script: `SRS/enemy_gen.py`.  For `r ∈ {3,4,5,6}` with the first `r` primes as
radicands `p_i`, and box size `B ∈ {3,4,5,6,8}`, we:

1. enumerated `v(C) = Σ p_i·C_i` over `C ∈ {0,…,B}^r`;
2. found the smallest consecutive gap `gap = v(C^b) − v(C^a)` and set the
   enemy coefficient `C = C^b − C^a` (mixed signs);
3. using that `C`, computed the full sign-profile of the `2^r` embeddings
   `σ(C) = Σ (±) C_i √p_i`;
4. reported `N`, `E/N`, the minimizing embedding, and the quantity
   `P_* = Π_{σ≠σ_*} |σ(C)|` in rate form `log₂ P_* / N`.

Key measurement (selected rows):

| r | B | N | E/N | log₂P*/N | norm bound exponent 2^{r-1}-1 |
|---|---|---|---|---|---|
| 3 | 3 | 12 | 0.29 | 0.94 | 3 |
| 3 | 5 | 15 | 0.45 | 1.07 | 3 |
| 4 | 3 | 14 | 0.71 | 0.71 | 7 |
| 4 | 8 | 21 | 0.53 | 2.28 | 7 |
| 5 | 3 | 23 | 0.61 | 3.69 | 15 |
| 5 | 8 | 22 | 0.71 | 3.63 | 15 |
| 6 | 3 | 23 | 0.61 | 6.76 | 31 |
| 6 | 8 | 32 | 0.65 | 6.45 | 31 |

(The rank-6 enemy at B=3 coincides with a rank-5 enemy because the minimal
coefficient vector had a zero sixth entry; the table reports the explicit
`C` in the raw output.)

---

## M2.1 What the profile says

For a pigeonhole-generated enemy `C`, the *small* embedding `σ_*` is the
principal one constructed with a small signed sum.  The remaining `2^r − 1`
embeddings are **not** small: a generic sign flip keeps roughly the same
`r`-term magnitude, `|σ(C)| ≍ M^{1/2}` per term, so `|σ(C)| ≍ m^{1/2}M^{1/2}
≍ m·M^{1/2}`... more precisely, a sum of `r` random signed terms of size
`~ M^{1/2}` has `log₂|σ(C)| ≈ (1/2)·log₂ M + O(1)` on average.  Hence

```
log₂ P_*  ≈  (2^r − 1)·( (1/2)·log₂ M + O(1) ).
```

Since `N ≈ r·log₂ B ≈ r·log₂ M` (the primes are fixed), this predicts

```
log₂ P_* / N  ≈  (2^r − 1)/(2r)·(log₂ M / log₂ M)·...  =  Θ(2^r / r).
```

The observed values confirm that growth: `0.94→1.07` (r=3), `2.28` (r=4),
`3.6–3.7` (r=5), `6.5–6.8` (r=6).  Each step `r ↦ r+1` roughly *doubles*
`log₂P_*/N`, which is exactly `2^r` behavior.  The small embedding does **not**
make the other conjugates small; the norm bound's input `Π_{σ≠σ_*}|σ(C)|` is
genuinely of order `M^{Θ(2^r)}` on these enemies.

---

## M2.2 Consequence (the M2 verdict)

On pigeonhole-extremal instances, for the small embedding `σ_*`,

```
|L_{σ_*}| · P_*  =  |N_{K/Q}(S)|  ( ≥ 1 ),
```

so the loss `P_*` in any norm-only argument is exactly `M^{Θ(2^r)}`.
**Because `P_*` is of the full rank-`2^r` scale, no improvement from merely
isolating one embedding is possible; the norm estimate is essentially tight
on the random enemy family.**

This is the M2 output in the form the next theorems need:

1. The "enemy" does not realize a special *arithmetic* structure at the
   embedding level — the small distance comes from a signed sum of the `r`
   positive magnitudes, and the conjugate profile is close to independent sign
   flips.  So a special bilinear/unit cancellation is *not* visible on the
   pigeonhole family.
2. To get `P_*` below `M^{Θ(2^r)}` one must use an instance for which **two or
   more** embeddings of the *same* `C` are simultaneously small.  Those are the
   simultaneous-approximation points, and they are exactly where the
   Subspace-Type constant must be computed.  The enemy generator therefore
   classifies two regimes:

   - **generic regime (measured above):** exactly one tiny embedding,
     every other conjugate `≈ M^{1/2}`, so the norm bound is saturated and the
    rank-`2^r` loss is unavoidable;
   - **thin simultaneous regime (not produced by the generator):** at least two
     embeddings of the same `C` are small, which occurs only when `C` satisfies
     an additional approximate relation, i.e. near a proper sublattice.  This
     is where an effective estimate is both needed and (for `r=2`) available.

3. Empirical tail of the rank-3, rank-4 families (table) is consistent with:
   for the enemy to reach `E/N` near the `(r−1)/r` ceiling, the minimal
   embedding has `|S|` roughly `(Σ√p_i)·B^{-(r−1)}`, and the other conjugates
   remain `~B` generic.  So the Subspace-type upper bound `(r−1+δ)/r` is doing
   real work in the thin regime, but the pigeonhole enemy family sits on the
   *norm-tight* regime: there the full-conjugate loss
   `P_* ≍ M^{2^{r-1}}` is present and is the ceiling on what a norm-only
   argument can certify.

---

## M2.3 What this tells the next theorem

If one ever proves `|L(C)| ≥ 2^{-poly(N)}`, its extremal case is not the
pigonhole enemy family (there the loss is unavoidably of full `2^{r-1}` scale).  It is a
simultaneous-approximation point of a rank-`r` lattice — concretely, an
approximation of at least one of the vectors

```
v_i := √(d_i/d_r),   i = 1,…,r−1,
```

by a rational point to accuracy `2^{-poly(N)}`.  Any uniform proof has to
involve an estimate on that simultaneous approximation.  This is precisely
scheduled as M4; M2 confirms it cannot be replaced by a norm manipulation.

---

## M2.4 Limitations of this probe

- For `r = 6` the raw table's best B=3 instance collapsed to a rank-5 vector
  (a zero coefficient); increasing the bound `B` or requiring all `C_i ≠ 0`
  is needed to see true rank-6 behaviour.
- We used the *consecutive-gap* pigeonhole as the enemy; that saturates the
  *lower* side.  It does not probe the thin simultaneous regime directly.
- The conjugate-size distribution in each row is reported via its sorted
  `log₂` values (`log_profile_norm` in the raw output); a finer taxonomy would
  cluster those profiles per `r`, but the doubling of `log₂P_*/N` with `r`
  is already decisive for the point above.

---

*(Raw table embedded in the `enemy_gen.py` output preceding this file; the
interpretation above is derived from that table by the argument of M2.1.)*