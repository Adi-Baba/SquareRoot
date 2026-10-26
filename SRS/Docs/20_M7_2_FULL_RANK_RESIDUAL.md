# M7.2 — Full-Rank Residual Regime Audit (protocol)

Status: EXECUTED.  9104 residuals extracted (r=3..7, B=2..8).  The
residual regime is **explained as near-antipodal one-small geometry** —
not a new mechanism.  See `m7_2_residuals.csv` in `SRS/` for the full data.

**Precise wording:** this is *explained in the tested residual population*,
not a proof that all simultaneous-small configurations are near-antipodal.

---

## 20.1 The regime

Use canonical representatives: fix the first independent generator's sign
to `+` (one representative per Taù antipodal class `{σ, −σ}`).  For each
`C`, sort the canonical `2^{r−1}` magnitudes of

```
T_σ(C) := Σ_i C_i σ_i √d_i,     σ_i ∈ {±1}, σ₁ = +1.
```

Let the two smallest canonical magnitudes correspond to `σ⁰, σ¹`.  Their
tai-distance set

```
J := { i : σ⁰_i ≠ σ¹_i },
```

and `r_J := rank_{F2}{d_i : i ∈ J}`.

**Residual candidate:**  the two canonical magnitudes are both substantially
small, `|T_{σ⁰}|,|T_{σ¹}| ≤ 2^{-E/4}` (i.e. `Delta_2 ≥ E/4`),
and the support is *not* collapsed:

```
|J| ≥ r − 2,   r_J = |J|.
```

This is the M7.1 category C.  `M4.1` already says
`ε` cannot drop below the norm-inverse of the form on `J`, so the only
interest is in registration `(C,d,σ⁰,σ¹)`.

---

## 20.2 The `X/Y` decomposition (the key geometric step)

For a residual candidate define

```
X := (T_{σ⁰} + T_{σ¹})/2 = Σ_{i∉J} C_i·(σ⁰_i+σ¹_i)/2·√d_i,
Y := (T_{σ⁰} − T_{σ¹})/2 = Σ_{i∈J} C_i·(σ⁰_i−σ¹_i)/2·√d_i.
```

`X` is supported on the *agreement* coordinates `i∉J`, `Y` is supported on
the *difference* coordinates `i∈J`, and both are sums of the same family
`C_i √d_i`.  Since `T_{σ⁰}=X+Y` and `T_{σ¹}=X−Y`,

```
|X+Y| ≤ 2^{-Δ},   |X−Y| ≤ 2^{-Δ}   (Δ = min(E,Delta_2)/1),
```

so both `|X|` and `|Y|` are `≤ 2^{-Δ}`.  The M7.2 working hypothesis is
that those two simultaneous small sums, sharing the same coefficients and
radicand field, expose structure beyond the norm, precisely because they are
two independent radical forms of controlled support being simultaneously
tiny on the same `C`.

---

## 20.3 Data-to-collect, per residual candidate

For each `(C,d,σ⁰,σ¹)` meeting §20.1 record:

```
(E, Δ₂, Δ₂/E);
C ∈ Z^m,  gcd(C₁,…,C_m),  support{i : C_i ≠ 0};
the two forms X and Y explicitly: their signed coefficient vectors;
all remaining canonical embedding magnitudes;
J with |J|, r_J;
the radicand multiset of the generated subfield Q(√d_i : i ∈ J);
for each pair (a,b) of distinct radicands from J, the product d_a·d_b
(helps compute subfield discriminants);
the two-by-two difference vector (σ⁰_i − σ¹_i)/2 ∈ {0,±1}^m;
```

Then test, in order:

1. `gcd(C)=1` vs nontrivial — i.e. is the coefficient vector already
   primitive relative to the radicands?
2. coefficient support vs `m` — is it full or can a subset-only enemy
   explain the configuration?
3. the `X`/`Y` word length `|supp X| + |supp Y|`: the "residual regime"
   should manifest near the balanced boundary `|J|≈r−2`, so ask if the
   *number of non-zero entries on each side* correlates with `Δ₂/E`.
4. whether `Y`'s squarefree part is degenerate — `Y` is nonzero exactly
   when the `d_i:i∈J` are `F₂`-independent (Paper A); so `Y ≠ 0` already
   holds in the residual regime and yields a nonzero algebraic integer.
5. the induced quadratic subfield behaviour: each pair of `Y`-radicands
   defines a quadratic subfield; does its regulator / Pell behaviour
   explain the pair at random?  (Characterise and store; do NOT conclude
   from a single example.)
6. discriminants of the generated subfield `Q(√{d_i}:i∈J)`: compute
   exactly and store modulo a table so the same `d`/different `C` pattern
   can be compared.
7. the recurrence: if `C` yields residual data `(E_{r=4}...)`, does the
   analogous surplus survive after *reducing each d_i to its fundamental
   representative*?

---

## 20.4 Execution results

Summary (qualifying canonical pairs, `Delta_2 >= E/4`, `|J| >= r-2`, `r_J=|J|`):

```
r=3:  4/124 (B=2),  16/728 (B=4),  48/4912 (B=8)     rate 3.2% -> 0.98%
r=4:  36/624 (B=2), 156/6560 (B=4), 812/83520 (B=8)  rate 5.8% -> 0.97%
r=5:  192/3124 (B=2), 1760/59048 (B=4)               rate 6.2% -> 3.0%
r=6:  1184/15624 (B=2)                               rate 7.6%
r=7:  4896/78124 (B=2)                               rate 6.3%
```

Support distribution of `X=(T_0+T_1)/2` and `Y=(T_0-T_1)/2`:

```
(suppX, suppY) ~ (1, r-1) or (2, r-2)  for >99% of residuals
```

Key diagnostics:

> **RETRACTED (v1.1).**  The original v1.0 conclusion of this section
> ("near-antipodal geometry explains the residual regime", and the median
> `E_Y/E = 1.000`) was **wrong**, and rested on a buggy filter: the script
> zipped a genuine-case list against a differently-filtered list when
> computing medians, and its `X = 0` / "artifact" classification was based on
> floats.  Corrected figures are in `Docs/21`, §21.5.

**Corrected conclusion.**  The residual regime splits into two genuinely
different algebraic classes, separated by hypothesis (H) of `Docs/21`:

* **Genuine two-small** (`|X| ≥ |Y|/2`, i.e. hypothesis (H) holds): 1616 of
  the 9104 residuals, all with `supp(X) = 2`, all with `gcd(C) = 1`.  These
  are *not* artifacts: `X` is a genuine two-term form.
* **Near-antipodal**: the remainder.

Crucially `supp(X) = 2` in **every** genuine case (3344/3344 rows satisfying
(H), over all `B`).  So the correct statement is not "the phenomenon is
illusory" but "the agreement part is always tiny-support", which is exactly
the hypothesis under which `Docs/21` proves a strict improvement.

* Rate decays with `B` (3.2% -> 0.98% for r=3): finite-box coincidence.
* Deepest genuine case: `r=6`, `d=(2,3,5,7,11,13)`,
  `C=(−2,−1,−1,−1,−1,−1)`, `E = 4.947`, `E_X = 2.453`, `E_Y = 2.217`,
  `|J| = 4`, `supp(X) = 2`, `supp(Y) = 4`.

**Conclusion.**  The residual regime is *not* explained away as an artifact.
It is explained — provably and quantitatively — by the support structure of
the agreement form `X`, via `Docs/21`.  What remains open is the isolated
one-small-embedding case, not the two-small case.

## 20.5 Decide M7.2-phase stop

The entire purpose is bookkeeping, not a theorem.  The audit ends when one
of the following is unambiguous:

* a **stable invariant** surviving `B→2B`, a change of the prime set, and
  a change of `r` (3→5→7), which would be a serious candidate for Bucket 3;
* a **proof that all surviving residuals lie on the one-small boundary**
  after quotienting by the obvious transformations (`C→−C`, the
  global-flip class, monodromy on the sign group), which would close the
  residual theme without new algebra;
* a measured phenomenon **unique to the M4.1 full-rank threshold**, which
  would re-open the M4.1 scale with sharper constants.

If none of these, the files remain open, and Bucket 3 continues as
unadopted.

---

## 20.5 Relation to existing documents

* `Docs/12` (M4.1) proves that whenever `r_J < |J|` the corresponding
  mixed form `T_{σ⁰} − T_{σ¹}` already has a smaller effective field
  (exponent `2^{r_J}`), so the *full-rank residual regime is the only one
  in which M4.1 does not already hand us a reason for `Δ₂` to be large.
* `Docs/13` (M4-theory) proves the sign-difference rank of the column
  pattern is `≤ k−1`, which is exactly why the two-small geometry is
  invisible to a rank-compression, supporting the idea that the *pair*
  `(X,Y)`, not the rank of `D`, is the right finer invariant.