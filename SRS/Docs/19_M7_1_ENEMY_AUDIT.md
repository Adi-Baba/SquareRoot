# M7.1 — Enemy audit, one record per closed mechanism

Version 1.0.  The aim is not a smallest SRS value; it is an instance that
maximises the information a given *closed* route consumes while that same
route fails to separate it.  Every row uses the uniform record

```
input family  ->  S  ->  smallness  ->  visible-to-route  ->  why insufficient
```

The seven mechanisms are indexed by the closed table in `TODO.md` M7.

| Closed mechanism (Docs ref) | Enemy construction | Smallness / E, N | P* or visible data | Why the route fails to certify |
|---|---|---|---|---|
| Full norm (Paper B)        | One exceptional sign-flip: take `C` from the r=3 consecutive-gap enemy with target embedding tiny, all other signs generic | `E/N≈0.3–0.71`, `E ≤ 6N+O(1)` already effective | `P_* = Π_{σ≠σ0}|S_σ|` with `log₂P_*/N ≈ Θ(2^r/r)` (Paper 11) | `P_*` sits at full rank-`2^r` scale, so `|S| ≥ 1/(27 d₃^{3/2}M^{9/2})` cannot be improved by pairing; exponent `2^{r-1}−1` is structural. |
| Relative norm tower (Docs/10) | Same family under successive quadratic norms `N_{K_i/K_{i-1}}` | Telescoping exposure of the same count of embeddings | `N_{K_i/K_{i-1}}(L) = L·τ_iL`, product over telescope = `N_{K/Q}(L)` | The telescope equals the single-step full norm; no extra restriction at any level. |
| Sign-pattern difference rank descent (Docs/13) | `C` concentrated by construction in a block sharing the *same* change pattern on the flipped coordinates | The pair-wise difference matrix D has the minimal possible nontrivial-pattern rank | column constant on blocks `V_c`; `B_c ≍ M^{1/2}` not small | Only the combination `PᵀB` through the column-pattern vectors is required to be small; the whole `B_c` is not forced small.  `D` of rank `q<m` leaves `m−q` free coordinates. |
| One-small pigeonhole (Docs/11) | Consecutive-gap gaps on the prime radicand family | `E/N` at the `(r−1)/r` ceiling with `N≈r log B` | the remaining embeddings are generic; no two small | It already saturates the lower side exactly: `|S| ≍ (Σ√p_i)·B^{-(r−1)}`.  Identifying "what a uniform effective theorem must reject" is: the *thin* simultaneous regime below, not this norm-tight one. |
| Two-small lemma (Docs/12) | pairs `σ,τ` with both `|T_σ|,|T_τ|≤ ε ≪ M^{3/2}` sought over `{d_i}` prime, small `C` | measure of any computed in `M7.1` | we found *all* such pairs (r=3 196/196, r=4 1564, r=5 400 samples to B=4 thr=1e-2) have difference support of full rank `r_J=#`diffs — in other words, the difference is supported on a single independent sub-basis of maximum size | This is the direct witness: there is **no** pair in the search whose sign-difference set is rank-deficient; the M4.1 threshold is never passed, so the induced "lower-rank form" does not simplify the estimate. |
| Generator/coefficient descent | sparse or cancellation-heavy `C` on a single radical subset | support rank evaluated directly by F₂ row-reduction | In every tested example the one-small embedding is the all-principal-roots one and is minimised by a pigeonhole-mixed signed sum; the coefficient vector has a cancellation pattern across the *full* prime set, not on a low-rank support | Any descent would have to impose a relation whose support rank is strictly less than full — which M4.1 already denies unless ε is at least the full-rank norm threshold. |
| Generic Baker / unit (Docs/16) | rank-3 `H(C)=±1` branch, form `ρ = d₃u²` | `d₃²F(C) ∈ Z`, dangerous branch `H=±1` | unit heights via Zassenhaus; regulator is fixed by `K` | `𝔇 = C₃·4⁵·(1+log B)·Πh(ε_i)` in the three-log estimate has `𝔇 ≤ 2^{O(N)}` only; cannot be combined with a single-coordinate equation and cannot supply a polynomial upper bound on `E_L`. |

---

## 19.1 The decisive two-small search result (corrected: antipodal flip removed)

Flask interpretation:

> In the tested prime-radicand/small-coefficient families, after removing the
> antipodal duplication, every observed pair of simultaneously small
> embeddings has an `F₂`-independent difference support (`r_J = |J|`).

Antipodal cleanup matters: without pinning one representative per
`(σ, −σ)` class via e.g. `ε₁ = +1`, the second-smallest magnitude in the
sorted embedding list is often the trivial global-flip copy of the first,
so such a `Δ₂ = −log₂|S_(2)|` does not measure independent simultaneous
smallness.  All numbers below are from the de-redundant profile (one
representative per antipodal class).

Refreshed table (qualifying canonical pairs with `Δ₂ ≥ E/4`):

```
r=3, B=4:  16/728  qualifying;   diff/rank histogram {(1,1):12, (2,2):4}
r=4, B=3: 244/2400 qualifying;  histogram {(1,1):152, (2,2):52, (3,3):40}
r=5, B=2: 900/3125? (3124);     histogram {(1,1):548, (2,2):160, (3,3):80, (4,4):112}
```

There are three observed regimes: (A) rank-1/Pell local  (`|J|=r_J=1`,
dominates); (B) short-but-independent multi-radical supports
(`|J|=r_J ∈ {2,3}`); (C) large/full residual support (`|J|≈r−1`, `r_J=|J|`,
e.g. the r=5 `(4,4)` entries).

Importantly:

* There is **no** observed rank-deficient difference support.  After the
  antipodal de-redundancy, M4.1 behaves exactly as the theory predicts:
  whenever the second small embedding is genuinely distinct, its difference
  support is independent in these samples.
* But this does **not** mean M4.1 explains the simultaneous smallness.
  It has only explained why the obvious *lower-rank escape route* does not
  appear.  Regime C's full-rank configuration is **not** disposed of by M4.1.

The high-frequency rate of category-C-type events (e.g. `112/900` of
qualifying r=5 pairs) shows that this residual is not a negligible exceptional
set; it recurs at nontrivial proportion.  No measure-theoretic "measure zero"
language is attached to it: the enumeration only establishes nontrivial
frequency at small `C`.

---

## 19.2 Why this reframes Bucket 1/2/3 (corrected)

Three buckets from `Docs/18`:

* Bucket 1 (algebraic beyond full norm): **narrowed, not refuted.**  The
  specific candidate mechanism
  `two simultaneous smalls ⇒ rank-deficient difference support ⇒ smaller-field
  norm` is empirically absent and already bounded below by M4.1.  But
  special-field information can still live in the full-rank regime.
* Bucket 2 (simultaneous approximation beyond one form): **open; elementary
  sign-pattern realisation exhausted.**  The sign-difference matrix does not
  compress the field degree (`Docs/13`), but the full-rank residual
  configuration is exactly simultaneous closeness of many independent
  radicals to rationals and to a single `C`.  That quantitative geometry is
  not produced by `M D` rank alone.
* Bucket 3 (special-field mechanism): **open.**  Now becomes a natural
  candidate *because* M7.1 eliminated one obvious low-rank route, not
  because the computation proved it to be the only possibility.

Strongest statement this record currently supports:

```
low-rank-support algebraic descent:        not observed; M4.1 controls it.
Full-rank simultaneous-small configs:      definitely observed at nontrivial rate.
Extra arithmetic structure:                unknown; see M7.2.
```

---

## 19.3 Limitations of the M7.1 compute job

The search was over the first few primes with small `C`, and interprets
every pair `T_{σ⁰},T_σ` of the same sample.  It is evidence, not proof;
its value is that, to leading order, it shows *what a uniform effective
argument cannot rely on*.  The only way to make such statements theorems is
to enhance the analytic lower bound; where M7.1 stops being a search and
starts being a proof is the statement M4.1 already proved.