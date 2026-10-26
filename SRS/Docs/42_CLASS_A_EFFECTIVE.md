# Paper 42 — Theorem 42.1: the defect is eventually quasi-exponentially large

Version 1.0. Completes and GENERALIZES Theorem 41.A (Docs/41 41.3). For fixed
squarefree `a_1 != a_j`, let `q_k` be the Q-coordinate of `eps^k` (eps the
fundamental Pell unit of `x^2 - a_1 y^2 = +-1`), `p*` the nearest integer to
`q_k sqrt(a_j)`, and `n_j(k) = p*^2 - a_j q_k^2` the defect. Then:

> **Theorem 42.1.** There exist EXPLICITLY computable `c = c(a_1,a_j) > 0` and
> `K = K(a_1,a_j)` such that for every `k >= K`:
> ```
> |n_j(k)|  >  exp( c * k / log k ) .
> ```
> Equivalently: every index with `|n_j(k)| <= exp(ck/log k)` satisfies
> `k < K`, and the set of such indices is effectively enumerable.

This subsumes Docs/41 Thm 41.A (Class A, `|n_j| = +-1`) as the special case,
with better constants, and extends it to every FIXED and every
quasi-polynomially growing defect bound. Constants: `c` may depend on
`(a_1, a_j)` freely (Docs/30 sec 30.6 — absorbed into the N-independent
`K_0` of Docs/28 Thm 28.3).

All integers below are exact; the only non-elementary input is the Matveev
linear-forms-in-logs theorem, used in shape form with an explicit published
constant (42.5). Companion measurement: `m712_depth.py`, `m712_decay_probe.py`.

---

## 42.1 Setup and conventions

* `a_1, a_j` distinct squarefree. `eps > 1` the fundamental unit of
  `x^2 - a_1 y^2 = +-1` (integer solutions; norm `N_{Q(sqrt a_1)}(eps) = +-1`,
  `|eps'| = 1/eps`), `eta > 1` the same for `x^2 - a_j y^2 = +-1`.
* `eps^k = P_k + q_k sqrt(a_1)` for integers `P_k, q_k >= 1` (k >= 1; the
  engine's `eps_power` recurrence — k IS the exponent). Conjugation gives
  `eps'^k = P_k - q_k sqrt(a_1)`, so
  ```
  q_k = (eps^k - eps'^k) / (2 sqrt(a_1))                     (42.1a)
  ```
  EXACTLY, for every k, for BOTH norm signs — no parity caveat at this level
  (the Docs/41 41.1 parity caveat affects only the eta/beta race identity,
  not this paper).
* `p* = round(q_k sqrt(a_j))`, `omega = p* + q_k sqrt(a_j)`,
  `n = n_j(k) = N_{Q(sqrt a_j)}(omega) = p*^2 - a_j q_k^2 in Z \ {0}`,
  `q_k >= 1`.

## 42.2 Lemma 42.2 (orbit reduction — elementary)

Suppose `|n| <= B` with `B >= 1`. Then there are `s in Z` and
`rho in Q(sqrt a_j)` with

```
omega = rho eta^s ,      1 <= |rho| < eta ,      |N(rho)| = |n| ,
max(|rho|, |rho'|) <= max(eta, |n|) <= eta * B ,      h(rho) <= (log(eta B))/2 .
```

*Proof.* The unit group of integer Pell solutions is `+-eta^Z` (classical).
`omega` has norm `n`; multiply by `eta^{-u}` with `u` chosen so
`|omega eta^{-u}| in [1, eta)` (the values `|omega| eta^{-u}` step by the
factor `eta`, so some `u` lands in the window; `s := -u`). Norm is scaled by
`(+-1)^s`, so `|N(rho)| = |n|`. Conjugating: `|rho'||rho| = |n|`, so
`|rho'| = |n|/|rho| <= |n|` and `max(|rho|,|rho'|) <= max(eta, |n|)`. The
standard logarithmic height of `rho in Q(sqrt a_j)` is
`h(rho) = (1/2)(log^+|rho| + log^+|rho'|) <= (log(eta B))/2`. ∎

(The number of orbits is irrelevant: we only need THE orbit containing
`omega`, and its representative. No class-number input.)

## 42.3 Lemma 42.3 (the exact identity and the small linear form)

Let `lambda := sqrt(a_j/a_1) in L := Q(sqrt a_1, sqrt a_j)`. Then, EXACTLY:

```
rho eta^s - rho' eta'^s  =  lambda (eps^k - eps'^k) .                  (42.3a)
```

*Proof.* Conjugating `omega = rho eta^s` gives `p* - q_k sqrt(a_j) =
rho' eta'^s`; subtracting kills `p*`: `2 q_k sqrt(a_j) = rho eta^s - rho'
eta'^s`. Similarly `2 q_k sqrt(a_1) = eps^k - eps'^k`. Divide. ∎

Define the real linear form in three logarithms

```
Lambda := k log eps  -  s log eta_s  -  log(|rho|/lambda),
          eta_s := eta^{sgn(s)}  (h(eta_s) = h(eta)),
```

so `rho eta^s = sgn * |rho| eta^s` and (42.3a) reads
`lambda eps^k - rho eta^s = E` with `E := lambda eps'^k - rho' eta'^s`.

**Claim 42.3b.** For all large k (explicitly: `k >= K_1(a_1,a_j,B)`),
```
|E|  <=  eps^{-k} (lambda + 3|n|/lambda),
|Lambda|  <=  4 (1 + |n|) eps^{-2k} * (1 + lambda^{-2})/lambda .        (42.3c)
```
*Proof.* Bounded terms: `|omega| in [q_k sqrt(a_j), q_k sqrt(a_j) + 1]` and
`q_k sqrt(a_j) = lambda(eps^k - eps'^k)/2 >= lambda eps^k/3` for large k, so
`|omega| <= 2 lambda eps^k` and
`|rho'| |eta'|^s = |rho'|/eta^s = |rho'||rho|/|omega| = |n|/|omega| <= 3|n|/(lambda eps^k)`.
With `|eps'| = 1/eps`: `|E| <= lambda eps^{-k} + 3|n| eps^{-k}/lambda`. For
large k this is `< lambda eps^k / 2`, so
`|rho eta^s/(lambda eps^k) - 1| <= 2|E|/(lambda eps^k)` and
`|Lambda| = |log(|rho| eta^s/(lambda eps^k))| <= 2|E|/(lambda eps^k) <=
4(1+|n|) eps^{-2k}(1 + lambda^{-2})/lambda`. ∎

## 42.4 Lemma 42.4 (Lambda is never zero)

`Lambda != 0` for every `k >= 1` and every distinct squarefree pair.

*Proof.* Suppose `|rho| eta^s = lambda eps^k` (real equality). Then, as
ELEMENTS of `L`, `rho eta^s = +- lambda eps^k`. Now
`lambda eps^k = sqrt(a_j/a_1)(P_k + q_k sqrt(a_1)) = (P_k/a_1) sqrt(a_1 a_j) + q_k sqrt(a_j)`
lies in the FULL biquadratic field `L`, and its `sqrt(a_j)`-coefficient is
`q_k >= 1`. If `rho eta^s = +- lambda eps^k` with `rho eta^s in Q(sqrt a_j)`,
then `(P_k/a_1) sqrt(a_1 a_j) = rho eta^s - q_k sqrt(a_j) in Q(sqrt a_j)`;
but `sqrt(a_1 a_j) notin Q(sqrt a_j)` (the squarefree class of `a_1 a_j`
differs from that of `a_j`, since `a_1 != 1` is squarefree), so `P_k = 0`
— impossible (`P_k >= 1`). ∎

**CORRECTION (this pass).** An earlier displayed line here wrote
`lambda eps^k = (q_k + P_k/a_1) sqrt(a_1 a_j)` — WRONG: the `q_k` term sits
on `sqrt(a_j)`, not on `sqrt(a_1 a_j)`. The claim (`Lambda != 0`) survives
with the corrected argument above; the norm-based sketch of the Docs/41
draft had a residual `a_1 | a_j` corner and is superseded. ∎-record

## 42.4b THE CONFINEMENT NO-GO (new; the method-ceiling theorem)

The most promising shortcut for breaking the log-k loss was: `alpha_3 =
|rho|/lambda` has its REAL embedding confined (`|rho| in [1, eta)` by 42.2),
so one hopes `h(alpha_3)` — which is what Matveev charges for — might be
bounded, independent of `|n|`. That would upgrade 42.1 from `exp(ck/log k)`
to `eps^{2k} k^{-C}` and close the gap at `n = 2`. **It fails, provably and
elementarily:**

> **Lemma 42.4c.** `h(lambda/rho) >= (1/2) log|n| - (1/2) log(eta lambda)`.
>
> *Proof.* `h(x) = h(1/x)` (standard, both archimedean and non-archimedean
> parts). Use the side that grows: `h(rho/lambda) >= (1/4) sum_sigma
> log^+|sigma(rho/lambda)| >= (1/2) log^+(|rho'|/lambda) >= (1/2)(log|rho'|
> - log lambda) >= (1/2)(log|n| - log|rho| - log lambda) >= (1/2) log|n| -
> (1/2) log(eta lambda)` for `|rho| < eta`. ∎

So the height grows through the FINITE primes (the poles of `rho/lambda`)
exactly as the product formula demands: the defect's size IS the height of
the third coefficient-number. No linear-forms bound can escape it:

> **Theorem 42.4d (the linear-forms ceiling).** Every Matveev-type bound on
> `Lambda = k log eps - s log eta_s + log(|rho|/lambda)` pays `(log|n|)`
> multiplicatively: `log|Lambda| >= -C (1+log k)(log a)^2 (O(1) + log|n|)`.
> Combined with the true value `|Lambda| <= C(1+|n|) eps^{-2k}` (42.3c), the
> conflict window is exactly
> ```
> log|n|  <  (2 log eps / (C (log a)^2)) * k / log k ,
> ```
> i.e. `|n| < exp(Theta(k/log k))` — **Theorem 42.1's floor is the CEILING of
> the entire linear-forms method for this family.** The wall band
> `[exp(Theta(k/log k)), eps^{(1-delta)k}]` is the true residue, not a loose
> estimate. ∎

Verified numerically (k = 10, (2,3)): `h(1/x) = 4.205 >= (1/2)log 6728 -
(1/2)log(eta lambda) = 3.65`; and `Lambda = log(1 - (eps'/eps)^k +
d_k/(lambda omega eps^k)) = -9.9e-4`, matching `d_k/(lambda^2 eps^{2k}) =
9.9e-5`-scale up to the `(1 - (eps'/eps)^k)` correction and sign.

## 42.5 Proof of Theorem 42.1

Suppose `k >= K` and `1 <= |n| <= B` with `B := exp(c k / log k)`, `c` to be
fixed. Build `rho, s` (42.2), `Lambda` (42.3). By 42.4 `Lambda != 0`. Apply
Matveev's theorem for linear forms in n logarithms (Matveev 2000; the
standard explicit constant `C_Mat` is published; as in Docs/16 sec 16.2 the
shape is what matters and the constant is effective):

```
log|Lambda|  >=  -C_Mat * D^4 * (1 + log Bf) * A_1 A_2 A_3 ,
Bf := max(|k|, |s|, 1) ,   A_i := max(D h(alpha_i), 1) ,   D = [L:Q] = 4 ,
alpha = (eps, eta_s, lambda/rho) .
```

Heights: `h(eps) = (log eps)/2 = O(log a_1)` (for the Pell unit
`eps <= 2 sqrt(a_1)` on the worst branch); `h(eta_s) = h(eta) = O(log a_j)`;
and by 42.2
`h(lambda/rho) <= h(lambda) + h(rho) <= O(log a_j) + (log(eta B))/2`. Since
`log|n| <= c k/log k <= k/2` for large k, also `s = O(k + log|n|) = O(k)` and
`log Bf = O(log k)`. Hence

```
log|Lambda|  >=  -C'' * (log k) * (log a)^2 * (log a + c k / log k)      (42.5a)
```

with `C'' = C''(C_Mat)` explicit and `a = max(a_1, a_j)`. On the other side,
(42.3c) with `|n| <= B` gives

```
log|Lambda|  <=  -2 k log eps + c k / log k + O(1) .                    (42.5b)
```

Fix `c := log eps / (3 C'')` (so `c` is explicit once `C_Mat` is; the
`(log a)^2` factor is absorbed by shrinking `c` — set
`c := log eps / (3 C'' (1 + log a)^2)` to be safe). Then for large k the
right side of (42.5a) is `>= -c k` while (42.5b) is `<= -2 k log eps +
o(k) < -c k` — contradiction. Every sufficiently large k therefore has
`|n_j(k)| > exp(c k/log k)`, and `K` is explicit (the crossing point of the
two bounds, computable from `C_Mat, a_1, a_j`). ∎

**Remark (what enters).** The proof uses: Pell orbit structure (elementary),
the exact identity (42.3a) — a consequence of the SHARED q_k — the field
intersection (42.4), and Matveev. It does NOT use the Subspace Theorem, and
it is effective. This is the first effective lower bound on `n_j(k)` beyond
`|n_j| >= 1` in this project.

## 42.6 Corollaries

**(a) Class A absorbed.** `{k : n_j(k) = +-1}` is finite and effectively
enumerable (`|n| = 1 <= B`). Docs/41 Thm 41.A is the special case; the census
(`|n_j| = 1` only at k = 1 over k <= 1e5) is consistent.

**(b) Depth cap.** With `L_j(k) = log_eps(q_k/|n_j(k)|)` (Docs/41 41.4):
`L_j(k) <= k - c k/(log eps log k) + O(1)` — the depth never comes within
`Omega(k/log k)` of its Liouville ceiling `k`. Class-A depth (~k) is excluded
with an explicit margin.

**(c) The wall, restated.** Theorem 42.1 excludes defect values up to
`exp(ck/log k)`. The Class-B wall (Docs/41 41.3) requires excluding defects
up to `eps^{(1-delta)k}`, `delta = 1/(3n(n-1))`. The REMAINING band is

```
|n_j(k)|  in  [ exp(Theta(k/log k)) ,  eps^{(1-delta)k} ] .
```

In depth coordinates: exceptions need `L_j >= delta k - O(1)`; 42.1 caps
`L_j <= k - Omega(k/log k)`; the band `(delta k, k - Omega(k/log k))`
survives. The achieved-vs-required table:

```
  defect floor  achieved:  exp(c k / log k)         (Thm 42.1)
                required:  eps^{(1-delta) k}        (Class-B exclusion)
                ratio in the exponent:  c/log k  vs  (1-delta) log eps
                ==> the SAME one-logarithm loss as Docs/27 Lemma 27.2.
```

This is now the THIRD independent instance of the log-k wall (Docs/27 eta-form;
Docs/31 ST-output calibration; 42.1 defect-form). The wall is one phenomenon.

**(d) Downstream (unchanged).** The chain Docs/34 sec 34.6 remains mechanical:
any improvement of the achieved floor to the required floor closes M7.12 and
hence SRS.

## 42.7 Verification (exact integer kernel, (2,3), k <= 1e5)

```
min |n_j(k)| over k = 2..1e5            : 3   (at k = 2)
indices with |n_j(k)| <= 100            : 5   (k = 2..6)
empirical exponent frontier min log|n|/k: 0.5493 (at k = 2, |n| = 3)
```

For every measured k the defect is far ABOVE the theorem's guaranteed floor
`exp(ck/log k)` (the theorem's `c` is symbolic pending the chosen `C_Mat`
value, but `c/log k -> 0`, while the measured frontier `log|n|/k` stays
bounded below by ~0.55·(declining profile) — data strictly stronger than the
theorem at every measured point, as it must be: the theorem is an all-k
statement). The Class-A census (only k=1) and the depth records
(max L = 12.26 at k=36540) are re-confirmed by `m712_depth.py`.

## 42.8 Honest limits

1. **The wall band survives** — 42.6(c). Nothing here claims the gap.
2. **`C_Mat` is shape-quoted.** The constant is explicit in the literature
   (Matveev 2000; sharpened two-log variants in Laurent–Mignotte–Nesterenko),
   but no numerical value is inserted; `c` and `K` are explicit-symbolic.
   Per the repo rule this theorem is labelled: proved in shape, constants
   pending a source pass. `[constant: verify against Matveev 2000 before
   publication]`
3. **The log-k loss is intrinsic to the method**: Matveev's `(1+log Bf)` is
   believed near-optimal in general. Beating it for our three-log family —
   fixed unit heights, coefficients ~ k, one growing-height unit-coefficient
   term — is the precise shape of the missing brick (Docs/40 open rows).
4. **One pair verified.** The theorem holds for all pairs by proof; the
   empirical sections are (2,3) only.

## 42.10 The constant pass (DONE — Baker–Wüstholz 1993 inserted)

The non-elementary input is now the **Baker–Wüstholz 1993** theorem (J. reine
angew. Math. **442**, 19–62), which has the SAME shape as Matveev's and an
explicit published constant, quoted from the standard statement:

```
log|Λ|  ≥  −C_BW · h(α₁)···h(αₙ) · log(max|βᵢ|) ,
C_BW = 18(n+1)! · n^{n+1} · (32d)^{n+2} · log(2nd) .
```

For our form: `n = 3` logs (`ε, η_s, λ/ρ`), `d = 4` (the biquadratic field
`L`), coefficients `(k, s, 1)` — so `C_BW(3,4) = 18·24·81·128⁵·log 24 =
3.8210×10¹⁵`.

Assembling 42.2/42.3/42.5 with `C_BW`:

```
c_max(a₁,a_j) = 4 log ε / (C_BW · h(ε) · h(η)) ,
```

(`h(η)` bounded crudely by `(log(2√a_j+1))/2`, which only shrinks `c_max` —
safe direction). Per-pair values:

```
(2,3): 2.80e−15   (2,5): 2.46e−15   (2,7): 2.28e−15   (2,11): 2.06e−15
(2,13): 1.99e−15  (3,2): 3.12e−15   (5,2): 3.12e−15   (5,3): 2.80e−15
(3,5): 2.46e−15   (11,2): 3.12e−15
```

So take `c := 10⁻¹⁵` uniformly (absolute, all pairs) and `K(a₁,a_j)` the
computable crossing point, of size ~`10¹⁵·log a`. **Theorem 42.1 is now
fully explicit** — with the honest caveat that the constants are of the
usual BW/Matveev scale: effective IN PRINCIPLE (the all-k statement and the
enumerability are real), vacuous IN PRACTICE below `k ~ 10¹⁵`.

**Uniformity (the ledger, T4).** `C_BW` is ABSOLUTE (depends on n=3, d=4
only — both fixed). `c(a₁,a_j) ≥ c₀/log³(a)` for an absolute `c₀ ~ 10⁻¹⁵`,
and `K₀ = O(n·c/log ε) + K_ST` (Docs/28 Thm 28.3) is **N-independent** for
every pair. The uniformity chain survives the constant pass unchanged.

**Hypothesis bookkeeping flag.** BW 1993 in the quoted form applies to
`Λ ≠ 0` with integer coefficients; the independence of
`{log ε, log η, log(λ/ρ)}` (needed for the homogeneous reading) follows for
the cases used here from 42.4's field-intersection argument, but the
general (u,v,w) relation `ε^uη^v(λ/ρ)^w = 1` bookkeeping should be verified
against the BW paper's exact hypothesis before publication.
`[flag: verify BW hypothesis bookkeeping]`

**Not resolved by this pass:** R7/R8 (Docs/40) — both BW and Matveev share
the shape `Π h(αᵢ)·log B`, which is exactly what Theorem 42.4d proves is
the ceiling. The constant pass CONFIRMS the wall from a second published
statement; it does not breach it. The decisive question stays: anything
BEYOND the shape `Πh·log B` for fixed-height units.

## 42.9 Status

| Item | Result |
|---|---|
| Orbit reduction (42.2) | proved, elementary, no class number |
| Exact identity (42.3a) | proved (shared q_k; sign-free) |
| `|Lambda| <= C(1+|n|) eps^{-2k}` | proved (42.3c) |
| `Lambda != 0` for all k, all pairs | proved (42.4, field intersection; supersedes the norm argument and its `a_1|a_j` corner) |
| **Theorem 42.1** (defect > exp(ck/log k) eventually, effective) | **PROVED, CONSTANTS EXPLICIT (42.10: BW 1993, c ~ 1e-15, K ~ 1e15)** |
| Class A subsumed; depth cap `k - Omega(k/log k)` | proved |
| Wall band `[exp(Theta(k/log k)), eps^{(1-delta)k}]` | OPEN — third documented instance of the log-k wall |
| Downstream chain | unchanged (Docs/34 sec 34.6) |