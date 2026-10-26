# Paper 41 — the class decomposition and the softness-depth function

Version 1.0. Produced by the L13 attack plan (TODO). Three things happen here:
the EHS/KMN exception set is DECOMPOSED into two classes with different
mechanisms (41.3); the deeper class (A) is effectively bounded by a classical
Baker argument that Docs/23-39 missed (41.4); and the entire remaining gap is
reformulated as a tail bound on a NEW computable function, the softness depth,
whose measured behaviour is metric-like with the exact predicted rate (41.5).

Companion module: `m712_depth.py` (MEASUREMENT ONLY). All integers exact.

---

## 41.1 The race identity (proved)

Let `q_k` be the Q-coordinate of `eps^k` (eps = fundamental Pell unit of
`x^2 - a_1 y^2 = +-1`), `p*` the nearest integer to `q_k sqrt(a_j)`, and

```
n_j(k) := p*^2 - a_j q_k^2          (nonzero integer, the defect).
```

From `q_k sqrt(a_j) = beta_j(eps^k - eps^{-k})` (Docs/34 corrected; sign
caveat below), `p* + q_k sqrt(a_j) = 2 q_k sqrt(a_j) (1 + O(q^{-2}))`, and
`delta_k := q_k sqrt(a_j) - p* = -n_j/(p* + q_k sqrt(a_j))`:

```
eta_k = ||beta_j eps^k||  =  | beta_j eps^{-k} - n_j(k)/(2 q_k sqrt(a_j)) |        (41.1)
        (exact whenever |delta_k -/+ beta_j eps^{-k}| <= 1/2; asymptotically exact)

2 q_k sqrt(a_j) * beta_j eps^{-k} = (a_j/(2 a_1)) (1 - eps^{-2k})                  (41.2, EXACT)
```

So `eta_k` is a RACE between two tiny terms: the universal correction
`beta_j eps^{-k}` and the defect term `n_j/(2 q_k sqrt(a_j))`. Verified
numerically to 1e-17 (even k) / float precision (odd k) at k = 20..41.

**Parity caveat (new; N(eps) = -1, i.e. a_1 = 2).** `eps^{-k} = (-1)^k (P_k -
Q_k sqrt(a_1))`, so the clean identity `q_k sqrt(a_j) = beta_j(eps^k -
eps^{-k})` holds for EVEN k only; odd k needs `beta_j(eps^k + eps^{-k})`.
Magnitudes and every O(eps^{-k}) statement are unaffected (the beta-term
enters with sign `(-1)^k`). For units of norm +1 (the generic fundamental
Pell unit) the identity holds for all k. This retroactively explains the
Docs/34 section 34.2 "index shift" pathology, and the m=9 "verification"
there compared two values differing at digit 13.

**Unit caveat (new).** The Pell connection requires eps = the fundamental
unit OF THE PELL EQUATION (norm +-1); `1 + sqrt(a_1)` is that unit only for
a_1 = 2. Docs/30 section 30.5's Pisot dichotomy is stated for `1 + sqrt(a_1)`.
Theorem 30.3 survives regardless (Prop 30.2's archimedean failure of
pseudo-Pisot is eps-independent). Cleanup flagged; the a_1 = 2 calibration
chain is unaffected.

## 41.2 The defect-window theorem (proved)

With `theta_0 = eps^{-1/(3n(n-1))}` (Docs/29), (41.1)-(41.2) give, for all
large k:

```
eta_k < theta_0^k   <=>   | n_j(k) - (a_j/(2a_1))(1 - eps^{-2k}) |
                            <~  2 q_k sqrt(a_j) theta_0^k
                          ~  C(a_1,a_j) eps^{(1 - 1/(3n(n-1))) k} .       (41.3)
```

The target `(a_j/(2a_1))(1 - eps^{-2k})` converges EXponentially to the fixed
rational `a_j/(2a_1) = O(1)`. Hence, for large k, (41.3) forces
`|n_j(k)| <~ C eps^{(1-delta)k}` with `delta = 1/(3n(n-1))`:

> **Exception of the reduced problem <=> the defect is mid-soft:
> `|n_j(k)| <~ eps^{(1-delta)k}`.** The generic defect is `Theta(q_k) =
> Theta(eps^k)`; the Pell floor is `|n_j| = 1`. Everything between
> `eps^{(1-delta)k}` and the generic scale is forbidden by SRS.

This is Docs/35 Thm 35.4's 9n(n-1) gap and Docs/36 Thm 36.2's no-go territory,
now restated as a statement about a single integer sequence `n_j(k)`.

## 41.3 The class decomposition

```
Class A:  n_j(k) = +-1   (q_k is simultaneously a Pell denominator of
                          sqrt(a_1) and of sqrt(a_j)).
          Here both race terms sit at eps^{-k} scale; cancellation corners
          exist (a_j = 2a_1 with n_j = -1 drives eta to eps^{-2k} scale).

Class B:  2 <= |n_j(k)| <~ eps^{(1-delta)k}   ("mid-soft" defects).
```

**Theorem 41.A (Class A is effectively bounded; proof sketch).** At a
Class-A index, `p_k + q_k sqrt(a_1) = zeta_1 eps^t` and
`p* + q_k sqrt(a_j) = zeta_2 eta^s` with eps, eta the fixed fundamental Pell
units (norm +-1) and t ~ s ~ k. Equality of the common Q-coordinate gives

```
zeta_1 eps^t - lambda zeta_2 eta^s = zeta_1' eps'^t - lambda zeta_2' eta'^s ,
lambda = sqrt(a_1/a_j),      |RHS| <= C eps'^t ,
```

so the three-log form `Lambda = t log eps - s log eta - log(zeta_1/(lambda
zeta_2))` satisfies `|Lambda| <= C rho^t`, `rho = eps'/eps < 1` —
EXPONENTIALLY small. All three algebraic numbers are FIXED (heights O(log a));
only the coefficients t, s grow, entering Matveev's bound as log B:

```
|Lambda| >= exp( - C' (log k) log eps log eta log(a_1 a_j) ) = exp(-C'' polylog(a) log k).
```

Since `e^{-Theta(k)} << exp(-C'' polylog(a) log k)` for large k, contradiction;
the threshold is explicit and tiny (`k* = O(polylog(a) log k*)`). This kills
the deepest possible dips with an explicit classical bound. NOT in Docs/23-39:
Docs/37 Track A searched LRS-count literature and missed the simultaneous-Pell
route. TO COMPLETE: write out Matveev with explicit constants; handle the
`a_j = 2a_1` cancellation corner (constant-level only, same kill).

**Empirical census (this pass, k <= 1e5, exact kernel): Class A = {k = 1}
ONLY.** (q_1 = 1: 4 - 3 = 1.) Consistent with Theorem 41.A.

**Class B is the entire remaining wall.** No mechanism is known; every
conjugate-product attempt is excluded by Docs/36 Thm 36.2; the empirical
margin is enormous (below).

## 41.4 The softness-depth function (the new attack object)

```
nu_j(k) := |n_j(k)| / q_k                    (relative defect)
L_j(k)  := log_eps( q_k / |n_j(k)| )         (SOFTNESS DEPTH, base eps)
```

Properties (all proved or exact-verified):

1. `eta_k ~ nu_j(k)/(2 sqrt(a_j))` (from 41.1). So M7.12 (`eta_k >= k^{-C}`)
   is EXACTLY the statement `L_j(k) <= C log_eps k + O(1)`.
2. Liouville cap: `L_j(k) <= log_eps q_k = k + O(1)` (n_j >= 1). Equality is
   approached only on Class A (depth ~ k).
3. The gap (Class-B exclusion) is the depth TAIL bound:
   `L_j(k) >= delta k - O(1)` forbidden for large k, `delta = 1/(3n(n-1))`.
4. For every threshold L, the candidate set `{k : L_j(k) >= L}` is EXACTLY
   computable from the integer kernel — a one-parameter family of shrinking-
   target problems interpolating between the full exception set (L = delta k)
   and the deep records (L ~ 12).

**Measured (k <= 1e5, (2,3)):**

```
max depth  L = 12.26 at k = 36540   (the wall there: delta*k = 6090)
depth tail counts vs metric prediction 2 K eps^{-L}:
  L :  2      4      6      8     10     12
  obs: 9911   1682   300    57    10     1
  obs/pred: 0.29 0.29 0.30 0.33 0.34 0.20   <- STABLE ratio, rate = eps^{-L}
```

Two readings, both important:

* The depth process decays at EXACTLY the metric rate `eps^{-L}` in L, with a
  stable constant ~0.3 (vs the uniform value 2 — the folded-ray orbit is not
  equidistributed; the CONSTANT is biased, the RATE is metric). This is the
  strongest evidence yet that the mechanism is the Docs/39(e) shrinking-target
  one, with target radius ~eps^{-L}.
* The depth is FLAT: bounded by 12.26 over four decades of k, while the wall
  sits at delta*k = 16667 at k = 1e5. The gap between observed (O(1)) and
  required (delta*k) is a factor ~1300 in the exponent and growing.

## 41.5 What the function buys, and the honest limit

**Buys:**
1. UNIFICATION: Class A (depth ~ k, now effectively bounded) and Class B (the
   wall) are two regimes of ONE function; the gap is one tail bound.
2. The stronger conjecture the data suggests: `L_j(k) = O(log k)`, maybe O(1).
   Any of these is far stronger than the required `delta k` — so the target
   theorem has enormous slack. The required statement is WEAK; that is new
   information about where to look.
3. The exact interface for the effective Borel-Cantelli step: metric theory
   gives, for the folded-ray orbit (Docs/39(e)) with target radius
   `psi(k) = eps^{-delta k}`, a.s. finiteness because `sum psi(k) < inf`.
   The measured rate `eps^{-L}` says the orbit behaves metrically per-k; the
   missing theorem is the EFFECTIVE Borel-Cantelli for this one algebraic
   orbit — and `L_j` is the quantity it must bound.

**Does not buy:** a proof. Nothing here bounds the tail; Docs/36 Thm 36.2
still excludes every conjugate-product mechanism, and no effective metric
theorem for a fixed algebraic orbit is known (Docs/40 first pass: none found).

## 41.6 Status

| Item | Result |
|---|---|
| Race identity (41.1)-(41.2) | proved; verified 1e-17; parity caveat recorded |
| Defect-window theorem (41.3) | proved |
| Class A effectively bounded | Thm 41.A sketch; completion = explicit Matveev constants |
| Class A census k <= 1e5 | {1} only — consistent |
| Softness-depth function | defined; properties proved |
| Depth flatness + metric rate | measured (12.26 max; rate eps^{-L} stable) |
| Class-B exclusion (the wall) | OPEN — now one sentence about one integer sequence |
| eps/unit consistency cleanup | flagged (41.1 caveats); needed for uniformity in a_1 |

**Next actions.** (i) Complete Theorem 41.A with explicit constants
(Docs/42 target). (ii) Extend the depth census to (a1,a_j) beyond (2,3) — the
flatness and the 0.3 constant are so far one-pair measurements. (iii) The
Class-B wall: any effective Borel-Cantelli / shrinking-target input for the
folded ray; the required bound is delta*k against an observed O(1) — the
slack is the invitation.