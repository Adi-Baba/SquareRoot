# Src

Verification suite for `Docs/Draft1.md`, `Docs/Paper1.md` (audit) and
`Docs/Paper2.md` (theory: Möbius conjugacy and its corollaries).

```
bracket_zoom.py        reference implementation of Algorithms 1-3
test_bracket_zoom.py   35 checks, all exact arithmetic
exp_u_aware.py         the u-aware bound: derivation + exact verification
exp_depth.py           recursion depth, full-interval vs u-aware bound
test_u_aware.py        24 checks pinning every bug found in the two above,
                      incl. the optimal-certificate null result (C vs Theorem 8)
exp_precision.py       precision scheduling vs Newton (experiment 1)
exp_hybrid.py          float-prefix variant (experiment 1b)
exp_compare.py         cubic vs plain vs Newton; where the time goes
exp_mirror.py          reflection across r; mirror + average = Newton
test_exp_mirror.py     11 checks pinning two bugs found in exp_mirror
exp_certified.py       §10.1 certified benchmark: bisection, interval
                       Newton, zoom, Newton+res — all exact-integer certified
test_certified.py      25 checks pinning the §10.1 bugs found in exp_certified
exact_decide.py        exact decision procedures (cubic on rational interval)
verify_theorems.py     proved / decided / measured ledger
mutation_test.py      14 injected faults; 12 killed, 2 analysed in the paper §8.1
opcount.py           exact executed-arithmetic counts (machine independent)
test_opcount.py       7 checks pinning the corrected op-count claims
mobius.py            Mobius structure behind Paper2 Thm 2-3 (exact rationals only)
                      incl. the exact sharpness-attainment factor
run_all.py           registry of every harness + one runner; --list to see it
test_mobius.py        29 checks: closed form vs iteration, no sampling; attainment;
                      Cor 7 slack identity, disjointness, 2/(27(2+t)) constant
mobius_n.py          the n-th root map: error identity, multiplier, 1/n degradation
test_mobius_n.py      19 checks: Thm 10 for n=2..5, Cor 12 sign flip, Thm 13 repelling
test_decimal_oracle.py  5 checks: Decimal vs math.isqrt at up to 6000 digits
test_exact_decide.py  19 checks: deciders vs brute force and direct expansion
```

```bash
cd Src

# Everything, with one command.  --list shows what exists and which paper
# section each harness reproduces; --all includes the slow sweeps.
python3 run_all.py --list
python3 run_all.py            # tests + the fast harnesses
python3 run_all.py --all      # adds exp_certified and mutation_test
python3 run_all.py certified  # just the harnesses matching a substring

# Or drive them individually, which is what the paper's Appendix A does:
python3 -m unittest test_bracket_zoom test_u_aware test_exp_mirror \
        test_certified test_exact_decide test_opcount test_mobius \
        test_mobius_n test_decimal_oracle -v      # 179 pass
python3 verify_theorems.py                       # proved/decided/measured
```

**179 checks total, all passing.** The tests use `Fraction` for every rigorous
claim (exact rationals, so no rounding can decide whether a bound holds) and
`mpmath` only where the integers would otherwise blow up — cubing a rational
every level makes the numerators grow without bound, which is why depth
counting lives in `exp_depth.py` and not in the Fraction module.

`bracket_zoom.py` is numeric-type agnostic — pass `int`, `Fraction`, or
`Decimal`. The tests use `Fraction` for every rigorous claim and `Decimal`
where exact rationals get too slow (past a few hundred digits).

## What the tests cover

| Paper section | Claim | Result |
|---|---|---|
| §2, §3 | Claim 1 weighted average `α = (q·x²+p·y²)/d` | holds for all α |
| §2, §3 | Claim 2 `r = √((x+p/d)² + pq/d²)` | holds exactly |
| §2 | pronic side fact (`x+p/d = HM`, `GM² = AM·HM`) | holds |
| §5 | **Theorem 1** `r − M(a,z) = (r−a)(z−r)/(a+z)` | holds on 150,000 exact cases |
| §5 | plain under-estimates, corrected over-estimates | 0 violations in 20,000 |
| §6 | **Theorem 2** `0 ≤ est − r ≤ w³/(16a²)` | valid |
| §6 | closed form `est − r = u²(w−u)/(2a(2a+w)+u(2a+u))` | exact |
| §8.1 | median 2.00 / 3.00 digits per level | reproduced |
| §8.1 | 2946 and 3152 sample sizes | reproduced |
| §8.1 | bracket invariant never broken | holds |
| §8.2 | `f''(c) = 0` ⇒ 3 digits/level | **confirmed** (open question #4) |
| §8.5 | recursion depth | reproduced within the paper's ±1 convention |
| Paper §3 | Corollary 1, $m_0+v/2m_0>r$ with error $v^2/(4m_0^2(\widehat m_0+r))$ | holds |
| Paper §5.3 | Theorem 4 valid for any $a>0$, so $0<\alpha<1$ is covered | holds |
| Paper §5.3 | Prop. 5 thresholds $w<4a$ and $w<3\sqrt3\,a$ | exact |

## Findings

Tests marked `[FINDING]` assert something *disagreeing* with the paper. Each is
a real defect in the draft, not a test artefact.

1. **Algorithm 3 diverges when `w > 4a`.** `ε = w³/(16a²) < w` ⟺ `w < 4a`.
   Not stated in the paper and not guarded in its code listing. With `a=1, w=5`
   the bracket *widens* 5 → 7.81 and the invariant `a < √α < b` is destroyed
   by the second step. → `TestAlgorithm3Precondition`, `zoom_status()`.

2. **Theorem 2 is loose by 27/16 ≈ 1.69×.** `u²(w−u)` peaks at `u = 2w/3` with
   value `4w³/27`, so the sharp constant is `4/27`, not `1/16`. The ratio to
   the paper's bound is exactly `64a²/(108a² + 90aw + 12w²)`, which approaches
   `16/27` as `w/a → 0` and is well below it when `w ~ a` (`a=w=1` gives 0.305).
   → `TestTheorem2.test_paper_bound_attained_at_two_thirds`.

3. **§8.1's explanation of the extreme-α case is wrong.** The paper says the
   level-0 error is `~1e-25` and blames the measurement. It is `~2.5e-13`, and
   the cause is structural: the digit scan cannot advance `a` at all, because
   `(10⁶ + 10⁻ᵏ)² > 10¹²+1` for every `k ≤ 6`. With `a` frozen the plain error
   `≈ u·w/(2a)` shrinks only linearly in `w`, so the method degrades to exactly
   1 digit/level for six levels. This is an argument *for* the recursive zoom
   (which re-centres `a ≈ r`), not against the method.
   → `TestSection81.test_extreme_alpha_mechanism`.

4. **§8.5's speed ratio is non-monotonic and contradicts its own cost model.**
   The stated model is `O(D log D)` vs `O(D)`, so zoom/Newton must increase
   with D. It goes 5.70 → 8.79 → **4.50**. The 10000-digit row is suspect.
   → `TestReportedNumbers.test_zoom_over_newton_ratio_is_non_monotonic`.

5. **The `isqrt` baseline is in §8.4 but absent from the abstract.** It is
   10–675× faster than the zoom; even the best recursive number is ~10× off the
   C reference at 3000 digits. The conclusion compares only to Newton.
   → `TestReportedNumbers.test_isqrt_baseline_gap_omitted_from_abstract`.

6. **"65× faster" is the best case only.** The §8.5 range is 1.8×–65×, and
   "2 to 9× slower than Newton" is contradicted by the table's own 4.5× at
   10000 digits. → `TestReportedNumbers.test_abstract_speedup_range`.

7. **§8.2 understates the non-π gains.** It reports "~2.0 to 2.1"; the true
   asymptotics are e: 2.74, dottie: 2.66, ln2: 2.30, φ: 2.35. The integers the
   paper reports are digit-scan quantisation plateaus, not the asymptotic rate.

## The u-aware bound (Aditya's observation)

The error is a function of `u` alone: `u = r −a` is the distance from the left,
`w −u` from the right, and

```
r −m   =  u(w−u) / (2a+w)
est −r  =  u²(w−u) / (2a(2a+w) + u(2a+u))
```

`u` itself cannot be evaluated — that would already be the answer — but after
computing `m` and `est` we know `m−a < u < est−a`, an interval far narrower
than `(0,w)`. Maximising the same numerator over the interval we actually know
gives a tighter bound for free.

**Rigour notes, each paid for by a bug:**

* **Same footing.** Both bounds must divide by `D`. An early version let the
  candidate divide by `2a(2a+w)` while the baseline returned bare `4w³/27`,
  which reported "100% tighter". The honest numbers, both from the same
  `_bound` routine, are in `exp_u_aware.py` §2.

* **Never search for the maximiser.** Ternary search on `n(u)/D(u)` returns the
  value at the search endpoint, which is *below* the true maximum — a bound
  that undershoots certifies nothing. Instead: `D` is strictly increasing
  (`D′ = 2a+2u > 0`), so on any sub-interval `n(u)` is at most that
  sub-interval's numerator max and `D(u)` is at least `D` at its left end.
  Splitting the interval gives an upper bound for every piece count, in exact
  rationals. Verified `bound ≥ error` at pieces ∈ {1, 4, 64}, 0 violations in
  30,000 cases.

* **Regime matters.** Tightness depends strongly on `w/a`:

  | regime | median `u-aware / published` |
  |---|---|
  | `w ≤ 100a` | 0.265 |
  | `w ≤ a` | 0.576 |
  | `w ≤ a/100` (realistic) | 0.616 |

* **Depth, measured not assumed.** Saving levels does **not** buy `3^depth`.
  Both schemes must pay for the final level in full, so the honest speedup is
  1.0× at ≤100 digits, 1.14–1.29× at 1000–30000 digits under fixed working
  precision (1.01–1.25× under per-level precision). Full table: `exp_depth.py`.

## The mirror experiment (and two bugs it caught)

The symmetry of `√α` is `x ↦ α/x`, a *reflection* — not `2r − x`. Seeding
`[a, α/a]` brackets `r` by `a ≤ r ≤ α/a` (from `a ≤ √α`, square both sides),
and because `r = GM(a, α/a) < AM(a, α/a)` the bracket halves once, free, with
no test. Measured (`exp_mirror.py`):

| mode | 1000 digits | 3000 digits | depth @3000 |
|---|---|---|---|
| plain `[a, a+1]` | 1.00× | 1.00× | 7.6 |
| mirror `[a, α/a]` | 1.02× | 1.00× | 7.5 |
| narrower of the two | 1.03× | 1.01× | 7.4 |
| mirror + free half | **1.11×** | **1.07×** | 7.2 |

Mirror bracket is narrower in 5 of 16 cases, equal in 2, **wider in 9** — the
`1/1` bracket is already fine, and no scheme saves more than a fraction of a
level, so this is not the missing 3.6×.

The punchline: `(x + α/x)/2` *is* Newton's step. The mirror idea, made
rigorous, converges to the competitor.

**Two bugs, both pinned by `test_exp_mirror.py`:**

1. **Perfect squares hang.** For `α = 4`, `α/a = a` exactly, so the mirror and
   narrow brackets collapse to width 0 → `p = 0`, `est = a`, `eps = 0`, and
   `Decimal(0).adjusted()` is `0`, never `≤ −digits`. The loop never exits
   (`exit=124` on timeout). Plain mode was immune only because it starts at
   `w = 1`. Fixed by returning the exact root at seeding when `a·a == α`.

2. **The free half may only be taken at seeding.** `ab = α` holds only before
   the first zoom; afterwards `a ← est − eps` breaks the invariant. Halving
   *every* iteration passed 0/16 by luck. It is unsound because halving keeps
   the left half `[a, a+w/2]`, valid only while `est − r ≥ eps/2`, which fails
    as `u → 0`. `test_root_can_sit_in_the_{right,left}_half` proves both
    directions in exact rationals.

## The certified benchmark (§10.1)

`exp_certified.py` answers a question §10 could not: what happens when the
*competitors* also carry certificates. Four bracket-carrying methods, identical
α set, identical sizes, adaptive per-step precision, certification outside the
timed region. `certify()` proves `lo² ≤ α ≤ hi²` and `hi − lo ≤ 10^-D` purely
from `Decimal.as_integer_ratio()` — **no reference square root in the check**.

| digits | bisection | int.Newton | zoom (cert.) | Newton+res | N+res(int) |
|---|---|---|---|---|---|
| 50 | 1.272 ms / 167 | 0.088 ms / 5.4 | 0.269 ms / 3.8 | 0.028 ms / 2.0 | 0.035 ms / 6.1 |
| 200 | 10.107 ms / 665 | 0.138 ms / 7.2 | 0.390 ms / 4.9 | 0.043 ms / 3.0 | 0.067 ms / 7.6 |
| 1000 | 219.475 ms / 3322 | 0.408 ms / 9.7 | 1.069 ms / 6.6 | 0.631 ms / 6.0 | 0.671 ms / 9.9 |
| 3000 | 2139.626 ms / 9966 | 1.967 ms / 11.2 | 5.131 ms / 7.6 | 3.005 ms / 7.0 | 5.001 ms / 11.6 |

Step counts land on theory, which is the check that this is not an artefact:
bisection at `3.322·D` to three digits, interval Newton at `log2 D`, the zoom
at `log3 D`. **320/320 runs certified** (4 sizes × 16 α × 5 methods);
**480/480** on a six-size sweep.

**The hoped-for claim is refuted.** On identical seeds the zoom is 2.61× slower
than interval Newton at 3000 digits — it takes *fewer* steps (7.6 vs 11.2) and
still loses, because cubic convergence buys fewer levels, not less work. It is
417× faster than bisection and 1.71× slower than Newton+res.

**Three bugs and two methodological errors, all pinned by `test_certified.py`:**

1. **The termination test was 10× too lax.** `eps.adjusted() <= -digits`
   accepts every `ε` in `[10^-D, 10^(1-D))` — a bracket up to ten times wider
   than promised. `α = 257` at 50 digits lands exactly there
   (`1.198e-50 > 1e-50`), and 1 of 64 points failed before the fix. Strict `<`
   in all five files that carried it: **1/64 → 0/64**.

2. **A bracket that certified before being printed and failed after.** The
   exact identity `r ∈ [est−eps, est]` holds for the *exact* `est`. Round the
   endpoints naively and `est`'s rounding uncertainty vanishes from the printed
   interval: **16 of 48 certificates failed** while the invariant itself held
   48/48 — the loss was purely in extraction. `certified.py` builds
   `_est_bracket`, an outward enclosure with directed rounding: 48/48.

3. **Certification was inside the timed region.** An intermediate version of
   `newton_certified` called `certify()` in-loop to decide whether it had
   converged, adding ~2.7 ms of a 9.8 ms row. Verification belongs in
   `verify()`, outside the timing; the in-loop check that is *allowed* is a
   single Decimal subtraction `hi_b - lo_b > 10**-digits`, which is what
   catches the pessimistic digit counter falling one step short (30 of 480
   runs failed without it).

4. **An under-scheduled competitor.** Interval Newton run with the `rate=1`
   precision schedule is rounding-limited to ~16 digits/step and needs **65
   iterations for 1000 digits**; at `rate=2` (next width ~w²) it needs 9.7. The
   first run of this benchmark made a sound method look absurd.

5. **A seeding asymmetry that moved the ranking.** Bisection, interval Newton
   and the zoom all start at `floor(√α)` with ~3 correct digits — the integer
   bracket is part of the zoom's definition. Newton is conventionally seeded
   from a 30-digit `sqrt`, a 27-digit head start the others never get. Worth
   **1.24×/1.55×/1.06×/1.66×** at 50/200/1000/3000 digits: it cuts Newton from
   11.6 steps to 7.0 but the removed steps are the cheap ones. On the shared
   seed the zoom is at *parity* with Newton (5.131 vs 5.001 ms) instead of
   1.71× behind — so the headline gap to Newton was a seeding artefact, not an
   algorithmic one.

The harness **interleaves methods per α** with best-of-N timing. An early run
measured every method ~1.8× slower than a later one purely from machine load;
sequential per-method timing would have charged that drift to whichever method
happened to run last. Ratios, not milliseconds, are what the paper cites.

Certifying the zoom costs `8.4×` at 50 digits and `5.2×` at 200 — a fixed
~0.4 ms of outward-rounded `Decimal` context setup — but only `1.2×` at 1000
and is within noise at 3000, where recursion depth is *identical* to the plain
§10 path (3.75, 4.94, 6.56, 7.56). At 3000 digits the certified path is
marginally **faster**, because the plain path recomputes the terminating level
at full target precision and the certified path does not.

## A bug the suite caught in the paper itself

`Corollary 1` was first written with the error term
$v^2/(4m_0(\widehat m_0+r))$, missing a factor of $m_0$. The suite's numeric
check failed by a factor of exactly $m_0$ and the paper was corrected to
$v^2/(4m_0^2(\widehat m_0+r))$.

## The missing theorem (open question #4, now settled)

The paper leaves "π wale `sin″(π) = 0` argument ko test karna" as future work.
It holds, and the same argument explains every row of §8.2 — including the two
the paper does not explain. For a simple root `c` with `0 < |f'(c)| < ∞`, and
`A = c−a`, `B = b−c`:

```
|c − x_L|  =  (A·B / 2) · |f''(ξ) / f'(c)|
```

* `f''(c) ≠ 0` ⇒ `A, B = O(w)` ⇒ `O(w²)` ⇒ **2 digits/level** — e, ln2, φ, ∛2, dottie
* `f''(c) = 0` ⇒ `f''(ξ) = f'''(c)(ξ−c) = O(w)` ⇒ `O(w³)` ⇒ **3 digits/level** — π, `tanh(30(x−π))`
* `|f'(c)| → 0` or `∞` ⇒ the `w²` term dies ⇒ `O(w)` ⇒ **1 digit/level** — flat and steep roots

This is worth promoting to Theorem 3. See
`TestInterpolationRateTheory` for the measurements.
