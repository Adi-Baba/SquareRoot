# SRS — Square-Root Sum Probe

A modular, exact-arithmetic toolkit for the **Square-Root Sum Problem (SRS)**:
decide the sign of
```
S = Σ_{i=1}^k √a_i − Σ_{j=1}^k √b_j
```
on a Turing machine (bit-complexity model), where each input integer has ≤ n bits
and `N = k·n` is the total input size.

## Documentation

| file | what it is | read it |
|---|---|---|
| [`STATE_OF_THE_ART.md`](STATE_OF_THE_ART.md) | the literature check, and the correction that changed the answer | **first** |
| [`THE_MISSING_THEOREM.md`](THE_MISSING_THEOREM.md) | the roadmap: problem, both proved sides, the exact open problem, seven routes, every retraction | second |
| [`PROOFS.md`](PROOFS.md) | the proofs | for verification |
| [`TODO.md`](TODO.md) | what is decided vs still open, with stopping rules | for status |

**The one-paragraph version.** `sign(S)` for `S ≠ 0` is open, but the obstacle
is **not** an exponential barrier. Eisenbrand–Haeberle–Singer (arXiv:2312.02057;
SoCG 2024) apply Schmidt's Subspace Theorem to this exact problem and get
`|E| ≥ γ·(n‖x‖_∞)^{−2n}`, exponent sharpenable to `n+ε` — **polynomial in the
number of terms**. The `2^r` in `PROOFS.md` Thm 2.6 is a property of the
**norm**, not of the problem. What remains open is the **effectivity of the
constant `γ`**, which depends on the radicands, and those are input.

## The problem, stripped

Written proofs for the theoretical claims below live in **`PROOFS.md`**:
the independence behind the zero test (Thm 1), the separation bounds (Thm 2 —
general Liouville in 2.1, the sharp two-term constant in 2.2, the orbit
second moment in 2.3, the quadratic tower in 2.4 with its support rank in
2.4′, the sharpest norm bound via the global-flip pairing in 2.6, the achieved
rate by counting in 2.5, and the Subspace-Theorem bound in 2.8 — the last
**quoted from the literature, not proved here**), and the exact F1 law with its
explicit residual (Thm 3).
Computation corroborates; it does not carry these claims.
**`THE_MISSING_THEOREM.md`** is the map: both proved sides, the exact open
problem, all seven routes, and every retraction.

| Variant | Complexity | Why |
|---|---|---|
| `S = 0` | **P** | Squarefree reduction + coefficient aggregation → hash map check (distinct squarefree parts are ℚ-linearly independent — **proved: `PROOFS.md` Thm 1**) |
| `sign(S)` for `S ≠ 0` | **Open** | Root separation. **The exponent is settled: `E/N ≤ (r−1+δ)/r`** (`PROOFS.md` Thm 2.9, single-form Subspace Theorem), matching the achieved rate `(r−1)/r` of Thm 2.5. So the sharp conjecture is the *asymptotic* `limsup_{N→∞} E/N = (r−1)/r` and both sides meet. (The `sup` form is **retracted**: a certified `r=3` instance gives `E/N = 1.2795 > 2/3` — `PROOFS.md` §2.10.2.) What is open is the **effectivity of the constant** `c(d,δ)`, which depends on the input radicands. The `k = 2` constant is closed and sharp at `1/2` **with an explicit constant** (this repo: Thm 2.2 + 2.6). Best known complexity class containing the sign problem is PSPACE. |

**Precondition for any meaningful statement:** input integers are bounded — `a_i, b_j ≤ 2^n − 1`, `k` terms, total bit-size `N = k·n`. Without this, `√(M²+1) − M ≈ 1/(2M)` makes the problem vacuously "small separation" for unbounded `M`. All instance generators in this repo return `N` explicitly.

## What we can and cannot do

| Can do (proved / implemented) | Cannot do (open / out of scope) |
|---|---|
| Exact zero test for `S = 0` (in P) | Prove `2^{−poly(N)}` lower bound for all `S ≠ 0` |
| Exact sign for small `k` via exact algebraic isolation | Prove SRS ∈ P on a Turing machine |
| Certified `√α` to `D` digits (exact interval, no reference sqrt) | Guarantee polynomial precision suffices for all instances |
| **Proved** `log₂(1/|S|) = N/4 + 1 − log₂|D| + ρ`, `0 < ρ ≤ 1.31e-3` (F1 — **`PROOFS.md` Thm 3 / Cor 3.3**, exact with constants) | Rule out existence of `2^{−N²}` or `2^{−2^N}` separations |
| **Proved** sharpest *norm* bound `log₂(1/|S|) ≤ ((2^{r−1}−1)/2)(log₂M+1)`, `r` = F₂-rank (**`PROOFS.md` Thm 2.6**) — the optimum of the norm class, and the best fully *effective* bound | Make the Subspace-Theorem constant `γ` effective: bound `Q₀(a_1..a_n)` with `log₂Q₀ ≤ poly(N)` (**`THE_MISSING_THEOREM.md` §3.2**) |
| **Proved** sharp `k = 2` constant: `E/N → 1/2`, tight from above and below (**Thm 2.2 + 2.6**) — **stronger than the literature at this rank**, which does not use the radicands | Prove `2^{−poly(N)}` — needs only the effectivity of one constant, not a new `ε`-improvement |
| **Proved here** `E/N ≤ (r−1+δ)/r` — the single-form Subspace Theorem, exponent `r−1`, **equal to the achieved rate at every rank** (**`PROOFS.md` Thm 2.9**). At `r=2` this is Roth, confirmed sharp at `1/2` by Thm 2.2/2.6 | Bound the constant `c(d,δ) ≤ 2^{−poly(N)}`, and describe the exceptional rational subspaces (**detectable**, so the generic case is already constant-free) |
| **Quoted from the literature** `\|E\| ≥ γ·(k‖x‖_∞)^{−2k}` (arXiv:2312.02057, SoCG 2024) — the *simultaneous* application, exponent `2n`, so `1.5×` looser than Thm 2.9 at `k=3` | Decide SRS in polynomial time — needs a *computable* constant to set a precision threshold |
| **Searched, not found** no instance with a *growing* `E/N`, certified signs and magnitudes, `r=2..5` (`refutation.py`) — a negative result, reported as such. The harness's ceiling is the *asymptotic* rate; finite-`N` instances above `(r−1)/r` exist and are expected (additive slack) | Prove `limsup_{N→∞} E/N = (r−1)/r`, the sharp form |
| Rigorous measurement of `log₂(1/|S|)/N` across families | Prove the sign problem is in NP (best known upper bound is PSPACE) |

**What this repo actually delivers:**
- A rigorous one-sided bound: `∃ instances with |S| = 2^{−Θ(N)}` (proved, exact — **`PROOFS.md` Cor 3.3**, with the `O(1)` explicit and one-sided)
- **A settled exponent.** Thm 2.9 puts the guaranteed rate at `(r−1+δ)/r`, equal to the achieved rate — so the sharp conjecture `limsup_{N→∞} E/N = (r−1)/r` is a *constant* question, not a growth-rate one (**`PROOFS.md` §2.9, §2.10**)
- **The constant is a Markov constant** (`PROOFS.md` §2.10.1). At `r=2` it is `1/M(α)` for the quadratic irrational `α=√(d₂/d₁)`, hence computable from a periodic continued fraction — that is *why* `r=2` closes. For `r≥3` the periodicity is gone.
- **The constant is effective for `r ≤ 3`** (`PROOFS.md` Thm 2.20): `c ≥ 1/(2√D)` at `r=3`, so `log₂(1/c) ≤ 1.5 log₂D + 1 = O(N)`; at `r=2` no constant is needed at all (Thm 2.2, below). At `r=3` the two-variable form factors through a **single** quadratic irrational, which is the whole mechanism. **`r ≥ 4` is OPEN** — an `r≤4` generalisation was claimed and **retracted** (§2.11.9): it rested on `dist(aW,Z) = a·dist(W,Z)`, false unless `aZ = Z`.
- **Is the conjecture proved? NO.** Proved for `r ≤ 3` and for `r ≤ (K−1)·log₂N` for each fixed `K`; the conjecture needs one fixed `K` for all `r ≤ N/2`. The region `4 ≤ r ≤ N/2` with `r ≿ log₂N` is open, and the only obstacle there is the single constant `log₂(1/c)`.
- **The exact obstruction** (`PROOFS.md` §2.10.3, §2.11): the elementary norm bound is `|S| ≥ (3T)^{−(2^{r−1}−1)}` — effective, but exponential in `r`, and `r = Θ(N)` is achievable. Thm 2.9 is linear in `r` but ineffective. The open step is `poly(r)` **and** effective at once.
- **The exponent is the algebraic degree** (`PROOFS.md` Thm 2.17). `[Q(S):Q] = 2^r` *exactly*, so no norm argument can ever reach `poly(r)` — not by a better estimate, not by different representatives. Sharp: zero one `Cᵢ` and the degree is exactly `2^{r−1}`. Progress needs a **third family of argument** that is neither a norm nor the Subspace Theorem.
- **The conjecture is proved for `r = O(log N)`** (`PROOFS.md` Cor 2.19), effectively and uniformly in `d`. With `r ≤ N/2` (Thm 2.18) the open band is `O(log N) ≤ r ≤ N/2`, and by the radicand-magnitude dichotomy the danger is exactly *small radicands at large `r`*.
- **That exponent is irreducible** (`PROOFS.md` §2.11.2–2.11.3). It is exactly the number of complement pairs `2^{r−1}−1`, forced by AM-GM; and squaring `S` to halve the field degree does not help, because the rational part `ΣCᵢ²dᵢ` always destroys the complement pairing. No variant of the norm argument reaches `poly(r)`.
- **But every fixed rank is already settled** (`PROOFS.md` §2.11.1, Thm 2.12): `E ≤ (2^{r−1}−1)(N + log₂(3r))` is effective, and `2^{r−1}−1` is a *constant* at fixed rank, so the conjectured polynomial separation is a theorem rank-by-rank. Only **uniformity in `r`** is open.
- **`r = 2` needs no constant at all** (`PROOFS.md` §2.11.4, Thm 2.2): the two-term bound `log₂(1/|S|) ≤ N/2 + 3/2` is effective outright, with zero slack in the constant. (It does **not** come from `M(α) ≤ 1` — that bound is false, see §2.11.4; `M` is `Θ(1)` and computable from the periodic CF, and it is not needed here.) At `r = 3` the constant measures `≤ 22` bits over `d ≤ 20000` with median `5.8` — evidence that `c = Ω(1)`.
- **An active refutation hunt** (`refutation.py`) that reports no counterexample with certified magnitudes, and treats "nothing found" as a reproducible negative rather than a shrug
- **A precise statement of what is missing** — the one constant that remains open, seven routes tried, and every error retracted (**`THE_MISSING_THEOREM.md`**, **`STATE_OF_THE_ART.md`**)
- A measured frontier: `log₂(1/|S|)/N` across 6 families, no counterexample found
- A two-sided experiment: F1–F3 map the frontier; F4/F5 hunt for super-polynomial separation; **F6 attains the `r=3` ceiling** and so calibrates what "at the ceiling" looks like
- **Not a solution.** If the conjecture is false, F4/F5 are the only path to a real result. Note that the *guaranteed* side is no longer the obstacle: the best known bound is polynomial modulo an ineffective constant, so the repo's effective `O(2^r)` is a limitation of the implementation, not of the theory.

## Module map (dependencies strictly downward)

```
SRS/
├── squarefree.py     # L0: exact a = c²·s decomposition, aggregation, zero test
├── field.py          # L0: multiquadratic field, XOR basis, Norm, independence assert
├── f1_theorem.py     # L0: PROVED separation law (exact, not measured)
├── signexact.py      # L1: exact sign via an integer isqrt bracket kernel
├── certsingle.py     # L2: certified √α bracket (wraps ../Src/exp_certified)
├── intervalsum.py    # L3: certified rational bracket for S, relative width
├── families.py       # L4: F1–F6 instance generators
├── refutation.py     # L5: certified hunt for a counterexample to the asymptotic rate
├── probe.py          # L5: harness, tables, growth fits
├── run.py            # registry + runner (mirrors ../Src/run_all.py)
├── _srclib.py        # loads ../Src/ modules by explicit file path
├── _deps.py          # re-export shim over _srclib (no production callers)
└── tests/            # one per module, + test_no_float.py and test_probe.py
```

`signexact.py` used to be described here as "exact sign via conjugate poly +
root isolation". It does neither: its kernel is `t = isqrt(s << 2b)`, which
certifies `t/2^b ≤ √s ≤ (t+1)/2^b` and lets both sign tests be integer
comparisons. No polynomial is formed and no root is isolated. `field.py` is
where minimal/characteristic polynomials live, and it is not on the sign path.

The Rust crate is a **sibling** of `SRS/`, not a child — the path is
`../Rust/srs/`, capital R, which matters on a case-sensitive filesystem:

```
../Rust/srs/
├── Cargo.toml
├── src/
│   ├── lib.rs        # re-exports
│   ├── error.rs      # error types
│   ├── field.rs      # multiquadratic over BigInt
│   ├── squarefree.rs # exact decomposition, Brent rho
│   ├── signexact.rs  # certified sign decision
│   ├── f1.rs         # Pell convergents at 10⁷ bits, exact identity assertion
│   └── bin/probe.rs  # CLI harness (mirrors python probe.py)
```

## Design constraints (baked in, not optional)

1. **No float in any decision path.** `test_no_float.py` AST-scans for `float(` / `math.sqrt` — same source-contract as `../Src/test_certified.py`.
2. **Input size is a precondition.** Every `families.py` return carries `N = k·n`. No `N` = no comparison.
3. **Two independent sign paths must agree.** `signexact` (algebraic, small `k`) vs `intervalsum` (certified brackets, any `k`). Disagreement = **bug**, not warning.
4. **Witnesses on every rejection** (per `../Src/exact_decide.py`).
5. **Squarefree honesty, stated per size.** `squarefree_decompose` proves
   `c·c·s == n` on every return path. Whether `s` is squarefree then depends
   on one thing: if the cofactor is below 2⁶⁴ the 7-base Miller–Rabin set is a
   proof and squarefreeness is **proven**; above 2⁶⁴ the test is only probable,
   and that is said in the docstring rather than glossed. Generated radicands
   avoid the question entirely — `random_squarefree` builds them as products of
   distinct primes each below 2⁶⁴, so an in-repo instance is certified at every
   width.
6. **No duplication of certified core.** `certsingle.py` wraps
   `../Src/exp_certified.py` (already mutation-tested), loaded by explicit file
   path through `_srclib.load_src_module` rather than `sys.path` mutation.
   `_deps.py` re-exports the same names for external callers; nothing in
   `SRS/` imports it.

## Run

```bash
cd /home/aditya/Documents/Irrational/SRS
python3 run.py --list
python3 run.py                  # tests + fast harnesses
python3 run.py --all            # all 14 Src harnesses too; budget ~15-20 min
python3 -m pytest tests -q      # 347 tests, 12093 subtests
python3 refutation.py            # refutation hunt: no counterexample found
python3 refutation.py --r 3 --budget 400 --bits 80   # deeper, at rank 3
python3 probe.py --max-bits 4096  # precision cap is threaded and enforced
```

`--all` used to run only 2 of the `Src/` harnesses while still printing
"SRS FULL SUITE"; it now runs all of them, which is why it takes minutes rather
than seconds. `python3 probe.py --max-bits 8` exits non-zero, because an
undecided run must not look like a clean one.

## Rust crate

Build: `cd ../Rust/srs && cargo run --release --bin probe -- --a 2 --n 1 --bits 1000000`
Tests: `cd ../Rust/srs && cargo test && cargo clippy --all-targets`

The crate has no F4/F5 generators and no general `deg(minpoly)` sweep, so it
covers a strict subset of what the Python side does.

## Stopping rule (from TODO.md)

> If F1–F6 all hold ratio `log₂(1/|S|)/N < 2` for all instances up to `N ≈ 2000`, declare **no counterexample found**, write up as labelled negative measurement, **stop**. Prevents an open-ended research hole.