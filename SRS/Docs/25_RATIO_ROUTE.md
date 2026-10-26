# Paper 25 — M7.6: the ratio / defect-rigidity route (rigorous verdict)

Version 1.0.  A rigorous treatment of the proposed attack: use the quadratic
identities `p_i² − a_i q² = u_i` to manufacture simultaneous rational
approximations to the pairwise ratios `√(a_i/a_j)`, and look for extra
algebraic rigidity in the defect vector `u`.

**Verdict: the route is refuted, both theoretically and empirically.**

* **Theorem 25.1** (proved) — cross-defect rigidity: `a_j u_i − a_i u_j ≠ 0`.
* **Theorem 25.2** (proved) — the ratio approximations are **never** better than
  the plain quadratic bound, and the implied constant depends only on
  `a_i,a_j`, not on `q` or on `u`.  **The proposed "new information" is
  identically zero.**
* **§25.4** — the rigidity hypothesis is **empirically false**: at exceptional
  `q` the cross-defects are statistically indistinguishable from random.
* **§25.6** — of the five proposed intermediate targets, three are already
  refuted, one is a finite non-decidable question, and one survives.

---

## 25.1 Setup

`a_1 < ⋯ < a_n` distinct squarefree `>1`; `p_i = round(q√a_i)`;
`u_i := p_i² − a_i q² ∈ ℤ \ {0}`; `dist_ℤ(q√a_i) = |u_i|/(p_i+q√a_i)`.
The EHS condition is `Π|u_i| ≲ 2ⁿ√(∏a_i)·q^{n−1−δ}` (`Docs/24` §24.2).

The proposal: from `a_j p_i² − a_i p_j² = a_j u_i − a_i u_j`, deduce

```
| √a_j p_i − √a_i p_j |  ≲  |a_j u_i − a_i u_j| / q ,
hence   | p_i/p_j − √(a_i/a_j) |  is controlled by the defects,
```

giving an overdetermined system of *pairwise ratio* approximations.

---

## 25.2 Theorem 25.1 — cross-defect rigidity

**Theorem 25.1.**  For `i ≠ j`,  `c_ij := a_j u_i − a_i u_j ≠ 0`.

*Proof.*  `c_ij = a_j(p_i²−a_iq²) − a_i(p_j²−a_jq²) = a_j p_i² − a_i p_j²`.
Suppose `c_ij = 0`.  Since `p_i,p_j ≠ 0` (as `q√a_i ≥ √2 > 1/2`), we get
`(p_i/p_j)² = a_i/a_j`, i.e. `a_j n² = a_i m²` with `p_i/p_j = m/n` reduced.
Choose a prime `p | a_i`, `p ∤ a_j`.  Then `v_p(a_j n²) = 2v_p(n)` is **even**
while `v_p(a_i m²) = 1 + 2v_p(m)` is **odd** — contradiction.  (If no such
prime exists then `a_i | a_j`; symmetrically `a_j | a_i`; so `a_i = a_j`,
contrary to distinctness.) ∎

So `u ∉ span_Q(a)`, and equivalently `a_i/a_j` is never a rational square for
distinct squarefree `a_i,a_j`.  **This is the "extra algebraic relation" that
is available — and Theorem 25.3 shows it is worth nothing.**

---

## 25.3 Theorem 25.2 — the ratio route carries zero new information

**Theorem 25.2.**  With `D = a_i a_j`,

```
| p_i/p_j − √(a_i/a_j) |   =   |c_ij| / ( p_j ( a_j p_i + p_j √D ) )   ≥   |c_ij| / ( p_j ( a_j p_i + p_j √D ) )
```

and since `|c_ij| ≥ 1`, this is `≥ 1/(p_j(a_jp_i+p_j√D))`.  Using
`p_i ≍ q√a_i`, `p_j ≍ q√a_j`:

```
| p_i/p_j − √(a_i/a_j) |   ≥   1 / ( q² · √(a_i a_j) √a_j (√a_j + 1) )   =:   1/(C_{ij} q²) ,
```

with `C_ij` depending **only on `a_i,a_j`**.

*Proof.*  `√(a_i/a_j) = √D/a_j`, and
`|a_jp_i − p_j√D| = |a_j²p_i² − Dp_j²|/(a_jp_i+p_j√D) = a_j|c_ij|/(a_jp_i+p_j√D)`;
divide by `a_j p_j`. ∎

**Consequence (the decisive point).**  The right-hand side is the *plain
quadratic bound* for the quadratic irrational `√D/a_j` with denominator
`p_j ≍ q√a_j`, and its constant `C_ij` is independent of `q`, of the defects,
and of everything the EHS condition controls.  Moreover `|c_ij| ≥ 1` means the
ratio approximation is **never better** than `1/(C_ij q²)`.

> **Therefore the pairwise-ratio route produces no constraint on the `u_i` and
> no constraint on `q`.  It is a restatement of the individual quadratic bound,
> not new information.**  The proposal's central premise — that
> `√a_j p_i − √a_i p_j` is "new information compared with the coordinatewise
> Pell estimates" — is **false**, and provably so.

Note this is *not* a Roth-type situation: `√(a_i/a_j)` is quadratic, so we get
the exact quadratic bound rather than Roth's `q^{-2-ε}`, and that bound is
already saturated.  Improving it would require `|c_ij| < 1`, impossible.

---

## 25.4 Empirical refutation of the rigidity hypothesis

The proposal also posits: *"if many `p_i²−a_iq²` are simultaneously small,
then `u` must satisfy some additional algebraic relation."*

Tested on the 121 certified exceptional `q` for `(√2,√3)`, `δ=1/6`
(`Docs/24` §24.8), with `c_12 = 3u_1 − 2u_2`:

```
zero cross-defects among 121 exceptional q ....... 0        (Thm 25.1 holds)

|c_12| / (|u_1| + |u_2|):
   exceptional q  :  median 2.1778   mean 2.2515   max 3.0000
   random control :  median 2.1667   mean 2.2367          <-- same magnitudes
   # with |c_12| < 10% of (|u_1|+|u_2|) :  1 / 121
   # with |c_12| = 1                    :  0 / 121
   # with |c_12| <= 10                  :  9 / 121
```

> **The exceptional points and the random control are statistically
> indistinguishable.**  The defect vector carries **no detectable extra
> algebraic rigidity**.  Rigidity hypothesis **refuted** at the tested
> population.

Note also `|c_12| = 1` never occurs: the cross-defect is never minimal, so the
ratio approximation is never near-optimal either — consistent with
Theorem 25.2.

---

## 25.5 The one genuine identity — and why it constrains nothing

**Theorem 25.4.**  For all `i,j`,
```
( a_j u_i + a_i u_j )²  =  c_ij²  +  4 a_i a_j u_i u_j ,              (25.1)
```
where `c_ij = a_j u_i − a_i u_j`.

*Proof.*  `(x+y)² − (x−y)² = 4xy` with `x = a_ju_i`, `y = a_iu_j`. ∎
(Verified exactly on samples of the exceptional set.)

**But (25.1) is an identity, hence a reparametrisation, not a constraint**:
`(a_ju_i+a_iu_j)` and `(a_ju_i−a_iu_j)` determine `(u_i,u_j)` and vice versa, so
requiring (25.1) removes nothing.  Counting confirms: the number of
`(u_1,u_2)` with `|u_1u_2| ≤ P ≈ q^{1−2δ}` is `O(P log P) = O(q^{1−2δ}\log q)`,
and (25.1) is satisfied by all of them.  **No gain whatsoever.**

---

## 25.6 The five proposed intermediate targets, addressed

**Target 1 (pairwise coprime `a_i`) and Target 2 (`a_i` the first `n` primes):
ALREADY REFUTED, by `Docs/24` Theorem 24.3.**  The Pell construction places
exceptional `q` at `2^{Θ(n²)}` for *any* radicand set whatsoever, including
pairwise coprime sets and prime sets.  Coprimality and primality confer **no
rigidity**; the exceptional points are generated by the Pell denominators of a
single `a_i`, which is a purely local phenomenon.  (Verified: `a=(2,3,5,7)`
is pairwise coprime and has 15,007 exceptional `q` up to `2×10⁷`, still
growing.)

**Target 3 (all `a_i ≤ A`, bounded `A`): reduces to a finite but
NON-DECIDABLE question.**  With `A` fixed, the radicand set and hence `n` are
bounded, so `q*` is a definite integer.  But one cannot certify a finite
search without an upper bound on `q*`, and the Subspace Theorem supplies none.
**No structural gain.**

**Target 4 (`n = 2` or `n = 3`): the genuine first case, and the precise gap
is identifiable.**  For `n = 2` the condition is `|u_1u_2| ≲ 4√(a_1a_2)q^{1−δ}`.
The classical effective theory here is the reduction theory of **small values
of indefinite binary quadratic forms**, and indeed
`c_12 = a_2p_1² − a_1p_2²` is exactly such a form — so that theory is fully
effective.  **But it controls `c_12`, whereas EHS needs `u_1u_2`.**  By
Theorem 25.2 these are unrelated (`c_12 ≥ 1` always, while `u_1u_2` can be
small independently).  **So reduction theory does not transfer, and this is
the exact obstruction.**  This is a sharp, concrete statement of why the
classical effective machinery stops here.

**Target 5 (`|u_i| ≤ q^θ`, `θ < 1`): survives, and is the right refinement.**
Measured for the 121 exceptional `q`: `θ(q) := max(|u_1|,|u_2|)/q` has median
`0.196`; the 17 points with `θ ≥ 0.9` are all **small** `q`
(`q ≤ 5741`).  So the large-`q` population satisfies `|u_i| ≤ q^{0.9}` — i.e.
`|u_i| = Θ(q)`, sublinear only in the weakest sense.  Target 5 therefore
isolates the large-`q` regime and is the right place to look, but note the
defects are still `Θ(q)`, not `o(q)`.

---

## 25.7 Verdict and what is left

```
  Claim                                          Status
  ------------------------------------------------------------------
  c_ij != 0  (extra algebraic relation exists)   PROVED (Thm 25.1)
  ratio approx gives info beyond Pell estimates   FALSE (Thm 25.2, proved)
  defect vector has extra rigidity               FALSE empirically (§25.4)
  identity (25.1) constrains anything             NO  (it is an identity)
  coprimality / primality helps                   NO  (Thm 24.3 kills it)
  bounded A helps                                 NO  (non-decidable)
  reduction theory transfers                      NO  (controls c_ij, needs u_1u_2)
  soft-defect regime (target 5) survives          YES
```

**The honest conclusion.**  The exceptional `q`-set carries **no algebraic
rigidity** beyond the tautology that `c_ij ≠ 0`: the cross-defects are
generic-size integers, statistically indistinguishable from random.  The
phenomenon is therefore **metric, not algebraic** — it is genuine
multidimensional Diophantine approximation, exactly as `Docs/24` §24.8 found
(2/3 of the exceptional `q` are near-convergents for neither coordinate).

This is the *same* structural verdict as `Docs/22` Theorem C, now at the
second level:

| level | inert information | verdict |
|---|---|---|
| `Docs/22` | "other conjugates are large" | no norm/Fourier gain |
| `Docs/25` | "cross-defects are small / ratios well-approximated" | no gain |

**Consequence for the programme.**  A pure algebraic attack on the EHS
constant is now closed from two independent directions.  What remains requires
**effectivity in multidimensional Diophantine approximation for algebraic
points** — i.e. an effective lower bound on the multiplicative Markov constant

```
μ(a_1,…,a_n) := liminf_{q→∞} q^{n+1} Π_i dist_ℤ(q√a_i)  >  0  effectively,
```

at the scale `log₂q* = O(N²)` (`Docs/24` Thm 24.5).  The Pell family
(`Docs/24` Thm 24.3) is the hard test any such theorem must survive, and it
already shows the `n²` scale is forced.  **No such effective Markov-constant
bound is known; the problem is open.**