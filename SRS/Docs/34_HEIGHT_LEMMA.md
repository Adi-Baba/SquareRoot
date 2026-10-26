# Paper 34 — M7.12: the height lemma, and a no-go for the Liouville route

Version 1.0.  Works through the proposed final lemma.  Two corrections to the
project record, a verification of the "everything downstream is mechanical"
claim, and a **provable no-go** for the height/Liouville route the final lemma
was expected to use.

---

## 34.1 Correction 1 — a factor `2` in `Docs/27` §27.2

`Docs/27` displayed
```
q_k √a_j  =  β_j ( ε^k − ε^{−k} ) / 2 ,        β_j := ½√(a_j/a_1) .        ✗
```
The correct identity is
```
q_k √a_j  =  β_j ( ε^k − ε^{−k} )  =  √(a_j/a_1)·Q_k√a_1 / 2 · … ,       ✓
```
verified at `m = 9`: `β_j(ε^9 − ε^{−9}) = 1706.0700454553433` versus
`Q_9√3 = 985√3 = 1706.0700454553441`.

**Consequences: none.**  The *substance* of `Docs/27` — Lemma 27.1,
`‖β_jε^k‖ = ‖q_k√a_j‖ + O(ε^{−k})` — was verified numerically and is
correct, as is everything downstream (`Docs/28`–`Docs/33`, the M7.7/M7.9/M7.11
chain).  Only the displayed formula was wrong.  `Docs/27` corrected in place.

## 34.2 Correction 2 — a pervasive index-shift hazard

The Pell list produced by the standard continued-fraction recursion satisfies
`ε^{j+1} = p_j + Q_j√a_1`, i.e. **list index `j` ↔ exponent `m = j+1`**.  This
offset caused repeated false conclusions in my own verification scripts
(including one spurious "BAD" identity check in this session).  **Convention
fixed henceforth: `q` in the list at index `j` is `Q_j`, and every Pell-index
identity must be written with `ε^{j+1}`.**

---

## 34.3 The `ζ` reformulation, corrected and verified

From `Docs/32` Thm 32.3, an exceptional `k` forces `e_k = ζ(ε')^k` with
`ζ = −Ψ_2/a_0`, hence (taking traces) `ζ` is **determined by `k`**:

> **`ζ_k = 2p_k / ( ε^k + (ε')^k )`**,   and, since `p_k = β_jε^k + e_k`,
> ```
> ζ_k − 2β_j  =  − 2β_j ρ^k/(1+ρ^k)  +  2e_k / ( ε^k (1+ρ^k) ) ,     (34.1)
> ρ := ε'/ε ,   |ρ| = (√a₁−1)/(√a₁+1) < 1 .
> ```

**Correction to the proposed reformulation.**  The rate in (34.1) is *not*
`|ρ|^k` but `θ/ε`: for an exceptional `k`, `|e_k| < θ^k`, so
```
|ζ_k − 2β_j|   ≤   2(θ/ε)^k  +  2β_j|ρ|^k   ≲   2(θ/ε)^k ,            (34.2)
```
and `|ρ| < θ/ε` in every case tested, so the **rounding-error term dominates**.
Verified numerically (exact integer arithmetic, 60-digit magnitudes,
`a_1=2, a_2=3`, `θ₀ = 0.863384`): at `k = 41` (exceptional, `q = 1.7469e15`),
```
ζ_k − 2β_j  =  −1.63691e−18 ,      2|e_k|/ε^k = 1.66e−18 ,      2β_j|ρ|^k = 1.3e−32 .
```
So the observed value is the `e_k`-term, confirming (34.1)–(34.2).  Also
confirmed: `h(ζ_k) = O(k)`, and `ζ_k → 2β_j` exponentially fast.

---

## 34.4 Theorem 34.1 — a NO-GO for the height/Liouville route

**Theorem 34.1.**  *No* argument using only (i) integrality, (ii) height/degree
bounds in `L`, or (iii) Liouville-type inequalities in `L`, can bound the
exceptional indices `k`.

*Proof.*  `ζ_k − 2β_j ∈ L`, `L` of degree `4`, and `h(ζ_k − 2β_j) = O(k)`.
Liouville in `L` gives only a **power law**:
```
|ζ_k − 2β_j|   ≥   c · H(ζ_k − 2β_j)^{1−4}   ≳   c·k^{−3} .            (34.3)
```
But the actual value is **exponentially small**, `≲ 2(θ/ε)^k` by (34.2).
Since `2(θ/ε)^k ≪ k^{-3}` for large `k`, the inequality (34.3) is **always
satisfied** — a power-law lower bound is simply *weaker* than an exponentially
small value, so no contradiction can ever arise, and (34.3) yields no upper
bound on `k`. ∎

**Consequences.**

1. **The proposed "minimal sufficient theorem" is not a shortcut.**
   ```
   |ζ_i − 2β_j|  ≥  exp(−poly(N))                                       (34.4)
   ```
   is, by (34.2), **equivalent to** `k ≤ poly(N)`.  It is not a consequence of
   heights, of Liouville, or of anything else standard — it *is* the statement
   to be proved.  So the route "bound the hyperplane height, then use height
   separation" is **circular**.
2. The obstruction is structural and worth stating: `ζ_k − 2β_j` is an element
   of `L` that is exponentially small **for every `k`** — exactly what units
   produce.  There is nothing arithmetically anomalous about the `ζ_k`.  So no
   *existence* obstruction can exist; the last lemma must be about **which
   `ζ_k` the Subspace Theorem selects**, not about whether `ζ_k` exists.
3. `Docs/31` §31.3's warning is confirmed from a second direction: `o(k)`
   (ST output) versus `O(log k)` (needed) was not merely an artefact — the
   same `k`-versus-`log k` mismatch reappears here as "exponentially small"
   versus "power-law lower bound".

---

## 34.5 Case A is ruled out

The proposed best case was: *the KMN/ST construction already contains a height
bound for its exceptional hyperplanes.*  **This is false, and provably so from
the proof as written.**

* In `Docs/32`/KMN Prop 3.4 the relation `a_0p + Σ_𝒃a_𝒃y_𝒃 = 0` is produced
  by **Theorem 2.1, the Subspace Theorem**.  The `a_𝒃 ∈ L` are the ST's output;
  the classical ST (and every quantitative variant) bounds the **number** `t`
  of hyperplanes, never their **coefficients' heights**.  KMN then discard the
  `a_𝒃` and retain only `η_{i,j}` with `h(η_{i,j}) = o(k)`.
* `ζ_k` is **determined by `k`** (§34.3).  Hence
  ```
  bounding  h(ζ)   ⟺   bounding  k  ⟺   bounding  q*   ⟺   SRS.
  ```
  So "bound the exceptional-hyperplane height" is **not a weaker, more
  accessible version of the problem — it is the same problem.**  This is the
  decisive reason Case A cannot hold in a form that helps.

---

## 34.6 The downstream chain is mechanical — verified

The claim "`k* ≤ P(N)` ⟹ SRS" checks out.  With `q_k ≍ ε^k`,

```
k* ≤ P(N)   ⟹   log₂ q*  =  k*·log₂ε  =  O(P(N))  =  poly(N) ,
log₂(1/γ)   =  (n+1)·log₂Q₀ + O(n)  ≥  (n+1)·log₂q* + O(n)  =  O(N·P(N)) ,
E = log₂(1/|S|)  ≤  O(N²)  +  log₂(1/γ)  =  O(N·P(N))  =  poly(N) .
```

(`Docs/24` Thm 28.5 / `Docs/33`: the `O(N²)` term is the EHS exponent
`2n·log₂(n‖C‖∞)` with `n ≤ N`, `log₂‖C‖∞ ≤ N`.)  So indeed **no further
exponential-rank obstacle waits behind `k*`** — provided the `O(poly(N))`
constant `C = C(a_1,a_j)` in `‖β_jε^k‖ ≥ k^{-C}` may depend on the data
without being controlled in `N` (`Docs/29` §29.5).

---

## 34.7 Where this leaves Case B

Everything special to our configuration is now exhausted:

| ingredient | status |
|---|---|
| ambient ST dimension `= 3` | `Docs/32` Thm 32.1 (effective) |
| explicit ST threshold `n_0` | `Docs/33` Thm 33.1 |
| one exception per hyperplane, computable | `Docs/33` Thm 33.2 |
| effective cardinality `\|ℳ_θ\| ≤ n_0 + t` | `Docs/33` Thm 33.3–33.4 |
| no persistent relation (full contradiction) | `Docs/32` Thm 32.3 |
| `ζ_k` explicit, `h = O(k)`, `ζ_k − 2β_j ≲ (θ/ε)^k` | §34.3 |
| height/Liouville route | **REFUTED** — Thm 34.1 |
| "ST already bounds the hyperplane heights" | **REFUTED** — §34.5 (circular) |

> **M7.12, final.**  The one remaining statement is
> ```
> ∃ C = C(a_1,a_j)  such that   ‖β_jε^k‖ ≥ k^{−C}   for all k ≥ K(a_1,a_j).   (34.5)
> ```
> Equivalently (by §34.3): bound the `k`'s at which the Subspace Theorem's
> exceptional hyperplane of parameter `ζ_k` occurs.  This is **Case B**, and
> §34.4–§34.5 show it cannot be reached by heights, by Liouville, or by
> extracting a height bound from the ST.

**Assessment.**  This is a genuine, sharply-formulated Diophantine-approximation
problem in a *fixed biquadratic field with a two-mode exponential orbit and a
single varying scalar*, with all surrounding infrastructure explicit.  It is
Case B, and on the evidence of `Docs/32`/`Docs/33` the special structure has
removed everything except the effectivity of the ST's own selection step.
Whether an existing low-dimensional effective-ST or `S`-unit theorem reaches
(34.5) is **not determined here**; I am not asserting one way or the other.