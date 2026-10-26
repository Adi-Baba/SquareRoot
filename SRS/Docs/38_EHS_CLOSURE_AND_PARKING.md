# Paper 38 — M7.13: EHS parameter-range closure, and the parking of M7.12

Version 1.0.  The last cheap check.  It closes the Track B lever from `Docs/37`
§37.5, and records the formal end of the M7.12 programme.

---

## 38.1 What had to be checked

`Docs/37` Theorem 37.3 gives the product estimate

```
prod_{i=2..n} || q sqrt(a_i) ||   >=   C_n · q^{-(n+1)} ,        C_n > 0 ,
```

and `Docs/37` Fact 37.2 gives `‖Q_k√a₁‖ = ε^{−k} ≍ q^{−1}` exactly.  Since the
EHS exceptional condition is

```
q^{1+delta} · prod_{i=1..n} || q sqrt(a_i) ||   <   1 ,
```

the `a₁` coordinate absorbs the `q^{1}`, leaving the requirement

```
prod_{i=2..n} || q sqrt(a_i) ||   <   q^{-delta} .              (38.1)
```

Theorem 37.3 contradicts (38.1) **iff `n+1 < δ`**.  So the whole of Track B
reduces to a single question about EHS: **is `δ > n+1` available?**

## 38.2 The check, and its answer

**No.**  Eisenbrand–Haeberle–Singer fix

```
delta  =  1 / (3n)
```

in the successive-minima step, and it is not a free parameter.  Their Theorem 4
assumes

```
Q  >=  ( 2n ||x||_inf )^{3/2} ,        lambda_1  >=  Q^{-1 - 1/(3n)} ,
```

and concludes `|E| >= Q^{-(n+1)}` (user's direct reading of the paper;
arXiv:2312.02057, SoCG 2024).  The `1/(3n)` is forced by the
geometry-of-numbers step: the first successive minimum must be separated from
`Q^{-1}` by that specific margin before the Subspace-Theorem machinery applies.
It is not an oversight available for optimisation.

Independently confirmed at the abstract level: the headline result is
`|E| >= gamma · (n · max_i |x_i|)^{-2n}`, improving the previous
`(n · max_i(|x_i|√a_i))^{-2^n}`, with `gamma` **non-explicit** and coming from
the Subspace Theorem.

> **M7.13 (EHS parameter-range closure).**  No available EHS parameter
> optimisation converts the exponent `n+1` of `Docs/37` Thm 37.3 into the
> required exponent `< 1`.  ∎

## 38.3 Why there was never any slack — the scale is forced

The user's sharper observation, which is the real reason this is a closure
rather than a missed opportunity.

`Docs/37` Thm 37.3 returns exponent `n+1 = 1 + n`.  But `q^{−1}` is *exactly*
what the first coordinate contributes (`Docs/37` Fact 37.2), and the remaining
`n−1` coordinates each contribute at best `q^{−1}` under a product-of-conjugates
argument (`Docs/36` Thm 36.1, exponent `1`).  So

```
q^{-1}  ×  q^{-n}   =   q^{-(n+1)}
```

is **precisely the natural arithmetic scale** of the problem.  Theorem 37.3 is
not a loose estimate that a better choice of `δ` could tighten; it is the
correct answer to the question the Pell structure poses.  `Docs/37`'s
`n+1` is not an artefact of a weak method — it is the shape of the arithmetic.

## 38.4 Formal status: M7.12 PARKED

> **M7.12 is PARKED** — *not* because the problem is experimentally exhausted,
> and *not* because M7.12 has been proved impossible.  It is parked because the
> project has isolated the remaining obstruction and excluded every major
> natural mechanism for reaching it.

Mechanisms excluded, with the reference that excludes each:

| mechanism | excluded by | how |
|---|---|---|
| norm / product-of-conjugates | `Docs/36` Thm 36.2 | structural: exponent `≥ 1` vs required `< 1`, all `n ≥ 2` |
| single-index improvement | `Docs/36` Thm 36.2 | degree-4 forced, all three surviving conjugates `~ ε^k` |
| two-index improvement | `Docs/36` Thm 36.1 | reaches exponent `1` but bounds a **sum**, not `d_k` |
| heights / Liouville | `Docs/34` Thm 34.1 | power law is *weaker* than the exponentially small actual value |
| ST hyperplane heights | `Docs/34` §34.5 | circular: `ζ_k` is determined by `k`, so bounding `h(ζ)` **is** bounding `k` |
| Ridout / ST power laws | `Docs/35` §35.6 | `exp(−(1+ε)k log ε)` vs `exp(−3k log ε)`; `1+ε < 3` |
| recurrence invariant | `Docs/36` §36.3 | `I_k ≡ 0` identically; `Q(u,v) = (v−εu)(v−ε'u)` |
| conjugate amplification | `Docs/36` §36.6 | measured flat: `A_k → 2β_j³(2a₁)^4`, ratio `1.00000` |
| product-level Pell coupling | `Docs/37` Thm 37.3 | moves `n−2` coordinates to exponent `1`, leaves **one** at `3` |
| LRS literature route | `Docs/37` §37.1 | Loxton–Van der Poorten family counts solutions, does not locate the last |
| EHS parameter optimisation | `Docs/38` §38.2 (M7.13) | `δ = 1/(3n)` is forced by geometry of numbers |

**What is *not* claimed.**  M7.12 is not proved impossible.  No theorem is
asserted not to exist.  The claim is exactly this: *the algebraic toolkit
developed here cannot reach it, the product and literature routes do not supply
the missing input, and the exponent it needs is below anything
product-of-conjugates can produce.*

## 38.5 The research conclusion

> **SRS is reduced to a genuinely new small-fractional-part problem for a
> trace-zero quadratic recurrence.**
>
> Namely: for distinct squarefree `a₁ ≠ a_j`, `ε = 1+√a₁`,
> `β_j = √(a_j/a₁)/2`, prove that
> ```
> { k >= 1 : || ( sqrt(a_j) / (2 sqrt(a_1)) ) (1+sqrt(a_1))^k ||
>                  <=  (1+sqrt(a_1))^{-k/(3n(n-1))} }
> ```
> is finite with `k*·log ε <= poly(N)`, where the sequence is a **trace-zero
> order-2 linear recurrence** (`Docs/35` Thm 35.3) living in a fixed biquadratic
> field, and where `β_j ∉ Q(√a₁)` (`Docs/35` Thm 35.2) removes the automatic-Pisot
> mechanism that would otherwise make the statement false.

This is a well-posed question about the fractional parts of a quadratic
recurrence with an algebraic coefficient outside the generating field — a
recognised subject, but not one reached by the tools in this repository.  The
surrounding infrastructure is explicit and reusable: the exceptional set is
enumerated (`Docs/33` Thm 33.2), bounded (`Docs/33` Thm 33.4), the threshold is
explicit (`Docs/33` Thm 33.1), and the obstruction is a single hyperbolic
quantity `‖β_jε^k‖` whose distance to the threshold has been measured over
`k = 1…80`.

## 38.6 Programme closure

`Docs/21`–`Docs/38` close the M7 / L-programme.  Final state:

* **M7.1–M7.4** — audits and the two-small theorem (`Docs/19`, `21`, `22`).
* **L1–L7** — literature, q-set, ratio, Markov, Pell-index, field structure
  (`Docs/23`–`29`).
* **L8–L11** — KMN audit, proof audit, relation enumeration, effective count
  (`Docs/30`–`33`).
* **L12** — height lemma; no-go for the Liouville route; Case A refuted
  (`Docs/34`).
* **L12′** — M7.12 sharpened to a quantified `9n(n−1)` gap (`Docs/35`).
* **L12″** — engine built; **structural** no-go for integrality arguments
  (`Docs/36`).
* **L12‴** — Tracks A and B; three new results; product attack does not close
  (`Docs/37`).
* **M7.13** — EHS parameter-range closure; **M7.12 parked** (`Docs/38`).

Norm, invariant, and height experimentation is **frozen**.  The next action for
this line of work, if any, is external: a literature or expert question about
small fractional parts of trace-zero quadratic recurrences — **not** further
algebraic manipulation of `‖β_jε^k‖` by norm, lattice, or height methods, which
§38.4 shows cannot succeed.