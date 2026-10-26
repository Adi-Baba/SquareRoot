# The Möbius-CF Exact Correspondence for Pronic Numbers

## Abstract

We prove an exact algebraic correspondence between the rational iterates of the Möbius map φ(z) = a + a/(a+z) and the continued fraction convergents of √α for pronic numbers α = a(a+1). For all n ≥ 0, the n-th Möbius iterate zₙ = φⁿ(a+1) satisfies:

- **Odd n:** denominator equals the CF numerator
- **Even n:** a × denominator equals the CF numerator
- **Odd n:** numerator equals a(a+1) × CF denominator
- **Even n:** numerator equals (a+1) × CF denominator

This identity provides certified rational approximations with geometric
convergence (ratio |λ| = (r−a)/(r+a) < 1) and exact error bounds via the
Möbius conjugacy.

---

## 1. Introduction

Let α = a(a+1) for integer a ≥ 2 (pronic numbers: 6, 12, 20, 30, 42, ...). The continued fraction of √α is periodic with period 2:

```
√α = [a; 2, 2a, 2, 2a, 2, 2a, ...]
```

The convergents pₙ/qₙ satisfy the recurrence:
```
p₀ = a,          q₀ = 1
p₁ = 2a+1,       q₁ = 2
pₙ = aₙpₙ₋₁ + pₙ₋₂,   qₙ = aₙqₙ₋₁ + qₙ₋₂
where aₙ = 2 if n odd, 2a if n even (for n ≥ 1)
```

Consider the Möbius map derived from the cubic bracket zoom method (Paper 1):
```
φ(z) = a + a/(a+z) = (a z + a²+a) / (z + a)
```
with matrix representation M = [[a, a(a+1)], [1, a]].

Starting from z₀ = a+1 (the right endpoint of the integer bracket [a, a+1]), define the exact rational iterates:
```
z₀ = a+1
zₙ = φ(zₙ₋₁)  for n ≥ 1
```

**Theorem 1 (Exact Möbius-CF Correspondence).** For α = a(a+1) and zₙ = φⁿ(a+1) with zₙ = Aₙ/Bₙ in lowest terms:
1. **Odd n = 2k+1:** B_{2k+1} = p_{2k+1}  (denominator = CF numerator)
2. **Even n = 2k:**   a·B_{2k} = p_{2k}   (a × denominator = CF numerator)
3. **Odd n = 2k+1:** A_{2k+1} = a(a+1) · q_{2k+1}  (numerator = α × CF denominator)
4. **Even n = 2k:**   A_{2k} = (a+1)·q_{2k}     (numerator = (a+1) × CF denominator)

Equivalently, the exact iterates are:
```
z_{2k+1} = a(a+1)·q_{2k+1} / p_{2k+1}
z_{2k}   = (a+1)·q_{2k}   / (p_{2k}/a)
```

---

## 2. Proof

### 2.1 Matrix Setup

The Möbius map φ corresponds to the matrix:
```
M = [[a, a(a+1)], [1, a]]
```
The n-th iterate is zₙ = Mⁿ · (a+1) in Möbius action:
```
zₙ = (Mⁿ₁₁(a+1) + Mⁿ₁₂) / (Mⁿ₂₁(a+1) + Mⁿ₂₂)
```

The continued fraction of √(a(a+1)) has period 2 with partial quotients:
```
a₀ = a,  a₁ = 2,  a₂ = 2a,  a₃ = 2,  a₄ = 2a, ...
```

The CF period matrix for one period (2, 2a) is:
```
C = [[2, 1], [1, 0]] · [[2a, 1], [1, 0]] = [[4a+1, 2], [2a, 1]]
```

### 2.2 The matrix identity (corrected)

An earlier version of this section claimed

```
M² = a·C + (a-1)I,   C = [[4a+1, 2], [2a, 1]]
```

**That identity is false.** At a = 2 it reads M² = [[10, 24], [4, 10]] against
a·C + (a-1)I = [[19, 4], [8, 3]]; the intermediate expansion below it was
correct but the last step does not follow. The correct identities are:

**Lemma 1.** M² = 2a·M + a·I (Cayley–Hamilton: tr M = 2a, det M = −a give the
characteristic polynomial λ² − 2aλ − a). Equivalently, folding α = a(a+1)
into the entries,

```
M² = a·K,   K = [[2a+1, 2a(a+1)], [2, 2a+1]]
```

*Proof.* Direct expansion:
```
M² = [[a, a²+a], [1, a]]² = [[a² + a²+a, a(a²+a) + a(a²+a)], [2a, a² + a²+a]]
   = [[a(2a+1), 2a²(a+1)], [2a, a(2a+1)]]
   = 2a·[[a, a²+a], [1, a]] + a·[[1, 0], [0, 1]]
   = a·[[2a+1, 2a(a+1)], [2, 2a+1]]
```
∎

K is the two-step matrix in (q, p) coordinates. Substituting zₙ = αqₙ/pₙ
into φ² (whose matrix is M²) gives

```
φ²(zₙ) = α·[(2a+1)qₙ + 2pₙ] / [2α·qₙ + (2a+1)pₙ],
```

so the pair (qₙ, pₙ) is advanced two CF steps by C̃ = [[2a+1, 2], [2a(a+1), 2a+1]]
— exact for every n, and consistent with the convergent recurrences
(`[q_{n+2}; p_{n+2}] = C̃·[q_n; p_n]`, verified alongside Theorem 1 below).

### 2.3 Eigenvalue relation

The eigenvalues of M are λ₊ = a + √(a²+a) and λ₋ = a − √(a²+a).
The eigenvalues of C̃ (and of C) are μ₊ = (2a+1) + 2√(a²+a) and
μ₋ = (2a+1) − 2√(a²+a).

**Lemma 2.** μ = 2λ + 1.

*Proof.* Direct: 2(a ± √(a²+a)) + 1 = (2a+1) ± 2√(a²+a) = μ.
∎

### 2.4 Three exact facts about the convergents

Write

```
Dₙ = pₙ² − αqₙ²   (the Pell residue),
eₙ = pₙpₙ₊₁ − αqₙqₙ₊₁   (the cross term).
```

**Lemma 3 (residue law).** Dₙ = +1 for odd n, −a for even n.

**Lemma 4 (cross term).** eₙ = −a·(−1)ⁿ.

**Lemma 5 (determinant).** pₙ₊₁qₙ − pₙqₙ₊₁ = (−1)ⁿ.

*Proof (joint strong induction; bases at n = 0 and 1).* Each quantity at
level n+1 expands, via pₙ₊₁ = aₙ₊₁pₙ + pₙ₋₁ and qₙ₊₁ = aₙ₊₁qₙ + qₙ₋₁, into
the same quantities at levels n and n−1 plus explicit constants:

```
Dₙ₊₁ = aₙ₊₁²·Dₙ + Dₙ₋₁ + 2aₙ₊₁·eₙ₋₁
eₙ   = aₙ₊₁·Dₙ + eₙ₋₁
pₙ₊₁qₙ − pₙqₙ₊₁ = pₙ₋₁qₙ − pₙqₙ₋₁   (the aₙ₊₁ terms cancel)
```

*Bases:* D₀ = a² − α = −a; D₁ = (2a+1)² − 4α = 1; e₀ = a(2a+1) − 2α = −a;
e₁ = p₁p₂ − αq₁q₂ with p₂ = 2a(2a+1)+a = 4a²+3a, q₂ = 4a+1, which expands to
(2a+1)(4a²+3a) − 2a(a+1)(4a+1) = a; and p₁q₀ − p₀q₁ = (2a+1) − 2a = 1 = (−1)⁰.

*Steps, by parity* (using the hypotheses at n and n−1; recall aₙ₊₁ = 2 on
even n, 2a on odd n):

- n even: Dₙ₊₁ = 4·(−a) + 1 + 4a = 1 ✓ (n+1 odd); eₙ = 2·(−a) + a = −a ✓;
  the determinant identity flips sign ✓.
- n odd: Dₙ₊₁ = 4a²·1 + (−a) + 4a·(−a) = −a ✓ (n+1 even); eₙ = 2a·1 − a = a ✓. ∎

(An equivalent route to Lemma 3 is the classical d-sequence of the CF
square-root algorithm: for α = a²+a it is periodic with d = 1 on even steps
and a on odd steps, giving Dₙ = (−1)ⁿ⁺¹dₙ₊₁. Both routes verified for
a = 2..29.)

### 2.5 Proof of Theorem 1

The four items are reductions of one unified identity:

**Claim.** zₙ = α·qₙ/pₙ for every n ≥ 0.

*Base.* z₀ = a+1 = αq₀/p₀ = a(a+1)/a ✓. z₁ = a + a/(2a+1) = 2a(a+1)/(2a+1)
= αq₁/p₁ ✓.

*Step.* Assume zₙ₊₁ = αqₙ₊₁/pₙ₊₁. Then

```
zₙ₊₂ = φ(zₙ₊₁) = a + a/(a + αqₙ₊₁/pₙ₊₁)
     = α(pₙ₊₁ + a·qₙ₊₁) / (a(pₙ₊₁ + (a+1)qₙ₊₁))
     = (a+1)(pₙ₊₁ + a·qₙ₊₁) / (pₙ₊₁ + (a+1)qₙ₊₁),
```

using α = a(a+1). It remains to show this equals αqₙ₊₂/pₙ₊₂, i.e.

```
(a+1)(pₙ₊₁ + a·qₙ₊₁)·pₙ₊₂ = α·qₙ₊₂·(pₙ₊₁ + (a+1)qₙ₊₁).
```

Substituting the recurrences pₙ₊₂ = aₙ₊₂pₙ₊₁ + pₙ, qₙ₊₂ = aₙ₊₂qₙ₊₁ + qₙ and
expanding, all terms cancel except

```
LHS − RHS = (a+1)·[aₙ₊₂·Dₙ₊₁ + eₙ − a·(−1)ⁿ].
```

Both parities give zero, by Lemmas 3 and 4:

- n even: aₙ₊₂ = 2a, Dₙ₊₁ = 1, eₙ = −a:  2a − a − a = 0.
- n odd:  aₙ₊₂ = 2,  Dₙ₊₁ = −a, eₙ = a:  −2a + a + a = 0.

So zₙ₊₂ = αqₙ₊₂/pₙ₊₂, and the claim holds for all n by induction. ∎

*Items 1–4 from the claim.* Let g = gcd(αqₙ, pₙ), so Aₙ = αqₙ/g and
Bₙ = pₙ/g in lowest terms. Writing pₙ = gP and αqₙ = gG, the residue
factors as Dₙ = pₙ² − αqₙ² = g(gP² − Gq), so **g divides Dₙ**.

- Odd n: Dₙ = 1 forces g = 1. Hence Bₙ = pₙ (item 1) and Aₙ = αqₙ (item 3).
- Even n: g | a. For the reverse divisibility, note q₂ₖ ≡ 1 (mod a): on even
  n the recurrence reads qₙ₊₂ = 2a·qₙ₊₁ + qₙ ≡ qₙ (mod a), and q₀ = 1. Then
  pₙ² = αqₙ² − a = a·[(a+1)qₙ² − 1] with (a+1)qₙ² − 1 ≡ qₙ² − 1 ≡ 0 (mod a),
  so a² | pₙ², hence a | pₙ, hence a | g, hence g = a. Therefore
  a·Bₙ = pₙ (item 2) and Aₙ = αqₙ/a = (a+1)qₙ (item 4). ∎

**Computational verification retained as a cross-check:** the formulas hold
exactly for all a ∈ [2, 100] and n ≤ 20 in `fractions.Fraction`
(`Src/mobius_cf_theorem.py`, `SRS/f1_theorem.py`), and the two-step identity
`[qₙ₊₂; pₙ₊₂] = C̃·[qₙ; pₙ]` holds for a = 2..29, n ≤ 25. These now
corroborate the proof above; they are no longer its only support.

---

## 3. Consequences

### 3.1 Certified Geometric Convergence

The Möbius conjugacy T(z) = (z - r)/(z + r) conjugates φ to scaling:
```
T(φ(z)) = λ T(z),  λ = -u/(2a+u),  u = √α - a
```
The exact error is:
```
zₙ - r = 2r λⁿ T(z₀) / (1 - λⁿ T(z₀))
```
Since |λ| < 1, the error decays geometrically with ratio |λ| =
(r−a)/(r+a). One application of the corrected-estimator construction has
cubic error in the bracket width w (Paper 1, Theorem 3); but iterating φ
*at a fixed anchor* is geometric, not cubic — the cubic rate in Paper 1
belongs to the zoom, which re-anchors at every level. For pronic α one has
u = r−a < 1/2 (u = 1/(√(1+1/a)+1), rising to 1/2 from below), so
|λ| = u/(2a+u) < 1/(4a): about log₁₀(4a) correct digits per iterate, exact
for every n via the closed form.

### 3.2 Exact Error Certification

The conjugacy provides an exact error identity:
```
|zₙ - r| = |2r λⁿ T(z₀)| / |1 - λⁿ T(z₀)|
```
The quantities r, λ, T(z₀) live in Q(√α), not Q, so the identity itself is
not a rational computation — but certified rational *bounds* are:

- r < a+1, so 2r < 2(a+1);
- |λ| = u/(2a+u) with u = r−a ∈ (0,1), so |λ| < 1/(2a) for a ≥ 2;
- |T(z₀)| < 1 since z₀ = a+1 > r.

Combining:

```
|zₙ − r| < 2(a+1)·(2a)⁻ⁿ / (1 − (2a)⁻ⁿ),
```

computable in exact rational arithmetic for every n. Stronger still, λ < 0,
so the iterates alternate strictly about r (Paper 2, Corollary 4): zₙ and
zₙ₊₁ bracket r, and |zₙ₊₁ − zₙ| = |zₙ − r| + |zₙ₊₁ − r| is a certified
two-sided error bound needing no bound on r at all.

### 3.3 Exact CF Convergents from Möbius Iterates

The CF convergents are recovered from the Möbius iterates:
```
p_{2k+1} = B_{2k+1},    q_{2k+1} = A_{2k+1} / a(a+1)
p_{2k}   = a·B_{2k},    q_{2k}   = A_{2k} / (a+1)
```

---

## 4. Why Only Pronic Numbers?

The exact correspondence fails for non-pronic α = a²+p with p|2a, p ≠ a.

**Proposition.** The exact denominator match B_{2k+1} = p_{2k+1} holds iff p = a.

*Proof sketch.* The proof of Theorem 1 (§2.5) uses α = a(a+1) at exactly one
point: the step φ(αq/p) = (a+1)(p+aq)/(p+(a+1)q) needs a²p + ap = a(a+1)p = αp,
which holds only when p = a. For α = a²+p with p ≠ a the step becomes
φ(αq/p) = [a(a+1)p + aαq]/(ap + αq) and the induction collapses. (The residue
law itself survives — for p | 2a the d-sequence gives Dₙ = +1 odd, −p even —
but the odd-n match Bₙ = pₙ then still fails, as the verification shows.)
Computational verification confirms the match fails for all tested (a,p) with
p|2a, p ≠ a, a ≤ 20.
∎

---

## 5. Algorithmic Implications

The correspondence yields a certified root-finding algorithm for pronic
numbers (geometric convergence, ratio |λ| = (r−a)/(r+a) < 1/(2a); certified
two-sided bracket from the alternation of §3.2):

```python
from fractions import Fraction as F

def certified_pronic_sqrt(a, target_digits, max_steps=1000):
    """Return (lo, hi, n): a certified rational bracket for sqrt(a(a+1)).

    The iterates alternate strictly about r (lambda < 0, Paper 2 Cor 4),
    so consecutive iterates bracket r and their difference is a certified
    two-sided error bound (Sec 3.2) -- no bound on r itself is needed.
    """
    if a < 2:
        raise ValueError("a must be >= 2")
    alpha = a * (a + 1)
    z_prev = F(a + 1)
    tol = F(1, 10 ** target_digits)
    for n in range(1, max_steps + 1):
        z = F(a) + F(a) / (F(a) + z_prev)
        if abs(z - z_prev) < tol:
            return (min(z, z_prev), max(z, z_prev), n)
        z_prev = z
    raise RuntimeError("no certified bracket within max_steps")
```

The algorithm produces exact rational brackets with geometric convergence
and a certified two-sided error bound at every step.

---

## 6. Conclusion

The Möbius-CF exact correspondence for pronic numbers reveals an algebraic unity: the geometrically convergent Möbius iteration (ratio u/(2a+u) < 1/(4a), exact closed form for every iterate) and the continued-fraction convergents (error O(1/qₙ²)) are two faces of the same algebraic structure in the quadratic field Q(√(a(a+1))). The Möbius conjugacy provides exact certified error bounds that the CF lacks, while the CF gives the optimal rational approximations.

This identity is the foundation for a new class of certified arbitrary-precision root-finding algorithms for pronic numbers, with potential extensions via Möbius transformations to other quadratic irrationals.

---

## Appendix: Verification Code

```python
from fractions import Fraction as F

def verify_pronic_theorem(a_max=100, n_max=15):
    """Verify Theorem 1 for all pronic α = a(a+1) up to a_max."""
    for a in range(2, a_max+1):
        alpha = a*(a+1)
        
        # CF convergents
        p = [a, 2*a+1]
        q = [1, 2]
        for n in range(2, n_max+1):
            an = 2 if n%2==1 else 2*a
            p.append(an*p[-1] + p[-2])
            q.append(an*q[-1] + q[-2])
        
        z = F(a+1)
        for n in range(n_max+1):
            A, B = z.numerator, z.denominator
            if n % 2 == 1:  # odd
                assert B == p[n], f"Odd n={n}: B={B}, p[n]={p[n]}"
                assert A == a*(a+1)*q[n], f"Odd n={n}: A={A}, α·q={a*(a+1)*q[n]}"
            else:  # even
                assert a*B == p[n], f"Even n={n}: a*B={a*B}, p[n]={p[n]}"
                assert A == (a+1)*q[n], f"Even n={n}: A={A}, (a+1)·q={(a+1)*q[n]}"
            z = F(a) + F(a) / (F(a) + z)
    print(f"All theorems verified for a ≤ {a_max}, n ≤ {n_max}")

if __name__ == "__main__":
    verify_pronic_theorem(100, 15)
    print("All 100 pronic numbers verified up to 15 iterations.")
```
