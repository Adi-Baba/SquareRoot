# Exact Dynamics of a Bracketed Zoom Method for Square Roots

**The Möbius structure, a closed form for every iterate, and the sharp
attainment of the classical error bound**

---

## Abstract

We analyse a recursive method for computing $\sqrt{\alpha}$ that begins from a
single observation: if $x^2<\alpha<(x+1)^2$, the position of $\alpha$ between
two consecutive perfect squares encodes an approximation to $\sqrt\alpha$,
and refining the same construction inside the resulting interval yields a
bracket that never breaks and contracts cubically per iteration. The method is
usually presented asymptotically. We show that no asymptotic treatment is
needed. The estimator family is a **Möbius transformation**; conjugating it by
$T(z)=\frac{z-r}{z+r}$ turns it into a pure scaling with multiplier
$\lambda=-u/(2a+u)$, and we obtain a closed form for the $n$-th iterate that
holds exactly for every $n\ge0$:

$$\varphi^{\,n}(z)-r=\frac{2r\,\lambda^{n}\,T(z)}{1-\lambda^{n}\,T(z)} .$$

Three results follow as corollaries rather than as separate arguments. First,
one-sidedness is exact for *every* iterate — the sequence provably alternates
about the root, and the corrected estimator is simply the even subsequence,
i.e. $\varphi^2$. Second, the classical sharp bound $w^3/(27a^2)$, together
with its constant $4/27$ and its maximiser $u^\star=2w/3$, is recovered as the
numerator of the closed form; and we make the word *sharp* precise by giving
the exact attainment factor
$\frac{18a^2}{18a^2+15aw+2w^2}$, which shows the bound is attained only in
the limit $w/a\to0$. Third, an existing $u$-aware certificate is recast in
these coordinates, where the interval $[\ell,h]$ of Theorem 8 is seen to arise
from a single unknown scalar $\lambda$ rather than from an interval of
positions. We also restate the general-$f$ interpolation rate, which settles a
conjecture left open in the source manuscript and explains two of its
unexplained numerical rows. We then sharpen what "sharp" means: an exact
factorisation (§5.4, Corollary 7) shows the constant $1/27$ is attained at no
bracket whatsoever, and that the sharp *position-blind* constant is
$2/(27(2+t))$, which the companion audit's baseline was already computing. A
final section closes the obvious follow-up: the optimal *position-aware*
certificate turns out to buy nothing over the one we already certify (§5.5). Finally we locate the boundary of the closed form.
For the $n$-th root the error polynomial is the complete homogeneous symmetric
polynomial $h_{n-2}(z,r,a)$, which is constant exactly when $n=2$; hence the
map has degree $n-1$ and **no** rational change of coordinates can linearise it
for $n\ge3$. The obstruction is not merely structural: the multiplier satisfies
$\lambda_n\to-(n-1)$ as the anchor weakens, so for $n\ge3$ the scalar iteration
is *repelling* unless the root is already known to within a relative
$\sim1.59/n$ -- a degradation with no counterpart at $n=2$, where the map
contracts from every bracket. All identities are confirmed symbolically and
re-checked in exact rational arithmetic; the accompanying test suite contains
no sampling of the kind the paper elsewhere withdraws as evidence.

**Keywords:** square root, Möbius transformation, functional iteration,
certified rounding, sharp error bound, interval refinement, interpolation
remainder.

---

## 1. Introduction

A family of elementary algorithms for $\sqrt\alpha$ proceeds by bracketing: if
$x^2<\alpha<(x+1)^2$ then $x<\sqrt\alpha<x+1$, and the residual
$p=\alpha-x^2$ measures how far $\alpha$ sits between the two squares. Writing
$r=\sqrt\alpha$ and $u=r-x$ for the offset of the root inside the bracket,
one has the exact factorisation

$$p=\alpha-x^2=(r-x)(r+x)=u(2x+u),\qquad r=x+\frac{p}{x+r},$$

so $r$ is a fixed point of the map $z\mapsto x+\frac{p}{x+z}$. This is the
starting point of the method studied here, and it has been developed in the
source manuscript [S1] with a one-sided error bound, a divergence condition,
and a recursion law. We take that construction as given and ask a question
that the asymptotic treatment leaves open:

> Since the iteration is a rational function of the current iterate, and the
> iterate enters only through one quotient, must it not be *exactly*
> solvable, rather than merely asymptotically analysable?

It must. The map is a linear fractional transformation, and the answer is
short. The purpose of this paper is to write that argument down cleanly,
because it changes what several later claims are entitled to say: the sharp
constant stops being a correction to someone else's derivation and becomes a
corollary of our own; the corrected estimator stops being a heuristic and
becomes the second iterate of a single map; and the one-sidedness guarantee
extends from two specific iterates to all of them.

### 1.1 Contributions

1. **Theorem 2 (conjugacy).** The estimator family is Möbius; $T$ conjugates it
   to a scaling with $\lambda=-u/(2a+u)$.
2. **Theorem 3 (closed form).** $\varphi^n(z)-r=\frac{2r\lambda^nT(z)}{1-\lambda^nT(z)}$
   exactly, for all $n\ge0$.
3. **Corollary 4 (alternation).** One-sidedness holds for every iterate; odd
   iterates of $z>r$ under-estimate, even ones over-estimate, and the contraction
   ratio is exactly $u/(2a+u)$.
4. **Corollary 5 (the corrected estimator).** The corrected estimate is
   $\varphi^2$, and its $k$-fold iteration is $\varphi^{2k}$.
5. **Corollary 6 (sharpness made precise).** The classical bound is recovered,
   with the exact attainment factor $18a^2/(18a^2+15aw+2w^2)$.
6. **Corollary 7 (strict slack).** $1/27$ is attained at no bracket: for every
   $t=w/a>0$ and $\eta=u/w\in[0,1]$ the slack against $\tfrac{2}{27(2+t)}$ is
   given by an explicit strictly-positive sum of two non-negative terms with
   disjoint zero sets. This replaces Corollary 6's asymptotic attainment factor
   with an exact identity, and identifies $2/(27(2+t))$ — not $1/27$ — as the
   sharp position-blind constant, which the §10.1 baseline already computes.
7. **Theorem 8 ($u$-aware certificate).** The position of the root inside the
   bracket yields a strictly tighter certificate at no extra cost. We also show
   (§5.5) that the *optimal* such certificate — the supremum
   $\max_{\eta\in[L,U]}F(t,\eta)$ of the exact error over the information
   possessed — beats Theorem 8 by a median of exactly zero and at most
   $5\times10^{-6}$, so the optimality is not worth pursuing.
8. **Theorem 9 (general $f$).** The per-decade gain of linear interpolation is
   $2$, $3$, or $1$ digits according as $f''(c)$ is nonzero, zero, or $f'(c)$
   degenerates; this settles an open conjecture of [S1].
9. **Theorem 10 (exact error for the $n$-th root).** The error polynomial is the
   complete homogeneous symmetric polynomial $h_{n-2}(z,r,a)$, constant only
   when $n=2$.
10. **Theorem 11 (the obstruction).** For $n\ge3$ no rational linearisation
    exists: $\deg(\Psi\circ\varphi_n)=e(n-1)\ne e=\deg(\lambda\Psi)$.
11. **Corollary 12 (alternation survives).** For every $n\ge2$ the sign of
    $r-\varphi_n(z)$ is the sign of $z-r$, so the iterates alternate — while
    *convergence* is a separate question and often fails (§9.4).
12. **Theorem 13 and its consequence.** $\lambda_n\to-(n-1)$ as $a\to0$, so the
    scalar $n$-th root iteration is repelling unless the root is already known to
    within a relative $\sim1.59/n$; $n=2$ is the only case that contracts
    unconditionally.

### 1.2 Relationship to prior work

**Numbering across the two papers.** [S2] and this paper are companion papers
read together, but each is numbered independently — monotonic in its own
document order, with no gaps and no duplicates. One result appears in *both*
and, after aligning the sequences, it carries the **same** number in both:

| result | this paper | [S2] |
|---|---|---|
| the $u$-aware bound (here §6) | **Theorem 8** | **Theorem 8** (§11) |

A second result appears in both under *different* numbers, because the papers
place it in different sections: the general-$f$ interpolation rate is Theorem 9
here (§7) and Theorem 7 in [S2] (§7). Cite the paper with the number in that
case.

The two papers also number *corollaries* by different conventions: [S2] gives a
corollary its parent theorem's number (Corollary 1 of Theorem 1, Corollary 2 of
Theorem 2), while this paper runs a single shared sequence in which corollaries
take their own numbers (Corollaries 4–7 all follow Theorem 3). Both sequences
are internally consistent — monotonic, no gaps, no duplicate numbers — so this
is a stylistic difference, not an error; it is recorded because "Corollary 7"
denotes different things in the two documents. Results unique to this paper
(Proposition 1, Theorems 2, 3 and 10–13, Corollaries 4–7 and 12) have no [S2]
counterpart.

None of the ingredients is claimed as new in isolation, and we do not claim
otherwise. The plain estimate $x+\frac{p}{x+z}$ is linear interpolation —
regula falsi — and its under/over-estimate structure is standard for that
family. The accelerated variant that iterates it with a decreasing bracket is
a known bracketing device [S1]. The general-$f$ analysis of Theorem 9 rests
on the classical interpolation remainder and a count of vanishing derivatives
at the root. What is new here is the observation that the *whole* family is a
single Möbius transformation and therefore exactly solvable, and the
consequences that follow from taking that seriously.

We note two scope boundaries. Section 9 concerns the *scalar* map of §1 for
arbitrary $n$; the companion-matrix and power-method literature addresses a
different construction in $n-1$ coordinates, and we do not treat it as a
remedy for the scalar method (see §9.5). And the method studied here is not
the fastest way to compute a square root in a working implementation, and we
make no such claim:
a companion audit [S2] benchmarks it against bisection, interval Newton and
Newton with a residual certificate, and reports that it loses to interval
Newton at 3000 digits. That is a statement about implementations and is
deliberately kept out of this paper, which is about the mathematics.

---

## 2. Notation and standing assumptions

| Symbol | Meaning |
|---|---|
| $\alpha$ | the radicand; $\alpha>0$ |
| $r$ | $\sqrt\alpha$ |
| $a$ | the bracket's left endpoint; $a>0$ |
| $b$ | the bracket's right endpoint, $b=a+w$ |
| $w$ | bracket width, $w>0$ |
| $u$ | $r-a$, the offset of the root inside the bracket; $0<u<w$ |
| $p$ | $\alpha-a^2=u(2a+u)$ |
| $\varphi$ | the plain estimator map, $z\mapsto a+\frac{p}{a+z}$ |
| $\lambda$ | the conjugate multiplier, $\frac{a-r}{a+r}=-\frac{u}{2a+u}$ |
| $T$ | the conjugacy, $z\mapsto\frac{z-r}{z+r}$ |

**Assumptions.** Throughout, $a>0$ and $0<u<w$, so $r=a+u>0$. The map $\varphi$
is non-degenerate: its representing matrix has determinant
$a^2-\alpha=-p$, which vanishes exactly when $a=r$. We assume $a\neq r$; the
case $a=r$ is the trivial fixed point. (The condition $a\neq r$ is what is
needed; the stronger hypothesis that $\alpha$ is not a perfect square,
adopted in an earlier version of this argument, is sufficient but not
necessary, and we correct it here.)

Notation $\varphi^n$ denotes the $n$-fold iterate, $\varphi^0=\mathrm{id}$.

---

## 3. The estimator family and its exact error

We first record the algebra the closed form will have to reproduce, so that
the later recovery is visible as a recovery rather than asserted.

**Proposition 1 (exact error).** For every admissible $z$ and every
$a>0,\ 0<u<w$,

$$M(a,z)-r=\frac{(r-a)(z-r)}{a+z},\qquad
\widehat r-r=\frac{u^2(w-u)}{D},\quad D=2a(2a+w)+u(2a+u),$$

where $M(a,z)=a+\frac{p}{a+z}$ is the plain estimate and
$\widehat r=M\bigl(a,M(a,b)\bigr)$ the corrected one.

*Proof.* Since $p=u(2a+u)$ and $r^2=a^2+p$,
$$M(a,z)-r=\frac{az+\alpha-r(z+a)}{z+a}=\frac{(a-r)(z-r)}{z+a}.$$
For the corrected estimate, substitute $z=M(a,b)$ in place of $b$ and combine;
the result is stated in [S1, Thm 3] and re-derived in §4.3 below from a
different direction, which is the point of the present argument. $\square$

---

## 4. The Möbius structure

### 4.1 The map is a linear fractional transformation

Since $a^2+p=\alpha$,
$$\varphi(z)=a+\frac{p}{a+z}=\frac{az+\alpha}{z+a},$$
represented by the matrix
$$M=\begin{pmatrix}a&\alpha\\1&a\end{pmatrix},\qquad \det M=a^2-\alpha=-p\neq0 .$$

Its fixed points solve $\varphi(z)=z$, i.e. $az+\alpha=z^2+az$, i.e.
$z^2=\alpha$. Hence exactly two: $r$ and $-r$. The second is spurious for our
purposes but is what makes the closed form possible.

### 4.2 Conjugacy to a scaling

**Theorem 2 (Möbius conjugacy).** With $T(z)=\dfrac{z-r}{z+r}$,
$$T\bigl(\varphi(z)\bigr)=\lambda\,T(z),\qquad
\lambda=\frac{a-r}{a+r}=-\frac{u}{2a+u},$$
for all $z\notin\{-r,-a\}$. In particular $|\lambda|<1$ for every $a>0$.

*Proof.* The structural reason is worth stating first, because it is what
makes $T$ the *right* conjugacy rather than a lucky guess: $T$ sends the two
fixed points to $0$ and $\infty$, and the only Möbius transformations fixing
$0$ and $\infty$ are pure scalings $w\mapsto\lambda w$. So $\varphi$ is
conjugate to *some* scaling, and computing the multiplier suffices.

Directly, using $\alpha=r^2$:
$$\varphi(z)-r=\frac{az+r^2-rz-ra}{z+a}=\frac{(a-r)(z-r)}{z+a},\qquad
\varphi(z)+r=\frac{(a+r)(z+r)}{z+a},$$
and taking the quotient gives the claim, since $a-r=-u$ and $a+r=2a+u$.
Finally $|\lambda|=\frac{u}{2a+u}<1$ iff $u<2a+u$ iff $a>0$. $\square$

### 4.3 The closed form

**Theorem 3 (exact iterate).** For every $n\ge0$ and $z\notin\{-r,-a\}$,
$$T\bigl(\varphi^{\,n}(z)\bigr)=\lambda^{n}T(z)
\quad\Longrightarrow\quad
\boxed{\ \varphi^{\,n}(z)-r=\frac{2r\,\lambda^{n}\,T(z)}{1-\lambda^{n}\,T(z)}\ } .$$

*Proof.* Induction on $n$ using Theorem 2. Since
$T^{-1}(v)=r\frac{1+v}{1-v}$ (substitute to check), writing
$v=\lambda^{n}T(z)$ gives $\varphi^n(z)=r\frac{1+v}{1-v}$, and

$$\varphi^n(z)-r=r\frac{(1+v)-(1-v)}{1-v}=\frac{2rv}{1-v}. \qquad\square$$

**Domain.** (The closed form) requires $z\neq-r$ and $|\lambda^nT(z)|<1$. For
$z>r>0$ we have $0<T(z)<1$ and $|\lambda|<1$, so both conditions hold for all
$n$: the iteration is globally well defined on the half-line $z>r$.

---

## 5. Consequences

### 5.1 One-sidedness for every iterate

**Corollary 4 (alternation).** For $z>r$, the iterates alternate strictly about
$r$: $\varphi^n(z)<r$ for odd $n$ and $\varphi^n(z)>r$ for even $n$. More
generally, for $z\neq r$ the sign of the error is
$$\operatorname{sign}\bigl(\varphi^n(z)-r\bigr)=\operatorname{sign}(\lambda^{n})\,\operatorname{sign}\bigl(T(z)\bigr),$$
and the signs reverse for $z<r$. The contraction ratio is exactly
$u/(2a+u)$, which equals $\varphi'(r)$.

*Proof.* $T(z)>0$ for $z>r$, the denominator of Theorem 3 is positive
(being $\ge 1-|T(z)|>0$ as $|\lambda^n|\le1$), and $\lambda<0$. For
$\varphi'(r)$: $\varphi'(z)=-p/(a+z)^2$ gives
$-u(2a+u)/(2a+u)^2=-u/(2a+u)$. $\square$

The classical under/over-estimate statement [S1, Cor 2] is the $n=1$ and $n=2$
case of this; the general statement was not previously available.

### 5.2 The corrected estimator is the even subsequence

**Corollary 5.** $m=\varphi(b)$ and $\widehat r=\varphi(m)=\varphi^2(b)$. The
$k$-fold corrected iteration is $\varphi^{2k}$, with exact error from
Theorem 3 at $n=2k$.

So the "correction" is not a new estimator: it is the second iterate of the map
already in hand. This is the paper's main interpretive claim — the corrected
estimate inherits its over-estimate property from $\lambda<0$, not from a
separate inequality.

### 5.3 Recovering the sharp bound, and what "sharp" means

Substituting $z=b=a+w$, $n=2$,
$$T(b)=\frac{w-u}{2a+w+u},\qquad \lambda^2=\frac{u^2}{(2a+u)^2}.$$
Theorem 3 then reads
$$\widehat r-r=\frac{2(a+u)\,u^2(w-u)\big/[(2a+u)^2(2a+w+u)]}
{1-u^2(w-u)\big/[(2a+u)^2(2a+w+u)]}.$$
The numerator of the denominator collapses:
$$(2a+u)^2(2a+w+u)-u^2(w-u)
=(2a+u)^3+u^3+w\bigl[(2a+u)^2-u^2\bigr]
=8a^3+12a^2u+6au^2+2u^3+4a^2w+4auw
=2(a+u)\bigl(4a^2+2au+u^2+2aw\bigr)
=2(a+u)\underbrace{\bigl[2a(2a+w)+u(2a+u)\bigr]}_{D}.$$
Cancelling gives exactly
$$\widehat r-r=\frac{u^2(w-u)}{D},$$
recovering Eq. (5) of [S1] from a different direction.

**Corollary 6 (sharpness, precisely).** Since $f(u)=u^2(w-u)$ satisfies
$f'(u)=u(2w-3u)$ with its unique interior critical point at $u^\star=2w/3$ and
$f(2w/3)=4w^3/27$, and $D\ge4a^2$,
$$0\le\widehat r-r\le\frac{4w^3/27}{4a^2}=\frac{w^3}{27a^2}.$$
The constant $4/27$ and the maximiser $2w/3$ are therefore not independent
discoveries: they are the numerator of Theorem 3.

We make the word *sharp* exact, because the bound is the product of two
*separately taken* suprema — $f$ at $u^\star=2w/3$ and $D$ as $u\to0$ — which
are not attained at the same $u$. The exact attainment factor is

$$\frac{f(2w/3)/D(2w/3)}{w^3/(27a^2)}=\frac{18a^2}{18a^2+15aw+2w^2}
=\frac{18a^2}{(3a+2w)(6a+w)}\ \xrightarrow[w/a\to0]{}\ 1 .$$

It equals $0.92$ at $w/a=0.1$, $0.51$ at $w/a=1$, and $0.22$ at $w/a=3$. Two
consequences, one a correction to how this result is usually stated:

* $u^\star=2w/3$ is the maximiser of the **numerator** $f$, not of the ratio
  $f(u)/D(u)$; the latter has a closed-form but unwielded maximiser, obtained
  by solving a quartic. The two coincide only in the limit $w/a\to0$. (Corollary 7
  below shows that quartic maximiser is $\tfrac23-\tfrac{t}{27}+O(t^2)$.)
* "The bound is sharp" is an **asymptotic** statement. We do not know of an
  earlier source that distinguishes this, and an earlier version of this
  argument did not.

Corollary 6 measures the failure to attain along one path, $u=2w/3$. The next
lemma settles it on *every* path at once, and upgrades the asymptotic to an
exact non-attainment statement.

### 5.4 Non-attainment: the constant is never reached

The attainment factor above is an asymptotic statement about a limit. We can
replace it with an exact identity, and the replacement is stronger than what the
limit alone gives.

Normalise the state by

$$t=\frac wa,\qquad \eta=\frac uw,\qquad 0<\eta<1 ,$$

in which case Proposition 1 reads

$$\boxed{\ \frac{e\,a^2}{w^3}=F(t,\eta)=\frac{\eta^2(1-\eta)}{4+2t+2t\eta+t^2\eta^2}\ }$$

and the published constant is recovered exactly as $F(0,\tfrac23)=\tfrac1{27}$.
That single evaluation locates the constant at a **corner** of the state space:
$t=0$ and $\eta=2/3$ are both excluded by $a>0$, $w>0$ and $0<u<w$.

**Corollary 7 (strict slack).** For every $t>0$ and every $\eta\in[0,1]$,

$$\frac{2}{27(2+t)}-F(t,\eta)=\frac{\underbrace{27(2+t)\bigl(\eta-\tfrac23\bigr)^2\bigl(\eta+\tfrac13\bigr)}_{\ge\,0,\;=0\iff\eta=2/3}+\underbrace{2t\eta(2+t\eta)}_{\ge\,0,\;=0\iff\eta=0}}{27(2+t)\bigl(4+2t+2t\eta+t^2\eta^2\bigr)}\;>\;0 .$$

*Proof.* Both numerator summands are non-negative on $[0,1]\times[0,\infty)$,
being a positive multiple of a square and a product of non-negative factors. The
first vanishes only at $\eta=2/3$; the second only at $\eta=0$ or $t=0$. Their
zero sets are therefore disjoint on $[0,1]\times(0,\infty)$, so their sum is
strictly positive there, and the denominator is strictly positive. $\square$

The identity is verified symbolically and re-checked in exact rational
arithmetic (`Src/mobius.py`, `normed_error`, `position_blind_constant`,
`strict_slack`; tests `test_corollary7_*`), which is the appropriate instrument per
§10 — it is an identity in $\mathbb{Q}$, not a measurement.

Three consequences.

1. **$1/27$ is never attained at any bracket.** For every $a>0$, $w>0$,
   $0<u<w$ one has strictly

   $$\sup_{\eta\in[0,1]}F(t,\eta)\;<\;\frac{2}{27(2+t)}\;<\;\frac1{27}.$$

   The published constant is a supremum over a boundary the recursion never
   reaches, not a value the error takes. Corollary 6's asymptotic statement is
   thus the weakest correct form of this fact; Corollary 7 is exact.

2. **The sharp position-blind constant is $2/(27(2+t))$, and the code already
   computes it.** Since $\tfrac{2}{27(2+t)}=\tfrac{1}{27}\bigl(1-\tfrac{t}{2+t}\bigr)$,
   a bound strictly better than $1/27$ is available with *no* position
   information whatsoever — by relative $t/(2+t)$, i.e. $33\%$ at $w/a=1$,
   $4.8\%$ at $w/a=1/10$, $0.05\%$ at $w/a=1/1000$. In unnormalised form this is

   $$\frac{4w^3/27}{D(0)}=\frac{2w^3}{27a(2a+w)},$$

   which is exactly what the certified benchmark's baseline returns. So $w^3/27a^2$
   was never the best constant the method could certify from position-blind
   information; §10.1's baseline was the sharper one all along, and neither paper
   recorded the gap.

3. **The "unwieldy quartic maximiser" of §8.3 is $2/3$ displaced by first order.**
   Expanding the stationary condition of $F$ about $\eta=2/3$,

   $$\eta^\star(t)=\frac23-\frac{t}{27}+\frac{t^2}{243}+O(t^3),\qquad
   \sup_{\eta}F(t,\eta)=\frac1{27}-\frac{5t}{162}+\frac{16t^2}{729}+O(t^3).$$

   So the true maximiser of the *ratio* and the maximiser $2w/3$ of the
   *numerator* agree to zeroth order and part company at $O(t)$; and the
   position-blind constant trails the true supremum by only $t/162$ at first
   order. The quartic is real but is a small perturbation of a known point, not
   an independent object.

A corollary of (2) worth one line, since it touches a stated guard in [S1]: the
progress condition $\varepsilon<w$ with $\varepsilon=\frac{2w^3}{27a(2a+w)}$ is
$2w^2<27a(2a+w)$, i.e. $w<\tfrac{3a}{4}\bigl(9+\sqrt{129}\bigr)\approx15.27\,a$
— wider than the $w<\sqrt{27}\,a\approx5.20\,a$ implied by the loose
$\varepsilon=\frac{w^3}{27a^2}$. So the sharper constant does relax the
bracket-shrinking guard.

It does not matter in practice, and that is *decided* rather than asserted:
running the actual recursion over its 16 test radicands, the largest $w/a$ at any
visited state is $1$ at level 0 and at most $10^{-2}$ thereafter — against
thresholds of $5.20$ and $15.27$. The guard is not the binding constraint at any
state the method reaches. (`Src/test_mobius.py`, `TestProgressGuard`: the
sharpness of $\varepsilon$, both guard roots compared exactly, and the bound on
$w/a$ over states the recursion really produces, in rationals.)

**What this does not do.** Corollary 7 bounds the error with *less* information. It
therefore does not compete with Theorem 8, which uses *more*: an actual interval
$[\ell,h]\ni u$. Nor does sharpening the constant compensate for the §10.1
finding that the method loses to interval Newton on precision-weighted volume.
It is a correction to how the sharp constant should be described, not a new
source of speed.

### 5.5 The optimal position-aware certificate, and why it buys nothing

The normalisation above poses the sharpest version of the position-aware
question. If the computation establishes $\eta=u/w\in[L,U]$, then the *optimal*
certificate is not an inequality but a supremum:

$$e\;\le\;C(t;L,U)\frac{w^3}{a^2},\qquad C(t;L,U)=\max_{\eta\in[L,U]}F(t,\eta),$$

the supremum of the exact error over precisely the information possessed.
Theorem 8 is an upper bound on that supremum; $C$ is it. It is worth asking
whether closing that gap is worth a paper. It is not, for a structural reason
rather than a numerical one. Theorem 8 is already tight on
the intervals the method produces: the bracket it hands the certificate is
$[\ell,h]=[m-a,\widehat r-a]$, whose width $h-\ell=\widehat r-m$ *is the current
error*. So the interval over which $C$ maximises is the quantity being bounded.
Refining the maximisation therefore buys resolution far below the quantity's own
size.

*Measured* (not proved — this is a measurement of two implementations, quoted
as ratios only, `Src/exp_u_aware.py` section (2b)): with $C$ evaluated by 64-fold
interval refinement against Theorem 8's closed form on the same interval and
through the same routine, the ratio $C/B$ is **1.0000 median and 1.0000 worst**
at $w\le a/100$, and $0.9979$ median at the far wider $w\le a$. No search is
involved: refining the piece count only ever returns a value of the bound at an
actual point, so it cannot under-shoot.

The instructive part is that a count of improved states would have called this a
result. In exact rationals, 95 of 120 visited states *do* improve — but by at
most $5\times10^{-6}$ relative, median exactly zero. The improvement rounds away
entirely in floating point, which is why the printed ratio reads $1.0000$. The
correct quantity is the magnitude, not the frequency
(`Src/test_u_aware.py`, `TestOptimalCertificate`, which asserts on the
magnitude and fails if it exceeds $10^{-4}$).

The consequence is that the "optimal certification" programme has nothing to
recover, and $1/27$ — corrected by Corollary 7 to $2/(27(2+t))$ — is the right
object to keep describing. What remains genuinely open in this direction is not
the constant but the *cost*: Theorem 8 buys its gain from two numbers the
iteration already computed, whereas recovering optimality would buy at most
$10^{-6}$ of a bound from the same two numbers.

---

## 6. A $u$-aware certificate

Because the exact error depends on $u$ alone, and $u$ is *not* known but
bracketed, the certificate can be tightened using only what the iteration
already computes.

**Theorem 8 ($u$-aware bound).** Let $a>0$, $0<u<w$, and let $[\ell,h]\subseteq[0,w]$
be any interval with $u\in[\ell,h]$ — in particular, the interval the method
produces for free, $\ell=m-a$, $h=\widehat r-a$. Then

$$\widehat r-r\ \le\ B(a,w;\ell,h):=\frac{\displaystyle\max_{v\in[\ell,h]}v^2(w-v)}
{\displaystyle 2a(2a+w)+\ell(2a+\ell)}\ \le\ \frac{w^3}{27a^2},$$

and $B$ is a valid error bound. The maximisation of the numerator is
elementary — its only critical point is $2w/3$ — so $B$ is a closed form with
no search. This is a consequence of Proposition 1, not an independent
construction: it uses that $D(u)=2a(2a+w)+u(2a+u)$ is strictly increasing
($D'(u)=2a+2u>0$), hence minimised at $\ell$.

**Remark (the conjugacy restated).** Theorem 3 needs only the *scalar*
$\lambda=-u/(2a+u)$, not the position $u$. A tighter certificate is therefore
equivalent to a tighter interval for $\lambda$, and the interval $[\ell,h]$
above is exactly the statement that $\lambda\in[\lambda(h),\lambda(\ell)]$
since $\lambda$ is a decreasing bijection of $(0,w)$ onto $(-1,0)$. In these
coordinates the "$u$-aware bound" is not an interval-of-positions device at
all: it is ordinary interval arithmetic on one scalar.

---

## 7. General root functions

**Theorem 9 (interpolation rate).** Let $f\in C^3([a,b])$ with a simple root
$c\in(a,b)$, $0<|f'(c)|<\infty$, and set $A=c-a$, $B=b-c$. The linear
interpolation error satisfies
$$|c-x_L|=\frac{A B}{2}\,\Bigl|\frac{f''(\xi)}{f'(c)}\Bigr|$$
for some $\xi\in(a,b)$, whence the per-decade gain of digits is

| condition | leading behaviour | digits per decade |
|---|---|---|
| $f''(c)\neq0$ | $A,B=O(w)$ | **2** |
| $f''(c)=0$ | $f''(\xi)=f'''(c)(\xi-c)=O(w)$ | **3** |
| $|f'(c)|\to0$ or $\to\infty$ | $w^2$ term dies | **1** |

*Proof sketch.* The interpolation remainder gives the identity directly; the
three rows are the three ways the factor $|f''(\xi)/f'(c)|$ can scale with $w$.
$\square$

This settles a conjecture left open in [S1] — that $\pi$ enjoys the
$3$-digit regime because $\sin''(\pi)=0$ — and it supplies the reason for two
rows of that manuscript's table which it could not explain.

---

## 8. Discussion

### 8.1 What the closed form buys

* **The recursion is not heuristic.** Theorem 3 replaces asymptotic reasoning
  with arithmetic; the classical bound is recovered as a corollary.
* **The "correction" is the even iterate.** Corollary 5 explains one-sidedness
  by the sign of $\lambda$ rather than by a separate inequality.
* **The certificate is one-dimensional.** §6 shows the position interval is
  equivalent to an interval on a single scalar.

### 8.2 What it does not buy

* **No speed.** The method remains slower than interval Newton at 3000 digits
  in the implementation audited in [S2]. This paper makes no performance
  claim; the closed form is a statement about mathematics, not about runtime.
* **No better order.** The closed form is *exact*, not higher-order. The
  natural next question is whether a rational map of $T$ squaring $\lambda$ —
  the degree-$4$ analogue, which is Halley's method in these coordinates —
  improves anything. Since Theorem 3 shows the plain iteration already
  contracts geometrically with a known factor, such a variant must earn its
  keep in operations per step, and [S2] measures that budget as the binding
  constraint.

### 8.3 Limitations

* $a\neq r$ is required for non-degeneracy; the trivial case $a=r$ is
  excluded.
* The exact maximiser of $f(u)/D(u)$ is a quartic root and is reported only as
  an attainment factor (§5.3), not in closed radical form. §5.4 locates it to
  $\tfrac23-\tfrac{t}{27}+\tfrac{t^2}{243}+O(t^3)$, so the residual uncertainty
  is $O(t)$ about a known point rather than an unlocated root.
* Theorem 8 is stated for the bracket produced by the iteration; its
  generalisation to arbitrary $f$ is open, and would need $D'(u)>0$-type
  monotonicity, which is an artefact of $x^2-\alpha$ and is not expected to
  survive.
* The $0<\alpha<1$ regime requires rescaling and is not analysed here. In
  exact arithmetic the recursion provably cannot start from $a=0$: the
  corrected estimator degenerates to $\widehat r=1$ for every such $\alpha$, and
  the progress condition becomes $w<0$.
* §9 treats the scalar map only. It says nothing about the cost or accuracy of
  any $n$-th-root algorithm, and in particular nothing about the
  higher-order Newton variants that §9.4 recommends on convergence grounds
  alone.

---

## 9. Where the structure stops: the $n$-th root

The closed form of §4–§5 is not an artefact of a particular proof strategy.
It is available for $n=2$ and for no other $n$, and this section locates the
boundary.

### 9.1 The $n$-th root map and its exact error

For $\alpha>0$, $r=\alpha^{1/n}$, and $0<a<r$, put
$$S_n(z,a)=\sum_{k=0}^{n-1}z^{\,n-1-k}a^k,\qquad
\varphi_n(z)=a+\frac{r^n-a^n}{S_n(z,a)} .$$
This is the direct analogue of §1, since $r=a+\frac{r^n-a^n}{S_n(r,a)}$.
As a rational function $\varphi_n$ has degree $n-1$, and its fixed points are
the $n$ numbers $r\zeta^j$ with $\zeta=e^{2\pi i/n}$: indeed
$\varphi_n(z)=z\iff(z-a)S_n(z,a)=r^n-a^n\iff z^n=r^n$.

**Theorem 10 (exact error, general $n$).** For $0<a<r$,
$$r-\varphi_n(z)=\frac{(r-a)(z-r)\,Q_n(z,r,a)}{S_n(z,a)},\qquad
Q_n(z,r,a)=\sum_{k=0}^{n-2}\sum_{j=0}^{n-2-k}z^{\,n-2-k-j}r^{j}a^{k},$$
and $Q_n$ is the **complete homogeneous symmetric polynomial**
$$Q_n(z,r,a)=h_{n-2}(z,r,a):=\!\!\sum_{i+j+k=n-2}\!\!z^{i}r^{j}a^{k}.$$

*Proof.* $S_n(z,a)-S_n(r,a)=\sum_{k=0}^{n-2}a^k\bigl(z^{n-1-k}-r^{n-1-k}\bigr)$,
and $z^{m}-r^{m}=(z-r)\sum_{j=0}^{m-1}z^{m-1-j}r^{j}$; the inner exponents
satisfy $(n-2-k-j)+j+k=n-2$, and the surviving monomials are exactly all
triples $(i,j,k)$ of non-negative integers summing to $n-2$. $\square$

For $n=2$ this gives $Q_2=1$ and recovers Prop 1 verbatim. For $n=3$,
$Q_3=z+r+a$.

### 9.2 The obstruction

**Theorem 11 (no rational linearisation).** For $n\ge3$ there is no rational
$\Psi$ of degree $e\ge1$ and no $\lambda$ with
$\Psi\circ\varphi_n=\lambda\,\Psi$ as rational maps.

*Proof.* Degrees of rational maps multiply under composition:
$$e(n-1)=\deg(\Psi\circ\varphi_n)=\deg(\lambda\Psi)=e .$$
Since $e\ge1$ this forces $n=2$. $\square$

The point of the symmetric-polynomial description is that it makes §9.1 and
§9.2 *the same fact seen twice*. A Möbius map has degree $1$, and it is
precisely the case $n-1=1$ — that is, $Q_2$ constant — in which the error
$\frac{(r-a)(z-r)Q_n}{S_n}$ is a ratio of linear forms. For every $n\ge3$,
$Q_n$ has degree $n-2\ge1$ and the obstruction is absolute rather than
methodological: no change of rational coordinates can help.

### 9.3 One-sidedness survives

**Corollary 12.** For $a,r,z>0$ and every $n\ge2$,
$$\operatorname{sign}\bigl(r-\varphi_n(z)\bigr)=\operatorname{sign}(z-r),$$
because $(r-a)$, $Q_n$ and $S_n$ are all strictly positive. Hence each step
moves to the opposite side of $r$, and the iterates alternate for **all** $n$ —
not only for $n=2$.

Note carefully what does and does not follow. Alternation is an algebraic fact
about signs and holds unconditionally; *convergence* is a separate question,
and for $n\ge3$ the answer is often no (§9.4). Alternation is not a convergence
guarantee.

### 9.4 The multiplier, and a $1/n$ degradation

**Theorem 13.** With $t=a/r\in(0,1)$,
$$\lambda_n:=\varphi_n'(r)=-\frac{(r-a)S_n'(r,a)}{S_n(r,a)}
=1-\frac{n}{\sum_{k=0}^{n-1}t^{k}} .$$
In particular $\lambda_n\to-(n-1)$ as $a\to0^{+}$, and $\lambda_2$ is the only
case for which $|\lambda_n|<1$ for every bracket.

*Proof.* From $\varphi_n'(z)=-p_nS_n'(z,a)/S_n(z,a)^2$ and
$S_n(r,a)=p_n/(r-a)$; for the second form, $S_n(z,a)=\frac{z^n-a^n}{z-a}$
gives $S_n'(r,a)=\frac{nr^{n-1}(r-a)-p_n}{(r-a)^2}$, and substitution yields
$\lambda_n=1-\frac{nr^{n-1}}{\sum_{k<n}r^{n-1-k}a^k}$. $\square$

For $n=2$, $|\lambda_2|=\frac{u}{2a+u}<1$ always: the square-root map
contracts from *every* bracket, which is why the method of [S1] needs no
guard on that point. For $n\ge3$ the fixed point is **repelling** once the
anchor is too weak. Writing $\lambda_n=-1$ gives the boundary
$\sum_{k=0}^{n-1}t^k=n/2$, and the required relative accuracy $1-t$ before the
iteration contracts at all is

| $n$ | 2 | 3 | 4 | 10 | 100 | 500 |
|---|---|---|---|----|-----|-----|
| required $u/r$ | none | $0.634$ | $0.456$ | $0.168$ | $0.0160$ | $0.0032$ |

and asymptotically $1-t_{\text{crit}}\sim 1.5936/n$, from the solution of
$1-e^{-y}=y/2$. This is a usable obstruction rather than a curiosity: it
predicts that a fixed-anchor regula-falsi scheme is unusable for high-order
roots, and that the appropriate family is the higher-order Newton variants
(Halley, Householder), whose attracting regions do not shrink with $n$.

### 9.5 What escapes, and what we do not claim

Locally, Koenigs' theorem supplies a conjugacy $\Psi_n$ with
$\Psi_n(\varphi_n(z))=\lambda_n\Psi_n(z)$ — but only where
$0<|\lambda_n|<1$, i.e. only inside the attracting region of §9.4, and
$\Psi_n$ is generally a transcendental power series rather than a rational
function. Iterating it gives an asymptotic law
$\varphi_n^m(z)-r\sim\Psi_n(z)\lambda_n^m$, which is the correct replacement
for the exact formula of §4.3, but a strictly weaker statement.

Two things are deliberately **not** claimed here.

First, this section is about the *scalar* map of §1. A different, classical
construction linearises the companion-matrix pencil $aI+C_n$ by eigenvector
decomposition and contracts for every $a>0$ and every $n$ — the power method
for eigenvalues. That is a different algorithm in $n-1$ coordinates, not a
remedy for this one: it does not carry the bracket, and its convergence is
unrelated to $\lambda_n$. Observing that it contracts precisely where the
scalar map is repelling says nothing about rescuing the scalar method.

Second, no performance claim is made anywhere in this section. §9.4 predicts
that higher-order methods are the right choice for $n\ge3$; it does not
measure them.

---

## 10. Reproducibility

Every identity in §4–§6 is confirmed symbolically and re-checked in exact
rational arithmetic. There is no sampling anywhere in the verification, and in
particular no claim of the form "no violation observed in $N$ random trials" —
that form of evidence is withdrawn in [S2] for the reason that it describes a
sample rather than a universally quantified statement.

| item | what it checks |
|---|---|
| `Src/mobius.py` | the algebraic definitions and the closed form, in `Fraction` |
| `Src/mobius_n.py` | the $n$-th root map, the symmetric polynomial $Q_n$, the multiplier and the attracting threshold, in `Fraction` |
| `Src/test_mobius_n.py` | Theorem 10 exactly for $n=2..5$; the $n=2$ recovery; the sign flip of Cor 12; the repelling region of Thm 13 |
| `Src/test_mobius.py` | closed form vs actual iteration, $n=0..6$; alternation; the attainment factor; the $\alpha<1$ degeneracy; Corollary 7 (identity, strictness, the $2/(27(2+t))$ constant, and its equality with `exp_u_aware.bound_published`) |
| `Src/test_u_aware.py` | Theorem 8's validity and same-footing; the §5.5 null result — `C` never under-shoots, never exceeds Theorem 8, and the relative gain over visited states stays under $10^{-4}$ |
| `Src/exact_decide.py` | decides cubic inequalities on rational intervals exactly, with witnesses |
| `Src/verify_theorems.py` | the proved / decided / measured ledger |

For §9 the corresponding claims are decided rather than illustrated: the
identity of Theorem 10 is checked at exact rational points, the sign flip of
Corollary 12 follows from the sign of the three positive factors rather than
from iterating, and the repelling region of Theorem 13 is decided by
`attracting(a, r, n)`, which is an exact comparison in $\mathbb{Q}$. The
asymptotic constant $1.5936$ is used only for the table in §9.4, which is
illustrative; every claim in that subsection other than the table's stated
asymptotics is exact.

The symbolic identities used for §4-§6 are recorded in `Src/mobius.py` as
(S1)–(S9): $T\circ\varphi=\lambda T$; the closed form at $n=2$ reproducing
Eq. (5); the numerator maximised at $2w/3$ with value $4w^3/27$;
$\varphi'(r)=\lambda$; the attainment factor of Corollary 6; the normalised
error $e\,a^2/w^3=F(t,\eta)$ and $F(0,\tfrac23)=\tfrac1{27}$; the strict-slack
identity of Corollary 7; and the expansions
$\eta^\star(t)=\tfrac23-\tfrac{t}{27}+\tfrac{t^2}{243}+O(t^3)$ and
$\sup_\eta F=\tfrac1{27}-\tfrac{5t}{162}+O(t^2)$.

Corollary 7 is the only claim of §5.4 whose proof is a *disjointness* argument rather
than a one-line identity, so its test does not merely re-evaluate the formula: it
asserts that each summand is non-negative, that the sum is strictly positive for
$t>0$, and that the two zero sets ($\eta=2/3$ and $\eta=0$) are distinct. That is
the content of the proof, and it is what the tests check. Eight single-edit
faults injected into the three supporting functions — a wrong factor
$\eta+\tfrac13\to\eta+\tfrac12$, the mixed term dropped, the two summands
swapped, $27\to26$, the constant reverted to $1/27$, and three slips inside $F$
itself — are all caught; a suite that only checked the formula numerically would
have passed at least the disjointness claim.

---

## References

- **[S1]** The source manuscript whose claims are corrected and extended here;
  see the companion audit [S2] for its full statement and a list of its
  specific errors.
- **[S2]** Companion audit of [S1]: certified benchmarking against bisection,
  interval Newton and Newton with a residual certificate, and the correction
  of the manuscript's performance hypotheses.

*Bibliography note.* The references above are placeholders identifying
internal documents. A submission-ready version requires a checked external
bibliography — the classical anchors for this material are the standard
treatments of linear fractional transformations and of the interpolation
remainder (e.g. the relevant chapters of standard texts on numerical analysis
and on approximation theory) and the classical papers on higher-order
root-finders. **No external citation has been verified here, and none should
be copied from this draft into a submission without checking.**
