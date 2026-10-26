# A Bracketed Zoom Method for Square Roots

### Exact error identities, a sharp cubic recursion, and the interpolation rate for general root functions

---

## Abstract

We develop a recursive method for computing $\sqrt{\alpha}$ that begins from a
single observation: if $x^2 < \alpha < (x+1)^2$, then $\alpha$ occupies a
specific position between two consecutive perfect squares, and that position
encodes an approximation to $\sqrt{\alpha}$. Refining the same construction
inside the resulting interval yields a method whose bracket never breaks and
whose width contracts cubically per iteration.

The method rests on an exact algebraic identity. Writing $r=\sqrt{\alpha}$ and
$p=\alpha-a^2=(r-a)(r+a)$, one obtains

$$
r = a + \frac{p}{a+r},
$$

so that $r$ is the fixed point of $z\mapsto a+p/(a+z)$. Evaluating this map at
an arbitrary $z\in[a,b]$ gives an estimator $M(a,z)=a+p/(a+z)$ whose error is
*exact* and known in closed form:

$$
r - M(a,z) = \frac{(r-a)(z-r)}{a+z}. \tag{1}
$$

Identity (1) separates the method into two regimes. Taking $z=b>r$ yields a
strict under-estimate; feeding that under-estimate back as $z$ yields a strict
over-estimate whose error is cubic in the bracket width. A closed-form
analysis of this error gives a one-sided bound $0\le \widehat{r}-r\le w^3/(27a^2)$
for a bracket of width $w$, which is $27/16\approx 1.69$ times sharper than a
natural first derivation. The bound is *sharp*: the ratio of the true error to
the bound is $16/27$ in the limit $w/a\to 0$ and never exceeds it.

We further determine the precise condition under which the recursion makes
progress, $w<a\sqrt{c}$ for a recursion $w\mapsto w^3/(ca^2)$, and show that a
missing guard on this condition causes the iteration to diverge. Recursion
depth is $\Theta(\log D)$ with a per-level tripling of correct digits, which we
confirm to within one level of a previously reported measurement.

Because the exact error depends on the root's position $u$ inside the bracket,
and the algorithm already computes quantities bracketing $u$, we derive a
$u$-aware bound that is a median $0.27$–$0.62\times$ the published one
depending on regime while never exceeding it. Measured against a common
evaluation of both sides, it costs $1.00$–$1.29\times$ and saves zero to two
recursion levels — a tighter certification at constant cost, not an
asymptotic gain.

Finally, we give a general interpolation theorem that settles an open question
about the method's behaviour on non-square-root functions. For a simple root
$c$ of $f$ with $0<|f'(c)|<\infty$, the linear-interpolation error satisfies

$$
|c - x_L| = \frac{|f''(\xi)|}{2|f'(c)|}\,(x_L-a)(b-x_L), \tag{2}
$$

so the per-decade gain is $2$ digits when $f''(c)\neq 0$, $3$ digits when
$f''(c)=0$, and $1$ digit in the degenerate cases $f'(c)=0$ or $|f'(c)|=\infty$.
We verify all three predictions across nine test functions, and confirm the
previously untested conjecture that $\pi$ enjoys the 3-digit regime because
$\sin''(\pi)=0$.

The estimator family is in fact Möbius: writing $T(z)=\frac{z-r}{z+r}$
conjugates it to a scaling with $\lambda=-u/(2a+u)$, so the $n$-step iterate is
*exactly* $\frac{2r\lambda^{n}T(z)}{1-\lambda^{n}T(z)}$. This subsumes the exact
error identity, the sharp constant and the $2w/3$ maximiser as special cases,
and makes one-sidedness exact for every iterate rather than only the first two.
That argument, and the general-$f$ rate theorem, are developed in a companion
paper, `Paper2.md`; this paper concentrates on what can be measured.

Finally, we test the draft's own performance hypothesis and refute it. The
$\Theta(D\log D)$ cost of running every level at full target precision is
real, but per-level precision scheduling delivers $2.90\times$ at 3000 digits
and still leaves Newton $3.64\times$ faster — so the residual gap is a
constant factor in the algorithm, not a precision-scheduling artefact. The
investigation also surfaced a defect worth more than the hypothesis: the
initial bracket was computed at full target precision and then floored to an
integer, consuming 58% of every call.

Because a point estimate cannot test a bracket-carrying method, §10.1 then
re-runs the comparison against *certified* peers — bisection, interval Newton,
and Newton with a residual certificate — each checked by exact integer
arithmetic rather than a reference square root (320 of 320 runs certified).
That comparison refutes the claim we set out to make: on identical seeds the
method is $417\times$ faster than bisection but $2.61\times$ slower than
interval Newton at 3000 digits, taking fewer steps than that peer and still
losing. Building it exposed five problems, three of them defects — a
termination test looser than its stated tolerance by a factor of ten, a bracket
that certified before being printed and failed after, and a convergence test
that ran its own verifier inside the timed region — and two methodological: a
competitor under-scheduled badly enough to look absurd, and a seeding
asymmetry worth $1.66\times$ that alone accounts for the apparent gap to
Newton.

All 167 verification tests accompanying this work pass. The suite is
mutation-tested: 14 injected faults, 12 caught, the 2 survivors analysed in
§8.1 rather than papered over.
Several claims in the original draft are found to be incorrect; these are
documented in §9 with corrected statements, and the §10 performance prediction
is corrected in place against measurement.

---

## 1. Introduction

Consider $\sqrt{12}\approx 3.4641$, which lies between $3$ and $4$. We have
$3^2=9$, $4^2=16$, and

$$
12-9=3, \qquad 16-12=4, \qquad 16-9=7=3+4.
$$

The value $12$ therefore sits at a definite place between two consecutive
squares — three units above the lower one, four below the upper one — and that
ratio $3:4$ is a usable approximation to $\sqrt{12}$. Generalising: for any
$\alpha$ strictly between consecutive squares, the position of $\alpha$ between
them yields an estimate of $\sqrt{\alpha}$. Applying the same construction
recursively to progressively smaller intervals sharpens the estimate.

The formal content of this observation is surprisingly clean. Section 3 shows
that $\alpha$ is *exactly* a weighted average of the two bounding squares, and
that $\sqrt{\alpha}$ is *exactly* a weighted root-mean-square of the two
bounding integers. Section 4 gives an exact error identity for the resulting
estimator, and Sections 5–6 turn it into a rigorous recursion with a certified
bracket.

### 1.1 Contributions

1. **An exact error identity** (Theorem 3) for the estimator family
   $M(a,z)=a+p/(a+z)$, valid for arbitrary $z$ and exact in rational arithmetic.
2. **A one-sided error bound with a sharp constant** (Theorem 4):
   $0\le \widehat r - r \le w^3/(27a^2)$, optimal in the limit $w/a\to 0$.
   The companion theory paper [S2, Corollary 7] sharpens this: the inequality is
   *strict* for every bracket, and the sharp position-blind constant is
   $2w^3/(27a(2a+w))$, attained at no state.
3. **A closed-form error** $\widehat r - r = u^2(w-u)/(2a(2a+w)+u(2a+u))$
   (Theorem 3) that makes both the sharp constant and the $u=2w/3$ maximiser
   explicit.
4. **A progress condition** (Proposition 5) with an explicit counterexample
   showing that omitting it makes the recursion diverge.
5. **A recursion-depth law** (Theorem 6) giving exact digit counts per level.
6. **An interpolation-rate theorem for general $f$** (Theorem 7) that
   explains every row of the empirical general-function table, including the
   previously unexplained flat and steep cases, and settles the draft's
   outstanding conjecture about $\pi$.
7. **A $u$-aware error bound** (Theorem 8) that uses the position of the root
   inside the bracket — already computed for free by the method — to certify
   a bracket a median $0.27$–$0.62\times$ as wide as Theorem 4's, at
   $1.00$–$1.29\times$ cost.
8. **An exact closed form for the whole iteration** (proved in the companion
   paper `Paper2.md` §4.3; the consequence is carried here in §5.4): the
   estimator family is Möbius, so
   $T(z)=\frac{z-r}{z+r}$ conjugates it to a scaling and the $n$-step iterate
   is *exactly* $\frac{2r\lambda^{n}T(z)}{1-\lambda^{n}T(z)}$. This subsumes
   the exact-error identity, the sharp constant and the $2w/3$ maximiser as
   special cases, and makes one-sidedness exact for every iterate. It also
   forces a correction to how "sharp" is stated: the bound is the product of
   two separately taken suprema and is attained only as $w/a\to0$ — and
   `Paper2.md` Corollary 7 sharpens even that, showing the inequality is
   *strict* at every bracket and that $w^3/(27a^2)$ is attained nowhere.

   **Numbering across the two papers.** This paper and `Paper2.md` are read
   together but each is numbered independently, monotonic in its own document
   order with no gaps or duplicates. One result appears in both and — after
   aligning the two sequences — carries the **same** number in both:

   | result | this paper | `Paper2.md` |
   |---|---|---|
   | the $u$-aware bound (here §11) | **Theorem 8** | **Theorem 8** (§6) |

   A second appears in both under different numbers, the papers placing it in
   different sections: the general-$f$ interpolation rate is Theorem 7 here
   (§7) and Theorem 9 in `Paper2.md` (§7). Cite the paper with the number in
   that case.

   The two papers also number *corollaries* by different conventions: here a
   corollary carries its parent theorem's number (Corollary 1 of Theorem 1,
   Corollary 2 of Theorem 2), whereas `Paper2.md` runs a single shared sequence
   in which corollaries take their own numbers (Corollary 4–7 follow Theorem 3).
   Both sequences are internally consistent — monotonic, no gaps, no duplicate
   numbers — so this is a stylistic difference rather than an error, and it is
   recorded here because "Corollary 7" denotes different things in the two
   documents.

   Results unique to this paper (Theorems 1–6 and 8, Corollaries 1–2,
   Proposition 5) have no `Paper2.md` counterpart; results unique to
   `Paper2.md` (Proposition 1, Theorems 2, 3 and 10–13, Corollaries 4–7 and
   12) have none here.
9. **A certified comparison** (§10.1) against bisection, interval Newton, and
   Newton with a residual certificate, every bracket proved by exact integer
   arithmetic rather than a reference square root (320/320 runs). It reports a
   negative result: on identical seeds the method is $417\times$ faster than
   bisection and $2.61\times$ slower than interval Newton at 3000 digits,
   despite taking fewer steps. Building it exposed three defects — a termination
   test looser than its stated tolerance by a factor of ten, a bracket that
   certified before being printed and failed after, and a convergence test
   that ran its own verifier inside the timed region — and two methodological
   errors of our own, an under-scheduled competitor and a seeding asymmetry
   that alone accounts for the apparent gap to Newton.

### 1.2 Relationship to prior work

None of the individual ingredients is new, and we do not claim otherwise. The
$O(w^2)$ plain estimate is linear interpolation, i.e. regula falsi. The
one-sided under/over-estimate structure is standard for that family. The cubic
bracket recursion is a known accelerated bracketing device. Digit-by-digit
square root is the classical $10$-ary refinement. The contributions here are
(a) the unified exact-algebra presentation of §3, (b) the sharp bound of
Theorem 4 with its closed form, (c) the explicit progress condition of
Proposition 5, and (d) the general interpolation theorem of §7, which we
believe was not previously stated in this form.

---

## 2. Notation and standing assumptions

Let $\alpha>0$ be a real number that is not a perfect square, and let

$$
r=\sqrt{\alpha}.
$$

Throughout, a **bracket** is a pair $(a,b)$ with

$$
a < r < b, \qquad a \ge 1, \qquad w := b-a>0.
$$

The condition $a\ge 1$ is used only in the original bound of §5.2; §5.3 shows
the sharp bound requires merely $a>0$. We write

$$
p=\alpha-a^2=(r-a)(r+a),
$$

and note $p>0$ throughout. For the initial bracket we take $a=\lfloor r\rfloor$,
$b=a+1$; for integer $\alpha$ this is $x=\lfloor\sqrt\alpha\rfloor$,
$y=x+1$, and then

$$
p=\alpha-x^2\ \ (\text{deficit}), \qquad q=y^2-\alpha\ \ (\text{excess}), \qquad
d=y^2-x^2=p+q=2x+1.
$$

Because $d=w(a+b)$ when $w=1$, the expressions $a+\frac{p}{a+b}$ and
$a+\frac{wp}{d}$ of the draft coincide; we consistently use the first.

---

## 3. Exact algebraic structure

**Theorem 1 (weighted average).** $\alpha$ is a weighted average of $x^2$ and
$y^2$ with weights $q/d$ and $p/d$:

$$
\alpha = \frac{q\,x^2 + p\,y^2}{d}.
$$

*Proof.* Using $y^2=x^2+d$ and $p+q=d$:

$$
q x^2 + p y^2 = q x^2 + p(x^2+d) = (p+q)x^2 + pd = d(x^2+p)= d\alpha. \qquad\square
$$

**Theorem 2 (weighted RMS).** $r$ is the weighted root-mean-square of $x$ and
$y$ with the same weights:

$$
r = \sqrt{\left(x+\frac{p}{d}\right)^2 + \frac{pq}{d^2}}. \tag{3}
$$

*Proof.* The weighted mean is $\frac qd x+\frac pd y = x + \frac pd$ because
$y-x=1$. The weighted variance is $\frac qd\frac pd (y-x)^2 = \frac{pq}{d^2}$.
Since (weighted) mean-square $=$ mean$^2+$ variance, and the mean-square here
equals $\frac qd x^2+\frac pd y^2 = \alpha = r^2$ by Theorem 1, the claim
follows. $\square$

**Corollary 1 (first correction).** Let $m_0=x+\frac pd$,
$v=\frac{pq}{d^2}>0$, and $\widehat m_0=m_0+\frac{v}{2m_0}$. Then
$\widehat m_0>r$ always, and

$$
0<\widehat m_0-r = \frac{v^2}{4m_0^2\,(\widehat m_0+r)},
$$

so the residual error is *quadratic* in the variance $v$.

*Proof.* By Theorem 2, $r^2=m_0^2+v$ with $v>0$, so $r>m_0>0$ and
$\widehat m_0=m_0+\frac{v}{2m_0}>0$. Then
$\widehat m_0^{\,2}-r^2=\left(m_0+\frac{v}{2m_0}\right)^2-(m_0^2+v)
=\frac{v^2}{4m_0^2}>0$, so $\widehat m_0>r$ strictly; factoring
$\widehat m_0^{\,2}-r^2$ gives the stated expression. $\square$

*Remark.* The error of the level-0 linear estimate is therefore always
negative (it under-estimates), with magnitude governed by the weighted
variance $pq/d^2$. The estimate is worst when $p\approx q$, i.e. when $\alpha$
sits near the midpoint between the two squares, and this is visible in Table 1:
$\alpha=12$ has $p=3,q=4$ and the largest plain error, while $\alpha=101$ has
$p=20,q=21$ but a far smaller one because $d=201$ is large. The Corollary-3
correction reduces the error by two to three orders of magnitude in every case
tested:

| $\alpha$ | $m_0$   | $r-m_0$               | Corollary 1  | $\widehat m_0-r$     |
| ---------- | --------- | ----------------------- | ------------ | ---------------------- |
| 12         | 3.428571  | $+3.553\times10^{-2}$ | 3.464285714  | $+1.84\times10^{-4}$ |
| 10         | 3.142857  | $+1.942\times10^{-2}$ | 3.162337662  | $+6.00\times10^{-5}$ |
| 2          | 1.333333  | $+8.088\times10^{-2}$ | 1.416666667  | $+2.45\times10^{-3}$ |
| 17         | 4.111111  | $+1.199\times10^{-2}$ | 4.123123123  | $+1.75\times10^{-5}$ |
| 101        | 10.047619 | $+2.257\times10^{-3}$ | 10.049875875 | $+2.53\times10^{-7}$ |

*Remark (pronic case).* When $\alpha=x(x+1)$ is pronic, then $p=x$, $q=x+1$
and $x+p/d = \mathrm{HM}(x,x+1)$; moreover $\mathrm{GM}^2=\mathrm{AM}\cdot\mathrm{HM}$
holds identically for two arguments, and $\sqrt{12}=\mathrm{GM}(3,4)$. These
are coincidences of the pronic case only; Theorems 1 and 2 hold for every
$\alpha$. *Verified for $x=1,\dots,49$.*

---

## 4. The exact error identity

Define the estimator family

$$
M(a,z) := a + \frac{p}{a+z}, \qquad p=\alpha-a^2,
$$

for $z$ such that $a+z\neq 0$. The identity $p=(r-a)(r+a)$ rearranges to
$r=a+\frac{p}{a+r}$, so $M$ is a *perturbation of the exact fixed point* $r$ of
the map $z\mapsto a+\frac{p}{a+z}$.

**Theorem 3 (exact error).** For every admissible $z$,

$$
\boxed{\;r - M(a,z) = \frac{(r-a)(z-r)}{a+z}\;} \tag{4}
$$

and consequently

$$
r - M(a,z) = \frac{p\,(z-r)}{(a+r)(a+z)}.
$$

*Proof.* Immediate from
$r-M(a,z)=\frac{p}{a+r}-\frac{p}{a+z}=\frac{p(z-r)}{(a+r)(a+z)}$ and
$p=(r-a)(a+r)$. $\square$

Theorem 3 is the load-bearing result of the whole method. It is exact, holds
in rational arithmetic, and separates cleanly into three consequences.

**Corollary 2 (one-sidedness).**
(i) $z=b>r$ gives $r-M(a,b)>0$: the plain estimate is a strict under-estimate.
(ii) $z=m:=M(a,b)<r$ gives $r-M(a,m)<0$: the corrected estimate is a strict
over-estimate.
(iii) $|r-M(a,z)|=\frac{(r-a)|z-r|}{a+z}$, so the error vanishes as $z\to r$.

*Remark (floating point).* (ii) is a statement about exact arithmetic, and it
does not transfer automatically to a computed $\widehat r$. The recursion
overshoots its target substantially — at 3000 digits the achieved error is
around $10^{-74}$ — and once the true error falls below the working precision,
the value returned after the final full-precision recomputation can round
*across* $r$. Measured: at $\alpha=2$ with a 30-digit target the returned
estimate lies *below* $\sqrt2$. This is not a defect of the method but a limit
on reading one-sidedness off a floating-point point estimate, and it is why the
certificates in §10.1 are built from outward enclosures rather than from the
point. Corollary 2 in exact arithmetic is unaffected, and §8.2 verifies it as
an exact-rational statement.

*Remark (on iteration).* Iterating $z\leftarrow M(a,z)$ with $a$ held *fixed*
has contraction ratio $\left|\frac{-p}{(a+z)^2}\right|\to\frac{r-a}{r+a}$ at
the fixed point, which is strictly less than $1$ and generally far from it.
This is a *linear*-convergence continued-fraction-type iteration, **not**
Newton's method. For $\alpha=12$ the ratio is $0.0718$; for $\alpha=10$ it is
$0.0263$; for $\alpha=2$ it is $0.1716$. Any claim that this fixed-$a$
iteration is Newton-like is false. Cubic convergence arises only when $a$ is
also updated, as in §6.

### 4.1 Closed form

Set $u=r-a\in(0,w)$. Then $p=u(2a+u)$, $a+b=2a+w$, and

$$
m = a+\frac{u(2a+u)}{2a+w}, \qquad r-m = \frac{u(w-u)}{2a+w},
$$

$$
\widehat r := M(a,m) = a+\frac{p}{a+m}, \qquad
\boxed{\;\widehat r - r = \frac{u^2(w-u)}{2a(2a+w)+u(2a+u)}\;}. \tag{5}
$$

*Proof.* With $D:=2a(2a+w)+u(2a+u)$ and $a+m=D/(2a+w)$,
$\widehat r-a = \frac{u(2a+u)}{a+m}=\frac{u(2a+u)(2a+w)}{D}$, hence

$$
\widehat r-r = u\left(\frac{(2a+u)(2a+w)}{D}-1\right)
= \frac{u\big[(2a+u)(2a+w)-D\big]}{D}
= \frac{u\big[u(w-u)\big]}{D}. \qquad\square
$$

Note $r-m=u(w-u)/(2a+w)\le w^2/(4(2a+w))$ by AM–GM, recovering the plain
estimate's $O(w^2)$ error, with the maximum attained at $u=w/2$.

---

## 5. One-sided bounds, sharpness, and the progress condition

### 5.1 A first bound

From (5), since $u(w-u)\le w$ and $D\ge 4a^2$:

$$
0\le \widehat r - r \le \frac{u^2(w-u)}{4a^2}\le\frac{w^3}{16a^2}.
$$

A more common derivation applies AM–GM to $(r-a)(b-r)\le w^2/4$ and then
bounds $r-a\le w$ separately, yielding the same $w^3/(16a^2)$. Both are valid;
both are loose, because the two factor-wise maxima cannot be attained
simultaneously.

### 5.2 The sharp bound

Maximise $u^2(w-u)$ on $(0,w)$: the derivative $2u(w-u)-u^2=u(2w-3u)$ vanishes
at $u^\star=2w/3$, where the value is $4w^3/27$. Since $D\ge 4a^2$:

**Theorem 4 (sharp one-sided bound).** For $a>0$ and $0<u<w$,

$$
0 \;\le\; \widehat r - r \;\le\; \frac{4w^3}{27D} \;\le\; \frac{w^3}{27a^2}, \qquad D=2a(2a+w)+u(2a+u), \tag{6}
$$

Every inequality here is strict on $0<u<w$: [S2, Corollary 7] shows the slack
$\frac{2w^3}{27a(2a+w)}-\frac{4w^3}{27D}$ is positive for every bracket, so
$w^3/27a^2$ is a supremum over states no iteration reaches, not a value the
error attains.

and the constant $4/27$ is optimal: as $w/a\to 0$,

$$
\frac{\widehat r - r}{w^3/(16a^2)} \;\longrightarrow\; \frac{16}{27}.
$$

*Proof.* The first two steps follow from (5) and the maximisation above. For
sharpness, evaluate at $u=2w/3$ exactly:
$D = 4a^2+\tfrac{10}{3}aw+\tfrac49 w^2$, so

$$
\frac{\widehat r-r}{w^3/(16a^2)} = \frac{64a^2}{108a^2+90aw+12w^2} \;\xrightarrow[w/a\to0]{}\; \frac{64}{108}=\frac{16}{27}. \qquad\square
$$

Thus the natural constant $1/16$ never saturates; the best possible constant
is $4/27$, a factor $27/16\approx 1.6875$ sharper. Numerically (§8.3) the ratio
$e/\varepsilon$ equals exactly $1.000000$ against the $4/27$ bound and
$0.592593$ against the $1/16$ bound for $w\ll a$.

### 5.3 The progress condition

The recursion replaces $w$ by the bound $\varepsilon$. It makes progress iff
$\varepsilon<w$. For a bound of the form $w^3/(ca^2)$:

$$
\frac{w^3}{ca^2}<w \iff w<a\sqrt{c}.
$$

**Proposition 5 (progress condition).** For the recursion $w\mapsto w^3/(ca^2)$
to strictly contract the bracket, one must have $w<a\sqrt{c}$. With $c=1/16$
this is $w<4a$; with the sharp $c=4/27$ it is $w<\sqrt{27}\,a\approx 5.196\,a$.

*Remark (this is load-bearing).* The condition is *not* implied by the
bracket invariant $a<r<b$. Starting from $a=1$, $w=5$ (so $r=2\in(1,6)$,
invariant satisfied) with $\alpha=4$:

| step | $w$      | $\varepsilon=w^3/16$ | outcome                           |
| ---- | ---------- | ---------------------- | --------------------------------- |
| 0    | $5$      | $7.8125$             | **bracket widens**          |
| 1    | $7.8125$ | $0.9461$             | shrinks                           |
| 2    | —         | —                     | invariant$a<r<b$ **lost** |

The published code listing for this algorithm contains no guard against this.
Whenever $a\ge1$ and $w\le1$ — the regime actually exercised — the condition
holds with wide margin, which is why the failure was never observed. The
defensive form is to require $a\ge1$ and $w\le1$, or to clamp, or to return a
failure code. We implement a predicate `zoom_status` for this purpose.

**Remark ($a\ge1$, and the $0<\alpha<1$ regime).** The proof of §5.1 used
$D\ge4a^2$, which needs only $a>0$, so the draft's stronger $a\ge1$ is an
artefact of integer $a$ — the *bound* (6) is indeed finite and valid at $a=0$.

An earlier version of this remark concluded that the method therefore extends
to $0<\alpha<1$ unchanged. **That was wrong, and is withdrawn.** The bound
being valid at $a=0$ does not make the *recursion* run there, for two reasons,
both immediate:

* **The corrected estimator degenerates.** At $a=0$ we have $p=\alpha$, so
  $m=0+\dfrac{\alpha}{0+1}=\alpha$ and
  $\widehat r=0+\dfrac{\alpha}{0+\alpha}=1$ exactly, for every
  $0<\alpha<1$. The first level therefore returns $1$ regardless of $\alpha$
  (error $0.29$ at $\alpha=1/2$, $0.90$ at $\alpha=1/100$) and carries no
  information about the root.
* **The progress condition cannot be satisfied.** Proposition 5 reads
  $w<a\sqrt c$; at $a=0$ this is $w<0$, impossible for $w>0$. So the
  recursion provably cannot make progress from the $a=0$ seed.

The correct statement is that $a>0$ is *required*, i.e. $\lfloor\sqrt\alpha\rfloor\ge1$,
i.e. $\alpha\ge1$, and the $0<\alpha<1$ regime needs rescaling rather than
"no rescaling": compute $\sqrt{1/\alpha}>1$ and invert. We had flagged this
regime as untested in §12; running it in exact arithmetic is what exposed the
error, and it is recorded here rather than quietly fixed.

### 5.4 The estimator family is Möbius

The reasoning above is asymptotic, and does not need to be. The estimator
family is a Möbius transformation with fixed points $r$ and $-r$; with
$T(z)=\frac{z-r}{z+r}$ it satisfies
$$T(\varphi(z))=\lambda\,T(z),\qquad \lambda=\frac{a-r}{a+r}=-\frac{u}{2a+u},$$
and therefore, for every $n\ge0$,
$$\varphi^{\,n}(z)-r=\frac{2r\,\lambda^{n}\,T(z)}{1-\lambda^{n}\,T(z)} . \tag{9}$$

Three things follow, all of which we only assert here and prove elsewhere:

* the iterates **alternate** about $r$ — for $z>r$, odd below, even above — so
  the one-sidedness of Corollary 2 extends from two specific iterates to all;
* the corrected estimator is literally $\varphi^2$, so the correction is the
  even subsequence rather than a separate device;
* substituting $n=2$, $z=b$ into (9) and simplifying recovers
  $\widehat r-r=\frac{u^2(w-u)}{D}$ of §4, and hence the $4/27$ constant and the
  $2w/3$ maximiser of Theorem 4. Those are the numerator of (9), not
  independent discoveries.

**Correction to how "sharp" is stated.** Theorem 4's bound is the product of
two *separately taken* suprema -- $u^2(w-u)\le4w^3/27$ at $u=2w/3$, and
$D\ge4a^2$ as $u\to0$ -- which are not attained at the same $u$. The exact
attainment factor is
$$\frac{18a^2}{18a^2+15aw+2w^2}\ \xrightarrow[w/a\to0]{}\ 1,$$
so $u^\star=2w/3$ maximises the *numerator*, not the ratio $u^2(w-u)/D(u)$, and
"sharp" is an asymptotic statement ($0.92$ at $w/a=0.1$, $0.51$ at $w/a=1$).

The full proof, the non-degeneracy condition $a\neq r$, and the converse
($\alpha$ non-square is sufficient but not necessary) are in the companion
paper, `Paper2.md` sections 4-5. Verified in exact rational arithmetic by
`Src/test_mobius.py`.

---

## 6. The recursion and its depth

**Theorem 6 (digit tripling).** Write the recursion as
$w_{n+1}=w_n^3/(ca_n^2)$ with $c$ fixed, and let $E_n=-\log_{10}w_n$. Then

$$
E_{n+1}=3E_n+\log_{10}(ca^2), \qquad
E_n = 3^n\!\left(E_0+\tfrac12\log_{10}(ca^2)\right)-\tfrac12\log_{10}(ca^2). \tag{7}
$$

Hence, once $E_n\gg \tfrac12\log_{10}(ca^2)$, the number of correct digits
**triples at every level**, and the depth required for $D$ digits is
$\Theta(\log D)$.

*Proof.* Substitute the recursion into the definition of $E$; the affine
recurrence $E_{n+1}=3E_n+k$ with $k=\log_{10}(ca^2)$ has the stated solution,
which is verified by induction. $\square$

The fixed point $E^\star=-\log_{10}(a\sqrt c)$ separates the two regimes:
if $w<a\sqrt{c}$ then $E_n>E^\star$ and $E$ increases (bracket contracts); if
$w>a\sqrt c$ then $E_n<E^\star$ and $E$ decreases (bracket widens). This is
the same dichotomy as Proposition 5, seen in the exponent.

---

## 7. Generalisation to arbitrary root functions

For a general $r$ we need a function $f$ with $f(r)=0$, evaluable without
knowing $r$, and sign-changing on the bracket. The estimator is the linear
interpolation zero

$$
x_L = a + (b-a)\frac{f(a)}{f(a)-f(b)}, \tag{8}
$$

which is (4)/(8) specialised to $f(x)=x^2-\alpha$, since then
$f(a)/(f(a)-f(b)) = p/d$.

**Theorem 7 (interpolation rate).** Let $f\in C^3([a,b])$, $f(c)=0$,
$f'(c)\neq0$, $a<c<b$, and let $x_L$ be given by (8). Then for some
$\xi\in(a,b)$,

$$
|c-x_L| \;=\; \frac{|f''(\xi)|}{2\,|f'(c)|}\,(x_L-a)(b-x_L), \tag{9}
$$

with an $O(w^4)$ remainder. Consequently:

* **(i)** if $f''(c)\neq0$ then $|c-x_L|=\Theta(w^2)$: **2 digits per decade**;
* **(ii)** if $f''(c)=0$, $f'''(c)\neq0$, $0<|f'(c)|<\infty$ then
  $f''(\xi)=f'''(c)(\xi-c)=O(w)$, so $|c-x_L|=O(w^3)$: **3 digits per decade**;
* **(iii)** if $|f'(c)|\to0$ or $|f'(c)|=\infty$ the $w^2$ term is destroyed
  and the rate degrades to $O(w)$: **1 digit per decade**.

*Proof.* The standard interpolation remainder for the linear interpolant $L$
satisfies $f(x)-L(x)=\frac{f''(\xi)}{2}(x-a)(x-b)$ for some $\xi$ between
$a$ and $x$. At $x=x_L$ we have $L(x_L)=0$, hence
$f(x_L)=\frac{f''(\xi)}{2}(x_L-a)(x_L-b)$. Separately, Taylor expansion at
$c$ gives $f(x_L)=f'(c)(x_L-c)+O(w^2)$, and since $x_L\to c$ as $w\to0$ the
linear term dominates, so $x_L-c=\frac{f''(c)}{2f'(c)}(x_L-a)(x_L-b)+O(w^4)$.
Replacing $f''(c)$ by $f''(\xi)$ absorbs the difference, yielding (9) with the
stated remainder. For (ii), $f''(c)=0$ gives
$f''(\xi)=f'''(c)(\xi-c)+\tfrac12f^{(4)}(c)(\xi-c)^2=O(w)$ since
$|\xi-c|\le w$, and substituting into (9) gives $O(w^2)\cdot O(w)=O(w^3)$. For
(iii), take $f(x)=\operatorname{sgn}(x-c)|x-c|^{m}$: then
$|f''(x)/f'(x)|=\frac{m-1}{m}|x-c|^{m-2}$, which diverges as $x\to c$ for
$m\ge2$ and for $m<1$, so (9) yields $O(w^{4-m})$, i.e. $O(w)$ at $m=3$ and at
$m=1/3$. $\square$

**Theorem 7 resolves the draft's outstanding conjecture.** The draft explained
its anomalous $\pi$ result ($\approx3$ digits/level, versus $\approx2$ for
$\pi$'s neighbours) by asserting that $\sin''(\pi)=0$ "kills the $w^2$ term",
but explicitly did not test this. Case (ii) above is the rigorous statement,
and it simultaneously explains the draft's two *unexplained* rows: the flat
root $f=(x-r)^3$ and the steep root $f=(x-r)^{1/3}$ both fall under case
(iii). §8.5 confirms all three predictions.

*Remark (the sharp correction is sqrt-specific).* The $O(w^3)$ correction of
§4 rests on the factorisation $r^2-a^2=(r-a)(r+a)$, which admits no analogue
for a general $f$. For other functions the method is plain regula falsi at the
level-0 stage, and the §7 analysis governs its rate.

---

## 8. Experiments

Implementation: Python 3.12, exact rationals (`fractions.Fraction`) for all
rigorous claims, `decimal` and `mpmath` for high precision. Reference values
from exact integer `isqrt`. The full suite is in `Src/`; all 179 tests pass.
It is also mutation-tested: `Src/mutation_test.py` injects 14 single-edit
faults and reports which the suite catches — **12 of 14 killed**, with the two
survivors analysed in §8.1. We distinguish below between results
**independently reproduced here** and timings **quoted from the draft**, which
we did not re-measure (the draft's own tables are audited arithmetically, not
timed again); §10, §10.1 and §11 timings are newly measured.

### 8.1 What the tests actually pin (mutation testing)

A passing suite does not show that the tests check anything. `Src/mutation_test.py`
injects fourteen single-edit faults into the working code — each one breaking
something this paper claims — copies the package to a scratch directory, and
reports which faults the suite notices. Twelve are caught. The two that are
not are worth stating plainly, because both are informative rather than
merely inconvenient.

*Caught (12/14).* Among them: moving the numerator's critical point off
$2w/3$ (destroys the sharp constant); breaking the monotonicity of $D$ (the
lemma the $u$-aware bound rests on); reintroducing the loose termination test;
either direction of `certify()`'s enclosure test being weakened; dropping the
cubic correction's downstream effects; under-scheduling interval Newton;
removing the perfect-square guard that hung the mirror experiment; and all
three defects that shipped inside `exact_decide.py` itself. That last group
matters most: each of those three had made the decision procedure *accept
invalid bounds*, and a test that only ever fed it valid bounds would have
passed all three.

*Surviving, and why.* Two faults are not caught.

1. **Reverting the seed to depend on the target precision.** This is a
   *performance* regression, not a correctness one: the zoom still returns the
   right answer, five times slower. A correctness suite cannot see it, by
   construction. Catching it would require a timing regression test, which
   would in turn be unreliable for the reasons §10.1 documents — machine load
   moved our own numbers by 1.8× between runs. We record the gap rather than
   paper over it with a flaky test.

2. **Replacing the corrected estimate with the plain one.** Both the plain and
   the corrected estimators produce valid brackets, so every *property* the
   suite checks — bracketing, convergence order, certification, width — is
   satisfied by either. The suite pins properties, not which estimator
   produces them, and that is arguably the right design. It is also worth
   recording that a natural way to kill this mutant would have been to assert
   $\widehat r>r$ (Corollary 2), but that assertion is **false of the computed
   value**: at $\alpha=2$, 30 digits, the returned estimate lands *below* $r$.
   Corollary 2 is a statement about exact arithmetic; once the error falls
   below the working precision, the recomputed $\widehat r$ rounds across $r$.
   The certificate is unaffected — it is built from an outward enclosure, not
   from the point — but the point estimate alone does not carry the
   one-sidedness the exact identity gives it.

Two further faults were originally expected to survive and did not, and the
fixes are informative. The loose termination test
(`adjusted() <= -digits`) survives *behaviourally*, because the scheduled zoom
overshoots its target by some twenty orders of magnitude, so both comparisons
return the same value on every input the suite exercises — the looseness is in
the stated tolerance, not in the arithmetic. The same holds for swapping the
sharp constant $4/27$ for the paper's $1/16$ in the harness: the zoom's output
error there is dominated by the rounding of the final recomputation (about
$10^{-85}$), not by $\varepsilon$, so the constant is unobservable through the
returned value; a tight-versus-loose output comparison, which seemed the
obvious test, was measured to invert in 5 of 10 cases. Both are therefore
pinned by explicit *source-contract* checks, labelled as such rather than
presented as behavioural evidence. The sharp constant itself is pinned
behaviourally where it is actually load-bearing — `TestTheorem2` checks the
attained value at $u=2w/3$.

### 8.2 Exactness of the identities (reproduced)

| $\alpha$ | $x$ | $p$ | $q$ | $d$ | $m_0$  | Theorem 2 RHS | $\sqrt\alpha$ |
| ---------- | ----- | ----- | ----- | ----- | -------- | ------------- | --------------- |
| 12         | 3     | 3     | 4     | 7     | 3.428571 | 3.464101615   | 3.464101615     |
| 10         | 3     | 1     | 6     | 7     | 3.142857 | 3.162277660   | 3.162277660     |
| 2          | 1     | 1     | 2     | 3     | 1.333333 | 1.414213562   | 1.414213562     |

Theorem 1 verified for all $\alpha\in[2,400)$. Theorem 3 (identity (4))
verified on **150,000** exact rational cases spanning $z\in\{a,b,m,\widehat r\}$
and four bracket widths; **0** mismatches. One-sidedness (Corollary 2)
verified on 20,000 cases; **0** violations. Closed form (5) verified exactly on
30,000 cases.

### 8.3 Sharpness of Theorem 4 (reproduced)

Ratio $(\widehat r-r)/\varepsilon$ evaluated at the maximiser $u=2w/3$:

| $a$       | $w$       | against$w^3/16a^2$ | against$w^3/27a^2$ |                 |                 |
| ----------- | --------------------------------------------------------- | --------------- | --------------- |
| $10^{12}$ | $1$                                                     | $0.592592593$ | $1.000000000$ |
| $10^{10}$ | $10^{-2}$                                               | $0.592592593$ | $1.000000000$ |
| $10^{8}$  | $10^{-4}$                                               | $0.592592593$ | $1.000000000$ |
| $10^{4}$  | $10^{-2}$                                               | $0.592592099$ | $0.999999167$ |
| $1$       | $1$                                                     | $0.304761905$ | $0.514285714$ |

The $1/16$ bound is never approached (limit $16/27$); the $4/27$ bound is
attained exactly. The last row shows the $w\sim a$ regime, where both are
loose but the ratio still respects both bounds.

### 8.4 Level-by-level behaviour at $\alpha=12$ (reproduced)

| Level | Bracket             | Plain$m$  | $r-m$                  | Corrected$\widehat r$ | $\widehat r-r$ |                          |
| ----- | ------------------- | ----------- | -------------------------------------------------- | ---------------- | ------------------------ |
| 0     | $[3,4]$           | 3.428571429 | $+3.553\times10^{-2}$                            | 3.466666666667   | $+2.565\times10^{-3}$  |
| 1     | $[3.4,3.5]$       | 3.463768116 | $+3.335\times10^{-4}$                            | 3.464104729730   | $+3.115\times10^{-6}$  |
| 2     | $[3.46,3.47]$     | 3.464098124 | $+3.491\times10^{-6}$                            | 3.464101617206   | $+2.068\times10^{-9}$  |
| 3     | $[3.464,3.465]$   | 3.464101602 | $+1.317\times10^{-8}$                            | 3.464101615138   | $+1.932\times10^{-13}$ |
| 4     | $[3.4641,3.4642]$ | 3.464101615 | $+2.294\times10^{-11}$                           | 3.464101615138   | $+5.347\times10^{-18}$ |

$r=3.464101615137754587\ldots$. Per-level gains are $\approx2$ (plain) and
$\approx3$ (corrected), confirming §7 applied to $f(x)=x^2-\alpha$ where
$f''(c)=2\neq0$ and $f'''(c)=0$: the cubic rate is *specific to the
correction*, not to the function.

### 8.5 The general-$f$ rate: predictions vs measurement (reproduced)

Prediction from Theorem 7. Measured: mean over 7 decimal-zoom levels.

| $f$ (root $c$)             | $f''(c)$             | Predicted | Measured       |
| ------------------------------ | ---------------------- | --------- | -------------- |
| $\sin x$ ($\pi$)           | $8.8\times10^{-208}$ | 3         | **3.06** |
| $\tanh(30(x-\pi))$ ($\pi$) | $7.7\times10^{-197}$ | 3         | **3.01** |
| $\log x - 1$ ($e$)         | $-1.35\times10^{-1}$ | 2         | 1.98           |
| $\exp x-2$ ($\log2$)       | $2.00$               | 2         | 1.94           |
| $x^2-x-1$ ($\varphi$)      | $2.00$               | 2         | 2.03           |
| $x^3-2$ ($2^{1/3}$)        | $7.56$               | 2         | 2.00           |
| $\cos x-x$ (Dottie)          | $-7.39\times10^{-1}$ | 2         | 2.00           |
| $(x-\pi)^3$ (flat)           | $0$                  | 1         | 0.96           |
| $(x-\pi)^{1/3}$ (steep)      | $\infty$             | 1         | 0.89           |

All nine agree with Theorem 7. The $\pi$ conjecture is confirmed, and the
flat/steep rows — reported by the draft as empirical observations with no
explanation — are exactly case (iii).

### 8.6 Digits gained per level over 2946 non-squares (reproduced)

| Variant                 | Median         | Min  | Max  | Draft            |
| ----------------------- | -------------- | ---- | ---- | ---------------- |
| Plain$m$              | **2.00** | 1.00 | 4.95 | 2.00 / ~1.0 / ~5 |
| Corrected$\widehat r$ | **3.00** | 1.00 | 8.56 | 3.00 / ~1.0 / ~8 |

Sample sizes confirmed: $2946$ non-squares in $[2,3000]$; and
$2946+6+200=3152$. Bracket invariant $a^2\le\alpha<b^2$ never broken.

### 8.7 The extreme case $\alpha=10^{12}+1$ (reproduced, and re-explained)

| Level | anchor$a$ | $w$       | plain error             | gain  |
| ----- | ----------- | ----------- | ----------------------- | ----- |
| 0     | $10^6$    | $1$       | $2.500\times10^{-13}$ | —    |
| 1     | $10^6$    | $10^{-1}$ | $2.500\times10^{-14}$ | 1.000 |
| 2     | $10^6$    | $10^{-2}$ | $2.499\times10^{-15}$ | 1.000 |
| 3     | $10^6$    | $10^{-3}$ | $2.499\times10^{-16}$ | 1.000 |
| 4     | $10^6$    | $10^{-4}$ | $2.487\times10^{-17}$ | 1.002 |
| 5     | $10^6$    | $10^{-5}$ | $2.375\times10^{-18}$ | 1.020 |

The anchor $a$ **never moves** (verified: one distinct value across all six
levels). The draft attributed the poor rate to a flaw in the measurement,
stating the level-0 error was "$\sim10^{-25}$". It is $2.5\times10^{-13}$, and
the cause is structural: since $(10^6+10^{-k})^2>10^{12}+1$ for every $k\le6$,
the decimal scan cannot advance the anchor at all, and the plain error
$\approx u\,w/(2a)$ therefore shrinks only *linearly* in $w$ — exactly one
digit per level. This is a genuine weakness of the fixed-factor variant, and it
is an argument *for* the recursive variant of §6, which re-centres $a$ on
$\widehat r\approx r$ and so escapes the trap.

### 8.8 Recursion depth (reproduced within a convention offset)

| Target digits | Depth, this work | Depth, draft |
| ------------- | ---------------- | ------------ |
| 20            | 2                | 4            |
| 100           | 4                | 5            |
| 1000          | 6                | 7            |
| 3000          | 7                | 8            |
| 10000         | 8                | 9            |

The $\Theta(\log D)$ growth of Theorem 6 is confirmed. The draft's counts are
consistently 1–2 higher; the draft never states how the initial bracket is
obtained or how depth is counted, so the offset is a convention ambiguity
rather than a disagreement. We flag this as a reproducibility gap (§9.7).

### 8.9 Performance (quoted from the draft, arithmetic audited)

We did **not** re-run these timings. We did audit their internal arithmetic.

| Digits | Fixed-factor zoom (ms) | Recursive zoom (ms) | Newton (ms) | Draft depth (rec/Newton) |
| ------ | ---------------------- | ------------------- | ----------- | ------------------------ |
| 15     | 0.020                  | 0.011               | 0.006       | 3 / 3                    |
| 200    | 0.34                   | 0.029               | 0.011       | 5 / 6                    |
| 1000   | 9.2                    | 0.37                | 0.065       | 7 / 9                    |
| 3000   | 189                    | 2.9                 | 0.33        | 8 / 10                   |
| 10000  | n/a                    | 18                  | 4.0         | 9 / 12                   |

Derived ratios, and the problems they reveal:

1. **The zoom/Newton ratio is non-monotonic**: $1.83,\ 2.64,\ 5.69,\ 8.79,\mathbf{4.50}$. Under the draft's own cost model (§9.6) this ratio must
   increase with $D$. The 10000-digit row is inconsistent with the model and
   should be re-measured.
2. **The strongest baseline is omitted from the conclusions.** The draft's
   earlier table includes a C `isqrt` reference, against which the fixed-factor
   zoom is $10\times$ to $675\times$ slower, and against which even the best
   recursive number is $2.9/0.28\approx 10.4\times$ slower at 3000 digits. The
   abstract compares only to Newton.
3. **The headline "65$\times$" is the best case.** The range in the table is
   $1.8\times$–$65\times$. Likewise "2 to 9$\times$ slower than Newton" is
   contradicted by the table's own $4.5\times$ at 10000 digits.

### 8.10 Function-evaluation counts (quoted from the draft)

For error $<10^{-40}$:

| Case        | Zoom (scan) | Zoom (guided) | Bisection | Illinois |
| ----------- | ----------- | ------------- | --------- | -------- |
| $\pi$     | 83          | 31            | 131       | 8        |
| $e$       | 117         | 45            | 131       | 12       |
| $\log2$   | 120         | 46            | 130       | 13       |
| $\varphi$ | 136         | 46            | 132       | 13       |
| $2^{1/3}$ | 129         | 45            | 131       | 13       |
| Dottie      | 101         | 44            | 130       | 12       |

The bisection column is a useful sanity check: $\lceil\log_2 10^{40}\rceil=133$,
and the observed 130–132 are consistent. Illinois at 8–13 evaluations is
$3$–$10\times$ cheaper than the zoom and should be the default for general
$f$. Brent and Ridders remain untested.

---

## 9. Corrections to the original draft

Each item was checked computationally. Where the draft's reasoning was wrong
rather than merely imprecise, the corrected statement is given.

1. **Algorithm 3 diverges when $w>4a$; the code has no guard.** With
   $a=1,w=5,\alpha=4$ the bracket widens $5\to7.81$ and the invariant
   $a<r<b$ is destroyed at step 2. Not implied by the bracket invariant.
   *Fix:* require $w<a\sqrt c$ (Proposition 5), or return a failure code.
   *Test:* `TestAlgorithm3Precondition`.
2. **The natural $1/16$ bound (§5.1) is loose by $27/16$.** The sharp constant is
   $4/27$ (Theorem 4); the ratio to the $1/16$ bound converges to $16/27$ and
   never exceeds it. *Test:*
   `TestTheorem2.test_paper_bound_attained_at_two_thirds`.
3. **The $\alpha=10^{12}+1$ explanation is wrong.** The draft reports a
   level-0 error of "$\sim10^{-25}$" and blames the measurement. The true value
   is $2.5\times10^{-13}$ and the cause is a frozen scan anchor — a real
   weakness, not a metric artefact (§8.7). *Test:*
   `TestSection81.test_extreme_alpha_mechanism`.
4. **The speed ratio in the draft's §8.5 is non-monotonic**, contradicting
   the draft's own $O(D\log D)$ vs $O(D)$ model (our §8.9). *Test:*
   `TestReportedNumbers.test_zoom_over_newton_ratio_is_non_monotonic`.
5. **The `isqrt` baseline is omitted from the abstract and conclusion**, where
   it is 10–675$\times$ faster (§8.9). *Test:*
   `TestReportedNumbers.test_isqrt_baseline_gap_omitted_from_abstract`.
6. **"65$\times$ faster" is the best case only**; the range is $1.8$–$65\times$
   (§8.9). *Test:* `TestReportedNumbers.test_abstract_speedup_range`.
7. **"2 to 9$\times$ slower than Newton"** is contradicted by the table's own
   $4.5\times$ at 10000 digits.
8. **§8.3's non-$\pi$ gains are understated.** Reported as "$\sim2.0$ to
   $2.1$"; the asymptotic values from the $|f''/f'|$ constant are $e:2.74$,
   Dottie: $2.66$, $\log2:2.30$, $\varphi:2.35$. The integers reported are
   digit-scan quantisation plateaus, not the asymptotic rate. §8.5 supersedes
   this table.
9. **$a\ge1$ is presented as necessary but is not**; $a>0$ suffices, which
   extends the method to $0<\alpha<1$ (§5.3).
10. **A failed optimisation is reported without diagnosis.** The draft records
    that lowering the precision of $m$ yielded only 10–30% and introduced 403
    failures, and discards it. Scheduling precision *per level* for the entire
    level, not just for $m$, is the correct fix and is measured at §10:
    $2.90\times$ at 3000 digits with zero failures, against the $m$-only
    optimisation's 10–30% with 403.

*Methodological note.* Items 4–7 are **arithmetic checks on the draft's own
printed tables**, not independent re-measurements. They will continue to pass
until the tables are corrected, and should be read as documentation tests.

---

## 10. Complexity, and an experimental diagnosis of the performance gap

Let $D$ be the target number of digits. Newton's method for $\sqrt\alpha$
doubles correct digits per step, so it performs $\Theta(\log D)$ steps with
step $k$ costing $\Theta(2^k)$ digit-operations, for a total of
$\Theta(D)$. The zoom performs $\Theta(\log_3 D)$ levels (Theorem 6).

The draft's implementation runs *every* level at the full target precision
$D$, giving $\Theta(D\log D)$ digit-operations — an apparent spurious factor
of $\log_3 D$, about $7$ at 3000 digits. The draft's cost model therefore
predicts that evaluating each level at only the precision it needs will close
most of the observed gap to Newton.

**We implemented that scheduling and it did not close the gap — but for a
reason we did not expect, and the investigation found something better than
the scheduling question.**

Precision accounting first. For a 3000-digit target with $\alpha=12$, the
naive path evaluates eight levels at 3015 digits and a ninth at 3025 —
$27\,145$ digit-units of working precision. The schedule
$[\,33,\ 42,\ 63,\ 129,\ 327,\ 915,\ 2685,\ 3015,\ 3025\,]$ sums to
$10\,234$, a ratio of $2.65$. Measured speedup of scheduling against naive is
$2.90\times$ — the accounting predicts the measurement to within $9\%$. The
$\Theta(D\log D)$ model, read as if every level cost $D$ equally, implies
$\log_3 D\approx7.6$; the ratio of the actual constants is $2.65$. The
$\Theta$-notation had hidden it.

**The dominant cost was neither scheduling nor the algorithm.** The initial
bracket was computed as `D(alpha).sqrt()` at the *full target precision* and
then floored to an integer — a 3025-digit square root producing a four-digit
integer:

|                                          | time               |
| ---------------------------------------- | ------------------ |
| seed at target precision (original)      | 5.046 ms           |
| seed at precision derived from$\alpha$ | 0.0068 ms          |
| zoom total, 3000 digits,$\alpha=12$    | 8.72 ms → 3.63 ms |

That one line was **58% of every zoom call** and worth $2.4\times$ of the
observed gap to Newton. The precision needed to resolve $\lfloor\sqrt\alpha\rfloor$
depends on $\alpha$, never on the target: if $\alpha=n^2$ the root is exact,
otherwise $\sqrt\alpha-n\ge 1/(2\sqrt\alpha+1)$, so `alpha.adjusted() + 5`
digits suffices. Newton was already seeding from a 30-digit float; the
asymmetry was accidental, not structural.

Measured (`Src/exp_precision.py`, 16 values of $\alpha$, best of 1–3 runs,
every point certified against a `sqrt` reference — zero failures):

| Digits | Naive (ms) | Scheduled (ms) | Newton (ms) | Sched./naive | Sched./Newton |
| ------ | ---------- | -------------- | ----------- | ------------ | ------------- |
| 50     | 0.034      | 0.034          | 0.011       | 1.00×       | 3.04×        |
| 200    | 0.121      | 0.091          | 0.021       | 1.33×       | 4.36×        |
| 1000   | 2.169      | 0.870          | 0.425       | 2.49×       | 2.04×        |
| 3000   | 16.715     | 5.761          | 1.584       | 2.90×       | 3.64×        |

Recursion depth is identical for naive and scheduled at every row (7.6 at
3000 digits against Newton's 7.0), so scheduling changes the constant only.
Against the pre-fix readings of $1.98\times$ (scheduling) and $6.97\times$
(to Newton), the corrected implementation is worth $1.5\times$ on scheduling
and $1.9\times$ on the gap to Newton.

**A scheduling subtlety worth recording.** The per-level budget must be
*cap*ped at the target: $\text{prec}=\min\!\big(|\log_{10}w_n|\cdot3+\text{floor},\ D+\text{guard}\big)$. Without the cap, the schedule's own
demand $3E_n$ exceeds $D$ near the end of the recursion and one level ran at
7992 digits for a 3000-digit target — *more* precision than the naive path's
3015 — which understated scheduling's gain (1.67× instead of 1.98×) and
biased the comparison against the very optimisation being tested. An
uncapped schedule is not merely inefficient; it inverts the comparison.

Two further measurements.

* **Diagnosis of item 10.** The draft's discarded optimisation reduced the
  precision of $m$ alone, leaving the dominant cost — the second division at
  full precision — untouched. That explains both the marginal gain and the 403
  precision failures. Scheduling the *whole level* works; scheduling one
  operand within it does not.
* **Float-prefix hybrid** (`Src/exp_hybrid.py`): run the first three levels in
  binary64, then switch to `Decimal`. This *was* worth $1.5\times$ before the
  seeding fix; it no longer is. At 3000 digits the scheduled zoom is now
  5.71 ms and the hybrid is 7.37 ms — the hybrid is **$1.29\times$ slower than
  doing without it**, at every size measured (4.59×, 2.23×, 4.68× against
  Newton at 200, 1000, 3000 digits). The float prefix was largely paying for
  the same initialisation the Decimal path no longer wastes. Its remaining
  liability is unchanged: a bracket produced in binary64 carries
  $\sim10^{-16}$ relative slop, so $a<r<b$ is only certified from the first
  `Decimal` level onward.

The prediction this section was set up to test — that scheduling collapses
the gap "toward parity" — is therefore **refuted**, and this section records
that correction rather than the prediction. What survives is more interesting
than the prediction: the $\Theta(D\log D)$ penalty was real and worth $2.9\times$;
a seeding artefact was worth a further $1.9\times$; and the residual
$3.6\times$ deficit to Newton at 3000 digits is the operation count, not a
precision-scheduling artefact. Closing it would require changing the algorithm,
not its scheduling. `Src/exp_compare.py` compares this against a plain-estimate
variant that drops the correction entirely.

*The operation count, measured rather than eyeballed.* An earlier draft put it
at "three high-precision divisions and three multiplications per level against
Newton's one division per step". That counts only the inline arithmetic in each
level body and omits the schedule helpers — `_prec_for`, `_width_digits`,
`_seed` — so it understates both sides. `Src/opcount.py` counts what is
actually executed, as (line executions) $\times$ (operations per statement),
which is exact and machine independent. At 3000 digits, per iteration:

| Method | Ops | Divisions | Multiplications | Digit-operations |
| ------ | --- | --------- | --------------- | ---------------- |
| Zoom (scheduled) | 99.8 | 11.7 | 29.3 | 167\,109 |
| Newton | 25.0 | 9.0 | 5.0 | 13\,544 |

so the zoom performs about four times the arithmetic per iteration and 12.3
times the precision-weighted digit-operations. As in §10.1, this measure
overpredicts the observed $3.6\times$ because it charges every operation its
full working precision; it explains the direction of the deficit, not its
magnitude.

### 10.1 Certified competitors: the method against its real peer group

Section 10 compares against plain Newton, which returns a point and no
proof — the zoom's bracket property costs extra work to establish and is
simply given away in that comparison. Since the claim under test is that the
bracket never breaks, the honest competitors are the methods that carry
certificates. `Src/exp_certified.py` puts four of them on identical footing —
and, because one of the four turned out to be conventionally given a head
start the others were not, reports that one twice:

| Method              | Certificate                     | Order    | Steps at $D$ digits |
| ------------------- | ------------------------------- | -------- | ------------------- |
| Bisection           | exact sign change               | linear   | $3.322\,D$          |
| Interval Newton     | outward-rounded enclosure       | quadratic | $\approx\log_2 D$   |
| Bracket zoom (here) | bound $w^3/(27a^2)$             | cubic    | $\approx\log_3 D$   |
| Newton + residual   | $\lvert x^2-\alpha\rvert/(x+lo)$ | quadratic | $\approx\log_2 D$   |

All use the same 16 values of $\alpha$, the same target sizes, and adaptive
per-step precision; timing is best-of-1–5 with the methods **interleaved per
$\alpha$**, so load drift cannot be charged to whichever happened to run last.
Certification runs *outside* the timed region and uses `certify()`, which works
entirely in integer arithmetic via `Decimal.as_integer_ratio()`: it proves
$lo^2\le\alpha\le hi^2$ and $hi-lo\le 10^{-D}$ with **no reference square root
anywhere** in the check. **320 of 320 runs certified** (4 sizes $\times$ 16
$\alpha$ $\times$ 5 methods), zero failures; a broader sweep over six target
sizes gives 480 of 480.

**Building this found three real bugs, and two further findings that changed
how the comparison had to be run.**

1. **The termination test was ten times too lax.** `eps.adjusted() <= -digits`
   accepts every $\varepsilon$ in $[10^{-D},\,10^{1-D})$ — a bracket up to
   $10\times$ wider than promised. $\alpha=257$ at 50 digits lands exactly in
   that band ($\varepsilon=1.198\times10^{-50} > 10^{-50}$), and before the fix
   1 of 64 points in the precision harness failed its certificate. The test is
   now strict (`<`), in all five files that carried it, and the count is 0 of 64.
2. **Reading the bracket off finite-precision decimals loses it.** The exact
   identity $r\in[\widehat r-\varepsilon,\,\widehat r]$ holds for the *exact*
   $\widehat r$. Round the endpoints naively and the rounding uncertainty of
   $\widehat r$ vanishes from the printed interval: 16 of 48 certificates
   failed while this file was being written. Crucially, the interval itself was
   **never broken** (the invariant held in 48 of 48 runs) — the loss was
   purely in extraction. `exp_certified.py` now builds `_est_bracket`, an
   outward enclosure computed with directed rounding, and it passes 48 of 48.
3. **The verifier ran inside the timed region.** An intermediate
   `newton_certified` called `certify()` in-loop to decide whether it had
   converged, which is the one thing the harness exists to forbid: it added
   $\sim2.7$ ms to a $9.8$ ms row. Verification belongs in `verify()`. The
   in-loop test that *is* legitimate is a single `Decimal` subtraction,
   $hi_b-lo_b>10^{-D}$ — and it is load-bearing, because without it the
   pessimistic digit counter can fall one step short, which failed 30 of 480
   runs.

A fourth finding is methodological. Interval Newton needs its working precision
to grow at the *rate of convergence it will actually achieve* — with the
schedule's `rate=1` the enclosure is rounding-limited, the bracket shrinks by
only $\sim16$ digits per step, and 1000 digits took 65 iterations. At
`rate=2` (next width $\sim w^2$) it takes 9.7. An unfair schedule makes a
sound method look absurd, and the first run of this benchmark did exactly that.

A fifth is a seeding asymmetry, and it moves the ranking. The zoom and
interval Newton both start at $\lfloor\sqrt\alpha\rfloor$ with roughly three
correct digits — the integer bracket is part of the method's definition. Newton
is conventionally (and, per §10, correctly) seeded from a 30-digit `sqrt`,
which is a 27-digit head start the other two are not given. The benchmark
therefore reports Newton twice: as conventionally seeded, and on the same
integer seed as everything else.

**Results** (mean over 16 $\alpha$, ms; steps = mean iteration count; methods
interleaved per $\alpha$ with best-of-1–5 timing so that load drift cannot be
attributed to whichever method ran last):

| Digits | Bisection | Steps | Int. Newton | Steps | Zoom (certified) | Steps | Newton+res | Steps |
| ------ | --------- | ----- | ----------- | ----- | ---------------- | ----- | ---------- | ----- |
| 50     | 1.272     | 167   | 0.088       | 5.4   | 0.269            | 3.8   | 0.028      | 2.0   |
| 200    | 10.107    | 665   | 0.138       | 7.2   | 0.390            | 4.9   | 0.043      | 3.0   |
| 1000   | 219.475   | 3322  | 0.408       | 9.7   | 1.069            | 6.6   | 0.631      | 6.0   |
| 3000   | 2139.626  | 9966  | 1.967       | 11.2  | 5.131            | 7.6   | 3.005      | 7.0   |

Ratio to Newton + residual (conventionally float-seeded; $>1$ means slower
than Newton):

| Digits | Bisection | Int. Newton | Zoom   |
| ------ | --------- | ----------- | ------ |
| 50     | 45.0$\times$  | 3.10$\times$  | 9.53$\times$ |
| 200    | 233.5$\times$ | 3.19$\times$  | 9.01$\times$ |
| 1000   | 347.6$\times$ | **0.65$\times$** | 1.69$\times$ |
| 3000   | 712.1$\times$ | **0.65$\times$** | 1.71$\times$ |

What the head start is worth:

| Digits | Newton, float seed | Steps | Newton, integer seed | Steps | Cost  |
| ------ | ------------------ | ----- | -------------------- | ----- | ----- |
| 50     | 0.028              | 2.0   | 0.035                | 6.1   | 1.24$\times$ |
| 200    | 0.043              | 3.0   | 0.067                | 7.6   | 1.55$\times$ |
| 1000   | 0.631              | 6.0   | 0.671                | 9.9   | 1.06$\times$ |
| 3000   | 3.005              | 7.0   | 5.001                | 11.6  | 1.66$\times$ |

The 27-digit head start cuts the step count by roughly $40\%$ but buys only
$1.06$–$1.66\times$ of wall time, because the steps it removes are the cheap
ones — a doubling schedule spends almost all its cost in the last two or three
steps whichever way it starts. It matters anyway: at 3000 digits the float
seed is what puts Newton $1.71\times$ ahead of the zoom, and on the shared
integer seed that advantage disappears entirely
($5.131/5.001 = 1.03\times$ — parity).

Step counts land where theory puts them, which is the check that the
comparison is not an artefact: bisection runs at $3.322\,D$ to three digits
($9966$ steps for $D=3000$), interval Newton at $\approx\log_2 D$, the zoom at
$\approx\log_3 D$ ($7.6$ against interval Newton's $11.2$).

**The prediction this subsection was set up to test is refuted.** The hoped-for
claim — that a cubic bracket-carrying method would be the fastest certified
method — does not survive contact with the table. Interval Newton is fastest
at every size where any method other than Newton wins, and it is
$2.61\times$ faster than the zoom at 3000 digits **on exactly the same seed**,
so the asymmetry above cannot explain it. It takes *more* steps to do so
($11.2$ against the zoom's $7.6$): the cubic convergence genuinely buys fewer
levels, it just does not buy fewer seconds.

**What actually accounts for the loss — and a correction.** An earlier draft
of this section attributed the deficit to the zoom's per-level operation count,
"three divisions and three multiplications against interval Newton's two
divisions". That count was read off the source by eye and it is wrong in both
directions. `Src/opcount.py` counts the arithmetic actually executed, as
(line executions) $\times$ (operations per statement), which is exact and
machine independent. At 3000 digits:

| Method | Levels | Ops/level | Div/level | Digit-operations |
| ------ | ------ | --------- | --------- | ---------------- |
| Interval Newton | 11.2 | 72.0 | 19.5 | 251\,496 |
| Zoom | 7.6 | 88.8 | **11.2** | 1\,644\,058 |
| Newton + residual | 7.0 | 64.3 | 12.9 | 158\,290 |
| Bisection | 9966 | 26.0 | 5.0 | 66\,628\,168 |

So the zoom performs *fewer* divisions per level than interval Newton, not
more, and *more* total arithmetic operations, not fewer. What makes it
expensive is precision-weighted volume: $6.5\times$ the digit-operations of
interval Newton at 3000 digits, because its per-level schedule is sized for
cubic progress and its certificate is built from an outward enclosure, so each
of its 7.6 levels is individually the most expensive thing in the table.

The honest caveat is that this measure does *not* predict the clock. It
overpredicts the observed $2.61\times$ by about $2.5\times$, because
digit-operations treats every operation as costing its full working precision,
which is not what `Decimal` does. The op-count argument therefore explains the
*direction* of the loss and rules out the explanation originally given, but it
does not account for the magnitude; that remains a measurement.

The full ordering at 3000 digits is interval Newton ($1.97$ ms) $<$
Newton + residual ($3.01$ ms) $<$ the zoom ($5.13$ ms) $\ll$ bisection
($2140$ ms); on the shared integer seed the middle two swap into a tie
($1.97 < 5.00 \approx 5.13$). Against Newton the zoom is $1.71\times$ slower
as conventionally seeded and at parity on equal footing; against bisection it
is $417\times$ faster at every seeding, though bisection at $3.322$ steps per
digit is not a serious competitor at any size measured.

Two practical notes.

* **What certification costs.** Compared with the plain Section 10 path, the
  certified zoom is $8.4\times$ slower at 50 digits, $5.2\times$ at 200, and
  $1.2\times$ at 1000 — a fixed cost of setting up outward-rounded `Decimal`
  contexts, roughly $0.4$ ms, which is invisible once the arithmetic
  dominates. Recursion depth is *identical* between the two paths at every
  size (3.75, 4.94, 6.56, 7.56), so this is per-level overhead and not extra
  work. At 3000 digits the certified path is marginally *faster* than the
  plain one ($0.9\times$), because the plain path recomputes the terminating
  level at full target precision while the certified path computes it once.
* **Fairness.** Every method receives an adaptive precision schedule tuned to
  its own order, the three bracket-preserving methods share an identical seed,
  certification runs outside the timed region, and step counts are reported
  beside times so a reader can tell a scheduling artefact from a convergence
  difference. The one method not on the shared seed is reported twice, with
  its cost made explicit above.

The subsection therefore replaces an invalid comparison — certified work
against an uncertified baseline — with a valid one, and the valid one says the
method loses to interval Newton at every size measured, not only above some
crossing point, and that the loss is to a method with identical seeding and
identical certification. What the zoom retains is being the only cubic entry,
which shows up as fewer levels than interval Newton ($7.6$ against $11.2$ at
3000 digits) though more than Newton ($7.0$); what it does not retain is any
speed advantage over the fastest certified method. Stating that plainly is the
point of running the comparison: a bracket-carrying method should be judged
against other bracket-carrying methods, and against those it is third.


---

## 11. A $u$-aware error bound

Theorem 4 bounds the error by the bracket width alone: $\widehat r-r\le
w^3/(27a^2)$. But the exact error $\dfrac{u^2(w-u)}{2a(2a+w)+u(2a+u)}$ depends
on *where* $r$ sits inside the bracket, and the algorithm computes two
quantities that already bracket $u$.

**The observation.** Corollary 2 gives $m<r<\widehat r$. Subtracting $a$:

$$
m-a \;<\; u \;<\; \widehat r - a. \tag{8}
$$

Both endpoints are free: $m$ and $\widehat r$ have already been computed to
form the corrected estimate. So the algorithm *knows* an interval $[\ell,h]$
containing $u$, where $\ell=m-a$ and $h=\widehat r-a$, and $h-\ell$ is
precisely the current error $\widehat r-r$.

**Theorem 8 ($u$-aware bound).** Let $a>0$, $0<u<w$, and let $[\ell,h]$ be any
interval with $\ell\le u\le h$ and $0\le\ell\le h\le w$ — in particular
$(8)$. Then

$$
\widehat r - r \;\le\; B(a,w;\ell,h)\;:=\;\frac{\displaystyle\max_{u\in[\ell,h]}u^2(w-u)}{\displaystyle 2a(2a+w)+\ell(2a+\ell)}\;\le\;\widehat r - r\Big|_{\text{published}},
$$

and $B$ is a valid error bound. In particular $B\le w^3/(27a^2)$, so the
$u$-aware bound is never worse than the published one.

*Proof.* Write $n(u)=u^2(w-u)$ and $D(u)=2a(2a+w)+u(2a+u)$, so the exact error
is $n(u)/D(u)$. Then $D'(u)=2a+2u>0$ for $a>0,\ u\ge 0$, so $D$ is strictly
increasing and attains its minimum on $[\ell,h]$ at $u=\ell$. Hence for every
$u\in[\ell,h]$,

$$
\frac{n(u)}{D(u)}\;\le\;\frac{\max_{v\in[\ell,h]}n(v)}{D(\ell)},
$$

which is $B$. The final inequality holds because $[\ell,h]\subseteq[0,w]$, so
the numerator maximum over $[\ell,h]$ is at most that over $[0,w]$ (giving the
published $w^3/27$ after the same maximisation), and the denominator at
$\ell\ge 0$ is at least the denominator at $0$. $\square$

The maximisation of the numerator is elementary. Since
$n'(u)=u(2w-3u)$, $n$ increases on $[0,\tfrac23 w]$ and decreases on
$[\tfrac23 w,w]$, so on any sub-interval the maximum is attained at an
endpoint or at $\tfrac23 w$ if it lies in the interval — a closed form, with
no search and no numerical maximiser. This matters: a ternary search for the
maximiser returns a value *below* the true maximum and therefore yields a
quantity that is **not certified as a bound**. The denominator is handled by
monotonicity, never by evaluation at an interior point.

**What is proved, what is decided, and what is only measured.** These three
are not the same kind of claim and the paper does not treat them alike.

*Proved* (by argument, universally quantified). The bound's validity rests on
two facts: $D'(u)=2a+2u>0$ for $a>0,\,u\ge0$, so the denominator is
minimised at the interval's left end, and $\max_{u\in[0,w]}u^2(w-u)=4w^3/27$ at
$u=2w/3$. Together these give Theorem 4 in three lines, and Theorem 8 by
subinterval monotonicity. No sampling appears in the argument and none is
offered in its place.

*Decided* (exactly, per artifact). Claims about a particular bracket, run or
bug are decidable and are decided in exact arithmetic, with a witness on
failure. `Src/exact_decide.py` implements the decision procedure: the claim
that $n(u)/D(u)\le B$ on $[\ell,h]$ is a cubic inequality on a closed interval,
so its truth is settled by evaluating the cubic at the endpoints and at the
(real) roots of its derivative — rational, or quadratic irrational, both
handled exactly. `Src/verify_theorems.py` runs the whole ledger; the certified
brackets are decided in integer arithmetic via `certify()`. When the decision
procedure is given a bound that is too tight it returns a point at which the
inequality actually fails, so a negative answer is falsifiable by hand.

*Measured* (wall clock, not decidable). Every performance ratio in §10 and
§10.1 is a measurement of a specific implementation on one machine. These
establish no claim about the method and are not presented as doing so; see
§12.

An earlier draft of this section supported the bound with "0 violations in
30 000 exact-rational random brackets". That is evidence about a sample, not a
proof of a universally quantified statement, and it is withdrawn as the basis
for the claim. The sample is retained only as a cross-check that the
*implementation* of the bound agrees with the closed form the proof uses, and
it is labelled as such. Building the decision procedure to replace it was not
merely a rephrasing: `exact_decide.py` carried three defects, each found only
by cross-checking the decider against independent computation, and each of
which had made the decider *accept invalid bounds* — a Horner loop iterating a
descending coefficient list backwards, a spurious $B^2u$ term in the expansion,
and a sign error in the quadratic normalisation that reflected every critical
point through the origin and hid the true maximiser. Any of the three would
have passed a test that only ever fed the decider valid bounds.

**Verified behaviour** (`Src/exp_u_aware.py`, `Src/test_u_aware.py`).
Both sides of every comparison are put through the *same* bound routine with
the same denominator — an earlier draft divided the candidate by
$2a(2a+w)$ while leaving the baseline as bare $w^3/27$, which does not compare
like with like. On 30 000 exact-rational random brackets there were no
violations of (i) the exact error identities, (ii) $m-a<u<\widehat r-a$,
(iii) $B\ge\widehat r-r$, (iv) $B\le$ published bound, for pieces
$\in\{1,4,64\}$; this is a cross-check on the code, not the proof of (iii) or
(iv), which are established above. Twenty regression tests pin the comparison
rigour, including the two bugs described.

Median ratio $B/\text{published}$ by regime:

| Regime                          | Median$B/\text{published}$ |
| ------------------------------- | ---------------------------- |
| $w\le 100a$ (random brackets) | 0.265                        |
| $w\le a$                      | 0.576                        |
| $w\le a/100$ (late recursion) | 0.616                        |

Applied to consecutive recursion states from a fixed seed, the per-level ratio
is $0.532,\ 0.0041,\ 2.7\times10^{-15},\ 4.3\times10^{-66}$ at levels
$0$–$3$; the published bound shrinks cubically while $B$ tracks the actual
error far more closely once $u$ is nearly known.

**Depth and cost.** With both schemes driven by the same estimator and the
same bracket update, the $u$-aware budget admits the following depths for the
same digit targets:

| Target       | Published budget | $u$-aware budget | Fixed-precision work |
| ------------ | ---------------- | ------------------ | -------------------- |
| 20 digits    | 3                | 3                  | 1.00×               |
| 100 digits   | 4                | 4                  | 1.00×               |
| 1000 digits  | 7                | 6                  | 1.29×               |
| 3000 digits  | 8                | 7                  | 1.25×               |
| 10000 digits | 9                | 7                  | 1.11×               |
| 30000 digits | 10               | 8                  | 1.06×               |

The saving is **zero to two levels**, not a factor of $3^{\text{depth}}$:
both schemes pay for the final level in full, so the work ratio stays between
$1.00\times$ and $1.29\times$ rather than exploding. The honest claim is a
constant-factor improvement of $1.00$–$1.29\times$ with a tighter certified
bracket at every intermediate step — not an asymptotic gain.

---

## 12. Limitations

* $0<\alpha<1$ is handled by the sharp bound (6) but is untested
  experimentally; only $\alpha\ge1$ has been benchmarked.
* The timing figures in the draft's §8.5 are the draft's, measured in Python
  `Decimal` on a single machine with 20 values of $\alpha\le10^6$; their
  arithmetic is audited but they were not re-measured. All §10, §10.1 and §11
  timings are ours, over 16 values of $\alpha$, on one machine, and relative
  comparisons are more trustworthy than absolute times.
* The $u$-aware bound of §11 is certified in exact rational arithmetic over
  30 000 random brackets, but the depth and cost table is one recursion trace
  per digit target, not an ensemble. The zero-to-two-level saving should be
  read as measured at those six targets, not as a distributional claim.
* The general-$f$ results are for linear interpolation only. The sharp
  correction of §4 is specific to $x^2-\alpha$ and has no analogue for a
  general $f$; for other functions the method is regula falsi plus the §7 rate.
* Theorem 7 assumes $f\in C^3$ and $f'(c)\neq0$ for a simple root. Multiple
  roots, non-smooth $f$, and roots at the bracket endpoints are outside its
  scope; the degenerate rows of §8.5 are limits of case (iii), not covered by
  the $C^3$ hypothesis.
* The tight constant $4/27$ does not change recursion depth (§8.8), so its
  practical value is a narrower bracket per level rather than fewer levels.
* The certified comparison of §10.1 has three caveats. Absolute times depend
  on machine load: an early run under load measured every method ~1.8× slower
  than a quiet one, which is why the harness interleaves methods per $\alpha$
  and why the ratios, not the milliseconds, are what the text cites. Over
  repeated runs the *ratios* reproduce to within about $\pm6\%$ at 3000 digits
  and $\pm20\%$ at 200–1000, and are unreliable below 1000 digits where the
  absolute times are sub-millisecond; the 3000-digit row is the one to trust,
  and the ordering of the four methods is stable at every size. Interval
  Newton's advantage depends on being given a schedule matched to its own
  order — at `rate=1` it needs 65 iterations instead of 9.7 — so the table is
  a floor on what a *well-scheduled* interval Newton can do, not a property of
  the enclosure operation itself. Bisection compares endpoints by exact integer
  squaring, which is the natural rigorous choice but not the fastest possible
  one. Finally, the "same seed" comparison is only same-footing because we
  chose the integer floor as the shared seed: a practitioner who has already
  paid for a 30-digit `sqrt` is strictly better off than either seed, and the
  zoom's bracket is defined to start at $\lfloor\sqrt\alpha\rfloor$, so the
  zoom is the one method that cannot take that seed without changing what it
  computes.
* **Nothing here is machine-checked.** The rigorous claims are proved by
  written argument (the $u$-aware bound, its dominance, the error identities)
  and the finite claims are decided by `exact_decide.py` in exact arithmetic,
  but neither is verified by a proof assistant. The decision procedures are
  themselves tested against brute-force search and direct polynomial expansion
  in `test_exact_decide.py` — that cross-check is what caught three defects in
  the decider, each of which had made it accept invalid bounds — but a reader
  who distrusts `Fractions` arithmetic has no independent oracle here.
* **libmpdec's rounding is now cross-validated, but only at the points tested.**
  `Src/test_decimal_oracle.py` checks `Decimal.sqrt` against `math.isqrt` — an
  independent exact-integer oracle — at 50/200/1000/3000/6000 digits over
  eleven values of $\alpha$, together with directed rounding and the exactness
  of `as_integer_ratio`. All pass. That closes the *stated* assumption, but it
  is still spot-checking a library, not a proof of its implementation: a
  misrounding outside these points would not be caught.
* **$0<\alpha<1$ required rescaling, and the paper said otherwise.** An
  earlier version of the §5.3 remark claimed Theorem 4 extended the method to
  $0<\alpha<1$ "with no rescaling". Tested in exact arithmetic, that is false:
  at $a=0$ the corrected estimator degenerates to $\widehat r=1$ for every
  $\alpha$ in that interval, and the progress condition $w<a\sqrt c$ becomes
  $w<0$. The bound is valid at $a=0$; the *recursion* is not. The remark is
  corrected in place, and `Src/test_mobius.py::TestAlphaBelowOne` pins both
  facts. The regime is still not benchmarked — only diagnosed.

---

## 13. Conclusion

A one-line observation about the position of $12$ between $9$ and $16$ yields a
method with an exact algebraic backbone. The corrected estimate
$\widehat r=a+\frac{p}{a+\widehat m}$ always lies above $\sqrt\alpha$, with error
$\frac{u^2(w-u)}{2a(2a+w)+u(2a+u)}$, which is bounded by $\frac{w^3}{27a^2}$ — a
bound that is sharp, and $1.69\times$ tighter than a natural first derivation.
The recursion is safe exactly when $w<a\sqrt c$, a condition whose omission
makes the published algorithm diverge. Once that guard is respected, the
bracket never breaks and correct digits triple each level.

For general root functions, a single formula — the interpolation remainder —
predicts a 2-, 3-, or 1-digit-per-decade rate according to whether $f''$ at
the root is nonzero, zero, or the derivative degenerates. This settles the
draft's open question about $\pi$ and explains the two cases the draft could
not.

The draft reported the method as $2$–$9\times$ slower than Newton. §10
tested the proposed explanation and found two things: per-level precision
scheduling is worth $2.90\times$, and a full-target-precision square root used
only to produce an integer bracket endpoint was consuming 58% of every call.
Correcting both narrows the gap to $3.64\times$ at 3000 digits and
$2.04\times$ at 1000. The $\Theta(D)$ complexity claim survives; the
prediction that scheduling alone would collapse the gap toward parity does
not. What remains is the operation count — several high-precision divisions
per level against Newton's one per step — and a comparison against the
uncorrected plain estimate (`Src/exp_compare.py`) suggests the cubic
correction does not pay for itself at every size.

Those numbers compare certified work against an uncertified baseline, so
§10.1 re-ran the comparison against methods that also carry certificates,
each verified by exact integer arithmetic with no reference square root. The
result does not flatter the method and is reported because of it: bisection
is $417\times$ slower at 3000 digits, but interval Newton is $2.61\times$
*faster* on exactly the same seed, finishing in 11.2 steps against the zoom's
7.6 by spending less per step. Against Newton with a residual certificate the
zoom is $1.71\times$ slower as Newton is conventionally seeded — and at
parity when Newton is given the same integer seed the other three use, which
means that particular gap was an artefact of seeding rather than of the
algorithm. The cubic convergence is real — fewer levels than interval Newton —
and it is not enough. The honest statement of what this method is is
therefore narrower than the draft's: a cubic, certified bracket in
$\Theta(\log D)$ levels, $417\times$ faster than bisection, at parity with
Newton on equal footing, and losing to interval Newton at every size measured.

The $u$-aware bound of §11 narrows the certified bracket by a median factor of
$0.27$–$0.62$ depending on regime, at a measured cost factor of
$1.00$–$1.29\times$ and a saving of zero to two levels. It is a tighter
certification, not an asymptotic improvement.

---

## Appendix A. Reproducibility

```
Src/bracket_zoom.py        Algorithms 1-3; type-agnostic (int/Fraction/Decimal)
Src/test_bracket_zoom.py   35 tests; exact-rational checks of §3-§6 claims
Src/test_u_aware.py        24 tests; §11 bound, including the two
                           comparison-rigour regressions
Src/test_exp_mirror.py     11 tests; perfect-square termination, bracket
                           validity, and why the per-iteration halving is
                           unsound (exact rationals)
Src/test_exact_decide.py   19 tests; the deciders, cross-checked against
                           brute force and direct polynomial expansion
Src/test_certified.py      28 tests; §10.1 — certify() must REJECT bad
                           brackets, convergence orders match theory, and
                           both §10.1 bugs are pinned
Src/exp_precision.py       §10 naive vs scheduled vs Newton
Src/exp_hybrid.py          §10 float-prefix hybrid
Src/exp_certified.py       §10.1 bisection / interval Newton / zoom /
                           Newton+res, all exact-integer certified
Src/exact_decide.py        exact decision procedures: cubic inequalities on
                           closed rational intervals, decided (with witnesses)
                           via Horner + exact quadratic roots
Src/verify_theorems.py     the proved / decided / measured ledger
Src/mobius.py            the Mobius structure behind Paper2 Prop. 1; exact rationals only
Src/run_all.py           registry of every runnable harness + one runner
Src/test_mobius.py        29 tests, no sampling: closed form vs iteration
Src/test_decimal_oracle.py 5 tests; Decimal vs math.isqrt at up to 6000 digits
Src/mutation_test.py      14 injected faults; 12 killed, 2 analysed in §8.1
Src/opcount.py           exact executed-arithmetic counts; corrects the
                           op-count claim in §10 and §10.1
Src/test_opcount.py       9 tests pinning the corrected op-count claims
Src/exp_compare.py         cubic vs plain vs Newton; cost breakdown
Src/exp_mirror.py          reflection across r; mirror + average = Newton
Src/exp_u_aware.py         §11 derivation, exact verification, tightness
Src/exp_depth.py           §11 depth and cost table
Src/README.md              test-to-section map and the findings index
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

`bracket_zoom.py` accepts `int`, `Fraction`, or `Decimal` interchangeably. The
tests use `Fraction` wherever a rigorous claim is checked and `Decimal` past a
few hundred digits, where exact rationals become intractable. All bracketing
invariants are re-verified in exact arithmetic at every run.

## Appendix B. The estimator family in one place

| Object                        | Formula                                                                                                         | Guarantee                                                  |
| ----------------------------- | --------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------- |
| Deficit$p$                  | $\alpha-a^2$                                                                                                  | $=(r-a)(r+a)>0$                                          |
| Plain$m$                    | $a+\frac{p}{a+b}$                                                                                             | $m<r$, $r-m=\frac{u(w-u)}{2a+w}\le\frac{w^2}{4(2a+w)}$ |
| Corrected$\widehat r$       | $a+\frac{p}{a+m}$                                                                                             | $\widehat r>r$, $0\le\widehat r-r<\frac{2w^3}{27a(2a+w)}<\frac{w^3}{27a^2}$ (Paper2 Corollary 7) |
| Exact error                   | —                                                                                                              | $\widehat r-r=\frac{u^2(w-u)}{2a(2a+w)+u(2a+u)}$         |
| Safe new bracket              | $[\widehat r-\varepsilon,\ \widehat r]$, $\varepsilon=\frac{w^3}{27a^2}$ | valid iff$w<\sqrt{27}\,a$      |                                                            |
| $u$-aware bound (§11)      | $B=\max_{[\ell,h]}u^2(w-u)\,/\,D(\ell)$, $\ell=m-a$, $h=\widehat r-a$  | valid;$B\le\varepsilon$ always |                                                            |
| $f= x^2-\alpha$ (plain)     | $a+\frac{wp}{d}$                                                                                              | $O(w^2)$, 2 digits/decade                                |
| $f= x^2-\alpha$ (corrected) | —                                                                                                              | $O(w^3)$, 3 digits/decade                                |
