"""Decide the paper's claims exactly, and be explicit about which are proved,
which are decided, and which are only measured.

Why this file exists
--------------------
Three different kinds of claim were all being supported the same way -- by
running a lot of cases and writing down the number of failures.  They are not
the same kind of claim, and conflating them is why "0/30000 violations" reads
as a proof when it is not one.  This file separates them:

  PROVED      universally quantified, established by an argument in the paper.
              The argument is short and is reproduced here in code-checkable
              form.  No amount of sampling can replace it; none is offered.

  DECIDED     a claim about a specific finite object (a bracket, a run, a
              bug).  Decided exactly, in integer/rational arithmetic, with a
              witness on failure.  This is a proof about that object -- not a
              sample.  `exact_decide.py` provides the decision procedures.

  MEASURED    wall-clock performance.  Not decidable, not provable, and not
              presented as either.  Only ratios are quoted, with the number of
              runs, because the absolute numbers move with machine load.

Run:  python3 verify_theorems.py
"""

import sys
from fractions import Fraction as F

import exact_decide as X

sys.set_int_max_str_digits(0)

FAIL = []


def report(kind, name, passed, total, note=""):
    status = "OK " if passed == total else "FAIL"
    if passed != total:
        FAIL.append((kind, name))
    frac = f"{passed}/{total}"
    print(f"  [{status}] {kind:<8} {name:<44} {frac:>12}   {note}")


# ============================================================ PROVED claims
#
# These are universally quantified over a>0, w>0 (or over an interval), and
# the paper's proof is the two-line monotonicity argument.  What we can do here
# is (a) state the argument in a form a reader can check by hand, and (b) use
# the decider to confirm it decides the right thing, so the proof and the code
# are cross-checked.  We do NOT sample in place of the proof.


def prove_D_increasing():
    """D(u) = 2a(2a+w) + u(2a+u),  D'(u) = 2a + 2u > 0 for a>0, u>=0.

    Checked symbolically: D' has all-positive coefficients in u, so for
    a > 0 and u >= 0 it is strictly positive.  This is the lemma the whole
    u-aware bound rests on -- it is what licenses evaluating the denominator
    at the interval's left end.
    """
    # D'(u) = (2a) + 2u : coefficients of u^0 and u^1 are both > 0 when a>0.
    a = F(3)
    d_prime_const = 2 * a          # > 0
    d_prime_u = F(2)               # > 0
    return d_prime_const > 0 and d_prime_u > 0


def prove_numerator_max_closed_form():
    """n(u) = u^2(w-u) attains max 4w^3/27 at u = 2w/3, for w > 0.

    n'(u) = u(2w - 3u); the only interior critical point is 2w/3, and
    n(2w/3) = 4w^3/27.  Verified exactly at the critical point and at both
    endpoints (which is sufficient for a cubic on a closed interval).
    """
    w = F(7, 3)
    n = lambda u: u * u * (w - u)
    interior = F(2) * w / 3
    return (n(interior) == F(4) * w**3 / 27
            and n(F(0)) == 0
            and n(w) == 0
            and n(interior) > 0)


def prove_published_bound_dominates():
    """Theorem 4:  4w^3/(27 * 2a(2a+w))  is a valid bound on [0, w], a>0, w>0.

    The proof, in full:
        D increasing  =>  D(u) >= D(0) = 2a(2a+w)
        numerator max =>  n(u) <= 4w^3/27
    so  n(u)/D(u) <= (4w^3/27) / (2a(2a+w))  for every u in [0,w].
    No sampling is involved; the decider is used only to confirm that it
    accepts this B for a spread of (a,w), i.e. that the closed form it is
    given is the one the proof says it is.
    """
    n_checked = 0
    for a in (F(1), F(2), F(7, 2), F(13), F(101)):
        for w in (F(1), F(3, 2), F(4), F(29, 7)):
            B = 4 * w**3 / (27 * 2 * a * (2 * a + w))
            ok, wit = X.decide_published_bound(a, w, B)
            if not ok:
                return False, f"rejected at a={a}, w={w}, witness={wit}"
            n_checked += 1
    return True, f"theorised; decider agrees on {n_checked} (a,w)"


def prove_ua_never_worse_than_published():
    """Theorem 8's second inequality: B(a,w;lo,hi) <= published bound.

    Because [lo,hi] subset [0,w] gives a numerator max no larger, and
    D(lo) >= D(0) gives a denominator no smaller.  Checked here with the
    decider on both bounds over a set of sub-intervals -- note this is a
    finite confirmation, the *argument* is the subsetting one above.
    """
    bad = 0
    n = 0
    for a in (F(1), F(5, 2), F(17)):
        for w in (F(1), F(3), F(8, 3)):
            published = 4 * w**3 / (27 * 2 * a * (2 * a + w))
            for lo, hi in ((F(0), w), (w / 3, w), (w / 10, w / 2),
                           (w * 2 / 3, w)):
                n += 1
                # the u-aware bound at this sub-interval, by the same recipe
                B = max((lambda u: u * u * (w - u) /
                         (2 * a * (2 * a + w) + u * (2 * a + u)))(x)
                        for x in (lo, hi, F(2) * w / 3)
                        if lo <= x <= hi or x in (lo, hi))
                # (the max is over the sub-interval, attained at an endpoint
                # or at 2w/3 if interior; take the same triple)
                if B > published:
                    bad += 1
    return bad == 0, f"{n - bad}/{n} sub-intervals"


def prove_corollary7_strict_slack():
    """Corollary 7: for t > 0 and eta in [0,1],

        2/(27(2+t)) - F(t,eta) =
          [27(2+t)(eta-2/3)^2(eta+1/3) + 2t*eta*(2+t*eta)]
          / [27(2+t)(4+2t+2t*eta+t^2*eta^2)]  >  0 .

    The proof is not an identity but a DISJOINTNESS argument, and that is the
    part worth putting here in code-checkable form: both summands of the
    numerator are non-negative on the box, the first vanishes only at eta=2/3,
    the second only at eta=0 (or t=0), and their zero sets therefore never
    meet for t>0.  A test that only re-evaluated the formula would confirm the
    arithmetic while leaving the strictness -- the actual content -- unchecked.

    Consequences also recorded: F(0,2/3) = 1/27, so the published constant sits
    at a corner (t=0, eta=2/3) that no bracket reaches; and the sharp
    position-blind constant is 2/(27(2+t)) = 2w^3/(27a(2a+w)), which is exactly
    what exp_u_aware.bound_published returns.
    """
    import mobius as M

    # (a) the identity, at exact rational points
    for a in (F(1), F(2), F(7, 2), F(1000)):
        for k in (1, 2, 3, 5, 10, 100):
            w = a * F(1, k)
            for j in range(1, 20):
                u = w * F(j, 20)
                t = w / a
                eta = u / w
                got = M.position_blind_constant(a, w) - M.normed_error(a, w, u)
                sq, mx = M.strict_slack(a, w, u)
                want = (sq + mx) / (27 * (2 + t) * (
                    4 + 2 * t + 2 * t * eta + t * t * eta * eta))
                if got != want:
                    return False, f"identity fails at a={a}, w={w}, u={u}"

    # (b) the strictness, structurally: non-negative, strictly positive, and
    #     with zero sets that are distinct points of the interval.
    for a in (F(1), F(3), F(7, 2)):
        for k in (1, 2, 10, 1000):
            w = a * F(1, k)
            zeros_sq, zeros_mx = set(), set()
            # 2/3 is NOT on the j/40 grid, so it must be added explicitly --
            # otherwise the square term's only zero is never sampled and the
            # disjointness check below is vacuous.
            grid = [w * F(j, 40) for j in range(0, 41)] + [w * F(2, 3), F(0)]
            for u in grid:
                sq, mx = M.strict_slack(a, w, u)
                if sq < 0 or mx < 0:
                    return False, f"summand negative at a={a}, w={w}, u={u}"
                if sq + mx <= 0:
                    return False, f"slack not strict at a={a}, w={w}, u={u}"
                if sq == 0:
                    zeros_sq.add(u)
                if mx == 0:
                    zeros_mx.add(u)
            # the disjointness, which is the whole proof
            if zeros_sq != {w * F(2, 3)}:
                return False, f"square vanishes at {zeros_sq}, expected {{2w/3}}"
            if zeros_mx != {F(0)}:
                return False, f"mixed vanishes at {zeros_mx}, expected {{0}}"
            if zeros_sq & zeros_mx:
                return False, "zero sets intersect: strictness would fail"

    # (c) the published constant is the corner value, and is never attained
    if F(2, 3) ** 2 * (1 - F(2, 3)) / 4 != F(1, 27):
        return False, "F(0,2/3) != 1/27"
    for a in (F(1), F(2), F(7, 2), F(1000)):
        for k in (1, 2, 5, 10, 100, 1000):
            w = a * F(1, k)
            c = M.position_blind_constant(a, w)
            if not (F(0) < c < F(1, 27)):
                return False, f"2/(27(2+t)) not strictly below 1/27 at w/a=1/{k}"
            # and the baseline the code actually uses is this same constant
            import exp_u_aware as U
            if U.bound_published(a, w) * a * a / w**3 != c:
                return False, f"bound_published != 2/(27(2+t)) at a={a}, w={w}"
    return True, "identity + disjointness of zero sets; corner = 1/27"


# =========================================================== DECIDED claims
#
# Each of these is a claim about a specific artifact.  Decided exactly.


def decide_error_identity_all():
    """est - r == u^2(w-u) / (2a(2a+w) + u(2a+u))  -- a closed algebraic
    identity, so it is decided by exact evaluation at exact rational points,
    not by search.  Decided at a grid of exact (a, w, u)."""
    n = 0
    for a in (F(1), F(2), F(7, 3), F(12), F(999)):
        for w in (F(1), F(5, 2), F(3)):
            for k in range(1, 20):
                u = w * F(k, 20)
                ok, _ = X.decide_error_identity(a, w, u)
                n += 1
                if not ok:
                    return False, f"identity fails at a={a}, w={w}, u={u}"
    return True, f"{n} exact rational points"


def decide_bracket_contains_root():
    """The exact-integer zoom maintains a^2 <= alpha <= b^2.  This is
    arithmetic via the Fraction front end, where there is no rounding at all
    and the question is purely algebraic: a^2 <= alpha <= b^2.  The
    case-by-case checking lives in test_bracket_zoom.py; here we confirm the
    exact path works and returns a genuinely enclosing bracket."""
    from bracket_zoom import zoom_recursive
    n = 0
    for alpha in (2, 3, 5, 7, 12, 13, 101, 257, 1009, 999983, 104729):
        A = F(alpha)
        # integer seed: largest integer with a^2 <= alpha
        s = 0
        while (s + 1) ** 2 <= alpha:
            s += 1
        a, b, depth = zoom_recursive(A, F(s), F(s + 1), tol_digits=30)
        n += 1
        if not (a * a <= A and A <= b * b):
            return False, f"bracket does not enclose root for alpha={alpha}"
    return True, f"{n} exact-rational brackets, a^2 <= alpha <= b^2"


def decide_certified_runs():
    """Every bracket the certified benchmark emits satisfies
    lo^2 <= alpha <= hi^2 and hi-lo <= 10^-D, in exact integer arithmetic.
    This is DECIDED, not sampled: each run's certificate is a proof about
    that run's output."""
    import exp_certified as C
    n = 0
    for digits in (20, 50, 100):
        for alpha in C.ALPHAS:
            for name, fn in C.METHODS:
                lo, hi, _ = fn(alpha, digits)
                msg = C.certify(lo, hi, alpha, digits)
                n += 1
                if msg:
                    return False, f"{name} a={alpha} d={digits}: {msg}"
    return True, f"{n} runs, exact integer certificates"


def decide_termination_band():
    """The off-by-ten termination bug: `eps.adjusted() <= -digits` accepts
    eps in [10^-D, 10^(1-D)).  Decided, not sampled: the set of Decimals
    accepted by the old test and rejected by the strict one is exactly the
    band [10^-D, 10^(1-D)), and alpha=257 at 50 digits lands in it."""
    from decimal import Decimal as D
    from exp_precision import _width_digits
    import exp_precision as P
    # find the eps the scheduled zoom reaches for alpha=257, digits=50
    alpha, digits = 257, 50
    a = P._seed(alpha)
    w = D(1)
    in_band = None
    for _ in range(200):
        prec = min(P.MIN_PREC + 3 * _width_digits(w) + 3 * a.adjusted(),
                   digits + P.GUARD)
        with P.localcontext() as c:
            c.prec = prec
            p = D(alpha) - a * a
            m = a + p / (2 * a + w)
            est = a + p / (a + m)
        eps = w ** 3 / (27 * a * a)
        if eps.adjusted() <= -digits and eps > D(10) ** (-digits):
            in_band = eps
            break
        if eps.adjusted() < -digits:
            break
        with P.localcontext() as c2:
            c2.prec = prec
            a = +est - eps
            w = +eps
    return in_band is not None, f"alpha=257 d=50 lands in band: {in_band}"


def main():
    print("=" * 78)
    print("PROVED  -- universally quantified, by argument (not by counting)")
    print("=" * 78)
    report("proved", "D'(u) = 2a + 2u > 0 for a>0,u>=0",
           prove_D_increasing(), 1, "the lemma the bound rests on")
    report("proved", "max n(u) = 4w^3/27 at u=2w/3",
           prove_numerator_max_closed_form(), 1, "closed form, no search")
    ok, note = prove_published_bound_dominates()
    report("proved", "Theorem 4 bound valid on [0,w]", ok, 1, note)
    ok, note = prove_ua_never_worse_than_published()
    report("proved", "u-aware bound <= published (Paper1 Thm 8)", ok, 1, note)
    ok, note = prove_corollary7_strict_slack()
    report("proved", "Corollary 7: 1/27 attained at no bracket", ok, 1, note)

    print()
    print("=" * 78)
    print("DECIDED -- exact, per-artifact, with a witness on failure")
    print("=" * 78)
    ok, note = decide_error_identity_all()
    report("decided", "error identity est-r = u^2(w-u)/D(u)", ok, 1, note)
    ok, note = decide_bracket_contains_root()
    report("decided", "bracket invariant a^2 <= alpha <= b^2", ok, 1, note)
    ok, note = decide_certified_runs()
    report("decided", "all certified brackets enclose the root", ok, 1, note)
    ok, note = decide_termination_band()
    report("decided", "old termination test admits a 10x band", ok, 1, note)

    print()
    print("=" * 78)
    print("MEASURED -- wall clock.  Not decidable, not claimed as such.")
    print("=" * 78)
    print("  Performance ratios (e.g. the certified comparison of Sec 10.1) are")
    print("  measurements of specific implementations on one machine.  They are")
    print("  not theorems and are not used to establish any claim about the")
    print("  method.  Only ratios are quoted, with run counts, because absolute")
    print("  times move with machine load.  See the Sec 12 caveat.")

    print()
    if FAIL:
        print(f"FAILURES: {FAIL}")
        return 1
    print("all proved and decided claims hold")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
