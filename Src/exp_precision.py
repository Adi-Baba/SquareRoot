"""
Precision scheduling, plus the shared reference implementations.

This module does double duty, deliberately: it is both the Sec 10 harness and
the library that `exp_certified`, `exp_hybrid` and `exp_mirror` import.  The
shared pieces are `_seed`, `_width_digits`, `_level`, `zoom`, `newton`,
`reference` and `verify`; the Sec 10-only parts are `bench` and `main`.

Hypothesis under test (Paper1.md Sec 10):

    The draft's O(D log D) cost is an artefact of running every level at the
    full target precision D.  If level n is instead evaluated at only the
    precision it needs -- roughly |log10 w| digits of *relative* accuracy on a,
    since the next bracket is a distance ~w^3 away -- then the total work is
    sum_n Theta(3^n) = Theta(D), matching Newton.

RESULT: REFUTED.  Scheduling is worth 2.90x at 3000 digits, not the ~7x the
cost model implies, and Newton remains 3.64x faster at that size.  The theta
classes do match; the constant gap is not a log D artefact.

  (These are the corrected figures.  An earlier version of this docstring
  quoted 1.98x and ~7x, which were measured before the seeding fix below and
  are no longer reproducible.)

Two things were found and fixed while testing it:

  * the initial bracket was computed as a full-target-precision sqrt and then
    floored to an integer -- 58% of every call spent on a four-digit quantity.
    See _seed().  Fixed: 8.72 ms -> 3.63 ms at 3000 digits.
  * the per-level precision budget must be capped at the target; uncapped, the
    schedule's own 3*E_n exceeded D near the end of the recursion and one
    level ran at 7992 digits for a 3000-digit target -- more than the naive
    path, which inverted the comparison under test.

Every result is verified against an exact reference before being timed.

NOTE ON `zoom` HERE: this `zoom(alpha, digits)` returns (estimate, depth) --
an uncertified point.  The *certified* zoom in `exp_certified` shares the name
but returns a rigorous bracket (lo, hi, steps).  They are different functions;
callers disambiguate by module.
"""

import time
from decimal import Decimal as D, localcontext, MAX_PREC

GUARD = 15
MIN_PREC = 30
SEED_PREC = 30


# --------------------------------------------------------------- helpers

def _width_digits(w):
    """Number of significant digits needed to resolve a quantity of size w.

    Decimal.adjusted() is the base-10 exponent of the leading digit, so a
    positive value w has about 1 - w.adjusted() significant digits.
    """
    return 1 - w.adjusted()


def _seed(alpha):
    """floor(sqrt(alpha)) as an exact integer.

    Only enough digits to resolve the integer part unambiguously are needed --
    a quantity that depends on alpha, never on the target precision.  This
    originally ran at the target precision and then floored the result, i.e.
    a full D-digit sqrt to produce a four-digit integer.  Measured at 3000
    digits: 5.05 ms of an 8.72 ms call, 58% of the runtime, for a value that
    costs 6.8 us here -- worth 2.4x of the entire observed gap to Newton.

    Justification.  If alpha = n^2 the sqrt is exact.  Otherwise
    sqrt(alpha) - n = (alpha - n^2)/(sqrt(alpha) + n) >= 1/(2*sqrt(alpha)+1),
    so once sqrt(alpha)*10^-prec < 1/(2*sqrt(alpha)+1) the correctly-rounded
    sqrt cannot fall on the wrong side of n, which holds once
    prec > log10(2*alpha).  alpha.adjusted() + 5 clears that with margin.
    """
    value = D(alpha)
    with localcontext() as ctx:
        ctx.prec = max(MIN_PREC, value.adjusted() + 5)
        return +value.sqrt().to_integral_value(rounding="ROUND_FLOOR")


def _level(a, w, alpha, prec):
    """One zoom level at working precision `prec`. Returns the corrected
    estimate and the plain estimate."""
    with localcontext() as ctx:
        ctx.prec = prec
        p = alpha - a * a
        m = a + p / (2 * a + w)
        est = a + p / (a + m)
        return +est, +m


# --------------------------------------------------------------- zoom

def zoom(alpha, digits, scheduled=True, constant="tight", max_depth=200):
    """Recursive cubic zoom (Algorithm 3).

    scheduled=False -> every level at `digits + GUARD` (the draft's approach)
    scheduled=True  -> level n at just enough precision for the NEXT bracket
    """
    with localcontext() as ctx:
        ctx.prec = digits + GUARD + 10
        a = _seed(alpha)
        w = D(1)
        depth = 0
        while depth < max_depth:
            # A level must resolve `a` to much better than the width it will
            # produce next.  The endpoint a_{n+1} = est_n - eps_n inherits
            # est_n's relative error 10^-prec_n, i.e. an ABSOLUTE error of
            # about r*10^-prec_n, and it must beat w_{n+1} ~ 10^-E_{n+1}.
            # With E_{n+1} = 3E_n + 2*log10(r) this needs prec_n ~ 3E_n + 3*log10(r).
            if scheduled:
                # Cap at the target: once the bracket is already narrow enough
                # to certify `digits`, resolving further is work the naive path
                # never does either.  Without this cap the schedule computed
                # 3*E_n digits even when 3*E_n > D, so a level ran at 7992
                # digits for a 3000-digit target -- more precision than the
                # naive version's 3015 -- which understated scheduling's gain
                # and biased the comparison against it.
                prec = min(MIN_PREC + 3 * _width_digits(w) + 3 * a.adjusted(),
                           digits + GUARD)
            else:
                prec = digits + GUARD
            est, _m = _level(a, w, alpha, prec)
            eps = (w ** 3 / (27 * a * a)) if constant == "tight" else (w ** 3 / (16 * a * a))
            if eps.adjusted() < -digits:
                # The bracket is narrow enough, but `est` was only computed to
                # this level's precision.  Recompute it at full target
                # precision -- one extra level, not one per level.
                est, _ = _level(a, w, alpha, digits + GUARD + 10)
                return est, depth + 1
            # store the bracket at a precision that will survive to the next
            # level; never below what the next level itself will demand
            store = max(prec, MIN_PREC + 3 * _width_digits(eps) + 3 * a.adjusted())
            with localcontext() as ctx2:
                ctx2.prec = store
                a = +est - eps
                w = +eps
            depth += 1
    raise RuntimeError("depth exceeded")


# --------------------------------------------------------------- newton

def newton(alpha, digits, max_iter=200):
    """Newton with doubling precision -- the fair adaptive baseline.

    A well-implemented Newton needs total Theta(M(D)) work, so it is the right
    comparison.  Two implementation details decide whether that is achieved:

    * the starting point must already have several correct digits.  Starting
      from floor(sqrt(a)) gives only ~3, so the doubling schedule spends its
      early (cheap) steps on almost nothing and reaches the target at a far
      higher working precision than necessary.  Here we seed from a float,
      which is essentially free.
    * the stopping test must use the residual |x^2-a| = 2r|x-r| directly rather
      than demanding 10^(-2*digits), which would double the working precision.
    """
    with localcontext() as ctx:
        ctx.prec = SEED_PREC
        a = D(alpha)
        # Seed with Decimal's own correctly-rounded sqrt at SEED_PREC digits.
        # Seeding from the integer part gives only ~2 correct digits, and since
        # Newton squares the error the schedule then needs log2(D/2) steps and
        # overshoots the working precision.  A fair baseline starts near its
        # final accuracy.
        x = a.sqrt()
        c = SEED_PREC                        # correct digits of x, tracked exactly
        guard = 20
        for it in range(1, max_iter + 1):
            # Forming x - a/x cancels catastrophically: the difference is
            # ~2r*10^-c while the operands are ~r, so 2c digits must be carried.
            # Measuring the step at a precision beyond that returns a spurious
            # zero, so the stopping rule tracks c rather than testing the step.
            c = 2 * c
            prec = 2 * c + guard
            with localcontext() as ctx2:
                ctx2.prec = prec
                x = (x + a / x) / 2
            if c >= digits + guard:
                with localcontext() as ctx3:
                    ctx3.prec = prec
                    return +x, it
    raise RuntimeError("no convergence")


# --------------------------------------------------------------- verification

def reference(alpha, digits):
    with localcontext() as ctx:
        ctx.prec = digits + 40
        return D(alpha).sqrt()


def verify(alpha, digits):
    ref = reference(alpha, digits)
    limit = D(10) ** (-digits)
    for name, fn in [
        ("naive", lambda: zoom(alpha, digits, False)),
        ("sched", lambda: zoom(alpha, digits, True)),
        ("newton", lambda: newton(alpha, digits)),
    ]:
        est, _ = fn()
        with localcontext() as ctx:
            ctx.prec = digits + 40
            if abs(est - ref) >= limit:
                return f"{name}: err={abs(est-ref):.3e} >= 1e-{digits}"
    return None


# --------------------------------------------------------------- timing

ALPHAS = [2, 3, 5, 7, 10, 12, 13, 17, 23, 31, 101, 257, 1009, 7919, 104729, 999983]
DIGITS = [50, 200, 1000, 3000]


def bench(fn, reps):
    best, out = float("inf"), None
    for _ in range(reps):
        t0 = time.perf_counter()
        out = fn()
        best = min(best, time.perf_counter() - t0)
    return best * 1e3, out


def main():
    print(f"max Decimal precision supported: {MAX_PREC}\n")
    print(f"{'digits':>7} | {'naive':>9} | {'sched':>9} | {'newton':>9} | "
          f"{'sched/naive':>12} | {'sched/newton':>12} | depth n/s/N")
    print("-" * 92)
    for dg in DIGITS:
        tn = ts = tN = 0.0
        dn = ds = dN = 0
        for al in ALPHAS:
            bad = verify(al, dg)
            if bad:
                print(f"  !! alpha={al} digits={dg}: {bad}")
                continue
            reps = 1 if dg >= 3000 else (2 if dg >= 1000 else 3)
            v1, (e1, d1) = bench(lambda: zoom(al, dg, False), reps)
            v2, (e2, d2) = bench(lambda: zoom(al, dg, True), reps)
            v3, (e3, d3) = bench(lambda: newton(al, dg), reps)
            tn += v1; ts += v2; tN += v3
            dn += d1; ds += d2; dN += d3
        n = len(ALPHAS)
        print(f"{dg:>7} | {tn/n:>9.3f} | {ts/n:>9.3f} | {tN/n:>9.3f} | "
              f"{tn/ts:>11.2f}x | {ts/tN:>11.2f}x | {dn/n:.1f} / {ds/n:.1f} / {dN/n:.1f}")
    print("\nms per sqrt, mean over", len(ALPHAS), "values of alpha")


if __name__ == "__main__":
    main()
