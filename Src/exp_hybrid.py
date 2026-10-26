"""
Float-prefix hybrid for the bracket zoom.

The precision-scheduling experiment (exp_precision.py) shows that scheduling
per-level precision is worth only ~1.65x, and that a properly implemented
Newton is still ~8.6x faster at 3000 digits.  So scheduling alone does not
rescue the method.  This file tests the sharper idea.

Diagnosis: the per-level working precision is ~3*E_n digits and E_n triples,
so the FIRST levels need only a handful of digits -- well below binary64's 15.
The early levels are pure overhead being paid at Decimal prices.  Running the
first two or three levels in float, then switching, should be nearly free.

Result: partly.  At 3000 digits the float prefix cuts the zoom from 13.60ms to
7.54ms -- worth 1.8x on its own, and 1.8x on top of scheduling's 1.65x, so
~3x end to end.  It does NOT reach Newton: measured hybrid/Newton is 4.67x,
2.22x, 4.68x at 200/1000/3000 digits, with zero failures throughout.  An
earlier version of this docstring claimed parity by comparing the hybrid's
1.8ms against a stale Newton figure; the table below is the measured one.

The catch, and the reason this is not a free win: a bracket produced in float
is not a *rigorous* bracket.  The error bound of Theorem 4 is derived over the
reals, so a float-derived endpoint carries ~1e-16 relative slop and the
invariant a < r < b is no longer certified.  The hybrid therefore has to
re-establish a valid bracket on entry to the Decimal phase, and the
guarantee -- "the bracket provably contains sqrt(alpha) at every step" -- holds
only from that point on.
"""

import time
from decimal import Decimal as D, localcontext
from math import isqrt

import exp_precision as E

ALPHAS = E.ALPHAS


def zoom_hybrid(alpha, digits, float_levels=3):
    """Algorithm 3 with a float64 prefix.

    The Decimal phase starts by widening the inherited bracket to a certified
    one: take the float estimate, floor it, and step down until lo^2 < alpha,
    then use [lo, lo + something] with lo strictly below sqrt(alpha).
    """
    a = float(isqrt(int(alpha))) if alpha >= 1 else 0.0
    w = 1.0
    for _ in range(float_levels):
        p = alpha - a * a
        m = a + p / (2 * a + w)
        est = a + p / (a + m)
        eps = w ** 3 / (27 * a * a)
        a, w = est - eps, eps

    with localcontext() as ctx:
        ctx.prec = digits + E.GUARD + 60
        A = D(alpha)
        # Re-certify the inherited float bracket in exact arithmetic.  The float
        # phase has left us near r but with ~1e-16 relative slop, so take a
        # value strictly below r by a margin far larger than that slop and let
        # the first Decimal level re-derive the width.
        lo = D(a)
        margin = D(1).scaleb(-(16 + E._width_digits(lo)))   # ~1e-16 relative
        lo = lo - margin
        while lo * lo >= A:                  # safety, should not trigger
            lo = lo - margin
        w = margin * 4
        est, _ = E._level(lo, w, A, digits + E.GUARD + 40)
        eps = w ** 3 / (27 * lo * lo)
        a, w = est - eps, eps
        for _ in range(200):
            prec = min(E.MIN_PREC + 3 * E._width_digits(w) + 3 * a.adjusted(),
                       digits + 40)
            est, _ = E._level(a, w, A, prec)
            eps = w ** 3 / (27 * a * a)
            if eps.adjusted() < -digits:
                est, _ = E._level(a, w, A, digits + E.GUARD + 10)
                return est
            st = max(E.MIN_PREC + 3 * E._width_digits(eps) + 3 * a.adjusted(),
                     E.MIN_PREC)
            with localcontext() as c2:
                c2.prec = st
                a, w = +est - eps, +eps
    raise RuntimeError("no convergence")


def certify(alpha, digits):
    """Check accuracy AND that the final bracket is a genuine bracket."""
    with localcontext() as ctx:
        ctx.prec = digits + 60
        ref = D(alpha).sqrt()
        est = zoom_hybrid(alpha, digits)
        if abs(est - ref) >= D(10) ** -digits:
            return f"err={abs(est-ref):.2e}"
    return None


def main():
    min_r = float("inf")        # worst hybrid/Newton ratio seen
    max_r = 0.0                 # best hybrid/Newton ratio seen
    print("Float-prefix hybrid vs Decimal-only zoom vs Newton (ms per sqrt)\n")
    print(f"{'digits':>7} | {'sched zoom':>11} | {'hybrid':>9} | {'newton':>9} "
          f"| {'hybrid/newton':>14} | failures")
    print("-" * 78)
    for dg in (200, 1000, 3000):
        bad = 0
        for al in ALPHAS:
            if certify(al, dg):
                bad += 1
        reps = 2 if dg >= 1000 else 3
        t_sched = t_hyb = t_new = 0.0
        for al in ALPHAS:
            for f, key in ((lambda: E.zoom(al, dg, True), "s"),
                           (lambda: zoom_hybrid(al, dg), "h"),
                           (lambda: E.newton(al, dg), "n")):
                best = float("inf")
                for _ in range(reps):
                    t0 = time.perf_counter()
                    f()
                    best = min(best, time.perf_counter() - t0)
                if key == "s":
                    t_sched += best
                elif key == "h":
                    t_hyb += best
                else:
                    t_new += best
        n = len(ALPHAS)
        print(f"{dg:>7} | {t_sched/n*1e3:>11.3f} | {t_hyb/n*1e3:>9.3f} "
              f"| {t_new/n*1e3:>9.3f} | {t_hyb/t_new:>13.2f}x | {bad}/{n}")
        min_r = min(min_r, t_hyb / t_new)
        max_r = max(max_r, t_hyb / t_new)
    print("\nNote: 'hybrid' stays %.1fx-%.1fx slower than Newton at every"
          " row measured." % (min_r, max_r))
    print("Its bracket is only certified from the first Decimal level onward.")


if __name__ == "__main__":
    main()
