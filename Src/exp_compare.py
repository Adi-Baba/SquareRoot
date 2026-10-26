"""Why does the zoom lose to Newton, and can it be made to win?

Three variants, identical scheduler, identical certification, cheap seed:

  cubic   corrected estimate a + p/(a+m);  E <- 3E + log10(16 a^2)
  plain   plain estimate m = a + p/(a+b);  E <- 2E + log10(4(2a+w))
  newton  x <- (x + alpha/x)/2

Measures precision-work (sum of working precisions) separately from wall time,
because they fail for different reasons.
"""
import time
from decimal import Decimal as D, localcontext

ALPHAS = [2, 3, 5, 7, 10, 12, 13, 17, 23, 31, 101, 257, 1009, 7919, 104729, 999983]
GUARD = 25
MIN_PREC = 25


def _width_digits(w):
    return 1 - w.adjusted()


def _seed(alpha, digits=None):
    """floor(sqrt(alpha)) as an exact integer.

    Precision comes from alpha alone (alpha.adjusted() + 5), never from the
    target `digits`.  Justification is in exp_precision._seed: for
    alpha = n^2 the sqrt is exact, otherwise sqrt(alpha) - n >= 1/(2 sqrt + 1),
    so prec > log10(2*alpha) makes the floor unambiguous.

    `digits` is accepted only so call sites written for the old signature
    keep working; it is deliberately ignored.
    """
    value = D(alpha)
    with localcontext() as ctx:
        ctx.prec = max(30, value.adjusted() + 5)
        return +value.sqrt().to_integral_value(rounding="ROUND_FLOOR")


def _seed_naive(alpha, digits):
    """The pre-fix seeding: a full-target-precision sqrt, then floored.

    Kept ONLY so the cost of the original bug can still be measured -- see
    the 'seed cost' section of main().  Nothing in the timed paths uses it.
    """
    with localcontext() as ctx:
        ctx.prec = digits + GUARD + 10
        return +D(alpha).sqrt().to_integral_value(rounding="ROUND_FLOOR")


def zoom(alpha, digits, variant="cubic", cheap_seed=True):
    """One recursion. Returns (estimate, depth, sum_of_precisions, n_div, n_mul)."""
    with localcontext() as ctx:
        ctx.prec = digits + GUARD + 10
        a = _seed(alpha) if cheap_seed else _seed_naive(alpha, digits)
        w = D(1)
        depth = 0
        prec_sum = ndiv = nmul = 0
        while True:
            E = _width_digits(w)
            prec = min(MIN_PREC + 3 * E + 3 * a.adjusted(), digits + GUARD) \
                if variant == "cubic" else \
                min(MIN_PREC + 2 * E + 3 * a.adjusted(), digits + GUARD)
            with localcontext() as c2:
                c2.prec = prec
                prec_sum += prec
                p = alpha - a * a                      # 1 mul
                nmul += 1
                if variant == "cubic":
                    m = a + p / (2 * a + w)            # 1 div
                    ndiv += 1
                    est = a + p / (a + m)              # 1 div
                    ndiv += 1
                    eps = w ** 3 / (27 * a * a)        # 2 mul + 1 div
                    ndiv += 1
                    nmul += 3
                else:                                  # plain
                    est = a + p / (2 * a + w)          # 1 div
                    ndiv += 1
                    eps = w ** 2 / (4 * (2 * a + w))   # 1 mul + 1 div
                    ndiv += 1
                    nmul += 1
            if eps.adjusted() < -digits:
                with localcontext() as c3:
                    c3.prec = digits + GUARD + 10
                    prec_sum += digits + GUARD + 10
                    p = alpha - a * a
                    if variant == "cubic":
                        m = a + p / (2 * a + w)
                        est = a + p / (a + m)
                    else:
                        est = a + p / (2 * a + w)
                    return +est, depth + 1, prec_sum, ndiv, nmul
            with localcontext() as c4:
                c4.prec = prec
                if variant == "cubic":
                    a = +est - eps        # est > r: left endpoint drops below r
                else:
                    a = +est              # est < r: est IS the left endpoint
                w = +eps
            depth += 1


def newton(alpha, digits):
    with localcontext() as ctx:
        ctx.prec = 30
        a = D(alpha)
        x = a.sqrt()
        c = 30
        guard = 20
        prec_sum = ndiv = 0
        for _ in range(200):
            c *= 2
            prec = 2 * c + guard
            with localcontext() as c2:
                c2.prec = prec
                prec_sum += prec
                ndiv += 1
                x = (x + a / x) / 2
            if c >= digits + guard:
                with localcontext() as c3:
                    c3.prec = prec
                    return +x, prec_sum, ndiv


def reference(alpha, digits):
    with localcontext() as ctx:
        ctx.prec = digits + 40
        return D(alpha).sqrt()


def check(est, alpha, digits):
    with localcontext() as ctx:
        ctx.prec = digits + 40
        return abs(est - reference(alpha, digits)) < D(10) ** (-digits)


def bench(fn, reps=3):
    best = 1e9
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best * 1e3


def main():
    DIGITS = [200, 1000, 3000]
    print(f"{'digits':>7} | {'variant':>7} | {'ms':>8} | {'depth':>5} | "
          f"{'sum prec':>9} | {'div':>4} | {'mul':>4} | ok")
    print("-" * 82)
    agg = {}
    for dg in DIGITS:
        for var in ("cubic", "plain"):
            ts = td = tp = nd = nm = ok = 0
            for al in ALPHAS:
                est, dep, sp, dv, ml = zoom(al, dg, var)
                if not check(est, al, dg):
                    print(f"  !! {var} alpha={al} d={dg} FAILED")
                    continue
                ok += 1; td += dep; tp += sp; nd += dv; nm += ml
                ts += bench(lambda: zoom(al, dg, var))
            n = len(ALPHAS)
            agg[(dg, var)] = ts / n
            print(f"{dg:>7} | {var:>7} | {ts/n:8.3f} | {td/n:5.1f} | "
                  f"{tp/n:9.0f} | {nd/n:4.1f} | {nm/n:4.1f} | {ok}/{n}")
        ts = tp = nd = ok = 0
        for al in ALPHAS:
            est, sp, dv = newton(al, dg)
            if not check(est, al, dg):
                print(f"  !! newton alpha={al} d={dg} FAILED"); continue
            ok += 1; tp += sp; nd += dv
            ts += bench(lambda: newton(al, dg))
        n = len(ALPHAS)
        agg[(dg, "newton")] = ts / n
        print(f"{dg:>7} | {'newton':>7} | {ts/n:8.3f} | {'':>5} | "
              f"{tp/n:9.0f} | {nd/n:4.1f} | {'':>4} | {ok}/{n}")
        print("-" * 82)

    print("\nratio to newton (wall time):")
    for dg in DIGITS:
        nv = agg[(dg, "newton")]
        print(f"  {dg:>5} digits: cubic {agg[(dg,'cubic')]/nv:5.2f}x   "
              f"plain {agg[(dg,'plain')]/nv:5.2f}x")

    # how much is the naive seed alone?
    print("\nseed cost, 3000 digits, alpha=12:")
    with localcontext() as ctx:
        ctx.prec = 3000 + GUARD + 10
        t0 = time.perf_counter(); _seed_naive(12, 3000); t_naive = time.perf_counter() - t0
        t0 = time.perf_counter(); _seed(12); t_cheap = time.perf_counter() - t0
    full = bench(lambda: zoom(12, 3000, "cubic", cheap_seed=False))
    cheap = bench(lambda: zoom(12, 3000, "cubic", cheap_seed=True))
    print(f"  naive seed {t_naive*1e3:.3f} ms   cheap seed {t_cheap*1e6:.1f} us")
    print(f"  zoom total: naive seed {full:.2f} ms   cheap seed {cheap:.2f} ms"
          f"   -> seed was {1 - cheap/full:.0%} of the call")


if __name__ == "__main__":
    main()
