"""Does Blender-style mirroring help the bracket zoom?

The reflection symmetry of sqrt across r = sqrt(alpha) is NOT additive --
it is multiplicative.  Since r^2 = alpha, r is the geometric mean of any
pair (x, alpha/x):

    x  ---reflect--->  alpha/x
    GM(x, alpha/x) = sqrt(x * alpha/x) = alpha^0.5 = r

so the mirror of a point across the root is alpha/x, not 2r - x.  Three
consequences are tested here:

  A. mirror seed  -- bracket [a, alpha/a] instead of [a, a+1]
  B. narrower     -- take whichever of the two brackets is narrower
  C. free half    -- for a mirror bracket (a*b = alpha) we have r = GM < AM,
                     so r lies in the LEFT half without testing anything.
                     That is a bisection bit for free.

Punchline to check: mirror + average is exactly Newton,
    x <- (x + alpha/x) / 2.
If that is what the idea reduces to, then the idea is the competitor.
"""
import time
from decimal import Decimal as D, localcontext
from exp_precision import _seed, _width_digits, GUARD, MIN_PREC

ALPHAS = [2, 3, 5, 7, 10, 12, 13, 17, 23, 31, 101, 257, 1009, 7919, 104729, 999983]


def zoom_mirror(alpha, digits, mode="plain"):
    """Cubic zoom, differing only in how the INITIAL bracket is chosen.

    mode="plain" -- [a, a+1]                       (what we do now)
    mode="mirror" -- [a, alpha/a]                  (the reflected bracket)
    mode="narrow" -- whichever is narrower         (free: one division)
    """
    with localcontext() as ctx:
        ctx.prec = digits + GUARD + 10
        A = D(alpha)
        a = _seed(alpha)
        # Guard: for a perfect square alpha/a == a exactly, so the mirror and
        # narrow brackets collapse to width 0.  Then p = 0, est = a, eps = 0,
        # and Decimal(0).adjusted() is 0 -- never <= -digits -- so the loop
        # below spins forever.  (plain mode is unaffected: it starts at w=1.)
        # a is already the exact root here, so nothing is lost by returning it.
        if a * a == A:
            return +a, 0
        if mode == "plain":
            w = D(1)
        else:
            # reflected endpoint: alpha / a.  r >= a and r <= alpha/a both hold
            # because a <= r <= alpha/a is equivalent to a^2 <= alpha <= a*r,
            # i.e. a <= r, which is exactly the seeding guarantee.
            inv = A / a
            w = inv - a
            if mode == "narrow" and w > D(1):
                w = D(1)
        depth = 0
        while True:
            prec = min(MIN_PREC + 3 * _width_digits(w) + 3 * a.adjusted(),
                       digits + GUARD)
            with localcontext() as c2:
                c2.prec = prec
                p = A - a * a
                m = a + p / (2 * a + w)
                est = a + p / (a + m)
                eps = w ** 3 / (27 * a * a)
            if eps.adjusted() < -digits:
                with localcontext() as c3:
                    c3.prec = digits + GUARD + 10
                    p = A - a * a
                    m = a + p / (2 * a + w)
                    return + (a + p / (a + m)), depth + 1
            with localcontext() as c4:
                c4.prec = prec
                a, w = +est - eps, +eps
            depth += 1


def zoom_freehalf(alpha, digits):
    """Mirror bracket + its free bisection bit, applied ONLY where valid.

    For a mirror bracket ab = alpha we have r = GM(a,b) < AM(a,b), so r lies
    in [a, AM] and the width halves with no test at all.  That deduction
    depends on ab = alpha, which holds only at seeding -- after the first zoom
    level a is est - eps and the invariant is gone, so the half must NOT be
    taken again.  (An earlier version halved every iteration and passed 0/16
    by luck: the bracket is only safe to halve while r sits below the
    arithmetic mean, which is not guaranteed once est - r < eps/2.)
    """
    with localcontext() as ctx:
        ctx.prec = digits + GUARD + 10
        A = D(alpha)
        a = _seed(alpha)
        # Same perfect-square guard as zoom_mirror: alpha/a == a gives w = 0
        # and an infinite loop below.  a is exact here, so return it.
        if a * a == A:
            return +a, 0
        w = (A / a - a) / 2                 # free half, mirror invariant holds here
        depth = 0
        while True:
            prec = min(MIN_PREC + 3 * _width_digits(w) + 3 * a.adjusted(),
                       digits + GUARD)
            with localcontext() as c2:
                c2.prec = prec
                p = A - a * a
                m = a + p / (2 * a + w)
                est = a + p / (a + m)
                eps = w ** 3 / (27 * a * a)
            if eps.adjusted() < -digits:
                with localcontext() as c3:
                    c3.prec = digits + GUARD + 10
                    p = A - a * a
                    m = a + p / (2 * a + w)
                    return + (a + p / (a + m)), depth + 1
            with localcontext() as c4:
                c4.prec = prec
                a, w = +est - eps, +eps
            depth += 1


def check(est, alpha, digits):
    with localcontext() as ctx:
        ctx.prec = digits + 40
        return abs(est - D(alpha).sqrt()) < D(10) ** (-digits)


def bench(fn, reps=3):
    best = 1e9
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best * 1e3


def main():
    print("=== A/B/C: seed strategy, 3000 digits ===\n")
    modes = [("plain", lambda a, d: zoom_mirror(a, d, "plain")),
             ("mirror", lambda a, d: zoom_mirror(a, d, "mirror")),
             ("narrow", lambda a, d: zoom_mirror(a, d, "narrow")),
             ("freehalf", zoom_freehalf)]
    for digits in (1000, 3000):
        print(f"--- {digits} digits ---")
        print(f"{'mode':>10} | {'ms':>7} | {'depth':>5} | {'fail':>5} | vs plain")
        base = None
        for name, fn in modes:
            tot = dep = bad = 0
            for al in ALPHAS:
                est, d = fn(al, digits)
                if not check(est, al, digits):
                    bad += 1
                tot += bench(lambda: fn(al, digits))
                dep += d
            n = len(ALPHAS)
            ms = tot / n
            if base is None:
                base = ms
            print(f"{name:>10} | {ms:7.3f} | {dep/n:5.1f} | {bad:>2}/{n} | {base/ms:5.2f}x")
        print()

    # does the mirror bracket actually start narrower?
    print("=== mirror bracket width vs width-1, by alpha ===")
    with localcontext() as ctx:
        ctx.prec = 40
        narrower = same = wider = 0
        for al in ALPHAS:
            a = _seed(al)
            w_mirror = D(al) / a - a
            verdict = "narrower" if w_mirror < 1 else ("same" if w_mirror == 1 else "wider")
            if verdict == "narrower":
                narrower += 1
            elif verdict == "same":
                same += 1
            else:
                wider += 1
        print(f"  mirror narrower: {narrower}   equal: {same}   wider: {wider}"
              f"   (of {len(ALPHAS)})")

    # the punchline: is mirror + average exactly newton?
    print("\n=== is (x + alpha/x)/2 Newton? ===")
    with localcontext() as ctx:
        ctx.prec = 60
        al, x = 12, D(3)
        for i in range(4):
            x = (x + D(al) / x) / 2
            print(f"  step {i+1}: x = {x:.40f}")
        print(f"  true sqrt   = {D(al).sqrt():.40f}")


if __name__ == "__main__":
    main()
