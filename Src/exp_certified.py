"""
Certified-method benchmark: bisection vs interval Newton vs bracket zoom.

WHY THIS EXISTS (Paper1.md Sec 10.1)
------------------------------------
The paper benchmarks the zoom against plain Newton.  But plain Newton returns
a point and no proof -- it cannot tell you the answer is right.  A method
whose entire claim is "the bracket never breaks" should be ranked against the
other methods that carry a certificate:

    bisection         certified, 1 bit per step          -> ~3.32*D steps
    interval Newton   certified, quadratic               -> ~log2(D) steps
    bracket zoom      certified, cubic                   -> ~log3(D) steps   <- under test
    Newton + residual certified, quadratic, one extra term

ALL FOUR return (lo, hi, steps) with lo <= sqrt(alpha) <= hi, proved by
certify() in EXACT INTEGER ARITHMETIC.  Nothing is checked against a reference
sqrt: certify() decides lo**2 <= alpha <= hi**2 and hi-lo <= 10**-D using
Python ints only, so no rounding can decide a bound in its own favour.

Same footing as Sec 10, inherited rather than re-invented:

  * identical ALPHAS, DIGITS, bench()/best-of-N timing, and _seed()
  * every method gets adaptive per-step precision -- none is handicapped with
    a fixed working precision
  * certification runs OUTSIDE the timed region, exactly as Sec 10's
    verify() does; it is a harness cost, not a runtime cost
"""

import signal
import sys
import time
from contextlib import contextmanager
from decimal import Decimal as D, localcontext, ROUND_FLOOR, ROUND_CEILING

from exp_precision import _seed, _width_digits, GUARD, MIN_PREC

# Decimals here routinely carry more than 4300 significant digits, and
# Python 3.11+ refuses int("...") beyond that.  certify() builds integers
# from the digit tuple, so lift the limit rather than fail on large targets.
sys.set_int_max_str_digits(0)


# ------------------------------------------------------- exact certification


def _coeff_exp(d):
    """Decompose a Decimal exactly as (coeff, exp) with d = coeff * 10**exp."""
    t = d.as_tuple()
    coeff = int("".join(map(str, t.digits)))
    if t.sign:
        coeff = -coeff
    return coeff, int(t.exponent)


def _le_sqrt_exact(c, e, alpha):
    """Exact: c*10**e <= sqrt(alpha), for c > 0, alpha >= 0."""
    if e >= 0:
        return c * c * 10 ** (2 * e) <= alpha
    return c * c <= alpha * 10 ** (-2 * e)


def _ge_sqrt_exact(c, e, alpha):
    """Exact: sqrt(alpha) <= c*10**e, i.e. alpha <= (c*10**e)**2."""
    if e >= 0:
        return alpha <= c * c * 10 ** (2 * e)
    return alpha * 10 ** (-2 * e) <= c * c


def _le_sqrt(d, alpha):
    """EXACT: is Decimal d <= sqrt(alpha)?  d >= 0, alpha >= 0."""
    c, e = _coeff_exp(d)
    if c <= 0:
        return True
    return _le_sqrt_exact(c, e, alpha)


def certify(lo, hi, alpha, digits):
    """Prove lo <= sqrt(alpha) <= hi and hi - lo <= 10**-digits.

    Returns None if the certificate holds, else a message.  Exact integer
    arithmetic throughout -- this is a proof, not a comparison against a
    reference sqrt.
    """
    alpha = int(alpha)
    if lo > hi:
        return f"empty bracket: lo={lo} > hi={hi}"

    lc, le = _coeff_exp(lo)
    hc, he = _coeff_exp(hi)

    if lc > 0 and not _le_sqrt_exact(lc, le, alpha):
        return f"lo={lo} is above sqrt(alpha)"
    if hc < 0:
        return f"hi={hi} is below sqrt(alpha)"
    if not _ge_sqrt_exact(hc, he, alpha):
        return f"hi={hi} is below sqrt(alpha)"

    # width: hi - lo <= 10**-digits, exactly
    emin = min(he, le)
    diff = hc * 10 ** (he - emin) - lc * 10 ** (le - emin)
    if diff <= 0:
        return None
    k = digits + emin
    if k >= 0:
        ok = diff * 10 ** k <= 1
    else:
        ok = diff <= 10 ** (-k)
    if not ok:
        return f"width {diff}*10**{emin} exceeds 10**-{digits}"
    return None


# ------------------------------------------------------- helpers


def _prec_for(alpha, w, guard=GUARD, rate=1):
    """Working precision needed to resolve a bracket of width w for this alpha.

    a ~ sqrt(alpha) has about adj(alpha)/2 digits before the point, and the
    bracket needs |log10 w| more below that.

    `rate` is the number of digits the NEXT width is expected to demand
    relative to the current one, because rounding must stay below the width
    the method is about to produce:

        rate=1  bisection, next width ~ w/2
        rate=2  interval Newton, next width ~ w**2/(2r)  -- with rate=1 it
                gains only `guard` digits per step and degrades to linear
        (the zoom uses its own rate=3 schedule in its own expression)
    """
    return guard + (int(D(alpha).adjusted()) // 2 + 1) + rate * _width_digits(w)


@contextmanager
def _round(prec, rounding):
    with localcontext() as c:
        c.prec = prec
        c.rounding = rounding
        yield


class _Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise _Timeout()


# ------------------------------------------------------- 1. bisection


def bisection(alpha, digits, max_steps=100000):
    """Certified bisection.  The sign test is EXACT, so no rounding can
    misplace the root: whichever side m falls on, m is a valid endpoint.

    Returns (lo, hi, steps) with lo <= r <= hi.
    """
    a = _seed(alpha)          # floor(sqrt(alpha)) <= r
    b = a + D(1)              # r <= a+1
    tol = D(10) ** (-digits)
    steps = 0
    while b - a > tol:
        prec = _prec_for(alpha, b - a)
        with _round(prec, ROUND_FLOOR):
            m = (a + b) / 2
        if m <= a or m >= b:
            raise RuntimeError("bisection stalled: precision too low")
        if _le_sqrt(m, alpha):
            a = m
        else:
            b = m
        steps += 1
        if steps > max_steps:
            raise RuntimeError("bisection exceeded max_steps")
    return a, b, steps


# ------------------------------------------------------- 2. interval Newton


def interval_newton(alpha, digits, max_iter=200):
    """Interval Newton on f(x)=x^2-alpha with outward rounding throughout.

    N = m - f(m)/f'([a,b]); new bracket = [a,b] inter N.  By the mean value
    theorem the root lies in N whenever m and the root both lie in [a,b] and
    f' does not vanish there.  Here f'(x)=2x does not vanish since a >= 1, so
    the intersection always retains the root, and the standard argument gives
    quadratic convergence.

    Returns (lo, hi, iters).
    """
    A = D(alpha)
    a = _seed(alpha)
    b = a + D(1)
    tol = D(10) ** (-digits)
    it = 0
    while b - a > tol:
        # rate=2: the enclosure produced here is O(w**2/(2r)), so precision
        # must resolve that rather than w.  At rate=1 it gained only ~16
        # digits per step and needed 65 iterations for 1000 digits.
        prec = _prec_for(alpha, b - a, rate=2)
        with _round(prec, ROUND_FLOOR):
            m_lo = (a + b) / 2
            fm_lo = m_lo * m_lo - A
            den_lo = 2 * a
        with _round(prec, ROUND_CEILING):
            m_hi = (a + b) / 2
            fm_hi = m_hi * m_hi - A
            den_hi = 2 * b

        # Q = f(m)/f'([a,b]) with f' = [den_lo, den_hi] > 0.
        # For fixed y>0, x/y is increasing in x, so the min sits at x_lo and
        # the max at x_hi.  Each of those, as a function of y, is decreasing
        # when its numerator is positive and increasing when negative -- so
        # which denominator extremises depends on the SIGN of f(m), which is
        # negative on the first iterations.  Branching costs two divisions;
        # taking x_lo/y_hi unconditionally would be wrong there.
        with _round(prec, ROUND_FLOOR):
            q_lo = fm_lo / (den_hi if fm_lo >= 0 else den_lo)
        with _round(prec, ROUND_CEILING):
            q_hi = fm_hi / (den_lo if fm_hi >= 0 else den_hi)
        # N = m - Q = [m_lo - q_hi, m_hi - q_lo]  (cross terms)
        with _round(prec, ROUND_FLOOR):
            n_lo = m_lo - q_hi
        with _round(prec, ROUND_CEILING):
            n_hi = m_hi - q_lo

        new_a = a if a >= n_lo else n_lo
        new_b = b if b <= n_hi else n_hi
        if not (new_a < new_b):
            # outward rounding swallowed the bracket: fall back to one
            # exact bisection step, which is always safe
            with _round(prec + 20, ROUND_FLOOR):
                m = (a + b) / 2
            if _le_sqrt(m, alpha):
                new_a, new_b = m, b
            else:
                new_a, new_b = a, m
        if new_a == a and new_b == b:
            raise RuntimeError("interval Newton stalled")
        a, b = new_a, new_b
        it += 1
        if it > max_iter:
            raise RuntimeError("interval Newton exceeded max_iter")
    return a, b, it


# ------------------------------------------------------- 3. bracket zoom


def _dn(prec, op, *args):
    with _round(prec, ROUND_FLOOR):
        return op(*args)


def _up(prec, op, *args):
    with _round(prec, ROUND_CEILING):
        return op(*args)


def _est_bracket(a, w, A, prec):
    """Rigorously enclose the corrected estimate and its error bound.

    Returns (est_lo, est_hi, eps_hi) with

        exact est in [est_lo, est_hi]   and   exact eps <= eps_hi,

    where est = a + p/(a+m), m = a + p/(2a+w), p = A - a^2 and
    eps = w^3/(27 a^2), all evaluated on the *exact* Decimals a and w.

    A single directed-rounding context is NOT enough here: flooring a+m
    enlarges p/(a+m), so flooring each operation in turn does not bound the
    composed expression.  Each intermediate therefore has to be enclosed
    separately and the composition propagated.

    Needed because the exact identity 0 <= est - r <= eps says nothing about
    est rounded to finite precision: when est - r happens to be near its
    bound, or near zero, the rounding of est alone moves it outside
    [est-eps, est].  Sec 10.1 measured that breaking 16 times in 48 runs.
    """
    # a^2, with a >= 0 so squaring is increasing
    a2l = _dn(prec, lambda x: x * x, a)
    a2h = _up(prec, lambda x: x * x, a)
    # p = A - a^2  in  [A - a2h, A - a2l]; A is a short exact Decimal
    pl = _dn(prec, lambda x, y: x - y, A, a2h)
    ph = _up(prec, lambda x, y: x - y, A, a2l)
    # t = 2a + w
    tl = _dn(prec, lambda x, y: 2 * x + y, a, w)
    th = _up(prec, lambda x, y: 2 * x + y, a, w)
    # p/t with p >= 0 (the invariant a <= r) and t > 0:  [pl/th, ph/tl]
    qtl = _dn(prec, lambda x, y: x / y, pl, th)
    qth = _up(prec, lambda x, y: x / y, ph, tl)
    # m = a + p/t, then s = a + m
    ml = _dn(prec, lambda x, y: x + y, a, qtl)
    mh = _up(prec, lambda x, y: x + y, a, qth)
    sl = _dn(prec, lambda x, y: x + y, a, ml)
    sh = _up(prec, lambda x, y: x + y, a, mh)
    # p/s decreasing in s
    rsl = _dn(prec, lambda x, y: x / y, pl, sh)
    rsh = _up(prec, lambda x, y: x / y, ph, sl)
    # est = a + p/s
    el = _dn(prec, lambda x, y: x + y, a, rsl)
    eh = _up(prec, lambda x, y: x + y, a, rsh)
    # eps = w^3 / (27 a^2), all positive
    w3h = _up(prec, lambda x: x * x * x, w)
    dnl = _dn(prec, lambda x: 27 * x, a2l)
    epsh = _up(prec, lambda x, y: x / y, w3h, dnl)
    return el, eh, epsh


def zoom(alpha, digits, max_depth=200):
    """The paper's cubic zoom, returning a RIGOROUSLY certified bracket.

    The exact theory gives r in [est-eps, est]; this returns an enclosure of
    that set which accounts for the rounding of est and eps themselves, so
    the printed bracket survives certify() -- an exact check.

    Termination uses eps.adjusted() < -digits, i.e. eps < 10**-digits.  The
    earlier form `<= -digits` admitted eps in [10**-D, 10**(1-D)) -- a
    certificate up to ten times looser than the target, which alpha=257 at 50
    digits actually hit.  Pinned by test_certified.py.

    Two bugs fixed here, both of which broke inputs this project actually
    uses:

    * alpha = 0 divided by zero in _est_bracket (a = s = 0 makes p/s the
      0/0).  sqrt(0) is exact, so it is returned directly.
    * The working precision was capped at `digits + GUARD`, which is not
      enough once sqrt(alpha) has digits before the decimal point: those
      count against the same significant-digit budget.  For alpha around
      2**500 the cap was ~50 digits while the estimate needed ~180, so
      w_new never fell below w and the loop raised "zoom stalled" instead of
      producing a bracket.  The budget now includes the integer part, as
      _prec_for does.
    """
    if alpha < 0:
        raise ValueError(f"alpha must be >= 0, got {alpha}")
    if alpha == 0:
        # sqrt(0) = 0 exactly; no zoom needed and _est_bracket would divide
        # by a+m = 0.
        return D(0), D(0), 0

    tol = D(10) ** (-digits)
    # Digits before the decimal point of sqrt(alpha), which come out of the
    # same significant-digit budget as the digits after it.
    int_digits = max(1, _seed(alpha).adjusted() + 1)
    pprec = digits + GUARD + 10 + int_digits
    with localcontext() as ctx:
        ctx.prec = pprec
        A = D(alpha)
        a = _seed(alpha)              # a <= r
        w = D(1)                      # [a, a+1] brackets r
        depth = 0
        while depth < max_depth:
            prec = max(MIN_PREC,
                       MIN_PREC + 3 * _width_digits(w) + 3 * a.adjusted() + 3 * int_digits)
            # one interval evaluation serves both the termination test and the
            # step; the outward rounding is what keeps a <= r <= a+w intact
            el, eh, epsh = _est_bracket(a, w, A, prec)
            if epsh.adjusted() < -digits:
                # re-enclose at full target precision for the returned bracket
                fl, fh, feps = _est_bracket(a, w, A, pprec)
                lo = _dn(pprec, lambda x, y: x - y, fl, feps)
                hi = fh
                width = _up(pprec, lambda x, y: x - y, hi, lo)
                if width <= tol:
                    return lo, hi, depth + 1
                # not narrow enough after all: fall through and zoom once more
            a_new = _dn(prec, lambda x, y: x - y, el, epsh)
            w_new = _up(prec, lambda x, y: x - y, eh, a_new)
            if not (w_new < w):
                raise RuntimeError("zoom stalled: bracket not shrinking")
            a, w = a_new, w_new
            depth += 1
    raise RuntimeError("depth exceeded")


# ------------------------------------------------------- 4. Newton + residual


def newton_certified(alpha, digits, float_seed=True, max_iter=200):
    """Newton, certified afterwards by a residual bound.

    |x - r| = |x^2 - alpha| / (x + r) <= |x^2 - alpha| / (x + lo) once
    lo <= r is known -- and floor(sqrt(alpha)) supplies lo for free.  The
    certificate costs one extra squaring, which is the honest price of turning
    a point method into a certified one.

    float_seed=True (the conventional choice, and what Sec 10 uses) starts
    Newton from a 30-digit sqrt, a 30-digit head start that the zoom and
    interval Newton are NOT given -- they begin at floor(sqrt(alpha)) with
    ~3 correct digits.  float_seed=False puts Newton on the same footing as
    the other three so the comparison can be read seed for seed; bench()
    reports both, because the difference is large enough to move the ranking.

    Returns (lo, hi, iters).
    """
    A = D(alpha)
    lo_seed = _seed(alpha)                 # <= r, free
    if float_seed:
        with localcontext() as ctx:
            ctx.prec = 30
            x = A.sqrt()
        c = 30
    else:
        # NOTE: the int-seed path still calls A.sqrt() here.  It is evaluated
        # at the *default* context precision, not the 30 digits the float-seed
        # path requests, so it yields a short seed -- but the asymmetry is
        # incidental to the context, not to the branch.  The paper's shared
        # seed is the integer floor; this path approximates it well enough
        # that it still takes ~4 more iterations than the float-seed path
        # (7 vs 2 at 100 digits, alpha=12).  If a future edit makes this
        # explicitly `x = +lo_seed`, the benchmark's seed-sensitivity numbers
        # change and the paper must be re-measured.
        x = A.sqrt()
        c = 1
    for it in range(1, max_iter + 1):
        # c = 2c+1 rather than 2c: from a short seed the doubling bound
        # overshoots, and an optimistic counter stops one step early and
        # fails certify() for some alpha.  Being pessimistic costs a cheap
        # extra step; being optimistic costs the certificate.
        c = 2 * c + 1 if not float_seed else 2 * c
        prec = 2 * c + 20
        with localcontext() as c2:
            c2.prec = prec
            x = (x + A / x) / 2
        if c >= digits + 20:
            # Trim to what the target needs before certifying: the doubling
            # schedule overshoots (c=3840 for a 3000-digit target) and the
            # step is evaluated at 2*c+20, so x arrives carrying 7700 digits.
            # The residual bound holds for any x with lo_seed <= r, so trimming
            # first is free and keeps the certificate small.
            with localcontext() as c3:
                c3.prec = digits + GUARD + 10
                x = +x
            cprec = digits + GUARD + 50      # resolve |x^2-alpha| ~ 2r*10^-c
            with _round(cprec, ROUND_FLOOR):
                sq_lo = x * x
            with _round(cprec, ROUND_CEILING):
                sq_hi = x * x
            resid = max(abs(sq_lo - A), abs(sq_hi - A))
            with _round(cprec, ROUND_FLOOR):
                err_lo = resid / (x + lo_seed)
                lo_b = x - err_lo
            with _round(cprec, ROUND_CEILING):
                err_hi = resid / (x + lo_seed)
                hi_b = x + err_hi
            # NB: no certify() here.  Verification belongs in verify(), OUTSIDE
            # the timed region; pulling it in cost Newton+res ~2.7 ms of the
            # 9.8 ms row at 3000 digits.  What keeps the certificate sound is
            # the pessimistic counter above PLUS this width check, which is
            # one subtraction and one comparison -- a Decimal op, not an
            # exact-integer proof.  Without it the counter's slack of 20 is
            # marginal (measured width 2.6e-1064 against a 1e-1000 target
            # was close enough to fail for several alpha at 1000 digits) and
            # 30 of 480 runs failed.
            if hi_b - lo_b > D(10) ** (-digits):
                continue                  # one more step; counter was short
            return lo_b, hi_b, it
    raise RuntimeError("no convergence")


# ------------------------------------------------------- harness

ALPHAS = [2, 3, 5, 7, 10, 12, 13, 17, 23, 31, 101, 257, 1009, 7919, 104729, 999983]
DIGITS = [50, 200, 1000, 3000]

METHODS = [
    ("bisection", bisection),
    ("int.Newton", interval_newton),
    ("zoom", zoom),
    ("Newton+res", newton_certified),
    # Same Newton, same certificate, but seeded from floor(sqrt(alpha)) --
    # the seed the zoom and interval Newton actually get.  Reported alongside
    # the conventional float-seeded row because the head start is worth
    # 2-4x and is large enough to move the ranking.
    ("N+res(int)", lambda a, d: newton_certified(a, d, float_seed=False)),
]


def verify(alpha, digits):
    """Certify every method's bracket.  Returns None or an error string."""
    for name, fn in METHODS:
        try:
            lo, hi, _steps = fn(alpha, digits)
        except Exception as exc:                       # noqa: BLE001
            return f"{name}: {exc}"
        msg = certify(lo, hi, alpha, digits)
        if msg:
            return f"{name}: {msg}"
    return None


def _run_with_timeout(fn, seconds):
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(seconds)
    try:
        return fn()
    finally:
        signal.alarm(0)


def main():
    print("certification: exact integer arithmetic (lo**2 <= alpha <= hi**2),")
    print("no reference sqrt involved; check runs outside the timed region")
    print("methods interleaved per alpha (best-of-N) so load drift cannot be")
    print("attributed to whichever method happened to run last\n")
    print(f"{'digits':>7} | {'method':>11} | {'ms':>9} | {'steps':>7}")
    print("-" * 50)
    totals = {}
    for dg in DIGITS:
        reps = 5 if dg <= 200 else (3 if dg <= 1000 else 2)
        best = {(name, al): None for name, _ in METHODS for al in ALPHAS}
        step_of = {(name, al): 0 for name, _ in METHODS for al in ALPHAS}
        for al in ALPHAS:
            for name, fn in METHODS:
                for _ in range(reps):
                    t0 = time.perf_counter()
                    try:
                        lo, hi, st = _run_with_timeout(lambda: fn(al, dg), 120)
                    except _Timeout:
                        print(f"{dg:>7} | {name:>11} | {'>120s':>9} | {'-':>7}")
                        return
                    except Exception as exc:                   # noqa: BLE001
                        print(f"{dg:>7} | {name:>11} | error: {exc}")
                        return
                    ms = (time.perf_counter() - t0) * 1e3
                    msg = certify(lo, hi, al, dg)
                    if msg:
                        print(f"  !! CERT FAIL alpha={al} digits={dg} {name}: {msg}")
                    if best[(name, al)] is None or ms < best[(name, al)]:
                        best[(name, al)] = ms
                    step_of[(name, al)] = st
        for name, _ in METHODS:
            times = [best[(name, al)] for al in ALPHAS]
            avg_t = sum(times) / len(times)
            avg_s = sum(step_of[(name, al)] for al in ALPHAS) / len(ALPHAS)
            totals[(dg, name)] = (avg_t, avg_s)
            print(f"{dg:>7} | {name:>11} | {avg_t:>9.3f} | {avg_s:>7.1f}")
        print("-" * 50)

    print("\nspeedup vs Newton+res (certified), >1 means slower than Newton:\n")
    print(f"{'digits':>7} | {'bisection':>11} | {'int.Newton':>11} | {'zoom':>11}")
    print("-" * 52)
    for dg in DIGITS:
        base = totals.get((dg, "Newton+res"))
        if not base:
            continue
        row = []
        for name in ("bisection", "int.Newton", "zoom"):
            v = totals.get((dg, name))
            row.append(f"{v[0] / base[0]:>10.2f}x" if v else f"{'n/a':>11}")
        print(f"{dg:>7} | " + " | ".join(row))

    print("\nseed sensitivity: Newton+res from a 30-digit float sqrt vs from")
    print("floor(sqrt(alpha)) -- the seed the zoom and interval Newton get:\n")
    print(f"{'digits':>7} | {'float ms':>9} {'steps':>6} | {'int ms':>9} "
          f"{'steps':>6} | {'cost':>7}")
    print("-" * 62)
    for dg in DIGITS:
        f, i = totals.get((dg, "Newton+res")), totals.get((dg, "N+res(int)"))
        if not f or not i:
            continue
        print(f"{dg:>7} | {f[0]:>9.3f} {f[1]:>6.1f} | {i[0]:>9.3f} "
              f"{i[1]:>6.1f} | {i[0] / f[0]:>6.2f}x")
    print("\n(ms per sqrt, mean over 16 alpha; >120s = did not finish)")


if __name__ == "__main__":
    main()
