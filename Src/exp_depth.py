"""
Does the u-aware bound reduce recursion depth, and by how much?

Tooling note: exact Fractions were the wrong instrument here.  They are right
for PROVING a bound (the quantities stay small), but for COUNTING iterations
each level cubes a rational, so integers grow without bound and pure-Python
integer arithmetic swamps the run.  The recursion needs only the order of
magnitude of each bound, so the right instrument is mpmath at fixed working
precision.

Both bounds maximise the SAME numerator u^2(w-u) and divide by the SAME
denominator lower bound 2a(2a+w); they differ only in the interval over which
the maximisation runs:

    published : u in (0, w)                -- nothing known about u
    u-aware   : u in (m-a, est-a)          -- known after one evaluation

so bound_published must divide by 2a(2a+w) too, or the comparison is rigged.

Convergence signal: p = alpha - a^2 underflows to exactly 0 once `a` becomes
the correctly-rounded mpf root, which makes b == 0.  That is SUCCESS -- we have
reached working precision, which the caller set to digits + guard -- not
failure.  Treating b <= 0 as an error sentinel (as an earlier draft did)
reports success as failure and poisons the median with -1.
"""

import math
import statistics as st

from mpmath import mp, mpf

SAMPLE = [2, 3, 5, 7, 10, 12, 13, 17, 23, 31, 101, 257, 1009, 7919, 104729, 999983]
MAX_LEVELS = 200
FAILURE = -1
LOG10_TO_LOG2 = 3.3219280948873623


def _p(a, alpha):
    """alpha - a^2.  Uses alpha rather than r^2: at level 0 alpha is an exact
    integer, so this avoids importing rounding error from the oracle sqrt into
    the very first evaluation."""
    return alpha - a * a


def corrected_from_m(a, m, p):
    """est from already-computed m and p.  A naive corrected(a,w,alpha) would
    call plain() again, costing a second full-precision a^2 per level -- that
    was the dominant term, and doubling it doubled the run time at 10000
    digits."""
    return a + p / (a + m)


def _num_interval_max(w, lo, hi):
    """max of u^2(w-u) over the closed interval [lo, hi]; critical point 2w/3."""
    cand = [lo, hi]
    c = 2 * w / 3
    if lo < c < hi:
        cand.append(c)
    return max(u * u * (w - u) for u in cand)


def bound_published(a, w):
    """Sharp bound using only the assumption u in (0, w): maximise the
    numerator over the whole interval, divide by the same denominator bound the
    u-aware version uses."""
    return _num_interval_max(w, 0, w) / (2 * a * (2 * a + w))


def bound_known_interval(a, w, m, est):
    """The same expression restricted to the interval we actually know."""
    return _num_interval_max(w, m - a, est - a) / (2 * a * (2 * a + w))


def _reached(b, digits):
    """True iff b <= 10^-digits.

    An earlier version derived this from mp.mag() to avoid mp.log10(), which at
    the ambient 30000-digit working precision is a transcendental evaluation
    and dominated the whole script.  But mpmath's mag() returns e for a
    mantissa in [0.5, 1), not [1, 2), so the floor(log2) inequalities built on
    the other convention misfired exactly at the threshold.

    _log10() already pins the precision at 50 digits, which is cheap and
    resolves an exponent of 30000 to far more digits than the comparison needs.
    """
    return -_log10(b) >= digits


def m_le_r_le_est(a, est, r):
    """Bracket invariant: the new left endpoint must not overshoot r, and r must
    not pass the new right endpoint.  Checked every level so that a bound which
    ever under-shot would surface as FAILURE instead of silent garbage."""
    return a <= r <= est


def _step(a, w, alpha, u_aware):
    """One level of either scheme.  Returns (est, b).

    Both schemes compute the SAME est and apply the SAME bracket update
    a' = est - b, w' = b; the only difference is which bound supplies b.  That
    is the entire experiment, so nothing else may differ.  An earlier version
    held `a` fixed on the published side while the u-aware side moved it, which
    confounds the comparison -- and, because `a` only ever increases, quietly
    favoured the baseline: smaller a gives smaller D, hence a larger b, hence
    fewer levels.
    """
    p = _p(a, alpha)
    m = a + p / (2 * a + w)
    est = corrected_from_m(a, m, p)
    b = bound_known_interval(a, w, m, est) if u_aware else bound_published(a, w)
    return est, b


def _run(a0, alpha, r, digits, u_aware, max_levels=MAX_LEVELS):
    a, w = a0, mpf(1)
    n = 0
    while n < max_levels:
        est, b = _step(a, w, alpha, u_aware)
        if b == 0:
            # p underflowed to exactly 0: `a` is the correctly-rounded mpf
            # root, so we reached working precision.  That is success, not the
            # failure sentinel -- an earlier version returned -1 here and
            # reported the winning scheme as having failed.
            return n + 1
        if _reached(b, digits):
            # Stop BEFORE the invariant check.  The invariant est-b <= r is a
            # precondition for the NEXT bracket, not a certificate of this
            # bound: once the true est-r falls below the working precision
            # (here b reached 1e-911 at dps 560), est - r is rounding noise
            # larger than b, and the check reports FAILURE for a run that has
            # already succeeded.  An earlier ordering did exactly that.
            return n + 1
        if b < 0 or not m_le_r_le_est(est - b, est, r):
            return FAILURE
        a, w = est - b, b
        n += 1
    return FAILURE


def depth_published(a, alpha, r, digits, max_levels=MAX_LEVELS):
    """Full-interval bound.  Same driver as depth_adaptive; only b differs."""
    return _run(a, alpha, r, digits, u_aware=False, max_levels=max_levels)


def depth_adaptive(a, alpha, r, digits, max_levels=MAX_LEVELS):
    """u-aware bound.  Same driver as depth_published; only b differs.

    `r` (the true sqrt) supplies the initial bracket on both sides and never
    enters a bound calculation.
    """
    return _run(a, alpha, r, digits, u_aware=True, max_levels=max_levels)


def _log10(x):
    """log10 at 50 digits regardless of ambient working precision.

    mp.log10() runs at mp.dps, which here can exceed 30000; that made it the
    dominant cost of the whole script.  A digit count needs far less.
    """
    if x == 0:
        return mp.ninf
    with mp.workdps(50):
        return mp.log10(x)


def trace(a, alpha, r, digits, u_aware, max_levels=MAX_LEVELS):
    """Digits certified after each level, for either scheme.  Same driver as
    the depth counters -- depth and trace are the same walk, reported
    differently."""
    w = mpf(1)
    out = []
    for _ in range(max_levels):
        est, b = _step(a, w, alpha, u_aware)
        if b == 0:
            out.append(mp.inf)
            break
        out.append(-_log10(b))
        if out[-1] >= digits:
            break
        a, w = est - b, b
    return out


def trace_full(a, alpha, r, digits, max_levels=MAX_LEVELS):
    return trace(a, alpha, r, digits, u_aware=False, max_levels=max_levels)


def trace_adaptive(a, alpha, r, digits, max_levels=MAX_LEVELS):
    return trace(a, alpha, r, digits, u_aware=True, max_levels=max_levels)


def cost_fixed_precision(depth, digits):
    """Model A: working precision is set once to digits + guard, as any mpmath
    caller does, so every level costs the same.  Work is proportional to depth
    and the saving is exactly the depth ratio."""
    return float(depth)


def cost_adaptive_precision(trace, digits):
    """Model B: each level is run at only the precision it needs, and the final
    level is clamped to what the caller asked for.  Work is the sum of the
    per-level precisions.  This is the model in which 'precision triples per
    level' is even relevant -- and it bounds the saving by the geometric sum,
    not by 3^depth: both schemes must pay for the final level in full, so
    saving two levels cannot buy 9x."""
    guard = 60
    total = 0.0
    for d in trace:
        if d == mp.inf or d >= digits:
            total += digits + guard
            break
        total += float(d) + guard
    return total


def run(digits, dps):
    """Both schemes at the SAME working precision, returning digit-traces.
    workdps keeps this from leaking precision into the caller (an earlier
    version assigned mp.dps globally, which made results depend on call
    order)."""
    with mp.workdps(dps):
        out = []
        for al in SAMPLE:
            alpha = mpf(al)
            a = mpf(math.isqrt(al))
            r = mp.sqrt(alpha)
            out.append((trace_full(a, alpha, r, digits),
                        trace_adaptive(a, alpha, r, digits)))
        return out


def _depth_of(trace):
    """Level count at which the target was first met."""
    return len(trace) if trace else FAILURE


def _assert_ok(traces, digits):
    """A short trace means the scheme never reached `digits`.  Letting that
    through would put FAILURE into st.median and silently skew it."""
    for label, tr in traces:
        assert tr, f"{label} produced no trace at {digits} digits"
        assert tr[-1] == mp.inf or tr[-1] >= digits, (
            f"{label} stopped short at {digits} digits: {tr}")


def main():
    print("Recursion depth: full-interval bound vs u-aware bound")
    print("median over 16 values of alpha; working precision = digits + 60\n")
    print(f"{'digits':>7} | {'full':>4} | {'u-aw':>4} | {'saved':>5} | "
          f"{'fixed-prec':>10} | {'adaptive-prec':>13} | stability")
    print("-" * 78)

    for digits in (20, 100, 1000, 3000, 10000, 30000):
        dps = digits + 60
        pairs = run(digits, dps)
        for pf, pa in pairs:
            _assert_ok([("full", pf), ("u-aware", pa)], digits)

        # Necessary-but-not-sufficient check: doubling the working precision
        # leaves every trace unchanged, so no run was sitting within the guard
        # of a resolution boundary at the precision used.  This rules out the
        # coarse failure modes; it is NOT a proof that none was.
        pairs2 = run(digits, dps * 2)
        stable = ([_depth_of(pf) for pf, _ in pairs]
                  == [_depth_of(pf) for pf, _ in pairs2]
                  and [_depth_of(pa) for _, pa in pairs]
                  == [_depth_of(pa) for _, pa in pairs2])

        med_pf = st.median([_depth_of(pf) for pf, _ in pairs])
        med_pa = st.median([_depth_of(pa) for _, pa in pairs])
        # Model A: fixed working precision -> work proportional to depth
        ca_f = st.median([cost_fixed_precision(_depth_of(pf), digits)
                          for pf, _ in pairs])
        ca_a = st.median([cost_fixed_precision(_depth_of(pa), digits)
                          for _, pa in pairs])
        # Model B: per-level precision -> work is the sum of per-level precisions
        cb_f = st.median([cost_adaptive_precision(pf, digits)
                          for pf, _ in pairs])
        cb_a = st.median([cost_adaptive_precision(pa, digits)
                          for _, pa in pairs])
        saved = med_pf - med_pa
        print(f"{digits:>7} | {med_pf:>4.0f} | {med_pa:>4.0f} | {saved:>5.0f} | "
              f"{ca_f / ca_a:>9.2f}x | {cb_f / cb_a:>12.2f}x | "
              f"{'consistent' if stable else 'UNSTABLE'}")

    print("\nPer-alpha depth at 3000 digits (a is NOT ignored: both sides vary):")
    with mp.workdps(3060):
        for al in SAMPLE:
            alpha = mpf(al)
            a = mpf(math.isqrt(al))
            r = mp.sqrt(alpha)
            f = trace_full(a, alpha, r, 3000)
            g = trace_adaptive(a, alpha, r, 3000)
            print(f"  alpha={al:>7}  a={int(a):>4}  full={len(f)}  u-aware={len(g)}")

    print("\nDigits certified per level (alpha=17, target 3500):")
    with mp.workdps(3600):
        alpha = mpf(17)
        a = mpf(4)
        r = mp.sqrt(alpha)
        f = trace_full(a, alpha, r, 3500)
        g = trace_adaptive(a, alpha, r, 3500)
        print("  full :", " ".join(f"{float(x):8.1f}" for x in f))
        print("  u-aw :", " ".join(f"{float(x):8.1f}" for x in g))
        def ratios(t):
            return [float(t[i + 1] / t[i]) for i in range(len(t) - 1)
                    if t[i] != mp.inf and t[i + 1] != mp.inf and t[i] > 0]
        gf, ga = ratios(f), ratios(g)
        print(f"  growth factor/level: full {gf[0]:.2f}..{gf[-1]:.2f}"
              f"   u-aware {ga[0]:.2f}..{ga[-1]:.2f}")


if __name__ == "__main__":
    main()
