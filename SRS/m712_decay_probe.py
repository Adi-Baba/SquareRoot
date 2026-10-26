"""
T3 decay probe (TODO L13).  MEASUREMENT ONLY -- decides nothing.

Measures the proxy

    eta_proxy(k) = || q_k sqrt(a_j) ||

for the calibration pair (a_1, a_j) = (2, 3), k = 1..KMAX, from exact integers:

    q_k = engine Pell denominator  (m712_engine.eps_power; k IS the exponent),
    p   = isqrt(a_j q^2)  (floor of q sqrt(a_j), exact),
    n   = min(a_j q^2 - p^2, (p+1)^2 - a_j q^2)   (exact distance-squared numerator),
    eta = n / (p* + q sqrt(a_j)),  p* the nearer integer,
          evaluated as an exact-int numerator over a positive mpmath denominator
          (no cancellation; 60 significant digits throughout).

Proxy validity: ||q_k sqrt(a_j)|| = ||beta_j eps^k|| + O(eps^{-k}) (Docs/34
corrected identity, verified).  eps^{-k} < 1e-100 for k >= 250 while eta is
Theta(1) empirically, so the correction is invisible at every scale probed.

Questions (TODO L13 T3):
  (i)  does eta ever dip below eta_min = 0.0036793 (attained at k=43, Docs/36 36.6)?
  (ii) does the running minimum trend downward?
  (iii) does G1(k) = -log(eta)/(k log eps) trend to 0 (no exponential decay)
        or to a positive constant (Cor 30.4's e^{-o(k)} is empirically loose)?
  (iv) observed power-law constant C_obs = max_k (-log eta_k)/log k.
All outputs are measurements; no claim is made or refuted.
"""

from math import isqrt, log

import m712_engine as E

A1, AJ = 2, 3
LOG_EPS = log(1.0 + 2.0 ** 0.5)


def eta_proxy(k: int, mp):
    pt = E.eps_power(A1, k)
    q = pt.qk
    v = AJ * q * q
    p = isqrt(v)
    d_lo = v - p * p
    d_hi = (p + 1) ** 2 - v
    if d_lo <= d_hi:
        n, pstar = d_lo, p
    else:
        n, pstar = d_hi, p + 1
    denom = mp.mpf(pstar) + mp.mpf(q) * mp.sqrt(AJ)
    return mp.mpf(n) / denom


LN2 = log(2.0)


def _log_int(x: int) -> float:
    """log(x) for a positive Python int, error < 1e-18 (64-bit mantissa trick)."""
    bl = x.bit_length()
    if bl <= 64:
        return log(x)
    m = x >> (bl - 64)
    return (bl - 64) * LN2 + log(m)


def eta_log(k: int) -> float:
    """log ||q_k sqrt(AJ)|| with the SAME exact integer kernel, log-only.

    Denominator: p* + q sqrt(AJ) = 2 q sqrt(AJ) (1 + O(1/q)) because
    |q sqrt(AJ) - p*| <= 1/2 for the nearest integer p*, so for q >= 10 the
    relative error of the approximation is < 1e-10 and its log error is far
    below the 1e-18 mantissa error.  No mpmath in the hot loop.
    """
    pt = E.eps_power(A1, k)
    q = pt.qk
    return _eta_log_from_q(q)


def _eta_log_from_q(q: int) -> float:
    v = AJ * q * q
    p = isqrt(v)
    d_lo = v - p * p
    d_hi = (p + 1) ** 2 - v
    n = d_lo if d_lo <= d_hi else d_hi
    return _log_int(n) - _log_int(q) - log(2.0 * (AJ ** 0.5))


def scan(kmax: int):
    """Incremental exact scan: one Pell recurrence walk, k IS the exponent.

    Same recurrence m712_engine.eps_power certifies: (P,Q) -> (P + a1 Q, P + Q)
    from (1,0), so Q after k steps IS q_k -- no index drift, no rebuild.
    Yields (k, log eta_k) for k = 1..kmax.
    """
    p_, q_ = 1, 0
    for k in range(1, kmax + 1):
        p_, q_ = p_ + A1 * q_, p_ + q_
        yield k, _eta_log_from_q(q_)


def main(kmax: int = 5000) -> None:
    from mpmath import mp
    mp.dps = 60
    print("T3 DECAY PROBE  (a1,a_j)=(%d,%d)  k=1..%d  MEASUREMENT ONLY" % (A1, AJ, kmax))
    print("threshold theta_0^k = eps^{-k/6}, log eps = %.6f" % LOG_EPS)

    run_min = None
    run_min_k = 0
    c_obs = 0.0
    g1_first = None
    g1_last = None
    exceptions = []
    records = []
    checkpoints = [10, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000]
    if kmax > 20000:
        checkpoints += [50000, 100000]
    # fast path validation: exact vs log formula at three k
    if kmax > 1000:
        for kv in (1000, 1586, 5000):
            exact = log(float(eta_proxy(kv, mp)))
            fast = eta_log(kv)
            print("validate k=%-5d exact %.12f  fast %.12f  diff %.2e"
                  % (kv, exact, fast, abs(exact - fast)))
    print("\n  k     eta_k          running-min     argmin   G1(k)     -log(eta)/log k")
    for k, le in scan(kmax):
        g1 = -le / (k * LOG_EPS)
        if k in (1, 2):
            g1_first = g1
        g1_last = g1
        pow_c = -le / log(k) if k >= 2 else 0.0
        if pow_c > c_obs:
            c_obs = pow_c
        if le < -k * LOG_EPS / 6.0:
            exceptions.append(k)
        if run_min is None or le < run_min:
            run_min, run_min_k = le, k
            records.append((k, le))
        if k in checkpoints:
            print("%4d  e^%-11.4f e^%-12.4f %-8d %-9.5f %.4f"
                  % (k, le, run_min, run_min_k, g1, pow_c))

    print("\nRECORD EVENTS (k, -log(eta)/log k) -- law test: slope of log(eta) vs log(k)")
    for k, le in records:
        print("   k=%-6d eta=e^%.3f   log-log exponent %.3f"
              % (k, le, -le / log(k) if k >= 2 else 0.0))
    if len(records) >= 3:
        ks = [log(k) for k, _ in records[1:]]
        ys = [le for _, le in records[1:]]
        mk = sum(ks) / len(ks)
        my = sum(ys) / len(ys)
        slope = (sum((a - mk) * (b - my) for a, b in zip(ks, ys))
                 / sum((a - mk) ** 2 for a in ks))
        print("   log-log slope over record events: %.4f  (=> running min ~ k^%.4f)"
              % (slope, -slope))
    print("\n(i)  final minimum e^%.4f = %.6g at k=%d" % (run_min, _e(run_min), run_min_k))
    print("(iii) G1: k=1 %.4f -> k=%d %.6f  (0 => no exponential decay trend)"
          % (g1_first, kmax, g1_last))
    print("(iv)  C_obs = max_k -log(eta_k)/log k = %.4f  (target needs SOME finite C; "
          "observed bound eta >= k^{-C_obs} over the range)" % c_obs)
    print("dips below the BARE scale eps^{-k/6} (constant factor of the exact EHS"
          " threshold not applied; certified exception list is Docs/27): %s"
          % (exceptions if exceptions else "none in range"))


def _e(x: float) -> float:
    from math import exp
    return exp(x)


if __name__ == "__main__":
    main()