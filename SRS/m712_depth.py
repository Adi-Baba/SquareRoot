"""
The softness-depth function (Docs/41).  MEASUREMENT ONLY -- decides nothing.

Race identity (Docs/41 41.1, proved): for large k,
    eta_k = || beta_j eps^k ||  =  | beta_j eps^{-k} - n_j/(p*+q sqrt(a_j)) | ,
    n_j = p*^2 - a_j q_k^2 ,  p* = nearest integer to q_k sqrt(a_j),
with the parity sign caveat for N(eps) = -1 (a_1 = 2): the beta-term enters
with sign (-1)^k; magnitudes are unaffected.

Define the RELATIVE DEFECT and SOFTNESS DEPTH
    nu_j(k) := |n_j(k)| / q_k ,
    L_j(k)  := log_eps( q_k / |n_j(k)| ) = -log_eps nu_j(k).
Then eta_k ~ nu_j(k)/(2 sqrt(a_j)) and the Class-B wall (Docs/41 41.2) is:
    exceptions <=> |n_j| <~ eps^{(1-delta)k}  <=>  L_j(k) >= delta*k - O(1),
delta = 1/(3n(n-1)).  Liouville cap: L <= log_eps q_k ~ k (attained only on
Class A, n_j = +-1).  Observed (this module): L is FLAT ~ O(10) over k <= 1e5.

Outputs: race-identity check, Class-A census (|n_j| = 1), depth records,
depth tail counts vs the metric prediction ~ 2 K eps^{-L}, Class-B census.
"""

from math import isqrt, log

A1, AJ = 2, 3
N = 2
DELTA = 1.0 / (3 * N * (N - 1))
LOG_EPS = log(1.0 + 2.0 ** 0.5)
LN2 = log(2.0)


def _log_int(x: int) -> float:
    bl = x.bit_length()
    if bl <= 64:
        return log(x)
    return (bl - 64) * LN2 + log(x >> (bl - 64))


def depth_row(q: int, aj: int, le: float):
    """(n_j, L) from the exact integer kernel; q = q_k, L = log_eps(q/|n_j|)."""
    v = aj * q * q
    p = isqrt(v)
    d_lo = v - p * p
    d_hi = (p + 1) ** 2 - v
    n = d_lo if d_lo <= d_hi else d_hi          # |n_j|, exact
    L = (_log_int(q) - _log_int(n)) / le if n > 0 else float("inf")
    return n, L


def fpu(a1: int):
    """Fundamental unit (p,q) of x^2 - a1 y^2 = +-1 via the CF of sqrt(a1).

    Returns the first convergent with |p^2 - a1 q^2| = 1 (norm -1 or +1 --
    either is a valid Pell unit for the orbit; Docs/27's convention).
    """
    a0 = isqrt(a1)
    if a0 * a0 + 1 == a1:          # 0-th convergent a0/1 has norm -1
        return a0, 1
    m, d, a = 0, 1, a0
    p2, p1, q2, q1 = 1, a0, 0, 1        # p_{-1}, p_0; q_{-1}, q_0
    while True:
        m = d * a - m
        d = (a1 - m * m) // d
        a = (a0 + m) // d
        p2, p1 = p1, a * p1 + p2
        q2, q1 = q1, a * q1 + q2
        if abs(p1 * p1 - a1 * q1 * q1) == 1:
            return p1, q1


def scan(a1: int, aj: int, kmax: int):
    """Incremental exact walk of eps^k for the TRUE Pell unit eps = p + q sqrt(a1)."""
    from math import sqrt
    pe, qe = fpu(a1)
    le = log(pe + qe * sqrt(a1))               # per-pair log eps (fixes unit caveat)
    P, Q = 1, 0
    for k in range(1, kmax + 1):
        P, Q = P * pe + a1 * qe * Q, P * qe + pe * Q
        yield k, Q, depth_row(Q, aj, le)


def census(a1: int, aj: int, kmax: int = 50000):
    """One-pair summary: Class-A events, max depth, frontier, tail profile."""
    cls_a = []
    rec_L, rec_k = -1.0, 0
    front = (1e18, 0, 0)          # min log|n|/k
    tail = {}
    for k, q, (n, L) in scan(a1, aj, kmax):
        if n == 1 and k > 1:
            cls_a.append(k)
        if L > rec_L:
            rec_L, rec_k = L, k
        e = _log_int(n) / k if k > 1 else 1e18
        if e < front[0]:
            front = (e, k, n)
        for ell in (4, 8, 12):
            if L >= ell:
                tail[ell] = tail.get(ell, 0) + 1
    return {"class_a": cls_a, "maxL": rec_L, "argmax": rec_k,
            "front": front, "tail": tail}


def main(kmax: int = 100000) -> None:
    from mpmath import mp
    mp.dps = 60
    e = mp.mpf(1) + mp.sqrt(A1)
    beta = mp.sqrt(mp.mpf(AJ) / A1) / 2
    print("SOFTNESS DEPTH  (a1,a_j)=(%d,%d)  n=%d  delta=1/%d  k=1..%d"
          % (A1, AJ, N, round(1 / DELTA), kmax))

    # --- race identity check at both parities (exact vs mpmath) ---
    print("\nrace identity eta_k = |(+/-)beta eps^-k - n_j/(p*+q sqrt aj)|:")
    for k in (20, 21, 25, 26, 40, 41):
        q = list(scan(A1, AJ, k))[-1][1]
        v = AJ * q * q
        p = isqrt(v)
        n = (v - p * p) if (v - p * p) <= ((p + 1) ** 2 - v) else ((p + 1) ** 2 - v)
        pstar = p if (v - p * p) <= ((p + 1) ** 2 - v) else p + 1
        x = beta * e ** k
        eta = abs(x - mp.floor(x + mp.mpf("0.5")))
        lhs = float(eta)
        sgn = 1 if k % 2 == 0 else -1     # beta-term sign by parity (N(eps) = -1)
        rhs = abs(sgn * float(beta * e ** (-k)) - n / (pstar + q * float(mp.sqrt(AJ))))
        print("  k=%-3d parity %s  eta=%.12f  race=%.12f  diff=%.2e"
              % (k, "even" if k % 2 == 0 else " odd", lhs, rhs, abs(lhs - rhs)))

    # --- full scan ---
    class_a = []
    class_b = []
    rec_L, rec_k = -1.0, 0
    tail = {}
    for k, q, (n, L) in scan(A1, AJ, kmax):
        if n == 1:
            class_a.append(k)
        if n >= 2 and L >= DELTA * k - 12.0:   # generous O(1) slack
            class_b.append((k, n, L))
        if L > rec_L:
            rec_L, rec_k = L, k
        for ell in (2, 4, 6, 8, 10, 12):
            if L >= ell:
                tail[ell] = tail.get(ell, 0) + 1

    print("\nClass-A census (|n_j| = 1, simultaneous Pell): k in %s"
          % (class_a if class_a else "NONE"))
    print("Class-B census (2 <= |n_j| and L >= delta*k - 12): %s"
          % (class_b[:10] if class_b else "NONE"))
    print("max depth L = %.3f at k=%d   (delta*k there = %.1f)"
          % (rec_L, rec_k, DELTA * rec_k))
    print("\ndepth tail counts vs metric prediction ~ 2*K*eps^-L:")
    print("  L   count   2*K*eps^-L")
    for ell in sorted(tail):
        pred = 2 * kmax * 2.718281828 ** (-ell * LOG_EPS)
        print("  %-3d %-7d %.3g" % (ell, tail[ell], pred))


if __name__ == "__main__":
    main()