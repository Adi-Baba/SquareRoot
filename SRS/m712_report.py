"""
Reporting layer for m712_engine.  MEASUREMENT ONLY.

This module exists so that m712_engine.py can stay exactly decidable.  Nothing
here decides anything: it takes logarithms of quantities the engine computed
exactly, and prints.  tests/test_no_float.py lists it as a MEASUREMENT module
for that reason, and m712_engine.py as a DECISION_PATH module.

The question this is built to answer is the research one from Docs/35:
which arithmetic quantity, if any, grows as d_k = ||x_k|| falls?  Answer, for
the record: NONE of them do, and the conjugate-amplification candidate is
provably flat at 2 beta_j^3 (2 a_1)^4 (Docs/36 section 36.6).
"""

import m712_engine as E

__all__ = ["_mp", "profile", "main"]


def _mp(a1, aj):
    from mpmath import mp
    mp.dps = 50
    e = mp.mpf(1) + mp.sqrt(a1)
    ep = mp.mpf(1) - mp.sqrt(a1)
    return mp, e, ep, mp.sqrt(mp.mpf(aj) / a1) / 2


def profile(a1: int, aj: int, kmax: int = 60, n: int = 2):
    """
    The diagnostic table: every layer's quantity, measured, one dict per k.
    This is the function the research question is about.
    """
    mp, e, ep, beta = _mp(a1, aj)
    th0 = e ** (-mp.mpf(1) / (3 * n * (n - 1)))
    le = mp.log(e)
    out = []
    for k in range(1, kmax + 1):
        pt = E.eps_power(a1, k)
        m = E.nearest_integer(a1, aj, k)
        x = (mp.mpf(pt.qk) / 2) * mp.sqrt(aj) + (mp.mpf(pt.pk) / (2 * a1)) * mp.sqrt(a1 * aj)
        dk = abs(x - mp.mpf(m))
        y = E.y_of(a1, aj, k, m)
        nrm = E.field_norm(y, a1, aj)
        nrm = mp.mpf(nrm.numerator) / mp.mpf(nrm.denominator)
        u, v, nc, _ = E.integrality_row(a1, aj, k)
        barrier = e ** (-3 * k) / (2 * beta ** 3)
        out.append({
            "k": k, "m": m,
            "d": dk,
            "G1": -mp.log(dk) / (k * le),
            "G_norm": 3 * k * le + mp.log(dk),
            "R_norm": dk / barrier,
            "A_conj": abs(nrm) / (dk * e ** (3 * k)),
            "N_y": nrm,
            "log_abs_Nc": mp.log(abs(nc)) if nc != 0 else mp.inf,
            "abs_c": abs(mp.mpf(u) + mp.mpf(v) * mp.sqrt(a1)),
            "threshold": th0 ** k,
            "exception": dk < th0 ** k,
        })
    return out


def main():
    mp, e, ep, beta = _mp(2, 3)
    print("=" * 78)
    print("M7.12 GAP PROFILER   a_1=2, a_j=3, n=2, theta_0 = eps^(-1/6) = %s"
          % mp.nstr(e ** (-mp.mpf(1) / 6), 8))
    print("=" * 78)
    print("layer 3 DEGENERATE: I_k == 0 for k in 1..39 ->",
          all(E.invariant_row(2, 3, k) for k in range(1, 40)))
    print("  Q(u,v) == (v-eps u)(v-eps' u) exactly; x_{k+1}=eps x_k kills factor 1.")
    print()
    rows = profile(2, 3, 60)
    print("{:>3} {:>11} {:>8} {:>10} {:>11} {:>11} {:>12} {:>4}".format(
        "k", "d_k", "G1", "G_norm", "R_norm", "A_conj", "|N(c)|^(1/k)", "exc"))
    for r in rows[:5] + [None] + rows[-3:]:
        if r is None:
            print("   ...")
            continue
        print("{:>3} {:>11} {:>8} {:>10} {:>11} {:>11} {:>12} {:>4}".format(
            r["k"], mp.nstr(r["d"], 6), mp.nstr(r["G1"], 4), mp.nstr(r["G_norm"], 5),
            mp.nstr(r["R_norm"], 4), mp.nstr(r["A_conj"], 4),
            mp.nstr(mp.exp(r["log_abs_Nc"] / r["k"]), 6),
            "YES" if r["exception"] else "-"))
    print()
    print("exceptions below theta_0^k, k=1..60 : %d"
          % len([r for r in rows if r["exception"]]))
    print("min d_k                            : %s at k=%d"
          % (mp.nstr(min(r["d"] for r in rows), 6),
             min(rows, key=lambda r: r["d"])["k"]))
    print()
    print("WHICH QUANTITY GROWS AS d_k FALLS?  corr with -log d_k")
    xs = [-mp.log(r["d"]) for r in rows]
    n_ = len(xs)
    mx = sum(xs) / n_

    def corr(ys):
        my = sum(ys) / n_
        cov = sum((xs[i] - mx) * (ys[i] - my) for i in range(n_))
        vx = sum((v - mx) ** 2 for v in xs) ** 0.5
        vy = sum((v - my) ** 2 for v in ys) ** 0.5
        return cov / (vx * vy) if vx and vy else mp.mpf(0)

    for name, ys in [
        ("G_norm", [r["G_norm"] for r in rows]),
        ("A_conj", [mp.log(r["A_conj"]) for r in rows]),
        ("log|N(c)|/k", [r["log_abs_Nc"] / r["k"] for r in rows]),
    ]:
        print("   corr( -log d_k , %-16s ) = %s" % (name, mp.nstr(corr(ys), 5)))
    print()
    print("=== CONJUGATE AMPLIFICATION IS EXACTLY FLAT ===")
    asym = 2 * beta ** 3 * (2 * 2) ** 4
    print("   A_conj -> %s ; predicted 2 beta^3 (2 a_1)^4 = %s ; ratio %s"
          % (mp.nstr(rows[-1]["A_conj"], 8), mp.nstr(asym, 8),
             mp.nstr(rows[-1]["A_conj"] / asym, 8)))
    print("   => the conjugate estimates in Thm 35.1 are asymptotically TIGHT;")
    print("      there is no slack to convert into a better exponent.")
    print()
    print("=== HOW CLOSE DOES d_k GET TO THE THRESHOLD? ===")
    for q, k, d, t in sorted((r["d"] / r["threshold"], r["k"], r["d"], r["threshold"])
                              for r in rows)[:4]:
        print("   k=%-3d d_k=%-12s theta_0^k=%-12s ratio=%s"
              % (k, mp.nstr(d, 5), mp.nstr(t, 5), mp.nstr(q, 5)))
    print()
    print("=== EXCEPTIONS: Theorem 36.1 does NOT isolate them ===")
    exc = sorted(r["k"] for r in rows if r["exception"])
    adj = [(exc[i], exc[i + 1]) for i in range(len(exc) - 1) if exc[i + 1] == exc[i] + 1]
    print("   exception indices:", exc)
    print("   adjacent pairs:", adj if adj else "none")


if __name__ == "__main__":
    main()
