"""M2: enemy generator for the uniform radical-sum separation.

For each rank r in 3..6 we take the prime radicands p_1,...,p_r and, for
increasing box bound B, produce a near-extremal coefficient vector by the
pigeonhole/gap construction (sort all v(C)=sum C_i sqrt(p_i), C in {0..B}^r,
and take the smallest consecutive gap).  This yields one specific embedding
sigma_+ with |L_{sigma_+}| small.  We then inspect the FULL profile of the
2^r sign-embeddings of the SAME coefficient vector:

    E = log2(1/|S|),  N = input_bits,  E/N,
    sigma_* = argmin |L_sigma|,
    the normalized profile {log2|L_sigma| / N},
    P_* = prod_{sigma != sigma_*} |L_sigma|,
    log2(P_*)/N.

The norm identity says |L_{sigma_*}| * P_* = |N_{K/Q}(S)| >= 1 after
clearing denominators, so |L_{sigma_*}| >= 1/P_*.  The experiment tells us
whether, at small embedding, P_* is much smaller than the trivial
(m M^{3/2})^{2^r-1}.
"""
import itertools
from mpmath import mp, mpf, sqrt, log, fabs
LN2 = log(2)
def log2(x): return log(x)/LN2

mp.dps = 80

def input_bits_vec(Cs, ds):
    t = 0
    for c, d in zip(Cs, ds):
        if c == 0:
            continue
        t += abs(c).bit_length() + d.bit_length()
    return t

def all_values(B, ds):
    out = []
    for C in itertools.product(range(B + 1), repeat=len(ds)):
        v = sum(C[i] * sqrt(ds[i]) for i in range(len(ds)))
        if v == 0:
            continue
        out.append((v, C))
    return out

def profile_and_gap(ds, B, want):
    vals = all_values(B, ds)
    vals.sort(key=lambda t: t[0])
    # smallest consecutive gap as the enemy S
    best = None
    for a in range(len(vals) - 1):
        va, Ca = vals[a]
        vb, Cb = vals[a + 1]
        gap = vb - va
        if gap <= 0:
            continue
        C = tuple(Cb[i] - Ca[i] for i in range(len(ds)))
        if all(c == 0 for c in C):
            continue
        if best is None or gap < best[0]:
            best = (gap, C)
    if best is None:
        return None
    gap, C = best
    S = sum(C[i] * sqrt(ds[i]) for i in range(len(ds)))
    # full sign-embedding profile of the SAME C
    embs = []
    for eps in itertools.product([1, -1], repeat=len(ds)):
        v = sum(C[i] * eps[i] * sqrt(ds[i]) for i in range(len(ds)))
        embs.append(v)
    i_star = min(range(len(embs)), key=lambda i: fabs(embs[i]))
    P_star = mpf(1)
    for i in range(len(embs)):
        if i == i_star:
            continue
        P_star *= fabs(embs[i]) if embs[i] != 0 else mpf(0)
    N = input_bits_vec(C, ds)
    if N == 0 or S == 0:
        return None
    E = -log2(fabs(S))
    prof = [log2(fabs(v)) / N if v != 0 else mpf('-inf') for v in embs]
    return {
        'C': C, 'N': N, 'E': float(E), 'EN': float(E) / N,
        'sigma_star': i_star, 'tiny_S': float(S),
        'P_star_over_N': float(log2(P_star) / N),
        'emb_abs_sorted': sorted(float(fabs(v)) for v in embs),
        'log_profile_norm': sorted(float(p) for p in prof if p != mpf('-inf')),
    }

def classify(res):
    # crude taxonomy by the magnitudes of the non-minimizing embeddings
    logs = [float(log2(x)) for x in res['emb_abs_sorted'] if x > 0]
    logs.sort()
    k = len(logs)
    if k >= 2:
        spread = logs[-1] - logs[k // 2]
    else:
        spread = 0.0
    if k >= 4 and logs[k - 2] - logs[k // 2] < spread * 0.2:
        return 'clustered-large'
    if spread > 3:
        return 'one-tiny rest-large'
    return 'comparable'

PRIMES = [2, 3, 5, 7, 11, 13, 17, 19]
for r in (3, 4, 5, 6):
    ds = PRIMES[:r]
    print(f'=== r={r}, d={ds} ===')
    for B in (3, 4, 5, 6, 8):
        res = profile_and_gap(ds, B, 1)
        if res is None:
            continue
        res['tax'] = classify(res)
        print(f"  B={B}: N={res['N']:>3} E/N={res['EN']:.3f} "
              f"log2P*/N={res['P_star_over_N']:.3f} tax={res['tax']} C={res['C']}")
    print()
