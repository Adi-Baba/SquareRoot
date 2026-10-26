import itertools, math, csv, sys
import numpy as np
from collections import defaultdict

def f2_rank(ds):
    if not ds: return 0
    prlist = sorted({p for d in ds for p in factor_set(d)})
    if not prlist: return 0
    M = []
    for d in ds:
        L = [0]*len(prlist)
        for p in factor_set(d): L[prlist.index(p)] = 1
        M.append(L)
    M = [row[:] for row in M]
    r = 0
    for c in range(len(prlist)):
        sel = next((i for i in range(r, len(M)) if M[i][c] == 1), None)
        if sel is not None:
            M[r], M[sel] = M[sel], M[r]
            for i in range(len(M)):
                if i != r and M[i][c] == 1:
                    M[i] = [(M[i][k] + M[r][k]) % 2 for k in range(len(prlist))]
            r += 1
    return r

def factor_set(n):
    s = set(); p = 2
    while p*p <= n:
        while n % p == 0: n //= p; s.add(p)
        p += 1
    if n > 1: s.add(n)
    return s

def canonical_signs(r):
    reps = []
    for eps in itertools.product([1,-1], repeat=r):
        if eps[0] == 1:
            reps.append(eps)
    return reps

def run(ds, B, theta=0.25):
    r = len(ds)
    sqrts = [math.sqrt(d) for d in ds]
    signs = canonical_signs(r)
    n_signs = len(signs)
    # Precompute sign matrix as numpy array
    S = np.array(signs, dtype=np.float64)  # (n_signs, r)
    SQ = np.array(sqrts, dtype=np.float64)  # (r,)
    
    residuals = []
    total = 0
    for C in itertools.product(range(-B, B+1), repeat=r):
        if all(c == 0 for c in C): continue
        total += 1
        C_arr = np.array(C, dtype=np.float64)
        vals = S @ (C_arr * SQ)  # (n_signs,)
        abs_vals = np.abs(vals)
        order = np.argsort(abs_vals)
        i0, i1 = order[0], order[1]
        v0, v1 = abs_vals[i0], abs_vals[i1]
        if v0 == 0 or v1 == 0: continue
        E = -math.log2(v0)
        D2 = -math.log2(v1)
        if not (E > 1 and D2 >= theta * E): continue
        s0, s1 = signs[i0], signs[i1]
        J = [i for i in range(r) if s0[i] != s1[i]]
        if len(J) < r - 2: continue
        rJ = f2_rank([ds[i] for i in J])
        if rJ != len(J): continue
        # Compute X and Y
        X_coeffs = []
        Y_coeffs = []
        for i in range(r):
            if s0[i] == s1[i]:
                X_coeffs.append((i, C[i] * s0[i]))
            else:
                Y_coeffs.append((i, C[i] * (s0[i] - s1[i]) // 2))
        X_val = sum(c * sqrts[i] for i, c in X_coeffs)
        Y_val = sum(c * sqrts[i] for i, c in Y_coeffs)
        EX = -math.log2(abs(X_val)) if X_val != 0 else float('inf')
        EY = -math.log2(abs(Y_val)) if Y_val != 0 else float('inf')
        # Remaining canonical magnitudes
        remaining = sorted([abs_vals[i] for i in range(n_signs) if i != i0 and i != i1])
        # Products d_a*d_b for a,b in J
        prods = [ds[a]*ds[b] for a in J for b in J if a < b]
        residuals.append({
            'C': C, 'd': ds, 's0': s0, 's1': s1,
            'E': E, 'D2': D2, 'D2/E': D2/E,
            'J': J, 'rJ': rJ,
            'gcdC': math.gcd(*[abs(c) for c in C]),
            'suppC': sum(1 for c in C if c != 0),
            'suppX': len(X_coeffs), 'suppY': len(Y_coeffs),
            'X': X_val, 'Y': Y_val, 'EX': EX, 'EY': EY,
            'remaining': remaining, 'prods': prods,
        })
    return residuals, total

def main():
    grids = [
        ((2,3,5), [2,4,8]),
        ((2,3,5,7), [2,4,8]),
        ((2,3,5,7,11), [2,4]),
        ((2,3,5,7,11,13), [2]),
        ((2,3,5,7,11,13,17), [2]),
    ]
    all_rows = []
    summary = []
    for ds, Bs in grids:
        for B in Bs:
            print(f"Running ds={ds} B={B}...", file=sys.stderr)
            residuals, total = run(ds, B)
            r = len(ds)
            rate = len(residuals) / total if total > 0 else 0
            if residuals:
                med_D2E = sorted([x['D2/E'] for x in residuals])[len(residuals)//2]
                ex_list = sorted([x['EX']/x['E'] for x in residuals if x['EX'] != float('inf')])
                ey_list = sorted([x['EY']/x['E'] for x in residuals if x['EY'] != float('inf')])
                med_EXE = ex_list[len(ex_list)//2] if ex_list else float('nan')
                med_EYE = ey_list[len(ey_list)//2] if ey_list else float('nan')
            else:
                med_D2E = med_EXE = med_EYE = float('nan')
            summary.append({
                'r': r, 'B': B, 'samples': total, 'residuals': len(residuals),
                'rate': rate, 'med_D2/E': med_D2E, 'med_EX/E': med_EXE, 'med_EY/E': med_EYE,
            })
            for x in residuals:
                all_rows.append({
                    'r': r, 'B': B, 'd': ';'.join(map(str,ds)),
                    'C': ';'.join(map(str,x['C'])),
                    's0': ';'.join(map(str,x['s0'])),
                    's1': ';'.join(map(str,x['s1'])),
                    'E': x['E'], 'D2': x['D2'], 'D2/E': x['D2/E'],
                    'J': ';'.join(map(str,x['J'])), 'rJ': x['rJ'],
                    'gcdC': x['gcdC'], 'suppC': x['suppC'],
                    'suppX': x['suppX'], 'suppY': x['suppY'],
                    'X': x['X'], 'Y': x['Y'], 'EX': x['EX'], 'EY': x['EY'],
                    'remaining': ';'.join(f"{v:.6e}" for v in x['remaining']),
                    'prods': ';'.join(map(str,x['prods'])),
                })
            print(f"  -> {len(residuals)} residuals out of {total} (rate {rate:.4f})", file=sys.stderr)
    
    # Write CSV
    with open('m7_2_residuals.csv', 'w', newline='') as f:
        if all_rows:
            w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
            w.writeheader()
            w.writerows(all_rows)
    
    # Print summary
    print("\n=== M7.2 Summary ===")
    print(f"{'r':>3} {'B':>3} {'samples':>10} {'residuals':>10} {'rate':>8} {'med_D2/E':>10} {'med_EX/E':>10} {'med_EY/E':>10}")
    for s in summary:
        print(f"{s['r']:>3} {s['B']:>3} {s['samples']:>10} {s['residuals']:>10} {s['rate']:>8.4f} {s['med_D2/E']:>10.3f} {s['med_EX/E']:>10.3f} {s['med_EY/E']:>10.3f}")
    
    # Grouped by (|J|, rJ)
    print("\n=== By (|J|, r_J) ===")
    groups = defaultdict(list)
    for x in all_rows:
        groups[(x['r'], len(x['J'].split(';')), x['rJ'])].append(x)
    for k in sorted(groups.keys()):
        v = groups[k]
        print(f"  r={k[0]}, |J|={k[1]}, rJ={k[2]}: {len(v)} residuals")
    
    # Grouped by gcd(C)
    print("\n=== By gcd(C) ===")
    groups2 = defaultdict(list)
    for x in all_rows:
        groups2[x['gcdC']].append(x)
    for k in sorted(groups2.keys()):
        print(f"  gcd={k}: {len(groups2[k])} residuals")
    
    # Grouped by (suppX, suppY)
    print("\n=== By (suppX, suppY) ===")
    groups3 = defaultdict(list)
    for x in all_rows:
        groups3[(x['suppX'], x['suppY'])].append(x)
    for k in sorted(groups3.keys()):
        print(f"  suppX={k[0]}, suppY={k[1]}: {len(groups3[k])} residuals")

main()
