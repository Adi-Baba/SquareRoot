"""M7.1 enemy: search for two genuinely small embeddings with no small-support relation.

For a coefficient vector C in a rank-r multiquadratic field, define
v_sigma(C) = sum_i C_i * sigma(sqrt(d_i)).  We ask:
  * is there a pair sigma != tau with both |v_sigma|,|v_tau| <= tau_scaled ?
  * for such a pair, does the difference delta = v_sigma - v_tau have support
    on many independent radicals (i.e. the sign-difference set has full rank) ?
This tests whether the "thin simultaneous regime" exists at all.
"""
import itertools
from sympy import sqrt

def all_v(C, ds):
    ds = list(ds)
    out = {}
    for eps in itertools.product([1,-1], repeat=len(ds)):
        out[eps] = sum(C[i]*eps[i]*sqrt(ds[i]) for i in range(len(ds)))
    return out

def f2_rank(ds):
    # parity vectors over primes
    vecs = []
    ps = set()
    for d in ds:
        f = {}
        n=d
        p=2
        while p*p<=n:
            c=0
            while n%p==0: n//=p; c+=1
            if c%2: f[p]=1
            p+=1
        if n>1: f[n]=1
        ps|=set(f)
        v=0
        for q in sorted(ps):
            if f.get(q): v |= 1<<sorted(ps).index(q)
        vecs.append((f))
    # rank over F2 of dict vectors
    rows=[]
    for f in vecs:
        rows.append(f)
    prlist=sorted({p for f in vecs for p in f})
    cur=[ {p:(f.get(p,0)&1) for p in prlist} for f in vecs]
    rank=0
    for col in range(len(prlist)):
        # find a row with 1 in this col nonzero below
        sel=None
        for r in range(rank,len(cur)):
            if cur[r][prlist[col]]==1:
                sel=r;break
        if sel!=None:
            cur[rank],cur[sel]=cur[sel],cur[rank]
            for r in range(len(cur)):
                if r!=rank and cur[r][prlist[col]]==1:
                    cur[r] = {p:(cur[r][p]+cur[rank][p])%2 for p in prlist}
            rank+=1
    return rank

def f2_rankv(ds):
    import numpy as np
    prlist = sorted({p for d in ds for p in factor_set(d)})
    idx={p:i for i,p in enumerate(prlist)}
    M=[]
    for d in ds:
        v=[0]*len(prlist)
        for p in factor_set(d): v[idx[p]]=1
        M.append(v)
    M=np.array(M,dtype=int)%2
    # rank over F2 via mod2 elimination
    M=M.tolist()
    r=0
    for c in range(len(prlist)):
        sel=None
        for i in range(r,len(M)):
            if M[i][c]==1: sel=i;break
        if sel is not None:
            M[r],M[sel]=M[sel],M[r]
            for i in range(len(M)):
                if i!=r and M[i][c]==1:
                    M[i]=[(M[i][k]+M[r][k])%2 for k in range(len(prlist))]
            r+=1
    return r

def factor_set(n):
    s=set(); p=2
    while p*p<=n:
        while n%p==0: n//=p; s.add(p)
        p+=1
    if n>1: s.add(n)
    return s

def runs():
    for (ds,B) in [((2,3,5),10),((2,3,5,7),6),((2,3,5,7,11),4)]:
        r=f2_rankv(ds)
        two_small=0; two_small_fullrank_diff=0; total=0
        worstpair=None
        for C in itertools.product(range(-B,B+1), repeat=len(ds)):
            if all(c==0 for c in C): continue
            total+=1
            vs=list(all_v(C,ds).items())
            # find smallest two absolute values among embeddings
            pairs=[]
            for a in range(len(vs)):
                for b in range(a+1,len(vs)):
                    ea,va=vs[a]; eb,vb=vs[b]
                    pairs.append((float(abs(va)),float(abs(vb)),ea,eb,va,vb))
            if not pairs: continue
            # threshold: small means both <= B^0.6 * max|sqrt d| * something; use absolute smallness proxy:
            small=[p for p in pairs if p[0]<1e-2 and p[1]<1e-2]
            if small:
                two_small+=1
                # measure difference support rank
                a0,b0,ea,eb,va,vb=small[0]
                diffs=[]
                for i in range(len(ds)):
                    if ea[i]!=eb[i]:
                        diffs.append(ds[i])
                rr=f2_rankv(diffs) if diffs else 0
                if rr>=len(diffs) and len(diffs)>=2 and rr==r:
                    two_small_fullrank_diff+=1
                if worstpair is None or (a0+b0)<(worstpair[0]+worstpair[1]):
                    worstpair=(a0,b0,ea,eb,diffs,rr)
        print(f"ds={ds} r={r} B={B} total={total} two_small(<1e-2)={two_small} fullrank_diff_pairs={two_small_fullrank_diff}")
        if worstpair:
            print(f"   best pair: |v1|={worstpair[0]:.2e} |v2|={worstpair[1]:.2e} diffs={worstpair[4]} support_rank={worstpair[5]}")

if __name__ == "__main__":
    # This used to be a bare module-level `runs()`, so `import two_small`
    # executed a full combinatorial search (B=7 over 4 primes) and never
    # returned.  Nothing in the repo imports this module -- it is a standalone
    # probe -- so gating it preserves behaviour for every consumer and makes
    # the file safe to import.
    runs()
