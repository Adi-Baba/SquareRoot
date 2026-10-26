import itertools, math
import numpy as np

def factor_set(n):
    s=set(); p=2
    while p*p<=n:
        while n%p==0: n//=p; s.add(p)
        p+=1
    if n>1: s.add(n)
    return s

def ses(ds):
    prlist=sorted({p for d in ds for p in factor_set(d)})
    M=[]
    for d in ds:
        v=[0]*len(prlist)
        for p in factor_set(d):
            k=prlist.index(p); v[k]=1
        M.append(v)
    M=np.array(M, dtype=int)%2
    M=M.tolist()
    r=0
    for c in range(len(prlist)):
        sel=None
        for i in range(r,len(M)):
            if M[i][c]==1: sel=i; break
        if sel is not None:
            M[r],M[sel]=M[sel],M[r]
            for i in range(len(M)):
                if i!=r and M[i][c]==1:
                    M[i]=[(M[i][k]+M[r][k])%2 for k in range(len(prlist))]
            r+=1
    return r

def smallness(B, ds):
    small_pairs=0
    classinfo={"full":0,"partial":0}
    best=1e18
    for C in itertools.product(range(-B,B+1), repeat=len(ds)):
        if all(c==0 for c in C): continue
        xs=[c*math.sqrt(d) for c,d in zip(C,ds)]
        vals=[]
        for eps in itertools.product([1,-1], repeat=len(ds)):
            v=sum(x*e for x,e in zip(xs,eps))
            vals.append((abs(v),eps))
        vals.sort()
        if vals[0][0]<1e-1 and vals[1][0]<1e-1:
            small_pairs+=1
            diffs=[d for i,d in enumerate(ds) if vals[0][1][i]!=vals[1][1][i]]
            if diffs:
                rr=ses(diffs)
                if rr==len(diffs):
                    classinfo["full"]+=1
                else:
                    classinfo["partial"]+=1
            b=vals[0][0]*vals[1][0]
            if b<best: best=b
    return small_pairs,classinfo,best

for ds,B in [((2,3,5),12),((2,3,5,7),7)]:
    r=ses(ds)
    sp,ci,best=smallness(B,ds)
    print(ds,"r=",r,"B=",B,"pairs_both<1e-1:",sp,ci,"best_prod:",best)
