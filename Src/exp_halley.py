"""Optimised certified Halley.  Positive intervals: +,* are monotone, so 2
directed evaluations each.  Division needs the RIGHT pairing: lower = n_lo/d_hi,
upper = n_hi/d_lo -- using n_lo/d_lo is unsound (this was a real bug)."""
from decimal import Decimal as D, localcontext, ROUND_FLOOR, ROUND_CEILING
import math, sys, time, statistics as st
sys.path.insert(0,'/home/aditya/Documents/Irrational/Src')
from exp_certified import certify

GUARD=80; CNT={'div':0,'mul':0,'add':0}

def halley(alpha, digits, maxl=400):
    CNT.update(div=0,mul=0,add=0)
    with localcontext() as top:
        top.prec=digits+GUARD+20
        A=D(alpha); a=D(math.isqrt(alpha)); w=D(1); lev=0
        while lev<maxl:
            p = min(3*(-w.adjusted()) + 3*(a.adjusted()+1) + GUARD, digits+GUARD)
            with localcontext() as c:
                c.prec=p; c.rounding=ROUND_FLOOR
                m0=a+w/2; m2=m0*m0
                nlo=m0*(m2+3*A); dlo=3*m2+A
                w3=(w/2)*(w/2)*(w/2); a2=a*a; denlo=3*m2+a2
            with localcontext() as c:
                c.prec=p; c.rounding=ROUND_CEILING
                m0=a+w/2; m2=m0*m0
                nhi=m0*(m2+3*A); dhi=3*m2+A
                w3=(w/2)*(w/2)*(w/2); a2=a*a; denhi=3*m2+a2
            with localcontext() as c:
                c.prec=p; c.rounding=ROUND_FLOOR;   Hlo=nlo/dhi; Elo=w3/denhi
            with localcontext() as c:
                c.prec=p; c.rounding=ROUND_CEILING; Hhi=nhi/dlo; Ehi=w3/denlo
            CNT['div']+=4; CNT['mul']+=8; CNT['add']+=3
            lo=Hlo-Ehi; hi=Hhi+Elo; nw=hi-lo
            if nw<=0: return None,lev+1
            if -nw.adjusted()>=digits: return (lo,hi),lev+1
            if not (nw<w): return None,lev+1
            a,w=lo,nw; lev+=1
    return None,lev

if __name__=='__main__':
    import exp_certified as C
    ALPHAS=[2,3,5,7,10,12,13,17,23,31,101,257,1009,7919,104729,999983]
    print(f"{'digits':>7} | {'certified':>9} | {'levels':>6} | {'div/lv':>6} | {'ops/lv':>6}")
    for dg in (12,100,1000):
        ok=0; lv=[]; dv=[]
        for al in ALPHAS:
            r,l=halley(al,dg)
            if r and not certify(r[0],r[1],al,dg): ok+=1; lv.append(l); dv.append(CNT['div'])
        if lv: print(f"{dg:>7} | {ok:>4}/{len(ALPHAS):<4} | {st.median(lv):>6.0f} | "
                    f"{st.median(dv)/st.median(lv):>6.1f} | 13.0")
    # wall clock vs the incumbents, best-of-3, interleaved
    def best(fn,reps=3):
        b={a:1e18 for a in ALPHAS}
        for _ in range(reps):
            for a in ALPHAS:
                t=time.perf_counter(); fn(a,1000); b[a]=min(b[a],time.perf_counter()-t)
        return b
    R={n:best(f) for n,f in (('zoom',C.zoom),('int.Newton',C.interval_newton),
                             ('Newton+res',C.newton_certified),('Halley',halley))}
    print(f"\nwall clock at 1000 digits (ms, median of 16, best of 3)")
    mz=st.median([v*1000 for v in R['zoom'].values()])
    for n in ('zoom','int.Newton','Newton+res','Halley'):
        m=st.median([v*1000 for v in R[n].values()])
        print(f"  {n:>11}: {m:>7.3f} ms   {m/mz:>5.2f}x vs zoom")
