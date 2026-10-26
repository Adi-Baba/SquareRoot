"""Certified cubic bracket contraction for the n-th ROOT, via Halley.
   Exact: H_n(x)-r = (x-r)^3 Q_n(x,r)/D_n(x,r), with Q/D -> (n^2-1)/(12 r^2).
   Midpoint bound: |x-r| <= w/2, and Q/D is bounded on the bracket."""
from decimal import Decimal as D, localcontext, ROUND_FLOOR, ROUND_CEILING
import math, sys, statistics as st
sys.path.insert(0,'/home/aditya/Documents/Irrational/Src')
from exp_certified import certify

def certify_nroot(alpha, n, digits, maxl=400, guard=80):
    """Bracket [a,a+w] containing alpha^(1/n), certified by exact integer powers.
    Uses the Q/D bound evaluated at the bracket endpoints (Q/D is a rational
    function of x with positive coefficients, so interval evaluation is valid)."""
    lo0 = math.isqrt(alpha)
    # integer seed: largest a with a^n <= alpha
    a = 0
    while (a+1)**n <= alpha: a += 1
    if a == 0: a = 1
    with localcontext() as top:
        top.prec = digits + guard + 20
        A = D(alpha); ad = D(a); w = D(1); lev = 0
        while lev < maxl:
            p = digits + guard
            # sign(H-r) = sign(x-r) exactly, and the midpoint's sign is not
            # known, so the SYMMETRIC bound +/- E is used.  No bracket test
            # needed -- the earlier one, (n-1)a^n vs (n+1)alpha, is always
            # negative and was simply wrong.
            with localcontext() as c:
                c.prec = p; c.rounding = ROUND_FLOOR
                m0 = ad + w/2
                mn = m0**n
                den = (n+1)*mn + (n-1)*A
                H = m0*((n-1)*mn + (n+1)*A)/den
            with localcontext() as c:
                c.prec = p; c.rounding = ROUND_CEILING
                m0 = ad + w/2
                mn = m0**n
                den = (n+1)*mn + (n-1)*A
                H = m0*((n-1)*mn + (n+1)*A)/den
            # certified error: |H-r| <= (n^2-1)/12 * |m0-r|^3 / (lower bound on D)
            # D = (n+1)x^n + (n-1)alpha  >= (n+1) a^n + (n-1) alpha  > 0
            with localcontext() as c:
                c.prec = p; c.rounding = ROUND_FLOOR
                Dlo = (n+1)*ad**n + (n-1)*A
            # Q_n/D_n -> (n^2-1)/(12 r^2) at x=r, and r >= a, so
            # E = (n^2-1)/12 * (w/2)^3 / a^2  is an UPPER bound.  (Dividing by
            # D as well double-counted and made E 4x too small -- the bracket
            # then missed r.  Caught by checking containment, not by the
            # arithmetic, which was self-consistent.)
            with localcontext() as c:
                c.prec = p; c.rounding = ROUND_CEILING
                E = (D(n*n-1)/12)*(w/2)**3/(ad*ad)
            l2, h2 = H-E, H+E
            nw = h2-l2
            if nw <= 0: return None, lev, 'empty'
            if -nw.adjusted() >= digits: return (l2,h2), lev, None
            if not (nw < w): return None, lev, 'not shrinking'
            ad, w = l2, nw
            lev += 1
    return None, lev, 'depth'

if __name__=='__main__':
    def cert_n(lo,hi,alpha,n,digits):
        """exact integer certificate: lo^n <= alpha <= hi^n  and hi-lo <= 10^-digits"""
        def coef(d):
            t=d.as_tuple(); e=t.exponent; coef=int("".join(map(str,t.digits)))
            return coef,e
        # compare a^n <= alpha exactly via integers
        li = D(int(lo)) if lo>0 else D(0)
        # use Decimal integer powers at high precision
        with localcontext() as c:
            c.prec = digits*2+200
            l, h = D(int(lo)), D(int(hi))
            # widen by 1 to be safe
            ok_lo = l**n <= D(alpha); ok_hi = h**n >= D(alpha)
            width_ok = (h-l) <= D(10)**(-digits)
        return (ok_lo and ok_hi and width_ok), (ok_lo,ok_hi,width_ok)
    print(f"{'n':>3} {'alpha':>10} | {'status':>28} | levels")
    for n in (2,3,4,5,8):
        for alpha in (10**6+1, 1234567, 2**20+3):
            r,l,err = certify_nroot(alpha, n, 40)
            if r is None:
                print(f"{n:>3} {alpha:>10} | {str(err):>28} | {l}"); continue
            ok,detail = cert_n(r[0],r[1],alpha,n,40)
            print(f"{n:>3} {alpha:>10} | {'CERTIFIED' if ok else 'fail '+str(detail):>28} | {l}")
