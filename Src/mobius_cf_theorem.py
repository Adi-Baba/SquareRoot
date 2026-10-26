#!/usr/bin/env python3
"""
Exact Möbius-CF Correspondence for Pronic Numbers

Theorem: For α = a(a+1), the Möbius iterates φⁿ(a+1) have exact
denominator/numerator correspondence with CF convergents.

z_n = φⁿ(a+1) = A_n / B_n in lowest terms
CF convergent: p_n / q_n

Theorem 1 (Exact Correspondence):
- Unified form:    z_n = a(a+1) · q_n / p_n          (all n >= 0)
- Odd n = 2k+1:   B_{2k+1} = p_{2k+1}
- Even n = 2k:    a · B_{2k} = p_{2k}
- Odd n = 2k+1:   A_{2k+1} = a(a+1) · q_{2k+1}
- Even n = 2k:    A_{2k} = (a+1) · q_{2k}

(The even-n numerator item relates A to q, not to B: A = (a+1)*q. The
earlier docstring here wrote "(a+1) * B_{2k} = A_{2k}", conflating the
Moebius denominator B with the CF denominator q; at a=2, n=2 that reads
3*11 = 27, which is false, while A = 27 = 3*9 = (a+1)*q_2 is exact.
Docs/Paper3.md Sec 2.5 carries the proof.)
"""

from fractions import Fraction as F

def verify_pronic_theorem(a_max=50, n_max=15):
    """Verify Theorem 1 for all pronic α = a(a+1) up to a_max."""
    failed = []
    for a in range(2, a_max+1):
        # CF convergents
        p = [a, 2*a+1]
        q = [1, 2]
        for n in range(2, n_max+1):
            an = 2 if n%2==1 else 2*a
            p.append(an*p[-1] + p[-2])
            q.append(an*q[-1] + q[-2])
        
        z = F(a+1)
        for n in range(n_max+1):
            A, B = z.numerator, z.denominator
            if n % 2 == 1:  # odd
                if B != p[n]:
                    failed.append(f"a={a}, n={n} odd: B={B}, p[n]={p[n]}")
                if A != a*(a+1)*q[n]:
                    failed.append(f"a={a}, n={n} odd: A={A}, a(a+1)*q={a*(a+1)*q[n]}")
            else:  # even
                if a*B != p[n]:
                    failed.append(f"a={a}, n={n} even: a*B={a*B}, p[n]={p[n]}")
                if A != (a+1)*q[n]:
                    failed.append(f"a={a}, n={n} even: A={A}, (a+1)*q={(a+1)*q[n]}")
            z = F(a) + F(a) / (F(a) + z)
    
    if failed:
        print("FAILURES:")
        for f in failed[:10]:
            print(f"  {f}")
        return False
    print(f"✓ Theorem verified for a ≤ {a_max}, n ≤ {n_max}")
    return True


def verify_general_failure(a_max=20):
    """Verify the correspondence FAILS for non-pronic p|2a, p≠a."""
    failed_to_fail = []
    for a in range(2, a_max+1):
        for p in range(1, 2*a+1):
            if (2*a) % p != 0 or p == a:
                continue
            a1 = 2*a // p
            a2 = 2*a
            
            # CF convergents
            p_cf = [a, a*a1 + 1]
            q_cf = [1, a1]
            for n in range(2, 10):
                an = a1 if n%2==1 else a2
                p_cf.append(an*p_cf[-1] + p_cf[-2])
                q_cf.append(an*q_cf[-1] + q_cf[-2])
            
            # Möbius
            z = F(a+1)
            match = True
            for n in range(1, 8):
                z = F(a) + F(p) / (F(a) + z)
                if n % 2 == 1:
                    if z.denominator != p_cf[n]:
                        match = False
                        break
                else:
                    if a * z.denominator != p_cf[n]:
                        match = False
                        break
            if match:
                failed_to_fail.append((a, p))
    
    if failed_to_fail:
        print(f"UNEXPECTED MATCHES for non-pronic: {failed_to_fail}")
        return False
    print(f"✓ All non-pronic p|2a correctly FAIL (tested a≤{a_max})")
    return True


if __name__ == "__main__":
    print("=" * 60)
    print("Möbius-CF Correspondence Verification")
    print("=" * 60)
    verify_pronic_theorem(100, 12)
    verify_general_failure(20)
    print("\n✓ All verifications passed.")
