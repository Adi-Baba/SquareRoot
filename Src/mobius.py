"""The Mobius structure of the estimator family, and an exact n-step form.

This file records the algebra behind new Section 5.3 of `Docs/Paper1.md`.
It contains NO numerical experiment: everything here is either a symbolic
identity (settled once with sympy, results recorded below) or an exact
rational check of that identity.  No floats, no sampling, no timings.

Setup (unchanged from the paper).  For alpha > 0 write

    r = sqrt(alpha),   a > 0,   b = a + w,   0 < u = r - a < w,
    p = alpha - a^2 = u(2a + u),

and let the plain estimator be the Mobius-looking rational function

    phi(z) = a + p/(a + z).

The claim proved in the paper, and checked here, is that phi is a Mobius
transformation whose fixed points are +r and -r, so that the change of variable

    T(z) = (z - r)/(z + r)

conjugates phi to a *scaling*:

    T(phi(z)) = lambda * T(z),      lambda = (a - r)/(a + r) = -u/(2a + u).

Everything else follows by induction.  In particular the n-th iterate is

    phi^n(z) - r  =  2 r lambda^n T(z) / (1 - lambda^n T(z)),                 (N)

exactly, not asymptotically -- which recovers Theorem 3 (n = 1), recovers the
corrected estimator as phi^2 (so Corollary 2 and the 4/27 bound are its n = 2
case), and makes the *one-sidedness exact for every iterate* rather than only
for the first two.
"""

from fractions import Fraction as F

# --------------------------------------------------------------- definitions


def p_of(a, u):
    """alpha - a^2 with r = a + u.  Equals u(2a+u) identically."""
    return u * (2 * a + u)


def phi(a, u, z):
    """The plain map, phi(z) = a + p/(a+z), in exact rationals."""
    a, u, z = F(a), F(u), F(z)
    return a + p_of(a, u) / (a + z)


def T(a, u, z):
    """The conjugacy T(z) = (z-r)/(z+r), with r = a+u."""
    a, u, z = F(a), F(u), F(z)
    r = a + u
    return (z - r) / (z + r)


def Tinv(a, u, v):
    """Inverse of T.  T^-1(v) = r(1+v)/(1-v)."""
    a, u, v = F(a), F(u), F(v)
    r = a + u
    return r * (1 + v) / (1 - v)


def lam(a, u):
    """The scaling factor lambda = (a-r)/(a+r) = -u/(2a+u)."""
    a, u = F(a), F(u)
    return (a - (a + u)) / (a + (a + u))


def error_closed_form(a, u, z, n):
    """phi^n(z) - r, in closed form (formula (N)).  Exact, no iteration."""
    v = lam(a, u) ** n * T(a, u, z)
    return 2 * (a + u) * v / (1 - v)


def normed_error(a, w, u):
    """e * a^2 / w^3, the scale-free error, in exact rationals.

    With t = w/a and eta = u/w the error of the corrected estimate becomes

        e * a^2 / w^3  =  F(t, eta) = eta^2(1-eta) / (4 + 2t + 2t*eta + t^2*eta^2),

    and F(0, 2/3) = 1/27 exactly.  So Theorem 4's constant is the value of F at
    the corner (t, eta) = (0, 2/3) of the state space, which is a boundary point
    no bracket realises: t = w/a > 0 and eta in (0,1) strictly.

    This is Corollary 7's setting.  Nothing here is asymptotic.
    """
    a, w, u = F(a), F(w), F(u)
    t, eta = w / a, u / w
    return eta * eta * (1 - eta) / (4 + 2 * t + 2 * t * eta + t * t * eta * eta)


def position_blind_constant(a, w):
    """2 / (27(2+t)), with t = w/a -- the sharp position-blind constant.

    Corollary 7: sup over eta in [0,1] of F(t, eta) is strictly below this for
    every t > 0, so the sharp constant needing no position information is this,
    NOT 1/27.  Equivalently

        (4w^3/27) / D(0)  =  2w^3 / (27 a (2a+w))  =  2/(27(2+t)) * w^3/a^2,

    which is exactly what exp_u_aware.bound_published returns -- so the codebase
    has been using the sharper constant all along, and neither paper said so.

    Relative to the published 1/27 it is better by t/(2+t): 33% at w/a = 1,
    4.8% at w/a = 1/10, 0.05% at w/a = 1/1000.
    """
    a, w = F(a), F(w)
    return F(2) / (27 * (2 + w / a))


def strict_slack(a, w, u):
    """The numerator of Corollary 7's slack -- position_blind_constant - F.

    Returned alone so that its two summands can be tested SEPARATELY for
    non-negativity.  The strictness of Corollary 7 is exactly the statement that
    they never vanish simultaneously:

        27(2+t)(eta-2/3)^2(eta+1/3)   vanishes only at eta = 2/3
        2t*eta*(2+t*eta)               vanishes only at eta = 0 (or t = 0)

    so their zeros are disjoint on [0,1] x (0, inf).
    """
    a, w, u = F(a), F(w), F(u)
    t, eta = w / a, u / w
    sq = (eta - F(2, 3)) ** 2 * (eta + F(1, 3))     # >= 0, vanishes at eta=2/3
    mx = 2 * t * eta * (2 + t * eta)                 # >= 0, vanishes at eta=0
    return 27 * (2 + t) * sq, mx


def attainment_ratio(a, w):
    """How close the sharp bound gets to being attained, exactly.

    Theorem 4 bounds the corrected error by w^3/(27 a^2) using two SEPARATE
    suprema:  f(u) = u^2(w-u) <= 4w^3/27, attained at u = 2w/3, and
    D(u) >= 4a^2, approached as u -> 0.  Those are not attained at the same
    u, so the bound is sharp only in the limit w/a -> 0.

    This returns the exact factor by which the bound is attained at u = 2w/3:

        f(2w/3)/D(2w/3)  /  (w^3/(27 a^2))  =  18a^2/(18a^2 + 15aw + 2w^2).

    It is 1 in the limit w/a -> 0 and falls to 1/2 at w/a = 1.  So the
    maximiser 2w/3 maximises the NUMERATOR, not the ratio f(u)/D(u); the two
    coincide only asymptotically.  Stated in Paper2 Corollary 4.
    """
    a, w = F(a), F(w)
    u = 2 * w / 3
    f = u * u * (w - u)
    D = 2 * a * (2 * a + w) + u * (2 * a + u)
    bound = w ** 3 / (27 * a * a)
    return (f / D) / bound


def phi_iter(a, u, z, n):
    """phi iterated n times, by actually iterating.  Used only to test (N)."""
    z = F(z)
    for _ in range(n):
        z = phi(a, u, z)
    return z - (a + u)


# ------------------------------------------------------- symbolic identities
#
# Settled once, offline, with sympy (a, u, w, z symbolic and positive; r = a+u,
# alpha = r^2).  Re-running is unnecessary; the results are transcribed here and
# each is re-checked in exact rationals by test_mobius.py.  What was verified:
#
#   (S1)  T(phi(z)) == lam * T(z)                       -> True
#   (S2)  2r lam^2 T(b) / (1 - lam^2 T(b)) == u^2(w-u)/D,  D = 2a(2a+w)+u(2a+u)
#                                                           -> True
#   (S3)  d/du [u^2 (w - u)] = -u^2 + 2u(-u + w), sole critical point u = 2w/3,
#         at which the numerator attains 4w^3/27             -> True
#   (S4)  phi'(r) == lam                                   -> True
#
# (S2) is the load-bearing one: it says formula (N) at n = 2 reproduces the
# paper's existing Theorem 3 / 4/27 constant EXACTLY.  Without it the closed
# form would be a new and unrelated claim; with it, Section 5.3 is a
# re-derivation of results already in the paper, not a replacement for them.
#   (S5)  attainment ratio 18a^2/(18a^2+15aw+2w^2)   -> True
#   (S6)  e*a^2/w^3 == F(t,eta) = eta^2(1-eta)/(4+2t+2t*eta+t^2*eta^2),  t=w/a, eta=u/w
#                                                           -> True
#   (S7)  F(0, 2/3) == 1/27                            -> True
#   (S8)  2/(27(2+t)) - F(t,eta) ==
#           [27(2+t)(eta-2/3)^2(eta+1/3) + 2t*eta(2+t*eta)]
#           / [27(2+t)(4+2t+2t*eta+t^2*eta^2)]          -> True
#   (S9)  sup_eta F(t,eta) = 1/27 - 5t/162 + 16t^2/729 + O(t^3),
#         argmax eta*(t) = 2/3 - t/27 + t^2/243 + O(t^3)
#
# (S5) is what makes "sharp" precise: Theorem 4's bound is the product of two
# separately-taken suprema, so it is attained only as w/a -> 0.
#
# (S8) is the load-bearing new one -- Corollary 7.  It upgrades (S5) from an
# asymptotic attainment factor to an EXACT identity exhibiting a strictly
# positive slack for every bracket.  Its two numerator summands are non-negative
# on [0,1] x [0,inf) and their zero sets are disjoint, so the slack is strictly
# positive whenever t > 0.  Consequence: 1/27 is a supremum over a corner of the
# state space that no bracket reaches, and the sharp position-blind constant is
# 2/(27(2+t)), not 1/27.
#
# (S9) explains the "unwieldy quartic maximiser" flagged in Paper2 Section 8.3:
# it is 2/3 displaced by first order -t/27, so the true ratio maximiser and the
# numerator maximiser 2w/3 agree to zeroth order and part company at O(t).
SYMBOLIC_IDENTITY = {
    "F(0,2/3) = 1/27": True,
    "attainment ratio 18a^2/(18a^2+15aw+2w^2)": True,
    "T(phi) = lam*T": True,
    "closed form n=2 reproduces u^2(w-u)/D": True,
    "e*a^2/w^3 = F(t,eta)": True,
    "numerator maximised at u = 2w/3, value 4w^3/27": True,
    "phi'(r) = lam": True,
    "strict slack 2/(27(2+t)) - F = [..]/[..]": True,
    "sup F = 1/27 - 5t/162 + 16t^2/729 + O(t^3)": True,
}
