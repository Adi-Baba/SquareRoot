//! Certified cubic bracket-zoom square roots, in exact integer arithmetic.
//!
//! Every method returns a bracket `(lo, hi)` proved by `certify()` to satisfy
//! `lo^2 <= alpha <= hi^2` and `hi - lo <= 10^-digits`. Nothing here is
//! floating point and nothing is sampled.
//!
//! # Why this crate exists
//!
//! The Python audit (`Docs/Paper1.md` §10.1) reported that the zoom loses to
//! interval Newton, and attributed it to precision-weighted volume: 11.67
//! divisions and 85 operations per iteration for the zoom against 19.62 and
//! 72.75 for interval Newton. Measuring that Python implementation directly
//! showed the arithmetic -- divisions included -- is about 2% of wall clock;
//! the rest is CPython's `localcontext` dispatch. So the Python timings cannot
//! separate the methods, and the op-count proxy is not a cost model.
//!
//! Rust removes that confound: the arithmetic *is* the program. Division counts
//! (`Counters::divs`) are exact and machine-independent, and wall clock here is
//! comparable across methods.
//!
//! # The identity everything rests on
//!
//! Halley's step for `f(x) = x^2 - alpha` is
//!
//! ```text
//! H(x) = x(x^2 + 3 alpha) / (3 x^2 + alpha)
//! ```
//!
//! and it satisfies, *exactly*,
//!
//! ```text
//! H(x) - r = (x - r)^3 / (3 x^2 + r^2),        r = sqrt(alpha).
//! ```
//!
//! First, `sign(H(x) - r) = sign(x - r)`, so Halley is one-sided in the same
//! sense as the paper's corrected estimator `phi^2`: from `x > r` it
//! overshoots. Second, `|H(x) - r| = |x-r|^3/(3x^2+r^2)`, so any upper bound on
//! the distance `|x - r|` gives a certified cubic error bound.
//!
//! Internally the state is `(aa, ww)`: the left endpoint and the width, both
//! scaled by `10^scale`. Working scaled keeps every intermediate an integer and
//! makes the width condition exact.

use num_bigint::BigUint;
use num_traits::{One, Zero};

/// Counting of executed integer operations, per method. Exact.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Counters {
    pub divs: u64,
    pub muls: u64,
    pub adds: u64,
    pub levels: u64,
}

impl Counters {
    pub fn ops(&self) -> u64 {
        self.divs + self.muls + self.adds
    }
}

fn cdiv(n: &BigUint, d: &BigUint, c: &mut Counters) -> BigUint {
    c.divs += 1;
    n / d
}

/// Division rounded UP. Needed wherever a quantity is used as an UPPER bound:
/// a truncated bound is not a bound.
fn cdiv_up(n: &BigUint, d: &BigUint, c: &mut Counters) -> BigUint {
    c.divs += 1;
    let (q, r) = (n / d, n % d);
    if r.is_zero() {
        q
    } else {
        q + BigUint::one()
    }
}

fn cmul(a: &BigUint, b: &BigUint, c: &mut Counters) -> BigUint {
    c.muls += 1;
    a * b
}

fn cadd(a: &BigUint, b: &BigUint, c: &mut Counters) -> BigUint {
    c.adds += 1;
    a + b
}

pub fn pow10(d: u64) -> BigUint {
    BigUint::from(10u32).pow(d as u32)
}

/// Largest `a` with `a^2 <= alpha`.
pub fn isqrt(alpha: &BigUint) -> BigUint {
    if alpha.is_zero() {
        return BigUint::zero();
    }
    let mut x = BigUint::one();
    while (&x * &x) <= *alpha {
        x <<= 1;
    }
    let mut lo = x.clone() >> 1;
    let mut hi = x;
    loop {
        if lo >= hi {
            return lo;
        }
        let mid = (&lo + &hi + BigUint::one()) >> 1;
        if (&mid * &mid) <= *alpha {
            lo = mid;
        } else {
            hi = mid - BigUint::one();
        }
    }
}

/// Largest integer `a` with `a^n <= alpha`.
///
/// A proper binary search. An earlier version started at `1 << (bits/n + 2)`
/// and halved while `x^n > alpha`, which returns a **power of two**, not the
/// integer root: `iroot(101, 2)` gave `8`, so the seed bracket `[8, 9]` did not
/// contain `sqrt(101) ~ 10.05` and every subsequent step was working on the
/// wrong interval. Symptom: n-th-root certification failing with `hi^n < alpha`
/// while the square-root methods were fine.
pub fn iroot(alpha: &BigUint, n: u32) -> BigUint {
    if alpha.is_zero() {
        return BigUint::zero();
    }
    if *alpha == BigUint::one() {
        return BigUint::one();
    }
    // Upper bound: 2^ceil(bits/n) has nth power >= 2^bits > alpha.
    let mut hi = BigUint::one() << ((alpha.bits() + n as u64 - 1) / n as u64);
    let mut lo = BigUint::one();
    // Invariant: lo^n <= alpha < hi^n.  Find the largest lo with lo^n <= alpha.
    while &hi - &lo > BigUint::one() {
        let mid: BigUint = (&lo + &hi) >> 1;
        if mid.pow(n) <= *alpha {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo
}

/// `floor(log10(v))` for `v >= 1`; the number of integer decades in `v`.
pub fn decades(v: &BigUint) -> u64 {
    if v.is_zero() {
        return 0;
    }
    let mut d = ((v.bits() as f64) * std::f64::consts::LOG10_2).floor() as u64;
    while d > 0 && *v < pow10(d) {
        d -= 1;
    }
    while *v >= pow10(d + 1) {
        d += 1;
    }
    d
}

/// Prove, in exact integers, that the scaled bracket `(lo, hi)` certifies
/// `sqrt(alpha)` to `digits` decimals. `lo`, `hi` are scaled by `10^scale`.
pub fn certify(
    lo: &BigUint,
    hi: &BigUint,
    alpha: &BigUint,
    digits: u64,
    scale: u64,
) -> Result<(), String> {
    if lo > hi {
        return Err("lo > hi".into());
    }
    let at = alpha * pow10(2 * scale);
    if lo * lo > at {
        return Err("lo^2 > alpha".into());
    }
    if hi * hi < at {
        return Err("hi^2 < alpha".into());
    }
    let width = hi - lo;
    let tol = if scale >= digits {
        pow10(scale - digits)
    } else {
        BigUint::zero()
    };
    if width > tol {
        return Err(format!("width exceeds 10^-{}", digits));
    }
    Ok(())
}

/// Prove the same for an n-th root bracket scaled by `10^scale`.
pub fn certify_nth(
    lo: &BigUint,
    hi: &BigUint,
    alpha: &BigUint,
    n: u32,
    digits: u64,
    scale: u64,
) -> Result<(), String> {
    if lo > hi {
        return Err("lo > hi".into());
    }
    let at = alpha * pow10(n as u64 * scale);
    if lo.pow(n) > at {
        return Err("lo^n > alpha".into());
    }
    if hi.pow(n) < at {
        return Err("hi^n < alpha".into());
    }
    let width = hi - lo;
    let tol = if scale >= digits {
        pow10(scale - digits)
    } else {
        BigUint::zero()
    };
    if width > tol {
        return Err(format!("width exceeds 10^-{}", digits));
    }
    Ok(())
}

/// Outcome of one run.
pub struct Run {
    pub lo: BigUint,
    pub hi: BigUint,
    pub scale: u64,
    pub levels: u64,
    pub counters: Counters,
}

/// Halley's step `H(x) = x(x^2+3A)/(3x^2+A)` scaled by `10^scale`, returned as
/// `(floor, ceil)`. Every intermediate is rounded outward, which is why this
/// costs two divisions rather than one: a truncated step is not a bound.
/// `as2` is `alpha * 10^(2*scale)`, hoisted out of the loop by the caller
/// because `10^(2*scale)` is a big-integer power and recomputing it per level
/// dominated the run time.
///
/// Derivation of that scaling: with `X = x*S`, `S = 10^scale`,
/// `H = X(x^2+3A)/(3x^2+A) = X(X^2 + 3A S^2)/(3X^2 + A S^2)`. Using unscaled
/// `alpha` here silently produces a bracket that misses `r`.
fn halley(x: &BigUint, as2: &BigUint, c: &mut Counters) -> (BigUint, BigUint) {
    let x2 = cmul(x, x, c);
    let num = cmul(x, &cadd(&x2, &cmul(&BigUint::from(3u32), &as2, c), c), c);
    let den = cadd(&cmul(&BigUint::from(3u32), &x2, c), &as2, c);
    (cdiv(&num, &den, c), cdiv_up(&num, &den, c))
}

/// Certified bracket for `sqrt(alpha)` by Halley at the bracket **midpoint**.
///
/// # Soundness
///
/// `r` lies in `[a, a+w]` by the loop invariant, so `|m0 - r| <= w/2` exactly.
/// The identity gives
/// `|H(m0) - r| = |m0-r|^3/(3 m0^2 + r^2) <= (w/2)^3/(3 m0^2 + a^2)`,
/// using `r >= a`. So `E = w^3/(8(3 m0^2 + a^2))` is a valid upper bound and
/// `[H-E, H+E]` contains `r`.
///
/// Divisions per level: 2 (one for `H`, one for `E`).
pub fn halley_midpoint(alpha: &BigUint, digits: u64) -> Run {
    let mut c = Counters::default();
    let scale = digits + 20;
    let s = pow10(scale);
    let as2 = alpha * &pow10(2 * scale);
    let tol = pow10(scale - digits);
    let (mut aa, mut ww) = (&isqrt(alpha) * &s, s.clone());

    for _ in 0..4000u64 {
        c.levels += 1;
        let bb = &aa + &ww;
        let mm = (&aa + &bb) >> 1;
        let (hl, hu) = halley(&mm, &as2, &mut c);

        let w3 = cmul(&cmul(&ww, &ww, &mut c), &ww, &mut c);
        let m2 = cmul(&mm, &mm, &mut c);
        let a2 = cmul(&aa, &aa, &mut c);
        let den = cadd(&cmul(&BigUint::from(3u32), &m2, &mut c), &a2, &mut c);
        let e = cdiv_up(&w3, &cmul(&BigUint::from(8u32), &den, &mut c), &mut c);

        let lo = if hl > e { &hl - &e } else { BigUint::zero() };
        let hi = &hu + &e;
        let nw = &hi - &lo;
        // nw is already scaled by 10^scale, so the width condition
        // nw <= 10^(scale-digits) is decades(nw) >= scale - digits.
        // real width = nw / 10^scale, so width <= 10^-digits  <=>  nw <= 10^(scale-digits)
        // <=> floor(log10(nw)) <= scale - digits.
        if nw <= tol || nw >= ww {
            return Run { lo, hi, scale, levels: c.levels, counters: c };
        }
        aa = lo;
        ww = nw;
    }
    let hi = &aa + &ww;
    Run { lo: aa, hi, scale, levels: c.levels, counters: c }
}

/// Certified bracket by Halley at the **right endpoint** `b = a + w`.
///
/// # Soundness
///
/// `H(b) - r = (b-r)^3/(3b^2+r^2) >= 0`, so `H(b) >= r` and the new bracket is
/// one-sided: `[H(b) - E, H(b)]`. The plain estimate `m = a + p/(2a+w)`
/// satisfies `m < r` (the paper's `r - m = u(w-u)/(2a+w) > 0`), hence
/// `b - r < b - m`. Both effects then favour validity: the numerator is
/// larger and the denominator `3b^2 + a^2 <= 3b^2 + r^2` is smaller, so
/// `E = (b-m)^3/(3b^2+a^2) >= (b-r)^3/(3b^2+r^2) = H(b) - r`.
///
/// Divisions per level: 3.
pub fn halley_right(alpha: &BigUint, digits: u64) -> Run {
    let mut c = Counters::default();
    let scale = digits + 20;
    let s = pow10(scale);
    let as2 = alpha * &pow10(2 * scale);
    let tol = pow10(scale - digits);
    let (mut aa, mut ww) = (&isqrt(alpha) * &s, s.clone());

    for _ in 0..4000u64 {
        c.levels += 1;
        let bb = &aa + &ww;
        let a2 = cmul(&aa, &aa, &mut c);
        let ps = alpha * &s;
        let p = if ps > a2 { &ps - &a2 } else { BigUint::zero() };
        let tden = cadd(&cmul(&BigUint::from(2u32), &aa, &mut c), &ww, &mut c);
        let mq = cdiv(&p, &tden, &mut c);
        let mm = cadd(&aa, &mq, &mut c);

        let (hl, hu) = halley(&bb, &as2, &mut c);

        let d1 = if bb > mm { &bb - &mm } else { BigUint::zero() };
        let d1c = cmul(&cmul(&d1, &d1, &mut c), &d1, &mut c);
        let den = cadd(&cmul(&BigUint::from(3u32), &cmul(&bb, &bb, &mut c), &mut c), &a2, &mut c);
        let e = cdiv_up(&d1c, &den, &mut c);

        let lo = if hl > e { &hl - &e } else { BigUint::zero() };
        let hi = hu;
        let nw = &hi - &lo;
        // nw is already scaled by 10^scale, so the width condition
        // nw <= 10^(scale-digits) is decades(nw) >= scale - digits.
        // real width = nw / 10^scale, so width <= 10^-digits  <=>  nw <= 10^(scale-digits)
        // <=> floor(log10(nw)) <= scale - digits.
        if nw <= tol || nw >= ww {
            return Run { lo, hi, scale, levels: c.levels, counters: c };
        }
        aa = lo;
        ww = nw;
    }
    let hi = &aa + &ww;
    Run { lo: aa, hi, scale, levels: c.levels, counters: c }
}

/// The paper's corrected zoom in exact integers, for comparison.
///
/// `m = a + p/(2a+w)`, `est = a + p/(a+m)`, certified width `eps = w^3/(27a^2)`
/// (Theorem 4). Divisions per level: 3.
pub fn paper_zoom(alpha: &BigUint, digits: u64) -> Run {
    let mut c = Counters::default();
    let scale = digits + 20;
    let s = pow10(scale);
    let as2 = alpha * &pow10(2 * scale);
    let tol = pow10(scale - digits);
    let (mut aa, mut ww) = (&isqrt(alpha) * &s, s.clone());

    for _ in 0..4000u64 {
        c.levels += 1;
        let a2 = cmul(&aa, &aa, &mut c);
        // p = alpha - a^2, scaled by S^2.  BOTH terms are exact (alpha is a
        // short integer, aa is an exact integer times S), so p is exact too and
        // needs no interval.  An earlier version propagated p over the full
        // aa^2, which made est uncertain by ~aa/2 -- half the bracket -- and the
        // method could not contract at all.
        let p = &as2 - &a2;
        let tden = cadd(&cmul(&BigUint::from(2u32), &aa, &mut c), &ww, &mut c);

        // m = a + p/(2a+w): one division, so bound it both ways.
        let m_lo = cadd(&aa, &cdiv(&p, &tden, &mut c), &mut c);
        let m_hi = cadd(&aa, &cdiv_up(&p, &tden, &mut c), &mut c);

        // est = a + p/(a+m) DECREASES in m, so max est comes from m_lo.
        let am_lo = cadd(&aa, &m_lo, &mut c);
        let am_hi = cadd(&aa, &m_hi, &mut c);
        let est_hi = cadd(&aa, &cdiv_up(&p, &am_lo, &mut c), &mut c);
        let est_lo = cadd(&aa, &cdiv(&p, &am_hi, &mut c), &mut c);

        // eps = w^3/(27 a^2) is an upper bound on est - r (Theorem 4).
        let w3 = cmul(&cmul(&ww, &ww, &mut c), &ww, &mut c);
        let eps = cdiv_up(&w3, &cmul(&BigUint::from(27u32), &a2, &mut c), &mut c);

        // lo must use the SMALLEST est and the LARGEST eps; inflating both
        // pushes lo above r (found by certify, not by inspection).
        let lo = if est_lo > eps { &est_lo - &eps } else { BigUint::zero() };
        let hi = est_hi.clone();
        let nw = &hi - &lo;
        if nw <= tol || nw >= ww {
            return Run { lo, hi, scale, levels: c.levels, counters: c };
        }
        aa = lo;
        ww = nw;
    }
    let hi = &aa + &ww;
    Run { lo: aa, hi, scale, levels: c.levels, counters: c }
}

/// Bisection: certified, no divisions (the work is squarings).
pub fn bisect(alpha: &BigUint, digits: u64) -> Run {
    let mut c = Counters::default();
    let scale = digits + 20;
    let s = pow10(scale);
    let mut aa = &isqrt(alpha) * &s;
    let mut ww = s;
    let tol = pow10(scale - digits);
    while ww > tol {
        c.levels += 1;
        let mid = &aa + (&ww >> 1);
        if (&mid * &mid) <= alpha * &pow10(2 * scale) {
            aa = mid;
        }
        ww >>= 1;
    }
    let hi = &aa + &ww;
    Run { lo: aa, hi, scale, levels: c.levels, counters: c }
}

/// Certified bracket for `alpha^(1/n)` by Halley at the midpoint.
///
/// # The identity
///
/// For `f(x) = x^n - alpha` the Halley step is
///
/// ```text
/// H_n(x) = x[(n-1)x^n + (n+1)alpha] / [(n+1)x^n + (n-1)alpha]
/// ```
///
/// and a symbolic factorisation gives, for every `n >= 2`,
///
/// ```text
/// H_n(x) - r = (x-r)^3 Q_n(x,r) / D_n(x,r),   D_n = (n+1)x^n + (n-1)r^n,
/// ```
///
/// where `Q_n` is homogeneous of degree `n-2` with **non-negative coefficients**
/// and `Q_n(r,r) <= n(n^2-1)/6 * r^(n-2)`. Both structural facts were checked
/// symbolically for `n = 2..16`.
///
/// # The bound
///
/// `r` lies in `[a, b]` with `b = a + w`, so `r >= a` and `r <= b`. Monotonicity
/// of `Q_n` gives `Q_n(x,r) <= Q_n(b,b) <= n(n^2-1)/6 * b^(n-2)`, and
/// `D_n(x,r) >= 2n a^n`. With `|x - r| <= w/2` at the midpoint,
///
/// ```text
/// |H_n(x) - r| <= (w/2)^3 (n^2-1) b^(n-2) / (12 a^n) =: E.
/// ```
///
/// In scaled units (`w = ww/S`, `b = bb/S`, `a = aa/S`) the `S` factors cancel:
///
/// ```text
/// E_scaled = (n^2-1) ww^3 bb^(n-2) / (96 aa^n).
/// ```
///
/// At `n = 2` this reduces to `w^3/(32 a^2)`, exactly the square-root midpoint
/// bound — a useful consistency check. Getting this wrong by an extra `S^(n-2)`
/// silently makes `E` too small for every `n >= 3`.
///
/// # Why this matters
///
/// `Paper2.md` §9.4 proves the paper's own map `phi_n` is **repelling** for
/// `n >= 3` (`lambda_n -> -(n-1)`). Halley's multiplier at the root is `0` for
/// every `n`, so a bracket-carried cubic contraction exists there regardless.
pub fn halley_nth_root(alpha: &BigUint, n: u32, digits: u64) -> Run {
    match halley_nth_root_inner(alpha, n, digits, false) {
        Ok(r) => r,
        Err(e) => panic!("halley_nth_root produced an invalid bracket: {}", e),
    }
}

/// As `halley_nth_root`, but returns the level at which the bracket first failed
/// to enclose `alpha^(1/n)` instead of panicking. Used by the tests: a certified
/// method that is wrong at level 3 but only *reported* wrong at the end is much
/// harder to debug than one that stops at level 3.
pub fn halley_nth_root_checked(alpha: &BigUint, n: u32, digits: u64) -> Result<Run, String> {
    halley_nth_root_inner(alpha, n, digits, true)
}

fn halley_nth_root_inner(alpha: &BigUint, n: u32, digits: u64, verify: bool) -> Result<Run, String> {
    if n < 2 {
        return Err(format!("n must be >= 2, got {}", n));
    }
    let mut c = Counters::default();
    let scale = digits + 20;
    let s = pow10(scale);
    let asn = alpha * &pow10(n as u64 * scale);
    let at = alpha * &pow10(n as u64 * scale);
    // H_n(x) = x[(n-1)x^n + (n+1)A] / [(n+1)x^n + (n-1)A]
    //          ^ n-1 ^ n+1        ^ n+1      ^ n-1
    // An earlier version used `n` in place of `n+1` in the numerator's alpha
    // term AND in the denominator's x^n term, which is not a Halley step at
    // all: it does not even fix r.  For n = 2 it gave the wrong denominator
    // 2x^2 + A instead of 3x^2 + A, and the bracket silently drifted off the
    // root from level 1.  Caught by the per-level invariant check.
    let np1 = BigUint::from(n + 1);
    let nm1 = BigUint::from(n - 1);
    let n2m1 = BigUint::from((n as u64) * (n as u64) - 1);

    let tol = pow10(scale - digits);
    let mut aa = &iroot(alpha, n) * &s;
    let mut ww = s.clone();

    // First, narrow the bracket with bisection until the Halley error bound
    // permits contraction. The bound requires E <= w/2, i.e. (n^2-1) b^(n-2) w^2
    // <= 48 a^n. Starting from [a, a+1], this fails for large n. Bisection is
    // certified, division-free, and only costs squarings.
    let mut a_pow_n = aa.pow(n);
    let mut b = &aa + &ww;
    while ww > tol {
        let n2m1 = BigUint::from((n as u64) * (n as u64) - 1);
        let b_pow_nm2 = b.pow(n.saturating_sub(2));
        let lhs = &n2m1 * &b_pow_nm2 * &ww * &ww;
        let rhs = BigUint::from(48u32) * &a_pow_n;
        if lhs <= rhs {
            break;
        }
        // One bisection step: m = (a+b)/2.  mid^n <= A?  If so a = m, else b = m.
        let mid: BigUint = (&aa + &b) >> 1;
        let mid_pow_n = mid.pow(n);
        if mid_pow_n <= at {
            aa = mid;
            a_pow_n = mid_pow_n;
        } else {
            b = mid;
        }
        ww = &b - &aa;
    }

    for lev in 0..4000u64 {
        c.levels += 1;
        let bb = &aa + &ww;
        let mm: BigUint = (&aa + &bb) >> 1;
        let mn = mm.pow(n);
        let num = &mm * (&(&nm1 * &mn) + &(&np1 * &asn));
        let den = &(&np1 * &mn) + &(&nm1 * &asn);
        let hl = cdiv(&num, &den, &mut c);
        let hu = cdiv_up(&num, &den, &mut c);

        let w3 = cmul(&cmul(&ww, &ww, &mut c), &ww, &mut c);
        let numer = cmul(&n2m1, &cmul(&w3, &bb.pow(n - 2), &mut c), &mut c);
        let denom = cmul(&BigUint::from(96u32), &aa.pow(n), &mut c);
        let e = cdiv_up(&numer, &denom, &mut c);

        let lo = if hl > e { &hl - &e } else { BigUint::zero() };
        let hi = &hu + &e;

        if verify {
            if lo.pow(n) > at {
                return Err(format!("level {}: lo^n > alpha (lo overshoots)", lev));
            }
            if hi.pow(n) < at {
                return Err(format!("level {}: hi^n < alpha (E or hu too small)", lev));
            }
        }

        let nw = &hi - &lo;
        if nw <= tol || nw >= ww {
            return Ok(Run { lo, hi, scale, levels: c.levels, counters: c });
        }
        aa = lo;
        ww = nw;
    }
    Err("depth exceeded".into())
}

#[cfg(test)]
mod tests {
    use super::*;

    const ALPHAS: [u64; 12] = [2, 3, 5, 7, 10, 13, 101, 257, 1009, 7919, 104729, 999983];

    #[test]
    fn isqrt_is_exact() {
        for a in ALPHAS {
            let al = BigUint::from(a);
            let r = isqrt(&al);
            assert!(r.clone() * r.clone() <= al, "isqrt({}) too large", a);
            assert!((r.clone() + BigUint::one()).pow(2) > al, "isqrt({}) too small", a);
        }
    }

    #[test]
    fn iroot_is_exact() {
        // Regression: an earlier iroot halved a power of two and returned 8 for
        // 101, so the seed bracket did not contain the root at all.
        for n in 2u32..=8 {
            for a in ALPHAS {
                let al = BigUint::from(a);
                let r = iroot(&al, n);
                assert!(r.pow(n) <= al, "iroot({},{}) too large", a, n);
                assert!(
                    (r.clone() + BigUint::one()).pow(n) > al,
                    "iroot({},{}) too small", a, n
                );
            }
        }
    }

    #[test]
    fn decades_is_floor_log10() {
        for d in 0u64..40 {
            assert_eq!(decades(&pow10(d)), d, "decades(10^{})", d);
            assert_eq!(decades(&(&pow10(d) + BigUint::one())), d);
        }
        assert_eq!(decades(&BigUint::one()), 0);
        assert_eq!(decades(&BigUint::zero()), 0);
    }

    fn check_sq(f: fn(&BigUint, u64) -> Run, name: &str, digits: u64) {
        for a in ALPHAS {
            let al = BigUint::from(a);
            let r = f(&al, digits);
            certify(&r.lo, &r.hi, &al, digits, r.scale)
                .unwrap_or_else(|e| panic!("{} alpha={} digits={}: {}", name, a, digits, e));
        }
    }

    #[test]
    fn square_root_methods_certify() {
        for d in [20u64, 60, 200] {
            check_sq(halley_midpoint, "halley_midpoint", d);
            check_sq(halley_right, "halley_right", d);
            check_sq(bisect, "bisect", d);
        }
    }

    /// Whether the certified bound can contract from the naive integer seed.
    ///
    /// `E <= w/2` requires `(n^2-1) b^(n-2) w^2 <= 48 a^n`. With the seed
    /// bracket `[a, a+w]`, `a = iroot(alpha, n)` and `w = 1`, this fails for
    /// large `n` and small `alpha`: for `n = 4, alpha = 2` it reads
    /// `15 * 2^2 * 1 <= 48 * 1^4`, i.e. `60 <= 48`, which is false. The bound
    /// is then wider than the bracket and the iteration cannot start. That is a
    /// real limitation of the *bound*, not of the implementation: it needs a
    /// seed of relative width well below `a^(2/n)`, which the naive integer root
    /// does not supply.
    fn can_contract(alpha: &BigUint, n: u32) -> bool {
        let a = iroot(alpha, n);
        let b = &a + BigUint::one();
        // (n^2-1) b^(n-2) * 1 <= 48 a^n, all integer-exact
        BigUint::from((n as u64) * (n as u64) - 1) * b.pow(n - 2) <= BigUint::from(48u32) * a.pow(n)
    }

    #[test]
    fn nth_root_bracket_always_contains_the_root() {
        // The invariant that must hold unconditionally: whatever the width,
        // the returned bracket encloses alpha^(1/n), checked at EVERY level.
        for n in [2u32, 3, 4, 5, 8, 16] {
            for a in ALPHAS {
                let al = BigUint::from(a);
                halley_nth_root_checked(&al, n, 100)
                    .unwrap_or_else(|e| panic!("n={} alpha={}: {}", n, a, e));
            }
        }
    }

    #[test]
    fn nth_root_reaches_full_width_where_the_bound_permits() {
        // Only where the bound can actually contract; see `can_contract`.
        let mut checked = 0usize;
        for n in [2u32, 3, 4, 5, 8, 16] {
            for a in ALPHAS {
                let al = BigUint::from(a);
                if !can_contract(&al, n) {
                    continue;
                }
                checked += 1;
                let r = halley_nth_root_checked(&al, n, 100)
                    .unwrap_or_else(|e| panic!("n={} alpha={}: {}", n, a, e));
                certify_nth(&r.lo, &r.hi, &al, n, 100, r.scale)
                    .unwrap_or_else(|e| panic!("n={} alpha={} certify: {}", n, a, e));
            }
        }
        assert!(checked > 0, "no (n, alpha) pair could contract -- bound too loose");
    }

    #[test]
    fn can_contract_matches_observation() {
        // The exclusion predicate must actually predict what the method can do,
        // or it is just an excuse.  Checked against the observed contraction.
        for n in [2u32, 3, 4, 5, 8, 16] {
            for a in ALPHAS {
                let al = BigUint::from(a);
                let r = halley_nth_root_checked(&al, n, 100).unwrap();
                let width = &r.hi - &r.lo;
                let reached = width <= pow10(r.scale - 100);
                let predicted = can_contract(&al, n);
                assert!(
                    !predicted || reached,
                    "n={} alpha={}: predicted contraction but width not reached",
                    n, a
                );
            }
        }
    }

    #[test]
    fn nth_root_rejects_bad_n() {
        let al = BigUint::from(2u64);
        assert!(halley_nth_root_checked(&al, 1, 50).is_err());
    }

    #[test]
    fn counters_are_plausible() {
        let al = BigUint::from(1000003u64);
        let r = halley_midpoint(&al, 200);
        assert!(r.levels >= 3, "cubic method should take several levels");
        // 3 divisions per level: two for H (a floor and a ceil, since a
        // truncated step is not a bound) and one for E.
        assert_eq!(r.counters.divs, 3 * r.counters.levels,
                   "halley_midpoint: 3 divisions/level, got {} over {} levels",
                   r.counters.divs, r.counters.levels);
    }
}
