//! Exact squarefree decomposition for the SRS sign path.
//!
//! Given `n`, find `c` and squarefree `s` with `n = c^2 * s`.
//!
//! WHY THIS EXISTS SEPARATELY FROM `field::is_squarefree`
//! -------------------------------------------------------
//! The field constructor only ever needs a yes/no answer, and it rejects
//! anything it cannot settle.  The sign path needs the decomposition itself:
//! aggregating coefficients by squarefree part is what makes the zero test
//! exact, because `{sqrt(d) : d distinct squarefree}` is linearly independent
//! over Q.  So this module returns the pair, and verifies the identity
//! `c * c * s == n` before returning, which means a misfactorisation cannot
//! pass silently.
//!
//! The cost profile matters here and drove the design.  The F1 family builds
//! radicands `p^2` and `alpha * q^2` with `p, q` of up to 10^7 bits.  Neither
//! is factorable in any useful time, and neither needs to be: `p^2` is
//! detected by a single integer square root, and `alpha * q^2` only needs the
//! (tiny) factorisation of `alpha`.  So the perfect-square test comes first
//! and the general path is only reached for genuinely awkward inputs.

use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

use crate::error::{Result, SrsError};

/// Trial-division bound. Generators and small cofactors live far below this.
const TRIAL_DIVISION_BOUND: u64 = 1 << 20;


/// `a` divides `b`. Named because num-bigint 0.4 has no `is_multiple_of`.
fn is_multiple_of(a: &BigInt, b: &BigInt) -> bool {
    !b.is_zero() && (a % b).is_zero()
}

#[inline]
fn rem_is_multiple_of(a: &BigInt, b: &BigInt) -> bool {
    is_multiple_of(a, b)
}

/// `n = c^2 * s` with `s` squarefree.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Decomposition {
    pub c: BigInt,
    pub s: BigInt,
}

/// Decompose `n` exactly.
///
/// Returns an error rather than a guess if the remaining cofactor is too
/// large to split.  The caller is deciding a sign, so an approximate answer
/// would be worse than none.
pub fn decompose(n: &BigInt) -> Result<Decomposition> {
    if *n < BigInt::one() {
        return Err(SrsError::InvalidInput(format!(
            "squarefree decomposition needs n >= 1, got {n}"
        )));
    }
    if *n == BigInt::one() {
        return Ok(Decomposition { c: BigInt::one(), s: BigInt::one() });
    }

    // Fast path: a perfect square needs no factoring at all. This is the F1
    // radicand p^2, and at 10^7 bits it is the difference between instant and
    // impossible.
    let root = n.sqrt();
    if &root * &root == *n {
        return Ok(Decomposition { c: root, s: BigInt::one() });
    }

    let (mut c, mut s, remaining) = trial_division(n)?;

    if remaining == BigInt::one() {
        return finish(c, s, n);
    }

    let r = remaining.sqrt();
    if &r * &r == remaining {
        c *= r;
        return finish(c, s, n);
    }

    if is_probable_prime(&remaining) {
        s *= remaining;
        return finish(c, s, n);
    }

    // General case: split the cofactor. Factoring is genuinely needed here.
    let mut exps: Vec<(BigInt, u32)> = Vec::new();
    factor(&remaining, &mut exps, 0)?;
    for (p, e) in exps {
        c *= p.pow(e / 2);
        if e % 2 == 1 {
            s *= p;
        }
    }
    finish(c, s, n)
}

/// Verify the identity and hand back the pair.
fn finish(c: BigInt, s: BigInt, n: &BigInt) -> Result<Decomposition> {
    if &c * &c * &s != *n {
        return Err(SrsError::ArithmeticError(format!(
            "squarefree decomposition is inconsistent for {n}: c^2*s != n"
        )));
    }
    if !is_squarefree(&s)? {
        return Err(SrsError::ArithmeticError(format!(
            "squarefree decomposition produced a non-squarefree part {s} for {n}"
        )));
    }
    Ok(Decomposition { c, s })
}

/// Strip small prime factors, returning (square part, squarefree part so far,
/// remaining cofactor).
fn trial_division(n: &BigInt) -> Result<(BigInt, BigInt, BigInt)> {
    let (c, s, remaining, _) = trial_division_report(n)?;
    Ok((c, s, remaining))
}

/// As `trial_division`, but also reports whether any small prime occurred with
/// an even exponent >= 2 during the scan.
///
/// The parity flag is what `is_squarefree` needs, and it cannot be recovered
/// from the returned triple: trial division puts p^2 into `c` and nothing into
/// `s`, so 18 reduces to (c=3, s=1, cofactor=1) and every output field looks
/// squarefree even though 9 divides it. The earlier `is_squarefree` read the
/// cofactor and got 18 wrong for exactly that reason.
fn trial_division_report(n: &BigInt) -> Result<(BigInt, BigInt, BigInt, bool)> {
    let mut c = BigInt::one();
    let mut s = BigInt::one();
    let mut remaining = n.clone();
    let mut repeated_small_prime = false;

    let mut p: u64 = 2;
    while p <= TRIAL_DIVISION_BOUND {
        let bp = BigInt::from(p);
        if &bp * &bp > remaining {
            return Ok((c, s, remaining, repeated_small_prime));
        }
        if rem_is_multiple_of(&remaining, &bp) {
            let mut e: u32 = 0;
            while rem_is_multiple_of(&remaining, &bp) {
                remaining /= &bp;
                e += 1;
            }
            if e >= 2 {
                repeated_small_prime = true;
            }
            c *= bp.pow(e / 2);
            if e % 2 == 1 {
                s *= bp;
            }
        }
        p += if p == 2 { 1 } else { 2 };   // 2, then odds only
    }

    // Out of trial-division range but the cofactor may already be small.
    let limit = BigInt::from(TRIAL_DIVISION_BOUND) * BigInt::from(TRIAL_DIVISION_BOUND);
    if remaining <= limit {
        return Ok((c, s, remaining, repeated_small_prime));
    }

    // Not a failure: hand the large cofactor back and let the caller decide.
    Ok((c, s, remaining, repeated_small_prime))
}

/// Accumulate the prime factorisation of `n` into `out`.
fn factor(n: &BigInt, out: &mut Vec<(BigInt, u32)>, depth: u32) -> Result<()> {
    if *n == BigInt::one() {
        return Ok(());
    }
    if depth > 64 {
        return Err(SrsError::ArithmeticError(
            "factorisation recursion too deep".into(),
        ));
    }
    if is_probable_prime(n) {
        push_factor(out, n.clone(), 1);
        return Ok(());
    }
    let d = pollard_rho(n)?;
    let other = n / &d;
    factor(&d, out, depth + 1)?;
    factor(&other, out, depth + 1)?;
    Ok(())
}

fn push_factor(out: &mut Vec<(BigInt, u32)>, p: BigInt, e: u32) {
    for slot in out.iter_mut() {
        if slot.0 == p {
            slot.1 += e;
            return;
        }
    }
    out.push((p, e));
}

fn gcd(a: &BigInt, b: &BigInt) -> BigInt {
    let mut x = a.clone();
    let mut y = b.clone();
    while !y.is_zero() {
        let t = x % &y;
        x = y;
        y = t;
    }
    x.abs()
}

/// Brent's variant of Pollard's rho, deterministic in `n`.
fn pollard_rho(n: &BigInt) -> Result<BigInt> {
    if is_multiple_of(n, &BigInt::from(2u8)) {
        return Ok(BigInt::from(2u8));
    }
    if is_multiple_of(n, &BigInt::from(3u8)) {
        return Ok(BigInt::from(3u8));
    }

    // rho needs about n^(1/4) iterations; scale the budget with the input so
    // a merely-large cofactor is not rejected as though it were hopeless.
    let budget: u64 = std::cmp::max(1 << 12, (n.bits() / 4) * 64);
    let two = BigInt::from(2u8);

    for attempt in 1..=4u64 {
        let one = BigInt::one();
        let mut y = one.clone() + BigInt::from(attempt);
        let c = one.clone() + BigInt::from(attempt);
        let mut g = one.clone();
        let mut q = one.clone();
        let mut r: u64 = 1;
        let mut ys = y.clone();
        let mut x = y.clone();
        let mut spent: u64 = 0;

        while g.is_one() && spent < budget {
            x = y.clone();
            for _ in 0..r {
                y = (&y * &y + &c) % n;
                spent += 1;
            }
            let mut k: u64 = 0;
            while k < r && g.is_one() {
                ys = y.clone();
                for _ in 0..std::cmp::min(64, r - k) {
                    y = (&y * &y + &c) % n;
                    q = (&q * (&x - &y).abs()) % n;
                }
                g = gcd(&q, n);
                k += 64;
            }
            r *= 2;
        }

        if g == *n {
            // The batch collapsed onto n; back up one difference at a time.
            g = BigInt::one();
            while g.is_one() {
                ys = (&ys * &ys + &c) % n;
                g = gcd(&(&x - &ys).abs(), n);
            }
        }

        if g > one && g < *n {
            return Ok(g);
        }
        let _ = &two;
    }

    Err(SrsError::FactorizationFailed(format!(
        "Pollard's rho could not split a {}-bit cofactor within {budget} rounds \
         x 4 restarts; it needs Pollard p-1, ECM, or a pre-factored input",
        n.bits()
    )))
}

/// Deterministic Miller-Rabin below 2^64, probabilistic above.
pub fn is_probable_prime(n: &BigInt) -> bool {
    if *n < BigInt::from(2u8) {
        return false;
    }
    const SMALL: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    for p in SMALL {
        let bp = BigInt::from(p);
        if is_multiple_of(n, &bp) {
            return *n == bp;
        }
    }

    let bases: Vec<u64> = if n.bits() <= 64 {
        // The standard deterministic set for n < 2^64.
        vec![2, 325, 9375, 28178, 450775, 9780504, 1795265022]
    } else {
        // Probabilistic above 2^64; the caller must treat primality as
        // probable there, which is why `decompose` re-verifies c^2*s == n and
        // s squarefree rather than trusting this alone.
        vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    };

    let one = BigInt::one();
    let mut d = n - 1u32;
    let mut s: u32 = 0;
    while !(&d % BigInt::from(2u8)).is_zero() {
        d >>= 1u32;
        s += 1;
    }

    for a in bases {
        let base = BigInt::from(a);
        if base >= *n {
            continue;
        }
        let mut x = base.modpow(&d, n);
        if x == one || x == n - 1u32 {
            continue;
        }
        let mut witness = true;
        for _ in 0..s.saturating_sub(1) {
            x = (&x * &x) % n;
            if x == n - 1u32 {
                witness = false;
                break;
            }
        }
        if witness {
            return false;
        }
    }
    true
}

/// Exact squarefree test.
pub fn is_squarefree(n: &BigInt) -> Result<bool> {
    if n.is_negative() {
        return Err(SrsError::InvalidInput(format!(
            "squarefree test needs n >= 0, got {n}"
        )));
    }
    if *n < BigInt::from(4u8) {
        return Ok(true);
    }
    if is_multiple_of(n, &BigInt::from(4u8)) {
        return Ok(false);
    }
    let root = n.sqrt();
    if &root * &root == *n {
        return Ok(false);
    }

    // A small prime seen with exponent >= 2 during the scan settles it: p^2
    // divides n. Reading the cofactor alone is not enough, and was a real
    // bug -- is_squarefree(18) returned true because trial division reduces 18
    // to cofactor 1 while quietly absorbing 3^2 into the square part.
    let (_, _, remaining, repeated_small_prime) = trial_division_report(n)?;
    if repeated_small_prime {
        return Ok(false);
    }
    if remaining == BigInt::one() {
        return Ok(true);
    }
    is_squarefree_core(&remaining)
}

/// Squarefreeness of a cofactor with no small prime factor.
fn is_squarefree_core(n: &BigInt) -> Result<bool> {
    if *n == BigInt::one() {
        return Ok(true);
    }
    let r = n.sqrt();
    if &r * &r == *n {
        return Ok(false);
    }
    if is_probable_prime(n) {
        return Ok(true);
    }
    let d = pollard_rho(n)?;
    let other = n / &d;
    if d == other || gcd(&d, &other) != BigInt::one() {
        return Ok(false);
    }
    Ok(is_squarefree_core(&d)? && is_squarefree_core(&other)?)
}

/// Reference decomposition by exhaustive trial division, used only to
/// cross-check the real one.
///
/// The subtlety this gets right, and an earlier version got wrong: after
/// trial division stops, whatever cofactor is left belongs to the squarefree
/// part, not the square part. Reading n = 2 as (c=1, s=1) is the error --
/// 2 is itself squarefree, so s = 2.
///
/// One pass over odd divisors is enough for the range this is used on
/// (n < 3000): the loop condition is q*q <= n, and when it exits whatever is
/// left is 1 or a prime greater than q.
#[cfg(test)]
fn decompose_brute(n: &BigInt) -> (BigInt, BigInt) {
    let mut c = BigInt::one();
    let mut s = BigInt::one();
    let mut rest = n.clone();

    let mut q: u64 = 2;
    loop {
        let bq = BigInt::from(q);
        if &bq * &bq > rest {
            break;
        }
        let mut e = 0u32;
        while is_multiple_of(&rest, &bq) {
            rest /= &bq;
            e += 1;
        }
        if e > 0 {
            c *= bq.pow(e / 2);
            if e % 2 == 1 {
                s *= bq;
            }
        }
        q = if q == 2 { 3 } else { q + 2 };
    }

    // Whatever remains is 1 or a single prime, hence squarefree.
    if !rest.is_one() {
        s *= rest;
    }
    (c, s)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bi(n: u64) -> BigInt {
        BigInt::from(n)
    }

    #[test]
    fn test_known_decompositions() {
        let cases: [(u64, u64, u64); 10] = [
            (1, 1, 1),
            (2, 1, 2),
            (4, 2, 1),
            (12, 2, 3),
            (18, 3, 2),
            (20, 2, 5),
            (72, 6, 2),
            (100, 10, 1),
            (999983, 1, 999983),
            (200, 10, 2),
        ];
        for (n, want_c, want_s) in cases {
            let got = decompose(&bi(n)).unwrap();
            assert_eq!(
                (got.c, got.s),
                (bi(want_c), bi(want_s)),
                "n={n}"
            );
        }
    }

    #[test]
    fn test_matches_brute_force_exhaustively() {
        for n in 1u64..3000 {
            let got = decompose(&bi(n)).unwrap();
            let (want_c, want_s) = decompose_brute(&bi(n));
            assert_eq!(got.c, want_c, "c wrong for n={n}");
            assert_eq!(got.s, want_s, "s wrong for n={n}");
            assert_eq!(&got.c * &got.c * &got.s, bi(n), "identity broken for n={n}");
        }
    }

    #[test]
    fn test_perfect_square_fast_path_at_scale() {
        // The F1 radicand shape: p^2 with p of 1000 bits. Detected by one
        // integer square root, with no factoring attempted.
        let p: BigInt = BigInt::from(2u32).pow(999) + BigInt::one();
        let sq = &p * &p;
        let got = decompose(&sq).unwrap();
        assert_eq!(got.c, p);
        assert_eq!(got.s, BigInt::one());
    }

    #[test]
    fn test_f1_second_radicand_shape() {
        // alpha * q^2: only alpha needs factoring, q is a perfect-square
        // multiplier. This is the other F1 radicand.
        for (a, alpha_sf) in [(2u64, 6u64), (3, 3), (7, 14), (10, 110)] {
            let q = bi(1234567);
            let n = BigInt::from(alpha_sf) * &q * &q;
            let got = decompose(&n).unwrap();
            assert_eq!(got.s, bi(alpha_sf), "a={a}");
            assert_eq!(&got.c * &got.c * &got.s, n, "a={a}: identity");
        }
    }

    #[test]
    fn test_large_prime_cofactor() {
        // p*q with both primes above the trial-division bound needs rho.
        let p = bi(1000003);
        let q = bi(1000033);
        let n = &p * &q;
        let got = decompose(&n).unwrap();
        assert_eq!(got.c, BigInt::one());
        assert_eq!(&got.c * &got.c * &got.s, n);
    }

    #[test]
    fn test_large_prime_square_cofactor() {
        let p = bi(1000003);
        let n = &p * &p;
        let got = decompose(&n).unwrap();
        assert_eq!(got.c, p);
        assert_eq!(got.s, BigInt::one());
    }

    #[test]
    fn test_hard_cofactor_reports_the_input() {
        // Two 61-bit primes: rho cannot split this in reasonable time, and it
        // must say so rather than return a wrong factorisation.
        let p = BigInt::from(2305843009213693951u64);
        let q = BigInt::from(2305843009213693921u64);
        let n = &p * &q;
        let err = decompose(&n).unwrap_err();
        let message = err.to_string();
        assert!(
            message.contains("Pollard") || message.contains("cofactor"),
            "unexpected error: {message}"
        );
    }

    #[test]
    fn test_rejects_non_positive() {
        assert!(decompose(&BigInt::zero()).is_err());
        assert!(decompose(&BigInt::from(-4i8)).is_err());
    }

    #[test]
    fn test_is_squarefree_matches_definition() {
        for n in 0u64..4000 {
            let got = is_squarefree(&bi(n)).unwrap();
            let mut want = true;
            let mut p = 2u64;
            while p * p <= n {
                if n % (p * p) == 0 {
                    want = false;
                    break;
                }
                p += 1;
            }
            assert_eq!(got, want, "n={n}");
        }
    }

    #[test]
    fn test_is_squarefree_large() {
        let p = bi(1000003);
        assert!(is_squarefree(&p).unwrap());
        assert!(!is_squarefree(&(&p * &p)).unwrap());
        assert!(is_squarefree(&(&p * bi(1000033))).unwrap());
    }

    #[test]
    fn test_miller_rabin_rejects_known_composites() {
        for n in [341u64, 561, 1105, 1729, 2047, 2465, 2821, 6601] {
            assert!(!is_probable_prime(&bi(n)), "{n} is composite");
        }
    }

    #[test]
    fn test_miller_rabin_accepts_known_primes() {
        for p in [2u64, 3, 5, 97, 1000003, 1000033, 999983, 104729] {
            assert!(is_probable_prime(&bi(p)), "{p} is prime");
        }
    }

    #[test]
    fn test_rho_returns_a_real_divisor() {
        for n in [bi(1000003) * bi(1000033), bi(2) * bi(3) * bi(1000003)] {
            let d = pollard_rho(&n).unwrap();
            assert!(d > BigInt::one() && d < n, "rho returned {d} for {n}");
            assert!(is_multiple_of(&n, &d), "rho returned a non-divisor for {n}");
        }
    }
}
