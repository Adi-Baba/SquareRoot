//! F1: the Pell family and its separation law, in exact integer arithmetic.
//!
//! Fix a >= 2 and set alpha = a(a+1).  The continued fraction of sqrt(alpha) is
//! periodic from the start:
//!
//! ```text
//! sqrt(a(a+1)) = [a; 2, 2a, 2, 2a, 2, 2a, ...]
//! ```
//!
//! and the SRS instance built from its convergents is
//!
//! ```text
//! S_n = sqrt(p_n^2) - sqrt(alpha q_n^2) = p_n - q_n sqrt(alpha)
//! ```
//!
//! The exact residue is
//!
//! ```text
//! p_n^2 - alpha q_n^2 = +1  for odd n
//! p_n^2 - alpha q_n^2 = -a  for even n
//! ```
//!
//! Both are exact integer identities, so the sign of S_n is fixed by parity
//! alone: +1 on odd n, -1 on even n.  A past revision's docstring said the
//! residue was "+-1" with even n "varying"; it does not, and the sign is
//! derivable rather than computed.
//!
//! The separation law is
//!
//! ```text
//! log2(1/|S_n|) = N/4 + 1 - log2(|D|) + o(1),   D = p_n^2 - alpha q_n^2
//! ```
//!
//! with N the bit size of the two radicands, so log2(1/|S_n|)/N -> 1/4.  The
//! float helpers below MEASURE that; the identities are checked with integers.

use num_bigint::BigInt;
use num_traits::{One, Signed, ToPrimitive};

use crate::error::{Result, SrsError};

/// Largest `a` accepted, so that `a*(a+1)` cannot overflow a `u64`.
const MAX_A: u64 = (1u64 << 31) - 2;

/// CF convergents `(p_n, q_n)` of `sqrt(a(a+1))` for `n = 0..=n_max`.
///
/// `a_1 = 2`, `a_2 = 2a`, `a_3 = 2`, ... per the period above; index 0 is the
/// integer part `a/1`.
pub fn pell_convergents(a: u64, n_max: usize) -> Result<Vec<(BigInt, BigInt)>> {
    check_a(a)?;
    let mut p = vec![BigInt::from(a), BigInt::from(2 * a + 1)];
    let mut q = vec![BigInt::one(), BigInt::from(2u8)];

    for i in 2..=n_max {
        let a_i = if i % 2 == 1 { 2 } else { 2 * a };
        let a_i_big = BigInt::from(a_i);
        p.push(&a_i_big * &p[i - 1] + &p[i - 2]);
        q.push(&a_i_big * &q[i - 1] + &q[i - 2]);
    }

    Ok(p.into_iter().zip(q).collect())
}

fn check_a(a: u64) -> Result<()> {
    if a < 2 {
        return Err(SrsError::InvalidInput(format!("a must be >= 2, got {a}")));
    }
    if a > MAX_A {
        return Err(SrsError::InvalidInput(format!(
            "a={a} exceeds {MAX_A}; a*(a+1) must fit in a u64"
        )));
    }
    Ok(())
}

/// The exact residue `p_n^2 - alpha q_n^2`: +1 on odd n, `-a` on even n.
pub fn pell_residue(a: u64, n: usize) -> Result<BigInt> {
    check_a(a)?;
    Ok(if n % 2 == 1 {
        BigInt::one()
    } else {
        -BigInt::from(a)
    })
}

/// The sign of `S_n`, fixed by parity: +1 on odd n, -1 on even n.
pub fn expected_sign(n: usize) -> i32 {
    if n % 2 == 1 {
        1
    } else {
        -1
    }
}

/// The F1 instance as `(coefficient, radicand)` pairs:
/// `S_n = sqrt(p_n^2) - sqrt(alpha q_n^2)`.
pub fn generate_f1(a: u64, n: usize) -> Result<Vec<(i64, BigInt)>> {
    let cv = pell_convergents(a, n)?;
    let (p, q) = cv
        .get(n)
        .ok_or_else(|| SrsError::InvalidInput(format!("n={n} is out of range")))?;
    let alpha = BigInt::from(a * (a + 1));
    Ok(vec![(1, p * p), (-1, alpha * q * q)])
}

/// Verify `p_n^2 - alpha q_n^2 = +1` for every odd `n` in `1..=n_max`.
///
/// Exact integer arithmetic; no floating point is involved.
pub fn verify_pell_identity(a: u64, n_max: usize) -> Result<()> {
    let cv = pell_convergents(a, n_max)?;
    let alpha = BigInt::from(a * (a + 1));

    for (n, (p, q)) in cv.iter().enumerate() {
        if n % 2 == 0 || n == 0 {
            continue;
        }
        let diff = p * p - &alpha * q * q;
        if diff != BigInt::one() {
            return Err(SrsError::ArithmeticError(format!(
                "a={a}, n={n}: p^2-alpha*q^2={diff}, expected 1"
            )));
        }
    }
    Ok(())
}

/// Verify the full residue law, both parities, in exact integers.
pub fn verify_pell_residues(a: u64, n_max: usize) -> Result<()> {
    let cv = pell_convergents(a, n_max)?;
    let alpha = BigInt::from(a * (a + 1));
    for (n, (p, q)) in cv.iter().enumerate() {
        let diff = p * p - &alpha * q * q;
        let want = pell_residue(a, n)?;
        if diff != want {
            return Err(SrsError::ArithmeticError(format!(
                "a={a}, n={n}: p^2-alpha*q^2={diff}, expected {want}"
            )));
        }
    }
    Ok(())
}

/// Input size `N` in bits: `bit_len(p^2) + bit_len(alpha q^2)`, exact.
pub fn input_bits(a: u64, n: usize) -> Result<u64> {
    let cv = pell_convergents(a, n)?;
    let (p, q) = cv
        .get(n)
        .ok_or_else(|| SrsError::InvalidInput(format!("n={n} is out of range")))?;
    let alpha = BigInt::from(a * (a + 1));
    Ok((p * p).bits() + (alpha * q * q).bits())
}

/// Largest integer an `f64` represents exactly.
const F64_EXACT_INT: f64 = 9_007_199_254_740_992.0; // 2^53

/// MEASUREMENT ONLY.  `(n, log2(1/|S_n|)/N, log2(1/|S_n|))` for `n = 1..=n_max`.
///
/// Float, and bounded on purpose: past 2^53 an `f64` no longer represents
/// `p_n` or `q_n` exactly, so the reported ratio is not a rounded version of
/// the truth but a different number.  This errors there rather than printing
/// it.  Use `input_bits` and the exact residue for anything that has to be
/// right, and expect this to stop around n = 17.
pub fn log2_ratio_law(a: u64, n_max: usize) -> Result<Vec<(usize, f64, f64)>> {
    check_a(a)?;
    let cv = pell_convergents(a, n_max)?;
    let alpha_int = BigInt::from(a * (a + 1));
    let alpha = a as f64 * (a + 1) as f64;
    let sqrt_alpha = alpha.sqrt();

    let mut results = Vec::with_capacity(n_max);
    for (n, (p, q)) in cv.iter().enumerate().skip(1) {
        let (p_f, q_f) = match (p.to_f64(), q.to_f64()) {
            (Some(pf), Some(qf)) if pf < F64_EXACT_INT && qf < F64_EXACT_INT => (pf, qf),
            _ => {
                return Err(SrsError::ArithmeticError(format!(
                    "a={a}, n={n}: p_n or q_n exceeds 2^53, so f64 no longer \
                     represents them exactly and the ratio would not be a \
                     rounded truth. Use input_bits and the exact residue."
                )))
            }
        };
        // |S_n| = |D| / (p + q sqrt(alpha)), so log2(1/|S_n|) is the log of
        // (p + q sqrt(alpha)) divided by |D|.  The earlier revision had this
        // sign reversed, which made every reported separation negative.
        let d = (p * p - &alpha_int * q * q)
            .abs()
            .to_f64()
            .unwrap_or(f64::INFINITY);
        let log2_inv_s = (p_f + q_f * sqrt_alpha).log2() - d.log2();
        let n_bits = 2.0 * p_f.log2() + 2.0 * q_f.log2() + alpha.log2();
        let ratio = if n_bits > 0.0 { log2_inv_s / n_bits } else { f64::INFINITY };
        results.push((n, ratio, log2_inv_s));
    }
    Ok(results)
}

/// `1 - log2(|D|)`: the exact constant in `log2(1/|S_n|) = N/4 + that`.
pub fn separation_constant(a: u64, n: usize) -> Result<f64> {
    let d = pell_residue(a, n)?
        .abs()
        .to_f64()
        .ok_or_else(|| SrsError::ArithmeticError("residue does not fit in f64".into()))?;
    Ok(1.0 - d.log2())
}

/// log2(1/|S_n|) for an F1 instance, from the closed form.
///
/// Measured, because N and |S| are only available in floating point. The
/// identity it checks is
///
/// ```text
///     log2(1/|S_n|) = N/4 + 1 - log2(|D|) + o(1),   D = p_n^2 - alpha q_n^2
/// ```
///
/// so the residual `log2(1/|S_n|) - N/4` should equal `1 - log2|D|`. The
/// o(1) comes from `p_n + q_n sqrt(alpha) = 2 p_n + o(1)` and is largest at
/// the smallest n, decaying from there; at these sizes it is under 1e-4.
pub fn separation_residual(a: u64, n: usize) -> Result<f64> {
    let cv = pell_convergents(a, n)?;
    let (p, q) = cv
        .get(n)
        .ok_or_else(|| SrsError::InvalidInput(format!("n={n} is out of range")))?;
    let alpha = BigInt::from(a * (a + 1));
    // Same 2^53 guard as log2_ratio_law, and for the same reason: past that
    // point an f64 no longer represents p_n or q_n exactly, so the residual
    // would not be a rounded version of the truth but a different number.
    let p_f = p.to_f64().ok_or_else(|| {
        SrsError::ArithmeticError(format!("p_n exceeds f64 range at n={n}"))
    })?;
    let q_f = q.to_f64().ok_or_else(|| {
        SrsError::ArithmeticError(format!("q_n exceeds f64 range at n={n}"))
    })?;
    if p_f >= F64_EXACT_INT || q_f >= F64_EXACT_INT {
        return Err(SrsError::ArithmeticError(format!(
            "a={a}, n={n}: p_n or q_n exceeds 2^53, so f64 no longer represents \
             them exactly and the residual would not be a rounded truth. Use \
             input_bits and the exact residue instead."
        )));
    }
    let alpha_f = (a * (a + 1)) as f64;
    let d = (p * p - &alpha * q * q).abs().to_f64().ok_or_else(|| {
        SrsError::ArithmeticError("residue exceeds f64 range".into())
    })?;

    let log2_inv_s = (p_f + q_f * alpha_f.sqrt()).log2() - d.log2();
    let n_bits = 2.0 * p_f.log2() + 2.0 * q_f.log2() + alpha_f.log2();
    Ok(log2_inv_s - n_bits / 4.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pell_convergents() {
        let cv = pell_convergents(3, 6).unwrap();
        assert_eq!(cv[0], (BigInt::from(3u8), BigInt::from(1u8)));
        assert_eq!(cv[1], (BigInt::from(7u8), BigInt::from(2u8)));
        assert_eq!(cv[2], (BigInt::from(45u8), BigInt::from(13u8)));
        assert_eq!(cv[5], (BigInt::from(1351u32), BigInt::from(390u32)));
    }

    #[test]
    fn test_pell_identity() {
        verify_pell_identity(3, 20).unwrap();
    }

    #[test]
    fn test_residue_law_both_parities() {
        for a in 2..20u64 {
            verify_pell_residues(a, 20).unwrap();
        }
    }

    #[test]
    fn test_residue_values() {
        // a=3, n=2: 45^2 - 12*13^2 = 2025 - 2028 = -3
        let cv = pell_convergents(3, 2).unwrap();
        let (p, q) = &cv[2];
        let diff = p * p - BigInt::from(12u8) * q * q;
        assert_eq!(diff, BigInt::from(-3i8));
        assert_eq!(pell_residue(3, 2).unwrap(), BigInt::from(-3i8));
    }

    #[test]
    fn test_expected_sign_follows_parity() {
        assert_eq!(expected_sign(0), -1);
        assert_eq!(expected_sign(1), 1);
        assert_eq!(expected_sign(2), -1);
        assert_eq!(expected_sign(7), 1);
    }

    #[test]
    fn test_rejects_bad_a() {
        for bad in [0u64, 1] {
            assert!(pell_convergents(bad, 3).is_err());
            assert!(pell_residue(bad, 1).is_err());
        }
        assert!(pell_convergents(MAX_A + 1, 1).is_err());
    }

    #[test]
    fn test_generate_f1_matches_construction() {
        let terms = generate_f1(3, 5).unwrap();
        let cv = pell_convergents(3, 5).unwrap();
        let (p, q) = &cv[5];
        let alpha = BigInt::from(12u8);
        assert_eq!(terms, vec![(1, p * p), (-1, alpha * q * q)]);
    }

    #[test]
    fn test_input_bits_is_exact_and_large() {
        // Past n=60 the radicands exceed f64 entirely, so a float-based
        // N would have overflowed or lost every low bit.
        let bits = input_bits(2, 200).unwrap();
        assert!(bits > 1000, "expected a large bit count, got {bits}");
    }

    #[test]
    fn test_log2_ratio_law_errors_instead_of_reporting_nonsense() {
        // Past n~17 the denominators stop being exact in f64, so the helper
        // must refuse rather than print a plausible wrong ratio.
        assert!(log2_ratio_law(2, 40).is_err());
        // In range it works and the separation is positive.
        let rows = log2_ratio_law(2, 10).unwrap();
        assert_eq!(rows.len(), 10);
        for (_, ratio, log2_inv_s) in rows {
            assert!(log2_inv_s > 0.0, "separation must be positive");
            assert!(ratio > 0.0, "ratio must be positive");
        }
    }

    #[test]
    fn test_log2_ratio_law_tends_to_one_quarter() {
        // The proved law: log2(1/|S_n|)/N -> 1/4, approached from above on
        // odd n because the correction there is exactly +1.
        for a in [2u64, 3, 4] {
            for (n, ratio, _) in log2_ratio_law(a, 12).unwrap() {
                if n >= 5 && n % 2 == 1 {
                    assert!(ratio > 0.25, "a={a}, n={n}: ratio={ratio}");
                    assert!(ratio < 0.25 + 0.05, "a={a}, n={n}: ratio={ratio}");
                }
            }
        }
    }

    /// E4 from TODO.md: the symbolic limit of the separation law, pinned.
    ///
    /// `log2(1/|S_n|) - N/4` must equal `1 - log2|D|`, where |D| is 1 on odd
    /// n and a on even n. This turns the documented limit into something the
    /// build fails on if it regresses, rather than a comment.
    #[test]
    fn test_separation_law_symbolic_limit() {
        for a in [2u64, 3, 4, 5, 6] {
            for n in 1..=20usize {
                let want = separation_constant(a, n).unwrap();
                let got = separation_residual(a, n).unwrap();
                assert!(
                    (got - want).abs() < 1e-3,
                    "a={a}, n={n}: residual {got:.6}, want 1-log2|D| = {want:.6}"
                );
            }
        }
    }

    #[test]
    fn test_separation_residual_refuses_beyond_f64() {
        // Past 2^53 the float ratio is not a rounded truth, so it must error
        // rather than print a plausible wrong number.
        let err = separation_residual(2, 40).unwrap_err();
        assert!(err.to_string().contains("2^53"), "unexpected error: {err}");
    }

    #[test]
    fn test_separation_constant() {
        assert_eq!(separation_constant(3, 5).unwrap(), 1.0);   // |D| = 1
        assert!((separation_constant(3, 6).unwrap() - (1.0 - 3f64.log2())).abs() < 1e-12);
        assert!((separation_constant(4, 6).unwrap() - (-1.0)).abs() < 1e-12); // 1 - 2
    }
}
