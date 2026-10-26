//! Exact sign decision for square-root sums, integers only.
//!
//! Decides the sign of `S = SUM_i c_i * sqrt(a_i)` exactly. This is the piece
//! the Python reference implementation (`../../SRS/signexact.py`) has and this
//! crate previously did not.
//!
//! WHY IT IS EXACT
//! ---------------
//! The whole decision path is integer arithmetic. For `s >= 0` and `b >= 0`,
//! `t = isqrt(s << 2b)` certifies
//!
//! ```text
//!     t^2 <= s * 2^(2b) < (t+1)^2
//! ```
//!
//! and since `t >= 0`, dividing by `2^(2b)` and taking square roots gives
//! `t/2^b <= sqrt(s) <= (t+1)/2^b`, a bracket of width exactly `2^-b`.
//! Combining those term by term, with the endpoints swapped for negative
//! coefficients, yields integers `L, H` with `S` in `[L/2^b, H/2^b]`. So
//! `L > 0` proves `S > 0` and `H < 0` proves `S < 0`, both by integer
//! comparison. No float is involved, and none can be: `f64` has 53 bits of
//! mantissa and the F1 family needs separation on the order of `2^-N`.
//!
//! WHY IT TERMINATES
//! -----------------
//! `S = 0` if and only if every aggregated coefficient vanishes, because
//! `{sqrt(d) : d distinct squarefree}` is linearly independent over `Q`. The
//! zero test is therefore exact, and once `S != 0` the bracket width
//! `C * 2^-b` (with `C` the sum of absolute coefficients) drops below `|S|`
//! for all sufficiently large `b`. If the caller-set cap arrives first, this
//! returns an error rather than a guess.
//!
//! WHAT IT COSTS
//! -------------
//! The number of fractional bits needed. That is governed by how close `S`
//! sits to zero, which is the separation problem itself, and it is
//! exponential in the input size for adversarial inputs. The bit count is
//! reported in the result so callers see the cost rather than being surprised
//! by it.

use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};
use std::collections::BTreeMap;

use crate::error::{Result, SrsError};
use crate::squarefree::decompose;

/// Refuse to spin forever on adversarial inputs.
pub const DEFAULT_MAX_BITS: u64 = 1 << 22;

/// A sum to be decided: `coeff * sqrt(radicand)` pairs.
pub type Terms = Vec<(i64, BigInt)>;

/// Outcome of a sign decision.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SignResult {
    /// -1, 0 or +1.
    pub sign: i8,
    /// "zero-test", "rational" or "certified-interval".
    pub method: &'static str,
    /// Fractional bits the decision needed; 0 for the exact-algebraic methods.
    pub bits: u64,
    /// Number of distinct squarefree parts after aggregation.
    pub parts: usize,
}

/// The decision could not be made within the precision cap.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SignUnknown {
    pub bits: u64,
    /// log2 of the bracket width still being carried, as a bit count.
    pub width_log2: i64,
}

impl std::fmt::Display for SignUnknown {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "sign not decided within {} fractional bits (bracket still ~2^{} wide)",
            self.bits, self.width_log2
        )
    }
}

impl std::error::Error for SignUnknown {}

/// Decide the sign of `S = SUM (coeff, radicand)` exactly.
pub fn exact_sign(terms: &[(i64, BigInt)]) -> Result<SignResult> {
    exact_sign_within(terms, DEFAULT_MAX_BITS)
}

/// As `exact_sign`, with an explicit precision cap.
pub fn exact_sign_within(terms: &[(i64, BigInt)], max_bits: u64) -> Result<SignResult> {
    if max_bits == 0 {
        return Err(SrsError::InvalidInput("max_bits must be positive".into()));
    }
    let agg = aggregate(terms)?;
    decide_aggregated(&agg, max_bits)
}

/// Aggregate coefficients by squarefree part.
///
/// This is what makes the zero test exact: two terms contribute to the same
/// bucket precisely when they share a squarefree part, and distinct buckets
/// are linearly independent over `Q`.
pub fn aggregate(terms: &[(i64, BigInt)]) -> Result<BTreeMap<BigInt, BigInt>> {
    let mut agg: BTreeMap<BigInt, BigInt> = BTreeMap::new();
    for (coeff, radicand) in terms {
        if *radicand <= BigInt::zero() {
            return Err(SrsError::InvalidInput(format!(
                "radicand must be positive, got {radicand}"
            )));
        }
        if *coeff == 0 {
            continue;
        }
        let d = decompose(radicand)?;
        let contribution = BigInt::from(*coeff) * d.c;
        *agg.entry(d.s).or_insert_with(BigInt::zero) += contribution;
    }
    // Drop buckets that cancelled exactly.
    agg.retain(|_, v| !v.is_zero());
    Ok(agg)
}

/// Certified integer bracket for `sqrt(s)`, scaled by `2^bits`.
///
/// Returns `(lo_scaled, hi_scaled)` with
/// `lo/2^bits <= sqrt(s) <= hi/2^bits` and `hi - lo == 1`.
pub fn sqrt_bracket_int(s: &BigInt, bits: u64) -> Result<(BigInt, BigInt)> {
    if s.is_negative() {
        return Err(SrsError::InvalidInput(format!(
            "sqrt bracket needs s >= 0, got {s}"
        )));
    }
    if *s == BigInt::zero() {
        return Ok((BigInt::zero(), BigInt::zero()));
    }
    let t = (s << (2 * bits)).sqrt();
    Ok((t.clone(), t + BigInt::one()))
}

/// Certified bracket for an aggregated sum, scaled by `2^bits`.
///
/// Endpoint selection is sign-aware: a positive coefficient maps `[lo, hi]`
/// to `[c*lo, c*hi]`, a negative one to `[c*hi, c*lo]`. Summing lower against
/// upper endpoints keeps the enclosure valid term by term, with no assumption
/// about cancellation.
pub fn sum_bracket_int(agg: &BTreeMap<BigInt, BigInt>, bits: u64) -> Result<(BigInt, BigInt)> {
    let mut lo = BigInt::zero();
    let mut hi = BigInt::zero();
    for (s, c) in agg {
        if c.is_zero() {
            continue;
        }
        if *s == BigInt::one() {
            let exact = c << bits;
            lo += &exact;
            hi += exact;
            continue;
        }
        let (t, t1) = sqrt_bracket_int(s, bits)?;
        if c.is_positive() {
            lo += c * t;
            hi += c * t1;
        } else {
            lo += c * t1;
            hi += c * t;
        }
    }
    Ok((lo, hi))
}

/// Decide the sign of an already-aggregated sum.
pub fn decide_aggregated(agg: &BTreeMap<BigInt, BigInt>, max_bits: u64) -> Result<SignResult> {
    if agg.is_empty() {
        return Ok(SignResult { sign: 0, method: "zero-test", bits: 0, parts: 0 });
    }

    // A sum with no irrational part is a plain integer comparison.
    if agg.len() == 1 && agg.contains_key(&BigInt::one()) {
        let total = agg[&BigInt::one()].clone();
        let sign = if total.is_positive() {
            1
        } else if total.is_negative() {
            -1
        } else {
            0
        };
        return Ok(SignResult { sign, method: "rational", bits: 0, parts: 1 });
    }

    let parts = agg.len();
    let c_sum: BigInt = agg.values().fold(BigInt::zero(), |a, b| a + b.abs());
    let start = start_bits(agg);

    let mut bits = start;
    loop {
        if bits > max_bits {
            return Err(SrsError::ArithmeticError(SignUnknown {
                bits: max_bits,
                width_log2: c_sum.bits() as i64 - bits as i64,
            }
            .to_string()));
        }
        let (lo, hi) = sum_bracket_int(agg, bits)?;
        if lo.is_positive() {
            return Ok(SignResult { sign: 1, method: "certified-interval", bits, parts });
        }
        if hi.is_negative() {
            return Ok(SignResult { sign: -1, method: "certified-interval", bits, parts });
        }
        bits *= 2;
    }
}

/// First precision to try.
///
/// The bracket width is `C * 2^-bits` with `C` the sum of absolute
/// coefficients, so log2(C) plus the scale of the radicands is the natural
/// starting point. This avoids wasted halvings without assuming anything
/// about how close `S` is to zero; the doubling loop still does the real work.
fn start_bits(agg: &BTreeMap<BigInt, BigInt>) -> u64 {
    let c_sum: BigInt = agg.values().fold(BigInt::zero(), |a, b| a + b.abs());
    let c_bits = c_sum.bits();
    let max_radicand_bits = agg.keys().map(|s| s.bits()).max().unwrap_or(0);
    16 + c_bits + 2 * max_radicand_bits
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::f1;
    use num_traits::ToPrimitive;

    fn bi(n: u64) -> BigInt {
        BigInt::from(n)
    }

    fn signs(terms: &[(i64, BigInt)]) -> i8 {
        exact_sign(terms).unwrap().sign
    }

    #[test]
    fn test_empty_sum_is_zero() {
        let r = exact_sign(&[]).unwrap();
        assert_eq!(r.sign, 0);
        assert_eq!(r.method, "zero-test");
    }

    #[test]
    fn test_all_zero_coefficients_is_zero() {
        let r = exact_sign(&[(0, bi(5)), (0, bi(7))]).unwrap();
        assert_eq!(r.sign, 0);
    }

    #[test]
    fn test_exact_zero_by_cancellation() {
        // sqrt(12) - 2 sqrt(3) = 2 sqrt(3) - 2 sqrt(3) = 0
        let r = exact_sign(&[(1, bi(12)), (-2, bi(3))]).unwrap();
        assert_eq!(r.sign, 0);
        assert_eq!(r.method, "zero-test");
    }

    #[test]
    fn test_rational_sum() {
        // 2 - 3 = -1
        let r = exact_sign(&[(1, bi(4)), (-1, bi(9))]).unwrap();
        assert_eq!(r.sign, -1);
        assert_eq!(r.method, "rational");
    }

    #[test]
    fn test_known_signs() {
        // sqrt(2) - sqrt(3) < 0
        assert_eq!(signs(&[(1, bi(2)), (-1, bi(3))]), -1);
        // sqrt(2) + sqrt(3) - sqrt(5) > 0
        assert_eq!(signs(&[(1, bi(2)), (1, bi(3)), (-1, bi(5))]), 1);
        // 2 sqrt(2) - 3 sqrt(3) < 0
        assert_eq!(signs(&[(2, bi(2)), (-3, bi(3))]), -1);
        // sqrt(12) - sqrt(3) = sqrt(3) > 0
        assert_eq!(signs(&[(1, bi(12)), (-1, bi(3))]), 1);
    }

    #[test]
    fn test_sqrt_bracket_certifies_by_squaring() {
        // The defining property, checked exactly: lo^2 <= s <= hi^2 after
        // rescaling. No float anywhere.
        for s in [2u64, 3, 5, 12, 17, 999, 1000003] {
            for bits in [1u64, 4, 17, 64] {
                let (lo, hi) = sqrt_bracket_int(&bi(s), bits).unwrap();
                let scale = BigInt::one() << bits;
                let s_scaled = bi(s) * &scale * &scale;
                assert!(lo.is_positive(), "lo should be positive for s={s}");
                assert!(&lo * &lo <= s_scaled, "lo^2 > s*2^2b for s={s}, bits={bits}");
                assert!(s_scaled <= &hi * &hi, "s*2^2b > hi^2 for s={s}, bits={bits}");
                assert_eq!(hi - lo, BigInt::one(), "width should be 1 ulp");
            }
        }
    }

    #[test]
    fn test_sqrt_bracket_edge_cases() {
        assert_eq!(sqrt_bracket_int(&BigInt::zero(), 8).unwrap(), (BigInt::zero(), BigInt::zero()));
        let (lo, hi) = sqrt_bracket_int(&BigInt::one(), 8).unwrap();
        let scale = BigInt::one() << 8;
        assert_eq!(lo, scale);
        assert_eq!(hi, scale + BigInt::one());
        assert!(sqrt_bracket_int(&BigInt::from(-1i8), 8).is_err());
    }

    #[test]
    fn test_deep_cancellation_is_exact() {
        // sqrt(N^2-1) - N is negative but tiny; a float surrogate would
        // return 0 or the wrong sign here.
        // N large enough that sqrt(N^2-1) and N are indistinguishable in f64.
        let big = BigInt::from(1_000_000_007u64);
        let nn = &big * &big;
        let r = exact_sign(&[(1, &nn - 1u32), (-1, nn.clone())]).unwrap();
        assert_eq!(r.sign, -1);
        assert!(r.bits > 32, "should have needed many bits, used {}", r.bits);
    }

    #[test]
    fn test_aggregate_merges_by_squarefree_part() {
        let agg = aggregate(&[(1, bi(12)), (2, bi(3))]).unwrap();
        // sqrt(12) + 2 sqrt(3) = 2 sqrt(3) + 2 sqrt(3) = 4 sqrt(3)
        assert_eq!(agg.len(), 1);
        assert_eq!(agg[&bi(3)], bi(4));
    }

    #[test]
    fn test_aggregate_drops_cancelled_buckets() {
        let agg = aggregate(&[(1, bi(12)), (-2, bi(3))]).unwrap();
        assert!(agg.is_empty());
    }

    #[test]
    fn test_rejects_bad_terms() {
        assert!(exact_sign(&[(1, BigInt::zero())]).is_err());
        assert!(exact_sign(&[(1, BigInt::from(-4i8))]).is_err());
        assert!(exact_sign_within(&[(1, bi(2))], 0).is_err());
    }

    #[test]
    fn test_cap_raises_rather_than_guessing() {
        // A cap far too small to decide must error, not return a sign.
        let big = BigInt::from(1_000_000_007u64);
        let nn = &big * &big;
        let err = exact_sign_within(&[(1, &nn - 1u32), (-1, nn.clone())], 8).unwrap_err();
        assert!(
            err.to_string().contains("not decided"),
            "unexpected error: {err}"
        );
    }

    /// THE CROSS-CHECK THAT MATTERS: the certified decider must agree with
    /// the closed-form sign of the F1 family, which comes from the exact
    /// Pell residue and needs no computation to know.
    #[test]
    fn test_f1_family_signs_match_the_pell_residue() {
        for a in 2..=8u64 {
            for n in 1..=25usize {
                let cv = f1::pell_convergents(a, n).unwrap();
                let (p, q) = &cv[n];
                let alpha = BigInt::from(a * (a + 1));
                let terms = vec![(1i64, p * p), (-1i64, alpha.clone() * q * q)];
                let got = exact_sign(&terms).unwrap();
                assert_eq!(
                    i32::from(got.sign),
                    f1::expected_sign(n),
                    "a={a}, n={n}: decider disagrees with the closed form"
                );
            }
        }
    }

    /// The bit cost must scale with the input, which is the separation law
    /// showing up in the implementation rather than only on paper.
    #[test]
    fn test_bit_cost_grows_with_instance_size() {
        let mut previous = 0u64;
        for n in [1usize, 5, 20, 60] {
            let cv = f1::pell_convergents(2, n).unwrap();
            let (p, q) = &cv[n];
            let alpha = BigInt::from(6u8);
            let r = exact_sign(&[(1i64, p * p), (-1i64, alpha.clone() * q * q)]).unwrap();
            assert!(
                r.bits > previous,
                "n={n}: bits did not grow ({} then {})",
                previous,
                r.bits
            );
            previous = r.bits;
        }
    }

    #[test]
    fn test_sign_is_always_in_range() {
        // Property: whatever the input, the reported sign is -1, 0 or +1.
        let mut state: u64 = 0x243F_6A88_85A3_08D3;
        let mut next = || {
            // xorshift64, so the cases are varied but reproducible
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for _ in 0..200 {
            let k = (next() % 3) as usize + 1;
            let terms: Terms = (0..k)
                .map(|_| {
                    let coeff = if next() % 2 == 0 { 1i64 } else { -1i64 };
                    let rad = BigInt::from(next() % 200 + 2);
                    (coeff, rad)
                })
                .collect();
            let sign = exact_sign(&terms).unwrap().sign;
            assert!((-1..=1).contains(&sign), "sign {sign} out of range");
        }
    }

    #[test]
    fn test_bracket_contains_the_true_value() {
        // Independent check of the enclosure: the bracket's endpoints must
        // bracket the sum to within the reported width, verified by
        // comparing the |S| estimate against a high-precision float at a size
        // where f64 still has plenty of bits. This is an oracle, not a
        // decision path.
        for (coeffs, rads) in [
            (vec![1i64, -1i64], vec![2u64, 3]),
            (vec![1, 1, -1], vec![2, 3, 5]),
            (vec![3, -1], vec![7, 11]),
        ] {
            let terms: Terms = coeffs
                .iter()
                .zip(rads.iter())
                .map(|(c, r)| (*c, bi(*r)))
                .collect();
            let agg = aggregate(&terms).unwrap();
            let bits = 64u64;
            let (lo, hi) = sum_bracket_int(&agg, bits).unwrap();
            let scale = 2f64.powi(bits as i32);
            let approx: f64 = rads
                .iter()
                .zip(coeffs.iter())
                .map(|(r, c)| *c as f64 * (*r as f64).sqrt())
                .sum();
            let lo_f = lo.to_f64().unwrap() / scale;
            let hi_f = hi.to_f64().unwrap() / scale;
            assert!(
                lo_f <= approx + 1e-9 && approx <= hi_f + 1e-9,
                "bracket [{lo_f}, {hi_f}] does not contain {approx}"
            );
        }
    }
}
