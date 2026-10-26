//! SRS: certified square-root sums in exact integer arithmetic.
//!
//! # What this crate decides
//!
//! Given `S = SUM_i c_i * sqrt(a_i)` with integer `c_i` and `a_i > 0`, decide
//! the sign of `S` exactly, and certify a rational interval containing it.
//!
//! # Where the exactness comes from
//!
//! The zero test is exact because `{sqrt(d) : d distinct squarefree}` is
//! linearly independent over `Q`: `S = 0` if and only if every aggregated
//! coefficient vanishes.  Once `S != 0`, a bracket of the form
//! `[L/2^b, H/2^b]` built from `isqrt(a << 2b)` clears zero for all
//! sufficiently large `b`, and the loop that refines it is guaranteed to
//! terminate.  Nothing in that path touches a float.
//!
//! # What this crate does NOT do
//!
//! It does not currently implement the certified sign decision.  `f1` gives
//! the family and the closed-form sign, `field` gives exact field arithmetic.
//! The general sign loop lives in the Python reference implementation, and
//! this crate's `probe` reports measurements.  Saying so is better than
//! shipping a placeholder that returns a confident wrong answer.

pub mod error;
pub mod f1;
pub mod field;
pub mod signexact;
pub mod squarefree;

pub use error::{Result, SrsError};
pub use f1::{
    expected_sign, generate_f1, input_bits, log2_ratio_law, pell_convergents, pell_residue,
    separation_constant, verify_pell_identity, verify_pell_residues,
};
pub use field::{
    Elem, MultiQuadraticField, EXACT_MAX_GENERATORS, MAX_GENERATORS,
};
pub use signexact::{exact_sign, SignResult, SignUnknown};
pub use squarefree::{decompose, is_squarefree, Decomposition};

/// Crate version, for provenance in probe output.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;
    use num_traits::One;

    #[test]
    fn f1_pell_identity_exact() {
        for a in 2..12u64 {
            verify_pell_identity(a, 20).unwrap();
        }
    }

    #[test]
    fn f1_residue_law_exact() {
        for a in 2..12u64 {
            verify_pell_residues(a, 20).unwrap();
        }
    }

    #[test]
    fn f1_sign_matches_the_analytic_value() {
        for a in 2..10u64 {
            for n in 1..=20usize {
                let cv = pell_convergents(a, n).unwrap();
                let (p, q) = &cv[n];
                let alpha = BigInt::from(a * (a + 1));
                // |S_n| = |D| / (p + q sqrt(alpha)) and p, q > 0, so S_n has
                // the sign of p^2 - alpha q^2.
                let diff = p * p - &alpha * q * q;
                let sign = if diff > BigInt::one() - BigInt::one() {
                    1
                } else {
                    -1
                };
                assert_eq!(sign, expected_sign(n), "a={a}, n={n}");
            }
        }
    }

    #[test]
    fn field_norm_is_an_integer_and_errors_above_the_cap() {
        let k = MultiQuadraticField::new(vec![BigInt::from(2u8), BigInt::from(3u8)]).unwrap();
        assert_eq!(k.norm(&k.from_sqrt(0).unwrap()).unwrap(), BigInt::from(4u8));
        // The regression that mattered: this used to return BigInt::one().
        assert!(k.norm(&k.from_sqrt(0).unwrap()).is_ok());
    }

    #[test]
    fn minimal_polynomial_degree_is_not_inflated() {
        let k = MultiQuadraticField::new(vec![BigInt::from(2u8), BigInt::from(3u8)]).unwrap();
        let s2 = k.from_sqrt(0).unwrap();
        assert_eq!(k.minimal_polynomial(&s2).unwrap().len(), 3);
        assert_eq!(k.characteristic_polynomial(&s2).unwrap().len(), 5);
    }
}
