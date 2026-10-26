//! Multiquadratic field arithmetic over the integers, exact.
//!
//! K = Q(sqrt(s_1), ..., sqrt(s_m)) with the s_i distinct, squarefree and
//! independent.  Elements are sparse maps from a bitmask to a BigInt
//! coefficient; bit i of a mask is set when the basis element
//! sqrt(prod of s_j for j in mask) appears.
//!
//! Two limits are enforced rather than left to run out of memory:
//!
//! * `MAX_GENERATORS` caps m.  The basis is a table of 2^m entries, so m = 40
//!   would ask for 2^40 BigInts before any work is done.  A past revision
//!   allocated that table with no check at all.
//! * `EXACT_MAX_GENERATORS` caps the operations that are genuinely
//!   exponential in m (`norm`, `trace`, `minimal_polynomial`,
//!   `characteristic_polynomial`).  Past revisions returned a placeholder
//!   value from `norm` above m = 6, which is worse than an error: a caller
//!   that routes a sign decision through it gets a confident wrong answer.
//!   Every one of them now returns `Result` and errors instead.

use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};
use std::collections::{HashMap, HashSet};

use crate::error::{Result, SrsError};

/// Hard cap on generators: the basis table has 2^m entries.
pub const MAX_GENERATORS: usize = 20;

/// Cap on the operations that cost 2^m field multiplications.
pub const EXACT_MAX_GENERATORS: usize = 6;

/// A field element: bitmask -> coefficient, with zero coefficients removed.
pub type Elem = HashMap<usize, BigInt>;

pub struct MultiQuadraticField {
    /// The squarefree generators, sorted and deduplicated.
    pub s: Vec<BigInt>,
    /// Number of generators.
    pub m: usize,
    /// `basis_radicand[mask]` is the product of `s[i]` over the bits of `mask`.
    pub basis_radicand: Vec<BigInt>,
}

impl MultiQuadraticField {
    /// Build Q(sqrt(s_1), ..., sqrt(s_m)) from squarefree generators.
    ///
    /// Errors if the list is empty, a generator is below 2, a generator is not
    /// squarefree, the generators are not independent, or `m` exceeds
    /// `MAX_GENERATORS`.
    pub fn new(squarefree_parts: Vec<BigInt>) -> Result<Self> {
        if squarefree_parts.is_empty() {
            return Err("at least one squarefree part required".into());
        }

        let mut s = squarefree_parts;
        s.sort();
        s.dedup();

        for v in &s {
            if *v < BigInt::from(2u8) {
                return Err(SrsError::InvalidInput(format!(
                    "squarefree parts must be >= 2, got {v}"
                )));
            }
            if !Self::is_squarefree(v)? {
                return Err(SrsError::InvalidInput(format!("{v} is not squarefree")));
            }
        }

        let m = s.len();
        if m > MAX_GENERATORS {
            return Err(SrsError::InvalidInput(format!(
                "too many generators: m={m} exceeds MAX_GENERATORS={MAX_GENERATORS} \
                 (the basis alone would need 2^{m} entries)"
            )));
        }

        Self::verify_independence(&s)?;

        let size = 1usize << m;
        let mut basis_radicand = vec![BigInt::one(); size];
        for mask in 1..size {
            let lsb = mask & mask.wrapping_neg();
            let i = lsb.trailing_zeros() as usize;
            let prev = mask ^ lsb;
            basis_radicand[mask] = &basis_radicand[prev] * &s[i];
        }

        Ok(Self { s, m, basis_radicand })
    }

    /// Exact squarefreeness test by trial division.
    ///
    /// O(sqrt(n)) integer work.  The generator magnitudes in this crate come
    /// from squarefree parts, which are small by construction, so a plain
    /// trial division is the honest tool.  `TRIAL_DIVISION_LIMIT` bounds it:
    /// past that the loop gives up and reports failure rather than grinding,
    /// so an oversized input is a clear error and not a hang.
    fn is_squarefree(n: &BigInt) -> Result<bool> {
        if *n < BigInt::zero() {
            return Err(SrsError::InvalidInput(format!(
                "squarefree test needs n >= 0, got {n}"
            )));
        }
        if *n < BigInt::from(4u8) {
            return Ok(true);
        }
        if (n % BigInt::from(4u8)).is_zero() {
            return Ok(false);
        }
        let root = n.sqrt();
        if &root * &root == *n {
            return Ok(false);
        }

        let mut rem = n.clone();
        let mut p: u64 = 2;
        loop {
            if BigInt::from(p) * BigInt::from(p) > rem {
                return Ok(true);
            }
            if p > TRIAL_DIVISION_LIMIT {
                return Err(SrsError::ArithmeticError(format!(
                    "squarefree test of {n} exceeded the trial-division limit of \
                     {TRIAL_DIVISION_LIMIT}; the input is too large for this check"
                )));
            }
            let bp = BigInt::from(p);
            if (&rem % (&bp * &bp)).is_zero() {
                return Ok(false);
            }
            while (&rem % &bp).is_zero() {
                rem /= &bp;
            }
            p += 1;
        }
    }

    /// Reject a generator set in which one generator is a product of the rest
    /// up to squares, e.g. {2, 3, 6}.
    fn verify_independence(s: &[BigInt]) -> Result<()> {
        let m = s.len();
        for i in 0..m {
            let others: Vec<usize> = (0..m).filter(|j| *j != i).collect();
            for mask in 1..(1usize << others.len()) {
                let mut prod = s[i].clone();
                for (bit_idx, j) in others.iter().enumerate() {
                    if mask & (1 << bit_idx) != 0 {
                        prod *= &s[*j];
                    }
                }
                if Self::is_perfect_square(&prod) {
                    return Err(SrsError::InvalidInput(format!(
                        "squarefree parts not independent: s[{i}]={} times a subset is a perfect square ({prod})",
                        s[i]
                    )));
                }
            }
        }
        Ok(())
    }

    fn is_perfect_square(n: &BigInt) -> bool {
        if *n < BigInt::zero() {
            return false;
        }
        let r = n.sqrt();
        &r * &r == *n
    }

    fn guard_exact(&self, what: &str) -> Result<()> {
        if self.m > EXACT_MAX_GENERATORS {
            return Err(SrsError::ArithmeticError(format!(
                "{what} needs 2^{} field multiplications; m={} exceeds \
                 EXACT_MAX_GENERATORS={EXACT_MAX_GENERATORS}. Returning a \
                 placeholder would be a wrong answer, so this is an error.",
                self.m, self.m
            )));
        }
        Ok(())
    }

    // -- constructors --------------------------------------------------------

    pub fn zero(&self) -> Elem {
        HashMap::new()
    }

    pub fn one(&self) -> Elem {
        [(0usize, BigInt::one())].into_iter().collect()
    }

    pub fn from_rational(&self, q: i64) -> Elem {
        [(0usize, BigInt::from(q))].into_iter().collect()
    }

    /// The rational element `q`, given exactly as a BigInt.
    pub fn from_bigint(&self, q: &BigInt) -> Elem {
        if q.is_zero() {
            HashMap::new()
        } else {
            [(0usize, q.clone())].into_iter().collect()
        }
    }

    pub fn from_sqrt(&self, idx: usize) -> Result<Elem> {
        if idx >= self.m {
            return Err(SrsError::InvalidInput(format!(
                "generator index {idx} out of range [0, {})",
                self.m
            )));
        }
        Ok([(1usize << idx, BigInt::one())].into_iter().collect())
    }

    fn clean(d: Elem) -> Elem {
        d.into_iter().filter(|(_, v)| !v.is_zero()).collect()
    }

    // -- arithmetic ----------------------------------------------------------

    pub fn add(&self, x: &Elem, y: &Elem) -> Elem {
        let mut res = x.clone();
        for (k, v) in y {
            *res.entry(*k).or_insert_with(BigInt::zero) += v;
        }
        Self::clean(res)
    }

    pub fn neg(&self, x: &Elem) -> Elem {
        x.iter().map(|(k, v)| (*k, -v)).collect()
    }

    pub fn sub(&self, x: &Elem, y: &Elem) -> Elem {
        let mut res = x.clone();
        for (k, v) in y {
            *res.entry(*k).or_insert_with(BigInt::zero) -= v;
        }
        Self::clean(res)
    }

    pub fn mul(&self, x: &Elem, y: &Elem) -> Elem {
        let mut res: Elem = HashMap::new();
        for (mask_x, coeff_x) in x {
            for (mask_y, coeff_y) in y {
                let mask_res = mask_x ^ mask_y;
                // sqrt(prod) * sqrt(prod) collapses to the integer product.
                let mut coeff_mult = BigInt::one();
                let mut inter = mask_x & mask_y;
                while inter != 0 {
                    let lsb = inter & inter.wrapping_neg();
                    let idx = lsb.trailing_zeros() as usize;
                    coeff_mult *= &self.s[idx];
                    inter ^= lsb;
                }
                *res.entry(mask_res).or_insert_with(BigInt::zero) +=
                    coeff_x * coeff_y * coeff_mult;
            }
        }
        Self::clean(res)
    }

    pub fn mul_scalar(&self, x: &Elem, scalar: i64) -> Elem {
        Self::clean(x.iter().map(|(k, v)| (*k, v * scalar)).collect())
    }

    pub fn sqr(&self, x: &Elem) -> Elem {
        self.mul(x, x)
    }

    pub fn pow(&self, x: &Elem, n: u32) -> Elem {
        if n == 0 {
            return self.one();
        }
        let mut res = self.one();
        let mut base = x.clone();
        let mut exp = n;
        while exp > 0 {
            if exp & 1 == 1 {
                res = self.mul(&res, &base);
            }
            base = self.mul(&base, &base);
            exp >>= 1;
        }
        res
    }

    /// True when two elements have the same expansion in the basis.
    pub fn eq_elem(&self, x: &Elem, y: &Elem) -> bool {
        let xs: Elem = Self::clean(x.clone());
        let ys: Elem = Self::clean(y.clone());
        xs == ys
    }

    // -- conjugates ----------------------------------------------------------

    /// The 2^m Galois images of `x`, one per sign vector, with multiplicity.
    fn conjugates(&self, x: &Elem) -> Vec<Elem> {
        let mut out = Vec::with_capacity(1usize << self.m);
        for mask in 0..(1usize << self.m) {
            let mut val: Elem = HashMap::new();
            for (basis_mask, coeff) in x {
                let mut sign = 1i8;
                let mut bm = *basis_mask;
                while bm != 0 {
                    let lsb = bm & bm.wrapping_neg();
                    let idx = lsb.trailing_zeros() as usize;
                    if (mask >> idx) & 1 == 1 {
                        sign = -sign;
                    }
                    bm ^= lsb;
                }
                if sign > 0 {
                    val.insert(*basis_mask, coeff.clone());
                } else {
                    val.insert(*basis_mask, -coeff);
                }
            }
            out.push(Self::clean(val));
        }
        out
    }

    /// Field norm: the product of all 2^m Galois conjugates, an integer.
    ///
    /// Errors above `EXACT_MAX_GENERATORS` instead of returning a placeholder.
    pub fn norm(&self, x: &Elem) -> Result<BigInt> {
        if x.is_empty() {
            return Ok(BigInt::zero());
        }
        self.guard_exact("norm")?;

        let mut prod = self.one();
        for conj in self.conjugates(x) {
            prod = self.mul(&prod, &conj);
        }
        self.rational_part(&prod, "norm")
    }

    /// Field trace: the sum of all 2^m Galois conjugates, an integer.
    pub fn trace(&self, x: &Elem) -> Result<BigInt> {
        if x.is_empty() {
            return Ok(BigInt::zero());
        }
        self.guard_exact("trace")?;

        let mut total = BigInt::zero();
        for mask in 0..(1usize << self.m) {
            let mut val = BigInt::zero();
            for (basis_mask, coeff) in x {
                let mut sign = 1i8;
                let mut bm = *basis_mask;
                while bm != 0 {
                    let lsb = bm & bm.wrapping_neg();
                    let idx = lsb.trailing_zeros() as usize;
                    if (mask >> idx) & 1 == 1 {
                        sign = -sign;
                    }
                    bm ^= lsb;
                }
                if sign > 0 {
                    val += coeff;
                } else {
                    val -= coeff;
                }
            }
            total += val;
        }
        Ok(total)
    }

    fn rational_part(&self, x: &Elem, what: &str) -> Result<BigInt> {
        for (mask, coeff) in x {
            if *mask != 0 && !coeff.is_zero() {
                return Err(SrsError::ArithmeticError(format!(
                    "{what} should be rational but has a non-rational part at mask {mask}"
                )));
            }
        }
        Ok(x.get(&0).cloned().unwrap_or_else(BigInt::zero))
    }

    /// The minimal polynomial of `x` over Q, ascending coefficients.
    ///
    /// Product of (t - sigma(x)) over the DISTINCT orbit values.  Getting
    /// "distinct" right is the whole point: multiplying over all 2^m images
    /// with multiplicity yields the characteristic polynomial, which is a
    /// power of this and reports the wrong degree.
    pub fn minimal_polynomial(&self, x: &Elem) -> Result<Vec<BigInt>> {
        if x.is_empty() {
            return Ok(vec![BigInt::zero(), BigInt::one()]);
        }
        self.guard_exact("minimal_polynomial")?;

        // Elem is a HashMap and so is not itself hashable; key the dedup on a
        // canonically ordered Vec of (mask, coefficient) pairs instead.
        // A polynomial coefficient may be a non-integer rational, so build
        // the constant term with Num::from and let the rationals convert
        // exactly; to_i64 would truncate.
        let mut seen: HashSet<Vec<(usize, BigInt)>> = HashSet::new();
        for conj in self.conjugates(x) {
            seen.insert(Self::canon(&conj));
        }
        let roots: Vec<Elem> = seen.into_iter().map(|k| Self::from_canon(&k)).collect();
        self.product_of_linear_factors(&roots)
    }

    /// Canonically ordered coefficient list, usable as a map key.
    fn canon(x: &Elem) -> Vec<(usize, BigInt)> {
        let mut v: Vec<(usize, BigInt)> = x.iter().map(|(k, c)| (*k, c.clone())).collect();
        v.sort_by_key(|(k, _)| *k);
        v
    }

    fn from_canon(v: &[(usize, BigInt)]) -> Elem {
        v.iter().cloned().collect()
    }

    /// The characteristic polynomial of multiplication by `x` on K.
    ///
    /// Product over all 2^m images, with multiplicity.  This is a power of
    /// the minimal polynomial; use `minimal_polynomial` for that.
    pub fn characteristic_polynomial(&self, x: &Elem) -> Result<Vec<BigInt>> {
        if x.is_empty() {
            return Ok(vec![-BigInt::one(), BigInt::one()]);
        }
        self.guard_exact("characteristic_polynomial")?;
        let roots = self.conjugates(x);
        self.product_of_linear_factors(&roots)
    }

    fn product_of_linear_factors(&self, roots: &[Elem]) -> Result<Vec<BigInt>> {
        // Coefficients stay field elements until the final extraction.
        let mut poly: Vec<Elem> = vec![self.one()];
        for conj in roots {
            let neg_conj = self.neg(conj);
            let mut new_poly: Vec<Elem> = vec![self.zero(); poly.len() + 1];
            for (i, coeff) in poly.iter().enumerate() {
                let scaled = self.mul(coeff, &neg_conj);
                new_poly[i] = self.add(&new_poly[i], &scaled);
                new_poly[i + 1] = self.add(&new_poly[i + 1], coeff);
            }
            poly = new_poly;
        }

        let mut out = Vec::with_capacity(poly.len());
        for coeff in &poly {
            out.push(self.rational_part(coeff, "polynomial coefficient")?);
        }
        Ok(primitive(&out))
    }
}

const TRIAL_DIVISION_LIMIT: u64 = 1 << 20;

/// Scale an integer polynomial to be primitive: integer coefficients, gcd 1.
fn primitive(poly: &[BigInt]) -> Vec<BigInt> {
    let mut g = BigInt::zero();
    let mut any = false;
    for c in poly {
        if c.is_zero() {
            continue;
        }
        any = true;
        g = num_integer::gcd(g, c.abs());
    }
    if !any || g <= BigInt::one() {
        return poly.to_vec();
    }
    poly.iter().map(|c| c / &g).collect()
}

impl std::fmt::Display for MultiQuadraticField {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "MultiQuadraticField(s={:?}, m={})", self.s, self.m)
    }
}

impl std::fmt::Debug for MultiQuadraticField {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "MultiQuadraticField(s={:?}, m={})", self.s, self.m)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::ToPrimitive;

    /// `count` distinct primes, which are squarefree and independent.
    fn first_primes(count: usize) -> Vec<BigInt> {
        const PRIMES: [u64; 25] = [
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
            79, 83, 89, 97,
        ];
        PRIMES[..count].iter().map(|p| BigInt::from(*p)).collect()
    }

    fn k2() -> MultiQuadraticField {
        MultiQuadraticField::new(vec![BigInt::from(2u8), BigInt::from(3u8)]).unwrap()
    }

    fn k3() -> MultiQuadraticField {
        MultiQuadraticField::new(vec![
            BigInt::from(2u8),
            BigInt::from(3u8),
            BigInt::from(5u8),
        ])
        .unwrap()
    }

    #[test]
    fn test_field_creation() {
        let k = k2();
        assert_eq!(k.m, 2);
        assert_eq!(
            k.basis_radicand,
            vec![
                BigInt::from(1u8),
                BigInt::from(2u8),
                BigInt::from(3u8),
                BigInt::from(6u8)
            ]
        );
    }

    #[test]
    fn test_basic_arithmetic() {
        let k = k2();
        let a = k.from_sqrt(0).unwrap();
        let b = k.from_sqrt(1).unwrap();
        assert_eq!(k.mul(&a, &b).get(&3), Some(&BigInt::one()));
        assert_eq!(k.mul(&a, &a).get(&0), Some(&BigInt::from(2u8)));
        assert_eq!(k.mul(&b, &b).get(&0), Some(&BigInt::from(3u8)));
    }

    #[test]
    fn test_add_sub_neg() {
        let k = k2();
        let a = k.from_sqrt(0).unwrap();
        let b = k.from_sqrt(1).unwrap();
        assert_eq!(k.add(&a, &b).len(), 2);
        assert!(k.sub(&a, &a).is_empty());
        assert_eq!(k.neg(&a).get(&1), Some(&BigInt::from(-1i8)));
        assert!(k.eq_elem(&k.sub(&a, &a), &k.zero()));
    }

    #[test]
    fn test_norm() {
        let k = k2();
        assert_eq!(k.norm(&k.from_sqrt(0).unwrap()).unwrap(), BigInt::from(4u8));
        assert_eq!(k.norm(&k.one()).unwrap(), BigInt::one());
        assert_eq!(k.norm(&k.zero()).unwrap(), BigInt::zero());
    }

    #[test]
    fn test_norm_errors_above_the_cap_instead_of_guessing() {
        // This is the regression that matters: the old code returned
        // BigInt::one() here, a confident wrong answer.
        let k = MultiQuadraticField::new(first_primes(EXACT_MAX_GENERATORS + 1)).unwrap();
        assert_eq!(k.m, EXACT_MAX_GENERATORS + 1);
        let err = k.norm(&k.one()).unwrap_err();
        assert!(
            err.to_string().contains("EXACT_MAX_GENERATORS"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn test_trace() {
        let k = k2();
        assert_eq!(k.trace(&k.one()).unwrap(), BigInt::from(4u8));
        assert_eq!(k.trace(&k.from_sqrt(0).unwrap()).unwrap(), BigInt::zero());
        // trace(sqrt2 - 1) = 4 * (-1)
        let x = k.sub(&k.from_sqrt(0).unwrap(), &k.one());
        assert_eq!(k.trace(&x).unwrap(), BigInt::from(-4i8));
    }

    #[test]
    fn test_minimal_polynomial_degree_is_right() {
        // sqrt2 in Q(sqrt2,sqrt3) has degree 2, not 4.
        let k = k2();
        let s2 = k.from_sqrt(0).unwrap();
        let min = k.minimal_polynomial(&s2).unwrap();
        assert_eq!(min.len(), 3, "minimal polynomial should be quadratic: {min:?}");
        assert_eq!(min, vec![
            BigInt::from(-2i8),
            BigInt::zero(),
            BigInt::one()
        ]);

        let char = k.characteristic_polynomial(&s2).unwrap();
        assert_eq!(char.len(), 5, "characteristic polynomial should be quartic");
        // (t^2 - 2)^2 = t^4 - 4t^2 + 4
        assert_eq!(char, vec![
            BigInt::from(4u8),
            BigInt::zero(),
            BigInt::from(-4i8),
            BigInt::zero(),
            BigInt::one()
        ]);
    }

    #[test]
    fn test_minimal_polynomial_of_sqrt2_plus_sqrt3() {
        let k = k2();
        let x = k.add(&k.from_sqrt(0).unwrap(), &k.from_sqrt(1).unwrap());
        // t^4 - 10t^2 + 1
        assert_eq!(
            k.minimal_polynomial(&x).unwrap(),
            vec![
                BigInt::from(1u8),
                BigInt::zero(),
                BigInt::from(-10i8),
                BigInt::zero(),
                BigInt::one()
            ]
        );
    }

    #[test]
    fn test_minimal_polynomial_vanishes_on_every_conjugate() {
        // The defining property, checked symbolically in the field.
        let k = k2();
        for x in [
            k.from_sqrt(0).unwrap(),
            k.from_sqrt(1).unwrap(),
            k.sub(&k.from_sqrt(0).unwrap(), &k.one()),
            k.add(&k.from_sqrt(0).unwrap(), &k.from_sqrt(1).unwrap()),
            k.one(),
        ] {
            let poly = k.minimal_polynomial(&x).unwrap();
            for conj in k.conjugates(&x) {
                let mut acc = k.zero();
                for coeff in poly.iter().rev() {
                    let scaled = k.mul(&acc, &conj);
                    acc = k.add(&scaled, &k.from_bigint(coeff));
                }
                assert!(acc.is_empty(), "polynomial does not vanish on a conjugate");
            }
        }
    }

    #[test]
    fn test_three_generators() {
        let k = k3();
        let a = k.from_sqrt(0).unwrap();
        let b = k.from_sqrt(1).unwrap();
        let c = k.from_sqrt(2).unwrap();
        let ab = k.mul(&a, &b);
        assert_eq!(ab.get(&3), Some(&BigInt::one()));
        let abc = k.mul(&ab, &c);
        assert_eq!(abc.get(&7), Some(&BigInt::one()));
    }

    #[test]
    fn test_three_generators_norm() {
        let k = k3();
        assert_eq!(
            k.norm(&k.from_sqrt(0).unwrap()).unwrap(),
            BigInt::from(16u8)
        );
        assert_eq!(k.norm(&k.one()).unwrap(), BigInt::one());
    }

    #[test]
    fn test_independence_rejected() {
        assert!(MultiQuadraticField::new(vec![BigInt::from(2u8), BigInt::from(8u8)]).is_err());
        assert!(
            MultiQuadraticField::new(vec![
                BigInt::from(2u8),
                BigInt::from(3u8),
                BigInt::from(6u8)
            ])
            .is_err()
        );
    }

    #[test]
    fn test_non_squarefree_rejected() {
        assert!(MultiQuadraticField::new(vec![BigInt::from(4u8)]).is_err());
        assert!(MultiQuadraticField::new(vec![BigInt::from(12u8)]).is_err());
    }

    #[test]
    fn test_squarefree_test_matches_definition() {
        // Definition: no prime square divides n.  Cross-checked against a
        // naive trial-division oracle over the whole small range, because
        // is_squarefree's own loop divides as it goes and that interaction
        // is exactly where an off-by-one would hide.
        fn oracle(n: u64) -> bool {
            let mut p = 2u64;
            while p * p <= n {
                if n.is_multiple_of(p * p) {
                    return false;
                }
                p += 1;
            }
            true
        }

        for n in 0u64..5000 {
            let got = MultiQuadraticField::is_squarefree(&BigInt::from(n)).unwrap();
            assert_eq!(got, oracle(n), "n={n}");
        }
    }

    #[test]
    fn test_squarefree_test_on_large_prime_squares() {
        // A square factor much larger than the trial-division reach must
        // still be caught, by the isqrt fast path.
        let p = 1_000_003u64;
        let sq = BigInt::from(p) * BigInt::from(p);
        assert!(MultiQuadraticField::is_squarefree(&BigInt::from(p)).unwrap());
        assert!(!MultiQuadraticField::is_squarefree(&sq).unwrap());
    }

    #[test]
    fn test_squarefree_test_rejects_negative() {
        assert!(MultiQuadraticField::is_squarefree(&BigInt::from(-1i8)).is_err());
    }

    #[test]
    fn test_every_exact_op_guards_above_the_cap() {
        // The regression: `norm` used to return BigInt::one() above m = 6.
        // Every exponential operation must error there, not answer.
        let m = EXACT_MAX_GENERATORS + 1;
        let primes: Vec<BigInt> = first_primes(m);
        let k = MultiQuadraticField::new(primes).unwrap();
        let x = k.one();

        for name in ["norm", "trace", "minimal_polynomial", "characteristic_polynomial"] {
            let result = match name {
                "norm" => k.norm(&x).map(|_| ()),
                "trace" => k.trace(&x).map(|_| ()),
                "minimal_polynomial" => k.minimal_polynomial(&x).map(|_| ()),
                "characteristic_polynomial" => k.characteristic_polynomial(&x).map(|_| ()),
                _ => unreachable!(),
            };
            let err = result.expect_err(&format!("{name} should refuse above the cap"));
            assert!(
                err.to_string().contains("EXACT_MAX_GENERATORS"),
                "{name} error was not the expected cap message: {err}"
            );
        }
    }

    #[test]
    fn test_exact_ops_work_at_the_cap() {
        let k = MultiQuadraticField::new(first_primes(EXACT_MAX_GENERATORS)).unwrap();
        assert_eq!(k.m, EXACT_MAX_GENERATORS);
        let x = k.one();
        assert_eq!(k.norm(&x).unwrap(), BigInt::one());
        assert_eq!(k.trace(&x).unwrap(), BigInt::from(64u8));
        assert!(k.minimal_polynomial(&x).is_ok());
    }

    #[test]
    fn test_generator_count_capped_before_allocating() {
        // 2^21 entries would be allocated by a revision with no cap.
        let err =
            MultiQuadraticField::new(first_primes(MAX_GENERATORS + 1)).unwrap_err();
        assert!(err.to_string().contains("MAX_GENERATORS"), "{err}");
    }

    #[test]
    fn test_from_sqrt_range_checked() {
        let k = k2();
        assert!(k.from_sqrt(2).is_err());
        assert!(k.from_sqrt(usize::MAX).is_err());
    }

    #[test]
    fn test_empty_and_small_inputs_rejected() {
        assert!(MultiQuadraticField::new(vec![]).is_err());
        assert!(MultiQuadraticField::new(vec![BigInt::from(1u8)]).is_err());
        assert!(MultiQuadraticField::new(vec![BigInt::from(0u8)]).is_err());
    }

    #[test]
    fn test_minimal_polynomial_degree_tracks_orbit_size() {
        // The defining relationship: the minimal polynomial's degree is the
        // number of distinct Galois images.  Characteristic polynomial degree
        // is always 2^m.  Getting this pair confused is what the old code did.
        for m in 1..=3usize {
            let k = MultiQuadraticField::new(first_primes(m)).unwrap();
            let x = k.one();
            assert_eq!(
                k.characteristic_polynomial(&x).unwrap().len(),
                (1usize << m) + 1,
                "m={m}: characteristic polynomial should have degree 2^m"
            );
            // 1 is rational, so its orbit is a single point regardless of m.
            assert_eq!(
                k.minimal_polynomial(&x).unwrap().len(),
                2,
                "m={m}: the minimal polynomial of 1 is linear"
            );
        }
    }

    #[test]
    fn test_minimal_polynomial_of_rational_has_degree_one() {
        let k = k3();
        for q in 1..=9i64 {
            let poly = k.minimal_polynomial(&k.from_rational(q)).unwrap();
            assert_eq!(poly.len(), 2, "minimal polynomial of {q} should be linear");
            assert_eq!(poly, vec![BigInt::from(-q), BigInt::one()]);
        }
    }

    #[test]
    fn test_arithmetic_agrees_with_numeric_evaluation() {
        // Check field algebra against real arithmetic on sqrt2 and sqrt3,
        // to 12 places.  This is the only float in the crate and it is a
        // test oracle, not part of any decision path.
        fn approx(k: &MultiQuadraticField, x: &Elem) -> f64 {
            let mut total = 0.0f64;
            for (mask, coeff) in x {
                let mut rad = 1.0f64;
                let mut bm = *mask;
                while bm != 0 {
                    let lsb = bm & bm.wrapping_neg();
                    let idx = lsb.trailing_zeros() as usize;
                    rad *= &k.s[idx]
                        .to_f64()
                        .expect("generator fits in f64 for this test");
                    bm ^= lsb;
                }
                total += coeff.to_f64().expect("coefficient fits in f64") * rad.sqrt();
            }
            total
        }

        let k = k2();
        let s2 = k.from_sqrt(0).unwrap();
        let s3 = k.from_sqrt(1).unwrap();
        let cases = vec![
            k.add(&s2, &s3),
            k.sub(&s2, &s3),
            k.mul(&s2, &s3),
            k.sqr(&s2),
            k.add(&k.mul(&s2, &s3), &k.from_rational(5)),
        ];
        for x in cases {
            let numeric = approx(&k, &x);
            // (numeric)^2 - 2 and (numeric)^2 - 3 are the relations that hold
            // for the basis elements; verify the element squares to what the
            // algebra says it should.
            let squared = approx(&k, &k.sqr(&x));
            assert!(
                (numeric * numeric - squared).abs() < 1e-9,
                "mul-then-square disagrees with squaring numerically"
            );
        }
    }

    #[test]
    fn test_primitive_normalises_content() {
        // (t^2-2)^2 = t^4-4t^2+4 has content 1 already; check a genuinely
        // non-primitive input is scaled down rather than left as-is.
        let k = k2();
        let x = k.mul_scalar(&k.from_sqrt(0).unwrap(), 2);
        let poly = k.minimal_polynomial(&x).unwrap();
        let leading = poly.last().unwrap();
        assert_eq!(leading, &BigInt::one(), "leading coefficient should be 1");
    }

    #[test]
    fn test_rejects_bad_a() {
        for bad in [0u64, 1] {
            assert!(pell_convergents_rejects(bad), "a={bad} should be rejected");
        }
    }

    fn pell_convergents_rejects(a: u64) -> bool {
        crate::f1::pell_convergents(a, 3).is_err()
    }

    #[test]
    fn test_pow() {
        let k = k2();
        let a = k.from_sqrt(0).unwrap();
        assert_eq!(k.pow(&a, 0), k.one());
        assert_eq!(k.pow(&a, 1), a);
        assert_eq!(k.pow(&a, 2).get(&0), Some(&BigInt::from(2u8)));
        assert_eq!(k.pow(&a, 3).get(&1), Some(&BigInt::from(2u8)));
    }

    #[test]
    fn test_mul_scalar() {
        let k = k2();
        let a = k.from_sqrt(0).unwrap();
        assert_eq!(k.mul_scalar(&a, 5).get(&1), Some(&BigInt::from(5u8)));
        assert!(k.mul_scalar(&a, 0).is_empty());
    }
}
