//! Precision schedules for iterative square-root algorithms.
//!
//! Two schedules are modelled, both returning the *working precision* (in
//! decimal digits) used by each level of the recursion.
//!
//! 1. **Cubic bracket-zoom** — the method of Papers 1 and 2. Each level gains
//!    three digits, so the per-level budget is
//!    `prec_n = MIN + 3*width_digits(w_n) + 3*adjusted(a)`, capped at the
//!    target. This is the schedule validated in Paper1 section 10; the values
//!    this module produces are checked against that reference in the tests.
//!
//! 2. **Radix-2^b restoring** — the standard digit-recurrence square root. Each
//!    step resolves `b` bits, so the number of steps is about `D*log2(10)/b`
//!    and every step runs at (roughly) the full width. This is the honest
//!    baseline for a hardware square root, and it is what shipped FPGA/ASIC
//!    SQR units actually implement.
//!
//! The whole point of the cost model is to compare these two on equal terms, so
//! both are expressed as a `Vec<usize>` of per-level working precisions.

/// The width in decimal digits of a quantity of size ~`w`. Mirrors
/// `exp_precision._width_digits`, which is `1 - w.adjusted()`.
///
/// `w_adjusted` is the base-10 exponent of the leading digit, tracked exactly
/// as an `i64` rather than computed from an `f64` logarithm. That distinction
/// is not cosmetic: an `f64` recurrence drifts, and at D = 3000 it drifts far
/// enough to change the schedule (327 -> 321, 915 -> 900, 2685 -> 2640), which
/// would silently de-synchronise this model from the Python implementation it
/// claims to mirror. The whole schedule is a chain of `floor(log10(.))`
/// decisions, so it is computed on exact exponents throughout.
pub fn width_digits(w_adjusted: i64) -> i64 {
    1 - w_adjusted
}

/// Configuration for the cubic zoom schedule.
#[derive(Debug, Clone, Copy)]
pub struct ZoomConfig {
    /// Target decimal digits, `D`.
    pub digits: u32,
    /// Guard bits/digits added to the target.
    pub guard: u32,
    /// Floor on the working precision.
    pub min_prec: u32,
    /// Number of radicands to average over.
    ///
    /// The anchor `a = floor(sqrt(alpha))` sets both the seed magnitude and
    /// the number of levels, so the cost is not a function of `D` alone.
    pub num_alphas: usize,
    /// Radicands to use.
    pub alphas: &'static [u64],
}

impl Default for ZoomConfig {
    fn default() -> Self {
        ZoomConfig {
            digits: 3000,
            guard: 15,
            min_prec: 30,
            num_alphas: 16,
            alphas: &[2, 3, 5, 7, 10, 12, 13, 17, 23, 31, 101, 257, 1009, 7919, 104729, 999983],
        }
    }
}

/// Per-level working precisions of the cubic zoom for one radicand.
///
/// Mirrors the loop in `exp_precision.zoom(scheduled=True)`:
///
/// ```text
/// prec  = min(MIN + 3*width_digits(w) + 3*adjusted(a), digits + guard)
/// eps   = w^3 / (27 a^2)
/// stop  when eps.adjusted() < -digits
/// ```
///
/// plus one final full-precision recompute (the terminating level), which
/// Paper1 section 10 counts as a separate step.
///
/// The bracket values themselves are not simulated. Only `adjusted(w)` is
/// propagated, because that is the only property of the bracket the schedule
/// reads. `a` is the seed `floor(sqrt(alpha))` and its `adjusted()` is stable
/// across the recursion, since every iterate stays within a factor of two of
/// the seed.
pub fn cubic_schedule(cfg: &ZoomConfig, alpha: u64) -> Vec<u32> {
    let a_seed = isqrt_floor(alpha);
    let adj_a = adjusted_int(a_seed);
    // The bracket width starts at 1: exponent 0, significand 1.
    let mut w_adj: i64 = 0;
    let mut w_sig: f64 = 1.0;
    let mut out = Vec::new();
    let cap = cfg.digits + cfg.guard;

    // The anchor `a` is *not* constant. Python's loop sets `a = +est - eps` each
    // level, so it migrates from the integer seed `floor(sqrt(alpha))` towards
    // `sqrt(alpha)`. That drift is not cosmetic: for alpha = 12 the seed is 3
    // but the limit is 3.4641, and since `eps = w^3/(27 a^2)` the denominator
    // grows by 27*12 - 27*9 = 81, which is 1.9 digits. Holding `a` fixed
    // therefore mis-schedules every level from the third onward.
    //
    // `a` is carried as a significand in [1,10) plus the integer exponent
    // `adj_a`, which is stable because `a` stays within a factor of two of the
    // seed. Every exponent in the loop below is exact integer arithmetic; only
    // significands are floating point, and each is renormalised into [1,10)
    // every step, so rounding cannot accumulate the way a running f64
    // logarithm does.
    let scale = 10f64.powi(adj_a as i32);
    let root_sig = (alpha as f64).sqrt() / scale;
    let mut a_sig = a_seed as f64 / scale;

    for _ in 0..200 {
        let prec = (cfg.min_prec as i64 + 3 * width_digits(w_adj) + 3 * adj_a)
            .clamp(0, cap as i64) as u32;
        out.push(prec);

        // eps = [w_sig^3 / (27 a_sig^2)] * 10^(3 w_adj - 2 adj_a)
        let m = w_sig * w_sig * w_sig / (27.0 * a_sig * a_sig);
        let eps_log10 = 3.0 * w_adj as f64 - 2.0 * adj_a as f64 + m.log10();
        let eps_adj = eps_log10.floor() as i64;

        if eps_adj < -(cfg.digits as i64) {
            out.push(cfg.digits + cfg.guard + 10);
            return out;
        }
        // Renormalise eps into a significand, then advance the anchor the way
        // Python does: a <- est - eps, with est -> sqrt(alpha).
        let eps_sig = 10f64.powf(eps_log10 - eps_adj as f64);
        // eps is eps_sig * 10^eps_adj in absolute scale. The anchor `a` is
        // measured in units of 10^adj_a, so dividing through by that scale is
        // all that is required here.
        let eps_in_a_scale = eps_sig * 10f64.powi((eps_adj - adj_a) as i32);
        a_sig = root_sig - eps_in_a_scale;
        // Guard against a_sig leaving [1,10) on the first step.
        if a_sig >= 10.0 {
            a_sig /= 10.0;
        } else if a_sig < 1.0 {
            a_sig *= 10.0;
        }
        w_sig = eps_sig;
        w_adj = eps_adj;
    }
    out
}

/// Leading decimal digit of an integer, as a base-10 exponent.
fn adjusted_int(n: u64) -> i64 {
    if n == 0 {
        return 0;
    }
    (n as f64).log10().floor() as i64
}

/// Per-level working precisions of a radix-`2^b` restoring digit recurrence.
///
/// Each step produces `b` bits, so `steps ~= ceil(D * log2(10) / b)`. Every
/// step is taken at the full target width, because the remainder register must
/// be able to hold the partial square root at all times.
pub fn restoring_schedule(cfg: &ZoomConfig, radix_bits: u32) -> Vec<u32> {
    let bits = (cfg.digits as f64) * std::f64::consts::LOG2_10;
    let steps = (bits / radix_bits as f64).ceil() as usize;
    let width = cfg.digits + cfg.guard;
    vec![width; steps]
}

/// Integer square root, `floor(sqrt(n))`. Used only for seeding; `u128`
/// intermediate keeps it exact up to 2^64.
pub fn isqrt_floor(n: u64) -> u64 {
    if n < 2 {
        return n;
    }
    let mut x = n as u128;
    let mut y = (x + 1) / 2;
    while y < x {
        x = y;
        y = (x + n as u128 / x) / 2;
    }
    x as u64
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Paper1 section 10 reports, for alpha = 12 at D = 3000, the schedule
    /// 33, 42, 63, 129, 327, 915, 2685, 3015, 3025. Our loop appends the
    /// terminating recompute as `digits + guard + 10` = 3025, so the expected
    /// vector is the first eight entries plus 3025. Guard against drift in the
    /// reference itself by asserting the head of the sequence.
    #[test]
    fn cubic_schedule_matches_paper1_reference() {
        let cfg = ZoomConfig { digits: 3000, ..Default::default() };
        let sched = cubic_schedule(&cfg, 12);
        assert_eq!(
            &sched[..8],
            &[33, 42, 63, 129, 327, 915, 2685, 3015],
            "head of the cubic schedule must match Paper1 section 10"
        );
        assert_eq!(sched.last().copied(), Some(3025));
    }

    #[test]
    fn cubic_depth_is_log3_like() {
        let cfg = ZoomConfig { digits: 3000, ..Default::default() };
        for alpha in [2, 12, 999983] {
            let d = cubic_schedule(&cfg, alpha).len();
            assert!(
                (6..=10).contains(&d),
                "alpha {alpha}: depth {d} outside the expected log3(D) range"
            );
        }
    }

    #[test]
    fn cubic_schedule_is_nondecreasing() {
        let cfg = ZoomConfig { digits: 3000, ..Default::default() };
        for alpha in [2, 3, 12, 101, 999983] {
            let s = cubic_schedule(&cfg, alpha);
            for w in s.windows(2) {
                assert!(w[1] >= w[0], "alpha {alpha}: schedule must not shrink");
            }
        }
    }

    #[test]
    fn restoring_depth_scales_inversely_with_radix() {
        let cfg = ZoomConfig { digits: 3000, ..Default::default() };
        let d2 = restoring_schedule(&cfg, 2).len();
        let d4 = restoring_schedule(&cfg, 4).len();
        let d8 = restoring_schedule(&cfg, 8).len();
        assert_eq!(d4, (d2 + 1) / 2, "radix 4 halves the step count");
        assert_eq!(d8, (d4 + 1) / 2, "radix 8 halves it again");
        // 3000 decimal digits is 9965.78 bits, and a radix-2 step resolves one
        // bit, so 4983 steps. Note this is still ~550x the zoom's depth, which
        // is the whole reason the zoom is interesting as hardware.
        assert_eq!(d2, 4983, "radix-2 restoring at D=3000");
        assert!(d2 > 400 * cubic_schedule(&cfg, 12).len());
    }

    #[test]
    fn isqrt_is_correct() {
        for n in [0u64, 1, 2, 3, 4, 12, 15, 16, 99, 101, 999983, 10u64.pow(12)] {
            let r = isqrt_floor(n);
            assert!(r * r <= n && (r + 1) * (r + 1) > n, "isqrt wrong at {n}");
        }
    }
}
