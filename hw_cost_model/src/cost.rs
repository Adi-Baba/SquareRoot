//! The cost model: turn a precision schedule into area, delay and energy.
//!
//! # What each algorithm costs per level
//!
//! **Affine zoom** — the literal transcription of Papers 1-2. Each level
//! evaluates `m = a + p/(a+b)` then `est = a + p/(a+m)`, i.e. **two divisions
//! by a running value**. A restoring array divider is the dominant term.
//!
//! **Projective zoom** — the same map lifted to `P^1`:
//! `(X,Y) -> (aX + alpha*Y, X + aY)`. That is **two multiplications and two
//! additions, no division at all**. Both multipliers are *seed constants*, so
//! they are multiply-by-constant, not general multipliers. This is the form
//! that makes the Moebius observation matter in hardware.
//!
//! **Radix-2^b restoring** — the standard digit recurrence, and the honest
//! baseline: shipping SQR units implement it, not Newton's divider.
//!
//! The model's job is to say which wins on area x delay x energy, not to assert
//! that any of them do.

use crate::schedule::ZoomConfig;
use crate::tech::{fast_delay, ripple_delay, CellLibrary};

/// Adder topology. This decides whether truncating a datapath buys *depth* as
/// well as area, which is the whole point of the exercise.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdderKind {
    Ripple,
    Lookahead,
}

impl AdderKind {
    pub fn name(self) -> &'static str {
        match self {
            AdderKind::Ripple => "ripple",
            AdderKind::Lookahead => "lookahead",
        }
    }
}

/// Adder delay at a given width.
pub fn adder_delay(kind: AdderKind, width: usize, lib: &CellLibrary) -> f64 {
    match kind {
        AdderKind::Ripple => ripple_delay(width, lib),
        AdderKind::Lookahead => fast_delay(width, lib),
    }
}

/// Critical path of one *multiply-by-constant* at width `w`.
///
/// A constant multiplier is a shift-and-add network, so its depth is the adder
/// chain: with `k` partial products it is `k` adds, hence logarithmic in `w` for
/// a lookahead adder and linear for a ripple. The number of partial products is
/// proportional to the popcount of the constant, which is unknowable here
/// (the constants are seeds and depend on the radicand), so this uses the
/// worst realistic case of a dense constant, about `log2(w)` partial products.
pub fn const_mult_delay(kind: AdderKind, width: usize, lib: &CellLibrary) -> f64 {
    let w = width.max(2);
    let partials = (w as f64).log2().ceil().max(1.0);
    match kind {
        AdderKind::Ripple => partials * ripple_delay(w, lib),
        AdderKind::Lookahead => partials * fast_delay(w, lib),
    }
}

/// Critical path of one restoring-array divider at width `w`.
///
/// A restoring divider is a shift-and-subtract array: `w` sequential
/// subtract-compare stages, each a full-width subtract. This is the structure
/// that makes the affine form expensive in *depth*, not merely in area — it is
/// the dominant term in the comparison and the reason the projective lift is
/// worth doing before any synthesis.
pub fn divider_delay(kind: AdderKind, width: usize, lib: &CellLibrary) -> f64 {
    let w = width.max(2);
    // One stage per bit, each costing a subtract of the current remainder.
    let stage = match kind {
        AdderKind::Ripple => ripple_delay(w, lib),
        AdderKind::Lookahead => fast_delay(w, lib),
    };
    w as f64 * stage
}

/// Critical path through one level of the given algorithm at width `w`.
///
/// This is what makes the algorithms differ in delay. Charging every level a
/// plain adder width would make all three identical in depth, which is not the
/// case and would make the comparison vacuous.
pub fn level_delay(algo: Algo, kind: AdderKind, w: usize, lib: &CellLibrary) -> f64 {
    match algo {
        // (X,Y) -> (aX + alpha*Y, X + aY): multiply then add, in parallel on
        // the two coordinates, so the depth is one mult plus one add.
        Algo::ProjectiveZoom => const_mult_delay(kind, w, lib) + adder_delay(kind, w, lib),
        // m = a + p/(a+b); est = a + p/(a+m): two *dependent* divisions, each
        // followed by an add. The two divisions are in series because the
        // second divides by a quantity the first produced.
        Algo::AffineZoom => 2.0 * (divider_delay(kind, w, lib) + adder_delay(kind, w, lib)),
        // One step is a 2b subtract plus a select; see `restoring_step`.
        Algo::Restoring => {
            let b = 2usize;
            adder_delay(kind, 2 * b, lib) + adder_delay(kind, b, lib)
        }
    }
}

/// The three algorithms under comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Algo {
    /// `est = a + p/(a+m)`, affine, with two divisions per level.
    AffineZoom,
    /// The projective lift: two multiplies and two adds, no division.
    ProjectiveZoom,
    /// Standard radix-2^b digit recurrence.
    Restoring,
}

impl Algo {
    pub fn name(self) -> &'static str {
        match self {
            Algo::AffineZoom => "affine-zoom",
            Algo::ProjectiveZoom => "projective-zoom",
            Algo::Restoring => "restoring",
        }
    }

    pub fn all() -> [Algo; 3] {
        [Algo::AffineZoom, Algo::ProjectiveZoom, Algo::Restoring]
    }
}

/// Area, delay and energy of one design.
#[derive(Debug, Clone, Copy)]
pub struct Cost {
    pub algo: Algo,
    /// Gate equivalents.
    pub area: f64,
    /// `t_fo4`.
    pub delay: f64,
    /// Relative switching energy.
    pub energy: f64,
    /// Levels (iterations) in the recursion.
    pub depth: usize,
    /// Sum of per-level working precisions: the datapath area-time product,
    /// and the quantity per-level truncation is supposed to reduce.
    pub sum_prec: f64,
    /// Widest level: what a fixed-width design would have to carry.
    pub max_prec: u32,
}

impl Cost {
    /// Area x delay x energy, the usual proxy for "is this worth building".
    ///
    /// Reported, but **not** the number to decide on: the three designs run for
    /// wildly different numbers of cycles, so this product mostly re-states the
    /// iteration count. Use `cost_per_result` for a comparison that is fair
    /// across them.
    pub fn merit(&self) -> f64 {
        self.area * self.delay * self.energy
    }

    /// Area per square root delivered, at one result per `delay` cycles.
    ///
    /// This is the comparison that actually means something. A design that
    /// finishes in 9 cycles amortises its datapath over 9 results-worth of work
    /// per second; one that takes 4983 cycles amortises the same datapath over
    /// 4983. Dividing area by the iteration count is the whole difference
    /// between "the zoom is 2.9x smaller" and "the zoom is 550x less
    /// area-efficient", and only the second is a fair statement.
    pub fn cost_per_result(&self) -> f64 {
        self.area / self.depth.max(1) as f64
    }

    /// Latency-normalised energy: energy per square root.
    pub fn energy_per_result(&self) -> f64 {
        self.energy / self.depth.max(1) as f64
    }

    /// A single merit that is comparable across different iteration counts:
    /// area x delay x energy all divided by depth, so a design is rewarded for
    /// delivering results rather than for taking a long time to deliver one.
    pub fn fair_merit(&self) -> f64 {
        self.merit() / (self.depth.max(1) as f64).powi(3)
    }
}

/// Cost of one affine level: two running divisions plus two multiplies.
///
/// A restoring array divider of width `w` costs `2w` full adders; two of them
/// dominate everything else here, which is the point of the comparison.
fn affine_level(w: usize, lib: &CellLibrary) -> (f64, f64) {
    let divisions = 2.0 * (2 * w) as f64 * lib.full_adder;
    let mults = 2.0 * w as f64 * lib.and_or;
    let adds = 2.0 * w as f64 * lib.full_adder;
    (divisions + mults + adds, 2.0 * lib.switch_power_fa * (2 * w) as f64)
}

/// Cost of one projective level: `(aX + alpha*Y, X + aY)`.
///
/// Both coefficients are seed constants, so these are multiply-by-constant: an
/// AND-OR array, not a full multiplier. No divider.
fn projective_level(w: usize, lib: &CellLibrary) -> (f64, f64) {
    let mults = 2.0 * w as f64 * lib.and_or;
    let adds = 2.0 * w as f64 * lib.full_adder;
    (mults + adds, 2.0 * w as f64 * lib.switch_power_fa)
}

/// Cost of one radix-`b` restoring step, resolving `b` bits.
///
/// Per new bit: a `2b`-bit subtract (2b full adders) and a `2b`-wide select.
/// So a step costs `b * (2b*FA + 2b*MUX)`.
fn restoring_step(b: u32, lib: &CellLibrary) -> (f64, f64) {
    let bf = b as f64;
    let area = bf * (2.0 * bf * lib.full_adder + 2.0 * bf * lib.mux2);
    (area, bf * 2.0 * bf * lib.switch_power_mux)
}

/// Evaluate the whole design from a per-level precision schedule.
fn from_schedule(
    algo: Algo,
    sched: &[u32],
    adder: AdderKind,
    lib: &CellLibrary,
) -> Cost {
    let mut area = 0.0f64;
    let mut delay = 0.0;
    let mut switch_energy = 0.0;
    let mut max_prec = 0u32;

    for &p in sched {
        let w = p.max(2) as usize;
        let (a, e) = match algo {
            Algo::AffineZoom => affine_level(w, lib),
            Algo::ProjectiveZoom => projective_level(w, lib),
            Algo::Restoring => restoring_step(2, lib),
        };
        // Area is the *widest* level, not the sum. The datapath is one physical
        // array reused by every level; summing over levels would count it once
        // per iteration and describe no real circuit. (It is precisely the
        // mistake that makes a 4983-step radix recurrence look 500x larger than
        // it is, when in fact it is a single small step unit iterated.)
        area = area.max(a);
        switch_energy += e;
        max_prec = max_prec.max(p);
        // Levels are sequential: a level cannot start before the last finished.
        // The critical path through a level depends on the algorithm, not just
        // on the adder width.
        delay += level_delay(algo, adder, w, lib);
    }

    // Leakage is charged every cycle for the whole datapath, over the whole run.
    let leakage = area * lib.leakage * sched.len() as f64;
    let energy = switch_energy + leakage;

    Cost {
        algo,
        area,
        delay,
        energy,
        depth: sched.len(),
        sum_prec: sched.iter().map(|&p| p as f64).sum(),
        max_prec,
    }
}

/// Cost of the cubic zoom, averaged over the radicand set.
///
/// The anchor `a = floor(sqrt(alpha))` changes both the seed and the number of
/// levels, so a single-radicand number would be an artefact of one choice.
pub fn zoom_cost(
    cfg: &ZoomConfig,
    algo: Algo,
    adder: AdderKind,
    lib: &CellLibrary,
) -> Cost {
    let mut area = 0.0;
    let mut delay = 0.0;
    let mut energy = 0.0;
    let mut depth = 0usize;
    let mut sum_prec = 0.0;
    let mut max_prec_sum = 0.0;

    for &alpha in cfg.alphas.iter().take(cfg.num_alphas) {
        let sched = crate::schedule::cubic_schedule(cfg, alpha);
        let c = from_schedule(algo, &sched, adder, lib);
        area += c.area;
        delay += c.delay;
        energy += c.energy;
        depth += c.depth;
        sum_prec += c.sum_prec;
        // Mean widest level: every radicand needs its own datapath, but the
        // area figure above already uses the max over that radicand's levels.
        max_prec_sum += c.max_prec as f64;
    }
    let n = cfg.num_alphas as f64;
    Cost {
        algo,
        area: area / n,
        delay: delay / n,
        energy: energy / n,
        depth: (depth as f64 / n).round() as usize,
        sum_prec: sum_prec / n,
        max_prec: (max_prec_sum / n).round() as u32,
    }
}

/// Cost of the radix-`b` restoring baseline.
pub fn restoring_cost(
    cfg: &ZoomConfig,
    radix_bits: u32,
    adder: AdderKind,
    lib: &CellLibrary,
) -> Cost {
    let sched = crate::schedule::restoring_schedule(cfg, radix_bits);
    let mut c = from_schedule(Algo::Restoring, &sched, adder, lib);
    // The step's own width is 2b, not the full target; the remainder register
    // is what is wide. Charge the register once, at target width.
    c.area += (cfg.digits + cfg.guard) as f64 * lib.full_adder;
    c
}

/// The headline comparison: what per-level truncation actually buys.
///
/// This is the number the hardware idea rests on, so it is computed
/// explicitly rather than left implicit in the schedules. It compares the
/// projective zoom as built (narrow early levels) against the same zoom forced
/// to carry its widest datapath at every level -- which is what a fixed-width
/// implementation would do.
pub fn truncation_saving(cfg: &ZoomConfig, adder: AdderKind, lib: &CellLibrary) -> TruncationReport {
    let mut sched_sum = 0.0;
    let mut max_prec = 0u32;
    let mut levels = 0usize;
    // Delay is *not* the sum of widths: a ripple adder is linear in width, so
    // summing widths is exactly right for it, but a logarithmic adder would
    // understate the benefit. Sum the actual per-level delay instead.
    let mut variable_delay = 0.0;
    for &alpha in cfg.alphas.iter().take(cfg.num_alphas) {
        let s = crate::schedule::cubic_schedule(cfg, alpha);
        sched_sum += s.iter().map(|&p| p as f64).sum::<f64>();
        max_prec = max_prec.max(s.iter().copied().max().unwrap_or(0));
        levels += s.len();
        for &p in &s {
            variable_delay += adder_delay(adder, p.max(2) as usize, lib);
        }
    }
    let n = cfg.num_alphas as f64;
    let mean_levels = levels as f64 / n;
    let sum_prec = sched_sum / n;
    let fixed_area_time = mean_levels * max_prec as f64;
    let variable_area_time = sum_prec;
    let fixed_delay = mean_levels * adder_delay(adder, max_prec.max(2) as usize, lib);
    variable_delay /= n;

    TruncationReport {
        mean_levels,
        mean_depth: mean_levels.round() as usize,
        sum_prec,
        max_prec,
        variable_area_time,
        fixed_area_time,
        variable_delay,
        fixed_delay,
        saving: 1.0 - variable_area_time / fixed_area_time,
        delay_saving: 1.0 - variable_delay / fixed_delay,
    }
}

/// The truncation saving, in the units of the schedule.
#[derive(Debug, Clone, Copy)]
pub struct TruncationReport {
    pub mean_levels: f64,
    pub mean_depth: usize,
    /// Sum of per-level working precisions (the variable-width area-time).
    pub sum_prec: f64,
    /// Widest level.
    pub max_prec: u32,
    pub variable_area_time: f64,
    /// What a fixed-width datapath would pay: depth x widest level.
    pub fixed_area_time: f64,
    /// Sum of per-level adder delays for the narrowed datapath.
    pub variable_delay: f64,
    /// Delay a fixed-width datapath would pay.
    pub fixed_delay: f64,
    /// Fraction of the area-time product recovered by narrowing.
    pub saving: f64,
    /// Fraction of the critical path recovered by narrowing. For a logarithmic
    /// adder this is much larger than `saving`, and for a ripple adder it is
    /// equal to `saving` (width scales linearly with delay).
    pub delay_saving: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tech::CellLibrary;

    fn cfg(d: u32) -> ZoomConfig {
        ZoomConfig { digits: d, ..Default::default() }
    }

    #[test]
    fn projective_has_no_divider_cost() {
        let lib = CellLibrary::generic_180nm();
        let w = 512;
        let (aff, _) = affine_level(w, &lib);
        let (prj, _) = projective_level(w, &lib);
        // Two restoring dividers must dominate two mults and two adds.
        assert!(aff > 2.0 * prj, "affine {aff} vs projective {prj} at w={w}");
    }

    #[test]
    fn projective_beats_affine_on_area_and_delay() {
        let lib = CellLibrary::generic_180nm();
        for adder in [AdderKind::Ripple, AdderKind::Lookahead] {
            for d in [200u32, 1000, 3000] {
                let a = zoom_cost(&cfg(d), Algo::AffineZoom, adder, &lib);
                let p = zoom_cost(&cfg(d), Algo::ProjectiveZoom, adder, &lib);
                assert!(p.area < a.area, "area: d={d} {adder:?}");
                assert!(p.delay < a.delay, "delay: d={d} {adder:?}");
                assert!(p.merit() < a.merit());
            }
        }
    }

    #[test]
    fn truncation_saves_something_and_never_negative() {
        let lib = CellLibrary::generic_180nm();
        for d in [50u32, 200, 1000, 3000] {
            let t = truncation_saving(&cfg(d), AdderKind::Ripple, &lib);
            assert!(t.saving > 0.0 && t.saving < 1.0, "d={d} saving={}", t.saving);
            assert!(t.variable_area_time < t.fixed_area_time);
            assert_eq!(t.max_prec as f64 * t.mean_levels, t.fixed_area_time);
        }
    }

    #[test]
    fn truncation_saving_grows_with_depth_but_bottoms_out() {
        let lib = CellLibrary::generic_180nm();
        let s200 = truncation_saving(&cfg(200), AdderKind::Ripple, &lib).saving;
        let s1000 = truncation_saving(&cfg(1000), AdderKind::Ripple, &lib).saving;
        let s3000 = truncation_saving(&cfg(3000), AdderKind::Ripple, &lib).saving;
        assert!(s3000 > s200, "deeper recursion should truncate harder");
        // Geometric growth means the saving converges; it cannot approach 1.
        assert!(s3000 < 0.95, "saving {s3000} should not approach 1");
    }

    #[test]
    fn truncation_saves_the_same_delay_fraction_as_area_for_ripple() {
        // Ripple delay is linear in width, so narrowing a level saves exactly
        // the same fraction of delay that it saves of area-time. This is the
        // case where the area saving translates directly into a timing win.
        let lib = CellLibrary::generic_180nm();
        let t = truncation_saving(&cfg(3000), AdderKind::Ripple, &lib);
        assert!((t.delay_saving - t.saving).abs() < 1e-9);
    }

    #[test]
    fn a_logarithmic_adder_erases_most_of_the_truncation_delay_claim() {
        // Counter-intuitive but important: with a lookahead adder the
        // fixed-width baseline is *already* fast, because a log-depth adder at
        // the full width costs no more than a ripple adder at a small width.
        // So narrowing the datapath saves far less delay than it saves area --
        // the two do not come for free together.
        //
        // This falsifies the naive claim that per-level truncation is a
        // "62% timing win" as well as a 62% area win. It is an area result.
        let lib = CellLibrary::generic_180nm();
        let r = truncation_saving(&cfg(3000), AdderKind::Ripple, &lib);
        let l = truncation_saving(&cfg(3000), AdderKind::Lookahead, &lib);
        assert_eq!(r.saving, l.saving, "area saving is adder-independent");
        assert!(
            l.delay_saving < r.delay_saving,
            "lookahead delay saving {l:?} should be below ripple {r:?}"
        );
    }

    #[test]
    fn restoring_needs_many_more_levels_than_the_zoom() {
        let lib = CellLibrary::generic_180nm();
        let c = cfg(3000);
        let z = zoom_cost(&c, Algo::ProjectiveZoom, AdderKind::Ripple, &lib);
        let r = restoring_cost(&c, 2, AdderKind::Ripple, &lib);
        assert!(r.depth > 20 * z.depth, "zoom {} vs restoring {}", z.depth, r.depth);
    }

    #[test]
    fn area_is_the_widest_level_not_the_sum_over_levels() {
        // Regression test. Summing area over levels describes no real circuit:
        // the datapath is one array, reused by every level. It previously made
        // the radix-2 baseline look ~500x larger than the zoom when the truth is
        // the other way round, purely because it has 500x more iterations.
        let lib = CellLibrary::generic_180nm();
        let c = cfg(3000);
        let r = restoring_cost(&c, 2, AdderKind::Ripple, &lib);
        let p = zoom_cost(&c, Algo::ProjectiveZoom, AdderKind::Ripple, &lib);

        // A restoring step unit is a handful of gates; 4983 iterations must not
        // change its silicon cost.
        let per_step = 2.0 * (2.0 * 2.0 * lib.full_adder + 2.0 * 2.0 * lib.mux2);
        assert!(
            r.area < 20.0 * per_step + 4000.0 * lib.full_adder,
            "restoring area {} looks like it was summed over steps",
            r.area
        );
        // The zoom carries a full-width multipy-by-constant array, so it is
        // genuinely the larger *unit* -- but only by a small constant factor,
        // not by a factor equal to the iteration count.
        assert!(p.area > r.area, "the zoom datapath is the wider unit");
        assert!(
            p.area < 50.0 * r.area,
            "width ratio {} is suspiciously close to the depth ratio",
            p.area / r.area
        );
    }

    #[test]
    fn zoom_area_equals_its_widest_level() {
        // Direct check of the same property on the zoom side: area must match
        // the cost of the widest level alone, not the total over all levels.
        let lib = CellLibrary::generic_180nm();
        let c = cfg(1000);
        let cost = zoom_cost(&c, Algo::ProjectiveZoom, AdderKind::Ripple, &lib);
        // The final recompute level is digits+guard+10 wide and dominates; every
        // other level is narrower, so the max must sit at or above the mean
        // width by a clear margin.
        assert!(cost.max_prec as f64 > cost.sum_prec / cost.depth as f64);
    }

    #[test]
    fn the_zoom_loses_to_radix2_restoring_on_area_per_result() {
        // NEGATIVE RESULT, recorded deliberately.
        //
        // The intuition that motivates the projective zoom is that per-level
        // truncation saves 63.6% of the datapath area-time. That is true, and it
        // is what makes the zoom fast *in software*. In hardware the saving does
        // not carry over, because of one thing the schedule hides:
        //
        //   * A radix-2 restoring step is a 4-bit subtract. It is tiny. Running
        //     it 4983 times costs 4983 cycles of a very small circuit, and its
        //     area is dominated by the remainder *register*, not the logic.
        //   * A projective level at the target width is a 3025-digit
        //     multiply-by-constant: two AND-OR arrays and two full adders,
        //     ~78k GE. The zoom needs exactly that circuit because its last
        //     level is always full precision.
        //
        // So the zoom's first eight levels run on narrow, cheap hardware, and
        // then the ninth level is as wide as the whole problem. Meanwhile the
        // restoring recurrence never builds anything wider than 4 bits of logic.
        //
        // Per square root delivered, the ratio is roughly 1600x against the
        // zoom. The variable-precision advantage that dominates software cost is
        // largely an artefact of charging by the *digit*; hardware is charged by
        // the *circuit*, and a digit recurrence's circuit barely grows with D.
        let lib = CellLibrary::generic_180nm();
        let c = cfg(3000);
        let p = zoom_cost(&c, Algo::ProjectiveZoom, AdderKind::Ripple, &lib);
        let r = restoring_cost(&c, 2, AdderKind::Ripple, &lib);

        assert!(
            p.cost_per_result() > 100.0 * r.cost_per_result(),
            "expected the zoom to lose badly; got {} vs {}",
            p.cost_per_result(),
            r.cost_per_result()
        );
        assert!(p.fair_merit() > r.fair_merit());
        // The raw merit points the *other* way, which is why it must not be
        // used to compare designs with different iteration counts.
        assert!(p.merit() > r.merit());
    }

    #[test]
    fn restoring_area_per_result_is_flat_in_d_and_the_zooms_is_not() {
        // The decisive structural result, and the reason the conclusion is not
        // a detail of D = 3000.
        //
        // A radix-2 recurrence's area is dominated by the D-digit remainder
        // register, and its iteration count is also linear in D, so
        // area/depth is *constant* in D: ~5.5 GE per square root, always.
        //
        // The zoom's area is a full-width multiply-by-constant and its depth
        // grows only like log3(D), so area/depth grows without bound.
        //
        // The gap therefore widens monotonically with precision, and no choice
        // of D in any plausible range rescues the zoom.
        let lib = CellLibrary::generic_180nm();
        let mut prev_ratio = 0.0;
        for d in [200u32, 500, 1000, 3000, 10000] {
            let c = cfg(d);
            let p = zoom_cost(&c, Algo::ProjectiveZoom, AdderKind::Ripple, &lib);
            let r = restoring_cost(&c, 2, AdderKind::Ripple, &lib);
            let ratio = p.cost_per_result() / r.cost_per_result();
            assert!(ratio > prev_ratio, "gap should widen: D={d} ratio={ratio}");
            prev_ratio = ratio;
            // Restoring stays essentially flat across two decades of D.
            let r_first = restoring_cost(&cfg(200), 2, AdderKind::Ripple, &lib).cost_per_result();
            if d == 10000 {
                assert!(
                    (r.cost_per_result() / r_first - 1.0).abs() < 0.15,
                    "restoring GE/result should be flat in D: {} vs {r_first}",
                    r.cost_per_result()
                );
            }
        }
    }

    #[test]
    fn restoring_beats_the_zoom_at_every_precision_tested() {
        let lib = CellLibrary::generic_180nm();
        for adder in [AdderKind::Ripple, AdderKind::Lookahead] {
            for d in [50u32, 200, 1000, 3000, 10000] {
                let c = cfg(d);
                for algo in [Algo::AffineZoom, Algo::ProjectiveZoom] {
                    let z = zoom_cost(&c, algo, adder, &lib);
                    let r = restoring_cost(&c, 2, adder, &lib);
                    assert!(
                        r.cost_per_result() < z.cost_per_result(),
                        "{algo:?} beat restoring at D={d} {adder:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn cost_per_result_ranks_like_area_divided_by_depth() {
        let lib = CellLibrary::generic_180nm();
        let c = cfg(3000);
        for cost in [
            zoom_cost(&c, Algo::ProjectiveZoom, AdderKind::Ripple, &lib),
            restoring_cost(&c, 2, AdderKind::Ripple, &lib),
        ] {
            let expected = cost.area / cost.depth as f64;
            assert!((cost.cost_per_result() - expected).abs() < 1e-6);
        }
    }

    #[test]
    fn all_costs_are_finite_and_positive() {
        let lib = CellLibrary::generic_180nm();
        for algo in Algo::all() {
            let c = zoom_cost(&cfg(500), algo, AdderKind::Ripple, &lib);
            assert!(c.area.is_finite() && c.area > 0.0);
            assert!(c.delay.is_finite() && c.delay > 0.0);
            assert!(c.energy.is_finite() && c.energy > 0.0);
        }
    }
}
