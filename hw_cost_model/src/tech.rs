//! Technology assumptions: what a gate costs.
//!
//! Every number here is an assumption, and the whole cost model is only as
//! trustworthy as this file. They are collected in one place so that a reader
//! can substitute their own PDK and re-run.
//!
//! Units. Area is in arbitrary "gate equivalents" (GE): one GE is the area of a
//! single 2-input NAND. Delay is in units of the fanout-of-four delay `t_fo4`.
//! Energy is a switching-energy proxy: dynamic energy per evaluation is
//! approximated as proportional to toggled capacitance, so we charge activity
//! factors per node type rather than joules.

/// A standard-cell library, in gate equivalents and `t_fo4`.
#[derive(Debug, Clone, Copy)]
pub struct CellLibrary {
    pub name: &'static str,
    /// 2-input NAND / inverter.
    pub inv: f64,
    pub nand2: f64,
    /// Full adder, the workhorse of any adder.
    pub full_adder: f64,
    /// Half adder.
    pub half_adder: f64,
    /// 2:1 multiplexer, one bit.
    pub mux2: f64,
    /// AND-OR pair (an alternative adder element; comparable to a FA).
    pub and_or: f64,
    /// Static CMOS power, relative units per node per evaluation.
    pub switch_power_inv: f64,
    pub switch_power_fa: f64,
    pub switch_power_mux: f64,
    /// Leakage per node (relative), charged per cycle whether or not it toggles.
    pub leakage: f64,
}

impl CellLibrary {
    /// Generic 180 nm style numbers. Deliberately conservative and old: a
    /// square-root datapath is dominated by adders and muxes, whose ratios
    /// barely move between nodes.
    pub fn generic_180nm() -> Self {
        CellLibrary {
            name: "generic-180nm",
            inv: 1.0,
            nand2: 1.0,
            full_adder: 9.0,
            half_adder: 5.0,
            mux2: 3.0,
            and_or: 4.0,
            switch_power_inv: 1.0,
            switch_power_fa: 12.0,
            switch_power_mux: 3.0,
            leakage: 0.02,
        }
    }

    /// A modern node, roughly 3x denser and 2.5x faster than the 180nm entry.
    pub fn modern_7nm() -> Self {
        CellLibrary {
            name: "modern-7nm",
            inv: 0.3,
            nand2: 0.3,
            full_adder: 2.6,
            half_adder: 1.5,
            mux2: 0.9,
            and_or: 1.2,
            switch_power_inv: 0.3,
            switch_power_fa: 3.0,
            switch_power_mux: 0.9,
            leakage: 0.10,
        }
    }
}

/// Delay of a ripple-carry `width`-bit adder, in `t_fo4`.
///
/// Ripple carry, `t = width * fa_delay`. The FA delay itself is estimated
/// from its area: on a standard-cell library delay tracks area roughly
/// linearly in the first order, and we use `fa_delay = full_adder / 3`.
pub fn ripple_delay(width: usize, lib: &CellLibrary) -> f64 {
    width as f64 * (lib.full_adder / 3.0)
}

/// A carry-select/lookahead adder: logarithmic in width, with a larger
/// constant. Used for the sensitivity sweep, since the whole point of
/// truncating a datapath is that it buys *depth* as well as area.
pub fn fast_delay(width: usize, lib: &CellLibrary) -> f64 {
    let w = width.max(2) as f64;
    let log_w = w.log2();
    // O(log w) stages of roughly constant cost.
    log_w * (lib.full_adder / 3.0) + 4.0 * lib.full_adder / 3.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn delay_is_monotone_in_width() {
        let lib = CellLibrary::generic_180nm();
        let mut prev = 0.0;
        for w in [8usize, 32, 128, 512, 2048, 8192] {
            let d = ripple_delay(w, &lib);
            assert!(d > prev, "ripple delay must grow with width");
            prev = d;
        }
        let mut prevf = 0.0;
        for w in [8usize, 32, 128, 512, 2048, 8192] {
            let d = fast_delay(w, &lib);
            assert!(d > prevf, "fast delay must grow with width");
            prevf = d;
        }
    }

    #[test]
    fn lookahead_beats_ripple_at_large_width() {
        let lib = CellLibrary::generic_180nm();
        for w in [512usize, 2048, 8192, 30000] {
            assert!(
                fast_delay(w, &lib) < ripple_delay(w, &lib),
                "logarithmic adder must beat ripple at width {w}"
            );
        }
    }

    #[test]
    fn modern_library_is_smaller_and_faster() {
        let a = CellLibrary::generic_180nm();
        let b = CellLibrary::modern_7nm();
        assert!(b.full_adder < a.full_adder);
        assert!(ripple_delay(2048, &b) < ripple_delay(2048, &a));
    }
}
