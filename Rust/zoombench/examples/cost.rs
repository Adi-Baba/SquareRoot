//! Is division the bottleneck?  Measure it directly in Rust.
//!
//! Two experiments:
//!  A. the price of one BigUint division vs multiplication vs addition, at the
//!     operand sizes the zoom actually uses (d digits each);
//!  B. the share of Halley's certified run spent in division, by timing the same
//!     computation with the divisions replaced by multiplications.

use num_bigint::BigUint;
use std::time::Instant;
use zoombench::*;

fn ms(f: impl Fn() -> ()) -> f64 {
    let t = Instant::now();
    f();
    t.elapsed().as_secs_f64() * 1e3
}

fn main() {
    println!("A. cost of ONE BigUint operation, operands of D digits each\n");
    println!("{:>7} {:>10} {:>10} {:>10} {:>12}", "digits", "mul ms", "add ms", "div ms", "div/mul");
    for d in [100u64, 500, 1000, 3000, 10000] {
        let a = pow10(d);          // ~d+1 digits
        let b = pow10(d / 2 + 3);  // different size, avoids trivial cancellation
        let n = 200;
        let tm = ms(|| { for _ in 0..n { let _ = &a * &b; } }) / n as f64 * 1e3;
        let ta = ms(|| { for _ in 0..n { let _ = &a + &b; } }) / n as f64 * 1e3;
        let td = ms(|| { for _ in 0..n { let _ = &a / &b; } }) / n as f64 * 1e3;
        println!("{:>7} {:>10.4} {:>10.4} {:>10.4} {:>12.1}", d, tm, ta, td, td / tm);
    }

    println!("\nB. share of a real certified run spent in division");
    println!("   Halley-mid at 1000 digits, 16 radicands:\n");
    let alphas: [u64; 16] = [2, 3, 5, 7, 10, 12, 13, 17, 23, 31, 101, 257, 1009, 7919, 104729, 999983];
    let mut levels = 0u64;
    let mut divs = 0u64;
    let mut muls = 0u64;
    let mut adds = 0u64;
    for a in alphas {
        let r = halley_midpoint(&BigUint::from(a), 1000);
        levels += r.levels;
        divs += r.counters.divs;
        muls += r.counters.muls;
        adds += r.counters.adds;
    }
    println!("   totals over 16 radicands: {} levels, {} divisions, {} muls, {} adds", levels, divs, muls, adds);
    println!("   divisions are {:.0}% of all counted operations\n", divs as f64 / (divs + muls + adds) as f64 * 100.0);

    // Price the two components with the measured per-op costs above at 1000 digits.
    let (tm, td) = {
        let a = pow10(1000);
        let b = pow10(503);
        let n = 200;
        let m = ms(|| { for _ in 0..n { let _ = &a * &b; } }) / n as f64 * 1e3;
        let d = ms(|| { for _ in 0..n { let _ = &a / &b; } }) / n as f64 * 1e3;
        (m, d)
    };
    println!("   at 1000 digits: div = {:.4} ms, mul = {:.4} ms", td, tm);
    println!("   predicted cost of the divisions alone : {:.3} ms", divs as f64 * td);
    println!("   predicted cost of the multiplications : {:.3} ms", muls as f64 * tm);
    println!("   => divisions are {:.0}% of the arithmetic, so removing them",
        divs as f64 * td / (divs as f64 * td + muls as f64 * tm) * 100.0);
    println!("      could save at most that share of the ARITHMETIC, not of the run.");
}
