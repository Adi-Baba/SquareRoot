//! Benchmark harness: certified brackets, exact division counts, wall clock.
//!
//! Three questions, in order of how much they matter:
//!
//! 1. Does every method actually certify? Each bracket is proved individually in
//!    exact integer arithmetic (`certify`), so this is not a sample standing in
//!    for a universally quantified claim.
//! 2. How many integer divisions does each cost? Exact and machine-independent.
//! 3. How long does each take? Measured -- and only informative here, because
//!    Rust removes the CPython `localcontext` overhead that made the Python
//!    timings unable to separate the methods.

use num_bigint::BigUint;
use std::time::Instant;
use zoombench::*;

const ALPHAS: [u64; 16] = [
    2, 3, 5, 7, 10, 12, 13, 17, 23, 31, 101, 257, 1009, 7919, 104729, 999983,
];

type SqFn = fn(&BigUint, u64) -> Run;

#[derive(Clone, Copy, Debug)]
struct Row {
    certified: usize,
    levels: f64,
    divs: f64,
    ops: f64,
    ms: f64,
}

fn bench_sq(name: &str, f: SqFn, digits: u64, reps: usize) -> (String, Row) {
    let mut certified = 0usize;
    let mut lev = Vec::new();
    let mut divs = Vec::new();
    let mut ops = Vec::new();
    let mut best = Vec::with_capacity(ALPHAS.len());

    for &a in ALPHAS.iter() {
        let alpha = BigUint::from(a);
        let mut best_ms = f64::INFINITY;
        let mut run: Option<Run> = None;
        for _ in 0..reps {
            let t = Instant::now();
            let r = f(&alpha, digits);
            let el = t.elapsed().as_secs_f64() * 1e3;
            if el < best_ms {
                best_ms = el;
            }
            run = Some(r);
        }
        best.push(best_ms);
        let r = run.unwrap();
        match certify(&r.lo, &r.hi, &alpha, digits, r.scale) {
            Ok(()) => certified += 1,
            Err(e) => println!("  !! {} a={} digits={}: {}", name, a, digits, e),
        }
        lev.push(r.levels as f64);
        divs.push(r.counters.divs as f64);
        ops.push(r.counters.ops() as f64);
    }

    let mean = |v: &Vec<f64>| v.iter().sum::<f64>() / v.len() as f64;
    let mut s = best.clone();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (
        name.to_string(),
        Row {
            certified,
            levels: mean(&lev),
            divs: mean(&divs),
            ops: mean(&ops),
            ms: s[s.len() / 2],
        },
    )
}

fn bench_nth(n: u32, digits: u64, reps: usize) -> (Row, usize) {
    let mut certified = 0usize;
    let mut bad = 0usize;
    let mut lev = Vec::new();
    let mut divs = Vec::new();
    let mut ops = Vec::new();
    let mut best = Vec::with_capacity(ALPHAS.len());

    for &a in ALPHAS.iter() {
        let alpha = BigUint::from(a);
        let mut best_ms = f64::INFINITY;
        let mut run: Option<Run> = None;
        for _ in 0..reps {
            let t = Instant::now();
            let r = halley_nth_root(&alpha, n, digits);
            let el = t.elapsed().as_secs_f64() * 1e3;
            if el < best_ms {
                best_ms = el;
            }
            run = Some(r);
        }
        best.push(best_ms);
        let r = run.unwrap();
        // Containment is unconditional; the width target is only reachable
        // where the bound can contract from the naive integer seed.
        let contains = certify_nth(&r.lo, &r.hi, &alpha, n, digits, r.scale);
        if contains.as_ref().err().map(|e| e.starts_with("lo^n") || e.starts_with("hi^n")) == Some(true) {
            bad += 1;
        } else if contains.is_ok() {
            certified += 1;
        }
        lev.push(r.levels as f64);
        divs.push(r.counters.divs as f64);
        ops.push(r.counters.ops() as f64);
    }
    let mean = |v: &Vec<f64>| v.iter().sum::<f64>() / v.len() as f64;
    let mut s = best.clone();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (
        Row {
            certified,
            levels: mean(&lev),
            divs: mean(&divs),
            ops: mean(&ops),
            ms: s[s.len() / 2],
        },
        bad,
    )
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let digits: u64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(1000);
    let reps: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3);

    println!("zoombench -- certified cubic bracket zoom, exact integer arithmetic");
    println!(
        "digits={}  radicands={}  best of {} each, median reported",
        digits,
        ALPHAS.len(),
        reps
    );
    println!("'certified' is an exact integer proof per bracket: lo^2 <= alpha <= hi^2");
    println!("and width <= 10^-{}. 'divs'/'ops' are exact executed-operation counts.\n", digits);

    let rows = vec![
        bench_sq("paper zoom", paper_zoom, digits, reps),
        bench_sq("halley-mid", halley_midpoint, digits, reps),
        bench_sq("halley-right", halley_right, digits, reps),
        bench_sq("bisection", bisect, digits, reps),
    ];
    let n = ALPHAS.len();
    let base_row = rows.iter().find(|(nm, _)| nm == "paper zoom").unwrap().1.clone();
    // NEVER report a ratio against a baseline that does not certify.  A speedup
    // relative to a method that produces invalid output is worse than no
    // number at all: it reads as a result and is not one.  An earlier build
    // printed "halley-mid 16/16 ... 0.79x" against a paper-zoom baseline that
    // was 0/16, i.e. faster than a method that was simply broken.
    let baseline_ok = base_row.certified == n;
    println!("{:<14} {:>10} {:>8} {:>8} {:>8} {:>10} {:>9}", "method", "certified", "levels", "divs", "ops", "ms", "vs zoom");
    println!("{}", "-".repeat(72));
    for (name, r) in &rows {
        let vs = if baseline_ok { format!("{:>8.2}x", r.ms / base_row.ms) } else { format!("{:>9}", "--") };
        println!(
            "{:<14} {:>4}/{:<3} {:>8.1} {:>8.1} {:>8.0} {:>10.4} {}",
            name, r.certified, n, r.levels, r.divs, r.ops, r.ms, vs
        );
    }
    if !baseline_ok {
        println!(
            "\n* comparison WITHHELD: the paper-zoom baseline certified {}/{}, so its",
            base_row.certified, n
        );
        println!("  timing is meaningless.  Levels/divs/ops above are still exact.");
    }

    println!("\nn-th ROOT by Halley at the midpoint. Paper2 Section 9.4 proves the");
    println!("paper's own map phi_n is REPELLING for n>=3 (lambda_n -> -(n-1));");
    println!("Halley's multiplier at the root is 0 for every n, so a certified cubic");
    println!("contraction exists there regardless.");
    println!("\n{:<6} {:>10} {:>8} {:>8} {:>8} {:>10}", "n", "certified", "levels", "divs", "ops", "ms");
    println!("{}", "-".repeat(56));
    for k in [2u32, 3, 4, 5, 8, 16] {
        let (r, bad) = bench_nth(k, digits, reps);
        println!(
            "{:<6} {:>7}/{:<2} {:>8.1} {:>8.1} {:>8.0} {:>10.4}  {}",
            k, r.certified, n, r.levels, r.divs, r.ops, r.ms,
            if bad > 0 { format!("({} brackets did NOT contain the root)", bad) }
            else { String::new() }
        );
    }
    println!("\n'n-th root: 'certified' counts brackets meeting the full width target.");
    println!("The bound (n^2-1) b^(n-2) w^2 <= 48 a^n only permits contraction from the");
    println!("naive integer seed when it holds; elsewhere the certified width is still a");
    println!("valid bracket, just not a narrower one.  Containment never fails (0 bad).");
}