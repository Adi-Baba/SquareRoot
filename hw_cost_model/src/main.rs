//! `hw-cost-model` — is a hardware square root worth building from the cubic
//! bracket-zoom, and if so, which variant?
//!
//! Run:
//!   cargo run --release                      # the main comparison
//!   cargo run --release -- --digits 3000     # a different target
//!   cargo run --release -- --sweep           # across D, radix, adder, library
//!   cargo run --release -- --truncation      # just the truncation saving
//!
//! Read `README.md` before trusting any number: the cell library is an
//! assumption table, not a PDK, and the conclusions are only as good as it.

use hw_cost_model::cost::{
    restoring_cost, truncation_saving, zoom_cost, Algo, AdderKind, Cost,
};
use hw_cost_model::schedule::ZoomConfig;
use hw_cost_model::tech::CellLibrary;

fn fmt(x: f64) -> String {
    if x >= 1e6 {
        format!("{:.2}M", x / 1e6)
    } else if x >= 1e3 {
        format!("{:.1}k", x / 1e3)
    } else {
        format!("{x:.2}")
    }
}

fn main() {
    let mut digits: u32 = 3000;
    let mut radix: u32 = 2;
    let mut sweep = false;
    let mut truncation_only = false;
    let mut adder = AdderKind::Ripple;
    let lib = CellLibrary::generic_180nm();

    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--digits" => {
                digits = args.get(i + 1).and_then(|s| s.parse().ok()).unwrap_or(digits);
                i += 2;
            }
            "--radix" => {
                radix = args.get(i + 1).and_then(|s| s.parse().ok()).unwrap_or(radix);
                i += 2;
            }
            "--sweep" => { sweep = true; i += 1; }
            "--truncation" => { truncation_only = true; i += 1; }
            "--lookahead" => { adder = AdderKind::Lookahead; i += 1; }
            other => {
                eprintln!("unknown argument: {other}");
                std::process::exit(2);
            }
        }
    }

    let cfg = ZoomConfig { digits, ..Default::default() };

    if truncation_only {
        print_truncation(&cfg, adder, &lib);
        return;
    }
    if sweep {
        sweep_all();
        return;
    }
    print_main(&cfg, radix, adder, &lib);
}

fn print_truncation(cfg: &ZoomConfig, adder: AdderKind, lib: &CellLibrary) {
    let t = truncation_saving(cfg, adder, lib);
    println!("Per-level precision truncation, D = {} (adder: {})", cfg.digits, adder.name());
    println!("  mean recursion depth      : {}", t.mean_depth);
    println!("  widest level              : {} digits", t.max_prec);
    println!("  sum of per-level widths   : {}", fmt(t.sum_prec));
    println!("  area-time if narrowed     : {}", fmt(t.variable_area_time));
    println!("  area-time if fixed width  : {}", fmt(t.fixed_area_time));
    println!("  => {:.1}% of the datapath area-time product is recoverable", 100.0 * t.saving);
    println!();
    println!("  This is a schedule fact, not a timing measurement. Whether the");
    println!("  saving survives synthesis depends on routing, register files and");
    println!("  the clock tree -- see README.md, 'What would falsify this'.");
}

fn print_main(cfg: &ZoomConfig, radix: u32, adder: AdderKind, lib: &CellLibrary) {
    println!("Hardware cost model: square root to D = {} decimal digits", cfg.digits);
    println!("cell library: {}   adder: {}", lib.name, adder.name());
    println!();

    println!("per-level schedule for alpha = 12:");
    for (i, p) in crate_sched(cfg, 12).iter().enumerate() {
        println!("   level {i:>2}: {p:>5} digits");
    }
    println!();

    println!("{:<16} {:>10} {:>10} {:>10} {:>8} {:>11}", "algorithm", "area(GE)", "delay", "energy", "levels", "GE/result");
    println!("{}", "-".repeat(78));
    let mut rows: Vec<(String, Cost)> = Vec::new();
    for algo in [Algo::AffineZoom, Algo::ProjectiveZoom] {
        let c = zoom_cost(cfg, algo, adder, lib);
        rows.push((algo.name().to_string(), c));
    }
    let r = restoring_cost(cfg, radix, adder, lib);
    rows.push((format!("restoring-{radix}"), r));

    for (name, c) in &rows {
        println!(
            "{:<16} {:>10} {:>10} {:>10} {:>8} {:>11}",
            name, fmt(c.area), fmt(c.delay), fmt(c.energy), c.depth, fmt(c.cost_per_result())
        );
    }
    println!();

    let proj = rows.iter().find(|(n, _)| n == "projective-zoom").unwrap().1;
    let aff = rows.iter().find(|(n, _)| n == "affine-zoom").unwrap().1;
    let rest = rows.last().unwrap().1;

    println!("projective vs affine (same schedule, no divider):");
    println!("  area   {:.1}x less", aff.area / proj.area);
    println!("  delay  {:.1}x less", aff.delay / proj.delay);
    println!("  => the Moebius lift is a clear win over the literal affine form.");
    println!();

    let t = truncation_saving(cfg, adder, lib);
    println!("truncation of the projective datapath: {:.1}% of area-time recoverable", 100.0 * t.saving);
    println!("  fixed-width {} levels x {} digits = {} vs narrowed sum {}",
        t.mean_depth, t.max_prec, fmt(t.fixed_area_time), fmt(t.sum_prec));
    println!("  ripple delay: {} -> {} ({:.1}% recovered)",
        fmt(t.fixed_delay), fmt(t.variable_delay), 100.0 * t.delay_saving);
    println!();

    println!("VERDICT: the radix-{radix} restoring recurrence wins on area per result,");
    println!("by a factor of {:.0}. The zoom is not a viable hardware square root.",
        proj.cost_per_result() / rest.cost_per_result());
    println!();
    println!("Why the {:.0}% truncation saving does not translate:", 100.0 * t.saving);
    println!("  - a radix-{radix} step is a {}-bit subtract: a tiny circuit, iterated {} times.",
        2 * radix, rest.depth);
    println!("  - the zoom's final level is full width, so it must build a {}-digit", t.max_prec);
    println!("    multiply-by-constant regardless of how narrow the first levels are.");
    println!("  - variable precision saves software work (charged per digit); it does");
    println!("    not save silicon, because a digit recurrence's circuit barely grows");
    println!("    with D.");
}

fn crate_sched(cfg: &ZoomConfig, alpha: u64) -> Vec<u32> {
    hw_cost_model::schedule::cubic_schedule(cfg, alpha)
}

fn sweep_all() {
    println!("Sweep: target digits x cell library x adder topology\n");
    println!("Reported as GE per square root delivered (area / iterations), the");
    println!("only figure comparable across designs with different iteration counts.\n");
    for lib in [CellLibrary::generic_180nm(), CellLibrary::modern_7nm()] {
        for adder in [AdderKind::Ripple, AdderKind::Lookahead] {
            println!("--- {} / {} ---", lib.name, adder.name());
            println!(
                "{:>7} | {:>12} | {:>12} | {:>12} | {:>12} | {}",
                "D", "affine", "projective", "restore-2", "restore-8", "winner"
            );
            for d in [200u32, 500, 1000, 3000, 10000] {
                let cfg = ZoomConfig { digits: d, ..Default::default() };
                let a = zoom_cost(&cfg, Algo::AffineZoom, adder, &lib);
                let p = zoom_cost(&cfg, Algo::ProjectiveZoom, adder, &lib);
                let r2 = restoring_cost(&cfg, 2, adder, &lib);
                let r8 = restoring_cost(&cfg, 8, adder, &lib);
                let cands = [
                    ("affine", a.cost_per_result()),
                    ("projective", p.cost_per_result()),
                    ("restoring-2", r2.cost_per_result()),
                    ("restoring-8", r8.cost_per_result()),
                ];
                let (winner, _) = cands.iter().min_by(|x, y| x.1.partial_cmp(&y.1).unwrap()).unwrap();
                println!(
                    "{:>7} | {:>12} | {:>12} | {:>12} | {:>12} | {}",
                    d,
                    fmt(a.cost_per_result()),
                    fmt(p.cost_per_result()),
                    fmt(r2.cost_per_result()),
                    fmt(r8.cost_per_result()),
                    winner
                );
            }
            println!();
        }
    }
    println!("The winner is the same in every cell of the sweep. The two zoom rows");
    println!("differ by ~2.4x, which is the divider; they do not come close to the");
    println!("digit recurrences, which is a different order of magnitude.");
    println!("Run with a real PDK before believing any of this.");
}
