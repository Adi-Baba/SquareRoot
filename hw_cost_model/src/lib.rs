//! A hardware cost model for square root by cubic bracket-zoom.
//!
//! The question this crate answers: **if the bracket-zoom of Papers 1 and 2
//! were built in hardware, would it beat the radix-`2^b` restoring digit
//! recurrence that shipped SQR units actually implement?** Three things have to
//! be true for that to be worth asking, and each is modelled separately:
//!
//! 1. The per-level operation must be cheap. In affine form it is *not* — it
//!    needs two divisions by a running value. The projective lift
//!    `(X,Y) -> (aX + alpha*Y, X + aY)` is two multiplies and two adds, with no
//!    division, and is exact, so the algebra is not a trick that breaks
//!    bit-for-bit.
//! 2. The per-level precision must be *truncatable*, otherwise the narrow early
//!    levels save nothing. Paper1 section 10 gives the schedule; the model
//!    reports how much datapath area-time that recovers.
//! 3. The cell assumptions must be honest. They are in `tech.rs`, isolated, and
//!    swappable for a real PDK.
//!
//! # Reading the output
//!
//! `merit = area * delay * energy` is a *proxy*. The parts are always printed
//! alongside it, because the product hides which term is doing the work, and
//! because a multiplier that is 3x worse on delay can be 10x better on area
//! with no clocking pressure to care. See `README.md` for what would falsify
//! the conclusions.
//!
//! # Scope
//!
//! This models arithmetic *area*, *depth* and *switching energy* for the
//! datapath only. It does not model memory for storing the radicand, the
//! register file between levels, routing congestion, or clock tree power —
//! and for a 3000-digit datapath those omissions are not small.

pub mod cost;
pub mod schedule;
pub mod tech;

pub use cost::{Algo, AdderKind, Cost};
pub use schedule::ZoomConfig;
pub use tech::CellLibrary;
