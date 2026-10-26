# hw-cost-model

An area–delay–energy cost model asking one question:

> The cubic bracket-zoom of Papers 1 and 2 is fast **in software**. Would it be
> a good **hardware** square root?

**Answer: no.** The radix-2 restoring digit recurrence — what shipped SQR
units actually implement — wins by roughly three orders of magnitude on area
per square root delivered, and the gap *widens* with precision. The reasoning
and the numbers are below, because a negative result is only worth as much as
the ability to check it.

## What was modelled

| | affine zoom | projective zoom | radix-`2^b` restoring |
|---|---|---|---|
| per-level operation | 2 divisions + 2 mults | 2 mults + 2 adds | one `b`-bit subtract |
| levels at D=3000 | 9 | 9 | 4983 |
| datapath width | 3025 digits | 3025 digits | 4 bits of logic |

The **projective** form is the Moebius lift from Paper 2: the map
`est = a + p/(a + m)` becomes the pair update `(X,Y) -> (aX + alpha*Y, X + aY)`,
which is exact, division-free, and worth doing — it is 2.4x smaller and 405x
faster in delay than the literal affine transcription. That part of the idea
works.

The **truncation schedule** also works, in the sense Paper 1 claims: per-level
precision at D=3000 is `33, 42, 63, 129, 327, 915, 2685, 3015, 3025`, summing
to 10234 against a fixed-width `9 x 3025 = 27225`, recovering **63.6%** of the
datapath area-time product. This is verified against the Python reference by
`cubic_schedule_matches_paper1_reference`.

Neither of those is the question that matters, and this is the point of the
exercise.

## Why the truncation saving does not become a hardware saving

```
algorithm          area(GE)      delay     energy   levels   GE/result
affine-zoom          187.6k    143.68M     484.8k        9       20.8k
projective-zoom       78.7k     354.9k     239.8k        9        8.7k
restoring-2           27.2k      89.7k     129.2k     4983        5.46
```

Three structural facts, none of them visible from the schedule:

1. **Area is the widest level, never the sum.** The datapath is one physical
   array reused by every level. Summing area over levels — the obvious first
   implementation, and one this crate originally shipped and tested wrong —
   makes a 4983-step recurrence look 500x larger than a 9-level one, when in
   fact it is a single small step unit iterated. See
   `area_is_the_widest_level_not_the_sum_over_levels`.

2. **The zoom's last level is full width.** Eight narrow levels and then one
   3025-digit multiply-by-constant. The saving applies to levels 0–7; level 8
   sets the area floor regardless. A digit recurrence never faces this: its
   logic is 4 bits wide no matter how large D is, and its area is dominated by a
   register that grows only linearly.

3. **Therefore the iteration counts do not cancel.** The restoring recurrence's
   area is dominated by a `D`-digit register and its iteration count is linear
   in `D`, so GE per result is *flat* at ~5.5 across two decades of precision.
   The zoom's area is full-width and its depth grows only like `log3(D)`, so GE
   per result grows without bound.

| D | affine | projective | restoring-2 |
|---|---|---|---|
| 200 | 2.3k | 975 | 6.10 |
| 1000 | 7.9k | 3.3k | 5.56 |
| 3000 | 20.8k | 8.7k | 5.46 |
| 10000 | 62.2k | 26.1k | 5.43 |

## The general lesson

Variable precision is a **software** optimisation. Software cost is charged per
digit, so computing 33 digits instead of 3000 is a ~90x saving. Hardware cost
is charged per circuit, and a digit recurrence's circuit barely grows with D.
The two cost models disagree, and for this particular algorithm the hardware
model is the one that decides the question.

This is the same shape of error as the ZKP proposal rejected in Paper 2's
sibling analysis: a property that holds in the arithmetic and fails to
transfer to the setting it is being used to justify.

## What would falsify this

The model omits, and would need before the conclusion could be revisited:

- **Register file and memory.** The zoom stores a growing `D`-digit bracket
  between levels; the recurrence stores one partial root. This *widens* the gap
  further, so it does not rescue the zoom.
- **Routing and clock tree.** At 3025 digits (~10k bits) the projective
  multiply-by-constant is routing-bound, not gate-bound, so its real area is
  worse than the GE estimate. Again, against the zoom.
- **A real PDK.** The cell library is an assumption table in `tech.rs`, not a
  foundry library. Sensitivity to it is swept and the verdict is unchanged, but
  the absolute numbers are not synthesis numbers.
- **Multi-precision parallelism.** If the workload is *many* square roots and
  the datapath can be time-shared across lanes, the fixed-width recurrence is
  far better at filling a pipeline. Against the zoom.
- **Fused or pipelined digit recurrences.** A pipelined radix-`2^b` unit is the
  standard high-throughput answer and is not separately modelled here.

The only credible route to a favourable result would be a way to avoid
building the full-width multiply, e.g. a Newton-style correction that reuses
narrow levels, or a digit-recurrence *transcription* of the zoom. The latter is
what the projective form hints at, and it is the obvious next thing to try.

## Running it

```
cargo test                        # 22 tests
cargo run --release               # the main comparison
cargo run --release -- --sweep    # D x cell library x adder topology
cargo run --release -- --truncation
cargo run --release -- --digits 10000 --radix 4
```

## Layout

| file | contents |
|---|---|
| `src/tech.rs` | cell library assumptions, adder delay models |
| `src/schedule.rs` | cubic zoom and restoring schedules |
| `src/cost.rs` | the cost model, level delays, figures of merit |
| `src/main.rs` | CLI and reporting |

`src/schedule.rs` reproduces `Src/exp_precision.py`'s schedule loop, including
the drift of the anchor `a` from the integer seed toward `sqrt(alpha)`. That
detail matters: holding `a` fixed mis-schedules every level from the third
onward, because `eps = w^3/(27a^2)` depends on it.
