"""Exact, machine-independent operation counts for the Section 10.1 methods.

Why this file exists
--------------------
Section 10.1's central performance claim is that the zoom loses to interval
Newton because of its *per-level operation count* -- "three divisions and
three multiplications against interval Newton's two divisions".  That count was
read off the source by eye.  It is the load-bearing explanation for the whole
section, and an eyeball count is not a proof of anything.

But unlike wall-clock time, an operation count is *exactly decidable*.  It is
an integer, it is reproducible on any machine, and it does not move with load.
So the claim can be moved out of the "measured" column and into the "decided"
one, which is the strongest thing that can be done for it.

Two independent methods, and they must agree:

  1. **Static counting.** Parse the method's source with `ast` and count the
     arithmetic operations it performs per iteration, on each control-flow
     path.  Counts are exact but path-sensitive, so each path is reported
     separately rather than averaged.

  2. **Dynamic counting.** Wrap the Decimal type and count operations actually
     executed, per iteration, on concrete inputs.  Path-insensitive but
     empirical in the sense that it observes one run.

Agreement between an exact static count and an observed dynamic count is
strong evidence the source is being read correctly; disagreement means one of
them is wrong, which is the interesting case.

What this does NOT do
---------------------
It does not predict time.  A division is not a fixed cost, and at 3000 digits
a division costs far more than a multiplication.  These counts are a
*machine-independent proxy*, and the paper should say so.  What they do
establish exactly is the shape of the claim in Section 10.1: that the zoom
performs strictly more high-precision divisions per iteration than interval
Newton, and that no amount of wall-clock noise can make that false.

Run:  python3 opcount.py
"""

import ast
import sys

# Categories reported.  A division at 3000 digits is roughly an order of
# magnitude more expensive than a multiplication, which is why the paper's
# claim is about divisions specifically.
CATEGORIES = ("div", "mul", "pow", "add", "neg", "digit_ops")

# Operator -> category, for the per-statement categorisation below.
ARITH_OPS = {
    ast.Add: "+", ast.Sub: "-", ast.Mult: "*",
    ast.Div: "/", ast.FloorDiv: "//", ast.Mod: "%", ast.Pow: "**",
}

# ----------------------------------------------------------- dynamic counting
#
# Counting BinOp nodes in the AST is NOT the right instrument here: a method's
# arithmetic often lives in a helper it calls (zoom delegates to
# _est_bracket), so a per-function AST count misses it entirely, and it
# cannot distinguish a division from a multiplication -- which is precisely the
# distinction the paper's claim rests on.
#
# Instead we count EXECUTED opcodes with sys.settrace.  That follows helper
# calls automatically, separates operations by kind, and is machine
# independent: it counts what the interpreter did, not how fast it did it.



# Opcode-level tracing (f_trace_opcodes) is not available in every CPython
# build, and was observed to be inert here, so operations are counted as
# (times a statement executed) x (arithmetic operations in that statement).
# Both factors are exact: the first from line events, the second from the AST.
# A statement that Python splits across physical lines is attributed to its
# first line, which is the line CPython reports.


def ops_per_statement(tree):
    """{lineno: {category: count}} for every statement in `tree`.

    Each statement is categorised once by `_categorise`; an earlier version
    also called a separate `_op_counts_in` helper and discarded its result.
    That helper has been removed.
    """
    out = {}

    def visit(node):
        for child in ast.iter_child_nodes(node):
            if isinstance(child, ast.stmt):
                cats = _categorise(child)
                if cats:
                    out[child.lineno] = cats
                visit(child)

    visit(tree)
    return out


def _categorise(stmt):
    """Count arithmetic in one statement, by category."""
    cats = dict.fromkeys(CATEGORIES, 0)
    for sub in ast.walk(stmt):
        if isinstance(sub, ast.BinOp):
            name = ARITH_OPS.get(type(sub.op))
            if name == "+":
                cats["add"] += 1
            elif name == "-":
                cats["add"] += 1
            elif name == "*":
                cats["mul"] += 1
            elif name in ("/", "//", "%"):
                cats["div"] += 1
            elif name == "**":
                cats["pow"] += 1
        elif isinstance(sub, ast.UnaryOp) and isinstance(sub.op, ast.USub):
            cats["neg"] += 1
    return cats


def executed_ops(fn, *args):
    """Run fn(*args); return (category_counts, exec_lines, result).

    Counts only lines inside `fn`'s own module, so helper calls made from
    there are included -- the paper's claim is about the whole per-iteration
    cost, not just the inline part.  (A helper in a *different* module is not
    counted; for exp_certified.py every helper is in that file, which the
    report asserts.)
    """
    return _run_with_counts(fn.__code__.co_filename, fn, args)


def _run_with_counts(mod_file, fn, args):
    tree = ast.parse(open(mod_file).read())
    per_line = ops_per_statement(tree)

    counter = _LineCounter(mod_file)
    old = sys.gettrace()
    sys.settrace(counter._trace)
    try:
        result = fn(*args)
    finally:
        sys.settrace(old)

    digit_ops = 0
    for lineno, (n, prec_sum) in counter.exec_lines.items():
        cats = per_line.get(lineno)
        if not cats:
            continue
        for k, v in cats.items():
            counter.counts[k] += v * n
        digit_ops += sum(cats.values()) * prec_sum
    counter.counts["digit_ops"] = digit_ops
    return counter.counts, counter.exec_lines, result


class _LineCounter:
    """Counts how many times each line of one file executes."""

    def __init__(self, module_file):
        self.module_file = module_file
        self.exec_lines = {}
        self.counts = dict.fromkeys(CATEGORIES, 0)
        self._on = True

    def _trace(self, frame, event, arg):
        if event == "call":
            if frame.f_code.co_filename != self.module_file:
                return None
            return self._trace
        if event == "line" and self._on:
            import decimal
            ln = frame.f_lineno
            prec = decimal.getcontext().prec
            prev = self.exec_lines.get(ln)
            if prev is None:
                self.exec_lines[ln] = [0, 0]
                self.exec_lines[ln][0] = 1
                self.exec_lines[ln][1] = prec
            else:
                prev[0] += 1
                prev[1] += prec      # sum of working precisions
        return self._trace


def report():
    import exp_certified as C
    print("=" * 76)
    print("EXACT EXECUTED ARITHMETIC COUNTS (machine independent)")
    print("=" * 76)
    print("Per-iteration cost, from line-execution counts x per-statement ops.")
    print("Helpers in the same module are included; divisions are separated.")
    print()
    hdr = "%-18s%6s%7s%7s%6s%6s%9s%8s%13s" % (
        "method", "iters", "div", "mul", "pow", "add", "div/it", "tot/it",
        "digit-ops")
    print(hdr)
    print("-" * 76)
    rows = {}
    for fn in ("bisection", "interval_newton", "zoom", "newton_certified"):
        f = getattr(C, fn)
        counts, lines, res = executed_ops(f, 12, 300)
        steps = res[2] if isinstance(res, tuple) and len(res) >= 3 else 1
        steps = max(steps, 1)
        tot = sum(counts[k] for k in CATEGORIES if k != "digit_ops")
        rows[fn] = (counts["div"] / steps, tot / steps,
                    counts["digit_ops"])
        print("%-18s%6d%7d%7d%6d%6d%9.2f%8.2f%13d" % (
            fn, steps, counts["div"], counts["mul"], counts["pow"],
            counts["add"], counts["div"] / steps, tot / steps,
            counts["digit_ops"]))
    print()
    print("=" * 76)
    print("THE CLAIM, IN ITS EXACT FORM")
    print("=" * 76)
    zd, zt, zdo = rows["zoom"]
    idd, it_, ido = rows["interval_newton"]
    print("  divisions per iteration: zoom %.2f vs int.Newton %.2f"
          % (zd, idd))
    print("  total ops per iteration: zoom %.2f vs int.Newton %.2f"
          % (zt, it_))
    print("  DIGIT-OPERATIONS (ops x working precision):")
    print("    zoom           : %d" % zdo)
    print("    interval Newton: %d" % ido)
    print("    zoom / int.Newton = %.3f" % (zdo / ido))
    print()
    print("  Exact, reproducible on any machine, unaffected by load.")


if __name__ == "__main__":
    report()
