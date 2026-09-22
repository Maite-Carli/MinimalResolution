# Charting the E2 page

`anss_chart.py` turns a finished `mr_BP` run into an SVG chart of the
Adams–Novikov E2 page. It reads only the text tables `mr_BP` already writes,
so it never recomputes anything and needs no rebuild — you can point it at
old output, including output from other machines.

```sh
./mr_st 35 31 && ./BPtab 35 && ./mr_BP 35 30    # the usual pipeline
./anss_chart.py 35                              # -> 35_anss_E2.svg
```

Requires Python 3 and nothing else. The argument is the same `<halfT>` you
gave the pipeline; that's how it finds the files. Run it in the directory
holding the run, or pass `--dir`.

### Charting a comodule other than the sphere

`mr_BP_comod` (see [`GENERAL_COMODULES.md`](GENERAL_COMODULES.md)) writes its
tables with a `<halfT>_<comodule>BP...` prefix rather than `<halfT>_BP...`.
Pass `-c/--comodule` to chart one:

```sh
./BPtab 20 && ./mr_BP_comod 20 4 alpha_1   # resolve S/alpha_1
./anss_chart.py 20 -c alpha_1              # -> 20_alpha_1anss_E2.svg
```

The comodule name lands in the output filename and the chart title, so
charts for different complexes don't overwrite each other. Nothing in the
script is specific to the sphere: the glyphs, groups and structure lines
below are all derived from whichever run's tables it is pointed at.

One caveat when reading such a chart: `mr_BP_comod` runs `mult_table()` (so
the `α₁` structure lines are drawn) but deliberately **not** `mult_theta()`,
which is specific to the Moore spectrum and carries a hardcoded table of
theta degrees. The `theta_i` tables are Bockstein-side and aren't drawn on an
ANSS chart anyway (see §4), so this costs nothing here.

## 1. Grading convention

The default (`--grading anss`) plots the Adams–Novikov bidegree:

- **x** = stem = `t - s`
- **y** = `s`, the homological degree

**This is not the pair `mr_BP` prints.** `algNov.cpp:118` emits
`|deg=(t-s, s+i)`, where `i` is the algebraic Novikov filtration, so
`v1^1[1-0]` (which is α₂, in stem 7 of Ext¹) is printed at height **2** and
belongs at height **1**. The chart takes the stem from the printed pair but
reads `s` off the `[s-n]` bracket in the class name instead. The internal
degree is recoverable as `t = x + y`.

`--grading algnov` plots `mr_BP`'s printed pair as-is, one dot per class,
keeps the 3-multiples as vertical towers, and draws the algebraic Novikov
differentials. That view is useful for checking a run against the raw
tables; it is not an ANSS chart.

## 2. One mark per cyclic summand

The classes in the table are an F₃-basis of `gr Ext`, so a `Z/9` appears as
*two* of them. Drawing one dot each would misrepresent the group, so the
ANSS view groups them: multiplication by 3 is exactly what
`<prefix>AANSS_a0.txt` records, and chaining those edges recovers each cyclic
summand. A class that is not 3 times any *surviving* class generates a
summand; the chain of its 3-multiples gives the order.

The glyphs follow Belmont's published 3-primary ANSS chart:

| glyph | meaning |
|---|---|
| filled square, marked ∞ | `Z_(3)` |
| filled dot | `Z/3` |
| dot inside `n-1` concentric rings | `Z/3^n` |
| dashed outermost ring | the tower runs off the end of the computed range, so the order drawn is a **lower bound** |

Several summands in one bidegree are drawn side by side. Each is labelled
with the name of its generator, typeset the way the published charts do —
`v₀v₁³[1-1]` for the table's `v0^1v1^3[1-1]`. `--raw-labels` keeps the
table's exact text, `--no-labels` drops labels entirely, and `--no-legend`
drops the key under the axis.

Two rules decide the glyph, and they come from different places:

- **Filtration 0 is always `Z_(3)`**, by mathematics rather than by the
  tables: `Ext^0 = Prim(M)` is a submodule of a free `BP_*`-module, hence
  torsion-free. (This is why `S/α₁` has a box in stem 4: `Ext^{0,4}` there is
  `Z_(3){v_1x_0 + 3x_4}` — see [`SES_CHECK.md`](SES_CHECK.md) §4.)
- **Above filtration 0** the order comes from the `a0` chain. Where that
  chain hits the truncation the mark gets the dashed ring: the printed
  second coordinate `s + i` is cut off at the resolution length (§6), and the
  `a0` table is pruned one step earlier still (`multiplication.cpp:169`), so
  near the top of the filtration range the top of a tower is simply not
  visible. At `halfT=185, L=14`, 60 of 344 summands are in that situation.

A class killed by an algebraic Novikov differential is dropped, along with
the tag that killed it: lines containing `<-` record a differential and both
ends die. (The tag never gets a line of its own — `SS_table::output`,
`SS.h:139-146`, only emits untagged entries — so dropping the `<-` lines
removes both.)

Note that 3-multiples must be found from the `a0` table, **not** by looking
for `v0` in the name. `v0^1[1-1]` carries a `v0` but is not 3 times anything
that survives — its predecessor `[1-1]` supports a d₂ — and it is the
generator of the `Z/9` in stem 11. A syntactic filter would hide it.

If the `a0` file is missing the script says so and draws every class as a
plain dot, since without it no isomorphism type is knowable.

## 3. Structure lines

Solid tan lines are multiplication by `h0 = (η_R(v1) - η_L(v1))/p` = α₁,
read from `<prefix>AANSS_h0.txt`, so they run `(+3, +1)` — slope 1/3, as in
the published charts. Products between classes inside the same summand are
suppressed; each line joins two marks.

Multiplication by α₂ (`(+7, +1)` lines) and β₁-divisibility (which Belmont's
chart indicates by colouring the classes) are not currently available:
`BPInit::mult_table` (`BP_init.cpp:180-187`) only computes the `h0` table for
the algebraic Novikov side, and the `theta_i` tables `mult_theta` writes are
Bockstein-side. Adding either means calling `mult_table(<class>, <degree>,
"<name>.txt")` with the appropriate `BPBP` element (`BP_oper.thetas()` has
β₁) and teaching this script the extra file — the parsing side already
handles any multiplication table. Massey products, which that chart also
draws, are out of reach here entirely.

## 4. Bockstein tables are not supported

`Boc.cpp:37-38` computes `(degree*2 - k, k + filtration)` but prints only the
first coordinate, and that coordinate is `2t - s`, not the stem — for α₁ it
prints `7` where the algebraic Novikov table prints stem `3`. So the
Bockstein `.txt` files are missing the second coordinate entirely and use a
different first one. Charting them needs either a fix at that line or a
reconstruction from the class names (`Boc_table::v_valuation`, `Boc.cpp:4-6`,
counts only powers of `v0`, so `y = s + v0-exponent`).

## 5. Verification

At `halfT=35, s=30` the chart's marks are

```
(0,0) box        (3,1) (7,1) (10,2) (11,1)°° (13,3) (15,1) (19,1) (20,4)
(23,1)°° (23,5) (26,2) (27,1) (29,3) (31,1)
```

i.e. 1, α₁, α₂, β₁, α_{3/2}, α₁β₁, α₄, α₅, β₁², α_{6/2}, α₁β₁², β₂, α₇,
α₁β₂, α₈ — agreeing with the published p=3 ANSS charts on every stem in
range. The two double-ringed marks (°°) are stems 11 and 23, which are
exactly the two α-family classes of order 9 in that range: `α_{i/j}` has
order `3^j`, and `j = 2` precisely when `3 | i`.

For `S/α₁` the same chart has boxes in stems 0 **and 4** and nothing in
stems 3, 13, 23, 29 — see [`SES_CHECK.md`](SES_CHECK.md), which checks the
whole page against the sphere's automatically.

## 6. Truncation

Classes near the top stems are missing because of the degree bound, not
because they are absent — `mr_BP` prints things like `out of range for beta1`
when a needed class falls outside the budget. The chart cannot tell the
difference, so treat the right-hand edge as unreliable and crop it with
`--max-stem` when showing the chart to anyone. The same applies upwards in
filtration: `s + i` is capped by `mr_BP`'s second argument, which is what
the dashed rings of §2 are about.

## 7. Options

| Flag | Effect |
|---|---|
| `-d, --dir` | directory holding the tables (default `.`) |
| `-c, --comodule` | chart an `mr_BP_comod` run for this comodule (e.g. `alpha_1`); default is a plain `mr_BP` run |
| `-o, --output` | output file (default `<halfT>_<comodule>anss_E2.svg`) |
| `--grading anss\|algnov` | see §1 |
| `--no-labels` | drop the class names |
| `--raw-labels` | label with the table text verbatim (`v1^3[1-0]`) |
| `--label-size` | font size for the names (default 5.5) |
| `--no-legend` | drop the glyph key under the axis |
| `--max-rings` | most rings drawn before the order is written beside the mark (default 4) |
| `--ring-gap` | px between rings (default 2.0) |
| `--omit-stem0` | drop `Ext^0` in stem 0, as the published charts do |
| `--max-stem`, `--max-filt` | crop |
| `--unit` | px per lattice step (default 31.2, matching published charts) |
| `--dot-radius`, `--tick`, `--title`, `--no-title` | cosmetics |

Every mark carries an SVG `<title>`, so hovering it in a browser shows the
class name, the group, and both gradings.
