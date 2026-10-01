# Where this work stands

A handoff note covering the work done on `S/α₁`, the charts, and the
degree-convention and `mod I` fixes. Written so that a fresh session (or a
reader six months from now) can pick up without re-deriving what cost time to
establish. Everything described here is on `prime3`, commits `74bbdf1` ..
`4ca4485`.

## 1. What was added

| file | what it is |
|---|---|
| `ses_check.py` | checks a two-cell complex's E2 page against the sphere's via the cofiber long exact sequence; exits 1 on a mismatch |
| [`SES_CHECK.md`](SES_CHECK.md) | the derivation behind it, the controls, and the limitations |
| `anss_chart.py` | rewritten: one mark per **cyclic summand** (box / dot / rings), labelled, in the convention of Belmont's published p=3 chart |
| [`CHARTS.md`](CHARTS.md) | rewritten around that, plus what the glyphs prove and what they assume |
| `BP_mod_I.cpp` | **bug fix**: the two exponent slots were transposed (§5) |
| `README.md`, [`BUILD_AND_RUN.md`](BUILD_AND_RUN.md), [`pipelines/BP.md`](pipelines/BP.md), [`pipelines/STEENROD.md`](pipelines/STEENROD.md) | the degree-convention correction (§3) |

Quick start after a pull (needs `g++` with OpenMP, `libgmp-dev`, `python3`):

```sh
sh st_compiling && sh BPtable_compile && sh BP_compile && sh BP_comod_compile
./mr_st 35 13 && ./BPtab 35 && ./mr_BP 35 12      # classic pipeline
./mr_BP_comod 35 12 alpha_1                       # S/α₁
./anss_chart.py 35 -c alpha_1                     # chart
./ses_check.py 35 -c alpha_1                      # 21/21 bidegrees agree
```

## 2. Gradings — the thing to get right first

Following Rognes's Definition 4.5 (`Adams_sseq.pdf`, on this branch): an
element of `E_r^{s,t}` has **filtration `s`**, **total degree `t-s`** (the
stem), **internal degree `t`**. Charts use `(t-s, s)` coordinates.

The program prints a *third* grading alongside these, and that is the usual
source of confusion:

```
v0^1v1^3[1-1]    |deg=(a, b)
                   a = t - s          the stem
                   b = s + i          NOT the filtration
```

- `i` = **algebraic Novikov filtration** = the total number of `v`'s in the
  name, `v0 = p` included (`algNov.cpp:10-16`; `MinimalResolution.pdf` p. 9).
- `[s-j]` = the `j`-th generator of the `s`-th term of the minimal
  resolution — exactly Rognes's `γ_{s,j}`. **`s` is read off this bracket**,
  never off `b`. `anss_chart.py:81-83` does precisely that.
- `t = a + s`, recoverable, never printed.

The classic trap: `v1^1[1-0]` is α₂, filtration 1, printed at height 2.

## 3. The command-line arguments

**`argv[1]` is the maximal internal degree `t`, not half of it.** The README
said "half of t"; that was inherited from Guozhen Wang's p=2 code and is
wrong here. `exponents.cpp:31` stores full topological degrees
(`|v_n| = |t_n| = 2(3^n-1)`, so `|v_1| = 4`) and `monomial_index` keeps every
monomial of degree `d <= argv[1]` (`mon_index.cpp:44`), on both the BP side
(`BP.cpp:8`) and the Steenrod side (`steenrod.cpp:5`). Measured, not inferred:
`./mr_st 24 6` stops at generator degree 24; `./mr_BP_comod 24 8` has β₁² at
`t = 24` while the same run at 23 stops at `t = 20`.

`argv[2]` is the resolution length `L`, and it caps the *other* direction:
only classes with `s + i < L` are printed (`algNov.cpp:201-205`,
`BP_init.cpp:34`). That is why the `v0`-towers stop where they do — 12 steps
in stem 0 at `L = 12`.

Consequence worth remembering: since `t = stem + s`, **the computed region is
cut off along a diagonal**, not at the top of either axis. At `halfT=35` the
sphere has a class in stem 31 (`t = 32`) but not in stem 30, because β₁³ sits
at `s = 6`, i.e. `t = 36`. The chart draws that diagonal now.

## 4. `BPBP`'s exponent slots

`BPBP = polynomial<BP>` is **not** "polynomial in the `t_i` with `BP_*`
coefficients". The **outer** exponent indexes the `v_i` and the **inner**
(coefficient) exponent indexes the `t_i` (`BP.cpp:69`, and the warning atop
`comodules.cpp`). So:

- `BP_oper.h0()` is `t_1` — outer 0, inner exponent 1. Use it.
- `BPBP_opers.monomial(singleVar(1,1), unit(1))` is `η_R(v_1)`, **not** `t_1`.

## 5. The bug that was found and fixed (`35c6370`)

`reduce_BPBP_mod_I` was written for the layout the type name suggests, so it
read both slots backwards. Probed against the real types:

```
before:  t_1 -> 0        η_R(v_1) -> t_1
after:   t_1 -> t_1      η_R(v_1) -> 0
```

Every off-diagonal coaction entry is a polynomial in the `t_i`, so the old
code reduced all of them to zero and handed phase 2 a **split** comodule —
the mod-`I` model of `S⁰ ∨ S⁴` in place of `S/α₁`.

Effect on results, measured: `sphere` and `triv_01` are byte-identical before
and after (their coaction is `1`, which reduces correctly either way, so the
documented `mr_BP`-vs-`mr_BP_comod` regression could never catch this).
`alpha_1`'s tables changed **only in generator indices** inside class names
(`[1-1]→[1-0]`, …); the degree profile, the summand structure and the
α₁-multiplication edges are identical, and `ses_check.py` still reports 21/21.
The BP side always used the true coaction and the model only guided the
generator search — so the damage was indexing, not answers. That is luck, not
a guarantee.

**Two values pin the convention down; re-check them after touching that
file**: `h0()` must reduce to `t_1`, and `η_R(v_1)` must reduce to `0`.

## 6. The mathematics established

### The cofiber sequence

`BP_*` is even, so `α₁` acts as 0 on BP-homology and

    0 --> BP_* --> BP_*(S/α₁) --> Σ⁴BP_* --> 0

(shift 4, not 3 — the top *cell* is `e⁴`). Split over `BP_*`, **not** over
`BP_*BP`. The `BP_*`-splitting is what keeps the cobar complexes exact, hence
gives a long exact sequence; the failure of a comodule splitting is the
connecting map, which is the Yoneda product with the extension class — the
off-diagonal `t_1`, i.e. `α₁`. In chart coordinates:

    0 -> coker(α₁ into (n,s)) -> E2(S/α₁)(n,s) -> ker(α₁ out of (n-4,s)) -> 0

α₁-multiplication is `(stem +3, s +1)`; the top cell shifts `(stem +4, s +0)`.

### The (4,0) class — the one that looks wrong and isn't

`ker(α₁ : Ext^{0,0} → Ext^{1,4}) = 3Z_(3) ≠ 0`, because `α₁` has order 3. So
there **is** a class in stem 4, filtration 0. Verified by hand, independently
of the LES and of the code: with `η_R(v_1) = v_1 + 3t_1`,

    Ext^{0,4}(BP_*(S/α₁)) = Prim(M)_4 = Z_(3){ v_1x_0 + 3x_4 }

and the same computation gives `Prim(M)_8 = 0`, matching the chart. Stem 0 is
the only place this can happen — it is the only bidegree where the sphere's
Ext is not 3-torsion.

### Group structures

A non-zero `a0` entry is sound (`3x ≠ 0`); a `-> o` entry is not (`3x` may
have jumped filtration); the class count is exact. So **one `a0` chain
exhausting a bidegree pins the group down** (length `n` plus an element of
order `3^n` forces cyclic). That is how `S/α₁` gets `Z/9` in (7,1) and `Z/27`
in (11,1). The Bockstein tables confirm both independently: `v0^2[1-0]` dies
by a `d2`, `v0^3[1-1]` by a `d3`.

Where this is *not* settled: bidegrees holding several chains. At `halfT=35`
there are none (all 15 sphere and 14 `S/α₁` bidegrees are exact); at
`halfT=185, L=14` the sphere has 168 exact, 37 lower bounds (truncated) and
64 open. `anss_chart.py` now reports that split after every run.

## 7. What is verified, and how to re-verify

| check | command | result |
|---|---|---|
| generic path == shipped `mr_BP` | `diff 35_BPAANSS_table.txt 35_sphereBPAANSS_table.txt` | byte-identical (algNov and Bockstein) |
| `S/α₁` against the LES | `./ses_check.py 35 -c alpha_1` | 21/21 bidegrees, all sharp |
| same at other ranges | `20 / 6`, `28 / 9` | 11/11, 17/17 |
| split control | `./ses_check.py 35 -c triv_01 --attaching none` | 28/28 |
| negative control | `./ses_check.py 35 -c alpha_1 --attaching none --top-cell 4` | 8 mismatches, as it should |
| sphere chart vs published charts | `./anss_chart.py 35` | matches, including `Z/9` at stems 11 and 23 |

A second negative control worth knowing about: adding a scratch comodule with
`rank = 2`, `degree = {0,4}` and the **identity** coaction (`S⁰ ∨ S⁴`) and
checking it against the `S/α₁` prediction fails in the same 8 bidegrees. Its
chart shows `[0-1]` — the top cell's own generator — where `S/α₁` shows
`v1^1[0-0]`.

## 8. Open threads

- **β₁-divisibility** (Belmont's chart colours those classes blue) and **α₂
  lines** are not available: `BPInit::mult_table` only builds the `h0` table
  on the algebraic Novikov side, and `mult_theta`'s `theta_i` tables are
  Bockstein-side. Both would be a `mult_table(<class>, <degree>, "<name>.txt")`
  call in `mr_BP_comod` plus one branch in the chart script, which already
  parses any multiplication table. `BP_oper.thetas()` has β₁.
- **The 64 ambiguous bidegrees at 185.** The Bockstein tables
  (`BocSS_table.txt`, with `B2A_table.txt` translating names) are the
  independent record of `p`-divisibility. The chart does not read them because
  `Boc.cpp:37-38` prints `2t - s` and drops the second coordinate — fixable
  there, or reconstructible from names.
- **A startup self-check** asserting `h0() -> t_1` and `η_R(v_1) -> 0` would
  turn a repeat of §5 into a loud crash instead of a silent wrong model. Not
  added; five lines in `mr_BP_comod.cpp` if wanted.
- **Massey products** (the grey dashed lines in Belmont's chart) are out of
  reach from these tables.
- **`claude/code-output-explanation-vy84sm`** is an unmerged branch from a
  different session, "Document how to build, run and chart from scratch". It
  likely overlaps the README and `CHARTS.md` edits merged here, so expect
  conflicts and review it before merging.
- Comodules are registered in `comodules.cpp`; `mr_BP_comod --list` shows
  them. Adding one is a builder function plus a table row. Nothing anywhere
  checks the comodule axioms.

## 9. Practical notes

- `libgmp-dev` is needed for `BPtab` only; the two Python scripts need no
  packages.
- Runtimes on a cloud container: `35 / 12` is seconds; `BPtab 185` ≈ 8 s and
  `mr_BP_comod 185 14 sphere` ≈ 209 s, giving 620 classes out to stem 182.
- The `<halfT>` placeholder in filenames and CLI flags is a misnomer kept for
  compatibility; its value is `t`.
