# The S/α₁ work: what was established

Companion to [`../CLAUDE.md`](../CLAUDE.md), which is the orientation note for
the repository as a whole — read that first for build/run, the two-phase
algorithm, the heights mechanism, the traps, and the checking recipe. This
page covers what a separate strand of work established about **`S/α₁`, the
gradings, and the group structures**, with the arguments rather than just the
conclusions, so they do not have to be re-derived. Commits `74bbdf1` ..
`4ca4485` on `prime3`.

## 0. What happened, in brief

The question that started it: **how do you check the code's `S/α₁` output,
when no E₂ page for it is published?** The answer is that one does not need a
published chart — the cofiber sequence determines the page from the sphere's,
which *is* known and which this repo computes. That became `ses_check.py` and
[`SES_CHECK.md`](SES_CHECK.md), and the run agrees in every bidegree in range
(§4).

Four corrections came out of the same work:

1. **The degree convention.** `argv[1]` is the maximal internal degree `t`,
   not half of it; the README's claim was inherited from the p=2 original.
   Measured against runs at 12, 23, 24 and 40, then fixed in the README,
   `BUILD_AND_RUN.md` and the two pipeline docs that carried it as an open
   `TODO(math)`.
2. **The charts.** `anss_chart.py` drew one dot per class, which
   misrepresents a `Z/9` as two dots and hides a `Z_(3)` entirely. It now
   draws one mark per cyclic summand — box, dot, concentric rings — in the
   convention of Belmont's published p=3 chart, labelled by generator name.
3. **A real bug**, raised by another agent and confirmed by a probe against
   the live types: `reduce_BPBP_mod_I` had `BPBP`'s exponent slots
   transposed, so every off-diagonal coaction entry reduced to zero and phase
   1 silently received a *split* comodule. Fixed in `35c6370`; see
   `../CLAUDE.md` trap 2. The sphere regression structurally could not catch
   it — the sphere's coaction is `1`, which reduces correctly either way.
4. **The stem-4 class**, which was challenged as obviously wrong and turned
   out to be right: `ker(α₁ : Ext^{0,0} → Ext^{1,4}) = 3Z_(3) ≠ 0` because α₁
   has order 3. Settled by computing the primitive by hand, independently of
   both the long exact sequence and the code (§2).

Two conclusions worth carrying forward. First, the checks that catch things
here are the ones with an *independently known* answer — the sphere
regression, the long exact sequence, a hand-computed `Ext^0`, the Bockstein
tables — because nothing in this pipeline validates its own input and almost
every failure mode is silent. Second, three of the four corrections above
were found by running something and comparing, not by reading the code.

## 1. The three gradings

Following Rognes's Definition 4.5 (`Adams_sseq.pdf`, in the repo root): an
element of `E_r^{s,t}` has **filtration `s`**, **total degree `t-s`** (the
stem), **internal degree `t`**. Charts use `(t-s, s)`.

The program prints a *third* grading beside those, which is the usual source
of confusion:

```
v0^1v1^3[1-1]    |deg=(a, b)        a = t - s      the stem
                                    b = s + i      NOT the filtration
```

| piece | meaning |
|---|---|
| `i` | algebraic Novikov filtration = total number of `v`'s in the name, `v0 = p` included (`algNov.cpp:10-16`, `MinimalResolution.pdf` p. 9) |
| `[s-j]` | the `j`-th generator of the `s`-th term of the minimal resolution — Rognes's `γ_{s,j}`. **`s` is read off this bracket**, never off `b` (`anss_chart.py:81-83`) |
| `t` | `= a + s`; recoverable, never printed |

The standing trap: `v1^1[1-0]` is α₂, filtration 1, printed at height 2.

Because `t = stem + s`, the computed region is **cut off along a diagonal**,
not at the top of either axis — at `t = 35` the sphere has a class in stem 31
(`t = 32`) but not in stem 30, where β₁³ sits at `s = 6`, i.e. `t = 36`. The
chart draws that diagonal. A class missing high up usually means raise the
*first* argument, not the second.

## 2. The cofiber sequence, and the class that looks wrong

`BP_*` is even, so α₁ acts as 0 on BP-homology and

    0 --> BP_* --> BP_*(S/α₁) --> Σ⁴BP_* --> 0

— shift 4, not 3: the top *cell* is `e⁴`. Split over `BP_*`, **not** over
`BP_*BP`. The `BP_*`-splitting keeps the cobar complexes exact, hence gives a
long exact sequence; the failure of a comodule splitting *is* the connecting
map, which is the Yoneda product with the extension class — the off-diagonal
`t_1`, i.e. α₁. See [`SES_CHECK.md`](SES_CHECK.md) for the full derivation and
for `ses_check.py`, which checks a run against it bidegree by bidegree.

The one result worth recording here, because it is the thing that looks like a
bug and is not: **there is a class in stem 4, filtration 0**, because
`ker(α₁ : Ext^{0,0} → Ext^{1,4}) = 3Z_(3) ≠ 0` — α₁ has order 3, so three
times the top cell survives even though the top cell itself does not. Checked
by hand, independently of the long exact sequence and of the code, using
`η_R(v_1) = v_1 + 3t_1` and `D(am) = η_L(a)D(m) + (η_L - η_R)(a) ⊗ m`:

    Ext^{0,4}(BP_*(S/α₁)) = Prim(M)_4 = Z_(3){ v_1x_0 + 3x_4 }
    Ext^{0,8}(BP_*(S/α₁)) = 0                   (same computation, degree 8)

both matching the run. The chart names that generator `v1^1[0-0]`, after its
leading term `v_1x_0` — the names are coordinates in the *resolution*, where
no `x_4` exists, not in `M`. For `S⁰ ∨ S⁴` the same bidegree instead shows
`[0-1]`, the top cell's own generator: that difference is exactly the `+3x_4`.

Stem 0 is the only bidegree where this can happen, since it is the only place
the sphere's Ext is not 3-torsion.

## 3. Which group structures are proved

| evidence | strength |
|---|---|
| a non-zero `a0` entry `x -> y` | **sound**: the leading term of `3x` is non-zero, so `3x ≠ 0` |
| an `a0` entry `x -> o` | **not sound**: `3x` may have jumped filtration — a hidden extension |
| the class count in a bidegree | **exact**: the classes are an F₃-basis of `gr Ext` |

So a chain gives a lower bound on the order and the count gives the exact
length, and **when one chain exhausts a bidegree the group is pinned down** —
length `n` plus an element of order `3^n` forces cyclic. That is why `S/α₁`
has `Z/9` in (7,1) and `Z/27` in (11,1), not `Z/3 ⊕ Z/3` and `Z/9 ⊕ Z/3`. The
Bockstein tables confirm both independently: `v0^2[1-0]` dies by a `d2`,
`v0^3[1-1]` by a `d3`.

Bidegrees holding several chains are *not* settled — a hidden extension could
merge them. `anss_chart.py` reports the census after every run: at `t=35` all
15 sphere and 14 `S/α₁` bidegrees are exact; at `t=185, L=14` the sphere has
168 exact, 37 lower bounds (truncated towers), 64 open. The Bockstein tables
(`BocSS_table.txt`, with `B2A_table.txt` translating names) are the way into
the open ones; the chart does not read them because `Boc.cpp:37-38` prints
`2t - s` and drops the second coordinate.

All of this is the **characteristic-0** story. Over `BP_*/I_n` every summand
is `Z/p`; the chart reads which case it is from `run_info.txt`.

## 4. Results on record

| check | expected |
|---|---|
| `./ses_check.py 35 -c alpha_1` | 21 bidegrees agree, all sharp |
| the same at `t=20` / `t=28` | 11/11, 17/17 |
| `./ses_check.py 35 -c triv_01 --attaching none` | 28/28 (split control) |
| `./ses_check.py 35 -c alpha_1 --attaching none --top-cell 4` | 8 mismatches — the negative control |
| sphere chart at `t=35` | matches the published p=3 charts, including `Z/9` at stems 11 and 23 |

A second negative control worth knowing: a scratch comodule with `rank = 2`,
`degree = {0,4}` and the **identity** coaction (`S⁰ ∨ S⁴`), checked against
the `S/α₁` prediction, fails in the same 8 bidegrees.

## 5. Still open from this strand

- **β₁-divisibility** (the blue classes in Belmont's published chart) and
  **α₂ lines**: both need a `mult_table(<class>, <degree>, "<name>.txt")` call
  in `mr_BP_comod` plus one branch in `anss_chart.py`, which already parses
  any multiplication table. `BP_oper.thetas()` has β₁.
- **Massey products** (her grey dashed lines) are out of reach from these
  tables.
- **A startup self-check** asserting `h0() -> t_1` and `η_R(v_1) -> 0` would
  turn a repeat of the `reduce_BPBP_mod_I` bug into a loud crash instead of a
  silent split comodule. Five lines in `mr_BP_comod.cpp`; not added.
