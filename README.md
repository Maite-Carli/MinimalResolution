# MinimalResolution

This is a p=3 fork of [Guozhen Wang's code](https://github.com/pouiyter/MinimalResolution)
for computing the Adams-Novikov E<sub>2</sub> page for the sphere (originally
at p=2), using the algebraic Novikov spectral sequence. See
[this repository](https://github.com/ebelmont/ANSS_data) for sample data and an
explanation of how to interpret the data files output by the program.
Guozhen's original instructions are preserved at the bottom of this file.

This fork adds the ability to resolve **any** finitely generated
`BP_*BP`-comodule whose underlying BP_*-module is free over a quotient `BP_*/I_n`.

## Computing an E<sub>2</sub> page

### Build

Needs `g++` (C++11, OpenMP) and GMP (`libgmp-dev` — for `BPtab` only):

```sh
sh BPtable_compile     # -> BPtab         structure maps of the BP Hopf algebroid
sh BP_comod_compile    # -> mr_BP_comod   the resolver: any comodule, any height
```

That is the whole build. The other scripts are the original sphere-only path
(`st_compiling` → `mr_st`, `BP_compile` → `mr_BP`) and unrelated pipelines
(`kos_compile`, `ex_compile`, `e2p_compile`, `mot*_compile`, `tau*_compile`,
all hardcoded to p=2, plus `trunc_compile` for a toy Hopf algebroid). See
[`docs/BUILD_AND_RUN.md`](docs/BUILD_AND_RUN.md).

### The commands, end to end

Everything needed to produce both the tables and the charts for the sphere,
`S/α₁` and `S/p`:

```sh
# once, to build
sh BPtable_compile
sh BP_comod_compile

# once per degree range
./BPtab 30

# one run per comodule: computes the E2 page and writes the tables
./mr_BP_comod 30 8 sphere
./mr_BP_comod 30 8 alpha_1
./mr_BP_comod 30 8 mod_p

# one run per comodule: reads those tables and draws the chart
./anss_chart.py 30 -c sphere
./anss_chart.py 30 -c alpha_1
./anss_chart.py 30 -c mod_p
```

Run them in that order, in one directory. Nothing else is needed — in
particular **no `mr_st` run**; see "Do you need `mr_st` and `mr_BP`?" below.

### What each command produces

| command | produces |
|---|---|
| `sh BPtable_compile` | the executable `BPtab` |
| `sh BP_comod_compile` | the executable `mr_BP_comod` |
| `./BPtab 30` | `30_etaL`, `30_R2L`, `30_delta` — the structure maps `η_L`, `R2L`, `Δ` of the Hopf algebroid `(BP_*, BP_*BP)`. They describe the algebroid itself, not any comodule, so **one run serves every comodule and every height** |
| `./mr_BP_comod 30 8 sphere` | the tables `30_sphereBP…` (listed below), plus `30_sphereP…` intermediates |
| `./mr_BP_comod 30 8 alpha_1` | the same for `S/α₁`, prefixed `30_alpha_1BP…` |
| `./mr_BP_comod 30 8 mod_p` | the same for `S/p`, prefixed `30_mod_pBP…` |
| `./anss_chart.py 30 -c sphere` | the chart `30_sphereanss_E2.svg` |
| `./anss_chart.py 30 -c alpha_1` | `30_alpha_1anss_E2.svg` |
| `./anss_chart.py 30 -c mod_p` | `30_mod_panss_E2.svg` |

Each `mr_BP_comod` run writes its answer under the prefix
`<t>_<comodule>BP` and its internal phase-1 model under `<t>_<comodule>P`,
so runs for different comodules coexist in one directory without clobbering
each other. The readable files are:

| file | contents |
|---|---|
| `…AANSS_table.txt` | **the E<sub>2</sub> page** — one line per class, with the algebraic Novikov differentials |
| `…AANSS_h0.txt` | multiplication by `h₀ = t₁` (α₁) |
| `…AANSS_a0.txt` / `…AANSS_a<n>.txt` | multiplication by `p` at height 0; by `v_n` at height `n` |
| `…BocSS_table.txt` | the `v₀`-Bockstein at height 0; the `v_n`-Bockstein at height `n` |
| `…BocSS_h0.txt`, `…BocSS_a*.txt` | the same products on the Bockstein page |
| `…B2A_table.txt` | dictionary from Bockstein names to algebraic Novikov names |
| `…run_info.txt` | height, characteristic, `t`, length — read by `anss_chart.py` |

Everything else (`…res`, `…maps*`, `…back*`, `…_binary`, and the `P`-prefixed
files) is intermediate.

`anss_chart.py` recomputes nothing: it only reads the `.txt` tables from the
matching `mr_BP_comod` run. So you can redraw with different options at any
time without re-resolving.

### Choosing the numbers

`30` and `8` are **your choice**; nothing in the code prefers them. They are

- **`30` = the maximum internal degree `t`.** It must be the *same number* in
  `BPtab`, in every `mr_BP_comod` call, and in every `anss_chart.py` call,
  because it is used directly in the filenames each program looks for —
  `./BPtab 40` writes `40_etaL`, so `./mr_BP_comod 30 8 sphere` would not
  find it. It cuts the page off along the diagonal `t = stem + s`.
- **`8` = the resolution length**, i.e. how many steps of the minimal
  resolution to compute. It caps the *printed* second coordinate: classes
  appear while `s + i < 8`, where `i` is the algebraic Novikov filtration —
  so the homological degree `s` alone typically only reaches 4–5 at `8`.
  It appears only in the `mr_BP_comod` calls; the chart reads it back from
  the tables.

Both can be raised freely: `./BPtab 60` then `./mr_BP_comod 60 14 sphere`
then `./anss_chart.py 60 -c sphere` is a perfectly good larger run. The
length is the expensive parameter — each step is another resolution step over
`BP_*BP` — so raising `t` is usually cheaper per unit of chart. Note also
that only the `v_n`/`t_n` with `|v_n| = 2(3ⁿ−1) ≤ t` exist at all, so `v₂`
first appears at `t = 16` and `v₃` at `t = 52`.

### Comodules already available

| name | complex | `n` | rank | degrees |
|---|---|---|---|---|
| `sphere` *(default)* | `S` | 0 | 1 | `0` |
| `alpha_1` | `S/α₁ = cofib(S³ → S⁰)` | 0 | 2 | `0, 4` |
| `triv_01` | `S ∨ S¹` | 0 | 2 | `0, 1` |
| `mod_p` | `S/p` | 1 | 1 | `0` |
| `mod_p_v1` | `S/(p,v₁) = V(1)` | 2 | 1 | `0` |
| `alpha_1_mod_p` | `S/α₁ ∧ S/p` | 1 | 2 | `0, 4` |

To compute with any of them, put its name where `sphere` stands above — the
syntax is otherwise identical, and the same `./BPtab 30` still serves:

```sh
./mr_BP_comod 30 8 mod_p_v1        # V(1)
./anss_chart.py 30 -c mod_p_v1     # -> 30_mod_p_v1anss_E2.svg
```

`./mr_BP_comod --list` prints the available comodules with their heights,
and `--help` the usage.

The comodule name is optional for `mr_BP_comod` and defaults to `sphere`, so
`./mr_BP_comod 30 8` is the sphere. **`anss_chart.py` is not the same**: with
no `-c` it looks for `<t>_BP…`, the prefix a plain `mr_BP` run writes, so
charting an `mr_BP_comod` sphere run needs `-c sphere` explicitly —
`./anss_chart.py 30` on its own will report `30_BPAANSS_table.txt not
found`.

`n` is the **height**: the comodule's underlying module is free over
`BP_*/I_n`, where `I_n = (p, v₁, …, v_{n−1})`. `n = 0` means free over `BP_*`
itself. Since `I_n` is invariant, `(BP_*/I_n, BP_*BP/I_n)` is again a Hopf
algebroid and change of rings gives the same `Ext`, which is what lets
`BP_*/p` — free over `BP_*` for no `n` at all — go through the same
machinery. [`docs/GENERAL_COMODULES.md`](docs/GENERAL_COMODULES.md) has the
argument and what it costs. You do not have to supply `n`: it is recorded
with each comodule, and picked up automatically.

### Do you need `mr_st` and `mr_BP`?

Not for anything above — not even for the sphere. This trips people up
because Guozhen's original instructions (at the bottom of this file) require
`mr_st` before `mr_BP`, so it looks like a property of the sphere. It is not;
it is a property of *which driver you use*.

Both routes run the same two phases: resolve `M/I` over the field-based
`P = F₃[t₁, t₂, …]`, then lift that model to `BP_*BP`. They differ in where
phase 1 happens.

- **`mr_BP` is a two-program pipeline.** `mr_st` performs phase 1 and writes
  the model to disk as `<t>_gens_data` and `<t>_BPtables`; `mr_BP` reads
  exactly those back. Hence the ordering rule, and the caveat that `mr_st`'s
  second argument must be at least one larger than `mr_BP`'s.
- **`mr_BP_comod` does both phases in one process**, at `length + 1`
  internally. `BPtab` is its only prerequisite.

It cannot simply reuse an `mr_st` run, because phase 1 resolves `M/I`, which
is *comodule-specific*: `mr_st` only ever resolves the sphere's reduction,
its comodule being hardcoded exactly as `mr_BP`'s is. So for `S/α₁` or `S/p`
there is nothing in `<t>_gens_data` to share.

Two things still favour the old pair, if you are working on the sphere:

- `mr_st`'s model **is** reusable across `mr_BP` runs — resolve once at
  length 9, then lift at 8, 7, … without redoing phase 1, whereas
  `mr_BP_comod` redoes phase 1 every run.
- `mr_BP` also calls `mult_theta()`, the `θ`-multiplication table for the
  Moore spectrum, which `mr_BP_comod` deliberately skips: it carries a
  hardcoded table of `θ` degrees that is meaningless for an arbitrary
  comodule.

Otherwise `./mr_BP_comod <t> <length> sphere` reaches the same answer by a
shorter route — its tables are byte-identical to an `mr_st`/`mr_BP` run's.

### Adding your own comodule

Two edits, both in `comodules.cpp`; nothing else in the program changes.

**1. A builder.** `coaction_rows(i)` returns generator `i`'s coaction as
sparse `(j, c)` pairs — `j` another generator, `c` an element of `BP_*BP`.
Rows need not be sorted; `set_comodule` sorts and reduces them.

```cpp
static void build_myComplex(BP_Op &BP_oper, int &rank, std::vector<int> &degree,
                             std::function<vectors<matrix_index,BPBP>(int)> &coaction_rows){
    rank   = 2;              // generators over BP_*/I_n
    degree = {0, 4};         // FULL topological degrees: |v_1| = |t_1| = 4 at p=3

    BPBP t1 = BP_oper.h0();  // t_1, taken from the loaded structure tables

    coaction_rows = [&BP_oper, t1](int i) -> vectors<matrix_index,BPBP>{
        vectors<matrix_index,BPBP> row;
        BPBP one = BP_oper.BPBP_opers.unit(1);
        if(i == 0) row.push({(matrix_index)0, one});      // psi(x_0) = 1 (x) x_0
        else { row.push({(matrix_index)0, t1});           // psi(x_4) = t_1 (x) x_0
               row.push({(matrix_index)1, one}); }        //          + 1 (x) x_4
        return row;
    };
}
```

**2. One row in `comodule_table`**, whose last field is the height:

```cpp
{"myComplex", "one-line description, shown by --list", build_myComplex, 0},
```

Then `sh BP_comod_compile && ./mr_BP_comod 30 8 myComplex`.

#### Remarks

- **`BPBP`'s exponent slots are the opposite of what the type name
  suggests.** The *outer* exponent indexes the `v_i` via `η_R`; the *inner*
  one indexes the `t_i`. So `BPBP_opers.monomial(singleVar(1,1), unit(1))` is
  `η_R(v₁)`, **not** `t₁`. Use `BP_oper.h0()` for `t₁` and
  `BP_oper.thetas()` for `β₁` and friends. Getting this backwards fails
  silently — the warning at the top of `comodules.cpp` has the details.
- **Nothing checks the comodule axioms.** A coaction violating
  coassociativity or counitality will not crash; it will quietly produce a
  wrong `Ext`. Check by hand.
- **Nothing checks the height either, in either direction.** Too small and
  the machinery treats a module that is not free over `BP_*/I_n` as though it
  were. Too large and you get a perfectly self-consistent run of the wrong
  computation — declaring `sphere` at `n = 1` quietly computes `S/p`. State
  `n` from a proof that the module is free over `BP_*/I_n` and over nothing
  larger.

At height `n > 0` write the coaction in `BP_*BP/I_n` — in practice, exactly
what you would have written over `BP_*BP`, since `BP_oper` is already
reducing mod `I_n` by the time your builder runs.

**Find a prediction before trusting the output.** For `S/p`,
`Ext⁰ = F₃[v₁]`, so filtration 0 must be exactly one class at each of the
stems 0, 4, 8, …. For a two-cell complex, `ses_check.py` checks the whole
page against the one-cell run automatically (next section but one).

### What is out of reach

Only the invariant *prime* ideals `I_n`. `S/p^k` and `S/(p, v₁^k)` for
`k > 1` do live over invariant quotients — `(p^k)` and `(p, v₁^k)` are
invariant ideals — but those quotients are not polynomial rings, so they
would need a truncated multiplication rather than just a reduction, closer to
what `trunc_hopf.cpp` does.

## Charts

`anss_chart.py` draws an SVG chart of the Adams-Novikov E<sub>2</sub> page
from a finished `mr_BP` run (Python 3, no dependencies, nothing recomputed):

```sh
./mr_st 35 31 && ./BPtab 35 && ./mr_BP 35 30
./anss_chart.py 35            # -> 35_anss_E2.svg
```

Marks are plotted at `(t-s, s)` — note this is *not* the pair `mr_BP` prints,
which is `(t-s, s+i)` for `i` the algebraic Novikov filtration. There is one
mark per **cyclic summand**, not per class, following the convention of the
published p=3 charts: a filled square marked ∞ for `Z_(3)`, a filled dot for
`Z/3`, a dot inside `n-1` rings for `Z/3^n` (the 3-multiplications come from
the `a0` table), and a dashed outer ring where the range cuts a tower off.
Each mark is labelled with its generator's name from the table
(`v₀v₁³[1-1]`), and tan slope-1/3 lines are multiplication by α₁. Works for
any comodule via `-c`. See [`docs/CHARTS.md`](docs/CHARTS.md) for the full
convention and the limitations.

## Checking a two-cell complex against the sphere

There is no published E<sub>2</sub> chart for `S/α₁` to compare
`mr_BP_comod`'s output with, but the cofiber sequence determines it from the
sphere's page: the connecting map of the long exact sequence is
multiplication by α₁, so the answer is `coker(α₁)` on the bottom cell plus
`ker(α₁)` shifted four stems for the top cell. `ses_check.py` carries that
comparison out from finished runs:

```sh
./BPtab 35 && ./mr_BP_comod 35 12 sphere && ./mr_BP_comod 35 12 alpha_1
./ses_check.py 35 -c alpha_1     # 21 bidegrees, all sharp, 0 mismatches
```

[`docs/SES_CHECK.md`](docs/SES_CHECK.md) derives this, and in particular
explains where the *non-split* comodule structure of `BP_*(S/α₁)` enters — it
is carried entirely by that connecting map.

## Documentation

This codebase had no architecture documentation beyond this README and
inline comments. [`docs/index.html`](docs/index.html) is a generated
documentation dashboard covering the class structure, module relationships,
and full build/run pipeline — start there, or jump straight to
[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) for the big picture. It's now
grounded directly in `MinimalResolution.pdf` (the algorithm writeup
referenced below), with citations to its definitions and propositions
throughout rather than guesswork. It also flags a few things worth knowing
before you dig in: several sub-pipelines (`kos`, `mr_ex`, and the whole
motivic pipeline) are hardcoded to p=2 — the docs explain why this is an
independent cross-check rather than an unfinished port — and there's some
dead/duplicate code left over from earlier refactors — see
`docs/ARCHITECTURE.md` §6–7 for specifics.

******************************************************************************************************

The algorithm is explained in the pdf file MinimalResolution.pdf

The codes can be compiled with GCC. The GNU Multiple Precision Arithmetic Library should be installed.

******************************************************************************************************

To compile, run the following batch files:

sh st_compiling

sh BPtable_complile

sh BP_compile

*******************************************************************************************************

To get the minimal resolution for BP/I, for t<=25, s<=21 (say), run

./mr_st 25 21

To get the structure maps of the BP Hopf algebroid for t<=25, run

./BPtab 25

To get the minimal resolution for BP, for t<=25, s<=20, run

./mr_BP 25 20

*******************************************************************************************************

Warning:

The first input parameter is the maximal internal degree t.

The three executalbes are dependent, and should be run in the above order. 

The s for the minimal resolution for BP/I should be at least one larger than the s for that of BP.

Any mistake of the input could result in unpredictible behaviour, usually a break-down of the program such as a segmentation error.

[EB asked here whether the p=3 version has different restrictions on the
degrees. It does: this fork's first parameter is the internal degree t
itself, not half of it, so the three examples above are truncated at t<=25
rather than t<=50. Guozhen's original text said "half of t"; that is a p=2
carry-over (the p=2 code is not in this repo to check against, but its
degrees 2(2^n - 1) are all even, so a halved convention there is plausible).
Here exponents.cpp:31 holds the full topological degrees
|v_n| = |t_n| = 2(3^n - 1), so |v_1| = 4, and monomial_index keeps every
monomial of degree d <= argv[1] (mon_index.cpp:44) on both the BP side
(BP.cpp:8) and the Steenrod side (steenrod.cpp:5).

Verified: ./mr_st 24 6 tops out at generator degree 24, and ./mr_BP_comod
24 8 produces classes of internal degree t = 24 (beta_1^2 = [4-0], printed
|deg=(20,4)), while the same run at 23 stops at t = 20. A degree-t class
needs argv[1] >= t, and note that only the v_n/t_n with |v_n| <= argv[1] are
built at all (mon_index.cpp:12-15, used at BPtable.cpp:12), so v_2/t_2 first
appear at 16 and v_3/t_3 at 52.

Classes are printed as |deg=(t-s, s+i): the stem t-s, then the homological
degree s plus the algebraic Novikov filtration i. See docs/CHARTS.md.]
