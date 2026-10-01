# Working on this repository

Context carried forward from the session that added the height mechanism
(September 2026). Written for whoever — human or Claude — picks this up next.

## What this is

A p=3 fork of Guozhen Wang's code for the Adams–Novikov E<sub>2</sub> page,
i.e. `Ext_{BP_*BP}(BP_*, M)`, computed through the algebraic Novikov spectral
sequence. The algorithm is `MinimalResolution.pdf`; `docs/` documents the
code. Originally sphere-only; it now resolves any finitely generated
`BP_*BP`-comodule whose underlying module is free over some `BP_*/I_n`.

Everything lives on branch **`prime3`**. There is no CI and no test suite —
verification is the recipe in "Checking a change" below, and it matters,
because almost every failure mode here is silent.

## Build and run

```sh
sh BPtable_compile     # -> BPtab        needs GMP (libgmp-dev); BPtab only
sh BP_comod_compile    # -> mr_BP_comod  the resolver
./BPtab 30             # structure maps; one run serves every comodule/height
./mr_BP_comod 30 8 mod_p
./anss_chart.py 30 -c mod_p
```

`30` = max internal degree `t`, and must be identical in every command (it is
part of the filenames). `8` = resolution length, which caps the **printed**
second coordinate `s + i < 8` (`i` = algebraic Novikov filtration), not `s`.
See the README for the full operational account.

`mr_st` and `mr_BP` are the original sphere-only pair and are **not** needed
by `mr_BP_comod`, not even for the sphere: `mr_BP` reads a phase-1 model that
`mr_st` wrote to disk (`<t>_gens_data`, `<t>_BPtables`), whereas
`mr_BP_comod` runs phase 1 itself at `length + 1`. It cannot reuse `mr_st`'s
output anyway, because phase 1 resolves `M/I`, which is comodule-specific.
`kos`, `mr_ex` and the motivic pipelines are hardcoded to p=2 and unrelated.

## The algorithm, compressed

A pivot search over `BP_* = Z_(3)[v_1, v_2, …]` is impossible:
`PolyOp::invertible` (`polynomial/6.h:27`) only recognises invertible bare
constants, so almost nothing a real coaction produces is ever a pivot. Hence
two phases, which is Wang's central idea (minimality is a mod-`I` condition):

1. resolve `M/I` over `P = BP_*BP/I = F_3[t_1, t_2, …]` from scratch —
   `pre_resolution_tab`, correct because `F_3` is a field;
2. lift that model to `BP_*BP` — `pre_resolution_modeled`, which replays the
   model's cogenerator choices and pivots and only lifts values.

Then `Hom(BP_*, BP_*BP ⊗ V_s) = V_s`, so taking primitives
(`BPComplex`/`primitive_data`) gives a cochain complex of free `BP_*`-modules
whose cohomology is `Ext`. Two filtrations on that one complex, through the
same `SS_table`/Curtis machinery: `I`-adic → `algNov_table`, `p`-adic →
`Boc_table`. `multiplication.cpp` adds products the Atiyah–Hirzebruch way.

## Heights: what the recent work added

`I_n = (p, v_1, …, v_{n-1})` is invariant, so `(BP_*/I_n, BP_*BP/I_n)` is
again a Hopf algebroid and change of rings gives
`Ext_{BP_*BP}(BP_*, M) = Ext_{BP_*BP/I_n}(BP_*/I_n, M)` for `I_n M = 0`. Each
comodule carries a **height** `n` (`ComoduleSpec::height`); `n = 0` is the
classical case and is bit-for-bit unchanged.

Two facts make this cheap, and are worth not re-deriving:

* `(BP_*BP/I_n)/(I/I_n) = P` for **every** `n`, so phase 1,
  `reduce_coaction_rows_mod_I` and the lift are all untouched;
* killing `v_1 … v_{n-1}` is reduction modulo a **monomial** ideal, so reduced
  elements are closed under `+` and `×` — reduce the structure tables and the
  input coaction once and nothing downstream needs overriding.

Where it lives: `Z3_mod_p_Op` (`Z3.h`/`.cpp`) gives `F_p` arithmetic on the
existing `Z3` representation; `BP_Op::height` plus `reduce_{v,t,right,left}_mod_I`
and height-aware `load_etaL`/`load_R2L`/`load_delta` (`BP.cpp`);
`primitive_data::make_primitives` skips killed monomials;
`algNov_table::cycle_pot` drops the `v_0` direction in characteristic p;
`Boc_table::v_valuation` becomes the `v_n`-Bockstein;
`multiplication::vn_extension` replaces the vacuous `a_0` table with
`AANSS_a<n>.txt`; `BP_Op::h0()` has a char-p branch because
`(η_R(v_1) − η_L(v_1))/p` is `0/p` there.

**Do not introduce a second coefficient type.** `polynomial<Fp>` is already
`P`, and `matrix<P>::moduleOper` is a static that both halves of
`mr_BP_comod` share in one process.

## Traps that have already cost time

1. **`BPBP`'s exponent slots are the reverse of the type name.** The *outer*
   exponent indexes the `v_i` via `η_R`; the *inner* one indexes the `t_i`.
   The source of truth is `BP_Op::etaR`'s comment "the right unit, vn is in
   the outer" (`BP.cpp`, search for that phrase — line numbers drift), with
   the full treatment in the `WARNING` block at the top of `comodules.cpp`.
   So `BPBP_opers.monomial(singleVar(1,1), unit(1))` is `η_R(v_1)`, **not**
   `t_1`. Use `BP_Op::h0()`/`t1()`.
2. **That exact confusion was a live bug.** `reduce_BPBP_mod_I` had the slots
   transposed (fixed in `35c6370`): it sent `t_1 ↦ 0` and `η_R(v_1) ↦ v_1`,
   so phase 1 silently received a *split* comodule. Invisible for the sphere
   and for split examples — only a nonzero off-diagonal coaction entry
   exposes it. Probe any change here with the two values `h0()` and
   `η_R(v_1)`.
3. **Nothing checks the comodule axioms, and nothing checks the height** — in
   either direction. Too small and a non-free module is treated as free; too
   large and you get a self-consistent run of `Ext(BP_*, M/I_n)` instead of
   `Ext(BP_*, M)`.
4. **`anss_chart.py` used to assert freeness.** `make_mark` classified every
   `s = 0` summand as `Z_(3)` unconditionally, which drew `S/p`'s
   `Ext^0 = F_3[v_1]` as a row of boxes. The ground ring is now *published* by
   each run in `<prefix>run_info.txt` and read, not assumed — and if that file
   is missing while no class carries a `v0` (which only happens in
   characteristic p) the script refuses to draw rather than guess.
5. **`-c` asymmetry.** With no `-c` the chart prefers the `mr_BP` prefix
   `<t>_BP…` and falls back to `<t>_sphereBP…`; an explicit `-c` never falls
   back.
6. **`.gitignore` vs `example_data/`.** Run output and binaries are ignored by
   pattern, which matches the deliberately committed examples too — hence the
   `!example_data/**` negation. Without it a *newly added* example is
   invisible to `git add`.
7. **The two docs HTML files were rendered by different passes and disagree.**
   `docs/CHARTS.html` joins a list item's soft line breaks into one line and
   leaves `'`/`"` raw; `docs/GENERAL_COMODULES.html` keeps the newlines and
   escapes them as `&#39;`/`&quot;`. Both escape `&<>` and rewrite `.md` links
   to `.html`. When editing an `.html`, render a block you have *not* changed
   first and assert it appears verbatim in the file — that calibrates the
   dialect before you touch anything.
8. **Git in this environment:** pushing works, but `git push --delete` of a
   ref returns **HTTP 403** — the session credentials cannot delete refs.
   Report it rather than retrying.
9. **`prime3` moved twice mid-session.** Always `git fetch` before assuming a
   fast-forward.

## Checking a change

Run all of this. The first item is the one that catches most mistakes.

```sh
# 1. height 0 must be byte-identical. Build the pre-change tree in a worktree
#    and diff every .txt output of mr_BP and of mr_BP_comod for sphere,
#    alpha_1 and triv_01 (22 files at t=20, length=4).

# 2. Ext^0 is the ring of invariants of the base ring, exactly:
#      mod_p    -> F_3[v_1]: filtration 0 at stems 0,4,8,...  and nothing else
#      mod_p_v1 -> F_3[v_2]: stems 0,16,...  and no v1 anywhere in the tables

# 3. the cofiber long exact sequence, bidegree by bidegree
./BPtab 30 && ./mr_BP_comod 30 8 mod_p && ./mr_BP_comod 30 8 alpha_1_mod_p
./ses_check.py 30 -c alpha_1_mod_p -s 30_mod_pBP   # 26 agree, 23 sharp, 0 mismatch
./ses_check.py 30 -c alpha_1                       # 17 agree, 15 sharp, 0 mismatch

# 4. the committed examples must still redraw themselves byte for byte
for c in sphere alpha_1 mod_p; do
  ./anss_chart.py 30 -c $c -d example_data -o /tmp/r.svg
  cmp /tmp/r.svg example_data/30_${c}anss_E2.svg
done

# 5. all four build scripts warning-free with -Wall
```

`example_data/` holds three complete runs at `t=30, length=8` precisely so
that (4) works; keep it complete if you add to it.

When documenting a command or a code sample, **run it verbatim from an empty
directory before committing it.** Doing that caught two wrong claims in this
session's README work.

## Open items

* The remote branch `claude/keen-edison-t29s2m` still exists — it is fully
  merged into `prime3` and should be deleted, but the session could not (403,
  trap 8). Three other stale `claude/*` branches are also behind `prime3`.
* Only the invariant **prime** ideals `I_n` are reachable. `S/p^k` and
  `S/(p, v_1^k)` for `k > 1` do live over invariant quotients, but those are
  truncated polynomial rings and would need a truncated multiplication,
  closer to what `trunc_hopf.cpp` does.
* `v_n`-multiplication lines are not drawn on the ANSS chart (only `α_1` is);
  the `--grading algnov` view does draw them.
* A `les_check.py` was written and passing in-session — the `p`-Bockstein long
  exact sequence check for `S/p` against the sphere (23 bidegrees, 22 sharp,
  0 mismatches) — but never committed. `ses_check.py` cannot do it, since it
  assumes a rank-2 two-cell comodule and `BP_*(S/p)` has rank 1.
* The README's "Charts" section still illustrates with the old
  `mr_st`/`BPtab`/`mr_BP` path.
