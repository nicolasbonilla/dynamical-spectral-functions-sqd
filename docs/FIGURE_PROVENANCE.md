# FIGURE DATA PROVENANCE — single source of truth

**Rule of order:** every figure's numbers come from ONE authoritative data file listed below.
No hand-typed values may drift from these files. Hardware figures carry the exact **job id**.
Before editing any figure, check this table. Before trusting any `.json`, check it is **not** in the
"STALE / DO NOT USE" list at the bottom.

> **This table is DERIVED, and the derivation is checked.** The label column is not maintained by
> hand: `python src/check_provenance.py` expands `\input` and `\includegraphics` from
> `paper/main.tex` and fails if the set of figures documented here is not exactly the set the
> manuscript typesets, or if a cell names a file that is not in the deposit and does not say so.
>
> **Why that check exists.** Until 2026-09-19 this file documented five figures the manuscript had
> stopped typesetting (`fig:amort`, `fig:bench`, `fig:gallery`, `fig:gflow`, `fig:molsuite`) and
> omitted two it had started (`fig:gapscaling`, `fig:thm1iii-violation`). Nothing noticed, because
> nothing derived the table from the manuscript. Run the checker before believing this page.

**Measured state of the manuscript — run `python src/check_provenance.py` for the current census;
do not quote this paragraph instead of running it.** At the time of the C3 deposit (2026-09-19) it
reported 50 source files reached from `paper/main.tex` and **34 floats — 17 figures and 17 tables,
all double-column** — and *failed*, because the preamble census comment of `paper/main.tex` still
read `33 of the 17+16 floats are starred`. *Re-measured 2026-09-26 on the v3 manuscript as
submitted: 60 source files reached (the Supplemental Material is ten further files), 34 floats (17
figures + 17 tables), all starred, and it **passes** — the census comment now reads `34 of the 17+17`.*
*Re-measured 2026-09-28: 61 source files (the new `app_proof.tex`, Appendix A), the same 34 floats; it passes.*

Paths: figure sources live in `paper/figs/` and their float wrappers in `paper/figs/float_*.tex`,
`paper/figs/captions.tex` and `paper/carried/`; generating scripts in `src/`, data in `data/`.

**Known gaps between what a script does and what the deposit contains are listed in
[`KNOWN_DISCREPANCIES.md`](KNOWN_DISCREPANCIES.md). Read it before re-running any compute script.**

---

## The figure table

`#` is the figure number **printed in v3 as submitted** (1–9 main text, S1–S8 Supplemental
Material); the `fig:` label is the stable key and the column `src/check_provenance.py` checks.
Documents written between 2026-09-19 and 2026-09-25 number the figures 1–17 continuously; the
conversion table opens [`KNOWN_DISCREPANCIES.md`](KNOWN_DISCREPANCIES.md).
`enters as` — **PDF** = `\includegraphics`, **frag** = a native pgfplots fragment compiled with the
paper. `guarded` — regenerated and diffed byte-for-byte by `src/check_figures.py`.
A cell that names a file which is not in this deposit says **NOT DEPOSITED** in that cell.

| # | Body ref | Assembled in | Figure body | enters as | Generator | AUTHORITATIVE data | guarded |
|---|----------|--------------|-------------|-----------|-----------|--------------------|---------|
| 1 | fig:circ | `paper/sec_2.tex` | `fig_circuit.pdf` | PDF | hand-drawn quantikz; standalone source `paper/figs/src/fig_circuit_native_FIXED.tex` (deposited 2026-09-26; recompiles pixel-identically, not script-guarded) | — (schematic; plots no data) | no |
| 2 | fig:method | `paper/carried/fig_method.tex` | `fig_method_native_frag.tex` + `aw_method.dat` | frag | `make_method_fig_max.py` | **`method_max.json`** | **yes** |
| 3 | fig:akwsampled | `paper/figs/float_f5.tex` | `fig_akw_sampled_honest.tex` + `fig5_caption.tex` | frag | `make_akw_sampled_honest_fig.py` (paths from `__file__` since the §A repair); the caption's finite-shot paragraph from `fig3_finite_shot.py` | **`sampled_akw_L8.json`** (panels a,b); panel (c) from `sampled_honest.json` + `honest_sampling.json`; caption finite-shot numbers from `fig3_finite_shot_T2.6e6.json` + `fig3_finite_shot_T5.2e6.json` | **yes** |
| 4 | fig:lattice | `paper/carried/fig_lattice.tex` | `fig_akw_native.pdf` | PDF | standalone source `paper/figs/src/fig_akw_fixed.tex` with its raster `paper/figs/src/akw_field_v2.png` (deposited 2026-09-26; recompiles pixel-identically). That raster is the output of `recolor_akw_v2.py`, whose own inputs (the viridis raster and `fig_akw_v2_L12.tex`) are **NOT DEPOSITED**, so the recolouring cannot run on a clean clone | **`akw_lanczos_L12.json`** | no |
| 5 | fig:scaling | `paper/figs/float_f12.tex` | `fig_scaling2_native.tex` + `caption_fig12.tex` | frag | ⚠ **none usable** — `make_scaling_fig.py` emits a DIFFERENT two-panel figure and would destroy panel (c); §3 of `KNOWN_DISCREPANCIES.md` | **`scaling_data.json`** | no |
| 6 | fig:gapscaling | `paper/figs/float_n2.tex` | `fig_gapscaling_native.tex` + `fig_gapscaling_caption.tex` | frag | `src/recovered/make_fig_gapscaling.py` (recovered 2026-09-26). The fragment regenerates except one provenance-timestamp comment line; the caption was edited by hand after generation and is the authoritative one | **`gap_scaling.json`** (computed by `gap_scaling.py`) | **yes** (fragment; one comment line normalised) |
| 7 | fig:sqw | `paper/carried/fig_sqw.tex` | `fig_sqw.pdf` | PDF | standalone source `paper/figs/src/fig_sqw_fixed.tex` (deposited 2026-09-26; recompiles pixel-identically). Its raster `sqw_field.png` is deposited as rendered: `sqw_field.py` draws the same field in viridis, and the recolouring to hot is a step that is not deposited | **`sqw_L12.json`** + `sqw_edges.dat` (by `make_sqw_edges.py`) | no |
| 8 | fig:spin | `paper/carried/fig_spin.tex` | `fig_spinqw.pdf` | PDF | standalone source `paper/figs/src/fig_spinqw_fixed.tex` (deposited 2026-09-26; recompiles pixel-identically); no deposited script writes its raster `paper/figs/src/spinqw_field.png`; the caption's deviations of the plotted peak from the des Cloizeaux–Pearson boundary (about 3 %, 10 %, 29 %) are printed by `src/spin_peak_deviation.py` → `data/spin_peak_deviation.json` (2026-09-28, late) and checked by `verify.py` | **`spinqw_L12.json`** + `spinqw_edges.dat` (by `spin_lanczos.py`) | no |
| 9 | fig:heron | `paper/carried/fig_heron.tex` | `fig_hardware_hero_frag.tex` + `heron_hot.dat` | frag | `make_hardware_hero.py` | ⭐ **`hw_lucj_n2_result.json` — job `da125f2ein7c73bcsqs0`** (panel b) + **`heron_spectral.json`** (panel a) | **yes** |
| S1 | fig:witness | `paper/sm_nogo.tex` | `fig_witness_native.tex` | frag | none — hand-maintained against the JSON (whose `chi_max`, the largest Schmidt rank over all cuts of the pairing ordering, is the χ = 2 the figure and the theorem print; recorded by `apsg_witness.py` since 2026-09-28, late, and recomputed by `verify.py`) | **`apsg_witness.json`** | no |
| S2 | fig:decoupling | `paper/figs/captions.tex` | `fig_decoupling_native.tex` | frag | `make_decoupling_native.py` | **`cost_vs_ent.json`** + `n19_suite.json` + `resource_master.json` + `stats_resource.json` | **yes** |
| S3 | fig:master | `paper/figs/captions.tex` | `fig_resource_master_native.tex` | frag | `make_decoupling_native.py` (emits both resource fragments). `make_resource_master_fig.py` writes the same file — **two writers, one artefact**; the guarded one is the first | **`resource_master.json`** + `stats_resource.json` | **yes** |
| S4 | fig:ladder | `paper/carried/fig_ladder.tex` | `fig_ladder_native.tex` | frag | none — hand-maintained against the JSON | **`ladder_vs_chain.json`** | no |
| S5 | fig:thm1iii-violation | `paper/figs/float_n3.tex` | `fig_thm1iii_violation_native.tex` + `caption_thm1iii_violation.tex` + eight `n3_*.dat` | frag | `src/recovered/build_n3.py` (recovered 2026-09-26) writes the ten `n3_*.dat`; `src/recovered/n3_fig5_whiskers.py` (2026-09-28, late) then writes `n3_fig5_summary.json` and the two whisker blocks of the fragment; `fix_n3_caption.py` and `fix_n3_fig5_series.py` are deposited as the record and are not run | the deposited **`data/cert_stress.json`, `cert_akw.json`, `cert_teqsci.json`** (since 2026-09-28, late; until then the earlier run of the same scans in `data/thm1iii_violation/cert_earlier_run/`, which stays deposited as a record and is no longer read) + `data/thm1iii_violation/thm3_results.json` + `thm3_L8b.json` | **tables, whiskers and summary yes** (`check_figures.py`); the plotted ratios and every printed count re-derived from the certificate rows by `verify.py` section 9.6c, §B |
| S6 | fig:hero | `paper/carried/fig_hero.tex` | `fig_hero.pdf` | PDF | `make_hero_fig.py` writes a `fig_hero.tex` into the **cwd** (that file is **NOT DEPOSITED**; measured 2026-09-26, it differs from the source of the committed PDF in one line, the position of the `R (Å)` label, adjusted by hand afterwards); the standalone source of the committed PDF is `paper/figs/src/fig_hero_STANDALONE.tex` (deposited 2026-09-26; recompiles pixel-identically) | **`n2_hero.json`** + `hero_aw.dat` + `hero_res.dat` | no |
| S7 | fig:noise | `paper/carried/fig_noise.tex` | `fig_noise_score.pdf` | PDF | `make_noise_fig.py` writes a `fig_noise_native.tex` into the **cwd** (that file is **NOT DEPOSITED**; measured 2026-09-26, it differs from the source of the committed PDF in two lines, the length of the ε=16% marker and the size and position of the left annotation, adjusted by hand afterwards with a comment saying why); the standalone source of the committed PDF is `paper/figs/src/fig_noise_score_STANDALONE.tex` (deposited 2026-09-26; recompiles pixel-identically) | **`noise_spectral.json`** + `noise.dat` | no |
| S8 | fig:noiserec | `paper/carried/fig_noiserec.tex` | `fig_noise_recovery_native.tex` | frag | `make_noise_recovery_native.py` | **`molecular_noise.json`** + `molecular_noise_sweep.json` | **yes** |

### Column "guarded", counted

`src/check_figures.py` regenerates **19 artefacts** covering **8 of the 17 figures**
(2, 3, 6, 9, S2, S3, S5 — its ten tables — and S8). The other nine have no byte-for-byte guard:
four have no deposited generator at all (1, 5, S1, S4), and five ship as a PDF (4, 7, 8, S6, S7).
Fig. 1 is also a PDF; it is counted once, among the four without a generator.

*Until 2026-09-26 the count was 8 artefacts and 6 figures (2, 3, 9, S2, S3, S8): the generators of
Figs. 6 and S5 were not in the deposit, and none of the six PDF figures had its standalone source
here. Since that day the generators are in `src/recovered/` and the six sources in `paper/figs/src/`,
where each recompiles to a pixel-identical rendering of the committed PDF — but no script guards
those six, so they are not counted as guarded (`KNOWN_DISCREPANCIES.md` §29).*

---

## §A — the one guard that did not fire, measured — **REPAIRED 2026-09-19**

> **Closed.** Both literals below are gone: `REPO` is now `os.path.dirname(HERE)` and the two
> JSONs are read from the deposited `data/`. `check_figures.py` guards 8 of 8 artefacts
> (sentinel test, positive control) and writes nothing into the repository (sha256 + mtime
> snapshot of 231 files, before and after). What follows is the record of the defect.

`src/make_akw_sampled_honest_fig.py` lines 17–18 set

    REPO = r"C:\Users\Nicolas\Downloads\Proyecto_SQD_ML\P2_dynamical_spectral_functions\release"
    CALC = r"C:\Users\...\scratchpad\p2_calc"

Both are absolute. Two consequences, both measured on 2026-09-19 rather than reasoned about:

1. **It cannot run on a clean clone.** `CALC` is a session scratch directory that is not part of
   this deposit, and panel (c) is read from it — although `data/sampled_honest.json` and
   `data/honest_sampling.json` are deposited, the script does not read those copies.
2. **`src/check_figures.py` cannot see a corrupted Fig. 3.** The checker copies `src/` and `data/`
   into a scratch tree and runs the generators there; this one ignores the scratch tree and writes
   through `REPO` into the real `paper/figs/`. Negative control: corrupting one plotted literal in
   the committed `fig_akw_sampled_honest.tex` and re-running the checker on that tree gives
   **exit 0, "all 8 regenerated artefacts are byte-identical"**. The same corruption applied to
   `fig_decoupling_native.tex`, whose generator resolves paths from the cwd, gives **exit 1 and
   `[DIFFER]`**. So the positive control fires and this one does not.
   A side effect of the same bug: running `make check-figures` **rewrites four files in the working
   tree** (`fig_akw_sampled_honest.tex`, its two caption files and
   `data/sampled_akw_figure_numbers.json`) with identical bytes, which contradicts the checker's
   own docstring ("It never writes into the repository").

Until those two constants are resolved from `__file__` and `data/`, read the checker's result as
**8 of 8 guarded since the §25 repair of 2026-09-19** — before it, 7 of 8. Measured by seeding each artefact with a sentinel line and checking it is gone after the run, with a positive control asserting the sentinel was in all eight first.

## §B — the figure the guardian declares it cannot anchor

`src/verify.py` registers `fig_thm1iii_violation_native.tex` in `COVERAGE_GAP_V3` and prints it as
`[XFAIL]` on every run: the fragment plots eight `n3_*.dat` tables (2432 lines), and no check in
`verify.py` anchors those ratios to the certificate rows. This is printed, not hidden, and it is the
honest status of Fig. S5 (`fig:thm1iii-violation`). *Since 2026-09-26 the builder and its inputs are
deposited (`src/recovered/build_n3.py`, `data/thm1iii_violation/`) and `check_figures.py`
regenerates all ten tables byte-for-byte; the coverage gap stays, because the check it names has
still not been written.* What remained to close it was that check: a test in `verify.py` that anchors
the plotted ratios to the certificate rows. *Closed 2026-09-28 (late): `verify.py` section 9.6c
re-classifies the pool from the five deposited JSON sources (`thm3_results.json`, `thm3_L8b.json`,
`data/cert_stress.json`, `cert_akw.json`, `cert_teqsci.json`) with code that shares nothing with
`build_n3.py`, requires every point of the eight plotted tables to be the error-to-bound ratio of a
classified row (to six significant figures), checks the whisker summary against the 64 rows of the
Fig. 3 setting, recomputes the 3.1e-8 control of those rows against `sampled_akw_L8.json`, and checks
every count the caption, Sec. III D, Sec. S2 and Sec. S9 print. The entry is removed from
`COVERAGE_GAP_V3` with a dated comment; the guardian now prints four items by name, the known open
defects. The fragment is diffed by `check_figures.py` since the same day (its whisker blocks are
written by `n3_fig5_whiskers.py`), and `build_n3.py`'s self-check still requires its legend counts to
match the tables.*

---

## ⭐ HARDWARE — the single authoritative run (do not confuse)

**N₂ / ibm_marrakesh energy = job `da125f2ein7c73bcsqs0`** (8-seed recovery sweep, "Run all").
File: `hw_lucj_n2_result.json`. The **full per-seed × per-iteration trajectories** are stored in
`e_hw_seed_hist` (8×5 energies) and `dE_hist_seeds` (8×5 mHa), with `dE_hist_mean`/`dE_hist_std` the
per-step 8-seed statistics — reconstructed from the run log so the per-point error bars are never lost
again. Converged dE = **0.59 ± 0.14 mHa** (best seed 0.40); noiseless sim 29.5 mHa.
Confirmed authoritative by Nicolás (2026-08-17); per-seed trajectories restored 2026-08-17.

**A(ω) / ibm_fez** = the `ibm_fez` L=6 Hubbard run (|S|=300/300, exact by coverage). File:
`heron_spectral.json` → `figs/heron_hot.dat`. The 300 determinants were retained by post-selection in
reversed bit order, so every one of them came from a device error; no configuration recovery and no
readout-error mitigation were applied (measurement twirling, Pauli twirling and dynamical decoupling
only) — `KNOWN_DISCREPANCIES.md` §30.

---

## Sec. V has no figure — and, since 2026-09-19, it has a source

This file is a figure map, and Sec. V typesets no figure, so nothing about it ever appeared here.
That silence was doing real damage: **Tables III, IV and V, the frontier, the calibration and the
convergence census had no deposited source at all**, and `sec_5_body.tex` said so in its own
provenance comment — *"this section is still the one part of the paper a reader cannot regenerate"*.

It is closed. The 204-row sweep those tables summarise is in **`data/c3_frontier/`** and the code
that produced it in **`src/frontier/`**, documented file by file in
[`C3_FRONTIER.md`](C3_FRONTIER.md). The map for Sec. V, in the same form as the figure table above:

| Sec. V object | authoritative file | producer |
|---|---|---|
| Tables III–V, the frontier, the calibration, the 204 evaluations | `data/c3_frontier/C3_grid.json` | `src/frontier/fsum.py`, from the 54 `data/c3_frontier/pt_L*_FR*_nl*.json` |
| the convergence census (187 / 6 / 11) | same | classifier in `fsum.conv_flag`, re-derived in `verify.py` §(10b) |
| the zero-Born / tie audit of the cut | `data/c3_frontier/ties.json` | `src/frontier/ftie.py` |
| the cross-validation against a dense Λ_S | `data/c3_frontier/validation.json`, `validation_full.json` | `src/frontier/fvalid.py` |
| the proof-ladder displacement (how much of the frontier is Cauchy–Schwarz) | `data/c3_frontier/verdict_frontier.json` ← `ladder_L{4,6,8}.json` | `src/frontier/verdict.py` ← `run_ladder.py`, `fron_lib.py` |
| the rows at the fractions the manuscript publishes, and L=14 | `data/c3_frontier/published_fraction/` | `src/frontier/pubsum.py`, `run_L14.py` |
| the four adversarial re-runs named in `sec_5_body.tex` | `data/c3_frontier/adversarial/` | `src/frontier/z_{fine,rank,e0,suppK}.py` |
| the null a referee will propose (`1−w_S`) and the out-of-sample calibration | `data/c3_frontier/null/` | `src/frontier/null_ws4-7.py`, `calib_oos.py` |

`src/check_provenance.py` does not police this table: it validates `fig:` labels against the
manuscript, and there are none here. What polices it is `src/verify.py` §(10b), which opens the 54
point files and re-derives all 204 rows — 2 387 assertions, seven negative controls.

---

## The exceptions to "every plotted value is deposited"

The Data Availability Statement ([`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md)) points at this file
as the figure → file map. That map resolves for sixteen of the seventeen figures. It does **not**
resolve for these, which are therefore named here rather than left implicit *(since 2026-09-26 only
Fig. 1 remains, see the second row)*:

| figure | what is not in `data/` |
|---|---|
| **Fig. 1** (fig:circ) | a circuit schematic. It plots no data at all, so there is nothing to deposit. |
| **Fig. S5** (fig:thm1iii-violation) | ~~two of the four JSONs behind the plotted ratios (`thm3_results.json`, `n3_fig5_summary.json`) are not deposited, and neither is the script that turned them into the `n3_*.dat` tables~~ — **resolved 2026-09-26**: all of them are in `data/thm1iii_violation/` and the script is `src/recovered/build_n3.py`. See §B. |

*Withdrawn from this table 2026-09-19: a row reading "Fig. 6, panel (c): the four free-fermion
support sizes 35, 336, 3496, 37361 are constants typed into `src/make_decoupling_native.py` and
written from there into the fragment." **That was a v2 row renumbered from Fig. 1 to Fig. 6
without being re-measured, and it had stopped being true.** Panel (c) of `fig_decoupling_native.tex`
was rebuilt for v3 as "(c) Hubbard: χ runs forwards"; a byte-level grep over every `.tex`, `.dat`
and `.py` under `paper/` and `src/` finds `37361` in exactly one place, the **prose** of
`paper/sec_6_body.tex:337`, and `make_decoupling_native.py` no longer contains it. Those four
numbers were still hand-typed constants with no deposited producer — that limit was real and is
recorded in `KNOWN_DISCREPANCIES.md` — but it is a sentence in Sec. VI, not a figure exception,
and printing it here would advertise an exception for a panel that no longer exists. *(2026-09-28,
late: they have a producer, `src/free_fermion_support.py` → `data/free_fermion_support.json` —
open chain, half filling, the smallest determinant set carrying all but ε_w = 2.5e-3 of the weight,
inclusive at the exact four-way tie of L = 4 — and `verify.py` recomputes the four counts
independently and checks the two sentences, Secs. VI and S13, that print them; §33.)*
`src/check_provenance.py` did not catch this because it validates the label SET against
`main.tex`, not the CONTENTS of a cell.*

---

## Orphans in `paper/figs/` — files with no producer and no consumer

Derived by `python src/check_provenance.py -v`, 2026-09-19; counts re-measured 2026-09-26, after
`src/recovered/` was deposited and added to the scan (10 / 46 / 6). The tiers are defined by two facts
that need no guessing: does any deposited file *name* the artefact, and does the *manuscript reach it*.

**Tier 1 — no deposited *code* names it and the manuscript does not use it (10).** Being written
down here does not give a file a producer, so documenting an orphan does not remove it from this
list:

- *(Left tier 1 on 2026-09-26: `paper/figs/n3_fig5_eta.dat` and `n3_fig5_frac.dat`, the 64 Fig.-3 rows
  that were summarised into the whiskers of Fig. S5. `src/recovered/build_n3.py` now writes them
  and `src/recovered/fix_n3_fig5_series.py` is the record of the summary; the whisker values were
  still literals in the fragment until 2026-09-28 (late), when `src/recovered/n3_fig5_whiskers.py`
  took over writing them, guarded by `check_figures.py`.)*
- `data/cert_frontier_6-8.json`, `data/cert_frontier_10.json` — written by
  `leakage_certificate_suite.py` under a *computed* filename
  (`"cert_frontier_%s.json" % "-".join(...)`), so the literal names appear in no script. Only a
  comment in `paper/sec_5_body.tex` mentions them, and nothing reads them back.
- the eight `docs/img/*.png` README thumbnails. Since 2026-09-26 they are rendered from the v3
  figures by `docs/img/make_thumbnails.py` (figure body only, no caption); that script lives next
  to them rather than in `src/`, so this checker does not count it as a producer. The ten earlier
  thumbnails, hand-cropped from the **v1** PDF — two of them (`bench.png`, `gallery.png`) of figures
  the manuscript withdrew, one (`scaling.png`) a block of withdrawn v1 text, one (`sqw.png`) showing
  the two fabricated `sqw_edges.dat` rows — are in `_superseded/v1_readme_img/`, with a README
  saying what each showed.

**Tier 3 — the manuscript typesets it and no deposited code names it (6):** the six PDFs
(`fig_akw_native.pdf`, `fig_circuit.pdf`, `fig_hero.pdf`, `fig_noise_score.pdf`, `fig_spinqw.pdf`,
`fig_sqw.pdf`). Their standalone LaTeX sources are in `paper/figs/src/` since 2026-09-26 and recompile
to pixel-identical renderings, but LaTeX is not code in this checker's sense, so they stay here. *(Until
2026-09-26 this tier also held seven of the eight `n3_*.dat` tables Fig. S5 plots; `build_n3.py` names
them now.)*

Tier 2 (produced by code, unused by the paper — 46 files, mostly datasets whose only remaining
reader is `src/verify.py`) is printed by the checker and itemised in
[`FILE_INDEX.md`](FILE_INDEX.md).

---

## ⛔ STALE / DO NOT USE (archive or ignore)

- **`notebooks/HW_LUCJ_N2_Heron_READY.ipynb` SAVED CELL OUTPUTS** = an OLDER run, job `d9qe731dsedc73af67d0` (single seed, final 1.2 mHa, noiseless 28.5). These are NOT the paper's numbers. The authoritative run is `da125…` in `hw_lucj_n2_result.json`. Re-run "Run all" and re-save the notebook to sync, or ignore the saved outputs.
- **The `rebuild/` tree** is a scratch build directory (git-ignored) holding `make_table.py`'s output. The **shipped** table is `paper/table_molecules_v3.tex`.
- `scaling_data_dense_backup.json`, `sqw_L12_eta020_backup.json` — backups, not the plotted data (git-ignored by `*_backup.json`).
- **`paper/figs/fig_akw_sampled_honest_caption.tex` and `..._caption_macro.tex` — SUPERSEDED.**
  They are written by `make_akw_sampled_honest_fig.py`, nothing reads them, and their finite-shot
  penalties (`1.46×` at L=6, `2.05×` at L=8) predate the ones the paper prints. The caption the
  manuscript typesets is `figs/fig5_caption.tex` (`1.50±0.25`, `1.96±0.12`, `2.05±0.07`), and those
  three numbers are exactly the `ratio` row of Table IX. They also still open with *“Declared gap: no
  finite-shot run exists for the configuration of (a)–(b)”*, which stopped being true on 2026-09-19,
  when that configuration was run at finite shots ($5.2×10⁶$ shots per channel; generator `src/fig3_finite_shot.py`, deposited 2026-09-28) and both Sec. 2 and
  `fig5_caption.tex` replaced the declaration with the measured budget. **Substituting the
  generator's caption because it is machine-written would reinstate withdrawn numbers and a gap that
  is now closed.** `KNOWN_DISCREPANCIES.md` §26.
- **`paper/figs/fig_amort_native.tex` — RETIRED.** The body of `fig:amort`, which Sec. 9.4 withdraws.
  It stays in `paper/figs/` only because `src/verify.py` names it; nothing typesets it.
- `figs/scaling.dat`, `figs/conv_method.dat`, `figs/resource_axis.dat` — all three **DELETED 2026-09-18** as stale orphans; rationale in `KNOWN_DISCREPANCIES.md` §8. Do not restore them from `_superseded/` or from the arXiv v2 bundle, which still contains two of them.
- **`_superseded/`** holds the v2-era manuscript as it stood in `paper/` on 2026-09-19 (`v2_paper/`; edited after the 2026-08-22 posting, so it is not the posted v2 source, which is `paper/arxiv-submission.tar.gz`), the v1 README thumbnails (`v1_readme_img/`), the retired v2 figure sources (`v2_figs/`) and the bodies of the five figures withdrawn in v3 (`v3_retired/figs/`: `fig_benchmark.pdf`, `fig_gflownet.pdf`, `fig_molecular_gallery.pdf`, `fig_molecular_suite.pdf`, `fig_hardware_hero.pdf`, `heron_spectral.pdf`, `_method_standalone.*`, `bench_U*.dat`). **RETIRED — nothing in `paper/` reads them.** Their generators (`make_bench_fig.py`, `make_gflow_fig.py`, `make_gallery_native.py`, `make_magic_suite_fig.py`) are still in `src/` and would recreate `paper/figs/bench_U*.dat` if run; that is the only reason to open that directory.
- **`paper/arxiv-submission.tar.gz`** is the frozen **arXiv v2** bundle and must never be re-uploaded as v3: it carries the seven v2 section files, the five retired figures and the two fabricated rows of `figs/sqw_edges.dat`. Run `python src/check_tarball.py` for the live comparison.
