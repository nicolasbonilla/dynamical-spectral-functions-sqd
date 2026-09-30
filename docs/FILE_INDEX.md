# FILE INDEX — every file in this repository, described

*Rebuilt 2026-09-19 against the v3 REVTeX manuscript; updated 2026-09-26 for v3 as submitted to arXiv
on 2026-09-25 (one column, 66 pp: a 26-page main text, the reference list, then the Supplemental
Material; "25-page" here until 2026-09-28, which miscounted by one) and for the recovery of that day (`KNOWN_DISCREPANCIES.md` §29). The version before
2026-09-19 described the **v2** file set in its entirety — `introduction.tex`, `method.tex`,
`results.tex`, `discussion.tex`, `conclusion.tex`, `hardware.tex`, `theorem_b2.tex`, twenty figures, a
40-page PDF. None of those files has been in `paper/` since commit `606f969`. See
[`KNOWN_DISCREPANCIES.md`](KNOWN_DISCREPANCIES.md) §28.*

Grouped by role. For *how to run* each script, see [`REPRODUCE.md`](REPRODUCE.md); for *which data feeds
which figure*, see [`FIGURE_PROVENANCE.md`](FIGURE_PROVENANCE.md); for *where the deposit and the code
disagree*, see [`KNOWN_DISCREPANCIES.md`](KNOWN_DISCREPANCIES.md).

Figure and table numbers below are those **printed in v3 as submitted**: Figs. 1–9 and Tables I–V in the
main text, Figs. S1–S8, Tables S1–S12 and Secs. S1–S16 in the Supplemental Material. The conversion from
the numbering of earlier builds opens `KNOWN_DISCREPANCIES.md`.

**Counted, not asserted.** 425 files tracked by git after the commits of 2026-09-26 (389 at
`f83c49f`, the GitHub head of 2026-09-25). `rebuild/`, `build/`, `__pycache__/`, `*.log`, the
`paper/main.{aux,log,out,bbl,blg}` build files and every `ckpt/` file except the three deposited ones
are present on disk and `.gitignore`d; they are build scratch and are described at the end rather than
listed.

| directory | tracked files |
|---|---|
| `paper/` | 101 — 33 at the root, 9 in `carried/`, 49 in `figs/` and 10 in `figs/src/` |
| `src/` | 104 — **66** at the root, **29 in `src/frontier/`** (the C3 sweep, deposited 2026-09-19) and **9 in `src/recovered/`** (recovered 2026-09-26) *(2026-09-28, late: `free_fermion_support.py`, `spin_peak_deviation.py`, `frontier/rerank_L14.py`, `recovered/n3_fig5_whiskers.py`)* |
| `data/` | 166 — **43** at the root, **96 in `data/c3_frontier/`**, **20 in `data/gflow_runs/`** and **7 in `data/thm1iii_violation/`** *(2026-09-28, late: `free_fermion_support.json`, `spin_peak_deviation.json`, `c3_frontier/published_fraction/L14_rerank.json`)* |
| `docs/` | 16 — 7 Markdown, 8 PNG thumbnails and the script that renders them |
| `ckpt/` | 3 — the L=14 ground state, its ranking and their README (deposited 2026-09-25) |
| `_superseded/` | 36 — `v2_paper/` 11, `v2_figs/` 2, `v3_retired/` 11, `v1_readme_img/` 11 (ten PNG and a README), and the folder's own `README.md` |
| `notebooks/` | 3 |
| root | 8 — `README.md`, `CITATION.cff`, `.zenodo.json`, `LICENSE`, `Makefile`, `requirements.txt`, `.gitignore`, `.github/workflows/ci.yml` |

---

## `paper/` — the v3 manuscript (REVTeX 4.2, one column for submission)

`main.tex` is the only hand-written driver. Everything it reaches is listed here, and the list is
derived: `python src/check_provenance.py` expands `\input` from `main.tex` and reports **61 source
files reached** and **34 floats — 17 figures and 17 tables, all starred** (re-measured 2026-09-28; 60 files on 2026-09-26,
before `app_proof.tex`).

### Main text

| File | Contents |
|---|---|
| `main.tex` | Preamble (REVTeX options, the `figladder` font restoration, float placement, theorem environments, hyperref), title, author block, abstract, and the `\input` list. The author block carries three affiliations (Universidad Nacional de Colombia first), set 2026-09-19. |
| `sec_1.tex` | **I. The operational question** (`sec:intro`). |
| `sec_2.tex` | **II. One sampling primitive, four channels** (`sec:method`) — holds Table I (`tab:provenance`) and Fig. 1 (`fig:circ`); Figs. 2–4 follow it from `main.tex`. |
| `sec_3_body.tex` | **III. Two budgets: the missed weight and the boundary** (`sec:leakage`) — the central theorem. |
| `sec_3_1.tex` | **III D. The captured weight is not the controlling magnitude** (`sec:retraction`), `\input` from inside Sec. III. |
| `sec_4_body.tex` | **IV. What the structure does certify: moment exactness** (`sec:moments`) — Table II, and Proposition 1 (`prop:shots`, sampling capture; in the SM until 2026-09-28). |
| `sec_5_body.tex` | **V. Where the leakage certificate is informative: its non-vacuity threshold, measured** (`sec:frontier`) — Table III. |
| `sec_6_body.tex` | **VI. No shortcut to the boundary: an obstruction and its witness** (`sec:nogo`). |
| `sec_7.tex` | **VII. Three prices: support, resolution, and shots** (`sec:prices`) — Table IV; Fig. 5 follows it. |
| `sec_8.tex` | **VIII. Physical validation: the two gaps scale in opposite directions** (`sec:physics`) — Table V; Figs. 6–8 follow it. |
| `sec_9.tex` | **IX. Hardware, scope, and the region where the method has no advantage** (`sec:scope`, `sec:honest`) — Fig. 9, and the data availability statement as Sec. IX D. Its header records that `fig:gflow` and `fig:amort` are deliberately absent. |
| `sec_10.tex` | **X. Conclusion** (`sec:conclusion`), and the acknowledgments. |
| `app_proof.tex` | **Appendix A. Proof of the leakage certificate and of its closed form** (`app:proof`), added 2026-09-28: items (1)–(5) of the proof of Theorem 1, moved verbatim from `sm_leakage.tex`, with Eq. (A1) (`eq:split`). `main.tex` inputs it after `sec_10.tex` and before the references, inside a group that issues `\appendix` (the one global effect, the per-section equation reset, is undone in `sm_begin.tex`). |
| `bibliography.tex` | Hand-written `thebibliography` — 142 entries, arXiv-safe, no external `.bib`/`.bbl`. In first-citation order (restored 2026-09-30, when six keys moved into Sec. I and Sec. I began citing `sierant2026`): the 103 references of the main text, among them `supplemental` (ref. [19]; ref. [12] from 2026-09-28 to the reordering: the Supplemental Material cited as a reference in the APS form, "See Supplemental Material at [URL will be inserted by publisher] for …, which includes Refs. [104–142]", first cited where the SM is first mentioned, at the end of Sec. I, and again at the end of Sec. X; it was ref. [102], cited only in Sec. X, until the review pass of 2026-09-28), then the 39 works cited only in the SM (the range is computed with `\citenum`, and `pra_split.py --build` checks it). |

### Supplemental Material

`sm_begin.tex` opens it (numbering restarts as S1, S2, …). The last eight files are the appendices of the
2026-09-19 build, now Secs. S9–S16.

| File | Contents |
|---|---|
| `sm_begin.tex` | opens the Supplemental Material; no section of its own. |
| `sm_method.tex` | **S1.** The primitive: provenance declarations in full. |
| `sm_leakage.tex` | **S2.** The leakage certificate and the weight-only bound: details (the proof of Theorem 1 moved to Appendix A, `app_proof.tex`, on 2026-09-28; the (H1)–(H3) list it repeated is kept once, in S11). |
| `sm_moments.tex` | **S3.** Moment exactness: the proof route, the conjecture and the shot bound. |
| `sm_frontier.tex` | **S4.** The non-vacuity threshold: protocol, tables and limits — Tables S1–S3. |
| `sm_nogo.tex` | **S5.** One-body magic, the geminal witness and the determinant support — Fig. S1, Table S4, Figs. S2–S3 (via `figs/captions.tex`) and Fig. S4. |
| `sm_prices.tex` | **S6.** The three prices: the resolution grid and finite shots — Tables S5–S6. |
| `sm_physics.tex` | **S7.** Physical validation: spin scaling, pole census, sum rule and size — Table S7. |
| `sm_scope.tex` | **S8.** Hardware and scope: full account (job identifiers, shot arithmetic, the "executed on" qualification). |
| `sec_3_app.tex` | **S9.** Counterexamples to the weight-only bound, and the pooled sweep. Fig. S5 (`figs/float_n3.tex`) follows it. |
| `sec_4_app.tex` | **S10.** Moment exactness: statement, proof, sharpness, and the shot count. |
| `sec_5_app.tex` | **S11.** The leakage certificate: hypotheses, protocol, and the anatomy of the threshold. |
| `sec_6_app.tex` | **S12.** The geminal witness: proof and angle sweep — Table S8. |
| `app_carried_theorem.tex` | **S13.** One-body magic cannot bound the determinant support (carried from v2's `theorem_b2.tex`). |
| `app_stats.tex` | **S14.** The resource statistics, stratified — Table S9. |
| `app_carried_repro.tex` | **S15.** Reproducibility data for the molecular suite — Tables S10 (`tab:repro`), S11 and S12, and Figs. S6–S8. **Corrected 2026-09-19:** `R_eq(H₂O)` now prints `0.958`, the value the geometry in `src/n19_suite.py:123` gives (`√(0.757² + 0.587²) = 0.957924 Å`); it printed `0.957`. |
| `app_dequant.tex` | **S16.** The dequantization test (`app:dequant`); its numbers come from `data/gflow.json`, written by `src/gflow_dequant.py`. |
| `table_molecules_v3.tex` | Table S12 (`tab:molecules`), the nineteen-molecule suite, `\input` from `app_carried_repro.tex`. The v2 `table_molecules.tex` is in `_superseded/v2_paper/`. |

### `paper/carried/` — float wrappers carried verbatim from the v2 tree

Nine one-float files, each `\begin{figure*}` + `\input` of a fragment (or `\includegraphics`) +
caption: `fig_method.tex`, `fig_lattice.tex`, `fig_ladder.tex`, `fig_sqw.tex`, `fig_spin.tex`,
`fig_hero.tex`, `fig_heron.tex`, `fig_noise.tex`, `fig_noiserec.tex`. They exist because v3 keeps the
v2 captions unchanged for those nine figures; the float and the `\input` share a line, which is why a
line-level `\input` expander finds two floats in this manuscript instead of thirty-three.

`carried/fig_amort.tex` and `carried/fig_gflow.tex` do **not** exist: their `\input` lines in
`main.tex` are commented out, with the reason (both figures are withdrawn in v3).

### `paper/figs/` — 49 files, plus `src/` (10)

| group | files | note |
|---|---|---|
| float wrappers | `float_f5.tex`, `float_f12.tex`, `float_n2.tex`, `float_n3.tex`, `captions.tex` | hand-written; `captions.tex` carries two floats (`fig:decoupling`, `fig:master`) |
| captions | `fig5_caption.tex`, `caption_fig12.tex`, `fig_gapscaling_caption.tex`, `caption_thm1iii_violation.tex` | the numbers in each are traced in `FIGURE_PROVENANCE.md` |
| native fragments (typeset) | `fig_method_native_frag.tex`, `fig_akw_sampled_honest.tex`, `fig_witness_native.tex`, `fig_decoupling_native.tex`, `fig_resource_master_native.tex`, `fig_ladder_native.tex`, `fig_scaling2_native.tex`, `fig_gapscaling_native.tex`, `fig_hardware_hero_frag.tex`, `fig_noise_recovery_native.tex`, `fig_thm1iii_violation_native.tex` | 11 fragments; `src/verify.py` opens all 11 and declares one coverage gap |
| compiled figure bodies | `fig_circuit.pdf`, `fig_akw_native.pdf`, `fig_sqw.pdf`, `fig_spinqw.pdf`, `fig_hero.pdf`, `fig_noise_score.pdf` | 6 PDFs. Their standalone sources are in `src/` below since 2026-09-26 (`KNOWN_DISCREPANCIES.md` §5, §29) |
| plotted data | 19 `.dat` files | 10 read by a deposited fragment, 9 not — table below |
| generator by-products | `figure_numbers.json`, `fig_akw_sampled_honest_caption.tex`, `fig_akw_sampled_honest_caption_macro.tex`, `fig_amort_native.tex` | **none of these four is read by anything.** The first three are written by `make_akw_sampled_honest_fig.py`; the two caption files carry **superseded** numbers and must not be substituted for `fig5_caption.tex` (`KNOWN_DISCREPANCIES.md` §26). `fig_amort_native.tex` is the body of a figure v3 withdrew, kept because `src/` still names it |
| `src/` | `fig_circuit_native_FIXED.tex`, `fig_akw_fixed.tex`, `fig_sqw_fixed.tex`, `fig_spinqw_fixed.tex`, `fig_hero_STANDALONE.tex`, `fig_noise_score_STANDALONE.tex`, the rasters `akw_field_v2.png`, `sqw_field.png`, `spinqw_field.png`, and a `README.md` | **deposited 2026-09-26**: the standalone sources of the six PDFs, recovered from the author's scratch directory. Each recompiles to a pixel-identical rendering of the committed PDF (216 dpi, same text). Nothing in the manuscript reads this folder |

### `paper/main.pdf`, the Physical Review A split, `paper/arxiv-submission.tar.gz`

| File | Contents |
|---|---|
| `main.pdf` | The compiled v3 manuscript — **68 pp**: main text pp. 1–28 (p. 28 holds the end of Sec. X, the acknowledgments and Appendix A, the proof of Theorem 1), the reference list pp. 29–33, then the Supplemental Material from p. 34; rebuilt on 2026-09-30 after the last source edit with the two `pdflatex` passes of `make paper` (0 errors, 0 undefined references or citations, no rerun request, 0 overfull `\hbox`, 0 overfull `\vbox` and no float-too-large warning after the layout pass of that day, which changed no text (before it, two overfull `\vbox`, 24.3 pt on p. 48 and 5.0 pt on p. 55, and a float-too-large warning of 30.2 pt for Fig. 3); `KNOWN_DISCREPANCIES.md` §34; the build of 2026-09-28, late, had one overfull `\vbox`, 12.8 pt on p. 48). Before the review pass of that day it was 66 pp (main text pp. 1–26, references pp. 27–31, SM from p. 32). It is **not** the arXiv v3 PDF: it carries the edits made after submission (README.md, "`paper/main.pdf` is the v3 manuscript"). Rebuild it at the end of any editing pass, then run `make pra`. |
| `main_pra.tex` | The **Physical Review A** main-text driver (added 2026-09-28). It holds no text of the paper: it runs `main.tex` itself and ends the document where `main.tex` opens the Supplemental Material (`\input{sm_begin.tex}`), i.e. after the reference list. It adds only `xr-hyper` (SM labels read from `sm_pra.aux`, `[nocite]`), the drop of REVTeX's duplicated bookkeeping labels, and four values written to `main_pra.aux` for `sm_pra.tex`. |
| `sm_pra.tex` | The **Physical Review A** Supplemental Material driver (added 2026-09-28). It runs `main.tex`, skips everything before `\input{sm_begin.tex}`, runs the SM `\input` list of `main.tex` (so the SM opens on p. 1 with the title block of `sm_begin.tex`), restores the theorem/proposition/footnote counters and REVTeX's table-strut state from `main_pra.aux`, and ends with a reference list of only the works the SM cites, in the order it first cites them, each entry taken as tokens from `bibliography.tex` at build time. |
| `pra_split.py` | Builds and checks the two drivers: `python paper/pra_split.py --build` (or `make pra`) builds `main.tex`, then `main_pra` and `sm_pra` alternately until neither asks for a rerun, in a temporary copy of `paper/`, and copies only `main_pra.pdf` and `sm_pra.pdf` back. It exits non-zero unless both logs have 0 errors, 0 undefined or multiply defined labels and citations and no overfull box or warning line the `main.tex` build lacks; the SM reference list equals the SM's citations in first-citation order; the APS reference to the SM names exactly the works cited only in the SM; and every page of `paper/main.pdf` is in exactly one of the two PDFs (main text pixel-identical below the page number; SM identical up to citation numbers, line for line). `--compare` runs only the page comparison. Needs `pdflatex`, `pdftotext`, `pdftoppm`, Pillow and an `xr-hyper` of 2023 or later (for `[nocite]`; built with v7.01o). Build `main.pdf` first: it is what the pages are compared against. |
| `main_pra.pdf`, `sm_pra.pdf` | The two builds of 2026-09-30: `main_pra.pdf` 33 pp. (= `main.pdf` pp. 1–33, pixel-identical below the page number), `sm_pra.pdf` 38 pp. (35 pp. of SM = `main.pdf` pp. 34–68, 20 pixel-identical and 15 differing in citation numbers only (19 and 16 before the layout pass of 2026-09-30), then 3 pp. of its own reference list, 92 entries; 18 / 17 and 90 entries in the builds of 2026-09-28). Regenerate with `make pra`. |
| `arxiv-submission.tar.gz` | The **frozen arXiv v2** bundle: 49 entries, 47 files. It is kept deliberately as the byte-exact record of what was posted, and must never be re-uploaded — it still carries the seven v2 body files, the withdrawn figures and the two fabricated rows of `figs/sqw_edges.dat`. Live comparison: `python src/check_tarball.py` (measured 2026-09-26: 9 byte-identical, 17 differ, 21 bundle-only). v3 was built separately, by `src/build_arxiv_bundle.py` (`ARXIV_SUBMISSION.md`). |

---

## `src/` — 66 scripts at the root, run from the repository root as `python src/<name>.py`; 29 more in `src/frontier/` and 9 files in `src/recovered/`

### The checkers

| File | Purpose |
|---|---|
| `verify.py` | **The adversarial guardian** (~15–35 s, numpy/scipy only). Recomputes the half-filled Hubbard ground states for both boundary conditions (open chain 9.6047, periodic ring 9.3851) and the geminal witness from first principles, then checks **934 checks / 11 674 numeric assertions** across `data/*.json`, `data/c3_frontier/*.json`, `paper/*.tex` (the abstract included), `paper/figs/*.tex` and `paper/figs/*.dat` against them. Its section **(10b)** opens the 54 C3 point files and re-derives all 204 grid rows — see [`C3_FRONTIER.md`](C3_FRONTIER.md) §7; seven negative controls for that section fire and name the defect. Exits non-zero and names the defect. `[XFAIL]` lines are registered open defects (`KNOWN_OPEN`), open claims (`CLAIMS_OPEN`) and declared coverage gaps (`COVERAGE_GAP_V3`); they do not set the exit code. `--fast` is for editing loops and **certifies nothing** (it exits 3). |
| `check_figures.py` | Regenerates eight generators into a throw-away tree — the six of 2026-09-19 and the two recovered into `src/recovered/` on 2026-09-26 — and diffs 19 artefacts against the committed ones; one provenance-timestamp comment line of the Fig. 6 fragment is normalised, and the run says so. It writes nothing into the repository (`KNOWN_DISCREPANCIES.md` §25, §29). |
| `check_provenance.py` | **New, 2026-09-19.** Derives the figure inventory from `paper/main.tex` and fails if `FIGURE_PROVENANCE.md`, the `main.tex` float census or `check_figures.py`'s artefact list no longer describes the manuscript. Prints the orphan inventory (it scans `src/`, `src/frontier/` and, since 2026-09-26, `src/recovered/`). Nine self-tests with negative controls run first; it exits 2 if one does not fire. |
| `check_coherence.py` | **New, 2026-09-19.** Checks the manuscript against itself: ten cross-section claims (the vacuity range, the published sizes, the largest sector, the hardware shot budget, the guardian counts quoted in the data statement, …) must be stated the same way wherever they are stated; two synthetic controls per check. |
| `check_trim.py` | **New, 2026-09-19.** The guardian of a length-reduction pass: against a git baseline (default `bc9670c`), no distinctive number and no concession, limitation, attribution or retraction may leave `paper/*.tex`. Items removed on purpose are listed in it with their proof. The proof for three of them — the literal `2608.16436` and two sentences on the retracted Theorem 1(iii) of v1–v2, moved out of the body into the arXiv Comments field — is read from `build/arxiv_comments.txt`, which is gitignored: on the author's tree the check passes, **on a fresh clone it reports FAIL** on those three (measured 2026-09-26). |
| `check_tarball.py` | Compares `paper/arxiv-submission.tar.gz` with the current `paper/` tree, so no document has to carry a hand-written count. |
| `build_arxiv_bundle.py` | **New, 2026-09-19.** Builds a submission bundle from the working tree (the v2 tarball is a record, not a source); it built the package submitted as v3 on 2026-09-25. |
| `abs_paste.py` | Writes the abstract exactly as pasted into the arXiv form (`build/abstract_arxiv.txt`) and counts it against arXiv's 1920-character limit. |
| `build_repro_notebook.py` | Assembles `notebooks/00_Reproduce_Everything.ipynb` from source. It is the v1–v2 pipeline (see `notebooks/` below). Measured 2026-09-26: the committed notebook differs from what the script writes in the wording of one markdown cell (the fifth), and in nothing else. |

### Shared exact-diagonalization engines

| File | Purpose |
|---|---|
| `akw_lanczos.py` | Sector-basis Hubbard engine: sparse ground state + Haydock continued fraction for `A(k,ω)` (PBC). Imported by others. |
| `scaling_lanczos.py` | Haydock/Lanczos scaling engine with explicit sparse H (good to L≈12). |
| `scaling_lanczos_mf.py` | Matrix-free, RAM-staged, **resumable** scaling to L=14 (28 qubits); checkpoints to `ckpt/`. |
| `_maxrigor.py` | Cross-checks: first-moment sum rule, peak dispersion, Haydock convergence. |
| `nonfreeness.py`, `rigor_floor.py` | Convention-free nonfreeness; shared-floor / CP-overhead structure. `cost_vs_entanglement.py` imports `rigor_floor`, so the module is load-bearing even though its JSON is not deposited (§13). |

### Spectral / dynamical response

| File | Purpose | → |
|---|---|---|
| `sqw_lanczos.py` | Exact `S(q,ω)` (density response, N sector). **Needs `--eta 0.18`**; it refuses to run at its 0.20 default (§1, §16). | `sqw_L12.json` |
| `spin_lanczos.py` | Exact `S^zz(q,ω)`, the spin channel (at L=12 a set of discrete poles; v3 calls it neither gapless nor a continuum, Sec. VIII). | `spinqw_L12.json`, `figs/spinqw_edges.dat` |
| `spin_peak_deviation.py` | Added 2026-09-28 (late): the deviation of the plotted spin peak (the grid argmax of each broadened line of `spinqw_L12.json`) from the des Cloizeaux–Pearson boundary, the percentages the Fig. 8 caption prints (about 3 % for q ≤ π/2, 10 % at 2π/3, 29 % at 5π/6); `verify.py` recomputes every row. | `spin_peak_deviation.json` |
| `sampled_akw.py` | `A(k,ω)` on the Born-ranked subspace — the **infinite-shot limit** of the protocol, not a finite-shot run. | `sampled_akw_L8.json` |
| `sampled_honest.py`, `honest_sampling_scaling.py` | The measured **finite-shot** sweeps that Table S6 and Fig. 3(c) report. | `sampled_honest.json`, `honest_sampling.json` |
| `fig3_finite_shot.py` | The finite-shot run of the literal Fig. 3(a,b) configuration (T = 2.6e6 and 5.2e6 shots per channel, 4 seeds, all 16 channels) and its post-processing: every finite-shot number of the Fig. 3 caption and Sec. S1. Written 2026-09-19, recovered from the session transcript and deposited 2026-09-28 (`KNOWN_DISCREPANCIES.md` §31); `--summary` prints the caption numbers. | `fig3_finite_shot_T2.6e6.json`, `fig3_finite_shot_T5.2e6.json` |
| `spectral_validation_max.py` | Max-level validation of the sampled `A(ω)` (Hubbard L=6, U/t=8). | `method_max.json` |
| `charge_gap_ed.py` | μ⁺, μ⁻ and the Mott gap Δ of the half-filled ring, L=6…12, by sector diagonalization. | `charge_gap.json` |
| `gap_scaling.py` | Finite-size scaling of the **two** gaps (spin and charge) with a Heisenberg control. | `gap_scaling.json` |
| `pole_structure.py` | Pole census behind Table S7. | `pole_structure.json` |
| `eta_sweep.py` | Resolution sweep behind the η prices of Sec. VII. | `eta_sweep.json` |
| `make_sqw_edges.py` | Generator of `figs/sqw_edges.dat` (band edges + peak). Until 2026-09-18 that file had **no generator anywhere**. | `figs/sqw_edges.dat`, `sqw_edges_provenance.json` |
| `sqw_field.py`, `recolor_akw_v2.py` | Rasters for the `S(q,ω)` / `A(k,ω)` colour maps. `sqw_field.py` writes the `S(q,ω)` field in viridis; the hot-coloured rasters the figures embed are deposited as rendered in `paper/figs/src/` (2026-09-26). `recolor_akw_v2.py` **cannot run here** — neither of its inputs is deposited (§13). | `build/` |

### The certificate (Theorem B and its falsification)

| File | Purpose | → |
|---|---|---|
| `leakage_certificate.py` | The certificate itself, with a self-test. | `leakage_certificate_selftest.json` |
| `leakage_certificate_suite.py` | The 606-subspace sweep behind Sec. V and Fig. S5. Writes the frontier files under a **computed** name, which is why `cert_frontier_*.json` appear orphaned. The Fig. S5 pool was built from an **earlier run** of its three scans, deposited in `data/thm1iii_violation/cert_earlier_run/` (2026-09-26). | `cert_akw.json`, `cert_stress.json`, `cert_teqsci.json`, `cert_frontier_6-8.json`, `cert_frontier_10.json` |
| `support_witness.py` | Support witness behind Sec. VI. | `support_witness.json` |
| `gflow_dequant.py` | Dequantization control of Sec. S16 (needs pyscf, torch, qiskit-ibm-runtime). | `gflow.json`, `gflow_runs/` |
| `moments_table.py` | Table II (`tab:moments`), its negative control and the two-level series of Sec. S10, rebuilt from scratch (2026-09-25). | `moments_table.json` |
| `apsg_witness.py` | The **provable** geminal witness: `F_k` large yet χ=2, `\|S\|=2^K`. | `apsg_witness.json` |

### Molecular suite and the resource theory

| File | Purpose | → |
|---|---|---|
| `n19_suite.py` | Nineteen molecules: FCI-verified `F₁`, `\|S\|`, χ at two geometries. **Needs `pyscf`.** | `n19_suite.json` |
| `n19_spectral.py` | Suite-wide sampled `A(ω)` reconstruction. **Needs `pyscf`.** | `n19_spectral.json` |
| `n2_hero_data.py` | N₂ dissociation: spectrum and resources together; the `F₁=2Nᵤ` identity. **Needs `pyscf`.** | `n2_hero.json`, `figs/hero_*.dat` |
| `cost_vs_entanglement.py` | The 30-point Hubbard sweep: the cost is governed by χ, not `F₁`. | `cost_vs_ent.json` |
| `free_fermion_support.py` | Added 2026-09-28 (late): the four free-fermion (U = 0) site-basis supports of Secs. VI and S13, `\|S\|_{ε_w}` = 35, 336, 3496, 37361 for L = 4, 6, 8, 10 — open chain, half filling, the smallest determinant set carrying all but ε_w = 2.5e-3 of the weight, inclusive at the exact four-way tie of L = 4 (a strict cut gives 36). Until then hand-typed constants with no producer; `verify.py` recomputes the counts independently and checks both sentences. | `free_fermion_support.json` |
| `ladder_vs_chain.py` | Chain vs two-leg ladder at matched `F₁`. | `ladder_vs_chain.json` |
| `gate1_ladder.py` | Response-targeted selection on the ladder. **Its output is not deposited** (§13). | — |
| `stats_resource.py` | Every published correlation, stratified; the exact permutation test. | `stats_resource.json` |
| `scaling_data.py` | Scaling data L=4,6,8 — `F₁` grows, subspace fraction falls. | `scaling_data.json` |
| `headtohead_ms.py` | Multi-seed selector benchmark (time-evolution vs Krylov vs CIPSI). | `headtohead_ms.json` |

### Hardware and noise

| File | Purpose | → |
|---|---|---|
| `molecular_noise_energy.py` | Noise-assisted SQD energy across the suite (simulated bit-flip + recovery). | `molecular_noise.json` |
| `molecular_noise_sweep.py` | The same, versus the per-qubit bit-flip rate ε. | `molecular_noise_sweep.json` |
| `noise_spectral.py` | Spectral rel-L1 vs ε (L=6 Hubbard). | `noise_spectral.json`, `figs/noise.dat` |
| `amortized_recovery.py` | Amortized configuration recovery, leave-one-out. **Its figure was withdrawn in v3**; the data stays. | `amortized_recovery.json` |

### Figure and table generators

**Live** — they produce something the v3 manuscript typesets:
`make_decoupling_native.py` (Figs. S2 **and** S3), `make_method_fig_max.py` (Fig. 2),
`make_akw_sampled_honest_fig.py` (Fig. 3; its absolute paths were repaired on 2026-09-19, §25),
`make_hardware_hero.py` (Fig. 9), `make_noise_recovery_native.py` (Fig. S8),
`make_hero_fig.py` (Fig. S6, into the cwd), `make_noise_fig.py` (Fig. S7, into the cwd),
`make_heron_fig.py` (writes `figs/heron*.dat`), `make_table.py` (Table S12, into `rebuild/`).

**Second writer, not the guarded one:** `make_resource_master_fig.py` writes
`figs/fig_resource_master_native.tex`, the same file `make_decoupling_native.py` writes.

**Hazard:** `make_scaling_fig.py` emits a *different* figure and would destroy panel (c) of Fig. 5. It
is deliberately out of `make figures` and out of `check_figures.py` (§3).

**Dead** — generators of the five figures v3 withdrew, kept as the record of how they were made:
`make_bench_fig.py` (`fig:bench`, and it would recreate `figs/bench_U*.dat`),
`make_gflow_fig.py` (`fig:gflow`), `make_gallery_native.py` (`fig:gallery`),
`make_magic_suite_fig.py` (`fig:molsuite`).

**Recovered on 2026-09-26** (absent until then, and named by their own output): `make_fig_gapscaling.py`
(Fig. 6) and `build_n3.py` with `fix_n3_caption.py` / `fix_n3_fig5_series.py` (Fig. S5), now in
`src/recovered/` below. See `KNOWN_DISCREPANCIES.md` §27 and §29.

---

### `src/frontier/` — 29 scripts, the C3 certificate-frontier sweep (deposited 2026-09-19; four more on 2026-09-25; `rerank_L14.py` on 2026-09-28)

Sec. V used to be the one part of the manuscript a reader could not regenerate. This is what closes
that. Full reproduction instructions, per-file costs and the declared gaps are in
[`C3_FRONTIER.md`](C3_FRONTIER.md); the module names are kept exactly as they were run, because
renaming them would have meant rewiring the import graph rather than moving path constants.

| group | files | role |
|---|---|---|
| **A — the engine** | `c0_lib.py`, `c0_big.py`, `c0_sweep.py`, `fcore.py`, `frun.py`, `frun12.py`, `fsum.py` | the Hubbard sector (matrix-free above L=8), the Born ranking at the repository's own protocol (K=18, dt=0.5, sub=16), the memory-lean Gram identity that makes the sweep affordable, the two point drivers (`frun12` stores the Krylov basis in `|S|` coordinates), and the aggregation to `C3_grid.json` |
| **B — the gates** | `fvalid.py`, `fproto.py`, `ftie.py` | six algebraic controls at L=6 and the dense cross-validation of Λ_S at L=8; the proof that the ranking *is* the repository's (set overlap 1.000000 against `L10_order.npy`); the census of zero-Born determinants and `argsort` ties at the cut |
| **C — the proof ladder** | `fron_lib.py`, `run_ladder.py`, `verdict.py` | the exact five-rung instrumentation of the Theorem-B proof — no Krylov anywhere — and the verdict on how much of the frontier is Cauchy–Schwarz and how much is physics |
| **D — adversarial re-runs** | `z_fine.py`, `z_rank.py`, `z_e0.py`, `z_suppK.py` | the four re-runs named in the provenance note of `paper/sec_5_body.tex`: fine-grained frontier, ranking sensitivity, an *estimated* E₀, and the support of the Krylov space |
| **D′ — the support of the Krylov space, measured (2026-09-25)** | `z_kblock.py`, `z_momseed_support.py`, `z_suppK_power.py`, `z_suppK_sym.py` | the momentum-block route of Sec. V E (the Krylov space of a momentum probe lies in one momentum block, dim/L determinants); the coordinate support of the momentum seeds of Fig. 3 at L=6 and 8; the support of K(H, φ) by powers (a lower bound, for the sizes where the projector cannot be stored); and the corrected estimator, exact in the reflection symmetry. Outputs in `data/c3_frontier/adversarial/` |
| **E — published fractions** | `pubsum.py`, `run_L14.py`, `rerank_L14.py` | the four points at the fractions the manuscript actually publishes, and the L=14 point at FR=0.08; `rerank_L14.py` (2026-09-28, late) re-runs the K=18 ranking protocol at L=14 on the deposited checkpoint and compares the selection with the stored `ckpt/L14_order.npy` as sets (`C3_FRONTIER.md` §5) |
| **F — post-processing** | `calib_oos.py`, `null_ws4.py`, `null_ws5.py`, `null_ws6.py`, `null_ws7.py` | zero new compute: the out-of-sample test of the calibration, and the null a referee will propose — does `1−w_S` locate the error as well as the certificate? |

### `src/recovered/` — 9 files: 8 recovered 2026-09-26, and `n3_fig5_whiskers.py` added 2026-09-28

Two figure generators and the scripts behind two of their inputs, written on 2026-09-18 in the
author's session scratch directory and never committed until they were recovered on 2026-09-26. The
folder's `README.md` gives the detail; `KNOWN_DISCREPANCIES.md` §29 the measurements.

| File | Status | Purpose |
|---|---|---|
| `make_fig_gapscaling.py` | ported, guarded | Fig. 6 (`fig:gapscaling`) → `paper/figs/fig_gapscaling_native.tex`; caption and facts to `build/` |
| `build_n3.py` | ported, guarded | the ten `paper/figs/n3_*.dat` of Fig. S5 from `data/thm1iii_violation/` (the two thm3 JSONs) and, since 2026-09-28 (late), the deposited `data/cert_{stress,akw,teqsci}.json` — no longer the earlier run in `cert_earlier_run/` |
| `n3_fig5_whiskers.py` | added 2026-09-28 (late), guarded | run after `build_n3.py`: the min/median/max whiskers of the 64 Fig.-3 rows, written to `data/thm1iii_violation/n3_fig5_summary.json` and into the two whisker blocks of the Fig. S5 fragment (what `fix_n3_fig5_series.py` did once by hand) |
| `fix_n3_fig5_series.py`, `fix_n3_caption.py` | record — do not run | the one-shot patchers of the Fig. S5 fragment and caption |
| `thm3_run.py`, `thm3_core.py`, `thm3_L8b.py` | record (numpy 1.x) | the producers of `thm3_results.json` and `thm3_L8b.json` |
| `README.md` | — | what was recovered, what was changed, what was measured |

---

## `data/c3_frontier/` — 96 files, the C3 sweep itself (`published_fraction/L14_rerank.json`, the K=18 ranking re-run at L=14, added 2026-09-28, late; `C3_FRONTIER.md` §5)

54 point files, the 204-row grid they aggregate to, three gate outputs, three proof ladders, the
verdict, and three sub-directories (`published_fraction/`, `adversarial/` — including the Krylov-support
measurements of 2026-09-25 — and `null/`). Every number in
Sec. V comes from here. Documented file by file in [`C3_FRONTIER.md`](C3_FRONTIER.md) §1; `data/` is
excluded from the arXiv package, so the deposit costs the submission nothing.

---

## `data/` — 43 datasets at the root, and three sub-directories

| file | feeds |
|---|---|
| `akw_lanczos_L12.json` | Fig. 4 `fig:lattice` (exact `A(k,ω)`, L=12) |
| `sampled_akw_L8.json` | Fig. 3(a,b) `fig:akwsampled` |
| `sampled_honest.json`, `honest_sampling.json` | Fig. 3(c) and Table S6 (the finite-shot sweeps) |
| `fig3_finite_shot_T2.6e6.json`, `fig3_finite_shot_T5.2e6.json` | the finite-shot paragraph of the Fig. 3 caption and Sec. S1 (the configuration of Fig. 3(a,b) at finite shots; each file carries its reproduction check on the former window) |
| `method_max.json` | Fig. 2 `fig:method` + `figs/aw_method.dat` |
| `sqw_L12.json` | Fig. 7 `fig:sqw` + `figs/sqw_edges.dat` |
| `spinqw_L12.json` | Fig. 8 `fig:spin` + `figs/spinqw_edges.dat` |
| `spin_peak_deviation.json` | the peak-dispersion percentages of the Fig. 8 caption (2026-09-28, late) |
| `free_fermion_support.json` | the four free-fermion supports printed in Secs. VI and S13 (2026-09-28, late) |
| `gap_scaling.json` | Fig. 6 `fig:gapscaling` (and Table V) |
| `charge_gap.json` | the Δ = 4.97 t of Sec. VIII and Table V |
| `pole_structure.json` | Table S7 `tab:poles` |
| `cost_vs_ent.json` | Fig. S2(b,c) `fig:decoupling` |
| `n19_suite.json` | Fig. S2(a), Table S12, Table S10 |
| `n19_spectral.json` | Table S12 (the reconstruction columns) |
| `resource_master.json` | Fig. S3 `fig:master` |
| `stats_resource.json` | Fig. S2(d), Sec. S14, Table S9 |
| `apsg_witness.json` | Fig. S1 `fig:witness` |
| `support_witness.json` | Secs. III and VI (the support witness; the 0.43 / 0.35 of the abstract) |
| `ladder_vs_chain.json` | Fig. S4 `fig:ladder` |
| `scaling_data.json` | Fig. 5 `fig:scaling` |
| `eta_sweep.json` | Sec. VII (the resolution price) |
| `cert_akw.json`, `cert_stress.json`, `cert_teqsci.json`, `cert_frontier_6-8.json`, `cert_frontier_10.json`, `leakage_certificate_selftest.json` | Secs. V and S9 (the re-audit of the 378 rows). **Not** the scans behind Fig. S5: those are an earlier run, in `thm1iii_violation/cert_earlier_run/` (on these files the pool gives 469 and 2.005 instead of 468 and 2.10, as Sec. S9 says) |
| `moments_table.json` | Table II `tab:moments`, written by `src/moments_table.py` |
| `hw_lucj_n2_result.json` | ⭐ Fig. 9(b) — **job `da125f2ein7c73bcsqs0`** (authoritative N₂ run) |
| `heron_spectral.json` | Fig. 9(a) — the `ibm_fez` L=6 `A(ω)` run → `figs/heron_hot.dat` |
| `noise_spectral.json` | Fig. S7 `fig:noise` |
| `molecular_noise.json`, `molecular_noise_sweep.json` | Fig. S8 `fig:noiserec` |
| `n2_hero.json` | Fig. S6 `fig:hero` |
| `headtohead_ms.json` | **no figure in v3** — `fig:bench` was withdrawn; the data stays |
| `gflow.json`, `gflow_runs/` (20 per-seed files) | Sec. S16 (the dequantization control), written by `src/gflow_dequant.py` |
| `amortized_recovery.json` | **not used by the paper** — no figure, no table and no number in the text |
| `sqw_edges_provenance.json` | the provenance record of `figs/sqw_edges.dat` |
| `sampled_akw_figure_numbers.json` | the number→source manifest of Fig. 3 |
| `thm1iii_violation/` (7 files, 2026-09-26) | the inputs of Fig. S5: `thm3_results.json`, `thm3_L8b.json`, `n3_fig5_summary.json` (rewritten by `n3_fig5_whiskers.py` since 2026-09-28, late), the three earlier-run certificate scans in `cert_earlier_run/` (a record since that day: the pool is built from the deposited `data/cert_*.json`), and a `README.md` |
| `c3_frontier/` (95 files) | Sec. V and Tables III, S1–S3 — see above and [`C3_FRONTIER.md`](C3_FRONTIER.md) |

Backups (`*_backup.json`) are not the plotted data and are `.gitignore`d.

---

## Orphan inventory

Produced by `python src/check_provenance.py -v`, re-measured 2026-09-26 over the 79 artefacts that are
*supposed* to have a producer (`paper/figs/*.dat|pdf|json`, machine-written `paper/figs/*.tex`,
`data/*.json`, `docs/img/*.png`), with `src/`, `src/frontier/` and `src/recovered/` counted as code.
Hand-written LaTeX is excluded — it has no producer by design.

The tiers use two facts that need no guessing: does any deposited **code** name the file, and does the
**manuscript** reach it. Writing a file down in a document does not give it a producer, so documenting
an orphan does not remove it from tier 1.

### Tier 1 — no deposited code names it, and the manuscript does not use it (10)

| file | why it is here |
|---|---|
| `data/cert_frontier_6-8.json`, `cert_frontier_10.json` | written by `leakage_certificate_suite.py:697` under a **computed** filename (`"cert_frontier_%s.json" % "-".join(...)`), so the literal names appear in no script. Only a comment in `paper/sec_5_body.tex:70` mentions them, and nothing reads them back. |
| the eight `docs/img/*.png` | the README thumbnails, rendered on 2026-09-26 from the v3 figures by `docs/img/make_thumbnails.py` (figure body only, no caption). That script sits next to them, not in `src/`, so the checker does not count it as code. The ten thumbnails they replaced were hand-cropped from the v1 PDF and showed withdrawn or corrected content; they are in `_superseded/v1_readme_img/`. |

*Left tier 1 on 2026-09-26: `paper/figs/n3_fig5_eta.dat` and `n3_fig5_frac.dat`, the 64 Fig.-3 rows
summarised into the whiskers of Fig. S5 — `src/recovered/build_n3.py` now writes them.*

### Tier 2 — code names it, the manuscript does not use it (51)

Mostly datasets whose figure was withdrawn or whose only remaining reader is the guardian:
`gflow.json`, `amortized_recovery.json`, `headtohead_ms.json`, the six certificate files,
`eta_sweep.json`, `support_witness.json`, `sampled_akw_figure_numbers.json`,
`figs/figure_numbers.json`, the seven `.dat` tables read only by standalone sources
(`hero_aw.dat`, `hero_res.dat`, `heron.dat`, `heron_hw.dat`, `noise.dat`, `spinqw_edges.dat`,
`sqw_edges.dat`), since 2026-09-26 the two `n3_fig5_*.dat` that `src/recovered/build_n3.py`
writes, and since 2026-09-28 the two `fig3_finite_shot_T*.json` behind the finite-shot paragraph of the Fig. 3
caption (a caption quotes them; no float reads them; `verify.py` checks every printed value against them), and since late that day `free_fermion_support.json` and `spin_peak_deviation.json`, read by `verify.py` and quoted by two sentences and a caption.
*(Counted by `python src/check_provenance.py -v` on 2026-09-28, late: 51; 49 earlier that day.)* **None of these is a defect** — a dataset with a producer and a guardian is doing
its job — but none of them is load-bearing for the v3 manuscript either.

### Tier 3 — the manuscript typesets it and no deposited code names it (6)

The six PDFs (`fig_akw_native.pdf`, `fig_circuit.pdf`, `fig_hero.pdf`, `fig_noise_score.pdf`,
`fig_spinqw.pdf`, `fig_sqw.pdf`). Their *numbers* are in committed JSONs, and since 2026-09-26 their
standalone LaTeX sources are in `paper/figs/src/`, where each recompiles to a pixel-identical rendering;
LaTeX is not code for this checker, and no script regenerates them, so they stay in this tier. *(Until
2026-09-26 this tier also held seven of the eight `n3_*.dat` tables Fig. S5 plots.)*

### `paper/figs/*.dat` — who reads them

| file | read by a fragment **in this deposit**? | read by | written by |
|---|---|---|---|
| `aw_method.dat` | **yes** | `fig_method_native_frag.tex` | `make_method_fig_max.py` |
| `heron_hot.dat` | **yes** | `fig_hardware_hero_frag.tex` | `make_hardware_hero.py` |
| `n3_pub_eta_hi/lo.dat`, `n3_pub_frac_hi/lo.dat`, `n3_thmB_eta/frac.dat`, `n3_inf_eta/frac.dat` | **yes** (8 files, 2432 lines) | `fig_thm1iii_violation_native.tex` | `src/recovered/build_n3.py` (recovered 2026-09-26; nothing deposited before that, §27) |
| `n3_fig5_eta.dat`, `n3_fig5_frac.dat` | no | nothing; `fix_n3_fig5_series.py` summarised them into literals of the Fig. S5 fragment | `src/recovered/build_n3.py` |
| `hero_aw.dat`, `hero_res.dat` | no | the standalone `paper/figs/src/fig_hero_STANDALONE.tex` (deposited 2026-09-26) — note it is *not* `paper/carried/fig_hero.tex`, which is the float wrapper | `make_hero_fig.py` |
| `noise.dat` | no | the standalone `paper/figs/src/fig_noise_score_STANDALONE.tex` (deposited 2026-09-26) | `make_noise_fig.py` |
| `heron.dat`, `heron_hw.dat` | no | `work/notes/fig_heron_native.tex`, a **superseded** version of the hardware figure — `work/` is the author's sibling working tree and is **NOT in this deposit** | `make_heron_fig.py` |
| `spinqw_edges.dat` | no | an older standalone `work/notes/fig_spinqw.tex`, **not in this deposit**; the deposited `paper/figs/src/fig_spinqw_fixed.tex` does not read it | `spin_lanczos.py` |
| `sqw_edges.dat` | no | the standalone `paper/figs/src/fig_sqw_fixed.tex` (deposited 2026-09-26) | `make_sqw_edges.py` (since 2026-09-18; before that, **nothing**) |

---

## `_superseded/` — 36 tracked files, kept so the history is legible

`_superseded/README.md` has one line per subfolder.

| subtree | contents |
|---|---|
| `v2_paper/` | the v2-era manuscript as it stood in `paper/` on 2026-09-19 (edited after the 2026-08-22 posting; the posted v2 source is `paper/arxiv-submission.tar.gz` — measured 2026-09-26, six of these ten source files differ from it, and `main.pdf` is a 40-page build, not the posted PDF): `main.tex`, `main.pdf`, `introduction.tex`, `method.tex`, `resource.tex`, `results.tex`, `hardware.tex`, `discussion.tex`, `conclusion.tex`, `theorem_b2.tex`, `table_molecules.tex`. Every one of them was deleted from `paper/` in `606f969` and is preserved here, so the deletion commit does not lose the trail. (`main.bbl` is on disk but untracked: `.gitignore`'s `*.bbl` excludes it, and the v2 bibliography was hand-written anyway.) |
| `v2_figs/` | `fig_akw_sampled_native.tex` and `make_akw_sampled_fig.py` — the Fig.-5 fragment retired on 2026-09-19 because it labelled an infinite-shot curve "sampled (85% sector)", and its generator. Replaced by `fig_akw_sampled_honest.tex`. |
| `v1_readme_img/` | the ten README thumbnails cropped from the v1 PDF on 2026-08-17 and superseded on 2026-09-26 (`scaling.png` was withdrawn v1 text, `sqw.png` showed the two fabricated `sqw_edges.dat` rows, `spinqw.png` the word GAPLESS, `hardware.png` the 50,000-shot caption); its `README.md` says what each showed |
| `v3_retired/figs/` | the bodies of the five figures v3 withdrew and their data: `fig_benchmark.pdf` + `bench_U4.0/8.0/12.0.dat` (`fig:bench`), `fig_gflownet.pdf` (`fig:gflow`), `fig_molecular_gallery.pdf` (`fig:gallery`), `fig_molecular_suite.pdf` (`fig:molsuite`), plus `fig_hardware_hero.pdf`, `heron_spectral.pdf` and `_method_standalone.tex/.pdf`, earlier standalone builds that the native fragments replaced. |

**Nothing in `paper/` reads any of it.** The one live consequence: `make_bench_fig.py` is still in
`src/` and would recreate `paper/figs/bench_U*.dat` if run.

The body of the fifth withdrawn figure, `fig_amort_native.tex`, was **not** moved here — it is still
in `paper/figs/`, because `src/verify.py` names it. That asymmetry is deliberate and is recorded in
`verify.py`'s `COVERED_FRAGMENTS` comment.

---

## `notebooks/`

| File | Contents |
|---|---|
| `00_Reproduce_Everything.ipynb` | **The v1–v2 narrated pipeline, kept as a record** — written for arXiv:2608.16436v1–v2, it reproduces the v1–v2 figure set, several of which v3 withdrew, and it is **not** the reproduction path of the v3 manuscript (use `src/verify.py` and `REPRODUCE.md`); its first cell says so since 2026-09-26. Assembled by `src/build_repro_notebook.py`. **Needs `pyscf`; not verified in this pass.** |
| `HW_LUCJ_N2_Heron_READY.ipynb` | N₂ energy on `ibm_marrakesh`: LUCJ circuit, sampling, self-consistent recovery, 8-seed error-bar cell. Credentials scrubbed. **Its saved cell outputs are an older run** (`d9qe731dsedc73af67d0`) and are not the paper's numbers. A first cell added on 2026-09-26 says that the manuscript reports this run as *executed on* the device, not as a demonstration, and that the record cannot say whether noise helps. |
| `Spectral_Heron.ipynb` | `A(ω)` on `ibm_fez`: the L=6 Hubbard spectral run (full-sector coverage). Kept as it ran: its post-processing reads the bitstrings in reversed order and applies post-selection only, with no configuration recovery, and the `# TREX` of its sampler cell is measurement twirling, not readout-error mitigation (`KNOWN_DISCREPANCIES.md` §30; a note in its first cell says so since 2026-09-26). Credentials scrubbed. |

Both hardware notebooks have the API **token and instance CRN removed** (placeholders), and neither
ships raw counts: they fetch them from the IBM Quantum service by job id.

---

## Support files

| File | Purpose |
|---|---|
| `README.md` | Front page, figure gallery, quick start. |
| `CITATION.cff` | Machine-readable citation (arXiv DOI `10.48550/arXiv.2608.16436` and the Zenodo concept DOI). |
| `.zenodo.json` | Zenodo deposit metadata (title, creator, licence, the `isSupplementTo` link to arXiv:2608.16436). |
| `LICENSE` | MIT (code) + CC-BY-4.0 (paper, figures and data). |
| `Makefile` | `make verify` / `check-figures` / `check-coherence` / `figures` / `data` / `paper` / `reproduce` / `ci`. Measured results per target in `KNOWN_DISCREPANCIES.md` §10. |
| `requirements.txt` | Python dependencies. |
| `.gitignore` | Includes `paper/*Notes.bib` (a TeXstudio artefact that would otherwise be committed) and `ckpt/` (except its three deposited files), `rebuild/`, `build/`, `*.log`. |
| `.github/workflows/ci.yml` | Three jobs: `python src/verify.py`, `make check-figures` and, since 2026-09-26, `make check-coherence` — the three commands of the `Makefile`'s `ci` target. Tracked since 2026-09-18; it has run, green, on every push since 2026-09-19 (`KNOWN_DISCREPANCIES.md` §24). |
| `docs/REPRODUCE.md` | Figure/number → script → command, and what a clean clone cannot reach. |
| `docs/FIGURE_PROVENANCE.md` | Figure → source → generator → data → job id, plus the stale-file list. Checked by `src/check_provenance.py`. |
| `docs/KNOWN_DISCREPANCIES.md` | Every place where running a committed script does **not** reproduce the deposit, repaired or not. **Opens with a banner about `paper/arxiv-submission.tar.gz`, which must never be uploaded as v3.** |
| `docs/DATA_AVAILABILITY.md` | The Data Availability Statement drafted for Physical Review A, plus the audit licensing each of its sentences. |
| `docs/ARXIV_SUBMISSION.md` | How the v3 source package was built and submitted (2026-09-25), the form fields, and the Comments filed. |
| `docs/C3_FRONTIER.md` | The C3 certificate-frontier sweep of Sec. V, file by file. |
| `docs/img/*.png`, `docs/img/make_thumbnails.py` | Eight README thumbnails, rendered on 2026-09-26 from the v3 figures (figure body only) by the script beside them, which needs PyMuPDF (tier-1 orphans above). |
| `ckpt/L14_gs.npz`, `ckpt/L14_order.npy`, `ckpt/README.txt` | The L=14 ground state (94 MB) and its Born ranking (41 MB), deposited 2026-09-25 with their SHA-256 in the README, because they take hours to regenerate; `src/frontier/run_L14.py` reuses them. |
| `src/recovered/README.md`, `data/thm1iii_violation/README.md`, `paper/figs/src/README.md`, `_superseded/README.md`, `_superseded/v1_readme_img/README.md` | The five READMEs of 2026-09-26, one per folder added or relabelled that day. |

### Present on disk, not in the deposit (`.gitignore`d)

`ckpt/` apart from its three deposited files (resumable checkpoints of `scaling_lanczos_mf.py`),
`rebuild/` (scratch output of `make_table.py`), `build/` (raster and log by-products, the caption and
summary files of the two recovered generators, and the arXiv v3 package with its filed Comments and
abstract), `src/__pycache__/`, `paper/main.{aux,log,out,bbl,blg}`, and two stray root logs
`texput.log` / `x2.log` left by a TeX run.
