# FILE INDEX — every file in this repository, described

*Rebuilt 2026-09-19 against the v3 REVTeX manuscript. The version before this one described the **v2**
file set in its entirety — `introduction.tex`, `method.tex`, `results.tex`, `discussion.tex`,
`conclusion.tex`, `hardware.tex`, `theorem_b2.tex`, twenty figures, a 40-page PDF. None of those files
has been in `paper/` since commit `606f969`. See [`KNOWN_DISCREPANCIES.md`](KNOWN_DISCREPANCIES.md) §28.*

Grouped by role. For *how to run* each script, see [`REPRODUCE.md`](REPRODUCE.md); for *which data feeds
which figure*, see [`FIGURE_PROVENANCE.md`](FIGURE_PROVENANCE.md); for *where the deposit and the code
disagree*, see [`KNOWN_DISCREPANCIES.md`](KNOWN_DISCREPANCIES.md).

**Counted, not asserted.** 225 files: 223 tracked by git at commit `28c12ea`, plus `src/check_provenance.py`
and `src/build_arxiv_bundle.py`, added in this pass. `ckpt/`, `rebuild/`, `build/`, `__pycache__/`,
`*.log` and `paper/main.{aux,log,out}` are present on disk and `.gitignore`d; they are build scratch and
are described at the end rather than listed.

| directory | tracked files |
|---|---|
| `paper/` | 81 — 23 at the root, 9 in `carried/`, 49 in `figs/` |
| `src/` | 55 tracked + 2 new = **57**, plus **24 in `src/frontier/`** (the C3 sweep, deposited 2026-09-19) |
| `data/` | 37, plus **82 in `data/c3_frontier/`** (the C3 sweep, 2.47 MB) |
| `docs/` | 17 — 7 Markdown, 10 PNG thumbnails |
| `_superseded/` | 24 |
| `notebooks/` | 3 |
| root | 7 — `README.md`, `CITATION.cff`, `LICENSE`, `Makefile`, `requirements.txt`, `.gitignore`, `.github/workflows/ci.yml` |

---

## `paper/` — the v3 manuscript (REVTeX 4.2, two-column, Physical Review A)

`main.tex` is the only hand-written driver. Everything it reaches is listed here, and the list is
derived: `python src/check_provenance.py` expands `\input` from `main.tex` and reports **50 source
files reached** and **33 floats — 17 figures and 16 tables, all double-column**.

### Body

| File | Contents |
|---|---|
| `main.tex` | Preamble (REVTeX options, the `figladder` font restoration, float placement, theorem environments, hyperref), title, author block, abstract, and the `\input` list. The author block is marked `*** NEEDS NICOLAS'S DECISION ***`: `\affiliation{Daita AI}` carries no postal address and the printed e-mail is on a different domain. |
| `sec_1.tex` | **I. The operational question** (`sec:intro`). |
| `sec_2.tex` | **II. One sampling primitive, four channels** (`sec:method`) — holds Table I (`tab:provenance`) and Fig. 1 (`fig:circ`). |
| `sec_3_body.tex` | **III. Two budgets: what the shots buy and what they do not** (`sec:leakage`) — the central theorem. |
| `sec_3_1.tex` | **III.1 Retraction of Theorem 1(iii)** (`sec:retraction`), `\input` from inside Sec. III. |
| `sec_4_body.tex` | **IV. What the structure does certify: moment exactness** (`sec:moments`) — Table II. |
| `sec_5_body.tex` | **V. What the certificate can and cannot assert** (`sec:frontier`) — Tables III–V. |
| `sec_6_body.tex` | **VI. No shortcut to the boundary: an obstruction and its witness** (`sec:nogo`) — Fig. 5, Table VI. |
| `sec_7.tex` | **VII. Three prices: support, resolution, and shots** (`sec:prices`) — Tables VII–IX. |
| `sec_8.tex` | **VIII. Physical validation: the two gaps scale in opposite directions** (`sec:physics`) — Tables X–XI. |
| `sec_9.tex` | **IX. Hardware, scope, and the region where the method has no advantage** (`sec:scope`, `sec:honest`) — Table XII. Its §9.4 withdraws two statistics *and the two figures that drew them*; its header records that `fig:gflow` and `fig:amort` are deliberately absent. |
| `sec_10.tex` | **X. Conclusion** (`sec:conclusion`). |
| `bibliography.tex` | Hand-written `thebibliography` — 133 entries, arXiv-safe, no external `.bib`/`.bbl`. |
| `table_molecules_v3.tex` | Table XIII (`tab:molecules`), the nineteen-molecule suite. The v2 `table_molecules.tex` is in `_superseded/v2_paper/`. |

### Appendices

| File | Contents |
|---|---|
| `sec_3_app.tex` | **A.** Counterexamples to the retracted bound, and the pooled sweep. `\input`s `figs/float_n3.tex` (Fig. 17). |
| `sec_4_app.tex` | **B.** Moment exactness: statement, proof, sharpness, shot count. |
| `sec_5_app.tex` | **C.** The leakage certificate: hypotheses, protocol, anatomy of the threshold. |
| `sec_6_app.tex` | **D.** The geminal witness: proof and angle sweep — Table XIV. |
| `app_carried_theorem.tex` | **E.** One-body magic cannot bound the determinant support (carried from v2's `theorem_b2.tex`). |
| `app_stats.tex` | **F.** The resource statistics, stratified — Table XV. |
| `app_carried_repro.tex` | **G.** Reproducibility data for the molecular suite — Table XVI (`tab:repro`). **Corrected 2026-09-19:** `R_eq(H₂O)` now prints `0.958`, the value the geometry in `src/n19_suite.py:123` gives (`√(0.757² + 0.587²) = 0.957924 Å`); it printed `0.957`. |

### `paper/carried/` — float wrappers carried verbatim from the v2 tree

Nine one-float files, each `\begin{figure*}` + `\input` of a fragment (or `\includegraphics`) +
caption: `fig_method.tex`, `fig_lattice.tex`, `fig_ladder.tex`, `fig_sqw.tex`, `fig_spin.tex`,
`fig_hero.tex`, `fig_heron.tex`, `fig_noise.tex`, `fig_noiserec.tex`. They exist because v3 keeps the
v2 captions unchanged for those nine figures; the float and the `\input` share a line, which is why a
line-level `\input` expander finds two floats in this manuscript instead of thirty-three.

`carried/fig_amort.tex` and `carried/fig_gflow.tex` do **not** exist: their `\input` lines in
`main.tex` are commented out, with the reason (Sec. 9.4 withdraws both figures).

### `paper/figs/` — 49 files

| group | files | note |
|---|---|---|
| float wrappers | `float_f5.tex`, `float_f12.tex`, `float_n2.tex`, `float_n3.tex`, `captions.tex` | hand-written; `captions.tex` carries two floats (`fig:decoupling`, `fig:master`) |
| captions | `fig5_caption.tex`, `caption_fig12.tex`, `fig_gapscaling_caption.tex`, `caption_thm1iii_violation.tex` | the numbers in each are traced in `FIGURE_PROVENANCE.md` |
| native fragments (typeset) | `fig_method_native_frag.tex`, `fig_akw_sampled_honest.tex`, `fig_witness_native.tex`, `fig_decoupling_native.tex`, `fig_resource_master_native.tex`, `fig_ladder_native.tex`, `fig_scaling2_native.tex`, `fig_gapscaling_native.tex`, `fig_hardware_hero_frag.tex`, `fig_noise_recovery_native.tex`, `fig_thm1iii_violation_native.tex` | 11 fragments; `src/verify.py` opens all 11 and declares one coverage gap |
| compiled figure bodies | `fig_circuit.pdf`, `fig_akw_native.pdf`, `fig_sqw.pdf`, `fig_spinqw.pdf`, `fig_hero.pdf`, `fig_noise_score.pdf` | 6 PDFs; none has a `.tex` source in this deposit (`KNOWN_DISCREPANCIES.md` §5) |
| plotted data | 19 `.dat` files | 10 read by a deposited fragment, 9 not — table below |
| generator by-products | `figure_numbers.json`, `fig_akw_sampled_honest_caption.tex`, `fig_akw_sampled_honest_caption_macro.tex`, `fig_amort_native.tex` | **none of these four is read by anything.** The first three are written by `make_akw_sampled_honest_fig.py`; the two caption files carry **superseded** numbers and must not be substituted for `fig5_caption.tex` (`KNOWN_DISCREPANCIES.md` §26). `fig_amort_native.tex` is the body of a figure v3 withdrew, kept because `src/` still names it |

### `paper/main.pdf`, `paper/arxiv-submission.tar.gz`

| File | Contents |
|---|---|
| `main.pdf` | The compiled v3 preprint — **59 pp**, body 1–47, appendices 48–55, references 56–59; 0 errors, 0 undefined references or citations. **It must be rebuilt at the end of any editing pass**: `paper/app_carried_repro.tex` and `paper/main.tex` were edited on 2026-09-19 and one of those edits changes a printed number. |
| `arxiv-submission.tar.gz` | The **frozen arXiv v2** bundle: 49 entries, 47 files. It is kept deliberately as the byte-exact record of what was posted, and must never be re-uploaded — it still carries the seven v2 body files, the withdrawn figures and the two fabricated rows of `figs/sqw_edges.dat`. Live comparison: `python src/check_tarball.py` (measured 2026-09-19: 9 byte-identical, 17 differ, 21 bundle-only). |

---

## `src/` — 57 scripts, run from the repository root as `python src/<name>.py`

### The three checkers

| File | Purpose |
|---|---|
| `verify.py` | **The adversarial guardian** (~15 s, numpy/scipy only). Recomputes the half-filled Hubbard ground states for both boundary conditions (open chain 9.6047, periodic ring 9.3851) and the geminal witness from first principles, then checks **831 checks / 9 846 numeric assertions** across `data/*.json`, `data/c3_frontier/*.json`, `paper/*.tex` (the abstract included), `paper/figs/*.tex` and `paper/figs/*.dat` against them. Its section **(10b)** opens the 54 C3 point files and re-derives all 204 grid rows — see [`C3_FRONTIER.md`](C3_FRONTIER.md) §7; seven negative controls for that section fire and name the defect. Exits non-zero and names the defect. `[XFAIL]` lines are registered open defects (`KNOWN_OPEN`), open claims (`CLAIMS_OPEN`) and declared coverage gaps (`COVERAGE_GAP_V3`); they do not set the exit code. `--fast` is for editing loops and **certifies nothing** (it exits 3). |
| `check_figures.py` | Regenerates six generators into a throw-away tree and diffs 8 artefacts against the committed ones. **Read its PASS as 7, not 8** — one artefact is compared with a copy of itself, and the run rewrites four files in the working tree; measured, with a negative control, in `KNOWN_DISCREPANCIES.md` §25. |
| `check_provenance.py` | **New, 2026-09-19.** Derives the figure inventory from `paper/main.tex` and fails if `FIGURE_PROVENANCE.md`, the `main.tex` float census or `check_figures.py`'s artefact list no longer describes the manuscript. Prints the orphan inventory. Nine self-tests with negative controls run first; it exits 2 if one does not fire. |
| `check_tarball.py` | Compares `paper/arxiv-submission.tar.gz` with the current `paper/` tree, so no document has to carry a hand-written count. |
| `build_arxiv_bundle.py` | **New, 2026-09-19.** Builds a submission bundle from the working tree (the v2 tarball is a record, not a source). |
| `build_repro_notebook.py` | Assembles `notebooks/00_Reproduce_Everything.ipynb` from source; the notebook is generated, not hand-edited. |

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
| `spin_lanczos.py` | Exact `S^zz(q,ω)` (gapless two-spinon continuum). | `spinqw_L12.json`, `figs/spinqw_edges.dat` |
| `sampled_akw.py` | `A(k,ω)` on the Born-ranked subspace — the **infinite-shot limit** of the protocol, not a finite-shot run. | `sampled_akw_L8.json` |
| `sampled_honest.py`, `honest_sampling_scaling.py` | The measured **finite-shot** sweeps that Table IX and Fig. 3(c) report. | `sampled_honest.json`, `honest_sampling.json` |
| `spectral_validation_max.py` | Max-level validation of the sampled `A(ω)` (Hubbard L=6, U/t=8). | `method_max.json` |
| `charge_gap_ed.py` | μ⁺, μ⁻ and the Mott gap Δ of the half-filled ring, L=6…12, by sector diagonalization. | `charge_gap.json` |
| `gap_scaling.py` | Finite-size scaling of the **two** gaps (spin and charge) with a Heisenberg control. | `gap_scaling.json` |
| `pole_structure.py` | Pole census behind Table XI. | `pole_structure.json` |
| `eta_sweep.py` | Resolution sweep behind the η prices of Sec. VII. | `eta_sweep.json` |
| `make_sqw_edges.py` | Generator of `figs/sqw_edges.dat` (band edges + peak). Until 2026-09-18 that file had **no generator anywhere**. | `figs/sqw_edges.dat`, `sqw_edges_provenance.json` |
| `sqw_field.py`, `recolor_akw_v2.py` | Perceptual "hot" rasters for the `S(q,ω)` / `A(k,ω)` colormaps. `recolor_akw_v2.py` **cannot run here** — neither of its inputs is deposited (§13). | `build/` |

### The certificate (Theorem B and its falsification)

| File | Purpose | → |
|---|---|---|
| `leakage_certificate.py` | The certificate itself, with a self-test. | `leakage_certificate_selftest.json` |
| `leakage_certificate_suite.py` | The 606-subspace sweep behind Sec. V and Fig. 17. Writes the frontier files under a **computed** name, which is why `cert_frontier_*.json` appear orphaned. | `cert_akw.json`, `cert_stress.json`, `cert_teqsci.json`, `cert_frontier_6-8.json`, `cert_frontier_10.json` |
| `support_witness.py` | Support witness behind Sec. VI. | `support_witness.json` |
| `apsg_witness.py` | The **provable** geminal witness: `F_k` large yet χ=2, `\|S\|=2^K`. | `apsg_witness.json` |

### Molecular suite and the resource theory

| File | Purpose | → |
|---|---|---|
| `n19_suite.py` | Nineteen molecules: FCI-verified `F₁`, `\|S\|`, χ at two geometries. **Needs `pyscf`.** | `n19_suite.json` |
| `n19_spectral.py` | Suite-wide sampled `A(ω)` reconstruction. **Needs `pyscf`.** | `n19_spectral.json` |
| `n2_hero_data.py` | N₂ dissociation: spectrum and resources together; the `F₁=2Nᵤ` identity. **Needs `pyscf`.** | `n2_hero.json`, `figs/hero_*.dat` |
| `cost_vs_entanglement.py` | The 30-point Hubbard sweep: the cost is governed by χ, not `F₁`. | `cost_vs_ent.json` |
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
`make_decoupling_native.py` (Figs. 6 **and** 7), `make_method_fig_max.py` (Fig. 2),
`make_akw_sampled_honest_fig.py` (Fig. 3 — ⚠ absolute paths, §25),
`make_hardware_hero.py` (Fig. 14), `make_noise_recovery_native.py` (Fig. 16),
`make_hero_fig.py` (Fig. 13, into the cwd), `make_noise_fig.py` (Fig. 15, into the cwd),
`make_heron_fig.py` (writes `figs/heron*.dat`), `make_table.py` (Table XIII, into `rebuild/`).

**Second writer, not the guarded one:** `make_resource_master_fig.py` writes
`figs/fig_resource_master_native.tex`, the same file `make_decoupling_native.py` writes.

**Hazard:** `make_scaling_fig.py` emits a *different* figure and would destroy panel (c) of Fig. 9. It
is deliberately out of `make figures` and out of `check_figures.py` (§3).

**Dead** — generators of the five figures v3 withdrew, kept as the record of how they were made:
`make_bench_fig.py` (`fig:bench`, and it would recreate `figs/bench_U*.dat`),
`make_gflow_fig.py` (`fig:gflow`), `make_gallery_native.py` (`fig:gallery`),
`make_magic_suite_fig.py` (`fig:molsuite`).

**Absent, and named by their own output:** `make_fig_gapscaling.py` (Fig. 10) and `build_n3.py` with
`fix_n3_caption.py` / `fix_n3_fig5_series.py` (Fig. 17). See `KNOWN_DISCREPANCIES.md` §27.

---

### `src/frontier/` — 24 scripts, the C3 certificate-frontier sweep (deposited 2026-09-19)

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
| **E — published fractions** | `pubsum.py`, `run_L14.py` | the four points at the fractions the manuscript actually publishes, and the L=14 point at FR=0.08 |
| **F — post-processing** | `calib_oos.py`, `null_ws4.py`, `null_ws5.py`, `null_ws6.py`, `null_ws7.py` | zero new compute: the out-of-sample test of the calibration, and the null a referee will propose — does `1−w_S` locate the error as well as the certificate? |

---

## `data/c3_frontier/` — 82 files, the C3 sweep itself

54 point files, the 204-row grid they aggregate to, three gate outputs, three proof ladders, the
verdict, and four sub-directories (`published_fraction/`, `adversarial/`, `null/`). Every number in
Sec. V comes from here. Documented file by file in [`C3_FRONTIER.md`](C3_FRONTIER.md) §1; `data/` is
excluded from the arXiv package, so the deposit costs the submission nothing.

---

## `data/` — 37 authoritative datasets

| file | feeds |
|---|---|
| `akw_lanczos_L12.json` | Fig. 4 `fig:lattice` (exact `A(k,ω)`, L=12) |
| `sampled_akw_L8.json` | Fig. 3(a,b) `fig:akwsampled` |
| `sampled_honest.json`, `honest_sampling.json` | Fig. 3(c) and Table IX (the finite-shot sweeps) |
| `method_max.json` | Fig. 2 `fig:method` + `figs/aw_method.dat` |
| `sqw_L12.json` | Fig. 11 `fig:sqw` + `figs/sqw_edges.dat` |
| `spinqw_L12.json` | Fig. 12 `fig:spin` + `figs/spinqw_edges.dat` |
| `gap_scaling.json` | Fig. 10 `fig:gapscaling` |
| `charge_gap.json` | the Δ = 4.97 t of Sec. VIII and Table X |
| `pole_structure.json` | Table XI `tab:poles` |
| `cost_vs_ent.json` | Fig. 6(b,c) `fig:decoupling` |
| `n19_suite.json` | Fig. 6(a), Table XIII, Table XVI |
| `n19_spectral.json` | Table XIII (the reconstruction columns) |
| `resource_master.json` | Fig. 7 `fig:master` |
| `stats_resource.json` | Fig. 6(d), Appendix F, Table XV |
| `apsg_witness.json` | Fig. 5 `fig:witness` |
| `support_witness.json` | Sec. VI (the support witness) |
| `ladder_vs_chain.json` | Fig. 8 `fig:ladder` |
| `scaling_data.json` | Fig. 9 `fig:scaling` |
| `eta_sweep.json` | Sec. VII (the resolution price) |
| `cert_akw.json`, `cert_stress.json`, `cert_teqsci.json`, `cert_frontier_6-8.json`, `cert_frontier_10.json`, `leakage_certificate_selftest.json` | Sec. V, Tables III–V, and Fig. 17 via the absent `build_n3.py` |
| `hw_lucj_n2_result.json` | ⭐ Fig. 14(b) — **job `da125f2ein7c73bcsqs0`** (authoritative N₂ run) |
| `heron_spectral.json` | Fig. 14(a) — the `ibm_fez` L=6 `A(ω)` run → `figs/heron_hot.dat` |
| `noise_spectral.json` | Fig. 15 `fig:noise` |
| `molecular_noise.json`, `molecular_noise_sweep.json` | Fig. 16 `fig:noiserec` |
| `n2_hero.json` | Fig. 13 `fig:hero` |
| `headtohead_ms.json` | **no figure in v3** — `fig:bench` was withdrawn; the data stays |
| `gflow.json`, `amortized_recovery.json` | **no figure in v3** — `fig:gflow` and `fig:amort` were withdrawn by Sec. 9.4. `gflow.json` additionally has **no generator anywhere in `src/`** (§11.2 item 7) |
| `sqw_edges_provenance.json` | the provenance record of `figs/sqw_edges.dat` |
| `sampled_akw_figure_numbers.json` | the number→source manifest of Fig. 3 |

Backups (`*_backup.json`) are not the plotted data and are `.gitignore`d.

---

## Orphan inventory

Produced by `python src/check_provenance.py -v`, 2026-09-19, over the 80 artefacts that are *supposed*
to have a producer (`paper/figs/*.dat|pdf|json`, machine-written `paper/figs/*.tex`, `data/*.json`,
`docs/img/*.png`). Hand-written LaTeX is excluded — it has no producer by design.

The tiers use two facts that need no guessing: does any deposited **code** name the file, and does the
**manuscript** reach it. Writing a file down in a document does not give it a producer, so documenting
an orphan does not remove it from tier 1.

### Tier 1 — no deposited code names it, and the manuscript does not use it (14)

| file | why it is here |
|---|---|
| `paper/figs/n3_fig5_eta.dat`, `n3_fig5_frac.dat` | the 64 Fig.-3 rows that `caption_thm1iii_violation.tex` says were summarised into Fig. 17's whiskers by `fix_n3_fig5_series.py`, which is not deposited. The whisker values are literals in the fragment. |
| `data/cert_frontier_6-8.json`, `cert_frontier_10.json` | written by `leakage_certificate_suite.py:697` under a **computed** filename (`"cert_frontier_%s.json" % "-".join(...)`), so the literal names appear in no script. Only a comment in `paper/sec_5_body.tex:70` mentions them, and nothing reads them back. |
| the ten `docs/img/*.png` | hand-exported README thumbnails; no script produces them. Two — `bench.png` and `witness.png` — are not embedded anywhere either, and `bench.png` is a thumbnail of a figure v3 withdrew. |

### Tier 2 — code names it, the manuscript does not use it (43)

Mostly datasets whose figure was withdrawn or whose only remaining reader is the guardian:
`gflow.json`, `amortized_recovery.json`, `headtohead_ms.json`, the six certificate files,
`eta_sweep.json`, `support_witness.json`, `sampled_akw_figure_numbers.json`,
`figs/figure_numbers.json`, and the seven `.dat` tables read only by the standalone sources of §5
(`hero_aw.dat`, `hero_res.dat`, `heron.dat`, `heron_hw.dat`, `noise.dat`, `spinqw_edges.dat`,
`sqw_edges.dat`). **None of these is a defect** — a dataset with a producer and a guardian is doing
its job — but none of them is load-bearing for the v3 manuscript either.

### Tier 3 — the manuscript typesets it and no deposited code names it (13)

The six PDFs (`fig_akw_native.pdf`, `fig_circuit.pdf`, `fig_hero.pdf`, `fig_noise_score.pdf`,
`fig_spinqw.pdf`, `fig_sqw.pdf`) and the seven `n3_*.dat` tables Fig. 17 plots. Their *numbers* are in
committed JSONs; their *renderings* are reproducible by nothing in this deposit.

### `paper/figs/*.dat` — who reads them

| file | read by a fragment **in this deposit**? | read by | written by |
|---|---|---|---|
| `aw_method.dat` | **yes** | `fig_method_native_frag.tex` | `make_method_fig_max.py` |
| `heron_hot.dat` | **yes** | `fig_hardware_hero_frag.tex` | `make_hardware_hero.py` |
| `n3_pub_eta_hi/lo.dat`, `n3_pub_frac_hi/lo.dat`, `n3_thmB_eta/frac.dat`, `n3_inf_eta/frac.dat` | **yes** (8 files, 2432 lines) | `fig_thm1iii_violation_native.tex` | **nothing deposited** — `build_n3.py`, §27 |
| `n3_fig5_eta.dat`, `n3_fig5_frac.dat` | no | nothing, anywhere | **nothing deposited** |
| `hero_aw.dat`, `hero_res.dat` | no | a standalone `fig_hero.tex` that is **NOT DEPOSITED** — note it is *not* `paper/carried/fig_hero.tex`, which is the float wrapper — but which `make_hero_fig.py` rewrites into the cwd | `make_hero_fig.py` |
| `noise.dat` | no | a standalone `fig_noise_native.tex`, **NOT DEPOSITED**, rewritten into the cwd by its generator | `make_noise_fig.py` |
| `heron.dat`, `heron_hw.dat` | no | `work/notes/fig_heron_native.tex`, a **superseded** version of the hardware figure — `work/` is the author's sibling working tree and is **NOT in this deposit** | `make_heron_fig.py` |
| `spinqw_edges.dat` | no | `work/notes/fig_spinqw.tex` — **not in this deposit** | `spin_lanczos.py` |
| `sqw_edges.dat` | no | `work/notes/fig_sqw.tex` — **not in this deposit** | `make_sqw_edges.py` (since 2026-09-18; before that, **nothing**) |

---

## `_superseded/` — 24 tracked files, kept so the history is legible

| subtree | contents |
|---|---|
| `v2_paper/` | the arXiv v2 manuscript: `main.tex`, `main.pdf`, `introduction.tex`, `method.tex`, `resource.tex`, `results.tex`, `hardware.tex`, `discussion.tex`, `conclusion.tex`, `theorem_b2.tex`, `table_molecules.tex`. Every one of them was deleted from `paper/` in `606f969` and is preserved here, so the deletion commit does not lose the trail. (`main.bbl` is on disk but untracked: `.gitignore`'s `*.bbl` excludes it, and the v2 bibliography was hand-written anyway.) |
| `v2_figs/` | `fig_akw_sampled_native.tex` and `make_akw_sampled_fig.py` — the Fig.-5 fragment retired on 2026-09-19 because it labelled an infinite-shot curve "sampled (85% sector)", and its generator. Replaced by `fig_akw_sampled_honest.tex`. |
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
| `00_Reproduce_Everything.ipynb` | **Master reproduction** — one narrated pipeline regenerating every figure and number from the committed `data/*.json`. Assembled by `src/build_repro_notebook.py`. **Needs `pyscf`; not verified in this pass.** |
| `HW_LUCJ_N2_Heron_READY.ipynb` | N₂ energy on `ibm_marrakesh`: LUCJ circuit, sampling, self-consistent recovery, 8-seed error-bar cell. Credentials scrubbed. **Its saved cell outputs are an older run** (`d9qe731dsedc73af67d0`) and are not the paper's numbers. |
| `Spectral_Heron.ipynb` | `A(ω)` on `ibm_fez`: the L=6 Hubbard spectral run (full-sector coverage). Credentials scrubbed. |

Both hardware notebooks have the API **token and instance CRN removed** (placeholders), and neither
ships raw counts: they fetch them from the IBM Quantum service by job id.

---

## Support files

| File | Purpose |
|---|---|
| `README.md` | Front page, figure gallery, quick start. |
| `CITATION.cff` | Machine-readable citation. |
| `LICENSE` | MIT (code) + CC-BY-4.0 (paper, figures and data). |
| `Makefile` | `make verify` / `check-figures` / `figures` / `data` / `paper` / `reproduce` / `ci`. Measured results per target in `KNOWN_DISCREPANCIES.md` §10. |
| `requirements.txt` | Python dependencies. |
| `.gitignore` | Includes `paper/*Notes.bib` (a TeXstudio artefact that would otherwise be committed) and `ckpt/`, `rebuild/`, `build/`, `*.log`. |
| `.github/workflows/ci.yml` | Two jobs: `python src/verify.py` and `make check-figures`. Tracked since 2026-09-18 — **it has still never executed**, because the commit carrying it has not been pushed (`KNOWN_DISCREPANCIES.md` §24). |
| `docs/REPRODUCE.md` | Figure/number → script → command, and what a clean clone cannot reach. |
| `docs/FIGURE_PROVENANCE.md` | Figure → source → generator → data → job id, plus the stale-file list. Checked by `src/check_provenance.py`. |
| `docs/KNOWN_DISCREPANCIES.md` | Every place where running a committed script does **not** reproduce the deposit, repaired or not. **Opens with a banner about `paper/arxiv-submission.tar.gz`, which must never be uploaded as v3.** |
| `docs/DATA_AVAILABILITY.md` | The Data Availability Statement drafted for Physical Review A, plus the audit licensing each of its sentences. |
| `docs/ARXIV_SUBMISSION.md` | Submission metadata and the step-by-step for posting. |
| `docs/img/*.png` | Ten README thumbnails; hand-exported, no generator (tier-1 orphans above). |

### Present on disk, not in the deposit (`.gitignore`d)

`ckpt/` (resumable checkpoints of `scaling_lanczos_mf.py`), `rebuild/` (scratch output of
`make_table.py`), `build/` (raster and log by-products), `src/__pycache__/`,
`paper/main.{aux,log,out}`, and two stray root logs `texput.log` / `x2.log` left by a TeX run.
