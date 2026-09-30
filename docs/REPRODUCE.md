# REPRODUCE — every figure and every number, traced to source

This is the referee's map, and it is meant to be **literally true**. Where a clean clone reaches the
paper's numbers, the row says how; where it cannot — because a figure step is not deposited, the
IBM jobs are closed, or `pyscf` is needed — the row says that instead, and says why.

*Rewritten 2026-09-19 for the v3 manuscript. The previous version described twenty figures and named
generators for five that no longer exist. Figure and table numbers updated 2026-09-26 to those printed
in v3 as submitted: 17 figures (1–9 in the main text, S1–S8 in the Supplemental Material) and 17
tables (I–V and S1–S12).*

All scripts live in [`../src/`](../src/) and all committed data in [`../data/`](../data/), with the
environment of [`requirements.txt`](../requirements.txt). The two-step pattern is:
`python src/<compute>.py` → writes `data/<data>.json` → `python src/make_<fig>.py` → writes
`paper/figs/<fragment>.tex` → `make paper`.

> **Where the output goes (changed 2026-09-18).** Every compute script resolves its output path from
> **its own location on disk**, so `data/` means *this clone's* `data/` no matter what directory you
> are standing in. Pass `--out PATH` to send it somewhere else, `--in PATH` to read from somewhere
> else. Until 2026-09-18 eighteen of these scripts wrote to `/w/`, the working directory of the
> container the original runs were made in: **following these instructions did not repopulate
> `data/`, and this document was not true.** See [`KNOWN_DISCREPANCIES.md`](KNOWN_DISCREPANCIES.md)
> §2 for the repair and the evidence that it works.
>
> (Repaired 2026-09-19, §25: every generator now derives its paths from `__file__`.)

---

## The four checks, and what each is worth

| command | what it proves | measured 2026-09-19 |
|---|---|---|
| `python src/verify.py` (`make verify`) | recomputes the physics — `F₁(L=6)` for the open chain (9.6047) and the periodic ring (9.3851), the geminal witness `F₁=F₂=4K`, `\|S\|=2^K`, `χ=2` — from first principles, then checks the deposited `.json`, `.dat` and `.tex` (the abstract included) against it. Since 2026-09-19 its section **(10b)** also opens the 54 C3 point files and re-derives every column of the 204-row frontier grid, including **Theorem B itself on all 204 rows** | **PASS**, 868 checks, 0 failures, **9 884 numeric assertions**, 5 xfail, 0 xpass, 2 skips *(2026-09-19)* |
| `python src/check_figures.py` (`make check-figures`) | re-runs eight generators in a scratch tree — the six of 2026-09-19 plus the two recovered into `src/recovered/` on 2026-09-26 — and diffs 19 artefacts against the committed ones | **PASS**, **19 of 19** *(2026-09-26; one comment line of the Fig. 6 fragment, a provenance timestamp, is normalised and the run says so — §29)*. Until the §25 repair one artefact was compared with a copy of itself and the run rewrote four files in the working tree; both are measured, with controls, in §25 |
| `python src/check_coherence.py` (`make check-coherence`) | expands `\input` from `paper/main.tex` and checks the manuscript against **itself**: the vacuity range, the number of published sizes and the largest sector are re-derived from `data/c3_frontier/published_fraction/` and must be stated the same way in the abstract, Sec. I, Sec. V, Sec. IX and Sec. X; the hardware shot budget must be one number in all four places; and the data statement (Sec. IX D of v3 as printed) must quote the guardian counts `verify.py` actually prints. Each check has an absence half and two synthetic controls | **PASS**, 10 checks, **20 of 20 controls fire** *(2026-09-19; written that day, after seven contradictions between sections survived a green run of the other three)* |
| `python src/check_provenance.py` | derives the figure inventory from `paper/main.tex` and fails if this document, `FIGURE_PROVENANCE.md` or the `main.tex` float census no longer describes the manuscript | **PASS** — 61 source files reached, 34 floats (17 figures + 17 tables), all starred *(2026-09-28, with `app_proof.tex`; 60 on 2026-09-26, 50 on 2026-09-19)* |

> **`--fast` certifies nothing.** It skips the `L=10` diagonalisations and the checks that depend on
> them, and it exits 3 even when nothing fails. Use it while editing; never quote it.

> **The `[XFAIL]` lines are real.** They are registered open defects (`KNOWN_OPEN`), open claims
> (`CLAIMS_OPEN`) and one declared coverage gap (`COVERAGE_GAP_V3`). They print on every run, they do
> not set the exit code, and they are the honest agenda. An `[XPASS]` means a registered defect has
> been repaired and the registry entry should now be deleted so the slot becomes a hard check.

> **No `[XPASS]` is expected, and the count of `[XFAIL]` is five.** `paper/app_carried_repro.tex`
> now prints the correct `R_eq(H₂O) = 0.958`; `KNOWN_OPEN['tab_repro.H2O.R_eq']` was deleted and
> that slot is a hard check again. The fifth registered defect is new on 2026-09-19 and belongs to
> the C3 deposit: `data/c3_frontier/ladder_L4.json` was written four minutes before
> `src/frontier/fron_lib.make_grid` widened its integration window, so re-running it today gives a
> different ω grid. The registry pins both values; the measured consequence — L=4 rows of
> `verdict_frontier.json` move by ≤ 2.8 % relative, L=6 and L=8 unchanged, no manuscript number
> affected — is in [`C3_FRONTIER.md`](C3_FRONTIER.md) §4.

---

## Figures — all seventeen

Figure numbers are those **printed in v3 as submitted**: 1–9 in the main text, S1–S8 in the
Supplemental Material. The `\label` is the stable key; documents written between 2026-09-19 and
2026-09-25 number the figures 1–17 continuously, and the conversion table opens
[`KNOWN_DISCREPANCIES.md`](KNOWN_DISCREPANCIES.md). `guard` = regenerated byte-for-byte by
`check_figures.py`. Full provenance, including the data behind every panel, is in
[`FIGURE_PROVENANCE.md`](FIGURE_PROVENANCE.md).

| printed (v3) | `\label` | Paper object | Compute (exact) → data | Figure generator | guard |
|:--:|---|---|---|---|:--:|
| 1 | `fig:circ` | Sampling-primitive circuit schematic | — (schematic, no data) | hand-drawn quantikz; standalone source `paper/figs/src/fig_circuit_native_FIXED.tex` (deposited 2026-09-26; recompiles pixel-identically, §29) | — |
| 2 | `fig:method` | Method validation `A(ω)` + Krylov convergence, L=6 | `spectral_validation_max.py` → `method_max.json` | `make_method_fig_max.py` → fragment + `figs/aw_method.dat` | ✔ |
| 3 | `fig:akwsampled` | `A(k,ω)` on a Born-ranked subspace, L=8, + the shot cost | `sampled_akw.py` → `sampled_akw_L8.json`; `sampled_honest.py` / `honest_sampling_scaling.py` → the finite-shot sweeps; `fig3_finite_shot.py --T 2600000` and `--T 5200000` → `fig3_finite_shot_T*.json`, the caption's finite-shot numbers (`--summary` prints them) | `make_akw_sampled_honest_fig.py` | ✔ |
| 4 | `fig:lattice` | Exact `A(k,ω)`, L=12 Mott map | `akw_lanczos.py` → `akw_lanczos_L12.json` | `recolor_akw_v2.py`; **neither of its two inputs is deposited**, so it cannot run here. The figure ships as `figs/fig_akw_native.pdf`; its standalone source and the recoloured raster it embeds are in `paper/figs/src/` (2026-09-26, recompiles pixel-identically) | — |
| 5 | `fig:scaling` | Scaling to 28 qubits, three panels | `scaling_data.py` (L≤8) / `scaling_lanczos_mf.py` (L≤14) → `scaling_data.json` | ⚠ **none usable.** `make_scaling_fig.py` emits a *different* two-panel figure and would destroy panel (c); it refuses to run without `--i-know-this-drops-panel-c --out <path>` (§3) | — |
| 6 | `fig:gapscaling` | The two gaps of the ring scale oppositely | `gap_scaling.py` → `gap_scaling.json` (spin gap by Lanczos with full reorthogonalisation; charge gap by sector ED; Heisenberg control) | `src/recovered/make_fig_gapscaling.py` — recovered and deposited 2026-09-26 (§27, §29); the fragment regenerates except one provenance-timestamp comment line, and the caption was hand-edited after generation | ✔ |
| 7 | `fig:sqw` | `S(q,ω)` charge structure factor | **`python src/sqw_lanczos.py 12 --eta 0.18`** → `sqw_L12.json` (the flag is required: the default is 0.20 and the script **refuses to run** without it); edges by `make_sqw_edges.py` → `figs/sqw_edges.dat` | `sqw_field.py` (writes the field in viridis; the hot-coloured raster the figure embeds was recoloured by a step that is not deposited) + the standalone `paper/figs/src/fig_sqw_fixed.tex` (deposited 2026-09-26, recompiles pixel-identically) | — |
| 8 | `fig:spin` | `S^zz(q,ω)` spin structure factor | `spin_lanczos.py` → `spinqw_L12.json`, `figs/spinqw_edges.dat` | the standalone `paper/figs/src/fig_spinqw_fixed.tex` (deposited 2026-09-26, recompiles pixel-identically); no deposited script writes its raster `spinqw_field.png` | — |
| 9 | `fig:heron` | **IBM Heron hardware** — `A(ω)` on `ibm_fez` + N₂ energy on `ibm_marrakesh` | the two notebooks; the deposited post-processing is `heron_spectral.json` and `hw_lucj_n2_result.json` (job **`da125f2ein7c73bcsqs0`**) | `make_hardware_hero.py` → fragment + `figs/heron_hot.dat` | ✔ |
| S1 | `fig:witness` | The geminal (APSG) decoupling witness | `apsg_witness.py` → `apsg_witness.json` | none — `figs/fig_witness_native.tex` is hand-maintained against that JSON | — |
| S2 | `fig:decoupling` | `F₁`–`\|S\|` decoupling, stratified | `n19_suite.py`, `cost_vs_entanglement.py`, `stats_resource.py` → `n19_suite.json`, `cost_vs_ent.json`, `stats_resource.json` | `make_decoupling_native.py` | ✔ |
| S3 | `fig:master` | The resource thesis over 74 deposited rows | (aggregated) → `resource_master.json` | `make_decoupling_native.py` (also writes this one). `make_resource_master_fig.py` is a **second writer** of the same file | ✔ |
| S4 | `fig:ladder` | Chain vs two-leg ladder at matched `F₁` | `ladder_vs_chain.py` → `ladder_vs_chain.json` | none — `figs/fig_ladder_native.tex` is hand-maintained | — |
| S5 | `fig:thm1iii-violation` | Theorem 1(iii) of v1–v2 fails on the subspaces the paper uses | `leakage_certificate_suite.py` → `cert_akw.json`, `cert_stress.json`, `cert_teqsci.json`, `cert_frontier_*.json` | `src/recovered/build_n3.py` on `data/thm1iii_violation/` — recovered and deposited 2026-09-26; all ten tables regenerate byte-for-byte. `src/verify.py` still registers the figure as its one declared coverage gap (no check there anchors the ratios) | ✔ (tables) |
| S6 | `fig:hero` | N₂ dissociation hero (`F₁ = 2Nᵤ`) | `n2_hero_data.py` → `n2_hero.json`, `figs/hero_*.dat` (**needs `pyscf`**) | `make_hero_fig.py` — writes `fig_hero.tex` into the **cwd**, so the figure is regenerable in two steps; the standalone source of the committed PDF is `paper/figs/src/fig_hero_STANDALONE.tex` (2026-09-26, recompiles pixel-identically) | — |
| S7 | `fig:noise` | Noise robustness (bit-flip, L=6) | `noise_spectral.py` → `noise_spectral.json`, `figs/noise.dat` | `make_noise_fig.py` — writes `fig_noise_native.tex` into the **cwd**; the standalone source of the committed PDF is `paper/figs/src/fig_noise_score_STANDALONE.tex` (2026-09-26, recompiles pixel-identically) | — |
| S8 | `fig:noiserec` | Noise-assisted energy across chemistry | `molecular_noise_energy.py`, `molecular_noise_sweep.py` → `molecular_noise*.json` (**need `pyscf`**) | `make_noise_recovery_native.py` | ✔ |

**Tables.** Table S12 (`tab:molecules`) is regenerated by `make_table.py` into `rebuild/`, a scratch
directory; the **shipped** table is `paper/table_molecules_v3.tex`, whose data rows are identical and
whose caption carries hand-added text. The other sixteen tables are hand-maintained inside their
section files (Table II, `tab:moments`, is rebuilt by `src/moments_table.py` since 2026-09-25) and are
checked value-by-value by `src/verify.py`.

**Withdrawn in v3, and no longer figures:** the selector benchmark (`headtohead_ms.json`), the
molecular gallery and magic suite (`n19_suite.json`), and the
amortized-recovery study (`amortized_recovery.json`). Their data and their generators are still in the
deposit; their bodies are in `_superseded/v3_retired/figs/`. The v3 manuscript does not typeset
them and does not discuss them: its version history is kept out of the body, in the arXiv Comments
field, `README.md` and `CITATION.cff`.

---

## Key numbers → where they come from

Unless the row says otherwise, every entry below is checked by `src/verify.py` against the file named,
on every run. The rows added on 2026-09-26 cover the headline numbers of the v3 abstract and name
their guard explicitly, including where there is none.

| Number in the paper | Value | Source |
|---|---|---|
| The bound is vacuous at every operating fraction of the resource scan, at all five sizes (abstract, Secs. I, V, X; Sec. IX C stopped restating it on 2026-09-30); Fig. 3, at 85% of its sector, is the certified exception (abstract, Secs. I, V, X) | vacuity factor **1.66–8.61**; largest sector **D = 10 306 296** (L=14) | `data/c3_frontier/published_fraction/PUB_rows.json` (16 rows, L=6–12; minimum 1.659 at L=6, η=0.25) and `L14_FR008.json` (maximum 8.613 at L=14, η=0.10, at its final Lanczos depth, 250). Guard: `python src/check_coherence.py` re-derives the range and the largest sector from these files and requires every section that states them to agree |
| The bound crosses the trivial one at a true relative error of … (abstract, Sec. V) | **1.2–7.3×10⁻³, median 2.8×10⁻³** over 16 `(L,η)` cells | `data/c3_frontier/null/CALIB_OOS.json` (1.227×10⁻³ … 7.268×10⁻³, median 2.746×10⁻³), written by `python src/frontier/calib_oos.py` from the frontier block of `data/c3_frontier/C3_grid.json`, which it regenerates exactly (last table of this page). Guard: `verify.py` §(10b) re-derives `C3_grid.json`; it does **not** assert these three numbers |
| A two-constant rule in `1−w`, fitted out of sample, is tighter on 11 of 14 splits (abstract, Secs. I, V, X) | **11** null wins, **3** certificate wins | `python src/frontier/null_ws7.py` (re-run 2026-09-26, ~2 min: 11 of the 14 train/test splits go to a fitted `w_S` rule, 3 to the certificate). It **prints to stdout only; no output file is deposited**, and no guardian asserts the count |
| On the probe's support, where `w = 1` exactly, the error is still … (abstract, Sec. VI) | relative `L₁` **0.43** at L=6 and **0.35** at L=8, η=0.18 t | `data/support_witness.json` (`src/support_witness.py`); asserted by `verify.py`, together with the 2237-determinant support at L=8 |
| Collapse identity `F₁ = 4 tr[γ(1−γ)] = 2Nᵤ` | to `<10⁻¹³` | `n19_suite.py`, `n2_hero_data.py` |
| Open-chain `F₁(L=6)` (the resource families) | 9.6047 | recomputed from first principles inside `verify.py` |
| Periodic-ring `F₁(L=6)` (the scaling figure) | 9.3851 | *a different system* — `scaling_data.json`; see `KNOWN_DISCREPANCIES.md` §4 |
| Subspaces in the certificate sweep (Fig. S5, Secs. III D and S9) | 621 pooled, 606 live, the weight-only bound violated on **469**, and on all **363** with `w_S ≥ 0.99`, the mildest by a factor **2.00**; below `w_S = 0.99` it survives on 137 of 243 | **the deposited `data/cert_*.json`, since 2026-09-28 (late).** Until that day the plotted pool read its three certificate scans from an earlier run of the same scans (kept in `data/thm1iii_violation/cert_earlier_run/` as a record), which gave 468 / 2.10 / 138; the two runs differ on 28 XXZ and t–V stress subspaces. `python src/recovered/build_n3.py` now prints the numbers above from the deposited files, `src/recovered/n3_fig5_whiskers.py` writes the whiskers, `check_figures.py` guards the ten tables, the fragment and the summary, and `verify.py` section 9.6c re-classifies the pool independently and checks every printed count and every plotted ratio — the coverage gap it used to declare is closed (§33) |
| **Everything in Sec. V** — Tables III–V, the frontier, the calibration, the convergence census | 204 evaluations, 187 CONVERGED / 6 MARGINAL / 11 NOT-CONVERGED | `data/c3_frontier/` (54 point files → `C3_grid.json`), re-derived in `verify.py` §(10b). **Deposited 2026-09-19**; before that date this row read *not deposited* |
| At L=12 and the **published** fraction 0.18, the certificate is vacuous | leak/trivial = 7.51, 4.86, 4.00, 2.81 at η = 0.10, 0.15, 0.18, 0.25 | `C3_grid.json`, asserted in `verify.py` §(10b) |
| Theorem B holds on every measured row | `(1−w_S) ≤ rel-L₁ ≤ min{leak, trivial}`, **0 violations in 204 rows** | `C3_grid.json`, checked in `verify.py` §(10b) |
| The cost of the C3 sweep | **2 873.9 s = 0.80 h**, peak **1.64 GiB** over 51 subspaces | recomputed from `data/c3_frontier/` on every `verify.py` run and compared with the stanza in [`C3_FRONTIER.md`](C3_FRONTIER.md) §2. The pre-run estimate was 20–40 h and 9 GB |
| Within-stratum Spearman `ρ(F₁,\|S\|)` | exactly −1 in all six strata | `cost_vs_ent.json`, recomputed in `verify.py` |
| Exact one-sided permutation p | `3.3×10⁻¹³` | `stats_resource.json` → `hubbard_stratified` |
| Mott gap Δ = μ⁺ − μ⁻ at L=12 | 4.97 t | `charge_gap.json` (and independently `gap_scaling.json`) |
| Hardware N₂ energy (`ibm_marrakesh`) | **0.59 ± 0.14 mHa** (8-seed mean ± s.d.) | `hw_lucj_n2_result.json`, job `da125f2ein7c73bcsqs0` |
| Noiseless statevector, same circuit | 29.5 mHa | same file (`e_sim`) |
| Hardware `A(ω)` (`ibm_fez`) coverage | `\|S\| = 300/300` (exact by coverage) | `heron_spectral.json` |

**Hardware provenance:** the authoritative N₂ run is job `da125f2ein7c73bcsqs0` (8-seed recovery
sweep). `hw_lucj_n2_result.json` stores the full per-seed × per-iteration trajectories
(`e_hw_seed_hist`, `dE_hist_seeds`) so the per-point error bars cannot be lost again. The **saved cell
outputs of `HW_LUCJ_N2_Heron_READY.ipynb` are an older single-seed run** and are not the paper's
numbers.

---

## Compute engines (shared)

- **`akw_lanczos.py`** — sector-basis (up-string ⊗ dn-string) Hubbard engine: sparse ground state,
  Haydock continued fraction for `G(z)`, `c†_k`/`c_k` maps. Reused by `sampled_akw.py`,
  `scaling_data.py`, `charge_gap_ed.py`.
- **`scaling_lanczos_mf.py`** — fully matrix-free, RAM-staged, resumable (checkpoints to `ckpt/`) to
  reach L=14 (28 qubits) on a memory-limited machine.
- Every `eigsh` call site passes a fixed start vector (`default_rng(20260918)`). Deposited values were
  produced by *unseeded* runs, so a re-run agrees to ARPACK tolerance (~1e-13), not in the last digit
  — and from now on it agrees with itself (`KNOWN_DISCREPANCIES.md` §14).

---

## What a clean clone can and cannot reach — measured, not asserted

Environment: Windows workstation, Python 3.11, numpy + scipy, **no `pyscf`**, no TeX in the loop, no
`/w` directory of any kind. Rows marked *(2026-09-19)* were re-run today; the others are carried from
the 2026-09-18 pass and say so.

### Reaches the paper's numbers

| command | result |
|---|---|
| `python src/verify.py` | **PASS** — 1280 checks, 0 failures, 14 563 numeric assertions, 4 registered xfail, 0 xpass, 2 skips *(2026-09-30, final review: +66 checks and +66 assertions, verify.py checks the Fig. 3 caption's panel (b) mean and maximum, window and branch thresholds (section 9.6b), the printed noise values of Figs. S7-S8, Sec. S8 and Sec. VIII (section 9.9b) and every restated Fig. 3 certificate number in the abstract and Secs. I, V, X and S3 (section 9.9c); 1214 / 14 497 earlier on 2026-09-30: +2 checks and +2 assertions, verify.py section 9.3 checks every printed exact permutation p (Sec. I, Sec. VI, Fig. 6 caption) against `exact_perm_p_one_sided_F1` at 1e-15, unified on the 3.35e-13 rounding; 2026-09-28, late, closing the review, +278 checks and +2 821 assertions: section 9.6c (the free-fermion supports, the spin-peak deviations and the Fig. S5 pool re-classified from the deposited certificate rows), the chi_max of the witness, the LF hash of the Fig. 3 generator and the L=14 re-ranking record; the Fig. S5 coverage gap is closed, xfail 5 → 4, §33; 934 / 11 674 earlier that day, after the Fig. 3 finite-shot run was deposited and its caption numbers checked from its data, +60 checks and +187 assertions; 874 / 11 487 earlier the same day, after the Fig. 3 frequency grid was widened to [-9t, 17t], +1 596 assertions; 874 / 9 891 earlier the same day, after two phrase guards on Sec. S8; 872 / 9 889 on 2026-09-26, after the hardware bit-order check; 868 / 9 884 on 2026-09-19, after section (10b) and the `ladder_L4` open defect arrived with the C3 frontier deposit; the 807 / 7 454 / 4 reading of earlier the same day is superseded)* |
| `python src/check_provenance.py` | **PASS** — the documentation describes the manuscript in the tree *(2026-09-19)* |
| `python src/check_tarball.py` | runs; reports 9 byte-identical / 17 differ / 21 bundle-only against the frozen v2 bundle *(2026-09-19)* |
| `python src/check_figures.py` | **PASS** on 8 artefacts, **all 8 genuinely guarded** since the §25 repair — verified by seeding each artefact with a sentinel line and checking it is gone after the run, and by snapshotting sha256 + mtime of all 231 files to confirm the run writes nothing into the repository *(2026-09-19)*. **PASS on 19 artefacts from 8 generators** since the recovery of 2026-09-26, with four negative controls and a 381-file snapshot that shows it still writes nothing (§29) |
| `python src/recovered/build_n3.py`, `python src/recovered/n3_fig5_whiskers.py`, `python src/recovered/make_fig_gapscaling.py` | recovered 2026-09-26 (§29). `build_n3.py` regenerates the ten `n3_*.dat` of Fig. S5 byte-for-byte from the two thm3 JSONs of `data/thm1iii_violation/` and, since 2026-09-28 (late), the deposited `data/cert_*.json`, and prints 621 / 606 / 469 / 363 / 2.005; `n3_fig5_whiskers.py`, run next, writes the whisker summary and the two whisker blocks of the fragment; `make_fig_gapscaling.py` regenerates the Fig. 6 fragment except one provenance-timestamp comment line. All write into `paper/figs/` (the caption and summary by-products go to `build/`); `check_figures.py` runs them in a scratch tree instead *(2026-09-26; nine generators, 21 artefacts since 2026-09-28)* |
| the six commands of the `figures` recipe, run in order | **exit 0** on all six, regenerated fragments byte-identical *(2026-09-18; `check_figures.py` re-runs the same six today)* |
| the two commands of the `data` recipe | **exit 0** — rebuilds `figs/sqw_edges.dat` and `data/charge_gap.json` from the deposited spectra, ~4 min of L=4…12 exact diagonalization *(2026-09-18)* |
| `python src/sqw_lanczos.py 12 --eta 0.18` | reproduces the **physics** of `data/sqw_L12.json` to `max|dS| = 8.2e-05`, `rel-L1 = 3.0e-06`, `|dE0| = 3.1e-15` — **not** byte for byte; the deposited file predates the seeded ARPACK start vector *(re-measured 2026-09-19, 425 s; §12, which withdraws the earlier bit-identity claim)*. Without `--eta 0.18` it **exits 2 and computes nothing**, instead of running ~700 s and overwriting the deposit at η = 0.20 |
| `python src/akw_lanczos.py <L>`, `spin_lanczos.py`, `spin_peak_deviation.py`, `free_fermion_support.py`, `sampled_akw.py`, `scaling_data.py`, `ladder_vs_chain.py`, `spectral_validation_max.py`, `gap_scaling.py` | run with numpy/scipy only; minutes at the deposited sizes *(2026-09-18)* |
| `python src/frontier/fsum.py` | **regenerates `data/c3_frontier/C3_grid.json` exactly** — deep-equal to the committed file, 204 grid rows and 16 frontier rows, in under a second from the 54 deposited point files *(2026-09-19)* |
| `python src/frontier/verdict.py 4 6 8`, `pubsum.py`, `calib_oos.py`, `null_ws4.py` | **regenerate `verdict_frontier.json`, `published_fraction/PUB_rows.json`, `null/CALIB_OOS.json` and `null/NULL_WS4.json` exactly** *(2026-09-19)* |
| `python src/frontier/frun.py`, `frun12.py`, `ftie.py`, `fvalid.py`, `run_ladder.py`, `z_fine.py`, `z_rank.py`, `z_e0.py`, `z_suppK.py`, `null_ws5/6/7.py` | all run; the recomputed points agree with the deposited ones to **3×10⁻¹⁶ … 6×10⁻⁹** on the headline columns (w_S, Λ̂/η, both bounds, rel-L₁). The residual is ARPACK's unseeded start vector. Per-file commands, arguments and costs: [`C3_FRONTIER.md`](C3_FRONTIER.md) §3 *(2026-09-19)* |

### Does **not** reach it, and why

| what | why not |
|---|---|
| the six figures that ship as PDF (Figs. 1, 4, 7, 8, S6, S7) | **their standalone sources are deposited since 2026-09-26** (`paper/figs/src/`), and each recompiles to a pixel-identical rendering of the committed PDF — but no script regenerates or guards them, and none of the three rasters they embed can be regenerated here: `akw_field_v2.png` is the output of the deposited `src/recolor_akw_v2.py`, whose viridis input is not deposited; `sqw_field.png` was recoloured by a step that is not deposited; and no deposited script writes `spinqw_field.png`. §5, §29. Their *numbers* are all in committed `.json`/`.dat`, and that is the claim `FIGURE_PROVENANCE.md` makes. |
| ~~**Fig. 3's generator**~~ | **repaired 2026-09-19 (§25).** It hard-coded two absolute Windows paths, one of them a scratch directory that is not deposited, and failed at the first `open` on any other machine. Both roots now come from `__file__`, and it regenerates the committed fragment in a tree that contains only `src/` and `data/`. |
| ~~**Fig. 6's generator**~~ (`fig:gapscaling`) | **recovered 2026-09-26** (`src/recovered/make_fig_gapscaling.py`, §29); what follows is the record. Its own header names `make_fig_gapscaling.py`, which was not in the deposit. The data (`gap_scaling.json`) and its compute script (`gap_scaling.py`) **are** here; the JSON→figure step is not. §27. |
| ~~**Fig. S5's generator and two of its inputs**~~ (`fig:thm1iii-violation`) | **recovered 2026-09-26** (`src/recovered/build_n3.py`, `data/thm1iii_violation/`, §29); `verify.py` had no check anchoring the ratios until 2026-09-28 (late), when its section 9.6c re-derived every plotted ratio from the certificate rows (§33). What follows is the record: `build_n3.py`, `fix_n3_caption.py`, `fix_n3_fig5_series.py`, `thm3_results.json` and `n3_fig5_summary.json` are all absent. The eight `n3_*.dat` tables it plots are committed and readable; nothing recomputes them. §27. |
| `src/recolor_akw_v2.py` | **neither of its two inputs is deposited.** It fails at the first `imread`, and nothing can be done about that from inside this repository. §13. |
| `src/gate1_ladder.py`, `src/rigor_floor.py` | they run, but their outputs were never deposited, so there is nothing to compare against. §13. |
| `src/frontier/fproto.py` — the C3 protocol gate | it needs `ckpt/L10_order.npy`, a by-product of `src/scaling_lanczos_mf.py`, which is not deposited: `ckpt/` is `.gitignore`d except for the L=14 checkpoint deposited on 2026-09-25. **It exits 3 with a message naming the file**, not a traceback. With `P2_CKPT` pointed at a tree that has it, the set overlap with the repository's own ranking is **1.000000** at every fraction *(measured 2026-09-19)*. [`C3_FRONTIER.md`](C3_FRONTIER.md) §5. |
| `src/frontier/run_L14.py` — the L=14 certificate point | it reuses `ckpt/L14_gs.npz` (94 MB) and `ckpt/L14_order.npy` (41 MB), both deposited since 2026-09-25, rather than spending three hours on a ground state. The **result** is deposited (`data/c3_frontier/published_fraction/L14_FR008.json` with its run log), and that file records that the ranking protocol was **not** re-run at L=14. Same exit-3 behaviour. |
| `src/frontier/run_ladder.py 4` | reproduces every field of `data/c3_frontier/ladder_L4.json` **except the ω grid**, which the library widened four minutes after that file was written. Registered as a `KNOWN_OPEN` defect, consequence measured: [`C3_FRONTIER.md`](C3_FRONTIER.md) §4. |
| everything requiring `pyscf` — `n19_suite.py`, `n19_spectral.py`, `n2_hero_data.py`, `molecular_noise_*.py` | `pyscf` is a hard dependency and is **not installed** in the environment these checks were run in. They are **untested here**; no claim is made about them beyond the fact that their output paths were repaired. The 269 molecular assertions of `verify.py` are therefore *transcription* checks against the deposited JSON, not recomputations. §11.2 item 3. |
| `make reproduce` | executes `notebooks/00_Reproduce_Everything.ipynb`, which needs `pyscf`. **Not verified in this pass.** The notebook is the v1–v2 pipeline, kept as a record: it reproduces the v1–v2 figure set, several of which v3 withdrew, and is not the reproduction path of the v3 manuscript. |
| `data/gflow.json` | `src/gflow_dequant.py` — the algorithm of `calculations/validate_realnoise.py` of the companion review (arXiv:2608.05314), unchanged, with every seed written to `data/gflow_runs/` and merged into `data/gflow.json` (`--merge`). Needs `pyscf`, `torch` and `qiskit-ibm-runtime` (FakeTorino readout rates): the Docker image `sqd-nb` plus `pip install qiskit-ibm-runtime`; twenty single-threaded runs (ten seeds at 500 and at 1000 shots, symmetry-pinned gauge, the companion's ladder configuration); run a few at a time, since twenty at once exhausted the memory of a 16-core laptop. |
| the two IBM Heron runs | the jobs are closed and a device run today has a different calibration. **The measured counts are not deposited** — the notebooks fetch them from the IBM Quantum service by job id — so neither the measurement nor the classical post-processing of the raw counts runs from the deposit alone. What is deposited, and can be checked offline, is the post-processed arrays that enter the figures, with backend, shots and job id. |
| `make paper` | needs a TeX distribution with `pgfplots ≥ 1.18` and REVTeX 4.2. Re-run on 2026-09-30 (its two `pdflatex` passes, MiKTeX): `paper/main.pdf` is that build (68 pp: main text pp. 1–28 with Appendix A, references pp. 29–33, Supplemental Material from p. 34), made on 2026-09-30 after the last source edit; 0 errors, 0 undefined references or citations, no rerun request, 0 overfull `\hbox` or `\vbox` and no float-too-large warning after the layout pass of that day, which changed no text (before it, two overfull `\vbox`, 24.3 pt, p. 48, and 5.0 pt, p. 55, and a float-too-large warning of 30.2 pt for Fig. 3; the p. 55 one came back with the final pass's Sec. S5 rewording and went again with a `\looseness=-1` on one SM paragraph, no word changed; `KNOWN_DISCREPANCIES.md` §34; late on 2026-09-28 one, 12.8 pt, p. 48; the build of that morning had two, 10.9 pt and 12.8 pt, and the one before the review pass 5.8 pt and 12.1 pt). It carries the edits made after the arXiv v3 submission (README.md lists them), so it is **not** the arXiv v3 PDF. |
| `make pra` | added 2026-09-28: runs `make paper`, then `python paper/pra_split.py --build`, which builds `paper/main_pra.pdf` (33 pp) and `paper/sm_pra.pdf` (38 pp) for Physical Review A and checks them against `paper/main.pdf` (exit 0 on 2026-09-30: 0 errors / undefined / multiply defined, no new overfull box or warning, the reference list in first-citation order with the SM entry naming exactly the 39 works cited only there, [104–142], main_pra pp. 1–33 pixel-identical to main.pdf pp. 1–33, sm_pra pp. 1–35 = main.pdf pp. 34–68 up to citation numbers — 20 pages pixel-identical, 15 differing in citation numbers only, line for line (19 and 16 before the layout pass of 2026-09-30); late on 2026-09-28 it was 38 works, 18 and 17 pages). Needs `pdftotext`, `pdftoppm`, Pillow and `xr-hyper` ≥ 2023. |
| the CI workflow | runs `make verify` and `make check-figures` on every push since 2026-09-19 (§24, closed); `make check-coherence` was added as a third job on 2026-09-26 and runs from the next push on. Listed here only because this row used to say it had never executed |

**The Data Availability Statement drafted from these facts is in
[`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md).** It must claim nothing this table does not support.
