# FIGURE DATA PROVENANCE — single source of truth

**Rule of order:** every figure's numbers come from ONE authoritative data file listed below.
No hand-typed values may drift from these files. Hardware figures carry the exact **job id**.
Before editing any figure, check this table. Before trusting any `.json`, check it is **not** in the
"STALE / DO NOT USE" list at the bottom.

Paths: figure sources live in `release/paper/figs/`; generating scripts and data live one level up in
`release/`. **The Aug-17 `rebuild/` tree is a SUPERSEDED scratch copy carrying pre-2026-09-05
values (it still plots the retracted `(9.385,147)` and records the reversed provenance note) --
DO NOT USE IT.** (The `(9.385,147)` there is wrong *in the resource-master figure*, which is the
OPEN chain and must read `(9.6047,147)`; the ring value 9.3851 is correct where it belongs, in
`fig_scaling2_native.tex` — see H3 in `KNOWN_DISCREPANCIES.md`.)

**Known gaps between what a script does and what the deposit contains are listed in
[`KNOWN_DISCREPANCIES.md`](KNOWN_DISCREPANCIES.md). Read it before re-running any compute script.**

## The two exceptions to "every plotted value is deposited"

The Data Availability Statement ([`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md)) claims that each
figure's plotted values are deposited as a JSON or plain-text table. That claim is true of eighteen of
the twenty figures. It is **not** true of these two, which are therefore named in the statement itself
rather than left inside a table row:

| figure | what is not in `data/` |
|---|---|
| **Fig. 1, panel (c)** | the four free-fermion support sizes `35, 336, 3496, 37361` are constants typed into `src/make_decoupling_native.py:34` and written from there into the fragment. **No file in `data/` holds them and no deposited script recomputes them.** The same four numbers are printed in `paper/resource.tex`. |
| **Fig. 14** | a circuit schematic. It plots no data at all, so there is nothing to deposit. |

Everything else in this table resolves to a file.

| # | Body ref | Figure source (`rebuild/figs/`) | Generator | AUTHORITATIVE data | Notes |
|---|----------|--------------------------------|-----------|--------------------|-------|
| 1 | fig:decoupling | fig_decoupling_native.tex | make_decoupling_native.py | n19_suite.json + cost_vs_ent.json | panel (b) shows 3 of 6 Hubbard sweeps (representative); panel (c) free-fermion \|S\| hard-typed [35,336,3496,37361] |
| 2 | fig:master | fig_resource_master_native.tex | make_resource_master_fig.py | **resource_master.json** | chain L=6 F1 = **9.6047** (OPEN chain, ED). 2026-08-17 had overwritten it with the PERIODIC-ring 9.3851 from scaling_data.json — a different system; reverted 2026-09-05, rank-preserving so ρ(χ,\|S\|)=0.72 pooled and ρ(F1,\|S\|)=0.60 are unchanged |
| 3 | fig:method | fig_method_native_frag.tex | make_method_fig_max.py | figs/aw_method.dat + **method_max.json** | panel (b) sampled floor = **~0.019** (method_max.json K=12=0.0185), NOT 0.025. `figs/conv_method.dat` **DELETED 2026-09-18**: stale single-seed orphan (its K=12 row read 0.024984 against the correct 0.018550 of method_max.json); it had no writer and no reader anywhere in the deposit |
| 4 | fig:lattice | fig_akw_native.pdf | akw_field*/recolor_akw_v2.py (path repaired 2026-09-18; **its two inputs, `akw_field_v2.png` and `fig_akw_v2_L12.tex`, are NOT in the deposit** — see KNOWN_DISCREPANCIES §5) | **akw_lanczos_L12.json** | exact A(k,ω), L=12, U/t=8, η=0.18t; peak split 4.99t vs exact gap 4.97t (η broadening) |
| 5 | fig:akwsampled | fig_akw_sampled_native.tex | make_akw_sampled_fig.py | **sampled_akw_L8.json** | L=8, 85% of the sector, mean rel-L1 0.0022. **CORRECTED 2026-09-18 — this row said "REAL sampled" and that was wrong.** `sampled_akw.py:54-57` accumulates the *exact* Born weights `wc` of the time-evolved seed over `K+1` steps and takes `np.argsort(wc)[::-1][:round(FR*nS)]`. There is no random number generator, no shot count and no seed anywhere in the file: this is the **infinite-shot limit of the declared protocol**, i.e. the subspace the device would select given unlimited measurements — not an oracle, but not a finite-shot run either. **No shot budget is quoted because none was ever drawn.** At matched subspace size the finite-shot number is worse (measured: 1.45x at L=6, 1.83x at L=8). Six other scripts in this deposit *do* sample with an RNG, seeds and +/-sigma bands (`spectral_validation_max.py`, `headtohead_ms.py`, `noise_spectral.py`, `molecular_noise_energy.py`, `molecular_noise_sweep.py`); this one does not. |
| 6 | fig:sqw | fig_sqw.pdf | sqw_lanczos.py / sqw_field.py; edges by **make_sqw_edges.py** | **sqw_L12.json** + figs/sqw_edges.dat | S(q,ω); reflection residual ~1e-4 (NOT machine precision). Until 2026-09-18 `figs/sqw_edges.dat` had **no generator anywhere in the repository**; `make_sqw_edges.py` now rebuilds it from `data/sqw_L12.json` |
| 6b | fig:spin | fig_spinqw.pdf | `fig_spinqw.tex` (hand standalone) + spinqw_field.png | **spinqw_L12.json** + figs/spinqw_edges.dat | 2026-08-17: added cyan "extracted peak" overlay from the `peak` column of spinqw_edges.dat |
| 7 | fig:spin | fig_spinqw.pdf | spin_lanczos.py | **spinqw_L12.json** + figs/spinqw_edges.dat | S^zz(q,ω); reflection residual 7e-13 (machine precision, genuine) |
| 8 | fig:gallery | fig_molecular_gallery.pdf | make_gallery_native.py | **n19_suite.json** (PyMOL orbitals) | 19 molecules; names in dark ink, magic bar carries regime colour |
| 9 | fig:molsuite | fig_molecular_suite.pdf | make_magic_suite_fig.py | **n19_suite.json** | F1/F1max eq vs diss, 19 molecules |
| 10 | fig:hero | fig_hero.pdf | make_hero_fig.py | **n2_hero.json** + figs/hero_aw.dat + figs/hero_res.dat | N2 dissociation; uses italic N_u (body uses roman N_u) |
| 11 | fig:bench | fig_benchmark.pdf | make_bench_fig.py | **headtohead_ms.json** + figs/bench_U*.dat | U/t=4,8,12; 0.042/0.729/0.431 at U/t=12,\|S\|=200 |
| 12 | fig:scaling | fig_scaling2_native.tex | ⚠ **hand-maintained** — `make_scaling_fig.py` emits a DIFFERENT two-panel figure and would DESTROY panel (c) if run; see [KNOWN_DISCREPANCIES.md](KNOWN_DISCREPANCIES.md) | **scaling_data.json** | F1 6.6→22.5, frac 0.92→0.08; 2026-08-17 added panel (c) |S|=frac·Np1_sector & D=Np1_sector, both exponential. |S| = round(frac·Np1_sector) = 22, 246, 2180, **19183**, 124407, **824503** (the last two were printed as 19184 / 824504 until 2026-09-18). `figs/scaling.dat` DELETED (orphan). **Same caveat as row 5:** the fractions come from `scaling_lanczos.py:66` / `scaling_lanczos_mf.py:101` / `scaling_data.py:56`, which rank determinants by the **exact** time-averaged Born weight and cut at a threshold — the infinite-shot limit, not a shot-limited run. |
| 13 | fig:ladder | fig_ladder_native.tex | ladder_vs_chain.py | **ladder_vs_chain.json** | chain rows are OBC; 2026-08-17 overwrote their F1 with the PERIODIC-ring values from scaling_data.json (L6 9.3851, L8 12.8138, L10 16.0183). Reverted 2026-09-05 to the OBC values 9.6047 / 12.8163 / 16.0306; caption now reads 9.60/9.69. The figure plots only χ and fractions, which were never affected |
| 14 | fig:circ | fig_circuit.pdf | fig_circuit_native.tex (quantikz) | — (schematic, no data) | interaction column drawn as ctrl dots, labelled R_zz |
| 15 | fig:heron | fig_hardware_hero_frag.tex | make_hardware_hero.py | **hw_lucj_n2_result.json — job `da125f2ein7c73bcsqs0`** | ⭐ AUTHORITATIVE HW RUN. Panel (b) = **8-seed recovery MEAN per step + shaded ±std envelope** (`dE_hist_mean`/`dE_hist_std` = [31.1±4.0, 6.1±0.8, 2.1±0.3, 1.0±0.2, **0.59±0.14**]); capped whisker at converged pt; noiseless 29.5 |
| 16 | fig:noise | fig_noise_score.pdf | make_noise_fig.py | **noise_spectral.json** + figs/noise.dat | naive vs S-CoRe; "naive" curve colour clashes with orange=learned elsewhere |
| 17 | fig:noiserec | fig_noise_recovery_native.tex | make_noise_recovery_native.py | **molecular_noise.json** + molecular_noise_sweep.json | panel (a) 5 mols (8 seeds), panel (b) 19 mols at ε=2% (5 seeds) — different ensembles |
| 18 | fig:gflow | fig_gflownet.pdf | make_gflow_fig.py | **gflow.json** | bars 43.2/39.2/27.0/26.0/22.9 |
| 19 | fig:amort | fig_amort_native.tex | amortized_recovery.py | **amortized_recovery.json** | classical vs amortized LOO |
| 20 | fig:witness | fig_witness_native.tex | (native) | **apsg_witness.json** | F1=F2=4K, \|S\|=2^K, χ=2; uses roman N_u (correct) |

## ⭐ HARDWARE — the single authoritative run (do not confuse)

**N₂ / ibm_marrakesh energy = job `da125f2ein7c73bcsqs0`** (8-seed recovery sweep, "Run all").
File: `hw_lucj_n2_result.json`. The **full per-seed × per-iteration trajectories** are stored in
`e_hw_seed_hist` (8×5 energies) and `dE_hist_seeds` (8×5 mHa), with `dE_hist_mean`/`dE_hist_std` the
per-step 8-seed statistics — reconstructed from the run log so the per-point error bars are never lost
again. Converged dE = **0.59 ± 0.14 mHa** (best seed 0.40); noiseless sim 29.5 mHa.
Confirmed authoritative by Nicolás (2026-08-17); per-seed trajectories restored 2026-08-17.

**A(ω) / ibm_fez** = the `ibm_fez` L=6 Hubbard run (|S|=300/300, exact by coverage). File: `heron_spectral.json` / figs/heron*.dat.

## ⛔ STALE / DO NOT USE (archive or ignore)

- **`notebooks/HW_LUCJ_N2_Heron_READY.ipynb` SAVED CELL OUTPUTS** = an OLDER run, job `d9qe731dsedc73af67d0` (single seed, final 1.2 mHa, noiseless 28.5). These are NOT the paper's numbers. The authoritative run is `da125…` in `hw_lucj_n2_result.json`. Re-run "Run all" and re-save the notebook to sync, or ignore the saved outputs.
- `_archive/figuras/hw_lucj_n2_result.json` — old (Aug 6), superseded.
- `scaling_data_dense_backup.json`, `sqw_L12_eta020_backup.json` — backups, not the plotted data.
- `figs/scaling.dat` — DELETED (was an orphan stale copy of scaling_data.json).
- `figs/conv_method.dat` — **DELETED 2026-09-18.** Stale single-seed convergence table (13 rows) whose
  K=12 entry read `240  0.024984` while the authoritative `data/method_max.json` (8 seeds) gives
  `|S|=242.5  rel-L1=0.018550`. No script in the deposit wrote it and no fragment read it; the paper
  prints the correct value. Superseded generator of the same family: `_superseded/_archive/`.
- `figs/resource_axis.dat` — **DELETED 2026-09-18.** Orphan of a figure that never entered the paper
  (header `t Fn Sn nn F S n`, 16 rows). No writer, no reader, and no committed `.json` carries that
  column set; the only surviving generator is `_superseded/_archive/resource_axis_v2.py`.
