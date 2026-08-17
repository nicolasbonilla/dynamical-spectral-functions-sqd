# FIGURE DATA PROVENANCE — single source of truth

**Rule of order:** every figure's numbers come from ONE authoritative data file listed below.
No hand-typed values may drift from these files. Hardware figures carry the exact **job id**.
Before editing any figure, check this table. Before trusting any `.json`, check it is **not** in the
"STALE / DO NOT USE" list at the bottom.

Paths: figure sources live in `rebuild/figs/`; generating scripts and data live one level up in
`02_Paper_Amortizacion/`.

| # | Body ref | Figure source (`rebuild/figs/`) | Generator | AUTHORITATIVE data | Notes |
|---|----------|--------------------------------|-----------|--------------------|-------|
| 1 | fig:decoupling | fig_decoupling_native.tex | make_decoupling_native.py | n19_suite.json + cost_vs_ent.json | panel (b) shows 3 of 6 Hubbard sweeps (representative); panel (c) free-fermion \|S\| hard-typed [35,336,3496,37361] |
| 2 | fig:master | fig_resource_master_native.tex | make_resource_master_fig.py | **resource_master.json** | chain L=6 F1 CORRECTED 9.6047→**9.3851** (ED-certified); ρ(χ,\|S\|)=0.72 pooled, ρ(F1,\|S\|)=0.60 |
| 3 | fig:method | fig_method_native_frag.tex | (native, inline) | figs/aw_method.dat + figs/conv_method.dat + method_max.json | panel (b) sampled floor = **~0.019** (method_max.json K=12=0.0185), NOT 0.025 |
| 4 | fig:lattice | fig_akw_native.pdf | akw_field*/recolor_akw_v2.py | **akw_lanczos_L12.json** | exact A(k,ω), L=12, U/t=8, η=0.18t; peak split 4.99t vs exact gap 4.97t (η broadening) |
| 5 | fig:akwsampled | fig_akw_sampled_native.tex | (script `_mk_akw_sampled`) | **sampled_akw_L8.json** | REAL sampled vs exact, L=8, 85% sector, mean rel-L1 0.0022 |
| 6 | fig:sqw | fig_sqw.pdf | sqw_lanczos.py / sqw_field.py | **sqw_L12.json** + figs/sqw_edges.dat | S(q,ω); reflection residual ~1e-4 (NOT machine precision) |
| 6b | fig:spin | fig_spinqw.pdf | `fig_spinqw.tex` (hand standalone) + spinqw_field.png | **spinqw_L12.json** + figs/spinqw_edges.dat | 2026-08-17: added cyan "extracted peak" overlay from the `peak` column of spinqw_edges.dat |
| 7 | fig:spin | fig_spinqw.pdf | spin_lanczos.py | **spinqw_L12.json** + figs/spinqw_edges.dat | S^zz(q,ω); reflection residual 7e-13 (machine precision, genuine) |
| 8 | fig:gallery | fig_molecular_gallery.pdf | make_gallery_native.py | **n19_suite.json** (PyMOL orbitals) | 19 molecules; names in dark ink, magic bar carries regime colour |
| 9 | fig:molsuite | fig_molecular_suite.pdf | make_magic_suite_fig.py | **n19_suite.json** | F1/F1max eq vs diss, 19 molecules |
| 10 | fig:hero | fig_hero.pdf | make_hero_fig.py | **n2_hero.json** + figs/hero_aw.dat + figs/hero_res.dat | N2 dissociation; uses italic N_u (body uses roman N_u) |
| 11 | fig:bench | fig_benchmark.pdf | make_bench_fig.py | **headtohead_ms.json** + figs/bench_U*.dat | U/t=4,8,12; 0.042/0.729/0.431 at U/t=12,\|S\|=200 |
| 12 | fig:scaling | fig_scaling2_native.tex | make_scaling_fig.py (⚠ panel c added by hand to the .tex) | **scaling_data.json** | F1 6.6→22.5, frac 0.92→0.08; 2026-08-17 added panel (c) |S|=frac·Np1_sector & D=Np1_sector, both exponential. `figs/scaling.dat` DELETED (orphan) |
| 13 | fig:ladder | fig_ladder_native.tex | ladder_vs_chain.py | **ladder_vs_chain.json** | ⚠ json chain-L6 F1 still 9.6047 (stale) — figure/caption use corrected 9.39; overwrite json to match |
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
