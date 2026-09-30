# `data/thm1iii_violation/` — the inputs of Fig. S5 (`fig:thm1iii-violation`)

Deposited 2026-09-26. These files were written on 2026-09-18 in the author's session scratch
directory and read from there by `build_n3.py` (now `src/recovered/build_n3.py`) to build the ten
`paper/figs/n3_*.dat` tables of Fig. S5. They were recovered unchanged; the untouched copies, with
SHA-256 checksums, are kept outside this repository.

| file | rows | written by | what it is |
|---|---|---|---|
| `thm3_results.json` | 239 | `src/recovered/thm3_run.py` | the original sweep: 3 Hubbard L=6 Hamiltonians × 8 subspace families |
| `thm3_L8b.json` | 4 | `src/recovered/thm3_L8b.py` | a **site-seeded** L=8 run — not the paper's Fig. 3 configuration, which seeds with c†_{k↑}\|GS⟩ |
| `cert_earlier_run/cert_stress.json` | 162 | `src/leakage_certificate_suite.py`, earlier run | the stress families |
| `cert_earlier_run/cert_akw.json` | 112 | same | the A(k,ω) certificate scan; its 64 L=8 rows are the paper's Fig. 3 setting |
| `cert_earlier_run/cert_teqsci.json` | 104 | same | the TE-QSCI scan |
| `n3_fig5_summary.json` | — | `src/recovered/n3_fig5_whiskers.py` (since 2026-09-28; `fix_n3_fig5_series.py` before) | min / median / max of the 64 Fig.-3 rows at each abscissa, the whiskers drawn in Fig. S5 |

**`cert_earlier_run/` is not the same as the deposited `data/cert_*.json`.** Both are runs of the same
three scans; the plotted pool was built from the earlier one, and the two differ on 28 XXZ and t–V
stress subspaces. `src/recovered/build_n3.py` on these files gives the numbers the paper prints —
621 pooled, 606 live, the weight-only bound violated on 468, on all 363 with `w_S ≥ 0.99`, mildest
factor 2.10 — and regenerates all ten `n3_*.dat` byte-for-byte. With `data/cert_*.json` in their
place it gives 469 and 2.005, which is what the paper's Sec. S9 says. Both were measured on
2026-09-26.

**2026-09-28 (late): the pool is now built from the deposited `data/cert_*.json`.** `build_n3.py`
reads them instead of `cert_earlier_run/`, Fig. S5 and every count quoted from it were regenerated
(469 violations of 606, all 363 with `w_S ≥ 0.99`, mildest factor 2.005, survival on 137 of 243 below
0.99), and the paper prints one set of numbers. `cert_earlier_run/` stays deposited as the record of the
run the figure was first built from and is no longer read by any deposited script.
`n3_fig5_summary.json` is now written by `src/recovered/n3_fig5_whiskers.py` from the two
`n3_fig5_*.dat` tables (`fix_n3_fig5_series.py` is the record of the one-shot patch), and
`src/verify.py` section 9.6c re-classifies the pool from the deposited sources and checks every plotted
ratio and printed count against it.
