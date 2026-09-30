# `src/recovered/` — two figure generators and their compute, recovered on 2026-09-26

Until 2026-09-26 the manuscript printed two figures whose own source files name a generator that was
not in this deposit (`docs/KNOWN_DISCREPANCIES.md` §27): Fig. 6 (`fig:gapscaling`, built by
`make_fig_gapscaling.py`) and Fig. S5 (`fig:thm1iii-violation`, built by `build_n3.py`). Both scripts,
and the scripts behind two of the inputs of the second, were written on 2026-09-18 in the author's
session scratch directory and never committed. They were found there on 2026-09-26, copied out of it
(the untouched copies, with SHA-256 checksums, are kept outside this repository), and deposited here.
The five record files are unchanged: four byte-identical, and `fix_n3_caption.py` identical except
for its line endings (CRLF → LF). The two generators were ported, as the table says.

They live in a subdirectory, not in `src/`, on purpose: they are a recovery, not part of the pipeline
the manuscript was submitted with, and `src/verify.py` counts the files at the top level of `src/`.

| file | status | what it does |
|---|---|---|
| `make_fig_gapscaling.py` | **ported, runs, guarded** | writes `paper/figs/fig_gapscaling_native.tex` (Fig. 6) from `data/gap_scaling.json`, `data/spinqw_L12.json` and `data/akw_lanczos_L12.json`. Paths now come from `__file__`; the caption and a facts file go to `build/`, because the committed caption was edited by hand after generation and is the authoritative one. From the deposited data the fragment matches the committed one except for one comment line, line 5, which records the timestamp of the `gap_scaling.json` it was built from (the numerical content of the two JSON files is identical). |
| `build_n3.py` | **ported, runs, guarded** | writes the ten `paper/figs/n3_*.dat` tables of Fig. S5 from `data/thm1iii_violation/` and prints the counts the paper quotes. Paths now come from `__file__`, and two search keys of its legend self-check were updated to the legend wording the fragment prints today. Byte-identical on all ten tables. |
| `fix_n3_fig5_series.py`, `fix_n3_caption.py` | record only — **do not run** | one-shot patchers that, on 2026-09-18, replaced the 64-diamond overlay of the Fig. S5 fragment by min/median/max whiskers (writing `n3_fig5_summary.json`) and traced two caption numbers. They edit the fragment and the caption in place, in the directory they sit in; the committed fragment and caption have been edited by hand since, so they are kept as the record of how those literals were made, with their original relative paths, not as tools. |
| `thm3_run.py`, `thm3_core.py` | record | the original Theorem 1(iii) sweep: 3 Hubbard L=6 Hamiltonians × 8 subspace families, 239 rows → `data/thm1iii_violation/thm3_results.json`. They write next to themselves, and use `np.trapz`, which numpy 2.0 removed: they run under numpy 1.x, as they did on 2026-09-18. |
| `thm3_L8b.py` | record (and read by `build_n3.py`) | a **site-seeded** L=8 run, 4 rows → `data/thm1iii_violation/thm3_L8b.json`. Not the paper's Fig. 3 configuration, which seeds with c†_{k↑}|GS⟩; its docstring says so. `build_n3.py` reads `L`, `U` and `η` out of this file's source. Same numpy 1.x note. |

Measured on 2026-09-26, in a scratch copy of the repository and again through `src/check_figures.py`:

- `python src/recovered/build_n3.py` regenerates all ten `n3_*.dat` byte-for-byte and prints 621
  pooled / 606 live / 468 violations of the weight-only bound / 363 of 363 at `w_S ≥ 0.99`, mildest
  factor 2.10 — the numbers of Secs. III D and S9 and of the Fig. S5 caption. Those come from an
  **earlier run** of the three certificate scans (`data/thm1iii_violation/cert_earlier_run/`). With the
  deposited `data/cert_*.json` in their place it gives **469** and **2.005** instead, and survival on
  137 of 243 below `w_S = 0.99` — the discrepancy the paper states in Sec. S9.
  *2026-09-28 (late): `build_n3.py` now reads the deposited `data/cert_*.json` (a one-line change of
  its `CERT` path, recorded in its docstring), Fig. S5 and every count quoted from it were regenerated
  (621 / 606 / 469 / 363 of 363 / 2.005 / 137 of 243), the paper prints one set of numbers, and
  `cert_earlier_run/` stays deposited as a record and is no longer read. The whiskers of the 64 Fig.-3
  rows are written by the new `n3_fig5_whiskers.py`, run after `build_n3.py`; `check_figures.py` runs
  both and compares the ten tables, the fragment and `n3_fig5_summary.json`.*
- `python src/recovered/make_fig_gapscaling.py` exits 0, all of its internal guards pass, and the
  fragment differs from the committed one in line 5 only.

Still open, and not closed by this deposit: no check in `src/verify.py` anchors the plotted ratios of
Fig. S5 to the certificate rows. `verify.py` keeps printing that as its declared coverage gap.
*Closed 2026-09-28 (late): `verify.py` section 9.6c re-classifies the pool from the five JSON sources
independently of `build_n3.py` and requires every plotted point to be the ratio of a classified row.*
