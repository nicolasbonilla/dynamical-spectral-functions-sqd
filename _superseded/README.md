# `_superseded/` — material the current manuscript no longer uses, kept for the trail

Nothing in `paper/` reads anything in this directory. It is kept so that every retraction and every
replaced file stays legible; see [`docs/FILE_INDEX.md`](../docs/FILE_INDEX.md) for the file-by-file list.

| subfolder | what it holds |
|---|---|
| `v2_paper/` | the v2-era manuscript as it stood in `paper/` on 2026-09-19 (edited after the 2026-08-22 posting; the posted v2 source is `paper/arxiv-submission.tar.gz`). Measured 2026-09-26: six of its ten source files (`main.tex`, `resource.tex`, `results.tex`, `discussion.tex`, `conclusion.tex`, `theorem_b2.tex`) differ from the posted bundle, and its `main.pdf` is a 40-page build, not the posted v2 PDF. |
| `v2_figs/` | `fig_akw_sampled_native.tex` and its generator `make_akw_sampled_fig.py`, retired on 2026-09-19 because the fragment labelled an infinite-shot curve "sampled (85% sector)". |
| `v3_retired/` | the bodies and data of the five figures v3 withdrew (`fig:bench`, `fig:gflow`, `fig:gallery`, `fig:molsuite`; the fifth, `fig_amort_native.tex`, stays in `paper/figs/` because `src/verify.py` names it) and earlier standalone builds the native fragments replaced. |
| `v1_readme_img/` | the ten README thumbnails cropped from the v1 PDF on 2026-08-17, superseded on 2026-09-26; its own `README.md` says what each showed. |
| `arxiv_v4/` | added 2026-11-18 (deposit 3.5.0): the files of deposit 3.4.1, the arXiv-v4 text, that the manuscript prepared for JCP does not use — the arXiv-v4 PDF (`paper/main.pdf`), its Physical Review A split (`main_pra.*`, `sm_pra.*`, `pra_split.py`) and the first-order tables of the blue series of Fig. S5 (`paper/figs/n3_thmB_*.dat`); its own `README.md` gives one line per file. The arXiv-v4 sources are at tag `v3.4.1`. |
