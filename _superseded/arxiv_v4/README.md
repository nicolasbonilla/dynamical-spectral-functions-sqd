# `_superseded/arxiv_v4/` — files of deposit 3.4.1 (arXiv:2608.16436v4) that the JCP text does not use

Moved or copied here on 2026-11-18, when deposit 3.5.0 put the manuscript prepared for The Journal of Chemical
Physics in `paper/`. The arXiv-v4 sources themselves are not kept here: they are the files of `paper/` at tag
`v3.4.1` (Zenodo 10.5281/zenodo.23088494). One line per file:

- `paper/main.pdf` — the arXiv-v4 build of 3.4.1 (68 pages, main text and Supplemental Material in one document); `docs/img/make_thumbnails.py` crops the README thumbnails from it.
- `paper/main_pra.tex` — the Physical Review A driver of 3.4.1 (main text and references); it runs the `main.tex` of 3.4.1, so it builds only at tag `v3.4.1`.
- `paper/sm_pra.tex` — the Physical Review A driver of the Supplemental Material of 3.4.1; same remark.
- `paper/main_pra.pdf` — its build (33 pages), as deposited in 3.4.1.
- `paper/sm_pra.pdf` — its build (38 pages), as deposited in 3.4.1.
- `paper/pra_split.py` — the script that built and checked the two PRA PDFs against `paper/main.pdf` of 3.4.1 (`make pra` of 3.4.1).
- `paper/figs/n3_thmB_frac.dat`, `paper/figs/n3_thmB_eta.dat` — the blue series of Fig. S5 of 3.4.1 (error over the first-order certificate); `src/make_fig_thmB_data.py --old` rebuilds both byte for byte, as its control. The JCP text plots the second-order certificate (`paper/figs/` of 3.5.0).
