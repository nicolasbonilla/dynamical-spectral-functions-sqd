# `paper/figs/src/` — standalone sources of the six figures that enter the paper as PDF

Six figures of the manuscript are included as compiled PDFs (`\includegraphics`), not as native
fragments. Until 2026-09-26 their standalone LaTeX sources were not in this deposit
(`docs/KNOWN_DISCREPANCIES.md` §5). They were recovered that day from the author's session scratch
directory of 2026-09-18/19, where each sat next to a byte-identical copy of the committed PDF, and are
deposited here unchanged except for line endings (CRLF → LF). Nothing in the manuscript reads this
folder; the PDFs in `paper/figs/` remain what the paper includes.

| source | builds | printed as | reads |
|---|---|---|---|
| `fig_circuit_native_FIXED.tex` | `fig_circuit.pdf` | Fig. 1, `fig:circ` | — (schematic, `quantikz`) |
| `fig_akw_fixed.tex` | `fig_akw_native.pdf` | Fig. 4, `fig:lattice` | `akw_field_v2.png` |
| `fig_sqw_fixed.tex` | `fig_sqw.pdf` | Fig. 7, `fig:sqw` | `sqw_field.png`, `sqw_edges.dat` |
| `fig_spinqw_fixed.tex` | `fig_spinqw.pdf` | Fig. 8, `fig:spin` | `spinqw_field.png` |
| `fig_hero_STANDALONE.tex` | `fig_hero.pdf` | Fig. S6, `fig:hero` | `paper/figs/hero_aw.dat`, `paper/figs/hero_res.dat` |
| `fig_noise_score_STANDALONE.tex` | `fig_noise_score.pdf` | Fig. S7, `fig:noise` | `paper/figs/noise.dat` |

**Measured 2026-09-26** (MiKTeX `pdflatex`, one pass each): all six compile with 0 errors, and each
result renders **pixel-identical** to the committed PDF at 216 dpi (PyMuPDF), with identical extracted
text and page size. For `fig:sqw` the source that does so is the last one written, at 00:28 on
2026-09-19; an earlier version of the same file differs from the committed PDF in one tick label
(`≈8` instead of `8.8`), and is kept outside this repository with the rest of the recovered material.
To repeat it:

```bash
mkdir -p /tmp/figsrc/paper/figs
cp paper/figs/src/* paper/figs/*.dat /tmp/figsrc/ && cp paper/figs/*.dat /tmp/figsrc/paper/figs/
cd /tmp/figsrc && for f in *.tex; do pdflatex -interaction=nonstopmode "$f"; done
```

**The three rasters are deposited as rendered, and are not regenerated here.** `akw_field_v2.png` is
the hot-coloured output of `src/recolor_akw_v2.py`, whose own input (the viridis raster) is still not
deposited (§13). `src/sqw_field.py` renders the same `S(q,ω)` field from `data/sqw_L12.json` at the
same size, but in viridis; the hot-coloured `sqw_field.png` embedded here was recoloured by a step that
is not deposited. No deposited script is known to write `spinqw_field.png`. What every raster plots is
in the committed JSON (`akw_lanczos_L12.json`, `sqw_L12.json`, `spinqw_L12.json`).
