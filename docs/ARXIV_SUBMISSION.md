# arXiv submission — the preprint is ready

> **Scope of this document (added 2026-09-18).** Everything below describes the bundle that was posted
> as **arXiv v2** (39 pp). It is *not* a description of the current working tree. Nine files in
> `paper/` have moved ahead of `paper/arxiv-submission.tar.gz` since v2 was posted and two more have
> been deleted from the tree altogether, and `paper/main.pdf` now builds to **40 pp**. Run
> `python src/check_tarball.py` for the current comparison rather than trusting this sentence — it was
> wrong once already. The bundle still carries the two fabricated rows of `figs/sqw_edges.dat`, so it
> must be rebuilt, and this page re-verified, before any new submission. See
> [`KNOWN_DISCREPANCIES.md`](KNOWN_DISCREPANCIES.md) §7.

The upload-ready source bundle is [`../paper/arxiv-submission.tar.gz`](../paper/arxiv-submission.tar.gz)
(~4.5 MB, well under arXiv's 50 MB limit). It has been verified to compile **from a clean extraction**
with `pdflatex` alone — the exact process arXiv's AutoTeX runs.

## Verification (what was checked)

| Check | Result |
|---|---|
| Clean-room `pdflatex main.tex` ×2 from the extracted tarball | ✅ 39 pages, exit 0 |
| Undefined references / citations | ✅ 0 |
| LaTeX / package errors | ✅ 0 |
| Sole `\documentclass` (unambiguous main file) | ✅ `main.tex` only (`article`, 11pt) |
| External `.bib` / `.bbl` needed | ✅ none — bibliography is an inline `thebibliography` |
| Absolute paths in `\input` / `\includegraphics` | ✅ none |
| All `\input` (21) and `\includegraphics` (10) targets present | ✅ all resolve |
| Fonts | ✅ Latin Modern + AMS only (embedded; standard on arXiv) |
| Embedded PDF metadata (title/author/keywords) | ✅ set via `\hypersetup` |

## What the bundle contains / omits

- **Contains:** the eleven section `*.tex` files, `figs/` with the native pgfplots `\input` fragments,
  the vector figure PDFs actually referenced, and the plotted `.dat` tables.
- **Omits (deliberately):** `main.pdf` (arXiv rebuilds it), `figs/_method_standalone.tex` (a standalone
  helper with its own `\documentclass` — would confuse AutoTeX about which file is primary) and its
  unused `_method_standalone.pdf`.

## Rebuild the bundle from source

```bash
cd paper
tar czf arxiv-submission.tar.gz \
  --exclude='main.pdf' --exclude='arxiv-submission.tar.gz' \
  --exclude='figs/_method_standalone.tex' --exclude='figs/_method_standalone.pdf' \
  *.tex figs/
```

## Submission form metadata

| Field | Value |
|---|---|
| **Primary category** | `quant-ph` (Quantum Physics) |
| **Cross-list** | `cond-mat.str-el` (strongly-correlated electrons; Hubbard spectra, Mott gap) |
| **Title** | Dynamical spectral functions from bitstring-sampled quantum subspaces: entanglement, not one-body magic, tracks the sampling cost |
| **Authors** | Nicolás Bonilla Vargas |
| **Abstract** | copy from the compiled PDF (`paper/main.pdf`), or from the `abstract` environment in `main.tex` |
| **License** | CC BY 4.0 (matches this repository's paper license) |
| **Comments** | e.g. "39 pages, 20 figures. Code and full reproduction: <this repo URL>" |

## Steps

1. Sign in at <https://arxiv.org> → **Submit**.
2. Upload `paper/arxiv-submission.tar.gz`. AutoTeX detects `main.tex` and compiles it.
3. Review arXiv's generated PDF against `paper/main.pdf` (should be identical, 39 pp).
4. Enter the metadata above; select the CC BY 4.0 license.
5. Add the repository URL in **Comments** for reproducibility.
6. Submit; note the assigned `arXiv:2608.16436` and update the badge/`CITATION.cff` in this repo.
