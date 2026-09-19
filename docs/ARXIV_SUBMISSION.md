# arXiv submission — the **v3** source package

> # ⛔ DO NOT UPLOAD `paper/arxiv-submission.tar.gz`
>
> That file is the **frozen arXiv v2 package**, kept deliberately as a historical record. Measured,
> not remembered (`python src/check_tarball.py`): **47 files**, and they are the wrong ones. It
> carries `introduction.tex`, `method.tex`, `results.tex`, `discussion.tex`, `conclusion.tex`,
> `hardware.tex`, `resource.tex`, `theorem_b2.tex` — none of which exist in this manuscript any more —
> plus **five withdrawn figures** (`fig_gflownet.pdf`, `fig_benchmark.pdf`, `fig_molecular_gallery.pdf`,
> `fig_molecular_suite.pdf`, `fig_amort_native.tex`), the **two fabricated rows** of
> `figs/sqw_edges.dat`, and the superseded `19184` / `824504`. Uploading it would re-post every one of
> them. See [`KNOWN_DISCREPANCIES.md`](KNOWN_DISCREPANCIES.md) §7.
>
> **The package to upload is built by a script and lands in `build/arxiv-v3.tar.gz`.**

```bash
python src/build_arxiv_bundle.py          # build + audit + clean-room compile + PDF equivalence
python src/build_arxiv_bundle.py --list   # just print the derived file list
python src/build_arxiv_bundle.py --selftest   # the five negative controls
```

`build/` is gitignored on purpose: the package is a **build product**, one command away at any time,
never a deposited artefact that can go stale behind the sources the way the v2 tarball did.

---

## 1. How the file list is derived — and why it is not written by hand

`src/build_arxiv_bundle.py` starts at `paper/main.tex` and expands, recursively and after stripping
comments, every

| form | example in this manuscript |
|---|---|
| `\input{...}` / `\include{...}` | `\input{sec_9.tex}`, `\input{carried/fig_sqw.tex}` |
| `\includegraphics[...]{...}` | `\includegraphics[width=\linewidth]{figs/fig_hero.pdf}` |
| `table{...}`, `table[...]{...}` | `table[x=w,y=Aex]{aw_method.dat}` |
| `table {...}` **with no brackets, on its own line** | `table {n3_pub_frac_lo.dat}` in `figs/fig_thm1iii_violation_native.tex` |

Search paths are read *out of* `main.tex` (`\graphicspath{{figs/}}` and
`\pgfplotsset{table/search path={figs}}`), not assumed.

**The trap, and the measurement.** An earlier attempt wrote the pattern as
`table\[[^]]*\]\{` — it *required* the optional bracket group — and lost `.dat` files without a
word. `--selftest` control 1 re-runs that broken derivation on today's tree and prints the damage:
**60 files instead of 66, six of the ten `.dat` files gone**
(`n3_pub_eta_hi/lo`, `n3_pub_frac_hi/lo`, `n3_thmB_eta`, `n3_thmB_frac`). The parser here is a
scanner, not a line regex: keyword → *optional* `[...]` → `{...}`, across newlines.

**An unresolved reference is a hard error**, never a silent omission. That is the whole difference.

### The closure is confirmed by two mechanisms that share no code with it

1. **Filename-shaped-literal sweep.** Every `{...}` group in every shipped `.tex` whose content looks
   like `name.dat|csv|tsv|txt|pdf|png|jpg|eps` must already be in the closure. Catches a macro form
   the scanner does not know about.
2. **TeX's own testimony.** The clean-room compile runs `pdflatex -recorder`; every `INPUT` line of the
   resulting `main.fls` that points inside the extraction directory must be in the closure, **and vice
   versa**. This is the oracle that cannot be fooled by a parser bug, because TeX wrote it.
   Current result: **66 opened, 0 opened-but-not-shipped, 0 shipped-but-never-opened**, out of 719
   `INPUT` lines (the rest are the class, the packages and the fonts).

   *And the room is checked for leaks.* An `INPUT` line outside the extraction directory is normally a
   texmf file — but if `TEXINPUTS` or a stray absolute path reached back into the repository, the
   compile would succeed on a file **the package does not contain**, and this whole test would pass
   for the wrong reason. Any such path is now named rather than skipped. Current result: **0 — the
   clean room is sealed.**

The script also prints the **15 files left behind** in `paper/` (`sqw_edges.dat`, `hero_aw.dat`,
`fig_amort_native.tex`, `main.pdf`, the v2 tarball, …) so that the *opposite* failure — shipping
extraneous files, which arXiv explicitly forbids — is visible too.

---

## 2. What the package contains

**66 files, ~2.5 MB uncompressed, ~2.0 MB gzipped.**

| kind | count |
|---|---|
| `.tex` | 50 — **21 at the root** (`main.tex`, 15 `sec_*`, `app_carried_repro`, `app_carried_theorem`, `app_stats`, `bibliography`, `table_molecules_v3`), **20 in `figs/`**, **9 in `carried/`** |
| `.pdf` (vector figures) | 6 — `fig_akw_native`, `fig_circuit`, `fig_hero`, `fig_noise_score`, `fig_spinqw`, `fig_sqw` |
| `.dat` (plotted tables) | 10 — `aw_method`, `heron_hot`, and the eight `n3_*` series of Fig. 5 |

The tarball is built **deterministically** (member order sorted, mtime/uid/gid/mode fixed, gzip
`mtime=0`), so two builds of the same tree are byte-identical — checked on every run.

---

## 3. Verification — what was actually measured, 2026-09-19

Everything below is one run of `python src/build_arxiv_bundle.py --passes 3`.

### 3.1 Compiles from a clean extraction, outside the repository

```
pdflatex exit 0 | pages 59 | errors 0 | undefined refs 0 | undefined citations 0 | missing files 0
```

Three passes, extraction into a fresh temporary directory that contains nothing but the package.

> **A detector that was blind, and is not any more.** This script asks `pdflatex` for
> `-file-line-error`, which makes errors print as `sec_9.tex:120: LaTeX Error: …` instead of
> `! LaTeX Error: …`. The first version of `log_stats()` only matched `! `, so it scored a run with six
> real errors as **“0 errors”** — caught on the DVI probe below. It now matches both shapes, and
> `--selftest` controls 2 and 2b create genuine errors and demand that *this very function* sees them.

### 3.2 The PDF is equivalent to the repository's

Snapshot of the 12:46 run of 2026-09-19 (the hashes move with every source edit — re-run, do not quote):

| compared | pages | chars of extracted text | text SHA-1 |
|---|---|---|---|
| package, clean room, 3 passes | 59 | 274 969 | `80715e7724440fcc` |
| full `paper/` tree, same 3 passes, same binary | 59 | 274 969 | `80715e7724440fcc` |

**Identical text flow on every one of the 62 pages.** The mirror recompilation is the comparison that
matters: same compiler, same pass count, same minute — so any difference would be the *package's*
fault and nothing else's.

The script also compares against the committed `paper/main.pdf`. When that differs, it does **not**
guess which side is wrong: it checks whether any packaged source is *newer* than `main.pdf` and names
them. At the time of writing it did differ, on pages 43, 44, 45 and 55, and the diagnosis was
unambiguous — five sources (`figs/fig_akw_sampled_honest.tex`, `sec_9.tex`, `main.tex`,
`app_carried_repro.tex`, `sec_10.tex`) had been edited **after** `main.pdf` was last built. The
committed PDF was stale; the package was right.
**Rebuild `paper/main.pdf`, then re-run this script, before uploading.**

Page geometry, measured from the produced PDF rather than assumed: **612 × 792 pt (US Letter) on every
page**, from REVTeX's `aps` substyle.

### 3.3 Nothing retired, fabricated or auxiliary is in the package

Run as code against **the tarball**, not against the tree — the bytes that would actually be uploaded:

| checked | result |
|---|---|
| 15 retired basenames (7 v2 body files + 5 withdrawn figures + `mainNotes.bib`, `main.pdf`, the v2 tarball) | none present |
| 22 forbidden extensions (`.aux .log .out .toc .fls .bbl .blg .bib .py .ipynb .gz …`) | none present |
| 7 forbidden directories (`_superseded`, `__pycache__`, `.git`, `rebuild`, `ckpt`, `notebooks`, `data`) | none present |
| the two fabricated `sqw_edges.dat` rows (`0.0000 4.969 5.000 5.000`, `2.0000 …`) | absent — `sqw_edges.dat` is not in the package at all |
| superseded `19184` / `824504`, matched on numeric word boundaries | 0 occurrences |
| stray `CR` (0x0D) bytes — the defect class that typeset a build note on page 43 | 0, checked at byte level |
| absolute paths (`C:\`, `/home/`, `/tmp/`, `/w/`) in any shipped `.tex` | 0 |
| every member byte-compared against its source in `paper/` | all 66 identical |
| `fig_amort` / `fig_gflow` reachable from *live* (uncommented) code in `main.tex` | no |

### 3.4 The five negative controls (`--selftest`)

A PASS with no negative control is not evidence. All five fire:

| control | result |
|---|---|
| 1. the old bracket-requiring parser | **FIRED** — 60 files instead of 66, 6 `.dat` lost |
| 2. delete `figs/aw_method.dat` | **FIRED** — 5 errors, first `figs/fig_method_native_frag.tex:19: Package pgfplots Error:` |
| 2b. delete `sec_10.tex` | **FIRED** — 3 errors, first `! LaTeX Error: File 'sec_10.tex' not found.` |
| 3. re-introduce `19184` into `fig_scaling2_native.tex` | **FIRED** — audit names the file and the line |
| 4. insert a sentinel paragraph into `sec_10.tex` | **FIRED** — the PDF text comparison sees it |

---

## 4. arXiv's requirements — re-read from the source, not recalled

**Checked against <https://info.arxiv.org/help/submit_tex.html> and
<https://info.arxiv.org/help/00README.html> on 2026-09-19.** This mattered, because one of the three
rules this repository had been carrying is **no longer true**.

| rule this project believed | status today | what the package does |
|---|---|---|
| “`\pdfoutput=1` in the first five lines” | ❌ **OBSOLETE, and now the opposite.** arXiv retired AutoTeX in **April 2025** (Submission System 1.5). The page now says verbatim: *“You should not use `\pdfoutput` to change the output format.”* The processor is chosen by the submitter at upload. | `main.tex` does **not** set `\pdfoutput` — which is now the correct state. **Do not add it.** |
| “explicit paper size” | ➖ **Not a rule on the 1.5 system** — it was a `dvips`-era concern. No paper-size requirement appears on the current pages. | Measured anyway: every page **612 × 792 pt**, set by `revtex4-2`'s `aps` substyle. Nothing to do. |
| “no `\today`” | ✅ **REAL and current.** Verbatim: *“arXiv recommends against using the `\today` macro in the standard `\date` field. Because pdf are occasionally rebuilt this date will change and may cause confusion.”* | ✅ **CLOSED 2026-09-19.** The `\date{\today}` line was deleted from `paper/main.tex` (a comment records why). The fix had been verified against the class rather than guessed -- `revtex4-2.cls:2521` is `\def\@date{}` -- and the composed PDF confirms it: `Dated:` occurs **0** times in the extracted text of all 57 pages. |

### What replaced the `\pdfoutput` rule, and why it still bites here

The engine is now selected in the upload form, and recorded in a `00README.json`
(`{"process": {"compiler": "pdflatex"}, "sources": [{"filename": "main.tex", "usage": "toplevel"}]}`).
arXiv says it is “not required (nor recommended)” to write that file by hand, so the package does not
ship one — **but the choice must be made correctly in the form**, because this package cannot survive
the other path. Measured, not assumed: the script re-compiles the extracted package through
`latex` (the DVI route), and it dies —

```
latex main.tex -> 6 error(s); first: carried/fig_hero.tex:2: LaTeX Error: Cannot determine size of graphic in figs/f...
```

one per embedded PDF figure, because `graphics` loads `dvips.def` in DVI mode. **Select “LaTeX with
PDFLaTeX”.**

### Further rules verified against the package

| rule (verbatim from `submit_tex.html`) | result |
|---|---|
| “do not include any associated auxiliary file or intermediate or resulting output file, e.g. `foo.aux, foo.log, foo.toc, foo.lot, foo.lof, foo.dvi, foo.pdf`” | PASS — none present |
| “do not include extraneous files (including unused figure files), leftover files, backup files” | PASS — the closure is the definition of *used*; 15 unused files are left behind, and listed |
| “we require that you use a unified figure format” | PASS — 6 PDF figures, 0 EPS/PS |
| “Hidden files will be deleted upon announcement” | PASS — no dotfiles in the package |
| “Submissions with embedded JavaScript are automatically rejected” | PASS — all 6 figure PDFs scanned for `/JavaScript`, `/JS`, `/RichMedia`, `/Movie`, `/Launch`; clean |
| “In case the `.bbl` file is not uploaded, and at least one necessary `.bib` file is missing, the submission system will block you” | N/A — the bibliography is an inline `thebibliography` in `bibliography.tex`; no `\bibliography{}`, no `.bib`, no `.bbl` |
| compilation runs from the root of the submission directory | PASS — `main.tex` is at the package root; `figs/` and `carried/` are relative to it |
| only one `\documentclass` in the package | PASS |

*Optional, not done:* arXiv notes that INSPIRE's reference extraction is more accurate when the
references live in a file named `main.bbl` rather than `bibliography.tex`. That would be an edit to
`paper/main.tex`, so it is left as a decision, not a change.

---

## 5. Submission form

This is a **replacement (v3)** of an existing entry, not a new submission.

| field | value |
|---|---|
| **Entry** | `arXiv:2608.16436` → **Replace** |
| **Processor** | **LaTeX with PDFLaTeX** (see §4 — the DVI route fails on all six PDF figures) |
| **Upload** | `build/arxiv-v3.tar.gz` |
| **Primary category** | `quant-ph` |
| **Cross-list** | `cond-mat.str-el` |
| **Title** | Captured weight and boundary leakage bound the error of sample-based spectral functions |
| **Authors** | Nicolás Bonilla Vargas |
| **Abstract** | copy the `abstract` environment from `paper/main.tex` (**it changed in v3**). **Type the ten Greek letters as Unicode characters (η, ρ, γ, Λ) and keep the exponents and subscripts (`10^-3`, `L_1`): so typed it is 1910 characters against arXiv's limit of 1920. Spelled out as `eta`, `rho`, ... it is 1940 and arXiv rejects it.** The margin depends on the convention, so decide the exact string before opening the form; `scratchpad/cierre/abs_paste.py` prints it and counts it both ways. Note that `mide_resumen.py` deletes `^` and `_` and therefore reads 8 characters short. |
| **License** | CC BY 4.0 |
| **Comments** | paste the block below verbatim (measured 2026-09-19: 62 pages, 17 figures, 17 tables). arXiv asks that a replacement merge the old comments with the new ones, and that the reason for the replacement appear in this field. |

```
62 pages, 17 figures, 17 tables. Code, data and full reproduction repository:
https://github.com/nicolasbonilla/dynamical-spectral-functions-sqd
v3: major revision, retitled. Theorem 1(iii) of v1-v2, an error bound indexed on the
captured Born weight alone, does not hold - a two-level counterexample forces its
constant to the trivial value - and is replaced by a proved two-sided leakage bound,
evaluated on five Hubbard sectors up to dimension 10306296 and vacuous on all of them
by factors 1.66-8.61. The two rank correlations that v1-v2 pooled over the thirty
Hubbard points are Simpson reversals of the six (L,N) strata and are replaced by a
stratified test. Further claims of v1-v2 are not made here and are listed, item by
item, in docs/KNOWN_DISCREPANCIES.md of the repository above.
```


### Steps

1. Rebuild `paper/main.pdf` from the current sources, then run
   `python src/build_arxiv_bundle.py --passes 3` and require **`RESULT: PASS`**.
2. Sign in at <https://arxiv.org>, open `arXiv:2608.16436`, choose **Replace**.
3. Upload **`build/arxiv-v3.tar.gz`** — *not* `paper/arxiv-submission.tar.gz`.
4. When asked for the processor, choose **LaTeX with PDFLaTeX**.
5. Review arXiv's generated PDF: **62 pages**, and the last body page must read `X. CONCLUSION`
   with no stray build note above it.
6. Update the abstract and the Comments field; keep CC BY 4.0.
7. After announcement, update the badge and `CITATION.cff` with the v3 date.

---

## 6. Open items that are *not* in this package's gift

| # | item | owner |
|---|---|---|
| ~~1~~ | ~~`\date{\today}` in `paper/main.tex` (§4)~~ — **done 2026-09-19**; the composed PDF carries no date line. | closed |
| 2 | `paper/main.pdf` must be rebuilt after the last source edit, or §3.2's second comparison stays red | whoever compiles last |
