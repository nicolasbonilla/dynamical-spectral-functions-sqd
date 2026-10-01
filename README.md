<div align="center">

# Captured weight and boundary leakage bound the error of sample-based spectral functions

### *a two-sided, a posteriori bound — and an honest account of where it is vacuous*

**Nicolás Bonilla Vargas** &nbsp;[![ORCID](https://img.shields.io/badge/ORCID-0009--0006--6155--4391-A6CE39?logo=orcid&logoColor=white)](https://orcid.org/0009-0006-6155-4391)

[![arXiv](https://img.shields.io/badge/arXiv-2608.16436-b31b1b.svg)](https://arxiv.org/abs/2608.16436)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22967220.svg)](https://doi.org/10.5281/zenodo.22967220)
[![Paper](https://img.shields.io/badge/paper-PDF-blue.svg)](paper/main.pdf)
[![License: MIT](https://img.shields.io/badge/code-MIT-green.svg)](LICENSE)
[![License: CC BY 4.0](https://img.shields.io/badge/paper%20%26%20data-CC--BY--4.0-lightgrey.svg)](LICENSE)
[![Repository](https://img.shields.io/badge/repository-github-181717?logo=github)](https://github.com/nicolasbonilla/dynamical-spectral-functions-sqd)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](requirements.txt)
[![Reproducible](https://img.shields.io/badge/reproducible-gaps%20declared-brightgreen.svg)](docs/REPRODUCE.md)

</div>

> **One sentence.** Sample-based quantum diagonalization is moving from energies to dynamics; given the
> subspace the shots populated, we prove a **two-sided bound on the relative `L₁` error** of the
> spectral function built on it — at least `1−w`, at most `(1+√w)(√(1−w) + Λ̂(η)/η)`, with `w` the
> captured Born weight and `Λ̂` the normalized coupling of `H` across the subspace boundary at
> resolution `η`, computable from the Ritz vectors in hand. **Shots buy the first term down, not the
> second.**

> **And the retraction, in the same breath.** Versions 1–2 of this preprint claimed a *weight-only*
> bound. An analytic counterexample forces its constant to the trivial value, and it is **withdrawn**.
> On Hubbard sectors up to dimension 10 306 296 the replacement is worse than the trivial bound at
> every operating fraction of the resource scan, at all five sizes, by factors of 1.66–8.61 — **so those
> reconstructions are uncertified**; the one reconstruction drawn at 85% of its sector (Fig. 3) is the
> certified exception, loosely. What survives is a calibration: over sixteen `(L,η)` cells whose true errors span five
> decades, the bound crosses the trivial one at a true relative error of `1.2–7.3×10⁻³`, median
> `2.8×10⁻³`, computed without the answer. **Yet a two-constant rule in `1−w`, fitted out of sample,
> localizes that error more tightly on 11 of 14 train/test splits, so no superiority is claimed for the
> bound as an error tracker.** **No quantum advantage is claimed.**

**Repository:** <https://github.com/nicolasbonilla/dynamical-spectral-functions-sqd> — this URL is the
one the manuscript's data-availability statement points to. It is also in
[`CITATION.cff`](CITATION.cff).

This repository is the **reproduction record** of the paper: every number and every figure traces to a
script and a data file, except the declared gaps below, listed in [`docs/REPRODUCE.md`](docs/REPRODUCE.md) and
[`docs/FILE_INDEX.md`](docs/FILE_INDEX.md), and `python src/verify.py` recomputes the physics and checks
the deposit against it. **It is not a claim that everything re-runs here.** What a clean clone reaches —
and the parts it does not, because they need `pyscf`, an IBM Quantum account, or a figure step that is
not deposited — is tabulated, measured rather than asserted, at the end of
[`docs/REPRODUCE.md`](docs/REPRODUCE.md), and the gaps are enumerated in
[`docs/KNOWN_DISCREPANCIES.md`](docs/KNOWN_DISCREPANCIES.md). The statement that may be published about
all this is [`docs/DATA_AVAILABILITY.md`](docs/DATA_AVAILABILITY.md).

> **Companion works.** Review — *Machine learning for sample-based quantum diagonalization: a review of
> generative configuration recovery and the classical-simulability frontier* ([arXiv:2608.05314](https://arxiv.org/abs/2608.05314)).
> Method — *Reconstruction-independent moment tests for quantum-computed dynamical spectra* (in preparation).

---

## Figure gallery

| Momentum-resolved `A(k,ω)` (exact Mott map) | The resource thesis (74 rows, 71 distinct exact states) |
|:---:|:---:|
| ![akw](docs/img/akw.png) | ![resource](docs/img/resource.png) |
| **Charge structure factor `S(q,ω)`, L=12** | **Spin structure factor `S^zz(q,ω)`, L=12** |
| ![sqw](docs/img/sqw.png) | ![spinqw](docs/img/spinqw.png) |
| **Scaling: fraction falls, absolute cost grows** | **IBM Heron hardware (two runs, stated at their size)** |
| ![scaling](docs/img/scaling.png) | ![hardware](docs/img/hardware.png) |
| **The geminal witness: magic large, `χ = 2`** | **N₂ dissociation hero (`F₁ = 2Nᵤ`)** |
| ![witness](docs/img/witness.png) | ![hero](docs/img/hero.png) |

<div align="center"><em>Both collective channels come from the <strong>same sampling primitive</strong>. At L=12 each is a set of discrete poles, so the paper uses the spin–charge contrast (spin scale πJ/2≈0.79t, far below the charge gap Δ≈4.97t) as a check of the method, not as a result, and calls neither channel a continuum or gapless (Sec. VIII).</em></div>

---

## What the paper shows

1. **A two-sided, a posteriori error bound for sampled spectral functions.** One sampling primitive,
   four channels — `A(ω)`, `A(k,ω)`, `S(q,ω)`, `S^zz(q,ω)` — with **no Hadamard tests or controlled
   unitaries**, and a certificate computed from the Ritz vectors the run already produced. The upper
   branch splits the error into what shots can fix (captured weight `1−w`) and what they cannot
   (boundary leakage `Λ̂(η)/η`).
2. **A retraction, and what replaces it.** The weight-only bound of versions 1–2 is withdrawn: an
   analytic counterexample forces its constant to the trivial value, and none of 606 test subspaces
   violates the replacement. The obstruction is structural — in the determinant basis the probe's
   Krylov space has coordinate support on every symmetry-allowed determinant, measured from `L=6` to
   `L=14`, so at full captured weight only a subspace holding all of them has zero leakage.
3. **An obstruction theorem for one-body magic.** The fermionic AntiFlatness collapses to a single
   1-RDM invariant `F₁ = 4 tr[γ(1−γ)] = 2Nᵤ`. Being an orbital-rotation (Gaussian) invariant, it is
   **provably decoupled** from the basis-dependent support `|S|`: constant along an orbit on which the
   support runs from `O(K)` to `2^K` at bond dimension 2. Within each of six `(L,N)` strata its rank
   correlation with the support is exactly `−1` (exact `p = 3.3×10⁻¹³`); **both pooled coefficients are
   withdrawn** as a Simpson reversal.
4. **Validation, measured prices, and hardware.** Exact diagonalization across Hubbard chains and a
   nineteen-molecule suite; the three prices (support, resolution, shots) measured rather than
   asserted; and two real **IBM Heron** runs.

**Hardware, stated at its size.** Two runs were executed on IBM Heron processors. On `ibm_fez` (L=6
Hubbard, 12 qubits), post-selection retained all 300 determinants of the sector, so the reconstructed
A(ω) equals the exact one *by coverage*: an identity, not a fidelity measurement. The notebook read the
bitstrings in reversed qubit order, so every retained shot is a device error: a noiseless simulation of
the same circuits keeps 0 of 350 000 shots under that reading and all of them under the correct one
(`src/hw_bitorder_check.py`; `docs/KNOWN_DISCREPANCIES.md` §30). On `ibm_marrakesh`
(stretched N₂, CAS(10e,12o), 24 qubits), recovery over 8 seeds reaches 0.59 ± 0.14 mHa of the exact
active-space energy, below the 29.5 mHa of a noiseless simulation of the same circuit. The recovered
subspace sizes on which that comparison turns were not stored, so the record cannot say whether device
noise helps. The `ibm_fez` run says nothing about accuracy at incomplete coverage.

---

## Repository structure

```
.
├── src/                       # 66 computation & figure scripts, plus 29 in src/frontier/ (the C3 sweep); run from repo root as `python src/X.py`
│   ├── verify.py              #   the adversarial guardian: recomputes the physics, then checks the deposit
│   ├── check_figures.py       #   do the committed generators still produce the committed fragments?
│   ├── check_provenance.py    #   does the documentation still describe the manuscript? (derives it from main.tex)
│   ├── akw_lanczos.py …       #   exact-diagonalization engines + spectral/resource/scaling compute
│   ├── make_*.py              #   one data-driven figure generator per figure (no hand-typed numbers)
│   └── recovered/             #   the generators of Figs. 6 and S5 and their compute, recovered 2026-09-26 (README inside)
├── data/                      # 38 datasets, plus data/c3_frontier/ (95 files), data/gflow_runs/ (20) and
│                              #   data/thm1iii_violation/ (the inputs of Fig. S5): every plotted number (*.json)
├── ckpt/                      # the L=14 ground-state checkpoint (94 MB) and ranking, with SHA-256 (deposited 2026-09-25)
├── notebooks/
│   ├── 00_Reproduce_Everything.ipynb   #   the v1–v2 narrated pipeline, kept as a record (needs pyscf)
│   ├── HW_LUCJ_N2_Heron_READY.ipynb    #   IBM Heron N₂ energy   (token SCRUBBED — see below)
│   └── Spectral_Heron.ipynb            #   IBM Heron A(ω)        (token SCRUBBED — see below)
├── paper/                     # the v3 manuscript: REVTeX 4.2 (aps, pra), one column for submission
│   ├── main.tex               #   preamble, abstract and the \input list
│   ├── sec_1.tex … sec_10.tex #   the ten sections of the main text; app_proof.tex is its Appendix A
│   ├── sm_*.tex               #   nine Supplemental Material files: Secs. S1–S8 (sm_begin.tex opens the SM)
│   ├── sec_*_app.tex, app_*.tex#  Secs. S9–S16 of the Supplemental Material (the former appendices)
│   ├── carried/               #   nine one-float wrappers carried verbatim from v2
│   ├── figs/                  #   11 native fragments, 6 figure PDFs, 19 plotted .dat tables, captions
│   │   └── src/               #   the standalone sources of the 6 figure PDFs (deposited 2026-09-26)
│   ├── main.pdf               #   the compiled preprint (68 pp)
│   ├── main_pra.tex, sm_pra.tex#  the Physical Review A split: main text + references, and the SM on its own
│   ├── main_pra.pdf, sm_pra.pdf#  their builds (33 pp and 38 pp), made by pra_split.py
│   ├── pra_split.py           #   builds the two drivers and checks them page by page against main.pdf
│   └── arxiv-submission.tar.gz#   the FROZEN arXiv v2 bundle — a record, never a source
├── _superseded/               # the v2-era manuscript, the five figures v3 withdrew and the v1 README thumbnails, kept for the trail
├── docs/
│   ├── REPRODUCE.md           #   figure/number → script → exact command, and what a clean clone cannot reach
│   ├── FILE_INDEX.md          #   every file, described, plus the orphan inventory
│   ├── FIGURE_PROVENANCE.md   #   figure → source → generator → data → job id
│   ├── KNOWN_DISCREPANCIES.md #   where running a script does NOT reproduce the deposit (read this)
│   ├── DATA_AVAILABILITY.md   #   the statement drafted for PRA, and the audit behind each sentence
│   ├── C3_FRONTIER.md         #   the C3 certificate-frontier sweep behind Sec. V, file by file
│   ├── ARXIV_SUBMISSION.md    #   how the v3 source package was built and submitted
│   └── img/                   #   README thumbnails, rendered from the v3 figures by make_thumbnails.py
├── .github/workflows/ci.yml   # CI: `make verify`, `make check-figures` and `make check-coherence` on every push
├── requirements.txt           # Python dependencies
├── Makefile                   # `make verify`, `make check-figures`, `make figures`, `make paper`, `make pra`, `make ci`
├── CITATION.cff               # citation metadata
├── .zenodo.json               # Zenodo deposit metadata
└── LICENSE                    # MIT (code) + CC-BY-4.0 (paper text, figures and data)
```

---

## Quick start

```bash
# 1. environment
python -m venv .venv && source .venv/bin/activate      # (Windows: .venv\Scripts\activate)
pip install -r requirements.txt

# 2. the adversarial guardian (~10 s, numpy/scipy only). It RECOMPUTES the physics
#    (Hubbard F₁/E₀ for both boundary conditions, the geminal witness) and then checks every
#    deposited artefact — figures, tables, .dat files and the abstract — against it.
python src/verify.py        # or: make verify   ->   PASS, plus a named list of open defects
#    It exits non-zero and names the defect if anything disagrees. The `[XFAIL]` lines it
#    prints are REAL open defects, claims and coverage gaps, deliberately visible.

# 2b. do the committed figures still come out of the committed generators? (see the caveat below)
python src/check_figures.py # or: make check-figures

# 2c. does the DOCUMENTATION still describe the manuscript? (derived from paper/main.tex)
python src/check_provenance.py          # add -v for the orphan inventory, --markdown for the table

# 3. the v1–v2 narrated pipeline, kept as a record (needs pyscf; NOT the v3 figure set)
jupyter notebook notebooks/00_Reproduce_Everything.ipynb    # or headless: make reproduce

# 4. regenerate a single data-driven figure fragment from its committed data
python src/make_decoupling_native.py    # -> paper/figs/fig_decoupling_native.tex

# 5. build the paper (needs a TeX distribution with REVTeX 4.2 and pgfplots ≥ 1.18)
make paper                              # -> paper/main.pdf

# 6. the two files Physical Review A takes (after step 5: it checks them against paper/main.pdf)
python paper/pra_split.py --build       # or: make pra  -> paper/main_pra.pdf, paper/sm_pra.pdf
```

**The Physical Review A split.** `paper/main.tex` builds the main text, the references and the
Supplemental Material as one document — the arXiv form. PRA takes the Supplemental Material as a
separate file, so two drivers in `paper/` split that same document without copying any of it:
`main_pra.tex` runs `main.tex` and stops where the SM begins (main text and full reference list,
33 pp.), and `sm_pra.tex` runs `main.tex`, skips to the SM and ends with a reference list of only the
works the SM cites (38 pp.: 35 of SM, 3 of references). Cross-references between the two print as in
`main.pdf` ("Sec. S3", "Fig. 9") through `xr-hyper`, whose `[nocite]` option needs an `xr-hyper` of 2023
or later. `python paper/pra_split.py --build` builds both in a temporary copy of `paper/` (it writes
only the two PDFs into `paper/`) and checks them: 0 errors, 0 undefined or multiply defined labels and
citations, no overfull box or warning that the `main.tex` build does not have, the SM reference list
against the SM's citations, the APS form of the Supplemental Material reference (ref. [19] of the main
list since the reordering of 2026-09-30, "See Supplemental Material at [URL will be inserted by publisher] …, which includes Refs.
[104–142]", the 39 works cited only in the SM), and every page of `paper/main.pdf` in exactly one of the
two PDFs — pixel-identical below the page number for the main text, identical up to citation numbers
for the SM. Build `paper/main.pdf` first: that is the PDF the pages are compared against.

**How v3 was built and submitted:** see **[`docs/ARXIV_SUBMISSION.md`](docs/ARXIV_SUBMISSION.md)** (v3 was
submitted on 2026-09-25 from a package built by `src/build_arxiv_bundle.py`, and was announced on
2026-09-28). **`paper/arxiv-submission.tar.gz` is the frozen snapshot of what was posted as
v2 and must never be re-uploaded**: 17 files in `paper/` have since moved ahead of it, 21 of its files
no longer exist in `paper/`, and the bundle still contains the two fabricated rows of `figs/sqw_edges.dat` and the
seven v2 body files the v3 sections replace. Run `python src/check_tarball.py` for the live comparison,
and see [`docs/KNOWN_DISCREPANCIES.md`](docs/KNOWN_DISCREPANCIES.md) §7.

---

## Reproducibility at a glance

*Measured 2026-09-19 unless the row says otherwise; the guardian rows were re-measured on 2026-09-26.*

| Layer | Reproducible here? | How |
|---|---|---|
| **The guardian** | ✅ **PASS** — 1280 checks, 0 failures, 14 563 numeric assertions, 4 xfail, 0 xpass, 2 skips *(re-measured 2026-09-30, final review: +66 checks and +66 assertions, verify.py checks the Fig. 3 caption's panel (b) mean and maximum, window and branch thresholds (section 9.6b), the printed noise values of Figs. S7-S8, Sec. S8 and Sec. VIII (section 9.9b) and every restated Fig. 3 certificate number in the abstract and Secs. I, V, X and S3 (section 9.9c); 1214 / 14 497 earlier on 2026-09-30: +2 checks and +2 assertions, verify.py section 9.3 now checks every printed exact permutation p (Sec. I, Sec. VI, Fig. 6 caption) against `exact_perm_p_one_sided_F1` at 1e-15, unified on the 3.35e-13 rounding; re-measured 2026-09-28, late, closing the review: +278 checks and +2 821 assertions for section 9.6c (the free-fermion supports, the spin-peak deviations and the Fig. S5 pool re-classified from the deposited certificate rows), the chi_max of the witness, the LF hash of the Fig. 3 generator and the L=14 re-ranking record, and the declared Fig. S5 coverage gap closed, xfail 5 → 4, `KNOWN_DISCREPANCIES.md` §33; 934 / 11 674 after the finite-shot numbers of the Fig. 3 caption were deposited with their generator and checked from its two data files: +60 checks, +187 assertions; 874 / 11 487 after the Fig. 3 frequency grid was widened to [-9t, 17t]: 874 pass, 0 fail, 11 487 assertions, 5 xfail, 2 skips; 874 / 9 891 after two phrase guards on Sec. S8 were added the same day; 872 / 9 889 on 2026-09-26, after the hardware bit-order check)* | `python src/verify.py` (or `make verify`) |
| **The documentation** | ✅ **PASS** — derived from `paper/main.tex`: 61 source files, 34 floats, 17 figures, 17 tables *(2026-09-28; 60 on 2026-09-26, before `app_proof.tex`)* | `python src/check_provenance.py` |
| **Manuscript coherence** | ✅ **PASS** — 10 cross-section claims, 20 synthetic controls, all 20 fire *(2026-09-19)*; since 2026-09-30 also 16 per-place checks (each vacuity statement in each section that carries it), each with its own control, all 16 fire. New that day, because the other three guardians were green while seven sections contradicted each other: they compare a printed number against a deposited file, not a claim against a claim | `python src/check_coherence.py` |
| **The figure generators** | ✅ **PASS on 21 artefacts from 9 generators** *(2026-09-28, late: the Fig. S5 whisker script `src/recovered/n3_fig5_whiskers.py`, with the fragment and its summary JSON; 19 from 8 on 2026-09-26)* — the six of 2026-09-19, the two recovered on 2026-09-26 into `src/recovered/`, and the whisker script; one provenance-timestamp comment line of the Fig. 6 fragment is normalised, and the run says so. Measured with four negative controls in `KNOWN_DISCREPANCIES.md` §29. (Until 2026-09-19 one of the first eight artefacts was compared with a copy of itself; §25.) | `python src/check_figures.py` |
| **Exact diagonalization, numpy/scipy only** (`A(ω)`, `A(k,ω)`, `S(q,ω)`, `S^zz`, the two gaps, the χ–\|S\| resource map, the geminal witness, the scaling) | ✅ **verified locally**, minutes | `python src/<script>.py`; see `docs/REPRODUCE.md` |
| **The 19-molecule suite** | ⚠️ runs, but needs `pyscf`; **untested in this pass**, so the molecular checks are transcription checks, not recomputations | `python src/n19_suite.py` etc. |
| **The v1–v2 narrated notebook** | ⚠️ **a record, not the v3 reproduction path** — it reproduces the v1–v2 figure set, several of which v3 withdrew, needs `pyscf` and was **not run**; use `src/verify.py` and `docs/REPRODUCE.md` | `notebooks/00_Reproduce_Everything.ipynb` |
| **Matrix-free scaling** to 28 qubits | ✅ (RAM-staged, resumable) | `src/scaling_lanczos_mf.py` |
| **Figures** | ⚠️ **8 of 17 are guarded by regeneration**, two of them only in part: Figs. 2, 3, 9, S2, S3 and S8 byte-for-byte; Fig. 6 byte-for-byte after one provenance-timestamp comment line is normalised; Fig. S5 in its ten `.dat` tables, its whisker blocks and its summary JSON (since 2026-09-28, late; the legend literals of the fragment are checked by `build_n3.py`'s self-check and every plotted ratio is re-derived from the certificate rows by `verify.py` section 9.6c). Nine are not guarded: four have no deposited generator (Figs. 1, S1, S4, 5), and five ship as a PDF (Figs. 4, 7, 8, S6, S7). The standalone sources of all six PDF figures (Fig. 1 included) are deposited since 2026-09-26 and recompile pixel-identically, but no script guards them, and none of the three rasters they embed can be regenerated here: the viridis input of `src/recolor_akw_v2.py` (which writes `akw_field_v2.png`) is not deposited, the step that recoloured `sqw_field.png` is not deposited, and no deposited script writes `spinqw_field.png`. *(Until 2026-09-26: 6 of 17, and the generators of Figs. 6 and S5 were not in the tree.)* | inventory in [`docs/FIGURE_PROVENANCE.md`](docs/FIGURE_PROVENANCE.md) |
| **IBM Heron hardware runs** | ⚠️ needs an IBM Quantum account, and the jobs are closed | `notebooks/` (tokens scrubbed; **no raw counts are deposited**) |

> **`paper/main.pdf` is the manuscript for Physical Review A**, the v3 text revised in the review pass of
> 2026-09-28 and in the passes of 2026-09-30 (listed at the end of this note), last recompiled on
> 2026-09-30 after the final source edit: 68 pages, one column, a 28-page
> main text (p. 28 holds the end of Sec. X, the acknowledgments and Appendix A, the proof of
> Theorem 1), the reference list on pp. 29–33, then the Supplemental Material from p. 34. Before that
> pass it was 66 pages, with a 26-page main text, references on pp. 27–31 and the SM from p. 32. *(This file said "a 25-page main text" until 2026-09-28, and
> so does the Comments field filed with arXiv v3, which is kept as that record; the v3 package itself,
> rebuilt from `build/arxiv-v3.tar.gz`, also ends its main text on p. 26.)* The build log shows 0 errors,
> 0 undefined references or citations, 0 overfull `\hbox` and 0 overfull `\vbox` (the layout pass of
> 2026-09-30, which changed no text, removed the two overfull `\vbox` of that day's earlier build, 24.3 pt
> on p. 48 and 5.0 pt on p. 55, and the report of Fig. 3 as 30.2 pt too large for its page,
> `KNOWN_DISCREPANCIES.md` §34; the p. 55 one came back when the final pass of that day reworded Sec. S5,
> and went again with a `\looseness=-1` on one SM paragraph, no word changed; before the 2026-09-30 rewording of Sec. S5 the p. 48 one was 12.8 pt; a second one of 10.9 pt on p. 43 went late on 2026-09-28 when the footnote of Sec. S5 was folded
> into its paragraph; 5.8 pt and 12.1 pt before the review pass). It differs from the source package submitted to arXiv as v3 in these later edits: the
> Zenodo DOI in the data-availability statement and in the repository reference (which now carries the
> current title), the sentence noting that the L=14 checkpoint is deposited, and the disclosure of
> 2026-09-26 that the L=6 `ibm_fez` run was post-selected in reversed bit order, without configuration
> recovery or readout-error mitigation (Secs. IX A and S8, Table I, the caption of Fig. 9;
> `docs/KNOWN_DISCREPANCIES.md` §30); and, on 2026-09-28, the sentence of Sec. S8 that still called
> that run a proof of principle with "recovery to the correct particle-number sector" (it now claims
> only that the steps executed, with a post-selection that was itself wrong), the known-gaps sentence of
> the data-availability statement (it now names what is still open, `KNOWN_DISCREPANCIES.md` §5, §13,
> §29), the guardian counts it quotes (934 / 11 674), and the Supplemental Material cited as a reference
> in the APS form (the works cited only in the SM were then the 38 of [103–140]; 39, [104–142], since
> 2026-09-30).
> **The review pass of 2026-09-28** then revised the text for Physical Review A: the abstract (the
> reconstructions are classical Born draws or the infinite-shot ranking, bar the device run; on 2026-09-30 the abstract was reworded to "one device run, an execution record, aside", since only the `ibm_fez` run enters a Hubbard spectral function; the bound
> presupposes the exact ground state), Secs. I and X to match; Theorem 1 gains hypothesis (H4), the exact
> probe, and its proof moves from Sec. S2 into a new Appendix A (`paper/app_proof.tex`); the proposition
> on sampling capture moves into Sec. IV C as Proposition 1 (the former Proposition 1 is Proposition 2);
> the coefficient-insensitivity statement of Sec. VI and Sec. S13 is corrected (the raw support is
> lower- but not upper-semicontinuous); Fig. 3 is redrawn and rescored on a window that holds both
> branches, its finite-shot run is deposited (`src/fig3_finite_shot.py`, recovered from the session that
> made it and checked to reproduce its original output exactly: all 32 rows at T = 2.6e6, the 28 recorded
> rows and every summary line at T = 5.2e6) and re-run on that window, and the Fig. S8
> bands are drawn unfloored (`docs/KNOWN_DISCREPANCIES.md` §31); the
> hardware, scope and noise captions describe the runs as executed (Sec. IX, Sec. S8, Figs. 9 and
> S6–S8); and the printed values that disagreed with a committed record are corrected (§32). The
> Supplemental Material is now first cited where it is first mentioned, at the end of Sec. I (ref. [12]
> then, ref. [19] since the reordering of 2026-09-30), which renumbered the main-text references first
> cited after Sec. I by one.
> **Closing that review, late on 2026-09-28** (`docs/KNOWN_DISCREPANCIES.md` §33 and the second table of
> §32): every number that was still declared rather than generated received a deposited generator and a
> `verify.py` recomputation — the four free-fermion supports of Secs. VI and S13
> (`src/free_fermion_support.py`), the spin-peak percentages of the Fig. 8 caption
> (`src/spin_peak_deviation.py`), the bond dimension of the geminal witness over all cuts
> (`src/apsg_witness.py`) — Fig. S5 was rebuilt from the deposited certificate scans (469 violations of
> 606, mildest factor 2.00, one set of numbers everywhere; `src/recovered/n3_fig5_whiskers.py` writes its
> whiskers, and the coverage gap `verify.py` had declared for it since 2026-09-19 is closed), the K=18
> ranking was re-run at L=14 and reproduces the stored selection as a set (`src/frontier/rerank_L14.py`),
> Secs. I and X claim of the geminal orbit only what Theorem 3 proves, and a dozen sentences were made
> exact against their records. Run `make paper` after any further source edit, then `make pra`.
> **The passes of 2026-09-30** (`docs/KNOWN_DISCREPANCIES.md` §34): the abstract, Sec. I and Sec. X were
> rewritten to state results rather than concessions (no result moved); a citation audit corrected what
> the text attributed to cited works (among them the priority for Gaussian-invariant measures of
> non-Gaussianity, now Sierant et al. for the antiflatness used here, with Tarabunga et al. and Leone and
> Bittel for what they prove; Benthien and Jeckelmann's sizes; the k = π cut of Nocera and Alvarez; the
> additive-approximation scope of ExtraFerm), and the reference list was put back in first-citation
> order (142 entries, the SM reference at [19], the SM-only works at [104–142]); the classical cost at
> generic rotated points of the geminal orbit is stated as not established, in the main text and the SM
> alike; and `verify.py` now also checks the Fig. 3 caption's panel (b) numbers, window and branch
> thresholds, the printed noise values and every restated Fig. 3 certificate number, while
> `check_coherence.py` requires each vacuity statement in each place that carries it.
> **The prose pass of 2026-10-01** (`docs/KNOWN_DISCREPANCIES.md` §35): the abstract, the opening of
> Sec. I and Sec. X were rewritten in short sentences, one idea per paragraph; no number, claim or limit
> changed in the paper, and three items moved from the abstract to the body (§35 lists them). After its
> two closing reviews the guardians pass, the three PDFs keep their page counts, and
> `build/arxiv-v4.tar.gz` was rebuilt locally from these sources (`RESULT: PASS`; nothing uploaded).
>
> **CI.** `.github/workflows/ci.yml` runs `make verify` and `make check-figures` on every push, and has
> done so since 2026-09-19 (see the Actions tab). `make check-coherence` was added as a third job on
> 2026-09-26 and runs from the next push on, so that the workflow runs what `make ci` runs.

**Data availability:** the statement drafted for Physical Review A, and the audit of what the deposit
does and does not support, are in [`docs/DATA_AVAILABILITY.md`](docs/DATA_AVAILABILITY.md).

The two hardware runs are the **entirety** of the quantum hardware used; every other result is an exact
classical simulation reproducible from this repository.

---

## 🔒 Security note (please read before pushing)

The hardware notebooks in `notebooks/` have had the author's **IBM Quantum API token and instance CRN
removed** and replaced with `<YOUR_IBM_QUANTUM_TOKEN>` / `<YOUR_IBM_QUANTUM_INSTANCE_CRN>` placeholders.
To run them, insert **your own** credentials (never commit them — use `QiskitRuntimeService.save_account`
locally or an environment variable). This repository contains **no secrets**.

---

## Citation

If you use this work, please cite it (see [`CITATION.cff`](CITATION.cff)):

```bibtex
@article{BonillaVargas2026CapturedWeight,
  title   = {Captured weight and boundary leakage bound the error of
             sample-based spectral functions},
  author  = {Bonilla Vargas, Nicol\'as},
  journal = {arXiv preprint arXiv:2608.16436},
  eprint        = {2608.16436},
  archivePrefix = {arXiv},
  primaryClass  = {quant-ph},
  doi     = {10.48550/arXiv.2608.16436},
  year    = {2026}
}
```

and, for the code and data, the archived deposit (this DOI always resolves to the latest version):

```bibtex
@software{BonillaVargas2026Repository,
  title     = {Captured weight and boundary leakage bound the error of
               sample-based spectral functions --- reproduction repository},
  author    = {Bonilla Vargas, Nicol\'as},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.22967220},
  url       = {https://doi.org/10.5281/zenodo.22967220},
  year      = {2026}
}
```

## License

- **Code** (`*.py`, `notebooks/`, `Makefile`): [MIT](LICENSE).
- **Paper text, figures and data** (`paper/`, `data/`): [CC-BY-4.0](LICENSE).
