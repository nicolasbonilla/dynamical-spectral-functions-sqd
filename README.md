<div align="center">

# Dynamical spectral functions from bitstring-sampled quantum subspaces
### *entanglement, not one-body magic, tracks the sampling cost*

**Nicolás Bonilla Vargas** &nbsp;[![ORCID](https://img.shields.io/badge/ORCID-0009--0006--6155--4391-A6CE39?logo=orcid&logoColor=white)](https://orcid.org/0009-0006-6155-4391)

[![arXiv](https://img.shields.io/badge/arXiv-2608.16436-b31b1b.svg)](https://arxiv.org/abs/2608.16436)
[![Paper](https://img.shields.io/badge/paper-PDF-blue.svg)](paper/main.pdf)
[![License: MIT](https://img.shields.io/badge/code-MIT-green.svg)](LICENSE)
[![License: CC BY 4.0](https://img.shields.io/badge/paper%20%26%20data-CC--BY--4.0-lightgrey.svg)](LICENSE)
[![Repository](https://img.shields.io/badge/repository-github-181717?logo=github)](https://github.com/nicolasbonilla/dynamical-spectral-functions-sqd)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](requirements.txt)
[![Reproducible](https://img.shields.io/badge/reproducible-every%20number-brightgreen.svg)](docs/REPRODUCE.md)

</div>

> **One sentence.** We lift sample-based quantum diagonalization (SQD/QSCI) from the ground-state *energy*
> to *dynamical* spectral functions — `A(ω)`, `A(k,ω)`, `S(q,ω)`, `S^zz(q,ω)` — reconstructed from purely
> bitstring-sampled subspaces, and prove that the one-body fermionic magic `F₁` is **decoupled** from the
> sampling cost `|S|`, which is instead **lower-bounded and empirically tracked** by the entanglement
> (minimal bond dimension `χ`).

**Repository:** <https://github.com/nicolasbonilla/dynamical-spectral-functions-sqd> — this URL is the
one the manuscript's data-availability statement points to. It is also in
[`CITATION.cff`](CITATION.cff).

This repository is the **reproduction record** of the paper: every number and every figure traces to a
script and a data file listed in [`docs/REPRODUCE.md`](docs/REPRODUCE.md) and
[`docs/FILE_INDEX.md`](docs/FILE_INDEX.md), and `python src/verify.py` recomputes the physics and checks
the deposit against it. **It is not a claim that everything re-runs here.** What a clean clone reaches —
and the parts it does not, because they need `pyscf`, an IBM Quantum account, or a figure source that is
not deposited — is tabulated, measured rather than asserted, at the end of
[`docs/REPRODUCE.md`](docs/REPRODUCE.md), and the gaps are enumerated in
[`docs/KNOWN_DISCREPANCIES.md`](docs/KNOWN_DISCREPANCIES.md). The statement that may be published about
all this is [`docs/DATA_AVAILABILITY.md`](docs/DATA_AVAILABILITY.md).

> **Companion works.** Review — *Machine learning for sample-based quantum diagonalization: generative
> configuration recovery and the classical-simulability frontier* ([arXiv:2608.05314](https://arxiv.org/abs/2608.05314)).
> Method — *Self-falsifying quantum spectroscopy: a transportable necessary-condition screen for
> quantum-computed dynamical spectra* (companion, arXiv posting in progress).

---

## Figure gallery

| Momentum-resolved `A(k,ω)` (exact Mott map) | The resource thesis (74 exact states) |
|:---:|:---:|
| ![akw](docs/img/akw.png) | ![resource](docs/img/resource.png) |
| **Charge structure factor `S(q,ω)` — gapped** | **Spin structure factor `S^zz(q,ω)` — gapless** |
| ![sqw](docs/img/sqw.png) | ![spinqw](docs/img/spinqw.png) |
| **Scaling: fraction falls, absolute cost grows** | **IBM Heron hardware (two real runs)** |
| ![scaling](docs/img/scaling.png) | ![hardware](docs/img/hardware.png) |
| **Nineteen-molecule magic gallery** | **N₂ dissociation hero (`F₁ = 2Nᵤ`)** |
| ![gallery](docs/img/gallery.png) | ![hero](docs/img/hero.png) |

<div align="center"><em>The two collective channels come from the <strong>same sampling primitive</strong>: the
<strong>spin–charge separation</strong> of the 1D Mott insulator — charge gapped (Δ≈5t), spin gapless (πJ/2≈0.79t).</em></div>

---

## What the paper shows

1. **A new use of the sampled-subspace toolkit** — frequency-resolved dynamical response from
   bitstring-sampled subspaces (charged `(N±1)` single-particle spectra and neutral-`N` structure
   factors), with **no Hadamard tests or controlled unitaries**. Validated against exact diagonalization
   on Hubbard chains and across a **nineteen-molecule suite** (energies FCI-verified to `<10⁻⁵` Ha), and
   demonstrated on the **IBM Heron** processor.
2. **The resource that sets the cost — and the one that does not.** The fermionic AntiFlatness collapses
   to a single 1-RDM invariant `F₁ = 4 tr[γ(1−γ)] = 2Nᵤ`. Being an orbital-rotation (Gaussian) invariant,
   it is **provably decoupled** from the basis-dependent support `|S|`; the cost is instead lower-bounded
   and empirically tracked by the bond dimension `χ` (Spearman `ρ = 0.90`). No quantum-advantage claim is
   made — the open question is relocated, not settled.

**Headline hardware number:** stretched N₂ on `ibm_marrakesh` reaches **0.59 ± 0.14 mHa** of exact FCI via
noise-assisted configuration recovery (the same shallow circuit gives 29.5 mHa in noiseless simulation).

---

## Repository structure

```
.
├── src/                       # 55 computation & figure scripts (run from repo root as `python src/X.py`)
│   ├── verify.py              #   fast numpy/scipy-only reproducibility check (what CI runs — see the CI note below)
│   ├── akw_lanczos.py …       #   exact-diagonalization engines + spectral/resource/scaling compute
│   └── make_*.py              #   one data-driven figure generator per figure (no hand-typed numbers)
├── data/                      # authoritative data: every plotted number, FCI-verified (*.json)
├── notebooks/
│   ├── 00_Reproduce_Everything.ipynb   #   ★ MASTER notebook: reproduces every figure & number, end to end
│   ├── HW_LUCJ_N2_Heron_READY.ipynb    #   IBM Heron N₂ energy   (token SCRUBBED — see below)
│   └── Spectral_Heron.ipynb            #   IBM Heron A(ω)        (token SCRUBBED — see below)
├── paper/                     # arXiv-ready LaTeX source + compiled PDF
│   ├── main.tex               #   master file (\input's the sections below)
│   ├── introduction.tex … conclusion.tex, theorem_b2.tex, bibliography.tex
│   ├── figs/                  #   figure fragments (.tex), vector PDFs, plotted data (.dat)
│   ├── main.pdf               #   the compiled preprint (40 pp; arXiv v2 was 39 pp)
│   └── arxiv-submission.tar.gz#   ready-to-upload source bundle
├── docs/
│   ├── REPRODUCE.md           #   figure/number → script → exact command
│   ├── FILE_INDEX.md          #   every file, described
│   ├── FIGURE_PROVENANCE.md   #   figure → data-source → job-id single source of truth
│   ├── KNOWN_DISCREPANCIES.md #   where running a script does NOT reproduce the deposit (read this)
│   └── img/                   #   README thumbnails
├── requirements.txt           # Python dependencies
├── Makefile                   # `make verify`, `make figures`, `make reproduce`, `make paper`
├── CITATION.cff               # citation metadata
└── LICENSE                    # MIT (code) + CC-BY-4.0 (paper text/figures)
```

---

## Quick start

```bash
# 1. environment
python -m venv .venv && source .venv/bin/activate      # (Windows: .venv\Scripts\activate)
pip install -r requirements.txt

# 2. run the adversarial guardian (~15 s, numpy/scipy only). It RECOMPUTES the physics
#    (Hubbard F₁/E₀ for both boundary conditions, the geminal witness) and then checks every
#    deposited artefact — figures, tables, .dat files and the abstract — against it.
python src/verify.py        # or: make verify   ->   PASS, plus a named list of open defects
#    It exits non-zero and names the defect if anything disagrees. The `[XFAIL]` lines it
#    prints are REAL open defects and claims, deliberately visible; they are the v3 agenda.

# 2b. do the committed figures still come out of the committed generators? (writes nothing)
python src/check_figures.py # or: make check-figures

# 3. reproduce EVERYTHING in one coherent, narrated pass — every figure and number, top to bottom
jupyter notebook notebooks/00_Reproduce_Everything.ipynb    # or headless: make reproduce

# 4. regenerate a single data-driven figure fragment from its committed data
python src/make_decoupling_native.py    # -> paper/figs/fig_decoupling_native.tex

# 5. build the paper (needs a TeX distribution, e.g. TeX Live / MiKTeX)
make paper                              # -> paper/main.pdf
```

**Posting to arXiv?** The upload-ready bundle `paper/arxiv-submission.tar.gz` compiles clean-room with
`pdflatex` alone (39 pp, 0 undefined refs) — see **[`docs/ARXIV_SUBMISSION.md`](docs/ARXIV_SUBMISSION.md)**
for the verification, metadata, and step-by-step. **Note:** that bundle is a frozen snapshot of what was
posted as **v2**; nine files in `paper/` have since moved ahead of it and two more have been deleted
from the tree. Run `python src/check_tarball.py` for the live comparison. **The bundle still contains the
two fabricated rows of `figs/sqw_edges.dat`**, so it must be rebuilt before any new submission — see
[`docs/KNOWN_DISCREPANCIES.md`](docs/KNOWN_DISCREPANCIES.md) §7.

The **master notebook** [`notebooks/00_Reproduce_Everything.ipynb`](notebooks/00_Reproduce_Everything.ipynb)
walks the whole pipeline end to end — the decoupling identity and witness, the sampled spectral functions,
the χ–\|S\| resource map, the scaling, the 19-molecule magic suite, the selector/noise/AI results, the
hardware run, and the figures — each number loaded from the same `data/*.json` that feeds the paper.
Every figure's one-line reproduce command is in **[`docs/REPRODUCE.md`](docs/REPRODUCE.md)**.

---

## Reproducibility at a glance

| Layer | Reproducible here? | How |
|---|---|---|
| **The guardian** | ✅ **verified 2026-09-18** — PASS, 801 checks, 5 014 numeric assertions, 11 registered xfail | `python src/verify.py` (or `make verify`) |
| **Everything, one coherent pass** | ⚠️ **not verified in the 2026-09-18 pass** — the notebook needs `pyscf`, which was not installed in the environment these checks were run in | `notebooks/00_Reproduce_Everything.ipynb` (or `make reproduce`) |
| **Exact diagonalization, numpy/scipy only** (`A(ω)`, `A(k,ω)`, `S(q,ω)`, `S^zz`, the χ–\|S\| resource map, the geminal witness, the scaling) | ✅ **verified locally**, minutes | `python src/<script>.py`; see `docs/REPRODUCE.md` |
| **The 19-molecule magic suite** | ⚠️ runs, but needs `pyscf`; **untested in the 2026-09-18 pass** | `python src/n19_suite.py` etc. |
| **Matrix-free scaling** to 28 qubits | ✅ (RAM-staged) | `src/scaling_lanczos_mf.py` (resumable, checkpointed) |
| **Figures** | ✅ the ten native fragments; ⚠️ the ten PDF figures need sources kept outside this deposit | `make figures` (native pgfplots from the `data/*.json`/`.dat`); inventory in [`docs/FILE_INDEX.md`](docs/FILE_INDEX.md) |
| **IBM Heron hardware runs** | ⚠️ needs an IBM Quantum account | `notebooks/` (cached counts + a re-run cell; tokens scrubbed) |

> **Where the output lands (changed 2026-09-18).** Every compute script now writes into *this clone's*
> `data/` by default — the path is derived from the script's own location, so the working directory does
> not matter — and `--out PATH` overrides it. Until that date **eighteen of them wrote to `/w/`**, the
> working directory of the container the original runs were made in, so running the documented pipeline
> did **not** repopulate `data/`. That is fixed and the evidence is in
> [`docs/KNOWN_DISCREPANCIES.md`](docs/KNOWN_DISCREPANCIES.md) §2.
>
> **CI note.** `.github/workflows/ci.yml` was listed in `.gitignore` until 2026-09-18, so the workflow
> was **not in the repository and had never run**. It is tracked now and both of its commands pass
> locally, but it will execute for the first time only when the branch carrying it is pushed. §15.

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
@article{BonillaVargas2026DynamicalSpectral,
  title   = {Dynamical spectral functions from bitstring-sampled quantum subspaces:
             entanglement, not one-body magic, tracks the sampling cost},
  author  = {Bonilla Vargas, Nicol\'as},
  journal = {arXiv preprint arXiv:2608.16436},
  year    = {2026}
}
```

## License

- **Code** (`*.py`, `notebooks/`, `Makefile`): [MIT](LICENSE).
- **Paper text and figures** (`paper/`): [CC-BY-4.0](LICENSE).
