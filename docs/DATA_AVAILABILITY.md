# Data Availability Statement — for Physical Review A

*Rewritten 2026-09-18, second pass, after an adversarial audit found **three false sentences in Part 1**
— the only part that reaches a reader. Every claim below was re-checked against the deposit on that date
by executing something, not by reading. Where the deposit does not support a claim, the claim is not
made; where a previous draft made it anyway, §3 records what it said and why it was wrong.*

Physical Review A requires that *"all published articles must include a Data Availability Statement (DAS)
that explains whether and how authors are sharing their data."* A DAS is a statement of fact about an
artefact, in the same class as a figure caption: **if it is not true of the deposit, it is a false
statement in a published paper.** This file therefore has three parts — the statement itself, the audit
that licenses each of its sentences, and the record of the sentences that were removed.

---

## Part 1 — the statement to paste into the manuscript

> **Data availability.** The data that underlie the figures and the quoted numbers of this work are
> openly available in the reproduction repository at
> `https://github.com/nicolasbonilla/dynamical-spectral-functions-sqd`, with the code under the MIT
> licence and the data and manuscript under CC BY 4.0. Each figure's plotted values are deposited as a
> JSON or plain-text table, and `docs/FIGURE_PROVENANCE.md` maps every figure to the file that is
> authoritative for it; the two exceptions — a schematic that plots no data, and one panel whose four
> free-fermion support sizes are written into the figure source itself — are named there.
>
> The classical data — exact-diagonalization spectral functions of the one-dimensional Hubbard model, the
> molecular suite, and the resource statistics — are additionally *regenerable*: the computation scripts
> that produced them are deposited in `src/`, each resolves its output path from its own location on disk
> and writes into the repository's own `data/` directory, and `docs/REPRODUCE.md` gives the command for
> each figure together with a measured table of what a clean clone does and does not reach. An automated
> guardian (`python src/verify.py`) recomputes the underlying physics by exact diagonalization and checks
> the deposited artefacts against it: 801 checks over 5 014 numeric assertions, exiting non-zero on any
> discrepancy. It covers the deposited data files and the numbers printed in the figure sources and the
> abstract; it is not a proof that every sentence of the manuscript is verified, and eleven known open
> defects and claims are printed by name on every run.
>
> The quantum-hardware results are deposited as the arrays that enter the figures, together with the
> backend name, the shot count and the IBM Quantum job identifier for each of the two runs
> (`data/heron_spectral.json`, `data/hw_lucj_n2_result.json`). **The raw measured bitstrings are not part
> of the deposit**: the notebooks in `notebooks/` retrieve the counts from the IBM Quantum service by job
> identifier, which requires an account and a job that is still retrievable, so the classical
> post-processing of the raw counts cannot be re-run from the deposited files alone. The hardware runs
> themselves cannot be regenerated in any case: the jobs are closed, and a device executed today has a
> different calibration, so new counts would be new data rather than a reproduction. For the
> generative-model control of Fig. 18 only the five summary rows are deposited, in `data/gflow.json`; no
> script in the repository regenerates them, and the per-seed raw traces are available from the authors
> on reasonable request.
>
> Known gaps between the deposited code and the deposited artefacts — including ten figures that enter as
> PDF without their standalone LaTeX source, three computation scripts whose output is not deposited, and
> the parts of the pipeline that require `pyscf` or an IBM Quantum account — are enumerated in
> `docs/KNOWN_DISCREPANCIES.md` rather than left for the reader to discover.

**Before this is signed, two things must be true and are not automatically true:**

1. **The repository URL must be live and public.** It is the URL already printed in `paper/hardware.tex`
   and it is now also in `README.md` and `CITATION.cff` — but *this deposit cannot verify that the remote
   exists.* Open it in a browser from a signed-out session before submitting.
2. **The deposit must be committed.** Everything described above is true of the working tree. Until it is
   committed and pushed, `git clone` delivers an earlier tree in which several of these sentences are
   false. **Signing the DAS while the repairs are uncommitted publishes a false statement.**

---

## Part 2 — the audit that licenses each sentence

### 2.1 "openly available … MIT / CC BY 4.0" — TRUE, after the licence was made explicit

`LICENSE` is dual: MIT for code, CC BY 4.0 for `paper/`. Until 2026-09-18 **`data/` was covered by
neither heading**, so the earlier draft's "the data … under CC BY 4.0" claimed a licence the file did not
grant. `LICENSE` now names `data/` explicitly under the CC BY 4.0 heading, and Part 1 states the split
instead of flattening it to one licence. `CITATION.cff` already declared `MIT AND CC-BY-4.0`.

### 2.2 "Each figure's plotted values are deposited" — TRUE, with two exceptions now named in Part 1

`docs/FIGURE_PROVENANCE.md` carries one row per figure naming the authoritative file. The two things that
sentence must not be allowed to cover, and which Part 1 therefore names:

* **Fig. 14** is a circuit schematic; it plots no data.
* **Panel (c) of Fig. 1** carries four free-fermion support sizes (35, 336, 3496, 37361). They are
  *constants typed into the generator* — `src/make_decoupling_native.py:34` — and from there into
  `paper/figs/fig_decoupling_native.tex`; **no file in `data/` holds them**, and no deposited script
  recomputes them. `docs/FIGURE_PROVENANCE.md` row 1 has always said "hard-typed"; the point here is
  that the DAS must not let "every plotted value is deposited" quietly cover them. The same four numbers
  are printed in `paper/resource.tex`.

This claim is deliberately weaker than "every figure can be recompiled": see 2.6.

### 2.3 "the computation scripts … each … writes into the repository's own `data/`" — TRUE AS OF 2026-09-18, IN TWO PASSES

**First pass (morning).** Eighteen scripts wrote to `/w/`, the working directory of the Docker container
the original runs were made in. Running the documented pipeline did not repopulate `data/`; the deposited
files had been moved there by hand. **A DAS signed on that tree would have been false in its central
sentence.**

**Second pass (this one).** The audit found that the repair had missed four computation scripts —
`amortized_recovery.py`, `ladder_vs_chain.py`, `sampled_akw.py`, `spin_lanczos.py` — which were broken in
a quieter way: they wrote to the *relative* path `data/…`, which lands in whatever directory the reader
happens to be standing in and raises `FileNotFoundError` from anywhere but the repository root. The word
*"each"* in Part 1 was therefore false for four files. All four now carry the same `DEPOSIT PATHS` block
as the other eighteen. The fifteen figure generators had the same defect and are anchored too.

Verified three ways, all today, all from a directory outside the repository: the path block of each
repaired script executed in isolation and its destination printed (**every destination absolute and
inside the repository**); `python src/sqw_lanczos.py 6 --eta 0.18 --out …` run from `C:\Users\…` to
completion, exit 0; `python src/ladder_vs_chain.py --out …` run to completion from the scratch directory,
exit 0, output written where `--out` asked and nowhere else.

### 2.4 The guardian sentence — TRUE, and it now states its size and its limits

`python src/verify.py`: **801 checks, 0 failures, 11 registered xfail, 2 skips, 5 014 numeric assertions,
12.3 s, exit 0**, measured on the final tree. An earlier draft quoted 4 951 assertions, which no longer
reproduced; quoting a stale count in a statement about verification is its own small defect.

The audit also asked the right question — *does the guardian actually bite?* Two things were done:

* An adversarial pass showed that **`data/charge_gap.json` had no coverage at all**: falsifying
  Δ(L=12) from 4.968759 to 4.900000 still gave PASS, 0 FAIL. That file backs the Δ = 4.97 t of the
  manuscript. Section 9.14 of the guardian now covers it — each row's μ⁺, μ⁻ and Δ rebuilt from its own
  three sector energies, particle–hole symmetry μ⁺+μ⁻ = U, the three sector dimensions against their
  binomials, the ground-state energies against the *same* energies deposited by three other scripts in
  three other files, and the Lieb–Wu limit recomputed here by quadrature (4.679517107464081 against the
  deposited 4.679517107464084). Both falsifications were re-run against the new code: the first now
  fails on "Δ does not follow from the chemical potentials in its own row"; a subtler one that moved
  E₀(L=10) by 10⁻⁶ *and kept the row internally consistent* fails on the cross-file check against
  `data/scaling_data.json`.
* Part 1 states the guardian's **scope** in its own sentence rather than letting "an automated guardian
  checks the deposited artefacts" be read as "the manuscript is verified". It is not: the eleven xfail
  lines include three sentences the manuscript still prints and cannot defend.

### 2.5 The hardware sentence — REWRITTEN; the previous version was false twice

The previous draft said *"the raw bitstring histograms … are included in `data/heron_spectral.json` and
`data/hw_lucj_n2_result.json`"* and *"all classical post-processing of those counts is reproducible from
the deposited files without further device access."* Both are false, and a referee would find the first
in under a minute by opening either file:

| file | what is actually in it |
|---|---|
| `data/heron_spectral.json` | `L, U, eta, nsector, grid[600], A_exact[600], A_hw[600], A_local[600], backend, shots, job_id, hw_relL1, hw_S` — **post-processed spectral arrays**. No bitstring, no count, no list of recovered determinants. |
| `data/hw_lucj_n2_result.json` | energies. The key that reads like a histogram, `hist_hw`, is **five SQD iteration energies** (−108.781…, −108.801…, −108.806…, −108.807…), not a measurement histogram. |

A scan of every `data/*.json` for a key containing *counts*, *bit* or *hist* in the measurement sense
returns nothing. And the notebooks are explicit about where the counts come from:
`HW_LUCJ_N2_Heron_READY.ipynb` runs `result_hw = job.result()[0]` under the comment
`job = service.job('PASTE_JOB_ID_HERE')`; `Spectral_Heron.ipynb` runs `res = job.result()` and then
`get_counts()`. The counts arrive over the network, from an account. Part 1 now says exactly that.

What *is* true and is claimed: the arrays that enter the figures are deposited as measured, with backend,
shots and job id, so every number plotted from the hardware run can be checked against the file it came
from. `docs/REPRODUCE.md` carried the same false sentence ("the notebooks ship with cached device
counts") and has been corrected.

One further honesty item, repeated here because a DAS is where a reader looks for it: the deposited
`heron_spectral.json` records `hw_S = 300` out of `nsector = 300`. The recovered subspace is the **whole**
sector, so `A_hw` matches `A_exact` **by coverage, not by fidelity** — and it matches to machine
precision (max|diff| = 7.9 × 10⁻¹⁵), not byte-for-byte as the file's own provenance string used to say.
That string has been corrected in the deposited file.

### 2.6 What the DAS does **not** claim, and why

| not claimed | why |
|---|---|
| that every figure can be **recompiled** from the deposit | ten of the twenty figures ship as PDF with no `.tex` source inside the repository (their standalone sources live in the author's working tree), and `fig_molecular_gallery.pdf` additionally needs 19 trimmed PyMOL orbital PNGs that are not deposited. `KNOWN_DISCREPANCIES.md` §5. |
| that `make reproduce` runs | it executes a notebook that requires `pyscf`, which is not installed in the verified environment. It was **not** run in this pass and no claim is made for it. |
| that the GFlowNet comparison is reproducible | `data/gflow.json` holds five summary rows; **no script in `src/` generates it.** Part 1 now says both halves of this — no generator *and* raw on request — rather than only the second. |
| that the hardware measurement can be re-obtained | see 2.5. |
| that the classical post-processing of the raw counts runs offline | see 2.5. It does not; the counts are not deposited. |
| that the guardian runs in continuous integration | `.github/workflows/ci.yml` exists and is valid, and `.gitignore` no longer excludes it, but **it has never executed**: it cannot until the commit containing it is pushed. The sentence *"The guardian is executed on every commit by continuous integration"* may be added to Part 1 **only after** a run has actually completed. |
| that `paper/arxiv-submission.tar.gz` is a source of data | it is the frozen v2 bundle and **still contains two fabricated data rows** and the superseded `19184`/`824504`; see the banner at the top of `KNOWN_DISCREPANCIES.md`. Part 1 no longer points at the arXiv source package at all — see §3 item 4. It is kept as a historical record and **must never be uploaded as v3**. |

### 2.7 The hardware runs, in one table

| what | where | backend | job id |
|---|---|---|---|
| N₂ energy (LUCJ + recovery, 8 seeds) | `data/hw_lucj_n2_result.json` | `ibm_marrakesh` | `da125f2ein7c73bcsqs0` |
| spectral `A(ω)` | `data/heron_spectral.json` | `ibm_fez` | `d9s16avpemts73ct6g8g` |

### 2.8 Three deposited scripts whose output is not deposited

`src/gate1_ladder.py`, `src/rigor_floor.py` and `src/recolor_akw_v2.py`. Detail in
`KNOWN_DISCREPANCIES.md` §13. They are kept because a deposited script with a missing output is a gap a
reader can see, whereas a deleted script is a gap nobody can see.

---

## Part 3 — the sentences that were removed, and why

Kept in full, because a DAS that quietly improves is a DAS nobody can audit.

1. **"the raw bitstring histograms … are included in `data/heron_spectral.json` and
   `data/hw_lucj_n2_result.json`"** — FALSE. Neither file contains a bitstring or a count (2.5).
2. **"all classical post-processing of those counts is reproducible from the deposited files without
   further device access"** — FALSE, and the more damaging of the two, because it is exactly the kind of
   promise PRA and SciPost check. The counts come from `job.result()` over the network (2.5).
3. **"the computation scripts … each writes into the repository's own `data/` directory"** — FALSE for
   four scripts, which wrote relative to the caller's working directory (2.3). Repaired, so the sentence
   is now kept rather than removed.
4. **"and in the source package of arXiv:2608.16436"** — REMOVED as a *data source*. The v2 source
   package is exactly where the two fabricated rows of `sqw_edges.dat` and the superseded 19184/824504
   still live. Pointing a reader there for data, with the warning buried in Part 2, would send them to
   the one copy that is wrong. The repository is the deposit; the e-print is the paper.
5. **"under CC BY 4.0"** applied to everything — REPLACED by the actual split (2.1).
6. **"`LICENSE` is CC BY 4.0"** (Part 2 of the previous draft) — it is dual, and `data/` was in neither
   half. Fixed in `LICENSE`, not just in prose.
7. **"4 951 numeric assertions, 12.9 s"** — stale; the measured values are 5 014 and 12.3 s (2.4).
8. **"(Replace `<OWNER>/<REPO>` … it is already in `README.md` and `CITATION.cff`)"** — it was in
   neither. The URL printed in `paper/hardware.tex` has been added to both, and Part 1 now carries the
   real URL with an explicit instruction to confirm the remote is live.

---

## Part 4 — what to do before submitting

1. **Commit and push.** Nothing above is true of a clone until it is. This is the single largest gap
   between what this document says and what a referee would download.
2. Open `https://github.com/nicolasbonilla/dynamical-spectral-functions-sqd` signed out and confirm it
   resolves. If the repository has a different name, change it in Part 1, `README.md`, `CITATION.cff`
   **and** `paper/hardware.tex`.
3. Confirm the workflow in `.github/workflows/ci.yml` has actually run at least once. Only then may the
   continuous-integration sentence of 2.6 be added.
4. Decide about the GFlowNet per-seed traces: deposit them, or keep "available on request". Either is
   acceptable to PRA; silence is not.
5. Re-run `python src/verify.py` against the final tree and require exit 0. If the numbers move, update
   the count in Part 1 — it is a factual claim like any other.
6. Build the v3 source package from the working tree. **Do not upload `paper/arxiv-submission.tar.gz`.**

---

## Appendix — the same statement for the other candidate venues

**SciPost Physics Core** publishes its acceptance criteria and requires that data and code be made
available; its report is public and permanent, so a DAS that over-promises stays on the record next to
the paper forever. The Part 1 text satisfies it as written, provided the repository is live at
submission time.

**IOP (Quantum Science and Technology)** asks for a data-availability statement in a fixed set of forms.
The matching form is:

> "The data that support the findings of this study are openly available at
> `https://github.com/nicolasbonilla/dynamical-spectral-functions-sqd`. The quantum-hardware results are
> included as the deposited arrays with their backend names, shot counts and job identifiers; the raw
> measured bitstrings are not deposited and the hardware runs cannot be re-executed, for the reasons
> given in the repository's `docs/DATA_AVAILABILITY.md`."
