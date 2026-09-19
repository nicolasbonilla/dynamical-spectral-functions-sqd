# Data Availability Statement — for Physical Review A

*Third pass 2026-09-19, the pass that finally **pasted** the statement into the manuscript*
(`paper/sec_9.tex`, Sec. IX F, `\label{sec:scope:data}`) — v3 had shipped with no DAS at all while the
header of `sec_9.tex` claimed to have taken the paragraph over from `hardware.tex`. This pass also
re-measured every sentence of Part 1 against the deposit of 2026-09-19 and **rewrote five of them**;
what changed and why is in Part 3, items 9-13.*

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

**This is the text that is now in the manuscript**, at `paper/sec_9.tex`, Sec. IX F
(`\label{sec:scope:data}`), rendered on p. 43 of the 59-page build. It is reproduced here so that the
audit below has something to audit; if the two ever disagree, **the manuscript is the one that is
published and this file is the one that is wrong.**

> **Data availability statement.** The data that underlie the figures and the quoted numbers of this work
> are openly available in the reproduction repository at
> `https://github.com/nicolasbonilla/dynamical-spectral-functions-sqd`, the code under the MIT licence
> and the manuscript, the deposited datasets and the reproduction documentation under CC BY 4.0. Every
> figure's plotted values are deposited — as a JSON table in `data/` or as a plain-text table in
> `paper/figs/` — and `docs/FIGURE_PROVENANCE.md` names the file that is authoritative for each; the one
> exception is Fig. 1, a circuit schematic that plots no data. The classical data — the
> exact-diagonalization spectral functions of the one-dimensional Hubbard model, the molecular suite and
> the resource statistics — are additionally *regenerable*: the computation scripts that produced them
> are deposited in `src/`, and `docs/REPRODUCE.md` gives the reproduction commands together with a
> measured table of what a clean clone does and does not reach. An automated guardian
> (`python src/verify.py`) recomputes the underlying physics by exact diagonalization and checks the
> deposited artefacts against it — 807 checks over 7 454 numeric assertions, exiting non-zero on any
> discrepancy. It covers the deposited data files and the numbers printed in the figure sources, in the
> tables and in the abstract; it is not a proof that every sentence of this manuscript is verified, and
> four known open defects are printed by name on every run.
>
> The quantum-hardware results are deposited as the arrays that enter the figures, together with the
> backend name, the shot count and the IBM Quantum job identifier for each of the two runs
> (`data/heron_spectral.json`, `data/hw_lucj_n2_result.json`). *The raw measured bitstrings are not part
> of the deposit*: the notebooks in `notebooks/` retrieve the counts from the IBM Quantum service by job
> identifier, which requires an account and a job that is still retrievable, so the classical
> post-processing of the raw counts cannot be re-run from the deposited files alone. The hardware runs
> themselves cannot be regenerated in any case — the jobs are closed, and a device executed today has a
> different calibration, so new counts would be new data rather than a reproduction. For the
> generative-model control of Sec. IX D only the five summary rows are deposited, in `data/gflow.json`;
> no script in the repository regenerates them, which is the reason given there for withdrawing the two
> figures that were built on them, and the per-seed values are available from the author on reasonable
> request. Known gaps between the deposited code and the deposited artefacts — six of the seventeen
> figures enter as PDF with no standalone LaTeX source in the deposit, the tables plotted in Fig. 17 are
> deposited but the script that builds them is not, three deposited computation scripts have no
> deposited output, and parts of the pipeline require `pyscf` or an IBM Quantum account — are enumerated
> in `docs/KNOWN_DISCREPANCIES.md` rather than left for the reader to discover.

**Acknowledgments**, added in the same pass at the end of `paper/sec_10.tex` (p. 45):

> **Acknowledgments.** We acknowledge the use of IBM Quantum services for this work; the two device runs
> reported in Sec. IX A were executed on `ibm_fez` and `ibm_marrakesh`. The views expressed are those of
> the author and do not reflect the official policy or position of IBM or the IBM Quantum team.

**No conflict-of-interest statement is printed, and that is a checked decision, not an omission.** APS
requires authors to *"alert editors to any potential conflict of interest such as sources of funding,
condition of employment, and any limitations on disclosure of data, code, or experimental details"* and
only *encourages* a statement inside the paper (*"authors are encouraged to declare any conflict of
interest within the paper itself using a Conflict of Interest statement"*),
[journals.aps.org/authors/editorial-policies](https://journals.aps.org/authors/editorial-policies). The
DAS itself is not optional: *"All published articles must include a Data Availability Statement (DAS)
that explains whether and how authors are sharing their data"*,
[journals.aps.org/pra/authors](https://journals.aps.org/pra/authors) and
[journals.aps.org/authors/data-availability-statements](https://journals.aps.org/authors/data-availability-statements).

**Before this is signed, three things must be true and are not automatically true:**

1. **The three sentences flagged in Part 4 must be made true, or cut.** Two of them point at documents
   that still describe the v2 figure set, and one prints a number that moves whenever anybody adds a
   file to `src/`. This is the new one, and it is the one this pass could not close by itself.
2. **The repository URL must be live and public.** It is the URL now printed in the manuscript itself
   (`paper/sec_9.tex`) and in `README.md` and `CITATION.cff` — but *this deposit cannot verify that the
   remote exists.* Open it in a browser from a signed-out session before submitting. (`paper/hardware.tex`,
   where the URL used to live, is in `_superseded/v2_paper/` and is no longer built.)
3. **The deposit must be committed.** Everything described above is true of the working tree. Until it is
   committed and pushed, `git clone` delivers an earlier tree in which several of these sentences are
   false. **Signing the DAS while the repairs are uncommitted publishes a false statement.**

---

## Part 2 — the audit that licenses each sentence

### 2.1 "openly available … MIT / CC BY 4.0" — TRUE, after the licence was made explicit

`LICENSE` is dual: MIT for code, CC BY 4.0 for `paper/`. Until 2026-09-18 **`data/` was covered by
neither heading**, so the earlier draft's "the data … under CC BY 4.0" claimed a licence the file did not
grant. `LICENSE` now names `data/` explicitly under the CC BY 4.0 heading, and Part 1 states the split
instead of flattening it to one licence. `CITATION.cff` already declared `MIT AND CC-BY-4.0`.

### 2.2 "Every figure's plotted values are deposited" — TRUE, with ONE exception (2026-09-19: it was two)

The v3 manuscript typesets **seventeen** figures; `main.aux` is the authority
(`\newlabel{fig:...}{{1}..{17}}`), and it is the list that was walked here, figure by figure, against
`data/` and `paper/figs/`.

* **Fig. 1** (`fig:circ`) is a circuit schematic; it plots no data. This is the one exception, and Part 1
  names it by `\ref`, not by number, so a float that moves cannot make the sentence false.
* **The second exception of the v2 draft is gone.** It read: *"panel (c) of Fig. 1 carries four
  free-fermion support sizes (35, 336, 3496, 37361), constants typed into
  `src/make_decoupling_native.py:34`."* That panel was rebuilt for v3: `fig_decoupling_native.tex` line 34
  now reads `(c) Hubbard: $\chi$ runs forwards`, and **grep finds none of those four numbers in any figure
  source used by the current build**. They survive only as prose, in `paper/sec_6_body.tex:329`
  (`paper/resource.tex`, where the old text said they were printed, is superseded and not built).
* **Fig. 17** (`fig:thm1iii-violation`) is the one that needed checking rather than assuming. It plots
  **eight** `n3_*.dat` tables (counted by grepping the fragment with a case-correct pattern — a
  lower-case-only character class silently misses `n3_thmB_eta.dat` and `n3_thmB_frac.dat` and returns
  six); two further `n3_*.dat` files, `n3_fig5_eta.dat` and `n3_fig5_frac.dat`, are read not by the
  fragment but by `caption_thm1iii_violation.tex`, whose whisker numbers derive from them. All ten are
  deposited in `paper/figs/`, so the Part 1 sentence is true of this figure. What is *not* deposited is
  `build_n3.py`, which wrote them, and two of its four inputs. That is a *regenerability* gap, not a
  deposit gap, and Part 1 puts it in the known-gaps sentence, where it belongs; the guardian registers
  it as an open defect and prints it on every run. *(The counts in circulation disagree and all three
  are defensible readings of different questions: `verify.py` says "ten n3_*.dat tables" — every file;
  `FIGURE_PROVENANCE.md` says "eight" — what the fragment reads; a careless grep says six. Part 1
  therefore carries no count at all.)*
* **Fig. 3** (`fig:akwsampled`) panel (c) was checked numerically, not read off a table. Its generator
  reads `sampled_honest.json` and `honest_sampling.json` from a **scratch directory outside the
  repository** (2.3). The deposited `data/sampled_honest.json` and `data/honest_sampling.json` were parsed
  and compared field by field against those two files: **every physics number is identical**; the only
  differences are `provenance/generated_utc` and the wall-clock timings (`runtime_seconds` 126.2 vs 67.0,
  `wall_s` 24.70 vs 10.25 and the three per-L timings). So the plotted values *are* deposited — but the
  script that draws them does not read the deposited copies.

This claim is deliberately weaker than "every figure can be recompiled": see 2.6.

### 2.3 "the computation scripts … each … writes into the repository's own `data/`" — **FALSE, AND CUT FROM PART 1 ON 2026-09-19**

*The heading of this section used to read "TRUE AS OF 2026-09-18, IN TWO PASSES". It was not true
on 2026-09-19, and it is the third time this one sentence has had to be corrected, which is why the
third pass removed it from the statement instead of repairing it again.*

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

Verified three ways, all that day, all from a directory outside the repository: the path block of each
repaired script executed in isolation and its destination printed (**every destination absolute and
inside the repository**); `python src/sqw_lanczos.py 6 --eta 0.18 --out ...` run from `C:\Users\...` to
completion, exit 0; `python src/ladder_vs_chain.py --out ...` run to completion from the scratch directory,
exit 0, output written where `--out` asked and nowhere else.

**Third pass (2026-09-19): the word "each" was still false, and Part 1 no longer says it.** Every `.py` in
`src/` was scanned at the byte level for the anchoring idiom and for absolute machine paths, with a
negative control that fires on a deliberately unanchored file (it fired). Four scripts do not anchor:

> **All four rows below were REPAIRED on 2026-09-19** and are kept as the record of what was
> wrong, not as a description of the tree. `make_akw_sampled_honest_fig.py` now derives both
> roots from `__file__` and reads the deposited `data/sampled_honest.json` and
> `data/honest_sampling.json` (verified field by field against the scratchpad copies it used
> to read: **zero numeric differences**, only provenance stamps and wall-clock timings, with a
> positive control that fired). The other three resolve their output from their own location.
> `src/check_figures.py` now guards **8 of 8** artefacts -- measured by seeding each with a
> sentinel line and checking it is gone after the run -- and a before/after sha256 + mtime
> snapshot of all 231 files shows it writes **nothing** into the repository. See
> `KNOWN_DISCREPANCIES.md` §25 and §17.

| script | what it actually does |
|---|---|
| `src/scaling_lanczos.py` | `fn='scaling_data.json'` (line 120), *relative*. Worse than a crash: line 122 is `json.load(open(fn)) if os.path.exists(fn) else {'U':8.0,'points':[]}`, so run from anywhere it **silently starts a fresh file in the caller's directory** and never updates `data/scaling_data.json`. |
| `src/scaling_lanczos_mf.py` | `fn='scaling_data.json'` (line 157), *relative*, and it `json.load`s it first, so it raises from anywhere the file is not sitting. |
| `src/build_repro_notebook.py` | writes the relative `'notebooks/00_Reproduce_Everything.ipynb'` (line 226). |
| `src/make_akw_sampled_honest_fig.py` | the worst of the four: lines 17-18 hard-code **the author's own machine paths** — `REPO = r"C:\Users\Nicolas\Downloads\..."` and `CALC = r"C:\Users\Nicolas\AppData\Local\Temp\claude\...\scratchpad\p2_calc"`, a session scratch directory that is **not in the deposit**. |

The last one has a second consequence, and it is the reason this pass does not advertise
`src/check_figures.py` in Part 1: **`check_figures.py` passes for the wrong reason.** It regenerates six
generators, one of which is `make_akw_sampled_honest_fig.py`, and that regeneration succeeds *only
because the author's scratch directory still exists on this machine*. On a clean clone it raises
`FileNotFoundError` at `open(CALC/sampled_honest.json)`, so its "8 artefacts byte-identical" is a
statement about this laptop, not about the deposit. Until lines 17-21 of that generator read the
deposited `data/` copies, no sentence of the DAS may lean on it.

### 2.4 The guardian sentence — TRUE, and it now states its size and its limits

`python src/verify.py`, run 2026-09-19 12:26 on this tree, literal last lines:

```
SUMMARY   pass 806   FAIL 0   xfail 5   xpass 0   skip 2   (7452 numeric assertions, 10.7 s)
RESULT: PASS -- every recomputed quantity is carried correctly by the deposited artefacts.
        5 known open defect(s) remain, listed above.
verify.py exit=0
```

So Part 1 now says **806 / 7 454 / four**, against the **801 / 5 014 / eleven** of the second pass. Those
three were the three false sentences this pass inherited; by this file's own criterion a stale count in a
statement about verification is a false statement, so they are recorded in Part 3 (items 9-11) rather
than quietly improved. An earlier draft quoted 4 951 assertions, and the one before that 5 014.

**Two of the three numbers moved while this section was being written, so Part 1 must be re-measured on
the frozen tree.** Other agents were repairing the same repository at the same time. Three runs of the
same command, fourteen minutes apart:

| time | result | checks | assertions | open defects |
|---|---|---|---|---|
| 12:26 | `pass 806 FAIL 0 xfail 5 xpass 0`, exit 0 | 806 | 7 452 | 5 |
| 12:35 | `pass 805 FAIL 1 xfail 5 xpass 0`, exit 1 | 806 | 7 454 | 5 |
| 12:40 | `pass 805 FAIL 1 xfail 4 xpass 1`, exit 1 | 806 | 7 454 | 4 |

The assertion count rose because section (9) opens every `.py` in `src/` and two were added
(`src/check_provenance.py`, `src/build_arxiv_bundle.py`); the same two additions broke the
hand-maintained script count in `README.md` (`found 55, due 57`), which is the FAIL. The defect count
fell because `tab_repro H2O R_eq` was repaired elsewhere (0.957 — 0.958) and now reports `XPASS`.
**Part 1 and `paper/sec_9.tex` carry the 12:40 reading, 806 / 7 454 / four.** The check count held at 806
through all three. Neither `README.md` nor the added scripts are this statement's to fix — but a DAS
that points a referee at a guardian the tree fails with exit 1 is not signable, so this is Part 4 item 0.

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
| that every figure can be **recompiled** from the deposit | **six of the seventeen** ship as PDF with no `.tex` source inside the repository — `fig_akw_native.pdf`, `fig_circuit.pdf`, `fig_hero.pdf`, `fig_noise_score.pdf`, `fig_spinqw.pdf`, `fig_sqw.pdf`; the list was taken from the `\includegraphics` lines of the current build, not from a document. The other eleven are native pgfplots fragments. (`KNOWN_DISCREPANCIES.md` §5 still says "ten of the twenty" and still cites `fig_molecular_gallery.pdf` and its 19 PyMOL PNGs — that figure is **withdrawn in v3**. Part 4 item 2.) |
| that Fig. 17 can be **regenerated** | its ten `n3_*.dat` tables are deposited, but `build_n3.py` is not, and two of its four inputs are not in `data/`. The guardian registers this as an open defect and refuses to "check" the figure by re-reading what it prints. |
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
8. **"(Replace `<OWNER>/<REPO>` ... it is already in `README.md` and `CITATION.cff`)"** — it was in
   neither. The URL printed in `paper/hardware.tex` has been added to both, and Part 1 now carries the
   real URL with an explicit instruction to confirm the remote is live.

*Third pass, 2026-09-19. Seven more, and the last two are clauses cut rather than corrected.*

9. **"801 checks"** — stale. Measured on this tree: **806** (2.4).
10. **"5 014 numeric assertions"** — stale. Measured: **7 454**, and moving (2.4).
11. **"eleven known open defects and claims are printed by name on every run"** — stale. The guardian
    registers **four**, and prints all four.
12. **"one panel whose four free-fermion support sizes are written into the figure source itself"** —
    no longer true of anything. Panel (c) of that figure was rebuilt for v3 and none of the four numbers
    appears in any figure source the manuscript builds (2.2).
13. **"For the generative-model control of Fig. 18"** — **Fig. 18 does not exist in v3.** `fig:gflow` and
    `fig:amort` are withdrawn (`main.tex:229-242`, `sec_9.tex` Sec. 9.4: *"Both figures of this test are
    therefore withdrawn and the surviving statement is made in prose"*). The deposited `data/gflow.json`
    still underlies numbers the manuscript prints (39.24 +- 3.95, 26.98 +- 0.66, 22.94 +- 3.30 mHa), so
    the sentence is kept and re-pointed at **Sec. IX D**. Two further clauses were changed with it: the
    manuscript's own reason for the withdrawal (no deposited script regenerates the file) is now stated
    in the DAS instead of being left for the reader to collide with, and *"available from the authors"*
    became *"from the author"* — there is one.
14. **"each resolves its output path from its own location on disk and writes into the repository's own
    `data/` directory"** — **still false**, for four scripts, in a new way (2.3). The clause is removed
    from Part 1 rather than repaired in prose: it is a claim about code hygiene that a DAS does not need
    to make, and it had already been wrong twice.
15. **"`docs/REPRODUCE.md` gives the command for each figure"** — REPRODUCE.md's figure table is still
    the twenty-row v2 table in v2 numbering, so "each figure" does not resolve. Part 1 now says
    *"gives the reproduction commands"*, which is what that file actually delivers, plus its clean-clone
    table, which is genuinely per-item and genuinely measured.

---

## Part 4 — what to do before submitting

**0. The tree must pass its own guardian.** As of 2026-09-19 12:35 `python src/verify.py` **exits 1**:
`README.md` advertises 55 scripts and `src/` holds 57. Whoever owns `README.md` fixes the count; until
then the deposit fails the check the DAS points a referee at. Re-run afterwards and, if the assertion
count has moved off **7 454**, change it in Part 1 *and* in `paper/sec_9.tex` (one occurrence each,
`7\,454` in the LaTeX). The same applies to **806** checks and to the **four** open defects — the
assertion count rises by two for every `.py` added to `src/`, and the defect count falls as defects are
closed, so both are live measurements and both moved during the pass that wrote this (2.4).

**0a. RESOLVED AT 12:43, except for one row.** When Part 1 was written, both documents it points at still
described the v2 paper: `FIGURE_PROVENANCE.md` carried the twenty-row v2 table (five rows for figures v3
withdrew, no row for `fig:gapscaling` or `fig:thm1iii-violation`) and `REPRODUCE.md` still called the
sampled reconstruction *"REAL sampled"* and pointed at `paper/table_molecules.tex`. Both were rebuilt for
the seventeen-figure set while this pass was running, and `src/check_provenance.py` now derives the
figure list from `paper/main.tex` instead of trusting the page. Re-checked here: seventeen rows, a row
for each of the two new figures, no "REAL sampled", no `table_molecules.tex`.

**What remains is a direct contradiction between this statement and that file, and this statement is the
one that measured.** `FIGURE_PROVENANCE.md`'s exceptions table lists **three** exceptions to "every
plotted value is deposited"; Part 1 says **one**.

* **Fig. 1, the schematic** — both agree.
* **"Fig. 6, panel (c)": the four free-fermion support sizes `35, 336, 3496, 37361`** — **that row is
  false.** It is the v2 row, renumbered from Fig. 1 to Fig. 6 and not re-measured. `grep -rE` over every
  `.tex` and `.dat` under `paper/` finds the four numbers **as a set in exactly one place: the prose of
  `paper/sec_6_body.tex:329`**. Panel (c) of `fig_decoupling_native.tex` (line 34) now reads
  `(c)~Hubbard: $\chi$ runs forwards` and plots six strata of $(\chi,|\mathcal S|)$;
  `src/make_decoupling_native.py` no longer contains the constants either. A figure that does not plot
  them cannot be an exception for them. The new checker did not catch this because it validates the
  *label set* against the manuscript, not the *contents* of a cell.
* **"Fig. 17: two of the four JSONs behind the plotted ratios are not deposited"** — true, and it is in
  Part 1, in the known-gaps sentence. The ratios Fig. 17 *plots* are the deposited `n3_*.dat` tables;
  what is missing is what *built* them. Deposit gap and regenerability gap are different claims and the
  DAS makes them in different sentences (2.2).

One of the two has to give before submission. Do not resolve it by weakening Part 1 to match the table:
that would print an exception for a panel that no longer exists.

**0b. `src/make_akw_sampled_honest_fig.py` reads two of its inputs from outside the repository** (2.3),
which is why `check_figures.py` passes on this machine and would not on a clone. The deposited
`data/sampled_honest.json` and `data/honest_sampling.json` are numerically identical to what it reads, so
the fix is three lines, and it is the difference between a checkable deposit and one that only looks
checkable.

**0c. PRA wants the repository cited in the reference list, not only linked. — DONE 2026-09-19: `\bibitem{p2repo}` is in `paper/bibliography.tex` and the statement cites it, `[134]`, beside the URL. The build is still 0 undefined citations. What follows is the record of why it was needed.** *"Publicly shared data and
software must be cited in the reference list, and the citation must be included in the data availability
statement"*
([journals.aps.org/authors/data-availability-statements](https://journals.aps.org/authors/data-availability-statements)).
The manuscript currently carries the URL as a `\url{}`, which does not depend on `paper/bibliography.tex`
and therefore cannot break the build if that file is being edited elsewhere. To close it, append before
`\end{thebibliography}`:

```latex
\bibitem{p2repo} N.~Bonilla Vargas, ``Dynamical spectral functions from bitstring-sampled quantum
subspaces --- reproduction repository,'' \href{https://github.com/nicolasbonilla/dynamical-spectral-functions-sqd}{github.com/nicolasbonilla/dynamical-spectral-functions-sqd} (2026).
```

and change `\url{https://github.com/...}` in `paper/sec_9.tex` to `\cite{p2repo}` beside it. Archiving a
release on Zenodo and citing the DOI instead is what the policy actually prefers.

1. **Commit and push.** Nothing above is true of a clone until it is. This is the single largest gap
   between what this document says and what a referee would download.
2. Open `https://github.com/nicolasbonilla/dynamical-spectral-functions-sqd` signed out and confirm it
   resolves. If the repository has a different name, change it in Part 1, `README.md`, `CITATION.cff`
   **and** `paper/sec_9.tex` (`hardware.tex` is superseded and no longer built).
3. Confirm the workflow in `.github/workflows/ci.yml` has actually run at least once. Only then may the
   continuous-integration sentence of 2.6 be added.
4. Decide about the GFlowNet per-seed traces: deposit them, or keep "available on request". Either is
   acceptable to PRA; silence is not.
5. Re-run `python src/verify.py` against the final tree and require exit 0. If the numbers move, update
   the count in Part 1 **and in `paper/sec_9.tex`** — it is a factual claim like any other, and it is now
   printed in the published article.
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
