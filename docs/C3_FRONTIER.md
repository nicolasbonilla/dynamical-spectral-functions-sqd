# C3 — THE CERTIFICATE FRONTIER SWEEP, deposited

*Deposited 2026-09-19. Until this date Sec. V was the one part of the manuscript a reader could not
regenerate: its three tables came from a 204-row sweep that existed only in a working tree, and
`sec_5_body.tex` said so about itself. The sweep is here now — 54 measured points, the grid they
aggregate to, the three gates that had to pass before any of it was believed, and the code that
produced all of it.*

The companion documents are [`REPRODUCE.md`](REPRODUCE.md) (how to run everything else),
[`FILE_INDEX.md`](FILE_INDEX.md) (what every file is) and
[`KNOWN_DISCREPANCIES.md`](KNOWN_DISCREPANCIES.md) (where the deposit and the code disagree).

---

## 1. What is here

| path | contents |
|---|---|
| `data/c3_frontier/C3_grid.json` | the sweep: **204 grid rows** (L × fraction × η) and a 16-row frontier block (one crossing per (L, η) cell) |
| `data/c3_frontier/pt_L*_FR*_nl*.json` | **54 point files** — one per (L, fraction, Krylov depth). Every number in the grid comes from these |
| `data/c3_frontier/ties.json` | gate: how many determinants in `S` carry exactly zero accumulated Born weight, and how the `argsort` ties at the cut are treated |
| `data/c3_frontier/validation.json` | gate: the cross-validation of the memory-lean algebra against a dense Rayleigh–Ritz Λ_S at L=8 (`fvalid.py B`) |
| `data/c3_frontier/validation_full.json` | the same gate re-run from this repository on 2026-09-19 with **no argument**, so it also carries part (A): the six algebraic controls A1–A6 at L=6, all PASS |
| `data/c3_frontier/ladder_L{4,6,8}.json` | the exact five-rung instrumentation of the Theorem-B proof (no Krylov anywhere): where the slack is lost, rung by rung |
| `data/c3_frontier/verdict_frontier.json` | 54 rows: where the frontier sits at each rung of the proof, i.e. how much of it is Cauchy–Schwarz and how much is physics |
| `data/c3_frontier/published_fraction/` | the four points measured at the fractions the manuscript actually publishes (`_PUB`), their summary `PUB_rows.json`, and the L=14 point `L14_FR008.json` with its run log |
| `data/c3_frontier/adversarial/` | the four adversarial re-runs named in the provenance note of `sec_5_body.tex` |
| `data/c3_frontier/null/` | the post-processing that tests the certificate against the null a referee will propose (`1−w_S`), and the out-of-sample test of the calibration |
| `data/c3_frontier/C3_tables.txt`, `L1*.txt` | the tables and the run logs of the sweep, as produced |
| `src/frontier/` | **24 Python files.** Group A the engine, B the gates, C the proof ladder, D the adversarial re-runs, E the published-fraction rows, F the post-processing |

Total: **2.47 MB** in 82 files. `data/` is excluded from the arXiv package by
`src/build_arxiv_bundle.py` (`FORBIDDEN_PATH_PARTS` contains `data`, `FORBIDDEN_EXTS` contains
`.py`), so this deposit adds **nothing** to the submission.

That was checked member by member rather than by comparing two file sizes, because the file sizes do
not match and the reason is not this deposit. Rebuilding the bundle after depositing gives **66
members — the same 66, none added, none removed, and not one of them from `data/` or `src/`**. The
archive did grow, from 2 082 244 to 2 087 913 bytes, and every byte of that is attributable to
`paper/*.tex` edits made in the same working tree while this deposit was being prepared:
`sec_5_body.tex` +10 959, `sec_5_app.tex` +3 941, `figs/fig5_caption.tex` +474 and eight smaller
ones, +15 661 uncompressed in total. Comparing the sizes alone would have let this deposit take
credit — or blame — for someone else's edit.

---

## 2. The measured cost — the estimate was 25–50× pessimistic

`P2_VERDAD_ESTABLECIDA.md` §11, `P2_PLAN.md` and `P2_VEREDICTO_PUERTAS.md` all budget C3 at
**“~20–40 h of CPU and ~9 GB of RAM”**. That figure was written *before the sweep was run*. The run
is now in the deposit and carries its own stopwatch: every point file records `t_point` and
`peak_vec_GiB`, and the grid carries them forward.

The numbers below are **recomputed from `data/c3_frontier/` by `src/verify.py` on every run** and
compared against this block. If the deposit changes and this document does not, the build fails.

```
C3_COST_SECONDS  = 2873.9
C3_COST_HOURS    = 0.798
C3_PEAK_GIB      = 1.6357
C3_UNIQUE_POINTS = 51
C3_GRID_ROWS     = 204
C3_POINT_FILES   = 54
C3_CONVERGED     = 187
C3_MARGINAL      = 6
C3_NOTCONVERGED  = 11
```

Split by size, in seconds of subspace work: **L=6 41.3 · L=8 288.6 · L=10 969.5 · L=12 1574.5**.
The peak of 1.64 GiB is a single L=12 point at FR ≥ 0.85 with `n_l = 500`; everything below L=12
stays under 0.3 GiB. **0.80 h and 1.6 GiB, not 20–40 h and 9 GB.**

Two honesty notes on that comparison, because it flatters this deposit:

* `t_point` is **subspace work only** — the Lanczos basis and the Gram matrix. It excludes the
  ground state, the ranking and the reference spectrum, which are computed once per `(L, n_ref)`
  and are the larger share at small L. The wall clock of the whole sweep was longer than 0.80 h.
* The η sweep really is free: `V` and the Gram do not depend on η, so the four resolutions of each
  row cost one subspace each, not four. That is why the grid has 204 rows and only 51 subspaces.

The same over-estimation happened in the other direction elsewhere and is recorded rather than
quietly dropped: for the four published-fraction points of `published_fraction/`, a prior estimate
of “41 s of core time” came out at **128.7 s of core and 471 s of wall clock** — 3× and 11×
optimistic.

---

## 3. How to regenerate every deposited artefact

Run from the repository root. Outputs go to `$C3_OUT`, which defaults to `data/c3_frontier` — **set
it to a scratch directory to compare against the deposit instead of overwriting it**, exactly as
`src/check_figures.py` does for the figure fragments:

```
set C3_OUT=%TEMP%\c3check          (Windows)      export C3_OUT=/tmp/c3check   (POSIX)
```

Large temporary Krylov memmaps go to `$C3_TMP`, default `build/` (which is `.gitignore`d).

| artefact | command | cost here |
|---|---|---|
| `C3_grid.json` | `python src/frontier/fsum.py` | < 1 s |
| the 54 `pt_*.json`, L ≤ 10 | `python src/frontier/frun.py <L> <n_l> <FR,FR,…> <n_ref> [n_ref_o] [tag]` | 3 s (L=6) … 2 min (L=10) per point |
| the 16 `pt_L12_*.json` | `python src/frontier/frun12.py <L> <n_l> <FR,…> <n_ref> <memmap> [tag]` | ≈ 110 s per point |
| `validation.json` | `python src/frontier/fvalid.py B` (or no argument for A **and** B) | 2 min 23 s |
| `ties.json` | `python src/frontier/ftie.py` | 38 s |
| protocol gate | `python src/frontier/fproto.py` | 46 s — **needs a checkpoint, see §5** |
| `ladder_L4.json` | `PER_ETA=12 python src/frontier/run_ladder.py 4 0.2,0.4,0.6,0.8,0.9` | 4 s |
| `ladder_L6.json` | `PER_ETA=16 python src/frontier/run_ladder.py 6 0.10,0.20,0.30,0.40,0.50,0.56,0.60,0.70,0.80,0.85,0.90,0.93,0.95,0.97,0.98,0.99` | 14 s |
| `ladder_L8.json` | `PER_ETA=16 python src/frontier/run_ladder.py 8 0.10,0.20,0.30,0.40,0.50,0.56,0.625,0.70,0.80,0.85,0.90,0.95,0.97,0.99` | ≈ 3 min |
| `verdict_frontier.json` | `python src/frontier/verdict.py 4 6 8` | < 1 s |
| `published_fraction/PUB_rows.json` | `python src/frontier/pubsum.py` | < 1 s |
| `published_fraction/L14_FR008.json` | `python src/frontier/run_L14.py` | 37 min — **needs two checkpoints, see §5** |
| `adversarial/z_fine_L6.json` | `python src/frontier/z_fine.py 6 255,285,291,294,295,296,297,298,299 0.18,0.15,0.05 check` | 6 s |
| `adversarial/z_rank_L6.json` | `python src/frontier/z_rank.py 6 0.18` | 5 s |
| `adversarial/z_e0_L6.json` | `python src/frontier/z_e0.py 6 0.85` | 10 s |
| `adversarial/z_suppK_L8.txt` | `python src/frontier/z_suppK.py 6 8` | 1 min 11 s |
| `null/CALIB_OOS.json` | `python src/frontier/calib_oos.py` | < 1 s |
| `null/NULL_WS4.json` | `python src/frontier/null_ws4.py` | < 1 s |
| the null’s out-of-sample tables | `python src/frontier/null_ws5.py`, `null_ws6.py`, `null_ws7.py` | < 1 s each, stdout only |

**`PER_ETA` matters and is recorded in each ladder file** (`per_eta`): the L=4 ladder was run at 12,
the L=6 and L=8 ladders at 16. Running them at the default 12 changes `ngrid` and therefore the
quadrature, and the files will not match.

### Measured, not asserted — every one of the 24 files was run from this repository on 2026-09-19

| what was re-run | against the deposit |
|---|---|
| `fsum.py` | **`C3_grid.json` regenerated exactly** — deep-equal, 204 rows, 16 frontier rows |
| `verdict.py 4 6 8` | **`verdict_frontier.json` regenerated exactly** |
| `calib_oos.py`, `null_ws4.py`, `pubsum.py` | **`CALIB_OOS.json`, `NULL_WS4.json`, `PUB_rows.json` regenerated exactly** |
| `frun.py 6 400 0.20 290 290` | `pt_L6_FR0.20_nl400.json`: headline columns (w_S, Λ̂/η, bounds, rel-L1) to **3.5 × 10⁻¹⁶** |
| `frun.py 8 500 0.56 500` | `pt_L8_FR0.56_nl500.json`: headline columns to **6.8 × 10⁻¹⁴** |
| `frun12.py 8 150 0.85 500 1 _coordcheck` | `pt_L8_FR0.85_nl150_coordcheck.json`: headline columns to **6.0 × 10⁻⁹** |
| `ftie.py` | `ties.json` to **2.1 × 10⁻¹⁰** |
| `fvalid.py` | part (A) **PASS** (A1–A6), part (B) reproduces Λ_S to 2 × 10⁻¹³ and Λ_K to 5 × 10⁻⁷ |
| `fproto.py` (with the checkpoint) | ranking set-overlap **1.000000** with the repository’s own `L10_order.npy` at every fraction, K=18 |
| `run_ladder.py 6` at `PER_ETA=16` | `ladder_L6.json`: `ngrid` identical everywhere; every other field to **1.6 × 10⁻³** or better |
| `run_ladder.py 4` at `PER_ETA=12` | **does not match** — see §4 |
| `z_rank.py`, `z_e0.py` | `z_rank_L6.json`, `z_e0_L6.json` regenerated **exactly** |
| `z_fine.py`, `z_suppK.py` | run; `z_fine_L6.json` reproduces, `z_suppK` reprints the L=6/L=8 support census |
| `run_L14.py` | declared gap (§5): exits 3 with a message naming the missing checkpoint |

The residual deviations are ARPACK’s: every `eigsh` in the C3 engine uses an unseeded start vector,
so a re-run agrees to solver tolerance, not in the last bit. The quantities where the deviation is
largest — `beta_n` at Lanczos breakdown, `Lambda_K_in` at 10⁻⁷, `plateau/rel_change` at 10⁻¹⁴ — are
truncation residuals, i.e. numbers whose whole content is that they are small.

---

## 4. The one stale artefact, declared

`data/c3_frontier/ladder_L4.json` was written at 17:32 on 2026-09-18. Four minutes later
`fron_lib.make_grid` widened its integration window from a constant `pad = 4.0` to
`pad = max(4, 40 η)`. The L=4 ladder therefore predates its own library; the L=6 and L=8 ladders
(17:40 and 17:59) do not.

Re-running it today reproduces every field except the ω grid: `ngrid` at η = 0.50 goes from **1497
to 2265**, and the finer, wider quadrature moves the numbers that depend on it.

**Measured consequence, so that nobody has to guess at it:** feeding the re-run `ladder_L4.json`
into `verdict.py` changes **only the L=4 rows** of `verdict_frontier.json` — 72 of 54 × … entries,
all at L=4 — by at most **2.8 % relative on `d_CS`** (7 × 10⁻⁴ absolute in fraction units). Every
L=6 and L=8 row is bit-identical. No number in the manuscript comes from the L=4 ladder.

This is registered in `src/verify.py` as `KNOWN_OPEN["ladder_L4.ngrid_eta050"]`, which pins 1497 and
names 2265 as the value the current code produces. Regenerating the file turns that `[XFAIL]` into
an `[XPASS]` telling you to delete the registry entry. It is not silently tolerated and it is not
silently repaired.

---

## 5. What is NOT deposited, and why

A declared gap is worth more than a covered one. Three remain.

1. **`work/ckpt/L10_order.npy`** — the ranking `src/scaling_lanczos_mf.py` checkpoints. `ckpt/` is
   `.gitignore`d. `fproto.py`, the gate that proves the sweep uses the repository’s own ranking
   protocol rather than merely claiming to, needs it. Without it the script **exits 3 with a message
   naming the file**, instead of a traceback. Point `P2_CKPT` at a tree that has it, or regenerate
   it with `python src/scaling_lanczos_mf.py 10`. With it: set-overlap 1.000000 at every fraction.

2. **`work/ckpt/L14_gs.npz` (94 MB) and `L14_order.npy` (39 MB)** — the L=14 ground state and
   ranking. `run_L14.py` reuses them rather than spending three hours recomputing a ground state;
   they are too large for the deposit. The **result** of that run is deposited
   (`published_fraction/L14_FR008.json`, with `L14_run.txt`), and the file itself records that the
   **ranking protocol was not re-run at L=14** — it is inherited from the checkpoint, verified at
   L=10 by set overlap. Same exit-3 behaviour as above.

3. **`z_fine_L8`** — the L=8 half of the fine-grained frontier scan. Its log is deposited
   (`adversarial/z_fine_L8.txt`, 1089 s of measurement) but it was not re-run in this pass; only the
   L=6 half was. Nothing in the manuscript depends on a number that appears only there.

Not deposited because it is not part of C3, and named so that the absence is not mistaken for an
oversight: the exploratory scripts of the same working tree (`tab.py`, `tighten.py`, `slackL.py`,
`s5model.py`, `selftest.py`, `selftest2.py`, `hyp_*.py`, `run_ladder`’s hypothesis probes, and the
`atk/` attack directory with its 5.8 MB `.npy` intermediates). They produced no number that reaches
the manuscript. `null_ws.py`, `null_ws2.py` and `null_ws3.py` are likewise absent: they are the
three earlier, weaker versions of the null test, superseded by `null_ws4.py`–`null_ws7.py`, which
are here.

---

## 6. The path repairs — three lines, not two, and the import graph untouched

The C3 working tree lived in two scratch directories and reached the repository by absolute path.
Depositing it required moving **four** path constants (the plan predicted two), plus routing the
outputs through `$C3_OUT`. No module was renamed, no import was rewired, no function moved:

| file | was | is |
|---|---|---|
| `c0_lib.py:18` | `SRC = r"C:\…\release\src"` | `SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))` |
| `fcore.py:27` | `PUERTAS = os.path.join(os.path.dirname(HERE), 'p2_puertas')` | `PUERTAS = HERE` |
| `fron_lib.py:36` | `PUERTAS = r"C:\…\scratchpad\p2_puertas"` | `PUERTAS = os.path.dirname(os.path.abspath(__file__))` |
| `fproto.py:16` | `CK = r"C:\…\work\ckpt"` | `CK = os.environ.get('P2_CKPT') or os.path.join(ROOT, 'ckpt')` + an exit-3 guard |
| `z_e0.py`, `z_fine.py`, `z_rank.py`, `z_suppK.py`, `run_L14.py`, `pubsum.py` | the same two scratch constants | `HERE` / `$C3_OUT` |

Two facts about that repair worth recording, both found by running the code rather than reading it:

* **`fvalid.py` already binds the name `OUT` to its results dict.** The deposit block therefore
  binds the directory to `OUTDIR` in that file alone. Injecting `OUT` blindly produced
  `TypeError: expected str … not dict` after 2 minutes 23 seconds of correct computation.
* **`fsum.py` now skips `*_PUB.json` explicitly.** Its glob is `pt_L*_nl*.json`, and the four
  published-fraction points measured on 2026-09-19 match it. Aggregating them would add rows and
  silently change the 204 the manuscript quotes, so they live in `published_fraction/` and the skip
  is written down in the file with the reason. This is what makes `python src/frontier/fsum.py`
  reproduce the deposited grid exactly rather than approximately.

A deposit-wide scan asserts that no file under `src/frontier/` contains `C:\Users`, `scratchpad`,
`p2_puertas` or `p2_frontera` — in a path *or in prose*; the three docstrings that named the old
scratch layout were rewritten, not left to mislead.

---

## 7. What the guardian checks, and what it cannot

`src/verify.py` section **(10b)** opens the 54 point files — not the grid — and:

* re-derives `bound_leak`, `bound_trivial`, `cert`, `in/out`, `slack`, the non-triviality flag and
  the convergence flag for all 204 rows from Theorem B’s algebra and the classifier. Measured
  agreement: **exact** (max relative deviation 0.0);
* checks **Theorem B itself** on every row that carries a measured error —
  `(1 − w_S) ≤ rel-L₁ ≤ min{leak, trivial}`. **0 violations in 204 rows.** This is the only place in
  the repository where the manuscript’s central inequality is checked against a *measured*
  reconstruction error rather than a synthetic sweep;
* re-derives `Λ̂ = Λ/‖φ‖` and `Λ² = Λ_in² + Λ_out²` from each point file’s own `norm_phi2`. The
  missing `‖φ‖` would be a silent factor √2 on every bound in Sec. V;
* recomputes the cost and the convergence census and compares them with the stanza in §2 of this
  document;
* asserts the Sec. V headline: at L = 12 and the **published** fraction 0.18, the certificate is
  vacuous at all four resolutions (leak/trivial = 7.51, 4.86, 4.00, 2.81);
* asserts closure both ways — no grid row without a deposited point file, and no deposited point
  file that is neither aggregated nor superseded by a deeper `n_l` at the same (L, FR). Three files
  are superseded and are named: `pt_L8_FR0.85_nl150_coordcheck.json`, `pt_L12_FR0.50_nl300.json`,
  `pt_L12_FR0.60_nl300.json`.

2 387 numeric assertions, with a coverage floor of the same number in `SECTION_FLOORS`: if a point
file is deleted the assertions vanish and the floor fires.

**A PASS is worth nothing without a control that fires, so seven were run** — each on a full copy of
the repository, never on the tree itself. All seven fired and all seven named the defect:

| control | result |
|---|---|
| delete one deposited point file | exit 1, *closure* |
| multiply one `bound_leak` in the grid by 1.05 | exit 1, *Theorem B's algebra* |
| set one row's `rel-L₁` to 1.5 × its certificate | exit 1, *upper bound* |
| drop the `‖φ‖` from one point's `Λ̂` (a factor √2) | exit 1, *normalisation* |
| let this document drift back to `C3_COST_HOURS = 20.0` | exit 1, *`C3_COST_HOURS`* |
| silently "repair" the stale `ladder_L4` grid to 2265 | exit 0 with `[XPASS]` asking for the registry entry to be deleted |
| delete `data/c3_frontier/` entirely | exit 1, *`c3_frontier`* |

**What it does not do:** it does not re-run the sweep. A grid whose *primitives* were all corrupted
consistently would still satisfy every relation above. What defends against that is §3 — the fact,
measured and printed there, that the code in `src/frontier/` reproduces the deposited primitives
from scratch.
