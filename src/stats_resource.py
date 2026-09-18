# -*- coding: utf-8 -*-
r"""stats_resource.py -- Correct statistical analysis of the resource thesis (Sec. "One-body
magic and its decoupling from the sampling cost |S|").

WHAT IT COMPUTES
----------------
The manuscript summarises the resource evidence with six pooled Spearman coefficients and one
significance statement.  Every one of those coefficients is reproduced here as a control, and then
re-analysed under the design that actually generated the data.  Six blocks:

  (1) HUBBARD, STRATIFIED.  data/cost_vs_ent.json is not 30 exchangeable points: it is a 6 x 5
      design, six (L,N) strata each swept over U in {1,2,4,8,16}.  Within every stratum
      rho(F1,|S|) = -1 exactly.  The pooled rho = -0.11 is Simpson's paradox, not a null result.
      Reported: per-stratum rho, exact permutation null of the summed statistic (obtained by
      complete enumeration of the 5! relabellings per stratum with branch-and-bound over the
      6-fold product), Fisher combination of the exact per-stratum p-values, and stratified
      Kendall concordance as an effect size.

  (2) MOLECULES, CLUSTERED + TIED.  The n=38 molecular points are 19 molecules x 2 geometries
      (dependent pairs) and 11 of them sit at the floor |S| = 1 (a 29% tie block).  Neither the
      asymptotic Spearman p-value nor n=38 is admissible.  Reported: the tie structure, a
      cluster-level (molecule-level) permutation test, the intracluster correlation, the
      resulting effective n, and a cluster-bootstrap confidence interval.

  (3) DEPENDENT-CORRELATION COMPARISONS the manuscript never performs.  chi (rho=0.90) is claimed
      to beat F1 (rho=0.86) on the SAME 38 points, and chi (0.72) to beat F1 (0.60) on the SAME
      74.  Two correlations sharing a variable must be compared with Steiger's / Williams' test,
      not read off side by side.  Reported: Williams' t, Steiger's z, and a cluster-bootstrap
      interval on the difference; and an explicit verdict on whether the difference is resolvable.

  (4) chain2 AND ladder (n=3 each).  In both, rho(F1,|S|) = +1.000 exactly while rho(chi,|S|) is
      0.50 and 0.87: in 2 of the 4 families F1 outranks chi.  The manuscript does not mention it.
      Reported: the exact permutation floor at n=3, and the system-size confound (both families
      are L-sweeps at fixed U, and F1/L is constant to 3 digits, so F1 is a proxy for L).

  (5) GFLOWNET paired comparison.  The manuscript contrasts "3.1 s.d." against "1.2 s.d." using
      DIFFERENT denominators.  With 5 paired seeds the correct statistic is a paired t on the
      differences.  data/gflow.json deposits only 5 summary rows, so the paired t is not
      computable; this script instead derives the rigorous INTERVAL of paired t consistent with
      the deposited means and s.d.'s, and states which of the two contrasts is decidable.

  (6) FORENSICS of the orphan constant pooled_spearman_F1_S = 0.5991983669238197 in
      data/resource_master.json, which does not reproduce from the 74 deposited points.

WHICH CLAIMS OF THE PAPER THIS SUPPORTS
---------------------------------------
  * resource.tex: "chi holds at rho=0.57 while F1 collapses to rho=-0.11 (n=30, p~0.56, not
    significant)" -- REFUTED in the direction that STRENGTHENS the paper's thesis: stratified,
    rho(F1,|S|) = -1 in all six strata, p ~ 3e-13.  F1 does not merely fail to predict the cost;
    it predicts it with a perfect inverted sign.
  * resource.tex: "both p < 0.01" for the two tracking coefficients -- the molecular half of that
    statement is not supported by the asymptotic test actually used (clustered + tied data); it
    IS supported by the cluster permutation test computed here.
  * resource.tex / fig:master: the side-by-side 0.90-vs-0.86 and 0.72-vs-0.60 comparisons -- the
    formal test is supplied here, with a verdict.

HOW TO RUN
----------
    python src/stats_resource.py [--root .] [--out data/stats_resource.json]
                                 [--boot 20000] [--perm 1000000]

Inputs (read-only, relative to --root):
    data/cost_vs_ent.json, data/resource_master.json, data/n19_suite.json, data/gflow.json

RUNTIME
-------
~20 s on one core of a laptop (numpy 1.26 / scipy 1.13) at the default resampling sizes.
Memory < 200 MB.  No QPU, no network, no compiled extension beyond numpy/scipy.

SELF-VERIFICATION (mandatory; the script aborts on any mismatch)
---------------------------------------------------------------
Before any new statistic is reported, the script recomputes the six coefficients the manuscript
prints and compares them to hard-coded literals at 1e-4 absolute tolerance:
    molecular n=38 : rho(chi,|S|) = 0.9011 , rho(F1,|S|) = 0.8618
    hubbard  n=30  : rho(chi,|S|) = 0.5692 , rho(F1,|S|) = -0.1088
    pooled   n=74  : rho(chi,|S|) = 0.7224 , rho(F1,|S|) = 0.5992
It additionally checks the 6x5 stratum design, the 19x2 molecular pairing (matched by exact F1
value against data/n19_suite.json), and the exact-permutation machinery against the closed-form
Spearman null for n=5.  Any failure raises SystemExit(1) and nothing is written.

Deterministic: every random number comes from numpy Generator(PCG64(seed)) with the seeds fixed
in SEEDS below.  Bit-for-bit reproducible on a given numpy version.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import platform
import sys
import time

import numpy as np
import scipy
from scipy import stats

# --- NUMPY_TRAPEZOID_BRIDGE ------------------------------------------------
# numpy 2.0 ADDED np.trapezoid and REMOVED np.trapz.  Files in this repository use
# both names, so without this bridge no single numpy version runs the whole deposit.
if not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz          # numpy < 2.0
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid          # numpy >= 2.0
# ---------------------------------------------------------------------------

# Repository layout: this file lives in <repo>/src/, the deposit in <repo>/data/.
# NO container paths: every default is resolved from __file__.
SRC_DIR = os.path.dirname(os.path.abspath(__file__))
REPO    = os.path.dirname(SRC_DIR)


def _repo_rel(p):
    """A path as it must appear in a DEPOSITED provenance block: relative to the
    repository, forward slashes, never the author's machine.  An absolute path in a
    deposited file is both a privacy leak and a dangling pointer -- it names a
    directory the reader does not have.  Added 2026-09-18 (second pass) after an
    audit found this file's own output carrying an absolute machine path."""
    if not p:
        return p
    a = os.path.abspath(str(p))
    r = os.path.abspath(REPO)
    if os.path.normcase(a).startswith(os.path.normcase(r) + os.sep):
        return os.path.relpath(a, r).replace(os.sep, "/")
    return "<not part of the deposit: %s>" % os.path.basename(a)


DEFAULT_OUT = os.path.join(REPO, "data", "stats_resource.json")

# ----------------------------------------------------------------------------------------------
# Fixed RNG seeds.  Do not change without also updating the provenance block consumers.
# ----------------------------------------------------------------------------------------------
SEEDS = {
    "cluster_permutation_molecular": 20260918,
    "cluster_bootstrap_molecular": 20260919,
    "family_bootstrap_pooled": 20260920,
}

# Control values printed by the manuscript.  Self-verification targets.
CONTROLS = {
    "molecular_rho_chi_S": 0.9011,
    "molecular_rho_F1_S": 0.8618,
    "hubbard_rho_chi_S": 0.5692,
    "hubbard_rho_F1_S": -0.1088,
    "pooled_rho_chi_S": 0.7224,
    "pooled_rho_F1_S": 0.5992,
}
CONTROL_TOL = 1e-4

U_GRID = (1.0, 2.0, 4.0, 8.0, 16.0)


def srho(a, b):
    """Spearman rho with mid-rank tie correction (scipy convention)."""
    return float(stats.spearmanr(np.asarray(a, float), np.asarray(b, float)).statistic)


def fail(msg):
    print("SELF-VERIFICATION FAILED: " + msg, file=sys.stderr)
    raise SystemExit(1)


# ----------------------------------------------------------------------------------------------
# Vectorised rank machinery.  Spearman's rho is Pearson's r on mid-ranks, so a whole batch of
# resamples can be ranked at once: encode every value by its ordinal in the global sorted-unique
# pool, count the ordinals present in each resample, and read the mid-rank off the cumulative
# counts.  This is EXACT -- it reproduces scipy.stats.spearmanr to machine precision, which the
# self-verification below checks explicitly -- and it is what makes B = 1e6 affordable on a
# laptop.  Naive looping over scipy would take ~90 minutes for the same numbers.
# ----------------------------------------------------------------------------------------------
def codes_of(values):
    """Ordinal encoding of `values` against their own sorted-unique pool.  Returns (codes, K)."""
    u = np.unique(np.asarray(values, float))
    return np.searchsorted(u, np.asarray(values, float)).astype(np.int64), int(u.size)


def batch_midranks(samp_codes, K):
    """Mid-ranks of each row of `samp_codes` (shape (B, n), ordinal codes in [0, K))."""
    B, n = samp_codes.shape
    flat = (samp_codes + (np.arange(B, dtype=np.int64)[:, None] * K)).ravel()
    counts = np.bincount(flat, minlength=B * K).reshape(B, K).astype(np.float64)
    cum = np.cumsum(counts, axis=1) - counts                 # values strictly smaller
    mid = cum + 0.5 * (counts + 1.0)
    return np.take_along_axis(mid, samp_codes, axis=1)


def row_rho(ra, rb):
    """Pearson r between corresponding rows of two rank matrices."""
    a = ra - ra.mean(axis=1, keepdims=True)
    b = rb - rb.mean(axis=1, keepdims=True)
    return (a * b).sum(axis=1) / np.sqrt((a * a).sum(axis=1) * (b * b).sum(axis=1))


# ==============================================================================================
# Data loading
# ==============================================================================================
def load(root):
    p = lambda *a: os.path.join(root, *a)
    with open(p("data", "cost_vs_ent.json")) as f:
        cve = json.load(f)
    with open(p("data", "resource_master.json")) as f:
        rmaster = json.load(f)
    with open(p("data", "n19_suite.json")) as f:
        n19 = json.load(f)
    with open(p("data", "gflow.json")) as f:
        gflow = json.load(f)
    return cve, rmaster, n19, gflow


# ==============================================================================================
# Exact permutation machinery for small strata
# ==============================================================================================
def stratum_null(x, y):
    """Exact permutation null of Spearman rho for one stratum.

    x is held fixed and every one of the n! relabellings of y is applied.  Tie corrections are
    handled by scipy, so this is valid for the tied chi column as well.  Returns the sorted array
    of n! rho values (with multiplicity).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    out = [srho(x, np.asarray(perm, float)) for perm in itertools.permutations(y)]
    return np.sort(np.asarray(out))


def exact_one_sided_p(null, obs, side):
    """P(rho <= obs) ('left') or P(rho >= obs) ('right') under the exact null, with a 1e-12 slack
    so the floating-point image of the observed permutation is counted."""
    n = null.size
    if side == "left":
        return float(np.count_nonzero(null <= obs + 1e-12)) / n
    return float(np.count_nonzero(null >= obs - 1e-12)) / n


def exact_sum_tail(nulls, t_obs, side):
    """EXACT tail probability of T = sum_s rho_s under the product of the per-stratum exact nulls.

    The full product has (n!)^S points (here 120^6 = 2.99e12), far too many to enumerate.  The
    tail is obtained exactly by depth-first search over the strata with a branch-and-bound cut:
    each stratum's values are visited in the favourable order and the branch is abandoned as soon
    as the partial sum plus the best achievable remainder can no longer reach t_obs.  No
    approximation is made -- every product point satisfying the inequality is counted.

    Returns (probability, number_of_product_points_in_the_tail).
    """
    S = len(nulls)
    if side == "right":
        vals = [np.sort(nl)[::-1] for nl in nulls]          # descending
    else:
        vals = [np.sort(nl) for nl in nulls]                # ascending
    best = [float(v[0]) for v in vals]
    suffix = [0.0] * (S + 1)
    for i in range(S - 1, -1, -1):
        suffix[i] = suffix[i + 1] + best[i]
    total = 1
    for v in vals:
        total *= v.size
    eps = 1e-12
    count = 0

    def rec(i, partial):
        nonlocal count
        v = vals[i]
        rest = suffix[i + 1]
        last = (i == S - 1)
        if side == "right":
            for j in range(v.size):
                if partial + v[j] + rest < t_obs - eps:
                    break                                    # descending: all later are worse
                if last:
                    count += 1
                else:
                    rec(i + 1, partial + v[j])
        else:
            for j in range(v.size):
                if partial + v[j] + rest > t_obs + eps:
                    break                                    # ascending: all later are worse
                if last:
                    count += 1
                else:
                    rec(i + 1, partial + v[j])

    rec(0, 0.0)
    return count / total, count


def fisher_combine(pvals):
    """Fisher's method on independent one-sided p-values."""
    pv = np.asarray(pvals, float)
    chi2 = float(-2.0 * np.sum(np.log(pv)))
    df = 2 * pv.size
    return chi2, df, float(stats.chi2.sf(chi2, df))


def stratified_kendall(groups):
    """Stratified Kendall concordance: concordant/discordant pairs accumulated WITHIN strata only
    (no cross-stratum pair is ever formed).  Returns (tau, nC, nD, nTies)."""
    nc = nd = nt = 0
    for gx, gy in groups:
        gx = np.asarray(gx, float)
        gy = np.asarray(gy, float)
        m = gx.size
        for i in range(m):
            for j in range(i + 1, m):
                s = np.sign(gx[i] - gx[j]) * np.sign(gy[i] - gy[j])
                if s > 0:
                    nc += 1
                elif s < 0:
                    nd += 1
                else:
                    nt += 1
    denom = nc + nd
    tau = (nc - nd) / denom if denom else float("nan")
    return float(tau), nc, nd, nt


# ==============================================================================================
# Dependent-correlation comparisons
# ==============================================================================================
def williams_t(r_xz, r_yz, r_xy, n):
    """Hotelling-Williams t for H0: rho_xz == rho_yz with x, y, z jointly observed (common z).
    Applied to rank-transformed data this is the standard test for two dependent Spearman
    coefficients that share a variable.  df = n - 3."""
    detR = (1.0 - r_xy ** 2 - r_xz ** 2 - r_yz ** 2 + 2.0 * r_xy * r_xz * r_yz)
    rbar = 0.5 * (r_xz + r_yz)
    num = (r_xz - r_yz) * math.sqrt((n - 1) * (1.0 + r_xy))
    den = math.sqrt(2.0 * ((n - 1) / (n - 3.0)) * detR + (rbar ** 2) * (1.0 - r_xy) ** 3)
    t = num / den
    df = n - 3
    return float(t), int(df), float(2.0 * stats.t.sf(abs(t), df))


def steiger_z(r_xz, r_yz, r_xy, n):
    """Steiger's (1980) z-test for two dependent correlations sharing one variable, in the
    Fisher-z metric (Dunn & Clark covariance)."""
    zf = lambda r: 0.5 * math.log((1 + r) / (1 - r))
    rbar2 = 0.5 * (r_xz ** 2 + r_yz ** 2)
    cov_num = (r_xy * (1.0 - 2.0 * rbar2)
               - 0.5 * rbar2 * (1.0 - 2.0 * rbar2 - r_xy ** 2))
    cov = cov_num / ((1.0 - rbar2) ** 2)
    z = (zf(r_xz) - zf(r_yz)) * math.sqrt((n - 3.0) / (2.0 - 2.0 * cov))
    return float(z), float(2.0 * stats.norm.sf(abs(z)))


# ==============================================================================================
# MAIN
# ==============================================================================================
def main():
    ap = argparse.ArgumentParser(description="Correct statistics for the resource section.")
    ap.add_argument("--root", default=REPO,
                    help="repository root to READ data/ from (default: the repository this file "
                         "lives in).  --root never decides where output is WRITTEN: see --out.")
    ap.add_argument("--out", default=DEFAULT_OUT,
                    help="output JSON (default: <repo>/data/stats_resource.json).  A "
                         "relative path resolves against the CURRENT DIRECTORY -- never "
                         "against --root, which is the repository itself.")
    ap.add_argument("--boot", type=int, default=20000, help="cluster-bootstrap replicates")
    ap.add_argument("--perm", type=int, default=1000000,
                    help="cluster-permutation replicates")
    args = ap.parse_args()

    t_start = time.time()
    cve, rmaster, n19, gflow = load(args.root)
    R = {}
    log = lambda *a: print(*a, flush=True)

    # ------------------------------------------------------------------------------------------
    # SELF-VERIFICATION  (part 1: the six published coefficients)
    # ------------------------------------------------------------------------------------------
    log("=" * 94)
    log("SELF-VERIFICATION")
    log("=" * 94)

    pts = rmaster["points"]
    fam = np.array([p["family"] for p in pts])
    F1_all = np.array([p["F1"] for p in pts], float)
    S_all = np.array([p["S"] for p in pts], float)
    CH_all = np.array([p["chi"] for p in pts], float)
    if len(pts) != 74:
        fail(f"resource_master.json has {len(pts)} points, expected 74")

    mmask = fam == "molecule"
    if mmask.sum() != 38:
        fail(f"{int(mmask.sum())} molecular points, expected 38")
    hmask = fam == "hubbard"
    if hmask.sum() != 30:
        fail(f"{int(hmask.sum())} hubbard points, expected 30")

    got = {
        "molecular_rho_chi_S": srho(CH_all[mmask], S_all[mmask]),
        "molecular_rho_F1_S": srho(F1_all[mmask], S_all[mmask]),
        "hubbard_rho_chi_S": srho(CH_all[hmask], S_all[hmask]),
        "hubbard_rho_F1_S": srho(F1_all[hmask], S_all[hmask]),
        "pooled_rho_chi_S": srho(CH_all, S_all),
        "pooled_rho_F1_S": srho(F1_all, S_all),
    }
    for k, target in CONTROLS.items():
        d = abs(got[k] - target)
        mark = "OK " if d <= CONTROL_TOL else "BAD"
        log(f"  [{mark}] {k:24s} recomputed {got[k]:+.6f}   paper {target:+.4f}   |d|={d:.2e}")
        if d > CONTROL_TOL:
            fail(f"{k}: recomputed {got[k]!r} vs paper {target!r}, |d|={d:.3e} > {CONTROL_TOL}")
    R["self_verification_controls"] = {
        k: {"recomputed": got[k], "paper": CONTROLS[k], "abs_diff": abs(got[k] - CONTROLS[k])}
        for k in CONTROLS}

    rows = cve["rows"]
    if len(rows) != 30:
        fail(f"cost_vs_ent.json has {len(rows)} rows, expected 30")
    if abs(srho([r["FAF"] for r in rows], [r["S"] for r in rows]) - cve["r_FAF_S"]) > 1e-12:
        fail("cost_vs_ent.json r_FAF_S does not reproduce from its own rows")
    if abs(srho([r["chi"] for r in rows], [r["S"] for r in rows]) - cve["r_chi_S"]) > 1e-12:
        fail("cost_vs_ent.json r_chi_S does not reproduce from its own rows")
    log("  [OK ] cost_vs_ent.json self-consistent (r_FAF_S, r_chi_S reproduce from its rows)")

    # SELF-VERIFICATION part 2: the exact-permutation machinery against the closed form.
    base = np.arange(5.0)
    nl = stratum_null(base, base)
    closed = np.sort(np.array([1.0 - 6.0 * np.sum((base - np.asarray(p, float)) ** 2) / 120.0
                               for p in itertools.permutations(base)]))
    if nl.size != 120 or np.max(np.abs(nl - closed)) > 1e-12:
        fail("exact permutation null for n=5 disagrees with the closed-form Spearman formula")
    log(f"  [OK ] exact n=5 permutation null: 120 points, "
        f"max|rho_perm - closed form| = {np.max(np.abs(nl - closed)):.2e}")

    # ------------------------------------------------------------------------------------------
    # BLOCK 1 -- HUBBARD, STRATIFIED
    # ------------------------------------------------------------------------------------------
    log("")
    log("=" * 94)
    log("BLOCK 1 -- HUBBARD n=30 : the design is 6 strata (L,N) x 5 values of U")
    log("=" * 94)

    strata = {}
    for r in rows:
        strata.setdefault((r["L"], r["N"]), []).append(r)
    if len(strata) != 6 or any(len(v) != 5 for v in strata.values()):
        fail(f"expected a 6x5 design, found {[(k, len(v)) for k, v in strata.items()]}")
    for k, v in strata.items():
        if tuple(sorted(r["U"] for r in v)) != tuple(sorted(U_GRID)):
            fail(f"stratum {k} does not carry the U grid {U_GRID}")
    log("  [OK ] design verified: 6 strata x 5 U values, identical U grid in each")

    strat_rows = []
    nulls_F1, nulls_chi = [], []
    groups_F1, groups_chi = [], []
    for key in sorted(strata, key=lambda k: (k[0], -k[1])):
        v = sorted(strata[key], key=lambda r: r["U"])
        F = np.array([r["FAF"] for r in v], float)
        S = np.array([r["S"] for r in v], float)
        C = np.array([r["chi"] for r in v], float)
        rF, rC = srho(F, S), srho(C, S)
        nF = stratum_null(S, F)
        nC = stratum_null(S, C)
        pF = exact_one_sided_p(nF, rF, "left")
        pC = exact_one_sided_p(nC, rC, "right")
        nulls_F1.append(nF)
        nulls_chi.append(nC)
        groups_F1.append((F, S))
        groups_chi.append((C, S))
        strat_rows.append({"L": key[0], "N": key[1], "n": 5,
                           "rho_F1_S": rF, "exact_p_one_sided_F1": pF,
                           "rho_chi_S": rC, "exact_p_one_sided_chi": pC,
                           "chi_values": [int(x) for x in C],
                           "n_distinct_chi": int(len(set(C.tolist())))})
        log(f"  (L={key[0]:2d},N={key[1]:2d})  rho(F1,|S|) = {rF:+.4f} [exact p_one={pF:.5f}]"
            f"   rho(chi,|S|) = {rC:+.4f} [exact p_one={pC:.5f}]  chi={[int(x) for x in C]}")
    R["hubbard_strata"] = strat_rows

    T_F1 = sum(s["rho_F1_S"] for s in strat_rows)
    T_chi = sum(s["rho_chi_S"] for s in strat_rows)
    log("")
    log(f"  combined statistic T = sum_s rho_s :  F1 {T_F1:+.6f}   chi {T_chi:+.6f}"
        f"   (attainable range [-6,+6])")

    t0 = time.time()
    pF_exact, nleafF = exact_sum_tail(nulls_F1, T_F1, "left")
    pC_exact, nleafC = exact_sum_tail(nulls_chi, T_chi, "right")
    log("  EXACT permutation p (product of the six 120-point nulls, 120^6 = 2.986e12 points):")
    log(f"      F1  : p_one-sided = {pF_exact:.6e}   ({nleafF} product points in the tail)")
    log(f"      chi : p_one-sided = {pC_exact:.6e}   ({nleafC} product points in the tail)")
    log(f"      [branch-and-bound enumeration, {time.time()-t0:.2f} s; exact, not sampled]")

    chi2F, dfF, pfishF = fisher_combine([s["exact_p_one_sided_F1"] for s in strat_rows])
    chi2C, dfC, pfishC = fisher_combine([s["exact_p_one_sided_chi"] for s in strat_rows])
    log("  Fisher combination of the six exact per-stratum p-values:")
    log(f"      F1  : X2 = {chi2F:.3f}  df = {dfF}  p = {pfishF:.6e}")
    log(f"      chi : X2 = {chi2C:.3f}  df = {dfC}  p = {pfishC:.6e}")

    tauF, ncF, ndF, ntF = stratified_kendall(groups_F1)
    tauC, ncC, ndC, ntC = stratified_kendall(groups_chi)
    log("  EFFECT SIZE (stratified Kendall over the 6x10 = 60 within-stratum pairs):")
    log(f"      F1  : tau = {tauF:+.4f}   concordant {ncF} / discordant {ndF} / tied {ntF}")
    log(f"      chi : tau = {tauC:+.4f}   concordant {ncC} / discordant {ndC} / tied {ntC}")
    log(f"  mean within-stratum rho:  F1 {np.mean([s['rho_F1_S'] for s in strat_rows]):+.4f}"
        f"   chi {np.mean([s['rho_chi_S'] for s in strat_rows]):+.4f}")

    pooled_p_F1 = float(stats.spearmanr(F1_all[hmask], S_all[hmask]).pvalue)
    pooled_p_chi = float(stats.spearmanr(CH_all[hmask], S_all[hmask]).pvalue)
    log("")
    log("  FOR CONTRAST, the pooled (design-ignoring) numbers the manuscript prints:")
    log(f"      F1  : rho = {got['hubbard_rho_F1_S']:+.4f}  asymptotic p = {pooled_p_F1:.4f}"
        f"   <-- reported as 'not significant'")
    log(f"      chi : rho = {got['hubbard_rho_chi_S']:+.4f}  asymptotic p = {pooled_p_chi:.6f}")
    log(f"  => Simpson's paradox: pooled rho(F1,|S|) = {got['hubbard_rho_F1_S']:+.4f} while EVERY")
    log(f"     stratum gives exactly -1.  The correct p is {pF_exact:.3e}, not {pooled_p_F1:.2f}.")

    R["hubbard_stratified"] = {
        "T_F1": T_F1, "T_chi": T_chi,
        "exact_perm_p_one_sided_F1": pF_exact, "exact_perm_tail_points_F1": nleafF,
        "exact_perm_p_one_sided_chi": pC_exact, "exact_perm_tail_points_chi": nleafC,
        "product_space_size": 120 ** 6,
        "fisher_F1": {"chi2": chi2F, "df": dfF, "p": pfishF},
        "fisher_chi": {"chi2": chi2C, "df": dfC, "p": pfishC},
        "stratified_kendall_F1": {"tau": tauF, "nC": ncF, "nD": ndF, "nTies": ntF},
        "stratified_kendall_chi": {"tau": tauC, "nC": ncC, "nD": ndC, "nTies": ntC},
        "mean_within_stratum_rho_F1": float(np.mean([s["rho_F1_S"] for s in strat_rows])),
        "mean_within_stratum_rho_chi": float(np.mean([s["rho_chi_S"] for s in strat_rows])),
        "paper_pooled_p_F1": pooled_p_F1, "paper_pooled_p_chi": pooled_p_chi,
        "effective_n_of_the_pooled_test": 6,
    }

    # ------------------------------------------------------------------------------------------
    # BLOCK 2 -- MOLECULES: 19 clusters of 2, and a 29% tie block at |S| = 1
    # ------------------------------------------------------------------------------------------
    log("")
    log("=" * 94)
    log("BLOCK 2 -- MOLECULES n=38 : 19 molecules x 2 geometries, dependent pairs, |S| floor ties")
    log("=" * 94)

    mols = n19["mols"]
    if len(mols) != 19:
        fail(f"n19_suite.json carries {len(mols)} molecules, expected 19")
    lookup = {}
    for m in mols:
        lookup.setdefault(m["FAF_eq"], []).append((m["mol"], "eq"))
        lookup.setdefault(m["FAF_diss"], []).append((m["mol"], "diss"))
    idx_m = np.flatnonzero(mmask)
    assign = []
    for i in idx_m:
        key = float(F1_all[i])
        if key not in lookup or len(lookup[key]) != 1:
            fail(f"molecular row {i} (F1={key!r}) does not match exactly one n19_suite entry")
        assign.append(lookup[key][0])
    if len(set(assign)) != 38:
        fail("molecular row -> (molecule, geometry) map is not injective")
    log("  [OK ] 38 molecular rows matched 1:1 to 19 molecules x 2 geometries by exact F1 value")

    mol_names = sorted({a[0] for a in assign})
    cid = np.array([mol_names.index(a[0]) for a in assign])
    Fm, Sm, Cm = F1_all[idx_m], S_all[idx_m], CH_all[idx_m]

    uniq, cnt = np.unique(Sm, return_counts=True)
    n_floor = int(cnt[uniq == 1][0]) if (uniq == 1).any() else 0
    n_tied = int(cnt[cnt > 1].sum())
    log(f"  |S| tie structure: {len(uniq)} distinct values over 38 points; "
        f"{n_tied} points ({100*n_tied/38:.0f}%) lie in a tied group")
    log(f"  |S| = 1 floor block: {n_floor} points ({100*n_floor/38:.0f}%) "
        f"share a single mid-rank")
    log(f"  tied groups (|S| value: multiplicity): "
        f"{ {int(u): int(c) for u, c in zip(uniq, cnt) if c > 1} }")

    # Ceiling imposed by the ties: the largest Spearman rho attainable between this tied |S|
    # column and ANY untied predictor is corr(sorted midranks of |S|, 1..38) < 1.  Coefficients
    # must be read against that ceiling, not against 1.
    mid_S = np.sort(stats.rankdata(Sm))
    rho_ceiling = float(np.corrcoef(mid_S, np.arange(1.0, 39.0))[0, 1])
    log(f"  TIE CEILING: with this |S| tie pattern no predictor can exceed rho = {rho_ceiling:.4f};")
    log(f"  every coefficient below should be read as a fraction of that ceiling, not of 1.")

    rho_chi_m = srho(Cm, Sm)
    rho_F1_m = srho(Fm, Sm)
    p_asym_chi = float(stats.spearmanr(Cm, Sm).pvalue)
    p_asym_F1 = float(stats.spearmanr(Fm, Sm).pvalue)
    log("  as published (asymptotic, n=38 assumed i.i.d. and untied):")
    log(f"      rho(chi,|S|) = {rho_chi_m:.4f}  p = {p_asym_chi:.3e}")
    log(f"      rho(F1,|S|)  = {rho_F1_m:.4f}  p = {p_asym_F1:.3e}")
    log(f"  as a fraction of the tie ceiling {rho_ceiling:.4f}:  chi {rho_chi_m/rho_ceiling:.4f}"
        f"   F1 {rho_F1_m/rho_ceiling:.4f}")
    log("      ^ neither p is admissible: the 38 points are 19 dependent pairs and 29% of |S|")
    log("        sits in one tied block.  The cluster permutation test below is.")

    # Intracluster correlation of log2|S| within molecule (one-way random effects, m = 2).
    y = np.log2(Sm)
    grp = np.array([y[cid == g] for g in range(19)])
    grand = y.mean()
    msb = 2.0 * np.sum((grp.mean(axis=1) - grand) ** 2) / (19 - 1)
    msw = np.sum((grp - grp.mean(axis=1, keepdims=True)) ** 2) / (19 * (2 - 1))
    icc = (msb - msw) / (msb + (2 - 1) * msw)
    deff = 1.0 + (2 - 1) * icc
    n_eff = 38.0 / deff
    log(f"  intracluster correlation of log2|S| within molecule: ICC = {icc:.4f}")
    log(f"  design effect = 1 + (m-1)*ICC = {deff:.4f}  ->  EFFECTIVE n = 38/{deff:.4f}"
        f" = {n_eff:.1f}   (not 38)")

    rng = np.random.default_rng(SEEDS["cluster_permutation_molecular"])
    ordidx = np.lexsort((np.array([a[1] for a in assign]), cid))
    Fp = Fm[ordidx].reshape(19, 2)
    Cp = Cm[ordidx].reshape(19, 2)
    Sp = Sm[ordidx].reshape(19, 2)
    if abs(srho(Cp.ravel(), Sp.ravel()) - rho_chi_m) > 1e-12:
        fail("cluster reshaping altered the molecular rho")

    # Molecule-level permutation, done in closed form.  Under a relabelling of the 19 molecules
    # the multiset of X values is unchanged, so the mid-rank vectors are fixed and only their
    # pairing with |S| moves.  Writing rx, ry for the fixed (19 x 2) mid-rank blocks, the whole
    # statistic is the assignment sum T(p) = sum_i A[p(i), i] with
    #     A[j, i] = rx[j,0]*ry[i,0] + rx[j,1]*ry[i,1],
    # and rho is an affine function of T.  Exhaustive enumeration would need 19! = 1.2e17 terms;
    # this evaluates B relabellings with two array lookups.
    def cluster_perm_p(Xp, Spp, obs, B, rng, chunk=200000):
        rx = stats.rankdata(Xp.ravel()).reshape(19, 2)
        ry = stats.rankdata(Spp.ravel()).reshape(19, 2)
        A = rx[:, 0:1] * ry[:, 0][None, :] + rx[:, 1:2] * ry[:, 1][None, :]   # (19 molecules, 19 slots)
        allr = stats.rankdata(Xp.ravel())
        mx, sx = allr.mean(), allr.std()
        ally = stats.rankdata(Spp.ravel())
        my, sy = ally.mean(), ally.std()
        n = 38.0
        to_rho = lambda T: (T - n * mx * my) / (n * sx * sy)
        if abs(to_rho(A[np.arange(19), np.arange(19)].sum()) - obs) > 1e-10:
            fail("vectorised cluster-permutation statistic disagrees with scipy on the identity")
        ge = 0
        cols = np.arange(19)
        done = 0
        while done < B:
            m = min(chunk, B - done)
            perms = np.argsort(rng.random((m, 19)), axis=1)
            T = A[perms, cols[None, :]].sum(axis=1)
            ge += int(np.count_nonzero(to_rho(T) >= obs - 1e-12))
            done += m
        return (ge + 1.0) / (B + 1.0), ge

    Bp = args.perm
    t0 = time.time()
    p_perm_chi, ge_chi = cluster_perm_p(Cp, Sp, rho_chi_m, Bp, rng)
    p_perm_F1, ge_F1 = cluster_perm_p(Fp, Sp, rho_F1_m, Bp, rng)
    log(f"  CLUSTER PERMUTATION (permute the 19 molecules, B = {Bp}, seed "
        f"{SEEDS['cluster_permutation_molecular']}, {time.time()-t0:.1f} s):")
    log(f"      rho(chi,|S|) = {rho_chi_m:.4f}   p = {p_perm_chi:.3e}   ({ge_chi} exceedances)")
    log(f"      rho(F1,|S|)  = {rho_F1_m:.4f}   p = {p_perm_F1:.3e}   ({ge_F1} exceedances)")
    log(f"      resolution floor of this test: 1/(B+1) = {1.0/(Bp+1):.2e}")

    rngb = np.random.default_rng(SEEDS["cluster_bootstrap_molecular"])
    Bb = args.boot
    codeC, KC = codes_of(Cp.ravel())
    codeF, KF = codes_of(Fp.ravel())
    codeS, KS = codes_of(Sp.ravel())
    codeC, codeF, codeS = codeC.reshape(19, 2), codeF.reshape(19, 2), codeS.reshape(19, 2)
    bc, bf = np.empty(Bb), np.empty(Bb)
    done, chunk = 0, 5000
    while done < Bb:
        m = min(chunk, Bb - done)
        pick = rngb.integers(0, 19, size=(m, 19))
        rS = batch_midranks(codeS[pick].reshape(m, 38), KS)
        bc[done:done + m] = row_rho(batch_midranks(codeC[pick].reshape(m, 38), KC), rS)
        bf[done:done + m] = row_rho(batch_midranks(codeF[pick].reshape(m, 38), KF), rS)
        done += m
    bd = bc - bf
    # self-check of the vectorised ranker against scipy on the unresampled data
    _chk = row_rho(batch_midranks(codeC.reshape(1, 38), KC),
                   batch_midranks(codeS.reshape(1, 38), KS))[0]
    if abs(_chk - rho_chi_m) > 1e-12:
        fail(f"vectorised batch ranker disagrees with scipy spearmanr ({_chk} vs {rho_chi_m})")
    ci = lambda v: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
    ci_c, ci_f, ci_d = ci(bc), ci(bf), ci(bd)
    log(f"  CLUSTER BOOTSTRAP (B = {Bb}, seed {SEEDS['cluster_bootstrap_molecular']}), "
        f"95% percentile CI:")
    log(f"      rho(chi,|S|) = {rho_chi_m:.4f}  CI [{ci_c[0]:.4f}, {ci_c[1]:.4f}]")
    log(f"      rho(F1,|S|)  = {rho_F1_m:.4f}  CI [{ci_f[0]:.4f}, {ci_f[1]:.4f}]")

    R["molecular"] = {
        "n_points": 38, "n_clusters": 19, "cluster_size": 2,
        "rho_chi_S": rho_chi_m, "rho_F1_S": rho_F1_m,
        "asymptotic_p_chi_as_published": p_asym_chi,
        "asymptotic_p_F1_as_published": p_asym_F1,
        "tie_structure_S": {str(int(u)): int(c) for u, c in zip(uniq, cnt)},
        "n_at_floor_S_eq_1": n_floor,
        "fraction_in_tied_group": n_tied / 38.0,
        "max_attainable_rho_given_S_ties": rho_ceiling,
        "rho_chi_S_as_fraction_of_ceiling": rho_chi_m / rho_ceiling,
        "rho_F1_S_as_fraction_of_ceiling": rho_F1_m / rho_ceiling,
        "ICC_log2S_within_molecule": float(icc),
        "design_effect": float(deff),
        "effective_n": float(n_eff),
        "cluster_permutation": {"B": Bp, "seed": SEEDS["cluster_permutation_molecular"],
                                "p_chi": p_perm_chi, "p_F1": p_perm_F1,
                                "exceedances_chi": ge_chi, "exceedances_F1": ge_F1,
                                "resolution_floor": 1.0 / (Bp + 1)},
        "cluster_bootstrap": {"B": Bb, "seed": SEEDS["cluster_bootstrap_molecular"],
                              "ci95_rho_chi_S": list(ci_c), "ci95_rho_F1_S": list(ci_f)},
        "molecule_order": mol_names,
    }

    # ------------------------------------------------------------------------------------------
    # BLOCK 3 -- THE COMPARISON THE PAPER NEVER MAKES
    # ------------------------------------------------------------------------------------------
    log("")
    log("=" * 94)
    log("BLOCK 3 -- formal comparison of the two DEPENDENT correlations (common variable |S|)")
    log("=" * 94)

    comps = {}

    r_xy_m = srho(Cm, Fm)
    tW, dfW, pW = williams_t(rho_chi_m, rho_F1_m, r_xy_m, 38)
    zS, pS = steiger_z(rho_chi_m, rho_F1_m, r_xy_m, 38)
    tW_eff, dfW_eff, pW_eff = williams_t(rho_chi_m, rho_F1_m, r_xy_m, 19)
    log(f"  (a) MOLECULAR, n=38 : rho(chi,|S|) = {rho_chi_m:.4f}  vs  rho(F1,|S|) = {rho_F1_m:.4f}")
    log(f"      rho(chi,F1) = {r_xy_m:.4f}   difference = {rho_chi_m - rho_F1_m:+.4f}")
    log(f"      Williams t = {tW:+.3f}  df = {dfW}   p = {pW:.4f}   (n=38, dependence IGNORED)")
    log(f"      Steiger  z = {zS:+.3f}            p = {pS:.4f}")
    log(f"      Williams at the EFFECTIVE n = 19 : t = {tW_eff:+.3f}  df = {dfW_eff}  "
        f"p = {pW_eff:.4f}")
    log(f"      cluster-bootstrap 95% CI on the difference: [{ci_d[0]:+.4f}, {ci_d[1]:+.4f}]"
        f"   contains 0: {ci_d[0] <= 0.0 <= ci_d[1]}")
    resolvable_a = (pW < 0.05) and (pW_eff < 0.05) and not (ci_d[0] <= 0.0 <= ci_d[1])
    log(f"      VERDICT: the 0.90-vs-0.86 gap is "
        f"{'RESOLVABLE' if resolvable_a else 'NOT RESOLVABLE'} at these sample sizes.")
    comps["molecular_38"] = {"r_chi_S": rho_chi_m, "r_F1_S": rho_F1_m, "r_chi_F1": r_xy_m,
                             "difference": rho_chi_m - rho_F1_m,
                             "williams_t": tW, "williams_df": dfW, "williams_p": pW,
                             "steiger_z": zS, "steiger_p": pS,
                             "williams_at_effective_n19": {"t": tW_eff, "df": dfW_eff,
                                                           "p": pW_eff},
                             "cluster_bootstrap_ci95_difference": list(ci_d),
                             "resolvable": bool(resolvable_a)}

    r_chi_p, r_F1_p = got["pooled_rho_chi_S"], got["pooled_rho_F1_S"]
    r_xy_p = srho(CH_all, F1_all)
    tWp, dfWp, pWp = williams_t(r_chi_p, r_F1_p, r_xy_p, 74)
    zSp, pSp = steiger_z(r_chi_p, r_F1_p, r_xy_p, 74)
    rngf = np.random.default_rng(SEEDS["family_bootstrap_pooled"])
    fams = {f: np.flatnonzero(fam == f) for f in sorted(set(fam.tolist()))}
    gC, KgC = codes_of(CH_all)
    gF, KgF = codes_of(F1_all)
    gS, KgS = codes_of(S_all)
    bdp = np.empty(args.boot)
    done, chunk = 0, 5000
    while done < args.boot:
        m = min(chunk, args.boot - done)
        sel = np.concatenate([rngf.choice(v, size=(m, v.size), replace=True)
                              for v in fams.values()], axis=1)
        rS = batch_midranks(gS[sel], KgS)
        bdp[done:done + m] = (row_rho(batch_midranks(gC[sel], KgC), rS)
                              - row_rho(batch_midranks(gF[sel], KgF), rS))
        done += m
    ci_dp = (float(np.percentile(bdp, 2.5)), float(np.percentile(bdp, 97.5)))
    log("")
    log(f"  (b) POOLED, n=74 : rho(chi,|S|) = {r_chi_p:.4f}  vs  rho(F1,|S|) = {r_F1_p:.4f}")
    log(f"      rho(chi,F1) = {r_xy_p:.4f}   difference = {r_chi_p - r_F1_p:+.4f}")
    log(f"      Williams t = {tWp:+.3f}  df = {dfWp}   p = {pWp:.5f}")
    log(f"      Steiger  z = {zSp:+.3f}            p = {pSp:.5f}")
    log(f"      within-family bootstrap 95% CI on the difference: "
        f"[{ci_dp[0]:+.4f}, {ci_dp[1]:+.4f}]   contains 0: {ci_dp[0] <= 0.0 <= ci_dp[1]}")
    resolvable_b = (pWp < 0.05) and not (ci_dp[0] <= 0.0 <= ci_dp[1])
    log(f"      VERDICT: the 0.72-vs-0.60 gap is "
        f"{'RESOLVABLE' if resolvable_b else 'NOT RESOLVABLE'}.")
    log("      CAVEAT: the 74 points mix four families with different designs; the pooled")
    log("      coefficient is a descriptive summary, not an estimate of one population parameter.")
    comps["pooled_74"] = {"r_chi_S": r_chi_p, "r_F1_S": r_F1_p, "r_chi_F1": r_xy_p,
                          "difference": r_chi_p - r_F1_p,
                          "williams_t": tWp, "williams_df": dfWp, "williams_p": pWp,
                          "steiger_z": zSp, "steiger_p": pSp,
                          "within_family_bootstrap_ci95_difference": list(ci_dp),
                          "resolvable": bool(resolvable_b)}

    Fh, Sh, Ch = F1_all[hmask], S_all[hmask], CH_all[hmask]
    r_xy_h = srho(Ch, Fh)
    tWh, dfWh, pWh = williams_t(got["hubbard_rho_chi_S"], got["hubbard_rho_F1_S"], r_xy_h, 30)
    tWh6, dfWh6, pWh6 = williams_t(got["hubbard_rho_chi_S"], got["hubbard_rho_F1_S"], r_xy_h, 9)
    log("")
    log(f"  (c) HUBBARD, n=30 : rho(chi,|S|) = {got['hubbard_rho_chi_S']:.4f}  vs  "
        f"rho(F1,|S|) = {got['hubbard_rho_F1_S']:.4f}")
    log(f"      Williams t = {tWh:+.3f}  df = {dfWh}   p = {pWh:.3e}   (n=30, strata IGNORED)")
    log(f"      at n = 9 (a deliberately conservative stand-in for the 6 strata): "
        f"t = {tWh6:+.3f}  df = {dfWh6}  p = {pWh6:.4f}")
    log("      NOTE: this contrast is the one the paper leans on, and it survives even the")
    log("      conservative reading -- but the correct statement is the STRATIFIED one of Block 1,")
    log("      where F1 is not 'uncorrelated' but PERFECTLY ANTI-correlated.")
    comps["hubbard_30"] = {"r_chi_S": got["hubbard_rho_chi_S"],
                           "r_F1_S": got["hubbard_rho_F1_S"], "r_chi_F1": r_xy_h,
                           "williams_t": tWh, "williams_df": dfWh, "williams_p": pWh,
                           "williams_conservative_n9": {"t": tWh6, "df": dfWh6, "p": pWh6}}
    R["dependent_correlation_comparisons"] = comps

    # ------------------------------------------------------------------------------------------
    # BLOCK 4 -- chain2 AND ladder (n = 3 each)
    # ------------------------------------------------------------------------------------------
    log("")
    log("=" * 94)
    log("BLOCK 4 -- chain2 and ladder (n=3 each): F1 outranks chi in 2 of the 4 families")
    log("=" * 94)

    small = {}
    n3_null = stratum_null(np.arange(3.0), np.arange(3.0))
    p_floor_one = exact_one_sided_p(n3_null, 1.0, "right")
    for f in ("chain2", "ladder"):
        sel = fam == f
        if sel.sum() != 3:
            fail(f"family {f} has {int(sel.sum())} points, expected 3")
        Ff, Sf, Cf = F1_all[sel], S_all[sel], CH_all[sel]
        o = np.argsort(Sf)
        Ff, Sf, Cf = Ff[o], Sf[o], Cf[o]
        rf, rc = srho(Ff, Sf), srho(Cf, Sf)
        pf = exact_one_sided_p(stratum_null(Sf, Ff), rf, "right")
        pc = exact_one_sided_p(stratum_null(Sf, Cf), rc, "right")
        Lgrid = np.array([6.0, 8.0, 10.0])          # the ladder_vs_chain design: U=8, L=6,8,10
        f_over_L = Ff / Lgrid
        log(f"  {f:7s}: F1 = {np.round(Ff,4).tolist()}   |S| = {Sf.astype(int).tolist()}   "
            f"chi = {Cf.astype(int).tolist()}")
        log(f"           rho(F1,|S|) = {rf:+.4f}  exact p_one = {pf:.4f}"
            f"      rho(chi,|S|) = {rc:+.4f}  exact p_one = {pc:.4f}")
        log(f"           F1/L = {np.round(f_over_L,4).tolist()}   spread "
            f"{f_over_L.max()-f_over_L.min():.4f}   (F1 = 1.60*L to three digits)")
        small[f] = {"F1": Ff.tolist(), "S": Sf.astype(int).tolist(),
                    "chi": Cf.astype(int).tolist(), "L": Lgrid.tolist(),
                    "F1_over_L": f_over_L.tolist(),
                    "rho_F1_S": rf, "exact_p_one_sided_F1": pf,
                    "rho_chi_S": rc, "exact_p_one_sided_chi": pc}
    chi2s, dfs, ps = fisher_combine([small["chain2"]["exact_p_one_sided_F1"],
                                     small["ladder"]["exact_p_one_sided_F1"]])
    log("")
    log(f"  EXACT FLOOR AT n=3: the smallest attainable one-sided p is P(rho=+1) = "
        f"{p_floor_one:.4f} = 1/6.")
    log(f"  No n=3 family can reach p < 0.05.  rho = +1.000 here carries log2(6) = "
        f"{math.log2(6):.2f} bits -- less evidence than three coin flips.")
    log(f"  Fisher over the two families: X2 = {chi2s:.3f}  df = {dfs}  p = {ps:.4f}  "
        f"-- and even this is inadmissible:")
    log("     (i)   the direction is chosen post hoc, after seeing rho = +1;")
    log("     (ii)  the two families are NOT independent: same U = 8, same L grid {6,8,10}, and")
    log("           all three chain2 rows are duplicates of three of the 30 Hubbard rows;")
    log("     (iii) THE CONFOUND: within these families the swept parameter is L, not U.  F1 =")
    log("           1.60*L to three digits while |S| ~ exp(L), so rho(F1,|S|) = +1 restates")
    log("           'both grow with system size' and carries no predictive content.  After")
    log("           conditioning on L, n=3 over a 3-point L grid leaves ZERO residual degrees of")
    log("           freedom: the partial correlation is undefined, not small.")
    log("  VERDICT: arithmetically true, statistically uninterpretable.  The manuscript should")
    log("  state it and disarm it rather than omit it -- an omitted +1.000 is what a referee finds.")
    R["small_families"] = {**small,
                           "exact_p_floor_n3_one_sided": p_floor_one,
                           "fisher_two_families": {"chi2": chi2s, "df": dfs, "p": ps},
                           "verdict": ("arithmetically true, statistically uninterpretable: n=3 "
                                       "floor p=1/6, families not independent, and the sweep "
                                       "variable is L with F1 = 1.60*L, so F1 is a proxy for "
                                       "system size.")}

    # ------------------------------------------------------------------------------------------
    # BLOCK 5 -- GFLOWNET
    # ------------------------------------------------------------------------------------------
    log("")
    log("=" * 94)
    log("BLOCK 5 -- GFlowNet: the correct paired test, and the deposit gap that blocks it")
    log("=" * 94)

    meth = {m["label"]: m for m in gflow["methods"]}
    nseeds = int(gflow["nseeds"])
    has_raw = any(("per_seed" in m) or ("seeds" in m) for m in gflow["methods"])
    log(f"  data/gflow.json deposits {len(gflow['methods'])} SUMMARY rows (mean, std) over "
        f"{nseeds} seeds.")
    log(f"  per-seed raw present in the deposit: {has_raw}   "
        f'(the json note says: "Per-seed raw available on request")')
    log("  => THE PAIRED t IS NOT COMPUTABLE FROM THE DEPOSIT.  That is a deposit defect and it")
    log('     contradicts hardware.tex\'s "complete reproduction repository".  What follows is the')
    log("     rigorous INTERVAL of paired t consistent with the deposited summaries.")

    def paired_bounds(a, b, n):
        """For paired samples with marginal s.d. sa, sb, the s.d. of the differences obeys
        |sa - sb| <= sd <= sa + sb (the pairing correlation lies in [-1,1]).  Returns the induced
        interval of the paired t statistic and its two-sided p."""
        d = a["mean_mHa"] - b["mean_mHa"]
        sa, sb = a["std_mHa"], b["std_mHa"]
        sd_lo, sd_hi = abs(sa - sb), sa + sb
        t_hi = d / (sd_lo / math.sqrt(n)) if sd_lo > 0 else float("inf")
        t_lo = d / (sd_hi / math.sqrt(n))
        pl = float(2.0 * stats.t.sf(abs(t_lo), n - 1))
        ph = float(2.0 * stats.t.sf(abs(t_hi), n - 1)) if math.isfinite(t_hi) else 0.0
        return {"delta_mHa": d, "sd_diff_lower": sd_lo, "sd_diff_upper": sd_hi,
                "t_worst_case": t_lo, "p_worst_case_two_sided": pl,
                "t_best_case": t_hi, "p_best_case_two_sided": ph, "df": n - 1}

    tcrit = float(stats.t.ppf(0.975, nseeds - 1))
    cmpA = paired_bounds(meth["ibm"], meth["ibm+cheap"], nseeds)
    cmpB = paired_bounds(meth["ibm+cheap"], meth["gfn(fused)"], nseeds)
    log("")
    log(f"  two-sided t critical value at df = {nseeds-1}: {tcrit:.3f}")
    verdicts = {}
    for tag, c, desc in (("A", cmpA, "ibm -> ibm+cheap        (the FREE classical tail)"),
                         ("B", cmpB, "ibm+cheap -> gfn(fused) (the LEARNED generator)")):
        dec = ("SIGNIFICANT for every admissible pairing" if c["t_worst_case"] > tcrit
               else ("NOT significant for any admissible pairing" if c["t_best_case"] < tcrit
                     else "UNDECIDABLE from the deposit"))
        verdicts[tag] = dec
        log(f"  contrast {tag}: {desc}")
        log(f"      delta = {c['delta_mHa']:.2f} mHa;  sd(diff) in "
            f"[{c['sd_diff_lower']:.2f}, {c['sd_diff_upper']:.2f}]")
        log(f"      paired t in [{c['t_worst_case']:.3f}, {c['t_best_case']:.3f}]  ->  p in "
            f"[{c['p_best_case_two_sided']:.4f}, {c['p_worst_case_two_sided']:.4f}]")
        log(f"      => {dec}")

    dA = meth["ibm"]["mean_mHa"] - meth["ibm+cheap"]["mean_mHa"]
    dB = meth["ibm+cheap"]["mean_mHa"] - meth["gfn(fused)"]["mean_mHa"]
    conv = {
        "paper_as_printed": (dA / meth["ibm"]["std_mHa"], dB / meth["gfn(fused)"]["std_mHa"]),
        "baseline_sd_both": (dA / meth["ibm"]["std_mHa"], dB / meth["ibm+cheap"]["std_mHa"]),
        "target_sd_both": (dA / meth["ibm+cheap"]["std_mHa"], dB / meth["gfn(fused)"]["std_mHa"]),
        "pooled_sd_both": (
            dA / math.sqrt(0.5 * (meth["ibm"]["std_mHa"] ** 2
                                  + meth["ibm+cheap"]["std_mHa"] ** 2)),
            dB / math.sqrt(0.5 * (meth["ibm+cheap"]["std_mHa"] ** 2
                                  + meth["gfn(fused)"]["std_mHa"] ** 2)),
        ),
    }
    log("")
    log("  THE DENOMINATOR INCONSISTENCY (s.d. units for contrast A and contrast B):")
    for k, (a, b) in conv.items():
        log(f"      {k:20s} :  A = {a:6.2f} s.d.   B = {b:6.2f} s.d.")
    log("  The manuscript's '3.1 s.d.' divides by the BASELINE dispersion; its '1.2 s.d.' divides")
    log("  by the NEW method's dispersion.  Applied consistently in the baseline convention the")
    log(f"  second contrast reads {conv['baseline_sd_both'][1]:.2f} s.d., not 1.2.  The")
    log("  inconsistency runs AGAINST the author's own result.  None of the four conventions is")
    log("  the right test; the paired t is, and it needs the per-seed raw.")

    p_sign_floor = 2.0 / (2 ** nseeds)
    need = int(math.ceil(math.log2(2.0 / 0.05)))
    log(f"  NONPARAMETRIC FLOOR: with {nseeds} paired seeds the exact sign-flip permutation test")
    log(f"  has a minimum two-sided p of 2/2^{nseeds} = {p_sign_floor:.4f} > 0.05.  No 5-seed")
    log(f"  paired experiment can be declared significant distribution-free.  {need} seeds put the")
    log("  floor below 0.05.  Report effect sizes with intervals, or run more seeds.")

    R["gflownet"] = {
        "per_seed_raw_in_deposit": bool(has_raw),
        "n_seeds": nseeds,
        "t_crit_two_sided_0p05": tcrit,
        "contrast_A_ibm_to_ibmcheap": cmpA,
        "contrast_B_ibmcheap_to_gfnfused": cmpB,
        "verdict_A": verdicts["A"],
        "verdict_B": verdicts["B"],
        "sd_unit_conventions": {k: list(v) for k, v in conv.items()},
        "nonparametric_min_two_sided_p": p_sign_floor,
        "seeds_needed_for_nonparametric_floor_below_0p05": need,
    }

    # ------------------------------------------------------------------------------------------
    # BLOCK 6 -- FORENSICS of the orphan pooled coefficient
    # ------------------------------------------------------------------------------------------
    log("")
    log("=" * 94)
    log("BLOCK 6 -- forensics: the orphan constant pooled_spearman_F1_S in resource_master.json")
    log("=" * 94)

    dep_F1 = rmaster["pooled_spearman_F1_S"]
    dep_chi = rmaster["pooled_spearman_chi_S"]
    rec_F1 = srho(F1_all, S_all)
    rec_chi = srho(CH_all, S_all)
    log(f"  deposited pooled_spearman_chi_S = {dep_chi!r}")
    log(f"  recomputed                      = {rec_chi!r}   -> reproduces exactly")
    log(f"  deposited pooled_spearman_F1_S  = {dep_F1!r}")
    log(f"  recomputed                      = {rec_F1!r}   -> does NOT reproduce "
        f"(|d| = {abs(dep_F1 - rec_F1):.3e})")

    tests = {}
    tgt = np.isclose(F1_all, 9.604689, atol=1e-9)
    hub_idx = np.flatnonzero(tgt & (fam == "hubbard"))
    ch_idx = np.flatnonzero(tgt & (fam == "chain2"))
    if hub_idx.size != 1 or ch_idx.size != 1:
        fail("expected exactly one hubbard and one chain2 row at F1 = 9.604689")
    f_break = F1_all.copy()
    f_break[ch_idx[0]] += 1e-9                        # break the tie, preserve the rank order
    r_break = srho(f_break, S_all)
    tests["tie_broken_chain2_vs_hubbard_at_9.604689"] = r_break
    f_pbc = F1_all.copy()
    f_pbc[ch_idx[0]] = 9.385082                       # the pre-fix PBC value
    tests["chain2_L6_set_to_PBC_9.385082"] = srho(f_pbc, S_all)
    f_hub = F1_all.copy()
    f_hub[hub_idx[0]] = 9.3851                        # the pre-fix committed value
    tests["hubbard_L6N6U8_set_to_committed_9.3851"] = srho(f_hub, S_all)
    f_both = F1_all.copy()
    f_both[[hub_idx[0], ch_idx[0]]] = 9.3851          # both moved: the tie is PRESERVED
    tests["both_rows_set_to_9.3851_tie_preserved"] = srho(f_both, S_all)

    log("")
    log("  HYPOTHESIS: the deposited constant predates the 2026-09-05 provenance fix, which set")
    log("  the Hubbard (L=6,N=6,U=8) row and the chain2 L=6 row to the same F1 = 9.604689 and")
    log("  thereby CREATED an F1 tie that did not previously exist.  Test -- re-break that tie:")
    for k, v in tests.items():
        mark = "   <== EXACT MATCH to the deposited constant" if abs(v - dep_F1) < 1e-15 else ""
        log(f"      {k:46s} -> {v!r}{mark}")
    explained = abs(r_break - dep_F1) < 1e-15
    log("")
    if explained:
        log("  RESOLVED.  The orphan constant is neither irreproducible nor wrong: it is the pooled")
        log("  Spearman coefficient computed with the F1 tie between the Hubbard L=6 row and the")
        log("  chain2 L=6 row ABSENT.  Any tie-breaking perturbation of either row reproduces it to")
        log("  machine precision.  The provenance fix restored both rows to 9.604689, created the")
        log("  tie, and moved the coefficient to 0.5992028038362396 -- but the stored constant was")
        log("  never recomputed.  chi is unaffected because chi = 5 in both rows already.")
        log("")
        log("  CONSEQUENCE: the provenance note in resource_master.json states 'Rank-preserving:")
        log("  no F1 value in the 74 lies between the two, so every published Spearman coefficient")
        log("  is unchanged (... F1-|S| 0.5992, verified after the edit)'.  That claim is FALSE as")
        log("  stated: the edit preserved rank ORDER but not TIE STRUCTURE, and the coefficient did")
        log("  change, in the 6th decimal (0.5991983669 -> 0.5992028038).  Numerically irrelevant")
        log("  -- both print as 0.60 -- but the audit trail is wrong, and a referee who recomputes")
        log("  the deposit finds a number that does not reproduce.")
        log("  FIX: set pooled_spearman_F1_S to 0.5992028038362396 and amend the note.")
    else:
        log("  NOT RESOLVED by the tie hypothesis.  Declare the constant irreproducible and replace")
        log("  it with the recomputed value.")

    R["orphan_constant_forensics"] = {
        "deposited_pooled_spearman_F1_S": dep_F1,
        "recomputed_pooled_spearman_F1_S": rec_F1,
        "abs_diff": abs(dep_F1 - rec_F1),
        "deposited_pooled_spearman_chi_S": dep_chi,
        "recomputed_pooled_spearman_chi_S": rec_chi,
        "chi_reproduces_exactly": bool(abs(dep_chi - rec_chi) < 1e-15),
        "hypothesis_tests": tests,
        "explained_by_F1_tie_creation": bool(explained),
        "recommended_value": rec_F1,
    }

    # ------------------------------------------------------------------------------------------
    # PROVENANCE + WRITE
    # ------------------------------------------------------------------------------------------
    elapsed = time.time() - t_start
    R["provenance"] = {
        "script": "stats_resource.py",
        "purpose": "correct statistical re-analysis of the resource section of the manuscript",
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": elapsed,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "seeds": SEEDS,
        "parameters": {"boot": args.boot, "perm": args.perm,
                       "root": _repo_rel(args.root),
                       "control_tolerance": CONTROL_TOL, "U_grid": list(U_GRID)},
        "inputs": ["data/cost_vs_ent.json", "data/resource_master.json",
                   "data/n19_suite.json", "data/gflow.json"],
        "control_values_reproduced": {k: {"paper": CONTROLS[k], "recomputed": got[k]}
                                      for k in CONTROLS},
        "determinism": ("all randomness from numpy.random.default_rng with the seeds above; "
                        "bit-for-bit reproducible on numpy " + np.__version__),
    }

    # FIXED 2026-09-18: a RELATIVE --out used to be joined onto --root, and --root is
    # the repository. A verification run therefore wrote its scratch output INTO the
    # deposit (it left release/bite/stats_rerun.json behind). --root says where to READ
    # the data from; it must never decide where results are WRITTEN. A relative --out
    # now resolves against the current working directory, like every other tool.
    outp = os.path.abspath(args.out)
    os.makedirs(os.path.dirname(os.path.abspath(outp)), exist_ok=True)
    with open(outp, "w", encoding="utf-8") as f:
        json.dump(R, f, indent=1, sort_keys=False)
    log("")
    log("=" * 94)
    log(f"wrote {outp}   ({elapsed:.1f} s)")
    log("=" * 94)


if __name__ == "__main__":
    main()
