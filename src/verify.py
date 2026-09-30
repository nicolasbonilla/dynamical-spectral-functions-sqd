# -*- coding: utf-8 -*-
r"""verify.py -- adversarial reproducibility guardian for this repository.

WHAT THIS IS
------------
A guardian that RECOMPUTES the load-bearing numbers of the manuscript from first
principles and then OPENS the deposited artefacts -- the JSON datasets, the .dat
tables that feed the figures, the native pgfplots figure sources, and the printed
numbers in paper/*.tex -- and checks that every one of them carries the value the
recomputation produces.  It exits non-zero, naming the defect, when any of them
does not.

WHY IT LOOKS LIKE THIS
----------------------
The previous version of this file did not open a single file of the repository.
Its geminal block compared the literal expression ``4*K`` against ``4.0*K``; eight
of its ten [PASS] lines were ``4*K == 4*K``.  It never built the geminal state,
never formed gamma, never counted a support and never computed chi, although its
docstring promised chi.  Consequences, measured: re-introducing the exact 2026-08-17
defect (the open-chain F1 9.604689 overwritten by the periodic-ring value 9.3851 in
data/resource_master.json, data/ladder_vs_chain.json,
paper/figs/fig_resource_master_native.tex and paper/results.tex) left it printing
``ALL PASS``; so did destroying the geminal witness (|S|=1, chi=99).

Three rules follow from that history and they shape everything below.

 1. TRUTH IS RECOMPUTED, NEVER STORED.  No headline value is hard-coded as an
    expected constant.  F1 for the open chain and for the periodic ring, the
    geminal invariants, gamma, the support, chi, the Schmidt spectra, the
    correlations and the spectral edges are all computed here, and the deposited
    artefacts are then compared against that computation.  What IS hard-coded is
    the assignment of each artefact slot to a physical system -- which is exactly
    the statement a guardian is supposed to make.

 2. OBC AND PBC ARE DIFFERENT SYSTEMS.  F1(L=6,U=8) is 9.604689 for the OPEN chain
    and 9.385082 for the PERIODIC ring; 12.816276 / 12.813821 at L=8; 16.030613 /
    16.018265 at L=10.  These are two systems, not two versions of one number.  The
    2026-08-17 edit substituted one for the other and was reverted on 2026-09-05.
    Every slot below is declared OBC or PBC and the check is two-sided: the value
    must match its own system AND must NOT match the other one.  In particular
    paper/figs/fig_scaling2_native.tex legitimately carries 9.3851 because its
    engine is the ring (data/scaling_data.json, akw_lanczos.hop, j=(i+1)%L);
    "correcting" it to 9.6047 is a NEW defect and this guardian fails on it.

 3. OPEN DEFECTS ARE NAMED, NOT HIDDEN.  A handful of numbers in the tree are known
    to be wrong and are not this file's to fix.  They are listed in KNOWN_OPEN below
    with the wrong value, the correct value and the reason.  Each is reported as
    [XFAIL] and does not set the exit code -- but the registry PINS the wrong value:
    if such a slot changes to any third value it is a hard [FAIL], and if it is
    repaired it is reported as [XPASS] with an instruction to delete the entry.  A
    registry entry can therefore excuse exactly one known number and nothing else.

RUN
---
    python src/verify.py                 full run   (~15 s; includes L=10 ED)
    python src/verify.py --fast          quick run  (~3 s; skips the L=10 ED)
    python src/verify.py --root PATH     check a tree other than the parent of src/
    python src/verify.py -v              also print the PASS lines

Dependencies: numpy + scipy only.  Nothing here imports pyscf, so the molecular
checks read the deposited data and the geometry table (parsed with `ast`, not
imported) rather than re-running quantum chemistry.

Exit status: 0 iff no hard check failed.  Any hard failure prints a line that names
the defect -- which slot, which system, which value was found and which was due.
"""

import argparse
import ast
import hashlib
import json
import math
import os
import re
import sys
import time
from itertools import combinations

import numpy as np
# --- NUMPY_TRAPEZOID_BRIDGE ------------------------------------------------
# numpy 2.0 ADDED np.trapezoid and REMOVED np.trapz.  Files in this repository use
# both names, so without this bridge no single numpy version runs the whole deposit
# -- this guardian calls np.trapz and would itself die on numpy 2.x.
if not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz          # numpy < 2.0
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid          # numpy >= 2.0
# ---------------------------------------------------------------------------
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.stats import pearsonr, spearmanr

# ---------------------------------------------------------------------------
# Registry of KNOWN OPEN DEFECTS.
# Each entry pins the wrong value that is currently in the tree.  Nothing else is
# excused: a third value in the same slot is a hard failure.
# ---------------------------------------------------------------------------
KNOWN_OPEN = {
    "sqw_lanczos.eta_default": dict(
        wrong=0.20, right=0.18, tol=1e-9,
        why="data/sqw_L12.json (and the S(q,w) figure built on it) was computed at eta=0.18, "
            "but src/sqw_lanczos.py:21 now defaults to eta=0.20, so running the deposited "
            "script as committed does not regenerate the deposited spectrum."),
    "resource_master.pooled_spearman_F1_S": dict(
        wrong=0.5991983669238197, right=0.5992028038362396, tol=1e-12,
        why="the stored pooled Spearman F1-|S| is not reproducible from the 74 deposited "
            "points with either set of F1 values; recomputation gives 0.5992028038362396. "
            "Both round to the 0.60 printed in resource.tex, so no printed number moves."),
    # "tab_repro.H2O.R_eq" was here until 2026-09-19.  paper/app_carried_repro.tex
    # printed R_eq(H2O) = 0.957 A while the geometry in src/n19_suite.py gives
    # sqrt(0.757^2 + 0.587^2) = 0.95792 -> 0.958.  The manuscript was corrected, this
    # guardian reported XPASS and asked for the entry to go, and it is gone: the slot
    # is now a HARD check that fails if the printed value drifts again.
    "method_max.sumrule": dict(
        wrong=0.5000000000000002, right=0.48658771750884894, tol=1e-9,
        why="data/method_max.json stores <Psi0|c c^dag|Psi0> = 1-<n_p>, an analytic "
            "ground-state identity, under the name 'sumrule'.  The integral of the "
            "deposited A_exact over the deposited window is 0.48659 (the reconstruction "
            "gives 0.48647): a finite-window effect, not machine precision.  The "
            "gap is 2.7 per cent of the sum rule, a finite-window effect.  The manuscript "
            "half of this defect is CLOSED: the figure now labels the windowed integral "
            "0.4865 and no caption says 'to machine precision' (see CLAIMS_RETRACTED). "
            "What remains open is the DATA file, which still stores the analytic identity "
            "under the name 'sumrule'."),
    "ladder_L4.ngrid_eta050": dict(
        wrong=1497.0, right=2265.0, tol=0.5, pin_tol=0.5,
        why="data/c3_frontier/ladder_L4.json was computed at 17:32 on 2026-09-18, four "
            "minutes before src/frontier/fron_lib.make_grid widened its integration "
            "window from pad=4.0 to pad=max(4, 40*eta).  Re-running "
            "`PER_ETA=12 python src/frontier/run_ladder.py 4 0.2,0.4,0.6,0.8,0.9` today "
            "reproduces every other field of that file but integrates over a wider "
            "window, giving ngrid=2265 at eta=0.50.  Measured consequence: only the L=4 "
            "rows of verdict_frontier.json move, by at most 2.8 per cent relative on d_CS "
            "(7e-4 absolute in fraction units); every L=6 and L=8 row is unchanged.  The "
            "L=6 and L=8 ladders postdate the change and reproduce exactly at PER_ETA=16."),
}

# Historical note, deliberately left here: until 2026-09-18 this registry also carried
# `sqw_edges.fabricated_rows` -- two rows of paper/figs/sqw_edges.dat (q/pi = 0.0000 and
# 2.0000) that no computation in this repository produced.  q=0 is inaccessible by
# construction (sqw_lanczos.py 'skip q=0'; sqw_field.py 'S(q=0,w>0)=0'), and the value
# 5.000 in them is not even on the frequency grid.  That table has since been regenerated
# from data/sqw_L12.json and all eleven surviving rows now reconstruct byte-for-byte, so
# the entry is gone and any row without a computational origin is a hard failure again.

# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

PASS, FAIL, XFAIL, XPASS, SKIP = "PASS", "FAIL", "XFAIL", "XPASS", "SKIP"


class Report(object):
    def __init__(self, verbose=False):
        self.verbose = verbose
        self.counts = {PASS: 0, FAIL: 0, XFAIL: 0, XPASS: 0, SKIP: 0}
        self.failures = []
        self.xfails = []
        self.xpasses = []
        self.skips = []
        self.numbers = 0          # scalar comparisons actually performed
        self.section = ""
        self.per_section = []

    def head(self, title):
        self.section = title
        self.per_section.append([title, 0, 0])
        print("\n" + title)
        print("-" * len(title))

    def emit(self, status, name, msg, n=1):
        self.counts[status] += 1
        if status in (PASS, FAIL, XFAIL, XPASS):
            self.numbers += n
            if self.per_section:
                self.per_section[-1][1] += n
                self.per_section[-1][2] += 1
        line = "  [%-5s] %s%s" % (status, name, (" -- " + msg) if msg else "")
        if status == FAIL:
            self.failures.append((self.section, name, msg))
            print(line)
        elif status == XFAIL:
            self.xfails.append((self.section, name, msg))
            if self.verbose:
                print(line)
        elif status == XPASS:
            self.xpasses.append((self.section, name, msg))
            print(line)
        elif status == SKIP:
            self.skips.append((self.section, name, msg))
            print(line)
        elif self.verbose:
            print(line)

    # -- primitive comparisons ------------------------------------------------
    def eq(self, name, got, want, tol=0.0, defect=None, n=1):
        """Hard numeric check.  `defect` names what is wrong if it fails."""
        try:
            good = abs(float(got) - float(want)) <= tol
        except (TypeError, ValueError):
            good = False
        if good:
            self.emit(PASS, name, "%.10g == %.10g" % (float(got), float(want)), n)
        else:
            d = defect or "value does not match the recomputed/derived quantity"
            self.emit(FAIL, name, "%s: found %r, due %.10g (tol %g)" % (d, got, float(want), tol), n)
        return good

    def eq_int(self, name, got, want, defect=None, n=1):
        good = (got == want)
        if good:
            self.emit(PASS, name, "%r == %r" % (got, want), n)
        else:
            d = defect or "value does not match the recomputed/derived quantity"
            self.emit(FAIL, name, "%s: found %r, due %r" % (d, got, want), n)
        return good

    def truth(self, name, cond, msg_ok="", msg_bad="", n=1):
        if cond:
            self.emit(PASS, name, msg_ok, n)
        else:
            self.emit(FAIL, name, msg_bad, n)
        return bool(cond)

    def known(self, key, name, got, n=1):
        """Check a slot listed in KNOWN_OPEN.  Pins the wrong value.

        The wrong value is pinned at `pin_tol` (1e-12 by default, i.e. exactly as the
        literal is written); only the REPAIRED branch uses the looser `tol`.  With a
        single shared tolerance a third value -- neither the defect nor the repair --
        could land inside the repaired window and be reported as fixed.  It cannot now.
        """
        e = KNOWN_OPEN[key]
        g = float(got)
        pin = float(e.get("pin_tol", 1e-12))
        if abs(g - float(e["wrong"])) <= pin:
            self.emit(XFAIL, name, "known open defect: %.10g present, %.10g correct. %s"
                      % (g, float(e["right"]), e["why"]), n)
            return
        if abs(g - float(e["right"])) <= e["tol"]:
            self.emit(XPASS, name,
                      "REPAIRED (now %.10g). Delete KNOWN_OPEN['%s'] so this becomes a hard check."
                      % (g, key), n)
        else:
            self.emit(FAIL, name,
                      "UNREGISTERED value in a slot with a known defect: found %.10g; the "
                      "registry pins %.10g (wrong) and %.10g (correct). %s"
                      % (g, float(e["wrong"]), float(e["right"]), e["why"]), n)

    def missing(self, name, msg):
        self.emit(FAIL, name, msg, n=0)

    def skip(self, name, msg):
        self.emit(SKIP, name, msg, n=0)


# ---------------------------------------------------------------------------
# Section 1a -- Hubbard ground states, recomputed.  OBC and PBC.
# ---------------------------------------------------------------------------

def strings(L, n):
    S = sorted(sum(1 << i for i in c) for c in combinations(range(L), n))
    return S, {m: a for a, m in enumerate(S)}


def hop_matrix(L, n, periodic, t=1.0):
    """Spin-resolved hopping in the fixed-particle-number determinant basis.

    periodic=False -> OPEN chain, bonds (0,1)...(L-2,L-1)  [the OBC family]
    periodic=True  -> PERIODIC ring, bond (L-1,0) closes   [the PBC family;
                      identical to akw_lanczos.hop, j=(i+1)%L]
    """
    S, idx = strings(L, n)
    bonds = [(i, (i + 1) % L) for i in range(L)] if periodic else [(i, i + 1) for i in range(L - 1)]
    r, c, v = [], [], []
    for a, m in enumerate(S):
        for (i, j) in bonds:
            for (p, q) in ((i, j), (j, i)):
                if (m >> q) & 1 and not (m >> p) & 1:
                    m2 = (m & ~(1 << q)) | (1 << p)
                    lo, hi = min(p, q), max(p, q)
                    mask = m & (((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))
                    sg = -1.0 if (bin(mask).count("1") & 1) else 1.0
                    r.append(idx[m2]); c.append(a); v.append(-t * sg)
    return sp.csr_matrix((v, (r, c)), shape=(len(S), len(S))), S, idx


def onerdm_spin(rho, S, idx, L):
    """Spin-resolved 1-RDM G[i,j] = <c^dag_i c_j> from the reduced density matrix
    `rho` of one spin species in the determinant basis."""
    G = np.zeros((L, L))
    for a, m in enumerate(S):
        for j in range(L):
            if not (m >> j) & 1:
                continue
            for i in range(L):
                if i == j:
                    m2, sg = m, 1.0
                elif (m >> i) & 1:
                    continue
                else:
                    m2 = (m & ~(1 << j)) | (1 << i)
                    lo, hi = min(i, j), max(i, j)
                    mask = m & (((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))
                    sg = -1.0 if (bin(mask).count("1") & 1) else 1.0
                G[i, j] += float(rho[idx[m2], a].real) * sg
    return G


_HUB_CACHE = {}


def hubbard_f1(L, U, periodic, t=1.0):
    """Ground state of the half-filled Hubbard model; returns (F1, E0).

    F1 = 4 tr[gamma(1-gamma)] over spin-orbitals -- the one-body fermionic magic,
    equal to 2N_u on number-conserving states.  gamma is built explicitly, not
    assumed.
    """
    key = (L, U, periodic, t)
    if key in _HUB_CACHE:
        return _HUB_CACHE[key]
    n = L // 2
    T, S, idx = hop_matrix(L, n, periodic, t)
    D = len(S)
    occ = np.array([[(m >> i) & 1 for i in range(L)] for m in S], float)
    H = (sp.kron(T, sp.identity(D), format="csr")
         + sp.kron(sp.identity(D), T, format="csr")
         + sp.diags(U * (occ @ occ.T).ravel())).tocsr()
    v0 = np.random.default_rng(20260918).standard_normal(D * D)
    # machine-precision Lanczos where it is cheap; at L=10 a 1e-6 residual with a wider
    # Krylov basis gives the same F1 to six decimals in a tenth of the time.
    tol = 0.0 if D * D <= 6000 else 1e-6
    w, v = eigsh(H, k=1, which="SA", v0=v0, ncv=min(60, D * D - 1), tol=tol)
    psi = v[:, 0].reshape(D, D)
    rho = psi @ psi.conj().T                       # up-spin reduced density matrix
    G = onerdm_spin(rho, S, idx, L)
    g = np.clip(np.linalg.eigvalsh(G), 0.0, 1.0)
    # 4 tr[gamma(1-gamma)] over ALL spin-orbitals; at Sz=0 the down block repeats the up block
    F1 = float(8.0 * np.sum(g * (1.0 - g)))
    out = (F1, float(w[0]))
    _HUB_CACHE[key] = out
    return out


# ---------------------------------------------------------------------------
# Section 1b -- the geminal (APSG) witness, recomputed.
# Construction follows src/apsg_witness.py; the Jordan-Wigner operators are
# rebuilt here so that verify.py stays numpy+scipy only (apsg_witness.py imports
# n19_suite.py, which imports pyscf).
# ---------------------------------------------------------------------------

def jw_ops(M):
    dim = 1 << M
    C = []
    for p in range(M):
        r, c, d = [], [], []
        for s in range(dim):
            if (s >> p) & 1:
                sg = -1.0 if (bin(s & ((1 << p) - 1)).count("1") & 1) else 1.0
                r.append(s & ~(1 << p)); c.append(s); d.append(sg)
        C.append(sp.csr_matrix((d, (r, c)), shape=(dim, dim)))
    return C, [x.T.conj().tocsr() for x in C]


def apsg_state(thetas):
    r"""Antisymmetrised product of K strongly orthogonal 2e singlet geminals:
        |Psi_K> = prod_i ( cos th |b_i.up b_i.dn> + sin th |a_i.up a_i.dn> ) |0>
    Geminal i occupies spatial orbitals (2i, 2i+1) -> spin-orbitals 4i..4i+3.
    Returns (psi, C, Cd, M, Ne)."""
    K = len(thetas)
    M = 4 * K
    C, Cd = jw_ops(M)
    psi = np.zeros(1 << M, complex)
    psi[0] = 1.0
    for i, th in enumerate(thetas):
        bu, bd, au, ad = 4 * i, 4 * i + 1, 4 * i + 2, 4 * i + 3
        psi = np.cos(th) * (Cd[bu] @ (Cd[bd] @ psi)) + np.sin(th) * (Cd[au] @ (Cd[ad] @ psi))
    psi /= np.linalg.norm(psi)
    return psi, C, Cd, M, 2 * K


def majorana_fk(psi, C, Cd, M, kmax=3):
    """F_k = M - 0.5 * sum(sv^{2k}) from the Williamson spectrum of the Majorana
    covariance matrix -- the definition used throughout the paper."""
    G = []
    for p in range(M):
        cp = C[p] @ psi
        dp = Cd[p] @ psi
        G.append(cp + dp)
        G.append(-1j * (cp - dp))
    G = np.array(G)
    Mc = np.zeros((2 * M, 2 * M))
    for a in range(2 * M):
        for b in range(a + 1, 2 * M):
            val = (-1j * np.vdot(G[a], G[b])).real
            Mc[a, b] = val
            Mc[b, a] = -val
    sv = np.linalg.svd(Mc, compute_uv=False)
    return {k: float(M - 0.5 * np.sum(sv ** (2 * k))) for k in range(1, kmax + 1)}


def gamma_so(psi, C, M):
    """Spin-orbital 1-RDM gamma_pq = <psi| c^dag_p c_q |psi>."""
    cp = [C[p] @ psi for p in range(M)]
    g = np.zeros((M, M), complex)
    for p in range(M):
        for q in range(M):
            g[p, q] = np.vdot(cp[p], cp[q])
    return g


def n_unpaired(psi, C, Cd, ncas):
    """N_u = sum_a n_a (2 - n_a) over spatial orbitals (gauge-invariant
    multireference index)."""
    tot = 0.0
    for a in range(ncas):
        na = (np.vdot(psi, Cd[2 * a] @ (C[2 * a] @ psi)).real
              + np.vdot(psi, Cd[2 * a + 1] @ (C[2 * a + 1] @ psi)).real)
        tot += na * (2.0 - na)
    return float(tot)


def support_eps(psi, thr=1.0 - 1e-9):
    """eps-support: smallest number of determinants carrying `thr` of the weight.
    Same recipe as apsg_witness.support."""
    w = np.sort(np.abs(psi) ** 2)[::-1]
    c = np.cumsum(w) / w.sum()
    return int(np.searchsorted(c, thr) + 1)


def schmidt_rank(psi, M, cut, tol=1e-10, ptol=1e-12):
    """Schmidt rank across the mode bipartition at position `cut`.  Identical
    reshape convention to apsg_witness.chi_cut."""
    A = psi.reshape(1 << cut, 1 << (M - cut))
    s = np.linalg.svd(A, compute_uv=False)
    s = s[s > tol]
    p = s ** 2 / np.sum(s ** 2)
    return int(np.sum(p > ptol))


def witness_row(K, theta):
    psi, C, Cd, M, Ne = apsg_state([theta] * K)
    F = majorana_fk(psi, C, Cd, M, kmax=3)
    g = gamma_so(psi, C, M)
    occ = np.clip(np.linalg.eigvalsh(g).real, 0.0, 1.0)
    f1_gamma = float(4.0 * np.sum(occ * (1.0 - occ)))
    nu = n_unpaired(psi, C, Cd, 2 * K)
    ranks = [schmidt_rank(psi, M, cut) for cut in range(1, M)]
    return dict(K=K, theta=theta, M=M, Ne=Ne,
                F1=F[1], F2=F[2], F3=F[3],
                F1_gamma=f1_gamma, twoNu=2.0 * nu,
                S_raw=int(np.sum(np.abs(psi) ** 2 > 1e-14)),
                S_eps=support_eps(psi),
                chi_central=schmidt_rank(psi, M, 2 * K),
                chi_max=max(ranks))


# ---------------------------------------------------------------------------
# Small helpers for reading the deposited artefacts
# ---------------------------------------------------------------------------

COORD = re.compile(r"\(\s*(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*,\s*(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*\)")
NUMTOK = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


# ---------------------------------------------------------------------------
# THE v3 MANUSCRIPT BODY
#
# v2 was eight body files; v3 is ten sections plus seven appendices, and the prose
# this guardian pins moved rather than vanished -- often split, so that one v2 file's
# sentences now live in two or three v3 files.  read_body() therefore hands the
# sentence-level checks the WHOLE body, expanded from paper/main.tex in driver order.
#
# This is not a relaxation.  Every check still has to find its sentence, and the number
# in that sentence still has to match the recomputation; only the haystack grew from
# one file to the body.  Deriving the list from main.tex also means a section added
# later is searched automatically, and one removed stops being searched -- neither can
# happen without this list changing.
# ---------------------------------------------------------------------------
_BODY_CACHE = {}


def body_files(root):
    """Every .tex paper/main.tex reaches, in the order the driver inputs them."""
    out, seen = [], set()

    def walk(rel):
        rel = rel.replace("\\", "/")
        if rel in seen:
            return
        seen.add(rel)
        path = os.path.join(root, "paper", rel.replace("/", os.sep))
        if not os.path.exists(path):
            return
        out.append("paper/" + rel)
        with open(path, encoding="utf-8", errors="replace") as fh:
            txt = fh.read()
        txt = "\n".join(l for l in txt.split("\n") if not l.lstrip().startswith("%"))
        for m in re.finditer(re.escape(chr(92)) + r"input\{([^}]+)\}", txt):
            s = m.group(1)
            if not s.endswith(".tex"):
                s += ".tex"
            for cand in (s, os.path.dirname(rel) + "/" + s if os.path.dirname(rel) else s):
                if os.path.exists(os.path.join(root, "paper", cand.replace("/", os.sep))):
                    walk(cand)
                    break

    walk("main.tex")
    return out


def read_body(root):
    """The concatenated manuscript body, comments stripped."""
    key = os.path.abspath(root)
    if key not in _BODY_CACHE:
        chunks = []
        for rel in body_files(root):
            with open(os.path.join(root, rel.replace("/", os.sep)),
                      encoding="utf-8", errors="replace") as fh:
                txt = fh.read()
            chunks.append("\n".join(l for l in txt.split("\n")
                                    if not l.lstrip().startswith("%")))
        _BODY_CACHE[key] = "\n".join(chunks)
    return _BODY_CACHE[key]


def lines_body(root):
    return read_body(root).split("\n")


def read(root, rel):
    path = os.path.join(root, rel.replace("/", os.sep))
    with open(path, encoding="utf-8", errors="replace") as fh:
        return fh.read()


def read_json(root, rel):
    return json.loads(read(root, rel))


def lines(root, rel):
    return read(root, rel).split("\n")


def chi_of(row, name="<row>"):
    """The central-cut Schmidt rank of a witness row, under either field name.

    C11's own report recommended renaming this field from `chi` to `chi_central_cut`
    (because it is NOT the minimal MPS bond dimension of the theorem, which is the maximum
    over all cuts).  Following that recommendation used to kill this guardian with a
    KeyError and exit status 2 -- a guardian that dies when the deposit improves is worse
    than no guardian.  Both names are accepted; neither value is excused."""
    for k in ("chi", "chi_central_cut"):
        if k in row:
            return row[k], k
    raise KeyError("%s carries neither 'chi' nor 'chi_central_cut'; the central-cut "
                   "Schmidt rank of the geminal witness has no field at all" % name)


def coords_of(line):
    """(x,y) pairs on a pgfplots \addplot line, excluding \addlegendentry text."""
    body = line.split("\\addlegendentry")[0]
    return [(float(a), float(b)) for a, b in COORD.findall(body)]


def find_line(ls, needle, start=0):
    for i in range(start, len(ls)):
        if needle in ls[i]:
            return i
    return -1


def comb(n, k):
    return math.comb(n, k)


# Declared slack on a RECOMPUTED quantity before it is compared with a deposited
# literal.  A bare rounding test is a knife edge: F1(L=8,OBC) = 12.8162762237 sits
# 2.8e-7 from the six-decimal rounding boundary, so a different BLAS/ARPACK build could
# make this guardian accuse an impeccable repository of carrying the wrong system --
# in a public, permanent refereeing record.  The slack is 2000x smaller than the
# smallest OBC/PBC separation it has to resolve (2.455e-3 at L=8), so it costs no
# discrimination whatsoever.
SOLVER_SLACK = 5e-7


def rounds_to(value, printed, ndp, slack=SOLVER_SLACK):
    """True iff the deposited literal `printed` is the correct rounding of `value` to
    ndp decimals, allowing SOLVER_SLACK of eigensolver jitter on `value`."""
    return abs(float(value) - float(printed)) <= 0.5 * 10 ** (-ndp) + slack


def infer_ndp(printed):
    """How many decimals the deposited literal actually carries.  The tree stores the
    same physical number at different precisions (9.604689 in one file, 12.8163 in
    another), so the comparison has to be made at the precision that is written down."""
    s = repr(float(printed))
    if "e" in s or "E" in s:
        return 12
    return len(s.split(".")[1]) if "." in s else 0


def panel_blocks(ls):
    """Split a pgfplots `groupplot` source into its panels.

    Returns a list of (title, [lines]) in plotting order.  Panel-aware parsing is not
    cosmetic: in fig_resource_master_native.tex the same `rmChn` series appears twice,
    once against chi and once against F1, and a parser that cannot tell them apart
    compares a bond dimension with a magic."""
    out, cur = [], None
    for ln in ls:
        if "\\nextgroupplot" in ln:
            if cur is not None:
                out.append(cur)
            cur = [ln]
        elif cur is not None:
            cur.append(ln)
    if cur is not None:
        out.append(cur)
    named = []
    for i, blk in enumerate(out):
        # the title may sit several lines below \nextgroupplot, inside its option list
        m = re.search(r"title=\{\(([a-z])\)", "\n".join(blk[:8]))
        named.append((m.group(1) if m else str(i), blk))
    return named


# ---------------------------------------------------------------------------
# The two-sided OBC/PBC check.  This is the check whose absence let the
# 2026-08-17 defect through.
# ---------------------------------------------------------------------------

def system_check(rep, name, printed, L, system, obc, pbc, ndp=None, where=""):
    """`printed` must carry the value of `system` at size L, and must NOT carry the
    value of the other system.  ndp=None compares at full precision; otherwise the
    comparison is against the correctly rounded value."""
    mine, other = (obc[L], pbc[L]) if system == "OBC" else (pbc[L], obc[L])
    myname, othername = (("OPEN chain (OBC)", "PERIODIC ring (PBC)") if system == "OBC"
                         else ("PERIODIC ring (PBC)", "OPEN chain (OBC)"))
    if ndp is None:
        ndp = min(infer_ndp(printed), 6)
    ok_mine = rounds_to(mine, printed, ndp)
    ok_other = rounds_to(other, printed, ndp)
    if ok_mine and not ok_other:
        rep.emit(PASS, name, "%s L=%d: %.6g is the %s value" % (where, L, float(printed), myname))
        rep.numbers += 0  # counted by emit
    elif ok_other and not ok_mine:
        rep.emit(FAIL, name,
                 "WRONG SYSTEM: %s carries %.6g, which is the %s value; this slot is the %s "
                 "and must carry %.6f. They are different systems, not two versions of one "
                 "number (see the provenance note in data/ladder_vs_chain.json)."
                 % (where, float(printed), othername, myname, mine))
    elif ok_mine and ok_other:
        rep.emit(FAIL, name,
                 "PRECISION TOO LOW TO DISCRIMINATE: %s L=%d carries %s, which rounds to both "
                 "the %s value %.6f and the %s value %.6f at ndp=%s. A literal written at a "
                 "precision that cannot tell the two systems apart disables the only check that "
                 "catches the 2026-08-17 defect; write it with enough digits."
                 % (where, L, printed, myname, mine, othername, other, ndp))
    else:
        rep.emit(FAIL, name,
                 "UNKNOWN VALUE: %s carries %.6g, which is neither the %s value %.6f nor the "
                 "%s value %.6f" % (where, float(printed), myname, mine, othername, other))


# ---------------------------------------------------------------------------
# Registry of RETRACTED CLAIMS -- sentences, not scalars.
#
# 2026-09-19: all seven were still printed by the v2 manuscript and none is printed by
# the v3 one, so the registry was inverted rather than deleted.  Each entry now asserts
# that its sentence is ABSENT; its return is a hard failure.  The text of each entry is
# left exactly as it was, because the `why` is the record of why the sentence went.
#
# KNOWN_OPEN pins a wrong NUMBER in one slot.  That is not enough: the sum-rule
# defect, for instance, lives in three places (a JSON field, a figure label and a
# caption), and repairing the JSON alone made the previous guardian print
# "REPAIRED -- delete the registry entry" while the false sentence was still printed
# twice in the PDF.  This registry pins the SENTENCES.  Each entry is an exact string
# that is currently in the tree and is known to be wrong or unsupported; it is reported
# as XFAIL while it is there and as XPASS (with an instruction to delete the entry) once
# it goes.  Nothing here sets the exit code -- but nothing here can be edited silently,
# and the list is the retraction agenda for v3, carried by the guardian rather than by a
# note somebody has to remember to read.
# ---------------------------------------------------------------------------
CLAIMS_RETRACTED = {
    "sumrule.figure_label": dict(
        where="paper/figs/fig_method_native_frag.tex",
        text=r"$\int\! A=0.500$",
        why="the figure LABELS 0.500 as an integral of A.  It is not one: it is the "
            "analytic identity 1-<n_p> stored in data/method_max.json as 'sumrule'.  The "
            "integral of the deposited A_exact over the plotted window is 0.48659 and of "
            "the reconstruction 0.48647 (2.7 per cent below, a finite-window effect).  "
            "Repairing data/method_max.json alone does NOT repair this."),
    "sumrule.results_caption": dict(
        where="body",
        text=r"to machine precision",
        why="the caption of the method figure says the reconstructed weights obey "
            "\\int A^+ dw = 0.500 'to machine precision'.  The deviation is 2.7 per cent "
            "(finite window).  The companion figure's caption already admits about 2 per "
            "cent: the two captions contradict each other."),
    "abstract.fci_1e-5": dict(
        where="paper/main.tex",
        text=r"FCI-verified to $<10^{-5}$~Ha",
        why="the molecular energies are CASCI in the same active space compared with "
            "themselves; the residual confessed in paper/theorem_b2.tex is 1.1e-13 Ha.  "
            "The phrase states an external verification that was not performed."),
    "abstract.runs_on_heron": dict(
        where="paper/main.tex",
        text=r"and runs on the",
        why="data/heron_spectral.json carries hw_S = 300 of nsector = 300: the retained "
            "subspace is the WHOLE sector, so A_hw == A_exact by coverage, not by "
            "fidelity.  The caption of the hardware figure already says 'exact by "
            "coverage'; the abstract and the conclusion do not.  'executed on' is the "
            "defensible verb."),
    "abstract.not_significant": dict(
        where="paper/main.tex",
        text=r"($n=30$, not significant)",
        why="pooled over the six (L,N) strata of data/cost_vs_ent.json the F1-|S| rank "
            "correlation reads -0.11; WITHIN each of the six strata it is exactly -1.000. "
            "This is Simpson's paradox: the stratified evidence is far stronger than the "
            "published null, in the opposite direction.  'Not significant' is wrong."),
    "results.qto0_charge_gap": dict(
        where="body",
        text=r"as $q\to0$ its lower edge approaches the charge gap",
        why="the smallest accessible momentum is q/pi = 1/6 and its onset is 5.600, which "
            "is 12.7 per cent ABOVE Delta = 4.968759 at the 3 per cent threshold the "
            "generator uses.  The approach was illustrated by two rows of "
            "paper/figs/sqw_edges.dat that no computation produced (q=0 and q=2, with the "
            "lower edge hand-set equal to Delta).  Those rows are gone; the sentence is "
            "not."),
    "hw.configuration_recovery": dict(
        where="body",
        text=r"configuration recovery returned",
        why="the L=6 ibm_fez run (job d9s16avpemts73ct6g8g) was plain post-selection on particle "
            "number and S_z, in reversed bit order (notebooks/Spectral_Heron.ipynb cell 15); no "
            "configuration recovery was applied, and every retained shot is a device error "
            "(data/hw_bitorder_check.json)."),
    "hw.trex": dict(
        where="body",
        text=r"TREX",
        why="SamplerV2 ran with measurement twirling (twirling.enable_measure); it returns twirled raw "
            "counts and applies no readout-error extinction, so the run had no readout mitigation."),
    "hw.proof_of_principle": dict(
        where="paper/carried/fig_heron.tex",
        text=r"a proof of principle",
        why="the retained shots of panel (a) are all device errors (reversed bit order); the panel "
            "is an execution record, not a proof of principle of the method."),
    # 2026-09-28: the same claim stood in Sec. S8, in the paragraph that itself discloses the
    # reversed bit order: "recovery to the correct particle-number sector ... which is what a
    # proof of principle is".  The sentence now claims executability only.  Both phrases are
    # refused in paper/sm_scope.tex; "recovery to the correct" is not refused body-wide,
    # because configuration recovery of the ibm_marrakesh run does restore particle number.
    "hw.sm_scope_proof_of_principle": dict(
        where="paper/sm_scope.tex",
        text=r"proof of principle",
        why="Sec. S8 discloses that the L=6 ibm_fez run was post-selected in reversed bit order, so "
            "every retained determinant is a device error (data/hw_bitorder_check.json); the run "
            "shows that the steps execute end to end, not a proof of principle of the method."),
    "hw.sm_scope_recovery_to_correct": dict(
        where="paper/sm_scope.tex",
        text=r"recovery to the correct",
        why="the L=6 ibm_fez run applied no configuration recovery, and its post-selection read the "
            "bitstrings in reversed bit order (notebooks/Spectral_Heron.ipynb cell 15; "
            "data/hw_bitorder_check.json: 0 of 350000 noiseless shots kept under that reading)."),
    "fig_scaling2.beyond_classical": dict(
        where="paper/figs/fig_scaling2_native.tex",
        text=r"{beyond-classical}",
        why="panel (b) shades |S|/dim in [0,0.5] and labels it 'beyond-classical', with no "
            "criterion, no mention in the caption and no mention in the text -- in a paper "
            "that disclaims any beyond-classical claim six times.  The figure says what "
            "the text does not."),
}


def claim_check(rep, root, key):
    """Check one entry of CLAIMS_RETRACTED.  The assertion is ABSENCE.

    These seven sentences were the v2 retraction agenda.  v3 removed all seven, so the
    useful check is no longer "is it still there" but "has it come back".  A registry
    entry that is merely deleted once the sentence goes protects nothing; this one keeps
    protecting.  `where` may be a file, or the sentinel "body", which means the whole
    manuscript with comment lines stripped -- so a retracted sentence cannot be
    reinstated by moving it into another section.
    """
    e = CLAIMS_RETRACTED[key]
    if e["where"] == "body":
        txt = read_body(root)
    else:
        try:
            txt = read(root, e["where"])
        except IOError:
            rep.missing("claim %s" % key, "cannot open %s" % e["where"])
            return
    name = "%s (%s)" % (key, e["where"])
    if e["text"] in txt:
        rep.emit(FAIL, name,
                 "RETRACTED CLAIM IS BACK: %r is printed again in %s. It was withdrawn "
                 "for this reason and the reason has not changed. %s"
                 % (e["text"], e["where"], e["why"]), 1)
    else:
        rep.emit(PASS, name,
                 "retracted and still absent: %r is not in %s" % (e["text"], e["where"]), 1)


# ---------------------------------------------------------------------------
# Section 9 helpers
# ---------------------------------------------------------------------------

BS = chr(92)              # a backslash, written this way so no escape can be mis-read
RBS = re.escape(BS)       # the same backslash, safe to paste into a regex
RE_INPUT = re.compile(re.escape(BS + "input{figs/") + r"([A-Za-z0-9_]+)\.tex\}")

# Figure sources that THIS guardian opens and anchors to a dataset.  The structural
# check below refuses to pass if paper/*.tex \input's a fragment that is not in here:
# coverage can then never silently shrink by adding an unwatched figure.
COVERED_FRAGMENTS = {
    "fig_decoupling_native.tex",
    "fig_resource_master_native.tex",
    "fig_scaling2_native.tex",
    "fig_witness_native.tex",
    "fig_ladder_native.tex",
    "fig_akw_sampled_honest.tex",
    "fig_method_native_frag.tex",
    "fig_hardware_hero_frag.tex",
    "fig_noise_recovery_native.tex",
    # v3: fig:gapscaling and fig:thm1iii-violation are new floats.  The first is checked
    # below against data/gap_scaling.json.  The second is NOT, and that is recorded in
    # COVERAGE_GAP_V3 rather than papered over.
    "fig_gapscaling_native.tex",
    "fig_thm1iii_violation_native.tex",
    # fig_amort_native.tex is no longer \input by anything: Sec. 9.4 withdraws that
    # figure in the body ("Both figures of this test are therefore withdrawn").  It stays
    # in paper/figs/ because src/ still names it; nothing typesets it.
}

# ---------------------------------------------------------------------------
# COVERAGE GAPS this guardian knows about and does NOT hide.  Unlike KNOWN_OPEN,
# which pins a wrong number, this pins a MISSING CHECK: an artefact the manuscript
# typesets and no check in this file anchors to a dataset.
#
# It is emitted as XFAIL, exactly the convention KNOWN_OPEN uses -- printed on every
# run, carried into the summary, and NOT setting the exit code.  (An earlier version of
# this comment said it set the exit code.  It did not; nothing read the dictionary at
# all.  Fixed 2026-09-19, and said out loud because a false comment in the guardian is
# the same defect as a false sentence in the paper.)
#
# Close an entry by WRITING THE CHECK, not by deleting the entry.
# ---------------------------------------------------------------------------
COVERAGE_GAP_V3 = {
    # "tab:moments (Sec. IV)" was declared here until 2026-09-25; it is closed by
    # src/moments_table.py and the assertions in section 9.2 below.
    # 2026-09-26: build_n3.py and its inputs were recovered and deposited, and
    # src/check_figures.py now regenerates the ten n3_*.dat tables byte-for-byte.  The entry
    # stayed, because what it named -- a check IN THIS FILE -- had not been written.
    # 2026-09-28 (late): "fig_thm1iii_violation_native.tex" is closed by section 9.6c below.
    # build_n3.py now reads the DEPOSITED data/cert_{stress,akw,teqsci}.json (until that day
    # it read an earlier run of the same scans, kept in data/thm1iii_violation/cert_earlier_run/),
    # and 9.6c re-classifies the pool from the five JSON sources with code that shares nothing
    # with build_n3.py, then checks the eight plotted tables, the whisker summary and every
    # count the caption, Sec. III D, Sec. S2 and Sec. S9 print against that classification.
}

# Minimum number of numeric assertions each section must actually perform.  A guardian
# that only checks what it happens to find cannot notice a DELETED artefact: erasing a
# whole \addplot series used to take the total from 996 to 984 and still print PASS.
# These floors are the fix.  They are the counts of the clean tree; raise them when you
# add checks, never lower them to make a run go green.
# measured on the clean tree of 2026-09-18, less a declared 2 % margin.  The margin exists
# so that a legitimately absent optional file (a CI workflow, say) cannot make the guardian
# cry wolf; it is far too small to hide a deleted artefact, which costs tens to hundreds.
SECTION_FLOORS = [
    ("(2) RECOMPUTED", 32),
    ("(3) WITNESS ARTEFACTS", 98),      # 2026-09-28 (late): 100 measured, chi_max checks added
    ("(4) OPEN CHAIN", 62),
    ("(5) DERIVED QUANTITIES", 440),
    ("(6) STATISTICS", 23),
    ("(7) MOLECULAR SUITE", 263),
    ("(8) SPECTRAL ARTEFACTS", 69),
    ("(9) THE FILES", 10850),           # 2026-09-30 (final review): 11071 measured, less 2 %,
                                        # after sections 9.6b (caption panel (b)), 9.9b and 9.9c;
                                        # 2026-09-28 (late): 11003 measured, less 2 %; sections
                                        # 9.6b/9.6c (Fig. 3 finite shots, the free-fermion
                                        # supports, the spin peaks, the Fig. S5 pool) added
    ("(10b) C3 FRONTIER", 2387),
    ("(10) RETRACTED CLAIMS", 7),   # the header was renamed when the registry was inverted;
                                   # the floor of 7 is unchanged and is met by the seven
                                   # absence assertions, not by lowering it.
    ("(11) COVERAGE FLOOR", 0),
]


def all_coords(text):
    """Every (x,y) pair in every `coordinates {...}` block of a pgfplots fragment."""
    out = []
    for blk in re.findall(r"coordinates\s*\{(.*?)\}\s*(?:;|\\|$)", text, re.S):
        out.extend([(float(a), float(b)) for a, b in COORD.findall(blk)])
    return out


def logical_lines(text):
    """pgfplots statements, re-joined.  A statement may span several source lines
    (fig_hardware_hero_frag.tex breaks one \addplot across two), and a parser that
    reads raw lines silently sees an empty series -- which is exactly how a deleted
    artefact hides."""
    _b = chr(92)
    starts = tuple([_b + w for w in ("addplot", "node", "nextgroupplot", "fill",
                                     "draw", "begin", "end", "def", "definecolor",
                                     "addlegendentry")] + ["%"])
    out, cur = [], ""
    for ln in text.split(chr(10)):
        st = ln.lstrip()
        if cur and (st.startswith(starts) or cur.rstrip().endswith(";")):
            out.append(cur)
            cur = ln
        else:
            cur = (cur + " " + ln) if cur else ln
    if cur:
        out.append(cur)
    return out


def series_coords(line):
    """(x,y) pairs on one \\addplot line, ignoring the `+- (0,sigma)` error bars."""
    body = line.split(BS + "addlegendentry")[0]
    body = re.sub(r"\+-\s*\([^)]*\)", " ", body)
    return [(float(a), float(b)) for a, b in COORD.findall(body)]


def errbars(line):
    """the sigma values of `(x,y) +- (0,sigma)` pairs, in order."""
    return [float(m.group(1)) for m in
            re.finditer(r"\+-\s*\(\s*0\s*,\s*(-?[\d.eE+]+)\s*\)", line)]


def table_cols(root, rel):
    """Read a whitespace-separated .dat with a header line.  Returns (names, columns)."""
    ls = [l for l in read(root, rel).split(chr(10)) if l.strip()]
    names = ls[0].split()
    cols = [[] for _ in names]
    for l in ls[1:]:
        parts = l.split()
        for i, p in enumerate(parts):
            cols[i].append(float(p))
    return names, cols


def curve_on_grid(rep, name, pairs, grid, ref, ndp, defect, offset=None):
    """Compare a plotted (x,y) series against a reference array sampled on `grid`.

    The figure plots a SUBSET of the deposited grid (every third point), so matching by
    position is wrong; each plotted abscissa is looked up in the deposited grid instead.
    An abscissa that is not on the grid is a fabricated point and says so."""
    idx = {}
    for i, g in enumerate(grid):
        idx.setdefault(round(float(g), ndp), i)
    off = offset
    bad, orphan = [], []
    for x, y in pairs:
        i = idx.get(round(float(x), ndp))
        if i is None:
            orphan.append(x)
            continue
        if off is None:
            off = round(float(y) - float(ref[i]), ndp)
        if abs(float(y) - (float(ref[i]) + off)) > 0.5 * 10 ** (-ndp) + 1e-9:
            bad.append((x, y, ref[i] + off))
    if orphan:
        rep.emit(FAIL, name, "FABRICATED ABSCISSAE: %d plotted point(s) sit at frequencies "
                             "that are not on the deposited grid, first %.4f"
                 % (len(orphan), orphan[0]), n=len(pairs))
        return False, off
    if bad:
        rep.emit(FAIL, name, "%s: %d of %d plotted points do not reproduce (offset %s); "
                             "first is x=%.4f, plotted %.4f, deposited %.4f"
                 % (defect, len(bad), len(pairs), off, bad[0][0], bad[0][1], bad[0][2]),
                 n=len(pairs))
        return False, off
    rep.emit(PASS, name, "%d plotted points reproduce from the dataset (offset %s)"
             % (len(pairs), off), n=len(pairs))
    return True, off


def curve_check(rep, name, plotted, ref, ndp, defect, offset=None):
    """Compare a plotted series with a reference array, allowing a constant offset.

    Returns (ok, offset).  `plotted` is a list of y-values already rounded to ndp by
    the generator; `ref` is the full-precision array.  The offset (curves 'offset
    vertically by k') is INFERRED and then required to be constant: a forged curve that
    is a shifted copy of the wrong array fails on the second point."""
    if len(plotted) != len(ref):
        rep.emit(FAIL, name, "%s: %d plotted points against %d deposited values"
                 % (defect, len(plotted), len(ref)), n=1)
        return False, None
    if offset is None:
        offset = round(float(plotted[0]) - float(ref[0]), ndp)
    tol = 0.5 * 10 ** (-ndp) + 1e-9
    bad = [(i, plotted[i], ref[i] + offset)
           for i in range(len(ref)) if abs(plotted[i] - (ref[i] + offset)) > tol]
    if bad:
        rep.emit(FAIL, name, "%s: %d of %d plotted points do not reproduce (offset %.4f); "
                             "first is index %d, plotted %.4f, deposited %.4f"
                 % (defect, len(bad), len(ref), offset, bad[0][0], bad[0][1], bad[0][2]),
                 n=len(ref))
        return False, offset
    rep.emit(PASS, name, "%d points reproduce from the dataset (offset %.4f)"
             % (len(ref), offset), n=len(ref))
    return True, offset


# ===========================================================================
#  SECTION 9 -- the files no earlier guardian opened
# ===========================================================================

def section9(rep, root, ctx):
    np_ = np
    rep.head("(9) THE FILES NO EARLIER GUARDIAN OPENED -- abstract, prose, figure sources")

    obc, pbc = ctx["obc"], ctx["pbc"]

    # ---------------------------------------------------------------- 9.1
    # The 2026-08-17 defect, hunted across EVERY source file rather than in the five
    # slots somebody happened to list.  A PBC literal is legal only inside the scaling
    # figure (whose engine is data/scaling_data.json, the periodic ring); an OBC literal
    # is illegal inside that figure's plotted lines.  This is what catches the defect
    # when it is reintroduced in prose (results.tex) rather than in a dataset.
    # The periodic-ring literals (6.6074, 9.3851, ...) may appear in the ring figure
    # AND in that figure's own caption, which v3 moved into its own file.  Both are
    # required to SAY they are the ring, so a file cannot join this set by accident:
    # the scan below refuses any home that does not carry the word "periodic".
    PBC_HOMES = ("paper/figs/fig_scaling2_native.tex", "paper/figs/caption_fig12.tex")
    for _h in PBC_HOMES:
        rep.truth("PBC home %s declares periodic boundary conditions" % _h,
                  "periodic" in read(root, _h).lower(),
                  "the file says so in words, so a ring literal in it is labelled",
                  "%s is allowed to carry periodic-ring literals but never says it is "
                  "the ring; a reader cannot tell which system the numbers belong to" % _h)
    pbc_lits, obc_lits = {}, {}
    for L in sorted(obc):
        if L not in pbc:
            continue
        for nd in (2, 3, 4, 5, 6):
            a = "%.*f" % (nd, obc[L])
            b = "%.*f" % (nd, pbc[L])
            if a != b:
                obc_lits.setdefault(a, (L, nd))
                pbc_lits.setdefault(b, (L, nd))
    srcs = []
    for sub in ("paper", "paper/figs"):
        d = os.path.join(root, sub.replace("/", os.sep))
        for f in sorted(os.listdir(d)):
            if f.endswith(".tex"):
                srcs.append(sub + "/" + f)
    misplaced = []
    nscan = 0
    for rel in srcs:
        txt = read(root, rel)
        for lit, (L, nd) in pbc_lits.items():
            for m in re.finditer(r"(?<![\d.])" + re.escape(lit) + r"(?![\d])", txt):
                nscan += 1
                if rel not in PBC_HOMES:
                    ln = txt[:m.start()].count(chr(10)) + 1
                    misplaced.append((rel, ln, lit, "PERIODIC ring (PBC)", L,
                                      "%.6f" % obc[L]))
        for lit, (L, nd) in obc_lits.items():
            for m in re.finditer(r"(?<![\d.])" + re.escape(lit) + r"(?![\d])", txt):
                nscan += 1
                if rel in PBC_HOMES:
                    ln = txt[:m.start()].count(chr(10)) + 1
                    # comments are the place where the note "9.6047 is the OTHER system"
                    # legitimately lives
                    line = txt.split(chr(10))[ln - 1]
                    if line.lstrip().startswith("%"):
                        continue
                    misplaced.append((rel, ln, lit, "OPEN chain (OBC)", L,
                                      "%.6f" % pbc[L]))
    rep.truth("no OBC/PBC literal sits in the wrong file (whole-tree scan of %d .tex sources)"
              % len(srcs),
              not misplaced,
              "%d system-discriminating literals scanned, all in the right file" % nscan,
              "WRONG SYSTEM: %s -- these slots must carry the other system's value. This is "
              "the 2026-08-17 defect, and it is reintroducible in prose, not only in the "
              "five slots previously listed." % misplaced,
              n=max(nscan, 1))

    # ---------------------------------------------------------------- 9.2
    # v3 inserts a float WRAPPER between the body and the fragment
    # (sec_X.tex -> figs/float_nN.tex -> figs/fig_..._native.tex), so the scan has to
    # follow \input through paper/figs/ as well.  Reading only the body files found
    # five wrappers and none of the fragments behind them, i.e. it would have reported
    # the wrappers as unwatched while the real figure sources went unnoticed.
    inputs, frontier, seen_in = set(), [r for r in srcs if not r.startswith("paper/figs/")], set()
    while frontier:
        rel = frontier.pop()
        if rel in seen_in:
            continue
        seen_in.add(rel)
        for m in RE_INPUT.finditer(read(root, rel)):
            name = m.group(1) + ".tex"
            sub_rel = "paper/figs/" + name
            if os.path.exists(os.path.join(root, sub_rel.replace("/", os.sep))) \
               and ("float_" in name or name.endswith("caption.tex")
                    or name.startswith("caption") or name == "captions.tex"
                    or name.endswith("_caption.tex")):
                frontier.append(sub_rel)       # a wrapper / caption file: keep walking
            else:
                inputs.add(name)
    for _gapkey in sorted(COVERAGE_GAP_V3):
        # The key used to be pasted behind a literal "paper/figs/", which was right while
        # every gap was a figure fragment and became a FALSE PATH the moment one was not:
        # the tab:moments entry printed "coverage gap: paper/figs/tab:moments (Sec. IV)",
        # naming a file that does not exist.  A guardian that prints a path it has not
        # checked is committing the defect it exists to catch.  Keys that name a real
        # fragment keep their directory; keys that name something else print as they are.
        _donde = ("paper/figs/%s" % _gapkey) if _gapkey.endswith(".tex") else _gapkey
        rep.emit(XFAIL, "coverage gap: %s" % _donde,
                 "the manuscript typesets this and no check in this file anchors "
                 "it to a dataset. %s" % COVERAGE_GAP_V3[_gapkey], 1)
    # ------------------------------------------------ tab:moments, regenerated (2026-09-25)
    # The four lost artefacts are replaced by src/moments_table.py, which rebuilds the
    # L = 6 open-chain sector from scratch (independently of src/frontier/c0_lib.py) and
    # writes data/moments_table.json.  The rows printed in sec_4_body.tex, the negative
    # control and the two-level series of sec_4_app.tex are asserted against it.  The
    # roundoff entries (max_{j<=2K+1} eps_j, "largest over six implementations") are not
    # asserted digit by digit: they are floating-point residues, and only their order is
    # claimed.
    def _tex_sci(x, nd):
        e = int(math.floor(math.log10(abs(x))))
        m = round(x / 10 ** e, nd)
        if m >= 10:
            m, e = m / 10, e + 1
        return ("%.*f" % (nd, m)) + "\\times10^{%d}" % e
    _mt = os.path.join(root, "data", "moments_table.json")
    if not os.path.exists(_mt):
        rep.missing("tab:moments regenerated",
                    "data/moments_table.json is absent; run python src/moments_table.py")
    else:
        with open(_mt, encoding="utf-8") as _fh:
            _M = json.load(_fh)
        _rows = {r["K"]: r for r in _M["table"]["rows"]}
        with open(os.path.join(root, "paper", "sec_4_body.tex"), encoding="utf-8") as _fh:
            _t4 = _fh.read()
        with open(os.path.join(root, "paper", "sec_4_app.tex"), encoding="utf-8") as _fh:
            _t4a = " ".join(_fh.read().split())
        for _K in (0, 1, 2):
            _m = re.search(r"\$%d\$ & \$(\d+)\$ & \$([\d.]+)\$" % _K, _t4)
            rep.eq_int("tab:moments |S| at K=%d" % _K, int(_m.group(1)) if _m else -1,
                       _rows[_K]["S"])
            rep.eq("tab:moments |S|/D at K=%d" % _K, float(_m.group(2)) if _m else -1.0,
                   _rows[_K]["frac"], tol=5e-4)
        for _K in (0, 1):
            _rel = _tex_sci(_rows[_K]["rel_eps_2K2"], 1)
            _abs = _tex_sci(_rows[_K]["abs_eps_2K2"], 2)
            rep.truth("tab:moments eps_{2K+2} at K=%d" % _K, _rel in _t4 and _abs in _t4,
                      "relative %s and absolute %s as printed" % (_rel, _abs),
                      "the table does not print the regenerated %s (relative) and %s "
                      "(absolute)" % (_rel, _abs))
        _nc = _M["negative_control"]
        rep.truth("tab:moments negative control", _nc["fail_at_zeroth_moment"] == _nc["draws"]
                  == 300 and "every single draw fails" in " ".join(_t4.split()),
                  "300 of 300 random size-200 subspaces fail at the zeroth moment",
                  "regenerated control: %d of %d fail" % (_nc["fail_at_zeroth_moment"],
                                                         _nc["draws"]))
        _tl = ",\\ ".join("%.4f" % x["L1"] for x in _M["two_level"])
        rep.truth("app:moments:conj two-level series", _tl in _t4a,
                  "the seven values printed are the converged ones",
                  "sec_4_app.tex does not print the converged series %s" % _tl)
    # ------------------------------------------------ the support witness S* = supp(phi)
    # (2026-09-25) supp(phi) is counted above 1e-12 of the largest amplitude, across a gap
    # that src/support_witness.py checks (W3) and against the independent k = 0 count of
    # z_suppK_sym (W8).  Every sentence that quotes the witness, its size, the support floor
    # of Sec. III or the out-of-support share of the Born cut is asserted here.
    _sw = read_json(root, "data/support_witness.json")
    _swr = {r["L"]: r for r in _sw["results"]}
    _rel = lambda L, e: [t["relL1"] for t in _swr[L]["table"] if abs(t["eta"] - e) < 1e-12][0]
    _body = " ".join(" ".join(open(os.path.join(root, "paper", f), encoding="utf-8").read()
                              .split()) for f in sorted(os.listdir(os.path.join(root, "paper")))
                     if f.endswith(".tex"))
    _q = "$0.43$ at $L=6$ and $%.2f$ at $L=8$" % _rel(8, 0.18)
    _all18 = re.findall(r"\$0\.\d\d\$ at \$L=6\$ and \$0\.\d\d\$ at \$L=8\$", _body)
    rep.truth("support witness: rel-L1(0.18) at L=6, 8", "%.2f" % _rel(6, 0.18) == "0.43"
              and _body.count(_q) >= 1 and all(x == _q for x in _all18),
              "%s, %d times" % (_q, _body.count(_q)),
              "the witness errors printed are not those of data/support_witness.json (%.4f, %.4f)"
              % (_rel(6, 0.18), _rel(8, 0.18)))
    _q = "rising to $%.2f$ and $%.2f$ at $\\eta=0.05" % (_rel(6, 0.05), _rel(8, 0.05))
    _all05 = re.findall(r"rising to \$[\d.]+\$ and \$[\d.]+\$ at \$\\eta=0\.05", _body)
    rep.truth("support witness: rel-L1(0.05) at L=6, 8", _body.count(_q) >= 1
              and all(x == _q for x in _all05),
              "%s, %d times" % (_q, _body.count(_q)),
              "the eta = 0.05 t witness errors printed are not %.4f, %.4f"
              % (_rel(6, 0.05), _rel(8, 0.05)))
    rep.eq("support witness: |S*|/D at L=8 == 0.571", _swr[8]["frac"], 0.571, tol=5e-4)
    rep.eq_int("support witness: |S*| at L=8 printed as 2237", _swr[8]["nS"], 2237,
               "the paper prints 2237 of 3920 for supp(phi) at L=8")
    rep.eq_int("support witness: |S*| at L=6 == (1/2+1/L)D", _swr[6]["nS"], 200)
    rep.truth("support witness: 426 of the 4900 ground-state amplitudes vanish",
              _swr[8]["ground_state_zeros"] == 426 and _swr[8]["ground_state_dimension"] == 4900
              and "$426$ of the $4900$ ground-state amplitudes" in _body,
              "%d of %d" % (_swr[8]["ground_state_zeros"], _swr[8]["ground_state_dimension"]),
              "the ground-state zero count printed is not the deposited one")
    _k0 = {L: read_json(root, "data/c3_frontier/adversarial/z_suppK_sym_L%d.json" % L)
           for L in (10, 12, 14)}
    rep.eq("support floor at L=10 == 0.600",
           _k0[10]["union_count_by_k"]["1e-12"][0] / _k0[10]["sector_dimension"], 0.600, tol=5e-4)
    _lo = min(_k0[L]["union_count_by_k"][th][0] / _k0[L]["sector_dimension"]
              for L in (12, 14) for th in ("1e-10", "1e-12"))
    rep.truth("support floor above 0.54 at L=12, 14 at both thresholds", _lo > 0.54,
              "smallest |supp(phi)|/D = %.4f" % _lo,
              "the floor printed for L = 12, 14 is not supported: %.4f" % _lo)
    # ------------------------------------------------ the dequantization control (App. dequant)
    # (2026-09-25) data/gflow.json is written by src/gflow_dequant.py with every seed (ten seeds,
    # symmetry-pinned gauge, the companion review's ladder configuration).  The statistics are
    # RE-DERIVED here from the per-seed values, and every number the appendix and Sec. IX print
    # is asserted.
    _gf = read_json(root, "data/gflow.json")["runs"]
    _gtxt = " ".join(open(os.path.join(root, "paper", "app_dequant.tex"), encoding="utf-8")
                     .read().split())
    _s9 = " ".join(open(os.path.join(root, "paper", "sec_9.tex"), encoding="utf-8").read().split())

    def _pt(a, b):
        _d = np.asarray(a, float) - np.asarray(b, float)
        return (_d.mean(), _d.std(ddof=1), _d.mean() / (_d.std(ddof=1) / np.sqrt(_d.size)),
                int((_d > 0).sum()))

    _A = {k: np.asarray(v["per_seed"], float) for k, v in _gf["1000"]["methods"].items()}
    _B = {k: np.asarray(v["per_seed"], float) for k, v in _gf["500"]["methods"].items()}
    rep.eq_int("gflow.json: ten paired seeds at 1000 and at 500 shots",
               len({len(v) for v in list(_A.values()) + list(_B.values())} | {10}), 1,
               "every method at both shot counts must carry the same ten seeds")
    _cl, _ge, _gn = (_pt(_A["ibm"], _A["ibm+cheap"]), _pt(_A["ibm+cheap"], _A["gfn-fused"]),
                     _pt(_A["ibm+cheap"], _A["gfn-nocheap"]))
    _cl5, _ge5, _gn5 = (_pt(_B["ibm"], _B["ibm+cheap"]), _pt(_B["ibm+cheap"], _B["gfn-fused"]),
                        _pt(_B["ibm+cheap"], _B["gfn-nocheap"]))
    _ha = _pt(_A["ibm+cheap"][:5], _A["gfn-fused"][:5])
    _hb = _pt(_A["ibm+cheap"][5:], _A["gfn-fused"][5:])
    _need = ["$%.2f\\pm%.2f$" % (_A[k].mean(), _A[k].std(ddof=1))
             for k in ("ibm", "ibm+cheap", "gfn-fused", "gfn-nocheap")]
    _need += ["$%.2f$~mHa, a $%.0f\\%%$ reduction, in all ten seeds, with paired $t_9=%.2f$"
              % (_cl[0], 100 * _cl[0] / _A["ibm"].mean(), _cl[2]),
              "$%.0f\\%%$ with $t_9=%.1f$" % (100 * _cl5[0] / _B["ibm"].mean(), _cl5[2]),
              "a further $%.2f\\pm%.2f$~mHa" % (_ge[0], _ge[1]),
              "with paired $t_9=%.2f$" % _ge[2],
              "($%.2f$~mHa, $\\lvert t_9\\rvert=%.2f$)" % (abs(_gn[0]), abs(_gn[2])),
              "$%.2f$~mHa ($\\lvert t_9\\rvert=%.2f$," % (abs(_ge5[0]), abs(_ge5[2])),
              "$%.2f$~mHa ($\\lvert t_9\\rvert=%.1f$)" % (abs(_gn5[0]), abs(_gn5[2])),
              "by $%.2f\\pm%.2f$~mHa in every seed, with $t_4=%.2f$" % (_ha[0], _ha[1], _ha[2]),
              "$%.2f\\pm%.2f$~mHa with $t_4=%.2f$" % (_hb[0], _hb[1], _hb[2])]
    for _q in _need:
        rep.truth("App. dequant prints %s" % _q[:52], _q in _gtxt, _q,
                  "the appendix does not print the re-derived %s" % _q)
    rep.truth("App. dequant verdicts: classical 10/10; fused unresolved at 1000; worse at 500",
              _cl[3] == 10 and abs(_ge[2]) < 2.262 and _ge5[0] < 0 and abs(_ge5[2]) > 2.262
              and _ha[3] == 5 and _hb[0] < 0,
              "t9 critical 2.262 (two-sided 5%)",
              "a verbal verdict of App. dequant no longer matches data/gflow.json")
    _q9 = ("from $%.2f$ to $%.2f$~mHa (paired $t_9=%.2f$, all ten seeds)"
           % (_A["ibm"].mean(), _A["ibm+cheap"].mean(), _cl[2]))
    rep.truth("Sec. IX dequantization paragraph matches gflow.json", _q9 in _s9 and
              ("a further $%.2f$~mHa, a difference ten paired seeds do not resolve ($t_9=%.2f$)"
               % (_ge[0], _ge[2])) in _s9, _q9,
              "Sec. IX prints numbers that are not those of data/gflow.json")
    _ms = {r["L"]: r for r in
           read_json(root, "data/c3_frontier/adversarial/z_momseed_support.json")["results"]}
    _q = "support on $%d$ to $%d$ of the $%d$ determinants" % (
        _ms[8]["support_min"], _ms[8]["support_max"], _ms[8]["D"])
    rep.truth("momentum-seed support at L=8 (z_momseed_support.json)", _q in _body,
              _q, "the momentum-seed support printed in Sec. III is not the deposited %s" % _q)
    _ti = {b["L"]: {p["FR"]: p for p in b["points"]}
           for b in read_json(root, "data/c3_frontier/ties.json")}
    for _L, _want in ((6, "30.2"), (8, "42.0")):
        _p = _ti[_L][0.85]
        _got = "%.1f" % (100.0 * _p["out_of_supp_in_S"] / _p["S"])
        rep.truth("ties.json: share of S outside supp(phi) at L=%d, FR=0.85" % _L,
                  _got == _want and ("$%s\\%%$" % _want) in _body,
                  "%s%% (%d of %d)" % (_got, _p["out_of_supp_in_S"], _p["S"]),
                  "the out-of-support share printed (%s%%) is not the deposited %s%%"
                  % (_want, _got))

    uncovered = sorted(inputs - COVERED_FRAGMENTS)
    rep.truth("every native figure source the paper \\input's is opened here",
              not uncovered,
              "%d fragments \\input'ed, all covered" % len(inputs),
              "UNWATCHED FIGURE SOURCES: %s are \\input'ed by the manuscript and no check "
              "in this file opens them. Every literal in them could be edited without this "
              "guardian noticing." % uncovered,
              n=max(len(inputs), 1))

    # ---------------------------------------------------------------- 9.3  the ABSTRACT
    # Never opened before.  Eight numbers live here and this is the text that is indexed.
    mt = re.sub(r"\s+", " ", read_body(root))
    mF, mS, mX = ctx["mF"], ctx["mS"], ctx["mX"]
    hF, hS, hX = ctx["hF"], ctx["hS"], ctx["hX"]
    # The v3 abstract no longer prints the four pooled coefficients: it prints the
    # stratified statement that replaced them ("exactly -1 within each of six strata,
    # exact p=3.3e-13; both pooled coefficients are withdrawn").  So the four checks are
    # replaced by checks on what is actually printed, plus an assertion that the four
    # retracted coefficients stay out of the abstract.
    _hub_rows = read_json(root, "data/cost_vs_ent.json")["rows"]
    _st = {}
    for _r in _hub_rows:
        _st.setdefault((_r["L"], _r["N"]), []).append(_r)
    _within = [spearmanr([r["FAF"] for r in rs], [r["S"] for r in rs]).correlation
               for _, rs in sorted(_st.items())]
    # 2026-09-30 (reader-value reframe): the stratified statistics left the abstract; their
    # home is Sec. I (sec_1.tex) and Sec. VI (sec_6_body.tex), with the Fig. 6 caption
    # (figs/captions.tex).  The two checks below were named for the abstract but searched
    # the whole body, so they kept passing on the Sec. I wording; they are now named for
    # where the statement lives, the regex admits every printed form ('exact $p=', 'exact
    # permutation $p=', 'exact one-sided permutation $p='), at least TWO places must print
    # the value (Sec. I and Sec. VI), and every printed value must be the 3.35 rounding of
    # exact_perm_p_one_sided_F1 = (1/120)^6 = 3.349e-13: tolerance 5e-15 -> 1e-15, so the
    # 3.3e-13 rounding that Sec. I and the Fig. 6 caption carried until today fails.
    rep.truth("Sec. I: 'exactly -1' within each of six strata",
              re.search(r"exactly \$-1\$", mt) is not None
              and len(_within) == 6 and all(abs(v + 1.0) < 1e-12 for v in _within),
              "six strata, every one at rho = -1.000000",
              "Sec. I claims an exact -1 within every stratum; recomputation gives %s"
              % [round(v, 6) for v in _within], n=6)
    _sr = read_json(root, "data/stats_resource.json")["hubbard_stratified"]
    _ms = list(re.finditer(r"exact(?:\s+one-sided)?(?:\s+permutation)?\s+\$p=([\d.]+)"
                           + RBS + r"times10\^\{-(\d+)\}\$", mt))
    if len(_ms) < 2:
        rep.missing("Secs. I/VI exact permutation p",
                    "fewer than two places (Sec. I and Sec. VI) state the exact permutation "
                    "p-value that replaced the withdrawn pooled coefficients (%d found)" % len(_ms))
    else:
        for _i, _m in enumerate(_ms):
            _pp = float(_m.group(1)) * 10 ** (-int(_m.group(2)))
            rep.eq("Secs. I/VI exact permutation p (occurrence %d of %d)" % (_i + 1, len(_ms)),
                   _pp, float(_sr["exact_perm_p_one_sided_F1"]), 1e-15,
                   "a printed exact permutation p is not the 3.35e-13 rounding of "
                   "exact_perm_p_one_sided_F1 of data/stats_resource.json")
    # v3, 2026-09-19.  The manuscript no longer announces a withdrawal anywhere; it
    # states positively what the two pooled coefficients ARE -- Simpson reversals of the
    # six (L,N) strata -- and that neither is used.  The old check here asserted the
    # presence of the sentence "pooled coefficients are withdrawn"; it was named for the
    # abstract but searched read_body(), and it was in fact satisfied by a paragraph
    # heading in app_stats.tex.  It could not detect the failure that matters now: a
    # pooled coefficient quietly deleted, leaving Fig. 6(d) and App. F printing a number
    # the text never explains.  Both halves are asserted instead.
    _simpson_named = "Simpson reversal" in mt
    _declared_unused = ("Neither reading is used here" in mt
                        and "Why neither pooled coefficient is used" in mt)
    rep.truth("the two pooled coefficients are named Simpson reversals and declared unused",
              _simpson_named and _declared_unused,
              "the body says what the two pooled coefficients are, and that neither is used",
              "the manuscript stopped explaining the two pooled coefficients "
              "(named as a Simpson reversal: %s; declared unused: %s). Deleting them "
              "instead of explaining them leaves Fig. 6(d) and App. F printing numbers "
              "the text does not account for." % (_simpson_named, _declared_unused))
    abstract_rhos = []
    for pat, value, label in abstract_rhos:
        m = re.search(pat, mt)
        if m is None:
            rep.missing(label, "the abstract sentence carrying this correlation could not be "
                               "located (pattern %r); the guardian cannot protect a number it "
                               "cannot find" % pat)
        else:
            rep.eq(label, float(m.group(1)), round(value, 2), 5e-3,
                   "the ABSTRACT prints a Spearman rho that disagrees with the value "
                   "recomputed from data/resource_master.json (%.6f). This is the text that "
                   "is indexed by arXiv" % value)
    for pat, want, label in (
            (r"\$n=38\$", len(ctx["mol_pts"]), "abstract n=38 molecular points"),
            (r"\$n=30\$", len(ctx["hub_pts"]), "abstract n=30 Hubbard points")):
        m = re.search(pat, mt)
        if m is None:
            rep.missing(label, "the abstract no longer states this sample size")
        else:
            rep.eq_int(label, want, int(re.search(r"\d+", m.group(0)).group(0)),
                       "the abstract states a sample size that is not the number of points "
                       "of that family in data/resource_master.json")

    # ---------------------------------------------------------------- 9.4
    # results.tex, matched-magic chain-vs-ladder, in PROSE.  The previous guardian
    # checked the caption thirty-six lines below and missed this one.
    rt = re.sub(r"\s+", " ", read_body(root))   # was paper/results.tex (v2)
    lvc = ctx["lvc"]
    m = re.search(r"\(\$([\d.]+)\$ vs \$([\d.]+)\$ at \$16\$ qubits", rt)
    if m is None:
        # v3 rewrote the sentence: it no longer prints the two F1 literals, it states the
        # relative match ("at one-body magic matched to within 1%").  Check THAT instead,
        # against the same two deposited numbers, so the claim is still pinned.
        _mm = re.search(r"one-body magic matched to within \$(\d+)" + RBS + r"%\$", rt)
        if _mm is None:
            rep.missing("matched-F1 chain vs ladder (prose)",
                        "neither the v2 sentence ('$12.82$ vs $12.91$ at $16$ qubits') nor "
                        "the v3 one ('one-body magic matched to within N%') is in the body; "
                        "the matched-magic comparison is unpinned")
        else:
            _tol = float(_mm.group(1)) / 100.0
            # inside section9 the two tables arrive through ctx, not as globals
            _ch = ctx["obc"].get(8)
            _la = lvc["results"].get("ladder-2x4", {}).get("F1")
            if _ch is None or _la is None:
                rep.missing("matched-F1 chain vs ladder (v3 prose)",
                            "the deposited chain or ladder F1 at 16 qubits is missing")
            else:
                _rel = abs(_la - _ch) / _ch
                rep.truth("matched-F1 chain vs ladder, within the stated %s%%"
                          % _mm.group(1), _rel <= _tol,
                          "|%.4f - %.4f| / %.4f = %.4f <= %.2f"
                          % (_la, _ch, _ch, _rel, _tol),
                          "the body says the two are matched to within %s%%; the deposited "
                          "values differ by %.2f%%" % (_mm.group(1), 100 * _rel))
    else:
        system_check(rep, "results.tex matched-F1 prose, chain L=8", float(m.group(1)),
                     8, "OBC", obc, pbc,
                     where="paper/results.tex (matched-magic sentence, chain)")
        rep.eq("results.tex matched-F1 prose, ladder 2x4", float(m.group(2)),
               round(lvc["results"]["ladder-2x4"]["F1"], 2), 5e-3,
               "the printed ladder F1 disagrees with data/ladder_vs_chain.json")

    # ---------------------------------------------------------------- 9.5  fig_ladder
    lt = read(root, "paper/figs/fig_ladder_native.tex")
    lls = logical_lines(lt)
    R = lvc["results"]
    for tag, key_fmt, field, ndp in (("2-leg ladder", "ladder-2x%d", "chi", 0),
                                     ("chain", "chain-L%d", "chi", 0)):
        ln = [l for l in lls if ("addlegendentry{" + tag + "}") in l]
        if not ln:
            rep.missing("fig_ladder (a) %s chi series" % tag,
                        "the series is gone from the figure source")
            continue
        pts_ = series_coords(ln[0])
        rep.eq_int("fig_ladder (a) %s chi series point count" % tag, len(pts_), 3,
                   "a point was added to or removed from the plotted chi series")
        for x, y in pts_:
            L = int(round(x / 2))
            k = key_fmt % (L // 2 if "ladder" in key_fmt else L)
            rep.eq("fig_ladder (a) %s chi at 2L=%d" % (tag, int(x)), y, R[k][field], 0.0,
                   "the plotted minimal bond dimension disagrees with "
                   "data/ladder_vs_chain.json")
    # the annotated ratio
    m = re.search(r"\{\$" + RBS + r"times([\d.]+)\$" + RBS + r"}", lt)
    if m:
        rep.eq("fig_ladder (a) annotated chi ratio at 16 qubits", float(m.group(1)),
               round(R["ladder-2x4"]["chi"] / float(R["chain-L8"]["chi"]), 1), 5e-2,
               "the annotated ladder/chain bond-dimension ratio is not the ratio of the "
               "deposited chi values")
    # panel (b): the required (N+-1) fraction
    # Panel (b) of the rebuilt fragment carries no legend of its own (the panel-(a)
    # legend serves both), so the series are identified by the colour macro that draws
    # them -- ldLad for the ladder, ldChain for the chain -- and by being the ones with
    # exactly two points.
    for tag, style, keys in (("2-leg ladder", "ldLad", ("ladder-2x3", "ladder-2x4")),
                             ("chain", "ldChain", ("chain-L6", "chain-L8"))):
        ln = [l for l in lls if style in l and len(series_coords(l)) == 2]
        if not ln:
            rep.missing("fig_ladder (b) %s fraction series" % tag,
                        "no two-point %s series found in the fragment" % style)
            continue
        pts_ = series_coords(ln[0])
        rep.eq_int("fig_ladder (b) %s fraction point count" % tag, len(pts_), 2,
                   "a point was added to or removed from the plotted fraction series")
        for (x, y), k in zip(pts_, keys):
            rep.eq("fig_ladder (b) %s fraction at 2L=%d" % (tag, int(x)), y,
                   R[k]["frac"], 5.1e-4,
                   "the plotted (N+-1) fraction disagrees with data/ladder_vs_chain.json")

    # ---------------------------------------------------------------- 9.6  fig_akw_sampled
    # The flagship sampled figure, 3562 literals, never opened before.  Two things are
    # checked: that every plotted point comes from the deposited reconstruction, and that
    # the "sampled" curve is NOT a copy of the exact one -- a perfect reconstruction is
    # exactly what a forger would draw.
    sk = read_json(root, "data/sampled_akw_L8.json")
    at = read(root, "paper/figs/fig_akw_sampled_honest.tex")
    als = logical_lines(at)
    ex_lines = [l for l in als if "addplot[akExact,line width" in l]
    sm_lines = [l for l in als if "addplot[akSamp,line width" in l]
    ov = sk["overlay"]
    ovk = sorted(ov.keys(), key=float)
    rep.eq_int("fig_akw_sampled: exact curves plotted", len(ex_lines), len(ovk),
               "the number of exact A(k,w) curves in the figure is not the number of "
               "momenta in data/sampled_akw_L8.json['overlay']")
    rep.eq_int("fig_akw_sampled: sampled curves plotted", len(sm_lines), len(ovk),
               "the number of sampled A(k,w) curves in the figure is not the number of "
               "momenta in data/sampled_akw_L8.json['overlay']")
    for i, k in enumerate(ovk):
        if i >= len(ex_lines) or i >= len(sm_lines):
            break
        ex_pairs = series_coords(ex_lines[i])
        sm_pairs = series_coords(sm_lines[i])
        ok, off = curve_on_grid(rep, "fig_akw_sampled exact curve k=%s" % k, ex_pairs,
                                sk["wg"], ov[k]["exact"], 3,
                                "the plotted EXACT A(k,w) is not the deposited exact array")
        curve_on_grid(rep, "fig_akw_sampled sampled curve k=%s" % k, sm_pairs,
                      sk["wg"], ov[k]["sampled"], 3,
                      "the plotted SAMPLED A(k,w) is not the deposited reconstruction",
                      offset=off)
        ex_plot = [y for _, y in ex_pairs]
        sm_plot = [y for _, y in sm_pairs]
        rep.truth("fig_akw_sampled: the sampled curve at k=%s is not a copy of the exact one"
                  % k,
                  any(a != b for a, b in zip(ex_plot, sm_plot)),
                  "the two curves differ, as a real reconstruction must",
                  "FORGERY SIGNATURE: the plotted 'sampled (85%% sector)' curve at k=%s is "
                  "identical to the exact curve. A reconstruction from a strict subspace "
                  "cannot be exact; a figure that shows it is claiming something the "
                  "deposited rel-L1 of %.4g contradicts" % (k, sk["mean_relL1"]))
    # panel (b): per-k rel-L1 markers
    mk = [l for l in als if "addplot[akSamp,only marks" in l]
    if not mk:
        rep.missing("fig_akw_sampled (b) per-k rel-L1 markers", "the series is gone")
    else:
        pts_ = series_coords(mk[0])
        per_k = sk["per_k"]
        rep.eq_int("fig_akw_sampled (b) marker count", len(pts_), len(per_k),
                   "the number of plotted per-k rel-L1 markers is not the number of "
                   "momenta in data/sampled_akw_L8.json['per_k']")
        for x, y in pts_:
            key = "%.3f" % x
            if key not in per_k:
                rep.emit(FAIL, "fig_akw_sampled (b) marker at k/pi=%s" % key,
                         "FABRICATED POINT: k/pi=%s has no entry in "
                         "data/sampled_akw_L8.json['per_k']" % key, n=1)
            else:
                rep.eq("fig_akw_sampled (b) rel-L1 at k/pi=%s" % key, y,
                       per_k[key]["relL1"], 5e-5,
                       "the plotted per-k rel-L1 disagrees with the deposited value")
    for pat, value, label in (
            (r"mean per-k rel-L1 ([\d.]+)", sk["mean_relL1"], "header mean rel-L1"),
            (r"max ([\d.]+)\)", sk["max_relL1"], "header max rel-L1")):
        m = re.search(pat, at)
        if m:
            rep.eq("fig_akw_sampled %s" % label, float(m.group(1)),
                   round(value, len(m.group(1).split(".")[1])), 5e-5,
                   "the figure header advertises an accuracy that is not the deposited one "
                   "(%.6g)" % value)

    # ---------------------------------------------------------------- 9.6b Fig. 3 finite shots
    # Added 2026-09-28.  The finite-shot paragraph of the Fig. 3 caption (and its echoes in
    # Sec. II, Sec. VII, Sec. S1 and Sec. S8) had no deposited generator until
    # src/fig3_finite_shot.py was recovered and deposited; its two outputs are
    # data/fig3_finite_shot_T2.6e6.json and _T5.2e6.json.  Every printed number is recomputed
    # HERE from the 32 per-(channel, seed) rows, not read from the file's own summary, and the
    # summary is checked against the rows too.  The two files must also be scored on the SAME
    # window as Fig. 3(a,b) (data/sampled_akw_L8.json), and each must carry an exact
    # reproduction of the 2026-09-19 output on the former window.
    import statistics as _st
    fs = {}
    for _T, _tag in ((2600000, "2.6e6"), (5200000, "5.2e6")):
        _rel = "data/fig3_finite_shot_T%s.json" % _tag
        if not os.path.exists(os.path.join(root, _rel)):
            rep.missing("Fig. 3 finite-shot data %s" % _rel,
                        "the file behind the finite-shot paragraph of the Fig. 3 caption is gone")
            continue
        d_ = read_json(root, _rel)
        c_ = d_["config"]
        rows_ = d_["rows"]
        rel_ = [r["relL1"] for r in rows_]
        fa_ = [r["frac_add"] for r in rows_]
        fr_ = [r["frac_rem"] for r in rows_]
        rep.eq_int("%s: T" % _rel, c_["T"], _T, "the file is not at the budget its name states")
        rep.truth("%s: configuration of Fig. 3(a,b)" % _rel,
                  (c_["L"], c_["U"], c_["eta"], c_["K"], c_["dt"]) == (sk["L"], sk["U"], sk["eta"], 16, 0.5)
                  and c_["NL"] == 260 and c_["seeds"] == [7001, 7002, 7003, 7004],
                  "L=8, U=8, eta=0.18, K=16, dt=0.5, NL=260, seeds 7001-7004",
                  "the finite-shot run is not at the configuration of Fig. 3(a,b): %r" % c_)
        rep.truth("%s: scored on the window of Fig. 3(a,b)" % _rel,
                  c_["window"] == sk["window"] and c_["npts"] == len(sk["wg"]),
                  "%s x %d" % (c_["window"], c_["npts"]),
                  "the finite-shot errors are scored on %s x %d while Fig. 3(a,b) uses %s x %d; "
                  "their ratio would compare two different integrals"
                  % (c_["window"], c_["npts"], sk["window"], len(sk["wg"])))
        rep.eq_int("%s: 32 (channel, seed) rows" % _rel, len(rows_), 32,
                   "8 momenta x 4 seeds, each row both branches summed")
        rep.eq_int("%s: rows are the 8 momenta x 4 seeds" % _rel,
                   sorted((round(r["k_over_pi"], 3), r["seed"]) for r in rows_),
                   sorted((round(2 * n / 8.0, 3), s_) for n in range(8) for s_ in c_["seeds"]),
                   "a (momentum, seed) pair is missing or duplicated")
        rep.truth("%s: fractions are |S|/3920" % _rel,
                  all(abs(r["frac_add"] - r["S_add"] / 3920.0) < 1e-15 and
                      abs(r["frac_rem"] - r["S_rem"] / 3920.0) < 1e-15 for r in rows_),
                  "64 fractions", "a stored fraction is not its |S| over the sector dimension",
                  n=2 * len(rows_))
        st_ = dict(mean=_st.mean(rel_), sd=_st.stdev(rel_), mx=max(rel_),
                   fa=_st.mean(fa_), fr=_st.mean(fr_))
        sm_ = d_["summary"]
        for k_, key_ in (("mean", "mean_relL1"), ("sd", "sd_relL1_ddof1"), ("mx", "max_relL1"),
                         ("fa", "mean_frac_add"), ("fr", "mean_frac_rem")):
            rep.eq("%s: summary.%s recomputed from the rows" % (_rel, key_), sm_[key_], st_[k_],
                   1e-12 * max(1.0, abs(st_[k_])),
                   "the stored summary is not the statistic of the stored rows")
        rep.eq("%s: mean_relL1 (top level) from the rows" % _rel, d_["mean_relL1"], st_["mean"], 1e-15)
        rep.eq("%s: max_relL1 (top level) from the rows" % _rel, d_["max_relL1"], st_["mx"], 0.0)
        rc_ = d_.get("reproduction_check") or {}
        pv_ = rc_.get("printed_vs_reproduced") or {}
        rep.truth("%s: exact reproduction of the 2026-09-19 output on [-9t, 9t]" % _rel,
                  rc_.get("exact") is True and rc_.get("n_rows_mismatched") == 0
                  and rc_.get("n_rows_compared") == (32 if _T == 2600000 else 28)
                  and rc_.get("window") == [-9.0, 9.0] and rc_.get("npts") == 600
                  and len(pv_) == 5 and all(a == b for a, b in pv_.values()),
                  "%s rows and all summary strings identical" % rc_.get("n_rows_compared"),
                  "the deposited generator does not carry an exact reproduction of the run it "
                  "replaces: %r" % {k: rc_.get(k) for k in ("exact", "n_rows_compared", "n_rows_mismatched")})
        fs[_T] = dict(st_, T=c_["T"], n=len(rows_), nseeds=len(c_["seeds"]))
    if len(fs) == 2:
        sh_ = read_json(root, "data/sampled_honest.json")
        pen_ = [m_["penalty_factor_topm"] for r_ in sh_["results"] if r_["L"] == 8
                for m_ in r_["matched_subspace"] if m_.get("S_matched") == 2180]
        rep.eq_int("Table S6 transferred factor (L=8, |S|=2180) found once", len(pen_), 1)
        pen0 = pen_[0] if pen_ else float("nan")
        ratio_ = fs[2600000]["mean"] / sk["mean_relL1"]
        under_ = 1.0 - pen0 / ratio_
        # 2026-09-28 (late): the generator writes the understatement only at T = 2.6e6 (the
        # budget whose fraction, 0.839, is the one next to the 0.850 of Fig. 3(a,b)); at
        # T = 5.2e6 the field is null (it used to hold a meaningless -1.24).  Each data file
        # also carries the sha256 of its generator with line endings normalised to LF; it is
        # recomputed here from src/fig3_finite_shot.py, so an edit to the algorithm block
        # cannot leave a stored reproduction check green.
        with open(os.path.join(root, "src", "fig3_finite_shot.py"), "rb") as _fh:
            _gen_lf = hashlib.sha256(_fh.read().replace(b"\r\n", b"\n")).hexdigest()
        for _T, _tag in ((2600000, "2.6e6"), (5200000, "5.2e6")):
            _d = read_json(root, "data/fig3_finite_shot_T%s.json" % _tag)
            sm_, pv_ = _d["summary"], _d.get("provenance") or {}
            rep.eq("fig3_finite_shot_T%s: ratio to Fig. 3(b) mean from the rows" % _tag,
                   sm_["ratio_mean_to_panel_b"], fs[_T]["mean"] / sk["mean_relL1"], 1e-12)
            rep.eq("fig3_finite_shot_T%s: ratio to Fig. 3(b) max from the rows" % _tag,
                   sm_["ratio_max_to_panel_b"], fs[_T]["mx"] / sk["max_relL1"], 1e-12)
            if _T == 2600000:
                rep.eq("fig3_finite_shot_T%s: understatement from the rows and Table S6" % _tag,
                       sm_["understatement_of_transferred_penalty"], under_, 1e-12)
            else:
                rep.truth("fig3_finite_shot_T%s: no understatement is written at this budget" % _tag,
                          sm_.get("understatement_of_transferred_penalty") is None, "null",
                          "a value is stored for a quantity the generator defines only at T=2.6e6")
            rep.truth("fig3_finite_shot_T%s: provenance.generator_sha256 is the LF hash of "
                      "src/fig3_finite_shot.py" % _tag, pv_.get("generator_sha256") == _gen_lf,
                      _gen_lf[:12], "the deposited data were not written by the deposited generator "
                      "(stored %s, file %s)" % (str(pv_.get("generator_sha256"))[:12], _gen_lf[:12]))
            rep.truth("fig3_finite_shot_T%s: provenance.argv carries no absolute path" % _tag,
                      all(not os.path.isabs(str(a)) for a in pv_.get("argv", [])), "relative",
                      "a local absolute path is recorded in the deposit")
        cap = " ".join(read(root, "paper/figs/fig5_caption.tex").split())
        f5, f26 = fs[5200000], fs[2600000]

        def _num(pat, where, text, name, want, tol, defect):
            m_ = re.search(pat, text)
            if not m_:
                rep.missing(name, "the printed %s is gone from %s (pattern %s)" % (name, where, pat))
                return
            rep.eq(name, float(m_.group(1)), want, tol, defect)

        D = "the caption no longer states what data/fig3_finite_shot_T*.json holds"
        m_ = re.search(r"\$T=([\d.]+)\\times10\^\{6\}\$ shots per channel---\$([\d.]+)\\times10\^\{7\}\$ "
                       r"over the (\d+) channels", cap)
        if not m_:
            rep.missing("Fig. 3 caption: budget sentence",
                        "'T=... shots per channel---... over the 16 channels' is gone")
        else:
            rep.eq("Fig. 3 caption: T per channel (x1e6)", float(m_.group(1)), f5["T"] / 1e6, 0.0, D)
            rep.eq("Fig. 3 caption: shots over the channels (x1e7)", float(m_.group(2)),
                   round(int(m_.group(3)) * f5["T"] / 1e7, 1), 1e-12, D)
            rep.eq_int("Fig. 3 caption: number of channels", int(m_.group(3)), 16, D)
        _num(r"reaches rel-\$L_1=\(([\d.]+)\\pm[\d.]+\)\\times10\^\{-3\}\$, maximum", "caption", cap,
             "Fig. 3 caption: T=5.2e6 mean rel-L1 (x1e-3)", round(f5["mean"] * 1e3, 2), 1e-12, D)
        _num(r"reaches rel-\$L_1=\([\d.]+\\pm([\d.]+)\)\\times10\^\{-3\}\$, maximum", "caption", cap,
             "Fig. 3 caption: T=5.2e6 s.d. (x1e-3)", round(f5["sd"] * 1e3, 2), 1e-12, D)
        _num(r"maximum \$([\d.]+)\\times10\^\{-3\}\$: panel", "caption", cap,
             "Fig. 3 caption: T=5.2e6 max rel-L1 (x1e-3)", round(f5["mx"] * 1e3, 2), 1e-12, D)
        _num(r"covering \$([\d.]+)\$ of the addition sector", "caption", cap,
             "Fig. 3 caption: T=5.2e6 addition fraction", round(f5["fa"], 3), 1e-12, D)
        _num(r"against the \$([\d.]+)\$ of \(a\)--\(b\)", "caption", cap,
             "Fig. 3 caption: fraction of (a)-(b)", round(sk["frac"], 3), 1e-12, D)
        # 2026-09-30 (final numeric review NR-1): the widened window moved the panel (b) mean
        # and maximum of the caption (2.19e-3 -> 2.36e-3, 3.28e-3 -> 3.27e-3) and added the two
        # branch thresholds and the window itself; until now only the figure header was tied to
        # data/sampled_akw_L8.json, so restoring the stale caption values left every guardian
        # green.  Each printed value is now its rounding of the deposited one.
        DB = "the Fig. 3 caption, panel (b), does not state what data/sampled_akw_L8.json holds"
        _num(r"removal branch summed over that window: mean \$([\d.]+)\\times10\^\{-3\}\$", "caption", cap,
             "Fig. 3 caption: panel (b) mean rel-L1 (x1e-3)", round(sk["mean_relL1"] * 1e3, 2), 1e-12, DB)
        _num(r"removal branch summed over that window: mean \$[\d.]+\\times10\^\{-3\}\$, maximum "
             r"\$([\d.]+)\\times10\^\{-3\}\$\.", "caption", cap,
             "Fig. 3 caption: panel (b) max rel-L1 (x1e-3)", round(sk["max_relL1"] * 1e3, 2), 1e-12, DB)
        m_ = re.search(r"window, \$\\omega-E_0\\in\[(-?[\d.]+)\\,t,(-?[\d.]+)\\,t\]\$, is symmetric", cap)
        if not m_:
            rep.missing("Fig. 3 caption: plotted and scored window", "'[-9 t,17 t]' is gone")
        else:
            rep.eq("Fig. 3 caption: window lower edge", float(m_.group(1)), float(sk["window"][0]), 0.0, DB)
            rep.eq("Fig. 3 caption: window upper edge", float(m_.group(2)), float(sk["window"][1]), 0.0, DB)
        _cg8 = [r_ for r_ in read_json(root, "data/charge_gap.json")["results"] if r_["L"] == sk["L"]]
        if len(_cg8) != 1:
            rep.missing("Fig. 3 caption: branch thresholds",
                        "data/charge_gap.json has no single row at L=%s" % sk["L"])
        else:
            DG = "the Fig. 3 caption's branch threshold is not data/charge_gap.json at L=8"
            _num(r"removal branch lies below \$\\omega-E_0=([\d.]+)\\,t\$", "caption", cap,
                 "Fig. 3 caption: removal threshold mu_- (t)", round(_cg8[0]["mu_minus"], 2), 1e-12, DG)
            _num(r"addition branch above \$([\d.]+)\\,t\$", "caption", cap,
                 "Fig. 3 caption: addition threshold mu_+ (t)", round(_cg8[0]["mu_plus"], 2), 1e-12, DG)
            rep.truth("Fig. 3 caption: the window holds both branches",
                      sk["window"][0] < _cg8[0]["mu_minus"] and sk["window"][1] > _cg8[0]["mu_plus"],
                      "holds", "the scored window cuts a branch the caption says it holds")
        _num(r"At half the budget, \$T=([\d.]+)\\times10\^\{6\}\$", "caption", cap,
             "Fig. 3 caption: half budget (x1e6)", f26["T"] / 1e6, 0.0, D)
        _num(r"subspace covers \$([\d.]+)\$---just under", "caption", cap,
             "Fig. 3 caption: T=2.6e6 addition fraction", round(f26["fa"], 3), 1e-12, D)
        _num(r"error is \$\(([\d.]+)\\pm[\d.]+\)\\times10\^\{-3\}\$", "caption", cap,
             "Fig. 3 caption: T=2.6e6 mean rel-L1 (x1e-3)", round(f26["mean"] * 1e3, 2), 1e-12, D)
        _num(r"error is \$\([\d.]+\\pm([\d.]+)\)\\times10\^\{-3\}\$", "caption", cap,
             "Fig. 3 caption: T=2.6e6 s.d. (x1e-3)", round(f26["sd"] * 1e3, 2), 1e-12, D)
        _num(r"a factor \$([\d.]+)\$ above the mean of \(b\)", "caption", cap,
             "Fig. 3 caption: ratio to the mean of (b)", round(ratio_, 2), 1e-12, D)
        _num(r"understated it by \$(\d+)\\%\$", "caption", cap,
             "Fig. 3 caption: understatement of the transferred factor (%)",
             round(100 * under_), 1e-12, D)
        m_ = re.search(r"\(\$1-([\d.]+)/([\d.]+)\$\)", cap)
        if not m_:
            rep.missing("Fig. 3 caption: definition of the understatement",
                        "the explicit '(1-2.05/ratio)' is gone")
        else:
            rep.eq("Fig. 3 caption: transferred factor in the definition", float(m_.group(1)),
                   round(pen0, 2), 1e-12, D)
            rep.eq("Fig. 3 caption: measured ratio in the definition", float(m_.group(2)),
                   round(ratio_, 2), 1e-12, D)
        # 2026-09-28 (late): the caption counts '(momentum, seed) pairs' (8 momenta x 4 seeds,
        # both branches summed in each), no longer '(channel, seed)', which read as 16 x 4.
        _num(r"(\d+) seeds, \$\d+\$ \(momentum, seed\) pairs, both branches summed", "caption", cap,
             "Fig. 3 caption: seed count", f5["nseeds"], 0, D)
        _num(r"seeds, \$(\d+)\$ \(momentum, seed\) pairs, both branches summed", "caption", cap,
             "Fig. 3 caption: pair count", f5["n"], 0, D)
        rep.eq_int("Fig. 3 caption: pairs = momenta x seeds", f5["n"],
                   len(set(r["k_over_pi"] for r in rows_)) * f5["nseeds"],
                   "the row count is not 8 momenta x 4 seeds")
        smm = " ".join(read(root, "paper/sm_method.tex").split())
        _num(r"penalty is a measured \$([\d.]+)\$ and not", "sm_method.tex", smm,
             "Sec. S1: measured finite-shot penalty", round(ratio_, 2), 1e-12,
             "Sec. S1 quotes a penalty that is not T=2.6e6 mean / Fig. 3(b) mean")
        # 2026-09-28 (late): 5.2e6 is the budget that SUFFICES (the two budgets measured are a
        # factor two apart), and 2.61 is measured at fraction 0.839, 'just under' 0.850.
        _num(r"a factor two apart: \$([\d.]+)\\times10\^\{6\}\$ per channel suffice", "sm_method.tex", smm,
             "Sec. S1: T per channel (x1e6)", f5["T"] / 1e6, 0.0, "Sec. S1 quotes another budget")
        _num(r"just under the announced one \(\$([\d.]+)\$ against \$[\d.]+\$\), the penalty",
             "sm_method.tex", smm, "Sec. S1: fraction at which 2.61 is measured",
             round(f26["fa"], 3), 1e-12, "Sec. S1 states the wrong fraction for the measured penalty")
        _num(r"just under the announced one \(\$[\d.]+\$ against \$([\d.]+)\$\), the penalty",
             "sm_method.tex", smm, "Sec. S1: the announced fraction", round(sk["frac"], 3), 1e-12,
             "Sec. S1 states the wrong announced fraction")
        _num(r"not the \$([\d.]+)\$ Table~\\ref\{tab:finiteshot\} would transfer", "sm_method.tex", smm,
             "Sec. S1: transferred factor", round(pen0, 2), 1e-12, "not Table S6's L=8 factor")
        s2 = " ".join(read(root, "paper/sec_2.tex").split())
        _num(r"at finite shots \$([\d.]+)\\times10\^\{6\}\$ per channel suffice and \$[\d.]+\\times10\^\{6\}\$ do not",
             "sec_2.tex", s2, "Table I: T per channel that suffices (x1e6)", f5["T"] / 1e6, 0.0,
             "Table I quotes a budget that is not the deposited run's")
        _num(r"at finite shots \$[\d.]+\\times10\^\{6\}\$ per channel suffice and \$([\d.]+)\\times10\^\{6\}\$ do not",
             "sec_2.tex", s2, "Table I: T per channel that does not (x1e6)", f26["T"] / 1e6, 0.0,
             "Table I quotes a budget that is not the deposited run's")
        _num(r"panel~\(b\) at \$([\d.]+)\\times10\^\{6\}\$ shots per channel and not at half that budget",
             "sec_2.tex", s2, "Sec. II B text: T per channel (x1e6)", f5["T"] / 1e6, 0.0,
             "Sec. II quotes a budget that is not the deposited run's")
        for rel__, pat_ in (("paper/sec_7.tex", r"accuracy with \$([\d.]+)\\times10\^\{7\}\$, and not with half as many"),
                            ("paper/sm_scope.tex", r"channels, \$([\d.]+)\\times10\^\{7\}\$ shots as measured")):
            _num(pat_, rel__, " ".join(read(root, rel__).split()),
                 "%s: total shots of the Fig. 3 run (x1e7)" % rel__, round(16 * f5["T"] / 1e7, 1), 1e-12,
                 "the total shot count is not 16 channels x the deposited T")
        # 2026-09-28 (late): Sec. S8 quotes the 8.3e7 with the accuracy it bought and says the
        # cost at the 0.05 threshold was not measured (half the budget gives 6.16e-3).
        sc_ = " ".join(read(root, "paper/sm_scope.tex").split())
        _num(r"to reach its relative \$L_1\$ of \$([\d.]+)\\times10\^\{-3\}\$", "sm_scope.tex", sc_,
             "Sec. S8: accuracy bought by 8.3e7 shots (x1e-3)", round(f5["mean"] * 1e3, 2), 1e-12,
             "Sec. S8 quotes an accuracy that is not the T=5.2e6 mean")
        _num(r"half that budget already gives \$([\d.]+)\\times10\^\{-3\}\$", "sm_scope.tex", sc_,
             "Sec. S8: accuracy at half the budget (x1e-3)", round(f26["mean"] * 1e3, 2), 1e-12,
             "Sec. S8 quotes an accuracy that is not the T=2.6e6 mean")
        rep.truth("Sec. S8: the cost at the 0.05 threshold is declared unmeasured",
                  "which was not measured" in sc_, "stated",
                  "Sec. S8 presents 8.3e7 as the cost of the observable without its accuracy qualifier")
        body_ = " ".join(read_body(root).split())
        for stale_ in (r"narrower window $\omega-E_0\in[-9\,t,9\,t]$",
                       "have not been recomputed on the window"):
            rep.truth("no stale finite-shot window disclosure: %s" % stale_[:40], stale_ not in body_,
                      "absent", "the manuscript still says the finite-shot numbers are on the old "
                      "window, but data/fig3_finite_shot_T*.json are scored on the window of (a)-(b)")

    # ---------------------------------------------------------------- 9.6c  three generators of 2026-09-28 (late)
    # (i)   the free-fermion determinant supports 35 / 336 / 3496 / 37361 of Secs. VI and S13,
    #       hand-typed constants with no producer until src/free_fermion_support.py wrote
    #       data/free_fermion_support.json.  Recomputed here from the open-chain orbitals with
    #       code that shares nothing with that script (Cauchy-Binet Slater weights, cumulative
    #       cut, inclusive tie convention).
    # (ii)  the spin-peak deviations of the Fig. 8 caption (src/spin_peak_deviation.py ->
    #       data/spin_peak_deviation.json), recomputed from data/spinqw_L12.json.
    # (iii) the pooled sweep of Fig. S5, rebuilt by src/recovered/build_n3.py from the DEPOSITED
    #       data/cert_{stress,akw,teqsci}.json (it read an earlier run of the same scans until
    #       2026-09-28).  The pool is re-classified here from the five JSON sources, the eight
    #       plotted n3_*.dat tables and the whisker summary are checked against that
    #       classification, and so is every count the caption, Sec. III D, Sec. S2 and Sec. S9
    #       print.  This closes the COVERAGE_GAP_V3 entry of 2026-09-19.
    def _pnum(pat, where, text, name, want, tol, defect):
        m_ = re.search(pat, text)
        if not m_:
            rep.missing(name, "the printed %s is gone from %s (pattern %s)" % (name, where, pat))
            return None
        rep.eq(name, float(m_.group(1)), want, tol, defect)
        return float(m_.group(1))

    def _tex(rel):
        """A .tex file as one line, comment lines dropped (captions carry numbers in comments)."""
        return " ".join(" ".join(ln for ln in read(root, rel).splitlines()
                                 if not ln.lstrip().startswith("%")).split())

    # ---- (i) free-fermion supports ----------------------------------------------------------
    try:
        ffs = read_json(root, "data/free_fermion_support.json")
    except Exception as exc:                                    # noqa: BLE001
        rep.missing("data/free_fermion_support.json", "cannot read: %s" % exc)
        ffs = None
    if ffs:
        _eps = float(ffs["eps_w"])
        rep.eq("free_fermion_support: eps_w is the lattice threshold eps^2 at eps = 0.05", _eps, 0.05 ** 2, 1e-15,
               "the discarded-weight threshold is not the 2.5e-3 the text prints")
        _ffrows = {int(r["L"]): r for r in ffs["rows"]}
        _ffgot = {}
        for L in (4, 6, 8, 10):
            jj = np.arange(1, L + 1)[:, None]
            mm = np.arange(1, L + 1)[None, :]
            Phi = np.sqrt(2.0 / (L + 1)) * np.sin(np.pi * jj * mm / (L + 1))
            Phi = Phi[:, np.argsort(-2.0 * np.cos(np.pi * np.arange(1, L + 1) / (L + 1)))][:, :L // 2]
            w1 = np.array([np.linalg.det(Phi[list(I), :]) ** 2 for I in combinations(range(L), L // 2)])
            w = np.sort(np.outer(w1, w1).ravel())[::-1]
            w = w / w.sum()
            tail = np.concatenate([np.cumsum(w[::-1])[::-1], [0.0]])    # tail[k] = sum_{i>=k} w_i
            n_in = int(np.argmax(tail <= _eps * (1 + 1e-12)))
            n_st = int(np.argmax(tail < _eps * (1 - 1e-12)))
            _ffgot[L] = (n_in, n_st, w, tail)
            r = _ffrows.get(L)
            if r is None:
                rep.missing("free_fermion_support L=%d" % L, "row missing from the deposited file")
                continue
            rep.eq("free_fermion_support L=%d: Cauchy-Binet norm of the Slater weights" % L,
                   float(w1.sum()), 1.0, 1e-12, "the up-spin determinant weights do not sum to one")
            rep.eq_int("free_fermion_support L=%d: sector dimension" % L, int(r["dim"]), int(w.size),
                       "the deposited dimension is not C(L,L/2)^2")
            rep.eq_int("free_fermion_support L=%d: |S|_eps_w, inclusive cut, recomputed" % L,
                       int(r["support_inclusive"]), n_in,
                       "the deposited inclusive support is not the smallest set whose discarded weight is <= eps_w")
            rep.eq_int("free_fermion_support L=%d: |S|_eps_w, strict cut, recomputed" % L,
                       int(r["support_strict"]), n_st,
                       "the deposited strict support is not the smallest set whose discarded weight is < eps_w")
            rep.eq("free_fermion_support L=%d: discarded weight at the inclusive cut" % L,
                   float(r["discarded_weight_inclusive"]), float(tail[n_in]), 1e-12,
                   "the recorded discarded weight is not the tail sum at the cut")
            rep.truth("free_fermion_support L=%d: discarded weight at the cut <= eps_w" % L,
                      tail[n_in] <= _eps * (1 + 1e-12), "%.3e" % tail[n_in],
                      "the kept set discards more than eps_w")
            rep.eq_int("free_fermion_support L=%d: published value == inclusive count" % L,
                       int(ffs["published"][str(L)]), n_in, "the value the paper prints is not the inclusive count")
        if all(L in _ffgot for L in (4, 6, 8, 10)):
            w4 = _ffgot[4][2]
            rep.truth("free_fermion_support L=4: exactly four weights equal eps_w = 1/400",
                      bool(np.all(np.abs(w4[-4:] - _eps) < 1e-15) and abs(w4[-5] - _eps) > 1e-6),
                      "1/400 x 4", "the L=4 tie the text declares is not in the recomputed weights")
            rep.eq_int("free_fermion_support L=4: the strict cut keeps one more determinant",
                       _ffgot[4][1], _ffgot[4][0] + 1, "the L=4 tie no longer separates the two conventions")
            rep.eq_int("free_fermion_support L=4: tie count recorded", int(_ffrows[4]["n_weights_tied_at_boundary"]), 4,
                       "the deposited tie count at L=4 is not four")
            for L in (6, 8, 10):
                rep.eq_int("free_fermion_support L=%d: inclusive and strict cuts agree" % L,
                           _ffgot[L][0], _ffgot[L][1], "a tie at the cut appeared where the text says there is none")
        for c_ in ffs.get("ed_crosscheck", []):
            rep.truth("free_fermion_support L=%s: ED cross-check recorded below 1e-12" % c_.get("L"),
                      float(c_.get("max_abs_weight_diff", 1.0)) < 1e-12, "%.1e" % float(c_.get("max_abs_weight_diff", 1.0)),
                      "the determinant weights do not match the U=0 exact diagonalization")
        rep.truth("free_fermion_support: the definition is the inclusive cumulative cut",
                  "inclusive" in str(ffs.get("definition", "")), "stated",
                  "the tie convention is not written into the deposit")
        smn = _tex("paper/sm_nogo.tex")
        m_ = re.search(r"\|\\mathcal S\|_\{\\varepsilon_w\}=(\d+),\\,(\d+),\\,(\d+),\\,(\d+)\$ for \$L=4,6,8,10\$ "
                       r"on the open chain at the lattice threshold \$\\varepsilon_w=([\d.]+)\\times10\^\{-3\}\$", smn)
        if not m_:
            rep.missing("Sec. S13: free-fermion supports",
                        "the sentence '|S|_eps_w = 35, 336, 3496, 37361 for L=4,6,8,10 on the open chain' is gone")
        else:
            for i_, L in enumerate((4, 6, 8, 10)):
                rep.eq_int("Sec. S13: |S|_eps_w at L=%d" % L, int(m_.group(i_ + 1)), _ffgot[L][0],
                           "Sec. S13 prints a support that is not the recomputed inclusive count")
            rep.eq("Sec. S13: eps_w", float(m_.group(5)) * 1e-3, _eps, 1e-15, "Sec. S13 prints another threshold")
        _pnum(r"a strict cut gives \$(\d+)\$", "sm_nogo.tex", smn, "Sec. S13: strict-cut value at L=4",
              _ffgot[4][1], 0, "the strict-cut value at L=4 is not the recomputed one")
        rep.truth("Sec. S13: the L=4 tie is declared",
                  "the four smallest determinant weights equal $\\varepsilon_w$ exactly" in smn, "stated",
                  "the tie convention behind 35 is no longer stated where the number is printed")
        rep.truth("Sec. S13: the generator is named", "free\\_fermion\\_support.py" in smn, "named",
                  "the four supports are printed without their generator")
        s6_ = _tex("paper/sec_6_body.tex")
        m_ = re.search(r"\(open chain, \$\\varepsilon_w=([\d.]+)\\times10\^\{-3\}\$\) of \$(\d+)\$, \$(\d+)\$, "
                       r"\$(\d+)\$ and \$(\d+)\$ at \$L=4\$--\$10\$", s6_)
        if not m_:
            rep.missing("Sec. VI: free-fermion supports",
                        "the sentence '(open chain, eps_w=2.5e-3) of 35, 336, 3496 and 37361 at L=4--10' is gone")
        else:
            rep.eq("Sec. VI: eps_w", float(m_.group(1)) * 1e-3, _eps, 1e-15, "Sec. VI prints another threshold")
            for i_, L in enumerate((4, 6, 8, 10)):
                rep.eq_int("Sec. VI: |S|_eps_w at L=%d" % L, int(m_.group(i_ + 2)), _ffgot[L][0],
                           "Sec. VI prints a support that is not the recomputed inclusive count")
        rep.truth("no 'no generator is deposited' disclaimer remains in the manuscript",
                  "no generator for these four values" not in " ".join(read_body(root).split()),
                  "absent", "the manuscript still publishes the four supports as generator-less constants")

    # ---- (ii) spin-peak deviations of the Fig. 8 caption -------------------------------------
    try:
        spd = read_json(root, "data/spin_peak_deviation.json")
        sqz = read_json(root, "data/spinqw_L12.json")
    except Exception as exc:                                    # noqa: BLE001
        rep.missing("data/spin_peak_deviation.json", "cannot read: %s" % exc)
        spd = None
    if spd:
        wg_ = np.array(sqz["wg"], dtype=float)
        J_ = float(sqz["J"])
        got_ = {}
        for key_, A_ in sqz["S"].items():
            q_ = float(key_)
            pk_ = float(wg_[int(np.argmax(np.array(A_, dtype=float)))])
            dcp_ = 0.5 * np.pi * J_ * abs(np.sin(np.pi * q_))
            got_[round(q_, 5)] = (pk_, dcp_, (pk_ - dcp_) / dcp_ if dcp_ > 1e-12 else None)
        rows_ = {round(float(r["q_over_pi"]), 5): r for r in spd["rows"]}
        rep.eq_int("spin_peak_deviation: one row per momentum of spinqw_L12.json", len(rows_), len(got_),
                   "the deposited rows are not the momenta of data/spinqw_L12.json")
        for q_, (pk_, dcp_, dev_) in sorted(got_.items()):
            r = rows_.get(q_)
            if r is None:
                rep.missing("spin_peak_deviation q/pi=%.5f" % q_, "row missing")
                continue
            rep.eq("spin_peak_deviation q/pi=%.5f: grid-argmax peak" % q_, float(r["peak"]), pk_, 1e-12,
                   "the deposited peak is not the argmax of the broadened line on its own grid")
            rep.eq("spin_peak_deviation q/pi=%.5f: des Cloizeaux-Pearson boundary" % q_, float(r["dCP_lower"]), dcp_, 1e-12,
                   "the deposited boundary is not (pi J/2)|sin q|")
            if dev_ is None:
                rep.truth("spin_peak_deviation q/pi=%.5f: no deviation where the boundary vanishes" % q_,
                          r["rel_dev"] is None, "null", "a deviation is stored at q=pi")
            else:
                rep.eq("spin_peak_deviation q/pi=%.5f: relative deviation" % q_, float(r["rel_dev"]), dev_, 1e-12,
                       "the deposited deviation is not (peak - boundary)/boundary")
        le_ = max(abs(d_) for q_, (_, _, d_) in got_.items() if d_ is not None and (q_ <= 0.5 + 1e-9 or q_ >= 1.5 - 1e-9))
        d23 = got_[round(2.0 / 3.0, 5)][2]
        d56 = got_[round(5.0 / 6.0, 5)][2]
        sp_ = spd["summary"]
        rep.eq("spin_peak_deviation summary: max |dev| for q <= pi/2", float(sp_["max_abs_rel_dev_q_le_half_pi"]), le_, 1e-12,
               "the summary maximum is not the maximum over the rows")
        rep.eq("spin_peak_deviation summary: deviation at 2pi/3", float(sp_["rel_dev_2pi_3"]), d23, 1e-12,
               "the summary value at 2pi/3 is not the row's")
        rep.eq("spin_peak_deviation summary: deviation at 5pi/6", float(sp_["rel_dev_5pi_6"]), d56, 1e-12,
               "the summary value at 5pi/6 is not the row's")
        rep.eq("spin_peak_deviation summary: grid spacing", float(sp_["grid_spacing"]), float(wg_[1] - wg_[0]), 1e-15,
               "the recorded grid spacing is not that of spinqw_L12.json")
        rep.eq("spin_peak_deviation summary: J", float(sp_["J"]), J_, 0.0, "another J")
        fsp = _tex("paper/carried/fig_spin.tex")
        v_ = _pnum(r"to within about \$(\d+)\\%\$ for \$q\\le\\pi/2\$", "fig_spin.tex", fsp,
                   "Fig. 8 caption: bound on |dev| for q <= pi/2 (%)", float(round(100 * le_)), 0.0,
                   "the caption's 'about N%' is not the rounded maximum deviation for q <= pi/2")
        if v_ is not None:
            rep.truth("Fig. 8 caption: 'about %d%%' holds to half a point" % int(v_), 100 * le_ <= v_ + 0.5,
                      "%.2f%%" % (100 * le_), "the printed percentage understates the measured deviation")
        _pnum(r"by \$(\d+)\\%\$ at \$q=2\\pi/3\$", "fig_spin.tex", fsp, "Fig. 8 caption: deviation at 2pi/3 (%)",
              float(round(100 * d23)), 0.0, "the caption's percentage at 2pi/3 is not the recomputed one")
        _pnum(r"and \$(\d+)\\%\$ at \$q=5\\pi/6\$", "fig_spin.tex", fsp, "Fig. 8 caption: deviation at 5pi/6 (%)",
              float(round(100 * d56)), 0.0, "the caption's percentage at 5pi/6 is not the recomputed one")
        _pnum(r"on the \$([\d.]+)\\,t\$ frequency grid", "fig_spin.tex", fsp, "Fig. 8 caption: frequency grid (t)",
              round(float(wg_[1] - wg_[0]), 3), 1e-12, "the caption's grid spacing is not that of spinqw_L12.json")
        rep.truth("Fig. 8 caption: the peak is declared as the grid argmax of the broadened line",
                  "the maximum of each broadened line on the" in fsp, "stated",
                  "the caption no longer says which peak its percentages refer to")

    # ---- (iii) the pooled sweep of Fig. S5, re-classified from the deposited sources --------
    try:
        P_ = []
        R_ = read_json(root, "data/thm1iii_violation/thm3_results.json")
        for r in R_:
            P_.append(dict(src="thm3", kind=r["kind"], rel=r["relL1"], err=r["LHS"], pub=r["RHS_pap"],
                           B=min(r["RHS_new"], r["RHS_triv"]), wS=r["wS"], frac=r["frac"], eta=r["eta"],
                           L=None, leak=bool(r["RHS_new"] < r["RHS_triv"])))
        L8_ = read_json(root, "data/thm1iii_violation/thm3_L8b.json")
        _l8src = read(root, "src/recovered/thm3_L8b.py")
        _m8 = re.search(r"^L,U,eta\s*=\s*([0-9.]+)\s*,\s*([0-9.]+)\s*,\s*([0-9.]+)\s*$", _l8src, re.M)
        _eta8 = float(_m8.group(3)) if _m8 else float("nan")
        for r in L8_:
            n2_ = r["LB"] / (1.0 - r["wS"])
            triv_ = n2_ * (1.0 + r["wS"])
            P_.append(dict(src="L8b", kind="topTE", rel=r["relL1"], err=r["LHS"], pub=r["RHS_pap"],
                           B=min(r["RHS_new"], triv_), wS=r["wS"], frac=r["FR"], eta=_eta8, L=8,
                           leak=bool(r["RHS_new"] < triv_)))
        nsrc_ = {"thm3": len(R_), "L8b": len(L8_)}
        shots_ = None
        for f_ in ("cert_stress.json", "cert_akw.json", "cert_teqsci.json"):
            d_ = read_json(root, "data/" + f_)
            nsrc_[f_] = len(d_["rows"])
            if f_ == "cert_teqsci.json":
                shots_ = d_.get("shots")
            for r in d_["rows"]:
                P_.append(dict(src=f_, kind="", rel=r["relL1_true"], err=r["L1_true"], pub=2.0 * r["lower"],
                               B=r["upper"], wS=r["w_S"], frac=r["frac"], eta=r["eta"], L=r.get("L"),
                               leak=bool(r["bound_nontrivial"])))
        cak_ = read_json(root, "data/cert_akw.json")
    except Exception as exc:                                    # noqa: BLE001
        rep.missing("Fig. S5 pool sources", "cannot read: %s" % exc)
        P_ = None
    if P_:
        NOISE_ = 1e-10
        live_ = [r for r in P_ if r["rel"] > NOISE_]
        inf_ = [r for r in live_ if r["pub"] <= 0.0]
        fin_ = [r for r in live_ if r["pub"] > 0.0]
        for r in fin_:
            r["ratio"] = r["err"] / r["pub"]
        for r in live_:
            r["ratioB"] = r["err"] / r["B"] if r["B"] > 0 else float("inf")
        n_pool, n_deg, n_live = len(P_), len(P_) - len(live_), len(live_)
        n_viol = len(inf_) + sum(1 for r in fin_ if r["ratio"] > 1.0)
        hi_ = [r for r in live_ if r["wS"] >= 0.99]
        hi_fin = [r for r in hi_ if r["pub"] > 0.0]
        hi_inf = [r for r in hi_ if r["pub"] <= 0.0]
        n_hiv = len(hi_inf) + sum(1 for r in hi_fin if r["ratio"] > 1.0)
        mild_ = min(r["ratio"] for r in hi_fin)
        fmax_ = max(r["ratio"] for r in fin_)
        lo_ = [r for r in live_ if r["wS"] < 0.99]
        lo_surv = sum(1 for r in lo_ if r["pub"] > 0.0 and r["ratio"] <= 1.0)
        n_violB = sum(1 for r in live_ if r["ratioB"] > 1.0)
        Bmax_ = max(r["ratioB"] for r in live_)
        Bmin_ = min(r["ratioB"] for r in live_)
        sat_ = [r for r in live_ if abs(r["ratioB"] - 1.0) < 1e-9]
        leak_n = sum(1 for r in live_ if r["leak"])
        n_fs = sum(1 for r in live_ if r["src"] == "cert_teqsci.json" or (r["src"] == "thm3" and r["kind"] == "teqsci"))
        f5_ = [r for r in fin_ if r["src"] == "cert_akw.json" and r["L"] == 8]
        etas_ = sorted(set(round(r["eta"], 4) for r in live_))
        print("   Fig. S5 pool: %d pooled, %d exact, %d live, %d zero-bound, %d violations, %d/%d at w_S>=0.99 "
              "(mildest %.4f), %d/%d survive below, ThmB %d violations, ratios %.3g..%.10f"
              % (n_pool, n_deg, n_live, len(inf_), n_viol, n_hiv, len(hi_), mild_, lo_surv, len(lo_), n_violB,
                 Bmin_, Bmax_))
        rep.eq_int("Fig. S5 pool: five sources present", len(nsrc_), 5, "a source of the pool is missing")
        rep.eq_int("Fig. S5 pool: Theorem 1 (the deposited bound) violated on none of the live rows", n_violB, 0,
                   "the deposited certificate is violated on the pooled sweep")
        rep.truth("Fig. S5 pool: every saturated row has w_S = 0 (A_S == 0, the trivial branch as an identity)",
                  all(r["wS"] == 0.0 for r in sat_) and len(sat_) > 0, "%d rows" % len(sat_),
                  "a row attains the bound with w_S > 0")
        rep.eq_int("Fig. S5 pool: no zero-bound row below w_S = 0.99", sum(1 for r in inf_ if r["wS"] < 0.99), 0,
                   "the caption's split of the infinite-ratio rows changed")
        rep.eq_int("Fig. S5 pool: every high-weight live row is a violation", n_hiv, len(hi_),
                   "the weight-only bound survives on a subspace with w_S >= 0.99")

        def _dat(name):
            out_ = []
            for ln in read(root, "paper/figs/" + name).splitlines()[1:]:
                if ln.strip():
                    a_, b_ = ln.split()
                    out_.append((a_, b_))
            return sorted(out_)

        def _fmt(pairs):
            return sorted(("%.6g" % x, "%.6g" % y) for x, y in pairs)

        for dat_, want_ in (("n3_pub_frac_hi.dat", [(r["frac"], r["ratio"]) for r in fin_ if r["wS"] >= 0.99]),
                            ("n3_pub_frac_lo.dat", [(r["frac"], r["ratio"]) for r in fin_ if r["wS"] < 0.99]),
                            ("n3_pub_eta_hi.dat", [(r["eta"], r["ratio"]) for r in fin_ if r["wS"] >= 0.99]),
                            ("n3_pub_eta_lo.dat", [(r["eta"], r["ratio"]) for r in fin_ if r["wS"] < 0.99]),
                            ("n3_thmB_frac.dat", [(r["frac"], r["ratioB"]) for r in live_]),
                            ("n3_thmB_eta.dat", [(r["eta"], r["ratioB"]) for r in live_]),
                            ("n3_inf_frac.dat", [(r["frac"], 1.0) for r in inf_]),
                            ("n3_inf_eta.dat", [(r["eta"], 1.0) for r in inf_]),
                            ("n3_fig5_frac.dat", [(r["frac"], r["ratio"]) for r in f5_]),
                            ("n3_fig5_eta.dat", [(r["eta"], r["ratio"]) for r in f5_])):
            got_d = _dat(dat_)
            rep.eq_int("Fig. S5 table %s: row count from the certificate rows" % dat_, len(got_d), len(want_),
                       "the plotted table does not hold one point per classified row", n=1)
            rep.truth("Fig. S5 table %s: every plotted point is a classified certificate row (6 s.f.)" % dat_,
                      got_d == _fmt(want_), "%d points" % len(got_d),
                      "the plotted ratios are not the error-to-bound ratios of the deposited rows", n=len(want_))
        # the whisker summary of the 64 rows of the Fig. 3 setting
        try:
            ws_ = read_json(root, "data/thm1iii_violation/n3_fig5_summary.json")
        except Exception as exc:                                # noqa: BLE001
            rep.missing("data/thm1iii_violation/n3_fig5_summary.json", "cannot read: %s" % exc)
            ws_ = None
        med_ = {}
        if ws_:
            f5d = [(float(a), float(b)) for a, b in _dat("n3_fig5_frac.dat")]
            rep.eq_int("Fig. S5 whiskers: 64 rows of the Fig. 3 setting", int(ws_["n_rows"]), len(f5d),
                       "the whisker summary does not describe the 64 plotted rows")
            for g_ in ws_["frac_groups"]:
                ys = sorted(y for x, y in f5d if abs(x - g_["x"]) < 1e-9)
                n_ = len(ys)
                m_ = 0.5 * (ys[n_ // 2 - 1] + ys[n_ // 2]) if n_ % 2 == 0 else ys[n_ // 2]
                med_[g_["x"]] = m_
                rep.eq("Fig. S5 whiskers at fraction %.2f: min" % g_["x"], float(g_["min"]), ys[0], 1e-9,
                       "the whisker minimum is not that of the plotted rows")
                rep.eq("Fig. S5 whiskers at fraction %.2f: max" % g_["x"], float(g_["max"]), ys[-1], 1e-9,
                       "the whisker maximum is not that of the plotted rows")
                rep.eq("Fig. S5 whiskers at fraction %.2f: median" % g_["x"], float(g_["median"]), m_, 1e-9,
                       "the whisker median is not that of the plotted rows")
                rep.eq_int("Fig. S5 whiskers at fraction %.2f: 16 rows" % g_["x"], int(g_["n"]), n_,
                           "a fraction group does not hold 8 momenta x 2 branches")
            rep.eq("Fig. S5 whiskers: overall min", float(ws_["overall_min"]), min(y for _, y in f5d), 1e-9,
                   "the overall minimum is not that of the plotted rows")
            rep.eq("Fig. S5 whiskers: overall max", float(ws_["overall_max"]), max(y for _, y in f5d), 1e-9,
                   "the overall maximum is not that of the plotted rows")
        # the 3.1e-8 control: cert_akw.json's own L=8, FR=0.85 per-momentum errors on the
        # figure window against the per-k errors of the deposited sampled_akw_L8.json
        try:
            mine_ = {kk: vv["0.85"]["relL1_figure_window"] for kk, vv in cak_["per_L"]["8"]["per_k"].items()}
            dev_ = max(abs(mine_[kk] - sk["per_k"][kk]["relL1"]) for kk in sk["per_k"])
        except Exception as exc:                                # noqa: BLE001
            rep.missing("cert_akw.json L=8 control", "cannot recompute: %s" % exc)
            dev_ = None
        if dev_ is not None:
            rep.eq("cert_akw.json: provenance control == max per-k deviation from sampled_akw_L8.json",
                   float(cak_["provenance"]["controls_reproduced"]["max_abs_dev_per_k"]), dev_, 1e-15,
                   "the stored control is not the deviation of the stored rows")
            rep.truth("cert_akw.json: the Fig. 3 setting reproduces the deposited sampled_akw_L8.json below 5e-5",
                      dev_ < 5e-5, "%.2e" % dev_, "the highlighted rows of Fig. S5 are not the Fig. 3 reconstruction")
        # ---- the printed numbers: caption, Sec. III D, Sec. S2, Sec. S9 ----
        capn = _tex("paper/figs/caption_thm1iii_violation.tex")
        DN = "the Fig. S5 caption prints a count that is not the re-classified pool's"
        _pnum(r"We pooled \$(\d+)\$ subspaces", "caption", capn, "Fig. S5 caption: pooled", n_pool, 0, DN)
        m_ = re.search(r"\(\$(\d+)\+(\d+)\$ from the original sweep; \$(\d+)\+(\d+)\+(\d+)\$", capn)
        if not m_:
            rep.missing("Fig. S5 caption: source counts", "the '239+4 ... 162+112+104' parenthesis is gone")
        else:
            for g_, k_ in zip(m_.groups(), ("thm3", "L8b", "cert_stress.json", "cert_akw.json", "cert_teqsci.json")):
                rep.eq_int("Fig. S5 caption: rows from %s" % k_, int(g_), nsrc_[k_], DN)
        _pnum(r"\$(\d+)\$ of them reproduce the target", "caption", capn, "Fig. S5 caption: exact rows", n_deg, 0, DN)
        _pnum(r"leaving \$(\d+)\$ live subspaces", "caption", capn, "Fig. S5 caption: live rows", n_live, 0, DN)
        _pnum(r"\$(\d+)\$ are finite-shot realizations", "caption", capn, "Fig. S5 caption: finite-shot rows", n_fs, 0, DN)
        _pnum(r"finite-shot realizations at \$(\d+)\$ shots per snapshot", "caption", capn,
              "Fig. S5 caption: shots per snapshot", float(shots_ or -1), 0, DN)
        _pnum(r"violated on \$(\d+)\$ of the \$\d+\$", "caption", capn, "Fig. S5 caption: violations", n_viol, 0, DN)
        _pnum(r"violated on \$\d+\$ of the \$(\d+)\$", "caption", capn, "Fig. S5 caption: live rows (again)", n_live, 0, DN)
        _pnum(r"On \$(\d+)\$ of them it is identically zero", "caption", capn, "Fig. S5 caption: zero-bound rows", len(inf_), 0, DN)
        _pnum(r"Among the \$(\d+)\$ finite-ratio points", "caption", capn, "Fig. S5 caption: finite high-weight rows", len(hi_fin), 0, DN)
        _pnum(r"the mildest by a factor \$([\d.]+)\$; with the", "caption", capn, "Fig. S5 caption: mildest factor",
              round(mild_, 2), 1e-12, DN)
        _pnum(r"with the \$(\d+)\$ infinite ones that is", "caption", capn, "Fig. S5 caption: infinite high-weight rows",
              len(hi_inf), 0, DN)
        _pnum(r"infinite ones that is \$(\d+)\$ of \$\d+\$", "caption", capn, "Fig. S5 caption: high-weight violations", n_hiv, 0, DN)
        _pnum(r"infinite ones that is \$\d+\$ of \$(\d+)\$", "caption", capn, "Fig. S5 caption: high-weight rows", len(hi_), 0, DN)
        m_ = re.search(r"the finite factors reach \$([\d.]+)\\times10\^\{(\d+)\}\$", capn)
        if not m_:
            rep.missing("Fig. S5 caption: largest finite factor", "'the finite factors reach' is gone")
        else:
            rep.eq("Fig. S5 caption: largest finite factor", float(m_.group(1)) * 10 ** int(m_.group(2)), fmax_,
                   0.005 * 10 ** int(m_.group(2)), DN)
        _pnum(r"the bound survives on \$(\d+)\$ of \$\d+\$ points", "caption", capn, "Fig. S5 caption: survivors below 0.99",
              lo_surv, 0, DN)
        _pnum(r"the bound survives on \$\d+\$ of \$(\d+)\$ points", "caption", capn, "Fig. S5 caption: rows below 0.99",
              len(lo_), 0, DN)
        if dev_ is not None:
            _pnum(r"to \$([\d.]+)\\times10\^\{-8\}\$ per momentum", "caption", capn,
                  "Fig. S5 caption: control of the plotted Fig. 3 rows (x1e-8)", round(dev_ * 1e8, 1), 1e-9,
                  "the caption quotes a control that is not the deviation of the rows Fig. S5 plots")
        if f5_:
            f5r = [r["ratio"] for r in f5_]
            f85 = [r["ratio"] for r in f5_ if abs(r["frac"] - 0.85) < 1e-9]
            m_ = re.search(r"violation factor is \$([\d.]+)\$--\$([\d.]+)\$, of which \$([\d.]+)\$--\$([\d.]+)\$ at the fraction \$0\.85\$", capn)
            if not m_:
                rep.missing("Fig. S5 caption: Fig. 3-setting range", "'violation factor is A--B, of which C--D at the fraction 0.85' is gone")
            else:
                for g_, w_, nm_ in zip(m_.groups(), (min(f5r), max(f5r), min(f85), max(f85)),
                                       ("min", "max", "min at 0.85", "max at 0.85")):
                    rep.eq("Fig. S5 caption: Fig. 3-setting violation factor %s" % nm_, float(g_), round(w_, 1), 1e-12, DN)
            if med_:
                m_ = re.search(r"\(medians \$([\d.]+), ([\d.]+), ([\d.]+), ([\d.]+)\$ at fractions \$0\.50, 0\.70, 0\.85, 0\.95\$\)", capn)
                if not m_:
                    rep.missing("Fig. S5 caption: medians", "the '(medians ... at fractions 0.50, 0.70, 0.85, 0.95)' parenthesis is gone")
                else:
                    for g_, x_ in zip(m_.groups(), (0.5, 0.7, 0.85, 0.95)):
                        rep.eq("Fig. S5 caption: median at fraction %.2f" % x_, float(g_), round(med_[x_], 1), 1e-12, DN)
        m_ = re.search(r"is violated on none of the \$(\d+)\$; its ratios span \$([\d.]+)\\times10\^\{-5\}\$ to \$([\d.]+)\$", capn)
        if not m_:
            rep.missing("Fig. S5 caption: Theorem 1 ratio range", "'violated on none of the 606; its ratios span' is gone")
        else:
            rep.eq_int("Fig. S5 caption: Theorem 1 tested on every live row", int(m_.group(1)), n_live, DN)
            rep.eq("Fig. S5 caption: smallest Theorem 1 ratio (x1e-5)", float(m_.group(2)), round(Bmin_ * 1e5, 2), 1e-12, DN)
            rep.eq("Fig. S5 caption: largest Theorem 1 ratio", float(m_.group(3)), round(Bmax_, 10), 1e-12, DN)
        rep.truth("Fig. S5 caption: the maximum is reached on four subspaces",
                  "reached on four subspaces" in capn and len(sat_) == 4, "%d" % len(sat_), DN)
        _pnum(r"is the smaller of the two on only \$(\d+)\$ of the \$\d+\$ rows", "caption", capn,
              "Fig. S5 caption: rows on which the leakage branch is the smaller", leak_n, 0, DN)
        m_ = re.search(r"Hubbard-chain rows exist only at \$\\eta/t=([\d.,\\ ]+)\$ and the stress[^$]*\$([\d.,\\ ]+)\$", capn)
        if m_:
            e_ = sorted(float(x) for x in (m_.group(1) + "," + m_.group(2)).replace("\\,", "").split(",") if x.strip())
            rep.truth("Fig. S5 caption: the eta values named are the eta values of the pool", e_ == etas_,
                      "%s" % etas_, "the caption's list of resolutions is not the pool's")
        # Sec. III D
        s31 = _tex("paper/sec_3_1.tex")
        DN3 = "Sec. III D prints a count that is not the re-classified pool's"
        _pnum(r"fails on \$(\d+)\$, and on all \$\d+\$ with \$w_S\\ge0\.99\$", "sec_3_1.tex", s31, "Sec. III D: violations", n_viol, 0, DN3)
        _pnum(r"fails on \$\d+\$, and on all \$(\d+)\$ with \$w_S\\ge0\.99\$", "sec_3_1.tex", s31, "Sec. III D: high-weight rows", len(hi_), 0, DN3)
        _pnum(r"with \$w_S\\ge0\.99\$, the mildest by a factor \$([\d.]+)\$", "sec_3_1.tex", s31, "Sec. III D: mildest factor",
              round(mild_, 2), 1e-12, DN3)
        if med_:
            _pnum(r"median violation factor rises from \$([\d.]+)\$ to \$[\d.]+\$ between sector fractions \$0\.50\$ and \$0\.95\$",
                  "sec_3_1.tex", s31, "Sec. III D: median at 0.50", round(med_[0.5], 1), 1e-12, DN3)
            _pnum(r"median violation factor rises from \$[\d.]+\$ to \$([\d.]+)\$ between sector fractions \$0\.50\$ and \$0\.95\$",
                  "sec_3_1.tex", s31, "Sec. III D: median at 0.95", round(med_[0.95], 1), 1e-12, DN3)
        rep.truth("Sec. III D prints one value per count (no 'on the deposited certificate files' aside)",
                  "on the deposited certificate files" not in s31 and "$468$" not in s31 and "$2.10$" not in s31,
                  "one set", "Sec. III D still carries the earlier-run values beside the regenerated ones")
        # Sec. S9 (sec_3_app)
        s3a = _tex("paper/sec_3_app.tex")
        DN9 = "Sec. S9 prints a count that is not the re-classified pool's"
        m_ = re.search(r"contribute \$(\d+)\$ subspaces: \$(\d+)\$ from the original sweep, and \$(\d+)\+(\d+)\+(\d+)\$", s3a)
        if not m_:
            rep.missing("Sec. S9: source counts", "'contribute 621 subspaces: 243 from the original sweep, and 162+112+104' is gone")
        else:
            rep.eq_int("Sec. S9: pooled", int(m_.group(1)), n_pool, DN9)
            rep.eq_int("Sec. S9: original sweep rows", int(m_.group(2)), nsrc_["thm3"] + nsrc_["L8b"], DN9)
            for g_, k_ in zip(m_.groups()[2:], ("cert_stress.json", "cert_akw.json", "cert_teqsci.json")):
                rep.eq_int("Sec. S9: rows from %s" % k_, int(g_), nsrc_[k_], DN9)
        _pnum(r"leaving \$(\d+)\$ live", "sec_3_app.tex", s3a, "Sec. S9: live rows", n_live, 0, DN9)
        _pnum(r"of which \$(\d+)\$ are finite-shot realizations at \$\d+\$ shots", "sec_3_app.tex", s3a, "Sec. S9: finite-shot rows", n_fs, 0, DN9)
        _pnum(r"the weight-only bound fails on \$(\d+)\$; among the", "sec_3_app.tex", s3a, "Sec. S9: violations", n_viol, 0, DN9)
        _pnum(r"among the \$(\d+)\$ with \$w_S\\ge0\.99\$", "sec_3_app.tex", s3a, "Sec. S9: high-weight rows", len(hi_), 0, DN9)
        _pnum(r"it fails on all \$(\d+)\$, the mildest by a factor", "sec_3_app.tex", s3a, "Sec. S9: high-weight violations", n_hiv, 0, DN9)
        _pnum(r"it fails on all \$\d+\$, the mildest by a factor \$([\d.]+)\$, with \$\d+\$ of them", "sec_3_app.tex", s3a,
              "Sec. S9: mildest factor", round(mild_, 2), 1e-12, DN9)
        _pnum(r"the mildest by a factor \$[\d.]+\$, with \$(\d+)\$ of them carrying \$w_S=1\$", "sec_3_app.tex", s3a,
              "Sec. S9: zero-bound rows", len(inf_), 0, DN9)
        _pnum(r"it survives on \$(\d+)\$ of \$\d+\$, in every case", "sec_3_app.tex", s3a, "Sec. S9: survivors below 0.99", lo_surv, 0, DN9)
        _pnum(r"it survives on \$\d+\$ of \$(\d+)\$, in every case", "sec_3_app.tex", s3a, "Sec. S9: rows below 0.99", len(lo_), 0, DN9)
        _pnum(r"the largest ratio of error to bound is \$([\d.]+)\$", "sec_3_app.tex", s3a, "Sec. S9: largest Theorem 1 ratio",
              round(Bmax_, 10), 1e-12, DN9)
        _pnum(r"is the smaller of the two on only \$(\d+)\$ of the \$\d+\$", "sec_3_app.tex", s3a,
              "Sec. S9: rows on which the leakage branch is the smaller", leak_n, 0, DN9)
        rep.truth("Sec. S9 no longer explains the plotted pool as an earlier run",
                  "earlier run" not in s3a, "absent", "Sec. S9 still says the plotted pool was read from an earlier run")
        # Sec. S2 (sm_leakage)
        sml = _tex("paper/sm_leakage.tex")
        DN2 = "Sec. S2 prints a count that is not the re-classified pool's"
        _pnum(r"contributes \$(\d+)\$ subspaces from five Hamiltonian families", "sm_leakage.tex", sml, "Sec. S2: pooled", n_pool, 0, DN2)
        _pnum(r"implementations; \$(\d+)\$ reproduce their target", "sm_leakage.tex", sml, "Sec. S2: exact rows", n_deg, 0, DN2)
        _pnum(r"On the \$(\d+)\$ that remain the weight-only bound fails on \$\d+\$", "sm_leakage.tex", sml, "Sec. S2: live rows", n_live, 0, DN2)
        _pnum(r"On the \$\d+\$ that remain the weight-only bound fails on \$(\d+)\$", "sm_leakage.tex", sml, "Sec. S2: violations", n_viol, 0, DN2)
        _pnum(r"it fails on \$(\d+)\$ of \$\d+\$, the mildest by a factor", "sm_leakage.tex", sml, "Sec. S2: high-weight violations", n_hiv, 0, DN2)
        _pnum(r"it fails on \$\d+\$ of \$(\d+)\$, the mildest by a factor", "sm_leakage.tex", sml, "Sec. S2: high-weight rows", len(hi_), 0, DN2)
        _pnum(r"of \$\d+\$, the mildest by a factor \$([\d.]+)\$ \(Sec", "sm_leakage.tex", sml, "Sec. S2: mildest factor",
              round(mild_, 2), 1e-12, DN2)
        if med_ and f5_:
            m_ = re.search(r"is \$([\d.]+)\$, \$([\d.]+)\$, \$([\d.]+)\$, \$([\d.]+)\$, the \$(\d+)\$ rows spanning \$([\d.]+)\$ to \$([\d.]+)\$ "
                           r"overall and \$([\d.]+)\$ to \$([\d.]+)\$ at the fraction \$0\.85\$", sml)
            if not m_:
                rep.missing("Sec. S2: the medians and ranges of the Fig. 3 setting", "the sentence is gone or reworded")
            else:
                for g_, x_ in zip(m_.groups()[:4], (0.5, 0.7, 0.85, 0.95)):
                    rep.eq("Sec. S2: median at fraction %.2f" % x_, float(g_), round(med_[x_], 1), 1e-12, DN2)
                rep.eq_int("Sec. S2: rows of the Fig. 3 setting", int(m_.group(5)), len(f5_), DN2)
                f5r = [r["ratio"] for r in f5_]
                f85 = [r["ratio"] for r in f5_ if abs(r["frac"] - 0.85) < 1e-9]
                for g_, w_, nm_ in zip(m_.groups()[5:], (min(f5r), max(f5r), min(f85), max(f85)),
                                       ("min", "max", "min at 0.85", "max at 0.85")):
                    rep.eq("Sec. S2: Fig. 3-setting violation factor %s" % nm_, float(g_), round(w_, 1), 1e-12, DN2)
            _pnum(r"rise\} monotonically, by a factor \$([\d.]+)\$", "sm_leakage.tex", sml, "Sec. S2: rise of the medians",
                  round(med_[0.95] / med_[0.5], 1), 1e-12, DN2)
        rep.truth("no manuscript sentence quotes the earlier-run counts 468 / 2.10 / 138",
                  all(x_ not in " ".join(read_body(root).split()) for x_ in ("on $468$", "$468$ of", "factor $2.10$", "on $138$ of")),
                  "absent", "an earlier-run value of the Fig. S5 pool is still printed")

    # ---------------------------------------------------------------- 9.7  fig_method
    mm = ctx["mm"]
    mtx = read(root, "paper/figs/fig_method_native_frag.tex")
    names, cols = table_cols(root, "paper/figs/aw_method.dat")
    rep.truth("aw_method.dat has the three columns the method figure reads",
              names == ["w", "Aex", "Asm"], "columns %s" % names,
              "the columns of the .dat the figure plots changed: %s" % names)
    for col, key in (("w", "grid"), ("Aex", "A_exact"), ("Asm", "A_final")):
        if col not in names:
            continue
        v = cols[names.index(col)]
        ref = mm[key]
        if len(v) != len(ref):
            rep.emit(FAIL, "aw_method.dat column %s" % col,
                     "the .dat the method figure actually plots has %d rows against %d "
                     "deposited values -- rows were added or deleted"
                     % (len(v), len(ref)), n=1)
            continue
        bad = sum(1 for a, b in zip(v, ref) if abs(a - b) > 5e-7 * max(1.0, abs(b)))
        rep.truth("aw_method.dat column %s reproduces from data/method_max.json" % col,
                  bad == 0, "%d values" % len(v),
                  "%d of %d values of the plotted .dat do not reproduce from "
                  "data/method_max.json" % (bad, len(v)), n=len(v))
    m = re.search(r"rel-\$L_1=([\d.]+)\$", mtx)
    if m:
        grid = np_.array(mm["grid"])
        rel = float(np_.trapz(np_.abs(np_.array(mm["A_final"]) - np_.array(mm["A_exact"])),
                              grid) / np_.trapz(np_.array(mm["A_exact"]), grid))
        rep.eq("fig_method (a) annotated rel-L1", float(m.group(1)), round(rel, 3), 5e-4,
               "the annotated rel-L1 is not the relative L1 distance between the deposited "
               "A_final and A_exact (%.6f)" % rel)
    # the geometric-fit annotation
    m = re.search(RBS + r"rho" + RBS + r"approx([\d.]+)", mtx)
    if m:
        rep.eq("fig_method (b) annotated geometric ratio rho", float(m.group(1)),
               round(mm["rho"], 1), 5e-2,
               "the annotated geometric ratio is not data/method_max.json['rho'] (%.6f)"
               % mm["rho"])

    # ---------------------------------------------------------------- 9.8  hardware hero
    hh = read(root, "paper/figs/fig_hardware_hero_frag.tex")
    hls = logical_lines(hh)
    lucj = read_json(root, "data/hw_lucj_n2_result.json")
    efci = lucj["e_fci"]
    sim = [l for l in hls if "simSlate,line width" in l]
    if not sim:
        rep.missing("fig_hardware_hero (b) noiseless-simulator trajectory", "series gone")
    else:
        pts_ = series_coords(sim[0])
        ref = [(e - efci) * 1000.0 for e in lucj["hist_sim"]]
        rep.eq_int("fig_hardware_hero (b) simulator point count", len(pts_), len(ref),
                   "the simulator trajectory changed length; it is hist_sim of "
                   "data/hw_lucj_n2_result.json")
        for (x, y), r in zip(pts_, ref):
            rep.eq("fig_hardware_hero (b) simulator dE at iter %d" % int(x), y, r, 5e-3,
                   "the plotted noiseless-simulator error is not (hist_sim - e_fci) of "
                   "data/hw_lucj_n2_result.json")
    hw = [l for l in hls if "hwWarm,line width" in l]
    hw = hw + [l for l in hls if "coordinates {(0,31.079)" in l]
    hwl = [l for l in hls if "+-" in l and "coordinates" in l]
    if hwl:
        pts_ = series_coords(hwl[0])
        sig = errbars(hwl[0])
        rep.eq("fig_hardware_hero (b) final hardware dE", pts_[-1][1],
               lucj["dE_mean_mHa"], 5e-4,
               "the final hardware point is not dE_mean_mHa of "
               "data/hw_lucj_n2_result.json (%.6f)" % lucj["dE_mean_mHa"])
        if sig:
            rep.eq("fig_hardware_hero (b) final hardware sigma", sig[-1],
                   lucj["dE_std_mHa"], 5e-4,
                   "the final hardware error bar is not dE_std_mHa of "
                   "data/hw_lucj_n2_result.json (%.6f)" % lucj["dE_std_mHa"])
        # the band must be the mean +- sigma it advertises
        hi = [l for l in hls if "name path=hwhi" in l]
        lo = [l for l in hls if "name path=hwlo" in l]
        if hi and lo:
            H = [y for _, y in series_coords(hi[0])]
            Lo = [y for _, y in series_coords(lo[0])]
            bad = [i for i in range(min(len(H), len(Lo), len(pts_), len(sig)))
                   if abs((H[i] - Lo[i]) / 2.0 - sig[i]) > 1e-3
                   or abs((H[i] + Lo[i]) / 2.0 - pts_[i][1]) > 1e-3]
            rep.truth("fig_hardware_hero (b) the shaded band is exactly mean +- sigma",
                      not bad, "%d band points consistent" % len(H),
                      "the shaded hardware band is not the plotted mean plus/minus the "
                      "plotted sigma at indices %s -- a band drawn independently of its "
                      "own series can be drawn anywhere" % bad, n=max(len(H), 1))
    # the .dat panel (a) actually plots
    names, cols = table_cols(root, "paper/figs/heron_hot.dat")
    her = ctx["her"]
    if names == ["w", "A"]:
        g, A = cols[0], cols[1]
        rep.eq_int("heron_hot.dat row count", len(g), len(her["grid"]),
                   "the .dat the hardware figure plots no longer has one row per "
                   "deposited grid point")
        n = min(len(g), len(her["grid"]))
        bad = sum(1 for i in range(n)
                  if abs(g[i] - her["grid"][i]) > 5e-5
                  or abs(A[i] - her["A_hw"][i]) > 5e-6)
        rep.truth("heron_hot.dat reproduces from data/heron_spectral.json",
                  bad == 0, "%d rows" % n,
                  "%d of %d rows of the plotted hardware .dat do not reproduce from the "
                  "deposited A_hw" % (bad, n), n=n)
    else:
        rep.missing("heron_hot.dat", "unexpected columns %s" % names)

    # ---------------------------------------------------------------- 9.9  noise recovery
    nr = read(root, "paper/figs/fig_noise_recovery_native.tex")
    nls = logical_lines(nr)
    sw = read_json(root, "data/molecular_noise_sweep.json")
    bymol = dict((m_["mol"], m_) for m_ in sw["mols"])
    tags = {"nrH6": "H6", "nrNH3": "NH3", "nrN2": "N2", "nrCO": "CO", "nrC2": "C2"}
    for tag, mol in sorted(tags.items()):
        if mol not in bymol:
            rep.missing("fig_noise_recovery %s" % mol,
                        "molecule not in data/molecular_noise_sweep.json")
            continue
        rows = bymol[mol]["rows"]
        solid = [l for l in nls if ("addplot[" + tag + ",line width=1.4pt") in l]
        dash = [l for l in nls if ("addplot[" + tag + ",line width=0.9pt") in l]
        lo = [l for l in nls if ("name path=" + mol + "lo") in l]
        hi = [l for l in nls if ("name path=" + mol + "hi") in l]
        for lines_, key, label in ((solid, "score", "S-CoRe"), (dash, "naive", "naive")):
            if not lines_:
                rep.missing("fig_noise_recovery %s %s series" % (mol, label),
                            "series gone from the figure source")
                continue
            pts_ = series_coords(lines_[0])
            rep.eq_int("fig_noise_recovery %s %s point count" % (mol, label),
                       len(pts_), len(rows),
                       "the plotted series has a different number of points from the "
                       "deposited sweep")
            for (x, y), r in zip(pts_, rows):
                rep.eq("fig_noise_recovery %s %s at eps=%g" % (mol, label, x), y,
                       r[key], 5e-5,
                       "the plotted error disagrees with data/molecular_noise_sweep.json")
        if lo and hi:
            Lo = [y for _, y in series_coords(lo[0])]
            Hi = [y for _, y in series_coords(hi[0])]
            # NOT CHECKABLE from the deposit: the band is a per-seed spread and
            # data/molecular_noise_sweep.json stores only the mean and its standard
            # deviation (the band is visibly asymmetric, so it is not mean +- std).
            # What is an invariant, and what a fabricated band would break, is that it
            # brackets its own curve and is contained in the deposited +-1 sigma scale.
            bad = [i for i in range(min(len(Lo), len(Hi), len(rows)))
                   if not (Lo[i] <= rows[i]["score"] <= Hi[i])
                   or Hi[i] - Lo[i] > 6.0 * max(rows[i]["score_std"], 1e-12)]
            rep.truth("fig_noise_recovery %s band brackets the deposited score" % mol,
                      not bad, "%d band points bracket their curve" % len(Lo),
                      "the shaded band for %s does not contain the deposited score, or is "
                      "wider than six deposited standard deviations, at indices %s -- a "
                      "band drawn independently of its own data can be drawn anywhere"
                      % (mol, bad), n=max(len(Lo), 1))

    # ---------------------------------------------------------------- 9.9b printed noise numbers
    # 2026-09-30 (final numeric review NR-2).  The captions of Figs. noise and noiserec, Sec. S8
    # and Sec. VIII print values of data/noise_spectral.json, data/molecular_noise_sweep.json and
    # data/pole_structure.json (generators src/noise_spectral.py, src/molecular_noise_sweep.py,
    # src/pole_structure.py) that no check read: a one-digit change in any of them left every
    # guardian green.  Each printed value is now its rounding of the deposited one, and the two
    # qualitative claims they carry are recomputed: 'beyond the seed bands' (mean+s.d. at the
    # quoted eps below mean-s.d. at eps=0) and 'near 0.05 for 4% <= eps <= 16%'.
    try:
        ns_ = read_json(root, "data/noise_spectral.json")
    except Exception as exc:                                    # noqa: BLE001
        rep.missing("data/noise_spectral.json", "cannot read: %s" % exc)
        ns_ = None
    fn_ = _tex("paper/carried/fig_noise.tex")
    fr_ = _tex("paper/carried/fig_noiserec.tex")
    ss_ = _tex("paper/sm_scope.tex")
    s8_ = _tex("paper/sec_8.tex")
    DN = "the printed noise value is not its rounding of the deposited sweep"
    if ns_:
        eps_ = [float(e) for e in ns_["eps"]]
        i0_, i13_ = eps_.index(0.0), eps_.index(0.13)
        rep.eq("noise_spectral.json: the two arms coincide at eps=0", ns_["score"][i0_], ns_["naive"][i0_],
               0.0, "recovery and naive post-selection differ on noiseless bitstrings")
        for rel_, txt_, pat_ in (
                ("carried/fig_noise.tex", fn_, r"the two arms coincide at \$([\d.]+)\$"),
                ("sm_scope.tex", ss_, r"from \$([\d.]+)\$ at \$\\varepsilon=0\$ to \$[\d.]+\$ at \$13\\%\$")):
            _pnum(pat_, rel_, txt_, "%s: spectral rel-L1 at eps=0" % rel_,
                  round(ns_["score"][i0_], 3), 1e-12, DN)
        _pnum(r"at \$\\varepsilon=0\$ to \$([\d.]+)\$ at \$13\\%\$", "sm_scope.tex", ss_,
              "sm_scope.tex: spectral rel-L1 at eps=13%", round(ns_["score"][i13_], 3), 1e-12, DN)
        rep.truth("noise_spectral.json: eps=13% lies below eps=0 beyond the seed bands",
                  ns_["score"][i13_] + ns_["score_std"][i13_] < ns_["score"][i0_] - ns_["score_std"][i0_],
                  "beyond", "Sec. S8 reports the eps=13% value 'beyond the seed bands' and it is not")
        m_ = re.search(r"keeps the error near \$([\d.]+)\$ for \$(\d+)\\%\\le\\varepsilon\\le(\d+)\\%\$", fn_)
        if not m_:
            rep.missing("carried/fig_noise.tex: 'near 0.05 for 4% <= eps <= 16%'", "the range is gone")
        else:
            c_, a_, b_ = float(m_.group(1)), float(m_.group(2)) / 100, float(m_.group(3)) / 100
            inr_ = [s for e, s in zip(eps_, ns_["score"]) if a_ - 1e-12 <= e <= b_ + 1e-12]
            rep.truth("carried/fig_noise.tex: every recovered score in [%g, %g] is within 0.01 of %g"
                      % (a_, b_, c_), len(inr_) >= 2 and all(abs(s - c_) <= 0.01 for s in inr_),
                      "%d scores" % len(inr_), "a recovered score in the quoted range is not near %g" % c_)
    try:
        ms_ = dict((m_["mol"], m_["rows"]) for m_ in read_json(root, "data/molecular_noise_sweep.json")["mols"])
    except Exception as exc:                                    # noqa: BLE001
        rep.missing("data/molecular_noise_sweep.json (9.9b)", "cannot read: %s" % exc)
        ms_ = {}
    for mol_, tex_, eps_q, nd_ in (("NH3", r"\\mathrm\{NH_3\}\$", 0.10, 2), ("CO", "CO", 0.20, 2),
                                  ("C2", r"\\mathrm\{C_2\}\$", 0.30, 3)):
        rows_m = ms_.get(mol_)
        if not rows_m:
            rep.missing("molecular_noise_sweep.json %s" % mol_, "molecule gone")
            continue
        r0_ = [r for r in rows_m if r["eps"] == 0.0][0]
        rq_ = [r for r in rows_m if abs(r["eps"] - eps_q) < 1e-12][0]
        rep.truth("molecular_noise_sweep.json: %s at eps=%g lies below eps=0 beyond the seed bands"
                  % (mol_, eps_q), rq_["score"] + rq_["score_std"] < r0_["score"] - r0_["score_std"],
                  "beyond", "the SM and the caption report %s 'beyond the seed bands' and it is not" % mol_)
        pat_ = (r"%s from \$([\d.]+)\$ to \$([\d.]+)\$~mHa at \$?(?:\\varepsilon=)?(\d+)\\%%" % tex_)
        m_ = re.search(pat_, ss_)
        if not m_:
            rep.missing("sm_scope.tex: %s noise pair" % mol_, "pattern %s gone" % pat_)
        else:
            rep.eq("sm_scope.tex: %s error at eps=0 (mHa)" % mol_, float(m_.group(1)),
                   round(r0_["score"], nd_), 1e-12, DN)
            rep.eq("sm_scope.tex: %s error at eps=%g (mHa)" % (mol_, eps_q), float(m_.group(2)),
                   round(rq_["score"], nd_), 1e-12, DN)
            rep.eq("sm_scope.tex: %s quoted eps (%%)" % mol_, float(m_.group(3)), round(100 * eps_q), 0.0, DN)
    if "NH3" in ms_:
        r0_ = [r for r in ms_["NH3"] if r["eps"] == 0.0][0]
        rq_ = [r for r in ms_["NH3"] if abs(r["eps"] - 0.10) < 1e-12][0]
        m_ = re.search(r"\\mathrm\{NH_3\}\$ \$([\d.]+)\\to([\d.]+)\$~mHa at \$\\varepsilon=(\d+)\\%\$", fr_)
        if not m_:
            rep.missing("carried/fig_noiserec.tex: NH3 pair", "'NH3 2.53 -> 1.24 mHa at eps=10%' is gone")
        else:
            rep.eq("carried/fig_noiserec.tex: NH3 at eps=0 (mHa)", float(m_.group(1)), round(r0_["score"], 2), 1e-12, DN)
            rep.eq("carried/fig_noiserec.tex: NH3 at eps=10% (mHa)", float(m_.group(2)), round(rq_["score"], 2), 1e-12, DN)
            rep.eq("carried/fig_noiserec.tex: NH3 eps (%)", float(m_.group(3)), 10.0, 0.0, DN)
    try:
        ps_ = read_json(root, "data/pole_structure.json")
        cg12_ = [r for r in read_json(root, "data/charge_gap.json")["results"] if r["L"] == 12][0]
    except Exception as exc:                                    # noqa: BLE001
        rep.missing("data/pole_structure.json / charge_gap.json L=12 (9.9b)", "cannot read: %s" % exc)
        ps_ = None
    if ps_:
        ch_ = [r for r in ps_["table"] if r["channel"] == "charge" and r["L"] == 12]
        wmin_ = min(r["w_retained_frac"] for r in ch_)
        _pnum(r"---\$([\d.]+)\\%\$ or more of the weight at each \$q\$---lies above it", "sec_8.tex", s8_,
              "Sec. VIII: retained charge weight at every q (%)", math.floor(1000 * wmin_) / 10, 1e-12,
              "Sec. VIII: the printed floor is not floor(100 x min w_retained_frac, 1 d.p.)")
        rep.eq_int("pole_structure.json: charge momenta at L=12", len(ch_), 6,
                   "the charge channel no longer carries the six nonzero momenta of L=12")
        rep.truth("Sec. VIII: every retained charge pole at L=12 lies above Delta(12)",
                  min(r["lowest_retained_pole"] for r in ch_) > cg12_["Delta"],
                  "%.3f > %.3f" % (min(r["lowest_retained_pole"] for r in ch_), cg12_["Delta"]),
                  "a retained charge pole lies below the charge gap Sec. VIII says it lies above")

    # ---------------------------------------------------------------- 9.9c restated certificate numbers
    # 2026-09-30 (final numeric review NR-3).  The Fig. 3 certificate (relative bound 0.77-1.75,
    # 350-900 times the true error), the momentum-block need (455-483 of 490 at L=8), the spread
    # of the crossing (5.9) and the 0.98 of the sector at L=8 are printed in the abstract and in
    # Secs. I, V, X and S3, and no check read any copy.  Every copy is now its rounding of the
    # deposited record, and each is required where the text carries it.
    DC = "a restated certificate number is not its rounding of the deposited record"
    try:
        ca_ = [r for r in read_json(root, "data/cert_akw.json")["rows"]
               if r["L"] == 8 and abs(r["FR"] - 0.85) < 1e-9]
        kb_ = read_json(root, "data/c3_frontier/adversarial/z_kblock_L8.json")["rows"]
        co_ = read_json(root, "data/c3_frontier/null/CALIB_OOS.json")
        zf_ = read(root, "data/c3_frontier/adversarial/z_fine_L8.txt")
    except Exception as exc:                                    # noqa: BLE001
        rep.missing("certificate records (9.9c)", "cannot read: %s" % exc)
        ca_ = None
    if ca_:
        rep.eq_int("cert_akw.json: L=8, FR=0.85 channels", len(ca_), 16,
                   "the Fig. 3 certificate no longer has its sixteen channels")
        rep.truth("cert_akw.json: every L=8, FR=0.85 channel is non-vacuous",
                  all(r["bound_nontrivial"] for r in ca_), "16 of 16",
                  "a Fig. 3 channel is vacuous, and the text says all sixteen are bounded")
        lo_u, hi_u = min(r["rel_upper"] for r in ca_), max(r["rel_upper"] for r in ca_)
        lo_s, hi_s = min(r["slack"] for r in ca_), max(r["slack"] for r in ca_)

        def _sf2(x):
            return float("%.2g" % x)
        kb_lo, kb_hi = min(r["first_nonvacuous"] for r in kb_), max(r["first_nonvacuous"] for r in kb_)
        kb_blk = set(r["block_size"] for r in kb_)
        spread_ = max(r["relL1"] for r in co_) / min(r["relL1"] for r in co_)
        fine_ = [(int(m_.group(1)), float(m_.group(2))) for m_ in
                 re.finditer(r"\|S\|=(\d+) FR=[\d.]+ eta=0\.18 \| .*?relR5=([\d.e+-]+)", zf_)]
        above_ = [s for s, v in fine_ if v >= 0.05]
        below_ = [s for s, v in fine_ if v < 0.05]
        frac_ = round(min(below_) / 3920.0, 2) if below_ else float("nan")
        A_ = ("sec_1.tex", "sec_5_body.tex", "sec_10.tex", "sm_frontier.tex")
        checks_ = (   # (pattern, deposited value at its printed rounding, label, files that must carry it)
            (r"\$([\d.]+)\$--\$[\d.]+\$ against the trivial \$2\$", round(lo_u, 2), "relative bound, low", A_),
            (r"\$[\d.]+\$--\$([\d.]+)\$ against the trivial \$2\$", round(hi_u, 2), "relative bound, high", A_),
            (r"\$(\d+)\$--\$\d+\$ times (?:its|the) true error", _sf2(lo_s), "slack, low (2 s.f.)",
             ("main.tex",) + A_),
            (r"\$\d+\$--\$(\d+)\$ times (?:its|the) true error", _sf2(hi_s), "slack, high (2 s.f.)",
             ("main.tex",) + A_),
            (r"\$(\d+)\$--\$\d+\$ of (?:its )?\$490\$", kb_lo, "momentum block, first non-vacuous low",
             ("sec_5_body.tex", "sec_10.tex", "sm_frontier.tex")),
            (r"\$\d+\$--\$(\d+)\$ of (?:its )?\$490\$", kb_hi, "momentum block, first non-vacuous high",
             ("sec_5_body.tex", "sec_10.tex", "sm_frontier.tex")),
            (r"a spread of \$([\d.]+)\$", round(spread_, 1), "calibration spread",
             ("sec_1.tex", "sec_5_body.tex", "sm_frontier.tex")),
            (r"(?:from|at least) \$([\d.]+)\$ of the sector at \$L=8\$", frac_,
             "fraction at which the certificate reaches 0.05 at L=8", A_),
        )
        for pat_, want_, lab_, need_files in checks_:
            for rel_ in ("main.tex", "sec_1.tex", "sec_5_body.tex", "sec_10.tex", "sm_frontier.tex"):
                txt_ = _tex("paper/" + rel_)
                if rel_ == "main.tex":
                    m0_ = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", txt_)
                    txt_ = m0_.group(1) if m0_ else ""
                hits_ = [float(m_.group(1)) for m_ in re.finditer(pat_, txt_)]
                if not hits_:
                    if rel_ in need_files:
                        rep.missing("%s: %s" % (rel_, lab_), "the restated value is gone (pattern %s)" % pat_)
                    continue
                for h_ in hits_:
                    rep.eq("%s: %s" % (rel_, lab_), h_, want_, 1e-12, DC)
        rep.truth("z_fine_L8.txt: the 0.05 crossing is bracketed at integer sizes",
                  bool(above_) and bool(below_) and max(above_) < min(below_),
                  "%d < |S| <= %d" % (max(above_ or [0]), min(below_ or [0])),
                  "the fine scan no longer brackets the crossing the text places at 0.98")
        rep.truth("z_kblock_L8.json: the momentum blocks hold 490 determinants", kb_blk == {490},
                  "490", "the block size printed as 490 is not the deposited one (%s)" % sorted(kb_blk))

    # ---------------------------------------------------------------- 9.10 amortized
    am = read(root, "paper/figs/fig_amort_native.tex")
    ar = read_json(root, "data/amortized_recovery.json")
    summ = ar["summary"]
    for tag, mean_key, std_key in (("amCls", "class_mean", "class_std"),
                                   ("amGen", "amort_mean", "amort_std")):
        ln = [l for l in logical_lines(am) if ("addplot[" + tag + ",line width") in l]
        if not ln:
            rep.missing("fig_amort %s series" % tag, "series gone from the figure source")
            continue
        pts_ = series_coords(ln[0])
        sig = errbars(ln[0])
        rep.eq_int("fig_amort %s point count" % tag, len(pts_), len(summ),
                   "the plotted series has a different number of fractions from "
                   "data/amortized_recovery.json['summary']")
        for i, (x, y) in enumerate(pts_):
            key = ("%g" % x)
            if key not in summ:
                rep.emit(FAIL, "fig_amort %s at frac=%s" % (tag, key),
                         "FABRICATED POINT: fraction %s has no entry in "
                         "data/amortized_recovery.json['summary']" % key, n=1)
                continue
            rep.eq("fig_amort %s mean at frac=%s" % (tag, key), y, summ[key][mean_key],
                   5e-4, "the plotted mean rel-L1 disagrees with the deposited summary")
            if i < len(sig):
                rep.eq("fig_amort %s sigma at frac=%s" % (tag, key), sig[i],
                       summ[key][std_key], 5e-4,
                       "the plotted error bar disagrees with the deposited summary")

    # ---------------------------------------------------------------- 9.10b panel (c)
    # The support counts of the scaling figure live in TWO places: the plotted series and
    # the provenance comment that tells a referee where they come from.  Only one of the
    # two was covered, so the comment could be reverted to the superseded 19184 / 824504
    # while the graph stayed right -- and the comment is what a referee reads.  Both are
    # checked here, and they must be integers: |S| is a count of determinants, and a
    # tolerance that accepts 19183.4 accepts a non-integer number of determinants.
    sc = read(root, "paper/figs/fig_scaling2_native.tex")
    sdj = read_json(root, "data/scaling_data.json")
    from fractions import Fraction
    due = []
    for p in sorted(sdj["points"], key=lambda q: q["L"]):
        v = Fraction(p["frac"]).limit_denominator(10 ** 15) * p["Np1_sector"] \
            if False else Fraction(p["frac"]) * p["Np1_sector"]
        due.append(int(round(float(v))))
    m = re.search(r"with frac from scaling_data\.json: ([\d, ]+)\.", sc)
    plotted = None
    for ln in logical_lines(sc):
        if "coordinates {(8,22)" in ln or ("(8,%d)" % due[0]) in ln and "824" in ln:
            plotted = [int(round(y)) for _, y in series_coords(ln)]
            break
    rep.truth("fig_scaling2 (c) plots the support counts derivable from scaling_data.json",
              plotted == due, "%s" % due,
              "the plotted absolute supports are %s; frac x Np1_sector gives %s. |S| is a "
              "count of determinants: it must be an integer and it must be THIS integer"
              % (plotted, due), n=max(len(due), 1))
    if plotted is not None:
        nonint = [y for _, y in series_coords(
            [ln for ln in logical_lines(sc) if ("(8,%d)" % due[0]) in ln][0])
            if abs(y - round(y)) > 0]
        rep.truth("fig_scaling2 (c) support counts are integers", not nonint,
                  "all integers",
                  "non-integer determinant counts plotted: %s" % nonint)
    if m is None:
        rep.missing("fig_scaling2 provenance comment",
                    "the comment that tells a referee how |S| is derived is gone")
    else:
        incomment = [int(x) for x in m.group(1).replace(" ", "").split(",") if x]
        rep.truth("fig_scaling2 provenance comment carries the same counts as the graph",
                  incomment == due, "%s" % incomment,
                  "the provenance COMMENT says %s while the graph plots %s. The comment is "
                  "what a referee reads to find out where the number comes from; half a "
                  "repair is not a repair." % (incomment, due), n=max(len(due), 1))

    # ---------------------------------------------------------------- 9.11 legality of q
    # The fabricated rows were removed from the .dat.  That is not enough: the previous
    # guardian compared the .dat with the JSON, so laundering the same two rows THROUGH
    # the JSON passed -- with the coverage counter going UP.  The grid is now checked
    # against the physics, not against a sibling artefact.
    for jsonrel, datrel, L, label in (
            ("data/sqw_L12.json", "paper/figs/sqw_edges.dat", 12, "S(q,w)"),
            ("data/spinqw_L12.json", "paper/figs/spinqw_edges.dat", 12, "S^zz(q,w)")):
        d = read_json(root, jsonrel)
        keys = sorted(d["S"].keys(), key=float) if "S" in d else []
        if not keys:
            keys = sorted([k for k in d.keys() if re.match(r"^[\d.]+$", str(k))], key=float)
        qs = [float(k) for k in keys]
        legal = [2.0 * n / L for n in range(1, L)]          # q/pi = 2n/L, n = 1..L-1
        illegal = [q for q in qs if min(abs(q - g) for g in legal) > 1e-4]
        rep.truth("%s: every momentum in %s is an accessible lattice momentum"
                  % (label, jsonrel), not illegal,
                  "%d momenta, all of the form q/pi = 2n/L with 1 <= n <= L-1" % len(qs),
                  "ILLEGAL MOMENTA in %s: %s. q=0 is inaccessible by construction "
                  "(sqw_lanczos.py skips it; sqw_field.py sets S(q=0,w>0)=0 by density "
                  "conservation) and q/pi=2 is its image. A row whose momentum cannot "
                  "exist is fabricated whether or not a sibling file also carries it."
                  % (jsonrel, illegal), n=max(len(qs), 1))
        rep.truth("%s: q=0 and q/pi=2 are absent from %s" % (label, jsonrel),
                  not any(abs(q) < 1e-6 or abs(q - 2.0) < 1e-6 for q in qs),
                  "neither momentum is present",
                  "q=0 or q/pi=2 is present in %s. These are exactly the two rows that "
                  "were fabricated in the .dat on 2026-08-17; putting them in the dataset "
                  "instead is the same defect laundered one level down." % jsonrel)
        names, cols = table_cols(root, datrel)
        rep.eq_int("%s edge table row count" % label, len(cols[0]), len(qs),
                   "the deposited edge table has a different number of rows from the "
                   "spectrum it is derived from -- a DELETED row is as much a defect as "
                   "an invented one, and used to be invisible")
        want = set(round(q, 4) for q in qs)
        got = set(round(q, 4) for q in cols[0])
        rep.truth("%s edge table carries exactly the computed momenta" % label,
                  want == got, "%d momenta match" % len(want),
                  "the momenta of %s are not the momenta of %s: missing %s, extra %s"
                  % (datrel, jsonrel, sorted(want - got), sorted(got - want)),
                  n=max(len(want), 1))

    # ---------------------------------------------------------------- 9.12 sum_rule block
    # Nine numbers were added to data/method_max.json and nothing watched them; a full
    # falsification of the block passed.
    if "sum_rule" in mm:
        sr = mm["sum_rule"]
        grid = np_.array(mm["grid"])
        Ie = float(np_.trapz(np_.array(mm["A_exact"]), grid))
        If = float(np_.trapz(np_.array(mm["A_final"]), grid))
        rep.eq("method_max.json sum_rule.sum_rule_exact_integrated",
               sr.get("sum_rule_exact_integrated"), Ie, 1e-9,
               "the deposited integral of A_exact is not the integral of the deposited "
               "A_exact (%.10f)" % Ie)
        rep.eq("method_max.json sum_rule.sum_rule_reconstructed_integrated",
               sr.get("sum_rule_reconstructed_integrated"), If, 1e-9,
               "the deposited integral of A_final is not the integral of the deposited "
               "A_final (%.10f)" % If)
        rep.eq("method_max.json sum_rule.sum_rule_analytic",
               sr.get("sum_rule_analytic"), mm["sumrule"], 1e-9,
               "the analytic identity in the sum_rule block disagrees with the legacy "
               "'sumrule' field of the same file")
        if "window_deficit_frac" in sr and "tail_outside_window" in sr:
            rep.eq("method_max.json sum_rule window deficit is the tail, not quadrature",
                   sr["tail_outside_window"],
                   float(sr["sum_rule_analytic"]) - Ie - float(sr.get("quad_error", 0.0)),
                   1e-6,
                   "the window/quadrature decomposition does not close: the deposited "
                   "tail plus the deposited quadrature residue is not the measured "
                   "deficit. The decomposition is the whole defence of the 2.7 per cent")
    else:
        rep.missing("method_max.json sum_rule block",
                    "the integrated sum rule added on 2026-09-18 is gone from "
                    "data/method_max.json; the figure caption's claim has no dataset again")

    # ---------------------------------------------------------------- 9.13 plumbing
    pyfiles = sorted(f for f in os.listdir(os.path.join(root, "src")) if f.endswith(".py"))
    rd = read(root, "README.md")
    m = re.search(r"(\d+)\s+computation & figure scripts", rd)
    if m:
        rep.eq_int("README.md script count", int(m.group(1)), len(pyfiles),
                   "README.md advertises a number of scripts that is not the number of "
                   ".py files in src/. A hand-maintained count goes stale the moment "
                   "anybody adds a script")
    else:
        rep.missing("README.md script count",
                    "the sentence carrying the script count could not be located")
    # every script that integrates must carry the numpy 1.x/2.x bridge, or the
    # documented environment cannot run the documented pipeline
    nobridge = []
    for f in pyfiles:
        t = read(root, "src/" + f)
        bridged = ("np.trapezoid = np.trapz" in t and "np.trapz = np.trapezoid" in t)
        if ("np.trapz(" in t or "np.trapezoid" in t) and not bridged:
            nobridge.append(f)
    rep.truth("every integrating script carries the numpy 1.x/2.x bridge",
              not nobridge, "%d scripts checked" % len(pyfiles),
              "THE DEPOSITED ENVIRONMENT CANNOT RUN THE DEPOSITED PIPELINE: %s call "
              "np.trapz or np.trapezoid without the compatibility bridge. np.trapezoid "
              "exists only from numpy 2.0 and np.trapz was REMOVED in numpy 2.0, so "
              "without the bridge no single numpy version runs the whole repository."
              % nobridge, n=max(len(pyfiles), 1))
    # the arXiv bundle: the counts printed in three documents must be the measured ones.
    # This sentence went stale once already and nobody noticed for a day.
    try:
        import tarfile
        same, diff, only = [], [], []
        with tarfile.open(os.path.join(root, "paper", "arxiv-submission.tar.gz")) as tf:
            for mem in tf.getmembers():
                if not mem.isfile():
                    continue
                relp = mem.name[2:] if mem.name.startswith("./") else mem.name
                loc = os.path.join(root, "paper", relp.replace("/", os.sep))
                data = tf.extractfile(mem).read()
                if not os.path.exists(loc):
                    only.append(relp)
                elif (open(loc, "rb").read().replace(chr(13).encode() + chr(10).encode(),
                                                     chr(10).encode())
                      == data.replace(chr(13).encode() + chr(10).encode(), chr(10).encode())):
                    same.append(relp)
                else:
                    diff.append(relp)
    except Exception as exc:                                        # noqa: BLE001
        rep.missing("arxiv bundle comparison", "cannot read the bundle: %s" % exc)
        same = diff = only = None
    if same is not None:
        WORDS = {2: "two", 8: "eight", 9: "nine", 36: "36", 39: "39"}
        kd = read(root, "docs/KNOWN_DISCREPANCIES.md")
        rep.truth("docs/KNOWN_DISCREPANCIES.md states the measured bundle comparison",
                  ("**%d are byte-identical**" % len(same)) in kd
                  and ("**%d differ**" % len(diff)) in kd
                  and (("**%d exist only in the bundle**" % len(only)) in kd
                       or ("**%s exist only in the bundle**"
                           % WORDS.get(len(only), len(only))) in kd),
                  "%d identical / %d differ / %d bundle-only" % (len(same), len(diff), len(only)),
                  "STALE DERIVED NUMBER in the document that records stale derived numbers: the "
                  "measured comparison is %d identical, %d differ, %d only in the bundle "
                  "(%s). Run python src/check_tarball.py."
                  % (len(same), len(diff), len(only), diff))
        rd2 = read(root, "README.md")
        rep.truth("README.md states the measured number of files ahead of the bundle",
                  ("%s files in `paper/` have since moved ahead"
                   % WORDS.get(len(diff), len(diff))) in rd2,
                  "%d files ahead" % len(diff),
                  "README.md states a number of files ahead of paper/arxiv-submission.tar.gz "
                  "that is not the measured %d" % len(diff))
        rep.truth("the arXiv bundle still carries the fabricated S(q,w) rows",
                  "figs/sqw_edges.dat" in diff,
                  "the bundle is behind on figs/sqw_edges.dat, as the documents say",
                  "figs/sqw_edges.dat now matches the bundle. Either the bundle was rebuilt "
                  "(delete this check and say so) or the repaired table was reverted.")

    # CI must invoke the guardian at the path where the guardian lives
    cirel = ".github/workflows/ci.yml"
    if os.path.exists(os.path.join(root, cirel.replace("/", os.sep))):
        ci = read(root, cirel)
        rep.truth("CI invokes the guardian at its real path",
                  "python src/verify.py" in ci,
                  "the workflow runs src/verify.py",
                  "the CI workflow does not run 'python src/verify.py'. Running "
                  "'python verify.py' from the repository root cannot work: the file "
                  "lives in src/. README.md calls this file 'what CI runs'.")
    else:
        rep.skip("CI workflow", "%s is not present in the working tree" % cirel)

    # ---------------------------------------------------------------- 9.14
    # data/charge_gap.json -- deposited 2026-09-18, and UNGUARDED until now.
    # An adversarial pass falsified Delta(L=12) from 4.968759 to 4.900000 and this
    # guardian still reported PASS, 0 FAIL: the only occurrence of "charge_gap" in
    # this file was an entry in the open-claims registry, which reads results.tex,
    # not the datum.  The file backs the Delta = 4.97 t that the manuscript prints,
    # so a silent corruption there is a corrupted number in a published paper.
    #
    # What is checked is NOT the printed value against itself -- that proves nothing.
    # Every row must be internally consistent (mu+ and mu- rebuilt from the three
    # sector energies, Delta rebuilt from mu+ and mu-, particle-hole symmetry), the
    # sector dimensions must be the binomials they claim to be, the ground-state
    # energies must agree with the SAME energies deposited by three other scripts in
    # three other files, and the Lieb-Wu L->infinity limit is recomputed here by
    # quadrature instead of copied.
    cg_rel = "data/charge_gap.json"
    if not os.path.exists(os.path.join(root, cg_rel.replace("/", os.sep))):
        rep.missing("data/charge_gap.json present",
                    "the file that backs Delta = 4.97 t is not in the deposit")
    else:
        cg = read_json(root, cg_rel)
        rows = cg["results"]
        rep.eq_int("charge_gap.json carries every size L=6,8,10,12",
                   [r["L"] for r in rows], [6, 8, 10, 12],
                   "the deposited charge gap no longer covers the four sizes the "
                   "scaling statement is made over")
        for r in rows:
            L, U = r["L"], r["U"]
            w = "charge_gap L=%d" % L
            rep.eq(w + ": mu+ = E0(N+1) - E0(N)", r["mu_plus"],
                   r["E0_Np1"] - r["E0_N"], 1e-12,
                   "the deposited mu+ is not the difference of the deposited sector "
                   "energies: the row was edited without recomputing it")
            rep.eq(w + ": mu- = E0(N) - E0(N-1)", r["mu_minus"],
                   r["E0_N"] - r["E0_Nm1"], 1e-12,
                   "the deposited mu- is not the difference of the deposited sector "
                   "energies")
            rep.eq(w + ": Delta = mu+ - mu-", r["Delta"],
                   r["mu_plus"] - r["mu_minus"], 1e-12,
                   "THE CHARGE GAP DOES NOT FOLLOW FROM THE CHEMICAL POTENTIALS IN ITS "
                   "OWN ROW -- which is exactly what a hand-edited Delta looks like")
            rep.eq(w + ": particle-hole symmetry mu+ + mu- = U",
                   r["mu_plus"] + r["mu_minus"], U, 1e-8,
                   "half filling on a bipartite ring forces mu+ + mu- = U exactly; a row "
                   "that violates it is not the half-filled Hubbard ring it claims to be")
            rep.eq(w + ": mu+_rel = mu+ - U/2", r["mu_plus_rel"],
                   r["mu_plus"] - U / 2.0, 1e-12,
                   "the particle-hole-symmetric level is not offset from mu+ by U/2")
            rep.eq(w + ": mu-_rel = mu- - U/2", r["mu_minus_rel"],
                   r["mu_minus"] - U / 2.0, 1e-12,
                   "the particle-hole-symmetric level is not offset from mu- by U/2")
            rep.eq(w + ": L*Delta", r["L_times_Delta"], L * r["Delta"], 1e-10,
                   "the tabulated L*Delta is not L times the tabulated Delta")
            rep.eq_int(w + ": dim(N) = C(L,L/2)^2", r["dim_N"], comb(L, L // 2) ** 2,
                       "the N-sector dimension is not the half-filled binomial: the row "
                       "does not describe the system it claims to")
            rep.eq_int(w + ": dim(N+1) = C(L,L/2+1)*C(L,L/2)", r["dim_Np1"],
                       comb(L, L // 2 + 1) * comb(L, L // 2),
                       "the (N+1)-sector dimension is not the binomial it must be")
            rep.eq_int(w + ": dim(N-1) = C(L,L/2-1)*C(L,L/2)", r["dim_Nm1"],
                       comb(L, L // 2 - 1) * comb(L, L // 2),
                       "the (N-1)-sector dimension is not the binomial it must be")
            if r.get("E0_N_dense") is not None:
                rep.eq(w + ": ARPACK E0(N) vs dense eigh", r["E0_N"], r["E0_N_dense"], 1e-10,
                       "the iterative and the dense ground-state energy of the same sector "
                       "disagree: one of the two is not converged")
        # cross-file: the SAME ground-state energy, written by four different scripts
        sd = dict((p["L"], p["E0"]) for p in read_json(root, "data/scaling_data.json")["points"])
        for r in rows:
            if r["L"] in sd:
                rep.eq("charge_gap L=%d: E0(N) == data/scaling_data.json" % r["L"],
                       r["E0_N"], sd[r["L"]], 1e-9,
                       "two deposited files disagree about the ground-state energy of the "
                       "same ring: charge_gap_ed.py and scaling_data.py must agree, or one "
                       "of them is not solving the stated Hamiltonian")
        for rel2 in ("data/sqw_L12.json", "data/akw_lanczos_L12.json"):
            r12 = [r for r in rows if r["L"] == 12]
            if r12 and os.path.exists(os.path.join(root, rel2.replace("/", os.sep))):
                rep.eq("charge_gap L=12: E0(N) == %s" % rel2, r12[0]["E0_N"],
                       read_json(root, rel2)["E0"], 1e-9,
                       "the charge-gap ground state is not the ground state the deposited "
                       "L=12 spectral functions were built on")
        # Lieb-Wu: recomputed here, not copied out of the file it is checking
        try:
            from scipy.integrate import quad as _quad
            from scipy.special import j1 as _j1
            Ulw, tlw = 8.0, 1.0
            def _lwint(w):
                x = w * Ulw / (2.0 * tlw)
                if x > 700.0:                 # exp overflows; the integrand is 0 there
                    return 0.0
                return _j1(w) / (w * (1.0 + math.exp(x)))
            integ = _quad(_lwint, 0.0, 200.0, limit=400)[0]
            lw = Ulw - 4.0 * tlw + 8.0 * tlw * integ
            rep.eq("Lieb-Wu Delta(L->inf) at U/t=8, recomputed by quadrature",
                   cg["lieb_wu"]["Delta_inf"], lw, 1e-7,
                   "the deposited Lieb-Wu limit is not what its own stated formula "
                   "produces: the finite-size gaps are being compared against a wrong "
                   "asymptote")
        except Exception as exc:                                    # noqa: BLE001
            rep.skip("Lieb-Wu quadrature", "scipy quadrature unavailable: %s" % exc)
        dl = [r["Delta"] for r in rows]
        rep.truth("charge gap decreases monotonically with L",
                  all(dl[i] > dl[i + 1] for i in range(len(dl) - 1)),
                  "Delta(L) = " + ", ".join("%.6f" % x for x in dl),
                  "the deposited charge gap is not monotone in L, which the finite-size "
                  "extrapolation assumes", n=max(len(dl) - 1, 1))
        rep.truth("every finite-size charge gap lies above the Lieb-Wu limit",
                  all(x > cg["lieb_wu"]["Delta_inf"] for x in dl),
                  "min Delta = %.6f > %.6f" % (min(dl), cg["lieb_wu"]["Delta_inf"]),
                  "a finite ring is reported with a SMALLER charge gap than the infinite "
                  "chain, which is the wrong sign for the 1/L correction", n=len(dl))
        cgtxt = read(root, cg_rel)
        rep.truth("charge_gap.json carries no absolute machine path",
                  ("C:" + chr(92)) not in cgtxt and "/Users/" not in cgtxt,
                  "provenance is repository-relative",
                  "a deposited file names a directory on the author's machine: both a "
                  "privacy leak and a pointer no reader can follow")


# ===========================================================================
#  SECTION 10 -- open claims, and the coverage floor
# ===========================================================================

def section10(rep, root):
    rep.head("(10) RETRACTED CLAIMS -- %d sentences that must stay out of the manuscript"
             % len(CLAIMS_RETRACTED))
    for key in sorted(CLAIMS_RETRACTED):
        claim_check(rep, root, key)


def section_c3(rep, root):
    r"""The C3 frontier sweep: 204 rows, 54 point files, deposited 2026-09-19.

    Until that date Sec. V of the manuscript was the one part of the paper a reader
    could not regenerate: its three tables came from a sweep that lived only in a
    working tree.  The sweep is now in data/c3_frontier/ and its code in src/frontier/.
    This section is what makes the deposit worth having.  It does NOT read the grid and
    echo it back.  It opens the 54 POINT files, re-derives every derived column of the
    204-row grid from the primitives measured in them (w_S, Lambda_in, Lambda_out,
    rel-L1), and then checks Theorem B itself on all 204 rows.  Six things are checked:

      (a) closure      -- every grid row is backed by a deposited point file, and every
                          deposited point file is either used or SUPERSEDED by a deeper
                          Krylov depth at the same (L, FR).  A file that is neither is a
                          FAIL: it would be a measurement in the deposit that the grid
                          silently ignores.
      (b) arithmetic   -- bound_leak, bound_trivial, cert, in/out, slack, the
                          non-triviality flag and the convergence flag are RECOMPUTED
                          here from Theorem B's own algebra and from the classifier, and
                          compared with what the grid stores.
      (c) Theorem B    -- ||phi||^2 (1-w_S) <= ||A-A_S||_1 <= min{leak, trivial}, in
                          relative form (1-w_S) <= rel-L1 <= cert, on every row that
                          carries a measured error.  This is the manuscript's central
                          inequality and this is the only place it is checked against a
                          measured error rather than against a synthetic sweep.
      (d) normalisation-- Lambda_hat = Lambda/||phi|| (NOT Lambda), re-derived from each
                          point file's own norm_phi2.  The factor sqrt(2) between the two
                          is exactly the sort of silent convention error that would make
                          every certificate in Sec. V wrong by 41 per cent.
      (e) the cost     -- the run cost is recomputed from the point files (sum of
                          t_point over the distinct (L, FR, n_l) points; max peak) and
                          checked against the numbers docs/C3_FRONTIER.md prints.  The
                          estimate this deposit replaces was 20-40 h and 9 GB; the
                          measurement is 0.8 h and 1.6 GiB, and a document that drifts
                          back towards the estimate now fails the build.
      (f) the vacuity  -- at the fraction the manuscript published for L=12 the
                          certificate must be VACUOUS (leak > trivial) at all four
                          resolutions.  That is Sec. V's headline negative result; if it
                          ever stops being true the section is wrong, not the data.

    The point files are opened, not stat'ed: a deleted or truncated one fails here.
    """
    import glob as _glob

    rep.head("(10b) C3 FRONTIER SWEEP -- 204 rows re-derived from the 54 deposited points")

    ddir = os.path.join(root, "data", "c3_frontier")
    if not os.path.isdir(ddir):
        rep.missing("data/c3_frontier",
                    "the C3 frontier deposit is absent. Sec. V of the manuscript is then "
                    "unreproducible again and its three tables have no source in this "
                    "repository. Restore data/c3_frontier/ and src/frontier/.")
        return

    grid_file = "data/c3_frontier/C3_grid.json"
    G = read_json(root, grid_file)
    grid, front = G["grid"], G["frontier"]

    # -- (a) closure ---------------------------------------------------------
    files = sorted(os.path.basename(p) for p in _glob.glob(os.path.join(ddir, "pt_*.json")))
    pts = {}
    for n in files:
        d = read_json(root, "data/c3_frontier/" + n)
        pts[n] = d
    keys_file = {}
    for n, d in pts.items():
        keys_file[n] = (d["meta"]["L"], round(d["FR"], 4), d["nl"])
    deepest = {}
    for n, (L, FR, nl) in keys_file.items():
        if (L, FR) not in deepest or nl > deepest[(L, FR)][1]:
            deepest[(L, FR)] = (n, nl)
    used = {v[0] for v in deepest.values()}
    keys_grid = {(g["L"], round(g["FR"], 4), g["nl"]) for g in grid}
    backed = {keys_file[n] for n in used}
    rep.truth("C3 closure: every grid row is backed by a deposited point file",
              keys_grid <= backed,
              "%d distinct (L,FR,n_l) rows, all present in data/c3_frontier/" % len(keys_grid),
              "grid rows with no deposited point file: %s" % sorted(keys_grid - backed),
              n=len(keys_grid))
    stray = []
    for n in files:
        if n in used:
            continue
        L, FR, nl = keys_file[n]
        win, wnl = deepest[(L, FR)]
        if not (win in used and wnl > nl):
            stray.append(n)
    rep.truth("C3 closure: no point file is orphaned", not stray,
              "%d point files: %d aggregated, %d superseded by a deeper n_l at the same "
              "(L,FR) -- %s" % (len(files), len(used), len(files) - len(used),
                                ", ".join(sorted(set(files) - used)) or "none"),
              "point files that are neither aggregated nor superseded: %s" % stray,
              n=len(files))
    # the row count is DERIVED from the etas each point actually carries, never typed
    etas = sorted({round(g["eta"], 4) for g in grid})
    expect = sum(1 for n in used for e in etas
                 if any(abs(r["eta"] - e) < 1e-9 and r["n"] == pts[n]["nl"]
                        for r in pts[n]["rows"]))
    rep.eq_int("C3_grid.json row count", len(grid), expect,
               "the grid no longer has one row per (deposited point) x (resolution). "
               "Either a point file was dropped or a row was added by hand")
    rep.eq_int("C3_grid.json frontier-block row count", len(front),
               len({(g["L"], round(g["eta"], 4)) for g in grid}),
               "the frontier block must carry exactly one crossing per (L, eta) cell")

    # -- (b) every derived column is the arithmetic of the primitives --------
    worst = (0.0, "")
    nbad = 0
    for g in grid:
        w = g["w_S"]
        B = (1.0 + math.sqrt(w)) * (math.sqrt(max(1.0 - w, 0.0)) + g["Lam_over_eta"])
        triv = 1.0 + w
        io = g["Lam_in"] / max(g["Lam_out"], 1e-300)
        want = dict(bound_leak=B, bound_trivial=triv, cert=min(B, triv), in_over_out=io)
        if g["relL1"]:
            want["slack"] = g["cert"] / g["relL1"]
            want["slack_leakbranch"] = g["bound_leak"] / g["relL1"]
        for k, v in want.items():
            d = abs(g[k] - v) / max(abs(v), 1e-300)
            if d > worst[0]:
                worst = (d, "L=%d FR=%.2f eta=%.2f %s: stored %.12g, recomputed %.12g"
                         % (g["L"], g["FR"], g["eta"], k, g[k], v))
            if d > 1e-12:
                nbad += 1
        if bool(B < triv) != bool(g["nontrivial"]):
            nbad += 1
        pr = g["plateau_relspread"]
        conv = ("CONVERGED" if (io < 0.05 and pr < 1e-3)
                else ("MARGINAL" if io < 0.30 else "NOT-CONVERGED"))
        if conv != g["conv"]:
            nbad += 1
            worst = (float("inf"), "L=%d FR=%.2f eta=%.2f conv: stored %s, classifier says %s"
                     % (g["L"], g["FR"], g["eta"], g["conv"], conv))
    rep.truth("C3 grid: every derived column recomputed from Theorem B's algebra",
              nbad == 0,
              "204 rows x 8 columns, max relative deviation %.2e" % worst[0],
              "%d column(s) do not follow from the stored primitives. Worst: %s"
              % (nbad, worst[1]), n=len(grid) * 8)

    # -- (c) Theorem B on every row with a measured error --------------------
    lo_bad, hi_bad, nrows = [], [], 0
    for g in grid:
        if not g["relL1"]:
            continue
        nrows += 1
        tag = "L=%d FR=%.2f eta=%.2f" % (g["L"], g["FR"], g["eta"])
        if g["relL1"] < (1.0 - g["w_S"]) - 1e-12:
            lo_bad.append(tag)
        if g["relL1"] > g["cert"] + 1e-12:
            hi_bad.append(tag)
    rep.truth("Theorem B lower bound (1-w_S) <= rel-L1 on every deposited row",
              not lo_bad, "%d rows, 0 violations" % nrows,
              "%d violation(s): %s" % (len(lo_bad), lo_bad[:4]), n=nrows)
    rep.truth("Theorem B upper bound rel-L1 <= min{leak, trivial} on every deposited row",
              not hi_bad, "%d rows, 0 violations" % nrows,
              "%d violation(s) -- THE THEOREM OF SEC. III FAILS ON THE DEPOSIT: %s"
              % (len(hi_bad), hi_bad[:4]), n=nrows)

    # -- (d) the normalisation of Lambda_hat ---------------------------------
    worstn = (0.0, "")
    nn = 0
    for n in sorted(used):
        d = pts[n]
        nphi = math.sqrt(d["meta"]["norm_phi2"])
        for r in d["rows"]:
            if r["n"] != d["nl"]:
                continue
            nn += 1
            lam = math.hypot(r["Lambda_K_in"], r["Lambda_K_out"])
            got = r["LamHat_over_eta"] * r["eta"] * nphi
            dev = abs(got - lam) / max(lam, 1e-300)
            if dev > worstn[0]:
                worstn = (dev, "%s eta=%.2f: Lambda_hat*eta*||phi|| = %.12g, "
                               "hypot(in,out) = %.12g" % (n, r["eta"], got, lam))
    rep.truth("Lambda_hat = Lambda/||phi|| and Lambda^2 = in^2 + out^2, on every point",
              worstn[0] < 1e-9, "%d (point, eta) pairs, max deviation %.2e" % (nn, worstn[0]),
              "the certificate's normalisation does not hold: %s. A missing ||phi|| is a "
              "factor sqrt(2) on every bound in Sec. V" % worstn[1], n=nn)

    # -- (e) the measured cost, and the document that prints it --------------
    uniq = {}
    for g in grid:
        uniq[(g["L"], g["FR"], g["nl"])] = (g["t_point"], g["peak_GiB"])
    tot = sum(v[0] for v in uniq.values())
    peak = max(v[1] for v in uniq.values())
    cens = {}
    for g in grid:
        cens[g["conv"]] = cens.get(g["conv"], 0) + 1
    measured = {
        "C3_COST_SECONDS": round(tot, 1),
        "C3_COST_HOURS": round(tot / 3600.0, 3),
        "C3_PEAK_GIB": round(peak, 4),
        "C3_UNIQUE_POINTS": len(uniq),
        "C3_GRID_ROWS": len(grid),
        "C3_POINT_FILES": len(files),
        "C3_CONVERGED": cens.get("CONVERGED", 0),
        "C3_MARGINAL": cens.get("MARGINAL", 0),
        "C3_NOTCONVERGED": cens.get("NOT-CONVERGED", 0),
    }
    doc_rel = "docs/C3_FRONTIER.md"
    doc_path = os.path.join(root, "docs", "C3_FRONTIER.md")
    if not os.path.isfile(doc_path):
        rep.missing(doc_rel, "the C3 deposit has no reproduction document; Sec. V would be "
                             "regenerable only by reading the source of src/frontier/")
    else:
        txt = read(root, doc_rel)
        for key, want in sorted(measured.items()):
            m = re.search(r"^\s*" + key + r"\s*=\s*([0-9.]+)\s*$", txt, re.M)
            if m is None:
                rep.missing("%s: %s" % (doc_rel, key),
                            "the machine-readable cost stanza of %s has no line %r. The "
                            "estimate this deposit replaced (20-40 h, 9 GB) was 25-50x "
                            "pessimistic; the measurement has to stay checkable."
                            % (doc_rel, key))
                continue
            rep.eq("%s: %s" % (doc_rel, key), float(m.group(1)), float(want), tol=5e-4,
                   defect="the document prints a cost/census the deposited points do not "
                          "support (recomputed from data/c3_frontier/ on this run)")

    # -- (f) the published fraction of L=12 is uncertified --------------------
    pub = [g for g in grid if g["L"] == 12 and abs(g["FR"] - 0.18) < 1e-9]
    rep.truth("Sec. V headline: at L=12, FR=0.18 the certificate is VACUOUS at every eta",
              bool(pub) and all(g["bound_leak"] > g["bound_trivial"] for g in pub),
              "4 resolutions, leak/trivial = %s"
              % ", ".join("%.2f" % (g["bound_leak"] / g["bound_trivial"]) for g in pub),
              "the negative result of Sec. V no longer holds on the deposited grid: %s"
              % [(g["eta"], g["bound_leak"], g["bound_trivial"]) for g in pub],
              n=len(pub))

    # -- the proof ladder, and the one stale artefact in the deposit ---------
    # ladder_L4.json predates a change to the integration window of
    # src/frontier/fron_lib.make_grid by four minutes; the registry pins that, so the
    # deposit declares its one stale artefact instead of hiding it.  ladder_L6 and
    # ladder_L8 postdate the change and are checked against the CURRENT convention.
    lad = {L: read_json(root, "data/c3_frontier/ladder_L%d.json" % L) for L in (4, 6, 8)}
    ng = {}
    for L, d in lad.items():
        per = d["per_eta"]
        rows = d["cases"][0]["rows"]
        ng[L] = {round(r["eta"], 4): r["ngrid"] for r in rows}
        # ngrid = 2*ntail + ceil((hi-lo)/(eta/per_eta)) + 1 is a pure function of the
        # pole span and eta, so it is the cheapest witness of which window was used.
        ratio = ng[L][0.5] / float(ng[L][0.05])
        rep.truth("ladder_L%d.json: per_eta = %g recorded and its grid is monotone in eta"
                  % (L, per), 0.0 < ratio < 1.0,
                  "ngrid(eta=0.05)=%d > ngrid(eta=0.50)=%d" % (ng[L][0.05], ng[L][0.5]),
                  "the omega grid does not coarsen with eta: %s" % ng[L], n=len(ng[L]))
    rep.known("ladder_L4.ngrid_eta050", "ladder_L4.json: integration grid at eta=0.50",
              float(ng[4][0.5]))
    vf = read_json(root, "data/c3_frontier/verdict_frontier.json")
    rep.eq_int("verdict_frontier.json rows", len(vf),
               len({r["L"] for r in vf}) * len({r["eta"] for r in vf})
               * len({r["thr"] for r in vf}),
               "the proof-ladder verdict must carry one row per (L, eta, threshold)")
    rep.truth("verdict_frontier.json covers the three deposited ladders",
              {r["L"] for r in vf} == set(lad), "L = %s" % sorted(lad),
              "verdict_frontier.json was computed from a different set of ladder files "
              "than the three deposited ones: %s vs %s" % (sorted({r["L"] for r in vf}),
                                                           sorted(lad)), n=len(lad))

    # -----------------------------------------------------------------------
    # (h) the L=14 ranking, re-run (2026-09-28, late).  The L=14 certificate point takes
    #     its subspace from the stored ranking ckpt/L14_order.npy; until this date the
    #     K=18 protocol had been re-run only at L=10.  src/frontier/rerank_L14.py re-runs
    #     it at L=14 from the stored ground state and compares the top-k selections as
    #     sets; the record is data/c3_frontier/published_fraction/L14_rerank.json.  This
    #     guardian cannot repeat the 49-minute ranking; it checks that the record is the
    #     one the text describes (set overlap 1 with nothing unshared at the published
    #     |S| and at every other cut, the same captured weight, an unambiguous cut at the
    #     published fraction, the deposited checkpoints by SHA-256) and that Sec. S3 says so.
    # -----------------------------------------------------------------------
    try:
        rr = read_json(root, "data/c3_frontier/published_fraction/L14_rerank.json")
    except Exception as exc:                                    # noqa: BLE001
        rep.missing("data/c3_frontier/published_fraction/L14_rerank.json", "cannot read: %s" % exc)
        rr = None
    if rr:
        rep.eq_int("L14_rerank.json: L", int(rr["L"]), 14, "not the L=14 record")
        rep.eq_int("L14_rerank.json: sector dimension", int(rr["D"]), 10306296,
                   "not the (N+1) sector of L=14, U/t=8, half filling")
        for k_, v_ in (("K", 18), ("dt", 0.5), ("sub", 16), ("m", 6)):
            rep.eq("L14_rerank.json: protocol %s" % k_, float(rr["protocol"][k_]), float(v_), 0.0,
                   "the ranking was re-run with a protocol other than the repository's K=18 one")
        rows_ = rr["rows"]
        rep.truth("L14_rerank.json: the published |S| = 824504 is one of the cuts",
                  any(int(r["k"]) == 824504 for r in rows_), "present",
                  "the re-run does not compare the selection at the published L=14 count")
        rep.truth("L14_rerank.json: at least four cuts compared", len(rows_) >= 4, "%d" % len(rows_),
                  "fewer cuts than the text reports")
        for r in rows_:
            rep.eq("L14_rerank.json k=%d: set overlap with the stored ranking" % int(r["k"]),
                   float(r["set_overlap"]), 1.0, 0.0,
                   "the re-run ranking does not reproduce the stored selection as a set")
            rep.eq_int("L14_rerank.json k=%d: determinants not shared" % int(r["k"]), int(r["n_only_stored"]), 0,
                       "the two selections differ")
            rep.eq_int("L14_rerank.json k=%d: n_common == k" % int(r["k"]), int(r["n_common"]), int(r["k"]),
                       "the overlap is not the whole selection")
            rep.eq("L14_rerank.json k=%d: captured weight, stored vs re-run" % int(r["k"]),
                   float(r["w_S_rerun"]), float(r["w_S_stored"]), 1e-12,
                   "the two selections carry different captured weight")
            if int(r["k"]) == 824504:
                rep.truth("L14_rerank.json k=824504: the cut is unambiguous", not r["cut_ambiguous"],
                          "%d/%d tied at the cut, all inside" % (r["tied_at_cut_inside"], r["tied_at_cut_total"]),
                          "a determinant tied at the cut weight lies outside the published selection")
                rep.eq("L14_rerank.json k=824504: captured weight is the L14_FR008 one",
                       float(r["w_S_stored"]), 0.998387504078, 1e-11,
                       "the stored captured weight at the published count moved")
        ck_ = read(root, "ckpt/README.txt")
        for key_, nm_ in (("L14_gs_sha256", "L14_gs.npz"), ("L14_order_sha256", "L14_order.npy")):
            rep.truth("L14_rerank.json: %s hashed as ckpt/README.txt records" % nm_,
                      str(rr["provenance"].get(key_, "")) in ck_ and len(str(rr["provenance"].get(key_, ""))) == 64,
                      str(rr["provenance"].get(key_, ""))[:12],
                      "the re-run did not read the deposited checkpoint (its SHA-256 is not the recorded one)")
        smf_ = " ".join(read(root, "paper/sm_frontier.tex").split())
        rep.truth("Sec. S3 reports the L=14 re-run and its set overlap",
                  "reproduces the stored selection as a set to $1.000000$" in smf_
                  and "$\\lvert S\\rvert=824\\,504$ and at the fractions $0.04$, $0.16$ and $0.32$" in smf_,
                  "stated", "Sec. S3 no longer reports the re-run ranking at L=14")
        for rel_, bad_ in (("paper/sm_frontier.tex", "was \\emph{not} re-run at that size"),
                           ("paper/sec_5_body.tex", "was not re-run"),
                           ("paper/sec_1.tex", "not re-ranked"), ("paper/sec_10.tex", "not re-ranked")):
            rep.truth("%s: the superseded 'not re-run' caveat is gone" % rel_,
                      bad_ not in " ".join(read(root, rel_).split()), "absent",
                      "the text still says the L=14 ranking was not re-run")


def coverage_floor(rep):
    """A guardian that only checks what it finds cannot see a deleted artefact."""
    rep.head("(11) COVERAGE FLOOR -- did every section actually do its work?")
    have = dict((t, n) for t, n, _ in rep.per_section)
    for prefix, floor in SECTION_FLOORS:
        got = None
        for t, n in have.items():
            if t.startswith(prefix):
                got = n
                break
        if got is None:
            rep.emit(FAIL, "section %s ran" % prefix,
                     "SECTION MISSING: no section beginning %r reported any assertion. "
                     "A guardian cannot certify a repository it did not read." % prefix, n=1)
        elif got < floor:
            rep.emit(FAIL, "section %s coverage floor" % prefix,
                     "COVERAGE COLLAPSED: %d numeric assertions against a floor of %d. "
                     "Checks do not vanish on their own -- an artefact (a whole plotted "
                     "series, a table, a dataset) has been deleted, and every check that "
                     "read it went quiet. Find what is missing; do not lower the floor."
                     % (got, floor), n=1)
        else:
            rep.emit(PASS, "section %s coverage floor" % prefix,
                     "%d assertions >= floor %d" % (got, floor), n=1)


# ===========================================================================
#  MAIN
# ===========================================================================

def main(argv=None):
    ap = argparse.ArgumentParser(description="Adversarial reproducibility guardian.")
    ap.add_argument("--root", default=None, help="repository root (default: parent of src/)")
    ap.add_argument("--fast", action="store_true", help="skip the L=10 exact diagonalisations")
    ap.add_argument("-v", "--verbose", action="store_true", help="print PASS and XFAIL lines too")
    args = ap.parse_args(argv)

    root = args.root or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    root = os.path.abspath(root)
    rep = Report(verbose=args.verbose)
    t0 = time.time()

    print("verify.py -- adversarial reproducibility guardian")
    print("repository root: %s" % root)
    print("mode: %s" % ("fast (no L=10 exact diagonalisation)" if args.fast else "full"))

    # -----------------------------------------------------------------------
    rep.head("(1) RECOMPUTED FROM FIRST PRINCIPLES -- half-filled Hubbard, U/t=8")
    # -----------------------------------------------------------------------
    sizes = [4, 6, 8] + ([] if args.fast else [10])
    OBC, PBC, E0_OBC, E0_PBC = {}, {}, {}, {}
    for L in sizes:
        OBC[L], E0_OBC[L] = hubbard_f1(L, 8.0, periodic=False)
        PBC[L], E0_PBC[L] = hubbard_f1(L, 8.0, periodic=True)
        print("  L=%-2d  OBC F1=%.6f E0=%.6f   |   PBC F1=%.6f E0=%.6f   (gap %.6f)"
              % (L, OBC[L], E0_OBC[L], PBC[L], E0_PBC[L], abs(OBC[L] - PBC[L])))
    if args.fast:
        rep.skip("ED L=10", "run without --fast to recompute F1(L=10) for both systems")

    # -----------------------------------------------------------------------
    rep.head("(2) RECOMPUTED FROM FIRST PRINCIPLES -- the geminal (APSG) witness")
    # -----------------------------------------------------------------------
    W = {K: witness_row(K, np.pi / 4) for K in (1, 2, 3, 4)}
    print("   K  M  Ne     F1      F1(gamma)   2N_u      F2       F3    |S|raw |S|eps  chi(2K) chi(max)")
    for K in (1, 2, 3, 4):
        w = W[K]
        print("  %2d %3d %3d  %8.5f %9.5f %9.5f %8.5f %8.5f %6d %6d %7d %8d"
              % (K, w["M"], w["Ne"], w["F1"], w["F1_gamma"], w["twoNu"], w["F2"], w["F3"],
                 w["S_raw"], w["S_eps"], w["chi_central"], w["chi_max"]))
    for K in (1, 2, 3, 4):
        w = W[K]
        rep.eq("witness K=%d  F1 == 4K" % K, w["F1"], 4.0 * K, 1e-9,
               "the Majorana antiflatness of the theta=pi/4 geminal product is not 4K")
        rep.eq("witness K=%d  F1 == 4 tr[gamma(1-gamma)]" % K, w["F1_gamma"], w["F1"], 1e-9,
               "the one-body collapse F1 = 4 tr[gamma(1-gamma)] fails on the witness state")
        rep.eq("witness K=%d  F1 == 2N_u" % K, w["twoNu"], w["F1"], 1e-9,
               "the identity F1 = 2N_u fails on the witness state")
        rep.eq("witness K=%d  F2 == 4K" % K, w["F2"], 4.0 * K, 1e-9,
               "second-order antiflatness is not 4K at theta=pi/4")
        rep.eq_int("witness K=%d  raw support == 2^K" % K, w["S_raw"], 2 ** K,
                   "the pairing-basis determinant support of the witness is not 2^K")
        rep.eq_int("witness K=%d  eps-support == 2^K" % K, w["S_eps"], 2 ** K,
                   "the eps-support of the witness is not 2^K")
        rep.eq_int("witness K=%d  minimal MPS bond dimension chi == 2" % K, w["chi_max"], 2,
                   "the minimal MPS bond dimension of the witness (max Schmidt rank over all "
                   "spin-orbital cuts) is not 2 -- the theorem's chi=2 is this quantity")
        rep.eq_int("witness K=%d  central-cut rank == %d" % (K, 2 if K % 2 else 1),
                   w["chi_central"], 2 if K % 2 else 1,
                   "the Schmidt rank at the single cut position 2K (what data/apsg_witness.json "
                   "stores under the name 'chi') changed")
    # theta sweep at K=3
    TH = [0.02, 0.10, 0.25, 0.40, 0.50]
    SW = {t: witness_row(3, t * np.pi) for t in TH}

    # -----------------------------------------------------------------------
    rep.head("(3) WITNESS ARTEFACTS -- data/apsg_witness.json, theorem display, figure")
    # -----------------------------------------------------------------------
    try:
        aw = read_json(root, "data/apsg_witness.json")["rows"]
    except Exception as exc:                                    # noqa: BLE001
        rep.missing("data/apsg_witness.json", "cannot read: %s" % exc)
        aw = []
    byK = {r["K"]: r for r in aw if r.get("theta") == "pi/4"}
    for K in (1, 2, 3, 4):
        if K not in byK:
            rep.missing("apsg_witness.json K=%d" % K, "row missing from the deposited witness")
            continue
        r, w = byK[K], W[K]
        rep.eq("apsg_witness.json K=%d FAF1" % K, r["FAF1"], w["F1"], 1e-9,
               "deposited FAF_1 disagrees with the recomputed geminal state")
        rep.eq("apsg_witness.json K=%d twoNu" % K, r["twoNu"], w["twoNu"], 1e-9,
               "deposited 2N_u disagrees with the recomputed geminal state")
        rep.eq("apsg_witness.json K=%d F2" % K, r["F2"], w["F2"], 1e-9,
               "deposited F_2 disagrees with the recomputed geminal state")
        rep.eq("apsg_witness.json K=%d F3" % K, r["F3"], w["F3"], 1e-9,
               "deposited F_3 disagrees with the recomputed geminal state")
        rep.eq_int("apsg_witness.json K=%d Sdet" % K, r["Sdet"], w["S_eps"],
                   "deposited determinant support disagrees with the recomputed geminal state")
        rep.eq_int("apsg_witness.json K=%d chi (central cut)" % K,
                   chi_of(r, "apsg_witness.json K=%d" % K)[0], w["chi_central"],
                   "deposited 'chi' disagrees with the Schmidt rank at the cut apsg_witness.py "
                   "actually takes (position 2K)")
        rep.eq_int("apsg_witness.json K=%d M" % K, r["M"], w["M"], "mode count changed")
        rep.eq_int("apsg_witness.json K=%d Ne" % K, r["Ne"], w["Ne"], "electron count changed")
        # 2026-09-28 (late): apsg_witness.py now records chi_max, the largest Schmidt rank over
        # all M-1 cuts of the pairing ordering -- the chi = 2 the theorem and Fig. S1 print --
        # so that number is generated by the deposited diagonalization and not declared.
        rep.eq_int("apsg_witness.json K=%d chi_max (largest Schmidt rank over all cuts)" % K,
                   int(r.get("chi_max", -1)), w["chi_max"],
                   "the deposited maximum over cuts is absent or disagrees with recomputation: the "
                   "printed chi = 2 would be declared, not generated")
    sweep_rows = {r["theta"]: r for r in aw if r.get("K") == 3 and r.get("theta") != "pi/4"}
    for t in TH:
        key = "%gpi" % t
        if key not in sweep_rows:
            rep.missing("apsg_witness.json theta=%s" % key, "sweep row missing")
            continue
        r, w = sweep_rows[key], SW[t]
        rep.eq("apsg_witness.json theta=%s FAF1" % key, r["FAF1"], w["F1"], 1e-9,
               "deposited theta-sweep FAF_1 disagrees with recomputation")
        rep.eq("apsg_witness.json theta=%s twoNu" % key, r["twoNu"], w["twoNu"], 1e-9,
               "deposited theta-sweep 2N_u disagrees with recomputation")
        rep.eq("apsg_witness.json theta=%s F2" % key, r["F2"], w["F2"], 1e-9,
               "deposited theta-sweep F_2 disagrees with recomputation")
        rep.eq_int("apsg_witness.json theta=%s Sdet" % key, r["Sdet"], w["S_eps"],
                   "deposited theta-sweep support disagrees with recomputation")
        rep.eq_int("apsg_witness.json theta=%s chi" % key,
                   chi_of(r, "apsg_witness.json theta=%s" % key)[0], w["chi_central"],
                   "deposited theta-sweep central-cut rank disagrees with recomputation")
        rep.eq_int("apsg_witness.json theta=%s chi_max (largest Schmidt rank over all cuts)" % key,
                   int(r.get("chi_max", -1)), w["chi_max"],
                   "the deposited maximum over cuts is absent or disagrees with recomputation "
                   "(2 at every interior angle, 1 at theta = pi/2)")
    # the two sentences that print the central-cut sequence beside the maximum (2026-09-28 (late))
    _cc = [str(W[K]["chi_central"]) for K in (1, 2, 3, 4)]
    for _rel, _pat in (("paper/sm_nogo.tex", r"at the central cut alone it is \$(\d)\$, \$(\d)\$, \$(\d)\$, \$(\d)\$ for \$K=1,\\dots,4\$"),
                       ("paper/sec_6_app.tex", r"central spin-orbital cut alone is \$(\d)\$, \$(\d)\$, \$(\d)\$, \$(\d)\$ for \$K=1,\\dots,4\$")):
        _t = " ".join(read(root, _rel).split())
        _m = re.search(_pat, _t)
        if not _m:
            rep.missing("%s: central-cut sequence" % _rel, "the sentence giving the central-cut ranks is gone")
        else:
            rep.truth("%s: central-cut ranks 2, 1, 2, 1 are the recomputed ones" % _rel, list(_m.groups()) == _cc,
                      ", ".join(_cc), "the printed central-cut sequence is not the recomputed one")
        rep.truth("%s: chi = 2 is named as the largest Schmidt rank over all cuts" % _rel,
                  "largest Schmidt rank over all cuts" in _t, "stated",
                  "the text no longer says which quantity the printed chi = 2 is")

    # --- the theorem display in paper/theorem_b2.tex ---
    try:
        th = read_body(root)                  # was paper/theorem_b2.tex (v2)
    except Exception as exc:                                    # noqa: BLE001
        rep.missing("paper/theorem_b2.tex", "cannot read: %s" % exc)
        th = ""
    thf = re.sub(r"\s+", " ", th)
    m = re.search(r"\\mathcal F_1\s*=\s*2N_\{\s*\\mathrm\s*\{?u\}?\s*\}\s*=\s*(\d+)\s*K", thf)
    if m:
        rep.eq("theorem display: F_1 = 2N_u = <c>K", float(m.group(1)), W[1]["F1"] / 1.0, 1e-9,
               "the theorem states F_1 = %sK but the recomputed witness gives F_1 = %g at K=1"
               % (m.group(1), W[1]["F1"]))
    else:
        rep.missing("theorem display: F_1 = 2N_u = 4K",
                    "the display could not be found in paper/theorem_b2.tex -- the theorem "
                    "statement this guardian is supposed to protect is gone or was reworded")
    m = re.search(r"\\mathcal F_2\s*=\s*(\d+)\s*K", thf)
    if m:
        rep.eq("theorem display: F_2 = <c>K", float(m.group(1)), W[1]["F2"], 1e-9,
               "the theorem states F_2 = %sK but the recomputed witness gives F_2 = %g at K=1"
               % (m.group(1), W[1]["F2"]))
    else:
        rep.missing("theorem display: F_2 = 4K", "display not found in paper/theorem_b2.tex")
    m = re.search(r"\|\\mathcal S\|_\{\\mathrm\{pair\}\}\s*=\s*(\d+)\^\{?K\}?", thf)
    if m:
        rep.eq_int("theorem display: |S|_pair = <b>^K", int(m.group(1)), W[1]["S_eps"],
                   "the theorem states |S|_pair = %s^K but the recomputed witness support at "
                   "K=1 is %d" % (m.group(1), W[1]["S_eps"]))
    else:
        rep.missing("theorem display: |S|_pair = 2^K", "display not found in paper/theorem_b2.tex")
    m = re.search(r"\|\\mathcal S\|_\{\\mathrm\{pair\}\}\s*=\s*\d+\^\{?K\}?\s*,\s*"
                  r"\\qquad\s*\\chi\s*=\s*(\d+)", thf)
    if m:
        rep.eq_int("theorem display: chi = <v>", int(m.group(1)), W[1]["chi_max"],
                   "the theorem states chi = %s; the minimal MPS bond dimension of the "
                   "recomputed witness (max Schmidt rank over all cuts) is %d"
                   % (m.group(1), W[1]["chi_max"]))
    else:
        rep.missing("theorem display: chi = 2",
                    "the chi entry of the witness theorem display could not be found -- the "
                    "single number that makes the witness a statement about cost is missing")

    # --- fig_witness_native.tex ---
    try:
        fw = lines(root, "paper/figs/fig_witness_native.tex")
    except Exception as exc:                                    # noqa: BLE001
        rep.missing("paper/figs/fig_witness_native.tex", "cannot read: %s" % exc)
        fw = []
    wp = dict(panel_blocks(fw))
    # panel (a): the four resources against K
    series_a = ((r"support \$\|\\mathcal\{S\}\|=2\^\{K\}\$", "S_eps", "determinant support"),
                (r"\$\\mathcal\{F\}_1=2N_\{\\mathrm u\}=4K\$", "F1", "one-body magic"),
                (r"\$\\mathcal\{F\}_2=4K\$", "F2", "second-order magic"),
                (r"bond dim\.\\? \$\\chi=2\$", "chi_max", "minimal MPS bond dimension"))
    blk = "\n".join(wp.get("a", []))
    if not blk:
        rep.missing("fig_witness_native.tex panel (a)", "panel not found")
    for pat, key, human in series_a:
        m = re.search(r"(coordinates *\{[^}]*\})[^\n]*\n?\s*\\addlegendentry\{" + pat + r"\}", blk)
        if not m:
            rep.missing("fig_witness (a) %s series" % human,
                        "the %s series is gone from the witness figure" % human)
            continue
        for (x, y) in coords_of(m.group(1)):
            K = int(round(x))
            if K not in W:
                continue
            rep.eq("fig_witness (a) %s at K=%d" % (human, K), y, float(W[K][key]), 5e-4,
                   "the witness figure plots %g for the %s at K=%d; the recomputed geminal "
                   "state gives %g" % (y, human, K, float(W[K][key])))
    # panel (b): the theta sweep at K=3
    blk = "\n".join(wp.get("b", []))
    for pat, key, human in ((r"\$\\mathcal\{F\}_1=2N_\{\\mathrm u\}\$", "F1", "one-body magic"),
                            (r"\$\\mathcal\{F\}_2\$", "F2", "second-order magic")):
        m = re.search(r"(coordinates *\{[^}]*\})[^\n]*\n?\s*\\addlegendentry\{" + pat + r"\}", blk)
        if not m:
            rep.missing("fig_witness (b) %s series" % human, "series not found")
            continue
        for (x, y) in coords_of(m.group(1)):
            t = round(x, 2)
            if t not in SW:
                continue
            rep.eq("fig_witness (b) %s at theta/pi=%.2f" % (human, t), y, SW[t][key], 5e-5,
                   "the witness figure plots %g; the recomputed theta sweep gives %g"
                   % (y, SW[t][key]))

    # -----------------------------------------------------------------------
    rep.head("(4) OPEN CHAIN vs PERIODIC RING -- the 2026-08-17 defect, both directions")
    # -----------------------------------------------------------------------
    rmj = read_json(root, "data/resource_master.json")
    lvc = read_json(root, "data/ladder_vs_chain.json")
    scal = read_json(root, "data/scaling_data.json")

    # 4.1 data/ladder_vs_chain.json -- declared OBC by its own provenance note
    for L in (6, 8, 10):
        key = "chain-L%d" % L
        row = lvc["results"].get(key)
        if row is None:
            rep.missing("ladder_vs_chain %s" % key, "chain row missing")
            continue
        if L in OBC:
            system_check(rep, "ladder_vs_chain.json %s F1" % key, row["F1"], L, "OBC",
                         OBC, PBC, ndp=None, where="data/ladder_vs_chain.json:%s.F1" % key)
        else:
            rep.skip("ladder_vs_chain.json %s F1" % key, "L=10 ED skipped (--fast)")

    # 4.2 data/resource_master.json -- the 'chain2' family is the OBC chain
    chain2 = [p for p in rmj["points"] if p["family"] == "chain2"]
    chain_by_S = {}
    for L in (6, 8, 10):
        row = lvc["results"].get("chain-L%d" % L)
        if row:
            chain_by_S[row["Sdet"]] = L
    for p in chain2:
        L = chain_by_S.get(p["S"])
        if L is None:
            rep.missing("resource_master chain2 point |S|=%r" % p["S"],
                        "no chain row in data/ladder_vs_chain.json has this support; the "
                        "chain family of the master dataset no longer matches the chain "
                        "dataset")
            continue
        if L in OBC:
            system_check(rep, "resource_master.json chain2 L=%d F1" % L, p["F1"], L, "OBC",
                         OBC, PBC, ndp=None,
                         where="data/resource_master.json points[family=chain2].F1")
        else:
            rep.skip("resource_master.json chain2 L=%d F1" % L, "L=10 ED skipped (--fast)")

    # 4.3 data/scaling_data.json -- the PERIODIC ring (akw_lanczos.hop, j=(i+1)%L)
    for pt in scal["points"]:
        L = pt["L"]
        if L in PBC:
            system_check(rep, "scaling_data.json L=%d FAF" % L, pt["FAF"], L, "PBC",
                         OBC, PBC, ndp=None, where="data/scaling_data.json points[L=%d].FAF" % L)
            rep.eq("scaling_data.json L=%d E0" % L, pt["E0"], E0_PBC[L], 1e-5,
                   "the deposited ring ground-state energy does not match exact "
                   "diagonalisation of the periodic ring")
        else:
            rep.skip("scaling_data.json L=%d FAF" % L,
                     "L=%d not recomputed here (%s)" % (L, "--fast" if L == 10 else
                                                        "beyond the seconds-scale ED budget"))
        rep.eq("scaling_data.json L=%d FAF_per_site" % L, pt["FAF_per_site"], pt["FAF"] / L, 1e-12,
               "FAF_per_site is not FAF/L")
        n = L // 2
        rep.eq_int("scaling_data.json L=%d (N+1) sector dim" % L, pt["Np1_sector"],
                   comb(L, n + 1) * comb(L, n),
                   "the deposited (N+1) sector dimension is not C(L,N/2+1)*C(L,N/2)")
        s_abs = pt["frac"] * pt["Np1_sector"]
        rep.eq("scaling_data.json L=%d frac*dim is an integer" % L, s_abs, round(s_abs), 1e-6,
               "frac is |S|/dim with |S| an integer count; frac*dim is not an integer")

    # 4.4 paper/figs/fig_resource_master_native.tex -- chain markers are OBC
    frm = lines(root, "paper/figs/fig_resource_master_native.tex")
    frm_panels = dict(panel_blocks(frm))
    chain_lines_b = [ln for ln in frm_panels.get("b", [])
                     if "rmChn" in ln and "coordinates" in ln]
    if not chain_lines_b:
        rep.missing("fig_resource_master_native.tex panel (b) chain series",
                    "the rmChn (open-chain) series is gone from panel (b) of the master figure")
    for ln in chain_lines_b:
        for (x, y) in coords_of(ln):
            L = chain_by_S.get(int(round(y)))
            if L is None or L not in OBC:
                continue
            system_check(rep, "fig_resource_master (b) chain L=%d F1" % L, x, L, "OBC",
                         OBC, PBC, ndp=3,
                         where="paper/figs/fig_resource_master_native.tex rmChn")

    # 4.5 paper/figs/fig_decoupling_native.tex -- Hubbard chain series are OBC
    fdc = lines(root, "paper/figs/fig_decoupling_native.tex")
    # The rebuilt fragment writes the stratum legends with a thin-space brace,
    # \addlegendentry{$(6{,}6)$}, and repeats each in panels (b) and (c).  The old
    # spelling matched nothing and the check went quiet -- which is exactly the failure
    # mode the coverage floor exists to catch, and it did catch it.
    # The rebuilt fragment puts the stratum LEGENDS in panel (c), where the abscissa is
    # chi, and repeats the same strata unlabelled in panel (b), where the abscissa is F1.
    # Taking the line next to the legend therefore read chi as if it were magic.  The
    # strata are identified here by their DATA instead: a panel-(b) series belongs to the
    # (L,N) stratum whose (FAF, S) pairs it reproduces.  Nothing is assumed about colours.
    _ce_rows = read_json(root, "data/cost_vs_ent.json")["rows"]
    _strata = {}
    for _r in _ce_rows:
        _strata.setdefault((_r["L"], _r["N"]), []).append(_r)
    _dec_b = {}
    for _ln in fdc:
        _cs = coords_of(_ln)
        if len(_cs) != 5:
            continue
        for _k, _rs in _strata.items():
            if all(any(abs(_r["FAF"] - x) < 6e-4 and abs(_r["S"] - y) < 0.51 for _r in _rs)
                   for (x, y) in _cs):
                _dec_b.setdefault(_k, _cs)
    rep.eq_int("fig_decoupling panel (b) plots all six (L,N) strata", len(_dec_b), 6,
               "a stratum series of the F1 panel could not be matched to any (L,N) group "
               "of data/cost_vs_ent.json; found %s" % sorted(_dec_b))
    for legend, L in (("$(6{,}6)$", 6), ("$(8{,}8)$", 8)):
        cs = _dec_b.get((L, L))
        if cs is None or L not in OBC:
            if L in OBC:
                rep.missing("fig_decoupling series %s" % legend,
                            "no panel-(b) series reproduces the (%d,%d) rows of "
                            "data/cost_vs_ent.json" % (L, L))
            continue
        # the U/t=8 point is the one whose F1 matches either system's value
        cand = [c for c in cs if min(abs(c[0] - OBC[L]), abs(c[0] - PBC[L])) < 5e-3]
        if not cand:
            rep.missing("fig_decoupling series %s U=8 point" % legend,
                        "no point of this series is near either the OBC (%.4f) or the PBC "
                        "(%.4f) value of F1(L=%d,U=8)" % (OBC[L], PBC[L], L))
            continue
        system_check(rep, "fig_decoupling (b) %s U=8 F1" % legend, cand[0][0], L, "OBC",
                     OBC, PBC, ndp=3, where="paper/figs/fig_decoupling_native.tex %s" % legend)
        # y is log2 of the support carried by data/ladder_vs_chain.json
        row = lvc["results"].get("chain-L%d" % L)
        if row:
            # the rebuilt panel plots |S| itself on a log axis, not log2|S|
            rep.eq("fig_decoupling (b) %s U=8 |S|" % legend, cand[0][1],
                   float(row["Sdet"]), 0.51,
                   "the plotted support does not match the deposited chain support %d"
                   % row["Sdet"])

    # 4.6 paper/figs/fig_scaling2_native.tex -- panel (a) is the PERIODIC ring.
    #     H3: 9.3851 here is CORRECT.  "Fixing" it to 9.6047 must fail.
    fsc = lines(root, "paper/figs/fig_scaling2_native.tex")
    fsc_panels = dict(panel_blocks(fsc))
    plotted_faf = {}
    for ln in fsc_panels.get("a", []):
        for (x, y) in coords_of(ln):
            if x in (8.0, 12.0, 16.0, 20.0, 24.0, 28.0):
                plotted_faf.setdefault(int(x) // 2, []).append(y)
    if not plotted_faf:
        rep.missing("fig_scaling2 panel (a)",
                    "the magic series could not be located -- panel (a) of the scaling figure "
                    "is the slot the 2026-08-17 edit touched, and it must stay PBC")
    for L, vals in sorted(plotted_faf.items()):
        for y in vals:
            if L in PBC:
                system_check(rep, "fig_scaling2 (a) L=%d F1" % L, y, L, "PBC", OBC, PBC, ndp=4,
                             where="paper/figs/fig_scaling2_native.tex panel (a)")
            else:
                pt = [p for p in scal["points"] if p["L"] == L]
                if pt:
                    rep.eq("fig_scaling2 (a) L=%d F1 (vs dataset)" % L, y, pt[0]["FAF"], 5e-5,
                           "the plotted ring magic does not match data/scaling_data.json")

    # 4.7 paper/results.tex -- the chain-vs-ladder caption
    res = read_body(root)                     # was paper/results.tex (v2)
    m = re.search(r"open chain vs\.?\\?\s*ladder:\s*\$([\d.]+)/([\d.]+)\$,\s*\n?\s*"
                  r"\$([\d.]+)/([\d.]+)\$,\s*\$([\d.]+)/([\d.]+)\$", res)
    if m:
        g = [float(x) for x in m.groups()]
        for k, L in enumerate((6, 8, 10)):
            if L in OBC:
                system_check(rep, "results.tex ladder caption chain L=%d F1" % L, g[2 * k], L,
                             "OBC", OBC, PBC, ndp=2,
                             where="paper/results.tex (fig:ladder caption, chain)")
            row = lvc["results"].get("ladder-2x%d" % (L // 2))
            if row:
                rep.eq("results.tex ladder caption ladder L=%d F1" % L, g[2 * k + 1],
                       round(row["F1"], 2), 5e-3,
                       "the printed ladder magic does not match data/ladder_vs_chain.json")
    else:
        rep.missing("results.tex ladder caption F1 pairs",
                    "the 'open chain vs. ladder' magic pairs could not be located in "
                    "paper/results.tex")
    m = re.search(r"\$\\times([\d.]+),\\times([\d.]+),\\times([\d.]+)\$", res)
    if m:
        for k, L in enumerate((6, 8, 10)):
            ch = lvc["results"].get("chain-L%d" % L)
            la = lvc["results"].get("ladder-2x%d" % (L // 2))
            if ch and la:
                rep.eq("results.tex ladder/chain chi ratio L=%d" % L, float(m.group(k + 1)),
                       round(la["chi"] / ch["chi"], 1), 5e-2,
                       "the printed chi ratio does not match the deposited chi values "
                       "(%d/%d)" % (la["chi"], ch["chi"]))
    m = re.search(r"\$(\d+),(\d+),(\d+)\$ at \$12/16/20\$ qubits", res)
    if m:
        for k, L in enumerate((6, 8, 10)):
            ch = lvc["results"].get("chain-L%d" % L)
            if ch:
                rep.eq_int("results.tex chain chi L=%d" % L, int(m.group(k + 1)), ch["chi"],
                           "the printed chain bond dimension does not match "
                           "data/ladder_vs_chain.json")

    # 4.8 data/method_max.json is an OBC L=6 run: its E0 must be the open-chain E0
    mm = read_json(root, "data/method_max.json")
    if 6 in E0_OBC:
        ok = abs(mm["E0"] - E0_OBC[6]) <= 1e-5
        other = abs(mm["E0"] - E0_PBC[6]) <= 1e-5
        if ok and not other:
            rep.emit(PASS, "method_max.json E0 is the OPEN chain L=6 energy",
                     "%.9f" % mm["E0"])
        elif other:
            rep.emit(FAIL, "method_max.json E0",
                     "WRONG SYSTEM: E0=%.9f is the periodic-ring energy; the method figure is "
                     "an open-chain run" % mm["E0"])
        else:
            rep.emit(FAIL, "method_max.json E0",
                     "E0=%.9f matches neither the open chain (%.9f) nor the ring (%.9f)"
                     % (mm["E0"], E0_OBC[6], E0_PBC[6]))

    # -----------------------------------------------------------------------
    rep.head("(5) DERIVED QUANTITIES -- figures against the datasets that generate them")
    # -----------------------------------------------------------------------
    # 5.1 fig_scaling2 panels (b) and (c) against data/scaling_data.json
    by_L = {p["L"]: p for p in scal["points"]}
    # Panel (b) now carries TWO curves.  name path=Bpub is the infinite-shot limit,
    # which is the fraction data/scaling_data.json deposits; name path=Bhon is the
    # measured finite-shot curve, which is a different quantity and must NOT be compared
    # with it (they differ by up to 0.07 and comparing them was the failure that first
    # exposed the rebuild).  The two are told apart by the pgfplots path name, not by
    # position in the file.
    hon = read_json(root, "data/honest_sampling.json")
    hon_by_L = {int(k): v for k, v in hon["per_L"].items()}
    n_pub = n_hon = 0
    # In the rebuilt fragment an \addplot spans three physical lines: the options
    # (which carry name path=) and the coordinates are not on the same line.  Reading
    # physical lines found the path names and the coordinates in different places and
    # matched neither -- 0 and 0 points, which the count assertion below caught.
    for ln in logical_lines("\n".join(fsc_panels.get("b", []))):
        is_pub = "name path=Bpub" in ln
        is_hon = "name path=Bhon" in ln
        if not (is_pub or is_hon):
            continue
        for (x, y) in coords_of(ln):
            L = int(x) // 2
            if L not in by_L:
                continue
            if is_pub:
                n_pub += 1
                rep.eq("fig_scaling2 (b) T->inf L=%d fraction" % L, y,
                       round(by_L[L]["frac"], 4), 5e-5,
                       "the plotted infinite-shot sector fraction does not match "
                       "data/scaling_data.json")
            else:
                n_hon += 1
                row = hon_by_L.get(L)
                if row is None:
                    rep.missing("fig_scaling2 (b) finite-T L=%d fraction" % L,
                                "data/honest_sampling.json has no per_L entry for L=%d, so "
                                "the measured curve of panel (b) is unanchored" % L)
                else:
                    rep.eq("fig_scaling2 (b) finite-T L=%d fraction" % L, y,
                           round(float(row["honest_frac_mean"]), 4), 5.1e-4,
                           "the plotted finite-shot fraction does not match "
                           "honest_frac_mean in data/honest_sampling.json")
                    # errbars() returns the sigmas in plotting order, so pair them with
                    # the coordinates by position -- the same order coords_of() returns.
                    _eb = errbars(ln)
                    _cs = coords_of(ln)
                    if len(_eb) == len(_cs):
                        _i = [k for k, c in enumerate(_cs) if abs(c[0] - x) < 1e-6]
                        if _i:
                            rep.eq("fig_scaling2 (b) finite-T L=%d s.d." % L, _eb[_i[0]],
                                   round(float(row["honest_frac_std"]), 4), 5.1e-5,
                                   "the plotted error bar is not honest_frac_std")
    rep.truth("fig_scaling2 (b) plots both curves the caption describes",
              n_pub >= 4 and n_hon >= 4,
              "%d infinite-shot and %d finite-shot points checked" % (n_pub, n_hon),
              "panel (b) should carry an infinite-shot series (name path=Bpub) and a "
              "measured finite-shot series (name path=Bhon); found %d and %d points"
              % (n_pub, n_hon), n=max(n_pub + n_hon, 1))
    # Panel (c) carries three series in v3 (T, dim, |S|); the old code guessed which
    # was which from the first y value, which silently mislabels the new shot-budget
    # curve.  The legend entry that FOLLOWS each \addplot names it.
    _c_lines = fsc_panels.get("c", [])
    for _i, ln in enumerate(_c_lines):                        # panel (c): T, dim and |S|
        cs = coords_of(ln)
        if not cs:
            continue
        _tail = " ".join(_c_lines[_i:_i + 3])
        if BS + "addlegendentry{$T$}" in _tail:
            continue                     # the shot budget: checked in section (8), not here
        is_dim = (BS + "addlegendentry{$" + BS + "dim$}" in _tail)
        for (x, y) in cs:
            L = int(x) // 2
            if L not in by_L:
                continue
            n = L // 2
            if is_dim:
                rep.eq_int("fig_scaling2 (c) L=%d sector dim" % L, int(round(y)),
                           comb(L, n + 1) * comb(L, n),
                           "the plotted sector dimension is not C(L,N/2+1)*C(L,N/2)")
            else:
                want = by_L[L]["frac"] * by_L[L]["Np1_sector"]
                key = "fig_scaling2.support.L%d" % L
                if key in KNOWN_OPEN:
                    rep.known(key, "fig_scaling2 (c) L=%d support" % L, y)
                else:
                    rep.eq("fig_scaling2 (c) L=%d support" % L, y, round(want), 0.5,
                           "the plotted absolute support is not frac*dim from "
                           "data/scaling_data.json (= %.4f)" % want)

    # 5.2 data/ladder_vs_chain.json internal consistency
    for key, row in sorted(lvc["results"].items()):
        L, N = row["L"], row["N"]
        nup = N // 2
        if key.startswith("chain"):
            rep.eq_int("ladder_vs_chain %s sector dim D" % key, row["D"], comb(L, nup) ** 2,
                       "D is not C(L,N/2)^2")
        if "nSector" in row:
            rep.eq_int("ladder_vs_chain %s (N+1) sector" % key, row["nSector"],
                       comb(L, nup + 1) * comb(L, nup),
                       "nSector is not C(L,N/2+1)*C(L,N/2)")
        if "frac" in row and "kR" in row and "nSector" in row:
            rep.eq("ladder_vs_chain %s frac == kR/nSector" % key, row["frac"],
                   round(row["kR"] / row["nSector"], 4), 5e-5,
                   "the reported response fraction is not kR/nSector")

    # 5.3 the two scatter figures against data/resource_master.json, point by point
    pts = rmj["points"]
    fam_of = {"rmMol": "molecule", "rmHub": "hubbard", "rmChn": "chain2", "rmLad": "ladder"}
    n_fam = {}          # (panel, family) -> list of the coordinate lists that draw it
    for panel, xkey in (("a", "chi"), ("b", "F1")):
        for ln in frm_panels.get(panel, []):
            tag = [k for k in fam_of if k in ln and "coordinates" in ln]
            if not tag:
                continue
            fam = fam_of[tag[0]]
            cs = coords_of(ln)
            family_pts = [p for p in pts if p["family"] == fam]
            unmatched = []
            for (x, y) in cs:
                hit = [p for p in family_pts
                       if abs(p["S"] - y) < 0.51 and abs(p[xkey] - x) < 6e-4]
                if not hit:
                    unmatched.append((x, y))
            rep.truth("fig_resource_master (%s) %s series (%s vs |S|), %d points"
                      % (panel, fam, xkey, len(cs)), not unmatched,
                      "every plotted point is a deposited point",
                      "%d plotted point(s) have no counterpart in data/resource_master.json: %s"
                      % (len(unmatched), unmatched[:4]), n=2 * len(cs))
            # The rebuilt figure draws the molecular family as TWO series -- the 27 with
            # |S| > 1 as filled marks and the 11 pinned at the floor |S| = 1 as open ones,
            # each labelled with its own count -- so the count assertion is on the union
            # and on the split, not on one series.  It is now stronger than it was: the
            # split itself has to be right.
            n_fam.setdefault((panel, fam), []).append(cs)

    for (panel, fam), groups in sorted(n_fam.items()):
        allpts = [c for g in groups for c in g]
        family_pts = [p for p in pts if p["family"] == fam]
        rep.eq_int("fig_resource_master (%s) %s point count (union of %d series)"
                   % (panel, fam, len(groups)), len(allpts), len(family_pts),
                   "the figure plots a different number of points than the dataset holds")
        if fam == "molecule" and len(groups) == 2:
            floor = sorted(len([c for c in g if abs(c[1] - 1.0) < 0.51]) for g in groups)
            due = sorted((len([p for p in family_pts if p["S"] <= 1]),
                          0))
            rep.eq_int("fig_resource_master (%s) molecules at the |S|=1 floor" % panel,
                       max(floor), max(due),
                       "the open-marker series is supposed to be exactly the molecules "
                       "pinned at |S| = 1 in data/resource_master.json")

    mol_pts = [p for p in pts if p["family"] == "molecule"]
    hub_pts = [p for p in pts if p["family"] == "hubbard"]
    dec_mol, dec_hub = [], []
    for ln in fdc:
        if "coordinates" not in ln:
            continue
        if any(t in ln for t in ("dcCov", "dcMrf", "dcIon")):
            dec_mol += coords_of(ln)
        elif any(t in ln for t in ("dcH1", "dcH2", "dcH3")):
            dec_hub += coords_of(ln)
    # v3 draws the six Hubbard strata twice -- once against F1 and once against chi -- so
    # a style-name sweep would feed the chi panel into an F1 comparison.  The panel-(b)
    # series were already identified above by matching the deposited (FAF, S) rows, and
    # those are the ones that belong here.
    if not dec_hub:
        dec_hub = [c for _cs in _dec_b.values() for c in _cs]
    for label, cs, src in (("molecules", dec_mol, mol_pts), ("Hubbard", dec_hub, hub_pts)):
        # the rebuild plots |S| on a log axis; v2 plotted log2|S| on a linear one
        unmatched = [(x, y) for (x, y) in cs
                     if not [p for p in src
                             if abs(p["F1"] - x) < 6e-4
                             and abs(max(p["S"], 1) - y) < 0.51]]
        rep.truth("fig_decoupling %s series (F1 vs |S|), %d points" % (label, len(cs)),
                  not unmatched, "every plotted point is a deposited point",
                  "%d plotted point(s) have no counterpart in data/resource_master.json: %s"
                  % (len(unmatched), unmatched[:4]), n=2 * len(cs))
    rep.eq_int("fig_decoupling molecule point count", len(dec_mol), len(mol_pts),
               "the molecular panel plots a different number of points than the dataset holds")

    # 5.4 fig_decoupling panel (c): the free-fermion labels are log2 of the plotted y
    labels = [int(m.group(1)) for m in re.finditer(r"inner sep=1\.5pt\] at \(axis cs:\d+,[\d.]+\)\{\$(\d+)\$\}",
                                                   "\n".join(fdc))]
    i = find_line(fdc, "dcFree")
    if i >= 0 and labels:
        cs = coords_of(fdc[i])
        for (x, y), nlab in zip(cs, labels):
            rep.eq("fig_decoupling (c) L=%d log2|S| label" % int(x), y,
                   round(math.log2(nlab), 3), 5e-4,
                   "the plotted ordinate is not log2 of the printed support %d" % nlab)

    # 5.5 the molecular trend line and the annotated Pearson r
    mF = np.array([p["F1"] for p in mol_pts])
    mS = np.array([max(p["S"], 1) for p in mol_pts], float)
    slope, intercept = np.polyfit(mF, np.log2(mS), 1)
    m = re.search(r"\{([\d.]+)\*x\+([\d.]+)\}", "\n".join(fdc))
    if m:
        rep.eq("fig_decoupling (a) trend slope", float(m.group(1)), round(slope, 4), 5e-5,
               "the drawn trend line is not the least-squares fit of log2|S| on F1 over the "
               "deposited molecular points")
        rep.eq("fig_decoupling (a) trend intercept", float(m.group(2)), round(intercept, 4), 5e-5,
               "the drawn trend line is not the least-squares fit over the deposited points")
    m = re.search(r"\{\$r=([\d.]+)\$\}", "\n".join(fdc))
    if m:
        rep.eq("fig_decoupling (a) annotated Pearson r", float(m.group(1)),
               round(pearsonr(mF, np.log2(mS))[0], 2), 5e-3,
               "the annotated r is not the Pearson correlation of F1 with log2|S| over the "
               "deposited molecular points")

    # -----------------------------------------------------------------------
    rep.head("(6) STATISTICS -- every published correlation, recomputed from the points")
    # -----------------------------------------------------------------------
    aF = np.array([p["F1"] for p in pts])
    aS = np.array([p["S"] for p in pts], float)
    aX = np.array([p["chi"] for p in pts], float)
    rep.eq_int("resource_master.json point count", len(pts), 74,
               "the master dataset no longer holds the 74 states the paper reports")
    for fam, want in (("molecule", 38), ("hubbard", 30), ("chain2", 3), ("ladder", 3)):
        rep.eq_int("resource_master.json family '%s' count" % fam,
                   sum(1 for p in pts if p["family"] == fam), want,
                   "the family composition of the master dataset changed")
    rep.eq("resource_master.json pooled Spearman chi-|S|", rmj["pooled_spearman_chi_S"],
           spearmanr(aX, aS).correlation, 1e-12,
           "the stored pooled Spearman is not reproducible from the deposited points")
    rep.known("resource_master.pooled_spearman_F1_S",
              "resource_master.json pooled Spearman F1-|S|", rmj["pooled_spearman_F1_S"])

    rsrc = re.sub(r"\s+", " ", read_body(root))  # was paper/resource.tex (v2)
    mX = np.array([p["chi"] for p in mol_pts], float)
    hF = np.array([p["F1"] for p in hub_pts])
    hS = np.array([p["S"] for p in hub_pts], float)
    hX = np.array([p["chi"] for p in hub_pts], float)
    # v3 prints the same six coefficients, in different sentences and in different
    # files.  Four of them now carry the word "withdrawn" beside them -- which does NOT
    # exempt them from being right: a retracted number quoted wrongly is still a wrong
    # number, and a reader checking the retraction is the person most likely to look it up.
    printed = [
        ("chi-|S| over n=38 molecules", spearmanr(mX, mS).correlation,
         r"\\rho=(-?[\d.]+)\$ within the \$n=38\$ molecular suite"),
        ("F1-|S| over n=38 molecules", spearmanr(mF, mS).correlation,
         r"comparable\s+\$\\rho=(-?[\d.]+)\$"),
        ("chi-|S| over n=30 Hubbard (withdrawn, still quoted)",
         spearmanr(hX, hS).correlation,
         r"over these same \$30\$ points: \$\\rho=(-?[\d.]+)\$"),
        ("F1-|S| over n=30 Hubbard (withdrawn, still quoted)",
         spearmanr(hF, hS).correlation,
         r"and \$\\rho=(-?[\d.]+)\$ for\s+\$\\mathcal F_1\$"),
        ("pooled chi-|S| over the 74", spearmanr(aX, aS).correlation,
         r"pooled Spearman\s+\$\\rho=(-?[\d.]+)\$"),
    ]
    # The sixth, the pooled F1-|S| over the 74, is no longer printed anywhere: v3 retracts
    # it instead of quoting it.  Assert that, so it cannot come back unnoticed.
    _f1_74 = re.search(r"Spearman \$\\rho=(-?[\d.]+)\$ for \$\\mathcal F_1\$ vs", rsrc)
    rep.truth("pooled F1-|S| over the 74 is no longer asserted in the body",
              _f1_74 is None,
              "the sentence that offered the pooled F1 coefficient as evidence is gone, "
              "which is what Sec. 6 says it did",
              "the retracted pooled F1-|S| coefficient is printed again (%s); it is a "
              "Simpson artefact and the manuscript retracts it"
              % (_f1_74.group(1) if _f1_74 else ""))
    for label, value, pattern in printed:
        m = re.search(pattern, rsrc)
        if m is None:
            rep.missing("resource.tex printed rho for %s" % label,
                        "the sentence carrying this correlation could not be located; the "
                        "guardian cannot protect a number it cannot find (pattern %r)" % pattern)
        else:
            rep.eq("resource.tex printed rho, %s" % label, float(m.group(1)), round(value, 2),
                   5e-3,
                   "the printed Spearman rho disagrees with the value recomputed from the "
                   "deposited points (%.6f)" % value)
    m = re.search(r"\$7\$ of the \$38\$ molecular points", rsrc)
    nviol = sum(1 for p in mol_pts if p["chi"] > p["S"])
    if m:
        rep.eq_int("resource.tex '7 of the 38' chi > |S|_eps violations", 7, nviol,
                   "the confessed number of points where chi exceeds the thresholded support "
                   "does not match the deposited points")
    else:
        rep.missing("resource.tex '7 of the 38'",
                    "the confession that chi > |S|_eps in 7 of 38 points is gone; recomputed "
                    "count is %d" % nviol)
    m = re.search(r"differing by up to \$\\sim\\!([\d.]+)\\times\$", rsrc)
    if m:
        ratios = [lvc["results"]["ladder-2x%d" % (L // 2)]["Sdet"]
                  / lvc["results"]["chain-L%d" % L]["Sdet"] for L in (6, 8, 10)]
        rep.eq("resource.tex chain/ladder |S| ratio", float(m.group(1)),
               round(max(ratios), 1), 5e-2,
               "the printed maximum chain-to-ladder support ratio does not match "
               "data/ladder_vs_chain.json (%.4f)" % max(ratios))
    m = re.search(r"\$p\\approx([\d.]+)\$", rsrc)
    if m:
        rep.eq("resource.tex printed p-value for the pooled F1-|S| null", float(m.group(1)),
               spearmanr(hF, hS).pvalue, 1e-2,
               "the printed p disagrees with the recomputed two-sided Spearman p (%.4f). "
               "NOTE: this pooled test is a Simpson's-paradox artefact -- within each of the "
               "six (L,N) strata the F1-|S| rank correlation is exactly -1"
               % spearmanr(hF, hS).pvalue)
    # the stratified statement the pooled null hides -- computed, not asserted
    ce = read_json(root, "data/cost_vs_ent.json")
    rows = ce["rows"]
    rep.eq("cost_vs_ent.json r_FAF_S vs its own rows", ce["r_FAF_S"],
           spearmanr([r["FAF"] for r in rows], [r["S"] for r in rows]).correlation, 1e-12,
           "the stored correlation is not reproducible from the rows in the same file")
    rep.eq("cost_vs_ent.json r_chi_S vs its own rows", ce["r_chi_S"],
           spearmanr([r["chi"] for r in rows], [r["S"] for r in rows]).correlation, 1e-12,
           "the stored correlation is not reproducible from the rows in the same file")
    strata = {}
    for r in rows:
        strata.setdefault((r["L"], r["N"]), []).append(r)
    rhos = []
    for k, rs in sorted(strata.items()):
        rhos.append(spearmanr([r["FAF"] for r in rs], [r["S"] for r in rs]).correlation)
    rep.truth("cost_vs_ent.json: F1-|S| is exactly -1 within every (L,N) stratum",
              len(rhos) == 6 and all(abs(r + 1.0) < 1e-12 for r in rhos),
              "six strata, all rho = -1.000000",
              "the stratified structure changed: %s" % [round(r, 4) for r in rhos],
              n=len(rhos))

    # -----------------------------------------------------------------------
    rep.head("(7) MOLECULAR SUITE -- both printed tables against the deposited data")
    # -----------------------------------------------------------------------
    suite = read_json(root, "data/n19_suite.json")["mols"]
    spec = {r["mol"]: r for r in read_json(root, "data/n19_spectral.json")["mols"]}
    by_mol = {r["mol"]: r for r in suite}
    rep.eq_int("n19_suite.json molecule count", len(suite), 19,
               "the suite no longer holds nineteen molecules")

    def clean(name):
        s = re.sub(r"\$\^?\{?\\dagger\}?\$?", "", name)
        s = s.replace("$_", "").replace("$", "").replace("{", "").replace("}", "")
        s = s.replace("\\", "").replace("_", "").strip()
        return s

    # 7.1 tab:molecules  (paper/table_molecules.tex)
    tm = lines(root, "paper/table_molecules_v3.tex")
    nrows = 0
    for ln in tm:
        if "&" not in ln or "multicolumn" in ln or "Molecule &" in ln:
            continue
        cells = [c.strip() for c in ln.split("\\\\")[0].split("&")]
        if len(cells) != 8:
            continue
        mol = clean(cells[0])
        r = by_mol.get(mol)
        if r is None:
            rep.missing("table_molecules row %r" % mol,
                        "printed row has no counterpart in data/n19_suite.json")
            continue
        nrows += 1
        mcas = re.search(r"\((\d+)\+(\d+),(\d+)\)", cells[2])
        if mcas:
            na, nb, no = (int(x) for x in mcas.groups())
            rep.eq_int("table_molecules %s CAS n_e" % mol, (na, nb), tuple(r["nelec"]),
                       "printed active-space electron count disagrees with n19_suite.json")
            rep.eq_int("table_molecules %s CAS n_o" % mol, no, r["ncas"],
                       "printed active-space orbital count disagrees with n19_suite.json")
        rep.eq("table_molecules %s F1_eq" % mol, float(cells[3]), round(r["FAF_eq"], 3), 5e-4,
               "printed equilibrium magic disagrees with n19_suite.json (%.6f)" % r["FAF_eq"])
        rep.eq("table_molecules %s F1_diss" % mol, float(cells[4]), round(r["FAF_diss"], 3), 5e-4,
               "printed dissociation magic disagrees with n19_suite.json (%.6f)" % r["FAF_diss"])
        rep.eq("table_molecules %s F1_diss/F1_max" % mol, float(cells[5]),
               round(r["FAF_diss"] / r["fafmax"], 2), 5e-3,
               "printed normalised magic is not FAF_diss/fafmax from n19_suite.json")
        sp_r = spec.get(mol)
        if sp_r is not None:
            cell = cells[6]
            mm2 = re.search(r"<\\!10\^\{(-?\d+)\}", cell)
            if mm2:
                thr = 10.0 ** int(mm2.group(1))
                rep.truth("table_molecules %s rel-L1 < 1e%s" % (mol, mm2.group(1)),
                          sp_r["relL1_final"] < thr,
                          "%.3g < %g" % (sp_r["relL1_final"], thr),
                          "printed bound is violated: n19_spectral.json gives %.6g"
                          % sp_r["relL1_final"])
            else:
                mm3 = re.search(r"(\d+)\\times10\^\{(-?\d+)\}", cell)
                if mm3:
                    want = float(mm3.group(1)) * 10.0 ** int(mm3.group(2))
                    got = sp_r["relL1_final"]
                    exp = math.floor(math.log10(got)) if got > 0 else 0
                    rounded = round(got / 10 ** exp) * 10 ** exp
                    rep.eq("table_molecules %s rel-L1" % mol, want, rounded,
                           abs(want) * 0.5 + 1e-30,
                           "printed rel-L1 does not round to the deposited %.4g" % got)
            rep.eq("table_molecules %s |S|/D" % mol, float(cells[7]),
                   round(sp_r["frac_final"], 2), 5e-3,
                   "printed sector fraction disagrees with n19_spectral.json (%.6f)"
                   % sp_r["frac_final"])
    rep.eq_int("table_molecules row count", nrows, 19,
               "the printed suite table no longer has nineteen data rows")

    # 7.2 tab:repro (inside paper/theorem_b2.tex), against n19_suite.json and the
    #     geometry table parsed out of src/n19_suite.py (parsed with ast, not imported:
    #     importing it would pull in pyscf).
    geom = {}
    try:
        tree = ast.parse(read(root, "src/n19_suite.py"))
        mols_node = None
        for node in tree.body:
            if isinstance(node, ast.Assign) and any(
                    isinstance(t, ast.Name) and t.id == "MOLS" for t in node.targets):
                mols_node = node.value
        def atoms_of(node):
            if isinstance(node, ast.Constant):
                return node.value
            if isinstance(node, ast.Call) and getattr(node.func, "id", "") == "diat":
                a, b, R = (ast.literal_eval(x) for x in node.args)
                return "%s 0 0 0; %s 0 0 %s" % (a, b, R)
            raise ValueError("unrecognised geometry entry")
        def first_bond(atoms):
            parts = [p.split() for p in atoms.split(";")]
            p0 = np.array([float(x) for x in parts[0][1:4]])
            p1 = np.array([float(x) for x in parts[1][1:4]])
            return float(np.linalg.norm(p1 - p0))
        for elt in mols_node.elts:
            label = ast.literal_eval(elt.elts[0])
            geom[label] = (first_bond(atoms_of(elt.elts[7])),
                           first_bond(atoms_of(elt.elts[8])))
    except Exception as exc:                                    # noqa: BLE001
        rep.skip("src/n19_suite.py geometry table",
                 "could not parse the MOLS table (%s); R_eq/R_diss not checked" % exc)

    trow = 0
    for ln in lines(root, "paper/app_carried_repro.tex"):
        if "&" not in ln or "multicolumn" in ln or "Molecule &" in ln:
            continue
        cells = [c.strip() for c in ln.split("\\\\")[0].split("&")]
        if len(cells) != 9:
            continue
        mol = clean(cells[0])
        r = by_mol.get(mol)
        if r is None:
            continue
        trow += 1
        if mol in geom:
            key = "tab_repro.%s.R_eq" % mol
            if key in KNOWN_OPEN:
                rep.known(key, "tab_repro %s R_eq" % mol, float(cells[2]))
            else:
                rep.eq("tab_repro %s R_eq" % mol, float(cells[2]), round(geom[mol][0], 3), 5e-4,
                       "printed equilibrium bond length is not the geometry in "
                       "src/n19_suite.py (%.5f A)" % geom[mol][0])
            rep.eq("tab_repro %s R_diss" % mol, float(cells[3]), round(geom[mol][1], 2), 5e-3,
                   "printed dissociation bond length is not the geometry in src/n19_suite.py "
                   "(%.5f A)" % geom[mol][1])
        na, nb = r["nelec"]
        rep.eq_int("tab_repro %s sector dim D" % mol, int(cells[4]),
                   comb(r["ncas"], na) * comb(r["ncas"], nb),
                   "printed sector dimension is not C(ncas,na)*C(ncas,nb)")
        rep.eq_int("tab_repro %s |S|_eq" % mol, int(cells[5]), r["n99_eq"],
                   "printed equilibrium support disagrees with n19_suite.json")
        rep.eq_int("tab_repro %s |S|_diss" % mol, int(cells[6]), r["n99_diss"],
                   "printed dissociation support disagrees with n19_suite.json")
        rep.eq_int("tab_repro %s chi_eq" % mol, int(cells[7]), r["chi_eq"],
                   "printed equilibrium bond dimension disagrees with n19_suite.json")
        rep.eq_int("tab_repro %s chi_diss" % mol, int(cells[8]), r["chi_diss"],
                   "printed dissociation bond dimension disagrees with n19_suite.json")
    rep.eq_int("tab_repro row count", trow, 19,
               "the reproducibility table no longer has nineteen data rows")

    # -----------------------------------------------------------------------
    rep.head("(8) SPECTRAL ARTEFACTS -- provenance of the .dat tables and the sum rule")
    # -----------------------------------------------------------------------
    def edges_from(jsonrel, thrfrac=0.03):
        """Reconstruct the edge table from the deposited spectra with the recipe of
        src/spin_lanczos.py:35-41."""
        d = read_json(root, jsonrel)
        wg = np.array(d["wg"])
        qs = sorted(float(k) for k in d["S"])
        M = np.array([d["S"]["%.5f" % q] for q in qs])
        out = {}
        for i, q in enumerate(qs):
            A = M[i]
            thr = thrfrac * M.max()
            sig = wg[A > thr]
            lo = sig.min() if len(sig) else 0.0
            hi = sig.max() if len(sig) else 0.0
            pk = wg[int(np.argmax(A))]
            out[round(q, 4)] = (float("%.3f" % lo), float("%.3f" % hi), float("%.3f" % pk))
        return out, np.array(d["wg"]), d

    for datrel, jsonrel, label in (("paper/figs/spinqw_edges.dat", "data/spinqw_L12.json",
                                    "S^zz(q,w) edges"),
                                   ("paper/figs/sqw_edges.dat", "data/sqw_L12.json",
                                    "S(q,w) edges")):
        try:
            recon, wg, d = edges_from(jsonrel)
        except Exception as exc:                                # noqa: BLE001
            rep.missing(label, "cannot reconstruct from %s: %s" % (jsonrel, exc))
            continue
        rows = [l.split() for l in read(root, datrel).strip().split("\n")[1:] if l.strip()]
        orphans = []
        for cells in rows:
            q = round(float(cells[0]), 4)
            vals = tuple(float(x) for x in cells[1:4])
            if q not in recon:
                orphans.append((q, vals))
                continue
            ok = all(abs(a - b) < 5e-4 for a, b in zip(vals, recon[q]))
            rep.truth("%s q/pi=%.4f" % (label, q), ok,
                      "row reproduces from %s" % jsonrel,
                      "row does NOT reproduce from %s: file has %s, the deposited spectrum "
                      "gives %s" % (jsonrel, vals, recon[q]), n=3)
        if orphans:
            rep.emit(FAIL, "%s rows with no computational origin" % label,
                     "FABRICATED ROWS: %d row(s) of %s carry momenta that do not exist in "
                     "%s: %s" % (len(orphans), datrel, jsonrel, orphans),
                     n=3 * len(orphans))
        else:
            rep.emit(PASS, "%s: every row has a computational origin" % label, "", n=0)

    # eta consistency between the deposited spectrum and its generator's default
    for jsonrel, srcrel, var in (("data/sqw_L12.json", "src/sqw_lanczos.py", "eta"),
                                 ("data/spinqw_L12.json", "src/spin_lanczos.py", "eta")):
        d = read_json(root, jsonrel)
        # strip comments first: a '# eta = 0.18' note above the code would otherwise be
        # read as the default and would silently repair this check (attack E2).
        _srctxt = re.sub(r"#[^" + chr(10) + "]*", "", read(root, srcrel))
        m = re.search(r"(?<![A-Za-z_])eta\s*=\s*([\d.]+)", _srctxt)
        if m is None:
            rep.missing("%s eta default" % os.path.basename(srcrel), "no eta default found")
            continue
        name = "%s eta matches %s default" % (os.path.basename(jsonrel),
                                              os.path.basename(srcrel))
        key = "sqw_lanczos.eta_default"
        if "sqw_lanczos" in srcrel and key in KNOWN_OPEN:
            rep.known(key, name, float(m.group(1)))
        else:
            rep.truth(name, abs(float(m.group(1)) - d["eta"]) < 1e-12,
                      "eta=%g" % d["eta"],
                      "the deposited spectrum was computed at eta=%g but %s now defaults to "
                      "eta=%s: the deposited file cannot be regenerated by running the script "
                      "as committed" % (d["eta"], srcrel, m.group(1)))

    # the sum rule of the method figure
    grid = np.array(mm["grid"])
    A_exact = np.array(mm["A_exact"])
    A_final = np.array(mm["A_final"])
    I_exact = float(np.trapz(A_exact, grid))
    I_final = float(np.trapz(A_final, grid))
    rep.known("method_max.sumrule", "method_max.json 'sumrule' vs the integral of A_exact",
              mm["sumrule"])
    print("        stored 'sumrule' (an analytic ground-state identity)  : %.8f" % mm["sumrule"])
    print("        integral of the deposited A_exact over the window     : %.8f  "
          "(%.2f%% below the sum rule -- finite window)"
          % (I_exact, 100 * (mm["sumrule"] - I_exact) / mm["sumrule"]))
    print("        integral of the deposited reconstruction A_final      : %.8f  "
          "(%.2f%% below the sum rule)"
          % (I_final, 100 * (mm["sumrule"] - I_final) / mm["sumrule"]))

    # hardware: "exact by coverage" must be an identity, not a fidelity claim
    her = read_json(root, "data/heron_spectral.json")
    hg = np.array(her["grid"])
    hw = np.array(her["A_hw"])
    hx = np.array(her["A_exact"])
    rel = float(np.trapz(np.abs(hw - hx), hg) / np.trapz(hx, hg))
    rep.eq("heron_spectral.json hw_relL1 reproduces from the deposited arrays",
           her["hw_relL1"], rel, 1e-9,
           "the stored hardware relative-L1 is not the integral of |A_hw - A_exact| over the "
           "deposited grid (%.3g)" % rel)
    rep.truth("heron_spectral.json: hw_relL1 ~ 0 only because the subspace is the whole sector",
              (rel < 1e-9) == (her["hw_S"] == her["nsector"]),
              "hw_S = %d of nsector = %d -- exact BY COVERAGE, not by fidelity"
              % (her["hw_S"], her["nsector"]),
              "the coverage identity broke: hw_relL1=%.3g with hw_S=%d of %d. If the retained "
              "subspace is no longer the full sector, a zero error is a fidelity claim and must "
              "be defended as one" % (rel, her["hw_S"], her["nsector"]))
    bo = read_json(root, "data/hw_bitorder_check.json")
    rep.truth("hw_bitorder_check.json: the notebook's reversed bit order retains none of the intended shots",
              bo["kept_notebook_reading"] == 0 and bo["kept_correct_reading"] == bo["total_shots"] == 350000
              and bo["support_correct_reading"] == bo["sector_dimension"] == her["nsector"],
              "0 of %d noiseless shots kept under int(bs[::-1],2), all %d (support %d/%d) under int(bs,2): "
              "every determinant retained from the device is a device error"
              % (bo["total_shots"], bo["kept_correct_reading"], bo["support_correct_reading"],
                 bo["sector_dimension"]),
              "the deposited bit-order check no longer says what Sec. IX A and Sec. S8 print: %r" % bo)

    # -----------------------------------------------------------------------
    section9(rep, root, dict(obc=OBC, pbc=PBC, mF=mF, mS=mS, mX=mX, hF=hF, hS=hS, hX=hX,
                             mol_pts=mol_pts, hub_pts=hub_pts, lvc=lvc, mm=mm, her=her))
    section10(rep, root)
    section_c3(rep, root)
    if not args.fast:
        coverage_floor(rep)
    else:
        rep.skip("coverage floor", "--fast cannot meet the floors; it certifies nothing")

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    dt = time.time() - t0
    print("\n" + "=" * 78)
    print("COVERAGE  checks / numeric assertions, by section")
    for title, nnum, nchk in rep.per_section:
        print("   %-70s %4d / %4d" % (title[:70], nchk, nnum))
    print("=" * 78)
    print("SUMMARY   pass %d   FAIL %d   xfail %d   xpass %d   skip %d   "
          "(%d numeric assertions, %.1f s)"
          % (rep.counts[PASS], rep.counts[FAIL], rep.counts[XFAIL], rep.counts[XPASS],
             rep.counts[SKIP], rep.numbers, dt))
    print("=" * 78)

    if rep.xfails:
        print("\nKNOWN OPEN DEFECTS (registered in KNOWN_OPEN; they do not set the exit code):")
        for sec, name, msg in rep.xfails:
            print("  - %s\n      %s" % (name, msg))
    if rep.xpasses:
        print("\nREPAIRED (delete the KNOWN_OPEN entry to turn these into hard checks):")
        for sec, name, msg in rep.xpasses:
            print("  - %s\n      %s" % (name, msg))
    if rep.skips:
        print("\nNOT CHECKED:")
        for sec, name, msg in rep.skips:
            print("  - %s: %s" % (name, msg))

    if args.fast:
        print("\nRESULT: INCOMPLETE -- --fast skipped the L=10 diagonalisations, the checks that "
              "depend on them and the coverage floor.\n"
              "        It CERTIFIES NOTHING and exits 3 by construction, pass or fail. An "
              "adversarial pass reintroduced the\n"
              "        2026-08-17 open-chain/periodic-ring defect at L=10 and watched --fast "
              "exit 0. Run without --fast.")
        if rep.failures:
            print("\nFAILURES (found even in fast mode):")
            for sec, name, msg in rep.failures:
                print("  * [%s] %s\n      %s" % (sec, name, msg))
        return 3

    if rep.failures:
        print("\nFAILURES:")
        for sec, name, msg in rep.failures:
            print("  * [%s] %s\n      %s" % (sec, name, msg))
        print("\nRESULT: FAIL -- %d check(s) failed. Each line above names the defect."
              % len(rep.failures))
        return 1

    print("\nRESULT: PASS -- every recomputed quantity is carried correctly by the deposited "
          "artefacts.")
    if rep.xfails:
        print("        %d known open defect(s) remain, listed above." % len(rep.xfails))
    return 0


def guarded_main(argv=None):
    """A guardian that crashes must still fail loudly and say what it could not read."""
    try:
        return main(argv)
    except Exception as exc:                                    # noqa: BLE001
        import traceback
        traceback.print_exc()
        print("\nRESULT: FAIL -- the guardian could not complete: %s: %s\n"
              "        A missing or malformed artefact is itself a defect; do not treat an "
              "aborted run as a pass." % (type(exc).__name__, exc))
        return 2


if __name__ == "__main__":
    sys.exit(guarded_main())
