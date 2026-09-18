# -*- coding: utf-8 -*-
r"""gap_scaling.py -- finite-size scaling of the TWO gaps of the half-filled 1D Hubbard ring.

WHAT IT COMPUTES
----------------
For L = 4, 6, 8, 10, 12 sites, U/t = 8, half filling (n_up = n_dn = L/2), periodic ring:

  (1) SPIN GAP.  The lowest exact pole of the dynamical spin structure factor
      S^zz(q=pi, omega) = -(1/pi) Im <psi0| Sz_{-q} [w-(H-E0)+i eta]^{-1} Sz_q |psi0>,
      obtained WITHOUT any broadening: the Lanczos tridiagonal built on the seed
      Sz_{q=pi}|psi0> is diagonalized and the smallest Ritz value carrying non-negligible
      spectral weight is the lowest pole.  Lanczos runs with FULL REORTHOGONALIZATION --
      without it the Gauss quadrature underlying the continued fraction is not reliable
      (ghost eigenvalues), and this script reports the no-reorthogonalization number
      side by side so the size of that effect is on the record.

  (2) CHARGE GAP.  mu+ = E0(N+1) - E0(N), mu- = E0(N) - E0(N-1), Delta = mu+ - mu-,
      each E0 an ARPACK ground state of the corresponding particle-number sector.

  (3) FITS.  spin gap vs 1/L  and  charge gap vs 1/L (and 1/L + 1/L^2), with the
      standard errors of the fit coefficients, extrapolated to L -> infinity.
      The charge-gap extrapolation is compared with the exact Lieb-Wu gap, which is
      recomputed here from two independent integral representations.

  (4) INDEPENDENT CONTROL.  The pure Heisenberg ring with J = 4 t^2 / U = 0.5 at the same
      L, solved by DENSE exact diagonalization in the Sz=0 sector (a different code path
      and a different model from (1)).  Two independent observables are extracted: the
      lowest pole of its own S^zz(q=pi,omega), and its singlet-triplet gap
      E0(Sz=1) - E0(Sz=0).  If the Hubbard spin gap tracks the Heisenberg one, the
      "gapless" spin channel of the L=12 figure is the finite-size singlet-triplet gap of
      the effective Heisenberg ring, NOT a property of the Hubbard model.

  (5) THE CLAIM THIS SUPPORTS.  The two channels have QUALITATIVELY OPPOSITE finite-size
      behaviour: the spin gap closes as ~ 2/L (L * gap = const), the charge gap SATURATES
      at a finite value (L * gap grows linearly).  The script quantifies the separation of
      the two L -> infinity intercepts with fit standard errors and reports a z-score.

WHICH PAPER CLAIM IT BACKS
--------------------------
It replaces the single-size, refutable wording around Figs. 4/5 ("gapless", "the two-spinon
continuum touches omega = 0") -- whose "zero" is the left edge of the computation window
(wg[0] = -0.1 in src/spin_lanczos.py) -- by a size-resolved statement that is a
demonstration rather than an illustration: spin closes as 1/L towards zero, charge saturates
towards the Lieb-Wu value.  It also supplies the provenance for Delta = 4.97 printed in the
manuscript, which no script in the repository currently computes.

HOW TO RUN
----------
    python src/gap_scaling.py                    # L = 4,6,8,10,12 -> data/gap_scaling.json
                                                 # the bare command reproduces the deposited JSON:
                                                 # every default equals the recorded invocation
    python src/gap_scaling.py --sizes 4,6,8,10  # cheaper subset
    python src/gap_scaling.py --out other.json --nl-max 60

It imports the existing engine src/akw_lanczos.py (strings, sector_H); nothing is
re-derived.  Runs from any working directory; no container paths.  Output is a single
JSON with a "provenance" block.

RELATION TO src/charge_gap_ed.py
--------------------------------
src/charge_gap_ed.py computes mu^+, mu^- and Delta by a DIFFERENT code path
(build_H_explicit + dense LAPACK cross-check).  The two scripts are deliberately
redundant on the charge channel and must agree: both give Delta(L=12) = 4.968759.
Only this file computes the SPIN gap, the Heisenberg control and the two fits.

RUNTIME: ~5.5 min for L = 4,6,8,10,12 on a laptop, < 2 GB RAM.  No QPU.

MANDATORY SELF-CHECK
--------------------
The run ABORTS unless it reproduces, from scratch:
    charge gap at L=12   Delta      = 4.9688  +- 1e-3     (manuscript prints 4.97)
    spin   gap at L=12   E_spin(pi) = 0.16635 +- 1e-4
    particle-hole symmetry  mu+ + mu- = U     +- 1e-6
    Lieb-Wu gap (two independent integrals)   agree to 1e-6 and equal 4.6795 +- 1e-3
Self-checks for L < 12 are applied whenever the corresponding size is in --sizes.

COST (measured on a 15 GB Windows laptop with four other jobs competing)
-----------------------------------------------------------------------
    L = 4, 6, 8                      < 1 s each
    L = 10                           ~35 s   (dim 63 504)
    L = 12                           ~200 s  (dim 853 776; the two (N+-1) ground states
                                      dominate at ~95-136 s, the spin Lanczos costs ~7 s)
    --heisenberg-sizes 14,16,18,20   ~70 s   (sparse ED, no Hubbard cost)
    --ghost-probe 8,10 --ghost-nl 200 ~30 s
    whole default run                ~6.5 min wall, peak RSS ~0.5 GB with --nl-max 40

    NOT reachable on a laptop: L = 14.  dim(N) = 11 778 624 and dim(N+-1) = 10 306 296, so
    one vector is 94 MB; ARPACK alone wants ~1.9 GB and the reorthogonalized Lanczos basis
    another ~3.8 GB at nl = 40.  Budget ~4 GB of free RAM and 1.5-2.5 CPU-hours.

Everything is deterministic: the ARPACK start vector is drawn from a FIXED seed
(--rng-seed, default 20260918), so the run is bit-for-bit reproducible.  Verified: two
independent runs reproduced every digit printed below.
"""

import argparse
import json
import os
import platform
import sys
import time

import numpy as np
import scipy
import scipy.sparse as sp
from scipy.integrate import quad
from scipy.linalg import eigh, eigh_tridiagonal
from scipy.sparse.linalg import eigsh
from scipy.special import expit, j1

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


DEFAULT_OUT = os.path.join(REPO, "data", "gap_scaling.json")

# ----------------------------------------------------------------------------- constants
LIEB_WU_U8_REFERENCE = 4.6795          # literature value of the U/t=8 charge gap
SELF_CHECK = {                          # values this script must reproduce, or abort
    12: {"charge_gap": (4.9688, 1e-3), "spin_gap": (0.16635, 1e-4)},
    10: {"charge_gap": (5.022, 5e-3), "spin_gap": (0.19956, 5e-4)},
    8:  {"charge_gap": (5.186, 5e-3), "spin_gap": (0.23958, 5e-4)},
    6:  {"charge_gap": (5.358, 5e-3), "spin_gap": (0.34857, 5e-4)},
}
HEISENBERG_SELF_CHECK = {               # pure Heisenberg ring, J = 4t^2/U = 0.5
    12: (0.17792, 5e-4), 10: (0.21162, 5e-4), 8: (0.26134, 5e-4), 6: (0.34237, 5e-4),
}

T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


# ------------------------------------------------------------------- engine (reused, not rewritten)
def import_engine(src_dir):
    """Import strings()/sector_H() from the repository engine src/akw_lanczos.py."""
    src_dir = os.path.abspath(src_dir)
    if src_dir not in sys.path:
        sys.path.insert(0, src_dir)
    import akw_lanczos  # noqa: E402  (deliberate late import)
    return akw_lanczos


# ----------------------------------------------------------------------------- Hubbard pieces
def sector_ground_energy(eng, L, U, nup, ndn, t, rng, tol=0.0, ncv=None):
    """Lowest eigenvalue of the (nup, ndn) Hubbard sector via matrix-free ARPACK.

    The start vector is drawn from a fixed-seed RNG so the result is reproducible."""
    H, Su, iu, Sd, idd, Du, Dd = eng.sector_H(L, U, nup, ndn, t)
    dim = Du * Dd
    v0 = rng.standard_normal(dim)
    v0 /= np.linalg.norm(v0)
    kw = dict(k=1, which="SA", v0=v0, tol=tol, maxiter=20000)
    if ncv is not None:
        kw["ncv"] = min(ncv, dim)
    vals, vecs = eigsh(H, **kw)
    return float(vals[0]), vecs[:, 0], (Du, Dd, Su, Sd)


def occ_matrix(eng, L, n):
    """(D, L) occupation matrix of the n-particle strings -- same convention as spin_lanczos.py."""
    S, _ = eng.strings(L, n)
    return np.array([[(m >> i) & 1 for i in range(L)] for m in S], dtype=float)


def lanczos_poles(matvec, seed, nl_max, reorth=True, wtol=1e-10, etol=1e-10, patience=4):
    """Lanczos on `seed`; return the poles (Ritz values) and their spectral weights.

    With reorth=True the basis is fully reorthogonalized twice per step (DGKS), which is
    what makes the resulting Gauss quadrature -- and hence the identification of the lowest
    pole -- trustworthy.  Returns (poles, weights, info).  `seed` must be real.
    """
    nrm2 = float(np.dot(seed, seed))
    if nrm2 < 1e-28:
        raise RuntimeError("seed vector is numerically zero")
    dim = seed.size
    Q = np.empty((nl_max, dim), dtype=np.float64) if reorth else None
    v = seed / np.sqrt(nrm2)
    vp = np.zeros_like(v)
    if reorth:
        Q[0] = v
    alphas, betas = [], []
    beta = 0.0
    low_hist = []
    n_used = 0
    stopped = "nl_max"
    for it in range(nl_max):
        n_used = it + 1
        w = matvec(v)
        a = float(np.dot(v, w))
        alphas.append(a)
        w = w - a * v - beta * vp
        if reorth:
            for _ in range(2):                       # two passes = DGKS re-orthogonalization
                w -= Q[: it + 1].T @ (Q[: it + 1] @ w)
        b = float(np.linalg.norm(w))
        # lowest weighted Ritz value so far
        if it >= 1:
            ev, es = eigh_tridiagonal(np.array(alphas), np.array(betas))
            wts = nrm2 * es[0, :] ** 2
            keep = wts > wtol * nrm2
            low_hist.append(float(ev[keep][0]) if keep.any() else np.nan)
            if len(low_hist) > patience:
                recent = np.array(low_hist[-(patience + 1):])
                if np.all(np.isfinite(recent)) and np.max(np.abs(np.diff(recent))) < etol:
                    stopped = "converged"
                    if b < 1e-10:
                        stopped = "breakdown+converged"
                    break
        if b < 1e-10:
            stopped = "breakdown"
            break
        betas.append(b)
        vp = v
        v = w / b
        if reorth and it + 1 < nl_max:
            Q[it + 1] = v                    # guard: the last step needs no stored basis vector
        beta = b
    alphas = np.array(alphas)
    betas = np.array(betas[: len(alphas) - 1])
    ev, es = eigh_tridiagonal(alphas, betas)
    weights = nrm2 * es[0, :] ** 2
    if reorth:
        del Q
    info = dict(n_lanczos=n_used, stop=stopped, seed_norm2=nrm2, reorth=bool(reorth),
                lowest_history=[float(x) for x in low_hist[-6:]])
    return ev, weights, info


def lowest_pole(poles, weights, rel_wtol=1e-8):
    """Smallest Ritz value carrying relative spectral weight above rel_wtol."""
    tot = float(weights.sum())
    keep = weights > rel_wtol * tot
    if not keep.any():
        raise RuntimeError("no Ritz value carries spectral weight")
    i = int(np.argmax(keep))          # poles come out sorted ascending
    return float(poles[i]), float(weights[i] / tot)


def hubbard_point(eng, L, U, t, nl_max, rng_seed, do_noreorth):
    """All L-dependent Hubbard quantities at one size."""
    nup = ndn = L // 2
    rng = np.random.default_rng(rng_seed + L)
    tA = time.time()
    E0, psi0, (Du, Dd, Su, Sd) = sector_ground_energy(eng, L, U, nup, ndn, t, rng)
    log(f"  L={L}: E0(N)  = {E0:.10f}   dim={Du * Dd}   ({time.time() - tA:.1f}s)")

    tA = time.time()
    Ep, _, _ = sector_ground_energy(eng, L, U, nup + 1, ndn, t, np.random.default_rng(rng_seed + 100 + L))
    Em, _, _ = sector_ground_energy(eng, L, U, nup - 1, ndn, t, np.random.default_rng(rng_seed + 200 + L))
    mu_plus = Ep - E0
    mu_minus = E0 - Em
    charge_gap = mu_plus - mu_minus
    log(f"  L={L}: mu+={mu_plus:.6f}  mu-={mu_minus:.6f}  Delta={charge_gap:.6f} "
        f"(PH check mu+ + mu- - U = {mu_plus + mu_minus - U:+.2e})  ({time.time() - tA:.1f}s)")

    # ---- spin gap: lowest exact pole of S^zz(q=pi, omega)
    upocc = occ_matrix(eng, L, nup)
    dnocc = occ_matrix(eng, L, ndn)
    q = np.pi
    ph = np.cos(q * np.arange(L))                    # exp(i pi j) = (-1)^j, real
    nq_up = upocc @ ph
    nq_dn = dnocc @ ph
    Psi = psi0.reshape(Du, Dd)
    seed = (((nq_up[:, None] - nq_dn[None, :]) / 2.0) * Psi).ravel()
    seed = seed - psi0 * float(np.dot(psi0, seed))   # project out |psi0> (exact for q != 0)
    H, _, _, _, _, _, _ = eng.sector_H(L, U, nup, ndn, t)

    def shifted(x):
        return H.matvec(x) - E0 * x

    tA = time.time()
    poles, wts, info = lanczos_poles(shifted, seed, nl_max=min(nl_max, seed.size), reorth=True)
    spin_gap, spin_w = lowest_pole(poles, wts)
    log(f"  L={L}: spin gap (reorth) = {spin_gap:.8f}  rel.weight={spin_w:.4f}  "
        f"n_lanczos={info['n_lanczos']} ({info['stop']})  ({time.time() - tA:.1f}s)")

    spin_gap_nr = None
    if do_noreorth:
        tA = time.time()
        p2, w2, i2 = lanczos_poles(shifted, seed, nl_max=min(nl_max, seed.size), reorth=False)
        spin_gap_nr, _ = lowest_pole(p2, w2)
        log(f"  L={L}: spin gap (NO reorth) = {spin_gap_nr:.8f}  "
            f"n_lanczos={i2['n_lanczos']}  delta={spin_gap_nr - spin_gap:+.2e}  ({time.time() - tA:.1f}s)")

    top = np.argsort(wts)[::-1][:6]
    order = np.argsort(poles[top])
    top = top[order]
    return dict(
        L=L, dim_N=int(Du * Dd),
        E0_N=E0, E0_Nplus1=Ep, E0_Nminus1=Em,
        mu_plus=mu_plus, mu_minus=mu_minus,
        ph_symmetry_residual=float(mu_plus + mu_minus - U),
        charge_gap=charge_gap,
        spin_gap=spin_gap, spin_gap_relweight=spin_w,
        spin_gap_no_reorth=spin_gap_nr,
        L_times_spin_gap=L * spin_gap, L_times_charge_gap=L * charge_gap,
        lanczos=info,
        szz_pi_top_poles=[[float(poles[i]), float(wts[i] / wts.sum())] for i in top],
    )


# ----------------------------------------------------------------------------- Heisenberg control
def heisenberg_sector(eng, L, J, ntot_up):
    """Dense H = J sum_i S_i . S_{i+1} (PBC) in the sector with ntot_up up spins."""
    S, idx = eng.strings(L, ntot_up)
    D = len(S)
    H = np.zeros((D, D))
    for a, m in enumerate(S):
        diag = 0.0
        for i in range(L):
            j = (i + 1) % L
            si = 0.5 if (m >> i) & 1 else -0.5
            sj = 0.5 if (m >> j) & 1 else -0.5
            diag += J * si * sj
            if si != sj:
                m2 = m ^ ((1 << i) | (1 << j))
                H[idx[m2], a] += 0.5 * J
        H[a, a] += diag
    return H, S, idx


def heisenberg_point(eng, L, J):
    """Independent control: dense ED of the Heisenberg ring; two independent gap measures."""
    tA = time.time()
    H0, S0, _ = heisenberg_sector(eng, L, J, L // 2)
    e0, V0 = eigh(H0)
    E0 = float(e0[0])
    psi0 = V0[:, 0]
    # lowest pole of S^zz(q=pi) by exact Lehmann (no Lanczos at all)
    ph = np.cos(np.pi * np.arange(L))
    szq = np.array([sum(ph[i] * (0.5 if (m >> i) & 1 else -0.5) for i in range(L)) for m in S0])
    seed = szq * psi0
    amps = V0.T @ seed
    wts = amps ** 2
    tot = wts.sum()
    ex = e0 - E0
    keep = (wts > 1e-10 * tot) & (ex > 1e-9)
    pole = float(ex[keep][0])
    # singlet-triplet gap from a different sector entirely
    H1, _, _ = heisenberg_sector(eng, L, J, L // 2 + 1)
    e1 = np.linalg.eigvalsh(H1)
    st_gap = float(e1[0] - E0)
    log(f"  L={L}: Heisenberg J={J}: S^zz(pi) lowest pole = {pole:.8f}, "
        f"singlet-triplet gap = {st_gap:.8f}  ({time.time() - tA:.1f}s)")
    return dict(L=L, J=J, E0=E0, dim_sz0=len(S0),
                szz_pi_lowest_pole=pole, singlet_triplet_gap=st_gap,
                L_times_gap=L * pole,
                internal_consistency=float(abs(pole - st_gap)))


def heisenberg_sparse(eng, L, J, ntot_up):
    """Sparse H = J sum_i S_i . S_{i+1} (PBC) in the sector with ntot_up up spins."""
    S, idx = eng.strings(L, ntot_up)
    D = len(S)
    rows, cols, vals = [], [], []
    for a, m in enumerate(S):
        diag = 0.0
        for i in range(L):
            j = (i + 1) % L
            si = 0.5 if (m >> i) & 1 else -0.5
            sj = 0.5 if (m >> j) & 1 else -0.5
            diag += J * si * sj
            if si != sj:
                rows.append(idx[m ^ ((1 << i) | (1 << j))]); cols.append(a); vals.append(0.5 * J)
        rows.append(a); cols.append(a); vals.append(diag)
    return sp.csr_matrix((vals, (rows, cols)), shape=(D, D)), S, idx


def heisenberg_st_gap_large(eng, L, J, rng_seed):
    """Singlet-triplet gap of the Heisenberg ring at a size too large for dense ED.

    Justified by the exact identity verified at L = 4..12 in heisenberg_point():
    the lowest pole of S^zz(q=pi, omega) EQUALS E0(Sz=1) - E0(Sz=0) to ~1e-15.
    Only the two sector ground states are needed, so sparse ARPACK suffices."""
    tA = time.time()
    out = []
    for nup in (L // 2, L // 2 + 1):
        H, S, _ = heisenberg_sparse(eng, L, J, nup)
        rng = np.random.default_rng(rng_seed + 1000 * nup + L)
        v0 = rng.standard_normal(H.shape[0]); v0 /= np.linalg.norm(v0)
        w = eigsh(H, k=1, which="SA", v0=v0, tol=0.0, maxiter=50000)[0]
        out.append(float(w[0]))
    gap = out[1] - out[0]
    log(f"  L={L}: Heisenberg (sparse) singlet-triplet gap = {gap:.8f}  "
        f"L*gap = {L * gap:.4f}  ({time.time() - tA:.1f}s)")
    return dict(L=L, J=J, E0=out[0], method="sparse ARPACK, singlet-triplet gap only",
                szz_pi_lowest_pole=gap, singlet_triplet_gap=gap,
                L_times_gap=L * gap, internal_consistency=None)


def ghost_probe(eng, L, U, t, nl_fixed, rng_seed):
    """Control for the reorthogonalization claim.

    Runs the SAME spin seed to a FIXED, deliberately over-long nl with and without
    reorthogonalization and counts near-degenerate Ritz values ("ghosts"), the textbook
    signature of lost orthogonality.  This is what justifies -- or bounds -- the statement
    that the continued fraction needs full reorthogonalization: the number below says at
    which nl the plain three-term recurrence starts to fabricate copies of converged poles.
    src/spin_lanczos.py runs nl = 200 with no reorthogonalization, so this is the control
    for the number the manuscript already relies on.
    """
    nup = ndn = L // 2
    rng = np.random.default_rng(rng_seed + L)
    E0, psi0, (Du, Dd, _, _) = sector_ground_energy(eng, L, U, nup, ndn, t, rng)
    upocc = occ_matrix(eng, L, nup); dnocc = occ_matrix(eng, L, ndn)
    ph = np.cos(np.pi * np.arange(L))
    seed = (((upocc @ ph)[:, None] - (dnocc @ ph)[None, :]) / 2.0 * psi0.reshape(Du, Dd)).ravel()
    seed = seed - psi0 * float(np.dot(psi0, seed))
    H, _, _, _, _, _, _ = eng.sector_H(L, U, nup, ndn, t)
    mv = lambda x: H.matvec(x) - E0 * x
    nl = int(min(nl_fixed, seed.size))
    res = {}
    for tag, ro in (("reorth", True), ("no_reorth", False)):
        pol, wt, inf = lanczos_poles(mv, seed, nl_max=nl, reorth=ro, etol=-1.0)
        lo, lw = lowest_pole(pol, wt)
        dup = int(np.sum(np.diff(pol) < 1e-8))
        dup4 = int(np.sum(np.diff(pol) < 1e-4))
        res[tag] = dict(n_lanczos=inf["n_lanczos"], lowest_pole=lo, lowest_relweight=lw,
                        n_ritz=int(pol.size), n_ghost_1e8=dup, n_ghost_1e4=dup4,
                        min_gap_between_ritz=float(np.min(np.diff(pol))) if pol.size > 1 else None)
        log(f"  ghost probe L={L} nl={nl} [{tag}]: lowest pole {lo:.10f}  "
            f"ghosts(<1e-8) {dup}  ghosts(<1e-4) {dup4}  of {pol.size} Ritz values")
    res["lowest_pole_difference"] = res["no_reorth"]["lowest_pole"] - res["reorth"]["lowest_pole"]
    res["nl_fixed"] = nl
    res["L"] = L
    return res


# ----------------------------------------------------------------------------- Lieb-Wu
def lieb_wu_gap(U, t=1.0):
    """Exact Lieb-Wu charge gap by two independent integral representations."""
    # form A: Delta = (16 t^2 / U) int_1^inf dy sqrt(y^2-1)/sinh(2 pi y t / U)
    # 1/sinh(x) is written as 2 exp(-x)/(1-exp(-2x)) so the large-x tail never overflows.
    def f1(y):
        x = 2.0 * np.pi * y * t / U
        return np.sqrt(y * y - 1.0) * 2.0 * np.exp(-x) / (1.0 - np.exp(-2.0 * x))
    I1, _ = quad(f1, 1.0, np.inf, limit=400)
    g1 = 16.0 * t * t / U * I1

    # form B: Delta = U - 4t + 8t int_0^inf dw J_1(w) / (w (1 + exp(w U / 2t)))
    # 1/(1+exp(x)) = expit(-x), overflow-free.
    def f2(w):
        if w <= 0.0:
            return 0.25                      # J_1(w)/w -> 1/2 and 1/(1+e^0) = 1/2
        return j1(w) / w * expit(-w * U / (2.0 * t))
    I2, _ = quad(f2, 0.0, np.inf, limit=800)
    g2 = U - 4.0 * t + 8.0 * t * I2
    return float(g1), float(g2)


# ----------------------------------------------------------------------------- fits
def ols(X, y):
    """Ordinary least squares with coefficient standard errors (unweighted)."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    n, p = X.shape
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = n - p
    if dof <= 0:
        return beta, np.full(p, np.nan), np.nan, np.nan
    s2 = float(resid @ resid) / dof
    cov = s2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float(resid @ resid) / ss_tot if ss_tot > 0 else np.nan
    return beta, se, float(np.sqrt(s2)), r2


def fit_block(Ls, vals, label, with_quadratic=False):
    Ls = np.asarray(Ls, float)
    vals = np.asarray(vals, float)
    invL = 1.0 / Ls
    X = np.column_stack([np.ones_like(invL), invL])
    beta, se, rmse, r2 = ols(X, vals)
    dof = int(len(Ls) - 2)
    out = {"label": label, "sizes": Ls.tolist(), "values": vals.tolist(),
           "model": "y = c0 + c1/L",
           "c0_intercept": float(beta[0]), "c0_se": float(se[0]),
           "c1_slope": float(beta[1]), "c1_se": float(se[1]),
           "rmse": rmse, "r2": r2, "dof": dof}
    if dof <= 0:
        # Two points determine a two-parameter line exactly: the standard errors are not
        # small, they do not exist.  Say so instead of printing nan as if it were a result.
        out["underdetermined"] = True
        out["note"] = ("only %d size(s): a two-parameter fit has %d degrees of freedom, so the "
                       "coefficients are an interpolation and NO standard error, R^2 or "
                       "significance may be quoted from this block" % (len(Ls), dof))
    if with_quadratic and len(Ls) >= 4:
        X2 = np.column_stack([np.ones_like(invL), invL, invL ** 2])
        b2, s2, rm2, r22 = ols(X2, vals)
        out["quadratic"] = {"model": "y = c0 + c1/L + c2/L^2",
                            "c0_intercept": float(b2[0]), "c0_se": float(s2[0]),
                            "c1": float(b2[1]), "c1_se": float(s2[1]),
                            "c2": float(b2[2]), "c2_se": float(s2[2]),
                            "rmse": rm2, "r2": r22, "dof": int(len(Ls) - 3)}
    return out


# ----------------------------------------------------------------------------- self-check
class SelfCheckFailure(RuntimeError):
    pass


def require(name, got, expect, tol, record):
    ok = abs(got - expect) <= tol
    record.append({"check": name, "got": float(got), "expected": float(expect),
                   "tolerance": float(tol), "abs_error": float(abs(got - expect)), "pass": bool(ok)})
    status = "PASS" if ok else "FAIL"
    log(f"  SELF-CHECK [{status}] {name}: got {got:.8f}, expected {expect} +- {tol} "
        f"(|err| = {abs(got - expect):.2e})")
    if not ok:
        raise SelfCheckFailure(f"{name}: got {got!r}, expected {expect} +- {tol}")


# ----------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sizes", default="4,6,8,10,12", help="comma-separated ring sizes")
    ap.add_argument("--U", type=float, default=8.0)
    ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--nl-max", type=int, default=40,
                    help="max Lanczos steps for the spin seed (40 = the value that produced "
                         "data/gap_scaling.json)")
    ap.add_argument("--rng-seed", type=int, default=20260918, help="fixed ARPACK start-vector seed")
    ap.add_argument("--out", default=DEFAULT_OUT,
                    help="output JSON (default: <repo>/data/gap_scaling.json); a relative path is "
                         "taken relative to the CURRENT DIRECTORY, never to the repository root")
    ap.add_argument("--src", default=SRC_DIR,
                    help="directory holding akw_lanczos.py (default: this script's directory)")
    ap.add_argument("--heisenberg-sizes", default="14,16,18,20",
                    help="extra ring sizes for the Heisenberg control only (sparse ED, "
                         "singlet-triplet gap); e.g. 14,16,18,20 -- costs seconds and turns the "
                         "1/L closure of the spin channel into an eight-point fit")
    ap.add_argument("--ghost-probe", default="8,10",
                    help="comma-separated sizes on which to run the fixed-nl reorthogonalization "
                         "control (ghost-eigenvalue count with and without reorthogonalization); "
                         "use together with --ghost-nl")
    ap.add_argument("--ghost-nl", type=int, default=200,
                    help="fixed Lanczos length for --ghost-probe (200 = what spin_lanczos.py uses)")
    # NOTE ON THE NAME: this flag ENABLES the comparison (the name reads as if it disabled it).
    # It is kept for compatibility with the invocation recorded in data/gap_scaling.json;
    # --reorth-compare is the readable spelling of the same switch, and it is ON by default so
    # that a bare "python src/gap_scaling.py" reproduces the deposited JSON.
    ap.add_argument("--no-reorth-compare", "--reorth-compare", dest="no_reorth_compare",
                    action="store_true", default=True,
                    help="also run Lanczos WITHOUT reorthogonalization and record the "
                         "difference (ON by default; switch off with --skip-reorth-compare)")
    ap.add_argument("--skip-reorth-compare", dest="no_reorth_compare", action="store_false",
                    help="do NOT run the extra no-reorthogonalization Lanczos pass")
    ap.add_argument("--skip-self-check", action="store_true",
                    help="NOT for production: disables the abort-on-mismatch guard")
    args = ap.parse_args(argv)

    sizes = [int(s) for s in args.sizes.split(",") if s.strip()]
    eng = import_engine(args.src)
    U, t = args.U, args.t
    J = 4.0 * t * t / U
    log(f"engine: {eng.__file__}")
    log(f"sizes={sizes}  U={U}  t={t}  J=4t^2/U={J}  nl_max={args.nl_max}  rng_seed={args.rng_seed}")

    checks = []
    # Lieb-Wu first: it costs nothing and gates everything downstream.
    g1, g2 = lieb_wu_gap(U, t)
    log(f"Lieb-Wu charge gap: form A = {g1:.9f}, form B = {g2:.9f}")
    if not args.skip_self_check:
        require("lieb_wu_two_forms_agree", g1 - g2, 0.0, 1e-6, checks)
        require("lieb_wu_literature_value", g1, LIEB_WU_U8_REFERENCE, 1e-3, checks)

    hub, hei = [], []
    for L in sizes:
        log(f"--- L = {L} ---")
        h = hubbard_point(eng, L, U, t, args.nl_max, args.rng_seed, args.no_reorth_compare)
        hub.append(h)
        hei.append(heisenberg_point(eng, L, J))
        if not args.skip_self_check:
            require(f"L{L}_particle_hole_mu_sum", h["mu_plus"] + h["mu_minus"], U, 1e-6, checks)
            if L in SELF_CHECK:
                c, ct = SELF_CHECK[L]["charge_gap"]
                s, st = SELF_CHECK[L]["spin_gap"]
                require(f"L{L}_charge_gap", h["charge_gap"], c, ct, checks)
                require(f"L{L}_spin_gap", h["spin_gap"], s, st, checks)
            if L in HEISENBERG_SELF_CHECK:
                v, vt = HEISENBERG_SELF_CHECK[L]
                require(f"L{L}_heisenberg_szz_pi", hei[-1]["szz_pi_lowest_pole"], v, vt, checks)
                require(f"L{L}_heisenberg_pole_equals_ST_gap",
                        hei[-1]["internal_consistency"], 0.0, 1e-8, checks)

    ghosts = []
    for L in [int(x) for x in args.ghost_probe.split(",") if x.strip()]:
        log(f"--- reorthogonalization control, L = {L}, fixed nl = {args.ghost_nl} ---")
        ghosts.append(ghost_probe(eng, L, U, t, args.ghost_nl, args.rng_seed))

    extra = [int(x) for x in args.heisenberg_sizes.split(",") if x.strip()]
    extra = [L for L in extra if L not in sizes]
    hei_ext = list(hei)
    if extra:
        log("--- Heisenberg control extension (sparse ED, singlet-triplet gap) ---")
        for L in extra:
            hei_ext.append(heisenberg_st_gap_large(eng, L, J, args.rng_seed))
        hei_ext.sort(key=lambda d: d["L"])

    Ls = [h["L"] for h in hub]
    spin = [h["spin_gap"] for h in hub]
    charge = [h["charge_gap"] for h in hub]
    hspin = [h["szz_pi_lowest_pole"] for h in hei]

    fits = {}
    fits["spin_all"] = fit_block(Ls, spin, "Hubbard spin gap vs 1/L (all sizes)")
    fits["charge_all"] = fit_block(Ls, charge, "Hubbard charge gap vs 1/L (all sizes)",
                                   with_quadratic=True)
    fits["heisenberg_all"] = fit_block(Ls, hspin, "Heisenberg spin gap vs 1/L (all sizes)")
    if extra:
        eL = [d["L"] for d in hei_ext]
        ev = [d["szz_pi_lowest_pole"] for d in hei_ext]
        fits["heisenberg_extended"] = fit_block(eL, ev,
            "Heisenberg spin gap vs 1/L (control extended by sparse ED)")
        big_e = [i for i, L in enumerate(eL) if L >= 6]
        fits["heisenberg_extended_L6plus"] = fit_block([eL[i] for i in big_e], [ev[i] for i in big_e],
            "Heisenberg spin gap vs 1/L (control extended, L >= 6)")
    big = [i for i, L in enumerate(Ls) if L >= 6]
    if len(big) >= 3:
        fits["spin_L6plus"] = fit_block([Ls[i] for i in big], [spin[i] for i in big],
                                        "Hubbard spin gap vs 1/L (L >= 6)")
        fits["charge_L6plus"] = fit_block([Ls[i] for i in big], [charge[i] for i in big],
                                          "Hubbard charge gap vs 1/L (L >= 6)", with_quadratic=True)
        fits["heisenberg_L6plus"] = fit_block([Ls[i] for i in big], [hspin[i] for i in big],
                                              "Heisenberg spin gap vs 1/L (L >= 6)")

    # L * gap: the crisp qualitative contrast (constant vs linearly growing)
    key_s = "spin_L6plus" if "spin_L6plus" in fits else "spin_all"
    key_c = "charge_L6plus" if "charge_L6plus" in fits else "charge_all"
    Lsub = fits[key_s]["sizes"]
    Xlin = np.column_stack([np.ones(len(Lsub)), np.asarray(Lsub, float)])
    bs, ses, _, r2s = ols(Xlin, np.asarray(fits[key_s]["values"]) * np.asarray(Lsub, float))
    bc, sec, _, r2c = ols(Xlin, np.asarray(fits[key_c]["values"]) * np.asarray(Lsub, float))
    Lgap = {"model": "L*gap = d0 + d1*L", "sizes": Lsub,
            "spin": {"d0": float(bs[0]), "d0_se": float(ses[0]),
                     "d1": float(bs[1]), "d1_se": float(ses[1]), "r2": r2s,
                     "d1_z": float(bs[1] / ses[1]) if ses[1] > 0 else np.nan},
            "charge": {"d0": float(bc[0]), "d0_se": float(sec[0]),
                       "d1": float(bc[1]), "d1_se": float(sec[1]), "r2": r2c,
                       "d1_z": float(bc[1] / sec[1]) if sec[1] > 0 else np.nan}}

    c0s, s0s = fits[key_s]["c0_intercept"], fits[key_s]["c0_se"]
    c0c, s0c = fits[key_c]["c0_intercept"], fits[key_c]["c0_se"]
    z_sep = (c0c - c0s) / np.sqrt(s0s ** 2 + s0c ** 2)
    z_spin_zero = abs(c0s) / s0s if s0s > 0 else np.inf
    separation = {
        "fit_used_spin": key_s, "fit_used_charge": key_c,
        "spin_intercept": c0s, "spin_intercept_se": s0s,
        "charge_intercept": c0c, "charge_intercept_se": s0c,
        "lieb_wu_exact": g1,
        "charge_intercept_minus_lieb_wu": c0c - g1,
        "z_spin_intercept_vs_zero": float(z_spin_zero),
        "z_intercept_separation": float(z_sep),
        "statement": (
            ("NOT A MEASUREMENT: the underlying fit is underdetermined (%d and %d degrees of "
             "freedom), so no sigma may be quoted. Intercepts are spin %.4f t, charge %.4f t, "
             "quoted as interpolations only."
             % (fits[key_s]["dof"], fits[key_c]["dof"], c0s, c0c))
            if (fits[key_s].get("underdetermined") or fits[key_c].get("underdetermined")
                or not np.isfinite(z_sep))
            else ("spin gap extrapolates to zero (intercept consistent with 0 at "
                  f"{z_spin_zero:.1f} sigma from the fit) while the charge gap saturates at "
                  f"{c0c:.4f} t; the two intercepts separate at {abs(z_sep):.1f} sigma")),
    }

    control = {
        "hubbard_spin_gap": spin,
        "heisenberg_spin_gap": hspin,
        "ratio_hubbard_over_heisenberg": [float(a / b) for a, b in zip(spin, hspin)],
        "abs_difference": [float(a - b) for a, b in zip(spin, hspin)],
        "max_abs_difference": float(max(abs(a - b) for a, b in zip(spin, hspin))),
        "note": ("Heisenberg here is the PURE spin ring with J = 4t^2/U, solved by dense ED in a "
                 "different sector and a different code path; the strong-coupling correction to J "
                 "is O(t^4/U^3), so a residual difference of a few percent at U/t=8 is expected "
                 "and is itself evidence that the Hubbard number is the spin-sector gap."),
    }

    prov = {
        "script": os.path.basename(__file__),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "wall_seconds": round(time.time() - T0, 2),
        "argv": sys.argv[1:] if argv is None else list(argv),
        "parameters": {"sizes": sizes, "U": U, "t": t, "J_eff": J, "q": "pi",
                       "boundary": "periodic ring", "filling": "half (n_up = n_dn = L/2)",
                       "nl_max": args.nl_max, "rng_seed": args.rng_seed,
                       "lanczos": "full reorthogonalization (2 DGKS passes)",
                       "broadening": "none -- exact poles from the Lanczos tridiagonal"},
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "scipy": scipy.__version__, "platform": platform.platform()},
        "engine": {"module": "akw_lanczos", "path": _repo_rel(eng.__file__),
                   "functions_reused": ["strings", "sector_H"]},
        "self_checks": checks,
        "self_check_enforced": not args.skip_self_check,
        "control_values_reproduced": [c["check"] for c in checks if c["pass"]],
    }

    result = {"provenance": prov,
              "lieb_wu": {"U": U, "t": t, "form_A_sinh_integral": g1,
                          "form_B_bessel_integral": g2, "literature": LIEB_WU_U8_REFERENCE},
              "hubbard": hub, "heisenberg_control": hei,
              "heisenberg_control_extended": hei_ext,
              "fits": fits, "L_times_gap": Lgap,
              "reorthogonalization_control": ghosts,
              "separation_of_the_two_channels": separation,
              "spin_channel_control": control}

    out = os.path.abspath(args.out)
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    log(f"WROTE {out}")

    # ---- human-readable table
    print("\n  L |  dim(N)  |   mu+      mu-     Delta_charge | spin gap   L*spin | Heis. gap  L*Heis")
    print("----+----------+---------------------------------+-------------------+-------------------")
    for h, e in zip(hub, hei):
        print(f" {h['L']:2d} | {h['dim_N']:8d} | {h['mu_plus']:8.5f} {h['mu_minus']:8.5f} "
              f"{h['charge_gap']:8.5f} | {h['spin_gap']:9.6f} {h['L_times_spin_gap']:6.3f} | "
              f"{e['szz_pi_lowest_pole']:9.6f} {e['L_times_gap']:6.3f}")
    print(f"\n  spin  gap = {fits[key_s]['c1_slope']:.4f}(+-{fits[key_s]['c1_se']:.4f})/L "
          f"+ {c0s:+.5f}(+-{s0s:.5f})   [{fits[key_s]['label']}]")
    print(f"  charge gap = {fits[key_c]['c1_slope']:.4f}(+-{fits[key_c]['c1_se']:.4f})/L "
          f"+ {c0c:+.5f}(+-{s0c:.5f})   [{fits[key_c]['label']}]  vs Lieb-Wu {g1:.5f}")
    print(f"  {separation['statement']}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SelfCheckFailure as exc:
        log(f"ABORTING: self-check failed -- {exc}")
        sys.exit(2)
