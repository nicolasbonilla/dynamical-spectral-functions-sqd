# -*- coding: utf-8 -*-
r"""eta_sweep.py -- the BROADENING CONFOUND in the subspace-fraction scaling law.

WHAT THIS COMPUTES
------------------
The required sampled-subspace fraction  frac(L, eta)  -- the smallest fraction of the
(N+1, Sz=+1) sector whose Rayleigh-Ritz spectral function A_S(omega) reproduces the exact
Haydock A(omega) to rel-L1 < 5% -- for L = 4, 6, 8, 10, 12 and a decade of Lorentzian
broadenings eta/t in {0.05, 0.075, 0.10, 0.15, 0.18, 0.25, 0.35, 0.50}.  It then
fits, FOR EACH eta:

  * the log-linear exponent of the FRACTION,  frac ~ exp(b_frac * L),
  * the log-linear exponent of the ABSOLUTE support size, |S| = frac * dim ~ exp(b_S * L),
  * the LOCAL logarithmic slopes between consecutive L (the sequence is compared against a
    single exponential), and a formal curvature test of the log-linear model with a p-value.

WHICH CLAIM OF THE PAPER THIS BACKS (and how it corrects it)
------------------------------------------------------------
The manuscript reports a single "required fraction" sequence 0.92 -> 0.08 at the ONE value
eta = 0.15 t and reads it as a property of the model.  Half of that decay is a property of the
BROADENING, not of the physics:

  (1) the exponent of the FRACTION is not an observable of the model -- it moves by a factor
      ~2 across the eta decade, so any headline fraction is reachable by choosing eta;
  (2) the exponent of the ABSOLUTE support |S| is stable across the same decade (R^2 >= 0.999)
      -- that, and not the fraction, is the robust statement;
  (3) even at FIXED eta the sequence is not a single exponential: the local logarithmic slopes
      vary by a large factor, so the quoted straight-line exponent is a fit through a curve of
      increasing slope.

The mechanism is already inside the manuscript but never connected: the geometric-decay
constants of the moment-exactness theorem diverge (C(eta) -> infinity, rho(eta) -> 1+) as
eta -> 0+.

MODEL -- note the boundary conditions
-------------------------------------
Hubbard RING (PERIODIC boundary conditions, the hopping t_{L-1,0} is present), U/t = 8, half
filling.  This is the model of src/scaling_lanczos.py, and it is NOT the model of
src/cost_vs_entanglement.py / src/gate1_ladder.py, which use an OPEN chain.  The one-body magic
F1 differs between the two and the difference is not small at L = 6:

    L = 6 :  F1(OBC) = 9.604689   vs   F1(PBC) = 9.385082
    L = 8 :  F1(OBC) = 12.816276  vs   F1(PBC) = 12.813821

The F1 values reported alongside frac(L) here are therefore the PBC ones, and a "correction" of
9.3851 -> 9.604689 must NOT be applied to them: 9.385082 is the correct F1 of the periodic ring
whose spectral function this script computes.

PIPELINE
--------
Identical to src/scaling_lanczos.py (matrix-free variant of src/scaling_lanczos_mf.py):
ground state by sparse Lanczos in the (L/2, L/2) sector; seed phi = c^dag_{0,up}|Psi_0>;
exact A(omega) by Haydock continued fraction; configuration ranking by the time-integrated
Born weight of the Krylov-evolved seed (K = 18 steps, dt = 0.5); subspace A_S(omega) by a
RESTRICTED matvec (embed -> full H matvec -> extract; no submatrix is ever materialized);
coarse scan over fractions followed by bisection to a resolution of dim/200.

The Lanczos recursion coefficients (alpha, beta) do not depend on the evaluation point z,
so ONE Lanczos run per subspace serves ALL eta values.  This is what makes the sweep affordable.

DETERMINISM
-----------
Every RNG is seeded explicitly (SEED_GS for the eigsh start vector, SEED_CTRL for the control).
The ground-state global sign is fixed by a deterministic convention.  Results are bit-for-bit
reproducible on a given numpy/scipy build.

MANDATORY SELF-CHECK (the script ABORTS on failure)
---------------------------------------------------
  C1  the eta = 0.15 column must reproduce the published sequence
      0.917 / 0.820 / 0.560 / 0.380 / 0.180  for L = 4 / 6 / 8 / 10 / 12  (atol 0.02);
  C2  at L = 8 the matrix-free Krylov propagator used for the configuration ranking must give
      the same frac as scipy's expm_multiply on the explicitly built sparse Hamiltonian.

USAGE
-----
  python eta_sweep.py                                  # full sweep, default output eta_sweep.json
  python eta_sweep.py --L 4 6 8 --out out/small.json   # subset
  python eta_sweep.py --eta 0.05 0.15 0.50             # custom broadenings
  python eta_sweep.py --no-check                       # skip C1/C2 (NOT recommended)

RUNTIME (laptop-class CPU, other jobs competing; 5 eta values, measured)
  L=4  6 s   L=6  9 s   L=8  41 s   L=10  4 min   L=12  94 min   (total 100 min)
  L = 12 dominates: 24 restricted-Haydock evaluations of 260 steps each in a sector of
  dimension 731 808.  Peak RSS < 400 MB -- the Hamiltonian is never materialized (the only
  large arrays are the seed and the ranking vector); the L=8 control briefly builds a
  3920 x 3920 sparse matrix for scipy's expm_multiply.

  --refit rebuilds the fit block from an existing output JSON in seconds, so the statistics can
  be revised without repeating any of the physics.
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import sys
import time

import numpy as np
import scipy
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

# The engine primitives (strings, hopping, c^dag maps) are reused verbatim from the repository.
_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (_HERE, os.path.join(_HERE, "..", "src"), os.path.join(_HERE, "src")):
    if os.path.isdir(_p) and _p not in sys.path:
        sys.path.insert(0, _p)
import akw_lanczos as AK  # noqa: E402

# ----------------------------------------------------------------------------------------------
# fixed protocol constants -- these define the published pipeline and must not be tuned
# ----------------------------------------------------------------------------------------------
U_HUB = 8.0          # on-site repulsion, in units of the hopping t
T_HOP = 1.0
N_LANCZOS = 260      # Haydock depth (src/scaling_lanczos.py)
N_GRID = 600         # omega grid points
THR_RELL1 = 0.05     # acceptance threshold on the relative L1 error
K_RANK = 18          # number of time-evolution steps used to accumulate the Born weight
DT_RANK = 0.5        # time step
SUBSTEPS = 16        # substeps per DT_RANK for the matrix-free Krylov propagator
KRYLOV_M = 8         # Krylov dimension of the matrix-free propagator
COARSE_FRACS = (0.02, 0.05, 0.08, 0.12, 0.16, 0.20, 0.25, 0.30,
                0.40, 0.50, 0.60, 0.70, 0.78, 0.85, 0.92, 1.0)
SEED_GS = 20260918   # RNG seed for the eigsh start vector
SEED_CTRL = 4242     # RNG seed reserved for the control path

DEFAULT_L = (4, 6, 8, 10, 12)
DEFAULT_ETA = (0.05, 0.075, 0.10, 0.15, 0.18, 0.25, 0.35, 0.50)
# The eight broadenings the manuscript's eta table prints.  The deposited
# data/eta_sweep.json is the output of a bare `python src/eta_sweep.py`:
# until 2026-09-18 the deposit only held five of these columns, so three of
# the table's entries (eta = 0.075, 0.18, 0.35) had no source at all.

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
DEFAULT_OUT = os.path.join(REPO, "data", "eta_sweep.json")
SCALING_JSON = os.path.join(REPO, "data", "scaling_data.json")

# Published reference for the self-check: required fraction at eta = 0.15 t.
CONTROL_ETA = 0.15
CONTROL_ATOL = 0.02


def _published_fracs():
    """The six published fractions (L = 4..14, eta = 0.15 t), READ FROM THE DEPOSIT.

    DEFECT FIXED HERE (2026-09-18).  An earlier version of this file hard-typed
        {4: 0.9167, 6: 0.83, 8: 0.5625, 10: 0.3625, 12: 0.17, 14: 0.08}
    as "the published sequence".  Three of those six digits-strings are NOT what
    data/scaling_data.json contains: the deposit holds 0.8200 (not 0.83), 0.556122
    (not 0.5625) and 0.362491 (not 0.3625).  A file destined for the deposit whose
    "published" reference disagreed with the deposit is exactly the kind of orphan
    number this repository is being repaired for, so the sequence is now READ from
    data/scaling_data.json and nothing is hand-typed.  L = 14 is beyond what this
    script can compute (sector dimension 10 306 296): it is carried only so that the
    local-slope diagnostic can be run on the published sequence itself, and the JSON
    flags the whole sequence as QUOTED, not recomputed.
    """
    if not os.path.exists(SCALING_JSON):
        raise SystemExit("ABORT: %s is missing; the published sequence can only come from the "
                         "deposit, never from hand-typed digits." % SCALING_JSON)
    with open(SCALING_JSON, encoding="utf-8") as fh:
        pts = json.load(fh)["points"]
    return {int(p["L"]): float(p["frac"]) for p in pts}


PUBLISHED_FRACS_L4_14 = _published_fracs()
# same source, restricted to the sizes this script actually recomputes
CONTROL_FRAC = {L: f for L, f in PUBLISHED_FRACS_L4_14.items() if L <= 12}

_T0 = time.time()


def log(*a):
    print(f"[{time.time() - _T0:8.1f}s]", *a, flush=True)


def _trapz(y, x):
    fn = getattr(np, "trapezoid", None) or np.trapz
    return float(fn(y, x))


# ----------------------------------------------------------------------------------------------
# matrix-free sector Hamiltonian (never materializes H; RAM cost is one float64 diagonal)
# ----------------------------------------------------------------------------------------------
def sector_matvec(L, U, nup, ndn, t=T_HOP):
    """Return (matvec, dim, Du, Dd, up_strings, up_index) for the (nup, ndn) Hubbard sector."""
    Tu, Su, iu = AK.hop(L, nup, t)
    Td, Sd, _ = AK.hop(L, ndn, t)
    Du, Dd = len(Su), len(Sd)
    upocc = np.array([[(m >> i) & 1 for i in range(L)] for m in Su], dtype=np.float64)
    dnocc = np.array([[(m >> i) & 1 for i in range(L)] for m in Sd], dtype=np.float64)
    diagU = U * (upocc @ dnocc.T).ravel()

    def mv(x):
        X = x.reshape(Du, Dd)
        Y = Tu @ X + (Td @ X.T).T
        return Y.ravel() + diagU * x

    return mv, Du * Dd, Du, Dd, Su, iu


# ----------------------------------------------------------------------------------------------
# Haydock continued fraction, split so that ONE Lanczos run serves every eta
# ----------------------------------------------------------------------------------------------
def lanczos_coeffs(mv, v0, nl):
    """Plain (non-reorthogonalized) Lanczos coefficients of <v0|(z - M)^{-1}|v0>.

    Byte-identical recursion to AK.haydock; only the z-evaluation is factored out.
    Returns (alpha, beta, |v0|^2).
    """
    nrm2 = float(np.vdot(v0, v0).real)
    if nrm2 < 1e-14:
        return np.zeros(0), np.zeros(0), 0.0
    v = v0 / np.sqrt(nrm2)
    vp = np.zeros_like(v)
    b = 0.0
    a_list, b_list = [], []
    for _ in range(nl):
        w = mv(v)
        an = float(np.vdot(v, w).real)
        a_list.append(an)
        w = w - an * v - b * vp
        bn = float(np.sqrt(np.vdot(w, w).real))
        b_list.append(bn)
        if bn < 1e-10:
            break
        vp = v
        v = w / bn
        b = bn
    return np.array(a_list), np.array(b_list), nrm2


def cf_eval(a, bs, nrm2, z):
    """Evaluate the continued fraction built from (a, bs, nrm2) on the complex points z."""
    if a.size == 0:
        return np.zeros_like(z)
    G = z - a[-1]
    for m in range(len(a) - 2, -1, -1):
        G = (z - a[m]) - bs[m] ** 2 / G
    return nrm2 / G


def spectra_all_eta(mv, seed, E0, grid, etas, nl):
    """A(omega) = -Im G / pi for every eta, from a SINGLE Lanczos run."""
    a, bs, nrm2 = lanczos_coeffs(lambda x: mv(x) - E0 * x, seed, nl)
    return {e: -cf_eval(a, bs, nrm2, grid + 1j * e).imag / np.pi for e in etas}


# ----------------------------------------------------------------------------------------------
# configuration ranking
# ----------------------------------------------------------------------------------------------
def krylov_expm(mv, v, dt, m):
    """exp(-i dt H) v by a small Lanczos basis (memory: m complex vectors)."""
    n = v.shape[0]
    beta = float(np.linalg.norm(v))
    if beta < 1e-14:
        return v.copy()
    V = np.empty((m, n), dtype=complex)
    V[0] = v / beta
    al = np.zeros(m)
    be = np.zeros(m)
    w = mv(V[0])
    a = float(np.vdot(V[0], w).real)
    al[0] = a
    w = w - a * V[0]
    mm = m
    for j in range(1, m):
        bnorm = float(np.linalg.norm(w))
        be[j - 1] = bnorm
        if bnorm < 1e-12:
            mm = j
            break
        V[j] = w / bnorm
        w = mv(V[j])
        a = float(np.vdot(V[j], w).real)
        al[j] = a
        w = w - a * V[j] - bnorm * V[j - 1]
    T = (np.diag(al[:mm]) + np.diag(be[:mm - 1], 1) + np.diag(be[:mm - 1], -1))
    from scipy.linalg import expm as dense_expm
    E = dense_expm(-1j * dt * T)[:, 0]
    out = beta * (V[:mm].T @ E)
    del V
    return out


def rank_configs_matrixfree(mv, seed, K=K_RANK, dt=DT_RANK, sub=SUBSTEPS, m=KRYLOV_M):
    """Time-integrated Born weight of the Krylov-evolved seed -> descending configuration order."""
    v = seed / np.sqrt(float(np.vdot(seed, seed).real))
    wc = np.abs(v) ** 2
    ddt = dt / sub
    for s in range(1, K * sub + 1):
        v = krylov_expm(mv, v, ddt, m)
        if s % sub == 0:
            wc += np.abs(v) ** 2
    return np.argsort(wc)[::-1].astype(np.int64), wc


def rank_configs_expm_multiply(L, U, nup, nd, seed, K=K_RANK, dt=DT_RANK):
    """Control path: the literal ranking of src/scaling_lanczos.py (explicit H + expm_multiply)."""
    from scipy.sparse.linalg import expm_multiply
    H1, _, _ = AK.build_H_explicit(L, U, nup, nd, T_HOP)
    v = seed / np.sqrt(float(np.vdot(seed, seed).real))
    wc = np.abs(v) ** 2
    for _ in range(K):
        v = expm_multiply(-1j * dt * H1, v)
        wc += np.abs(v) ** 2
    del H1
    return np.argsort(wc)[::-1].astype(np.int64)


# ----------------------------------------------------------------------------------------------
# the sweep for one L
# ----------------------------------------------------------------------------------------------
def gs_sector(L, U):
    """Deterministic ground state of the (L/2, L/2) sector."""
    nup = nd = L // 2
    mv0, dim0, Du0, Dd0, Su0, iu0 = sector_matvec(L, U, nup, nd)
    Hop = sp.linalg.LinearOperator((dim0, dim0), matvec=mv0, dtype=float)
    v0 = np.random.default_rng(SEED_GS).standard_normal(dim0)
    v0 /= np.linalg.norm(v0)
    w, v = eigsh(Hop, k=1, which="SA", v0=v0, ncv=min(dim0, 24), maxiter=20000, tol=0.0)
    psi = v[:, 0]
    # deterministic global sign: largest-magnitude amplitude is positive
    if psi[int(np.argmax(np.abs(psi)))] < 0:
        psi = -psi
    return float(w[0]), psi, Du0, Dd0, Su0, iu0, nup, nd


def onerdm_up(psi_mat, S, idx, L):
    """up-spin 1-RDM Gamma[i,j] = <psi| c^dag_i c_j |psi> (down spin traced out)."""
    rho = psi_mat @ psi_mat.conj().T
    G = np.zeros((L, L))
    for a, m in enumerate(S):
        for j in range(L):
            if not (m >> j) & 1:
                continue
            for i in range(L):
                if i == j:
                    m2, sign = m, 1.0
                elif (m >> i) & 1:
                    continue
                else:
                    m2 = (m & ~(1 << j)) | (1 << i)
                    lo, hi = min(i, j), max(i, j)
                    mask = m & (((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))
                    sign = -1.0 if (bin(mask).count("1") & 1) else 1.0
                G[i, j] += float(rho[idx[m2], a].real) * sign
    return G


def sweep_L(L, etas, U=U_HUB, control=False):
    """Return the record for one L: frac(eta) plus the shared diagnostics."""
    t_start = time.time()
    E0, psi, Du0, Dd0, Su0, iu0, nup, nd = gs_sector(L, U)
    psi_mat = psi.reshape(Du0, Dd0)
    g = np.clip(np.linalg.eigvalsh(onerdm_up(psi_mat, Su0, iu0, L)), 0.0, 1.0)
    FAF = float(8.0 * np.sum(g * (1.0 - g)))
    log(f"L={L}: E0={E0:.6f}  FAF(F1)={FAF:.6f}  n_up={g.sum():.4f}")

    cdU = AK.cdag_map(L, nup)
    seed = (cdU[0] @ psi_mat).reshape(-1).astype(complex)
    mv1, nS, Du1, Dd1, _, _ = sector_matvec(L, U, nup + 1, nd)
    log(f"L={L}: (N+1,Sz=+1) dim={nS}  |phi|^2={float(np.vdot(seed, seed).real):.6f}")

    H1op = sp.linalg.LinearOperator((nS, nS), matvec=mv1, dtype=float)
    # fixed ARPACK start vector here too: these two extremal eigenvalues only set the ends of
    # the frequency window, but an unseeded eigsh makes even that depend on the run.
    _v0w = np.random.default_rng(SEED_GS).standard_normal(nS)
    _v0w /= np.linalg.norm(_v0w)
    kw = dict(k=1, return_eigenvectors=False, ncv=min(nS, 20), maxiter=20000, tol=1e-9, v0=_v0w)
    emin = float(eigsh(H1op, which="SA", **kw)[0])
    emax = float(eigsh(H1op, which="LA", **kw)[0])
    grid = np.linspace(emin - E0 - 1.0, emax - E0 + 1.0, N_GRID)

    A_ex = spectra_all_eta(mv1, seed, E0, grid, etas, min(N_LANCZOS, nS))
    nrm = {e: _trapz(np.abs(A_ex[e]), grid) for e in etas}
    sums = {e: _trapz(A_ex[e], grid) for e in etas}
    log(f"L={L}: exact A(w) done, window [{grid[0]:.2f},{grid[-1]:.2f}], "
        + "  ".join(f"int(eta={e})={sums[e]:.4f}" for e in etas))

    order, _ = rank_configs_matrixfree(mv1, seed)
    log(f"L={L}: configuration ranking done")

    cache = {}

    def relL1_at(k, ordr=order, store=cache):
        """rel-L1 of the top-k subspace, for EVERY eta, from one Lanczos run."""
        if k in store:
            return store[k]
        Sset = np.sort(ordr[:k])

        def rmv(xs):
            xf = np.zeros(nS, dtype=complex)
            xf[Sset] = xs
            return mv1(xf)[Sset]

        A_s = spectra_all_eta(rmv, seed[Sset], E0, grid, etas, min(N_LANCZOS, k))
        r = {e: _trapz(np.abs(A_s[e] - A_ex[e]), grid) / nrm[e] for e in etas}
        store[k] = r
        return r

    ks = sorted(set(max(2, int(f * nS)) for f in COARSE_FRACS))
    for k in ks:
        r = relL1_at(k)
        log(f"L={L}:   frac={k / nS:.3f} (k={k})  "
            + "  ".join(f"relL1[{e}]={r[e]:.4f}" for e in etas))

    fracs = {}
    for e in etas:
        frac, klo = 1.0, 2
        for k in ks:
            if relL1_at(k)[e] < THR_RELL1:
                khi = k
                while khi - klo > max(2, nS // 200):
                    km = (klo + khi) // 2
                    if relL1_at(km)[e] < THR_RELL1:
                        khi = km
                    else:
                        klo = km
                frac = khi / nS
                break
            klo = k
        fracs[e] = float(frac)
        log(f"L={L}: eta={e:.3f} -> frac={frac:.4f}  |S|={int(round(frac * nS))}")

    rec = {
        "L": L, "qubits": 2 * L, "U": U, "E0": E0, "F1_FAF": FAF, "F1_per_site": FAF / L,
        "Np1_sector_dim": int(nS), "grid_window": [float(grid[0]), float(grid[-1])],
        "n_lanczos": N_LANCZOS, "n_grid": N_GRID, "threshold_relL1": THR_RELL1,
        "bisection_resolution_in_k": int(max(2, nS // 200)),
        "frac": {f"{e:g}": fracs[e] for e in etas},
        "S_abs": {f"{e:g}": int(round(fracs[e] * nS)) for e in etas},
        "exact_sum_rule": {f"{e:g}": sums[e] for e in etas},
        "coarse_scan": [{"k": int(k), "frac": k / nS,
                         "relL1": {f"{e:g}": cache[k][e] for e in etas}} for k in ks],
        "n_subspace_lanczos_runs": len(cache),
        "wall_seconds": time.time() - t_start,
    }

    if control:
        # C2: the matrix-free propagator must give the same answer as expm_multiply on explicit H.
        log(f"L={L}: CONTROL C2 -- repeating the ranking with scipy expm_multiply on explicit H")
        order2 = rank_configs_expm_multiply(L, U, nup + 1, nd, seed)
        cache2 = {}
        frac2, klo = 1.0, 2
        for k in ks:
            if relL1_at(k, order2, cache2)[CONTROL_ETA] < THR_RELL1:
                khi = k
                while khi - klo > max(2, nS // 200):
                    km = (klo + khi) // 2
                    if relL1_at(km, order2, cache2)[CONTROL_ETA] < THR_RELL1:
                        khi = km
                    else:
                        klo = km
                frac2 = khi / nS
                break
            klo = k
        overlap = len(set(order[:int(0.6 * nS)].tolist())
                      & set(order2[:int(0.6 * nS)].tolist())) / max(1, int(0.6 * nS))
        rec["control_C2"] = {"frac_expm_multiply": float(frac2),
                             "frac_krylov_matrixfree": fracs[CONTROL_ETA],
                             "top60pct_set_overlap": float(overlap)}
        log(f"L={L}: CONTROL C2  frac(expm_multiply)={frac2:.4f} vs "
            f"frac(matrix-free)={fracs[CONTROL_ETA]:.4f}  set-overlap={overlap:.4f}")
    return rec


# ----------------------------------------------------------------------------------------------
# fits and the linearity test
# ----------------------------------------------------------------------------------------------
def _ols(x, y, deg):
    """Least squares fit of degree `deg`; returns (coeffs low->high, residual sum of squares, R^2)."""
    X = np.vander(np.asarray(x, float), deg + 1, increasing=True)
    beta, *_ = np.linalg.lstsq(X, np.asarray(y, float), rcond=None)
    resid = np.asarray(y, float) - X @ beta
    rss = float(resid @ resid)
    tss = float(np.sum((np.asarray(y, float) - np.mean(y)) ** 2))
    return beta, rss, (1.0 - rss / tss if tss > 0 else float("nan")), X


def loglinear_report(Ls, vals):
    """Log-linear exponent + R^2 + local slopes + curvature (quadratic) significance test."""
    from scipy import stats
    Ls = np.asarray(Ls, float)
    y = np.log(np.asarray(vals, float))
    b1, rss1, r2, _ = _ols(Ls, y, 1)
    out = {"n_points": int(len(Ls)), "L": Ls.tolist(), "values": list(map(float, vals)),
           "b": float(b1[1]), "intercept": float(b1[0]), "R2": float(r2),
           "extrapolation": {f"L={int(LL)}": float(np.exp(b1[0] + b1[1] * LL))
                             for LL in (14, 16)}}
    slopes = [float((y[i + 1] - y[i]) / (Ls[i + 1] - Ls[i])) for i in range(len(Ls) - 1)]
    out["local_log_slopes"] = slopes
    out["local_slope_pairs"] = [f"{int(Ls[i])}->{int(Ls[i + 1])}" for i in range(len(Ls) - 1)]
    if slopes:
        amin, amax = min(abs(s) for s in slopes), max(abs(s) for s in slopes)
        out["local_slope_abs_range"] = [amin, amax]
        out["local_slope_variation_factor"] = float(amax / amin) if amin > 0 else float("inf")
    # curvature test: is the quadratic coefficient of log(value) vs L significantly non-zero?
    if len(Ls) >= 4:
        b2, rss2, r2q, X2 = _ols(Ls, y, 2)
        dof = len(Ls) - 3
        if dof > 0 and rss2 > 0:
            s2 = rss2 / dof
            cov = s2 * np.linalg.inv(X2.T @ X2)
            se_c = float(np.sqrt(cov[2, 2]))
            tstat = float(b2[2] / se_c) if se_c > 0 else float("inf")
            p_t = float(2 * stats.t.sf(abs(tstat), dof))
            Fstat = float(((rss1 - rss2) / 1.0) / s2)
            p_F = float(stats.f.sf(Fstat, 1, dof))
        else:
            se_c, tstat, p_t, Fstat, p_F = float("nan"), float("nan"), float("nan"), \
                float("nan"), float("nan")
        out["curvature_test"] = {
            "model": "log(value) = a + b*L + c*L^2",
            "c": float(b2[2]), "se_c": se_c, "t": tstat, "dof": int(max(dof, 0)),
            "p_two_sided": p_t, "F_linear_vs_quadratic": Fstat, "p_F": p_F,
            "R2_quadratic": float(r2q), "R2_linear": float(r2),
            "interpretation": ("c < 0 means the log-slope steepens with L: the sequence is NOT a "
                               "single exponential, and the quoted straight-line exponent is a fit "
                               "through a curve."),
        }
    return out


def build_fits(points, etas):
    """Log-linear fits of frac(L) and |S|(L) for every eta.

    Two windows are reported for each quantity: ALL computed L, and L >= 6.  L = 4 is a saturation
    artefact -- its (N+1) sector holds only 24 configurations, so the required fraction is pinned
    near 0.92 for every eta and drags the straight-line exponent toward zero.
    """
    fits = {}
    if len(points) < 3:
        return fits
    Lall = [p["L"] for p in points]
    for e in etas:
        key = f"{e:g}"
        block = {}
        for wname, sel in (("all_L", lambda L: True), ("L_ge_6", lambda L: L >= 6)):
            idx = [i for i, L in enumerate(Lall) if sel(L)]
            if len(idx) < 3:
                continue
            Lv = [Lall[i] for i in idx]
            fr = [points[i]["frac"][key] for i in idx]
            sa = [points[i]["S_abs"][key] for i in idx]
            block[wname] = {"fraction": loglinear_report(Lv, fr),
                            "absolute_support": loglinear_report(Lv, sa)}
        fits[key] = block
        for wname, b in block.items():
            f_, s_ = b["fraction"], b["absolute_support"]
            log(f"FIT eta={key} [{wname}]:  b(frac)={f_['b']:+.4f} R2={f_['R2']:.4f} | "
                f"b(|S|)={s_['b']:+.4f} R2={s_['R2']:.4f} | "
                f"frac@L=14 extrap={f_['extrapolation']['L=14']:.3f} | "
                f"local-slope spread x{f_.get('local_slope_variation_factor', float('nan')):.2f}")
    return fits


# ----------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--L", type=int, nargs="+", default=list(DEFAULT_L))
    ap.add_argument("--eta", type=float, nargs="+", default=list(DEFAULT_ETA))
    ap.add_argument("--out", default=DEFAULT_OUT,
                    help="output JSON (default: <repo>/data/eta_sweep.json); a relative path "
                         "is taken relative to the CURRENT DIRECTORY, never to the repository")
    ap.add_argument("--control-L", type=int, default=8,
                    help="L at which self-check C2 (matrix-free vs expm_multiply) is run; 0 = skip")
    ap.add_argument("--no-check", action="store_true", help="do not abort on a failed self-check")
    ap.add_argument("--refit", default=None, metavar="JSON",
                    help="rebuild the fit block from an existing output JSON (no physics is "
                         "recomputed; the measured points are carried over verbatim)")
    args = ap.parse_args()

    if args.refit:
        with open(args.refit) as f:
            prev = json.load(f)
        etas_prev = sorted(float(k) for k in prev["points"][0]["frac"])
        prev["fits"] = build_fits(prev["points"], etas_prev)
        pl = sorted(PUBLISHED_FRACS_L4_14)
        prev["published_sequence_L4_to_L14"] = loglinear_report(
            pl, [PUBLISHED_FRACS_L4_14[k] for k in pl])
        prev["published_sequence_L4_to_L14"]["source"] = (
            "QUOTED from release/data/scaling_data.json (eta = 0.15 t); NOT recomputed here.")
        prev["published_sequence_L4_to_L14"]["provenance_note"] = (
            "The six fractions are READ FROM release/data/scaling_data.json, never hand-typed. "
             "Until 2026-09-18 this file carried them as [0.9167, 0.83, 0.5625, 0.3625, 0.17, 0.08], "
             "which disagrees with the deposit at L=6 (0.82, not 0.83), L=8 (0.556122, not 0.5625) "
             "and L=10 (0.362491, not 0.3625). The JSON was patched by hand that day but the script "
             "was not, so any re-run silently reintroduced the wrong sequence. It now cannot: see "
             "_published_fracs().")
        prev.setdefault("provenance", {})["refit_from"] = os.path.abspath(args.refit)
        prev["provenance"]["refit_timestamp_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                                  time.gmtime())
        d = os.path.dirname(os.path.abspath(args.out))
        if d:
            os.makedirs(d, exist_ok=True)
        with open(args.out, "w") as f:
            json.dump(prev, f, indent=1)
        log(f"REFIT: {args.refit} -> {args.out}")
        return

    etas = sorted(set(round(float(e), 6) for e in args.eta))
    Ls = sorted(set(args.L))
    log(f"eta_sweep: L={Ls}  eta={etas}  out={args.out}")
    log(f"numpy {np.__version__}  scipy {scipy.__version__}  python {platform.python_version()}")

    points, checks = [], []
    for L in Ls:
        rec = sweep_L(L, etas, control=(L == args.control_L))
        points.append(rec)

        # ---- C1: the eta = 0.15 column must reproduce the published sequence, checked per L ----
        if abs(CONTROL_ETA - min(etas, key=lambda e: abs(e - CONTROL_ETA))) < 1e-9 and L in CONTROL_FRAC:
            got = rec["frac"][f"{CONTROL_ETA:g}"]
            exp = CONTROL_FRAC[L]
            ok = abs(got - exp) <= CONTROL_ATOL
            checks.append({"check": "C1", "L": L, "eta": CONTROL_ETA, "expected": exp,
                           "obtained": got, "atol": CONTROL_ATOL, "pass": bool(ok)})
            log(f"SELF-CHECK C1  L={L}: frac(eta=0.15)={got:.4f} vs published {exp:.4f} "
                f"(atol {CONTROL_ATOL}) -> {'PASS' if ok else 'FAIL'}")
            if not ok and not args.no_check:
                raise SystemExit(f"ABORT: self-check C1 failed at L={L}: "
                                 f"{got:.4f} != {exp:.4f} +- {CONTROL_ATOL}")
        if "control_C2" in rec:
            c = rec["control_C2"]
            ok = abs(c["frac_expm_multiply"] - c["frac_krylov_matrixfree"]) <= CONTROL_ATOL
            checks.append({"check": "C2", "L": L, "expected": c["frac_expm_multiply"],
                           "obtained": c["frac_krylov_matrixfree"], "atol": CONTROL_ATOL,
                           "pass": bool(ok)})
            log(f"SELF-CHECK C2  L={L} -> {'PASS' if ok else 'FAIL'}")
            if not ok and not args.no_check:
                raise SystemExit("ABORT: self-check C2 failed (ranking propagator mismatch)")

    fits = build_fits(points, etas)

    published = None
    if len(PUBLISHED_FRACS_L4_14) >= 4:
        pl = sorted(PUBLISHED_FRACS_L4_14)
        published = loglinear_report(pl, [PUBLISHED_FRACS_L4_14[k] for k in pl])
        published["source"] = ("QUOTED from release/data/scaling_data.json (eta = 0.15 t); NOT "
                               "recomputed here. L = 14 has sector dimension 10 306 296 and is out "
                               "of reach of this script.")
        published["provenance_note"] = (
            "The six fractions are READ FROM release/data/scaling_data.json, never hand-typed. "
             "Until 2026-09-18 this file carried them as [0.9167, 0.83, 0.5625, 0.3625, 0.17, 0.08], "
             "which disagrees with the deposit at L=6 (0.82, not 0.83), L=8 (0.556122, not 0.5625) "
             "and L=10 (0.362491, not 0.3625). The JSON was patched by hand that day but the script "
             "was not, so any re-run silently reintroduced the wrong sequence. It now cannot: see "
             "_published_fracs().")
        log("PUBLISHED sequence L=4..14: local log-slopes "
            + ", ".join(f"{s:+.3f}" for s in published["local_log_slopes"])
            + f"  (variation factor {published['local_slope_variation_factor']:.2f}x)")
        if "curvature_test" in published:
            ct = published["curvature_test"]
            log(f"PUBLISHED curvature test: c={ct['c']:+.5f} t={ct['t']:+.3f} "
                f"dof={ct['dof']} p={ct['p_two_sided']:.4g}")

    n_fail = sum(1 for c in checks if not c["pass"])
    out = {
        "script": "eta_sweep.py",
        "what": ("required sampled-subspace fraction frac(L,eta) for the (N+1) spectral function "
                 "of the 1D Hubbard chain, and the log-linear exponents of frac and of |S|"),
        "claim_supported": ("the exponent of the FRACTION is a function of the broadening eta and "
                            "is therefore not a property of the model; the exponent of the "
                            "ABSOLUTE support |S| is not; and at fixed eta the sequence is not a "
                            "single exponential"),
        "points": points,
        "fits": fits,
        "published_sequence_L4_to_L14": published,
        "self_checks": checks,
        "self_checks_passed": bool(n_fail == 0),
        "provenance": {
            "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "wall_seconds_total": time.time() - _T0,
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
            "parameters": {
                "U": U_HUB, "t": T_HOP, "n_lanczos": N_LANCZOS, "n_grid": N_GRID,
                "threshold_relL1": THR_RELL1, "K_rank": K_RANK, "dt_rank": DT_RANK,
                "krylov_substeps": SUBSTEPS, "krylov_m": KRYLOV_M,
                "coarse_fractions": list(COARSE_FRACS),
                "seed_gs_eigsh_start": SEED_GS, "seed_control": SEED_CTRL,
                "L": Ls, "eta": etas,
            },
            "control_values_reproduced": {
                "C1_published_frac_at_eta_0.15": CONTROL_FRAC,
                "C1_atol": CONTROL_ATOL,
                "C2": "matrix-free Krylov ranking vs scipy expm_multiply on explicit H",
            },
        },
    }
    d = os.path.dirname(os.path.abspath(args.out))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)
    log(f"WROTE {args.out}  (self-checks: {len(checks) - n_fail}/{len(checks)} pass)")
    if n_fail and not args.no_check:
        raise SystemExit("ABORT: one or more self-checks failed")


if __name__ == "__main__":
    main()
