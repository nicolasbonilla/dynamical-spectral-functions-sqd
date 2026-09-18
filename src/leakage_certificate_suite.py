# -*- coding: utf-8 -*-
r"""leakage_certificate_suite.py -- evaluates the Theorem-B certificate on the subspaces the paper ACTUALLY uses.

WHAT THIS COMPUTES
------------------
Four independent blocks, all built on src/leakage_certificate.py:

  teqsci    the 13 finite-shot TE-QSCI points of spectral_validation_max.py (L=6, U/t=8, eta=0.15,
            8 seeds, 4000 shots/step).  For every (seed, K) subspace it reports
            |S|, |S|/D, w_S, K_S, Lambda_S/eta, the lower bound, the upper bound, the true rel-L1
            and the slack.  CONTROL: the per-K mean |S| and mean rel-L1 must reproduce
            data/method_max.json bit for bit (the protocol is deterministic given the seeds).

  akw       the FR fractions of sampled_akw.py at L=6 and L=8 (U/t=8, eta=0.18), FR in
            {0.50,0.70,0.85,0.95}, all L momenta, both (N+1) and (N-1) branches.
            CONTROL: at L=8, FR=0.85 the per-k combined rel-L1 must reproduce
            data/sampled_akw_L8.json (mean 0.0021858473952153625, max 0.0032826276455945177).

  frontier  the UTILITY FRONTIER: Lambda_S(eta)/eta as a function of the retained sector fraction,
            on the exact ranking protocol of scaling_lanczos.py (seed c^dag_{0,up}|GS>, Born weight
            of the time-evolved state accumulated over K=18 steps of dt=0.5, eta=0.15), for
            L = 6, 8, 10.  The certificate says something non-vacuous only where Lambda_S/eta < 1.
            CONTROL: E0 and the published frac@rel-L1<0.05 of data/scaling_data.json.

  stress    a fresh falsification battery for the two-sided bound: subspaces from Hamiltonian and
            selection families that are NOT Hubbard-chain / not top-Born, and eta values outside
            the paper's range.  Any single violation is the headline result.

WHICH PAPER CLAIM THIS SUPPORTS
-------------------------------
Everything here supports the REPLACEMENT of Theorem 1(iii) by Theorem B (see src/leakage_certificate.py).  It
also produces, for the first time in this repository, the numbers w_S and K_S that the published
"genuine certificate" is written in terms of but that no script ever evaluated.

HOW TO RUN
----------
    python src/leakage_certificate_suite.py selftest
    python src/leakage_certificate_suite.py teqsci   [outdir]      #  ~2.5 min
    python src/leakage_certificate_suite.py akw      [outdir]      #  ~22 min (two dense eigh of 3920)
    python src/leakage_certificate_suite.py frontier [outdir] [Ls] #  L=6,8 ~4.5 min; L=10 ~16 min, 0.9 GB
    python src/leakage_certificate_suite.py stress   [outdir]      #  ~8.5 min
    python src/leakage_certificate_suite.py all      [outdir]

outdir defaults to <repo>/data/; every block writes <outdir>/cert_<block>.json with a
"provenance" section (parameters, library versions, runtime, controls reproduced).

The Hubbard sector engine is imported from src/akw_lanczos.py and the reference JSONs are read
from <repo>/data/, both resolved from __file__ -- the script runs from any working directory and
contains no container path.  P2_SRC / P2_DATA override the two locations if ever needed.

All ARPACK start vectors are fixed (ARPACK_SEED below), so every block is bit-for-bit
reproducible; without that the stress-battery slacks moved between runs.

MEASURED 2026-09-18 by running every block from the repository root (numpy 1.26.4 /
scipy 1.13.1, 16-core laptop, other jobs competing):
selftest 5 s, teqsci 138 s, akw 1339 s (22 min), frontier(6,8) 264 s, frontier(10) 973 s
(16 min), stress 513 s (8.5 min).  Peak RSS ~0.9 GB (frontier L=10).
Every block reproduces its control and is bit-for-bit reproducible: the stress block was run
three independent times and all 5730 JSON leaves were identical.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

# --- NUMPY_TRAPEZOID_BRIDGE ------------------------------------------------
# numpy 2.0 ADDED np.trapezoid and REMOVED np.trapz.  Files in this repository use
# both names, so without this bridge no single numpy version runs the whole deposit.
if not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz          # numpy < 2.0
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid          # numpy >= 2.0
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
if os.environ.get("P2_SRC"):
    sys.path.insert(0, os.environ["P2_SRC"])

import akw_lanczos as AK                                    # noqa: E402
from leakage_certificate import (certify, selftest, provenance,   # noqa: E402
                                 spectral_from_pairs, _TRAPZ)

# Repository layout: this file lives in <repo>/src/, the deposit in <repo>/data/.
# NO container paths: the default output directory is resolved from __file__.
REPO = os.path.dirname(_HERE)
DEFAULT_OUTDIR = os.path.join(REPO, "data")

# FIXED 2026-09-18.  Every eigsh() below used to start from the RANDOM vector scipy draws for
# itself, so the last ~3 digits of E0 -- and with them every rel-L1, every Lambda and the whole
# stress-battery slack table -- moved from run to run.  Two runs of the stress block differed in
# min_slack (0.0 vs 4.7e-40) and in median_slack (1.8188 vs 1.8317) for no reason but that.
# A deposit that does not reproduce bit for bit is not a deposit, so the start vector is now
# drawn once from a fixed seed.
ARPACK_SEED = 20260918


def arpack_v0(n):
    """Deterministic ARPACK start vector -> bit-reproducible eigenpairs."""
    v = np.random.default_rng(ARPACK_SEED).standard_normal(n)
    return v / np.linalg.norm(v)

T0 = time.time()


def log(*a):
    print("[%7.1fs]" % (time.time() - T0), *a, flush=True)


def _dump(obj, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, name)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=float)
    log("wrote", path)
    return path


# ==================================================================================
# BLOCK 1 -- the 13 TE-QSCI points of spectral_validation_max.py
# ==================================================================================
def _build_L6_fock(L=6, U=8.0, t_hop=1.0):
    """Verbatim Fock-space builder of spectral_validation_max.py (OBC chain)."""
    M = 2 * L
    dim = 1 << M

    def c_op(p):
        r, c, d = [], [], []
        for s in range(dim):
            if (s >> p) & 1:
                sg = (-1) ** bin(s & ((1 << p) - 1)).count('1')
                r.append(s & ~(1 << p)); c.append(s); d.append(float(sg))
        return sp.csr_matrix((d, (r, c)), shape=(dim, dim))

    C = [c_op(p) for p in range(M)]
    Cd = [c.T.conj() for c in C]
    Nn = np.array([bin(s).count('1') for s in range(dim)])
    Sz = np.array([sum(((s >> (2 * i)) & 1) - ((s >> (2 * i + 1)) & 1) for i in range(L))
                   for s in range(dim)])
    H = sp.csr_matrix((dim, dim))
    for i in range(L - 1):
        for spn in (0, 1):
            a = 2 * i + spn; b = 2 * (i + 1) + spn
            H = H - t_hop * (Cd[a] @ C[b] + Cd[b] @ C[a])
    for i in range(L):
        H = H + U * (Cd[2 * i] @ C[2 * i]) @ (Cd[2 * i + 1] @ C[2 * i + 1])
    return H.tocsr(), C, Cd, Nn, Sz, M, dim


def run_teqsci(outdir=".", L=6, U=8.0, K=12, dt=0.4, shots=4000, eta=0.15, NSEEDS=8,
               ref_json=None):
    t_blk = time.time()
    log("TEQSCI: building the L=%d Fock space ..." % L)
    H, C, Cd, Nn, Sz, M, dim = _build_L6_fock(L, U)
    gi = np.where((Nn == L) & (Sz == 0))[0]
    Hg = H[gi][:, gi].toarray()
    w, v = np.linalg.eigh(Hg)
    E0 = float(w[0])
    psi0 = np.zeros(dim); psi0[gi] = v[:, 0]
    phi_full = Cd[0] @ psi0
    si = np.where((Nn == L + 1) & (Sz == 1))[0]
    H7 = sp.csr_matrix(H[si][:, si])
    En, Vn = np.linalg.eigh(H7.toarray())
    phi = phi_full[si]
    wnorm = float(np.vdot(phi, phi).real)
    amp = Vn.conj().T @ phi
    poles = En - E0
    wts = np.abs(amp) ** 2
    log("TEQSCI: E0=%.12f sumrule=%.6f dim(N+1,Sz+1)=%d" % (E0, wnorm, len(si)))

    # --- the paper's own window recipe, needed to reproduce method_max.json exactly ---
    pgrid = np.linspace(poles.min() - 1, poles.max() + 1, 600)
    A_ex_w = spectral_from_pairs(poles, wts, pgrid, eta)
    normA = float(_TRAPZ(A_ex_w, pgrid))

    def relL1_paper(Sset):
        Ss = np.array(sorted(Sset))
        Em, Um = np.linalg.eigh(H7[Ss][:, Ss].toarray())
        aS = Um.conj().T @ phi[Ss]
        A_S = spectral_from_pairs(Em - E0, np.abs(aS) ** 2, pgrid, eta)
        return float(_TRAPZ(np.abs(A_S - A_ex_w), pgrid) / normA)

    def evolve(tk):
        return Vn @ (np.exp(-1j * En * tk) * (Vn.conj().T @ phi))

    rows = []
    per_seed_S = np.zeros((NSEEDS, K + 1))
    per_seed_l1 = np.zeros((NSEEDS, K + 1))
    for sd in range(NSEEDS):
        rng = np.random.default_rng(1000 + sd)                 # verbatim seed of the paper script
        seen = {int(np.argmax(np.abs(phi) ** 2))}
        for k in range(K + 1):
            vk = evolve(k * dt)
            pr = np.abs(vk) ** 2; pr /= pr.sum()
            for d in np.unique(rng.choice(len(si), size=shots, p=pr)):
                seen.add(int(d))
            Ss = np.array(sorted(seen))
            per_seed_S[sd, k] = len(Ss)
            per_seed_l1[sd, k] = relL1_paper(seen)
            cert = certify(H7, phi, Ss, eta, E0=E0, ref_poles=En, ref_weights=wts,
                           dense_max=10 ** 6)
            cert.update(seed=int(1000 + sd), K=int(k), relL1_paper=per_seed_l1[sd, k])
            rows.append(cert)
        log("TEQSCI: seed %d done (|S|_final=%d)" % (1000 + sd, per_seed_S[sd, -1]))

    Smean = per_seed_S.mean(0); l1mean = per_seed_l1.mean(0); l1std = per_seed_l1.std(0)

    # ---- CONTROL: reproduce data/method_max.json -----------------------------------
    controls = {}
    if ref_json and os.path.exists(ref_json):
        ref = json.load(open(ref_json))
        dS = max(abs(Smean[r['K']] - r['S']) for r in ref['sampled'])
        dl = max(abs(l1mean[r['K']] - r['l1']) for r in ref['sampled'])
        dE = abs(E0 - ref['E0']); dsr = abs(wnorm - ref['sumrule'])
        assert dS < 1e-9, "CONTROL FAILED: mean |S| differs from method_max.json by %g" % dS
        assert dl < 1e-9, "CONTROL FAILED: mean rel-L1 differs from method_max.json by %g" % dl
        assert dE < 1e-10 and dsr < 1e-10, "CONTROL FAILED: E0/sum rule"
        controls = dict(source=os.path.basename(ref_json), max_abs_dev_meanS=dS,
                        max_abs_dev_meanRelL1=dl, dev_E0=dE, dev_sumrule=dsr,
                        tol="1e-9 absolute")
        log("TEQSCI CONTROL OK: method_max.json reproduced (max dev |S| %.2e, rel-L1 %.2e)" % (dS, dl))
    else:
        log("TEQSCI: reference json not found -- control SKIPPED")

    # per-K aggregate over the 8 seeds
    agg = []
    for k in range(K + 1):
        sub = [r for r in rows if r['K'] == k]
        agg.append(dict(K=k,
                        S_mean=float(Smean[k]), frac_mean=float(Smean[k] / len(si)),
                        w_S_mean=float(np.mean([r['w_S'] for r in sub])),
                        w_S_min=float(np.min([r['w_S'] for r in sub])),
                        K_S=int(sub[0]['K_S']),
                        K_S_all_equal=bool(len(set(r['K_S'] for r in sub)) == 1),
                        Lambda_over_eta_mean=float(np.mean([r['Lambda_over_eta'] for r in sub])),
                        Lambda_hat_over_eta_mean=float(np.mean([r['Lambda_hat_over_eta'] for r in sub])),
                        rel_lower_mean=float(np.mean([r['rel_lower'] for r in sub])),
                        rel_upper_mean=float(np.mean([r['rel_upper'] for r in sub])),
                        relL1_true_mean=float(np.mean([r['relL1_true'] for r in sub])),
                        relL1_paper_mean=float(l1mean[k]), relL1_paper_std=float(l1std[k]),
                        slack_mean=float(np.mean([r['upper'] / r['L1_true'] for r in sub])),
                        slack_min=float(np.min([r['upper'] / r['L1_true'] for r in sub])),
                        n_viol_lower=int(sum(not r['lower_ok'] for r in sub)),
                        n_viol_upper=int(sum(not r['upper_ok'] for r in sub))))
    nv = sum(a['n_viol_lower'] + a['n_viol_upper'] for a in agg)
    log("TEQSCI: %d subspaces certified, %d bound violations" % (len(rows), nv))
    out = dict(block="teqsci", L=L, U=U, eta=eta, K=K, dt=dt, shots=shots, NSEEDS=NSEEDS,
               dim_sector=int(len(si)), E0=E0, sumrule=wnorm,
               n_subspaces=len(rows), n_violations=int(nv), per_K=agg, rows=rows,
               provenance=provenance(dict(L=L, U=U, eta=eta, K=K, dt=dt, shots=shots,
                                          NSEEDS=NSEEDS, seeds="np.random.default_rng(1000+sd)"),
                                     controls, t_blk, script=os.path.basename(__file__)))
    return out


# ==================================================================================
# BLOCK 2 -- the FR fractions of sampled_akw.py
# ==================================================================================
def run_akw(outdir=".", Ls=(6, 8), U=8.0, eta=0.18, K=16, dt=0.5,
            FRs=(0.50, 0.70, 0.85, 0.95), dense_max=400, n_kry=500, ref_json=None):
    t_blk = time.time()
    out_L = {}
    controls = {}
    rows = []
    for L in Ls:
        nup = nd = L // 2
        Hexpl, Du, Dd = AK.build_H_explicit(L, U, nup, nd)
        E0v, V0 = eigsh(Hexpl, k=1, which='SA', v0=arpack_v0(Hexpl.shape[0]),
                        ncv=min(Hexpl.shape[0], 40), tol=1e-11, maxiter=100000)
        E0 = float(E0v[0]); Psi = V0[:, 0].reshape(Du, Dd)
        wg = np.linspace(-9, 9, 600)                     # the figure's own window

        cdU = AK.cdag_map(L, nup); cUr = AK.cdag_map(L, nup - 1)
        Hadd, _, _ = AK.build_H_explicit(L, U, nup + 1, nd)
        Hrem, _, _ = AK.build_H_explicit(L, U, nup - 1, nd)
        log("AKW L=%d: E0=%.12f  dim(add)=%d dim(rem)=%d -- dense eigh ..."
            % (L, E0, Hadd.shape[0], Hrem.shape[0]))
        Ea, Ua = np.linalg.eigh(Hadd.toarray())
        Er, Ur = np.linalg.eigh(Hrem.toarray())
        log("AKW L=%d: eigh done" % L)

        ks = [2 * np.pi * n / L for n in range(L)]
        per_k = {}
        for n, k in enumerate(ks):
            ph = np.exp(1j * k * np.arange(L)) / np.sqrt(L)
            add = np.zeros((len(AK.strings(L, nup + 1)[0]), Dd), dtype=complex)
            for j in range(L):
                add += ph[j] * (cdU[j] @ Psi)
            rem = np.zeros((len(AK.strings(L, nup - 1)[0]), Dd), dtype=complex)
            for j in range(L):
                rem += np.conj(ph[j]) * (cUr[j].T @ Psi)
            seeds = {'add': (add.ravel(), Hadd, Ea, Ua, +1.0), 'rem': (rem.ravel(), Hrem, Er, Ur, -1.0)}

            # ranking by accumulated Born weight of the time-evolved state (verbatim protocol)
            rank = {}
            A_ex_tot = np.zeros_like(wg)
            for tag, (sd, Hs, Em, Um, sgn) in seeds.items():
                coef = Um.conj().T @ sd
                A_ex_tot += spectral_from_pairs(sgn * (Em - E0), np.abs(coef) ** 2, wg, eta)
                wc = np.zeros(Hs.shape[0])
                for kk in range(K + 1):
                    vk = Um @ (np.exp(-1j * Em * kk * dt) * coef)
                    wc += np.abs(vk) ** 2
                rank[tag] = np.argsort(wc)[::-1]
            nrm_w = float(_TRAPZ(np.abs(A_ex_tot), wg))

            for FR in FRs:
                A_S_tot = np.zeros_like(wg)
                bound_tot = 0.0; norm2_tot = 0.0; lo_tot = 0.0
                for tag, (sd, Hs, Em, Um, sgn) in seeds.items():
                    nSfull = Hs.shape[0]
                    kk = max(1, int(round(FR * nSfull)))
                    Sset = np.sort(rank[tag][:kk])
                    cert = certify(Hs, sd, Sset, eta, E0=E0, ref_poles=Em,
                                   ref_weights=np.abs(Um.conj().T @ sd) ** 2,
                                   dense_max=dense_max, n_kry=n_kry,
                                   kry_check=(n == 0 and FR == 0.85), return_pairs=True)
                    th = cert.pop("_ritz_theta"); rw = cert.pop("_ritz_weight")
                    cert.update(L=L, branch=tag, FR=FR, k_over_pi=2.0 * n / L)
                    rows.append(cert)
                    bound_tot += cert['upper']; norm2_tot += cert['norm_phi'] ** 2
                    lo_tot += cert['lower']
                    # the figure's own combined A_S on its own window, from the SAME Ritz pairs
                    A_S_tot += spectral_from_pairs(sgn * (th - E0), rw, wg, eta)
                r_fig = float(_TRAPZ(np.abs(A_S_tot - A_ex_tot), wg) / nrm_w)
                per_k.setdefault("%.3f" % (2 * n / L), {})["%.2f" % FR] = dict(
                    relL1_figure_window=r_fig, combined_upper=bound_tot, combined_lower=lo_tot,
                    combined_rel_upper=bound_tot / norm2_tot,
                    combined_rel_lower=lo_tot / norm2_tot)
            log("AKW L=%d: k/pi=%.3f done" % (L, 2 * n / L))
        out_L[str(L)] = dict(E0=E0, dim_add=int(Hadd.shape[0]), dim_rem=int(Hrem.shape[0]),
                             per_k=per_k)

    # ---- CONTROL: reproduce data/sampled_akw_L8.json at FR=0.85 ---------------------
    if ref_json and os.path.exists(ref_json) and "8" in out_L:
        ref = json.load(open(ref_json))
        mine = {kk: vv["0.85"]["relL1_figure_window"] for kk, vv in out_L["8"]["per_k"].items()}
        dev = max(abs(mine[kk] - ref['per_k'][kk]['relL1']) for kk in ref['per_k'])
        mmean = float(np.mean(list(mine.values())))
        # tolerance: A_S is built from a 500-step Krylov compression, not from the dense
        # eigendecomposition the original script used; 5e-5 absolute on a quantity of 2.2e-3.
        assert dev < 5e-5, "CONTROL FAILED: sampled_akw_L8.json per-k rel-L1 dev %g" % dev
        controls = dict(source=os.path.basename(ref_json), max_abs_dev_per_k=dev,
                        mean_relL1_mine=mmean, mean_relL1_published=ref['mean_relL1'],
                        tol="5e-5 absolute")
        log("AKW CONTROL OK: mean rel-L1 %.10f vs published %.10f (max per-k dev %.2e)"
            % (mmean, ref['mean_relL1'], dev))
    else:
        log("AKW: reference json not found -- control SKIPPED")

    nv = sum((not r['lower_ok']) + (not r['upper_ok']) for r in rows)
    log("AKW: %d channel subspaces certified, %d bound violations" % (len(rows), nv))
    return dict(block="akw", U=U, eta=eta, K=K, dt=dt, FRs=list(FRs), per_L=out_L,
                n_subspaces=len(rows), n_violations=int(nv), rows=rows,
                provenance=provenance(dict(Ls=list(Ls), U=U, eta=eta, K=K, dt=dt, FRs=list(FRs),
                                           dense_max=dense_max, n_kry=n_kry),
                                      controls, t_blk, script=os.path.basename(__file__)))


# ==================================================================================
# BLOCK 3 -- the utility frontier  Lambda_S/eta  vs  retained fraction
# ==================================================================================
def _krylov_expm(mv, v, dt, m=8):
    """e^{-i dt H} v via a small Lanczos basis (memory = m complex vectors)."""
    from scipy.linalg import expm as dense_expm
    n = v.shape[0]
    beta = float(np.linalg.norm(v))
    if beta < 1e-14:
        return v.copy()
    V = np.empty((m, n), dtype=complex)
    V[0] = v / beta
    al = np.zeros(m); be = np.zeros(m); mm = m
    w = mv(V[0]); a = np.vdot(V[0], w).real; al[0] = a; w = w - a * V[0]
    for j in range(1, m):
        b = float(np.linalg.norm(w)); be[j - 1] = b
        if b < 1e-12:
            mm = j; break
        V[j] = w / b
        w = mv(V[j]); a = np.vdot(V[j], w).real; al[j] = a
        w = w - a * V[j] - b * V[j - 1]
    T = np.diag(al[:mm]) + np.diag(be[:mm - 1], 1) + np.diag(be[:mm - 1], -1)
    E = dense_expm(-1j * dt * T)[:, 0]
    return beta * (V[:mm].T @ E)


def run_frontier(outdir=".", Ls=(6, 8), U=8.0, eta=0.15, K=18, dt=0.5,
                 fracs=(0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80,
                        0.85, 0.90, 0.95, 0.98, 1.00),
                 n_kry=260, ref_json=None):
    t_blk = time.time()
    curves = {}
    controls = {}
    for L in Ls:
        nup = nd = L // 2
        H0, Du0, Dd0 = AK.build_H_explicit(L, U, nup, nd)
        wv, vv = eigsh(H0, k=1, which='SA', tol=1e-12, v0=arpack_v0(H0.shape[0]),
                       ncv=min(H0.shape[0], 40), maxiter=100000)
        E0 = float(wv[0]); psi_mat = vv[:, 0].reshape(Du0, Dd0)
        cdU = AK.cdag_map(L, nup)
        H1, Du1, Dd1 = AK.build_H_explicit(L, U, nup + 1, nd)
        nS = H1.shape[0]
        seed = (cdU[0] @ psi_mat).reshape(-1).astype(complex)
        log("FRONTIER L=%d: E0=%.12f  (N+1) dim=%d  ||phi||^2=%.6f"
            % (L, E0, nS, float(np.vdot(seed, seed).real)))

        # ranking: accumulated Born weight of the time-evolved state (scaling_lanczos.py protocol)
        vt = seed / np.sqrt(np.vdot(seed, seed).real)
        wc = np.abs(vt) ** 2
        for _ in range(1, K + 1):
            vt = _krylov_expm(lambda x: H1.dot(x), vt, dt, m=10)
            wc += np.abs(vt) ** 2
        order = np.argsort(wc)[::-1]
        del vt, wc
        log("FRONTIER L=%d: ranking done" % L)

        # exact reference: Lehmann pairs if dense is affordable, else Haydock on the figure window
        dense_ref = nS <= 4200
        if dense_ref:
            Em, Um = np.linalg.eigh(H1.toarray())
            ref_poles, ref_w = Em, np.abs(Um.conj().T @ seed) ** 2
            del Em, Um
        else:
            ref_poles = ref_w = None
            _v0 = arpack_v0(H1.shape[0])
            emin = float(eigsh(H1, k=1, which='SA', return_eigenvectors=False, tol=1e-8,
                               v0=_v0, maxiter=100000)[0])
            emax = float(eigsh(H1, k=1, which='LA', return_eigenvectors=False, tol=1e-8,
                               v0=_v0, maxiter=100000)[0])
            grid = np.linspace(emin - E0 - 1.0, emax - E0 + 1.0, 600)
            A_ex = (-AK.haydock(lambda x: H1.dot(x) - E0 * x, seed, min(260, nS),
                                grid + 1j * eta).imag / np.pi)
            nrm = float(_TRAPZ(np.abs(A_ex), grid))

        pts = []
        for f in fracs:
            kk = max(2, int(round(f * nS)))
            Sset = np.sort(order[:kk])
            cert = certify(H1, seed, Sset, eta, E0=E0, ref_poles=ref_poles, ref_weights=ref_w,
                           dense_max=1200, n_kry=n_kry)
            if not dense_ref:
                def rmv(xs, Sset=Sset):
                    xf = np.zeros(nS, dtype=complex); xf[Sset] = xs
                    return H1.dot(xf)[Sset]
                A_s = (-AK.haydock(lambda x: rmv(x) - E0 * x, seed[Sset], min(260, kk),
                                   grid + 1j * eta).imag / np.pi)
                cert['relL1_figure_window'] = float(_TRAPZ(np.abs(A_s - A_ex), grid) / nrm)
            cert.update(L=L, frac_target=f)
            pts.append(cert)
            log("FRONTIER L=%d frac=%.3f |S|=%d  w_S=%.6f  Lam/eta=%.4g  relUB=%.4g  relL1=%s"
                % (L, kk / nS, kk, cert['w_S'], cert['Lambda_hat_over_eta'], cert['rel_upper'],
                   ("%.4g" % cert['relL1_true']) if 'relL1_true' in cert
                   else ("%.4g" % cert.get('relL1_figure_window', float('nan')))))
        curves[str(L)] = dict(E0=E0, dim=int(nS), points=pts,
                              ref="dense-Lehmann" if dense_ref else "Haydock nl=260, figure window")
        if ref_json and os.path.exists(ref_json):
            ref = json.load(open(ref_json))
            p = [q for q in ref['points'] if q['L'] == L]
            if p:
                dE = abs(E0 - p[0]['E0'])
                assert dE < 1e-8, "CONTROL FAILED: E0(L=%d) dev %g vs scaling_data.json" % (L, dE)
                controls["E0_L%d" % L] = dict(mine=E0, published=p[0]['E0'], dev=dE,
                                              dim_mine=int(nS), dim_published=p[0]['Np1_sector'],
                                              published_frac=p[0]['frac'], tol="1e-8 absolute")
                log("FRONTIER CONTROL OK: E0(L=%d) reproduces scaling_data.json (dev %.2e)" % (L, dE))
    return dict(block="frontier", U=U, eta=eta, K=K, dt=dt, curves=curves,
                provenance=provenance(dict(Ls=list(Ls), U=U, eta=eta, K=K, dt=dt,
                                           fracs=list(fracs), n_kry=n_kry), controls, t_blk,
                                      script=os.path.basename(__file__)))


# ==================================================================================
# BLOCK 4 -- fresh falsification battery (families and eta NOT used in the earlier sweep)
# ==================================================================================
def _H_heisenberg_xxz(N, Delta, Jz_seed=None):
    """XXZ chain, PBC, Sz=0 sector.  Completely different model class from the Hubbard chain."""
    from itertools import combinations
    states = [sum(1 << b for b in c) for c in combinations(range(N), N // 2)]
    states.sort()
    idx = {m: i for i, m in enumerate(states)}
    D = len(states)
    rows, cols, val = [], [], []
    for a, m in enumerate(states):
        diag = 0.0
        for i in range(N):
            j = (i + 1) % N
            si = 1 if (m >> i) & 1 else -1
            sj = 1 if (m >> j) & 1 else -1
            diag += 0.25 * Delta * si * sj
            if si != sj:
                m2 = m ^ (1 << i) ^ (1 << j)
                rows.append(idx[m2]); cols.append(a); val.append(0.5)
        rows.append(a); cols.append(a); val.append(diag)
    return sp.csr_matrix((val, (rows, cols)), shape=(D, D)), states


def _H_tV_ring(L, N, V, t=1.0):
    """Spinless fermions with nearest-neighbour repulsion on a ring (t-V model)."""
    from itertools import combinations
    states = [sum(1 << b for b in c) for c in combinations(range(L), N)]
    states.sort()
    idx = {m: i for i, m in enumerate(states)}
    D = len(states)
    rows, cols, val = [], [], []
    for a, m in enumerate(states):
        diag = 0.0
        for i in range(L):
            j = (i + 1) % L
            if ((m >> i) & 1) and ((m >> j) & 1):
                diag += V
            for (p, q) in ((i, j), (j, i)):
                if (m >> q) & 1 and not (m >> p) & 1:
                    m2 = (m & ~(1 << q)) | (1 << p)
                    lo, hi = min(p, q), max(p, q)
                    mask = m & (((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))
                    sgn = -1.0 if (bin(mask).count('1') & 1) else 1.0
                    rows.append(idx[m2]); cols.append(a); val.append(-t * sgn)
        rows.append(a); cols.append(a); val.append(diag)
    return sp.csr_matrix((val, (rows, cols)), shape=(D, D)), states


def _H_anderson(D, W, seed):
    """1D Anderson chain with box disorder (single particle, open boundaries)."""
    r = np.random.default_rng(seed)
    diag = r.uniform(-W, W, D)
    off = -np.ones(D - 1)
    return sp.diags([off, diag, off], [-1, 0, 1], format='csr')


def _H_sparse_random(D, deg, seed):
    """Sparse symmetric random matrix, average degree `deg` (no lattice structure at all)."""
    r = np.random.default_rng(seed)
    nnz = int(deg * D)
    rows = r.integers(0, D, nnz)
    cols = r.integers(0, D, nnz)
    vals = r.normal(size=nnz)
    M = sp.coo_matrix((vals, (rows, cols)), shape=(D, D)).tocsr()
    M = (M + M.T) * 0.5
    return (M + sp.diags(r.normal(size=D))).tocsr()


def _select(kind, H, phi, frac, seed):
    """Subspace selection families, none of which is the paper's top-Born ranking."""
    D = H.shape[0]
    n = max(2, int(round(frac * D)))
    r = np.random.default_rng(seed)
    if kind == "top_phi":
        return np.sort(np.argsort(np.abs(phi))[::-1][:n])
    if kind == "random":
        return np.sort(r.choice(D, n, replace=False))
    if kind == "antiphi":                      # adversarial: drop the LARGEST amplitudes
        return np.sort(np.argsort(np.abs(phi))[:n])
    if kind == "diag_window":                  # configurations closest in diagonal energy to <phi|H|phi>
        d = np.asarray(H.diagonal())
        e0 = float(np.vdot(phi, H @ phi).real / np.vdot(phi, phi).real)
        return np.sort(np.argsort(np.abs(d - e0))[:n])
    if kind == "half_phi_half_random":
        a = np.argsort(np.abs(phi))[::-1][:n // 2]
        rest = np.setdiff1d(np.arange(D), a)
        b = r.choice(rest, n - len(a), replace=False)
        return np.sort(np.concatenate([a, b]))
    if kind == "block":                        # a contiguous block of configurations
        s = int(r.integers(0, D - n + 1))
        return np.arange(s, s + n)
    raise ValueError(kind)


def run_stress(outdir="."):
    t_blk = time.time()
    rows = []
    fams = []

    # --- family A: XXZ chain, N=12 spins, Sz=0 (D=924); seed = S^+_0 |GS> projected back ---
    for Delta in (0.5, 1.0, 2.0):
        H, states = _H_heisenberg_xxz(12, Delta)
        ev, V = np.linalg.eigh(H.toarray())
        gs = V[:, 0]
        # seed: S^z_0 |GS| (stays in the same Sz sector) -- a genuine dynamical structure-factor seed
        sz0 = np.array([0.5 if (m >> 0) & 1 else -0.5 for m in states])
        phi = sz0 * gs
        fams.append(("xxz_D%.1f" % Delta, H, phi, ev, np.abs(V.conj().T @ phi) ** 2))

    # --- family B: t-V ring, L=12, N=6 (D=924) ---
    for V_ in (1.0, 4.0):
        H, states = _H_tV_ring(12, 6, V_)
        ev, Vv = np.linalg.eigh(H.toarray())
        gs = Vv[:, 0]
        n0 = np.array([1.0 if (m >> 0) & 1 else 0.0 for m in states])
        phi = (n0 - n0.mean()) * gs
        fams.append(("tV_L12_V%.0f" % V_, H, phi, ev, np.abs(Vv.conj().T @ phi) ** 2))

    # --- family C: Anderson chain, D=600 ---
    for W in (1.0, 5.0):
        H = _H_anderson(600, W, seed=31415 + int(W))
        ev, V = np.linalg.eigh(H.toarray())
        phi = np.zeros(600); phi[300] = 1.0
        fams.append(("anderson_W%.0f" % W, H, phi, ev, np.abs(V.conj().T @ phi) ** 2))

    # --- family D: structureless sparse random matrix, D=800 ---
    for deg in (6, 20):
        H = _H_sparse_random(800, deg, seed=2718 + deg)
        ev, V = np.linalg.eigh(H.toarray())
        r = np.random.default_rng(99 + deg)
        phi = r.normal(size=800)
        fams.append(("sprand_deg%d" % deg, H, phi, ev, np.abs(V.conj().T @ phi) ** 2))

    kinds = ["top_phi", "random", "antiphi", "diag_window", "half_phi_half_random", "block"]
    etas = [0.03, 0.07, 0.22, 0.40, 0.90]            # deliberately outside the paper's 0.10-0.18
    fr = [0.10, 0.35, 0.75]
    idx = 0
    for fname, H, phi, ev, wts in fams:
        for kind in kinds:
            for f in fr:
                idx += 1
                eta = etas[idx % len(etas)]
                S = _select(kind, H, phi, f, seed=7000 + idx)
                cert = certify(H, phi, S, eta, E0=0.0, ref_poles=ev, ref_weights=wts,
                               dense_max=1200, n_kry=300)
                cert.update(family=fname, selection=kind, frac_target=f, case=idx)
                rows.append(cert)
        log("STRESS: family %-16s done (%d cases so far)" % (fname, len(rows)))

    viol = [r for r in rows if not (r['lower_ok'] and r['upper_ok'])]
    slacks = np.array([r['upper'] / r['L1_true'] for r in rows])
    # The raw minimum is NOT the tightest case: it is dominated by rows where BOTH sides are
    # numerically zero (a strongly localised state whose subspace already contains everything,
    # e.g. bound 8.8e-55 against a true error of 1.9e-15).  A ratio of two machine zeros says
    # nothing about how tight the bound is, so the honest statistic is the minimum restricted to
    # rows whose true error is actually resolved.  Both are reported; neither is hidden.
    TRUE_FLOOR = 1e-10
    resolved = np.array([r['upper'] / r['L1_true'] for r in rows
                         if r['relL1_true'] > TRUE_FLOOR])
    min_nd = float(resolved.min()) if resolved.size else float('nan')
    # And the two ramps must never be mixed.  The upper bound is a MINIMUM of two expressions:
    # the trivial ||phi||^2 (1 + w_S) and the Theorem-B leakage expression.  A slack of exactly
    # 1 is attained on the TRIVIAL branch at w_S = 0, where A_S vanishes identically and the
    # bound is an identity -- that says nothing about the leakage branch.  Report them apart.
    leak = np.array([r['upper'] / r['L1_true'] for r in rows
                     if r['relL1_true'] > TRUE_FLOOR
                     and r['upper_unclipped'] < r['upper_trivial']])
    min_leak = float(leak.min()) if leak.size else float('nan')
    log("STRESS: %d NEW subspaces, %d violations; min slack %.4g raw / %.4f over the %d rows "
        "with a resolved true error / %.4f over the %d rows where the LEAKAGE branch binds; "
        "median slack %.2f"
        % (len(rows), len(viol), slacks.min(), min_nd, resolved.size,
           min_leak, leak.size, np.median(slacks)))
    return dict(block="stress", n_subspaces=len(rows), n_violations=len(viol),
                violations=viol, min_slack=float(slacks.min()),
                min_slack_note=("the raw minimum is a ratio of two numerical zeros (both bound "
                                "and true error below 1e-12) and is NOT a tightness statement; "
                                "use min_slack_resolved"),
                min_slack_resolved=min_nd,
                n_rows_resolved=int(resolved.size),
                min_slack_resolved_leakage_branch=min_leak,
                n_rows_leakage_branch=int(leak.size),
                branch_note=("min_slack_resolved is attained on the TRIVIAL branch "
                             "||phi||^2 (1+w_S), where it is an identity at w_S = 0. The only "
                             "figure that says anything about Theorem B is "
                             "min_slack_resolved_leakage_branch, restricted to the rows where "
                             "the leakage expression is the one that binds."),
                true_error_floor=TRUE_FLOOR,
                median_slack=float(np.median(slacks)),
                families=sorted(set(r['family'] for r in rows)),
                selections=kinds, etas=etas, fracs=fr, rows=rows,
                provenance=provenance(dict(families=[f[0] for f in fams], selections=kinds,
                                           etas=etas, fracs=fr, rng="default_rng(7000+case)"),
                                      dict(note="novelty by construction: no Hubbard-chain "
                                                "Hamiltonian and no top-Born selection here"),
                                      t_blk, script=os.path.basename(__file__)))


# ==================================================================================
if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "selftest"
    outdir = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_OUTDIR
    DATA = os.environ.get("P2_DATA", os.path.join(REPO, "data"))
    mm = os.path.join(DATA, "method_max.json")
    sa = os.path.join(DATA, "sampled_akw_L8.json")
    sd = os.path.join(DATA, "scaling_data.json")

    log("SELF-VERIFICATION of src/leakage_certificate.py (mandatory) ...")
    selftest()
    log("SELF-VERIFICATION PASSED")

    if cmd in ("teqsci", "all"):
        _dump(run_teqsci(outdir, ref_json=mm), outdir, "cert_teqsci.json")
    if cmd in ("akw", "all"):
        _dump(run_akw(outdir, ref_json=sa), outdir, "cert_akw.json")
    if cmd in ("frontier", "all"):
        Ls = tuple(int(x) for x in sys.argv[3].split(",")) if len(sys.argv) > 3 else (6, 8)
        _dump(run_frontier(outdir, Ls=Ls, ref_json=sd), outdir,
              "cert_frontier_%s.json" % "-".join(str(x) for x in Ls))
    if cmd in ("stress", "all"):
        _dump(run_stress(outdir), outdir, "cert_stress.json")
    log("DONE")
