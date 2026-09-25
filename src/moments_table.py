# -*- coding: utf-8 -*-
"""Table tab:moments (Sec. IV), its negative control, and the two-level series of
App. app:moments:conj -- regenerated from nothing but this file.

Replaces the four lost artefacts (final_out.txt, bc2_out.txt, control_out.txt,
toy_out.txt) that verify.py declared as a coverage gap.  The Hamiltonian is built here
from scratch, independently of src/frontier/c0_lib.py, so the table is also a
cross-check of that builder.

Open Hubbard chain, L = 6, U/t = 8, half filling; seed phi = c^dag_{0,up} |Psi_0> in the
(N+1) sector (4 up, 3 down), dimension D = 300; S_K = union_{a<=K} supp(H^a phi).
Moments mu_j = <phi|(H - E0)^j|phi>, reconstructed moments on span(S) with the same shift.

What is reproducible exactly and what is not.  |S|, |S|/D, w_S, K_S and the order-(2K+2)
discrepancies are properties of the construction and must agree with the table.  The
entries max_{j<=2K+1} eps_j are floating-point residues (~1e-15); the table quotes the
largest over six implementations, so this single implementation must give a value of the
same order and no more is claimed.

    python src/moments_table.py        ->  data/moments_table.json
"""
from __future__ import print_function
import json
import os
from itertools import combinations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "moments_table.json")


def strings(L, n):
    S = sorted(sum(1 << b for b in c) for c in combinations(range(L), n))
    return S, {m: i for i, m in enumerate(S)}


def hop(L, n, pbc, t=1.0):
    S, idx = strings(L, n)
    r, c, v = [], [], []
    bonds = [(i, (i + 1) % L) for i in range(L if pbc else L - 1)]
    for a, m in enumerate(S):
        for (i, j) in bonds:
            for (p, q) in ((i, j), (j, i)):
                if (m >> q) & 1 and not (m >> p) & 1:
                    m2 = (m & ~(1 << q)) | (1 << p)
                    lo, hi = min(p, q), max(p, q)
                    mask = m & (((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))
                    s = -1.0 if bin(mask).count("1") & 1 else 1.0
                    r.append(idx[m2]); c.append(a); v.append(-t * s)
    return sp.csr_matrix((v, (r, c)), shape=(len(S), len(S))), S


def hamiltonian(L, U, nu, nd, pbc):
    Tu, Su = hop(L, nu, pbc)
    Td, Sd = hop(L, nd, pbc)
    uo = np.array([[(m >> i) & 1 for i in range(L)] for m in Su], float)
    do = np.array([[(m >> i) & 1 for i in range(L)] for m in Sd], float)
    H = sp.kron(Tu, sp.identity(len(Sd))) + sp.kron(sp.identity(len(Su)), Td) \
        + sp.diags(U * (uo @ do.T).ravel())
    return H.tocsr(), Su, Sd


def cdag_up(L, n, j):
    Sa, _ = strings(L, n)
    _, ib = strings(L, n + 1)
    r, c, v = [], [], []
    for a, m in enumerate(Sa):
        if not (m >> j) & 1:
            s = -1.0 if bin(m & ((1 << j) - 1)).count("1") & 1 else 1.0
            r.append(ib[m | (1 << j)]); c.append(a); v.append(s)
    return sp.csr_matrix((v, (r, c)), shape=(len(ib), len(Sa)))


def moments(Hs, v, jmax):
    out, x = [], v.copy()
    for _ in range(jmax + 1):
        out.append(float(v @ x))
        x = Hs @ x
    return out


def table_rows(L=6, U=8.0, pbc=False, jmax=7, tol=1e-12):
    H0, Su, Sd = hamiltonian(L, U, L // 2, L // 2, pbc)
    E0, V = eigsh(H0, k=1, which="SA", v0=np.random.default_rng(1).standard_normal(H0.shape[0]))
    E0 = float(E0[0])
    Psi = V[:, 0].reshape(len(Su), len(Sd))
    phi = np.asarray((cdag_up(L, L // 2, 0) @ Psi)).ravel()
    phi = phi if np.sum(phi) >= 0 else -phi
    H1 = hamiltonian(L, U, L // 2 + 1, L // 2, pbc)[0].toarray()
    D = H1.shape[0]
    Hs = H1 - E0 * np.eye(D)
    mu = moments(Hs, phi, jmax)
    krylov = [phi]
    for _ in range(4):
        krylov.append(H1 @ krylov[-1])
    rows, S = [], set()
    for K in range(3):
        S |= set(np.nonzero(np.abs(krylov[K]) > tol * np.linalg.norm(krylov[K]))[0].tolist())
        idx = np.array(sorted(S))
        out = np.setdiff1d(np.arange(D), idx)
        KS = -1
        for a, x in enumerate(krylov):
            if np.linalg.norm(x[out]) <= 1e-12 * np.linalg.norm(x):
                KS = a
            else:
                break
        pS = phi[idx]
        muS = moments(Hs[np.ix_(idx, idx)], pS, jmax)
        rel = [abs(mu[j] - muS[j]) / abs(mu[j]) for j in range(jmax + 1)]
        ab = [abs(mu[j] - muS[j]) for j in range(jmax + 1)]
        full = len(idx) == D
        rows.append(dict(
            K=K, S=len(idx), frac=round(len(idx) / D, 3), w_S=float(pS @ pS / (phi @ phi)),
            K_S=("full sector" if full else KS),
            max_rel_eps_upto_2K1=max(rel[:2 * K + 2]),
            rel_eps_2K2=(None if full else rel[2 * K + 2]),
            abs_eps_2K2=(None if full else ab[2 * K + 2])))
    return dict(D=D, E0=E0, rows=rows), Hs, phi


def negative_control(Hs, phi, size=200, draws=300, seed=20260918):
    rng = np.random.default_rng(seed)
    D = len(phi)
    mu0 = float(phi @ phi)
    fails = 0
    worst_pass = 0.0
    for _ in range(draws):
        idx = np.sort(rng.choice(D, size=size, replace=False))
        e0 = abs(mu0 - float(phi[idx] @ phi[idx])) / mu0
        if e0 > 1e-10:
            fails += 1
        else:
            worst_pass = max(worst_pass, e0)
    return dict(size=size, draws=draws, seed=seed, fail_at_zeroth_moment=fails)


def two_level_series(eta=0.15, gs=(1, 3, 10, 30, 100, 300, 1000)):
    """||A - A_{K_0}||_1 for H = [[0,g],[g,0]], phi = e0: two poles at +-g (weight 1/2 each)
    against one pole at 0 (weight 1); Lorentzians of half width eta, integrated adaptively
    over the WHOLE real line.  A finite window drops the Lorentzian tails and understates
    the norm by up to 1e-3 at g = 1000; that is why the values printed before 2026-09-25
    were low from g = 30 on."""
    from scipy.integrate import quad
    lor = lambda x: (eta / np.pi) / (x * x + eta * eta)
    out = []
    for g in gs:
        f = lambda w: abs(0.5 * lor(w - g) + 0.5 * lor(w + g) - lor(w))
        edges = [-np.inf, -g, -g / 2.0, 0.0, g / 2.0, g, np.inf]
        tot = sum(quad(f, a, b, limit=500, epsabs=1e-13, epsrel=1e-12)[0]
                  for a, b in zip(edges, edges[1:]))
        out.append(dict(g=g, L1=float(tot)))
    return out


if __name__ == "__main__":
    tab, Hs, phi = table_rows()
    res = dict(table=tab, negative_control=negative_control(Hs, phi),
               two_level=two_level_series())
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=1)
    for r in tab["rows"]:
        print("K=%d |S|=%d %.3f w=%.12f K_S=%s  max eps<=2K+1 %.1e  eps_2K+2 %s (abs %s)"
              % (r["K"], r["S"], r["frac"], r["w_S"], r["K_S"], r["max_rel_eps_upto_2K1"],
                 "%.2e" % r["rel_eps_2K2"] if r["rel_eps_2K2"] is not None else "-",
                 "%.3e" % r["abs_eps_2K2"] if r["abs_eps_2K2"] is not None else "-"))
    nc = res["negative_control"]
    print("control: %d of %d random size-%d subspaces fail at the zeroth moment"
          % (nc["fail_at_zeroth_moment"], nc["draws"], nc["size"]))
    print("two-level eta=0.15:", " ".join("%.4f" % x["L1"] for x in res["two_level"]))
