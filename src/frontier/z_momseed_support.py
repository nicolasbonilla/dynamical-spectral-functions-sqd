# -*- coding: utf-8 -*-
r"""z_momseed_support.py -- coordinate support of the momentum seeds of Fig. akwsampled.

    phi_k = L^{-1/2} sum_j e^{ikj} c^dag_{j,up} |Psi_0>,   k = 2 pi m / L,

on the half-filled Hubbard ring (PBC, U/t = 8) at L = 6 and 8, in the (N+1, Sz=+1) sector.
Sec. III quotes the L = 8 range ("support on 3848 to 3920 of the 3920 determinants") to say
that the support floor of the LOCAL seed does not apply to the momentum sweep.

The support is counted above 1e-10 and above 1e-8 of the largest |phi_k| amplitude; the two
counts must agree for every k (the amplitudes split across a gap, so the count is not a
choice of threshold) or the script exits non-zero and writes nothing.  The ground state is an
ARPACK eigsh from a fixed-seed start vector, as in src/support_witness.py.

    python src/frontier/z_momseed_support.py      # -> data/c3_frontier/adversarial/z_momseed_support.json

Runtime: a few seconds.  No QPU.
"""
import json
import os
import sys

import numpy as np
from scipy.sparse.linalg import eigsh

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
REPO = os.path.dirname(SRC)
sys.path.insert(0, SRC)
import akw_lanczos as AK  # noqa: E402

U, T, SEED = 8.0, 1.0, 20260918
OUT = os.path.join(REPO, "data", "c3_frontier", "adversarial", "z_momseed_support.json")


def run(L):
    nup = nd = L // 2
    H0, Du, Dd = AK.build_H_explicit(L, U, nup, nd, T)
    v0 = np.random.default_rng(SEED).standard_normal(Du * Dd)
    _, V = eigsh(H0, k=1, which="SA", v0=v0 / np.linalg.norm(v0), tol=1e-12, maxiter=100000)
    Psi = V[:, 0].reshape(Du, Dd)
    cd = AK.cdag_map(L, nup)
    ph = [np.asarray((cd[j] @ Psi).reshape(-1)) for j in range(L)]
    D = int(ph[0].size)
    rows = []
    for m in range(L):
        k = 2 * np.pi * m / L
        v = sum(np.exp(1j * k * j) * ph[j] for j in range(L)) / np.sqrt(L)
        a = np.abs(v)
        n10 = int((a > 1e-10 * a.max()).sum())
        n8 = int((a > 1e-8 * a.max()).sum())
        if n10 != n8:
            sys.exit("L=%d m=%d: support %d above 1e-10 but %d above 1e-8 -- no gap, "
                     "nothing written" % (L, m, n10, n8))
        rows.append(dict(m=m, k_over_pi=2.0 * m / L, support=n10))
    s = [r["support"] for r in rows]
    print("L=%d  D=%d  support per m: %s  (min %d, max %d)" % (L, D, s, min(s), max(s)))
    return dict(L=L, D=D, rows=rows, support_min=min(s), support_max=max(s))


def main():
    res = [run(L) for L in (6, 8)]
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(dict(generated_by="src/frontier/z_momseed_support.py",
                       model="1D Hubbard ring, PBC, half filling, U/t=8",
                       seed="phi_k = L^-1/2 sum_j e^{ikj} c^dag_{j,up}|Psi_0>",
                       thresholds="1e-10 and 1e-8 of max|phi_k|, required to agree",
                       arpack_seed=SEED, results=res), fh, indent=1)
    print("WROTE", os.path.relpath(OUT, REPO))


if __name__ == "__main__":
    main()
