# -*- coding: utf-8 -*-
"""Coordinate support of K(H, phi), exact in the reflection symmetry -- the corrected estimator
of Sec. V.E.

Why this file exists.  The reflection R: j -> -j about the seed site commutes with H, and the
seed phi = c^dag_{0,up}|Psi_0> is an eigenvector of it (parity s = +1 or -1).  Every vector of
K(H, phi) therefore vanishes EXACTLY on the R-fixed determinants of parity -s: 4 of 300 at L = 6,
18 of 3920 at L = 8, 48 of 52 920 at L = 10.  Repeated multiplication by H amplifies round-off in
that forbidden sector until it crosses any amplitude threshold, which is why the unprojected
estimator (z_suppK_power.py) and the dense projector both reported "the whole sector" -- a
floating-point artefact.  Here every power is projected back onto the parity sector of phi, so
the forbidden determinants stay at zero and what is counted is the true support.

Checks, each fatal: [R, H] = 0 on a random vector; R phi = s phi.

    python src/frontier/z_suppK_sym.py 6 8 10 12          (L = 12: about a minute)
    python src/frontier/z_suppK_sym.py 14        (reads ckpt/L14_gs.npz; P2_CKPT overrides)
writes data/c3_frontier/adversarial/z_suppK_sym_L<L>.json
"""
from __future__ import print_function
import json
import os
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import c0_lib as C  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(ROOT, "data", "c3_frontier", "adversarial")
THRESHOLDS = (1e-10, 1e-12)
KMAX = 36


def reflection_species(L, n):
    """R on single-species strings (sorted by integer, as akw_lanczos.strings):
    c^dag_{o1}...c^dag_{on}|0> -> c^dag_{-o1}...c^dag_{-on}|0>, reordered ascending."""
    S, idx = C.AK.strings(L, n)
    perm = np.zeros(len(S), dtype=np.int64)
    sign = np.zeros(len(S))
    for a, m in enumerate(S):
        img = [(-i) % L for i in range(L) if (m >> i) & 1]
        arr, sg = list(img), 1.0
        for i in range(len(arr)):
            for j in range(len(arr) - 1 - i):
                if arr[j] > arr[j + 1]:
                    arr[j], arr[j + 1] = arr[j + 1], arr[j]
                    sg = -sg
        perm[a] = idx[sum(1 << i for i in img)]
        sign[a] = sg
    return perm, sign


def reflection(L, nup, ndn):
    pu, gu = reflection_species(L, nup)
    pd, gd = reflection_species(L, ndn)
    Dd = len(pd)
    P = (pu[:, None] * Dd + pd[None, :]).ravel()
    G = (gu[:, None] * gd[None, :]).ravel()

    def R(v):
        out = np.empty_like(v)
        out[P] = G * v
        return out
    fixed = np.nonzero(P == np.arange(P.size))[0]
    return R, fixed, G[fixed]


def system(L):
    if L != 14:
        s = C.system(L)
        H = s["H"]
        return s["phi"], (lambda x: H @ x), s["D"]
    from c0_big import sector_mv
    ck = os.environ.get("P2_CKPT") or os.path.join(ROOT, "ckpt")
    psi = np.load(os.path.join(ck, "L14_gs.npz"))["psi"].astype(float)
    _, _, Du0, Dd0 = sector_mv(14, 8.0, 7, 7)
    phi = np.asarray((C.AK.cdag_map(14, 7)[0] @ psi.reshape(Du0, Dd0)).reshape(-1), dtype=float)
    mv, D, _, _ = sector_mv(14, 8.0, 8, 7)
    return phi, mv, D


def run(L):
    t0 = time.time()
    phi, mv, D = system(L)
    R, fixed, fixed_sign = reflection(L, L // 2 + 1, L // 2)
    x = np.random.default_rng(7).standard_normal(D)
    comm = np.linalg.norm(R(mv(x)) - mv(R(x))) / np.linalg.norm(mv(x))
    if comm > 1e-12:
        raise SystemExit("L=%d: [R,H] = %.1e != 0 -- the reflection does not commute here" % (L, comm))
    s = float(np.sign(phi @ R(phi)))
    par = np.linalg.norm(R(phi) - s * phi) / np.linalg.norm(phi)
    # The ground state is non-degenerate, hence exactly reflection-symmetric; what is tested
    # is the eigensolver's accuracy.  The L = 14 checkpoint has residual 2.7e-9 (run_L14.py),
    # which puts the parity defect of its seed at 1.6e-9, so the tolerance is 1e-7; the
    # projection below then removes that component.
    if par > 1e-7:
        raise SystemExit("L=%d: phi is not a reflection eigenvector (defect %.1e)" % (L, par))
    forbidden = int(np.sum(fixed_sign == -s))
    proj = lambda v: 0.5 * (v + s * R(v))
    v = proj(phi / np.linalg.norm(phi))
    mask = {t: np.abs(v) > t for t in THRESHOLDS}
    counts = {t: [int(mask[t].sum())] for t in THRESHOLDS}
    for _ in range(KMAX):
        v = proj(mv(v))
        v = v / np.linalg.norm(v)
        for t in THRESHOLDS:
            mask[t] |= np.abs(v) > t
            counts[t].append(int(mask[t].sum()))
    allowed = D - forbidden
    res = dict(L=L, sector_dimension=D, reflection_parity_of_phi=int(s), commutator=comm,
               parity_defect=par, reflection_fixed=int(fixed.size), forbidden=forbidden,
               allowed=allowed, kmax=KMAX,
               union_count_by_k={"%g" % t: counts[t] for t in THRESHOLDS},
               at_k30={"%g" % t: counts[t][30] for t in THRESHOLDS},
               at_k36={"%g" % t: counts[t][36] for t in THRESHOLDS},
               all_allowed_at_k36={"%g" % t: counts[t][36] == allowed for t in THRESHOLDS},
               seconds=round(time.time() - t0, 1))
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "z_suppK_sym_L%d.json" % L), "w") as f:
        json.dump(res, f, indent=1)
    print("L=%2d D=%9d parity %+d forbidden %d allowed %d | k<=36: 1e-10 -> %d, 1e-12 -> %d | "
          "k<=30: 1e-10 -> %d  [comm %.0e, %.0f s]"
          % (L, D, s, forbidden, allowed, counts[1e-10][36], counts[1e-12][36], counts[1e-10][30],
             comm, res["seconds"]))


if __name__ == "__main__":
    for L in [int(a) for a in sys.argv[1:]] or [6, 8, 10]:
        run(L)
