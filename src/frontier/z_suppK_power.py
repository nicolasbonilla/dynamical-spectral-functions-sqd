# -*- coding: utf-8 -*-
"""Coordinate support of the Krylov space K(H, phi), by powers -- the estimator of Sec. V.E
for the sizes at which the projector onto K cannot be stored.

Every normalized power H^k phi / ||H^k phi|| lies in K = K(H, phi), so the union over
k <= kmax of the determinants on which it exceeds an amplitude threshold is contained in
the coordinate support of K.  It is a LOWER bound on that support, and a count of this
kind means nothing without its threshold, so every count is written with it.

What it establishes and what it does not (Sec. III and V.E): Lambda_S = 0 iff span(S)
contains K(H, phi_S), the Krylov space of the RETAINED part of the probe.  At full
captured weight (phi_S = phi) that is K(H, phi) itself, so a coordinate support equal to
the whole sector means that every proper coordinate subspace holding the whole probe
leaks.  Below full captured weight the space that decides is K(H, phi_S), which this
script does not measure.

Uses the deposited builder c0_lib.system(L) (the (N+1) sector of the half-filled ring,
U/t = 8, seed c^dag_{0,up}|Psi_0>) -- the same H and phi as every other frontier script.

    python src/frontier/z_suppK_power.py 8 10 12
writes data/c3_frontier/adversarial/z_suppK_power_L<L>.json for each L.
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


def _system14():
    """L = 14 from the ground-state checkpoint, exactly as run_L14.py builds it (the
    half-filled ground state of 11.8 million determinants is too slow to recompute here).
    The checkpoint is deposited in ckpt/ (94 MB); P2_CKPT overrides the folder."""
    from c0_big import sector_mv
    ck = os.environ.get("P2_CKPT") or os.path.join(ROOT, "ckpt")
    g = np.load(os.path.join(ck, "L14_gs.npz"))
    psi = g["psi"].astype(float)
    nup = nd = 7
    _, _, Du0, Dd0 = sector_mv(14, 8.0, nup, nd)
    phi = np.asarray((C.AK.cdag_map(14, nup)[0] @ psi.reshape(Du0, Dd0)).reshape(-1), dtype=float)
    mv, D, _, _ = sector_mv(14, 8.0, nup + 1, nd)
    return dict(U=8.0, phi=phi, D=D, H=None, mv=mv)


def run(L):
    t0 = time.time()
    s = _system14() if L == 14 else C.system(L)
    D = s["D"]
    H = s["H"] if s.get("H") is not None else None
    step = (lambda x: H @ x) if H is not None else s["mv"]
    v = s["phi"] / np.linalg.norm(s["phi"])
    mask = {t: np.abs(v) > t for t in THRESHOLDS}
    counts = {t: [int(mask[t].sum())] for t in THRESHOLDS}     # k = 0: the probe itself
    for _ in range(KMAX):
        v = step(v)
        v = v / np.linalg.norm(v)
        for t in THRESHOLDS:
            mask[t] |= np.abs(v) > t
            counts[t].append(int(mask[t].sum()))
    res = dict(
        L=L, U=s["U"], sector_dimension=D, probe_support=int(np.sum(np.abs(s["phi"]) > 0)),
        kmax=KMAX,
        union_count_by_k={"%g" % t: counts[t] for t in THRESHOLDS},
        at_k30={"%g" % t: counts[t][30] for t in THRESHOLDS},
        at_k36={"%g" % t: counts[t][36] for t in THRESHOLDS},
        whole_sector_at_k30={"%g" % t: counts[t][30] == D for t in THRESHOLDS},
        seconds=round(time.time() - t0, 1),
        note="lower bound on the coordinate support of K(H,phi); see module docstring",
    )
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "z_suppK_power_L%d.json" % L), "w") as f:
        json.dump(res, f, indent=1)
    print("L=%2d  D=%9d  k<=30: 1e-10 -> %d, 1e-12 -> %d   k<=36: 1e-10 -> %d   (%.0f s)"
          % (L, D, counts[1e-10][30], counts[1e-12][30], counts[1e-10][36], res["seconds"]))
    return res


if __name__ == "__main__":
    for L in [int(x) for x in sys.argv[1:]] or [8, 10]:
        run(L)
