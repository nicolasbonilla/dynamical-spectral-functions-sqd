# -*- coding: utf-8 -*-
r"""Re-run of the K=18 ranking protocol at L=14, against the stored ranking ckpt/L14_order.npy.

Added 2026-09-28.  The L=14 certificate point (run_L14.py, Table III, Sec. S4) selects its
subspace from the STORED ranking ckpt/L14_order.npy, which had been verified only to be a
permutation of the (N+1) sector; the ranking protocol itself had been re-run at L=10 (fproto.py,
set overlap 1.000000) but not at L=14.  This script re-runs it at L=14 with the same code path
fproto.py uses at L=10 -- fcore.ranking(): K=18 Born snapshots, dt=0.5, sub=16, Krylov m=6,
the accumulated Born weight of the time-evolved seed -- on the seed c^dag_{0 up}|Psi_0> built
from the stored ground state ckpt/L14_gs.npz exactly as run_L14.py builds it, and compares the
top-k selection with the stored one AS SETS at the published fraction (k = 824504) and at a few
other cuts, with the tie census of fcore.subspace() at each cut.

Run:  python src/frontier/rerank_L14.py   -> data/c3_frontier/published_fraction/L14_rerank.json
(about an hour single-threaded; ~2 GB of memory).  Needs ckpt/ (or P2_CKPT).
"""
import os, sys, time, gc, json, hashlib, datetime
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
ROOT = os.path.dirname(SRC)
OUT = os.environ.get('C3_OUT') or os.path.join(ROOT, 'data', 'c3_frontier')
PUB = os.path.join(OUT, 'published_fraction')
sys.dont_write_bytecode = True
sys.path.insert(0, HERE)
import fcore as F                                    # noqa: E402
from fcore import log                                # noqa: E402
from c0_big import sector_mv                         # noqa: E402
from c0_lib import AK                                # noqa: E402

CK = os.environ.get('P2_CKPT') or os.path.join(ROOT, 'ckpt')
L = 14; U = 8.0; nup = nd = 7
K_PUB = 824504                                       # |S| of the L=14 certificate point
_gs = os.path.join(CK, 'L14_gs.npz'); _ord = os.path.join(CK, 'L14_order.npy')
if not (os.path.isfile(_gs) and os.path.isfile(_ord)):
    print("SKIP: ckpt/L14_gs.npz and ckpt/L14_order.npy are needed (set P2_CKPT)."); sys.exit(3)


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


g = np.load(_gs)
E0 = float(g['E0']); psi = g['psi'].astype(float)
mv0, D0, Du0, Dd0 = sector_mv(L, U, nup, nd)
Psi = psi.reshape(Du0, Dd0); cdU = AK.cdag_map(L, nup)
phi = np.asarray((cdU[0] @ Psi).reshape(-1), dtype=float)
del Psi, cdU, psi, g, mv0; gc.collect()
mv, D, Du, Dd = sector_mv(L, U, nup + 1, nd)
log("L=%d D=%d E0=%.10f |phi|^2=%.12f  [stored ground state]" % (L, D, E0, float(phi @ phi)))

order_new, wc, t_rank = F.ranking(dict(mv=mv, phi=phi, D=D), K=18)
log("ranking re-run (K=18, dt=0.5, sub=16, m=6): %.0fs" % t_rank)
order_old = np.load(_ord)
assert order_old.size == D and order_new.size == D

rows = []
for k in (K_PUB, int(round(0.04 * D)), int(round(0.16 * D)), int(round(0.32 * D))):
    a = np.zeros(D, bool); a[order_old[:k]] = True
    b = np.zeros(D, bool); b[order_new[:k]] = True
    inter = int(np.count_nonzero(a & b))
    _, diag = F.subspace(order_new, wc, D, k / D)
    diag = {kk: (v if not isinstance(v, (np.floating, np.integer)) else v.item()) for kk, v in diag.items()}
    # the cut of subspace() is round(FR*D); re-derive with the exact k to be safe
    wsel = wc[order_new[:k]]; cut = float(wc[order_new[k - 1]])
    n_at_cut_tot = int(np.count_nonzero(wc == cut)); n_at_cut_in = int(np.count_nonzero(wsel == cut))
    rows.append(dict(k=int(k), FR=k / D, set_overlap=inter / k, n_common=inter,
                     n_only_stored=int(k - inter),
                     tied_at_cut_inside=n_at_cut_in, tied_at_cut_total=n_at_cut_tot,
                     cut_ambiguous=bool(n_at_cut_tot > n_at_cut_in),
                     w_S_stored=float(np.sum(phi[order_old[:k]] ** 2) / float(phi @ phi)),
                     w_S_rerun=float(np.sum(phi[order_new[:k]] ** 2) / float(phi @ phi))))
    log("k=%d (FR %.5f): set overlap %.6f (%d not shared); ties at cut %d/%d; w_S %.9f vs %.9f" % (
        k, k / D, inter / k, k - inter, n_at_cut_in, n_at_cut_tot, rows[-1]['w_S_stored'], rows[-1]['w_S_rerun']))

with open(os.path.abspath(__file__), 'rb') as f:
    gen = hashlib.sha256(f.read().replace(b'\r\n', b'\n')).hexdigest()
out = dict(what="K=18 ranking protocol re-run at L=14 vs the stored ranking ckpt/L14_order.npy, set overlap of the "
                "top-k selections", L=L, U=U, D=int(D), protocol=dict(K=18, dt=0.5, sub=16, m=6),
           rows=rows, t_rank_s=round(t_rank, 1),
           provenance=dict(script="src/frontier/rerank_L14.py", generator_sha256_lf=gen,
                           L14_gs_sha256=sha256(_gs), L14_order_sha256=sha256(_ord),
                           generated_utc=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")))
os.makedirs(PUB, exist_ok=True)
json.dump(out, open(os.path.join(PUB, 'L14_rerank.json'), 'w'), indent=1)
log("WROTE L14_rerank.json")
