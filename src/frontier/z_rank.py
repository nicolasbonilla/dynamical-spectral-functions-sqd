# -*- coding: utf-8 -*-
"""ATK2 / the decisive test of 'whose property is the frontier?'

If the certificate frontier sits near 1 because of the SAMPLING PROTOCOL (the Born ranking
picks the wrong determinants), then an ORACLE ranking should push it far down.  If it sits
near 1 because of the GEOMETRY of the Krylov space in the determinant basis, no ranking helps.
Three rankings at L=6, same certificate:
   (a) repo protocol  : accumulated Born weight of the 19 time snapshots   (what the device does)
   (b) |phi_i|^2      : Born weight of the probe alone
   (c) (P_K)_ii       : the diagonal of the projector onto the Krylov space K(H,phi) -- a pure
                        ORACLE: the only determinants that matter for A(w), ranked by how much
                        of K they carry.  No sampler can compute this.
"""
import os, sys, time, json, numpy as np
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
# --- deposit paths (repaired 2026-09-19; see docs/C3_FRONTIER.md) ----------
SRC = os.path.dirname(HERE)                       # release/src
ROOT = os.path.dirname(SRC)                       # repository root
OUT = os.environ.get('C3_OUT') or os.path.join(ROOT, 'data', 'c3_frontier')
TMP = os.environ.get('C3_TMP') or os.path.join(ROOT, 'build')
# --------------------------------------------------------------------------
ADV = os.path.join(OUT, 'adversarial')
if not os.path.isdir(ADV):
    os.makedirs(ADV)
PUERTAS = HERE
sys.path.insert(0, PUERTAS)
import c0_lib as C
import z_fine as Z

L = int(sys.argv[1]); etas = [float(x) for x in sys.argv[2].split(',')]
s = C.system(L); H = s['H']; phi = s['phi']; E0 = s['E0']; D = s['D']
Em, Um = np.linalg.eigh(H.toarray()); c = Um.T@phi; nph = np.linalg.norm(phi)
order_repo, wc = Z.ranking(H, phi)
# Krylov projector diagonal
live = np.abs(c) > 1e-12*nph
E = Em[live]; Uv = Um[:, live]; cv = c[live]; cols = []
i = 0
while i < len(E):
    j = i
    while j+1 < len(E) and E[j+1]-E[i] < 1e-9: j += 1
    v = Uv[:, i:j+1]@cv[i:j+1]; n = np.linalg.norm(v)
    if n > 1e-12*nph: cols.append(v/n)
    i = j+1
K = np.array(cols).T; diagK = (K*K).sum(1)
RANK = dict(repo=order_repo, phi=np.argsort(np.abs(phi)**2)[::-1], krylov=np.argsort(diagK)[::-1])
print("L=%d D=%d dimK=%d  |supp K|=%d" % (L, D, K.shape[1], int((diagK > 1e-14).sum())))
print("%-8s %-7s %-8s %-14s %-13s %-13s %-10s" % ("ranking", "|S|", "FR", "w_S", "rel err R0", "ThmB R5", "trivial"))
ks = [int(round(f*D)) for f in (0.30, 0.50, 0.56, 0.70, 0.85, 0.95, 0.98)]
res = []
for name, order in RANK.items():
    for k in ks:
        Sset = np.sort(order[:k])
        rows = Z.evaluate(L, Em, Um, H, phi, E0, Sset, etas)
        for r in rows:
            r['ranking'] = name; res.append(r)
            print("%-8s %-7d %-8.4f %-14.12f %-13.4e %-13.4e %-10.4f%s" %
                  (name, k, k/D, r['w_S'], r['relR0'], r['relR5'], r['trivial_rel'],
                   "   <-- NON-VACUOUS" if r['relR5'] < r['trivial_rel'] else ""))
    print()
json.dump(res, open(os.path.join(ADV, 'z_rank_L%d.json' % L), 'w'), indent=1)
