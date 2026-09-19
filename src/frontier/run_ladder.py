# -*- coding: utf-8 -*-
"""C6: run the exact 5-rung ladder over (L, FR, eta).  Usage: run_ladder.py L FR1,FR2,... """
import os, sys, gc, json, time, numpy as np
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
# --- deposit paths (repaired 2026-09-19; see docs/C3_FRONTIER.md) ----------
SRC = os.path.dirname(HERE)                       # release/src
ROOT = os.path.dirname(SRC)                       # repository root
OUT = os.environ.get('C3_OUT') or os.path.join(ROOT, 'data', 'c3_frontier')
TMP = os.environ.get('C3_TMP') or os.path.join(ROOT, 'build')
# --------------------------------------------------------------------------
import fron_lib as F
from fron_lib import log

L = int(sys.argv[1])
FRS = [float(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else [0.10,0.30,0.50,0.70,0.85,0.95]
ETAS = [float(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else [0.05,0.10,0.15,0.18,0.25,0.50]
PER_ETA = float(os.environ.get('PER_ETA', '12'))

s = F.build(L); H = s['H']; phi = s['phi']; E0 = s['E0']; D = s['D']
log("L=%d D=%d E0=%.12f |phi|^2=%.6f" % (L, D, E0, phi @ phi))
t = time.time(); order, wc = F.repo_ranking(H, phi, K=18, dt=0.5, sub=16)
log("repo ranking (K=18,dt=0.5,sub=16 -> 19 snapshots) done in %.1fs" % (time.time()-t))

# ---- condition 5 bookkeeping: zero-amplitude determinants & argsort ties
supp = np.where(np.abs(phi) > 0)[0]
log("supp(phi) = %d/%d = %.6f   (1/2+1/L = %.6f)" % (len(supp), D, len(supp)/D, 0.5+1.0/L))
wcs = np.sort(wc)[::-1]

t = time.time(); Em, Um = np.linalg.eigh(H.toarray()); log("dense eigh(H) %.1fs" % (time.time()-t))

res = dict(L=L, D=int(D), E0=float(E0), norm_phi2=float(phi@phi), supp=int(len(supp)),
           protocol=dict(K=18, dt=0.5, sub=16, snapshots=19, ranking="born_order_mf on time-evolved seed"),
           per_eta=PER_ETA, cases=[])
cases = [('FR=%.4f' % fr, np.sort(order[:max(1, int(round(fr*D)))])) for fr in FRS]
cases.append(('S*=supp(phi)', np.sort(supp)))
for tag, Sset in cases:
    k = len(Sset)
    nzero = int(np.sum(np.abs(phi[Sset]) == 0.0))
    # exact tie at the cut of the wc ranking?
    tie = float('nan')
    if k < D:
        tie = float(wcs[k-1] - wcs[k])
    t = time.time()
    R = F.ladder(Em, Um, phi, E0, H, Sset, ETAS, nb=128, per_eta=PER_ETA, tag=tag)
    R['zero_amp_in_S'] = nzero; R['S_cap_supp'] = int(k - nzero)
    R['wc_gap_at_cut'] = tie
    log("%-14s |S|=%5d FR=%.4f w_S=%.12f  zero-amp in S=%d (%.1f%%)  wc gap at cut=%.3e  [%.1fs]"
        % (tag, k, k/D, R['w_S'], nzero, 100*nzero/max(k,1), tie, time.time()-t))
    hdr = "      %5s %11s %11s %11s %11s %11s %11s | %8s %8s %8s %8s %8s | %7s"
    log(hdr % ("eta","R0 truth","R1 pre-CS","R2 post-CS","R3 Mink-F","R4 Mink-G","R5 ThmB",
               "R1/R0","R2/R1","R3/R2","R4/R3","R5/R4","R5/R0"))
    def rr(a,b):
        return (a/b) if b > 1e-13 else float('nan')
    for r in R['rows']:
        log("      %5.2f %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e | %8.3f %8.3f %8.3f %8.3f %8.3f | %7.4g"
            % (r['eta'], r['R0'], r['R1'], r['R2'], r['R3'], r['R4'], r['R5'],
               rr(r['R1'],r['R0']), rr(r['R2'],r['R1']), rr(r['R3'],r['R2']), rr(r['R4'],r['R3']), rr(r['R5'],r['R4']), rr(r['R5'],r['R0'])))
        log("            cos(F,G)=%.4f  N[F]=%.4f (|phi|+|phi_S|=%.4f)  N[G]=%.4e  N[b]=%.4e  Lam/eta=%.4e"
            "  | quad chk |Rphi|^2->1: %.8f  NGquad/NG=%.8f"
            % (r['cos_FG'], r['NF'], R['norm_phi2']**0.5 + R['norm_phiS'], r['NG'], r['Nb'], r['Lam']/r['eta'],
               r['chk_normu'], r['chk_NG']))
    res['cases'].append(R)
    gc.collect()
out = os.path.join(OUT, 'ladder_L%d.json' % L)
json.dump(res, open(out, 'w'), indent=1)
log("WROTE " + out)
