# -*- coding: utf-8 -*-
"""DELIVERABLE (f): where does the FRONTIER sit at each rung of the proof?
The repo defines its published fraction by rel-L1 < thr (scaling_lanczos_mf.stage_frac, thr=0.05).
FR*_k(thr) = smallest fraction at which rung k of the ladder is below thr.
  FR*_0 = what the system actually needs        (PHYSICS)
  FR*_1 = after the pointwise step
  FR*_2 = after Cauchy-Schwarz in omega         (FR*_2 - FR*_1  IS the Cauchy-Schwarz displacement)
  FR*_5 = what the certificate can guarantee    (what Theorem B would force you to publish)
The displacement FR*_5 - FR*_0 is the proof-induced shift of the frontier."""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
# --- deposit paths (repaired 2026-09-19; see docs/C3_FRONTIER.md) ----------
SRC = os.path.dirname(HERE)                       # release/src
ROOT = os.path.dirname(SRC)                       # repository root
OUT = os.environ.get('C3_OUT') or os.path.join(ROOT, 'data', 'c3_frontier')
TMP = os.environ.get('C3_TMP') or os.path.join(ROOT, 'build')
# --------------------------------------------------------------------------
THRS = [0.10, 0.05, 0.02]
def frstar(frs, vals, thr):
    """smallest FR at which vals <= thr, by linear interpolation in log(val) vs FR (vals decreasing)."""
    frs = np.asarray(frs, float); vals = np.asarray(vals, float)
    o = np.argsort(frs); frs = frs[o]; vals = vals[o]
    for i in range(len(frs)):
        if vals[i] <= thr:
            if i == 0: return float(frs[0]), '<='
            x0, x1 = frs[i-1], frs[i]; y0, y1 = np.log(vals[i-1]), np.log(vals[i])
            return float(x0 + (np.log(thr)-y0)*(x1-x0)/(y1-y0)), 'ok'
    return float('nan'), '>'
print("%-3s %-6s %-6s | %-8s %-8s %-8s %-8s %-8s %-8s | %-9s %-9s" %
      ("L","eta","thr","FR*_0","FR*_1","FR*_2","FR*_3","FR*_4","FR*_5","d(CS)","d(TOTAL)"))
rows = []
for L in [int(x) for x in sys.argv[1:]]:
    R = json.load(open(os.path.join(OUT, 'ladder_L%d.json' % L)))
    n2 = R['norm_phi2']
    cases = [c for c in R['cases'] if c['tag'].startswith('FR=')]
    etas = [r['eta'] for r in cases[0]['rows']]
    for ie, eta in enumerate(etas):
        frs = [c['FR'] for c in cases]
        v = {k: [c['rows'][ie]['R%d' % k]/n2 for c in cases] for k in range(6)}
        for thr in THRS:
            f = {k: frstar(frs, v[k], thr) for k in range(6)}
            dcs = f[2][0]-f[1][0]; dtot = f[5][0]-f[0][0]
            print("%-3d %-6.2f %-6.2f | %-8s %-8s %-8s %-8s %-8s %-8s | %-9s %-9s" %
                  (L, eta, thr,
                   *["%.4f" % f[k][0] if f[k][1] == 'ok' else ("<%.2f" % f[k][0] if f[k][1] == '<=' else ">0.95")
                     for k in range(6)],
                   "%.4f" % dcs if np.isfinite(dcs) else "n/a",
                   "%.4f" % dtot if np.isfinite(dtot) else "n/a"))
            rows.append(dict(L=L, eta=eta, thr=thr, **{('FRstar%d' % k): f[k][0] for k in range(6)},
                             d_CS=dcs, d_total=dtot))
    print("-"*110)
json.dump(rows, open(os.path.join(OUT, 'verdict_frontier.json'), 'w'), indent=1)
