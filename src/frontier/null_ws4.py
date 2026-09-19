# -*- coding: utf-8 -*-
"""Final, careful version of the null test.  Reports the whole operating curve,
not one point, so the answer cannot be cherry-picked."""
import json, math, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
# --- deposit paths (repaired 2026-09-19; see docs/C3_FRONTIER.md) ----------
SRC = os.path.dirname(HERE)                       # release/src
ROOT = os.path.dirname(SRC)                       # repository root
OUT = os.environ.get('C3_OUT') or os.path.join(ROOT, 'data', 'c3_frontier')
TMP = os.environ.get('C3_TMP') or os.path.join(ROOT, 'build')
# --------------------------------------------------------------------------
NULL = os.path.join(OUT, 'null')
if not os.path.isdir(NULL):
    os.makedirs(NULL)
G = json.load(open(os.path.join(OUT, 'C3_grid.json')))
grid, front = G['grid'], G['frontier']
cells = [(f['L'], f['eta']) for f in front]
def series(L, eta):
    rows = sorted([r for r in grid if r['L'] == L and abs(r['eta']-eta) < 1e-9
                   and r['conv'] == 'CONVERGED' and r['relL1'] and r['relL1'] > 0],
                  key=lambda r: r['FR'])
    return (np.array([r['FR'] for r in rows]), np.array([1.0-r['w_S'] for r in rows]),
            np.array([r['relL1'] for r in rows]))
S = {c: series(*c) for c in cells}
def errs(fn):
    v = []
    for (L, eta) in cells:
        fr, om, rl = S[(L, eta)]
        lvl = fn(L, eta)
        if lvl > om[0] or lvl < om[-1]: return None
        FRx = float(np.interp(-math.log(lvl), -np.log(om), fr))
        v.append(math.exp(float(np.interp(FRx, fr, np.log(rl)))))
    return np.array(v)
Ec = np.array([math.exp(float(np.interp(f['FR_star'], S[(f['L'],f['eta'])][0],
              np.log(S[(f['L'],f['eta'])][2])))) for f in front])
print("CERTIFICATE (0 free parameters): median %.3e  band %.3e..%.3e  spread x%.2f\n"
      % (np.median(Ec), Ec.min(), Ec.max(), Ec.max()/Ec.min()))
print("NULL A: stop at a fixed (1-w_S).  Whole reachable curve, one line per level:")
print("  %-11s %-11s %-11s %-11s %s" % ("1-w_S","median err","min","max","spread"))
rows=[]
for lg in np.arange(-5.4, -2.0, 0.2):
    e = errs(lambda L, eta, l=10.0**lg: l)
    if e is None:
        print("  %-11.2e  NOT REACHABLE on the measured grid" % 10.0**lg); continue
    rows.append((10.0**lg, np.median(e), e.min(), e.max(), e.max()/e.min()))
    print("  %-11.2e %-11.3e %-11.3e %-11.3e x%.2f" % rows[-1])
print("\nNULL B: (1-w_S) < c*eta^p, 2 free parameters, scanned wide (c in 1e-7..1e-1, p in 0..5):")
best=[]
for p in np.linspace(0,5,101):
    for lg in np.arange(-7,-1.0,0.05):
        e = errs(lambda L,eta,c=10.0**lg,p=p: c*eta**p)
        if e is None: continue
        best.append((e.max()/e.min(), float(np.median(e)), 10.0**lg, p, e.min(), e.max()))
best.sort()
print("  best overall: spread x%.2f at median %.3e (c=%.2e, p=%.2f)" % (best[0][0],best[0][1],best[0][2],best[0][3]))
# best null constrained to the certificate's own operating decade
tgt = float(np.median(Ec))
near = [b for b in best if 0.5*tgt <= b[1] <= 2.0*tgt]
if near:
    print("  best with median within a factor 2 of the certificate's (%.2e): spread x%.2f at median %.3e (c=%.2e, p=%.2f)"
          % (tgt, near[0][0], near[0][1], near[0][2], near[0][3]))
else:
    print("  NO 2-parameter w_S rule reaches within a factor 2 of the certificate's operating point on this grid")
json.dump(dict(cert=dict(median=float(np.median(Ec)), spread=float(Ec.max()/Ec.min())),
               nullA=[list(map(float,r)) for r in rows],
               nullB_best=list(map(float,best[0])),
               nullB_matched=(list(map(float,near[0])) if near else None)),
          open(os.path.join(NULL,'NULL_WS4.json'),'w'), indent=1)
