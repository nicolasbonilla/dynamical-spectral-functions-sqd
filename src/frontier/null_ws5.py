# -*- coding: utf-8 -*-
"""OUT-OF-SAMPLE null test.  Fit the w_S rule on the SMALL sizes only, test on the
large ones, and score the certificate on the SAME test cells.  This removes the
in-sample advantage the null had in null_ws4."""
import json, math, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
# --- deposit paths (repaired 2026-09-19; see docs/C3_FRONTIER.md) ----------
SRC = os.path.dirname(HERE)                       # release/src
ROOT = os.path.dirname(SRC)                       # repository root
OUT = os.environ.get('C3_OUT') or os.path.join(ROOT, 'data', 'c3_frontier')
TMP = os.environ.get('C3_TMP') or os.path.join(ROOT, 'build')
# --------------------------------------------------------------------------
G = json.load(open(os.path.join(OUT, 'C3_grid.json')))
grid, front = G['grid'], G['frontier']
cells = [(f['L'], f['eta']) for f in front]
FRS = {(f['L'], f['eta']): f['FR_star'] for f in front}
def series(L, eta):
    rows = sorted([r for r in grid if r['L'] == L and abs(r['eta']-eta) < 1e-9
                   and r['conv'] == 'CONVERGED' and r['relL1'] and r['relL1'] > 0],
                  key=lambda r: r['FR'])
    return (np.array([r['FR'] for r in rows]), np.array([1.0-r['w_S'] for r in rows]),
            np.array([r['relL1'] for r in rows]))
S = {c: series(*c) for c in cells}
def err_at_level(c, lvl):
    fr, om, rl = S[c]
    if lvl > om[0] or lvl < om[-1]: return None
    FRx = float(np.interp(-math.log(lvl), -np.log(om), fr))
    return math.exp(float(np.interp(FRx, fr, np.log(rl))))
def cert_err(c):
    fr, om, rl = S[c]
    return math.exp(float(np.interp(FRS[c], fr, np.log(rl))))

TRAIN = [c for c in cells if c[0] in (6, 8)]
TEST  = [c for c in cells if c[0] in (10, 12)]
print("TRAIN = L=6,8  (%d cells)      TEST = L=10,12  (%d cells)\n" % (len(TRAIN), len(TEST)))

ct_tr = np.array([cert_err(c) for c in TRAIN]); ct_te = np.array([cert_err(c) for c in TEST])
print("CERTIFICATE (no fitting at all):")
print("   train band %.3e..%.3e x%.2f   TEST band %.3e..%.3e  spread x%.2f"
      % (ct_tr.min(), ct_tr.max(), ct_tr.max()/ct_tr.min(),
         ct_te.min(), ct_te.max(), ct_te.max()/ct_te.min()))

def fit_and_test(p=None, tag=""):
    best = None
    for lg in np.arange(-6.0, -1.5, 0.02):
        c0 = 10.0**lg
        tr = [err_at_level(c, c0 * (c[1]**p if p is not None else 1.0)) for c in TRAIN]
        if any(v is None for v in tr): continue
        tr = np.array(tr)
        sc = tr.max()/tr.min()
        if best is None or sc < best[0]: best = (sc, c0, tr)
    if best is None: return None
    sc, c0, tr = best
    te = [err_at_level(c, c0 * (c[1]**p if p is not None else 1.0)) for c in TEST]
    if any(v is None for v in te):
        print("   %s level %.3e fitted on TRAIN (x%.2f) -> NOT REACHABLE on some TEST cell" % (tag, c0, sc))
        return None
    te = np.array(te)
    print("   %s level c=%.3e%s   train x%.2f   TEST band %.3e..%.3e  spread x%.2f"
          % (tag, c0, ("" if p is None else " * eta^%.2f" % p), sc, te.min(), te.max(), te.max()/te.min()))
    return te.max()/te.min(), te

print("\nNULL A  (1 free parameter, eta-blind), fitted on TRAIN:")
a = fit_and_test(None, "(1-w_S) <")
print("\nNULL B  (2 free parameters, eta-aware), exponent scanned, level fitted on TRAIN:")
bb = None
for p in np.linspace(0, 5, 51):
    r = fit_and_test(p, "(1-w_S) <") if False else None
best_p = None
for p in np.linspace(0, 5, 101):
    best = None
    for lg in np.arange(-6.0, -1.5, 0.02):
        c0 = 10.0**lg
        tr = [err_at_level(c, c0*c[1]**p) for c in TRAIN]
        if any(v is None for v in tr): continue
        tr = np.array(tr); sc = tr.max()/tr.min()
        if best is None or sc < best[0]: best = (sc, c0)
    if best is None: continue
    te = [err_at_level(c, best[1]*c[1]**p) for c in TEST]
    if any(v is None for v in te): continue
    te = np.array(te)
    if best_p is None or best[0] < best_p[0]: best_p = (best[0], best[1], p, te)
if best_p:
    sc, c0, p, te = best_p
    print("   best TRAIN fit: c=%.3e * eta^%.2f (train x%.2f)"
          "  ->  TEST band %.3e..%.3e  spread x%.2f" % (c0, p, sc, te.min(), te.max(), te.max()/te.min()))
print("\nSCORE on the SAME eight test cells (lower = locates the error more tightly):")
print("   certificate           x%.2f  (zero parameters, threshold given by the theorem)" % (ct_te.max()/ct_te.min()))
if a: print("   w_S, 1 param fitted   x%.2f" % a[0])
if best_p: print("   w_S, 2 params fitted  x%.2f" % (best_p[3].max()/best_p[3].min()))
