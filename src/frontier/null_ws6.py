# -*- coding: utf-8 -*-
"""The protocol, stated before looking: everything is chosen on TRAIN (L=6,8) and
evaluated on TEST (L=10,12), and the null is CONSTRAINED to the certificate's own
operating point (train median within a factor 2), so the two rules are compared at
the same error level and not at whichever level makes the null's curve flattest."""
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
def err_at(c, lvl):
    fr, om, rl = S[c]
    if lvl > om[0] or lvl < om[-1]: return None
    return math.exp(float(np.interp(float(np.interp(-math.log(lvl), -np.log(om), fr)), fr, np.log(rl))))
def cert(c):
    fr, om, rl = S[c]
    return math.exp(float(np.interp(FRS[c], fr, np.log(rl))))
TRAIN = [c for c in cells if c[0] in (6, 8)]
TEST  = [c for c in cells if c[0] in (10, 12)]
ct_tr = np.array([cert(c) for c in TRAIN]); ct_te = np.array([cert(c) for c in TEST])
m_tr = float(np.median(ct_tr))
print("CERTIFICATE  train median %.3e  spread x%.2f   |   TEST band %.3e..%.3e  spread x%.2f"
      % (m_tr, ct_tr.max()/ct_tr.min(), ct_te.min(), ct_te.max(), ct_te.max()/ct_te.min()))
print("Null constrained to train median in [%.2e, %.2e]  (factor 2 of the certificate's)\n"
      % (m_tr/2, m_tr*2))
for name, ps in (("NULL A  eta-blind, 1 free param", [None]),
                 ("NULL B  eta-aware, 2 free params", list(np.linspace(0, 5, 201)))):
    best = None
    for p in ps:
        for lg in np.arange(-6.0, -1.5, 0.01):
            c0 = 10.0**lg
            tr = [err_at(c, c0*(1.0 if p is None else c[1]**p)) for c in TRAIN]
            if any(v is None for v in tr): continue
            tr = np.array(tr); med = float(np.median(tr))
            if not (m_tr/2 <= med <= m_tr*2): continue
            sc = tr.max()/tr.min()
            if best is None or sc < best[0]: best = (sc, c0, p, med)
    if best is None:
        print("%s: NO setting reaches the certificate's operating point on TRAIN" % name); continue
    sc, c0, p, med = best
    te = [err_at(c, c0*(1.0 if p is None else c[1]**p)) for c in TEST]
    if any(v is None for v in te):
        print("%s: fitted (c=%.3e%s, train x%.2f, median %.2e) but NOT REACHABLE on a TEST cell"
              % (name, c0, "" if p is None else " eta^%.2f" % p, sc, med)); continue
    te = np.array(te)
    print("%s:\n   fitted c=%.3e%s  train median %.3e  train spread x%.2f"
          % (name, c0, "" if p is None else " * eta^%.2f" % p, med, sc))
    print("   -> TEST median %.3e  band %.3e..%.3e  spread x%.2f"
          % (np.median(te), te.min(), te.max(), te.max()/te.min()))
print("\nSCORE on the same eight held-out cells (L=10,12), matched operating point:")
print("   certificate  TEST spread x%.2f  (zero parameters)" % (ct_te.max()/ct_te.min()))
