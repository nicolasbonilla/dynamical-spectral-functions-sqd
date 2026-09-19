# -*- coding: utf-8 -*-
"""Robustness of the out-of-sample null test across every 2-vs-2 and 1-vs-3 split."""
import json, math, os, itertools, numpy as np
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
sizes = [6, 8, 10, 12]
print("%-14s %-14s %-9s %-9s %-9s %s" % ("TRAIN L", "TEST L", "cert", "nullA", "nullB", "winner"))
tally = {'cert': 0, 'null': 0, 'na': 0}
for k in (1, 2, 3):
    for tr_s in itertools.combinations(sizes, k):
        te_s = [s for s in sizes if s not in tr_s]
        TRAIN = [c for c in cells if c[0] in tr_s]; TEST = [c for c in cells if c[0] in te_s]
        ct_tr = np.array([cert(c) for c in TRAIN]); ct_te = np.array([cert(c) for c in TEST])
        m_tr = float(np.median(ct_tr)); cs = ct_te.max()/ct_te.min()
        res = {}
        for name, ps in (('A', [None]), ('B', list(np.linspace(0, 5, 101)))):
            best = None
            for p in ps:
                for lg in np.arange(-6.0, -1.5, 0.02):
                    c0 = 10.0**lg
                    tr = [err_at(c, c0*(1.0 if p is None else c[1]**p)) for c in TRAIN]
                    if any(v is None for v in tr): continue
                    tr = np.array(tr); med = float(np.median(tr))
                    if not (m_tr/2 <= med <= m_tr*2): continue
                    sc = tr.max()/tr.min()
                    if best is None or sc < best[0]: best = (sc, c0, p)
            if best is None: res[name] = None; continue
            te = [err_at(c, best[1]*(1.0 if best[2] is None else c[1]**best[2])) for c in TEST]
            res[name] = None if any(v is None for v in te) else (np.max(te)/np.min(te))
        vals = [v for v in res.values() if v is not None]
        win = 'cert' if (vals and cs < min(vals)) else ('null' if vals else 'n/a')
        tally['cert' if win == 'cert' else ('null' if win == 'null' else 'na')] += 1
        print("%-14s %-14s x%-8.2f %-9s %-9s %s" %
              (tr_s, tuple(te_s), cs,
               ("x%.2f" % res['A']) if res['A'] else "unreach",
               ("x%.2f" % res['B']) if res['B'] else "unreach", win.upper()))
print("\nsplits where the ZERO-parameter certificate beats BOTH fitted w_S rules: %d / %d"
      % (tally['cert'], sum(tally.values())))
