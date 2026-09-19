# -*- coding: utf-8 -*-
"""OUT-OF-SAMPLE test of the paper's ONE positive result (the calibration).

Paper claim (sec_5_body.tex:317-322): across the sixteen (L,eta) cells the crossing of
the trivial bound occurs at a TRUE relative L1 of 1.2e-3 .. 7.3e-3, median 2.8e-3.

Question nobody asked: is that band PREDICTIVE across size, or only descriptive?
Test: declare the band from the small sizes only, then check the large ones.
Pure post-processing of out/C3_grid.json -- zero new compute.
"""
import json, math, numpy as np, os
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

def interp_relL1(L, eta, FRx):
    """log-linear interpolation of the measured rel-L1 in FR, CONVERGED rows only."""
    rows = sorted([r for r in grid if r['L'] == L and abs(r['eta']-eta) < 1e-9
                   and r['conv'] == 'CONVERGED' and r['relL1'] is not None
                   and r['relL1'] > 0], key=lambda r: r['FR'])
    if len(rows) < 2: return None, 0, None
    x = np.array([r['FR'] for r in rows]); y = np.log(np.array([r['relL1'] for r in rows]))
    if FRx < x[0] or FRx > x[-1]:
        return None, len(rows), (float(x[0]), float(x[-1]))
    return float(math.exp(np.interp(FRx, x, y))), len(rows), (float(x[0]), float(x[-1]))

print("cell-by-cell: true rel-L1 AT the certificate's crossing of the trivial bound")
print("%-3s %-5s %-9s %-11s %-6s %-22s" % ("L","eta","FR*","relL1@FR*","nrows","FR range of CONV rows"))
cells = []
for f in front:
    L, eta, frs = f['L'], f['eta'], f['FR_star']
    v, n, rng = interp_relL1(L, eta, frs)
    cells.append(dict(L=L, eta=eta, FR_star=frs, relL1=v, nrows=n, rng=rng))
    print("%-3d %-5.2f %-9.4f %-11s %-6d %s" %
          (L, eta, frs, ("%.3e" % v) if v else "OUT-OF-RANGE", n, rng))

ok = [c for c in cells if c['relL1']]
vals = np.array([c['relL1'] for c in ok])
print("\nALL %d cells reproduced: min=%.2e  max=%.2e  median=%.2e  spread=%.2f"
      % (len(ok), vals.min(), vals.max(), np.median(vals), vals.max()/vals.min()))

for split in ([6], [6, 8], [6, 8, 10]):
    tr = np.array([c['relL1'] for c in ok if c['L'] in split])
    te = [c for c in ok if c['L'] not in split]
    if len(tr) == 0 or len(te) == 0: continue
    lo, hi = tr.min(), tr.max()
    inside = [c for c in te if lo <= c['relL1'] <= hi]
    print("\nDECLARE from L=%s (%d cells): band [%.2e, %.2e], spread x%.2f"
          % (split, len(tr), lo, hi, hi/lo))
    print("   TEST on L=%s (%d cells): %d/%d inside the band"
          % (sorted({c['L'] for c in te}), len(te), len(inside), len(te)))
    for c in te:
        fac = c['relL1']/np.median(tr)
        print("      L=%-3d eta=%.2f  relL1@FR*=%.3e   x%.2f of the declared median   %s"
              % (c['L'], c['eta'], c['relL1'], fac, "IN" if lo <= c['relL1'] <= hi else "OUT"))
json.dump(cells, open(os.path.join(NULL, 'CALIB_OOS.json'), 'w'), indent=1, default=float)
