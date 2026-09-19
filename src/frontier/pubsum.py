# -*- coding: utf-8 -*-
"""Summarise the published-fraction points with the C3 classifier (fsum), verbatim."""
import os, sys, glob, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
# --- deposit paths (repaired 2026-09-19; see docs/C3_FRONTIER.md) ----------
SRC = os.path.dirname(HERE)                       # release/src
ROOT = os.path.dirname(SRC)                       # repository root
OUT = os.environ.get('C3_OUT') or os.path.join(ROOT, 'data', 'c3_frontier')
TMP = os.environ.get('C3_TMP') or os.path.join(ROOT, 'build')
# --------------------------------------------------------------------------
PUB = os.path.join(OUT, 'published_fraction')
if not os.path.isdir(PUB):
    os.makedirs(PUB)
sys.path.insert(0, HERE)
sys.dont_write_bytecode = True
import fsum

ETAS = [0.10, 0.15, 0.18, 0.25]
rows = []
for f in sorted(glob.glob(os.path.join(PUB, 'pt_L*_PUB.json'))):
    d = json.load(open(f))
    tr = fsum.top_rows(d)
    for eta in ETAS:
        r = tr.get(round(eta, 4))
        if r is None:
            continue
        pl = fsum.plateau_of(d, eta)
        rows.append(dict(L=d['meta']['L'], D=d['meta']['D'], S=d['S'],
                         FR_real=d['FR_real'], eta=eta, nl=d['nl'], w_S=d['w_S'],
                         B=r['bound_leak'], triv=r['bound_trivial'],
                         vac=r['bound_leak'] / r['bound_trivial'],
                         relL1=r.get('relL1'), io=r['in_over_out'],
                         relspread=pl['rel_spread'],
                         conv=fsum.conv_flag(r, pl['rel_spread']),
                         nontrivial=bool(r['bound_leak'] < r['bound_trivial']),
                         orth=d.get('orth_defect_sampled'), dw=d.get('dw'),
                         t=d['t_point'], src=os.path.basename(f)))
json.dump(rows, open(os.path.join(PUB, 'PUB_rows.json'), 'w'), indent=1, default=float)
print("%-3s %-8s %-6s %-5s %-9s %-9s %-8s %-10s %-9s %-9s %s" %
      ("L", "|S|", "FR", "eta", "B", "triv", "VAC", "relL1", "in/out", "relspr", "conv"))
for r in rows:
    print("%-3d %-8d %-6.4f %-5.2f %-9.4f %-9.4f x%-7.3f %-10.4e %-9.3e %-9.2e %s%s" %
          (r['L'], r['S'], r['FR_real'], r['eta'], r['B'], r['triv'], r['vac'],
           r['relL1'], r['io'], r['relspread'], r['conv'],
           "" if not r['nontrivial'] else "  <-- NON-VACUOUS"))
