# -*- coding: utf-8 -*-
"""Viridis raster of the dynamical structure factor S(q,omega) (native axes added in LaTeX).
S(q=0,w>0)=0 (density conservation) and S(q)=S(2pi-q); periodic cubic interpolation in q."""
# --- DEPOSIT PATHS (repaired 2026-09-18) -----------------------------------
# This script used to hard-code its output under '/w/'.  /w was the working
# directory of the Docker container the published runs were made in; outside that
# container the documented pipeline wrote nothing a reader could find, and data/
# was in fact repopulated BY HAND.  That made `make data` and the README recipe
# untrue.  Repaired: every path is now an ARGUMENT with a default RELATIVE TO THIS
# REPOSITORY, so a clean clone reproduces into its own tree.
#     read   <repo>/data/<name>      override with  --in  PATH
#     write  <repo>/data/<name>      override with  --out PATH
#     write  <repo>/build/<name>     for by-products that are NOT part of the deposit
# Paths only -- no physics and no computational default was changed here.
import os as _os, sys as _sys
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))

def _flag(name):
    """Value of a `--name VALUE` (or `-n VALUE`) command-line flag, else None."""
    a = _sys.argv[1:]
    for f in ('--' + name, '-' + name[0]):
        if f in a and a.index(f) + 1 < len(a):
            return a[a.index(f) + 1]
    return None

def _argv_positional():
    """argv[1:] with the --in/--out flags and their values removed."""
    a, keep, i = _sys.argv[1:], [], 0
    while i < len(a):
        if a[i] in ('--out', '-o', '--in', '-i'):
            i += 2
            continue
        keep.append(a[i]); i += 1
    return keep

def _outpath(name, sub='data'):
    """Absolute path to write `name` to: --out if given, else <repo>/<sub>/<name>."""
    p = _os.path.abspath(_flag('out') or _os.path.join(_REPO, sub, name))
    _os.makedirs(_os.path.dirname(p), exist_ok=True)
    return p

def _inpath(name, sub='data'):
    """Absolute path to read `name` from: --in if given, else <repo>/<sub>/<name>."""
    return _os.path.abspath(_flag('in') or _os.path.join(_REPO, sub, name))
# ---------------------------------------------------------------------------
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
L=12
d=json.load(open(_inpath(f'sqw_L{L}.json')))
wg=np.array(d['wg']); qs=np.array(sorted(float(x) for x in d['S'].keys()))
M=np.array([d['S'][f'{q:.5f}'] for q in qs])                 # (nq, nw)
qs_f=np.concatenate([[0.0],qs,[2.0]])                         # add q=0 and q=2 (=0)
M_f=np.vstack([np.zeros_like(M[0]),M,np.zeros_like(M[0])])
sel=(wg>=0)&(wg<=12); wgc=wg[sel]; M_f=M_f[:,sel]
cs=CubicSpline(qs_f,M_f,axis=0,bc_type='periodic')
qg=np.linspace(0,2,600); Z=np.clip(cs(qg),0,None).T          # (nw, nq)
VMAX=float(np.percentile(Z[Z>1e-6],99.5))
fig=plt.figure(figsize=(6.6,3.9),dpi=340); ax=fig.add_axes([0,0,1,1]); ax.axis('off')
ax.imshow(Z,origin='lower',extent=[0,2,0,12],aspect='auto',cmap='viridis',
          vmin=0,vmax=VMAX,interpolation='bilinear')
fig.savefig(_outpath('sqw_field.png', sub='build'),dpi=340)
print(f"L={L}  Z.max={Z.max():.4f}  VMAX={VMAX:.4f}")
