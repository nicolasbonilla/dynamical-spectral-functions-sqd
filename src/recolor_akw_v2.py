# -*- coding: utf-8 -*-
r"""Recolor akw_field_v2.png from viridis -> hot WITHOUT changing the data (invert viridis to scalar via
nearest-LUT, then apply hot). Also inject a smooth 32-point 'hot' colormap into fig_akw_v2_L12.tex colorbar."""
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
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np
from scipy.spatial import cKDTree

img=plt.imread(_inpath('akw_field_v2.png', sub='build'))[...,:3]                    # H x W x 3, float [0,1] (viridis-coloured)
vir=matplotlib.colormaps['viridis'](np.linspace(0,1,256))[:,:3]
h,w,_=img.shape
_,idx=cKDTree(vir).query(img.reshape(-1,3))                      # nearest viridis index per pixel -> scalar
val=(idx/255.0).reshape(h,w)
hot=matplotlib.colormaps['hot'](val)[...,:3]
plt.imsave(_outpath('akw_field_v2.png', sub='build'), hot)
print('recoloured akw_field_v2.png viridis -> hot (data unchanged)')

# smooth 32-point hot colormap for the native colorbar (no quarter-steps)
cols=(np.array([matplotlib.colormaps['hot'](x)[:3] for x in np.linspace(0,1,32)])*255).round().astype(int)
cmstr=' '.join(f'rgb255=({r},{g},{b})' for r,g,b in cols)
p=_outpath('fig_akw_v2_L12.tex', sub='build'); L=open(p).read().splitlines()
out=[]
for l in L:
    if 'colormap/viridis' in l:
        out.append(r'\pgfplotsset{colormap={mplhot}{'+cmstr+'}}')       # define smooth hot
        out.append(l.replace('colormap/viridis','colormap name=mplhot'))# use it
    else: out.append(l)
open(p,'w').write('\n'.join(out)+'\n')
print('fig_akw_v2_L12.tex colorbar -> smooth hot')
