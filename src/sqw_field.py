# -*- coding: utf-8 -*-
"""Viridis raster of the dynamical structure factor S(q,omega) (native axes added in LaTeX).
S(q=0,w>0)=0 (density conservation) and S(q)=S(2pi-q); periodic cubic interpolation in q."""
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
L=12
d=json.load(open(f'/w/sqw_L{L}.json'))
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
fig.savefig('/w/sqw_field.png',dpi=340)
print(f"L={L}  Z.max={Z.max():.4f}  VMAX={VMAX:.4f}")
