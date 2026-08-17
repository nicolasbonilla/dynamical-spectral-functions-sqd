# -*- coding: utf-8 -*-
r"""Recolor akw_field_v2.png from viridis -> hot WITHOUT changing the data (invert viridis to scalar via
nearest-LUT, then apply hot). Also inject a smooth 32-point 'hot' colormap into fig_akw_v2_L12.tex colorbar."""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np
from scipy.spatial import cKDTree

img=plt.imread('/w/akw_field_v2.png')[...,:3]                    # H x W x 3, float [0,1] (viridis-coloured)
vir=matplotlib.colormaps['viridis'](np.linspace(0,1,256))[:,:3]
h,w,_=img.shape
_,idx=cKDTree(vir).query(img.reshape(-1,3))                      # nearest viridis index per pixel -> scalar
val=(idx/255.0).reshape(h,w)
hot=matplotlib.colormaps['hot'](val)[...,:3]
plt.imsave('/w/akw_field_v2.png', hot)
print('recoloured akw_field_v2.png viridis -> hot (data unchanged)')

# smooth 32-point hot colormap for the native colorbar (no quarter-steps)
cols=(np.array([matplotlib.colormaps['hot'](x)[:3] for x in np.linspace(0,1,32)])*255).round().astype(int)
cmstr=' '.join(f'rgb255=({r},{g},{b})' for r,g,b in cols)
p='/w/fig_akw_v2_L12.tex'; L=open(p).read().splitlines()
out=[]
for l in L:
    if 'colormap/viridis' in l:
        out.append(r'\pgfplotsset{colormap={mplhot}{'+cmstr+'}}')       # define smooth hot
        out.append(l.replace('colormap/viridis','colormap name=mplhot'))# use it
    else: out.append(l)
open(p,'w').write('\n'.join(out)+'\n')
print('fig_akw_v2_L12.tex colorbar -> smooth hot')
