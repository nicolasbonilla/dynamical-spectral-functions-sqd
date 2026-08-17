# -*- coding: utf-8 -*-
r"""Native fig_heron from the REAL IBM Heron spectral run (Spectral_Heron_READY.ipynb):
L=6 Hubbard A(omega) on ibm_fez, 50000 shots, job d9s16avpemts73ct6g8g.
HONEST framing: the recovered configurations span the COMPLETE (N+-1) sector (|S|=300/300), so the
sampled-subspace reconstruction is exact BY CONSTRUCTION (rel-L1=0). This is a coverage proof-of-
principle that the spectral pipeline runs end-to-end on a real device -- NOT a fidelity benchmark and
NOT a beyond-classical claim; at this system size full-sector coverage makes rel-L1=0 automatic.
We now plot the REAL saved hardware array A_hw (from the committed json, sourced from Google Drive)
as discrete markers on the exact target -- they coincide because of full-sector coverage."""
import json, numpy as np
d=json.load(open('data/heron_spectral.json'))
g=np.array(d['grid']); Aex=np.array(d['A_exact']); Ahw=np.array(d['A_hw'])
job=d['job_id']; shots=d['shots']; back=d['backend']; S=d['hw_S']; nsec=d['nsector']
# honest cross-check
assert np.allclose(Ahw,Aex,atol=1e-9), "A_hw should equal A_exact (full-sector coverage)"
with open('paper/figs/heron.dat','w') as f:
    f.write('w Aex Ahw\n')
    for x,ye,yh in zip(g,Aex,Ahw): f.write(f'{x:.4f} {ye:.5f} {yh:.5f}\n')
# subsample the hardware markers so they read as data points on the exact curve (not a solid line)
step=max(1,len(g)//34)
with open('paper/figs/heron_hw.dat','w') as f:
    f.write('w Ahw\n')
    for i in range(0,len(g),step): f.write(f'{g[i]:.4f} {Ahw[i]:.5f}\n')

tex=r"""\documentclass[border=6pt]{standalone}
\usepackage[T1]{fontenc}\usepackage{lmodern}\usepackage{amsmath}
\usepackage{tikz}\usetikzlibrary{arrows.meta}
\usepackage{pgfplots}\pgfplotsset{compat=1.18}
\definecolor{inkPrim}{HTML}{1B1B1B}\definecolor{inkMute}{HTML}{6E6E6E}
\definecolor{clExact}{HTML}{1F3A5F}\definecolor{coral}{HTML}{D55E00}\definecolor{gridClr}{HTML}{ECEAE5}
\definecolor{hwbox}{HTML}{F4EEE7}
\begin{document}\begin{tikzpicture}[font=\rmfamily]
\begin{axis}[width=12.8cm,height=7.8cm,
  xlabel={frequency \; $\omega-E_0$},ylabel={$A(\omega)$},
  xmin=__XMIN__,xmax=__XMAX__,ymin=0,ymajorgrids,y grid style={gridClr},
  axis line style={inkPrim,line width=0.7pt},tick style={inkPrim,line width=0.7pt},
  tick label style={font=\normalsize},label style={font=\large},
  legend style={font=\normalsize,draw=none,fill=white,fill opacity=0.7,text opacity=1,
    at={(0.985,0.52)},anchor=north east},legend cell align=left]
  \addplot[draw=none,fill=clExact,fill opacity=0.14,forget plot] table[x=w,y=Aex]{paper/figs/heron.dat} \closedcycle;
  \addplot[clExact,line width=1.3pt] table[x=w,y=Aex]{paper/figs/heron.dat};
  \addlegendentry{exact target $A(\omega)$}
  \addplot[only marks,mark=*,mark size=1.7pt,draw=coral,fill=coral]
     table[x=w,y=Ahw]{paper/figs/heron_hw.dat};
  \addlegendentry{IBM Heron \texttt{ibm\_fez} (measured)}
  \node[anchor=north east,align=left,font=\small,draw=coral,fill=hwbox,rounded corners=2.5pt,
        line width=0.8pt,inner sep=6pt] at (rel axis cs:0.985,0.965)
    {\textbf{Real quantum hardware} --- IBM Heron \texttt{ibm\_fez} (12 qubits),\\
     $50\,000$ shots, TREX $+$ twirling $+$ dynamical decoupling.\\
     Recovered configs span the \emph{complete} $(N{+}1)$ sector\\
     ($|\mathcal S|=300/300$), so the reconstruction is exact\\
     \emph{by construction} (rel-$L_1=0$): a coverage proof of\\
     principle, not a fidelity or beyond-classical claim.};
\end{axis}\end{tikzpicture}\end{document}
"""
xmin,xmax=g[0],g[-1]
open('fig_heron_native.tex','w',encoding='utf-8').write(
    tex.replace('__XMIN__',f'{xmin:.2f}').replace('__XMAX__',f'{xmax:.2f}'))
print(f'wrote fig_heron_native.tex ; REAL A_hw from job {job} ({back}, {shots} shots, |S|={S}/{nsec})')
print(f'  hardware markers: {len(range(0,len(g),step))} pts ; A_hw==A_exact (allclose): {np.allclose(Ahw,Aex)}')
