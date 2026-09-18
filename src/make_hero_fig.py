# -*- coding: utf-8 -*-
r"""Generate fig_hero.tex (native pgfplots, 2 panels) from the REAL FCI-verified n2_hero.json.
(a) the addition spectral function A(w) along the N2 triple-bond dissociation -- a single quasiparticle
peak at equilibrium splits into a strongly-correlated multiplet on stretch (waterfall, coloured by R);
(b) TWO genuinely independent quantities along the same coordinate: the fermionic magic F1 = 2 N_u
(EXACT identity, Eq. eq:collapse -- one-body magic IS the unpaired-electron count) RISES, while the
coherent quasiparticle weight Z (spectral sum rule) COLLAPSES. No tautology (F1 and N_u are the same
quantity by identity, so only ONE is drawn); the second curve Z is independent and ties (b) to (a).
Paper identity: vibrant molecular palette (covalent vermilion / ionic blue), neutral dark text, larger
fonts, no text over data. No hand-typed numbers."""
# --- DEPOSIT PATHS (repaired 2026-09-18, second pass) -----------------------
# This figure generator addresses every file it reads and writes by a path relative
# to the repository root ('data/...', 'paper/figs/...'), so it only ever worked when
# launched from that root and raised FileNotFoundError from anywhere else.  Rather
# than rewrite every literal -- which risks changing an output byte -- the process's
# working directory is anchored to the repository THIS FILE lives in.  Run from a
# copy of the tree (as src/check_figures.py does), it anchors to that copy, which is
# the behaviour that check wants.  Paths only; no figure content changed.
import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
# ---------------------------------------------------------------------------
import json, numpy as np
d=json.load(open('data/n2_hero.json'))
grid=np.array(d['grid']); P=d['points']
R=[p['R'] for p in P]; FAF=[p['FAF'] for p in P]; NU=[p['Nu'] for p in P]
SENT=[p['S_ent'] for p in P]; NPK=[p['npeaks'] for p in P]; ZW=[p['sumrule'] for p in P]
# sanity: the paper's central identity F1 = 2 N_u, verified to machine precision
assert all(abs(f-2*n)<1e-9 for f,n in zip(FAF,NU)), "F1 = 2 N_u identity broken!"
# waterfall data: grid + one A column per R, each normalized to unit peak for shape comparison
A=[np.array(p['A']) for p in P]
Amax=max(a.max() for a in A)
off=1.05                      # vertical offset between stacked curves (in normalized-A units)
with open('paper/figs/hero_aw.dat','w') as f:
    f.write('w '+' '.join(f'A{i}' for i in range(len(P)))+'\n')
    for j,wv in enumerate(grid):
        f.write(f'{wv:.4f} '+' '.join(f'{A[i][j]/Amax + i*off:.5f}' for i in range(len(P)))+'\n')
# resources vs R: magic F1 = 2 N_u (normalized to F1max = 12), and coherent weight Z (sum rule)
with open('paper/figs/hero_res.dat','w') as f:
    f.write('R magic z fabs nuabs\n')
    for i in range(len(P)):
        f.write(f'{R[i]:.3f} {NU[i]/6.0:.4f} {ZW[i]:.4f} {FAF[i]:.3f} {NU[i]:.3f}\n')

# sequential warm ramp (light -> dark maroon), one colour per R -- echoes the hot spectral identity
def ramp(t):  # t in [0,1] -> (r,g,b) light-orange to dark-maroon
    c0=np.array([0.99,0.85,0.62]); c1=np.array([0.40,0.0,0.12])
    return tuple((c0+(c1-c0)*t))
cols="".join(f"\\definecolor{{wf{i}}}{{rgb}}{{{r:.3f},{g:.3f},{b:.3f}}}\n"
             for i,(r,g,b) in enumerate(ramp(i/(len(P)-1)) for i in range(len(P))))
wfplots="".join(
  f"  \\addplot[wf{i},line width=1.1pt] table[x=w,y=A{i}]{{paper/figs/hero_aw.dat}};\n"
  f"  \\node[inkPrim,font=\\footnotesize,anchor=west] at (axis cs:1.18,{i*off+0.35:.3f}) {{{R[i]:.2f}}};\n"
  for i in range(len(P)))
Rmin,Rmax=min(R),max(R)
tmpl=r"""\documentclass[11pt,border=6pt]{standalone}
\usepackage[T1]{fontenc}\usepackage{lmodern}\usepackage{amsmath}
\usepackage{tikz}\usetikzlibrary{arrows.meta}
\usepackage{pgfplots}\pgfplotsset{compat=1.18}\usepgfplotslibrary{groupplots}
\definecolor{inkPrim}{HTML}{1B1B1B}\definecolor{inkMute}{HTML}{6E6E6E}
\definecolor{cMag}{HTML}{E8543A}  % vibrant vermilion -- one-body magic F1 = 2 N_u (covalent bond-breaking)
\definecolor{cZ}{HTML}{2E6FE0}    % vibrant blue      -- coherent quasiparticle weight Z
\definecolor{gridClr}{HTML}{ECEAE5}
__COLS__
\begin{document}
\begin{tikzpicture}[font=\rmfamily]
\begin{groupplot}[group style={group size=2 by 1, horizontal sep=1.9cm},
  width=7.6cm,height=8.6cm,
  axis line style={inkPrim,line width=0.7pt},tick style={inkPrim,line width=0.7pt},
  tick label style={font=\normalsize},label style={font=\normalsize},
  title style={font=\footnotesize,yshift=2pt}]
% ---- (a) spectral waterfall ----
\nextgroupplot[xlabel={addition energy \; $\omega-E_0$ \; (Ha)},
  ylabel={$A(\omega)$ \; (offset by $R$)}, ytick=\empty,
  title={(a)\quad quasiparticle $\to$ multiplet},
  xmin=-0.2,xmax=1.62, clip=false, ymin=-1.15, ymax=__YMAXA__,
  xtick={0,0.4,0.8,1.2}]
__WFPLOTS__
  \node[inkPrim,font=\footnotesize,anchor=west] at (axis cs:1.20,__TOPLBL__) {$R\,(\text{\AA})$};
  \node[font=\footnotesize\itshape,inkPrim,anchor=west]
     at (axis cs:-0.08,-0.72) {single quasiparticle (eq)};
  \node[font=\footnotesize\itshape,inkPrim,anchor=west]
     at (axis cs:-0.08,__MULTY__) {correlated multiplet (diss.)};
% ---- (b) two independent axes: magic F1=2N_u rises, coherent weight Z collapses ----
\nextgroupplot[xlabel={N--N bond length \; $R$ \; (\AA)},
  ylabel={fraction of maximum}, ylabel style={yshift=-2pt},
  title={(b)\quad magic $\mathcal F_1\!=\!2N_{\mathrm u}$ rises, weight $1{-}\langle n_p\rangle$ falls},
  xmin=__RMIN__,xmax=__RMAX__, ymin=0,ymax=1.15, ymajorgrids,y grid style={gridClr},
  ytick={0,0.2,0.4,0.6,0.8,1.0},
  legend style={font=\scriptsize,draw=inkMute,line width=0.4pt,fill=white,fill opacity=0.9,text opacity=1,
    at={(0.97,0.045)},anchor=south east,row sep=1.2pt},legend cell align=left]
  \draw[inkMute,line width=0.5pt,dashed] (axis cs:__RMIN__,1.0) -- (axis cs:__RMAX__,1.0);
  \node[inkPrim,font=\footnotesize,anchor=east] at (axis cs:__RMAX__,1.065)
     {$N_{\mathrm u}\!\to\!6$: two N($^{4}S$) atoms};
  % magic F1 = 2 N_u (identity) -- rising
  \addplot[cMag,line width=1.7pt,mark=*,mark size=2.6pt,mark options={fill=cMag,draw=white,line width=0.6pt}]
     table[x=R,y=magic]{paper/figs/hero_res.dat};
  \addlegendentry{fermionic magic $\mathcal F_1\!=\!2N_{\mathrm u}$}
  % coherent quasiparticle weight Z (spectral sum rule) -- collapsing
  \addplot[cZ,line width=1.7pt,mark=square*,mark size=2.4pt,mark options={fill=cZ,draw=white,line width=0.6pt}]
     table[x=R,y=z]{paper/figs/hero_res.dat};
  \addlegendentry{addition weight $1{-}\langle n_p\rangle$}
\end{groupplot}
\end{tikzpicture}
\end{document}
"""
NP=len(P)
out=(tmpl.replace('__COLS__',cols).replace('__WFPLOTS__',wfplots)
        .replace('__TOPLBL__',f'{(NP-1)*off+0.7:.2f}')
        .replace('__YMAXA__',f'{(NP-1)*off+2.0:.2f}')
        .replace('__MULTY__',f'{(NP-1)*off+1.5:.2f}')
        .replace('__RMIN__',f'{Rmin-0.05:.2f}').replace('__RMAX__',f'{Rmax+0.05:.2f}'))
open('fig_hero.tex','w',encoding='utf-8').write(out)
print('wrote fig_hero.tex ; FAF',round(min(FAF),2),'->',round(max(FAF),2),
      '| Nu',round(min(NU),2),'->',round(max(NU),2),
      '| Z',round(ZW[0],2),'->',round(ZW[-1],2),'| peaks',NPK)
print('IDENTITY F1 = 2 N_u verified to machine precision (only ONE drawn); Z is independent.')
