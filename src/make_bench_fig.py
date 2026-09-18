# -*- coding: utf-8 -*-
r"""Native benchmark figure from the REAL multi-seed headtohead_ms.json (L=6 Hubbard ED).
Three panels U/t=4,8,12: spectral rel-L1 vs subspace size |S| for the quantum time-evolution selector
(coral, seed-averaged with spread band) vs classical polynomial Krylov (blue) and CIPSI-GF (pink).
At strong coupling the quantum selector pulls decisively ahead."""
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
import json
d=json.load(open('data/headtohead_ms.json'))
Us=[('U4.0','4'),('U8.0','8'),('U12.0','12')]
import numpy as np
for key,_ in Us:
    r=d[key]; ks=r['Ks']; Q=r['Q_mean']; Qs=r['Q_std']; K=r['K']; S=r['SCI']
    with open(f'paper/figs/bench_{key}.dat','w') as f:
        f.write('S Q Qlo Qhi K SCI\n')
        for i,k in enumerate(ks):
            qlo=max(Q[i]-Qs[i],1e-4); qhi=Q[i]+Qs[i]
            kk=K[i] if i<len(K) else ''
            ss=S[i] if i<len(S) else ''
            f.write(f'{k} {Q[i]:.4f} {qlo:.4f} {qhi:.4f} {kk} {ss}\n')

LEGK = r"\addlegendentry{classical Krylov}"
LEGS = r"\addlegendentry{classical CIPSI-GF}"
LEGQ = r"\addlegendentry{time-evolution selector}"
def panel(key,letter,legend=False):
    ylab = r'ylabel={relative-$L_1$ error},' if letter=='a' else ''
    ltn = r'legend to name=benchleg,' if legend else ''
    lk = LEGK if legend else ''; ls = LEGS if legend else ''; lq = LEGQ if legend else ''
    U=d[key]['U']
    s=rf"""\nextgroupplot[title={{({letter})\quad $U/t={U:.0f}$}},
  xlabel={{subspace size \; $|\mathcal{{S}}|$}}, {ylab} {ltn}
  ymode=log,ymin=0.03,ymax=1.15,xmin=0,xmax=270,
  ytick={{0.03,0.1,0.3,1}},yticklabels={{$0.03$,$0.1$,$0.3$,$1$}},
  ymajorgrids,y grid style={{gridClr}}]
  \addplot[name path=lo,draw=none,forget plot] table[x=S,y=Qlo]{{paper/figs/bench_{key}.dat}};
  \addplot[name path=hi,draw=none,forget plot] table[x=S,y=Qhi]{{paper/figs/bench_{key}.dat}};
  \addplot[coral,opacity=0.16,forget plot] fill between[of=lo and hi];
  \addplot[blue,line width=1.2pt,mark=square*,mark size=1.6pt,mark options={{fill=blue,draw=white}}]
     table[x=S,y=K]{{paper/figs/bench_{key}.dat}}; {lk}
  \addplot[pink,line width=1.2pt,mark=triangle*,mark size=2pt,mark options={{fill=pink,draw=white}}]
     table[x=S,y=SCI]{{paper/figs/bench_{key}.dat}}; {ls}
  \addplot[coral,line width=1.5pt,mark=*,mark size=1.8pt,mark options={{fill=coral,draw=white}}]
     table[x=S,y=Q]{{paper/figs/bench_{key}.dat}}; {lq}"""
    return s

LETTERS={'U4.0':'a','U8.0':'b','U12.0':'c'}
body="\n".join(panel(k,LETTERS[k],legend=(k=='U4.0')) for k,_ in Us)
tex=r"""\documentclass[11pt,border=6pt]{standalone}
\usepackage[T1]{fontenc}\usepackage{lmodern}\usepackage{amsmath}
\usepackage{pgfplots}\pgfplotsset{compat=1.18}\usepgfplotslibrary{groupplots,fillbetween}
\definecolor{inkPrim}{HTML}{1B1B1B}\definecolor{coral}{HTML}{D55E00}
\definecolor{blue}{HTML}{0072B2}\definecolor{pink}{HTML}{CC79A7}\definecolor{gridClr}{HTML}{ECEAE5}
\begin{document}\begin{tikzpicture}[font=\rmfamily]
\begin{groupplot}[group style={group size=3 by 1,horizontal sep=1.43cm},
  width=5.21cm,height=5.55cm,
  axis line style={inkPrim,line width=0.7pt},tick style={inkPrim,line width=0.7pt},
  tick label style={font=\normalsize},label style={font=\normalsize},title style={font=\normalsize},
  legend columns=-1,
  legend style={font=\normalsize,draw=none,/tikz/every even column/.append style={column sep=12pt}},
  legend cell align=left]
""" + body + r"""
\end{groupplot}
\node[anchor=north,yshift=-1.55cm] at (group c2r1.south) {\ref{benchleg}};
\end{tikzpicture}\end{document}
"""
open('fig_bench_native.tex','w',encoding='utf-8').write(tex)
print("wrote fig_bench_native.tex")
