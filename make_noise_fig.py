# -*- coding: utf-8 -*-
r"""Native fig_noise from the REAL noise_spectral.json (L=6 Hubbard, per-qubit bit-flip channel):
spectral rel-L1 vs bit-flip rate eps, naive post-selection vs S-CoRe configuration recovery."""
import json
d=json.load(open('noise_spectral.json'))
eps=d['eps']
with open('paper/figs/noise.dat','w') as f:
    f.write('eps na na_lo na_hi sc sc_lo sc_hi\n')
    for i,e in enumerate(eps):
        na,ns=d['naive'][i],d['naive_std'][i]; sc,ss=d['score'][i],d['score_std'][i]
        f.write(f'{e:.3f} {na:.4f} {max(na-ns,0):.4f} {na+ns:.4f} {sc:.4f} {max(sc-ss,0):.4f} {sc+ss:.4f}\n')

tex=r"""\documentclass[11pt,border=6pt]{standalone}
\usepackage[T1]{fontenc}\usepackage{lmodern}\usepackage{amsmath}
\usepackage{tikz}\usetikzlibrary{arrows.meta}
\usepackage{pgfplots}\pgfplotsset{compat=1.18}\usepgfplotslibrary{fillbetween}
\definecolor{inkPrim}{HTML}{1B1B1B}\definecolor{inkMute}{HTML}{6E6E6E}
\definecolor{coral}{HTML}{D1495B}\definecolor{teal}{HTML}{009E73}\definecolor{gridClr}{HTML}{ECEAE5}
\begin{document}\begin{tikzpicture}[font=\rmfamily]
\begin{axis}[width=11.2cm,height=7.2cm,
  xlabel={per-qubit bit-flip rate \; $\varepsilon$},ylabel={spectral rel-$L_1$ error},
  xmin=0,xmax=0.205,ymin=0,ymax=0.50,
  xtick={0,0.05,0.10,0.15,0.20},xticklabels={$0$,$5\%$,$10\%$,$15\%$,$20\%$},
  ymajorgrids,y grid style={gridClr},axis line style={inkPrim,line width=0.6pt},tick style={inkPrim},
  tick label style={font=\normalsize},label style={font=\normalsize},
  legend style={font=\normalsize,draw=none,fill=white,fill opacity=0.75,text opacity=1,
    at={(0.03,0.97)},anchor=north west},legend cell align=left]
  \addplot[name path=nlo,draw=none,forget plot] table[x=eps,y=na_lo]{paper/figs/noise.dat};
  \addplot[name path=nhi,draw=none,forget plot] table[x=eps,y=na_hi]{paper/figs/noise.dat};
  \addplot[coral,opacity=0.15,forget plot] fill between[of=nlo and nhi];
  \addplot[name path=slo,draw=none,forget plot] table[x=eps,y=sc_lo]{paper/figs/noise.dat};
  \addplot[name path=shi,draw=none,forget plot] table[x=eps,y=sc_hi]{paper/figs/noise.dat};
  \addplot[teal,opacity=0.15,forget plot] fill between[of=slo and shi];
  \addplot[coral,line width=1.5pt,mark=square*,mark size=2pt,mark options={fill=coral,draw=white}]
     table[x=eps,y=na]{paper/figs/noise.dat}; \addlegendentry{naive post-selection}
  \addplot[teal,line width=1.5pt,mark=*,mark size=2pt,mark options={fill=teal,draw=white}]
     table[x=eps,y=sc]{paper/figs/noise.dat}; \addlegendentry{S-CoRe configuration recovery}
  \draw[inkMute,dashed,line width=0.6pt] (axis cs:0.16,0) -- (axis cs:0.16,0.47);
  \node[inkPrim,font=\normalsize,anchor=west,align=left] at (axis cs:0.085,0.15)
     {recovery holds rel-$L_1\!\approx\!0.05$\\ up to $\varepsilon\approx16\%$};
  \node[inkPrim,font=\footnotesize,anchor=south east,align=right] at (axis cs:0.204,0.255)
     {naive post-selection\\ degrades with noise};
\end{axis}\end{tikzpicture}\end{document}
"""
open('fig_noise_native.tex','w',encoding='utf-8').write(tex)
print("wrote fig_noise_native.tex")
