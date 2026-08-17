# -*- coding: utf-8 -*-
r"""G7 (decisive resource figure): the whole thesis in one plot, over 74 EXACT ground states
(38 molecular eq+diss, 30 Hubbard-sweep, 3 chain + 3 two-leg-ladder). (a) |S| vs chi: the cost tracks
the bond dimension in every family (chi lower-bounds |S|; the geminal witness chi=2,|S|=2^K is the
deliberately loose lower bound). (b) |S| vs F1: one-body magic does NOT predict the cost -- the Hubbard
branch runs the WRONG way (U up -> F1 up, |S| down) and the ladder pairs sit at the SAME F1 with very
different |S|. Reads resource_master.json (built from n19_suite/cost_vs_ent/ladder_vs_chain)."""
import json, numpy as np
d=json.load(open('resource_master.json'))
P=d['points']
fams=['molecule','hubbard','chain2','ladder']
col={'molecule':'rmMol','hubbard':'rmHub','chain2':'rmChn','ladder':'rmLad'}
mk ={'molecule':'*','hubbard':'square*','chain2':'triangle*','ladder':'diamond*'}
def coords(fam,xk,yk):
    s=[]
    for p in P:
        if p['family']!=fam: continue
        x=p['chi'] if xk=='chi' else p['F1']; y=max(p['S'],1)
        s.append(f"({x:.3f},{y})")
    return " ".join(s)
plotsA="".join(
 f"\\addplot[only marks,mark={mk[f]},mark size=1.7pt,{col[f]},mark options={{fill={col[f]},draw=white,line width=0.3pt}}] coordinates {{{coords(f,'chi','S')}}};\n"
 for f in fams)
plotsB="".join(
 f"\\addplot[only marks,mark={mk[f]},mark size=1.7pt,{col[f]},mark options={{fill={col[f]},draw=white,line width=0.3pt}}] coordinates {{{coords(f,'F1','S')}}};\n"
 for f in fams)
wit=" ".join(f"(2,{2**k})" for k in range(1,6))
rXS=d['pooled_spearman_chi_S']; rFS=d['pooled_spearman_F1_S']
tmpl=r"""% G7 master resource figure -- 74 exact ground states. make_resource_master_fig.py (resource_master.json).
\definecolor{rmMol}{HTML}{0072B2}\definecolor{rmHub}{HTML}{E69F00}\definecolor{rmChn}{HTML}{009E73}\definecolor{rmLad}{HTML}{D55E00}
\definecolor{rmInk}{HTML}{1B1B1B}
\begin{tikzpicture}
\begin{groupplot}[group style={group size=2 by 1, horizontal sep=1.9cm},
  width=8.0cm, height=6.2cm, paperaxis,
  title style={at={(0,1)},anchor=south west,font=\normalsize\bfseries,yshift=2pt},
  label style={font=\normalsize}, tick label style={font=\footnotesize},
  legend style={font=\scriptsize,draw=black!25,fill=white,fill opacity=0.9,text opacity=1},legend cell align=left]
% ---- (a) |S| vs chi : the cost governor ----
\nextgroupplot[title={(a)~cost tracks bond dimension $\chi$}, xmode=log,ymode=log,
  xlabel={minimal bond dimension $\chi$}, ylabel={determinant support $|\mathcal S|$},
  xmin=1.4,xmax=60,ymin=1,ymax=60000, ymajorgrids,xmajorgrids, grid style={draw=black!8},
  legend pos=north west]
__PLOTSA__\addlegendentry{molecules ($\times38$)}\addlegendentry{Hubbard sweep ($\times30$)}\addlegendentry{chain}\addlegendentry{2-leg ladder}
\addplot[rmInk,densely dashed,line width=0.9pt,mark=none] coordinates {__WIT__};
\node[rmInk,font=\scriptsize,anchor=west] at (axis cs:2.1,40) {witness: $\chi{=}2,\,|\mathcal S|{=}2^{K}$};
\node[rmInk,font=\scriptsize,anchor=south east] at (axis cs:55,1.5) {Spearman $\rho=__RXS__$};
% ---- (b) |S| vs F1 : magic does NOT predict cost ----
\nextgroupplot[title={(b)~one-body magic $\mathcal F_1$ does not}, ymode=log,
  xlabel={one-body magic $\mathcal F_1$}, ylabel={determinant support $|\mathcal S|$},
  xmin=-0.6,xmax=23,ymin=1,ymax=60000, ymajorgrids, grid style={draw=black!8},
  legend pos=north west]
__PLOTSB__\node[rmInk,font=\scriptsize,anchor=north east] at (axis cs:22.5,45000) {$\rho=__RFS__$};
\node[rmHub,font=\scriptsize,anchor=west] at (axis cs:6.2,150) {Hubbard: $U\!\uparrow\Rightarrow\mathcal F_1\!\uparrow,|\mathcal S|\!\downarrow$};
\node[rmLad,font=\scriptsize,anchor=south west] at (axis cs:13,1800) {ladder \& chain:};
\node[rmLad,font=\scriptsize,anchor=north west] at (axis cs:13,1400) {same $\mathcal F_1$, $|\mathcal S|$ jumps};
\end{groupplot}
\end{tikzpicture}
"""
out=(tmpl.replace('__PLOTSA__',plotsA).replace('__PLOTSB__',plotsB).replace('__WIT__',wit)
        .replace('__RXS__',f'{rXS:.2f}').replace('__RFS__',f'{rFS:.2f}'))
open('paper/figs/fig_resource_master_native.tex','w',encoding='utf-8').write(out)
print('wrote fig_resource_master_native.tex ; pooled rho chi/F1 =',round(rXS,2),round(rFS,2))
